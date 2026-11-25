#!/usr/bin/env python3
"""MPFR golden generator for zenith-float-num/tests/mpfr (r6 item 6, Paul call 2026-11-24).

Run through tools/mpfr_golden/gen.sh, which checks the pinned team env (oracle-2026-10-30:
python-flint 0.9.0) and drives the tests with `--features mpfr-golden-gen`. Each test spawns
`gen.py --serve <test>` and sends one query per line:

    <op> <prec> <rnd> <arg> ...      rnd in N Z U D A (MPFR_RNDN, RNDZ, RNDU, RNDD, RNDA)

args are values (below), `i<int>` / `u<int>` integers, `s<len>:<text>` strings, or `q<n>/<d>`
fractions. The answer is `<value> <ternary>` (two pairs for sin_cos and sinh_cosh): the result
correctly rounded to <prec> bits with an UNBOUNDED exponent range, and the sign of
(rounded - exact). The test applies MPFR's exponent-range rules (mpfr_check_range) itself.

Values: nan, +inf, -inf, +0, -0, [+-]<odd hex>p<exp> (= hex * 2^exp exactly), and [+-]huge /
[+-]tiny for results whose exponent is beyond +-(2^62 + 2^32) (overflow / underflow in every range).

How results are computed (FLINT first, standing rule):
- add, sub, mul, div, fma, sqr, *_ui, remainder, parse, q_*: exact integer / Fraction arithmetic,
  then one rounding (a far-away addend becomes a sticky bit, which rounds identically);
- sqrt, cbrt, rootn_ui, hypot: exact integer roots with a sticky bit;
- pow: exact when the result is a dyadic number (integer powers, perfect roots), else arb;
- everything else: python-flint arb in a Ziv loop: the working precision doubles until both ends
  of the enclosure round to the same value, which lies strictly outside the enclosure (so the
  ternary is known). Exact results (MPFR's special values, exp(0), log(1), ...) are listed per
  function from the MPFR 4.2.2 manual and ISO C99 Annex F.
"""
import hashlib
import os
import re
import time
import sys
from fractions import Fraction

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, 'zenith-float-bound', 'scripts'))
import goldlib  # noqa: E402  (at import: popguard memory scope + the pinned team env, fail closed)

from flint import arb, arf, ctx, fmpz  # noqa: E402

HUGE = (1 << 62) + (1 << 32)
MAXWP = 1 << 22


class GenError(Exception):
    pass


# ------------------------------------------------------------------ values
# ('nan',) ('inf', neg) ('zero', neg) ('num', neg, M, E): (-1)^neg * M * 2^E, M > 0

NAN = ('nan',)


def inf(neg):
    return ('inf', neg)


def zero(neg):
    return ('zero', neg)


def num(neg, m, e):
    assert m > 0
    tz = (m & -m).bit_length() - 1
    return ('num', neg, m >> tz, e + tz)


def parse_val(s):
    if s == 'nan':
        return NAN
    neg = s[0] == '-'
    b = s[1:]
    if b == 'inf':
        return inf(neg)
    if b == '0':
        return zero(neg)
    m, e = b.split('p')
    return num(neg, int(m, 16), int(e))


def fmt(v, prec=None):
    k = v[0]
    if k == 'nan':
        return 'nan'
    s = '-' if v[1] else '+'
    if k == 'inf':
        return s + 'inf'
    if k == 'zero':
        return s + '0'
    if k in ('huge', 'tiny'):
        return s + k
    _, neg, m, e = v
    top = e + m.bit_length()
    if top > HUGE:
        return s + 'huge'
    if top < -HUGE:
        return s + 'tiny'
    return '%s%xp%d' % (s, m, e)


def sgn(neg, d):
    return -d if neg else d


def is_int(v):
    return v[0] == 'zero' or (v[0] == 'num' and v[3] >= 0)


def to_int(v):
    if v[0] == 'zero':
        return 0
    _, neg, m, e = v
    assert e >= 0
    return -(m << e) if neg else (m << e)


def is_odd_int(v):
    return v[0] == 'num' and v[3] == 0  # M is odd after normalization


def to_frac(v):
    _, neg, m, e = v
    f = Fraction(m << e) if e >= 0 else Fraction(m, 1 << -e)
    return -f if neg else f


# ------------------------------------------------------------------ rounding

def round_mag(neg, n, e, prec, rnd, sticky=False):
    """|v| = n * 2^e (exact; with sticky: strictly between n and n+1 units). -> (value, ternary)."""
    if n == 0:
        assert not sticky
        return None
    L = n.bit_length()
    if sticky:
        assert L >= prec + 2, 'sticky rounding needs prec+2 bits'
    if L <= prec:
        return num(neg, n, e), 0
    k = L - prec
    q = n >> k
    rem = n & ((1 << k) - 1)
    half = 1 << (k - 1)
    rb = rem >= half
    sb = sticky or (rem & (half - 1)) != 0
    e2 = e + k
    if not rb and not sb:
        return num(neg, q, e2), 0
    if rnd == 'N':
        up = rb and (sb or (q & 1))
    elif rnd == 'Z':
        up = False
    elif rnd == 'U':
        up = not neg
    elif rnd == 'D':
        up = neg
    else:
        up = True
    if up:
        q += 1
    return num(neg, q, e2), sgn(neg, 1 if up else -1)


def round_signed(n, e, prec, rnd, sticky=False, zero_neg=False):
    if n == 0:
        assert not sticky
        return zero(zero_neg), 0
    return round_mag(n < 0, abs(n), e, prec, rnd, sticky)


def round_frac(f, prec, rnd, zero_neg=False):
    if f == 0:
        return zero(zero_neg), 0
    neg = f < 0
    a, b = abs(f.numerator), f.denominator
    # scale so the quotient has >= prec+2 bits
    s = prec + 2 - (a.bit_length() - b.bit_length()) + 1
    if s >= 0:
        q, r = divmod(a << s, b)
    else:
        q, r = divmod(a, b << -s)
    return round_mag(neg, q, -s, prec, rnd, r != 0)


def round_val(v, prec, rnd):
    if v[0] != 'num':
        return v, 0
    return round_mag(v[1], v[2], v[3], prec, rnd)


# ------------------------------------------------------------------ exact arithmetic

def v2(n):
    """2-adic valuation of a nonzero integer."""
    return (n & -n).bit_length() - 1


def top(v):
    return v[3] + v[2].bit_length()


def add_num(x, y, prec, rnd):
    """x + y for two 'num' values (exact, one rounding)."""
    nx = -x[2] if x[1] else x[2]
    ny = -y[2] if y[1] else y[2]
    big, small = (x, y) if top(x) >= top(y) else (y, x)
    g = max(0, prec + 3 - big[2].bit_length())
    if top(small) < big[3] - g:
        # small lies strictly below big's extended last unit: a sticky bit
        n = (big[2] << (g + 1)) + (1 if big[1] == small[1] else -1)
        return round_mag(big[1], n, big[3] - g - 1, prec, rnd)
    e = min(x[3], y[3])
    s = (nx << (x[3] - e)) + (ny << (y[3] - e))
    return round_signed(s, e, prec, rnd, zero_neg=(rnd == 'D'))


def op_add(prec, rnd, x, y):
    if x[0] == 'nan' or y[0] == 'nan':
        return NAN, 0
    if x[0] == 'inf' or y[0] == 'inf':
        if x[0] == 'inf' and y[0] == 'inf' and x[1] != y[1]:
            return NAN, 0
        return (x if x[0] == 'inf' else y), 0
    if x[0] == 'zero' and y[0] == 'zero':
        neg = (x[1] and y[1]) or (x[1] != y[1] and rnd == 'D')
        return zero(neg), 0
    if x[0] == 'zero':
        return round_val(y, prec, rnd)
    if y[0] == 'zero':
        return round_val(x, prec, rnd)
    return add_num(x, y, prec, rnd)


def neg_v(v):
    return v if v[0] == 'nan' else (v[0], not v[1]) + tuple(v[2:])


def op_sub(prec, rnd, x, y):
    return op_add(prec, rnd, x, neg_v(y))


def mul_exact(x, y):
    """exact product of two non-NaN values, or NAN for 0*inf."""
    neg = x[1] != y[1]
    if 'inf' in (x[0], y[0]):
        if 'zero' in (x[0], y[0]):
            return NAN
        return inf(neg)
    if 'zero' in (x[0], y[0]):
        return zero(neg)
    return num(neg, x[2] * y[2], x[3] + y[3])


def op_mul(prec, rnd, x, y):
    if x[0] == 'nan' or y[0] == 'nan':
        return NAN, 0
    return round_val(mul_exact(x, y), prec, rnd)


def op_sqr(prec, rnd, x):
    return op_mul(prec, rnd, x, x)


def op_div(prec, rnd, x, y):
    if x[0] == 'nan' or y[0] == 'nan':
        return NAN, 0
    neg = x[1] != y[1]
    if x[0] == 'inf':
        return (NAN, 0) if y[0] == 'inf' else (inf(neg), 0)
    if y[0] == 'inf':
        return zero(neg), 0
    if y[0] == 'zero':
        return (NAN, 0) if x[0] == 'zero' else (inf(neg), 0)
    if x[0] == 'zero':
        return zero(neg), 0
    a, b = x[2], y[2]
    s = max(0, prec + 3 + b.bit_length() - a.bit_length())
    q, r = divmod(a << s, b)
    return round_mag(neg, q, x[3] - y[3] - s, prec, rnd, r != 0)


def op_fma(prec, rnd, a, b, c):
    if NAN in (a, b, c):
        return NAN, 0
    p = mul_exact(a, b)
    if p == NAN:
        return NAN, 0
    if p[0] == 'inf':
        if c[0] == 'inf' and c[1] != p[1]:
            return NAN, 0
        return p, 0
    if c[0] == 'inf':
        return c, 0
    return op_add(prec, rnd, p, c)


def uval(n):
    return zero(False) if n == 0 else num(n < 0, abs(n), 0)


def op_add_ui(prec, rnd, x, n):
    return op_add(prec, rnd, x, uval(n))


def op_div_ui(prec, rnd, x, n):
    return op_div(prec, rnd, x, uval(n))


def op_ui_div(prec, rnd, n, x):
    return op_div(prec, rnd, uval(n), x)


def op_remainder(prec, rnd, x, y):
    if x[0] == 'nan' or y[0] == 'nan' or x[0] == 'inf' or y[0] == 'zero':
        return NAN, 0
    if y[0] == 'inf' or x[0] == 'zero':
        return round_val(x, prec, rnd)
    mx, ex, my, ey = x[2], x[3], y[2], y[3]
    if top(x) < top(y) - 1:
        return round_val(x, prec, rnd)  # |x| < |y|/2: quotient 0
    if ex >= ey:
        X = (mx * pow(2, ex - ey, 2 * my)) % (2 * my)
        Y, e = my, ey
    else:
        X, Y, e = mx, my << (ey - ex), ex
        X %= 2 * Y
    odd = X >= Y
    r = X % Y
    if 2 * r > Y or (2 * r == Y and odd):
        r -= Y
    if x[1]:
        r = -r
    if r == 0:
        return zero(x[1]), 0
    return round_signed(r, e, prec, rnd)


def iroot(n, k):
    """floor(n^(1/k)) for n >= 0 (FLINT fmpz_root)."""
    return int(fmpz(n).root(k))


def root_num(neg, m, e, k, prec, rnd):
    """(-1)^neg * (m 2^e)^(1/k), m > 0 odd, correctly rounded."""
    if k >= m.bit_length() and k > 2:
        # An odd m > 1 that is a k-th power has at least 1.58 k bits, so the root is exact only for m = 1 and k | e,
        # and it is never a midpoint (a midpoint's k-th power is an odd k-th power too): the Ziv loop terminates.
        if m == 1 and e % k == 0:
            return num(neg, 1, e // k), 0
        x = num(neg, m, e)
        return ziv(lambda: (-1 if neg else 1) * (abs(to_arb(x)).log() / k).exp(), prec, rnd, 'rootn k=%d' % k)
    if k * (prec + 3) > (1 << 24):
        raise GenError('rootn: k=%d with a %d-bit mantissa not tabulated' % (k, m.bit_length()))
    s = max(0, k * (prec + 3) - m.bit_length())
    s += (e - s) % k
    n = m << s
    r = iroot(n, k)
    return round_mag(neg, r, (e - s) // k, prec, rnd, r ** k != n)


def op_sqrt(prec, rnd, x):
    if x[0] == 'nan':
        return NAN, 0
    if x[0] == 'zero':
        return x, 0
    if x[1]:
        return NAN, 0
    if x[0] == 'inf':
        return x, 0
    return root_num(False, x[2], x[3], 2, prec, rnd)


def op_rootn_ui(prec, rnd, x, k):
    if k == 0 or x[0] == 'nan':
        return NAN, 0
    if x[0] == 'zero':
        return zero(x[1] and k % 2 == 1), 0
    if x[1] and k % 2 == 0:
        return NAN, 0
    if x[0] == 'inf':
        return x, 0
    return root_num(x[1], x[2], x[3], k, prec, rnd)


def op_cbrt(prec, rnd, x):
    if x[0] == 'zero':
        return x, 0
    return op_rootn_ui(prec, rnd, x, 3)


def op_hypot(prec, rnd, x, y):
    if x[0] == 'inf' or y[0] == 'inf':
        return inf(False), 0
    if x[0] == 'nan' or y[0] == 'nan':
        return NAN, 0
    if x[0] == 'zero':
        return round_val((y[0], False) + tuple(y[2:]), prec, rnd)
    if y[0] == 'zero':
        return round_val((x[0], False) + tuple(x[2:]), prec, rnd)
    big, small = (x, y) if top(x) >= top(y) else (y, x)
    g = max(0, prec + 3 - big[2].bit_length())
    if 2 * top(small) < big[3] - g + top(big) - 2:
        # sqrt(B^2 + s^2) = B (1 + s^2/(2 B^2) + ...): s^2/(2B) is below B's extended last unit
        return round_mag(False, (big[2] << (g + 1)) + 1, big[3] - g - 1, prec, rnd)
    e = min(x[3], y[3])
    s = (x[2] << (x[3] - e)) ** 2 + (y[2] << (y[3] - e)) ** 2
    return root_num(False, s, 2 * e, 2, prec, rnd)


# ------------------------------------------------------------------ arb + Ziv

def to_arb(v):
    if v[0] == 'zero':
        return arb(0)
    _, neg, m, e = v
    a = arb(arf((fmpz(m), fmpz(e))))
    return -a if neg else a


def arb_ends(y):
    """exact (lo, hi) of an arb as (n, e) pairs with a common exponent e."""
    mm, me = y.mid().man_exp()
    rm, re_ = y.rad().man_exp()
    mm, me, rm, re_ = int(mm), int(me), int(rm), int(re_)
    if rm == 0:
        return (mm, me), (mm, me)
    if mm == 0:
        me = re_
    e = min(me, re_)
    c = mm << (me - e)
    r = rm << (re_ - e)
    return (c - r, e), (c + r, e)


def cmp_ne(a, b):
    """compare a=(n, e) with b=(n, e) exactly."""
    (na, ea), (nb, eb) = a, b
    e = min(ea, eb)
    x, y = na << (ea - e), nb << (eb - e)
    return (x > y) - (x < y)


def decide(y, prec, rnd):
    if not y.is_finite():
        return None
    lo, hi = arb_ends(y)
    if lo == hi:
        n, e = lo
        return round_signed(n, e, prec, rnd)
    if lo[0] <= 0 <= hi[0]:
        return None
    neg = hi[0] < 0
    ra, ta = round_mag(neg, abs(lo[0]), lo[1], prec, rnd)
    rb, tb = round_mag(neg, abs(hi[0]), hi[1], prec, rnd)
    if ra != rb:
        return None
    rv = ((-ra[2]) if neg else ra[2], ra[3])
    if cmp_ne(rv, hi) > 0:
        return ra, 1
    if cmp_ne(rv, lo) < 0:
        return ra, -1
    return None


def ziv(f, prec, rnd, what):
    wp = prec + 40
    while wp <= MAXWP:
        ctx.prec = wp
        r = decide(f(), prec, rnd)
        if r is not None:
            return r
        wp *= 2
    raise GenError('Ziv loop did not converge for %s at %d bits' % (what, MAXWP))


def near(neg, v_is_one_below, prec, rnd, m, e):
    """(-1)^neg * (m 2^e -+ epsilon): epsilon far below the last unit (sticky representative)."""
    g = max(0, prec + 3 - m.bit_length())
    n = (m << (g + 1)) + (-1 if v_is_one_below else 1)
    return round_mag(neg, n, e - g - 1, prec, rnd)


ONE = num(False, 1, 0)


def magnitude_log2(v):
    return top(v)


# ------------------------------------------------------------------ elementary functions

def arb1(name):
    def f(x):
        return getattr(to_arb(x), name)()
    return f


# f(x) for |x| < 2^-(prec+8): f(x) = x (1 +- O(x^2)) or 1 +- O(x): the correction is far below the last unit
# of the extended result, so the sticky representative rounds exactly as f(x). 'below'/'above' is the
# magnitude of f(x) against |x| (or against 1 for cos/cosh); for expm1/log1p it depends on the sign of x.
TINY = {
    'sin': 'below', 'tan': 'above', 'atan': 'below', 'asin': 'above', 'asinh': 'below',
    'atanh': 'above', 'sinh': 'above', 'tanh': 'below',
    'expm1': ('above', 'below'), 'log1p': ('below', 'above'),
}


def tz_exact(prec, rnd, x, fn, name):
    if x[0] == 'num' and top(x) < -(prec + 8):
        # x-relative cases: only when x lies on the (prec + 1)-bit grid of rounding boundaries, so that the
        # open interval between x and f(x) holds no boundary. A wider x can sit within the correction of a
        # boundary (tan of -(1 - 2^-7935) 2^-3628 at 768 bits, Up: f(x) passes -2^-3628), so it takes the Ziv loop.
        if name in TINY and x[2].bit_length() <= prec + 1:
            t = TINY[name]
            if isinstance(t, tuple):
                t = t[1] if x[1] else t[0]
            return near(x[1], t == 'below', prec, rnd, x[2], x[3])
        if name in ('cos', 'cosh'):
            return near(False, name == 'cos', prec, rnd, 1, 0)
        if name == 'exp':
            return near(False, x[1], prec, rnd, 1, 0)
    return ziv(lambda: fn(x), prec, rnd, '%s(%s)' % (name, fmt(x)))


def op_exp(prec, rnd, x):
    if x[0] == 'nan':
        return NAN, 0
    if x[0] == 'inf':
        return (zero(False), 0) if x[1] else (inf(False), 0)
    if x[0] == 'zero':
        return ONE, 0
    if top(x) > 64:
        return (('tiny', False), -1) if x[1] else (('huge', False), 1)
    return tz_exact(prec, rnd, x, arb1('exp'), 'exp')


def op_exp2(prec, rnd, x):
    if x[0] == 'nan':
        return NAN, 0
    if x[0] == 'inf':
        return (zero(False), 0) if x[1] else (inf(False), 0)
    if x[0] == 'zero':
        return ONE, 0
    if top(x) > 64:
        return (('tiny', False), -1) if x[1] else (('huge', False), 1)
    if is_int(x):
        return num(False, 1, to_int(x)), 0
    if top(x) < -(prec + 8):
        # |x ln 2| < 2^-(prec+8): 1 +- a sticky unit (sign of x), as for exp
        return near(False, x[1], prec, rnd, 1, 0)
    return ziv(lambda: (to_arb(x) * arb.const_log2()).exp(), prec, rnd, 'exp2')


def op_exp10(prec, rnd, x):
    if x[0] == 'nan':
        return NAN, 0
    if x[0] == 'inf':
        return (zero(False), 0) if x[1] else (inf(False), 0)
    if x[0] == 'zero':
        return ONE, 0
    if top(x) > 64:
        return (('tiny', False), -1) if x[1] else (('huge', False), 1)
    if is_int(x) and not x[1] and to_int(x) * 2.33 < prec + 64:
        k = to_int(x)
        return round_val(num(False, 5 ** k, k), prec, rnd)
    if top(x) < -(prec + 10):
        # |x ln 10| < 2^(top(x) + 2) < 2^-(prec+8): 1 +- a sticky unit (sign of x)
        return near(False, x[1], prec, rnd, 1, 0)
    return ziv(lambda: (to_arb(x) * arb(10).log()).exp(), prec, rnd, 'exp10')


def op_expm1(prec, rnd, x):
    if x[0] == 'nan':
        return NAN, 0
    if x[0] == 'inf':
        return (num(True, 1, 0), 0) if x[1] else (inf(False), 0)
    if x[0] == 'zero':
        return x, 0
    if top(x) > 64:
        if x[1]:
            return near(True, True, prec, rnd, 1, 0)  # -1 + tiny
        return ('huge', False), 1
    return tz_exact(prec, rnd, x, arb1('expm1'), 'expm1')


def log_common(prec, rnd, x, fn, name, exact=None):
    if x[0] == 'nan':
        return NAN, 0
    if x[0] == 'zero':
        return inf(True), 0
    if x[1]:
        return NAN, 0
    if x[0] == 'inf':
        return x, 0
    if x[2] == 1 and x[3] == 0:
        return zero(False), 0
    if exact is not None:
        r = exact(x)
        if r is not None:
            return round_val(uval(r), prec, rnd)
    return ziv(lambda: fn(to_arb(x)), prec, rnd, '%s(%s)' % (name, fmt(x)))


def op_log(prec, rnd, x):
    return log_common(prec, rnd, x, lambda a: a.log(), 'log')


def op_log2(prec, rnd, x):
    return log_common(prec, rnd, x, lambda a: a.log() / arb.const_log2(), 'log2',
                      lambda v: v[3] if v[2] == 1 else None)


def log10_exact(v):
    m, e = v[2], v[3]
    # x = 10^e = 5^e 2^e (e > 0): m = 5^e has about 2.32 e bits (checked before building 5^e)
    if 0 < e < (1 << 20) and abs(m.bit_length() - 2.321928 * e) < 4 and 5 ** e == m:
        return e
    return None


def op_log10(prec, rnd, x):
    return log_common(prec, rnd, x, lambda a: a.log() / arb(10).log(), 'log10', log10_exact)


def op_log1p(prec, rnd, x):
    if x[0] == 'nan':
        return NAN, 0
    if x[0] == 'zero':
        return x, 0
    if x[0] == 'inf':
        return (NAN, 0) if x[1] else (x, 0)
    if x[1] and x[2] == 1 and x[3] == 0:
        return inf(True), 0
    if x[1] and top(x) > 0 and to_frac(x) < -1:
        return NAN, 0
    return tz_exact(prec, rnd, x, arb1('log1p'), 'log1p')


def trig(name, zero_is_one=False):
    def op(prec, rnd, x):
        if x[0] in ('nan', 'inf'):
            return NAN, 0
        if x[0] == 'zero':
            return (ONE, 0) if zero_is_one else (x, 0)
        return tz_exact(prec, rnd, x, arb1(name), name)
    return op


op_sin = trig('sin')
op_cos = trig('cos', True)
op_tan = trig('tan')


def half_pi(neg, prec, rnd, mult=1, div=2):
    return ziv(lambda: (-1 if neg else 1) * arb.pi() * mult / div, prec, rnd, 'pi*%d/%d' % (mult, div))


def op_asin(prec, rnd, x):
    if x[0] in ('nan', 'inf'):
        return NAN, 0
    if x[0] == 'zero':
        return x, 0
    if top(x) > 1 or (top(x) == 1 and not (x[2] == 1 and x[3] == 0)):
        return NAN, 0
    if x[2] == 1 and x[3] == 0:
        return half_pi(x[1], prec, rnd)
    return tz_exact(prec, rnd, x, arb1('asin'), 'asin')


def op_acos(prec, rnd, x):
    if x[0] in ('nan', 'inf'):
        return NAN, 0
    if x[0] == 'zero':
        return half_pi(False, prec, rnd)
    if top(x) > 1 or (top(x) == 1 and not (x[2] == 1 and x[3] == 0)):
        return NAN, 0
    if x[2] == 1 and x[3] == 0:
        return (zero(False), 0) if not x[1] else half_pi(False, prec, rnd, 1, 1)
    return tz_exact(prec, rnd, x, arb1('acos'), 'acos')


def op_atan(prec, rnd, x):
    if x[0] == 'nan':
        return NAN, 0
    if x[0] == 'inf':
        return half_pi(x[1], prec, rnd)
    if x[0] == 'zero':
        return x, 0
    return tz_exact(prec, rnd, x, arb1('atan'), 'atan')


def op_atan2(prec, rnd, y, x):
    if x[0] == 'nan' or y[0] == 'nan':
        return NAN, 0
    yn = y[1]
    if y[0] == 'zero':
        if (x[0] in ('num', 'inf') and not x[1]) or (x[0] == 'zero' and not x[1]):
            return zero(yn), 0
        return half_pi(yn, prec, rnd, 2, 2)
    if y[0] == 'inf':
        if x[0] == 'inf':
            return half_pi(yn, prec, rnd, 1, 4) if not x[1] else half_pi(yn, prec, rnd, 3, 4)
        return half_pi(yn, prec, rnd)
    if x[0] == 'zero':
        return half_pi(yn, prec, rnd)
    if x[0] == 'inf':
        return (zero(yn), 0) if not x[1] else half_pi(yn, prec, rnd, 2, 2)
    return ziv(lambda: arb.atan2(to_arb(y), to_arb(x)), prec, rnd, 'atan2(%s,%s)' % (fmt(y), fmt(x)))


def op_sinh(prec, rnd, x):
    if x[0] in ('nan', 'inf', 'zero'):
        return x, 0
    if top(x) > 64:
        return ('huge', x[1]), sgn(x[1], 1)
    return tz_exact(prec, rnd, x, arb1('sinh'), 'sinh')


def op_cosh(prec, rnd, x):
    if x[0] == 'nan':
        return NAN, 0
    if x[0] == 'inf':
        return inf(False), 0
    if x[0] == 'zero':
        return ONE, 0
    if top(x) > 64:
        return ('huge', False), 1
    return tz_exact(prec, rnd, x, arb1('cosh'), 'cosh')


def op_tanh(prec, rnd, x):
    if x[0] in ('nan', 'zero'):
        return x, 0
    if x[0] == 'inf':
        return num(x[1], 1, 0), 0
    # 1 - |tanh x| < 2 e^(-2|x|) < 2^-(prec+8) once |x| >= 2^(top(x)-1) >= 0.35 (prec + 10) (ln 2 / 2 < 0.35)
    if top(x) > 64 or 20 * (1 << (top(x) - 1) if top(x) >= 1 else 0) >= 7 * (prec + 10):
        return near(x[1], True, prec, rnd, 1, 0)
    return tz_exact(prec, rnd, x, arb1('tanh'), 'tanh')


def op_asinh(prec, rnd, x):
    if x[0] in ('nan', 'inf', 'zero'):
        return x, 0
    return tz_exact(prec, rnd, x, arb1('asinh'), 'asinh')


def op_acosh(prec, rnd, x):
    if x[0] == 'nan':
        return NAN, 0
    if x[0] == 'inf':
        return (NAN, 0) if x[1] else (x, 0)
    if x[0] == 'zero' or x[1] or top(x) < 1:
        return NAN, 0
    if x[2] == 1 and x[3] == 0:
        return zero(False), 0
    return tz_exact(prec, rnd, x, arb1('acosh'), 'acosh')


def op_atanh(prec, rnd, x):
    if x[0] in ('nan', 'zero'):
        return x, 0
    if x[0] == 'inf':
        return NAN, 0
    if x[2] == 1 and x[3] == 0:
        return inf(x[1]), 0
    if top(x) >= 1:
        return NAN, 0
    return tz_exact(prec, rnd, x, arb1('atanh'), 'atanh')


# ------------------------------------------------------------------ pow

def pow_class(x, y):
    """'huge' / 'tiny' when |x|^y is certainly beyond every exponent range (|y log2|x|| > HUGE + 64, a 128-bit arb
    enclosure decides), else None. x and y are finite, nonzero, |x| != 1."""
    ctx.prec = 128
    m, e = x[2], x[3]
    if top(x) in (0, 1) and e < 0:
        # |x| in [1/2, 2): log1p of the exact |x| - 1 (a 128-bit log of |x| itself is only good to 2^-128 absolute)
        lx = arb(arf((fmpz(m - (1 << -e)), fmpz(e)))).log1p()
    else:
        lx = abs(to_arb(x)).log()
    t = to_arb(y) * lx / arb.const_log2()
    lim = arb(HUGE + 64)
    if t > lim:
        return 'huge'
    if t < -lim:
        return 'tiny'
    return None


def int_pow(prec, rnd, x, yv):
    """x^y for finite nonzero x and an integer value y != 0 (y = +-M 2^E, E >= 0, possibly huge)."""
    yodd = yv[3] == 0
    neg = x[1] and yodd
    m, e = x[2], x[3]
    if m == 1 and e == 0:
        return num(neg, 1, 0), 0
    c = pow_class(x, yv)
    if c == 'huge':
        return ('huge', neg), sgn(neg, 1)
    if c == 'tiny':
        return ('tiny', neg), sgn(neg, -1)
    n = to_int(yv)
    if m == 1:
        return num(neg, 1, e * n), 0
    if abs(n) * m.bit_length() <= (1 << 22):
        p = m ** abs(n)
        if n > 0:
            return round_val(num(neg, p, e * n), prec, rnd)
        s = prec + 3 + p.bit_length()
        q, r = divmod(1 << s, p)
        return round_mag(neg, q, -s - e * (-n), prec, rnd, r != 0)
    return ziv(lambda: to_arb(x) ** fmpz(n), prec, rnd, 'pow_int')


def op_pow(prec, rnd, x, y):
    if y[0] == 'zero':
        return ONE, 0
    if x == ONE:
        return ONE, 0
    if x[0] == 'nan' or y[0] == 'nan':
        return NAN, 0
    yodd = is_odd_int(y)
    yneg = y[1]
    if x[0] == 'zero':
        if yneg:
            return inf(x[1] and yodd), 0
        return zero(x[1] and yodd), 0
    if y[0] == 'inf':
        if x[1] and x[2] == 1 and x[3] == 0:
            return ONE, 0
        absx_lt1 = top(x) <= 0
        if yneg:
            return (inf(False), 0) if absx_lt1 else (zero(False), 0)
        return (zero(False), 0) if absx_lt1 else (inf(False), 0)
    if x[0] == 'inf':
        if x[1]:
            if yneg:
                return zero(yodd), 0
            return inf(yodd), 0
        return (zero(False), 0) if yneg else (inf(False), 0)
    if is_int(y):
        return int_pow(prec, rnd, x, y)
    if x[1]:
        return NAN, 0
    # y = my 2^ey, ey < 0
    my = -y[2] if y[1] else y[2]
    d = -y[3]
    m, e = x[2], x[3]
    if m == 1:
        # 2^(e y) with y = my 2^-d is dyadic iff 2^d divides e my (valuations: no 2^d is built, d may be ~2^62)
        if e == 0 or v2(e) + v2(abs(my)) >= d:
            return num(False, 1, (e * my) >> d), 0
    elif d <= 30 and e % (1 << d) == 0:
        r = iroot(m, 1 << d)
        if r ** (1 << d) == m:
            return int_pow(prec, rnd, num(False, r, e >> d), uval(my))
    c = pow_class(x, y)
    if c == 'huge':
        return ('huge', False), 1
    if c == 'tiny':
        return ('tiny', False), -1
    # x^y = exp(y ln x) with |y ln x| < 2^-(prec+8): 1 +- a sticky unit (the sign of y ln x), as for exp.
    # |ln x| <= (|top(x)| + 2) for every positive x (<= 1 on [1/2, 2)), so log2|y ln x| <= top(y) + bits(|top(x)| + 2).
    if top(y) + (abs(top(x)) + 2).bit_length() < -(prec + 8):
        z_pos = (not y[1]) == (top(x) >= 1)
        return near(False, not z_pos, prec, rnd, 1, 0)
    return ziv(lambda: to_arb(x) ** to_arb(y), prec, rnd, 'pow(%s,%s)' % (fmt(x), fmt(y)))


def op_pow_ui(prec, rnd, x, n):
    if n == 0:
        return ONE, 0
    return op_pow(prec, rnd, x, uval(n))


def op_pow_si(prec, rnd, x, n):
    if n == 0:
        return ONE, 0
    return op_pow(prec, rnd, x, uval(n))


# ------------------------------------------------------------------ special functions

def is_nonpos_int(x):
    return x[0] == 'zero' or (x[0] == 'num' and x[1] and x[3] >= 0)


def op_gamma(prec, rnd, x):
    if x[0] == 'nan':
        return NAN, 0
    if x[0] == 'inf':
        return (NAN, 0) if x[1] else (x, 0)
    if x[0] == 'zero':
        return inf(x[1]), 0
    if is_nonpos_int(x):
        return NAN, 0
    if is_int(x) and to_int(x) <= 100000:
        f = 1
        for k in range(2, to_int(x)):
            f *= k
        return round_val(num(False, f, 0), prec, rnd)
    if not x[1] and top(x) > 64:
        return ('huge', False), 1
    return tz_exact(prec, rnd, x, arb1('gamma'), 'gamma')


def op_lngamma(prec, rnd, x):
    if x[0] == 'nan':
        return NAN, 0
    if x[0] == 'inf' or is_nonpos_int(x):
        return inf(False), 0
    if x in (ONE, num(False, 1, 1)):
        return zero(False), 0
    if x[1]:
        f = to_frac(x)
        k = (-f).__floor__()  # -x in (k, k+1)
        if k % 2 == 0:  # -2k'-1 < x < -2k': Gamma(x) < 0
            return NAN, 0
        return ziv(lambda: abs(to_arb(x).gamma()).log(), prec, rnd, 'lngamma')
    return tz_exact(prec, rnd, x, arb1('lgamma'), 'lngamma')


def op_digamma(prec, rnd, x):
    if x[0] == 'nan':
        return NAN, 0
    if x[0] == 'inf':
        return (NAN, 0) if x[1] else (x, 0)
    if x[0] == 'zero':
        return inf(not x[1]), 0
    if is_nonpos_int(x):
        return NAN, 0
    return tz_exact(prec, rnd, x, arb1('digamma'), 'digamma')


def op_erf(prec, rnd, x):
    if x[0] in ('nan', 'zero'):
        return x, 0
    if x[0] == 'inf':
        return num(x[1], 1, 0), 0
    # 1 - |erf x| = erfc|x| <= e^(-x^2) < 2^-(prec+8) once x^2 >= 2^(2 top(x) - 2) >= 0.7 (prec + 8) (ln 2 < 0.7)
    if top(x) > 40 or (top(x) >= 1 and 10 * (1 << (2 * top(x) - 2)) >= 7 * (prec + 8)):
        return near(x[1], True, prec, rnd, 1, 0)
    return tz_exact(prec, rnd, x, arb1('erf'), 'erf')


def op_erfc(prec, rnd, x):
    if x[0] == 'nan':
        return NAN, 0
    if x[0] == 'inf':
        return (num(False, 1, 1), 0) if x[1] else (zero(False), 0)
    if x[0] == 'zero':
        return ONE, 0
    if top(x) > 40:
        if x[1]:
            return near(False, True, prec, rnd, 1, 1)  # 2 - tiny
        return ('tiny', False), -1
    return tz_exact(prec, rnd, x, arb1('erfc'), 'erfc')


def op_eint(prec, rnd, x):
    if x[0] == 'nan':
        return NAN, 0
    if x[0] == 'inf':
        return (zero(True), 0) if x[1] else (x, 0)
    if x[0] == 'zero':
        return inf(True), 0
    if not x[1] and top(x) > 64:
        return ('huge', False), 1
    return tz_exact(prec, rnd, x, arb1('ei'), 'eint')


def op_ai(prec, rnd, x):
    if x[0] == 'nan':
        return NAN, 0
    if x[0] == 'inf':
        return zero(False), 0
    return tz_exact(prec, rnd, x, arb1('airy_ai'), 'ai')


def op_jn(prec, rnd, n, x):
    if x[0] == 'nan':
        return NAN, 0
    if x[0] == 'inf':
        return zero(False), 0
    if x[0] == 'zero':
        if n == 0:
            return ONE, 0
        return zero(abs(n) % 2 == 1 and (x[1] != (n < 0))), 0
    return ziv(lambda: to_arb(x).bessel_j(n), prec, rnd, 'jn(%d,%s)' % (n, fmt(x)))


def op_yn(prec, rnd, n, x):
    if x[0] == 'nan':
        return NAN, 0
    if x[0] == 'zero':
        return inf(not (n < 0 and n % 2 == 1)), 0
    if x[1]:
        return NAN, 0
    if x[0] == 'inf':
        return zero(False), 0
    return ziv(lambda: to_arb(x).bessel_y(n), prec, rnd, 'yn(%d,%s)' % (n, fmt(x)))


def op_gamma_inc(prec, rnd, a, x):
    if a[0] == 'nan' or x[0] == 'nan':
        return NAN, 0
    if x[0] == 'zero':
        return op_gamma(prec, rnd, a)
    if a[0] != 'num' or x[0] != 'num':
        raise GenError('gamma_inc special arguments %s %s not tabulated' % (fmt(a), fmt(x)))
    return ziv(lambda: to_arb(x).gamma_upper(to_arb(a)), prec, rnd, 'gamma_inc')


def op_const_pi(prec, rnd):
    return ziv(lambda: arb.pi(), prec, rnd, 'pi')


def op_const_log2(prec, rnd):
    return ziv(lambda: arb.const_log2(), prec, rnd, 'log2')


def op_const_euler(prec, rnd):
    return ziv(lambda: arb.const_euler(), prec, rnd, 'euler')


def pair(f, g):
    def op(prec, rnd, x):
        return f(prec, rnd, x) + g(prec, rnd, x)
    return op


op_sin_cos = pair(op_sin, op_cos)
op_sinh_cosh = pair(op_sinh, op_cosh)


def parse_text(text, base):
    t = text.strip()
    neg = t.startswith('-')
    t = t.lstrip('+-')
    if t.lower() in ('inf', '@inf@'):
        return inf(neg), None
    if t.lower() in ('nan', '@nan@'):
        return NAN, None
    sep = 'e' if base == 10 else 'p'
    mant, _, ex = t.lower().partition(sep) if sep in t.lower() else (t, '', '0')
    if base == 2 and 'e' in t.lower():
        mant, _, ex = t.lower().partition('e')
    ip, _, fp = mant.partition('.')
    digits = (ip + fp) or '0'
    n = int(digits, base)
    f = Fraction(n) * Fraction(base) ** (int(ex) - len(fp))
    if f == 0:
        return zero(neg), None
    return None, (-f if neg else f)


def op_parse(base):
    def op(prec, rnd, s):
        v, f = parse_text(s, base)
        if v is not None:
            return v, 0
        return round_frac(f, prec, rnd)
    return op


op_parse10 = op_parse(10)
op_parse2 = op_parse(2)

OPS = {k[3:]: v for k, v in globals().items() if k.startswith('op_') and callable(v)}


def qparse(a):
    n, d = a[1:].split('/')
    return Fraction(int(n), int(d))


def answer(line):
    parts = line.split(' ')
    op, prec, rnd = parts[0], int(parts[1]), parts[2]
    if op.startswith('q_'):
        a, b = qparse(parts[3]), qparse(parts[4])
        r = {'q_add': a + b, 'q_sub': a - b, 'q_mul': a * b, 'q_div': a / b}[op]
        return 'q%d/%d' % (r.numerator, r.denominator)
    args = []
    i = 3
    rest = line.split(' ', 3)[3] if len(parts) > 3 else ''
    while rest:
        if rest[0] == 's':
            ln, _, tail = rest[1:].partition(':')
            ln = int(ln)
            args.append(tail[:ln])
            rest = tail[ln + 1:]
            continue
        tok, _, rest = rest.partition(' ')
        if tok[0] in 'iu':
            args.append(int(tok[1:]))
        else:
            args.append(parse_val(tok))
    fn = OPS.get(op)
    if fn is None:
        raise GenError('unknown op ' + op)
    res = fn(prec, rnd, *args)
    out = []
    for k in range(0, len(res), 2):
        v, t = res[k], res[k + 1]
        if v[0] == 'num':
            assert v[2].bit_length() <= prec, (line[:200], v)
        out.append('%s %d' % (fmt(v), t))
    return ' '.join(out)


def serve(name):
    # Optional progress file (MPFR_GOLDEN_PROGRESS=<dir>): "<queries answered> <last op>", rewritten every 5 s, so a
    # supervisor can tell a long stream from a stalled one.
    pdir = os.environ.get('MPFR_GOLDEN_PROGRESS')
    count, last = 0, 0.0
    cache = {}  # a repeated query gets the same answer (the stream stores it as an @k back-reference)
    for line in sys.stdin:
        line = line.rstrip('\n')
        count += 1
        if pdir and time.time() - last > 5:
            last = time.time()
            with open(os.path.join(pdir, name), 'w') as f:
                f.write('%d %s\n' % (count, line.split(' ', 1)[0]))
        try:
            a = cache.get(line)
            if a is None:
                a = cache[line] = answer(line)
        except Exception as e:  # reported to the test, which fails loudly
            a = 'ERROR %s: %s' % (type(e).__name__, str(e)[:300])
            sys.stderr.write('gen.py %s: %s on %s\n' % (name, a, line[:300]))
            if pdir:
                with open(os.path.join(pdir, name + '.errors'), 'a') as f:
                    f.write(line + '\n')
        sys.stdout.write(a + '\n')
        sys.stdout.flush()


def finish(d):
    """Writes SHA256SUMS and PROVENANCE for every stream file in d: <name>.gold, or its parts <name>.gold.<k> (sorted, LF)."""
    names = sorted(f for f in os.listdir(d) if re.fullmatch(r'.+\.gold(\.[0-9]+)?', f))
    if not names:
        sys.exit('gen.py --finish: no .gold files in %s' % d)
    sums = []
    for n in names:
        with open(os.path.join(d, n), 'rb') as f:
            sums.append('%s  %s' % (hashlib.sha256(f.read()).hexdigest(), n))
        goldlib.record_provenance(os.path.join(d, n))
    with open(os.path.join(d, 'SHA256SUMS'), 'wb') as f:
        f.write(''.join(l + '\n' for l in sums).encode('ascii'))
    print('gen.py: %d goldens, team-env=%s' % (len(names), goldlib.TEAM_ENV))


if __name__ == '__main__':
    if len(sys.argv) >= 2 and sys.argv[1] == '--check-env':
        print('gen.py: team-env=%s' % goldlib.TEAM_ENV)
    elif len(sys.argv) >= 3 and sys.argv[1] == '--finish':
        finish(sys.argv[2])
    elif len(sys.argv) >= 3 and sys.argv[1] == '--serve':
        serve(sys.argv[2])
    elif len(sys.argv) >= 2 and sys.argv[1] == '--answer':
        print(answer(sys.argv[2]))
    else:
        sys.exit('usage: gen.py --check-env | --serve <test> | --answer "<query>" | --finish <dir>')
