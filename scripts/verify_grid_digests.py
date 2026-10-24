#!/usr/bin/env python3
"""Verify every entry of a grid-digest dump (ZF_DIGEST_DUMP, see `dd` in ext.rs grid_digest_tests) with mpmath.

Each line: op, p, rm, ints, operands (repr|inexact), '=', result repr, inexact, inexact_accumulated.
Finite results are checked against the exact value (exact integer arithmetic with an infinitesimal
sticky for unbounded exponent gaps) or against mpmath at p+256 and p+512 bits (both must decide the
same rounding), correctly rounded at p bits under the mode, with the exponent range [EMIN, EMAX]
(subnormals at EMIN, quantum 2^(EMIN-p)). Special results are tabulated per operand class for review.
usage: verify_grid_digests.py DUMPFILE [max_examples]
"""
import re, sys, collections
import mpmath
from mpmath import mp, mpf
from mpmath.libmp import from_man_exp

EMIN = -(2**62 - 1)
MAXP = (1 << 60) - 1
EMAX = 2**62 - 1
NUM = re.compile(r'Value\(ExactNumNumber \{ e: (-?\d+), s: (Pos|Neg), m: Mantissa \{ m: WordBuf \{ inner: \[([^\]]*)\] \}, n: (\d+) \}, inexact: (true|false), op_inexact: (None|Some\((?:true|false)\))')

class V:
    __slots__ = ('kind', 'neg', 'M', 'E', 'acc', 'n', 'raw')
def parse(s):
    v = V(); v.raw = s; v.acc = False; v.n = 0
    if 'Inf(' in s:
        v.kind = 'inf'; v.neg = 'Inf(Neg)' in s; return v
    if 'NaN(' in s:
        v.kind = 'nan'; v.neg = False; v.M = s[s.index('NaN('):].split(')')[0] + ')'; return v
    m = NUM.search(s)
    if not m:
        raise ValueError(s)
    e, sg, ws, n, acc, _ = m.groups()
    words = [int(w.strip().strip('"')) for w in ws.split(',') if w.strip()]
    M = 0
    for i, w in enumerate(words):
        M |= w << (64 * i)
    v.neg = sg == 'Neg'; v.acc = acc == 'true'; v.n = int(n)
    if M == 0:
        v.kind = 'zero'; v.M = 0; v.E = 0; return v
    v.kind = 'fin'; tz = (M & -M).bit_length() - 1
    v.M = M >> tz; v.E = int(e) - 64 * len(words) + tz  # n is the bit length; the words are left-aligned
    return v

def cls(v):
    if v.kind == 'nan': return 'N'
    return {'zero': 'Z', 'fin': 'F', 'inf': 'I'}[v.kind] + ('-' if v.neg else '+')

# exact dyadic with infinitesimal: value = sgn*(M*2^E) + st*eps, st in {-1,0,1} (absolute direction)
class D:
    __slots__ = ('neg', 'M', 'E', 'st')
    def __init__(s, neg, M, E, st=0):
        if M:
            tz = (M & -M).bit_length() - 1; M >>= tz; E += tz
        s.neg, s.M, s.E, s.st = neg, M, E, st
def dv(v):
    return D(v.neg, v.M if v.kind == 'fin' else 0, v.E if v.kind == 'fin' else 0)
def top(d): return d.E + d.M.bit_length()
GAP = 20000
CUR_RM = ['ToEven']
class Ambig(Exception): pass
def dadd(a, b):
    if a.M == 0 and b.M == 0:
        if a.st and b.st and a.st != b.st: raise Ambig('eps')
        return D(a.neg if a.neg == b.neg else CUR_RM[0] == 'Down', 0, 0, a.st or b.st)
    if a.M == 0:
        if a.st and a.st != b.st and b.st: raise Ambig('eps')
        return D(b.neg, b.M, b.E, b.st or a.st)
    if b.M == 0: return dadd(b, a)
    if top(a) < top(b): a, b = b, a
    if top(a) - top(b) > GAP:
        st = -1 if b.neg else 1
        if a.st and a.st != st: raise Ambig('eps')
        return D(a.neg, a.M, a.E, st)
    if a.st or b.st:
        if a.st and b.st and a.st != b.st: raise Ambig('eps')
    E = min(a.E, b.E)
    x = (-1 if a.neg else 1) * (a.M << (a.E - E)) + (-1 if b.neg else 1) * (b.M << (b.E - E))
    return D(x < 0 if x else CUR_RM[0] == 'Down', abs(x), E, a.st or b.st)
def dneg(a): return D(not a.neg, a.M, a.E, -a.st)
def dmul(a, b):
    if a.st or b.st: raise Ambig('eps-mul')
    return D(a.neg != b.neg, a.M * b.M, a.E + b.E)
def ddiv(a, b, bits):
    if b.st: raise Ambig('eps-div')
    k = max(0, bits + b.M.bit_length() - a.M.bit_length() + 4)
    q, r = divmod(a.M << k, b.M)
    neg = a.neg != b.neg
    st = (-1 if neg else 1) if r else 0
    if a.st:
        if st and st != a.st * (-1 if b.neg else 1): raise Ambig('eps-div2')
        st = st or a.st * (-1 if b.neg else 1)
    return D(neg, q, a.E - b.E - k, st)
def iroot(x, n):
    if x == 0: return 0
    r = 1 << ((x.bit_length() + n - 1) // n)
    while True:
        y = ((n - 1) * r + x // r ** (n - 1)) // n
        if y >= r: break
        r = y
    while r ** n > x: r -= 1
    while (r + 1) ** n <= x: r += 1
    return r
def droot(a, n, bits):
    if a.st: raise Ambig('eps-root')
    s = n * (bits + 4) - a.M.bit_length()
    s = max(s, 0); E = a.E - s
    s += E % n; E -= E % n
    X = a.M << s
    r = iroot(X, n)
    st = (-1 if a.neg else 1) if r ** n != X else 0
    return D(a.neg, r, E // n, st)

def rnd(d, p, rm):
    """-> ('zero',neg) | ('fin',neg,M,E) | ('inf',neg), inexact"""
    neg = d.neg
    if d.M == 0:
        if d.st: raise Ambig('pure-eps')
        return ('zero', neg), False
    st = d.st * (-1 if neg else 1)  # relative to magnitude
    e = top(d)
    qe = (e - p) if e >= EMIN else (EMIN - p)
    sh = qe - d.E
    if sh <= 0:
        R = d.M << (-sh)
        c = 'exact' if st == 0 else 'lo' if st > 0 else 'below'
    elif sh > d.M.bit_length() + 1:
        R, c = 0, 'lo'
    else:
        R = d.M >> sh; rem = d.M & ((1 << sh) - 1); half = 1 << (sh - 1)
        if rem == 0: c = 'exact' if st == 0 else ('lo' if st > 0 else 'below')
        elif rem < half: c = 'lo'
        elif rem > half: c = 'hi'
        else: c = 'half' if st == 0 else ('hi' if st > 0 else 'lo')
    if c == 'below':
        if d.M == 1 and e - 1 >= EMIN and sh <= 0:
            qe = e - 1 - p; R = 1 << p
        R -= 1; c = 'hi_full'
    inexact = c != 'exact'
    up = False
    if inexact:
        if rm == 'Up': up = not neg
        elif rm == 'Down': up = neg
        elif rm == 'ToZero' or rm == 'None': up = False
        elif rm == 'FromZero': up = True
        elif rm in ('ToEven', 'ToOdd'):
            if c in ('hi', 'hi_full'): up = True
            elif c == 'half': up = (R & 1 == 1) if rm == 'ToEven' else (R & 1 == 0)
            else: up = False
        else: raise ValueError(rm)
    if up: R += 1
    if R == 0:
        return ('zero', neg), inexact
    if qe + R.bit_length() > EMAX:
        away = {'Up': not neg, 'Down': neg, 'ToZero': False, 'None': False}.get(rm, True)
        if away: return ('inf', neg), True
        return ('fin', neg, (1 << p) - 1, EMAX - p), True
    tz = (R & -R).bit_length() - 1
    return ('fin', neg, R >> tz, qe + tz), inexact

def tle(x, y):
    # x <= y for ('zero',neg) / ('fin',neg,M,E) tuples, without materialising 2^E
    sx = 0 if x[0] == 'zero' else (-1 if x[1] else 1)
    sy = 0 if y[0] == 'zero' else (-1 if y[1] else 1)
    if sx != sy or sx == 0: return sx <= sy
    tx, ty = x[3] + x[2].bit_length(), y[3] + y[2].bit_length()
    if tx != ty: mag = tx < ty
    else:
        k = min(x[3], y[3]); mx, my = x[2] << (x[3] - k), y[2] << (y[3] - k)
        if mx == my: return True
        mag = mx < my
    return mag if sx > 0 else not mag
def got_t(v):
    if v.kind == 'fin': return ('fin', v.neg, v.M, v.E)
    if v.kind == 'zero': return ('zero', v.neg)
    if v.kind == 'inf': return ('inf', v.neg)
    return ('nan',)

def to_mpf(v):
    if v.kind == 'fin': return mp.make_mpf(from_man_exp(-v.M if v.neg else v.M, v.E))
    if v.kind == 'zero': return mpf(0)
    if v.kind == 'inf': return mpf('-inf') if v.neg else mpf('inf')
    return mpf('nan')

def mp_to_D(y):
    man, ex = mpmath.libmp.normalize(*y._mpf_[:3], 1 << 30, 'n')[1:3] if False else (None, None)
    s, m, e, bc = y._mpf_
    return D(bool(s), int(m), int(e))

# taylor coefficients cache for tiny-argument functions
TAY = {}
def tiny_eval(fname, f, u, p, rm):
    key = fname
    if key not in TAY:
        with mp.workprec(400):
            cs = mpmath.taylor(f, 0, 7)
        TAY[key] = [c if abs(c) > mpf(2) ** -300 else mpf(0) for c in cs]
    cs = TAY[key]
    nz = [k for k, c in enumerate(cs) if c != 0]
    a = nz[0]; b = nz[1] if len(nz) > 1 else None
    ud = dv(u)
    sgn_u = -1 if u.neg else 1
    def dyadic(c):
        s, m, e, bc = c._mpf_
        return (bc <= 64), D(bool(s), int(m), int(e))
    ok, cd = dyadic(cs[a])
    if not ok:
        raise Ambig('tiny-nondyadic')
    T = cd if a == 0 else dmul(cd, ud) if a == 1 else None
    if T is None: raise Ambig('tiny-a>1')
    stsign = 0 if b is None else (1 if cs[b] > 0 else -1) * (sgn_u ** b)
    T.st = stsign
    return T

def mpfun(op, xs, ints):
    x = xs[0]
    F = {
        'sin': mpmath.sin, 'cos': mpmath.cos, 'tan': mpmath.tan, 'atan': mpmath.atan,
        'sinh': mpmath.sinh, 'cosh': mpmath.cosh, 'tanh': mpmath.tanh, 'atanh': mpmath.atanh,
        'exp': mpmath.exp, 'exp2': lambda t: mpmath.power(2, t), 'exp10': lambda t: mpmath.power(10, t),
        'ln': mpmath.log, 'log2': lambda t: mpmath.log(t, 2), 'log10': lambda t: mpmath.log10(t),
        'log1p': mpmath.log1p,
    }
    if op in F: return F[op](x)
    if op == 'logb': return mpmath.log(x) / mpmath.log(xs[1])
    if op == 'pow': return mpmath.power(x, xs[1])
    if op in JAC:
        return jac(op, x, xs[1])
    raise KeyError(op)
JAC = {'am', 'sn', 'cn', 'dn', 'cd', 'ns', 'nc', 'nd', 'sc', 'sd', 'cs', 'ds', 'dc'}
def jac(op, u, m):
    if op == 'am':
        if m == 0: return u
        return mpmath.quad(lambda t: mpmath.ellipfun('dn', t, m), [0, u])
    return mpmath.ellipfun(op, u, m)

def decide(op, xs, ints, p, rm):
    last = None
    for extra in (256, 512):
        wp = p + extra
        with mp.workprec(wp + 64):
            try:
                y = mpfun(op, xs, ints)
            except (ValueError, ZeroDivisionError):
                return ('nan',), None
            if isinstance(y, mpmath.mpc):
                if abs(y.imag) > mpf(2) ** (-wp) * (abs(y.real) + mpf(2) ** EMIN):
                    return ('nan',), None
                y = y.real
            if mpmath.isnan(y): return ('nan',), None
            if mpmath.isinf(y): return ('inf', y < 0), False
            if y == 0:
                return ('zero', None), None
            if op == "tanh" and abs(y) == 1 and xs[0] != 0:
                # |tanh x| < 1 for every finite x: a saturated value lies just inside 1 (1 - |tanh x| ~ 2e^-2|x|)
                return rnd(D(y < 0, 1, 0, 1 if y < 0 else -1), p, rm)
            yd = mp_to_D(y)
            if yd.M.bit_length() <= 128:
                return rnd(yd, p, rm)  # exact: 128+ trailing zero bits at p+256 only for an exact value
            lo = D(yd.neg, yd.M, yd.E, -1 if not yd.neg else 1)
            hi = D(yd.neg, yd.M, yd.E, 1 if not yd.neg else -1)
            # error bound: shift interval by 2 ulps at wp
            ulp_e = top(yd) - wp
            e1 = dadd(D(yd.neg, yd.M, yd.E), D(True, 1, ulp_e + 1))
            e2 = dadd(D(yd.neg, yd.M, yd.E), D(False, 1, ulp_e + 1))
            r1, i1 = rnd(e1, p, rm); r2, i2 = rnd(e2, p, rm)
            if r1 == r2:
                return r1, True
            last = (r1, r2)
    raise Ambig('interval %s' % (last,))

def exact_op(op, a, ints, p, rm, raw):
    """exact D for algebraic ops, or None"""
    if op in ('add', 'op+', 'add_full'): return dadd(dv(a[0]), dv(a[1]))
    if op in ('sub', 'op-', 'sub_full'): return dadd(dv(a[0]), dneg(dv(a[1])))
    if op in ('mul', 'op*', 'mul_full', 'sqr'):
        return dmul(dv(a[0]), dv(a[1] if len(a) > 1 else a[0]))
    if op in ('div', 'op/'): return ddiv(dv(a[0]), dv(a[1]), p + 8)
    if op == 'reciprocal': return ddiv(D(False, 1, 0), dv(a[0]), p + 8)
    if op == 'ldexp': d = dv(a[0]); d.E += ints[0]; return d
    if op == 'powi':
        b = dv(a[0]); n = ints[0]
        if b.st: return None
        return D(b.neg and n % 2 == 1, b.M ** n, b.E * n)
    if op == 'nth_root':
        # a zero root keeps the sign for odd n and for the square root (MPFR sqrt and rootn_ui)
        if a[0].kind == 'zero': return D(a[0].neg and (ints[0] % 2 == 1 or ints[0] == 2), 0, 0)
        if a[0].neg and ints[0] % 2 == 0: return 'nan'
        return droot(dv(a[0]), ints[0], p + 8)
    if op == 'fma': return dadd(dmul(dv(a[0]), dv(a[1])), dv(a[2]))
    if op == 'hypot':
        x, y = dv(a[0]), dv(a[1]); x.neg = y.neg = False
        if x.M == 0 or y.M == 0:
            return x if y.M == 0 else y
        if abs(top(x) - top(y)) > GAP:
            big = x if top(x) > top(y) else y
            return D(False, big.M, big.E, 1)
        s = dadd(dmul(x, x), dmul(y, y))
        return droot(s, 2, p + 8)
    return None

def expected(op, p, rm, ints, args, res_t):
    """-> (expected tuple, expected op-inexact or None) ; raises Ambig"""
    ws = max((x.n for x in args if x.kind in ('fin', 'zero')), default=64)
    if op.startswith('op'): p, rm = 128, 'ToEven'
    if op in ('add_full', 'sub_full', 'mul_full'):
        d = exact_op(op, args, ints, p, rm, None)
        if d.st: raise Ambig('full-eps')
        if d.M == 0: return ('zero', None), False
        if top(d) < EMIN and (EMIN - top(d)) + 64 * ((d.M.bit_length() + 63) // 64) > MAXP:
            # beyond MAX_PRECISION the library cannot widen: it rounds at the word width onto the subnormal grid
            return rnd(d, 64 * ((d.M.bit_length() + 63) // 64), rm)
        return ('fin', d.neg, d.M, d.E), False
    if op == 'round':
        d = dv(args[0])
        if d.M == 0: return ('zero', args[0].neg), False
        n = ints[0]; qe = -n
        if d.E >= qe: return ('fin', d.neg, d.M, d.E), False
        # emulate rounding at fixed quantum by choosing p so that e-p == qe
        pp = top(d) - qe
        if pp <= 0:
            pp = 1; d2 = d
            # value below half quantum or between: handle via integer
            R = 0; rem_ratio = d
        r, ix = rnd_fixed(d, qe, rm)
        return r, ix
    d = exact_op(op, args, ints, p, rm, None)
    if d == 'nan': return ('nan',), None
    if d is not None:
        r, ix = rnd(d, p, rm)
        return r, ix
    # transcendental
    xs = [to_mpf(x) for x in args]
    if op in TINY and args[0].kind == 'fin' and top(dv(args[0])) < -5000:
        f = tiny_fun(op, xs)
        if f == 'pole':
            pass
        else:
            T = tiny_eval(op + str(xs[1:]), f, args[0], p, rm)
            return rnd(T, p, rm)
    return decide(op, xs, ints, p, rm)

def rnd_fixed(d, qe, rm):
    neg = d.neg
    sh = qe - d.E
    if sh > d.M.bit_length() + 1:
        R, rem, half = 0, 1, 2
    else:
        R = d.M >> sh; rem = d.M & ((1 << sh) - 1); half = 1 << (sh - 1)
    c = 'exact' if rem == 0 else 'lo' if rem < half else 'hi' if rem > half else 'half'
    ix = c != 'exact'; up = False
    if ix:
        if rm == 'Up': up = not neg
        elif rm == 'Down': up = neg
        elif rm in ('ToZero', 'None'): up = False
        elif rm == 'FromZero': up = True
        else: up = c == 'hi' or (c == 'half' and R & 1 == 1)
    if up: R += 1
    if R == 0: return ('zero', neg), ix
    tz = (R & -R).bit_length() - 1
    return ('fin', neg, R >> tz, qe + tz), ix

TINY = {'sin', 'cos', 'tan', 'atan', 'sinh', 'cosh', 'tanh', 'atanh', 'exp', 'exp2', 'exp10', 'log1p'} | JAC
def tiny_fun(op, xs):
    if op in ('ns', 'cs', 'ds'): return 'pole'
    if op in JAC:
        m = xs[1]
        if op == 'am':
            return lambda t: mpmath.quad(lambda s: mpmath.ellipfun('dn', s, m), [0, t]) if t != 0 else mpf(0)
        return lambda t: mpmath.ellipfun(op, t, m)
    return {'sin': mpmath.sin, 'cos': mpmath.cos, 'tan': mpmath.tan, 'atan': mpmath.atan,
            'sinh': mpmath.sinh, 'cosh': mpmath.cosh, 'tanh': mpmath.tanh, 'atanh': mpmath.atanh,
            'exp': mpmath.exp, 'exp2': lambda t: mpmath.power(2, t), 'exp10': lambda t: mpmath.power(10, t),
            'log1p': mpmath.log1p}[op]

def cplx(op, p, rm, ints, a, exact=False):
    k, part = ints
    A, B, C, Dd = (dv(x) for x in a)
    if any(x.kind in ('nan', 'inf') for x in a): return None
    if op == 'cplx_as':
        s = dadd(A if part == 0 else B, (C if part == 0 else Dd) if k == 0 else dneg(C if part == 0 else Dd))
        return rnd(s, p, rm)
    if k == 0:  # mul: componentwise rounded products then rounded sum
        if rm == "None" or exact:  # None bounds come from the exact value, not the per-mode componentwise model
            x = dadd(dmul(A, C), dneg(dmul(B, Dd))) if part == 0 else dadd(dmul(A, Dd), dmul(B, C))
            return rnd(x, p, rm)
        def rr(d):
            t, _ = rnd(d, p, rm)
            if t[0] == 'inf': raise Ambig('cplx-inf')
            return D(t[1], t[2], t[3]) if t[0] == 'fin' else D(t[1], 0, 0)
        if part == 0: x = dadd(rr(dmul(A, C)), dneg(rr(dmul(B, Dd))))
        else: x = dadd(rr(dmul(A, Dd)), rr(dmul(B, C)))
        return rnd(x, p, rm)
    den = dadd(dmul(C, C), dmul(Dd, Dd))
    num = dadd(dmul(A, C), dmul(B, Dd)) if part == 0 else dadd(dmul(B, C), dneg(dmul(A, Dd)))
    if den.M == 0: return None
    if num.M == 0:
        if num.st: raise Ambig('cdiv-eps')
        return None
    return rnd(ddiv(num, den, p + 8), p, rm)

def cplx_md_none_interval(p, ints, a):
    """None-mode complex mul: each product may land anywhere between its ToZero and FromZero
    roundings (None keeps extra bits), then the sum is rounded: -> (lo, hi) tuples or None"""
    k, part = ints
    if k != 0: return None
    A, B, C, Dd = (dv(x) for x in a)
    def toD(t):
        if t[0] == 'inf': raise Ambig('cplx-inf')
        return D(t[1], t[2], t[3]) if t[0] == 'fin' else D(t[1], 0, 0)
    def pb(d):
        t1, _ = rnd(d, p, 'ToZero'); t2, _ = rnd(d, p, 'FromZero')
        return (toD(t1), toD(t2)) if tle(t1, t2) else (toD(t2), toD(t1))
    if part == 0:
        (l1, h1), (l2, h2) = pb(dmul(A, C)), pb(dmul(B, Dd))
        lo, hi = dadd(l1, dneg(h2)), dadd(h1, dneg(l2))
    else:
        (l1, h1), (l2, h2) = pb(dmul(A, Dd)), pb(dmul(B, C))
        lo, hi = dadd(l1, l2), dadd(h1, h2)
    return rnd(lo, p, 'Down')[0], rnd(hi, p, 'Up')[0]

def main():
    path = sys.argv[1]; maxex = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    stats = collections.Counter(); ex = collections.defaultdict(list); spec = collections.Counter()
    flags = collections.Counter()
    for ln in open(path):
        f = ln.rstrip('\n').split('\t')
        op, p, rm, ints = f[0], int(f[1]), f[2], eval(f[3])
        CUR_RM[0] = 'ToEven' if op.startswith('op') else rm
        i = f.index('=')
        args = []
        for s in f[4:i]:
            r, ixf = s.rsplit('|', 1); v = parse(r); args.append(v)
        res = parse(f[i + 1]); r_ix = f[i + 2] == 'true'; r_acc = f[i + 3] == 'true'
        key = op
        got = got_t(res)
        if any(x.kind in ('nan', 'inf') for x in args) or res.kind in ('nan', 'inf') and op not in ('exp','exp2','exp10','sinh','cosh','mul','sqr','op*','div','op/','ldexp','powi','add','sub','fma','hypot','pow','reciprocal','ns','cs','ds','nc','nd','dc','sc','sd','cd','tan'):
            spec[(op, rm if rm in ('Down',) else '', tuple(cls(x) for x in args), cls(res) if res.kind != 'nan' else res.M)] += 1
            stats[key, 'special'] += 1
            continue
        try:
            if op.startswith('cplx'):
                e = cplx(op, p, rm, ints, args)
                if e is None:
                    spec[(op, rm if rm == 'Down' else '', tuple(cls(x) for x in args), cls(res), str(ints))] += 1
                    stats[key, 'special'] += 1; continue
                exp_t, exp_ix = e
            else:
                exp_t, exp_ix = expected(op, p, rm, ints, args, got)
        except ZeroDivisionError:
            ok_z = res.kind in ("inf", "nan")
            stats[key, "divzero-ok" if ok_z else "BAD"] += 1
            if not ok_z and len(ex[key, "BAD"]) < maxex: ex[key, "BAD"].append(("divzero", ln[:600]))
            continue
        except Ambig as a:
            stats[key, 'ambig'] += 1
            if len(ex[key, 'ambig']) < maxex: ex[key, 'ambig'].append((str(a), ln[:600]))
            continue
        if rm == 'None' and exp_ix:
            # 'skip rounding' may keep extra bits: the result must lie between the ToZero and FromZero roundings at p
            try:
                lo_t, _ = cplx(op, p, "ToZero", ints, args, True) if op.startswith("cplx") else expected(op, p, "ToZero", ints, args, got)
                hi_t, _ = cplx(op, p, "FromZero", ints, args, True) if op.startswith("cplx") else expected(op, p, "FromZero", ints, args, got)
            except Ambig:
                stats[key, 'ambig'] += 1; continue
            fin = all(t[0] in ("fin", "zero") for t in (lo_t, hi_t, got))
            if not fin:
                ok_n = got in (lo_t, hi_t) or (hi_t[0] == "inf" and got[0] == "fin" and got[1] == hi_t[1] and (tle(lo_t, got) if not got[1] else tle(got, lo_t)))
            else:
                ok_n = (tle(lo_t, got) and tle(got, hi_t)) or (tle(hi_t, got) and tle(got, lo_t))
                if ok_n and not r_ix and not op.startswith("cplx"):
                    p2 = max(got[2].bit_length() + 2, EMIN - got[3] + 2) if got[0] == "fin" else 3
                    ok_n = expected(op, p2, "ToZero", ints, args, got)[0] == got == expected(op, p2, "FromZero", ints, args, got)[0]
            if not ok_n and op == "cplx_md" and fin:
                # complex mul in None mode may round the component products: accept the componentwise bounds too
                lo2, _ = cplx(op, p, "ToZero", ints, args); hi2, _ = cplx(op, p, "FromZero", ints, args)
                ok_n = all(t[0] in ("fin", "zero") for t in (lo2, hi2)) and ((tle(lo2, got) and tle(got, hi2)) or (tle(hi2, got) and tle(got, lo2)))
                if not ok_n:
                    iv = cplx_md_none_interval(p, ints, args)
                    ok_n = iv is not None and all(t[0] in ("fin", "zero") for t in iv) and tle(iv[0], got) and tle(got, iv[1])
            stats[key, 'none-ok' if ok_n else 'BAD'] += 1
            if not ok_n and len(ex[key, 'BAD']) < maxex: ex[key, 'BAD'].append(('none', ln[:900]))
            continue
        if exp_t == ('nan',):
            ok = res.kind == 'nan'
        elif exp_t[0] == 'zero' and exp_t[1] is None:
            ok = res.kind == 'zero'
        else:
            ok = exp_t == got
        if ok:
            stats[key, 'ok'] += 1
            if res.kind in ("inf", "nan"):
                pass  # Inf and NaN carry no flags
            elif exp_ix is not None and exp_ix != r_ix and not op.startswith("cplx") and "op_inexact: None" not in f[i + 1]:
                flags[key, 'inexact-flag'] += 1
                if len(ex[key, 'flag']) < maxex: ex[key, 'flag'].append((exp_ix, ln[:600]))
            acc_exp = r_ix or any(x.acc for x in args)
            if op in ("mul", "sqr", "op*", "mul_full", "div", "op/") and res.kind == "zero" and any(x.kind == "zero" for x in args):
                acc_exp = r_ix  # an exact zero factor makes an exact zero: no history
            if r_acc != acc_exp and res.kind not in ("inf", "nan") and not op.startswith("cplx"):
                flags[key, 'acc-flag'] += 1
                if len(ex[key, 'acc']) < maxex: ex[key, 'acc'].append((acc_exp, ln[:600]))
        else:
            stats[key, 'BAD'] += 1
            if len(ex[key, 'BAD']) < maxex: ex[key, 'BAD'].append((str(exp_t)[:300], ln[:900]))
    for k in sorted(stats): print('STAT', k, stats[k])
    for k in sorted(flags): print('FLAG', k, flags[k])
    for k in sorted(spec): print('SPEC', k, spec[k])
    for k in sorted(ex):
        for e in ex[k]: print('EX', k, e)
main()
