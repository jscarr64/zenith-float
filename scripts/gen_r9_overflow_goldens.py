#!/usr/bin/env python3
"""R9 goldens: overflow and underflow at the i64 exponent limits, every directed mode.

Writes zenith-float-num/tests/oracle/goldens/r9_overflow.jsonl in the canon of
gen_s9_specials_goldens.py (little-endian 64-bit limbs, value = 0.1m * 2^e,
exponent range +-(2^62 - 1), subnormal grid 2^(EMIN - p)).

Oracles:
  * "model" ops (add, sub, mul, sqr, fma, ldexp, powi, powsi, exp2/pow of exact powers
    of two): an exact integer model rounded once with s9.round_rat. Operands whose
    exponents are far apart are folded into a sticky unit (see `fold`) so no huge
    shifts happen. Each model row is also checked against mpmath at p+256 and p+512:
    when both error intervals round to one result it must equal the model.
  * transcendental ops (exp, exp10, sinh, cosh, pow, parse): mpmath at p+256 and
    p+512 (Ziv, escalating to p+16384) whose rounded intervals must agree.
Rows stream to disk.  Run with ~/.local/oracle-venv/bin/python (mpmath 1.3.0).
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("s9", HERE / "gen_s9_specials_goldens.py")
s9 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s9)

EMAX, EMIN, MODES = s9.EMAX, s9.EMIN, s9.MODES
round_rat, zero, Arg, interval_round = s9.round_rat, s9.zero, s9.Arg, s9.interval_round
OUT = s9.ROOT / "zenith-float-num" / "tests" / "oracle" / "goldens" / "r9_overflow.jsonl"
EXTRAS = [256, 512, 1024, 2048, 4096, 8192, 16384]
PS = [24, 53, 64, 113, 256, 1024]


def A(sign, n, ez):
    return Arg(sign, n, ez)


def maxv(p, sign="+"):
    return A(sign, (1 << p) - 1, EMAX - p)


def minsub(p, sign="+"):
    return A(sign, 1, EMIN - p)


def minnorm(sign="+"):
    return A(sign, 1, EMIN - 1)


def enc_exact(a: Arg, pa: int) -> dict:
    d = a.enc(pa)
    if d["k"] != "f":
        raise ValueError(f"arg {a.sign}{a.n}*2^{a.ez} not finite at {pa}")
    if a.n != 0:
        back = round_rat(a.sign, a.n, 1, a.ez, pa, "Up")
        if back != d:
            raise ValueError(f"arg {a.sign}{a.n}*2^{a.ez} not exact at {pa} bits")
    return d


# ------------------------------------------------------------ exact model

def top(a):
    return a[2] + a[1].bit_length()


def fold(big, small, guard):
    """Replace `small` by a sticky unit far below big's last bit (no breakpoint between)."""
    k = guard
    s = 1 if small[0] == big[0] else -1
    return (big[0], (big[1] << k) + s, big[2] - k)


def add_model(a, b, rm, guard):
    """a, b: (sign, n, ez) exact; returns (sign, n, ez) of a+b or ('0', zero sign)."""
    if a[1] == 0 and b[1] == 0:
        if a[0] == b[0]:
            return ("z", a[0])
        return ("z", "-" if rm == "Down" else "+")
    if a[1] == 0:
        return b
    if b[1] == 0:
        return a
    if top(a) - top(b) > 2 * guard:
        return fold(a, b, guard)
    if top(b) - top(a) > 2 * guard:
        return fold(b, a, guard)
    ez = min(a[2], b[2])
    va = (a[1] << (a[2] - ez)) * (1 if a[0] == "+" else -1)
    vb = (b[1] << (b[2] - ez)) * (1 if b[0] == "+" else -1)
    v = va + vb
    if v == 0:
        return ("z", "-" if rm == "Down" else "+")
    return ("-" if v < 0 else "+", abs(v), ez)


def tri(a: Arg):
    return (a.sign, a.n, a.ez)


def neg(a: Arg):
    return Arg("-" if a.sign == "+" else "+", a.n, a.ez)


def rnd(t, p, rm):
    if t[0] == "z":
        return zero(t[1], p)
    return round_rat(t[0], t[1], 1, t[2], p, rm)


def xsign(a, b):
    return "+" if a.sign == b.sign else "-"


def model(fn, args, k, p, rm, guard, s=None):
    if fn == "add":
        return rnd(add_model(tri(args[0]), tri(args[1]), rm, guard), p, rm)
    if fn == "sub":
        return rnd(add_model(tri(args[0]), tri(neg(args[1])), rm, guard), p, rm)
    if fn in ("mul", "sqr"):
        a = args[0]
        b = args[1] if fn == "mul" else args[0]
        if a.n == 0 or b.n == 0:
            return zero(xsign(a, b), p)
        return round_rat(xsign(a, b), a.n * b.n, 1, a.ez + b.ez, p, rm)
    if fn == "fma":
        a, b, c = args
        prod = (xsign(a, b), a.n * b.n, a.ez + b.ez)
        return rnd(add_model(prod, tri(c), rm, guard), p, rm)
    if fn == "ldexp":
        a = args[0]
        if a.n == 0:
            return zero(a.sign, p)
        return round_rat(a.sign, a.n, 1, a.ez + k, p, rm)
    if fn == "powi":
        a = args[0]
        s = "-" if (a.sign == "-" and k % 2 == 1) else "+"
        return round_rat(s, a.n ** k, 1, a.ez * k, p, rm)
    if fn == "powsi":  # k < 0
        a = args[0]
        m = -k
        s = "-" if (a.sign == "-" and m % 2 == 1) else "+"
        return round_rat(s, 1, a.n ** m, -a.ez * m, p, rm)
    if fn == "exp2":  # integer argument only
        a = args[0]
        return round_rat("+", 1, 1, (a.n << a.ez) * (1 if a.sign == "+" else -1), p, rm)
    if fn == "pow":  # x = 2^j, y integer
        x, y = args
        j = x.ez + x.n.bit_length() - 1
        yi = (y.n << y.ez) * (1 if y.sign == "+" else -1)
        return round_rat("+", 1, 1, j * yi, p, rm)
    if fn == "parse":
        txt = s
        sign = "-" if txt.startswith("-") else "+"
        txt = txt.lstrip("+-")
        mant, _, ex = txt.partition("e")
        ip, _, fp = mant.partition(".")
        digits = int(ip + fp)
        dexp = int(ex or "0") - len(fp)
        if digits == 0:
            return zero(sign, p)
        if dexp >= 0:
            return round_rat(sign, digits * 10 ** dexp, 1, 0, p, rm)
        return round_rat(sign, digits, 10 ** -dexp, 0, p, rm)
    raise ValueError(fn)


# ------------------------------------------------------------- mpmath

def mpf_of(fn, args, k, s):
    xs = [a.mpf() for a in args]
    if fn == "add":
        return mp.fadd(xs[0], xs[1])
    if fn == "sub":
        return mp.fsub(xs[0], xs[1])
    if fn == "mul":
        return mp.fmul(xs[0], xs[1])
    if fn == "sqr":
        return mp.fmul(xs[0], xs[0])
    if fn == "fma":
        return mp.fadd(mp.fmul(xs[0], xs[1], exact=True), xs[2])
    if fn == "ldexp":
        return mp.ldexp(xs[0], k)
    if fn == "powi":
        return mp.power(xs[0], k)
    if fn == "powsi":
        return mp.power(xs[0], k)
    if fn == "exp":
        return mp.exp(xs[0])
    if fn == "exp2":
        return mp.power(2, xs[0])
    if fn == "exp10":
        return mp.power(10, xs[0])
    if fn == "sinh":
        return mp.sinh(xs[0])
    if fn == "cosh":
        return mp.cosh(xs[0])
    if fn == "pow":
        return mp.power(xs[0], xs[1])
    if fn == "parse":
        return mp.mpf(s)
    raise ValueError(fn)


def ziv(fn, args, k, s, p, rm, cache):
    for i in range(len(EXTRAS) - 1):
        got = []
        for extra in (EXTRAS[i], EXTRAS[i + 1]):
            if extra not in cache:
                with mp.workprec(p + extra):
                    cache[extra] = mpf_of(fn, args, k, s)
            got.append(interval_round(cache[extra], p + extra, p, rm))
        if got[0] is not None and got[0] == got[1]:
            return got[0], EXTRAS[i]
        if fn in MODEL_FNS and i >= 1:
            return None
    return None


MODEL_FNS = {"add", "sub", "mul", "sqr", "fma", "ldexp", "powi", "powsi"}


# --------------------------------------------------------------- cases

CASES: list[dict] = []


def case(fn, mk, k=0, s=None, ps=None, pa=lambda p: p, exact=False, note=""):
    CASES.append({"fn": fn, "mk": mk, "k": k, "s": s, "ps": ps or PS, "pa": pa,
                  "exact": exact, "note": note})


def plus64(p):
    return p + 64


def build():
    for sg in "+-":
        ng = "-" if sg == "+" else "+"
        # ---- add/sub at the overflow edge
        case("add", lambda p, sg=sg: [maxv(p, sg), A(sg, 1, EMAX - p - 1)], note="MAX+half ulp (tie)")
        case("add", lambda p, sg=sg: [maxv(p, sg), A(sg, 1, EMAX - p - 2)], note="MAX+quarter ulp")
        case("add", lambda p, sg=sg: [maxv(p, sg), A(sg, 3, EMAX - p - 2)], note="MAX+3/4 ulp")
        case("add", lambda p, sg=sg: [maxv(p, sg), A(sg, 1, EMAX - p)], note="MAX+1 ulp")
        case("add", lambda p, sg=sg: [maxv(p, sg), maxv(p, sg)], note="MAX+MAX far past")
        case("add", lambda p, sg=sg: [maxv(p, sg), minsub(p, sg)], note="MAX+min subnormal")
        case("add", lambda p, sg=sg: [maxv(p, sg), A(sg, 1, EMAX - p - 1)], pa=plus64, note="MAX+half ulp wide args")
        case("add", lambda p, sg=sg: [maxv(p + 64, sg), A(sg, 1, EMAX - p - 70)], pa=plus64,
             note="wide MAX rounds up past MAX")
        case("sub", lambda p, sg=sg: [maxv(p, sg), A(sg, 1, EMAX - p - 1)], note="MAX-half ulp")
        case("sub", lambda p, sg=sg, ng=ng: [maxv(p, sg), A(ng, 1, EMAX - p)], note="MAX-(-ulp)")
        case("sub", lambda p, sg=sg: [maxv(p, sg), maxv(p, sg)], note="exact cancel at MAX")
        # ---- add/sub at the underflow edge
        case("sub", lambda p, sg=sg: [minnorm(sg), minsub(p, sg)], note="min normal - 1 subnormal ulp")
        case("sub", lambda p, sg=sg: [minnorm(sg), A(sg, 1, EMIN - p - 1)], pa=plus64,
             note="min normal - half subnormal ulp")
        case("sub", lambda p, sg=sg: [minnorm(sg), A(sg, 1, EMIN - p - 40)], pa=plus64,
             note="min normal - tiny")
        case("sub", lambda p, sg=sg: [A(sg, 3, EMIN - p - 1), A(sg, 1, EMIN - p - 1)], pa=plus64,
             note="exact one subnormal ulp")
        case("add", lambda p, sg=sg, ng=ng: [A(sg, 1, EMIN - p), A(ng, 1, EMIN - p - 1)], pa=plus64,
             note="half min subnormal (tie to 0)")
        case("add", lambda p, sg=sg, ng=ng: [A(sg, 1, EMIN - p), A(ng, 3, EMIN - p - 2)], pa=plus64,
             note="quarter min subnormal")
        case("add", lambda p, sg=sg, ng=ng: [A(sg, 1, EMIN - p - 20), A(ng, 1, EMIN - p - 30)],
             pa=lambda p: p + 64, note="far below min subnormal")
        case("add", lambda p, sg=sg, ng=ng: [A(sg, 3, EMIN - p - 1), A(ng, 1, EMIN - p)], pa=plus64,
             note="1.5-1 subnormal ulps")
        case("add", lambda p, sg=sg, ng=ng: [A(sg, 5, EMIN - p - 1), A(ng, 1, EMIN - p)], pa=plus64,
             note="tie 1.5 subnormal ulps -> 2")
        case("add", lambda p, sg=sg: [minsub(p, sg), minsub(p, sg)], note="2 subnormal ulps exact")
        case("add", lambda p, sg=sg: [A(sg, (1 << p) - 1, EMIN - p), minsub(p, sg)],
             note="largest subnormal + ulp = min normal")
        case("add", lambda p, sg=sg, ng=ng: [minsub(p, sg), minsub(p, ng)], note="exact cancel at min subnormal")
        case("add", lambda p, sg=sg, ng=ng: [A(sg, (1 << (p + 63)) + 1, EMIN - p - 63),
                                             A(ng, 1, EMIN - 1)], pa=plus64,
             note="double-rounding probe: below min normal + 2^-63 ulp")
        # ---- mul/sqr
        case("mul", lambda p, sg=sg: [maxv(p, sg), A("+", 1, 1)], note="MAX*2")
        case("mul", lambda p, sg=sg: [maxv(p, sg), A("+", 3, -1)], note="MAX*1.5")
        case("mul", lambda p, sg=sg: [maxv(p, sg), A("+", (1 << p) - 1, -p)], note="MAX*(1-ulp)")
        case("mul", lambda p, sg=sg: [A(sg, 1, (1 << 61)), A("+", 1, (1 << 61) - 1)], note="2^(2^62) exactly past")
        case("mul", lambda p, sg=sg: [A(sg, 1, (1 << 61)), A("+", 1, (1 << 61) - 2)], note="2^(EMAX-1) at limit")
        case("mul", lambda p, sg=sg: [A(sg, 3, (1 << 61) - 1), A("+", 1, (1 << 61) - 2)], note="1.5*2^(EMAX-1)")
        case("mul", lambda p, sg=sg: [A(sg, 1, -(1 << 61)), A("+", 1, -(1 << 61) + 1)], note="min normal exactly")
        case("mul", lambda p, sg=sg: [A(sg, 1, -(1 << 61)), A("+", 1, -(1 << 61) + 1 - p)], note="min subnormal exactly")
        case("mul", lambda p, sg=sg: [A(sg, 1, -(1 << 61)), A("+", 1, -(1 << 61) - p)], note="half min subnormal")
        case("mul", lambda p, sg=sg: [A(sg, 3, -(1 << 61) - 1), A("+", 1, -(1 << 61) - p)], note="3/4 min subnormal")
        case("mul", lambda p, sg=sg: [A(sg, 1, -(1 << 61)), A("+", 1, -(1 << 61) - p - 50)], note="far below")
        case("mul", lambda p, sg=sg: [minsub(p, sg), A("+", 3, -1)], note="1.5 subnormal ulps tie")
        case("mul", lambda p, sg=sg: [A(sg, 3, EMIN - p), A("+", 1, -1)], note="1.5 subnormal ulps tie to 2")
        case("mul", lambda p, sg=sg: [A(sg, (1 << (p - 1)) + 1, EMIN + 10 - p),
                                      A("+", (1 << (p - 1)) + 1, -13 - p + 1)],
             note="double-rounding probe (p-bit round then grid would tie-even)")
        case("mul", lambda p, sg=sg: [A(sg, (1 << p) - 1, EMIN + 5 - p), A("+", (1 << p) - 1, -9 - p)],
             note="subnormal product rounding")
        case("sqr", lambda p, sg=sg: [A(sg, 1, (1 << 61) - 1)], note="(2^(2^61))^2 = 2^(2^62) past")
        case("sqr", lambda p, sg=sg: [A(sg, 3, (1 << 61) - 3)], note="sqr just under")
        case("sqr", lambda p, sg=sg: [A(sg, 1, -(1 << 61))], note="sqr below min normal")
        case("sqr", lambda p, sg=sg: [A(sg, 1, -(1 << 61) - p // 2)], note="sqr into subnormal")
        case("sqr", lambda p, sg=sg: [maxv(p, sg)], note="sqr MAX far past")
        case("sqr", lambda p, sg=sg: [minsub(p, sg)], note="sqr min subnormal far below")
        # ---- fma
        case("fma", lambda p, sg=sg: [maxv(p, sg), A("+", 1, 0), A(sg, 1, EMAX - p - 1)], note="MAX*1/2*2 + half ulp")
        case("fma", lambda p, sg=sg: [A(sg, 1, (1 << 61)), A("+", 1, (1 << 61) - 1), maxv(p, "-" if sg == "+" else "+")],
             note="2^(2^62) - MAX exact intermediate overflow")
        case("fma", lambda p, sg=sg: [maxv(p, sg), A("+", 1, 1), maxv(p, "-" if sg == "+" else "+")],
             note="2MAX - MAX = MAX")
        case("fma", lambda p, sg=sg: [A(sg, 1, -(1 << 61)), A("+", 1, -(1 << 61) - p), A(sg, 1, EMIN - p)],
             note="half min sub + min sub tie")
        case("fma", lambda p, sg=sg, ng=ng: [A(sg, 1, -(1 << 61)), A("+", 1, -(1 << 61) + 1), A(ng, 1, EMIN - 1)],
             note="exact cancel at min normal")
        case("fma", lambda p, sg=sg, ng=ng: [minsub(p, sg), A("+", 1, -1), A(ng, 1, EMIN - p - 40)],
             pa=plus64, note="half min sub - tiny")
        case("fma", lambda p, sg=sg: [minsub(p, sg), minsub(p, "+"), A(sg, 1, EMIN - p)],
             note="tiny product + min sub")
        # ---- ldexp
        case("ldexp", lambda p, sg=sg: [maxv(p, sg)], k=1, note="ldexp MAX,1")
        case("ldexp", lambda p, sg=sg: [A(sg, 1, 0)], k=EMAX - 1, note="ldexp 1,EMAX-1 at limit")
        case("ldexp", lambda p, sg=sg: [A(sg, 1, 0)], k=EMAX, note="ldexp 1,EMAX past")
        case("ldexp", lambda p, sg=sg: [A(sg, 1, 0)], k=(1 << 63) - 1, note="ldexp i64::MAX")
        case("ldexp", lambda p, sg=sg: [A(sg, 1, 0)], k=-(1 << 63), note="ldexp i64::MIN")
        case("ldexp", lambda p, sg=sg: [maxv(p, sg)], k=-(1 << 63), note="ldexp MAX, i64::MIN")
        case("ldexp", lambda p, sg=sg: [A(sg, 1, -1)], k=EMIN, note="ldexp 1/2,EMIN = min normal")
        case("ldexp", lambda p, sg=sg: [A(sg, 1, 0)], k=lambda p: EMIN - p, note="min subnormal")
        case("ldexp", lambda p, sg=sg: [A(sg, 1, 0)], k=lambda p: EMIN - p - 1, note="half min subnormal")
        case("ldexp", lambda p, sg=sg: [A(sg, 3, 0)], k=lambda p: EMIN - p - 2, note="3/4 min subnormal")
        case("ldexp", lambda p, sg=sg: [A(sg, 3, 0)], k=lambda p: EMIN - p - 1, note="1.5 min subnormal")
        case("ldexp", lambda p, sg=sg: [maxv(p, sg)], k=EMIN - EMAX + 3, note="MAX into subnormal rounding")
        case("ldexp", lambda p, sg=sg: [A(sg, (1 << p) - 1, -p)], k=EMIN - 2, note="(1-ulp)*2^(EMIN-2)")
        case("ldexp", lambda p, sg=sg: [minsub(p, sg)], k=EMAX, note="min sub scaled up to normal")
        case("ldexp", lambda p, sg=sg: [minsub(p, sg)], k=2 * EMAX, note="min sub scaled by 2*EMAX")
        # ---- powi / powsi
        case("powi", lambda p, sg=sg: [A(sg, 1, (1 << 61) - 1)], k=2, note="powi 2 past")
        case("powi", lambda p, sg=sg: [A(sg, 3, (1 << 60) - 2)], k=3, note="powi 3 past")
        case("powi", lambda p, sg=sg: [A(sg, 3, 1537228672809129300)], k=3, note="powi 3 near")
        case("powi", lambda p, sg=sg: [A(sg, 3, -(1 << 60) - 2)], k=3, note="powi 3 subnormal")
        case("powi", lambda p, sg=sg: [A(sg, 1, -(1 << 61) - 1)], k=2, note="powi 2 far below")
        case("powi", lambda p, sg=sg: [A(sg, 3, -1537228672809129302)], k=3, note="powi 3 near underflow")
        case("powsi", lambda p, sg=sg: [A(sg, 1, -(1 << 61))], k=-2, note="powsi -2 past")
        case("powsi", lambda p, sg=sg: [A(sg, 1, (1 << 61))], k=-2, note="powsi -2 below min normal")
        case("powsi", lambda p, sg=sg: [A(sg, 3, (1 << 60) - 3)], k=-3, note="powsi -3 underflow rounding")
        case("powsi", lambda p, sg=sg: [minsub(p, sg)], k=-1, note="1/min subnormal past")
        case("powsi", lambda p, sg=sg: [maxv(p, sg)], k=-1, note="1/MAX subnormal")
    # ---- exp family (positive results)
    T = 3196577161300663808  # ~ EMAX*ln2
    case("exp", lambda p: [A("+", T, 0)], note="exp near overflow")
    case("exp", lambda p: [A("+", T + (1 << 40), 0)], note="exp past")
    case("exp", lambda p: [A("+", 3196577161300663915, 0)], note="exp just past threshold")
    case("exp", lambda p: [A("+", 3196577161300663914, 0)], note="exp just under threshold")
    case("exp", lambda p: [A("+", 1, 70)], note="exp far past")
    case("exp", lambda p: [A("-", 3196577161300663914, 0)], note="exp near min normal")
    case("exp", lambda p: [A("-", 3196577161300663915 + 20, 0)], note="exp into subnormal")
    case("exp", lambda p: [A("-", 3196577161300663915 + 2 * p, 0)], note="exp below subnormal")
    case("exp", lambda p: [A("-", 1, 70)], note="exp far below")
    for k in (EMAX - 1, EMAX, EMAX + 1, EMIN - 1, EMIN - 2, EMIN - 30, -(1 << 63)):
        case("exp2", lambda p, k=k: [A("+" if k > 0 else "-", abs(k), 0)], exact=True, note=f"exp2 {k}")
    for kk in (lambda p: EMIN - p, lambda p: EMIN - p - 1, lambda p: EMIN - p - 2):
        case("exp2", lambda p, kk=kk: [A("-", -kk(p), 0)], exact=True, note="exp2 subnormal grid")
    case("exp2", lambda p: [A("+", 2 * EMAX - 1, -1)], note="exp2 EMAX-1/2")
    case("exp2", lambda p: [A("+", 2 * EMAX + 1, -1)], note="exp2 EMAX+1/2")
    case("exp2", lambda p: [A("-", 2 * (-(EMIN - p)) - 1, -1)], note="exp2 EMIN-p+1/2")
    case("exp2", lambda p: [A("-", 2 * (-(EMIN - p)) + 1, -1)], note="exp2 EMIN-p-1/2 (below min sub, above half)")
    case("exp2", lambda p: [A("-", 2 * (-(EMIN - p - 1)) + 1, -1)], note="exp2 EMIN-p-3/2 (below half)")
    D = 1388255822130839282  # floor(EMAX*log10 2)
    case("exp10", lambda p: [A("+", D, 0)], note="exp10 under")
    case("exp10", lambda p: [A("+", D + 1, 0)], note="exp10 past")
    case("exp10", lambda p: [A("-", D, 0)], note="exp10 near min normal")
    case("exp10", lambda p: [A("-", D + 2, 0)], note="exp10 below")
    case("exp10", lambda p: [A("-", D + 1 + (p * 3) // 10, 0)], note="exp10 subnormal")
    for sg in "+-":
        case("sinh", lambda p, sg=sg: [A(sg, 3196577161300663914 + 1, 0)], note="sinh near overflow")
        case("sinh", lambda p, sg=sg: [A(sg, 3196577161300663915 + 1, 0)], note="sinh past")
        case("sinh", lambda p, sg=sg: [A(sg, 1, 70)], note="sinh far past")
        case("cosh", lambda p, sg=sg: [A(sg, 3196577161300663915 + 1, 0)], note="cosh past")
        case("cosh", lambda p, sg=sg: [A(sg, 3196577161300663914 + 1, 0)], note="cosh near")
    # ---- pow
    case("pow", lambda p: [A("+", 1, 1), A("+", EMAX - 1, 0)], exact=True, note="2^(EMAX-1)")
    case("pow", lambda p: [A("+", 1, 1), A("+", EMAX, 0)], exact=True, note="2^EMAX past")
    case("pow", lambda p: [A("+", 1, 1), A("-", -(EMIN - p), 0)], exact=True, note="2^(EMIN-p)")
    case("pow", lambda p: [A("+", 1, 1), A("-", -(EMIN - p - 1), 0)], exact=True, note="2^(EMIN-p-1) tie")
    case("pow", lambda p: [A("+", 1, 1), A("-", -(EMIN - p - 5), 0)], exact=True, note="2^(EMIN-p-5)")
    case("pow", lambda p: [A("+", 3, 0), A("+", 2909155196620405700, 0)], note="3^y near overflow")
    case("pow", lambda p: [A("+", 3, 0), A("+", 2909155196620405800, 0)], note="3^y past")
    case("pow", lambda p: [A("+", 3, 0), A("-", 2909155196620405800, 0)], note="3^-y below")
    case("pow", lambda p: [A("+", 10, 0), A("+", 1, 70)], note="10^huge far past")
    case("pow", lambda p: [A("+", 3, -1), A("+", 7, 59)], note="1.5^(7*2^59)")
    # ---- parse (decimal), incl cosmology
    for s in [f"1e{D}", f"1e{D + 1}", f"9.9e{D}", f"1.1897e{D}", f"1.1898e{D}", f"5e{D}",
              f"-1e{D + 1}", f"-1.1898e{D}", f"1e{D + 100}", f"-1e{D + 100}",
              f"1e-{D}", f"1e-{D + 1}", f"-1e-{D + 1}", f"1e-{D + 30}", f"-1e-{D + 30}",
              f"1e-{D + 400}", f"-1e-{D + 400}",
              "6.02214076e23", "1.616255e-35", "1e120", "-1e-120", "1.9884e30", "4.35e17",
              "1e100000", "1e-100000"]:
        small = abs(int(s.partition("e")[2] or "0")) <= 200000
        case("parse", lambda p: [], s=s, exact=small, note=f"parse {s}")


def main():
    build()
    t0 = time.time()
    rows = mism = 0
    hows: dict = {}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w") as fh:
        for c in CASES:
            fn = c["fn"]
            for p in c["ps"]:
                pa = c["pa"](p) if fn in MODEL_FNS else max(c["pa"](p), 128)
                try:
                    args = c["mk"](p)
                    encs = [enc_exact(a, pa) for a in args]
                except ValueError as e:
                    print("skip", fn, c["note"], p, e, file=sys.stderr)
                    continue
                kv = c["k"](p) if callable(c["k"]) else c["k"]
                cache: dict = {}
                guard = 2 * (pa + p) + 16
                for rm in MODES:
                    if fn in MODEL_FNS or c["exact"]:
                        want = model(fn, args, kv, p, rm, guard, c["s"])
                        how = "model"
                        if fn in MODEL_FNS or fn == "parse":
                            z = ziv(fn, args, kv, c["s"], p, rm, cache)
                            if z is not None:
                                if z[0] != want:
                                    mism += 1
                                    print("MISMATCH", fn, c["note"], p, rm, z[0], want, file=sys.stderr)
                                how = "model+ziv"
                    else:
                        z = ziv(fn, args, kv, c["s"], p, rm, cache)
                        if z is None:
                            print("UNRESOLVED", fn, c["note"], p, rm, file=sys.stderr)
                            mism += 1
                            continue
                        want, how = z[0], "ziv"
                    hows[how] = hows.get(how, 0) + 1
                    row = {"fn": fn, "p": p, "rm": rm, "args": encs, "expect": want,
                           "how": how, "note": c["note"]}
                    if fn in ("ldexp", "powi", "powsi"):
                        row["n"] = kv
                    if fn == "parse":
                        row["s"] = c["s"]
                    fh.write(json.dumps(row) + "\n")
                    rows += 1
            fh.flush()
    print(f"rows={rows} hows={hows} mismatches={mism} secs={time.time() - t0:.1f}")
    if mism:
        sys.exit(1)


if __name__ == "__main__":
    main()
