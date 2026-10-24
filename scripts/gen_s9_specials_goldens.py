#!/usr/bin/env python3
"""S9 directed special-function goldens (R7 + R7b).

Writes zenith-float-num/tests/oracle/goldens/s9_specials.jsonl.

Functions: gamma, ln_gamma, erf, erfc, expm1, log1p (R7); digamma, J0, J1,
zeta, beta (R7b). Modes: ToEven, Up, Down, ToZero, FromZero.

Every value is rounded at exactly p bits onto the i64 exponent range
EXPONENT_MIN/MAX = -/+(2^62 - 1) (canon of gen_oracle_goldens.py: little-endian
64-bit limbs, value = 0.1m * 2^e, round_p(p) limbs, subnormal grid 2^(EMIN-p)).

Ziv check: each value is computed at p+256 and at p+512 bits. The error
interval of each approximation must round to one result in the requested
mode, and the two results must agree. Otherwise the extra precision doubles
(1024 .. 16384) and the row is marked hard. Rows whose value sits within
2^-(p+16384) of a p-bit grid point carry a documented residual sign
("resid": the sign of f - a for the p-bit point a); they are marked hard too.
Exact rows (gamma(n), zeta(-n), beta of integers, ...) use a rational oracle.

Run with ~/.local/oracle-venv/bin/python (mpmath 1.3.0, sympy 1.14.0).
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from math import factorial
from pathlib import Path

import mpmath as mp
import sympy as sp

EMAX = (1 << 62) - 1
EMIN = -EMAX
WORD = 64
MODES = ["ToEven", "Up", "Down", "ToZero", "FromZero"]
EXTRAS = [256, 512, 1024, 2048, 4096, 8192, 16384]
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "zenith-float-num" / "tests" / "oracle" / "goldens" / "s9_specials.jsonl"
HARD_OUT = ROOT / "zenith-float-num" / "tests" / "oracle" / "goldens" / "s9_specials_hard.txt"


def round_p(p: int) -> int:
    return ((p + WORD - 1) // WORD) * WORD


def words(mant: int, prec: int) -> list[str]:
    mask = (1 << WORD) - 1
    return [f"{(mant >> (WORD * i)) & mask:016x}" for i in range(prec // WORD)]


def zero(sign: str, p: int) -> dict:
    return {"k": "f", "s": sign, "e": 0, "m": words(0, round_p(p))}


def nan(err: str = "InvalidArgument") -> dict:
    return {"k": "nan", "err": err}


def inf(sign: str) -> dict:
    return {"k": "inf", "s": sign}


def away(rm: str, positive: bool, rem_cmp: int, q_odd: bool) -> bool:
    """rem_cmp: sign of (remainder - half ulp); called only for a nonzero remainder."""
    if rm == "ToZero":
        return False
    if rm == "FromZero":
        return True
    if rm == "Up":
        return positive
    if rm == "Down":
        return not positive
    if rem_cmp > 0:
        return True
    if rem_cmp < 0:
        return False
    return q_odd


def overflow(sign: str, p: int, rm: str) -> dict:
    toward_zero = rm == "ToZero" or (rm == "Down" and sign == "+") or (rm == "Up" and sign == "-")
    if toward_zero:
        mant = ((1 << p) - 1) << (round_p(p) - p)
        return {"k": "f", "s": sign, "e": EMAX, "m": words(mant, round_p(p))}
    return inf(sign)


def floor_log2(n: int, d: int) -> int:
    fl = n.bit_length() - d.bit_length()
    if fl >= 0:
        if n < (d << fl):
            fl -= 1
    elif (n << -fl) < d:
        fl -= 1
    return fl


def div_round(n: int, d: int, s: int, rm: str, positive: bool) -> int:
    """round(n/d * 2^s) to an integer in mode rm."""
    if s >= 0:
        num, den = n << s, d
    else:
        num, den = n, d << -s
    q, r = divmod(num, den)
    if r == 0:
        return q
    two_r = 2 * r
    cmp = (two_r > den) - (two_r < den)
    return q + 1 if away(rm, positive, cmp, q & 1 == 1) else q


def round_rat(sign: str, n: int, d: int, ez: int, p: int, rm: str) -> dict:
    """Round sign * n/d * 2^ez (n, d > 0) at exactly p bits."""
    positive = sign == "+"
    prec = round_p(p)
    fl = floor_log2(n, d) + ez
    e = fl + 1
    if e > EMAX:
        return overflow(sign, p, rm)
    if e >= EMIN:
        q = div_round(n, d, ez + p - e, rm, positive)
        if q == (1 << p):
            q >>= 1
            e += 1
            if e > EMAX:
                return overflow(sign, p, rm)
        return {"k": "f", "s": sign, "e": e, "m": words(q << (prec - p), prec)}
    # subnormal grid 2^(EMIN - p), e = EMIN
    if EMIN - p - fl > 4 * p + 8:
        q = 1 if away(rm, positive, -1, False) else 0  # far below half the grid step
    else:
        q = div_round(n, d, ez - (EMIN - p), rm, positive)
    if q == 0:
        return zero(sign, p)
    e = EMIN
    if q == (1 << p):
        q >>= 1
        e += 1
    return {"k": "f", "s": sign, "e": e, "m": words(q << (prec - p), prec)}


def round_fraction(v: Fraction, p: int, rm: str, zero_sign: str = "+") -> dict:
    if v == 0:
        return zero(zero_sign, p)
    sign = "-" if v < 0 else "+"
    return round_rat(sign, abs(v.numerator), v.denominator, 0, p, rm)


def mpf_parts(x) -> tuple[str, int, int]:
    s, man, exp, _bc = x._mpf_
    return ("-" if s else "+", int(man), int(exp))


# ---------------------------------------------------------------- arguments

class Arg:
    """An exact dyadic argument: sign * n * 2^ez, or a signed zero."""

    def __init__(self, sign: str, n: int, ez: int):
        self.sign, self.n, self.ez = sign, n, ez

    def mpf(self):
        if self.n == 0:
            return mp.mpf(0) if self.sign == "+" else -mp.mpf(0)
        return mp.mpf((-self.n if self.sign == "-" else self.n, self.ez))

    def sym(self):
        v = sp.Integer(self.n) * sp.Integer(2) ** self.ez
        return -v if self.sign == "-" else v

    def frac(self) -> Fraction:
        v = Fraction(self.n) * (Fraction(2) ** self.ez)
        return -v if self.sign == "-" else v

    def enc(self, p: int) -> dict:
        if self.n == 0:
            return zero(self.sign, p)
        d = round_rat(self.sign, self.n, 1, self.ez, p, "ToZero")
        return d

    def tiny_bits(self) -> int:
        return self.n.bit_length()


def arg_at(spec, p: int) -> Arg:
    """Exact p-bit argument from a spec.

    ("d", "1.5")         decimal, rounded ToEven at p bits
    ("pow2", k)          2^k, optional sign via ("pow2", k, "-")
    ("off", base, k)     base + 2^k with base a decimal string, rounded at p bits
    ("mp", fn)           mpmath constant fn() rounded at p bits
    ("z", sign)          signed zero
    """
    kind = spec[0]
    if kind == "z":
        return Arg(spec[1], 0, 0)
    if kind == "pow2":
        return Arg(spec[2] if len(spec) > 2 else "+", 1, spec[1])
    with mp.workprec(round_p(p) + 128):
        if kind == "d":
            v = mp.mpf(spec[1])
        elif kind == "off":
            v = mp.mpf(spec[1])
            with mp.workprec(round_p(p) + 4 * p + 256):
                v = mp.mpf(spec[1]) + mp.ldexp(1, spec[2])
        elif kind == "mp":
            v = spec[1]()
        else:
            raise ValueError(spec)
    if v == 0:
        return Arg("+", 0, 0)
    s, man, exp = mpf_parts(v)
    r = round_rat(s, man, 1, exp, p, "ToEven")
    limbs = r["m"]
    mant = 0
    for i, w in enumerate(limbs):
        mant |= int(w, 16) << (WORD * i)
    return Arg(s, mant, r["e"] - round_p(p))


# ---------------------------------------------------------------- functions

def J0(x):
    return mp.besselj(0, x)


def J1(x):
    return mp.besselj(1, x)


MPF = {
    "gamma": lambda a: mp.gamma(a[0]),
    # ln|Gamma(x)| (ruling 2): the real part of the principal loggamma.
    "ln_gamma": lambda a: mp.re(mp.loggamma(a[0])),
    "erf": lambda a: mp.erf(a[0]),
    "erfc": lambda a: mp.erfc(a[0]),
    "expm1": lambda a: mp.expm1(a[0]),
    "log1p": lambda a: mp.log1p(a[0]),
    # psi(1/2) = -euler - 2 ln 2: mpmath digamma is very slow at 16k+ bits.
    "digamma": lambda a: (-mp.euler - 2 * mp.ln2) if a[0] == mp.mpf(0.5) else mp.digamma(a[0]),
    "bessel_j0": lambda a: J0(a[0]),
    "bessel_j1": lambda a: J1(a[0]),
    "zeta": lambda a: mp.zeta(a[0]),
    "beta": lambda a: mp.beta(a[0], a[1]),
}

SYM = {
    "gamma": lambda a: sp.gamma(a[0]),
    "ln_gamma": lambda a: sp.log(sp.Abs(sp.gamma(a[0]))),
    "erf": lambda a: sp.erf(a[0]),
    "erfc": lambda a: sp.erfc(a[0]),
    "expm1": lambda a: sp.exp(a[0]) - 1,
    "log1p": lambda a: sp.log(1 + a[0]),
    "digamma": lambda a: sp.digamma(a[0]),
    "bessel_j0": lambda a: sp.besselj(0, a[0]),
    "bessel_j1": lambda a: sp.besselj(1, a[0]),
    "zeta": lambda a: sp.zeta(a[0]),
    "beta": lambda a: sp.beta(a[0], a[1]),
}


def interval_round(x, wp: int, p: int, rm: str):
    """Round the error interval of approximation x (rel. error < 2^-(wp-16)); None if ambiguous."""
    sign, man, exp = mpf_parts(x)
    if man == 0:
        return None
    k = wp - 16
    lo = round_rat(sign, (man << k) - man, 1, exp - k, p, rm)
    hi = round_rat(sign, (man << k) + man, 1, exp - k, p, rm)
    return lo if lo == hi else None


def ziv(fn: str, args: list[Arg], p: int, rm: str, cache: dict):
    """Return (expect, extra, hard) or None if unresolved up to the top extra."""
    for i in range(len(EXTRAS) - 1):
        got = []
        for extra in (EXTRAS[i], EXTRAS[i + 1]):
            if extra not in cache:
                with mp.workprec(p + extra):
                    cache[extra] = MPF[fn]([a.mpf() for a in args])
            v = cache[extra]
            if not isinstance(v, mp.mpf):
                got.append(None)
                continue
            got.append(interval_round(v, p + extra, p, rm))
        if got[0] is not None and got[0] == got[1]:
            return got[0], EXTRAS[i], i > 0
    return None


def resid_round(fn, args, p, rm, resid: int, cache: dict):
    """f = a + resid * eps with a the p-bit point at the top precision and 0 < eps << ulp."""
    top = EXTRAS[-1]
    if top not in cache:
        with mp.workprec(p + top):
            cache[top] = MPF[fn]([x.mpf() for x in args])
    v = cache[top]
    sign, man, exp = mpf_parts(v)
    near = round_rat(sign, man, 1, exp, p, "ToEven")
    if near["k"] != "f":
        raise AssertionError("resid row off the finite range")
    mant = 0
    for i, w in enumerate(near["m"]):
        mant |= int(w, 16) << (WORD * i)
    ez = near["e"] - round_p(p)
    # check that the top approximation is within 2^-(p+top-32) of the p-bit point
    with mp.workprec(p + top + 64):
        point = mp.mpf((-mant if sign == "-" else mant, ez))
        gap = abs(v - point)
        if gap != 0 and gap > abs(point) * mp.ldexp(1, -(p + top - 32)):
            raise AssertionError(f"{fn} resid row is not near a grid point")
    k = 64
    signed = mant if sign == "+" else -mant
    total = (signed << k) + resid
    s = "-" if total < 0 else "+"
    return round_rat(s, abs(total), 1, ez - k, p, rm)


# ---------------------------------------------------------------- cases

P_LOW = [24, 53, 64, 113, 256]
P_MID = [24, 53, 113, 256, 1024]
P_ALL = [24, 53, 64, 113, 256, 1024, 4096]

CASES: list[dict] = []


def add(fn, family, specs, ps=None, exact=None, outcome=None, resid=None, note=""):
    CASES.append({
        "fn": fn, "family": family, "specs": specs, "ps": ps or P_MID,
        "exact": exact, "outcome": outcome, "resid": resid, "note": note,
    })


def frac_gamma_int(n):
    return Fraction(factorial(n - 1))


E_TOP = EMAX - 1          # 2^(EMAX-1) is the largest power of two (0.1 * 2^EMAX)
E_BOT = EMIN - 1          # 2^(EMIN-1) = 0.1 * 2^EMIN, the smallest normal


def build_cases():
    d = lambda s: ("d", s)
    # ---------------- gamma
    for s in ["0.5", "1.5", "2.5", "3.25", "7.75", "-0.5", "-2.5", "-5.75", "0.001", "33.5", "170.6", "1000.25"]:
        add("gamma", "moderate", [d(s)])
    for n in [1, 2, 3, 5, 10, 21, 30]:
        add("gamma", "exact", [d(str(n))], exact=lambda a, n=n: frac_gamma_int(n), ps=P_LOW)
    for s in ["1e-60", "1e-120", "-1e-90", "1e16"]:
        add("gamma", "cosmology", [d(s)], ps=[24, 53, 113, 256])
    add("gamma", "overflow", [d("1e60")], outcome="overflow+", ps=P_LOW)
    add("gamma", "overflow", [("pow2", 57)], outcome="overflow+", ps=P_LOW)
    add("gamma", "extreme", [("pow2", EMIN + 10)], resid=-1, ps=[24, 53, 113])
    add("gamma", "extreme", [("pow2", EMIN + 10, "-")], resid=-1, ps=[24, 53, 113])
    add("gamma", "overflow", [("pow2", E_BOT)], outcome="overflow+", ps=[24, 53])
    for n in [0, 1, 2, 10]:
        for sgn in ["+", "-"]:
            add("gamma", "near-pole", [("off", f"-{n}", None, sgn)], ps=[24, 53, 113, 256])
        add("gamma", "pole", [d(f"-{n}")], outcome="nan", ps=[53])
    add("gamma", "pole", [("z", "+")], outcome="nan", ps=[53])
    add("gamma", "pole", [("z", "-")], outcome="nan", ps=[53])
    add("gamma", "large-prec", [d("2.5")], ps=[16384, 65536])
    # ---------------- ln_gamma (positive arguments)
    for s in ["0.5", "1.5", "2.5", "3.25", "10.5", "100.125", "1e4", "0.01"]:
        add("ln_gamma", "moderate", [d(s)])
    for s in ["1", "2"]:
        add("ln_gamma", "exact", [d(s)], exact=lambda a: Fraction(0), ps=P_LOW)
    for base in ["1", "2"]:
        for sgn in ["+", "-"]:
            add("ln_gamma", "near-zero", [("off", base, None, sgn)], ps=[24, 53, 113, 256])
    for s in ["1e60", "1e120", "1e-60", "1e-120"]:
        add("ln_gamma", "cosmology", [d(s)], ps=[24, 53, 113, 256])
    add("ln_gamma", "extreme", [("pow2", (1 << 61))], ps=[24, 53])
    add("ln_gamma", "overflow", [("pow2", E_TOP)], outcome="overflow+", ps=[24, 53])
    add("ln_gamma", "extreme", [("pow2", EMIN + 10)], ps=[24, 53])
    # ruling 2: ln|Gamma| has a +Inf pole at 0 and at the negative integers
    add("ln_gamma", "pole", [d("-1")], outcome="+inf", ps=[53])
    add("ln_gamma", "pole", [("z", "+")], outcome="+inf", ps=[53])
    add("ln_gamma", "large-prec", [d("1.5")], ps=[16384])
    # ---------------- erf
    for s in ["0.5", "-0.5", "1", "2.75", "-3.5", "6", "0.0625"]:
        add("erf", "moderate", [d(s)])
    add("erf", "exact", [("z", "+")], exact=lambda a: Fraction(0), ps=P_LOW)
    add("erf", "exact", [("z", "-")], exact=lambda a: Fraction(0), ps=P_LOW, note="-0")
    for k in [-10, -40, -200]:
        add("erf", "near-zero", [("pow2", k)], ps=[24, 53, 113, 256])
    for s in ["1e-60", "-1e-120", "1e60", "-1e120"]:
        add("erf", "cosmology", [d(s)], ps=[24, 53, 113, 256],
            resid=(-1 if s == "1e60" else (1 if s == "-1e120" else None)))
    add("erf", "extreme", [("pow2", E_BOT + 4)], ps=[24, 53, 113])
    add("erf", "extreme", [("pow2", E_TOP)], resid=-1, ps=[24, 53])
    add("erf", "large-prec", [d("0.5")], ps=[16384, 65536])
    # ---------------- erfc
    for s in ["0.5", "-0.5", "1", "3", "10", "27.5", "-4"]:
        add("erfc", "moderate", [d(s)])
    add("erfc", "exact", [("z", "+")], exact=lambda a: Fraction(1), ps=P_LOW)
    for s in ["1e9", "1e-60", "-1e-120"]:
        add("erfc", "cosmology", [d(s)], ps=[24, 53, 113],
            resid=(-1 if s == "1e-60" else (1 if s == "-1e-120" else None)))
    add("erfc", "underflow", [("pow2", 32)], outcome="underflow+", ps=P_LOW)
    add("erfc", "underflow", [d("1e60")], outcome="underflow+", ps=[24, 53, 113])
    add("erfc", "extreme", [("pow2", E_BOT + 4)], resid=-1, ps=[24, 53])
    add("erfc", "extreme", [("pow2", E_TOP, "-")], resid=-1, ps=[24, 53], note="2 - tiny")
    add("erfc", "large-prec", [d("2")], ps=[16384])
    # ---------------- expm1
    for s in ["0.75", "-0.75", "1", "10.5", "-20", "700.25", "1e-5"]:
        add("expm1", "moderate", [d(s)])
    add("expm1", "exact", [("z", "+")], exact=lambda a: Fraction(0), ps=P_LOW)
    add("expm1", "exact", [("z", "-")], exact=lambda a: Fraction(0), ps=P_LOW, note="-0")
    for k in [-30, -100, -300]:
        add("expm1", "near-zero", [("pow2", k)], ps=[24, 53, 113, 256])
        add("expm1", "near-zero", [("pow2", k, "-")], ps=[24, 53, 113, 256])
    for s in ["1e-60", "-1e-120", "1e18"]:
        add("expm1", "cosmology", [d(s)], ps=[24, 53, 113])
    add("expm1", "overflow", [("pow2", 62)], outcome="overflow+", ps=P_LOW)
    add("expm1", "extreme", [d("-1e60")], resid=1, ps=[24, 53], note="-1 + tiny")
    add("expm1", "extreme", [("pow2", E_BOT + 4)], resid=1, ps=[24, 53, 113])
    add("expm1", "extreme", [("pow2", E_BOT + 4, "-")], resid=1, ps=[24, 53, 113])
    add("expm1", "large-prec", [d("0.75")], ps=[16384, 65536])
    # ---------------- log1p
    for s in ["0.75", "-0.75", "1", "10.5", "1e10", "-0.999", "1e-5"]:
        add("log1p", "moderate", [d(s)])
    add("log1p", "exact", [("z", "+")], exact=lambda a: Fraction(0), ps=P_LOW)
    add("log1p", "exact", [("z", "-")], exact=lambda a: Fraction(0), ps=P_LOW, note="-0")
    for k in [-30, -100, -300]:
        add("log1p", "near-zero", [("pow2", k)], ps=[24, 53, 113, 256])
        add("log1p", "near-zero", [("pow2", k, "-")], ps=[24, 53, 113, 256])
    for s in ["1e60", "1e120", "1e-60", "-1e-120"]:
        add("log1p", "cosmology", [d(s)], ps=[24, 53, 113])
    add("log1p", "near-pole", [("off", "-1", None, "+")], ps=[24, 53, 113])
    add("log1p", "extreme", [("pow2", E_TOP)], ps=[24, 53])
    add("log1p", "extreme", [("pow2", E_BOT + 4)], resid=-1, ps=[24, 53, 113])
    add("log1p", "pole", [d("-1")], outcome="-inf", ps=[53])
    add("log1p", "domain", [d("-2")], outcome="nan", ps=[53])
    add("log1p", "large-prec", [d("0.75")], ps=[16384, 65536])
    # ---------------- digamma
    for s in ["0.5", "1", "2.5", "10.25", "-0.5", "-2.75", "100.5", "0.01"]:
        add("digamma", "moderate", [d(s)])
    add("digamma", "near-zero", [("mp", lambda: mp.findroot(mp.digamma, 1.4616))], ps=P_ALL,
        note="positive root x0 rounded at p bits")
    for n in [0, 1, 3]:
        for sgn in ["+", "-"]:
            add("digamma", "near-pole", [("off", f"-{n}", None, sgn)], ps=[24, 53, 113, 256])
        add("digamma", "pole", [d(f"-{n}")], outcome="nan", ps=[53])
    for s in ["1e-60", "1e120", "-1e-90"]:
        add("digamma", "cosmology", [d(s)], ps=[24, 53, 113])
    add("digamma", "extreme", [("pow2", 1 << 61)], ps=[24, 53])
    add("digamma", "extreme", [("pow2", EMIN + 10)], resid=-1, ps=[24, 53], note="-1/x - gamma")
    add("digamma", "large-prec", [d("0.5")], ps=[16384, 65536])
    # ---------------- J0 / J1
    for fn in ["bessel_j0", "bessel_j1"]:
        for s in ["0.5", "1", "-2.25", "7.5", "30", "100.5", "1e-3"]:
            add(fn, "moderate", [d(s)])
        for s in ["1e-60", "-1e-120", "1e20"]:
            add(fn, "cosmology", [d(s)], ps=[24, 53, 113],
                resid=(None if s == "1e20" else (-1 if fn == "bessel_j0" else (-1 if s == "1e-60" else 1))))
        add(fn, "extreme", [("pow2", E_BOT + 8)], resid=-1, ps=[24, 53, 113])
        add(fn, "large-prec", [d("1")], ps=[16384])
    add("bessel_j0", "exact", [("z", "+")], exact=lambda a: Fraction(1), ps=P_LOW)
    add("bessel_j1", "exact", [("z", "+")], exact=lambda a: Fraction(0), ps=P_LOW)
    add("bessel_j1", "exact", [("z", "-")], exact=lambda a: Fraction(0), ps=P_LOW, note="-0")
    for k in [1, 2]:
        add("bessel_j0", "near-zero", [("mp", lambda k=k: mp.besseljzero(0, k))], ps=P_ALL,
            note=f"j0,{k} rounded at p bits")
        add("bessel_j1", "near-zero", [("mp", lambda k=k: mp.besseljzero(1, k))], ps=P_ALL,
            note=f"j1,{k} rounded at p bits")
    # ---------------- zeta
    for s in ["0.5", "1.5", "3", "4", "2.5", "-0.5", "-2.5", "50", "-30.5"]:
        add("zeta", "moderate", [d(s)])
    add("zeta", "exact-closed", [d("2")], ps=P_ALL, note="pi^2/6")
    add("zeta", "exact", [d("0")], exact=lambda a: Fraction(-1, 2), ps=P_LOW)
    for n in [1, 3, 11, 21]:
        b = Fraction(int(sp.bernoulli(n + 1).p), int(sp.bernoulli(n + 1).q))
        add("zeta", "exact", [d(f"-{n}")], exact=lambda a, b=b, n=n: -b / (n + 1), ps=P_LOW)
    add("zeta", "exact", [d("-2")], exact=lambda a: Fraction(0), ps=P_LOW, note="trivial zero, +0")
    for sgn in ["+", "-"]:
        add("zeta", "near-pole", [("off", "1", None, sgn)], ps=[24, 53, 113, 256])
    for s in ["1e-60", "-1e-120"]:
        add("zeta", "cosmology", [d(s)], ps=[24, 53, 113], resid=(-1 if s == "1e-60" else 1))
    add("zeta", "cosmology", [d("1e60")], resid=1, ps=[24, 53], note="1 + 2^-s")
    add("zeta", "extreme", [("pow2", E_TOP)], resid=1, ps=[24, 53])
    add("zeta", "pole", [d("1")], outcome="nan", ps=[53])
    add("zeta", "large-prec", [d("2")], ps=[16384, 65536])
    # ---------------- beta
    for a, b in [("0.5", "0.5"), ("1.5", "2.25"), ("3.5", "0.25"), ("-0.5", "0.25"), ("10.5", "20.25")]:
        add("beta", "moderate", [d(a), d(b)])
    for a, b in [(1, 1), (2, 3), (5, 7), (10, 12)]:
        add("beta", "exact", [d(str(a)), d(str(b))],
            exact=lambda x, a=a, b=b: Fraction(factorial(a - 1) * factorial(b - 1), factorial(a + b - 1)),
            ps=P_LOW)
    add("beta", "exact", [("pow2", -100), d("1")], exact=lambda x: Fraction(2) ** 100, ps=P_LOW,
        note="B(a,1) = 1/a")
    add("beta", "exact", [("pow2", EMIN + 10), d("1")], exact=None, outcome="pow2:%d" % (-(EMIN + 10)),
        ps=[24, 53], note="B(a,1) = 1/a = 2^-(EMIN+10)")
    add("beta", "cosmology", [d("1e10"), d("1e10")], ps=[24, 53])
    add("beta", "underflow", [d("1e60"), d("1e60")], outcome="underflow+", ps=[24, 53])
    add("beta", "cosmology", [d("1e-60"), d("1e-60")], ps=[24, 53, 113])
    add("beta", "pole", [("z", "+"), d("1")], outcome="nan", ps=[53])
    add("beta", "large-prec", [d("0.5"), d("1.5")], ps=[16384])
    # ---------------- ln_gamma, negative arguments (ruling 2; appended so earlier ids stay put)
    for s in ["-0.5", "-2.5", "-5.75", "-100.25"]:
        add("ln_gamma", "negative", [d(s)])
    for n in [1, 3]:
        for sgn in ["+", "-"]:
            add("ln_gamma", "near-pole", [("off", f"-{n}", None, sgn)], ps=[24, 53, 113, 256])
    for s in ["-2", "-10"]:
        add("ln_gamma", "pole", [d(s)], outcome="+inf", ps=[24, 53])
    add("ln_gamma", "pole", [("z", "-")], outcome="+inf", ps=[24, 53])


def resolve_specs(specs, p):
    out = []
    for s in specs:
        if s[0] == "off":
            k = -(p // 2)
            # base +/- 2^k: exact in p bits for |base| <= 2^(p/2 - 4)
            sgn = s[3]
            b = Fraction(int(s[1]))
            v = b + (Fraction(2) ** k if sgn == "+" else -(Fraction(2) ** k))
            sign = "-" if v < 0 else "+"
            av = abs(v)
            n, den = av.numerator, av.denominator
            ez = -(den.bit_length() - 1)
            out.append(Arg(sign, n, ez))
        else:
            out.append(arg_at(s, p))
    return out


def outcome_expect(outcome: str, p: int, rm: str) -> dict:
    if outcome == "nan":
        return nan()
    if outcome == "-inf":
        return inf("-")
    if outcome == "+inf":
        return inf("+")
    if outcome == "overflow+":
        return overflow("+", p, rm)
    if outcome == "underflow+":
        # 0 < f < 2^(EMIN - p - 8): rounds to +0 or to the smallest subnormal
        return round_rat("+", 1, 1, EMIN - 4 * p - 1024, p, rm)
    if outcome.startswith("pow2:"):
        k = int(outcome[5:])
        return round_rat("+", 1, 1, k, p, rm)
    raise ValueError(outcome)


def main() -> None:
    build_cases()
    rows = []
    hard = []
    unresolved = []
    sym_checked = 0
    sym_failed = []
    t0 = time.time()
    for ci, c in enumerate(CASES):
        for p in c["ps"]:
            args = resolve_specs(c["specs"], p)
            cache: dict = {}
            sym_done = False
            for rm in MODES:
                extra, is_hard, how = None, False, "ziv"
                if c["outcome"] is not None:
                    exp_ = outcome_expect(c["outcome"], p, rm)
                    how = "outcome"
                elif c["exact"] is not None:
                    zs = args[0].sign if c["fn"] in ("erf", "expm1", "log1p", "bessel_j1") else "+"
                    exp_ = round_fraction(c["exact"](args), p, rm, zero_sign=zs)
                    how = "exact"
                else:
                    r = ziv(c["fn"], args, p, rm, cache)
                    if r is None:
                        if c["resid"] is None:
                            unresolved.append((c["fn"], c["specs"], p, rm))
                            continue
                        exp_ = resid_round(c["fn"], args, p, rm, c["resid"], cache)
                        extra, is_hard, how = EXTRAS[-1], True, "resid"
                    else:
                        exp_, extra, is_hard = r
                    if not sym_done and p in (53, 113) and c["family"] in ("moderate", "near-zero", "exact-closed") and ci % 3 == 0:
                        sym_done = True
                        sym_checked += 1
                        dps = int(p * 0.302) + 200
                        try:
                            sv = sp.N(SYM[c["fn"]]([a.sym() for a in args]), dps)
                            with mp.workprec(p + 512):
                                mv = cache[min(cache)]
                                sv_m = mp.mpf(str(sp.re(sv)))
                                if abs(sv_m - mv) > abs(mv) * mp.ldexp(1, -(p + 64)):
                                    sym_failed.append((c["fn"], c["specs"], p))
                        except Exception as ex:  # sympy cannot evaluate: record, do not hide
                            sym_failed.append((c["fn"], c["specs"], p, repr(ex)))
                row = {
                    "id": f"s9-{c['fn']}-{p}-{rm}-{len(rows)}",
                    "fn": c["fn"], "family": c["family"], "p": p, "rm": rm,
                    "args": [a.enc(p) for a in args],
                    "expect": exp_, "how": how,
                    "extra": extra, "hard": is_hard,
                }
                if c["note"]:
                    row["note"] = c["note"]
                rows.append(row)
                if is_hard:
                    hard.append(f"{row['id']} fn={c['fn']} family={c['family']} p={p} rm={rm} "
                                f"how={how} extra={extra} specs={c['specs'][0] if c['specs'][0][0] != 'mp' else c['note']}")
        print(f"[{time.time() - t0:7.1f}s] case {ci + 1}/{len(CASES)} {c['fn']} {c['family']} rows={len(rows)}",
              file=sys.stderr, flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w") as f:
        for r in rows:
            f.write(json.dumps(r, separators=(",", ":")) + "\n")
    with HARD_OUT.open("w") as f:
        f.write("# S9 hard cases: the p+256/p+512 Ziv pair disagreed or straddled a grid point\n")
        for h in hard:
            f.write(h + "\n")
    print(f"rows={len(rows)} hard={len(hard)} unresolved={len(unresolved)} "
          f"sympy_checked={sym_checked} sympy_failed={len(sym_failed)}")
    for u in unresolved:
        print("UNRESOLVED", u)
    for s in sym_failed:
        print("SYMPY", s)
    if unresolved or sym_failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
