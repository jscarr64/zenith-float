#!/usr/bin/env python3
"""Borrowed unary goldens: ExactNum::abs_ref / neg_ref / floor_ref and ldexp_ref in all 7 rounding modes.

Writes zenith-float-num/tests/oracle/goldens/unary_ref.jsonl in the canon of gen_s9_specials_goldens.py
(little-endian 64-bit limbs, value = 0.1m * 2^e, exponent range +-(2^62 - 1), subnormal grid 2^(EMIN - p)).

Oracles, in this order:
  * FLINT (python-flint arf/arb): the exact result (|x|, -x, floor(x), x * 2^n) as (man, exp).
  * mpmath: the same exact result independently (fabs, fneg, floor, ldexp at enough bits); it must
    equal FLINT's.
  * an exact integer model rounds that value at p bits (s9.round_rat, with ties-to-odd added and
    RoundingMode::None as truncation at p rounded up to whole words, with an overflow to Inf). For in-range normal results in the
    five IEEE modes, mpmath's own rounding (libmp.normalize) must give the same bits.
Rows: every mode, both signs, signed zeros, Inf, NaN, subnormal operands and results, exact and inexact
ldexp, ties at p bits, overflow at EXPONENT_MAX and underflow below EXPONENT_MIN (Exponent extremes too).

Run with the pinned team oracle env (python-flint 0.9.0, mpmath 1.3.0).
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import flint
import mpmath as mp
from mpmath import libmp

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("r13", HERE / "gen_r13_goldens.py")
r13 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(r13)
s9 = r13.s9

EMAX, EMIN, WORD = s9.EMAX, s9.EMIN, s9.WORD
MODES7 = ["ToEven", "ToOdd", "Up", "Down", "ToZero", "FromZero", "None"]
OUT = s9.ROOT / "zenith-float-num" / "tests" / "oracle" / "goldens" / "unary_ref.jsonl"
I64_MAX, I64_MIN = (1 << 63) - 1, -(1 << 63)

_away5 = s9.away


def away7(rm: str, positive: bool, rem_cmp: int, q_odd: bool) -> bool:
    if rm == "ToOdd":
        return rem_cmp > 0 or (rem_cmp == 0 and not q_odd)
    return _away5(rm, positive, rem_cmp, q_odd)


s9.away = away7  # round_rat looks `away` up in its module


def model_round(sign: str, n: int, ez: int, p: int, rm: str) -> dict:
    if rm == "None":
        # No rounding: the stored words are truncated, and an overflow is Inf (RoundingMode::None does
        # not round toward zero; see ExactNum::result_to_ext_rm).
        if n.bit_length() + ez > EMAX:
            return s9.inf(sign)
        return s9.round_rat(sign, n, 1, ez, s9.round_p(p), "ToZero")
    return s9.round_rat(sign, n, 1, ez, p, rm)


def flint_exact(fn: str, x: r13.Arg, k: int):
    """(sign, n, ez) of the exact result from FLINT, n == 0 for a zero (sign then from the model)."""
    if x.n == 0:
        return None
    a = flint.arf((-x.n if x.sign == "-" else x.n, x.ez))
    old = flint.ctx.prec
    flint.ctx.prec = x.n.bit_length() + 64
    try:
        if fn == "abs":
            r = abs(a)
        elif fn == "neg":
            r = -a
        elif fn == "ldexp":
            r = a * flint.arf((1, k))
        else:
            b = flint.arb(a).floor()
            if b.rad() != 0:
                raise AssertionError("FLINT floor is not exact")
            r = b.mid()
    finally:
        flint.ctx.prec = old
    man, exp = (int(v) for v in r.man_exp())
    return ("-" if man < 0 else "+", abs(man), exp)


def mp_exact(fn: str, x: r13.Arg, k: int):
    if x.n == 0:
        return None
    with mp.workprec(x.n.bit_length() + 64):
        v = x.mpf()
        r = {"abs": lambda: mp.fabs(v), "neg": lambda: -v, "ldexp": lambda: mp.ldexp(v, k),
             "floor": lambda: mp.floor(v)}[fn]()
    if r == 0:
        return ("+", 0, 0)
    return s9.mpf_parts(r)


def canon_odd(t):
    sign, n, ez = t
    if n == 0:
        return (0,)
    tz = (n & -n).bit_length() - 1
    return (sign, n >> tz, ez + tz)


def mp_round_check(sign: str, n: int, ez: int, p: int, rm: str, expect: dict) -> str:
    rnd = {"ToEven": "n", "Up": "c", "Down": "f", "ToZero": "d", "FromZero": "u"}.get(rm)
    if rnd is None or expect["k"] != "f" or r13.num(expect) == 0 or expect["e"] in (EMIN, EMAX):
        return "model"
    s, man, exp, bc = libmp.normalize(1 if sign == "-" else 0, n, ez, n.bit_length(), p, rnd)
    got = s9.round_rat("-" if s else "+", int(man), 1, int(exp), p, "ToZero")
    if got != expect:
        raise AssertionError(f"mpmath rounding {got} vs model {expect}")
    return "mpmath"


def expect_of(fn: str, x: r13.Arg, k: int, p: int, pa: int, rm: str):
    fl, mpv = flint_exact(fn, x, k), mp_exact(fn, x, k)
    if fl is not None and canon_odd(fl) != canon_odd(mpv):
        raise AssertionError(f"{fn} FLINT {fl} vs mpmath {mpv}")
    if fn == "ldexp":
        if x.n == 0:
            return s9.zero(x.sign, p), "exact"
        sign, n, ez = fl
        e = model_round(sign, n, ez, p, rm)
        return e, mp_round_check(sign, n, ez, p, rm, e)
    # abs / neg / floor are exact and keep the operand's storage
    if fl is None or fl[1] == 0:
        zs = {"abs": "+", "neg": "-" if x.sign == "+" else "+", "floor": x.sign}[fn]
        if fn == "floor" and fl is not None:
            zs = x.sign  # floor of a positive fraction is +0
        return s9.zero(zs, pa), "exact"
    sign, n, ez = fl
    return s9.round_rat(sign, n, 1, ez, pa, "ToZero"), "exact"


def top(a: r13.Arg) -> int:
    return a.ez + a.n.bit_length()


D = lambda s: ("d", s)
VALUES = [D("1.5"), D("-2.75"), D("0.1"), D("-0.1"), D("1"), D("-1"), D("-0.5"), D("3.25"), D("-3.75"),
          D("1e60"), D("-1e-30"), D("6.02214076e23"), D("-123456.789"), ("z", "+"), ("z", "-"),
          ("sub", 1), ("sub", 3, "-"), ("max",), ("max", "-"), ("pow2", EMIN + 1), ("pow2", EMAX - 1, "-")]
TIES = [("tie",), ("tie", "-"), ("tie3",), ("tie3", "-")]
P_LD = [24, 64, 256]
P_UN = [24, 53, 64, 113, 256]


def ldexp_ns(x: r13.Arg, p: int) -> list[int]:
    ns = [0, 1, -1, 37, -1000, 1 << 40, I64_MAX, I64_MIN]
    if x.n:
        t = top(x)
        ns += [EMAX - t, EMAX - t + 1, EMIN - t, EMIN - t - 5, EMIN - t - p - 1, EMIN - t - 4 * p - 100]
    return [n for n in dict.fromkeys(ns) if I64_MIN <= n <= I64_MAX]


def specials(fn: str):
    """Inf and NaN operands: the result kind follows the by-value op (no rounding)."""
    out = []
    for s in "+-":
        r = {"abs": s9.inf("+"), "neg": s9.inf("-" if s == "+" else "+"), "floor": s9.inf(s),
             "ldexp": s9.inf(s)}[fn]
        out.append(({"k": "inf", "s": s}, r))
    out.append(({"k": "nan", "err": "InvalidArgument"}, s9.nan("InvalidArgument")))
    return out


def main() -> None:
    n = 0
    how_n: dict[str, int] = {}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w") as f:
        def emit(row):
            nonlocal n
            row["id"] = f"unary-{row['fn']}-{n}"
            f.write(json.dumps(row, separators=(",", ":")) + "\n")
            how_n[row["how"]] = how_n.get(row["how"], 0) + 1
            n += 1

        for fn in ("abs", "neg", "floor", "ldexp"):
            for arg, want in specials(fn):
                for rm in (MODES7 if fn == "ldexp" else ["None"]):
                    emit({"fn": fn, "family": "special", "p": 64, "pa": 64, "rm": rm, "k": 3,
                          "args": [arg], "expect": want, "how": "exact"})
        for fn in ("abs", "neg", "floor"):
            for p in P_UN:
                for extra in (0, 70):
                    pa = p + extra
                    for spec in VALUES:
                        x = r13.resolve(spec, p, pa)
                        e, how = expect_of(fn, x, 0, p, pa, "None")
                        emit({"fn": fn, "family": "value", "p": p, "pa": pa, "rm": "None", "k": 0,
                              "args": [x.enc(pa)], "expect": e, "how": how})
        for p in P_LD:
            for extra, specs in ((0, VALUES), (70, VALUES), (1, TIES)):
                pa = p + extra
                for spec in specs:
                    x = r13.resolve(spec, p, pa)
                    ns = ldexp_ns(x, p) if extra != 1 else [0, -3, 5]
                    for k in ns:
                        for rm in MODES7:
                            e, how = expect_of("ldexp", x, k, p, pa, rm)
                            emit({"fn": "ldexp", "family": "tie" if extra == 1 else "value", "p": p,
                                  "pa": pa, "rm": rm, "k": k, "args": [x.enc(pa)], "expect": e, "how": how})
    print(f"rows={n} how={how_n}")


if __name__ == "__main__":
    main()
