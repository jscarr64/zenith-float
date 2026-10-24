#!/usr/bin/env python3
"""R13 goldens: div, reciprocal, hypot (directed, all modes) and div_exact.

Writes zenith-float-num/tests/oracle/goldens/r13_div_hypot.jsonl in the canon of
gen_s9_specials_goldens.py (little-endian 64-bit limbs, value = 0.1m * 2^e,
i64 exponent range +-(2^62 - 1), subnormal grid 2^(EMIN - p)).

Two oracles per row:
  * mpmath at p+256 and at p+512 bits (mp.fdiv / mp.hypot). The error interval of
    each must round to one p-bit result and the two must agree ("how": "ziv").
    When the p+256 and p+512 values are the same number the value is exact
    (checked with exact=True mpmath arithmetic: v*v == x*x + y*y, v*y == x) (an
    exact quotient, a Pythagorean triple, a tie, a zero); rounding it at p bits
    must give the model's result ("how": "exact").
  * an exact integer model (rational rounding for div, isqrt with a sticky bit for
    hypot). Every row's expect comes from it, and it must equal the mpmath result.
Rows where neither mpmath check settles the result (a value within 2^-(p+256) of
a p-bit number or midpoint without being one, such as x >> y in hypot) are "how": "model"; the
script reports their count and fails if any row's two oracles disagree.

Arguments may carry more bits than p ("pa"), which is how exact ties are built.

Run with ~/.local/oracle-venv/bin/python (mpmath 1.3.0).
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from math import isqrt
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("s9", HERE / "gen_s9_specials_goldens.py")
s9 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s9)

EMAX, EMIN, MODES = s9.EMAX, s9.EMIN, s9.MODES
round_rat, round_p, zero, Arg = s9.round_rat, s9.round_p, s9.zero, s9.Arg
OUT = s9.ROOT / "zenith-float-num" / "tests" / "oracle" / "goldens" / "r13_div_hypot.jsonl"


# ------------------------------------------------------------- arguments

def arg(spec, pa: int) -> Arg:
    """("i", n, ez[, sign]) exact n*2^ez; ("max",) largest finite at pa; ("sub", k) k*2^(EMIN-pa);
    otherwise a gen_s9 spec ("d", "1.5"), ("pow2", k[, sign]), ("z", sign)."""
    kind = spec[0]
    if kind == "i":
        a = Arg(spec[3] if len(spec) > 3 else "+", spec[1], spec[2])
        if a.n.bit_length() > pa:
            raise ValueError(f"{spec} needs more than {pa} bits")
        return a
    if kind == "max":
        return Arg(spec[1] if len(spec) > 1 else "+", (1 << pa) - 1, EMAX - pa)
    if kind == "sub":
        return Arg(spec[2] if len(spec) > 2 else "+", spec[1], EMIN - pa)
    return s9.arg_at(spec, pa)


def mpf(a: Arg):
    return a.mpf()


def num(d: dict) -> int:
    m = 0
    for i, w in enumerate(d["m"]):
        m |= int(w, 16) << (64 * i)
    return m


# ------------------------------------------------------------- exact models

def div_model(x: Arg, y: Arg, p: int, rm: str) -> dict:
    s = "+" if x.sign == y.sign else "-"
    if x.n == 0:
        return zero(s, p)
    return round_rat(s, x.n, y.n, x.ez - y.ez, p, rm)


def recip_model(x: Arg, p: int, rm: str) -> dict:
    return round_rat(x.sign, 1, x.n, -x.ez, p, rm)


def hypot_model(x: Arg, y: Arg, p: int, rm: str) -> dict:
    if x.n == 0 and y.n == 0:
        return zero("+", p)
    if x.n == 0 or y.n == 0:
        a = y if x.n == 0 else x
        return round_rat("+", a.n, 1, a.ez, p, rm)
    if x.ez + x.n.bit_length() < y.ez + y.n.bit_length():
        x, y = y, x
    g_bits = p + 128
    if (x.ez + x.n.bit_length()) - (y.ez + y.n.bit_length()) > g_bits:
        # y^2 lies below 2^-(2 g_bits) of x^2: a stand-in 1 below x^2 * 4^g_bits rounds the same,
        # since |x| (at most p + 64 bits) and the midpoints next to it are its only nearby breaks.
        big = (x.n * x.n) << (2 * g_bits)
        s_val, base = big + 1, x.ez - g_bits
    else:
        k = min(x.ez, y.ez)
        s_val = ((x.n * x.n) << (2 * (x.ez - k))) + ((y.n * y.n) << (2 * (y.ez - k)))
        base = k
    g = max(0, p + 70 - s_val.bit_length() // 2)
    t = s_val << (2 * g)
    f = isqrt(t)
    sticky = 1 if f * f != t else 0
    return round_rat("+", 2 * f + sticky, 1, base - g - 1, p, rm)


def exact_model(x: Arg, y: Arg, p: int):
    """div_exact: the quotient if every mode gives one finite value, else None."""
    if x.n == 0:
        return zero("+" if x.sign == y.sign else "-", p)
    lo = div_model(x, y, p, "ToZero")
    hi = div_model(x, y, p, "FromZero")
    if lo != hi or lo["k"] != "f" or num(lo) == 0:
        return None
    return lo


# ------------------------------------------------------------- mpmath oracle

MPF = {
    "div": lambda a: mp.fdiv(a[0], a[1]),
    "reciprocal": lambda a: mp.fdiv(1, a[0]),
    "hypot": lambda a: mp.hypot(a[0], a[1]),
}


def is_exact(fn, args, v) -> bool:
    """mpmath's value v satisfies the defining identity exactly (exact=True arithmetic)."""
    with mp.workprec(max(x.n.bit_length() for x in args) + 64):
        a = [mpf(x) for x in args]  # exact: the context precision covers every argument bit
    if fn == "hypot":
        if v == 0:
            return a[0] == 0 and a[1] == 0
        nz = [x for x in args if x.n != 0]
        tops = [x.ez + x.n.bit_length() for x in nz]
        if len(tops) == 2 and abs(tops[0] - tops[1]) > 8 * 65536:
            return False  # x^2 + y^2 would need ~2^63 bits: not exact anyway at these p
        lhs = mp.fmul(v, v, exact=True)
        rhs = mp.fadd(mp.fmul(a[0], a[0], exact=True), mp.fmul(a[1], a[1], exact=True), exact=True)
        return lhs == rhs
    if fn == "div":
        return mp.fmul(v, a[1], exact=True) == a[0]
    return mp.fmul(v, a[0], exact=True) == 1


def mp_check(fn, args, p, rm, expect):
    """'ziv' | 'exact' | 'model'; raises on a disagreement with the model."""
    vals = []
    for extra in (256, 512):
        with mp.workprec(p + extra):
            vals.append(MPF[fn]([mpf(a) for a in args]))
    got = [s9.interval_round(v, p + e, p, rm) if v != 0 else None for v, e in zip(vals, (256, 512))]
    if got[0] is not None and got[0] == got[1]:
        if got[0] != expect:
            raise AssertionError(f"{fn} {p} {rm}: mpmath {got[0]} model {expect}")
        return "ziv"
    # exact: the p+256 and p+512 values are the same number (mpmath rounds each operation
    # correctly, so the value is exact); rounding it directly must give the model's result.
    if vals[0] == vals[1] and is_exact(fn, args, vals[0]):
        if vals[0] == 0:
            ok = expect["k"] == "f" and num(expect) == 0
        else:
            sign, man, exp = s9.mpf_parts(vals[0])
            ok = round_rat(sign, man, 1, exp, p, rm) == expect
        if not ok:
            raise AssertionError(f"{fn} {p} {rm}: exact mpmath value rounds unlike model {expect}")
        return "exact"
    if got[0] is not None and got[1] is not None and got[0] != got[1]:
        raise AssertionError(f"{fn} {p} {rm}: mpmath p+256 {got[0]} vs p+512 {got[1]}")
    return "model"


# ------------------------------------------------------------- cases

P_STD = [24, 53, 64, 113, 256, 1024, 4096]
P_BIG = [16384, 65536]
CASES: list[dict] = []


def add(fn, family, specs, ps=None, pa_extra=0):
    CASES.append({"fn": fn, "family": family, "specs": specs, "ps": ps or P_STD, "pa_extra": pa_extra})


def build_cases():
    d = lambda s: ("d", s)
    i = lambda n, ez=0, s="+": ("i", n, ez, s)
    E_TOP = EMAX - 1
    cosmo = [("6.02214076e23", "1.380649e-23"), ("1e60", "1e-60"), ("1.989e30", "5.972e24"),
             ("4.4e26", "8.8e26"), ("1.616255e-35", "5.391247e-44"), ("1e120", "1e120"),
             ("1e-120", "3e-120"), ("-2.99792458e8", "6.6743e-11")]
    # ---------------- hypot
    for a, b in [(3, 4), (5, 12), (8, 15), (20, 21), (119, 120), (4961, 6480)]:
        add("hypot", "triple", [i(a), i(b)])
        add("hypot", "triple", [i(a, 0, "-"), i(b, 0, "-")], ps=[24, 113])
    for k in [EMAX - 3, EMAX - 4, EMIN + 1, EMIN - 40]:
        add("hypot", "triple-extreme", [i(3, k), i(4, k)])
    add("hypot", "triple-subnormal", [("sub", 3), ("sub", 4)])
    for a, b in [("1", "1"), ("1", "2"), ("0.1", "0.2"), ("-1.5", "2.25"), ("7", "1e-3")]:
        add("hypot", "moderate", [d(a), d(b)])
    for a, b in cosmo:
        add("hypot", "cosmology", [d(a), d(b)])
    for k in [-40, -200, -5000]:
        add("hypot", "x>>y", [d("1"), ("pow2", k)])
    add("hypot", "x>>y", [("pow2", E_TOP), ("pow2", EMIN + 2)])
    add("hypot", "x>>y", [("pow2", EMIN + 2, "-"), ("pow2", E_TOP)])
    add("hypot", "near-limit", [("pow2", E_TOP), ("pow2", E_TOP)])
    add("hypot", "near-limit", [("pow2", EMIN - 1), ("pow2", EMIN - 1)])
    add("hypot", "near-limit", [("pow2", E_TOP - 1), ("pow2", E_TOP - 3)])
    add("hypot", "overflow", [("max",), ("max",)])
    add("hypot", "overflow", [("max",), ("pow2", EMIN + 5)])
    add("hypot", "overflow", [("max", "-"), ("pow2", E_TOP - 70)])
    add("hypot", "subnormal", [("sub", 1), ("sub", 1)])
    add("hypot", "subnormal", [("sub", 7), ("sub", 2)])
    add("hypot", "subnormal", [("pow2", EMIN - 1), ("sub", 1)])
    add("hypot", "zero", [("z", "+"), ("z", "-")])
    add("hypot", "zero", [("z", "-"), d("-2.5")])
    add("hypot", "zero", [("sub", 1, "-"), ("z", "+")])
    # ties: |x| has p+1 bits and is a midpoint at p
    add("hypot", "tie", [("tie",), ("z", "+")], pa_extra=1)
    add("hypot", "tie", [("tie", "-"), ("pow2", -4000)], pa_extra=1)
    add("hypot", "tie", [("tie3",), ("z", "+")], pa_extra=1)
    for a, b in [("1", "1"), ("3", "4"), ("1e60", "7")]:
        add("hypot", "large-prec", [d(a), d(b)] if a != "3" else [i(3), i(4)], ps=P_BIG)
    # ---------------- div
    for a, b in [(1, 3), (22, 7), (2, -3), (6, 3), (1, 1024), ((1 << 64) - 1, 3), (10, 4), (7, 7)]:
        sa, sb = ("-" if a < 0 else "+"), ("-" if b < 0 else "+")
        add("div", "moderate", [i(abs(a), 0, sa), i(abs(b), 0, sb)],
            ps=P_STD if abs(a) < (1 << 24) else [64, 113, 256, 1024, 4096])
    for a, b in cosmo:
        add("div", "cosmology", [d(a), d(b)])
    add("div", "overflow", [("max",), ("pow2", -1)])
    add("div", "overflow", [("pow2", E_TOP, "-"), ("pow2", EMIN + 2)])
    add("div", "underflow", [("pow2", EMIN + 2), ("pow2", E_TOP)])
    add("div", "subnormal", [("sub", 1), d("3")])
    add("div", "tie", [("sub", 1), d("2")])
    add("div", "tie", [("sub", 3), d("-2")])
    add("div", "tie", [("tie",), d("1")], pa_extra=1)
    add("div", "tie", [("tie3",), d("-4")], pa_extra=1)
    add("div", "near-limit", [("max",), d("1.5")])
    add("div", "near-limit", [("pow2", EMIN - 1), d("0.75")])
    add("div", "zero", [("z", "-"), d("3")])
    for a, b in [("1", "3"), ("1e60", "7")]:
        add("div", "large-prec", [d(a), d(b)], ps=P_BIG)
    # ---------------- reciprocal
    for a in ["3", "-7", "0.1", "1024", "1e60", "-1e-120", "6.6743e-11"]:
        add("reciprocal", "moderate", [d(a)])
    add("reciprocal", "overflow", [("sub", 1)])
    add("reciprocal", "overflow", [("pow2", EMIN - 1, "-")])
    add("reciprocal", "underflow", [("pow2", E_TOP)])
    add("reciprocal", "underflow", [("max",)])
    add("reciprocal", "large-prec", [d("3")], ps=P_BIG)


def resolve(spec, p, pa):
    kind = spec[0]
    if kind in ("tie", "tie3"):
        # (2^p + 1) and 3*(2^(p-1)) + 1... both p+1 bits, a midpoint at p bits
        n = (1 << p) + 1 if kind == "tie" else (3 << (p - 1)) + 1
        return Arg(spec[1] if len(spec) > 1 else "+", n, -p)
    return arg(spec, pa)


def main() -> None:
    build_cases()
    rows, how_n = [], {"ziv": 0, "exact": 0, "model": 0}
    t0 = time.time()
    for ci, c in enumerate(CASES):
        for p in c["ps"]:
            pa = p + c["pa_extra"]
            args = [resolve(s, p, pa) for s in c["specs"]]
            for rm in MODES:
                if c["fn"] == "div":
                    exp_ = div_model(args[0], args[1], p, rm)
                elif c["fn"] == "reciprocal":
                    exp_ = recip_model(args[0], p, rm)
                else:
                    exp_ = hypot_model(args[0], args[1], p, rm)
                how = mp_check(c["fn"], args, p, rm, exp_)
                how_n[how] += 1
                rows.append({
                    "id": f"r13-{c['fn']}-{p}-{rm}-{len(rows)}", "fn": c["fn"], "family": c["family"],
                    "p": p, "pa": pa, "rm": rm, "args": [a.enc(pa) for a in args], "expect": exp_, "how": how,
                })
            if c["fn"] == "div":
                ex = exact_model(args[0], args[1], p)
                rows.append({
                    "id": f"r13-div_exact-{p}-{len(rows)}", "fn": "div_exact", "family": c["family"],
                    "p": p, "pa": pa, "rm": "-", "args": [a.enc(pa) for a in args],
                    "expect": ex if ex is not None else {"k": "none"}, "how": "model",
                })
        print(f"[{time.time() - t0:7.1f}s] case {ci + 1}/{len(CASES)} {c['fn']} {c['family']} rows={len(rows)}",
              file=sys.stderr, flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w") as f:
        for r in rows:
            f.write(json.dumps(r, separators=(",", ":")) + "\n")
    per_fn = {}
    for r in rows:
        per_fn[r["fn"]] = per_fn.get(r["fn"], 0) + 1
    print(f"rows={len(rows)} per_fn={per_fn} how={how_n}")


if __name__ == "__main__":
    main()
