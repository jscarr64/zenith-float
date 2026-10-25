#!/usr/bin/env python3
"""Borrowing-ops goldens: ExactNum::add_ref / sub_ref / mul_ref / div_ref in every rounding mode.

Writes zenith-float-num/tests/oracle/goldens/borrow_ops.jsonl (canon of gen_s9_specials_goldens.py),
one row per line, streamed to disk as it is built.

Two oracles per row, as in gen_r13_goldens.py:
  * mpmath at p+256 and p+512 bits (mp.fadd / fsub / fmul / fdiv). Both error intervals must round
    to the same p-bit result ("ziv"), or the two values are the same exact number ("exact",
    checked with exact=True arithmetic) and rounding it gives the model's result.
  * an exact integer model (rational rounding). It must equal every mpmath result; rows mpmath cannot
    settle (x + y with |y| far below 2^-(p+512)|x|) are "model" and are counted.
Rows with "alias": true pass the same argument as both operands (the test calls op(&x, &x)).

Run with ~/.local/oracle-venv/bin/python (mpmath 1.3.0).
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("r13", HERE / "gen_r13_goldens.py")
r13 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(r13)
s9 = r13.s9

EMAX, EMIN, MODES = s9.EMAX, s9.EMIN, s9.MODES
round_rat, zero, Arg, num = s9.round_rat, s9.zero, s9.Arg, r13.num
OUT = s9.ROOT / "zenith-float-num" / "tests" / "oracle" / "goldens" / "borrow_ops.jsonl"


def neg(a: Arg) -> Arg:
    return Arg("-" if a.sign == "+" else "+", a.n, a.ez)


def top(a: Arg) -> int:
    return a.ez + a.n.bit_length()


def sum_model(x: Arg, y: Arg, p: int, rm: str) -> dict:
    if x.n == 0 and y.n == 0:
        s = x.sign if x.sign == y.sign else ("-" if rm == "Down" else "+")
        return zero(s, p)
    if x.n == 0 or y.n == 0:
        a = y if x.n == 0 else x
        return round_rat(a.sign, a.n, 1, a.ez, p, rm)
    if top(x) < top(y):
        x, y = y, x
    g = 4 * p + 512
    if top(y) < x.ez - g:
        # |y| < 2^-g ulp(x): a stand-in of one unit at 2^(x.ez - g) on y's side rounds the same.
        n = (x.n << g) + (1 if x.sign == y.sign else -1)
        return round_rat(x.sign, n, 1, x.ez - g, p, rm)
    k = min(x.ez, y.ez)
    v = (x.n if x.sign == "+" else -x.n) << (x.ez - k)
    v += (y.n if y.sign == "+" else -y.n) << (y.ez - k)
    if v == 0:
        return zero("-" if rm == "Down" else "+", p)
    return round_rat("+" if v > 0 else "-", abs(v), 1, k, p, rm)


def model(fn, x: Arg, y: Arg, p: int, rm: str) -> dict:
    if fn == "add":
        return sum_model(x, y, p, rm)
    if fn == "sub":
        return sum_model(x, neg(y), p, rm)
    s = "+" if x.sign == y.sign else "-"
    if fn == "mul":
        if x.n == 0 or y.n == 0:
            return zero(s, p)
        return round_rat(s, x.n * y.n, 1, x.ez + y.ez, p, rm)
    if y.n == 0:
        # IEEE 754: x / ±0 is a signed Inf (the division-by-zero error is in div_exact), 0 / 0 is NaN.
        return s9.nan("InvalidArgument") if x.n == 0 else s9.inf(s)
    return r13.div_model(x, y, p, rm)


MPF = {"add": mp.fadd, "sub": mp.fsub, "mul": mp.fmul, "div": mp.fdiv}


def is_exact(fn, x: Arg, y: Arg, v) -> bool:
    a, b = x.mpf(), y.mpf()
    if fn in ("add", "sub"):
        if x.n and y.n and abs(top(x) - top(y)) > 1 << 20:
            return False  # an exact sum would need ~2^62 bits; not exact at these p anyway
        return mp.fsub(v, MPF[fn](a, b, exact=True), exact=True) == 0
    if fn == "mul":
        return v == mp.fmul(a, b, exact=True)
    return mp.fmul(v, b, exact=True) == a


def mp_check(fn, x: Arg, y: Arg, p: int, rm: str, expect: dict) -> str:
    if expect["k"] in ("nan", "inf") and y.n == 0:
        return "model"
    vals = []
    for extra in (256, 512):
        with mp.workprec(p + extra):
            vals.append(MPF[fn](x.mpf(), y.mpf()))
    got = [s9.interval_round(v, p + e, p, rm) if v != 0 else None for v, e in zip(vals, (256, 512))]
    if got[0] is not None and got[0] == got[1]:
        if got[0] != expect:
            raise AssertionError(f"{fn} {p} {rm}: mpmath {got[0]} model {expect}")
        return "ziv"
    if vals[0] == vals[1] and is_exact(fn, x, y, vals[0]):
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


P_STD = [24, 53, 64, 113, 256, 1024]
CASES: list[dict] = []


def add(fns, family, specs, ps=None, alias=False, pa_extra=0):
    for fn in fns:
        CASES.append({"fn": fn, "family": family, "specs": specs, "ps": ps or P_STD,
                      "alias": alias, "pa_extra": pa_extra})


ALL = ("add", "sub", "mul", "div")


def build_cases():
    d = lambda s: ("d", s)
    E_TOP = EMAX - 1
    for a, b in [("1.5", "0.1"), ("-2.75", "1e-3"), ("0.1", "-0.3"), ("-7", "-11"), ("1e60", "3e-40"),
                 ("6.02214076e23", "-1.380649e-23")]:
        add(ALL, "moderate", [d(a), d(b)])
    add(ALL, "moderate-big-p", [d("0.1"), d("-3")], ps=[4096])
    # cancellation, exact zero sums (sign by mode), signed zeros
    add(("add",), "cancel", [d("1.25"), d("-1.25")])
    add(("sub",), "cancel", [d("-0.1"), d("-0.1")])
    add(("add", "sub"), "near-cancel", [d("1"), ("off", "-1", -30)])
    for sa, sb in [("+", "+"), ("+", "-"), ("-", "+"), ("-", "-")]:
        add(ALL[:3], "zero-zero", [("z", sa), ("z", sb)])
        add(ALL, "zero-x", [("z", sa), ("pow2", -3, sb)])
        add(("add", "sub", "mul"), "x-zero", [("pow2", 5, sa), ("z", sb)])
    add(("div",), "div-by-zero", [d("2"), ("z", "-")])
    add(("div",), "div-by-zero", [("z", "+"), ("z", "-")])
    # subnormal operands and results
    add(ALL, "subnormal", [("sub", 1), ("sub", 3, "-")])
    add(ALL, "subnormal", [("sub", 5), d("0.75")])
    add(("add", "sub"), "subnormal", [("pow2", EMIN - 1), ("sub", 1, "-")])
    add(("mul", "div"), "underflow", [("pow2", EMIN + 2), ("pow2", -5, "-")])
    add(("add", "sub"), "subnormal-x-normal", [("sub", 1), d("1")])
    # near EXPONENT_MAX: overflow is directional
    add(ALL, "near-max", [("max",), ("max",)])
    add(ALL, "near-max", [("max", "-"), ("pow2", E_TOP - 70, "-")])
    add(("mul", "div"), "overflow", [("pow2", E_TOP), ("pow2", 3, "-")])
    add(("div",), "overflow", [("max",), ("pow2", -2)])
    add(("add", "sub"), "huge-gap", [("pow2", E_TOP), ("pow2", EMIN + 2, "-")])
    add(("add", "sub"), "huge-gap", [("sub", 1, "-"), ("max",)])
    # ties at p bits
    add(("add", "sub"), "tie", [("tie",), ("z", "+")], pa_extra=1)
    add(("mul", "div"), "tie", [("tie3", "-"), d("1")], pa_extra=1)
    # aliasing: the same &ExactNum as both operands
    for a in [d("1.5"), d("-0.1"), ("z", "-"), ("sub", 3), ("max",), ("max", "-"), ("pow2", EMIN + 1)]:
        add(ALL, "alias", [a, a], alias=True)


def main() -> None:
    build_cases()
    how_n = {"ziv": 0, "exact": 0, "model": 0}
    per_fn: dict[str, int] = {}
    n = 0
    t0 = time.time()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w") as f:
        for ci, c in enumerate(CASES):
            for p in c["ps"]:
                pa = p + c["pa_extra"]
                x, y = [r13.resolve(s, p, pa) for s in c["specs"]]
                for rm in MODES:
                    exp_ = model(c["fn"], x, y, p, rm)
                    how = mp_check(c["fn"], x, y, p, rm, exp_)
                    how_n[how] += 1
                    per_fn[c["fn"]] = per_fn.get(c["fn"], 0) + 1
                    row = {"id": f"borrow-{c['fn']}-{p}-{rm}-{n}", "fn": c["fn"], "family": c["family"],
                           "p": p, "pa": pa, "rm": rm, "alias": c["alias"],
                           "args": [x.enc(pa), y.enc(pa)], "expect": exp_, "how": how}
                    f.write(json.dumps(row, separators=(",", ":")) + "\n")
                    n += 1
            f.flush()
            print(f"[{time.time() - t0:7.1f}s] case {ci + 1}/{len(CASES)} {c['fn']} {c['family']} rows={n}",
                  file=sys.stderr, flush=True)
    print(f"rows={n} per_fn={per_fn} how={how_n}")


if __name__ == "__main__":
    main()
