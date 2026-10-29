#!/usr/bin/env python3
"""R15 goldens: directed real exp2 and cbrt, every mode, p 24..65536.

Writes zenith-float-num/tests/oracle/goldens/r15_exp2_cbrt.jsonl in the canon of
gen_s9_specials_goldens.py (little-endian 64-bit limbs, value = 0.1m * 2^e, i64
exponent range +-(2^62 - 1), subnormal grid 2^(EMIN - p)).

Only the gaps left by arith.jsonl (cbrt of +-8), elementary.jsonl (exp2 of small
integers at p 64/128) and r9_overflow.jsonl (exp2 range edges at p <= 1024):
non-exact arguments at every p, exact cases (exp2 of integers, cbrt of perfect
cubes of both signs), and range edges at p >= 4096.

  * "exact": the value is a known dyadic (2^k, or n * 2^k for a cube n^3 * 2^(3k));
    rounding it at p bits (round_rat, overflow and subnormal grid included) is the
    expect, and mpmath at p+512 must agree whenever its interval decides.
  * "inexact": the flag the result must carry (false only for an exact value that
    fits p bits and the range; overflow, underflow and rounding set it).
  * "ziv": mpmath at p+256 and p+512 bits; both error intervals must round to the
    same p-bit result. A row neither settles is an error (none are expected).
Tiers: fast p <= 1024, deep p in {4096, 16384}, full p = 65536.

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
_spec = importlib.util.spec_from_file_location("s9", HERE / "gen_s9_specials_goldens.py")
s9 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s9)

EMAX, EMIN, MODES, Arg = s9.EMAX, s9.EMIN, s9.MODES, s9.Arg
OUT = s9.ROOT / "zenith-float-num" / "tests" / "oracle" / "goldens" / "r15_exp2_cbrt.jsonl"
P_STD = [24, 53, 64, 113, 256, 1024, 4096]
P_BIG = [16384, 65536]
P_EDGE = [4096, 16384, 65536]
SMALL = "--small" in sys.argv


def tier(p: int) -> str:
    return "fast" if p <= 1024 else ("deep" if p < 65536 else "full")


def fits(a: Arg, p: int) -> bool:
    if a.n == 0:
        return True
    tz = (a.n & -a.n).bit_length() - 1
    if (a.n >> tz).bit_length() > p:
        return False
    top = a.ez + a.n.bit_length()
    return EMIN - p <= a.ez + tz and top <= EMAX


def build(spec, p):
    kind = spec[0]
    if kind == "i":  # sign * n * 2^ez exactly
        return Arg(spec[3] if len(spec) > 3 else "+", spec[1], spec[2])
    if kind == "max":
        return Arg("+", (1 << p) - 1, EMAX - p)
    if kind == "sub":
        return Arg(spec[2] if len(spec) > 2 else "+", spec[1], EMIN - p)
    return s9.arg_at(spec, p)


MPF = {"exp2": lambda x: mp.power(2, x), "cbrt": lambda x: mp.sign(x) * mp.cbrt(abs(x))}  # real root


def exact_value(fn, a: Arg):
    """(sign, n, ez) of the exact result, or None."""
    if fn == "exp2":
        if a.n == 0:
            return ("+", 1, 0)
        tz = (a.n & -a.n).bit_length() - 1
        if a.ez + tz < 0:
            return None
        k = (a.n << a.ez) if a.ez >= 0 else (a.n >> -a.ez)
        return ("+", 1, -k if a.sign == "-" else k)
    if a.n == 0:
        return None
    n, ez = a.n, a.ez
    while ez % 3:
        n <<= 1
        ez -= 1
    r = round(n ** (1 / 3)) if n < (1 << 900) else int(mp.cbrt(mp.mpf(n)))
    for c in (r - 1, r, r + 1):
        if c > 0 and c ** 3 == n:
            return (a.sign, c, ez // 3)
    lo, hi = 1, 1 << (n.bit_length() // 3 + 2)
    while lo < hi:
        mid = (lo + hi) // 2
        if mid ** 3 < n:
            lo = mid + 1
        else:
            hi = mid
    return (a.sign, lo, ez // 3) if lo ** 3 == n else None


def oracle(fn, a: Arg, p: int, rm: str):
    ex = exact_value(fn, a)
    got = []
    for extra in (256, 512):
        with mp.workprec(p + extra):
            v = MPF[fn](a.mpf())
        got.append(s9.interval_round(v, p + extra, p, rm) if v != 0 else None)
    if ex is not None:
        want = s9.round_rat(ex[0], ex[1], 1, ex[2], p, rm)
        if got[1] is not None and got[1] != want:
            raise AssertionError(f"{fn} {p} {rm}: exact {want} vs mpmath {got[1]}")
        lo = s9.round_rat(ex[0], ex[1], 1, ex[2], p, "ToZero")
        hi = s9.round_rat(ex[0], ex[1], 1, ex[2], p, "FromZero")
        return want, "exact", lo != hi or want["k"] != "f"
    if got[0] is not None and got[0] == got[1]:
        return got[0], "ziv", True
    raise AssertionError(f"{fn} p={p} {rm}: unresolved ({got})")


CASES = []


def add(fn, family, spec, ps):
    CASES.append((fn, family, spec, ps))


def build_cases():
    d = lambda s: ("d", s)
    i = lambda n, ez=0, s="+": ("i", n, ez, s)
    ALL = P_STD + P_BIG
    for s in ["0.5", "-0.5", "0.1", "-7.3", "1.4426950408889634", "100.25", "-1000.125", "300000.7"]:
        add("exp2", "moderate", d(s), ALL)
    for k in [0, 1, -1, 63, -64, 1000, -1074, 123457, -987655]:
        add("exp2", "integer", i(abs(k), 0, "-" if k < 0 else "+") if k else ("z", "+"), ALL)
    for k, fam in [(EMAX - 1, "edge-top"), (EMIN - 1, "edge-bottom")]:
        add("exp2", fam, i(abs(k), 0, "-" if k < 0 else "+"), P_EDGE)
    add("exp2", "edge-overflow", i(EMAX), P_EDGE)
    add("exp2", "edge-overflow", ("i", 2 * EMAX + 1, -1), P_EDGE)      # EMAX + 1/2
    add("exp2", "edge-top", ("i", 2 * EMAX - 1, -1), P_EDGE)          # EMAX - 1/2
    for p in P_EDGE:  # subnormal floor EMIN - p and half of it (a tie), EMIN - p - 1/2
        add("exp2", "edge-subnormal", i(-(EMIN - p), 0, "-"), [p])
        add("exp2", "edge-underflow", i(-(EMIN - p - 1), 0, "-"), [p])
        add("exp2", "edge-underflow", ("i", -(2 * (EMIN - p) - 1), -1, "-"), [p])
    for s in ["2", "3", "-2", "0.1", "-1e60", "1e-300", "7.5", "-0.001"]:
        add("cbrt", "moderate", d(s), ALL)
    for n, k, sg in [(27, 0, "+"), (27, 0, "-"), (1, -3, "+"), (343, -9, "-"), (12345 ** 3, 0, "-"),
                     (1000001 ** 3, 30, "+"), ((2 ** 40 + 1) ** 3, -60, "-"), (3 ** 120, 0, "+")]:
        add("cbrt", "cube", i(n, k, sg), ALL)
    E3 = EMAX - 1 - ((EMAX - 1) % 3)  # 2^E3 is a cube, largest such power below 2^EMAX
    add("cbrt", "edge-top", ("i", 1, E3 - 1), ALL)
    add("cbrt", "edge-top", ("i", 1, EMAX - 2, "-"), ALL)
    add("cbrt", "edge-top", ("max",), ALL)
    add("cbrt", "edge-top", ("i", 27, E3 - 6, "-"), ALL)
    add("cbrt", "edge-bottom", ("i", 1, EMIN - 1), ALL)
    add("cbrt", "edge-bottom", ("i", 1, (EMIN - 3) - ((EMIN - 3) % 3), "-"), ALL)
    add("cbrt", "edge-subnormal", ("sub", 1), ALL)
    add("cbrt", "edge-subnormal", ("sub", 27, "-"), ALL)


def main():
    build_cases()
    rows, how_n, skipped = [], {"exact": 0, "ziv": 0}, 0
    t0 = time.time()
    for ci, (fn, fam, spec, ps) in enumerate(CASES):
        for p in ps:
            if SMALL and p > 256:
                continue
            a = build(spec, p)
            if not fits(a, p):
                skipped += 1
                continue
            for rm in MODES:
                exp_, how, inx = oracle(fn, a, p, rm)
                how_n[how] += 1
                rows.append({"id": f"r15-{fn}-{p}-{rm}-{len(rows)}", "tier": tier(p), "fn": fn,
                             "family": fam, "p": p, "rm": rm, "args": [a.enc(p)], "expect": exp_,
                             "how": how, "inexact": inx})
        print(f"[{time.time() - t0:7.1f}s] case {ci + 1}/{len(CASES)} {fn} {fam} rows={len(rows)}",
              file=sys.stderr, flush=True)
    if not SMALL:
        with OUT.open("w") as f:
            for r in rows:
                f.write(json.dumps(r, separators=(",", ":")) + "\n")
    tiers = {}
    for r in rows:
        tiers[(r["fn"], r["tier"])] = tiers.get((r["fn"], r["tier"]), 0) + 1
    print(f"rows={len(rows)} skipped_p={skipped} how={how_n} tiers={tiers}")


if __name__ == "__main__":
    main()
