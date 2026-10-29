#!/usr/bin/env python3
"""R15b goldens: exp, exp2 and pow at p in {8, 24, 32} of integer-valued arguments whose
exponent is 33..61 (integer part past 32 bits, mantissa of one 32-bit word).

Writes zenith-float-num/tests/oracle/goldens/r15b_int_exp.jsonl in the canon of
gen_s9_specials_goldens.py. exp2(n) and pow(4, n) are exact powers of two (rounded at p bits,
range edges included); exp(n) and pow(1.5, n) are Ziv rows: mpmath at p+256 and p+1024 bits
must round to the same p-bit result.

Run with ~/.local/oracle-venv/bin/python (mpmath 1.3.0).
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("s9", HERE / "gen_s9_specials_goldens.py")
s9 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s9)

OUT = s9.ROOT / "zenith-float-num" / "tests" / "oracle" / "goldens" / "r15b_int_exp.jsonl"
PS = [8, 24, 32]
# (odd mantissa of <= 8 bits, bit length of n): exponents 33..61.
NS = [(1, 33), (1, 34), (3, 40), (181, 41), (5, 48), (255, 56), (1, 60), (3, 61), (181, 61)]
BASES = {"exp": None, "exp2": None, "pow15": s9.Arg("+", 3, -1), "pow4": s9.Arg("+", 1, 2)}


def ziv(f, p, rm):
    got = []
    for extra in (256, 1024):
        with mp.workprec(p + extra + 64):
            v = f()
        got.append(s9.interval_round(v, p + extra, p, rm))
    if got[0] is None or got[0] != got[1]:
        raise AssertionError(f"unresolved {got}")
    return got[0]


def main():
    rows = []
    for fn, base in BASES.items():
        for m, bl in NS:
            n = m << (bl - m.bit_length())
            for sign in "+-":
                a = s9.Arg(sign, n, 0)
                k = -n if sign == "-" else n
                for p in PS:
                    for rm in s9.MODES:
                        if fn == "exp2":
                            want, how = s9.round_rat("+", 1, 1, k, p, rm), "exact"
                        elif fn == "pow4":
                            want, how = s9.round_rat("+", 1, 1, 2 * k, p, rm), "exact"
                        elif fn == "exp":
                            want, how = ziv(lambda: mp.exp(a.mpf()), p, rm), "ziv"
                        else:
                            want, how = ziv(lambda: mp.power(base.mpf(), a.mpf()), p, rm), "ziv"
                        args = [a.enc(p)] if base is None else [base.enc(p), a.enc(p)]
                        rows.append({"id": f"r15b-{fn}-{p}-{rm}-{len(rows)}", "fn": fn,
                                     "family": f"int-e{bl}", "p": p, "rm": rm, "args": args,
                                     "expect": want, "how": how, "inexact": how == "ziv"})
    with OUT.open("w") as f:
        for r in rows:
            f.write(json.dumps(r, separators=(",", ":")) + "\n")
    print(f"rows={len(rows)}")


if __name__ == "__main__":
    main()
