#!/usr/bin/env python3
"""ln, log2 and log10 at the edges of the i64 exponent range.

Writes zenith-float-num/tests/oracle/goldens/log_extreme.tsv in the canon of
gen_s9_specials_goldens.py (little-endian 64-bit limbs, i64 exponents +-(2^62 - 1)).
One row per (function, operand, precision, mode): id, fn, p, rm, argument, expect, how.

The reference is mpmath at p+256 and at p+512 bits; the error interval of each must
round to one p-bit result and the two must agree ("how": "ziv"). Exact results
(ln 1, log2 2^k, log10 10^k) come from their integer value ("how": "exact").
A value within 2^-(p+240) of a p-bit number or midpoint (ln(1+2^-(p-1)) and kin) is settled by
the next pair (p+512, p+1024) and so on ("how": "ziv+512"); any pair that disagrees, or a row
nothing settles, fails the script.

Run with ~/.local/oracle-venv/bin/python (mpmath 1.3.0).
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("s9", HERE / "gen_s9_specials_goldens.py")
s9 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s9)

EMAX, EMIN, MODES = s9.EMAX, s9.EMIN, s9.MODES
OUT = s9.ROOT / "zenith-float-num" / "tests" / "oracle" / "goldens" / "log_extreme.tsv"
PS = [24, 53, 113, 256, 1024]
FNS = {"ln": mp.log, "log2": lambda t: mp.log(t, 2), "log10": mp.log10}


def operands(p):
    """33 positive p-bit operands: range ends, subnormals, near one, powers of 2 and 10."""
    d = lambda s: s9.arg_at(("d", s), p)
    pw = lambda k: s9.Arg("+", 1, k)
    out = [
        pw(EMAX - 1), pw(EMAX - 2), s9.Arg("+", (1 << p) - 1, EMAX - p), pw(EMIN), pw(EMIN + 1),
        pw(EMIN - 1), pw(EMIN - p), s9.Arg("+", 3, EMIN - p), s9.Arg("+", (1 << p) - 1, EMIN - p),
        pw(1 << 40), pw(-(1 << 40)), s9.Arg("+", 3, (1 << 61)), s9.Arg("+", 5, -(1 << 61)),
        d("1"), d("2"), d("0.5"), d("10"), d("1000"), d("1e-3"), d("0.1"), d("3"), d("7"),
        s9.arg_at(("off", "1", -(p // 2)), p), s9.arg_at(("off", "1", -(p - 1)), p),
        s9.arg_at(("off", "1", -p - 1), p),
        s9.Arg("+", (1 << p) - 1, -p), s9.Arg("+", (1 << p) - 1, -p + 1),
        d("1e300"), d("1e-300"), d("6.02214076e23"), d("1.616255e-35"),
        s9.arg_at(("mp", mp.e), p), s9.arg_at(("mp", mp.pi), p),
    ]
    assert len(out) == 33, len(out)
    return out


def exact_value(fn, a: s9.Arg):
    """The integer result when it is exact, else None."""
    n, ez = a.n, a.ez
    if n == 0:
        return None
    tz = (n & -n).bit_length() - 1
    n, ez = n >> tz, ez + tz
    if fn == "ln":
        return 0 if (n, ez) == (1, 0) else None
    if fn == "log2":
        return ez if n == 1 else None
    if ez < 0:
        return None
    # 10^k = 5^k * 2^k: the odd part is 5^ez (checked only while 5^ez can fit in n)
    if 2 * ez > n.bit_length():
        return None
    return ez if n == 5 ** ez else None


def canon(r: dict) -> str:
    if r["k"] == "inf":
        return f"inf:{r['s']}"
    return f"f:{r['s']}:{r['e']}:{','.join(r['m'])}"


def main() -> None:
    rows, how_n = [], {"ziv": 0, "exact": 0}
    for p in PS:
        for oi, a in enumerate(operands(p)):
            argc = canon(a.enc(p))
            for fn, f in FNS.items():
                ex = exact_value(fn, a)
                vals = {}
                if ex is None:
                    for extra in s9.EXTRAS:
                        with mp.workprec(p + extra):
                            vals[extra] = f(a.mpf())
                for rm in MODES:
                    if ex is not None:
                        r = s9.zero("+", p) if ex == 0 else s9.round_rat("-" if ex < 0 else "+", abs(ex), 1, 0, p, rm)
                        how = "exact"
                    else:
                        r, how = None, None
                        for i in range(len(s9.EXTRAS) - 1):
                            e0, e1 = s9.EXTRAS[i], s9.EXTRAS[i + 1]
                            g0 = s9.interval_round(vals[e0], p + e0, p, rm)
                            g1 = s9.interval_round(vals[e1], p + e1, p, rm)
                            if g0 is not None and g1 is not None and g0 != g1:
                                raise AssertionError(f"{fn} op{oi} p={p} {rm}: p+{e0} {g0} p+{e1} {g1}")
                            if g0 is not None and g0 == g1:
                                r, how = g0, "ziv" if i == 0 else f"ziv+{e0}"
                                break
                        if r is None:
                            raise AssertionError(f"{fn} op{oi} p={p} {rm}: unresolved")
                    how_n[how] = how_n.get(how, 0) + 1
                    rows.append(f"lx-{fn}-{p}-{rm}-op{oi}\t{fn}\t{p}\t{rm}\t{argc}\t{canon(r)}\t{how}")
        print(f"p={p} rows={len(rows)}", file=sys.stderr, flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(rows) + "\n")
    print(f"rows={len(rows)} how={how_n}")


if __name__ == "__main__":
    main()
