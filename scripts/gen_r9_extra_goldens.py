#!/usr/bin/env python3
"""R9 extra goldens, APPENDED to tests/oracle/goldens/r9_overflow.jsonl (existing rows untouched).

  * radix-2 parse underflow (the `ps.out_of_range() < 0 && rdx == Radix::Bin` arm of
    `parse_word`): exact digits placed on the subnormal grid and rounded once by `rm`;
  * hex and octal parse whose radix exponent is below EXPONENT_MIN (`parse_out_of_range`):
    the minimum subnormal or a signed zero by `rm`;
  * exp at |x| just below 2^63, both signs (`exp_neg_scaled`'s loop path / overflow).

Oracle: mpmath at p+256 and p+512. Parse values are exact dyadics; the two evaluations
must be the identical exact value, which is then rounded once (how="mpmath-exact"). exp
uses the Ziv pair of gen_r9_overflow_goldens.ziv (how="ziv"). Rows stream to disk;
run with ~/.local/oracle-venv/bin/python (mpmath 1.3.0). Run once: it appends.
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("r9", HERE / "gen_r9_overflow_goldens.py")
r9 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(r9)

EMIN, MODES, PS, OUT, A = r9.EMIN, r9.MODES, r9.PS, r9.OUT, r9.A
round_rat = r9.round_rat
DIG = "0123456789abcdef"


def in_radix(n: int, b: int) -> str:
    out = ""
    while True:
        n, d = divmod(n, b)
        out = DIG[d] + out
        if n == 0:
            return out


def numeral(sign: str, mant: str, w: int, rdx: int) -> str:
    """sign mant * rdx^w, exponent digits written in radix rdx."""
    sep = "_e" if rdx > 10 else "e"
    es = ("-" if w < 0 else "") + in_radix(abs(w), rdx)
    return ("-" if sign == "-" else "") + mant + sep + es


def value_of(sign: str, mant: str, w: int, rdx: int):
    ip, _, fp = mant.partition(".")
    digits = int(ip + fp, rdx)
    bits = {2: 1, 8: 3, 16: 4}[rdx]
    x = mp.ldexp(mp.mpf(digits), bits * (w - len(fp)))
    return -x if sign == "-" else x


def parse_cases(p):
    """(sign, mantissa, written exponent, radix, note); every radix-digit exponent is below EMIN."""
    out = []
    for sg in "+-":
        # radix 2: 0.digits * 2^e with e = w + int_len < EMIN
        out += [
            (sg, "1", EMIN - 2, 2, "bin parse 2^(EMIN-2) subnormal exact"),
            (sg, "1011011", EMIN - p // 2 - 6, 2, "bin parse mid subnormal"),
            (sg, "1.0110111", EMIN - p + 3, 2, "bin parse near min subnormal, frac digits"),
            (sg, "1", EMIN - p - 1, 2, "bin parse half min subnormal (tie)"),
            (sg, "11", EMIN - p - 2, 2, "bin parse 3/4 min subnormal"),
            (sg, "1", EMIN - p - 10, 2, "bin parse far below min subnormal"),
            # hex / octal: radix exponent below EXPONENT_MIN (flush)
            (sg, "1", EMIN - 2, 16, "hex parse below range, flush"),
            (sg, "f.f", EMIN - 5, 16, "hex parse f.f below range, flush"),
            (sg, "1", EMIN - 2, 8, "oct parse below range, flush"),
            (sg, "7.7", EMIN - 5, 8, "oct parse 7.7 below range, flush"),
        ]
    return out


def exact_round(x_lo, x_hi, p, rm):
    if x_lo != x_hi:
        return None
    s, man, exp, _ = x_lo._mpf_
    return round_rat("-" if s else "+", int(man), 1, int(exp), p, rm)


def main():
    t0 = time.time()
    rows = bad = 0
    with OUT.open("a") as fh:
        for p in PS:
            for sg, mant, w, rdx, note in parse_cases(p):
                s = numeral(sg, mant, w, rdx)
                with mp.workprec(p + 256):
                    lo = value_of(sg, mant, w, rdx)
                with mp.workprec(p + 512):
                    hi = value_of(sg, mant, w, rdx)
                for rm in MODES:
                    want = exact_round(lo, hi, p, rm)
                    if want is None:
                        print("DISAGREE", note, p, rm, file=sys.stderr)
                        bad += 1
                        continue
                    row = {"fn": "parse", "p": p, "rm": rm, "args": [], "expect": want,
                           "how": "mpmath-exact", "note": note, "s": s, "rdx": rdx}
                    fh.write(json.dumps(row) + "\n")
                    rows += 1
                fh.flush()
            pa = max(p, 128)
            for sg in "+-":
                for a, note in [(A(sg, (1 << p) - 1, 63 - p), "exp |x| = 2^63 - 2^(63-p)"),
                                (A(sg, (1 << 128) - 1, 63 - 128), "exp |x| = 2^63 - 2^-65")]:
                    enc = r9.enc_exact(a, pa)
                    cache: dict = {}
                    for rm in MODES:
                        z = r9.ziv("exp", [a], 0, None, p, rm, cache)
                        if z is None:
                            print("UNRESOLVED exp", note, p, rm, file=sys.stderr)
                            bad += 1
                            continue
                        row = {"fn": "exp", "p": p, "rm": rm, "args": [enc], "expect": z[0],
                               "how": "ziv", "note": note}
                        fh.write(json.dumps(row) + "\n")
                        rows += 1
                    fh.flush()
    print(f"appended rows={rows} bad={bad} secs={time.time() - t0:.1f}")
    if bad:
        sys.exit(1)


if __name__ == "__main__":
    main()
