#!/usr/bin/env python3
"""S3 goldens: the Euler-Mascheroni constant gamma at exactly p bits (`Consts::euler_gamma`).

Writes zenith-float-num/tests/oracle/goldens/s3_euler_gamma.jsonl in the canon of
gen_s9_specials_goldens.py (little-endian 64-bit limbs at round_p(p) bits, value =
0.1m * 2^e). One row per (constant, p, mode), streamed to disk as it is computed.

Each row is settled twice: mpmath at p+256 and at p+512 bits (round to nearest), and
the error interval [v - ulp, v + ulp] of each must round to one p-bit result in the
mode; the two results must agree, or the script fails. Five modes (no ToOdd: an
irrational constant is never a tie, so it equals ToEven).

Run with ~/.local/oracle-venv/bin/python3 (mpmath 1.3.0).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from mpmath import libmp as L

WORD = 64
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "zenith-float-num" / "tests" / "oracle" / "goldens" / "s3_euler_gamma.jsonl"
PS = [24, 53, 64, 113, 256, 1024, 4096, 65536, 131072]
MODES = {"ToEven": "n", "Up": "c", "Down": "f", "ToZero": "d", "FromZero": "u"}
CONSTS = {"EulerGamma": L.mpf_euler}


def round_p(p: int) -> int:
    return ((p + WORD - 1) // WORD) * WORD


def settle(name: str, p: int, extra: int, rnd: str):
    """p-bit result in `rnd` from the value at p+extra bits, or None if the interval straddles."""
    q = p + extra
    v = CONSTS[name](q, "n")
    sign, man, exp, bc = v
    ulp = L.from_man_exp(1, exp + bc - q)  # >= one ulp at q bits
    lo = L.mpf_pos(L.mpf_sub(v, ulp, q + 8, "f"), p, rnd)
    hi = L.mpf_pos(L.mpf_add(v, ulp, q + 8, "c"), p, rnd)
    return lo if lo == hi else None


def enc(v, p: int) -> dict:
    sign, man, exp, bc = v
    prec = round_p(p)
    m = man << (prec - bc)
    mask = (1 << WORD) - 1
    return {"k": "f", "s": "-" if sign else "+", "e": exp + bc,
            "m": [f"{(m >> (WORD * i)) & mask:016x}" for i in range(prec // WORD)]}


def main() -> int:
    rows = 0
    with OUT.open("w") as f:
        for name in CONSTS:
            for p in PS:
                for mode, rnd in MODES.items():
                    a, b = settle(name, p, 256, rnd), settle(name, p, 512, rnd)
                    if a is None or b is None or a != b:
                        print(f"unsettled {name} p={p} {mode}", file=sys.stderr)
                        return 1
                    row = {"id": f"s3-{name}-{p}-{mode}", "c": name, "p": p, "rm": mode,
                           "expect": enc(a, p), "how": "mpmath-p+256-p+512"}
                    f.write(json.dumps(row, separators=(",", ":")) + "\n")
                    rows += 1
            f.flush()
            print(f"{name} done ({rows} rows)", flush=True)
    print(f"wrote {rows} rows to {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
