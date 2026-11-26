#!/usr/bin/env python3
"""Red test (official be63db9 oracle-deep, 2026-11-17): pack/round_mpf must
keep a residue far below dest ulp so directed modes can see it.

oracle-deep fail-fast killed `test elliptic` (2 failed / 1372, lines 776 and
782): elliptic_e_phi at p=512, Up and FromZero, limb mismatch.

At p=512, ToEven fl(π/2) is slightly above π/2. E(φ|1) is then
1 + ~2^{-1030}, which is above 1 and far below one ulp (2^{-511}).
Correctly rounded Up and FromZero are 1+ulp; Down, ToZero, ToEven,
ToOdd, and None stay 1. pack() collapsed the residue to 1, so the golds
matched the wrong value.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import gen_oracle_goldens as g  # noqa: E402
import gen_specials_goldens as s  # noqa: E402
import mpmath as mp  # noqa: E402

fail = 0


def bad(msg: str) -> None:
    global fail
    print(f"FAIL: {msg}")
    fail = 1


def one_plus_ulp_512() -> dict:
    # p=512, exp 1, leading bit plus the lowest limb bit set.
    return {
        "k": "f",
        "s": "+",
        "e": 1,
        "m": [
            "0000000000000001",
            "0000000000000000",
            "0000000000000000",
            "0000000000000000",
            "0000000000000000",
            "0000000000000000",
            "0000000000000000",
            "8000000000000000",
        ],
    }


def exact_one_512() -> dict:
    z = one_plus_ulp_512()
    z = dict(z)
    z["m"] = ["0000000000000000"] * 7 + ["8000000000000000"]
    return z


def main() -> int:
    # A value 1 + 2^{-1030} is the same shape as E(fl_512(π/2)|1): above 1,
    # many bits past dest ulp. Construct it at a precision that can hold it.
    mp.mp.prec = 2048
    x = 1 + mp.power(2, -1030)

    got_up = s.pack(x, 512, "Up")
    got_fz = s.pack(x, 512, "FromZero")
    got_dn = s.pack(x, 512, "Down")
    got_tz = s.pack(x, 512, "ToZero")
    got_te = s.pack(x, 512, "ToEven")
    got_exact_up = s.pack(mp.mpf(1), 512, "Up")

    want_up = one_plus_ulp_512()
    want_one = exact_one_512()
    if got_up != want_up:
        bad(f"pack(1+2^{{-1030}}, 512, Up) = {got_up}, expected 1+ulp {want_up}")
    if got_fz != want_up:
        bad(f"pack(1+2^{{-1030}}, 512, FromZero) = {got_fz}, expected 1+ulp {want_up}")
    if got_dn != want_one:
        bad(f"pack(1+2^{{-1030}}, 512, Down) = {got_dn}, expected 1")
    if got_tz != want_one:
        bad(f"pack(1+2^{{-1030}}, 512, ToZero) = {got_tz}, expected 1")
    if got_te != want_one:
        bad(f"pack(1+2^{{-1030}}, 512, ToEven) = {got_te}, expected 1")
    if got_exact_up != want_one:
        bad(f"pack(1, 512, Up) = {got_exact_up}, expected exact 1 (no false sticky)")

    # The official failing rows themselves.
    row_up = s.eval_phi("elliptic_e_phi", [s.HALF_PI, "1"], 512, "Up")
    row_fz = s.eval_phi("elliptic_e_phi", [s.HALF_PI, "1"], 512, "FromZero")
    row_dn = s.eval_phi("elliptic_e_phi", [s.HALF_PI, "1"], 512, "Down")
    if row_up != want_up:
        bad(f"eval_phi E(fl512(pi/2)|1) Up = {row_up}, expected 1+ulp")
    if row_fz != want_up:
        bad(f"eval_phi E(fl512(pi/2)|1) FromZero = {row_fz}, expected 1+ulp")
    if row_dn != want_one:
        bad(f"eval_phi E(fl512(pi/2)|1) Down = {row_dn}, expected 1")

    if fail == 0:
        print("test_round_mpf_directed_residue: ok")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
