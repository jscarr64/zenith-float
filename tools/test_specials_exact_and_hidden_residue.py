#!/usr/bin/env python3
"""Red test (gold regen at ceef8ae, 2026-11-17): the fixed packer plus
sufficient dps must not invent a directed residue on a mathematically
exact value, and must not flush a genuine component past dest ulp.

A high-dps regen of every pack/round_mpf gold found 7 expect diffs.

Barnes G(1)=G(2)=1 and G(4)=2 are exact factorial products. mp.barnesg
leaves a residual that shrinks with precision (sign-changing, ~2^{-prec}).
At dest p=53, Up would then false-sticky to 1+ulp / 2+ulp. The library
already uses the exact product and returns 1, 1, 2 in every mode. The
committed golds are right; the generator must emit the exact integer.

c.bessel_j on a real argument with integer order is real. A complex
mpmath path at dest p+2048 bits left a ghost imag at ~2^{-2401}. The
library and the committed golds have imag 0.

c.bessel_k(-500+0i, 1) is the branch value -K_1(500) - i π I_1(500)
(DLMF 10.34.2, ν=1, upper side). The real part is ~2^{-725}, representable
at p=128 and p=256. Official golds had flushed it to +0 because mpmath
itself returns re=0 below ~2000 bits and the old packer would have
dropped it. The library already matches the connection formula.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import gen_oracle_goldens as g  # noqa: E402
import gen_specials_goldens as s  # noqa: E402
import mpmath as mp  # noqa: E402

GOLD = ROOT / "zenith-float-num" / "tests" / "oracle" / "goldens" / "specials"

fail = 0


def bad(msg: str) -> None:
    global fail
    print(f"FAIL: {msg}")
    fail = 1


def work_prec(p: int) -> int:
    return g.round_p(p) + 2048


def exact_one_53() -> dict:
    return {"k": "f", "s": "+", "e": 1, "m": ["8000000000000000"]}


def exact_two_53() -> dict:
    return {"k": "f", "s": "+", "e": 2, "m": ["8000000000000000"]}


def one_plus_ulp_53() -> dict:
    return {"k": "f", "s": "+", "e": 1, "m": ["8000000000000800"]}


def two_plus_ulp_53() -> dict:
    return {"k": "f", "s": "+", "e": 2, "m": ["8000000000000800"]}


# Library / connection-formula pack of K_1(-500+0i) at ToEven (ceef8ae probe).
K_128 = {
    "k": "c",
    "re": {"k": "f", "s": "-", "e": -725, "m": ["37fa5994f40746e6", "b491fbfaf7bfa68f"]},
    "im": {"k": "f", "s": "-", "e": 718, "m": ["3fc4e4f30e9ffb70", "91f23cf96fd5252a"]},
}
K_256 = {
    "k": "c",
    "re": {
        "k": "f",
        "s": "-",
        "e": -725,
        "m": [
            "d7dea1dbf0bf506e",
            "36b8bd963224bcd6",
            "37fa5994f40746e6",
            "b491fbfaf7bfa68f",
        ],
    },
    "im": {
        "k": "f",
        "s": "-",
        "e": 718,
        "m": [
            "3a0f827c80f9983a",
            "e6f8f97407adf027",
            "3fc4e4f30e9ffb6f",
            "91f23cf96fd5252a",
        ],
    },
}


def gold_row(path: Path, row_id: str) -> dict:
    for line in path.read_text().splitlines():
        row = json.loads(line)
        if row["id"] == row_id:
            return row
    raise SystemExit(f"missing gold row {row_id} in {path}")


def eval_at(op: str, args: list[str], p: int, rm: str, extra: dict | None = None) -> dict:
    mp.mp.prec = work_prec(p)
    return s.eval_row(op, args, p, rm, extra or {})


def is_zero_component(v: dict) -> bool:
    if v.get("k") == "zero":
        return True
    if v.get("k") == "f" and all(limb == "0000000000000000" for limb in v.get("m", [])):
        return True
    return False


def main() -> int:
    # --- Barnes integers must stay exact, including Up at p=53. ---
    for tok, want, sticky in (
        ("1", exact_one_53(), one_plus_ulp_53()),
        ("2", exact_one_53(), one_plus_ulp_53()),
        ("4", exact_two_53(), two_plus_ulp_53()),
    ):
        got = eval_at("barnes_g", [tok], 53, "Up")
        if got == sticky:
            bad(f"barnes_g({tok}) p=53 Up false-stickied to {got}")
        elif got != want:
            bad(f"barnes_g({tok}) p=53 Up = {got}, expected exact {want}")

        gold = gold_row(GOLD / "barnes.jsonl", f"barnes_g-53-Up-" + {"1": "16", "2": "23", "4": "30"}[tok])
        if gold["expect"] != want:
            bad(f"committed barnes_g({tok}) p=53 Up gold is {gold['expect']}, expected exact {want}")

    # --- J_n(-20000) is real; no ghost imag at dest p+2048. ---
    for row_id, nu in (
        ("c.bessel_j-256-ToEven-133", "1"),
        ("c.bessel_j-256-ToEven-135", "-3"),
    ):
        got = eval_at("c.bessel_j", [f"c:-20000,0", f"c:{nu},0"], 256, "ToEven")
        if got.get("k") != "c" or not is_zero_component(got.get("im", {})):
            bad(f"{row_id} generator imag is {got.get('im')}, expected 0")
        gold = gold_row(GOLD / "complex.jsonl", row_id)
        if not is_zero_component(gold["expect"].get("im", {})):
            bad(f"{row_id} committed gold imag is {gold['expect'].get('im')}, expected 0")
        if got["re"] != gold["expect"]["re"]:
            bad(f"{row_id} generator re {got['re']} != gold {gold['expect']['re']}")

    # --- K_1(-500+0i) keeps the -K_1(500) real part. ---
    for row_id, p, want in (
        ("c.bessel_k-128-ToEven-125", 128, K_128),
        ("c.bessel_k-256-ToEven-234", 256, K_256),
    ):
        got = eval_at("c.bessel_k", ["c:-500,0", "c:1,0"], p, "ToEven")
        if got != want:
            bad(f"{row_id} generator = {got}, expected {want}")
        gold = gold_row(GOLD / "complex.jsonl", row_id)
        if gold["expect"] != want:
            bad(f"{row_id} committed gold = {gold['expect']}, expected {want}")

    if fail:
        print("test_specials_exact_and_hidden_residue: FAIL")
        return 1
    print("test_specials_exact_and_hidden_residue: ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
