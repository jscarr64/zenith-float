#!/usr/bin/env python3
"""Goldens for complex Γ and ψ at dest p 1920..2240 (`ExactComplex::gamma` / `digamma`).

Prints the CASES rows of zenith-float-num/tests/complex_gamma_high_p_mpmath.rs. Inputs are exact
dyadics n·2^e. Each part is settled twice: mpmath at p+256 and at p+512 bits (round to nearest);
the interval [v - ulp, v + ulp] of each must round to one p-bit result in the mode and the two
must agree, or the script fails. p is a multiple of 64 (the dest p is round_p(p)). No ToOdd: the
values are irrational, so ToOdd equals the truncation with the low bit set.

Run with ~/team-env/mpmath-1.3.0/bin/python (verify_env.sh mpmath-1.3.0 first).
"""
from __future__ import annotations

import sys
import mpmath
from mpmath import libmp as L

MODES = {"ToEven": "n", "Up": "c", "Down": "f", "ToZero": "d", "FromZero": "u"}

# (name, function, p, mode, (re_n, re_e), (im_n, im_e)); value = n * 2^e
CASES = [
    ("gamma_three_quarters_quarter_i_p1920", "gamma", 1920, "ToEven", (3, -2), (1, -2)),
    ("gamma_three_quarters_quarter_i_p1984", "gamma", 1984, "ToEven", (3, -2), (1, -2)),
    ("gamma_three_quarters_quarter_i_p2048", "gamma", 2048, "ToEven", (3, -2), (1, -2)),
    ("gamma_three_quarters_minus_quarter_i_p1984", "gamma", 1984, "ToEven", (3, -2), (-1, -2)),
    ("gamma_tiny_real_p2048_to_even", "gamma", 2048, "ToEven", (1, -42), (0, 0)),
    ("gamma_tiny_real_p2048_up", "gamma", 2048, "Up", (1, -42), (0, 0)),
    ("gamma_tiny_real_p2048_down", "gamma", 2048, "Down", (1, -42), (0, 0)),
    ("gamma_tiny_real_p2048_to_zero", "gamma", 2048, "ToZero", (1, -42), (0, 0)),
    ("gamma_tiny_real_p2048_from_zero", "gamma", 2048, "FromZero", (1, -42), (0, 0)),
    ("gamma_tiny_real_p2240", "gamma", 2240, "ToEven", (1, -42), (0, 0)),
    ("gamma_tiny_imag_p1984", "gamma", 1984, "ToEven", (0, 0), (1, -10)),
    ("gamma_two_pow_20_minus_1_p2048", "gamma", 2048, "ToEven", (1048575, 0), (1, -1)),
    ("gamma_two_pow_20_plus_1_p2048", "gamma", 2048, "ToEven", (1048577, 0), (1, -1)),
    ("gamma_two_pow_21_minus_1_p1984", "gamma", 1984, "ToEven", (2097151, 0), (1, -1)),
    ("gamma_two_pow_21_plus_1_p1984", "gamma", 1984, "ToEven", (2097153, 0), (1, -1)),
    ("digamma_three_quarters_quarter_i_p2048", "digamma", 2048, "ToEven", (3, -2), (1, -2)),
    ("digamma_tiny_real_p2048", "digamma", 2048, "ToEven", (1, -42), (0, 0)),
    ("digamma_two_pow_20_minus_1_p2048", "digamma", 2048, "ToEven", (1048575, 0), (1, -1)),
]


def settle_part(f, z, p: int, part: str, rnd: str):
    res = []
    for extra in (256, 512):
        q = p + extra
        with mpmath.workprec(q + 64):
            w = f(z)
            v = (w.real if part == "re" else w.imag)
        with mpmath.workprec(q):
            v = +v
        t = v._mpf_
        if t == L.fzero:
            res.append(None)
            continue
        sign, man, exp, bc = t
        ulp = L.from_man_exp(1, exp + bc - q)
        lo = L.mpf_pos(L.mpf_sub(t, ulp), p, rnd)
        hi = L.mpf_pos(L.mpf_add(t, ulp), p, rnd)
        if lo != hi:
            sys.exit(f"interval straddles at {p} {part} extra={extra}")
        res.append(lo)
    if res[0] != res[1]:
        sys.exit(f"settle disagreement at {p} {part}")
    return res[0]


def rust_part(t) -> str:
    if t is None:
        return '(false, "0", 0)'
    sign, man, exp, bc = t
    return f'({"true" if sign else "false"}, "{man:x}", {exp})'


def main() -> None:
    rows = []
    for name, fn, p, mode, (rn, re_), (inn, ie) in CASES:
        assert p % 64 == 0
        with mpmath.workprec(p + 1024):
            z = mpmath.mpc(mpmath.ldexp(rn, re_), mpmath.ldexp(inn, ie))
        f = mpmath.gamma if fn == "gamma" else mpmath.digamma
        r = settle_part(f, z, p, "re", MODES[mode])
        i = settle_part(f, z, p, "im", MODES[mode])
        rows.append(f'    ("{name}", "{fn}", {p}, RoundingMode::{mode}, ({rn}, {re_}), ({inn}, {ie}), {rust_part(r)}, {rust_part(i)}),')
        print(name, "ok", file=sys.stderr, flush=True)
    sys.stdout.write("\n".join(rows) + "\n")


if __name__ == "__main__":
    main()
