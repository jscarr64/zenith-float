#!/usr/bin/env python3
"""Goldens for quadrature, ODE, and balls.

True values come from mpmath. The Rust harness checks
|returned − true| <= the error the call reports.
"""

import json
from pathlib import Path

import mpmath as mp

mp.mp.dps = 150

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "zenith-float-num" / "tests" / "oracle" / "goldens"


def dec(x):
    return mp.nstr(mp.mpf(x), 120, strip_zeros=False)


def row(**fields):
    return fields


# ε-extrapolation candidates. Names match `sample()` in methods_oracle.rs.
# Truth is evaluated at 600 digits and printed with 200.
EPS_CASES = [
    ("xpow_m0p9", lambda: 1 / (1 - mp.mpf("0.9"))),
    ("log_invsqrt", lambda: mp.mpf(-4)),
    ("xpow_m0p99", lambda: 1 / (1 - mp.mpf("0.99"))),
    ("log2_invsqrt", lambda: mp.mpf(16)),
    ("log3_invsqrt", lambda: mp.mpf(-96)),
    ("xpow_m0p75", lambda: mp.mpf(4)),
]

# Must equal CHOSEN in zenith-float-num/tests/gk_eps_probe.rs.
EPS_CHOSEN = [
    "xpow_m0p9",
    "log_invsqrt",
    "xpow_m0p99",
    "log2_invsqrt",
    "log3_invsqrt",
    "xpow_m0p75",
]


def eps_truth(make):
    with mp.workdps(600):
        return mp.nstr(make(), 200, strip_zeros=False)


def quad_rows():
    e = mp.e
    pi = mp.pi
    cases = [
        ("exp", "0", "1", "1e-8", 16, e - 1, "smooth"),
        ("sin", "0", "1", "1e-6", 16, 1 - mp.cos(1), "smooth"),
        ("sin", "0", dec(pi), "1e-6", 8, mp.mpf(2), "smooth"),
        ("one", "0", "1", "1e-8", 4, mp.mpf(1), "smooth"),
        ("sin40", "0", "1", "1e-4", 32, (1 - mp.cos(40)) / 40, "oscillatory"),
        ("invsqrt", "0", "1", "1e-6", 32, mp.mpf(2), "endpoint"),
        ("sq", "-1", "1", "1e-8", 4, mp.mpf(2) / 3, "smooth"),
    ]
    rows = []
    fast = [(53, "ToEven"), (53, "Up"), (113, "ToEven")]
    deep = [(64, "ToEven"), (128, "ToEven"), (256, "ToEven")]
    for p, rm in fast + deep:
        tier = "fast" if (p, rm) in fast else "deep"
        for integ, a, b, tol, panels, truth, shape in cases:
            rows.append(
                row(
                    id=f"gk-{integ}-{a}-{b}-{p}-{rm}",
                    tier=tier,
                    op="gk_adaptive",
                    shape=shape,
                    p=p,
                    rm=rm,
                    integrand=integ,
                    a=a,
                    b=b,
                    tol=tol,
                    panels=panels,
                    true=dec(truth),
                )
            )
        for integ, a, b, truth, shape in (
            ("exp", "0", "1", e - 1, "smooth"),
            ("sin", "0", "1", 1 - mp.cos(1), "smooth"),
        ):
            rows.append(
                row(
                    id=f"gk1-{integ}-{p}-{rm}",
                    tier=tier,
                    op="gk_interval",
                    shape=shape,
                    p=p,
                    rm=rm,
                    integrand=integ,
                    a=a,
                    b=b,
                    true=dec(truth),
                )
            )
        rows.append(
            row(
                id=f"gk-inv-{p}-{rm}",
                tier=tier,
                op="gk_adaptive",
                shape="singular",
                p=p,
                rm=rm,
                integrand="inv",
                a="0",
                b="1",
                tol="1e-6",
                panels=8,
                status="NotConverged",
            )
        )
        rows.append(
            row(
                id=f"gl-sq-{p}-{rm}",
                tier=tier,
                op="gauss_legendre",
                shape="smooth",
                p=p,
                rm=rm,
                integrand="sq",
                a="-1",
                b="1",
                n=8,
                true=dec(mp.mpf(2) / 3),
            )
        )
        rows.append(
            row(
                id=f"ts-exp-{p}-{rm}",
                tier=tier,
                op="tanh_sinh",
                shape="endpoint",
                p=p,
                rm=rm,
                integrand="invsqrt",
                a="0",
                b="1",
                true=dec(2),
            )
        )
        rows.append(
            row(
                id=f"laguerre-1-{p}-{rm}",
                tier=tier,
                op="gauss_laguerre",
                shape="infinite",
                p=p,
                rm=rm,
                integrand="one",
                n=12,
                true=dec(1),
            )
        )
        rows.append(
            row(
                id=f"hermite-1-{p}-{rm}",
                tier=tier,
                op="gauss_hermite",
                shape="infinite",
                p=p,
                rm=rm,
                integrand="one",
                n=12,
                true=dec(mp.sqrt(pi)),
            )
        )
        rows.append(
            row(
                id=f"wynn-half-{p}-{rm}",
                tier=tier,
                op="wynn",
                shape="epsilon",
                p=p,
                rm=rm,
                true=dec(2),
            )
        )
        # Infinite-range (7, 15). The rational map is part of the row.
        for kind, integ, a, b, truth in (
            ("pos", "exp_neg", "0", "inf", mp.mpf(1)),
            ("alg", "invsq", "1", "inf", mp.mpf(1)),
            ("bi", "gauss", "-inf", "inf", mp.sqrt(pi)),
        ):
            rows.append(
                row(
                    id=f"gk-inf-{kind}-{p}-{rm}",
                    tier=tier,
                    op="gk_infinite",
                    shape="infinite",
                    p=p,
                    rm=rm,
                    integrand=integ,
                    a=a,
                    b=b,
                    tol="1e-8",
                    est_tol="1e-8",
                    panels=24,
                    true=dec(truth),
                    note="infinite-range GK",
                )
            )
        # Endpoint singularity. The ε table supplies the estimate. The row
        # checks that estimate against the truth, and that abs_err_est covers it.
        rows.append(
            row(
                id=f"gk-eps-invsqrt-{p}-{rm}",
                tier=tier,
                op="gk_adaptive",
                shape="epsilon",
                p=p,
                rm=rm,
                integrand="invsqrt",
                a="0",
                b="1",
                tol="1e-8",
                est_tol="1e-6",
                panels=32,
                true=dec(2),
                note="epsilon-extrapolation GK",
            )
        )
    return rows


def ode_rows():
    rows = []
    specs = [
        ("euler", "exp", 100, None, None, "5e-2"),
        ("rk4", "exp", 40, None, None, "1e-5"),
        ("rk45", "exp", None, "1e-8", "1e-6", None),
        ("euler", "t", 100, None, None, "5e-2"),
        ("rk4", "t", 20, None, None, "1e-6"),
        ("rk45", "t", None, "1e-8", "1e-6", None),
    ]
    fast_p = [53, 128]
    deep_p = [256, 1024]
    for p in fast_p + deep_p:
        tier = "fast" if p in fast_p else "deep"
        for solver, problem, steps, atol, rtol, tol in specs:
            truth = mp.e if problem == "exp" else mp.mpf("0.5")
            item = row(
                id=f"{solver}-{problem}-{p}",
                tier=tier,
                op=solver,
                problem=problem,
                p=p,
                rm="ToEven",
                t0="0",
                t1="1",
                y0="1" if problem == "exp" else "0",
                true=dec(truth),
            )
            if steps is not None:
                item["steps"] = steps
                item["tol"] = tol
            else:
                item["atol"] = atol
                item["rtol"] = rtol
            rows.append(item)
    return rows


def ball_targets(op, mid, rad, other, other_rad):
    """Points the image of the input ball must contain."""
    a = mp.mpf(mid)
    r = mp.mpf(rad)
    b = mp.mpf(other)
    s = mp.mpf(other_rad)
    if op == "add":
        vals = [a + b, (a - r) + (b - s), (a + r) + (b + s)]
    elif op == "mul":
        ends = [a - r, a + r]
        oends = [b - s, b + s]
        vals = [u * v for u in ends for v in oends]
        vals.append(a * b)
    elif op == "exp":
        vals = [mp.exp(a - r), mp.exp(a), mp.exp(a + r)]
    elif op == "sin":
        vals = [mp.sin(a - r), mp.sin(a), mp.sin(a + r)]
    else:
        raise SystemExit(op)
    return [dec(v) for v in vals]


def ball_rows():
    rows = []
    ops = ("add", "mul", "exp", "sin")
    layouts = (
        ("point", "0.5", "0", "0.25", "0"),
        ("interval", "0.5", "0.05", "0.25", "0.02"),
    )
    fast_p = [53, 128]
    deep_p = [256]
    for p in fast_p + deep_p:
        tier = "fast" if p in fast_p else "deep"
        for op in ops:
            for where, mid, rad, other, other_rad in layouts:
                rows.append(
                    row(
                        id=f"ball-{op}-{where}-{p}",
                        tier=tier,
                        op=f"ball_{op}",
                        where=where,
                        p=p,
                        rm="ToEven",
                        mid=mid,
                        rad=rad,
                        other=other,
                        other_rad=other_rad,
                        true=ball_targets(op, mid, rad, other, other_rad),
                    )
                )
    return rows


def dump(name, rows):
    path = OUT / name
    text = "\n".join(json.dumps(r, separators=(",", ":")) for r in rows) + "\n"
    path.write_text(text)
    print(f"{name}: {len(rows)}")


def main():
    dump("quad.jsonl", quad_rows())
    dump("ode.jsonl", ode_rows())
    dump("balls.jsonl", ball_rows())


if __name__ == "__main__":
    main()
