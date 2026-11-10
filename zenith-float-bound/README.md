# zenith-bound

**Pre-release (0.0.0). Nothing is guaranteed yet.** No enclosure, rounding or
determinism property described anywhere in this repository holds until the
release that states it, together with its tests. Do not rely on this crate.

zenith-bound will be a guaranteed-bounds (ball arithmetic) crate built on
[zenith-float](https://github.com/jscarr64/zenith-float). Every result will
either carry a defensible enclosure or say explicitly that the guarantee could
not be established.

- 100% Rust, `#![forbid(unsafe_code)]`, no IEEE `f32`/`f64` in the math.
- MSRV: Rust 1.99 (frozen pin; changes only in a release that says so).
- The plan lives in [`docs/PLAN.md`](docs/PLAN.md). Contributing:
  [`CONTRIBUTING.md`](CONTRIBUTING.md). Security: [`SECURITY.md`](SECURITY.md).

## Function status

Words used in this table:

- **Guaranteed**: the true result is always inside the returned ball, or the call returns an error or `Indeterminate`. A function is listed as Guaranteed (or called "rigorous") **only after its contract tests pass in CI**.
- **Heuristic**: an estimate with no proof, under `heuristic::*` only.
- **Refused**: returns `Indeterminate` or an infinite radius by design.
- **Not yet**: not implemented, or its contract tests are not yet passing.

| Area | Functions | Planned | Status |
|---|---|---|---|
| L0 arithmetic | add, sub, neg, abs, mul, sqr, powi, div, sqrt, hypot, sum, dot | Guaranteed | Not yet |
| L0 queries | Tri comparisons, contains, overlaps, union, intersect, floor, ceil, sgn, is_unique_at | Guaranteed | Not yet |
| L0 I/O | decimal parse/print, canonical encoding | Guaranteed | Not yet |
| L1 elementary | exp, exp2, expm1, log, log2, log10, log1p, sin, cos, tan, asin, acos, atan, atan2, sinh, cosh, tanh, asinh, acosh, atanh, pow; pi, e, ln2, ln10 | Guaranteed | Not yet |
| L2 | complex balls, polynomials, series, matrices; gamma, erf/erfc, Bessel J0/J1, zeta | Guaranteed | Not yet |
| L2 | Meijer G, Fox H, Heun | Refused (v1) | Not yet |
| L3 | rigorous quadrature | Guaranteed | Not yet |
| L4 | root isolation, verified linear solve, Lohner ODE | Guaranteed | Not yet |
| heuristic | error estimates, point solvers | Heuristic | Not yet |

## Credit

Portions will be derived from [puremp](https://github.com/KarpelesLab/puremp)
by Karpelès Lab Inc., MIT License. See `NOTICE-puremp`. puremp's `f64` lattice
path is not ported.

## Test oracles are dev-only

FLINT/Arb (via python-flint), MPFR and mpmath are used only as development test
oracles. They are never dependencies and never shipped, and goldens derived
from them are dev-only test data. FLINT/Arb and MPFR are LGPL-licensed, so
zenith-bound implements algorithms from published papers only, never from
their source (see `CONTRIBUTING.md` and `NOTICE`).

## License

Licensed under either of

- Apache License, Version 2.0 (`LICENSE-APACHE`)
- MIT license (`LICENSE-MIT`)

at your option. Files derived from puremp also keep puremp's MIT notice.
