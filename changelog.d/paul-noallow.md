### Removed (Jeff 2026-11-24 1:08/1:17 PM CT: zero #[allow] in zenith-float, paul)

- zenith-float-num: the historical public constants `CHEBYSHEV_MAX_DEGREE`, `JACOBI_AGM_MAX`, `QUADRATURE_MAX_NODES`
  and `TANH_SINH_LEVELS_MAX`. None of them limits anything since the 2.0 cap removals (Chebyshev degree, Jacobi AGM
  steps, Gauss node count and tanh-sinh levels are bounded by the error bound or allocatable precision); they were kept
  only under `#[allow(dead_code)]`. Also removed the unused private caps (HYPERGEOM_TRANSFORM_MAX, GK_NEWTON_MAX,
  SVD_ITER_MAX, EIGEN_ITER_MAX, LOG_EXACT_MAX, QUADRATURE_NEWTON_MAX, QUADRATURE_BISECT_STEPS) and helpers
  (series_term_budget, series_term_cap; series_n_max moved into the one regression test that shows the old cap).

### Changed (paul, same ruling)

- Tests: the wall-clock budget helper `test_time` moved from `zenith-float-num/src/test_time.rs` (included by six test
  binaries through `#[path]` under `#[allow(dead_code)]`) into the unpublished dev-only workspace crate
  `zenith-float-testkit`. Its 5 unit tests now run once there instead of once per including binary (num lib 481 -> 476,
  tier1 117 -> 112, c2f1 7 -> 2, oracle 12 -> 7, r9 20 -> 15, s9 17 -> 12, specials_oracle 27 -> 22; every other test
  unchanged). check_zero_test_targets: 70 -> 71 workspace targets.

### Changed — public API (Jeff 2026-11-24 2:35 PM CT: no public fn over clippy's default 7 arguments)

- `integrate_adaptive_gk`, `integrate_adaptive_gk_estimate`, `integrate_adaptive_gk_infinite` and
  `integrate_adaptive_gk_infinite_estimate` take `GkOptions { tol_rel, max_subintervals }` in place of the two
  arguments: `integrate_adaptive_gk(f, a, b, GkOptions { tol_rel: &tol, max_subintervals: 64 }, p, rm, cc)`.
  Same arithmetic and results.
- `rk4(f, OdeSpan { t0, y0, t1 }, n_steps, p, rm, cc)` and
  `rk45_adaptive(f, OdeSpan { t0, y0, t1 }, OdeTol { atol, rtol }, p, rm, cc)`: the interval/initial value and the
  tolerances are grouped. Same arithmetic and results.
- `newton(f, df, x0, NewtonOptions { tol, max_iter }, p, rm, cc)`: the stop tolerance and iteration cap are grouped.
  Same arithmetic and results.
- `ExactNum::meijer_g(MeijerParams { an, ap, bm, bq }, p, rm, cc)`, `ExactNum::fox_h(FoxHParams { an, ap, bm, bq }, p, rm, cc)`,
  `ExactNum::heun_g(HeunGParams { a, q, alpha, beta, gamma, delta }, p, rm, cc)` and
  `ExactNum::heun_c(HeunCParams { alpha, gamma, delta, eps, q }, p, rm, cc)`: the parameter groups are one struct, fields
  in the old argument order. Same arithmetic and results.
