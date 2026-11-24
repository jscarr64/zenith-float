### Removed (Jeff 2026-11-24 1:08/1:17 PM CT: zero #[allow] in zenith-float, paul)

- zenith-float-num: the historical public constants `CHEBYSHEV_MAX_DEGREE`, `JACOBI_AGM_MAX`, `QUADRATURE_MAX_NODES`
  and `TANH_SINH_LEVELS_MAX`. None of them limits anything since the 2.0 cap removals (Chebyshev degree, Jacobi AGM
  steps, Gauss node count and tanh-sinh levels are bounded by the error bound or allocatable precision); they were kept
  only under `#[allow(dead_code)]`. Also removed the unused private caps (HYPERGEOM_TRANSFORM_MAX, GK_NEWTON_MAX,
  SVD_ITER_MAX, EIGEN_ITER_MAX, LOG_EXACT_MAX, QUADRATURE_NEWTON_MAX, QUADRATURE_BISECT_STEPS) and helpers
  (series_term_budget, series_term_cap; series_n_max moved into the one regression test that shows the old cap).
