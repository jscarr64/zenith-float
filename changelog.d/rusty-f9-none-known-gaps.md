### Mode matrices: None everywhere it passes the strict oracle, explicit NONE-GAP rows elsewhere

- The named mode matrices of exact_bit_functions, wide_primitives, retry_budget, subnormal_ops (ToOdd too), exponent_range,
  r13_props (ToOdd too), sinh_cosh_overflow and add_sub_zero_tests (ToOdd too) now hold all 7 modes. None is checked against
  MPFR toward zero at rp(p, 64) (FIXLIST ZF-M27).
- The families that fail that oracle loop over `NONE_GAP_MODES` (every mode but None). Each use carries a `NONE-GAP F-x:` label,
  and `tools/test_rounding_mode_matrices.sh` holds the KNOWN-GAP list, queued for the ZF-M27 conformance PR:
  - F-A: exact pow, nth_root and log-of-a-power results come back as their predecessor.
  - F-B: None overflow gives ±Inf; Jeff ruled ±MAX.
  - F-C: log with a base near a grid point.
