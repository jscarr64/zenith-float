### Changed (bound layer, CBall)

- `CBall::mul`, `CBall::sqr` and `CBall::div` are now a faithful port of Arb `acb_mul` / `acb_sqr` / `acb_div` /
  `acb_inv` (FLINT 3.6.0), including the `arb_mul`, `arb_addmul`, `arb_div` (`div_arf`, `div_wide`, reduced working
  precision) and `acb_inv` radius bounds they use. An exactly-zero part no longer makes the other part all reals, an
  axis divisor gives an exact 0 where Arb does, and every radius is within 2x FLINT's on the 646 golden rows and the
  4,006 Tier 2 c_mul/c_div rows (all contain the exact image). Arb's indeterminate results are the all-reals part.
  Midpoints stay round-to-nearest-even with half-ulp error (Paul call per CONTRACT O19).
