### Added: `Consts::get_exact_prec`, `Consts::is_exact`, `ConstantId`

- `get_exact_prec(&mut self, c: ConstantId, p: usize, rm: RoundingMode) -> Result<ExactNum, Error>`
  returns π, e, ln 2, ln 10, √2, φ or γ correctly rounded at exactly `p` bits in every
  mode: `Down ≤ c ≤ Up` for every `p ≥ 1` (1, 2, 63, 64, 65 and 65536 tested), both
  inexact flags set, bit-identical whatever the cache state or request order.
- Method: a Ziv bracket around a cached approximation trusted only to 2^64 of its ulps,
  both ends rounded at `p`; the working precision doubles until they agree.
- `is_exact(c: ConstantId) -> bool` is false for all seven. `p = 0`, `p > MAX_PRECISION`
  and `RoundingMode::None` are `InvalidArgument`.
