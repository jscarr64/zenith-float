### Changed — public API (Jeff 2026-11-24 4:40 PM CT item 7: refusals are Result, not NaN; paul)

- `ExactNum::meijer_g`, `ExactNum::fox_h`, `ExactNum::heun_g` and `ExactNum::heun_c` now return
  `Result<ExactNum, Error>`. Every input they used to answer with a NaN (no Mellin-Barnes contour, coincident poles,
  a Heun recurrence pole, z outside the convergence disk, a non-finite argument, ...) is now
  `Err(Error::InvalidArgument)` (or the carried error of a NaN argument); a returned `Ok` is never NaN. Numeric values
  of every `Ok` result are unchanged (tier1 113/0/0, zf42/zf43 goldens unchanged). The catalog entry for Meijer G maps
  `Err` to its own NaN cell exactly as before.
