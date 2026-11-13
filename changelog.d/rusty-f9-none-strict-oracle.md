### MPFR harness: a strict reference for RoundingMode::None

- `tests/mpfr/common.rs`: `mpfr_reference_shifted` answers `None` with MPFR toward zero at the storage precision rp(p, 64)
  (FIXLIST ZF-M27), instead of an error. `storage_bits` gives rp(p, 64).
- `tests/mpfr/shifted.rs`: every Term reference maps `None` to `ToZero` at rp(p, 64), so its sticky bits and stand-in grids sit
  below the storage precision.
- None now joins the mode matrices of fuzz_round_modes, the shifted self-test, inexact_per_op, r14_raw_parts_props (ToOdd too)
  and r14_consts. In r14_consts, None is a refusal row (`Err(InvalidArgument)`) with no golds.
