### R14 review fixes

- `RawParts::significant_bits` now returns `Result<usize, Error>`: lossless `leading_zeros`
  conversion and checked arithmetic (`InvalidArgument` on overflow) instead of a cast.
- `from_raw_parts` / `value_parts` match every `RawKind` / flavor explicitly.
- `float!` literals: raw parts are validated at compile time (compile error if rejected); the
  runtime rebuild no longer falls back to a bare `NAN` but yields a NaN carrying the error.
- Docs: a Zero may carry `op_inexact` (underflowed zero), matching the code.
- Test: `get_exact_prec` with `p > MAX_PRECISION` returns `InvalidArgument`.
