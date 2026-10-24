### Breaking: `RawParts` carries the per-operation inexact flag

- New field `RawParts::op_inexact` (`ExactNum::inexact()`); `inexact` stays the accumulated
  flag. `as_raw_parts`/`from_raw_parts` now round-trip both bit-identically.
- `from_raw_parts` rejects `op_inexact` on Inf/NaN and `op_inexact && !inexact`
  (`InvalidArgument`). A Zero may carry `op_inexact` (an underflowed zero is inexact).
- Precision: no requested precision is stored apart from the mantissa word count, so a value
  built at p = 100 round-trips as 128 bits and behaves identically (documented on `RawParts`).
- The binary format keeps one flag; decoded values set both from it. `float!` literals emit
  the new field.
