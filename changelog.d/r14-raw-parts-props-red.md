### Tests (red): RawParts precision round trip and every-kind proptest

- `tests/r14_raw_parts_props.rs`: values built at p = 1, 63, 64, 65, 100, 128, 129 round-trip
  to the same representation and behave identically afterwards (add/mul/div/sqrt at p and
  p+1, all modes); a proptest over Finite, ±Zero, ±Inf, NaN and subnormals at `EXPONENT_MIN`
  with p in 1..=4096 checks the round trip is bit-identical, both inexact flags included.
- Red: `RawParts` does not carry the per-operation inexact flag, so a value with
  `inexact() == false` and `inexact_accumulated() == true` comes back with `inexact() == true`.
