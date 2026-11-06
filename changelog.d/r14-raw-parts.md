### Breaking: `RawParts` replaces the raw tuple API

- `ExactNum::as_raw_parts(&self) -> RawParts<'_>` now returns a struct that borrows the
  mantissa words (no copy) and covers every value: `kind: RawKind { Finite, Zero, Inf, NaN }`,
  `words` (little-endian, most significant word last), `precision` (= `words.len() * 64`),
  `negative`, `exponent`, `inexact`. Value = (-1)^negative · M · 2^(exponent − precision).
- `ExactNum::from_raw_parts(&RawParts) -> Result<ExactNum, Error>` validates its input and
  returns `InvalidArgument` for: a non-normalized mantissa (top bit clear unless subnormal at
  `EXPONENT_MIN`), precision 0, a word count not matching the precision, an exponent outside
  `[EXPONENT_MIN, EXPONENT_MAX]` (e.g. ±2^62), and malformed Zero/Inf/NaN parts.
- Inf/NaN report empty words, precision 0, exponent 0; NaN's error kind is not carried.
- The `fp!` macro and internal callers were migrated.
