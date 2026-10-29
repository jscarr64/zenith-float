### Fixed

- `RoundingMode::None` at the entry points that round the requested precision themselves (`ldexp`,
  `nth_root`, `cbrt`, `sqrt`, `reciprocal`, the elementary and special functions, the constants, ...)
  now keeps `rp(p, 64)` bits on every target, like the core operations. On 32-bit targets these
  truncated at `rp(p, 32)` bits (a 224-bit `nth_root` at `p = 200`, an inexact `ldexp(x, 0, 1, None)`
  of a 64-bit `x`). No value changes on 64-bit targets.
