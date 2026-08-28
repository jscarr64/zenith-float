# zenith-float

Arbitrary-precision software floating-point numbers in Rust.

All arithmetic runs on integer limbs. The library does not use hardware floating-point registers for calculations. Construct numbers from integers or from binary, octal, decimal, or hexadecimal strings.

The library can work without `std` if a memory allocator is available.

Besides the usual `+ − × ÷` and the elementary functions (`sqrt`, `exp`, `exp2`, `exp10`, `ln`, trig, hyperbolic), the public API includes `hypot`, `atan2`, `log1p`, `expm1`, and `rem_pi`.

License: MIT OR Apache-2.0.

## Crate layout

- `zenith-float` — public crate: `ExactNum`, `expr!`, constants, rounding.
- `zenith-float-num` — numeric kernel (you normally depend on `zenith-float` only).
- `zenith-float-macro` — the `expr!` procedural macro.

## Documentation

- [Build checklist](doc/BUILD_CHECKLIST.md) — what is implemented, tested, and still required (including vs astro-float / dashu-float).
- [Error bounds](doc/README.md) — ulp and series error theory used by `expr!` and property tests.

## Features

| Feature | Default | Purpose |
| --- | --- | --- |
| `std` | yes | Formatting, `FromStr`, serde (when enabled) |
| `random` | no | `ExactNum::random_normal` for tests and fuzzing |
| `serde` | no | Serialize/deserialize as a decimal string (or integer) |
| `mpfr-tests` | no | Optional MPFR comparison tests (Linux x86_64, needs `rug`) |

`no_std` (allocator required):

```toml
[dependencies]
zenith-float = { version = "0.1.0", default-features = false }
```

## Rounding

- `ExactNum` methods that take a rounding mode other than `RoundingMode::None` round the result to the requested precision.
- `RoundingMode::None` skips that final rounding and may keep extra bits.
- `expr!` raises working precision to compensate cancellation. It does **not** itself perform correct rounding; the completion of a rounding loop in finite time depends on the expression.
- Overflow of the configured exponent range becomes `±Inf`. Allocation and invalid arguments become `NaN`; `ExactNum::err()` returns the associated `Error`.
- Subnormals, `NaN`, and infinities are software values, not hardware floating-point registers.

## Example

```rust
use zenith_float::Consts;
use zenith_float::RoundingMode;
use zenith_float::ctx::Context;
use zenith_float::expr;

let mut ctx = Context::new(
    1024,
    RoundingMode::ToEven,
    Consts::new().expect("Constants cache initialized"),
    -10000,
    10000,
);

let pi = expr!(6 * atan(1 / sqrt(3)), &mut ctx);
let pi_lib = ctx.const_pi();
assert_eq!(pi.cmp(&pi_lib), Some(0));
```

Set the smallest exponent range that still holds your results. Internal working precision can grow with the exponent.

## Tests

From this directory:

```bash
bash scripts/ci.sh
```

That runs the default workspace tests, `std`-only tests, and tests with `random` and `serde` enabled. MPFR goldens are not part of that script:

```bash
cargo test -p zenith-float-num --features mpfr-tests -- --test-threads=1
```
