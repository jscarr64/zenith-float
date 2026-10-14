# zenith-float-num

Numeric kernel of [zenith-float](https://github.com/jscarr64/zenith-float). Depend on `zenith-float`; you do not normally need this crate directly. The text below is the zenith-float README.

> **2.0.0 is a breaking release.** It raises the MSRV to Rust 1.99, makes
> `parse` consume the whole string, tags invalid operations
> `NaN(InvalidArgument)`, changes `GkResult::abs_err_est`, adds angle-form
> elliptic functions, and removes every platform intrinsic and `unsafe` block
> from the library. Arithmetic behaviour changes too: `add`, `sub`, `mul`, and
> `div` round at exactly the requested bit precision, an overflow rounded toward
> zero returns the largest finite value, `inexact()` describes the last operation
> only (`inexact_accumulated()` keeps the old meaning), `from_raw_parts` validates
> its parts, and the exponent range is `i32::MIN..=i32::MAX` on 32-bit targets too.
> See [CHANGELOG.md](https://github.com/jscarr64/zenith-float/blob/main/CHANGELOG.md) for each change and its
> migration. 2.0.0 is not yet published; the current crates.io release is 1.0.12.

Arbitrary-precision software floating-point numbers in Rust, plus software IEEE-754 binary32/binary64 (`Ieee32` / `Ieee64`) and 1-D arrays. Version **2.0.0**. First stable release: **1.0.0** on [crates.io](https://crates.io/crates/zenith-float).

All arithmetic runs on integer limbs. The library does not use hardware floating-point registers for calculations. Construct `ExactNum` from integers or from binary, octal, decimal, or hexadecimal strings; construct IEEE widths from integer bit patterns (`from_bits`).

The library can work without `std` if a memory allocator is available.

MSRV is **Rust 1.99**: `rust-version = "1.99"` in every crate, and `rust-toolchain.toml` pins channel `1.99` with clippy and rustfmt.

The library crates forbid `unsafe_code` and use no platform intrinsics, SIMD, inline assembly, or `f32`/`f64`; `scripts/ci_no_intrinsics.sh` enforces this.

Besides the usual `+ − × ÷` and the elementary functions (`sqrt`, `exp`, `exp2`, `exp10`, `ln`, trig, hyperbolic), the public API includes `hypot`, `atan2`, `log1p`, `expm1`, `rem_pi`, specials (`erf`, `gamma`, `bessel_j`, Jacobi `sn`/`cn`/`dn`), IEEE split (`frexp`, `ldexp`, `logb`), extra constants (√2, φ, γ), `cexpr!` for complex expressions, and `SharedConsts` for sharing a constant cache across threads.

License: MIT OR Apache-2.0. See [CONTRIBUTING.md](https://github.com/jscarr64/zenith-float/blob/main/CONTRIBUTING.md) if you want to report a bug or send a patch.

## Credits

Parts of this kernel are derived from [astro-float](https://github.com/stencillogic/astro-float) (`astro-float-num` 0.3.7), Copyright (c) 2022 stencillogic, MIT License. The MIT copyright and permission notice is in [NOTICE](NOTICE). Original work is under Apache-2.0 and is also offered as MIT. MIT-derived code stays under the MIT terms, with that notice preserved. The file list is in the `zenith-float` package at `doc/PROVENANCE.md`.

## Acknowledgments

zenith-float was developed with the help of AI assistants from SpaceXAI, Anthropic, OpenAI, Google, and Meta.

## Crate layout

- `zenith-float` — public crate: `ExactNum`, `Ieee32`/`Ieee64`, arrays, `expr!`, `cexpr!`, constants, rounding.
- `zenith-float-num` — numeric kernel (you normally depend on `zenith-float` only).
- `zenith-float-macro` — `expr!` / `cexpr!` procedural macros.

## Documentation

- [Capabilities](https://github.com/jscarr64/zenith-float/blob/main/doc/ZENITH_FLOAT_CAPABILITIES.md) — what the crate does and does not do
- [Getting started](https://github.com/jscarr64/zenith-float/blob/main/doc/GETTING_STARTED.md)
- [Help](https://github.com/jscarr64/zenith-float/blob/main/doc/HELP.md) — recipes, mistakes, FAQ
- [Library](https://github.com/jscarr64/zenith-float/blob/main/doc/LIBRARY.md) — public API
- [`expr!` rounding](https://github.com/jscarr64/zenith-float/blob/main/doc/EXPR.md)
- [Precision](https://github.com/jscarr64/zenith-float/blob/main/doc/PRECISION.md)
- [Reproducibility](https://github.com/jscarr64/zenith-float/blob/main/doc/REPRODUCIBILITY.md)
- [Contributing](https://github.com/jscarr64/zenith-float/blob/main/CONTRIBUTING.md)
- rustdoc: <https://docs.rs/zenith-float>

## Features

| Feature | Default | Purpose |
| --- | --- | --- |
| `std` | yes | Formatting, `FromStr`, serde (when enabled) |
| `random` | no | `ExactNum::random_normal` for tests and fuzzing |
| `serde` | no | Serialize/deserialize as a decimal string or integer (`std` required) |
| `mpfr-tests` | no | Optional MPFR comparison tests (Linux x86_64, needs `rug`) |

`std` is on by default (formatting and `FromStr`):

```toml
[dependencies]
zenith-float = "2.0"
```

`no_std` (allocator required):

```toml
[dependencies]
zenith-float = { version = "2.0", default-features = false }
```

`no_std` is tested: `scripts/ci_nostd.sh` runs the mpmath reference cases in a `#![no_std]` harness (`nostd-tests/`) on a libc-only Linux binary and bare-metal on `thumbv7em-none-eabihf` under QEMU.

## Known issues

What 2.0.0 guarantees, and what it does not:

- Rounding is done in software on integer limbs, so the same inputs, precision, and rounding mode give bit-identical results on every platform.
- The Gauss–Kronrod `abs_err_est` is a conservative error estimate, validated against oracles. It is not a guaranteed bound.
- Guaranteed error bounds are not part of 2.0.0. Ball arithmetic that provides them is planned in a separate crate.

No other open issues. The items listed in the 1.0.5 audit (`doc/AUDIT.md`) were closed in 1.0.6 and 1.0.11. A Heun argument whose real segment from 0 meets a singular point is `NaN(InvalidArgument)`; that is the domain of the real solution, not an open defect.

## Rounding

- `RoundingMode` has seven modes. `ToEven` and `ToOdd` round to nearest and break ties to the even or odd neighbour. `Up`, `Down`, `ToZero`, and `FromZero` are directed: toward +∞, toward −∞, toward zero, and away from zero, with no tie rule. `None` skips the final rounding.
- `add`, `sub`, `mul`, and `div` round at exactly `p` bits, for any `p`; `ExactNum::add(&b, 24, RoundingMode::ToEven)` is a correctly rounded binary32-precision sum, subnormals included. The storage is still whole 64-bit words (32-bit on 32-bit targets), and the bits below `p` are zero.
- The other operations (`sqrt`, `cbrt`, `hypot`, `fma`, `rem`, parsing, and the transcendental and special functions) round to `p` rounded up to a whole number of words. Call `set_precision(p, rm)` afterwards for exactly `p` bits; that is a second rounding.
- `RoundingMode::None` skips that final rounding and keeps the storage precision, with possible extra bits.
- `inexact()` is true when the last `add`, `sub`, `mul`, or `div` rounded. `inexact_accumulated()` is true when anything upstream rounded. Other functions report the accumulated flag through both.
- `expr!` raises working precision to compensate cancellation. It does **not** itself perform correct rounding; the completion of a rounding loop in finite time depends on the expression.
- Overflow of the exponent range `i32::MIN..=i32::MAX` (the same on every target) becomes `±Inf`, except that `add`, `sub`, `mul`, `div`, `fma`, `pow`, `powi`, `hypot`, `reciprocal`, `ldexp`, `exp`, `exp2`, `exp10`, `sinh`, and `cosh` return `±max_value(p)` when the mode rounds toward zero for the result's sign (`ToZero`, `Down` for positive, `Up` for negative). Allocation and invalid arguments become `NaN`; `ExactNum::err()` returns the associated `Error`.
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

From the [repository](https://github.com/jscarr64/zenith-float) root:

```bash
bash scripts/ci.sh
```

Debug tests must finish in 10 minutes and the MPFR gate in 30 minutes (`CI_DEBUG_SECS` / `CI_MPFR_SECS`). Replay a random failure with `ZENITH_TEST_SEED=<seed> cargo test <name> -- --test-threads=1` (the seed is printed on panic).

That runs the default workspace tests, `std`-only tests, and tests with `random` and `serde` enabled. On Linux x86_64 the script also runs MPFR goldens in release:

```bash
cargo test -p zenith-float-num --features mpfr-tests -- --test-threads=1
```

Benchmarks (Criterion, release profile, deterministic fixtures):

```bash
./scripts/bench.sh              # quick smoke
./scripts/bench.sh --full       # full samples + HTML report under target/criterion/
./scripts/bench-compare.sh      # run benches, compare to doc/bench-baselines.tsv
./scripts/bench-compare.sh --update   # refresh the TSV baseline after intentional changes
cargo bench -p zenith-float-num --bench transcendentals -- ln/1024
cargo bench -p zenith-float-num --bench arithmetic -- mul_fft   # FFT-scale mul (slow)
```

Regression baselines live in `doc/bench-baselines.tsv` (tab-separated). `REGRESSION_PCT` defaults to 15.

Cross-library comparison (bigfloat-bench compatible workloads, TSV output):

```bash
./scripts/compare-bench.sh --quick    # zenith vs astro @ 132 bits, core tasks → doc/compare-results.tsv
./scripts/compare-bench.sh --full       # 132 / 1000 / 10000 bits, full task matrix
./scripts/compare-bench.sh --dashu    # include dashu-float (FBig)
```
