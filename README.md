# zenith-float

Arbitrary-precision software floating-point numbers in Rust, plus software IEEE-754 binary32/binary64 (`Ieee32` / `Ieee64`) and 1-D arrays.

All arithmetic runs on integer limbs. The library does not use hardware floating-point registers for calculations. Construct `ExactNum` from integers or from binary, octal, decimal, or hexadecimal strings; construct IEEE widths from integer bit patterns (`from_bits`).

The library can work without `std` if a memory allocator is available.

Besides the usual `+ − × ÷` and the elementary functions (`sqrt`, `exp`, `exp2`, `exp10`, `ln`, trig, hyperbolic), the public API includes `hypot`, `atan2`, `log1p`, `expm1`, `rem_pi`, specials (`erf`, `gamma`, `bessel_j`), IEEE split (`frexp`, `ldexp`, `logb`), extra constants (√2, φ, γ), `cexpr!` for complex expressions, and `SharedConsts` for sharing a constant cache across threads.

License: MIT OR Apache-2.0.

## Crate layout

- `zenith-float` — public crate: `ExactNum`, `Ieee32`/`Ieee64`, arrays, `expr!`, `cexpr!`, constants, rounding.
- `zenith-float-num` — numeric kernel (you normally depend on `zenith-float` only).
- `zenith-float-macro` — `expr!` / `cexpr!` procedural macros.

## Documentation

- [Getting started](doc/GETTING_STARTED.md) — construct `ExactNum`, `Consts`, `expr!`, format, scoped rounding.
- [Help](doc/HELP.md) — longer tutorial: precision, rounding, recipes, mistakes, what is not in the crate.
- [Library inventory](doc/LIBRARY.md) — complete public API: types, every `ExactNum` / `ExactComplex` method, macros, constants, features, rounding, I/O.
- [Error bounds](doc/README.md) — ulp and series error theory used by `expr!` and property tests.
- [`expr!` rounding contract](doc/EXPR.md) — per-op working precision and final `set_precision`.
- [Reproducibility](doc/REPRODUCIBILITY.md) — what determines a result; how to replay; citation line.

## Features

| Feature | Default | Purpose |
| --- | --- | --- |
| `std` | yes | Formatting, `FromStr`, serde (when enabled) |
| `random` | no | `ExactNum::random_normal` for tests and fuzzing |
| `serde` | no | Serialize/deserialize as a decimal string or integer (`std` required) |
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
