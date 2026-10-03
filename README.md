# zenith-float

> **Please upgrade to zenith-float 1.0.8.** Versions 1.0.0 through 1.0.4
> return wrong values from `erf`, `erfc`, `normal_cdf`, `hypergeom_2f1`,
> `betainc`, the gamma family, Bessel and several other special functions for
> some inputs, with no error reported, and their quadrature, root-finding and
> ODE routines return `None` for intervals such as `[2, 3]`. 1.0.5 fixes these
> with no API changes:
> `cargo update -p zenith-float -p zenith-float-num -p zenith-float-macro`.
> Versions 1.0.0 through 1.0.4 have been yanked from crates.io.
> Details: [1.0.5 release notes](https://github.com/jscarr64/zenith-float/releases/tag/v1.0.5). 1.0.6 is an additive SoftFloat pack: Jacobi \(P_n^{(\alpha,\beta)}\), real large-argument Bessel Hankel, complex Bessel cut/speed fixes, catalog specials (Scorer, Kelvin, Struve, Anger–Weber, Clausen, Barnes \(G\), polygamma / Hurwitz, inverse Jacobi, \({}_pF_q\), Lambert \(W\)), plus Meijer \(G\), Fox \(H\), and local / confluent Heun. 1.0.7 is a compile-only fix: the 1.85-required unsafe wrappers around `_addcarry_u64`/`_subborrow_u64` stay, and `unused_unsafe` is allowed so the crate builds on rustc 1.93+ (those intrinsics became safe in 1.93, not 1.87). No math change.

Arbitrary-precision software floating-point numbers in Rust, plus software IEEE-754 binary32/binary64 (`Ieee32` / `Ieee64`) and 1-D arrays. Current release: **1.0.8**. First stable release: **1.0.0** on [crates.io](https://crates.io/crates/zenith-float).

All arithmetic runs on integer limbs. The library does not use hardware floating-point registers for calculations. Construct `ExactNum` from integers or from binary, octal, decimal, or hexadecimal strings; construct IEEE widths from integer bit patterns (`from_bits`).

The library can work without `std` if a memory allocator is available.

Besides the usual `+ − × ÷` and the elementary functions (`sqrt`, `exp`, `exp2`, `exp10`, `ln`, trig, hyperbolic), the public API includes `hypot`, `atan2`, `log1p`, `expm1`, `rem_pi`, specials (`erf`, `gamma`, `bessel_j`, Jacobi `sn`/`cn`/`dn`), IEEE split (`frexp`, `ldexp`, `logb`), extra constants (√2, φ, γ), `cexpr!` for complex expressions, and `SharedConsts` for sharing a constant cache across threads.

License: MIT OR Apache-2.0. See [CONTRIBUTING.md](CONTRIBUTING.md) if you want to report a bug or send a patch.

## Crate layout

- `zenith-float` — public crate: `ExactNum`, `Ieee32`/`Ieee64`, arrays, `expr!`, `cexpr!`, constants, rounding.
- `zenith-float-num` — numeric kernel (you normally depend on `zenith-float` only).
- `zenith-float-macro` — `expr!` / `cexpr!` procedural macros.

## Documentation

- [Capabilities](doc/ZENITH_FLOAT_CAPABILITIES.md) — what the crate does and does not do
- [Getting started](doc/GETTING_STARTED.md)
- [Help](doc/HELP.md) — recipes, mistakes, FAQ
- [Library](doc/LIBRARY.md) — public API
- [`expr!` rounding](doc/EXPR.md)
- [Precision](doc/PRECISION.md)
- [Reproducibility](doc/REPRODUCIBILITY.md)
- [Contributing](CONTRIBUTING.md)
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
zenith-float = "1.0"
```

`no_std` (allocator required):

```toml
[dependencies]
zenith-float = { version = "1.0", default-features = false }
```

`no_std` is tested: `scripts/ci_nostd.sh` runs the mpmath reference cases in a `#![no_std]` harness (`nostd-tests/`) on a libc-only Linux binary and bare-metal on `thumbv7em-none-eabihf` under QEMU.

## Known issues

The 1.0.5 complex-Bessel cut hang and the real large-\(|x|\) `InvalidArgument` are closed in 1.0.6. Complex `bessel_k` can still be slow for non-integer \(\nu\) when \(|z|\) is just inside the series regime (the Temme \({}_2F_0\) path covers the large-\(|z|\) cases). Meijer \(G\) returns `NaN` at coincident poles (no logarithmic residue). Local / confluent Heun are series for \(\lvert z\rvert<1\) only. See `CHANGELOG.md`.

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
