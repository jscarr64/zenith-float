# Changelog

## Unreleased

- Native `+`, `-`, `*`, `/` operators on `ExactNum` (default precision, round-to-even).
- `exact!` / `fbig!` compile-time decimal float literals.
- Parse/format for radices 2–36 (`Radix::try_new`); `_e` exponent separator for bases > 10.
- `RadixFloat` wrapper and `ExactNum::with_radix` for radix-tagged I/O.
- `ExactNum::nth_root` / `expr!` `root(x, n)` — general n-th root (`sqrt`/`cbrt` delegation, composite factors, Newton for primes).
- `ExactNum::sinh_cosh` — paired hyperbolic evaluation with a single `exp(|x|)` path.
- `ExactNum::copysign`, `ExactNum::next_after` — software IEEE sign and successor.
- `ExactNum::fma` / `expr!` `fma(a, b, c)` — fused multiply-add with single final rounding.
- Criterion benchmarks in `zenith-float-num/benches/` (arithmetic, transcendentals, composite); `./scripts/bench.sh`.
- FFT-scale multiply benchmark tier (`arithmetic/mul_fft`, 346k and 524k bits).
- TSV benchmark baselines (`doc/bench-baselines.tsv`) and `scripts/bench-compare.sh` for regression checks without JSON.
- Cross-library compare harness (`zenith-float-compare/`, `scripts/compare-bench.sh`) for release vs astro/dashu.
- Newton reciprocal at three words and up.
- `hypot`, `atan2`, `log1p`, and `expm1` on `ExactNum` and in `expr!`.
- `exp2`, `exp10`, and public `rem_pi` on `ExactNum` and in `expr!`.

## 0.1.0

First public release.

- Software-limb `ExactNum` (no hardware floating-point in calculations).
- Elementary functions, constants cache, radix parse/format, `expr!`.
- Optional `std`, `random`, and `serde` features.
- Dual license MIT OR Apache-2.0.
