# Changelog

## Unreleased

- Criterion benchmarks in `zenith-float-num/benches/` (arithmetic, transcendentals, composite); `./scripts/bench.sh`.
- FFT-scale multiply benchmark tier (`arithmetic/mul_fft`, 346k and 524k bits).
- Newton reciprocal at three words and up.
- `hypot`, `atan2`, `log1p`, and `expm1` on `ExactNum` and in `expr!`.
- `exp2`, `exp10`, and public `rem_pi` on `ExactNum` and in `expr!`.

## 0.1.0

First public release.

- Software-limb `ExactNum` (no hardware floating-point in calculations).
- Elementary functions, constants cache, radix parse/format, `expr!`.
- Optional `std`, `random`, and `serde` features.
- Dual license MIT OR Apache-2.0.
