# Changelog

## Unreleased

- `cexpr!` — per-part `errs[]` for real vs imaginary cancellation; principal branch cuts documented; leaves match the complex-capable `expr!` set (`cbrt`/`root`, logs/exps, `hypot`/`fma`, `abs`/`arg`/`conj`, `ldexp`/`scalb`/`logb`). No `atan2` or `rem_pi` (complex trig uses `x+iy` identities).
- `ExactComplex::{log2, log10, log, log1p, exp2, exp10, expm1, ldexp, scalb, logb, cbrt, nth_root, hypot, fma, mul_add}`.
- Getting-started guide (`doc/GETTING_STARTED.md`) and help tutorial (`doc/HELP.md`).
- Maintainer build checklist is not packaged on crates.io.
- `Context::with_rounding_mode`.
- `ExactNum::{two_sum, two_product, fused_sum, fused_dot, polyval}`.
- `ExactComplex`: `tan`, `sinh`, `cosh`, `tanh`, `sqrt`, `pow`, `asin`, `acos`, `atan`, `asinh`, `acosh`, `atanh`.
- `ExactComplex` rectangular complex type (`+ − × ÷`, `exp`, `ln`, `sin`, `cos`, `abs`, `arg`).
- `ExactNum::sin_cos`; `mul_add` alias of `fma`.
- Special functions: `erf` / `erfc`, `gamma` / `ln_gamma`, integer-order `bessel_j`; `expr!` leaves.
- MPFR compare coverage for `fma`, paired `sin_cos` / `sinh_cosh`.
- MPFR compare coverage for `exp2`, `exp10`, and `rem_pi`; release gate on Linux x86_64 in `scripts/ci.sh`.
- MPFR fuzz: `add` / `mul` / `sqrt` bit-exact under all five IEEE rounding modes (`tests/mpfr/fuzz_round_modes.rs`).
- Progressive `ConstCache` (`cache_info`, √2, φ) and `CachedFBig` extra-precision wrapper.
- Public `Ball` enclosures and `ziv_round` correct-rounding retry helper.
- Stack-inlined mantissas (`INLINE_WORDS` = 2 limbs) before heap promotion.
- `doc/EXPR.md` — `expr!` per-op rounding contract (working precision, cancellation slots, final `set_precision`).
- `fma` skips the full-width product when `|a*b|` and `c` differ by more than `p + 2` words.
- MPFR fuzz covers `nth_root` for n=4,5 at 1-ULP (`mpfr_rootn_ui`).
- Optional nightly bench gate: `CI_BENCH=1 ./scripts/ci.sh`.
- `Consts::euler_gamma` (Euler–Mascheroni γ) on the progressive constant cache.
- `ExactNum::{frexp, ldexp, scalb, logb, ilogb}` — IEEE-style exponent split without hardware floats.
- `LowerExp` / `UpperExp` (`{:e}` / `{:E}`) and `LowerHex` formatting.
- `SharedConsts` — mutex-wrapped constant cache for parallel batch evaluation (`std`).
- MPFR oracles for `erf`/`erfc`, `gamma`/`ln_gamma`, `bessel_j`, extra constants, and `ExactComplex` add.
- `expr!` constants `sqrt2`, `phi`, `euler_gamma`.
- `expr!` leaves `ldexp` / `scalb` / `logb`.
- `serde`: JSON round-trip for `ExactNum` (decimal string / integer, Inf/NaN) and `ExactComplex`; feature implies `std`.
- Structured radix 2–36 parse/format round-trip stress test (`tests/radix_roundtrip.rs`).
- Native `+`, `-`, `*`, `/` operators on `ExactNum` (default precision, round-to-even).
- `exact!` / `fbig!` compile-time decimal float literals.
- Parse/format for radices 2–36 (`Radix::try_new`); `_e` exponent separator for bases > 10.
- `RadixFloat` wrapper and `ExactNum::with_radix` for radix-tagged I/O.
- `ExactNum::nth_root` / `expr!` `root(x, n)` — general n-th root (`sqrt`/`cbrt` delegation, composite factors, Newton for primes).
- `ExactNum::sinh_cosh` — paired hyperbolic evaluation with a single `exp(|x|)` path.
- `ExactNum::copysign`, `ExactNum::next_after` — software IEEE sign and successor.
- MPFR bit-oracle for `copysign` at high precision (`compare_copysign_test`); sign is applied before rounding so directed modes match MPFR.
- Seeded test RNG (`reseed_random`, `ZENITH_TEST_SEED`; panic hook prints the seed).
- `scripts/ci.sh` enforces debug / MPFR wall-time budgets (10 min / 30 min).
- Allocation failure at huge precision returns `NaN` (`Error::MemoryAllocation`), with a large-precision add/mul smoke test.
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
