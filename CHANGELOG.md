# Changelog

## Unreleased

- Quadrature: `gauss_legendre`, `tanh_sinh`, `gauss_laguerre`, `gauss_hermite`. Golds: \(x^2\) on \([-1,1]\) is \(2/3\); 20-point \(x^{38}\) is \(2/39\); \(1/\sqrt{1-x^2}=\pi\); Laguerre \(x^2=\Gamma(3)=2\); Hermite \(1=\sqrt{\pi}\).
- Orthogonal polynomials: `hermite_he`/`hermite_h`, `laguerre`/`gen_laguerre`, `chebyshev_t`/`chebyshev_u`, `gegenbauer`. Golds: `He_4(0)=3`, `L_3(0)=1`, `T_5(cos(π/5))=-1`, `C_2^{(1)}=4x²−1=U_2`, `2C_2^{(1/2)}=3x²−1`, `T_6=2xT_5−T_4`. `ORTHOPOLY_N_MAX=256`.
- Chebyshev interpolation: `chebyshev_coeffs` / `chebyshev_eval` / `clenshaw` / `chebyshev_error_bound`. `exp` on `[-1,1]` with 20 terms error `< 10^{-15}`; Clenshaw of `[1,2,3]` at `1/2` is `1/2`; `CHEBYSHEV_MAX_DEGREE=256`.
- `ExactNumPoly` — dense univariate (low-to-high coeffs); `div_rem(x²−1, x−1)=(x+1, 0)`; `gcd=x−1`; `compose(x², x+1)=x²+2x+1`; `∂(x³)=3x²`; `roots_real(x²−2)=±√2` (companion closed form through `POLY_COMPANION_CLOSED_DEG=2`).
- RNG: `random_uniform`, `random_gaussian` (Box–Muller), `random_exponential`, `ExactNumArray::random_fill` / `RandomDist`. Existing `random_normal(p, exp_from, exp_to)` is unchanged.
- Distribution kernels: `normal_pdf`/`cdf`, `gamma_pdf`, `beta_pdf`, `poisson_pmf`, `binomial_pmf`, `chi_squared_cdf`, `student_t_pdf`. Golds: \(1/\sqrt{2\pi}\), \(1/2\), \(e^{-1}\), \(\chi^2_2(2\ln 20)=19/20\).
- `ziv_round_vec` — same Ziv loop on a slice; `hypot(3,4)=5`; `atan2(1,1)=π/4` at 256 bits.
- `parse_exact` / `format_exact` — `0.1` is exact `1/10`; dyadic `0.5` / `0.125`; `parse("1.5e3")=1500`.
- `ExactInt` — limb integer; `20!`, `gcd(48,18)=6`, `2^100`, `div_rem(17,5)=(3,2)`.
- `ExactRational` — exact `num/den` with integer GCD reduction; `1/3+1/6=1/2`; `2/4=1/2`; sign onto the numerator; 256-bit `1/3`.
- `doc/REPRODUCIBILITY.md` — what determines a result, how to replay tests, citation line.
- `# Precision` rustdoc on `ExactNum` and `ExactComplex` specials (algorithm, thresholds, Ziv / ULP, MPFR oracle).
- `ExactNumArray::fft` / `ifft` — radix-2 Cooley–Tukey; impulse / cosine / Parseval golds; `FFT_MAX_POINTS=4096`.
- `ExactNumArray::eigen_decomp` — symmetric QR; \(Av=\lambda v\), \(V\Lambda V^T=A\); \(\begin{pmatrix}2&1\\1&2\end{pmatrix}\to(3,1)\); `EIGEN_ITER_MAX=64`.
- `ExactNumArray::svd_decomp` — Golub–Reinsch; \(U\Sigma V^T=A\), \(U^\top U=V^\top V=I\); \(\operatorname{diag}(3,2)\to\sigma=(3,2)\); `SVD_ITER_MAX=64`.
- Living docs are the build plan, capabilities, and build checklist only. `ZENITH_FLOAT_TODO.md` is retired.
- `ExactNumArray::qr_decomp` — modified Gram–Schmidt; \(QR=A\), \(Q^\top Q=I\); rank-deficient zero diagonal.
- `ExactNumArray::lu_decomp` — partial pivoting; \(PA=LU\) gold; singular → `None`.
- `ComplexBall` disk enclosures: `add`/`mul`/`exp`/`ln`/`sin`/`cos`. Unit-disk `exp` and \(\sin^2+\cos^2=1\) golds.
- `Ball::{cos,ln,sqrt,erf,bessel_j0,bessel_j1}` certified enclosures (`BALL_TRANSCENDENTAL_ERROR_TERMS`). `ExactNumArray` `(2×3)` `sin` and `bessel_j_nu` golds; `signum` ufunc.
- `ExactComplex::hypergeom_2f1` — series / Euler / Pfaff / Kummer in \(\mathbb{C}\). `cexpr!` leaf `hypergeom_2f1`. Cut on \([1,+\infty)\); non-positive integer \(c\) → NaN.
- `ExactComplex` elliptic \(K,E,\Pi\) (complete and incomplete) via Carlson \(R_F,R_C,R_D,R_J\) in \(\mathbb{C}\). `cexpr!` leaves `elliptic_k` / `elliptic_e` / `elliptic_f` / `elliptic_e_inc` / `elliptic_pi` / `elliptic_pi_inc`. \(K(1)=+\infty\); cut of \(K\) on \([1,+\infty)\).
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
