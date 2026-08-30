# zenith-float — Complete Build Plan

Every item on this list must be implemented, golded, and passing CI before zenith-float is complete. No item is optional. No item is deferred. Build in the order listed — each section's dependencies are satisfied by the sections above it.

## Rules for every item

- Software limb arithmetic only — no hardware floating-point arithmetic, `libm`, or compiler float codegen
- Software `Ieee32` / `Ieee64` are permitted as interchange types using integer bit manipulation only
- Working precision `p_wrk = p + WORD_BIT_SIZE`; `MAX_PREC_RETRY` bounds retries → `PrecisionRetryExhausted`
- Per-op precision at caller-chosen bits — never global `SOFT_PREC`
- Domain guards → `NaN` (`InvalidArgument`); overflow → `±Inf`; never invent values outside the stated domain
- No `unwrap` without a proof comment
- No dead code, no commented-out code, no magic numbers — named constants only
- Every new method has a doc comment stating its mathematical contract
- Every item needs a gold on the method callers will use — not `is_some()`, not a string check, a locked expected value
- `expr!` or `cexpr!` leaf added where arity fits
- MPFR oracle gold under `mpfr-tests` where MPFR has the function
- Derivative identity golded alongside the function
- `scripts/ci_full.sh` green before marking done

---

## Section 1 — Complex specials completion (ExactComplex)

Elementary `+−×÷`, `exp`/`ln`, trig/hyperbolic, principal `sqrt`/`pow`/inverses already exist. The items below complete the complex special function set. Do not wrap the real series on `|z|` and call it complex — use proper complex algorithms throughout.

### 1.1 Complex `Ei` / `Si` / `Ci` / `li` / Fresnel

**What:** `ExactComplex::ei`, `si`, `ci`, `li`, `fresnel_s`, `fresnel_c`

**Prompt:**
Implement `ExactComplex::ei(p, rm, cc)`, `si`, `ci`, `li`, `fresnel_s`, `fresnel_c` on `ExactComplex` in zenith-float. Software limbs only. No hardware float.

Algorithm:
- `Ei(z)`: power series for `|z| < EI_SERIES_THRESHOLD` (named constant); asymptotic expansion for large `|z|`, stopped at the smallest term. Cut on `(-∞, 0]` — `Ei(-1+0i)` and `Ei(-1-0i)` must differ by `2πi`.
- `Si(z)` and `Ci(z)`: use `Si(z) = (Ei(iz) - Ei(-iz))/(2i) - π/2` and `Ci(z) = -(Ei(iz) + Ei(-iz))/2`. Inherits `Ei` branch cut. `Ci(0)` is a pole → `NaN`.
- `li(z)`: `Ei(ln z)`. Cut on `(-∞, 1]`. `li(1)` is a pole → `NaN`.
- `fresnel_s(z)`, `fresnel_c(z)`: use `erf` which is already implemented. `fresnel_s(z) = (1+i)/4 · (erf((1+i)√π z/2) - erf((1-i)√π z/2))` and analogous for `fresnel_c`. Both entire.

`cexpr!` leaves: `ei(z)`, `si(z)`, `ci(z)`, `li(z)`, `fresnel_s(z)`, `fresnel_c(z)`.

Golds:
- `ei(1+0i)` matches real `Ei(1)` to working precision
- `si(0+0i) = 0`; `ci` pole at `z=0` → `NaN`
- `Si(-z) = -Si(z)` symmetry gold
- `fresnel_s(0) = 0`, `fresnel_c(0) = 0`
- `D(Ei)(z) = e^z/z` derivative gold
- `D(Si)(z) = sin(z)/z` derivative gold
- `li(e+0i)` matches real `li(e)` to working precision
- Cut gold: `ei(-1+0i)` vs `ei(-1-0i)` differ by `2πi`

Named constants: `EI_SERIES_THRESHOLD`. No magic numbers. `scripts/ci_full.sh` green.

---

### 1.2 Complex Bessel `J_ν` / `Y_ν` / `I_ν` / `K_ν`

**What:** `ExactComplex::bessel_j_nu`, `bessel_y`, `bessel_i`, `bessel_k`

**Prompt:**
Implement `ExactComplex::bessel_j_nu(nu, p, rm, cc)`, `bessel_y(nu, p, rm, cc)`, `bessel_i(nu, p, rm, cc)`, `bessel_k(nu, p, rm, cc)` on `ExactComplex`. Software limbs only.

Algorithm:
- Build on Hankel functions: `H_ν^(1)(z) = J_ν(z) + iY_ν(z)`, `H_ν^(2)(z) = J_ν(z) - iY_ν(z)`
- For `|z| < BESSEL_SERIES_THRESHOLD` (named): power series for `J_ν`; `Y_ν` via the Wronskian relation
- For large `|z|`: uniform asymptotic expansion, stopped at the smallest term
- `I_ν(z) = i^{-ν} J_ν(iz)`; `K_ν(z) = (π/2) i^{ν+1} H_ν^(1)(iz)`
- Branch cuts: `J_ν`, `I_ν` entire for integer `ν`; cut on `(-∞, 0]` for non-integer `ν`. `Y_ν`, `K_ν`: cut on `(-∞, 0]`.

`cexpr!` leaves: `bessel_j_nu(z, ν)`, `bessel_y(z, ν)`, `bessel_i(z, ν)`, `bessel_k(z, ν)`.

Golds:
- `J_0(1+0i)`, `J_1(1+0i)` match real values to working precision
- Wronskian `J_ν(z)Y_{ν+1}(z) - J_{ν+1}(z)Y_ν(z) = -2/(πz)` at a non-trivial complex point
- `I_0(1+0i)` matches real `I_0(1)`; `K_0(1+0i)` matches real `K_0(1)`
- `I_ν(z) = i^{-ν} J_ν(iz)` identity gold
- `D(J_0)(z) = -J_1(z)` derivative gold
- `z = 0` for non-integer `ν` → `NaN`
- Cut: `J_{1/2}(-1+0i)` vs `J_{1/2}(-1-0i)` differ correctly

Named constants: `BESSEL_SERIES_THRESHOLD`. No magic numbers.

---

### 1.3 Complex elliptic integrals

**What:** `ExactComplex::elliptic_k`, `elliptic_e_complete`, `elliptic_f`, `elliptic_e`, `elliptic_pi_complete`, `elliptic_pi`

**Prompt:**
Implement complex elliptic integrals on `ExactComplex` by extending the existing real Carlson symmetric form kernels `R_F`, `R_C`, `R_D`, `R_J` to complex arguments. The duplication algorithm works unchanged in ℂ when arguments avoid the branch cuts. `CARLSON_DUPE_MAX = 128` unchanged.

Convention: parameter `m = k²`, incomplete argument `x = sin φ` — same as the real implementation.

Branch cuts: `K(m)` has a cut on `[1, +∞)` in `m`. Document incomplete form cuts explicitly in doc comments and pin by gold.

`cexpr!` leaves: same names as real leaves — `elliptic_k(m)`, `elliptic_e(m)`, `elliptic_f(x, m)`, `elliptic_e_inc(x, m)`, `elliptic_pi(n, m)`, `elliptic_pi_inc(n, x, m)`.

Golds:
- `K(0+0i) = π/2`, `E(0+0i) = π/2` match real values
- `K(0.5+0i)` matches real `K(0.5)` to working precision
- Legendre relation `E(m)K'(m) + E'(m)K(m) - K(m)K'(m) = π/2` at a complex `m`
- `F(x, 0+0i) = arcsin(x)` identity gold
- Cut: `K(2+0i)` vs `K(2-0i)` differ correctly

---

### 1.4 Complex `₂F₁`

**What:** `ExactComplex::hypergeom_2f1(a, b, c, z, p, rm, cc)` where all parameters are `ExactComplex`

**Prompt:**
Implement `ExactComplex::hypergeom_2f1` where `a`, `b`, `c`, `z` are all `ExactComplex`. Software limbs only.

Algorithm:
- Series `Σ (a)_n(b)_n/(c)_n · z^n/n!` for `|z| < 1`; cap `HYPERGEOM_SERIES_MAX_TERMS` (named)
- Euler transformation `(1-z)^{c-a-b} ₂F₁(c-a,c-b;c;z)` for `|1-z| < 1`
- Pfaff transformation `(1-z)^{-a} ₂F₁(a,c-b;c;z/(z-1))` for `Re(z) < 1/2`
- Kummer evaluation at `z = 1` when `Re(c-a-b) > 0`
- Cut on `[1, +∞)` in `z`; `c` a non-positive integer → `NaN`
- `NaN` when no transformation brings `z` into the convergence region within `HYPERGEOM_TRANSFORM_MAX` (named) attempts

`cexpr!` leaf: `hypergeom_2f1(a, b, c, z)`.

Golds:
- `₂F₁(1,1;2;0.5+0i)` matches real value `2ln2`
- `₂F₁(0.5,0.5;1;0.5+0i)` matches real value `2K(0.5)/π`
- `₂F₁(a,b;c;0+0i) = 1` for any complex `a`, `b`, `c`
- Euler transformation identity at a non-trivial complex `z`
- `c = 0+0i` → `NaN`
- Cut: `₂F₁(1,1;2;2+0i)` vs `₂F₁(1,1;2;2-0i)` differ correctly

---

## Section 2 — Array ufuncs

Every scalar `ExactNum` special available elementwise on `ExactNumArray`, `Ieee32Array`, `Ieee64Array`.

### 2.1 ExactNumArray ufuncs — elementary

**Prompt:**
Add elementwise methods to `ExactNumArray` for every elementary transcendental already on `ExactNum`: `sqrt`, `cbrt`, `ln`, `log2`, `log10`, `log1p`, `exp`, `exp2`, `exp10`, `expm1`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`, `abs`, `signum`, `ceil`, `floor`, `int`, `fract`.

Each method signature: `fn name(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> ExactNumArray` or without `cc` where the scalar method does not need it. Shape is preserved. A `NaN` element propagates to `NaN` in the output at that position.

Golds: elementwise `sin` on a `(2×3)` array matches applying `sin` to each element individually. Shape mismatch between two arrays in a binary ufunc returns `None`.

---

### 2.2 ExactNumArray ufuncs — specials

**Prompt:**
Add elementwise methods to `ExactNumArray` for every special function on `ExactNum`: `erf`, `erfc`, `gamma`, `ln_gamma`, `digamma`, `gammainc`, `gammainc_upper`, `ei`, `si`, `ci`, `li`, `fresnel_s`, `fresnel_c`, `bessel_j_nu`, `bessel_y`, `bessel_i`, `bessel_k`, `elliptic_k`, `elliptic_e_complete`, `elliptic_f`, `elliptic_e`, `elliptic_pi_complete`, `elliptic_pi`, `legendre_p`, `assoc_legendre_p`, `hypergeom_2f1`, `betainc`.

Same rules as §2.1. Functions that take an extra parameter (e.g. `bessel_j_nu(nu, p, rm, cc)`) take that parameter as a scalar applied uniformly across the array.

Golds: elementwise `bessel_j_nu(0.5)` on a 1D array matches scalar results. Domain violations (`NaN` inputs) propagate correctly.

---

### 2.3 Ieee32Array / Ieee64Array ufuncs

**Prompt:**
Add elementwise special function methods to `Ieee32Array` and `Ieee64Array` via conversion through `ExactNum` at a working precision sufficient to round correctly to binary32/binary64. The conversion path: `Ieee64::to_exact(p_wrk)` → apply `ExactNum` special → round back to `Ieee64`. No hardware float at any point.

Cover: `sqrt`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `exp`, `ln`, `erf`, `erfc`, `gamma`, `bessel_j_nu`, `bessel_y`, `bessel_i`, `bessel_k`.

Golds: `Ieee64Array::sin` matches `ExactNum::sin` rounded to binary64 for a representative set of values.

---

### 2.4 Integer SIMD for IEEE arrays

**Prompt:**
Add SIMD-accelerated paths for `Ieee32Array` and `Ieee64Array` elementwise `add`, `sub`, `mul`, `div`, `sqrt`, and `fma` using integer SIMD lanes (`u32x8` / `u64x4` or equivalent via `std::simd` or `packed_simd2`). The arithmetic is software IEEE — integer lanes carrying bit patterns, arithmetic implemented in software, no hardware FPU instructions. The scalar and SIMD paths must produce bit-identical results.

Named constant: `IEEE_SIMD_LANE_WIDTH`. Use `cfg` feature gates so the scalar fallback compiles everywhere.

Golds: SIMD `add`/`mul`/`div`/`sqrt` on a 1000-element `Ieee64Array` produces bit-identical output to the scalar path. Results match `ExactNum` rounded to binary64.

---

## Section 3 — Certified interval arithmetic

### 3.1 Ball arithmetic for transcendentals

**Prompt:**
Extend `Ball { mid, rad }` to support certified enclosures for: `sin`, `cos`, `exp`, `ln`, `sqrt`, `erf`, `bessel_j(0, ...)`, `bessel_j(1, ...)`.

For each function `f`, implement `Ball::f(p, rm, cc) -> Ball` such that the returned ball provably contains `f(mid)` — i.e. `|f(mid) - result.mid| ≤ result.rad`. Use interval extensions of the existing `ExactNum` kernels plus a certified error bound derived from the Taylor remainder or the function's known Lipschitz constant on the ball.

Named constant: `BALL_TRANSCENDENTAL_ERROR_TERMS` — number of Taylor terms used in the error bound.

Golds:
- `Ball::sin` on `[π/6 ± 2^{-p}]` contains `1/2`
- `Ball::exp` on `[1 ± 2^{-p}]` contains `e`
- `Ball::erf` on `[1 ± 2^{-p}]` contains real `erf(1)`
- Composing two `Ball` ops produces a ball that still contains the true value — verify by checking `mid ± rad` encloses the `ExactNum` result

---

### 3.2 Certified interval for complex functions

**Prompt:**
Extend `Ball` to `ComplexBall { mid: ExactComplex, rad: ExactNum }` (a disk in ℂ). Implement `add`, `mul`, `exp`, `ln`, `sin`, `cos` with certified enclosures on `ComplexBall`. The radius grows by the Lipschitz constant of the operation on the disk.

Golds:
- `ComplexBall::exp` on the unit disk contains `e^{mid}`
- Composition of `ComplexBall::sin` and `ComplexBall::cos` on a small disk satisfies the Pythagorean identity enclosure

---

## Section 4 — Arbitrary-precision linear algebra

### 4.1 LU decomposition

**Prompt:**
Implement LU decomposition with partial pivoting on `ExactNumArray` (2D, square or rectangular). Return `(L, U, P)` where `P` is a permutation vector. Use exact `ExactNum` arithmetic with explicit `(p, rm)` at each step. A singular matrix returns `None` for the `U` factor — do not invent a result.

Golds:
- `[[2,1],[4,3]]` → `L`, `U`, `P` such that `P·A = L·U` exactly
- Singular matrix `[[1,2],[2,4]]` → `None`
- `L·U` reconstructs `P·A` to working precision

---

### 4.2 QR decomposition (Gram-Schmidt)

**Prompt:**
Implement QR decomposition via modified Gram-Schmidt on `ExactNumArray`. Return `(Q, R)` where `Q` is orthogonal and `R` is upper triangular, all at explicit precision `(p, rm)`. Rank-deficient input: zero column in `R` at the deficient position, not a panic.

Golds:
- `Q·R` reconstructs `A` to working precision
- `Q^T · Q = I` to working precision
- Rank-deficient matrix handled correctly

---

### 4.3 SVD (Golub-Reinsch)

**Prompt:**
Implement SVD on `ExactNumArray` via the Golub-Reinsch bidiagonalization algorithm at explicit precision `(p, rm)`. Return `(U, Σ, V^T)`. Singular values in descending order. Software limbs only — no LAPACK, no hardware float.

Cap: `SVD_ITER_MAX` named constant for the QR iteration convergence check. Return `None` if convergence is not reached within the cap.

Golds:
- `U · Σ · V^T` reconstructs `A` to working precision
- `U^T · U = I`, `V^T · V = I` to working precision
- Known singular values of `[[3,0],[0,2]]` are `3` and `2`

---

### 4.4 Eigenvalue decomposition

**Prompt:**
Implement eigenvalue decomposition for real symmetric matrices via the symmetric QR algorithm (tridiagonalization + QR iteration) at explicit precision `(p, rm)`. Return `(eigenvalues, eigenvectors)` as `(ExactNumArray, ExactNumArray)`. Non-symmetric input returns `None` — do not attempt a general eigendecomposition here.

Cap: `EIGEN_ITER_MAX` named constant.

Golds:
- `[[2,1],[1,2]]` eigenvalues are `1` and `3`; eigenvectors are orthogonal
- `A · v = λ · v` for each eigenpair to working precision
- `V · Λ · V^T = A` reconstruction to working precision

---

### 4.5 Multiprecision FFT

**Prompt:**
Implement a radix-2 Cooley-Tukey FFT on `ExactNumArray` (complex, length must be a power of 2) at explicit precision `(p, rm, cc)`. Twiddle factors computed from `ExactNum::sin_cos` at working precision. Software limbs only.

Named constant: `FFT_MAX_POINTS` — maximum supported length.

Golds:
- FFT of `[1, 0, 0, 0]` is `[1, 1, 1, 1]`
- FFT of a pure cosine at frequency `k` has energy only at bin `k` and `N-k`
- IFFT(FFT(x)) = x to working precision
- Parseval's theorem: `sum |x_n|² = (1/N) sum |X_k|²` to working precision

---

## Section 5 — Precision documentation per function

### 5.1 Error bound documentation

**Prompt:**
For every special function in §12 of `LIBRARY.md` and `ZENITH_FLOAT_CAPABILITIES.md`, add a doc comment section "Precision" that states:
- The error bound in ULP for the primary series/algorithm region
- The argument range where each algorithm is used (series threshold, asymptotic threshold)
- The named constants that control those thresholds
- Whether an MPFR oracle gold exists and what domain it covers

This is documentation only — no code changes. The doc comments go on the `ExactNum` and `ExactComplex` methods. Format: follow the existing `erf` doc comment as a template.

---

## Section 6 — no_std certified builds

### 6.1 no_std compliance audit and fix

**Prompt:**
Audit every item in zenith-float for `no_std` compatibility. For each item that currently requires `std` beyond formatting traits: either provide an allocator-only alternative or document explicitly in the doc comment that the `std` feature is required and why.

Specifically:
- All special functions must compile with `default-features = false` and a global allocator
- `Ball` and `ComplexBall` must compile without `std`
- Array types must compile without `std`
- `SharedConsts` correctly gates on `std` feature
- CI adds a `no_std` build target to `scripts/ci_full.sh`

Golds: `cargo build --no-default-features --target thumbv7em-none-eabihf` succeeds for the kernel crate.

---

## Section 7 — Reproducibility story

### 7.1 Reproducibility documentation

**Prompt:**
Add a `doc/REPRODUCIBILITY.md` to zenith-float with the following content:

- Statement that all results are determined by the source code, input values, precision `p`, and rounding mode `rm` — no platform-dependent behavior, no hardware float
- Instructions for running the full gold suite: `cargo test --features mpfr-tests`
- Statement that `ZENITH_TEST_SEED` allows replay of any random test failure
- Statement that the gold values in `golds/` are locked expected outputs, not tolerances
- Citation format: "Computed with zenith-float {version}, reproducible by running `cargo test` at commit {hash}"

This document is for researchers who need to cite numerical results in papers.

---

## Section 8 — Version and publish checklist

### 8.1 Pre-publish checklist

**Prompt:**
Create `scripts/prepublish.sh` that runs the following and fails loudly on any error:

1. `cargo test --all-features` — all tests pass
2. `cargo test --features mpfr-tests` — MPFR oracle golds pass
3. `cargo build --no-default-features` — no_std compiles
4. `scripts/compare-bench.sh --quick` — benchmark baseline refreshed
5. `grep -r "f32\|f64\|f128" zenith-float-num/src/ --include="*.rs"` — kernel purity check (hardware float type names absent; software `Ieee32`/`Ieee64` are permitted as they use integer arithmetic)
6. Verify `ZENITH_FLOAT_CAPABILITIES.md` version and date match `Cargo.toml`
7. Verify `doc/LIBRARY.md` `expr!` leaf list matches the macro expansion

Script exits 0 only when all checks pass. No publish without a green prepublish run.
