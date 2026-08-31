# zenith-float — Complete Build Plan

Every item on this list must be implemented, golded, and passing CI before zenith-float is complete. No item is optional. No item is deferred. Build in the order listed — each section's dependencies are satisfied by the sections above it.

**Living docs (same trio as Accumath — there is no TODO file):**

| Document | Role |
| --- | --- |
| This file | Walk list. Status table at the top. Flip a row when a gold lands. |
| [`ZENITH_FLOAT_CAPABILITIES.md`](ZENITH_FLOAT_CAPABILITIES.md) | What the crate can do. Must match the object. |
| [`BUILD_CHECKLIST.md`](BUILD_CHECKLIST.md) | Maintainer: what is in CI, crates, and the public API today. |

After every slice: update all three. Do not keep a fourth inventory.

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

## Status vs tree (2026-08-30)

Compared to `zenith-float-num` / macros / docs. ✅ = method + object gold. 🟡 = present but short of the prompt. ⬜ = not on the object. Section 11 duplicates Section 4.

| Item | Status | On the object |
| --- | --- | --- |
| 1.1 Complex Ei/Si/Ci/li/Fresnel | ✅ | `complex_ei.rs`; `cexpr!` leaves; cut / identity golds. Plan derivative golds not separately locked |
| 1.2 Complex Bessel | ✅ | `complex_bessel.rs`; Wronskian / cut golds |
| 1.3 Complex elliptic | ✅ | `complex_elliptic.rs`; Legendre + cut golds |
| 1.4 Complex `_2F1` | ✅ | `complex_hypergeom.rs`; `2ln2`, `2K/π`, Euler, `c=0`, cut |
| 2.1 ExactNumArray elementary ufuncs | ✅ | `(2×3)` `sin`; shape mismatch → `None`; `signum` |
| 2.2 ExactNumArray specials | ✅ | `bessel_j_nu(1/2)` matches scalar; NaN propagates |
| 2.3 Ieee32/64 array ufuncs | ✅ | via `ExactNum`; `bin64_array_exp_sin` |
| 2.4 Integer SIMD | 🟡 | add/mul lanes only; no named `IEEE_SIMD_LANE_WIDTH`; no SIMD div/sqrt/fma |
| 3.1 Ball transcendentals | ✅ | `sin`/`cos`/`exp`/`ln`/`sqrt`/`erf`/`J0`/`J1` |
| 3.2 ComplexBall | ✅ | disk add/mul/exp/ln/sin/cos |
| 4.1 / 11.1 LU | ✅ | `lu_decomp`; \(PA=LU\); singular → `None` |
| 4.2 / 11.2 QR | ✅ | `qr_decomp` modified Gram–Schmidt; \(QR=A\); \(Q^\top Q=I\); rank-deficient zero diagonal |
| 4.3 / 11.3 SVD | ✅ | `svd_decomp`; \(U\Sigma V^T=A\); \(U^\top U=V^\top V=I\); \(\operatorname{diag}(3,2)\); wide \(\sigma=(2,1)\) |
| 4.4 / 11.4 eigen | ✅ | `eigen_decomp`; \(Av=\lambda v\); \(V\Lambda V^T=A\); \(\begin{pmatrix}2&1\\1&2\end{pmatrix}\to(3,1)\); non-symmetric `None` |
| 4.5 / 11.5 FFT | ✅ | `fft` / `ifft`; impulse `[1,0,0,0]→[1,1,1,1]`; cosine bins; IFFT; Parseval; `FFT_MAX_POINTS=4096` |
| 5.1 / 18.1 Precision doc comments | ✅ | `# Precision` on `ExactNum` / `ExactComplex` specials (algorithm, thresholds, ULP/Ziv, MPFR) |
| 6.1 no_std / thumb | 🟡 | allocator `no_std` compiles; no `thumbv7em-none-eabihf` CI gold |
| 7.1 / 20.1 Reproducibility.md | ✅ | `doc/REPRODUCIBILITY.md`; unit tests lock values (no `golds/` tree) |
| 8 | — | skipped by plan |
| 9.1 ExactRational | ⬜ | |
| 9.2 ExactInt | ⬜ | |
| 9.3 parse_exact / format_exact | ⬜ | scientific parse exists; exact-min-prec APIs do not |
| 10.1 Per-op precision | ✅ | specials take `(p,rm,cc)`; no `SOFT_PREC` in the kernel |
| 10.2 expr! / cexpr! `p_wrk` | 🟡 | leaves use `p_wrk`; not every listed composite gold is locked |
| 10.3 `ziv_round_vec` | ⬜ | |
| 12–17, 18.2–18.3, 19–20 | ⬜ | distributions, poly, quadrature, roots, ODE, DSP, crypto, serde-all, binary I/O, HDF5, HELP rewrite, MPFR extend, proptest, prepublish, hex CI |

Walk this table top to bottom. Do not start a later ⬜ while an earlier ⬜ remains. Next implementation slice: **§9.1 ExactRational**.

When a row flips, add `**Status:** done YYYY-MM-DD` under that section heading and update CAPABILITIES + BUILD_CHECKLIST in the same session.

---

## Section 1 — Complex specials completion (ExactComplex)

Elementary `+−×÷`, `exp`/`ln`, trig/hyperbolic, principal `sqrt`/`pow`/inverses already exist. The items below complete the complex special function set. Do not wrap the real series on `|z|` and call it complex — use proper complex algorithms throughout.

### 1.1 Complex `Ei` / `Si` / `Ci` / `li` / Fresnel

**Status:** done 2026-08-30 — `complex_ei.rs`; `cexpr!` leaves; cut / identity golds.

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

**Status:** done 2026-08-30 — `complex_bessel.rs`; Wronskian / cut golds.

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

**Status:** done 2026-08-30 — `complex_elliptic.rs`; Legendre + cut golds.

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

**Status:** done 2026-08-30 — `complex_hypergeom.rs`; `2ln2`, `2K/π`, Euler, `c=0`, cut.

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

**Status:** done 2026-08-30 — `(2×3)` `sin`; shape mismatch → `None`; `signum`.

**Prompt:**
Add elementwise methods to `ExactNumArray` for every elementary transcendental already on `ExactNum`: `sqrt`, `cbrt`, `ln`, `log2`, `log10`, `log1p`, `exp`, `exp2`, `exp10`, `expm1`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`, `abs`, `signum`, `ceil`, `floor`, `int`, `fract`.

Each method signature: `fn name(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> ExactNumArray` or without `cc` where the scalar method does not need it. Shape is preserved. A `NaN` element propagates to `NaN` in the output at that position.

Golds: elementwise `sin` on a `(2×3)` array matches applying `sin` to each element individually. Shape mismatch between two arrays in a binary ufunc returns `None`.

---

### 2.2 ExactNumArray ufuncs — specials

**Status:** done 2026-08-30 — `bessel_j_nu(1/2)` matches scalar; NaN propagates.

**Prompt:**
Add elementwise methods to `ExactNumArray` for every special function on `ExactNum`: `erf`, `erfc`, `gamma`, `ln_gamma`, `digamma`, `gammainc`, `gammainc_upper`, `ei`, `si`, `ci`, `li`, `fresnel_s`, `fresnel_c`, `bessel_j_nu`, `bessel_y`, `bessel_i`, `bessel_k`, `elliptic_k`, `elliptic_e_complete`, `elliptic_f`, `elliptic_e`, `elliptic_pi_complete`, `elliptic_pi`, `legendre_p`, `assoc_legendre_p`, `hypergeom_2f1`, `betainc`.

Same rules as §2.1. Functions that take an extra parameter (e.g. `bessel_j_nu(nu, p, rm, cc)`) take that parameter as a scalar applied uniformly across the array.

Golds: elementwise `bessel_j_nu(0.5)` on a 1D array matches scalar results. Domain violations (`NaN` inputs) propagate correctly.

---

### 2.3 Ieee32Array / Ieee64Array ufuncs

**Status:** done 2026-08-30 — via `ExactNum`; `bin64_array_exp_sin`.

**Prompt:**
Add elementwise special function methods to `Ieee32Array` and `Ieee64Array` via conversion through `ExactNum` at a working precision sufficient to round correctly to binary32/binary64. The conversion path: `Ieee64::to_exact(p_wrk)` → apply `ExactNum` special → round back to `Ieee64`. No hardware float at any point.

Cover: `sqrt`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `exp`, `ln`, `erf`, `erfc`, `gamma`, `bessel_j_nu`, `bessel_y`, `bessel_i`, `bessel_k`.

Golds: `Ieee64Array::sin` matches `ExactNum::sin` rounded to binary64 for a representative set of values.

---

### 2.4 Integer SIMD for IEEE arrays

**Status:** partial 2026-08-30 — integer-lane add/mul (SSE2/NEON); no named `IEEE_SIMD_LANE_WIDTH`; no SIMD div/sqrt/fma.

**Prompt:**
Add SIMD-accelerated paths for `Ieee32Array` and `Ieee64Array` elementwise `add`, `sub`, `mul`, `div`, `sqrt`, and `fma` using integer SIMD lanes (`u32x8` / `u64x4` or equivalent via `std::simd` or `packed_simd2`). The arithmetic is software IEEE — integer lanes carrying bit patterns, arithmetic implemented in software, no hardware FPU instructions. The scalar and SIMD paths must produce bit-identical results.

Named constant: `IEEE_SIMD_LANE_WIDTH`. Use `cfg` feature gates so the scalar fallback compiles everywhere.

Golds: SIMD `add`/`mul`/`div`/`sqrt` on a 1000-element `Ieee64Array` produces bit-identical output to the scalar path. Results match `ExactNum` rounded to binary64.

---

## Section 3 — Certified interval arithmetic

### 3.1 Ball arithmetic for transcendentals

**Status:** done 2026-08-30 — `sin`/`cos`/`exp`/`ln`/`sqrt`/`erf`/`J0`/`J1`.

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

**Status:** done 2026-08-30 — `ComplexBall` disk add/mul/exp/ln/sin/cos.

**Prompt:**
Extend `Ball` to `ComplexBall { mid: ExactComplex, rad: ExactNum }` (a disk in ℂ). Implement `add`, `mul`, `exp`, `ln`, `sin`, `cos` with certified enclosures on `ComplexBall`. The radius grows by the Lipschitz constant of the operation on the disk.

Golds:
- `ComplexBall::exp` on the unit disk contains `e^{mid}`
- Composition of `ComplexBall::sin` and `ComplexBall::cos` on a small disk satisfies the Pythagorean identity enclosure

---

## Section 4 — Arbitrary-precision linear algebra

### 4.1 LU decomposition

**Status:** done 2026-08-30 — `lu_decomp`; \(PA=LU\); singular → `None`.

**Prompt:**
Implement LU decomposition with partial pivoting on `ExactNumArray` (2D, square or rectangular). Return `(L, U, P)` where `P` is a permutation vector. Use exact `ExactNum` arithmetic with explicit `(p, rm)` at each step. A singular matrix returns `None` for the `U` factor — do not invent a result.

Golds:
- `[[2,1],[4,3]]` → `L`, `U`, `P` such that `P·A = L·U` exactly
- Singular matrix `[[1,2],[2,4]]` → `None`
- `L·U` reconstructs `P·A` to working precision

---

### 4.2 QR decomposition (Gram-Schmidt)

**Status:** done 2026-08-30 — `qr_decomp`; \(QR=A\); \(Q^\top Q=I\); rank-deficient zero diagonal.

**Prompt:**
Implement QR decomposition via modified Gram-Schmidt on `ExactNumArray`. Return `(Q, R)` where `Q` is orthogonal and `R` is upper triangular, all at explicit precision `(p, rm)`. Rank-deficient input: zero column in `R` at the deficient position, not a panic.

Golds:
- `Q·R` reconstructs `A` to working precision
- `Q^T · Q = I` to working precision
- Rank-deficient matrix handled correctly

---

### 4.3 SVD (Golub-Reinsch)

**Status:** done 2026-08-30 — `svd_decomp`; \(U\Sigma V^T=A\); \(U^\top U=V^\top V=I\); \(\operatorname{diag}(3,2)\to\sigma=(3,2)\). Cap `SVD_ITER_MAX=64` sweeps per value.

**Prompt:**
Implement SVD on `ExactNumArray` via the Golub-Reinsch bidiagonalization algorithm at explicit precision `(p, rm)`. Return `(U, Σ, V^T)`. Singular values in descending order. Software limbs only — no LAPACK, no hardware float.

Cap: `SVD_ITER_MAX` named constant for the QR iteration convergence check. Return `None` if convergence is not reached within the cap.

Golds:
- `U · Σ · V^T` reconstructs `A` to working precision
- `U^T · U = I`, `V^T · V = I` to working precision
- Known singular values of `[[3,0],[0,2]]` are `3` and `2`

---

### 4.4 Eigenvalue decomposition

**Status:** done 2026-08-30 — `eigen_decomp`; \(Av=\lambda v\); \(V\Lambda V^T=A\); \(\lambda=(3,1)\); non-symmetric `None`. Cap `EIGEN_ITER_MAX=64`.

**Prompt:**
Implement eigenvalue decomposition for real symmetric matrices via the symmetric QR algorithm (tridiagonalization + QR iteration) at explicit precision `(p, rm)`. Return `(eigenvalues, eigenvectors)` as `(ExactNumArray, ExactNumArray)`. Non-symmetric input returns `None` — do not attempt a general eigendecomposition here.

Cap: `EIGEN_ITER_MAX` named constant.

Golds:
- `[[2,1],[1,2]]` eigenvalues are `1` and `3`; eigenvectors are orthogonal
- `A · v = λ · v` for each eigenpair to working precision
- `V · Λ · V^T = A` reconstruction to working precision

---

### 4.5 Multiprecision FFT

**Status:** done 2026-08-30 — `fft` / `ifft`; impulse, cosine bins, IFFT round-trip, Parseval. Cap `FFT_MAX_POINTS=4096`.

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

**Status:** done 2026-08-30 — `# Precision` on `ExactNum` and `ExactComplex` specials.

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

**Status:** done 2026-08-30 — `doc/REPRODUCIBILITY.md`. Locked values are unit tests, not a `golds/` tree.

**Prompt:**
Add a `doc/REPRODUCIBILITY.md` to zenith-float with the following content:

- Statement that all results are determined by the source code, input values, precision `p`, and rounding mode `rm` — no platform-dependent behavior, no hardware float
- Instructions for running the full gold suite: `cargo test --features mpfr-tests`
- Statement that `ZENITH_TEST_SEED` allows replay of any random test failure
- Statement that the gold values in `golds/` are locked expected outputs, not tolerances
- Citation format: "Computed with zenith-float {version}, reproducible by running `cargo test` at commit {hash}"

This document is for researchers who need to cite numerical results in papers.

---

## Section 8 — SKIPPEDS INTENTIONALLY

---

## Section 9 — Number representation completeness

### 9.1 Rational arithmetic (`ExactRational`)

**Prompt:**
Implement `ExactRational { num: ExactNum, den: ExactNum }` as an exact rational type where both numerator and denominator are `ExactNum` values. This is distinct from `ExactNum` which is a floating-point type — `ExactRational` is exact.

- `new(num, den)` — reduces to lowest terms; `den < 0` normalizes sign to numerator
- `from_i64(n, d)` — from integer ratio
- `add`, `sub`, `mul`, `div` — exact rational arithmetic
- `to_exact_num(p, rm)` — convert to `ExactNum` at given precision
- `is_integer()` — true if denominator is 1
- `floor()`, `ceil()`, `round()` — exact integer parts
- `partial_cmp` — exact comparison without conversion to float

Golds:
- `ExactRational::new(1, 3).add(ExactRational::new(1, 6)) = ExactRational::new(1, 2)`
- `ExactRational::new(2, 4) = ExactRational::new(1, 2)` (reduction)
- `ExactRational::new(1, 3).to_exact_num(256, ToEven)` matches `ExactNum::parse("0.333...")`
- `ExactRational::new(-3, -4) = ExactRational::new(3, 4)` (sign normalization)

---

### 9.2 Arbitrary-precision integer (`ExactInt`)

**Prompt:**
Implement `ExactInt` — an arbitrary-precision signed integer using the existing limb infrastructure. This is not `ExactNum` truncated — it is a proper big-integer type.

- `from_i64`, `from_u64`, `from_i128`, `from_u128`
- `add`, `sub`, `mul` — exact
- `div_rem(a, b)` → `(quotient, remainder)` — truncated division
- `gcd(a, b)` — Euclidean
- `pow(base, exp: u64)` — binary exponentiation
- `to_exact_num(p, rm)` — convert to floating point
- `from_exact_num(x)` — truncate a finite `ExactNum` to integer; `None` if not finite
- `bit_length()` — number of bits in the representation
- `is_zero()`, `is_negative()`, `signum()`

Golds:
- `ExactInt::from_i64(factorial(20))` matches `2432902008176640000`
- `gcd(ExactInt::from_i64(48), ExactInt::from_i64(18)) = 6`
- `pow(ExactInt::from_i64(2), 100)` matches known value
- `div_rem(ExactInt::from_i64(17), ExactInt::from_i64(5)) = (3, 2)`

---

### 9.3 Decimal string I/O at full precision

**Prompt:**
Extend `ExactNum::parse` and `ExactNum::format` to handle:
- Scientific notation input: `1.23e-45`, `1.23E+100`
- Arbitrary-precision decimal strings with more digits than `SOFT_PREC`
- Binary, octal, hex input with exact conversion
- `parse_exact(s)` — parse to the minimum precision needed to represent the value exactly (for terminating decimals)
- `format_exact(radix)` — format to the minimum digits needed for exact round-trip

Golds:
- `parse_exact("0.1")` is exactly `1/10` — verify `10 * parse_exact("0.1") == 1` exactly
- `parse_exact("0.5")` is exactly representable in binary — bit pattern gold
- `format_exact(Radix::Dec, parse_exact("0.125"))` = `"0.125"` (exact round-trip)
- Scientific notation: `parse("1.5e3", Dec, 256, ToEven)` = `1500` exactly

---

## Section 10 — Extended precision infrastructure

### 10.1 Per-operation precision tracking

**Status:** done 2026-08-30 — specials take `(p, rm, cc)`; no `SOFT_PREC` in the kernel.

**Prompt:**
Implement per-operation precision at caller-chosen bits for ALL unary special functions — not global `SOFT_PREC`. This is the near-term item that has been open since the beginning.

Every special function (`erf`, `erfc`, `gamma`, `ln_gamma`, `digamma`, `gammainc`, `ei`, `si`, `ci`, `li`, `fresnel_s`, `fresnel_c`, `bessel_j_nu`, `bessel_y`, `bessel_i`, `bessel_k`, `elliptic_k`, `elliptic_e_complete`, `elliptic_f`, `elliptic_e`, `elliptic_pi_complete`, `elliptic_pi`, `legendre_p`, `assoc_legendre_p`, `hypergeom_2f1`, `betainc`, `ai`, `bi`) must accept `(p: usize, rm: RoundingMode, cc: &mut Consts)` and use that precision throughout — not a global constant.

Named constant `SOFT_PREC` must only be used as a default when the caller does not specify precision — it must not be hardcoded inside any kernel function.

Golds:
- `erf(1.0, 64, ToEven)` correct to 64 bits
- `erf(1.0, 512, ToEven)` correct to 512 bits — different result, more digits
- `bessel_j_nu(0.5, 1.0, 256, ToEven)` matches MPFR oracle at 256 bits
- No global `SOFT_PREC` reference inside any kernel function — verified by grep

---

### 10.2 Precision propagation through `expr!`

**Prompt:**
Extend `expr!` and `cexpr!` to propagate precision through special function leaves correctly. Currently the context precision `p` may not be correctly passed to all leaf functions.

Verify and fix:
- Every `expr!` leaf calls its underlying function at `p_wrk` (working precision), not at a fixed precision
- `cexpr!` leaves for complex specials use `p_wrk` for both real and imaginary parts
- The final `set_precision(p, rm)` is applied after all leaves are evaluated at `p_wrk`
- `errs[]` slots correctly detect cancellation in special function results

Golds:
- `expr!(erf(x) + erfc(x), ctx)` at 256 bits = `1` to full precision (no cancellation loss)
- `expr!(bessel_j(0, x)^2 + bessel_y(0, x)^2, ctx)` computes at working precision throughout
- `cexpr!(erf(z), ctx)` for complex `z` uses `p_wrk` inside the Faddeeva computation

---

### 10.3 `ziv_round` extension for vector operations

**Prompt:**
Extend `ziv_round` to support vectorized correct-rounding: given a function `f: &[ExactNum] -> ExactNum`, find the precision `p_wrk` such that the result rounds uniquely to `p` bits.

`ziv_round_vec(p, rm, inputs, compute)` — same Ziv loop as `ziv_round` but `compute` takes a slice of inputs all evaluated at `p_wrk`.

This is needed for functions like `hypot(x, y)` and `atan2(y, x)` where all inputs must be at the same working precision.

Golds:
- `ziv_round_vec` for `hypot(3, 4)` returns `5` exactly at any precision
- `ziv_round_vec` for `atan2(1, 1)` returns `π/4` correctly rounded at 256 bits
- `MAX_PREC_RETRY` still bounds the retry count

---

## Section 11 — Matrix and linear algebra

Duplicate of Section 4. Do not implement twice. Status is §4.1–§4.5.

### 11.1 LU decomposition for `ExactNumArray`

**Prompt:**
Implement LU decomposition with partial pivoting on `ExactNumArray` (2D square). Return `(L, U, P)` where `P` is a permutation vector. Use exact `ExactNum` arithmetic with explicit `(p, rm)` at each step. Singular matrix → `None`.

Golds:
- `[[2,1],[4,3]]` → `L`, `U`, `P` such that `P·A = L·U` to working precision
- Singular `[[1,2],[2,4]]` → `None`
- `L·U` reconstructs `P·A` to working precision
- `L` is unit lower triangular, `U` is upper triangular

---

### 11.2 QR decomposition

**Prompt:**
Implement QR decomposition via modified Gram-Schmidt on `ExactNumArray`. Return `(Q, R)` at explicit `(p, rm)`. Rank-deficient: zero column in `R`, not a panic.

Golds:
- `Q·R` reconstructs `A` to working precision
- `Q^T · Q = I` to working precision
- Rank-deficient matrix: zero column in `R` at the deficient position

---

### 11.3 SVD

**Status:** done 2026-08-30 — same method as §4.3 (`svd_decomp`). Do not implement twice.

**Prompt:**
Implement SVD via Golub-Reinsch bidiagonalization on `ExactNumArray` at explicit `(p, rm)`. Return `(U, Σ, V^T)`. Singular values in descending order. Cap `SVD_ITER_MAX` named constant. `None` if not converged.

Golds:
- `U · Σ · V^T` reconstructs `A` to working precision
- `U^T · U = I`, `V^T · V = I` to working precision
- Known singular values of `[[3,0],[0,2]]` are `3` and `2`

---

### 11.4 Eigenvalue decomposition

**Status:** done 2026-08-30 — same method as §4.4 (`eigen_decomp`). Do not implement twice.

**Prompt:**
Implement eigenvalue decomposition for real symmetric matrices via symmetric QR algorithm at explicit `(p, rm)`. Return `(eigenvalues, eigenvectors)`. Non-symmetric → `None`. Cap `EIGEN_ITER_MAX`.

Golds:
- `[[2,1],[1,2]]` eigenvalues are `1` and `3`; eigenvectors orthogonal
- `A · v = λ · v` for each eigenpair to working precision
- `V · Λ · V^T = A` reconstruction

---

### 11.5 Multiprecision FFT

**Status:** done 2026-08-30 — same methods as §4.5 (`fft` / `ifft`). Do not implement twice.

**Prompt:**
Implement radix-2 Cooley-Tukey FFT on complex `ExactNumArray` (length power of 2) at explicit `(p, rm, cc)`. Twiddle factors from `ExactNum::sin_cos`. Cap `FFT_MAX_POINTS`.

Golds:
- FFT of `[1, 0, 0, 0]` is `[1, 1, 1, 1]`
- FFT of pure cosine at frequency `k` has energy only at bins `k` and `N-k`
- IFFT(FFT(x)) = x to working precision
- Parseval: `Σ|x_n|² = (1/N)Σ|X_k|²` to working precision

---

## Section 12 — Statistics and distributions in zenith-float

### 12.1 Statistical distribution kernels

**Prompt:**
Implement the following distribution PDF/CDF kernels at arbitrary precision in zenith-float. These are the numeric primitives that Accumath's symbolic layer calls for numeric evaluation.

- `normal_pdf(x, mu, sigma, p, rm, cc)` — `exp(-(x-mu)²/(2sigma²)) / (sigma·sqrt(2π))`
- `normal_cdf(x, mu, sigma, p, rm, cc)` — `(1 + erf((x-mu)/(sigma·sqrt(2))))/2`; uses existing `erf`
- `gamma_pdf(x, alpha, beta, p, rm, cc)` — `x^(alpha-1) · exp(-x/beta) / (beta^alpha · Γ(alpha))`
- `beta_pdf(x, alpha, beta, p, rm, cc)` — `x^(alpha-1)(1-x)^(beta-1) / B(alpha,beta)`; uses `betainc`
- `poisson_pmf(k, lambda, p, rm)` — `lambda^k · exp(-lambda) / k!` for integer `k`
- `binomial_pmf(k, n, prob, p, rm)` — `C(n,k) · prob^k · (1-prob)^(n-k)`
- `chi_squared_cdf(x, k, p, rm, cc)` — via regularized incomplete gamma `γ(k/2, x/2)/Γ(k/2)`
- `student_t_pdf(x, nu, p, rm, cc)` — via `Γ` and `betainc`

Golds:
- `normal_pdf(0, 0, 1)` = `1/sqrt(2π)` to 256 bits
- `normal_cdf(0, 0, 1)` = `1/2` exactly
- `gamma_pdf(1, 1, 1)` = `exp(-1)` to 256 bits
- `poisson_pmf(0, 1)` = `exp(-1)` to 256 bits
- `chi_squared_cdf` at known quantiles matches standard tables to 10 decimal places

---

### 12.2 Random number generation at arbitrary precision

**Prompt:**
Extend the existing `random` feature to support:
- `ExactNum::random_uniform(a, b, p, rm)` — uniform on `[a, b]` at precision `p`
- `ExactNum::random_normal(mu, sigma, p, rm, cc)` — normal via Box-Muller at precision `p`
- `ExactNum::random_exponential(lambda, p, rm)` — exponential via inverse CDF
- `ExactNumArray::random_fill(shape, dist, p, rm, cc)` — fill a 2D array with samples
- Deterministic seeding: `reseed_random(seed)` makes all subsequent samples reproducible

Named constant: `DEFAULT_RANDOM_SEED` (already exists — unchanged).

Golds:
- `random_uniform(0, 1, 256)` produces values in `[0, 1]`
- With same seed, two calls to `random_normal` produce identical sequences
- `random_exponential(1.0)` sample mean over 10000 samples is within 0.05 of 1.0
- `ExactNumArray::random_fill((100, 100), Normal(0,1))` has correct shape

---

## Section 13 — Polynomial arithmetic

### 13.1 Dense polynomial arithmetic

**Prompt:**
Implement `ExactNumPoly { coeffs: Vec<ExactNum>, p: usize, rm: RoundingMode }` — a univariate polynomial with `ExactNum` coefficients at explicit precision.

- `eval(x, p, rm, cc)` — Horner evaluation
- `add`, `sub`, `mul` — exact at precision `p`
- `div_rem(f, g, p, rm)` → `(quotient, remainder)`
- `gcd(f, g, p, rm)` — Euclidean with content removal
- `compose(f, g, p, rm, cc)` — `f(g(x))`
- `derivative(p, rm)` — formal derivative
- `integral(p, rm)` — indefinite integral (no constant)
- `roots_real(p, rm, cc)` — numerical real roots via companion matrix + eigenvalues

Golds:
- `(x² - 1).div_rem(x - 1) = (x + 1, 0)`
- `gcd(x² - 1, x - 1) = x - 1`
- `compose(x², x + 1) = x² + 2x + 1`
- `roots_real(x² - 2)` contains `±√2` to working precision
- `derivative(x³) = 3x²`

---

### 13.2 Chebyshev approximation

**Prompt:**
Implement Chebyshev polynomial approximation:
- `chebyshev_coeffs(f, n, a, b, p, rm, cc)` — compute `n` Chebyshev coefficients of `f` on `[a,b]` at precision `p`
- `chebyshev_eval(coeffs, x, a, b, p, rm)` — evaluate the Chebyshev expansion
- `chebyshev_error_bound(coeffs, p)` — bound on the truncation error from the tail coefficients
- `clenshaw(coeffs, x, p, rm)` — Clenshaw recurrence for stable evaluation

Named constant: `CHEBYSHEV_MAX_DEGREE`.

Golds:
- Chebyshev approximation of `exp(x)` on `[-1,1]` with 20 terms: error < `10^{-15}` at 256 bits
- `chebyshev_eval` at the nodes equals the function exactly (up to rounding)
- `clenshaw` matches direct evaluation of the Chebyshev polynomial to working precision

---

### 13.3 Orthogonal polynomials

**Prompt:**
Implement the following orthogonal polynomial families as `ExactNum` computations:

- `hermite_he(n, x, p, rm)` — probabilist's Hermite polynomial `He_n(x)` via recurrence
- `hermite_h(n, x, p, rm)` — physicist's Hermite polynomial `H_n(x)` via recurrence
- `laguerre(n, x, p, rm)` — Laguerre polynomial `L_n(x)`
- `gen_laguerre(n, alpha, x, p, rm)` — generalized Laguerre `L_n^(α)(x)`
- `chebyshev_t(n, x, p, rm)` — Chebyshev polynomial of first kind `T_n(x)`
- `chebyshev_u(n, x, p, rm)` — Chebyshev polynomial of second kind `U_n(x)`
- `gegenbauer(n, lambda, x, p, rm)` — Gegenbauer (ultraspherical) polynomial

Named cap: `ORTHOPOLY_N_MAX` for all of the above.

Golds:
- `hermite_he(4, 0) = 3` (known value)
- `laguerre(3, 0) = 1` (L_n(0) = 1 for all n)
- `chebyshev_t(5, cos(π/5)) = cos(π)= -1` (Chebyshev property)
- `gegenbauer(2, 1, x) = 3x² - 1` (matches Legendre for λ=1 up to normalization)
- Recurrence: `T_{n+1}(x) = 2x T_n(x) - T_{n-1}(x)` verified at `n=5`

---

## Section 14 — Numerical methods

### 14.1 Numerical integration (quadrature)

**Prompt:**
Implement high-precision numerical quadrature:
- `gauss_legendre(f, a, b, n, p, rm, cc)` — Gauss-Legendre quadrature with `n` nodes; nodes and weights computed at precision `p` from Legendre polynomial roots
- `tanh_sinh(f, a, b, p, rm, cc)` — tanh-sinh (double exponential) quadrature; excellent for integrands with endpoint singularities
- `gauss_laguerre(f, n, p, rm, cc)` — Gauss-Laguerre for `∫₀^∞ f(x) e^{-x} dx`
- `gauss_hermite(f, n, p, rm, cc)` — Gauss-Hermite for `∫_{-∞}^∞ f(x) e^{-x²} dx`

Named constants: `QUADRATURE_MAX_NODES`, `TANH_SINH_LEVELS_MAX`.

Golds:
- `gauss_legendre(|x| x², -1, 1, 10)` = `2/3` to 256 bits
- `tanh_sinh(|x| 1/sqrt(1-x²), -1, 1, ...)` = `π` to 256 bits (singular integrand)
- `gauss_laguerre(|x| x², 10)` = `Γ(3) = 2` to 256 bits
- 20-point Gauss-Legendre integrates degree-39 polynomials exactly

---

### 14.2 Root finding at arbitrary precision

**Prompt:**
Implement root finding algorithms at arbitrary precision:
- `bisect(f, a, b, tol, p, rm, cc)` — bisection; guaranteed convergence; `None` if `f(a)f(b) > 0`
- `newton(f, df, x0, tol, max_iter, p, rm, cc)` — Newton-Raphson with exact derivative
- `brent(f, a, b, tol, p, rm, cc)` — Brent's method; combines bisection + secant + inverse quadratic
- `illinois(f, a, b, tol, p, rm, cc)` — Illinois algorithm; superlinear convergence

Named constants: `ROOT_MAX_ITER`, `ROOT_DEFAULT_TOL`.

Golds:
- `bisect(sin, 3, 4)` finds `π` to 256 bits
- `newton(|x| x²-2, |x| 2x, 1.0)` finds `√2` to 256 bits in a few iterations
- `brent(sin, 3, 4)` finds `π` faster than bisection
- `bisect(sin, 0, 1)` → `None` (no sign change)

---

### 14.3 ODE solvers at arbitrary precision

**Prompt:**
Implement ODE solvers at arbitrary precision for `ExactNum`:
- `rk4(f, t0, y0, t1, n_steps, p, rm, cc)` — classical RK4; fixed step; returns `(t_values, y_values)` as `ExactNumArray`
- `rk45_adaptive(f, t0, y0, t1, atol, rtol, p, rm, cc)` — Dormand-Prince RK5(4) at arbitrary precision; step control based on `ExactNum` comparison
- `euler(f, t0, y0, t1, n_steps, p, rm)` — explicit Euler; educational use

Named constants: `ODE_MAX_STEPS`, `ODE_MIN_STEP`.

Golds:
- `rk4(|t,y| -y, 0, 1, 1, 1000)` matches `exp(-1)` to 50 decimal places at 256 bits
- `rk45_adaptive` for `y' = -y` meets `atol = 1e-50` at 256 bits
- `euler` with 1000 steps for `y' = y` has error bounded by `O(h)` — verified

---

## Section 15 — Signal processing primitives

### 15.1 Discrete transforms

**Prompt:**
Implement discrete signal processing transforms at arbitrary precision:
- `dct(signal, p, rm, cc)` — Discrete Cosine Transform (Type II) via FFT
- `idct(signal, p, rm, cc)` — Inverse DCT
- `dst(signal, p, rm, cc)` — Discrete Sine Transform
- `idst(signal, p, rm, cc)` — Inverse DST
- `fft_real(signal, p, rm, cc)` — FFT of a real signal (returns complex output)
- `ifft_real(spectrum, p, rm, cc)` — IFFT back to real

Golds:
- `idct(dct(x)) = x` to working precision
- `dct` of a constant signal has energy only in DC bin
- `fft_real` of a pure cosine at frequency `k` has energy at bins `k` and `N-k` only
- Parseval: `Σ|x_n|² = (1/N)Σ|X_k|²`

---

### 15.2 Window functions

**Prompt:**
Implement windowing functions as `ExactNumArray` generators:
- `hann_window(n, p, rm, cc)` — Hann window: `0.5(1 - cos(2πk/N))`
- `hamming_window(n, p, rm, cc)` — Hamming: `0.54 - 0.46cos(2πk/N)`
- `blackman_window(n, p, rm, cc)` — Blackman: three-term cosine
- `kaiser_window(n, beta, p, rm, cc)` — Kaiser: uses `I_0` from the Bessel functions
- `rectangular_window(n, p, rm)` — all ones (identity)

Golds:
- `hann_window(4)` = `[0, 0.75, 0.75, 0]` (known values)
- `hamming_window` endpoints are `0.08` not `0` (distinguishes from Hann)
- `kaiser_window` with `beta=0` equals rectangular window
- All windows sum to positive values — no negative sums

---

## Section 16 — Cryptographic primitives in zenith-float

### 16.1 Modular arithmetic on `ExactInt`

**Prompt:**
Implement modular arithmetic on `ExactInt` (from §9.2):
- `mod_pow(base, exp, modulus)` — fast modular exponentiation; `exp` is `ExactInt`
- `mod_inv(a, modulus)` — extended Euclidean; `None` if not invertible
- `miller_rabin(n, witnesses)` — primality test; deterministic for `n < 3·10^{18}`
- `pollard_rho(n)` — Brent's variant; `None` if `n` is prime; cap `POLLARD_RHO_ITER_MAX`

These are the number-theoretic primitives that the Accumath number theory layer calls.

Golds:
- `mod_pow(2, 100, 1000000007)` matches known value
- `mod_inv(3, 7) = 5`
- `miller_rabin(2^{31} - 1, [2,3,5,7])` = `true` (Mersenne prime)
- `pollard_rho(8051)` finds a factor (known: `83 × 97`)

---

### 16.2 Hash functions

**Prompt:**
Implement cryptographic hash functions operating on byte arrays, needed for license key validation and signed package verification:
- `sha256(data: &[u8]) -> [u8; 32]` — pure Rust, no external crates
- `sha512(data: &[u8]) -> [u8; 64]`
- `hmac_sha256(key: &[u8], data: &[u8]) -> [u8; 32]`
- `constant_time_eq(a: &[u8], b: &[u8]) -> bool` — timing-safe comparison

These are pure integer operations — no floating point involved.

Golds:
- `sha256(b"")` matches the known SHA-256 of the empty string
- `sha256(b"abc")` matches the known test vector
- `hmac_sha256` matches NIST test vectors
- `constant_time_eq` returns the same result regardless of where the first difference is

---

## Section 17 — Serialization and interchange

### 17.1 Serde support for all types

**Prompt:**
Extend serde support (already present for `ExactNum` and `ExactComplex`) to all new types:
- `ExactRational` — serializes as `{"num": "...", "den": "..."}` decimal strings
- `ExactInt` — serializes as a decimal string
- `ExactNumArray` — serializes as `{"shape": [rows, cols], "data": ["...", ...], "p": n, "rm": "ToEven"}`
- `Ieee32Array`, `Ieee64Array` — serialize as `{"shape": [rows, cols], "data": [f32/f64 values]}`
- `Ball` — serializes as `{"mid": "...", "rad": "..."}`

All deserializers reject mismatched shapes and return `Err`. Round-trip: `deserialize(serialize(x)) == x` for all types.

Golds:
- `ExactNumArray` round-trip through JSON preserves all values to full precision
- `Ieee64Array` round-trip through JSON preserves all bit patterns
- `ExactRational` round-trip: `1/3` serializes and deserializes correctly
- Mismatched shape in `ExactNumArray` JSON returns `Err`

---

### 17.2 Binary format

**Prompt:**
Implement a compact binary format for `ExactNum` and `ExactNumArray` that is more efficient than JSON for large arrays:
- `ExactNum::to_bytes() -> Vec<u8>` — little-endian limb array with sign, exponent, precision header
- `ExactNum::from_bytes(bytes: &[u8]) -> Result<Self, Error>`
- `ExactNumArray::to_bytes() -> Vec<u8>` — shape header + packed limb arrays
- `ExactNumArray::from_bytes(bytes: &[u8]) -> Result<Self, Error>`
- Version byte in the header for forward compatibility

Golds:
- Round-trip: `from_bytes(to_bytes(x)) == x` for all finite `ExactNum` values
- `to_bytes(ExactNum::NAN)` produces a recognizable NaN encoding
- Array round-trip preserves shape and all element values
- Invalid bytes return `Err` not a panic

---

### 17.3 HDF5 and CSV I/O for arrays

**Prompt:**
Implement I/O for `ExactNumArray` and `Ieee64Array`:
- `Ieee64Array::from_csv(path) -> Result<Self, Error>` — parse CSV; missing → `NaN`
- `Ieee64Array::to_csv(path) -> Result<(), Error>` — write CSV at full `f64` precision
- `ExactNumArray::from_hdf5(path, dataset) -> Result<Self, Error>` via `hdf5-rust` crate
- `ExactNumArray::to_hdf5(path, dataset) -> Result<(), Error>`
- `Ieee64Array::from_hdf5` / `to_hdf5` — same via native HDF5 double type

Named constant: `CSV_MAX_ROWS`.

Golds:
- CSV round-trip: write 100×3 `Ieee64Array`, read back, values match
- HDF5 round-trip: write `ExactNumArray`, read back at same precision, values match
- Missing value in CSV → `NaN` at that position
- Wrong dataset name in HDF5 → `Err`

---

## Section 18 — Documentation completeness

### 18.1 Per-function precision documentation

**Status:** done 2026-08-30 — same `# Precision` comments as §5.1. Do not write twice.

**Prompt:**
For every special function in zenith-float add a `# Precision` section to the doc comment stating:
- Algorithm used (series, AGM, Carlson, asymptotic)
- The argument range where each algorithm applies
- Named threshold constants
- Error bound in ULP for the primary region
- Whether MPFR oracle exists and what domain it covers
- Known slow cases (e.g. large `|e|` for trig, near-zero for `Ci`)

This is documentation only — no code changes. Format follows the existing `erf` doc comment as a template.

Golds (documentation audit):
- Every function in §12 of `LIBRARY.md` has a `# Precision` section
- Every named threshold constant is referenced in the doc comment for its function
- No function says "arbitrary precision" without specifying the algorithm

---

### 18.2 `GETTING_STARTED.md` update

**Prompt:**
Update `doc/GETTING_STARTED.md` to cover all new capabilities added since the original version:
- Software `Ieee32` / `Ieee64` — how to use them for interchange
- 2D array operations — `from_shape`, `matmul`, elementwise ops
- Special functions at caller-chosen precision — the `(p, rm, cc)` pattern
- `ExactRational` and `ExactInt` — when to use each
- Complex specials via `cexpr!` — branch cut conventions
- `Ball` interval arithmetic — what it guarantees and what it doesn't

Each new section follows the pattern of the existing sections: prose explanation, working code example, link to the detailed doc.

---

### 18.3 `HELP.md` complete rewrite

**Prompt:**
Rewrite `doc/HELP.md` as a comprehensive user guide — not a reference (that is `LIBRARY.md`) and not a tutorial (that is `GETTING_STARTED.md`) but an explanation of why the API looks the way it does, with recipes for common tasks and a complete FAQ.

Structure:
- **The precision model** — why you specify `p` on every operation and what happens if you don't
- **Rounding modes** — when to use each; why `ToEven` is the default
- **Constants cache** — why `Consts` exists; how to reuse it; `SharedConsts` for threads
- **`expr!` vs methods** — when to use each; what `expr!` guarantees and what it doesn't
- **Complex arithmetic** — branch cuts; why `cexpr!` has different leaves than `expr!`
- **Special functions** — algorithm overview; what to do when `PrecisionRetryExhausted`
- **Arrays** — shape semantics; why elementwise requires matching shapes
- **Software IEEE** — what it is; when to use `Ieee64` vs `ExactNum`
- **Recipes** — 20+ worked examples for common numerical tasks
- **Common mistakes** — the 15 most common errors and how to fix them
- **FAQ** — 30+ questions with direct answers

---

## Section 19 — CI and quality infrastructure

### 19.1 MPFR oracle extension

**Prompt:**
Extend the `mpfr-tests` feature to cover all new special functions:
- Complex `erf`/`erfc` — MPFR has `mpc_erf` via the MPC library
- Complex `Γ`/`ln_gamma` — via MPC
- Real and complex Airy `Ai`/`Bi` — MPFR has `mpfr_ai`
- Complex Bessel — via MPC/ARB
- `ExactRational` arithmetic — compare with GMP rational

For functions where MPFR/MPC has no equivalent (`si`, `ci`, `li`, Fresnel, complex elliptic, `₂F₁`): use identity/series golds instead. Document clearly which category each function falls in.

Golds:
- Every new function either has an MPFR oracle gold or an explicitly documented identity gold
- The `mpfr-tests` CI pass rate stays at 100%

---

### 19.2 Property-based testing with `proptest`

**Prompt:**
Add `proptest` integration for key properties that should hold for all inputs:
- `ExactNum` arithmetic: `a + b = b + a` (commutativity) at random precisions
- `ExactNum` rounding monotonicity: rounding at precision `p` then `q < p` equals rounding at `q` directly (when rounding mode is the same)
- Special functions: `erf(-x) = -erf(x)` for all real `x` (odd symmetry)
- `ExactNumArray` matmul: `(AB)C = A(BC)` for random matrices at fixed precision
- `ExactRational` arithmetic: `(a + b) - b = a` for all rational `a`, `b`

Named constant: `PROPTEST_CASES = 1000`.

Golds:
- All property tests pass with `PROPTEST_CASES = 1000` random inputs
- Any failure produces a minimal counterexample via shrinking
- Property tests run in `cargo test` without the `mpfr-tests` feature

---

### 19.3 Benchmark suite

**Prompt:**
Extend `scripts/bench.sh` to cover all new capabilities:
- Special functions at 64, 128, 256, 512, 1024 bits — one benchmark per function family
- Array matmul at various sizes: `(10×10)`, `(100×100)`, `(1000×1000)`
- FFT at sizes `256`, `1024`, `4096`
- LU decomposition for `(50×50)`, `(200×200)`
- Comparison against `astro-float` and `dashu-float` for the functions they implement

`scripts/compare-bench.sh --quick` must run in under 5 minutes. Full benchmarks may take longer.

Golds:
- Benchmark results recorded in `doc/bench-baselines.tsv`
- No regression of more than 10% in any existing benchmark without a documented reason
- New benchmarks added to `doc/compare-results.tsv` for library comparison

---

### 19.4 Pre-publish checklist

**Prompt:**
Create `scripts/zenith_prepublish.sh` that runs and fails loudly on any error:

1. `cargo test --all-features` — all tests pass
2. `cargo test --features mpfr-tests` — MPFR oracle golds pass
3. `cargo build --no-default-features` — no_std compiles
4. `cargo build --no-default-features --target thumbv7em-none-eabihf` — embedded target compiles
5. `scripts/compare-bench.sh --quick` — benchmark baseline refreshed
6. Kernel purity grep: no hardware float type tokens in `zenith-float-num/src/` — software `Ieee32`/`Ieee64` permitted
7. Every function in §12 of `LIBRARY.md` has a `# Precision` doc section — `scripts/check_precision_docs.sh`
8. `ZENITH_FLOAT_CAPABILITIES.md` version and date match `Cargo.toml`
9. `doc/LIBRARY.md` `expr!` and `cexpr!` leaf lists match the macro expansion — `scripts/check_leaves.sh`
10. All comparison table cells in `ZENITH_FLOAT_CAPABILITIES.md` §24 are verified — no `<!-- verify -->` markers
11. `proptest` property tests pass with `PROPTEST_CASES = 1000`
12. No `TODO`, `FIXME`, `HACK`, or `unwrap()` without proof comment in source

Script exits 0 only when all 12 checks pass.

---

## Section 20 — Reproducibility and citation support

### 20.1 Reproducibility documentation

**Status:** done 2026-08-30 — same file as §7.1. Do not write twice.

**Prompt:**
Write `doc/REPRODUCIBILITY.md` with the following content:

- Statement that all results are determined by: source code, input values, precision `p`, rounding mode `rm`, and (for constants) the `Consts` cache state — no platform-dependent behavior, no hardware float
- Instructions for running the full gold suite: `cargo test --features mpfr-tests`
- Statement that `ZENITH_TEST_SEED` allows replay of any random test failure
- Statement that the gold values in `golds/` are locked expected outputs, not tolerances
- Citation format: `"Computed with zenith-float {version} (commit {hash}), reproducible by running cargo test at that commit"`
- Note that software `Ieee32`/`Ieee64` results are also reproducible — they use integer arithmetic, not the hardware FPU
- Statement that results are bit-identical across platforms (x86_64, aarch64, wasm32) — verified by the hex limb CI job

---

### 20.2 Platform verification CI

**Prompt:**
Extend CI to verify bit-identical results across platforms:
- `scripts/ci_hex_arm.sh` — already exists for aarch64; extend to cover new special functions
- Add `scripts/ci_hex_wasm.sh` — run the hex limb gold suite compiled to WebAssembly via `wasmtime`
- Add `scripts/ci_hex_32bit.sh` — run on 32-bit target to verify `WORD_BIT_SIZE = 32` path

Each script: compile for target, run the hex gold suite, compare output byte-for-byte with x86_64 reference.

Golds:
- All three platform scripts exit 0 on the existing gold suite
- New special function hex golds added to `golds/hex/` for at least one representative value per new function

