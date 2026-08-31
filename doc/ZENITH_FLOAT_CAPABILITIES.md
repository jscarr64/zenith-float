# zenith-float Capability Reference

**Version:** 0.1.0  
**Date:** 2026-08-31  
**License:** MIT OR Apache-2.0  
**Status:** Public crate  

This document states exactly what `zenith-float` provides, what it does not provide, and what the named constants and caps mean in practice. It is not a tutorial. For a working example, see `doc/GETTING_STARTED.md`. For full method signatures, see `doc/LIBRARY.md`. For `expr!` rounding contracts, see `doc/EXPR.md`.

Engineering walk list: [`ZENITH_FLOAT_BUILD_PLAN.md`](ZENITH_FLOAT_BUILD_PLAN.md) (status table at the top). Maintainer CI inventory: [`BUILD_CHECKLIST.md`](BUILD_CHECKLIST.md). There is no TODO file — keep these three current after every slice.

A row marked ✅ is backed by a unit test or MPFR oracle gold. A row marked 🟡 has a named bound or condition; read the notes. Every `NaN` return has an associated `Error`; no operation silently discards a failure.

---

## How to read this document

| Symbol | Meaning |
| --- | --- |
| ✅ | Implemented, tested |
| 🟡 | Partial — named cap, condition, or platform constraint; read the notes |
| ⬜ | Not implemented |
| 🚫 | Explicitly out of scope; will not be added to this crate |

---

## 1. What this crate is

`zenith-float` is an arbitrary-precision software floating-point library. Every calculation uses integer limbs. Hardware floating-point arithmetic is not used. IEEE widths are software types `Ieee32` / `Ieee64` (`u32` / `u64` bits, `from_bits` / `to_bits`). Arbitrary precision is `ExactNum`.

Depend on `zenith-float`, not `zenith-float-num`. The kernel crate is an implementation detail.

---

## 2. Cargo features

| Feature | Default | What it enables |
| --- | --- | --- |
| `std` | yes | `Display`, `LowerExp`, `UpperExp`, `Binary`, `Octal`, `UpperHex`, `LowerHex`, `FromStr`, `std::error::Error` for `Error`, `SharedConsts`, serde when `serde` is also on |
| `random` | no | `ExactNum::random_normal`, `random_seed`, `reseed_random`, `seeded_random`, `DEFAULT_RANDOM_SEED` |
| `serde` | no | `Serialize` / `Deserialize` for `ExactNum` / `ExactComplex` / `ExactRational` / `ExactInt` / arrays / `Ball`; decimal strings carry `@p=`; implies `std` |
| `mpfr-tests` | no | Optional MPFR bit-oracle tests; Linux x86_64 + `rug` only; not a runtime dependency |

`no_std` is supported when a global allocator is available (`default-features = false`). Formatting traits require `std`.

---

## 3. Numeric model

| Item | Value |
| --- | --- |
| `Word` | `u64` on 64-bit targets; `u32` on 32-bit |
| `WORD_BIT_SIZE` | 64 or 32 accordingly |
| `INLINE_WORDS` | `2` — mantissas of ≤ 128 bits live on the stack |
| `EXPONENT_MAX` | `i32::MAX` on 64-bit; `i32::MAX / 4` on 32-bit |
| `EXPONENT_MIN` | `i32::MIN` on 64-bit; `i32::MIN / 4` on 32-bit |
| Default operator precision | 128 bits, `ToEven` — used by `+` `−` `×` `÷` operator traits only |
| `MAX_PREC_RETRY` | `256` — extra word-sized retry steps beyond `p` for correct-rounding loops |
| `SVD_ITER_MAX` | `64` — QR sweeps per singular value in `svd_decomp`; then `None` |
| `EIGEN_ITER_MAX` | `64` — QR sweeps per eigenvalue in `eigen_decomp`; then `None` |
| `FFT_MAX_POINTS` | `4096` — max length of `fft` / `ifft`; longer → `None` |
| `POLY_COMPANION_CLOSED_DEG` | `2` — `ExactNumPoly::roots_real` closed-form companion eigenvalues; higher degree → `None` |
| `CHEBYSHEV_MAX_DEGREE` | `256` — max coefficient count for `chebyshev_coeffs`; larger `n` → `None` |
| `ORTHOPOLY_N_MAX` | `256` — max degree for Hermite / Laguerre / Chebyshev T,U / Gegenbauer; larger `n` → `NaN` |
| `QUADRATURE_MAX_NODES` | `64` — max Gauss nodes; larger `n` → `None` |
| `TANH_SINH_LEVELS_MAX` | `8` — max tanh–sinh step halvings after `h = 2π/(p ln 2)` |
| `ROOT_MAX_ITER` | `256` — cap on bisection / Newton / Brent / Illinois steps |
| `ROOT_DEFAULT_TOL` | `−256` — default absolute tolerance exponent (`2^{ROOT_DEFAULT_TOL}`) |
| `ODE_MAX_STEPS` | `65536` — accepted-step cap for RK4 / RK45 / Euler |
| `ODE_MIN_STEP` | `−256` — minimum RK45 step exponent (`h_min = 2^{ODE_MIN_STEP}`) |
| `DSP_MAX_POINTS` | `2048` — max real length for `dct`/`idct`/`dst`/`idst` (`2N`-point FFT) |
| `IEEE_SIMD_LANE_WIDTH` | `4` — `u32` lanes per integer SIMD vector; binary64 uses 2 `u64` lanes |
| `BINARY_INLINE_LEN` | `16` — stack record: flags, version, `n_sig`, inexact, `i32` exponent BE, two `u32` limbs BE |
| `BINARY_INLINE_MANT_BITS` | `64` — max mantissa bits for [`ExactNum::to_inline_bytes`]; wider → `MemoryAllocation` |
| `BINARY_FORMAT_VERSION` | `1` — first version byte of every record |
| `BINARY_MAX_U32` | `65536` — max `u32` limbs in a heap record |
| `BINARY_MAX_ELEMS` | `1048576` — max array elements in `from_bytes` |
| `CSV_MAX_ROWS` | `1048576` — max data rows in `from_csv_str` |
| `CSV_MAX_COLS` | `4096` — max columns in one CSV row |
| `PROPTEST_CASES` | `1000` — cases per property in `tests/proptest_props.rs` |
| `POLLARD_RHO_ITER_MAX` | `1048576` — `f` evaluations per `c` in Brent Pollard ρ |

**Special values:** `+Inf`, `−Inf`, `NaN` (with optional `Error`), subnormals at `EXPONENT_MIN`. Public sentinels: `INF_POS`, `INF_NEG`, `NAN`.

**Precision rounding:** requested bit count is always rounded up to the next word boundary. Asking for 1 bit gives 64 bits on a 64-bit target.

---

## 4. Rounding modes

| Mode | Meaning |
| --- | --- |
| `None` | Skip the final round; extra bits may remain |
| `Up` | Toward +∞ |
| `Down` | Toward −∞ |
| `ToZero` | Toward zero |
| `FromZero` | Away from zero |
| `ToEven` | Nearest, ties to even (recommended default) |
| `ToOdd` | Nearest, ties to odd |

Use `RoundingMode::None` for intermediate steps; round once at the end with `set_precision`. That is the practical rule for any multi-step computation.

---

## 5. Error handling

| Variant | Typical result |
| --- | --- |
| `ExponentOverflow(Sign)` | `±Inf` |
| `DivisionByZero` | `NaN` |
| `InvalidArgument` | `NaN` |
| `PrecisionRetryExhausted` | `NaN` — Ziv budget, not a domain error |
| `MemoryAllocation` | `NaN` |

`ExactNum::err()` returns `Option<Error>` on `NaN`. `Error` implements `Display`; with `std` it implements `std::error::Error`. `From<TryReserveError>` and `From<LayoutError>` convert to `MemoryAllocation`.

No operation panics on a numeric domain error. `lu_decomp` / `qr_decomp` / `svd_decomp` return `None` when a workspace `try_reserve_exact` fails. Panics are possible only when the allocator itself panics.

---

## 6. Arithmetic

| Operation | Status | Notes |
| --- | --- | --- |
| `add` / `sub` / `mul` / `div` | ✅ | Explicit `(p, rm)` |
| `add_full_prec` / `sub_full_prec` / `mul_full_prec` | ✅ | No precision reduction; exact when finite values fit |
| `fma` / `mul_add` | ✅ | Fused `a*b+c`; one final round; full-width product when magnitudes overlap |
| `rem` | ✅ | fmod-style; no `p` / `rm` |
| `reciprocal` | ✅ | Schoolbook (< 3 words); Newton (≥ 3 words) |
| `powi` | ✅ | `usize` exponent; binary exponentiation |
| `powsi` | ✅ | `isize` exponent |
| `pow` | ✅ | Arbitrary `ExactNum` exponent; needs `Consts` |
| `two_sum` / `two_product` | ✅ | `(hi, lo)`; `p` / `rm` round `hi` only; `hi + lo` equals the exact result |
| `fused_sum` / `fused_dot` | ✅ | Extra-precision accumulation; one final round |
| `polyval` | ✅ | Horner; `coeffs[0]` is the constant term (lowest degree first) |
| `+=` / `-=` / `*=` / `/=` operator traits | 🚫 | Not provided; use named methods with explicit `p` and `rm` |
| `%` operator trait | 🚫 | Use `rem` method or `expr!` `%` |
| `Hash` / total `Ord` including `NaN` | 🚫 | `NaN` makes these semantically wrong |

**Operator traits** (`Add`, `Sub`, `Mul`, `Div`) exist for all combinations of `ExactNum` and `&ExactNum` at **128 bits, `ToEven`**. Use named methods for any other precision or rounding mode.

---

## 7. Comparison, clamp, sign

| Method | Notes |
| --- | --- |
| `cmp` | `Option<SignedWord>` — `None` if either is NaN |
| `abs_cmp` | Compare absolute values |
| `min` / `max` / `clamp` | |
| `signum` | −1, 0, or +1 as `ExactNum` |
| `abs` | |
| `copysign` | Sign of `sign`, magnitude of `self` |
| `next_after` | Next representable value toward `toward` at precision `p` |

---

## 8. Exponent and mantissa

| Method | Notes |
| --- | --- |
| `exponent` | `Option<Exponent>` |
| `precision` | `Option<usize>` mantissa bit length |
| `sign` | `Option<Sign>` |
| `frexp` | `(significand in [0.5, 1), exponent)` — method only, not an `expr!` leaf |
| `ldexp` / `scalb` | `self × 2^n` |
| `logb` | `floor(log2(\|self\|))` as `ExactNum` |
| `ilogb` | Same as `Option<Exponent>` — method only, not an `expr!` leaf |
| `set_precision` | Round to `p`; `Result<(), Error>` |
| `try_set_precision` | Ziv-style: round only if the rounding is unique given working precision `s` |

---

## 9. Integer / rounding of the value

| Method | Notes |
| --- | --- |
| `int` | Integer part |
| `fract` | Fractional part |
| `ceil` / `floor` | |
| `round` | Round with `n` binary fractional bits and mode `rm` |

---

## 10. Roots, logs, exponentials

All take `(p, rm)`. Those marked `cc` also need a `Consts` cache.

| Method | `expr!` leaf | Notes |
| --- | --- | --- |
| `sqrt` | `sqrt(x)` | |
| `cbrt` | `cbrt(x)` | |
| `nth_root(n, p, rm)` | `root(x, n)` | `n=0` → NaN; `n=2`/`3` delegate to `sqrt`/`cbrt`; even `n` of negative → NaN |
| `ln` / `log2` / `log10` | yes | `cc` |
| `log(base, p, rm, cc)` | `log(x, b)` | `cc` |
| `log1p` | `log1p(x)` | `cc` |
| `exp` / `exp2` / `exp10` | yes | `cc` |
| `expm1` | yes | `cc` |
| `pow` | `pow(b, x)` | `cc` |

---

## 11. Circular and hyperbolic

All take `(p, rm, cc)` except `hypot` (no cache needed).

| Method | `expr!` leaf | Notes |
| --- | --- | --- |
| `sin` / `cos` / `tan` | yes | `rem_pi` reduction internally |
| `sin_cos` | **no** | Paired; shared argument reduction |
| `asin` / `acos` / `atan` | yes | |
| `atan2(y=self, x, …)` | `atan2(y, x)` | |
| `hypot` | `hypot(x, y)` | No cache |
| `rem_pi` | `rem_pi(x)` | Reduces into `(−2π, 2π)` |
| `sinh` / `cosh` / `tanh` | yes | |
| `sinh_cosh` | **no** | Paired |
| `asinh` / `acosh` / `atanh` | yes | Large `\|x\|` adds `2\|e\|` extra bits |

---

## 12. Special functions

| Method | `expr!` leaf | Notes |
| --- | --- | --- |
| `erf` / `erfc` | yes | MPFR 1-ULP on `\|x\| ≲ 4`. Rustdoc `# Precision` on every row in this table |
| `gamma` | yes | Poles at non-positive integers → NaN |
| `ln_gamma` | yes | Positive `self` only |
| `digamma` | yes | Reflection for \(z<0\); poles at non-positive integers → NaN |
| `gammainc` | `gammainc(s, x)` | Lower \(\gamma(s,x)\); \(s>0\), \(x\ge 0\) |
| `gammainc_upper` | `gammainc_upper(s, x)` | Upper \(\Gamma(s,x)=\Gamma(s)-\gamma(s,x)\) |
| `ei` | yes | Cauchy PV for \(x<0\); \(x=0\) is a pole |
| `si` | yes | Odd; series or auxiliary \(f,g\). \(+\infty\to\pi/2\) |
| `ci` | yes | `self > 0`; series or auxiliary \(f,g\). \(+\infty\to 0\) |
| `li` | yes | `self > 0`, `self ≠ 1`; `Ei(ln self)` |
| `fresnel_s` / `fresnel_c` | yes | Odd; series or auxiliary \(f,g\). \(\pm\infty\to\pm 1/2\) |
| `bessel_j(n, p, rm, cc)` | `bessel_j(x, n)` | Integer order; Miller recurrence for large \(n\); \(J_n(0)=\delta_{n0}\) |
| `bessel_j_nu` / `bessel_y` / `bessel_i` / `bessel_k` | `bessel_j_nu(x, ν)` etc. | Real order. \(K\): \(x>0\) |
| `elliptic_k` / `elliptic_e_complete` | `elliptic_k` / `elliptic_e` | Complete; \(m=k^2\); \(K(1)=+\infty\); \(K(m>1)=m^{-1/2}K(1/m)\); \(E\) for \(m\le 1\) |
| `elliptic_f` / `elliptic_e` | `elliptic_f` / `elliptic_e_inc` | Incomplete; \(x=\sin\varphi\), \(\lvert x\rvert\le 1\) |
| `elliptic_pi_complete` / `elliptic_pi` | `elliptic_pi` / `elliptic_pi_inc` | \(n<1\), \(m<1\) complete |
| `legendre_p` / `assoc_legendre_p` | `legendre_p(x, n)` / `legendre_p_assoc(x, n, m)` | Integer \(n\); recurrence at \(p+O(n)\) bits; Condon–Shortley |
| `hypergeom_2f1` | `hypergeom_2f1(a,b,c,z)` | Series / Gauss / Pfaff; real continuation for \(z\le -1\) when defined; non-real \(z>1\) → NaN |
| `betainc` | `betainc(a,b,x)` | Regularized \(I_x(a,b)\); \(a>0\), \(b>0\), \(x\in[0,1]\) |
| `normal_pdf` / `normal_cdf` | `normal_pdf(x,μ,σ)` | \(\varphi\); \(\Phi=(1+\mathrm{erf}(z/\sqrt{2}))/2\); \(\sigma>0\). \(\varphi(0,0,1)=1/\sqrt{2\pi}\); \(\Phi(0,0,1)=1/2\) |
| `gamma_pdf` / `beta_pdf` | yes | Scale \(\beta\); \(B\) via \(\Gamma\). \(\mathrm{gamma\_pdf}(1,1,1)=e^{-1}\) |
| `poisson_pmf` / `binomial_pmf` | yes | Non-negative integer \(k\); \(\mathrm{poisson}(0,1)=e^{-1}\) |
| `chi_squared_cdf` / `student_t_pdf` | yes | Regularized lower gamma; \(t\) via \(\Gamma\). \(\chi^2_2(2\ln 20)=19/20\) |

---

## 12b. Software IEEE and arrays

| Item | Status | Notes |
| --- | --- | --- |
| `Ieee32` / `Ieee64` | ✅ | Integer IEEE-754 binary32/binary64; `from_bits` / `to_bits`; add/mul/div/sqrt/FMA |
| `Ieee32Array` / `Ieee64Array` | ✅ | Row-major; elementwise, `sum`/`dot`, software `matmul`; integer SIMD add/sub/mul/div/sqrt/fma; specials via `ExactNum` |
| `ExactNumArray` | ✅ | Shared `p`; row-major elementwise, software `matmul`, `lu_decomp`, `qr_decomp`, `svd_decomp`, `eigen_decomp`, `fft`/`ifft`; `ExactNum` specials; `(2×3)` `sin` matches scalar; shape mismatch → `None` |
| Integer SIMD (IEEE) | ✅ | `u32`/`u64` lanes; SSE2/NEON; `IEEE_SIMD_LANE_WIDTH=4`; add/sub/mul/div/sqrt/fma bit-identical to scalar kernel; 1000-element `Ieee64Array` gold; not an FPU |
| `lu_decomp` / `qr_decomp` | ✅ | Partial-pivot LU; modified Gram–Schmidt QR; singular or failed workspace reserve → `None`; rank-deficient QR → zero \(R_{kk}\) |
| `svd_decomp` | ✅ | Golub–Reinsch; \((U,\Sigma,V^T)\); \(\sigma\) descending; `SVD_ITER_MAX=64` sweeps/value → `None`; empty/NaN/Inf → `None` |
| `eigen_decomp` | ✅ | Symmetric QR; \((\Lambda,V)\) with \(\lambda\) descending; `EIGEN_ITER_MAX=64`; non-symmetric / empty / non-finite → `None` |
| `fft` / `ifft` | ✅ | Radix-2 Cooley–Tukey; `(1,n)`/`(n,1)` real or `(2,n)` complex; unnormalized DFT; `ifft` divides by `n`; `FFT_MAX_POINTS=4096` |
| `ExactRational` | ✅ | `num/den` reduced; `den>0`; `from_i64` / `new`; add/sub/mul/div; `to_exact_num`; `1/3+1/6=1/2`; `2/4=1/2`; `(-3)/(-4)=3/4`; 256-bit `1/3` |
| `ExactInt` | ✅ | Little-endian `Word` limbs; `from_i64`/`u64`/`i128`/`u128`; add/sub/mul; `div_rem`; `gcd`; `pow`; `20!`; `2^100`; `gcd(48,18)=6` |
| `ExactNumPoly` | ✅ | Dense univariate, low-to-high coeffs; `eval`/`add`/`sub`/`mul`/`div_rem`/`gcd`/`compose`/`derivative`/`integral`; `roots_real` through degree `POLY_COMPANION_CLOSED_DEG=2`; `(x²−1)÷(x−1)=(x+1,0)`; `gcd=x−1`; compose; `∂(x³)=3x²`; `±√2` |
| Chebyshev approx | ✅ | `chebyshev_coeffs`/`chebyshev_eval`/`clenshaw`/`chebyshev_error_bound`; Gauss nodes; `exp` on `[-1,1]` 20 terms `<10^{-15}`; Clenshaw `[1,2,3](1/2)=1/2` |
| Orthogonal polys | ✅ | `hermite_he`/`hermite_h`/`laguerre`/`gen_laguerre`/`chebyshev_t`/`chebyshev_u`/`gegenbauer`; `He_4(0)=3`; `L_3(0)=1`; `T_5(\cos(\pi/5))=-1`; `C_2^{(1)}=4x^2-1` |
| Quadrature | ✅ | `gauss_legendre`/`tanh_sinh`/`gauss_laguerre`/`gauss_hermite`; \(x^2\) on \([-1,1]\) is \(2/3\); 20-point \(x^{38}\) is \(2/39\); \(1/\sqrt{1-x^2}=\pi\); Laguerre \(x^2=2\) |
| Root finding | ✅ | `bisect`/`newton`/`brent`/`illinois`; `bisect(sin,[3,4])=π`; `newton(x²−2)=√2`; Brent fewer iters than bisection; no sign change → `None` |
| ODE solvers | ✅ | `rk4` / `rk45_adaptive` (Dormand–Prince 5(4)) / `euler`; `y'=-y` RK4 1000-step error `<10^{-12}`; RK45 `atol=10^{-12}`; Euler error shrinks when `h` halves |
| Discrete transforms | ✅ | `dct`/`idct` (type II / III); `dst`/`idst`; `fft_real`/`ifft_real`; `idct(dct(x))=x`; constant → DC only; cosine energy at `k` and `N-k`; Parseval |
| Window functions | ✅ | Symmetric Hann / Hamming / Blackman; Kaiser `I_0`; rectangular. `hann(4)=[0,3/4,3/4,0]`; Hamming ends `0.08`; Kaiser `β=0` is ones |
| Modular `ExactInt` | ✅ | `mod_pow`/`mod_inv`/`miller_rabin`/`pollard_rho`; `2^{100} ≡ 976371285 (mod 10^9+7)`; `3^{-1}≡5 (mod 7)`; `2^{31}−1` prime; `8051=83×97` |
| Hash / HMAC | ✅ | `sha256`/`sha512`/`hmac_sha256`/`constant_time_eq`; empty and `abc` FIPS vectors; RFC 4231 HMAC TC1 |
| Binary interchange | ✅ | `to_inline_bytes` / `write_bytes` / `to_bytes` / `from_bytes`; 16-byte BE inline; heap `u32` limbs; array shape; invalid → `Err` |
| CSV | ✅ | `Ieee64Array` / `ExactNumArray` `to_csv` / `from_csv`; cells are binary64 bit integers or `Display@p=`; empty → `NAN` |
| HDF5 | ⬜ | Leftover. No C `libhdf5`. Own contiguous subset only, if ever |
| BLAS / blocked / FFT matmul | ⬜ | Not this crate |

Hardware IEEE arithmetic stays forbidden.

---

## 13. `expr!` macro (real)

`expr!(expression, context)` evaluates a tree of `ExactNum` operations at extra working precision, then rounds **once** to the context precision `p` with rounding mode `rm`.

**Context:** `Context`, or tuples `(p, rm, &mut Consts)` / `(p, rm, &mut Consts, emin, emax)`.

**Working precision:** starts at `p + WORD_BIT_SIZE`. Per-subexpression `errs[]` slots grow when cancellation is detected. `p_wrk = p_rnd + sum(errs)`.

**What `expr!` does not guarantee:**
- Bit-identity with MPFR for composite trees — leaves are 1-ULP oracle tested; composite trees use the cancellation heuristic
- Correct rounding of the entire expression in the IEEE fused sense — only the root is rounded to `p`

**Function leaves (complete list):**

`recip`, `sqrt`, `cbrt`, `root`, `ln`, `log2`, `log10`, `log`, `log1p`, `exp`, `exp2`, `exp10`, `expm1`, `pow`, `rem_pi`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`, `hypot`, `fma`, `mul_add`, `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`, `erf`, `erfc`, `gamma`, `ln_gamma`, `digamma`, `gammainc`, `gammainc_upper`, `ei`, `si`, `ci`, `li`, `fresnel_s`, `fresnel_c`, `bessel_j`, `bessel_j_nu`, `bessel_y`, `bessel_i`, `bessel_k`, `elliptic_k`, `elliptic_e`, `elliptic_e_inc`, `elliptic_f`, `elliptic_pi`, `elliptic_pi_inc`, `legendre_p`, `legendre_p_assoc`, `hypergeom_2f1`, `betainc`, `normal_pdf`, `normal_cdf`, `gamma_pdf`, `beta_pdf`, `poisson_pmf`, `binomial_pmf`, `chi_squared_cdf`, `student_t_pdf`, `ldexp`, `scalb`, `logb`.

**Named constants in the expression:** `pi`, `e`, `ln_2`, `ln_10`, `sqrt2`, `phi`, `euler_gamma`.

**Not in `expr!`:** `frexp`, `ilogb`, `sin_cos`, `sinh_cosh`, `int`/`fract`/`ceil`/`floor`/`round`, `min`/`max`/`clamp`, `cmp`, `copysign`, `next_after`, `powi`/`powsi`, `abs`/`signum`, parse/format, raw parts, `nth_root` by that name (use `root`). Complex values need `cexpr!`.

**Scoped rounding mode:**

```rust
let down = ctx.with_rounding_mode(RoundingMode::Down, |ctx| expr!(1 / 3, ctx));
```

The previous mode is restored when the closure returns.

---

## 14. `cexpr!` macro (complex)

`cexpr!(expression, context)` is the complex analogue of `expr!`. Same extra-precision loop; cancellation tracked independently on real and imaginary parts via two `errs[]` slots (`complex_cancel_bits` returns `(Option<usize>, Option<usize>)`).

**`I`:** expands to `ExactComplex::i(p_wrk)`. Lowercase `i` remains a variable. Do not name a variable `e` — it is Euler's number in `cexpr!` as in `expr!`.

**Operators:** `+`, `−`, `×`, `/`. No `%`.

**Function leaves (complete list):**

`recip`, `sqrt`, `cbrt`, `root`, `ln`, `log2`, `log10`, `log`, `log1p`, `exp`, `exp2`, `exp10`, `expm1`, `pow`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `hypot`, `fma`, `mul_add`, `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`, `abs`, `arg`, `conj`, `ldexp`, `scalb`, `logb`, `erf`, `erfc`, `gamma`, `ln_gamma`, `digamma`, `ei`, `si`, `ci`, `li`, `fresnel_s`, `fresnel_c`, `bessel_j_nu`, `bessel_y`, `bessel_i`, `bessel_k`, `elliptic_k`, `elliptic_e`, `elliptic_e_inc`, `elliptic_f`, `elliptic_pi`, `elliptic_pi_inc`, `hypergeom_2f1`.

**Named constants:** `I`, `pi`, `e`, `ln_2`, `ln_10`, `sqrt2`, `phi`, `euler_gamma` (reals lifted as `x + 0i`).

**Branch cuts (observable, pinned by tests):**

| Family | Cut / range |
| --- | --- |
| `arg` | `atan2(im, re)` in `(−π, π]` |
| `ln`, `log*`, `log1p` | Cut on `(−∞, 0]`; `ln(−1) = iπ` |
| `sqrt`, `cbrt`, `root` | Cut on `(−∞, 0]`; `Re(sqrt) ≥ 0`; `sqrt(−1) = +i` |
| `pow` | `exp(w · ln(z))`; `ln` cut on the base |
| `asin` / `acos` / `atan` / `asinh` / `acosh` / `atanh` | Principal branches |
| `ei` / `si` / `ci` | `Ei` cut on (−∞, 0]; `Si`/`Ci` via `Ei(±iz)`; `Ci(0)` → NaN |
| `li` | `Ei(ln z)`; cut on (−∞, 1]; pole at 1 → NaN |
| `fresnel_s` / `fresnel_c` | via `erf`; entire |
| `bessel_j_nu` / `bessel_i` | entire for integer ν; cut on (−∞, 0] otherwise |
| `bessel_y` / `bessel_k` | cut on (−∞, 0]; \(z=0\) → NaN |
| `elliptic_k` | cut on \([1,+\infty)\); \(K(1)=+\infty\) |
| `hypergeom_2f1` | cut on \([1,+\infty)\) in \(z\) |
| incomplete \(F,E,\Pi\) | Carlson cuts on \(1-x^2\), \(1-mx^2\), \(1-nx^2\) along \((-\infty,0]\) |

**Trig identity:** `sin(x + iy) = sin(x)cosh(y) + i·cos(x)sinh(y)`. The complex argument is never passed through `rem_pi`.

**Not in `cexpr!`:**

| Leaf | Why |
| --- | --- |
| `%` / `rem_pi` | Remainder and π-reduction are real-valued |
| `atan2` | No standard two-complex analogue; use `arg` |

---

## 15. `exact!` and `fbig!`

Compile-time parse of a string literal into an exact `ExactNum`. `fbig!` is an alias. Both take a decimal string literal.

```rust
let x = exact!("1.25");
let y = fbig!("3.14159");
```

---

## 16. Constants cache (`Consts`)

`Consts::new() -> Result<Self, Error>` — allocates progressive caches for π, e, ln 2, ln 10, √2, φ, and Euler–Mascheroni γ. Extend on demand; never recompute from scratch.

| Method | Constant |
| --- | --- |
| `pi(p, rm)` | π |
| `e(p, rm)` | e |
| `ln_2(p, rm)` | ln 2 |
| `ln_10(p, rm)` | ln 10 |
| `sqrt2(p, rm)` | √2 |
| `phi(p, rm)` | φ = (1 + √5) / 2 |
| `euler_gamma(p, rm)` | Euler–Mascheroni γ |
| `cache_info()` | Bit lengths of each cached constant |

**Reuse one `Consts` (or one `Context`) across calls.** Creating a new cache per call recomputes everything from scratch.

**`SharedConsts` (`std`):** mutex-wrapped `Consts` for multi-threaded use. `new()`, `with(|cc: &mut Consts| …)`.

**`CachedFBig`:** wrap a precomputed `ExactNum`; `round(p, rm)` for lazy precision extension.

---

## 17. `Context`

`zenith_float::ctx::Context` bundles precision, rounding mode, `Consts`, `emin`, and `emax`.

| Method | Notes |
| --- | --- |
| `new(p, rm, cc, emin, emax)` | |
| `with_rounding_mode(rm, closure)` | Scoped rounding mode; restores on return |
| `set_precision` / `set_rounding_mode` / `set_consts` / `set_emin` / `set_emax` | |
| `precision` / `rounding_mode` / `consts` / `emin` / `emax` | |
| `const_pi` / `const_e` / `const_ln2` / `const_ln10` / `const_sqrt2` / `const_phi` / `const_euler_gamma` | At context `p` / `rm` |
| `clone` | `Result<Self, Error>` — deep-clones the cache |

**Exponent window:** results outside `[emin, emax]` become `0` or `Inf`. Use the narrowest window that still holds real results — a large unused range can increase internal working precision for no benefit.

---

## 18. `ExactComplex`

Cartesian `re + i·im` as two `ExactNum` values.

| API | Notes |
| --- | --- |
| `new(re, im)` / `zero(p)` / `one(p)` / `i(p)` | |
| `re` / `im` | References |
| `is_nan` | Either part NaN |
| `conj` | |
| `abs(p, rm)` | Modulus (`hypot`) |
| `arg(p, rm, cc)` | Argument (`atan2`) |
| `add` / `sub` / `mul` / `div` | With `(p, rm)` |
| `exp` / `ln` / `sin` / `cos` / `tan` | With `(p, rm, cc)` |
| `sinh` / `cosh` / `tanh` | |
| `sqrt` / `pow` | Principal branch |
| `asin` / `acos` / `atan` / `asinh` / `acosh` / `atanh` | Principal branches |
| `log2` / `log10` / `log` / `log1p` | Via `ln` ratios |
| `exp2` / `exp10` / `expm1` | Via `exp` |
| `ldexp` / `scalb` | Scale both parts by `2^n` |
| `cbrt` / `nth_root` | Principal; `exp(ln / n)` |
| `hypot` | Principal `sqrt(z² + w²)` |
| `fma` / `mul_add` | Extra-precision `a*b+c` |
| `logb` | `logb(\|z\|)` as a real |
| `Add` `Sub` `Mul` `Div` | 128-bit `ToEven` like reals |
| `erf` / `erfc` | ✅ Faddeeva \(w(z)\); entire |
| `gamma` / `ln_gamma` / `digamma` | ✅ Stirling + reflection; \(\ln\Gamma\) cut on \((-\infty,0]\); poles → NaN |
| `ei` / `si` / `ci` / `li` | ✅ series or asymptotic `Ei`; `Si`/`Ci` via `Ei(±iz)`; `li=Ei(ln z)` |
| `fresnel_s` / `fresnel_c` | ✅ via `erf`; entire |
| `bessel_j_nu` / `bessel_y` / `bessel_i` / `bessel_k` | ✅ series or Hankel; \(I_ν=i^{-ν}J_ν(iz)\); \(K_ν=(\pi/2)i^{ν+1}H_ν^{(1)}(iz)\) |
| `elliptic_k` / `elliptic_e_complete` / `elliptic_f` / `elliptic_e` / `elliptic_pi_*` | ✅ Carlson in \(\mathbb{C}\); \(K(0)=E(0)=\pi/2\); \(K(1)=+\infty\); Legendre and cut golds |
| `hypergeom_2f1` | ✅ Series / Euler / Pfaff / Kummer; \(2\ln 2\), \(2K/\pi\), Euler identity, \(c=0\) NaN, cut golds |

No `expr!` for complex — use `cexpr!`. Serde: struct `{ "re", "im" }` of decimal strings.

---

## 19. Parse, format, radix conversion

| API | Notes |
| --- | --- |
| `parse(s, rdx, p, rm, cc)` | Returns `ExactNum` (not `Result`); failure → NaN; `err()` recovers the `Error` |
| `parse_exact` / `format_exact` | `ExactRational::parse_exact("0.1")` is `1/10`; `ExactNum::parse_exact` is dyadic-only; `format_exact(Dec, 1/8)="0.125"` |
| `format(rdx, rm, cc)` | `Result<String, Error>` |
| `convert_from_radix` / `convert_to_radix` | Low-level digit-byte API |
| `FromStr` (`std`) | Decimal, `ToEven`, unbounded precision |
| `Display` (`std`) | Decimal |
| `LowerExp` / `UpperExp` (`std`) | Decimal scientific |
| `Binary` / `Octal` / `UpperHex` / `LowerHex` (`std`) | |

**`parse` returns `ExactNum`, not `Result`.** A failed parse is NaN. **`format` returns `Result<String, Error>`.** That asymmetry is intentional and easy to miss.

Radix 2–36. For bases > 10 the exponent uses `_e` so `e` can be a digit.

---

## 20. Correct-rounding infrastructure

| Item | Notes |
| --- | --- |
| `ziv_round(p, rm, compute)` | Call `compute(p_wrk)`, then `try_set_precision` until uniquely rounded or `MAX_PREC_RETRY` exhausted (→ NaN / `PrecisionRetryExhausted`) |
| `ziv_round_vec(p, rm, inputs, compute)` | Same loop; inputs lifted to `p_wrk`; `hypot(3,4)=5`; `atan2(1,1)=π/4` at 256 bits |
| `Ball { mid, rad }` | Certified `add` / `mul` / `exp` / `sin` / `cos` / `ln` / `sqrt` / `erf` / `bessel_j0` / `bessel_j1`; Lipschitz + `BALL_TRANSCENDENTAL_ERROR_TERMS` ulps; `contains(x, p)` |
| `ComplexBall { mid, rad }` | Disk in \(\mathbb{C}\); `add`/`mul`/`exp`/`ln`/`sin`/`cos`; unit-disk `exp` and \(\sin^2+\cos^2=1\) golds |
| `MAX_PREC_RETRY = 256` | Extra word-sized budget per operation; caps at `256 × WORD_BIT_SIZE` bits above `p` |

---

## 21. Serde (`serde` feature)

- **`ExactNum`:** decimal string `"<Display>@p=<bits>"`; JSON integers use `DEFAULT_P`. Rehydration parses at the stored bit count so limbs match across targets
- **`ExactComplex`:** struct with fields `re` and `im` (each encoded as above)
- **`ExactRational`:** `{"num": "...@p=", "den": "...@p="}`
- **`ExactInt`:** signed decimal string
- **`ExactNumArray`:** `{"shape", "data", "p", "rm"}`; shape product must equal `data.len()`
- **`Ieee32Array` / `Ieee64Array`:** `{"shape", "data"}` where `data` is integer IEEE bit patterns (not hardware `f32`/`f64` JSON numbers)
- **`Ball`:** `{"mid", "rad"}` encoded strings

---

## 22. Random (`random` feature)

| Item | Notes |
| --- | --- |
| `DEFAULT_RANDOM_SEED` | `0x5EED_CAFE_BADC_0D00` |
| `random_seed()` / `reseed_random(seed)` | |
| `seeded_random::<T>()` | Draw from the crate RNG |
| `ExactNum::random_normal` | Random finite with exponent in `[exp_from, exp_to]` |
| `ZENITH_TEST_SEED` | Environment variable to replay test failures |

Unit tests default to the deterministic seed. Without reseeding, production code uses OS entropy.

---

## 23. What is not in this crate (permanent)

These are design decisions, not a backlog:

- Hardware floating-point arithmetic, `libm`, or compiler float codegen
- Rust hardware IEEE type tokens in kernel sources; convert at the boundary with `Ieee32::from_bits` / `Ieee64::from_bits` (or pack from `frexp` / `ilogb` on `ExactNum`)
- `+=` / `-=` / `*=` / `/=` assigning operators
- `Hash` or total `Ord` that includes NaN
- `%` operator trait — use `rem` method or `expr!` `%`
- Symbolic CAS, formula rewriting, or expression IR — this crate is numeric only
- MPFR as a runtime engine — it is an oracle in `mpfr-tests` only

---

## 24. Comparison vs astro-float and dashu-float

| Capability | zenith-float | astro-float 0.9.6 | dashu-float 0.6.0 |
| --- | :---: | :---: | :---: |
| Pure Rust kernel (no MPFR in lib) | ✅ | ✅ | ✅ |
| `no_std` + allocator | ✅ | ✅ | ✅ |
| Hardware IEEE arithmetic | 🚫 forbidden | converters | API-edge literals |
| Software `Ieee32` / `Ieee64` + arrays / matmul | ✅ | ⬜ | ⬜ |
| `expr!` / `cexpr!` + `Context` | ✅ both | ✅ real only | 🟡 macros elsewhere |
| `hypot`, `atan2`, `log1p`, `expm1` | ✅ | ⬜ | ✅ |
| `exp2`, `exp10`, `rem_pi` | ✅ | ⬜ | 🟡 partial |
| Native `+ − × ÷` operators | ✅ | ⬜ | ✅ |
| `fma` / `mul_add` | ✅ | ⬜ | 🟡 |
| `two_sum` / `two_product` / `fused_sum` / `fused_dot` / `polyval` | ✅ | ⬜ | ⬜ |
| `fbig!` / compile-time float literals | ✅ | ⬜ | ✅ |
| Parse/format bases 2–36 | ✅ | ⬜ (bin/oct/dec/hex) | ✅ |
| Complex (`ExactComplex` + `cexpr!`) | ✅ elementary | ⬜ | ✅ `CBig` |
| Complex specials (`erf`, `Γ`, `ψ`, `Ei`, Bessel, elliptic, `_2F1`) | ✅ | ⬜ | ⬜ |
| General `nth_root(n)` | ✅ | ✅ | ✅ |
| `sin_cos` / `sinh_cosh` paired | ✅ | ✅ | ✅ |
| Special functions (`erf`, `Γ`, `J_n`) | ✅ | ⬜ | 🟡 / separate |
| `euler_gamma` constant | ✅ | ⬜ | <!-- verify --> |
| `with_rounding_mode` scoped closure | ✅ | ⬜ | <!-- verify --> |
| Ziv + `Ball` correct-rounding proof | ✅ | ⬜ | ✅ |
| Stack-inlined small values | ✅ `INLINE_WORDS = 2` | ⬜ | ✅ |
| MPFR bit-oracle tests (optional) | ✅ | ✅ | fuzz + unit |

<!-- Verify the dashu-float rows marked with comments before publishing. The astro-float rows were current at 0.9.6; re-check at publish time. -->

---

## 25. Leftovers (this crate — walk the build plan)

Not a second product. First open implementation slice is **§19.3 Benchmark suite**.

| Plan | Item |
| --- | --- |
| §19.3 | bench specials / matmul / FFT / LU |
| §17.3 leftover | HDF5: own contiguous subset, not `libhdf5` / not a general crate |
| §19.4 | `scripts/zenith_prepublish.sh` |
| §20.2 | hex limb CI (`arm` / `wasm` / `32bit`) |

---

## 26. Version history of this document

| Version | Date | Changes |
| --- | --- | --- |
| 0.1.0 | 2026-08-31 | SIMD IEEE div/sqrt/fma; thumb `no_std`; `expr!`/`cexpr!` composite golds. Prior: complex specials through `_2F1`; arrays; `Ball`/`ComplexBall`; LU/QR/SVD; FFT; Precision rustdoc; REPRO; `ExactRational`; `ExactInt`; `parse_exact`/`format_exact`; `ziv_round_vec`; distribution kernels; RNG; `ExactNumPoly`; Chebyshev; orthogonal polynomials; quadrature; root finding; ODE solvers; DCT/DST/`fft_real`; windows; modular `ExactInt`; SHA-2 / HMAC; serde `@p=` + IEEE bits; 16-byte BE binary interchange; CSV. TODO file retired; walk `ZENITH_FLOAT_BUILD_PLAN.md` |

---

## 27. Related documents

| Document | Location | Audience |
| --- | --- | --- |
| Build plan (walk list + status) | `doc/ZENITH_FLOAT_BUILD_PLAN.md` | Maintainers; next slice |
| Build checklist | `doc/BUILD_CHECKLIST.md` | Maintainers; CI and crate inventory |
| Getting started | `doc/GETTING_STARTED.md` | New users; first working example |
| Help | `doc/HELP.md` | Why the API looks this way; recipes; common mistakes |
| Full method inventory | `doc/LIBRARY.md` | Every public type and method |
| `expr!` rounding contract | `doc/EXPR.md` | Per-op working precision; what `expr!` guarantees and does not |
| Precision and retry budget | `doc/PRECISION.md` | `MAX_PREC_RETRY`; exponent scaling; `expr!` bounds |
| Reproducibility and citation | `doc/REPRODUCIBILITY.md` | What determines a result; how to replay; citation line |
| Error-bound theory | `doc/README.md` | ULP / series error bounds used by `expr!`; contributors |
