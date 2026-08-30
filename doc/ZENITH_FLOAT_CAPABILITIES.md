# zenith-float Capability Reference

**Version:** <!-- e.g. 0.1.0 -->  
**Date:** <!-- release date -->  
**License:** MIT OR Apache-2.0  
**Status:** Public crate  

This document states exactly what `zenith-float` provides, what it does not provide, and what the named constants and caps mean in practice. It is not a tutorial. For a working example, see `doc/GETTING_STARTED.md`. For full method signatures, see `doc/LIBRARY.md`. For `expr!` rounding contracts, see `doc/EXPR.md`.

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
| `serde` | no | `Serialize` / `Deserialize` for `ExactNum` and `ExactComplex`; implies `std` |
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

`ExactNum::err()` returns `Option<Error>` on `NaN`. `Error` implements `Display`; with `std` it implements `std::error::Error`. `From<TryReserveError>` converts to `MemoryAllocation`.

No operation panics on a numeric domain error. Panics are possible only on allocation failure when the allocator itself panics.

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
| `erf` / `erfc` | yes | MPFR 1-ULP on `\|x\| ≲ 4` |
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

---

## 12b. Software IEEE and arrays

| Item | Status | Notes |
| --- | --- | --- |
| `Ieee32` / `Ieee64` | ✅ | Integer IEEE-754 binary32/binary64; `from_bits` / `to_bits`; add/mul/div/sqrt/FMA |
| `Ieee32Array` / `Ieee64Array` | ✅ | Row-major; elementwise, `sum`/`dot`, software `matmul`; integer SIMD add/mul; specials via `ExactNum` |
| `ExactNumArray` | ✅ | Shared `p`; row-major elementwise, software `matmul`, `ExactNum` specials |
| Integer SIMD (IEEE add/mul) | ✅ | `u32`/`u64` lanes; SSE2/NEON; bit-identical to scalar kernel; not an FPU |
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

`recip`, `sqrt`, `cbrt`, `root`, `ln`, `log2`, `log10`, `log`, `log1p`, `exp`, `exp2`, `exp10`, `expm1`, `pow`, `rem_pi`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`, `hypot`, `fma`, `mul_add`, `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`, `erf`, `erfc`, `gamma`, `ln_gamma`, `digamma`, `gammainc`, `gammainc_upper`, `ei`, `si`, `ci`, `li`, `fresnel_s`, `fresnel_c`, `bessel_j`, `bessel_j_nu`, `bessel_y`, `bessel_i`, `bessel_k`, `elliptic_k`, `elliptic_e`, `elliptic_e_inc`, `elliptic_f`, `elliptic_pi`, `elliptic_pi_inc`, `legendre_p`, `legendre_p_assoc`, `hypergeom_2f1`, `betainc`, `ldexp`, `scalb`, `logb`.

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

`recip`, `sqrt`, `cbrt`, `root`, `ln`, `log2`, `log10`, `log`, `log1p`, `exp`, `exp2`, `exp10`, `expm1`, `pow`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `hypot`, `fma`, `mul_add`, `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`, `abs`, `arg`, `conj`, `ldexp`, `scalb`, `logb`.

**Named constants:** `I`, `pi`, `e`, `ln_2`, `ln_10`, `sqrt2`, `phi`, `euler_gamma` (reals lifted as `x + 0i`).

**Branch cuts (observable, pinned by tests):**

| Family | Cut / range |
| --- | --- |
| `arg` | `atan2(im, re)` in `(−π, π]` |
| `ln`, `log*`, `log1p` | Cut on `(−∞, 0]`; `ln(−1) = iπ` |
| `sqrt`, `cbrt`, `root` | Cut on `(−∞, 0]`; `Re(sqrt) ≥ 0`; `sqrt(−1) = +i` |
| `pow` | `exp(w · ln(z))`; `ln` cut on the base |
| `asin` / `acos` / `atan` / `asinh` / `acosh` / `atanh` | Principal branches |

**Trig identity:** `sin(x + iy) = sin(x)cosh(y) + i·cos(x)sinh(y)`. The complex argument is never passed through `rem_pi`.

**Not in `cexpr!`:**

| Leaf | Why |
| --- | --- |
| `%` / `rem_pi` | Remainder and π-reduction are real-valued |
| `atan2` | No standard two-complex analogue; use `arg` |
| `erf` / `erfc` / `gamma` / `ln_gamma` / `ei` / `si` / `ci` / `li` / `fresnel_*` / `bessel_j` | No complex kernel in this crate |

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

No `expr!` for complex — use `cexpr!`. Serde: struct `{ "re", "im" }` of decimal strings.

---

## 19. Parse, format, radix conversion

| API | Notes |
| --- | --- |
| `parse(s, rdx, p, rm, cc)` | Returns `ExactNum` (not `Result`); failure → NaN; `err()` recovers the `Error` |
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
| `Ball { mid, rad }` | First-order interval arithmetic; `add` / `mul` with rounding ulp in radius; `contains(x, p)` |
| `MAX_PREC_RETRY = 256` | Extra word-sized budget per operation; caps at `256 × WORD_BIT_SIZE` bits above `p` |

---

## 21. Serde (`serde` feature)

- **`ExactNum`:** serializes as a decimal string (`Display`); deserializes from that string or JSON integers (`i64` / `u64` / `i128` / `u128`)
- **`ExactComplex`:** struct with fields `re` and `im`

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
- Complex `erf`, `erfc`, `gamma`, `ln_gamma`, `ei`, `si`, `ci`, `li`, `fresnel_s`, `fresnel_c`, `bessel_j` — no complex kernel
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
| Complex (`ExactComplex` + `cexpr!`) | ✅ | ⬜ | ✅ `CBig` |
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

## 25. Version history of this document

| Version | Date | Changes |
| --- | --- | --- |
| 0.1.0 | <!-- date --> | Initial release |

---

## 26. Related documents

| Document | Location | Audience |
| --- | --- | --- |
| Getting started | `doc/GETTING_STARTED.md` | New users; first working example |
| Help | `doc/HELP.md` | Why the API looks this way; recipes; common mistakes |
| Full method inventory | `doc/LIBRARY.md` | Every public type and method |
| `expr!` rounding contract | `doc/EXPR.md` | Per-op working precision; what `expr!` guarantees and does not |
| Precision and retry budget | `doc/PRECISION.md` | `MAX_PREC_RETRY`; exponent scaling; `expr!` bounds |
| Error-bound theory | `doc/README.md` | ULP / series error bounds used by `expr!`; contributors |
