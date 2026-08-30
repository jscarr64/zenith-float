# zenith-float library inventory

This is a complete inventory of what the **`zenith-float` crate** exposes as a software math library. Arithmetic uses integer limbs only. Hardware IEEE binary interchange formats are not part of this crate; convert to those formats in the caller if needed.

Related docs: [getting started](GETTING_STARTED.md), [help](HELP.md), [error bounds](README.md), [`expr!` rounding](EXPR.md).

License: MIT OR Apache-2.0.

---

## 1. Crate layout

| Crate | Role |
| --- | --- |
| **`zenith-float`** | Public package. Re-exports the numeric kernel and the macros `expr!`, `cexpr!`, `exact!`, `fbig!`. |
| **`zenith-float-num`** | Numeric kernel. Applications should depend on `zenith-float`, not this crate. |
| **`zenith-float-macro`** | Procedural macros. Not a direct application dependency. |
| **`zenith-float-compare`** | Workspace comparison benches against other float crates. Not part of the public math API. |

`no_std` is supported when a memory allocator is available (`default-features = false` to drop `std`).

---

## 2. Cargo features (`zenith-float`)

| Feature | Default | What it enables |
| --- | --- | --- |
| `std` | yes | `Display` / radix format traits, `FromStr`, `std::error::Error` for `Error`, `SharedConsts`, serde when `serde` is on. |
| `random` | no | `ExactNum::random_normal`, `seeded_random`, `random_seed`, `reseed_random`, `DEFAULT_RANDOM_SEED`. |
| `serde` | no | `Serialize` / `Deserialize` for `ExactNum` and `ExactComplex`. Implies `std`. |
| `mpfr-tests` | no | Optional MPFR comparison tests in the kernel crate (Linux x86_64, `rug`). Not a runtime math engine. |

---

## 3. Numeric model

A finite `ExactNum` is a **sign**, a **binary mantissa** stored as an array of words, and an **`Exponent` (`i32`)**. Precision is a bit count, always rounded **up to the word size**.

| Item | Value |
| --- | --- |
| `Word` | `u64` except on 32-bit targets, where it is `u32` |
| `WORD_BIT_SIZE` | 64 or 32 accordingly |
| `WORD_MAX` | `Word::MAX` |
| `WORD_BASE` | `WORD_MAX + 1` as a double-width integer |
| `WORD_SIGNIFICANT_BIT` | high bit of a word |
| `INLINE_WORDS` | `2` — small mantissas stored inline |
| `EXPONENT_BIT_SIZE` | 32 |
| `EXPONENT_MAX` | `i32::MAX` on non-32-bit pointers; `i32::MAX / 4` on 32-bit |
| `EXPONENT_MIN` | `i32::MIN` on non-32-bit pointers; `i32::MIN / 4` on 32-bit |
| Default operator precision | 128 bits (`ToEven`) for `+ − × ÷` on `ExactNum` (internal default; not a public constant) |
| `MAX_PREC_RETRY` | `256` — extra-precision budget for correct-rounding retries (`try_set_precision`, `ziv_round`, transcendentals) |

Special values: **`+Inf`**, **`-Inf`**, **`NaN`** (optional associated `Error`). Finite values may be **subnormal** at `EXPONENT_MIN`. Results may be marked **inexact** when bits were rounded or an argument was inexact.

Public sentinels: `NAN`, `INF_POS`, `INF_NEG`.

---

## 4. Rounding

`RoundingMode` (copy, eq, debug):

| Variant | Meaning |
| --- | --- |
| `None` | Skip the final round; extra bits may remain. |
| `Up` | Round half toward +∞ |
| `Down` | Round half toward −∞ |
| `ToZero` | Round half toward zero |
| `FromZero` | Round half away from zero |
| `ToEven` | Round half to even |
| `ToOdd` | Round half to odd |

Methods that take a mode other than `None` round to the requested precision. `expr!` raises working precision for cancellation; it does not itself guarantee correct rounding. See [EXPR.md](EXPR.md).

---

## 5. Errors

`Error`:

| Variant | Typical public mapping |
| --- | --- |
| `ExponentOverflow(Sign)` | `±Inf` |
| `DivisionByZero` | `NaN` |
| `InvalidArgument` | `NaN` |
| `PrecisionRetryExhausted` | `NaN` — Ziv/`MAX_PREC_RETRY` budget, not a domain error |
| `MemoryAllocation` | `NaN` |

`ExactNum::err()` returns the associated error on `NaN`. `Error` implements `Display`; with `std` it implements `std::error::Error`. `From<TryReserveError>` → `MemoryAllocation`.

---

## 6. Sign and radix

**`Sign`:** `Pos`, `Neg`. Methods: `invert`, `is_positive`, `is_negative`, `to_int` (`i8`: +1 / −1). `Hash`, `Eq`, `Copy`.

**`Radix`:** bases **2 through 36**. Named constants `Bin` (2), `Oct` (8), `Dec` (10), `Hex` (16). Methods: `try_new(u8)`, `value`, `commensurable_shift` (power-of-two bases), `uses_underscore_exponent` (bases > 10 use `_e` so `e` can be a digit), `bits_per_digit`. `From<Radix> for u8`.

---

## 7. Public types and re-exports

From `zenith_float` / `zenith_float_num`:

- `ExactNum`, `ExactComplex`
- `Consts`, `ConstCache` (alias of `Consts`), `ConstCacheInfo`, `CachedFBig`, `SharedConsts` (`std` only)
- `Context` (module `zenith_float::ctx`), trait `Contextable`
- `Ball`, `ziv_round`
- `RadixFloat`
- `FromExt`
- `RoundingMode`, `Radix`, `Sign`, `Error`, `Exponent`, `Word`
- Word/exponent constants listed in §3
- `MAX_PREC_RETRY`, `INLINE_WORDS`
- `NAN`, `INF_POS`, `INF_NEG`
- Feature `random`: `random_seed`, `reseed_random`, `seeded_random`, `DEFAULT_RANDOM_SEED`

Module `ctx` is public. `macro_util` is `#[doc(hidden)]` and exists for `expr!` / `cexpr!` expansion (`check_exponent_range`, `check_complex_exponent_range`, `complex_cancel_bits`, `compute_added_err`, `ErrAlgo`, `TrigFun`, …). Do not treat it as application API.

Not re-exported: internal `Mantissa`, `WordBuf`, series helpers, `DEFAULT_P`.

---

## 8. `ExactNum` — constructors and classification

| Method / item | Signature / notes |
| --- | --- |
| `new` | `new(p: usize)` — zero at precision `p` |
| `nan` | `nan(err: Option<Error>)` |
| `from_word` | `from_word(d: Word, p)` |
| `from_i8`, `from_i16`, `from_i32`, `from_i64`, `from_i128` | signed integers + precision |
| `from_u8`, `from_u16`, `from_u32`, `from_u64`, `from_u128` | unsigned integers + precision |
| `From` | `i8`/`u8`/`i16`/`u16`/`i32`/`u32`/`i64`/`u64`/`i128`/`u128` at default 128-bit precision |
| `Default` | `new` at default 128-bit precision |
| `from_raw_parts` | mantissa slice, used length, sign, exponent, inexact flag |
| `from_words` | mantissa slice, sign, exponent |
| `max_value` / `min_value` | largest / most-negative finite at precision `p` |
| `min_positive` / `min_positive_normal` | smallest positive (incl. subnormal) / smallest normal |
| `random_normal` | feature `random`: random finite with exponent in `[exp_from, exp_to]` |
| `is_inf_pos` / `is_inf_neg` / `is_inf` | infinities |
| `is_nan` | |
| `is_int` | integer-valued finite |
| `is_zero` / `is_positive` / `is_negative` / `is_subnormal` | |
| `err` | `Option<Error>` on NaN |
| `classify` | `core::num::FpCategory` (`Zero`, `Subnormal`, `Normal`, `Infinite`, `Nan`) |
| `inexact` / `set_inexact` | exactness flag |
| `Clone`, `Debug` | |
| `PartialEq`, `Eq` | equality via `cmp`; NaN is not equal |
| `PartialOrd` | `None` if either is NaN |
| `Neg` | `ExactNum` and `&ExactNum` |

No `Hash`, no total `Ord`, no `AddAssign` / `MulAssign`.

---

## 9. `ExactNum` — arithmetic (explicit precision)

Unless noted, results round with `(p, rm)`.

| Method | Meaning |
| --- | --- |
| `add` / `sub` / `mul` / `div` | four operations |
| `add_full_prec` / `sub_full_prec` / `mul_full_prec` | no precision reduction (exact product/sum of finite values when it fits) |
| `fma` / `mul_add` | fused `a*b+c` with one final round; `mul_add` is an alias of `fma` |
| `two_sum` / `two_product` | `(hi, lo)`; `p`/`rm` round **hi** only; finite `hi + lo` equals the exact limb sum/product |
| `fused_sum` / `fused_dot` | extra-precision accumulation, one final round |
| `polyval` | Horner; `coeffs[0]` is the constant term (lowest degree first) |
| `rem` | remainder (`self` rem `d2`); no `p`/`rm` |
| `reciprocal` | `1/self` |
| `neg` | copy with inverted sign (also `Neg` / `inv_sign` in place) |
| `pow` | `self^n` (arbitrary `ExactNum` exponent); needs `Consts` |
| `powi` | `self^n` for `usize` `n` |
| `powsi` | `self^n` for `isize` `n` |

Operator traits **`Add` `Sub` `Mul` `Div`** for all combinations of `ExactNum` and `&ExactNum`, always at **128 bits, `ToEven`**. Use the named methods for other precision/rounding.

There is **no** `%` operator trait; remainder is `rem` or `expr!` `%`.

---

## 10. `ExactNum` — comparison, clamp, sign

| Method | Meaning |
| --- | --- |
| `cmp` | `Option<SignedWord>` — negative / zero / positive; `None` if unordered (NaN). `SignedWord` is `i128` on 64-bit targets and `i64` on 32-bit; it is the return type of `cmp`/`abs_cmp` even though it is not re-exported at the crate root. |
| `abs_cmp` | compare absolute values |
| `min` / `max` | |
| `clamp` | |
| `signum` | −1, 0, or +1 as `ExactNum` |
| `abs` | absolute value |
| `copysign` | magnitude of `self`, sign of `sign`; sign applied before precision reduction |
| `next_after` | next representable toward `toward` at precision `p` |

---

## 11. `ExactNum` — exponent, mantissa, precision

| Method | Meaning |
| --- | --- |
| `exponent` | `Option<Exponent>` |
| `precision` | `Option<usize>` mantissa bit length |
| `sign` | `Option<Sign>` |
| `as_raw_parts` | `Option<(&[Word], usize, Sign, Exponent, bool)>` |
| `mantissa_digits` | `Option<&[Word]>` |
| `mantissa_max_bit_len` | |
| `is_inline` | mantissa fits in `INLINE_WORDS` |
| `set_exponent` | |
| `set_sign` | |
| `set_precision` | round to `p`; `Result<(), Error>` |
| `try_set_precision` | Ziv-style: round if uniquely determined given working precision `s`; `bool` |
| `frexp` | `(significand in [0.5, 1), exponent)` as `(ExactNum, Exponent)` |
| `ldexp` / `scalb` | `self × 2^n` (`n: Exponent`) |
| `logb` | `floor(log2(|self|))` as a float |
| `ilogb` | same as integer `Option<Exponent>` |

`frexp` and `ilogb` are **methods only** (not `expr!` leaves).

---

## 12. `ExactNum` — integer / rounding of the value

| Method | Meaning |
| --- | --- |
| `int` | integer part |
| `fract` | fractional part |
| `ceil` / `floor` | |
| `round` | round with `n` binary fractional bits and mode `rm` |

---

## 13. `ExactNum` — roots, logs, exponentials

Need `Consts` where marked `cc`.

| Method | `expr!` leaf? |
| --- | --- |
| `sqrt(p, rm)` | `sqrt(x)` |
| `cbrt(p, rm)` | `cbrt(x)` |
| `nth_root(n, p, rm)` | `root(x, n)` — `n=0` → NaN; `n=2`/`3` delegate to sqrt/cbrt; even `n` of negative → NaN |
| `ln` / `log2` / `log10` | yes |
| `log(self, base, p, rm, cc)` | `log(x, b)` |
| `log1p` | `log1p(x)` |
| `exp` / `exp2` / `exp10` | yes |
| `expm1` | yes |
| `pow` | `pow(b, x)` |

---

## 14. `ExactNum` — circular and hyperbolic

All take `(p, rm, cc)` except `hypot` (no cache).

| Method | `expr!` |
| --- | --- |
| `sin` / `cos` / `tan` | yes |
| `sin_cos` → `(sin, cos)` | **no** (shared argument reduction) |
| `asin` / `acos` / `atan` | yes |
| `atan2(y=self, x, …)` | `atan2(y, x)` |
| `hypot` | `hypot(x, y)` |
| `rem_pi` | reduce modulo `2π` into `(-2π, 2π)` |
| `sinh` / `cosh` / `tanh` | yes |
| `sinh_cosh` → `(sinh, cosh)` | **no** |
| `asinh` / `acosh` / `atanh` | yes |

---

## 15. `ExactNum` — special functions

| Method | Notes | `expr!` |
| --- | --- | --- |
| `erf` / `erfc` | `erfc = 1 − erf` | yes |
| `gamma` | poles at non-positive integers → NaN; +Inf at 0 | yes |
| `ln_gamma` | log-gamma for positive `self` | yes |
| `digamma` | \(\psi(z)\); reflection for \(z<0\); poles at non-positive integers | yes |
| `gammainc(s, x)` | lower \(\gamma(s,x)\); \(s>0\), \(x\ge 0\) | `gammainc(s, x)` |
| `gammainc_upper(s, x)` | upper \(\Gamma(s,x)\) | `gammainc_upper(s, x)` |
| `ei` | Cauchy PV for \(x<0\); \(x=0\) is a pole | yes |
| `si` | odd; all real | yes |
| `ci` | `self > 0`; otherwise NaN | yes |
| `li` | `self > 0`, `self ≠ 1`; `Ei(ln self)` | yes |
| `fresnel_s` / `fresnel_c` | odd; series or auxiliary \(f,g\) | yes |
| `bessel_j(n, p, rm, cc)` | `J_n(self)`, integer order; Miller for large \(n\) | `bessel_j(x, n)` |
| `bessel_j_nu` / `bessel_y` / `bessel_i` / `bessel_k` | Real order; \(K\): \(x>0\) | yes |
| `elliptic_k` / `elliptic_e_complete` / `elliptic_f` / `elliptic_e` / `elliptic_pi_complete` / `elliptic_pi` | Carlson; \(m=k^2\), \(x=\sin\varphi\); \(K(1)=+\infty\); \(K(m>1)=m^{-1/2}K(1/m)\) | `elliptic_k`, `elliptic_e`, `elliptic_f`, `elliptic_e_inc`, `elliptic_pi`, `elliptic_pi_inc` |
| `legendre_p` / `assoc_legendre_p` | Integer \(n\); \(p+O(n)\) recurrence; Condon–Shortley | `legendre_p(x, n)`, `legendre_p_assoc(x, n, m)` |
| `hypergeom_2f1` | Series / Gauss / Pfaff; real \(z\le -1\) when defined; non-real \(z>1\) → NaN | `hypergeom_2f1(a,b,c,z)` |
| `betainc` | Regularized \(I_x(a,b)\) | `betainc(a,b,x)` |

---

## 16. `ExactNum` — parse, format, radix conversion

| API | Notes |
| --- | --- |
| `parse(s, rdx, p, rm, cc)` | returns `ExactNum` (not `Result`); Inf / NaN / `err()` on failure |
| `format(rdx, rm, cc)` | `Result<String, Error>`; Inf / −Inf / NaN / Err strings |
| `with_radix(self, radix)` | wrap as `RadixFloat` |
| `convert_from_radix(sign, digits, e, rdx, p, rm, cc)` | digit bytes, scientific exponent in that radix |
| `convert_to_radix(rdx, rm, cc)` | `(Sign, Vec<u8>, Exponent)` |
| `FromStr` (`std`) | decimal parse, `ToEven`, unbounded precision then value |
| `Display` (`std`) | decimal |
| `LowerExp` / `UpperExp` (`std`) | decimal scientific (`e` / `E`) |
| `Binary` / `Octal` (`std`) | |
| `UpperHex` / `LowerHex` (`std`) | hex; specials stay `Inf` / `NaN` |

For radices where `e` is a digit, the exponent is written with `_e`.

---

## 17. Macros

Public macros (crate root): `expr!`, `cexpr!`, `exact!`, `fbig!`. Import them like any other item, e.g. `use zenith_float::{expr, cexpr, exact};`, or call them as `zenith_float::exact!("1.25")`.

### `expr!(expression, context)`

**Context** (`Contextable`): `Context`, or tuples

- `(precision, RoundingMode, &mut Consts)` — exponent limits `EXPONENT_MIN`/`MAX`
- `(precision, RoundingMode, &mut Consts, emin, emax)` — `emin` clamped to `[EXPONENT_MIN, 0]`, `emax` to `[0, EXPONENT_MAX]`

**Syntax allowed:** paths (variables), integer / float / string literals, unary `-`, binary `+ − * / %`, calls, parentheses. Inputs are treated as exact (`set_inexact(false)`). Working precision is raised; final result uses the context rounding mode. Overflow of `emin`/`emax` → `±Inf`.

**Binary operators:** `+`, `-`, `*`, `/`, `%` (maps to `ExactNum::rem`).

**Function leaves (complete list):**

`recip`, `sqrt`, `cbrt`, `root`, `ln`, `log2`, `log10`, `log`, `log1p`, `exp`, `exp2`, `exp10`, `expm1`, `pow`, `rem_pi`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`, `hypot`, `fma`, `mul_add`, `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`, `erf`, `erfc`, `gamma`, `ln_gamma`, `digamma`, `gammainc`, `ei`, `si`, `ci`, `li`, `fresnel_s`, `fresnel_c`, `bessel_j`, `bessel_j_nu`, `bessel_y`, `bessel_i`, `bessel_k`, `elliptic_k`, `elliptic_e`, `elliptic_e_inc`, `elliptic_f`, `elliptic_pi`, `elliptic_pi_inc`, `legendre_p`, `legendre_p_assoc`, `hypergeom_2f1`, `betainc`, `ldexp`, `scalb`, `logb`.

**Named constants in the expression:** `pi`, `e`, `ln_2`, `ln_10`, `sqrt2`, `phi`, `euler_gamma`.

### `cexpr!(expression, context)`

Same **context** as `expr!` (`Context` or the two tuples above). Same extra-precision loop: working precision starts at `p + WORD_BIT_SIZE`, internals use `RoundingMode::None`, the root is rounded once with `(p, rm)`, then each part is clamped with `check_complex_exponent_range`.

**`I`:** expands to `ExactComplex::i(p_wrk)` (`0 + 1i`). Lowercase `i` remains a variable. Do not name a variable `e` (`e` is Euler’s number, lifted as `e + 0i`).

**Cancellation:** add/sub use **two** `errs[]` slots, one for the real part and one for the imaginary part (`complex_cancel_bits` returns `(Option<usize>, Option<usize>)`). A subtraction can lose bits on both parts at different rates; a shared max would underestimate working precision. `p_wrk` is `p_rnd + sum(errs)`.

**Operators:** `+`, `-`, `*`, `/`. No `%`.

**Function leaves (complete list):**

`recip`, `sqrt`, `cbrt`, `root`, `ln`, `log2`, `log10`, `log`, `log1p`, `exp`, `exp2`, `exp10`, `expm1`, `pow`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `hypot`, `fma`, `mul_add`, `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`, `abs`, `arg`, `conj`, `ldexp`, `scalb`, `logb`.

**Named constants:** `I`, `pi`, `e`, `ln_2`, `ln_10`, `sqrt2`, `phi`, `euler_gamma` (reals except `I` are `x + 0i`).

**Reals lifted as `x + 0i`:** `abs`, `arg`, `logb`, literals, and the real constants. `ldexp` / `scalb` scale **both** parts by `2^n`. `hypot(z, w)` is principal `sqrt(z² + w²)`. `fma` / `mul_add` are extra-precision `a*b+c`, then one round. `root(z, n)` is `nth_root` (integer `n`).

**Branch cuts (principal values, observable):**

| Family | Cut / range |
| --- | --- |
| `arg` | `atan2(im, re)`, values in (−π, π] |
| `ln`, `log2`, `log10`, `log`, `log1p` | cut on (−∞, 0]; `ln(−1) = iπ` |
| `sqrt`, `cbrt`, `root` | cut on (−∞, 0]; `Re(sqrt) ≥ 0`; `sqrt(−1) = +i` |
| `pow` | `exp(w · ln(z))`, so the `ln` cut on the **base** |
| `asin` / `acos` / `atan` / `asinh` / `acosh` / `atanh` | principal branches of the usual identities |

**Trig / hyperbolic:** `sin(x+iy) = sin(x)cosh(y) + i cos(x)sinh(y)` (and the matching identities). The **complex** argument is never passed to `rem_pi`. Only a real component uses real `sin_cos` / `sinh_cosh` (those may reduce that real).

**Explicitly not in `cexpr!`:**

| Leaf | Why |
| --- | --- |
| `%` / `rem_pi` | remainder and π-reduction are real; complex trig uses the identities above |
| `atan2` | no standard two-complex analogue; use `arg` for `atan2(im, re)` |
| `erf` / `erfc` / `gamma` / `ln_gamma` / `bessel_j` | no complex kernel in this crate |

Literals and real constants are lifted as `x + 0i`. Variables may be `ExactComplex` or anything `FromExt` can wrap as a real.

**Not in `expr!`:** `frexp`, `ilogb`, `sin_cos`, `sinh_cosh`, `int`/`fract`/`ceil`/`floor`/`round`, `min`/`max`/`clamp`, `cmp`, `copysign`, `next_after`, `powi`/`powsi`, `abs`/`signum`, parse/format, raw parts, `nth_root` under the name `nth_root` (use `root`). Complex values use `cexpr!`, not `expr!`.

### `exact!("…")` and `fbig!("…")`

Compile-time parse of a **string literal** into an exact `ExactNum` (same expansion). Decimal literal syntax as implemented by the macro parser.

---

## 18. Constants cache

**`Consts::new() -> Result<Self, Error>`** — allocate progressive caches.

| Method | Constant |
| --- | --- |
| `pi(p, rm)` | π |
| `e(p, rm)` | e |
| `ln_2(p, rm)` | ln 2 |
| `ln_10(p, rm)` | ln 10 |
| `sqrt2(p, rm)` | √2 |
| `phi(p, rm)` | φ = (1+√5)/2 |
| `euler_gamma(p, rm)` | Euler–Mascheroni γ |
| `cache_info()` | `ConstCacheInfo` { `pi`, `e`, `ln2`, `ln10`, `sqrt2`, `phi`, `euler` } bit lengths |

**`CachedFBig`:** wrap a precomputed `ExactNum`; `new`, `cached_bit_len`, `inner`, `round(p, rm)`.

**`SharedConsts` (`std`):** mutex around `Consts`. `new()`, `with(|cc: &mut Consts| …)` (recovers poisoned mutex).

---

## 19. `Context` and `Contextable`

**`zenith_float::ctx::Context`:** holds precision, rounding mode, `Consts`, `emin`, `emax`.

| Method | |
| --- | --- |
| `new(p, rm, cc, emin, emax)` | |
| `to_raw_parts` | `(p, rm, Consts, emin, emax)` |
| `with_rounding_mode` | run a closure with a temporary rounding mode, then restore |
| `set_precision` / `set_rounding_mode` / `set_consts` / `set_emin` / `set_emax` | |
| `precision` / `rounding_mode` / `consts` / `emin` / `emax` | |
| `const_pi` / `const_e` / `const_ln2` / `const_ln10` / `const_sqrt2` / `const_phi` / `const_euler_gamma` | at context `p`/`rm` |
| `clone` | `Result<Self, Error>` (deep-clone cache) |

**`Contextable`:** `precision`, `rounding_mode`, `consts`, the seven `const_*` methods, `emin`, `emax`. Implemented for `Context` and the two tuples in §17.

---

## 20. `ExactComplex`

Cartesian `re + i·im` as two `ExactNum`s.

| API | |
| --- | --- |
| `new(re, im)` | |
| `re` / `im` | references |
| `zero(p)` / `one(p)` / `i(p)` | |
| `is_nan` | either part NaN |
| `conj` | |
| `abs(p, rm)` | modulus (`hypot`) |
| `arg(p, rm, cc)` | argument (`atan2`) |
| `add` / `sub` / `mul` / `div` | with `(p, rm)` |
| `exp` / `ln` / `sin` / `cos` / `tan` | with `(p, rm, cc)` |
| `sinh` / `cosh` / `tanh` | |
| `sqrt` / `pow` | principal branch |
| `asin` / `acos` / `atan` / `asinh` / `acosh` / `atanh` | principal branches |
| `log2` / `log10` / `log` / `log1p` | principal; `log2`/`log10`/`log` via `ln` ratios |
| `exp2` / `exp10` / `expm1` | via `exp` |
| `ldexp` / `scalb` | scale both parts by `2^n` |
| `cbrt` / `nth_root` | principal; `exp(ln / n)` |
| `hypot` | principal `sqrt(z² + w²)` |
| `fma` / `mul_add` | extra-precision `a*b+c` |
| `logb` | `logb(\|z\|)` as a real |
| `Add` `Sub` `Mul` `Div` | 128-bit `ToEven` like reals |

No `expr!` for complexes — use `cexpr!`. Serde: struct `{ "re", "im" }` of decimal strings.

---

## 21. `Ball` and `ziv_round`

**`Ball`:** enclosure `mid ± rad` (`rad` stored non-negative).

| Method | |
| --- | --- |
| `new(mid, rad)` | |
| `mid` / `rad` | |
| `add` / `mul` | interval arithmetic with an extra rounding ulp in the radius |
| `contains(x, p)` | `x` in `[mid−rad, mid+rad]`; NaN/Inf never contained |

**`ziv_round(p, rm, compute)`:** call `compute(working_p)` and `try_set_precision` until the rounding is unique or `MAX_PREC_RETRY` is exhausted (then NaN / `PrecisionRetryExhausted`).

---

## 22. `RadixFloat`

Binary kernel plus a **parse/format radix** tag.

| Method | |
| --- | --- |
| `new(value, radix)` | validates 2..=36 |
| `with_radix` | no extra check (for named radices) |
| `value` / `radix` / `into_inner` | |
| `parse(s, radix, p, rm, cc)` | |
| `format(rm, cc)` | using the stored radix |
| `Deref` / `DerefMut` | to `ExactNum` |
| `ExactNum::with_radix` | constructor |

---

## 23. `FromExt`

`FromExt<T>::from_ext(v, p, rm, cc)`:

- For types with `From<T> for ExactNum`: convert then `set_precision`.
- For `&str`: `parse` in decimal.

Used by `expr!` to lift variables and literals.

---

## 24. Serde (`serde` feature)

- **`ExactNum`:** serialize as a **decimal string** (`Display`). Deserialize from that string, or from JSON integers (`i64`/`u64`/`i128`/`u128`).
- **`ExactComplex`:** serialize struct with fields `re` and `im`; deserialize the same map.

---

## 25. Random (`random` feature)

| Item | |
| --- | --- |
| `DEFAULT_RANDOM_SEED` | `0x5EED_CAFE_BADC_0D00` |
| `random_seed()` | current seed |
| `reseed_random(seed)` | |
| `seeded_random::<T>()` | draw from the crate RNG |
| `ExactNum::random_normal` | |
| `ZENITH_TEST_SEED` | env override for tests / replay |

Unit tests default to the deterministic seed. Without reseeding, non-test `random` feature uses OS entropy.

---

## 26. Implementation notes (not extra public functions)

The kernel uses integer add/mul (schoolbook, Karatsuba/Toom, FFT at large sizes), Newton division, series and argument reduction for elementary/special functions, and a progressive constant cache. None of those algorithms are separate public types.

MPFR/`rug` appear only in **tests** (`mpfr-tests`), not as the evaluation engine.

---

## 27. What is not in this crate

These are load-bearing product choices, not a backlog:

- Hardware binary interchange types or converters (including a feature-gated module). Callers or a separate crate pack bits from `frexp` / `ilogb`.
- Using `expr!` for complex values (that macro is real-valued). Complex expressions are `cexpr!`.
- Symbolic CAS, formula rewriting, or host-application IR
- Remainder as a Rust `%` operator on `ExactNum`
- `Hash` / total order including NaN
- Assigning operators (`+=` …)

Callers choose precision `p` and rounding; this crate supplies the mechanism, not per-domain “formula class” policies.
