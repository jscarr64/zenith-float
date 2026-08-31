# Help: zenith-float

This is the long guide. If you only want a working snippet, start with [GETTING_STARTED.md](GETTING_STARTED.md). If you need every method name, use [LIBRARY.md](LIBRARY.md).

You do not need to be a numerical analyst to use this crate. You do need to be willing to say **how many bits** you want and **how to round** leftovers. Those two choices are the price of extra accuracy.

---

## How these files fit together

| File | What it is for |
| --- | --- |
| [GETTING_STARTED.md](GETTING_STARTED.md) | Short path: depend, build a number, run `expr!`, print |
| **This file** | Why the API looks this way, recipes, mistakes, “what is not here” |
| [LIBRARY.md](LIBRARY.md) | Complete public inventory |
| [EXPR.md](EXPR.md) | How `expr!` raises working precision and rounds once at the end |
| [PRECISION.md](PRECISION.md) | How much extra work the library is allowed to do |
| [README.md](../README.md) | Crate overview and error-bound notes for contributors |

---

## The idea in plain language

A pocket calculator stores a **fixed** number of digits. If you ask it for π, it gives you something like 3.14159265 and **stops**. Multiply that by a huge number and the error grows with it.

**zenith-float** still stores a finite number of digits, but **you** pick how many (as **bits**, not decimal places). The arithmetic is ordinary integer arithmetic on those bits. The CPU’s hardware floating-point unit is **not** used. IEEE binary32/binary64 are software types `Ieee32` / `Ieee64` (`from_bits` / `to_bits`). Arbitrary precision is `ExactNum`.

Think of an `ExactNum` as a signed scientific-notation number in **base 2**:

> sign × mantissa × 2^exponent

The mantissa is a list of machine **words** (64-bit on a typical PC). Precision is a **bit count**, then rounded **up** to a whole number of words. Ask for 1 bit and you get a full word (64 bits on 64-bit targets). Ask for 65 bits and you get two words (128 bits).

Almost every operation also takes a **rounding mode**: after the library has more bits than you asked for, it **cuts** the result back to your precision using that rule.

---

## Words you will see

| Word | Meaning |
| --- | --- |
| **Precision `p`** | How many bits of mantissa you asked for (then rounded up to the word size). |
| **Rounding mode `rm`** | What to do with leftover bits (nearest-even is `ToEven`, a good default). |
| **Limb / word** | One integer chunk of the mantissa. |
| **Exponent** | The power of two in scientific notation. Too big → infinity; too small → zero or a **subnormal**. |
| **Context** | A bundle: `p`, `rm`, a constants cache, and allowed exponent range. `expr!` needs this. |
| **`Consts`** | A cache so π, e, ln 2, … are not recomputed from scratch every call. |
| **`inexact`** | A sticky flag: some bits were already rounded, or an input was already inexact. |
| **NaN** | “Not a number.” Errors at the API (bad arguments, division by zero, allocation failure) usually become this. |
| **`ExactComplex`** | Two `ExactNum`s: real part and imaginary part. |

---

## Install

In `Cargo.toml`:

```toml
[dependencies]
zenith-float = "0.1.0"
```

`std` is on by default (pretty-printing and `FromStr`). For `no_std` you still need an allocator:

```toml
zenith-float = { version = "0.1.0", default-features = false }
```

Depend on **`zenith-float`**, not `zenith-float-num`. The kernel crate is an implementation detail.

---

## Precision: how many bits is “enough”?

There is no universal answer. A useful picture:

- **53 bits** is roughly “calculator double” width (about 15–16 decimal digits of *binary* precision).
- **256 bits** is already far past everyday decimal homework.
- **1024 bits** is “I want π printed as a wall of digits.”

More bits means **more memory and more time**. Multiplication and `sin` get more expensive as `p` grows. Pick the smallest `p` that makes your printed result stable, then stop.

The exponent range on `Context` is a **safety window**. Results outside `[emin, emax]` become 0 or infinity. Use the **narrowest** window that still holds real answers. A huge unused range can make internal working precision grow for nothing.

---

## Rounding: what happens to leftover bits

Imagine you must keep 2 decimal digits and you have 1.255. Different rules pick 1.26 or 1.25.

| Mode | Everyday reading |
| --- | --- |
| `ToEven` | Nearest; ties go to even (usual scientific default). |
| `Up` / `Down` | Toward +∞ / −∞ (directed). |
| `ToZero` / `FromZero` | Toward zero / away from zero. |
| `ToOdd` | Nearest; ties to odd. |
| `None` | Do **not** cut back; extra bits may remain (faster internals, not a “final answer” mode). |

`expr!` computes with extra working bits, then rounds **once** to the context’s `p` and `rm`. That is **not** a promise that a long formula is correctly rounded in the IEEE fused sense. See [EXPR.md](EXPR.md).

---

## Context: a settings sheet for `expr!`

`expr!` is a macro that looks like math: `sin(pi / 6)`. Under the hood it still needs:

1. Precision  
2. Rounding mode  
3. A `Consts` cache (π and friends)  
4. Exponent limits  

```rust
use zenith_float::{expr, Consts, RoundingMode};
use zenith_float::ctx::Context;

let mut ctx = Context::new(
    256,
    RoundingMode::ToEven,
    Consts::new().expect("constants cache"),
    -10_000,
    10_000,
);

let y = expr!(sin(pi / 6), &mut ctx);
```

A one-shot tuple also works: `(p, rm, &mut cc)` or `(p, rm, &mut cc, emin, emax)`.

To try a **different rounding mode** for one sub-expression without losing the rest of the context:

```rust
let down = ctx.with_rounding_mode(RoundingMode::Down, |ctx| expr!(1 / 3, ctx));
```

The previous mode comes back when the closure returns.

---

## Two styles of code

**Style A — `expr!` / `cexpr!`.** Write the formula. Good for `sin`, `ln`, `pow`, nested arithmetic. Inputs are treated as exact at the start of the tree.

**Style B — methods on `ExactNum`.** You pass `p` and `rm` on every call: `a.add(&b, p, rm)`. Good when you need a **pair** of results (`sin_cos`, `sinh_cosh`), IEEE-style split (`frexp`, `ilogb`), or compensated helpers (`two_sum`, `fused_sum`, `polyval`).

You can mix them: compute a value with `expr!`, then call a method on the result.

Rust `+ − * /` on `ExactNum` exist too, with a **fixed default** precision (128 bits, `ToEven`). They are convenient, not a substitute for choosing `p` when accuracy matters.

---

## Making numbers

```rust
use zenith_float::{exact, Consts, ExactNum, RoundingMode};

let p = 256;
let rm = RoundingMode::ToEven;
let mut cc = Consts::new().expect("constants cache");

let a = ExactNum::from_u32(3, p);
let b = ExactNum::parse("0.5", zenith_float::Radix::Dec, p, rm, &mut cc);
let c = a.add(&b, p, rm); // 3.5
let x = exact!("1.25");
```

- Integers: `from_u8` … `from_u128` and matching `from_i*`.  
- Text: `parse` with a `Radix` (2 through 36). Named: `Bin`, `Oct`, `Dec`, `Hex`.  
- Compile-time decimals: `exact!("…")` or the alias `fbig!("…")`. You can write `zenith_float::exact!("1.25")` without importing the macro.

**`parse` does not return `Result`.** A bad string becomes **NaN**. Call `ExactNum::err()` to get the `Error`. **`format` does return `Result<String, Error>`.** That mismatch is easy to miss.

---

## Reading and writing text

After `expr!` returns, you own an `ExactNum`. The `&mut ctx` borrow from the macro is over.

```rust
println!("{}", y);      // decimal (needs `std`)
println!("{:e}", y);    // scientific

let s = y
    .format(zenith_float::Radix::Hex, RoundingMode::ToEven, ctx.consts())
    .expect("format");
```

`ctx.consts()` is the same cache you put in `Context`. Formatting π-related conversions needs it.

For bases where `e` is a digit (above 10), exponents use `_e` so `e` can still be a digit character. See [LIBRARY.md](LIBRARY.md) §16.

---

## When a result is not a “normal” number

These are **software** values, not hardware registers:

| Value | Typical cause |
| --- | --- |
| `+Inf` / `−Inf` | Exponent overflow (`Error::ExponentOverflow`). |
| `NaN` | Division by zero, invalid argument, allocation failure. |
| Subnormal | Tiny numbers at the smallest allowed exponent; leading mantissa bits can be zero. |

`ExactNum::err()` is for **NaN** with an associated error. Check `is_nan()`, `is_inf()`, `is_finite()` as you would in any numerical library.

---

## The constants cache

π, e, ln 2, ln 10, √2, φ, and Euler–Mascheroni γ are computed to the bits you ask for and **remembered**. Reuse one `Consts` (or one `Context`) in a loop. Creating a new cache on every call wastes work.

With `std`, `SharedConsts` wraps the cache in a mutex for threads. Each `Context` still owns its own `Consts` unless you design otherwise.

`expr!` names inside the formula: `pi`, `e`, `ln_2`, `ln_10`, `sqrt2`, `phi`, `euler_gamma`.

---

## What you can write inside `expr!`

Allowed: variables, integer/string literals, unary minus, `+ − * / %`, function calls, parentheses.

**Functions (real `expr!`):**  
`recip`, `sqrt`, `cbrt`, `root`, `ln`, `log2`, `log10`, `log`, `log1p`, `exp`, `exp2`, `exp10`, `expm1`, `pow`, `rem_pi`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`, `hypot`, `fma`, `mul_add`, `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`, `erf`, `erfc`, `gamma`, `ln_gamma`, `digamma`, `gammainc`, `gammainc_upper`, `ei`, `si`, `ci`, `li`, `fresnel_s`, `fresnel_c`, `bessel_j`, `bessel_j_nu`, `bessel_y`, `bessel_i`, `bessel_k`, `elliptic_k`, `elliptic_e`, `elliptic_e_inc`, `elliptic_f`, `elliptic_pi`, `elliptic_pi_inc`, `legendre_p`, `legendre_p_assoc`, `hypergeom_2f1`, `betainc`, `ldexp`, `scalb`, `logb`.

**Not in `expr!` (use methods):** `frexp`, `ilogb`, `sin_cos`, `sinh_cosh`, `ceil` / `floor` / `round`, `min` / `max`, `copysign`, `next_after`, `polyval`, `two_sum`, and friends.

**`%` in `expr!`** is remainder (`rem`). There is no Rust `%=` and no `%` **operator** on `ExactNum` values outside the macro.

Bessel \(J,Y,I,K\) (real order), elliptic \(K,E,\Pi\), Legendre \(P_n/P_n^m\), \({}_2F_1\), and regularized incomplete beta are in the crate.

---

## Methods when the formula macro is the wrong tool

Paired trig/hyperbolic, exponent split, and compensated primitives take **explicit** `p` and `rm`.

**`polyval`:** Horner scheme. **`coeffs[0]` is the constant term** (lowest degree first). `[1, 2, 3]` at `x = 2` is `1 + 2x + 3x² = 17`. That is the opposite of “highest power first” arrays in some other libraries.

**`two_sum` / `two_product`:** not a hardware Dekker two-sum. The library forms the **exact** limb sum or product, then **`p` and `rm` round only `hi`**. `lo` is the remainder. Put them back together with `hi.add(&lo, p, rm)` (not `add_full_prec` with precision 0).

**`fused_sum` / `fused_dot`:** extra working precision, **one** final round. Better than a long chain of ordinary `add` when many terms cancel.

```rust
use zenith_float::{ExactNum, RoundingMode};

let p = 256;
let rm = RoundingMode::ToEven;
let one = ExactNum::from_u32(1, p);
let two = ExactNum::from_u32(2, p);
let three = ExactNum::from_u32(3, p);

let six = ExactNum::fused_sum(&[one.clone(), two.clone(), three.clone()], p, rm);
let seventeen = ExactNum::polyval(&[one.clone(), two.clone(), three.clone()], &two, p, rm);
let (hi, lo) = one.two_sum(&two, p, rm);
let three_again = hi.add(&lo, p, rm);
```

---

## Complex numbers

Real formulas stay in `expr!`. Complex formulas use **`cexpr!`**. Same extra-precision idea; cancellation is watched on **both** the real and imaginary parts.

The imaginary unit in the expression is **`I`** (capital i). That leaves lowercase `i` free as a variable name.

**Do not name a variable `e` inside `cexpr!`.** In that macro, `e` is Euler’s number, the same as in `expr!`.

```rust
use zenith_float::{cexpr, Consts, RoundingMode};
use zenith_float::ctx::Context;

let mut ctx = Context::new(
    256,
    RoundingMode::ToEven,
    Consts::new().expect("constants cache"),
    -10_000,
    10_000,
);
let z = cexpr!(I * I, &mut ctx); // −1 + 0i
```

`cexpr!` has no `%`, no `atan2`, and no `rem_pi`. Leaves include roots (`sqrt`, `cbrt`, `root`), logs, exps, `pow`, circular and hyperbolic functions and inverses, `hypot` / `fma` / `mul_add`, `abs` / `arg` / `conj`, `ldexp` / `scalb` / `logb`, `erf` / `erfc` / `gamma` / `ln_gamma` / `digamma`, `ei` / `si` / `ci` / `li` / Fresnel, `bessel_j_nu` / `bessel_y` / `bessel_i` / `bessel_k`, and elliptic `elliptic_k` / `elliptic_e` / `elliptic_f` / `elliptic_e_inc` / `elliptic_pi` / `elliptic_pi_inc`. Principal branch cuts: `ln` and `sqrt` on the non-positive real axis (`ln(−1) = iπ`, `sqrt(−1) = +i`); \(K(m)\) on \([1,+\infty)\). Add/sub cancellation is tracked **separately** on the real and imaginary parts. See [LIBRARY.md](LIBRARY.md) §17.

---

## What this library does not do

These are **choices**, not missing tickets:

- Use hardware floating-point arithmetic. IEEE widths are `Ieee32` / `Ieee64` (`from_bits` / `to_bits`).  
- Treat `expr!` as complex (use `cexpr!`).  
- Computer algebra (rewrite formulas, Risch, expression IR).  
- Decide “how many bits this physics formula needs.” You choose `p` and `rm`.  
- Total order / `Hash` that includes NaN.  

---

## Common mistakes

1. **`parse` vs `format`.** Parse → `ExactNum` (NaN on failure). Format → `Result`.  
2. **Highest-first polynomial arrays.** `polyval` is lowest-first.  
3. **Classical two-sum.** Here `p` and `rm` round `hi`.  
4. **`cexpr!` variable `e`.** Use another name; `e` is the constant. Imaginary unit is `I`.  
5. **Huge `emin`/`emax` “just in case.”** Tight bounds keep working precision honest.  
6. **New `Consts` every call.** Reuse the cache.  
7. **Expecting `expr!` to be bit-identical to MPFR** for a whole tree. Leaves are tested; composite expressions use a cancellation heuristic.  
8. **Expecting hardware IEEE types.** Use `Ieee32::from_bits` / `Ieee64::to_bits`, or pack from `ExactNum` via `frexp` / `ilogb`.

---

## If something looks wrong

**The answer is NaN.** Print `y.err()` if it is `Some`. Division by zero, a domain error (for example `sqrt` of a negative real), or running out of memory all land here.

**The answer is Inf.** The exponent window is too tight, or the value really overflowed. Widen `emax` only if the mathematics needs it.

**Digits change when I raise `p`.** That can be normal: you asked for more bits. If they **never** settle, the formula may be cancelling (close numbers subtracting). `expr!` tries to add working bits; a method chain with `RoundingMode::None` until a final `set_precision` is the other pattern.

**Compile error about `e` in `cexpr!`.** Rename your variable. `e` is Euler’s number.

**I cloned `Context` and it returned `Result`.** Cloning copies the constants cache; that can fail allocation.

---

## FAQ

**Is this “exact” arithmetic?** Finite precision is still finite. The type is named `ExactNum` because the **limb representation** is exact for what it stores; rounding is explicit. It is not infinite-precision CAS.

**Why bits instead of decimal digits?** The mantissa is binary. Decimal digits are a display choice (`Radix::Dec`).

**Can I plot this?** Not inside this crate. Format a string, or in another crate pack bits from `frexp` / `ilogb`.

**Do I need `expr!`?** No. Methods are the full API. The macros are for readable formulas.

**What about `no_std`?** Allocator required. Formatting traits need the `std` feature.

**Where is the full method list?** [LIBRARY.md](LIBRARY.md).

---

## Next reading

- [GETTING_STARTED.md](GETTING_STARTED.md) if you have not run a first example  
- [LIBRARY.md](LIBRARY.md) for every public name  
- [EXPR.md](EXPR.md) if you need the `expr!` rounding contract  
- [PRECISION.md](PRECISION.md) if working precision or retries matter  

