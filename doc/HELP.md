# Help: zenith-float

This is the user guide: **why** the API looks this way, recipes, mistakes, and FAQ. It is not a first-run tutorial ([GETTING_STARTED.md](GETTING_STARTED.md)) and not a method inventory ([LIBRARY.md](LIBRARY.md)).

You pick **how many bits** (`p`) and **how to round leftovers** (`rm`). Those two choices are the price of extra accuracy. The CPU’s hardware floating-point unit is not used.

---

## How these files fit together

| File | Role |
| --- | --- |
| [GETTING_STARTED.md](GETTING_STARTED.md) | Short path: depend, build a number, run `expr!`, print |
| **This file** | Why the API is shaped this way; recipes; mistakes; FAQ |
| [LIBRARY.md](LIBRARY.md) | Every public name |
| [EXPR.md](EXPR.md) | How `expr!` / `cexpr!` raise working precision and round once |
| [PRECISION.md](PRECISION.md) | Retry budget vs exponent |
| [README.md](../README.md) | Crate overview for contributors |

---

## The precision model

An `ExactNum` is a signed binary scientific number:

> sign × mantissa × 2^exponent

The mantissa is a list of machine **words** (`Word`: `u64` on a 64-bit host, `u32` on Thumb). You ask for a bit count `p`. The kernel **rounds `p` up to a whole number of words**. Ask for 1 bit and you get 64 bits on a 64-bit target. Ask for 65 bits and you get two words (128 bits).

**Why `p` is an argument, not a global.** A global `SOFT_PREC` would silently mix widths. A product of a 64-bit factor and a 1024-bit factor must say which width the *result* is. So `add`, `sin`, `gamma`, and the rest take `(p, rm)` (and `cc` when they need π or e).

**What happens if you do not pass `p`.** Rust `+ − * /` on `ExactNum` exist. They use a **fixed** default of 128 bits and `ToEven`. They are convenient for sketches. They are not “the current session precision.” If you need 256 bits, call `a.add(&b, 256, rm)` or put 256 in a `Context` and use `expr!`.

**Working precision.** Leaves that must be correctly rounded (or close) compute at `p_wrk = p + WORD_BIT_SIZE` and then cut back. `expr!` does the same for a whole tree and rounds **once** at the end. See [EXPR.md](EXPR.md).

**Exponent window.** `Context` also has `[emin, emax]`. Values outside that window become 0 or Inf. A huge unused range can make internal working precision grow for nothing. Use the narrowest window that still holds real answers.

---

## Rounding modes

After the kernel has more bits than you asked for, it **cuts** using `rm`.

| Mode | When to use it |
| --- | --- |
| `ToEven` | Default. Nearest; ties to even. Matches usual scientific / IEEE “round to nearest, ties to even.” |
| `Up` / `Down` | Directed rounding toward +∞ / −∞. Interval endpoints, conservative bounds. |
| `ToZero` / `FromZero` | Toward zero / away from zero. Integer conversion, some error analyses. |
| `ToOdd` | Nearest; ties to odd. Rare; some proofs prefer it to `ToEven`. |
| `None` | Do **not** cut. Extra bits may remain. Use on a method *chain*, then one `set_precision(p, ToEven)` at the end. Not a “final answer” mode. |

`ToEven` is the default on operators and `expr!` because it does not systematically bias a long sum the way “always up” does, and it matches the software IEEE kernel.

`expr!` computes with extra bits and `RoundingMode::None` internally, then one `set_precision(p, rm)`. That is **not** a promise that a long formula is correctly rounded in the IEEE fused sense.

---

## The constants cache

π, e, ln 2, ln 10, √2, φ, and Euler–Mascheroni γ are computed to the bits you ask for and **remembered**. A fresh `Consts::new()` on every `sin` recomputes π. Reuse one `Consts` (or one `Context`) in a loop.

```rust
use zenith_float::{Consts, ExactNum, RoundingMode};

let p = 256;
let rm = RoundingMode::ToEven;
let mut cc = Consts::new().expect("constants cache");
let x = ExactNum::from_u8(1, p);
let _ = x.sin(p, rm, &mut cc);
let _ = x.cos(p, rm, &mut cc); // same cache
```

With `std`, `SharedConsts` is a mutex around `Consts`. Recover a poisoned lock and run a closure:

```rust
use zenith_float::SharedConsts;
use zenith_float::{ExactNum, RoundingMode};

let sc = SharedConsts::new().expect("cache");
let y = sc.with(|cc| {
    ExactNum::from_u8(1, 128).sin(128, RoundingMode::ToEven, cc)
});
```

Each `Context` still owns its own `Consts` unless you build it from a cache you already have. `expr!` names: `pi`, `e`, `ln_2`, `ln_10`, `sqrt2`, `phi`, `euler_gamma`.

---

## `expr!` vs methods

**`expr!` / `cexpr!`.** Write the formula. Leaves are treated as exact. Working precision grows when add/sub cancel. One final round to the context `p` and `rm`. Good for `sin(pi / 6)`, nested `ln`/`exp`, complex `I * I`.

**What `expr!` does not guarantee.** Bit-identity with MPFR for a *composite* tree. Leaves have golds; the tree uses a cancellation heuristic. It does not rewrite algebra (no CAS). It does not pick `p` for you.

**Methods.** You pass `p` and `rm` every time. Required for pairs (`sin_cos`, `sinh_cosh`), IEEE split (`frexp`, `ilogb`), compensated helpers (`two_sum`, `fused_sum`, `polyval`), and anything not a macro leaf (`ceil`, `min`, `next_after`).

Mix them: `expr!` for the formula, then a method on the result.

---

## Complex arithmetic

Real formulas stay in `expr!`. Complex formulas use **`cexpr!`**. Same extra-precision loop; cancellation is tracked on the **real and imaginary parts separately**.

- Imaginary unit: **`I`** (capital). Lowercase `i` is a free variable.
- Do **not** name a variable `e` inside `cexpr!` — that token is Euler’s number.

**Why the leaf lists differ.** Some real operations have no honest complex form in this crate (`atan2`, `rem_pi`, `%`, real-order `bessel_j`, distributions). Complex `sin`/`cos` use the `x+iy` identities; they do **not** call `rem_pi` on a complex value.

**Principal cuts** (locked by golds):

- `ln` and `sqrt` on the non-positive real axis: `ln(−1) = iπ`, `sqrt(−1) = +i`
- \(K(m)\) and \({}_2F_1(\ldots;z)\) on \([1,+\infty)\)

See [LIBRARY.md](LIBRARY.md) §17.

---

## Special functions

Each public special has a `# Precision` rustdoc: algorithm (series, AGM, Carlson, …), the region it applies, named thresholds, ULP notes, and whether an MPFR oracle exists.

**`PrecisionRetryExhausted`.** Ziv / `MAX_PREC_RETRY` ran out of extra words. This is **not** a domain error. The value is `NaN` with that `Error`. Raise `p`, tighten the argument, or treat it as “could not certify.” Do not invent a number.

**Slow cases (typical).** Large `|e|` for `sin`/`cos` (argument reduction). Near-zero `Ci`. Growing Airy / Bessel on the wrong half-line. Read the method’s `# Precision` before blaming the cache.

Domain errors (`sqrt` of a negative *real*, `gamma` at a non-positive integer, …) are `NaN` + `InvalidArgument`.

---

## Arrays

`ExactNumArray` stores a shared `p` and a row-major buffer. A 1-D vector is shape `(1, n)`.

- `from_shape` is `None` if `rows * cols != vals.len()`.
- Elementwise `add` / `mul` require **equal shapes**. There is no NumPy broadcast.
- `matmul` is a sequential triple loop. Not BLAS. Not SIMD-float.
- `lu_decomp` / `svd_decomp` return `None` if singular, empty, non-finite, not converged, or a workspace `try_reserve_exact` fails (no abort).

`Ieee32Array` / `Ieee64Array` store integer bit patterns. Add/sub/mul/div/sqrt/fma use integer SIMD lanes and must match the scalar kernel bit-for-bit.

---

## Software IEEE

`Ieee32` / `Ieee64` are **interchange widths**: 24- and 53-bit significands stored as `u32` / `u64`. Arithmetic is the integer kernel, not a hardware add.

| Use | Type |
| --- | --- |
| Caller-chosen bits, specials, `expr!` | `ExactNum` |
| FFI / file / “this was a binary64 bit pattern” | `Ieee64::from_bits` / `to_bits` |
| Widen to work, then pack | `to_exact(p)` / `from_exact` |

A host that already has hardware bits converts **once** at the boundary. Do not sprinkle hardware types through the kernel.

---

## Recipes

Snippets assume `use zenith_float::{…}` as shown. They are meant to be copied, not to replace [LIBRARY.md](LIBRARY.md).

### 1. `sin(π/6)` via `expr!`

```rust
use zenith_float::{expr, Consts, RoundingMode};
use zenith_float::ctx::Context;
let mut ctx = Context::new(256, RoundingMode::ToEven, Consts::new().unwrap(), -10_000, 10_000);
let y = expr!(sin(pi / 6), &mut ctx);
```

### 2. `1/3` under two rounding modes

```rust
let down = ctx.with_rounding_mode(RoundingMode::Down, |ctx| expr!(1 / 3, ctx));
let up = ctx.with_rounding_mode(RoundingMode::Up, |ctx| expr!(1 / 3, ctx));
```

### 3. `hypot(3, 4)` is 5

```rust
use zenith_float::{expr, Consts, RoundingMode};
use zenith_float::ctx::Context;
let mut ctx = Context::new(256, RoundingMode::ToEven, Consts::new().unwrap(), -1000, 1000);
let h = expr!(hypot(3, 4), &mut ctx);
```

### 4. Fused sum `1+2+3`

```rust
use zenith_float::{ExactNum, RoundingMode};
let p = 256;
let rm = RoundingMode::ToEven;
let one = ExactNum::from_u8(1, p);
let six = ExactNum::fused_sum(&[one.clone(), ExactNum::from_u8(2, p), ExactNum::from_u8(3, p)], p, rm);
```

### 5. Horner `1 + 2x + 3x²` at `x = 2` is 17

```rust
let two = ExactNum::from_u8(2, p);
let seventeen = ExactNum::polyval(&[one.clone(), two.clone(), ExactNum::from_u8(3, p)], &two, p, rm);
```

### 6. `two_sum` then reconstruct

```rust
let (hi, lo) = one.two_sum(&two, p, rm);
let three = hi.add(&lo, p, rm);
```

### 7. Exact `1/3 + 1/6 = 1/2`

```rust
use zenith_float::ExactRational;
let half = ExactRational::from_i64(1, 3).add(&ExactRational::from_i64(1, 6));
assert_eq!(half, ExactRational::from_i64(1, 2));
```

### 8. `gcd(48, 18) = 6` and `17 ÷ 5`

```rust
use zenith_float::ExactInt;
assert_eq!(ExactInt::from_i64(48).gcd(&ExactInt::from_i64(18)), ExactInt::from_i64(6));
let (q, r) = ExactInt::from_i64(17).div_rem(&ExactInt::from_i64(5)).unwrap();
assert_eq!((q, r), (ExactInt::from_i64(3), ExactInt::from_i64(2)));
```

### 9. Text `0.1` is the rational `1/10`

```rust
use zenith_float::ExactRational;
let t = ExactRational::parse_exact("0.1").unwrap();
assert_eq!(t, ExactRational::from_i64(1, 10));
```

### 10. Software binary64 bits

```rust
use zenith_float::Ieee64;
let x = Ieee64::from_bits(0x3FF0_0000_0000_0000);
assert_eq!(x.to_bits(), Ieee64::from_i32(1).to_bits());
```

### 11. `2×2` `matmul`

```rust
use zenith_float::{ExactNum, ExactNumArray};
let p = 64;
let n = |k: u8| ExactNum::from_u8(k, p);
let a = ExactNumArray::from_shape(p, 2, 2, &[n(1), n(2), n(3), n(4)]).unwrap();
let b = ExactNumArray::from_shape(p, 2, 2, &[n(5), n(6), n(7), n(8)]).unwrap();
let c = a.matmul(&b).unwrap();
assert_eq!(c.get2(0, 0).unwrap().cmp(&n(19)), Some(0));
```

### 12. Shape mismatch is `None`

```rust
assert!(a.add(&ExactNumArray::from_shape(p, 1, 2, &[n(1), n(2)]).unwrap()).is_none());
```

### 13. LU of a small nonsingular matrix

```rust
let (l, u, _perm) = a.lu_decomp(p, rm).unwrap();
```

Singular or failed reserve → `None`, not a panic.

### 14. `I * I = −1`

```rust
use zenith_float::{cexpr, Consts, RoundingMode};
use zenith_float::ctx::Context;
let mut ctx = Context::new(256, RoundingMode::ToEven, Consts::new().unwrap(), -10_000, 10_000);
let z = cexpr!(I * I, &mut ctx);
```

### 15. Principal `sqrt(−1) = +i`

```rust
let z = cexpr!(sqrt(-1), &mut ctx);
```

### 16. A `Ball` around 0 contains `sin` of that ball at 0

```rust
use zenith_float::{Ball, Consts, ExactNum, RoundingMode};
let p = 128;
let rm = RoundingMode::ToEven;
let mut cc = Consts::new().unwrap();
let b = Ball::new(ExactNum::from_u8(0, p), ExactNum::from_u8(1, p));
assert!(b.sin(p, rm, &mut cc).contains(&ExactNum::from_u8(0, p), p));
```

### 17. Bisection: `sin` on `[3, 4]` meets π

```rust
use zenith_float::{bisect, root_default_tol, Consts, ExactNum, RoundingMode};
let p = 256;
let rm = RoundingMode::ToEven;
let mut cc = Consts::new().unwrap();
let tol = root_default_tol(p);
let root = bisect(|x, p, rm, cc| x.sin(p, rm, cc), &ExactNum::from_u8(3, p), &ExactNum::from_u8(4, p), &tol, p, rm, &mut cc);
```

### 18. Gauss–Legendre of `x²` on `[-1, 1]`

```rust
use zenith_float::{gauss_legendre, Consts, ExactNum, RoundingMode};
let p = 128;
let rm = RoundingMode::ToEven;
let mut cc = Consts::new().unwrap();
let i = gauss_legendre(|x, p, rm, _| x.mul(x, p, rm), &ExactNum::from_i8(-1, p), &ExactNum::from_u8(1, p), 8, p, rm, &mut cc);
```

### 19. Symmetric Hann window of length 4

```rust
use zenith_float::{hann_window, Consts, RoundingMode};
let w = hann_window(4, 64, RoundingMode::ToEven, &mut Consts::new().unwrap()).unwrap();
// w = [0, 3/4, 3/4, 0]
```

### 20. Inline binary record (≤ 64 mantissa bits)

```rust
use zenith_float::ExactNum;
let x = ExactNum::from_u8(2, 64);
let buf = x.to_inline_bytes().unwrap();
let y = ExactNum::from_inline_bytes(buf.as_bytes()).unwrap();
assert_eq!(x.cmp(&y), Some(0));
```

### 21. CSV: empty cell is `NAN`

```rust
use zenith_float::{Ieee64, Ieee64Array};
let a = Ieee64Array::from_csv_str("1,,3\n").unwrap();
assert_eq!(a.get(1).unwrap().to_bits(), Ieee64::NAN.to_bits());
```

### 22. SHA-256 of the empty slice

```rust
use zenith_float::sha256;
let h = sha256(b"");
```

### 23. Modular inverse `3⁻¹ ≡ 5 (mod 7)`

```rust
use zenith_float::{mod_inv, ExactInt};
assert_eq!(mod_inv(&ExactInt::from_i64(3), &ExactInt::from_i64(7)).unwrap(), ExactInt::from_i64(5));
```

### 24. `erf` at caller `p`

```rust
use zenith_float::{Consts, ExactNum, RoundingMode};
let mut cc = Consts::new().unwrap();
let y = ExactNum::from_u8(1, 256).erf(256, RoundingMode::ToEven, &mut cc);
```

### 25. `RoundingMode::None` then one `set_precision`

```rust
let mut acc = ExactNum::from_u8(0, 256);
acc = acc.add(&ExactNum::from_u8(1, 256), 256, RoundingMode::None);
acc.set_precision(128, RoundingMode::ToEven).unwrap();
```

### 26. Hex format after `expr!`

```rust
let s = y.format(zenith_float::Radix::Hex, RoundingMode::ToEven, ctx.consts()).unwrap();
```

### 27. `2^10` as `ExactInt`

```rust
use zenith_float::ExactInt;
assert_eq!(ExactInt::from_i64(2).pow(10), ExactInt::from_i64(1024));
```

### 28. FFT of an impulse

```rust
// ExactNumArray::fft: real row (1, n), n a power of two ≤ FFT_MAX_POINTS.
// Impulse [1,0,0,0] → [1,1,1,1] (unnormalized). See LIBRARY.md.
```

### 29. Newton for `x² − 2 = 0`

```rust
use zenith_float::{newton, root_default_tol, Consts, ExactNum, RoundingMode};
let p = 256;
let rm = RoundingMode::ToEven;
let mut cc = Consts::new().unwrap();
let two = ExactNum::from_u8(2, p);
let x0 = ExactNum::from_u8(1, p);
let root = newton(
    |x, p, rm, _| x.mul(x, p, rm).sub(&two, p, rm),
    |x, p, rm, _| x.add(x, p, rm),
    &x0,
    &root_default_tol(p),
    64,
    p,
    rm,
    &mut cc,
);
```

### 30. Parse a decimal string (NaN on failure)

```rust
use zenith_float::{Consts, ExactNum, Radix, RoundingMode};
let mut cc = Consts::new().unwrap();
let x = ExactNum::parse("1.25", Radix::Dec, 128, RoundingMode::ToEven, &mut cc);
assert!(x.err().is_none());
```

---

## Common mistakes

1. **`parse` vs `format`.** `parse` returns `ExactNum` (NaN on failure). `format` returns `Result<String, Error>`.
2. **Highest-first polynomial arrays.** `polyval` / `ExactNumPoly` are **lowest degree first**.
3. **Classical Dekker two-sum.** Here `p` and `rm` round `hi`; `lo` is the remainder. Reconstruct with `hi.add(&lo, p, rm)`.
4. **Variable named `e` in `cexpr!`.** That token is Euler’s number. Imaginary unit is `I`.
5. **Huge `emin` / `emax` “just in case.”** Tight windows keep working precision honest.
6. **New `Consts` on every call.** Reuse the cache (or `SharedConsts` across threads).
7. **Expecting `expr!` to match MPFR bit-for-bit on a whole tree.** Leaves are golded; composites use a cancellation heuristic.
8. **Hardware IEEE types in the kernel.** Convert at the edge with `from_bits` / `to_bits`.
9. **Using `+` and thinking you chose `p`.** Operators use 128 bits and `ToEven`.
10. **NumPy-style broadcast on arrays.** Shapes must match; otherwise `None`.
11. **CSV cells as human `1.0`.** `Ieee64Array` CSV is **unsigned bit patterns**. Empty / `nan` → `NAN`.
12. **`Hash` / total `Ord` including NaN.** Not provided. `cmp` returns `None` if unordered.
13. **`expr!` for a complex formula.** Use `cexpr!`.
14. **Treating `PrecisionRetryExhausted` as a domain error.** It is a Ziv budget. Raise `p` or accept `NaN`.
15. **`%` as a Rust operator on `ExactNum`.** Only `expr!` `%` (remainder) or the `rem` method.
16. **Assuming `lu_decomp` panics on a bad matrix.** It returns `None` (singular or reserve failure).
17. **Cloning `Context` without expecting `Result`.** Copying the constants cache can fail allocation.
18. **`to_inline_bytes` on a 256-bit number.** The 16-byte record holds 64 mantissa bits. Wider values use `to_bytes` (heap record).
19. **A new `ExactNumArray` cell precision vs the array `p`.** `from_shape` rounds every entry to the array `p`.
20. **Linking `libhdf5`.** Not in this crate. CSV and the binary record are the interchange paths today.

---

## FAQ

**Is this exact arithmetic?** Finite precision is still finite. `ExactNum` means the **limbs** are exact for what they store; rounding is explicit. It is not infinite-precision CAS.

**Why bits instead of decimal digits?** The mantissa is binary. Decimal is a display / parse radix.

**Do I need `expr!`?** No. Methods are the full API. Macros are for readable formulas.

**What does `expr!` guarantee?** Extra working precision, one final round to `(p, rm)`, exponent clamping. Not fused IEEE correctness of an arbitrary tree.

**What is `p_wrk`?** `p + WORD_BIT_SIZE`. Leaves and `expr!` start there.

**Why is `ToEven` the default?** Unbiased ties; matches the software IEEE kernel.

**When do I use `RoundingMode::None`?** Intermediate steps of a method chain; then one `set_precision`.

**Why does every special take `cc`?** So π and friends are cached. There is no global cache.

**What is `SharedConsts`?** `std` mutex around `Consts` for threads. `with(|cc| …)` recovers a poisoned lock.

**Why is `cexpr!` a different macro?** Real and complex cancellation, cuts, and leaf sets are not the same.

**Why no `atan2` in `cexpr!`?** No honest complex `atan2` leaf in this crate. Use `arg` or real methods on parts.

**What is a principal cut?** A chosen jump of a multi-valued function. Ours are listed above and in LIBRARY §17.

**`sqrt` of a negative real?** `ExactNum`: `NaN`. `cexpr!`: `+i` times the positive sqrt.

**What is `PrecisionRetryExhausted`?** Ziv ran out of extra words. Not a domain error.

**How do I recover an error from `NaN`?** `ExactNum::err()` → `Option<Error>`.

**Why Inf?** Exponent overflow or a tight `emax`. Widen only if the math needs it.

**Why subnormals?** Tiny values at `EXPONENT_MIN`; leading mantissa bits can be zero.

**`Ieee64` vs `ExactNum`?** Fixed 53-bit interchange vs caller `p`.

**Are SIMD lanes hardware floats?** No. Integer `u32`/`u64` lanes. Must match the scalar kernel.

**Why does elementwise add return `None`?** Shapes differ.

**Why is `matmul` slow?** Triple loop. No BLAS.

**What does `lu_decomp` return on OOM?** `None` (failed `try_reserve_exact`).

**How do I serialize?** `serde` feature: strings carry `@p=`. IEEE arrays serialize **bits**. Or `to_bytes` / `to_inline_bytes`. Or CSV.

**Why `@p=` in JSON?** So a 256-bit value does not come back at 128 bits on another target.

**CSV `1` is not the number one?** For `Ieee64Array`, `1` is bit pattern 1 (a subnormal), not `1.0`.

**Is HDF5 supported?** No C `libhdf5`. A tiny owned subset is leftover. Use CSV or `to_bytes`.

**`no_std`?** Allocator required. Formatting traits and `SharedConsts` need `std`. Thumb CI is leftover (`lazy_static` still wants `std`).

**Can I plot?** Not in this crate. Format a string or pack bits in another crate.

**Is there computer algebra?** No. That is Accumath, not zenith-float.

**`Hash` of a float?** Not if NaN is involved. No total order including NaN.

**Assigning `+=`?** Not implemented. Use `add` and assign.

**`%` outside `expr!`?** Use `rem`.

**Default operator precision?** 128 bits, `ToEven`.

**Compile-time literals?** `exact!("1.25")` or `fbig!("1.25")`.

**Bases above 10?** Exponent marker is `_e` so `e` can be a digit.

**`Ball` vs correctly rounded `ExactNum`?** A ball is an enclosure. `ziv_round` uses balls to certify a unique midpoint.

**`contains` prove rounding?** No. It tests membership at precision `p`.

**MPFR at runtime?** Never. `mpfr-tests` is an optional oracle, Linux x86_64 only.

**Integer SIMD leftover?** Closed. Add/sub/mul/div/sqrt/fma are integer-lane SIMD and bit-identical to the scalar kernel.

**Where is the method list?** [LIBRARY.md](LIBRARY.md).

**Where is the first example?** [GETTING_STARTED.md](GETTING_STARTED.md).

---

## What this library does not do

These are choices, not missing tickets:

- Hardware floating-point arithmetic
- `expr!` for complex values (`cexpr!` instead)
- Symbolic CAS, formula rewriting, or an expression IR
- Choosing `p` for a physics formula
- Total order / `Hash` that includes NaN
- A C HDF5 dependency

---

## If something looks wrong

**The answer is NaN.** Print `y.err()`. Division by zero, a domain error, allocation failure, or `PrecisionRetryExhausted` all land here.

**The answer is Inf.** The exponent window is too tight, or the value overflowed.

**Digits change when I raise `p`.** You asked for more bits. If they never settle, the formula is cancelling. `expr!` adds working bits; a `None` chain plus final `set_precision` is the other pattern.

**Compile error about `e` in `cexpr!`.** Rename the variable.

**I cloned `Context` and got `Result`.** The constants cache copy can fail.

---

## Next reading

- [GETTING_STARTED.md](GETTING_STARTED.md) — first run
- [LIBRARY.md](LIBRARY.md) — every public name
- [EXPR.md](EXPR.md) — `expr!` / `cexpr!` rounding contract
- [PRECISION.md](PRECISION.md) — retries and working precision
