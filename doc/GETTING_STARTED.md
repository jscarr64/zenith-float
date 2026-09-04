# Getting started with zenith-float

This is a short path from a new crate dependency to a computed, formatted result. For a longer explanation (precision, rounding, recipes, mistakes), see [HELP.md](HELP.md). For every public name, see [LIBRARY.md](LIBRARY.md). For how `expr!` rounds, see [EXPR.md](EXPR.md).

All arithmetic is software integer limbs. There are no hardware binary interchange types in this crate.

## 1. Depend on the crate

```toml
[dependencies]
zenith-float = "1.0"
```

`std` is on by default (formatting and `FromStr`). For `no_std` plus an allocator:

```toml
zenith-float = { version = "1.0", default-features = false }
```

## 2. Construct numbers

Precision is a bit count, rounded up to the word size (64 bits on 64-bit targets). Rounding is explicit on almost every operation.

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

`parse` returns `ExactNum`, not `Result`. A failed parse is `NaN`; `ExactNum::err()` recovers the `Error`. `from_u32` is one of `from_u8` / `from_u16` / `from_u32` / `from_u64` / `from_u128` (and the signed `from_i*` counterparts). Compile-time decimals: `exact!("…")` or the alias `fbig!("…")`. You can also write `zenith_float::exact!("1.25")` without importing the macro.

Infinities and NaN are software values: `INF_POS`, `INF_NEG`, `NAN`. Errors (overflow, division by zero, bad arguments, allocation) become NaN; `ExactNum::err()` recovers the `Error`.

## 3. Wire a `Context` and use `expr!`

`expr!` needs a context: precision, rounding mode, a `Consts` cache, and an exponent range. Use the **narrowest** exponent range that still holds your results; working precision can grow with the exponent.

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

A temporary tuple works too: `(p, rm, &mut cc)` or `(p, rm, &mut cc, emin, emax)`.

`expr!` raises working precision to fight cancellation. It does **not** by itself guarantee a correctly rounded final bit pattern; see [EXPR.md](EXPR.md).

To evaluate a sub-expression with a different rounding mode:

```rust
let down = ctx.with_rounding_mode(RoundingMode::Down, |ctx| expr!(1 / 3, ctx));
```

The previous mode is restored when the closure returns.

## 4. Format output

`y` below is an owned `ExactNum` after `expr!` returns, so the `&mut ctx` borrow from the macro has ended. `format` takes `&self` and `&mut Consts`; `ctx.consts()` is that cache. `format` returns `Result<String, Error>` (unlike `parse`).

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

println!("{}", y);
println!("{:e}", y);

let s = y
    .format(zenith_float::Radix::Hex, RoundingMode::ToEven, ctx.consts())
    .expect("format");
```

## 5. When to call methods instead of `expr!`

Use `ExactNum` methods when you need paired results (`sin_cos`, `sinh_cosh`), IEEE-style split (`frexp`, `ilogb`), or compensated primitives (`two_sum`, `two_product`, `fused_sum`, `fused_dot`, `polyval`).

`fused_sum` / `fused_dot` add (or multiply-then-add) at extra working precision and round **once**. `polyval` is Horner with **lowest degree first**: `coeffs[0]` is `a0` in `a0 + a1 x + a2 x² + …`. `two_sum` / `two_product` compute the exact limb sum or product, then round that value to `p` bits with `rm` to get `hi`; `lo` is the remainder. Reconstruct with `hi.add(&lo, p, rm)` (not a hardware-float Dekker two-sum, which has no `p`/`rm`).

```rust
use zenith_float::{ExactNum, RoundingMode};

let p = 256;
let rm = RoundingMode::ToEven;
let one = ExactNum::from_u32(1, p);
let two = ExactNum::from_u32(2, p);
let three = ExactNum::from_u32(3, p);

let six = ExactNum::fused_sum(&[one.clone(), two.clone(), three.clone()], p, rm);
// polyval: [a0, a1, a2] = 1 + 2x + 3x²; at x = 2 that is 17
let seventeen = ExactNum::polyval(&[one.clone(), two.clone(), three.clone()], &two, p, rm);

let (hi, lo) = one.two_sum(&two, p, rm); // hi is 1+2 rounded to p bits
let three_again = hi.add(&lo, p, rm);
```

For **complex** expressions use `cexpr!` (documented with `expr!` in [LIBRARY.md](LIBRARY.md) §17). Same extra-precision loop; cancellation on both parts; imaginary unit `I`. `expr!` stays real-valued.

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

Hardware IEEE registers are not used. For a fixed-width interchange type see §6. For arrays, rationals, and balls see the sections after that.

## 6. Software `Ieee32` / `Ieee64` (interchange)

These store IEEE-754 binary32 / binary64 as integer bits. Arithmetic is the same integer kernel as `ExactNum`, cut to 24- or 53-bit significands. Use them at an FFI or file boundary, not as a second calculator.

```rust
use zenith_float::{ExactNum, Ieee64};

let bits: u64 = 0x3FF0_0000_0000_0000; // 1.0
let x = Ieee64::from_bits(bits);
assert_eq!(x.to_bits(), bits);

let wide = x.to_exact(128);
let back = Ieee64::from_exact(&wide);
assert_eq!(back.to_bits(), bits);

let one = ExactNum::from_u8(1, 64);
assert_eq!(Ieee64::from_exact(&one).to_bits(), Ieee64::from_i32(1).to_bits());
```

`from_bits` / `to_bits` are the only way in or out of a host that already has hardware bits. Details: [LIBRARY.md](LIBRARY.md) §7b.

## 7. 2-D arrays

A 1-D vector is shape `(1, n)`. `from_shape` returns `None` if `rows * cols` is not `vals.len()`. Elementwise add/mul require matching shapes. `matmul` is a sequential triple loop (not BLAS).

```rust
use zenith_float::{ExactNum, ExactNumArray, Ieee64, Ieee64Array};

let p = 64;
let n = |k: u8| ExactNum::from_u8(k, p);
let a = ExactNumArray::from_shape(p, 2, 2, &[n(1), n(2), n(3), n(4)]).unwrap();
let b = ExactNumArray::from_shape(p, 2, 2, &[n(5), n(6), n(7), n(8)]).unwrap();
let c = a.matmul(&b).unwrap();
assert_eq!(c.get2(0, 0).unwrap().cmp(&n(19)), Some(0));
assert!(a
    .add(&ExactNumArray::from_shape(p, 2, 3, &[n(1), n(1), n(1), n(1), n(1), n(1)]).unwrap())
    .is_none());

let u = Ieee64Array::from_shape(2, 2, &[Ieee64::from_i32(1); 4]).unwrap();
let v = u.matmul(&u).unwrap();
assert_eq!(v.get2(0, 0).unwrap().to_bits(), Ieee64::from_i32(2).to_bits());
```

CSV and the 16-byte binary record are in [LIBRARY.md](LIBRARY.md) §7b. Decompositions (`lu_decomp`, `svd_decomp`) return `None` if singular or a workspace reserve fails.

## 8. Specials at caller-chosen `(p, rm, cc)`

There is no process-wide precision. Every special takes the bit count, the rounding mode, and a `Consts` cache (when the algorithm needs π, e, or a similar constant).

```rust
use zenith_float::{Consts, ExactNum, RoundingMode};

let p = 256;
let rm = RoundingMode::ToEven;
let mut cc = Consts::new().expect("constants cache");
let x = ExactNum::from_u8(1, p);
let y = x.sin(p, rm, &mut cc); // sin(1)
```

Reuse `cc` in a loop. A new cache on every call recomputes π. Domain errors become `NaN`; `y.err()` is the `Error`. The same `(p, rm, cc)` arguments appear on `ExactNumArray` elementwise specials. Jacobi `sn`/`cn`/`dn` take parameter \(m=k^2\in[0,1]\); the inverse of `sn` is `elliptic_f`. See [HELP.md](HELP.md) recipe 32 and `# Precision` on each method.

## 9. `ExactRational` and `ExactInt`

`ExactRational` is reduced `num/den` (not a float). `0.1` parsed with `parse_exact` is `1/10`. `ExactInt` is a signed limb integer: factorials, `gcd`, `div_rem`. Convert to `ExactNum` only when you need a rounded float.

```rust
use zenith_float::{ExactInt, ExactRational};

let half = ExactRational::from_i64(1, 3).add(&ExactRational::from_i64(1, 6));
assert_eq!(half, ExactRational::from_i64(1, 2));

let g = ExactInt::from_i64(48).gcd(&ExactInt::from_i64(18));
assert_eq!(g, ExactInt::from_i64(6));
let (q, r) = ExactInt::from_i64(17).div_rem(&ExactInt::from_i64(5)).unwrap();
assert_eq!((q, r), (ExactInt::from_i64(3), ExactInt::from_i64(2)));
```

Inventory: [LIBRARY.md](LIBRARY.md) §7b.

## 10. Complex specials in `cexpr!`

`expr!` is real. `cexpr!` is complex: cancellation is tracked on each part. The imaginary unit is `I`. Do not name a variable `e` (that token is Euler’s number).

Principal cuts: `ln` and `sqrt` on the non-positive real axis (`ln(−1) = iπ`, `sqrt(−1) = +i`); \(K(m)\) and \({}_2F_1(\ldots;z)\) on \([1,+\infty)\).

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
let z = cexpr!(sqrt(-1), &mut ctx); // +i
let w = cexpr!(erf(I), &mut ctx);
```

Leaves include `erf` / `Γ` / `Ei` / Bessel / elliptic / `{}_2F_1`. Full list and cuts: [LIBRARY.md](LIBRARY.md) §17.

## 11. `Ball` enclosures

A `Ball` is `mid ± rad`. Transcendentals return a wider ball that still contains the image of the input disk (Lipschitz padding, not a tight interval library). `contains` is a software test at precision `p`. It does **not** prove a correctly rounded `ExactNum` by itself; `ziv_round` uses balls to certify a unique midpoint.

```rust
use zenith_float::{Ball, Consts, ExactNum, RoundingMode};

let p = 128;
let rm = RoundingMode::ToEven;
let mut cc = Consts::new().expect("constants cache");
let mid = ExactNum::from_u8(0, p);
let rad = ExactNum::from_u8(1, p);
let b = Ball::new(mid, rad);
let s = b.sin(p, rm, &mut cc);
assert!(s.contains(&ExactNum::from_u8(0, p), p));
```

`ComplexBall` is a disk in \(\mathbb{C}\). Contract: [LIBRARY.md](LIBRARY.md) (Ball / `ziv_round`).

## 12. Next documents

| Doc | Audience |
| --- | --- |
| [HELP.md](HELP.md) | Longer tutorial, recipes, FAQ |
| [LIBRARY.md](LIBRARY.md) | Every public type and method |
| [EXPR.md](EXPR.md) | Per-op working precision in `expr!` |
| [PRECISION.md](PRECISION.md) | Retry budget vs exponent |
| [README.md](README.md) | Error-bound theory (contributors) |
