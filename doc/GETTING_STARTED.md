# Getting started with zenith-float

This is a short path from a new crate dependency to a computed, formatted result. For the full public surface, see [LIBRARY.md](LIBRARY.md). For how `expr!` rounds, see [EXPR.md](EXPR.md).

All arithmetic is software integer limbs. There are no hardware binary interchange types in this crate.

## 1. Depend on the crate

```toml
[dependencies]
zenith-float = "0.1.0"
```

`std` is on by default (formatting and `FromStr`). For `no_std` plus an allocator:

```toml
zenith-float = { version = "0.1.0", default-features = false }
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

Hardware binary interchange (plotting, FFI) is not in this crate. A small downstream crate can call `frexp` / `ilogb` and pack bits.

## 6. Next documents

| Doc | Audience |
| --- | --- |
| [LIBRARY.md](LIBRARY.md) | Every public type and method |
| [EXPR.md](EXPR.md) | Per-op working precision in `expr!` |
| [PRECISION.md](PRECISION.md) | Retry budget vs exponent |
| [BUILD_CHECKLIST.md](BUILD_CHECKLIST.md) | Implemented vs out of scope |
| [README.md](README.md) | Error-bound theory (contributors) |
