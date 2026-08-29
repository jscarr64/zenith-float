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
use zenith_float::{Consts, ExactNum, RoundingMode};

let p = 256;
let rm = RoundingMode::ToEven;
let mut cc = Consts::new().expect("constants cache");

let a = ExactNum::from_u32(3, p);
let b = ExactNum::parse("0.5", zenith_float::Radix::Dec, p, rm, &mut cc);
let c = a.add(&b, p, rm); // 3.5
```

Compile-time decimals:

```rust
use zenith_float::exact;

let x = exact!("1.25");
```

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

With `std`:

```rust
#[cfg(feature = "std")]
println!("{}", y);           // decimal Display
#[cfg(feature = "std")]
println!("{:e}", y);         // scientific
```

Any radix 2–36:

```rust
let s = y.format(zenith_float::Radix::Hex, RoundingMode::ToEven, ctx.consts())
    .expect("format");
```

## 5. When to call methods instead of `expr!`

Use `ExactNum` methods when you need paired results (`sin_cos`, `sinh_cosh`), IEEE-style split (`frexp`, `ilogb`), or compensated primitives (`two_sum`, `two_product`, `fused_sum`, `fused_dot`, `polyval`).

For **complex** expressions use `cexpr!` (same working-precision loop as `expr!`, cancellation on both parts). The imaginary unit in the expression is `I`. `expr!` stays real-valued. This crate already ships macros (`expr!`, `cexpr!`, `exact!`, `fbig!`); the “no macros” rule is Accumath, not zenith-float.

Hardware binary interchange (plotting, FFI) is not in this crate. A small downstream crate can call `frexp` / `ilogb` and pack bits.

## 6. Next documents

| Doc | Audience |
| --- | --- |
| [LIBRARY.md](LIBRARY.md) | Every public type and method |
| [EXPR.md](EXPR.md) | Per-op working precision in `expr!` |
| [PRECISION.md](PRECISION.md) | Retry budget vs exponent |
| [BUILD_CHECKLIST.md](BUILD_CHECKLIST.md) | Implemented vs out of scope |
| [README.md](README.md) | Error-bound theory (contributors) |
