

zenith-float is an arbitrary precision floating-point numbers library designed for performance, portability, and implemented purely in Rust.

All arithmetic is software big-float on integer words. The library does not use hardware `f32` or `f64` for calculations. Construct numbers from integers or decimal/binary strings.

The library implements the basic operations and functions. It uses classical algorithms such as Karatsuba, Toom-Cook, Schönhage-Strassen algorithm, and others.

The library can work without the standard library provided there is a memory allocator.




Calculate Pi with 1024 bit precision rounded to the nearest even number.

``` rust
use zenith_float::Consts;
use zenith_float::RoundingMode;
use zenith_float::ctx::Context;
use zenith_float::expr;

// Create a context with precision 1024, and rounding to even.
let mut ctx = Context::new(1024, RoundingMode::ToEven, 
    Consts::new().expect("Constants cache initialized"),
    -10000, 10000);

// Compute pi: pi = 6*arctan(1/sqrt(3))
let pi = expr!(6 * atan(1 / sqrt(3)), &mut ctx);

// Use library's constant value for verifying the result.
let pi_lib = ctx.const_pi();

// Compare computed constant with library's constant
assert_eq!(pi.cmp(&pi_lib), Some(0));
```
