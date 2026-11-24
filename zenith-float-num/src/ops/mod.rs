// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! High-level operations on the numbers.

mod acos;
mod acosh;
mod airy;
mod asin;
mod asinh;
mod atan;
mod atan2;
mod atanh;
mod cbrt;
pub mod consts;
mod cos;
mod cosh;
pub(crate) mod exp_bb;
mod expm1;
mod fma;
mod hypot;
mod ieee;
mod jacobi;
mod log;
mod log1p;
mod log_newton;
pub(crate) mod log_primes;
mod nroot;
mod pow;
mod recip;
mod series;
mod sign;
mod sin;
mod sin_cos;
mod sinh;
pub(crate) mod special;
mod sqrt;
mod tan;
mod tanh;
mod util;
pub(crate) mod zeta_ui;

#[cfg(test)]
mod tests;
