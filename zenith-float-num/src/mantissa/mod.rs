// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! Everything related to mantissa.

mod base;
mod cbrt;
mod conv;
mod div;
mod fft;
mod mul;
mod sqrt;
mod toom2;
mod toom3;
mod util;

pub use base::AddSubSpec;
pub use base::Mantissa;
#[cfg(test)]
pub(crate) use base::CLONE_CALLS;
