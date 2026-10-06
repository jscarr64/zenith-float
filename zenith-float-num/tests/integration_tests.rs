// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! Integration tests.

#[cfg(test)]
#[cfg(feature = "mpfr-tests")]
#[cfg(all(target_arch = "x86_64", target_os = "linux"))]
mod mpfr;
