//! Integration tests.

#[cfg(test)]
#[cfg(feature = "mpfr-tests")]
#[cfg(all(target_arch = "x86_64", target_os = "linux"))]
mod mpfr;
