//! zenith-bound: guaranteed-bounds (ball arithmetic) on zenith-float.
//!
//! Pre-release. Nothing in this crate is guaranteed yet.
#![forbid(unsafe_code)]
#![deny(
    missing_docs,
    clippy::allow_attributes,
    clippy::arithmetic_side_effects
)]

pub mod array;
pub mod ball;
pub mod enclosure;
pub mod error;
pub mod mag;
pub mod policy;
pub mod precision;
pub mod tri;

// One module per parallel slice (S0, 2026-10-24). Each file belongs to its slice alone, so no
// two slices edit the same file; this list is final.
pub mod arith; // S1
pub mod bridge; // S2
pub mod carlson; // S4b
pub mod cball; // CB
pub mod cholesky; // LA1
pub mod eigen; // LA5
pub mod elementary; // S3
pub mod gls; // LA4
pub mod hypergeom; // S4a
pub mod lu; // LA2
pub mod matrix; // S6
pub mod qr; // LA3
pub mod quadrature; // S5
pub mod zf_special; // S4c

pub use precision::{Precision, MAX_PRECISION, MIN_PRECISION};
