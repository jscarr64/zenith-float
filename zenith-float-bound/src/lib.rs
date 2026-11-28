//! zenith-float-bound: guaranteed-bounds (ball arithmetic) on zenith-float.
//!
//! Pre-release. Nothing in this crate is guaranteed yet.
#![forbid(unsafe_code)]
#![deny(
    missing_docs,
    clippy::allow_attributes,
    clippy::arithmetic_side_effects
)]

// P6 (Jeff, 2026-10-30): 64-bit targets only.
#[cfg(not(target_pointer_width = "64"))]
compile_error!("zenith-float-bound requires a 64-bit target (target_pointer_width = 64); 32-bit targets are not supported.");

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
pub mod cball; // CB
pub mod lu; // LA2
pub mod matrix; // S6
pub mod qr; // LA3

// Held slices (HOLD.md at the repo root): not compiled; the mod line stays here, commented.
// HOLD H05: pub mod carlson; // S4b
// HOLD H06: pub mod cholesky; // LA1
// HOLD H08: pub mod eigen; // LA5
// HOLD H01: pub mod elementary; // S3
// HOLD H07: pub mod gls; // LA4
// HOLD H02: pub mod hypergeom; // S4a
// HOLD H04: pub mod quadrature; // S5
// HOLD H03: pub mod zf_special; // S4c

pub use precision::{Precision, MAX_PRECISION, MIN_PRECISION};
