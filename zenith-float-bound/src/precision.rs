//! The precision type and its limits (API-freeze unit 5; Paul and Rusty, 2026-10-23).
//!
//! These are zenith-float's own items, re-exported, so there is one definition (zenith-float-num `defs.rs`):
//! `Precision` is `usize`, `MIN_PRECISION` is 1 and `MAX_PRECISION` is the fixed `2^60 - 1`, so an exponent plus or
//! minus a precision stays inside `±2^62`.

pub use zenith_float_num::{Precision, MAX_PRECISION, MIN_PRECISION};
