//! The precision type and its limits (API-freeze unit 5; Paul and Rusty, 2026-10-23).
//!
//! zb reads these from zenith-float-num (ZB-TD02, ZF-M32): there is one definition, in
//! zenith-float-num defs.rs, and zb never duplicates the formula.

/// A precision in bits: zenith-float-num precision type.
pub type Precision = zenith_float_num::Precision;

/// Smallest accepted precision in bits (zenith-float-num MIN_PRECISION).
pub const MIN_PRECISION: Precision = zenith_float_num::MIN_PRECISION;

/// Largest accepted precision in bits: zenith-float-num MAX_PRECISION, read, not re-derived.
/// ZF-M32 sets it to 2^32 - 1 on every supported (64-bit) target.
pub const MAX_PRECISION: Precision = zenith_float_num::MAX_PRECISION;
