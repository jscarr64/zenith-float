//! The precision type and its limits (API-freeze unit 5; Paul and Rusty, 2026-10-23).
//!
//! zf's precision type is `usize` (`zenith_float::Precision` after Rusty's i64 work). The current pin (fe5c15d)
//! exports neither the alias nor the limits, so zb mirrors them here.
//! TODO(SP-i64): at the i64 re-pin, re-export `zenith_float::Precision`, `MIN_PRECISION` and `MAX_PRECISION`
//! instead of defining them.

/// A precision in bits: zf's precision type.
pub type Precision = usize;

/// Smallest accepted precision in bits.
pub const MIN_PRECISION: Precision = 1;

/// Largest accepted precision in bits: `min(2^60, usize::MAX / 4) - 1` (zf w6 defs.rs), so an exponent plus or minus
/// a precision stays inside `±2^62`. That is `2^60 - 1` on 64-bit targets and `2^30 - 1` on 32-bit targets.
pub const MAX_PRECISION: Precision = max_precision();

const fn max_precision() -> Precision {
    let by_size = usize::MAX >> 2;
    let lim = if let Some(by_exponent) = 1usize.checked_shl(60) {
        if by_exponent < by_size {
            by_exponent
        } else {
            by_size
        }
    } else {
        by_size
    };
    if let Some(v) = lim.checked_sub(1) {
        v
    } else {
        0
    }
}
