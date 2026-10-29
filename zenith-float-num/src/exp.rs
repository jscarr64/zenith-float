//! Overflow-free exponent arithmetic.
//!
//! Exponents lie in `EXPONENT_MIN..=EXPONENT_MAX` (±(2^62 − 1)) and precisions are at most
//! `MAX_PRECISION` (< 2^60), so a sum of two exponents fits `i64` but a sum of two exponents and a
//! shift, or a difference of extreme exponents, may not. Every such expression is evaluated here in
//! `i128` and narrowed with [`narrow`], which saturates to `±WIDE_LIMIT`: far outside the exponent
//! range (so overflow and underflow are still detected), yet leaving 2^61 of `i64` headroom for the
//! `± precision` adjustments the subnormal paths make afterwards.

use crate::defs::{Exponent, EXPONENT_MAX};

/// Saturation bound of [`narrow`]: `EXPONENT_MAX + 2^61`.
pub(crate) const WIDE_LIMIT: i64 = EXPONENT_MAX + (1 << 61);

const _: () = assert!(WIDE_LIMIT.checked_add(1 << 61).is_some());

/// Widens an exponent-like value for exact arithmetic.
#[inline]
pub(crate) fn wide(v: i64) -> i128 {
    i128::from(v)
}

/// Widens a bit count or shift.
#[inline]
pub(crate) fn wide_u(v: usize) -> i128 {
    // usize is at most 64 bits on supported targets, so this is lossless.
    i128::try_from(v).unwrap_or(i128::MAX)
}

/// Narrows an exact exponent to `i64`, saturating to `±WIDE_LIMIT`.
#[inline]
pub(crate) fn narrow(v: i128) -> i64 {
    let lim = i128::from(WIDE_LIMIT);
    let c = v.clamp(-lim, lim);
    i64::try_from(c).unwrap_or(if c < 0 { -WIDE_LIMIT } else { WIDE_LIMIT })
}

/// `|a − b|` as a shift count, saturating at `usize::MAX` (a shift that large discards the operand).
#[inline]
pub(crate) fn abs_diff_shift(a: i64, b: i64) -> usize {
    usize::try_from(a.abs_diff(b)).unwrap_or(usize::MAX)
}

/// `e · k` for an exponent and an integer multiplier (`pow`), narrowed with [`narrow`].
#[inline]
pub(crate) fn mul(e: Exponent, k: i64) -> i64 {
    narrow(wide(e).saturating_mul(wide(k)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::defs::EXPONENT_MIN;

    #[test]
    fn narrow_saturates_outside_the_range() {
        assert_eq!(narrow(wide(EXPONENT_MAX) * 4), WIDE_LIMIT);
        assert_eq!(narrow(wide(EXPONENT_MIN) * 4), -WIDE_LIMIT);
        assert_eq!(narrow(wide(EXPONENT_MIN) - 1), EXPONENT_MIN - 1);
        assert_eq!(mul(EXPONENT_MAX, -3), -WIDE_LIMIT);
        assert_eq!(
            abs_diff_shift(EXPONENT_MAX, EXPONENT_MIN),
            // 2 * EXPONENT_MAX in u64, saturated to `usize` like `abs_diff_shift` (32-bit: `usize::MAX`).
            usize::try_from(EXPONENT_MAX.unsigned_abs().saturating_mul(2)).unwrap_or(usize::MAX)
        );
    }
}
