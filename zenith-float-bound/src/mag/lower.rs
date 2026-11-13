//! Round-down products and differences.

use super::mul::mul_parts;
use super::round::norm;
use super::round::settle;
use super::Mag;
use core::num::NonZeroU128;
use core::num::NonZeroU64;

impl Mag {
    /// Greatest `Mag` at or below `self * other` (CONTRACT 14.2).
    ///
    /// An `Infinity` operand yields zero. A zero operand yields zero. Exponent
    /// underflow yields zero (CONTRACT 3.2). Exponent overflow yields the largest
    /// finite `Mag`.
    pub fn mul_lower(&self, other: &Mag) -> Mag {
        // `nonzero()` is `None` for `Infinity`, so an `Infinity` operand falls to zero.
        match (self.nonzero(), other.nonzero()) {
            (Some((ea, ma)), Some((eb, mb))) => mul_parts(ea, ma, eb, mb, false),
            _ => Mag::zero(),
        }
    }

    /// Greatest `Mag` at or below `self - other` (CONTRACT 14.2).
    ///
    /// The result is zero when `self <= other`, and when either operand is `Infinity`.
    /// Subtracting zero from a finite value gives that value (CONTRACT 2.3), built from its
    /// integer parts by `Mag::from_canonical`. A positive result whose exponent would fall
    /// below `i64::MIN` is zero.
    pub fn sub_lower(&self, other: &Mag) -> Mag {
        match (self.nonzero(), other.nonzero()) {
            (Some((ea, ma)), Some((eb, mb))) => sub_parts(ea, ma, eb, mb),
            (Some((e, m)), None) if !other.is_infinity() => Mag::from_canonical(e, m.get()),
            _ => Mag::zero(),
        }
    }
}

/// `a - b` rounded down, for canonical finite nonzero `a = (ea, ma)` and `b = (eb, mb)`.
fn sub_parts(ea: i64, ma: NonZeroU64, eb: i64, mb: NonZeroU64) -> Mag {
    settle(sub_checked(ea, ma, eb, mb), false)
}

fn sub_checked(ea: i64, ma: NonZeroU64, eb: i64, mb: NonZeroU64) -> Option<Mag> {
    // |ea - eb| < 2^65, so the checked subtraction cannot fail.
    let gap = i128::from(ea).checked_sub(i128::from(eb))?;
    match u32::try_from(gap) {
        Ok(gap) if gap < 64 => {
            // `ma << gap < 2^127`. `checked_sub` is `None` for a < b, and a zero difference is
            // a == b: both mean the exact difference is zero or less, so the answer is zero.
            let shifted = u128::from(ma.get()) << gap;
            match shifted
                .checked_sub(u128::from(mb.get()))
                .and_then(NonZeroU128::new)
            {
                Some(diff) => norm(diff, i128::from(eb).checked_sub(64)?, false, false),
                None => Some(Mag::zero()),
            }
        }
        // gap < 0: the minuend has the smaller exponent, so a < b and the difference is below zero.
        Err(_) if gap < 0 => Some(Mag::zero()),
        // gap >= 64.
        _ => sub_wide(ea, ma, gap, mb),
    }
}

/// Difference when `other` sits at least 64 exponents below `self`.
///
/// The smaller value is strictly below one ulp of `self`, so the round-down result is
/// `self` minus a fraction of that ulp. The subtraction is carried at 63 extra bits and
/// then truncated (CONTRACT 14.2). `gap >= 64` here.
fn sub_wide(ea: i64, ma: NonZeroU64, gap: i128, mb: NonZeroU64) -> Option<Mag> {
    let wide = u128::from(ma.get()) << 63;
    // Units of 2^(ea - 127): the subtrahend is mb / 2^(gap - 63), rounded up. shift >= 1 (gap >= 64).
    // A shift of 64 or more leaves only the fraction (mb < 2^64), which still rounds up to one unit:
    // for 64 <= shift < 128 the shifted arm would give `mb >> shift == 0` and, because the mask keeps
    // the low `shift >= 64` bits and mb is a nonzero value below 2^64, `frac == true`, i.e. exactly
    // the `(0, true)` of the other arm. So the guard `< 64` is the tightest exact one.
    let (subtrahend, frac) = match u32::try_from(gap.checked_sub(63)?) {
        Ok(shift) if shift < 64 => {
            let mb = u128::from(mb.get());
            (mb >> shift, (mb & !(u128::MAX << shift)) != 0)
        }
        _ => (0, true),
    };
    // wide >= 2^126 and subtrahend < 2^64, so both checked subtractions cannot fail.
    let diff = wide
        .checked_sub(subtrahend)?
        .checked_sub(u128::from(frac))?;
    let scale = i128::from(ea).checked_sub(127)?;
    // Range proof that `diff` is nonzero: ma >= 2^63 (canonical form), so wide >= 2^126. shift >= 1,
    // so subtrahend = mb >> shift < 2^63 (or 0), and frac <= 1. Hence diff >= 2^126 - 2^63 > 0 and
    // `NonZeroU128::new` is always `Some`; no zero-result arm is written. (Its `None` could only
    // reach `settle`, whose fallback is the sound bound.)
    norm(NonZeroU128::new(diff)?, scale, false, false)
}
