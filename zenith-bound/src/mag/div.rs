//! `div_up` and `div_lower` on a 128/64 division with a sticky remainder.

use super::round::loose;
use super::round::norm;
use super::round::settle;
use super::Mag;
use core::num::NonZeroU128;
use core::num::NonZeroU64;

impl Mag {
    /// Least `Mag` at or above `self / other` (CONTRACT 14.1, 14.4).
    ///
    /// `Infinity` is checked before a zero divisor, and both are checked before the
    /// zero-numerator shortcut, so `0 / 0` is `Infinity`. Any `Infinity` operand or a
    /// zero divisor yields `Infinity`. A zero numerator with a finite nonzero divisor
    /// yields zero. Exponent underflow rounds up to the smallest positive `Mag`;
    /// exponent overflow is `Infinity`.
    pub fn div_up(&self, other: &Mag) -> Mag {
        div_parts(self, other, true)
    }

    /// Greatest `Mag` at or below `self / other` (CONTRACT 14.2, 14.4).
    ///
    /// An `Infinity` operand or a zero divisor yields zero. A zero numerator yields
    /// zero. Exponent underflow yields zero. Exponent overflow yields the largest
    /// finite `Mag`.
    pub fn div_lower(&self, other: &Mag) -> Mag {
        div_parts(self, other, false)
    }
}

fn div_parts(num: &Mag, den: &Mag, up: bool) -> Mag {
    if num.is_infinity() || den.is_infinity() {
        return loose(up);
    }
    match (num.nonzero(), den.nonzero()) {
        // A zero divisor is decided before a zero numerator (CONTRACT 14.6): 0 / 0 is Infinity up.
        (_, None) => loose(up),
        (None, Some(_)) => Mag::zero(),
        (Some((ea, ma)), Some((eb, mb))) => settle(div_checked(ea, ma, eb, mb, up), up),
    }
}

fn div_checked(ea: i64, ma: NonZeroU64, eb: i64, mb: NonZeroU64, up: bool) -> Option<Mag> {
    let numerator = u128::from(ma.get()) << 64;
    let divisor = NonZeroU128::from(mb);
    // quot >= 2^127 / 2^64 = 2^63, so it is nonzero; division does not carry that in its type,
    // so a zero quotient would be `None` here and settle to the sound bound.
    let quot = NonZeroU128::new(numerator / divisor)?;
    let sticky = numerator % divisor != 0;
    // |ea - eb - 64| < 2^66, so the checked steps cannot fail.
    let scale = i128::from(ea)
        .checked_sub(i128::from(eb))?
        .checked_sub(64)?;
    norm(quot, scale, sticky, up)
}
