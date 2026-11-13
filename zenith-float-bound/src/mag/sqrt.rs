//! `sqrt_up` and `sqrt_lower` by an integer square root, exponent parity first.

use super::round::loose;
use super::round::norm;
use super::round::settle;
use super::Mag;
use core::num::NonZeroU128;
use core::num::NonZeroU64;

/// `2^64` as a nonzero `u128`. `isolate_highest_one` keeps only the top set bit of a value whose
/// top bit is bit 64, so there is no `None` arm to settle and no value to get wrong.
fn p64() -> NonZeroU128 {
    (NonZeroU128::MIN | (1u128 << 64)).isolate_highest_one()
}

/// `2^63` as a nonzero `u128`, built the same way.
fn p63() -> NonZeroU128 {
    (NonZeroU128::MIN | (1u128 << 63)).isolate_highest_one()
}

impl Mag {
    /// Least `Mag` at or above `sqrt(self)` (CONTRACT 14.1, 14.4).
    ///
    /// `sqrt(Infinity)` is `Infinity`. `sqrt(0)` is zero. The exponent parity is fixed
    /// first so the radicand is a `u128`, then the integer square root is rounded up
    /// when the square is not exact.
    pub fn sqrt_up(&self) -> Mag {
        sqrt(self, true)
    }

    /// Greatest `Mag` at or below `sqrt(self)` (CONTRACT 14.2, 14.4).
    ///
    /// `sqrt(Infinity)` and `sqrt(0)` are zero. A nonzero finite value uses the same
    /// integer square root as [`Mag::sqrt_up`], truncated toward zero.
    pub fn sqrt_lower(&self) -> Mag {
        sqrt(self, false)
    }
}

fn sqrt(mag: &Mag, up: bool) -> Mag {
    if mag.is_infinity() {
        return loose(up);
    }
    match mag.nonzero() {
        Some((exp, mant)) => settle(sqrt_checked(exp, mant, up), up),
        None => Mag::zero(),
    }
}

fn sqrt_checked(exp: i64, mant: NonZeroU64, up: bool) -> Option<Mag> {
    let exp = i128::from(exp);
    // value = mant * 2^(exp - 64). Even exp: N = mant * 2^64 with scale exp - 128. Odd exp:
    // N = mant * 2^63 with scale exp - 127. Either way `scale` is even, so `scale / 2` is exact.
    // |exp| < 2^63, so the checked steps cannot fail.
    let (factor, scale) = if exp.checked_rem_euclid(2)? == 0 {
        (p64(), exp.checked_sub(128)?)
    } else {
        (p63(), exp.checked_sub(127)?)
    };
    // mant < 2^64 and factor <= 2^64, so the product fits in a u128 and is nonzero by type.
    let radicand = NonZeroU128::from(mant).checked_mul(factor)?;
    let root = radicand.isqrt();
    // root < 2^64, so root * root < 2^128 and the checked multiply cannot fail.
    let sticky = root.get().checked_mul(root.get())? != radicand.get();
    norm(root, scale.checked_div(2)?, sticky, up)
}
