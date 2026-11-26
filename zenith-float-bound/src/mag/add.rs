//! `add_up`: the least `Mag` at or above the sum of two magnitudes.

use super::round::norm;
use super::round::settle;
use super::Mag;
use super::TOP;
use core::num::NonZeroU128;
use core::num::NonZeroU64;

impl Mag {
    /// Least `Mag` at or above `self + other` (CONTRACT 14.1, 14.3).
    ///
    /// An `Infinity` operand yields `Infinity` (CONTRACT 14.6). A zero operand yields the
    /// other value, before any alignment (CONTRACT 2.3); the result is built from that
    /// operand's integer parts by `Mag::from_canonical`, as for every other result. The
    /// operands are aligned to the larger exponent; bits shifted out become sticky, and a
    /// sticky sum rounds up one unit in the last place. When the exponent gap is at least 64
    /// and the smaller operand is nonzero, the result is the larger operand plus one unit.
    pub fn add_up(&self, other: &Mag) -> Mag {
        if self.is_infinity() || other.is_infinity() {
            Mag::infinity()
        } else {
            match (self.nonzero(), other.nonzero()) {
                (Some((ea, ma)), Some((eb, mb))) => add_parts(ea, ma, eb, mb),
                (Some((e, m)), None) | (None, Some((e, m))) => Mag::from_canonical(e, m.get()),
                (None, None) => Mag::zero(),
            }
        }
    }
}

fn add_parts(ea: i64, ma: NonZeroU64, eb: i64, mb: NonZeroU64) -> Mag {
    settle(add_checked(ea, ma, eb, mb), true)
}

fn add_checked(ea: i64, ma: NonZeroU64, eb: i64, mb: NonZeroU64) -> Option<Mag> {
    let (ea, ma, eb, mb) = if ea >= eb {
        (ea, ma, eb, mb)
    } else {
        (eb, mb, ea, ma)
    };
    // |ea - eb| < 2^65, so the checked subtraction cannot fail.
    let gap = i128::from(ea).checked_sub(i128::from(eb))?;
    match u32::try_from(gap) {
        // gap < 64: `ma << gap < 2^127`, so the aligned sum is exact in 128 bits (checked_add
        // cannot fail) and nonzero because `mb` is.
        Ok(gap) if gap < 64 => {
            let sum = NonZeroU128::from(mb).checked_add(u128::from(ma.get()) << gap)?;
            norm(sum, i128::from(eb).checked_sub(64)?, false, true)
        }
        // gap >= 64: the smaller operand is below one unit of the larger one's last place.
        _ => Some(successor(ea, ma)),
    }
}

/// The next `Mag` strictly above the canonical finite value `(exp, mant)`.
fn successor(exp: i64, mant: NonZeroU64) -> Mag {
    match mant.get().checked_add(1) {
        Some(next) => Mag::from_canonical(exp, next),
        // mant was u64::MAX: the next value is 2^63 at the next exponent, or Infinity (reachable).
        None => match exp.checked_add(1) {
            Some(exp) => Mag::from_canonical(exp, TOP),
            None => Mag::infinity(),
        },
    }
}
