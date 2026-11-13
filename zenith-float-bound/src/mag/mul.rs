//! `mul_up` and the shared product used by `mul_lower`.

use super::round::norm;
use super::round::settle;
use super::Mag;
use core::num::NonZeroU128;
use core::num::NonZeroU64;

impl Mag {
    /// Least `Mag` at or above `self * other` (CONTRACT 1.3, 14.1).
    ///
    /// An `Infinity` operand yields `Infinity`, including zero times `Infinity`
    /// (CONTRACT 2.4, 14.6). A zero operand with a finite partner yields zero
    /// (CONTRACT 2.3). The 128-bit product is normalised to 64 bits; a nonzero
    /// residue rounds up. Exponent underflow rounds up to the smallest positive
    /// `Mag` (CONTRACT 3.2). Exponent overflow is `Infinity` (CONTRACT 3.3).
    pub fn mul_up(&self, other: &Mag) -> Mag {
        if self.is_infinity() || other.is_infinity() {
            Mag::infinity()
        } else {
            match (self.nonzero(), other.nonzero()) {
                (Some((ea, ma)), Some((eb, mb))) => mul_parts(ea, ma, eb, mb, true),
                _ => Mag::zero(),
            }
        }
    }
}

/// Product of two canonical finite nonzero values.
pub(super) fn mul_parts(ea: i64, ma: NonZeroU64, eb: i64, mb: NonZeroU64, up: bool) -> Mag {
    settle(mul_checked(ea, ma, eb, mb, up), up)
}

fn mul_checked(ea: i64, ma: NonZeroU64, eb: i64, mb: NonZeroU64, up: bool) -> Option<Mag> {
    // `ma * mb < 2^128` and both factors are nonzero, so the product is nonzero and the checked
    // multiply cannot fail; |ea + eb - 128| < 2^66, so the exponent steps cannot fail either.
    let prod = NonZeroU128::from(ma).checked_mul(NonZeroU128::from(mb))?;
    let scale = i128::from(ea)
        .checked_add(i128::from(eb))?
        .checked_sub(128)?;
    norm(prod, scale, false, up)
}
