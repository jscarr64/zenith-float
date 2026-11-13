//! Unit in the last place, half an ulp, and the half-ulp predicate.
//!
//! `half_ulp` is a sound rounding-error term only for nearest modes. Directed modes
//! need a full ulp from [`Mag::ulp_mag`].

use super::Mag;
use super::TOP;
use crate::precision::Precision;

impl Mag {
    /// `2^(pre_round_exp - prec)`, placed as a canonical power of two (CONTRACT 14.7).
    ///
    /// The exponent field is computed in `i128` as `pre_round_exp - prec + 1` and then
    /// placed. `prec` is the full `Precision` value. An exponent field above `i64::MAX`
    /// (only `ulp_mag(i64::MAX, 0)`) is `Infinity`. An exponent field below `i64::MIN`
    /// rounds up to the smallest positive `Mag` (CONTRACT 3.2).
    pub fn ulp_mag(pre_round_exp: i64, prec: Precision) -> Mag {
        // |pre_round_exp - prec + 1| < 2^65, so the checked steps cannot fail (a `usize` prec always fits `i128`). `Infinity` is the
        // one fallback: an error term must never be too small.
        match i128::try_from(prec)
            .ok()
            .and_then(|p| i128::from(pre_round_exp).checked_sub(p))
            .and_then(|diff| diff.checked_add(1))
        {
            Some(mag_exp) => place_pow2(mag_exp),
            None => Mag::infinity(),
        }
    }

    /// `2^(mid_exp - prec - 1)`, placed as a canonical power of two (CONTRACT 14.7).
    ///
    /// At `prec == 0` this is `2^(mid_exp - 1)`. The exponent field is `mid_exp - prec`,
    /// computed in `i128` over the full `Precision` range. Underflow rounds up to the
    /// smallest positive `Mag`; overflow is `Infinity`.
    ///
    /// This value is the right midpoint-rounding term only under a nearest mode. A
    /// directed mode needs the full ulp from [`Mag::ulp_mag`].
    pub fn half_ulp(mid_exp: i64, prec: Precision) -> Mag {
        // |mid_exp - prec| < 2^64, so the checked step cannot fail. `Infinity` is the one fallback.
        match i128::try_from(prec)
            .ok()
            .and_then(|p| i128::from(mid_exp).checked_sub(p))
        {
            Some(mag_exp) => place_pow2(mag_exp),
            None => Mag::infinity(),
        }
    }

    /// Whether `self` is strictly below half an ulp of a midpoint at `mid_exp` and `prec`.
    ///
    /// For a finite nonzero `Mag` this is `exp <= mid_exp - prec - 1`, evaluated in
    /// `i128` so a threshold below `i64::MIN` is not rounded up (CONTRACT 1.4, 14.7).
    /// The canonical zero is below every half ulp. `Infinity` is not.
    pub fn smaller_than_half_ulp(&self, mid_exp: i64, prec: Precision) -> bool {
        // |mid_exp - prec - 1| < 2^65, so the checked steps cannot fail; `false` ("not known to be
        // below") is the one fallback.
        let Some(threshold) = i128::try_from(prec)
            .ok()
            .and_then(|p| i128::from(mid_exp).checked_sub(p))
            .and_then(|diff| diff.checked_sub(1))
        else {
            return false;
        };
        match self.nonzero() {
            Some((exp, _)) => i128::from(exp) <= threshold,
            None => !self.is_infinity(),
        }
    }
}

/// `2^(mag_exp - 1)` as `Finite { exp: mag_exp, mant: 2^63 }`, or the range clamp: above
/// `i64::MAX` is `Infinity`, below `i64::MIN` is `min_pos` (round-up underflow, CONTRACT 3.2).
fn place_pow2(mag_exp: i128) -> Mag {
    match i64::try_from(mag_exp) {
        Ok(exp) => Mag::from_canonical(exp, TOP),
        Err(_) => {
            if mag_exp > i128::from(i64::MAX) {
                Mag::infinity()
            } else {
                Mag::min_pos()
            }
        }
    }
}
