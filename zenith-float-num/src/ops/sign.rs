//! Sign manipulation and successor operations.

use crate::common::util::{round_p, sub_borrow};
use crate::defs::{Error, Word, EXPONENT_MAX, EXPONENT_MIN, WORD_MAX, WORD_SIGNIFICANT_BIT};
use crate::num::ExactNumNumber;
use crate::{RoundingMode, Sign};

impl ExactNumNumber {
    /// One ulp up in magnitude. `ExponentOverflow` past the largest finite magnitude.
    fn bump_up(x: &mut Self) -> Result<(), Error> {
        let was_subnormal = x.is_subnormal();
        if x.mantissa_mut().add_ulp() {
            if x.exponent() >= EXPONENT_MAX {
                return Err(Error::ExponentOverflow(x.sign()));
            }
            let ml = x.mantissa_mut().len() - 1;
            x.mantissa_mut().digits_mut()[ml] = WORD_SIGNIFICANT_BIT;
            x.set_exponent(x.exponent() + 1);
        }
        if was_subnormal {
            // The step can carry into the leading bit (largest subnormal -> least normal).
            x.mantissa_mut().update_bit_len();
        }
        Ok(())
    }

    /// One ulp down in magnitude, on the p-bit grid of `x`'s exponent; a nonzero `x` steps to
    /// zero (same sign) only from the least subnormal.
    fn bump_down(x: &mut Self, p: usize) -> Result<(), Error> {
        let mut borrow: Word = 1;
        for v in x.mantissa_mut().digits_mut().iter_mut() {
            borrow = sub_borrow(*v, 0, borrow, v);
            if borrow == 0 {
                break;
            }
        }
        let ml = x.mantissa_mut().len() - 1;
        if x.mantissa_mut().digits_mut()[ml] & WORD_SIGNIFICANT_BIT != 0 {
            return Ok(());
        }
        if x.exponent() > EXPONENT_MIN && !x.is_subnormal() {
            // From a power of two 2^(e-1) the predecessor is 0.11…1 · 2^(e-1): all p bits set
            // one binade lower, not the unnormalized 0.011…1 · 2^e.
            x.mantissa_mut().digits_mut().fill(WORD_MAX);
            x.set_exponent(x.exponent() - 1);
            return Ok(());
        }
        // On the subnormal grid (exponent field EXPONENT_MIN): fewer significant bits, or zero.
        if x.mantissa_mut().is_all_zero() {
            let s = x.sign();
            *x = Self::new(p)?;
            x.set_sign(s);
        } else {
            x.mantissa_mut().update_bit_len();
        }
        Ok(())
    }

    /// Returns a value with the magnitude of `self` and the sign of `sign`.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory.
    ///  - InvalidArgument: the precision is incorrect.
    pub fn copysign(&self, sign: &Self, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;

        let mut ret = self.clone()?;
        if sign.is_negative() {
            ret.set_sign(Sign::Neg);
        } else if sign.is_zero() {
            ret.set_sign(sign.sign());
        } else {
            ret.set_sign(Sign::Pos);
        }
        ret.set_precision(p, rm)?;

        ret.set_inexact(ret.inexact() | self.inexact() | sign.inexact());
        Ok(ret)
    }

    /// Returns the next representable value from `self` toward `toward` at precision `p`.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory.
    ///  - InvalidArgument: the precision is incorrect.
    pub fn next_after(&self, toward: &Self, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        // Compare the exact values: rounding both to p first could make distinct values equal.
        match self.cmp(toward) {
            0 => {
                let mut x = self.clone()?;
                x.set_precision(p, rm)?;
                Ok(x)
            }
            c => self.next_toward(c < 0, p),
        }
    }

    /// The next p-bit value after `self` going up (`upward`) or down. A `self` that is not
    /// representable at `p` goes to the nearest p-bit value strictly beyond it in that
    /// direction; past the largest finite magnitude the result is `ExponentOverflow` (±Inf).
    pub(crate) fn next_toward(&self, upward: bool, p: usize) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        let mut x = self.clone()?;
        x.set_inexact(false);
        x.set_precision(p, if upward { RoundingMode::Up } else { RoundingMode::Down })?;
        if x.inexact() || x.cmp(self) != 0 {
            // Rounding in the step direction already moved strictly past `self`.
            x.set_inexact(true);
            return Ok(x);
        }
        if x.is_zero() {
            let mut v = Self::min_positive(p)?;
            if !upward {
                v.set_sign(Sign::Neg);
            }
            v.set_inexact(true);
            return Ok(v);
        }
        // Up for a negative value (or down for a positive one) shrinks the magnitude.
        if upward == x.is_positive() {
            Self::bump_up(&mut x)?;
        } else {
            Self::bump_down(&mut x, p)?;
        }
        x.set_inexact(true);
        Ok(x)
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::Consts;

    #[test]
    fn test_copysign() -> Result<(), crate::Error> {
        let p = 128;
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;

        let pos = ExactNumNumber::from_word(5, p)?;
        let neg = ExactNumNumber::from_i8(-3, p)?;

        let r = pos.copysign(&neg, p, rm)?;
        assert!(r.is_negative());
        assert!(r.cmp(&ExactNumNumber::from_word(5, p)?.neg()?) == 0);

        let r = neg.copysign(&pos, p, rm)?;
        assert!(r.is_positive());
        assert!(r.cmp(&ExactNumNumber::from_word(3, p)?) == 0);

        let zero = ExactNumNumber::new(p)?;
        let neg_zero = zero.copysign(&neg, p, rm)?;
        let pos_zero = zero.copysign(&pos, p, rm)?;
        assert!(neg_zero.is_zero() && neg_zero.is_negative());
        assert!(pos_zero.is_zero() && pos_zero.is_positive());

        let x = ExactNumNumber::parse("1.8", crate::Radix::Hex, p, RoundingMode::None, &mut cc)?;
        let y = ExactNumNumber::parse("-1", crate::Radix::Hex, p, RoundingMode::None, &mut cc)?;
        let r = x.copysign(&y, p, rm)?;
        assert!(r.is_negative());
        Ok(())
    }

    #[test]
    fn test_next_after() -> Result<(), crate::Error> {
        let p = 128;
        let rm = RoundingMode::ToEven;

        let one = ExactNumNumber::from_word(1, p)?;
        let two = ExactNumNumber::from_word(2, p)?;
        let next = one.next_after(&two, p, rm)?;
        assert!(next.cmp(&one) > 0);
        assert!(next.cmp(&two) < 0);

        let prev = two.next_after(&one, p, rm)?;
        assert!(prev.cmp(&two) < 0);
        assert!(prev.cmp(&one) > 0);

        let zero = ExactNumNumber::new(p)?;
        let min_pos = ExactNumNumber::min_positive(p)?;
        assert!(zero.next_after(&two, p, rm)?.cmp(&min_pos) == 0);

        let neg_min = min_pos.neg()?;
        assert!(zero.next_after(&neg_min, p, rm)?.cmp(&neg_min) == 0);

        assert!(one.next_after(&one, p, rm)?.cmp(&one) == 0);
        Ok(())
    }
}
