//! Sign manipulation and successor operations.

use crate::common::util::{round_p, sub_borrow};
use crate::defs::{Error, Word, WORD_SIGNIFICANT_BIT, WORD_MAX};
use crate::num::ExactNumNumber;
use crate::{RoundingMode, Sign};

impl ExactNumNumber {
    fn bump_up(x: &mut Self) {
        if x.mantissa_mut().add_ulp() {
            let ml = x.mantissa_mut().len() - 1;
            x.mantissa_mut().digits_mut()[ml] = WORD_SIGNIFICANT_BIT;
            x.set_exponent(x.exponent() + 1);
        }
    }

    fn bump_down(x: &mut Self) {
        let mut borrow: Word = 1;
        for v in x.mantissa_mut().digits_mut().iter_mut() {
            borrow = sub_borrow(*v, 0, borrow, v);
            if borrow == 0 {
                return;
            }
        }
        let ml = x.mantissa_mut().len() - 1;
        x.mantissa_mut().digits_mut().fill(WORD_MAX);
        x.mantissa_mut().digits_mut()[ml] = WORD_SIGNIFICANT_BIT;
        x.set_exponent(x.exponent() - 1);
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

        let mut ret = self.abs()?;
        ret.set_precision(p, rm)?;

        if sign.is_negative() {
            ret.set_sign(Sign::Neg);
        } else if sign.is_zero() {
            ret.set_sign(sign.sign());
        } else {
            ret.set_sign(Sign::Pos);
        }

        ret.set_inexact(ret.inexact() | self.inexact() | sign.inexact());
        Ok(ret)
    }

    /// Returns the next representable value from `self` toward `toward` at precision `p`.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory.
    ///  - InvalidArgument: the precision is incorrect.
    pub fn next_after(
        &self,
        toward: &Self,
        p: usize,
        _rm: RoundingMode,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;

        let mut x = self.clone()?;
        x.set_precision(p, RoundingMode::None)?;

        let mut t = toward.clone()?;
        t.set_precision(p, RoundingMode::None)?;

        if x.cmp(&t) == 0 {
            return Ok(x);
        }

        let upward = x.cmp(&t) < 0;

        if x.is_zero() {
            let mut v = Self::min_positive(p)?;
            if !upward {
                v.set_sign(Sign::Neg);
            }
            return Ok(v);
        }

        if upward {
            Self::bump_up(&mut x);
        } else {
            Self::bump_down(&mut x);
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
    fn test_copysign() {
        let p = 128;
        let mut cc = Consts::new().unwrap();
        let rm = RoundingMode::ToEven;

        let pos = ExactNumNumber::from_word(5, p).unwrap();
        let neg = ExactNumNumber::from_i8(-3, p).unwrap();

        let r = pos.copysign(&neg, p, rm).unwrap();
        assert!(r.is_negative());
        assert!(r.cmp(&ExactNumNumber::from_word(5, p).unwrap().neg().unwrap()) == 0);

        let r = neg.copysign(&pos, p, rm).unwrap();
        assert!(r.is_positive());
        assert!(r.cmp(&ExactNumNumber::from_word(3, p).unwrap()) == 0);

        let zero = ExactNumNumber::new(p).unwrap();
        let neg_zero = zero.copysign(&neg, p, rm).unwrap();
        let pos_zero = zero.copysign(&pos, p, rm).unwrap();
        assert!(neg_zero.is_zero() && neg_zero.is_negative());
        assert!(pos_zero.is_zero() && pos_zero.is_positive());

        let x = ExactNumNumber::parse("1.8p+0", crate::Radix::Hex, p, RoundingMode::None, &mut cc)
            .unwrap();
        let y = ExactNumNumber::parse("-1p+0", crate::Radix::Hex, p, RoundingMode::None, &mut cc)
            .unwrap();
        let r = x.copysign(&y, p, rm).unwrap();
        assert!(r.is_negative());
    }

    #[test]
    fn test_next_after() {
        let p = 128;
        let rm = RoundingMode::ToEven;

        let one = ExactNumNumber::from_word(1, p).unwrap();
        let two = ExactNumNumber::from_word(2, p).unwrap();
        let next = one.next_after(&two, p, rm).unwrap();
        assert!(next.cmp(&one) > 0);
        assert!(next.cmp(&two) < 0);

        let prev = two.next_after(&one, p, rm).unwrap();
        assert!(prev.cmp(&two) < 0);
        assert!(prev.cmp(&one) > 0);

        let zero = ExactNumNumber::new(p).unwrap();
        let min_pos = ExactNumNumber::min_positive(p).unwrap();
        assert!(
            zero.next_after(&two, p, rm)
                .unwrap()
                .cmp(&min_pos)
                == 0
        );

        let neg_min = min_pos.neg().unwrap();
        assert!(
            zero.next_after(&neg_min, p, rm)
                .unwrap()
                .cmp(&neg_min)
                == 0
        );

        assert!(one.next_after(&one, p, rm).unwrap().cmp(&one) == 0);
    }
}
