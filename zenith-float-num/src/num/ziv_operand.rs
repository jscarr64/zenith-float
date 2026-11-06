//! The one approved Ziv working-copy helper (Jeff's exception, tools/NO_CLONE_RULES.md "Exception: Ziv working copy").
//!
//! A correct-rounding retry loop reads the same exact input at a larger precision on every pass, and the kernels
//! consume their argument, so each pass needs a fresh owned value of the input at that pass's precision. This file holds
//! the only function that makes one: [`ExactNumNumber::ziv_pass_operand`], one allocation plus one limb copy per pass.
//! The scan (`tools/rules_scan.sh`) accepts exactly one such helper, in this file, protected by CODEOWNERS.

use super::ExactNumNumber;
use crate::defs::{Error, WORD_BIT_SIZE};
use crate::mantissa::Mantissa;

impl ExactNumNumber {
    /// A fresh owned copy of `self` with room for `p` bits, for one retry pass: one allocation (the new mantissa, zero
    /// filled) and one limb copy (the words of `self` into the top of it). The value, the sign, the exponent and both
    /// inexact flags are carried over unchanged, and nothing is rounded, so the precision may only grow: `p` below the
    /// precision of `self` is refused. A precision between two word sizes rounds up to the word size, as everywhere.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for the mantissa.
    ///  - InvalidArgument: `p` is not a valid precision, or is smaller than the precision of `self`.
    pub(crate) fn ziv_pass_operand(&self, p: usize) -> Result<Self, Error> {
        Self::p_assertion(p)?;
        let src = self.m.digits();
        let mut m = Mantissa::new(p)?;
        let lift = m
            .digits()
            .len()
            .checked_sub(src.len())
            .ok_or(Error::InvalidArgument)?;
        m.digits_mut()[lift..].copy_from_slice(src);
        if self.m.bit_len() != 0 {
            m.set_bit_len(self.m.bit_len() + lift * WORD_BIT_SIZE);
        }
        Ok(ExactNumNumber {
            e: self.e,
            s: self.s,
            m,
            inexact: self.inexact,
            op_inexact: self.op_inexact,
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::defs::{Error, Sign, Word, EXPONENT_MIN, WORD_BIT_SIZE, WORD_SIGNIFICANT_BIT};
    use crate::mantissa::CLONE_CALLS;
    use crate::num::ExactNumNumber;

    const TOP: Word = WORD_SIGNIFICANT_BIT;
    /// A word value that fits the 32-bit word of the 32-bit targets too.
    const LOW: Word = 0x0123_4567;
    const W: usize = WORD_BIT_SIZE;

    fn clones() -> usize {
        CLONE_CALLS.with(|c| c.get())
    }

    /// A two-word number widened to four words: the two words move to the top of a four-word mantissa, the exponent, the
    /// sign and the bit length follow, nothing is rounded. The expected words are written out by hand.
    #[test]
    fn widening_moves_the_words_to_the_top() -> Result<(), Error> {
        let x = ExactNumNumber::from_words(&[LOW, TOP | 3], Sign::Neg, 5)?;
        let y = x.ziv_pass_operand(4 * W)?;
        assert_eq!(y.mantissa().digits(), &[0, 0, LOW, TOP | 3]);
        assert_eq!(y.mantissa().bit_len(), 4 * W);
        assert_eq!(y.exponent(), 5);
        assert!(y.is_negative());
        assert!(!y.inexact());
        Ok(())
    }

    /// The same precision gives the same words; a precision between two word sizes rounds up to the word size.
    #[test]
    fn same_precision_and_in_between_precisions() -> Result<(), Error> {
        let x = ExactNumNumber::from_words(&[7, TOP | 9], Sign::Pos, -3)?;
        let same = x.ziv_pass_operand(2 * W)?;
        assert_eq!(same.mantissa().digits(), &[7, TOP | 9]);
        assert_eq!(same.exponent(), -3);
        let between = x.ziv_pass_operand(2 * W + 1)?;
        assert_eq!(between.mantissa().digits(), &[0, 7, TOP | 9]);
        assert_eq!(between.mantissa().bit_len(), 3 * W);
        Ok(())
    }

    /// A subnormal operand keeps its short bit length (it only gains the new low words), at the smallest exponent.
    #[test]
    fn a_subnormal_operand_keeps_its_exponent_and_gains_low_words() -> Result<(), Error> {
        let x = ExactNumNumber::from_words(&[1, 1], Sign::Pos, EXPONENT_MIN)?;
        assert_eq!(x.mantissa().bit_len(), W + 1);
        let y = x.ziv_pass_operand(4 * W)?;
        assert_eq!(y.mantissa().digits(), &[0, 0, 1, 1]);
        assert_eq!(y.mantissa().bit_len(), W + 1 + 2 * W);
        assert_eq!(y.exponent(), EXPONENT_MIN);
        Ok(())
    }

    /// Zero stays zero (bit length 0, exponent 0) at any precision.
    #[test]
    fn zero_stays_zero() -> Result<(), Error> {
        let x = ExactNumNumber::from_words(&[0, 0], Sign::Pos, 17)?;
        let y = x.ziv_pass_operand(3 * W)?;
        assert_eq!(y.mantissa().digits(), &[0, 0, 0]);
        assert_eq!(y.mantissa().bit_len(), 0);
        assert!(y.is_zero());
        Ok(())
    }

    /// Both inexact flags are carried over, as a copy carries them.
    #[test]
    fn the_inexact_flags_are_carried_over() -> Result<(), Error> {
        let mut x = ExactNumNumber::from_words(&[1, TOP], Sign::Pos, 1)?;
        x.set_inexact(true);
        assert!(x.ziv_pass_operand(3 * W)?.inexact());
        x.set_op_inexact(false);
        let y = x.ziv_pass_operand(3 * W)?;
        assert!(y.inexact());
        assert!(!y.op_inexact());
        Ok(())
    }

    /// Narrowing is not what a retry pass does (precision only grows); it is refused, not rounded.
    #[test]
    fn narrowing_is_refused() -> Result<(), Error> {
        let x = ExactNumNumber::from_words(&[1, TOP], Sign::Pos, 1)?;
        assert!(matches!(x.ziv_pass_operand(W), Err(Error::InvalidArgument)));
        assert!(matches!(x.ziv_pass_operand(0), Err(Error::InvalidArgument)));
        Ok(())
    }

    /// The helper is its own copy primitive: it does not go through `Mantissa::clone`.
    #[test]
    fn it_does_not_call_mantissa_clone() -> Result<(), Error> {
        let x = ExactNumNumber::from_words(&[1, TOP], Sign::Pos, 1)?;
        let before = clones();
        let _y = x.ziv_pass_operand(1024)?;
        assert_eq!(clones(), before);
        Ok(())
    }
}
