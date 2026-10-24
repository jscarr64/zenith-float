//! `ln(1 + x)` without cancellation for small `x`.

use crate::common::util::{bump_prec_retry, round_p};
use crate::defs::{Error, Exponent, WORD_BIT_SIZE};
use crate::num::ExactNumNumber;
use crate::ops::consts::Consts;
use crate::ops::util::compute_small_exp;
use crate::RoundingMode;

impl ExactNumNumber {
    /// Computes `ln(1 + self)` with precision `p`. The result is rounded using `rm`.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - InvalidArgument: `self <= -1`, or the precision is incorrect.
    ///  - DivisionByZero: `self == -1` (the result is `-Inf` at the `ExactNum` layer).
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    pub fn log1p(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;

        if self.is_zero() {
            return Self::new2(p, self.sign(), self.inexact());
        }

        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len());

        // `ln(1 + x) = x - x^2/2 + ...`: the true value is below `x` in magnitude for `x > 0` and above it for
        // `x < 0`, so the one-bit correction goes against the sign of `x` only when `x` is positive.
        compute_small_exp!(self, self.exponent(), self.is_positive(), p_wrk, p, rm);

        p_wrk += p_inc;

        loop {
            let p_x = p_wrk + 3;
            let one = Self::from_word(1, p_x)?;

            // Huge x: `1 + x` can overflow, so use `ln(1 + x) = ln(x) + 1/x - ...`, whose tail
            // `1/(2x²)` is far below the working precision.
            let p_exp = Exponent::try_from(p_x).map_err(|_| Error::InvalidArgument)?;
            if self.is_positive() && self.exponent() > p_exp + 2 {
                let r = one.div(self, p_x, RoundingMode::None)?;
                let mut ret =
                    self.ln(p_x, RoundingMode::None, cc)?
                        .add(&r, p_x, RoundingMode::None)?;
                if ret.try_set_precision(p, rm, p_wrk)? {
                    ret.set_inexact(true);
                    return Ok(ret);
                }
                bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
                continue;
            }

            let u = one.add(self, p_x, RoundingMode::None)?;

            if u.is_zero() {
                return Err(Error::DivisionByZero);
            }
            if u.is_negative() {
                return Err(Error::InvalidArgument);
            }

            let mut ret = if u.cmp(&one) == 0 {
                let mut x = self.clone()?;
                x.set_precision(p_x, RoundingMode::None)?;
                x
            } else {
                let ln_u = u.ln(p_x, RoundingMode::None, cc)?;
                let um1 = u.sub(&one, p_x, RoundingMode::None)?;
                if um1.cmp(self) == 0 {
                    ln_u
                } else {
                    ln_u.mul(self, p_x, RoundingMode::None)?
                        .div(&um1, p_x, RoundingMode::None)?
                }
            };

            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }

            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::defs::WORD_BIT_SIZE;
    use crate::RoundingMode;

    #[test]
    fn test_log1p() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;
        let p = WORD_BIT_SIZE * 5;

        let z = ExactNumNumber::new(p)?;
        assert!(z.log1p(p, rm, &mut cc)?.is_zero());

        let one = ExactNumNumber::small_const(1);
        let ln2 = one.log1p(p, rm, &mut cc)?;
        let ln2b = ExactNumNumber::from_i8(2, p)?.ln(p, rm, &mut cc)?;
        assert_eq!(ln2.cmp(&ln2b), 0);

        assert!(matches!(
            one.neg()?.log1p(p, rm, &mut cc),
            Err(Error::DivisionByZero)
        ));
        assert!(matches!(
            ExactNumNumber::from_i8(-2, p)?.log1p(p, rm, &mut cc),
            Err(Error::InvalidArgument)
        ));

        for _ in 0..40 {
            let x = ExactNumNumber::random_normal(p, -30, -2)?;
            if x.is_negative() {
                continue;
            }
            let l = x.log1p(p, rm, &mut cc)?;
            let u = one.add(&x, p, rm)?;
            let back = l.exp(p, rm, &mut cc)?;
            let mut eps = ExactNumNumber::small_const(1);
            eps.set_exponent(u.exponent() - p as crate::Exponent + 5);
            assert!(back.sub(&u, p, rm)?.abs()?.cmp(&eps) < 0);
        }
        Ok(())
    }
}
