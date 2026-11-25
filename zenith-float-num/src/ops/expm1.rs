//! `e^x − 1` without cancellation for small `x`.

use crate::common::util::{bump_prec_retry, round_p};
use crate::defs::{Error, Sign, WORD_BIT_SIZE};
use crate::num::ExactNumNumber;
use crate::ops::consts::Consts;
use crate::ops::util::compute_small_exp;
use crate::Exponent;
use crate::RoundingMode;

impl ExactNumNumber {
    /// Computes `exp(self) - 1` with precision `p`. The result is rounded using `rm`.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - ExponentOverflow: `exp(self)` overflows the exponent range.
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - InvalidArgument: the precision is incorrect.
    pub fn expm1(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;

        if self.is_zero() {
            return Self::new2(p, self.sign(), self.inexact());
        }

        let mut p_inc = WORD_BIT_SIZE;
        // The retry budget counts from the working start: an argument wider than `p` must not exhaust it at once.
        let p_base = p.max(self.mantissa_max_bit_len());
        let mut p_wrk = p_base;

        // `exp(x) - 1 = x + x^2/2 + ...`: the true value is above `x` in magnitude for `x > 0` and below it for
        // `x < 0`, so the one-bit correction goes against the sign of `x` only when `x` is negative.
        compute_small_exp!(self, self.exponent(), self.is_negative(), p_wrk, p, rm);

        p_wrk += p_inc;

        loop {
            let p_x = p_wrk + 3;
            let e = match self.exp(p_x, RoundingMode::None, cc) {
                Err(Error::ExponentOverflow(s)) => {
                    if s == Sign::Neg {
                        // e^x underflowed: the true value is just above -1, round once.
                        let mut m1 = Self::from_word(1, p)?;
                        m1.set_sign(Sign::Neg);
                        let mut r = Self::round_offset(&m1, Sign::Pos, p, rm)?;
                        r.set_inexact(true);
                        return Ok(r);
                    }
                    return Err(Error::ExponentOverflow(s));
                }
                Err(e) => return Err(e),
                Ok(v) => v,
            };

            let one = Self::from_word(1, p_x)?;

            // ZF-37: e^x < 2^-(p+1) is below half an ulp of the result's binade (-1, -1/2), so the
            // value -1 + e^x rounds like "just above -1"; (e^x - 1)·x/ln(e^x) cannot decide it at
            // any width (e.g. x = -4.7e10: e^x = 2^-6.8e10).
            if self.is_negative() && !e.is_zero() && e.exponent() < -(p as Exponent) {
                let mut m1 = Self::from_word(1, p)?;
                m1.set_sign(Sign::Neg);
                let mut r = Self::round_offset(&m1, Sign::Pos, p, rm)?;
                r.set_inexact(true);
                return Ok(r);
            }

            if e.is_zero() {
                // e^x is below every p_x-bit number: the true value is just above -1.
                let mut m1 = Self::from_word(1, p)?;
                m1.set_sign(Sign::Neg);
                let mut r = Self::round_offset(&m1, Sign::Pos, p, rm)?;
                r.set_inexact(true);
                return Ok(r);
            }

            let mut ret = if e.cmp(&one) == 0 {
                self.try_rounded_to(p_x, RoundingMode::None)?
            } else {
                let em1 = e.sub(&one, p_x, RoundingMode::None)?;
                let ln_e = e.ln(p_x, RoundingMode::None, cc)?;
                if ln_e.cmp(self) == 0 {
                    em1
                } else {
                    em1.mul(self, p_x, RoundingMode::None)?
                        .div(&ln_e, p_x, RoundingMode::None)?
                }
            };

            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }

            bump_prec_retry(&mut p_wrk, &mut p_inc, p_base)?;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::defs::WORD_BIT_SIZE;
    use crate::RoundingMode;

    #[test]
    fn test_expm1() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;
        let p = WORD_BIT_SIZE * 5;

        let z = ExactNumNumber::new(p)?;
        assert!(z.expm1(p, rm, &mut cc)?.is_zero());

        let one = ExactNumNumber::small_const(1);
        let e1 = one.expm1(p, rm, &mut cc)?;
        assert!(e1.is_positive());
        assert!(!e1.is_zero());

        for _ in 0..40 {
            let x = ExactNumNumber::random_normal(p, -12, -2)?;
            let y = x.expm1(p, rm, &mut cc)?;
            let back = y.log1p(p, rm, &mut cc)?;
            let mut eps = ExactNumNumber::small_const(1);
            eps.set_exponent(x.exponent() - p as crate::Exponent + 6);
            assert!(back.sub(&x, p, rm)?.abs()?.cmp(&eps) < 0);
        }
        Ok(())
    }
}
