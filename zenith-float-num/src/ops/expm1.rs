//! `e^x − 1` without cancellation for small `x`.

use crate::common::util::{bump_prec_retry, round_p};
use crate::defs::{Error, Sign, WORD_BIT_SIZE};
use crate::num::ExactNumNumber;
use crate::ops::consts::Consts;
use crate::ops::util::compute_small_exp;
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
        let mut p_wrk = p.max(self.mantissa_max_bit_len());

        compute_small_exp!(self, self.exponent() as isize, false, p_wrk, p, rm);

        p_wrk += p_inc;

        loop {
            let p_x = p_wrk + 3;
            let e = match self.exp(p_x, RoundingMode::None, cc) {
                Err(Error::ExponentOverflow(s)) => {
                    if s == Sign::Neg {
                        let mut m1 = Self::from_word(1, p)?;
                        m1.set_sign(Sign::Neg);
                        m1.set_inexact(true | self.inexact());
                        return Ok(m1);
                    }
                    return Err(Error::ExponentOverflow(s));
                }
                Err(e) => return Err(e),
                Ok(v) => v,
            };

            let one = Self::from_word(1, p_x)?;

            if e.is_zero() {
                let mut m1 = Self::from_word(1, p)?;
                m1.set_sign(Sign::Neg);
                m1.set_inexact(m1.inexact() | self.inexact());
                return Ok(m1);
            }

            let mut ret = if e.cmp(&one) == 0 {
                let mut x = self.clone()?;
                x.set_precision(p_x, RoundingMode::None)?;
                x
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

            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::consts::ONE;
    use crate::defs::WORD_BIT_SIZE;
    use crate::RoundingMode;

    #[test]
    fn test_expm1() {
        let mut cc = Consts::new().unwrap();
        let rm = RoundingMode::ToEven;
        let p = WORD_BIT_SIZE * 5;

        let z = ExactNumNumber::new(p).unwrap();
        assert!(z.expm1(p, rm, &mut cc).unwrap().is_zero());

        let one = ONE.clone().unwrap();
        let e1 = one.expm1(p, rm, &mut cc).unwrap();
        assert!(e1.is_positive());
        assert!(!e1.is_zero());

        for _ in 0..40 {
            let x = ExactNumNumber::random_normal(p, -12, -2).unwrap();
            let y = x.expm1(p, rm, &mut cc).unwrap();
            let back = y.log1p(p, rm, &mut cc).unwrap();
            let mut eps = ONE.clone().unwrap();
            eps.set_exponent(x.exponent() - p as crate::Exponent + 6);
            assert!(back.sub(&x, p, rm).unwrap().abs().unwrap().cmp(&eps) < 0);
        }
    }
}
