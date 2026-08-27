//! Two-argument arctangent.

use crate::common::util::{bump_prec_retry, round_p};
use crate::defs::{Error, WORD_BIT_SIZE};
use crate::num::ExactNumNumber;
use crate::ops::consts::Consts;
use crate::{RoundingMode, Sign};

impl ExactNumNumber {
    /// Computes `atan2(self, x)` (quadrant-aware arctangent of `self/x`) with precision `p`.
    /// The result is rounded using the rounding mode `rm`.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - ExponentOverflow: the resulting exponent is out of range.
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - InvalidArgument: the precision is incorrect.
    pub fn atan2(
        &self,
        x: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;

        let y = self;
        let inexact = y.inexact() | x.inexact();

        if x.is_zero() {
            if y.is_zero() {
                return if x.is_negative() {
                    Self::signed_pi(y.sign(), p, rm, cc, inexact)
                } else {
                    Self::new2(p, y.sign(), inexact)
                };
            }
            return Self::signed_half_pi(y.sign(), p, rm, cc, inexact);
        }

        if y.is_zero() {
            return if x.is_positive() {
                Self::new2(p, y.sign(), inexact)
            } else {
                Self::signed_pi(y.sign(), p, rm, cc, inexact)
            };
        }

        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p + p_inc;

        loop {
            let p_x = p_wrk + 3;

            let mut ret = match y.div(x, p_x, RoundingMode::None) {
                Err(Error::ExponentOverflow(_)) => {
                    if x.is_positive() {
                        Self::signed_half_pi(y.sign(), p_x, RoundingMode::None, cc, inexact)?
                    } else {
                        Self::signed_pi(y.sign(), p_x, RoundingMode::None, cc, inexact)?
                    }
                }
                Err(e) => return Err(e),
                Ok(q) => {
                    let a = q.atan(p_x, RoundingMode::None, cc)?;
                    if x.is_positive() {
                        a
                    } else {
                        let pi = Self::signed_pi(y.sign(), p_x, RoundingMode::None, cc, false)?;
                        if y.is_positive() {
                            a.add(&pi, p_x, RoundingMode::None)?
                        } else {
                            a.sub(&pi, p_x, RoundingMode::None)?
                        }
                    }
                }
            };

            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | inexact);
                return Ok(ret);
            }

            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn signed_half_pi(
        s: Sign,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
        inexact: bool,
    ) -> Result<Self, Error> {
        let mut h = cc.pi_num(p, rm)?;
        h.set_exponent(1);
        h.set_sign(s);
        h.set_inexact(h.inexact() | inexact);
        Ok(h)
    }

    fn signed_pi(
        s: Sign,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
        inexact: bool,
    ) -> Result<Self, Error> {
        let mut pi = cc.pi_num(p, rm)?;
        pi.set_sign(s);
        pi.set_inexact(pi.inexact() | inexact);
        Ok(pi)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::consts::ONE;
    use crate::defs::WORD_BIT_SIZE;
    use crate::RoundingMode;

    #[test]
    fn test_atan2() {
        let mut cc = Consts::new().unwrap();
        let rm = RoundingMode::ToEven;
        let p = WORD_BIT_SIZE * 5;

        let z = ExactNumNumber::new(p).unwrap();
        let one = ONE.clone().unwrap();
        assert!(z.atan2(&one, p, rm, &mut cc).unwrap().is_zero());

        let y = ExactNumNumber::from_i8(1, p).unwrap();
        let x = ExactNumNumber::from_i8(1, p).unwrap();
        let a = y.atan2(&x, p, rm, &mut cc).unwrap();
        let at = y.atan(p, rm, &mut cc).unwrap();
        // atan2(1,1) = π/4 = atan(1)
        assert_eq!(a.cmp(&at), 0);

        let nx = x.neg().unwrap();
        let a2 = y.atan2(&nx, p, rm, &mut cc).unwrap();
        let mut pi = cc.pi_num(p, rm).unwrap();
        pi.set_exponent(1); // π/2
        let three_q = pi.add(&at, p, rm).unwrap();
        let mut eps = ONE.clone().unwrap();
        eps.set_exponent(a2.exponent() - p as crate::Exponent + 4);
        assert!(a2.sub(&three_q, p, rm).unwrap().abs().unwrap().cmp(&eps) < 0);

        for _ in 0..30 {
            let yy = ExactNumNumber::random_normal(p, -8, 8).unwrap();
            if yy.is_zero() {
                continue;
            }
            let a = yy.atan2(&one, p, rm, &mut cc).unwrap();
            let b = yy.atan(p, rm, &mut cc).unwrap();
            eps.set_exponent(a.exponent() - p as crate::Exponent + 4);
            assert!(a.sub(&b, p, rm).unwrap().abs().unwrap().cmp(&eps) < 0);
        }
    }
}
