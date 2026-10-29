//! Two-argument arctangent.

use crate::common::util::round_p_rm;
use crate::common::util::{bump_prec_retry, invert_rm_for_sign};
use crate::defs::{Error, EXPONENT_MIN, WORD_BIT_SIZE};
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
        let p = round_p_rm(p, rm);
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

        // `|y / x| < 2^(ey - ex + 1)`. When that is below half the smallest number of precision `p`, the angle
        // (`y / x` less a part in 2^(2^31)) is rounded like the quotient itself, which is in turn rounded to zero or
        // to the smallest number by the mode. The retry loop cannot resolve such a value: nothing it computes
        // has a nonzero part to show at any working precision.
        // In `i128`: with the i64 exponent range `ey - ex + 2 + p` reaches past `i64::MAX` for a huge `y` over a
        // tiny `x`, and a wrapped sum took this branch (an overflowing quotient, `Inf` for `+-pi/2`).
        if x.is_positive()
            && i128::from(y.exponent()) - i128::from(x.exponent()) + 2 + p as i128
                <= i128::from(EXPONENT_MIN)
        {
            let mut q = y.div(x, p, rm)?;
            q.set_inexact(true);
            return Ok(q);
        }

        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p + p_inc;

        loop {
            let p_x = p_wrk + 3;

            let mut ret = match y.div(x, p_x, RoundingMode::None) {
                // |y / x| is beyond the largest number: the angle is `+-pi/2` for either sign of `x`.
                Err(Error::ExponentOverflow(_)) => {
                    Self::signed_half_pi(y.sign(), p_x, RoundingMode::None, cc, inexact)?
                }
                Err(e) => return Err(e),
                Ok(q) => {
                    if x.is_positive() && !q.op_inexact() {
                        // An exact quotient: `atan(q)` is `q` less a correction that no working precision can
                        // show, so rounding `q` toward `atan(q)` needs `atan`'s own treatment of a small
                        // argument. (The retry loop below never decides a directed rounding of a number that is
                        // exactly representable.)
                        let mut a = q.atan(p, rm, cc)?;
                        a.set_inexact(a.inexact() | inexact);
                        return Ok(a);
                    }
                    let a = q.atan(p_x, RoundingMode::None, cc)?;
                    if x.is_positive() {
                        a
                    } else {
                        // `signed_pi` already carries `y`'s sign: Q2 is `atan + π`, Q3 is `atan − π`.
                        let pi = Self::signed_pi(y.sign(), p_x, RoundingMode::None, cc, false)?;
                        a.add(&pi, p_x, RoundingMode::None)?
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
        let rm = if s.is_negative() { invert_rm_for_sign(rm) } else { rm };
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
        let rm = if s.is_negative() { invert_rm_for_sign(rm) } else { rm };
        let mut pi = cc.pi_num(p, rm)?;
        pi.set_sign(s);
        pi.set_inexact(pi.inexact() | inexact);
        Ok(pi)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::defs::WORD_BIT_SIZE;
    use crate::RoundingMode;

    #[test]
    fn test_atan2() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;
        let p = WORD_BIT_SIZE * 5;

        let z = ExactNumNumber::new(p)?;
        let one = ExactNumNumber::small_const(1);
        assert!(z.atan2(&one, p, rm, &mut cc)?.is_zero());

        let y = ExactNumNumber::from_i8(1, p)?;
        let x = ExactNumNumber::from_i8(1, p)?;
        let a = y.atan2(&x, p, rm, &mut cc)?;
        let at = y.atan(p, rm, &mut cc)?;
        // atan2(1,1) = π/4 = atan(1)
        assert_eq!(a.cmp(&at), 0);

        let nx = x.neg()?;
        let a2 = y.atan2(&nx, p, rm, &mut cc)?;
        let mut pi = cc.pi_num(p, rm)?;
        pi.set_exponent(1); // π/2
        let three_q = pi.add(&at, p, rm)?;
        let mut eps = ExactNumNumber::small_const(1);
        eps.set_exponent(a2.exponent() - p as crate::Exponent + 4);
        assert!(a2.sub(&three_q, p, rm)?.abs()?.cmp(&eps) < 0);

        let ny = y.neg()?;
        let a3 = ny.atan2(&nx, p, rm, &mut cc)?;
        let three_q_neg = three_q.neg()?;
        eps.set_exponent(a3.exponent() - p as crate::Exponent + 4);
        assert!(a3.sub(&three_q_neg, p, rm)?.abs()?.cmp(&eps) < 0);
        assert!(a3.is_negative());

        for _ in 0..30 {
            let yy = ExactNumNumber::random_normal(p, -8, 8)?;
            if yy.is_zero() {
                continue;
            }
            let a = yy.atan2(&one, p, rm, &mut cc)?;
            let b = yy.atan(p, rm, &mut cc)?;
            eps.set_exponent(a.exponent() - p as crate::Exponent + 4);
            assert!(a.sub(&b, p, rm)?.abs()?.cmp(&eps) < 0);
        }
        Ok(())
    }
}
