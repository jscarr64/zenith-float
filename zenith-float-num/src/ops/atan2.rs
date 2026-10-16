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

        // A quotient this small is `atan(q) = q - q³/3 + ...`, a hair below `q`: no working precision below
        // about `-2 * exponent(q)` bits can tell `q` from the result, so the retry loop below never ends in
        // the directed modes (`PrecisionRetryExhausted`). The result is `q` rounded, or `q` minus a hair.
        if x.is_positive() {
            let p_q = p + 2 * WORD_BIT_SIZE;
            let e_est = (y.exponent() as i64).saturating_sub(x.exponent() as i64);
            let bits = y
                .mantissa_max_bit_len()
                .saturating_add(x.mantissa_max_bit_len());
            if e_est.saturating_mul(2) < -((p_q + bits + 2 * WORD_BIT_SIZE) as i64) {
                return Self::atan2_tiny_quotient(y, x, p, rm, p_q);
            }
        }

        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p + p_inc;

        loop {
            let p_x = p_wrk + 3;

            let mut ret = match y.div(x, p_x, RoundingMode::None) {
                // `|y/x|` is beyond the exponent range: the angle is `±π/2` for either sign of `x`
                // (`x` is negligible next to `y`), with the sign of `y`.
                Err(Error::ExponentOverflow(_)) => {
                    Self::signed_half_pi(y.sign(), p_x, RoundingMode::None, cc, inexact)?
                }
                Err(e) => return Err(e),
                Ok(q) => {
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

    /// `atan2(y, x)` for `x > 0` and `|y / x|` below `2^-(p_q + bits of y and x + 2 words) / 2`.
    fn atan2_tiny_quotient(
        y: &Self,
        x: &Self,
        p: usize,
        rm: RoundingMode,
        p_q: usize,
    ) -> Result<Self, Error> {
        let q = y.div(x, p_q, RoundingMode::None)?;
        if q.is_zero() || q.inexact() {
            // `q` underflowed, or `y / x` is not a `p_q`-bit number: its distance to a rounding boundary is
            // far above `q³/3`, so the result rounds like `y / x`.
            let mut ret = y.div(x, p, rm)?;
            ret.set_inexact(true);
            return Ok(ret);
        }
        // `y / x` is exactly `q`: the result is below it by less than the last bit of `q`. One bit down and
        // "something more" (the inexact flag is the sticky bit of the final rounding) is that open interval.
        let k = (q.exponent() as i64).saturating_sub(q.mantissa_max_bit_len() as i64);
        let unit = Self::pow2_unit(q.sign(), k, p_q)?;
        let mut ret = q.sub(&unit, p_q, RoundingMode::None)?;
        ret.set_inexact(true);
        ret.set_precision(p, rm)?;
        // the result differs from `q`, so it is inexact whatever the operands' flags
        ret.set_inexact(true);
        Ok(ret)
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

        let ny = y.neg().unwrap();
        let a3 = ny.atan2(&nx, p, rm, &mut cc).unwrap();
        let three_q_neg = three_q.neg().unwrap();
        eps.set_exponent(a3.exponent() - p as crate::Exponent + 4);
        assert!(
            a3.sub(&three_q_neg, p, rm)
                .unwrap()
                .abs()
                .unwrap()
                .cmp(&eps)
                < 0
        );
        assert!(a3.is_negative());

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
