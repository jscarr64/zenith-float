// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! Hyperbolic sine.

use crate::common::util::bump_prec_retry;
use crate::common::util::round_p;
use crate::defs::Error;
use crate::defs::RoundingMode;
use crate::num::ExactNumNumber;
use crate::ops::util::compute_small_exp;
use crate::Consts;
use crate::Sign;
use crate::WORD_BIT_SIZE;

impl ExactNumNumber {
    /// Computes the hyperbolic sine of a number with precision `p`. The result is rounded using the rounding mode `rm`.
    /// This function requires constants cache cc for computing the result.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory.
    ///  - InvalidArgument: the precision is incorrect.
    pub fn sinh(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        let p = round_p(p);

        if self.is_zero() {
            return Self::new2(p, self.sign(), self.inexact());
        }

        let mut p_inc = WORD_BIT_SIZE;
        // The retry budget counts from the working start: an argument wider than `p` must not exhaust it at once.
        let p_base = p.max(self.mantissa_max_bit_len());
        let mut p_wrk = p_base;

        compute_small_exp!(self, self.exponent() * 2 - 2, false, p_wrk, p, rm);

        p_wrk += p_inc;

        let mut x = self.clone()?;
        x.set_inexact(false);
        x.set_sign(Sign::Pos);

        // In u64: on 32-bit the doubled exponent can exceed `usize`, and then it exceeds any bit length.
        let ethres = (x.exponent().unsigned_abs().max(1) - 1).checked_mul(2);

        loop {
            let p_x = p_wrk + 4;
            x.set_precision(p_x, RoundingMode::None)?;

            let mut ret = if x.exponent() <= 0 {
                Self::sinh_series(x.clone()?, p_x, RoundingMode::None)?
            } else {
                let mut val = if ethres
                    .and_then(|t| usize::try_from(t).ok())
                    .is_none_or(|t| t > x.mantissa_max_bit_len() + 2)
                {
                    // e^|x| / 2 * signum(self)

                    x.exp(p_x, RoundingMode::None, cc)
                } else {
                    // (e^x - e^(-x)) / 2

                    let ex = match x.exp(p_x, RoundingMode::None, cc) {
                        Ok(val) => val,
                        Err(Error::ExponentOverflow(_)) => {
                            return Err(Error::ExponentOverflow(self.sign()));
                        }
                        Err(err) => return Err(err),
                    };

                    let xe = ex.reciprocal(p_x, RoundingMode::None)?;

                    ex.sub(&xe, p_x, RoundingMode::None)
                }
                .map_err(|e| -> Error {
                    if let Error::ExponentOverflow(_) = e {
                        Error::ExponentOverflow(self.sign())
                    } else {
                        e
                    }
                })?;

                val.div_by_2(RoundingMode::None);

                val
            };

            ret.set_sign(self.sign());

            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                break Ok(ret);
            }

            bump_prec_retry(&mut p_wrk, &mut p_inc, p_base)?;
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::{common::util::random_subnormal, Sign};

    use super::*;

    #[test]
    fn test_sinh() -> Result<(), crate::Error> {
        let p = 32000;
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;
        let mut n1 = ExactNumNumber::from_word(1, 1)?;
        n1.set_exponent(0);
        let _n2 = n1.sinh(p, rm, &mut cc)?;
        //println!("{:?}", n2.fp3(crate::Radix::Dec, rm).unwrap());

        let d1 = ExactNumNumber::max_value(p)?;
        let d2 = ExactNumNumber::min_value(p)?;

        assert!(matches!(
            d1.sinh(p, rm, &mut cc),
            Err(Error::ExponentOverflow(Sign::Pos))
        ));
        assert!(matches!(
            d2.sinh(p, rm, &mut cc),
            Err(Error::ExponentOverflow(Sign::Neg))
        ));

        let d3 = ExactNumNumber::min_positive(p)?;
        let zero = ExactNumNumber::new(1)?;
        let n1 = random_subnormal(p)?;

        assert!(d3.sinh(p, rm, &mut cc)?.cmp(&d3) == 0);
        assert!(zero.sinh(p, rm, &mut cc)?.is_zero());
        assert!(n1.sinh(p, rm, &mut cc)?.cmp(&n1) == 0);

        let mut large_pos = ExactNumNumber::from_word(1, 256)?;
        large_pos.set_exponent(150);
        assert!(matches!(
            large_pos.sinh(p, rm, &mut cc),
            Err(Error::ExponentOverflow(Sign::Pos))
        ));
        let mut large_neg = ExactNumNumber::from_word(1, 256)?;
        large_neg.set_exponent(150);
        large_neg.set_sign(Sign::Neg);
        assert!(matches!(
            large_neg.sinh(p, rm, &mut cc),
            Err(Error::ExponentOverflow(Sign::Neg))
        ));
        Ok(())
    }
}
