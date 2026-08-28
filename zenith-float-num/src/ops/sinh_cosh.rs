//! Paired hyperbolic sine and cosine.

use crate::common::util::bump_prec_retry;
use crate::common::util::round_p;
use crate::defs::Error;
use crate::defs::RoundingMode;
use crate::num::ExactNumNumber;
use crate::Consts;
use crate::Sign;
use crate::WORD_BIT_SIZE;

impl ExactNumNumber {
    /// Computes `(sinh(self), cosh(self))` with precision `p`.
    ///
    /// Both results share a single `exp(|x|)` evaluation; each output is rounded independently at `p`.
    ///
    /// ## Errors
    ///
    ///  - ExponentOverflow: a result is too large or too small.
    ///  - MemoryAllocation: failed to allocate memory.
    ///  - InvalidArgument: the precision is incorrect.
    pub fn sinh_cosh(
        &self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<(Self, Self), Error> {
        let p = round_p(p);

        if self.is_zero() {
            let sinh = Self::new2(p, self.sign(), self.inexact())?;
            let mut cosh = Self::from_word(1, p)?;
            cosh.set_inexact(self.inexact());
            return Ok((sinh, cosh));
        }

        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len());

        if p_wrk as isize + 1 < -(self.exponent() as isize * 2 - 2) {
            let mut x = self.clone()?;
            if p > x.mantissa_max_bit_len() {
                x.set_precision(p, RoundingMode::None)?;
            }
            let mut sinh = x.add_correction(false)?;
            sinh.set_precision(p, rm)?;
            sinh.set_sign(self.sign());
            let mut cosh = Self::from_word(1, p)?;
            cosh.set_inexact(self.inexact());
            return Ok((sinh, cosh));
        }

        p_wrk += p_inc;

        let mut x = self.clone()?;
        x.set_inexact(false);
        x.set_sign(Sign::Pos);

        let ethres = (x.exponent().unsigned_abs().max(1) as usize - 1) * 2;

        loop {
            let p_x = p_wrk + 4;
            x.set_precision(p_x, RoundingMode::None)?;

            let (mut sinh, mut cosh) = if x.exponent() <= 0 {
                let sinh = Self::sinh_series(x.clone()?, p_x, RoundingMode::None)?;
                let sh2 = sinh.mul(&sinh, p_x, RoundingMode::None)?;
                let one = Self::from_word(1, p_x)?;
                let cosh = one
                    .add(&sh2, p_x, RoundingMode::None)?
                    .sqrt(p_x, RoundingMode::None)?;
                (sinh, cosh)
            } else if ethres > x.mantissa_max_bit_len() + 2 {
                let mut ex = match x.exp(p_x, RoundingMode::None, cc) {
                    Ok(val) => val,
                    Err(Error::ExponentOverflow(_)) => {
                        return Err(Error::ExponentOverflow(self.sign()));
                    }
                    Err(err) => return Err(err),
                };
                let mut sinh = ex.clone()?;
                sinh.div_by_2(RoundingMode::None);
                ex.div_by_2(RoundingMode::None);
                (sinh, ex)
            } else {
                let ex = match x.exp(p_x, RoundingMode::None, cc) {
                    Ok(val) => val,
                    Err(Error::ExponentOverflow(_)) => {
                        return Err(Error::ExponentOverflow(self.sign()));
                    }
                    Err(err) => return Err(err),
                };

                let xe = ex.reciprocal(p_x, RoundingMode::None)?;

                let mut sinh = ex.sub(&xe, p_x, RoundingMode::None).map_err(|e| -> Error {
                    if let Error::ExponentOverflow(_) = e {
                        Error::ExponentOverflow(self.sign())
                    } else {
                        e
                    }
                })?;
                sinh.div_by_2(RoundingMode::None);

                let mut cosh = ex.add(&xe, p_x, RoundingMode::None)?;
                cosh.div_by_2(RoundingMode::None);

                (sinh, cosh)
            };

            sinh.set_sign(self.sign());

            let sinh_ok = sinh.try_set_precision(p, rm, p_wrk)?;
            let cosh_ok = cosh.try_set_precision(p, rm, p_wrk)?;

            if sinh_ok && cosh_ok {
                sinh.set_inexact(sinh.inexact() | self.inexact());
                cosh.set_inexact(cosh.inexact() | self.inexact());
                break Ok((sinh, cosh));
            }

            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::{common::util::random_subnormal, Sign};

    #[test]
    fn test_sinh_cosh() {
        let p = 320;
        let mut cc = Consts::new().unwrap();
        let rm = RoundingMode::ToEven;

        let d1 = ExactNumNumber::from_word(1, p).unwrap();
        let (s, c) = d1.sinh_cosh(p, rm, &mut cc).unwrap();
        assert!(s.cmp(&d1.sinh(p, rm, &mut cc).unwrap()) == 0);
        assert!(c.cmp(&d1.cosh(p, rm, &mut cc).unwrap()) == 0);

        let d1 = ExactNumNumber::max_value(p).unwrap();
        assert!(d1.sinh_cosh(p, rm, &mut cc).unwrap_err() == Error::ExponentOverflow(Sign::Pos));

        let d2 = ExactNumNumber::min_value(p).unwrap();
        assert!(d2.sinh_cosh(p, rm, &mut cc).unwrap_err() == Error::ExponentOverflow(Sign::Neg));

        let zero = ExactNumNumber::new(1).unwrap();
        let (s, c) = zero.sinh_cosh(p, rm, &mut cc).unwrap();
        assert!(s.is_zero());
        assert!(c.cmp(&ExactNumNumber::from_word(1, p).unwrap()) == 0);

        let n1 = random_subnormal(p);
        let (s, c) = n1.sinh_cosh(p, rm, &mut cc).unwrap();
        assert!(s.cmp(&n1.sinh(p, rm, &mut cc).unwrap()) == 0);
        assert!(c.cmp(&n1.cosh(p, rm, &mut cc).unwrap()) == 0);
    }
}
