//! Paired sine and cosine.

use crate::common::consts::TRIG_EXP_THRES;
use crate::common::util::bump_prec_retry;
use crate::common::util::round_p;
use crate::defs::Error;
use crate::defs::RoundingMode;
use crate::num::ExactNumNumber;
use crate::ops::consts::Consts;
use crate::WORD_BIT_SIZE;

impl ExactNumNumber {
    /// Computes `(sin(self), cos(self))` with precision `p`.
    ///
    /// Argument reduction is shared; each output is rounded independently at `p`.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory.
    ///  - InvalidArgument: the precision is incorrect.
    pub fn sin_cos(
        &self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<(Self, Self), Error> {
        let p = round_p(p);

        if self.is_zero() {
            let sin = Self::new2(p, self.sign(), self.inexact())?;
            let mut cos = Self::from_word(1, p)?;
            cos.set_inexact(self.inexact());
            return Ok((sin, cos));
        }

        let mut p_inc = WORD_BIT_SIZE;
        // The retry budget counts from the working start: an argument wider than `p` must not exhaust it at once.
        let p_base = p.max(self.mantissa_max_bit_len());
        let mut p_wrk = p_base;

        if (p_wrk as i128) + 1 < -(i128::from(self.exponent()) * 2 - 2) {
            let mut x = self.clone()?;
            if p > x.mantissa_max_bit_len() {
                x.set_precision(p, RoundingMode::None)?;
            }
            let mut sn = x.add_correction(true)?;
            sn.set_precision(p, rm)?;
            let mut cs = Self::from_word(1, p)?;
            cs.set_inexact(self.inexact());
            return Ok((sn, cs));
        }

        p_wrk += p_inc;

        let mut add_p = (3 - TRIG_EXP_THRES) as usize;
        loop {
            let mut x = self.clone()?;

            let p_x = p_wrk + add_p;
            x.set_precision(p_x, RoundingMode::None)?;

            x = x.reduce_trig_arg(cc, RoundingMode::None)?;

            let (t, _q) = x.trig_arg_pi_proximity(cc, RoundingMode::None)?;
            if add_p < t {
                add_p = t;
            } else {
                let mut sn = x.clone()?.sin_series(RoundingMode::None)?;
                let mut cs = x.cos_series(RoundingMode::None)?;

                let sn_ok = sn.try_set_precision(p, rm, p_wrk)?;
                let cs_ok = cs.try_set_precision(p, rm, p_wrk)?;

                if sn_ok && cs_ok {
                    sn.set_inexact(sn.inexact() | self.inexact());
                    cs.set_inexact(cs.inexact() | self.inexact());
                    break Ok((sn, cs));
                }

                bump_prec_retry(&mut p_wrk, &mut p_inc, p_base)?;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::util::random_subnormal;

    #[test]
    fn test_sin_cos() -> Result<(), crate::Error> {
        let p = 320;
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;

        let d1 = ExactNumNumber::from_word(1, p)?;
        let (s, c) = d1.sin_cos(p, rm, &mut cc)?;
        assert!(s.cmp(&d1.sin(p, rm, &mut cc)?) == 0);
        assert!(c.cmp(&d1.cos(p, rm, &mut cc)?) == 0);

        let zero = ExactNumNumber::new(1)?;
        let (s, c) = zero.sin_cos(p, rm, &mut cc)?;
        assert!(s.is_zero());
        assert!(c.cmp(&ExactNumNumber::from_word(1, p)?) == 0);

        let n1 = random_subnormal(p)?;
        let (s, c) = n1.sin_cos(p, rm, &mut cc)?;
        assert!(s.cmp(&n1.sin(p, rm, &mut cc)?) == 0);
        assert!(c.cmp(&n1.cos(p, rm, &mut cc)?) == 0);
        Ok(())
    }
}
