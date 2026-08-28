//! Auxiliary items.

use crate::common::consts::TRIG_EXP_THRES;
use crate::common::util::round_p;
use crate::{num::ExactNumNumber, Consts, Error, RoundingMode};

impl ExactNumNumber {
    /// Reduce `self` to the interval `(-2π, 2π)` for trigonometric evaluation.
    ///
    /// Uses the same working precision as `sin`/`cos`/`tan` before modular reduction.
    pub fn rem_pi(&self, p: usize, _rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        let p = round_p(p);
        let add_p = (3 - TRIG_EXP_THRES) as usize;
        let p_x = p.max(self.mantissa_max_bit_len()).saturating_add(add_p);
        let mut x = self.clone()?;
        x.set_precision(p_x, RoundingMode::None)?;
        x.reduce_trig_arg(cc, RoundingMode::None)
    }

    /// Reduce `self` to interval (-2*pi; 2*pi)
    pub(crate) fn reduce_trig_arg(self, cc: &mut Consts, rm: RoundingMode) -> Result<Self, Error> {
        if self.exponent() > 2 {
            let extra = (self.exponent() as usize).min(1 << 16);
            let mut pi = cc.pi_num(self.mantissa_max_bit_len().saturating_add(extra), rm)?;

            pi.set_exponent(pi.exponent() + 1);

            self.rem(&pi)
        } else {
            Ok(self)
        }
    }

    /// Determine how close `self` to pi or pi/2
    pub(crate) fn trig_arg_pi_proximity(
        &self,
        cc: &mut Consts,
        rm: RoundingMode,
    ) -> Result<(usize, usize), Error> {
        let mut q = 0;

        if self.exponent() > 0 {
            let mut pi = cc.pi_num(self.mantissa_max_bit_len() + self.exponent() as usize, rm)?;
            pi.set_sign(self.sign());

            if self.exponent() < 2 {
                pi.set_exponent(1);
            }

            let mut ret = self.sub(&pi, self.mantissa_max_bit_len(), rm)?;
            q += pi.exponent() as usize;

            if ret.exponent() > 0 && self.sign() == ret.sign() {
                if ret.exponent() < 2 {
                    pi.set_exponent(1);
                }

                ret = ret.sub(&pi, self.mantissa_max_bit_len(), rm)?;
                q += pi.exponent() as usize;
            }

            Ok((ret.exponent().unsigned_abs() as usize, q))
        } else {
            Ok((self.exponent().unsigned_abs() as usize, q))
        }
    }
}

/// Compute result for argument with small exponent.
macro_rules! compute_small_exp {
    ($arg:ident, $exp:expr, $sign_inverse:expr, $p_wrk:ident, $p:ident, $rm:ident) => {
        if $p_wrk as isize + 1 < -($exp) {
            let mut x = $arg.clone()?;
            if $p > x.mantissa_max_bit_len() {
                x.set_precision($p, RoundingMode::None)?;
            }
            let mut ret = x.add_correction($sign_inverse)?;
            ret.set_precision($p, $rm)?;
            return Ok(ret);
        }
    };
}

pub(super) use compute_small_exp;
