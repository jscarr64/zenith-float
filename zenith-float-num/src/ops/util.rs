// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! Auxiliary items.

use crate::common::consts::TRIG_EXP_THRES;
use crate::common::util::round_p;
use crate::{num::ExactNumNumber, Consts, Error, RoundingMode};

impl ExactNumNumber {
    /// Reduce `self` to the interval `(-2π, 2π)` for trigonometric evaluation.
    ///
    /// Uses the same working precision as `sin`/`cos`/`tan` before modular reduction.
    pub fn rem_pi(&self, p: usize, _rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        if p < crate::defs::MIN_PRECISION {
            return Err(Error::InvalidArgument);
        }
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
            let extra = trig_reduce_extra(self.exponent());
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
        if i128::from($p_wrk as i64) + 1 < -i128::from($exp) {
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

/// Extra bits of pi for reducing an argument with exponent `e` (`e > 2`), capped at 2^16.
fn trig_reduce_extra(e: crate::Exponent) -> usize {
    // `e > 2`; an exponent past `usize` saturates to the cap instead of truncating on 32-bit.
    usize::try_from(e).map_or(1 << 16, |x| x.min(1 << 16))
}

#[cfg(test)]
mod tests {
    use super::trig_reduce_extra;

    /// An exponent past `u32` must hit the 2^16 cap; a 32-bit `usize` truncation would give 5.
    #[test]
    fn trig_reduce_extra_exponent_beyond_u32() {
        assert_eq!(trig_reduce_extra((1i64 << 32) + 5), 1 << 16);
        assert_eq!(trig_reduce_extra(3), 3);
        assert_eq!(trig_reduce_extra(1i64 << 16), 1 << 16);
        assert_eq!(trig_reduce_extra(crate::EXPONENT_MAX), 1 << 16);
    }
}
