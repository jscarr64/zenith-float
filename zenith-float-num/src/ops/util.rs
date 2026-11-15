// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! Auxiliary items.

use crate::common::consts::TRIG_EXP_THRES;
use crate::common::util::{reduction_width, round_p};
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
            // π carries exponent(self) extra bits, so the remainder keeps every bit of `self`;
            // the width is checked (`MemoryAllocation` past what the machine holds), never capped.
            let width = reduction_width(self.mantissa_max_bit_len(), self.exponent())?;
            let mut pi = cc.pi_num(width, rm)?;

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

#[cfg(test)]
mod tests {
    use crate::common::util::reduction_width;
    use crate::Error;

    /// The π width follows the argument exponent with no limit: 2^16 (main's old silent cap),
    /// 2^22 (Tier 1's old refusal) and one past each are all carried in full. Only a width the
    /// exponent arithmetic cannot hold is refused, and as `MemoryAllocation`, not a precision limit.
    #[test]
    fn trig_reduction_width_has_no_cap() {
        for e in [
            3i64,
            (1 << 16) - 1,
            1 << 16,
            (1 << 16) + 1,
            (1 << 22) - 1,
            1 << 22,
            (1 << 22) + 1,
            (1i64 << 32) + 5,
        ] {
            assert_eq!(reduction_width(128, e), Ok(128 + e as usize), "e = {e}");
        }
        assert_eq!(
            reduction_width(128, crate::EXPONENT_MAX),
            Err(Error::MemoryAllocation)
        );
        assert_eq!(reduction_width(usize::MAX, 1), Err(Error::MemoryAllocation));
    }
}
