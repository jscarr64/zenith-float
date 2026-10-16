// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! Auxiliary items.

use crate::common::consts::TRIG_EXP_THRES;
use crate::common::util::round_p;
use crate::{num::ExactNumNumber, Consts, Error, Exponent, RoundingMode, Sign, WORD_BIT_SIZE};

impl ExactNumNumber {
    /// `±2^k` for `k >= EXPONENT_MIN - p` (`p` a multiple of the word size): the last bit of a `p`-bit number with the
    /// exponent `k + p`, including the subnormal range, where it is formed by an exact product of two powers of two.
    pub(crate) fn pow2_unit(sign: Sign, k: i64, p: usize) -> Result<Self, Error> {
        let min = crate::defs::EXPONENT_MIN as i64;
        if k >= min - 1 {
            return Self::sticky_term(sign, (k + 1) as Exponent);
        }
        let hi = Self::sticky_term(sign, crate::defs::EXPONENT_MIN)?;
        let lo = Self::sticky_term(Sign::Pos, (k - min + 2) as Exponent)?;
        hi.mul(&lo, p, RoundingMode::None)
    }

    /// A one-bit number `±2^(e - 1)`: the stand-in for a term that is far below the rounding position of a sum.
    /// Adding it exactly (`add_full_prec`) keeps the direction of the lost term, so that one final rounding
    /// by `set_precision` gives the correctly rounded result of the sum in every rounding mode.
    pub(crate) fn sticky_term(sign: Sign, e: Exponent) -> Result<Self, Error> {
        let mut t = Self::from_word(1, WORD_BIT_SIZE)?;
        t.set_exponent(e);
        t.set_sign(sign);
        Ok(t)
    }

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
