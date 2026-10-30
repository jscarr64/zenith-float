// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! Auxiliary items.

use crate::common::util::round_p_rm;
use crate::{num::ExactNumNumber, Consts, Error, RoundingMode};

/// Initial guard bits of the `rem_pi` reduction (doubled on each Ziv retry); the same on every target.
const REM_PI_GUARD: usize = 64;

impl ExactNumNumber {
    /// Reduce `self` to the interval `(-2π, 2π)`: `x - trunc(x / 2π) * 2π`, with the sign of `x`.
    ///
    /// The result is rounded at `p` bits with `rm`; `RoundingMode::None` truncates at
    /// `rp(p, 64)` bits, so the value is the same on every target. The reduction uses `2π` with
    /// `exp(x) + p + guard` bits: `x` has fewer than `2^(exp(x) - 2)` periods, and each carries an
    /// error below `2^(3 - prec(2π))`, so the absolute error is below `2^(exp(x) + 1 - prec(2π))`.
    /// The bits lost to cancellation, `exp(x) - exp(result)`, are charged against that bound and the
    /// guard is doubled (Ziv) until the result rounds correctly.
    pub fn rem_pi(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        if p < crate::defs::MIN_PRECISION {
            return Err(Error::InvalidArgument);
        }
        // `None` is a truncation at whole 64-bit words; truncation toward zero is `ToZero`.
        let (p, rm) = if rm == RoundingMode::None {
            (round_p_rm(p, rm), RoundingMode::ToZero)
        } else {
            (p, rm)
        };
        let ex = self.exponent();
        if self.is_zero() || ex <= 2 {
            // |x| < 4 < 2π: the reduction is the identity.
            let mut x = self.clone()?;
            x.set_precision(p, rm)?;
            return Ok(x);
        }
        let ex_bits = usize::try_from(ex).map_err(|_| Error::InvalidArgument)?;
        let mut guard = REM_PI_GUARD;
        loop {
            let pi_p = self
                .mantissa_max_bit_len()
                .max(p)
                .checked_add(ex_bits)
                .and_then(|v| v.checked_add(guard))
                .ok_or(Error::InvalidArgument)?;
            let mut two_pi = cc.pi_num(pi_p, RoundingMode::None)?;
            let e2 = two_pi
                .exponent()
                .checked_add(1)
                .ok_or(Error::InvalidArgument)?;
            two_pi.set_exponent(e2);
            let mut r = self.rem(&two_pi)?;
            if !r.is_zero() {
                // correct leading bits: exp(r) - (exp(x) + 1 - prec(2π))
                let prec = i64::try_from(two_pi.mantissa_max_bit_len())
                    .map_err(|_| Error::InvalidArgument)?;
                let s = r
                    .exponent()
                    .checked_add(prec)
                    .and_then(|v| v.checked_sub(ex))
                    .and_then(|v| v.checked_sub(1))
                    .ok_or(Error::InvalidArgument)?;
                let s = usize::try_from(s)
                    .unwrap_or(0)
                    .min(r.mantissa_max_bit_len());
                if s > p && r.try_set_precision(p, rm, s)? {
                    return Ok(r);
                }
            }
            guard = guard.checked_mul(2).ok_or(Error::InvalidArgument)?;
        }
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
