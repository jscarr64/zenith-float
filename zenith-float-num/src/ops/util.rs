// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! Auxiliary items.

use crate::common::consts::TRIG_EXP_THRES;
use crate::common::util::{reduction_width, round_p};
use crate::{num::ExactNumNumber, Consts, Error, RoundingMode};

impl ExactNumNumber {
    /// Reduce `self` modulo `2π` to the interval `(-2π, 2π)` (the sign of `self`), correctly rounded to `p`
    /// bits with `rm`. π carries `exponent(self)` extra bits, and the width grows until the remainder's own
    /// leading bits (it may lie close to a multiple of 2π) are certain.
    pub fn rem_pi(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        if p < crate::defs::MIN_PRECISION {
            return Err(Error::InvalidArgument);
        }
        let p = round_p(p);
        // The reduction is the identity, exact, below 2π: |self| < 4, or |self| in [4, 8) (exponent 3) below 2π
        // (π rounded down, so 2π here <= the true 2π). Dividing by the inexact 2π would flag the unchanged value
        // inexact, and a directed mode could then never round its all-zero tail (the Ziv loop below would widen π
        // without end).
        let identity = self.is_zero()
            || self.exponent() <= 2
            || (self.exponent() == 3 && {
                let w = p.max(self.mantissa_max_bit_len()) + crate::defs::WORD_BIT_SIZE;
                let mut two_pi = cc.pi_num(w, RoundingMode::Down)?;
                two_pi.set_exponent(two_pi.exponent() + 1);
                self.abs_cmp(&two_pi) < 0
            });
        if identity {
            return self.try_rounded_to(p, rm);
        }
        let add_p = (3 - TRIG_EXP_THRES) as usize;
        let mut p_x = p.max(self.mantissa_max_bit_len()).saturating_add(add_p);
        loop {
            let x = self.try_rounded_to(p_x, RoundingMode::None)?;
            let mut r = x.reduce_trig_arg(cc, RoundingMode::None)?;
            if r.is_zero() {
                p_x = round_p(p_x.saturating_mul(2));
                continue;
            }
            // |error| <= 2^(4 - p_x) absolute (|self| / 2π multiples of a p_x + e(self) bit 2π).
            let lost = (4 - r.exponent()).max(0) as usize;
            if let Some(s) = p_x.checked_sub(lost + 4) {
                let s = s / crate::defs::WORD_BIT_SIZE * crate::defs::WORD_BIT_SIZE;
                if s > p && r.try_set_precision(p, rm, s)? {
                    return Ok(r);
                }
            }
            p_x = round_p(
                p_x.saturating_add(lost)
                    .saturating_add(p_x / 4)
                    .saturating_add(64),
            );
        }
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
            let mut ret = $arg.add_correction_at($p, $sign_inverse)?;
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

    /// |x| < 2π (exponent 3: |x| in [4, 8)) reduces to x itself, exactly. The remainder by the
    /// inexact 2π used to come back flagged inexact with an all-zero tail, so a directed mode could
    /// never round it and the Ziv loop widened π without end (the full-tier MPFR `mpfr_compare_ops`
    /// hang, r6: x = -4.937..., p = 3392, Up). Gold: x itself (the definition), not inexact.
    #[cfg(feature = "std")]
    #[test]
    fn rem_pi_below_two_pi_is_the_identity_in_every_mode() -> Result<(), String> {
        use crate::{Consts, ExactNum, Exponent, RoundingMode};
        let p = 3392usize;
        // -(79/16) * (1 + 2^-3000) = -4.9375...: below 2π, exponent 3, a long mantissa
        let m = ExactNum::from_u8(1, p).add(
            &ExactNum::from_u8(1, p).ldexp(-3000 as Exponent, p, RoundingMode::ToEven),
            p,
            RoundingMode::ToEven,
        );
        let mut x = -ExactNum::from_u8(79, p)
            .ldexp(-4, p, RoundingMode::ToEven)
            .mul(&m, p, RoundingMode::ToEven);
        for rm in [
            RoundingMode::ToEven,
            RoundingMode::Up,
            RoundingMode::Down,
            RoundingMode::ToZero,
            RoundingMode::FromZero,
        ] {
            // x moves into the worker (a hang must not hang the test) and comes back with the result
            let (tx, rx) = std::sync::mpsc::channel();
            std::thread::spawn(move || {
                let r = Consts::new().map(|mut cc| ExactNum::rem_pi(&x, p, rm, &mut cc));
                let _sent = tx.send((r, x));
            });
            let (r, back) = rx
                .recv_timeout(std::time::Duration::from_secs(60))
                .map_err(|_| format!("rem_pi {rm:?}: no result in 60 s (Ziv loop does not end)"))?;
            x = back;
            let r = r.map_err(|e| format!("{e:?}"))?;
            if r != x || r.inexact() {
                return Err(format!(
                    "rem_pi {rm:?}: got {r:?} (inexact {}), want x exactly",
                    r.inexact()
                ));
            }
        }
        Ok(())
    }

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
