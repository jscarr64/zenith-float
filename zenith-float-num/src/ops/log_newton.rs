//! Newton logarithm, after FLINT 3.6 `arb_log_newton` (Fredrik Johansson,
//! LGPL-3.0-or-later).
//!
//! Below [`ARB_LOG_NEWTON_PREC`] `/ 2` bits the existing atanh series is
//! used (and is the recursion floor, so Newton cannot loop). Above that,
//! a low-precision logarithm is exponentiated once and a short integer
//! polynomial corrects `log(1+w)` with `w = x e^{-t} - 1`. The binary
//! exponent is factored by the caller as `n log 2` from the cached
//! constant, as in the existing `ln` factoring path. Arguments near 1
//! use the same short polynomial without a Newton step.

use crate::common::consts::ONE;
use crate::defs::Error;
use crate::defs::RoundingMode;
use crate::defs::Sign;
use crate::num::ExactNumNumber;
use crate::ops::consts::Consts;
use crate::WORD_BIT_SIZE;

/// FLINT 3.6 `ARB_LOG_NEWTON_PREC`. Newton is used above half of this
/// (1400 bits) so the seed logarithm never recurses into itself.
pub(crate) const ARB_LOG_NEWTON_PREC: usize = 2800;

/// Numerator coefficients of the truncated `log(1+w)` polynomial, and the
/// common denominator 232792560, from FLINT 3.6 `log_newton.c`.
const LOG_COEFFS: [i32; 21] = [
    0, 232792560, -116396280, 77597520, -58198140, 46558512, -38798760, 33256080, -29099070,
    25865840, -23279256, 21162960, -19399380, 17907120, -16628040, 15519504, -14549535, 13693680,
    -12932920, 12252240, -11639628,
];
const LOG_DEN: u64 = 232792560;

fn coeff_term(c: i32, p: usize) -> Result<ExactNumNumber, Error> {
    let mut v = ExactNumNumber::from_word(c.unsigned_abs() as u64, p)?;
    if c < 0 {
        v.set_sign(Sign::Neg);
    }
    Ok(v)
}

/// Horner evaluation of `sum_{k=0}^{n-1} c_k w^k / LOG_DEN`.
fn eval_log_poly(w: &ExactNumNumber, n: usize, p: usize) -> Result<ExactNumNumber, Error> {
    let n = n.clamp(1, LOG_COEFFS.len());
    let mut acc = coeff_term(LOG_COEFFS[n - 1], p)?;
    for k in (0..n - 1).rev() {
        acc = acc.mul(w, p, RoundingMode::None)?;
        if LOG_COEFFS[k] != 0 {
            let c = coeff_term(LOG_COEFFS[k], p)?;
            acc = acc.add(&c, p, RoundingMode::None)?;
        }
    }
    let den = ExactNumNumber::from_u64_internal(LOG_DEN, p)?;
    acc.div(&den, p, RoundingMode::None)
}

fn newton_term_count(prec: usize) -> usize {
    let n = if prec <= 3200 {
        4
    } else if prec <= 6000 {
        6
    } else if prec <= 300_000 {
        7
    } else if prec <= 1_000_000 {
        9
    } else {
        12
    };
    n.clamp(2, 16)
}

/// Magnitude bound `k` with `|x| < 2^k`, equal to `x.exponent()` for a
/// nonzero finite `x`. Zero is treated as a very negative magnitude.
fn abs_bound_lt_2exp(x: &ExactNumNumber) -> i64 {
    if x.is_zero() {
        i64::MIN / 4
    } else {
        x.exponent()
    }
}

impl ExactNumNumber {
    /// `|self - 1| < 2^{-p/16}`, the FLINT near-1 shortcut threshold.
    pub(crate) fn ln_is_near_one(&self, p: usize) -> Result<bool, Error> {
        let wp = p.saturating_add(8).max(WORD_BIT_SIZE);
        let d = self.sub(&ONE, wp, RoundingMode::None)?;
        let x1mag = abs_bound_lt_2exp(&d);
        Ok(x1mag < -((p as i64) / 16))
    }

    /// Short `log(1+w)` polynomial for an argument near 1, after
    /// `arb_log_newton`.
    pub(crate) fn ln_near_one(&self, p: usize) -> Result<Self, Error> {
        let wp = p.saturating_add(16).max(WORD_BIT_SIZE);
        let t = self.sub(&ONE, wp, RoundingMode::None)?;
        let x1mag = abs_bound_lt_2exp(&t);
        if t.is_zero() {
            return Self::from_word(0, p);
        }
        let n = if x1mag >= 0 {
            16
        } else {
            let den = (-x1mag) as usize;
            (p / den) + 2
        };
        eval_log_poly(&t, n, wp)
    }

    /// Newton logarithm of a positive argument that is not 1, at working
    /// precision `p`. The argument should already be reduced to
    /// `[1/2, 1)` by the caller when the binary exponent is large; this
    /// routine still accepts a general positive `x`.
    pub(crate) fn ln_newton(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        if self.is_zero() || self.is_negative() {
            return Err(Error::InvalidArgument);
        }

        let d = self.sub(&ONE, p.saturating_add(8).max(WORD_BIT_SIZE), RoundingMode::None)?;
        let x1mag = abs_bound_lt_2exp(&d);

        if d.is_zero() {
            return Self::from_word(0, p);
        }

        if x1mag < -((p as i64) / 16) {
            let mut ret = self.ln_near_one(p)?;
            ret.set_precision(p, RoundingMode::None)?;
            return Ok(ret);
        }

        if p <= ARB_LOG_NEWTON_PREC / 2 {
            let mut x = self.clone()?;
            x.set_precision(p, RoundingMode::None)?;
            return Self::ln_series(x, RoundingMode::None);
        }

        let n = newton_term_count(p);
        let extra = 10usize;
        let mut wp = p + 10;
        if x1mag < 0 {
            wp += (-x1mag) as usize;
        }
        let wp2 = wp * (n - 1) / n;
        let seed_p = (wp / n).saturating_add(extra);

        let t = if seed_p <= ARB_LOG_NEWTON_PREC / 2 {
            let mut x = self.clone()?;
            x.set_precision(seed_p, RoundingMode::None)?;
            Self::ln_series(x, RoundingMode::None)?
        } else {
            self.ln_newton(seed_p, cc)?
        };

        // w = x * exp(-t) - 1
        let mut neg_t = t.clone()?;
        if !t.is_zero() {
            neg_t.set_sign(if t.is_negative() { Sign::Pos } else { Sign::Neg });
        }
        // exp of a seed-sized argument at wp bits: uses bit-burst when wp is large.
        let e = neg_t.exp(wp, RoundingMode::None, cc)?;
        let mut xw = self.clone()?;
        xw.set_precision(wp, RoundingMode::None)?;
        let w = xw.mul(&e, wp, RoundingMode::None)?.sub(&ONE, wp2, RoundingMode::None)?;

        let corr = eval_log_poly(&w, n, wp2)?;
        let mut ret = t.add(&corr, p, RoundingMode::None)?;
        ret.set_precision(p, RoundingMode::None)?;
        Ok(ret)
    }
}
