//! Interval enclosures (`Ball`) and Ziv correct-rounding retries.

use crate::common::util::bump_prec_retry;
use crate::common::util::round_p;
use crate::defs::WORD_BIT_SIZE;
use crate::Consts;
use crate::Error;
use crate::ExactNum;
use crate::RoundingMode;

/// Enclosure `mid ± rad` used to certify that a rounded midpoint is unique.
#[derive(Clone, Debug)]
pub struct Ball {
    mid: ExactNum,
    rad: ExactNum,
}

impl Ball {
    /// `mid ± rad`. The radius is taken in absolute value.
    pub fn new(mid: ExactNum, rad: ExactNum) -> Self {
        Ball {
            mid,
            rad: rad.abs(),
        }
    }

    /// Midpoint of the enclosure.
    pub fn mid(&self) -> &ExactNum {
        &self.mid
    }

    /// Non-negative radius.
    pub fn rad(&self) -> &ExactNum {
        &self.rad
    }

    fn rounding_ulp(x: &ExactNum, p: usize) -> ExactNum {
        if x.is_nan() || x.is_inf() || x.is_zero() {
            return ExactNum::new(p);
        }
        let e = x.exponent().unwrap_or(0);
        let bits = x.mantissa_max_bit_len().unwrap_or(p) as i32;
        let mut u = ExactNum::from_word(1, p);
        u.set_exponent(e.saturating_sub(bits.saturating_sub(2)));
        u
    }

    /// Sum of two balls: midpoint add at precision `p`, radius `r1+r2` plus a rounding ulp.
    pub fn add(&self, other: &Self, p: usize, rm: RoundingMode) -> Self {
        let mid = self.mid.add(&other.mid, p, rm);
        let rad = self.rad.add(&other.rad, p, RoundingMode::Up).add(
            &Self::rounding_ulp(&mid, p),
            p,
            RoundingMode::Up,
        );
        Ball { mid, rad }
    }

    /// Product of two balls with a first-order radius bound plus a rounding ulp.
    pub fn mul(&self, other: &Self, p: usize, rm: RoundingMode) -> Self {
        let mid = self.mid.mul(&other.mid, p, rm);
        let a = self.mid.abs().mul(&other.rad, p, RoundingMode::Up);
        let b = other.mid.abs().mul(&self.rad, p, RoundingMode::Up);
        let c = self.rad.mul(&other.rad, p, RoundingMode::Up);
        let rad = a
            .add(&b, p, RoundingMode::Up)
            .add(&c, p, RoundingMode::Up)
            .add(&Self::rounding_ulp(&mid, p), p, RoundingMode::Up);
        Ball { mid, rad }
    }

    /// Exponential of a ball. `exp` is increasing; the radius uses
    /// \(\lvert\exp(m)\rvert(e^{r}-1)\) plus a rounding ulp (same `Up` convention as `add`/`mul`).
    pub fn exp(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let mid = self.mid.exp(p, rm, cc);
        let em1 = self.rad.expm1(p, RoundingMode::Up, cc);
        let rad = mid.abs().mul(&em1, p, RoundingMode::Up).add(
            &Self::rounding_ulp(&mid, p),
            p,
            RoundingMode::Up,
        );
        Ball { mid, rad }
    }

    /// Sine of a ball. \(\lvert\sin'\rvert\le 1\), so the image radius is at most `rad`
    /// plus a rounding ulp.
    pub fn sin(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let mid = self.mid.sin(p, rm, cc);
        let rad = self
            .rad
            .add(&Self::rounding_ulp(&mid, p), p, RoundingMode::Up);
        Ball { mid, rad }
    }

    /// True when `x` lies in `[mid − rad, mid + rad]` (NaN / Inf never contained).
    pub fn contains(&self, x: &ExactNum, p: usize) -> bool {
        if x.is_nan() || self.mid.is_nan() || self.rad.is_nan() {
            return false;
        }
        let d = self.mid.sub(x, p, RoundingMode::None).abs();
        matches!(d.cmp(&self.rad), Some(c) if c <= 0)
    }
}

/// Evaluate `compute` at increasing working precision until the result rounds uniquely
/// to `p` bits (`try_set_precision`). Same retry budget as the transcendental kernel
/// ([`crate::MAX_PREC_RETRY`]).
pub fn ziv_round<F>(p: usize, rm: RoundingMode, mut compute: F) -> ExactNum
where
    F: FnMut(usize) -> ExactNum,
{
    let mut p_inc = WORD_BIT_SIZE;
    let mut p_wrk = match round_p(p).checked_add(p_inc) {
        Some(v) => v,
        None => return ExactNum::nan(Some(Error::InvalidArgument)),
    };
    loop {
        let mut v = compute(p_wrk);
        if v.try_set_precision(p, rm, p_wrk) {
            return v;
        }
        if bump_prec_retry(&mut p_wrk, &mut p_inc, p).is_err() {
            return ExactNum::nan(Some(Error::PrecisionRetryExhausted));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ball_add_contains_true_sum() {
        let p = 128;
        let rm = RoundingMode::ToEven;
        let a = ExactNum::from(3);
        let b = ExactNum::from(4);
        let u = ExactNum::from_word(1, p);
        let ba = Ball::new(a.clone(), u.clone());
        let bb = Ball::new(b.clone(), u.clone());
        let sum = ba.add(&bb, p, rm);
        let true_sum = a.add(&b, p, rm);
        assert!(sum.contains(&true_sum, p));
    }

    #[test]
    fn ziv_round_sqrt_matches_direct() {
        let p = 128;
        let rm = RoundingMode::ToEven;
        let two = ExactNum::from(2);
        let via_ziv = ziv_round(p, rm, |pw| two.sqrt(pw, RoundingMode::None));
        let direct = two.sqrt(p, rm);
        assert_eq!(via_ziv.cmp(&direct), Some(0));
    }

    #[test]
    fn ball_exp_contains_one_and_two() {
        let p = 128;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new().unwrap();
        let two = ExactNum::from_u8(2, p);
        let rad = two.powsi(-20, p, rm);
        let z = Ball::new(ExactNum::from_u8(0, p), rad.clone());
        let ez = z.exp(p, rm, &mut cc);
        assert!(ez.contains(&ExactNum::from_u8(1, p), p));

        let ln2 = cc.ln_2(p, rm);
        let bln = Ball::new(ln2, rad);
        let e2 = bln.exp(p, rm, &mut cc);
        assert!(e2.contains(&two, p));
    }

    #[test]
    fn ball_sin_contains_zero_at_origin_and_pi() {
        let p = 128;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new().unwrap();
        let two = ExactNum::from_u8(2, p);
        let rad = two.powsi(-20, p, rm);
        let z = Ball::new(ExactNum::from_u8(0, p), rad.clone());
        let sz = z.sin(p, rm, &mut cc);
        assert!(sz.contains(&ExactNum::from_u8(0, p), p));

        let pi = cc.pi(p, rm);
        let bpi = Ball::new(pi, rad);
        let sp = bpi.sin(p, rm, &mut cc);
        assert!(sp.contains(&ExactNum::from_u8(0, p), p));
    }

    #[test]
    fn ball_exp_sin_contain_scalar_at_a_point() {
        let p = 128;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new().unwrap();
        let x = ExactNum::from_u8(1, p);
        let rad = ExactNum::from_u8(2, p).powsi(-12, p, rm);
        let b = Ball::new(x.clone(), rad);
        let hi = x.exp(256, rm, &mut cc);
        assert!(b.exp(p, rm, &mut cc).contains(&hi, p));
        let hs = x.sin(256, rm, &mut cc);
        assert!(b.sin(p, rm, &mut cc).contains(&hs, p));
    }
}
