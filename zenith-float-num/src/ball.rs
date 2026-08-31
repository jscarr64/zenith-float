//! Interval enclosures (`Ball`) and Ziv correct-rounding retries.

use crate::common::util::bump_prec_retry;
use crate::common::util::round_p;
use crate::defs::WORD_BIT_SIZE;
use crate::Consts;
use crate::Error;
use crate::ExactNum;
use crate::RoundingMode;

/// Extra rounding-ulp multiples added to a transcendental Lipschitz radius.
const BALL_TRANSCENDENTAL_ERROR_TERMS: u32 = 8;

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

    /// Cosine of a ball. \(\lvert\cos'\rvert\le 1\).
    pub fn cos(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let mid = self.mid.cos(p, rm, cc);
        let rad = self
            .rad
            .add(&Self::rounding_ulp(&mid, p), p, RoundingMode::Up);
        Ball { mid, rad }
    }

    /// Natural log of a ball. Domain: the ball must lie in \((0,+\infty)\).
    /// Lipschitz \(\lvert\ln'\rvert=1/x\le 1/(m-r)\).
    pub fn ln(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !self.strictly_positive(p) {
            return Self::nan_ball(p);
        }
        let mid = self.mid.ln(p, rm, cc);
        let den = self.mid.sub(&self.rad, p, RoundingMode::Down);
        let lip = ExactNum::from_u8(1, p).div(&den, p, RoundingMode::Up);
        let rad = lip.mul(&self.rad, p, RoundingMode::Up).add(
            &Self::transcendental_slack(&mid, p),
            p,
            RoundingMode::Up,
        );
        Ball { mid, rad }
    }

    /// Square root of a ball. Domain: the ball must lie in \((0,+\infty)\).
    /// Lipschitz \(1/(2\sqrt{x})\le 1/(2\sqrt{m-r})\).
    pub fn sqrt(&self, p: usize, rm: RoundingMode) -> Self {
        if !self.strictly_positive(p) {
            return Self::nan_ball(p);
        }
        let mid = self.mid.sqrt(p, rm);
        let lo = self
            .mid
            .sub(&self.rad, p, RoundingMode::Down)
            .sqrt(p, RoundingMode::Down);
        let two = ExactNum::from_u8(2, p);
        let den = two.mul(&lo, p, RoundingMode::Down);
        let lip = ExactNum::from_u8(1, p).div(&den, p, RoundingMode::Up);
        let rad = lip.mul(&self.rad, p, RoundingMode::Up).add(
            &Self::transcendental_slack(&mid, p),
            p,
            RoundingMode::Up,
        );
        Ball { mid, rad }
    }

    /// Error function of a ball. \(\lvert\mathrm{erf}'\rvert\le 2/\sqrt{\pi}\).
    pub fn erf(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let mid = self.mid.erf(p, rm, cc);
        let two = ExactNum::from_u8(2, p);
        let s = cc.pi(p, RoundingMode::Down).sqrt(p, RoundingMode::Down);
        let lip = two.div(&s, p, RoundingMode::Up);
        let rad = lip.mul(&self.rad, p, RoundingMode::Up).add(
            &Self::transcendental_slack(&mid, p),
            p,
            RoundingMode::Up,
        );
        Ball { mid, rad }
    }

    /// \(J_0\) of a ball. \(\lvert J_0'\rvert=\lvert J_1\rvert\le 1\).
    pub fn bessel_j0(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let mid = self.mid.bessel_j(0, p, rm, cc);
        let rad = self
            .rad
            .add(&Self::transcendental_slack(&mid, p), p, RoundingMode::Up);
        Ball { mid, rad }
    }

    /// \(J_1\) of a ball. \(\lvert J_1'\rvert=\lvert J_0-J_1/x\rvert\le 1+1/(\lvert m\rvert-r)\)
    /// when the ball excludes \(0\).
    pub fn bessel_j1(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !self.excludes_zero(p) {
            return Self::nan_ball(p);
        }
        let mid = self.mid.bessel_j(1, p, rm, cc);
        let den = self.mid.abs().sub(&self.rad, p, RoundingMode::Down);
        let extra = ExactNum::from_u8(1, p).div(&den, p, RoundingMode::Up);
        let lip = ExactNum::from_u8(1, p).add(&extra, p, RoundingMode::Up);
        let rad = lip.mul(&self.rad, p, RoundingMode::Up).add(
            &Self::transcendental_slack(&mid, p),
            p,
            RoundingMode::Up,
        );
        Ball { mid, rad }
    }

    fn transcendental_slack(mid: &ExactNum, p: usize) -> ExactNum {
        let u = Self::rounding_ulp(mid, p);
        u.mul(
            &ExactNum::from_u8(BALL_TRANSCENDENTAL_ERROR_TERMS as u8, p),
            p,
            RoundingMode::Up,
        )
    }

    fn nan_ball(p: usize) -> Self {
        let n = ExactNum::nan(Some(Error::InvalidArgument));
        let _ = p;
        Ball {
            mid: n.clone(),
            rad: n,
        }
    }

    fn strictly_positive(&self, p: usize) -> bool {
        if self.mid.is_nan() || self.rad.is_nan() || !self.mid.is_positive() {
            return false;
        }
        matches!(self.mid.cmp(&self.rad), Some(c) if c > 0)
            && !self.mid.sub(&self.rad, p, RoundingMode::Down).is_negative()
            && !self.mid.sub(&self.rad, p, RoundingMode::Down).is_zero()
    }

    fn excludes_zero(&self, p: usize) -> bool {
        if self.mid.is_nan() || self.rad.is_nan() || self.mid.is_zero() {
            return false;
        }
        matches!(self.mid.abs().cmp(&self.rad), Some(c) if c > 0)
            && !self
                .mid
                .abs()
                .sub(&self.rad, p, RoundingMode::Down)
                .is_zero()
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

    #[test]
    fn ball_plan_transcendental_golds() {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new().unwrap();
        let rad = ExactNum::from_u8(1, p).ldexp(-(p as i32), p, RoundingMode::None);
        let six = ExactNum::from_u8(6, p);
        let pi6 = cc.pi(p, rm).div(&six, p, rm);
        let bsin = Ball::new(pi6, rad.clone());
        let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, rm);
        assert!(bsin.sin(p, rm, &mut cc).contains(&half, p));

        let one = ExactNum::from_u8(1, p);
        let be = Ball::new(one.clone(), rad.clone());
        let e = cc.e(p, rm);
        assert!(be.exp(p, rm, &mut cc).contains(&e, p));

        let erf1 = one.erf(p, rm, &mut cc);
        assert!(be.erf(p, rm, &mut cc).contains(&erf1, p));

        let ln1 = one.ln(256, rm, &mut cc);
        assert!(be.ln(p, rm, &mut cc).contains(&ln1, p));
        let sq = one.sqrt(256, rm);
        assert!(be.sqrt(p, rm).contains(&sq, p));
        assert!(be
            .cos(p, rm, &mut cc)
            .contains(&one.cos(256, rm, &mut cc), p));
        assert!(be
            .bessel_j0(p, rm, &mut cc)
            .contains(&one.bessel_j(0, 256, rm, &mut cc), p));
        assert!(be
            .bessel_j1(p, rm, &mut cc)
            .contains(&one.bessel_j(1, 256, rm, &mut cc), p));

        let composed = be.sin(p, rm, &mut cc).exp(p, rm, &mut cc);
        let true_c = one.sin(256, rm, &mut cc).exp(256, rm, &mut cc);
        assert!(composed.contains(&true_c, p));
    }
}
