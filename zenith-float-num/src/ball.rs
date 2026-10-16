//! Interval enclosures (`Ball`) and Ziv correct-rounding retries.

use crate::common::util::bump_prec_retry;
use crate::common::util::round_p;
use crate::defs::WORD_BIT_SIZE;
use crate::Consts;
use crate::Error;
use crate::ExactComplex;
use crate::ExactNum;
use crate::RoundingMode;
use core::ops::Neg;

/// Extra rounding-ulp multiples added to a transcendental Lipschitz radius.
const BALL_TRANSCENDENTAL_ERROR_TERMS: u32 = 8;

/// Enclosure `mid ± rad` used to certify that a rounded midpoint is unique.
#[derive(Debug)]
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
    ///
    /// Accuracy: enclosure
    pub fn rad(&self) -> &ExactNum {
        &self.rad
    }

    fn rounding_ulp(x: &ExactNum, p: usize) -> ExactNum {
        if x.is_nan() || x.is_inf() || x.is_zero() {
            return ExactNum::new(p);
        }
        let e = x.exponent().unwrap_or(0);
        // The midpoint is rounded at `p` bits (2.0.0: exact-bit rounding), not at its storage width.
        let bits = x.mantissa_max_bit_len().unwrap_or(p).min(p.max(1)) as i32;
        let mut u = ExactNum::from_word(1, p);
        u.set_exponent(e.saturating_sub(bits.saturating_sub(2)));
        u
    }

    /// Sum of two balls: midpoint add at precision `p`, radius `r1+r2` plus a rounding ulp.
    pub fn add(self, other: &Self, p: usize, rm: RoundingMode) -> Self {
        let Ball { mid: m1, rad: r1 } = self;
        let mid = m1.add(&other.mid, p, rm);
        let rad = r1.add(&other.rad, p, RoundingMode::Up).add(
            &Self::rounding_ulp(&mid, p),
            p,
            RoundingMode::Up,
        );
        Ball { mid, rad }
    }

    /// Product of two balls with a first-order radius bound plus a rounding ulp.
    pub fn mul(&self, other: &Self, p: usize, rm: RoundingMode) -> Self {
        let mid = ExactNum::mul_ref(&self.mid, &other.mid, p, rm);
        let a = if self.mid.is_negative() {
            self.mid.neg_via_mul().mul(&other.rad, p, RoundingMode::Up)
        } else {
            ExactNum::mul_ref(&self.mid, &other.rad, p, RoundingMode::Up)
        };
        let b = if other.mid.is_negative() {
            other.mid.neg_via_mul().mul(&self.rad, p, RoundingMode::Up)
        } else {
            ExactNum::mul_ref(&other.mid, &self.rad, p, RoundingMode::Up)
        };
        let c = ExactNum::mul_ref(&self.rad, &other.rad, p, RoundingMode::Up);
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
        let rad = (if mid.is_negative() {
            mid.neg_via_mul().mul(&em1, p, RoundingMode::Up)
        } else {
            em1.mul(&mid, p, RoundingMode::Up)
        })
        .add(&Self::rounding_ulp(&mid, p), p, RoundingMode::Up);
        Ball { mid, rad }
    }

    /// Sine of a ball. \(\lvert\sin'\rvert\le 1\), so the image radius is at most `rad`
    /// plus a rounding ulp.
    pub fn sin(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let mid = self.mid.sin(p, rm, cc);
        let rad = Self::rounding_ulp(&mid, p).add(&self.rad, p, RoundingMode::Up);
        Ball { mid, rad }
    }

    /// Cosine of a ball. \(\lvert\cos'\rvert\le 1\).
    pub fn cos(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let mid = self.mid.cos(p, rm, cc);
        let rad = Self::rounding_ulp(&mid, p).add(&self.rad, p, RoundingMode::Up);
        Ball { mid, rad }
    }

    /// Natural log of a ball. Domain: the ball must lie in \((0,+\infty)\).
    /// Lipschitz \(\lvert\ln'\rvert=1/x\le 1/(m-r)\).
    pub fn ln(self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let Ball { mid: m, rad: r } = self;
        if !Self::above_zero(&m, &r) {
            return Self::nan_ball(p);
        }
        let mid = m.ln(p, rm, cc);
        let den = m.sub(&r, p, RoundingMode::Down);
        if den.is_negative() || den.is_zero() {
            return Self::nan_ball(p);
        }
        let lip = ExactNum::from_u8(1, p).div(&den, p, RoundingMode::Up);
        let rad = lip.mul(&r, p, RoundingMode::Up).add(
            &Self::transcendental_slack(&mid, p),
            p,
            RoundingMode::Up,
        );
        Ball { mid, rad }
    }

    /// Square root of a ball. Domain: the ball must lie in \((0,+\infty)\).
    /// Lipschitz \(1/(2\sqrt{x})\le 1/(2\sqrt{m-r})\).
    pub fn sqrt(self, p: usize, rm: RoundingMode) -> Self {
        let Ball { mid: m, rad: r } = self;
        if !Self::above_zero(&m, &r) {
            return Self::nan_ball(p);
        }
        let mid = m.sqrt(p, rm);
        let lo = m.sub(&r, p, RoundingMode::Down);
        if lo.is_negative() || lo.is_zero() {
            return Self::nan_ball(p);
        }
        let lo = lo.sqrt(p, RoundingMode::Down);
        let two = ExactNum::from_u8(2, p);
        let den = two.mul(&lo, p, RoundingMode::Down);
        let lip = ExactNum::from_u8(1, p).div(&den, p, RoundingMode::Up);
        let rad = lip.mul(&r, p, RoundingMode::Up).add(
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
        let rad = Self::transcendental_slack(&mid, p).add(&self.rad, p, RoundingMode::Up);
        Ball { mid, rad }
    }

    /// \(J_1\) of a ball. \(\lvert J_1'\rvert=\lvert J_0-J_1/x\rvert\le 1+1/(\lvert m\rvert-r)\)
    /// when the ball excludes \(0\).
    pub fn bessel_j1(self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let Ball { mid: m, rad: r } = self;
        if m.is_nan() || r.is_nan() || m.is_zero() || !matches!(m.abs_cmp(&r), Some(c) if c > 0) {
            return Self::nan_ball(p);
        }
        let mid = m.bessel_j(1, p, rm, cc);
        let den = if m.is_negative() {
            m.neg().sub(&r, p, RoundingMode::Down)
        } else {
            m.sub(&r, p, RoundingMode::Down)
        };
        if den.is_zero() {
            return Self::nan_ball(p);
        }
        let extra = ExactNum::from_u8(1, p).div(&den, p, RoundingMode::Up);
        let lip = ExactNum::from_u8(1, p).add(&extra, p, RoundingMode::Up);
        let rad = lip.mul(&r, p, RoundingMode::Up).add(
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

    /// `mid` is a positive number strictly above `rad` (the lower endpoint is then positive).
    fn above_zero(mid: &ExactNum, rad: &ExactNum) -> bool {
        if mid.is_nan() || rad.is_nan() || !mid.is_positive() {
            return false;
        }
        matches!(mid.compare(rad), Some(c) if c > 0)
    }

    /// True when `x` lies in `[mid − rad, mid + rad]` (NaN / Inf never contained).
    pub fn contains(&self, x: &ExactNum, p: usize) -> bool {
        if x.is_nan() || self.mid.is_nan() || self.rad.is_nan() {
            return false;
        }
        let d = ExactNum::sub_ref(&self.mid, x, p, RoundingMode::None).abs();
        matches!(d.compare(&self.rad), Some(c) if c <= 0)
    }
}

/// Disk enclosure in \(\mathbb{C}\): center `mid`, radius `rad`.
#[derive(Clone, Debug)]
pub struct ComplexBall {
    mid: ExactComplex,
    rad: ExactNum,
}

impl ComplexBall {
    /// Disk `mid` with radius `rad` (taken in absolute value).
    pub fn new(mid: ExactComplex, rad: ExactNum) -> Self {
        ComplexBall {
            mid,
            rad: rad.abs(),
        }
    }

    /// Center.
    pub fn mid(&self) -> &ExactComplex {
        &self.mid
    }

    /// Non-negative radius.
    ///
    /// Accuracy: enclosure
    pub fn rad(&self) -> &ExactNum {
        &self.rad
    }

    fn slack(mid: &ExactComplex, p: usize) -> ExactNum {
        let u_re = Ball::rounding_ulp(mid.re(), p);
        let u_im = Ball::rounding_ulp(mid.im(), p);
        let u = if matches!(u_re.compare(&u_im), Some(c) if c >= 0) { u_re } else { u_im };
        u.mul(
            &ExactNum::from_u8(BALL_TRANSCENDENTAL_ERROR_TERMS as u8, p),
            p,
            RoundingMode::Up,
        )
    }

    fn nan_disk() -> Self {
        let n = ExactNum::nan(Some(Error::InvalidArgument));
        ComplexBall {
            mid: ExactComplex::new(n.clone(), n.clone()),
            rad: n,
        }
    }

    /// Sum of two disks: radii add, plus rounding slack.
    pub fn add(self, other: &Self, p: usize, rm: RoundingMode) -> Self {
        let ComplexBall { mid: m1, rad: r1 } = self;
        let mid = m1.add(&other.mid, p, rm);
        let rad =
            r1.add(&other.rad, p, RoundingMode::Up)
                .add(&Self::slack(&mid, p), p, RoundingMode::Up);
        ComplexBall { mid, rad }
    }

    /// Product of two disks: \(\lvert m_1\rvert r_2+\lvert m_2\rvert r_1+r_1 r_2\) plus slack.
    pub fn mul(&self, other: &Self, p: usize, rm: RoundingMode) -> Self {
        let mid = self.mid.mul(&other.mid, p, rm);
        let a = self
            .mid
            .abs(p, RoundingMode::Up)
            .mul(&other.rad, p, RoundingMode::Up);
        let b = other
            .mid
            .abs(p, RoundingMode::Up)
            .mul(&self.rad, p, RoundingMode::Up);
        let c = ExactNum::mul_ref(&self.rad, &other.rad, p, RoundingMode::Up);
        let rad = a
            .add(&b, p, RoundingMode::Up)
            .add(&c, p, RoundingMode::Up)
            .add(&Self::slack(&mid, p), p, RoundingMode::Up);
        ComplexBall { mid, rad }
    }

    /// Exponential. Lipschitz \(\lvert e^z\rvert\le \exp(\mathrm{Re}\,m+r)\).
    pub fn exp(self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let ComplexBall { mid: m, rad: r } = self;
        let mid = m.exp(p, rm, cc);
        let (m_re, _m_im) = m.into_parts();
        let re_hi = m_re.add(&r, p, RoundingMode::Up);
        let lip = re_hi.exp(p, RoundingMode::Up, cc);
        let rad = lip
            .mul(&r, p, RoundingMode::Up)
            .add(&Self::slack(&mid, p), p, RoundingMode::Up);
        ComplexBall { mid, rad }
    }

    /// Principal logarithm. Domain: the disk must exclude \(0\).
    /// Lipschitz \(1/(\lvert m\rvert-r)\).
    pub fn ln(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let am = self.mid.abs(p, RoundingMode::Down);
        if matches!(am.compare(&self.rad), Some(c) if c <= 0) || am.is_zero() || self.rad.is_nan() {
            return Self::nan_disk();
        }
        let mid = self.mid.ln(p, rm, cc);
        let den = am.sub(&self.rad, p, RoundingMode::Down);
        let lip = ExactNum::from_u8(1, p).div(&den, p, RoundingMode::Up);
        let rad =
            lip.mul(&self.rad, p, RoundingMode::Up)
                .add(&Self::slack(&mid, p), p, RoundingMode::Up);
        ComplexBall { mid, rad }
    }

    /// Sine. \(\lvert\cos(z)\rvert\le\cosh(\lvert\mathrm{Im}\,m\rvert+r)\).
    pub fn sin(self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let ComplexBall { mid: m, rad: r } = self;
        let mid = m.sin(p, rm, cc);
        let (_m_re, im) = m.into_parts();
        // |im| + rad
        let im_hi = if im.is_negative() {
            im.neg().add(&r, p, RoundingMode::Up)
        } else {
            im.add(&r, p, RoundingMode::Up)
        };
        let lip = im_hi.sinh_cosh(p, RoundingMode::Up, cc).1;
        let rad = lip
            .mul(&r, p, RoundingMode::Up)
            .add(&Self::slack(&mid, p), p, RoundingMode::Up);
        ComplexBall { mid, rad }
    }

    /// Cosine. Same Lipschitz bound as `sin`.
    pub fn cos(self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let ComplexBall { mid: m, rad: r } = self;
        let mid = m.cos(p, rm, cc);
        let (_m_re, im) = m.into_parts();
        // |im| + rad
        let im_hi = if im.is_negative() {
            im.neg().add(&r, p, RoundingMode::Up)
        } else {
            im.add(&r, p, RoundingMode::Up)
        };
        let lip = im_hi.sinh_cosh(p, RoundingMode::Up, cc).1;
        let rad = lip
            .mul(&r, p, RoundingMode::Up)
            .add(&Self::slack(&mid, p), p, RoundingMode::Up);
        ComplexBall { mid, rad }
    }

    /// True when \(\lvert z-\mathrm{mid}\rvert\le\mathrm{rad}\).
    pub fn contains(&self, z: &ExactComplex, p: usize) -> bool {
        if z.is_nan() || self.mid.is_nan() || self.rad.is_nan() {
            return false;
        }
        let d = ExactComplex::sub_ref(z, &self.mid, p, RoundingMode::None).abs(p, RoundingMode::Up);
        matches!(d.compare(&self.rad), Some(c) if c <= 0)
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

/// Same Ziv loop as [`ziv_round`], with every input lifted to `p_wrk` before `compute`.
///
/// `MAX_PREC_RETRY` still bounds the number of bumps via [`bump_prec_retry`].
pub fn ziv_round_vec<F>(p: usize, rm: RoundingMode, inputs: &[ExactNum], mut compute: F) -> ExactNum
where
    F: FnMut(usize, &[ExactNum]) -> ExactNum,
{
    let mut p_inc = WORD_BIT_SIZE;
    let mut p_wrk = match round_p(p).checked_add(p_inc) {
        Some(v) => v,
        None => return ExactNum::nan(Some(Error::InvalidArgument)),
    };
    loop {
        let mut xs = alloc::vec::Vec::with_capacity(inputs.len());
        for x in inputs {
            let mut y = x.clone();
            if y.set_precision(p_wrk, RoundingMode::None).is_err() {
                y = x.clone();
            }
            xs.push(y);
        }
        let mut v = compute(p_wrk, &xs);
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
        let mk_u = || ExactNum::from_word(1, p);
        let ba = Ball::new(ExactNum::from(3), mk_u());
        let bb = Ball::new(ExactNum::from(4), mk_u());
        let sum = ba.add(&bb, p, rm);
        let true_sum = ExactNum::from(3).add(&ExactNum::from(4), p, rm);
        assert!(sum.contains(&true_sum, p));
    }

    #[test]
    fn ziv_round_vec_hypot_atan2() -> Result<(), crate::Error> {
        let rm = RoundingMode::ToEven;
        for p in [64usize, 128, 256] {
            let a = ExactNum::from_u8(3, p);
            let b = ExactNum::from_u8(4, p);
            let got = ziv_round_vec(p, rm, &[a, b], |pw, xs| {
                xs[0].hypot(&xs[1], pw, RoundingMode::None)
            });
            assert_eq!(got.compare(&ExactNum::from_u8(5, p)), Some(0));
        }

        let p = 256;
        let mut cc = Consts::new()?;
        let got = ziv_round_vec(
            p,
            rm,
            &[ExactNum::from_u8(1, p), ExactNum::from_u8(1, p)],
            |pw, xs| xs[0].atan2(&xs[1], pw, RoundingMode::None, &mut cc),
        );
        let quarter = cc.pi(p, rm).div(&ExactNum::from_u8(4, p), p, rm);
        assert_eq!(got.compare(&quarter), Some(0));
        Ok(())
    }

    #[test]
    fn ziv_round_sqrt_matches_direct() {
        let p = 128;
        let rm = RoundingMode::ToEven;
        let two = ExactNum::from(2);
        let via_ziv = ziv_round(p, rm, |pw| two.sqrt(pw, RoundingMode::None));
        let direct = two.sqrt(p, rm);
        assert_eq!(via_ziv.compare(&direct), Some(0));
    }

    #[test]
    fn ball_exp_contains_one_and_two() -> Result<(), crate::Error> {
        let p = 128;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let two = ExactNum::from_u8(2, p);
        let mk_rad = || two.powsi(-20, p, rm);
        let z = Ball::new(ExactNum::from_u8(0, p), mk_rad());
        let ez = z.exp(p, rm, &mut cc);
        assert!(ez.contains(&ExactNum::from_u8(1, p), p));

        let ln2 = cc.ln_2(p, rm);
        let bln = Ball::new(ln2, mk_rad());
        let e2 = bln.exp(p, rm, &mut cc);
        assert!(e2.contains(&two, p));
        Ok(())
    }

    #[test]
    fn ball_sin_contains_zero_at_origin_and_pi() -> Result<(), crate::Error> {
        let p = 128;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let two = ExactNum::from_u8(2, p);
        let mk_rad = || two.powsi(-20, p, rm);
        let z = Ball::new(ExactNum::from_u8(0, p), mk_rad());
        let sz = z.sin(p, rm, &mut cc);
        assert!(sz.contains(&ExactNum::from_u8(0, p), p));

        let pi = cc.pi(p, rm);
        let bpi = Ball::new(pi, mk_rad());
        let sp = bpi.sin(p, rm, &mut cc);
        assert!(sp.contains(&ExactNum::from_u8(0, p), p));
        Ok(())
    }

    #[test]
    fn ball_exp_sin_contain_scalar_at_a_point() -> Result<(), crate::Error> {
        let p = 128;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let x = ExactNum::from_u8(1, p);
        let rad = ExactNum::from_u8(2, p).powsi(-12, p, rm);
        let b = Ball::new(ExactNum::from_u8(1, p), rad);
        let hi = x.exp(256, rm, &mut cc);
        assert!(b.exp(p, rm, &mut cc).contains(&hi, p));
        let hs = x.sin(256, rm, &mut cc);
        assert!(b.sin(p, rm, &mut cc).contains(&hs, p));
        Ok(())
    }

    #[test]
    fn ball_plan_transcendental_golds() -> Result<(), crate::Error> {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let mk_rad = || ExactNum::from_u8(1, p).ldexp(-(p as i32), p, RoundingMode::None);
        let six = ExactNum::from_u8(6, p);
        let pi6 = cc.pi(p, rm).div(&six, p, rm);
        let bsin = Ball::new(pi6, mk_rad());
        let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, rm);
        assert!(bsin.sin(p, rm, &mut cc).contains(&half, p));

        let one = ExactNum::from_u8(1, p);
        let mk_be = || Ball::new(ExactNum::from_u8(1, p), mk_rad());
        let be = mk_be();
        let e = cc.e(p, rm);
        assert!(be.exp(p, rm, &mut cc).contains(&e, p));

        let erf1 = one.erf(p, rm, &mut cc);
        assert!(be.erf(p, rm, &mut cc).contains(&erf1, p));

        let ln1 = one.ln(256, rm, &mut cc);
        assert!(mk_be().ln(p, rm, &mut cc).contains(&ln1, p));
        let sq = one.sqrt(256, rm);
        assert!(mk_be().sqrt(p, rm).contains(&sq, p));
        assert!(be
            .cos(p, rm, &mut cc)
            .contains(&one.cos(256, rm, &mut cc), p));
        assert!(be
            .bessel_j0(p, rm, &mut cc)
            .contains(&one.bessel_j(0, 256, rm, &mut cc), p));
        assert!(mk_be()
            .bessel_j1(p, rm, &mut cc)
            .contains(&one.bessel_j(1, 256, rm, &mut cc), p));

        let composed = be.sin(p, rm, &mut cc).exp(p, rm, &mut cc);
        let true_c = one.sin(256, rm, &mut cc).exp(256, rm, &mut cc);
        assert!(composed.contains(&true_c, p));
        Ok(())
    }

    #[test]
    fn complex_ball_exp_and_pythagoras() -> Result<(), crate::Error> {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let mk_mid = || {
            ExactComplex::new(
                ExactNum::from_u8(3, p).div(&ExactNum::from_u8(10, p), p, rm),
                ExactNum::from_u8(2, p).div(&ExactNum::from_u8(10, p), p, rm),
            )
        };
        let mid = mk_mid();
        let rad = ExactNum::from_u8(1, p);
        let unit = ComplexBall::new(ExactComplex::zero(p), rad);
        let e_mid = mid.exp(p, rm, &mut cc);
        assert!(unit.exp(p, rm, &mut cc).contains(&e_mid, p));

        let mk_d = || {
            ComplexBall::new(
                mk_mid(),
                ExactNum::from_u8(1, p).ldexp(-40, p, RoundingMode::None),
            )
        };
        let d = mk_d();
        let s = mk_d().sin(p, rm, &mut cc);
        let c = d.cos(p, rm, &mut cc);
        let ss = s.mul(&s, p, rm);
        let cc2 = c.mul(&c, p, rm);
        let py = ss.add(&cc2, p, rm);
        assert!(py.contains(&ExactComplex::one(p), p));
        Ok(())
    }
}
