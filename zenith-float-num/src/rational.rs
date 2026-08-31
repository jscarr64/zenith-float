//! Exact rationals `num/den` with `ExactNum` limbs, reduced to lowest terms.

use crate::defs::RoundingMode;
use crate::defs::DEFAULT_P;
use crate::defs::WORD_BIT_SIZE;
use crate::ExactNum;
use crate::NAN;
use core::cmp::Ordering;

/// Exact rational `num/den`. Distinct from [`ExactNum`], which is floating-point.
///
/// Both parts are integer-valued `ExactNum`s after [`Self::new`]. The
/// denominator is positive. A zero denominator or a non-finite part is stored
/// as `NaN/1`.
#[derive(Clone, Debug)]
pub struct ExactRational {
    num: ExactNum,
    den: ExactNum,
}

fn rat_one() -> ExactNum {
    ExactNum::from_u8(1, DEFAULT_P)
}

fn work_p(xs: &[&ExactNum]) -> usize {
    let mut p = DEFAULT_P;
    for x in xs {
        if let Some(px) = x.precision() {
            p = p.max(px);
        }
    }
    p.saturating_add(WORD_BIT_SIZE)
}

fn is_int(x: &ExactNum) -> bool {
    if x.is_nan() || x.is_inf() {
        return false;
    }
    x.fract().is_zero()
}

fn gcd_int(mut a: ExactNum, mut b: ExactNum) -> ExactNum {
    a = a.abs();
    b = b.abs();
    while !b.is_zero() {
        let r = a.rem(&b);
        a = b;
        b = r;
    }
    if a.is_zero() {
        rat_one()
    } else {
        a
    }
}

fn reduce(mut num: ExactNum, mut den: ExactNum) -> ExactRational {
    if num.is_nan() || den.is_nan() || num.is_inf() || den.is_inf() || den.is_zero() {
        return ExactRational {
            num: NAN,
            den: rat_one(),
        };
    }
    if den.is_negative() {
        num = num.neg();
        den = den.neg();
    }
    if is_int(&num) && is_int(&den) {
        let g = gcd_int(num.clone(), den.clone());
        if !g.is_zero() && g.cmp(&rat_one()) != Some(0) {
            let p = work_p(&[&num, &den, &g]);
            num = num.div(&g, p, RoundingMode::ToEven);
            den = den.div(&g, p, RoundingMode::ToEven);
        }
    }
    ExactRational { num, den }
}

impl ExactRational {
    /// `num/den` reduced to lowest terms. Negative `den` moves the sign to `num`.
    /// Zero `den` or a non-finite part is `NaN`.
    pub fn new(num: ExactNum, den: ExactNum) -> Self {
        reduce(num, den)
    }

    /// Integer ratio `n/d` at default limb precision.
    pub fn from_i64(n: i64, d: i64) -> Self {
        Self::new(
            ExactNum::from_i64(n, DEFAULT_P),
            ExactNum::from_i64(d, DEFAULT_P),
        )
    }

    /// Numerator (integer-valued after construction).
    pub fn num(&self) -> &ExactNum {
        &self.num
    }

    /// Denominator (positive integer-valued after construction).
    pub fn den(&self) -> &ExactNum {
        &self.den
    }

    /// True if a part is NaN or the denominator was zero.
    pub fn is_nan(&self) -> bool {
        self.num.is_nan() || self.den.is_nan()
    }

    /// True if the denominator is `1`.
    pub fn is_integer(&self) -> bool {
        !self.is_nan() && self.den.cmp(&rat_one()) == Some(0)
    }

    /// `self + rhs` as an exact rational.
    pub fn add(&self, rhs: &Self) -> Self {
        if self.is_nan() || rhs.is_nan() {
            return Self::new(NAN, rat_one());
        }
        let p = work_p(&[&self.num, &self.den, &rhs.num, &rhs.den]);
        let rm = RoundingMode::None;
        let ad = self.num.mul(&rhs.den, p, rm);
        let bc = rhs.num.mul(&self.den, p, rm);
        let num = ad.add(&bc, p, rm);
        let den = self.den.mul(&rhs.den, p, rm);
        Self::new(num, den)
    }

    /// `self - rhs`.
    pub fn sub(&self, rhs: &Self) -> Self {
        self.add(&Self::new(rhs.num.neg(), rhs.den.clone()))
    }

    /// `self * rhs`.
    pub fn mul(&self, rhs: &Self) -> Self {
        if self.is_nan() || rhs.is_nan() {
            return Self::new(NAN, rat_one());
        }
        let p = work_p(&[&self.num, &self.den, &rhs.num, &rhs.den]);
        let rm = RoundingMode::None;
        Self::new(self.num.mul(&rhs.num, p, rm), self.den.mul(&rhs.den, p, rm))
    }

    /// `self / rhs`. Zero `rhs` is `NaN`.
    pub fn div(&self, rhs: &Self) -> Self {
        if self.is_nan() || rhs.is_nan() || rhs.num.is_zero() {
            return Self::new(NAN, rat_one());
        }
        let p = work_p(&[&self.num, &self.den, &rhs.num, &rhs.den]);
        let rm = RoundingMode::None;
        Self::new(self.num.mul(&rhs.den, p, rm), self.den.mul(&rhs.num, p, rm))
    }

    /// Convert to an `ExactNum` at `(p, rm)`.
    pub fn to_exact_num(&self, p: usize, rm: RoundingMode) -> ExactNum {
        if self.is_nan() {
            return NAN;
        }
        self.num.div(&self.den, p, rm)
    }

    fn quot(&self) -> ExactNum {
        let p = work_p(&[&self.num, &self.den]);
        self.num.div(&self.den, p, RoundingMode::None)
    }

    /// Greatest integer `≤ self` as a rational with denominator `1`.
    pub fn floor(&self) -> Self {
        if self.is_nan() {
            return Self::new(NAN, rat_one());
        }
        Self::new(self.quot().floor(), rat_one())
    }

    /// Least integer `≥ self` as a rational with denominator `1`.
    pub fn ceil(&self) -> Self {
        if self.is_nan() {
            return Self::new(NAN, rat_one());
        }
        Self::new(self.quot().ceil(), rat_one())
    }

    /// Nearest integer, ties away from zero, as a rational with denominator `1`.
    pub fn round(&self) -> Self {
        if self.is_nan() {
            return Self::new(NAN, rat_one());
        }
        Self::new(self.quot().round(0, RoundingMode::FromZero), rat_one())
    }

    /// Exact comparison via cross-multiplication. `None` if either value is NaN.
    pub fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.is_nan() || other.is_nan() {
            return None;
        }
        let p = work_p(&[&self.num, &self.den, &other.num, &other.den]);
        let rm = RoundingMode::None;
        let left = self.num.mul(&other.den, p, rm);
        let right = other.num.mul(&self.den, p, rm);
        match left.cmp(&right) {
            Some(c) if c > 0 => Some(Ordering::Greater),
            Some(c) if c < 0 => Some(Ordering::Less),
            Some(_) => Some(Ordering::Equal),
            None => None,
        }
    }
}

impl PartialEq for ExactRational {
    fn eq(&self, other: &Self) -> bool {
        matches!(self.partial_cmp(other), Some(Ordering::Equal))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Consts;
    use crate::Radix;

    fn r(n: i64, d: i64) -> ExactRational {
        ExactRational::new(
            ExactNum::from_i64(n, DEFAULT_P),
            ExactNum::from_i64(d, DEFAULT_P),
        )
    }

    #[test]
    fn rational_reduce_add_sign_to_float() {
        let half = r(1, 2);
        let two_four = r(2, 4);
        assert_eq!(two_four, half);
        assert_eq!(ExactRational::from_i64(2, 4), ExactRational::from_i64(1, 2));
        assert!(!two_four.is_integer());
        assert_eq!(
            two_four.num().cmp(&ExactNum::from_u8(1, DEFAULT_P)),
            Some(0)
        );
        assert_eq!(
            two_four.den().cmp(&ExactNum::from_u8(2, DEFAULT_P)),
            Some(0)
        );

        let s = r(1, 3).add(&r(1, 6));
        assert_eq!(s, half);

        let pos = ExactRational::from_i64(-3, -4);
        let want = ExactRational::from_i64(3, 4);
        assert_eq!(pos, want);
        assert!(pos.den().is_positive());

        let p = 256;
        let rm = RoundingMode::ToEven;
        let got = ExactRational::from_i64(1, 3).to_exact_num(p, rm);
        let want_f = ExactNum::from_i64(1, p).div(&ExactNum::from_i64(3, p), p, rm);
        assert_eq!(got.cmp(&want_f), Some(0));
        let mut cc = Consts::new().unwrap();
        let parsed = ExactNum::parse(
            "0.33333333333333333333333333333333333333333333333333333333333333333333333333333333",
            Radix::Dec,
            p,
            rm,
            &mut cc,
        );
        assert_eq!(got.cmp(&parsed), Some(0));

        assert!(ExactRational::from_i64(1, 3).partial_cmp(&half) == Some(Ordering::Less));
        assert!(ExactRational::from_i64(2, 2).is_integer());
        assert_eq!(
            ExactRational::from_i64(5, 3).floor(),
            ExactRational::from_i64(1, 1)
        );
        assert_eq!(
            ExactRational::from_i64(5, 3).ceil(),
            ExactRational::from_i64(2, 1)
        );
        assert!(ExactRational::from_i64(1, 0).is_nan());
    }
}
