//! Rectangular complex numbers over [`ExactNum`].

use crate::defs::DEFAULT_P;
use crate::Consts;
use crate::ExactNum;
use crate::RoundingMode;

/// Complex value `re + i·im` with software-limb real and imaginary parts.
#[derive(Debug, Clone)]
pub struct ExactComplex {
    re: ExactNum,
    im: ExactNum,
}

impl ExactComplex {
    /// Constructs `re + i·im`.
    pub fn new(re: ExactNum, im: ExactNum) -> Self {
        Self { re, im }
    }

    /// Real part.
    pub fn re(&self) -> &ExactNum {
        &self.re
    }

    /// Imaginary part.
    pub fn im(&self) -> &ExactNum {
        &self.im
    }

    /// `0 + 0i` at precision `p`.
    pub fn zero(p: usize) -> Self {
        Self::new(ExactNum::new(p), ExactNum::new(p))
    }

    /// `1 + 0i` at precision `p`.
    pub fn one(p: usize) -> Self {
        Self::new(ExactNum::from_u8(1, p), ExactNum::new(p))
    }

    /// `0 + 1i` at precision `p`.
    pub fn i(p: usize) -> Self {
        Self::new(ExactNum::new(p), ExactNum::from_u8(1, p))
    }

    /// True if either part is NaN.
    pub fn is_nan(&self) -> bool {
        self.re.is_nan() || self.im.is_nan()
    }

    /// Complex conjugate.
    pub fn conj(&self) -> Self {
        Self::new(self.re.clone(), self.im.neg())
    }

    /// Modulus `|z| = hypot(re, im)` at precision `p`.
    pub fn abs(&self, p: usize, rm: RoundingMode) -> ExactNum {
        self.re.hypot(&self.im, p, rm)
    }

    /// Argument `atan2(im, re)` at precision `p`.
    pub fn arg(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> ExactNum {
        self.im.atan2(&self.re, p, rm, cc)
    }

    /// `self + rhs` at precision `p`.
    pub fn add(&self, rhs: &Self, p: usize, rm: RoundingMode) -> Self {
        Self::new(self.re.add(&rhs.re, p, rm), self.im.add(&rhs.im, p, rm))
    }

    /// `self - rhs` at precision `p`.
    pub fn sub(&self, rhs: &Self, p: usize, rm: RoundingMode) -> Self {
        Self::new(self.re.sub(&rhs.re, p, rm), self.im.sub(&rhs.im, p, rm))
    }

    /// `self * rhs` at precision `p`.
    pub fn mul(&self, rhs: &Self, p: usize, rm: RoundingMode) -> Self {
        let ac = self.re.mul(&rhs.re, p, RoundingMode::None);
        let bd = self.im.mul(&rhs.im, p, RoundingMode::None);
        let ad = self.re.mul(&rhs.im, p, RoundingMode::None);
        let bc = self.im.mul(&rhs.re, p, RoundingMode::None);
        Self::new(ac.sub(&bd, p, rm), ad.add(&bc, p, rm))
    }

    /// `self / rhs` at precision `p`.
    pub fn div(&self, rhs: &Self, p: usize, rm: RoundingMode) -> Self {
        let ac = self.re.mul(&rhs.re, p, RoundingMode::None);
        let bd = self.im.mul(&rhs.im, p, RoundingMode::None);
        let bc = self.im.mul(&rhs.re, p, RoundingMode::None);
        let ad = self.re.mul(&rhs.im, p, RoundingMode::None);
        let den = rhs.re.mul(&rhs.re, p, RoundingMode::None).add(
            &rhs.im.mul(&rhs.im, p, RoundingMode::None),
            p,
            RoundingMode::None,
        );
        Self::new(
            ac.add(&bd, p, rm).div(&den, p, rm),
            bc.sub(&ad, p, rm).div(&den, p, rm),
        )
    }

    /// `e^self` using `exp(re) (cos(im) + i sin(im))`.
    pub fn exp(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let er = self.re.exp(p, RoundingMode::None, cc);
        let (s, c) = self.im.sin_cos(p, RoundingMode::None, cc);
        Self::new(er.mul(&c, p, rm), er.mul(&s, p, rm))
    }

    /// Principal logarithm `ln|z| + i Arg(z)`.
    pub fn ln(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let mag = self.abs(p, RoundingMode::None);
        Self::new(mag.ln(p, rm, cc), self.arg(p, rm, cc))
    }

    /// `sin(self)` via `sin(re)cosh(im) + i cos(re)sinh(im)`.
    pub fn sin(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let (sn, cs) = self.re.sin_cos(p, RoundingMode::None, cc);
        let (sh, ch) = self.im.sinh_cosh(p, RoundingMode::None, cc);
        Self::new(sn.mul(&ch, p, rm), cs.mul(&sh, p, rm))
    }

    /// `cos(self)` via `cos(re)cosh(im) - i sin(re)sinh(im)`.
    pub fn cos(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let (sn, cs) = self.re.sin_cos(p, RoundingMode::None, cc);
        let (sh, ch) = self.im.sinh_cosh(p, RoundingMode::None, cc);
        Self::new(cs.mul(&ch, p, rm), sn.mul(&sh, p, rm).neg())
    }
}

macro_rules! impl_cplx_binop {
    ($trait:ident, $method:ident) => {
        impl core::ops::$trait<&ExactComplex> for &ExactComplex {
            type Output = ExactComplex;
            fn $method(self, rhs: &ExactComplex) -> ExactComplex {
                ExactComplex::$method(self, rhs, DEFAULT_P, RoundingMode::ToEven)
            }
        }
        impl core::ops::$trait<ExactComplex> for &ExactComplex {
            type Output = ExactComplex;
            fn $method(self, rhs: ExactComplex) -> ExactComplex {
                ExactComplex::$method(self, &rhs, DEFAULT_P, RoundingMode::ToEven)
            }
        }
        impl core::ops::$trait<&ExactComplex> for ExactComplex {
            type Output = ExactComplex;
            fn $method(self, rhs: &ExactComplex) -> ExactComplex {
                ExactComplex::$method(&self, rhs, DEFAULT_P, RoundingMode::ToEven)
            }
        }
        impl core::ops::$trait<ExactComplex> for ExactComplex {
            type Output = ExactComplex;
            fn $method(self, rhs: ExactComplex) -> ExactComplex {
                ExactComplex::$method(&self, &rhs, DEFAULT_P, RoundingMode::ToEven)
            }
        }
    };
}

impl_cplx_binop!(Add, add);
impl_cplx_binop!(Sub, sub);
impl_cplx_binop!(Mul, mul);
impl_cplx_binop!(Div, div);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::NAN;

    #[test]
    fn test_complex_arith() {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new().unwrap();

        let a = ExactComplex::new(ExactNum::from_u8(3, p), ExactNum::from_u8(4, p));
        let mag = a.abs(p, rm);
        let five = ExactNum::from_u8(5, p);
        assert_eq!(mag.cmp(&five), Some(0));

        let i = ExactComplex::i(p);
        let i2 = i.mul(&i, p, rm);
        assert_eq!(i2.re().cmp(&ExactNum::from_i8(-1, p)), Some(0));
        assert!(i2.im().is_zero());

        let z = ExactComplex::one(p);
        let e = z.exp(p, rm, &mut cc);
        let back = e.ln(p, rm, &mut cc);
        let d = back.re().sub(z.re(), p, RoundingMode::None).abs();
        assert!(d.exponent().unwrap_or(0) < -((p as i32) / 4));
    }

    #[test]
    fn test_complex_nan() {
        let n = ExactComplex::new(NAN.clone(), ExactNum::new(64));
        assert!(n.is_nan());
    }
}
