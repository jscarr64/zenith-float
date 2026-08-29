//! Rectangular complex numbers over [`ExactNum`].

use crate::defs::DEFAULT_P;
use crate::Consts;
use crate::ExactNum;
use crate::RoundingMode;
use crate::WORD_BIT_SIZE;

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

    /// True if either part is inexact.
    pub fn inexact(&self) -> bool {
        self.re.inexact() || self.im.inexact()
    }

    /// Sets the inexact flag on both parts.
    pub fn set_inexact(&mut self, inexact: bool) {
        self.re.set_inexact(inexact);
        self.im.set_inexact(inexact);
    }

    /// Rounds both parts to precision `p`.
    pub fn set_precision(&mut self, p: usize, rm: RoundingMode) -> Result<(), crate::Error> {
        self.re.set_precision(p, rm)?;
        self.im.set_precision(p, rm)
    }

    /// Real `x` as `x + 0i`. Imaginary zero uses precision `p`.
    pub fn from_real(re: ExactNum, p: usize) -> Self {
        Self::new(re, ExactNum::new(p))
    }

    /// `1 / self`.
    pub fn reciprocal(&self, p: usize, rm: RoundingMode) -> Self {
        Self::one(p).div(self, p, rm)
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

    fn finish(self, p: usize, rm: RoundingMode) -> Self {
        let mut re = self.re;
        let mut im = self.im;
        if let Err(e) = re.set_precision(p, rm) {
            return Self::new(ExactNum::nan(Some(e)), ExactNum::nan(Some(e)));
        }
        if let Err(e) = im.set_precision(p, rm) {
            return Self::new(ExactNum::nan(Some(e)), ExactNum::nan(Some(e)));
        }
        Self::new(re, im)
    }

    fn work_p(p: usize) -> usize {
        p.saturating_add(WORD_BIT_SIZE)
    }

    /// `tan(self) = sin(self) / cos(self)`.
    pub fn tan(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        self.sin(p_x, RoundingMode::None, cc)
            .div(
                &self.cos(p_x, RoundingMode::None, cc),
                p_x,
                RoundingMode::None,
            )
            .finish(p, rm)
    }

    /// `sinh(self)` via `sinh(re)cos(im) + i cosh(re)sin(im)`.
    pub fn sinh(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        let (sh, ch) = self.re.sinh_cosh(p_x, RoundingMode::None, cc);
        let (sn, cs) = self.im.sin_cos(p_x, RoundingMode::None, cc);
        Self::new(
            sh.mul(&cs, p_x, RoundingMode::None),
            ch.mul(&sn, p_x, RoundingMode::None),
        )
        .finish(p, rm)
    }

    /// `cosh(self)` via `cosh(re)cos(im) + i sinh(re)sin(im)`.
    pub fn cosh(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        let (sh, ch) = self.re.sinh_cosh(p_x, RoundingMode::None, cc);
        let (sn, cs) = self.im.sin_cos(p_x, RoundingMode::None, cc);
        Self::new(
            ch.mul(&cs, p_x, RoundingMode::None),
            sh.mul(&sn, p_x, RoundingMode::None),
        )
        .finish(p, rm)
    }

    /// `tanh(self) = sinh(self) / cosh(self)`.
    pub fn tanh(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        self.sinh(p_x, RoundingMode::None, cc)
            .div(
                &self.cosh(p_x, RoundingMode::None, cc),
                p_x,
                RoundingMode::None,
            )
            .finish(p, rm)
    }

    /// Principal square root: `√r (cos(θ/2) + i sin(θ/2))`.
    pub fn sqrt(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        let r = self.abs(p_x, RoundingMode::None);
        let th = self.arg(p_x, RoundingMode::None, cc);
        let sr = r.sqrt(p_x, RoundingMode::None);
        let half_th = th.ldexp(-1, p_x, RoundingMode::None);
        let (s, c) = half_th.sin_cos(p_x, RoundingMode::None, cc);
        Self::new(
            sr.mul(&c, p_x, RoundingMode::None),
            sr.mul(&s, p_x, RoundingMode::None),
        )
        .finish(p, rm)
    }

    /// `self^rhs` as `exp(rhs * ln(self))` (principal branch).
    pub fn pow(&self, rhs: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        let ln = self.ln(p_x, RoundingMode::None, cc);
        rhs.mul(&ln, p_x, RoundingMode::None)
            .exp(p_x, RoundingMode::None, cc)
            .finish(p, rm)
    }

    /// Principal `asin`: `-i ln(i z + √(1 − z²))`.
    pub fn asin(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        let z2 = self.mul(self, p_x, RoundingMode::None);
        let one = Self::one(p_x);
        let rad = one
            .sub(&z2, p_x, RoundingMode::None)
            .sqrt(p_x, RoundingMode::None, cc);
        let iz = Self::i(p_x).mul(self, p_x, RoundingMode::None);
        let ln = iz
            .add(&rad, p_x, RoundingMode::None)
            .ln(p_x, RoundingMode::None, cc);
        Self::i(p_x)
            .mul(&ln, p_x, RoundingMode::None)
            .neg()
            .finish(p, rm)
    }

    /// Principal `acos`: `π/2 − asin(z)`.
    pub fn acos(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        let half_pi = cc
            .pi(p_x, RoundingMode::None)
            .ldexp(-1, p_x, RoundingMode::None);
        let asinv = self.asin(p_x, RoundingMode::None, cc);
        Self::new(half_pi, ExactNum::new(p_x))
            .sub(&asinv, p_x, RoundingMode::None)
            .finish(p, rm)
    }

    /// Principal `atan`: `(i/2) ln((i+z)/(i−z))`.
    pub fn atan(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        let i = Self::i(p_x);
        let num = i.add(self, p_x, RoundingMode::None);
        let den = i.sub(self, p_x, RoundingMode::None);
        let ln = num
            .div(&den, p_x, RoundingMode::None)
            .ln(p_x, RoundingMode::None, cc);
        let half_i = i.ldexp_parts(-1, p_x, RoundingMode::None);
        half_i.mul(&ln, p_x, RoundingMode::None).finish(p, rm)
    }

    /// Principal `asinh`: `ln(z + √(z² + 1))`.
    pub fn asinh(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        let z2 = self.mul(self, p_x, RoundingMode::None);
        let rad =
            z2.add(&Self::one(p_x), p_x, RoundingMode::None)
                .sqrt(p_x, RoundingMode::None, cc);
        self.add(&rad, p_x, RoundingMode::None)
            .ln(p_x, RoundingMode::None, cc)
            .finish(p, rm)
    }

    /// Principal `acosh`: `ln(z + √(z−1)√(z+1))`.
    pub fn acosh(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        let one = Self::one(p_x);
        let zm = self
            .sub(&one, p_x, RoundingMode::None)
            .sqrt(p_x, RoundingMode::None, cc);
        let zp = self
            .add(&one, p_x, RoundingMode::None)
            .sqrt(p_x, RoundingMode::None, cc);
        self.add(
            &zm.mul(&zp, p_x, RoundingMode::None),
            p_x,
            RoundingMode::None,
        )
        .ln(p_x, RoundingMode::None, cc)
        .finish(p, rm)
    }

    /// Principal `atanh`: `(1/2) ln((1+z)/(1−z))`.
    pub fn atanh(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        let one = Self::one(p_x);
        let num = one.add(self, p_x, RoundingMode::None);
        let den = one.sub(self, p_x, RoundingMode::None);
        num.div(&den, p_x, RoundingMode::None)
            .ln(p_x, RoundingMode::None, cc)
            .ldexp_parts(-1, p_x, RoundingMode::None)
            .finish(p, rm)
    }

    fn ldexp_parts(&self, n: crate::Exponent, p: usize, rm: RoundingMode) -> Self {
        Self::new(self.re.ldexp(n, p, rm), self.im.ldexp(n, p, rm))
    }

    fn neg(&self) -> Self {
        Self::new(self.re.neg(), self.im.neg())
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

impl crate::FromExt<ExactComplex> for ExactComplex {
    fn from_ext(mut v: ExactComplex, p: usize, rm: RoundingMode, _cc: &mut Consts) -> Self {
        if let Err(e) = v.set_precision(p, rm) {
            return Self::new(ExactNum::nan(Some(e)), ExactNum::nan(Some(e)));
        }
        v.set_inexact(false);
        v
    }
}

impl crate::FromExt<&ExactComplex> for ExactComplex {
    fn from_ext(v: &ExactComplex, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        <Self as crate::FromExt<ExactComplex>>::from_ext(v.clone(), p, rm, cc)
    }
}

impl<T> crate::FromExt<T> for ExactComplex
where
    ExactNum: crate::FromExt<T>,
{
    fn from_ext(v: T, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        Self::from_real(<ExactNum as crate::FromExt<T>>::from_ext(v, p, rm, cc), p)
    }
}

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

        let z0 = ExactComplex::zero(p);
        let t = z0.tan(p, rm, &mut cc);
        assert!(t.re().is_zero() || t.re().exponent().unwrap_or(0) < -((p as i32) / 4));
        assert!(t.im().is_zero() || t.im().exponent().unwrap_or(0) < -((p as i32) / 4));
        let sh = z0.sinh(p, rm, &mut cc);
        assert!(sh.re().is_zero() || sh.re().exponent().unwrap_or(0) < -((p as i32) / 4));
    }

    #[test]
    fn test_complex_nan() {
        let n = ExactComplex::new(NAN.clone(), ExactNum::new(64));
        assert!(n.is_nan());
    }
}
