//! Rectangular complex numbers over [`ExactNum`].

use crate::defs::DEFAULT_P;
use crate::Consts;
use crate::ExactNum;
use crate::Exponent;
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
    ///
    /// Branch: same as real `atan2`; values lie in (−π, π]. The cut of `ln` / `sqrt` /
    /// `pow` is the non-positive real axis, approached from above as +π and from below as −π.
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
    ///
    /// Branch cut: (−∞, 0] on the real axis. `ln(−1)` is `iπ` (argument +π).
    pub fn ln(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let mag = self.abs(p, RoundingMode::None);
        Self::new(mag.ln(p, rm, cc), self.arg(p, rm, cc))
    }

    /// `sin(self)` via `sin(re)cosh(im) + i cos(re)sinh(im)`.
    ///
    /// The complex value is **not** passed to `rem_pi`. Only the real (resp. imaginary)
    /// *component* uses real `sin_cos` / `sinh_cosh`.
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
    ///
    /// Branch cut: (−∞, 0]. Real part of the result is ≥ 0. `sqrt(−1)` is `+i`.
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

    /// `log2(self) = ln(self) / ln 2` (principal branch).
    pub fn log2(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        let ln = self.ln(p_x, RoundingMode::None, cc);
        let base = Self::from_real(cc.ln_2(p_x, RoundingMode::None), p_x);
        ln.div(&base, p_x, RoundingMode::None).finish(p, rm)
    }

    /// `log10(self) = ln(self) / ln 10` (principal branch).
    pub fn log10(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        let ln = self.ln(p_x, RoundingMode::None, cc);
        let base = Self::from_real(cc.ln_10(p_x, RoundingMode::None), p_x);
        ln.div(&base, p_x, RoundingMode::None).finish(p, rm)
    }

    /// `log_base(self) = ln(self) / ln(base)` (principal branch).
    pub fn log(&self, base: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        let ln = self.ln(p_x, RoundingMode::None, cc);
        let lnb = base.ln(p_x, RoundingMode::None, cc);
        ln.div(&lnb, p_x, RoundingMode::None).finish(p, rm)
    }

    /// `ln(1 + self)` (principal branch).
    pub fn log1p(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        Self::one(p_x)
            .add(self, p_x, RoundingMode::None)
            .ln(p_x, RoundingMode::None, cc)
            .finish(p, rm)
    }

    /// `exp2(self) = exp(self · ln 2)`.
    pub fn exp2(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        let ln2 = Self::from_real(cc.ln_2(p_x, RoundingMode::None), p_x);
        self.mul(&ln2, p_x, RoundingMode::None)
            .exp(p_x, RoundingMode::None, cc)
            .finish(p, rm)
    }

    /// `exp10(self) = exp(self · ln 10)`.
    pub fn exp10(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        let ln10 = Self::from_real(cc.ln_10(p_x, RoundingMode::None), p_x);
        self.mul(&ln10, p_x, RoundingMode::None)
            .exp(p_x, RoundingMode::None, cc)
            .finish(p, rm)
    }

    /// `exp(self) − 1`.
    pub fn expm1(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        self.exp(p_x, RoundingMode::None, cc)
            .sub(&Self::one(p_x), p_x, RoundingMode::None)
            .finish(p, rm)
    }

    /// Scale both parts by `2^n` (`ldexp` on re and im).
    pub fn ldexp(&self, n: Exponent, p: usize, rm: RoundingMode) -> Self {
        self.ldexp_parts(n, p, rm)
    }

    /// Same as [`ldexp`](Self::ldexp).
    pub fn scalb(&self, n: Exponent, p: usize, rm: RoundingMode) -> Self {
        self.ldexp(n, p, rm)
    }

    /// `logb(|z|)` as a real (`x + 0i`).
    pub fn logb(&self, p: usize, rm: RoundingMode) -> ExactNum {
        self.abs(p, rm).logb(p, rm)
    }

    /// Principal `n`-th root via `exp(ln(z) / n)`. Inherits the `ln` branch cut.
    pub fn nth_root(&self, n: usize, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if n == 0 {
            return Self::new(
                ExactNum::nan(Some(crate::Error::InvalidArgument)),
                ExactNum::nan(Some(crate::Error::InvalidArgument)),
            );
        }
        if n == 1 {
            let mut z = self.clone();
            if let Err(e) = z.set_precision(p, rm) {
                return Self::new(ExactNum::nan(Some(e)), ExactNum::nan(Some(e)));
            }
            return z;
        }
        if n == 2 {
            return self.sqrt(p, rm, cc);
        }
        let p_x = Self::work_p(p);
        let ln = self.ln(p_x, RoundingMode::None, cc);
        let ninv = Self::from_real(
            ExactNum::from_u32(1, p_x).div(
                &ExactNum::from_u32(n as u32, p_x),
                p_x,
                RoundingMode::None,
            ),
            p_x,
        );
        ln.mul(&ninv, p_x, RoundingMode::None)
            .exp(p_x, RoundingMode::None, cc)
            .finish(p, rm)
    }

    /// Principal cube root. Same branch as [`nth_root`](Self::nth_root) with `n = 3`.
    pub fn cbrt(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.nth_root(3, p, rm, cc)
    }

    /// Principal `sqrt(self² + other²)` (analytic continuation of real `hypot`).
    pub fn hypot(&self, other: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        self.mul(self, p_x, RoundingMode::None)
            .add(
                &other.mul(other, p_x, RoundingMode::None),
                p_x,
                RoundingMode::None,
            )
            .sqrt(p_x, RoundingMode::None, cc)
            .finish(p, rm)
    }

    /// `self * b + c` at extra working precision, then one round (not a fused complex hardware op).
    pub fn fma(&self, b: &Self, c: &Self, p: usize, rm: RoundingMode) -> Self {
        let p_x = Self::work_p(p);
        self.mul(b, p_x, RoundingMode::None)
            .add(c, p_x, RoundingMode::None)
            .finish(p, rm)
    }

    /// Alias of [`fma`](Self::fma).
    pub fn mul_add(&self, b: &Self, c: &Self, p: usize, rm: RoundingMode) -> Self {
        self.fma(b, c, p, rm)
    }

    /// `self^rhs` as `exp(rhs * ln(self))` (principal branch).
    ///
    /// Inherits the `ln` cut on `self`: non-positive real base uses Arg = ±π.
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
        assert!(d.is_zero() || d.exponent().unwrap_or(0) < -((p as i32) / 4));

        let z0 = ExactComplex::zero(p);
        let t = z0.tan(p, rm, &mut cc);
        assert!(t.re().is_zero() || t.re().exponent().unwrap_or(0) < -((p as i32) / 4));
        assert!(t.im().is_zero() || t.im().exponent().unwrap_or(0) < -((p as i32) / 4));
        let sh = z0.sinh(p, rm, &mut cc);
        assert!(sh.re().is_zero() || sh.re().exponent().unwrap_or(0) < -((p as i32) / 4));

        let two = ExactComplex::from_real(ExactNum::from_u8(2, p), p);
        let four = ExactComplex::from_real(ExactNum::from_u8(4, p), p);
        let lg = four.log2(p, rm, &mut cc);
        let d = lg.re().sub(two.re(), p, RoundingMode::None).abs();
        assert!(d.is_zero() || d.exponent().unwrap_or(0) < -((p as i32) / 8));
        assert!(lg.im().is_zero() || lg.im().exponent().unwrap_or(0) < -((p as i32) / 4));
        let e2 = two.exp2(p, rm, &mut cc);
        let d2 = e2.re().sub(four.re(), p, RoundingMode::None).abs();
        assert!(d2.is_zero() || d2.exponent().unwrap_or(0) < -((p as i32) / 8));
    }

    #[test]
    fn test_complex_branch_cuts() {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new().unwrap();
        let m1 = ExactComplex::from_real(ExactNum::from_i8(-1, p), p);

        let s = m1.sqrt(p, rm, &mut cc);
        assert!(s.re().is_zero() || s.re().exponent().unwrap_or(0) < -((p as i32) / 4));
        assert_eq!(s.im().cmp(&ExactNum::from_u8(1, p)), Some(0));

        let l = m1.ln(p, rm, &mut cc);
        assert!(l.re().is_zero() || l.re().exponent().unwrap_or(0) < -((p as i32) / 4));
        let pi = cc.pi(p, rm);
        let d = l.im().abs().sub(&pi, p, RoundingMode::None).abs();
        assert!(d.is_zero() || d.exponent().unwrap_or(0) < -((p as i32) / 8));
        assert!(l.im().is_positive());
    }

    #[test]
    fn test_complex_nan() {
        let n = ExactComplex::new(NAN.clone(), ExactNum::new(64));
        assert!(n.is_nan());
    }
}
