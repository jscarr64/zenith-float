//! Rectangular complex numbers over [`ExactNum`].

use crate::common::util::round_p;
use crate::defs::DEFAULT_P;
use crate::Consts;
use crate::ExactNum;
use crate::Exponent;
use crate::RoundingMode;
use crate::WORD_BIT_SIZE;
use core::ops::Neg;

fn abs_part(x: ExactNum) -> ExactNum {
    if x.is_negative() {
        x.neg()
    } else {
        x
    }
}

/// `a + b` (or `a - b`) of two exact products, exact whenever its width stays bounded.
///
/// When the leading exponents are further apart than both widths plus `p` and two guard words, the
/// smaller term only acts as a sticky bit below the final rounding position: the sum is then formed at
/// that bounded width without rounding (`RoundingMode::None`), so the full-width sum never has to span
/// the whole exponent range (`2^EXPONENT_MAX + 2^EXPONENT_MIN`).
fn exact_sum(a: &ExactNum, b: &ExactNum, sub: bool, p: usize) -> ExactNum {
    let wide = match (
        a.exponent(),
        b.exponent(),
        a.mantissa_max_bit_len(),
        b.mantissa_max_bit_len(),
    ) {
        (Some(ea), Some(eb), Some(la), Some(lb)) if !a.is_zero() && !b.is_zero() => {
            let gap = (i128::from(ea) - i128::from(eb)).unsigned_abs();
            let bound = (la as u128) + (lb as u128) + (p as u128) + 2 * (WORD_BIT_SIZE as u128);
            gap > bound
        }
        _ => false,
    };
    match (wide, sub) {
        (false, false) => a.add_full_prec(b),
        (false, true) => a.sub_full_prec(b),
        // The far smaller term only decides which side of a grid point the sum lies on. It is replaced
        // by a proxy of the same sign, a power of two just below the `w`-bit grid of the larger term, and
        // the sum is then formed exactly: the one later rounding at `p` sees the same side (and never a
        // tie), which a sum rounded at `w` without a sticky bit lost.
        (true, sub) => {
            let w = a.mantissa_max_bit_len().map_or(p, |la| la.max(b.mantissa_max_bit_len().map_or(p, |lb| lb)))
                .saturating_add(p)
                .saturating_add(2 * WORD_BIT_SIZE);
            let a_big = a.exponent() > b.exponent();
            let (big, small) = if a_big { (a, b) } else { (b, a) };
            let proxy = big
                .exponent()
                .and_then(|e| i64::try_from(w).ok().and_then(|w| w.checked_add(2)).and_then(|w2| e.checked_sub(w2)))
                .and_then(|e| ExactNum::pow2(e, WORD_BIT_SIZE).ok())
                .map(|t| if small.is_negative() { -t } else { t });
            match (proxy, a_big, sub) {
                (Some(t), true, false) => a.add_full_prec(&t),
                (Some(t), true, true) => a.sub_full_prec(&t),
                (Some(t), false, false) => t.add_full_prec(b),
                (Some(t), false, true) => t.sub_full_prec(b),
                (None, _, true) => ExactNum::sub_ref(a, b, w, RoundingMode::None),
                (None, _, false) => ExactNum::add_ref(a, b, w, RoundingMode::None),
            }
        }
    }
}

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

    /// Consumes the value and returns `(re, im)` by move; no copy of either part.
    pub fn into_parts(self) -> (ExactNum, ExactNum) {
        (self.re, self.im)
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

    /// True if this operation rounded either part. Operations on `ExactComplex` report the
    /// operation alone, like [`ExactNum::inexact`]; see [`Self::inexact_accumulated`].
    pub fn inexact(&self) -> bool {
        self.re.inexact() || self.im.inexact()
    }

    /// True if either part, or any value it was computed from, was rounded.
    pub fn inexact_accumulated(&self) -> bool {
        self.re.inexact_accumulated() || self.im.inexact_accumulated()
    }

    /// Copy of `self` without inexact history: the operand of a per-operation flag.
    pub(crate) fn exact_copy(&self) -> Self {
        Self::new(self.re.exact_copy_ext(), self.im.exact_copy_ext())
    }

    /// Flags `r`, computed from clean operands: the operation rounded if any step rounded
    /// either part; the accumulated flag adds the operands' history.
    pub(crate) fn per_op(inputs: bool, mut r: Self) -> Self {
        let op = r.re.inexact_accumulated() || r.im.inexact_accumulated();
        r.re.finish_part(op, inputs);
        r.im.finish_part(op, inputs);
        r
    }

    /// Sets the inexact flag on both parts.
    pub fn set_inexact(&mut self, inexact: bool) {
        self.re.set_inexact(inexact);
        self.im.set_inexact(inexact);
    }

    /// Rounds both parts to precision `p`.
    ///
    /// Accuracy: correctly-rounded
    pub fn set_precision(&mut self, p: usize, rm: RoundingMode) -> Result<(), crate::Error> {
        self.re.set_precision(p, rm)?;
        self.im.set_precision(p, rm)
    }

    /// Real `x` as `x + 0i`. Imaginary zero uses precision `p`.
    pub fn from_real(re: ExactNum, p: usize) -> Self {
        Self::new(re, ExactNum::new(p))
    }

    /// `1 / self`.
    ///
    /// Accuracy: correctly-rounded
    pub fn reciprocal(&self, p: usize, rm: RoundingMode) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.reciprocal_word(p, rm))
    }

    /// Implementation of [`Self::reciprocal`] on operands without inexact history.
    fn reciprocal_word(&self, p: usize, rm: RoundingMode) -> Self {
        Self::one(p).div(self, p, rm)
    }

    /// True if either part is NaN.
    pub fn is_nan(&self) -> bool {
        self.re.is_nan() || self.im.is_nan()
    }

    /// Complex conjugate.
    pub fn conj(&self) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.conj_word())
    }

    /// Implementation of [`Self::conj`] on operands without inexact history.
    fn conj_word(self) -> Self {
        Self::new(self.re, self.im.neg())
    }

    /// Modulus `|z| = hypot(re, im)` at precision `p`.
    ///
    /// Accuracy: correctly-rounded
    pub fn abs(&self, p: usize, rm: RoundingMode) -> ExactNum {
        self.re.hypot(&self.im, p, rm)
    }

    /// Argument `atan2(im, re)` at precision `p`.
    ///
    /// Branch: same as real `atan2`; values lie in (−π, π]. The cut of `ln` / `sqrt` /
    /// `pow` is the non-positive real axis, approached from above as +π and from below as −π.
    ///
    /// Accuracy: correctly-rounded
    pub fn arg(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> ExactNum {
        self.im.atan2(&self.re, p, rm, cc)
    }

    /// `a + b` with both operands borrowed: the sanctioned form for a site where both operands are
    /// still needed afterwards (see [`ExactNum::add_ref`]). Parts and flags are those of [`Self::add`];
    /// no operand is copied, a subnormal one included (the kernels read it as stored and renormalize the exact result; only `div` shifts a subnormal divisor into its own working buffer, because the division kernel needs a normalized divisor).
    ///
    /// Accuracy: correctly-rounded
    pub fn add_ref(a: &Self, b: &Self, p: usize, rm: RoundingMode) -> Self {
        let inputs = a.inexact_accumulated() || b.inexact_accumulated();
        let mut r = Self::new(
            ExactNum::add_ref(&a.re, &b.re, p, rm),
            ExactNum::add_ref(&a.im, &b.im, p, rm),
        );
        r.re.set_inexact(r.re.inexact());
        r.im.set_inexact(r.im.inexact());
        Self::per_op(inputs, r)
    }

    /// `a - b` with both operands borrowed. Same contract as [`Self::add_ref`].
    ///
    /// Accuracy: correctly-rounded
    pub fn sub_ref(a: &Self, b: &Self, p: usize, rm: RoundingMode) -> Self {
        let inputs = a.inexact_accumulated() || b.inexact_accumulated();
        let mut r = Self::new(
            ExactNum::sub_ref(&a.re, &b.re, p, rm),
            ExactNum::sub_ref(&a.im, &b.im, p, rm),
        );
        r.re.set_inexact(r.re.inexact());
        r.im.set_inexact(r.im.inexact());
        Self::per_op(inputs, r)
    }

    /// `self + rhs` at precision `p`. Takes `self` by value: each part of the result is written by
    /// the corresponding part of `self` (see [`ExactNum::add`]).
    ///
    /// Accuracy: correctly-rounded
    pub fn add(self, rhs: &Self, p: usize, rm: RoundingMode) -> Self {
        let inputs = self.inexact_accumulated() || rhs.inexact_accumulated();
        let (re, im) = self.into_parts();
        let mut r = Self::new(re.add(&rhs.re, p, rm), im.add(&rhs.im, p, rm));
        r.re.set_inexact(r.re.inexact());
        r.im.set_inexact(r.im.inexact());
        Self::per_op(inputs, r)
    }

    /// `self - rhs` at precision `p`. Takes `self` by value.
    ///
    /// Accuracy: correctly-rounded
    pub fn sub(self, rhs: &Self, p: usize, rm: RoundingMode) -> Self {
        let inputs = self.inexact_accumulated() || rhs.inexact_accumulated();
        let (re, im) = self.into_parts();
        let mut r = Self::new(re.sub(&rhs.re, p, rm), im.sub(&rhs.im, p, rm));
        r.re.set_inexact(r.re.inexact());
        r.im.set_inexact(r.im.inexact());
        Self::per_op(inputs, r)
    }

    /// `self * rhs` at precision `p`.
    ///
    /// Each of `ac`, `bd`, `ad`, `bc` is rounded at `(p, rm)`, then `ac−bd` and
    /// `ad+bc` are rounded at `(p, rm)`. That matches the MPFR componentwise gold.
    /// Callers that pass [`RoundingMode::None`] still keep full products (series paths).
    ///
    /// Accuracy: correctly-rounded
    pub fn mul(&self, rhs: &Self, p: usize, rm: RoundingMode) -> Self {
        let this = self.exact_copy();
        let rhs_c = rhs.exact_copy();
        let inputs = self.inexact_accumulated() || rhs.inexact_accumulated();
        Self::per_op(inputs, this.mul_word(&rhs_c, p, rm))
    }

    /// Implementation of [`Self::mul`] on operands without inexact history.
    fn mul_word(self, rhs: &Self, p: usize, rm: RoundingMode) -> Self {
        let ac = ExactNum::mul_ref(&self.re, &rhs.re, p, rm);
        let bd = ExactNum::mul_ref(&self.im, &rhs.im, p, rm);
        let ad = self.re.mul(&rhs.im, p, rm);
        let bc = self.im.mul(&rhs.re, p, rm);
        // A product that rounded to a zero part is dropped by the zero-aware sum's flags.
        let products_rounded = ac.inexact() || bd.inexact() || ad.inexact() || bc.inexact();
        let mut r = Self::new(ac.sub(&bd, p, rm), ad.add(&bc, p, rm));
        if products_rounded {
            r.set_inexact(true);
        }
        r
    }

    /// `self / rhs` at precision `p`.
    ///
    /// Accuracy: correctly-rounded
    pub fn div(&self, rhs: &Self, p: usize, rm: RoundingMode) -> Self {
        let this = self.exact_copy();
        let rhs_c = rhs.exact_copy();
        let inputs = self.inexact_accumulated() || rhs.inexact_accumulated();
        Self::per_op(inputs, this.div_word(&rhs_c, p, rm))
    }

    /// Implementation of [`Self::div`] on operands without inexact history.
    ///
    /// The products and the sums `ac + bd`, `bc - ad`, `c^2 + d^2` are formed exactly, so each part
    /// is the exact quotient rounded once by the final division. Rounding the numerator first
    /// rounded twice: `(-7i)/(-7i)` at `p = 1` gave `32/49` rounded, not `1`.
    fn div_word(self, rhs: &Self, p: usize, rm: RoundingMode) -> Self {
        let ac = self.re.mul_full_prec(&rhs.re);
        let bd = self.im.mul_full_prec(&rhs.im);
        let bc = self.im.mul_full_prec(&rhs.re);
        let ad = self.re.mul_full_prec(&rhs.im);
        let cc = rhs.re.mul_full_prec(&rhs.re);
        let dd = rhs.im.mul_full_prec(&rhs.im);
        let den = exact_sum(&cc, &dd, false, p);
        Self::new(
            exact_sum(&ac, &bd, false, p).div(&den, p, rm),
            exact_sum(&bc, &ad, true, p).div(&den, p, rm),
        )
    }

    /// `e^self` using `exp(re) (cos(im) + i sin(im))`.
    ///
    /// Accuracy: correctly-rounded
    pub fn exp(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.exp_word(p, rm, cc))
    }

    /// Implementation of [`Self::exp`] on operands without inexact history.
    fn exp_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let er = self.re.exp(p, RoundingMode::None, cc);
        let (s, c) = self.im.sin_cos(p, RoundingMode::None, cc);
        Self::new(ExactNum::mul_ref(&er, &c, p, rm), er.mul(&s, p, rm))
    }

    /// Principal logarithm `ln|z| + i Arg(z)`.
    ///
    /// Branch cut: (−∞, 0] on the real axis. `ln(−1)` is `iπ` (argument +π).
    /// `|z|` is `ln(max(|re|,|im|)) + ½ ln(1 + (min/max)²)`, so a point just off
    /// the cut keeps the small real part.
    ///
    /// Accuracy: correctly-rounded
    pub fn ln(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.ln_word(p, rm, cc))
    }

    /// Implementation of [`Self::ln`] on operands without inexact history.
    fn ln_word(self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.re.is_nan()
            || self.im.is_nan()
            || self.re.is_inf()
            || self.im.is_inf()
            || (self.re.is_zero() && self.im.is_zero())
        {
            return self.ln_mod_arg(p, rm, cc);
        }
        let prec = round_p(p);
        let work = prec.saturating_mul(2).saturating_add(256);
        // The argument only borrows `self`; compute it before the parts are moved.
        let arg = self.arg(p, rm, cc);
        let re_abs = abs_part(self.re);
        let im_abs = abs_part(self.im);
        let (big, small) = match re_abs.abs_cmp(&im_abs) {
            Some(ord) if ord >= 0 => (re_abs, im_abs),
            Some(_) => (im_abs, re_abs),
            None => {
                return Self::new(
                    ExactNum::nan(Some(crate::Error::InvalidArgument)),
                    ExactNum::nan(Some(crate::Error::InvalidArgument)),
                );
            }
        };
        let ln_big = big.ln(work, RoundingMode::None, cc);
        let ratio = small.div(&big, work, RoundingMode::None);
        let corr = ratio
            .sqr(work, RoundingMode::None)
            .log1p(work, RoundingMode::None, cc)
            .ldexp(-1, work, RoundingMode::None);
        let mut real = ln_big.add(&corr, work, RoundingMode::None);
        if let Err(err) = real.set_precision(prec, rm) {
            return Self::new(ExactNum::nan(Some(err)), ExactNum::nan(Some(err)));
        }
        Self::new(real, arg)
    }

    fn ln_mod_arg(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let mag = self.abs(p, RoundingMode::None);
        Self::new(mag.ln(p, rm, cc), self.arg(p, rm, cc))
    }

    /// `sin(self)` via `sin(re)cosh(im) + i cos(re)sinh(im)`.
    ///
    /// The complex value is **not** passed to `rem_pi`. Only the real (resp. imaginary)
    /// *component* uses real `sin_cos` / `sinh_cosh`.
    ///
    /// Accuracy: correctly-rounded
    pub fn sin(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.sin_word(p, rm, cc))
    }

    /// Implementation of [`Self::sin`] on operands without inexact history.
    fn sin_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let (sn, cs) = self.re.sin_cos(p, RoundingMode::None, cc);
        let (sh, ch) = self.im.sinh_cosh(p, RoundingMode::None, cc);
        Self::new(sn.mul(&ch, p, rm), cs.mul(&sh, p, rm))
    }

    /// `cos(self)` via `cos(re)cosh(im) - i sin(re)sinh(im)`.
    ///
    /// Accuracy: correctly-rounded
    pub fn cos(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.cos_word(p, rm, cc))
    }

    /// Implementation of [`Self::cos`] on operands without inexact history.
    fn cos_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let (sn, cs) = self.re.sin_cos(p, RoundingMode::None, cc);
        let (sh, ch) = self.im.sinh_cosh(p, RoundingMode::None, cc);
        Self::new(cs.mul(&ch, p, rm), sn.mul(&sh, p, rm).neg())
    }

    /// The square root when it is exact at `p_x` bits (`sqrt(3+4i) = 2+i`), found with the
    /// algebraic form `t = sqrt((|z|+|re|)/2)`, `im / (2t)`; `None` when any step rounds.
    fn sqrt_exact(&self, p_x: usize) -> Option<Self> {
        let finite = |x: &ExactNum| !x.is_nan() && !x.is_inf();
        if !finite(&self.re) || !finite(&self.im) || self.im.is_zero() {
            return None;
        }
        let rm = RoundingMode::ToEven;
        let exact = |x: &ExactNum| finite(x) && !x.inexact_accumulated();
        let a = self.re.hypot(&self.im, p_x, rm);
        let a_exact = exact(&a);
        let h = if self.re.is_negative() {
            a.sub(&self.re, p_x, rm)
        } else {
            a.add(&self.re, p_x, rm)
        }
        .ldexp(-1, p_x, rm);
        let t = h.sqrt(p_x, rm);
        // 2t, recomputed from `h` because `t` is still read below
        let t2 = h.sqrt(p_x, rm).ldexp(1, p_x, rm);
        let u = if self.im.is_negative() {
            self.im.neg_via_mul().div(&t2, p_x, rm)
        } else {
            ExactNum::div_ref(&self.im, &t2, p_x, rm)
        };
        if !(a_exact && exact(&h) && exact(&t) && exact(&u)) {
            return None;
        }
        Some(if self.re.is_negative() {
            let t = if self.im.is_negative() { t.neg() } else { t };
            Self::new(u, t)
        } else {
            let u = if self.im.is_negative() { u.neg() } else { u };
            Self::new(t, u)
        })
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
    ///
    /// Accuracy: correctly-rounded
    pub fn tan(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.tan_word(p, rm, cc))
    }

    /// Implementation of [`Self::tan`] on operands without inexact history.
    fn tan_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
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
    ///
    /// Accuracy: correctly-rounded
    pub fn sinh(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.sinh_word(p, rm, cc))
    }

    /// Implementation of [`Self::sinh`] on operands without inexact history.
    fn sinh_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
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
    ///
    /// Accuracy: correctly-rounded
    pub fn cosh(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.cosh_word(p, rm, cc))
    }

    /// Implementation of [`Self::cosh`] on operands without inexact history.
    fn cosh_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
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
    ///
    /// Accuracy: correctly-rounded
    pub fn tanh(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.tanh_word(p, rm, cc))
    }

    /// Implementation of [`Self::tanh`] on operands without inexact history.
    fn tanh_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
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
    ///
    /// Accuracy: correctly-rounded
    pub fn sqrt(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.sqrt_word(p, rm, cc))
    }

    /// Implementation of [`Self::sqrt`] on operands without inexact history.
    fn sqrt_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        if let Some(z) = self.sqrt_exact(p_x) {
            return z.finish(p, rm);
        }
        let r = self.abs(p_x, RoundingMode::None);
        let th = self.arg(p_x, RoundingMode::None, cc);
        let sr = r.sqrt(p_x, RoundingMode::None);
        let half_th = th.ldexp(-1, p_x, RoundingMode::None);
        let (s, c) = half_th.sin_cos(p_x, RoundingMode::None, cc);
        Self::new(
            ExactNum::mul_ref(&sr, &c, p_x, RoundingMode::None),
            sr.mul(&s, p_x, RoundingMode::None),
        )
        .finish(p, rm)
    }

    /// `log2(self) = ln(self) / ln 2` (principal branch).
    ///
    /// Accuracy: correctly-rounded
    pub fn log2(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.log2_word(p, rm, cc))
    }

    /// Implementation of [`Self::log2`] on operands without inexact history.
    fn log2_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        let ln = self.ln(p_x, RoundingMode::None, cc);
        let base = Self::from_real(cc.ln_2(p_x, RoundingMode::None), p_x);
        ln.div(&base, p_x, RoundingMode::None).finish(p, rm)
    }

    /// `log10(self) = ln(self) / ln 10` (principal branch).
    ///
    /// Accuracy: correctly-rounded
    pub fn log10(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.log10_word(p, rm, cc))
    }

    /// Implementation of [`Self::log10`] on operands without inexact history.
    fn log10_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        let ln = self.ln(p_x, RoundingMode::None, cc);
        let base = Self::from_real(cc.ln_10(p_x, RoundingMode::None), p_x);
        ln.div(&base, p_x, RoundingMode::None).finish(p, rm)
    }

    /// `log_base(self) = ln(self) / ln(base)` (principal branch).
    ///
    /// Accuracy: correctly-rounded
    pub fn log(&self, base: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let base_c = base.exact_copy();
        let inputs = self.inexact_accumulated() || base.inexact_accumulated();
        Self::per_op(inputs, this.log_word(&base_c, p, rm, cc))
    }

    /// Implementation of [`Self::log`] on operands without inexact history.
    fn log_word(&self, base: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        let ln = self.ln(p_x, RoundingMode::None, cc);
        let lnb = base.ln(p_x, RoundingMode::None, cc);
        ln.div(&lnb, p_x, RoundingMode::None).finish(p, rm)
    }

    /// `ln(1 + self)` (principal branch).
    ///
    /// Accuracy: correctly-rounded
    pub fn log1p(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.log1p_word(p, rm, cc))
    }

    /// Implementation of [`Self::log1p`] on operands without inexact history.
    fn log1p_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        Self::one(p_x)
            .add(self, p_x, RoundingMode::None)
            .ln(p_x, RoundingMode::None, cc)
            .finish(p, rm)
    }

    /// `exp2(self) = exp(self · ln 2)`.
    ///
    /// Accuracy: correctly-rounded
    pub fn exp2(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.exp2_word(p, rm, cc))
    }

    /// Implementation of [`Self::exp2`] on operands without inexact history.
    fn exp2_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        let ln2 = Self::from_real(cc.ln_2(p_x, RoundingMode::None), p_x);
        self.mul(&ln2, p_x, RoundingMode::None)
            .exp(p_x, RoundingMode::None, cc)
            .finish(p, rm)
    }

    /// `exp10(self) = exp(self · ln 10)`.
    ///
    /// Accuracy: correctly-rounded
    pub fn exp10(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.exp10_word(p, rm, cc))
    }

    /// Implementation of [`Self::exp10`] on operands without inexact history.
    fn exp10_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        let ln10 = Self::from_real(cc.ln_10(p_x, RoundingMode::None), p_x);
        self.mul(&ln10, p_x, RoundingMode::None)
            .exp(p_x, RoundingMode::None, cc)
            .finish(p, rm)
    }

    /// `exp(self) − 1`.
    ///
    /// Accuracy: correctly-rounded
    pub fn expm1(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.expm1_word(p, rm, cc))
    }

    /// Implementation of [`Self::expm1`] on operands without inexact history.
    fn expm1_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        self.exp(p_x, RoundingMode::None, cc)
            .sub(&Self::one(p_x), p_x, RoundingMode::None)
            .finish(p, rm)
    }

    /// Scale both parts by `2^n` (`ldexp` on re and im).
    ///
    /// Accuracy: correctly-rounded
    pub fn ldexp(&self, n: Exponent, p: usize, rm: RoundingMode) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.ldexp_word(n, p, rm))
    }

    /// Implementation of [`Self::ldexp`] on operands without inexact history.
    fn ldexp_word(self, n: Exponent, p: usize, rm: RoundingMode) -> Self {
        self.ldexp_parts(n, p, rm)
    }

    /// Same as [`ldexp`](Self::ldexp).
    ///
    /// Accuracy: correctly-rounded
    pub fn scalb(&self, n: Exponent, p: usize, rm: RoundingMode) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.scalb_word(n, p, rm))
    }

    /// Implementation of [`Self::scalb`] on operands without inexact history.
    fn scalb_word(self, n: Exponent, p: usize, rm: RoundingMode) -> Self {
        self.ldexp(n, p, rm)
    }

    /// `logb(|z|)` as a real (`x + 0i`).
    ///
    /// Accuracy: correctly-rounded
    pub fn logb(&self, p: usize, rm: RoundingMode) -> ExactNum {
        self.abs(p, rm).logb(p, rm)
    }

    /// Principal `n`-th root via `exp(ln(z) / n)`. Inherits the `ln` branch cut.
    ///
    /// Accuracy: correctly-rounded
    pub fn nth_root(&self, n: usize, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.nth_root_word(n, p, rm, cc))
    }

    /// Implementation of [`Self::nth_root`] on operands without inexact history.
    fn nth_root_word(&self, n: usize, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
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
    ///
    /// Accuracy: correctly-rounded
    pub fn cbrt(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.cbrt_word(p, rm, cc))
    }

    /// Implementation of [`Self::cbrt`] on operands without inexact history.
    fn cbrt_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.nth_root(3, p, rm, cc)
    }

    /// Principal `sqrt(self² + other²)` (analytic continuation of real `hypot`).
    ///
    /// Accuracy: correctly-rounded
    pub fn hypot(&self, other: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let other_c = other.exact_copy();
        let inputs = self.inexact_accumulated() || other.inexact_accumulated();
        Self::per_op(inputs, this.hypot_word(&other_c, p, rm, cc))
    }

    /// Implementation of [`Self::hypot`] on operands without inexact history.
    fn hypot_word(&self, other: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
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
    ///
    /// Accuracy: correctly-rounded
    pub fn fma(&self, b: &Self, c: &Self, p: usize, rm: RoundingMode) -> Self {
        let this = self.exact_copy();
        let b_c = b.exact_copy();
        let c_c = c.exact_copy();
        let inputs =
            self.inexact_accumulated() || b.inexact_accumulated() || c.inexact_accumulated();
        Self::per_op(inputs, this.fma_word(&b_c, &c_c, p, rm))
    }

    /// Implementation of [`Self::fma`] on operands without inexact history.
    fn fma_word(&self, b: &Self, c: &Self, p: usize, rm: RoundingMode) -> Self {
        let p_x = Self::work_p(p);
        self.mul(b, p_x, RoundingMode::None)
            .add(c, p_x, RoundingMode::None)
            .finish(p, rm)
    }

    /// Alias of [`fma`](Self::fma).
    ///
    /// Accuracy: correctly-rounded
    pub fn mul_add(&self, b: &Self, c: &Self, p: usize, rm: RoundingMode) -> Self {
        let this = self.exact_copy();
        let b_c = b.exact_copy();
        let c_c = c.exact_copy();
        let inputs =
            self.inexact_accumulated() || b.inexact_accumulated() || c.inexact_accumulated();
        Self::per_op(inputs, this.mul_add_word(&b_c, &c_c, p, rm))
    }

    /// Implementation of [`Self::mul_add`] on operands without inexact history.
    fn mul_add_word(&self, b: &Self, c: &Self, p: usize, rm: RoundingMode) -> Self {
        self.fma(b, c, p, rm)
    }

    /// `self^rhs` as `exp(rhs * ln(self))` (principal branch).
    ///
    /// Inherits the `ln` cut on `self`: non-positive real base uses Arg = ±π.
    ///
    /// Accuracy: correctly-rounded
    pub fn pow(&self, rhs: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let rhs_c = rhs.exact_copy();
        let inputs = self.inexact_accumulated() || rhs.inexact_accumulated();
        Self::per_op(inputs, this.pow_word(&rhs_c, p, rm, cc))
    }

    /// Implementation of [`Self::pow`] on operands without inexact history.
    fn pow_word(&self, rhs: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        let ln = self.ln(p_x, RoundingMode::None, cc);
        rhs.mul(&ln, p_x, RoundingMode::None)
            .exp(p_x, RoundingMode::None, cc)
            .finish(p, rm)
    }

    /// Principal `asin`: `-i ln(i z + √(1 − z²))`.
    ///
    /// On the real cut `|re| > 1`, an imaginary part of `+0` is the principal
    /// value and an imaginary part of `-0` is its conjugate.
    ///
    /// Accuracy: correctly-rounded
    pub fn asin(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.asin_word(p, rm, cc))
    }

    /// Implementation of [`Self::asin`] on operands without inexact history.
    fn asin_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let prec = round_p(p);
        let p_x = prec.saturating_add(WORD_BIT_SIZE);
        let z2 = self.mul(self, p_x, RoundingMode::None);
        let one = Self::one(p_x);
        let rad = one
            .sub(&z2, p_x, RoundingMode::None)
            .sqrt(p_x, RoundingMode::None, cc);
        let iz = Self::i(p_x).mul(self, p_x, RoundingMode::None);
        let ln = iz
            .add(&rad, p_x, RoundingMode::None)
            .ln_mod_arg(p_x, RoundingMode::None, cc);
        let mut v = Self::i(p_x).mul(&ln, p_x, RoundingMode::None).neg();
        if self.asin_below_cut() {
            v = v.conj();
        }
        v.finish(prec, rm)
    }

    /// `|re| > 1` and the imaginary part is negative zero.
    fn asin_below_cut(&self) -> bool {
        if !(self.im.is_zero() && self.im.is_negative()) {
            return false;
        }
        let one = ExactNum::from_u8(1, WORD_BIT_SIZE);
        matches!(self.re.abs_cmp(&one), Some(ord) if ord > 0)
    }

    /// Principal `acos`: `π/2 − asin(z)`.
    ///
    /// Accuracy: correctly-rounded
    pub fn acos(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.acos_word(p, rm, cc))
    }

    /// Implementation of [`Self::acos`] on operands without inexact history.
    fn acos_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
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
    ///
    /// Accuracy: correctly-rounded
    pub fn atan(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.atan_word(p, rm, cc))
    }

    /// Implementation of [`Self::atan`] on operands without inexact history.
    fn atan_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        let num = Self::i(p_x).add(self, p_x, RoundingMode::None);
        let den = Self::i(p_x).sub(self, p_x, RoundingMode::None);
        let ln = num
            .div(&den, p_x, RoundingMode::None)
            .ln(p_x, RoundingMode::None, cc);
        let half_i = Self::i(p_x).ldexp_parts(-1, p_x, RoundingMode::None);
        half_i.mul(&ln, p_x, RoundingMode::None).finish(p, rm)
    }

    /// Principal `asinh`: `ln(z + √(z² + 1))`.
    ///
    /// Accuracy: correctly-rounded
    pub fn asinh(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.asinh_word(p, rm, cc))
    }

    /// Implementation of [`Self::asinh`] on operands without inexact history.
    fn asinh_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        let z2 = self.mul(self, p_x, RoundingMode::None);
        let rad =
            z2.add(&Self::one(p_x), p_x, RoundingMode::None)
                .sqrt(p_x, RoundingMode::None, cc);
        rad.add(self, p_x, RoundingMode::None)
            .ln(p_x, RoundingMode::None, cc)
            .finish(p, rm)
    }

    /// Principal `acosh`: `ln(z + √(z−1)√(z+1))`.
    ///
    /// On the cut `z = x ± 0i`, `x < 1`, the result is `i·acos(x)` for `-1 ≤ x < 1` and
    /// `acosh(-x) + iπ` for `x < -1`, conjugated when the zero imaginary part is negative
    /// (C99 `cacosh` convention; the real part is never negative).
    ///
    /// Accuracy: correctly-rounded
    pub fn acosh(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.acosh_word(p, rm, cc))
    }

    /// Implementation of [`Self::acosh`] on operands without inexact history.
    fn acosh_word(self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        if self.im.is_zero() && !self.re.is_nan() && !self.re.is_inf() {
            let one = ExactNum::from_word(1, p_x);
            if matches!(self.re.compare(&one), Some(c) if c < 0) {
                let neg_one = one.neg();
                let im_negative = self.im.is_negative();
                let (re, im) = if matches!(self.re.compare(&neg_one), Some(c) if c >= 0) {
                    (
                        ExactNum::new(p_x),
                        self.re.acos(p_x, RoundingMode::None, cc),
                    )
                } else {
                    (
                        self.re.neg().acosh(p_x, RoundingMode::None, cc),
                        cc.pi(p_x, RoundingMode::None),
                    )
                };
                let im = if im_negative { im.neg() } else { im };
                return Self::new(re, im).finish(p, rm);
            }
        }
        // `z - 1` is `-1 + z` and `z + 1` is `1 + z`: the constants own the sums, `z` stays borrowed.
        let zm = Self::one(p_x)
            .neg()
            .add(&self, p_x, RoundingMode::None)
            .sqrt(p_x, RoundingMode::None, cc);
        let zp =
            Self::one(p_x)
                .add(&self, p_x, RoundingMode::None)
                .sqrt(p_x, RoundingMode::None, cc);
        zm.mul(&zp, p_x, RoundingMode::None)
            .add(&self, p_x, RoundingMode::None)
            .ln(p_x, RoundingMode::None, cc)
            .finish(p, rm)
    }

    /// Principal `atanh`: `(1/2) ln((1+z)/(1−z))`.
    ///
    /// `atanh(±1 + 0i) = ±∞ + 0i` (C99 `catanh`).
    ///
    /// Accuracy: correctly-rounded
    pub fn atanh(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.atanh_word(p, rm, cc))
    }

    /// Implementation of [`Self::atanh`] on operands without inexact history.
    fn atanh_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let p_x = Self::work_p(p);
        if self.im.is_zero() && !self.re.is_nan() {
            let one_r = ExactNum::from_word(1, p_x);
            let im0 = || {
                let mut z = ExactNum::new(p);
                if self.im.is_negative() {
                    z = z.neg();
                }
                z
            };
            if self.re.compare(&one_r) == Some(0) {
                return Self::new(crate::ext::INF_POS, im0());
            }
            if self.re.compare(&one_r.neg()) == Some(0) {
                return Self::new(crate::ext::INF_NEG, im0());
            }
        }
        let num = Self::one(p_x).add(self, p_x, RoundingMode::None);
        let den = Self::one(p_x).sub(self, p_x, RoundingMode::None);
        num.div(&den, p_x, RoundingMode::None)
            .ln(p_x, RoundingMode::None, cc)
            .ldexp_parts(-1, p_x, RoundingMode::None)
            .finish(p, rm)
    }

    fn ldexp_parts(self, n: crate::Exponent, p: usize, rm: RoundingMode) -> Self {
        Self::new(self.re.ldexp(n, p, rm), self.im.ldexp(n, p, rm))
    }

    pub(crate) fn neg(self) -> Self {
        Self::new(self.re.neg(), self.im.neg())
    }

    /// `-self` for an operand that has to stay borrowed (see [`ExactNum::neg_via_mul`]).
    pub(crate) fn neg_via_mul(&self) -> Self {
        Self::new(self.re.neg_via_mul(), self.im.neg_via_mul())
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

// `add` and `sub` take their first operand by value, so there is no `&ExactComplex + &ExactComplex`
// (see the operators on `ExactNum`): `&a + b` is `b + a`, `&a - b` is `-b + a`.
impl core::ops::Add<&ExactComplex> for ExactComplex {
    type Output = ExactComplex;
    fn add(self, rhs: &ExactComplex) -> ExactComplex {
        ExactComplex::add(self, rhs, DEFAULT_P, RoundingMode::ToEven)
    }
}
impl core::ops::Add<ExactComplex> for ExactComplex {
    type Output = ExactComplex;
    fn add(self, rhs: ExactComplex) -> ExactComplex {
        ExactComplex::add(self, &rhs, DEFAULT_P, RoundingMode::ToEven)
    }
}
impl core::ops::Add<ExactComplex> for &ExactComplex {
    type Output = ExactComplex;
    fn add(self, rhs: ExactComplex) -> ExactComplex {
        ExactComplex::add(rhs, self, DEFAULT_P, RoundingMode::ToEven)
    }
}
impl core::ops::Sub<&ExactComplex> for ExactComplex {
    type Output = ExactComplex;
    fn sub(self, rhs: &ExactComplex) -> ExactComplex {
        ExactComplex::sub(self, rhs, DEFAULT_P, RoundingMode::ToEven)
    }
}
impl core::ops::Sub<ExactComplex> for ExactComplex {
    type Output = ExactComplex;
    fn sub(self, rhs: ExactComplex) -> ExactComplex {
        ExactComplex::sub(self, &rhs, DEFAULT_P, RoundingMode::ToEven)
    }
}
impl core::ops::Sub<ExactComplex> for &ExactComplex {
    type Output = ExactComplex;
    fn sub(self, rhs: ExactComplex) -> ExactComplex {
        ExactComplex::add(rhs.neg(), self, DEFAULT_P, RoundingMode::ToEven)
    }
}
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
    fn div_rounds_the_exact_quotient_once() {
        // (-7i) / (-7i) = 1 exactly; rounding the numerator 49 to one bit first gave 32/49.
        let m7i = || ExactComplex::new(ExactNum::new(64), -ExactNum::from_u8(7, 64));
        for rm in [RoundingMode::ToZero, RoundingMode::Down, RoundingMode::Up, RoundingMode::ToEven] {
            let q = m7i().div(&m7i(), 1, rm);
            assert_eq!(q.re().compare(&ExactNum::from_u8(1, 64)), Some(0), "{rm:?}");
            assert!(q.im().is_zero(), "{rm:?}");
            assert!(!q.inexact(), "{rm:?}");
        }
    }

    #[test]
    fn div_keeps_the_direction_of_a_far_smaller_numerator_term() -> Result<(), crate::Error> {
        // (-7 + s i) / (-7 - 7i) with s = 1.5 * 2^(EMIN - 18): im = (-49 - 7s) / 98 = -1/2 - 7s/98, just below -1/2.
        // Rounding the wide sum `bc - ad` without a sticky bit gave exactly -1/2, so Down kept -1/2.
        let m7 = || -ExactNum::from_u8(7, 64);
        let s = ExactNum::from_mag_parts(0xC000_0000_0000_0000, i64::from(crate::EXPONENT_MIN) - 17, 64)?;
        let q = ExactComplex::new(m7(), s).div(&ExactComplex::new(m7(), m7()), 1, RoundingMode::Down);
        assert_eq!(q.im().compare(&-ExactNum::from_u8(1, 64)), Some(0));
        assert!(q.inexact());
        Ok(())
    }

    #[test]
    fn mul_reports_a_product_that_rounded_to_a_zero_part() {
        // im = (-7)(tiny) + 0 * 0: the product underflows to -0 at p = 1, so the operation rounded.
        let tiny = ExactNum::from_u8(3, 64)
            .ldexp(crate::EXPONENT_MIN + 1, 64, RoundingMode::None)
            .ldexp(-20, 64, RoundingMode::None);
        let x = ExactComplex::new(-ExactNum::from_u8(7, 64), ExactNum::new(64));
        let y = ExactComplex::new(ExactNum::new(64), tiny);
        let z = x.mul(&y, 1, RoundingMode::ToEven);
        assert!(z.im().is_zero());
        assert!(z.inexact());
        assert!(z.inexact_accumulated());
    }

    #[test]
    fn test_complex_arith() -> Result<(), crate::Error> {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;

        let a = ExactComplex::new(ExactNum::from_u8(3, p), ExactNum::from_u8(4, p));
        let mag = a.abs(p, rm);
        let five = ExactNum::from_u8(5, p);
        assert_eq!(mag.compare(&five), Some(0));

        let i = ExactComplex::i(p);
        let i2 = i.mul(&i, p, rm);
        assert_eq!(i2.re().compare(&ExactNum::from_i8(-1, p)), Some(0));
        assert!(i2.im().is_zero());

        let z = ExactComplex::one(p);
        let e = z.exp(p, rm, &mut cc);
        let back = e.ln(p, rm, &mut cc);
        let d = ExactNum::sub_ref(back.re(), z.re(), p, RoundingMode::None).abs();
        assert!(d.is_zero() || d.exponent().unwrap_or(0) < -((p as Exponent) / 4));

        let z0 = ExactComplex::zero(p);
        let t = z0.tan(p, rm, &mut cc);
        assert!(t.re().is_zero() || t.re().exponent().unwrap_or(0) < -((p as Exponent) / 4));
        assert!(t.im().is_zero() || t.im().exponent().unwrap_or(0) < -((p as Exponent) / 4));
        let sh = z0.sinh(p, rm, &mut cc);
        assert!(sh.re().is_zero() || sh.re().exponent().unwrap_or(0) < -((p as Exponent) / 4));

        let two = ExactComplex::from_real(ExactNum::from_u8(2, p), p);
        let four = ExactComplex::from_real(ExactNum::from_u8(4, p), p);
        let lg = four.log2(p, rm, &mut cc);
        let d = ExactNum::sub_ref(lg.re(), two.re(), p, RoundingMode::None).abs();
        assert!(d.is_zero() || d.exponent().unwrap_or(0) < -((p as Exponent) / 8));
        assert!(lg.im().is_zero() || lg.im().exponent().unwrap_or(0) < -((p as Exponent) / 4));
        let e2 = two.exp2(p, rm, &mut cc);
        let d2 = ExactNum::sub_ref(e2.re(), four.re(), p, RoundingMode::None).abs();
        assert!(d2.is_zero() || d2.exponent().unwrap_or(0) < -((p as Exponent) / 8));
        Ok(())
    }

    #[test]
    fn test_complex_branch_cuts() -> Result<(), crate::Error> {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let m1 = ExactComplex::from_real(ExactNum::from_i8(-1, p), p);

        let s = m1.sqrt(p, rm, &mut cc);
        assert!(s.re().is_zero() || s.re().exponent().unwrap_or(0) < -((p as Exponent) / 4));
        assert_eq!(s.im().compare(&ExactNum::from_u8(1, p)), Some(0));

        let l = m1.ln(p, rm, &mut cc);
        assert!(l.re().is_zero() || l.re().exponent().unwrap_or(0) < -((p as Exponent) / 4));
        let pi = cc.pi(p, rm);
        let d = (if l.im().is_negative() {
            l.im().neg_via_mul().sub(&pi, p, RoundingMode::None)
        } else {
            ExactNum::sub_ref(l.im(), &pi, p, RoundingMode::None)
        })
        .abs();
        assert!(d.is_zero() || d.exponent().unwrap_or(0) < -((p as Exponent) / 8));
        assert!(l.im().is_positive());

        let below = ExactComplex::new(ExactNum::from_i8(-1, p), ExactNum::new(p).neg());
        let lb = below.ln(p, rm, &mut cc);
        let db = (if lb.im().is_negative() {
            lb.im().neg_via_mul().sub(&pi, p, RoundingMode::None)
        } else {
            pi.neg().add(lb.im(), p, RoundingMode::None)
        })
        .abs();
        assert!(db.is_zero() || db.exponent().unwrap_or(0) < -((p as Exponent) / 8));
        assert!(lb.im().is_negative());
        Ok(())
    }

    #[test]
    fn test_complex_nan() {
        let n = ExactComplex::new(NAN, ExactNum::new(64));
        assert!(n.is_nan());
    }

    #[test]
    fn into_parts_moves_both_parts() {
        let re = ExactNum::from_u8(3, 64);
        let im = ExactNum::from_u8(5, 64);
        let (r, i) = ExactComplex::new(re, im).into_parts();
        assert!(matches!(r.compare(&ExactNum::from_u8(3, 64)), Some(0)));
        assert!(matches!(i.compare(&ExactNum::from_u8(5, 64)), Some(0)));
    }

    #[test]
    fn exponent_range_check_on_complex_parts() {
        let p = 64;
        let big = ExactNum::from_u8(1, p).ldexp(100, p, RoundingMode::ToEven);
        let small = ExactNum::from_u8(1, p).ldexp(-100, p, RoundingMode::ToEven);
        let ok = ExactNum::from_u8(3, p);
        let z = crate::macro_util::check_complex_exponent_range(
            ExactComplex::new(big.neg(), small),
            -50,
            50,
        );
        assert!(z.re().is_inf() && z.re().is_negative());
        assert!(z.im().is_zero());
        let z = crate::macro_util::check_complex_exponent_range(
            ExactComplex::new(ok.neg_via_mul(), ok),
            -50,
            50,
        );
        assert!(matches!(
            z.re().compare(&ExactNum::from_u8(3, p).neg()),
            Some(0)
        ));
        assert!(matches!(z.im().compare(&ExactNum::from_u8(3, p)), Some(0)));
    }

    #[test]
    fn conj_flips_only_the_imaginary_sign() {
        let p = 128;
        let z = ExactComplex::new(ExactNum::from_u8(3, p), ExactNum::from_u8(4, p));
        let c = z.conj();
        assert_eq!(c.re().compare(&ExactNum::from_u8(3, p)), Some(0));
        assert_eq!(c.im().compare(&ExactNum::from_i8(-4, p)), Some(0));
        let back = c.conj();
        assert_eq!(back.im().compare(&ExactNum::from_u8(4, p)), Some(0));
    }

    #[test]
    fn acosh_of_real_below_minus_one_has_pi_imaginary_part_with_signed_sign(
    ) -> Result<(), crate::Error> {
        let p = 128;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let pi = cc.pi(p, rm);
        let up =
            ExactComplex::new(ExactNum::from_i8(-2, p), ExactNum::new(p)).acosh(p, rm, &mut cc);
        assert_eq!(up.im().compare(&pi), Some(0));
        let down = ExactComplex::new(ExactNum::from_i8(-2, p), ExactNum::new(p).neg())
            .acosh(p, rm, &mut cc);
        assert_eq!(down.im().compare(&pi.neg_via_mul()), Some(0));
        assert_eq!(up.re().compare(down.re()), Some(0));
        // |x| < 1: purely imaginary, im = acos(x).
        let inside =
            ExactComplex::new(ExactNum::from_i8(-1, p), ExactNum::new(p)).acosh(p, rm, &mut cc);
        assert!(inside.re().is_zero());
        assert_eq!(inside.im().compare(&pi), Some(0));
        Ok(())
    }

    #[test]
    fn ln_of_second_quadrant_point_uses_the_larger_part_for_the_modulus() -> Result<(), crate::Error>
    {
        let p = 128;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        // ln(-3 + 4i): re = ln 5, |re part| < |im part| exercises the swapped (big, small) arm.
        let z = ExactComplex::new(ExactNum::from_i8(-3, p), ExactNum::from_u8(4, p));
        let l = z.ln(p, rm, &mut cc);
        let ln5 = ExactNum::from_u8(5, p).ln(p, rm, &mut cc);
        let d = ln5.neg().add(l.re(), p, RoundingMode::None).abs();
        assert!(d.is_zero() || d.exponent().unwrap_or(0) < -((p as Exponent) / 2));
        assert!(l.im().is_positive());
        Ok(())
    }
}
