//! ExactNum including finite numbers, NaN, and `Inf`.

use crate::common::util::log2_ceil;
use crate::defs::SignedWord;
use crate::defs::DEFAULT_P;
use crate::num::ExactNumNumber;
use crate::Consts;
use crate::Error;
use crate::Exponent;
use crate::Radix;
use crate::RoundingMode;
use crate::Sign;
use crate::Word;
use crate::WORD_BIT_SIZE;
use core::num::FpCategory;
use lazy_static::lazy_static;

#[cfg(feature = "std")]
use core::fmt::Write;

#[cfg(not(feature = "std"))]
use alloc::{string::String, vec::Vec};

/// Not a number.
pub const NAN: ExactNum = ExactNum {
    inner: Flavor::NaN(None),
};

/// Positive infinity.
pub const INF_POS: ExactNum = ExactNum {
    inner: Flavor::Inf(Sign::Pos),
};

/// Negative infinity.
pub const INF_NEG: ExactNum = ExactNum {
    inner: Flavor::Inf(Sign::Neg),
};

lazy_static! {

    /// 1
    pub static ref ONE: ExactNum = ExactNum { inner: Flavor::Value(ExactNumNumber::from_word(1, DEFAULT_P).expect("Constant ONE initialized")) };

    /// 2
    pub static ref TWO: ExactNum = ExactNum { inner: Flavor::Value(ExactNumNumber::from_word(2, DEFAULT_P).expect("Constant TWO initialized")) };
}

/// A floating point number of arbitrary precision.
#[derive(Debug)]
pub struct ExactNum {
    inner: Flavor,
}

#[derive(Debug)]
enum Flavor {
    Value(ExactNumNumber),
    NaN(Option<Error>),
    Inf(Sign), // signed Inf
}

impl ExactNum {
    /// Returns a new number with value of 0 and precision of `p` bits. Precision is rounded upwards to the word size.
    /// The function returns NaN if the precision `p` is incorrect.
    pub fn new(p: usize) -> Self {
        Self::result_to_ext(ExactNumNumber::new(p), false, true)
    }

    /// Constructs not-a-number with an associated error `err`.
    pub fn nan(err: Option<Error>) -> Self {
        ExactNum {
            inner: Flavor::NaN(err),
        }
    }

    /// Returns true if `self` is positive infinity.
    pub fn is_inf_pos(&self) -> bool {
        matches!(self.inner, Flavor::Inf(Sign::Pos))
    }

    /// Returns true if `self` is negative infinity.
    pub fn is_inf_neg(&self) -> bool {
        matches!(self.inner, Flavor::Inf(Sign::Neg))
    }

    /// Returns true if `self` is infinite.
    pub fn is_inf(&self) -> bool {
        matches!(self.inner, Flavor::Inf(_))
    }

    /// Return true if `self` is not a number.
    pub fn is_nan(&self) -> bool {
        matches!(self.inner, Flavor::NaN(_))
    }

    /// Return true if `self` is an integer number.
    pub fn is_int(&self) -> bool {
        match &self.inner {
            Flavor::Value(v) => v.is_int(),
            Flavor::NaN(_) => false,
            Flavor::Inf(_) => false,
        }
    }

    /// Returns the associated with NaN error, if any.
    pub fn err(&self) -> Option<Error> {
        match &self.inner {
            Flavor::NaN(Some(e)) => Some(*e),
            _ => None,
        }
    }

    /// Adds `d2` to `self` and returns the result of the operation with precision `p` rounded according to `rm`.
    /// Precision is rounded upwards to the word size.
    /// The function returns NaN if the precision `p` is incorrect.
    pub fn add(&self, d2: &Self, p: usize, rm: RoundingMode) -> Self {
        self.add_op(d2, p, rm, false)
    }

    /// Adds `d2` to `self` and returns the result of the operation.
    /// The resulting precision is equal to the full precision of the result.
    /// This operation can be used to emulate integer addition.
    pub fn add_full_prec(&self, d2: &Self) -> Self {
        self.add_op(d2, 0, RoundingMode::None, true)
    }

    fn add_op(&self, d2: &Self, p: usize, rm: RoundingMode, full_prec: bool) -> Self {
        match &self.inner {
            Flavor::Value(v1) => match &d2.inner {
                Flavor::Value(v2) => Self::result_to_ext(
                    if full_prec { v1.add_full_prec(v2) } else { v1.add(v2, p, rm) },
                    v1.is_zero(),
                    v1.sign() == v2.sign(),
                ),
                Flavor::Inf(s2) => ExactNum {
                    inner: Flavor::Inf(*s2),
                },
                Flavor::NaN(err) => Self::nan(*err),
            },
            Flavor::Inf(s1) => match &d2.inner {
                Flavor::Value(_) => ExactNum {
                    inner: Flavor::Inf(*s1),
                },
                Flavor::Inf(s2) => {
                    if *s1 != *s2 {
                        NAN
                    } else {
                        ExactNum {
                            inner: Flavor::Inf(*s2),
                        }
                    }
                }
                Flavor::NaN(err) => Self::nan(*err),
            },
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    /// Subtracts `d2` from `self` and returns the result of the operation with precision `p` rounded according to `rm`.
    /// Precision is rounded upwards to the word size.
    /// The function returns NaN if the precision `p` is incorrect.
    pub fn sub(&self, d2: &Self, p: usize, rm: RoundingMode) -> Self {
        self.sub_op(d2, p, rm, false)
    }

    /// Subtracts `d2` from `self` and returns the result of the operation.
    /// The resulting precision is equal to the full precision of the result.
    /// This operation can be used to emulate integer subtraction.
    pub fn sub_full_prec(&self, d2: &Self) -> Self {
        self.sub_op(d2, 0, RoundingMode::None, true)
    }

    fn sub_op(&self, d2: &Self, p: usize, rm: RoundingMode, full_prec: bool) -> Self {
        match &self.inner {
            Flavor::Value(v1) => match &d2.inner {
                Flavor::Value(v2) => Self::result_to_ext(
                    if full_prec { v1.sub_full_prec(v2) } else { v1.sub(v2, p, rm) },
                    v1.is_zero(),
                    v1.sign() == v2.sign(),
                ),
                Flavor::Inf(s2) => {
                    if s2.is_positive() {
                        INF_NEG
                    } else {
                        INF_POS
                    }
                }
                Flavor::NaN(err) => Self::nan(*err),
            },
            Flavor::Inf(s1) => match &d2.inner {
                Flavor::Value(_) => ExactNum {
                    inner: Flavor::Inf(*s1),
                },
                Flavor::Inf(s2) => {
                    if *s1 == *s2 {
                        NAN
                    } else {
                        ExactNum {
                            inner: Flavor::Inf(*s1),
                        }
                    }
                }
                Flavor::NaN(err) => Self::nan(*err),
            },
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    /// Multiplies `d2` by `self` and returns the result of the operation with precision `p` rounded according to `rm`.
    /// Precision is rounded upwards to the word size.
    /// The function returns NaN if the precision `p` is incorrect.
    pub fn mul(&self, d2: &Self, p: usize, rm: RoundingMode) -> Self {
        self.mul_op(d2, p, rm, false)
    }

    /// Multiplies `d2` by `self` and returns the result of the operation.
    /// The resulting precision is equal to the full precision of the result.
    /// This operation can be used to emulate integer multiplication.
    pub fn mul_full_prec(&self, d2: &Self) -> Self {
        self.mul_op(d2, 0, RoundingMode::None, true)
    }

    /// Computes `self * b + c` with precision `p`, rounded once with `rm`.
    ///
    /// Unlike `mul` followed by `add`, the product is not rounded to `p` before the addition.
    pub fn fma(&self, b: &Self, c: &Self, p: usize, rm: RoundingMode) -> Self {
        if self.is_nan() {
            return self.clone();
        }
        if b.is_nan() {
            return b.clone();
        }
        if c.is_nan() {
            return c.clone();
        }
        match (&self.inner, &b.inner, &c.inner) {
            (Flavor::Value(a), Flavor::Value(bv), Flavor::Value(cv)) => {
                Self::result_to_ext(a.fma(bv, cv, p, rm), false, true)
            }
            _ => {
                let prod = self.mul(b, p, RoundingMode::None);
                prod.add(c, p, rm)
            }
        }
    }

    /// Knuth–Dekker two-sum: `(hi, lo)` with `hi` rounded to `p` bits using `rm` and
    /// `hi + lo` equal to the exact sum of finite operands (via [`add_full_prec`](Self::add_full_prec)).
    /// Unlike a hardware-float Dekker two-sum, this takes `(p, rm)` because the high part is an
    /// `ExactNum` at a chosen precision, not an implicit machine word.
    ///
    /// Inf / NaN: `hi` is `self.add(b, p, rm)`; `lo` is zero (or NaN if `hi` is NaN).
    /// Reconstruct with `hi.add(&lo, p, rm)` (not `add_full_prec`, which uses internal precision 0).
    pub fn two_sum(&self, b: &Self, p: usize, rm: RoundingMode) -> (Self, Self) {
        if self.is_nan() {
            return (self.clone(), Self::nan(self.err()));
        }
        if b.is_nan() {
            return (b.clone(), Self::nan(b.err()));
        }
        if self.is_inf() || b.is_inf() {
            return (self.add(b, p, rm), Self::new(p));
        }
        let exact = self.add_full_prec(b);
        let mut hi = exact.clone();
        if let Err(err) = hi.set_precision(p, rm) {
            return (Self::nan(Some(err)), Self::nan(Some(err)));
        }
        let lo = exact.sub_full_prec(&hi);
        (hi, Self::normalize_eft_lo(lo, p))
    }

    fn normalize_eft_lo(lo: Self, p: usize) -> Self {
        if lo.is_nan() {
            return lo;
        }
        if lo.is_zero() {
            let mut z = Self::new(p);
            z.set_inexact(lo.inexact());
            return z;
        }
        lo
    }

    /// Two-product: `(hi, lo)` with `hi` rounded to `p` bits using `rm` and `hi + lo` equal to the
    /// exact product of finite operands (via [`mul_full_prec`](Self::mul_full_prec)).
    pub fn two_product(&self, b: &Self, p: usize, rm: RoundingMode) -> (Self, Self) {
        if self.is_nan() {
            return (self.clone(), Self::nan(self.err()));
        }
        if b.is_nan() {
            return (b.clone(), Self::nan(b.err()));
        }
        if self.is_inf() || b.is_inf() {
            return (self.mul(b, p, rm), Self::new(p));
        }
        let exact = self.mul_full_prec(b);
        let mut hi = exact.clone();
        if let Err(err) = hi.set_precision(p, rm) {
            return (Self::nan(Some(err)), Self::nan(Some(err)));
        }
        let lo = exact.sub_full_prec(&hi);
        (hi, Self::normalize_eft_lo(lo, p))
    }

    /// Sum `xs` at extra working precision and round once to `p` bits.
    pub fn fused_sum(xs: &[Self], p: usize, rm: RoundingMode) -> Self {
        if xs.is_empty() {
            return Self::new(p);
        }
        let extra = log2_ceil(xs.len().max(1)).saturating_add(2);
        let p_wrk = match p
            .checked_add(WORD_BIT_SIZE)
            .and_then(|v| v.checked_add(extra))
        {
            Some(v) => v,
            None => return Self::nan(Some(Error::InvalidArgument)),
        };
        let mut acc = Self::new(p_wrk);
        for x in xs {
            acc = acc.add(x, p_wrk, RoundingMode::None);
        }
        if let Err(err) = acc.set_precision(p, rm) {
            return Self::nan(Some(err));
        }
        acc
    }

    /// Dot product of equal-length slices: extra-precision `∑ xs[i]*ys[i]`, then one round to `p`.
    /// Length mismatch yields NaN (`InvalidArgument`).
    pub fn fused_dot(xs: &[Self], ys: &[Self], p: usize, rm: RoundingMode) -> Self {
        if xs.len() != ys.len() {
            return Self::nan(Some(Error::InvalidArgument));
        }
        if xs.is_empty() {
            return Self::new(p);
        }
        let extra = log2_ceil(xs.len().max(1)).saturating_add(2);
        let p_wrk = match p
            .checked_add(WORD_BIT_SIZE)
            .and_then(|v| v.checked_add(extra))
        {
            Some(v) => v,
            None => return Self::nan(Some(Error::InvalidArgument)),
        };
        let mut acc = Self::new(p_wrk);
        for (x, y) in xs.iter().zip(ys.iter()) {
            let prod = x.mul(y, p_wrk, RoundingMode::None);
            acc = acc.add(&prod, p_wrk, RoundingMode::None);
        }
        if let Err(err) = acc.set_precision(p, rm) {
            return Self::nan(Some(err));
        }
        acc
    }

    /// Horner evaluation `a₀ + x(a₁ + x(a₂ + …))` with fused multiply-add at extra working precision,
    /// then one round to `p`. `coeffs[0]` is the constant term (lowest degree first).
    /// Empty `coeffs` yields zero.
    pub fn polyval(coeffs: &[Self], x: &Self, p: usize, rm: RoundingMode) -> Self {
        if coeffs.is_empty() {
            return Self::new(p);
        }
        let extra = log2_ceil(coeffs.len().max(1)).saturating_add(2);
        let p_wrk = match p
            .checked_add(WORD_BIT_SIZE)
            .and_then(|v| v.checked_add(extra))
        {
            Some(v) => v,
            None => return Self::nan(Some(Error::InvalidArgument)),
        };
        let mut acc = coeffs[coeffs.len() - 1].clone();
        if let Err(err) = acc.set_precision(p_wrk, RoundingMode::None) {
            return Self::nan(Some(err));
        }
        for a in coeffs.iter().rev().skip(1) {
            acc = acc.fma(x, a, p_wrk, RoundingMode::None);
        }
        if let Err(err) = acc.set_precision(p, rm) {
            return Self::nan(Some(err));
        }
        acc
    }

    /// Alias of [`Self::fma`].
    pub fn mul_add(&self, b: &Self, c: &Self, p: usize, rm: RoundingMode) -> Self {
        if self.is_nan() {
            return self.clone();
        }
        if b.is_nan() {
            return b.clone();
        }
        if c.is_nan() {
            return c.clone();
        }
        match (&self.inner, &b.inner, &c.inner) {
            (Flavor::Value(a), Flavor::Value(bv), Flavor::Value(cv)) => {
                Self::result_to_ext(a.mul_add(bv, cv, p, rm), false, true)
            }
            _ => self.fma(b, c, p, rm),
        }
    }

    fn mul_op(&self, d2: &Self, p: usize, rm: RoundingMode, full_prec: bool) -> Self {
        match &self.inner {
            Flavor::Value(v1) => {
                match &d2.inner {
                    Flavor::Value(v2) => Self::result_to_ext(
                        if full_prec { v1.mul_full_prec(v2) } else { v1.mul(v2, p, rm) },
                        v1.is_zero(),
                        v1.sign() == v2.sign(),
                    ),
                    Flavor::Inf(s2) => {
                        if v1.is_zero() {
                            // 0*inf
                            NAN
                        } else {
                            let s = if v1.sign() == *s2 { Sign::Pos } else { Sign::Neg };
                            ExactNum {
                                inner: Flavor::Inf(s),
                            }
                        }
                    }
                    Flavor::NaN(err) => Self::nan(*err),
                }
            }
            Flavor::Inf(s1) => {
                match &d2.inner {
                    Flavor::Value(v2) => {
                        if v2.is_zero() {
                            // inf*0
                            NAN
                        } else {
                            let s = if v2.sign() == *s1 { Sign::Pos } else { Sign::Neg };
                            ExactNum {
                                inner: Flavor::Inf(s),
                            }
                        }
                    }
                    Flavor::Inf(s2) => {
                        let s = if s1 == s2 { Sign::Pos } else { Sign::Neg };
                        ExactNum {
                            inner: Flavor::Inf(s),
                        }
                    }
                    Flavor::NaN(err) => Self::nan(*err),
                }
            }
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    /// Divides `self` by `d2` and returns the result of the operation with precision `p` rounded according to `rm`.
    /// Precision is rounded upwards to the word size.
    /// The function returns NaN if the precision `p` is incorrect.
    pub fn div(&self, d2: &Self, p: usize, rm: RoundingMode) -> Self {
        match &self.inner {
            Flavor::Value(v1) => match &d2.inner {
                Flavor::Value(v2) => {
                    Self::result_to_ext(v1.div(v2, p, rm), v1.is_zero(), v1.sign() == v2.sign())
                }
                Flavor::Inf(_) => Self::new(v1.mantissa_max_bit_len()),
                Flavor::NaN(err) => Self::nan(*err),
            },
            Flavor::Inf(s1) => match &d2.inner {
                Flavor::Value(v) => {
                    if *s1 == v.sign() {
                        INF_POS
                    } else {
                        INF_NEG
                    }
                }
                Flavor::Inf(_) => NAN,
                Flavor::NaN(err) => Self::nan(*err),
            },
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    /// Returns the remainder of division of `|self|` by `|d2|`. The sign of the result is set to the sign of `self`.
    pub fn rem(&self, d2: &Self) -> Self {
        match &self.inner {
            Flavor::Value(v1) => match &d2.inner {
                Flavor::Value(v2) => {
                    Self::result_to_ext(v1.rem(v2), v1.is_zero(), v1.sign() == v2.sign())
                }
                Flavor::Inf(_) => self.clone(),
                Flavor::NaN(err) => Self::nan(*err),
            },
            Flavor::Inf(_) => NAN,
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    /// Compares `self` to `d2`.
    /// Returns positive if `self` > `d2`, negative if `self` < `d2`, zero if `self` == `d2`, None if `self` or `d2` is NaN.
    #[allow(clippy::should_implement_trait)]
    pub fn cmp(&self, d2: &ExactNum) -> Option<SignedWord> {
        match &self.inner {
            Flavor::Value(v1) => match &d2.inner {
                Flavor::Value(v2) => Some(v1.cmp(v2)),
                Flavor::Inf(s2) => {
                    if *s2 == Sign::Pos {
                        Some(-1)
                    } else {
                        Some(1)
                    }
                }
                Flavor::NaN(_) => None,
            },
            Flavor::Inf(s1) => match &d2.inner {
                Flavor::Value(_) => Some(*s1 as SignedWord),
                Flavor::Inf(s2) => Some(*s1 as SignedWord - *s2 as SignedWord),
                Flavor::NaN(_) => None,
            },
            Flavor::NaN(_) => None,
        }
    }

    /// Compares the absolute value of `self` to the absolute value of `d2`.
    /// Returns positive if `|self|` is greater than `|d2|`, negative if `|self|` is smaller than `|d2|`, 0 if `|self|` equals to `|d2|`, None if `self` or `d2` is NaN.
    pub fn abs_cmp(&self, d2: &Self) -> Option<SignedWord> {
        match &self.inner {
            Flavor::Value(v1) => match &d2.inner {
                Flavor::Value(v2) => Some(v1.cmp(v2)),
                Flavor::Inf(_) => Some(-1),
                Flavor::NaN(_) => None,
            },
            Flavor::Inf(_) => match &d2.inner {
                Flavor::Value(_) => Some(1),
                Flavor::Inf(_) => Some(0),
                Flavor::NaN(_) => None,
            },
            Flavor::NaN(_) => None,
        }
    }

    /// Reverses the sign of `self`.
    pub fn inv_sign(&mut self) {
        match &mut self.inner {
            Flavor::Value(v1) => v1.inv_sign(),
            Flavor::Inf(s) => self.inner = Flavor::Inf(s.invert()),
            Flavor::NaN(_) => {}
        }
    }

    /// Compute the power of `self` to the `n` with precision `p`. The result is rounded using the rounding mode `rm`.
    /// This function requires constants cache `cc` for computing the result.
    /// Precision is rounded upwards to the word size.
    /// The function returns NaN if the precision `p` is incorrect.
    pub fn pow(&self, n: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match &self.inner {
            Flavor::Value(v1) => {
                match &n.inner {
                    Flavor::Value(v2) => Self::result_to_ext(
                        v1.pow(v2, p, rm, cc),
                        v1.is_zero(),
                        v1.sign() == v2.sign(),
                    ),
                    Flavor::Inf(s2) => {
                        // v1^inf
                        let val = v1.cmp(&crate::common::consts::ONE);
                        if val > 0 {
                            ExactNum {
                                inner: Flavor::Inf(*s2),
                            }
                        } else if val < 0 {
                            Self::new(p)
                        } else {
                            Self::from_u8(1, p)
                        }
                    }
                    Flavor::NaN(err) => Self::nan(*err),
                }
            }
            Flavor::Inf(s1) => {
                match &n.inner {
                    Flavor::Value(v2) => {
                        // inf ^ v2
                        if v2.is_zero() {
                            Self::from_u8(1, p)
                        } else if v2.is_positive() {
                            if s1.is_negative() && v2.is_odd_int() {
                                // v2 is odd and has no fractional part.
                                INF_NEG
                            } else {
                                INF_POS
                            }
                        } else {
                            Self::new(p)
                        }
                    }
                    Flavor::Inf(s2) => {
                        // inf^inf
                        if s2.is_positive() {
                            INF_POS
                        } else {
                            Self::new(p)
                        }
                    }
                    Flavor::NaN(err) => Self::nan(*err),
                }
            }
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    /// Compute the power of `self` to the integer `n` with precision `p`. The result is rounded using the rounding mode `rm`.
    /// Precision is rounded upwards to the word size.
    /// The function returns NaN if the precision `p` is incorrect.
    pub fn powi(&self, n: usize, p: usize, rm: RoundingMode) -> Self {
        match &self.inner {
            Flavor::Value(v1) => Self::result_to_ext(v1.powi(n, p, rm), false, true),
            Flavor::Inf(s1) => {
                // inf ^ v2
                if n == 0 {
                    Self::from_u8(1, p)
                } else if s1.is_negative() && (n & 1 == 1) {
                    INF_NEG
                } else {
                    INF_POS
                }
            }
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    /// Compute the power of `self` to the signed integer `n` with precision `p`. The result is rounded using the rounding mode `rm`.
    /// Precision is rounded upwards to the word size.
    /// Negative `n` is a reciprocal of the corresponding positive power.
    /// The function returns NaN if the precision `p` is incorrect, or Inf if `self` is zero and `n` is negative.
    pub fn powsi(&self, n: isize, p: usize, rm: RoundingMode) -> Self {
        match &self.inner {
            Flavor::Value(v1) => Self::result_to_ext(v1.powsi(n, p, rm), false, true),
            Flavor::Inf(s1) => {
                if n == 0 {
                    Self::from_u8(1, p)
                } else if n < 0 {
                    Self::new(p)
                } else if s1.is_negative() && (n & 1 == 1) {
                    INF_NEG
                } else {
                    INF_POS
                }
            }
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    /// Computes the logarithm base `n` of a number with precision `p`. The result is rounded using the rounding mode `rm`.
    /// This function requires constants cache `cc` for computing the result.
    /// Precision is rounded upwards to the word size.
    /// The function returns NaN if the precision `p` is incorrect.
    pub fn log(&self, n: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match &self.inner {
            Flavor::Value(v1) => {
                match &n.inner {
                    Flavor::Value(v2) => {
                        if v2.is_zero() {
                            return INF_NEG;
                        }
                        Self::result_to_ext(v1.log(v2, p, rm, cc), false, true)
                    }
                    Flavor::Inf(s2) => {
                        // v1.log(inf)
                        if s2.is_positive() {
                            Self::new(p)
                        } else {
                            NAN
                        }
                    }
                    Flavor::NaN(err) => Self::nan(*err),
                }
            }
            Flavor::Inf(s1) => {
                if *s1 == Sign::Neg {
                    // -inf.log(any)
                    NAN
                } else {
                    match &n.inner {
                        Flavor::Value(v2) => {
                            // +inf.log(v2)
                            if v2.exponent() <= 0 {
                                INF_NEG
                            } else {
                                INF_POS
                            }
                        }
                        Flavor::Inf(_) => NAN, // +inf.log(inf)
                        Flavor::NaN(err) => Self::nan(*err),
                    }
                }
            }
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    /// Returns true if `self` is positive.
    /// The function returns false if `self` is NaN.
    pub fn is_positive(&self) -> bool {
        match &self.inner {
            Flavor::Value(v) => v.is_positive(),
            Flavor::Inf(s) => *s == Sign::Pos,
            Flavor::NaN(_) => false,
        }
    }

    /// Returns true if `self` is negative.
    /// The function returns false if `self` is NaN.
    pub fn is_negative(&self) -> bool {
        match &self.inner {
            Flavor::Value(v) => v.is_negative(),
            Flavor::Inf(s) => *s == Sign::Neg,
            Flavor::NaN(_) => false,
        }
    }

    /// Returns true if `self` is subnormal. A number is subnormal if the most significant bit of the mantissa is not equal to 1.
    pub fn is_subnormal(&self) -> bool {
        if let Flavor::Value(v) = &self.inner {
            return v.is_subnormal();
        }
        false
    }

    /// Returns true if `self` is zero.
    pub fn is_zero(&self) -> bool {
        match &self.inner {
            Flavor::Value(v) => v.is_zero(),
            Flavor::Inf(_) => false,
            Flavor::NaN(_) => false,
        }
    }

    /// Restricts the value of `self` to an interval determined by the values of `min` and `max`.
    /// The function returns `max` if `self` is greater than `max`, `min` if `self` is less than `min`, and `self` otherwise.
    /// If either argument is NaN or `min` is greater than `max`, the function returns NaN.
    pub fn clamp(&self, min: &Self, max: &Self) -> Self {
        if self.is_nan() || min.is_nan() || max.is_nan() || max.cmp(min).unwrap() < 0 {
            // call to unwrap() is unreacheable
            NAN
        } else if self.cmp(min).unwrap() < 0 {
            // call to unwrap() is unreacheable
            min.clone()
        } else if self.cmp(max).unwrap() > 0 {
            // call to unwrap() is unreacheable
            max.clone()
        } else {
            self.clone()
        }
    }

    /// Returns the value of `d1` if `d1` is greater than `self`, or the value of `self` otherwise.
    /// If either argument is NaN, the function returns NaN.
    pub fn max(&self, d1: &Self) -> Self {
        if self.is_nan() || d1.is_nan() {
            NAN
        } else if self.cmp(d1).unwrap() < 0 {
            // call to unwrap() is unreacheable
            d1.clone()
        } else {
            self.clone()
        }
    }

    /// Returns value of `d1` if `d1` is less than `self`, or the value of `self` otherwise.
    /// If either argument is NaN, the function returns NaN.
    pub fn min(&self, d1: &Self) -> Self {
        if self.is_nan() || d1.is_nan() {
            NAN
        } else if self.cmp(d1).unwrap() > 0 {
            // call to unwrap() is unreacheable
            d1.clone()
        } else {
            self.clone()
        }
    }

    /// Returns a ExactNum with the value -1 if `self` is negative, 1 if `self` is positive, zero otherwise.
    /// The function returns NaN If `self` is NaN.
    pub fn signum(&self) -> Self {
        if self.is_nan() {
            NAN
        } else if self.is_negative() {
            let mut ret = Self::from_u8(1, DEFAULT_P);
            ret.inv_sign();
            ret
        } else {
            Self::from_u8(1, DEFAULT_P)
        }
    }

    /// Parses a number from the string `s`.
    /// The function expects `s` to be a number in scientific format in radix `rdx`, or +-Inf, or NaN.
    /// if `p` equals to usize::MAX then the precision of the resulting number is determined automatically from the input.
    ///
    /// ## Examples
    ///
    /// ```
    /// # use zenith_float_num::ExactNum;
    /// # use zenith_float_num::Radix;
    /// # use zenith_float_num::RoundingMode;
    /// # use zenith_float_num::Consts;
    /// let mut cc = Consts::new().expect("Constants cache initialized.");
    ///
    /// let n = ExactNum::parse("0.0", Radix::Bin, 64, RoundingMode::ToEven, &mut cc);
    /// assert!(n.is_zero());
    ///
    /// let n = ExactNum::parse("-Inf", Radix::Hex, 1, RoundingMode::None, &mut cc);
    /// assert!(n.is_inf_neg());
    ///
    /// let n = ExactNum::parse("NaN", Radix::Oct, 2, RoundingMode::None, &mut cc);
    /// assert!(n.is_nan());
    /// ```
    pub fn parse(s: &str, rdx: Radix, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match crate::parser::parse(s, rdx) {
            Ok(ps) => {
                if ps.is_inf() {
                    if ps.sign() == Sign::Pos {
                        INF_POS
                    } else {
                        INF_NEG
                    }
                } else if ps.is_nan() {
                    NAN
                } else {
                    let (m, s, e) = ps.raw_parts();
                    Self::result_to_ext(
                        ExactNumNumber::convert_from_radix(s, m, e, rdx, p, rm, cc),
                        false,
                        true,
                    )
                }
            }
            Err(e) => Self::nan(Some(e)),
        }
    }

    #[cfg(feature = "std")]
    pub(crate) fn write_str<T: Write>(
        &self,
        w: &mut T,
        rdx: Radix,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<(), core::fmt::Error> {
        match &self.inner {
            Flavor::Value(v) => match v.format(rdx, rm, cc) {
                Ok(s) => w.write_str(&s),
                Err(e) => match e {
                    Error::ExponentOverflow(s) => {
                        if s.is_positive() {
                            w.write_str("Inf")
                        } else {
                            w.write_str("-Inf")
                        }
                    }
                    _ => w.write_str("Err"),
                },
            },
            Flavor::Inf(sign) => {
                let s = if sign.is_negative() { "-Inf" } else { "Inf" };
                w.write_str(s)
            }
            crate::ext::Flavor::NaN(_) => w.write_str("NaN"),
        }
    }

    /// Formats the number using radix `rdx` and rounding mode `rm`.
    /// Note, since hexadecimal digits include the character "e", the exponent part is separated
    /// from the mantissa by "_".
    /// For example, a number with mantissa `123abcdef` and exponent `123` would be formatted as `123abcdef_e+123`.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - ExponentOverflow: the resulting exponent becomes greater than the maximum allowed value for the exponent.
    pub fn format(&self, rdx: Radix, rm: RoundingMode, cc: &mut Consts) -> Result<String, Error> {
        let s = match &self.inner {
            Flavor::Value(v) => match v.format(rdx, rm, cc) {
                Ok(s) => return Ok(s),
                Err(e) => match e {
                    Error::ExponentOverflow(s) => {
                        if s.is_positive() {
                            "Inf"
                        } else {
                            "-Inf"
                        }
                    }
                    _ => "Err",
                },
            },
            Flavor::Inf(sign) => {
                if sign.is_negative() {
                    "-Inf"
                } else {
                    "Inf"
                }
            }
            crate::ext::Flavor::NaN(_) => "NaN",
        };

        let mut ret = String::new();
        ret.try_reserve_exact(s.len())?;
        ret.push_str(s);

        Ok(ret)
    }

    /// Wraps `self` in a [`crate::RadixFloat`] tagged with `radix` for parse/format.
    pub fn with_radix(self, radix: Radix) -> crate::radix_float::RadixFloat {
        crate::radix_float::RadixFloat::with_radix(self, radix)
    }

    /// Returns a random normalized (not subnormal) ExactNum number with exponent in the range
    /// from `exp_from` to `exp_to` inclusive. The sign can be positive and negative. Zero is excluded.
    /// Precision is rounded upwards to the word size.
    /// Function does not follow any specific distribution law.
    /// The intended use of this function is for testing.
    /// The function returns NaN if the precision `p` is incorrect or when `exp_from` is less than EXPONENT_MIN or `exp_to` is greater than EXPONENT_MAX.
    #[cfg(feature = "random")]
    pub fn random_normal(p: usize, exp_from: Exponent, exp_to: Exponent) -> Self {
        Self::result_to_ext(
            ExactNumNumber::random_normal(p, exp_from, exp_to),
            false,
            true,
        )
    }

    /// Returns category of `self`.
    pub fn classify(&self) -> FpCategory {
        match &self.inner {
            Flavor::Value(v) => {
                if v.is_subnormal() {
                    FpCategory::Subnormal
                } else if v.is_zero() {
                    FpCategory::Zero
                } else {
                    FpCategory::Normal
                }
            }
            Flavor::Inf(_) => FpCategory::Infinite,
            Flavor::NaN(_) => FpCategory::Nan,
        }
    }

    /// Computes the arctangent of a number with precision `p`. The result is rounded using the rounding mode `rm`.
    /// This function requires constants cache `cc` for computing the result.
    /// Precision is rounded upwards to the word size.
    /// The function returns NaN if the precision `p` is incorrect.
    pub fn atan(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match &self.inner {
            Flavor::Value(v) => Self::result_to_ext(v.atan(p, rm, cc), v.is_zero(), true),
            Flavor::Inf(s) => Self::result_to_ext(Self::half_pi(*s, p, rm, cc), false, true),
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    /// Computes `atan2(self, x)` with precision `p` (quadrant-aware arctangent of `self / x`).
    /// The result is rounded using the rounding mode `rm`.
    /// This function requires constants cache `cc`.
    /// Precision is rounded upwards to the word size.
    pub fn atan2(&self, x: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return self.clone();
        }
        if x.is_nan() {
            return x.clone();
        }

        match (&self.inner, &x.inner) {
            (Flavor::Inf(sy), Flavor::Inf(sx)) => {
                let mut q = cc.pi(p, rm);
                q = q.div(&ExactNum::from_word(4, p), p, rm);
                if sx.is_negative() {
                    let three = ExactNum::from_word(3, p);
                    q = three.mul(&q, p, rm);
                }
                if sy.is_negative() {
                    q.neg()
                } else {
                    q
                }
            }
            (Flavor::Inf(sy), Flavor::Value(_)) => {
                Self::result_to_ext(Self::half_pi(*sy, p, rm, cc), false, true)
            }
            (Flavor::Value(y), Flavor::Inf(sx)) => {
                if sx.is_positive() {
                    Self::result_to_ext(ExactNumNumber::new2(p, y.sign(), y.inexact()), false, true)
                } else {
                    let mut pi = cc.pi(p, rm);
                    pi.set_sign(y.sign());
                    pi
                }
            }
            (Flavor::Value(y), Flavor::Value(xv)) => {
                Self::result_to_ext(y.atan2(xv, p, rm, cc), false, false)
            }
            _ => NAN,
        }
    }

    /// Computes `sqrt(self² + other²)` with precision `p`.
    /// The result is rounded using the rounding mode `rm`.
    /// Precision is rounded upwards to the word size.
    /// `hypot(±Inf, y)` and `hypot(x, ±Inf)` are `+Inf`, including when the other argument is NaN.
    pub fn hypot(&self, other: &Self, p: usize, rm: RoundingMode) -> Self {
        if self.is_inf() || other.is_inf() {
            return INF_POS;
        }
        if self.is_nan() {
            return self.clone();
        }
        if other.is_nan() {
            return other.clone();
        }
        match (&self.inner, &other.inner) {
            (Flavor::Value(a), Flavor::Value(b)) => {
                Self::result_to_ext(a.hypot(b, p, rm), false, true)
            }
            _ => NAN,
        }
    }

    /// Computes `ln(1 + self)` with precision `p`.
    /// The result is rounded using the rounding mode `rm`.
    /// This function requires constants cache `cc`.
    /// Returns `-Inf` for `self == -1`, and NaN if `self < -1`.
    pub fn log1p(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match &self.inner {
            Flavor::Value(v) => Self::result_to_ext(v.log1p(p, rm, cc), false, false),
            Flavor::Inf(s) => {
                if s.is_positive() {
                    INF_POS
                } else {
                    NAN
                }
            }
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    /// Computes `exp(self) - 1` with precision `p`.
    /// The result is rounded using the rounding mode `rm`.
    /// This function requires constants cache `cc`.
    pub fn expm1(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match &self.inner {
            Flavor::Value(v) => Self::result_to_ext(v.expm1(p, rm, cc), false, true),
            Flavor::Inf(s) => {
                if s.is_positive() {
                    INF_POS
                } else {
                    ExactNum::from_i8(-1, p)
                }
            }
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    /// Computes the hyperbolic tangent of a number with precision `p`. The result is rounded using the rounding mode `rm`.
    /// This function requires constants cache `cc` for computing the result.
    /// Precision is rounded upwards to the word size.
    /// The function returns NaN if the precision `p` is incorrect.
    pub fn tanh(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match &self.inner {
            Flavor::Value(v) => Self::result_to_ext(v.tanh(p, rm, cc), v.is_zero(), true),
            Flavor::Inf(s) => Self::from_i8(s.to_int(), p),
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    fn half_pi(
        s: Sign,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<ExactNumNumber, Error> {
        let mut half_pi = cc.pi_num(p, rm)?;

        half_pi.set_exponent(1);
        half_pi.set_sign(s);

        Ok(half_pi)
    }

    fn result_to_ext(
        res: Result<ExactNumNumber, Error>,
        is_dividend_zero: bool,
        is_same_sign: bool,
    ) -> ExactNum {
        match res {
            Err(e) => match e {
                Error::ExponentOverflow(s) => {
                    if s.is_positive() {
                        INF_POS
                    } else {
                        INF_NEG
                    }
                }
                Error::DivisionByZero => {
                    if is_dividend_zero {
                        NAN
                    } else if is_same_sign {
                        INF_POS
                    } else {
                        INF_NEG
                    }
                }
                Error::MemoryAllocation => Self::nan(Some(Error::MemoryAllocation)),
                Error::InvalidArgument => Self::nan(Some(Error::InvalidArgument)),
                Error::PrecisionRetryExhausted => Self::nan(Some(Error::PrecisionRetryExhausted)),
            },
            Ok(v) => ExactNum {
                inner: Flavor::Value(v),
            },
        }
    }

    /// Returns the exponent of `self`, or None if `self` is Inf or NaN.
    pub fn exponent(&self) -> Option<Exponent> {
        match &self.inner {
            Flavor::Value(v) => Some(v.exponent()),
            _ => None,
        }
    }

    /// Returns the number of significant bits used in the mantissa, or None if `self` is Inf or NaN.
    /// Normal numbers use all bits of the mantissa.
    /// Subnormal numbers use fewer bits than the mantissa can hold.
    pub fn precision(&self) -> Option<usize> {
        match &self.inner {
            Flavor::Value(v) => Some(v.precision()),
            _ => None,
        }
    }

    /// Returns the maximum value for the specified precision `p`: all bits of the mantissa are set to 1,
    /// the exponent has the maximum possible value, and the sign is positive.
    /// Precision is rounded upwards to the word size.
    /// The function returns NaN if the precision `p` is incorrect.
    pub fn max_value(p: usize) -> Self {
        Self::result_to_ext(ExactNumNumber::max_value(p), false, true)
    }

    /// Returns the minimum value for the specified precision `p`: all bits of the mantissa are set to 1, the exponent has the maximum possible value, and the sign is negative. Precision is rounded upwards to the word size.
    /// The function returns NaN if the precision `p` is incorrect.
    pub fn min_value(p: usize) -> Self {
        Self::result_to_ext(ExactNumNumber::min_value(p), false, true)
    }

    /// Returns the minimum positive subnormal value for the specified precision `p`:
    /// only the least significant bit of the mantissa is set to 1, the exponent has
    /// the minimum possible value, and the sign is positive.
    /// Precision is rounded upwards to the word size.
    /// The function returns NaN if the precision `p` is incorrect.
    pub fn min_positive(p: usize) -> Self {
        Self::result_to_ext(ExactNumNumber::min_positive(p), false, true)
    }

    /// Returns the minimum positive normal value for the specified precision `p`:
    /// only the most significant bit of the mantissa is set to 1, the exponent has
    /// the minimum possible value, and the sign is positive.
    /// Precision is rounded upwards to the word size.
    /// The function returns NaN if the precision `p` is incorrect.
    pub fn min_positive_normal(p: usize) -> Self {
        Self::result_to_ext(ExactNumNumber::min_positive_normal(p), false, true)
    }

    /// Returns a new number with value `d` and the precision `p`. Precision is rounded upwards to the word size.
    /// The function returns NaN if the precision `p` is incorrect.
    pub fn from_word(d: Word, p: usize) -> Self {
        Self::result_to_ext(ExactNumNumber::from_word(d, p), false, true)
    }

    /// Returns a copy of the number with the sign reversed.
    pub fn neg(&self) -> Self {
        let mut ret = self.clone();
        ret.inv_sign();
        ret
    }

    /// Decomposes `self` into raw parts.
    /// The function returns a reference to a slice of words representing mantissa,
    /// numbers of significant bits in the mantissa, sign, exponent,
    /// and a bool value which specify whether the number is inexact.
    pub fn as_raw_parts(&self) -> Option<(&[Word], usize, Sign, Exponent, bool)> {
        if let Flavor::Value(v) = &self.inner {
            Some(v.as_raw_parts())
        } else {
            None
        }
    }

    /// Constructs a number from the raw parts:
    ///
    ///  - `m` is the mantisaa.
    ///  - `n` is the number of significant bits in mantissa.
    ///  - `s` is the sign.
    ///  - `e` is the exponent.
    ///  - `inexact` specify whether number is inexact.
    ///
    /// This function returns NaN in the following situations:
    ///
    /// - `n` is larger than the number of bits in `m`.
    /// - `n` is smaller than the number of bits in `m`, but `m` does not represent corresponding subnormal number mantissa.
    /// - `n` is smaller than the number of bits in `m`, but `e` is not the minimum possible exponent.
    /// - `n` or the size of `m` is too large (larger than isize::MAX / 2 + EXPONENT_MIN).
    /// - `e` is less than EXPONENT_MIN or greater than EXPONENT_MAX.
    pub fn from_raw_parts(m: &[Word], n: usize, s: Sign, e: Exponent, inexact: bool) -> Self {
        Self::result_to_ext(
            crate::mantissa::Mantissa::from_raw_parts(m, n)
                .map(|mantissa| ExactNumNumber::from_raw_unchecked(mantissa, s, e, inexact)),
            false,
            true,
        )
    }

    /// Constructs a number from the slice of words:
    ///
    ///  - `m` is the mantissa.
    ///  - `s` is the sign.
    ///  - `e` is the exponent.
    ///
    /// The function returns NaN if `e` is less than EXPONENT_MIN or greater than EXPONENT_MAX.
    pub fn from_words(m: &[Word], s: Sign, e: Exponent) -> Self {
        Self::result_to_ext(ExactNumNumber::from_words(m, s, e), false, true)
    }

    /// Returns the sign of `self`, or None if `self` is NaN.
    pub fn sign(&self) -> Option<Sign> {
        match &self.inner {
            Flavor::Value(v) => Some(v.sign()),
            Flavor::Inf(s) => Some(*s),
            Flavor::NaN(_) => None,
        }
    }

    /// Sets the exponent of `self`.
    /// Note that if `self` is subnormal, the exponent may not change, but the mantissa will shift instead.
    /// `e` will be clamped to the range from EXPONENT_MIN to EXPONENT_MAX if it's outside of the range.
    /// See example below.
    ///
    /// ## Examples
    ///
    /// ```
    /// # use zenith_float_num::ExactNum;
    /// # use zenith_float_num::EXPONENT_MIN;
    /// // construct a subnormal value.
    /// let mut n = ExactNum::min_positive(128);
    ///
    /// assert_eq!(n.exponent(), Some(EXPONENT_MIN));
    /// assert_eq!(n.precision(), Some(1));
    ///
    /// // increase exponent.
    /// let n_exp = n.exponent().expect("n is not NaN");
    /// n.set_exponent(n_exp + 1);
    ///
    /// // the outcome for subnormal number.
    /// assert_eq!(n.exponent(), Some(EXPONENT_MIN));
    /// assert_eq!(n.precision(), Some(2));
    /// ```
    pub fn set_exponent(&mut self, e: Exponent) {
        if let Flavor::Value(v) = &mut self.inner {
            v.set_exponent(e)
        }
    }

    /// Returns the maximum mantissa length of `self` in bits regardless of whether `self` is normal or subnormal.
    pub fn mantissa_max_bit_len(&self) -> Option<usize> {
        if let Flavor::Value(v) = &self.inner {
            Some(v.mantissa_max_bit_len())
        } else {
            None
        }
    }

    /// True when a finite value stores its mantissa on the stack (at most [`crate::INLINE_WORDS`] limbs).
    /// Inf and NaN return `false`.
    pub fn is_inline(&self) -> bool {
        match &self.inner {
            Flavor::Value(v) => v.is_inline(),
            Flavor::Inf(_) | Flavor::NaN(_) => false,
        }
    }

    /// Sets the precision of `self` to `p`.
    /// If the new precision is smaller than the existing one, the number is rounded using specified rounding mode `rm`.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - InvalidArgument: the precision is incorrect.
    pub fn set_precision(&mut self, p: usize, rm: RoundingMode) -> Result<(), Error> {
        if let Flavor::Value(v) = &mut self.inner {
            v.set_precision(p, rm)
        } else {
            Ok(())
        }
    }

    /// Computes the reciprocal of a number with precision `p`.
    /// The result is rounded using the rounding mode `rm`.
    /// Precision is rounded upwards to the word size.
    /// The function returns NaN if the precision `p` is incorrect.
    pub fn reciprocal(&self, p: usize, rm: RoundingMode) -> Self {
        match &self.inner {
            Flavor::Value(v) => Self::result_to_ext(v.reciprocal(p, rm), false, v.is_positive()),
            Flavor::Inf(s) => {
                let mut ret = Self::new(p);
                ret.set_sign(*s);
                ret
            }
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    /// Sets the sign of `self`.
    pub fn set_sign(&mut self, s: Sign) {
        match &mut self.inner {
            Flavor::Value(v) => v.set_sign(s),
            Flavor::Inf(_) => self.inner = Flavor::Inf(s),
            Flavor::NaN(_) => {}
        };
    }

    /// Returns the raw mantissa words of a number.
    pub fn mantissa_digits(&self) -> Option<&[Word]> {
        if let Flavor::Value(v) = &self.inner {
            Some(v.mantissa().digits())
        } else {
            None
        }
    }

    /// Converts an array of digits in radix `rdx` to ExactNum with precision `p`.
    /// `digits` represents mantissa and is interpreted as a number smaller than 1 and greater or equal to 1/`rdx`.
    /// The first element in `digits` is the most significant digit.
    /// `e` is the exponent part of the number, such that the number can be represented as `digits` * `rdx` ^ `e`.
    /// Precision is rounded upwards to the word size.
    /// if `p` equals usize::MAX then the precision of the resulting number is determined automatically from the input.
    ///
    /// ## Examples
    ///
    /// Code below converts `-0.1234567₈ × 10₈^3₈` given in radix 8 to ExactNum.
    ///
    /// ``` rust
    /// # use zenith_float_num::{ExactNum, Sign, RoundingMode, Radix, Consts};
    /// let mut cc = Consts::new().expect("Constants cache initialized.");
    ///
    /// let n = ExactNum::convert_from_radix(
    ///     Sign::Neg,
    ///     &[1, 2, 3, 4, 5, 6, 7, 0],
    ///     3,
    ///     Radix::Oct,
    ///     64,
    ///     RoundingMode::None,
    ///     &mut cc);
    /// assert!(!n.is_nan());
    /// assert!(n.is_negative());
    /// ```
    ///
    /// ## Errors
    ///
    /// On error, the function returns NaN with the following associated error:
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - ExponentOverflow: the resulting exponent becomes greater than the maximum allowed value for the exponent.
    ///  - InvalidArgument: the precision is incorrect, or `digits` contains unacceptable digits for given radix,
    ///    or when `e` is less than EXPONENT_MIN or greater than EXPONENT_MAX.
    pub fn convert_from_radix(
        sign: Sign,
        digits: &[u8],
        e: Exponent,
        rdx: Radix,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        Self::result_to_ext(
            ExactNumNumber::convert_from_radix(sign, digits, e, rdx, p, rm, cc),
            false,
            true,
        )
    }

    /// Converts `self` to radix `rdx` using rounding mode `rm`.
    /// The function returns sign, mantissa digits in radix `rdx`, and exponent such that the converted number
    /// can be represented as `mantissa digits` * `rdx` ^ `exponent`.
    /// The first element in the mantissa is the most significant digit.
    ///
    /// ## Examples
    ///
    /// ``` rust
    /// # use zenith_float_num::{ExactNum, Sign, RoundingMode, Radix, Consts};
    ///
    /// let mut cc = Consts::new().expect("Constants cache initialized.");
    /// let n = ExactNum::parse("123.45678", Radix::Dec, 64, RoundingMode::None, &mut cc);
    /// let (s, m, _e) = n.convert_to_radix(Radix::Dec, RoundingMode::None, &mut cc).expect("Conversion failed");
    /// assert_eq!(s, Sign::Pos);
    /// assert!(!m.is_empty());
    /// ```
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - ExponentOverflow: the resulting exponent becomes greater than the maximum allowed value for the exponent.
    ///  - InvalidArgument: `self` is Inf or NaN.
    pub fn convert_to_radix(
        &self,
        rdx: Radix,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<(Sign, Vec<u8>, Exponent), Error> {
        match &self.inner {
            Flavor::Value(v) => v.convert_to_radix(rdx, rm, cc),
            Flavor::NaN(_) => Err(Error::InvalidArgument),
            Flavor::Inf(_) => Err(Error::InvalidArgument),
        }
    }

    /// Returns true if `self` is inexact. The function returns false if `self` is Inf or NaN.
    pub fn inexact(&self) -> bool {
        if let Flavor::Value(v) = &self.inner {
            v.inexact()
        } else {
            false
        }
    }

    /// Marks `self` as inexact if `inexact` is true, or exact otherwise.
    /// The function has no effect if `self` is Inf or NaN.
    pub fn set_inexact(&mut self, inexact: bool) {
        if let Flavor::Value(v) = &mut self.inner {
            v.set_inexact(inexact);
        }
    }

    /// Try to round and then set the precision to `p`, given `self` has `s` correct digits in mantissa.
    /// The function returns true if rounding succeeded, or if `self` is Inf or NaN.
    /// If the fuction returns `false`, `self` is still modified, and should be discarded.
    /// In case of an error, `self` will be set to NaN with an associated error.
    /// If the precision `p` is incorrect `self` will be set to NaN.
    pub fn try_set_precision(&mut self, p: usize, rm: RoundingMode, s: usize) -> bool {
        if let Flavor::Value(v) = &mut self.inner {
            v.try_set_precision(p, rm, s).unwrap_or_else(|e| {
                self.inner = Flavor::NaN(Some(e));
                true
            })
        } else {
            true
        }
    }

    /// Split `self = m · 2^e` with `m` in `[0.5, 1)` (zeros return `(0, 0)`; Inf/NaN return `(self, 0)`).
    pub fn frexp(&self) -> (Self, Exponent) {
        match &self.inner {
            Flavor::Value(v) => match v.frexp() {
                Ok((m, e)) => (m.into(), e),
                Err(err) => (Self::nan(Some(err)), 0),
            },
            Flavor::Inf(_) | Flavor::NaN(_) => (self.clone(), 0),
        }
    }

    /// `self · 2^n`. Alias of [`Self::scalb`].
    pub fn ldexp(&self, n: Exponent, p: usize, rm: RoundingMode) -> Self {
        match &self.inner {
            Flavor::Value(v) => Self::result_to_ext(v.ldexp(n, p, rm), v.is_zero(), true),
            Flavor::Inf(_) | Flavor::NaN(_) => self.clone(),
        }
    }

    /// `self · 2^n` (IEEE `scalbn`).
    pub fn scalb(&self, n: Exponent, p: usize, rm: RoundingMode) -> Self {
        self.ldexp(n, p, rm)
    }

    /// `floor(log2(|self|))` as a float. Zero becomes `-Inf`; Inf/NaN unchanged in kind.
    pub fn logb(&self, p: usize, rm: RoundingMode) -> Self {
        match &self.inner {
            Flavor::Value(v) => Self::result_to_ext(v.logb(p, rm), v.is_zero(), true),
            Flavor::Inf(_) => INF_POS,
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    /// Integer `floor(log2(|self|))`. `None` for zero, Inf, or NaN.
    pub fn ilogb(&self) -> Option<Exponent> {
        match &self.inner {
            Flavor::Value(v) => v.ilogb().ok(),
            _ => None,
        }
    }
}

impl Clone for ExactNum {
    fn clone(&self) -> Self {
        match &self.inner {
            Flavor::Value(v) => Self::result_to_ext(v.clone(), false, true),
            Flavor::Inf(s) => {
                if s.is_positive() {
                    INF_POS
                } else {
                    INF_NEG
                }
            }
            Flavor::NaN(err) => Self::nan(*err),
        }
    }
}

macro_rules! gen_wrapper_arg {
    // function requires self as argument
    ($comment:literal, $fname:ident, $ret:ty, $pos_inf:block, $neg_inf:block, $($arg:ident, $arg_type:ty),*) => {
        #[doc=$comment]
        pub fn $fname(&self$(,$arg: $arg_type)*) -> $ret {
            match &self.inner {
                Flavor::Value(v) => Self::result_to_ext(v.$fname($($arg,)*), v.is_zero(), true),
                Flavor::Inf(s) => if s.is_positive() $pos_inf else $neg_inf,
                Flavor::NaN(err) => Self::nan(*err),
            }
        }
    };
}

macro_rules! gen_wrapper_arg_rm {
    // unwrap error, function requires self as argument
    ($comment:literal, $fname:ident, $ret:ty, $pos_inf:block, $neg_inf:block, $($arg:ident, $arg_type:ty),*) => {
        #[doc=$comment]
        pub fn $fname(&self$(,$arg: $arg_type)*, rm: RoundingMode) -> $ret {
            match &self.inner {
                Flavor::Value(v) => {
                    Self::result_to_ext(v.$fname($($arg,)* rm), v.is_zero(), true)
                },
                Flavor::Inf(s) => if s.is_positive() $pos_inf else $neg_inf,
                Flavor::NaN(err) => Self::nan(*err),
            }
        }
    };
}

macro_rules! gen_wrapper_arg_rm_cc {
    // unwrap error, function requires self as argument
    ($comment:literal, $fname:ident, $ret:ty, $pos_inf:block, $neg_inf:block, $($arg:ident, $arg_type:ty),*) => {
        #[doc=$comment]
        pub fn $fname(&self$(,$arg: $arg_type)*, rm: RoundingMode, cc: &mut Consts) -> $ret {
            match &self.inner {
                Flavor::Value(v) => {
                    Self::result_to_ext(v.$fname($($arg,)* rm, cc), v.is_zero(), true)
                },
                Flavor::Inf(s) => if s.is_positive() $pos_inf else $neg_inf,
                Flavor::NaN(err) => Self::nan(*err),
            }
        }
    };
}

macro_rules! gen_wrapper_log {
    ($comment:literal, $fname:ident, $ret:ty, $pos_inf:block, $neg_inf:block, $($arg:ident, $arg_type:ty),*) => {
        #[doc=$comment]
        pub fn $fname(&self$(,$arg: $arg_type)*, rm: RoundingMode, cc: &mut Consts) -> $ret {
            match &self.inner {
                Flavor::Value(v) => {
                    if v.is_zero() {
                        return INF_NEG;
                    }
                    Self::result_to_ext(v.$fname($($arg,)* rm, cc), v.is_zero(), true)
                },
                Flavor::Inf(s) => if s.is_positive() $pos_inf else $neg_inf,
                Flavor::NaN(err) => Self::nan(*err),
            }
        }
    };
}

impl ExactNum {
    gen_wrapper_arg!(
        "Returns the absolute value of `self`.",
        abs,
        Self,
        { INF_POS },
        { INF_POS },
    );
    /// Returns a value with the magnitude of `self` and the sign of `sign`.
    pub fn copysign(&self, sign: &Self, p: usize, rm: RoundingMode) -> Self {
        if self.is_nan() {
            return self.clone();
        }
        let sign_num = match &sign.inner {
            Flavor::Value(v) => v.clone(),
            Flavor::Inf(s) => ExactNumNumber::from_i8(s.to_int(), p),
            Flavor::NaN(_) => ExactNumNumber::new(p),
        };
        let sign_num = match sign_num {
            Ok(v) => v,
            Err(e) => return Self::nan(Some(e)),
        };
        match &self.inner {
            Flavor::Value(v) => Self::result_to_ext(v.copysign(&sign_num, p, rm), false, true),
            Flavor::Inf(_) => {
                if sign.is_negative() || (sign.is_zero() && sign_num.is_negative()) {
                    INF_NEG
                } else {
                    INF_POS
                }
            }
            Flavor::NaN(err) => Self::nan(*err),
        }
    }
    /// Returns the next representable value from `self` toward `toward` at precision `p`.
    pub fn next_after(&self, toward: &Self, p: usize, rm: RoundingMode) -> Self {
        if self.is_nan() {
            return self.clone();
        }
        if toward.is_nan() {
            return toward.clone();
        }
        match (&self.inner, &toward.inner) {
            (Flavor::Value(v), Flavor::Value(t)) => {
                Self::result_to_ext(v.next_after(t, p, rm), false, true)
            }
            (Flavor::Inf(_), _) | (_, Flavor::Inf(_)) => self.clone(),
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
        }
    }
    gen_wrapper_arg!("Returns the integer part of `self`.", int, Self, { NAN }, {
        NAN
    },);
    gen_wrapper_arg!(
        "Returns the fractional part of `self`.",
        fract,
        Self,
        { NAN },
        { NAN },
    );
    gen_wrapper_arg!(
        "Returns the smallest integer greater than or equal to `self`.",
        ceil,
        Self,
        { INF_POS },
        { INF_NEG },
    );
    gen_wrapper_arg!(
        "Returns the largest integer less than or equal to `self`.",
        floor,
        Self,
        { INF_POS },
        { INF_NEG },
    );
    gen_wrapper_arg_rm!("Returns the rounded number with `n` binary positions in the fractional part of the number using rounding mode `rm`.", 
        round,
        Self,
        { INF_POS },
        { INF_NEG },
        n,
        usize
    );
    gen_wrapper_arg_rm!(
        "Computes the square root of a number with precision `p`. The result is rounded using the rounding mode `rm`.
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.",
        sqrt,
        Self,
        { INF_POS },
        { NAN },
        p,
        usize
    );
    gen_wrapper_arg_rm!(
        "Computes the cube root of a number with precision `p`. The result is rounded using the rounding mode `rm`.
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.",
        cbrt,
        Self,
        { INF_POS },
        { INF_NEG },
        p,
        usize
    );
    /// Computes the `n`-th root of `self` with precision `p`. `n = 2` and `n = 3` delegate to [`sqrt`](Self::sqrt) and [`cbrt`](Self::cbrt).
    pub fn nth_root(&self, n: usize, p: usize, rm: RoundingMode) -> Self {
        if n == 0 {
            return Self::nan(Some(Error::InvalidArgument));
        }
        match &self.inner {
            Flavor::Value(v) => Self::result_to_ext(v.nth_root(n, p, rm), v.is_zero(), true),
            Flavor::Inf(s) => {
                if n % 2 == 0 {
                    if s.is_negative() {
                        NAN
                    } else {
                        INF_POS
                    }
                } else if s.is_negative() {
                    INF_NEG
                } else {
                    INF_POS
                }
            }
            Flavor::NaN(err) => Self::nan(*err),
        }
    }
    gen_wrapper_log!(
        "Computes the natural logarithm of a number with precision `p`. The result is rounded using the rounding mode `rm`.
        This function requires constants cache `cc` for computing the result.
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.",
        ln,
        Self,
        { INF_POS },
        { NAN },
        p,
        usize
    );
    gen_wrapper_log!(
        "Computes the logarithm base 2 of a number with precision `p`. The result is rounded using the rounding mode `rm`.
        This function requires constants cache `cc` for computing the result.
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.",
        log2,
        Self,
        { INF_POS },
        { NAN },
        p,
        usize
    );
    gen_wrapper_log!(
        "Computes the logarithm base 10 of a number with precision `p`. The result is rounded using the rounding mode `rm`.
        This function requires constants cache `cc` for computing the result.
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.",
        log10,
        Self,
        { INF_POS },
        { NAN },
        p,
        usize
    );
    gen_wrapper_arg_rm_cc!(
        "Computes `e` to the power of `self` with precision `p`. The result is rounded using the rounding mode `rm`.
        This function requires constants cache `cc` for computing the result.
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.",
        exp,
        Self,
        { INF_POS },
        { Self::new(p) },
        p,
        usize
    );
    gen_wrapper_arg_rm_cc!(
        "Computes `2` to the power of `self` with precision `p`. The result is rounded using the rounding mode `rm`.
        This function requires constants cache `cc` for computing the result.
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.",
        exp2,
        Self,
        { INF_POS },
        { Self::new(p) },
        p,
        usize
    );
    gen_wrapper_arg_rm_cc!(
        "Computes `10` to the power of `self` with precision `p`. The result is rounded using the rounding mode `rm`.
        This function requires constants cache `cc` for computing the result.
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.",
        exp10,
        Self,
        { INF_POS },
        { Self::new(p) },
        p,
        usize
    );
    gen_wrapper_arg_rm_cc!(
        "Reduces `self` modulo `2π` into the interval `(-2π, 2π)` using precision `p` and rounding mode `rm`.
        This function requires constants cache `cc` for computing the result.",
        rem_pi,
        Self,
        { NAN },
        { NAN },
        p,
        usize
    );

    gen_wrapper_arg_rm_cc!(
        "Computes the sine of a number with precision `p`. The result is rounded using the rounding mode `rm`.
        This function requires constants cache `cc` for computing the result.
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.",
        sin,
        Self,
        { NAN },
        { NAN },
        p,
        usize
    );
    gen_wrapper_arg_rm_cc!(
        "Computes the cosine of a number with precision `p`. The result is rounded using the rounding mode `rm`.
        This function requires constants cache `cc` for computing the result.
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.",
        cos,
        Self,
        { NAN },
        { NAN },
        p,
        usize
    );
    /// Computes `(sin(self), cos(self))` with precision `p` using a shared argument reduction.
    pub fn sin_cos(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> (Self, Self) {
        match &self.inner {
            Flavor::Value(v) => match v.sin_cos(p, rm, cc) {
                Ok((s, c)) => (
                    Self::result_to_ext(Ok(s), false, true),
                    Self::result_to_ext(Ok(c), false, true),
                ),
                Err(e) => (Self::nan(Some(e)), Self::nan(Some(e))),
            },
            Flavor::Inf(_) => (NAN, NAN),
            Flavor::NaN(err) => (Self::nan(*err), Self::nan(*err)),
        }
    }
    gen_wrapper_arg_rm_cc!(
        "Computes the tangent of a number with precision `p`. The result is rounded using the rounding mode `rm`.
        This function requires constants cache `cc` for computing the result.
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.",
        tan,
        Self,
        { NAN },
        { NAN },
        p,
        usize
    );
    gen_wrapper_arg_rm_cc!(
        "Computes the arcsine of a number with precision `p`. The result is rounded using the rounding mode `rm`.
        This function requires constants cache `cc` for computing the result.
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.", 
        asin,
        Self,
        {NAN},
        {NAN},
        p,
        usize
    );
    gen_wrapper_arg_rm_cc!(
        "Computes the arccosine of a number with precision `p`. The result is rounded using the rounding mode `rm`.
        This function requires constants cache `cc` for computing the result.
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.",
        acos,
        Self,
        { NAN },
        { NAN },
        p,
        usize
    );

    gen_wrapper_arg_rm_cc!(
        "Computes the hyperbolic sine of a number with precision `p`. The result is rounded using the rounding mode `rm`.
        This function requires constants cache cc for computing the result. 
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.",
        sinh,
        Self,
        { INF_POS },
        { INF_NEG },
        p,
        usize
    );
    gen_wrapper_arg_rm_cc!(
        "Computes the hyperbolic cosine of a number with precision `p`. The result is rounded using the rounding mode `rm`.
        This function requires constants cache cc for computing the result. 
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.",
        cosh,
        Self,
        { INF_POS },
        { INF_POS },
        p,
        usize
    );
    /// Computes `(sinh(self), cosh(self))` with precision `p` using a single `exp(|x|)` evaluation.
    pub fn sinh_cosh(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> (Self, Self) {
        match &self.inner {
            Flavor::Value(v) => match v.sinh_cosh(p, rm, cc) {
                Ok((s, c)) => (
                    Self::result_to_ext(Ok(s), false, true),
                    Self::result_to_ext(Ok(c), false, true),
                ),
                Err(Error::ExponentOverflow(s)) => {
                    if s.is_positive() {
                        (INF_POS, INF_POS)
                    } else {
                        (INF_NEG, INF_POS)
                    }
                }
                Err(e) => (Self::nan(Some(e)), Self::nan(Some(e))),
            },
            Flavor::Inf(s) => {
                if s.is_positive() {
                    (INF_POS, INF_POS)
                } else {
                    (INF_NEG, INF_POS)
                }
            }
            Flavor::NaN(err) => (Self::nan(*err), Self::nan(*err)),
        }
    }
    gen_wrapper_arg_rm_cc!(
        "Error function `erf(self)` with precision `p`.

# Precision

- Algorithm: Taylor series when `|x|.exponent() ≤ 2`; complementary asymptotic otherwise. Saturates to `±1` when `2|e| > p+4`.
- Bound: Ziv correct-rounding (`MAX_PREC_RETRY`). 1 ULP vs MPFR on `|x| ≲ 4`.
- Thresholds: exponent cut `≤ 2` (not a named constant).
- MPFR oracle: yes, `|x| ≲ 4` under `mpfr-tests`.",
        erf,
        Self,
        { ExactNum::from_u8(1, p) },
        { ExactNum::from_i8(-1, p) },
        p,
        usize
    );
    gen_wrapper_arg_rm_cc!(
        "Complementary error function `erfc(self) = 1 - erf(self)` with precision `p`.

# Precision

- Algorithm: `1 - erf` at extra working precision (same series / asymptotic as `erf`).
- Bound: Ziv correct-rounding (`MAX_PREC_RETRY`). 1 ULP vs MPFR on `|x| ≲ 4`.
- MPFR oracle: yes, `|x| ≲ 4` under `mpfr-tests`.",
        erfc,
        Self,
        { Self::new(p) },
        { ExactNum::from_u8(2, p) },
        p,
        usize
    );
    gen_wrapper_arg_rm_cc!(
        "Gamma function `Γ(self)` with precision `p`. Poles at non-positive integers yield NaN (or +Inf at 0).

# Precision

- Algorithm: Stirling series for `ln Γ` then `exp`; reflection across the negative axis. Integer factorials for small positive integers.
- Bound: Ziv correct-rounding (`MAX_PREC_RETRY`). 1 ULP vs MPFR on the oracle domain.
- MPFR oracle: yes, under `mpfr-tests`.",
        gamma,
        Self,
        { INF_POS },
        { NAN },
        p,
        usize
    );
    gen_wrapper_arg_rm_cc!(
        "`ln Γ(self)` for positive `self` with precision `p`.

# Precision

- Algorithm: Stirling series (Bernoulli) at working precision `p + WORD_BIT_SIZE`.
- Bound: Ziv correct-rounding (`MAX_PREC_RETRY`). 1 ULP vs MPFR on the oracle domain.
- MPFR oracle: yes, under `mpfr-tests`.",
        ln_gamma,
        Self,
        { INF_POS },
        { NAN },
        p,
        usize
    );
    gen_wrapper_arg_rm_cc!(
        "Digamma `ψ(self)`. Poles at non-positive integers. Reflection for z < 0.

# Precision

- Algorithm: recurrence to a large argument, then Bernoulli series; reflection for `z < 0`.
- Bound: Ziv correct-rounding (`MAX_PREC_RETRY`). Accumath `ψ` identities are the gold; MPFR `digamma` is not an oracle here.
- MPFR oracle: no.",
        digamma,
        Self,
        { INF_POS },
        { NAN },
        p,
        usize
    );
    /// Lower incomplete gamma `γ(self, x)` for `self > 0`, `x ≥ 0`.
    ///
    /// # Precision
    ///
    /// - Algorithm: power series in `x` at working precision; `+∞` in `x` returns `Γ(self)`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no (Accumath `γ(s,x)` golds).
    pub fn gammainc(&self, x: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match (&self.inner, &x.inner) {
            (Flavor::Value(s), Flavor::Value(xv)) => {
                Self::result_to_ext(s.gammainc(xv, p, rm, cc), xv.is_zero(), true)
            }
            (Flavor::Value(s), Flavor::Inf(sx)) => {
                if sx.is_positive() {
                    Self::result_to_ext(s.gamma(p, rm, cc), false, true)
                } else {
                    NAN
                }
            }
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
            (Flavor::Inf(_), _) => NAN,
        }
    }
    /// Upper incomplete gamma `Γ(self, x)` for `self > 0`, `x ≥ 0`.
    ///
    /// # Precision
    ///
    /// - Algorithm: `Γ(self) - γ(self, x)` at working precision; `+∞` in `x` returns 0.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    pub fn gammainc_upper(&self, x: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match (&self.inner, &x.inner) {
            (Flavor::Value(s), Flavor::Value(xv)) => {
                Self::result_to_ext(s.gammainc_upper(xv, p, rm, cc), xv.is_zero(), true)
            }
            (Flavor::Value(_), Flavor::Inf(sx)) => {
                if sx.is_positive() {
                    Self::new(p)
                } else {
                    NAN
                }
            }
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
            (Flavor::Inf(_), _) => NAN,
        }
    }
    gen_wrapper_arg_rm_cc!(
        "Exponential integral `Ei(self)` (principal value for `self < 0`). `0` is a pole.

# Precision

- Algorithm: power series for moderate `|x|`; factorial asymptotic when `|x|` is large (`exponent() > 6` and `|x| ≳ 0.7 p`).
- Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
- MPFR oracle: open (Accumath `Ei` golds).",
        ei,
        Self,
        { INF_POS },
        { Self::new(p) },
        p,
        usize
    );
    /// Sine integral `Si(self)`. `+∞ → π/2`, `−∞ → −π/2`.
    ///
    /// # Precision
    ///
    /// - Algorithm: series, or auxiliary `f,g` asymptotic on the same cut as `Ei`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: open.
    pub fn si(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match &self.inner {
            Flavor::Value(v) => Self::result_to_ext(v.si(p, rm, cc), v.is_zero(), true),
            Flavor::Inf(s) => Self::result_to_ext(Self::half_pi(*s, p, rm, cc), false, true),
            Flavor::NaN(err) => Self::nan(*err),
        }
    }
    gen_wrapper_arg_rm_cc!(
        "Cosine integral `Ci(self)` for `self > 0`.

# Precision

- Algorithm: series, or auxiliary `f,g` asymptotic (same `|x|` cut as `Ei`). Near-zero is a pole.
- Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
- MPFR oracle: open.",
        ci,
        Self,
        { Self::new(p) },
        { NAN },
        p,
        usize
    );
    gen_wrapper_arg_rm_cc!(
        "Logarithmic integral `li(self) = Ei(ln self)` for `self > 0`, `self ≠ 1`.

# Precision

- Algorithm: `Ei(ln self)` at extra working precision (inherits `Ei` series / asymptotic).
- Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
- MPFR oracle: open.",
        li,
        Self,
        { INF_POS },
        { NAN },
        p,
        usize
    );
    /// Fresnel sine integral `S(self)`. `±∞ → ±1/2`.
    ///
    /// # Precision
    ///
    /// - Algorithm: series, or auxiliary `f,g` when `|x|.exponent() ≥ 8` (or `3 x² > p`).
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: open.
    pub fn fresnel_s(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.fresnel_sc_ext(true, p, rm, cc)
    }

    /// Fresnel cosine integral `C(self)`. `±∞ → ±1/2`.
    ///
    /// # Precision
    ///
    /// - Algorithm: same series / auxiliary `f,g` split as [`Self::fresnel_s`].
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: open.
    pub fn fresnel_c(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.fresnel_sc_ext(false, p, rm, cc)
    }

    fn fresnel_sc_ext(&self, sine: bool, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match &self.inner {
            Flavor::Value(v) => {
                let inner = if sine { v.fresnel_s(p, rm, cc) } else { v.fresnel_c(p, rm, cc) };
                Self::result_to_ext(inner, v.is_zero(), true)
            }
            Flavor::Inf(s) => {
                let mut half = ExactNum::from_u8(1, p);
                half.set_exponent(0);
                half.set_sign(*s);
                half
            }
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    /// Airy \(\mathrm{Ai}(\mathrm{self})\). \(+\infty\to 0\); \(-\infty\) has no limit → NaN.
    ///
    /// # Precision
    ///
    /// - Algorithm: power series for `|x| < AIRY_SERIES_THRESHOLD` (`8`); asymptotic otherwise.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no (identity / ODE golds).
    pub fn ai(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match &self.inner {
            Flavor::Value(v) => Self::result_to_ext(v.ai(p, rm, cc), v.is_zero(), true),
            Flavor::Inf(s) => {
                if s.is_positive() {
                    Self::new(p)
                } else {
                    NAN
                }
            }
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    /// Airy \(\mathrm{Bi}(\mathrm{self})\). \(+\infty\to+\infty\); \(-\infty\) has no limit → NaN.
    ///
    /// # Precision
    ///
    /// - Algorithm: same `AIRY_SERIES_THRESHOLD = 8` split as [`Self::ai`].
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    pub fn bi(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match &self.inner {
            Flavor::Value(v) => Self::result_to_ext(v.bi(p, rm, cc), v.is_zero(), true),
            Flavor::Inf(s) => {
                if s.is_positive() {
                    INF_POS
                } else {
                    NAN
                }
            }
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    /// \(\mathrm{Ai}'(\mathrm{self})\). \(+\infty\to 0\); \(-\infty\) → NaN.
    ///
    /// # Precision
    ///
    /// - Algorithm: differentiated series / asymptotic; `AIRY_SERIES_THRESHOLD = 8`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    pub fn ai_prime(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match &self.inner {
            Flavor::Value(v) => Self::result_to_ext(v.ai_prime(p, rm, cc), v.is_zero(), true),
            Flavor::Inf(s) => {
                if s.is_positive() {
                    Self::new(p)
                } else {
                    NAN
                }
            }
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    /// \(\mathrm{Bi}'(\mathrm{self})\). \(+\infty\to+\infty\); \(-\infty\) → NaN.
    ///
    /// # Precision
    ///
    /// - Algorithm: differentiated series / asymptotic; `AIRY_SERIES_THRESHOLD = 8`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    pub fn bi_prime(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match &self.inner {
            Flavor::Value(v) => Self::result_to_ext(v.bi_prime(p, rm, cc), v.is_zero(), true),
            Flavor::Inf(s) => {
                if s.is_positive() {
                    INF_POS
                } else {
                    NAN
                }
            }
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    /// Bessel function of the first kind `J_n(self)` for integer order `n`.
    ///
    /// # Precision
    ///
    /// - Algorithm: power series; Miller recurrence for large `n` (`n ≤ 1024`).
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`). 1 ULP vs MPFR `jn` for `n = 0,1,2`.
    /// - MPFR oracle: yes, `n = 0,1,2` under `mpfr-tests`.
    pub fn bessel_j(&self, n: usize, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match &self.inner {
            Flavor::Value(v) => Self::result_to_ext(v.bessel_j(n, p, rm, cc), v.is_zero(), true),
            Flavor::Inf(_) => NAN,
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    fn bessel_nu_ext(
        &self,
        nu: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
        f: fn(
            &ExactNumNumber,
            &ExactNumNumber,
            usize,
            RoundingMode,
            &mut Consts,
        ) -> Result<ExactNumNumber, Error>,
    ) -> Self {
        match (&self.inner, &nu.inner) {
            (Flavor::Value(x), Flavor::Value(n)) => {
                Self::result_to_ext(f(x, n, p, rm, cc), x.is_zero(), true)
            }
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
            _ => NAN,
        }
    }

    /// \(J_ν(\mathrm{self})\) for real order `nu`.
    ///
    /// # Precision
    ///
    /// - Algorithm: series in `x`; integer `ν` delegates to [`Self::bessel_j`].
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no (identity golds).
    pub fn bessel_j_nu(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.bessel_nu_ext(nu, p, rm, cc, ExactNumNumber::bessel_j_nu)
    }

    /// \(Y_ν(\mathrm{self})\) for `self > 0`.
    ///
    /// # Precision
    ///
    /// - Algorithm: Wronskian / series from \(J_ν\); cut on \((-\infty, 0]\).
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    pub fn bessel_y(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.bessel_nu_ext(nu, p, rm, cc, ExactNumNumber::bessel_y)
    }

    /// \(I_ν(\mathrm{self})\).
    ///
    /// # Precision
    ///
    /// - Algorithm: series; \(I_ν(z) = i^{-ν} J_ν(iz)\) for the complex path.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    pub fn bessel_i(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.bessel_nu_ext(nu, p, rm, cc, ExactNumNumber::bessel_i)
    }

    /// \(K_ν(\mathrm{self})\) for `self > 0`. \(K_{-ν}=K_ν\).
    ///
    /// # Precision
    ///
    /// - Algorithm: series / Temme; large-`|x|` asymptotic `k_asymptotic`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    pub fn bessel_k(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.bessel_nu_ext(nu, p, rm, cc, ExactNumNumber::bessel_k)
    }
    gen_wrapper_arg_rm_cc!(
        "Complete elliptic `K(self)`. Parameter `m = k²`. `m = 1` is `+∞`; `m > 1` uses the reciprocal-modulus transform.

# Precision

- Algorithm: Carlson `R_F` duplication; cap `CARLSON_DUPE_MAX = 128`.
- Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
- MPFR oracle: no (Accumath `K` identities).",
        elliptic_k,
        Self,
        { NAN },
        { NAN },
        p,
        usize
    );
    gen_wrapper_arg_rm_cc!(
        "Complete elliptic `E(self)` for `self ≤ 1`. `E(1) = 1`.

# Precision

- Algorithm: Carlson `R_F` / `R_D`; `CARLSON_DUPE_MAX = 128`.
- Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
- MPFR oracle: no.",
        elliptic_e_complete,
        Self,
        { NAN },
        { NAN },
        p,
        usize
    );
    /// Incomplete `F(self | m)` for `|self| ≤ 1`. `self = sin φ`, `m = k²`.
    ///
    /// # Precision
    ///
    /// - Algorithm: Carlson `R_F`; `CARLSON_DUPE_MAX = 128`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    pub fn elliptic_f(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match (&self.inner, &m.inner) {
            (Flavor::Value(x), Flavor::Value(mv)) => {
                Self::result_to_ext(x.elliptic_f(mv, p, rm, cc), x.is_zero(), true)
            }
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
            _ => NAN,
        }
    }
    /// Incomplete `E(self | m)` for `|self| ≤ 1`.
    ///
    /// # Precision
    ///
    /// - Algorithm: Carlson `R_F` / `R_D`; `CARLSON_DUPE_MAX = 128`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    pub fn elliptic_e(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match (&self.inner, &m.inner) {
            (Flavor::Value(x), Flavor::Value(mv)) => {
                Self::result_to_ext(x.elliptic_e(mv, p, rm, cc), x.is_zero(), true)
            }
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
            _ => NAN,
        }
    }
    /// Complete `Π(self, m)` for `self < 1`, `m < 1`.
    ///
    /// # Precision
    ///
    /// - Algorithm: Carlson `R_J`; `CARLSON_DUPE_MAX = 128`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    pub fn elliptic_pi_complete(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        match (&self.inner, &m.inner) {
            (Flavor::Value(n), Flavor::Value(mv)) => {
                Self::result_to_ext(n.elliptic_pi_complete(mv, p, rm, cc), false, true)
            }
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
            _ => NAN,
        }
    }
    /// Incomplete `Π(self; x | m)`.
    ///
    /// # Precision
    ///
    /// - Algorithm: Carlson `R_J`; `CARLSON_DUPE_MAX = 128`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    pub fn elliptic_pi(
        &self,
        x: &Self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        match (&self.inner, &x.inner, &m.inner) {
            (Flavor::Value(n), Flavor::Value(xv), Flavor::Value(mv)) => {
                Self::result_to_ext(n.elliptic_pi(xv, mv, p, rm, cc), xv.is_zero(), true)
            }
            (Flavor::NaN(err), _, _) | (_, Flavor::NaN(err), _) | (_, _, Flavor::NaN(err)) => {
                Self::nan(*err)
            }
            _ => NAN,
        }
    }
    /// Legendre \(P_n(\mathrm{self})\) for integer `n`.
    ///
    /// # Precision
    ///
    /// - Algorithm: three-term recurrence at `p + O(n)` bits. Cap via caller (`LEGENDRE_N_MAX` on Accumath).
    /// - Bound: working-precision recurrence (not a Ziv leaf).
    /// - MPFR oracle: no.
    pub fn legendre_p(&self, n: u32, p: usize, rm: RoundingMode) -> Self {
        match &self.inner {
            Flavor::Value(v) => Self::result_to_ext(v.legendre_p(n, p, rm), false, true),
            Flavor::Inf(_) => NAN,
            Flavor::NaN(err) => Self::nan(*err),
        }
    }
    /// Associated \(P_n^m(\mathrm{self})\) (Condon–Shortley).
    ///
    /// # Precision
    ///
    /// - Algorithm: recurrence from \(P_n\); Condon–Shortley phase.
    /// - Bound: working-precision recurrence (not a Ziv leaf).
    /// - MPFR oracle: no.
    pub fn assoc_legendre_p(&self, n: u32, m: i32, p: usize, rm: RoundingMode) -> Self {
        match &self.inner {
            Flavor::Value(v) => Self::result_to_ext(v.assoc_legendre_p(n, m, p, rm), false, true),
            Flavor::Inf(_) => NAN,
            Flavor::NaN(err) => Self::nan(*err),
        }
    }
    /// Gaussian \({}_2F_1(\mathrm{self}, b; c; z)\).
    ///
    /// # Precision
    ///
    /// - Algorithm: series for `|z| < 1`; Gauss at `z = 1`; Pfaff / continuation. Cap `HYPERGEOM_TERM_MAX = 10_000`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`) when the series converges.
    /// - MPFR oracle: no.
    pub fn hypergeom_2f1(
        &self,
        b: &Self,
        c: &Self,
        z: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        match (&self.inner, &b.inner, &c.inner, &z.inner) {
            (Flavor::Value(a), Flavor::Value(bv), Flavor::Value(cv), Flavor::Value(zv)) => {
                Self::result_to_ext(a.hypergeom_2f1(bv, cv, zv, p, rm, cc), zv.is_zero(), true)
            }
            (Flavor::NaN(err), _, _, _)
            | (_, Flavor::NaN(err), _, _)
            | (_, _, Flavor::NaN(err), _)
            | (_, _, _, Flavor::NaN(err)) => Self::nan(*err),
            _ => NAN,
        }
    }
    /// Regularized incomplete beta \(I_x(a=\mathrm{self}, b)\).
    ///
    /// # Precision
    ///
    /// - Algorithm: series / continued fraction in `x ∈ [0, 1]` for `a > 0`, `b > 0`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    pub fn betainc(&self, b: &Self, x: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match (&self.inner, &b.inner, &x.inner) {
            (Flavor::Value(a), Flavor::Value(bv), Flavor::Value(xv)) => {
                Self::result_to_ext(a.betainc(bv, xv, p, rm, cc), xv.is_zero(), true)
            }
            (Flavor::NaN(err), _, _) | (_, Flavor::NaN(err), _) | (_, _, Flavor::NaN(err)) => {
                Self::nan(*err)
            }
            _ => NAN,
        }
    }
    gen_wrapper_arg_rm_cc!(
        "Computes the hyperbolic arcsine of a number with precision `p`. The result is rounded using the rounding mode `rm`.
        This function requires constants cache `cc` for computing the result.
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.",
        asinh,
        Self,
        { INF_POS },
        { INF_NEG },
        p,
        usize
    );
    gen_wrapper_arg_rm_cc!(
        "Computes the hyperbolic arccosine of a number with precision `p`. The result is rounded using the rounding mode `rm`.
        This function requires constants cache `cc` for computing the result.
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.",
        acosh,
        Self,
        { INF_POS },
        { NAN },
        p,
        usize
    );
    gen_wrapper_arg_rm_cc!(
        "Computes the hyperbolic arctangent of a number with precision `p`. The result is rounded using the rounding mode `rm`.
        This function requires constants cache `cc` for computing the result.
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.",
        atanh,
        Self,
        { NAN },
        { NAN },
        p,
        usize
    );
}

macro_rules! impl_int_conv {
    ($s:ty, $from_s:ident) => {
        impl ExactNum {
            /// Constructs ExactNum with precision `p` from an integer value `i`.
            /// Precision is rounded upwards to the word size.
            /// The function returns NaN if the precision `p` is incorrect.
            pub fn $from_s(i: $s, p: usize) -> Self {
                Self::result_to_ext(ExactNumNumber::$from_s(i, p), false, true)
            }
        }
    };
}

impl_int_conv!(i8, from_i8);
impl_int_conv!(i16, from_i16);
impl_int_conv!(i32, from_i32);
impl_int_conv!(i64, from_i64);
impl_int_conv!(i128, from_i128);

impl_int_conv!(u8, from_u8);
impl_int_conv!(u16, from_u16);
impl_int_conv!(u32, from_u32);
impl_int_conv!(u64, from_u64);
impl_int_conv!(u128, from_u128);

impl From<ExactNumNumber> for ExactNum {
    fn from(x: ExactNumNumber) -> Self {
        ExactNum {
            inner: Flavor::Value(x),
        }
    }
}

#[cfg(feature = "std")]
use core::{
    fmt::{Binary, Display, Formatter, Octal, UpperHex},
    str::FromStr,
};

use core::{cmp::Eq, cmp::Ordering, cmp::PartialEq, cmp::PartialOrd, ops::Neg};

impl Neg for ExactNum {
    type Output = ExactNum;
    fn neg(mut self) -> Self::Output {
        self.inv_sign();
        self
    }
}

impl Neg for &ExactNum {
    type Output = ExactNum;
    fn neg(self) -> Self::Output {
        let mut ret = self.clone();
        ret.inv_sign();
        ret
    }
}

//
// ordering traits
//

impl PartialEq for ExactNum {
    fn eq(&self, other: &Self) -> bool {
        let cmp_result = ExactNum::cmp(self, other);
        matches!(cmp_result, Some(0))
    }
}

impl<'a> PartialEq<&'a ExactNum> for ExactNum {
    fn eq(&self, other: &&'a ExactNum) -> bool {
        let cmp_result = ExactNum::cmp(self, other);
        matches!(cmp_result, Some(0))
    }
}

impl Eq for ExactNum {}

impl PartialOrd for ExactNum {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let cmp_result = ExactNum::cmp(self, other);
        match cmp_result {
            Some(v) => {
                if v > 0 {
                    Some(Ordering::Greater)
                } else if v < 0 {
                    Some(Ordering::Less)
                } else {
                    Some(Ordering::Equal)
                }
            }
            None => None,
        }
    }
}

impl<'a> PartialOrd<&'a ExactNum> for ExactNum {
    fn partial_cmp(&self, other: &&'a ExactNum) -> Option<Ordering> {
        let cmp_result = ExactNum::cmp(self, other);
        match cmp_result {
            Some(v) => {
                if v > 0 {
                    Some(Ordering::Greater)
                } else if v < 0 {
                    Some(Ordering::Less)
                } else {
                    Some(Ordering::Equal)
                }
            }
            None => None,
        }
    }
}

impl Default for ExactNum {
    fn default() -> ExactNum {
        ExactNum::new(DEFAULT_P)
    }
}

#[cfg(feature = "std")]
impl FromStr for ExactNum {
    type Err = Error;

    /// Returns parsed number or NAN in case of error.
    /// The implementation is not available in no_std environment.
    fn from_str(src: &str) -> Result<ExactNum, Self::Err> {
        let bf = crate::common::consts::TENPOWERS.with(|tp| {
            let cc = &mut tp.borrow_mut();
            ExactNum::parse(src, Radix::Dec, usize::MAX, RoundingMode::ToEven, cc)
        });

        if bf.is_nan() {
            if let Some(err) = bf.err() {
                return Err(err);
            }
        }

        Ok(bf)
    }
}

macro_rules! impl_from {
    ($tt:ty, $fn:ident) => {
        impl From<$tt> for ExactNum {
            fn from(v: $tt) -> Self {
                ExactNum::$fn(v, DEFAULT_P)
            }
        }
    };
}

impl_from!(i8, from_i8);
impl_from!(i16, from_i16);
impl_from!(i32, from_i32);
impl_from!(i64, from_i64);
impl_from!(i128, from_i128);
impl_from!(u8, from_u8);
impl_from!(u16, from_u16);
impl_from!(u32, from_u32);
impl_from!(u64, from_u64);
impl_from!(u128, from_u128);

#[cfg(feature = "std")]
macro_rules! impl_format_rdx {
    ($trait:ty, $rdx:path) => {
        impl $trait for ExactNum {
            /// Formats the number.
            /// The implementation is not available in no_std environment.
            fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), core::fmt::Error> {
                crate::common::consts::TENPOWERS.with(|tp| {
                    let cc = &mut tp.borrow_mut();
                    self.write_str(f, $rdx, RoundingMode::ToEven, cc)
                })
            }
        }
    };
}

#[cfg(feature = "std")]
impl_format_rdx!(Binary, Radix::Bin);
#[cfg(feature = "std")]
impl_format_rdx!(Octal, Radix::Oct);
#[cfg(feature = "std")]
impl_format_rdx!(Display, Radix::Dec);
#[cfg(feature = "std")]
impl_format_rdx!(core::fmt::LowerExp, Radix::Dec);
#[cfg(feature = "std")]
impl_format_rdx!(UpperHex, Radix::Hex);
#[cfg(feature = "std")]
impl core::fmt::UpperExp for ExactNum {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> Result<(), core::fmt::Error> {
        crate::common::consts::TENPOWERS.with(|tp| {
            let cc = &mut tp.borrow_mut();
            let mut s = String::new();
            self.write_str(&mut s, Radix::Dec, RoundingMode::ToEven, cc)?;
            f.write_str(&s.replace('e', "E"))
        })
    }
}
#[cfg(feature = "std")]
impl core::fmt::LowerHex for ExactNum {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> Result<(), core::fmt::Error> {
        crate::common::consts::TENPOWERS.with(|tp| {
            let cc = &mut tp.borrow_mut();
            let mut s = String::new();
            self.write_str(&mut s, Radix::Hex, RoundingMode::ToEven, cc)?;
            if matches!(s.as_str(), "Inf" | "-Inf" | "NaN" | "Err") {
                f.write_str(&s)
            } else {
                f.write_str(&s.to_ascii_lowercase())
            }
        })
    }
}

macro_rules! impl_exact_binop {
    ($trait:ident, $method:ident, $op:ident) => {
        impl core::ops::$trait<&ExactNum> for &ExactNum {
            type Output = ExactNum;

            fn $method(self, rhs: &ExactNum) -> ExactNum {
                ExactNum::$op(self, rhs, DEFAULT_P, RoundingMode::ToEven)
            }
        }

        impl core::ops::$trait<ExactNum> for &ExactNum {
            type Output = ExactNum;

            fn $method(self, rhs: ExactNum) -> ExactNum {
                ExactNum::$op(self, &rhs, DEFAULT_P, RoundingMode::ToEven)
            }
        }

        impl core::ops::$trait<&ExactNum> for ExactNum {
            type Output = ExactNum;

            fn $method(self, rhs: &ExactNum) -> ExactNum {
                ExactNum::$op(&self, rhs, DEFAULT_P, RoundingMode::ToEven)
            }
        }

        impl core::ops::$trait<ExactNum> for ExactNum {
            type Output = ExactNum;

            fn $method(self, rhs: ExactNum) -> ExactNum {
                ExactNum::$op(&self, &rhs, DEFAULT_P, RoundingMode::ToEven)
            }
        }
    };
}

impl_exact_binop!(Add, add, add);
impl_exact_binop!(Sub, sub, sub);
impl_exact_binop!(Mul, mul, mul);
impl_exact_binop!(Div, div, div);

/// A trait for conversion with additional arguments.
pub trait FromExt<T> {
    /// Converts `v` to ExactNum with precision `p` using rounding mode `rm`.
    fn from_ext(v: T, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self;
}

impl<T> FromExt<T> for ExactNum
where
    ExactNum: From<T>,
{
    fn from_ext(v: T, p: usize, rm: RoundingMode, _cc: &mut Consts) -> Self {
        let mut ret = ExactNum::from(v);
        if let Err(err) = ret.set_precision(p, rm) {
            ExactNum::nan(Some(err))
        } else {
            ret
        }
    }
}

impl FromExt<&str> for ExactNum {
    fn from_ext(v: &str, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        ExactNum::parse(v, crate::Radix::Dec, p, rm, cc)
    }
}

#[cfg(test)]
mod tests {

    use crate::common::util::rand_p;
    use crate::defs::DEFAULT_P;
    use crate::ext::ONE;
    use crate::ext::TWO;
    use crate::Consts;
    use crate::Error;
    use crate::ExactNum;
    use crate::Radix;
    use crate::Sign;
    use crate::Word;
    use crate::INF_NEG;
    use crate::INF_POS;
    use crate::NAN;
    use crate::{defs::RoundingMode, WORD_BIT_SIZE};

    use core::num::FpCategory;
    #[cfg(feature = "std")]
    use std::str::FromStr;

    #[cfg(not(feature = "std"))]
    use alloc::format;

    #[cfg(target_pointer_width = "32")]
    #[test]
    fn test_decimal_formatting_round_trip() {
        // Regression test: on 32-bit targets, decimal formatting must produce a string
        // that parses back to the original value without losing precision.
        let mut cc = Consts::new().unwrap();
        let p = 53;
        let rm = RoundingMode::ToEven;
        let value = ExactNum::parse("1.0", Radix::Dec, p, rm, &mut cc);

        let formatted = value.format(Radix::Dec, rm, &mut cc).unwrap();
        assert_eq!(formatted, "1.e+0");

        let reparsed = ExactNum::parse(&formatted, Radix::Dec, p, rm, &mut cc);
        assert_eq!(reparsed, value);
    }

    #[test]
    fn test_ext() {
        let rm = RoundingMode::ToOdd;
        let mut cc = Consts::new().unwrap();

        // Inf & NaN
        let d1 = ExactNum::from_u8(1, rand_p());
        assert!(!d1.is_inf());
        assert!(!d1.is_nan());
        assert!(!d1.is_inf_pos());
        assert!(!d1.is_inf_neg());
        assert!(d1.is_positive());

        let mut d1 = d1.div(&ExactNum::new(rand_p()), rand_p(), rm);
        assert!(d1.is_inf());
        assert!(!d1.is_nan());
        assert!(d1.is_inf_pos());
        assert!(!d1.is_inf_neg());
        assert!(d1.is_positive());

        d1.inv_sign();
        assert!(d1.is_inf());
        assert!(!d1.is_nan());
        assert!(!d1.is_inf_pos());
        assert!(d1.is_inf_neg());
        assert!(d1.is_negative());

        let d1 = ExactNum::new(rand_p()).div(&ExactNum::new(rand_p()), rand_p(), rm);
        assert!(!d1.is_inf());
        assert!(d1.is_nan());
        assert!(!d1.is_inf_pos());
        assert!(!d1.is_inf_neg());
        assert!(d1.sign().is_none());

        for _ in 0..1000 {
            let i = crate::common::test_rng::random::<i64>();
            let d1 = ExactNum::from_i64(i, rand_p());
            let n1 = ExactNum::parse(&format!("{}", i), Radix::Dec, rand_p(), rm, &mut cc);
            assert!(d1.cmp(&n1) == Some(0));

            let i = crate::common::test_rng::random::<u64>();
            let d1 = ExactNum::from_u64(i, rand_p());
            let n1 = ExactNum::parse(&format!("{}", i), Radix::Dec, rand_p(), rm, &mut cc);
            assert!(d1.cmp(&n1) == Some(0));

            let i = crate::common::test_rng::random::<i128>();
            let d1 = ExactNum::from_i128(i, rand_p());
            let n1 = ExactNum::parse(&format!("{}", i), Radix::Dec, rand_p(), rm, &mut cc);
            assert!(d1.cmp(&n1) == Some(0));

            let i = crate::common::test_rng::random::<u128>();
            let d1 = ExactNum::from_u128(i, rand_p());
            let n1 = ExactNum::parse(&format!("{}", i), Radix::Dec, rand_p(), rm, &mut cc);
            assert!(d1.cmp(&n1) == Some(0));
        }

        assert!(ONE.exponent().is_some());
        assert!(INF_POS.exponent().is_none());
        assert!(INF_NEG.exponent().is_none());
        assert!(NAN.exponent().is_none());

        assert!(ONE.as_raw_parts().is_some());
        assert!(INF_POS.as_raw_parts().is_none());
        assert!(INF_NEG.as_raw_parts().is_none());
        assert!(NAN.as_raw_parts().is_none());

        assert!(ONE.add(&ONE, rand_p(), rm).cmp(&TWO) == Some(0));
        assert!(ONE.add(&INF_POS, rand_p(), rm).is_inf_pos());
        assert!(INF_POS.add(&ONE, rand_p(), rm).is_inf_pos());
        assert!(ONE.add(&INF_NEG, rand_p(), rm).is_inf_neg());
        assert!(INF_NEG.add(&ONE, rand_p(), rm).is_inf_neg());
        assert!(INF_POS.add(&INF_POS, rand_p(), rm).is_inf_pos());
        assert!(INF_POS.add(&INF_NEG, rand_p(), rm).is_nan());
        assert!(INF_NEG.add(&INF_NEG, rand_p(), rm).is_inf_neg());
        assert!(INF_NEG.add(&INF_POS, rand_p(), rm).is_nan());

        assert!(ONE.add_full_prec(&ONE).cmp(&TWO) == Some(0));
        assert!(ONE.add_full_prec(&INF_POS).is_inf_pos());
        assert!(INF_POS.add_full_prec(&ONE).is_inf_pos());
        assert!(ONE.add_full_prec(&INF_NEG).is_inf_neg());
        assert!(INF_NEG.add_full_prec(&ONE).is_inf_neg());
        assert!(INF_POS.add_full_prec(&INF_POS).is_inf_pos());
        assert!(INF_POS.add_full_prec(&INF_NEG).is_nan());
        assert!(INF_NEG.add_full_prec(&INF_NEG).is_inf_neg());
        assert!(INF_NEG.add_full_prec(&INF_POS).is_nan());

        assert!(ONE.sub_full_prec(&ONE).is_zero());
        assert!(ONE.sub_full_prec(&INF_POS).is_inf_neg());
        assert!(INF_POS.sub_full_prec(&ONE).is_inf_pos());
        assert!(ONE.sub_full_prec(&INF_NEG).is_inf_pos());
        assert!(INF_NEG.sub_full_prec(&ONE).is_inf_neg());
        assert!(INF_POS.sub_full_prec(&INF_POS).is_nan());
        assert!(INF_POS.sub_full_prec(&INF_NEG).is_inf_pos());
        assert!(INF_NEG.sub_full_prec(&INF_NEG).is_nan());
        assert!(INF_NEG.sub_full_prec(&INF_POS).is_inf_neg());

        assert!(ONE.mul_full_prec(&ONE).cmp(&ONE) == Some(0));
        assert!(ONE.mul_full_prec(&INF_POS).is_inf_pos());
        assert!(INF_POS.mul_full_prec(&ONE).is_inf_pos());
        assert!(ONE.mul_full_prec(&INF_NEG).is_inf_neg());
        assert!(INF_NEG.mul_full_prec(&ONE).is_inf_neg());
        assert!(INF_POS.mul_full_prec(&INF_POS).is_inf_pos());
        assert!(INF_POS.mul_full_prec(&INF_NEG).is_inf_neg());
        assert!(INF_NEG.mul_full_prec(&INF_NEG).is_inf_pos());
        assert!(INF_NEG.mul_full_prec(&INF_POS).is_inf_neg());

        assert!(TWO.sub(&ONE, rand_p(), rm).cmp(&ONE) == Some(0));
        assert!(ONE.sub(&INF_POS, rand_p(), rm).is_inf_neg());
        assert!(INF_POS.sub(&ONE, rand_p(), rm).is_inf_pos());
        assert!(ONE.sub(&INF_NEG, rand_p(), rm).is_inf_pos());
        assert!(INF_NEG.sub(&ONE, rand_p(), rm).is_inf_neg());
        assert!(INF_POS.sub(&INF_POS, rand_p(), rm).is_nan());
        assert!(INF_POS.sub(&INF_NEG, rand_p(), rm).is_inf_pos());
        assert!(INF_NEG.sub(&INF_NEG, rand_p(), rm).is_nan());
        assert!(INF_NEG.sub(&INF_POS, rand_p(), rm).is_inf_neg());

        assert!(TWO.mul(&ONE, rand_p(), rm).cmp(&TWO) == Some(0));
        assert!(ONE.mul(&INF_POS, rand_p(), rm).is_inf_pos());
        assert!(INF_POS.mul(&ONE, rand_p(), rm).is_inf_pos());
        assert!(ONE.mul(&INF_NEG, rand_p(), rm).is_inf_neg());
        assert!(INF_NEG.mul(&ONE, rand_p(), rm).is_inf_neg());
        assert!(ONE.neg().mul(&INF_POS, rand_p(), rm).is_inf_neg());
        assert!(ONE.neg().mul(&INF_NEG, rand_p(), rm).is_inf_pos());
        assert!(INF_POS.mul(&ONE.neg(), rand_p(), rm).is_inf_neg());
        assert!(INF_NEG.mul(&ONE.neg(), rand_p(), rm).is_inf_pos());
        assert!(INF_POS.mul(&INF_POS, rand_p(), rm).is_inf_pos());
        assert!(INF_POS.mul(&INF_NEG, rand_p(), rm).is_inf_neg());
        assert!(INF_NEG.mul(&INF_NEG, rand_p(), rm).is_inf_pos());
        assert!(INF_NEG.mul(&INF_POS, rand_p(), rm).is_inf_neg());
        assert!(INF_POS.mul(&ExactNum::new(rand_p()), rand_p(), rm).is_nan());
        assert!(INF_NEG.mul(&ExactNum::new(rand_p()), rand_p(), rm).is_nan());
        assert!(ExactNum::new(rand_p()).mul(&INF_POS, rand_p(), rm).is_nan());
        assert!(ExactNum::new(rand_p()).mul(&INF_NEG, rand_p(), rm).is_nan());

        assert!(TWO.div(&TWO, rand_p(), rm).cmp(&ONE) == Some(0));
        assert!(TWO.div(&INF_POS, rand_p(), rm).is_zero());
        assert!(INF_POS.div(&TWO, rand_p(), rm).is_inf_pos());
        assert!(TWO.div(&INF_NEG, rand_p(), rm).is_zero());
        assert!(INF_NEG.div(&TWO, rand_p(), rm).is_inf_neg());
        assert!(TWO.neg().div(&INF_POS, rand_p(), rm).is_zero());
        assert!(TWO.neg().div(&INF_NEG, rand_p(), rm).is_zero());
        assert!(INF_POS.div(&TWO.neg(), rand_p(), rm).is_inf_neg());
        assert!(INF_NEG.div(&TWO.neg(), rand_p(), rm).is_inf_pos());
        assert!(INF_POS.div(&INF_POS, rand_p(), rm).is_nan());
        assert!(INF_POS.div(&INF_NEG, rand_p(), rm).is_nan());
        assert!(INF_NEG.div(&INF_NEG, rand_p(), rm).is_nan());
        assert!(INF_NEG.div(&INF_POS, rand_p(), rm).is_nan());
        assert!(INF_POS
            .div(&ExactNum::new(rand_p()), rand_p(), rm)
            .is_inf_pos());
        assert!(INF_NEG
            .div(&ExactNum::new(rand_p()), rand_p(), rm)
            .is_inf_neg());
        assert!(ExactNum::new(rand_p())
            .div(&INF_POS, rand_p(), rm)
            .is_zero());
        assert!(ExactNum::new(rand_p())
            .div(&INF_NEG, rand_p(), rm)
            .is_zero());

        assert!(TWO.rem(&TWO).is_zero());
        assert!(TWO.rem(&INF_POS).cmp(&TWO) == Some(0));
        assert!(INF_POS.rem(&TWO).is_nan());
        assert!(TWO.rem(&INF_NEG).cmp(&TWO) == Some(0));
        assert!(INF_NEG.rem(&TWO).is_nan());
        assert!(TWO.neg().rem(&INF_POS).cmp(&TWO.neg()) == Some(0));
        assert!(TWO.neg().rem(&INF_NEG).cmp(&TWO.neg()) == Some(0));
        assert!(INF_POS.rem(&TWO.neg()).is_nan());
        assert!(INF_NEG.rem(&TWO.neg()).is_nan());
        assert!(INF_POS.rem(&INF_POS).is_nan());
        assert!(INF_POS.rem(&INF_NEG).is_nan());
        assert!(INF_NEG.rem(&INF_NEG).is_nan());
        assert!(INF_NEG.rem(&INF_POS).is_nan());
        assert!(INF_POS.rem(&ExactNum::new(rand_p())).is_nan());
        assert!(INF_NEG.rem(&ExactNum::new(rand_p())).is_nan());
        assert!(ExactNum::new(rand_p()).rem(&INF_POS).is_zero());
        assert!(ExactNum::new(rand_p()).rem(&INF_NEG).is_zero());

        for op in [ExactNum::add, ExactNum::sub, ExactNum::mul, ExactNum::div] {
            assert!(op(&NAN, &ONE, rand_p(), rm).is_nan());
            assert!(op(&ONE, &NAN, rand_p(), rm).is_nan());
            assert!(op(&NAN, &INF_POS, rand_p(), rm).is_nan());
            assert!(op(&INF_POS, &NAN, rand_p(), rm).is_nan());
            assert!(op(&NAN, &INF_NEG, rand_p(), rm).is_nan());
            assert!(op(&INF_NEG, &NAN, rand_p(), rm).is_nan());
            assert!(op(&NAN, &NAN, rand_p(), rm).is_nan());
        }

        assert!(ExactNum::rem(&NAN, &ONE).is_nan());
        assert!(ExactNum::rem(&ONE, &NAN).is_nan());
        assert!(ExactNum::rem(&NAN, &INF_POS).is_nan());
        assert!(ExactNum::rem(&INF_POS, &NAN).is_nan());
        assert!(ExactNum::rem(&NAN, &INF_NEG).is_nan());
        assert!(ExactNum::rem(&INF_NEG, &NAN).is_nan());
        assert!(ExactNum::rem(&NAN, &NAN).is_nan());

        for op in [ExactNum::add_full_prec, ExactNum::sub_full_prec, ExactNum::mul_full_prec] {
            assert!(op(&NAN, &ONE).is_nan());
            assert!(op(&ONE, &NAN).is_nan());
            assert!(op(&NAN, &INF_POS).is_nan());
            assert!(op(&INF_POS, &NAN).is_nan());
            assert!(op(&NAN, &INF_NEG).is_nan());
            assert!(op(&INF_NEG, &NAN).is_nan());
            assert!(op(&NAN, &NAN).is_nan());
        }

        assert!(ONE.cmp(&ONE).unwrap() == 0);
        assert!(ONE.cmp(&INF_POS).unwrap() < 0);
        assert!(INF_POS.cmp(&ONE).unwrap() > 0);
        assert!(INF_POS.cmp(&INF_POS).unwrap() == 0);
        assert!(ONE.cmp(&INF_NEG).unwrap() > 0);
        assert!(INF_NEG.cmp(&ONE).unwrap() < 0);
        assert!(INF_NEG.cmp(&INF_NEG).unwrap() == 0);
        assert!(INF_POS.cmp(&INF_NEG).unwrap() > 0);
        assert!(INF_NEG.cmp(&INF_POS).unwrap() < 0);
        assert!(INF_POS.cmp(&INF_POS).unwrap() == 0);
        assert!(ONE.cmp(&NAN).is_none());
        assert!(NAN.cmp(&ONE).is_none());
        assert!(INF_POS.cmp(&NAN).is_none());
        assert!(NAN.cmp(&INF_POS).is_none());
        assert!(INF_NEG.cmp(&NAN).is_none());
        assert!(NAN.cmp(&INF_NEG).is_none());
        assert!(NAN.cmp(&NAN).is_none());

        assert!(ONE.abs_cmp(&ONE).unwrap() == 0);
        assert!(ONE.abs_cmp(&INF_POS).unwrap() < 0);
        assert!(INF_POS.abs_cmp(&ONE).unwrap() > 0);
        assert!(INF_POS.abs_cmp(&INF_POS).unwrap() == 0);
        assert!(ONE.abs_cmp(&INF_NEG).unwrap() < 0);
        assert!(INF_NEG.abs_cmp(&ONE).unwrap() > 0);
        assert!(INF_NEG.abs_cmp(&INF_NEG).unwrap() == 0);
        assert!(INF_POS.abs_cmp(&INF_NEG).unwrap() == 0);
        assert!(INF_NEG.abs_cmp(&INF_POS).unwrap() == 0);
        assert!(INF_POS.abs_cmp(&INF_POS).unwrap() == 0);
        assert!(ONE.abs_cmp(&NAN).is_none());
        assert!(NAN.abs_cmp(&ONE).is_none());
        assert!(INF_POS.abs_cmp(&NAN).is_none());
        assert!(NAN.abs_cmp(&INF_POS).is_none());
        assert!(INF_NEG.abs_cmp(&NAN).is_none());
        assert!(NAN.abs_cmp(&INF_NEG).is_none());
        assert!(NAN.abs_cmp(&NAN).is_none());

        assert!(ONE.is_positive());
        assert!(!ONE.is_negative());

        assert!(ONE.neg().is_negative());
        assert!(!ONE.neg().is_positive());
        assert!(!INF_POS.is_negative());
        assert!(INF_POS.is_positive());
        assert!(INF_NEG.is_negative());
        assert!(!INF_NEG.is_positive());
        assert!(!NAN.is_positive());
        assert!(!NAN.is_negative());

        assert!(ONE.pow(&ONE, rand_p(), rm, &mut cc).cmp(&ONE) == Some(0));
        assert!(ExactNum::new(DEFAULT_P)
            .pow(&INF_POS, rand_p(), rm, &mut cc)
            .is_zero());
        assert!(ExactNum::new(DEFAULT_P)
            .pow(&INF_NEG, rand_p(), rm, &mut cc)
            .is_zero());
        assert!(ONE.pow(&INF_POS, rand_p(), rm, &mut cc).cmp(&ONE) == Some(0));
        assert!(ONE.pow(&INF_NEG, rand_p(), rm, &mut cc).cmp(&ONE) == Some(0));
        assert!(TWO.pow(&INF_POS, rand_p(), rm, &mut cc).is_inf_pos());
        assert!(TWO.pow(&INF_NEG, rand_p(), rm, &mut cc).is_inf_neg());
        assert!(INF_POS.pow(&ONE, rand_p(), rm, &mut cc).is_inf_pos());
        assert!(INF_NEG.pow(&ONE, rand_p(), rm, &mut cc).is_inf_neg());
        assert!(INF_NEG.pow(&TWO, rand_p(), rm, &mut cc).is_inf_pos());
        assert!(INF_POS.pow(&ONE.neg(), rand_p(), rm, &mut cc).is_zero());
        assert!(INF_NEG.pow(&ONE.neg(), rand_p(), rm, &mut cc).is_zero());
        assert!(
            INF_POS
                .pow(&ExactNum::new(DEFAULT_P), rand_p(), rm, &mut cc)
                .cmp(&ONE)
                == Some(0)
        );
        assert!(
            INF_NEG
                .pow(&ExactNum::new(DEFAULT_P), rand_p(), rm, &mut cc)
                .cmp(&ONE)
                == Some(0)
        );
        assert!(INF_POS.pow(&INF_POS, rand_p(), rm, &mut cc).is_inf_pos());
        assert!(INF_NEG.pow(&INF_POS, rand_p(), rm, &mut cc).is_inf_pos());
        assert!(INF_POS.pow(&INF_NEG, rand_p(), rm, &mut cc).is_zero());
        assert!(INF_NEG.pow(&INF_NEG, rand_p(), rm, &mut cc).is_zero());

        let half = ONE.div(&TWO, rand_p(), rm);
        assert!(TWO.log(&TWO, rand_p(), rm, &mut cc).cmp(&ONE) == Some(0));
        assert!(TWO.log(&INF_POS, rand_p(), rm, &mut cc).is_zero());
        assert!(TWO.log(&INF_NEG, rand_p(), rm, &mut cc).is_nan());
        assert!(INF_POS.log(&TWO, rand_p(), rm, &mut cc).is_inf_pos());
        assert!(INF_NEG.log(&TWO, rand_p(), rm, &mut cc).is_nan());
        assert!(half.log(&half, rand_p(), rm, &mut cc).cmp(&ONE) == Some(0));
        assert!(half.log(&INF_POS, rand_p(), rm, &mut cc).is_zero());
        assert!(half.log(&INF_NEG, rand_p(), rm, &mut cc).is_nan());
        assert!(INF_POS.log(&half, rand_p(), rm, &mut cc).is_inf_neg());
        assert!(INF_NEG.log(&half, rand_p(), rm, &mut cc).is_nan());
        assert!(INF_POS.log(&INF_POS, rand_p(), rm, &mut cc).is_nan());
        assert!(INF_POS.log(&INF_NEG, rand_p(), rm, &mut cc).is_nan());
        assert!(INF_NEG.log(&INF_POS, rand_p(), rm, &mut cc).is_nan());
        assert!(INF_NEG.log(&INF_NEG, rand_p(), rm, &mut cc).is_nan());
        assert!(TWO.log(&ONE, rand_p(), rm, &mut cc).is_inf_pos());
        assert!(half.log(&ONE, rand_p(), rm, &mut cc).is_inf_pos());
        assert!(ONE.log(&ONE, rand_p(), rm, &mut cc).is_nan());

        assert!(ONE.pow(&NAN, rand_p(), rm, &mut cc).is_nan());
        assert!(NAN.pow(&ONE, rand_p(), rm, &mut cc).is_nan());
        assert!(INF_POS.pow(&NAN, rand_p(), rm, &mut cc).is_nan());
        assert!(NAN.pow(&INF_POS, rand_p(), rm, &mut cc).is_nan());
        assert!(INF_NEG.pow(&NAN, rand_p(), rm, &mut cc).is_nan());
        assert!(NAN.pow(&INF_NEG, rand_p(), rm, &mut cc).is_nan());
        assert!(NAN.pow(&NAN, rand_p(), rm, &mut cc).is_nan());

        assert!(NAN.powi(2, rand_p(), rm).is_nan());
        assert!(NAN.powi(0, rand_p(), rm).is_nan());
        assert!(INF_POS.powi(2, rand_p(), rm).is_inf_pos());
        assert!(INF_POS.powi(3, rand_p(), rm).is_inf_pos());
        assert!(INF_NEG.powi(4, rand_p(), rm).is_inf_pos());
        assert!(INF_NEG.powi(5, rand_p(), rm).is_inf_neg());
        assert!(INF_POS.powi(0, rand_p(), rm).cmp(&ONE) == Some(0));
        assert!(INF_NEG.powi(0, rand_p(), rm).cmp(&ONE) == Some(0));

        assert!(TWO.log(&NAN, rand_p(), rm, &mut cc).is_nan());
        assert!(NAN.log(&TWO, rand_p(), rm, &mut cc).is_nan());
        assert!(INF_POS.log(&NAN, rand_p(), rm, &mut cc).is_nan());
        assert!(NAN.log(&INF_POS, rand_p(), rm, &mut cc).is_nan());
        assert!(INF_NEG.log(&NAN, rand_p(), rm, &mut cc).is_nan());
        assert!(NAN.log(&INF_NEG, rand_p(), rm, &mut cc).is_nan());
        assert!(NAN.log(&NAN, rand_p(), rm, &mut cc).is_nan());

        assert!(INF_NEG.abs().is_inf_pos());
        assert!(INF_POS.abs().is_inf_pos());
        assert!(NAN.abs().is_nan());

        assert!(INF_NEG.int().is_nan());
        assert!(INF_POS.int().is_nan());
        assert!(NAN.int().is_nan());

        assert!(INF_NEG.fract().is_nan());
        assert!(INF_POS.fract().is_nan());
        assert!(NAN.fract().is_nan());

        assert!(INF_NEG.ceil().is_inf_neg());
        assert!(INF_POS.ceil().is_inf_pos());
        assert!(NAN.ceil().is_nan());

        assert!(INF_NEG.floor().is_inf_neg());
        assert!(INF_POS.floor().is_inf_pos());
        assert!(NAN.floor().is_nan());

        for rm in [
            RoundingMode::Up,
            RoundingMode::Down,
            RoundingMode::ToZero,
            RoundingMode::FromZero,
            RoundingMode::ToEven,
            RoundingMode::ToOdd,
        ] {
            assert!(INF_NEG.round(0, rm).is_inf_neg());
            assert!(INF_POS.round(0, rm).is_inf_pos());
            assert!(NAN.round(0, rm).is_nan());
        }

        assert!(INF_NEG.sqrt(rand_p(), rm).is_nan());
        assert!(INF_POS.sqrt(rand_p(), rm).is_inf_pos());
        assert!(NAN.sqrt(rand_p(), rm).is_nan());

        assert!(INF_NEG.cbrt(rand_p(), rm).is_inf_neg());
        assert!(INF_POS.cbrt(rand_p(), rm).is_inf_pos());
        assert!(NAN.cbrt(rand_p(), rm).is_nan());

        for op in [ExactNum::ln, ExactNum::log2, ExactNum::log10] {
            assert!(op(&INF_NEG, rand_p(), rm, &mut cc).is_nan());
            assert!(op(&INF_POS, rand_p(), rm, &mut cc).is_inf_pos());
            assert!(op(&NAN, rand_p(), rm, &mut cc).is_nan());
        }

        assert!(INF_NEG.exp(rand_p(), rm, &mut cc).is_zero());
        assert!(INF_POS.exp(rand_p(), rm, &mut cc).is_inf_pos());
        assert!(NAN.exp(rand_p(), rm, &mut cc).is_nan());

        assert!(INF_NEG.sin(rand_p(), rm, &mut cc).is_nan());
        assert!(INF_POS.sin(rand_p(), rm, &mut cc).is_nan());
        assert!(NAN.sin(rand_p(), rm, &mut cc).is_nan());

        assert!(INF_NEG.cos(rand_p(), rm, &mut cc).is_nan());
        assert!(INF_POS.cos(rand_p(), rm, &mut cc).is_nan());
        assert!(NAN.cos(rand_p(), rm, &mut cc).is_nan());

        assert!(INF_NEG.tan(rand_p(), rm, &mut cc).is_nan());
        assert!(INF_POS.tan(rand_p(), rm, &mut cc).is_nan());
        assert!(NAN.tan(rand_p(), rm, &mut cc).is_nan());

        assert!(INF_NEG.asin(rand_p(), rm, &mut cc).is_nan());
        assert!(INF_POS.asin(rand_p(), rm, &mut cc).is_nan());
        assert!(NAN.asin(rand_p(), rm, &mut cc).is_nan());

        assert!(INF_NEG.acos(rand_p(), rm, &mut cc).is_nan());
        assert!(INF_POS.acos(rand_p(), rm, &mut cc).is_nan());
        assert!(NAN.acos(rand_p(), rm, &mut cc).is_nan());

        let p = rand_p();
        let mut half_pi: ExactNum = cc.pi_num(p, rm).unwrap().into();
        half_pi.set_exponent(1);
        assert!(INF_NEG.atan(p, rm, &mut cc).cmp(&half_pi.neg()) == Some(0));
        assert!(INF_POS.atan(p, rm, &mut cc).cmp(&half_pi) == Some(0));
        assert!(NAN.atan(rand_p(), rm, &mut cc).is_nan());

        assert!(INF_NEG.sinh(rand_p(), rm, &mut cc).is_inf_neg());
        assert!(INF_POS.sinh(rand_p(), rm, &mut cc).is_inf_pos());
        assert!(NAN.sinh(rand_p(), rm, &mut cc).is_nan());

        assert!(INF_NEG.cosh(rand_p(), rm, &mut cc).is_inf_pos());
        assert!(INF_POS.cosh(rand_p(), rm, &mut cc).is_inf_pos());
        assert!(NAN.cosh(rand_p(), rm, &mut cc).is_nan());

        assert!(INF_NEG.tanh(rand_p(), rm, &mut cc).cmp(&ONE.neg()) == Some(0));
        assert!(INF_POS.tanh(rand_p(), rm, &mut cc).cmp(&ONE) == Some(0));
        assert!(NAN.tanh(rand_p(), rm, &mut cc).is_nan());

        assert!(INF_NEG.asinh(rand_p(), rm, &mut cc).is_inf_neg());
        assert!(INF_POS.asinh(rand_p(), rm, &mut cc).is_inf_pos());
        assert!(NAN.asinh(rand_p(), rm, &mut cc).is_nan());

        assert!(INF_NEG.acosh(rand_p(), rm, &mut cc).is_nan());
        assert!(INF_POS.acosh(rand_p(), rm, &mut cc).is_inf_pos());
        assert!(NAN.acosh(rand_p(), rm, &mut cc).is_nan());

        assert!(INF_NEG.atanh(rand_p(), rm, &mut cc).is_nan());
        assert!(INF_POS.atanh(rand_p(), rm, &mut cc).is_nan());
        assert!(NAN.atanh(rand_p(), rm, &mut cc).is_nan());

        assert!(INF_NEG.reciprocal(rand_p(), rm).is_zero());
        assert!(INF_POS.reciprocal(rand_p(), rm).is_zero());
        assert!(NAN.reciprocal(rand_p(), rm).is_nan());

        assert!(TWO.signum().cmp(&ONE) == Some(0));
        assert!(TWO.neg().signum().cmp(&ONE.neg()) == Some(0));
        assert!(INF_POS.signum().cmp(&ONE) == Some(0));
        assert!(INF_NEG.signum().cmp(&ONE.neg()) == Some(0));
        assert!(NAN.signum().is_nan());

        let d1 = ONE.clone();
        assert!(d1.exponent() == Some(1));
        let words: &[Word] = {
            #[cfg(not(target_pointer_width = "32"))]
            {
                &[0, 0x8000000000000000]
            }
            #[cfg(target_pointer_width = "32")]
            {
                &[0, 0, 0, 0x80000000]
            }
        };

        assert!(d1.mantissa_digits() == Some(words));
        assert!(d1.is_inline());
        assert!(d1.mantissa_max_bit_len() == Some(DEFAULT_P));
        assert!(d1.precision() == Some(DEFAULT_P));
        assert!(d1.sign() == Some(Sign::Pos));

        assert!(INF_POS.exponent().is_none());
        assert!(INF_POS.mantissa_digits().is_none());
        assert!(INF_POS.mantissa_max_bit_len().is_none());
        assert!(INF_POS.precision().is_none());
        assert!(INF_POS.sign() == Some(Sign::Pos));

        assert!(INF_NEG.exponent().is_none());
        assert!(INF_NEG.mantissa_digits().is_none());
        assert!(INF_NEG.mantissa_max_bit_len().is_none());
        assert!(INF_NEG.precision().is_none());
        assert!(INF_NEG.sign() == Some(Sign::Neg));

        assert!(NAN.exponent().is_none());
        assert!(NAN.mantissa_digits().is_none());
        assert!(NAN.mantissa_max_bit_len().is_none());
        assert!(NAN.precision().is_none());
        assert!(NAN.sign().is_none());

        INF_POS.clone().set_exponent(1);
        INF_POS.clone().set_precision(1, rm).unwrap();
        INF_POS.clone().set_sign(Sign::Pos);

        INF_NEG.clone().set_exponent(1);
        INF_NEG.clone().set_precision(1, rm).unwrap();
        INF_NEG.clone().set_sign(Sign::Pos);

        NAN.clone().set_exponent(1);
        NAN.clone().set_precision(1, rm).unwrap();
        NAN.clone().set_sign(Sign::Pos);

        assert!(INF_POS.min(&ONE).cmp(&ONE) == Some(0));
        assert!(INF_NEG.min(&ONE).is_inf_neg());
        assert!(NAN.min(&ONE).is_nan());
        assert!(ONE.min(&INF_POS).cmp(&ONE) == Some(0));
        assert!(ONE.min(&INF_NEG).is_inf_neg());
        assert!(ONE.min(&NAN).is_nan());
        assert!(NAN.min(&INF_POS).is_nan());
        assert!(NAN.min(&INF_NEG).is_nan());
        assert!(NAN.min(&NAN).is_nan());
        assert!(INF_NEG.min(&INF_POS).is_inf_neg());
        assert!(INF_POS.min(&INF_NEG).is_inf_neg());
        assert!(INF_POS.min(&INF_POS).is_inf_pos());
        assert!(INF_NEG.min(&INF_NEG).is_inf_neg());

        assert!(INF_POS.max(&ONE).is_inf_pos());
        assert!(INF_NEG.max(&ONE).cmp(&ONE) == Some(0));
        assert!(NAN.max(&ONE).is_nan());
        assert!(ONE.max(&INF_POS).is_inf_pos());
        assert!(ONE.max(&INF_NEG).cmp(&ONE) == Some(0));
        assert!(ONE.max(&NAN).is_nan());
        assert!(NAN.max(&INF_POS).is_nan());
        assert!(NAN.max(&INF_NEG).is_nan());
        assert!(NAN.max(&NAN).is_nan());
        assert!(INF_NEG.max(&INF_POS).is_inf_pos());
        assert!(INF_POS.max(&INF_NEG).is_inf_pos());
        assert!(INF_POS.max(&INF_POS).is_inf_pos());
        assert!(INF_NEG.max(&INF_NEG).is_inf_neg());

        assert!(ONE.clamp(&ONE.neg(), &TWO).cmp(&ONE) == Some(0));
        assert!(ONE.clamp(&TWO, &ONE).is_nan());
        assert!(ONE.clamp(&INF_POS, &ONE).is_nan());
        assert!(ONE.clamp(&TWO, &INF_NEG).is_nan());
        assert!(ONE.neg().clamp(&ONE, &TWO).cmp(&ONE) == Some(0));
        assert!(TWO.clamp(&ONE.neg(), &ONE).cmp(&ONE) == Some(0));
        assert!(INF_POS.clamp(&ONE, &TWO).cmp(&TWO) == Some(0));
        assert!(INF_POS.clamp(&ONE, &INF_POS).is_inf_pos());
        assert!(INF_POS.clamp(&INF_NEG, &ONE).cmp(&ONE) == Some(0));
        assert!(INF_POS.clamp(&NAN, &INF_POS).is_nan());
        assert!(INF_POS.clamp(&ONE, &NAN).is_nan());
        assert!(INF_POS.clamp(&NAN, &NAN).is_nan());
        assert!(INF_NEG.clamp(&ONE, &TWO).cmp(&ONE) == Some(0));
        assert!(INF_NEG.clamp(&ONE, &INF_POS).cmp(&ONE) == Some(0));
        assert!(INF_NEG.clamp(&INF_NEG, &ONE).is_inf_neg());
        assert!(INF_NEG.clamp(&NAN, &INF_POS).is_nan());
        assert!(INF_NEG.clamp(&ONE, &NAN).is_nan());
        assert!(INF_NEG.clamp(&NAN, &NAN).is_nan());
        assert!(NAN.clamp(&ONE, &TWO).is_nan());
        assert!(NAN.clamp(&NAN, &TWO).is_nan());
        assert!(NAN.clamp(&ONE, &NAN).is_nan());
        assert!(NAN.clamp(&NAN, &NAN).is_nan());
        assert!(NAN.clamp(&INF_NEG, &INF_POS).is_nan());

        assert!(ExactNum::min_positive(DEFAULT_P).classify() == FpCategory::Subnormal);
        assert!(INF_POS.classify() == FpCategory::Infinite);
        assert!(INF_NEG.classify() == FpCategory::Infinite);
        assert!(NAN.classify() == FpCategory::Nan);
        assert!(ONE.classify() == FpCategory::Normal);

        assert!(!INF_POS.is_subnormal());
        assert!(!INF_NEG.is_subnormal());
        assert!(!NAN.is_subnormal());
        assert!(ExactNum::min_positive(DEFAULT_P).is_subnormal());
        assert!(!ExactNum::min_positive_normal(DEFAULT_P).is_subnormal());
        assert!(!ExactNum::max_value(DEFAULT_P).is_subnormal());
        assert!(!ExactNum::min_value(DEFAULT_P).is_subnormal());

        let n1 = ExactNum::convert_from_radix(
            Sign::Pos,
            &[],
            0,
            Radix::Dec,
            usize::MAX - 1,
            RoundingMode::None,
            &mut cc,
        );
        assert!(n1.is_nan());
        assert!(n1.err() == Some(Error::InvalidArgument));

        assert!(
            n1.convert_to_radix(Radix::Dec, RoundingMode::None, &mut cc)
                == Err(Error::InvalidArgument)
        );
        assert!(
            INF_POS.convert_to_radix(Radix::Dec, RoundingMode::None, &mut cc)
                == Err(Error::InvalidArgument)
        );
        assert!(
            INF_NEG.convert_to_radix(Radix::Dec, RoundingMode::None, &mut cc)
                == Err(Error::InvalidArgument)
        );
    }

    #[cfg(feature = "std")]
    #[test]
    fn test_ops_std() {
        let mut cc = Consts::new().unwrap();

        let d1 = ExactNum::parse(
            "0.0123456789012345678901234567890123456789",
            Radix::Dec,
            DEFAULT_P,
            RoundingMode::None,
            &mut cc,
        );

        let d1str = format!("{}", d1);
        assert_eq!(&d1str, "1.23456789012345678901234567890123456789e-2");
        assert_eq!(format!("{:e}", d1), d1str);
        assert_eq!(
            format!("{:E}", d1),
            "1.23456789012345678901234567890123456789E-2"
        );
        let mut d2 = ExactNum::from_str(&d1str).unwrap();
        d2.set_precision(DEFAULT_P, RoundingMode::ToEven).unwrap();
        assert_eq!(d2, d1);

        let d1 = ExactNum::parse(
            "-123.456789012345678901234567890123456789",
            Radix::Dec,
            DEFAULT_P,
            RoundingMode::None,
            &mut cc,
        );
        let d1str = format!("{}", d1);
        assert_eq!(&d1str, "-1.23456789012345678901234567890123456789e+2");
        let mut d2 = ExactNum::from_str(&d1str).unwrap();
        d2.set_precision(DEFAULT_P, RoundingMode::ToEven).unwrap();
        assert_eq!(d2, d1);

        let d1str = format!("{}", INF_POS);
        assert_eq!(d1str, "Inf");

        let d1str = format!("{}", INF_NEG);
        assert_eq!(d1str, "-Inf");

        let d1str = format!("{}", NAN);
        assert_eq!(d1str, "NaN");

        assert!(ExactNum::from_str("abc").is_ok());
        assert!(ExactNum::from_str("abc").unwrap().is_nan());
    }

    #[test]
    pub fn test_ops() {
        let mut cc = Consts::new().unwrap();

        let d1 = -&(TWO.clone());
        assert!(d1.is_negative());

        let p = DEFAULT_P;
        let rm = RoundingMode::ToEven;
        let two = ExactNum::from_u8(2, p);
        let eighth = two.powsi(-3, p, rm);
        let expected = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(8, p), p, rm);
        assert_eq!(eighth.cmp(&expected), Some(0));
        assert_eq!(two.powsi(3, p, rm).cmp(&ExactNum::from_u8(8, p)), Some(0));
        assert!(
            ExactNum::from_i8(-123, p) == ExactNum::parse("-1.23e+2", Radix::Dec, p, rm, &mut cc)
        );
        assert!(
            ExactNum::from_u8(123, p) == ExactNum::parse("1.23e+2", Radix::Dec, p, rm, &mut cc)
        );
        assert!(
            ExactNum::from_i16(-12312, p)
                == ExactNum::parse("-1.2312e+4", Radix::Dec, p, rm, &mut cc)
        );
        assert!(
            ExactNum::from_u16(12312, p)
                == ExactNum::parse("1.2312e+4", Radix::Dec, p, rm, &mut cc)
        );
        assert!(
            ExactNum::from_i32(-123456789, p)
                == ExactNum::parse("-1.23456789e+8", Radix::Dec, p, rm, &mut cc)
        );
        assert!(
            ExactNum::from_u32(123456789, p)
                == ExactNum::parse("1.23456789e+8", Radix::Dec, p, rm, &mut cc)
        );
        assert!(
            ExactNum::from_i64(-1234567890123456789, p)
                == ExactNum::parse("-1.234567890123456789e+18", Radix::Dec, p, rm, &mut cc)
        );
        assert!(
            ExactNum::from_u64(1234567890123456789, p)
                == ExactNum::parse("1.234567890123456789e+18", Radix::Dec, p, rm, &mut cc)
        );
        assert!(
            ExactNum::from_i128(-123456789012345678901234567890123456789, p)
                == ExactNum::parse(
                    "-1.23456789012345678901234567890123456789e+38",
                    Radix::Dec,
                    p,
                    rm,
                    &mut cc
                )
        );
        assert!(
            ExactNum::from_u128(123456789012345678901234567890123456789, p)
                == ExactNum::parse(
                    "1.23456789012345678901234567890123456789e+38",
                    Radix::Dec,
                    p,
                    rm,
                    &mut cc
                )
        );

        let neg = ExactNum::from_i8(-3, WORD_BIT_SIZE);
        let pos = ExactNum::from_i8(5, WORD_BIT_SIZE);

        assert!(pos > neg);
        assert!(neg < pos);
        assert!(!(pos < neg));
        assert!(!(neg > pos));
        assert!(INF_NEG < neg);
        assert!(INF_NEG < pos);
        assert!(INF_NEG < INF_POS);
        assert!(!(INF_NEG > neg));
        assert!(!(INF_NEG > pos));
        assert!(!(INF_NEG > INF_POS));
        assert!(INF_POS > neg);
        assert!(INF_POS > pos);
        assert!(INF_POS > INF_NEG);
        assert!(!(INF_POS < neg));
        assert!(!(INF_POS < pos));
        assert!(!(INF_POS < INF_NEG));
        assert!(!(INF_POS > INF_POS));
        assert!(!(INF_POS < INF_POS));
        assert!(!(INF_NEG > INF_NEG));
        assert!(!(INF_NEG < INF_NEG));
        assert!(!(INF_POS > NAN));
        assert!(!(INF_POS < NAN));
        assert!(!(INF_NEG > NAN));
        assert!(!(INF_NEG < NAN));
        assert!(!(NAN > INF_POS));
        assert!(!(NAN < INF_POS));
        assert!(!(NAN > INF_NEG));
        assert!(!(NAN < INF_NEG));
        assert!(!(NAN > NAN));
        assert!(!(NAN < NAN));
        assert!(!(neg > NAN));
        assert!(!(neg < NAN));
        assert!(!(pos > NAN));
        assert!(!(pos < NAN));
        assert!(!(NAN > neg));
        assert!(!(NAN < neg));
        assert!(!(NAN > pos));
        assert!(!(NAN < pos));

        assert!(!(NAN == NAN));
        assert!(!(NAN == INF_POS));
        assert!(!(NAN == INF_NEG));
        assert!(!(INF_POS == NAN));
        assert!(!(INF_NEG == NAN));
        assert!(!(INF_NEG == INF_POS));
        assert!(!(INF_POS == INF_NEG));
        assert!(!(INF_POS == neg));
        assert!(!(INF_POS == pos));
        assert!(!(INF_NEG == neg));
        assert!(!(INF_NEG == pos));
        assert!(!(neg == INF_POS));
        assert!(!(pos == INF_POS));
        assert!(!(neg == INF_NEG));
        assert!(!(pos == INF_NEG));
        assert!(!(pos == neg));
        assert!(!(neg == pos));
        assert!(neg == neg);
        assert!(pos == pos);
        assert!(INF_NEG == INF_NEG);
        assert!(INF_POS == INF_POS);
    }

    #[test]
    fn test_oom_and_large_precision() {
        let oom = ExactNum::nan(Some(Error::MemoryAllocation));
        assert!(oom.is_nan());
        assert_eq!(oom.err(), Some(Error::MemoryAllocation));

        let n = ExactNum::new(usize::MAX);
        assert!(n.is_nan());
        assert_eq!(n.err(), Some(Error::InvalidArgument));

        let p = 128 * WORD_BIT_SIZE;
        let a = ExactNum::from_word(3, p);
        let b = ExactNum::from_word(5, p);
        let s = a.add(&b, p, RoundingMode::ToEven);
        let m = a.mul(&b, p, RoundingMode::ToEven);
        assert!(!s.is_nan(), "large-prec add hung or failed");
        assert!(!m.is_nan(), "large-prec mul hung or failed");
        assert_eq!(s.cmp(&ExactNum::from_word(8, p)), Some(0));
    }

    #[test]
    fn test_two_sum_fused_polyval() {
        let p = 128;
        let rm = RoundingMode::ToEven;
        let one = ExactNum::from_u8(1, p);
        let two = ExactNum::from_u8(2, p);
        let three = ExactNum::from_u8(3, p);

        let (hi, lo) = one.two_sum(&two, p, rm);
        let rec = hi.add(&lo, p, rm);
        assert_eq!(rec.cmp(&ExactNum::from_u8(3, p)), Some(0));

        let (ph, pl) = two.two_product(&three, p, rm);
        let pr = ph.add(&pl, p, rm);
        assert_eq!(pr.cmp(&ExactNum::from_u8(6, p)), Some(0));

        let sum = ExactNum::fused_sum(&[one.clone(), two.clone(), three.clone()], p, rm);
        assert_eq!(sum.cmp(&ExactNum::from_u8(6, p)), Some(0));

        let dot = ExactNum::fused_dot(
            &[one.clone(), two.clone()],
            &[three.clone(), one.clone()],
            p,
            rm,
        );
        assert_eq!(dot.cmp(&ExactNum::from_u8(5, p)), Some(0));

        // 1 + 2x + 3x² at x = 2 → 17
        let pv = ExactNum::polyval(&[one, two, three], &ExactNum::from_u8(2, p), p, rm);
        assert_eq!(pv.cmp(&ExactNum::from_u8(17, p)), Some(0));
    }
}

#[cfg(feature = "random")]
#[cfg(test)]
mod rand_tests {

    use super::*;
    use crate::common::util::TEST_EXP_BOUND;

    #[test]
    fn test_rand() {
        for _ in 0..100 {
            let p = crate::common::test_rng::random::<usize>() % 192 + DEFAULT_P;
            let exp_from = crate::common::test_rng::random::<Exponent>().abs() % TEST_EXP_BOUND;
            let span = (TEST_EXP_BOUND - exp_from).max(1);
            let exp_shift = crate::common::test_rng::random::<Exponent>().abs() % span;
            let exp_to = exp_from + exp_shift;

            let n = ExactNum::random_normal(p, exp_from, exp_to);

            assert!(!n.is_subnormal());
            assert!(n.exponent().unwrap() >= exp_from && n.exponent().unwrap() <= exp_to);
            assert!(n.precision().unwrap() >= p);
        }
    }
}
