// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! ExactNum including finite numbers, NaN, and `Inf`.

use crate::common::util::log2_ceil;
use crate::common::util::round_p;
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

#[cfg(feature = "std")]
use core::fmt::Write;

#[cfg(not(feature = "std"))]
use alloc::{string::String, vec::Vec};

/// Not a number.
pub const NAN: ExactNum = ExactNum {
    inner: Flavor::NaN(None),
};

/// The NaN that `clamp`, `min` and `max` hand out by reference.
static NAN_VALUE: ExactNum = NAN;

/// Positive infinity.
pub const INF_POS: ExactNum = ExactNum {
    inner: Flavor::Inf(Sign::Pos),
};

/// Negative infinity.
pub const INF_NEG: ExactNum = ExactNum {
    inner: Flavor::Inf(Sign::Neg),
};

fn invalid_op() -> ExactNum {
    ExactNum::nan(Some(Error::InvalidArgument))
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

    /// Positive or negative infinity, built fresh (no copy of an existing value).
    fn inf_with_sign(positive: bool) -> Self {
        if positive {
            INF_POS
        } else {
            INF_NEG
        }
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
    /// The result is rounded at exactly `p` bits (storage is still a whole number of words);
    /// with `RoundingMode::None` the result keeps the storage precision.
    /// The function returns NaN if the precision `p` is incorrect.
    ///
    /// Accuracy: correctly-rounded
    pub fn add(&self, d2: &Self, p: usize, rm: RoundingMode) -> Self {
        self.add_op(d2, p, rm, false)
    }

    /// Adds `d2` to `self` and returns the result of the operation.
    /// The resulting precision is equal to the full precision of the result.
    /// This operation can be used to emulate integer addition.
    ///
    /// Accuracy: correctly-rounded
    pub fn add_full_prec(&self, d2: &Self) -> Self {
        self.add_op(d2, 0, RoundingMode::None, true)
    }

    fn add_op(&self, d2: &Self, p: usize, rm: RoundingMode, full_prec: bool) -> Self {
        match &self.inner {
            Flavor::Value(v1) => match &d2.inner {
                Flavor::Value(v2) => Self::result_to_ext_rm(
                    if full_prec {
                        v1.add_full_prec(v2)
                    } else {
                        ExactNumNumber::exact_bits(p, rm, |q, r| v1.add_per_op(v2, q, r))
                    },
                    v1.is_zero(),
                    v1.sign() == v2.sign(),
                    p,
                    rm,
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
                        invalid_op()
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
    /// The result is rounded at exactly `p` bits (storage is still a whole number of words);
    /// with `RoundingMode::None` the result keeps the storage precision.
    /// The function returns NaN if the precision `p` is incorrect.
    ///
    /// Accuracy: correctly-rounded
    pub fn sub(&self, d2: &Self, p: usize, rm: RoundingMode) -> Self {
        self.sub_op(d2, p, rm, false)
    }

    /// Subtracts `d2` from `self` and returns the result of the operation.
    /// The resulting precision is equal to the full precision of the result.
    /// This operation can be used to emulate integer subtraction.
    ///
    /// Accuracy: correctly-rounded
    pub fn sub_full_prec(&self, d2: &Self) -> Self {
        self.sub_op(d2, 0, RoundingMode::None, true)
    }

    fn sub_op(&self, d2: &Self, p: usize, rm: RoundingMode, full_prec: bool) -> Self {
        match &self.inner {
            Flavor::Value(v1) => match &d2.inner {
                Flavor::Value(v2) => Self::result_to_ext_rm(
                    if full_prec {
                        v1.sub_full_prec(v2)
                    } else {
                        ExactNumNumber::exact_bits(p, rm, |q, r| v1.sub_per_op(v2, q, r))
                    },
                    v1.is_zero(),
                    v1.sign() == v2.sign(),
                    p,
                    rm,
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
                        invalid_op()
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
    /// The result is rounded at exactly `p` bits (storage is still a whole number of words);
    /// with `RoundingMode::None` the result keeps the storage precision.
    /// The function returns NaN if the precision `p` is incorrect.
    ///
    /// Accuracy: correctly-rounded
    pub fn mul(&self, d2: &Self, p: usize, rm: RoundingMode) -> Self {
        self.mul_op(d2, p, rm, false)
    }

    /// Multiplies `d2` by `self` and returns the result of the operation.
    /// The resulting precision is equal to the full precision of the result.
    /// This operation can be used to emulate integer multiplication.
    ///
    /// Accuracy: correctly-rounded
    pub fn mul_full_prec(&self, d2: &Self) -> Self {
        self.mul_op(d2, 0, RoundingMode::None, true)
    }

    /// Computes `self * b + c` with precision `p`, rounded once with `rm`.
    ///
    /// Unlike `mul` followed by `add`, the product is not rounded to `p` before the addition.
    ///
    /// Accuracy: correctly-rounded
    pub fn fma(&self, b: &Self, c: &Self, p: usize, rm: RoundingMode) -> Self {
        let x = self.exact_copy_ext();
        let b_c = b.exact_copy_ext();
        let c_c = c.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || b.inexact_accumulated() || c.inexact_accumulated(),
            p,
            rm,
            |q, r| x.fma_word(&b_c, &c_c, q, r),
        )
    }

    /// Word-precision implementation of [`Self::fma`].
    fn fma_word(&self, b: &Self, c: &Self, p: usize, rm: RoundingMode) -> Self {
        if self.is_nan() {
            return Self::nan(self.err());
        }
        if b.is_nan() {
            return Self::nan(b.err());
        }
        if c.is_nan() {
            return Self::nan(c.err());
        }
        match (&self.inner, &b.inner, &c.inner) {
            (Flavor::Value(a), Flavor::Value(bv), Flavor::Value(cv)) => {
                Self::result_to_ext_rm(a.fma(bv, cv, p, rm), false, true, p, rm)
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
    ///
    /// Accuracy: correctly-rounded
    pub fn two_sum(&self, b: &Self, p: usize, rm: RoundingMode) -> (Self, Self) {
        if self.is_nan() {
            return (Self::nan(self.err()), Self::nan(self.err()));
        }
        if b.is_nan() {
            return (Self::nan(b.err()), Self::nan(b.err()));
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
            z.set_inexact(lo.inexact_accumulated());
            return z;
        }
        lo
    }

    /// Two-product: `(hi, lo)` with `hi` rounded to `p` bits using `rm` and `hi + lo` equal to the
    /// exact product of finite operands (via [`mul_full_prec`](Self::mul_full_prec)).
    ///
    /// Accuracy: correctly-rounded
    pub fn two_product(&self, b: &Self, p: usize, rm: RoundingMode) -> (Self, Self) {
        if self.is_nan() {
            return (Self::nan(self.err()), Self::nan(self.err()));
        }
        if b.is_nan() {
            return (Self::nan(b.err()), Self::nan(b.err()));
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
    ///
    /// Accuracy: correctly-rounded
    pub fn fused_sum(xs: &[Self], p: usize, rm: RoundingMode) -> Self {
        let xs_c: Vec<Self> = xs.iter().map(|v| v.exact_copy_ext()).collect();
        Self::round_exact(xs.iter().any(|v| v.inexact_accumulated()), p, rm, |q, r| {
            Self::fused_sum_word(&xs_c, q, r)
        })
    }

    /// Word-precision implementation of [`Self::fused_sum`].
    fn fused_sum_word(xs: &[Self], p: usize, rm: RoundingMode) -> Self {
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
    ///
    /// Accuracy: correctly-rounded
    pub fn fused_dot(xs: &[Self], ys: &[Self], p: usize, rm: RoundingMode) -> Self {
        let xs_c: Vec<Self> = xs.iter().map(|v| v.exact_copy_ext()).collect();
        let ys_c: Vec<Self> = ys.iter().map(|v| v.exact_copy_ext()).collect();
        Self::round_exact(
            xs.iter().any(|v| v.inexact_accumulated())
                || ys.iter().any(|v| v.inexact_accumulated()),
            p,
            rm,
            |q, r| Self::fused_dot_word(&xs_c, &ys_c, q, r),
        )
    }

    /// Word-precision implementation of [`Self::fused_dot`].
    fn fused_dot_word(xs: &[Self], ys: &[Self], p: usize, rm: RoundingMode) -> Self {
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
    ///
    /// Accuracy: correctly-rounded
    pub fn polyval(coeffs: &[Self], x: &Self, p: usize, rm: RoundingMode) -> Self {
        let x_c = x.exact_copy_ext();
        let coeffs_c: Vec<Self> = coeffs.iter().map(|v| v.exact_copy_ext()).collect();
        Self::round_exact(
            x.inexact_accumulated() || coeffs.iter().any(|v| v.inexact_accumulated()),
            p,
            rm,
            |q, r| Self::polyval_word(&coeffs_c, &x_c, q, r),
        )
    }

    /// Word-precision implementation of [`Self::polyval`].
    fn polyval_word(coeffs: &[Self], x: &Self, p: usize, rm: RoundingMode) -> Self {
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
    ///
    /// Accuracy: correctly-rounded
    pub fn mul_add(&self, b: &Self, c: &Self, p: usize, rm: RoundingMode) -> Self {
        let x = self.exact_copy_ext();
        let b_c = b.exact_copy_ext();
        let c_c = c.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || b.inexact_accumulated() || c.inexact_accumulated(),
            p,
            rm,
            |q, r| x.mul_add_word(&b_c, &c_c, q, r),
        )
    }

    /// Word-precision implementation of [`Self::mul_add`].
    fn mul_add_word(&self, b: &Self, c: &Self, p: usize, rm: RoundingMode) -> Self {
        if self.is_nan() {
            return Self::nan(self.err());
        }
        if b.is_nan() {
            return Self::nan(b.err());
        }
        if c.is_nan() {
            return Self::nan(c.err());
        }
        match (&self.inner, &b.inner, &c.inner) {
            (Flavor::Value(a), Flavor::Value(bv), Flavor::Value(cv)) => {
                Self::result_to_ext_rm(a.mul_add(bv, cv, p, rm), false, true, p, rm)
            }
            _ => self.fma(b, c, p, rm),
        }
    }

    fn mul_op(&self, d2: &Self, p: usize, rm: RoundingMode, full_prec: bool) -> Self {
        match &self.inner {
            Flavor::Value(v1) => {
                match &d2.inner {
                    Flavor::Value(v2) => Self::result_to_ext_rm(
                        if full_prec {
                            v1.mul_full_prec(v2)
                        } else {
                            ExactNumNumber::exact_bits(p, rm, |q, r| v1.mul_per_op(v2, q, r))
                        },
                        v1.is_zero(),
                        v1.sign() == v2.sign(),
                        p,
                        rm,
                    ),
                    Flavor::Inf(s2) => {
                        if v1.is_zero() {
                            // 0*inf
                            invalid_op()
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
                            invalid_op()
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
    /// The result is rounded at exactly `p` bits (storage is still a whole number of words);
    /// with `RoundingMode::None` the result keeps the storage precision.
    /// The function returns NaN if the precision `p` is incorrect.
    ///
    /// Accuracy: correctly-rounded
    pub fn div(&self, d2: &Self, p: usize, rm: RoundingMode) -> Self {
        match &self.inner {
            Flavor::Value(v1) => match &d2.inner {
                Flavor::Value(v2) => Self::result_to_ext_rm(
                    ExactNumNumber::exact_bits(p, rm, |q, r| v1.div_per_op(v2, q, r)),
                    v1.is_zero(),
                    v1.sign() == v2.sign(),
                    p,
                    rm,
                ),
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
                Flavor::Inf(_) => invalid_op(),
                Flavor::NaN(err) => Self::nan(*err),
            },
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    /// Returns the remainder of division of `|self|` by `|d2|`. The sign of the result is set to the sign of `self`.
    ///
    /// The remainder is always exact: its last bit is no finer than the finer of the two
    /// operands' last bits and it is smaller than both, so it fits without rounding. It takes
    /// no precision or rounding mode, and `inexact()` is false.
    ///
    /// Accuracy: exact
    pub fn rem(&self, d2: &Self) -> Self {
        match &self.inner {
            Flavor::Value(v1) => match &d2.inner {
                Flavor::Value(v2) => {
                    // The remainder is exact; its flag reports this operation, which never
                    // rounds, and the accumulated flag keeps both operands' history.
                    let inputs = v1.inexact() || v2.inexact();
                    let r = match (v1.exact_copy(), v2.exact_copy()) {
                        (Ok(a), Ok(b)) => a.rem(&b).map(|mut r| {
                            r.finish_flags(inputs, true);
                            r
                        }),
                        (Err(e), _) | (_, Err(e)) => Err(e),
                    };
                    Self::result_to_ext(r, v1.is_zero(), v1.sign() == v2.sign())
                }
                Flavor::Inf(_) => self.clone(),
                Flavor::NaN(err) => Self::nan(*err),
            },
            Flavor::Inf(_) => invalid_op(),
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    /// Compares `self` to `d2`.
    /// Returns positive if `self` > `d2`, negative if `self` < `d2`, zero if `self` == `d2`, None if `self` or `d2` is NaN.
    /// The non-zero results are exactly `1` and `-1`.
    pub fn compare(&self, d2: &ExactNum) -> Option<SignedWord> {
        self.cmp_raw(d2).map(SignedWord::signum)
    }

    /// [`compare`](Self::compare) before normalization: equal exponents return the mantissa word difference.
    fn cmp_raw(&self, d2: &ExactNum) -> Option<SignedWord> {
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
    /// The non-zero results are exactly `1` and `-1`.
    pub fn abs_cmp(&self, d2: &Self) -> Option<SignedWord> {
        self.abs_cmp_raw(d2).map(SignedWord::signum)
    }

    fn abs_cmp_raw(&self, d2: &Self) -> Option<SignedWord> {
        match &self.inner {
            Flavor::Value(v1) => match &d2.inner {
                Flavor::Value(v2) => Some(v1.abs_cmp(v2)),
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
    ///
    /// Accuracy: correctly-rounded
    pub fn pow(&self, n: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let x = self.exact_copy_ext();
        let n_c = n.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || n.inexact_accumulated(),
            p,
            rm,
            |q, r| x.pow_word(&n_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::pow`].
    fn pow_word(&self, n: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match &self.inner {
            Flavor::Value(v1) => {
                match &n.inner {
                    // 0^(-n) is a pole (Inf in every mode); any other overflow follows `rm`.
                    Flavor::Value(v2) if v1.is_zero() => {
                        Self::result_to_ext(v1.pow(v2, p, rm, cc), true, v1.sign() == v2.sign())
                    }
                    Flavor::Value(v2) => Self::result_to_ext_rm(
                        v1.pow(v2, p, rm, cc),
                        false,
                        v1.sign() == v2.sign(),
                        p,
                        rm,
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
    ///
    /// Accuracy: correctly-rounded
    pub fn powi(&self, n: usize, p: usize, rm: RoundingMode) -> Self {
        let x = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            x.powi_word(n, q, r)
        })
    }

    /// Word-precision implementation of [`Self::powi`].
    fn powi_word(&self, n: usize, p: usize, rm: RoundingMode) -> Self {
        match &self.inner {
            Flavor::Value(v1) => Self::result_to_ext_rm(v1.powi(n, p, rm), false, true, p, rm),
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
    ///
    /// Accuracy: correctly-rounded
    pub fn powsi(&self, n: i64, p: usize, rm: RoundingMode) -> Self {
        let x = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            x.powsi_word(n, q, r)
        })
    }

    /// Word-precision implementation of [`Self::powsi`].
    fn powsi_word(&self, n: i64, p: usize, rm: RoundingMode) -> Self {
        match &self.inner {
            Flavor::Value(v1) => Self::result_to_ext_rm(v1.powsi(n, p, rm), false, true, p, rm),
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
    ///
    /// Accuracy: correctly-rounded
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
    pub fn clamp<'a>(&'a self, min: &'a Self, max: &'a Self) -> &'a Self {
        if self.is_nan()
            || min.is_nan()
            || max.is_nan()
            || matches!(max.compare(min), Some(c) if c < 0)
        {
            &NAN_VALUE
        } else if matches!(self.compare(min), Some(c) if c < 0) {
            min
        } else if matches!(self.compare(max), Some(c) if c > 0) {
            max
        } else {
            self
        }
    }

    /// Returns `d1` if `d1` is greater than `self`, or `self` otherwise (a borrow of one of the two
    /// operands, no copy). If either argument is NaN, the function returns NaN.
    pub fn max<'a>(&'a self, d1: &'a Self) -> &'a Self {
        if self.is_nan() || d1.is_nan() {
            &NAN_VALUE
        } else if matches!(self.compare(d1), Some(c) if c < 0) {
            d1
        } else {
            self
        }
    }

    /// Returns `d1` if `d1` is less than `self`, or `self` otherwise (a borrow of one of the two
    /// operands, no copy). If either argument is NaN, the function returns NaN.
    pub fn min<'a>(&'a self, d1: &'a Self) -> &'a Self {
        if self.is_nan() || d1.is_nan() {
            &NAN_VALUE
        } else if matches!(self.compare(d1), Some(c) if c > 0) {
            d1
        } else {
            self
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
    /// The whole string is consumed. Trailing characters, a bare exponent such as `1e` or `1e+`,
    /// and any other malformed text produce `NaN(InvalidArgument)`.
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
        Self::round_exact(false, p, rm, |q, r| Self::parse_word(s, rdx, q, r, cc))
    }

    /// Word-precision implementation of [`Self::parse`].
    fn parse_word(s: &str, rdx: Radix, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
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
    ///
    /// Accuracy: correctly-rounded
    pub fn atan(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match &self.inner {
            Flavor::Value(v) => match Self::exact_operand(v) {
                Ok((x, acc)) => Self::round_exact(acc, p, rm, |q, r| {
                    Self::result_to_ext(x.atan(q, r, cc), x.is_zero(), true)
                }),
                Err(e) => Self::result_to_ext(Err(e), false, true),
            },
            Flavor::Inf(s) => Self::result_to_ext(Self::half_pi(*s, p, rm, cc), false, true),
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    /// Computes `atan2(self, x)` with precision `p` (quadrant-aware arctangent of `self / x`).
    /// The result is rounded using the rounding mode `rm`.
    /// This function requires constants cache `cc`.
    /// Precision is rounded upwards to the word size.
    ///
    /// Accuracy: correctly-rounded
    pub fn atan2(&self, x: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let y_c = self.exact_copy_ext();
        let x_c = x.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || x.inexact_accumulated(),
            p,
            rm,
            |q, r| y_c.atan2_word(&x_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::atan2`].
    fn atan2_word(&self, x: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return Self::nan(self.err());
        }
        if x.is_nan() {
            return Self::nan(x.err());
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
    ///
    /// Accuracy: correctly-rounded
    pub fn hypot(&self, other: &Self, p: usize, rm: RoundingMode) -> Self {
        let x = self.exact_copy_ext();
        let other_c = other.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || other.inexact_accumulated(),
            p,
            rm,
            |q, r| x.hypot_word(&other_c, q, r),
        )
    }

    /// Word-precision implementation of [`Self::hypot`].
    fn hypot_word(&self, other: &Self, p: usize, rm: RoundingMode) -> Self {
        if self.is_inf() || other.is_inf() {
            return INF_POS;
        }
        if self.is_nan() {
            return Self::nan(self.err());
        }
        if other.is_nan() {
            return Self::nan(other.err());
        }
        match (&self.inner, &other.inner) {
            (Flavor::Value(a), Flavor::Value(b)) => {
                Self::result_to_ext_rm(a.hypot(b, p, rm), false, true, p, rm)
            }
            _ => NAN,
        }
    }

    /// Computes `ln(1 + self)` with precision `p`.
    /// The result is rounded using the rounding mode `rm`.
    /// This function requires constants cache `cc`.
    /// Returns `-Inf` for `self == -1`, and NaN if `self < -1`.
    ///
    /// Accuracy: correctly-rounded
    pub fn log1p(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match &self.inner {
            Flavor::Value(v) => match Self::exact_operand(v) {
                Ok((x, acc)) => Self::round_exact(acc, p, rm, |q, r| {
                    Self::result_to_ext(x.log1p(q, r, cc), false, false)
                }),
                Err(e) => Self::result_to_ext(Err(e), false, true),
            },
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
    ///
    /// Accuracy: correctly-rounded
    pub fn expm1(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match &self.inner {
            Flavor::Value(v) => match Self::exact_operand(v) {
                Ok((x, acc)) => Self::round_exact(acc, p, rm, |q, r| {
                    Self::result_to_ext_rm(x.expm1(q, r, cc), false, true, q, r)
                }),
                Err(e) => Self::result_to_ext(Err(e), false, true),
            },
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
    ///
    /// Accuracy: correctly-rounded
    pub fn tanh(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match &self.inner {
            Flavor::Value(v) => match Self::exact_operand(v) {
                Ok((x, acc)) => Self::round_exact(acc, p, rm, |q, r| {
                    Self::result_to_ext(x.tanh(q, r, cc), x.is_zero(), true)
                }),
                Err(e) => Self::result_to_ext(Err(e), false, true),
            },
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

    /// Evaluates `f(precision, rounding)`, a function computed on operands whose inexact flags
    /// were cleared, so that the result is rounded at exactly `p` bits and carries a
    /// per-operation `inexact` flag; `inputs_inexact` is folded into the accumulated flag.
    ///
    /// Word-aligned precisions with a normal result call `f(p, rm)` directly. Otherwise `f`
    /// runs toward zero with a guard word, and the result is rounded to odd and then to `p`
    /// bits with `rm`, which is a single correct rounding. A rounded subnormal result at an
    /// aligned precision takes the same path, so it is not rounded twice.
    pub(crate) fn round_exact<F>(inputs_inexact: bool, p: usize, rm: RoundingMode, mut f: F) -> Self
    where
        F: FnMut(usize, RoundingMode) -> Self,
    {
        let direct = rm == RoundingMode::None || p == 0 || p >= usize::MAX / 4;
        let mut r = if direct {
            f(p, rm)
        } else if p.is_multiple_of(WORD_BIT_SIZE) {
            let r = f(p, rm);
            let subnormal_rounded = matches!(
                &r.inner,
                Flavor::Value(v) if v.exponent() == crate::EXPONENT_MIN && v.inexact()
            );
            if subnormal_rounded {
                Self::round_odd_path(p, rm, &mut f)
            } else {
                r
            }
        } else {
            Self::round_odd_path(p, rm, &mut f)
        };
        if let Flavor::Value(v) = &mut r.inner {
            v.finish_flags(inputs_inexact, true);
        }
        r
    }

    fn round_odd_path<F>(p: usize, rm: RoundingMode, f: &mut F) -> Self
    where
        F: FnMut(usize, RoundingMode) -> Self,
    {
        let wide = f(round_p(p) + WORD_BIT_SIZE, RoundingMode::ToZero);
        match wide.inner {
            Flavor::Value(v) => {
                if v.is_zero() && !v.inexact() {
                    // Exact zero: its sign depends on the rounding mode.
                    return f(p, rm);
                }
                let is_zero = v.is_zero();
                Self::result_to_ext_rm(v.round_odd_to(p, rm), is_zero, true, p, rm)
            }
            // Inf and NaN do not depend on the precision; recompute for the mode's own value.
            Flavor::Inf(_) | Flavor::NaN(_) => f(p, rm),
        }
    }

    /// A copy of `self` with both inexact flags cleared, for [`Self::round_exact`].
    pub(crate) fn exact_copy_ext(&self) -> Self {
        let mut c = self.clone();
        c.set_inexact(false);
        c
    }

    /// The operand of a function evaluated through [`Self::round_exact`]: a copy with both
    /// flags cleared, and its accumulated flag.
    fn exact_operand(v: &ExactNumNumber) -> Result<(ExactNumNumber, bool), Error> {
        Ok((v.exact_copy()?, v.inexact()))
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

    /// Like `result_to_ext`, but an exponent overflow follows `rm` as in IEEE 754 and MPFR:
    /// the largest finite value of precision `p` with the overflow sign when `rm` rounds toward
    /// zero for that sign (`ToZero`, `Down` for a positive result, `Up` for a negative one),
    /// otherwise Inf. The largest finite value is marked inexact.
    fn result_to_ext_rm(
        res: Result<ExactNumNumber, Error>,
        is_dividend_zero: bool,
        is_same_sign: bool,
        p: usize,
        rm: RoundingMode,
    ) -> ExactNum {
        match res {
            Err(Error::ExponentOverflow(s)) => {
                let toward_zero = matches!(
                    (s.is_positive(), rm),
                    (_, RoundingMode::ToZero)
                        | (true, RoundingMode::Down)
                        | (false, RoundingMode::Up)
                );
                if toward_zero {
                    match ExactNumNumber::max_value(p) {
                        Ok(mut m) => {
                            // The largest value representable with exactly `p` bits.
                            if let Err(e) = m.set_precision(p, RoundingMode::ToZero) {
                                return Self::result_to_ext(Err(e), is_dividend_zero, is_same_sign);
                            }
                            m.set_sign(s);
                            m.set_inexact(true);
                            ExactNum {
                                inner: Flavor::Value(m),
                            }
                        }
                        Err(e) => Self::result_to_ext(Err(e), is_dividend_zero, is_same_sign),
                    }
                } else {
                    Self::result_to_ext(res, is_dividend_zero, is_same_sign)
                }
            }
            _ => Self::result_to_ext(res, is_dividend_zero, is_same_sign),
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

    /// `-self` for an operand that has to stay borrowed: fresh `-1` times `self` through the
    /// general multiplication kernel at the operand's own precision. The product is exact, keeps
    /// the sign of zero, Inf and NaN (with its error) and the accumulated inexact flag, and sets
    /// the per-operation flag to false (nothing was rounded); nothing is cloned.
    pub(crate) fn neg_via_mul(&self) -> Self {
        let p = self.mantissa_max_bit_len().unwrap_or(crate::WORD_BIT_SIZE);
        Self::from_i8(-1, p).mul(self, p, RoundingMode::None)
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
    /// - `n` or the size of `m` is too large (larger than i64::MAX / 2 + EXPONENT_MIN).
    /// - `e` is less than EXPONENT_MIN or greater than EXPONENT_MAX.
    pub fn from_raw_parts(m: &[Word], n: usize, s: Sign, e: Exponent, inexact: bool) -> Self {
        Self::result_to_ext(
            ExactNumNumber::from_raw_parts(m, n, s, e, inexact),
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
    ///
    /// Accuracy: correctly-rounded
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
    ///
    /// Accuracy: correctly-rounded
    pub fn reciprocal(&self, p: usize, rm: RoundingMode) -> Self {
        let x = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            x.reciprocal_word(q, r)
        })
    }

    /// Word-precision implementation of [`Self::reciprocal`].
    fn reciprocal_word(&self, p: usize, rm: RoundingMode) -> Self {
        match &self.inner {
            Flavor::Value(v) => {
                Self::result_to_ext_rm(v.reciprocal(p, rm), false, v.is_positive(), p, rm)
            }
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
    ///
    /// For `add`, `sub`, `mul`, and `div` the flag is per operation: it is true only if that
    /// operation rounded, whatever the operands' flags were. Other functions report
    /// [`ExactNum::inexact_accumulated`].
    pub fn inexact(&self) -> bool {
        if let Flavor::Value(v) = &self.inner {
            v.op_inexact()
        } else {
            false
        }
    }

    /// Returns true if `self`, or any value it was computed from, was rounded. The function
    /// returns false if `self` is Inf or NaN. This accumulated flag is what the `expr!` macro
    /// and the correct-rounding loops use to know whether a chain of operations is exact.
    pub fn inexact_accumulated(&self) -> bool {
        if let Flavor::Value(v) = &self.inner {
            v.inexact()
        } else {
            false
        }
    }

    /// Sets the per-operation flag to `op` and adds `inputs` to the accumulated flag, which
    /// otherwise stays this part's own. A zero part that this operation did not round keeps
    /// no history, as with a real product by an exact zero (Ziv loops rely on it).
    pub(crate) fn finish_part(&mut self, op: bool, inputs: bool) {
        if let Flavor::Value(v) = &mut self.inner {
            let inherit = inputs && !(v.is_zero() && !v.inexact());
            v.finish_flags(inherit, true);
            v.set_op_inexact(op);
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
    pub fn frexp(self) -> (Self, Exponent) {
        match self.inner {
            Flavor::Value(v) => match v.frexp() {
                Ok((m, e)) => (m.into(), e),
                Err(err) => (Self::nan(Some(err)), 0),
            },
            Flavor::Inf(s) => (Self::inf_with_sign(s.is_positive()), 0),
            Flavor::NaN(err) => (Self::nan(err), 0),
        }
    }

    /// `self · 2^n`. Alias of [`Self::scalb`]. Takes ownership; the exponent is changed in the operand's own
    /// storage.
    ///
    /// Accuracy: correctly-rounded
    pub fn ldexp(self, n: Exponent, p: usize, rm: RoundingMode) -> Self {
        let inputs = self.inexact_accumulated();
        match self.inner {
            Flavor::Inf(s) => Self::inf_with_sign(s.is_positive()),
            Flavor::NaN(err) => Self::nan(err),
            Flavor::Value(mut v) => {
                v.set_inexact(false);
                Self::ldexp_value(v, inputs, n, p, rm)
            }
        }
    }

    /// `ldexp` of a finite operand whose flags are cleared; `inputs` is its accumulated flag.
    /// Zero is handled first (its sign depends on the mode); every other operand takes one call at `p`
    /// bits, or one call at a wider width rounded toward zero followed by a round-to-odd step (see
    /// [`Self::round_exact`]).
    fn ldexp_value(
        v: ExactNumNumber,
        inputs: bool,
        n: Exponent,
        p: usize,
        rm: RoundingMode,
    ) -> Self {
        let zero = v.is_zero();
        let direct = rm == RoundingMode::None || p == 0 || p >= usize::MAX / 4;
        let mut r = if direct || zero {
            Self::result_to_ext_rm(v.ldexp(n, p, rm), zero, true, p, rm)
        } else {
            let wp = round_p(p) + WORD_BIT_SIZE;
            let wide = Self::result_to_ext_rm(
                v.ldexp(n, wp, RoundingMode::ToZero),
                false,
                true,
                wp,
                RoundingMode::ToZero,
            );
            match wide.inner {
                Flavor::Value(w) => {
                    let is_zero = w.is_zero();
                    Self::result_to_ext_rm(w.round_odd_to(p, rm), is_zero, true, p, rm)
                }
                other => ExactNum { inner: other },
            }
        };
        if let Flavor::Value(v) = &mut r.inner {
            v.finish_flags(inputs, true);
        }
        r
    }

    /// `self · 2^n` (IEEE `scalbn`). Takes ownership, like [`Self::ldexp`].
    ///
    /// Accuracy: correctly-rounded
    pub fn scalb(self, n: Exponent, p: usize, rm: RoundingMode) -> Self {
        self.ldexp(n, p, rm)
    }

    /// `floor(log2(|self|))` as a float. Zero becomes `-Inf`; Inf/NaN unchanged in kind.
    ///
    /// Accuracy: correctly-rounded
    pub fn logb(&self, p: usize, rm: RoundingMode) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.logb_word(q, r)
        })
    }

    /// Word-precision implementation of [`Self::logb`].
    fn logb_word(&self, p: usize, rm: RoundingMode) -> Self {
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

macro_rules! gen_wrapper_exact_rm {
    // rounded at exactly `p` bits through `round_exact`, with a per-operation flag
    ($comment:literal, $fname:ident, $ret:ty, $pos_inf:block, $neg_inf:block, $p:ident, usize) => {
        #[doc=$comment]
        pub fn $fname(&self, $p: usize, rm: RoundingMode) -> $ret {
            match &self.inner {
                Flavor::Value(v) => match Self::exact_operand(v) {
                    Ok((x, acc)) => Self::round_exact(acc, $p, rm, |q, r| {
                        Self::result_to_ext(x.$fname(q, r), x.is_zero(), true)
                    }),
                    Err(e) => Self::result_to_ext(Err(e), false, true),
                },
                Flavor::Inf(s) => if s.is_positive() $pos_inf else $neg_inf,
                Flavor::NaN(err) => Self::nan(*err),
            }
        }
    };
}

macro_rules! gen_wrapper_arg_rm_cc {
    // rounded at exactly `p` bits through `round_exact`, with a per-operation flag
    ($comment:literal, $fname:ident, $ret:ty, $pos_inf:block, $neg_inf:block, $p:ident, usize) => {
        #[doc=$comment]
        pub fn $fname(&self, $p: usize, rm: RoundingMode, cc: &mut Consts) -> $ret {
            match &self.inner {
                Flavor::Value(v) => match Self::exact_operand(v) {
                    Ok((x, acc)) => Self::round_exact(acc, $p, rm, |q, r| {
                        Self::result_to_ext(x.$fname(q, r, cc), x.is_zero(), true)
                    }),
                    Err(e) => Self::result_to_ext(Err(e), false, true),
                },
                Flavor::Inf(s) => if s.is_positive() $pos_inf else $neg_inf,
                Flavor::NaN(err) => Self::nan(*err),
            }
        }
    };
}

macro_rules! gen_wrapper_word_rm_cc {
    // a reduction that is not a rounded function value: word precision, accumulated flag
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

macro_rules! gen_wrapper_ovf_rm_cc {
    // like gen_wrapper_arg_rm_cc, for functions whose only ExponentOverflow away from zero is a
    // finite result past the range: that overflow follows the rounding mode
    ($comment:literal, $fname:ident, $ret:ty, $pos_inf:block, $neg_inf:block, $p:ident, usize) => {
        #[doc=$comment]
        pub fn $fname(&self, $p: usize, rm: RoundingMode, cc: &mut Consts) -> $ret {
            match &self.inner {
                Flavor::Value(v) => match Self::exact_operand(v) {
                    Ok((x, acc)) => Self::round_exact(acc, $p, rm, |q, r| {
                        if x.is_zero() {
                            // a pole at zero (Γ, Ei) is an exact infinity in every mode
                            Self::result_to_ext(x.$fname(q, r, cc), true, true)
                        } else {
                            Self::result_to_ext_rm(x.$fname(q, r, cc), false, true, q, r)
                        }
                    }),
                    Err(e) => Self::result_to_ext(Err(e), false, true),
                },
                Flavor::Inf(s) => if s.is_positive() $pos_inf else $neg_inf,
                Flavor::NaN(err) => Self::nan(*err),
            }
        }
    };
}

macro_rules! gen_wrapper_log {
    ($comment:literal, $fname:ident, $ret:ty, $pos_inf:block, $neg_inf:block, $p:ident, usize) => {
        #[doc=$comment]
        pub fn $fname(&self, $p: usize, rm: RoundingMode, cc: &mut Consts) -> $ret {
            match &self.inner {
                Flavor::Value(v) => {
                    if v.is_zero() {
                        return INF_NEG;
                    }
                    match Self::exact_operand(v) {
                        Ok((x, acc)) => Self::round_exact(acc, $p, rm, |q, r| {
                            Self::result_to_ext(x.$fname(q, r, cc), false, true)
                        }),
                        Err(e) => Self::result_to_ext(Err(e), false, true),
                    }
                },
                Flavor::Inf(s) => if s.is_positive() $pos_inf else $neg_inf,
                Flavor::NaN(err) => Self::nan(*err),
            }
        }
    };
}

impl ExactNum {
    /// Returns the absolute value of `self`. Takes ownership; nothing is copied.
    ///
    /// Accuracy: correctly-rounded
    pub fn abs(mut self) -> Self {
        if self.is_negative() {
            self.inv_sign();
        }
        self
    }

    /// Returns a value with the magnitude of `self` and the sign of `sign`.
    ///
    /// Accuracy: correctly-rounded
    pub fn copysign(&self, sign: &Self, p: usize, rm: RoundingMode) -> Self {
        if self.is_nan() {
            return Self::nan(self.err());
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
    ///
    /// Accuracy: correctly-rounded
    pub fn next_after(&self, toward: &Self, p: usize, rm: RoundingMode) -> Self {
        if self.is_nan() {
            return Self::nan(self.err());
        }
        if toward.is_nan() {
            return Self::nan(toward.err());
        }
        match (&self.inner, &toward.inner) {
            (Flavor::Value(v), Flavor::Value(t)) => {
                Self::result_to_ext(v.next_after(t, p, rm), false, true)
            }
            (Flavor::Inf(_), _) | (_, Flavor::Inf(_)) => self.clone(),
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
        }
    }
    gen_wrapper_arg!(
        "Returns the integer part of `self`.
Accuracy: correctly-rounded",
        int,
        Self,
        { NAN },
        { NAN },
    );
    gen_wrapper_arg!(
        "Returns the fractional part of `self`.
Accuracy: correctly-rounded",
        fract,
        Self,
        { NAN },
        { NAN },
    );
    gen_wrapper_arg!(
        "Returns the smallest integer greater than or equal to `self`.
Accuracy: correctly-rounded",
        ceil,
        Self,
        { INF_POS },
        { INF_NEG },
    );
    gen_wrapper_arg!(
        "Returns the largest integer less than or equal to `self`.
Accuracy: correctly-rounded",
        floor,
        Self,
        { INF_POS },
        { INF_NEG },
    );
    gen_wrapper_arg_rm!("Returns the rounded number with `n` binary positions in the fractional part of the number using rounding mode `rm`.
Accuracy: correctly-rounded", 
        round,
        Self,
        { INF_POS },
        { INF_NEG },
        n,
        usize
    );
    /// Computes the square root of a number with precision `p`. The result is rounded using
    /// the rounding mode `rm` at exactly `p` bits; `inexact()` reports this operation only.
    /// The function returns NaN if the precision `p` is incorrect.
    ///
    /// Accuracy: correctly-rounded
    pub fn sqrt(&self, p: usize, rm: RoundingMode) -> Self {
        match &self.inner {
            Flavor::Value(v) => match Self::exact_operand(v) {
                Ok((x, acc)) => Self::round_exact(acc, p, rm, |q, r| {
                    Self::result_to_ext(x.sqrt(q, r), x.is_zero(), true)
                }),
                Err(e) => Self::result_to_ext(Err(e), false, true),
            },
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
    gen_wrapper_exact_rm!(
        "Computes the cube root of a number with precision `p`. The result is rounded using the rounding mode `rm`.
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.
Accuracy: correctly-rounded",
        cbrt,
        Self,
        { INF_POS },
        { INF_NEG },
        p,
        usize
    );
    /// Computes the `n`-th root of `self` with precision `p`. `n = 2` and `n = 3` delegate to [`sqrt`](Self::sqrt) and [`cbrt`](Self::cbrt).
    ///
    /// Accuracy: correctly-rounded
    pub fn nth_root(&self, n: usize, p: usize, rm: RoundingMode) -> Self {
        let x = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            x.nth_root_word(n, q, r)
        })
    }

    /// Word-precision implementation of [`Self::nth_root`].
    fn nth_root_word(&self, n: usize, p: usize, rm: RoundingMode) -> Self {
        if n == 0 {
            return Self::nan(Some(Error::InvalidArgument));
        }
        match &self.inner {
            Flavor::Value(v) => Self::result_to_ext(v.nth_root(n, p, rm), v.is_zero(), true),
            Flavor::Inf(s) => {
                if n.is_multiple_of(2) {
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
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.
Accuracy: correctly-rounded",
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
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.
Accuracy: correctly-rounded",
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
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.
Accuracy: correctly-rounded",
        log10,
        Self,
        { INF_POS },
        { NAN },
        p,
        usize
    );
    gen_wrapper_ovf_rm_cc!(
        "Computes `e` to the power of `self` with precision `p`. The result is rounded using the rounding mode `rm`.
        This function requires constants cache `cc` for computing the result.
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.
Accuracy: correctly-rounded",
        exp,
        Self,
        { INF_POS },
        { Self::new(p) },
        p,
        usize
    );
    gen_wrapper_ovf_rm_cc!(
        "Computes `2` to the power of `self` with precision `p`. The result is rounded using the rounding mode `rm`.
        This function requires constants cache `cc` for computing the result.
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.
Accuracy: correctly-rounded",
        exp2,
        Self,
        { INF_POS },
        { Self::new(p) },
        p,
        usize
    );
    gen_wrapper_ovf_rm_cc!(
        "Computes `10` to the power of `self` with precision `p`. The result is rounded using the rounding mode `rm`.
        This function requires constants cache `cc` for computing the result.
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.
Accuracy: correctly-rounded",
        exp10,
        Self,
        { INF_POS },
        { Self::new(p) },
        p,
        usize
    );
    gen_wrapper_word_rm_cc!(
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
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.
Accuracy: correctly-rounded",
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
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.
Accuracy: correctly-rounded",
        cos,
        Self,
        { NAN },
        { NAN },
        p,
        usize
    );
    /// Computes `(sin(self), cos(self))` with precision `p` using a shared argument reduction.
    ///
    /// Accuracy: correctly-rounded
    pub fn sin_cos(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> (Self, Self) {
        // Each result goes through its own exact-p-bit rounding and per-operation flag.
        (self.sin(p, rm, cc), self.cos(p, rm, cc))
    }
    gen_wrapper_arg_rm_cc!(
        "Computes the tangent of a number with precision `p`. The result is rounded using the rounding mode `rm`.
        This function requires constants cache `cc` for computing the result.
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.
Accuracy: correctly-rounded",
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
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.
Accuracy: correctly-rounded", 
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
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.
Accuracy: correctly-rounded",
        acos,
        Self,
        { NAN },
        { NAN },
        p,
        usize
    );

    gen_wrapper_ovf_rm_cc!(
        "Computes the hyperbolic sine of a number with precision `p`. The result is rounded using the rounding mode `rm`.
        This function requires constants cache cc for computing the result. 
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.
Accuracy: correctly-rounded",
        sinh,
        Self,
        { INF_POS },
        { INF_NEG },
        p,
        usize
    );
    gen_wrapper_ovf_rm_cc!(
        "Computes the hyperbolic cosine of a number with precision `p`. The result is rounded using the rounding mode `rm`.
        This function requires constants cache cc for computing the result. 
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.
Accuracy: correctly-rounded",
        cosh,
        Self,
        { INF_POS },
        { INF_POS },
        p,
        usize
    );
    /// Computes `(sinh(self), cosh(self))` with precision `p` using a single `exp(|x|)` evaluation.
    ///
    /// Accuracy: correctly-rounded
    pub fn sinh_cosh(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> (Self, Self) {
        // Each result goes through its own exact-p-bit rounding, per-operation flag and
        // mode-aware overflow.
        (self.sinh(p, rm, cc), self.cosh(p, rm, cc))
    }
    gen_wrapper_arg_rm_cc!(
        "Error function `erf(self)` with precision `p`.

# Precision

- Algorithm: Taylor series when `|x|.exponent() ≤ 2`; complementary asymptotic otherwise. Saturates to `±1` when `2|e| > p+4`.
- Bound: Ziv correct-rounding (`MAX_PREC_RETRY`). 1 ULP vs MPFR on `|x| ≲ 4`.
- Thresholds: exponent cut `≤ 2` (not a named constant).
- MPFR oracle: yes, `|x| ≲ 4` under `mpfr-tests`. Complex `erf` on the real axis uses the same oracle; GNU MPC has no `mpc_erf`.
Accuracy: correctly-rounded",
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
- MPFR oracle: yes, `|x| ≲ 4` under `mpfr-tests`.
Accuracy: correctly-rounded",
        erfc,
        Self,
        { Self::new(p) },
        { ExactNum::from_u8(2, p) },
        p,
        usize
    );
    gen_wrapper_ovf_rm_cc!(
        "Gamma function `Γ(self)` with precision `p`. Poles at non-positive integers yield NaN (or +Inf at 0).

# Precision

- Algorithm: Stirling series for `ln Γ` then `exp`; reflection across the negative axis. Integer factorials for small positive integers.
- Bound: Ziv correct-rounding (`MAX_PREC_RETRY`). 1 ULP vs MPFR on the oracle domain.
- MPFR oracle: yes, under `mpfr-tests`.
Accuracy: correctly-rounded",
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
- MPFR oracle: yes, under `mpfr-tests`.
Accuracy: correctly-rounded",
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
- Bound: Ziv correct-rounding (`MAX_PREC_RETRY`). 1 ULP vs MPFR `digamma` on `z > 0`.
- MPFR oracle: yes, `z > 0` under `mpfr-tests`.
Accuracy: correctly-rounded",
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
    /// - MPFR oracle: identity `γ(s,x)=Γ(s)−Γ(s,x)`; upper uses `mpfr_gamma_inc`.
    ///
    /// Accuracy: correctly-rounded
    pub fn gammainc(&self, x: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let x_c = x.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || x.inexact_accumulated(),
            p,
            rm,
            |q, r| this.gammainc_word(&x_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::gammainc`].
    fn gammainc_word(&self, x: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
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
    /// - MPFR oracle: yes, `mpfr_gamma_inc` under `mpfr-tests`.
    ///
    /// Accuracy: correctly-rounded
    pub fn gammainc_upper(&self, x: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let x_c = x.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || x.inexact_accumulated(),
            p,
            rm,
            |q, r| this.gammainc_upper_word(&x_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::gammainc_upper`].
    fn gammainc_upper_word(&self, x: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
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
    gen_wrapper_ovf_rm_cc!(
        "Exponential integral `Ei(self)` (principal value for `self < 0`). `0` is a pole.

# Precision

- Algorithm: power series for moderate `|x|`; factorial asymptotic when `|x|` is large (`exponent() > 6` and `|x| ≳ 0.7 p`).
- Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
- MPFR oracle: yes, `mpfr_eint` under `mpfr-tests`.
Accuracy: correctly-rounded",
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
    /// - MPFR oracle: no (`Si` odd, `Si(0)=0`, `Si(+∞)=π/2`).
    ///
    /// Accuracy: correctly-rounded
    pub fn si(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.si_word(q, r, cc)
        })
    }

    /// Word-precision implementation of [`Self::si`].
    fn si_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
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
- MPFR oracle: no (identity / series golds; GNU MPFR has no `Si`/`Ci`).
Accuracy: correctly-rounded",
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
- MPFR oracle: no (`li(e)=Ei(1)` identity).
Accuracy: correctly-rounded",
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
    /// - MPFR oracle: no (`S` odd, `S(0)=0`, `S(+∞)=1/2`).
    ///
    /// Accuracy: correctly-rounded
    pub fn fresnel_s(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.fresnel_s_word(q, r, cc)
        })
    }

    /// Word-precision implementation of [`Self::fresnel_s`].
    fn fresnel_s_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.fresnel_sc_ext(true, p, rm, cc)
    }

    /// Fresnel cosine integral `C(self)`. `±∞ → ±1/2`.
    ///
    /// # Precision
    ///
    /// - Algorithm: same series / auxiliary `f,g` split as [`Self::fresnel_s`].
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no (same identities as [`Self::fresnel_s`]).
    ///
    /// Accuracy: correctly-rounded
    pub fn fresnel_c(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.fresnel_c_word(q, r, cc)
        })
    }

    /// Word-precision implementation of [`Self::fresnel_c`].
    fn fresnel_c_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
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
    /// - MPFR oracle: yes, `mpfr_ai` under `mpfr-tests`.
    ///
    /// Accuracy: correctly-rounded
    pub fn ai(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.ai_word(q, r, cc)
        })
    }

    /// Word-precision implementation of [`Self::ai`].
    fn ai_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
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
    /// - MPFR oracle: no (`Ai Bi' − Ai' Bi = 1/π`; GNU MPFR has no `Bi`).
    ///
    /// Accuracy: correctly-rounded
    pub fn bi(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.bi_word(q, r, cc)
        })
    }

    /// Word-precision implementation of [`Self::bi`].
    fn bi_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
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
    ///
    /// Accuracy: correctly-rounded
    pub fn ai_prime(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.ai_prime_word(q, r, cc)
        })
    }

    /// Word-precision implementation of [`Self::ai_prime`].
    fn ai_prime_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
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
    ///
    /// Accuracy: correctly-rounded
    pub fn bi_prime(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.bi_prime_word(q, r, cc)
        })
    }

    /// Word-precision implementation of [`Self::bi_prime`].
    fn bi_prime_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
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
    /// - Algorithm: power series; Miller recurrence for large `n` (`n ≤ 1024`);
    ///   Hankel expansion (DLMF 10.17, stopped at its smallest term) when
    ///   `|x| ≥ max(16, 0.35(p+112))` and `|n|² < 2|x|`, and whenever the series
    ///   guard would exceed `BESSEL_GUARD_MAX` under that order condition.
    ///   `J_n(-x)=(-1)^n J_n(x)`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`). 1 ULP vs MPFR `jn` for `n = 0,1,2`.
    /// - MPFR oracle: yes, `n = 0,1,2` under `mpfr-tests`.
    ///
    /// Accuracy: correctly-rounded
    pub fn bessel_j(&self, n: usize, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.bessel_j_word(n, q, r, cc)
        })
    }

    /// Word-precision implementation of [`Self::bessel_j`].
    fn bessel_j_word(&self, n: usize, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
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

    /// \(J_ν(\mathrm{self})\) for real order `nu`. Negative `self` only for integer `ν`.
    ///
    /// # Precision
    ///
    /// - Algorithm: series in `x`; integer `ν` delegates to [`Self::bessel_j`];
    ///   Hankel expansion (same switch as [`Self::bessel_j`]) for large `|x|`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no (identity golds).
    ///
    /// Accuracy: correctly-rounded
    pub fn bessel_j_nu(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let nu_c = nu.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || nu.inexact_accumulated(),
            p,
            rm,
            |q, r| this.bessel_j_nu_word(&nu_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::bessel_j_nu`].
    fn bessel_j_nu_word(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.bessel_nu_ext(nu, p, rm, cc, ExactNumNumber::bessel_j_nu)
    }

    /// \(Y_ν(\mathrm{self})\) for `self > 0`.
    ///
    /// # Precision
    ///
    /// - Algorithm: Wronskian / series from \(J_ν\); Hankel expansion (DLMF 10.17)
    ///   for large `x` (same switch as [`Self::bessel_j`]). Cut on \((-\infty, 0]\).
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: yes, `mpfr_yn` for `n = 0,1` under `mpfr-tests`.
    ///
    /// Accuracy: correctly-rounded
    pub fn bessel_y(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let nu_c = nu.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || nu.inexact_accumulated(),
            p,
            rm,
            |q, r| this.bessel_y_word(&nu_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::bessel_y`].
    fn bessel_y_word(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.bessel_nu_ext(nu, p, rm, cc, ExactNumNumber::bessel_y)
    }

    /// \(I_ν(\mathrm{self})\).
    ///
    /// # Precision
    ///
    /// - Algorithm: series; \(I_ν(z) = i^{-ν} J_ν(iz)\) for the complex path.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn bessel_i(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let nu_c = nu.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || nu.inexact_accumulated(),
            p,
            rm,
            |q, r| this.bessel_i_word(&nu_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::bessel_i`].
    fn bessel_i_word(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.bessel_nu_ext(nu, p, rm, cc, ExactNumNumber::bessel_i)
    }

    /// \(K_ν(\mathrm{self})\) for `self > 0`. \(K_{-ν}=K_ν\).
    ///
    /// # Precision
    ///
    /// - Algorithm: series / Temme; large-`|x|` asymptotic `k_asymptotic`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn bessel_k(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let nu_c = nu.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || nu.inexact_accumulated(),
            p,
            rm,
            |q, r| this.bessel_k_word(&nu_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::bessel_k`].
    fn bessel_k_word(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.bessel_nu_ext(nu, p, rm, cc, ExactNumNumber::bessel_k)
    }
    gen_wrapper_arg_rm_cc!(
        "Complete elliptic `K(self)`. Parameter `m = k²`. `m = 1` is `+∞`; `m > 1` uses the reciprocal-modulus transform.

# Precision

- Algorithm: Carlson `R_F` duplication; cap `max(128, p/4 + 64)` steps.
- Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
- MPFR oracle: no (identity golds; GNU MPFR has no Carlson `K`).
Accuracy: correctly-rounded",
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

- Algorithm: Carlson `R_F` / `R_D`; cap `max(128, p/4 + 64)` steps.
- Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
- MPFR oracle: no.
Accuracy: correctly-rounded",
        elliptic_e_complete,
        Self,
        { NAN },
        { NAN },
        p,
        usize
    );
    /// Incomplete elliptic integral of the first kind, `F(φ|m)`, with
    /// `self = sin φ` and `m = k²`.
    ///
    /// `self` is the sine of the amplitude, not the amplitude `φ`. An angle is
    /// not accepted in this argument: `|self| > 1` is `NaN(InvalidArgument)`,
    /// so a value past `π/2` (for example `2`) cannot be passed here.
    /// `elliptic_f(sin φ, m)` is `F(φ|m)` for `φ` in `[-π/2, π/2]`, and the
    /// function is odd in `self`. `elliptic_f(sin θ, m)` for `θ` outside that
    /// range is `F(arcsin(sin θ)|m)`, which is not `F(θ|m)`.
    /// `m` is the parameter `k²`, not the elliptic modulus `k`.
    /// `F(±1 | 1) = ±∞` (the logarithmic singularity; consistent with `elliptic_k(1) = +∞`).
    /// The amplitude `φ` itself belongs to [`Self::elliptic_f_phi`], not here.
    ///
    /// # Precision
    ///
    /// - Algorithm: Carlson `R_F`; cap `max(128, p/4 + 64)` steps.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn elliptic_f(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let m_c = m.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || m.inexact_accumulated(),
            p,
            rm,
            |q, r| this.elliptic_f_word(&m_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::elliptic_f`].
    fn elliptic_f_word(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match (&self.inner, &m.inner) {
            (Flavor::Value(x), Flavor::Value(mv)) => {
                if let Ok(one) = ExactNumNumber::from_word(1, WORD_BIT_SIZE) {
                    if mv.cmp(&one) == 0 && x.abs_cmp(&one) == 0 {
                        return if x.is_negative() { INF_NEG } else { INF_POS };
                    }
                }
                Self::result_to_ext(x.elliptic_f(mv, p, rm, cc), x.is_zero(), true)
            }
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
            _ => NAN,
        }
    }
    /// Incomplete elliptic integral of the second kind, `E(φ|m)`, with
    /// `self = sin φ` and `m = k²`.
    ///
    /// Same argument convention as [`Self::elliptic_f`]: `self` is `sin φ`,
    /// not the amplitude `φ`, and `|self| > 1` is `NaN(InvalidArgument)`.
    /// `m` is the parameter `k²`, not the modulus `k`. The function is odd in
    /// `self`. The amplitude `φ` itself belongs to [`Self::elliptic_e_phi`],
    /// not here.
    ///
    /// # Precision
    ///
    /// - Algorithm: Carlson `R_F` / `R_D`; cap `max(128, p/4 + 64)` steps.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn elliptic_e(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let m_c = m.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || m.inexact_accumulated(),
            p,
            rm,
            |q, r| this.elliptic_e_word(&m_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::elliptic_e`].
    fn elliptic_e_word(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
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
    /// - Algorithm: Carlson `R_J`; cap `max(128, p/4 + 64)` steps.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn elliptic_pi_complete(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        let this = self.exact_copy_ext();
        let m_c = m.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || m.inexact_accumulated(),
            p,
            rm,
            |q, r| this.elliptic_pi_complete_word(&m_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::elliptic_pi_complete`].
    fn elliptic_pi_complete_word(
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
    /// Incomplete elliptic integral of the third kind, `Π(n; φ|m)`, with
    /// `self = n`, `x = sin φ`, and `m = k²`.
    ///
    /// `x` is the sine of the amplitude, not the amplitude `φ`, same as
    /// [`Self::elliptic_f`]. `|x| > 1` is `NaN(InvalidArgument)`. `m` is the
    /// parameter `k²`, not the modulus `k`. `self` is the characteristic `n`,
    /// not an angle. The amplitude `φ` itself belongs to
    /// [`Self::elliptic_pi_phi`], not here.
    ///
    /// # Precision
    ///
    /// - Algorithm: Carlson `R_J`; cap `max(128, p/4 + 64)` steps.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn elliptic_pi(
        &self,
        x: &Self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        let this = self.exact_copy_ext();
        let x_c = x.exact_copy_ext();
        let m_c = m.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || x.inexact_accumulated() || m.inexact_accumulated(),
            p,
            rm,
            |q, r| this.elliptic_pi_word(&x_c, &m_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::elliptic_pi`].
    fn elliptic_pi_word(
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
    /// Incomplete elliptic integral of the first kind, `F(φ|m)`, with
    /// `self = φ` the amplitude and `m = k²`.
    ///
    /// This takes the angle `φ`. [`Self::elliptic_f`] takes `sin φ`, not the
    /// angle. An angle passed to `elliptic_f` is not `F(φ|m)`.
    ///
    /// Any real `φ` is reduced as `φ = r + kπ` with `|r| ≤ π/2`, at a working
    /// precision that grows with the exponent of `φ`, and
    /// `F(φ|m) = F(r|m) + 2k K(m)`. The function is odd in `φ`.
    /// `φ = ±π/2`, when that equality holds at the requested precision, is
    /// `±K(m)`.
    ///
    /// `m > 1` is `NaN(InvalidArgument)` when the integrand is not real on the
    /// path from `0` to `φ` (it reaches `|sin t| > 1/√m`, including every
    /// `|φ| ≥ π/2`). `m = 1` and `|φ| ≥ π/2` is signed infinity, the same
    /// singularity as `F(±1|1)` and `K(1)`. A `NaN` argument keeps its error
    /// tag. An infinity is `NaN(InvalidArgument)`.
    ///
    /// # Precision
    ///
    /// - Algorithm: argument reduction, then Carlson `R_F` and `K(m)`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn elliptic_f_phi(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let m_c = m.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || m.inexact_accumulated(),
            p,
            rm,
            |q, r| this.elliptic_f_phi_word(&m_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::elliptic_f_phi`].
    fn elliptic_f_phi_word(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match (&self.inner, &m.inner) {
            (Flavor::Value(phi), Flavor::Value(mv)) => {
                Self::result_to_ext(phi.elliptic_f_phi(mv, p, rm, cc), phi.is_zero(), true)
            }
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
            _ => Self::nan(Some(Error::InvalidArgument)),
        }
    }
    /// Incomplete elliptic integral of the second kind, `E(φ|m)`, with
    /// `self = φ` the amplitude and `m = k²`.
    ///
    /// This takes the angle `φ`. [`Self::elliptic_e`] takes `sin φ`, not the
    /// angle. Reduction and the oddness in `φ` match [`Self::elliptic_f_phi`]:
    /// `E(φ|m) = E(r|m) + 2k E(m)`. `φ = ±π/2` is `±E(m)`.
    ///
    /// `m > 1` is `NaN(InvalidArgument)` when the integrand is not real on the
    /// path. `m = 1` stays real for every real `φ` (`E(φ|1) = ∫|cos|`). A
    /// `NaN` argument keeps its error tag. An infinity is
    /// `NaN(InvalidArgument)`.
    ///
    /// # Precision
    ///
    /// - Algorithm: argument reduction, then Carlson `R_F` / `R_D` and `E(m)`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn elliptic_e_phi(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let m_c = m.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || m.inexact_accumulated(),
            p,
            rm,
            |q, r| this.elliptic_e_phi_word(&m_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::elliptic_e_phi`].
    fn elliptic_e_phi_word(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match (&self.inner, &m.inner) {
            (Flavor::Value(phi), Flavor::Value(mv)) => {
                Self::result_to_ext(phi.elliptic_e_phi(mv, p, rm, cc), phi.is_zero(), true)
            }
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
            _ => Self::nan(Some(Error::InvalidArgument)),
        }
    }
    /// Incomplete elliptic integral of the third kind, `Π(n; φ|m)`, with
    /// `self = n`, `phi` the amplitude, and `m = k²`.
    ///
    /// `phi` is the angle. [`Self::elliptic_pi`] takes `sin φ` in that slot,
    /// not the angle. `self` is the characteristic `n`, not an angle.
    /// Reduction matches [`Self::elliptic_f_phi`]:
    /// `Π(n; φ|m) = Π(n; r|m) + 2k Π(n|m)`, and the result is odd in `φ`.
    /// `φ = ±π/2` is `±Π(n|m)` when that complete integral is real.
    ///
    /// `n = 1` is `NaN(InvalidArgument)` at the poles where `sin² t = 1`
    /// (`φ = π/2 + kπ`). `n > 1` is the same error once the path reaches
    /// `sin² t = 1/n`. `n < 1` on either side of those poles is the real
    /// value, including negative `n`. `m > 1` is `NaN(InvalidArgument)` when
    /// the integrand is not real. A `NaN` argument keeps its error tag. An
    /// infinity is `NaN(InvalidArgument)`.
    ///
    /// # Precision
    ///
    /// - Algorithm: argument reduction, then Carlson `R_J` and `Π(n|m)`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn elliptic_pi_phi(
        &self,
        phi: &Self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        let this = self.exact_copy_ext();
        let phi_c = phi.exact_copy_ext();
        let m_c = m.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || phi.inexact_accumulated() || m.inexact_accumulated(),
            p,
            rm,
            |q, r| this.elliptic_pi_phi_word(&phi_c, &m_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::elliptic_pi_phi`].
    fn elliptic_pi_phi_word(
        &self,
        phi: &Self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        match (&self.inner, &phi.inner, &m.inner) {
            (Flavor::Value(n), Flavor::Value(ph), Flavor::Value(mv)) => {
                Self::result_to_ext(n.elliptic_pi_phi(ph, mv, p, rm, cc), ph.is_zero(), true)
            }
            (Flavor::NaN(err), _, _) | (_, Flavor::NaN(err), _) | (_, _, Flavor::NaN(err)) => {
                Self::nan(*err)
            }
            _ => Self::nan(Some(Error::InvalidArgument)),
        }
    }
    /// Jacobi amplitude `am(self | m)`. Parameter `m = k² ∈ [0, 1]`.
    ///
    /// # Precision
    ///
    /// - Algorithm: AGM / descending Landen; cap `JACOBI_AGM_MAX = 128`. \(m=0\) is trig; \(m=1\) is hyperbolic.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no (identity golds; GNU MPFR has no Jacobi `sn`).
    ///
    /// Accuracy: correctly-rounded
    pub fn jacobi_am(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let m_c = m.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || m.inexact_accumulated(),
            p,
            rm,
            |q, r| this.jacobi_am_word(&m_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::jacobi_am`].
    fn jacobi_am_word(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match (&self.inner, &m.inner) {
            (Flavor::Value(u), Flavor::Value(mv)) => {
                Self::result_to_ext(u.jacobi_am(mv, p, rm, cc), u.is_zero(), true)
            }
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
            _ => NAN,
        }
    }
    /// `sn(self | m)`. Parameter `m = k² ∈ [0, 1]`.
    /// The inverse on `(-K, K)` is [`Self::elliptic_f`]: `F(sn(u|m)|m) = u`.
    ///
    /// # Precision
    ///
    /// - Algorithm: AGM amplitude; cap `JACOBI_AGM_MAX = 128`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn jacobi_sn(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let m_c = m.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || m.inexact_accumulated(),
            p,
            rm,
            |q, r| this.jacobi_sn_word(&m_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::jacobi_sn`].
    fn jacobi_sn_word(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match (&self.inner, &m.inner) {
            (Flavor::Value(u), Flavor::Value(mv)) => {
                Self::result_to_ext(u.jacobi_sn(mv, p, rm, cc), u.is_zero(), true)
            }
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
            _ => NAN,
        }
    }
    /// `cn(self | m)`. Parameter `m = k² ∈ [0, 1]`.
    ///
    /// # Precision
    ///
    /// - Algorithm: AGM amplitude; cap `JACOBI_AGM_MAX = 128`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn jacobi_cn(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let m_c = m.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || m.inexact_accumulated(),
            p,
            rm,
            |q, r| this.jacobi_cn_word(&m_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::jacobi_cn`].
    fn jacobi_cn_word(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match (&self.inner, &m.inner) {
            (Flavor::Value(u), Flavor::Value(mv)) => {
                Self::result_to_ext(u.jacobi_cn(mv, p, rm, cc), u.is_zero(), true)
            }
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
            _ => NAN,
        }
    }
    /// `dn(self | m)`. Parameter `m = k² ∈ [0, 1]`.
    ///
    /// # Precision
    ///
    /// - Algorithm: AGM amplitude; cap `JACOBI_AGM_MAX = 128`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn jacobi_dn(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let m_c = m.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || m.inexact_accumulated(),
            p,
            rm,
            |q, r| this.jacobi_dn_word(&m_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::jacobi_dn`].
    fn jacobi_dn_word(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match (&self.inner, &m.inner) {
            (Flavor::Value(u), Flavor::Value(mv)) => {
                Self::result_to_ext(u.jacobi_dn(mv, p, rm, cc), false, true)
            }
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
            _ => NAN,
        }
    }
    /// `cd(self | m) = cn / dn`. Parameter `m = k² ∈ [0, 1]`.
    ///
    /// # Precision
    ///
    /// - Algorithm: AGM amplitude; cap `JACOBI_AGM_MAX = 128`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn jacobi_cd(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let m_c = m.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || m.inexact_accumulated(),
            p,
            rm,
            |q, r| this.jacobi_cd_word(&m_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::jacobi_cd`].
    fn jacobi_cd_word(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match (&self.inner, &m.inner) {
            (Flavor::Value(u), Flavor::Value(mv)) => {
                Self::result_to_ext(u.jacobi_cd(mv, p, rm, cc), false, true)
            }
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
            _ => NAN,
        }
    }
    /// `ns(self | m) = 1/sn`.
    ///
    /// # Precision
    ///
    /// - Algorithm: AGM amplitude; cap `JACOBI_AGM_MAX = 128`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn jacobi_ns(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let m_c = m.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || m.inexact_accumulated(),
            p,
            rm,
            |q, r| this.jacobi_ns_word(&m_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::jacobi_ns`].
    fn jacobi_ns_word(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match (&self.inner, &m.inner) {
            (Flavor::Value(u), Flavor::Value(mv)) => {
                Self::result_to_ext(u.jacobi_ns(mv, p, rm, cc), false, true)
            }
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
            _ => NAN,
        }
    }
    /// `nc(self | m) = 1/cn`.
    ///
    /// # Precision
    ///
    /// - Algorithm: AGM amplitude; cap `JACOBI_AGM_MAX = 128`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn jacobi_nc(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let m_c = m.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || m.inexact_accumulated(),
            p,
            rm,
            |q, r| this.jacobi_nc_word(&m_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::jacobi_nc`].
    fn jacobi_nc_word(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match (&self.inner, &m.inner) {
            (Flavor::Value(u), Flavor::Value(mv)) => {
                Self::result_to_ext(u.jacobi_nc(mv, p, rm, cc), false, true)
            }
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
            _ => NAN,
        }
    }
    /// `nd(self | m) = 1/dn`.
    ///
    /// # Precision
    ///
    /// - Algorithm: AGM amplitude; cap `JACOBI_AGM_MAX = 128`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn jacobi_nd(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let m_c = m.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || m.inexact_accumulated(),
            p,
            rm,
            |q, r| this.jacobi_nd_word(&m_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::jacobi_nd`].
    fn jacobi_nd_word(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match (&self.inner, &m.inner) {
            (Flavor::Value(u), Flavor::Value(mv)) => {
                Self::result_to_ext(u.jacobi_nd(mv, p, rm, cc), false, true)
            }
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
            _ => NAN,
        }
    }
    /// `sc(self | m) = sn/cn`.
    ///
    /// # Precision
    ///
    /// - Algorithm: AGM amplitude; cap `JACOBI_AGM_MAX = 128`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn jacobi_sc(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let m_c = m.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || m.inexact_accumulated(),
            p,
            rm,
            |q, r| this.jacobi_sc_word(&m_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::jacobi_sc`].
    fn jacobi_sc_word(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match (&self.inner, &m.inner) {
            (Flavor::Value(u), Flavor::Value(mv)) => {
                Self::result_to_ext(u.jacobi_sc(mv, p, rm, cc), false, true)
            }
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
            _ => NAN,
        }
    }
    /// `sd(self | m) = sn/dn`.
    ///
    /// # Precision
    ///
    /// - Algorithm: AGM amplitude; cap `JACOBI_AGM_MAX = 128`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn jacobi_sd(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let m_c = m.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || m.inexact_accumulated(),
            p,
            rm,
            |q, r| this.jacobi_sd_word(&m_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::jacobi_sd`].
    fn jacobi_sd_word(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match (&self.inner, &m.inner) {
            (Flavor::Value(u), Flavor::Value(mv)) => {
                Self::result_to_ext(u.jacobi_sd(mv, p, rm, cc), false, true)
            }
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
            _ => NAN,
        }
    }
    /// `cs(self | m) = cn/sn`.
    ///
    /// # Precision
    ///
    /// - Algorithm: AGM amplitude; cap `JACOBI_AGM_MAX = 128`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn jacobi_cs(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let m_c = m.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || m.inexact_accumulated(),
            p,
            rm,
            |q, r| this.jacobi_cs_word(&m_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::jacobi_cs`].
    fn jacobi_cs_word(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match (&self.inner, &m.inner) {
            (Flavor::Value(u), Flavor::Value(mv)) => {
                Self::result_to_ext(u.jacobi_cs(mv, p, rm, cc), false, true)
            }
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
            _ => NAN,
        }
    }
    /// `ds(self | m) = dn/sn`.
    ///
    /// # Precision
    ///
    /// - Algorithm: AGM amplitude; cap `JACOBI_AGM_MAX = 128`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn jacobi_ds(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let m_c = m.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || m.inexact_accumulated(),
            p,
            rm,
            |q, r| this.jacobi_ds_word(&m_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::jacobi_ds`].
    fn jacobi_ds_word(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match (&self.inner, &m.inner) {
            (Flavor::Value(u), Flavor::Value(mv)) => {
                Self::result_to_ext(u.jacobi_ds(mv, p, rm, cc), false, true)
            }
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
            _ => NAN,
        }
    }
    /// `dc(self | m) = dn/cn`.
    ///
    /// # Precision
    ///
    /// - Algorithm: AGM amplitude; cap `JACOBI_AGM_MAX = 128`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn jacobi_dc(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let m_c = m.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || m.inexact_accumulated(),
            p,
            rm,
            |q, r| this.jacobi_dc_word(&m_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::jacobi_dc`].
    fn jacobi_dc_word(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match (&self.inner, &m.inner) {
            (Flavor::Value(u), Flavor::Value(mv)) => {
                Self::result_to_ext(u.jacobi_dc(mv, p, rm, cc), false, true)
            }
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
            _ => NAN,
        }
    }
    /// Legendre \(P_n(\mathrm{self})\) for integer `n`.
    ///
    /// # Precision
    ///
    /// - Algorithm: three-term recurrence at `p + O(n)` bits. Cap `ORTHOPOLY_N_MAX`.
    /// - Bound: working-precision recurrence (not a Ziv leaf).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: working-precision, 32 bits
    pub fn legendre_p(&self, n: u32, p: usize, rm: RoundingMode) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.legendre_p_word(n, q, r)
        })
    }

    /// Word-precision implementation of [`Self::legendre_p`].
    fn legendre_p_word(&self, n: u32, p: usize, rm: RoundingMode) -> Self {
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
    ///
    /// Accuracy: working-precision, 32 bits
    pub fn assoc_legendre_p(&self, n: u32, m: i32, p: usize, rm: RoundingMode) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.assoc_legendre_p_word(n, m, q, r)
        })
    }

    /// Word-precision implementation of [`Self::assoc_legendre_p`].
    fn assoc_legendre_p_word(&self, n: u32, m: i32, p: usize, rm: RoundingMode) -> Self {
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
    /// - Domain: real `z ≤ 1` (`z = 1` needs `c − a − b > 0`); `z > 1` (branch cut) and non-positive
    ///   integer `c` without an earlier terminating numerator parameter → NaN.
    /// - Algorithm: terminating series (any `z`); series for `|z| ≤ 1/2`; Pfaff `z/(z−1)` for `z < −1/2`;
    ///   `1 − z` connection (DLMF 15.8.4, logarithmic form A&S 15.3.10–11 when `c − a − b` is an
    ///   integer) for `1/2 < z < 1`; Gauss at `z = 1`. Cancellation is measured and the evaluation
    ///   repeated with more guard bits. Series cap `HYPERGEOM_TERM_MAX`.
    /// - Bound: Ziv correct-rounding (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn hypergeom_2f1(
        &self,
        b: &Self,
        c: &Self,
        z: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        let this = self.exact_copy_ext();
        let b_c = b.exact_copy_ext();
        let c_c = c.exact_copy_ext();
        let z_c = z.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated()
                || b.inexact_accumulated()
                || c.inexact_accumulated()
                || z.inexact_accumulated(),
            p,
            rm,
            |q, r| this.hypergeom_2f1_word(&b_c, &c_c, &z_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::hypergeom_2f1`].
    fn hypergeom_2f1_word(
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
    ///
    /// Accuracy: correctly-rounded
    pub fn betainc(&self, b: &Self, x: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let b_c = b.exact_copy_ext();
        let x_c = x.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || b.inexact_accumulated() || x.inexact_accumulated(),
            p,
            rm,
            |q, r| this.betainc_word(&b_c, &x_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::betainc`].
    fn betainc_word(
        &self,
        b: &Self,
        x: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
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
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.
Accuracy: correctly-rounded",
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
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.
Accuracy: correctly-rounded",
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
        Precision is rounded upwards to the word size. The function returns NaN if the precision `p` is incorrect.
Accuracy: correctly-rounded",
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

//
// ordering traits
//

impl PartialEq for ExactNum {
    fn eq(&self, other: &Self) -> bool {
        let cmp_result = ExactNum::compare(self, other);
        matches!(cmp_result, Some(0))
    }
}

impl<'a> PartialEq<&'a ExactNum> for ExactNum {
    fn eq(&self, other: &&'a ExactNum) -> bool {
        let cmp_result = ExactNum::compare(self, other);
        matches!(cmp_result, Some(0))
    }
}

impl Eq for ExactNum {}

impl PartialOrd for ExactNum {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let cmp_result = ExactNum::compare(self, other);
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
        let cmp_result = ExactNum::compare(self, other);
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

    /// Returns the parsed number.
    /// Malformed text is `Ok(NaN(InvalidArgument))`.
    /// Other failures, such as exponent overflow, are `Err`.
    /// The implementation is not available in no_std environment.
    fn from_str(src: &str) -> Result<ExactNum, Self::Err> {
        let bf = crate::common::consts::with_tenpowers(|cc| {
            Ok(ExactNum::parse(
                src,
                Radix::Dec,
                usize::MAX,
                RoundingMode::ToEven,
                cc,
            ))
        })?;

        if bf.is_nan() {
            match bf.err() {
                Some(Error::InvalidArgument) | None => {}
                Some(err) => return Err(err),
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
                crate::common::consts::with_tenpowers_fmt(|cc| {
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
        crate::common::consts::with_tenpowers_fmt(|cc| {
            let mut s = String::new();
            self.write_str(&mut s, Radix::Dec, RoundingMode::ToEven, cc)?;
            f.write_str(&s.replace('e', "E"))
        })
    }
}
#[cfg(feature = "std")]
impl core::fmt::LowerHex for ExactNum {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> Result<(), core::fmt::Error> {
        crate::common::consts::with_tenpowers_fmt(|cc| {
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
    use crate::Consts;
    use crate::Error;

    static ONE: ExactNum = ExactNum {
        inner: super::Flavor::Value(crate::num::ExactNumNumber::small_const(1)),
    };
    static TWO: ExactNum = ExactNum {
        inner: super::Flavor::Value(crate::num::ExactNumNumber::small_const(2)),
    };
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
            assert!(d1.compare(&n1) == Some(0));

            let i = crate::common::test_rng::random::<u64>();
            let d1 = ExactNum::from_u64(i, rand_p());
            let n1 = ExactNum::parse(&format!("{}", i), Radix::Dec, rand_p(), rm, &mut cc);
            assert!(d1.compare(&n1) == Some(0));

            let i = crate::common::test_rng::random::<i128>();
            let d1 = ExactNum::from_i128(i, rand_p());
            let n1 = ExactNum::parse(&format!("{}", i), Radix::Dec, rand_p(), rm, &mut cc);
            assert!(d1.compare(&n1) == Some(0));

            let i = crate::common::test_rng::random::<u128>();
            let d1 = ExactNum::from_u128(i, rand_p());
            let n1 = ExactNum::parse(&format!("{}", i), Radix::Dec, rand_p(), rm, &mut cc);
            assert!(d1.compare(&n1) == Some(0));
        }

        assert!(ONE.exponent().is_some());
        assert!(INF_POS.exponent().is_none());
        assert!(INF_NEG.exponent().is_none());
        assert!(NAN.exponent().is_none());

        assert!(ONE.as_raw_parts().is_some());
        assert!(INF_POS.as_raw_parts().is_none());
        assert!(INF_NEG.as_raw_parts().is_none());
        assert!(NAN.as_raw_parts().is_none());

        assert!(ONE.add(&ONE, rand_p(), rm).compare(&TWO) == Some(0));
        assert!(ONE.add(&INF_POS, rand_p(), rm).is_inf_pos());
        assert!(INF_POS.add(&ONE, rand_p(), rm).is_inf_pos());
        assert!(ONE.add(&INF_NEG, rand_p(), rm).is_inf_neg());
        assert!(INF_NEG.add(&ONE, rand_p(), rm).is_inf_neg());
        assert!(INF_POS.add(&INF_POS, rand_p(), rm).is_inf_pos());
        assert!(INF_POS.add(&INF_NEG, rand_p(), rm).is_nan());
        assert!(INF_NEG.add(&INF_NEG, rand_p(), rm).is_inf_neg());
        assert!(INF_NEG.add(&INF_POS, rand_p(), rm).is_nan());

        assert!(ONE.add_full_prec(&ONE).compare(&TWO) == Some(0));
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

        assert!(ONE.mul_full_prec(&ONE).compare(&ONE) == Some(0));
        assert!(ONE.mul_full_prec(&INF_POS).is_inf_pos());
        assert!(INF_POS.mul_full_prec(&ONE).is_inf_pos());
        assert!(ONE.mul_full_prec(&INF_NEG).is_inf_neg());
        assert!(INF_NEG.mul_full_prec(&ONE).is_inf_neg());
        assert!(INF_POS.mul_full_prec(&INF_POS).is_inf_pos());
        assert!(INF_POS.mul_full_prec(&INF_NEG).is_inf_neg());
        assert!(INF_NEG.mul_full_prec(&INF_NEG).is_inf_pos());
        assert!(INF_NEG.mul_full_prec(&INF_POS).is_inf_neg());

        assert!(TWO.sub(&ONE, rand_p(), rm).compare(&ONE) == Some(0));
        assert!(ONE.sub(&INF_POS, rand_p(), rm).is_inf_neg());
        assert!(INF_POS.sub(&ONE, rand_p(), rm).is_inf_pos());
        assert!(ONE.sub(&INF_NEG, rand_p(), rm).is_inf_pos());
        assert!(INF_NEG.sub(&ONE, rand_p(), rm).is_inf_neg());
        assert!(INF_POS.sub(&INF_POS, rand_p(), rm).is_nan());
        assert!(INF_POS.sub(&INF_NEG, rand_p(), rm).is_inf_pos());
        assert!(INF_NEG.sub(&INF_NEG, rand_p(), rm).is_nan());
        assert!(INF_NEG.sub(&INF_POS, rand_p(), rm).is_inf_neg());

        assert!(TWO.mul(&ONE, rand_p(), rm).compare(&TWO) == Some(0));
        assert!(ONE.mul(&INF_POS, rand_p(), rm).is_inf_pos());
        assert!(INF_POS.mul(&ONE, rand_p(), rm).is_inf_pos());
        assert!(ONE.mul(&INF_NEG, rand_p(), rm).is_inf_neg());
        assert!(INF_NEG.mul(&ONE, rand_p(), rm).is_inf_neg());
        assert!(ExactNum::from_i8(-1, DEFAULT_P)
            .mul(&INF_POS, rand_p(), rm)
            .is_inf_neg());
        assert!(ExactNum::from_i8(-1, DEFAULT_P)
            .mul(&INF_NEG, rand_p(), rm)
            .is_inf_pos());
        assert!(INF_POS
            .mul(&ExactNum::from_i8(-1, DEFAULT_P), rand_p(), rm)
            .is_inf_neg());
        assert!(INF_NEG
            .mul(&ExactNum::from_i8(-1, DEFAULT_P), rand_p(), rm)
            .is_inf_pos());
        assert!(INF_POS.mul(&INF_POS, rand_p(), rm).is_inf_pos());
        assert!(INF_POS.mul(&INF_NEG, rand_p(), rm).is_inf_neg());
        assert!(INF_NEG.mul(&INF_NEG, rand_p(), rm).is_inf_pos());
        assert!(INF_NEG.mul(&INF_POS, rand_p(), rm).is_inf_neg());
        assert!(INF_POS.mul(&ExactNum::new(rand_p()), rand_p(), rm).is_nan());
        assert!(INF_NEG.mul(&ExactNum::new(rand_p()), rand_p(), rm).is_nan());
        assert!(ExactNum::new(rand_p()).mul(&INF_POS, rand_p(), rm).is_nan());
        assert!(ExactNum::new(rand_p()).mul(&INF_NEG, rand_p(), rm).is_nan());

        assert!(TWO.div(&TWO, rand_p(), rm).compare(&ONE) == Some(0));
        assert!(TWO.div(&INF_POS, rand_p(), rm).is_zero());
        assert!(INF_POS.div(&TWO, rand_p(), rm).is_inf_pos());
        assert!(TWO.div(&INF_NEG, rand_p(), rm).is_zero());
        assert!(INF_NEG.div(&TWO, rand_p(), rm).is_inf_neg());
        assert!(ExactNum::from_i8(-2, DEFAULT_P)
            .div(&INF_POS, rand_p(), rm)
            .is_zero());
        assert!(ExactNum::from_i8(-2, DEFAULT_P)
            .div(&INF_NEG, rand_p(), rm)
            .is_zero());
        assert!(INF_POS
            .div(&ExactNum::from_i8(-2, DEFAULT_P), rand_p(), rm)
            .is_inf_neg());
        assert!(INF_NEG
            .div(&ExactNum::from_i8(-2, DEFAULT_P), rand_p(), rm)
            .is_inf_pos());
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
        assert!(TWO.rem(&INF_POS).compare(&TWO) == Some(0));
        assert!(INF_POS.rem(&TWO).is_nan());
        assert!(TWO.rem(&INF_NEG).compare(&TWO) == Some(0));
        assert!(INF_NEG.rem(&TWO).is_nan());
        assert!(
            ExactNum::from_i8(-2, DEFAULT_P)
                .rem(&INF_POS)
                .compare(&ExactNum::from_i8(-2, DEFAULT_P))
                == Some(0)
        );
        assert!(
            ExactNum::from_i8(-2, DEFAULT_P)
                .rem(&INF_NEG)
                .compare(&ExactNum::from_i8(-2, DEFAULT_P))
                == Some(0)
        );
        assert!(INF_POS.rem(&ExactNum::from_i8(-2, DEFAULT_P)).is_nan());
        assert!(INF_NEG.rem(&ExactNum::from_i8(-2, DEFAULT_P)).is_nan());
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

        assert!(ONE.compare(&ONE).unwrap() == 0);
        assert!(ONE.compare(&INF_POS).unwrap() < 0);
        assert!(INF_POS.compare(&ONE).unwrap() > 0);
        assert!(INF_POS.compare(&INF_POS).unwrap() == 0);
        assert!(ONE.compare(&INF_NEG).unwrap() > 0);
        assert!(INF_NEG.compare(&ONE).unwrap() < 0);
        assert!(INF_NEG.compare(&INF_NEG).unwrap() == 0);
        assert!(INF_POS.compare(&INF_NEG).unwrap() > 0);
        assert!(INF_NEG.compare(&INF_POS).unwrap() < 0);
        assert!(INF_POS.compare(&INF_POS).unwrap() == 0);
        assert!(ONE.compare(&NAN).is_none());
        assert!(NAN.compare(&ONE).is_none());
        assert!(INF_POS.compare(&NAN).is_none());
        assert!(NAN.compare(&INF_POS).is_none());
        assert!(INF_NEG.compare(&NAN).is_none());
        assert!(NAN.compare(&INF_NEG).is_none());
        assert!(NAN.compare(&NAN).is_none());

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
        // Finite values compare by magnitude, whatever their signs (1.0.4 compared signed values).
        let two = ExactNum::from_u8(2, 64);
        let neg_two = two.neg_via_mul();
        assert!(neg_two.abs_cmp(&ONE).unwrap() > 0);
        assert!(ONE.abs_cmp(&neg_two).unwrap() < 0);
        assert!(neg_two.abs_cmp(&two).unwrap() == 0);
        assert!(ExactNum::from_i8(-1, DEFAULT_P).abs_cmp(&INF_POS).unwrap() < 0);
        // Results are exactly -1/0/1 (1.0.4 returned the mantissa word difference, e.g.
        // 3.cmp(2) = 2^62, so callers testing `== Some(1)` / `== Some(-1)` misfired).
        let three = ExactNum::from_u8(3, 64);
        assert_eq!(three.compare(&two), Some(1));
        assert_eq!(two.compare(&three), Some(-1));
        assert_eq!(three.abs_cmp(&neg_two), Some(1));
        assert_eq!(neg_two.abs_cmp(&three), Some(-1));
        assert_eq!(INF_POS.compare(&INF_NEG), Some(1));
        assert_eq!(INF_NEG.compare(&INF_POS), Some(-1));

        assert!(ONE.is_positive());
        assert!(!ONE.is_negative());

        assert!(ExactNum::from_i8(-1, DEFAULT_P).is_negative());
        assert!(!ExactNum::from_i8(-1, DEFAULT_P).is_positive());
        assert!(!INF_POS.is_negative());
        assert!(INF_POS.is_positive());
        assert!(INF_NEG.is_negative());
        assert!(!INF_NEG.is_positive());
        assert!(!NAN.is_positive());
        assert!(!NAN.is_negative());

        assert!(ONE.pow(&ONE, rand_p(), rm, &mut cc).compare(&ONE) == Some(0));
        assert!(ExactNum::new(DEFAULT_P)
            .pow(&INF_POS, rand_p(), rm, &mut cc)
            .is_zero());
        assert!(ExactNum::new(DEFAULT_P)
            .pow(&INF_NEG, rand_p(), rm, &mut cc)
            .is_zero());
        assert!(ONE.pow(&INF_POS, rand_p(), rm, &mut cc).compare(&ONE) == Some(0));
        assert!(ONE.pow(&INF_NEG, rand_p(), rm, &mut cc).compare(&ONE) == Some(0));
        assert!(TWO.pow(&INF_POS, rand_p(), rm, &mut cc).is_inf_pos());
        assert!(TWO.pow(&INF_NEG, rand_p(), rm, &mut cc).is_inf_neg());
        assert!(INF_POS.pow(&ONE, rand_p(), rm, &mut cc).is_inf_pos());
        assert!(INF_NEG.pow(&ONE, rand_p(), rm, &mut cc).is_inf_neg());
        assert!(INF_NEG.pow(&TWO, rand_p(), rm, &mut cc).is_inf_pos());
        assert!(INF_POS
            .pow(&ExactNum::from_i8(-1, DEFAULT_P), rand_p(), rm, &mut cc)
            .is_zero());
        assert!(INF_NEG
            .pow(&ExactNum::from_i8(-1, DEFAULT_P), rand_p(), rm, &mut cc)
            .is_zero());
        assert!(
            INF_POS
                .pow(&ExactNum::new(DEFAULT_P), rand_p(), rm, &mut cc)
                .compare(&ONE)
                == Some(0)
        );
        assert!(
            INF_NEG
                .pow(&ExactNum::new(DEFAULT_P), rand_p(), rm, &mut cc)
                .compare(&ONE)
                == Some(0)
        );
        assert!(INF_POS.pow(&INF_POS, rand_p(), rm, &mut cc).is_inf_pos());
        assert!(INF_NEG.pow(&INF_POS, rand_p(), rm, &mut cc).is_inf_pos());
        assert!(INF_POS.pow(&INF_NEG, rand_p(), rm, &mut cc).is_zero());
        assert!(INF_NEG.pow(&INF_NEG, rand_p(), rm, &mut cc).is_zero());

        let half = ONE.div(&TWO, rand_p(), rm);
        assert!(TWO.log(&TWO, rand_p(), rm, &mut cc).compare(&ONE) == Some(0));
        assert!(TWO.log(&INF_POS, rand_p(), rm, &mut cc).is_zero());
        assert!(TWO.log(&INF_NEG, rand_p(), rm, &mut cc).is_nan());
        assert!(INF_POS.log(&TWO, rand_p(), rm, &mut cc).is_inf_pos());
        assert!(INF_NEG.log(&TWO, rand_p(), rm, &mut cc).is_nan());
        assert!(half.log(&half, rand_p(), rm, &mut cc).compare(&ONE) == Some(0));
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
        assert!(INF_POS.powi(0, rand_p(), rm).compare(&ONE) == Some(0));
        assert!(INF_NEG.powi(0, rand_p(), rm).compare(&ONE) == Some(0));

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
        assert!(INF_NEG.atan(p, rm, &mut cc).compare(&half_pi.neg_via_mul()) == Some(0));
        assert!(INF_POS.atan(p, rm, &mut cc).compare(&half_pi) == Some(0));
        assert!(NAN.atan(rand_p(), rm, &mut cc).is_nan());

        assert!(INF_NEG.sinh(rand_p(), rm, &mut cc).is_inf_neg());
        assert!(INF_POS.sinh(rand_p(), rm, &mut cc).is_inf_pos());
        assert!(NAN.sinh(rand_p(), rm, &mut cc).is_nan());

        assert!(INF_NEG.cosh(rand_p(), rm, &mut cc).is_inf_pos());
        assert!(INF_POS.cosh(rand_p(), rm, &mut cc).is_inf_pos());
        assert!(NAN.cosh(rand_p(), rm, &mut cc).is_nan());

        assert!(
            INF_NEG
                .tanh(rand_p(), rm, &mut cc)
                .compare(&ExactNum::from_i8(-1, DEFAULT_P))
                == Some(0)
        );
        assert!(INF_POS.tanh(rand_p(), rm, &mut cc).compare(&ONE) == Some(0));
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

        assert!(TWO.signum().compare(&ONE) == Some(0));
        assert!(
            ExactNum::from_i8(-2, DEFAULT_P)
                .signum()
                .compare(&ExactNum::from_i8(-1, DEFAULT_P))
                == Some(0)
        );
        assert!(INF_POS.signum().compare(&ONE) == Some(0));
        assert!(INF_NEG.signum().compare(&ExactNum::from_i8(-1, DEFAULT_P)) == Some(0));
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
        // DEFAULT_P = 128 bits is INLINE_WORDS limbs on every target (2 on 64-bit, 4 on 32-bit).
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

        assert!(INF_POS.min(&ONE).compare(&ONE) == Some(0));
        assert!(INF_NEG.min(&ONE).is_inf_neg());
        assert!(NAN.min(&ONE).is_nan());
        assert!(ONE.min(&INF_POS).compare(&ONE) == Some(0));
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
        assert!(INF_NEG.max(&ONE).compare(&ONE) == Some(0));
        assert!(NAN.max(&ONE).is_nan());
        assert!(ONE.max(&INF_POS).is_inf_pos());
        assert!(ONE.max(&INF_NEG).compare(&ONE) == Some(0));
        assert!(ONE.max(&NAN).is_nan());
        assert!(NAN.max(&INF_POS).is_nan());
        assert!(NAN.max(&INF_NEG).is_nan());
        assert!(NAN.max(&NAN).is_nan());
        assert!(INF_NEG.max(&INF_POS).is_inf_pos());
        assert!(INF_POS.max(&INF_NEG).is_inf_pos());
        assert!(INF_POS.max(&INF_POS).is_inf_pos());
        assert!(INF_NEG.max(&INF_NEG).is_inf_neg());

        assert!(
            ONE.clamp(&ExactNum::from_i8(-1, DEFAULT_P), &TWO)
                .compare(&ONE)
                == Some(0)
        );
        assert!(ONE.clamp(&TWO, &ONE).is_nan());
        assert!(ONE.clamp(&INF_POS, &ONE).is_nan());
        assert!(ONE.clamp(&TWO, &INF_NEG).is_nan());
        assert!(
            ExactNum::from_i8(-1, DEFAULT_P)
                .clamp(&ONE, &TWO)
                .compare(&ONE)
                == Some(0)
        );
        assert!(
            TWO.clamp(&ExactNum::from_i8(-1, DEFAULT_P), &ONE)
                .compare(&ONE)
                == Some(0)
        );
        assert!(INF_POS.clamp(&ONE, &TWO).compare(&TWO) == Some(0));
        assert!(INF_POS.clamp(&ONE, &INF_POS).is_inf_pos());
        assert!(INF_POS.clamp(&INF_NEG, &ONE).compare(&ONE) == Some(0));
        assert!(INF_POS.clamp(&NAN, &INF_POS).is_nan());
        assert!(INF_POS.clamp(&ONE, &NAN).is_nan());
        assert!(INF_POS.clamp(&NAN, &NAN).is_nan());
        assert!(INF_NEG.clamp(&ONE, &TWO).compare(&ONE) == Some(0));
        assert!(INF_NEG.clamp(&ONE, &INF_POS).compare(&ONE) == Some(0));
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

        let bad = ExactNum::from_str("abc").unwrap();
        assert!(bad.is_nan());
        assert_eq!(bad.err(), Some(Error::InvalidArgument));
    }

    #[test]
    pub fn test_ops() {
        use core::ops::Not;
        let mut cc = Consts::new().unwrap();

        let d1 = -ExactNum::from_u8(2, DEFAULT_P);
        assert!(d1.is_negative());

        let p = DEFAULT_P;
        let rm = RoundingMode::ToEven;
        let two = ExactNum::from_u8(2, p);
        let eighth = two.powsi(-3, p, rm);
        let expected = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(8, p), p, rm);
        assert_eq!(eighth.compare(&expected), Some(0));
        assert_eq!(
            two.powsi(3, p, rm).compare(&ExactNum::from_u8(8, p)),
            Some(0)
        );
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
        assert!((pos < neg).not());
        assert!((neg > pos).not());
        assert!(INF_NEG < neg);
        assert!(INF_NEG < pos);
        assert!(INF_NEG < INF_POS);
        assert!((INF_NEG > neg).not());
        assert!((INF_NEG > pos).not());
        assert!((INF_NEG > INF_POS).not());
        assert!(INF_POS > neg);
        assert!(INF_POS > pos);
        assert!(INF_POS > INF_NEG);
        assert!((INF_POS < neg).not());
        assert!((INF_POS < pos).not());
        assert!((INF_POS < INF_NEG).not());
        assert!((INF_POS > INF_POS).not());
        assert!((INF_POS < INF_POS).not());
        assert!((INF_NEG > INF_NEG).not());
        assert!((INF_NEG < INF_NEG).not());
        assert!((INF_POS > NAN).not());
        assert!((INF_POS < NAN).not());
        assert!((INF_NEG > NAN).not());
        assert!((INF_NEG < NAN).not());
        assert!((NAN > INF_POS).not());
        assert!((NAN < INF_POS).not());
        assert!((NAN > INF_NEG).not());
        assert!((NAN < INF_NEG).not());
        assert!((NAN > NAN).not());
        assert!((NAN < NAN).not());
        assert!((neg > NAN).not());
        assert!((neg < NAN).not());
        assert!((pos > NAN).not());
        assert!((pos < NAN).not());
        assert!((NAN > neg).not());
        assert!((NAN < neg).not());
        assert!((NAN > pos).not());
        assert!((NAN < pos).not());

        assert!((NAN == NAN).not());
        assert!((NAN == INF_POS).not());
        assert!((NAN == INF_NEG).not());
        assert!((INF_POS == NAN).not());
        assert!((INF_NEG == NAN).not());
        assert!((INF_NEG == INF_POS).not());
        assert!((INF_POS == INF_NEG).not());
        assert!((INF_POS == neg).not());
        assert!((INF_POS == pos).not());
        assert!((INF_NEG == neg).not());
        assert!((INF_NEG == pos).not());
        assert!((neg == INF_POS).not());
        assert!((pos == INF_POS).not());
        assert!((neg == INF_NEG).not());
        assert!((pos == INF_NEG).not());
        assert!((pos == neg).not());
        assert!((neg == pos).not());
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
        assert_eq!(s.compare(&ExactNum::from_word(8, p)), Some(0));
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
        assert_eq!(rec.compare(&ExactNum::from_u8(3, p)), Some(0));

        let (ph, pl) = two.two_product(&three, p, rm);
        let pr = ph.add(&pl, p, rm);
        assert_eq!(pr.compare(&ExactNum::from_u8(6, p)), Some(0));

        let sum = ExactNum::fused_sum(&[one.clone(), two.clone(), three.clone()], p, rm);
        assert_eq!(sum.compare(&ExactNum::from_u8(6, p)), Some(0));

        let dot = ExactNum::fused_dot(
            &[one.clone(), two.clone()],
            &[three.clone(), one.clone()],
            p,
            rm,
        );
        assert_eq!(dot.compare(&ExactNum::from_u8(5, p)), Some(0));

        // 1 + 2x + 3x² at x = 2 → 17
        let pv = ExactNum::polyval(&[one, two, three], &ExactNum::from_u8(2, p), p, rm);
        assert_eq!(pv.compare(&ExactNum::from_u8(17, p)), Some(0));
    }

    #[test]
    fn test_jacobi_sn_public() {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new().unwrap();
        let one = ExactNum::from_u8(1, p);
        let zero = ExactNum::from_u8(0, p);
        let half = one.div(&ExactNum::from_u8(2, p), p, rm);
        let sn0 = zero.jacobi_sn(&half, p, rm, &mut cc);
        assert!(sn0.is_zero(), "sn(0)");
        let cn0 = zero.jacobi_cn(&half, p, rm, &mut cc);
        assert_eq!(cn0.compare(&one), Some(0), "cn(0)");
        let sn_m0 = one.jacobi_sn(&zero, p, rm, &mut cc);
        let sin1 = one.sin(p, rm, &mut cc);
        let d = sn_m0.sub(&sin1, p, rm).abs();
        assert!(d.is_zero() || d.exponent().unwrap() < -80, "sn(1|0)=sin 1");
        assert!(one
            .jacobi_sn(&ExactNum::from_i8(2, p), p, rm, &mut cc)
            .is_nan());
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

/// Characterization of the NaN / Inf argument paths: the result of an operation with a NaN operand is a NaN
/// carrying that operand's error (first NaN operand in argument order wins), Inf operands keep their sign.
#[cfg(test)]
mod special_arg_tests {
    use crate::defs::{RoundingMode, DEFAULT_P};
    use crate::{Consts, Error, ExactNum, INF_NEG, INF_POS};
    use core::ops::Neg;

    const P: usize = DEFAULT_P;
    const RM: RoundingMode = RoundingMode::ToEven;

    fn nans() -> [ExactNum; 2] {
        [ExactNum::nan(Some(Error::InvalidArgument)), ExactNum::nan(None)]
    }

    fn same_nan(r: &ExactNum, want: &ExactNum) {
        assert!(r.is_nan(), "expected NaN, got {r:?}");
        assert_eq!(r.err(), want.err());
    }

    fn one() -> ExactNum {
        ExactNum::from_u8(1, P)
    }

    #[test]
    fn nan_operand_wins_in_ternary_and_binary_ops() {
        let mut cc = match Consts::new() {
            Ok(c) => c,
            Err(e) => panic!("consts: {e:?}"),
        };
        for n in nans() {
            let v = one();
            for (a, b, c) in [(&n, &v, &v), (&v, &n, &v), (&v, &v, &n)] {
                same_nan(&a.fma(b, c, P, RM), &n);
                same_nan(&a.mul_add(b, c, P, RM), &n);
            }
            for (a, b) in [(&n, &v), (&v, &n)] {
                same_nan(&a.atan2(b, P, RM, &mut cc), &n);
                same_nan(&a.next_after(b, P, RM), &n);
            }
            same_nan(&n.copysign(&v, P, RM), &n);
            same_nan(&n.hypot(&v, P, RM), &n);
            same_nan(&v.hypot(&n, P, RM), &n);
            let (hi, lo) = n.two_sum(&v, P, RM);
            same_nan(&hi, &n);
            same_nan(&lo, &n);
            let (hi, lo) = v.two_sum(&n, P, RM);
            same_nan(&hi, &n);
            same_nan(&lo, &n);
            let (hi, lo) = n.two_product(&v, P, RM);
            same_nan(&hi, &n);
            same_nan(&lo, &n);
            let (hi, lo) = v.two_product(&n, P, RM);
            same_nan(&hi, &n);
            same_nan(&lo, &n);
        }
    }

    #[test]
    fn min_max_clamp_borrow_an_operand() {
        let (lo, hi, mid) = (one(), ExactNum::from_u8(3, P), ExactNum::from_u8(2, P));
        assert!(core::ptr::eq(lo.max(&hi), &hi));
        assert!(core::ptr::eq(hi.max(&lo), &hi));
        assert!(core::ptr::eq(lo.min(&hi), &lo));
        assert!(core::ptr::eq(hi.min(&lo), &lo));
        assert!(core::ptr::eq(mid.clamp(&lo, &hi), &mid));
        assert!(core::ptr::eq(lo.clamp(&mid, &hi), &mid));
        assert!(core::ptr::eq(hi.clamp(&lo, &mid), &mid));
        for n in nans() {
            assert!(lo.max(&n).is_nan() && n.min(&lo).is_nan() && lo.clamp(&n, &hi).is_nan());
        }
        assert!(hi.clamp(&hi, &lo).is_nan());
    }

    #[test]
    fn neg_flips_the_sign_and_keeps_everything_else() {
        let mk = || ExactNum::from_u32(0xDEAD_BEEF, 40);
        let n = mk().neg();
        assert!(n.is_negative());
        assert_eq!(n.mantissa_max_bit_len(), mk().mantissa_max_bit_len());
        assert_eq!(n.abs_cmp(&mk()), Some(0));
        assert!(matches!(n.neg().compare(&mk()), Some(0)));
        let z = ExactNum::new(64).neg();
        assert!(z.is_zero() && z.is_negative());
        assert!(INF_POS.neg().is_inf() && INF_POS.neg().is_negative());
        assert!(INF_NEG.neg().is_inf() && INF_NEG.neg().is_positive());
        for n in nans() {
            same_nan(&n.neg_via_mul(), &n);
        }
        let third = || ExactNum::from_u8(1, P).div(&ExactNum::from_u8(3, P), 64, RM);
        assert!(third().inexact());
        let t = third().neg();
        assert!(t.inexact() && t.is_negative());
        assert!(matches!(
            (-ExactNum::from_u8(2, P)).compare(&ExactNum::from_u8(2, P).neg()),
            Some(0)
        ));
    }

    #[test]
    fn frexp_and_ldexp_keep_specials() {
        for n in nans() {
            let (m, e) = ExactNum::nan(n.err()).frexp();
            same_nan(&m, &n);
            assert_eq!(e, 0);
            let want = ExactNum::nan(n.err());
            same_nan(&n.ldexp(3, P, RM), &want);
        }
        let (m, e) = INF_POS.frexp();
        assert!(m.is_inf() && m.is_positive() && e == 0);
        let (m, e) = INF_NEG.frexp();
        assert!(m.is_inf() && m.is_negative() && e == 0);
        assert!(INF_POS.ldexp(3, P, RM).is_positive() && INF_POS.ldexp(3, P, RM).is_inf());
        assert!(INF_NEG.ldexp(3, P, RM).is_negative() && INF_NEG.ldexp(3, P, RM).is_inf());
    }
    #[test]
    fn neg_via_mul_equals_neg_for_every_kind_of_operand() {
        let same = |mk: &dyn Fn() -> ExactNum| {
            let x = mk();
            let want = mk().neg();
            let got = x.neg_via_mul();
            if want.is_nan() {
                same_nan(&got, &want);
                return;
            }
            assert_eq!(got.is_inf(), want.is_inf(), "{x:?}");
            assert_eq!(got.is_negative(), want.is_negative(), "{x:?}");
            assert_eq!(got.is_zero(), want.is_zero(), "{x:?}");
            // The accumulated flag is kept; the per-operation flag is false because the
            // product is exact (neg keeps the operand's stale per-operation flag).
            assert_eq!(
                got.inexact_accumulated(),
                want.inexact_accumulated(),
                "{x:?}"
            );
            assert!(!got.inexact(), "{x:?}");
            assert_eq!(got.compare(&want), want.compare(&got), "{x:?}");
            if !want.is_inf() {
                assert!(matches!(got.compare(&want), Some(0)), "{x:?}");
                assert_eq!(
                    got.mantissa_max_bit_len(),
                    want.mantissa_max_bit_len(),
                    "{x:?}"
                );
                assert_eq!(got.exponent(), want.exponent(), "{x:?}");
            }
        };
        for p in [1usize, 40, 64, 65, 128, 300] {
            same(&|| ExactNum::new(p));
            same(&|| ExactNum::new(p).neg());
            same(&|| ExactNum::from_u32(0xDEAD_BEEF, p));
            same(&|| ExactNum::from_u32(0xDEAD_BEEF, p).neg());
            same(&|| ExactNum::from_u8(1, p).div(&ExactNum::from_u8(3, p), p, RM));
            same(&|| {
                ExactNum::from_u8(1, p)
                    .div(&ExactNum::from_u8(3, p), p, RM)
                    .neg()
            });
            same(&|| ExactNum::from_u8(1, p).ldexp(-10_000, p, RM));
        }
        same(&|| INF_POS);
        same(&|| INF_NEG);
        for n in nans() {
            same_nan(&n.neg_via_mul(), &n);
        }
    }
    // ---- unary derivations that will take `self`: every call below is on a temporary, so the
    // same test compiles before and after the signature change.

    /// `k / 4`: an exact binary fraction.
    fn q(k: i64) -> ExactNum {
        ExactNum::from_i64(k, P).div(&ExactNum::from_u8(4, P), P, RM)
    }

    fn eq(a: &ExactNum, b: &ExactNum) -> bool {
        matches!(a.compare(b), Some(0))
    }

    #[test]
    fn abs_characterization() {
        let third = || ExactNum::from_u8(1, P).div(&ExactNum::from_u8(3, P), 64, RM);
        assert!(eq(&q(-15).abs(), &q(15)));
        assert!(eq(&q(15).abs(), &q(15)));
        assert!(!q(-15).abs().is_negative());
        let z = ExactNum::new(64).neg().abs();
        assert!(z.is_zero() && !z.is_negative());
        assert!(INF_NEG.abs().is_inf() && INF_NEG.abs().is_positive());
        assert!(INF_POS.abs().is_inf() && INF_POS.abs().is_positive());
        for n in nans() {
            let want = ExactNum::nan(n.err());
            same_nan(&n.abs(), &want);
        }
        let t = third().neg().abs();
        assert!(t.is_positive() && t.inexact_accumulated() && t.inexact() == third().inexact());
        assert_eq!(t.mantissa_max_bit_len(), third().mantissa_max_bit_len());
        assert_eq!(t.compare(&third()), Some(0));
    }

    #[test]
    fn int_and_fract_characterization() {
        assert!(eq(&q(15).int(), &q(12)));
        assert!(eq(&q(-15).int(), &q(-12)));
        assert!(eq(&q(15).fract(), &q(3)));
        assert!(eq(&q(-15).fract(), &q(-3)));
        assert!(q(3).int().is_zero());
        assert!(eq(&q(3).fract(), &q(3)));
        assert!(eq(&q(-1).fract(), &q(-1)));
        assert!(q(16).fract().is_zero());
        assert!(eq(&q(16).int(), &q(16)));
        let big = ExactNum::from_u8(1, P).ldexp(300, P, RM);
        assert!(big.fract().is_zero());
        assert!(eq(&big.int(), &ExactNum::from_u8(1, P).ldexp(300, P, RM)));
        assert!(ExactNum::new(P).int().is_zero() && ExactNum::new(P).fract().is_zero());
        assert!(INF_POS.int().is_nan() && INF_NEG.fract().is_nan());
        for n in nans() {
            same_nan(&n.int(), &n);
            same_nan(&n.fract(), &n);
        }
        let x = q(23);
        assert_eq!(x.fract().mantissa_max_bit_len(), x.mantissa_max_bit_len());
    }

    #[test]
    fn round_characterization() {
        assert!(eq(&q(7).round(1, RM), &q(8)));
        assert!(eq(&q(7).round(2, RM), &q(7)));
        assert!(eq(&q(5).round(1, RM), &q(4)));
        assert!(eq(&q(-7).round(1, RM), &q(-8)));
        assert!(eq(&q(7).round(0, RoundingMode::ToZero), &q(4)));
        assert!(INF_POS.round(3, RM).is_inf() && INF_POS.round(3, RM).is_positive());
        assert!(INF_NEG.round(3, RM).is_inf() && INF_NEG.round(3, RM).is_negative());
        for n in nans() {
            same_nan(&n.round(3, RM), &n);
        }
    }

    #[test]
    fn frexp_ldexp_characterization() {
        let (m, e) = q(48).frexp();
        assert!(eq(&m, &q(3)) && e == 4);
        let (m, e) = q(-48).frexp();
        assert!(eq(&m, &q(-3)) && e == 4);
        let (m, e) = q(1).frexp();
        assert!(eq(&m, &q(2)) && e == -1);
        let (m, e) = ExactNum::new(P).frexp();
        assert!(m.is_zero() && e == 0);
        assert!(eq(&q(12).ldexp(5, P, RM), &q(384)));
        assert!(eq(&q(12).ldexp(-2, P, RM), &q(3)));
        assert!(eq(&q(-12).ldexp(1, P, RM), &q(-24)));
        assert!(ExactNum::new(P).ldexp(7, P, RM).is_zero());
        assert!(q(4).ldexp(crate::Exponent::MAX, P, RM).is_inf());
        assert!(q(-4).ldexp(crate::Exponent::MAX, P, RM).is_negative());
        let third = || ExactNum::from_u8(1, P).div(&ExactNum::from_u8(3, P), 64, RM);
        let t = third().ldexp(3, 64, RM);
        assert!(t.inexact_accumulated());
        assert!(eq(&t, &third().mul(&q(32), 64, RM)));
        // rounding to a smaller precision rounds once
        let r = ExactNum::from_u32(0xFFFF_FFFF, P).ldexp(1, 8, RoundingMode::ToZero);
        assert!(r
            .compare(&ExactNum::from_u32(0xFFFF_FFFF, P).ldexp(1, P, RM))
            .is_some_and(|c| c < 0));
    }
}

/// Characterization of `add` / `sub` when an operand is zero (the "return the other operand" shortcut that rule R
/// deletes): value, sign of zero, rounding to `p`, flags, Inf / NaN. Every call is on a temporary or a borrowed
/// second operand, so the same test compiles before and after the signature change.
#[cfg(test)]
mod add_sub_zero_tests {
    use crate::defs::RoundingMode;
    use crate::{Error, ExactNum, INF_NEG, INF_POS};
    use core::ops::Neg;

    const RMS: [RoundingMode; 5] = [
        RoundingMode::ToEven,
        RoundingMode::ToZero,
        RoundingMode::Up,
        RoundingMode::Down,
        RoundingMode::FromZero,
    ];

    fn zero(neg: bool) -> ExactNum {
        let z = ExactNum::new(64);
        if neg {
            z.neg()
        } else {
            z
        }
    }

    fn third(p: usize) -> ExactNum {
        ExactNum::from_u8(1, p).div(&ExactNum::from_u8(3, p), p, RoundingMode::ToEven)
    }

    /// Full-precision sum with a zero operand is the other operand, exactly. (Before the zero-aware
    /// kernel branch this returned a value with an empty mantissa, flagged inexact.)
    #[test]
    fn full_prec_with_zero_operand_is_the_other_operand() {
        for x in [
            ExactNum::from_u8(1, 64),
            ExactNum::from_u8(7, 64).neg(),
            third(64),
            third(192).neg(),
        ] {
            for neg in [false, true] {
                let sum = zero(neg).add_full_prec(&x);
                let diff = zero(neg).sub_full_prec(&x);
                let sum2 = x.add_full_prec(&zero(neg));
                assert!(sum == x, "0 + x");
                assert_eq!(sum.inexact(), x.inexact(), "0 + x flag");
                assert!(diff == x.neg_via_mul(), "0 - x");
                assert!(sum2 == x, "x + 0");
                assert_eq!(sum2.exponent(), x.exponent());
            }
        }
    }

    #[test]
    fn zero_plus_zero_sign_rule() {
        for rm in RMS {
            for (an, bn) in [(false, false), (false, true), (true, false), (true, true)] {
                let want_add = if an == bn { an } else { rm == RoundingMode::Down };
                let r = zero(an).add(&zero(bn), 64, rm);
                assert!(r.is_zero());
                assert_eq!(r.is_negative(), want_add, "{an} + {bn} {rm:?}");
                let want_sub = if an != bn { an } else { rm == RoundingMode::Down };
                let r = zero(an).sub(&zero(bn), 64, rm);
                assert!(r.is_zero());
                assert_eq!(r.is_negative(), want_sub, "{an} - {bn} {rm:?}");
            }
        }
    }

    #[test]
    fn one_zero_operand_returns_the_other_rounded_to_p() {
        for rm in RMS {
            for p in [64, 128, 192] {
                for neg in [false, true] {
                    let x = if neg { third(192).neg() } else { third(192) };
                    // reference: the same value rounded to p by set_precision
                    let mut want = ExactNum::new(64).add(&x, 192, RoundingMode::None);
                    want.set_precision(p, rm).unwrap();
                    for zn in [false, true] {
                        let z = zero(zn);
                        let a = x.add(&z, p, rm);
                        let b = z.add(&x, p, rm);
                        let c = x.sub(&z, p, rm);
                        let d = z.sub(&x, p, rm);
                        for r in [&a, &b, &c] {
                            assert_eq!(r.compare(&want), Some(0), "p={p} {rm:?} neg={neg}");
                            assert_eq!(r.mantissa_max_bit_len(), want.mantissa_max_bit_len());
                        }
                        let mut want_neg = ExactNum::new(64).sub(&x, 192, RoundingMode::None);
                        want_neg.set_precision(p, rm).unwrap();
                        assert_eq!(d.compare(&want_neg), Some(0), "p={p} {rm:?}");
                        assert_eq!(a.is_negative(), neg);
                        assert_eq!(d.is_negative(), !neg);
                        // 192 bits of 1/3 do not fit 64 or 128 bits
                        if p < 192 {
                            assert!(a.inexact_accumulated() && b.inexact_accumulated());
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn zero_with_exact_operand_is_exact() {
        let x = ExactNum::from_u8(12, 64);
        for rm in RMS {
            let r = x.add(&zero(false), 64, rm);
            assert!(!r.inexact_accumulated());
            assert_eq!(r.compare(&x), Some(0));
            let r = zero(true).sub(&x, 64, rm);
            assert_eq!(r.compare(&x.neg_via_mul()), Some(0));
        }
    }

    #[test]
    fn zero_with_inf_and_nan() {
        let nan = ExactNum::nan(Some(Error::InvalidArgument));
        for rm in RMS {
            assert!(zero(false).add(&INF_POS, 64, rm).is_positive());
            assert!(zero(true).add(&INF_NEG, 64, rm).is_inf() && INF_NEG.is_negative());
            assert!(INF_NEG.sub(&zero(false), 64, rm).is_negative());
            assert!(zero(false).sub(&INF_POS, 64, rm).is_negative());
            assert!(zero(false).sub(&INF_POS, 64, rm).is_inf());
            let r = zero(false).add(&nan, 64, rm);
            assert!(r.is_nan() && r.err() == nan.err());
            let r = nan.sub(&zero(true), 64, rm);
            assert!(r.is_nan() && r.err() == nan.err());
        }
    }
}

/// Characterization over grids of operands: `ldexp` / `scalb` (one operand) and `add` / `sub` / `add_full_prec` /
/// `sub_full_prec` (pairs, including every zero, Inf and NaN combination). Original text: `ldexp` / `scalb` over a grid of operands (zero, normal, near both ends of the exponent
/// range), shifts, precisions and rounding modes: the digest of every result (value, sign, exponent, flags)
/// is pinned (storage representation, inline or heap, is not part of the value), so the signature change to `self`
/// must reproduce every bit. Both the `&self` and the `self` implementation give this digest.
#[cfg(test)]
mod grid_digest_tests {
    use crate::defs::{RoundingMode, EXPONENT_MAX, EXPONENT_MIN};
    use crate::{Error, ExactNum, INF_NEG, INF_POS};
    use core::hash::{Hash, Hasher};
    use core::ops::Neg;
    use std::collections::hash_map::DefaultHasher;

    const RMS: [RoundingMode; 6] = [
        RoundingMode::None,
        RoundingMode::ToEven,
        RoundingMode::ToZero,
        RoundingMode::Up,
        RoundingMode::Down,
        RoundingMode::FromZero,
    ];

    const OPERANDS: usize = 21;
    /// Operands 21 and 22 are subnormal.
    const PAIR_OPERANDS: usize = 23;

    /// Operand `i` of the grid, built fresh on every call (`ldexp` consumes it).
    fn operand(i: usize) -> ExactNum {
        let third =
            |p| ExactNum::from_u8(1, p).div(&ExactNum::from_u8(3, p), p, RoundingMode::ToEven);
        let ks = [
            EXPONENT_MIN + 3,
            EXPONENT_MIN + 70,
            EXPONENT_MIN + 200,
            EXPONENT_MAX - 3,
            EXPONENT_MAX - 70,
        ];
        match i {
            0 => ExactNum::new(64),
            1 => ExactNum::new(64).neg(),
            2 => ExactNum::from_u8(1, 64),
            3 => ExactNum::from_u8(7, 64).neg(),
            4 => third(64),
            5 => third(192).neg(),
            6 => third(320),
            7 => INF_POS,
            8 => INF_NEG,
            9 => ExactNum::nan(Some(Error::InvalidArgument)),
            10 => ExactNum::nan(None),
            21 => ExactNum::from_u8(3, 64)
                .ldexp(EXPONENT_MIN + 1, 64, RoundingMode::None)
                .ldexp(-20, 64, RoundingMode::None),
            22 => third(192)
                .ldexp(EXPONENT_MIN + 1, 192, RoundingMode::None)
                .ldexp(-100, 192, RoundingMode::None)
                .neg(),
            _ => {
                let k = ks[(i - 11) / 2];
                if (i - 11).is_multiple_of(2) {
                    ExactNum::from_u8(3, 64).ldexp(k, 64, RoundingMode::None)
                } else {
                    third(192).ldexp(k, 192, RoundingMode::None).neg()
                }
            }
        }
    }

    /// `{r:?}` without the storage representation: inline capacity and heap buffers print differently for the
    /// same words, and which one a result lives in is not part of its value.
    fn value_repr(r: &ExactNum) -> String {
        let d = format!("{r:?}");
        let mut out = String::new();
        let mut rest = d.as_str();
        while let Some(i) = rest.find("inner: ") {
            let (head, tail) = rest.split_at(i + "inner: ".len());
            out.push_str(head);
            if let Some(t) = tail.strip_prefix("Inline { data: [") {
                let close = t.find("], len: ").unwrap();
                let words: Vec<&str> = t[..close].split(", ").collect();
                let after = &t[close + "], len: ".len()..];
                let digits = after.find(|c: char| !c.is_ascii_digit()).unwrap();
                let len: usize = after[..digits].parse().unwrap();
                out.push_str(&format!("{:?}", &words[..len]));
                rest = &after[digits + " }".len()..];
            } else if let Some(t) = tail.strip_prefix("Heap([") {
                let close = t.find("])").unwrap();
                let words: Vec<&str> = t[..close].split(", ").collect();
                out.push_str(&format!("{words:?}"));
                rest = &t[close + "])".len()..];
            } else {
                rest = tail;
            }
        }
        out.push_str(rest);
        out
    }

    fn digest(r: &ExactNum, h: &mut DefaultHasher) {
        value_repr(r).hash(h);
        r.inexact().hash(h);
        r.inexact_accumulated().hash(h);
    }

    #[test]
    fn ldexp_and_scalb_grid_digest() {
        let shifts = [
            0,
            1,
            -1,
            5,
            -64,
            64,
            -200,
            200,
            1 << 20,
            -(1 << 20),
            EXPONENT_MAX / 2,
            EXPONENT_MIN / 2,
            EXPONENT_MAX,
            EXPONENT_MIN,
            EXPONENT_MAX - 1,
            EXPONENT_MIN + 1,
        ];
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        for i in 0..OPERANDS {
            for n in shifts {
                for p in [1usize, 53, 64, 65, 128, 192, 200] {
                    for rm in RMS {
                        digest(&operand(i).ldexp(n, p, rm), &mut h);
                        digest(&operand(i).scalb(n, p, rm), &mut h);
                        count += 2;
                    }
                }
            }
        }
        assert_eq!(count, OPERANDS * 16 * 7 * 6 * 2);
        println!("ldexp grid digest = {}", h.finish());
        assert_eq!(h.finish(), 3_896_020_040_264_525_112);
    }

    #[test]
    fn add_sub_grid_digest() {
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        for i in 0..PAIR_OPERANDS {
            for j in 0..PAIR_OPERANDS {
                for p in [1usize, 53, 64, 65, 128, 192, 200] {
                    for rm in RMS {
                        digest(&operand(i).add(&operand(j), p, rm), &mut h);
                        digest(&operand(i).sub(&operand(j), p, rm), &mut h);
                        count += 2;
                    }
                }
                // full precision of operands 2^32 apart would need a 2^32-bit mantissa
                let (ei, ej) = (operand(i).exponent(), operand(j).exponent());
                let near = match (ei, ej) {
                    (Some(a), Some(b)) => a.abs_diff(b) < 4096,
                    _ => true,
                };
                if near {
                    digest(&operand(i).add_full_prec(&operand(j)), &mut h);
                    digest(&operand(i).sub_full_prec(&operand(j)), &mut h);
                    count += 2;
                }
            }
        }
        assert_eq!(count, 45_030);
        println!("add/sub grid digest = {}", h.finish());
        // Differs from the digest pinned on the pre-kernel-branch code in exactly 40 results: full-precision
        // sums with a zero operand (see `full_prec_with_zero_operand_is_the_other_operand`).
        assert_eq!(h.finish(), 7_232_426_507_586_278_165);
    }
}
