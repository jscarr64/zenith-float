// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! ExactNum including finite numbers, NaN, and `Inf`.

use crate::common::util::log2_ceil;
use crate::common::util::{round_p, round_p_rm};
use crate::defs::Precision;
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
        if p < crate::defs::MIN_PRECISION {
            return Self::nan(Some(Error::InvalidArgument));
        }
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
    pub fn add(self, d2: &Self, p: usize, rm: RoundingMode) -> Self {
        self.add_op(d2, p, rm, false)
    }

    /// `a + b` rounded once to `p` bits under `rm`, with both operands borrowed: the sanctioned
    /// borrowing form of [`Self::add`], for any site where an operand is still needed afterwards
    /// (`ExactNum::add_ref(&x, &x, ..)` included). Both operands are read in place, never copied;
    /// the only buffers allocated are the result and the kernels' own working space (as for
    /// `add`). The result is bit-identical to `add`: value, sign of zero in every rounding mode,
    /// both inexact flags, Inf/NaN, subnormal results, and directional overflow (`±Inf` when `rm`
    /// rounds away from zero for the result's sign: `ToEven`, `FromZero`, `Up` for a positive and
    /// `Down` for a negative result; otherwise the largest finite value of `p` bits, signed).
    /// Errors are NaN values, as for every `ExactNum` operation. `tools/rules_scan.sh` (BORROW
    /// rule) fails on any copy in the call graph of `add_ref`/`sub_ref`/`mul_ref`/`div_ref`.
    ///
    /// Accuracy: correctly-rounded
    pub fn add_ref(a: &Self, b: &Self, p: usize, rm: RoundingMode) -> Self {
        a.add_op(b, p, rm, false)
    }

    /// `a - b` with both operands borrowed. Same contract as [`Self::add_ref`], for subtraction.
    ///
    /// Accuracy: correctly-rounded
    pub fn sub_ref(a: &Self, b: &Self, p: usize, rm: RoundingMode) -> Self {
        a.sub_op(b, p, rm, false)
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
    pub fn sub(self, d2: &Self, p: usize, rm: RoundingMode) -> Self {
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
    pub fn mul(self, d2: &Self, p: usize, rm: RoundingMode) -> Self {
        self.mul_op(d2, p, rm, false)
    }

    /// `a * b` rounded once with both operands borrowed. Same contract as [`Self::add_ref`], for
    /// multiplication; bit-identical to [`Self::mul`]. (For a square prefer [`Self::sqr`].)
    ///
    /// Accuracy: correctly-rounded
    pub fn mul_ref(a: &Self, b: &Self, p: usize, rm: RoundingMode) -> Self {
        a.mul_op(b, p, rm, false)
    }

    /// `self * self`, rounded to `p` bits: the square of a value that is still needed afterwards.
    ///
    /// A plain function, not a copy primitive: it reads `self` twice, builds a fresh result in its own buffer
    /// and does not consume `self`; `self` is not copied, a subnormal one included (the product kernel reads it as stored). Value, sign of zero in every rounding mode, both inexact
    /// flags, `op_inexact`, Inf, subnormal results and the error kind of a NaN operand are those of
    /// `ExactNum::mul_ref(&x, &x, p, rm)` (and of the old `x * x`), pinned by `sqr_grid_digest`. It replaces
    /// every `mul_ref(&x, &x, ..)` in the library.
    ///
    /// The function returns NaN if the precision `p` is incorrect.
    ///
    /// Accuracy: correctly-rounded
    pub fn sqr(&self, p: usize, rm: RoundingMode) -> Self {
        self.mul_op(self, p, rm, false)
    }

    /// `a / b` rounded once with both operands borrowed. Same contract as [`Self::add_ref`], for
    /// division; bit-identical to [`Self::div`] (`x / ±0` is a signed Inf, `0 / 0` a NaN).
    ///
    /// Accuracy: correctly-rounded
    pub fn div_ref(a: &Self, b: &Self, p: usize, rm: RoundingMode) -> Self {
        a.div_op(b, p, rm)
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
                let prod = ExactNum::mul_ref(self, b, p, RoundingMode::None);
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
    pub fn two_sum(self, b: &Self, p: usize, rm: RoundingMode) -> (Self, Self) {
        if self.is_nan() {
            return (Self::nan(self.err()), Self::nan(self.err()));
        }
        if b.is_nan() {
            return (Self::nan(b.err()), Self::nan(b.err()));
        }
        if self.is_inf() || b.is_inf() {
            return (self.add(b, p, rm), Self::new(round_p_rm(p, rm)));
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
            return (
                ExactNum::mul_ref(self, b, p, rm),
                Self::new(round_p_rm(p, rm)),
            );
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
            return Self::new(round_p_rm(p, rm));
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
            return Self::new(round_p_rm(p, rm));
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
            let prod = ExactNum::mul_ref(x, y, p_wrk, RoundingMode::None);
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
            return Self::new(round_p_rm(p, rm));
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
    pub fn div(self, d2: &Self, p: usize, rm: RoundingMode) -> Self {
        self.div_op(d2, p, rm)
    }

    /// Divides `self` by `d2` exactly: `Ok(Some(q))` when the quotient is representable with
    /// `p` bits (word-rounded) inside the exponent range, including the subnormal grid, and
    /// `Ok(None)` when it is not (it would round, overflow, or underflow). `q` is the value
    /// [`Self::div`] returns in every rounding mode, and `q.inexact()` is false.
    ///
    /// `±Inf / finite` is `Some(±Inf)` and `finite / ±Inf` is `Some(±0)`.
    ///
    /// ## Errors
    ///
    ///  - InvalidArgument: an incorrect precision (`p == 0` or too large). The precision is
    ///    checked first, so `x / 0` with `p == 0` is InvalidArgument, not DivisionByZero.
    ///  - Then, in this order: a NaN operand gives its own error kind (InvalidArgument when it
    ///    has none); `0 / 0` is InvalidArgument; `x / 0` with `x` non-zero (Inf included) is
    ///    DivisionByZero; `±Inf / ±Inf` (any signs) is InvalidArgument.
    ///
    /// Accuracy: exact
    pub fn div_exact(&self, d2: &Self, p: usize) -> Result<Option<Self>, Error> {
        ExactNumNumber::p_assertion(round_p(p))?;
        if let Flavor::NaN(err) = &self.inner {
            return Err(err.unwrap_or(Error::InvalidArgument));
        }
        if let Flavor::NaN(err) = &d2.inner {
            return Err(err.unwrap_or(Error::InvalidArgument));
        }
        if d2.is_zero() {
            return Err(if self.is_zero() {
                Error::InvalidArgument
            } else {
                Error::DivisionByZero
            });
        }
        // Toward zero: an overflow saturates to the largest finite value (marked inexact) and an
        // underflow rounds to zero or the subnormal grid (marked inexact), never to Inf.
        let q = self.div_op(d2, p, RoundingMode::ToZero);
        if let Flavor::NaN(err) = &q.inner {
            return Err(err.unwrap_or(Error::InvalidArgument));
        }
        let finite = matches!(self.inner, Flavor::Value(_)) && matches!(d2.inner, Flavor::Value(_));
        if finite && q.inexact() {
            Ok(None)
        } else {
            Ok(Some(q))
        }
    }

    fn div_op(&self, d2: &Self, p: usize, rm: RoundingMode) -> Self {
        match &self.inner {
            Flavor::Value(v1) => match &d2.inner {
                Flavor::Value(v2) => Self::result_to_ext_rm(
                    ExactNumNumber::exact_bits(p, rm, |q, r| v1.div_per_op(v2, q, r)),
                    v1.is_zero(),
                    v1.sign() == v2.sign(),
                    p,
                    rm,
                ),
                Flavor::Inf(s2) => {
                    // IEEE 754: the zero quotient takes the sign of the quotient.
                    let z = Self::new(v1.mantissa_max_bit_len());
                    if v1.sign() == *s2 {
                        z
                    } else {
                        -z
                    }
                }
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
    /// Accuracy: correctly-rounded (the remainder is exact)
    pub fn rem(&self, d2: &Self) -> Self {
        match &self.inner {
            Flavor::Value(v1) => match &d2.inner {
                Flavor::Value(v2) => {
                    // The remainder is exact; its flag reports this operation, which never
                    // rounds, and the accumulated flag keeps both operands' history.
                    let inputs = v1.inexact() || v2.inexact();
                    // `rem` never rounds: the flag it leaves on its result is the operands' history.
                    let r = v1.rem(v2).map(|mut r| {
                        r.set_inexact(false);
                        r.finish_flags(inputs, true);
                        r
                    });
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
                        // v1^(+-inf) as in IEEE 754 / MPFR: decided by |v1| against 1, whatever the sign of v1.
                        let val = v1.abs_cmp(&crate::common::consts::ONE);
                        if val == 0 {
                            Self::from_u8(1, round_p_rm(p, rm))
                        } else if (val > 0) == s2.is_positive() {
                            INF_POS
                        } else {
                            Self::new(round_p_rm(p, rm))
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
                            Self::from_u8(1, round_p_rm(p, rm))
                        } else if v2.is_positive() {
                            if s1.is_negative() && v2.is_odd_int() {
                                // v2 is odd and has no fractional part.
                                INF_NEG
                            } else {
                                INF_POS
                            }
                        } else {
                            // (+-inf)^(negative): a zero, negative only for -inf to an odd integer.
                            let mut z = Self::new(round_p_rm(p, rm));
                            if s1.is_negative() && v2.is_odd_int() {
                                z.set_sign(Sign::Neg);
                            }
                            z
                        }
                    }
                    Flavor::Inf(s2) => {
                        // inf^inf
                        if s2.is_positive() {
                            INF_POS
                        } else {
                            Self::new(round_p_rm(p, rm))
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
                    Self::from_u8(1, round_p_rm(p, rm))
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
            // (-0)^(negative odd n) is -Inf: the pole takes the sign of the base only for odd n.
            Flavor::Value(v1) => Self::result_to_ext_rm(
                v1.powsi(n, p, rm),
                false,
                !(v1.is_negative() && n & 1 == 1),
                p,
                rm,
            ),
            Flavor::Inf(s1) => {
                if n == 0 {
                    Self::from_u8(1, round_p_rm(p, rm))
                } else if n < 0 {
                    Self::new(round_p_rm(p, rm))
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
        if p < crate::defs::MIN_PRECISION {
            return Self::nan(Some(Error::InvalidArgument));
        }
        let mut r = self.log_word(n, p, rm, cc);
        // Two finite operands round once at `p` through `exact_bits` (per-operation flag); the other
        // results are fixed values that only carry the operands' history.
        let both_finite = matches!(
            (&self.inner, &n.inner),
            (Flavor::Value(_), Flavor::Value(_))
        );
        if let (false, Flavor::Value(v)) = (both_finite, &mut r.inner) {
            v.finish_flags(self.inexact_accumulated() || n.inexact_accumulated(), true);
        }
        r
    }

    /// Word-precision implementation of [`Self::log`].
    fn log_word(&self, n: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match &self.inner {
            Flavor::Value(v1) => {
                match &n.inner {
                    Flavor::Value(v2) => {
                        if v2.is_zero() {
                            return INF_NEG;
                        }
                        Self::result_to_ext(
                            ExactNumNumber::exact_bits(p, rm, |q, r| v1.log(v2, q, r, cc)).map(
                                |mut r| {
                                    // Pin the per-operation flag (a word-multiple `p` leaves it implicit).
                                    let rounded = r.op_inexact();
                                    r.set_op_inexact(rounded);
                                    r
                                },
                            ),
                            false,
                            true,
                        )
                    }
                    Flavor::Inf(s2) => {
                        // v1.log(inf)
                        if s2.is_positive() {
                            Self::new(round_p_rm(p, rm))
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

    /// A finite numeral whose radix exponent left the range: `MAX`/`Inf` by `rm` when above, the minimum
    /// subnormal or a signed zero by `rm` when below (any radix digit exponent below `EXPONENT_MIN` is
    /// below 2^(2 * EXPONENT_MIN), far under the smallest subnormal).
    fn parse_out_of_range(s: Sign, above: bool, p: usize, rm: RoundingMode) -> Self {
        if above {
            return Self::result_to_ext_rm(Err(Error::ExponentOverflow(s)), false, true, p, rm);
        }
        let tiny = ExactNumNumber::from_word(1, round_p(p) + WORD_BIT_SIZE).and_then(|mut t| {
            t.set_sign(s);
            let far = i128::from(crate::EXPONENT_MIN) * 2;
            t.place_scaled(far, true, p, rm)
        });
        Self::result_to_ext(tiny, false, true)
    }

    /// Word-precision implementation of [`Self::parse`].
    fn parse_word(s: &str, rdx: Radix, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        match crate::parser::parse(s, rdx) {
            Ok(ps) => {
                if ps.out_of_range() < 0 && rdx == Radix::Bin {
                    // 0.digits * 2^e with e below the range: exact digits, placed and rounded once.
                    let (m, s, e) = ps.raw_parts();
                    let placed =
                        ExactNumNumber::convert_from_radix(s, m, 0, rdx, usize::MAX, rm, cc)
                            .and_then(|mut v| {
                                let w = v.mantissa_max_bit_len().max(round_p(p) + WORD_BIT_SIZE);
                                v.set_precision(w, RoundingMode::None)?;
                                let at = i128::from(e) + i128::from(v.exponent());
                                v.place_scaled(at, false, p, rm)
                            });
                    Self::result_to_ext(placed, false, true)
                } else if ps.out_of_range() != 0 {
                    // A radix exponent past the range: overflow, or a magnitude far below the subnormals.
                    Self::parse_out_of_range(ps.sign(), ps.out_of_range() > 0, p, rm)
                } else if ps.is_inf() {
                    if ps.sign() == Sign::Pos {
                        INF_POS
                    } else {
                        INF_NEG
                    }
                } else if ps.is_nan() {
                    NAN
                } else {
                    let (m, s, e) = ps.raw_parts();
                    Self::result_to_ext_rm(
                        ExactNumNumber::convert_from_radix(s, m, e, rdx, p, rm, cc),
                        false,
                        true,
                        p,
                        rm,
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
            Flavor::Inf(s) => Self::from_i8(s.to_int(), round_p_rm(p, rm)),
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
        // Every public entry point that goes through here rejects `p = 0` (2.0: `MIN_PRECISION` is 1).
        if p < crate::defs::MIN_PRECISION {
            return Self::nan(Some(Error::InvalidArgument));
        }
        let direct = rm == RoundingMode::None || p >= usize::MAX / 4;
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
            // Inf does not depend on the precision; recompute for the mode's own overflow value.
            Flavor::Inf(_) => f(p, rm),
            // The wide pass failed. At a word-aligned p, f(p, rm) is already rounded at p
            // (the contract the direct branch above relies on). Otherwise it may carry more
            // than p bits (128 at p=113), so its value is used only when it is exact and fits
            // p bits exactly (setting the precision drops no bit); else the error propagates.
            Flavor::NaN(err) => match f(p, rm).inner {
                Flavor::Value(v) if p.is_multiple_of(WORD_BIT_SIZE) => Self {
                    inner: Flavor::Value(v),
                },
                Flavor::Value(mut v) if !v.op_inexact() => {
                    let fits = v.set_precision(p, rm).is_ok() && !v.inexact();
                    if fits {
                        Self {
                            inner: Flavor::Value(v),
                        }
                    } else {
                        Self::nan(err)
                    }
                }
                Flavor::Value(_) | Flavor::Inf(_) | Flavor::NaN(_) => Self::nan(err),
            },
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
    /// `ζ(self)` for real `self` at `p` bits with `rm`; see [`ExactNumNumber::zeta_real`].
    pub(crate) fn zeta_real_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> ExactNum {
        match &self.inner {
            Flavor::Value(v) => Self::result_to_ext_rm(v.zeta_real(p, rm, cc), false, true, p, rm),
            Flavor::Inf(s) => {
                if s.is_positive() {
                    match ExactNumNumber::from_word(1, p) {
                        Ok(v) => ExactNum {
                            inner: Flavor::Value(v),
                        },
                        Err(e) => Self::nan(Some(e)),
                    }
                } else {
                    NAN
                }
            }
            Flavor::NaN(err) => Self::nan(*err),
        }
    }

    /// `B(self, b)` at `p` bits with `rm`; see [`ExactNumNumber::beta`].
    pub(crate) fn beta_word(
        &self,
        b: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> ExactNum {
        match (&self.inner, &b.inner) {
            (Flavor::Value(x), Flavor::Value(y)) => {
                Self::result_to_ext_rm(x.beta(y, p, rm, cc), false, true, p, rm)
            }
            (Flavor::NaN(err), _) | (_, Flavor::NaN(err)) => Self::nan(*err),
            (Flavor::Inf(_), _) | (_, Flavor::Inf(_)) => NAN,
        }
    }

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

    /// Number of leading zero bits of the mantissa words (least significant word first), or `None` when they are all
    /// zero.
    fn mantissa_leading_zero_bits(words: &[Word]) -> Option<usize> {
        words
            .iter()
            .rev()
            .enumerate()
            .find(|(_, w)| **w != 0)
            .map(|(i, w)| i * WORD_BIT_SIZE + w.leading_zeros() as usize)
    }

    /// The exponent of the leading (most significant) set bit: `floor(log2(|self|))` for a finite nonzero number,
    /// that is `exponent() - 1` for a normal number. A subnormal number has leading zero bits in its mantissa, so its
    /// leading bit is below `EXPONENT_MIN - 1` (hence `i64`: `min_positive(p)` gives `EXPONENT_MIN - p`, rounded up
    /// to a word). `None` for zero, Inf and NaN.
    ///
    /// The sign is ignored. This is the convention of `exponent()` made explicit: a nonzero finite number is
    /// `0.m * 2^e = m * 2^(e - k)` with `m` of `k` bits and its top bit set, and its value lies in
    /// `[2^(e - 1), 2^e)`; the leading-bit exponent is `e - 1`.
    ///
    /// ```
    /// # use zenith_float_num::ExactNum;
    /// assert_eq!(ExactNum::from_word(3, 64).leading_bit_exponent(), Some(1));
    /// assert_eq!(ExactNum::from_word(1, 64).leading_bit_exponent(), Some(0));
    /// assert_eq!(ExactNum::new(64).leading_bit_exponent(), None);
    /// ```
    pub fn leading_bit_exponent(&self) -> Option<i64> {
        let (words, _, _, e, _) = self.value_parts()?;
        let lz = Self::mantissa_leading_zero_bits(words)?;
        e.checked_sub(1)?.checked_sub(i64::try_from(lz).ok()?)
    }

    /// The top 64 bits of the magnitude, the sticky bit below them, and the exponent, for a finite nonzero number:
    /// `(mant, sticky, exp)` with `mant` in `[2^63, 2^64)` (a subnormal number is shifted left until its top bit is
    /// set), `sticky` true when any bit below those 64 is set, and `exp` the exponent of the normalized value in the
    /// `0.m * 2^e` convention of `exponent()`, so `|self| = 0.mant * 2^exp` plus the bits below (`exp` equals
    /// `leading_bit_exponent() + 1`, and `exponent()` for a normal number). Independent of the word size. The sign is
    /// ignored; `None` for zero, Inf and NaN.
    ///
    /// ```
    /// # use zenith_float_num::ExactNum;
    /// assert_eq!(ExactNum::from_word(3, 64).top_u64_mantissa_and_sticky(), Some((3 << 62, false, 2)));
    /// ```
    pub fn top_u64_mantissa_and_sticky(&self) -> Option<(u64, bool, i64)> {
        let (words, _, _, e, _) = self.value_parts()?;
        let lz = Self::mantissa_leading_zero_bits(words)?;
        let skip = lz / WORD_BIT_SIZE;
        // enough words to hold 64 bits after the shift that removes the leading zeros of the top nonzero word
        let k = (64 + WORD_BIT_SIZE).div_ceil(WORD_BIT_SIZE);
        let window = words.iter().rev().skip(skip).take(k);
        let taken = window.len();
        let acc = window.fold(0u128, |a, w| (a << WORD_BIT_SIZE) | u128::from(*w));
        let acc = acc << ((k - taken) * WORD_BIT_SIZE);
        let norm = acc << acc.leading_zeros();
        let below = words.iter().rev().skip(skip + k).any(|w| *w != 0);
        let sticky = (norm as u64) != 0 || below;
        let exp = e.checked_sub(i64::try_from(lz).ok()?)?;
        Some(((norm >> 64) as u64, sticky, exp))
    }

    /// The exact number `0.mant * 2^exp = mant * 2^(exp - 64)` (positive), the `Mag` convention of zenith-bound:
    /// `mant` in `[2^63, 2^64)` gives a value in `[2^(exp - 1), 2^exp)`, as `exponent()` reads it. Any nonzero
    /// `mant` is accepted; `mant == 0` gives positive zero whatever `exp` is. The result has the precision `p`
    /// (rounded up to the word size, at least 64 bits) and is never rounded.
    ///
    /// ## Errors
    ///
    ///  - `ExponentOverflow(Sign::Pos)`: `exp` is above `EXPONENT_MAX` (the value would not be finite).
    ///  - `InvalidArgument`: the value is not representable exactly: it has bits below the smallest subnormal of
    ///    precision `p` (`2^(EXPONENT_MIN - p)`), or `p` is not a valid precision.
    ///
    /// ```
    /// # use zenith_float_num::ExactNum;
    /// let x = ExactNum::from_mag_parts(3 << 62, 2, 64).expect("in range");
    /// assert_eq!(x.leading_bit_exponent(), Some(1)); // 3
    /// ```
    pub fn from_mag_parts(mant: u64, exp: i64, p: usize) -> Result<Self, Error> {
        let pp = round_p(p.max(64));
        if mant == 0 {
            return Ok(Self::new(pp));
        }
        // `Exponent` is `i64`, so `exp` is checked against the 62-bit range directly; `grid` and
        // `total` stay inside `i64` because `pp <= MAX_PRECISION < 2^60`.
        if exp > crate::EXPONENT_MAX {
            return Err(Error::ExponentOverflow(Sign::Pos));
        }
        let grid = crate::EXPONENT_MIN - i64::try_from(pp).map_err(|_| Error::InvalidArgument)?;
        if exp < grid {
            // the highest bit of the value is below the smallest subnormal
            return Err(Error::InvalidArgument);
        }
        let total = exp - 64;
        // the first shift stays in the exponent range (the value is then normal), the second goes below it
        let n1: Exponent = total.max(crate::EXPONENT_MIN);
        let n2: Exponent = total - n1;
        let base = Self::from_u64_exact(mant, pp);
        let r = base
            .ldexp(n1, pp, RoundingMode::ToZero)
            .ldexp(n2, pp, RoundingMode::ToZero);
        if r.is_nan() || r.is_inf() || r.inexact() || r.is_zero() {
            return Err(Error::InvalidArgument);
        }
        Ok(r)
    }

    /// The exact number `2^exp` (positive) with precision `p` (see [`Self::from_mag_parts`]).
    ///
    /// ## Errors
    ///
    ///  - `ExponentOverflow(Sign::Pos)`: `exp` is above `EXPONENT_MAX - 1`.
    ///  - `InvalidArgument`: `exp` is below the smallest subnormal `2^(EXPONENT_MIN - p)`, or `p` is invalid.
    ///
    /// ```
    /// # use zenith_float_num::ExactNum;
    /// let x = ExactNum::pow2(-3, 64).expect("in range");
    /// assert_eq!(x.leading_bit_exponent(), Some(-3));
    /// ```
    pub fn pow2(exp: i64, p: usize) -> Result<Self, Error> {
        let e = exp
            .checked_add(1)
            .ok_or(Error::ExponentOverflow(Sign::Pos))?;
        Self::from_mag_parts(1 << 63, e, p)
    }

    /// `mant` as an exact number of precision `pp` (a multiple of the word size, at least 64).
    fn from_u64_exact(mant: u64, pp: usize) -> Self {
        if WORD_BIT_SIZE >= 64 {
            Self::from_word(mant as Word, pp)
        } else {
            let hi = Self::from_word((mant >> 32) as Word, pp).ldexp(32, pp, RoundingMode::ToZero);
            let lo = Self::from_word(mant as Word, pp);
            hi.add(&lo, pp, RoundingMode::ToZero)
        }
    }

    /// Returns the number of significant bits used in the mantissa, or None if `self` is Inf or NaN.
    /// Normal numbers use all bits of the mantissa.
    /// Subnormal numbers use fewer bits than the mantissa can hold.
    pub fn precision(&self) -> Option<Precision> {
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
        if p < crate::defs::MIN_PRECISION {
            return Self::nan(Some(Error::InvalidArgument));
        }
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

    /// Mantissa words, significant bits, sign, exponent and inexact flag of a finite value.
    pub(crate) fn value_parts(&self) -> Option<(&[Word], usize, Sign, Exponent, bool)> {
        match &self.inner {
            Flavor::Value(v) => Some(v.as_raw_parts()),
            Flavor::Inf(_) | Flavor::NaN(_) => None,
        }
    }

    /// Borrows the raw representation of `self` without copying; see [`RawParts`] for the
    /// contract. `ExactNum::from_raw_parts(&x.as_raw_parts())` rebuilds `x` bit for bit (a NaN
    /// comes back as a NaN without its error kind).
    pub fn as_raw_parts(&self) -> RawParts<'_> {
        match &self.inner {
            Flavor::Value(v) => {
                let (words, n, s, exponent, inexact) = v.as_raw_parts();
                RawParts {
                    kind: if n == 0 { RawKind::Zero } else { RawKind::Finite },
                    words,
                    precision: words.len() * WORD_BIT_SIZE,
                    negative: s.is_negative(),
                    exponent,
                    inexact,
                    op_inexact: v.op_inexact(),
                }
            }
            Flavor::Inf(s) => RawParts::special(RawKind::Inf, s.is_negative()),
            Flavor::NaN(_) => RawParts::special(RawKind::NaN, false),
        }
    }

    /// Builds a number from validated raw parts; see [`RawParts`] for the contract.
    ///
    /// ## Errors
    ///
    ///  - InvalidArgument: `precision` is 0 or differs from `words.len() * WORD_BIT_SIZE`;
    ///    `op_inexact` is set without `inexact`; `exponent` is outside `[EXPONENT_MIN, EXPONENT_MAX]` = `±(2^62 − 1)`; a `Finite`
    ///    mantissa is all zero, or is not normalized (top bit of the last word clear) while
    ///    `exponent != EXPONENT_MIN`; a `Zero` has a non-zero word, no words or a non-zero
    ///    exponent; an `Inf` / `NaN` has words, a precision, an exponent, `inexact` or `op_inexact` set, or a
    ///    `NaN` has `negative` set. A `Zero` may carry `inexact` and `op_inexact`.
    ///  - MemoryAllocation: the mantissa could not be allocated.
    pub fn from_raw_parts(parts: &RawParts<'_>) -> Result<Self, Error> {
        let sign = if parts.negative { Sign::Neg } else { Sign::Pos };
        if matches!(parts.kind, RawKind::Inf | RawKind::NaN) {
            if !parts.words.is_empty()
                || parts.precision != 0
                || parts.exponent != 0
                || parts.inexact
                || parts.op_inexact
            {
                return Err(Error::InvalidArgument);
            }
            return match parts.kind {
                RawKind::NaN if parts.negative => Err(Error::InvalidArgument),
                RawKind::NaN => Ok(Self::nan(None)),
                RawKind::Inf => Ok(Self {
                    inner: Flavor::Inf(sign),
                }),
                RawKind::Finite | RawKind::Zero => Err(Error::InvalidArgument),
            };
        }
        let width = parts
            .words
            .len()
            .checked_mul(WORD_BIT_SIZE)
            .ok_or(Error::InvalidArgument)?;
        if parts.precision == 0 || parts.precision != width {
            return Err(Error::InvalidArgument);
        }
        if parts.op_inexact && !parts.inexact {
            return Err(Error::InvalidArgument);
        }
        if parts.exponent < crate::EXPONENT_MIN || parts.exponent > crate::EXPONENT_MAX {
            return Err(Error::InvalidArgument);
        }
        let n = parts.significant_bits()?;
        let ok = match parts.kind {
            RawKind::Zero => n == 0 && parts.exponent == 0,
            RawKind::Finite => n > 0 && (n == width || parts.exponent == crate::EXPONENT_MIN),
            RawKind::Inf | RawKind::NaN => false,
        };
        if !ok {
            return Err(Error::InvalidArgument);
        }
        let mut v =
            ExactNumNumber::from_raw_parts(parts.words, n, sign, parts.exponent, parts.inexact)?;
        if parts.op_inexact != parts.inexact {
            v.set_op_inexact(parts.op_inexact);
        }
        Ok(Self::result_to_ext(Ok(v), false, true))
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
    pub fn mantissa_max_bit_len(&self) -> Option<Precision> {
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
                let mut ret = Self::new(round_p_rm(p, rm));
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

macro_rules! gen_wrapper_own {
    // function takes self by value and hands its storage to the result
    ($comment:literal, $fname:ident, $ret:ty, $pos_inf:block, $neg_inf:block $(,)?) => {
        #[doc=$comment]
        pub fn $fname(self) -> $ret {
            match self.inner {
                Flavor::Value(v) => {
                    let is_zero = v.is_zero();
                    Self::result_to_ext(v.$fname(), is_zero, true)
                }
                Flavor::Inf(s) => if s.is_positive() $pos_inf else $neg_inf,
                Flavor::NaN(err) => Self::nan(err),
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
                Flavor::Inf(s) => {
                    // Results from an infinite operand are built at the rm-aware precision (64-bit words for None).
                    let $p = round_p_rm($p, rm);
                    let _ = $p;
                    if s.is_positive() $pos_inf else $neg_inf
                }
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
                Flavor::Inf(s) => {
                    // Results from an infinite operand are built at the rm-aware precision (64-bit words for None).
                    let $p = round_p_rm($p, rm);
                    let _ = $p;
                    if s.is_positive() $pos_inf else $neg_inf
                }
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
                Flavor::Inf(s) => {
                    // Results from an infinite operand are built at the rm-aware precision (64-bit words for None).
                    let $p = round_p_rm($p, rm);
                    let _ = $p;
                    if s.is_positive() $pos_inf else $neg_inf
                }
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
                Flavor::Inf(s) => {
                    // Results from an infinite operand are built at the rm-aware precision (64-bit words for None).
                    let $p = round_p_rm($p, rm);
                    let _ = $p;
                    if s.is_positive() $pos_inf else $neg_inf
                }
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
    gen_wrapper_own!(
        "Returns the integer part of `self`.
Accuracy: correctly-rounded",
        int,
        Self,
        { NAN },
        { NAN },
    );
    gen_wrapper_own!(
        "Returns the fractional part of `self`.
Accuracy: correctly-rounded",
        fract,
        Self,
        { NAN },
        { NAN },
    );
    gen_wrapper_own!(
        "Returns the smallest integer greater than or equal to `self`.
Accuracy: correctly-rounded",
        ceil,
        Self,
        { INF_POS },
        { INF_NEG },
    );
    gen_wrapper_own!(
        "Returns the largest integer less than or equal to `self`.
Accuracy: correctly-rounded",
        floor,
        Self,
        { INF_POS },
        { INF_NEG },
    );
    /// Returns the rounded number with `n` binary positions in the fractional part of the number using rounding mode `rm`.
    /// Takes ownership; the result is written in the operand's storage.
    ///
    /// Accuracy: correctly-rounded
    pub fn round(self, n: usize, rm: RoundingMode) -> Self {
        match self.inner {
            Flavor::Value(v) => {
                let is_zero = v.is_zero();
                Self::result_to_ext(v.round(n, rm), is_zero, true)
            }
            Flavor::Inf(s) => {
                if s.is_positive() {
                    INF_POS
                } else {
                    INF_NEG
                }
            }
            Flavor::NaN(err) => Self::nan(err),
        }
    }
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
                    Self::nan(Some(Error::InvalidArgument))
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
    gen_wrapper_ovf_rm_cc!(
        "`ln |Γ(self)|` with precision `p`. Poles at zero and the negative integers yield +Inf.

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
                    Self::new(round_p_rm(p, rm))
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
                    Self::new(round_p_rm(p, rm))
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
                    Self::new(round_p_rm(p, rm))
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

// `add` and `sub` take their first operand by value, so there is no `&ExactNum + &ExactNum`: nothing could
// own the result's storage without copying an operand. `&a + b` is `b + a`; `&a - b` is `-b + a`.
impl core::ops::Add<&ExactNum> for ExactNum {
    type Output = ExactNum;

    fn add(self, rhs: &ExactNum) -> ExactNum {
        ExactNum::add(self, rhs, DEFAULT_P, RoundingMode::ToEven)
    }
}

impl core::ops::Add<ExactNum> for ExactNum {
    type Output = ExactNum;

    fn add(self, rhs: ExactNum) -> ExactNum {
        ExactNum::add(self, &rhs, DEFAULT_P, RoundingMode::ToEven)
    }
}

impl core::ops::Add<ExactNum> for &ExactNum {
    type Output = ExactNum;

    fn add(self, rhs: ExactNum) -> ExactNum {
        // `rhs + self` is `self + rhs` for every pair except two NaNs, where the first operand's error wins.
        if self.is_nan() && rhs.is_nan() {
            return ExactNum::nan(self.err());
        }
        ExactNum::add(rhs, self, DEFAULT_P, RoundingMode::ToEven)
    }
}

impl core::ops::Sub<&ExactNum> for ExactNum {
    type Output = ExactNum;

    fn sub(self, rhs: &ExactNum) -> ExactNum {
        ExactNum::sub(self, rhs, DEFAULT_P, RoundingMode::ToEven)
    }
}

impl core::ops::Sub<ExactNum> for ExactNum {
    type Output = ExactNum;

    fn sub(self, rhs: ExactNum) -> ExactNum {
        ExactNum::sub(self, &rhs, DEFAULT_P, RoundingMode::ToEven)
    }
}

impl core::ops::Sub<ExactNum> for &ExactNum {
    type Output = ExactNum;

    fn sub(self, rhs: ExactNum) -> ExactNum {
        // `-rhs + self` is `self - rhs` for every pair except two NaNs, where the first operand's error wins.
        if self.is_nan() && rhs.is_nan() {
            return ExactNum::nan(self.err());
        }
        ExactNum::add(-rhs, self, DEFAULT_P, RoundingMode::ToEven)
    }
}

// `mul` and `div` follow the same rule: `a * &b`, `a * b`, `a / &b`, `a / b` and `&a * b` exist; `&a * &b` and
// `&a / &b` do not. `&a * b` is `b * a` (the first operand's NaN error wins, as before). `&a / b` cannot be
// reordered, so it is the one place where the operator uses `div_ref`.
impl core::ops::Mul<&ExactNum> for ExactNum {
    type Output = ExactNum;

    fn mul(self, rhs: &ExactNum) -> ExactNum {
        ExactNum::mul(self, rhs, DEFAULT_P, RoundingMode::ToEven)
    }
}

impl core::ops::Mul<ExactNum> for ExactNum {
    type Output = ExactNum;

    fn mul(self, rhs: ExactNum) -> ExactNum {
        ExactNum::mul(self, &rhs, DEFAULT_P, RoundingMode::ToEven)
    }
}

impl core::ops::Mul<ExactNum> for &ExactNum {
    type Output = ExactNum;

    fn mul(self, rhs: ExactNum) -> ExactNum {
        // `rhs * self` is `self * rhs` for every pair except two NaNs, where the first operand's error wins.
        if self.is_nan() && rhs.is_nan() {
            return ExactNum::nan(self.err());
        }
        ExactNum::mul(rhs, self, DEFAULT_P, RoundingMode::ToEven)
    }
}

impl core::ops::Div<&ExactNum> for ExactNum {
    type Output = ExactNum;

    fn div(self, rhs: &ExactNum) -> ExactNum {
        ExactNum::div(self, rhs, DEFAULT_P, RoundingMode::ToEven)
    }
}

impl core::ops::Div<ExactNum> for ExactNum {
    type Output = ExactNum;

    fn div(self, rhs: ExactNum) -> ExactNum {
        ExactNum::div(self, &rhs, DEFAULT_P, RoundingMode::ToEven)
    }
}

impl core::ops::Div<ExactNum> for &ExactNum {
    type Output = ExactNum;

    fn div(self, rhs: ExactNum) -> ExactNum {
        ExactNum::div_ref(self, &rhs, DEFAULT_P, RoundingMode::ToEven)
    }
}

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

    /// round_exact at p=113: when the wide (192-bit) pass fails, the fallback must not
    /// return f(113)'s 128-bit value (1 + 2^-120 here) unrounded; it propagates the Err.
    #[test]
    fn test_round_exact_wide_err_never_unrounded() -> Result<(), crate::Error> {
        let w = 128;
        let one = crate::num::ExactNumNumber::from_word(1, w)?;
        let mut tiny = crate::num::ExactNumNumber::from_word(1, w)?;
        tiny.set_exponent(-119);
        let v = one.add(&tiny, w, crate::RoundingMode::None)?;
        let mut once = Some(v);
        let f = move |q: usize, _rm: crate::RoundingMode| {
            if q > w {
                return ExactNum::nan(Some(Error::PrecisionRetryExhausted));
            }
            match once.take() {
                Some(x) => ExactNum {
                    inner: super::Flavor::Value(x),
                },
                None => ExactNum::nan(None),
            }
        };
        let r = ExactNum::round_exact(false, 113, crate::RoundingMode::ToEven, f);
        assert_eq!(r.err(), Some(Error::PrecisionRetryExhausted));
        Ok(())
    }

    use crate::common::util::rand_p;
    use crate::defs::DEFAULT_P;
    use crate::Consts;
    use crate::Error;

    const ONE: ExactNum = ExactNum {
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
    fn test_decimal_formatting_round_trip() -> Result<(), crate::Error> {
        // Regression test: on 32-bit targets, decimal formatting must produce a string
        // that parses back to the original value without losing precision.
        let mut cc = Consts::new()?;
        let p = 53;
        let rm = RoundingMode::ToEven;
        let value = ExactNum::parse("1.0", Radix::Dec, p, rm, &mut cc);

        let formatted = value.format(Radix::Dec, rm, &mut cc)?;
        assert_eq!(formatted, "1.e+0");

        let reparsed = ExactNum::parse(&formatted, Radix::Dec, p, rm, &mut cc);
        assert_eq!(reparsed, value);
        Ok(())
    }

    #[test]
    fn finite_divided_by_infinity_is_a_zero_with_the_sign_of_the_quotient() {
        // IEEE 754 / MPFR: x / ±Inf = ±0 with sign(x) xor sign(divisor).
        let xs = [
            ExactNum::new(64),
            -ExactNum::new(64),
            ExactNum::from_u8(2, 64),
            -ExactNum::from_u8(7, 64),
        ];
        for x in &xs {
            for d in [INF_POS, INF_NEG] {
                let want_neg = x.is_negative() != d.is_negative();
                for q in [
                    ExactNum::div_ref(x, &d, 64, RoundingMode::ToEven),
                    ExactNum::div_ref(x, &d, 64, RoundingMode::Down),
                    x / if d.is_negative() { INF_NEG } else { INF_POS },
                ] {
                    assert!(q.is_zero(), "{x:?} / {d:?} = {q:?}");
                    assert_eq!(q.is_negative(), want_neg, "{x:?} / {d:?} = {q:?}");
                }
            }
        }
    }

    #[test]
    fn test_ext() -> Result<(), crate::Error> {
        let rm = RoundingMode::ToOdd;
        let mut cc = Consts::new()?;

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

        assert!(ONE.as_raw_parts().kind == crate::RawKind::Finite);
        assert!(INF_POS.as_raw_parts().kind == crate::RawKind::Inf);
        assert!(INF_NEG.as_raw_parts().kind == crate::RawKind::Inf);
        assert!(NAN.as_raw_parts().kind == crate::RawKind::NaN);

        assert!(ExactNum::add_ref(&ONE, &ONE, rand_p(), rm).compare(&TWO) == Some(0));
        assert!(ExactNum::add_ref(&ONE, &INF_POS, rand_p(), rm).is_inf_pos());
        assert!(ExactNum::add_ref(&INF_POS, &ONE, rand_p(), rm).is_inf_pos());
        assert!(ExactNum::add_ref(&ONE, &INF_NEG, rand_p(), rm).is_inf_neg());
        assert!(ExactNum::add_ref(&INF_NEG, &ONE, rand_p(), rm).is_inf_neg());
        assert!(ExactNum::add_ref(&INF_POS, &INF_POS, rand_p(), rm).is_inf_pos());
        assert!(ExactNum::add_ref(&INF_POS, &INF_NEG, rand_p(), rm).is_nan());
        assert!(ExactNum::add_ref(&INF_NEG, &INF_NEG, rand_p(), rm).is_inf_neg());
        assert!(ExactNum::add_ref(&INF_NEG, &INF_POS, rand_p(), rm).is_nan());

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

        assert!(ExactNum::sub_ref(&TWO, &ONE, rand_p(), rm).compare(&ONE) == Some(0));
        assert!(ExactNum::sub_ref(&ONE, &INF_POS, rand_p(), rm).is_inf_neg());
        assert!(ExactNum::sub_ref(&INF_POS, &ONE, rand_p(), rm).is_inf_pos());
        assert!(ExactNum::sub_ref(&ONE, &INF_NEG, rand_p(), rm).is_inf_pos());
        assert!(ExactNum::sub_ref(&INF_NEG, &ONE, rand_p(), rm).is_inf_neg());
        assert!(ExactNum::sub_ref(&INF_POS, &INF_POS, rand_p(), rm).is_nan());
        assert!(ExactNum::sub_ref(&INF_POS, &INF_NEG, rand_p(), rm).is_inf_pos());
        assert!(ExactNum::sub_ref(&INF_NEG, &INF_NEG, rand_p(), rm).is_nan());
        assert!(ExactNum::sub_ref(&INF_NEG, &INF_POS, rand_p(), rm).is_inf_neg());

        assert!(ExactNum::mul_ref(&TWO, &ONE, rand_p(), rm).compare(&TWO) == Some(0));
        assert!(ExactNum::mul_ref(&ONE, &INF_POS, rand_p(), rm).is_inf_pos());
        assert!(ExactNum::mul_ref(&INF_POS, &ONE, rand_p(), rm).is_inf_pos());
        assert!(ExactNum::mul_ref(&ONE, &INF_NEG, rand_p(), rm).is_inf_neg());
        assert!(ExactNum::mul_ref(&INF_NEG, &ONE, rand_p(), rm).is_inf_neg());
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
        assert!(ExactNum::mul_ref(&INF_POS, &INF_POS, rand_p(), rm).is_inf_pos());
        assert!(ExactNum::mul_ref(&INF_POS, &INF_NEG, rand_p(), rm).is_inf_neg());
        assert!(ExactNum::mul_ref(&INF_NEG, &INF_NEG, rand_p(), rm).is_inf_pos());
        assert!(INF_POS.sqr(rand_p(), rm).is_inf_pos());
        assert!(INF_NEG.sqr(rand_p(), rm).is_inf_pos());
        assert!(ExactNum::mul_ref(&INF_NEG, &INF_POS, rand_p(), rm).is_inf_neg());
        assert!(INF_POS.mul(&ExactNum::new(rand_p()), rand_p(), rm).is_nan());
        assert!(INF_NEG.mul(&ExactNum::new(rand_p()), rand_p(), rm).is_nan());
        assert!(ExactNum::new(rand_p()).mul(&INF_POS, rand_p(), rm).is_nan());
        assert!(ExactNum::new(rand_p()).mul(&INF_NEG, rand_p(), rm).is_nan());

        assert!(ExactNum::div_ref(&TWO, &TWO, rand_p(), rm).compare(&ONE) == Some(0));
        assert!(ExactNum::div_ref(&TWO, &INF_POS, rand_p(), rm).is_zero());
        assert!(ExactNum::div_ref(&INF_POS, &TWO, rand_p(), rm).is_inf_pos());
        assert!(ExactNum::div_ref(&TWO, &INF_NEG, rand_p(), rm).is_zero());
        assert!(ExactNum::div_ref(&INF_NEG, &TWO, rand_p(), rm).is_inf_neg());
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
        assert!(ExactNum::div_ref(&INF_POS, &INF_POS, rand_p(), rm).is_nan());
        assert!(ExactNum::div_ref(&INF_POS, &INF_NEG, rand_p(), rm).is_nan());
        assert!(ExactNum::div_ref(&INF_NEG, &INF_NEG, rand_p(), rm).is_nan());
        assert!(ExactNum::div_ref(&INF_NEG, &INF_POS, rand_p(), rm).is_nan());
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

        for op in [ExactNum::add_ref, ExactNum::sub_ref, ExactNum::mul_ref, ExactNum::div_ref] {
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

        assert!(ONE.compare(&ONE).ok_or(crate::Error::InvalidArgument)? == 0);
        assert!(ONE.compare(&INF_POS).ok_or(crate::Error::InvalidArgument)? < 0);
        assert!(INF_POS.compare(&ONE).ok_or(crate::Error::InvalidArgument)? > 0);
        assert!(
            INF_POS
                .compare(&INF_POS)
                .ok_or(crate::Error::InvalidArgument)?
                == 0
        );
        assert!(ONE.compare(&INF_NEG).ok_or(crate::Error::InvalidArgument)? > 0);
        assert!(INF_NEG.compare(&ONE).ok_or(crate::Error::InvalidArgument)? < 0);
        assert!(
            INF_NEG
                .compare(&INF_NEG)
                .ok_or(crate::Error::InvalidArgument)?
                == 0
        );
        assert!(
            INF_POS
                .compare(&INF_NEG)
                .ok_or(crate::Error::InvalidArgument)?
                > 0
        );
        assert!(
            INF_NEG
                .compare(&INF_POS)
                .ok_or(crate::Error::InvalidArgument)?
                < 0
        );
        assert!(
            INF_POS
                .compare(&INF_POS)
                .ok_or(crate::Error::InvalidArgument)?
                == 0
        );
        assert!(ONE.compare(&NAN).is_none());
        assert!(NAN.compare(&ONE).is_none());
        assert!(INF_POS.compare(&NAN).is_none());
        assert!(NAN.compare(&INF_POS).is_none());
        assert!(INF_NEG.compare(&NAN).is_none());
        assert!(NAN.compare(&INF_NEG).is_none());
        assert!(NAN.compare(&NAN).is_none());

        assert!(ONE.abs_cmp(&ONE).ok_or(crate::Error::InvalidArgument)? == 0);
        assert!(ONE.abs_cmp(&INF_POS).ok_or(crate::Error::InvalidArgument)? < 0);
        assert!(INF_POS.abs_cmp(&ONE).ok_or(crate::Error::InvalidArgument)? > 0);
        assert!(
            INF_POS
                .abs_cmp(&INF_POS)
                .ok_or(crate::Error::InvalidArgument)?
                == 0
        );
        assert!(ONE.abs_cmp(&INF_NEG).ok_or(crate::Error::InvalidArgument)? < 0);
        assert!(INF_NEG.abs_cmp(&ONE).ok_or(crate::Error::InvalidArgument)? > 0);
        assert!(
            INF_NEG
                .abs_cmp(&INF_NEG)
                .ok_or(crate::Error::InvalidArgument)?
                == 0
        );
        assert!(
            INF_POS
                .abs_cmp(&INF_NEG)
                .ok_or(crate::Error::InvalidArgument)?
                == 0
        );
        assert!(
            INF_NEG
                .abs_cmp(&INF_POS)
                .ok_or(crate::Error::InvalidArgument)?
                == 0
        );
        assert!(
            INF_POS
                .abs_cmp(&INF_POS)
                .ok_or(crate::Error::InvalidArgument)?
                == 0
        );
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
        assert!(neg_two.abs_cmp(&ONE).ok_or(crate::Error::InvalidArgument)? > 0);
        assert!(ONE.abs_cmp(&neg_two).ok_or(crate::Error::InvalidArgument)? < 0);
        assert!(neg_two.abs_cmp(&two).ok_or(crate::Error::InvalidArgument)? == 0);
        assert!(
            ExactNum::from_i8(-1, DEFAULT_P)
                .abs_cmp(&INF_POS)
                .ok_or(crate::Error::InvalidArgument)?
                < 0
        );
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
        // IEEE 754 / MPFR: pow(+0, -Inf) = +Inf (|x| < 1); checked against mpfr_pow in
        // tests/mpfr/pow_specials.rs (`pow(+0, -Inf)`).
        assert!(ExactNum::new(DEFAULT_P)
            .pow(&INF_NEG, rand_p(), rm, &mut cc)
            .is_inf_pos());
        assert!(ONE.pow(&INF_POS, rand_p(), rm, &mut cc).compare(&ONE) == Some(0));
        assert!(ONE.pow(&INF_NEG, rand_p(), rm, &mut cc).compare(&ONE) == Some(0));
        assert!(TWO.pow(&INF_POS, rand_p(), rm, &mut cc).is_inf_pos());
        // IEEE 754 / MPFR: pow(2, -Inf) = +0 (|x| > 1).
        let r = TWO.pow(&INF_NEG, rand_p(), rm, &mut cc);
        assert!(r.is_zero() && r.is_positive());
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

        let half = ExactNum::div_ref(&ONE, &TWO, rand_p(), rm);
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
        let mut half_pi: ExactNum = cc.pi_num(p, rm)?.into();
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

        let d1 = ONE;
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

        for mut special in [INF_POS, INF_NEG, NAN] {
            special.set_exponent(1);
            special.set_precision(1, rm)?;
            special.set_sign(Sign::Pos);
        }

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
        Ok(())
    }

    #[cfg(feature = "std")]
    #[test]
    fn test_ops_std() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;

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
        let mut d2 = ExactNum::from_str(&d1str)?;
        d2.set_precision(DEFAULT_P, RoundingMode::ToEven)?;
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
        let mut d2 = ExactNum::from_str(&d1str)?;
        d2.set_precision(DEFAULT_P, RoundingMode::ToEven)?;
        assert_eq!(d2, d1);

        let d1str = format!("{}", INF_POS);
        assert_eq!(d1str, "Inf");

        let d1str = format!("{}", INF_NEG);
        assert_eq!(d1str, "-Inf");

        let d1str = format!("{}", NAN);
        assert_eq!(d1str, "NaN");

        let bad = ExactNum::from_str("abc")?;
        assert!(bad.is_nan());
        assert_eq!(bad.err(), Some(Error::InvalidArgument));
        Ok(())
    }

    #[test]
    pub fn test_ops() -> Result<(), crate::Error> {
        use core::ops::Not;
        let mut cc = Consts::new()?;

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
        Ok(())
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
        let s = ExactNum::from_word(3, p).add(&b, p, RoundingMode::ToEven);
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

        let (hi, lo) = ExactNum::from_u8(1, p).two_sum(&two, p, rm);
        let rec = hi.add(&lo, p, rm);
        assert_eq!(rec.compare(&ExactNum::from_u8(3, p)), Some(0));

        let (ph, pl) = two.two_product(&three, p, rm);
        let pr = ph.add(&pl, p, rm);
        assert_eq!(pr.compare(&ExactNum::from_u8(6, p)), Some(0));

        let sum = ExactNum::fused_sum(
            &[ExactNum::from_u8(1, p), ExactNum::from_u8(2, p), ExactNum::from_u8(3, p)],
            p,
            rm,
        );
        assert_eq!(sum.compare(&ExactNum::from_u8(6, p)), Some(0));

        let dot = ExactNum::fused_dot(
            &[ExactNum::from_u8(1, p), ExactNum::from_u8(2, p)],
            &[ExactNum::from_u8(3, p), ExactNum::from_u8(1, p)],
            p,
            rm,
        );
        assert_eq!(dot.compare(&ExactNum::from_u8(5, p)), Some(0));

        // 1 + 2x + 3x² at x = 2 → 17
        let pv = ExactNum::polyval(&[one, two, three], &ExactNum::from_u8(2, p), p, rm);
        assert_eq!(pv.compare(&ExactNum::from_u8(17, p)), Some(0));
    }

    #[test]
    fn test_jacobi_sn_public() -> Result<(), crate::Error> {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let one = ExactNum::from_u8(1, p);
        let zero = ExactNum::from_u8(0, p);
        let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, rm);
        let sn0 = zero.jacobi_sn(&half, p, rm, &mut cc);
        assert!(sn0.is_zero(), "sn(0)");
        let cn0 = zero.jacobi_cn(&half, p, rm, &mut cc);
        assert_eq!(cn0.compare(&one), Some(0), "cn(0)");
        let sn_m0 = one.jacobi_sn(&zero, p, rm, &mut cc);
        let sin1 = one.sin(p, rm, &mut cc);
        let d = sn_m0.sub(&sin1, p, rm).abs();
        assert!(
            d.is_zero() || d.exponent().ok_or(crate::Error::InvalidArgument)? < -80,
            "sn(1|0)=sin 1"
        );
        assert!(one
            .jacobi_sn(&ExactNum::from_i8(2, p), p, rm, &mut cc)
            .is_nan());
        Ok(())
    }

    #[test]
    fn rm_none_results_from_infinite_operands_keep_64_bit_words_on_every_target(
    ) -> Result<(), crate::Error> {
        let mut cc = crate::Consts::new()?;
        let rm = crate::RoundingMode::None;
        let p = 200;
        let zero = ExactNum::from_u8(0, 64);
        let w = |x: &ExactNum| x.mantissa_max_bit_len();
        assert_eq!(w(&INF_NEG.exp(p, rm, &mut cc)), Some(256), "exp(-inf)");
        assert_eq!(w(&INF_NEG.exp2(p, rm, &mut cc)), Some(256), "exp2(-inf)");
        assert_eq!(w(&INF_POS.tanh(p, rm, &mut cc)), Some(256), "tanh(inf)");
        assert_eq!(w(&INF_POS.reciprocal(p, rm)), Some(256), "1/inf");
        assert_eq!(w(&INF_POS.powi(0, p, rm)), Some(256), "inf^0 (powi)");
        assert_eq!(w(&INF_POS.pow(&zero, p, rm, &mut cc)), Some(256), "inf^0");
        Ok(())
    }
}

#[cfg(feature = "random")]
#[cfg(test)]
mod rand_tests {

    use super::*;
    use crate::common::util::TEST_EXP_BOUND;

    #[test]
    fn test_rand() -> Result<(), crate::Error> {
        for _ in 0..100 {
            let p = crate::common::test_rng::random::<usize>() % 192 + DEFAULT_P;
            let exp_from = crate::common::test_rng::random::<Exponent>().abs() % TEST_EXP_BOUND;
            let span = (TEST_EXP_BOUND - exp_from).max(1);
            let exp_shift = crate::common::test_rng::random::<Exponent>().abs() % span;
            let exp_to = exp_from + exp_shift;

            let n = ExactNum::random_normal(p, exp_from, exp_to);

            assert!(!n.is_subnormal());
            assert!(
                n.exponent().ok_or(crate::Error::InvalidArgument)? >= exp_from
                    && n.exponent().ok_or(crate::Error::InvalidArgument)? <= exp_to
            );
            assert!(n.precision().ok_or(crate::Error::InvalidArgument)? >= p);
        }
        Ok(())
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

    /// `sqrt` of a negative number is NaN with the `InvalidArgument` kind whether the negative operand is
    /// finite or -Inf (it used to be a kind-less NaN for -Inf).
    #[test]
    fn sqrt_of_negative_inf_has_the_same_nan_kind_as_a_negative_finite() {
        let finite = ExactNum::from_u8(3, 64)
            .neg()
            .sqrt(64, RoundingMode::ToEven);
        assert!(finite.is_nan());
        assert_eq!(finite.err(), Some(Error::InvalidArgument));
        for rm in [RoundingMode::ToEven, RoundingMode::Up, RoundingMode::Down] {
            let inf = INF_NEG.sqrt(64, rm);
            assert!(inf.is_nan());
            assert_eq!(inf.err(), finite.err(), "sqrt(-Inf) NaN kind");
        }
        assert!(INF_POS.sqrt(64, RoundingMode::ToEven).is_inf_pos());
    }

    fn one() -> ExactNum {
        ExactNum::from_u8(1, P)
    }

    #[test]
    fn nan_operand_wins_in_ternary_and_binary_ops() -> Result<(), Box<dyn std::error::Error>> {
        let mut cc = Consts::new()?;
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
            let (hi, lo) = ExactNum::nan(n.err()).two_sum(&v, P, RM);
            same_nan(&hi, &n);
            same_nan(&lo, &n);
            let (hi, lo) = one().two_sum(&n, P, RM);
            same_nan(&hi, &n);
            same_nan(&lo, &n);
            let (hi, lo) = n.two_product(&v, P, RM);
            same_nan(&hi, &n);
            same_nan(&lo, &n);
            let (hi, lo) = v.two_product(&n, P, RM);
            same_nan(&hi, &n);
            same_nan(&lo, &n);
        }
        Ok(())
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
                // A None product keeps whole 64-bit words (rp(p, 64)) on every target.
                assert_eq!(
                    got.mantissa_max_bit_len(),
                    want.mantissa_max_bit_len()
                        .map(|b| crate::common::util::round_p_rm(b, RoundingMode::None)),
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
        let big = || ExactNum::from_u8(1, P).ldexp(300, P, RM);
        assert!(big().fract().is_zero());
        assert!(eq(&big().int(), &big()));
        assert!(ExactNum::new(P).int().is_zero() && ExactNum::new(P).fract().is_zero());
        assert!(INF_POS.int().is_nan() && INF_NEG.fract().is_nan());
        for (n, m) in nans().into_iter().zip(nans()) {
            same_nan(&n.int(), &m);
        }
        for (n, m) in nans().into_iter().zip(nans()) {
            same_nan(&n.fract(), &m);
        }
        assert_eq!(
            q(23).fract().mantissa_max_bit_len(),
            q(23).mantissa_max_bit_len()
        );
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
        for (n, m) in nans().into_iter().zip(nans()) {
            same_nan(&n.round(3, RM), &m);
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
    fn one_zero_operand_returns_the_other_rounded_to_p() -> Result<(), crate::Error> {
        for rm in RMS {
            for p in [64, 128, 192] {
                for neg in [false, true] {
                    let x = if neg { third(192).neg() } else { third(192) };
                    // reference: the same value rounded to p by set_precision
                    let mut want = ExactNum::new(64).add(&x, 192, RoundingMode::None);
                    want.set_precision(p, rm)?;
                    for zn in [false, true] {
                        let z = zero(zn);
                        let a = ExactNum::add_ref(&x, &z, p, rm);
                        let b = ExactNum::add_ref(&z, &x, p, rm);
                        let c = ExactNum::sub_ref(&x, &z, p, rm);
                        let d = z.sub(&x, p, rm);
                        for r in [&a, &b, &c] {
                            assert_eq!(r.compare(&want), Some(0), "p={p} {rm:?} neg={neg}");
                            assert_eq!(r.mantissa_max_bit_len(), want.mantissa_max_bit_len());
                        }
                        let mut want_neg = ExactNum::new(64).sub(&x, 192, RoundingMode::None);
                        want_neg.set_precision(p, rm)?;
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
        Ok(())
    }

    #[test]
    fn zero_with_exact_operand_is_exact() {
        let x = ExactNum::from_u8(12, 64);
        for rm in RMS {
            let r = ExactNum::from_u8(12, 64).add(&zero(false), 64, rm);
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
            let r = ExactNum::nan(nan.err()).sub(&zero(true), 64, rm);
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
    fn value_repr(r: &ExactNum) -> Result<String, crate::Error> {
        let d = format!("{r:?}");
        let mut out = String::new();
        let mut rest = d.as_str();
        while let Some(i) = rest.find("inner: ") {
            let (head, tail) = rest.split_at(i + "inner: ".len());
            out.push_str(head);
            if let Some(t) = tail.strip_prefix("Inline { data: [") {
                let close = t.find("], len: ").ok_or(crate::Error::InvalidArgument)?;
                let words: Vec<&str> = t[..close].split(", ").collect();
                let after = &t[close + "], len: ".len()..];
                let digits = after
                    .find(|c: char| !c.is_ascii_digit())
                    .ok_or(crate::Error::InvalidArgument)?;
                let len: usize = after[..digits]
                    .parse()
                    .map_err(|_| crate::Error::InvalidArgument)?;
                out.push_str(&format!("{:?}", &words[..len]));
                rest = &after[digits + " }".len()..];
            } else if let Some(t) = tail.strip_prefix("Heap([") {
                let close = t.find("])").ok_or(crate::Error::InvalidArgument)?;
                let words: Vec<&str> = t[..close].split(", ").collect();
                out.push_str(&format!("{words:?}"));
                rest = &t[close + "])".len()..];
            } else {
                rest = tail;
            }
        }
        out.push_str(rest);
        Ok(out)
    }

    fn digest(r: &ExactNum, h: &mut DefaultHasher) -> Result<(), Box<dyn std::error::Error>> {
        value_repr(r)?.hash(h);
        r.inexact().hash(h);
        r.inexact_accumulated().hash(h);
        Ok(())
    }

    /// `digest` that also appends the entry to `$ZF_DIGEST_DUMP/<test name>` when that variable is set: the
    /// operation, precision, mode, integer arguments, each operand (repr and inexact flag) and the result with
    /// both flags, one line per entry, for the mpmath verifier `scripts/verify_grid_digests.py`.
    fn dd(
        h: &mut DefaultHasher,
        op: &str,
        args: &[&ExactNum],
        ints: &[i128],
        p: usize,
        rm: RoundingMode,
        r: &ExactNum,
    ) -> Result<(), Box<dyn std::error::Error>> {
        if let Some(dir) = std::env::var_os("ZF_DIGEST_DUMP") {
            use std::io::Write;
            let t = std::thread::current();
            let name = t.name().map_or("main", |s| s).replace("::", ".");
            let mut line = format!("{op}\t{p}\t{rm:?}\t{ints:?}");
            for a in args {
                line.push_str(&format!("\t{}|{}", value_repr(a)?, a.inexact()));
            }
            line.push_str(&format!(
                "\t=\t{}\t{}\t{}\n",
                value_repr(r)?,
                r.inexact(),
                r.inexact_accumulated()
            ));
            std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(std::path::Path::new(&dir).join(name))?
                .write_all(line.as_bytes())?;
        }
        digest(r, h)
    }

    #[test]
    fn ldexp_and_scalb_grid_digest() -> Result<(), Box<dyn std::error::Error>> {
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
                        dd(
                            &mut h,
                            "ldexp",
                            &[&operand(i)],
                            &[i128::from(n)],
                            p,
                            rm,
                            &operand(i).ldexp(n, p, rm),
                        )?;
                        dd(
                            &mut h,
                            "ldexp",
                            &[&operand(i)],
                            &[i128::from(n)],
                            p,
                            rm,
                            &operand(i).scalb(n, p, rm),
                        )?;
                        count += 2;
                    }
                }
            }
        }
        assert_eq!(count, OPERANDS * 16 * 7 * 6 * 2);
        println!("ldexp grid digest = {}", h.finish());
        assert_eq!(h.finish(), 1_196_293_272_222_150_608);
        Ok(())
    }

    #[test]
    fn round_grid_digest() -> Result<(), Box<dyn std::error::Error>> {
        // 2^40 on a 64-bit target (the pinned digest); the same distance from the top of a narrower `usize`.
        let huge = 1usize << (usize::BITS - 24);
        let positions = [0usize, 1, 5, 63, 64, 65, 100, 200, 1 << 20, huge, usize::MAX / 2];
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        for i in 0..OPERANDS {
            for n in positions {
                for rm in RMS {
                    dd(
                        &mut h,
                        "round",
                        &[&operand(i)],
                        &[i128::try_from(n)?],
                        0,
                        rm,
                        &operand(i).round(n, rm),
                    )?;
                    count += 1;
                }
            }
        }
        assert_eq!(count, OPERANDS * 11 * 6);
        println!("round grid digest = {}", h.finish());
        // pinned on a 64-bit target; a narrower `usize` has a different position grid
        if cfg!(target_pointer_width = "64") {
            assert_eq!(h.finish(), 10_921_324_020_573_550_304);
        }
        Ok(())
    }

    /// The `+` and `-` operators on every combination of owned and borrowed operands that exists.
    #[test]
    fn operator_grid_digest() -> Result<(), Box<dyn std::error::Error>> {
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        for i in 0..PAIR_OPERANDS {
            for j in 0..PAIR_OPERANDS {
                dd(
                    &mut h,
                    "op+",
                    &[&operand(i), &operand(j)],
                    &[],
                    0,
                    RoundingMode::ToEven,
                    &(&operand(i) + operand(j)),
                )?;
                dd(
                    &mut h,
                    "op+",
                    &[&operand(i), &operand(j)],
                    &[],
                    0,
                    RoundingMode::ToEven,
                    &(operand(i) + &operand(j)),
                )?;
                dd(
                    &mut h,
                    "op+",
                    &[&operand(i), &operand(j)],
                    &[],
                    0,
                    RoundingMode::ToEven,
                    &(operand(i) + operand(j)),
                )?;
                dd(
                    &mut h,
                    "op-",
                    &[&operand(i), &operand(j)],
                    &[],
                    0,
                    RoundingMode::ToEven,
                    &(&operand(i) - operand(j)),
                )?;
                dd(
                    &mut h,
                    "op-",
                    &[&operand(i), &operand(j)],
                    &[],
                    0,
                    RoundingMode::ToEven,
                    &(operand(i) - &operand(j)),
                )?;
                dd(
                    &mut h,
                    "op-",
                    &[&operand(i), &operand(j)],
                    &[],
                    0,
                    RoundingMode::ToEven,
                    &(operand(i) - operand(j)),
                )?;
                count += 6;
            }
        }
        assert_eq!(count, PAIR_OPERANDS * PAIR_OPERANDS * 6);
        println!("operator grid digest = {}", h.finish());
        // 32-bit: generated by logs/dg/gen32.py (mpmath values, Word=u32 layout; 64-bit control reproduces the pin above)
        #[cfg(target_pointer_width = "64")]
        const WANT: u64 = 17_513_242_264_762_901_032;
        #[cfg(target_pointer_width = "32")]
        const WANT: u64 = 17_551_922_019_031_792_118;
        assert_eq!(h.finish(), WANT);
        Ok(())
    }

    /// `ExactComplex::add` and `sub` over parts that include zeros, Inf, NaN, inexact and subnormal values.
    #[test]
    fn complex_add_sub_grid_digest() -> Result<(), Box<dyn std::error::Error>> {
        use crate::ExactComplex;
        const COP: &str = "cplx_as";
        let parts = [0usize, 1, 3, 4, 5, 7, 9, 21];
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        for a in parts {
            for b in parts {
                for c in parts {
                    for d in parts {
                        for p in [1usize, 64, 128] {
                            for rm in RMS {
                                let z1 = || ExactComplex::new(operand(a), operand(b));
                                let z2 = || ExactComplex::new(operand(c), operand(d));
                                for (k, r) in [z1().add(&z2(), p, rm), z1().sub(&z2(), p, rm)]
                                    .into_iter()
                                    .enumerate()
                                {
                                    dd(
                                        &mut h,
                                        COP,
                                        &[&operand(a), &operand(b), &operand(c), &operand(d)],
                                        &[i128::try_from(k)?, 0],
                                        p,
                                        rm,
                                        r.re(),
                                    )?;
                                    dd(
                                        &mut h,
                                        COP,
                                        &[&operand(a), &operand(b), &operand(c), &operand(d)],
                                        &[i128::try_from(k)?, 1],
                                        p,
                                        rm,
                                        r.im(),
                                    )?;
                                    count += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(count, 8 * 8 * 8 * 8 * 3 * 6 * 2);
        println!("complex add/sub grid digest = {}", h.finish());
        #[cfg(target_pointer_width = "64")]
        let want = 18_040_698_040_155_693_760;
        #[cfg(target_pointer_width = "32")]
        let want = 10_425_817_710_396_386_987;
        assert_eq!(h.finish(), want);
        Ok(())
    }

    /// `ExactNum::add_ref` / `sub_ref` (both operands borrowed) hash to the same digest as the
    /// `add`/`sub` grid pinned on the pre-rule-R code (`add_sub_grid_digest`): same grid, same order.
    #[test]
    fn add_ref_sub_ref_grid_digest() -> Result<(), Box<dyn std::error::Error>> {
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        for i in 0..PAIR_OPERANDS {
            for j in 0..PAIR_OPERANDS {
                for p in [1usize, 53, 64, 65, 128, 192, 200] {
                    for rm in RMS {
                        dd(
                            &mut h,
                            "add",
                            &[&operand(i), &operand(j)],
                            &[],
                            p,
                            rm,
                            &ExactNum::add_ref(&operand(i), &operand(j), p, rm),
                        )?;
                        dd(
                            &mut h,
                            "sub",
                            &[&operand(i), &operand(j)],
                            &[],
                            p,
                            rm,
                            &ExactNum::sub_ref(&operand(i), &operand(j), p, rm),
                        )?;
                        count += 2;
                    }
                }
                let (ei, ej) = (operand(i).exponent(), operand(j).exponent());
                let near = match (ei, ej) {
                    (Some(a), Some(b)) => a.abs_diff(b) < 4096,
                    _ => true,
                };
                if near {
                    dd(
                        &mut h,
                        "add_full",
                        &[&operand(i), &operand(j)],
                        &[],
                        0,
                        RoundingMode::None,
                        &operand(i).add_full_prec(&operand(j)),
                    )?;
                    dd(
                        &mut h,
                        "sub_full",
                        &[&operand(i), &operand(j)],
                        &[],
                        0,
                        RoundingMode::None,
                        &operand(i).sub_full_prec(&operand(j)),
                    )?;
                    count += 2;
                }
            }
        }
        assert_eq!(count, 45_030);
        #[cfg(target_pointer_width = "64")]
        let want = 8_129_670_429_050_627_241;
        #[cfg(target_pointer_width = "32")]
        let want = 7_814_455_472_458_603_533;
        assert_eq!(h.finish(), want);
        Ok(())
    }

    /// `ExactComplex::add_ref` / `sub_ref` against the complex add/sub grid pinned on the old code.
    #[test]
    fn complex_add_ref_sub_ref_grid_digest() -> Result<(), Box<dyn std::error::Error>> {
        use crate::ExactComplex;
        const COP: &str = "cplx_as";
        let parts = [0usize, 1, 3, 4, 5, 7, 9, 21];
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        for a in parts {
            for b in parts {
                for c in parts {
                    for d in parts {
                        for p in [1usize, 64, 128] {
                            for rm in RMS {
                                let z1 = ExactComplex::new(operand(a), operand(b));
                                let z2 = ExactComplex::new(operand(c), operand(d));
                                for (k, r) in [
                                    ExactComplex::add_ref(&z1, &z2, p, rm),
                                    ExactComplex::sub_ref(&z1, &z2, p, rm),
                                ]
                                .into_iter()
                                .enumerate()
                                {
                                    dd(
                                        &mut h,
                                        COP,
                                        &[&operand(a), &operand(b), &operand(c), &operand(d)],
                                        &[i128::try_from(k)?, 0],
                                        p,
                                        rm,
                                        r.re(),
                                    )?;
                                    dd(
                                        &mut h,
                                        COP,
                                        &[&operand(a), &operand(b), &operand(c), &operand(d)],
                                        &[i128::try_from(k)?, 1],
                                        p,
                                        rm,
                                        r.im(),
                                    )?;
                                    count += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(count, 8 * 8 * 8 * 8 * 3 * 6 * 2);
        #[cfg(target_pointer_width = "64")]
        let want = 18_040_698_040_155_693_760;
        #[cfg(target_pointer_width = "32")]
        let want = 10_425_817_710_396_386_987;
        assert_eq!(h.finish(), want);
        Ok(())
    }

    #[test]
    fn add_sub_grid_digest() -> Result<(), Box<dyn std::error::Error>> {
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        for i in 0..PAIR_OPERANDS {
            for j in 0..PAIR_OPERANDS {
                for p in [1usize, 53, 64, 65, 128, 192, 200] {
                    for rm in RMS {
                        dd(
                            &mut h,
                            "add",
                            &[&operand(i), &operand(j)],
                            &[],
                            p,
                            rm,
                            &operand(i).add(&operand(j), p, rm),
                        )?;
                        dd(
                            &mut h,
                            "sub",
                            &[&operand(i), &operand(j)],
                            &[],
                            p,
                            rm,
                            &operand(i).sub(&operand(j), p, rm),
                        )?;
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
                    dd(
                        &mut h,
                        "add_full",
                        &[&operand(i), &operand(j)],
                        &[],
                        0,
                        RoundingMode::None,
                        &operand(i).add_full_prec(&operand(j)),
                    )?;
                    dd(
                        &mut h,
                        "sub_full",
                        &[&operand(i), &operand(j)],
                        &[],
                        0,
                        RoundingMode::None,
                        &operand(i).sub_full_prec(&operand(j)),
                    )?;
                    count += 2;
                }
            }
        }
        assert_eq!(count, 45_030);
        println!("add/sub grid digest = {}", h.finish());
        // Differs from the digest pinned on the pre-kernel-branch code in exactly 40 results: full-precision
        // sums with a zero operand (see `full_prec_with_zero_operand_is_the_other_operand`).
        #[cfg(target_pointer_width = "64")]
        let want = 8_129_670_429_050_627_241;
        #[cfg(target_pointer_width = "32")]
        let want = 7_814_455_472_458_603_533;
        assert_eq!(h.finish(), want);
        Ok(())
    }

    /// `mul` and `div` over every pair of grid operands, all precisions and rounding modes. Pinned on the
    /// code before `mul`/`div` took `self` and before `mul_ref`/`div_ref` existed.
    #[test]
    fn mul_div_grid_digest() -> Result<(), Box<dyn std::error::Error>> {
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        for i in 0..PAIR_OPERANDS {
            for j in 0..PAIR_OPERANDS {
                for p in [1usize, 53, 64, 65, 128, 192, 200] {
                    for rm in RMS {
                        dd(
                            &mut h,
                            "mul",
                            &[&operand(i), &operand(j)],
                            &[],
                            p,
                            rm,
                            &operand(i).mul(&operand(j), p, rm),
                        )?;
                        dd(
                            &mut h,
                            "div",
                            &[&operand(i), &operand(j)],
                            &[],
                            p,
                            rm,
                            &operand(i).div(&operand(j), p, rm),
                        )?;
                        count += 2;
                    }
                }
                dd(
                    &mut h,
                    "mul_full",
                    &[&operand(i), &operand(j)],
                    &[],
                    0,
                    RoundingMode::None,
                    &operand(i).mul_full_prec(&operand(j)),
                )?;
                count += 1;
            }
        }
        assert_eq!(count, PAIR_OPERANDS * PAIR_OPERANDS * (7 * 6 * 2 + 1));
        println!("mul/div grid digest = {}", h.finish());
        assert_eq!(h.finish(), 8_827_555_224_262_690_728);
        Ok(())
    }

    /// `ExactComplex::mul` and `div` over parts that include zeros, Inf, NaN, inexact and subnormal values.
    #[test]
    fn complex_mul_div_grid_digest() -> Result<(), Box<dyn std::error::Error>> {
        use crate::ExactComplex;
        const COP: &str = "cplx_md";
        let parts = [0usize, 1, 3, 4, 5, 7, 9, 21];
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        for a in parts {
            for b in parts {
                for c in parts {
                    for d in parts {
                        for p in [1usize, 64, 128] {
                            for rm in RMS {
                                let z1 = || ExactComplex::new(operand(a), operand(b));
                                let z2 = || ExactComplex::new(operand(c), operand(d));
                                for (k, r) in [z1().mul(&z2(), p, rm), z1().div(&z2(), p, rm)]
                                    .into_iter()
                                    .enumerate()
                                {
                                    dd(
                                        &mut h,
                                        COP,
                                        &[&operand(a), &operand(b), &operand(c), &operand(d)],
                                        &[i128::try_from(k)?, 0],
                                        p,
                                        rm,
                                        r.re(),
                                    )?;
                                    dd(
                                        &mut h,
                                        COP,
                                        &[&operand(a), &operand(b), &operand(c), &operand(d)],
                                        &[i128::try_from(k)?, 1],
                                        p,
                                        rm,
                                        r.im(),
                                    )?;
                                    count += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(count, 8 * 8 * 8 * 8 * 3 * 6 * 2);
        println!("complex mul/div grid digest = {}", h.finish());
        #[cfg(target_pointer_width = "64")]
        let want = 11_266_967_239_943_518_429;
        #[cfg(target_pointer_width = "32")]
        let want = 6_964_102_913_802_569_992;
        assert_eq!(h.finish(), want);
        Ok(())
    }

    /// `ExactNum::mul_ref` / `div_ref` (both operands borrowed) against the digest pinned on the old code.
    #[test]
    fn mul_ref_div_ref_grid_digest() -> Result<(), Box<dyn std::error::Error>> {
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        for i in 0..PAIR_OPERANDS {
            for j in 0..PAIR_OPERANDS {
                for p in [1usize, 53, 64, 65, 128, 192, 200] {
                    for rm in RMS {
                        dd(
                            &mut h,
                            "mul",
                            &[&operand(i), &operand(j)],
                            &[],
                            p,
                            rm,
                            &ExactNum::mul_ref(&operand(i), &operand(j), p, rm),
                        )?;
                        dd(
                            &mut h,
                            "div",
                            &[&operand(i), &operand(j)],
                            &[],
                            p,
                            rm,
                            &ExactNum::div_ref(&operand(i), &operand(j), p, rm),
                        )?;
                        count += 2;
                    }
                }
                dd(
                    &mut h,
                    "mul_full",
                    &[&operand(i), &operand(j)],
                    &[],
                    0,
                    RoundingMode::None,
                    &operand(i).mul_full_prec(&operand(j)),
                )?;
                count += 1;
            }
        }
        assert_eq!(count, PAIR_OPERANDS * PAIR_OPERANDS * (7 * 6 * 2 + 1));
        assert_eq!(h.finish(), 8_827_555_224_262_690_728);
        Ok(())
    }

    /// Squares through `mul_ref(&x, &x)`: the 23 grid operands (NaN of both kinds, +-Inf, +-0, the subnormal
    /// band, the exponent extremes) at 7 precisions in all 6 rounding modes. Pinned on the code before
    /// `sqr` existed; `sqr_grid_digest` below must hash to the same value.
    #[test]
    fn square_via_mul_ref_grid_digest() -> Result<(), Box<dyn std::error::Error>> {
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        for i in 0..PAIR_OPERANDS {
            for p in [1usize, 24, 53, 64, 65, 128, 192, 200, 320] {
                for rm in RMS {
                    let x = operand(i);
                    dd(
                        &mut h,
                        "sqr",
                        &[&x],
                        &[],
                        p,
                        rm,
                        &ExactNum::mul_ref(&x, &x, p, rm),
                    )?;
                    count += 1;
                }
            }
        }
        assert_eq!(count, PAIR_OPERANDS * 9 * 6);
        println!("square grid digest = {}", h.finish());
        assert_eq!(h.finish(), 17_269_539_713_433_390_012);
        Ok(())
    }

    /// `ExactNum::sqr(&x, p, rm)` against the digest pinned on `mul_ref(&x, &x, p, rm)`: value, sign of zero,
    /// both inexact flags, `op_inexact`, NaN error kind, Inf, subnormal results.
    #[test]
    fn sqr_grid_digest() -> Result<(), Box<dyn std::error::Error>> {
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        for i in 0..PAIR_OPERANDS {
            for p in [1usize, 24, 53, 64, 65, 128, 192, 200, 320] {
                for rm in RMS {
                    dd(
                        &mut h,
                        "sqr",
                        &[&operand(i)],
                        &[],
                        p,
                        rm,
                        &operand(i).sqr(p, rm),
                    )?;
                    count += 1;
                }
            }
        }
        assert_eq!(count, PAIR_OPERANDS * 9 * 6);
        assert_eq!(h.finish(), 17_269_539_713_433_390_012);
        Ok(())
    }

    /// Layer 0a pin: the trigonometric and hyperbolic kernels (argument reduction / restoration included),
    /// hashed over the operand grid (NaN kinds, infinities, zeros, subnormals and a few large arguments),
    /// recorded on the code that still copied the working value. A refactor that consumes `self` must keep it.
    #[test]
    fn trig_hyperbolic_grid_digest() -> Result<(), Box<dyn std::error::Error>> {
        let mut cc = crate::Consts::new()?;
        let extra = |k: usize| match k {
            0 => ExactNum::from_u32(100_000, 64),
            1 => ExactNum::from_u8(5, 64).ldexp(-9, 64, RoundingMode::None),
            2 => ExactNum::from_u8(3, 192)
                .ldexp(6, 192, RoundingMode::None)
                .neg(),
            _ => ExactNum::from_u8(1, 64)
                .div(&ExactNum::from_u8(7, 64), 64, RoundingMode::ToEven)
                .ldexp(-4, 64, RoundingMode::None),
        };
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        for i in (0..11).chain(21..23).chain(100..104) {
            for p in [53usize, 128, 200] {
                for rm in RMS {
                    let x = || if i < 100 { operand(i) } else { extra(i - 100) };
                    dd(&mut h, "sin", &[&x()], &[], p, rm, &x().sin(p, rm, &mut cc))?;
                    dd(&mut h, "cos", &[&x()], &[], p, rm, &x().cos(p, rm, &mut cc))?;
                    dd(&mut h, "tan", &[&x()], &[], p, rm, &x().tan(p, rm, &mut cc))?;
                    dd(
                        &mut h,
                        "atan",
                        &[&x()],
                        &[],
                        p,
                        rm,
                        &x().atan(p, rm, &mut cc),
                    )?;
                    dd(
                        &mut h,
                        "sinh",
                        &[&x()],
                        &[],
                        p,
                        rm,
                        &x().sinh(p, rm, &mut cc),
                    )?;
                    dd(
                        &mut h,
                        "cosh",
                        &[&x()],
                        &[],
                        p,
                        rm,
                        &x().cosh(p, rm, &mut cc),
                    )?;
                    dd(
                        &mut h,
                        "tanh",
                        &[&x()],
                        &[],
                        p,
                        rm,
                        &x().tanh(p, rm, &mut cc),
                    )?;
                    count += 7;
                }
            }
        }
        assert_eq!(count, 17 * 3 * 6 * 7);
        println!("trig hyperbolic grid digest = {}", h.finish());
        assert_eq!(h.finish(), 2_882_922_101_969_868_672);
        Ok(())
    }

    /// Layer 0a pin: `exp`, `exp2`, `exp10`, `pow` (integer and fractional exponent) and `powi`, hashed over the
    /// operand grid on the code that still copied the constants `ONE`/`TWO`/`TEN`.
    #[test]
    fn exp_pow_grid_digest() -> Result<(), Box<dyn std::error::Error>> {
        let mut cc = crate::Consts::new()?;
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        let exps = || {
            [
                ExactNum::from_u8(3, 64),
                ExactNum::from_u8(5, 64).ldexp(-1, 64, RoundingMode::None),
                ExactNum::from_u8(7, 64)
                    .ldexp(-5, 64, RoundingMode::None)
                    .neg(),
                ExactNum::new(64),
            ]
        };
        for i in (0..11).chain(21..23) {
            for p in [53usize, 128, 200] {
                for rm in RMS {
                    dd(
                        &mut h,
                        "exp",
                        &[&operand(i)],
                        &[],
                        p,
                        rm,
                        &operand(i).exp(p, rm, &mut cc),
                    )?;
                    dd(
                        &mut h,
                        "exp2",
                        &[&operand(i)],
                        &[],
                        p,
                        rm,
                        &operand(i).exp2(p, rm, &mut cc),
                    )?;
                    dd(
                        &mut h,
                        "exp10",
                        &[&operand(i)],
                        &[],
                        p,
                        rm,
                        &operand(i).exp10(p, rm, &mut cc),
                    )?;
                    for e in exps() {
                        dd(
                            &mut h,
                            "pow",
                            &[&operand(i), &e],
                            &[],
                            p,
                            rm,
                            &operand(i).pow(&e, p, rm, &mut cc),
                        )?;
                        count += 1;
                    }
                    for n in [0usize, 1, 2, 3, 7, 12] {
                        dd(
                            &mut h,
                            "powi",
                            &[&operand(i)],
                            &[i128::try_from(n)?],
                            p,
                            rm,
                            &operand(i).powi(n, p, rm),
                        )?;
                        count += 1;
                    }
                    count += 3;
                }
            }
        }
        assert_eq!(count, 13 * 3 * 6 * 13);
        println!("exp pow grid digest = {}", h.finish());
        assert_eq!(h.finish(), 2_398_075_299_501_098_093);
        Ok(())
    }

    /// Layer 0a pin: high-precision series (the rectangular scheme caches powers of the step) for `sin`, `exp`,
    /// `atan` and `sinh` at 1024 and 4096 bits.
    #[test]
    fn series_high_precision_digest() -> Result<(), Box<dyn std::error::Error>> {
        let mut cc = crate::Consts::new()?;
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        for i in [2usize, 3, 4, 5, 21] {
            for p in [1024usize, 4096] {
                for rm in [RoundingMode::ToEven, RoundingMode::Up] {
                    dd(
                        &mut h,
                        "sin",
                        &[&operand(i)],
                        &[],
                        p,
                        rm,
                        &operand(i).sin(p, rm, &mut cc),
                    )?;
                    dd(
                        &mut h,
                        "exp",
                        &[&operand(i)],
                        &[],
                        p,
                        rm,
                        &operand(i).exp(p, rm, &mut cc),
                    )?;
                    dd(
                        &mut h,
                        "atan",
                        &[&operand(i)],
                        &[],
                        p,
                        rm,
                        &operand(i).atan(p, rm, &mut cc),
                    )?;
                    dd(
                        &mut h,
                        "sinh",
                        &[&operand(i)],
                        &[],
                        p,
                        rm,
                        &operand(i).sinh(p, rm, &mut cc),
                    )?;
                    count += 4;
                }
            }
        }
        assert_eq!(count, 5 * 2 * 2 * 4);
        println!("series high precision digest = {}", h.finish());
        // 32-bit: generated by logs/dg/gen32.py (mpmath values, Word=u32 layout; 64-bit control reproduces the pin above)
        #[cfg(target_pointer_width = "64")]
        const WANT: u64 = 2_612_837_872_790_074_070;
        #[cfg(target_pointer_width = "32")]
        const WANT: u64 = 2_481_807_751_546_904_152;
        assert_eq!(h.finish(), WANT);
        Ok(())
    }

    /// Layer 0a pin for `fma`: the addend is shifted far above and far below the product so that both
    /// stand-in ("tiny") branches and the plain branch are all taken, over every addend operand
    /// (zeros, infinities and NaNs included).
    #[test]
    fn fma_stand_in_branches_digest() -> Result<(), Box<dyn std::error::Error>> {
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        for a in [2usize, 21] {
            for b in [3usize, 9] {
                for c in 0..OPERANDS {
                    for shift in [-3000i32, 0, 3000] {
                        for p in [64usize, 200] {
                            for rm in [RoundingMode::ToEven, RoundingMode::Down] {
                                let cs = operand(c).ldexp(
                                    shift as crate::Exponent,
                                    4096,
                                    RoundingMode::ToEven,
                                );
                                dd(
                                    &mut h,
                                    "fma",
                                    &[&operand(a), &operand(b), &cs],
                                    &[],
                                    p,
                                    rm,
                                    &operand(a).fma(&operand(b), &cs, p, rm),
                                )?;
                                count += 1;
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(count, 2 * 2 * OPERANDS * 3 * 2 * 2);
        println!("fma stand-in digest = {}", h.finish());
        // 32-bit: generated by logs/dg/gen32.py (mpmath values, Word=u32 layout; 64-bit control reproduces the pin above)
        #[cfg(target_pointer_width = "64")]
        const WANT: u64 = 6_276_170_101_586_996_179;
        #[cfg(target_pointer_width = "32")]
        const WANT: u64 = 1_593_121_079_249_968_054;
        assert_eq!(h.finish(), WANT);
        Ok(())
    }

    /// Layer 0a pin for `hypot`: the second argument is shifted far below, a little below, level with and far
    /// above the first, so the stand-in (`dy`) path and both orderings are taken, over every first operand
    /// (zeros, infinities and NaNs included).
    /// A product below the subnormal grid and an addend near `EXPONENT_MAX` are almost `2^63` places apart: the
    /// gap must not wrap in `i64`, and the product is only a sticky bit below the addend.
    #[test]
    fn fma_of_a_subnormal_product_and_a_huge_addend_is_the_addend(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let c = operand(17);
        for p in [64usize, 200] {
            let got = operand(21).fma(&operand(3), &c, p, RoundingMode::ToEven);
            assert_eq!(got.compare(&c), Some(0), "p={p}: {got:?}");
            assert!(got.inexact(), "p={p}");
            let down = operand(21).fma(&operand(3), &c, p, RoundingMode::Down);
            assert_eq!(down.compare(&c), Some(-1), "p={p}: {down:?}");
        }
        Ok(())
    }

    #[test]
    fn hypot_scaling_digest() -> Result<(), Box<dyn std::error::Error>> {
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        for a in 0..OPERANDS {
            for b in [2usize, 3, 9, 21] {
                for shift in [-3000i32, -200, 0, 3000] {
                    for p in [64usize, 200] {
                        for rm in [RoundingMode::ToEven, RoundingMode::Up] {
                            let bs = operand(b).ldexp(
                                shift as crate::Exponent,
                                256,
                                RoundingMode::ToEven,
                            );
                            dd(
                                &mut h,
                                "hypot",
                                &[&operand(a), &bs],
                                &[],
                                p,
                                rm,
                                &operand(a).hypot(&bs, p, rm),
                            )?;
                            dd(
                                &mut h,
                                "hypot",
                                &[&bs, &operand(a)],
                                &[],
                                p,
                                rm,
                                &bs.hypot(&operand(a), p, rm),
                            )?;
                            count += 2;
                        }
                    }
                }
            }
        }
        assert_eq!(count, OPERANDS * 4 * 4 * 2 * 2 * 2);
        println!("hypot scaling digest = {}", h.finish());
        // 32-bit: generated by logs/dg/gen32.py (mpmath values, Word=u32 layout; 64-bit control reproduces the pin above)
        #[cfg(target_pointer_width = "64")]
        const WANT: u64 = 17_199_362_010_888_151_510;
        #[cfg(target_pointer_width = "32")]
        const WANT: u64 = 4_048_660_911_403_563_556;
        assert_eq!(h.finish(), WANT);
        Ok(())
    }

    /// Layer 0a pin for decimal parsing (`conv_from_dec`): short and long digit strings, positive and negative
    /// exponents (including ones that need several powers-of-ten steps), and near-tie inputs that make the
    /// first working precision fail to round and force a retry.
    #[test]
    fn decimal_parse_grid_digest() -> Result<(), Box<dyn std::error::Error>> {
        let mut cc = crate::Consts::new()?;
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        let inputs = [
            "0.1",
            "-2.5e3",
            "123456789012345678901234567890e-250",
            "7.000000000000000000000000000001e500",
            "1e-400",
            "-9007199254740993.00000000000000000000001",
            "9007199254740993.00000000000000000000001",
            "9007199254740993.000000000000000000000000000000000000001",
            "4503599627370497.5000000000000000000000000000000000000000001",
            "0.30000000000000000000000000000000000000000000000000000000000000000000000000000000000001",
            "8.5e-1000",
            "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899e7",
        ];
        for s in inputs {
            for p in [53usize, 64, 128, 500] {
                for rm in
                    [RoundingMode::ToEven, RoundingMode::Up, RoundingMode::Down, RoundingMode::None]
                {
                    digest(
                        &ExactNum::parse(s, crate::Radix::Dec, p, rm, &mut cc),
                        &mut h,
                    )?;
                    count += 1;
                }
            }
        }
        assert_eq!(count, 12 * 4 * 4);
        println!("decimal parse digest = {}", h.finish());
        assert_eq!(h.finish(), 13_825_435_214_421_359_920);
        Ok(())
    }

    /// Layer 0a pin for `atanh` on zeros (positive, negative, inexact from an underflow, of several widths) and
    /// on the whole operand grid: value, sign, inexact flags and NaN kind all enter the digest.
    #[test]
    fn atanh_zero_and_grid_digest() -> Result<(), Box<dyn std::error::Error>> {
        let mut cc = crate::Consts::new()?;
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        let zeros = || {
            [
                ExactNum::new(64),
                ExactNum::new(64).neg(),
                ExactNum::new(256),
                ExactNum::new(256).neg(),
                ExactNum::from_u8(1, 64).ldexp(crate::EXPONENT_MIN, 64, RoundingMode::ToZero),
                ExactNum::from_u8(1, 64)
                    .ldexp(crate::EXPONENT_MIN, 64, RoundingMode::ToZero)
                    .neg(),
                ExactNum::from_u8(1, 192).ldexp(crate::EXPONENT_MIN, 192, RoundingMode::ToZero),
            ]
        };
        for i in 0..zeros().len() {
            for p in [53usize, 64, 128, 300] {
                for rm in [RoundingMode::ToEven, RoundingMode::Down] {
                    let z = zeros();
                    dd(
                        &mut h,
                        "atanh",
                        &[&zeros()[i]],
                        &[],
                        p,
                        rm,
                        &z[i].atanh(p, rm, &mut cc),
                    )?;
                    count += 1;
                }
            }
        }
        for i in 0..OPERANDS {
            for p in [64usize, 200] {
                for rm in [RoundingMode::ToEven, RoundingMode::Down] {
                    dd(
                        &mut h,
                        "atanh",
                        &[&operand(i)],
                        &[],
                        p,
                        rm,
                        &operand(i).atanh(p, rm, &mut cc),
                    )?;
                    dd(
                        &mut h,
                        "atanh",
                        &[&operand(i).ldexp(-3, 256, rm)],
                        &[],
                        p,
                        rm,
                        &operand(i).ldexp(-3, 256, rm).atanh(p, rm, &mut cc),
                    )?;
                    count += 2;
                }
            }
        }
        assert_eq!(count, 7 * 4 * 2 + OPERANDS * 2 * 2 * 2);
        println!("atanh zero digest = {}", h.finish());
        // 32-bit: generated by logs/dg/gen32.py (mpmath values, Word=u32 layout; 64-bit control reproduces the pin above)
        #[cfg(target_pointer_width = "64")]
        const WANT: u64 = 8_560_969_362_905_363_186;
        #[cfg(target_pointer_width = "32")]
        const WANT: u64 = 5_803_567_897_899_187_373;
        assert_eq!(h.finish(), WANT);
        Ok(())
    }

    /// Layer 0a pin for the Jacobi functions: all thirteen, at `u` = 0 and several nonzero (and special) operands,
    /// for the three parameter shapes the kernel special-cases (`m` = 0, `m` = 1, and the general AGM path), and
    /// for an out-of-range `m` (the NaN / error kind enters the digest through the repr).
    #[test]
    fn jacobi_all_functions_digest() -> Result<(), Box<dyn std::error::Error>> {
        type J = fn(&ExactNum, &ExactNum, usize, RoundingMode, &mut crate::Consts) -> ExactNum;
        const JNAMES: [&str; 13] =
            ["am", "sn", "cn", "dn", "cd", "ns", "nc", "nd", "sc", "sd", "cs", "ds", "dc"];
        let fns: [J; 13] = [
            ExactNum::jacobi_am,
            ExactNum::jacobi_sn,
            ExactNum::jacobi_cn,
            ExactNum::jacobi_dn,
            ExactNum::jacobi_cd,
            ExactNum::jacobi_ns,
            ExactNum::jacobi_nc,
            ExactNum::jacobi_nd,
            ExactNum::jacobi_sc,
            ExactNum::jacobi_sd,
            ExactNum::jacobi_cs,
            ExactNum::jacobi_ds,
            ExactNum::jacobi_dc,
        ];
        let mut cc = crate::Consts::new()?;
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        let ms = || {
            [
                ExactNum::new(64),
                ExactNum::from_u8(1, 64),
                ExactNum::from_u8(1, 64).ldexp(-1, 64, RoundingMode::None),
                ExactNum::from_u8(3, 64).ldexp(-3, 64, RoundingMode::None),
                ExactNum::from_u8(3, 64),
            ]
        };
        for (fi, f) in fns.into_iter().enumerate() {
            for i in [0usize, 1, 2, 3, 5, 7, 21] {
                for mi in 0..5 {
                    for p in [64usize, 128] {
                        for rm in [RoundingMode::ToEven, RoundingMode::Down] {
                            let m = ms();
                            dd(
                                &mut h,
                                JNAMES[fi],
                                &[&operand(i), &m[mi]],
                                &[],
                                p,
                                rm,
                                &f(&operand(i), &m[mi], p, rm, &mut cc),
                            )?;
                            count += 1;
                        }
                    }
                }
            }
        }
        assert_eq!(count, 13 * 7 * 5 * 2 * 2);
        // oracle: am(u|0) = u. A 192-bit inexact `u` made the working-precision check overflow before the fix.
        let u = operand(5);
        let zero_m = ExactNum::new(64);
        let am = u.jacobi_am(&zero_m, 300, RoundingMode::ToEven, &mut cc);
        assert!(ExactNum::sub_ref(&am, &u, 300, RoundingMode::ToEven).is_zero());
        println!("jacobi digest = {}", h.finish());
        // Re-pinned (was 8_344_934_497_909_943_058): `am` dropped the \(2\pi\) multiples of the \(4K\) period
        // reduction (am(-7|1/2) gave +0.4105, mpmath -5.8727); 8 `am` entries (u = -7, m = 1/2, 3/8) changed.
        // Re-pinned (was 638_880_846_716_710_340): `m sn^2` underflowing to 0 made dn/nd of the subnormal
        // operand (m = 1/2, 3/8) an exact 1. 16 dn/nd entries changed: `ToEven` is now 1 flagged inexact
        // (mpmath: 1 - dn = m u^2/2 > 0), `Down` is now `NaN(PrecisionRetryExhausted)` like `cn` there.
        // Earlier re-pin (was 13_695_875_014_094_095_401): 23 entries x 2 precisions (all at the tiny operand 21, `ToEven`)
        // were `NaN(PrecisionRetryExhausted)`, a wrong answer. `am`, `sn`, `sc`, `sd` of `u = 3 * 2^(-2^31 - 19)`
        // are `u` and the other nine functions are 1, because `f(u|m) = t (1 - O(u^2))` is far closer to the
        // representable `t` than half an ulp, and a value just below a representable number rounds to it
        // to nearest; the roundability check had treated that as undecidable. Oracle: mpmath, see the CHANGELOG.
        // 32-bit: generated by logs/dg/gen32.py (mpmath values, Word=u32 layout; 64-bit control reproduces the pin above)
        #[cfg(target_pointer_width = "64")]
        const WANT: u64 = 8_659_318_462_419_821_720;
        #[cfg(target_pointer_width = "32")]
        const WANT: u64 = 6_139_828_297_236_387_860;
        assert_eq!(h.finish(), WANT);
        Ok(())
    }

    /// `f(u|m)` for `u = 3 * 2^(-2^31 - 19)` (a subnormal operand) is `t (1 - O(u^2))` with `t = u` for `am`, `sn`,
    /// `sc`, `sd` and `t = 1` for `cn`, `dn`, `cd`, `nc`, `nd`, `dc` (series in `u`; mpmath at a surrogate `u`
    /// confirms the constants). The deviation is below `2^-(2^32)`, so to nearest the result is `t` itself,
    /// not an error. `ToEven` used to return `NaN(PrecisionRetryExhausted)` for most of these: a value just
    /// below a representable number is not a hard case to round to nearest.
    #[test]
    fn jacobi_of_a_tiny_argument_rounds_to_the_leading_term(
    ) -> Result<(), Box<dyn std::error::Error>> {
        type J = fn(&ExactNum, &ExactNum, usize, RoundingMode, &mut crate::Consts) -> ExactNum;
        let leading_u: [J; 4] = [
            ExactNum::jacobi_am,
            ExactNum::jacobi_sn,
            ExactNum::jacobi_sc,
            ExactNum::jacobi_sd,
        ];
        let leading_one: [J; 6] = [
            ExactNum::jacobi_cn,
            ExactNum::jacobi_dn,
            ExactNum::jacobi_cd,
            ExactNum::jacobi_nc,
            ExactNum::jacobi_nd,
            ExactNum::jacobi_dc,
        ];
        let mut cc = crate::Consts::new()?;
        // m = 0, 1, 1/2, 3/16 (am(u|0) = u exactly; the others by the series)
        let ms = [
            ExactNum::new(64),
            ExactNum::from_u8(1, 64),
            ExactNum::from_u8(1, 64).ldexp(-1, 64, RoundingMode::None),
            ExactNum::from_u8(3, 64).ldexp(-4, 64, RoundingMode::None),
        ];
        for p in [64usize, 128] {
            for m in &ms {
                for f in leading_u {
                    let got = f(&operand(21), m, p, RoundingMode::ToEven, &mut cc);
                    assert_eq!(
                        got.compare(&operand(21)),
                        Some(0),
                        "u p={p} m={m:?}: {got:?}"
                    );
                }
                for f in leading_one {
                    let got = f(&operand(21), m, p, RoundingMode::ToEven, &mut cc);
                    assert_eq!(
                        got.compare(&ExactNum::from_u8(1, 64)),
                        Some(0),
                        "1 p={p} m={m:?}: {got:?}"
                    );
                }
            }
        }
        Ok(())
    }

    /// `dn(u|m) = 1 - m u^2/2 + O(u^4)` and `nd = 1/dn` for the subnormal `u = operand(21)` and `m = 1/2, 3/8`:
    /// the value is strictly below (`nd`: above) 1, so to nearest it is 1 and inexact (mpmath at a surrogate
    /// `u = 2^-200` gives `(1 - dn) / (m u^2/2) = 1` and `(1/dn - 1) / (m u^2/2) = 1` to 60 digits). `sn^2` underflows to an exact 0, which used to make
    /// `1 - m sn^2` an exact 1 with the inexact flag off.
    #[test]
    fn jacobi_dn_nd_of_a_subnormal_argument_are_inexact() -> Result<(), Box<dyn std::error::Error>>
    {
        type J = fn(&ExactNum, &ExactNum, usize, RoundingMode, &mut crate::Consts) -> ExactNum;
        let fns: [(&str, J); 2] = [("dn", ExactNum::jacobi_dn), ("nd", ExactNum::jacobi_nd)];
        let mut cc = crate::Consts::new()?;
        let ms = [
            ExactNum::from_u8(1, 64).ldexp(-1, 64, RoundingMode::None),
            ExactNum::from_u8(3, 64).ldexp(-3, 64, RoundingMode::None),
        ];
        for p in [64usize, 128] {
            for m in &ms {
                for (name, f) in fns {
                    let got = f(&operand(21), m, p, RoundingMode::ToEven, &mut cc);
                    assert_eq!(
                        got.compare(&ExactNum::from_u8(1, 64)),
                        Some(0),
                        "{name} p={p} m={m:?}: {got:?}"
                    );
                    assert!(
                        got.inexact(),
                        "{name} p={p} m={m:?} must be inexact: {got:?}"
                    );
                }
            }
        }
        Ok(())
    }

    /// Layer 0a pin for the four Airy functions: zero (closed form), small arguments (Taylor pair), large
    /// positive and negative arguments (asymptotic expansions), and every special operand.
    #[test]
    fn airy_grid_digest() -> Result<(), Box<dyn std::error::Error>> {
        let mut cc = crate::Consts::new()?;
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        let extra = |k: usize| match k {
            0 => ExactNum::from_u8(3, 64).ldexp(-2, 64, RoundingMode::None),
            1 => ExactNum::from_u8(5, 64).neg(),
            2 => ExactNum::from_u8(30, 64),
            3 => ExactNum::from_u8(30, 64).neg(),
            _ => ExactNum::from_u8(9, 64).ldexp(-1, 64, RoundingMode::None),
        };
        for p in [64usize, 128] {
            for rm in [RoundingMode::ToEven, RoundingMode::Down] {
                for i in (0..11).chain(21..23) {
                    let x = operand(i);
                    digest(&x.ai(p, rm, &mut cc), &mut h)?;
                    digest(&x.bi(p, rm, &mut cc), &mut h)?;
                    digest(&x.ai_prime(p, rm, &mut cc), &mut h)?;
                    digest(&x.bi_prime(p, rm, &mut cc), &mut h)?;
                    count += 4;
                }
                for k in 0..5 {
                    let x = extra(k);
                    digest(&x.ai(p, rm, &mut cc), &mut h)?;
                    digest(&x.bi(p, rm, &mut cc), &mut h)?;
                    digest(&x.ai_prime(p, rm, &mut cc), &mut h)?;
                    digest(&x.bi_prime(p, rm, &mut cc), &mut h)?;
                    count += 4;
                }
            }
        }
        assert_eq!(count, 2 * 2 * (13 + 5) * 4);
        println!("airy digest = {}", h.finish());
        assert_eq!(h.finish(), 1_394_724_982_563_781_065);
        Ok(())
    }

    /// Layer 0a pin for `powi` with exponents long enough to use the windowed square-and-multiply table of odd
    /// powers (and the plain loop for short ones), over positive and negative bases and all three modes.
    #[test]
    fn powi_windowed_digest() -> Result<(), Box<dyn std::error::Error>> {
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        for i in [2usize, 3, 4, 5, 6] {
            for n in [7usize, 255, 1000, 65537, 1_000_003] {
                for p in [64usize, 200] {
                    for rm in [RoundingMode::ToEven, RoundingMode::Up, RoundingMode::Down] {
                        dd(
                            &mut h,
                            "powi",
                            &[&operand(i)],
                            &[i128::try_from(n)?],
                            p,
                            rm,
                            &operand(i).powi(n, p, rm),
                        )?;
                        count += 1;
                    }
                }
            }
        }
        assert_eq!(count, 5 * 5 * 2 * 3);
        println!("powi windowed digest = {}", h.finish());
        // 32-bit: generated by logs/dg/gen32.py (mpmath values, Word=u32 layout; 64-bit control reproduces the pin above)
        #[cfg(target_pointer_width = "64")]
        const WANT: u64 = 12_976_918_398_664_607_263;
        #[cfg(target_pointer_width = "32")]
        const WANT: u64 = 14_780_506_644_587_673_446;
        assert_eq!(h.finish(), WANT);
        Ok(())
    }

    /// Layer 1a pin: the entry points whose Ziv / precision-retry loops keep a working copy of the operand
    /// across iterations (`ln`, `log` with a base, `log1p`, `sin_cos`, `sinh_cosh`, `nth_root`), over the whole
    /// operand grid (specials included, so error kinds and NaN kinds enter the digest) and several
    /// precisions and modes. Recorded on the code that still copies, before the owned-operand wrapper.
    #[test]
    fn ziv_entry_points_digest() -> Result<(), Box<dyn std::error::Error>> {
        let mut cc = crate::Consts::new()?;
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        for i in (0..11).chain(21..23) {
            for p in [64usize, 200] {
                for rm in [RoundingMode::ToEven, RoundingMode::Down] {
                    dd(
                        &mut h,
                        "ln",
                        &[&operand(i)],
                        &[],
                        p,
                        rm,
                        &operand(i).ln(p, rm, &mut cc),
                    )?;
                    dd(
                        &mut h,
                        "log1p",
                        &[&operand(i)],
                        &[],
                        p,
                        rm,
                        &operand(i).log1p(p, rm, &mut cc),
                    )?;
                    let (s, c) = operand(i).sin_cos(p, rm, &mut cc);
                    dd(&mut h, "sin", &[&operand(i)], &[], p, rm, &s)?;
                    dd(&mut h, "cos", &[&operand(i)], &[], p, rm, &c)?;
                    let (s, c) = operand(i).sinh_cosh(p, rm, &mut cc);
                    dd(&mut h, "sinh", &[&operand(i)], &[], p, rm, &s)?;
                    dd(&mut h, "cosh", &[&operand(i)], &[], p, rm, &c)?;
                    dd(
                        &mut h,
                        "nth_root",
                        &[&operand(i)],
                        &[3],
                        p,
                        rm,
                        &operand(i).nth_root(3, p, rm),
                    )?;
                    dd(
                        &mut h,
                        "nth_root",
                        &[&operand(i)],
                        &[5],
                        p,
                        rm,
                        &operand(i).nth_root(5, p, rm),
                    )?;
                    count += 8;
                    for b in [2usize, 4, 6] {
                        dd(
                            &mut h,
                            "logb",
                            &[&operand(i), &operand(b)],
                            &[],
                            p,
                            rm,
                            &operand(i).log(&operand(b), p, rm, &mut cc),
                        )?;
                        count += 1;
                    }
                }
            }
        }
        assert_eq!(count, 13 * 2 * 2 * 11);
        println!("ziv entry points digest = {}", h.finish());
        // Re-recorded: the old value (2_735_357_469_467_868_262) pinned a wrong `log1p` result for the positive
        // subnormal operand 21 with `Down` (the argument itself instead of the argument minus one grid unit);
        // only those two entries (p = 64 and p = 200) differ. See the CHANGELOG.
        // Re-recorded again for `nth_root` of a negative number with `Up`/`Down` (they were swapped for odd n >= 5):
        // the old value was 6_378_838_253_543_692_452; only `nth_root(5)` entries of negative operands differ,
        // each confirmed by an exact integer root. See the CHANGELOG.
        // 32-bit: generated by logs/dg/gen32.py (mpmath values, Word=u32 layout; 64-bit control reproduces the pin above)
        #[cfg(target_pointer_width = "64")]
        const WANT: u64 = 15_733_382_221_161_312_340;
        #[cfg(target_pointer_width = "32")]
        const WANT: u64 = 642_497_529_214_412_721;
        assert_eq!(h.finish(), WANT);
        Ok(())
    }

    /// Layer 1a pin for the exact-result branch of `log(x, base)` (the raw working value is only close to an
    /// exact result, so the nearest candidate is checked with integer powers): known exact results are
    /// asserted, and every case enters the digest.
    #[test]
    fn log_exact_results_digest() -> Result<(), Box<dyn std::error::Error>> {
        let mut cc = crate::Consts::new()?;
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        // (x, base, expected exact result or None when the result is not exact)
        let cases: [(ExactNum, ExactNum, Option<ExactNum>); 8] = [
            (
                ExactNum::from_u8(8, 64),
                ExactNum::from_u8(2, 64),
                Some(ExactNum::from_u8(3, 64)),
            ),
            (
                ExactNum::from_u8(100, 64),
                ExactNum::from_u8(10, 64),
                Some(ExactNum::from_u8(2, 64)),
            ),
            (
                ExactNum::from_u8(81, 64),
                ExactNum::from_u8(3, 64),
                Some(ExactNum::from_u8(4, 64)),
            ),
            (
                ExactNum::from_u8(125, 64),
                ExactNum::from_u8(5, 64),
                Some(ExactNum::from_u8(3, 64)),
            ),
            (
                ExactNum::from_u8(1, 64).ldexp(-10, 64, RoundingMode::None),
                ExactNum::from_u8(2, 64),
                Some(ExactNum::from_u8(10, 64).neg()),
            ),
            (ExactNum::from_u8(10, 64), ExactNum::from_u8(2, 64), None),
            (ExactNum::from_u8(7, 64), ExactNum::from_u8(3, 64), None),
            (ExactNum::from_u8(2, 64), ExactNum::from_u8(8, 64), None),
        ];
        for (x, b, want) in cases {
            for p in [53usize, 64, 128, 200] {
                for rm in [RoundingMode::ToEven, RoundingMode::Down, RoundingMode::Up] {
                    let got = x.log(&b, p, rm, &mut cc);
                    if let Some(w) = &want {
                        assert!(ExactNum::sub_ref(&got, w, 64, RoundingMode::ToEven).is_zero());
                    }
                    digest(&got, &mut h)?;
                    count += 1;
                }
            }
        }
        assert_eq!(count, 8 * 4 * 3);
        println!("log exact digest = {}", h.finish());
        assert_eq!(h.finish(), 11_246_996_610_564_981_343);
        Ok(())
    }

    /// Loop-conversion pin for `ln`, `log2`, `log10` (the per-pass working copy becomes `ziv_pass_operand`): the 23
    /// grid operands (NaN kinds, +-Inf, +-0, the subnormal band, the exponent extremes), exact powers 2^k and 10^n over
    /// a wide exponent range, 1 +- 2^-j near the cancellation point, and inexact-flagged inputs (including an
    /// inexact-flagged exact power of ten for the log10 exact branch), at 7 precisions in all 6 rounding modes.
    /// Recorded on the code that still clones the input.
    #[test]
    fn ln_log2_log10_conversion_digest() -> Result<(), Box<dyn std::error::Error>> {
        let mut cc = crate::Consts::new()?;
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        let mut xs: Vec<ExactNum> = (0..PAIR_OPERANDS).map(operand).collect();
        for k in [1i64, 7, 64, 1000, 100_000, 1_000_000_000, -1, -64, -100_000, -1_000_000_000] {
            xs.push(ExactNum::from_u8(1, 64).ldexp(k, 64, RoundingMode::None));
            xs.push(ExactNum::from_u8(3, 64).ldexp(k, 64, RoundingMode::None));
        }
        for n in [1usize, 2, 3, 19, 20, 40, 80] {
            xs.push(ExactNum::from_u8(10, 320).powi(n, 320, RoundingMode::None));
        }
        for j in [1i64, 10, 63, 64, 65, 200] {
            let d = ExactNum::from_u8(1, 64).ldexp(-j, 64, RoundingMode::None);
            xs.push(ExactNum::add_ref(
                &ExactNum::from_u8(1, 320),
                &d,
                320,
                RoundingMode::None,
            ));
            xs.push(ExactNum::sub_ref(
                &ExactNum::from_u8(1, 320),
                &d,
                320,
                RoundingMode::None,
            ));
        }
        let mut hundred = ExactNum::from_u8(100, 64);
        hundred.set_inexact(true);
        xs.push(hundred);
        let mut tenth_ish =
            ExactNum::from_u8(1, 64).div(&ExactNum::from_u8(10, 64), 64, RoundingMode::ToEven);
        tenth_ish.set_inexact(true);
        xs.push(tenth_ish);
        let nx = xs.len();
        for x in &xs {
            for p in [1usize, 24, 53, 64, 100, 200, 1024] {
                for rm in RMS {
                    dd(&mut h, "ln", &[x], &[], p, rm, &x.ln(p, rm, &mut cc))?;
                    dd(&mut h, "log2", &[x], &[], p, rm, &x.log2(p, rm, &mut cc))?;
                    dd(&mut h, "log10", &[x], &[], p, rm, &x.log10(p, rm, &mut cc))?;
                    count += 3;
                }
            }
        }
        assert_eq!(count, nx * 7 * 6 * 3);
        println!(
            "ln/log2/log10 conversion digest = {} over {} entries",
            h.finish(),
            count
        );
        assert_eq!(h.finish(), 2_640_369_975_324_367_864);
        Ok(())
    }

    /// Layer 0a pin for the reciprocal (Newton step) and `nth_root` (exact-candidate check) kernels.
    #[test]
    fn reciprocal_nth_root_grid_digest() -> Result<(), Box<dyn std::error::Error>> {
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        for i in 0..PAIR_OPERANDS {
            for p in [1usize, 53, 64, 128, 200, 320] {
                for rm in RMS {
                    dd(
                        &mut h,
                        "reciprocal",
                        &[&operand(i)],
                        &[],
                        p,
                        rm,
                        &operand(i).reciprocal(p, rm),
                    )?;
                    count += 1;
                    for n in [2usize, 3, 4, 5, 7] {
                        dd(
                            &mut h,
                            "nth_root",
                            &[&operand(i)],
                            &[i128::try_from(n)?],
                            p,
                            rm,
                            &operand(i).nth_root(n, p, rm),
                        )?;
                        count += 1;
                    }
                }
            }
        }
        assert_eq!(count, PAIR_OPERANDS * 6 * 6 * 6);
        println!("reciprocal nth_root grid digest = {}", h.finish());
        // Re-recorded for `nth_root` of a negative number with `Up`/`Down` (swapped for odd n >= 5; old value
        // 1_369_146_293_907_659_395): only the n = 5 and 7 entries of negative operands in those modes differ,
        // each confirmed by an exact integer root. See the CHANGELOG.
        assert_eq!(h.finish(), 2_411_520_695_892_384_189);
        Ok(())
    }

    /// Field-for-field agreement of `sqr` with `mul_ref(&x, &x)` (not just a hash), including the NaN kind of
    /// a NaN operand, the sign of a zero result in every mode, and that `x` is left untouched.
    #[test]
    fn sqr_equals_mul_ref_of_the_same_operand() -> Result<(), Box<dyn std::error::Error>> {
        for i in 0..PAIR_OPERANDS {
            for p in [1usize, 24, 53, 64, 65, 128, 192, 200, 320] {
                for rm in RMS {
                    let x = operand(i);
                    let before = value_repr(&x)?;
                    let want = ExactNum::mul_ref(&x, &x, p, rm);
                    let got = x.sqr(p, rm);
                    assert_eq!(value_repr(&got)?, value_repr(&want)?, "i={i} p={p} {rm:?}");
                    assert_eq!(got.inexact(), want.inexact(), "i={i} p={p} {rm:?}");
                    assert_eq!(
                        got.inexact_accumulated(),
                        want.inexact_accumulated(),
                        "i={i} p={p} {rm:?}"
                    );
                    assert_eq!(got.err(), want.err(), "NaN kind, i={i} p={p} {rm:?}");
                    assert_eq!(
                        value_repr(&x)?,
                        before,
                        "operand changed, i={i} p={p} {rm:?}"
                    );
                }
            }
        }
        Ok(())
    }

    /// The `*` and `/` operators on every combination of owned and borrowed operands that remains after
    /// Rule R (`&a * &b` and `&a / &b` do not exist). Pinned on the code that still had them, hashing the
    /// same three forms.
    #[test]
    fn operator_mul_div_grid_digest() -> Result<(), Box<dyn std::error::Error>> {
        let mut h = DefaultHasher::new();
        let mut count = 0usize;
        for i in 0..PAIR_OPERANDS {
            for j in 0..PAIR_OPERANDS {
                dd(
                    &mut h,
                    "op*",
                    &[&operand(i), &operand(j)],
                    &[],
                    0,
                    RoundingMode::ToEven,
                    &(&operand(i) * operand(j)),
                )?;
                dd(
                    &mut h,
                    "op*",
                    &[&operand(i), &operand(j)],
                    &[],
                    0,
                    RoundingMode::ToEven,
                    &(operand(i) * &operand(j)),
                )?;
                dd(
                    &mut h,
                    "op*",
                    &[&operand(i), &operand(j)],
                    &[],
                    0,
                    RoundingMode::ToEven,
                    &(operand(i) * operand(j)),
                )?;
                dd(
                    &mut h,
                    "op/",
                    &[&operand(i), &operand(j)],
                    &[],
                    0,
                    RoundingMode::ToEven,
                    &(&operand(i) / operand(j)),
                )?;
                dd(
                    &mut h,
                    "op/",
                    &[&operand(i), &operand(j)],
                    &[],
                    0,
                    RoundingMode::ToEven,
                    &(operand(i) / &operand(j)),
                )?;
                dd(
                    &mut h,
                    "op/",
                    &[&operand(i), &operand(j)],
                    &[],
                    0,
                    RoundingMode::ToEven,
                    &(operand(i) / operand(j)),
                )?;
                count += 6;
            }
        }
        assert_eq!(count, PAIR_OPERANDS * PAIR_OPERANDS * 6);
        println!("operator mul/div grid digest = {}", h.finish());
        assert_eq!(h.finish(), 11_238_293_337_447_225_181);
        Ok(())
    }

    /// All seven rounding modes, `ToOdd` included (the pinned grids above use six).
    const ALL_RMS: [RoundingMode; 7] = [
        RoundingMode::None,
        RoundingMode::ToEven,
        RoundingMode::ToZero,
        RoundingMode::Up,
        RoundingMode::Down,
        RoundingMode::FromZero,
        RoundingMode::ToOdd,
    ];

    /// Two results are the same bit for bit: the value as stored (sign, exponent, words, NaN kind), and both inexact flags.
    fn same(a: &ExactNum, b: &ExactNum) -> Result<bool, Box<dyn std::error::Error>> {
        Ok(value_repr(a)? == value_repr(b)?
            && a.inexact() == b.inexact()
            && a.inexact_accumulated() == b.inexact_accumulated())
    }

    /// `add_ref`, `sub_ref`, `mul_ref` and `div_ref` give the same result as `add`, `sub`, `mul` and `div` on a fresh
    /// owned operand, bit for bit, on the pinned pair grid (the 23 operands: NaN of both kinds, +-Inf, +-0, inexact thirds,
    /// the subnormal band, the exponent extremes) at the 7 precisions of the digest grids, in all 7 rounding modes.
    #[test]
    fn ref_forms_equal_owned_forms_bit_for_bit() -> Result<(), Box<dyn std::error::Error>> {
        let mut count = 0usize;
        for i in 0..PAIR_OPERANDS {
            for j in 0..PAIR_OPERANDS {
                for p in [1usize, 53, 64, 65, 128, 192, 200] {
                    for rm in ALL_RMS {
                        let (a, b) = (operand(i), operand(j));
                        let pairs = [
                            (
                                "add",
                                ExactNum::add_ref(&a, &b, p, rm),
                                operand(i).add(&b, p, rm),
                            ),
                            (
                                "sub",
                                ExactNum::sub_ref(&a, &b, p, rm),
                                operand(i).sub(&b, p, rm),
                            ),
                            (
                                "mul",
                                ExactNum::mul_ref(&a, &b, p, rm),
                                operand(i).mul(&b, p, rm),
                            ),
                            (
                                "div",
                                ExactNum::div_ref(&a, &b, p, rm),
                                operand(i).div(&b, p, rm),
                            ),
                        ];
                        for (name, by_ref, owned) in &pairs {
                            assert!(
                                same(by_ref, owned)?,
                                "{name}_ref differs from {name}: operands {i}, {j}, p = {p}, {rm:?}: {by_ref:?} vs {owned:?}"
                            );
                            count += 1;
                        }
                        // the operands are untouched by the ref forms
                        assert!(same(&a, &operand(i))? && same(&b, &operand(j))?);
                    }
                }
            }
        }
        assert_eq!(count, PAIR_OPERANDS * PAIR_OPERANDS * 7 * 7 * 4);
        Ok(())
    }

    /// The same operand passed twice by reference (`add_ref(&x, &x)`), as the squares and `x / x` sites do.
    #[test]
    fn ref_forms_with_one_operand_twice_equal_the_owned_forms(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut count = 0usize;
        for i in 0..PAIR_OPERANDS {
            for p in [1usize, 53, 64, 65, 128, 192, 200] {
                for rm in ALL_RMS {
                    let x = operand(i);
                    let pairs = [
                        (
                            "add",
                            ExactNum::add_ref(&x, &x, p, rm),
                            operand(i).add(&operand(i), p, rm),
                        ),
                        (
                            "sub",
                            ExactNum::sub_ref(&x, &x, p, rm),
                            operand(i).sub(&operand(i), p, rm),
                        ),
                        (
                            "mul",
                            ExactNum::mul_ref(&x, &x, p, rm),
                            operand(i).mul(&operand(i), p, rm),
                        ),
                        (
                            "div",
                            ExactNum::div_ref(&x, &x, p, rm),
                            operand(i).div(&operand(i), p, rm),
                        ),
                    ];
                    for (name, by_ref, owned) in &pairs {
                        assert!(
                            same(by_ref, owned)?,
                            "{name}_ref(&x, &x) differs from {name}: operand {i}, p = {p}, {rm:?}: {by_ref:?} vs {owned:?}"
                        );
                        count += 1;
                    }
                }
            }
        }
        assert_eq!(count, PAIR_OPERANDS * 7 * 7 * 4);
        Ok(())
    }

    /// `ExactComplex::add_ref` / `sub_ref` against `add` / `sub`, on the complex grid (8 parts per component, 3 precisions),
    /// all 7 rounding modes, both parts and both inexact flags of each part.
    #[test]
    fn complex_ref_forms_equal_owned_forms_bit_for_bit() -> Result<(), Box<dyn std::error::Error>> {
        use crate::ExactComplex;
        let parts = [0usize, 1, 3, 4, 5, 7, 9, 21];
        let mut count = 0usize;
        for a in parts {
            for b in parts {
                for c in parts {
                    for d in parts {
                        for p in [1usize, 64, 128] {
                            for rm in ALL_RMS {
                                let z1 = ExactComplex::new(operand(a), operand(b));
                                let z2 = ExactComplex::new(operand(c), operand(d));
                                let own1 = || ExactComplex::new(operand(a), operand(b));
                                let pairs = [
                                    (
                                        "add",
                                        ExactComplex::add_ref(&z1, &z2, p, rm),
                                        own1().add(&z2, p, rm),
                                    ),
                                    (
                                        "sub",
                                        ExactComplex::sub_ref(&z1, &z2, p, rm),
                                        own1().sub(&z2, p, rm),
                                    ),
                                ];
                                for (name, by_ref, owned) in &pairs {
                                    assert!(
                                        same(by_ref.re(), owned.re())? && same(by_ref.im(), owned.im())?,
                                        "complex {name}_ref differs from {name}: parts {a}, {b}, {c}, {d}, p = {p}, {rm:?}"
                                    );
                                    count += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(count, 8 * 8 * 8 * 8 * 3 * 7 * 2);
        Ok(())
    }
}

#[cfg(test)]
mod to_odd_tests {
    use crate::{ExactNum, RoundingMode};

    /// `ToOdd` is nearest with ties to odd, not sticky round-to-odd. Hand-derived values at `p = 4`:
    /// 17/16 = 1.0001b is an exact tie (1.000b = 1 for `ToEven`, 1.001b = 9/8 for `ToOdd`); 33/32 = 1.00001b is below
    /// the tie (1.000b under both nearest modes; a sticky round-to-odd would give 1.001b).
    #[test]
    fn to_odd_is_nearest_with_ties_to_odd() {
        let p = 4;
        let q = |a: u8, b: u8, rm| ExactNum::from_u8(a, 64).div(&ExactNum::from_u8(b, 64), p, rm);
        let one = ExactNum::from_u8(1, 64);
        let nine_eighths =
            ExactNum::from_u8(9, 64).div(&ExactNum::from_u8(8, 64), 64, RoundingMode::ToEven);

        assert_eq!(q(17, 16, RoundingMode::ToEven).compare(&one), Some(0));
        assert_eq!(
            q(17, 16, RoundingMode::ToOdd).compare(&nine_eighths),
            Some(0)
        );
        assert_eq!(q(33, 32, RoundingMode::ToEven).compare(&one), Some(0));
        assert_eq!(q(33, 32, RoundingMode::ToOdd).compare(&one), Some(0));
        // 33/32 is just above 1.000b, so `Up` gives 1.001b: the value a sticky round-to-odd would also give.
        assert_eq!(q(33, 32, RoundingMode::Up).compare(&nine_eighths), Some(0));
    }
}

#[cfg(test)]
mod subnormal_copy_tests {
    use crate::mantissa::CLONE_CALLS;
    use crate::{ExactNum, Exponent, RoundingMode};

    fn clones() -> usize {
        CLONE_CALLS.with(|c| c.get())
    }

    /// `3 * 2^(EXPONENT_MIN + 1) * 2^-20`: a subnormal (its mantissa has leading zero bits).
    fn subnormal(p: usize) -> ExactNum {
        ExactNum::from_u8(3, p)
            .ldexp(Exponent::MIN + 1, p, RoundingMode::ToEven)
            .ldexp(-20, p, RoundingMode::ToEven)
    }

    /// No operation copies an operand, subnormal or not: the number of `Mantissa::clone` calls made while
    /// an operation runs is zero, for the four borrowed-borrowed forms, `sqr`, and the owned forms, with the
    /// subnormal operand on either side. (A subnormal operand used to be cloned once by `normalize`; a subnormal divisor is shifted into a kernel working buffer, which is not a `Mantissa::clone`.)
    #[test]
    fn subnormal_operands_are_not_copied() {
        let p = 192;
        let rm = RoundingMode::ToEven;
        let two = ExactNum::from_u8(2, p);
        let sub = subnormal(p);

        let before = clones();
        let _ = ExactNum::add_ref(&sub, &two, p, rm);
        let _ = ExactNum::add_ref(&two, &sub, p, rm);
        let _ = ExactNum::sub_ref(&sub, &two, p, rm);
        let _ = ExactNum::sub_ref(&two, &sub, p, rm);
        let _ = ExactNum::mul_ref(&sub, &two, p, rm);
        let _ = ExactNum::mul_ref(&two, &sub, p, rm);
        let _ = ExactNum::div_ref(&sub, &two, p, rm);
        let _ = ExactNum::div_ref(&two, &sub, p, rm);
        let _ = sub.sqr(p, rm);
        assert_eq!(clones() - before, 0, "borrowed forms and sqr");

        let before = clones();
        let _ = subnormal(p).add(&two, p, rm);
        let _ = ExactNum::from_u8(2, p).add(&sub, p, rm);
        let _ = subnormal(p).sub(&two, p, rm);
        let _ = ExactNum::from_u8(2, p).sub(&sub, p, rm);
        let _ = subnormal(p).mul(&two, p, rm);
        let _ = ExactNum::from_u8(2, p).mul(&sub, p, rm);
        let _ = subnormal(p).div(&two, p, rm);
        let _ = ExactNum::from_u8(2, p).div(&sub, p, rm);
        assert_eq!(clones() - before, 0, "owned forms");
    }
}

#[cfg(test)]
mod round_op_flag_tests {
    use crate::{ExactNum, RoundingMode};

    // `round` used to keep the operand's per-operation flag, so `inexact()` reported the
    // operation that made the operand instead of this rounding.
    #[test]
    fn round_reports_its_own_rounding() -> Result<(), crate::Error> {
        let x = ExactNum::from_mag_parts(0xC000_0000_0000_0000, -1, 64)?; // 0.375, exact
        assert!(!x.inexact());
        for rm in [RoundingMode::ToEven, RoundingMode::ToZero, RoundingMode::Up] {
            let y = ExactNum::from_mag_parts(0xC000_0000_0000_0000, -1, 64)?.round(0, rm);
            assert!(y.inexact(), "{rm:?}");
            assert!(y.inexact_accumulated(), "{rm:?}");
            let z = y.round(64, rm); // nothing left to round
            assert!(!z.inexact(), "{rm:?}");
            assert!(z.inexact_accumulated(), "{rm:?}");
        }
        for rm in [RoundingMode::ToEven, RoundingMode::ToZero, RoundingMode::None] {
            // a subnormal operand rounded to an integer is an inexact zero
            let t = ExactNum::from_mag_parts(0xC000_0000_0000_0000, crate::EXPONENT_MIN + 5, 64)?
                .round(0, rm);
            assert!(t.is_zero() && t.inexact(), "{rm:?}");
        }
        let w =
            ExactNum::from_mag_parts(0xC000_0000_0000_0000, 2, 64)?.round(1, RoundingMode::ToEven); // 3.0 exact
        assert!(!w.inexact());
        assert!(!w.inexact_accumulated());
        Ok(())
    }
}

#[cfg(test)]
mod underflow_zero_flag_tests {
    use crate::{Consts, ExactNum, RoundingMode, EXPONENT_MIN};

    // `log(x, base)` rounded at the word size (`round_p(p)`), not at `p`, unlike every other
    // correctly rounded function: at `p = 200` it returned a 256-bit value.
    #[test]
    fn log_with_a_base_rounds_at_the_requested_precision() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let r = ExactNum::from_u8(3, 64).log(
            &ExactNum::from_u8(2, 64),
            200,
            RoundingMode::ToZero,
            &mut cc,
        );
        let on_grid = ExactNum::add_ref(&r, &ExactNum::new(64), 200, RoundingMode::ToZero);
        assert!(r == on_grid && r.inexact());
        Ok(())
    }

    // The exact squaring behind an integer power used to flush a tiny partial product to a fresh
    // exact zero, so an integer power of a tiny subnormal number was an exact zero in
    // every mode, even `Up`, where the result is the smallest subnormal number.
    #[test]
    fn power_of_a_tiny_subnormal_is_not_an_exact_zero() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let x = ExactNum::from_mag_parts(0xC000_0000_0000_0000, EXPONENT_MIN - 17, 64)?;
        let three = ExactNum::from_u8(3, 64);
        let up = x.pow(&three, 53, RoundingMode::Up, &mut cc);
        assert!(up.is_positive() && !up.is_zero());
        assert!(up.inexact());
        let tz = x.pow(&three, 53, RoundingMode::ToZero, &mut cc);
        assert!(tz.is_zero());
        assert!(tz.inexact());
        // the least positive value is the quantum of the result precision, and far below half of it
        // a nearest mode gives zero
        let up128 = x.pow(&three, 128, RoundingMode::Up, &mut cc);
        assert!(up128 == ExactNum::pow2(EXPONENT_MIN - 128, 128)? && up128.inexact());
        let ne = x.pow(&three, 200, RoundingMode::ToEven, &mut cc);
        assert!(ne.is_zero() && ne.inexact());
        // a power `N / 2` takes the exact square root, then the same integer power
        let r = ExactNum::from_mag_parts(0x9000_0000_0000_0000, EXPONENT_MIN - 17, 64)?; // its root is exact
        let half5 = ExactNum::from_mag_parts(0xA000_0000_0000_0000, 2, 64)?; // 2.5
        let up = r.pow(&half5, 53, RoundingMode::Up, &mut cc);
        assert!(up.is_positive() && !up.is_zero() && up.inexact());
        let tz = r.pow(&half5, 53, RoundingMode::ToZero, &mut cc);
        assert!(tz.is_zero() && tz.inexact());
        // a root that is not exact takes the logarithm path, which also has to round the underflow
        let s = ExactNum::from_mag_parts(0xC000_0000_0000_0000, EXPONENT_MIN - 17, 64)?; // 1.5 * 2^(EMIN - 18)
        let up = s.pow(&half5, 53, RoundingMode::Up, &mut cc);
        assert!(up.is_positive() && !up.is_zero() && up.inexact());
        let tz = s.pow(&half5, 53, RoundingMode::ToZero, &mut cc);
        assert!(tz.is_zero() && tz.inexact());
        let sq = ExactNum::mul_ref(&x, &x, 53, RoundingMode::ToZero); // flushed to zero
        assert!(sq.is_zero() && sq.inexact());
        Ok(())
    }
}

// The raw limbs below are 64-bit words.
#[cfg(all(test, not(target_pointer_width = "32")))]
mod full_prec_subnormal_tests {
    use super::{RawKind, RawParts};
    use crate::{Error, ExactNum, Exponent, Word, EXPONENT_MIN};

    fn raw(words: &[Word], negative: bool, exponent: Exponent) -> Result<ExactNum, Error> {
        ExactNum::from_raw_parts(&RawParts {
            kind: RawKind::Finite,
            words,
            precision: words.len() * crate::WORD_BIT_SIZE,
            negative,
            exponent,
            inexact: false,
            op_inexact: false,
        })
    }

    // A full-precision sum or product below the normal range was cut to the words of its
    // normalized mantissa before the subnormal shift, which dropped the low word.
    #[test]
    fn full_prec_keeps_every_bit_of_a_subnormal_result() -> Result<(), Error> {
        // a 192-bit subnormal number whose top word is zero: 0xAAAAAAA_AAAAAAAAAAAAAAAA * 2^(EMIN - 192)
        let x = raw(&[0xAAAA_AAAA_AAAA_AAAA, 0x0AAA_AAAA, 0], true, EXPONENT_MIN)?;
        let one = raw(&[0x8000_0000_0000_0000], false, 1)?;
        let twice = raw(&[0x5555_5555_5555_5554, 0x1555_5555, 0], true, EXPONENT_MIN)?;
        assert!(x.is_negative() && !x.is_zero());
        assert!(twice.is_negative() && !twice.is_zero());
        assert!(ExactNum::mul_full_prec(&x, &one) == x);
        assert!(ExactNum::add_full_prec(&x, &x) == twice);
        Ok(())
    }
}

/// What a [`RawParts`] describes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RawKind {
    /// A non-zero finite value (normal, or subnormal at `EXPONENT_MIN`).
    Finite,
    /// Zero; the sign is in `negative`.
    Zero,
    /// Infinity; the sign is in `negative`.
    Inf,
    /// Not a number.
    NaN,
}

/// Borrowed raw representation of an [`ExactNum`], produced by [`ExactNum::as_raw_parts`] and
/// consumed by [`ExactNum::from_raw_parts`].
///
/// ## Contract
///
///  - `words` is the mantissa, least-significant word first (little-endian word order, native
///    bit order inside a word). `Word` is `u64` (`u32` on 32-bit targets), `WORD_BIT_SIZE` bits.
///  - `precision` is the stored mantissa width in bits and always equals
///    `words.len() * WORD_BIT_SIZE` (a number created at precision `p` stores
///    `ceil(p / WORD_BIT_SIZE)` words, so `p = 1` reports `WORD_BIT_SIZE`). No requested
///    precision is kept apart from the word count: a value rounded at `p = 100` is stored in
///    128 bits whose low 28 are zero, so the round trip loses nothing and later operations,
///    which take their precision as an argument, behave identically.
///  - Finite value `= (-1)^negative · M · 2^(exponent − precision)`, where `M` is the unsigned
///    integer formed by `words`; equivalently `0.1m… · 2^exponent` for a normal value.
///  - Normal: the top bit of the last word is set and
///    `EXPONENT_MIN <= exponent <= EXPONENT_MAX`, i.e. `|exponent| <= 2^62 − 1`.
///  - Subnormal: top bit clear (leading zero bits or words allowed), mantissa non-zero,
///    `exponent == EXPONENT_MIN`.
///  - Zero: `kind == Zero`, every word 0, `exponent == 0`, sign of zero in `negative`;
///    a zero may carry `inexact` and `op_inexact` (an underflow to zero is inexact).
///  - Inf / NaN: `words` empty, `precision == 0`, `exponent == 0`, `inexact == false`;
///    `op_inexact == false`; `negative` is the sign of the infinity and is `false` for NaN. A
///    NaN's error kind is not carried.
///  - `inexact` is the accumulated inexact flag (see `ExactNum::inexact_accumulated`) and
///    `op_inexact` the per-operation flag (see `ExactNum::inexact`); `op_inexact` implies
///    `inexact`.
#[derive(Debug)]
pub struct RawParts<'a> {
    /// Finite, Zero, Inf or NaN.
    pub kind: RawKind,
    /// Mantissa words, least-significant first.
    pub words: &'a [Word],
    /// Stored width in bits, `words.len() * WORD_BIT_SIZE`.
    pub precision: usize,
    /// Sign.
    pub negative: bool,
    /// Binary exponent of `0.1m…` form.
    pub exponent: Exponent,
    /// Accumulated inexact flag.
    pub inexact: bool,
    /// Per-operation inexact flag; implies `inexact`.
    pub op_inexact: bool,
}

impl RawParts<'_> {
    fn special(kind: RawKind, negative: bool) -> Self {
        RawParts {
            kind,
            words: &[],
            precision: 0,
            negative,
            exponent: 0,
            inexact: false,
            op_inexact: false,
        }
    }

    /// Number of bits from the lowest bit up to and including the highest set bit of `words`;
    /// equals `precision` for a normal value and 0 for zero.
    ///
    /// ## Errors
    ///
    ///  - InvalidArgument: the bit count does not fit in `usize`.
    pub fn significant_bits(&self) -> Result<usize, Error> {
        match self.words.iter().rposition(|w| *w != 0) {
            Some(i) => match self.words.get(i) {
                Some(w) => {
                    let lz =
                        usize::try_from(w.leading_zeros()).map_err(|_| Error::InvalidArgument)?;
                    i.checked_add(1)
                        .and_then(|k| k.checked_mul(WORD_BIT_SIZE))
                        .and_then(|bits| bits.checked_sub(lz))
                        .ok_or(Error::InvalidArgument)
                }
                None => Ok(0),
            },
            None => Ok(0),
        }
    }
}
