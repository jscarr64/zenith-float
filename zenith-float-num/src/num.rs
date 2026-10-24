// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! ExactNumNumber definition and basic arithmetic, comparison, and number manipulation operations.

use crate::common::consts::ONE;
use crate::common::util::count_leading_zeroes_skip_first;
use crate::common::util::round_p;
use crate::defs::Error;
use crate::defs::Exponent;
use crate::defs::Precision;
use crate::defs::RoundingMode;
use crate::defs::Sign;
use crate::defs::SignedWord;
use crate::defs::Word;
use crate::defs::EXPONENT_MAX;
use crate::defs::EXPONENT_MIN;
use crate::defs::MAX_PRECISION;
use crate::defs::WORD_BIT_SIZE;
use crate::defs::WORD_SIGNIFICANT_BIT;
use crate::exp;
use crate::mantissa::AddSubSpec;
use crate::mantissa::Mantissa;

// The one approved Ziv working-copy helper (see its file); used only at the pinned sites in tools/ziv_shape_a_sites.txt.
mod ziv_operand;

/// A finite floating point number with mantissa of an arbitrary size, an exponent, and the sign.
#[derive(Debug, Hash)]
pub(crate) struct ExactNumNumber {
    e: Exponent,
    s: Sign,
    m: Mantissa,
    // Accumulated: set when this value or any value it was computed from was rounded.
    // Kernel code (Ziv loops, `try_set_precision`, repeated squaring) relies on it.
    inexact: bool,
    // Per operation: `Some(r)` when the operation that produced this value records whether
    // it rounded (`add_per_op` and friends). `None` falls back to `inexact`. Any later
    // in-place change of the value resets it.
    op_inexact: Option<bool>,
}

impl ExactNumNumber {
    #[inline]
    pub(crate) fn is_subnormal(&self) -> bool {
        // Leading bit unset (and not zero): same as the public ExactNum docs.
        !self.is_zero() && self.m.bit_len() < self.m.max_bit_len()
    }

    pub(crate) fn as_raw_parts(&self) -> (&[Word], usize, Sign, Exponent, bool) {
        let (m_slice, n) = self.m.as_raw_parts();
        (m_slice, n, self.s, self.e, self.inexact)
    }

    // Check the precision so it does not cause arithmetic overflows anywhere.
    pub(super) fn p_assertion(p: usize) -> Result<(), Error> {
        // Exponent arithmetic is done in i64; sizes in usize, which is 32 bits wide on some targets.
        if p == 0 || p > crate::defs::MAX_PRECISION {
            Err(Error::InvalidArgument)
        } else {
            Ok(())
        }
    }

    /// Returns a new number with value of 0 and precision of `p` bits. Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - InvalidArgument: the precision is incorrect.
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    pub fn new(p: usize) -> Result<Self, Error> {
        Self::p_assertion(p)?;
        Ok(ExactNumNumber {
            m: Mantissa::new(p)?,
            e: 0,
            s: Sign::Pos,
            inexact: false,
            op_inexact: None,
        })
    }

    /// Like [`Self::new`], but also accepts `p = 0` (an empty mantissa). Crate-private: the full-precision
    /// kernels (`add_full_prec`, `sub_full_prec`, `mul_full_prec`) use `p = 0` to mean "as wide as the exact
    /// result", and the constants cache builds on them. Public entry points reject `p = 0`.
    pub(crate) fn new_internal(p: usize) -> Result<Self, Error> {
        if p > crate::defs::MAX_PRECISION {
            return Err(Error::InvalidArgument);
        }
        Ok(ExactNumNumber {
            m: Mantissa::new(p)?,
            e: 0,
            s: Sign::Pos,
            inexact: false,
            op_inexact: None,
        })
    }

    /// Returns a new number with value of 0, precision of `p` bits, sign `s`, and marked as inexact if `inexact` is true.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - InvalidArgument: the precision is incorrect.
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    pub fn new2(p: usize, s: Sign, inexact: bool) -> Result<Self, Error> {
        Self::p_assertion(p)?;
        Ok(ExactNumNumber {
            m: Mantissa::new(p)?,
            e: 0,
            s,
            inexact,
            op_inexact: None,
        })
    }

    /// `d` at [`DEFAULT_P`](crate::defs::DEFAULT_P) bits, built without allocating or failing.
    /// Equal to `from_word(d, DEFAULT_P)`.
    pub(crate) const fn small_const(d: Word) -> Self {
        if d == 0 {
            return ExactNumNumber {
                m: Mantissa::zero_default(),
                e: 0,
                s: Sign::Pos,
                inexact: false,
                op_inexact: None,
            };
        }
        let shift = d.leading_zeros();
        ExactNumNumber {
            m: Mantissa::from_normalized_top_word(d << shift),
            e: (WORD_BIT_SIZE - shift as usize) as Exponent,
            s: Sign::Pos,
            inexact: false,
            op_inexact: None,
        }
    }

    /// Returns the maximum value for the specified precision `p`: all bits of the mantissa are set to 1,
    /// the exponent has the maximum possible value, and the sign is positive.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - InvalidArgument: the precision is incorrect.
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    pub fn max_value(p: usize) -> Result<Self, Error> {
        Self::p_assertion(p)?;
        Ok(ExactNumNumber {
            m: Mantissa::oned_mantissa(p)?,
            e: EXPONENT_MAX,
            s: Sign::Pos,
            inexact: false,
            op_inexact: None,
        })
    }

    /// Returns the minimum value for the specified precision `p`: all bits of the mantissa are set to 1, the exponent has the maximum possible value, and the sign is negative. Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - InvalidArgument: the precision is incorrect.
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    pub fn min_value(p: usize) -> Result<Self, Error> {
        Self::p_assertion(p)?;
        Ok(ExactNumNumber {
            m: Mantissa::oned_mantissa(p)?,
            e: EXPONENT_MAX,
            s: Sign::Neg,
            inexact: false,
            op_inexact: None,
        })
    }

    /// Returns the minimum positive subnormal value for the specified precision `p`:
    /// only the least significant bit of the mantissa is set to 1, the exponent has
    /// the minimum possible value, and the sign is positive.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - InvalidArgument: the precision is incorrect.
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    pub fn min_positive(p: usize) -> Result<Self, Error> {
        Self::p_assertion(p)?;
        Ok(ExactNumNumber {
            m: Mantissa::min(p)?,
            e: EXPONENT_MIN,
            s: Sign::Pos,
            inexact: false,
            op_inexact: None,
        })
    }

    /// Returns the minimum positive normal value for the specified precision `p`:
    /// only the most significant bit of the mantissa is set to 1, the exponent has
    /// the minimum possible value, and the sign is positive.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - InvalidArgument: the precision is incorrect.
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    pub fn min_positive_normal(p: usize) -> Result<Self, Error> {
        Self::p_assertion(p)?;
        Ok(ExactNumNumber {
            m: Mantissa::from_word(p, WORD_SIGNIFICANT_BIT)?,
            e: EXPONENT_MIN,
            s: Sign::Pos,
            inexact: false,
            op_inexact: None,
        })
    }

    /// Returns a new number with value `d` and the precision `p`. Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - InvalidArgument: the precision is incorrect.
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    pub fn from_word(mut d: Word, p: usize) -> Result<Self, Error> {
        Self::p_assertion(p)?;

        if d == 0 {
            Self::new(p)
        } else {
            let mut shift = 0;
            while d & WORD_SIGNIFICANT_BIT == 0 {
                d <<= 1;
                shift += 1;
            }
            Ok(ExactNumNumber {
                m: Mantissa::from_word(p, d)?,
                e: (WORD_BIT_SIZE - shift) as Exponent,
                s: Sign::Pos,
                inexact: false,
                op_inexact: None,
            })
        }
    }

    pub(crate) fn from_u64_internal(d: u64, p: usize) -> Result<Self, Error> {
        #[cfg(not(target_pointer_width = "32"))]
        {
            Self::from_word(d, p)
        }

        #[cfg(target_pointer_width = "32")]
        {
            Self::p_assertion(p)?;

            if d == 0 {
                Self::new(p)
            } else {
                let mut d = d;
                let mut shift = 0;
                while d < 0x8000000000000000u64 {
                    d <<= 1;
                    shift += 1;
                }
                let words = [d as Word, (d >> WORD_BIT_SIZE) as Word];
                if p < WORD_BIT_SIZE * 2 {
                    if words[0] != 0 {
                        Err(Error::InvalidArgument)
                    } else {
                        Ok(ExactNumNumber {
                            m: Mantissa::from_word(p, words[1])?,
                            e: (64 - shift) as Exponent,
                            s: Sign::Pos,
                            inexact: false,
                            op_inexact: None,
                        })
                    }
                } else {
                    Ok(ExactNumNumber {
                        m: Mantissa::from_words(p, &[d as Word, (d >> WORD_BIT_SIZE) as Word])?,
                        e: (64 - shift) as Exponent,
                        s: Sign::Pos,
                        inexact: false,
                        op_inexact: None,
                    })
                }
            }
        }
    }

    /// Returns a copy of the number with the sign reversed.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    pub fn neg(&self) -> Result<Self, Error> {
        let mut ret = self.clone()?;
        ret.s = ret.s.invert();
        Ok(ret)
    }

    /// Adds `d2` to `self` and returns the result of the operation with precision `p` rounded according to `rm`.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - ExponentOverflow: the resulting exponent becomes greater than the maximum allowed value for the exponent.
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - InvalidArgument: the precision is incorrect.
    #[inline]
    pub fn add(&self, d2: &Self, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        self.add_sub(d2, p, 1, rm, false, false)
    }

    /// Subtracts `d2` from `self` and returns the result of the operation with precision `p` rounded according to `rm`.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - ExponentOverflow: the resulting exponent becomes greater than the maximum allowed value for the exponent.
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - InvalidArgument: the precision is incorrect.
    #[inline]
    pub fn sub(&self, d2: &Self, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        self.add_sub(d2, p, -1, rm, false, false)
    }

    /// Adds `d2` to `self` and returns the result of the operation.
    /// The resulting precision is equal to the full precision of the result.
    /// This operation can be used to emulate integer addition.
    ///
    /// ## Errors
    ///
    ///  - ExponentOverflow: the resulting exponent becomes greater than the maximum allowed value for the exponent.
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    #[inline]
    pub fn add_full_prec(&self, d2: &Self) -> Result<Self, Error> {
        self.add_sub(d2, 0, 1, RoundingMode::None, true, false)
    }

    /// Subtracts `d2` from `self` and returns the result of the operation.
    /// The resulting precision is equal to the full precision of the result.
    /// This operation can be used to emulate integer subtraction.
    ///
    /// ## Errors
    ///
    ///  - ExponentOverflow: the resulting exponent becomes greater than the maximum allowed value for the exponent.
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    #[inline]
    pub fn sub_full_prec(&self, d2: &Self) -> Result<Self, Error> {
        self.add_sub(d2, 0, -1, RoundingMode::None, true, false)
    }

    /// Multiplies `d2` by `self` and returns the result of the operation with precision `p` rounded according to `rm`.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - ExponentOverflow: the resulting exponent becomes greater than the maximum allowed value for the exponent.
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - InvalidArgument: the precision is incorrect.
    #[inline]
    pub fn mul(&self, d2: &Self, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        self.mul_general_case(d2, p, rm, false, false)
    }

    /// Multiplies `d2` by `self` and returns the result of the operation.
    /// The resulting precision is equal to the full precision of the result.
    /// This operation can be used to emulate integer multiplication.
    ///
    /// ## Errors
    ///
    ///  - ExponentOverflow: the resulting exponent becomes greater than the maximum allowed value for the exponent.
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    #[inline]
    pub fn mul_full_prec(&self, d2: &Self) -> Result<Self, Error> {
        self.mul_general_case(d2, 0, RoundingMode::None, true, false)
    }

    /// Adds like `add`. The result also records whether this addition rounded (`op_inexact`).
    #[inline]
    pub(crate) fn add_per_op(&self, d2: &Self, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        self.add_sub(d2, p, 1, rm, false, true)
    }

    /// Subtracts like `sub`. The result also records whether this subtraction rounded.
    #[inline]
    pub(crate) fn sub_per_op(&self, d2: &Self, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        self.add_sub(d2, p, -1, rm, false, true)
    }

    /// Multiplies like `mul`. The result also records whether this multiplication rounded.
    #[inline]
    pub(crate) fn mul_per_op(&self, d2: &Self, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        self.mul_general_case(d2, p, rm, false, true)
    }

    /// Divides like `div`. The result also records whether this division rounded.
    #[inline]
    pub(crate) fn div_per_op(&self, d2: &Self, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        self.div_general_case(d2, p, rm, true)
    }

    fn mul_general_case(
        &self,
        d2: &Self,
        p: usize,
        rm: RoundingMode,
        full_prec: bool,
        per_op: bool,
    ) -> Result<Self, Error> {
        let p = round_p(p);

        if !full_prec {
            Self::p_assertion(p)?;
        }

        let s = if self.s == d2.s { Sign::Pos } else { Sign::Neg };

        if self.m.is_zero() || d2.m.is_zero() {
            let mut ret = Self::new_internal(p)?;
            ret.set_sign(s);
            return Ok(ret);
        }

        // The product kernel takes the operands as stored, a subnormal one included (it is exact and is
        // renormalized by the kernel), so no operand is copied or shifted.
        let inputs = self.inexact || d2.inexact;
        let mut inexact = false;

        let (e_shift, m3) = self
            .m
            .mul(&d2.m, p, rm, s == Sign::Pos, full_prec, &mut inexact)?;

        let e = exp::narrow(exp::wide(self.e) + exp::wide(d2.e) - exp::wide(e_shift));

        if e > EXPONENT_MAX as i64 {
            return Err(Error::ExponentOverflow(s));
        }

        if e < EXPONENT_MIN as i64 {
            let mut ret = ExactNumNumber {
                m: m3,
                s,
                e: EXPONENT_MIN,
                inexact,
                op_inexact: None,
            };

            if full_prec {
                ret.widen_for_subnormal(e)?;
            }
            ret.subnormalize(e, rm);
            ret.finish_flags(inputs, per_op);

            Ok(ret)
        } else {
            let mut ret = ExactNumNumber {
                m: m3,
                s,
                e: e as Exponent,
                inexact,
                op_inexact: None,
            };
            ret.finish_flags(inputs, per_op);

            Ok(ret)
        }
    }

    /// Divides `self` by `d2` and returns the result of the operation with precision `p` rounded according to `rm`.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - DivisionByZero: `d2` is zero.
    ///  - ExponentOverflow: the resulting exponent becomes greater than the maximum allowed value for the exponent.
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - InvalidArgument: both `self` and `d2` are zero or precision is incorrect.
    pub fn div(&self, d2: &Self, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        self.div_general_case(d2, p, rm, false)
    }

    fn div_general_case(
        &self,
        d2: &Self,
        p: usize,
        rm: RoundingMode,
        per_op: bool,
    ) -> Result<Self, Error> {
        if d2.m.is_zero() {
            return if self.is_zero() {
                Err(Error::InvalidArgument)
            } else {
                Err(Error::DivisionByZero)
            };
        }

        let p = round_p(p);
        Self::p_assertion(p)?;

        let s = if self.s == d2.s { Sign::Pos } else { Sign::Neg };

        if self.m.is_zero() {
            let mut ret = Self::new(p)?; // self / d2 = 0
            ret.set_sign(s);
            return Ok(ret);
        }

        // The division kernel takes the operands as stored (it accounts for a subnormal operand's leading
        // zero bits itself), so no operand is cloned.
        let inputs = self.inexact || d2.inexact;
        let mut inexact = false;

        let (e_shift, m3) = self.m.div(&d2.m, p, rm, s == Sign::Pos, &mut inexact)?;

        let e = exp::narrow(exp::wide(self.e) - exp::wide(d2.e) + exp::wide(e_shift));

        if e > EXPONENT_MAX as i64 {
            return Err(Error::ExponentOverflow(s));
        }

        if e < EXPONENT_MIN as i64 {
            let mut ret = ExactNumNumber {
                m: m3,
                s,
                e: EXPONENT_MIN,
                inexact,
                op_inexact: None,
            };

            ret.subnormalize(e, rm);
            ret.finish_flags(inputs, per_op);

            Ok(ret)
        } else {
            let mut ret = ExactNumNumber {
                m: m3,
                s,
                e: e as Exponent,
                inexact,
                op_inexact: None,
            };
            ret.finish_flags(inputs, per_op);

            Ok(ret)
        }
    }

    /// Returns the remainder of division of `|self|` by `|d2|`. The sign of the result is set to the sign of `self`.
    ///
    /// ## Errors
    ///
    ///  - ExponentOverflow: the resulting exponent becomes greater than the maximum allowed value for the exponent.
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - InvalidArgument: `d2` is zero.
    pub fn rem(&self, d2: &Self) -> Result<Self, Error> {
        if d2.m.is_zero() {
            return Err(Error::InvalidArgument);
        }

        if self.m.is_zero() {
            return self.clone();
        }

        // represent as mN * 2 ^ eNeff
        let e1eff = self.e as i64 - self.m.max_bit_len() as i64;

        let (e, m2_opt) = d2.normalize()?;
        let m2_normalized = m2_opt.as_ref().unwrap_or(&d2.m);
        let e2eff = e - m2_normalized.max_bit_len() as i64;

        let finalize = |mut m3: Mantissa, mut e: i64| -> Result<ExactNumNumber, Error> {
            let (m3, e) = if m3.bit_len() > 0 {
                if e > EXPONENT_MAX as i64 {
                    return Err(Error::ExponentOverflow(self.s));
                }

                if e < EXPONENT_MIN as i64 {
                    // first, see if m3 is already subnormal, then shift accordingly
                    let mut excess = m3.max_bit_len() - m3.bit_len();

                    // `usize::MAX as i64` is -1, so the old range check flushed every
                    // subnormal remainder to zero; convert the shift instead.
                    if let Ok(shift) = usize::try_from(EXPONENT_MIN as i64 - e) {
                        if excess < shift {
                            m3.extend_subnormal(shift)?;
                            excess = m3.max_bit_len() - m3.bit_len();
                        }

                        let shift = excess - shift;
                        m3.shift_left(shift);
                        m3.set_bit_len(m3.bit_len() + shift);

                        e = EXPONENT_MIN as i64;
                    } else {
                        m3.set_zero();
                        e = 0;
                    }
                } else {
                    m3.normilize2();
                }

                (m3, e)
            } else {
                (m3, 0)
            };

            let ret = ExactNumNumber {
                m: m3,
                s: self.s,
                e: e as Exponent,
                inexact: self.inexact || d2.inexact,
                op_inexact: None,
            };

            Ok(ret)
        };

        if e1eff > e2eff {
            // (m1 * 2^(e1eff - e2eff)) mod m2 = ((m1 mod m2) * ( 2^(e1eff - e2eff) mod m2 )) mod m2

            let mut two = Mantissa::from_raw_parts(&[2], 2)?;
            let powm = if (e1eff - e2eff) as usize >= m2_normalized.max_bit_len() {
                two.pow_mod((e1eff - e2eff) as usize, m2_normalized)?
            } else {
                two.pow2((e1eff - e2eff - 1) as usize)?;
                two // m3 still becomes 0 in rare case when two == m2_normalized
            };

            let m3 = self.m.rem(m2_normalized)?;

            let m3 = m3.mul_mod(&powm, m2_normalized)?;

            let e = e2eff + m3.bit_len() as i64;

            finalize(m3, e)
        } else if (self.m.bit_len() as i64 + self.e as i64) < (m2_normalized.bit_len() as i64 + e) {
            // self < d2, remainder = self
            self.clone()
        } else {
            // since self.e >= d2.e and e1eff <= e2eff, then e2eff - e1eff < m1.len()
            // (m1 * 2 ^ e1eff) mod (m2 * 2 ^ e2eff) = m1 mod (m2 * 2 ^ (e2eff - e1eff))

            let mut m2_normalized = if let Some(m2) = m2_opt { m2 } else { d2.m.clone()? };

            let e2eff = e - m2_normalized.max_bit_len() as i64;
            let ediff = (e2eff - e1eff) as usize;

            let m2l = m2_normalized.max_bit_len() + ediff;
            m2_normalized.set_length(m2l)?;

            let m3 = if m2_normalized.max_bit_len() > m2l {
                let mut m = self.m.clone()?;
                m.pow2(m2_normalized.max_bit_len() - m2l)?;
                m.rem(&m2_normalized)
            } else {
                self.m.rem(&m2_normalized)
            }?;

            let e = e1eff + m3.bit_len() as i64 - m2_normalized.max_bit_len() as i64 + m2l as i64;

            finalize(m3, e)
        }
    }

    // Return normilized mantissa and exponent with corresponding shift.
    pub(super) fn normalize(&self) -> Result<(i64, Option<Mantissa>), Error> {
        if self.is_subnormal() {
            let mut mantissa = self.m.clone()?;
            let shift = mantissa.normilize2();

            #[cfg(not(target_pointer_width = "32"))]
            {
                // A shift is at most the mantissa length (< 2^60), so the result stays far inside i64.
                debug_assert!(shift <= MAX_PRECISION);
            }

            Ok((
                exp::narrow(exp::wide(self.e) - exp::wide_u(shift)),
                Some(mantissa),
            ))
        } else {
            Ok((self.e as i64, None))
        }
    }

    // Normalize mantissa and return exponent shift of `self`.
    pub(crate) fn normalize2(&mut self) -> usize {
        self.m.normilize2()
    }

    // Combined add and sub operations.
    fn add_sub(
        &self,
        d2: &Self,
        p: usize,
        op: i8,
        rm: RoundingMode,
        full_prec: bool,
        per_op: bool,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        if !full_prec {
            Self::p_assertion(p)?;
        }

        let mut d3 = Self::new_internal(0)?;

        // Both operands are zero: the result is a zero of precision `p`, signed by the rules for zeros.
        if self.m.is_zero() && d2.m.is_zero() {
            let other = if op < 0 { d2.s.invert() } else { d2.s };
            let negative =
                if self.s == other { self.s == Sign::Neg } else { rm == RoundingMode::Down };
            let mut ret = Self::new_internal(p)?;
            ret.finish_flags(d2.inexact, per_op);
            ret.set_sign(if negative { Sign::Neg } else { Sign::Pos });
            return Ok(ret);
        }

        // One operand is zero: the general kernel sums the other operand with a zero mantissa, so the
        // result is written (and rounded to `p`) by the same code as every other sum.
        if self.m.is_zero() || d2.m.is_zero() {
            let (src, sign) = if self.m.is_zero() {
                (d2, if op < 0 { d2.s.invert() } else { d2.s })
            } else {
                (self, self.s)
            };
            let (e0, m0) = (src.e as i64, &src.m);
            let inputs = src.inexact;
            let mut inexact = false;
            d3.s = sign;
            let zero = Mantissa::new(WORD_BIT_SIZE)?;
            let (shift, m3) = m0.abs_add(
                &zero,
                AddSubSpec {
                    p,
                    m2_shift: 1,
                    rm,
                    is_positive: d3.is_positive(),
                    full_prec,
                },
                &mut inexact,
            )?;
            let mut ret =
                Self::finish_sum(d3, e0 - shift, m3, (inexact, inputs), rm, per_op, full_prec)?;
            if ret.is_zero() {
                // A subnormal operand rounded to zero keeps the minimal exponent, as `set_precision` left it.
                ret.e = EXPONENT_MIN;
            }
            return Ok(ret);
        }

        // The sum kernels take the operands as stored (a subnormal operand has leading zero bits and is
        // renormalized by the kernel together with the result), so no operand is cloned.
        let (e1, m1) = (self.e as i64, &self.m);
        let (e2, m2) = (d2.e as i64, &d2.m);

        let inputs = self.inexact || d2.inexact;
        let mut inexact = false;

        let e;
        let (shift, m3) = if (self.s != d2.s && op >= 0) || (op < 0 && self.s == d2.s) {
            // subtract
            let cmp = self.abs_cmp(d2);
            if cmp > 0 {
                d3.s = self.s;
                e = e1;
                m1.abs_sub(
                    m2,
                    AddSubSpec {
                        p,
                        m2_shift: exp::abs_diff_shift(e1, e2),
                        rm,
                        is_positive: d3.is_positive(),
                        full_prec,
                    },
                    &mut inexact,
                )
            } else if cmp < 0 {
                d3.s = if op >= 0 { d2.s } else { d2.s.invert() };
                e = e2;
                m2.abs_sub(
                    m1,
                    AddSubSpec {
                        p,
                        m2_shift: exp::abs_diff_shift(e2, e1),
                        rm,
                        is_positive: d3.is_positive(),
                        full_prec,
                    },
                    &mut inexact,
                )
            } else {
                let mut ret = Self::new_internal(p)?;
                ret.inexact = inexact;
                ret.finish_flags(inputs, per_op);
                if rm == RoundingMode::Down {
                    ret.set_sign(Sign::Neg);
                }
                return Ok(ret);
            }
        } else {
            // add
            d3.s = self.s;
            if e1 >= e2 {
                e = e1;
                m1.abs_add(
                    m2,
                    AddSubSpec {
                        p,
                        m2_shift: exp::abs_diff_shift(e1, e2),
                        rm,
                        is_positive: d3.is_positive(),
                        full_prec,
                    },
                    &mut inexact,
                )
            } else {
                e = e2;
                m2.abs_add(
                    m1,
                    AddSubSpec {
                        p,
                        m2_shift: exp::abs_diff_shift(e2, e1),
                        rm,
                        is_positive: d3.is_positive(),
                        full_prec,
                    },
                    &mut inexact,
                )
            }
        }?;

        #[cfg(not(target_pointer_width = "32"))]
        {
            debug_assert!(shift <= i64::MAX / 2 && e >= -exp::WIDE_LIMIT);
        }
        Self::finish_sum(
            d3,
            exp::narrow(exp::wide(e) - exp::wide(shift)),
            m3,
            (inexact, inputs),
            rm,
            per_op,
            full_prec,
        )
    }

    /// Writes the exponent and mantissa of a finished sum into `d3` (sign already set) and sets the flags.
    /// `e` is the exponent after the normalization shift; `flags` is (this sum rounded, an input was inexact).
    fn finish_sum(
        mut d3: Self,
        e: i64,
        m3: Mantissa,
        flags: (bool, bool),
        rm: RoundingMode,
        per_op: bool,
        full_prec: bool,
    ) -> Result<Self, Error> {
        let (inexact, inputs) = flags;
        d3.inexact |= inexact;

        if e > EXPONENT_MAX as i64 {
            return Err(Error::ExponentOverflow(d3.s));
        }

        d3.m = m3;

        if e < EXPONENT_MIN as i64 {
            d3.e = EXPONENT_MIN;
            if full_prec {
                d3.widen_for_subnormal(e)?;
            }
            d3.subnormalize(e, rm);
        } else {
            d3.e = e as Exponent;
        }
        d3.finish_flags(inputs, per_op);

        Ok(d3)
    }

    /// A full-precision result below the normal range keeps every bit: the mantissa grows by the
    /// subnormal shift `EXPONENT_MIN - e` before `subnormalize`, which then drops only zeros. A shift
    /// past `MAX_PRECISION` cannot be kept and is rounded by `subnormalize` as before.
    fn widen_for_subnormal(&mut self, e: i64) -> Result<(), Error> {
        let sh = i128::from(EXPONENT_MIN) - i128::from(e);
        if let Ok(sh) = usize::try_from(sh) {
            let p2 = self.mantissa_max_bit_len().saturating_add(sh);
            if sh > 0 && p2 <= MAX_PRECISION {
                self.set_precision(p2, RoundingMode::None)?;
            }
        }
        Ok(())
    }

    /// Make `self` subnormal
    pub(crate) fn subnormalize(&mut self, e: i64, rm: RoundingMode) {
        debug_assert_eq!(self.exponent(), EXPONENT_MIN);

        if self.is_zero() {
            return;
        }

        let is_positive = self.is_positive();

        if (self.m.max_bit_len() as i64) + e > EXPONENT_MIN as i64 {
            // subnormal

            let mut shift = (EXPONENT_MIN as i64 - e) as usize;

            let rounded_before = self.op_inexact();
            let mut inexact = false;

            if self.m.round_mantissa(
                shift,
                rm,
                is_positive,
                &mut false,
                self.m.max_bit_len(),
                &mut inexact,
            ) {
                shift -= 1;
            }

            self.inexact |= inexact;
            self.op_inexact = Some(rounded_before || inexact);

            if shift > 0 {
                self.m.shift_right(shift);
                self.m.update_bit_len();
            }
        } else {
            // Every significand bit sits strictly below the minimum positive.
            // Compare the magnitude with half an ulp of that minimum.
            let pbits = self.m.max_bit_len() as i64;
            let below_half = e < EXPONENT_MIN as i64 - pbits;
            let class = if below_half { -1 } else { Self::mantissa_vs_half(&self.m) };
            self.inexact = true;
            self.op_inexact = Some(true);
            if Self::flush_to_min(class, rm, is_positive) {
                self.m.set_zero();
                self.m.digits_mut()[0] = 1;
                self.m.set_bit_len(1);
            } else {
                self.m.set_zero();
                self.e = 0;
            }
        }
    }

    /// `-1` below half, `0` exact tie, `1` above half of the minimum positive.
    fn mantissa_vs_half(m: &Mantissa) -> i8 {
        let words = m.digits();
        if words.is_empty() {
            return -1;
        }
        let last = words.len() - 1;
        let top = words[last];
        if top > WORD_SIGNIFICANT_BIT {
            return 1;
        }
        if top < WORD_SIGNIFICANT_BIT {
            return -1;
        }
        for word in &words[..last] {
            if *word != 0 {
                return 1;
            }
        }
        0
    }

    /// True when a value in `(0, min_positive)` rounds to `min_positive`.
    fn flush_to_min(class: i8, rm: RoundingMode, is_positive: bool) -> bool {
        let away = rm == RoundingMode::FromZero
            || (is_positive && rm == RoundingMode::Up)
            || (!is_positive && rm == RoundingMode::Down);
        match class {
            1 => {
                !matches!(rm, RoundingMode::ToZero | RoundingMode::None)
                    && (rm != RoundingMode::Down || !is_positive)
                    && (rm != RoundingMode::Up || is_positive)
            }
            0 => match rm {
                RoundingMode::ToOdd | RoundingMode::FromZero => true,
                RoundingMode::Up => is_positive,
                RoundingMode::Down => !is_positive,
                _ => false,
            },
            _ => away,
        }
    }

    /// Compares `self` to `d2`.
    /// Returns positive if `self` is greater than `d2`, negative if `self` is smaller than `d2`, 0 otherwise.
    pub fn cmp(&self, d2: &Self) -> SignedWord {
        if self.is_zero() && d2.is_zero() {
            return 0;
        }

        if self.s != d2.s {
            return self.s as SignedWord;
        }

        self.abs_cmp(d2) * self.s as SignedWord
    }

    /// Compares the absolute value of `self` to the absolute value of `d2`.
    /// Returns positive if `|self|` is greater than `|d2|`, negative if `|self|` is smaller than `|d2|`, 0 otherwise.
    pub fn abs_cmp(&self, d2: &Self) -> SignedWord {
        if self.m.is_zero() {
            if d2.m.is_zero() {
                return 0;
            } else {
                return -1;
            }
        }

        if d2.m.is_zero() {
            return 1;
        }

        let n1 = self.mantissa_max_bit_len() as i64 - self.precision() as i64;

        let n2 = d2.mantissa_max_bit_len() as i64 - d2.precision() as i64;

        let e = exp::wide(self.e) - exp::wide(n1) - exp::wide(d2.e) + exp::wide(n2);
        if e > 0 {
            return 1;
        } else if e < 0 {
            return -1;
        }

        self.m.abs_cmp(&d2.m, n1 == n2)
    }

    /// Returns the absolute value of a number.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    pub fn abs(&self) -> Result<Self, Error> {
        let mut ret = self.clone()?;
        ret.s = Sign::Pos;
        Ok(ret)
    }

    /// Constructs a number from the raw parts:
    ///
    ///  - `m` is the mantissa.
    ///  - `n` is the number of significant bits in mantissa.
    ///  - `s` is the sign.
    ///  - `e` is the exponent.
    ///  - `inexact` specify whether number is inexact.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - InvalidArgument: `n` is larger than the number of bits in `m`;
    ///    `n` is smaller than the number of bits in `m`, but `m` does not represent corresponding subnormal number mantissa;
    ///    `n` is smaller than the number of bits in `m`, but `e` is not the minimum possible exponent;
    ///    `n` or the size of `m` is too large (larger than i64::MAX / 2 + EXPONENT_MIN);
    ///    `e` is less than EXPONENT_MIN or greater than EXPONENT_MAX.
    pub fn from_raw_parts(
        m: &[Word],
        n: usize,
        s: Sign,
        mut e: Exponent,
        inexact: bool,
    ) -> Result<Self, Error> {
        let p = m.len() * WORD_BIT_SIZE;
        Self::p_assertion(p)?;

        if n > p {
            Err(Error::InvalidArgument)
        } else {
            let m = Mantissa::from_words(p, m)?;

            if m.bit_len() != n {
                return Err(Error::InvalidArgument);
            }

            if n > 0 && n < p && e > EXPONENT_MIN {
                return Err(Error::InvalidArgument);
            }

            if n == 0 {
                e = 0;
            }

            Ok(ExactNumNumber {
                e,
                s,
                m,
                inexact,
                op_inexact: None,
            })
        }
    }

    /// Build ExactNumNumber from raw parts unchecked.
    pub(super) fn from_raw_unchecked(m: Mantissa, s: Sign, e: Exponent, inexact: bool) -> Self {
        ExactNumNumber {
            e,
            s,
            m,
            inexact,
            op_inexact: None,
        }
    }

    /// Constructs a number from the slice of words:
    ///
    ///  - `m` is the mantissa.
    ///  - `s` is the sign.
    ///  - `e` is the exponent.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - InvalidArgument: `m` is empty, or its size is larger than i64::MAX / 2 + EXPONENT_MIN, or
    ///    when `e` is less than EXPONENT_MIN or greater than EXPONENT_MAX.
    pub fn from_words(m: &[Word], s: Sign, mut e: Exponent) -> Result<Self, Error> {
        let p = m.len() * WORD_BIT_SIZE;
        Self::p_assertion(p)?;

        let mut m = Mantissa::from_words(p, m)?;
        let bl = m.bit_len();

        if m.is_zero() {
            e = 0;
        } else if bl < p {
            let mut shift = p - bl;
            if (e as i64) < EXPONENT_MIN as i64 + shift as i64 {
                shift = (e as i64 - EXPONENT_MIN as i64) as usize;
                e = EXPONENT_MIN;
            } else {
                e -= shift as Exponent;
            }
            m.shift_left(shift);
            m.set_bit_len(bl + shift);
        }

        Ok(ExactNumNumber {
            e,
            s,
            m,
            inexact: false,
            op_inexact: None,
        })
    }

    /// Returns the sign of a number.
    #[inline]
    pub fn sign(&self) -> Sign {
        self.s
    }

    /// Returns true if `self` is positive.
    #[inline]
    pub fn is_positive(&self) -> bool {
        self.s == Sign::Pos
    }

    /// Returns true if `self` is negative.
    #[inline]
    pub fn is_negative(&self) -> bool {
        self.s == Sign::Neg
    }

    /// Returns the exponent of `self`.
    #[inline]
    pub fn exponent(&self) -> Exponent {
        self.e
    }

    /// Returns true if `self` is zero.
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.m.is_zero()
    }

    /// Returns the largest integer less than or equal to `self`.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    pub fn floor(self) -> Result<Self, Error> {
        let below = self.is_negative() && !self.is_int();
        let int = self.int()?;
        if below {
            return int.sub(&ONE, int.mantissa_max_bit_len(), RoundingMode::ToZero);
        }
        Ok(int)
    }

    /// Returns the smallest integer greater than or equal to `self`.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    pub fn ceil(self) -> Result<Self, Error> {
        let above = self.is_positive() && !self.is_int();
        let int = self.int()?;
        if above {
            return int.add(&ONE, int.mantissa_max_bit_len(), RoundingMode::ToZero);
        }
        Ok(int)
    }

    /// Returns fractional part of a number.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    pub fn fract(mut self) -> Result<Self, Error> {
        if self.e > 0 {
            if (self.e as usize) < self.m.max_bit_len() {
                // remove integer part of mantissa & normalize at the same time
                if let Some(shift) = self.m.find_one_from(self.e as usize) {
                    let mut e = self.e as i64 - shift as i64;
                    if e < EXPONENT_MIN as i64 {
                        let shift1 = (self.e as i64 - EXPONENT_MIN as i64) as usize;
                        self.m.shift_left(shift1);
                        self.m.mask_bits(shift - shift1, true);
                        e = EXPONENT_MIN as i64;
                    } else {
                        self.m.shift_left(shift);
                    }
                    self.e = e as Exponent;
                } else {
                    self.m.set_zero();
                    self.e = 0;
                }
            } else {
                self.m.set_zero();
                self.e = 0;
            }
        }
        Ok(self)
    }

    /// Returns integer part of a number.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    pub fn int(mut self) -> Result<Self, Error> {
        if self.e > 0 {
            if (self.e as usize) < self.m.max_bit_len() {
                self.m
                    .mask_bits(self.m.max_bit_len() - self.e as usize, false)
            }
        } else {
            self.m.set_zero();
            self.e = 0;
        }

        Ok(self)
    }

    /// Returns true if `self` is odd integer number.
    pub(crate) fn is_odd_int(&self) -> bool {
        if self.e > 0 {
            if (self.e as usize) <= self.m.max_bit_len() {
                self.m.is_odd_int(self.e as usize)
            } else {
                false
            }
        } else {
            self.is_zero()
        }
    }

    /// Returns true if `self` is a nonzero power of two (a single set bit in the mantissa).
    pub(crate) fn is_power_of_two(&self) -> bool {
        !self.is_zero() && self.m.find_one_from(1).is_none()
    }

    /// Returns true if `self` is an integer number.
    pub fn is_int(&self) -> bool {
        if self.e > 0 {
            if (self.e as usize) < self.m.max_bit_len() {
                self.m.find_one_from(self.e as usize).is_none()
            } else {
                true
            }
        } else {
            self.is_zero()
        }
    }

    /// Returns integer part of a number as built-in integer.
    pub(super) fn int_as_usize(&self) -> Result<usize, Error> {
        if self.e > 0 {
            debug_assert!(core::mem::size_of::<usize>() >= core::mem::size_of::<Word>());
            if (self.e as usize) <= WORD_BIT_SIZE {
                let shift = WORD_BIT_SIZE - self.e as usize;
                let ret = self.m.most_significant_word() as usize;
                Ok(ret >> shift)
            } else {
                Err(Error::InvalidArgument)
            }
        } else {
            Ok(0)
        }
    }

    /// Sets the exponent of `self`.
    /// Note that if `self` is subnormal, the exponent may not change, but the mantissa will shift instead.
    /// `e` will be clamped to the range from EXPONENT_MIN to EXPONENT_MAX if it's outside of the range.
    pub fn set_exponent(&mut self, e: Exponent) {
        if !self.is_zero() {
            if self.is_subnormal() && e > EXPONENT_MIN {
                let ediff = e as i64 - EXPONENT_MIN as i64;

                let n = self.mantissa_max_bit_len() - self.precision();
                if n as i64 >= ediff {
                    self.m.shift_left(ediff as usize);
                    self.m.set_bit_len(self.m.bit_len() + ediff as usize);
                } else {
                    self.m.shift_left(n);
                    self.m.set_bit_len(self.m.max_bit_len());
                    self.e = e - n as Exponent;
                }
            } else {
                self.e = e;
            }
        }
    }

    /// Returns the maximum mantissa length of `self` in bits regardless of whether `self` is normal or subnormal.
    #[inline]
    pub fn mantissa_max_bit_len(&self) -> Precision {
        self.m.max_bit_len()
    }

    /// Returns the number of significant bits used in the mantissa. Normal numbers use all bits of the mantissa.
    /// Subnormal numbers use fewer bits than the mantissa can hold.
    #[inline]
    pub fn precision(&self) -> Precision {
        self.m.bit_len()
    }

    /// Returns the rounded number with `n` binary positions in the fractional part of the number using rounding mode `rm`.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - ExponentOverflow: rounding causes exponent overflow.
    ///  - InvalidArgument: `n` is too large.
    pub fn round(mut self, n: usize, rm: RoundingMode) -> Result<Self, Error> {
        // `n` reaches `|EXPONENT_MIN| + p` when rounding onto the subnormal grid; `EXPONENT_MAX + 2 * MAX_PRECISION`
        // fits `i64`, so `-(n - 1)` below cannot wrap. The old `i64::MAX / 2` cap rejected those `n` at the i64 exponent.
        if n as u64 > EXPONENT_MAX as u64 + 2 * MAX_PRECISION as u64 {
            return Err(Error::InvalidArgument);
        }
        let mut rounded = false;

        let max_bits = self.mantissa_max_bit_len();
        let e = max_bits as i64 - self.e as i64;

        if e > 0 && e as usize > n && !self.is_zero() {
            let m = e as usize - n;

            if m >= max_bits {
                let is_positive = self.is_positive();
                let msb_set = self.m.most_significant_word() & WORD_SIGNIFICANT_BIT != 0;
                let below_half =
                    count_leading_zeroes_skip_first(self.mantissa().digits()) != max_bits;

                self.m.set_zero();

                if rm == RoundingMode::FromZero
                    || (is_positive && rm == RoundingMode::Up)
                    || (!is_positive && rm == RoundingMode::Down)
                    || (m == max_bits
                        && ((msb_set && rm == RoundingMode::ToOdd)
                            || (rm == RoundingMode::ToEven && msb_set && below_half)))
                {
                    // non zero for directed rounding modes,
                    // non zero for rounding to even/odd when msb of self is the rounding bit
                    let Some(top) = self.m.digits_mut().last_mut() else {
                        return Err(Error::InvalidArgument); // a rounded result has at least one word
                    };
                    *top = WORD_SIGNIFICANT_BIT;
                    self.m.set_bit_len(self.mantissa_max_bit_len());

                    let e = -(n as i64 - 1);
                    if e < EXPONENT_MIN as i64 {
                        self.e = EXPONENT_MIN;
                        self.subnormalize(e, rm);
                    } else {
                        self.e = e as Exponent;
                    }
                } else {
                    self.e = 0;
                }

                self.inexact |= true; // self was not zero
                rounded = true;
            } else {
                let mut inexact = false;

                let ovf = self.m.round_mantissa(
                    m,
                    rm,
                    self.is_positive(),
                    &mut false,
                    self.m.max_bit_len(),
                    &mut inexact,
                );

                self.inexact |= inexact;
                rounded = inexact;

                if ovf {
                    if self.e == EXPONENT_MAX {
                        return Err(Error::ExponentOverflow(self.s));
                    }

                    self.e += 1;
                } else if self.m.is_all_zero() {
                    self.m.set_bit_len(0);
                    self.e = 0;
                } else if self.is_subnormal() {
                    self.m.update_bit_len();
                }
            }
        }

        self.set_op_inexact(rounded);
        Ok(self)
    }

    #[cfg(any(test, feature = "random"))]
    /// Returns a random normalized (not subnormal) ExactNum number with exponent in the range
    /// from `exp_from` to `exp_to` inclusive. The sign can be positive and negative. Zero is excluded.
    /// Precision is rounded upwards to the word size.
    /// Function does not follow any specific distribution law.
    /// The intended use of this function is for testing.
    ///
    /// ## Errors
    ///
    ///  - InvalidArgument: the precision is incorrect or when `exp_from` is less than EXPONENT_MIN or `exp_to` is greater than EXPONENT_MAX.
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    pub fn random_normal(p: usize, exp_from: Exponent, exp_to: Exponent) -> Result<Self, Error> {
        Self::p_assertion(p)?;

        let m = Mantissa::random_normal(p)?;
        let e = if exp_from < exp_to {
            (crate::common::test_rng::random::<i64>().abs() % (exp_to as i64 - exp_from as i64)
                + exp_from as i64) as Exponent
        } else {
            exp_from
        };
        let s = if crate::common::test_rng::random::<u8>() & 1 == 0 {
            Sign::Pos
        } else {
            Sign::Neg
        };
        Ok(ExactNumNumber {
            e,
            s,
            m,
            inexact: false,
            op_inexact: None,
        })
    }

    /// Clones the number.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    pub fn clone(&self) -> Result<Self, Error> {
        Ok(ExactNumNumber {
            e: self.e,
            s: self.s,
            m: self.m.clone()?,
            inexact: self.inexact,
            op_inexact: self.op_inexact,
        })
    }

    /// Sets the precision of `self` to `p`.
    /// If the new precision is smaller than the existing one, the number is rounded using specified rounding mode `rm`.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - InvalidArgument: the precision is incorrect.
    pub fn set_precision(&mut self, p: usize, rm: RoundingMode) -> Result<(), Error> {
        self.set_precision_internal(p, rm, false, self.mantissa_max_bit_len())
            .map(|_| {})
    }

    /// Try to round and then set the precision to `p`, given `self` has `s` correct digits in mantissa.
    /// Returns true if rounding succeeded. If the fuction returns `false`, `self` is still modified, and should be discarded.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - InvalidArgument: the precision is incorrect.
    pub(crate) fn try_set_precision(
        &mut self,
        p: usize,
        rm: RoundingMode,
        s: usize,
    ) -> Result<bool, Error> {
        self.set_precision_internal(p, rm, self.inexact(), s)
    }

    /// Evaluates a correctly-rounded operation `op(precision, rounding)` so that the result is
    /// rounded at exactly `p` bits. Word-aligned precisions and `RoundingMode::None` call `op`
    /// directly; otherwise `op` runs toward zero with a guard word and the result is rounded to
    /// odd and then to `p` bits with `rm`.
    pub(crate) fn exact_bits<F>(p: usize, rm: RoundingMode, op: F) -> Result<Self, Error>
    where
        F: Fn(usize, RoundingMode) -> Result<Self, Error>,
    {
        Self::p_assertion(p)?;
        if rm == RoundingMode::None {
            return op(p, rm);
        }
        if p.is_multiple_of(WORD_BIT_SIZE) {
            let r = op(p, rm)?;
            // A rounded subnormal result was rounded twice (at `p` bits, then onto the
            // subnormal grid); recompute it through the round-to-odd path.
            if r.e != EXPONENT_MIN || !r.op_inexact() {
                return Ok(r);
            }
        }
        let wide = op(round_p(p) + WORD_BIT_SIZE, RoundingMode::ToZero)?;
        if wide.is_zero() && !wide.op_inexact() {
            // Exact zero: its sign depends on the rounding mode.
            return op(p, rm);
        }
        wide.round_odd_to(p, rm)
    }

    /// Rounds a result that was computed toward zero at a precision of at least `p + 2` bits
    /// to exactly `p` bits with `rm`: the sticky bit is folded into the least significant bit
    /// (round to odd) before the final rounding, so no double-rounding error can occur.
    /// The per-operation flag must say whether the toward-zero computation rounded.
    pub(crate) fn round_odd_to(mut self, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        let rounded = self.op_inexact();
        let acc = self.inexact;

        if rounded {
            if self.m.is_zero() {
                // A nonzero magnitude below the finest grid of the wide format.
                self.m = Mantissa::min(self.mantissa_max_bit_len())?;
                self.e = EXPONENT_MIN;
            } else if let Some(w) = self.m.digits_mut().first_mut() {
                *w |= 1;
            }
        }

        self.inexact = false;
        self.set_precision(p, rm)?;
        let rounded_again = self.inexact;

        self.inexact = acc || rounded_again;
        self.op_inexact = Some(rounded || rounded_again);
        Ok(self)
    }

    fn set_precision_internal(
        &mut self,
        p: usize,
        rm: RoundingMode,
        mut check_roundable: bool,
        s: usize,
    ) -> Result<bool, Error> {
        Self::p_assertion(p)?;

        if self.mantissa_max_bit_len() > p && p > 0 {
            let mut inexact = self.inexact;

            let ovf = self.m.round_mantissa(
                self.mantissa_max_bit_len() - p,
                rm,
                self.is_positive(),
                &mut check_roundable,
                s,
                &mut inexact,
            );

            self.inexact |= inexact;
            self.op_inexact = None;

            if check_roundable && self.inexact {
                return Ok(false);
            }

            if ovf {
                if self.e == EXPONENT_MAX {
                    return Err(Error::ExponentOverflow(self.s));
                }

                self.e += 1;
            } else if self.m.is_all_zero() {
                self.m.set_bit_len(0);
                self.e = 0;
            } else if self.is_subnormal() {
                self.m.update_bit_len();
            }
        } else if p == 0 {
            if self.is_zero() {
                if check_roundable && self.inexact {
                    return Ok(false);
                }
            } else {
                self.inexact |= true;
                self.op_inexact = None;
            }
        }

        self.m.set_length(p)?;

        Ok(true)
    }

    /// Sets the sign of `self`.
    pub fn set_sign(&mut self, s: Sign) {
        self.s = s;
    }

    /// Inverts the sign of `self`.
    pub fn inv_sign(&mut self) {
        if self.is_negative() {
            self.set_sign(Sign::Pos);
        } else {
            self.set_sign(Sign::Neg);
        }
    }

    /// Create float from usize value.
    pub(crate) fn from_usize(u: usize) -> Result<Self, Error> {
        let (shift, m) = Mantissa::from_usize(u)?;
        let s = Sign::Pos;
        let e = (m.max_bit_len() - shift) as Exponent;

        Ok(ExactNumNumber {
            m,
            s,
            e,
            inexact: false,
            op_inexact: None,
        })
    }

    /// Returns the raw mantissa words of a number.
    pub fn mantissa_mut(&mut self) -> &mut Mantissa {
        &mut self.m
    }

    /// Returns the raw mantissa words of a number.
    pub fn mantissa(&self) -> &Mantissa {
        &self.m
    }

    /// True when the mantissa is stored without a heap allocation.
    #[inline]
    pub fn is_inline(&self) -> bool {
        self.m.is_inline()
    }

    /// Constructs ExactNumNumber with precision `p` from a signed integer value `i`.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - InvalidArgument: the precision is incorrect.
    pub fn from_i128(i: i128, p: usize) -> Result<Self, Error> {
        let sign = if i < 0 { Sign::Neg } else { Sign::Pos };
        let mut ret = Self::from_u128(i.unsigned_abs(), p)?;
        ret.set_sign(sign);
        Ok(ret)
    }

    /// Constructs ExactNumNumber with precision `p` from an unsigned integer value `u`.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - InvalidArgument: the precision is incorrect.        
    pub fn from_u128(mut v: u128, p: usize) -> Result<Self, Error> {
        const SZ: usize = core::mem::size_of::<u128>() * 8;

        if p < SZ {
            return Err(Error::InvalidArgument);
        }

        Self::p_assertion(p)?;

        if v == 0 {
            Self::new(p)
        } else {
            let mut shift = 0;
            while v <= (u128::MAX >> 1) {
                v <<= 1;
                shift += 1;
            }

            let mut words = [0; SZ.div_ceil(WORD_BIT_SIZE)];
            for w in &mut words {
                *w = v as Word;
                v >>= WORD_BIT_SIZE;
            }

            Ok(ExactNumNumber {
                m: Mantissa::from_words(p, &words)?,
                e: (SZ - shift) as Exponent,
                s: Sign::Pos,
                inexact: false,
                op_inexact: None,
            })
        }
    }

    // Add correction to x for flooring and ceiling rounding modes.
    pub(crate) fn add_correction(&self, inv_corr_sign: bool) -> Result<Self, Error> {
        let p = self.mantissa_max_bit_len() + 1;
        let mut corr = ExactNumNumber::min_positive(p)?;
        corr.set_exponent(self.exponent());
        corr.set_sign(if inv_corr_sign { self.sign().invert() } else { self.sign() });
        self.add(&corr, p, RoundingMode::None)
    }

    /// Divide `self` by 2.
    pub(crate) fn div_by_2(&mut self, rm: RoundingMode) {
        let e = self.exponent();
        if e == EXPONENT_MIN {
            self.subnormalize(e as i64 - 1, rm);
        } else {
            self.set_exponent(e - 1);
        }
    }

    /// Returns true if self is `inexact`.
    #[inline]
    pub fn inexact(&self) -> bool {
        self.inexact
    }

    /// Marks self as inexact if `inexact` is true, or exact otherwise.
    #[inline]
    pub fn set_inexact(&mut self, inexact: bool) {
        self.inexact = inexact;
        self.op_inexact = None;
    }

    /// A copy of `self` with both inexact flags cleared, so that a function evaluated on it
    /// reports only its own rounding.
    pub(crate) fn exact_copy(&self) -> Result<Self, Error> {
        let mut c = self.clone()?;
        c.set_inexact(false);
        Ok(c)
    }

    /// Records whether the operation that produced `self` rounded it.
    pub(crate) fn set_op_inexact(&mut self, op: bool) {
        self.op_inexact = Some(op);
    }

    /// Returns true if the operation that produced `self` rounded. Values whose producing
    /// operation does not record that report the accumulated flag, `inexact()`.
    #[inline]
    pub fn op_inexact(&self) -> bool {
        self.op_inexact.unwrap_or(self.inexact)
    }

    // Splits the flag after an operation that started from `inexact = false`: `self.inexact`
    // holds what this operation rounded; `inputs` is the operands' accumulated flag.
    pub(crate) fn finish_flags(&mut self, inputs: bool, per_op: bool) {
        let rounded = self.inexact;
        self.inexact = rounded || inputs;
        self.op_inexact = if per_op { Some(rounded) } else { None };
    }
}

macro_rules! impl_int_conv {
    ($s:ty, $u:ty, $from_s:ident, $from_u:ident) => {
        impl ExactNumNumber {
            /// Constructs ExactNumNumber with precision `p` from a signed integer value `i`.
            /// Precision is rounded upwards to the word size.
            ///
            /// ## Errors
            ///
            ///  - MemoryAllocation: failed to allocate memory for mantissa.
            ///  - InvalidArgument: the precision is incorrect.
            pub fn $from_s(i: $s, p: usize) -> Result<Self, Error> {
                let sign = if i < 0 { Sign::Neg } else { Sign::Pos };
                let mut ret = Self::$from_u(i.unsigned_abs(), p)?;
                ret.set_sign(sign);
                Ok(ret)
            }

            /// Constructs ExactNumNumber with precision `p` from an unsigned integer value `u`.
            /// Precision is rounded upwards to the word size.
            ///
            /// ## Errors
            ///
            ///  - MemoryAllocation: failed to allocate memory for mantissa.
            ///  - InvalidArgument: the precision is incorrect.
            pub fn $from_u(u: $u, p: usize) -> Result<Self, Error> {
                const SZ: usize = core::mem::size_of::<$u>() * 8;

                if p < SZ {
                    return Err(Error::InvalidArgument);
                }

                Self::from_u64_internal(u as u64, p)
            }
        }
    };
}

impl_int_conv!(i8, u8, from_i8, from_u8);
impl_int_conv!(i16, u16, from_i16, from_u16);
impl_int_conv!(i32, u32, from_i32, from_u32);
impl_int_conv!(i64, u64, from_i64, from_u64);

#[cfg(test)]
mod tests {

    use crate::{common::util::random_subnormal, defs::WORD_MAX, Consts};

    use super::*;
    use crate::common::test_rng::random;

    #[cfg(not(feature = "std"))]
    use alloc::format;

    #[test]
    fn test_number() -> Result<(), crate::Error> {
        let p_rng = 10;
        let p_min = 5;
        let p1 = (random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let p2 = (random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let p = (random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;

        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;

        let mut d1;
        let mut d2;
        let mut d3;
        let mut ref_num;
        let mut eps = ExactNumNumber::small_const(1);

        for _ in 0..1000 {
            let p = (random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;

            let i1: i8 = crate::common::test_rng::random::<i8>();
            let i2: i16 = crate::common::test_rng::random::<i16>();
            let i3: i32 = crate::common::test_rng::random::<i32>();
            let i4: i64 = crate::common::test_rng::random::<i64>();
            let i5: i128 = crate::common::test_rng::random::<i128>();

            let n1 = ExactNumNumber::from_i8(i1, p)?;
            let n2 = ExactNumNumber::from_i16(i2, p)?;
            let n3 = ExactNumNumber::from_i32(i3, p)?;
            let n4 = ExactNumNumber::from_i64(i4, p)?;
            let n5 = ExactNumNumber::from_i128(i5, p)?;

            let p1 = ExactNumNumber::parse(
                &format!("{}", i1),
                crate::Radix::Dec,
                p,
                RoundingMode::None,
                &mut cc,
            )?;
            let p2 = ExactNumNumber::parse(
                &format!("{}", i2),
                crate::Radix::Dec,
                p,
                RoundingMode::None,
                &mut cc,
            )?;
            let p3 = ExactNumNumber::parse(
                &format!("{}", i3),
                crate::Radix::Dec,
                p,
                RoundingMode::None,
                &mut cc,
            )?;
            let p4 = ExactNumNumber::parse(
                &format!("{}", i4),
                crate::Radix::Dec,
                p,
                RoundingMode::None,
                &mut cc,
            )?;
            let p5 = ExactNumNumber::parse(
                &format!("{}", i5),
                crate::Radix::Dec,
                p,
                RoundingMode::None,
                &mut cc,
            )?;

            assert!(p1.cmp(&n1) == 0);
            assert!(p2.cmp(&n2) == 0);
            assert!(p3.cmp(&n3) == 0);
            assert!(p4.cmp(&n4) == 0);
            assert!(p5.cmp(&n5) == 0);
        }

        // 0 * 0
        d1 = ExactNumNumber::new(p1)?;
        d2 = ExactNumNumber::new(p2)?;
        ref_num = ExactNumNumber::new(p)?;
        d3 = d1.mul(&d2, p, rm)?;
        assert!(d3.cmp(&ref_num) == 0);

        // 0 / d2
        d2 = ExactNumNumber::from_i8(1, p2)?; // non-zero divisor
        d3 = d1.div(&d2, p, rm)?;
        ref_num = ExactNumNumber::new(p)?;
        assert!(d3.cmp(&ref_num) == 0);

        // 0 / -d2
        d2 = d2.neg()?;
        d3 = d1.div(&d2, p, rm)?;
        assert!(d3.cmp(&ref_num) == 0);

        // add & sub & cmp
        for _ in 0..10000 {
            let p1 = (random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
            let p2 = (random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
            let p = (random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;

            d1 = ExactNumNumber::random_normal(p1, -(p1 as Exponent) / 2, p1 as Exponent)?;
            d2 = ExactNumNumber::random_normal(p2, -(p2 as Exponent) / 2, p2 as Exponent)?;

            let d3 = d1.sub(&d2, p, RoundingMode::ToEven)?;
            let d4 = d3.add(&d2, p, RoundingMode::ToEven)?;

            eps.set_exponent(d1.exponent().max(d2.exponent()) - p as Exponent + 2);

            assert!(d1.sub(&d4, p, RoundingMode::ToEven)?.abs()?.cmp(&eps) < 0);

            // subnormal
            d1 = random_subnormal(p1)?;
            d2 = random_subnormal(p2)?;

            let d3 = d1.sub(&d2, p, RoundingMode::ToEven)?;
            let d4 = d3.add(&d2, p, RoundingMode::ToEven)?;

            d1.set_precision(p, RoundingMode::ToEven)?;

            assert!(
                d1.sub(&d4, p, RoundingMode::ToEven)?
                    .abs()?
                    .cmp(&ExactNumNumber::min_positive(p)?)
                    <= 0
            );
        }

        // cancellation
        let w = WORD_MAX ^ 1 ^ 4 ^ 16;
        d1 = ExactNumNumber::from_raw_parts(
            &[WORD_MAX, WORD_MAX],
            WORD_BIT_SIZE * 2,
            Sign::Pos,
            0,
            false,
        )?;
        d2 = ExactNumNumber::from_raw_parts(
            &[w, WORD_MAX, WORD_MAX],
            WORD_BIT_SIZE * 3,
            Sign::Pos,
            0,
            false,
        )?;

        d3 = d1.sub(&d2, WORD_BIT_SIZE, RoundingMode::None)?;

        assert_eq!(w, d3.mantissa().digits()[0]);

        // increase precision
        d1 = ExactNumNumber::from_raw_parts(
            &[0, WORD_MAX - 10, WORD_MAX],
            WORD_BIT_SIZE * 3,
            Sign::Pos,
            0,
            false,
        )?;
        d2 = ExactNumNumber::from_raw_parts(
            &[0, WORD_MAX, WORD_MAX],
            WORD_BIT_SIZE * 3,
            Sign::Pos,
            -16,
            false,
        )?;

        d3 = d1.add(&d2, WORD_BIT_SIZE * 3, RoundingMode::ToEven)?;

        d1 = ExactNumNumber::from_raw_parts(
            &[WORD_MAX - 10, WORD_MAX],
            WORD_BIT_SIZE * 2,
            Sign::Pos,
            0,
            false,
        )?;
        d2 = ExactNumNumber::from_raw_parts(
            &[WORD_MAX, WORD_MAX],
            WORD_BIT_SIZE * 2,
            Sign::Pos,
            -16,
            false,
        )?;
        let d4 = d1.add(&d2, WORD_BIT_SIZE * 4, RoundingMode::ToEven)?;

        assert!(d4.mantissa_max_bit_len() == WORD_BIT_SIZE * 4);
        assert!(d3.cmp(&d4) == 0);

        d1 = ExactNumNumber::from_raw_parts(
            &[0, WORD_MAX - 10, WORD_MAX],
            WORD_BIT_SIZE * 3,
            Sign::Pos,
            0,
            false,
        )?;
        d2 = ExactNumNumber::from_raw_parts(
            &[0, WORD_MAX, WORD_MAX],
            WORD_BIT_SIZE * 3,
            Sign::Pos,
            -16,
            false,
        )?;

        d3 = d1.sub(&d2, WORD_BIT_SIZE * 3, RoundingMode::ToEven)?;

        d1 = ExactNumNumber::from_raw_parts(
            &[WORD_MAX - 10, WORD_MAX],
            WORD_BIT_SIZE * 2,
            Sign::Pos,
            0,
            false,
        )?;
        d2 = ExactNumNumber::from_raw_parts(
            &[WORD_MAX, WORD_MAX],
            WORD_BIT_SIZE * 2,
            Sign::Pos,
            -16,
            false,
        )?;
        let d4 = d1.sub(&d2, WORD_BIT_SIZE * 4, RoundingMode::ToEven)?;

        assert!(d4.mantissa_max_bit_len() == WORD_BIT_SIZE * 4);
        assert!(d3.cmp(&d4) == 0);

        // decrease precision
        d1 = ExactNumNumber::from_raw_parts(
            &[0, WORD_MAX - 10, WORD_MAX],
            WORD_BIT_SIZE * 3,
            Sign::Pos,
            0,
            false,
        )?;
        d2 = ExactNumNumber::from_raw_parts(
            &[0, WORD_MAX, WORD_MAX],
            WORD_BIT_SIZE * 3,
            Sign::Pos,
            -16,
            false,
        )?;
        d3 = d1.add(&d2, WORD_BIT_SIZE, RoundingMode::ToEven)?;
        let mut d4 = d1.add(&d2, WORD_BIT_SIZE * 3, RoundingMode::ToEven)?;
        d4.set_precision(WORD_BIT_SIZE, RoundingMode::ToEven)?;

        assert!(d4.mantissa_max_bit_len() == WORD_BIT_SIZE);
        assert!(d3.cmp(&d4) == 0);

        d3 = d1.sub(&d2, WORD_BIT_SIZE, RoundingMode::ToEven)?;
        let mut d4 = d1.sub(&d2, WORD_BIT_SIZE * 3, RoundingMode::ToEven)?;
        d4.set_precision(WORD_BIT_SIZE, RoundingMode::ToEven)?;

        assert!(d4.mantissa_max_bit_len() == WORD_BIT_SIZE);
        assert!(d3.cmp(&d4) == 0);

        // full prec
        for _ in 0..10000 {
            let p1 = (random::<usize>() % 3 + 1) * WORD_BIT_SIZE;
            let p2 = (random::<usize>() % 3 + 1) * WORD_BIT_SIZE;

            d1 = ExactNumNumber::random_normal(p1, -(p1 as Exponent) / 2, p1 as Exponent)?;
            d2 = ExactNumNumber::random_normal(p2, -(p2 as Exponent) / 2, p2 as Exponent)?;
            let d3 = d1.sub_full_prec(&d2)?;
            let d4 = d3.add_full_prec(&d2)?;
            //println!("\n=== res \n{:?} \n{:?} \n{:?} \n{:?} \n{:?} \n{:?}", d1, d2, d3, d4, d1.sub(&d4, RoundingMode::ToEven).unwrap().abs().unwrap(), eps);
            assert!(d1.cmp(&d4) == 0);
        }

        // mul & div
        for i in 0..10000 {
            let p1 = (random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
            let p2 = (random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
            let p = (random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;

            d1 = ExactNumNumber::random_normal(
                p1,
                EXPONENT_MIN / 2 + p as Exponent,
                EXPONENT_MAX / 2,
            )?;
            d2 = ExactNumNumber::random_normal(p2, EXPONENT_MIN / 2, EXPONENT_MAX / 2)?;

            if !d2.is_zero() {
                let d4 = if i & 1 == 0 {
                    let d3 = d1.div(&d2, p, RoundingMode::ToEven)?;
                    d3.mul(&d2, p, RoundingMode::ToEven)?
                } else {
                    let d3 = d1.mul(&d2, p, RoundingMode::ToEven)?;
                    d3.div(&d2, p, RoundingMode::ToEven)?
                };

                eps.set_exponent(d1.exponent() - p as Exponent + 2);

                //println!("\n{:?}\n{:?}\n{:?}\n{:?}", d1,d2,d3,d4);

                assert!(d1.sub(&d4, p, RoundingMode::ToEven)?.abs()?.cmp(&eps) < 0);
            }

            // subnormal
            d1 = random_subnormal(p1)?;
            d2 = random_subnormal(p2)?;

            let d3 = d1.div(&d2, p, RoundingMode::ToEven)?;
            let d4 = d3.mul(&d2, p, RoundingMode::ToEven)?;

            assert!(
                d1.sub(&d4, p, RoundingMode::ToEven)?
                    .abs()?
                    .cmp(&ExactNumNumber::min_positive(p)?)
                    <= 0
            );
        }

        // full prec
        for _ in 0..10000 {
            let p1 = (random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
            let p2 = (random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
            let p = (random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;

            d1 = ExactNumNumber::random_normal(
                p1,
                EXPONENT_MIN / 2 + p as Exponent, // avoid subnormal numbers
                EXPONENT_MAX / 2,
            )?;
            d2 = ExactNumNumber::random_normal(p2, EXPONENT_MIN / 2, EXPONENT_MAX / 2)?;

            if !d2.is_zero() {
                let d3 = d1.mul_full_prec(&d2)?;
                let d4 = d1.mul(&d2, p1 + p2, rm)?;
                //println!("\n{:?}\n{:?}\n{:?}\n{:?}", d1,d2,d3,d4);
                assert!(d3.cmp(&d4) == 0);
            }
        }

        // large mantissa
        let pp = 32000;
        for i in 0..20 {
            let p1 = (random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE + pp;
            let p2 = (random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE + pp;
            let p = (random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE + pp;

            d1 = ExactNumNumber::random_normal(
                p1,
                EXPONENT_MIN / 2 + p as Exponent, // avoid subnormal numbers
                EXPONENT_MAX / 2,
            )?;
            d2 = ExactNumNumber::random_normal(p2, EXPONENT_MIN / 2, EXPONENT_MAX / 2)?;

            if !d2.is_zero() {
                let d4 = if i & 1 == 0 {
                    let d3 = d1.div(&d2, p, RoundingMode::ToEven)?;
                    d3.mul(&d2, p, RoundingMode::ToEven)?
                } else {
                    let d3 = d1.mul(&d2, p, RoundingMode::ToEven)?;
                    d3.div(&d2, p, RoundingMode::ToEven)?
                };

                eps.set_exponent(d1.exponent() - p as Exponent + 2);

                assert!(d1.sub(&d4, p, RoundingMode::ToEven)?.abs()?.cmp(&eps) < 0);
            }
        }

        // variable precision
        d1 = ExactNumNumber::from_raw_parts(
            &[WORD_MAX - 10, WORD_MAX],
            WORD_BIT_SIZE * 2,
            Sign::Pos,
            0,
            false,
        )?;
        d2 = ExactNumNumber::from_raw_parts(
            &[WORD_MAX, WORD_MAX - 10],
            WORD_BIT_SIZE * 2,
            Sign::Pos,
            -16,
            false,
        )?;

        d3 = d1.mul_full_prec(&d2)?;

        let d4 = d1.mul(&d2, WORD_BIT_SIZE * 5, RoundingMode::ToEven)?;

        assert!(d4.mantissa_max_bit_len() == WORD_BIT_SIZE * 5);
        assert!(d3.cmp(&d4) == 0);

        let d4 = d1.mul(&d2, WORD_BIT_SIZE * 4, RoundingMode::ToEven)?;

        assert!(d4.mantissa_max_bit_len() == WORD_BIT_SIZE * 4);
        assert!(d3.cmp(&d4) == 0);

        d3.set_precision(WORD_BIT_SIZE * 3, RoundingMode::ToEven)?;
        let d4 = d1.mul(&d2, WORD_BIT_SIZE * 3, RoundingMode::ToEven)?;

        assert!(d4.mantissa_max_bit_len() == WORD_BIT_SIZE * 3);
        assert!(d3.cmp(&d4) == 0);

        d3.set_precision(WORD_BIT_SIZE, RoundingMode::ToEven)?;
        let d4 = d1.mul(&d2, WORD_BIT_SIZE, RoundingMode::ToEven)?;

        assert!(d4.mantissa_max_bit_len() == WORD_BIT_SIZE);
        assert!(d3.cmp(&d4) == 0);

        d1 = ExactNumNumber::from_i8(2, WORD_BIT_SIZE * 2)?;
        d2 = ExactNumNumber::from_i8(3, WORD_BIT_SIZE * 2)?;

        d3 = d1.div(&d2, WORD_BIT_SIZE * 5, RoundingMode::ToEven)?;

        let words = {
            #[cfg(not(target_pointer_width = "32"))]
            {
                [
                    12297829382473034411,
                    12297829382473034410,
                    12297829382473034410,
                    12297829382473034410,
                    12297829382473034410,
                ]
            }
            #[cfg(target_pointer_width = "32")]
            {
                [2863311531, 2863311530, 2863311530, 2863311530, 2863311530]
            }
        };

        assert!(d3.mantissa_max_bit_len() == WORD_BIT_SIZE * 5);
        assert!(d3.mantissa().digits() == words);

        d3 = d1.div(&d2, WORD_BIT_SIZE * 3, RoundingMode::ToEven)?;

        let words = {
            #[cfg(not(target_pointer_width = "32"))]
            {
                [12297829382473034411, 12297829382473034410, 12297829382473034410]
            }
            #[cfg(target_pointer_width = "32")]
            {
                [2863311531, 2863311530, 2863311530]
            }
        };

        assert!(d3.mantissa_max_bit_len() == WORD_BIT_SIZE * 3);
        assert!(d3.mantissa().digits() == words);

        d3 = d1.div(&d2, WORD_BIT_SIZE, RoundingMode::ToEven)?;

        let words = {
            #[cfg(not(target_pointer_width = "32"))]
            {
                [12297829382473034411]
            }
            #[cfg(target_pointer_width = "32")]
            {
                [2863311531]
            }
        };

        assert!(d3.mantissa_max_bit_len() == WORD_BIT_SIZE);
        assert!(d3.mantissa().digits() == words);

        // reciprocal
        for _ in 0..1000 {
            let p1 = (random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
            let p = (random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;

            d1 = ExactNumNumber::random_normal(
                p1,
                EXPONENT_MIN + p as Exponent, // avoid exponent overflow error
                EXPONENT_MAX,
            )?;

            if !d1.is_zero() {
                let d3 = d1.reciprocal(p, rm)?;
                let d4 = ONE.div(&d3, p, rm)?;

                eps.set_exponent(d1.exponent() - p as Exponent + 2);

                //println!("\n{:?}\n{:?}\n{:?}", d1, d3, d4);

                assert!(d1.sub(&d4, p, rm)?.abs()?.cmp(&eps) < 0);
            }
        }

        // reciprocal near 1
        d1 = ExactNumNumber::parse(
            "F.FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF2DC85F7E77EC487_e-1",
            crate::Radix::Hex,
            320,
            RoundingMode::None,
            &mut cc,
        )?;

        d2 = d1.reciprocal(320, RoundingMode::ToEven)?;
        d3 = ExactNumNumber::parse("1.00000000000000000000000000000000000000000000000000000000000000000D237A0818813B7A_e+0", crate::Radix::Hex, 320, RoundingMode::None, &mut cc)?;

        //println!("{:?}", d2.format(crate::Radix::Hex, rm).unwrap());

        assert!(d2.cmp(&d3) == 0);

        // MAX
        d1 = ExactNumNumber::max_value(p1)?;
        d2 = d1.reciprocal(p, rm)?;
        d3 = ONE.div(&d1, p, rm)?;
        assert!(d2.cmp(&d3) == 0);

        // variable precision
        d1 = ExactNumNumber::from_i8(3, WORD_BIT_SIZE * 2)?;
        d2 = d1.reciprocal(WORD_BIT_SIZE * 5, RoundingMode::ToEven)?;

        let words = {
            #[cfg(not(target_pointer_width = "32"))]
            {
                [
                    12297829382473034411,
                    12297829382473034410,
                    12297829382473034410,
                    12297829382473034410,
                    12297829382473034410,
                ]
            }
            #[cfg(target_pointer_width = "32")]
            {
                [2863311531, 2863311530, 2863311530, 2863311530, 2863311530]
            }
        };

        assert!(d2.mantissa().digits() == words);

        d2 = d1.reciprocal(WORD_BIT_SIZE, RoundingMode::ToEven)?;

        let words = {
            #[cfg(not(target_pointer_width = "32"))]
            {
                [12297829382473034411]
            }
            #[cfg(target_pointer_width = "32")]
            {
                [2863311531]
            }
        };

        assert!(d2.mantissa().digits() == words);

        // subnormal numbers basic sanity
        d1 = ExactNumNumber::min_positive(p)?;
        d2 = ExactNumNumber::min_positive(p)?;
        ref_num = ExactNumNumber::min_positive(p)?;
        let one = ExactNumNumber::from_word(1, p)?;
        ref_num = ref_num.mul(&one.add(&one, p, rm)?, p, rm)?;

        // min_positive + min_positive = 2*min_positive
        assert!(d1.add(&d2, p, rm)?.cmp(&ref_num) == 0);
        assert!(d1.add(&d2, p, rm)?.cmp(&d1) > 0);
        assert!(d1.cmp(&d1.add(&d2, p, rm)?) < 0);

        // min_positive - min_positive = 0
        ref_num = ExactNumNumber::new(p)?;
        assert!(d1.sub(&d2, p, rm)?.cmp(&ref_num) == 0);

        // 1 * min_positive = min_positive
        assert!(one.mul(&d2, p, rm)?.cmp(&d2) == 0);

        // min_positive / 1 = min_positive
        assert!(d2.div(&one, p, rm)?.cmp(&d2) == 0);

        // min_positive / 1 = min_positive
        assert!(d2.div(&one, p, rm)?.cmp(&d2) == 0);

        // normal -> subnormal -> normal
        d1 = ExactNumNumber::from_word(1, p)?;
        d1.e = EXPONENT_MIN;
        d2 = ExactNumNumber::min_positive(p)?;
        assert!(!d1.is_subnormal());
        assert!(d1.sub(&d2, p, rm)?.cmp(&d1) < 0);
        assert!(d1.cmp(&d1.sub(&d2, p, rm)?) > 0);
        d1 = d1.sub(&d2, p, rm)?;
        assert!(d1.is_subnormal());
        d1 = d1.add(&d2, p, rm)?;
        assert!(!d1.is_subnormal());

        // overflow basic sanity
        d1 = ExactNumNumber::from_word(1, p)?;
        d1.e = EXPONENT_MAX - (d1.m.max_bit_len() - 1) as Exponent;
        assert!(matches!(
            ExactNumNumber::max_value(p)?.add(&d1, p, rm),
            Err(Error::ExponentOverflow(Sign::Pos))
        ));
        assert!(matches!(
            ExactNumNumber::min_value(p)?.sub(&d1, p, rm),
            Err(Error::ExponentOverflow(Sign::Neg))
        ));
        assert!(matches!(
            ExactNumNumber::max_value(p)?.mul(&ExactNumNumber::max_value(p)?, p, rm),
            Err(Error::ExponentOverflow(Sign::Pos))
        ));
        d1 = ExactNumNumber::from_word(1, p)?;
        d1.e = EXPONENT_MIN;
        assert!(matches!(
            ExactNumNumber::max_value(p)?.div(&d1, p, rm),
            Err(Error::ExponentOverflow(Sign::Pos))
        ));

        // decompose and compose
        let (m, n, s, e, inexact) = d1.as_raw_parts();
        d2 = ExactNumNumber::from_raw_parts(m, n, s, e, inexact)?;
        assert!(d1.cmp(&d2) == 0);
        assert!(ExactNumNumber::from_raw_parts(
            &[1, WORD_SIGNIFICANT_BIT],
            WORD_BIT_SIZE * 3,
            Sign::Pos,
            123,
            false,
        )
        .is_err());
        assert!(ExactNumNumber::from_raw_parts(
            &[1, WORD_SIGNIFICANT_BIT],
            WORD_BIT_SIZE * 2,
            Sign::Pos,
            123,
            false,
        )
        .is_ok());
        assert!(ExactNumNumber::from_raw_parts(
            &[1, WORD_SIGNIFICANT_BIT],
            WORD_BIT_SIZE,
            Sign::Pos,
            123,
            false,
        )
        .is_err());
        assert!(ExactNumNumber::from_raw_parts(&[0, 0], 0, Sign::Pos, 0, false).is_ok());
        assert!(ExactNumNumber::from_raw_parts(&[1, 0], 0, Sign::Pos, 0, false).is_err());
        assert!(ExactNumNumber::from_raw_parts(
            &[1, 2],
            WORD_BIT_SIZE,
            Sign::Pos,
            EXPONENT_MIN,
            false
        )
        .is_err());
        assert!(
            ExactNumNumber::from_raw_parts(&[1, 2], WORD_BIT_SIZE + 2, Sign::Pos, 123, false)
                .is_err()
        );
        assert!(ExactNumNumber::from_raw_parts(
            &[1, 2],
            WORD_BIT_SIZE + 2,
            Sign::Pos,
            EXPONENT_MIN,
            false
        )
        .is_ok());

        // sign and exponent
        d1 = ExactNumNumber::from_word(1, p)?;
        assert!(d1.sign() == Sign::Pos);
        assert!(d1.is_positive());
        d1 = d1.neg()?;
        assert!(d1.sign() == Sign::Neg);
        assert!(d1.is_negative());
        assert!(d1.exponent() == 1);

        assert!(
            ExactNumNumber::min_positive(p)?
                .fract()?
                .cmp(&ExactNumNumber::min_positive(p)?)
                == 0
        );
        assert!(ExactNumNumber::min_positive(p)?.int()?.is_zero());

        d1 = ExactNumNumber::new(p)?;
        assert!(d1.fract()?.is_zero());
        assert!(ExactNumNumber::new(p)?.int()?.is_zero());

        // ceil & floor
        // rem
        for (prec1, prec2) in [(128, 128), (128, 320), (320, 128)] {
            for sign_inv in [false, true] {
                for (s1, s2) in [
                    ("9.2", "1.2"),
                    ("2.0", "4.0"),
                    ("7.0", "4.0"),
                    ("8.0", "4.0"),
                    ("9.0", "4.0"),
                    ("14", "15"),
                    ("15", "15"),
                    ("16", "15"),
                    ("33", "15"),
                    ("44", "15"),
                    ("0.000000007", "0.000000004"),
                    ("0.00000000000000009", "0.000000000000000004"),
                ] {
                    d1 = ExactNumNumber::parse(
                        s1,
                        crate::Radix::Dec,
                        prec1,
                        RoundingMode::None,
                        &mut cc,
                    )?;
                    d2 = ExactNumNumber::parse(
                        s2,
                        crate::Radix::Dec,
                        prec2,
                        RoundingMode::None,
                        &mut cc,
                    )?;

                    let p = prec1.max(prec2);

                    d3 = d1.sub(
                        &d1.div(&d2, p + 1, RoundingMode::ToEven)?.floor()?.mul(
                            &d2,
                            p + 1,
                            RoundingMode::ToEven,
                        )?,
                        p + 1,
                        RoundingMode::ToEven,
                    )?;

                    if sign_inv {
                        d1.inv_sign();
                        d3.inv_sign();
                    }

                    let ret = d1.rem(&d2)?;

                    // println!("\n{:?}", d1.format(crate::Radix::Dec, RoundingMode::None).unwrap());
                    // println!("{:?}", d2.format(crate::Radix::Dec, RoundingMode::None).unwrap());
                    // println!("{:?}", ret.format(crate::Radix::Dec, RoundingMode::None).unwrap());
                    // println!("{:?}", d3.format(crate::Radix::Dec, RoundingMode::None).unwrap());

                    assert!(ret.sub(&d3, p, RoundingMode::ToEven)?.is_zero());
                }
            }
        }

        for _ in 0..1000 {
            let p1 = (random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
            let p2 = (random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
            let p = p1.max(p2);
            let e_rng = (p_rng * WORD_BIT_SIZE) as Exponent;

            d1 = ExactNumNumber::random_normal(
                p1,
                EXPONENT_MIN / 2 - p1 as Exponent,
                EXPONENT_MAX / 2,
            )?
            .abs()?;
            d2 = ExactNumNumber::random_normal(p2, d1.exponent() - e_rng, d1.exponent() + e_rng)?
                .abs()?;

            let ret = d1.rem(&d2)?;

            d3 = d1.sub(
                &d1.div(&d2, p + p_rng * WORD_BIT_SIZE, RoundingMode::None)?
                    .floor()?
                    .mul(&d2, p + p_rng * WORD_BIT_SIZE, RoundingMode::None)?,
                ret.mantissa_max_bit_len(),
                RoundingMode::ToEven,
            )?;

            //println!("\n{:?}", d1.format(crate::Radix::Dec, RoundingMode::None).unwrap());
            //println!("{:?}", d2.format(crate::Radix::Dec, RoundingMode::None).unwrap());
            //println!("{:?}", ret.format(crate::Radix::Dec, RoundingMode::None).unwrap());
            //println!("{:?}", d3.format(crate::Radix::Dec, RoundingMode::None).unwrap());

            assert!(ret.sub(&d3, p, RoundingMode::ToEven)?.is_zero());

            // subnormal
            d1 = random_subnormal(p1)?;
            d2 = random_subnormal(p2)?;

            let mut ret = d1.rem(&d2)?;

            // check with normalized numbers.
            d1.set_exponent(0);
            d2.set_exponent(0);

            let mut d3 = d1.rem(&d2)?;

            if !ret.is_zero() {
                // check exponent
                assert!(ret.exponent() == EXPONENT_MIN);
            }

            eps.set_exponent(1 - (ret.precision() as Exponent));

            // compare mantissas
            ret.set_exponent((ret.mantissa_max_bit_len() - ret.precision()) as Exponent);
            d3.set_exponent(0);

            assert!(ret.sub(&d3, p, RoundingMode::ToEven)?.abs()?.cmp(&eps) <= 0);
        }

        // special case rem
        let d1 = ExactNumNumber::from_words(&[WORD_SIGNIFICANT_BIT], Sign::Pos, 2)?;
        let d2 = ExactNumNumber::from_words(&[WORD_SIGNIFICANT_BIT], Sign::Pos, 65)?;

        let d3 = d2.rem(&d1)?;
        assert!(d3.is_zero());

        // is_odd_int
        let d1 =
            ExactNumNumber::parse("3.0", crate::Radix::Dec, 128, RoundingMode::ToEven, &mut cc)?;
        assert!(d1.is_odd_int());
        let d1 = ExactNumNumber::parse(
            "3.01",
            crate::Radix::Dec,
            128,
            RoundingMode::ToEven,
            &mut cc,
        )?;
        assert!(!d1.is_odd_int());
        let d1 = ExactNumNumber::parse(
            "32.0",
            crate::Radix::Dec,
            128,
            RoundingMode::ToEven,
            &mut cc,
        )?;
        assert!(!d1.is_odd_int());
        let d1 = ExactNumNumber::parse(
            "32.01",
            crate::Radix::Dec,
            128,
            RoundingMode::ToEven,
            &mut cc,
        )?;
        assert!(!d1.is_odd_int());
        let d1 = ExactNumNumber::parse(
            "0.00000000000000000000001",
            crate::Radix::Dec,
            256,
            RoundingMode::ToEven,
            &mut cc,
        )?;
        assert!(!d1.is_odd_int());
        let d1 = ExactNumNumber::parse(
            "5.00000000000000000000001",
            crate::Radix::Dec,
            256,
            RoundingMode::ToEven,
            &mut cc,
        )?;
        assert!(!d1.is_odd_int());
        let d1 = ExactNumNumber::parse(
            "10000000000000000000001.0",
            crate::Radix::Dec,
            256,
            RoundingMode::ToEven,
            &mut cc,
        )?;
        assert!(d1.is_odd_int());
        let d1 = ExactNumNumber::parse(
            "10000000000000000000001.0000000000000000001",
            crate::Radix::Dec,
            256,
            RoundingMode::ToEven,
            &mut cc,
        )?;
        assert!(!d1.is_odd_int());
        let d1 = ExactNumNumber::parse("1.1", crate::Radix::Bin, 256, RoundingMode::None, &mut cc)?;
        assert!(!d1.is_odd_int());

        // build from words
        // An empty mantissa is precision 0, which is rejected (2.0: `MIN_PRECISION` is 1).
        assert!(matches!(
            ExactNumNumber::from_words(&[], Sign::Pos, EXPONENT_MAX),
            Err(Error::InvalidArgument)
        ));

        let d1 = ExactNumNumber::from_words(&[0, 0], Sign::Neg, EXPONENT_MIN)?;
        assert!(d1.is_zero());
        assert_eq!(d1.exponent(), 0);
        assert_eq!(d1.precision(), 0);
        assert_eq!(d1.mantissa_max_bit_len(), WORD_BIT_SIZE * 2);
        assert_eq!(d1.sign(), Sign::Neg);

        let d1 = ExactNumNumber::from_words(&[3, 1], Sign::Pos, EXPONENT_MAX)?;
        #[cfg(not(target_pointer_width = "32"))]
        {
            assert_eq!(
                d1.mantissa().digits(),
                [0x8000000000000000u64, 0x8000000000000001u64]
            );
        }
        #[cfg(target_pointer_width = "32")]
        {
            assert_eq!(d1.mantissa().digits(), [0x80000000u32, 0x80000001u32]);
        }
        assert_eq!(d1.exponent(), EXPONENT_MAX - WORD_BIT_SIZE as Exponent + 1);
        assert_eq!(d1.precision(), WORD_BIT_SIZE * 2);
        assert_eq!(d1.mantissa_max_bit_len(), WORD_BIT_SIZE * 2);
        assert_eq!(d1.sign(), Sign::Pos);

        let d1 = ExactNumNumber::from_words(&[3, 1], Sign::Neg, EXPONENT_MIN)?;
        assert_eq!(d1.mantissa().digits(), [3, 1]);
        assert_eq!(d1.exponent(), EXPONENT_MIN);
        assert_eq!(d1.precision(), WORD_BIT_SIZE + 1);
        assert_eq!(d1.mantissa_max_bit_len(), WORD_BIT_SIZE * 2);
        assert_eq!(d1.sign(), Sign::Neg);

        let d1 = ExactNumNumber::from_words(&[3, 1], Sign::Pos, EXPONENT_MIN + 5)?;
        assert_eq!(d1.mantissa().digits(), [3 << 5, 1 << 5]);
        assert_eq!(d1.exponent(), EXPONENT_MIN);
        assert_eq!(d1.precision(), WORD_BIT_SIZE + 6);
        assert_eq!(d1.mantissa_max_bit_len(), WORD_BIT_SIZE * 2);
        assert_eq!(d1.sign(), Sign::Pos);

        let words = {
            #[cfg(not(target_pointer_width = "32"))]
            {
                [3, 0x8000000000000000u64]
            }
            #[cfg(target_pointer_width = "32")]
            {
                [3, 0x80000000u32]
            }
        };
        let d1 = ExactNumNumber::from_words(&words, Sign::Pos, EXPONENT_MIN + 5)?;
        assert_eq!(d1.mantissa().digits(), words);
        assert_eq!(d1.exponent(), EXPONENT_MIN + 5);
        assert_eq!(d1.precision(), WORD_BIT_SIZE * 2);
        assert_eq!(d1.mantissa_max_bit_len(), WORD_BIT_SIZE * 2);
        assert_eq!(d1.sign(), Sign::Pos);
        Ok(())
    }

    #[test]
    fn test_rounding() -> Result<(), crate::Error> {
        // trailing bits
        // 0 0000
        // 0 0001
        // 0 1000
        // 0 1001
        // 1 1000
        // 1 1001

        let mantissas = {
            #[cfg(not(target_pointer_width = "32"))]
            {
                [
                    [0x8000000000000000u64, 0x8000000000000000u64],
                    [0x8000000000000001u64, 0x8000000000000000u64],
                    [0x8000000000000008u64, 0x8000000000000000u64],
                    [0x8000000000000009u64, 0x8000000000000000u64],
                    [0x8000000000000018u64, 0x8000000000000000u64],
                    [0x8000000000000019u64, 0x8000000000000000u64],
                ]
            }
            #[cfg(target_pointer_width = "32")]
            {
                [
                    [0x80000000u32, 0x80000000u32],
                    [0x80000001u32, 0x80000000u32],
                    [0x80000008u32, 0x80000000u32],
                    [0x80000009u32, 0x80000000u32],
                    [0x80000018u32, 0x80000000u32],
                    [0x80000019u32, 0x80000000u32],
                ]
            }
        };

        let rounding_results_posnum = {
            #[cfg(not(target_pointer_width = "32"))]
            {
                [
                    (RoundingMode::None, mantissas),
                    (
                        RoundingMode::Down,
                        [
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                        ],
                    ),
                    (
                        RoundingMode::Up,
                        [
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000020u64, 0x8000000000000000u64],
                            [0x8000000000000020u64, 0x8000000000000000u64],
                        ],
                    ),
                    (
                        RoundingMode::FromZero,
                        [
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000020u64, 0x8000000000000000u64],
                            [0x8000000000000020u64, 0x8000000000000000u64],
                        ],
                    ),
                    (
                        RoundingMode::ToZero,
                        [
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                        ],
                    ),
                    (
                        RoundingMode::ToEven,
                        [
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000020u64, 0x8000000000000000u64],
                            [0x8000000000000020u64, 0x8000000000000000u64],
                        ],
                    ),
                    (
                        RoundingMode::ToOdd,
                        [
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000020u64, 0x8000000000000000u64],
                        ],
                    ),
                ]
            }
            #[cfg(target_pointer_width = "32")]
            {
                [
                    (RoundingMode::None, mantissas),
                    (
                        RoundingMode::Down,
                        [
                            [0x80000000u32, 0x80000000u32],
                            [0x80000000u32, 0x80000000u32],
                            [0x80000000u32, 0x80000000u32],
                            [0x80000000u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                        ],
                    ),
                    (
                        RoundingMode::Up,
                        [
                            [0x80000000u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000020u32, 0x80000000u32],
                            [0x80000020u32, 0x80000000u32],
                        ],
                    ),
                    (
                        RoundingMode::FromZero,
                        [
                            [0x80000000u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000020u32, 0x80000000u32],
                            [0x80000020u32, 0x80000000u32],
                        ],
                    ),
                    (
                        RoundingMode::ToZero,
                        [
                            [0x80000000u32, 0x80000000u32],
                            [0x80000000u32, 0x80000000u32],
                            [0x80000000u32, 0x80000000u32],
                            [0x80000000u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                        ],
                    ),
                    (
                        RoundingMode::ToEven,
                        [
                            [0x80000000u32, 0x80000000u32],
                            [0x80000000u32, 0x80000000u32],
                            [0x80000000u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000020u32, 0x80000000u32],
                            [0x80000020u32, 0x80000000u32],
                        ],
                    ),
                    (
                        RoundingMode::ToOdd,
                        [
                            [0x80000000u32, 0x80000000u32],
                            [0x80000000u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000020u32, 0x80000000u32],
                        ],
                    ),
                ]
            }
        };

        let rounding_results_negnum = {
            #[cfg(not(target_pointer_width = "32"))]
            {
                [
                    (RoundingMode::None, mantissas),
                    (
                        RoundingMode::Down,
                        [
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000020u64, 0x8000000000000000u64],
                            [0x8000000000000020u64, 0x8000000000000000u64],
                        ],
                    ),
                    (
                        RoundingMode::Up,
                        [
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                        ],
                    ),
                    (
                        RoundingMode::FromZero,
                        [
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000020u64, 0x8000000000000000u64],
                            [0x8000000000000020u64, 0x8000000000000000u64],
                        ],
                    ),
                    (
                        RoundingMode::ToZero,
                        [
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                        ],
                    ),
                    (
                        RoundingMode::ToEven,
                        [
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000020u64, 0x8000000000000000u64],
                            [0x8000000000000020u64, 0x8000000000000000u64],
                        ],
                    ),
                    (
                        RoundingMode::ToOdd,
                        [
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000000u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000010u64, 0x8000000000000000u64],
                            [0x8000000000000020u64, 0x8000000000000000u64],
                        ],
                    ),
                ]
            }
            #[cfg(target_pointer_width = "32")]
            {
                [
                    (RoundingMode::None, mantissas),
                    (
                        RoundingMode::Down,
                        [
                            [0x80000000u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000020u32, 0x80000000u32],
                            [0x80000020u32, 0x80000000u32],
                        ],
                    ),
                    (
                        RoundingMode::Up,
                        [
                            [0x80000000u32, 0x80000000u32],
                            [0x80000000u32, 0x80000000u32],
                            [0x80000000u32, 0x80000000u32],
                            [0x80000000u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                        ],
                    ),
                    (
                        RoundingMode::FromZero,
                        [
                            [0x80000000u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000020u32, 0x80000000u32],
                            [0x80000020u32, 0x80000000u32],
                        ],
                    ),
                    (
                        RoundingMode::ToZero,
                        [
                            [0x80000000u32, 0x80000000u32],
                            [0x80000000u32, 0x80000000u32],
                            [0x80000000u32, 0x80000000u32],
                            [0x80000000u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                        ],
                    ),
                    (
                        RoundingMode::ToEven,
                        [
                            [0x80000000u32, 0x80000000u32],
                            [0x80000000u32, 0x80000000u32],
                            [0x80000000u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000020u32, 0x80000000u32],
                            [0x80000020u32, 0x80000000u32],
                        ],
                    ),
                    (
                        RoundingMode::ToOdd,
                        [
                            [0x80000000u32, 0x80000000u32],
                            [0x80000000u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000010u32, 0x80000000u32],
                            [0x80000020u32, 0x80000000u32],
                        ],
                    ),
                ]
            }
        };

        for (sign, rr) in
            [(Sign::Pos, rounding_results_posnum), (Sign::Neg, rounding_results_negnum)]
        {
            for (rm, expected_mantissas) in rr.iter() {
                for (m1, m2) in mantissas.iter().zip(expected_mantissas.iter()) {
                    // rounding
                    let d1 = ExactNumNumber::from_raw_parts(
                        m1,
                        WORD_BIT_SIZE * 2,
                        sign,
                        WORD_BIT_SIZE as Exponent,
                        false,
                    )?;
                    let d2 = d1.round(WORD_BIT_SIZE - 4, *rm)?;
                    let d3 = ExactNumNumber::from_raw_parts(
                        m2,
                        WORD_BIT_SIZE * 2,
                        sign,
                        WORD_BIT_SIZE as Exponent,
                        false,
                    )?;

                    //println!("\n{:?} {:?}\n{:?}\nresult {:?}\nexpect {:?}", sign, rm, d1.m, d2.m, d3.m);

                    assert!(d2.cmp(&d3) == 0);
                }
            }
        }

        // special cases
        let testset = [
            (
                [0, WORD_SIGNIFICANT_BIT],
                Sign::Pos,
                -99,
                99,
                RoundingMode::None,
                [0, 0],
                0,
            ),
            (
                [0, WORD_SIGNIFICANT_BIT],
                Sign::Pos,
                -99,
                99,
                RoundingMode::Up,
                [0, WORD_SIGNIFICANT_BIT],
                -98,
            ),
            (
                [0, WORD_SIGNIFICANT_BIT],
                Sign::Pos,
                -99,
                99,
                RoundingMode::ToOdd,
                [0, WORD_SIGNIFICANT_BIT],
                -98,
            ),
            (
                [0, WORD_SIGNIFICANT_BIT],
                Sign::Pos,
                -99,
                99,
                RoundingMode::ToEven,
                [0, 0],
                0,
            ),
            (
                [1, WORD_SIGNIFICANT_BIT],
                Sign::Pos,
                -99,
                99,
                RoundingMode::ToEven,
                [0, WORD_SIGNIFICANT_BIT],
                -98,
            ),
            (
                [1, WORD_SIGNIFICANT_BIT],
                Sign::Pos,
                -99,
                99,
                RoundingMode::Down,
                [0, 0],
                0,
            ),
            (
                [0, WORD_SIGNIFICANT_BIT],
                Sign::Neg,
                -99,
                99,
                RoundingMode::Down,
                [0, WORD_SIGNIFICANT_BIT],
                -98,
            ),
            (
                [1, WORD_SIGNIFICANT_BIT],
                Sign::Pos,
                -100,
                99,
                RoundingMode::ToEven,
                [0, 0],
                0,
            ),
            (
                [1, WORD_SIGNIFICANT_BIT],
                Sign::Pos,
                -100,
                99,
                RoundingMode::ToOdd,
                [0, 0],
                0,
            ),
            (
                [1, WORD_SIGNIFICANT_BIT],
                Sign::Pos,
                -100,
                99,
                RoundingMode::Down,
                [0, 0],
                0,
            ),
            (
                [0, WORD_SIGNIFICANT_BIT],
                Sign::Pos,
                -100,
                99,
                RoundingMode::Up,
                [0, WORD_SIGNIFICANT_BIT],
                -98,
            ),
            (
                [1, 1],
                Sign::Pos,
                EXPONENT_MIN,
                EXPONENT_MIN.unsigned_abs() as usize + WORD_BIT_SIZE - 1,
                RoundingMode::Up,
                [0, 2],
                EXPONENT_MIN,
            ),
            (
                [1, 1],
                Sign::Pos,
                EXPONENT_MIN,
                EXPONENT_MIN.unsigned_abs() as usize + WORD_BIT_SIZE - 1,
                RoundingMode::Down,
                [0, 0],
                0,
            ),
            (
                [1, 1],
                Sign::Pos,
                EXPONENT_MIN,
                EXPONENT_MIN.unsigned_abs() as usize + WORD_BIT_SIZE - 1,
                RoundingMode::ToEven,
                [0, 2],
                EXPONENT_MIN,
            ),
            (
                [0, 1],
                Sign::Pos,
                EXPONENT_MIN,
                EXPONENT_MIN.unsigned_abs() as usize + WORD_BIT_SIZE - 1,
                RoundingMode::ToEven,
                [0, 0],
                0,
            ),
            (
                [1, 1],
                Sign::Pos,
                EXPONENT_MIN,
                EXPONENT_MIN.unsigned_abs() as usize + WORD_BIT_SIZE - 1,
                RoundingMode::ToOdd,
                [0, 2],
                EXPONENT_MIN,
            ),
            (
                [0, 1],
                Sign::Pos,
                EXPONENT_MIN,
                EXPONENT_MIN.unsigned_abs() as usize + WORD_BIT_SIZE - 1,
                RoundingMode::ToOdd,
                [0, 2],
                EXPONENT_MIN,
            ),
        ];
        for (m1, s, e1, n, rm, m2, e2) in testset {
            //println!("{:?}", (m1, s, e1, n, rm, m2, e2));
            let d1 = ExactNumNumber::from_words(&m1, s, e1)?;
            let d1_inexact = d1.inexact();
            let d2 = d1.round(n, rm)?;

            let m2m = Mantissa::from_words(m2.len() * WORD_BIT_SIZE, &m2)?;
            assert_eq!(d2.mantissa().digits(), m2m.digits());
            assert_eq!(d2.mantissa().bit_len(), m2m.bit_len());
            assert_eq!(d2.sign(), s);
            assert_eq!(d2.exponent(), e2);
            assert!(!d1_inexact);
            assert!(d2.inexact());
        }
        Ok(())
    }

    #[test]
    fn test_inexact() -> Result<(), crate::Error> {
        // any arg is inexact
        let mut d1 = ExactNumNumber::from_words(&[0, 18 << (WORD_BIT_SIZE - 5)], Sign::Pos, 5)?;
        let mut d2 = ExactNumNumber::from_words(&[0, 6 << (WORD_BIT_SIZE - 3)], Sign::Pos, 3)?;

        for op in [
            ExactNumNumber::add,
            ExactNumNumber::sub,
            ExactNumNumber::mul,
            ExactNumNumber::div,
        ] {
            let d3 = op(&d1, &d2, 256, RoundingMode::None)?;
            assert!(!d3.inexact());

            d1.set_inexact(true);
            let d3 = op(&d1, &d2, 256, RoundingMode::None)?;
            assert!(d3.inexact());

            d2.set_inexact(true);
            let d3 = op(&d1, &d2, 256, RoundingMode::None)?;
            assert!(d3.inexact());

            d1.set_inexact(false);
            let d3 = op(&d1, &d2, 256, RoundingMode::None)?;
            assert!(d3.inexact());

            d2.set_inexact(false);
        }

        let d3 = d1.rem(&d2)?;
        assert!(!d3.inexact());

        d1.set_inexact(true);
        let d3 = d1.rem(&d2)?;
        assert!(d3.inexact());

        d2.set_inexact(true);
        let d3 = d1.rem(&d2)?;
        assert!(d3.inexact());

        d1.set_inexact(false);
        let d3 = d1.rem(&d2)?;
        assert!(d3.inexact());

        // inexact because of rounding
        let d1 = ExactNumNumber::from_words(&[2, 18 << (WORD_BIT_SIZE - 5)], Sign::Pos, 5)?;
        let d2 = ExactNumNumber::from_words(&[1, 6 << (WORD_BIT_SIZE - 3)], Sign::Pos, 3)?;

        assert!(!d1.inexact());
        assert!(!d2.inexact());

        for op in [
            ExactNumNumber::add,
            ExactNumNumber::sub,
            ExactNumNumber::mul,
            ExactNumNumber::div,
        ] {
            let d3 = op(&d1, &d2, WORD_BIT_SIZE * 2, RoundingMode::None)?;
            assert!(d3.inexact());
            let d3 = op(&d2, &d1, WORD_BIT_SIZE * 2, RoundingMode::None)?;
            assert!(d3.inexact());
        }

        // inexact because of subnormalization
        let d1 = ExactNumNumber::from_words(
            &[1, WORD_SIGNIFICANT_BIT],
            Sign::Pos,
            EXPONENT_MIN - EXPONENT_MIN / 2 - 1,
        )?;
        let d2 =
            ExactNumNumber::from_words(&[1, WORD_SIGNIFICANT_BIT], Sign::Pos, EXPONENT_MIN / 2)?;

        assert!(!d1.inexact());
        assert!(!d2.inexact());

        let d3 = d1.mul(&d2, WORD_BIT_SIZE * 4, RoundingMode::None)?;
        assert!(d3.inexact());

        let d1 = ExactNumNumber::from_words(&[1, WORD_SIGNIFICANT_BIT], Sign::Pos, EXPONENT_MIN)?;
        let d2 = ExactNumNumber::from_words(&[1, WORD_SIGNIFICANT_BIT], Sign::Pos, EXPONENT_MIN)?;

        assert!(!d1.inexact());
        assert!(!d2.inexact());

        let d3 = d1.mul(&d2, 256, RoundingMode::Up)?;
        assert!(d3.inexact());

        let d3 = d1.mul(&d2, 256, RoundingMode::None)?;
        assert!(d3.inexact());

        // set precision
        let mut d1 = ExactNumNumber::from_words(&[0, WORD_MAX], Sign::Pos, 0)?;
        assert!(!d1.inexact());
        d1.set_precision(WORD_BIT_SIZE, RoundingMode::ToEven)?;
        assert!(!d1.inexact());

        let mut d1 = ExactNumNumber::from_words(&[0, WORD_MAX], Sign::Pos, 0)?;
        d1.set_inexact(true);
        assert!(d1.inexact());
        d1.set_precision(WORD_BIT_SIZE, RoundingMode::ToEven)?;
        assert!(d1.inexact());

        let mut d1 = ExactNumNumber::from_words(&[WORD_SIGNIFICANT_BIT, WORD_MAX], Sign::Pos, 0)?;
        assert!(!d1.inexact());
        d1.set_precision(WORD_BIT_SIZE, RoundingMode::ToEven)?;
        assert!(d1.inexact());

        // round
        let d1 = ExactNumNumber::from_words(&[0, WORD_MAX], Sign::Pos, 0)?;
        assert!(!d1.inexact());
        let d2 = d1.round(WORD_BIT_SIZE, RoundingMode::ToEven)?;
        assert!(!d2.inexact());

        let mut d1 = ExactNumNumber::from_words(&[0, WORD_MAX], Sign::Pos, 0)?;
        d1.set_inexact(true);
        assert!(d1.inexact());
        let d2 = d1.round(WORD_BIT_SIZE, RoundingMode::ToEven)?;
        assert!(d2.inexact());

        let d1 = ExactNumNumber::from_words(&[WORD_SIGNIFICANT_BIT, WORD_MAX], Sign::Pos, 0)?;
        assert!(!d1.inexact());
        let d2 = d1.round(WORD_BIT_SIZE, RoundingMode::None)?;
        assert!(d2.inexact());

        let d1 = ExactNumNumber::from_words(&[0, WORD_MAX], Sign::Pos, -2)?;
        assert!(!d1.inexact());
        let d2 = d1.round(1, RoundingMode::None)?;
        assert!(d2.inexact());
        Ok(())
    }
}
