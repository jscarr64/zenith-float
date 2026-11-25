// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! Mantissa of a number.

use crate::common::buf::{WordBuf, INLINE_WORDS};
use crate::common::util::add_carry;
use crate::common::util::find_one_from;
use crate::common::util::shift_slice_left;
use crate::common::util::shift_slice_right;
use crate::common::util::sub_borrow;
use crate::defs::DoubleWord;
use crate::defs::Error;
use crate::defs::RoundingMode;
use crate::defs::SignedWord;
use crate::defs::Word;
use crate::defs::WORD_BASE;
use crate::defs::WORD_BIT_SIZE;
use crate::defs::WORD_MAX;
use crate::defs::WORD_SIGNIFICANT_BIT;
use crate::mantissa::util::ExtendedSlice;
use crate::mantissa::util::NormalizedView;
use crate::mantissa::util::RightShiftedSlice;
use core::mem::size_of;
use itertools::izip;

/// Mantissa representation.
#[derive(Debug, Hash)]
pub struct Mantissa {
    m: WordBuf,
    n: usize, // number of bits, 0 is for number 0
}

/// Parameters of [`Mantissa::abs_add`] and [`Mantissa::abs_sub`].
#[derive(Copy, Clone, Debug)]
pub struct AddSubSpec {
    /// Result precision in bits.
    pub p: usize,
    /// Right shift of the second operand, in bits.
    pub m2_shift: usize,
    /// Rounding mode.
    pub rm: RoundingMode,
    /// Sign of the result.
    pub is_positive: bool,
    /// Keep full precision instead of rounding to `p`.
    pub full_prec: bool,
}

impl Mantissa {
    // bit lenth to length in words.
    #[inline]
    fn bit_len_to_word_len(p: usize) -> usize {
        (p + WORD_BIT_SIZE - 1) / size_of::<Word>() / 8
    }

    // reserve a buffer for mantissa.
    fn reserve_new(sz: usize) -> Result<WordBuf, Error> {
        WordBuf::new(sz)
    }

    /// A [`DEFAULT_P`](crate::defs::DEFAULT_P)-bit mantissa whose top word is `top`, which must have its
    /// most significant bit set. Const and allocation-free.
    pub(crate) const fn from_normalized_top_word(top: Word) -> Self {
        Mantissa {
            m: WordBuf::inline_with_top_word(top),
            n: WORD_BIT_SIZE * INLINE_WORDS,
        }
    }

    /// The zero mantissa of [`DEFAULT_P`](crate::defs::DEFAULT_P) bits. Const and allocation-free.
    pub(crate) const fn zero_default() -> Self {
        Mantissa {
            m: WordBuf::inline_with_top_word(0),
            n: 0,
        }
    }

    /// New mantissa with length of at least `p` bits filled with zeroes.
    pub fn new(p: usize) -> Result<Self, Error> {
        let mut m = Self::reserve_new(Self::bit_len_to_word_len(p))?;

        m.fill(0);

        Ok(Mantissa { m, n: 0 })
    }

    /// New mantissa with length of at least `p` bits filled with 1.
    pub fn oned_mantissa(p: usize) -> Result<Self, Error> {
        let mut m = Self::reserve_new(Self::bit_len_to_word_len(p))?;

        let n = WORD_BIT_SIZE * m.len();

        m.fill(WORD_MAX);

        Ok(Mantissa { m, n })
    }

    /// New mantissa with length of at least `p` bits for the value of minimum positive value.
    pub fn min(p: usize) -> Result<Self, Error> {
        let mut m = Self::reserve_new(Self::bit_len_to_word_len(p))?;

        m.fill(0);

        m[0] = 1;

        Ok(Mantissa { m, n: 1 })
    }

    /// New mantissa with length of at least `p` bits for the value of 1.
    pub fn from_word(p: usize, mut d: Word) -> Result<Self, Error> {
        let mut m = Self::reserve_new(Self::bit_len_to_word_len(p))?;

        if m.is_empty() {
            return Err(Error::InvalidArgument);
        }

        m.fill(0);

        let l = m.len();

        if d > 0 {
            while d & WORD_SIGNIFICANT_BIT == 0 {
                d <<= 1;
            }
        }

        m[l - 1] = d;

        let n = WORD_BIT_SIZE * m.len();

        Ok(Mantissa { m, n })
    }

    /// New mantissa with length of at least `p` bits prefilled with `words`.
    pub fn from_words(p: usize, w: &[Word]) -> Result<Self, Error> {
        let mut m = Self::reserve_new(Self::bit_len_to_word_len(p))?;

        if m.len() < w.len() {
            return Err(Error::InvalidArgument);
        }

        let b = m.len() - w.len();
        (&mut m)[..b].fill(0);
        (&mut m)[b..].copy_from_slice(w);

        let n = Self::find_bit_len(&m);

        Ok(Mantissa { m, n })
    }

    /// Return true if mantissa represents zero.
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.n == 0
    }

    /// Fill with zeroes, ad set bit length to 0.
    #[inline]
    pub fn set_zero(&mut self) {
        self.m.fill(0);
        self.n = 0;
    }

    #[inline]
    pub fn set_bit_len(&mut self, n: usize) {
        self.n = n;
    }

    /// Recalculate bit length.
    pub fn update_bit_len(&mut self) {
        let n = Self::find_bit_len(&self.m);
        self.n = n;
    }

    /// Shift right by n bits.
    pub fn shift_right(&mut self, n: usize) {
        shift_slice_right(&mut self.m, n)
    }

    /// Shift self left by n digits.
    #[inline]
    pub fn shift_left(&mut self, n: usize) {
        shift_slice_left(&mut self.m, n)
    }

    /// Shift to the left, returns exponent shift as positive value.
    pub(crate) fn maximize(m: &mut [Word]) -> usize {
        let mut shift = 0;
        let mut d = 0;

        for v in m.iter().rev() {
            d = *v;

            if d != 0 {
                break;
            }

            shift += WORD_BIT_SIZE;
        }

        if d != 0 {
            while WORD_SIGNIFICANT_BIT & d == 0 {
                d <<= 1;
                shift += 1;
            }

            shift_slice_left(m, shift);
            shift
        } else {
            0
        }
    }

    /// Find the position of the first occurence of "1" starting from start_pos.
    pub fn find_one_from(&self, start_pos: usize) -> Option<usize> {
        find_one_from(&self.m, start_pos)
    }

    /// Compare to m2 and return positive if self > m2, negative if self < m2, 0 otherwise.
    pub fn abs_cmp(&self, m2: &Self, normalized: bool) -> SignedWord {
        let n = self.len().min(m2.len());

        if normalized {
            Self::cmp_iters(self.m.iter().rev(), m2.m.iter().rev(), n)
        } else {
            Self::cmp_iters(
                NormalizedView::new(self.m.iter().rev()),
                NormalizedView::new(m2.m.iter().rev()),
                n,
            )
        }
    }

    // The words are compared where they are stored (or as computed by NormalizedView): nothing is collected.
    fn cmp_iters<A, T>(mut iter1: T, mut iter2: T, n: usize) -> SignedWord
    where
        A: core::borrow::Borrow<Word>,
        T: Iterator<Item = A>,
    {
        for (a, b) in core::iter::zip(iter1.by_ref(), iter2.by_ref()).take(n) {
            let diff = *a.borrow() as SignedWord - *b.borrow() as SignedWord;
            if diff != 0 {
                return diff;
            }
        }

        for v in iter1 {
            if *v.borrow() != 0 {
                return 1;
            }
        }

        for v in iter2 {
            if *v.borrow() != 0 {
                return -1;
            }
        }

        0
    }

    /// Subtracts m2 from self. m2 is supposed to be shifted right by m2_shift bits.
    /// `inexact` is set to true if the result is not exact.
    pub fn abs_sub(
        &self,
        m2: &Self,
        spec: AddSubSpec,
        inexact: &mut bool,
    ) -> Result<(i64, Self), Error> {
        let AddSubSpec {
            p,
            m2_shift,
            rm,
            is_positive,
            full_prec,
        } = spec;
        // The operands are as stored (a subnormal operand has leading zero bits): the difference is exact
        // and `maximize` below renormalizes it before it is rounded to `p`.

        let p = Self::bit_len_to_word_len(p);
        let mut c = 0;

        let l2 = (m2_shift + WORD_BIT_SIZE - 2) / WORD_BIT_SIZE;
        let mut m3 = if l2 > p && l2 > self.len() && !full_prec {
            // m2 lays outside of m1's mantissa range and outside the desired precision
            //
            // m1 1XXXX00000000   - m1 and trailing zeroes
            // m2 000000001XXXX   - m2_shift, m2

            *inexact |= true;

            let m3l = self.len().max(p) + 1;
            let mut m3 = Mantissa::new(m3l * WORD_BIT_SIZE)?;

            // subtract 1
            let m1iter = ExtendedSlice::new(self.m.iter(), m3l - self.len(), &0);
            let m3iter = m3.m.iter_mut();

            c = 1;
            for (d, &a) in izip!(m3iter, m1iter) {
                if c > 0 {
                    if a == 0 {
                        *d = WORD_MAX;
                    } else {
                        *d = a - c;
                        c = 0;
                    }
                } else {
                    *d = a;
                }
            }

            m3
        } else {
            // m2_shift is smaller than the precision of m1
            // or m2_shift is smaller than the desired precision
            //
            // m1 1XXXXXX00000   - m1 and trailing zeroes
            // m2 0000001XXXXX   - m2_shift, m2

            let l = self.len().max(m2_shift.div_ceil(WORD_BIT_SIZE) + m2.len());

            let mut m3 = Mantissa::new(l * WORD_BIT_SIZE)?;

            let m1iter = ExtendedSlice::new(self.m.iter(), l - self.len(), &0);
            let m2iter = RightShiftedSlice::new(&m2.m, m2_shift, 0, l - m2.len());
            let m3iter = m3.m.iter_mut();

            for (d, a, b) in izip!(m3iter, m1iter, m2iter) {
                c = sub_borrow(*a, b, c, d);
            }

            m3
        };

        if !(c == 0) {
            return Err(Error::InvalidArgument);
        }

        let mut shift = Self::maximize(&mut m3.m) as i64;

        if full_prec {
            m3.m.trunc_trailing_zeroes();
        } else if m3.len() < p {
            let n = m3.len();
            m3.m.try_extend(p * WORD_BIT_SIZE)?;
            m3.m[..p - n].fill(0);
        } else if m3.len() > p {
            if m3.round_mantissa(
                (m3.len() - p) * WORD_BIT_SIZE,
                rm,
                is_positive,
                &mut false,
                m3.max_bit_len(),
                inexact,
            ) {
                shift -= 1;
            }
            m3.m.trunc_to(p * WORD_BIT_SIZE);
        }

        m3.n = m3.max_bit_len();

        Ok((shift, m3))
    }

    /// Returns exponent shift, and self + m2.
    pub fn abs_add(
        &self,
        m2: &Self,
        spec: AddSubSpec,
        inexact: &mut bool,
    ) -> Result<(i64, Self), Error> {
        let AddSubSpec {
            p,
            m2_shift,
            rm,
            is_positive,
            full_prec,
        } = spec;
        // The operands are as stored: a subnormal operand has leading zero bits and is not normalized
        // first. The sum is exact, and `maximize` below renormalizes it before it is rounded to `p`.
        // A zero `m2` (with `m2_shift` > 0) sums `self` with nothing: the result is `self` written out and
        // rounded by the code below, which is how `add` and `sub` handle a zero operand.
        if !(!m2.is_zero() || m2_shift > 0) {
            return Err(Error::InvalidArgument);
        }

        let p = Self::bit_len_to_word_len(p);
        let mut c = 0;

        let (mut m3, c) = if m2_shift == 0 {
            let l = self.len().max(m2.len()) + 1;

            let mut m3 = Mantissa::new(l * WORD_BIT_SIZE)?;

            let m1iter = RightShiftedSlice::new(&self.m, 1, 0, l - self.len());
            let m2iter = RightShiftedSlice::new(&m2.m, 1, 0, l - m2.len());
            let m3iter = m3.m.iter_mut();

            for (d, a, b) in izip!(m3iter, m1iter, m2iter) {
                c = add_carry(a, b, c, d);
            }

            if !(c == 0) {
                return Err(Error::InvalidArgument);
            }
            c = 1;

            (m3, c)
        } else {
            let l2 = (m2_shift + WORD_BIT_SIZE - 2) / WORD_BIT_SIZE;
            if l2 > p && l2 > self.len() && !full_prec {
                // m2 lays outside of m1's mantissa range and outside the desired precision
                //
                // m1 1XXXX00000000   - m1 and trailing zeroes
                // m2 000000001XXXX   - m2_shift, m2

                *inexact |= true;

                let mut m3 = Mantissa::new((p + 1) * WORD_BIT_SIZE)?;

                if m3.len() >= self.m.len() {
                    let l = m3.len() - self.len();
                    m3.m[l..].copy_from_slice(&self.m);
                } else {
                    let l = self.len() - m3.len();
                    m3.m.copy_from_slice(&self.m[l..]);
                }

                m3.m[0] |= 1; // add sticky for rounding

                (m3, 0)
            } else {
                // m2_shift is smaller than the precision of m1
                // or m2_shift is smaller than the desired precision
                //
                // m1 1XXXXXX00000   - m1 and trailing zeroes
                // m2 0000001XXXXX   - m2_shift, m2

                let l = self.len().max(m2_shift.div_ceil(WORD_BIT_SIZE) + m2.len()) + 1;

                let mut m3 = Mantissa::new(l * WORD_BIT_SIZE)?;

                let m1iter = ExtendedSlice::new(self.m.iter(), l - self.len(), &0);
                let m2iter = RightShiftedSlice::new(&m2.m, m2_shift, 0, l - m2.len());
                let m3iter = m3.m.iter_mut();

                for (d, a, b) in izip!(m3iter, m1iter, m2iter) {
                    c = add_carry(*a, b, c, d);
                }

                if c > 0 {
                    m3.shift_right(1);
                    let l = m3.m.len();
                    m3.m[l - 1] |= WORD_SIGNIFICANT_BIT;
                }

                (m3, c)
            }
        };

        let mut shift = -(c as i64);

        // Only a subnormal first operand (with a subnormal second one) leaves the top bit clear.
        if m3.m[m3.len() - 1] & WORD_SIGNIFICANT_BIT == 0 {
            shift += Self::maximize(&mut m3.m) as i64;
        }

        if full_prec {
            m3.m.trunc_trailing_zeroes();
        } else if m3.len() < p {
            let n = m3.len();
            m3.m.try_extend(p * WORD_BIT_SIZE)?;
            m3.m[..p - n].fill(0);
        } else if m3.len() > p {
            if m3.round_mantissa(
                (m3.len() - p) * WORD_BIT_SIZE,
                rm,
                is_positive,
                &mut false,
                m3.max_bit_len(),
                inexact,
            ) {
                shift -= 1;
            }
            m3.m.trunc_to(p * WORD_BIT_SIZE);
        }

        m3.n = m3.max_bit_len();

        Ok((shift, m3))
    }

    /// Multiply two mantissas, return result, exponent shift, and inexact flag.
    pub fn mul(
        &self,
        m2: &Self,
        p: usize,
        rm: RoundingMode,
        is_positive: bool,
        full_prec: bool,
        inexact: &mut bool,
    ) -> Result<(i64, Self), Error> {
        // The operands need not be normalized (a subnormal operand has leading zero bits): the product is
        // exact and `maximize` below renormalizes it, so the returned shift is the product's leading zero
        // count (one less when rounding carried), and the caller's exponent is `e1 + e2 - shift` with each
        // operand's own exponent.
        let mut m3 = Self::reserve_new(self.len() + m2.len())?;

        Self::mul_unbalanced(&self.m, &m2.m, &mut m3)?;

        let mut shift = Self::maximize(&mut m3) as i64;

        let mut m3 = Mantissa { m: m3, n: 0 };

        let p = Self::bit_len_to_word_len(p);

        if full_prec {
            m3.m.trunc_trailing_zeroes();
        } else if m3.len() < p {
            let n = m3.len();
            m3.m.try_extend(p * WORD_BIT_SIZE)?;
            m3.m[..p - n].fill(0);
        } else if m3.len() > p {
            if m3.round_mantissa(
                (m3.len() - p) * WORD_BIT_SIZE,
                rm,
                is_positive,
                &mut false,
                m3.max_bit_len(),
                inexact,
            ) {
                shift -= 1;
            }
            m3.m.trunc_to(p * WORD_BIT_SIZE);
        }

        if !(shift >= -1) {
            // prevent exponent overflow
            return Err(Error::InvalidArgument);
        }

        m3.n = m3.max_bit_len();

        Ok((shift, m3))
    }

    /// Number of leading zero bits of the mantissa (zero for a normalized mantissa).
    pub(crate) fn leading_zero_bits(&self) -> usize {
        let mut n = 0;
        for &w in self.m.iter().rev() {
            if w != 0 {
                return n + w.leading_zeros() as usize;
            }
            n += WORD_BIT_SIZE;
        }
        n
    }

    /// The word of `src` at index `i`, read in place; 0 past either end.
    fn word_at(src: &[Word], i: Option<usize>) -> Word {
        let Some(&w) = i.and_then(|j| src.get(j)) else {
            return 0;
        };
        w
    }

    /// Writes `src << shift` (bits shifted out at the top are dropped, zeros come in at the bottom) into
    /// `dst`, which is as long as `src`: each destination word is computed from the source words, no
    /// slice is copied.
    fn write_shifted_left(dst: &mut [Word], src: &[Word], shift: usize) {
        let ws = shift / WORD_BIT_SIZE;
        let bs = shift % WORD_BIT_SIZE;
        for (i, d) in dst.iter_mut().enumerate() {
            let hi = Self::word_at(src, i.checked_sub(ws));
            let lo = Self::word_at(src, i.checked_sub(ws + 1));
            *d = if bs == 0 { hi } else { (hi << bs) | (lo >> (WORD_BIT_SIZE - bs)) };
        }
    }

    /// Divide mantissa by mantissa, return result and exponent ajustment.
    pub fn div(
        &self,
        m2: &Self,
        p: usize,
        rm: RoundingMode,
        is_positive: bool,
        inexact: &mut bool,
    ) -> Result<(i64, Self), Error> {
        let p = Self::bit_len_to_word_len(p);

        // A mantissa is a fraction, so the divisor's zero low words (a short value held at a wide
        // precision) do not change its value: divide by the words above them. Without this a 1-word
        // divisor at 65600 bits ran the dense 1025-word division (r2nd/zeta root cause, with `mul`).
        let z2 =
            m2.m.iter()
                .position(|w| w != &0)
                .unwrap_or(m2.len())
                .min(m2.len().saturating_sub(1));
        let m2w: &[Word] = &m2.m[z2..];

        let k = (p + m2w.len()).max(self.len()) + 1;

        // The operands are as stored. A subnormal dividend is shifted left (normalized) while it is written
        // into the padded working buffer `m1` below, which exists for every dividend; a subnormal divisor
        // is shifted into its own working buffer `d`, because the division kernel needs a normalized
        // divisor. The caller's exponents are the operands' own, and the shifts are accounted for here.
        let lz1 = self.leading_zero_bits();
        let lz2 = m2.leading_zero_bits();

        let mut e_shift =
            (m2w.len() as i64 - k as i64) * WORD_BIT_SIZE as i64 + lz2 as i64 - lz1 as i64;

        let mut m1 = Self::reserve_new(k)?;
        let l = k - self.len();
        m1[..l].fill(0);
        Self::write_shifted_left(&mut m1[l..], &self.m, lz1);

        let d = if lz2 == 0 {
            None
        } else {
            let mut d = Self::reserve_new(m2w.len())?;
            Self::write_shifted_left(&mut d, m2w, lz2);
            Some(d)
        };
        let divisor: &[Word] = match &d {
            Some(d) => d,
            None => m2w,
        };

        let (q, r) = Self::div_unbalanced(&m1, divisor)?;

        e_shift += q.len() as i64 * WORD_BIT_SIZE as i64;

        let mut m3 = Mantissa { m: q, n: 0 };

        let r_sticky = r.iter().any(|x| x != &0);

        // rounding
        if r_sticky {
            *inexact |= true;

            if rm as u32 & 0b1100000 != 0 {
                m3.m[0] |= 1;
            } else if (rm == RoundingMode::FromZero
                || (is_positive && rm == RoundingMode::Up)
                || (!is_positive && rm == RoundingMode::Down))
                && m3.add_ulp()
            {
                let m3l = m3.len() - 1;
                m3.m[m3l] = WORD_SIGNIFICANT_BIT;
                e_shift += 1;
            }
        }

        e_shift -= Self::maximize(&mut m3.m) as i64;

        if m3.round_mantissa(
            (m3.len() - p) * WORD_BIT_SIZE,
            rm,
            is_positive,
            &mut false,
            m3.max_bit_len(),
            inexact,
        ) {
            e_shift += 1;
        }

        m3.m.trunc_to(p * WORD_BIT_SIZE);
        m3.n = m3.max_bit_len();

        Ok((e_shift, m3))
    }

    // Multiply d1 by word d and put result to d3 with overflow.
    pub(super) fn mul_by_word(d1: &[Word], d: DoubleWord, d3: &mut [Word]) {
        let mut m: DoubleWord = 0;
        for (v1, v2) in d1.iter().zip(d3.iter_mut()) {
            m = *v1 as DoubleWord * d + m / WORD_BASE;
            *v2 = (m % WORD_BASE) as Word;
        }
        d3[d1.len()] = (m / WORD_BASE) as Word;
    }

    /// `self` to the power `i` mod `n` (`i ≥ 1`; `i = 0` is `InvalidArgument`), by left-to-right
    /// square-and-multiply over every bit of `i`. The scan starts at the top set bit found with
    /// `leading_zeros`: the old scan shifted before testing, so it dropped bit 63 and an exponent
    /// at or past 2^63 (the gap of `MAX rem minpos`) was reduced mod 2^63 — a silent wrong remainder.
    pub fn pow_mod(&self, i: usize, n: &Self) -> Result<Self, Error> {
        if i == 0 {
            return Err(Error::InvalidArgument);
        }
        let top = usize::BITS - 1 - i.leading_zeros();
        let mut x = self.clone()?;
        for b in (0..top).rev() {
            x = x.mul_mod(&x, n)?;
            if (i >> b) & 1 != 0 {
                x = x.mul_mod(self, n)?;
            }
        }
        Ok(x)
    }

    /// Multiply `self` by 2 ^ i
    pub fn pow2(&mut self, i: usize) -> Result<(), Error> {
        if self.n + i > self.max_bit_len() {
            self.m.try_extend_2(self.n + i)?;
        }

        self.shift_left(i);

        Ok(())
    }

    /// Add at least n bits of precision, data is not moved
    pub fn extend_subnormal(&mut self, n: usize) -> Result<(), Error> {
        self.m.try_extend_2(self.n + n)
    }

    // self * m1 mod n
    pub fn mul_mod(&self, m1: &Self, n: &Self) -> Result<Self, Error> {
        if !(n.m[n.len() - 1] & WORD_SIGNIFICANT_BIT != 0) {
            return Err(Error::InvalidArgument);
        }

        let mut m = Self::reserve_new(self.len() + m1.len())?;

        Self::mul_unbalanced(&self.m, &m1.m, &mut m)?;

        m.trunc_leading_zeroes();

        let (_q, r) = Self::div_unbalanced(&m, &n.m)?;

        let n = Self::find_bit_len(&r);

        Ok(Mantissa { m: r, n })
    }

    // Returns remainder of division of `self` by `n`.
    pub fn rem(&self, n: &Self) -> Result<Self, Error> {
        let (_q, r) = Self::div_unbalanced(&self.m, &n.m)?;
        let n = Self::find_bit_len(&r);
        Ok(Mantissa { m: r, n })
    }

    fn find_bit_len(m: &[Word]) -> usize {
        let mut n = 0;

        for &v in m.iter().rev() {
            if v != 0 {
                let mut v = v;
                while v & WORD_SIGNIFICANT_BIT == 0 {
                    n += 1;
                    v <<= 1;
                }
                break;
            }
            n += WORD_BIT_SIZE;
        }

        m.len() * WORD_BIT_SIZE - n
    }

    pub fn from_usize(u: usize) -> Result<(usize, Self), Error> {
        let mut m;

        if u == 0 {
            m = Self::reserve_new(1)?;
            m[0] = 0 as Word;

            let ret = Mantissa { m, n: 0 };

            return Ok((0, ret));
        }

        m = Self::reserve_new(1)?;
        m[0] = u as Word;

        let shift = Self::maximize(&mut m);
        let mut ret = Mantissa { m, n: 0 };

        ret.n = ret.max_bit_len();

        Ok((shift, ret))
    }

    /// Shift `self` to the left and return shift value.
    pub fn normilize2(&mut self) -> usize {
        self.n = self.max_bit_len();
        Self::maximize(&mut self.m)
    }

    /// Set n bits to 0 from the left/right.
    pub fn mask_bits(&mut self, mut n: usize, from_left: bool) {
        if from_left {
            for v in self.m.iter_mut().rev() {
                if n >= WORD_BIT_SIZE {
                    *v = 0;
                    n -= WORD_BIT_SIZE;
                } else {
                    let mask = WORD_MAX >> n;
                    *v &= mask;
                    break;
                }
            }
        } else {
            for v in self.m.iter_mut() {
                if n >= WORD_BIT_SIZE {
                    *v = 0;
                    n -= WORD_BIT_SIZE;
                } else {
                    let mask = WORD_MAX << n;
                    *v &= mask;
                    break;
                }
            }
        }
    }

    /// Decompose to raw parts.
    pub fn as_raw_parts(&self) -> (&[Word], usize) {
        (&self.m, self.n)
    }

    /// Construct from raw parts.
    pub fn from_raw_parts(m: &[Word], n: usize) -> Result<Self, Error> {
        let mut mm = Self::reserve_new(m.len())?;
        mm.copy_from_slice(m);
        Ok(Mantissa { m: mm, n })
    }

    /// Returns true if all digits are equal to 0.
    pub fn is_all_zero(&self) -> bool {
        self.m.iter().all(|v| v == &0)
    }

    /// Returns length of the mantissa in words.
    #[inline]
    pub fn len(&self) -> usize {
        self.m.len()
    }

    /// Returns maximum possible length of the mantissa in digits of the mantissa's base.
    #[inline]
    pub fn max_bit_len(&self) -> usize {
        self.len() * WORD_BIT_SIZE
    }

    /// The roundability verdict of [`Self::round_mantissa`] called with `check_roundable = true`, read
    /// without changing `self`: true when rounding off `n` positions is decided by `s` correct digits (then
    /// `round_mantissa` with or without the check leaves the same kept bits), false when it is not.
    pub(crate) fn is_roundable(&self, n: usize, rm: RoundingMode, s: usize) -> bool {
        let s = s - s % WORD_BIT_SIZE;
        if rm.skips_rounding() {
            return true;
        }
        if !(n > 0 && n < self.max_bit_len()) {
            return false;
        }
        let mut check = true;
        let mut cc = (n - (self.max_bit_len() - s)) / WORD_BIT_SIZE;
        if rm == RoundingMode::ToEven || rm == RoundingMode::ToOdd {
            let known_from = self.max_bit_len() - s;
            let round_bit_known = n > known_from;
            let n = n - 1;
            let i = n / WORD_BIT_SIZE;
            let t = n % WORD_BIT_SIZE;
            let num = (self.m[i] >> t) & 1;
            let mut uniform: Option<Word> = None;
            if t > 0 && round_bit_known {
                let rb = self.m[i] << (WORD_BIT_SIZE - t);
                if rb == 0 {
                    uniform = Some(0);
                } else if rb == WORD_MAX << (WORD_BIT_SIZE - t) {
                    uniform = Some(WORD_MAX);
                } else {
                    check = false;
                }
            }
            if check && round_bit_known {
                let known_word = known_from / WORD_BIT_SIZE;
                for &v in self.m[known_word.min(i)..i].iter().rev() {
                    match uniform {
                        None if v == 0 => uniform = Some(0),
                        None if v == WORD_MAX => uniform = Some(WORD_MAX),
                        Some(u) if v == u => {}
                        _ => {
                            check = false;
                            break;
                        }
                    }
                }
            }
            if check
                && round_bit_known
                && ((uniform == Some(0) && num == 0) || (uniform == Some(WORD_MAX) && num == 1))
            {
                check = false;
            }
        } else {
            let mut td = 0;
            let i = n / WORD_BIT_SIZE;
            let mut z = i;
            let t = n % WORD_BIT_SIZE;
            if t > 0 {
                let rb = self.m[i] << (WORD_BIT_SIZE - t);
                if rb != 0 {
                    let msk = rb & (WORD_MAX << (WORD_BIT_SIZE - t));
                    if msk == rb {
                        td = WORD_MAX;
                    } else if msk != 0 {
                        check = false;
                    }
                }
            } else if cc > 0 && z > 0 {
                z -= 1;
                if self.m[z] == WORD_MAX {
                    td = WORD_MAX;
                } else if self.m[z] != 0 {
                    check = false;
                }
                cc -= 1;
            }
            for &v in self.m[..z].iter().rev() {
                if check && cc > 0 {
                    if v != td {
                        check = false;
                    } else {
                        cc -= 1;
                    }
                }
            }
        }
        !check
    }

    /// Round `n` positions, return true if exponent is to be incremented.
    /// If `check_roundable` is true on input, the function verifies whether the mantissa is roundable, given it contains `s` correct digits.
    /// If `check_roundable` is set to false on return, in any case it means rounding was successful.
    /// If some bits set to 1 during rounding were set to 0, `inexact` will take value true.
    pub fn round_mantissa(
        &mut self,
        n: usize,
        rm: RoundingMode,
        is_positive: bool,
        check_roundable: &mut bool,
        s: usize,
        inexact: &mut bool,
    ) -> bool {
        // The roundability scan below works on whole words of known bits. Claiming fewer correct bits is always
        // sound (the check only gets stricter), so an unaligned `s` is rounded down to a word multiple, in every
        // build (Jeff 2026-11-24 4:40 PM CT item 5; was a debug_assert that `s` is aligned).
        let s = s - s % WORD_BIT_SIZE;

        // This function is so complex, because it combines rounding for all rounding modes
        // and checks for roundability at the same time.

        if rm.skips_rounding() {
            *check_roundable = false;

            // inexact?
            if !(*inexact)
                && n > 0
                && n < self.max_bit_len()
                && self.find_one_from(self.max_bit_len() - n).is_some()
            {
                *inexact |= true;
            }

            return false;
        }

        let self_len = self.m.len();
        if n > 0 && n < self.max_bit_len() {
            let mut i;
            let mut t;
            let mut c = false;
            let mut cc = (n - (self.max_bit_len() - s)) / WORD_BIT_SIZE; // num of roundable significant words

            if rm == RoundingMode::ToEven || rm == RoundingMode::ToOdd {
                // Nearest: `n` is now the index of the round bit.
                let known_from = self.max_bit_len() - s; // bits below this index are unknown
                let round_bit_known = n > known_from;
                let n = n - 1;
                let mut rem_zero = true;

                // analyze bits at n and at n+1
                // to decide if we need to add 1 or not.
                let np1 = n + 1;
                i = n / WORD_BIT_SIZE;
                t = n % WORD_BIT_SIZE;
                let ip1 = np1 / WORD_BIT_SIZE;
                let tp1 = np1 % WORD_BIT_SIZE;
                let num = (self.m[i] >> t) & 1;

                // Roundability: the known bits below the round bit (the `t` bits of word `i` and the whole
                // words `known_from / WORD_BIT_SIZE ..i`) must not all be zeros after a round bit of one, or
                // all be ones after a round bit of zero (checked below), or the unknown bits could move the
                // value across the half-way point. `uniform` is the pattern
                // seen so far (zeros or ones); a bit group that breaks the pattern decides the rounding.
                let mut uniform: Option<Word> = None;
                if t > 0 {
                    let rb = self.m[i] << (WORD_BIT_SIZE - t);
                    if rb != 0 {
                        rem_zero = false;
                    }
                    if round_bit_known {
                        if rb == 0 {
                            uniform = Some(0);
                        } else if rb == WORD_MAX << (WORD_BIT_SIZE - t) {
                            uniform = Some(WORD_MAX);
                        } else {
                            *check_roundable = false; // self is roundable
                        }
                    }
                }
                if *check_roundable && round_bit_known {
                    let known_word = known_from / WORD_BIT_SIZE;
                    for &v in self.m[known_word.min(i)..i].iter().rev() {
                        match uniform {
                            None if v == 0 => uniform = Some(0),
                            None if v == WORD_MAX => uniform = Some(WORD_MAX),
                            Some(u) if v == u => {}
                            _ => {
                                *check_roundable = false; // self is roundable
                                break;
                            }
                        }
                    }
                }

                // Nearest only has a half-way point to be close to: a round bit of one followed by known
                // zeros, or a round bit of zero followed by known ones. A round bit of one followed by ones
                // (just below a representable number) rounds up, and a round bit of zero followed by zeros
                // (just above one) rounds down, whatever the unknown bits are.
                if *check_roundable
                    && round_bit_known
                    && ((uniform == Some(0) && num == 0) || (uniform == Some(WORD_MAX) && num == 1))
                {
                    *check_roundable = false; // self is roundable
                }

                let nump1 = if ip1 < self_len { (self.m[ip1] >> tp1) & 1 } else { 0 };

                // anything before n'th word becomes 0
                for v in &mut self.m[..i] {
                    if *v != 0 {
                        rem_zero = false;
                    }
                    *v = 0;
                }

                if *check_roundable {
                    return false; // self is not roundable
                }

                // inexactness
                *inexact |= !rem_zero || num == 1;

                // rounding
                let eq1 = num == 1 && rem_zero;
                let gt1 = num == 1 && !rem_zero;

                match rm {
                    RoundingMode::ToEven => {
                        if gt1 || (eq1 && (nump1 & 1) != 0) {
                            // add 1
                            c = true;
                        }
                    }
                    // Nearest, ties to odd (not sticky round-to-odd): only the exact tie `eq1` looks at the
                    // neighbour's parity.
                    RoundingMode::ToOdd => {
                        if gt1 || (eq1 && (nump1 & 1) == 0) {
                            // add 1
                            c = true;
                        }
                    }
                    RoundingMode::None
                    | RoundingMode::Work
                    | RoundingMode::Up
                    | RoundingMode::Down
                    | RoundingMode::ToZero
                    | RoundingMode::FromZero => {}
                };

                if c {
                    // add 1 at (n+1)'th position
                    if ip1 > i {
                        self.m[i] = 0;
                    }
                    i = ip1;
                    t = tp1;
                } else {
                    t += 1;
                }
            } else {
                let mut rem_zero = true;
                let mut td = 0; // either all bits set, or zero (default)

                i = n / WORD_BIT_SIZE;
                let mut z = i;
                t = n % WORD_BIT_SIZE;
                if t > 0 {
                    // bits coming before the least significant bit
                    let rb = self.m[i] << (WORD_BIT_SIZE - t);
                    if rb != 0 {
                        rem_zero = false;

                        if *check_roundable {
                            // try to get information from preceding bits
                            let msk = rb & (WORD_MAX << (WORD_BIT_SIZE - t));
                            if msk == rb {
                                td = WORD_MAX;
                            } else if msk != 0 {
                                *check_roundable = false; // self is roundable
                            }
                        }
                    }
                } else if *check_roundable && cc > 0 && z > 0 {
                    // try to get information from preceding word
                    z -= 1;
                    if self.m[z] == WORD_MAX {
                        td = WORD_MAX;
                        rem_zero = false; // sticky is set
                    } else if self.m[z] != 0 {
                        *check_roundable = false; // self is roundable
                        rem_zero = false; // sticky is set
                    }
                    cc -= 1;
                }

                // anything before n'th word becomes 0
                for v in self.m[..z].iter_mut().rev() {
                    if *check_roundable && cc > 0 {
                        if *v != td {
                            *check_roundable = false; // self is roundable
                        } else {
                            cc -= 1;
                        }
                    }
                    if *v != 0 {
                        rem_zero = false; // sticky is set
                    }
                    *v = 0;
                }

                if *check_roundable {
                    return false; // self is not roundable
                }

                // inexactness
                *inexact |= !rem_zero;

                // rounding
                match rm {
                    RoundingMode::ToZero => {}
                    RoundingMode::FromZero => {
                        if !rem_zero {
                            // add 1
                            c = true;
                        }
                    }
                    RoundingMode::Up => {
                        if !rem_zero && is_positive {
                            // add 1
                            c = true;
                        }
                    }
                    RoundingMode::Down => {
                        if !rem_zero && !is_positive {
                            // add 1
                            c = true;
                        }
                    }
                    RoundingMode::None
                    | RoundingMode::Work
                    | RoundingMode::ToEven
                    | RoundingMode::ToOdd => {}
                };
            }

            // `*check_roundable` is false here: both branches above end in `if *check_roundable { return false; }`
            // (Jeff 2026-11-24 4:40 PM CT item 5; was a debug_assert).

            if c {
                if i < self_len {
                    if (self.m[i] >> t) as DoubleWord + 1 < (WORD_BASE >> t) {
                        self.m[i] = ((self.m[i] >> t) + 1) << t;
                        return false;
                    } else {
                        self.m[i] = 0;
                    }
                }

                // process overflows
                i += 1;
                for v in &mut self.m[i..] {
                    if *v < WORD_MAX {
                        *v += 1;
                        return false;
                    } else {
                        *v = 0;
                    }
                }
                self.m[self_len - 1] = WORD_SIGNIFICANT_BIT;
                return true;
            } else {
                // just remove trailing bits
                self.m[i] = if t >= WORD_BIT_SIZE { 0 } else { (self.m[i] >> t) << t };
            }
        }
        false
    }

    /// Sets the precision to `p`.
    pub fn set_length(&mut self, p: usize) -> Result<(), Error> {
        let sz = Self::bit_len_to_word_len(p);
        let orig_len = self.m.len();
        if sz < orig_len {
            self.m.trunc_to(p);
            let nn = (orig_len - sz) * WORD_BIT_SIZE;
            if self.n >= nn {
                self.n -= nn;
            } else {
                self.n = 0;
            }
        } else if sz > orig_len {
            self.m.try_extend(p)?;
            if self.n != 0 {
                self.n += (sz - orig_len) * WORD_BIT_SIZE;
            }
        }
        Ok(())
    }

    /// Shift to the left by n bits and increase or decrease precision if needed just to fit the data.
    #[cfg(test)]
    pub fn shift_left_resize(&mut self, n: usize) -> Result<(), Error> {
        if self.max_bit_len() >= self.bit_len() + n {
            self.shift_left(n);
            self.m.trunc_leading_zeroes();
        } else {
            self.m.try_extend_3(self.bit_len() + n, n)?;
        }

        self.n += n;

        Ok(())
    }

    pub fn most_significant_word(&self) -> Word {
        if self.n > 0 {
            self.m[(self.n - 1) / WORD_BIT_SIZE]
        } else {
            0
        }
    }

    #[cfg(any(test, feature = "random"))]
    /// Returns randomized mantissa with at least p bits of length.
    pub fn random_normal(p: usize) -> Result<Self, Error> {
        let mut m = Self::reserve_new(Self::bit_len_to_word_len(p))?;
        for v in m.iter_mut() {
            *v = crate::common::test_rng::random::<Word>();
        }
        let mut ret = Mantissa { m, n: 0 };
        if !ret.is_all_zero() {
            Self::maximize(&mut ret.m);
            ret.n = WORD_BIT_SIZE * ret.m.len();
            ret.m[0] ^= crate::common::test_rng::random::<Word>() >> 1;
        }
        Ok(ret)
    }

    /// Clones the mantissa.
    pub fn clone(&self) -> Result<Self, Error> {
        #[cfg(test)]
        CLONE_CALLS.with(|c| c.set(c.get() + 1));
        let mut m = Self::reserve_new(self.m.len())?;
        m.copy_from_slice(&self.m);
        Ok(Mantissa { m, n: self.n })
    }

    /// The mantissa of a new result, read straight from `self`: `self` shifted left by `shl` bits (at most
    /// its leading zero bits), stored in `dst_bits` rounded up to whole words, and rounded with `rm` while
    /// the words are written. Gives the same words and bit length as `shift_left(shl)` (with the bit length
    /// raised by `shl`), `round_mantissa(.., check_roundable = false, s = max_bit_len)` and
    /// `set_length(dst_bits)` on an owned mantissa; `self` is only read.
    /// Returns the mantissa and whether rounding carried out of the top word (the exponent then grows by
    /// one). `inexact` is set when a dropped bit was one.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for the result.
    pub(crate) fn rounded_from(
        &self,
        shl: usize,
        dst_bits: usize,
        rm: RoundingMode,
        is_positive: bool,
        inexact: &mut bool,
    ) -> Result<(Self, bool), Error> {
        let src: &[Word] = &self.m;
        let (ws, bs) = (shl / WORD_BIT_SIZE, shl % WORD_BIT_SIZE);
        // Word `i` of `self << shl`.
        let word = |i: usize| -> Word {
            if i < ws {
                return 0;
            }
            let j = i - ws;
            let hi = src[j] << bs;
            if bs > 0 && j > 0 {
                hi | (src[j - 1] >> (WORD_BIT_SIZE - bs))
            } else {
                hi
            }
        };
        let (sl, dl) = (src.len(), Self::bit_len_to_word_len(dst_bits));
        let mut m = Self::reserve_new(dl)?;
        let mut n = self.n + shl;
        if dl >= sl {
            let low = dl - sl;
            for (i, d) in m[low..].iter_mut().enumerate() {
                *d = word(i);
            }
            if n != 0 {
                n += low * WORD_BIT_SIZE;
            }
            return Ok((Mantissa { m, n }, false));
        }
        let k = sl - dl;
        for (i, d) in m.iter_mut().enumerate() {
            *d = word(k + i);
        }
        let dropped = |to: usize| (0..to).any(|i| word(i) != 0);
        let up = match rm {
            RoundingMode::None | RoundingMode::Work => {
                *inexact |= dropped(k);
                false
            }
            RoundingMode::ToEven | RoundingMode::ToOdd => {
                let r = word(k - 1);
                let half = r >> (WORD_BIT_SIZE - 1) == 1;
                let rest = (r << 1) != 0 || dropped(k - 1);
                *inexact |= half || rest;
                let odd = m[0] & 1 == 1;
                half && (rest || odd == (rm == RoundingMode::ToEven))
            }
            RoundingMode::ToZero
            | RoundingMode::FromZero
            | RoundingMode::Up
            | RoundingMode::Down => {
                let nz = dropped(k);
                *inexact |= nz;
                nz && match rm {
                    RoundingMode::FromZero => true,
                    RoundingMode::Up => is_positive,
                    RoundingMode::Down => !is_positive,
                    _ => false,
                }
            }
        };
        let mut ovf = false;
        if up {
            ovf = true;
            for d in m.iter_mut() {
                if *d < WORD_MAX {
                    *d += 1;
                    ovf = false;
                    break;
                }
                *d = 0;
            }
            if ovf {
                m[dl - 1] = WORD_SIGNIFICANT_BIT;
            }
        }
        let nn = k * WORD_BIT_SIZE;
        n = n.saturating_sub(nn);
        Ok((Mantissa { m, n }, ovf))
    }

    /// The mantissa `set_precision(p, rm)` leaves, read straight from `self`, for `0 < p < max_bit_len`:
    /// the words of `round_mantissa(max_bit_len - p, rm, is_positive, &mut false, max_bit_len, inexact)`
    /// followed by `set_length(p)` on an owned mantissa. Only the `ceil(p / W)` kept words are written; the
    /// dropped words are only read for the round and sticky bits. The bit length is `self`'s, less the
    /// dropped words (as `set_length`). Returns the mantissa, whether the rounding carried out of the top
    /// word, and whether the rounded mantissa at full length would be all zeros.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for the result.
    pub(crate) fn precision_from(
        &self,
        p: usize,
        rm: RoundingMode,
        is_positive: bool,
        inexact: &mut bool,
    ) -> Result<(Self, bool, bool), Error> {
        let src: &[Word] = &self.m;
        let len = self.max_bit_len();
        if !(0 < p && p < len) {
            return Err(Error::InvalidArgument);
        }
        let drop = len - p;
        let (k, t) = (drop / WORD_BIT_SIZE, drop % WORD_BIT_SIZE);
        // Bit `i` (from the bottom) of `self`, and whether every bit below `i` is zero.
        let bit = |i: usize| (src[i / WORD_BIT_SIZE] >> (i % WORD_BIT_SIZE)) & 1 == 1;
        let zero_below = |i: usize| {
            let (w, b) = (i / WORD_BIT_SIZE, i % WORD_BIT_SIZE);
            src[..w].iter().all(|v| v == &0) && (b == 0 || src[w] << (WORD_BIT_SIZE - b) == 0)
        };
        let up = match rm {
            RoundingMode::None | RoundingMode::Work => {
                // `round_mantissa` only records inexactness: the low `t` bits of word `k` stay.
                *inexact |= !zero_below(drop);
                None
            }
            RoundingMode::ToEven | RoundingMode::ToOdd => {
                let half = bit(drop - 1);
                let rest_zero = zero_below(drop - 1);
                *inexact |= !rest_zero || half;
                let odd = bit(drop);
                let tie_up = if rm == RoundingMode::ToEven { odd } else { !odd };
                Some(half && (!rest_zero || tie_up))
            }
            RoundingMode::ToZero
            | RoundingMode::FromZero
            | RoundingMode::Up
            | RoundingMode::Down => {
                let nz = !zero_below(drop);
                *inexact |= nz;
                Some(
                    nz && match rm {
                        RoundingMode::FromZero => true,
                        RoundingMode::Up => is_positive,
                        RoundingMode::Down => !is_positive,
                        _ => false,
                    },
                )
            }
        };
        let mut m = Self::reserve_new(src.len() - k)?;
        for (d, w) in m.iter_mut().zip(src[k..].iter()) {
            *d = *w;
        }
        let mut ovf = false;
        match up {
            None => {}
            Some(false) => m[0] = (m[0] >> t) << t,
            Some(true) => {
                if ((m[0] >> t) as DoubleWord) + 1 < (WORD_BASE >> t) {
                    m[0] = ((m[0] >> t) + 1) << t;
                } else {
                    m[0] = 0;
                    ovf = true;
                    for v in m[1..].iter_mut() {
                        if *v < WORD_MAX {
                            *v += 1;
                            ovf = false;
                            break;
                        }
                        *v = 0;
                    }
                    if ovf {
                        let top = m.len() - 1;
                        m[top] = WORD_SIGNIFICANT_BIT;
                    }
                }
            }
        }
        let full_zero =
            m.iter().all(|v| v == &0) && (up.is_some() || zero_below(k * WORD_BIT_SIZE));
        let n = self.n.saturating_sub(k * WORD_BIT_SIZE);
        Ok((Mantissa { m, n }, ovf, full_zero))
    }

    /// A new mantissa holding the top `keep` bits of `self` and zeros below them, written straight from
    /// `self`: the words and bit length of `mask_bits(max_bit_len - keep, false)` on an owned mantissa (a
    /// `keep` of zero gives a zero, as `set_zero`).
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for the result.
    pub(crate) fn top_bits_from(&self, keep: usize) -> Result<Self, Error> {
        let mut m = Self::reserve_new(self.m.len())?;
        if keep == 0 {
            return Ok(Mantissa { m, n: 0 });
        }
        let mut drop = self.max_bit_len().saturating_sub(keep);
        for (d, &v) in m.iter_mut().zip(self.m.iter()) {
            if drop >= WORD_BIT_SIZE {
                drop -= WORD_BIT_SIZE;
            } else {
                *d = v & (WORD_MAX << drop);
                drop = 0;
            }
        }
        Ok(Mantissa { m, n: self.n })
    }

    pub fn digits(&self) -> &[Word] {
        &self.m
    }

    /// True when the digit buffer is stored without a heap allocation.
    #[inline]
    pub fn is_inline(&self) -> bool {
        self.m.is_inline()
    }

    pub fn digits_mut(&mut self) -> &mut [Word] {
        &mut self.m
    }

    pub fn bit_len(&self) -> usize {
        self.n
    }

    /// returns true if `self` represents odd integer.
    pub fn is_odd_int(&self, e: usize) -> bool {
        // e == 0: the value is 0.m in [1/2, 1), not an integer (Jeff 2026-11-24 4:40 PM CT item 5; was a debug_assert that e > 0).
        if e == 0 {
            return false;
        }
        let n = self.max_bit_len() - e;

        if n == 0 {
            return self.m[0] & 1 > 0;
        }

        if ((self.m[n / WORD_BIT_SIZE]) >> (n % WORD_BIT_SIZE)) & 1 != 0 {
            self.find_one_from(e).is_none()
        } else {
            false
        }
    }

    // Add 1 ulp to mantissa, return true if carry occures
    pub(crate) fn add_ulp(&mut self) -> bool {
        let mut c = 1;
        for v in self.m.iter_mut() {
            c = add_carry(*v, 0, c, v);
            if c == 0 {
                break;
            }
        }
        c > 0
    }

    /// Compute the square root.
    pub fn sqrt(
        &self,
        p: usize,
        rm: RoundingMode,
        is_positive: bool,
        inexact: &mut bool,
        odd_exp: bool,
    ) -> Result<(i64, Self), Error> {
        let p = Self::bit_len_to_word_len(p);

        let k = (p.max(self.len()) + 1) * 2;

        let mut e_shift = 0;

        let mut m1 = Self::reserve_new(k)?;

        let l = k - self.len();
        m1[l..].copy_from_slice(&self.m);
        m1[..l].fill(0);

        if odd_exp {
            shift_slice_right(&mut m1, 1);
        }

        let (q, r) = Self::sqrt_rem(&m1)?;

        let mut m3 = Mantissa { m: q, n: 0 };

        let r_sticky = r.iter().any(|x| x != &0);

        // rounding
        if r_sticky {
            *inexact |= true;

            if rm as u32 & 0b1100000 != 0 {
                m3.m[0] |= 1;
            } else if (rm == RoundingMode::FromZero
                || (is_positive && rm == RoundingMode::Up)
                || (!is_positive && rm == RoundingMode::Down))
                && m3.add_ulp()
            {
                let m3l = m3.len() - 1;
                m3.m[m3l] = WORD_SIGNIFICANT_BIT;
                e_shift += 1;
            }
        }

        _ = Self::maximize(&mut m3.m) as i64;

        if m3.round_mantissa(
            (m3.len() - p) * WORD_BIT_SIZE,
            rm,
            is_positive,
            &mut false,
            m3.max_bit_len(),
            inexact,
        ) {
            e_shift += 1;
        }

        m3.m.trunc_to(p * WORD_BIT_SIZE);
        m3.n = m3.max_bit_len();

        Ok((e_shift, m3))
    }

    /// Compute the cube root.
    pub fn cbrt(
        &self,
        p: usize,
        rm: RoundingMode,
        is_positive: bool,
        inexact: &mut bool,
        add_exp: i64,
    ) -> Result<(i64, Self), Error> {
        let p = Self::bit_len_to_word_len(p);

        let k = p.max(self.len()) + 1;
        let k3 = k * 3;

        let mut e_shift = -((k * WORD_BIT_SIZE) as i64);

        let mut m1 = Self::reserve_new(k3)?;

        let l = k3 - self.len();
        m1[l..].copy_from_slice(&self.m);
        m1[..l].fill(0);

        if !(add_exp >= 0) {
            return Err(Error::InvalidArgument);
        }

        shift_slice_right(&mut m1, add_exp as usize);

        let (q, r) = Self::cbrt_rem(m1)?;

        let mut m3 = Mantissa { m: q, n: 0 };

        let r_sticky = r.iter().any(|x| x != &0);

        // rounding
        if r_sticky {
            *inexact |= true;

            if rm as u32 & 0b1100000 != 0 {
                m3.m[0] |= 1;
            } else if (rm == RoundingMode::FromZero
                || (is_positive && rm == RoundingMode::Up)
                || (!is_positive && rm == RoundingMode::Down))
                && m3.add_ulp()
            {
                let m3l = m3.len() - 1;
                m3.m[m3l] = WORD_SIGNIFICANT_BIT;
                e_shift += 1;
            }
        }

        e_shift += (m3.len() * WORD_BIT_SIZE) as i64;

        e_shift -= Self::maximize(&mut m3.m) as i64;

        if m3.round_mantissa(
            (m3.len() - p) * WORD_BIT_SIZE,
            rm,
            is_positive,
            &mut false,
            m3.max_bit_len(),
            inexact,
        ) {
            e_shift += 1;
        }

        m3.m.trunc_to(p * WORD_BIT_SIZE);
        m3.n = m3.max_bit_len();

        Ok((e_shift, m3))
    }

    pub fn from_word_buf(m: WordBuf) -> Self {
        let n = Self::find_bit_len(&m);
        Mantissa { m, n }
    }
}

#[cfg(test)]
std::thread_local! {
    /// Number of `Mantissa::clone` calls made by this thread (unit tests count them to show that an
    /// operation does not copy its operands).
    pub(crate) static CLONE_CALLS: core::cell::Cell<usize> = const { core::cell::Cell::new(0) };
}
