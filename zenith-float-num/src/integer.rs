//! Arbitrary-precision signed integers on `Word` limbs.

use crate::common::util::add_carry;
use crate::common::util::sub_borrow;
use crate::defs::DoubleWord;
use crate::defs::Error;
use crate::defs::Exponent;
use crate::defs::RoundingMode;
use crate::defs::Sign;
use crate::defs::Word;
use crate::defs::DEFAULT_P;
use crate::defs::WORD_BIT_SIZE;
use crate::ExactNum;
use alloc::vec;
use alloc::vec::Vec;
use core::cmp::Ordering;
use core::ops::Neg;

/// Arbitrary-precision signed integer. Distinct from [`ExactNum`], which is floating-point.
///
/// Limbs are little-endian `Word`s with no leading zeros. Zero is an empty limb
/// vector and a positive sign.
#[derive(Clone, Debug, Eq)]
pub struct ExactInt {
    sign: Sign,
    limbs: Vec<Word>,
}

impl ExactInt {
    /// The integer `0`.
    pub fn zero() -> Self {
        Self {
            sign: Sign::Pos,
            limbs: Vec::new(),
        }
    }

    /// The integer `1`.
    pub fn one() -> Self {
        Self::from_u64(1)
    }

    fn normalize(&mut self) {
        while self.limbs.last().is_some_and(|w| *w == 0) {
            self.limbs.pop();
        }
        if self.limbs.is_empty() {
            self.sign = Sign::Pos;
        }
    }

    fn from_limbs(sign: Sign, limbs: Vec<Word>) -> Self {
        let mut s = Self { sign, limbs };
        s.normalize();
        s
    }

    /// Little-endian unsigned limbs; leading zeros are stripped.
    #[cfg(test)]
    pub(crate) fn from_unsigned_limbs(limbs: Vec<Word>) -> Self {
        Self::from_limbs(Sign::Pos, limbs)
    }

    /// Integer from `u64`.
    pub fn from_u64(n: u64) -> Self {
        Self::from_u128(n as u128)
    }

    /// Integer from `u128`.
    pub fn from_u128(mut n: u128) -> Self {
        if n == 0 {
            return Self::zero();
        }
        let mut limbs = Vec::new();
        while n > 0 {
            limbs.push(n as Word);
            n >>= WORD_BIT_SIZE;
        }
        Self {
            sign: Sign::Pos,
            limbs,
        }
    }

    /// Integer from `i64`.
    pub fn from_i64(n: i64) -> Self {
        Self::from_i128(n as i128)
    }

    /// Integer from `i128`.
    pub fn from_i128(n: i128) -> Self {
        if n == i128::MIN {
            let mut v = Self::from_u128(1u128 << 127);
            v.sign = Sign::Neg;
            return v;
        }
        let mut v = Self::from_u128(n.unsigned_abs());
        if n < 0 {
            v.sign = Sign::Neg;
        }
        v
    }

    /// True if the value is `0`.
    pub fn is_zero(&self) -> bool {
        self.limbs.is_empty()
    }

    /// `Some(n)` if `self` is a non-negative integer that fits in `u128`.
    pub(crate) fn to_u128(&self) -> Option<u128> {
        if self.is_negative() || self.bit_length() > 128 {
            return None;
        }
        let mut n = 0u128;
        for (i, &w) in self.limbs.iter().enumerate() {
            n |= (w as u128) << (i * WORD_BIT_SIZE);
        }
        Some(n)
    }

    /// True if the value is strictly negative.
    pub fn is_negative(&self) -> bool {
        self.sign == Sign::Neg && !self.is_zero()
    }

    /// `-1`, `0`, or `+1`.
    pub fn signum(&self) -> Self {
        if self.is_zero() {
            Self::zero()
        } else if self.sign == Sign::Neg {
            Self::from_i64(-1)
        } else {
            Self::one()
        }
    }

    /// Number of bits in `|self|`. Zero has bit length `0`.
    pub fn bit_length(&self) -> usize {
        match self.limbs.last() {
            None => 0,
            Some(last) => {
                (self.limbs.len() - 1) * WORD_BIT_SIZE
                    + (WORD_BIT_SIZE - last.leading_zeros() as usize)
            }
        }
    }

    fn cmp_abs(a: &[Word], b: &[Word]) -> Ordering {
        match a.len().cmp(&b.len()) {
            Ordering::Equal => {
                for (x, y) in a.iter().rev().zip(b.iter().rev()) {
                    match x.cmp(y) {
                        Ordering::Equal => {}
                        o => return o,
                    }
                }
                Ordering::Equal
            }
            o => o,
        }
    }

    /// Compare as signed integers.
    pub fn soft_cmp(&self, other: &Self) -> Ordering {
        if self.is_zero() && other.is_zero() {
            return Ordering::Equal;
        }
        match (self.sign, other.sign) {
            (Sign::Pos, Sign::Neg) => Ordering::Greater,
            (Sign::Neg, Sign::Pos) => Ordering::Less,
            (Sign::Pos, Sign::Pos) => Self::cmp_abs(&self.limbs, &other.limbs),
            (Sign::Neg, Sign::Neg) => Self::cmp_abs(&other.limbs, &self.limbs),
        }
    }

    fn add_abs(a: &[Word], b: &[Word]) -> Vec<Word> {
        let n = a.len().max(b.len());
        let mut out = vec![0; n + 1];
        let mut c = 0;
        for (i, out_i) in out.iter_mut().enumerate().take(n) {
            let x = a.get(i).copied().unwrap_or(0);
            let y = b.get(i).copied().unwrap_or(0);
            c = add_carry(x, y, c, out_i);
        }
        out[n] = c;
        out
    }

    fn sub_abs(a: &[Word], b: &[Word]) -> Vec<Word> {
        debug_assert!(Self::cmp_abs(a, b) != Ordering::Less);
        let mut out = vec![0; a.len()];
        let mut c = 0;
        for i in 0..a.len() {
            let y = b.get(i).copied().unwrap_or(0);
            c = sub_borrow(a[i], y, c, &mut out[i]);
        }
        debug_assert!(c == 0);
        out
    }

    /// `self + rhs`.
    pub fn add(&self, rhs: &Self) -> Self {
        if self.is_zero() {
            return rhs.clone();
        }
        if rhs.is_zero() {
            return self.clone();
        }
        if self.sign == rhs.sign {
            Self::from_limbs(self.sign, Self::add_abs(&self.limbs, &rhs.limbs))
        } else {
            match Self::cmp_abs(&self.limbs, &rhs.limbs) {
                Ordering::Equal => Self::zero(),
                Ordering::Greater => {
                    Self::from_limbs(self.sign, Self::sub_abs(&self.limbs, &rhs.limbs))
                }
                Ordering::Less => {
                    Self::from_limbs(rhs.sign, Self::sub_abs(&rhs.limbs, &self.limbs))
                }
            }
        }
    }

    /// In-place `self += rhs`. Same-sign adds reuse `self`'s limb buffer
    /// (the Bernoulli reverse walk sums hundreds of odd powers every step).
    pub fn add_assign(&mut self, rhs: &Self) {
        if rhs.is_zero() {
            return;
        }
        if self.is_zero() {
            *self = rhs.clone();
            return;
        }
        if self.sign != rhs.sign {
            *self = self.add(rhs);
            return;
        }
        let n = self.limbs.len().max(rhs.limbs.len());
        self.limbs.resize(n + 1, 0);
        let mut c = 0;
        for i in 0..n {
            let y = rhs.limbs.get(i).copied().unwrap_or(0);
            c = add_carry(self.limbs[i], y, c, &mut self.limbs[i]);
        }
        self.limbs[n] = c;
        self.normalize();
    }

    /// `self - rhs`.
    pub fn sub(&self, rhs: &Self) -> Self {
        self.add(&rhs.neg())
    }

    /// `-self`.
    pub fn neg(&self) -> Self {
        if self.is_zero() {
            return Self::zero();
        }
        Self {
            sign: self.sign.invert(),
            limbs: self.limbs.clone(),
        }
    }

    /// `self * rhs` for a `u64` factor (FLINT `fmpz_mul_ui`). The even-Bernoulli
    /// reverse walk multiplies each odd power by `j²` every step; those factors
    /// fit in a word for all `j` we store.
    pub fn mul_u64(&self, rhs: u64) -> Self {
        if self.is_zero() || rhs == 0 {
            return Self::zero();
        }
        if rhs == 1 {
            return self.clone();
        }
        let r = rhs as DoubleWord;
        let mut out = vec![0; self.limbs.len() + 1];
        let mut k: DoubleWord = 0;
        for (i, &d) in self.limbs.iter().enumerate() {
            let m = (d as DoubleWord) * r + k;
            out[i] = m as Word;
            k = m >> WORD_BIT_SIZE;
        }
        out[self.limbs.len()] = k as Word;
        Self::from_limbs(self.sign, out)
    }

    /// In-place `self *= rhs` for a `u64` factor.
    pub fn mul_u64_assign(&mut self, rhs: u64) {
        if self.is_zero() {
            return;
        }
        if rhs == 0 {
            *self = Self::zero();
            return;
        }
        if rhs == 1 {
            return;
        }
        let r = rhs as DoubleWord;
        let mut k: DoubleWord = 0;
        for d in self.limbs.iter_mut() {
            let m = (*d as DoubleWord) * r + k;
            *d = m as Word;
            k = m >> WORD_BIT_SIZE;
        }
        if k != 0 {
            self.limbs.push(k as Word);
        }
    }

    /// `self * rhs`.
    pub fn mul(&self, rhs: &Self) -> Self {
        if self.is_zero() || rhs.is_zero() {
            return Self::zero();
        }
        let mut out = vec![0; self.limbs.len() + rhs.limbs.len()];
        for (i, &d1) in self.limbs.iter().enumerate() {
            let d1 = d1 as DoubleWord;
            if d1 == 0 {
                continue;
            }
            let mut k = 0;
            for (j, &d2) in rhs.limbs.iter().enumerate() {
                let m = d1 * (d2 as DoubleWord) + out[i + j] as DoubleWord + k;
                out[i + j] = m as Word;
                k = m >> WORD_BIT_SIZE;
            }
            out[i + rhs.limbs.len()] = k as Word;
        }
        let sign = if self.sign == rhs.sign { Sign::Pos } else { Sign::Neg };
        Self::from_limbs(sign, out)
    }

    /// Product using the mantissa Toom/FFT multipliers once either factor
    /// exceeds 32 words. The schoolbook [`Self::mul`] path is kept for small
    /// operands and for callers that want an infallible multiply.
    pub(crate) fn mul_fast(&self, rhs: &Self) -> Result<Self, Error> {
        if self.is_zero() || rhs.is_zero() {
            return Ok(Self::zero());
        }
        if self.limbs.len() <= 32 && rhs.limbs.len() <= 32 {
            return Ok(self.mul(rhs));
        }
        let mut out = vec![0; self.limbs.len() + rhs.limbs.len()];
        crate::mantissa::Mantissa::mul_unbalanced(&self.limbs, &rhs.limbs, &mut out)?;
        let sign = if self.sign == rhs.sign { Sign::Pos } else { Sign::Neg };
        Ok(Self::from_limbs(sign, out))
    }

    fn shl_one(a: &mut Vec<Word>) {
        let mut c = 0;
        for w in a.iter_mut() {
            let n = (*w << 1) | c;
            c = *w >> (WORD_BIT_SIZE - 1);
            *w = n;
        }
        if c != 0 {
            a.push(c);
        }
    }

    fn bit(a: &[Word], i: usize) -> bool {
        let w = i / WORD_BIT_SIZE;
        let b = i % WORD_BIT_SIZE;
        a.get(w).is_some_and(|v| (*v >> b) & 1 == 1)
    }

    fn set_bit(a: &mut [Word], i: usize) {
        let w = i / WORD_BIT_SIZE;
        let b = i % WORD_BIT_SIZE;
        if w < a.len() {
            a[w] |= 1 << b;
        }
    }

    fn div_rem_abs(a: &[Word], b: &[Word]) -> (Vec<Word>, Vec<Word>) {
        debug_assert!(!b.is_empty());
        if Self::cmp_abs(a, b) == Ordering::Less {
            return (Vec::new(), a.to_vec());
        }
        let Some(top) = a.last() else {
            // a and b both empty: 0 / 0 has no quotient bits; the zero remainder is empty.
            return (Vec::new(), Vec::new());
        };
        let bits = (a.len() - 1) * WORD_BIT_SIZE + (WORD_BIT_SIZE - top.leading_zeros() as usize);
        let mut rem: Vec<Word> = Vec::new();
        let mut quot = vec![0; a.len()];
        for i in (0..bits).rev() {
            Self::shl_one(&mut rem);
            if rem.is_empty() {
                rem.push(0);
            }
            if Self::bit(a, i) {
                rem[0] |= 1;
            }
            while rem.last().is_some_and(|w| *w == 0) {
                rem.pop();
            }
            if Self::cmp_abs(&rem, b) != Ordering::Less {
                rem = Self::sub_abs(&rem, b);
                while rem.last().is_some_and(|w| *w == 0) {
                    rem.pop();
                }
                Self::set_bit(&mut quot, i);
            }
        }
        (quot, rem)
    }

    /// Truncated division: `(quotient, remainder)` with remainder sign matching `self`.
    /// Zero divisor is `None`.
    pub fn div_rem(&self, rhs: &Self) -> Option<(Self, Self)> {
        if rhs.is_zero() {
            return None;
        }
        if self.is_zero() {
            return Some((Self::zero(), Self::zero()));
        }
        let (q, r) = Self::div_rem_abs(&self.limbs, &rhs.limbs);
        let qsign = if self.sign == rhs.sign { Sign::Pos } else { Sign::Neg };
        Some((Self::from_limbs(qsign, q), Self::from_limbs(self.sign, r)))
    }

    /// Non-negative GCD. `gcd(0, 0) = 0`.
    pub fn gcd(&self, rhs: &Self) -> Self {
        let mut a = Self::from_limbs(Sign::Pos, self.limbs.clone());
        let mut b = Self::from_limbs(Sign::Pos, rhs.limbs.clone());
        while !b.is_zero() {
            let r = a.div_rem(&b).map(|(_, r)| r).unwrap_or_else(Self::zero);
            a = b;
            b = r;
        }
        a
    }

    /// `self.pow(exp)` by binary exponentiation. `self^0 = 1`.
    pub fn pow(&self, mut exp: u64) -> Self {
        if exp == 0 {
            return Self::one();
        }
        let mut base = self.clone();
        let mut acc = Self::one();
        while exp > 0 {
            if exp & 1 == 1 {
                acc = acc.mul(&base);
            }
            exp >>= 1;
            if exp > 0 {
                base = base.mul(&base);
            }
        }
        acc
    }

    /// True if the value is `+1`.
    pub(crate) fn is_one(&self) -> bool {
        self.sign == Sign::Pos && self.limbs.len() == 1 && self.limbs[0] == 1
    }

    /// Least-significant limb, or `0`.
    pub(crate) fn low_word(&self) -> Word {
        self.limbs.first().copied().unwrap_or(0)
    }

    /// Little-endian limbs with the given sign.
    pub(crate) fn from_le_words(sign: Sign, words: &[Word]) -> Self {
        Self::from_limbs(sign, words.to_vec())
    }

    /// `self << bits`.
    pub(crate) fn shl(&self, bits: usize) -> Self {
        if self.is_zero() || bits == 0 {
            return self.clone();
        }
        let woff = bits / WORD_BIT_SIZE;
        let boff = bits % WORD_BIT_SIZE;
        let mut out = vec![0; self.limbs.len() + woff + 1];
        if boff == 0 {
            out[woff..woff + self.limbs.len()].copy_from_slice(&self.limbs);
        } else {
            let mut c = 0;
            for (i, &w) in self.limbs.iter().enumerate() {
                out[woff + i] = (w << boff) | c;
                c = w >> (WORD_BIT_SIZE - boff);
            }
            out[woff + self.limbs.len()] = c;
        }
        Self::from_limbs(self.sign, out)
    }

    /// Convert to an `ExactNum` at `(p, rm)`.
    ///
    /// Packs the little-endian limbs into a mantissa in linear time. The
    /// previous limb-by-limb `ldexp`/`add` walk was quadratic in the limb
    /// count and dominated bit-burst `exp` at dest widths past ~2^16.
    pub fn to_exact_num(&self, p: usize, rm: RoundingMode) -> ExactNum {
        match self.to_exact_num_number(p, rm) {
            Ok(n) => ExactNum::from(n),
            Err(_) => self.to_exact_num_accum(p, rm),
        }
    }

    /// Same conversion, as an `ExactNumNumber`.
    pub(crate) fn to_exact_num_number(
        &self,
        p: usize,
        rm: RoundingMode,
    ) -> Result<crate::num::ExactNumNumber, Error> {
        if self.is_zero() {
            return crate::num::ExactNumNumber::from_word(0, p);
        }
        // `from_words` takes a little-endian (LSW-first) full-width mantissa.
        // Passing exponent = work (the unnormalized bit width) makes the
        // post-normalization exponent equal `bit_length()`, i.e. `|x| < 2^{bit_len}`.
        let work = self.limbs.len().saturating_mul(WORD_BIT_SIZE);
        let mut x =
            crate::num::ExactNumNumber::from_words(&self.limbs, self.sign, work as Exponent)?;
        x.set_precision(p, rm)?;
        Ok(x)
    }

    /// Slow limb-by-limb walk; fallback if the packed conversion cannot allocate.
    fn to_exact_num_accum(&self, p: usize, rm: RoundingMode) -> ExactNum {
        if self.is_zero() {
            return ExactNum::from_u8(0, p);
        }
        let work = p.max(self.bit_length()).saturating_add(WORD_BIT_SIZE);
        let mut acc = ExactNum::from_u8(0, work);
        for &limb in self.limbs.iter().rev() {
            acc = acc.ldexp(WORD_BIT_SIZE as Exponent, work, RoundingMode::None);
            acc = acc.add(&ExactNum::from_word(limb, work), work, RoundingMode::None);
        }
        if self.sign == Sign::Neg {
            acc = acc.neg();
        }
        let _ = acc.set_precision(p, rm);
        acc
    }

    /// Decimal string of `self` (sign plus digits, no exponent).
    pub fn to_dec_string(&self) -> alloc::string::String {
        if self.is_zero() {
            return alloc::string::String::from("0");
        }
        let ten = Self::from_u64(10);
        let mut n = if self.is_negative() { self.neg() } else { self.clone() };
        let mut digits = alloc::vec::Vec::new();
        while !n.is_zero() {
            let (q, r) = n.div_rem(&ten).unwrap_or((Self::zero(), Self::zero()));
            digits.push(b'0' + (r.low_word() as u8));
            n = q;
        }
        if self.is_negative() {
            digits.push(b'-');
        }
        digits.reverse();
        alloc::string::String::from_utf8(digits)
            .unwrap_or_else(|_| alloc::string::String::from("0"))
    }

    /// Parse a decimal integer string (`-?[0-9]+`). Empty or non-digits is `None`.
    pub fn from_dec_string(s: &str) -> Option<Self> {
        let s = s.trim();
        if s.is_empty() {
            return None;
        }
        let (neg, digits) = if let Some(rest) = s.strip_prefix('-') {
            (true, rest)
        } else {
            (false, s.strip_prefix('+').unwrap_or(s))
        };
        if digits.is_empty() || !digits.bytes().all(|c| c.is_ascii_digit()) {
            return None;
        }
        let ten = Self::from_u64(10);
        let mut acc = Self::zero();
        for b in digits.bytes() {
            acc = acc.mul(&ten).add(&Self::from_u64((b - b'0') as u64));
        }
        Some(if neg { acc.neg() } else { acc })
    }

    /// Truncate a finite `ExactNum` toward zero. `None` if Inf or NaN.
    pub fn from_exact_num(x: &ExactNum) -> Option<Self> {
        if x.is_nan() || x.is_inf() {
            return None;
        }
        // |x| < 1 truncates to zero
        if x.exponent().is_none_or(|e| e <= 0) {
            return Some(Self::zero());
        }
        let p = x.precision().unwrap_or(DEFAULT_P).max(WORD_BIT_SIZE + 8);
        let two_w = ExactNum::from_u8(1, p).ldexp(WORD_BIT_SIZE as Exponent, p, RoundingMode::None);
        let sign = if x.is_negative() { Sign::Neg } else { Sign::Pos };
        // The lowest limb is the integer part of `x mod 2^w` (the remainder keeps the sign of `x`, and the
        // limb is read from its magnitude); the rest comes from the truncated quotient.
        let limbs = truncated_limbs(x, &two_w, p)?;
        Some(Self::from_limbs(sign, limbs))
    }

    /// Little-endian limbs, no leading zeros.
    pub(crate) fn as_limbs(&self) -> &[Word] {
        &self.limbs
    }

    /// `self >> bits`, truncated toward zero.
    pub fn shr_towards_zero(&self, bits: usize) -> Self {
        if bits == 0 || self.is_zero() {
            return self.clone();
        }
        let w = bits / WORD_BIT_SIZE;
        let b = bits % WORD_BIT_SIZE;
        if w >= self.limbs.len() {
            return Self::zero();
        }
        let src = &self.limbs[w..];
        // Each result limb is written once from the two source limbs it straddles.
        let out: Vec<Word> = (0..src.len())
            .map(|i| {
                let hi = match src.get(i + 1) {
                    Some(h) if b > 0 => h << (WORD_BIT_SIZE - b),
                    _ => 0,
                };
                (src[i] >> b) | hi
            })
            .collect();
        Self::from_limbs(self.sign, out)
    }
}

/// The base-`two_w` limbs (lowest first) of `trunc(|x|)` for a finite `x` with `|x| >= 1`; `two_w` is
/// `2^w`. `|x| < 2^e` fits in `ceil(e / w)` limbs, so the loop has that fixed trip count: arithmetic
/// that never brings the quotient to zero (a `two_w` of 0 or 1, say) cannot grow the limbs without
/// end. A non-zero quotient left after the last limb means that arithmetic failed: `None`.
fn truncated_limbs(x: &ExactNum, two_w: &ExactNum, p: usize) -> Option<Vec<Word>> {
    let n_limbs = usize::try_from(x.exponent()?).ok()?.div_ceil(WORD_BIT_SIZE);
    let mut limbs = vec![exact_word_lt_base(&x.rem(two_w))];
    let mut cur = ExactNum::div_ref(x, two_w, p, RoundingMode::ToZero)
        .int()
        .abs();
    for _ in 1..n_limbs {
        let r = cur.rem(two_w);
        limbs.push(exact_word_lt_base(&r));
        cur = cur.div(two_w, p, RoundingMode::ToZero).int();
    }
    if !cur.is_zero() {
        return None;
    }
    Some(limbs)
}

fn exact_word_lt_base(x: &ExactNum) -> Word {
    if x.is_zero() {
        return 0;
    }
    let Some((m, _n, _s, e, _)) = x.value_parts() else {
        return 0;
    };
    let pbuf = m.len() * WORD_BIT_SIZE;
    let shift = e - pbuf as i64;
    shift_le_low_word(m, shift)
}

fn shift_le_low_word(m: &[Word], shift: i64) -> Word {
    if m.is_empty() {
        return 0;
    }
    if shift >= 0 {
        let s = shift as usize;
        let woff = s / WORD_BIT_SIZE;
        let boff = s % WORD_BIT_SIZE;
        let lo = m.get(woff).copied().unwrap_or(0);
        if boff == 0 {
            lo
        } else {
            let hi = m.get(woff + 1).copied().unwrap_or(0);
            (lo << boff) | (hi >> (WORD_BIT_SIZE - boff))
        }
    } else {
        let r = (-shift) as usize;
        let woff = r / WORD_BIT_SIZE;
        let boff = r % WORD_BIT_SIZE;
        let lo = m.get(woff).copied().unwrap_or(0);
        if boff == 0 {
            lo
        } else {
            let hi = m.get(woff + 1).copied().unwrap_or(0);
            (lo >> boff) | (hi << (WORD_BIT_SIZE - boff))
        }
    }
}

impl PartialEq for ExactInt {
    fn eq(&self, other: &Self) -> bool {
        self.soft_cmp(other) == Ordering::Equal
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shl_shr_roundtrip() {
        let a = ExactInt::from_u64(0x1234_5678_9abc_def0);
        let s = a.shl(70).mul(&ExactInt::from_u64(3));
        let back = s.shr_towards_zero(70);
        assert_eq!(back, a.mul(&ExactInt::from_u64(3)));
        assert!(ExactInt::from_u64(7).shr_towards_zero(3).is_zero());
        assert_eq!(
            ExactInt::from_i64(-32).shr_towards_zero(2),
            ExactInt::from_i64(-8)
        );
        assert_eq!(a.shl(0), a);
        assert!(a.shr_towards_zero(10_000).is_zero());
    }

    /// Golds: Python big-int `>>` on |x| with sign restored, cross-checked by mpmath 1.3.0 `int(ldexp(x, -k))` at 2000 bits.
    #[test]
    fn shr_towards_zero_matches_mpmath_golds() {
        const ROWS: &[(&str, usize, &str)] = &[
            ("265613988875874769338781322035779626829233452653394495974574961739092490901302182994384699056346", 0, "265613988875874769338781322035779626829233452653394495974574961739092490901302182994384699056346"),
            ("265613988875874769338781322035779626829233452653394495974574961739092490901302182994384699056346", 1, "132806994437937384669390661017889813414616726326697247987287480869546245450651091497192349528173"),
            ("265613988875874769338781322035779626829233452653394495974574961739092490901302182994384699056346", 5, "8300437152371086541836916313618113338413545395418577999205467554346640340665693218574521845510"),
            ("265613988875874769338781322035779626829233452653394495974574961739092490901302182994384699056346", 63, "28797926378176403265074915749429212066859228703256175227484055087293913848900"),
            ("265613988875874769338781322035779626829233452653394495974574961739092490901302182994384699056346", 64, "14398963189088201632537457874714606033429614351628087613742027543646956924450"),
            ("265613988875874769338781322035779626829233452653394495974574961739092490901302182994384699056346", 65, "7199481594544100816268728937357303016714807175814043806871013771823478462225"),
            ("265613988875874769338781322035779626829233452653394495974574961739092490901302182994384699056346", 127, "1561138717114823538298347261418894943943004057247603756397"),
            ("265613988875874769338781322035779626829233452653394495974574961739092490901302182994384699056346", 128, "780569358557411769149173630709447471971502028623801878198"),
            ("265613988875874769338781322035779626829233452653394495974574961739092490901302182994384699056346", 129, "390284679278705884574586815354723735985751014311900939099"),
            ("265613988875874769338781322035779626829233452653394495974574961739092490901302182994384699056346", 191, "84629499432356251143987332790052357727"),
            ("265613988875874769338781322035779626829233452653394495974574961739092490901302182994384699056346", 300, "130392"),
            ("-5817092933824343165432524003391691164919859649719340532627567207607656859034356995566589707894210757866827613621721127496191250", 0, "-5817092933824343165432524003391691164919859649719340532627567207607656859034356995566589707894210757866827613621721127496191250"),
            ("-5817092933824343165432524003391691164919859649719340532627567207607656859034356995566589707894210757866827613621721127496191250", 1, "-2908546466912171582716262001695845582459929824859670266313783603803828429517178497783294853947105378933413806810860563748095625"),
            ("-5817092933824343165432524003391691164919859649719340532627567207607656859034356995566589707894210757866827613621721127496191250", 5, "-181784154182010723919766375105990348903745614053729391644611475237739276844823656111455928371694086183338362925678785234255976"),
            ("-5817092933824343165432524003391691164919859649719340532627567207607656859034356995566589707894210757866827613621721127496191250", 63, "-630690479640242953884043866315223070633472261989798500635874112372131472681406811174098266016878632592637537"),
            ("-5817092933824343165432524003391691164919859649719340532627567207607656859034356995566589707894210757866827613621721127496191250", 64, "-315345239820121476942021933157611535316736130994899250317937056186065736340703405587049133008439316296318768"),
            ("-5817092933824343165432524003391691164919859649719340532627567207607656859034356995566589707894210757866827613621721127496191250", 65, "-157672619910060738471010966578805767658368065497449625158968528093032868170351702793524566504219658148159384"),
            ("-5817092933824343165432524003391691164919859649719340532627567207607656859034356995566589707894210757866827613621721127496191250", 127, "-34189799409593810678020064275276503335837544375575702730719652325997428637418337151588095"),
            ("-5817092933824343165432524003391691164919859649719340532627567207607656859034356995566589707894210757866827613621721127496191250", 128, "-17094899704796905339010032137638251667918772187787851365359826162998714318709168575794047"),
            ("-5817092933824343165432524003391691164919859649719340532627567207607656859034356995566589707894210757866827613621721127496191250", 129, "-8547449852398452669505016068819125833959386093893925682679913081499357159354584287897023"),
            ("-5817092933824343165432524003391691164919859649719340532627567207607656859034356995566589707894210757866827613621721127496191250", 191, "-1853432739836261319302466666837631586389970146892710989193121000318357"),
            ("-5817092933824343165432524003391691164919859649719340532627567207607656859034356995566589707894210757866827613621721127496191250", 300, "-2855665290846666351956751064424833016"),
            ("6277101735386680763835789423207666416102355444464034512895", 0, "6277101735386680763835789423207666416102355444464034512895"),
            ("6277101735386680763835789423207666416102355444464034512895", 1, "3138550867693340381917894711603833208051177722232017256447"),
            ("6277101735386680763835789423207666416102355444464034512895", 5, "196159429230833773869868419475239575503198607639501078527"),
            ("6277101735386680763835789423207666416102355444464034512895", 63, "680564733841876926926749214863536422911"),
            ("6277101735386680763835789423207666416102355444464034512895", 64, "340282366920938463463374607431768211455"),
            ("6277101735386680763835789423207666416102355444464034512895", 65, "170141183460469231731687303715884105727"),
            ("6277101735386680763835789423207666416102355444464034512895", 127, "36893488147419103231"),
            ("6277101735386680763835789423207666416102355444464034512895", 128, "18446744073709551615"),
            ("6277101735386680763835789423207666416102355444464034512895", 129, "9223372036854775807"),
            ("6277101735386680763835789423207666416102355444464034512895", 191, "1"),
            ("6277101735386680763835789423207666416102355444464034512895", 300, "0"),
            ("-1361129467683753853853498429727072845824", 0, "-1361129467683753853853498429727072845824"),
            ("-1361129467683753853853498429727072845824", 1, "-680564733841876926926749214863536422912"),
            ("-1361129467683753853853498429727072845824", 5, "-42535295865117307932921825928971026432"),
            ("-1361129467683753853853498429727072845824", 63, "-147573952589676412928"),
            ("-1361129467683753853853498429727072845824", 64, "-73786976294838206464"),
            ("-1361129467683753853853498429727072845824", 65, "-36893488147419103232"),
            ("-1361129467683753853853498429727072845824", 127, "-8"),
            ("-1361129467683753853853498429727072845824", 128, "-4"),
            ("-1361129467683753853853498429727072845824", 129, "-2"),
            ("-1361129467683753853853498429727072845824", 191, "0"),
            ("-1361129467683753853853498429727072845824", 300, "0"),
            ("338770000845734292516042252062085074415", 0, "338770000845734292516042252062085074415"),
            ("338770000845734292516042252062085074415", 1, "169385000422867146258021126031042537207"),
            ("338770000845734292516042252062085074415", 5, "10586562526429196641126320376940158575"),
            ("338770000845734292516042252062085074415", 63, "36729517088986129440"),
            ("338770000845734292516042252062085074415", 64, "18364758544493064720"),
            ("338770000845734292516042252062085074415", 65, "9182379272246532360"),
            ("338770000845734292516042252062085074415", 127, "1"),
            ("338770000845734292516042252062085074415", 128, "0"),
            ("338770000845734292516042252062085074415", 129, "0"),
            ("338770000845734292516042252062085074415", 191, "0"),
            ("338770000845734292516042252062085074415", 300, "0"),
        ];
        for &(x, k, want) in ROWS {
            let want = ExactInt::from_dec_string(want);
            assert!(want.is_some(), "gold row {x} >> {k} parses");
            let got = ExactInt::from_dec_string(x).map(|v| v.shr_towards_zero(k));
            assert_eq!(got, want, "{x} >> {k}");
        }
    }

    #[test]
    fn to_exact_num_matches_limb_accum() {
        let p = 256;
        let cases = [
            ExactInt::zero(),
            ExactInt::from_u64(1),
            ExactInt::from_u64(3),
            ExactInt::from_u64(u64::MAX),
            ExactInt::from_u128(1u128 << 64),
            ExactInt::from_u128(u128::MAX),
            ExactInt::from_u64(5).shl(200),
            ExactInt::from_u64(0x1234_5678_9abc_def0)
                .shl(70)
                .mul(&ExactInt::from_u64(3)),
        ];
        for n in cases {
            for neg in [false, true] {
                let x = if neg && !n.is_zero() { n.neg() } else { n.clone() };
                let fast = x.to_exact_num(p, RoundingMode::ToEven);
                let slow = x.to_exact_num_accum(p, RoundingMode::ToEven);
                assert_eq!(fast, slow, "n={x:?} neg={neg}");
            }
        }
    }

    #[test]
    fn to_exact_num_wide_is_linear() {
        // A ~2^18-bit integer. The accum walk is thousands of wide adds;
        // the packed conversion should finish in a few milliseconds.
        let n = ExactInt::from_u64(0xdead_beef_cafe_babe).shl(1 << 18);
        let t0 = std::time::Instant::now();
        let y = n.to_exact_num(1 << 18, RoundingMode::ToEven);
        let ms = t0.elapsed().as_millis();
        assert!(!y.is_nan() && !y.is_zero(), "wide conversion failed");
        assert!(
            crate::test_time::within_ms(ms, 200),
            "to_exact_num of a 2^18-bit integer took {ms} ms"
        );
        let slow = n.to_exact_num_accum(256, RoundingMode::ToEven);
        let fast = n.to_exact_num(256, RoundingMode::ToEven);
        assert_eq!(fast, slow);
    }

    #[test]
    fn div_rem_with_zero_operands_no_panic() {
        // Zero divisor is None; zero dividend is (0, 0); the limb kernel with an empty
        // dividend and a non-empty divisor returns an empty quotient and remainder.
        assert!(ExactInt::one().div_rem(&ExactInt::zero()).is_none());
        let zero_quot = ExactInt::zero().div_rem(&ExactInt::one());
        assert!(matches!(zero_quot, Some((q, r)) if q.is_zero() && r.is_zero()));
        assert_eq!(ExactInt::div_rem_abs(&[], &[1]), (Vec::new(), Vec::new()));
    }

    #[test]
    fn from_exact_num_truncates_toward_zero() {
        let p = 192;
        let half = ExactNum::from_u8(1, p).ldexp(-1, p, RoundingMode::None);
        let mags: [u128; 7] =
            [0, 1, 7, u64::MAX as u128, 1u128 << 64, (1u128 << 64) + 5, 1u128 << 100];
        for m in mags {
            for neg in [false, true] {
                for frac in [false, true] {
                    let mut x = ExactInt::from_u128(m).to_exact_num(p, RoundingMode::ToEven);
                    if frac {
                        x = x.add(&half, p, RoundingMode::None);
                    }
                    let want =
                        if neg { ExactInt::from_u128(m).neg() } else { ExactInt::from_u128(m) };
                    let x = if neg { x.neg() } else { x };
                    assert_eq!(
                        ExactInt::from_exact_num(&x),
                        Some(want),
                        "m={m} neg={neg} frac={frac}"
                    );
                }
            }
        }
        assert_eq!(
            ExactInt::from_exact_num(&ExactNum::new(p)),
            Some(ExactInt::zero())
        );
        assert!(ExactInt::from_exact_num(&crate::INF_POS).is_none());
    }

    /// `from_exact_num` at the limb boundaries: 2^64 − 1 is one limb, 2^64 two, 2^128 − 1 two,
    /// 2^128 three, 3·2^190 (< 2^192) three, 2^192 four; each with and without a fraction ½ and with either sign, against
    /// the `ExactInt` it came from.
    #[test]
    fn from_exact_num_limb_count_boundaries() -> Result<(), crate::Error> {
        let p = 320;
        let two = || ExactInt::from_i64(2);
        let mk = |i: usize| match i {
            0 => ExactInt::from_u128(u128::from(u64::MAX)),
            1 => ExactInt::from_u128(1u128 << 64),
            2 => ExactInt::from_u128(u128::MAX),
            3 => ExactInt::from_u128(1u128 << 127).mul(&two()),
            4 => ExactInt::from_i64(3).mul(&two().pow(190)),
            _ => two().pow(192),
        };
        let limbs = [1usize, 2, 2, 3, 3, 4];
        for (i, &want_limbs) in limbs.iter().enumerate() {
            for neg in [false, true] {
                for frac in [false, true] {
                    let mut x = mk(i).to_exact_num(p, RoundingMode::ToEven);
                    if frac {
                        x = x.add(
                            &ExactNum::from_u8(1, p).ldexp(-1, p, RoundingMode::None),
                            p,
                            RoundingMode::None,
                        );
                    }
                    let x = if neg { -x } else { x };
                    let want = if neg { mk(i).neg() } else { mk(i) };
                    let got = ExactInt::from_exact_num(&x).ok_or(crate::Error::InvalidArgument)?;
                    assert_eq!(
                        got.limbs.len(),
                        want_limbs,
                        "case {i} neg={neg} frac={frac}"
                    );
                    assert_eq!(got, want, "case {i} neg={neg} frac={frac}");
                }
            }
        }
        Ok(())
    }

    /// `truncated_limbs` takes exactly `ceil(e / 64)` limbs for `|x| < 2^e` (2^64 + ½ is two limbs,
    /// 2^64 − ½ one), and gives `None` instead of looping when the divisor never brings the
    /// quotient to zero (0 or 1 for `two_w`: what a broken `ldexp` hands it).
    #[test]
    fn truncated_limbs_has_a_fixed_trip_count() {
        let p = 128;
        let tw = || ExactNum::from_u8(1, p).ldexp(64, p, RoundingMode::None);
        let two_w = tw();
        let half = ExactNum::from_u8(1, p).ldexp(-1, p, RoundingMode::None);
        let x = tw().add(&half, p, RoundingMode::None);
        assert_eq!(truncated_limbs(&x, &two_w, p), Some(vec![0, 1]));
        let y = tw().sub(&half, p, RoundingMode::None);
        assert_eq!(truncated_limbs(&y, &two_w, p), Some(vec![Word::MAX]));
        let z = ExactNum::from_u8(5, p);
        assert_eq!(truncated_limbs(&z, &two_w, p), Some(vec![5]));
        assert_eq!(truncated_limbs(&x, &ExactNum::from_u8(1, p), p), None);
        assert_eq!(truncated_limbs(&x, &ExactNum::new(p), p), None);
        assert_eq!(truncated_limbs(&z, &ExactNum::from_u8(1, p), p), None);
    }

    #[test]
    fn exact_int_fact20_gcd_pow2_divrem() -> Result<(), crate::Error> {
        let mut fact = ExactInt::one();
        for k in 2..=20 {
            fact = fact.mul(&ExactInt::from_i64(k));
        }
        assert_eq!(fact, ExactInt::from_i64(2_432_902_008_176_640_000));

        let g = ExactInt::from_i64(48).gcd(&ExactInt::from_i64(18));
        assert_eq!(g, ExactInt::from_i64(6));

        let p2 = ExactInt::from_i64(2).pow(100);
        assert_eq!(p2, ExactInt::from_u128(1u128 << 100));
        assert_eq!(
            p2,
            ExactInt::from_u128(1_267_650_600_228_229_401_496_703_205_376)
        );
        assert_eq!(p2.bit_length(), 101);

        let (q, r) = ExactInt::from_i64(17)
            .div_rem(&ExactInt::from_i64(5))
            .ok_or(crate::Error::InvalidArgument)?;
        assert_eq!(q, ExactInt::from_i64(3));
        assert_eq!(r, ExactInt::from_i64(2));

        let n = ExactInt::from_i64(-3);
        let f = n.to_exact_num(64, RoundingMode::ToEven);
        assert_eq!(ExactInt::from_exact_num(&f), Some(n));
        assert!(ExactInt::from_exact_num(&crate::NAN).is_none());
        assert!(ExactInt::from_i64(1).div_rem(&ExactInt::zero()).is_none());
        Ok(())
    }
}
