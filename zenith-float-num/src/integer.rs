//! Arbitrary-precision signed integers on `Word` limbs.

use crate::common::util::add_carry;
use crate::common::util::sub_borrow;
use crate::defs::DoubleWord;
use crate::defs::RoundingMode;
use crate::defs::Sign;
use crate::defs::Word;
use crate::defs::DEFAULT_P;
use crate::defs::WORD_BIT_SIZE;
use crate::ExactNum;
use alloc::vec;
use alloc::vec::Vec;
use core::cmp::Ordering;

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
    pub fn cmp(&self, other: &Self) -> Ordering {
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
        for i in 0..n {
            let x = a.get(i).copied().unwrap_or(0);
            let y = b.get(i).copied().unwrap_or(0);
            c = add_carry(x, y, c, &mut out[i]);
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
        let bits = (a.len() - 1) * WORD_BIT_SIZE
            + (WORD_BIT_SIZE - a.last().unwrap().leading_zeros() as usize);
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

    /// Convert to an `ExactNum` at `(p, rm)`.
    pub fn to_exact_num(&self, p: usize, rm: RoundingMode) -> ExactNum {
        if self.is_zero() {
            return ExactNum::from_u8(0, p);
        }
        let work = p.max(self.bit_length()).saturating_add(WORD_BIT_SIZE);
        let mut acc = ExactNum::from_u8(0, work);
        for &limb in self.limbs.iter().rev() {
            acc = acc.ldexp(WORD_BIT_SIZE as i32, work, RoundingMode::None);
            acc = acc.add(&ExactNum::from_word(limb, work), work, RoundingMode::None);
        }
        if self.sign == Sign::Neg {
            acc = acc.neg();
        }
        let _ = acc.set_precision(p, rm);
        acc
    }

    /// Truncate a finite `ExactNum` toward zero. `None` if Inf or NaN.
    pub fn from_exact_num(x: &ExactNum) -> Option<Self> {
        if x.is_nan() || x.is_inf() {
            return None;
        }
        let v = x.int();
        if v.is_zero() {
            return Some(Self::zero());
        }
        let p = v.precision().unwrap_or(DEFAULT_P).max(WORD_BIT_SIZE + 8);
        let two_w = ExactNum::from_u8(1, p).ldexp(WORD_BIT_SIZE as i32, p, RoundingMode::None);
        let mut cur = v.abs();
        let mut limbs = Vec::new();
        while !cur.is_zero() {
            let r = cur.rem(&two_w);
            limbs.push(exact_word_lt_base(&r));
            cur = cur.div(&two_w, p, RoundingMode::ToZero).int();
        }
        let sign = if v.is_negative() { Sign::Neg } else { Sign::Pos };
        Some(Self::from_limbs(sign, limbs))
    }
}

fn exact_word_lt_base(x: &ExactNum) -> Word {
    if x.is_zero() {
        return 0;
    }
    let Some((m, _n, _s, e, _)) = x.as_raw_parts() else {
        return 0;
    };
    let pbuf = m.len() * WORD_BIT_SIZE;
    let shift = e as isize - pbuf as isize;
    shift_le_low_word(m, shift)
}

fn shift_le_low_word(m: &[Word], shift: isize) -> Word {
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
        self.cmp(other) == Ordering::Equal
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_int_fact20_gcd_pow2_divrem() {
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
            .unwrap();
        assert_eq!(q, ExactInt::from_i64(3));
        assert_eq!(r, ExactInt::from_i64(2));

        let n = ExactInt::from_i64(-3);
        let f = n.to_exact_num(64, RoundingMode::ToEven);
        assert_eq!(ExactInt::from_exact_num(&f), Some(n));
        assert!(ExactInt::from_exact_num(&crate::NAN).is_none());
        assert!(ExactInt::from_i64(1).div_rem(&ExactInt::zero()).is_none());
    }
}
