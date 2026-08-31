//! Exact rationals `num/den` with `ExactNum` limbs, reduced to lowest terms.

use crate::defs::RoundingMode;
use crate::defs::DEFAULT_P;
use crate::defs::WORD_BIT_SIZE;
use crate::ExactInt;
use crate::ExactNum;
use crate::Radix;
use crate::NAN;
use alloc::string::String;
use alloc::vec::Vec;
use core::cmp::Ordering;

/// Exact rational `num/den`. Distinct from [`ExactNum`], which is floating-point.
///
/// Both parts are integer-valued `ExactNum`s after [`Self::new`]. The
/// denominator is positive. A zero denominator or a non-finite part is stored
/// as `NaN/1`.
#[derive(Clone, Debug)]
pub struct ExactRational {
    num: ExactNum,
    den: ExactNum,
}

fn rat_one() -> ExactNum {
    ExactNum::from_u8(1, DEFAULT_P)
}

fn work_p(xs: &[&ExactNum]) -> usize {
    let mut p = DEFAULT_P;
    for x in xs {
        if let Some(px) = x.precision() {
            p = p.max(px);
        }
    }
    p.saturating_add(WORD_BIT_SIZE)
}

fn is_int(x: &ExactNum) -> bool {
    if x.is_nan() || x.is_inf() {
        return false;
    }
    x.fract().is_zero()
}

fn gcd_int(mut a: ExactNum, mut b: ExactNum) -> ExactNum {
    a = a.abs();
    b = b.abs();
    while !b.is_zero() {
        let r = a.rem(&b);
        a = b;
        b = r;
    }
    if a.is_zero() {
        rat_one()
    } else {
        a
    }
}

fn reduce(mut num: ExactNum, mut den: ExactNum) -> ExactRational {
    if num.is_nan() || den.is_nan() || num.is_inf() || den.is_inf() || den.is_zero() {
        return ExactRational {
            num: NAN,
            den: rat_one(),
        };
    }
    if den.is_negative() {
        num = num.neg();
        den = den.neg();
    }
    if is_int(&num) && is_int(&den) {
        let g = gcd_int(num.clone(), den.clone());
        if !g.is_zero() && g.cmp(&rat_one()) != Some(0) {
            let p = work_p(&[&num, &den, &g]);
            num = num.div(&g, p, RoundingMode::ToEven);
            den = den.div(&g, p, RoundingMode::ToEven);
        }
    }
    ExactRational { num, den }
}

impl ExactRational {
    /// `num/den` reduced to lowest terms. Negative `den` moves the sign to `num`.
    /// Zero `den` or a non-finite part is `NaN`.
    pub fn new(num: ExactNum, den: ExactNum) -> Self {
        reduce(num, den)
    }

    /// Integer ratio `n/d` at default limb precision.
    pub fn from_i64(n: i64, d: i64) -> Self {
        Self::new(
            ExactNum::from_i64(n, DEFAULT_P),
            ExactNum::from_i64(d, DEFAULT_P),
        )
    }

    /// Numerator (integer-valued after construction).
    pub fn num(&self) -> &ExactNum {
        &self.num
    }

    /// Denominator (positive integer-valued after construction).
    pub fn den(&self) -> &ExactNum {
        &self.den
    }

    /// True if a part is NaN or the denominator was zero.
    pub fn is_nan(&self) -> bool {
        self.num.is_nan() || self.den.is_nan()
    }

    /// True if the denominator is `1`.
    pub fn is_integer(&self) -> bool {
        !self.is_nan() && self.den.cmp(&rat_one()) == Some(0)
    }

    /// `self + rhs` as an exact rational.
    pub fn add(&self, rhs: &Self) -> Self {
        if self.is_nan() || rhs.is_nan() {
            return Self::new(NAN, rat_one());
        }
        let p = work_p(&[&self.num, &self.den, &rhs.num, &rhs.den]);
        let rm = RoundingMode::None;
        let ad = self.num.mul(&rhs.den, p, rm);
        let bc = rhs.num.mul(&self.den, p, rm);
        let num = ad.add(&bc, p, rm);
        let den = self.den.mul(&rhs.den, p, rm);
        Self::new(num, den)
    }

    /// `self - rhs`.
    pub fn sub(&self, rhs: &Self) -> Self {
        self.add(&Self::new(rhs.num.neg(), rhs.den.clone()))
    }

    /// `self * rhs`.
    pub fn mul(&self, rhs: &Self) -> Self {
        if self.is_nan() || rhs.is_nan() {
            return Self::new(NAN, rat_one());
        }
        let p = work_p(&[&self.num, &self.den, &rhs.num, &rhs.den]);
        let rm = RoundingMode::None;
        Self::new(self.num.mul(&rhs.num, p, rm), self.den.mul(&rhs.den, p, rm))
    }

    /// `self / rhs`. Zero `rhs` is `NaN`.
    pub fn div(&self, rhs: &Self) -> Self {
        if self.is_nan() || rhs.is_nan() || rhs.num.is_zero() {
            return Self::new(NAN, rat_one());
        }
        let p = work_p(&[&self.num, &self.den, &rhs.num, &rhs.den]);
        let rm = RoundingMode::None;
        Self::new(self.num.mul(&rhs.den, p, rm), self.den.mul(&rhs.num, p, rm))
    }

    /// Convert to an `ExactNum` at `(p, rm)`.
    pub fn to_exact_num(&self, p: usize, rm: RoundingMode) -> ExactNum {
        if self.is_nan() {
            return NAN;
        }
        self.num.div(&self.den, p, rm)
    }

    fn quot(&self) -> ExactNum {
        let p = work_p(&[&self.num, &self.den]);
        self.num.div(&self.den, p, RoundingMode::None)
    }

    /// Greatest integer `≤ self` as a rational with denominator `1`.
    pub fn floor(&self) -> Self {
        if self.is_nan() {
            return Self::new(NAN, rat_one());
        }
        Self::new(self.quot().floor(), rat_one())
    }

    /// Least integer `≥ self` as a rational with denominator `1`.
    pub fn ceil(&self) -> Self {
        if self.is_nan() {
            return Self::new(NAN, rat_one());
        }
        Self::new(self.quot().ceil(), rat_one())
    }

    /// Nearest integer, ties away from zero, as a rational with denominator `1`.
    pub fn round(&self) -> Self {
        if self.is_nan() {
            return Self::new(NAN, rat_one());
        }
        Self::new(self.quot().round(0, RoundingMode::FromZero), rat_one())
    }

    /// Build `num/den` from limb integers.
    pub fn from_ints(num: ExactInt, den: ExactInt) -> Self {
        if den.is_zero() {
            return Self::new(NAN, rat_one());
        }
        let p = num
            .bit_length()
            .max(den.bit_length())
            .max(DEFAULT_P)
            .saturating_add(WORD_BIT_SIZE);
        Self::new(
            num.to_exact_num(p, RoundingMode::None),
            den.to_exact_num(p, RoundingMode::None),
        )
    }

    /// Parse a decimal (optional `e`/`E` exponent) as an exact rational.
    ///
    /// `0.1` is `1/10`, not a binary float. Invalid syntax is `None`.
    pub fn parse_exact(s: &str) -> Option<Self> {
        let s = s.trim();
        if s.is_empty() {
            return None;
        }
        let (neg, rest) = match s.as_bytes()[0] {
            b'+' => (false, &s[1..]),
            b'-' => (true, &s[1..]),
            _ => (false, s),
        };
        if rest.is_empty() {
            return None;
        }
        let (mant, exp) = split_dec_exp(rest)?;
        let (int_part, frac_part) = split_dot(mant);
        if int_part.is_empty() && frac_part.is_empty() {
            return None;
        }
        if !int_part.bytes().all(|c| c.is_ascii_digit())
            || !frac_part.bytes().all(|c| c.is_ascii_digit())
        {
            return None;
        }
        let mut digits = String::new();
        digits.push_str(int_part);
        digits.push_str(frac_part);
        if digits.is_empty() || digits.bytes().all(|c| c == b'0') {
            return Some(Self::from_i64(0, 1));
        }
        let num = dec_digits_to_int(&digits)?;
        let exp_adj = exp - frac_part.len() as i32;
        let ten = ExactInt::from_i64(10);
        let (n, d) = if exp_adj >= 0 {
            (num.mul(&ten.pow(exp_adj as u64)), ExactInt::one())
        } else {
            (num, ten.pow((-exp_adj) as u64))
        };
        let mut r = Self::from_ints(n, d);
        if neg {
            r = Self::new(r.num.neg(), r.den.clone());
        }
        Some(r)
    }

    /// Minimum digits in `rdx` for an exact terminating representation. `None` if the
    /// denominator has a prime factor that does not divide the radix.
    pub fn format_exact(&self, rdx: Radix) -> Option<String> {
        if self.is_nan() {
            return None;
        }
        let num = ExactInt::from_exact_num(&self.num)?;
        let den = ExactInt::from_exact_num(&self.den)?;
        format_terminating(&num, &den, rdx)
    }

    /// Exact comparison via cross-multiplication. `None` if either value is NaN.
    pub fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.is_nan() || other.is_nan() {
            return None;
        }
        let p = work_p(&[&self.num, &self.den, &other.num, &other.den]);
        let rm = RoundingMode::None;
        let left = self.num.mul(&other.den, p, rm);
        let right = other.num.mul(&self.den, p, rm);
        match left.cmp(&right) {
            Some(c) if c > 0 => Some(Ordering::Greater),
            Some(c) if c < 0 => Some(Ordering::Less),
            Some(_) => Some(Ordering::Equal),
            None => None,
        }
    }
}

fn split_dec_exp(s: &str) -> Option<(&str, i32)> {
    let bytes = s.as_bytes();
    let mut epos = None;
    for (i, c) in bytes.iter().enumerate() {
        if *c == b'e' || *c == b'E' {
            epos = Some(i);
            break;
        }
    }
    match epos {
        None => Some((s, 0)),
        Some(i) => {
            if i == 0 {
                return None;
            }
            let exp_s = &s[i + 1..];
            if exp_s.is_empty() {
                return None;
            }
            let exp: i32 = exp_s.parse().ok()?;
            Some((&s[..i], exp))
        }
    }
}

fn split_dot(s: &str) -> (&str, &str) {
    match s.find('.') {
        Some(i) => (&s[..i], &s[i + 1..]),
        None => (s, ""),
    }
}

fn dec_digits_to_int(s: &str) -> Option<ExactInt> {
    let ten = ExactInt::from_i64(10);
    let mut v = ExactInt::zero();
    for c in s.bytes() {
        if !c.is_ascii_digit() {
            return None;
        }
        v = v.mul(&ten).add(&ExactInt::from_u64((c - b'0') as u64));
    }
    Some(v)
}

fn format_terminating(num: &ExactInt, den: &ExactInt, rdx: Radix) -> Option<String> {
    if den.is_zero() {
        return None;
    }
    let radix = ExactInt::from_u64(rdx.value() as u64);
    let mut pk = ExactInt::one();
    let max_k = den.bit_length().saturating_add(8);
    for k in 0..=max_k {
        if let Some((scale, rem)) = pk.div_rem(den) {
            if rem.is_zero() {
                let digits = num.mul(&scale);
                return Some(format_fixed(&digits, k, rdx));
            }
        }
        pk = pk.mul(&radix);
    }
    None
}

fn format_fixed(n: &ExactInt, frac_digits: usize, rdx: Radix) -> String {
    let neg = n.is_negative();
    let abs = if neg { n.neg() } else { n.clone() };
    let mut digits = to_radix_digits(&abs, rdx);
    let mut k = frac_digits;
    while k > 0 && digits.ends_with('0') {
        digits.pop();
        k -= 1;
    }
    if digits.is_empty() {
        digits.push('0');
    }
    let mut out = String::new();
    if neg && digits != "0" {
        out.push('-');
    }
    if k == 0 {
        out.push_str(&digits);
        return out;
    }
    if digits.len() <= k {
        out.push('0');
        out.push('.');
        for _ in 0..(k - digits.len()) {
            out.push('0');
        }
        out.push_str(&digits);
    } else {
        let split = digits.len() - k;
        out.push_str(&digits[..split]);
        out.push('.');
        out.push_str(&digits[split..]);
    }
    out
}

fn to_radix_digits(n: &ExactInt, rdx: Radix) -> String {
    if n.is_zero() {
        return String::from("0");
    }
    let r = ExactInt::from_u64(rdx.value() as u64);
    let mut v = if n.is_negative() { n.neg() } else { n.clone() };
    let mut digits = Vec::new();
    while !v.is_zero() {
        let (q, rem) = v
            .div_rem(&r)
            .unwrap_or((ExactInt::zero(), ExactInt::zero()));
        let d = rem.low_word() as u8;
        digits.push(if d < 10 { b'0' + d } else { b'a' + (d - 10) });
        v = q;
    }
    digits.reverse();
    String::from_utf8(digits).unwrap_or_else(|_| String::from("0"))
}

fn exact_num_to_rational(x: &ExactNum) -> Option<ExactRational> {
    if x.is_nan() || x.is_inf() {
        return None;
    }
    if x.is_zero() {
        return Some(ExactRational::from_i64(0, 1));
    }
    let (m, _n, s, e, _) = x.as_raw_parts()?;
    let pbuf = m.len() * WORD_BIT_SIZE;
    let mut num = ExactInt::from_le_words(s, m);
    let shift = e as i64 - pbuf as i64;
    if shift >= 0 {
        num = num.shl(shift as usize);
        Some(ExactRational::from_ints(num, ExactInt::one()))
    } else {
        let den = ExactInt::one().shl((-shift) as usize);
        Some(ExactRational::from_ints(num, den))
    }
}

impl ExactNum {
    /// Parse a decimal that is exact in binary (a dyadic rational). `None` if the
    /// value is not a finite dyadic (so `0.1` is `None` here; use [`ExactRational::parse_exact`]).
    pub fn parse_exact(s: &str) -> Option<Self> {
        let r = ExactRational::parse_exact(s)?;
        let den = ExactInt::from_exact_num(r.den())?;
        if !den.is_one() && !is_pow2(&den) {
            return None;
        }
        let p = den.bit_length().max(WORD_BIT_SIZE);
        Some(r.to_exact_num(p, RoundingMode::None))
    }

    /// Minimum digits in `rdx` that recover `self` exactly when the value is a
    /// terminating expansion in that radix.
    pub fn format_exact(&self, rdx: Radix) -> Option<String> {
        exact_num_to_rational(self)?.format_exact(rdx)
    }
}

fn is_pow2(n: &ExactInt) -> bool {
    if n.is_zero() || n.is_negative() {
        return false;
    }
    n.sub(&ExactInt::one()).bit_length() < n.bit_length()
}

impl PartialEq for ExactRational {
    fn eq(&self, other: &Self) -> bool {
        matches!(self.partial_cmp(other), Some(Ordering::Equal))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Consts;
    use crate::Radix;

    fn r(n: i64, d: i64) -> ExactRational {
        ExactRational::new(
            ExactNum::from_i64(n, DEFAULT_P),
            ExactNum::from_i64(d, DEFAULT_P),
        )
    }

    #[test]
    fn rational_reduce_add_sign_to_float() {
        let half = r(1, 2);
        let two_four = r(2, 4);
        assert_eq!(two_four, half);
        assert_eq!(ExactRational::from_i64(2, 4), ExactRational::from_i64(1, 2));
        assert!(!two_four.is_integer());
        assert_eq!(
            two_four.num().cmp(&ExactNum::from_u8(1, DEFAULT_P)),
            Some(0)
        );
        assert_eq!(
            two_four.den().cmp(&ExactNum::from_u8(2, DEFAULT_P)),
            Some(0)
        );

        let s = r(1, 3).add(&r(1, 6));
        assert_eq!(s, half);

        let pos = ExactRational::from_i64(-3, -4);
        let want = ExactRational::from_i64(3, 4);
        assert_eq!(pos, want);
        assert!(pos.den().is_positive());

        let p = 256;
        let rm = RoundingMode::ToEven;
        let got = ExactRational::from_i64(1, 3).to_exact_num(p, rm);
        let want_f = ExactNum::from_i64(1, p).div(&ExactNum::from_i64(3, p), p, rm);
        assert_eq!(got.cmp(&want_f), Some(0));
        let mut cc = Consts::new().unwrap();
        let parsed = ExactNum::parse(
            "0.33333333333333333333333333333333333333333333333333333333333333333333333333333333",
            Radix::Dec,
            p,
            rm,
            &mut cc,
        );
        assert_eq!(got.cmp(&parsed), Some(0));

        assert!(ExactRational::from_i64(1, 3).partial_cmp(&half) == Some(Ordering::Less));
        assert!(ExactRational::from_i64(2, 2).is_integer());
        assert_eq!(
            ExactRational::from_i64(5, 3).floor(),
            ExactRational::from_i64(1, 1)
        );
        assert_eq!(
            ExactRational::from_i64(5, 3).ceil(),
            ExactRational::from_i64(2, 1)
        );
        assert!(ExactRational::from_i64(1, 0).is_nan());
    }

    #[test]
    fn parse_exact_tenth_half_eighth_sci() {
        let tenth = ExactRational::parse_exact("0.1").unwrap();
        assert_eq!(
            tenth.mul(&ExactRational::from_i64(10, 1)),
            ExactRational::from_i64(1, 1)
        );

        let half = ExactNum::parse_exact("0.5").unwrap();
        let p = half.precision().unwrap();
        let want = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, RoundingMode::ToEven);
        assert_eq!(half.cmp(&want), Some(0));
        assert_eq!(half.ilogb(), Some(-1));
        let (_m, _n, s, _e, _) = half.as_raw_parts().unwrap();
        assert_eq!(s, crate::Sign::Pos);

        assert_eq!(
            ExactRational::parse_exact("0.125")
                .unwrap()
                .format_exact(Radix::Dec)
                .as_deref(),
            Some("0.125")
        );
        assert_eq!(
            ExactNum::parse_exact("0.125")
                .unwrap()
                .format_exact(Radix::Dec)
                .as_deref(),
            Some("0.125")
        );
        assert!(ExactNum::parse_exact("0.1").is_none());

        let mut cc = Consts::new().unwrap();
        let sci = ExactNum::parse("1.5e3", Radix::Dec, 256, RoundingMode::ToEven, &mut cc);
        assert_eq!(sci.cmp(&ExactNum::from_i32(1500, 256)), Some(0));
    }
}
