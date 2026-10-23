//! IEEE-style exponent decomposition without hardware floats.

use crate::common::util::round_p;
use crate::defs::{Error, EXPONENT_MAX, EXPONENT_MIN, WORD_BIT_SIZE};
use crate::num::ExactNumNumber;
use crate::{Exponent, RoundingMode, Sign};

impl ExactNumNumber {
    /// Split into a mantissa in `[0.5, 1)` and a power of two: `self = m · 2^e`.
    pub fn frexp(mut self) -> Result<(Self, Exponent), Error> {
        if self.is_zero() {
            return Ok((self, 0));
        }
        let e = self.exponent();
        self.set_exponent(0);
        Ok((self, e))
    }

    /// `self · 2^n` (same as `ldexp` / integer `scalb`).
    pub fn ldexp(mut self, n: Exponent, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        if self.is_zero() {
            self.set_precision(p, rm)?;
            return Ok(self);
        }
        let e = self.exponent();
        let sum = crate::exp::wide(e) + crate::exp::wide(n);
        match Exponent::try_from(sum) {
            Ok(sum) if (EXPONENT_MIN..=EXPONENT_MAX).contains(&sum) => {
                self.set_exponent(sum);
                self.set_precision(p, rm)?;
                Ok(self)
            }
            _ => {
                if n >= 0 {
                    Err(Error::ExponentOverflow(self.sign()))
                } else {
                    // Below the range: scale by 2^EXPONENT_MIN in a multiplication, which
                    // rounds onto the subnormal grid (or to zero / the least value) by mode.
                    // A scale below EXPONENT_MIN would underflow even before that product, so it
                    // is clamped (the result is then below the grid by far more than p bits).
                    let s1 = crate::exp::narrow(sum - crate::exp::wide(EXPONENT_MIN));
                    let s1 = s1.max(EXPONENT_MIN);
                    self.set_exponent(s1);
                    let mut b = Self::from_word(1, WORD_BIT_SIZE)?;
                    b.set_exponent(EXPONENT_MIN + 1);
                    self.mul(&b, p, rm)
                }
            }
        }
    }

    /// Binary `logb`: `floor(log2(|self|))` as an integer-valued float.
    pub fn logb(&self, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        if self.is_zero() {
            return Err(Error::ExponentOverflow(crate::Sign::Neg));
        }
        let e = self.exponent().saturating_sub(1);
        let mut v = Self::from_word(e.unsigned_abs() as crate::Word, p)?;
        if e < 0 {
            v.set_sign(Sign::Neg);
        }
        v.set_precision(p, rm)?;
        Ok(v)
    }

    /// Integer `logb`.
    pub fn ilogb(&self) -> Result<Exponent, Error> {
        if self.is_zero() {
            return Err(Error::InvalidArgument);
        }
        Ok(self.exponent().saturating_sub(1))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::defs::DEFAULT_P;
    use crate::ext::ExactNum;

    #[test]
    fn frexp_ldexp_roundtrip() {
        let rm = RoundingMode::ToEven;
        let x = ExactNum::from_word(6, DEFAULT_P);
        let (m, e) = ExactNum::from_word(6, DEFAULT_P).frexp();
        let m_exp = m.exponent();
        let y = m.ldexp(e, DEFAULT_P, rm);
        assert_eq!(x.compare(&y), Some(0));
        assert_eq!(
            ExactNum::from_word(6, DEFAULT_P)
                .scalb(1, DEFAULT_P, rm)
                .compare(&ExactNum::from_word(12, DEFAULT_P)),
            Some(0)
        );
        assert_eq!(x.ilogb(), Some(2));
        let lb = x.logb(DEFAULT_P, rm);
        assert_eq!(lb.compare(&ExactNum::from_word(2, DEFAULT_P)), Some(0));
        assert!(ExactNum::from_word(0, DEFAULT_P)
            .logb(DEFAULT_P, rm)
            .is_inf_neg());
        assert_eq!(m_exp, Some(0));
    }
}
