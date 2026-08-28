//! IEEE-style exponent decomposition without hardware floats.

use crate::common::util::round_p;
use crate::defs::{Error, EXPONENT_MAX, EXPONENT_MIN};
use crate::num::ExactNumNumber;
use crate::{Exponent, RoundingMode, Sign};

impl ExactNumNumber {
    /// Split into a mantissa in `[0.5, 1)` and a power of two: `self = m · 2^e`.
    pub fn frexp(&self) -> Result<(Self, Exponent), Error> {
        if self.is_zero() {
            return Ok((self.clone()?, 0));
        }
        let e = self.exponent();
        let mut m = self.clone()?;
        m.set_exponent(0);
        Ok((m, e))
    }

    /// `self · 2^n` (same as `ldexp` / integer `scalb`).
    pub fn ldexp(&self, n: Exponent, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        if self.is_zero() {
            let mut z = self.clone()?;
            z.set_precision(p, rm)?;
            return Ok(z);
        }
        let e = self.exponent();
        match e.checked_add(n) {
            Some(sum) if Self::exponent_in_range(sum) => {
                let mut ret = self.clone()?;
                ret.set_exponent(sum);
                ret.set_precision(p, rm)?;
                Ok(ret)
            }
            _ => {
                if n >= 0 {
                    Err(Error::ExponentOverflow(self.sign()))
                } else {
                    Ok(Self::new2(p, self.sign(), true)?)
                }
            }
        }
    }

    #[cfg_attr(
        not(target_pointer_width = "32"),
        allow(clippy::absurd_extreme_comparisons)
    )]
    fn exponent_in_range(e: Exponent) -> bool {
        (EXPONENT_MIN..=EXPONENT_MAX).contains(&e)
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
        let (m, e) = x.frexp();
        let y = m.ldexp(e, DEFAULT_P, rm);
        assert_eq!(x.cmp(&y), Some(0));
        assert_eq!(
            x.scalb(1, DEFAULT_P, rm)
                .cmp(&ExactNum::from_word(12, DEFAULT_P)),
            Some(0)
        );
        assert_eq!(x.ilogb(), Some(2));
        let lb = x.logb(DEFAULT_P, rm);
        assert_eq!(lb.cmp(&ExactNum::from_word(2, DEFAULT_P)), Some(0));
        assert!(ExactNum::from_word(0, DEFAULT_P)
            .logb(DEFAULT_P, rm)
            .is_inf_neg());
        assert_eq!(m.exponent(), Some(0));
    }
}
