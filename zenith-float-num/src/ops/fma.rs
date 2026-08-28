//! Fused multiply-add: `a * b + c` with a single final rounding.

use crate::common::util::{bump_prec_retry, round_p};
use crate::defs::{Error, RoundingMode};
use crate::num::ExactNumNumber;
use crate::WORD_BIT_SIZE;

impl ExactNumNumber {
    /// Computes `self * b + c` with precision `p`, rounded once with `rm`.
    ///
    /// Unlike `mul` followed by `add`, the product is not rounded to `p` before the addition.
    ///
    /// ## Errors
    ///
    ///  - ExponentOverflow: the result exceeds the exponent range.
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - InvalidArgument: the precision is incorrect.
    pub fn fma(&self, b: &Self, c: &Self, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;

        if self.is_zero() || b.is_zero() {
            let mut ret = c.clone()?;
            ret.set_precision(p, rm)?;
            return Ok(ret);
        }

        if c.is_zero() {
            return self.mul(b, p, rm);
        }

        // When |a*b| and c are separated by more than p plus two words, the
        // smaller term cannot change the rounded p-bit result. Skip the
        // full-width product in that case (Horner / dot-product tails).
        let e_prod = (self.exponent() as i64).saturating_add(b.exponent() as i64);
        let e_c = c.exponent() as i64;
        let sep = e_prod.abs_diff(e_c);
        if sep > (p as u64).saturating_add((2 * WORD_BIT_SIZE) as u64) {
            if e_prod > e_c {
                return self.mul(b, p, rm);
            }
            let mut ret = c.clone()?;
            ret.set_precision(p, rm)?;
            ret.set_inexact(true);
            return Ok(ret);
        }

        let inexact = self.inexact() | b.inexact() | c.inexact();
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p + p_inc;

        // Exact product + add, then a single round to `p` (IEEE FMA). A
        // truncated product would drop sticky bits and fail the MPFR oracle.
        // Exact product + add, then a single round to `p` (IEEE FMA). A
        // truncated product would drop sticky bits and fail the MPFR oracle.
        let prod = self.mul_full_prec(b)?;
        let mut sum = prod.add_full_prec(c)?;

        loop {
            if sum.is_zero() {
                let mut ret = Self::new(p)?;
                ret.set_inexact(inexact | sum.inexact());
                return Ok(ret);
            }
            if sum.try_set_precision(p, rm, p_wrk)? {
                sum.set_inexact(sum.inexact() | inexact);
                return Ok(sum);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    /// Alias of [`Self::fma`].
    pub fn mul_add(&self, b: &Self, c: &Self, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        self.fma(b, c, p, rm)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::defs::WORD_BIT_SIZE;

    #[test]
    fn test_fma_zero_factor() {
        let rm = RoundingMode::ToEven;
        let p = WORD_BIT_SIZE * 4;
        let a = ExactNumNumber::from_i8(5, p).unwrap();
        let b = ExactNumNumber::new(p).unwrap();
        let c = ExactNumNumber::from_i8(3, p).unwrap();
        assert_eq!(a.fma(&b, &c, p, rm).unwrap().cmp(&c), 0);
    }

    #[test]
    fn test_fma_disjoint_exponents_skips_full_product() {
        let rm = RoundingMode::ToEven;
        let p = WORD_BIT_SIZE * 2;
        let a = ExactNumNumber::from_i8(3, p).unwrap();
        let b = ExactNumNumber::from_i8(5, p).unwrap();
        let mut c = ExactNumNumber::from_i8(1, p).unwrap();
        c.set_exponent(c.exponent() + 10_000);
        let got = a.fma(&b, &c, p, rm).unwrap();
        let mut want = c.clone().unwrap();
        want.set_precision(p, rm).unwrap();
        assert_eq!(got.cmp(&want), 0);
    }
}
