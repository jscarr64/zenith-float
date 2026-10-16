//! Fused multiply-add: `a * b + c` with a single final rounding.

use crate::common::util::{bump_prec_retry, round_p};
use crate::defs::{Error, RoundingMode, EXPONENT_MIN};
use crate::num::ExactNumNumber;
use crate::Exponent;
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

        let inexact = self.inexact() | b.inexact() | c.inexact();
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p + p_inc;

        // Exact product + add, then a single round to `p` (IEEE FMA). A
        // truncated product would drop sticky bits and fail the MPFR oracle.
        let prod = self.mul_full_prec(b)?;

        // When one term lies entirely below the other's last bit by more than two words, only
        // its sign can affect the rounding. Replace it by a same-signed stand-in just below
        // that bit, so the full-width sum stays small (Horner / dot-product tails) and the
        // rounding and the inexact flag stay correct in every mode.
        // Exponents of the leading set bits: mantissas may carry leading zeros.
        let e_prod =
            prod.exponent() as i64 - (prod.mantissa_max_bit_len() as i64 - prod.precision() as i64);
        let e_c = c.exponent() as i64 - (c.mantissa_max_bit_len() as i64 - c.precision() as i64);
        // The stand-in must lie below the rounding position of the result as well, not only below the
        // other term's last bit: a term of fewer than `p` bits has its last bit above that position,
        // and a stand-in there would cut the sum (1 - tiny keeps only the bits down to the stand-in).
        let w_prod = (prod.mantissa_max_bit_len() as i64).max(p as i64);
        let w_c = (c.mantissa_max_bit_len() as i64).max(p as i64);
        let guard = 2 * WORD_BIT_SIZE as i64;
        // A stand-in whose exponent would fall below the smallest exponent cannot be stored (the exponent
        // would wrap around to a huge positive one): both terms are then within a mantissa's width of the
        // bottom of the range and the full-width sum is small, so it is formed without a stand-in.
        let min = EXPONENT_MIN as i64;
        let mut sum = if e_prod - e_c > w_prod + guard && e_prod - w_prod - 2 >= min {
            let mut tiny = c.clone()?;
            tiny.set_exponent((e_prod - w_prod - 2) as Exponent);
            prod.add_full_prec(&tiny)?
        } else if e_c - e_prod > w_c + guard && e_c - w_c - 2 >= min {
            let mut tiny = prod;
            tiny.set_exponent((e_c - w_c - 2) as Exponent);
            tiny.add_full_prec(c)?
        } else {
            prod.add_full_prec(c)?
        };

        loop {
            if sum.is_zero() {
                // An exact cancellation is +0, or -0 when rounding down.
                let mut ret = Self::new(p)?;
                if rm == RoundingMode::Down {
                    ret.set_sign(crate::Sign::Neg);
                }
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
    fn test_fma_zero_factor() -> Result<(), crate::Error> {
        let rm = RoundingMode::ToEven;
        let p = WORD_BIT_SIZE * 4;
        let a = ExactNumNumber::from_i8(5, p)?;
        let b = ExactNumNumber::new(p)?;
        let c = ExactNumNumber::from_i8(3, p)?;
        assert_eq!(a.fma(&b, &c, p, rm)?.cmp(&c), 0);
        Ok(())
    }

    #[test]
    fn test_fma_disjoint_exponents_skips_full_product() -> Result<(), crate::Error> {
        let rm = RoundingMode::ToEven;
        let p = WORD_BIT_SIZE * 2;
        let a = ExactNumNumber::from_i8(3, p)?;
        let b = ExactNumNumber::from_i8(5, p)?;
        let mut c = ExactNumNumber::from_i8(1, p)?;
        c.set_exponent(c.exponent() + 10_000);
        let got = a.fma(&b, &c, p, rm)?;
        let mut want = c;
        want.set_precision(p, rm)?;
        assert_eq!(got.cmp(&want), 0);
        Ok(())
    }
}
