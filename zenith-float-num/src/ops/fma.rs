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

        let inexact = self.inexact() | b.inexact() | c.inexact();
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p + p_inc;

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
}
