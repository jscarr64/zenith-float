//! Hypotenuse: `sqrt(x² + y²)` without intermediate overflow.

use crate::common::util::{bump_prec_retry, round_p};
use crate::defs::{Error, WORD_BIT_SIZE};
use crate::num::ExactNumNumber;
use crate::{RoundingMode, Sign};

impl ExactNumNumber {
    /// Computes `sqrt(self² + other²)` with precision `p`, rounded with `rm`.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - ExponentOverflow: the result exceeds the exponent range.
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - InvalidArgument: the precision is incorrect.
    pub fn hypot(&self, other: &Self, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;

        let mut ax = self.abs()?;
        let mut ay = other.abs()?;
        if ax.abs_cmp(&ay) < 0 {
            core::mem::swap(&mut ax, &mut ay);
        }

        let inexact = self.inexact() | other.inexact();

        if ax.is_zero() {
            return Self::new2(p, Sign::Pos, inexact);
        }

        if ay.is_zero() {
            ax.set_precision(p, rm)?;
            ax.set_inexact(ax.inexact() | inexact);
            return Ok(ax);
        }

        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p + p_inc;

        loop {
            let p_x = p_wrk + 2;
            let r = ay.div(&ax, p_x, RoundingMode::None)?;
            let one = Self::from_word(1, p_x)?;
            let r2 = r.mul(&r, p_x, RoundingMode::None)?;
            let t = one.add(&r2, p_x, RoundingMode::None)?;
            let s = t.sqrt(p_x, RoundingMode::None)?;
            let mut ret = ax.mul(&s, p_x, RoundingMode::None)?;

            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | inexact);
                return Ok(ret);
            }

            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::consts::ONE;
    use crate::defs::WORD_BIT_SIZE;
    use crate::RoundingMode;

    #[test]
    fn test_hypot() {
        let rm = RoundingMode::ToEven;
        let p = WORD_BIT_SIZE * 5;

        let z = ExactNumNumber::new(p).unwrap();
        let h = z.hypot(&z, p, rm).unwrap();
        assert!(h.is_zero());

        let x = ExactNumNumber::from_i8(3, p).unwrap();
        let y = ExactNumNumber::from_i8(4, p).unwrap();
        let h = x.hypot(&y, p, rm).unwrap();
        let five = ExactNumNumber::from_i8(5, p).unwrap();
        assert_eq!(h.cmp(&five), 0);

        let nx = x.neg().unwrap();
        assert_eq!(nx.hypot(&y, p, rm).unwrap().cmp(&five), 0);

        for _ in 0..40 {
            let a = ExactNumNumber::random_normal(p, -10, 10).unwrap();
            let b = ExactNumNumber::random_normal(p, -10, 10).unwrap();
            let h = a.hypot(&b, p, rm).unwrap();
            let p2 = p + WORD_BIT_SIZE;
            let aa = a.mul(&a, p2, rm).unwrap();
            let bb = b.mul(&b, p2, rm).unwrap();
            let s = aa.add(&bb, p2, rm).unwrap().sqrt(p, rm).unwrap();
            let mut eps = ONE.clone().unwrap();
            eps.set_exponent(h.exponent() - p as crate::Exponent + 6);
            assert!(h.sub(&s, p, rm).unwrap().abs().unwrap().cmp(&eps) < 0);
        }
    }
}
