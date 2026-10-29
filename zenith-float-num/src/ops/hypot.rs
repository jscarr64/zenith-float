//! Hypotenuse: `sqrt(x² + y²)` without intermediate overflow.
use crate::common::util::round_p_rm;
use crate::defs::{Error, Exponent, EXPONENT_MAX, EXPONENT_MIN, WORD_BIT_SIZE};
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
        let p = round_p_rm(p, rm);
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

        // sqrt(x^2 + y^2) with the sum of squares formed exactly after scaling x to
        // [0.5, 1): one correct rounding by sqrt, exact results included (3, 4 -> 5).
        // Mantissas may carry leading zeros (a subnormal operand has many); normalize so exponents
        // compare magnitudes. The exponents are kept in an `i64`, so a subnormal operand needs no
        // special case: the result is moved onto the subnormal grid below when it lies under the range.
        // ax and ay are not used again: move them.
        let ax_exp = ax.exponent();
        let mut sx = ax;
        let ex = ax_exp - sx.normalize2() as i64;
        let ay_exp = ay.exponent();
        let mut sy = ay;
        let ey = ay_exp - sy.normalize2() as i64;
        let wx = sx.mantissa_max_bit_len() as i64;
        sx.set_exponent(0);
        // y only adds a sticky bit once y^2 lies below both x^2's last bit and the
        // p-bit result's guard bits; past that it is replaced by a smaller stand-in.
        let lim = wx.max(p as i64 / 2 + 1) + WORD_BIT_SIZE as i64;
        let gap = crate::exp::wide(ex) - crate::exp::wide(ey);
        let dy = if gap > crate::exp::wide(lim + WORD_BIT_SIZE as i64) { -lim } else { ey - ex };
        sy.set_exponent(dy as Exponent);
        let s2 = sx
            .mul_full_prec(&sx)?
            .add_full_prec(&sy.mul_full_prec(&sy)?)?;
        let mut ret = s2.sqrt(p, rm)?;
        let e = crate::exp::narrow(crate::exp::wide(ret.exponent()) + crate::exp::wide(ex));
        if e > EXPONENT_MAX {
            return Err(Error::ExponentOverflow(Sign::Pos));
        }
        if e < EXPONENT_MIN {
            // Under the range (the operands were subnormal): the same step `sqrt`, `mul` and `div` take.
            ret.set_exponent(EXPONENT_MIN);
            ret.subnormalize(e, rm);
        } else {
            ret.set_exponent(e as Exponent);
        }
        ret.set_inexact(ret.inexact() | inexact);
        Ok(ret)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::defs::WORD_BIT_SIZE;
    use crate::RoundingMode;

    #[test]
    fn test_hypot() -> Result<(), crate::Error> {
        let rm = RoundingMode::ToEven;
        let p = WORD_BIT_SIZE * 5;

        let z = ExactNumNumber::new(p)?;
        let h = z.hypot(&z, p, rm)?;
        assert!(h.is_zero());

        let x = ExactNumNumber::from_i8(3, p)?;
        let y = ExactNumNumber::from_i8(4, p)?;
        let h = x.hypot(&y, p, rm)?;
        let five = ExactNumNumber::from_i8(5, p)?;
        assert_eq!(h.cmp(&five), 0);

        let nx = x.neg()?;
        assert_eq!(nx.hypot(&y, p, rm)?.cmp(&five), 0);

        for _ in 0..40 {
            let a = ExactNumNumber::random_normal(p, -10, 10)?;
            let b = ExactNumNumber::random_normal(p, -10, 10)?;
            let h = a.hypot(&b, p, rm)?;
            let p2 = p + WORD_BIT_SIZE;
            let aa = a.mul(&a, p2, rm)?;
            let bb = b.mul(&b, p2, rm)?;
            let s = aa.add(&bb, p2, rm)?.sqrt(p, rm)?;
            let mut eps = ExactNumNumber::small_const(1);
            eps.set_exponent(h.exponent() - p as crate::Exponent + 6);
            assert!(h.sub(&s, p, rm)?.abs()?.cmp(&eps) < 0);
        }
        Ok(())
    }
}
