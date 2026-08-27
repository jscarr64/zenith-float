//! Reciprocal via Newton iteration.

use crate::common::consts::ONE;
use crate::common::util::{bump_prec_retry, round_p};
use crate::defs::{Error, WORD_BIT_SIZE};
use crate::num::ExactNumNumber;
use crate::{RoundingMode, Sign};

/// Schoolbook `1/x` stays cheaper than Newton setup for a couple of words.
const NEWTON_MIN_WORDS: usize = 3;

impl ExactNumNumber {
    /// Computes the reciprocal of a number with precision `p`. The result is rounded using the rounding mode `rm`.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - DivisionByZero: `self` is zero.
    ///  - ExponentOverflow: the resulting exponent becomes greater than the maximum allowed value for the exponent.
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    pub fn reciprocal(&self, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;

        if self.is_zero() {
            return Err(Error::DivisionByZero);
        }

        if p <= WORD_BIT_SIZE * (NEWTON_MIN_WORDS - 1) {
            return ONE.div(self, p, rm);
        }

        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p + p_inc;

        loop {
            let mut ret = self.reciprocal_newton(p_wrk)?;

            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }

            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    /// Newton: `x ← x * (2 − a x)`, doubling working precision each step.
    fn reciprocal_newton(&self, p: usize) -> Result<Self, Error> {
        let rm = RoundingMode::None;
        let p = round_p(p);

        let mut a = self.clone()?;
        a.set_sign(Sign::Pos);
        a.set_inexact(false);
        a.set_precision(p, rm)?;

        let two = Self::from_word(2, p)?;

        let mut p_cur = WORD_BIT_SIZE;
        let mut a_k = a.clone()?;
        a_k.set_precision(p_cur, rm)?;
        let mut x = ONE.div(&a_k, p_cur, rm)?;
        x.set_inexact(false);

        while p_cur < p {
            // Stay within quadratic reach: at most double. `round_p` on
            // `2*p_cur + ε` would jump a whole extra word and poison x.
            let p_next = p.min(p_cur.saturating_mul(2));

            x = Self::reciprocal_newton_step(&a, &two, x, p_next)?;
            p_cur = p_next;
        }

        // One more full-precision step mops up rounding from the last double.
        x = Self::reciprocal_newton_step(&a, &two, x, p)?;

        // Exact dyadics (powers of two) leave a zero residual; marking those
        // inexact makes try_set_precision loop forever (all sticky bits zero).
        let ax = a.mul(&x, p, rm)?;
        let one = Self::from_word(1, p)?;
        let exact = ax.cmp(&one) == 0;

        x.set_sign(self.sign());
        x.set_inexact(!exact);
        Ok(x)
    }

    fn reciprocal_newton_step(
        a: &Self,
        two: &Self,
        mut x: Self,
        p_next: usize,
    ) -> Result<Self, Error> {
        let rm = RoundingMode::None;
        x.set_precision(p_next, rm)?;
        let mut a_k = a.clone()?;
        a_k.set_precision(p_next, rm)?;
        a_k.set_inexact(false);

        let ax = a_k.mul(&x, p_next, rm)?;
        let t = two.sub(&ax, p_next, rm)?;
        let mut x = x.mul(&t, p_next, rm)?;
        x.set_inexact(false);
        Ok(x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::consts::ONE;
    use crate::defs::WORD_BIT_SIZE;
    use crate::{Consts, RoundingMode};

    #[test]
    fn test_reciprocal_newton_matches_div() {
        let mut cc = Consts::new().unwrap();
        let rm = RoundingMode::ToEven;

        let d1 = ExactNumNumber::parse(
            "F.FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF2DC85F7E77EC487_e-1",
            crate::Radix::Hex,
            320,
            RoundingMode::None,
            &mut cc,
        )
        .unwrap();
        let newton = d1.reciprocal(320, rm).unwrap();
        let div = ONE.div(&d1, 320, rm).unwrap();
        assert_eq!(newton.cmp(&div), 0);

        let p = WORD_BIT_SIZE * 8;
        for _ in 0..80 {
            let x = ExactNumNumber::random_normal(p, -64, 64).unwrap();
            if x.is_zero() {
                continue;
            }
            let r = x.reciprocal(p, rm).unwrap();
            let d = ONE.div(&x, p, rm).unwrap();
            assert_eq!(r.cmp(&d), 0);
        }

        let p = WORD_BIT_SIZE * 5;
        let t = ExactNumNumber::from_i8(3, p).unwrap();
        assert_eq!(
            t.reciprocal(p, rm)
                .unwrap()
                .cmp(&ONE.div(&t, p, rm).unwrap()),
            0
        );

        // Exact dyadic reciprocals must not NaN / error under ToEven.
        for v in [1i8, 2, 4, 8, 16, 32] {
            let x = ExactNumNumber::from_i8(v, p).unwrap();
            let r = x.reciprocal(p, rm).unwrap();
            let d = ONE.div(&x, p, rm).unwrap();
            assert_eq!(r.cmp(&d), 0, "1/{v}");
            assert!(!r.inexact(), "1/{v} should be exact");
        }

        let z = ExactNumNumber::new(p).unwrap();
        assert!(matches!(z.reciprocal(p, rm), Err(Error::DivisionByZero)));
    }
}
