//! General n-th root computation.

use crate::common::util::{bump_prec_retry, round_p};
use crate::ops::consts::Consts;
use crate::{
    defs::{Error, WORD_BIT_SIZE},
    num::ExactNumNumber,
    RoundingMode, Sign,
};

impl ExactNumNumber {
    /// Computes the `n`-th root of a number with precision `p`. The result is rounded using the rounding mode `rm`.
    /// `n = 2` and `n = 3` delegate to [`sqrt`](Self::sqrt) and [`cbrt`](Self::cbrt).
    ///
    /// ## Errors
    ///
    ///  - InvalidArgument: `n` is 0, or an even root of a negative number, or the precision is incorrect.
    ///  - MemoryAllocation: failed to allocate memory.
    pub fn nth_root(&self, n: usize, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;

        if n == 0 {
            return Err(Error::InvalidArgument);
        }

        if n == 1 {
            let mut ret = self.clone()?;
            ret.set_precision(p, rm)?;
            return Ok(ret);
        }

        if n == 2 {
            return self.sqrt(p, rm);
        }

        if n == 3 {
            return self.cbrt(p, rm);
        }

        if self.is_zero() {
            let sign = if n % 2 == 0 { Sign::Pos } else { self.sign() };
            return Self::new2(p, sign, self.inexact());
        }

        let mut sign = Sign::Pos;
        let mut x = self.clone()?;
        if self.is_negative() {
            if n % 2 == 0 {
                return Err(Error::InvalidArgument);
            }
            sign = Sign::Neg;
            x = self.abs()?;
        }

        let mut n_rem = n;
        while n_rem % 2 == 0 {
            x = x.sqrt(p, rm)?;
            n_rem /= 2;
        }
        while n_rem % 3 == 0 {
            x = x.cbrt(p, rm)?;
            n_rem /= 3;
        }

        if n_rem > 1 {
            x = Self::nth_root_newton(&x, n_rem, p, rm)?;
        }

        x.set_sign(sign);
        Ok(x)
    }

    fn nth_root_newton(
        abs_self: &Self,
        n: usize,
        p: usize,
        rm: RoundingMode,
    ) -> Result<Self, Error> {
        let mut cc = Consts::new()?;
        let n_num = Self::from_word(n as crate::Word, p)?;
        let nm1 = n - 1;

        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(abs_self.mantissa_max_bit_len());

        let ln = abs_self.ln(p_wrk + 12, RoundingMode::None, &mut cc)?;
        let inv_n = n_num.reciprocal(p_wrk + 12, RoundingMode::None)?;
        let mut x = ln.mul(&inv_n, p_wrk + 12, RoundingMode::None)?.exp(
            p_wrk + 12,
            RoundingMode::None,
            &mut cc,
        )?;

        loop {
            let p_x = p_wrk + 12;
            x.set_precision(p_x, RoundingMode::None)?;

            let x_nm1 = x.powi(nm1, p_x, RoundingMode::None)?;
            let quot = abs_self.div(&x_nm1, p_x, RoundingMode::None)?;
            let scaled = x.mul(
                &Self::from_word(nm1 as crate::Word, p_x)?,
                p_x,
                RoundingMode::None,
            )?;
            let mut x_new = scaled.add(&quot, p_x, RoundingMode::None)?;
            x_new = x_new.div(&n_num, p_x, RoundingMode::None)?;

            if x_new.try_set_precision(p, rm, p_wrk)? {
                x_new.set_inexact(x_new.inexact() | abs_self.inexact());
                break Ok(x_new);
            }

            x = x_new;
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::Consts;

    #[test]
    fn test_nth_root() {
        let p = 320;
        let mut cc = Consts::new().unwrap();
        let rm = RoundingMode::ToEven;

        let sixteen = ExactNumNumber::from_word(16, p).unwrap();
        assert!(
            sixteen
                .nth_root(4, p, rm)
                .unwrap()
                .cmp(&ExactNumNumber::from_word(2, p).unwrap())
                == 0
        );

        let neg = ExactNumNumber::from_i8(-8, p).unwrap();
        let root = neg.nth_root(3, p, rm).unwrap();
        assert!(root.cmp(&ExactNumNumber::from_i8(-2, p).unwrap()) == 0);

        assert!(neg.nth_root(4, p, rm).is_err());
        assert!(sixteen.nth_root(0, p, rm).is_err());

        let d1 = ExactNumNumber::parse(
            "F.FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF2DC85F7E77EC4872DC85F7E77EC487_e-1",
            crate::Radix::Hex,
            p,
            RoundingMode::None,
            &mut cc,
        )
        .unwrap();
        let sqrt = d1.sqrt(p, rm).unwrap();
        assert!(d1.nth_root(2, p, rm).unwrap().cmp(&sqrt) == 0);
        let cbrt = d1.cbrt(p, rm).unwrap();
        assert!(d1.nth_root(3, p, rm).unwrap().cmp(&cbrt) == 0);
    }
}
