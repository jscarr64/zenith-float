//! General n-th root computation.

use crate::common::util::{bump_prec_retry, invert_rm_for_sign, round_p};
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
            let sign = if n.is_multiple_of(2) { Sign::Pos } else { self.sign() };
            return Self::new2(p, sign, self.inexact());
        }

        let mut sign = Sign::Pos;
        let mut x = self.clone()?;
        if self.is_negative() {
            if n.is_multiple_of(2) {
                return Err(Error::InvalidArgument);
            }
            sign = Sign::Neg;
            x = self.abs()?;
        }

        // One correctly rounded root. Chaining sqrt and cbrt rounded twice (x^(1/4) was
        // sqrt(sqrt(x)), each rounded at p bits) and could be one unit off.
        // The root of the magnitude is rounded in the direction `rm` has for the result's sign: `Up` on a negative
        // radicand rounds the magnitude down.
        let rm_abs = if sign == Sign::Neg { invert_rm_for_sign(rm) } else { rm };
        x = Self::nth_root_newton(&x, n, p, rm_abs)?;

        x.set_sign(sign);
        Ok(x)
    }

    /// Returns `cand` rounded to `p` bits if its `n`-th power equals `abs_self` exactly.
    fn nth_root_exact_candidate(
        abs_self: &Self,
        cand: &Self,
        n: usize,
        p: usize,
    ) -> Result<Option<Self>, Error> {
        let mut r = cand.clone()?;
        r.set_precision(p, RoundingMode::ToEven)?;
        // r = odd * 2^k with `t` significant bits: odd^n exceeds 2^((t - 1) * n), so r^n needs more than
        // (t - 1) * n significant bits and cannot equal `abs_self` once that reaches its significant bits.
        // This also bounds the power's precision below (a huge `n` used to ask for r^n at p * n bits).
        let t = Self::significant_bits(&r);
        if t > 1 && (t - 1).saturating_mul(n) >= Self::significant_bits(abs_self) {
            return Ok(None);
        }
        let bits = t.saturating_mul(n).saturating_add(WORD_BIT_SIZE);
        // A candidate rounded up past the largest number whose power fits (the root of a number near the top of
        // the range) has a power above the range: it is not the exact root.
        let pw = match r.powi(n, bits, RoundingMode::None) {
            Ok(pw) => pw,
            Err(Error::ExponentOverflow(_)) => return Ok(None),
            Err(e) => return Err(e),
        };
        if pw.cmp(abs_self) == 0 {
            r.set_inexact(abs_self.inexact());
            Ok(Some(r))
        } else {
            Ok(None)
        }
    }

    /// Number of bits from the leading one to the lowest one of the mantissa (0 for zero).
    fn significant_bits(v: &Self) -> usize {
        let mut tz = 0usize;
        for w in v.mantissa().digits() {
            if *w == 0 {
                tz += WORD_BIT_SIZE;
            } else {
                tz += w.trailing_zeros() as usize;
                break;
            }
        }
        v.mantissa_max_bit_len().saturating_sub(tz)
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

        // Error bound. Let p_x = p_wrk + 12 and u = 2^-p_x. The start x0 = exp(ln(a) / n) is formed at p_s bits:
        // |ln a| < 2^b with b = bitlen(|e|) + 1, so ln, 1/n and the product carry an absolute error below
        // 3 * 2^b * 2^-p_s / n + 2^-p_s, and exp adds one relative ulp: x0 = r * (1 + d0), |d0| < 2^(b + 3 - p_s).
        // One Newton step maps d0 to d1 = (n - 1) / 2 * d0^2 * (1 + O(d0)), below 2^(log2 n + 2(b + 3 - p_s)).
        // With p_s = p_x + b + bitlen(n) + 4 that is below 2^-(p_x + 6) = u / 64. The step itself
        // (powi < 2u, div, mul, add, div < u each) adds below 6u, so the total is below 7u < 2^12 * 2^-p_x,
        // which is the 2^-p_wrk relative error `try_set_precision(p, rm, p_wrk)` is told. Later iterations
        // start from the previous step at a higher p_wrk, with an even smaller d0.
        let e_abs = (abs_self.exponent() as i64).unsigned_abs();
        let b = (u64::BITS - e_abs.leading_zeros()) as usize + 1;
        let n_bits = (usize::BITS - n.leading_zeros()) as usize;
        let p_s = p_wrk + 12 + b + n_bits + 4;
        let ln = abs_self.ln(p_s, RoundingMode::None, &mut cc)?;
        let inv_n = n_num.reciprocal(p_s, RoundingMode::None)?;
        let mut x =
            ln.mul(&inv_n, p_s, RoundingMode::None)?
                .exp(p_s, RoundingMode::None, &mut cc)?;

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
                // The Newton value is exact only up to its working noise, so rounding it flags the result
                // inexact even when the root is exact (`nth_root(2^7, 7)` is 2): raise the rounded value to
                // the `n`-th power to see whether it is.
                if let Some(r) = Self::nth_root_exact_candidate(abs_self, &x_new, n, p)? {
                    break Ok(r);
                }
                x_new.set_inexact(x_new.inexact() | abs_self.inexact());
                break Ok(x_new);
            }
            // An exact root (e.g. 1e10^(1/5) = 100) can never be certified by the rounding test;
            // check the nearest p-bit candidate by exact powering instead.
            // `try_set_precision` modified `x_new`: the candidate is the same Newton step formed again from the
            // still-live operands.
            let cand =
                scaled
                    .add(&quot, p_x, RoundingMode::None)?
                    .div(&n_num, p_x, RoundingMode::None)?;
            if let Some(r) = Self::nth_root_exact_candidate(abs_self, &cand, n, p)? {
                break Ok(r);
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
    fn test_nth_root() -> Result<(), crate::Error> {
        let p = 320;
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;

        let sixteen = ExactNumNumber::from_word(16, p)?;
        assert!(
            sixteen
                .nth_root(4, p, rm)?
                .cmp(&ExactNumNumber::from_word(2, p)?)
                == 0
        );

        let neg = ExactNumNumber::from_i8(-8, p)?;
        let root = neg.nth_root(3, p, rm)?;
        assert!(root.cmp(&ExactNumNumber::from_i8(-2, p)?) == 0);

        assert!(neg.nth_root(4, p, rm).is_err());
        assert!(sixteen.nth_root(0, p, rm).is_err());

        let d1 = ExactNumNumber::parse(
            "F.FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF2DC85F7E77EC4872DC85F7E77EC487_e-1",
            crate::Radix::Hex,
            p,
            RoundingMode::None,
            &mut cc,
        )?;
        let sqrt = d1.sqrt(p, rm)?;
        assert!(d1.nth_root(2, p, rm)?.cmp(&sqrt) == 0);
        let cbrt = d1.cbrt(p, rm)?;
        assert!(d1.nth_root(3, p, rm)?.cmp(&cbrt) == 0);
        Ok(())
    }
}
