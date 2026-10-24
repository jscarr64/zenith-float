// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! Hyperbolic arccosine.

use crate::common::consts::ONE;
use crate::common::util::bump_prec_retry;
use crate::common::util::count_leading_zeroes_skip_first;
use crate::common::util::round_p;
use crate::defs::Error;
use crate::defs::RoundingMode;
use crate::num::ExactNumNumber;
use crate::Consts;
use crate::Sign;
use crate::EXPONENT_MAX;
use crate::WORD_BIT_SIZE;

impl ExactNumNumber {
    /// Computes the hyperbolic arccosine of a number with precision `p`. The result is rounded using the rounding mode `rm`.
    /// This function requires constants cache `cc` for computing the result.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - ExponentOverflow: the result is too large or too small number.
    ///  - MemoryAllocation: failed to allocate memory.
    ///  - InvalidArgument: when `self` < 1, or the precision is incorrect.
    pub fn acosh(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        let p = round_p(p);

        let cmpone = self.cmp(&ONE);
        if cmpone == 0 {
            return Self::new2(p, Sign::Pos, self.inexact());
        } else if cmpone < 0 {
            return Err(Error::InvalidArgument);
        }

        let p_max = p.max(self.mantissa_max_bit_len());

        let mut x = self.clone()?;
        x.set_inexact(false);

        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p_max + p_inc;

        loop {
            let mut p_x = p_wrk;

            let mut ret = if (self.exponent() - 1) / 2 > p_x as i64 + 2 {
                // acosh(x) = ln(2*x)
                if self.exponent() == EXPONENT_MAX {
                    // ln(2) + ln(x)
                    p_x += 4; // ln2, ln x and their sum: below 3 ulps of p_x

                    let lnx = x.ln(p_x, RoundingMode::None, cc)?;

                    let ln2 = cc.ln_2_num(p_x, RoundingMode::None)?;

                    ln2.add(&lnx, p_x, RoundingMode::None)?
                } else {
                    // ln at p_x is within one p_x ulp and ln(2x) - acosh(x) = 1/(4x^2) + ... is far below it;
                    // 3 guard bits keep the sum under the 2^-p_wrk the rounding test is told (it was 0).
                    p_x += 3;
                    x.set_exponent(self.exponent() + 1);

                    let lnx = x.ln(p_x, RoundingMode::None, cc)?;

                    x.set_exponent(self.exponent());

                    lnx
                }
            } else {
                // ln(x + sqrt(x*x - 1))
                let mut additional_prec = 0;
                if self.exponent() == 1 {
                    additional_prec = count_leading_zeroes_skip_first(self.mantissa().digits());
                }

                p_x += 6 + additional_prec;

                let xx = x.mul(&x, p_x, RoundingMode::None)?;

                let d1 = xx.sub(&ONE, p_x, RoundingMode::None)?;

                let d2 = d1.sqrt(p_x, RoundingMode::None)?;

                let d3 = d2.add(&x, p_x, RoundingMode::FromZero)?;

                d3.ln(p_x, RoundingMode::None, cc)?
            };

            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }

            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::Exponent;

    use super::*;

    #[test]
    fn test_acosh() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;

        // near 1
        let p = 448;
        let n1 = ExactNumNumber::parse("1.0000000000000000000000000000000000000000000B56A0EBA6F7D47E21A7B2A7806A698BABAF2F05BC61E2F8FB50FE0B98F55B181AC9C8_e+0", crate::Radix::Hex, p, RoundingMode::None, &mut cc)?;
        let n2 = n1.acosh(p, rm, &mut cc)?;
        let n3 = ExactNumNumber::parse("4.C31368910963B1A1BCFC0EDBD393FB7A5E876F9751D93A20E7E48EC0D16090ADA5F46DF2184D32A19C500088EA09CBD4F23DF713113D8A58_e-16", crate::Radix::Hex, p, RoundingMode::None, &mut cc)?;

        // println!("{:?}", n1.format(crate::Radix::Bin, rm).unwrap());
        // println!("{:?}", n2.format(crate::Radix::Hex, rm).unwrap());

        assert!(n2.cmp(&n3) == 0);

        // large exp
        let p = 320;
        let n1 = ExactNumNumber::parse("1.921FB54442D18469898CC51701B839A200000000000000004D3C337F7C8D419EBBFC39B4BEC14AF6_e+1000", crate::Radix::Hex, p, RoundingMode::None, &mut cc)?;
        let n2 = n1.acosh(p, rm, &mut cc)?;
        let n3 = ExactNumNumber::parse("2.C5DAB0AF9025886C3364C7B6D6741EB19D4FB009D3F92CA21B77498D9F0666363C665F2F324EAEC8_e+3", crate::Radix::Hex, p, RoundingMode::None, &mut cc)?;

        //println!("{:?}", n2.format(crate::Radix::Hex, rm).unwrap());

        assert!(n2.cmp(&n3) == 0);

        let d1 = ExactNumNumber::max_value(p)?;
        let d2 = ExactNumNumber::min_value(p)?;
        let d3 = ExactNumNumber::min_positive(p)?;

        // MAX
        let mut d4 = d1.acosh(p, RoundingMode::Down, &mut cc)?;

        d4.set_exponent(d4.exponent() - 1);
        let d5 = d4.cosh(p, rm, &mut cc)?;
        let mut d5 = d5.mul(&d5, p, rm)?;
        d5.set_exponent(d5.exponent() + 1);

        let mut eps = ExactNumNumber::small_const(1);
        eps.set_exponent(
            d1.exponent() - p as Exponent + core::mem::size_of::<Exponent>() as Exponent * 8,
        );

        assert!(d1.sub(&d5, p, RoundingMode::ToEven)?.abs()?.cmp(&eps) < 0);

        assert!(matches!(
            d2.acosh(p, rm, &mut cc),
            Err(Error::InvalidArgument)
        ));
        assert!(matches!(
            d3.acosh(p, rm, &mut cc),
            Err(Error::InvalidArgument)
        ));

        assert!(ONE.acosh(p, rm, &mut cc)?.is_zero());
        Ok(())
    }
}
