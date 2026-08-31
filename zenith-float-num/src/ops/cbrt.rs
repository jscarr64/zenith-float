//! Cube root computation.

use crate::{
    common::util::round_p,
    defs::{Error, EXPONENT_MIN},
    num::ExactNumNumber,
    Exponent, RoundingMode,
};

impl ExactNumNumber {
    /// Computes the cube root of a number with precision `p`. The result is rounded using the rounding mode `rm`.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory.
    ///  - InvalidArgument: the precision is incorrect.
    pub fn cbrt(&self, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;

        if self.is_zero() {
            return Self::new2(p, self.sign(), self.inexact());
        }

        let (e1, m1_opt) = self.normalize()?;
        let m1_normalized = m1_opt.as_ref().unwrap_or_else(|| self.mantissa());

        let exp_add = if e1 < 0 {
            -e1 % 3
        } else if e1 > 0 {
            3 - e1 % 3
        } else {
            0
        };

        let mut inexact = self.inexact();

        let (e_shift, m3) = m1_normalized.cbrt(p, rm, self.is_positive(), &mut inexact, exp_add)?;

        let e = (e1 + exp_add) / 3 + e_shift;

        if e < EXPONENT_MIN as isize {
            let mut ret =
                ExactNumNumber::from_raw_unchecked(m3, self.sign(), EXPONENT_MIN, inexact);

            ret.subnormalize(e, rm);

            Ok(ret)
        } else {
            Ok(ExactNumNumber::from_raw_unchecked(
                m3,
                self.sign(),
                e as Exponent,
                inexact,
            ))
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::{
        common::{consts::ONE, util::random_subnormal},
        defs::{EXPONENT_MAX, EXPONENT_MIN, WORD_BIT_SIZE},
        Consts, Exponent,
    };

    #[test]
    fn test_cbrt() {
        let mut cc = Consts::new().unwrap();
        let mut eps = ONE.clone().unwrap();

        for _ in 0..1000 {
            let prec = (crate::common::test_rng::random::<usize>() % 5 + 1) * WORD_BIT_SIZE;
            let d1 = ExactNumNumber::random_normal(prec, EXPONENT_MIN, EXPONENT_MAX).unwrap();
            let d2 = d1.cbrt(prec, RoundingMode::ToEven).unwrap();
            let d3 = d2
                .mul(&d2, prec, RoundingMode::ToEven)
                .unwrap()
                .mul(&d2, prec, RoundingMode::ToEven)
                .unwrap();

            eps.set_exponent(d1.exponent() - prec as Exponent + 3);

            //println!("d1 {:?}", d1.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            //println!("d2 {:?}", d2.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            //println!("d3 {:?}", d3.format(crate::Radix::Bin, RoundingMode::None).unwrap());

            assert!(
                d1.sub(&d3, prec, RoundingMode::ToEven)
                    .unwrap()
                    .abs()
                    .unwrap()
                    .cmp(&eps)
                    < 0
            );
        }

        // MAX
        let prec = 320;
        let d1 = ExactNumNumber::max_value(prec).unwrap();
        let d2 = d1.cbrt(prec, RoundingMode::ToEven).unwrap();
        let d3 = d2
            .mul(&d2, prec, RoundingMode::ToEven)
            .unwrap()
            .mul(&d2, prec, RoundingMode::ToEven)
            .unwrap();
        eps.set_exponent(d1.exponent() - prec as Exponent + 3);
        assert!(
            d1.sub(&d3, prec, RoundingMode::ToEven)
                .unwrap()
                .abs()
                .unwrap()
                .cmp(&eps)
                < 0
        );

        // MIN
        let d1 = ExactNumNumber::min_value(prec).unwrap();
        let d2 = d1.cbrt(prec, RoundingMode::ToEven).unwrap();
        let d3 = d2
            .mul(&d2, prec, RoundingMode::ToEven)
            .unwrap()
            .mul(&d2, prec, RoundingMode::ToEven)
            .unwrap();
        eps.set_exponent(d1.exponent() - prec as Exponent + 3);
        assert!(
            d1.sub(&d3, prec, RoundingMode::ToEven)
                .unwrap()
                .abs()
                .unwrap()
                .cmp(&eps)
                < 0
        );

        // MIN POSITIVE
        let d1 = ExactNumNumber::min_positive(prec).unwrap();
        let d2 = d1.cbrt(prec, RoundingMode::ToEven).unwrap();
        let d3 = d2
            .mul(&d2, prec, RoundingMode::ToEven)
            .unwrap()
            .mul(&d2, prec, RoundingMode::ToEven)
            .unwrap();
        let eps = ExactNumNumber::min_positive(prec).unwrap();
        assert!(
            d1.sub(&d3, prec, RoundingMode::ToEven)
                .unwrap()
                .abs()
                .unwrap()
                .cmp(&eps)
                <= 0
        );

        // near 1
        let p = 320;
        let d1 = ExactNumNumber::parse(
            "F.FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF2DC85F7E77EC4872DC85F7E77EC487_e-1",
            crate::Radix::Hex,
            p,
            RoundingMode::None,
            &mut cc,
        )
        .unwrap();
        let d2 = d1.cbrt(p, RoundingMode::ToEven).unwrap();
        let d3 = ExactNumNumber::parse(
            "F.FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFB9ED752A27F96D7B9ED752A27F96D8_e-1",
            crate::Radix::Hex,
            p,
            RoundingMode::None,
            &mut cc,
        )
        .unwrap();

        //println!("{:?}", d2.format(crate::Radix::Hex, RoundingMode::None).unwrap());

        assert!(d2.cmp(&d3) == 0);

        let d1 = ExactNumNumber::parse(
            "1.00000000000000000000000000000000000000000000000000000000000000002DC85F7E77EC487C",
            crate::Radix::Hex,
            p,
            RoundingMode::None,
            &mut cc,
        )
        .unwrap();
        let d2 = d1.cbrt(p, RoundingMode::ToEven).unwrap();
        let d3 = ExactNumNumber::parse(
            "1.00000000000000000000000000000000000000000000000000000000000000000F42CA7F7D4EC2D4",
            crate::Radix::Hex,
            p,
            RoundingMode::None,
            &mut cc,
        )
        .unwrap();

        // println!("{:?}", d2.format(crate::Radix::Hex, RoundingMode::None).unwrap());

        assert!(d2.cmp(&d3) == 0);

        // random subnormal arg
        for _ in 0..1000 {
            let d1 = random_subnormal(prec);

            let d2 = d1.cbrt(prec, RoundingMode::ToEven).unwrap();
            let d3 = d2
                .mul(&d2, prec, RoundingMode::ToEven)
                .unwrap()
                .mul(&d2, prec, RoundingMode::ToEven)
                .unwrap();

            let eps = ExactNumNumber::min_positive(prec).unwrap();
            // eps.set_exponent(eps.get_exponent() + 3);

            // println!("{:?}", d1.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("{:?}", d3.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("{:?}", eps.format(crate::Radix::Bin, RoundingMode::None).unwrap());

            assert!(
                d1.sub(&d3, prec, RoundingMode::ToEven)
                    .unwrap()
                    .abs()
                    .unwrap()
                    .cmp(&eps)
                    <= 0
            );
        }
    }
}
