// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! Cosine.

use crate::common::consts::ONE;
use crate::common::consts::TRIG_EXP_THRES;
use crate::common::consts::TWO;
use crate::common::util::bump_prec_retry;
use crate::common::util::calc_add_cost;
use crate::common::util::calc_mul_cost;
use crate::common::util::round_p_rm;
use crate::defs::Error;
use crate::defs::RoundingMode;
use crate::defs::EXPONENT_MIN;
use crate::defs::WORD_BIT_SIZE;
use crate::num::ExactNumNumber;
use crate::ops::consts::Consts;
use crate::ops::series::series_cost_optimize;
use crate::ops::series::series_prec;
use crate::ops::series::series_run;
use crate::ops::series::ArgReductionEstimator;
use crate::ops::series::FactPolycoeffGen;
use crate::ops::util::compute_small_exp;
use crate::Exponent;

struct CosArgReductionEstimator {}

impl ArgReductionEstimator for CosArgReductionEstimator {
    /// Estimates cost of reduction n times for number with precision p.
    fn reduction_cost(n: usize, p: usize) -> u64 {
        // n * (cost(add) + cost(mul))
        let cost_mul = calc_mul_cost(p);
        let cost_add = calc_add_cost(p);
        n as u64 * (cost_mul + cost_add) as u64
    }

    /// Given m, the negative power of 2 of a number, returns the negative power of 2 if reduction is applied n times.
    #[inline]
    fn reduction_effect(n: usize, m: i64) -> i64 {
        // n + m
        (n as i64) + m
    }
}

impl ExactNumNumber {
    /// Computes the cosine of a number with precision `p`. The result is rounded using the rounding mode `rm`.
    /// This function requires constants cache `cc` for computing the result.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory.
    ///  - InvalidArgument: the precision is incorrect.
    pub fn cos(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        let p = round_p_rm(p, rm);

        if self.is_zero() {
            let mut ret = Self::from_word(1, p)?;
            ret.set_inexact(self.inexact());
            return Ok(ret);
        }

        let mut p_inc = WORD_BIT_SIZE;
        // The retry budget counts from the working start: an argument wider than `p` must not exhaust it at once.
        let p_base = p.max(self.mantissa_max_bit_len());
        let mut p_wrk = p_base;

        compute_small_exp!(ONE, self.exponent() * 2 - 1, true, p_wrk, p, rm);

        p_wrk += p_inc;

        let mut add_p = (1 - TRIG_EXP_THRES) as usize;
        loop {
            let mut x = self.clone()?;

            let p_x = p_wrk + add_p;
            x.set_precision(p_x, RoundingMode::None)?;

            x = x.reduce_trig_arg(cc, RoundingMode::None)?;

            let (t, q) = x.trig_arg_pi_proximity(cc, RoundingMode::None)?;
            if q & 1 == 1 && add_p < t {
                add_p = t;
            } else {
                let mut ret = x.cos_series(RoundingMode::None)?;
                if ret.try_set_precision(p, rm, p_wrk)? {
                    ret.set_inexact(ret.inexact() | self.inexact());
                    break Ok(ret);
                }

                bump_prec_retry(&mut p_wrk, &mut p_inc, p_base)?;
            }
        }
    }

    /// cosine series
    pub(super) fn cos_series(mut self, rm: RoundingMode) -> Result<Self, Error> {
        // cos:  1 - x^2/2! + x^4/4! - x^6/6! + ...

        let p = self.mantissa_max_bit_len();
        let mut polycoeff_gen = FactPolycoeffGen::for_cos(p)?;
        let (reduction_times, niter, e_eff) = series_cost_optimize::<CosArgReductionEstimator>(
            p,
            &polycoeff_gen,
            -self.exponent(),
            2,
            false,
        );

        // Reduction gives 2^(-p+1) per iteration.
        // Restore gives 2^(-p+5) per iteration.
        // First terms of the series for any e_eff >= 0 give 2^(-p+6) at most.
        // The error of the remaining terms of the series is compensated.
        let add_prec = reduction_times as i64 * 6 + 6 - e_eff;
        let p_arg = series_prec(p, add_prec, 0)?;
        self.set_precision(p_arg, rm)?;

        let arg = if reduction_times > 0 { self.cos_arg_reduce(reduction_times)? } else { self };

        let acc = ExactNumNumber::from_word(1, p_arg)?; // 1
        let x_step = arg.mul(&arg, p_arg, rm)?; // x^2
        let x_first = arg.mul(&arg, p_arg, rm)?; // x^2 (the same product as `x_step`, formed again: one owner each)

        let ret = series_run(acc, x_first, x_step, niter, &mut polycoeff_gen)?;

        if reduction_times > 0 {
            ret.cos_arg_restore(reduction_times)
        } else {
            Ok(ret)
        }
    }

    // reduce argument n times.
    // cost: n * O(add)
    fn cos_arg_reduce(self, n: usize) -> Result<Self, Error> {
        // cos(2*x) = 2*cos(x)^2 - 1
        let mut ret = self;
        let p = ret.mantissa_max_bit_len();
        if ret.exponent() < EXPONENT_MIN + n as Exponent {
            ret.set_exponent(EXPONENT_MIN);
            for _ in 0..n - (ret.exponent() - EXPONENT_MIN) as usize {
                ret = ret.div(&TWO, p, RoundingMode::FromZero)?;
            }
        } else {
            ret.set_exponent(ret.exponent() - n as Exponent);
        }
        Ok(ret)
    }

    // restore value for the argument reduced n times.
    // cost: n * (4*O(mul) + O(add))
    fn cos_arg_restore(self, n: usize) -> Result<Self, Error> {
        // cos(2*x) = 2*cos(x)^2 - 1
        let mut cos = self;
        let p = cos.mantissa_max_bit_len();

        for _ in 0..n {
            let mut cos2 = cos.mul(&cos, p, RoundingMode::None)?;
            cos2.set_exponent(cos2.exponent() + 1);
            cos = cos2.sub(&ONE, p, RoundingMode::None)?;
        }

        Ok(cos)
    }
}

#[cfg(test)]
mod tests {

    use crate::common::util::random_subnormal;

    use super::*;

    #[test]
    fn test_cosine() -> Result<(), crate::Error> {
        let p = 320;
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;
        let mut n1 = ExactNumNumber::from_word(1, 320)?;
        n1.set_exponent(0);
        let _n2 = n1.cos(p, rm, &mut cc)?;
        //println!("{:?}", n2.format(crate::Radix::Dec, rm).unwrap());

        // asymptotic & extrema testing
        let mut half_pi = cc.pi_num(128, RoundingMode::None)?;
        half_pi.set_exponent(1);
        half_pi.set_precision(p, RoundingMode::None)?;

        let n2 = half_pi.cos(p, rm, &mut cc)?;
        let n3 = ExactNumNumber::parse("5.2049C1114CF98E804177D4C76273644A29410F31C6809BBDF2A33679A7486365EEEE1A43A7D13E58_e-21", crate::Radix::Hex, 640, RoundingMode::None, &mut cc)?;

        //println!("{:?}", n2.format(crate::Radix::Bin, rm).unwrap());

        assert!(n2.cmp(&n3) == 0);

        // large exponent
        half_pi.set_exponent(256);
        let n2 = half_pi.cos(p, rm, &mut cc)?;
        let n3 = ExactNumNumber::parse("3.2F00069261A9FFC022D09F662F2E2DDBEFD1AF138813F2A71D7601C58E793299EA052E4EBC59107C_e-1", crate::Radix::Hex, 640, RoundingMode::None, &mut cc)?;

        //println!("{:?}", n2.format(crate::Radix::Hex, rm).unwrap());

        assert!(n2.cmp(&n3) == 0);

        // small exponent
        let p = 384;
        let n1 = ExactNumNumber::parse("1.992EF09686C3DC782C05BFD7863A715ECBDAED45DBAEE3ADFEF1AB8F74DB393D8CD6EAF9CA8443A6C28CF59D35B7FF56_e-20", crate::Radix::Hex, p, RoundingMode::None, &mut cc)?;
        let n2 = n1.cos(p, RoundingMode::ToEven, &mut cc)?;
        let n3 = ExactNumNumber::parse("F.FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEB8FC7D51D69792F9AB263F754D161F6A_e-1", crate::Radix::Hex, p, RoundingMode::None, &mut cc)?;

        // println!("{:?}", n1.format(crate::Radix::Hex, rm).unwrap());

        assert!(n2.cmp(&n3) == 0);

        let d3 = ExactNumNumber::min_positive(p)?;
        let zero = ExactNumNumber::new(1)?;
        let d4 = random_subnormal(p)?;

        assert!(d3.cos(p, rm, &mut cc)?.cmp(&ONE) == 0);
        assert!(zero.cos(p, rm, &mut cc)?.cmp(&ONE) == 0);
        assert!(d4.cos(p, rm, &mut cc)?.cmp(&ONE) == 0);
        Ok(())
    }

    /// `cos` of a 9984-bit argument near 100 at 256 bits toward `-Inf` (an `mpfr_compare_special` case):
    /// the result is a plain normal number, never NaN.
    #[test]
    #[cfg(target_pointer_width = "64")]
    fn cos_wide_argument_converges() -> Result<(), crate::Error> {
        use crate::{defs::Word, Sign};
        let mut cc = Consts::new()?;
        let words: [Word; 156] = [
            14952479658564103063,
            12531007309131022797,
            6494135004533635233,
            5703182228362271784,
            9392849452679478258,
            11016009685198680855,
            5996340738141207241,
            6001976938276149104,
            12138510522634603581,
            6549352291948679157,
            10364232064643934041,
            7500257890532724367,
            5065846620764132106,
            15391355858284472248,
            13486771376285451775,
            7438706533472103577,
            16959335378302138335,
            3953388940761123208,
            10583067949125126737,
            4768909473997734362,
            1507592445337527272,
            13820704256080041369,
            7897300851900710509,
            18205162792456255045,
            157234480653204039,
            781419742571896064,
            8467846813433917881,
            9715460753583295045,
            7841272119379176018,
            4932011525869982352,
            6974247187197031402,
            14420340088455655025,
            11398699590452942047,
            10761602304646310314,
            14490628286885031495,
            9843180136729294891,
            4614293558707121367,
            3859662898394013586,
            18082915414602392168,
            12814166720816565619,
            17990552410248806049,
            7274189455706457373,
            5014732456614780353,
            16082652101509605942,
            7866121656090294317,
            444516014204549286,
            8748298368127327546,
            18058673092389013480,
            5461679451765179869,
            3779642020548194411,
            1435337285752320199,
            2459578943964150791,
            16428490146137279193,
            2562291777950016813,
            3757044925203597688,
            4227582198408344936,
            2699616139681735782,
            8367591073790052333,
            4307748803221246720,
            15618136957341217210,
            4068857912317108108,
            8404722335361197540,
            16615178834362110297,
            1350848282291893462,
            305839494018916590,
            4071228600131257448,
            15732984364502521394,
            16940796536648095068,
            17676815718716260336,
            477897831538350710,
            12076145041833976616,
            1498677032303736800,
            14740120496549454296,
            5464000709718529194,
            1089103341458840261,
            2569719744895420606,
            6478368139494549747,
            13089782898566569118,
            15374849526239595650,
            17618881973843197686,
            17308153697367366915,
            6450221200295301868,
            8050452619988882111,
            3684033904258821954,
            15816382715395843403,
            4958019789371259675,
            15726248754671195373,
            1699871023477459224,
            12712065177997531185,
            16562893830193414574,
            17942222196870894350,
            14516542467089250237,
            13967031003058161374,
            3946273353039114278,
            5617173773515130002,
            9727695878091096207,
            10643389760314269121,
            14989828010024131241,
            13266440626786813980,
            2261134052400748463,
            2538800413216073709,
            11075079702397895586,
            2917304821422473309,
            3354895935799741094,
            15833463160886431481,
            2703261255863651560,
            1915270317141724122,
            11079742494088684092,
            9831808858219764518,
            1905651643148920535,
            5441141238023653377,
            4889602962105372814,
            640161647614602033,
            8577537326653720290,
            13538127619335216492,
            5917495817861275660,
            15597657053442108004,
            17379384394372483172,
            14907990252627095147,
            2201298008914878877,
            12391001869773255639,
            12940829120977298631,
            10009937741821119613,
            17076388665822539530,
            12129638193458625124,
            12480344337048173107,
            1545454141666149421,
            1572361106993317136,
            4149303432552213221,
            16009113560346531608,
            13097978378517762761,
            11180049335738409615,
            16401677924195796483,
            3643515754058796603,
            17398650045445185916,
            7425368496004700164,
            11445099139962410605,
            2045464732017971899,
            9468076200223288726,
            7572309818504171359,
            11014195235928789914,
            13979310375781515013,
            5271575865889938237,
            12582815541414797286,
            17165588707022577573,
            864511594326308845,
            17603518614767922539,
            16466767132611215046,
            5755940542857986629,
            3470879405153129527,
            17263733006627652379,
            5857503583518590173,
            147421033984662306,
            2955010104097229940,
            14179128828124470481,
            14488038916154245684,
        ];
        let x = ExactNumNumber::from_words(&words, Sign::Pos, 7)?;
        let r = x.cos(256, RoundingMode::Down, &mut cc)?;
        assert!(r.exponent() <= 0 && !r.is_zero());
        Ok(())
    }
}
