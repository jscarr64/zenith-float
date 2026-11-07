// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! Cosine.

use crate::common::consts::ONE;
use crate::common::consts::TRIG_EXP_THRES;
use crate::common::consts::TWO;
use crate::common::util::bump_prec_retry;
use crate::common::util::calc_add_cost;
use crate::common::util::calc_mul_cost;
use crate::common::util::cost_from_u128;
use crate::common::util::cost_product;
use crate::common::util::cost_sum;
use crate::common::util::round_p;
use crate::common::util::wide;
use crate::defs::Error;
use crate::defs::RoundingMode;
use crate::defs::EXPONENT_MIN;
use crate::defs::WORD_BIT_SIZE;
use crate::num::ExactNumNumber;
use crate::ops::consts::Consts;
use crate::ops::series::series_cost_optimize;
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
        let per_reduction =
            cost_sum(&[Some(u128::from(calc_mul_cost(p))), Some(u128::from(calc_add_cost(p)))]);
        cost_from_u128(cost_product(Some(wide(n)), per_reduction))
    }

    /// Given m, the negative power of 2 of a number, returns the negative power of 2 if reduction is applied n times.
    #[inline]
    fn reduction_effect(n: usize, m: i64) -> usize {
        // n + m
        ((n as i64) + m) as usize
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
        let p = round_p(p);

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
        let add_prec = reduction_times as i64 * 6 + 6 - e_eff as i64;
        let p_arg = p + if add_prec > 0 { add_prec as usize } else { 0 };
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

#[cfg(test)]
mod cost_overflow_tests {
    use crate::ops::series::ArgReductionEstimator;

    use super::CosArgReductionEstimator;

    /// Unclamped `f(p)`. The `u64` clamp is applied to `n * (f(p) + p)`, not to `f`.
    fn f_mul(p: u128) -> u128 {
        if p < 70 {
            p * p
        } else if p < 1625 {
            (p * p * p).isqrt()
        } else {
            let q = p.isqrt();
            q * q * q
        }
    }

    #[test]
    fn reduction_cost_matches_integer_oracle() {
        // 1: one reduction; the product is f(p) + p before the clamp.
        // 15: fifteen reductions; the product is 15 * (f(p) + p) before the clamp.
        let ns = [1_usize, 15_usize];
        // 1625: first p of the isqrt(p)^3 branch; 40^2 = 1_600 <= 1625 < 1_681 = 41^2; <= MAX_PRECISION on every target.
        // 2_643_876 = 1626^2: isqrt = 1626, 1626^3 = 4_298_942_376 > u32::MAX; still <= 2^30 - 2.
        // MAX_PRECISION: 2^30 - 2 on 32-bit (isqrt = 32_767, 32_767^3 > u32::MAX) and 2^60 - 1 on 64-bit (isqrt = 2^30 - 1, cube > u64::MAX).
        let ps = [1625_usize, 2_643_876_usize, crate::defs::MAX_PRECISION];
        for n in ns {
            for p in ps {
                let n_u = crate::common::util::wide(n);
                let p_u = crate::common::util::wide(p);
                let oracle = (n_u * (f_mul(p_u) + p_u)).min(u128::from(u64::MAX));
                let result = CosArgReductionEstimator::reduction_cost(n, p);
                assert_eq!(u128::from(result), oracle, "n = {n}, p = {p}");
            }
        }
    }

    /// Fixed iter cost. `next` is never called by `series_cost_optimize`.
    struct FixedCostGen(u64);

    impl crate::ops::series::PolycoeffGen for FixedCostGen {
        fn next(
            &mut self,
            _rm: crate::RoundingMode,
        ) -> Result<&crate::num::ExactNumNumber, crate::defs::Error> {
            Err(crate::defs::Error::InvalidArgument)
        }

        fn iter_cost(&self) -> u64 {
            self.0
        }
    }

    #[test]
    fn series_cost_optimize_terminates_at_max_precision() {
        // On 64-bit the cost plateau (all costs clamped at u64::MAX) must stop the loop via `cost2 < cost1` being false.
        let p = crate::defs::MAX_PRECISION;
        // 2 * MAX_PRECISION: add cost of two numbers. 2*(2^60 - 1) and 2*(2^30 - 2) both fit in usize.
        // 2 * p <= 2^61 always fits u64 on every target.
        let gen = FixedCostGen(u64::try_from(2 * p).unwrap_or(u64::MAX));
        // m = 0: negative exponent is zero, so the first candidate reduction count is floor(log2(p))/2.
        // pwr_step = 2: cosine powers advance by x^2, the same arguments `cos_series` passes.
        // ext = false: the series part is `series_cost`.
        let (reduction_times, niter, _) = crate::ops::series::series_cost_optimize::<
            CosArgReductionEstimator,
        >(p, &gen, 0, 2, false);
        // floor(log2(p))/2 is the optimizer step; a finished count stays <= p plus that step.
        let step = crate::common::util::log2_floor(p) / 2;
        assert!(
            reduction_times <= p + step,
            "reduction_times {reduction_times} p {p}"
        );
        // niter > 0: the plateau stop happens while `reduction_effect` is still below p. The guard at p + step returns niter == 0.
        assert!(niter > 0, "niter {niter} p {p}");
    }

    /// The returned triple must be a candidate the loop evaluated.
    ///
    /// `initial` is the first reduction count: `step - m` when `step > m`, else 0.
    /// Evaluated counts are `initial + k * step`. `m_eff` is that count plus `m`,
    /// and `niter` is the series length at that `m_eff`.
    fn assert_evaluated_candidate(p: usize, m: i64, ext: bool, gen: &FixedCostGen) {
        let step = crate::common::util::log2_floor(p) / 2;
        assert!(step > 0, "reduction_num_step p {p}");
        let (rt, niter, m_eff) =
            crate::ops::series::series_cost_optimize::<CosArgReductionEstimator>(p, gen, m, 2, ext);
        // step fits in i64 for every accepted precision (log2_floor(p) / 2 <= 30).
        let step_i = i64::try_from(step).unwrap_or(i64::MAX);
        let initial =
            if step_i > m { usize::try_from(step_i - m).unwrap_or(usize::MAX) } else { 0 };
        assert!(
            rt >= initial,
            "rt {rt} initial {initial} p {p} m {m} ext {ext}"
        );
        let advanced = rt.checked_sub(initial);
        assert!(
            advanced.is_some(),
            "rt {rt} initial {initial} p {p} m {m} ext {ext}"
        );
        assert!(
            advanced.is_some_and(|delta| delta.is_multiple_of(step)),
            "rt {rt} initial {initial} step {step} p {p} m {m} ext {ext}"
        );
        let m_usize = usize::try_from(m).unwrap_or(usize::MAX);
        assert_eq!(
            rt.checked_add(m_usize),
            Some(m_eff),
            "m_eff {m_eff} rt {rt} m {m} p {p} ext {ext}"
        );
        let ln = crate::common::util::log2_floor(p);
        let lln = crate::common::util::log2_floor(ln);
        let denom = ln
            .checked_sub(lln)
            .and_then(|v| v.checked_add(m_eff))
            .and_then(|v| v.checked_sub(2));
        assert!(
            denom.is_some_and(|d| d > 0),
            "denom ln {ln} lln {lln} m_eff {m_eff} p {p}"
        );
        assert_eq!(
            denom.filter(|d| *d > 0).map(|d| p / d / 2),
            Some(niter),
            "niter {niter} m_eff {m_eff} p {p} m {m} ext {ext}"
        );
        assert!(
            p.checked_add(step).is_some_and(|limit| rt <= limit),
            "rt {rt} p {p} step {step} m {m} ext {ext}"
        );
    }

    #[test]
    fn series_cost_optimize_first_candidate_capped_does_not_underflow() {
        // p = 1 << 28 = 268_435_456, accepted on 32- and 64-bit. log2_floor(p) = 28, so reduction_num_step = 14.
        // For m = 14 the initial reduction_times is 0, m_eff = 14,
        // niter = 268_435_456 / (28 - 4 + 14 - 2) / 2 = 3_728_270,
        // calc_mul_cost(p) = (2^14)^3 = 2^42,
        // exact series cost 19_688_316_708_440_768_512 > u64::MAX, so the first candidate is capped.
        let p = 1_usize << 28;
        assert_eq!(p, 268_435_456);
        assert_eq!(crate::common::util::log2_floor(p) / 2, 14);
        // 2 * p = 2^29: add cost of two numbers. It fits in u64.
        let gen = FixedCostGen(u64::try_from(2 * p).unwrap_or(u64::MAX));
        // m = 1: initial = 13, m_eff = 14, niter = 3_728_270, series cost 19_688_316_708_440_768_512 > u64::MAX.
        assert_evaluated_candidate(p, 1, false, &gen);
        // m = 13: initial = 1, m_eff = 14, same niter 3_728_270 and the same series cost above u64::MAX.
        assert_evaluated_candidate(p, 13, false, &gen);
        // m = 14: initial = 0, m_eff = 14, niter = 3_728_270, series cost 19_688_316_708_440_768_512 > u64::MAX.
        assert_evaluated_candidate(p, 14, false, &gen);
        // m = 15: initial = 0, m_eff = 15, niter = 268_435_456 / 37 / 2 = 3_627_506, series cost 19_156_310_110_238_670_848 > u64::MAX.
        assert_evaluated_candidate(p, 15, false, &gen);
        // ext: series part is iter_cost * niter with iter_cost = u64::MAX, capped whenever niter > 1.
        let gen_max = FixedCostGen(u64::MAX);
        // m = 1, ext: m_eff = 14, niter = 3_728_270, and u64::MAX * 3_728_270 does not fit in u64.
        assert_evaluated_candidate(p, 1, true, &gen_max);
        // m = 13, ext: m_eff = 14, niter = 3_728_270, and u64::MAX * 3_728_270 does not fit in u64.
        assert_evaluated_candidate(p, 13, true, &gen_max);
        // m = 14, ext: m_eff = 14, niter = 3_728_270, and u64::MAX * 3_728_270 does not fit in u64.
        assert_evaluated_candidate(p, 14, true, &gen_max);
        // m = 15, ext: m_eff = 15, niter = 3_627_506, and u64::MAX * 3_627_506 does not fit in u64.
        assert_evaluated_candidate(p, 15, true, &gen_max);

        let p_max = crate::defs::MAX_PRECISION;
        // 2 * MAX_PRECISION fits u64: 2*(2^30 - 2) on 32-bit and 2*(2^60 - 1) on 64-bit.
        let gen_max_p = FixedCostGen(u64::try_from(2 * p_max).unwrap_or(u64::MAX));
        // m = 0: initial = step. 32-bit step = 14, m_eff = 14, niter = (2^30 - 2) / 37 / 2 = 14_510_024, and 14_510_024 * 32_767^3 = 510_479_344_801_353_209_912 > u64::MAX. 64-bit step = 29, m_eff = 29, denom = 81, niter = (2^60 - 1) / 81 / 2 > 1, and (2^30 - 1)^3 > u64::MAX so calc_mul_cost is already capped.
        assert_evaluated_candidate(p_max, 0, false, &gen_max_p);
        // m = 1: initial = step - 1, m_eff = step, the series cost matches the capped m = 0 cost on both widths.
        assert_evaluated_candidate(p_max, 1, false, &gen_max_p);
        // m = 0, ext: niter > 1 on both widths (14_510_024 on 32-bit; (2^60 - 1) / 81 / 2 on 64-bit), so u64::MAX * niter does not fit in u64.
        assert_evaluated_candidate(p_max, 0, true, &gen_max);
        // m = 1, ext: m_eff and niter match m = 0, so u64::MAX * niter does not fit in u64.
        assert_evaluated_candidate(p_max, 1, true, &gen_max);
    }

    fn widen_usize(n: usize) -> u128 {
        match u128::try_from(n) {
            Ok(v) => v,
            Err(_) => u128::MAX,
        }
    }

    fn clamp_to_u64(v: u128) -> u128 {
        let cap = u128::from(u64::MAX);
        if v > cap {
            cap
        } else {
            v
        }
    }

    /// Schoolbook `p^2`, then `isqrt(p^3)`, then `isqrt(p)^3`, then the `u64` clamp.
    /// Matches `calc_mul_cost` without calling it.
    fn mul_cost(p: usize) -> u128 {
        let p_u = widen_usize(p);
        let raw = if p_u < 70 {
            p_u.checked_mul(p_u)
        } else if p_u < 1625 {
            p_u.checked_mul(p_u)
                .and_then(|s| s.checked_mul(p_u))
                .map(|v| v.isqrt())
        } else {
            let q = p_u.isqrt();
            q.checked_mul(q).and_then(|s| s.checked_mul(q))
        };
        match raw {
            Some(v) => clamp_to_u64(v),
            None => u128::from(u64::MAX),
        }
    }

    fn add_cost(p: usize) -> u128 {
        clamp_to_u64(widen_usize(p))
    }

    /// `n * per`, with `None` and `u128` overflow both becoming `u64::MAX`. `n = 0` and `per = Some` is 0.
    fn scale_cost(n: usize, per: Option<u128>) -> u64 {
        let prod = per.and_then(|v| widen_usize(n).checked_mul(v));
        match prod {
            Some(v) => match u64::try_from(v) {
                Ok(c) => c,
                Err(_) => u64::MAX,
            },
            None => u64::MAX,
        }
    }

    /// Cosine reduction: `n * (mul(p) + add(p))`.
    fn cos_per(p: usize) -> Option<u128> {
        mul_cost(p).checked_add(add_cost(p))
    }

    fn assert_cap(p: usize, per: Option<u128>, cost: impl Fn(usize) -> u64) {
        let cap = u128::from(u64::MAX);
        let Some(per) = per else {
            assert_eq!(cost(0), u64::MAX, "none per p={p}");
            assert_eq!(cost(1), u64::MAX, "none per p={p}");
            return;
        };
        if per == 0 {
            assert_eq!(cost(0), 0, "p={p}");
            assert_eq!(cost(1), 0, "p={p}");
            return;
        }
        if per > cap {
            // One step already exceeds `u64`. Zero reductions stay exact.
            assert_eq!(cost(0), 0, "p={p}");
            assert_eq!(cost(1), u64::MAX, "per above u64 p={p}");
            return;
        }
        let n_below = cap / per;
        let exact = match n_below.checked_mul(per) {
            Some(v) => v,
            None => u128::MAX,
        };
        assert!(exact <= cap, "p={p}");
        let crossed = match n_below.checked_add(1).and_then(|n| n.checked_mul(per)) {
            Some(v) => v > cap,
            None => true,
        };
        assert!(crossed, "successor crosses u64 p={p}");
        match usize::try_from(n_below) {
            Ok(n) => {
                assert_eq!(u128::from(cost(n)), exact, "below cap p={p} n={n}");
                match n.checked_add(1) {
                    Some(n_hi) => assert_eq!(cost(n_hi), u64::MAX, "at cap p={p} n={n_hi}"),
                    None => assert!(crossed, "successor is not a usize p={p}"),
                }
            }
            Err(_) => {
                // The first capped count does not fit in `usize`, so the cap is unreachable.
                match widen_usize(usize::MAX).checked_mul(per) {
                    Some(v) => {
                        assert!(v <= cap, "unreachable cap still below u64 p={p}");
                        assert_eq!(u128::from(cost(usize::MAX)), v, "p={p}");
                    }
                    None => assert_eq!(cost(usize::MAX), u64::MAX, "p={p}"),
                }
            }
        }
    }

    #[test]
    fn reduction_cost_grid_matches_integer_oracle() {
        // Counts: 0 (empty product), 1, 2, and 2^20. The cap pair is the last `n` whose
        // product fits in `u64` and the next `n`, when that `n` fits in `usize`.
        let ns = [0usize, 1, 2, 1usize << 20];
        let ps = [
            1usize,
            2,
            69,
            70,
            1624,
            1625,
            2_643_876,
            1usize << 28,
            crate::defs::MAX_PRECISION,
        ];
        for p in ps {
            let per = cos_per(p);
            for n in ns {
                assert_eq!(
                    CosArgReductionEstimator::reduction_cost(n, p),
                    scale_cost(n, per),
                    "n={n} p={p}"
                );
            }
            assert_cap(p, per, |n| CosArgReductionEstimator::reduction_cost(n, p));
        }
    }
}
