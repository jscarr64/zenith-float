// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! Functons used by macros

use crate::{
    common::util::{count_leading_ones, count_leading_zeroes_skip_first},
    defs::DEFAULT_P,
    Consts, ExactNum, Exponent, RoundingMode, Sign, EXPONENT_BIT_SIZE, INF_NEG, INF_POS,
};

/// Computes error for ExactNum values near 1. This function is for internal use by macro `expr`.
pub fn compute_added_err_near_one(arg: &ExactNum, emin: Exponent) -> usize {
    if arg.is_zero() {
        return 0;
    }

    let n = match arg.exponent() {
        Some(1) => count_leading_zeroes_skip_first(arg.mantissa_digits().unwrap_or(&[])),
        Some(0) => count_leading_ones(arg.mantissa_digits().unwrap_or(&[])),
        _ => 0,
    };

    if n > emin.unsigned_abs() as usize {
        0
    } else {
        n
    }
}

#[derive(Debug)]
pub enum TrigFun {
    Sin,
    Cos,
    Tan,
}

/// Algorithm of error computation.
#[derive(Debug)]
pub enum ErrAlgo<'a> {
    Log(&'a ExactNum, usize, Exponent),
    Log2(&'a ExactNum, &'a ExactNum, Exponent),
    Pow(&'a ExactNum, &'a ExactNum, Exponent),
    Trig(&'a ExactNum, usize, TrigFun, &'a mut Consts, Exponent),
    Asin(&'a ExactNum, Exponent),
    Acos(&'a ExactNum, Exponent),
    Acosh(&'a ExactNum, Exponent),
    Atanh(&'a ExactNum, Exponent),
}

/// Computes the precision increment of an arguments to cover the error for a given algorithm.
#[inline]
pub fn compute_added_err(algo: ErrAlgo<'_>) -> usize {
    match algo {
        ErrAlgo::Log(arg, c, emin) => {
            if arg.inexact_accumulated() && arg.is_positive() {
                c + compute_added_err_near_one(arg, emin)
            } else {
                0
            }
        }
        ErrAlgo::Log2(base, arg, emin) => {
            if (arg.inexact_accumulated() || base.inexact_accumulated())
                && arg.is_positive()
                && base.is_positive()
            {
                let y = if base.inexact_accumulated() {
                    compute_added_err_near_one(base, emin)
                } else {
                    0
                };
                let x = if arg.inexact_accumulated() {
                    compute_added_err_near_one(arg, emin)
                } else {
                    0
                };
                5 + x.max(y)
            } else {
                0
            }
        }
        ErrAlgo::Pow(base, arg, emin) => {
            if arg.inexact_accumulated() || base.inexact_accumulated() {
                if let Some(earg) = arg.exponent() {
                    let v = if !base.inexact_accumulated() || earg <= 0 {
                        0
                    } else {
                        let n = compute_added_err_near_one(base, emin);
                        // In u64: `emin` is near -2^62, which does not fit a 32-bit `usize`.
                        let lim = u64::try_from(EXPONENT_BIT_SIZE)
                            .ok()
                            .and_then(|b| emin.unsigned_abs().checked_add(b));
                        match lim {
                            Some(l) if earg.unsigned_abs() <= l => {
                                usize::try_from(earg).map_or(n, |e| n.min(e))
                            }
                            _ => 0,
                        }
                    };
                    5 + v + EXPONENT_BIT_SIZE
                } else {
                    0
                }
            } else {
                0
            }
        }
        ErrAlgo::Trig(arg, p, fun, cc, emin) => {
            if arg.inexact_accumulated() {
                if let Some(e) = arg.exponent() {
                    let err = if e < 1 {
                        0
                    } else {
                        let err = match fun {
                            TrigFun::Sin => {
                                let pi = cc.pi(p, RoundingMode::None);
                                let mut argrem = arg.rem(&pi);
                                argrem.set_sign(Sign::Pos);
                                if argrem.exponent().unwrap_or(0) > 1 {
                                    argrem = argrem.sub(&pi, p, RoundingMode::None);
                                }
                                if argrem.is_zero() {
                                    if p > emin.unsigned_abs() as usize {
                                        0
                                    } else {
                                        p
                                    }
                                } else {
                                    argrem.exponent().map_or(0, |e| {
                                        if e < 0 && e >= emin {
                                            e.unsigned_abs() as usize
                                        } else {
                                            0
                                        }
                                    })
                                }
                            }
                            TrigFun::Cos => {
                                let mut pi = cc.pi(p, RoundingMode::None);
                                let mut argrem = arg.rem(&pi);
                                if argrem.exponent().unwrap_or(0) == 1 {
                                    argrem.set_sign(Sign::Pos);
                                    pi.set_exponent(1);
                                    argrem = argrem.sub(&pi, p, RoundingMode::None);
                                    if argrem.is_zero() {
                                        if p > emin.unsigned_abs() as usize {
                                            0
                                        } else {
                                            p
                                        }
                                    } else {
                                        argrem.exponent().map_or(0, |e| {
                                            if e < 0 && e >= emin {
                                                e.unsigned_abs() as usize
                                            } else {
                                                0
                                            }
                                        })
                                    }
                                } else {
                                    0
                                }
                            }
                            TrigFun::Tan => {
                                let mut pi = cc.pi(p, RoundingMode::None);
                                pi.set_exponent(1);
                                let mut argrem = arg.rem(&pi);
                                argrem.set_sign(Sign::Pos);
                                if argrem.exponent().unwrap_or(0) > 0 {
                                    argrem = argrem.sub(&pi, p, RoundingMode::None);
                                }
                                if argrem.is_zero() {
                                    if p > emin.unsigned_abs() as usize {
                                        0
                                    } else {
                                        p
                                    }
                                } else {
                                    argrem.exponent().map_or(0, |e| {
                                        if e < emin {
                                            0
                                        } else {
                                            e.unsigned_abs() as usize
                                        }
                                    })
                                }
                            }
                        };

                        err + e.unsigned_abs() as usize
                    };

                    err + 2
                } else {
                    0
                }
            } else {
                0
            }
        }
        ErrAlgo::Asin(arg, emin) => {
            if arg.inexact_accumulated() && arg.exponent().unwrap_or(1) < 1 {
                let n = compute_added_err_near_one(arg, emin);
                2 + n.div_ceil(2)
            } else {
                0
            }
        }
        ErrAlgo::Acos(arg, emin) => {
            if arg.inexact_accumulated() && arg.exponent().unwrap_or(1) < 1 {
                let n = compute_added_err_near_one(arg, emin);
                2 + if arg.is_positive() { n } else { n / 2 }
            } else {
                0
            }
        }
        ErrAlgo::Acosh(arg, emin) => {
            if arg.inexact_accumulated() && arg.is_positive() && arg.exponent().unwrap_or(0) >= 1 {
                2 + compute_added_err_near_one(arg, emin)
            } else {
                0
            }
        }
        ErrAlgo::Atanh(arg, emin) => {
            if arg.inexact_accumulated() && arg.exponent().unwrap_or(1) < 1 {
                2 + compute_added_err_near_one(arg, emin)
            } else {
                0
            }
        }
    }
}

/// Checks if the number's exponent is in the given exponent range.
/// Returns Inf if the exponent of `n` is larger than `emax`.
/// The sign of the infinity is determined by the sign of `n`.
/// Returns 0 if the exponent of `n` is smaller than `emin`.
/// Return `n` itself otherwise.
#[inline]
pub fn check_exponent_range(n: ExactNum, emin: Exponent, emax: Exponent) -> ExactNum {
    if let Some(e) = n.exponent() {
        if e > emax {
            if n.is_positive() {
                INF_POS
            } else {
                INF_NEG
            }
        } else if e < emin {
            ExactNum::new(n.mantissa_max_bit_len().unwrap_or(DEFAULT_P))
        } else {
            n
        }
    } else {
        n
    }
}

/// Applies [`check_exponent_range`] to both parts of a complex value.
#[inline]
pub fn check_complex_exponent_range(
    z: crate::ExactComplex,
    emin: Exponent,
    emax: Exponent,
) -> crate::ExactComplex {
    let (re, im) = z.into_parts();
    crate::ExactComplex::new(
        check_exponent_range(re, emin, emax),
        check_exponent_range(im, emin, emax),
    )
}

/// `Ok(None)`: no cancellation. `Err(InvalidArgument)`: the cancelled-bit count is not a representable
/// working precision (32-bit `usize`), so the result cannot be retried and must not be accepted.
fn part_cancel(
    a: &ExactNum,
    b: &ExactNum,
    r: &ExactNum,
    is_add: bool,
) -> Result<Option<usize>, crate::Error> {
    let (Some(e1), Some(e2), Some(e3)) = (a.exponent(), b.exponent(), r.exponent()) else {
        return Ok(None);
    };
    let (Some(s1), Some(s2)) = (a.sign(), b.sign()) else {
        return Ok(None);
    };
    let opp = if is_add { s1 != s2 } else { s1 == s2 };
    if (e1 - e2).abs() <= 1 && opp {
        // A zero result has exponent 0, so the gap can exceed a 32-bit `usize`: that is an error.
        e1.max(e2)
            .checked_sub(e3)
            .and_then(|gap| usize::try_from(gap.unsigned_abs()).ok())
            .and_then(|g| g.checked_add(1))
            .map(Some)
            .ok_or(crate::Error::InvalidArgument)
    } else {
        Ok(None)
    }
}

/// Extra working bits from cancellation in the **real** and **imaginary** parts,
/// independently. Either side may be `None` when that part did not cancel. An error means a count that
/// is not a representable working precision; `cexpr!` then returns NaN with that error.
///
/// `cexpr!` stores these in two `errs[]` slots so a large imaginary cancellation
/// cannot hide a smaller real one (or the reverse).
#[inline]
pub fn complex_cancel_bits(
    arg1: &crate::ExactComplex,
    arg2: &crate::ExactComplex,
    ret: &crate::ExactComplex,
    is_add: bool,
) -> Result<(Option<usize>, Option<usize>), crate::Error> {
    if !arg1.inexact() && !arg2.inexact() {
        return Ok((None, None));
    }
    Ok((
        part_cancel(arg1.re(), arg2.re(), ret.re(), is_add)?,
        part_cancel(arg1.im(), arg2.im(), ret.im(), is_add)?,
    ))
}

#[cfg(test)]
mod tests {

    fn one() -> crate::ExactNum {
        crate::ExactNum::from_word(1, crate::defs::DEFAULT_P)
    }

    use crate::{
        mantissa::Mantissa, num::ExactNumNumber, Consts, Exponent, RoundingMode, Sign, Word,
        EXPONENT_MAX, EXPONENT_MIN, INF_NEG, INF_POS, NAN, WORD_BIT_SIZE, WORD_MAX,
        WORD_SIGNIFICANT_BIT,
    };

    use super::*;
    use crate::common::util::TEST_EXP_BOUND;

    #[cfg(not(feature = "std"))]
    use alloc::vec;

    fn assert(expected: usize, words: &[Word], s: Sign, e: Exponent) {
        let d = ExactNum::from_words(words, s, e);
        assert_eq!(
            compute_added_err_near_one(&d, EXPONENT_MIN),
            expected,
            "Expected error {} for d = {:?}",
            expected,
            d
        );
    }

    #[test]
    fn test_compute_err_near_one() {
        for s in [Sign::Pos, Sign::Neg] {
            // 0
            assert(0, &[0, 0, 0], s, 0);

            // 1 and 0.5
            assert(1, &[0, 0, WORD_SIGNIFICANT_BIT], s, 0);
            assert(WORD_BIT_SIZE * 3, &[0, 0, WORD_SIGNIFICANT_BIT], s, 1);

            // close to 1
            assert(WORD_BIT_SIZE, &[0, 0, WORD_MAX], s, 0);
            assert(WORD_BIT_SIZE, &[0, WORD_MAX, WORD_SIGNIFICANT_BIT], s, 1);

            // exponent is neither 0 nor 1
            assert(0, &[0, 123, WORD_MAX], s, -1);
            assert(0, &[0, 123, WORD_SIGNIFICANT_BIT], s, 2);
        }

        // inf and nan
        assert_eq!(
            compute_added_err_near_one(&INF_POS, EXPONENT_MIN),
            0,
            "Expected error 0 for INF_POS"
        );
        assert_eq!(
            compute_added_err_near_one(&INF_NEG, EXPONENT_MIN),
            0,
            "Expected error 0 for INF_NEG"
        );
        assert_eq!(
            compute_added_err_near_one(&NAN, EXPONENT_MIN),
            0,
            "Expected error 0 for NAN"
        );
    }

    fn gen_pair(m1: Mantissa, mut e: Exponent) -> (ExactNum, ExactNum) {
        let s = if crate::common::test_rng::random::<i8>() & 1 == 0 {
            Sign::Pos
        } else {
            Sign::Neg
        };

        let mut m2 = m1.clone().unwrap();

        let n1 = ExactNumNumber::from_raw_unchecked(m1, s, e, true);

        if m2.add_ulp() {
            assert!(e < EXPONENT_MAX);
            if let Some(top) = m2.digits_mut().iter_mut().last() {
                *top = WORD_SIGNIFICANT_BIT;
            }
            e += 1;
        }

        let n2 = ExactNumNumber::from_raw_unchecked(m2, s, e, true);

        (n1.into(), n2.into())
    }

    fn gen_num_pair(
        p: usize,
        e: Exponent,
        near1: i32,
    ) -> Result<(ExactNum, ExactNum), crate::Error> {
        let mut m1 = Mantissa::random_normal(p)?;

        if near1 != 0 {
            let mut bits = crate::common::test_rng::random::<usize>() % (p - 1) + 1;
            let mut iter = m1.digits_mut().iter_mut().rev();

            while bits >= WORD_BIT_SIZE {
                if near1.is_negative() {
                    *iter.next().ok_or(crate::Error::InvalidArgument)? = WORD_MAX;
                } else {
                    *iter.next().ok_or(crate::Error::InvalidArgument)? = 0;
                }
                bits -= WORD_BIT_SIZE;
            }

            if bits > 0 {
                if near1.is_negative() {
                    *iter.next().ok_or(crate::Error::InvalidArgument)? |=
                        WORD_MAX << (WORD_BIT_SIZE - bits);
                } else {
                    *iter.next().ok_or(crate::Error::InvalidArgument)? &= WORD_MAX >> bits;
                }
            }

            *m1.digits_mut()
                .iter_mut()
                .last()
                .ok_or(crate::Error::InvalidArgument)? |= WORD_SIGNIFICANT_BIT;
        }

        Ok(gen_pair(m1, e))
    }

    fn gen_num_pair_trig(
        p: usize,
        e: Exponent,
        near_pi: i32,
        cc: &mut Consts,
        add_half_pi: bool,
    ) -> Result<(ExactNum, ExactNum), crate::Error> {
        let mut m1 = Mantissa::random_normal(p)?;

        if near_pi != 0 {
            let bits = crate::common::test_rng::random::<usize>() % (p - 1) + 1;

            let pirm = if near_pi < 0 { RoundingMode::ToZero } else { RoundingMode::FromZero };

            let mut pi = cc.pi(p, pirm);
            pi.set_exponent(0);
            let pi_rounded = pi.round(bits, pirm);
            let mut pi = cc.pi(p, pirm);
            let pi_digits = pi_rounded
                .mantissa_digits()
                .ok_or(crate::Error::InvalidArgument)?;

            let l = m1.len();
            m1.digits_mut()[l - pi_digits.len() + 1..].copy_from_slice(&pi_digits[1..]);
            let shift = bits % WORD_BIT_SIZE;
            if shift > 0 {
                m1.digits_mut()[l - pi_digits.len()] &= WORD_MAX >> shift;
                m1.digits_mut()[l - pi_digits.len()] |= pi_digits[0];
            }

            if add_half_pi {
                let n = ExactNumNumber::from_raw_unchecked(m1, Sign::Pos, e, true);
                pi.set_exponent(1);
                let n = pi.add(&n.into(), p, RoundingMode::ToZero);
                let w = n.value_parts().ok_or(crate::Error::InvalidArgument)?.0;
                m1 = Mantissa::from_words(w.len() * WORD_BIT_SIZE, w)?;
            }
        }

        Ok(gen_pair(m1, e))
    }

    fn calc_err(mut d1: ExactNum, mut d2: ExactNum, p: usize) -> Result<usize, crate::Error> {
        Ok(if d1 == d2 {
            0
        } else {
            if let (Some(e1), Some(e2)) = (d1.exponent(), d2.exponent()) {
                let emax = e1.max(e2);
                d1.set_exponent(e1 - emax);
                d2.set_exponent(e2 - emax);
                let a = d1
                    .sub(&d2, p, RoundingMode::None)
                    .exponent()
                    .ok_or(crate::Error::InvalidArgument)?
                    .unsigned_abs() as usize;
                p.saturating_sub(a)
            } else {
                0
            }
        })
    }

    fn err_exp_samples(ernd: Exponent) -> [Exponent; 5] {
        [
            0,
            1,
            2,
            EXPONENT_BIT_SIZE as Exponent + ernd,
            TEST_EXP_BOUND.min(256 + ernd.min(64)),
        ]
    }

    const ERR_TEST_OUTER: usize = 10;

    #[test]
    fn test_compute_added_err() -> Result<(), crate::Error> {
        // Own stream: random_dist tests reseed the shared generator in parallel.
        crate::common::test_rng::reseed_random_this_thread(
            crate::common::test_rng::DEFAULT_RANDOM_SEED,
        );
        let mut cc = Consts::new()?;
        let emin = EXPONENT_MIN;

        for _ in 0..ERR_TEST_OUTER {
            let ernd =
                crate::common::test_rng::random::<Exponent>() % (EXPONENT_BIT_SIZE as Exponent - 5);

            for e in err_exp_samples(ernd) {
                for esign in [1, -1] {
                    let near1set = if e <= 1 { vec![0, -1, 1] } else { vec![0] };

                    for near1 in near1set {
                        let e = e * esign;
                        let p =
                            (crate::common::test_rng::random::<usize>() % 10 + 1) * WORD_BIT_SIZE;

                        let (n1, n2) = gen_num_pair(p, e, near1)?;

                        //println!("n1 {:?}", n1);
                        //println!("n2 {:?}", n2);

                        // reciprocal
                        let d1 = ExactNum::reciprocal(&n1, p, RoundingMode::None);
                        let d2 = ExactNum::reciprocal(&n2, p, RoundingMode::None);
                        let err = calc_err(d1, d2, p)?;

                        //println!("recip {:?}", err);
                        assert!(err <= 2);

                        // sqrt
                        let d1 = ExactNum::sqrt(&n1, p, RoundingMode::None);
                        let d2 = ExactNum::sqrt(&n2, p, RoundingMode::None);
                        let err = calc_err(d1, d2, p)?;

                        //println!("sqrt {:?}", err);
                        assert!(err <= 1);

                        // cbrt
                        let d1 = ExactNum::cbrt(&n1, p, RoundingMode::None);
                        let d2 = ExactNum::cbrt(&n2, p, RoundingMode::None);
                        let err = calc_err(d1, d2, p)?;

                        //println!("cbrt {:?}", err);
                        assert!(err <= 1);

                        // ln
                        let d1 = ExactNum::ln(&n1, p, RoundingMode::None, &mut cc);
                        let d2 = ExactNum::ln(&n2, p, RoundingMode::None, &mut cc);

                        let err_estimate = compute_added_err(ErrAlgo::Log(&n1, 2, emin));
                        let err = calc_err(d1, d2, p)?;

                        //println!("ln {:?} {:?}", err, err_estimate);
                        assert!(err <= err_estimate);

                        // log base 2
                        let d1 = ExactNum::log2(&n1, p, RoundingMode::None, &mut cc);
                        let d2 = ExactNum::log2(&n2, p, RoundingMode::None, &mut cc);

                        let err_estimate = compute_added_err(ErrAlgo::Log(&n1, 3, emin));
                        let err = calc_err(d1, d2, p)?;

                        //println!("log2 {:?} {:?}", err, err_estimate);
                        assert!(err <= err_estimate);

                        // log base 10
                        let d1 = ExactNum::log10(&n1, p, RoundingMode::None, &mut cc);
                        let d2 = ExactNum::log10(&n2, p, RoundingMode::None, &mut cc);

                        let err_estimate = compute_added_err(ErrAlgo::Log(&n1, 6, emin));
                        let err = calc_err(d1, d2, p)?;

                        //println!("log10 {:?} {:?}", err, err_estimate);
                        assert!(err <= err_estimate);

                        // log base b, pow
                        let mut nc = [n1.clone(), n2.clone()];

                        for e2 in err_exp_samples(ernd) {
                            for esign2 in [1, -1] {
                                let near1set2 = if e2 <= 1 { vec![0, -1, 1] } else { vec![0] };

                                for near12 in near1set2 {
                                    let e2 = e2 * esign2;
                                    let (b1, b2) = gen_num_pair(p, e2, near12)?;
                                    let mut bc = [b1, b2];

                                    for (i1, i2) in [(0, 0), (0, 1), (1, 0), (1, 1)] {
                                        for (j1, j2) in [(0, 0), (0, 1), (1, 0), (1, 1)] {
                                            if i1 == i2 && j1 == j2 {
                                                continue;
                                            }

                                            if i1 == i2 {
                                                nc[0].set_inexact(false);
                                                nc[1].set_inexact(false);
                                            } else {
                                                nc[0].set_inexact(true);
                                                nc[1].set_inexact(true);
                                            }

                                            if j1 == j2 {
                                                bc[0].set_inexact(false);
                                                bc[1].set_inexact(false);
                                            } else {
                                                bc[0].set_inexact(true);
                                                bc[1].set_inexact(true);
                                            }

                                            let n = &nc[i1];
                                            let nn = &nc[i2];

                                            let b = &bc[j1];
                                            let bb = &bc[j2];

                                            let d1 =
                                                ExactNum::log(n, b, p, RoundingMode::None, &mut cc);
                                            let d2 = ExactNum::log(
                                                nn,
                                                bb,
                                                p,
                                                RoundingMode::None,
                                                &mut cc,
                                            );

                                            let err_estimate =
                                                compute_added_err(ErrAlgo::Log2(b, n, emin));
                                            let err = calc_err(d1, d2, p)?;

                                            let err_estimate = 1.max(err_estimate);

                                            //println!("log {:?} {:?}", err, err_estimate);
                                            assert!(err <= err_estimate);

                                            let d1 =
                                                ExactNum::pow(b, n, p, RoundingMode::None, &mut cc);
                                            let d2 = ExactNum::pow(
                                                bb,
                                                nn,
                                                p,
                                                RoundingMode::None,
                                                &mut cc,
                                            );

                                            let err_estimate =
                                                compute_added_err(ErrAlgo::Pow(b, n, emin));
                                            let err = calc_err(d1, d2, p)?;

                                            if !(b.abs_cmp(&one()) == Some(0)
                                                || bb.abs_cmp(&one()) == Some(0))
                                            {
                                                let err_estimate = 1.max(err_estimate);
                                                //println!("pow {:?} {:?}", err, err_estimate);
                                                assert!(err <= err_estimate);
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        // exp
                        let d1 = ExactNum::exp(&n1, p, RoundingMode::None, &mut cc);
                        let d2 = ExactNum::exp(&n2, p, RoundingMode::None, &mut cc);

                        let err = calc_err(d1, d2, p)?;

                        //println!("exp {:?}", err);
                        assert!(err <= EXPONENT_BIT_SIZE + 1);

                        // sinh
                        let d1 = ExactNum::sinh(&n1, p, RoundingMode::None, &mut cc);
                        let d2 = ExactNum::sinh(&n2, p, RoundingMode::None, &mut cc);

                        let err = calc_err(d1, d2, p)?;

                        //println!("sinh {:?}", err);
                        assert!(err <= EXPONENT_BIT_SIZE + 1);

                        // cosh
                        let d1 = ExactNum::cosh(&n1, p, RoundingMode::None, &mut cc);
                        let d2 = ExactNum::cosh(&n2, p, RoundingMode::None, &mut cc);

                        let err = calc_err(d1, d2, p)?;

                        //println!("cosh {:?}", err);
                        assert!(err <= EXPONENT_BIT_SIZE + 1);

                        // tanh
                        let d1 = ExactNum::tanh(&n1, p, RoundingMode::None, &mut cc);
                        let d2 = ExactNum::tanh(&n2, p, RoundingMode::None, &mut cc);

                        let err = calc_err(d1, d2, p)?;

                        //println!("tanh {:?}", err);
                        assert!(err <= 2);

                        // asinh
                        let d1 = ExactNum::asinh(&n1, p, RoundingMode::None, &mut cc);
                        let d2 = ExactNum::asinh(&n2, p, RoundingMode::None, &mut cc);

                        let err = calc_err(d1, d2, p)?;

                        //println!("asinh {:?}", err);
                        assert!(err <= 2);

                        // acosh
                        let d1 = ExactNum::acosh(&n1, p, RoundingMode::None, &mut cc);
                        let d2 = ExactNum::acosh(&n2, p, RoundingMode::None, &mut cc);

                        let err_estimate = compute_added_err(ErrAlgo::Acosh(&n1, emin));
                        let err = calc_err(d1, d2, p)?;

                        //println!("acosh {:?} {:?}", err, err_estimate);
                        assert!(err <= err_estimate);

                        // atanh
                        let d1 = ExactNum::atanh(&n1, p, RoundingMode::None, &mut cc);
                        let d2 = ExactNum::atanh(&n2, p, RoundingMode::None, &mut cc);

                        let err_estimate = compute_added_err(ErrAlgo::Atanh(&n1, emin));
                        let err = calc_err(d1, d2, p)?;

                        //println!("atanh {:?} {:?}", err, err_estimate);
                        assert!(err <= err_estimate);

                        // atan
                        let d1 = ExactNum::atan(&n1, p, RoundingMode::None, &mut cc);
                        let d2 = ExactNum::atan(&n2, p, RoundingMode::None, &mut cc);

                        let err = calc_err(d1, d2, p)?;

                        //println!("atan {:?}", err);
                        assert!(err <= 2);
                    }
                }
            }

            // sin, cos, tan
            for e in err_exp_samples(ernd) {
                for esign in [1, -1] {
                    for near_pi in [0, 1, -1] {
                        for add_half_pi in [false, true] {
                            let e = e * esign;
                            let p = (crate::common::test_rng::random::<usize>() % 10 + 1)
                                * WORD_BIT_SIZE;

                            let (n1, n2) = gen_num_pair_trig(p, e, near_pi, &mut cc, add_half_pi)?;

                            // sin
                            let d1 = ExactNum::sin(&n1, p, RoundingMode::None, &mut cc);
                            let d2 = ExactNum::sin(&n2, p, RoundingMode::None, &mut cc);
                            let err_estimate = compute_added_err(ErrAlgo::Trig(
                                &n1,
                                p,
                                TrigFun::Sin,
                                &mut cc,
                                emin,
                            ));
                            let err = calc_err(d1, d2, p)?;

                            //println!("sin {:?} {:?}", err, err_estimate);
                            assert!(err <= err_estimate);

                            // cos
                            let d1 = ExactNum::cos(&n1, p, RoundingMode::None, &mut cc);
                            let d2 = ExactNum::cos(&n2, p, RoundingMode::None, &mut cc);
                            let err_estimate = compute_added_err(ErrAlgo::Trig(
                                &n1,
                                p,
                                TrigFun::Cos,
                                &mut cc,
                                emin,
                            ));
                            let err = calc_err(d1, d2, p)?;

                            //println!("cos {:?} {:?}", err, err_estimate);
                            assert!(err <= err_estimate);

                            // tan
                            let d1 = ExactNum::tan(&n1, p, RoundingMode::None, &mut cc);
                            let d2 = ExactNum::tan(&n2, p, RoundingMode::None, &mut cc);

                            let err_estimate = compute_added_err(ErrAlgo::Trig(
                                &n1,
                                p,
                                TrigFun::Tan,
                                &mut cc,
                                emin,
                            ));
                            let err = calc_err(d1, d2, p)?;

                            //println!("tan {:?} {:?}", err, err_estimate);
                            assert!(err <= err_estimate);
                        }
                    }
                }
            }

            // asin, acos
            for e in err_exp_samples(ernd).into_iter().map(|v| -v) {
                let near1set = if e == 0 { vec![0, -1] } else { vec![0] };

                for near1 in near1set {
                    let p = (crate::common::test_rng::random::<usize>() % 10 + 1) * WORD_BIT_SIZE;

                    let (n1, n2) = gen_num_pair(p, e, near1)?;

                    // asin
                    let d1 = ExactNum::asin(&n1, p, RoundingMode::None, &mut cc);
                    let d2 = ExactNum::asin(&n2, p, RoundingMode::None, &mut cc);

                    let err_estimate = compute_added_err(ErrAlgo::Asin(&n1, emin));
                    let err = calc_err(d1, d2, p)?;

                    //println!("asin {:?} {:?}", err, err_estimate);
                    assert!(err <= err_estimate);

                    // acos
                    let d1 = ExactNum::acos(&n1, p, RoundingMode::None, &mut cc);
                    let d2 = ExactNum::acos(&n2, p, RoundingMode::None, &mut cc);

                    let err_estimate = compute_added_err(ErrAlgo::Acos(&n1, emin));
                    let err = calc_err(d1, d2, p)?;

                    //println!("acos {:?} {:?}", err, err_estimate);
                    assert!(err <= err_estimate);
                }
            }
        }
        Ok(())
    }

    #[test]
    fn test_complex_cancel_independent_parts() -> Result<(), crate::Error> {
        let p = 256;
        let rm = RoundingMode::None;
        let mk_inexact = |n: u8| {
            let mut v = ExactNum::from_u8(n, p);
            v.set_inexact(true);
            v
        };
        let tiny = ExactNum::from_u8(1, p).ldexp(-80, p, rm);
        let mut near_m1 = ExactNum::from_i8(-1, p).add(&tiny, p, rm);
        let re_pos = mk_inexact(1);
        near_m1.set_inexact(true);
        let im = mk_inexact(2);

        let a = crate::ExactComplex::new(re_pos, mk_inexact(2));
        let b = crate::ExactComplex::new(near_m1, im);
        let ret = crate::ExactComplex::add_ref(&a, &b, p, rm);
        let (re_err, im_err) = complex_cancel_bits(&a, &b, &ret, true)?;
        assert!(re_err.is_some(), "real part should report cancellation");
        assert!(
            im_err.is_none(),
            "imaginary part should not share the real slot"
        );

        let im_pos = mk_inexact(2);
        let mut im_neg = ExactNum::from_i8(-2, p).add(&tiny, p, rm);
        im_neg.set_inexact(true);
        let c = crate::ExactComplex::new(mk_inexact(1), im_pos);
        let d = crate::ExactComplex::new(mk_inexact(1), im_neg);
        let ret2 = crate::ExactComplex::add_ref(&c, &d, p, rm);
        let (re_err2, im_err2) = complex_cancel_bits(&c, &d, &ret2, true)?;
        assert!(re_err2.is_none(), "real part should stay independent");
        assert!(
            im_err2.is_some(),
            "imaginary part should report cancellation"
        );
        Ok(())
    }
}
