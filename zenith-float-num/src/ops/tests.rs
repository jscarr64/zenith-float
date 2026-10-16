// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! tests

use crate::common::consts::ONE;
use crate::common::test_rng::random;
use crate::common::util::TEST_EXP_BOUND;
use crate::common::util::{count_leading_ones, count_leading_zeroes_skip_first, log2_floor};
use crate::defs::{RoundingMode, EXPONENT_MIN, WORD_BIT_SIZE};
use crate::num::ExactNumNumber;
use crate::ops::consts::Consts;
use crate::{Exponent, Sign};

const TEST_ITERS: usize = 256;

const fn get_prec_rng() -> usize {
    157
}

#[test]
fn test_ln_exp() -> Result<(), crate::Error> {
    let prec_rng = get_prec_rng();
    let mut eps = ExactNumNumber::small_const(1);

    let mut cc = Consts::new()?;

    for i in 0..TEST_ITERS {
        let p1 = (crate::common::test_rng::random::<usize>() % prec_rng + 1) * WORD_BIT_SIZE;
        let prec = (crate::common::test_rng::random::<usize>() % prec_rng + 1) * WORD_BIT_SIZE;

        if i & 1 == 0 {
            let mut d1 = ExactNumNumber::random_normal(p1, -TEST_EXP_BOUND, TEST_EXP_BOUND)?;
            d1.set_sign(Sign::Pos);

            let d2 = d1.ln(prec, RoundingMode::ToEven, &mut cc)?;
            let d3 = d2.exp(prec, RoundingMode::ToEven, &mut cc)?;

            eps.set_exponent(d2.exponent() - prec as Exponent);
            let err = eps.exp(prec, RoundingMode::Up, &mut cc)?;

            //println!("{}", d1.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            //println!("{}", d2.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            //println!("{}", d3.format(crate::Radix::Bin, RoundingMode::None).unwrap());

            // d2 - ulp(d2)/2 <= ln(d1) <= d2 + ulp(d2)/2  ->  d3 / e^(ulp(d2)/2) <= d1 <= d3 * e^(ulp(d2)/2)
            assert!(
                d1.cmp(&d3.mul(&err, prec, RoundingMode::Up)?) <= 0,
                "{} {:?}",
                prec,
                d1
            );
            assert!(
                d1.cmp(&d3.div(&err, prec, RoundingMode::Down)?) >= 0,
                "{} {:?}",
                prec,
                d1
            );
        } else {
            let emax = log2_floor(TEST_EXP_BOUND as usize) as Exponent;
            let emin = -emax;
            let d1 = ExactNumNumber::random_normal(p1, emin, emax)?;

            let d2 = d1.exp(prec, RoundingMode::ToEven, &mut cc)?;

            let d3 = d2.ln(prec, RoundingMode::ToEven, &mut cc)?;

            if d1.exponent() < 1 {
                let addexp = if d1.is_negative() {
                    count_leading_ones
                } else {
                    count_leading_zeroes_skip_first
                }(d2.mantissa().digits()) as Exponent;

                eps.set_exponent(d1.exponent() - prec as Exponent + addexp + 1);
            } else {
                eps.set_exponent(d1.exponent() - prec as Exponent + 1);
            }

            // println!("{}", d1.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("{}", d2.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("{}", d3.format(crate::Radix::Bin, RoundingMode::None).unwrap());

            assert!(
                d1.sub(&d3, prec, RoundingMode::ToEven)?.abs()?.cmp(&eps) <= 0,
                "{} {:?}",
                prec,
                d1
            );
        }
    }
    Ok(())
}

#[test]
fn test_powi() -> Result<(), crate::Error> {
    let prec_rng = get_prec_rng();

    for _ in 0..TEST_ITERS {
        let i = random::<usize>() % 1000 + 1;
        let p1 = (crate::common::test_rng::random::<usize>() % prec_rng + 1) * WORD_BIT_SIZE;
        let prec = (crate::common::test_rng::random::<usize>() % prec_rng + 1) * WORD_BIT_SIZE;
        // d1 is built twice from one seed (the reference value needs its own owner); the stream then resumes from a
        // fresh seed so later draws do not repeat.
        let build_d1 = |seed: u64| -> Result<ExactNumNumber, crate::Error> {
            let resume = random::<u64>();
            crate::common::test_rng::pin_thread_seed(seed);
            let built = ExactNumNumber::random_normal(
                p1,
                -TEST_EXP_BOUND / i as Exponent,
                TEST_EXP_BOUND / i as Exponent,
            );
            crate::common::test_rng::pin_thread_seed(resume);
            let mut d = built?;
            d.set_sign(Sign::Pos);
            Ok(d)
        };
        let seed = random::<u64>();
        let d1 = build_d1(seed)?;

        let d2 = d1.powi(i, prec, RoundingMode::ToEven)?;

        let mut d3 = build_d1(seed)?;
        d3.set_precision(prec + core::mem::size_of::<usize>(), RoundingMode::None)?;
        for _ in 1..i {
            d3 = d3.mul(&d1, d3.mantissa_max_bit_len(), RoundingMode::None)?;
        }
        d3.set_precision(prec, RoundingMode::ToEven)?;

        // println!("d1 {}", d1.format(crate::Radix::Bin, RoundingMode::None).unwrap());
        // println!("d2 {}", d2.format(crate::Radix::Bin, RoundingMode::None).unwrap());
        // println!("d3 {}", d3.format(crate::Radix::Bin, RoundingMode::None).unwrap());
        // println!("i {}", i);

        assert!(d3.cmp(&d2) == 0, "{} {} {:?}", i, prec, d1);
    }
    Ok(())
}

#[test]
fn test_log2_log10_pow() -> Result<(), crate::Error> {
    let prec_rng = get_prec_rng();
    let mut eps = ExactNumNumber::small_const(1);

    let mut cc = Consts::new()?;

    for i in 0..TEST_ITERS {
        let p1 = (crate::common::test_rng::random::<usize>() % prec_rng + 1) * WORD_BIT_SIZE;
        let prec = (crate::common::test_rng::random::<usize>() % prec_rng + 1) * WORD_BIT_SIZE;

        if i & 1 == 0 {
            let mut d1 = ExactNumNumber::random_normal(p1, -TEST_EXP_BOUND, TEST_EXP_BOUND)?;
            d1.set_sign(Sign::Pos);

            let d2 = d1.log2(prec, RoundingMode::ToEven, &mut cc)?;
            let two = ExactNumNumber::from_word(2, prec)?;
            let d3 = two.pow(&d2, prec, RoundingMode::ToEven, &mut cc)?;

            let d4 = d1.log10(prec, RoundingMode::ToEven, &mut cc)?;
            let ten = ExactNumNumber::from_word(10, prec)?;
            let d5 = ten.pow(&d4, prec, RoundingMode::ToEven, &mut cc)?;

            //println!("{}", d1.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            //println!("{}", d2.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            //println!("{}", d3.format(crate::Radix::Bin, RoundingMode::None).unwrap());

            // d2 - ulp(d2)/2 <= log2(d1) <= d2 + ulp(d2)/2  ->  d3 / 2^(ulp(d2)/2) <= d1 <= d3 * 2^(ulp(d2)/2)
            eps.set_exponent(d2.exponent() - prec as Exponent);
            let err = two.pow(&eps, prec, RoundingMode::Up, &mut cc)?;

            assert!(
                d1.cmp(&d3.mul(&err, prec, RoundingMode::Up)?) <= 0,
                "{} {:?}",
                prec,
                d1
            );
            assert!(
                d1.cmp(&d3.div(&err, prec, RoundingMode::Down)?) >= 0,
                "{} {:?}",
                prec,
                d1
            );

            eps.set_exponent(d4.exponent() - prec as Exponent);
            let err = ten.pow(&eps, prec, RoundingMode::Up, &mut cc)?;

            assert!(
                d1.cmp(&d5.mul(&err, prec, RoundingMode::Up)?) <= 0,
                "{} {:?}",
                prec,
                d1
            );
            assert!(
                d1.cmp(&d5.div(&err, prec, RoundingMode::Down)?) >= 0,
                "{} {:?}",
                prec,
                d1
            );
        } else {
            let emax = log2_floor(TEST_EXP_BOUND as usize) as Exponent;
            let emin = -emax;
            let d1 = ExactNumNumber::random_normal(p1, emin, emax)?;

            let two = ExactNumNumber::from_word(2, prec)?;
            let d2 = two.pow(&d1, prec, RoundingMode::ToEven, &mut cc)?;
            let d3 = d2.log2(prec, RoundingMode::ToEven, &mut cc)?;

            let ten = ExactNumNumber::from_word(10, prec)?;
            let d4 = ten.pow(&d1, prec, RoundingMode::ToEven, &mut cc)?;
            let d5 = d4.log10(prec, RoundingMode::ToEven, &mut cc)?;

            let set_eps = |x, eps: &mut ExactNumNumber| {
                if d1.exponent() < 1 {
                    let addexp = if d1.is_negative() {
                        count_leading_ones
                    } else {
                        count_leading_zeroes_skip_first
                    }(x) as Exponent;

                    eps.set_exponent(d1.exponent() - prec.min(p1) as Exponent + addexp + 2);
                } else {
                    eps.set_exponent(d1.exponent() - prec.min(p1) as Exponent + 2);
                }
            };

            // println!("{}", d1.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("{}", d2.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("{}", d3.format(crate::Radix::Bin, RoundingMode::None).unwrap());

            set_eps(d2.mantissa().digits(), &mut eps);

            assert!(
                d1.sub(&d3, prec, RoundingMode::ToEven)?.abs()?.cmp(&eps) < 0,
                "{} {:?}",
                prec,
                d1
            );

            set_eps(d4.mantissa().digits(), &mut eps);

            assert!(
                d1.sub(&d5, prec, RoundingMode::ToEven)?.abs()?.cmp(&eps) < 0,
                "{} {:?}",
                prec,
                d1
            );
        }
    }
    Ok(())
}

#[test]
fn test_log_pow() -> Result<(), crate::Error> {
    let prec_rng = get_prec_rng();
    let mut eps = ExactNumNumber::small_const(1);

    let mut cc = Consts::new()?;

    for i in 0..TEST_ITERS {
        let p1 = (crate::common::test_rng::random::<usize>() % prec_rng + 1) * WORD_BIT_SIZE;
        let p2 = (crate::common::test_rng::random::<usize>() % prec_rng + 1) * WORD_BIT_SIZE;
        let prec = (crate::common::test_rng::random::<usize>() % prec_rng + 1) * WORD_BIT_SIZE;
        if i & 1 == 0 {
            let mut d1 = ExactNumNumber::random_normal(p1, -TEST_EXP_BOUND, TEST_EXP_BOUND)?;
            let mut b = ExactNumNumber::random_normal(p2, -TEST_EXP_BOUND, TEST_EXP_BOUND)?;
            d1.set_sign(Sign::Pos);
            b.set_sign(Sign::Pos);

            let d2 = d1.log(&b, prec, RoundingMode::ToEven, &mut cc)?;
            let d3 = b.pow(&d2, prec, RoundingMode::ToEven, &mut cc)?;

            // println!("d1 {}", d1.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("d2 {}", d2.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("d3 {}", d3.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("b  {}", b.format(crate::Radix::Bin, RoundingMode::None).unwrap());

            // d2 - ulp(d2)/2 <= log_b(d1) <= d2 + ulp(d2)/2  ->  d3 / b^(ulp(d2)/2) <= d1 <= d3 * b^(ulp(d2)/2)
            eps.set_exponent(d2.exponent() - prec as Exponent);
            if b.exponent() <= 0 {
                let b_exp = b.exponent().unsigned_abs() as usize;
                let ln_b_bits = log2_floor(b_exp.max(2)) as Exponent + 1;
                eps.set_exponent(eps.exponent() + ln_b_bits);
            }
            let err = b.pow(&eps, prec, RoundingMode::Up, &mut cc)?;
            let err = if err.is_zero() { ExactNumNumber::min_positive(prec)? } else { err };

            if b.exponent() > 0 {
                assert!(
                    d1.cmp(&d3.mul(&err, prec, RoundingMode::Up)?) <= 0,
                    "{} {:?} {:?}",
                    prec,
                    d1,
                    b
                );
                assert!(
                    d1.cmp(&d3.div(&err, prec, RoundingMode::Down)?) >= 0,
                    "{} {:?} {:?}",
                    prec,
                    d1,
                    b
                );
            } else {
                assert!(
                    d1.cmp(&d3.div(&err, prec, RoundingMode::Down)?) <= 0,
                    "{} {:?} {:?}",
                    prec,
                    d1,
                    b
                );
                assert!(
                    d1.cmp(&d3.mul(&err, prec, RoundingMode::Up)?) >= 0,
                    "{} {:?} {:?}",
                    prec,
                    d1,
                    b
                );
            }
        } else {
            let mut b = ExactNumNumber::random_normal(p2, -TEST_EXP_BOUND, TEST_EXP_BOUND)?;
            b.set_sign(Sign::Pos);

            // if b close to 1, error increases significantly, i.e.
            // let d2 - err <= b^d1 <= d2 + err, then log_b(d2 - err) <= d1 <= log_b(d2 + err),
            // let d2 - err <= b^d1 <= d2 + err, then log_b(d2 - err) <= d1 <= log_b(d2 + err),
            // let d2 - err <= b^d1 <= d2 + err, then log_b(d2 - err) <= d1 <= log_b(d2 + err),
            // and log_b(x) has steep derivative 1 / x / ln(b).
            let mut berr = 0;
            if b.exponent() == 0 {
                berr = count_leading_ones(b.mantissa().digits()) as Exponent;
            } else if b.exponent() == 1 {
                berr = count_leading_zeroes_skip_first(b.mantissa().digits()) as Exponent;
            }

            let n = b.exponent().unsigned_abs() as usize;
            let emax = log2_floor(TEST_EXP_BOUND as usize / if n == 0 { 1 } else { n }) as Exponent;
            let emin = -emax;
            let d1 = ExactNumNumber::random_normal(p1, emin, emax)?;

            let d2 = b.pow(&d1, prec, RoundingMode::ToEven, &mut cc)?;
            let d3 = d2.log(&b, prec, RoundingMode::ToEven, &mut cc)?;

            // println!("\nb  {}", b.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("d1 {}", d1.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("d2 {}", d2.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("d3 {}", d3.format(crate::Radix::Bin, RoundingMode::None).unwrap());

            if d1.exponent() < 1 {
                let addexp = if (b.exponent() > 0 && d1.is_negative())
                    || b.exponent() <= 0 && d1.is_positive()
                {
                    count_leading_ones
                } else {
                    count_leading_zeroes_skip_first
                }(d2.mantissa().digits()) as Exponent;

                eps.set_exponent(d1.exponent() - prec.min(p1) as Exponent + addexp + berr + 2);
            } else {
                eps.set_exponent(d1.exponent() - prec.min(p1) as Exponent + berr + 2);
            }

            assert!(
                d1.sub(&d3, prec, RoundingMode::ToEven)?.abs()?.cmp(&eps) < 0,
                "{} {:?} {:?}",
                prec,
                d1,
                b
            );
        }
    }
    Ok(())
}

#[test]
fn test_sin_asin() -> Result<(), crate::Error> {
    let prec_rng = get_prec_rng();
    let mut eps = ExactNumNumber::small_const(1);
    let mut thres = ExactNumNumber::small_const(1);
    let thres_exp = -8;
    thres.set_exponent(thres_exp);

    let mut cc = Consts::new()?;

    let pi = cc.pi_num((prec_rng + 1) * WORD_BIT_SIZE, RoundingMode::None)?;

    let mut half_pi = cc.pi_num((prec_rng + 1) * WORD_BIT_SIZE, RoundingMode::None)?;
    half_pi.set_exponent(1);

    // argument between -pi/2, pi/2
    for i in 0..TEST_ITERS {
        let p1 = (crate::common::test_rng::random::<usize>() % prec_rng + 1) * WORD_BIT_SIZE;
        let prec = (crate::common::test_rng::random::<usize>() % prec_rng + 1) * WORD_BIT_SIZE;

        if i & 1 == 0 {
            let mut d1 = ExactNumNumber::random_normal(p1, -100, 2)?;

            // -pi/2, pi/2
            while d1.abs()?.cmp(&half_pi) > 0 {
                if d1.is_positive() {
                    d1 = d1.sub(&half_pi, p1, RoundingMode::None)?;
                }
                if d1.is_negative() {
                    d1 = d1.add(&half_pi, p1, RoundingMode::None)?;
                }
            }

            let d2 = d1.sin(prec, RoundingMode::ToEven, &mut cc)?;
            let d3 = d2.asin(prec, RoundingMode::ToEven, &mut cc)?;

            // println!("{}", d1.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("{}", d2.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("{}", d3.format(crate::Radix::Bin, RoundingMode::None).unwrap());

            eps.set_exponent(
                d1.exponent() - prec as Exponent
                    + 1
                    + count_leading_ones(d2.mantissa().digits()) as Exponent,
            );

            assert!(
                d1.sub(&d3, prec, RoundingMode::ToEven)?.abs()?.cmp(&eps) < 0,
                "{} {:?}",
                prec,
                d1
            );
        } else {
            let d1 = ExactNumNumber::random_normal(p1, -(prec as Exponent), 0)?;

            let d2 = d1.asin(prec, RoundingMode::ToEven, &mut cc)?;
            let d3 = d2.sin(prec, RoundingMode::ToEven, &mut cc)?;

            // println!("{}", d1.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("{}", d2.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("{}", d3.format(crate::Radix::Bin, RoundingMode::None).unwrap());

            eps.set_exponent(d1.exponent() - prec.min(p1) as Exponent + 2);

            assert!(
                d1.sub(&d3, prec, RoundingMode::ToEven)?.abs()?.cmp(&eps) < 0,
                "{} {:?}",
                prec,
                d1
            );
        }
    }

    // argument between -pi, -pi/2 and between pi/2, pi
    for _ in 0..TEST_ITERS {
        let p1 = (crate::common::test_rng::random::<usize>() % prec_rng + 1) * WORD_BIT_SIZE;
        let prec = (crate::common::test_rng::random::<usize>() % prec_rng + 1) * WORD_BIT_SIZE;

        let mut d1 = ExactNumNumber::random_normal(p1, -100, 2)?;

        // -pi, -pi/2 and pi/2, pi
        while d1.abs()?.cmp(&half_pi) < 0 {
            if d1.is_positive() {
                d1 = d1.add(&half_pi, p1, RoundingMode::None)?;
            }
            if d1.is_negative() {
                d1 = d1.sub(&half_pi, p1, RoundingMode::None)?;
            }
        }

        let arg = if d1.is_positive() {
            pi.sub(&d1, p1, RoundingMode::ToEven)?
        } else {
            pi.add(&d1, p1, RoundingMode::ToEven)?
        };

        let d2 = arg.sin(prec, RoundingMode::ToEven, &mut cc)?;
        let d3 = d2.asin(prec, RoundingMode::ToEven, &mut cc)?;

        // avoid values of sin close to 1 because of limited precision
        if ONE.sub(&d2.abs()?, prec, RoundingMode::None)?.cmp(&thres) >= 0 {
            // println!("{}", arg.format(Radix::Dec, RoundingMode::None).unwrap());
            // println!("{}", d1.format(Radix::Dec, RoundingMode::None).unwrap());
            // println!("{}", d2.format(Radix::Dec, RoundingMode::None).unwrap());
            // println!("{}", d3.format(Radix::Dec, RoundingMode::None).unwrap());

            eps.set_exponent(d1.exponent() - prec as Exponent - thres_exp);

            assert!(
                d1.abs()?
                    .sub(&d3, prec, RoundingMode::ToEven)?
                    .abs()?
                    .cmp(&eps)
                    < 0
                    || arg
                        .abs()?
                        .sub(&d3, prec, RoundingMode::ToEven)?
                        .abs()?
                        .cmp(&eps)
                        < 0,
                "{} {:?}",
                prec,
                d1
            );
        }
    }
    Ok(())
}

#[test]
fn test_cos_acos() -> Result<(), crate::Error> {
    let prec_rng = get_prec_rng();
    let mut eps = ExactNumNumber::small_const(1);

    let mut cc = Consts::new()?;

    let pi = cc.pi_num((prec_rng + 1) * WORD_BIT_SIZE, RoundingMode::None)?;

    for i in 0..TEST_ITERS {
        let p1 = (crate::common::test_rng::random::<usize>() % prec_rng + 1) * WORD_BIT_SIZE;
        let prec = (crate::common::test_rng::random::<usize>() % prec_rng + 1) * WORD_BIT_SIZE;

        if i & 1 == 0 {
            let mut d1 = ExactNumNumber::random_normal(p1, -(prec as Exponent), 3)?;

            // -pi, pi
            while d1.abs()?.cmp(&pi) > 0 {
                if d1.is_positive() {
                    d1 = d1.sub(&pi, p1, RoundingMode::None)?;
                }
                if d1.is_negative() {
                    d1 = d1.add(&pi, p1, RoundingMode::None)?;
                }
            }

            let d2 = d1.cos(prec, RoundingMode::ToEven, &mut cc)?;
            let d3 = d2.acos(prec, RoundingMode::ToEven, &mut cc)?;

            // println!("d1 {}", d1.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("d1 {:?}", d1);
            // println!("d2 {}", d2.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("d3 {}", d3.format(crate::Radix::Bin, RoundingMode::None).unwrap());

            if d2.cmp(&ONE) != 0 {
                eps.set_exponent(
                    d1.exponent() - prec as Exponent
                        + 1
                        + count_leading_ones(d2.mantissa().digits()) as Exponent,
                );

                assert!(
                    d1.abs()?
                        .sub(&d3, prec, RoundingMode::ToEven)?
                        .abs()?
                        .cmp(&eps)
                        < 0,
                    "{} {:?}",
                    prec,
                    d1
                );
            }
        } else {
            let d1 = ExactNumNumber::random_normal(p1, -(prec as Exponent), 0)?;

            let d2 = d1.acos(prec, RoundingMode::ToEven, &mut cc)?;

            let mut hp = cc.pi_num(prec, RoundingMode::ToZero)?;
            hp.set_exponent(1);

            if d2.abs_cmp(&hp) < 0 {
                let d3 = d2.cos(prec, RoundingMode::ToEven, &mut cc)?;

                eps.set_exponent(d1.exponent() - prec.min(p1) as Exponent - d1.exponent() + 2);

                // println!("d1 {}", d1.format(crate::Radix::Bin, RoundingMode::None).unwrap());
                // println!("d2 {}", d2.format(crate::Radix::Bin, RoundingMode::None).unwrap());
                // println!("d3 {}", d3.format(crate::Radix::Bin, RoundingMode::None).unwrap());

                assert!(
                    d1.abs()?
                        .sub(&d3, prec, RoundingMode::ToEven)?
                        .abs()?
                        .cmp(&eps)
                        < 0,
                    "{} {:?}",
                    prec,
                    d1
                );
            }
        }
    }
    Ok(())
}

#[test]
fn test_tan_atan() -> Result<(), crate::Error> {
    let prec_rng = get_prec_rng();
    let mut eps = ExactNumNumber::small_const(1);

    let mut cc = Consts::new()?;

    let mut half_pi = cc.pi_num((prec_rng + 1) * WORD_BIT_SIZE, RoundingMode::None)?;
    half_pi.set_exponent(1);

    for i in 0..TEST_ITERS {
        let p1 = (crate::common::test_rng::random::<usize>() % prec_rng + 1) * WORD_BIT_SIZE;
        let prec = (crate::common::test_rng::random::<usize>() % prec_rng + 1) * WORD_BIT_SIZE;

        if i & 1 == 0 {
            let mut d1 = ExactNumNumber::random_normal(p1, -100, 2)?;

            // -pi/2, pi/2
            while d1.abs()?.cmp(&half_pi) > 0 {
                if d1.is_positive() {
                    d1 = d1.sub(&half_pi, p1, RoundingMode::None)?;
                }
                if d1.is_negative() {
                    d1 = d1.add(&half_pi, p1, RoundingMode::None)?;
                }
            }

            let d2 = d1.tan(prec, RoundingMode::ToEven, &mut cc)?;
            let d3 = d2.atan(prec, RoundingMode::ToEven, &mut cc)?;

            // println!("d1 {}", d1.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("d2 {}", d2.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("d3 {}", d3.format(crate::Radix::Bin, RoundingMode::None).unwrap());

            eps.set_exponent(d1.exponent() - prec as Exponent + 2);

            assert!(
                d1.sub(&d3, prec, RoundingMode::ToEven)?.abs()?.cmp(&eps) < 0,
                "{} {:?}",
                prec,
                d1
            );
        } else {
            let d1 = ExactNumNumber::random_normal(p1, -TEST_EXP_BOUND, TEST_EXP_BOUND)?;

            let d2 = d1.atan(prec, RoundingMode::ToZero, &mut cc)?;

            let mut hp = cc.pi_num(prec, RoundingMode::ToZero)?;
            hp.set_exponent(1);

            if d2.abs_cmp(&hp) < 0 {
                let d3 = d2.tan(prec, RoundingMode::ToEven, &mut cc)?;

                // Large |x|: tan(atan(x)) error scales like x^2 * ulp(atan(x)).
                let err_exp = if d1.exponent() > 0 {
                    d1.exponent() * 2 - prec as Exponent + 2
                } else {
                    d1.exponent() - prec as Exponent + 2
                };
                eps.set_exponent(err_exp);

                // println!("d1 {}", d1.format(crate::Radix::Bin, RoundingMode::None).unwrap());
                // println!("d2 {}", d2.format(crate::Radix::Bin, RoundingMode::None).unwrap());
                // println!("d3 {}", d3.format(crate::Radix::Bin, RoundingMode::None).unwrap());

                assert!(
                    d1.sub(&d3, prec, RoundingMode::ToEven)?.abs()?.cmp(&eps) < 0,
                    "{} {:?}",
                    prec,
                    d1
                );
            }
        }
    }
    Ok(())
}

#[test]
fn test_sinh_asinh() -> Result<(), crate::Error> {
    let prec_rng = get_prec_rng();
    let mut eps = ExactNumNumber::small_const(1);

    let mut cc = Consts::new()?;

    for i in 0..TEST_ITERS {
        let p1 = (crate::common::test_rng::random::<usize>() % prec_rng + 1) * WORD_BIT_SIZE;
        let prec = (crate::common::test_rng::random::<usize>() % prec_rng + 1) * WORD_BIT_SIZE;

        if i & 1 == 0 {
            let d1 = ExactNumNumber::random_normal(p1, -100, 10)?;

            let d2 = d1.sinh(prec, RoundingMode::ToEven, &mut cc)?;
            let d3 = d2.asinh(prec, RoundingMode::ToEven, &mut cc)?;

            // println!("d1 {}", d1.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("d2 {}", d2.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("d3 {}", d3.format(crate::Radix::Bin, RoundingMode::None).unwrap());

            eps.set_exponent(d1.exponent() - prec as Exponent + 2);

            assert!(
                d1.sub(&d3, prec, RoundingMode::ToEven)?.abs()?.cmp(&eps) < 0,
                "{} {:?}",
                prec,
                d1
            );
        } else {
            let mut d1 = ExactNumNumber::random_normal(p1, -TEST_EXP_BOUND, TEST_EXP_BOUND)?;

            let d2 = d1.asinh(prec, RoundingMode::ToEven, &mut cc)?;
            let d3 = d2.sinh(prec, RoundingMode::ToEven, &mut cc)?;

            if d1.exponent() < -(prec as Exponent) {
                // both asinh and sinh are linear near 0
                d1.set_precision(prec, RoundingMode::ToEven)?;

                assert!(d1.cmp(&d2) == 0, "{} {:?}", prec, d1);
                assert!(d2.cmp(&d3) == 0, "{} {:?}", prec, d1);
            } else {
                // |sinh(asinh(x)) - x| ~ |asinh(x)| * |x| * 2^(-p) when asinh has p relative bits.
                eps.set_exponent(
                    d1.exponent() - prec as Exponent + d2.exponent().unsigned_abs() as Exponent + 2,
                );

                assert!(
                    d1.sub(&d3, prec, RoundingMode::ToEven)?.abs()?.cmp(&eps) < 0,
                    "{} {:?}",
                    prec,
                    d1
                );
            }
        }
    }
    Ok(())
}

#[test]
fn test_cosh_acosh() -> Result<(), crate::Error> {
    let prec_rng = get_prec_rng();
    let mut eps = ExactNumNumber::small_const(1);

    let mut cc = Consts::new()?;

    for i in 0..TEST_ITERS {
        let p1 = (crate::common::test_rng::random::<usize>() % prec_rng + 1) * WORD_BIT_SIZE;
        let prec = (crate::common::test_rng::random::<usize>() % prec_rng + 1) * WORD_BIT_SIZE;

        if i & 1 == 0 {
            let d1 = ExactNumNumber::random_normal(p1, -100, 10)?;

            let d2 = d1.cosh(prec, RoundingMode::ToEven, &mut cc)?;
            let d3 = d2.acosh(prec, RoundingMode::ToEven, &mut cc)?;

            eps.set_exponent(
                d1.exponent() - prec as Exponent
                    + 2
                    + count_leading_zeroes_skip_first(d2.mantissa().digits()) as Exponent,
            );

            assert!(
                d1.abs()?
                    .sub(&d3, prec, RoundingMode::ToEven)?
                    .abs()?
                    .cmp(&eps)
                    < 0,
                "{} {:?}",
                prec,
                d1
            );
        } else {
            let mut d1 = ExactNumNumber::random_normal(p1, 1, TEST_EXP_BOUND)?;
            d1.set_sign(Sign::Pos);

            let d2 = d1.acosh(prec, RoundingMode::ToEven, &mut cc)?;
            let d3 = d2.cosh(prec, RoundingMode::ToEven, &mut cc)?;

            let exp = d2.exp(prec + 1, RoundingMode::ToEven, &mut cc)?;
            let expr = exp.reciprocal(prec + 1, RoundingMode::ToEven)?;
            let mut d4 = exp.add(&expr, prec + 1, RoundingMode::ToEven)?;
            d4.set_exponent(d4.exponent() - 1);

            d4.set_precision(prec, RoundingMode::ToEven)?;

            // println!("d1 {}", d1.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("d2 {}", d2.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("d3 {}", d3.format(crate::Radix::Bin, RoundingMode::None).unwrap());
            // println!("d4 {}", d4.format(crate::Radix::Bin, RoundingMode::None).unwrap());

            assert!(d3.cmp(&d4) == 0, "{} {:?}", prec, d1);
        };
    }
    Ok(())
}

#[test]
fn test_tanh_atanh() -> Result<(), crate::Error> {
    let prec_rng = get_prec_rng();
    let mut eps = ExactNumNumber::small_const(1);

    let mut cc = Consts::new()?;

    let exp_to;
    #[cfg(not(target_pointer_width = "32"))]
    {
        exp_to = 5;
    }
    #[cfg(target_pointer_width = "32")]
    {
        exp_to = 3;
    }

    for i in 0..TEST_ITERS {
        let p1 = (crate::common::test_rng::random::<usize>() % prec_rng + 1) * WORD_BIT_SIZE;
        let prec = (crate::common::test_rng::random::<usize>() % prec_rng + 1) * WORD_BIT_SIZE;

        let (d1, d3) = if i & 1 == 0 {
            let d1 = ExactNumNumber::random_normal(p1, -100, exp_to)?;

            let d2 = d1.tanh(prec, RoundingMode::ToEven, &mut cc)?;
            let d3 = d2.atanh(prec, RoundingMode::ToEven, &mut cc)?;

            eps.set_exponent(
                d1.exponent() - prec as Exponent
                    + 1
                    + count_leading_ones(d2.mantissa().digits()) as Exponent,
            );

            (d1, d3)
        } else {
            let d1 = ExactNumNumber::random_normal(p1, EXPONENT_MIN, 0)?;

            let d2 = d1.atanh(prec, RoundingMode::ToEven, &mut cc)?;
            let d3 = d2.tanh(prec, RoundingMode::ToEven, &mut cc)?;

            eps.set_exponent(d1.exponent() - prec as Exponent + 1);

            (d1, d3)
        };

        // println!("d1 {}", d1.format(crate::Radix::Bin, RoundingMode::None).unwrap());
        // println!("d2 {}", d2.format(crate::Radix::Bin, RoundingMode::None).unwrap());
        // println!("d3 {}", d3.format(crate::Radix::Bin, RoundingMode::None).unwrap());
        // println!("e {}", eps.format(crate::Radix::Bin, RoundingMode::None).unwrap());

        assert!(
            d1.sub(&d3, prec, RoundingMode::ToEven)?.abs()?.cmp(&eps) < 0,
            "{} {:?}",
            prec,
            d1
        );
    }
    Ok(())
}
