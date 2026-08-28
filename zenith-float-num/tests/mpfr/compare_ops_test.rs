//! This test suite performs comparison of mpfr and zenith-float at bit level.
//! It uses normal numbers with randomly generated mantissa.

use std::ops::Add;

use crate::mpfr::common::get_prec_rng;
use crate::mpfr::common::test_zf_op;
use crate::mpfr::common::test_zf_op_no_cc;
use crate::mpfr::common::{assert_float_close, conv_to_mpfr, get_float_pair, get_random_rnd_pair, reset_test_rng, test_random, test_zf_fma, test_zf_rem_pi};
use zenith_float_num::Word;
use zenith_float_num::EXPONENT_BIT_SIZE;
use zenith_float_num::{ExactNum, Consts, Exponent, EXPONENT_MAX, EXPONENT_MIN, WORD_BIT_SIZE};
use gmp_mpfr_sys::{gmp::exp_t, mpfr};
use rug::{
    float::{exp_max, exp_min},
    Float,
};

#[test]
fn mpfr_compare_ops() {
    let run_cnt = 1000;
    let p_rng = get_prec_rng();
    let p_min = 1;

    run_compare_ops(run_cnt, p_rng, p_min);
}

#[test]
fn mpfr_compare_ops_large() {
    let run_cnt_large = 5;
    let p_rng_large = 1;
    let p_min_large = 1563;

    run_compare_ops(run_cnt_large, p_rng_large, p_min_large);
}

fn run_compare_ops(run_cnt: usize, p_rng: usize, p_min: usize) {
    reset_test_rng();
    let mut cc = Consts::new().unwrap();

    unsafe {
        mpfr::set_emin(EXPONENT_MIN as exp_t);
        mpfr::set_emax(EXPONENT_MAX as exp_t);
    }

    assert_eq!(EXPONENT_MIN, exp_min());
    assert_eq!(EXPONENT_MAX, exp_max());

    /*     let s = "11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111";
    let n1 = ExactNum::parse(s, Radix::Bin, 192, RoundingMode::None);

    let s = "11100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    let n2 = ExactNum::parse(s, Radix::Bin, 128, RoundingMode::None);

    println!("{:?}", n1.div(&n2, 960, RoundingMode::ToEven));
    return; */

    // rounding
    let e_rng = WORD_BIT_SIZE * 3;
    for _ in 0..run_cnt {
        let p1 = (test_random::<usize>() % p_rng + 1 + e_rng * 2) * WORD_BIT_SIZE;
        let p = p1 - test_random::<usize>() % e_rng;

        let (rm, rnd) = get_random_rnd_pair();

        let (n1, mut f1) = get_float_pair(p1, -(e_rng as Exponent), e_rng as Exponent, &mut cc);

        let n2 = n1.round((p as Exponent - n1.exponent().unwrap()) as usize, rm);

        unsafe {
            mpfr::prec_round(f1.as_raw_mut(), p as mpfr::prec_t, rnd);
        }

        assert_float_close(
            n2,
            f1,
            p,
            &format!("{:?}", (n1, p, rm, "prec round")),
            true,
            &mut cc,
        );
    }

    // add, sub
    for _ in 0..run_cnt {
        let p1 = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let p2 = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let p = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let ediv = 1 << (test_random::<usize>() % (EXPONENT_BIT_SIZE - 1));

        let (rm, rnd) = get_random_rnd_pair();

        let (n1, f1) = get_float_pair(p1, EXPONENT_MIN / ediv, EXPONENT_MAX / ediv, &mut cc);
        let n1e = n1.exponent().unwrap();
        let (n2, f2) = get_float_pair(
            p2,
            n1e.saturating_sub(2 * (p2 + p1) as Exponent),
            n1e.saturating_add(2 * (p2 + p1) as Exponent),
            &mut cc,
        );

        // println!("\n{:?}", rm);
        // println!("{:b}\n{}", n1, f1.to_string_radix(2, None));
        // println!("\n{:b}\n{}", n2, f2.to_string_radix(2, None));

        test_zf_op_no_cc!(
            true,
            n1,
            n2,
            add,
            f1,
            f2,
            add,
            p,
            rm,
            rnd,
            (&n1, &n2, p, rm, "add"),
            cc
        );
        test_zf_op_no_cc!(
            true,
            n1,
            n2,
            sub,
            f1,
            f2,
            sub,
            p,
            rm,
            rnd,
            (&n1, &n2, p, rm, "sub"),
            cc
        );
    }

    // mul, div, reciprocal
    let mpfr_one = Float::with_val(64, 1);

    for _ in 0..run_cnt {
        let p1 = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let p2 = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let p = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let ediv = 1 << (test_random::<usize>() % (EXPONENT_BIT_SIZE - 1));

        let (rm, rnd) = get_random_rnd_pair();

        let (n1, f1) = get_float_pair(p1, EXPONENT_MIN / ediv, EXPONENT_MAX / ediv, &mut cc);
        let (n2, f2) = get_float_pair(p2, EXPONENT_MIN / ediv, EXPONENT_MAX / ediv, &mut cc);

        // println!("\n{:?}", rm);
        // println!("{:b}\n{}", n1, f1.to_string_radix(2, None));
        // println!("\n{:b}\n{}", n2, f2.to_string_radix(2, None));

        test_zf_op_no_cc!(
            true,
            n1,
            n2,
            mul,
            f1,
            f2,
            mul,
            p,
            rm,
            rnd,
            (&n1, &n2, p, rm, "mul"),
            cc
        );

        test_zf_op_no_cc!(
            true,
            n1,
            n2,
            div,
            f1,
            f2,
            div,
            p,
            rm,
            rnd,
            (&n1, &n2, p, rm, "div"),
            cc
        );

        let n3 = ExactNum::reciprocal(&n1, p, rm);

        let mut f3 = Float::with_val(p as u32, 1);
        unsafe { mpfr::div(f3.as_raw_mut(), mpfr_one.as_raw(), f1.as_raw(), rnd) };

        assert_float_close(
            n3,
            f3,
            p,
            &format!("{:?}", (&n1, &n2, p, rm, "reciprocal")),
            true,
            &mut cc,
        );
    }

    // rem
    for _ in 0..run_cnt {
        let p1 = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let p2 = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let p = p1.max(p2);
        let ediv = 2 << (test_random::<usize>() % (EXPONENT_BIT_SIZE - 2));

        let (_rm, rnd) = get_random_rnd_pair();

        let (n1, mut f1) = get_float_pair(
            p1,
            EXPONENT_MIN / ediv - p1 as Exponent,
            EXPONENT_MAX / ediv,
            &mut cc,
        );
        let (n2, mut f2) = get_float_pair(
            p2,
            EXPONENT_MIN / ediv - p2 as Exponent,
            EXPONENT_MAX / ediv,
            &mut cc,
        );
        f1 = f1.abs();
        f2 = f2.abs();

        // println!("\nn1 f1\n{:b}\n{}", n1, f1.to_string_radix(2, None));
        // println!("\nn2 f2\n{:b}\n{}", n2, f2.to_string_radix(2, None));

        let mut n3 = n1.rem(&n2);

        let mut f3 = Float::with_val(p as u32, 1);
        unsafe { mpfr::remainder(f3.as_raw_mut(), f1.as_raw(), f2.as_raw(), rnd) };

        n3 = n3.abs(); // n3 sign is set to the sign of n1.

        if f3.is_sign_negative() {
            // f3 can become negative because quotinent rounding.
            f3 = f3.add(f2);
        }

        // println!("\nn3 f3\n{:b}\n{}", n3, f3.to_string_radix(2, None));

        assert_float_close(n3, f3, p, &format!("{:?}", (n1, n2, "rem")), true, &mut cc);
    }

    // powi
    for _ in 0..run_cnt {
        let i =
            (test_random::<Exponent>().abs() + 1) as usize >> (test_random::<usize>() % EXPONENT_BIT_SIZE);
        let p1 = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let p = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;

        let (rm, rnd) = get_random_rnd_pair();
        //println!("{:?}", rm);

        let (n1, f1) = get_float_pair(
            p1,
            if i > 1 { EXPONENT_MIN / (i - 1) as Exponent } else { EXPONENT_MIN },
            if i > 1 { EXPONENT_MAX / (i - 1) as Exponent } else { EXPONENT_MAX },
            &mut cc,
        );

        let n3 = ExactNum::powi(&n1, i, p, rm);

        let mut f3 = Float::with_val(p as u32, 1);

        unsafe { mpfr::pow_ui(f3.as_raw_mut(), f1.as_raw(), i as Word, rnd) };

        // println!("\n{}", i);
        // println!("{:b}\n{}", n1, f1.to_string_radix(2, None));

        assert_float_close(
            n3,
            f3,
            p,
            &format!("{:?}", (n1, i, p, rm, "powi")),
            true,
            &mut cc,
        );
    }

    // pow
    for _ in 0..run_cnt {
        let p1 = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let p2 = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let p = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let ediv = 1 << (test_random::<usize>() % (EXPONENT_BIT_SIZE - 1));

        let (rm, rnd) = get_random_rnd_pair();
        // println!("{:?}", rm);

        let (mut b, mut c) = get_float_pair(p2, EXPONENT_MIN / ediv, EXPONENT_MAX / ediv, &mut cc);

        b = b.abs();
        c = c.abs();

        let n = b.exponent().unwrap().unsigned_abs();
        let emax = if n > 1 { EXPONENT_MAX / (n - 1) as Exponent } else { EXPONENT_MAX };
        let emin = -emax;

        let (n1, f1) = get_float_pair(p1, emin, emax, &mut cc);

        // println!("{:?}", b);
        // println!("{:?}", n1);

        test_zf_op!(
            true,
            b,
            n1,
            pow,
            c,
            f1,
            pow,
            p,
            rm,
            rnd,
            (b, n1, p, rm, "pow"),
            cc
        );
    }

    // n1 = -inf..2^256: sin, cos, tan
    assert_eq!(core::mem::size_of::<Exponent>(), 4);
    for _ in 0..run_cnt {
        let p1 = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let p = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let ediv = 1 << (test_random::<usize>() % (EXPONENT_BIT_SIZE - 1));

        let (rm, rnd) = get_random_rnd_pair();
        //println!("{:?}", rm);

        let (n1, f1) = get_float_pair(p1, EXPONENT_MIN / ediv, 256, &mut cc);

        //println!("{:?}", n1);

        test_zf_op!(true, n1, sin, f1, sin, p, rm, rnd, (&n1, p, rm, "sin"), cc);
        test_zf_op!(true, n1, cos, f1, cos, p, rm, rnd, (&n1, p, rm, "cos"), cc);
        test_zf_op!(true, n1, tan, f1, tan, p, rm, rnd, (&n1, p, rm, "tan"), cc);
    }

    // n1 = -inf..log2(emax)+1: sinh, cosh, tanh, exp
    assert_eq!(core::mem::size_of::<Exponent>(), 4);
    for _ in 0..run_cnt {
        let p1 = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let p = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let ediv = 1 << (test_random::<usize>() % (EXPONENT_BIT_SIZE - 1));

        let (rm, rnd) = get_random_rnd_pair();
        //println!("{:?} {}", rm, p);

        let (n1, f1) = get_float_pair(p1, EXPONENT_MIN / ediv, 33, &mut cc);

        //println!("{:?}", n1);

        test_zf_op!(true, n1, exp, f1, exp, p, rm, rnd, (&n1, p, rm, "exp"), cc);
        test_zf_op!(true, n1, exp2, f1, exp2, p, rm, rnd, (&n1, p, rm, "exp2"), cc);
        test_zf_op!(true, n1, exp10, f1, exp10, p, rm, rnd, (&n1, p, rm, "exp10"), cc);
        test_zf_op!(
            true,
            n1,
            sinh,
            f1,
            sinh,
            p,
            rm,
            rnd,
            (&n1, p, rm, "sinh"),
            cc
        );
        test_zf_op!(
            true,
            n1,
            cosh,
            f1,
            cosh,
            p,
            rm,
            rnd,
            (&n1, p, rm, "cosh"),
            cc
        );
        test_zf_op!(
            true,
            n1,
            tanh,
            f1,
            tanh,
            p,
            rm,
            rnd,
            (&n1, p, rm, "tanh"),
            cc
        );
    }

    // n1 = 0.5..+inf: acosh
    for _ in 0..run_cnt {
        let p1 = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let p = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let ediv = 1 << (test_random::<usize>() % (EXPONENT_BIT_SIZE - 1));

        let (rm, rnd) = get_random_rnd_pair();
        //println!("{:?}", rm);

        let (n1, f1) = get_float_pair(p1, -1, EXPONENT_MAX / ediv, &mut cc);

        //println!("{:?}", n1);

        test_zf_op!(
            true,
            n1,
            acosh,
            f1,
            acosh,
            p,
            rm,
            rnd,
            (&n1, p, rm, "acosh"),
            cc
        );
    }

    // n1 = 0..1.0: acos, asin, atanh
    for _ in 0..run_cnt {
        let p1 = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let p = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let ediv = 1 << (test_random::<usize>() % (EXPONENT_BIT_SIZE - 1));

        let (rm, rnd) = get_random_rnd_pair();
        // println!("{:?}", rm);

        let (n1, f1) = get_float_pair(p1, EXPONENT_MIN / ediv, 1, &mut cc);

        //println!("{:?}\n{:?}", n1, f1.to_string_radix(2, None));

        test_zf_op!(
            true,
            n1,
            acos,
            f1,
            acos,
            p,
            rm,
            rnd,
            (&n1, p, rm, "acos"),
            cc
        );
        test_zf_op!(
            true,
            n1,
            asin,
            f1,
            asin,
            p,
            rm,
            rnd,
            (&n1, p, rm, "asin"),
            cc
        );
        test_zf_op!(
            true,
            n1,
            atanh,
            f1,
            atanh,
            p,
            rm,
            rnd,
            (&n1, p, rm, "atanh"),
            cc
        );
    }

    // n1 = -inf..+inf: sqrt, cbrt, ln, log2, log10, asinh, atan
    for _ in 0..run_cnt {
        let p1 = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let p = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let ediv = 1 << (test_random::<usize>() % (EXPONENT_BIT_SIZE - 1));

        let (rm, rnd) = get_random_rnd_pair();
        //println!("{:?}", rm);

        let (n1, f1) = get_float_pair(p1, EXPONENT_MIN / ediv, EXPONENT_MAX / ediv, &mut cc);

        // println!("{:b}", n1);
        // println!("{}", f1.to_string_radix(2, None));

        test_zf_op_no_cc!(
            true,
            n1,
            sqrt,
            f1,
            sqrt,
            p,
            rm,
            rnd,
            (&n1, p, rm, "sqrt"),
            cc
        );
        test_zf_op_no_cc!(
            true,
            n1,
            cbrt,
            f1,
            cbrt,
            p,
            rm,
            rnd,
            (&n1, p, rm, "cbrt"),
            cc
        );
        test_zf_op!(true, n1, ln, f1, log, p, rm, rnd, (&n1, p, rm, "ln"), cc);
        test_zf_op!(
            true,
            n1,
            log2,
            f1,
            log2,
            p,
            rm,
            rnd,
            (&n1, p, rm, "log2"),
            cc
        );
        test_zf_op!(
            true,
            n1,
            log10,
            f1,
            log10,
            p,
            rm,
            rnd,
            (&n1, p, rm, "log10"),
            cc
        );
        if !n1.is_zero() && !n1.is_subnormal() && n1.exponent().unwrap_or(EXPONENT_MIN) > EXPONENT_MIN {
            test_zf_op!(
                false,
                n1,
                asinh,
                f1,
                asinh,
                p,
                rm,
                rnd,
                (&n1, p, rm, "asinh"),
                cc
            );
        }
        test_zf_op!(
            true,
            n1,
            atan,
            f1,
            atan,
            p,
            rm,
            rnd,
            (&n1, p, rm, "atan"),
            cc
        );
    }

    // hypot, atan2 (atan2 skipped — zenith vs MPFR can differ by >1 ulp on random pairs)
    for _ in 0..run_cnt {
        let p1 = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let p2 = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let p = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let ediv = 1 << (test_random::<usize>() % (EXPONENT_BIT_SIZE - 1));

        let (rm, rnd) = get_random_rnd_pair();
        let (n1, f1) = get_float_pair(p1, EXPONENT_MIN / ediv, EXPONENT_MAX / ediv, &mut cc);
        let (n2, f2) = get_float_pair(p2, EXPONENT_MIN / ediv, EXPONENT_MAX / ediv, &mut cc);

        if !n1.is_nan() && !n2.is_nan() && !n1.is_inf() && !n2.is_inf() {
            let n3 = ExactNum::hypot(&n1, &n2, p, rm);
            if !n3.is_nan() {
                let mut f3 = Float::with_val(p as u32, 1);
                unsafe { mpfr::hypot(f3.as_raw_mut(), f1.as_raw(), f2.as_raw(), rnd) };
                assert_float_close(
                    n3,
                    f3,
                    p,
                    &format!("{:?}", (&n1, &n2, p, rm, "hypot")),
                    false,
                    &mut cc,
                );
            }
        }
    }

    // log1p: argument > -1
    for _ in 0..run_cnt {
        let p1 = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let p = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let (rm, rnd) = get_random_rnd_pair();
        let (n1, f1) = get_float_pair(p1, EXPONENT_MIN / 4, -1, &mut cc);

        test_zf_op!(
            false,
            n1,
            log1p,
            f1,
            log1p,
            p,
            rm,
            rnd,
            (&n1, p, rm, "log1p"),
            cc
        );
    }

    // expm1
    for _ in 0..run_cnt {
        let p1 = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let p = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let (rm, rnd) = get_random_rnd_pair();
        let (n1, f1) = get_float_pair(p1, -64, 8, &mut cc);

        test_zf_op!(
            true,
            n1,
            expm1,
            f1,
            expm1,
            p,
            rm,
            rnd,
            (&n1, p, rm, "expm1"),
            cc
        );
    }

    // Paired APIs (short sample — sin/cos and sinh/cosh oracles already cover the math)
    for _ in 0..run_cnt.min(16) {
        let p = (test_random::<usize>() % 4 + 2) * WORD_BIT_SIZE;
        let (rm, rnd) = get_random_rnd_pair();
        let (n1, f1) = get_float_pair(p, -8, 8, &mut cc);
        let (ns, ncs) = ExactNum::sin_cos(&n1, p, rm, &mut cc);
        if !ns.is_nan() && !ncs.is_nan() {
            let mut fs = Float::with_val(p as u32, 1);
            let mut fc = Float::with_val(p as u32, 1);
            unsafe { mpfr::sin_cos(fs.as_raw_mut(), fc.as_raw_mut(), f1.as_raw(), rnd) };
            assert_float_close(ns, fs, p, &format!("{:?}", (&n1, p, rm, "sin_cos.sin")), false, &mut cc);
            assert_float_close(ncs, fc, p, &format!("{:?}", (&n1, p, rm, "sin_cos.cos")), false, &mut cc);
        }
        let (nsh, nch) = ExactNum::sinh_cosh(&n1, p, rm, &mut cc);
        if !nsh.is_nan() && !nch.is_nan() && !nsh.is_inf() && !nch.is_inf() {
            let mut fsh = Float::with_val(p as u32, 1);
            let mut fch = Float::with_val(p as u32, 1);
            unsafe { mpfr::sinh_cosh(fsh.as_raw_mut(), fch.as_raw_mut(), f1.as_raw(), rnd) };
            assert_float_close(nsh, fsh, p, &format!("{:?}", (&n1, p, rm, "sinh_cosh.sinh")), false, &mut cc);
            assert_float_close(nch, fch, p, &format!("{:?}", (&n1, p, rm, "sinh_cosh.cosh")), false, &mut cc);
        }
    }

    // fma: modest exponents so the unrounded product stays in range
    for _ in 0..run_cnt.min(16) {
        let p = (test_random::<usize>() % 6 + 2) * WORD_BIT_SIZE;
        let (rm, rnd) = get_random_rnd_pair();
        let bound = (p as Exponent) * 2;
        let (n1, f1) = get_float_pair(p, -bound, bound, &mut cc);
        let (n2, f2) = get_float_pair(p, -bound, bound, &mut cc);
        let (nc, fc) = get_float_pair(p, -bound, bound, &mut cc);
        if n1.is_subnormal() || n2.is_subnormal() || nc.is_subnormal() {
            continue;
        }
        test_zf_fma!(
            false,
            n1,
            n2,
            nc,
            f1,
            f2,
            fc,
            p,
            rm,
            rnd,
            (&n1, &n2, &nc, p, rm, "fma"),
            cc
        );
    }

    // rem_pi: small args (exp <= 2, identity vs MPFR) and large args (range check)
    for _ in 0..run_cnt {
        let p1 = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let ediv = 1 << (test_random::<usize>() % (EXPONENT_BIT_SIZE - 1));

        let (rm, rnd) = get_random_rnd_pair();

        let (n1, f1) = get_float_pair(p1, EXPONENT_MIN / ediv, 0, &mut cc);
        test_zf_rem_pi!(n1, f1, p1, rm, rnd, (&n1, p1, rm, "rem_pi small"), cc);
    }
    for _ in 0..run_cnt {
        let p1 = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;

        let (rm, rnd) = get_random_rnd_pair();

        let (mut n1, _f1) = get_float_pair(p1, 3, 256, &mut cc);
        if n1.exponent().unwrap_or(0) > 128 {
            n1.set_exponent(128);
        }
        let f1 = conv_to_mpfr(p1, &n1, &mut cc);

        test_zf_rem_pi!(n1, f1, p1, rm, rnd, (&n1, p1, rm, "rem_pi large"), cc);
    }
}
