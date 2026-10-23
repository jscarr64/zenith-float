// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! This test suite performs comparison of mpfr and zenith-float at bit level.
//! It uses special cases of numbers, like zero, one, maximum possible value, a number with many trailing zeroes in mantissa, etc.

use crate::mpfr::common::{RAND_EMAX, RAND_EMIN};
use std::ops::Add;
use std::ops::Neg;

use crate::mpfr::common::{
    assert_float_close, conv_to_mpfr, get_last_zero, get_near_one, get_oned_sides, get_oned_zeroed,
    get_periodic, get_random_rnd_pair, reset_test_rng, test_random, test_zf_op_no_cc,
    test_zf_rem_pi,
};
use crate::mpfr::common::{draw_seed, get_prec_rng, rebuild_from_seed, test_zf_op};
use gmp_mpfr_sys::{gmp::exp_t, mpfr};
use rug::Float;
use zenith_float_num::{
    Consts, ExactNum, Exponent, EXPONENT_BIT_SIZE, EXPONENT_MAX, EXPONENT_MIN, WORD_BIT_SIZE,
};

#[test]
fn mpfr_compare_special() -> Result<(), Box<dyn std::error::Error>> {
    let run_cnt = crate::mpfr::common::tier_runs(500);
    let p_rng = get_prec_rng();
    let p_min = 1;

    run_compare_special(run_cnt, p_rng, p_min)?;
    Ok(())
}

#[test]
fn mpfr_compare_special_large() -> Result<(), Box<dyn std::error::Error>> {
    let run_cnt_large = crate::mpfr::common::tier_runs(3);
    let p_rng_large = 1;
    let p_min_large = 1563;

    run_compare_special(run_cnt_large, p_rng_large, p_min_large)?;
    Ok(())
}

fn run_compare_special(
    run_cnt: usize,
    p_rng: usize,
    p_min: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    reset_test_rng();
    let mut cc = Consts::new()?;

    unsafe {
        mpfr::set_emin(EXPONENT_MIN as exp_t);
        mpfr::set_emax(EXPONENT_MAX as exp_t);
    }

    assert_eq!(EXPONENT_MIN, unsafe { mpfr::get_emin() });
    assert_eq!(EXPONENT_MAX, unsafe { mpfr::get_emax() });

    let mpfr_one = Float::with_val(64, 1);

    for _ in 0..run_cnt {
        let p1 = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let p2 = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let p = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let ediv = 1 << (test_random::<usize>() % (EXPONENT_BIT_SIZE - 1));
        let emin = RAND_EMIN / ediv;
        let emax = RAND_EMAX / ediv;

        let (rm, rnd) = get_random_rnd_pair();

        // `n` is consumed by add / sub / mul / div below: each list element has its own seed and is rebuilt from
        // it (one element, not the whole list) for every use.
        let mk_nth = |idx: usize, seed: u64| {
            rebuild_from_seed(seed, || match idx {
                0 => ExactNum::new(p1),
                1 => ExactNum::new(p1).neg(),
                2 => ExactNum::from_i8(1, p1),
                3 => ExactNum::from_i8(-1, p1),
                4 => ExactNum::max_value(p1),
                5 => ExactNum::min_value(p1),
                6 => ExactNum::min_positive_normal(p1),
                7 => ExactNum::min_positive_normal(p1).neg(),
                8 => get_oned_zeroed(p1, emin, emax),
                9 => get_oned_sides(p1, emin, emax),
                10 => get_periodic(p1, emin, emax),
                11 => get_last_zero(p1, emin, emax),
                _ => get_near_one(p1),
            })
        };
        let seeds: [u64; 13] = std::array::from_fn(|_| draw_seed());
        let nn: [ExactNum; 13] = std::array::from_fn(|idx| mk_nth(idx, seeds[idx]));

        for (idx, n) in nn.iter().enumerate() {
            //println!("{:?}", n);

            let f = conv_to_mpfr(p1, n, &mut cc)?;

            let nn1 = [
                ExactNum::new(p2),
                ExactNum::new(p2).neg(),
                ExactNum::from_i8(1, p2),
                ExactNum::from_i8(-1, p2),
                ExactNum::max_value(p2),
                ExactNum::min_value(p2),
                ExactNum::min_positive_normal(p2),
                ExactNum::min_positive_normal(p2).neg(),
                get_oned_zeroed(p2, emin, emax),
                get_oned_sides(p2, emin, emax),
                get_periodic(p2, emin, emax),
                get_last_zero(p2, emin, emax),
                get_near_one(p2),
            ];

            for n1 in nn1.iter() {
                let f1 = conv_to_mpfr(p2, n1, &mut cc)?;

                //println!("rm {:?}", rm);
                //println!("\n--{:?}\n{:?}", n, n1);
                //println!("\n--{:?}\n{:?}", f, f1);
                test_zf_op_no_cc!(
                    true,
                    mk_nth(idx, seeds[idx]),
                    n1,
                    add,
                    f,
                    f1,
                    add,
                    p,
                    rm,
                    rnd,
                    (n, n1, p, rm, "add"),
                    cc
                );
                test_zf_op_no_cc!(
                    true,
                    mk_nth(idx, seeds[idx]),
                    n1,
                    sub,
                    f,
                    f1,
                    sub,
                    p,
                    rm,
                    rnd,
                    (n, n1, p, rm, "sub"),
                    cc
                );
                test_zf_op_no_cc!(
                    true,
                    mk_nth(idx, seeds[idx]),
                    n1,
                    mul,
                    f,
                    f1,
                    mul,
                    p,
                    rm,
                    rnd,
                    (n, n1, p, rm, "mul"),
                    cc
                );
                test_zf_op_no_cc!(
                    true,
                    mk_nth(idx, seeds[idx]),
                    n1,
                    div,
                    f,
                    f1,
                    div,
                    p,
                    rm,
                    rnd,
                    (n, n1, p, rm, "div"),
                    cc
                );

                test_zf_op!(
                    false,
                    n,
                    n1,
                    pow,
                    f,
                    f1,
                    pow,
                    p,
                    rm,
                    rnd,
                    (n, n1, p, rm, "pow"),
                    cc
                );

                // rem
                let p = p1.max(p2);
                let f = Float::with_val(f.prec(), f.abs_ref());
                let f1 = f1.abs();

                let mut n3 = n.rem(n1);

                let mut f3 = Float::with_val(p as u32, 0);
                unsafe { mpfr::remainder(f3.as_raw_mut(), f.as_raw(), f1.as_raw(), rnd) };

                n3 = n3.abs(); // n3 sign is set to the sign of n1.

                if f3.is_sign_negative() && !f3.is_zero() {
                    // f3 can become negative because of quotinent rounding.
                    f3 = f3.add(f1);
                }

                //println!("\n{:?}\n{:?}", n3, f3);

                assert_float_close(n3, f3, p, &format!("{:?}", (n, n1, "rem")), true, &mut cc)?;
            }

            test_zf_op_no_cc!(true, n, sqrt, f, sqrt, p, rm, rnd, (n, p, rm, "sqrt"), cc);
            test_zf_op_no_cc!(true, n, cbrt, f, cbrt, p, rm, rnd, (n, p, rm, "cbrt"), cc);
            test_zf_op!(true, n, ln, f, log, p, rm, rnd, (n, p, rm, "ln"), cc);
            test_zf_op!(true, n, log2, f, log2, p, rm, rnd, (n, p, rm, "log2"), cc);
            test_zf_op!(
                true,
                n,
                log10,
                f,
                log10,
                p,
                rm,
                rnd,
                (n, p, rm, "log10"),
                cc
            );
            // MPFR's small-argument shortcuts compute -2 EXP(x), which wraps for exponents near its
            // emin_min (= EXPONENT_MIN), and then size a buffer from it. Below the former 32-bit
            // band the unary functions that follow take their operand at the band's bottom instead;
            // i64_edges covers the 62-bit edges of the arithmetic.
            let n_low;
            let f_low;
            let (n, f) = if n.exponent().is_some_and(|e| e < RAND_EMIN) {
                let mut m = n.clone();
                m.set_exponent(RAND_EMIN);
                f_low = conv_to_mpfr(p1, &m, &mut cc)?;
                n_low = m;
                (&n_low, f_low)
            } else {
                (n, f)
            };
            if !n.is_zero()
                && !n.is_subnormal()
                && !n.is_inf()
                && !n.is_nan()
                && n.exponent().unwrap_or(EXPONENT_MIN) > EXPONENT_MIN
            {
                if n.exponent() == Some(EXPONENT_MAX) {
                    // MPFR's asinh loses ln 2 at its largest exponent (its 2x overflows):
                    // asinh(x) = sign(x) (ln 2 + ln |x|) there, to far below the last bit.
                    let r = if n.is_negative() {
                        let mirrored = match rnd {
                            mpfr::rnd_t::RNDU => mpfr::rnd_t::RNDD,
                            mpfr::rnd_t::RNDD => mpfr::rnd_t::RNDU,
                            other => other,
                        };
                        -big_log_ref(&Float::with_val(f.prec(), f.abs_ref()), p, mirrored)
                    } else {
                        big_log_ref(&f, p, rnd)
                    };
                    let n3 = n.asinh(p, rm, &mut cc);
                    assert_float_close(
                        n3,
                        r,
                        p,
                        &format!("{:?}", (&n, p, rm, "asinh max")),
                        false,
                        &mut cc,
                    )?;
                } else {
                    test_zf_op!(
                        false,
                        n,
                        asinh,
                        f,
                        asinh,
                        p,
                        rm,
                        rnd,
                        (n, p, rm, "asinh"),
                        cc
                    );
                }
            }
            test_zf_op!(true, n, atan, f, atan, p, rm, rnd, (n, p, rm, "atan"), cc);

            let mut n_trig = n.clone();
            let f_trig = if n.exponent().ok_or("none value")? > 128 {
                n_trig.set_exponent(128); // large exponent causes very long computation.
                conv_to_mpfr(p1, &n_trig, &mut cc)?
            } else {
                f.clone()
            };

            test_zf_op!(
                true,
                n_trig,
                sin,
                f_trig,
                sin,
                p,
                rm,
                rnd,
                (n, p, rm, "sin"),
                cc
            );
            test_zf_op!(
                true,
                n_trig,
                cos,
                f_trig,
                cos,
                p,
                rm,
                rnd,
                (n, p, rm, "cos"),
                cc
            );
            test_zf_op!(
                true,
                n_trig,
                tan,
                f_trig,
                tan,
                p,
                rm,
                rnd,
                (n, p, rm, "tan"),
                cc
            );
            test_zf_rem_pi!(n_trig, f_trig, p1, rm, rnd, (n, p1, rm, "rem_pi"), cc);

            test_zf_op!(true, n, exp, f, exp, p, rm, rnd, (n, p, rm, "exp"), cc);
            test_zf_op!(true, n, exp2, f, exp2, p, rm, rnd, (n, p, rm, "exp2"), cc);
            test_zf_op!(
                true,
                n,
                exp10,
                f,
                exp10,
                p,
                rm,
                rnd,
                (n, p, rm, "exp10"),
                cc
            );
            test_zf_op!(true, n, sinh, f, sinh, p, rm, rnd, (n, p, rm, "sinh"), cc);
            test_zf_op!(true, n, cosh, f, cosh, p, rm, rnd, (n, p, rm, "cosh"), cc);
            test_zf_op!(true, n, tanh, f, tanh, p, rm, rnd, (n, p, rm, "tanh"), cc);

            if n.exponent() == Some(EXPONENT_MAX) && n.is_positive() {
                // MPFR's acosh loses ln 2 at its largest exponent (its 2x overflows):
                // acosh(x) = ln 2 + ln x there, to far below the last bit.
                let r = big_log_ref(&f, p, rnd);
                let n3 = n.acosh(p, rm, &mut cc);
                assert_float_close(
                    n3,
                    r,
                    p,
                    &format!("{:?}", (&n, p, rm, "acosh max")),
                    true,
                    &mut cc,
                )?;
            } else {
                test_zf_op!(
                    true,
                    n,
                    acosh,
                    f,
                    acosh,
                    p,
                    rm,
                    rnd,
                    (n, p, rm, "acosh"),
                    cc
                );
            }

            test_zf_op!(true, n, acos, f, acos, p, rm, rnd, (n, p, rm, "acos"), cc);
            test_zf_op!(true, n, asin, f, asin, p, rm, rnd, (n, p, rm, "asin"), cc);
            test_zf_op!(
                true,
                n,
                atanh,
                f,
                atanh,
                p,
                rm,
                rnd,
                (n, p, rm, "atanh"),
                cc
            );

            // powi: covered by compare_ops; special values (max/min/zero) diverge from MPFR.

            // reciprocal (skip zero/inf/nan — MPFR oracle differs on non-finite cases)
            if !n.is_zero() && !n.is_inf() && !n.is_nan() {
                let n3 = ExactNum::reciprocal(n, p, rm);

                let mut f3 = Float::with_val(p as u32, 1);
                unsafe { mpfr::div(f3.as_raw_mut(), mpfr_one.as_raw(), f.as_raw(), rnd) };

                assert_float_close(
                    n3,
                    f3,
                    p,
                    &format!("{:?}", (n, p, rm, "reciprocal")),
                    false,
                    &mut cc,
                )?;
            }
        }
    }

    // grades of pi
    for _ in 0..run_cnt {
        let p1 = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        let p = (test_random::<usize>() % p_rng + p_min) * WORD_BIT_SIZE;
        // The former 32-bit band: MPFR's sin reduction sizes pi by the exponent, so 62-bit
        // exponents would ask it for exabytes.
        let ediv: Exponent = 1 << (test_random::<u32>() % 31);

        let (rm, rnd) = get_random_rnd_pair();
        let (rm2, _) = get_random_rnd_pair();

        let mut n = cc.pi(p1, rm2);
        let e = test_random::<u64>() % (10 + (RAND_EMAX / ediv) as u64);
        n.set_exponent(RAND_EMIN / ediv + e as Exponent);

        let f = conv_to_mpfr(p1, &n, &mut cc)?;

        test_zf_op!(true, n, sin, f, sin, p, rm, rnd, (&n, p, rm, "sin"), cc);

        test_zf_op!(true, n, cos, f, cos, p, rm, rnd, (&n, p, rm, "cos"), cc);

        test_zf_op!(true, n, tan, f, tan, p, rm, rnd, (n, p, rm, "tan"), cc);
    }
    Ok(())
}

/// `ln 2 + ln x` rounded to `p` bits by `rnd`, from a 256-bit-wider evaluation.
fn big_log_ref(f: &Float, p: usize, rnd: mpfr::rnd_t) -> Float {
    let wp = p as u32 + 256;
    let mut l = Float::with_val(wp, 0);
    let mut two = Float::with_val(wp, 2);
    unsafe {
        mpfr::log(l.as_raw_mut(), f.as_raw(), mpfr::rnd_t::RNDN);
        mpfr::log(two.as_raw_mut(), two.as_raw(), mpfr::rnd_t::RNDN);
    }
    let mut r = Float::with_val(p as u32, 0);
    unsafe { mpfr::add(r.as_raw_mut(), l.as_raw(), two.as_raw(), rnd) };
    r
}
