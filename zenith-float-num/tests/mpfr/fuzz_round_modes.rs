//! Short MPFR differential: add / sub / mul / div bit-exact in every mode including ToOdd,
//! sqrt and nth_root under the five MPFR modes, exact-bit precision sweeps, and
//! exponent overflow of add / sub / mul / div (max finite toward zero, Inf otherwise).

use crate::mpfr::common::assert_exact_eq;
use crate::mpfr::common::assert_float_close;
use crate::mpfr::common::conv_to_mpfr;
use crate::mpfr::common::get_float_pair;
use crate::mpfr::common::get_oned_zeroed;
use crate::mpfr::common::mpfr_reference;
use crate::mpfr::common::reset_test_rng;
use crate::mpfr::common::test_random;
use crate::mpfr::common::test_zf_op_no_cc;
use crate::mpfr::common::zenith_binop;
use crate::mpfr::common::BinOp;
use gmp_mpfr_sys::mpfr::{self, rnd_t};
use rug::float::{exp_max, exp_min};
use rug::Float;
use std::ops::Neg;
use zenith_float_num::{
    Consts, ExactNum, Exponent, RoundingMode, EXPONENT_MAX, EXPONENT_MIN, WORD_BIT_SIZE,
};

/// `-x` for an operand that stays borrowed: fresh `-1` times `x` at `x`'s own precision.
fn neg_of(x: &ExactNum) -> ExactNum {
    let p = x
        .mantissa_max_bit_len()
        .unwrap_or(zenith_float_num::WORD_BIT_SIZE);
    ExactNum::from_i8(-1, p).mul(x, p, RoundingMode::None)
}

const ROUNDS: [(RoundingMode, rnd_t); 5] = [
    (RoundingMode::ToEven, rnd_t::RNDN),
    (RoundingMode::Up, rnd_t::RNDU),
    (RoundingMode::Down, rnd_t::RNDD),
    (RoundingMode::FromZero, rnd_t::RNDA),
    (RoundingMode::ToZero, rnd_t::RNDZ),
];

/// Every mode with a bit-exact reference, `ToOdd` included (emulated in `mpfr_reference`).
const ALL_MODES: [RoundingMode; 6] = [
    RoundingMode::ToEven,
    RoundingMode::ToOdd,
    RoundingMode::Up,
    RoundingMode::Down,
    RoundingMode::ToZero,
    RoundingMode::FromZero,
];

#[test]
fn fuzz_mpfr_bit_exact_all_round_modes() {
    reset_test_rng();
    let mut cc = Consts::new().unwrap();

    unsafe {
        mpfr::set_emin(EXPONENT_MIN as gmp_mpfr_sys::gmp::exp_t);
        mpfr::set_emax(EXPONENT_MAX as gmp_mpfr_sys::gmp::exp_t);
    }
    assert_eq!(EXPONENT_MIN, exp_min());
    assert_eq!(EXPONENT_MAX, exp_max());

    let emin: Exponent = -64;
    let emax: Exponent = 64;
    let ops = [BinOp::Add, BinOp::Sub, BinOp::Mul, BinOp::Div];

    for p in [64usize, 128, 192, 320, 1024] {
        for i in 0..40 {
            // Half of the operands have long runs of ones then zeros, so sums and
            // differences land on exact ties and on carries out of the last word.
            let (n1, f1, n2, f2) = if i % 2 == 0 {
                let (n1, f1) = get_float_pair(p, emin, emax, &mut cc);
                let (n2, f2) = get_float_pair(p, emin, emax, &mut cc);
                (n1, f1, n2, f2)
            } else {
                let mut n1 = get_oned_zeroed(p, emin, emax);
                let mut n2 = get_oned_zeroed(p, emin, emax);
                if test_random::<bool>() {
                    n1 = n1.neg();
                }
                if test_random::<bool>() {
                    n2 = n2.neg();
                }
                let f1 = conv_to_mpfr(p, &n1, &mut cc);
                let f2 = conv_to_mpfr(p, &n2, &mut cc);
                (n1, f1, n2, f2)
            };
            for op in ops {
                for rm in ALL_MODES {
                    let n3 = match op {
                        BinOp::Add => n1.clone().add(&n2, p, rm),
                        BinOp::Sub => n1.clone().sub(&n2, p, rm),
                        BinOp::Mul => ExactNum::mul_ref(&n1, &n2, p, rm),
                        BinOp::Div => ExactNum::div_ref(&n1, &n2, p, rm),
                    };
                    let f3 = mpfr_reference(&f1, &f2, op, p as u32, rm);
                    let info = format!("#{i} p={p} {op:?} {rm:?} {n1:?} {n2:?}");
                    assert_exact_eq(&n3, &f3, p, &info, &mut cc);
                }
            }
            for (rm, rnd) in ROUNDS {
                let a = if n1.is_negative() { neg_of(&n1) } else { n1.clone() };
                let fa = Float::with_val(p as u32, f1.abs_ref());
                test_zf_op_no_cc!(
                    true,
                    a,
                    sqrt,
                    fa,
                    sqrt,
                    p,
                    rm,
                    rnd,
                    (i, &a, p, rm, "fuzz sqrt"),
                    cc
                );
            }
            if !n1.is_zero() && i < 8 {
                let a = if n1.is_negative() { n1.neg() } else { n1.clone() };
                let fa = Float::with_val(p as u32, f1.abs_ref());
                for n in [4usize, 5] {
                    for (rm, rnd) in ROUNDS {
                        let n3 = ExactNum::nth_root(&a, n, p, rm);
                        let mut f3 = Float::with_val(p as u32, 1);
                        unsafe {
                            mpfr::rootn_ui(f3.as_raw_mut(), fa.as_raw(), n as u64, rnd);
                        }
                        assert_float_close(
                            n3,
                            f3,
                            p,
                            &format!("{:?}", (i, n, &a, p, rm, "fuzz nth_root")),
                            false,
                            &mut cc,
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn fuzz_mpfr_overflow_all_round_modes() {
    reset_test_rng();
    let mut cc = Consts::new().unwrap();
    unsafe {
        mpfr::set_emin(EXPONENT_MIN as gmp_mpfr_sys::gmp::exp_t);
        mpfr::set_emax(EXPONENT_MAX as gmp_mpfr_sys::gmp::exp_t);
    }

    for &p in &[64, 128, 192, 256, 512] {
        let max = ExactNum::max_value(p);
        let seconds = vec![
            max.clone(),
            neg_of(&max),
            ExactNum::from_i8(2, p),
            ExactNum::from_i8(-2, p),
            {
                let mut x = ExactNum::from_i8(3, p);
                x.set_exponent(1);
                x
            },
            {
                let mut x = ExactNum::from_i8(1, p);
                x.set_exponent(0);
                x
            },
            {
                let mut x = ExactNum::from_i8(-1, p);
                x.set_exponent(-1);
                x
            },
        ];

        let neg_max = neg_of(&max);
        for n1 in [&max, &neg_max] {
            for n2 in &seconds {
                let f1 = conv_to_mpfr(p, n1, &mut cc);
                let f2 = conv_to_mpfr(p, n2, &mut cc);
                let a = n1.clone();
                let b = n2.clone();

                for &(rm, rnd) in ROUNDS.iter() {
                    unsafe { mpfr::clear_flags() };
                    test_zf_op_no_cc!(
                        true,
                        a,
                        b,
                        add,
                        f1,
                        f2,
                        add,
                        p,
                        rm,
                        rnd,
                        (p, &a, &b, rm, "overflow add"),
                        cc
                    );

                    unsafe { mpfr::clear_flags() };
                    test_zf_op_no_cc!(
                        true,
                        a,
                        b,
                        sub,
                        f1,
                        f2,
                        sub,
                        p,
                        rm,
                        rnd,
                        (p, &a, &b, rm, "overflow sub"),
                        cc
                    );

                    unsafe { mpfr::clear_flags() };
                    test_zf_op_no_cc!(
                        true,
                        a,
                        b,
                        mul,
                        f1,
                        f2,
                        mul,
                        p,
                        rm,
                        rnd,
                        (p, &a, &b, rm, "overflow mul"),
                        cc
                    );

                    unsafe { mpfr::clear_flags() };
                    test_zf_op_no_cc!(
                        true,
                        a,
                        b,
                        div,
                        f1,
                        f2,
                        div,
                        p,
                        rm,
                        rnd,
                        (p, &a, &b, rm, "overflow div"),
                        cc
                    );
                }
            }
        }
    }
}

/// Results are correctly rounded at exactly `p` bits, not at the next multiple of the word
/// size: every mode (ToOdd included), both signs, four operations, word-aligned and
/// unaligned precisions, and exponents in the middle, near `EXPONENT_MAX` (overflow), and
/// near `EXPONENT_MIN` (subnormal results).
#[test]
fn fuzz_mpfr_exact_bit_precision() {
    reset_test_rng();
    let mut cc = Consts::new().unwrap();
    unsafe {
        mpfr::set_emin(EXPONENT_MIN as gmp_mpfr_sys::gmp::exp_t);
        mpfr::set_emax(EXPONENT_MAX as gmp_mpfr_sys::gmp::exp_t);
    }
    let modes = [
        RoundingMode::ToEven,
        RoundingMode::ToOdd,
        RoundingMode::Up,
        RoundingMode::Down,
        RoundingMode::ToZero,
        RoundingMode::FromZero,
    ];
    let ops = [BinOp::Add, BinOp::Sub, BinOp::Mul, BinOp::Div];
    // (band, first exponent range, second exponent range)
    let bands: [(&str, Exponent, Exponent, Exponent, Exponent); 3] = [
        ("mid", -40, 40, -40, 40),
        ("high", EXPONENT_MAX - 3, EXPONENT_MAX, -2, 3),
        ("low", EXPONENT_MIN, EXPONENT_MIN + 3, -3, 2),
    ];
    for p in [24usize, 53, 64, 65, 100, 113, 127, 128, 130, 200, 255] {
        let storage = p.div_ceil(WORD_BIT_SIZE) * WORD_BIT_SIZE;
        for (band, e1lo, e1hi, e2lo, e2hi) in bands {
            for i in 0..60 {
                let mut n1 = ExactNum::random_normal(storage, e1lo, e1hi);
                let mut n2 = ExactNum::random_normal(storage, e2lo, e2hi);
                n1.set_precision(p, RoundingMode::ToZero).unwrap();
                n2.set_precision(p, RoundingMode::ToZero).unwrap();
                if test_random::<bool>() {
                    n1 = n1.neg();
                }
                if test_random::<bool>() {
                    n2 = n2.neg();
                }
                let f1 = conv_to_mpfr(p, &n1, &mut cc);
                let f2 = conv_to_mpfr(p, &n2, &mut cc);
                for op in ops {
                    for rm in modes {
                        // the ToOdd emulation is for normal-range results only
                        if rm == RoundingMode::ToOdd && band == "low" {
                            continue;
                        }
                        let n3 = zenith_binop(&n1, &n2, op, p, rm);
                        let f3 = mpfr_reference(&f1, &f2, op, p as u32, rm);
                        let info = format!("{band} #{i} p={p} {op:?} {rm:?} {n1:?} {n2:?}");
                        assert_exact_eq(&n3, &f3, p, &info, &mut cc);
                    }
                }
            }
        }
    }
}

/// Regression for the f7de24b rewrite of `fuzz_mpfr_exact_bit_precision`: the mul/div arms built a fresh random
/// first operand instead of using `n1`, so zenith computed a different product than MPFR (`mid #0 p=24 Mul ToEven`).
/// The operands are the exact ones of that failure; the oracle is MPFR (`mpfr::mul`, RNDN, p = 24).
#[test]
fn zenith_binop_uses_both_operands_mpfr_mul_p24() {
    let mut cc = Consts::new().unwrap();
    unsafe {
        mpfr::set_emin(EXPONENT_MIN as gmp_mpfr_sys::gmp::exp_t);
        mpfr::set_emax(EXPONENT_MAX as gmp_mpfr_sys::gmp::exp_t);
    }
    let p = 24usize;
    let n1 = ExactNum::from_raw_parts(
        &[9994389172669906944],
        64,
        zenith_float_num::Sign::Pos,
        -5,
        true,
    );
    let n2 = ExactNum::from_raw_parts(
        &[16400126323906838528],
        64,
        zenith_float_num::Sign::Neg,
        38,
        true,
    );
    let f1 = conv_to_mpfr(p, &n1, &mut cc);
    let f2 = conv_to_mpfr(p, &n2, &mut cc);
    let f3 = mpfr_reference(&f1, &f2, BinOp::Mul, p as u32, RoundingMode::ToEven);
    // the value MPFR reports for this product (an outside oracle, independent of zenith)
    let want = -4.13765094e9_f64;
    assert!(
        ((f3.to_f64() - want) / want).abs() < 1e-8,
        "MPFR product {f3}"
    );
    for rm in [RoundingMode::ToEven, RoundingMode::ToZero, RoundingMode::Up] {
        let f3 = mpfr_reference(&f1, &f2, BinOp::Mul, p as u32, rm);
        let n3 = zenith_binop(&n1, &n2, BinOp::Mul, p, rm);
        assert_exact_eq(&n3, &f3, p, &format!("mul {rm:?}"), &mut cc);
        let f4 = mpfr_reference(&f1, &f2, BinOp::Div, p as u32, rm);
        let n4 = zenith_binop(&n1, &n2, BinOp::Div, p, rm);
        assert_exact_eq(&n4, &f4, p, &format!("div {rm:?}"), &mut cc);
    }
}
