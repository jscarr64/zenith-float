//! Short MPFR differential: add / mul / sqrt under every IEEE rounding mode, and
//! exponent overflow of add / sub / mul / div (max finite toward zero, Inf otherwise).

use crate::mpfr::common::assert_exact_eq;
use crate::mpfr::common::assert_float_close;
use crate::mpfr::common::conv_to_mpfr;
use crate::mpfr::common::get_float_pair;
use crate::mpfr::common::mpfr_reference;
use crate::mpfr::common::reset_test_rng;
use crate::mpfr::common::test_random;
use crate::mpfr::common::test_zf_op_no_cc;
use crate::mpfr::common::BinOp;
use gmp_mpfr_sys::mpfr::{self, rnd_t};
use rug::float::{exp_max, exp_min};
use rug::Float;
use zenith_float_num::{
    Consts, ExactNum, Exponent, RoundingMode, EXPONENT_MAX, EXPONENT_MIN, WORD_BIT_SIZE,
};

const ROUNDS: [(RoundingMode, rnd_t); 5] = [
    (RoundingMode::ToEven, rnd_t::RNDN),
    (RoundingMode::Up, rnd_t::RNDU),
    (RoundingMode::Down, rnd_t::RNDD),
    (RoundingMode::FromZero, rnd_t::RNDA),
    (RoundingMode::ToZero, rnd_t::RNDZ),
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

    let p = 2 * WORD_BIT_SIZE;
    let emin: Exponent = -64;
    let emax: Exponent = 64;

    for i in 0..24 {
        let (n1, f1) = get_float_pair(p, emin, emax, &mut cc);
        let (n2, f2) = get_float_pair(p, emin, emax, &mut cc);
        for (rm, rnd) in ROUNDS {
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
                (i, &n1, &n2, p, rm, "fuzz add"),
                cc
            );
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
                (i, &n1, &n2, p, rm, "fuzz mul"),
                cc
            );
            if !n1.is_negative() {
                test_zf_op_no_cc!(
                    true,
                    n1,
                    sqrt,
                    f1,
                    sqrt,
                    p,
                    rm,
                    rnd,
                    (i, &n1, p, rm, "fuzz sqrt"),
                    cc
                );
            }
        }
        if !n1.is_negative() && !n1.is_zero() && i < 8 {
            for n in [4usize, 5] {
                for (rm, rnd) in ROUNDS {
                    let n3 = ExactNum::nth_root(&n1, n, p, rm);
                    let mut f3 = Float::with_val(p as u32, 1);
                    unsafe {
                        mpfr::rootn_ui(f3.as_raw_mut(), f1.as_raw(), n as u64, rnd);
                    }
                    assert_float_close(
                        n3,
                        f3,
                        p,
                        &format!("{:?}", (i, n, &n1, p, rm, "fuzz nth_root")),
                        false,
                        &mut cc,
                    );
                }
            }
        }
        let _ = (f1, f2, i);
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
            max.neg(),
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

        for n1 in [&max, &max.neg()] {
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
                        let n3 = match op {
                            BinOp::Add => n1.add(&n2, p, rm),
                            BinOp::Sub => n1.sub(&n2, p, rm),
                            BinOp::Mul => n1.mul(&n2, p, rm),
                            BinOp::Div => n1.div(&n2, p, rm),
                        };
                        let f3 = mpfr_reference(&f1, &f2, op, p as u32, rm);
                        let info = format!("{band} #{i} p={p} {op:?} {rm:?} {n1:?} {n2:?}");
                        assert_exact_eq(&n3, &f3, p, &info, &mut cc);
                    }
                }
            }
        }
    }
}
