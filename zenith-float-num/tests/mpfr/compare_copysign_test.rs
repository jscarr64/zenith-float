//! MPFR oracle for `copysign`, including extreme precision.

use crate::mpfr::common::assert_float_close;
use crate::mpfr::common::get_float_pair;
use crate::mpfr::common::get_random_rnd_pair;
use crate::mpfr::common::reset_test_rng;
use crate::mpfr::common::test_zf_op_no_cc;
use gmp_mpfr_sys::{gmp::exp_t, mpfr};
use rug::float::{exp_max, exp_min, Special};
use rug::Float;
use zenith_float_num::{
    Consts, ExactNum, Exponent, RoundingMode, Sign, EXPONENT_MAX, EXPONENT_MIN, INF_NEG, INF_POS,
    WORD_BIT_SIZE,
};

fn set_mpfr_exp_range() {
    unsafe {
        mpfr::set_emin(EXPONENT_MIN as exp_t);
        mpfr::set_emax(EXPONENT_MAX as exp_t);
    }
    assert_eq!(EXPONENT_MIN, exp_min());
    assert_eq!(EXPONENT_MAX, exp_max());
}

fn run_copysign_pairs(run_cnt: usize, p: usize, emin: Exponent, emax: Exponent) {
    reset_test_rng();
    let mut cc = Consts::new().unwrap();
    set_mpfr_exp_range();

    for i in 0..run_cnt {
        let (n1, f1) = get_float_pair(p, emin, emax, &mut cc);
        let (n2, f2) = get_float_pair(p, emin, emax, &mut cc);
        let (rm, rnd) = get_random_rnd_pair();
        test_zf_op_no_cc!(
            true,
            n1,
            n2,
            copysign,
            f1,
            f2,
            copysign,
            p,
            rm,
            rnd,
            (i, p, rm, "copysign"),
            cc
        );
    }
}

#[test]
fn mpfr_compare_copysign() {
    run_copysign_pairs(64, 4 * WORD_BIT_SIZE, -64, 64);
}

#[test]
fn mpfr_compare_copysign_extreme() {
    // Same order of magnitude as `mpfr_compare_ops_large` (~10k+ bits).
    let p = 256 * WORD_BIT_SIZE;
    run_copysign_pairs(8, p, -128, 128);

    reset_test_rng();
    let mut cc = Consts::new().unwrap();
    set_mpfr_exp_range();
    let p1 = p;
    let p_out = 160 * WORD_BIT_SIZE;
    let (n1, f1) = get_float_pair(p1, -32, 32, &mut cc);
    let (n2, f2) = get_float_pair(p1, -32, 32, &mut cc);
    for _ in 0..5 {
        let (rm, rnd) = get_random_rnd_pair();
        test_zf_op_no_cc!(
            true,
            n1,
            n2,
            copysign,
            f1,
            f2,
            copysign,
            p_out,
            rm,
            rnd,
            (p1, p_out, rm, "copysign round down"),
            cc
        );
    }
}

#[test]
fn mpfr_compare_copysign_specials() {
    let mut cc = Consts::new().unwrap();
    set_mpfr_exp_range();
    let p = 2 * WORD_BIT_SIZE;
    let rm = RoundingMode::ToEven;
    let rnd = gmp_mpfr_sys::mpfr::rnd_t::RNDN;

    let mut pz = ExactNum::new(p);
    pz.set_sign(Sign::Pos);
    let mut nz = ExactNum::new(p);
    nz.set_sign(Sign::Neg);
    let fz = Float::with_val(p as u32, Special::Zero);
    let fnz = Float::with_val(p as u32, Special::NegZero);
    let finf = Float::with_val(p as u32, Special::Infinity);
    let fninf = Float::with_val(p as u32, Special::NegInfinity);
    reset_test_rng();
    let (nmag, fmag) = get_float_pair(p, -8, 8, &mut cc);

    test_zf_op_no_cc!(true, nmag, pz, copysign, fmag, fz, copysign, p, rm, rnd, "sign +0", cc);
    test_zf_op_no_cc!(true, nmag, nz, copysign, fmag, fnz, copysign, p, rm, rnd, "sign -0", cc);
    test_zf_op_no_cc!(
        true,
        nmag,
        INF_POS,
        copysign,
        fmag,
        finf,
        copysign,
        p,
        rm,
        rnd,
        "sign +inf",
        cc
    );
    test_zf_op_no_cc!(
        true,
        nmag,
        INF_NEG,
        copysign,
        fmag,
        fninf,
        copysign,
        p,
        rm,
        rnd,
        "sign -inf",
        cc
    );
    test_zf_op_no_cc!(
        true,
        INF_POS,
        nz,
        copysign,
        finf,
        fnz,
        copysign,
        p,
        rm,
        rnd,
        "inf copysign -0",
        cc
    );
}
