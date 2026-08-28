//! MPFR oracles for rectangular `ExactComplex` (re/im via real MPFR ops).

use crate::mpfr::common::{
    assert_float_close, get_float_pair, get_random_rnd_pair, reset_test_rng,
};
use gmp_mpfr_sys::{gmp::exp_t, mpfr};
use rug::{
    float::{exp_max, exp_min},
    Float,
};
use zenith_float_num::{
    Consts, ExactComplex, EXPONENT_MAX, EXPONENT_MIN, WORD_BIT_SIZE,
};

#[test]
fn mpfr_compare_complex() {
    reset_test_rng();
    let mut cc = Consts::new().unwrap();

    unsafe {
        mpfr::set_emin(EXPONENT_MIN as exp_t);
        mpfr::set_emax(EXPONENT_MAX as exp_t);
    }
    assert_eq!(EXPONENT_MIN, exp_min());
    assert_eq!(EXPONENT_MAX, exp_max());

    let p = 2 * WORD_BIT_SIZE;

    for _ in 0..12 {
        let (rm, rnd) = get_random_rnd_pair();
        let (ar, far) = get_float_pair(p, -8, 8, &mut cc);
        let (ai, fai) = get_float_pair(p, -8, 8, &mut cc);
        let (br, fbr) = get_float_pair(p, -8, 8, &mut cc);
        let (bi, fbi) = get_float_pair(p, -8, 8, &mut cc);
        if ar.is_nan() || ai.is_nan() || br.is_nan() || bi.is_nan() {
            continue;
        }

        let a = ExactComplex::new(ar, ai);
        let b = ExactComplex::new(br, bi);

        let s = a.add(&b, p, rm);
        let mut fre = Float::with_val(p as u32, 1);
        let mut fim = Float::with_val(p as u32, 1);
        unsafe {
            mpfr::add(fre.as_raw_mut(), far.as_raw(), fbr.as_raw(), rnd);
            mpfr::add(fim.as_raw_mut(), fai.as_raw(), fbi.as_raw(), rnd);
        }
        assert_float_close(s.re().clone(), fre, p, "cplx add re", true, &mut cc);
        assert_float_close(s.im().clone(), fim, p, "cplx add im", true, &mut cc);
    }
}
