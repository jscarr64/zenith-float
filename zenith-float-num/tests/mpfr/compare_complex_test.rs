//! MPFR oracles for rectangular `ExactComplex` (re/im via real MPFR ops).

use crate::mpfr::common::{
    assert_float_close, get_float_pair, get_random_rnd_pair, reset_test_rng,
};
use gmp_mpfr_sys::{gmp::exp_t, mpfr};
use rug::{
    float::{exp_max, exp_min},
    Float,
};
use zenith_float_num::{Consts, ExactComplex, ExactNum, EXPONENT_MAX, EXPONENT_MIN, WORD_BIT_SIZE};

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
        assert_float_close(s.re().clone(), fre.clone(), p, "cplx add re", true, &mut cc);
        assert_float_close(s.im().clone(), fim.clone(), p, "cplx add im", true, &mut cc);

        let pr = a.mul(&b, p, rm);
        let mut fac = Float::with_val(p as u32, 1);
        let mut fbd = Float::with_val(p as u32, 1);
        let mut fad = Float::with_val(p as u32, 1);
        let mut fbc = Float::with_val(p as u32, 1);
        unsafe {
            mpfr::mul(fac.as_raw_mut(), far.as_raw(), fbr.as_raw(), rnd);
            mpfr::mul(fbd.as_raw_mut(), fai.as_raw(), fbi.as_raw(), rnd);
            mpfr::mul(fad.as_raw_mut(), far.as_raw(), fbi.as_raw(), rnd);
            mpfr::mul(fbc.as_raw_mut(), fai.as_raw(), fbr.as_raw(), rnd);
            mpfr::sub(fre.as_raw_mut(), fac.as_raw(), fbd.as_raw(), rnd);
            mpfr::add(fim.as_raw_mut(), fad.as_raw(), fbc.as_raw(), rnd);
        }
        assert_float_close(pr.re().clone(), fre, p, "cplx mul re", true, &mut cc);
        assert_float_close(pr.im().clone(), fim, p, "cplx mul im", true, &mut cc);
    }
}

/// Complex specials on the real axis vs MPFR. GNU MPC has no erf / Γ / Ai / J_n.
#[test]
fn mpfr_compare_complex_real_axis() {
    reset_test_rng();
    let mut cc = Consts::new().unwrap();

    unsafe {
        mpfr::set_emin(EXPONENT_MIN as exp_t);
        mpfr::set_emax(EXPONENT_MAX as exp_t);
    }

    let p = 2 * WORD_BIT_SIZE;

    for _ in 0..8 {
        let (n1, f1) = get_float_pair(p, -1, 2, &mut cc);
        if n1.is_nan() || n1.is_inf() {
            continue;
        }
        let (rm, rnd) = get_random_rnd_pair();
        let z = ExactComplex::from_real(n1.clone(), p);

        let er = z.erf(p, rm, &mut cc);
        let mut fe = Float::with_val(p as u32, 1);
        unsafe {
            mpfr::erf(fe.as_raw_mut(), f1.as_raw(), rnd);
        }
        assert_float_close(er.re().clone(), fe, p, "cplx erf|R re", false, &mut cc);
        assert!(
            er.im().is_zero() || er.im().exponent().is_some_and(|e| e < -((p as i32) / 4)),
            "cplx erf|R im"
        );

        if n1.is_positive() {
            let g = z.gamma(p, rm, &mut cc);
            let mut fg = Float::with_val(p as u32, 1);
            unsafe {
                mpfr::gamma(fg.as_raw_mut(), f1.as_raw(), rnd);
            }
            assert_float_close(g.re().clone(), fg, p, "cplx gamma|R re", false, &mut cc);
        }

        let a = z.ai(p, rm, &mut cc);
        let mut fa = Float::with_val(p as u32, 1);
        unsafe {
            mpfr::ai(fa.as_raw_mut(), f1.as_raw(), rnd);
        }
        assert_float_close(a.re().clone(), fa, p, "cplx ai|R re", false, &mut cc);
    }

    for n_ord in [0i64, 1] {
        for _ in 0..4 {
            let (n1, f1) = get_float_pair(p, -1, 2, &mut cc);
            if n1.is_nan() || n1.is_inf() {
                continue;
            }
            let (rm, rnd) = get_random_rnd_pair();
            let z = ExactComplex::from_real(n1.clone(), p);
            let nu = ExactComplex::from_real(ExactNum::from_i64(n_ord, p), p);
            let j = z.bessel_j_nu(&nu, p, rm, &mut cc);
            let mut fj = Float::with_val(p as u32, 1);
            unsafe {
                mpfr::jn(fj.as_raw_mut(), n_ord, f1.as_raw(), rnd);
            }
            assert_float_close(j.re().clone(), fj, p, &format!("cplx J_{n_ord}|R"), false, &mut cc);
        }
    }
}
