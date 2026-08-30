//! MPFR oracles for erf / erfc / gamma / ln_gamma / bessel_j and extra constants.
//! Bounded domains so Stirling / series stay inside CI time.

use crate::mpfr::common::{
    assert_float_close, conv_to_mpfr, get_float_pair, get_random_rnd_pair, reset_test_rng,
    test_random,
};
use gmp_mpfr_sys::{gmp::exp_t, mpfr};
use rug::{
    float::{exp_max, exp_min},
    Float,
};
use zenith_float_num::{Consts, ExactNum, EXPONENT_MAX, EXPONENT_MIN, WORD_BIT_SIZE};

#[test]
fn mpfr_compare_special_fns() {
    reset_test_rng();
    let mut cc = Consts::new().unwrap();

    unsafe {
        mpfr::set_emin(EXPONENT_MIN as exp_t);
        mpfr::set_emax(EXPONENT_MAX as exp_t);
    }
    assert_eq!(EXPONENT_MIN, exp_min());
    assert_eq!(EXPONENT_MAX, exp_max());

    let p = 2 * WORD_BIT_SIZE;
    let (rm, rnd) = get_random_rnd_pair();

    // Integer Γ(n) = (n-1)! — factorial path, 1-ULP vs MPFR.
    for k in 2u8..=12 {
        let n = ExactNum::from_u8(k, p);
        let f = conv_to_mpfr(p, &n, &mut cc);
        let g = n.gamma(p, rm, &mut cc);
        let mut fg = Float::with_val(p as u32, 1);
        unsafe {
            mpfr::gamma(fg.as_raw_mut(), f.as_raw(), rnd);
        }
        assert_float_close(g, fg, p, &format!("gamma({k})"), false, &mut cc);

        let lg = n.ln_gamma(p, rm, &mut cc);
        let mut flg = Float::with_val(p as u32, 1);
        unsafe {
            mpfr::lngamma(flg.as_raw_mut(), f.as_raw(), rnd);
        }
        assert_float_close(lg, flg, p, &format!("ln_gamma({k})"), false, &mut cc);
    }

    // Γ(1/2) = √π
    let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, rm);
    let fhalf = conv_to_mpfr(p, &half, &mut cc);
    let g = half.gamma(p, rm, &mut cc);
    let mut fg = Float::with_val(p as u32, 1);
    unsafe {
        mpfr::gamma(fg.as_raw_mut(), fhalf.as_raw(), rnd);
    }
    assert_float_close(g, fg, p, "gamma(1/2)", false, &mut cc);

    // erf / erfc on a modest interval
    for _ in 0..12 {
        let (n1, f1) = get_float_pair(p, -2, 2, &mut cc);
        if n1.is_nan() || n1.is_inf() {
            continue;
        }
        let (rm, rnd) = get_random_rnd_pair();
        let er = n1.erf(p, rm, &mut cc);
        let mut fe = Float::with_val(p as u32, 1);
        unsafe {
            mpfr::erf(fe.as_raw_mut(), f1.as_raw(), rnd);
        }
        assert_float_close(er, fe, p, "erf", false, &mut cc);

        let erc = n1.erfc(p, rm, &mut cc);
        let mut fec = Float::with_val(p as u32, 1);
        unsafe {
            mpfr::erfc(fec.as_raw_mut(), f1.as_raw(), rnd);
        }
        assert_float_close(erc, fec, p, "erfc", false, &mut cc);
    }

    // J_0, J_1, J_2 on |x| < 4
    for n_ord in [0i64, 1, 2] {
        for _ in 0..8 {
            let (n1, f1) = get_float_pair(p, -2, 2, &mut cc);
            if n1.is_nan() || n1.is_inf() {
                continue;
            }
            let (rm, rnd) = get_random_rnd_pair();
            let j = n1.bessel_j(n_ord as usize, p, rm, &mut cc);
            let mut fj = Float::with_val(p as u32, 1);
            unsafe {
                mpfr::jn(fj.as_raw_mut(), n_ord, f1.as_raw(), rnd);
            }
            assert_float_close(j, fj, p, &format!("bessel_j n={n_ord}"), false, &mut cc);
        }
    }

    // High integer order at modest |x|
    for n_ord in [20i64, 40] {
        let (n1, f1) = get_float_pair(p, -1, 1, &mut cc);
        if n1.is_nan() || n1.is_inf() {
            continue;
        }
        let (rm, rnd) = get_random_rnd_pair();
        let j = n1.bessel_j(n_ord as usize, p, rm, &mut cc);
        let mut fj = Float::with_val(p as u32, 1);
        unsafe {
            mpfr::jn(fj.as_raw_mut(), n_ord, f1.as_raw(), rnd);
        }
        assert_float_close(j, fj, p, &format!("bessel_j n={n_ord}"), false, &mut cc);
    }

    // Ei (MPFR eint) on x > 0 and Cauchy PV for x < 0
    for emin_emax in [(0, 2), (-2, 1)] {
        for _ in 0..8 {
            let (n1, f1) = get_float_pair(p, emin_emax.0, emin_emax.1, &mut cc);
            if n1.is_nan() || n1.is_inf() || n1.is_zero() {
                continue;
            }
            let (rm, rnd) = get_random_rnd_pair();
            let e = n1.ei(p, rm, &mut cc);
            let mut fe = Float::with_val(p as u32, 1);
            unsafe {
                mpfr::eint(fe.as_raw_mut(), f1.as_raw(), rnd);
            }
            assert_float_close(e, fe, p, "ei", false, &mut cc);
        }
    }

    // Y_0, Y_1 (MPFR yn) on 0 < x < 4
    for n_ord in [0i64, 1] {
        for _ in 0..8 {
            let (n1, f1) = get_float_pair(p, -1, 2, &mut cc);
            if n1.is_nan() || n1.is_inf() || !n1.is_positive() {
                continue;
            }
            let nu = ExactNum::from(n_ord as i32);
            let (rm, rnd) = get_random_rnd_pair();
            let y = n1.bessel_y(&nu, p, rm, &mut cc);
            let mut fy = Float::with_val(p as u32, 1);
            unsafe {
                mpfr::yn(fy.as_raw_mut(), n_ord, f1.as_raw(), rnd);
            }
            assert_float_close(y, fy, p, &format!("bessel_y n={n_ord}"), false, &mut cc);
        }
    }

    // digamma on z > 0
    for _ in 0..8 {
        let (n1, f1) = get_float_pair(p, 0, 2, &mut cc);
        if n1.is_nan() || n1.is_inf() || !n1.is_positive() {
            continue;
        }
        let (rm, rnd) = get_random_rnd_pair();
        let d = n1.digamma(p, rm, &mut cc);
        let mut fd = Float::with_val(p as u32, 1);
        unsafe {
            mpfr::digamma(fd.as_raw_mut(), f1.as_raw(), rnd);
        }
        assert_float_close(d, fd, p, "digamma", false, &mut cc);
    }

    let _ = test_random::<u8>();
}

#[test]
fn mpfr_compare_extra_constants() {
    reset_test_rng();
    let mut cc = Consts::new().unwrap();

    unsafe {
        mpfr::set_emin(EXPONENT_MIN as exp_t);
        mpfr::set_emax(EXPONENT_MAX as exp_t);
    }

    for _ in 0..24 {
        let p = (test_random::<usize>() % 4 + 2) * WORD_BIT_SIZE;
        let (rm, rnd) = get_random_rnd_pair();

        let n = cc.sqrt2(p, rm);
        let two = Float::with_val(p as u32, 2);
        let mut f = Float::with_val(p as u32, 1);
        unsafe {
            mpfr::sqrt(f.as_raw_mut(), two.as_raw(), rnd);
        }
        assert_float_close(n, f, p, "const sqrt2", true, &mut cc);

        let n = cc.phi(p, rm);
        let five = Float::with_val(p as u32, 5);
        let mut s5 = Float::with_val(p as u32, 1);
        unsafe {
            mpfr::sqrt(s5.as_raw_mut(), five.as_raw(), rnd);
        }
        let mut f = Float::with_val(p as u32, 1);
        unsafe {
            mpfr::add_ui(f.as_raw_mut(), s5.as_raw(), 1, rnd);
            mpfr::div_ui(f.as_raw_mut(), f.as_raw(), 2, rnd);
        }
        assert_float_close(n, f, p, "const phi", true, &mut cc);

        let n = cc.euler_gamma(p, rm);
        let mut f = Float::with_val(p as u32, 1);
        unsafe {
            mpfr::const_euler(f.as_raw_mut(), rnd);
        }
        assert_float_close(n, f, p, "const euler_gamma", false, &mut cc);
    }
}
