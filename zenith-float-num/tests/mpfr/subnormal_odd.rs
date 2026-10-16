//! Results in the subnormal band round correctly in every mode, including `ToOdd`
//! (round half to odd) at exact ties on the subnormal grid (2.0.0 audit, round 2: item d).

use crate::mpfr::common::{assert_exact_eq, conv_to_mpfr, mpfr_reference_with, reset_test_rng};
use crate::mpfr::exact_bit_functions::{input, MODES, PRECISIONS};
use gmp_mpfr_sys::mpfr::{self, rnd_t};
use rug::Float;
use zenith_float_num::{Consts, ExactNum, RoundingMode, EXPONENT_MIN, WORD_BIT_SIZE};

/// `conv_to_mpfr` for values near the bottom of the range: converts `x / 2^e`, then scales.
fn conv_small(w: usize, x: &ExactNum, cc: &mut Consts) -> Float {
    let e = x.exponent().unwrap_or(0);
    let s = x.clone().ldexp(-e, w, RoundingMode::None);
    let f = conv_to_mpfr(w, &s, cc);
    let old = unsafe { mpfr::get_emin() };
    unsafe { mpfr::set_emin(mpfr::get_emin_min()) };
    let f = f << e;
    unsafe { mpfr::set_emin(old) };
    f
}

type MpOp<'a> = Box<dyn Fn(&mut Float, rnd_t) -> i32 + 'a>;

/// Reference for an operation whose exact value `f` computes at `wide` bits: `ToOdd`
/// is `ToEven` except at an exact tie, where it takes the other neighbour.
fn reference(
    p: usize,
    rm: RoundingMode,
    wide: u32,
    f: &dyn Fn(&mut Float, rnd_t) -> i32,
) -> (Float, bool) {
    if rm != RoundingMode::ToOdd {
        return mpfr_reference_with(p as u32, rm, f);
    }
    let (rn, inexact) = mpfr_reference_with(p as u32, RoundingMode::ToEven, f);
    let (lo, _) = mpfr_reference_with(p as u32, RoundingMode::ToZero, f);
    let (hi, _) = mpfr_reference_with(p as u32, RoundingMode::FromZero, f);
    let old = unsafe { (mpfr::get_emin(), mpfr::get_emax()) };
    unsafe {
        mpfr::set_emin(mpfr::get_emin_min());
        mpfr::set_emax(mpfr::get_emax_max());
    }
    let mut h = Float::new(wide);
    let t = f(&mut h, rnd_t::RNDN);
    assert_eq!(t, 0, "wide reference must be exact");
    let mut sum = Float::new(wide + 2);
    let mut mid = Float::new(wide + 2);
    unsafe {
        mpfr::add(sum.as_raw_mut(), lo.as_raw(), hi.as_raw(), rnd_t::RNDN);
        mpfr::div_2ui(mid.as_raw_mut(), sum.as_raw(), 1, rnd_t::RNDN);
        mpfr::set_emin(old.0);
        mpfr::set_emax(old.1);
    }
    if lo != hi && mid == h {
        (if rn == lo { hi } else { lo }, inexact)
    } else {
        (rn, inexact)
    }
}

#[test]
fn subnormal_band_every_mode() {
    reset_test_rng();
    let mut cc = Consts::new().unwrap();
    let lo = EXPONENT_MIN;
    for p in PRECISIONS {
        let w = p.div_ceil(WORD_BIT_SIZE) * WORD_BIT_SIZE;
        let wide = (4 * w + 256) as u32;
        for i in 0..24 {
            let x = input(p, lo + 1, lo + 2 * p as i32, true);
            let y = input(p, -(p as i32) - 2, -1, true);
            let z = input(p, lo + 1, lo + p as i32, true);
            // shifts that land exactly on and between subnormal grid points
            let k = i % 7 + 1;
            let fx = conv_small(w, &x, &mut cc);
            let fy = conv_small(w, &y, &mut cc);
            let fz = conv_small(w, &z, &mut cc);
            for rm in MODES {
                let cases: [(&str, ExactNum, MpOp); 5] = [
                    (
                        "mul",
                        ExactNum::mul_ref(&x, &y, p, rm),
                        Box::new(|r: &mut Float, rnd| unsafe {
                            mpfr::mul(r.as_raw_mut(), fx.as_raw(), fy.as_raw(), rnd)
                        }),
                    ),
                    (
                        "add",
                        ExactNum::add_ref(&x, &z, p, rm),
                        Box::new(|r: &mut Float, rnd| unsafe {
                            mpfr::add(r.as_raw_mut(), fx.as_raw(), fz.as_raw(), rnd)
                        }),
                    ),
                    (
                        "sub",
                        ExactNum::sub_ref(&x, &z, p, rm),
                        Box::new(|r: &mut Float, rnd| unsafe {
                            mpfr::sub(r.as_raw_mut(), fx.as_raw(), fz.as_raw(), rnd)
                        }),
                    ),
                    (
                        "fma",
                        x.fma(&y, &z, p, rm),
                        Box::new(|r: &mut Float, rnd| unsafe {
                            mpfr::fma(r.as_raw_mut(), fx.as_raw(), fy.as_raw(), fz.as_raw(), rnd)
                        }),
                    ),
                    (
                        "ldexp",
                        x.clone().ldexp(-k, p, rm),
                        Box::new(|r: &mut Float, rnd| unsafe {
                            mpfr::div_2ui(r.as_raw_mut(), fx.as_raw(), k as _, rnd)
                        }),
                    ),
                ];
                for (name, got, f) in cases.iter() {
                    let info = format!("{name} p={p} {rm:?} x={x:?}");
                    let (want, inexact) = reference(p, rm, wide, f.as_ref());
                    assert_exact_eq(got, &want, p, &info, &mut cc);
                    if !want.is_zero() {
                        assert_eq!(got.inexact(), inexact, "{info}: inexact flag, MPFR {want}");
                    }
                }
            }
        }
    }
}
