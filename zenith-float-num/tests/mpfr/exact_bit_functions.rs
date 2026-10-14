//! Functions beyond add / sub / mul / div round at exactly `p` bits in every mode, and their
//! `inexact()` reports that one operation (2.0.0 audit, round 2: items a and c).

use crate::mpfr::common::assert_exact_eq;
use crate::mpfr::common::conv_to_mpfr;
use crate::mpfr::common::mpfr_reference_with;
use crate::mpfr::common::reset_test_rng;
use crate::mpfr::common::test_random;
use gmp_mpfr_sys::mpfr::{self, rnd_t};
use rug::Float;
use zenith_float_num::{Consts, ExactNum, Exponent, RoundingMode, WORD_BIT_SIZE};

pub const PRECISIONS: [usize; 11] = [24, 53, 64, 65, 100, 113, 127, 128, 130, 200, 255];

pub const MODES: [RoundingMode; 6] = [
    RoundingMode::ToEven,
    RoundingMode::ToOdd,
    RoundingMode::Up,
    RoundingMode::Down,
    RoundingMode::ToZero,
    RoundingMode::FromZero,
];

/// A `p`-bit input in `2^elo ..= 2^ehi`, stored at the word precision zenith uses for `p`.
/// Its accumulated `inexact` flag is set so that a result cannot pass by inheriting it.
pub fn input(p: usize, elo: Exponent, ehi: Exponent, signed: bool) -> ExactNum {
    let storage = p.div_ceil(WORD_BIT_SIZE) * WORD_BIT_SIZE;
    let mut n = ExactNum::random_normal(storage, elo, ehi);
    n.set_precision(p, RoundingMode::ToZero).unwrap();
    // random_normal picks the sign at random
    if n.is_negative() != (signed && test_random::<bool>()) {
        n = n.neg();
    }
    n.set_inexact(true);
    n
}

/// `p`-bit value of a small integer ratio `num / 2^k`, stored at the word precision.
pub fn exact_input(num: i64, k: i32, p: usize) -> ExactNum {
    let storage = p.div_ceil(WORD_BIT_SIZE) * WORD_BIT_SIZE;
    let mut n = ExactNum::from_i64(num, storage);
    n = n.ldexp(-k, storage, RoundingMode::ToEven);
    n.set_inexact(true);
    n
}

/// Checks one result against MPFR: value bit for bit and the per-operation flag.
pub fn check(
    got: &ExactNum,
    p: usize,
    rm: RoundingMode,
    info: &str,
    cc: &mut Consts,
    f: impl Fn(&mut Float, rnd_t) -> i32,
) {
    let (want, inexact) = mpfr_reference_with(p as u32, rm, f);
    assert_exact_eq(got, &want, p, info, cc);
    if !want.is_infinite() && !want.is_nan() {
        assert_eq!(got.inexact(), inexact, "{info}: inexact flag, MPFR {want}");
    }
}

/// Sweeps a unary function over the precisions and modes on random inputs in each exponent
/// band and on the given exact inputs (`num / 2^k`).
pub fn sweep_unary(
    name: &str,
    bands: &[(Exponent, Exponent, bool)],
    exact: &[(i64, i32)],
    per_band: usize,
    zf: impl Fn(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
    mp: impl Fn(&mut Float, &Float, rnd_t) -> i32,
) {
    reset_test_rng();
    let mut cc = Consts::new().unwrap();
    for p in PRECISIONS {
        let mut inputs = Vec::new();
        for &(elo, ehi, signed) in bands {
            for _ in 0..per_band {
                inputs.push(input(p, elo, ehi, signed));
            }
        }
        for &(num, k) in exact {
            inputs.push(exact_input(num, k, p));
        }
        for x in &inputs {
            let fx = conv_to_mpfr(p.div_ceil(WORD_BIT_SIZE) * WORD_BIT_SIZE, x, &mut cc);
            for rm in MODES {
                let got = zf(x, p, rm, &mut cc);
                let info = format!("{name} p={p} {rm:?} x={x:?}");
                check(&got, p, rm, &info, &mut cc, |r, rnd| mp(r, &fx, rnd));
            }
        }
    }
}

#[test]
fn exact_bit_sqrt() {
    sweep_unary(
        "sqrt",
        &[(-40, 40, false), (-3, 3, false)],
        &[(1, 0), (4, 0), (9, 0), (1, 2), (49, 6), (2, 0), (3, 0)],
        12,
        |x, p, rm, _| x.sqrt(p, rm),
        |r, a, rnd| unsafe { mpfr::sqrt(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_cbrt() {
    sweep_unary(
        "cbrt",
        &[(-40, 40, true)],
        &[(27, 0), (-8, 0), (1, 3), (2, 0)],
        6,
        |x, p, rm, _| x.cbrt(p, rm),
        |r, a, rnd| unsafe { mpfr::cbrt(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_exp() {
    sweep_unary(
        "exp",
        &[(-6, 4, true)],
        &[(0, 0), (1, 0), (-1, 1)],
        6,
        |x, p, rm, cc| x.exp(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::exp(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_exp2() {
    sweep_unary(
        "exp2",
        &[(-6, 4, true)],
        &[(0, 0), (3, 0), (-5, 0), (1, 1)],
        6,
        |x, p, rm, cc| x.exp2(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::exp2(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_exp10() {
    sweep_unary(
        "exp10",
        &[(-6, 3, true)],
        &[(0, 0), (2, 0), (-1, 0)],
        6,
        |x, p, rm, cc| x.exp10(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::exp10(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_ln() {
    sweep_unary(
        "ln",
        &[(-40, 40, false), (0, 1, false)],
        &[(1, 0), (2, 0), (1, 1)],
        6,
        |x, p, rm, cc| x.ln(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::log(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_log2() {
    sweep_unary(
        "log2",
        &[(-40, 40, false)],
        &[(1, 0), (8, 0), (1, 5), (3, 0)],
        6,
        |x, p, rm, cc| x.log2(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::log2(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_log10() {
    sweep_unary(
        "log10",
        &[(-40, 40, false)],
        &[(1, 0), (100, 0), (10, 0), (2, 0)],
        6,
        |x, p, rm, cc| x.log10(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::log10(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_log1p() {
    sweep_unary(
        "log1p",
        &[(-8, 0, true), (0, 10, false)],
        &[(0, 0), (1, 0), (-1, 1)],
        6,
        |x, p, rm, cc| x.log1p(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::log1p(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_expm1() {
    sweep_unary(
        "expm1",
        &[(-8, 3, true)],
        &[(0, 0), (1, 0)],
        6,
        |x, p, rm, cc| x.expm1(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::expm1(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_sin() {
    sweep_unary(
        "sin",
        &[(-6, 4, true)],
        &[(0, 0), (1, 0), (3, 0)],
        6,
        |x, p, rm, cc| x.sin(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::sin(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_cos() {
    sweep_unary(
        "cos",
        &[(-6, 4, true)],
        &[(0, 0), (1, 0), (3, 0)],
        6,
        |x, p, rm, cc| x.cos(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::cos(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_tan() {
    sweep_unary(
        "tan",
        &[(-6, 3, true)],
        &[(0, 0), (1, 0)],
        6,
        |x, p, rm, cc| x.tan(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::tan(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_asin() {
    sweep_unary(
        "asin",
        &[(-8, 0, true)],
        &[(0, 0), (1, 0), (-1, 0), (1, 1)],
        6,
        |x, p, rm, cc| x.asin(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::asin(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_acos() {
    sweep_unary(
        "acos",
        &[(-8, 0, true)],
        &[(0, 0), (1, 0), (-1, 0), (1, 1)],
        6,
        |x, p, rm, cc| x.acos(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::acos(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_atan() {
    sweep_unary(
        "atan",
        &[(-6, 6, true)],
        &[(0, 0), (1, 0), (-1, 0)],
        6,
        |x, p, rm, cc| x.atan(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::atan(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_sinh() {
    sweep_unary(
        "sinh",
        &[(-6, 4, true)],
        &[(0, 0), (1, 0)],
        6,
        |x, p, rm, cc| x.sinh(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::sinh(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_cosh() {
    sweep_unary(
        "cosh",
        &[(-6, 4, true)],
        &[(0, 0), (1, 0)],
        6,
        |x, p, rm, cc| x.cosh(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::cosh(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_tanh() {
    sweep_unary(
        "tanh",
        &[(-6, 4, true)],
        &[(0, 0), (1, 0)],
        6,
        |x, p, rm, cc| x.tanh(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::tanh(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_asinh() {
    sweep_unary(
        "asinh",
        &[(-6, 6, true)],
        &[(0, 0), (1, 0)],
        6,
        |x, p, rm, cc| x.asinh(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::asinh(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_acosh() {
    sweep_unary(
        "acosh",
        &[(1, 20, false)],
        &[(1, 0), (2, 0)],
        6,
        |x, p, rm, cc| x.acosh(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::acosh(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_atanh() {
    sweep_unary(
        "atanh",
        &[(-8, 0, true)],
        &[(0, 0), (1, 1), (-1, 2)],
        6,
        |x, p, rm, cc| x.atanh(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::atanh(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_gamma() {
    sweep_unary(
        "gamma",
        &[(-3, 4, false)],
        &[(1, 0), (3, 0), (5, 0), (1, 1)],
        6,
        |x, p, rm, cc| x.gamma(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::gamma(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_ln_gamma() {
    sweep_unary(
        "ln_gamma",
        &[(-3, 6, false)],
        &[(1, 0), (2, 0), (3, 0)],
        6,
        |x, p, rm, cc| x.ln_gamma(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::lngamma(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_digamma() {
    sweep_unary(
        "digamma",
        &[(-3, 6, false)],
        &[(1, 0), (2, 0)],
        6,
        |x, p, rm, cc| x.digamma(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::digamma(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_erf() {
    sweep_unary(
        "erf",
        &[(-6, 2, true)],
        &[(0, 0), (1, 0)],
        6,
        |x, p, rm, cc| x.erf(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::erf(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_erfc() {
    sweep_unary(
        "erfc",
        &[(-6, 2, true)],
        &[(0, 0), (1, 0)],
        6,
        |x, p, rm, cc| x.erfc(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::erfc(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

#[test]
fn exact_bit_ei() {
    sweep_unary(
        "ei",
        &[(-6, 4, false)],
        &[(1, 0), (2, 0)],
        6,
        |x, p, rm, cc| x.ei(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::eint(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}

// Pair functions: both results round at exactly `p` bits with their own flag.

#[test]
fn exact_bit_sin_cos_pair() {
    for (part, mp) in [(0usize, mpfr::sin as unsafe extern "C" fn(_, _, _) -> _), (1, mpfr::cos)] {
        sweep_unary(
            if part == 0 { "sin_cos.0" } else { "sin_cos.1" },
            &[(-6, 4, true)],
            &[(0, 0), (1, 0), (3, 0)],
            6,
            |x, p, rm, cc| {
                let (s, c) = x.sin_cos(p, rm, cc);
                if part == 0 {
                    s
                } else {
                    c
                }
            },
            |r, a, rnd| unsafe { mp(r.as_raw_mut(), a.as_raw(), rnd) },
        );
    }
}

#[test]
fn exact_bit_sinh_cosh_pair() {
    for (part, mp) in [(0usize, mpfr::sinh as unsafe extern "C" fn(_, _, _) -> _), (1, mpfr::cosh)]
    {
        sweep_unary(
            if part == 0 { "sinh_cosh.0" } else { "sinh_cosh.1" },
            &[(-6, 4, true), (32, 33, true)],
            &[(0, 0), (1, 0)],
            6,
            |x, p, rm, cc| {
                let (s, c) = x.sinh_cosh(p, rm, cc);
                if part == 0 {
                    s
                } else {
                    c
                }
            },
            |r, a, rnd| unsafe { mp(r.as_raw_mut(), a.as_raw(), rnd) },
        );
    }
}
