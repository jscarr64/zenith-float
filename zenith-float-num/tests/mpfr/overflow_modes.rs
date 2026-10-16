//! Overflow rounds by mode in every operation: `±max` under `ToZero`, under `Down` for a
//! positive result and under `Up` for a negative result; `±inf` otherwise
//! (2.0.0 audit, round 2: item b).

use crate::mpfr::common::conv_to_mpfr;
use crate::mpfr::common::mpfr_reference_with;
use crate::mpfr::common::reset_test_rng;
use crate::mpfr::exact_bit_functions::{check, input, MODES, PRECISIONS};
use gmp_mpfr_sys::mpfr::{self, rnd_t};
use rug::Float;
use zenith_float_num::{Consts, ExactNum, RoundingMode, EXPONENT_MAX, WORD_BIT_SIZE};

/// `-x` for an operand that stays borrowed: fresh `-1` times `x` at `x`'s own precision.
fn neg_of(x: &ExactNum) -> ExactNum {
    let p = x.mantissa_max_bit_len().unwrap_or(WORD_BIT_SIZE);
    ExactNum::from_i8(-1, p).mul(x, p, RoundingMode::None)
}

/// Like `check`, except that `ToOdd` (round half to odd) is compared with the `ToEven`
/// reference when that one overflowed: the common reference emulates `ToOdd` at `p + 1`
/// bits, which cannot overflow.
fn check_ovf(
    got: &ExactNum,
    p: usize,
    rm: RoundingMode,
    info: &str,
    cc: &mut Consts,
    f: impl Fn(&mut Float, rnd_t) -> i32,
) {
    if rm == RoundingMode::ToOdd {
        let (rn, _) = mpfr_reference_with(p as u32, RoundingMode::ToEven, &f);
        if rn.is_infinite() {
            check(got, p, RoundingMode::ToEven, info, cc, f);
            return;
        }
    }
    check(got, p, rm, info, cc, f);
}

/// `conv_to_mpfr` for values near the exponent limit: converts `x / 2^e`, then scales.
fn conv_big(w: usize, x: &ExactNum, cc: &mut Consts) -> Float {
    let e = x.exponent().unwrap_or(0);
    let s = x.clone().ldexp(-e, w, RoundingMode::None);
    let f = conv_to_mpfr(w, &s, cc);
    // MPFR's default exponent range is narrower than zenith's.
    let old = unsafe { mpfr::get_emax() };
    unsafe { mpfr::set_emax(mpfr::get_emax_max()) };
    let f = f << e;
    unsafe { mpfr::set_emax(old) };
    f
}

type Binop = fn(&ExactNum, &ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum;
type Mpbin = fn(&mut Float, &Float, &Float, rnd_t) -> i32;

#[test]
fn overflow_by_mode_binary() {
    let ops: [(&str, Binop, Mpbin); 5] = [
        (
            "add",
            |a, b, p, rm, _| ExactNum::add_ref(a, b, p, rm),
            |r, a, b, rnd| unsafe { mpfr::add(r.as_raw_mut(), a.as_raw(), b.as_raw(), rnd) },
        ),
        (
            "mul",
            |a, b, p, rm, _| ExactNum::mul_ref(a, b, p, rm),
            |r, a, b, rnd| unsafe { mpfr::mul(r.as_raw_mut(), a.as_raw(), b.as_raw(), rnd) },
        ),
        (
            "div",
            |a, b, p, rm, _| ExactNum::div_ref(a, b, p, rm),
            |r, a, b, rnd| unsafe { mpfr::div(r.as_raw_mut(), a.as_raw(), b.as_raw(), rnd) },
        ),
        (
            "hypot",
            |a, b, p, rm, _| a.hypot(b, p, rm),
            |r, a, b, rnd| unsafe { mpfr::hypot(r.as_raw_mut(), a.as_raw(), b.as_raw(), rnd) },
        ),
        (
            "pow",
            |a, b, p, rm, cc| {
                if a.is_negative() {
                    neg_of(a).pow(b, p, rm, cc)
                } else {
                    a.pow(b, p, rm, cc)
                }
            },
            |r, a, b, rnd| unsafe {
                let aa = Float::with_val(a.prec(), a.abs_ref());
                mpfr::pow(r.as_raw_mut(), aa.as_raw(), b.as_raw(), rnd)
            },
        ),
    ];
    reset_test_rng();
    let mut cc = Consts::new().unwrap();
    for p in PRECISIONS {
        let w = p.div_ceil(WORD_BIT_SIZE) * WORD_BIT_SIZE;
        for _ in 0..3 {
            let x = input(p, EXPONENT_MAX - 1, EXPONENT_MAX, true);
            let y = input(p, EXPONENT_MAX - 1, EXPONENT_MAX, true);
            let small = input(p, -3, -1, true);
            let two = input(p, 3, 4, false);
            for (name, zf, mp) in ops.iter() {
                let b = match *name {
                    "div" => &small,
                    "pow" => &two,
                    _ => &y,
                };
                let fx = conv_big(w, &x, &mut cc);
                let fb = conv_big(w, b, &mut cc);
                for rm in MODES {
                    let got = zf(&x, b, p, rm, &mut cc);
                    let info = format!("{name} p={p} {rm:?} {x:?} {b:?}");
                    check_ovf(&got, p, rm, &info, &mut cc, |r, rnd| mp(r, &fx, &fb, rnd));
                }
            }
            // fma(x, y, x), ldexp(x, 4), powi(x, 3)
            let fx = conv_big(w, &x, &mut cc);
            let fy = conv_big(w, &y, &mut cc);
            for rm in MODES {
                let got = x.fma(&y, &x, p, rm);
                check_ovf(
                    &got,
                    p,
                    rm,
                    &format!("fma p={p} {rm:?} {x:?}"),
                    &mut cc,
                    |r, rnd| unsafe {
                        mpfr::fma(r.as_raw_mut(), fx.as_raw(), fy.as_raw(), fx.as_raw(), rnd)
                    },
                );
                let got = x.clone().ldexp(4, p, rm);
                check_ovf(
                    &got,
                    p,
                    rm,
                    &format!("ldexp p={p} {rm:?} {x:?}"),
                    &mut cc,
                    |r, rnd| unsafe { mpfr::mul_2si(r.as_raw_mut(), fx.as_raw(), 4, rnd) },
                );
                let got = x.powi(3, p, rm);
                check_ovf(
                    &got,
                    p,
                    rm,
                    &format!("powi p={p} {rm:?} {x:?}"),
                    &mut cc,
                    |r, rnd| unsafe { mpfr::pow_ui(r.as_raw_mut(), fx.as_raw(), 3, rnd) },
                );
            }
        }
    }
}

type Unop = fn(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum;
type Mpun = fn(&mut Float, &Float, rnd_t) -> i32;

#[test]
fn overflow_by_mode_unary() {
    let ops: [(&str, i32, i32, bool, Unop, Mpun); 7] = [
        (
            "exp",
            32,
            33,
            false,
            |x, p, rm, cc| x.exp(p, rm, cc),
            |r, a, rnd| unsafe { mpfr::exp(r.as_raw_mut(), a.as_raw(), rnd) },
        ),
        (
            "exp2",
            32,
            33,
            false,
            |x, p, rm, cc| x.exp2(p, rm, cc),
            |r, a, rnd| unsafe { mpfr::exp2(r.as_raw_mut(), a.as_raw(), rnd) },
        ),
        (
            "exp10",
            32,
            33,
            false,
            |x, p, rm, cc| x.exp10(p, rm, cc),
            |r, a, rnd| unsafe { mpfr::exp10(r.as_raw_mut(), a.as_raw(), rnd) },
        ),
        (
            "sinh",
            32,
            33,
            true,
            |x, p, rm, cc| x.sinh(p, rm, cc),
            |r, a, rnd| unsafe { mpfr::sinh(r.as_raw_mut(), a.as_raw(), rnd) },
        ),
        (
            "cosh",
            32,
            33,
            true,
            |x, p, rm, cc| x.cosh(p, rm, cc),
            |r, a, rnd| unsafe { mpfr::cosh(r.as_raw_mut(), a.as_raw(), rnd) },
        ),
        (
            "ei",
            32,
            33,
            false,
            |x, p, rm, cc| x.ei(p, rm, cc),
            |r, a, rnd| unsafe { mpfr::eint(r.as_raw_mut(), a.as_raw(), rnd) },
        ),
        (
            "gamma",
            29,
            30,
            false,
            |x, p, rm, cc| x.gamma(p, rm, cc),
            |r, a, rnd| unsafe { mpfr::gamma(r.as_raw_mut(), a.as_raw(), rnd) },
        ),
    ];
    reset_test_rng();
    let mut cc = Consts::new().unwrap();
    for p in PRECISIONS {
        let w = p.div_ceil(WORD_BIT_SIZE) * WORD_BIT_SIZE;
        for (name, elo, ehi, signed, zf, mp) in ops.iter() {
            for _ in 0..2 {
                let x = input(p, *elo, *ehi, *signed);
                let fx = conv_big(w, &x, &mut cc);
                for rm in MODES {
                    let got = zf(&x, p, rm, &mut cc);
                    let info = format!("{name} p={p} {rm:?} {x:?}");
                    check_ovf(&got, p, rm, &info, &mut cc, |r, rnd| mp(r, &fx, rnd));
                }
            }
        }
    }
}
