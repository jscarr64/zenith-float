//! Seeded `ToOdd` and `None` passes for the call-site families of `get_random_rnd_pair` (common.rs). The picker draws
//! only the five modes MPFR has a `rnd_t` for, and its `% 5` streams stay as they are. Each family gets one extra
//! pass here with its own seed, through [`extra_mode_reference`]:
//! - `None`: MPFR toward zero at rp(p, 64) (the storage precision), with the inexact flag from the ternary value.
//! - `ToOdd`: nearest with ties to odd, as the library defines it (defs.rs `ToOdd`), derived from MPFR toward zero at
//!   p + 1 bits and the ternary value (`mpfr_reference_with`). Normal-range results only.
//!
//! A family that fails is left out of these passes unless the failure is F-A, F-B or F-C shaped. Any other failure
//! is a new finding, reported and not listed as an exception.
//!
//! Families: ops, special_fn (positive arguments), const, special values, complex add. Not here:
//! - copysign: NEW FINDING N-4. At p = 65 and 100 it rounds at rp(p, 64), not p, in every mode (round_p in
//!   ops/sign.rs), so the family is stopped, not listed.
//! - complex mul and the complex functions: no MPFR oracle in that family for a single rounding (it checks mul with a
//!   tolerance), so the add pass carries the family.

use crate::mpfr::common::{
    assert_exact_eq, conv_to_mpfr, draw_seed, mpfr_reference_with, normal_from_seed,
    reset_test_rng, test_random,
};
use gmp_mpfr_sys::mpfr::{self, rnd_t};
use rug::Float;
use zenith_float_num::{Consts, ExactNum, RoundingMode};

type Res = Result<(), Box<dyn std::error::Error>>;

/// The modes this file adds: the two `get_random_rnd_pair` never draws.
const EXTRA_MODES: [RoundingMode; 2] = [RoundingMode::ToOdd, RoundingMode::None];

/// `f(result, rnd)` rounded as `rm` at `p` bits for a mode with no MPFR `rnd_t` (`ToOdd`, `None`), with whether it is
/// inexact. Any other mode is an error: the five MPFR modes go through `get_random_rnd_pair`.
pub fn extra_mode_reference<F>(
    p: usize,
    rm: RoundingMode,
    f: F,
) -> Result<(Float, bool), Box<dyn std::error::Error>>
where
    F: Fn(&mut Float, rnd_t) -> i32,
{
    match rm {
        RoundingMode::ToOdd | RoundingMode::None => mpfr_reference_with(u32::try_from(p)?, rm, f),
        _ => Err(format!("{rm:?} has an MPFR rnd_t: use get_random_rnd_pair").into()),
    }
}

/// `m / d` as MPFR computes it, for the helper tests (exact when `d` is a power of two).
fn ratio(m: u64, d: u64) -> impl Fn(&mut Float, rnd_t) -> i32 {
    move |r, rnd| {
        let (a, b) = (Float::with_val(64, m), Float::with_val(64, d));
        unsafe { mpfr::div(r.as_raw_mut(), a.as_raw(), b.as_raw(), rnd) }
    }
}

fn is_val(f: &Float, m: u64, d: u64) -> bool {
    let mut w = Float::with_val(256, m);
    w /= d;
    *f == w
}

#[test]
fn helper_exact_is_exact_in_both_modes() -> Res {
    for rm in EXTRA_MODES {
        let (f, ix) = extra_mode_reference(4, rm, ratio(3, 4))?;
        assert!(is_val(&f, 3, 4) && !ix, "{rm:?}: {f}");
    }
    Ok(())
}

#[test]
fn helper_none_truncates_at_the_storage_precision() -> Res {
    // 1/3 at p = 65: storage 128 bits, toward zero, inexact (not the 65-bit truncation)
    let (f, ix) = extra_mode_reference(65, RoundingMode::None, ratio(1, 3))?;
    let mut want = Float::new(128);
    unsafe {
        mpfr::div_ui(
            want.as_raw_mut(),
            Float::with_val(64, 1).as_raw(),
            3,
            rnd_t::RNDZ,
        )
    };
    assert!(ix && f == want && f.prec() >= 128, "None 1/3 at 65: {f}");
    Ok(())
}

#[test]
fn helper_to_odd_tie_with_lsb_already_1_stays() -> Res {
    // 19/16 = 1.0011b at p = 4: a tie between 1.001b (odd) and 1.010b; ToOdd keeps 1.001b = 9/8
    let (f, ix) = extra_mode_reference(4, RoundingMode::ToOdd, ratio(19, 16))?;
    assert!(is_val(&f, 9, 8) && ix, "ToOdd 19/16: {f}");
    Ok(())
}

#[test]
fn helper_to_odd_tie_with_lsb_0_goes_to_the_odd_neighbour() -> Res {
    // 17/16 = 1.0001b at p = 4: a tie between 1.000b (even) and 1.001b (odd); ToOdd gives 9/8
    let (f, ix) = extra_mode_reference(4, RoundingMode::ToOdd, ratio(17, 16))?;
    assert!(is_val(&f, 9, 8) && ix, "ToOdd 17/16: {f}");
    // 33/32 = 1.00001b is below the tie: nearest gives 1.000b (a sticky round-to-odd would give 1.001b)
    let (f, ix) = extra_mode_reference(4, RoundingMode::ToOdd, ratio(33, 32))?;
    assert!(is_val(&f, 1, 1) && ix, "ToOdd 33/32: {f}");
    Ok(())
}

#[test]
fn helper_binade_edge() -> Res {
    // 2 - 2^-200 at p = 64: None truncates to 2 - 2^-63 (all ones), ToOdd rounds to nearest, 2
    let below_two = |r: &mut Float, rnd: rnd_t| {
        let mut tiny = Float::with_val(64, 1);
        tiny >>= 200;
        let two = Float::with_val(64, 2);
        unsafe { mpfr::sub(r.as_raw_mut(), two.as_raw(), tiny.as_raw(), rnd) }
    };
    let (f, ix) = extra_mode_reference(64, RoundingMode::None, below_two)?;
    let mut ones = Float::with_val(256, 2);
    ones -= Float::with_val(64, 1) >> 63;
    assert!(ix && f == ones, "None 2 - 2^-200: {f}");
    let (f, ix) = extra_mode_reference(64, RoundingMode::ToOdd, below_two)?;
    assert!(ix && f == 2, "ToOdd 2 - 2^-200: {f}");
    // the five MPFR modes are refused
    assert!(extra_mode_reference(64, RoundingMode::ToEven, below_two).is_err());
    Ok(())
}

/// The precisions a pass draws from: word multiples (storage = p) and not (storage = rp(p, 64) > p).
fn draw_p() -> usize {
    [64usize, 65, 100, 128, 192, 200, 256, 320][test_random::<usize>() % 8]
}

fn check_exact(n: &ExactNum, (f, ix): &(Float, bool), p: usize, msg: &str, cc: &mut Consts) -> Res {
    assert_exact_eq(n, f, p, msg, cc)?;
    if !n.is_nan() && !n.is_inf() {
        assert_eq!(n.inexact(), *ix, "{msg}: inexact flag");
    }
    Ok(())
}

type Un = fn(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum;
type Bin = fn(&ExactNum, &ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum;
type MpUn = fn(&mut Float, &Float, rnd_t) -> i32;
type MpBin = fn(&mut Float, &Float, &Float, rnd_t) -> i32;

/// One seeded pass: `runs` draws of `p` and operands in `[2^lo, 2^hi]` (both signs), every unary and binary op, both
/// extra modes, bit-exact against [`extra_mode_reference`] with the inexact flag.
fn pass(
    runs: usize,
    lo: i64,
    hi: i64,
    positive: bool,
    un: &[(&str, Un, MpUn)],
    bin: &[(&str, Bin, MpBin)],
) -> Res {
    reset_test_rng();
    let mut cc = Consts::new()?;
    for _ in 0..runs {
        let p = draw_p();
        // operands at a word multiple (a random normal fills its words), results at `p`
        let po = 64 * (1 + test_random::<usize>() % 4);
        let (s1, s2) = (draw_seed(), draw_seed());
        let mut n1 = normal_from_seed(s1, po, lo, hi);
        let mut n2 = normal_from_seed(s2, po, lo, hi);
        if positive {
            n1 = n1.abs();
            n2 = n2.abs();
        }
        let (f1, f2) = (
            conv_to_mpfr(po, &n1, &mut cc)?,
            conv_to_mpfr(po, &n2, &mut cc)?,
        );
        for rm in EXTRA_MODES {
            for (name, zf, mp) in un {
                let r = extra_mode_reference(p, rm, |r, d| mp(r, &f1, d))?;
                let msg = format!("{name} p={p} {rm:?} x={n1:?}");
                check_exact(&zf(&n1, p, rm, &mut cc), &r, p, &msg, &mut cc)?;
            }
            for (name, zf, mp) in bin {
                let r = extra_mode_reference(p, rm, |r, d| mp(r, &f1, &f2, d))?;
                let msg = format!("{name} p={p} {rm:?} x={n1:?} y={n2:?}");
                check_exact(&zf(&n1, &n2, p, rm, &mut cc), &r, p, &msg, &mut cc)?;
            }
        }
    }
    Ok(())
}

/// compare_ops_test family (`get_random_rnd_pair` there: round, add, sub, mul, div, pow_ui, hypot, sin_cos,
/// sinh_cosh and the elementary functions): one extra seeded pass in ToOdd and None.
#[test]
fn ops_family_extra_modes() -> Res {
    let un: [(&str, Un, MpUn); 8] = [
        (
            "sqrt",
            |x, p, rm, _| x.sqrt(p, rm),
            |r, a, d| unsafe { mpfr::sqrt(r.as_raw_mut(), a.as_raw(), d) },
        ),
        (
            "sin",
            |x, p, rm, cc| x.sin(p, rm, cc),
            |r, a, d| unsafe { mpfr::sin(r.as_raw_mut(), a.as_raw(), d) },
        ),
        (
            "cos",
            |x, p, rm, cc| x.cos(p, rm, cc),
            |r, a, d| unsafe { mpfr::cos(r.as_raw_mut(), a.as_raw(), d) },
        ),
        (
            "sinh",
            |x, p, rm, cc| x.sinh(p, rm, cc),
            |r, a, d| unsafe { mpfr::sinh(r.as_raw_mut(), a.as_raw(), d) },
        ),
        (
            "cosh",
            |x, p, rm, cc| x.cosh(p, rm, cc),
            |r, a, d| unsafe { mpfr::cosh(r.as_raw_mut(), a.as_raw(), d) },
        ),
        (
            "exp",
            |x, p, rm, cc| x.exp(p, rm, cc),
            |r, a, d| unsafe { mpfr::exp(r.as_raw_mut(), a.as_raw(), d) },
        ),
        (
            "ln",
            |x, p, rm, cc| x.ln(p, rm, cc),
            |r, a, d| unsafe { mpfr::log(r.as_raw_mut(), a.as_raw(), d) },
        ),
        (
            "powi3",
            |x, p, rm, _| x.powi(3, p, rm),
            |r, a, d| unsafe { mpfr::pow_ui(r.as_raw_mut(), a.as_raw(), 3, d) },
        ),
    ];
    let bin: [(&str, Bin, MpBin); 5] = [
        (
            "add",
            |x, y, p, rm, _| ExactNum::add_ref(x, y, p, rm),
            |r, a, b, d| unsafe { mpfr::add(r.as_raw_mut(), a.as_raw(), b.as_raw(), d) },
        ),
        (
            "sub",
            |x, y, p, rm, _| ExactNum::sub_ref(x, y, p, rm),
            |r, a, b, d| unsafe { mpfr::sub(r.as_raw_mut(), a.as_raw(), b.as_raw(), d) },
        ),
        (
            "mul",
            |x, y, p, rm, _| ExactNum::mul_ref(x, y, p, rm),
            |r, a, b, d| unsafe { mpfr::mul(r.as_raw_mut(), a.as_raw(), b.as_raw(), d) },
        ),
        (
            "div",
            |x, y, p, rm, _| ExactNum::div_ref(x, y, p, rm),
            |r, a, b, d| unsafe { mpfr::div(r.as_raw_mut(), a.as_raw(), b.as_raw(), d) },
        ),
        (
            "hypot",
            |x, y, p, rm, _| x.hypot(y, p, rm),
            |r, a, b, d| unsafe { mpfr::hypot(r.as_raw_mut(), a.as_raw(), b.as_raw(), d) },
        ),
    ];
    pass(400, -20, 20, false, &un, &bin)
}

/// compare_special_fn_test family (erf, erfc, gamma, lngamma, digamma, ai there), on positive arguments: `ln_gamma`
/// is documented for positive arguments only.
#[test]
fn special_fn_family_extra_modes() -> Res {
    let un: [(&str, Un, MpUn); 6] = [
        (
            "erf",
            |x, p, rm, cc| x.erf(p, rm, cc),
            |r, a, d| unsafe { mpfr::erf(r.as_raw_mut(), a.as_raw(), d) },
        ),
        (
            "erfc",
            |x, p, rm, cc| x.erfc(p, rm, cc),
            |r, a, d| unsafe { mpfr::erfc(r.as_raw_mut(), a.as_raw(), d) },
        ),
        (
            "gamma",
            |x, p, rm, cc| x.gamma(p, rm, cc),
            |r, a, d| unsafe { mpfr::gamma(r.as_raw_mut(), a.as_raw(), d) },
        ),
        (
            "ln_gamma",
            |x, p, rm, cc| x.ln_gamma(p, rm, cc),
            |r, a, d| unsafe { mpfr::lngamma(r.as_raw_mut(), a.as_raw(), d) },
        ),
        (
            "digamma",
            |x, p, rm, cc| x.digamma(p, rm, cc),
            |r, a, d| unsafe { mpfr::digamma(r.as_raw_mut(), a.as_raw(), d) },
        ),
        (
            "ai",
            |x, p, rm, cc| x.ai(p, rm, cc),
            |r, a, d| unsafe { mpfr::ai(r.as_raw_mut(), a.as_raw(), d) },
        ),
    ];
    pass(60, -3, 3, true, &un, &[])
}

/// compare_const_test family (pi, ln 2, e there).
#[test]
fn const_family_extra_modes() -> Res {
    reset_test_rng();
    let mut cc = Consts::new()?;
    for _ in 0..60 {
        let p = draw_p();
        for rm in EXTRA_MODES {
            let r =
                extra_mode_reference(p, rm, |r, d| unsafe { mpfr::const_pi(r.as_raw_mut(), d) })?;
            check_exact(&cc.pi(p, rm), &r, p, &format!("pi p={p} {rm:?}"), &mut cc)?;
            let r =
                extra_mode_reference(p, rm, |r, d| unsafe { mpfr::const_log2(r.as_raw_mut(), d) })?;
            check_exact(
                &cc.ln_2(p, rm),
                &r,
                p,
                &format!("ln2 p={p} {rm:?}"),
                &mut cc,
            )?;
            let one = Float::with_val(64, 1);
            let r = extra_mode_reference(p, rm, |r, d| unsafe {
                mpfr::exp(r.as_raw_mut(), one.as_raw(), d)
            })?;
            check_exact(&cc.e(p, rm), &r, p, &format!("e p={p} {rm:?}"), &mut cc)?;
        }
    }
    Ok(())
}

/// compare_special_test family: the special-value operands drawn there (zeros, +-1, tiny and large normals), add, sub,
/// mul, div and ln. No pair overflows: a None overflow is NONE-GAP F-B, queued for ZF-M27.
#[test]
fn special_values_family_extra_modes() -> Res {
    let mut cc = Consts::new()?;
    let neg = |mut x: ExactNum| {
        x.inv_sign();
        x
    };
    let mk = |i: usize, p: usize| match i {
        0 => ExactNum::new(p),
        1 => neg(ExactNum::new(p)),
        2 => ExactNum::from_i8(1, p),
        3 => ExactNum::from_i8(-1, p),
        4 => ExactNum::from_u8(3, p).ldexp(-1000, p, RoundingMode::None),
        5 => ExactNum::from_u8(5, p).ldexp(1000, p, RoundingMode::None),
        _ => neg(ExactNum::from_u8(7, p).ldexp(-3, p, RoundingMode::None)),
    };
    for po in [64usize, 128] {
        for p in [64usize, 65, 128, 200] {
            for i in 0..7 {
                for j in 0..7 {
                    let (x, y) = (mk(i, po), mk(j, po));
                    let (fx, fy) = (
                        conv_to_mpfr(po, &x, &mut cc)?,
                        conv_to_mpfr(po, &y, &mut cc)?,
                    );
                    for rm in EXTRA_MODES {
                        let msg = |op: &str| format!("{op} p={p} {rm:?} i={i} j={j}");
                        let r = extra_mode_reference(p, rm, |r, d| unsafe {
                            mpfr::add(r.as_raw_mut(), fx.as_raw(), fy.as_raw(), d)
                        })?;
                        check_exact(
                            &ExactNum::add_ref(&x, &y, p, rm),
                            &r,
                            p,
                            &msg("add"),
                            &mut cc,
                        )?;
                        let r = extra_mode_reference(p, rm, |r, d| unsafe {
                            mpfr::sub(r.as_raw_mut(), fx.as_raw(), fy.as_raw(), d)
                        })?;
                        check_exact(
                            &ExactNum::sub_ref(&x, &y, p, rm),
                            &r,
                            p,
                            &msg("sub"),
                            &mut cc,
                        )?;
                        let r = extra_mode_reference(p, rm, |r, d| unsafe {
                            mpfr::mul(r.as_raw_mut(), fx.as_raw(), fy.as_raw(), d)
                        })?;
                        check_exact(
                            &ExactNum::mul_ref(&x, &y, p, rm),
                            &r,
                            p,
                            &msg("mul"),
                            &mut cc,
                        )?;
                        let r = extra_mode_reference(p, rm, |r, d| unsafe {
                            mpfr::div(r.as_raw_mut(), fx.as_raw(), fy.as_raw(), d)
                        })?;
                        check_exact(
                            &ExactNum::div_ref(&x, &y, p, rm),
                            &r,
                            p,
                            &msg("div"),
                            &mut cc,
                        )?;
                    }
                }
                let x = mk(i, po);
                let fx = conv_to_mpfr(po, &x, &mut cc)?;
                for rm in EXTRA_MODES {
                    let r = extra_mode_reference(p, rm, |r, d| unsafe {
                        mpfr::log(r.as_raw_mut(), fx.as_raw(), d)
                    })?;
                    check_exact(
                        &x.ln(p, rm, &mut cc),
                        &r,
                        p,
                        &format!("ln p={p} {rm:?} i={i}"),
                        &mut cc,
                    )?;
                }
            }
        }
    }
    Ok(())
}

/// compare_complex_test family, add: componentwise MPFR, as the family checks it.
#[test]
fn complex_add_family_extra_modes() -> Res {
    use zenith_float_num::ExactComplex;
    reset_test_rng();
    let mut cc = Consts::new()?;
    for _ in 0..100 {
        let p = draw_p();
        let s: [u64; 4] = [draw_seed(), draw_seed(), draw_seed(), draw_seed()];
        let part = |k: usize| normal_from_seed(s[k], 128, -8, 8);
        let f: Vec<Float> = (0..4)
            .map(|k| conv_to_mpfr(128, &part(k), &mut cc))
            .collect::<Result<_, _>>()?;
        for rm in EXTRA_MODES {
            let a = ExactComplex::new(part(0), part(1));
            let b = ExactComplex::new(part(2), part(3));
            let c = ExactComplex::add_ref(&a, &b, p, rm);

            let rr = extra_mode_reference(p, rm, |r, d| unsafe {
                mpfr::add(r.as_raw_mut(), f[0].as_raw(), f[2].as_raw(), d)
            })?;
            let ri = extra_mode_reference(p, rm, |r, d| unsafe {
                mpfr::add(r.as_raw_mut(), f[1].as_raw(), f[3].as_raw(), d)
            })?;
            // the flag is value-wide: ExactComplex::inexact is "this operation rounded either part"
            let msg = format!("cplx add p={p} {rm:?}");
            assert_exact_eq(c.re(), &rr.0, p, &msg, &mut cc)?;
            assert_exact_eq(c.im(), &ri.0, p, &msg, &mut cc)?;
            assert_eq!(c.inexact(), rr.1 || ri.1, "{msg}: inexact flag");
        }
    }
    Ok(())
}
