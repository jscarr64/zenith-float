//! Two- and three-operand functions, roots, powers, and parsing round at exactly `p` bits in
//! every mode with a per-operation `inexact` flag (2.0.0 audit, round 2: items a and c).

use crate::mpfr::common::assert_exact_eq;
use crate::mpfr::common::conv_to_mpfr;
use crate::mpfr::common::reset_test_rng;
use crate::mpfr::exact_bit_functions::{
    check, exact_input, input, MODES, P14_PRECISIONS, PRECISIONS,
};
use gmp_mpfr_sys::mpfr;
use rug::Float;
use std::ops::Neg;
use zenith_float_num::{Consts, ExactNum, Radix, WORD_BIT_SIZE};

fn storage(p: usize) -> usize {
    p.div_ceil(WORD_BIT_SIZE) * WORD_BIT_SIZE
}

fn pairs(
    p: usize,
    a: (i32, i32, bool),
    b: (i32, i32, bool),
    exact: &[(i64, i32, i64, i32)],
) -> Result<Vec<(ExactNum, ExactNum)>, Box<dyn std::error::Error>> {
    let mut v = Vec::new();
    for _ in 0..6 {
        v.push((input(p, a.0, a.1, a.2)?, input(p, b.0, b.1, b.2)?));
    }
    for &(n1, k1, n2, k2) in exact {
        v.push((exact_input(n1, k1, p), exact_input(n2, k2, p)));
    }
    Ok(v)
}

#[test]
fn exact_bit_hypot() -> Result<(), Box<dyn std::error::Error>> {
    reset_test_rng();
    let mut cc = Consts::new()?;
    for p in PRECISIONS {
        for (x, y) in pairs(
            p,
            (-20, 20, true),
            (-20, 20, true),
            &[(3, 0, 4, 0), (5, 0, 12, 0), (1, 0, 1, 0)],
        )? {
            let (fx, fy) = (
                conv_to_mpfr(storage(p), &x, &mut cc)?,
                conv_to_mpfr(storage(p), &y, &mut cc)?,
            );
            for rm in MODES {
                let got = x.hypot(&y, p, rm);
                let info = format!("hypot p={p} {rm:?} {x:?} {y:?}");
                check(&got, p, rm, &info, &mut cc, |r, rnd| unsafe {
                    mpfr::hypot(r.as_raw_mut(), fx.as_raw(), fy.as_raw(), rnd)
                })?;
            }
        }
    }
    Ok(())
}

#[test]
fn exact_bit_fma() -> Result<(), Box<dyn std::error::Error>> {
    reset_test_rng();
    let mut cc = Consts::new()?;
    for p in PRECISIONS {
        for (i, (x, y)) in pairs(
            p,
            (-20, 20, true),
            (-20, 20, true),
            &[(3, 0, 4, 0), (1, 1, 1, 1)],
        )?
        .into_iter()
        .enumerate()
        {
            // c cancels a*b exactly for the last exact case; otherwise random
            let z = if i % 3 == 0 { x.mul_full_prec(&y).neg() } else { input(p, -30, 30, true)? };
            let mut z = z;
            z.set_inexact(true);
            let zs = z.mantissa_max_bit_len().unwrap_or(storage(p));
            let (fx, fy, fz) = (
                conv_to_mpfr(storage(p), &x, &mut cc)?,
                conv_to_mpfr(storage(p), &y, &mut cc)?,
                conv_to_mpfr(zs, &z, &mut cc)?,
            );
            for rm in MODES {
                let got = x.fma(&y, &z, p, rm);
                let info = format!("fma p={p} {rm:?} {x:?} {y:?} {z:?}");
                check(&got, p, rm, &info, &mut cc, |r, rnd| unsafe {
                    mpfr::fma(r.as_raw_mut(), fx.as_raw(), fy.as_raw(), fz.as_raw(), rnd)
                })?;
            }
        }
    }
    Ok(())
}

#[test]
fn exact_bit_pow_atan2() -> Result<(), Box<dyn std::error::Error>> {
    reset_test_rng();
    let mut cc = Consts::new()?;
    for p in PRECISIONS {
        for (x, y) in pairs(
            p,
            (-4, 4, false),
            (-3, 3, true),
            &[(2, 0, 3, 0), (9, 0, 1, 1), (4, 0, -1, 1), (1, 0, 7, 0)],
        )? {
            let (fx, fy) = (
                conv_to_mpfr(storage(p), &x, &mut cc)?,
                conv_to_mpfr(storage(p), &y, &mut cc)?,
            );
            for rm in MODES {
                let got = x.pow(&y, p, rm, &mut cc);
                let info = format!("pow p={p} {rm:?} {x:?} {y:?}");
                check(&got, p, rm, &info, &mut cc, |r, rnd| unsafe {
                    mpfr::pow(r.as_raw_mut(), fx.as_raw(), fy.as_raw(), rnd)
                })?;
                let got = y.atan2(&x, p, rm, &mut cc);
                let info = format!("atan2 p={p} {rm:?} {y:?} {x:?}");
                check(&got, p, rm, &info, &mut cc, |r, rnd| unsafe {
                    mpfr::atan2(r.as_raw_mut(), fy.as_raw(), fx.as_raw(), rnd)
                })?;
            }
        }
    }
    Ok(())
}

#[test]
fn exact_bit_powi_root_reciprocal() -> Result<(), Box<dyn std::error::Error>> {
    reset_test_rng();
    let mut cc = Consts::new()?;
    for p in PRECISIONS {
        let mut xs: Vec<ExactNum> = (0..6)
            .map(|_| input(p, -6, 6, true))
            .collect::<Result<_, _>>()?;
        for (n, k) in [(3, 0), (-2, 0), (1, 1), (27, 0), (16, 0)] {
            xs.push(exact_input(n, k, p));
        }
        for x in &xs {
            let fx = conv_to_mpfr(storage(p), x, &mut cc)?;
            for rm in MODES {
                for n in [0usize, 1, 2, 3, 7] {
                    let got = x.powi(n, p, rm);
                    let info = format!("powi({n}) p={p} {rm:?} {x:?}");
                    check(&got, p, rm, &info, &mut cc, |r, rnd| unsafe {
                        mpfr::pow_ui(r.as_raw_mut(), fx.as_raw(), n as _, rnd)
                    })?;
                }
                for n in [-3i64, -1, 2] {
                    let got = x.powsi(n, p, rm);
                    let info = format!("powsi({n}) p={p} {rm:?} {x:?}");
                    check(&got, p, rm, &info, &mut cc, |r, rnd| unsafe {
                        mpfr::pow_si(r.as_raw_mut(), fx.as_raw(), n as _, rnd)
                    })?;
                }
                let got = x.reciprocal(p, rm);
                let info = format!("reciprocal p={p} {rm:?} {x:?}");
                check(&got, p, rm, &info, &mut cc, |r, rnd| unsafe {
                    mpfr::ui_div(r.as_raw_mut(), 1, fx.as_raw(), rnd)
                })?;
                if !x.is_negative() {
                    for n in [4usize, 5] {
                        let got = x.nth_root(n, p, rm);
                        let info = format!("nth_root({n}) p={p} {rm:?} {x:?}");
                        check(&got, p, rm, &info, &mut cc, |r, rnd| unsafe {
                            mpfr::rootn_ui(r.as_raw_mut(), fx.as_raw(), n as _, rnd)
                        })?;
                    }
                }
            }
        }
    }
    Ok(())
}

#[test]
fn exact_bit_parse() -> Result<(), Box<dyn std::error::Error>> {
    let mut cc = Consts::new()?;
    let decimals = [
        "0.1",
        "-0.1",
        "1.5",
        "3",
        "1e10",
        "-2.5e-7",
        "123456789.123456789",
        "0.333333333333333333333",
        "6.02214076e23",
        "1e-300",
        "-7.25",
        "0.0",
    ];
    for p in PRECISIONS {
        for s in decimals {
            for rm in MODES {
                let got = ExactNum::parse(s, Radix::Dec, p, rm, &mut cc);
                let info = format!("parse({s}) p={p} {rm:?}");
                let c = std::ffi::CString::new(s)?;
                let exact = Float::with_val(4096, Float::parse(s)?);
                check(&got, p, rm, &info, &mut cc, |r, rnd| {
                    unsafe { mpfr::set_str(r.as_raw_mut(), c.as_ptr(), 10, rnd) };
                    // set_str's return value is not the ternary value: recompute it
                    match (*r).partial_cmp(&exact) {
                        Some(core::cmp::Ordering::Less) => -1,
                        Some(core::cmp::Ordering::Greater) => 1,
                        _ => 0,
                    }
                })?;
            }
        }
    }
    Ok(())
}

#[test]
fn exact_bit_hypot_wide_spread() -> Result<(), Box<dyn std::error::Error>> {
    reset_test_rng();
    let mut cc = Consts::new()?;
    for p in [53usize, 128, 8512] {
        for _ in 0..4 {
            let x = input(p, -20, 20, true)?;
            let y = input(p, -400_000, -300_000, true)?;
            for (a, b) in [(&x, &y), (&y, &x)] {
                let (fa, fb) = (
                    conv_to_mpfr(storage(p), a, &mut cc)?,
                    conv_to_mpfr(storage(p), b, &mut cc)?,
                );
                for rm in MODES {
                    let got = a.hypot(b, p, rm);
                    let info = format!("hypot spread p={p} {rm:?}");
                    check(&got, p, rm, &info, &mut cc, |r, rnd| unsafe {
                        mpfr::hypot(r.as_raw_mut(), fa.as_raw(), fb.as_raw(), rnd)
                    })?;
                }
            }
        }
    }
    Ok(())
}

/// zenith-bound P14: `add_full_prec` and `sub_full_prec` never round. The result equals the exact MPFR sum or
/// difference (computed wide enough to be exact) bit for bit, and `inexact()` is false.
#[test]
fn add_sub_full_prec_are_exact() -> Result<(), Box<dyn std::error::Error>> {
    reset_test_rng();
    let mut cc = Consts::new()?;
    for p in PRECISIONS.into_iter().chain(P14_PRECISIONS) {
        let mut v = Vec::new();
        for _ in 0..8 {
            v.push((input(p, -40, 40, true)?, input(p, -200, 200, true)?));
        }
        // cancellation, a far-smaller operand, equal magnitudes
        v.push((exact_input(1, 0, p), exact_input(1, 0, p)));
        v.push((exact_input(3, 0, p), exact_input(-3, 0, p)));
        v.push((exact_input(1, 0, p), exact_input(-1, 300, p)));
        v.push((exact_input(-5, 2, p), exact_input(7, 100, p)));
        for (a, b) in &mut v {
            a.set_inexact(false);
            b.set_inexact(false);
            let wide = 2 * storage(p) + 1024;
            let fa = conv_to_mpfr(storage(p), a, &mut cc)?;
            let fb = conv_to_mpfr(storage(p), b, &mut cc)?;
            for (name, got, sub) in [
                ("add_full_prec", a.add_full_prec(b), false),
                ("sub_full_prec", a.sub_full_prec(b), true),
            ] {
                let mut want = Float::new(wide as u32);
                let t = unsafe {
                    if sub {
                        mpfr::sub(
                            want.as_raw_mut(),
                            fa.as_raw(),
                            fb.as_raw(),
                            mpfr::rnd_t::RNDN,
                        )
                    } else {
                        mpfr::add(
                            want.as_raw_mut(),
                            fa.as_raw(),
                            fb.as_raw(),
                            mpfr::rnd_t::RNDN,
                        )
                    }
                };
                let info = format!("{name} p={p} a={a:?} b={b:?}");
                assert_eq!(t, 0, "{info}: reference not exact");
                assert_exact_eq(&got, &want, wide, &info, &mut cc)?;
                assert!(!got.inexact(), "{info}: inexact");
            }
        }
    }
    Ok(())
}
