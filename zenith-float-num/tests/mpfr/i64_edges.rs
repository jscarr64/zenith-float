//! zenith-float 2.0 exponent range: `Exponent = i64`, `EXPONENT_MAX = 2^62 - 1 = -EXPONENT_MIN`
//! (MPFR's `emax_max` / `emin_min` on LP64). Every rounding mode is checked against MPFR within 0, 1,
//! p - 1, p and p + 1 of each limit: carry-out of an all-ones mantissa at `EXPONENT_MAX`, overflow,
//! and subnormals rounding up into the normal range at `EXPONENT_MIN`. Results near the bottom are
//! compared in a space scaled by `2^SHIFT` (see `mpfr_reference_shifted`), because MPFR cannot hold
//! the subnormal grid `2^(EXPONENT_MIN - p)` below its own `emin_min`.

use crate::mpfr::common::{emulated_range_ok, exact_mpfr_shifted, mpfr_reference_shifted};
use crate::mpfr::common::{show_mp, with_full_range};
use crate::mpfr::exact_bit_functions::MODES;
use gmp_mpfr_sys::mpfr::{self, rnd_t};
use rug::Float;
use std::ops::Neg;
use zenith_float_num::{
    ExactNum, Exponent, RoundingMode, EXPONENT_MAX, EXPONENT_MIN, MAX_PRECISION, WORD_BIT_SIZE,
};

type Res = Result<(), Box<dyn std::error::Error>>;

const PS: [usize; 5] = [64, 113, 16384, 16385, 65536];

/// Scale used for results near `EXPONENT_MIN`.
const SHIFT: i64 = 1 << 40;

fn offsets(p: usize) -> [i64; 5] {
    let p = p as i64;
    [0, 1, p - 1, p, p + 1]
}

fn storage(p: usize) -> usize {
    p.div_ceil(WORD_BIT_SIZE) * WORD_BIT_SIZE
}

/// `1 - 2^-p` (all `p` mantissa bits set) times `2^(e - 0)`: `m * 2^e` with `0.5 <= m < 1`.
fn all_ones(p: usize, e: i64) -> Result<ExactNum, Box<dyn std::error::Error>> {
    let one = ExactNum::from_word(1, storage(p));
    let ulp = ExactNum::from_word(1, storage(p)).ldexp(
        -(p as Exponent),
        storage(p),
        RoundingMode::ToEven,
    );
    let mut x = ExactNum::sub_ref(&one, &ulp, p, RoundingMode::ToZero);
    if x.inexact() {
        return Err("1 - 2^-p is not exact".into());
    }
    x.set_exponent(Exponent::try_from(e)?);
    Ok(x)
}

/// `0.5 * 2^e`.
fn half_at(p: usize, e: i64) -> Result<ExactNum, Box<dyn std::error::Error>> {
    let mut x = ExactNum::from_word(1, storage(p));
    x.set_exponent(Exponent::try_from(e)?);
    Ok(x)
}

/// zenith result vs MPFR reference, both in the space scaled by `2^shift`.
fn show_zf(x: &ExactNum) -> String {
    match x.as_raw_parts() {
        Some((w, n, s, e, inx)) => format!("{s:?} {w:x?} n={n} e={e} inexact={inx}"),
        None => format!("special nan={} inf={}", x.is_nan(), x.is_inf()),
    }
}

fn same(got: &ExactNum, want: &Float, inexact: bool, shift: i64, info: &str) -> Res {
    if want.is_infinite() {
        assert!(
            got.is_inf(),
            "{info}: MPFR {}, zenith {}",
            show_mp(want),
            show_zf(got)
        );
        assert_eq!(
            got.is_negative(),
            want.is_sign_negative(),
            "{info}: sign of Inf"
        );
        return Ok(());
    }
    assert!(
        !got.is_nan() && !got.is_inf(),
        "{info}: MPFR {}, zenith {}",
        show_mp(want),
        show_zf(got)
    );
    if want.is_zero() {
        assert!(got.is_zero(), "{info}: MPFR 0, zenith {}", show_zf(got));
        assert_eq!(
            got.is_negative(),
            want.is_sign_negative(),
            "{info}: sign of zero"
        );
    } else {
        let (old_min, old_max) = unsafe { (mpfr::get_emin(), mpfr::get_emax()) };
        unsafe {
            mpfr::set_emin(mpfr::get_emin_min());
            mpfr::set_emax(mpfr::get_emax_max());
        }
        let g = exact_mpfr_shifted(got, shift)?;
        let eq = g == *want;
        unsafe {
            mpfr::set_emin(old_min);
            mpfr::set_emax(old_max);
        }
        assert!(
            eq,
            "{info}: MPFR {}, zenith {}",
            show_mp(want),
            show_zf(got)
        );
    }
    assert_eq!(got.inexact(), inexact, "{info}: inexact flag");
    Ok(())
}

#[test]
fn limits_are_62_bit_and_match_mpfr() {
    assert_eq!(EXPONENT_MAX as i128, (1i128 << 62) - 1);
    assert_eq!(EXPONENT_MIN as i128, -((1i128 << 62) - 1));
    assert_eq!(size_of::<Exponent>(), 8);
    let by_size = (usize::MAX / 4) as u128;
    let lim = (1u128 << 60).min(by_size);
    assert_eq!(MAX_PRECISION as u128, lim - 1);
    unsafe {
        assert_eq!(mpfr::get_emax_max() as i128, EXPONENT_MAX as i128);
        assert_eq!(mpfr::get_emin_min() as i128, EXPONENT_MIN as i128);
    }
    assert!(emulated_range_ok(64, SHIFT));
}

/// x = all-ones at `EXPONENT_MAX - k`: x + half ulp (tie, carries out of the mantissa), x + x, x * 2 and
/// x * x, in every mode.
#[test]
fn top_edge_carry_out_and_overflow() -> Res {
    for p in PS {
        for k in offsets(p) {
            let e = EXPONENT_MAX as i64 - k;
            let two = ExactNum::from_word(2, 64);
            let f2 = Float::with_val(64, 2);
            for neg in [false, true] {
                let x = all_ones(p, e)?;
                let half_ulp = half_at(p, e - p as i64)?;
                let (x, h) = if neg { (x.neg(), half_ulp.neg()) } else { (x, half_ulp) };
                let fx = with_full_range(|| exact_mpfr_shifted(&x, 0))?;
                let fh = with_full_range(|| exact_mpfr_shifted(&h, 0))?;
                for rm in MODES {
                    let pu = p as u32;
                    let cases: [(&str, ExactNum, Box<dyn Fn(&mut Float, rnd_t) -> i32>); 4] = [
                        (
                            "x+h",
                            ExactNum::add_ref(&x, &h, p, rm),
                            Box::new(|r: &mut Float, rnd| unsafe {
                                mpfr::add(r.as_raw_mut(), fx.as_raw(), fh.as_raw(), rnd)
                            }),
                        ),
                        (
                            "x+x",
                            ExactNum::add_ref(&x, &x, p, rm),
                            Box::new(|r: &mut Float, rnd| unsafe {
                                mpfr::add(r.as_raw_mut(), fx.as_raw(), fx.as_raw(), rnd)
                            }),
                        ),
                        (
                            "x*2",
                            ExactNum::mul_ref(&x, &two, p, rm),
                            Box::new(|r: &mut Float, rnd| unsafe {
                                mpfr::mul(r.as_raw_mut(), fx.as_raw(), f2.as_raw(), rnd)
                            }),
                        ),
                        (
                            "x*x",
                            ExactNum::mul_ref(&x, &x, p, rm),
                            Box::new(|r: &mut Float, rnd| unsafe {
                                mpfr::mul(r.as_raw_mut(), fx.as_raw(), fx.as_raw(), rnd)
                            }),
                        ),
                    ];
                    for (name, got, f) in cases {
                        let (want, inexact) = mpfr_reference_shifted(pu, rm, 0, f)?;
                        let info = format!("top {name} p={p} k={k} neg={neg} {rm:?}");
                        same(&got, &want, inexact, 0, &info)?;
                    }
                }
            }
        }
    }
    Ok(())
}

/// Near `EXPONENT_MIN`: an all-ones normal at `EXPONENT_MIN + k` scaled down by `2^-(k + j)` lands on
/// the subnormal grid and rounds (up into the normal range for the all-ones pattern); the largest
/// subnormal plus the smallest positive is the smallest normal; and x - y for neighbouring normals.
#[test]
fn bottom_edge_subnormals_round_into_normal() -> Res {
    for p in PS {
        let pu = p as u32;
        if !emulated_range_ok(pu, SHIFT) {
            return Err(format!("MPFR cannot emulate the range at p={p}").into());
        }
        for k in offsets(p) {
            let e = EXPONENT_MIN as i64 + k;
            let x = all_ones(p, e)?;
            let fx = with_full_range(|| exact_mpfr_shifted(&x, SHIFT))?;
            let fxn = with_full_range(|| exact_mpfr_shifted(&all_ones(p, e)?.neg(), SHIFT))?;
            for j in [1i64, (p / 2) as i64, p as i64 - 1, p as i64, p as i64 + 1, p as i64 + 2] {
                let d = k + j;
                for rm in MODES {
                    let got = all_ones(p, e)?.ldexp(-(d as Exponent), p, rm);
                    let (want, inexact) = mpfr_reference_shifted(pu, rm, SHIFT, |r, rnd| unsafe {
                        mpfr::mul_2si(r.as_raw_mut(), fx.as_raw(), -d as _, rnd)
                    })?;
                    same(
                        &got,
                        &want,
                        inexact,
                        SHIFT,
                        &format!("bottom ldexp p={p} k={k} j={j} {rm:?}"),
                    )?;
                    let got = all_ones(p, e)?.neg().ldexp(-(d as Exponent), p, rm);
                    let (want, inexact) = mpfr_reference_shifted(pu, rm, SHIFT, |r, rnd| unsafe {
                        mpfr::mul_2si(r.as_raw_mut(), fxn.as_raw(), -d as _, rnd)
                    })?;
                    same(
                        &got,
                        &want,
                        inexact,
                        SHIFT,
                        &format!("bottom -ldexp p={p} k={k} j={j} {rm:?}"),
                    )?;
                }
            }
        }
        // min_positive / min_positive_normal round the precision up to whole words.
        let ps = storage(p);
        let psu = ps as u32;
        let tiny = ExactNum::min_positive(p);
        let big_sub = ExactNum::sub_ref(
            &ExactNum::min_positive_normal(p),
            &tiny,
            ps,
            RoundingMode::ToZero,
        );
        assert!(
            !big_sub.inexact() && big_sub.is_subnormal(),
            "largest subnormal p={p}"
        );
        let (ft, fb) = (
            with_full_range(|| exact_mpfr_shifted(&tiny, SHIFT))?,
            with_full_range(|| exact_mpfr_shifted(&big_sub, SHIFT))?,
        );
        for rm in MODES {
            let got = ExactNum::add_ref(&big_sub, &tiny, ps, rm);
            let (want, inexact) = mpfr_reference_shifted(psu, rm, SHIFT, |r, rnd| unsafe {
                mpfr::add(r.as_raw_mut(), fb.as_raw(), ft.as_raw(), rnd)
            })?;
            same(
                &got,
                &want,
                inexact,
                SHIFT,
                &format!("bottom max_sub+min p={p} {rm:?}"),
            )?;
        }
    }
    Ok(())
}

/// `f(x) = g + d` for `x = +-min_positive_normal(p)`, with `g` exactly `x` or `1` and `|d|` near
/// `2^(2 EXPONENT_MIN)`, far below any rounding bit; only the sign of `d` matters. MPFR cannot be the
/// oracle at `EXPONENT_MIN` (its small-argument shortcuts compute `2 EXP(x)`, which wraps at its
/// `emin_min`), so the reference is this expansion. Before the fix, the small-argument test negated
/// `2 e - 2 = i64::MIN` and every one of these returned NaN or a wrong value.
#[test]
fn small_argument_shortcut_at_exponent_min() -> Res {
    type F = fn(&ExactNum, usize, RoundingMode, &mut zenith_float_num::Consts) -> ExactNum;
    // (name, f, g is 1 (else x), sign of d relative to sign(x) (for g = 1: the sign of d itself))
    let fs: [(&str, F, bool, i32); 12] = [
        ("sin", |x, p, rm, cc| x.sin(p, rm, cc), false, -1),
        ("tan", |x, p, rm, cc| x.tan(p, rm, cc), false, 1),
        ("asin", |x, p, rm, cc| x.asin(p, rm, cc), false, 1),
        ("atan", |x, p, rm, cc| x.atan(p, rm, cc), false, -1),
        ("sinh", |x, p, rm, cc| x.sinh(p, rm, cc), false, 1),
        ("tanh", |x, p, rm, cc| x.tanh(p, rm, cc), false, -1),
        ("asinh", |x, p, rm, cc| x.asinh(p, rm, cc), false, -1),
        ("atanh", |x, p, rm, cc| x.atanh(p, rm, cc), false, 1),
        ("expm1", |x, p, rm, cc| x.expm1(p, rm, cc), false, 0),
        ("log1p", |x, p, rm, cc| x.log1p(p, rm, cc), false, 2),
        ("cos", |x, p, rm, cc| x.cos(p, rm, cc), true, -1),
        ("cosh", |x, p, rm, cc| x.cosh(p, rm, cc), true, 1),
    ];
    let mut cc = zenith_float_num::Consts::new()?;
    let mut bad = Vec::new();
    for p in [64usize, 128, 192] {
        for neg in [false, true] {
            let mk = || {
                if neg {
                    ExactNum::min_positive_normal(p).neg()
                } else {
                    ExactNum::min_positive_normal(p)
                }
            };
            let x = mk();
            for &(name, f, to_one, rel) in &fs {
                let g = if to_one { ExactNum::from_u8(1, p) } else { mk() };
                // expm1: d = x^2/2 > 0; log1p: d = -x^2/2 < 0 (whatever the sign of x).
                let d: i32 = match rel {
                    0 => 1,
                    2 => -1,
                    r if to_one => r,
                    r if neg => -r,
                    r => r,
                };
                for rm in MODES {
                    let r = f(&x, p, rm, &mut cc);
                    let toward_d = match rm {
                        RoundingMode::Up => d > 0,
                        RoundingMode::Down => d < 0,
                        RoundingMode::FromZero => (d > 0) != g.is_negative(),
                        RoundingMode::ToZero => (d > 0) == g.is_negative(),
                        RoundingMode::ToEven => false,
                        // Not pinned here (g itself or its neighbour on d's side): the shortcut's
                        // to-odd handling is covered against MPFR inside its range.
                        RoundingMode::ToOdd => {
                            let c = r.compare(&g);
                            if r.is_nan() || c.is_some_and(|s| s != 0 && (s > 0) != (d > 0)) {
                                bad.push(format!("{name} p={p} neg {neg} ToOdd: cmp {c:?}"));
                            }
                            continue;
                        }
                        RoundingMode::None => continue,
                    };
                    let c = r.compare(&g);
                    let ok = !r.is_nan()
                        && match c {
                            Some(0) => !toward_d,
                            Some(s) => toward_d && (s > 0) == (d > 0),
                            None => false,
                        };
                    if !ok {
                        bad.push(format!(
                            "{name}({}{} p={p}) {rm:?}: cmp {c:?} nan {}",
                            if neg { "-" } else { "+" },
                            "min_positive_normal",
                            r.is_nan()
                        ));
                    }
                }
            }
        }
    }
    assert!(bad.is_empty(), "{} wrong:\n{}", bad.len(), bad.join("\n"));
    Ok(())
}
