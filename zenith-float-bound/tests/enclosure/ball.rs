//! PLAN C5 (2026-10-23): sealed `Enclosure` with `prec`, `identical -> bool`, `contains -> Tri`;
//! CONTRACT 6.5 b (identical includes prec), O20. Oracle: hand-built balls of small dyadics.
//! COVERS: Enclosure for Ball, Sealed for Ball

use zenith_bound::ball::Ball;
use zenith_bound::enclosure::Enclosure;
use zenith_bound::mag::Mag;
use zenith_bound::tri::Tri;
use zenith_float_num::{ExactNum, INF_POS, NAN};

const TOP: u64 = 1u64 << 63;

fn n(v: i64) -> ExactNum {
    ExactNum::from(v)
}

/// `v ± 2^k` at `prec`.
fn ball(v: i64, k: i64, prec: usize) -> Result<Ball, String> {
    let rad = Mag::finite(TOP, k.checked_add(1).ok_or("k")?).map_err(|e| format!("{:?}", e))?;
    Ball::from_mid_rad(n(v), rad, prec).map_err(|e| format!("{:?}", e))
}

fn point(v: i64, prec: usize) -> Result<Ball, String> {
    Ball::exact(n(v), prec).map_err(|e| format!("{:?}", e))
}

fn generic_prec<E: Enclosure>(e: &E) -> usize {
    e.prec()
}

#[test]
fn prec_through_the_trait() -> Result<(), String> {
    if generic_prec(&point(3, 77)?) != 77 {
        return Err("prec".into());
    }
    Ok(())
}

#[test]
fn identical_compares_mid_rad_and_prec() -> Result<(), String> {
    let a = ball(5, 0, 64)?;
    let cases = [
        (ball(5, 0, 64)?, true),
        (ball(5, 0, 65)?, false),
        (ball(5, 1, 64)?, false),
        (ball(6, 0, 64)?, false),
        (point(5, 64)?, false),
    ];
    for (i, (b, want)) in cases.iter().enumerate() {
        if a.identical(b) != *want || b.identical(&a) != *want {
            return Err(format!("case {}: identical want {}", i, want));
        }
    }
    let z = Ball::exact(-n(0), 64).map_err(|e| format!("{:?}", e))?;
    if !z.identical(&point(0, 64)?) {
        return Err("-0 and +0 points must be identical after normalisation".into());
    }
    if !Ball::all_reals(9)
        .map_err(|e| format!("{:?}", e))?
        .identical(&Ball::all_reals(9).map_err(|e| format!("{:?}", e))?)
    {
        return Err("all-reals".into());
    }
    Ok(())
}

#[test]
fn contains_decides_exact_cases() -> Result<(), String> {
    // 5 ± 2 = [3, 7]
    let b = ball(5, 1, 64)?;
    let cases: [(i64, Tri); 7] = [
        (5, Tri::True),
        (3, Tri::True),
        (7, Tri::True),
        (2, Tri::False),
        (8, Tri::False),
        (-7, Tri::False),
        (1 << 40, Tri::False),
    ];
    for (x, want) in cases {
        let got = b.contains(&n(x));
        if got != want {
            return Err(format!(
                "5±2 contains {}: want {:?}, got {:?}",
                x, want, got
            ));
        }
    }
    Ok(())
}

#[test]
fn contains_on_points_and_all_reals() -> Result<(), String> {
    let p = point(-9, 64)?;
    if p.contains(&n(-9)) != Tri::True || p.contains(&n(9)) != Tri::False {
        return Err("point ball".into());
    }
    let wide = Ball::exact(ExactNum::from(u64::MAX), 8).map_err(|e| format!("{:?}", e))?;
    if wide.contains(&ExactNum::from(u64::MAX)) != Tri::True {
        return Err("exact wide point must contain itself".into());
    }
    let all = Ball::all_reals(64).map_err(|e| format!("{:?}", e))?;
    if all.contains(&n(i64::MIN)) != Tri::True {
        return Err("all-reals contains every real".into());
    }
    if all.contains(&NAN) != Tri::False || p.contains(&INF_POS) != Tri::False {
        return Err("NaN/Inf are never contained".into());
    }
    Ok(())
}

/// `sign * 2^e` as an exact zf value (1 bit).
fn pow2(e: i64, negative: bool) -> ExactNum {
    let one = if negative { n(-1) } else { n(1) };
    one.ldexp(e, 64, zenith_float_num::RoundingMode::Down)
}

#[test]
fn contains_is_sound_when_the_difference_is_inexact() -> Result<(), String> {
    // 1 ± 1 against x = ∓2^-300: |mid - x| = 1 + 2^-300 > 1, which needs ~300 bits, so `contains` must not say True.
    // L0b mutants 132-C-hiadd, 134-C-loadd, 136-C-upmax (a bracket end rounded the wrong way or the min taken).
    let one_one = ball(1, 0, 64)?;
    if one_one.contains(&pow2(-300, true)) == Tri::True {
        return Err("1±1 must not contain -2^-300 (distance 1 + 2^-300)".into());
    }
    let m_one = ball(-1, 0, 64)?;
    if m_one.contains(&pow2(-300, false)) == Tri::True {
        return Err("-1±1 must not contain 2^-300 (distance 1 + 2^-300)".into());
    }
    // and the inside cases stay decided: distance 1 - 2^-300 <= 1
    if one_one.contains(&pow2(-300, false)) != Tri::True
        || m_one.contains(&pow2(-300, true)) != Tri::True
    {
        return Err("1±1 contains 2^-300 and -1±1 contains -2^-300".into());
    }
    Ok(())
}

#[test]
fn contains_works_at_the_wider_operand_precision() -> Result<(), String> {
    // mid = 2^200 (exact at 64 bits), x = 2^200 - 1 + 2^-300 at 600 bits, rad 1: distance 1 - 2^-300 <= 1, so True.
    // At the narrower operand's precision the first step already rounds by 2^8 (L0b mutant 143-C-wpmin).
    let p = 600;
    let rm = zenith_float_num::RoundingMode::Down;
    let x = pow2(200, false)
        .sub(&n(1), p, rm)
        .add(&pow2(-300, false), p, rm);
    let mid = pow2(200, false);
    let b = Ball::from_mid_rad(
        mid,
        Mag::finite(TOP, 1).map_err(|e| format!("{:?}", e))?,
        64,
    )
    .map_err(|e| format!("{:?}", e))?;
    match b.contains(&x) {
        Tri::True => Ok(()),
        t @ (Tri::False | Tri::Unknown) => Err(format!(
            "2^200 ± 1 contains 2^200 - 1 + 2^-300: want True, got {:?}",
            t
        )),
    }
}

#[test]
fn contains_is_tight_at_both_edges() -> Result<(), String> {
    // Paul's ruling 2026-10-23: contains is tight. |mid - x| = 1 + 2^-300 > 1 = rad is decided, so the answer is False,
    // at the lower edge (1 ± 1 vs -2^-300, mirror -1 ± 1 vs 2^-300) and the upper edge (1 ± 1 vs 2 + 2^-300, mirror
    // -1 ± 1 vs -2 - 2^-300). Formerly Unknown (L0b survivors 138-C-lowgt, 141-C-lozero, 142-C-lohi).
    let rm = zenith_float_num::RoundingMode::Down;
    let up_edge = n(2).add(&pow2(-300, false), 600, rm);
    let down_edge = n(-2).sub(&pow2(-300, false), 600, rm);
    let cases = [
        ("1±1 vs -2^-300", ball(1, 0, 64)?, pow2(-300, true)),
        ("-1±1 vs 2^-300", ball(-1, 0, 64)?, pow2(-300, false)),
        ("1±1 vs 2 + 2^-300", ball(1, 0, 64)?, up_edge),
        ("-1±1 vs -2 - 2^-300", ball(-1, 0, 64)?, down_edge),
    ];
    for (what, b, x) in cases {
        let t = b.contains(&x);
        if t != Tri::False {
            return Err(format!("{}: want False, got {:?}", what, t));
        }
    }
    Ok(())
}

#[test]
fn contains_stays_unknown_when_rad_is_past_zenith_float_range() -> Result<(), String> {
    // d = max - (-max) = 2 max overflows to +Inf in zenith-float, and rad = 2 max is no zf value either, so the
    // comparison cannot be decided there; the true answer is True, so the one wrong answer is False.
    let m = Mag::from_exact_num(&ExactNum::max_value(64));
    let exp = m.exp().ok_or("exp")?.checked_add(1).ok_or("exp+1")?;
    let rad = Mag::finite(m.mant().ok_or("mant")?, exp).map_err(|e| format!("{:?}", e))?;
    let b = Ball::from_mid_rad(ExactNum::max_value(64), rad, 64).map_err(|e| format!("{:?}", e))?;
    let t = b.contains(&-ExactNum::max_value(64));
    if t == Tri::False {
        return Err("max ± 2 max vs -max: False is wrong".into());
    }
    Ok(())
}

#[test]
fn infinite_points_are_outside_even_a_rad_past_zenith_float_range() -> Result<(), String> {
    // Kills 127-C-naninf under the tight rule: with an in-range rad, +-Inf already falls out as False, so the
    // explicit is_inf test shows only when rad is no zf value (the overflowing bracket alone would say Unknown).
    let m = Mag::from_exact_num(&ExactNum::max_value(64));
    let exp = m.exp().ok_or("exp")?.checked_add(1).ok_or("exp+1")?;
    let rad = Mag::finite(m.mant().ok_or("mant")?, exp).map_err(|e| format!("{:?}", e))?;
    let b = Ball::from_mid_rad(n(0), rad, 64).map_err(|e| format!("{:?}", e))?;
    for x in [INF_POS, -INF_POS] {
        let t = b.contains(&x);
        if t != Tri::False {
            return Err(format!("0 ± 2 max vs {:?}: want False, got {:?}", x, t));
        }
    }
    Ok(())
}
