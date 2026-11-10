//! Harper cosmology-range rows for `Ball` construction and `Enclosure` (2026-10-29): midpoints at zf's EXPONENT_MAX and
//! EXPONENT_MIN, radii beyond zf's range on both sides (Mag holds all of i64), 100_000-bit points at the edges, and
//! `MAX_PRECISION`. Hand-worked oracles: every value is a power of two or a sum of two (value beside each row).
//! A `contains` answer is wrong when it is decided against the exact value; `Unknown` is allowed where noted.

use zenith_float_bound::ball::Ball;
use zenith_float_bound::enclosure::Enclosure;
use zenith_float_bound::error::BoundError;
use zenith_float_bound::mag::Mag;
use zenith_float_bound::tri::Tri;
use zenith_float_bound::MAX_PRECISION;
use zenith_float_num::{ExactNum, EXPONENT_MAX, EXPONENT_MIN};

const TOP: u64 = 0x8000_0000_0000_0000;

fn ok<T>(what: &str, r: Result<T, BoundError>) -> Result<T, String> {
    r.map_err(|e| format!("{what}: unexpected {e:?}"))
}

fn add(a: i64, b: i64) -> Result<i64, String> {
    a.checked_add(b)
        .ok_or_else(|| format!("{a} + {b} overflows"))
}

/// `2^k` exactly.
fn pow2(k: i64) -> Result<ExactNum, String> {
    ExactNum::from_mag_parts(TOP, add(k, 1)?, 64).map_err(|e| format!("pow2({k}): {e:?}"))
}

/// The Mag `2^k`.
fn mag2(k: i64) -> Result<Mag, String> {
    ok("Mag::finite", Mag::finite(TOP, add(k, 1)?))
}

fn ball(mid: ExactNum, rad: Mag, prec: usize) -> Result<Ball, String> {
    ok("Ball::from_mid_rad", Ball::from_mid_rad(mid, rad, prec))
}

fn want(what: &str, b: &Ball, x: &ExactNum, w: Tri) -> Result<(), String> {
    let got = b.contains(x);
    if got == w {
        Ok(())
    } else {
        Err(format!("{what}: contains gave {got:?}, want {w:?}"))
    }
}

/// `Unknown` or the decided answer `w`, never the opposite decision.
fn want_or_unknown(what: &str, b: &Ball, x: &ExactNum, w: Tri) -> Result<(), String> {
    match (b.contains(x), &w) {
        (Tri::Unknown, _) | (Tri::True, Tri::True) | (Tri::False, Tri::False) => Ok(()),
        (got @ (Tri::True | Tri::False), _) => Err(format!(
            "{what}: contains decided {got:?}, exact answer is {w:?}"
        )),
    }
}

/// Point balls at the top binade: `2^(EMAX-1) +- 0` contains itself (True) and not `-2^(EMAX-1)` (|d| = 2^EMAX, past zf).
#[test]
fn cosmo_point_ball_at_emax() -> Result<(), String> {
    let e = EXPONENT_MAX.checked_sub(1).ok_or("emax-1")?;
    for prec in [1, 64, 100_000] {
        let b = ball(pow2(e)?, Mag::zero(), prec)?;
        want("2^(EMAX-1) in itself", &b, &pow2(e)?, Tri::True)?;
        want_or_unknown("-2^(EMAX-1)", &b, &(-pow2(e)?), Tri::False)?;
        want_or_unknown("0", &b, &ExactNum::from(0), Tri::False)?;
        want_or_unknown("2^(EMAX-2)", &b, &pow2(add(e, -1)?)?, Tri::False)?;
    }
    Ok(())
}

/// `2^(EMAX-1) +- 2^EMAX`: the far edge `-2^(EMAX-1)` is at distance exactly 2^EMAX = rad (True or Unknown, never False);
/// `-2^EMAX`... is not an ExactNum, so the outside probe is `-2^(EMAX-1) - 2^(EMAX-2)` at distance 1.25 * 2^EMAX (False or Unknown).
#[test]
fn cosmo_radius_reaching_past_emax() -> Result<(), String> {
    let e = EXPONENT_MAX.checked_sub(1).ok_or("emax-1")?;
    for prec in [1, 64, 100_000] {
        let b = ball(pow2(e)?, mag2(EXPONENT_MAX)?, prec)?;
        want_or_unknown("far edge -2^(EMAX-1)", &b, &(-pow2(e)?), Tri::True)?;
        want("0 inside", &b, &ExactNum::from(0), Tri::True)?;
        let out = (-pow2(e)?).sub_full_prec(&pow2(add(e, -1)?)?);
        want_or_unknown("outside -1.5 * 2^(EMAX-1)", &b, &out, Tri::False)?;
    }
    Ok(())
}

/// `0 +- 2^(EMAX+100)` and `0 +- 2^(2^62 + 2^61)`: a radius far beyond zf's range still holds every finite ExactNum.
#[test]
fn cosmo_radius_beyond_zf_range() -> Result<(), String> {
    let e = EXPONENT_MAX.checked_sub(1).ok_or("emax-1")?;
    for k in [add(EXPONENT_MAX, 100)?, add(1 << 62, 1 << 61)?] {
        for prec in [1, 1000] {
            let b = ball(ExactNum::from(0), mag2(k)?, prec)?;
            want(
                &format!("2^(EMAX-1) in 0 +- 2^{k}"),
                &b,
                &pow2(e)?,
                Tri::True,
            )?;
            want(
                &format!("-2^(EMAX-1) in 0 +- 2^{k}"),
                &b,
                &(-pow2(e)?),
                Tri::True,
            )?;
            want(
                &format!("2^EMIN in 0 +- 2^{k}"),
                &b,
                &pow2(EXPONENT_MIN)?,
                Tri::True,
            )?;
        }
    }
    Ok(())
}

/// The bottom: `2^EMIN +- 0` contains itself, not 0 and not `-2^EMIN`; `0 +- 2^(EMIN-100)` (below zf's range) contains 0
/// and must not claim `2^EMIN` (|x| = 2^EMIN > rad); `0 +- 2^(i64::MIN + 64)` likewise.
#[test]
fn cosmo_edges_at_emin() -> Result<(), String> {
    for prec in [1, 64, 100_000] {
        let b = ball(pow2(EXPONENT_MIN)?, Mag::zero(), prec)?;
        want("2^EMIN in itself", &b, &pow2(EXPONENT_MIN)?, Tri::True)?;
        want_or_unknown("0 not in 2^EMIN", &b, &ExactNum::from(0), Tri::False)?;
        want_or_unknown("-2^EMIN", &b, &(-pow2(EXPONENT_MIN)?), Tri::False)?;
        for k in [add(EXPONENT_MIN, -100)?, add(i64::MIN, 64)?] {
            let t = ball(ExactNum::from(0), mag2(k)?, prec)?;
            want(
                &format!("0 in 0 +- 2^{k}"),
                &t,
                &ExactNum::from(0),
                Tri::True,
            )?;
            want_or_unknown(
                &format!("2^EMIN in 0 +- 2^{k}"),
                &t,
                &pow2(EXPONENT_MIN)?,
                Tri::False,
            )?;
            want_or_unknown(
                &format!("-2^EMIN in 0 +- 2^{k}"),
                &t,
                &(-pow2(EXPONENT_MIN)?),
                Tri::False,
            )?;
        }
    }
    Ok(())
}

/// 100_000-bit points at both edges: `m = 2^E`, `x = 2^E + 2^(E-100000)`. Rad `2^(E-100000)` holds x on its boundary
/// (True); rad `2^(E-100001)` does not (False or Unknown); and `x - 2^(E-100000)`... the mirror `2^E - 2^(E-100000)` too.
#[test]
fn cosmo_hundred_thousand_bit_points_at_edges() -> Result<(), String> {
    let hi = EXPONENT_MAX.checked_sub(2).ok_or("emax-2")?;
    let lo = add(EXPONENT_MIN, 100_002)?;
    for e in [hi, 997, -997, lo] {
        let g = add(e, -100_000)?;
        let x = pow2(e)?.add_full_prec(&pow2(g)?);
        let y = pow2(e)?.sub_full_prec(&pow2(g)?);
        for prec in [64, 100_001] {
            let b = ball(pow2(e)?, mag2(g)?, prec)?;
            want(&format!("E={e} x on the edge"), &b, &x, Tri::True)?;
            want(&format!("E={e} y on the edge"), &b, &y, Tri::True)?;
            let n = ball(pow2(e)?, mag2(add(g, -1)?)?, prec)?;
            want_or_unknown(&format!("E={e} x past half rad"), &n, &x, Tri::False)?;
            want_or_unknown(&format!("E={e} y past half rad"), &n, &y, Tri::False)?;
            let m = ball(x.add_full_prec(&ExactNum::from(0)), Mag::zero(), prec)?;
            want(&format!("E={e} x point in itself"), &m, &x, Tri::True)?;
            want_or_unknown(
                &format!("E={e} 2^E not in x point"),
                &m,
                &pow2(e)?,
                Tri::False,
            )?;
        }
    }
    Ok(())
}

/// Construction and identity at the edges: prec `MAX_PRECISION` is accepted (stored, not rounded), `MAX_PRECISION + 1`
/// is InvalidArgument; `identical` separates equal-valued balls by prec and radius at EMAX/EMIN; a -0 midpoint is +0.
#[test]
fn cosmo_construct_and_identity_at_edges() -> Result<(), String> {
    let e = EXPONENT_MAX.checked_sub(1).ok_or("emax-1")?;
    let b = ok(
        "exact at MAX_PRECISION",
        Ball::exact(pow2(e)?, MAX_PRECISION),
    )?;
    if Enclosure::prec(&b) != MAX_PRECISION || !b.rad().is_zero() {
        return Err(format!("MAX_PRECISION ball {b:?}"));
    }
    let over = MAX_PRECISION.checked_add(1).ok_or("max+1")?;
    match Ball::exact(pow2(e)?, over) {
        Err(BoundError::InvalidArgument { .. }) => {}
        other => {
            return Err(format!(
                "MAX_PRECISION + 1: want InvalidArgument, got {other:?}"
            ))
        }
    }
    for (m, k) in [(e, EXPONENT_MAX), (EXPONENT_MIN, add(EXPONENT_MIN, -100)?)] {
        let a1 = ball(pow2(m)?, mag2(k)?, 64)?;
        let a2 = ball(pow2(m)?, mag2(k)?, 64)?;
        let a3 = ball(pow2(m)?, mag2(k)?, 65)?;
        let a4 = ball(pow2(m)?, mag2(add(k, 1)?)?, 64)?;
        if !a1.identical(&a2) || a1.identical(&a3) || a1.identical(&a4) {
            return Err(format!("identical at 2^{m} +- 2^{k}: {a1:?} {a3:?} {a4:?}"));
        }
    }
    let z = ball(-ExactNum::from(0), mag2(EXPONENT_MIN)?, 64)?;
    if z.mid().is_negative() {
        return Err(format!("-0 midpoint kept negative: {z:?}"));
    }
    want(
        "2^(EMIN-1) in +0 +- 2^EMIN",
        &z,
        &pow2(add(EXPONENT_MIN, -1)?)?,
        Tri::True,
    )
}
