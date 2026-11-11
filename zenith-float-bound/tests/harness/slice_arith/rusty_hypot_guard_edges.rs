//! Rusty (2026-11-10, fold PR #39, edgegate MISSED arith.rs:297:24 x2 and 297:26): the `Ball::hypot` guard bits on the
//! endpoint magnitudes `|m| ± r` (HYPOT_GUARD_BITS = 64; `wp = prec + 64` while that is `<= MAX_PRECISION`, else `prec`).
//! Tightness: the endpoints are `hypot↓` and `hypot↑` of the exact `|m| ± r` (16.S1.8), so for point balls they are the
//! two `prec`-bit neighbours of the exact `hypot`, and the ball is their rounded midpoint with radius one ulp. Without the
//! guard, `|m|` is rounded to `prec` first and the upper endpoint moves one ulp further out. Hand oracle for
//! `x = y = 1 + 2^-60` at prec 53: `hypot = 2^0.5 (1 + 2^-60)`, mantissa (units of 2^-63)
//! `isqrt(2^7 (2^60 + 1)^2) = 13043817825332782223` (not exact), so `lo = 0xb504f333f9de6000`, `hi = lo + 2^11`; the
//! midpoint is a tie, ToEven keeps `lo` (even at 53 bits), radius `hi - lo = 2^-52`. All four sign combinations agree.
//! Top of the range: a valid `prec` is never `InvalidArgument` (only `Resource` here, zf cannot allocate the mantissa);
//! guard dropped at `prec + 64 = MAX_PRECISION + 1`, kept at `MAX_PRECISION - 63`. The exact boundary `prec + 64 ==
//! MAX_PRECISION` and `prec = MAX_PRECISION` are the ignored red (the zf-num InvalidArgument bug, see the powi file).
//! COVERS: Ball::hypot

use zenith_float_bound::ball::Ball;
use zenith_float_bound::enclosure::Enclosure;
use zenith_float_bound::error::BoundError;
use zenith_float_bound::tri::Tri;
use zenith_float_bound::MAX_PRECISION;
use zenith_float_num::ExactNum;

use super::{mag, num, pow2, want_exact, want_invalid, TOP};

#[test]
fn hypot_guard_bits_keep_the_endpoints_tight() -> Result<(), String> {
    let lo = ExactNum::from_mag_parts(0xb504_f333_f9de_6000, 1, 64)
        .map_err(|e| format!("lo: {:?}", e))?;
    let rad = mag(TOP, -51)?; // 2^-52
    for (sx, sy) in [(1i64, 1i64), (-1, 1), (1, -1), (-1, -1)] {
        let x = num(sx).mul_full_prec(&num(1).add_full_prec(&pow2(-60)?));
        let y = num(sy).mul_full_prec(&num(1).add_full_prec(&pow2(-60)?));
        let a = Ball::exact(x, 53).map_err(|e| format!("x: {:?}", e))?;
        let b = Ball::exact(y, 53).map_err(|e| format!("y: {:?}", e))?;
        let h = a.hypot(&b, 53).map_err(|e| format!("hypot: {:?}", e))?;
        want_exact(
            &format!("hypot(({})(1+2^-60), ({})(1+2^-60)) at 53", sx, sy),
            &h,
            &lo,
            &rad,
        )?;
    }
    Ok(())
}

/// `hypot` at a valid `prec`: an enclosure of 5 = hypot(3, 4), or `Resource`; never `InvalidArgument`.
fn valid_prec_row(sx: i64, sy: i64, prec: usize) -> Result<(), String> {
    let what = format!(
        "hypot({}, {}) at prec MAX_PRECISION - {}",
        3 * sx,
        4 * sy,
        MAX_PRECISION.checked_sub(prec).ok_or("prec")?
    );
    let a = Ball::exact(num(3 * sx), prec).map_err(|e| format!("{}: exact: {:?}", what, e))?;
    let b = Ball::exact(num(4 * sy), prec).map_err(|e| format!("{}: exact: {:?}", what, e))?;
    match a.hypot(&b, prec) {
        Ok(r) => match r.contains(&num(5)) {
            Tri::True => Ok(()),
            t @ (Tri::False | Tri::Unknown) => Err(format!(
                "{}: result {:?} does not contain 5 ({:?})",
                what, r, t
            )),
        },
        Err(BoundError::Resource { .. }) => Ok(()),
        // `BoundError` is `#[non_exhaustive]`, so a match from outside the crate needs this arm.
        Err(e) => Err(format!(
            "{}: want an enclosure or Resource for a valid prec, got {:?}",
            what, e
        )),
    }
}

#[test]
fn hypot_guard_dropped_above_max_precision() -> Result<(), String> {
    let prec = MAX_PRECISION.checked_sub(63).ok_or("prec")?; // prec + 64 = MAX_PRECISION + 1
    for (sx, sy) in [(1i64, 1i64), (-1, 1), (1, -1), (-1, -1)] {
        valid_prec_row(sx, sy, prec)?;
    }
    Ok(())
}

#[test]
fn hypot_guard_kept_below_the_boundary() -> Result<(), String> {
    let prec = MAX_PRECISION.checked_sub(127).ok_or("prec")?; // prec + 64 = MAX_PRECISION - 63
    for (sx, sy) in [(1i64, 1i64), (-1, 1), (1, -1), (-1, -1)] {
        valid_prec_row(sx, sy, prec)?;
    }
    Ok(())
}

#[test]
fn hypot_prec_out_of_range_is_invalid() -> Result<(), String> {
    let over = MAX_PRECISION.checked_add(1).ok_or("prec")?;
    let a = Ball::exact(num(3), 64).map_err(|e| format!("exact: {:?}", e))?;
    let b = Ball::exact(num(-4), 64).map_err(|e| format!("exact: {:?}", e))?;
    want_invalid("hypot prec MAX_PRECISION + 1", a.hypot(&b, over))?;
    want_invalid("hypot prec 0", a.hypot(&b, 0))
}

#[test]
#[ignore = "red: zf-num: ExactNum mul/add at a valid prec in (MAX_PRECISION - 64, MAX_PRECISION] that is not a multiple of 64, and at MAX_PRECISION - 64, return InvalidArgument, so Ball::hypot reports InvalidArgument for a valid prec; un-ignore when zf-num reports Resource there"]
fn hypot_valid_prec_at_the_top_is_not_invalid() -> Result<(), String> {
    for d in [64usize, 1, 0] {
        let prec = MAX_PRECISION.checked_sub(d).ok_or("prec")?;
        for (sx, sy) in [(1i64, 1i64), (-1, -1)] {
            valid_prec_row(sx, sy, prec)?;
        }
    }
    Ok(())
}
