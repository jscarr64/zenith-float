//! Rusty (2026-11-10, fold PR #39, edgegate MISSED arith.rs:145:24): the `Ball::powi` working-precision guard at the top
//! of the precision range. `wp = prec + POWI_GUARD_BITS` only while that is `<= MAX_PRECISION`; above it `wp = prec`.
//! A valid `prec` must never be reported as `InvalidArgument` (CONTRACT: only a `prec` outside
//! `MIN_PRECISION..=MAX_PRECISION` is). At these precisions zenith-float cannot allocate a mantissa, so the honest result is
//! `Resource` (or, were it possible, an enclosing ball). Rows: guard dropped (`prec + 64 = MAX_PRECISION + 1`), guard
//! kept (`prec + 64 = MAX_PRECISION - 63`), out of range on both sides, `n` = 0 to 3, both signs of the midpoint.
//! The guard-kept boundary `prec + 64 == MAX_PRECISION` and `prec = MAX_PRECISION` itself are the ignored red below:
//! zenith-float-num reports `InvalidArgument` there (a zf-num bug, not this guard).
//! COVERS: Ball::powi

use zenith_float_bound::ball::Ball;
use zenith_float_bound::enclosure::Enclosure;
use zenith_float_bound::error::BoundError;
use zenith_float_bound::tri::Tri;
use zenith_float_bound::MAX_PRECISION;

use super::{num, want_invalid};

/// 3^n with the sign of `v` for odd n.
fn power(v: i64, n: u32) -> Result<i64, String> {
    v.checked_pow(n).ok_or_else(|| format!("{}^{}", v, n))
}

/// `powi` at a valid `prec`: an enclosing ball, or `Resource`; never `InvalidArgument` or another error.
fn valid_prec_row(v: i64, n: u32, prec: usize) -> Result<(), String> {
    let what = format!(
        "({})^{} at prec MAX_PRECISION - {}",
        v,
        n,
        MAX_PRECISION.checked_sub(prec).ok_or("prec")?
    );
    let b = Ball::exact(num(v), prec).map_err(|e| format!("{}: exact: {:?}", what, e))?;
    match b.powi(u64::from(n), prec) {
        Ok(r) => match r.contains(&num(power(v, n)?)) {
            Tri::True => Ok(()),
            t @ (Tri::False | Tri::Unknown) => Err(format!(
                "{}: result {:?} does not contain the power ({:?})",
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
fn powi_guard_dropped_above_max_precision() -> Result<(), String> {
    // prec + 64 = MAX_PRECISION + 1: the guard must be dropped (a word-multiple prec, which zf-num accepts).
    let prec = MAX_PRECISION.checked_sub(63).ok_or("prec")?;
    for v in [3i64, -3] {
        for n in 0..=3u32 {
            valid_prec_row(v, n, prec)?;
        }
    }
    Ok(())
}

#[test]
fn powi_guard_kept_below_the_boundary() -> Result<(), String> {
    // prec + 64 = MAX_PRECISION - 63 <= MAX_PRECISION: the guard is kept.
    let prec = MAX_PRECISION.checked_sub(127).ok_or("prec")?;
    for v in [3i64, -3] {
        for n in 0..=3u32 {
            valid_prec_row(v, n, prec)?;
        }
    }
    Ok(())
}

#[test]
fn powi_prec_out_of_range_is_invalid() -> Result<(), String> {
    let over = MAX_PRECISION.checked_add(1).ok_or("prec")?;
    for v in [3i64, -3] {
        let b = Ball::exact(num(v), 64).map_err(|e| format!("exact: {:?}", e))?;
        for n in [0u64, 1, 2] {
            want_invalid("powi prec MAX_PRECISION + 1", b.powi(n, over))?;
            want_invalid("powi prec 0", b.powi(n, 0))?;
        }
    }
    Ok(())
}

#[test]
#[ignore = "red: zf-num: ExactNum mul/add at a valid prec in (MAX_PRECISION - 64, MAX_PRECISION] that is not a multiple of 64, and at MAX_PRECISION - 64, return InvalidArgument (Ball::mul maps it), so Ball::powi reports InvalidArgument for a valid prec; un-ignore when zf-num reports Resource there"]
fn powi_valid_prec_at_the_top_is_not_invalid() -> Result<(), String> {
    // prec + 64 == MAX_PRECISION (guard kept at the exact boundary), prec = MAX_PRECISION - 1 and MAX_PRECISION.
    for d in [64usize, 1, 0] {
        let prec = MAX_PRECISION.checked_sub(d).ok_or("prec")?;
        for v in [3i64, -3] {
            for n in 1..=2u32 {
                valid_prec_row(v, n, prec)?;
            }
        }
    }
    Ok(())
}
