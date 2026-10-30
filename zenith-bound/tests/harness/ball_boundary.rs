//! CONTRACT 3.2 (round-up underflow to the smallest positive `Mag`), 6.4 and 15.3 (`is_informative`), 14.7 (ulp terms);
//! TEST_PLAN L0b-14 (degenerate and all-reals) at the zf 2.0 exponent bounds `ZF_EXPONENT_MIN`, `MIN + 1`, `MAX - 1`, `ZF_EXPONENT_MAX`
//! (Rusty: `Exponent = i64` with range +-(2^62 - 1), so an exponent sum or difference of two in-range exponents fits in i64). The bounds
//! are defined once in `common::generators`. Spence's rule (2026-10-23 7:51 PM CT): an ulp exponent below what `Mag` holds rounds
//! outward to the smallest positive `Mag`, never to zero. Inside the zf range that underflow happens only at precisions above 2^62
//! (`ulp_mag(ZF_EXPONENT_MIN, p)` first clamps at p = 2^62 + 3, `half_ulp` at 2^62 + 2), so those precisions are listed explicitly.
//! Oracle: the exponent field `mid_exp - prec + 1` (full ulp) or `mid_exp - prec` (half ulp) computed in i128 and clamped up to
//! `i64::MIN`; the mantissa is always the top bit.
//! COVERS: Mag::ulp_mag, Mag::half_ulp, Mag::smaller_than_half_ulp, Mag::mul_up, Mag::div_up, Ball::from_mid_rad, Ball::is_informative,
//! Ball::is_exact

use crate::common::generators::{TOP, ZF_EXPONENT_MAX, ZF_EXPONENT_MIN};
use zenith_bound::ball::Ball;
use zenith_bound::mag::Mag;
use zenith_float_num::ExactNum;

fn bound(base: i64, step: i64) -> Result<i64, String> {
    base.checked_add(step)
        .ok_or_else(|| format!("{} + {} overflows i64", base, step))
}

/// The four zf exponent bounds: MIN, MIN + 1, MAX - 1, MAX.
fn bounds() -> Result<[i64; 4], String> {
    Ok([
        ZF_EXPONENT_MIN,
        bound(ZF_EXPONENT_MIN, 1)?,
        bound(ZF_EXPONENT_MAX, -1)?,
        ZF_EXPONENT_MAX,
    ])
}

/// 2^62 + k as a precision. zf and zb are 64-bit only (P6), so this is always a usize.
fn above_2_62(k: usize) -> Result<usize, String> {
    (1usize << 62).checked_add(k).ok_or_else(|| String::from("2^62 + k overflows usize"))
}

/// One precision list (P6, 64-bit only; ZB-TD03). MAX_PRECISION - 1, MAX_PRECISION = u32::MAX and MAX_PRECISION + 1 bracket
/// the ZF-M32 limit 2^32 - 1; 2^62 + k reaches the ulp underflow inside the zf range; usize::MAX is the far end.
fn precisions() -> Result<Vec<usize>, String> {
    let mut ps = vec![0, 1, 2, 53, 4096, (1 << 32) - 2, u32::MAX as usize, 1 << 32];
    for k in 1..=3 {
        ps.push(above_2_62(k)?);
    }
    ps.push(usize::MAX);
    Ok(ps)
}

fn pi(p: usize) -> Result<i128, String> {
    i128::try_from(p).map_err(|_| format!("prec {} as i128", p))
}

/// The power of two with exponent field `field`, clamped up to `i64::MIN` (outward rounding).
fn expect_pow2(what: &str, got: &Mag, field: i128) -> Result<(), String> {
    let want = i64::try_from(field.max(i128::from(i64::MIN)))
        .map_err(|_| format!("{}: field {} above i64::MAX", what, field))?;
    if got.exp() == Some(want) && got.mant() == Some(TOP) {
        Ok(())
    } else {
        Err(format!(
            "{}: want 2^({} - 1) (field {}), got {:?}",
            what, want, field, got
        ))
    }
}

#[test]
fn ulp_mag_at_the_zf_bounds_rounds_an_underflow_out_to_the_smallest_positive_mag(
) -> Result<(), String> {
    for e in bounds()? {
        for p in precisions()? {
            let field = i128::from(e) - pi(p)? + 1;
            expect_pow2(
                &format!("ulp_mag({}, {})", e, p),
                &Mag::ulp_mag(e, p),
                field,
            )?;
        }
    }
    Ok(())
}

#[test]
fn half_ulp_at_the_zf_bounds_rounds_an_underflow_out_to_the_smallest_positive_mag(
) -> Result<(), String> {
    for e in bounds()? {
        for p in precisions()? {
            let field = i128::from(e) - pi(p)?;
            expect_pow2(
                &format!("half_ulp({}, {})", e, p),
                &Mag::half_ulp(e, p),
                field,
            )?;
        }
    }
    Ok(())
}

/// A second route to the same outward clamp: a pre-round exponent of `i64::MIN + 1` at small precisions. The ulp field
/// `e - p + 1` is `i64::MIN` at p = 2 and clamps from p = 3; the half field `e - p` is `i64::MIN` at p = 1 and clamps from p = 2.
#[test]
fn the_ulp_underflow_from_i64_min_plus_1_clamps_to_the_smallest_positive_mag() -> Result<(), String> {
    let min_pos = Mag::finite(TOP, i64::MIN).map_err(|err| format!("min_pos: {:?}", err))?;
    let e = bound(i64::MIN, 1)?;
    let cases = [
        ("ulp_mag p=1", Mag::ulp_mag(e, 1), bound(i64::MIN, 1)?),
        ("ulp_mag p=2", Mag::ulp_mag(e, 2), i64::MIN),
        ("ulp_mag p=3", Mag::ulp_mag(e, 3), i64::MIN),
        ("ulp_mag p=4", Mag::ulp_mag(e, 4), i64::MIN),
        ("half_ulp p=1", Mag::half_ulp(e, 1), i64::MIN),
        ("half_ulp p=2", Mag::half_ulp(e, 2), i64::MIN),
        ("half_ulp p=3", Mag::half_ulp(e, 3), i64::MIN),
    ];
    for (what, got, want) in cases {
        if got.exp() != Some(want) || got.mant() != Some(TOP) {
            return Err(format!(
                "{} at i64::MIN + 1: want field {}, got {:?}",
                what, want, got
            ));
        }
    }
    let wide = usize::try_from(u32::MAX).map_err(|_| String::from("u32::MAX as usize"))?;
    for (what, got) in [
        ("ulp_mag", Mag::ulp_mag(ZF_EXPONENT_MIN, wide)),
        ("half_ulp", Mag::half_ulp(ZF_EXPONENT_MIN, wide)),
        ("ulp_mag i64::MIN", Mag::ulp_mag(i64::MIN, wide)),
    ] {
        let field = if what == "ulp_mag i64::MIN" {
            i128::from(i64::MIN) - pi(wide)? + 1
        } else if what == "ulp_mag" {
            i128::from(ZF_EXPONENT_MIN) - pi(wide)? + 1
        } else {
            i128::from(ZF_EXPONENT_MIN) - pi(wide)?
        };
        expect_pow2(&format!("{}(.., u32::MAX)", what), &got, field)?;
    }
    if Mag::ulp_mag(i64::MIN, wide) != min_pos {
        return Err(String::from(
            "ulp_mag(i64::MIN, u32::MAX) is not the smallest positive Mag",
        ));
    }
    Ok(())
}

#[test]
fn the_ulp_underflow_inside_the_zf_range_starts_just_above_precision_2_62() -> Result<(), String> {
    // At ZF_EXPONENT_MIN = -(2^62 - 1): ulp field -(2^62 - 1) - p + 1 reaches i64::MIN at p = 2^62 + 2 (still exact) and clamps from
    // p = 2^62 + 3; the half ulp field reaches i64::MIN at p = 2^62 + 1 and clamps from p = 2^62 + 2.
    let min_pos = Mag::finite(TOP, i64::MIN).map_err(|err| format!("min_pos: {:?}", err))?;
    let cases = [
        ("ulp_mag", Mag::ulp_mag(ZF_EXPONENT_MIN, above_2_62(2)?), i64::MIN),
        ("ulp_mag", Mag::ulp_mag(ZF_EXPONENT_MIN, above_2_62(3)?), i64::MIN),
        ("half_ulp", Mag::half_ulp(ZF_EXPONENT_MIN, above_2_62(1)?), i64::MIN),
        ("half_ulp", Mag::half_ulp(ZF_EXPONENT_MIN, above_2_62(2)?), i64::MIN),
        ("ulp_mag", Mag::ulp_mag(ZF_EXPONENT_MIN, above_2_62(1)?), bound(i64::MIN, 1)?),
    ];
    for (what, got, want) in cases {
        if got.exp() != Some(want) || got.mant() != Some(TOP) {
            return Err(format!(
                "{} at ZF_EXPONENT_MIN: want field {}, got {:?}",
                what, want, got
            ));
        }
    }
    if Mag::ulp_mag(ZF_EXPONENT_MIN, usize::MAX) != min_pos {
        return Err(String::from(
            "ulp_mag(ZF_EXPONENT_MIN, usize::MAX) is not the smallest positive Mag",
        ));
    }
    Ok(())
}

/// ZF-M32: u32::MAX is exactly MAX_PRECISION (2^32 - 1), so Ball::exact accepts it and keeps the precision.
#[test]
fn ball_exact_accepts_u32_max_precision() -> Result<(), String> {
    let wide = usize::try_from(u32::MAX).map_err(|_| String::from("u32::MAX as usize"))?;
    match zenith_bound::ball::Ball::exact(ExactNum::from_word(1, 64), wide) {
        Ok(b) if b.prec() == wide => Ok(()),
        other => Err(format!("Ball::exact(1, u32::MAX): want Ok at prec u32::MAX, got {:?}", other)),
    }
}

/// ZF-M32: MAX_PRECISION + 1 = 2^32 is refused with InvalidArgument.
#[test]
fn ball_exact_refuses_max_precision_plus_1() -> Result<(), String> {
    let over = 1usize << 32;
    match zenith_bound::ball::Ball::exact(ExactNum::from_word(1, 64), over) {
        Err(zenith_bound::error::BoundError::InvalidArgument { .. }) => Ok(()),
        other => Err(format!("Ball::exact(1, 2^32): want InvalidArgument, got {:?}", other)),
    }
}

#[test]
fn the_clamped_radius_is_not_below_half_an_ulp() -> Result<(), String> {
    // smaller_than_half_ulp compares against the unclamped threshold, so a radius rounded up to the smallest positive Mag is never
    // reported as below half an ulp (that would let a caller drop it).
    for e in bounds()? {
        for p in precisions()? {
            if Mag::half_ulp(e, p).smaller_than_half_ulp(e, p) {
                return Err(format!(
                    "half_ulp({}, {}) reported as smaller than half an ulp",
                    e, p
                ));
            }
        }
    }
    Ok(())
}

#[test]
fn balls_with_an_ulp_radius_at_the_zf_bounds_are_informative_and_not_exact() -> Result<(), String> {
    // Inside +-(2^62 - 1) no ulp overflows to Infinity (the largest field is ZF_EXPONENT_MAX + 1 = 2^62), so every ulp radius at every
    // bound and precision, clamped or not, keeps the ball informative.
    for e in bounds()? {
        for p in precisions()? {
            for rad in [Mag::ulp_mag(e, p), Mag::half_ulp(e, p)] {
                let b = Ball::from_mid_rad(ExactNum::from(1i64), rad, 53)
                    .map_err(|err| format!("from_mid_rad at {}, {}: {:?}", e, p, err))?;
                if !b.is_informative() || b.is_exact() {
                    return Err(format!(
                        "ulp radius at mid exponent {}, prec {}: is_informative {} is_exact {}",
                        e,
                        p,
                        b.is_informative(),
                        b.is_exact()
                    ));
                }
            }
        }
    }
    Ok(())
}

#[test]
fn an_exponent_sum_and_difference_at_both_zf_edges_do_not_overflow() -> Result<(), String> {
    // The point of +-(2^62 - 1): a + b and a - b of in-range exponents fit in i64, and a Mag product or quotient of powers of two at the
    // edges is the exact power (field a + b - 1, a - b + 1), neither Infinity nor clamped.
    let pairs = [
        (ZF_EXPONENT_MAX, ZF_EXPONENT_MAX),
        (ZF_EXPONENT_MIN, ZF_EXPONENT_MIN),
        (ZF_EXPONENT_MAX, ZF_EXPONENT_MIN),
        (ZF_EXPONENT_MIN, ZF_EXPONENT_MAX),
    ];
    for (a, b) in pairs {
        let sum = a
            .checked_add(b)
            .ok_or_else(|| format!("exponent sum {} + {} overflows i64", a, b))?;
        let diff = a
            .checked_sub(b)
            .ok_or_else(|| format!("exponent difference {} - {} overflows i64", a, b))?;
        let x = Mag::finite(TOP, a).map_err(|err| format!("Mag 2^({} - 1): {:?}", a, err))?;
        let y = Mag::finite(TOP, b).map_err(|err| format!("Mag 2^({} - 1): {:?}", b, err))?;
        expect_pow2(
            &format!("mul_up at {}, {}", a, b),
            &x.mul_up(&y),
            i128::from(sum) - 1,
        )?;
        expect_pow2(
            &format!("div_up at {}, {}", a, b),
            &x.div_up(&y),
            i128::from(diff) + 1,
        )?;
        let field_sum = i128::from(sum) - 1;
        let field_diff = i128::from(diff) + 1;
        if field_sum < i128::from(i64::MIN) || field_diff > i128::from(i64::MAX) {
            return Err(format!("Mag field out of i64 at {}, {}", a, b));
        }
    }
    Ok(())
}
