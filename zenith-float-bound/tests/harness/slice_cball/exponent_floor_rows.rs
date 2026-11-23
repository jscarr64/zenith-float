//! Slice CB rows at zenith-float's exponent floor and ceiling (Spence PR #40 review, B3 and B4). Arb has no exponent range,
//! so FLINT 3.6.0 cannot give these rows; the oracle is the exact image, worked in the comment beside each row.
//!
//! zf facts the rows rely on: `mul_full_prec` of two in-range values whose product lies 2^60 or more binades below
//! `EXPONENT_MIN` returns zero with the inexact flag (the value is lost), and a sum rounded to `p` bits in the bottom
//! binade keeps only the bits at or above `2^(EXPONENT_MIN − p)` (gradual underflow).

use super::*;
use zenith_float_num::{EXPONENT_MAX, EXPONENT_MIN};

const TIER_PRECS: [usize; 5] = [64, 128, 256, 300, 1024];

fn xball(x: ExactNum, prec: usize) -> Result<Ball, String> {
    ok("exact", Ball::exact(x, prec))
}

fn xc(re: ExactNum, im: ExactNum, prec: usize) -> Result<CBall, String> {
    ok("CBall::new", CBall::new(xball(re, prec)?, xball(im, prec)?))
}

fn neg(x: &ExactNum) -> ExactNum {
    num(0).sub_full_prec(x)
}

/// `[mid − rad, mid + rad]` as exact values (`None` for an all-reals ball).
fn ends(b: &Ball) -> Result<Option<(ExactNum, ExactNum)>, String> {
    if b.rad().is_infinity() {
        return Ok(None);
    }
    let r = if b.rad().is_zero() {
        num(0)
    } else {
        b.rad()
            .to_exact_num()
            .map_err(|e| format!("rad to_exact_num: {:?}", e))?
    };
    Ok(Some((b.mid().sub_full_prec(&r), b.mid().add_full_prec(&r))))
}

fn lt(x: &ExactNum, y: &ExactNum) -> bool {
    x.sub_full_prec(y).is_negative()
}

/// `2^p + 1`: a tie at `p` bits whose even neighbour is `2^p` (nearest-even rounds it toward 0).
fn tie(p: usize) -> Result<ExactNum, String> {
    let p = i64::try_from(p).map_err(|_| "prec")?;
    Ok(pow2(p)?.add_full_prec(&num(1)))
}

/// `2^(−2^61)`: times anything at or below `2^(EXPONENT_MIN + 10)` the product is lost (more than 2^60 binades below).
fn deep() -> Result<ExactNum, String> {
    pow2(-(1i64 << 61))
}

// ---- B3: a product lost below the range must reach the radius ---------------------------------------------------

/// re = a·c − b·d with a·c lost: a = −2^(EMIN+10), c = 2^(−2^61) (`mul_full_prec` gives 0, inexact), b = 2^p + 1, d = 1.
/// The exact re is −(2^p + 1) − δ with δ = 2^(EMIN+10−2^61) > 0. The midpoint sum is the tie −(2^p + 1), which rounds
/// to −2^p with half_ulp 1, so a ball −2^p ± 1 ends at −(2^p + 1) and misses the exact value by δ. The ball must reach
/// strictly past −(2^p + 1). im = a·d + b·c (both products in range) must be contained.
#[test]
fn cball_mul_lost_product_reaches_radius_re() -> Result<(), String> {
    for p in TIER_PRECS {
        let a = neg(&pow2(EXPONENT_MIN + 10)?);
        let z = ok("mul", xc(a, tie(p)?, p)?.mul(&xc(deep()?, num(1), p)?, p))?;
        let Some((lo, _)) = ends(z.re())? else {
            return Err(format!("p={} re: all reals for a finite product", p));
        };
        if !lt(&lo, &neg(&tie(p)?)) {
            return Err(format!(
                "p={} re: ball {:?} ends at −(2^p+1); the exact re −(2^p+1) − 2^(EMIN+10−2^61) is outside",
                p,
                z.re()
            ));
        }
        // im = (2^p+1)·2^(−2^61) − 2^(EMIN+10): the offset is ~2^61 binades below the ball's bits (and too far apart
        // for one exact ExactNum), so the image is inside iff lo < (2^p+1)·2^(−2^61) ≤ hi.
        let x = tie(p)?.mul_full_prec(&deep()?);
        let Some((lo, hi)) = ends(z.im())? else {
            return Err(format!("p={} im: all reals for a finite product", p));
        };
        if !lt(&lo, &x) || lt(&hi, &x) {
            return Err(format!(
                "p={} im: ball {:?} misses (2^p+1)·2^(−2^61) − 2^(EMIN+10)",
                p,
                z.im()
            ));
        }
    }
    Ok(())
}

/// im = a·d + b·c with b·c lost: a = 2^p + 1, d = 1, b = 2^(EMIN+10), c = 2^(−2^61). The exact im is (2^p + 1) + δ; the
/// midpoint tie rounds to 2^p, so a ball 2^p ± 1 misses it. The ball must reach strictly past 2^p + 1.
#[test]
fn cball_mul_lost_product_reaches_radius_im() -> Result<(), String> {
    for p in TIER_PRECS {
        let z = ok(
            "mul",
            xc(tie(p)?, pow2(EXPONENT_MIN + 10)?, p)?.mul(&xc(deep()?, num(1), p)?, p),
        )?;
        let Some((_, hi)) = ends(z.im())? else {
            return Err(format!("p={} im: all reals for a finite product", p));
        };
        if !lt(&tie(p)?, &hi) {
            return Err(format!(
                "p={} im: ball {:?} ends at 2^p+1; the exact im (2^p+1) + 2^(EMIN+10−2^61) is outside",
                p,
                z.im()
            ));
        }
    }
    Ok(())
}

/// sqr cannot hit the B3 tie: re = a² − b² with a² lost adds a positive δ to −b², and an odd square's tie bit always
/// rounds toward 0 (an odd square is 1 mod 8), so the exact value moves into the ball; im = 2ab lost gives a zero
/// midpoint with the bottom ε. Pinned anyway: a = 2^(−0.7·2^62) (a² lost), b = 2^32 + 1 (b² = 2^64 + 2^33 + 1, a tie at
/// p = 64; 2ab stays in range). The exact re −b² + δ and im 2ab must be inside. A guard row: it passes before and after
/// the B3 fix.
#[test]
fn cball_sqr_lost_square_contains_image() -> Result<(), String> {
    let ka = -((1i64 << 62) / 10 * 7);
    let b = pow2(32)?.add_full_prec(&num(1));
    for p in [64usize, 128, 256, 1024] {
        let z = ok("sqr", xc(pow2(ka)?, b.add_full_prec(&num(0)), p)?.sqr(p))?;
        let bb = neg(&b.mul_full_prec(&b));
        // The exact re is −b² + δ with 0 < δ = 2^(2ka) far below any bit of the ball's ends: contained iff lo ≤ −b² < hi.
        let Some((lo, hi)) = ends(z.re())? else {
            return Err(format!("p={} re: all reals", p));
        };
        if lt(&bb, &lo) || !lt(&bb, &hi) {
            return Err(format!(
                "p={} re: ball {:?} misses −b² + 2^(2ka)",
                p,
                z.re()
            ));
        }
        let ab2 = pow2(ka)?.mul_full_prec(&b).mul_full_prec(&num(2));
        want_contains(&format!("p={} im 2ab", p), z.im(), &ab2)?;
    }
    Ok(())
}

// ---- B4: |b|² in acb_inv / the exact-divisor branch of acb_div must not leave the range ---------------------------

/// The ball is finite and contains ±2^(−k)/6, decided exactly as 6·lo ≤ ±2^(−k) ≤ 6·hi.
fn want_sixth(what: &str, b: &Ball, k: i64, negative: bool) -> Result<(), String> {
    let Some((lo, hi)) = ends(b)? else {
        return Err(format!(
            "{}: all reals, but 2^{}/6 is representable",
            what, -k
        ));
    };
    let target = if negative { neg(&pow2(-k)?) } else { pow2(-k)? };
    let (lo6, hi6) = (lo.mul_full_prec(&num(6)), hi.mul_full_prec(&num(6)));
    if lt(&target, &lo6) || lt(&hi6, &target) {
        return Err(format!(
            "{}: ball {:?} misses {}2^{}/6",
            what,
            b,
            if negative { "−" } else { "" },
            -k
        ));
    }
    Ok(())
}

fn ball_pair(m: &ExactNum, rad: Mag, p: usize) -> Result<CBall, String> {
    let r2 = Mag::from_exact_num(&rad.to_exact_num().map_err(|e| format!("{:?}", e))?);
    ok(
        "CBall::new",
        CBall::new(
            ok("ball", Ball::from_mid_rad(m.add_full_prec(&num(0)), rad, p))?,
            ok("ball", Ball::from_mid_rad(m.add_full_prec(&num(0)), r2, p))?,
        ),
    )
}

/// 1/(c + i d) with c = d = 3·2^k: the exact inverse is (1 − i)/(6·2^k). k puts c² + d² = 18·2^(2k) a few binades inside
/// the bottom binade's rounding window at hprec = max(p+3, 30), so the old |b|² kept 3 bits of 18 (16·2^(2k), rounded
/// down): a 1/9 relative error charged as Arb's ulp at prec − 1. For p > 256 the exact-divisor branch of acb_div built
/// c² with `Ball::mul`, whose bottom ε swamps it. The result must contain the image and be finite (the image is about
/// 2^(2^61), representable). The ball divisor (radius 2^-40 relative) must contain the centre's image.
#[test]
fn cball_div_divisor_square_underflow_contains_image() -> Result<(), String> {
    for p in TIER_PRECS {
        let hprec = i64::try_from(p.max(27) + 3).map_err(|_| "prec")?;
        // The exponent of 18·2^(2k) is 2k + 5; put it near EMIN − hprec + 3 (EMIN is odd, so fix the parity).
        let two_k = EXPONENT_MIN - hprec - 2;
        let two_k = if two_k % 2 == 0 { two_k } else { two_k - 1 };
        let k = two_k / 2;
        let c = pow2(k)?.mul_full_prec(&num(3));
        let z = ok(
            "div",
            xc(num(1), num(0), p)?.div(
                &xc(c.add_full_prec(&num(0)), c.add_full_prec(&num(0)), p)?,
                p,
            ),
        )?;
        want_sixth(&format!("p={} exact re", p), z.re(), k, false)?;
        want_sixth(&format!("p={} exact im", p), z.im(), k, true)?;
        let z = ok(
            "div",
            xc(num(1), num(0), p)?.div(&ball_pair(&c, mag(TOP, k + 3 - 40)?, p)?, p),
        )?;
        want_sixth(&format!("p={} ball re", p), z.re(), k, false)?;
        want_sixth(&format!("p={} ball im", p), z.im(), k, true)?;
    }
    Ok(())
}

/// 1/(c + i d) with c = d = 2^k, k = 2^61 + 5: c² overflows zf's range (exponent 2k + 1 > EXPONENT_MAX) but the exact
/// inverse (1 − i)·2^(−k−1) is representable. The result must be finite and contain it, for exact and ball divisors.
#[test]
fn cball_div_divisor_square_overflow_contains_image() -> Result<(), String> {
    let k = (1i64 << 61) + 5;
    if 2 * k < EXPONENT_MAX {
        return Err("row setup: c² must overflow".into());
    }
    for p in TIER_PRECS {
        let want = pow2(-k - 1)?;
        let exact = ok(
            "div",
            xc(num(1), num(0), p)?.div(&xc(pow2(k)?, pow2(k)?, p)?, p),
        )?;
        let ball = ok(
            "div",
            xc(num(1), num(0), p)?.div(&ball_pair(&pow2(k)?, mag(TOP, k + 1 - 40)?, p)?, p),
        )?;
        for (what, z) in [("exact", &exact), ("ball", &ball)] {
            for (part, b, x) in [
                ("re", z.re(), want.add_full_prec(&num(0))),
                ("im", z.im(), neg(&want)),
            ] {
                if !b.is_informative() {
                    return Err(format!(
                        "p={} {} {}: all reals, but ±2^(−k−1) is representable",
                        p, what, part
                    ));
                }
                want_contains(&format!("p={} {} {}", p, what, part), b, &x)?;
            }
        }
    }
    Ok(())
}
