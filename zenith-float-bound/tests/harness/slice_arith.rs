//! Slice S1 harness rows (`src/arith.rs`). Owned by slice S1 alone.
//! Goldens: tests/goldens/slices/arith/ (its own SHA256SUMS and GENERATORS). Mutants: tools/mutants/manifest.d/arith.txt, ids 200..=299.
//! Golden file headers for this slice are declared here, not in tests/common/goldens.rs.
//! CONTRACT 2.2, 3.2, 6.3, 9.3, 9.4, 10.2b, 10.2e, 10.4, 15.1 and 16.S1.1 to 16.S1.8. Oracles: hand-verified exact rows
//! (each value worked by hand in the comment beside it); containment of exact images computed with zenith-float's exact
//! `add_full_prec`/`mul_full_prec` (and directed `sqrt` brackets) at sample points of the input balls; the committed L0b
//! enclosure rows tests/goldens/ball_enclosure_rows.txt (exact images rounded outward by scripts/gen_ball_enclosure.py).
//! COVERS: Ball::add, Ball::sub, Ball::mul, Ball::mul_exact, Ball::sqr, Ball::powi, Ball::sqrt, Neg for Ball, Ball::div, Ball::inv, Ball::hypot

use zenith_bound::ball::Ball;
use zenith_bound::enclosure::Enclosure;
use zenith_bound::error::BoundError;
use zenith_bound::mag::Mag;
use zenith_bound::policy::DomainPolicy;
use zenith_bound::tri::Tri;
use zenith_float_num::{ExactNum, RoundingMode, EXPONENT_MAX, EXPONENT_MIN};

use crate::common::goldens::golden_dir;
use crate::common::rng::SplitMix64;

const TOP: u64 = 0x8000_0000_0000_0000;
const PRECS: [usize; 9] = [1, 2, 53, 64, 65, 113, 1024, 4096, 16384];
/// Mid exponents: 0, small, 1e±300 (2^±997), and far beyond (to 2^±2^28, so a product stays inside zf's i32 range).
const EXPS: [i64; 13] = [
    0,
    10,
    -10,
    997,
    -997,
    1000,
    -1000,
    100_000,
    -100_000,
    1 << 20,
    -(1 << 20),
    1 << 28,
    -(1 << 28),
];

fn num(v: i64) -> ExactNum {
    ExactNum::from(v)
}

fn mag(mant: u64, exp: i64) -> Result<Mag, String> {
    Mag::finite(mant, exp).map_err(|e| format!("Mag::finite({:#x}, {}): {:?}", mant, exp, e))
}

/// `2^k` exactly.
fn pow2(k: i64) -> Result<ExactNum, String> {
    let e = k.checked_add(1).ok_or("pow2 exponent")?;
    ExactNum::from_mag_parts(TOP, e, 64).map_err(|e| format!("pow2({}): {:?}", k, e))
}

/// `v * 2^k` exactly.
fn dy(v: i64, k: i64) -> Result<ExactNum, String> {
    Ok(num(v).mul_full_prec(&pow2(k)?))
}

fn ok(what: &str, r: Result<Ball, BoundError>) -> Result<Ball, String> {
    r.map_err(|e| format!("{}: unexpected {:?}", what, e))
}

fn ball(mid: ExactNum, rad: Mag, prec: usize) -> Result<Ball, String> {
    ok("from_mid_rad", Ball::from_mid_rad(mid, rad, prec))
}

fn pt(v: i64, prec: usize) -> Result<Ball, String> {
    ok("exact", Ball::exact(num(v), prec))
}

fn want_exact(what: &str, b: &Ball, mid: &ExactNum, rad: &Mag) -> Result<(), String> {
    if b.mid() != mid || b.rad() != rad {
        return Err(format!(
            "{}: got {:?}, want mid {:?} rad {:?}",
            what, b, mid, rad
        ));
    }
    Ok(())
}

fn want_domain(what: &str, r: Result<Ball, BoundError>) -> Result<(), String> {
    match r {
        Err(BoundError::Domain { .. }) => Ok(()),
        other => Err(format!("{}: want Domain, got {:?}", what, other)),
    }
}

fn want_contains(what: &str, b: &Ball, x: &ExactNum) -> Result<(), String> {
    match b.contains(x) {
        Tri::True => Ok(()),
        t @ (Tri::False | Tri::Unknown) => Err(format!(
            "{}: {:?} does not contain {:?} ({:?})",
            what, b, x, t
        )),
    }
}

// ---- hand-verified exact rows ------------------------------------------------------------------------------------

#[test]
fn add_exact_points() -> Result<(), String> {
    // 3 + 5 = 8, both exact: rad 0 (CONTRACT 2.2).
    let b = ok("add", pt(3, 53)?.add(&pt(5, 53)?, 53))?;
    want_exact("3+5", &b, &num(8), &Mag::zero())?;
    // mixed signs: -3 + 5 = 2.
    let b = ok("add", pt(-3, 53)?.add(&pt(5, 53)?, 53))?;
    want_exact("-3+5", &b, &num(2), &Mag::zero())
}

#[test]
fn add_radii_sum_exactly() -> Result<(), String> {
    // (1 ± 2^-10) + (2 ± 2^-10) = 3 ± 2^-9; Mag 2^-10 is finite(TOP, -9), 2^-9 is finite(TOP, -8).
    let a = ball(num(1), mag(TOP, -9)?, 53)?;
    let b = ball(num(2), mag(TOP, -9)?, 53)?;
    want_exact("radii", &ok("add", a.add(&b, 53))?, &num(3), &mag(TOP, -8)?)
}

#[test]
fn add_rounding_adds_half_ulp() -> Result<(), String> {
    // 1 + 2^-60 at 53 bits rounds to 1 (exponent 1); ε = half_ulp(1, 53) = 2^-53 = finite(TOP, -52).
    let b = ok(
        "add",
        pt(1, 53)?.add(&ok("exact", Ball::exact(pow2(-60)?, 53))?, 53),
    )?;
    want_exact("1+2^-60", &b, &num(1), &mag(TOP, -52)?)?;
    want_contains("1+2^-60", &b, &num(1).add_full_prec(&pow2(-60)?))
}

#[test]
fn add_at_prec_one() -> Result<(), String> {
    // 1 + 2 = 3 at 1 bit: ties-to-even of 3 = 0.11b·2^2 gives 4 (0.1b·2^3); exponent 3, ε = half_ulp(3, 1) = 2^1 = finite(TOP, 2).
    let b = ok("add", pt(1, 1)?.add(&pt(2, 1)?, 1))?;
    want_exact("1+2 @1", &b, &num(4), &mag(TOP, 2)?)?;
    want_contains("1+2 @1", &b, &num(3))
}

#[test]
fn sub_rows() -> Result<(), String> {
    // (5 ± 2^-10) - (3 ± 2^-10) = 2 ± 2^-9.
    let a = ball(num(5), mag(TOP, -9)?, 53)?;
    let b = ball(num(3), mag(TOP, -9)?, 53)?;
    want_exact("5-3", &ok("sub", a.sub(&b, 53))?, &num(2), &mag(TOP, -8)?)?;
    // 3 - 3 = 0 exactly, +0 midpoint.
    let z = ok("sub", pt(3, 53)?.sub(&pt(3, 53)?, 53))?;
    want_exact("3-3", &z, &num(0), &Mag::zero())?;
    if z.mid().is_negative() {
        return Err("3-3 gave -0".into());
    }
    Ok(())
}

#[test]
fn neg_rows() -> Result<(), String> {
    // -(3 ± 2^-10) = -3 ± 2^-10; -(0 ± 0) = +0 ± 0 (CONTRACT 6.3); prec kept.
    let b = -ball(num(3), mag(TOP, -9)?, 53)?;
    want_exact("-(3)", &b, &num(-3), &mag(TOP, -9)?)?;
    let z = -pt(0, 7)?;
    if !z.mid().is_zero() || z.mid().is_negative() || z.prec() != 7 || !z.rad().is_zero() {
        return Err(format!("-(0) = {:?}", z));
    }
    let back = -(-ball(num(-5), mag(TOP, 1)?, 53)?);
    want_exact("--(-5)", &back, &num(-5), &mag(TOP, 1)?)
}

#[test]
fn mul_rows() -> Result<(), String> {
    // (3 ± 2^-10)(5 ± 2^-10) = 15 ± (3·2^-10 + 5·2^-10 + 2^-20) = 15 ± (2^-7 + 2^-20).
    // 2^-7 + 2^-20 = 2^-7 (1 + 2^-13): mant TOP | TOP >> 13, exp -6.
    let a = ball(num(3), mag(TOP, -9)?, 53)?;
    let b = ball(num(5), mag(TOP, -9)?, 53)?;
    let want_rad = mag(TOP | (TOP >> 13), -6)?;
    want_exact("3*5", &ok("mul", a.mul(&b, 53))?, &num(15), &want_rad)?;
    // mixed signs: (-3)(5) = -15 exact.
    want_exact(
        "-3*5",
        &ok("mul", pt(-3, 53)?.mul(&pt(5, 53)?, 53))?,
        &num(-15),
        &Mag::zero(),
    )?;
    // zero times a huge radius (beyond zf's range) with no radius on the zero: 0 ± (0·R + 0·0 + 0·R) = 0 ± 0.
    let huge = ball(num(1), mag(TOP, 1 << 40)?, 53)?;
    want_exact(
        "0*huge",
        &ok("mul", pt(0, 53)?.mul(&huge, 53))?,
        &num(0),
        &Mag::zero(),
    )
}

#[test]
fn mul_exact_rows() -> Result<(), String> {
    // (3 ± 2^-10) · (-4) = -12 ± 2^-8.
    let a = ball(num(3), mag(TOP, -9)?, 53)?;
    want_exact(
        "3*-4",
        &ok("mul_exact", a.mul_exact(&num(-4), 53))?,
        &num(-12),
        &mag(TOP, -7)?,
    )?;
    // NaN and Inf scalars are Representation.
    for c in [
        ExactNum::nan(Some(zenith_float_num::Error::InvalidArgument)),
        ExactNum::from(1i64).div(&num(0), 64, RoundingMode::ToEven),
    ] {
        match pt(1, 53)?.mul_exact(&c, 53) {
            Err(BoundError::Representation { .. }) => {}
            other => return Err(format!("mul_exact by {:?}: {:?}", c, other)),
        }
    }
    Ok(())
}

#[test]
fn sqr_rows() -> Result<(), String> {
    // (3 ± 2^-10)^2 = 9 ± (2·3·2^-10 + 2^-20) = 9 ± (3·2^-9 + 2^-20); value 2^-7 (0.75 + 2^-13): mant (3<<62)|(1<<51), exp -7.
    let a = ball(num(3), mag(TOP, -9)?, 53)?;
    want_exact(
        "3^2",
        &ok("sqr", a.sqr(53))?,
        &num(9),
        &mag(0xC000_0000_0000_0000 | (1 << 51), -7)?,
    )?;
    // (-3)^2 = 9 exact.
    want_exact(
        "(-3)^2",
        &ok("sqr", pt(-3, 53)?.sqr(53))?,
        &num(9),
        &Mag::zero(),
    )?;
    // 3^2 = 9 at 2 bits rounds to 8 (ties: 9 = 1001b, nearest 2-bit values 8 and 12): ε = half_ulp(4, 2) = 2^1.
    let b = ok("sqr", pt(3, 2)?.sqr(2))?;
    want_exact("3^2 @2", &b, &num(8), &mag(TOP, 2)?)?;
    want_contains("3^2 @2", &b, &num(9))
}

#[test]
fn powi_rows() -> Result<(), String> {
    // 3^5 = 243, 2^10 = 1024, (-2)^3 = -8: exact.
    want_exact(
        "3^5",
        &ok("powi", pt(3, 53)?.powi(5, 53))?,
        &num(243),
        &Mag::zero(),
    )?;
    want_exact(
        "2^10",
        &ok("powi", pt(2, 53)?.powi(10, 53))?,
        &num(1024),
        &Mag::zero(),
    )?;
    want_exact(
        "(-2)^3",
        &ok("powi", pt(-2, 53)?.powi(3, 53))?,
        &num(-8),
        &Mag::zero(),
    )?;
    want_exact(
        "7^1",
        &ok("powi", pt(7, 53)?.powi(1, 53))?,
        &num(7),
        &Mag::zero(),
    )?;
    // x^0 = 1 ± 0 for every ball, including 0 and a ball containing 0 (CONTRACT 10.2e).
    want_exact(
        "0^0",
        &ok("powi", pt(0, 53)?.powi(0, 53))?,
        &num(1),
        &Mag::zero(),
    )?;
    let wide = ball(num(0), mag(TOP, 5)?, 53)?;
    want_exact(
        "(0±16)^0",
        &ok("powi", wide.powi(0, 53))?,
        &num(1),
        &Mag::zero(),
    )?;
    // 3^41 > 2^53: inexact at 53 bits, still encloses 3^41 exactly.
    let mut p = num(1);
    for _ in 0..41 {
        p = p.mul_full_prec(&num(3));
    }
    let b = ok("powi", pt(3, 53)?.powi(41, 53))?;
    if b.is_exact() {
        return Err("3^41 at 53 bits cannot be exact".into());
    }
    want_contains("3^41", &b, &p)?;
    // (1 ± 2^-10)^2 = 1 ± (2^-9 + 2^-20) by squaring: contains 0.999^2 and 1.001^2 endpoints exactly.
    let a = ball(num(1), mag(TOP, -9)?, 53)?;
    let b = ok("powi", a.powi(2, 53))?;
    let lo = num(1).add_full_prec(&-pow2(-10)?);
    let hi = num(1).add_full_prec(&pow2(-10)?);
    want_contains("(1-r)^2", &b, &lo.mul_full_prec(&lo))?;
    want_contains("(1+r)^2", &b, &hi.mul_full_prec(&hi))
}

#[test]
fn sqrt_rows() -> Result<(), String> {
    let strict = DomainPolicy::Strict;
    let clip = DomainPolicy::Clip;
    // sqrt(4) = 2, sqrt(0) = 0: exact.
    want_exact(
        "sqrt 4",
        &ok("sqrt", pt(4, 53)?.sqrt(53, &strict))?,
        &num(2),
        &Mag::zero(),
    )?;
    want_exact(
        "sqrt 0",
        &ok("sqrt", pt(0, 53)?.sqrt(53, &strict))?,
        &num(0),
        &Mag::zero(),
    )?;
    // sqrt(2) at 53 bits brackets the true root.
    let b = ok("sqrt", pt(2, 53)?.sqrt(53, &strict))?;
    want_contains("sqrt2 lo", &b, &num(2).sqrt(256, RoundingMode::Down))?;
    want_contains("sqrt2 hi", &b, &num(2).sqrt(256, RoundingMode::Up))?;
    // (1 ± 2): lo = -1 < 0: Strict is Domain; Clip encloses sqrt over [0, 3].
    want_domain(
        "strict 1±2",
        ball(num(1), mag(TOP, 2)?, 53)?.sqrt(53, &strict),
    )?;
    let c = ok("clip 1±2", ball(num(1), mag(TOP, 2)?, 53)?.sqrt(53, &clip))?;
    want_contains("clip 0", &c, &num(0))?;
    want_contains("clip sqrt3", &c, &num(3).sqrt(256, RoundingMode::Up))?;
    // (1 ± 1): lo = 0 exactly, inside the closed domain: the endpoint ball [0, sqrt_up(2)] (D15d).
    let e = ok(
        "strict 1±1",
        ball(num(1), mag(TOP, 1)?, 53)?.sqrt(53, &strict),
    )?;
    want_contains("endpoint 0", &e, &num(0))?;
    want_contains("endpoint sqrt2", &e, &num(2).sqrt(256, RoundingMode::Up))?;
    if !e.is_informative() {
        return Err("endpoint ball must be finite".into());
    }
    // (-3 ± 1): misses [0, inf): Domain under both. (-1 ± 1): touches 0 only: Strict Domain, Clip exactly 0.
    want_domain(
        "strict -3±1",
        ball(num(-3), mag(TOP, 1)?, 53)?.sqrt(53, &strict),
    )?;
    want_domain(
        "clip -3±1",
        ball(num(-3), mag(TOP, 1)?, 53)?.sqrt(53, &clip),
    )?;
    want_domain(
        "strict -1±1",
        ball(num(-1), mag(TOP, 1)?, 53)?.sqrt(53, &strict),
    )?;
    want_exact(
        "clip -1±1",
        &ok(
            "clip -1±1",
            ball(num(-1), mag(TOP, 1)?, 53)?.sqrt(53, &clip),
        )?,
        &num(0),
        &Mag::zero(),
    )?;
    // (0 ± 2^-10): Strict Domain, Clip [0, 2^-5].
    want_domain(
        "strict 0±r",
        ball(num(0), mag(TOP, -9)?, 53)?.sqrt(53, &strict),
    )?;
    let z = ok("clip 0±r", ball(num(0), mag(TOP, -9)?, 53)?.sqrt(53, &clip))?;
    want_contains("clip 0±r top", &z, &pow2(-5)?)?;
    // the all-reals ball: Strict Domain, Clip all reals (10.2a).
    want_domain(
        "strict all",
        ok("all", Ball::all_reals(53))?.sqrt(53, &strict),
    )?;
    if ok("clip all", ok("all", Ball::all_reals(53))?.sqrt(53, &clip))?.is_informative() {
        return Err("sqrt Clip of all reals must be all reals".into());
    }
    // a radius below zf's range on a positive mid: lo > 0, decided through the Mag bounds.
    let t = ok(
        "tiny r",
        ball(num(4), mag(TOP, -(1 << 40))?, 53)?.sqrt(53, &strict),
    )?;
    want_contains("tiny r", &t, &num(2))
}

#[test]
fn sqrt_lo_zero_is_the_endpoint_ball() -> Result<(), String> {
    // 4 ± 4 has lo == 0: Strict accepts it, and the result is the endpoint ball [0, sqrt_up(8)] (16.S1.7), whose top
    // stays below 3. The lo > 0 formula would give 2 ± 2, reaching 4.
    let b = ok(
        "sqrt",
        ball(num(4), mag(TOP, 3)?, 53)?.sqrt(53, &DomainPolicy::Strict),
    )?;
    want_contains("sqrt lo 0", &b, &num(0))?;
    let top = Mag::from_exact_num(b.mid()).add_up(b.rad());
    if top >= Mag::from_exact_num(&num(3)) {
        return Err(format!(
            "sqrt(4 ± 4) reaches {:?}, want below 3: {:?}",
            top, b
        ));
    }
    Ok(())
}

#[test]
fn prec_out_of_range_is_invalid_argument() -> Result<(), String> {
    let checks = [
        pt(1, 53)?.add(&pt(1, 53)?, 0),
        pt(1, 53)?.sub(&pt(1, 53)?, 0),
        pt(1, 53)?.mul(&pt(1, 53)?, 0),
        pt(1, 53)?.mul_exact(&num(1), 0),
        pt(1, 53)?.sqr(0),
        pt(1, 53)?.powi(0, 0),
        pt(1, 53)?.sqrt(0, &DomainPolicy::Strict),
        pt(1, 53)?.add(
            &pt(1, 53)?,
            zenith_bound::MAX_PRECISION.checked_add(1).ok_or("max+1")?,
        ),
    ];
    for (i, r) in checks.into_iter().enumerate() {
        match r {
            Err(BoundError::InvalidArgument { .. }) => {}
            other => return Err(format!("case {}: want InvalidArgument, got {:?}", i, other)),
        }
    }
    Ok(())
}

#[test]
fn overflow_is_all_reals() -> Result<(), String> {
    // 2^(EMAX-1) + 2^(EMAX-1) = 2^EMAX overflows zf; so do its product and square, and powi of 2^(EMAX/2+1).
    let big = || -> Result<Ball, String> { ok("big", Ball::exact(pow2(EXPONENT_MAX - 1)?, 53)) };
    let results = [
        big()?.add(&big()?, 53),
        big()?.sub(&-big()?, 53),
        big()?.mul(&big()?, 53),
        big()?.mul_exact(&num(4), 53),
        big()?.sqr(53),
        ok("half", Ball::exact(pow2(EXPONENT_MAX / 2 + 1)?, 53))?.powi(2, 53),
    ];
    for (i, r) in results.into_iter().enumerate() {
        let b = ok("overflow", r)?;
        if b.is_informative() {
            return Err(format!("case {}: want all reals, got {:?}", i, b));
        }
    }
    Ok(())
}

#[test]
fn underflow_is_exact_or_enclosed() -> Result<(), String> {
    // 2^(2k) for k near EMIN/2 lies in or below zf's bottom binade. zf keeps a shallow one exactly (EMIN with leading
    // zero bits: rad 0 is then right, and the midpoint must be exactly 2^(2k)); a deeper one is lost, and then the radius
    // must be positive and the ball must reach from 0 up past its midpoint (16.S1.2: bottom ε). Never a silent flush.
    let half = EXPONENT_MIN / 2;
    for k in [half - 2, half - 40, half - 200, half - 5000] {
        let x = || -> Result<Ball, String> { ok("tiny", Ball::exact(pow2(k)?, 53)) };
        for (what, r) in [("mul", x()?.mul(&x()?, 53)), ("sqr", x()?.sqr(53))] {
            let what = format!("{} k {}", what, k);
            let b = ok(&what, r)?;
            if b.rad().is_zero() {
                let mut exact = pow2(k)?;
                exact.set_inexact(false);
                let sq = exact.mul_full_prec(&exact);
                if sq.inexact_accumulated()
                    || sq.is_zero()
                    || b.mid().abs_cmp(&sq) != Some(0)
                    || b.mid().is_negative()
                {
                    return Err(format!(
                        "{}: rad 0 but the midpoint is not exactly 2^(2k): {:?}",
                        what, b
                    ));
                }
            } else {
                want_contains(&what, &b, &num(0))?;
                want_contains(&what, &b, b.mid())?;
            }
        }
    }
    Ok(())
}

#[test]
fn huge_and_infinite_radius() -> Result<(), String> {
    // a radius far beyond zf's range survives add (no overflow of the midpoint) and contains 1e300-scale points.
    let h = ball(num(1), mag(TOP, 1 << 40)?, 53)?;
    let b = ok("huge add", h.add(&pt(1, 53)?, 53))?;
    if !b.is_informative() || b.rad() < &mag(TOP, 1 << 40)? {
        return Err(format!("huge add: {:?}", b));
    }
    want_contains("huge add", &b, &pow2(997)?)?;
    // an all-reals operand gives all reals.
    for r in [
        ok("all", Ball::all_reals(53))?.add(&pt(1, 53)?, 53),
        pt(1, 53)?.mul(&ok("all", Ball::all_reals(53))?, 53),
        ok("all", Ball::all_reals(53))?.sqr(53),
        ok("all", Ball::all_reals(53))?.powi(3, 53),
    ] {
        if ok("all", r)?.is_informative() {
            return Err("all-reals operand must give all reals".into());
        }
    }
    Ok(())
}

#[test]
fn zero_radius_inexact_mid_gets_eps() -> Result<(), String> {
    // 1/3-like: (2^60 + 1) · 3 at 53 bits is inexact, so rad > 0 even with exact inputs (CONTRACT 2.2).
    let x = num(1).add_full_prec(&pow2(60)?);
    let b = ok("mul", ok("x", Ball::exact(x, 53))?.mul(&pt(3, 53)?, 53))?;
    if b.rad().is_zero() {
        return Err("inexact product with rad 0".into());
    }
    want_contains("inexact product", &b, &num(3).add_full_prec(&dy(3, 60)?))
}

// ---- containment property tests ----------------------------------------------------------------------------------

/// A random midpoint near `2^e`: a 64-bit top chunk plus a 64-bit chunk `k` bits further down (so wide mids round at
/// high precisions too), random sign.
fn rand_mid(rng: &mut SplitMix64, e: i64, prec: usize) -> Result<ExactNum, String> {
    let span = if prec > 1000 { 20_000 } else { 200 };
    let k = 64i64
        .checked_add(i64::try_from(rng.below(span)).map_err(|_| "span")?)
        .ok_or("k")?;
    let hi =
        ExactNum::from_mag_parts(rng.next_u64() | TOP, e, 64).map_err(|e| format!("{:?}", e))?;
    let lo = ExactNum::from_mag_parts(rng.next_u64(), e.checked_sub(k).ok_or("e-k")?, 64)
        .map_err(|e| format!("{:?}", e))?;
    let m = hi.add_full_prec(&lo);
    Ok(if rng.below(2) == 1 { -m } else { m })
}

/// A random radius: zero, tiny, comparable, or larger than the midpoint.
fn rand_rad(rng: &mut SplitMix64, e: i64) -> Result<Mag, String> {
    let off = match rng.below(4) {
        0 => return Ok(Mag::zero()),
        1 => -300 - i64::try_from(rng.below(300)).map_err(|_| "off")?,
        2 => -i64::try_from(rng.below(60)).map_err(|_| "off")?,
        _ => i64::try_from(rng.below(8)).map_err(|_| "off")?,
    };
    mag(rng.next_u64() | TOP, e.checked_add(off).ok_or("rad exp")?)
}

fn rand_ball(rng: &mut SplitMix64, prec: usize) -> Result<Ball, String> {
    let i = usize::try_from(rng.below(13)).map_err(|_| "i")?;
    let e = *EXPS.get(i).ok_or("EXPS index")?;
    let m = rand_mid(rng, e, prec)?;
    let r = rand_rad(rng, e)?;
    ball(m, r, prec)
}

/// A random ball whose midpoint exponent lies within 200 of `e` (add/sub operands, see run_binary).
fn near_ball(rng: &mut SplitMix64, e: i64, prec: usize) -> Result<Ball, String> {
    let off = i64::try_from(rng.below(401))
        .map_err(|_| "off")?
        .checked_sub(200)
        .ok_or("off")?;
    let e = e.checked_add(off).ok_or("e")?;
    let m = rand_mid(rng, e, prec)?;
    let r = rand_rad(rng, e)?;
    ball(m, r, prec)
}

/// Exact sample points of a ball: m, m - r, m + r, m + r/2 (only m when r is outside zf's range).
fn points(b: &Ball) -> Result<Vec<ExactNum>, String> {
    let m = b.mid();
    let mut v = vec![m.add_full_prec(&num(0))];
    if let Ok(r) = b.rad().to_exact_num() {
        v.push(m.add_full_prec(&-r.add_full_prec(&num(0))));
        v.push(m.add_full_prec(&r));
        v.push(m.add_full_prec(&r.mul_full_prec(&pow2(-1)?)));
    }
    Ok(v)
}

fn check(what: &str, res: &Ball, x: &ExactNum) -> Result<(), String> {
    if x.is_nan() || x.is_inf() {
        return Ok(()); // the exact image left zf's range; nothing to compare
    }
    want_contains(what, res, x)
}

#[derive(Debug)]
enum Op {
    Add,
    Sub,
    Mul,
    MulExact,
    Sqr,
    Powi,
}

fn run_binary(op: &Op, seed: u64, cases: usize) -> Result<(), String> {
    let mut rng = SplitMix64::new(seed);
    for case in 0..cases {
        let prec = *PRECS
            .get(usize::try_from(rng.below(9)).map_err(|_| "p")?)
            .ok_or("PRECS index")?;
        let a = rand_ball(&mut rng, prec)?;
        // add/sub: b's exponent stays within 200 of a's. A gap of 2^29 made the exact oracle sum a 2^29-bit number
        // (minutes per case) without testing anything new: past prec + 64 bits the smaller operand only feeds ε.
        let b = match op {
            Op::Add | Op::Sub => {
                near_ball(&mut rng, a.mid().exponent().map_or(0, i64::from), prec)?
            }
            Op::Mul | Op::MulExact | Op::Sqr | Op::Powi => rand_ball(&mut rng, prec)?,
        };
        let (pa, pb) = (points(&a)?, points(&b)?);
        let n = rng.below(7).checked_add(1).ok_or("n")?;
        let what = format!(
            "{:?} case {} prec {} a {:?} b {:?} n {}",
            op, case, prec, a, b, n
        );
        let bm = b.mid().add_full_prec(&num(0));
        let res = ok(
            &what,
            match op {
                Op::Add => a.add(&b, prec),
                Op::Sub => a.sub(&b, prec),
                Op::Mul => a.mul(&b, prec),
                Op::MulExact => a.mul_exact(&bm, prec),
                Op::Sqr => a.sqr(prec),
                Op::Powi => a.powi(n, prec),
            },
        )?;
        for x in &pa {
            match op {
                Op::Add | Op::Sub | Op::Mul => {
                    for y in &pb {
                        let img = match op {
                            Op::Add => x.add_full_prec(y),
                            Op::Sub => x.add_full_prec(&-y.add_full_prec(&num(0))),
                            Op::Mul | Op::MulExact | Op::Sqr | Op::Powi => x.mul_full_prec(y),
                        };
                        check(&what, &res, &img)?;
                    }
                }
                Op::MulExact => check(&what, &res, &x.mul_full_prec(&bm))?,
                Op::Sqr => check(&what, &res, &x.mul_full_prec(x))?,
                Op::Powi => {
                    let mut p = num(1);
                    for _ in 0..n {
                        p = p.mul_full_prec(x);
                    }
                    check(&what, &res, &p)?;
                }
            }
        }
    }
    Ok(())
}

#[test]
fn containment_add() -> Result<(), String> {
    run_binary(&Op::Add, 0x5101, 300)
}

#[test]
fn containment_sub() -> Result<(), String> {
    run_binary(&Op::Sub, 0x5102, 300)
}

#[test]
fn containment_mul() -> Result<(), String> {
    run_binary(&Op::Mul, 0x5103, 300)
}

#[test]
fn containment_mul_exact() -> Result<(), String> {
    run_binary(&Op::MulExact, 0x5104, 300)
}

#[test]
fn containment_sqr() -> Result<(), String> {
    run_binary(&Op::Sqr, 0x5105, 300)
}

#[test]
fn containment_powi() -> Result<(), String> {
    run_binary(&Op::Powi, 0x5106, 200)
}

#[test]
fn containment_neg() -> Result<(), String> {
    let mut rng = SplitMix64::new(0x5107);
    for case in 0..200 {
        let prec = *PRECS
            .get(usize::try_from(rng.below(9)).map_err(|_| "p")?)
            .ok_or("PRECS index")?;
        let a = rand_ball(&mut rng, prec)?;
        let pa = points(&a)?;
        let res = -a;
        for x in &pa {
            check(
                &format!("neg case {}", case),
                &res,
                &-x.add_full_prec(&num(0)),
            )?;
        }
    }
    Ok(())
}

#[test]
fn containment_sqrt() -> Result<(), String> {
    let mut rng = SplitMix64::new(0x5108);
    for case in 0..300 {
        let prec = *PRECS
            .get(usize::try_from(rng.below(9)).map_err(|_| "p")?)
            .ok_or("PRECS index")?;
        let a = rand_ball(&mut rng, prec)?;
        let pa = points(&a)?;
        let lo_neg = match a.rad().to_exact_num() {
            Ok(r) => a.mid().add_full_prec(&-r).is_negative(),
            Err(_) => a.mid().is_negative() || a.rad() > &Mag::from_exact_num(a.mid()),
        };
        let strict_pick = rng.below(2) == 1;
        let policy = if strict_pick {
            DomainPolicy::Strict
        } else {
            DomainPolicy::Clip
        };
        let what = format!(
            "sqrt case {} prec {} {:?} strict {}",
            case, prec, a, strict_pick
        );
        let mid_was_zero = a.mid().is_zero();
        match a.sqrt(prec, &policy) {
            Err(BoundError::Domain { .. }) => {
                // Strict refuses iff m - r < 0; Clip only when no point is >= 0.
                if !lo_neg
                    || (!strict_pick
                        && pa.iter().any(|x| !x.is_negative() || x.is_zero())
                        && !mid_was_zero)
                {
                    return Err(format!("{}: unexpected Domain", what));
                }
            }
            Err(e) => return Err(format!("{}: {:?}", what, e)),
            Ok(res) => {
                if strict_pick && lo_neg {
                    return Err(format!("{}: Strict accepted m - r < 0", what));
                }
                for x in pa.iter().filter(|x| !x.is_negative() || x.is_zero()) {
                    // The oracle brackets sqrt(x) 128 bits finer than BOTH the input and the target precision (case 69: a
                    // 3264-bit x at prec 16384 was bracketed at only 3392 bits, far coarser than the ball's half ulp).
                    let p = x
                        .precision()
                        .unwrap_or(64)
                        .max(prec)
                        .checked_add(128)
                        .ok_or("p")?;
                    check(&what, &res, &x.sqrt(p, RoundingMode::Down))?;
                    check(&what, &res, &x.sqrt(p, RoundingMode::Up))?;
                }
            }
        }
    }
    Ok(())
}

// ---- the L0b enclosure rows (tests/goldens/ball_enclosure_rows.txt) ------------------------------------------------

/// A row value `SIGN EXP MANT` (mantissa hex, odd): `(-1)^SIGN * 0.MANT * 2^EXP`.
fn row_val(s: &str, e: &str, m: &str) -> Result<ExactNum, String> {
    if m == "0" {
        return Ok(num(0));
    }
    let mut x = num(0);
    let digits = m.as_bytes();
    for chunk in digits.chunks(15) {
        let text = core::str::from_utf8(chunk).map_err(|e| e.to_string())?;
        let v = u64::from_str_radix(text, 16).map_err(|e| format!("{}: {}", m, e))?;
        let shift = i64::try_from(chunk.len())
            .map_err(|_| "len")?
            .checked_mul(4)
            .ok_or("shift")?;
        x = x
            .mul_full_prec(&pow2(shift)?)
            .add_full_prec(&ExactNum::from(v));
    }
    let exp: i64 = e.parse().map_err(|_| format!("exp {}", e))?;
    x.set_exponent(exp);
    Ok(if s == "1" { -x } else { x })
}

#[test]
fn enclosure_rows_contained() -> Result<(), String> {
    let path = golden_dir().join("ball_enclosure_rows.txt");
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{:?}: {}", path, e))?;
    let mut checked = 0usize;
    for line in text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
    {
        let f: Vec<&str> = line.split_whitespace().collect();
        let op: &str = f.first().ok_or("empty row")?;
        let prec: usize = f.get(1).ok_or("prec")?.parse().map_err(|_| "prec")?;
        let v = |i: usize| -> Result<ExactNum, String> {
            match (f.get(i), f.get(i + 1), f.get(i + 2)) {
                (Some(s), Some(e), Some(m)) => row_val(s, e, m),
                _ => Err(format!("short row: {}", line)),
            }
        };
        let rad = |i: usize| -> Result<Mag, String> { Ok(Mag::from_exact_num(&v(i)?)) };
        let (res, lo_at) = match op {
            "sqrt" => (
                ball(v(3)?, rad(6)?, prec)?.sqrt(prec, &DomainPolicy::Clip),
                9,
            ),
            "add" | "sub" | "mul" | "div" => {
                let a = ball(v(3)?, rad(6)?, prec)?;
                let b = ball(v(9)?, rad(12)?, prec)?;
                let r = match op {
                    "add" => a.add(&b, prec),
                    "sub" => a.sub(&b, prec),
                    "mul" => a.mul(&b, prec),
                    _ => a.div(&b, prec),
                };
                (r, 15)
            }
            _ => continue,
        };
        let res = ok(line, res)?;
        want_contains(line, &res, &v(lo_at)?)?;
        want_contains(line, &res, &v(lo_at + 3)?)?;
        checked += 1;
    }
    if checked != 1185 {
        return Err(format!(
            "checked {} rows, want 1185 (237 each of add, sub, mul, div, sqrt)",
            checked
        ));
    }
    Ok(())
}

// ---- Harper L0b harness (child modules under tests/harness/slice_arith/) -------------------------------------------

mod harper_cmp;
mod harper_edges;
mod harper_far_side_edge;
mod harper_oracle;
mod harper_props;

// ---- div, inv, hypot (16.S1.8, zenith-float R13) -----------------------------------------------------------------

/// `[lo, hi]` brackets of an exact quotient or hypot at 256 bits, directed by zenith-float (Down / Up).
fn brackets(f: &dyn Fn(RoundingMode) -> ExactNum) -> [ExactNum; 2] {
    [f(RoundingMode::Down), f(RoundingMode::Up)]
}

fn want_all_reals(what: &str, r: Result<Ball, BoundError>) -> Result<(), String> {
    let b = ok(what, r)?;
    if b.is_informative() {
        return Err(format!("{}: want all reals, got {:?}", what, b));
    }
    Ok(())
}

fn want_invalid(what: &str, r: Result<Ball, BoundError>) -> Result<(), String> {
    match r {
        Err(BoundError::InvalidArgument { .. }) => Ok(()),
        other => Err(format!("{}: want InvalidArgument, got {:?}", what, other)),
    }
}

#[test]
fn div_rows() -> Result<(), String> {
    // 6 / 3 = 2 exactly; -7 / 2 = -3.5 exactly; radius 0.
    want_exact(
        "6/3",
        &ok("6/3", pt(6, 53)?.div(&pt(3, 53)?, 53))?,
        &num(2),
        &Mag::zero(),
    )?;
    want_exact(
        "-7/2",
        &ok("-7/2", pt(-7, 53)?.div(&pt(2, 53)?, 53))?,
        &dy(-7, -1)?,
        &Mag::zero(),
    )?;
    // 1/3 at 53 bits: inexact, so the radius is positive and the ball holds both 256-bit brackets of 1/3.
    let b = ok("1/3", pt(1, 53)?.div(&pt(3, 53)?, 53))?;
    if b.rad().is_zero() {
        return Err(format!("1/3: zero radius {:?}", b));
    }
    for x in brackets(&|rm| num(1).div(&num(3), 256, rm)) {
        want_contains("1/3", &b, &x)?;
    }
    // (6 ± 1) / (3 ± 1) = [5/4, 7/2]: both endpoints are exact dyadics and must be inside.
    let b = ok(
        "(6±1)/(3±1)",
        ball(num(6), mag(TOP, 1)?, 53)?.div(&ball(num(3), mag(TOP, 1)?, 53)?, 53),
    )?;
    want_contains("5/4", &b, &dy(5, -2)?)?;
    want_contains("7/2", &b, &dy(7, -1)?)?;
    // -(6 ± 1) / (3 ± 1) = [-7/2, -5/4] and (6 ± 1) / -(3 ± 1) the same.
    let b = ok(
        "neg num",
        ball(num(-6), mag(TOP, 1)?, 53)?.div(&ball(num(3), mag(TOP, 1)?, 53)?, 53),
    )?;
    want_contains("-7/2", &b, &dy(-7, -1)?)?;
    want_contains("-5/4", &b, &dy(-5, -2)?)?;
    let b = ok(
        "neg den",
        ball(num(6), mag(TOP, 1)?, 53)?.div(&ball(num(-3), mag(TOP, 1)?, 53)?, 53),
    )?;
    want_contains("-7/2 b", &b, &dy(-7, -1)?)?;
    want_contains("-5/4 b", &b, &dy(-5, -2)?)?;
    // 1 / (4 ± 1) = [1/5, 1/3]: brackets of both endpoints.
    let b = ok(
        "1/(4±1)",
        pt(1, 53)?.div(&ball(num(4), mag(TOP, 1)?, 53)?, 53),
    )?;
    for x in brackets(&|rm| num(1).div(&num(5), 256, rm)) {
        want_contains("1/5", &b, &x)?;
    }
    for x in brackets(&|rm| num(1).div(&num(3), 256, rm)) {
        want_contains("1/3 end", &b, &x)?;
    }
    Ok(())
}

#[test]
fn div_by_zero_and_zero_straddling() -> Result<(), String> {
    want_domain("1/0", pt(1, 53)?.div(&pt(0, 53)?, 53))?;
    want_domain("0/0", pt(0, 53)?.div(&pt(0, 53)?, 53))?;
    want_all_reals(
        "1/(0±1)",
        pt(1, 53)?.div(&ball(num(0), mag(TOP, 1)?, 53)?, 53),
    )?;
    // |n| == s: the denominator reaches 0 exactly.
    want_all_reals(
        "1/(1±1)",
        pt(1, 53)?.div(&ball(num(1), mag(TOP, 1)?, 53)?, 53),
    )?;
    want_all_reals(
        "1/(-1±1)",
        pt(1, 53)?.div(&ball(num(-1), mag(TOP, 1)?, 53)?, 53),
    )?;
    want_all_reals(
        "1/all",
        pt(1, 53)?.div(&ok("all", Ball::all_reals(53))?, 53),
    )?;
    want_all_reals(
        "all/2",
        ok("all", Ball::all_reals(53))?.div(&pt(2, 53)?, 53),
    )?;
    // 0 / 2 = 0 exactly.
    want_exact(
        "0/2",
        &ok("0/2", pt(0, 53)?.div(&pt(2, 53)?, 53))?,
        &num(0),
        &Mag::zero(),
    )?;
    Ok(())
}

#[test]
fn inv_rows() -> Result<(), String> {
    want_exact(
        "1/4",
        &ok("inv 4", pt(4, 53)?.inv(53))?,
        &dy(1, -2)?,
        &Mag::zero(),
    )?;
    want_exact(
        "1/-2",
        &ok("inv -2", pt(-2, 53)?.inv(53))?,
        &dy(-1, -1)?,
        &Mag::zero(),
    )?;
    let b = ok("inv 3", pt(3, 53)?.inv(53))?;
    if b.rad().is_zero() {
        return Err(format!("inv 3: zero radius {:?}", b));
    }
    for x in brackets(&|rm| num(1).div(&num(3), 256, rm)) {
        want_contains("inv 3", &b, &x)?;
    }
    // 1 / (2 ± 1) = [1/3, 1].
    let b = ok("inv 2±1", ball(num(2), mag(TOP, 1)?, 53)?.inv(53))?;
    want_contains("1", &b, &num(1))?;
    for x in brackets(&|rm| num(1).div(&num(3), 256, rm)) {
        want_contains("inv 2±1 lo", &b, &x)?;
    }
    // 1 / -(4 ± 1) = [-1/3, -1/5].
    let b = ok("inv -4±1", ball(num(-4), mag(TOP, 1)?, 53)?.inv(53))?;
    for x in brackets(&|rm| num(-1).div(&num(5), 256, rm)) {
        want_contains("-1/5", &b, &x)?;
    }
    for x in brackets(&|rm| num(-1).div(&num(3), 256, rm)) {
        want_contains("-1/3", &b, &x)?;
    }
    want_domain("inv 0", pt(0, 53)?.inv(53))?;
    want_all_reals("inv 0±1", ball(num(0), mag(TOP, 1)?, 53)?.inv(53))?;
    want_all_reals("inv 1±1", ball(num(1), mag(TOP, 1)?, 53)?.inv(53))?;
    want_all_reals("inv all", ok("all", Ball::all_reals(53))?.inv(53))?;
    Ok(())
}

#[test]
fn hypot_rows() -> Result<(), String> {
    want_exact(
        "hypot 3 4",
        &ok("h34", pt(3, 53)?.hypot(&pt(4, 53)?, 53))?,
        &num(5),
        &Mag::zero(),
    )?;
    want_exact(
        "hypot -3 4",
        &ok("h-34", pt(-3, 53)?.hypot(&pt(4, 53)?, 53))?,
        &num(5),
        &Mag::zero(),
    )?;
    want_exact(
        "hypot 0 0",
        &ok("h00", pt(0, 53)?.hypot(&pt(0, 53)?, 53))?,
        &num(0),
        &Mag::zero(),
    )?;
    let b = ok("hypot 1 1", pt(1, 53)?.hypot(&pt(1, 53)?, 53))?;
    if b.rad().is_zero() {
        return Err(format!("hypot 1 1: zero radius {:?}", b));
    }
    for x in brackets(&|rm| num(1).hypot(&num(1), 256, rm)) {
        want_contains("sqrt 2", &b, &x)?;
    }
    // hypot(3 ± 1, 4 ± 1) = [hypot(2, 3), hypot(4, 5)].
    let b = ok(
        "h ball",
        ball(num(3), mag(TOP, 1)?, 53)?.hypot(&ball(num(4), mag(TOP, 1)?, 53)?, 53),
    )?;
    for x in brackets(&|rm| num(2).hypot(&num(3), 256, rm)) {
        want_contains("hypot 2 3", &b, &x)?;
    }
    for x in brackets(&|rm| num(4).hypot(&num(5), 256, rm)) {
        want_contains("hypot 4 5", &b, &x)?;
    }
    // hypot(0 ± 1, 0 ± 1) reaches sqrt 2 at the corner.
    let b = ok(
        "h zero balls",
        ball(num(0), mag(TOP, 1)?, 53)?.hypot(&ball(num(0), mag(TOP, 1)?, 53)?, 53),
    )?;
    for x in brackets(&|rm| num(1).hypot(&num(1), 256, rm)) {
        want_contains("corner", &b, &x)?;
    }
    want_all_reals(
        "hypot all",
        ok("all", Ball::all_reals(53))?.hypot(&pt(1, 53)?, 53),
    )?;
    want_all_reals(
        "hypot , all",
        pt(1, 53)?.hypot(&ok("all", Ball::all_reals(53))?, 53),
    )?;
    Ok(())
}

#[test]
fn div_inv_hypot_overflow_and_prec() -> Result<(), String> {
    let big = || -> Result<Ball, String> { ok("big", Ball::exact(pow2(EXPONENT_MAX - 1)?, 53)) };
    want_all_reals(
        "big/2^-2",
        big()?.div(&ok("q", Ball::exact(pow2(-2)?, 53))?, 53),
    )?;
    want_all_reals(
        "inv tiny",
        ok("tiny", Ball::exact(pow2(EXPONENT_MIN)?, 53))?.inv(53),
    )?;
    want_all_reals("hypot big big", big()?.hypot(&big()?, 53))?;
    // A quotient in the subnormal range is exact: 2^(EMIN+1) / 2^10, checked by multiplying back.
    let b = ok(
        "under",
        ok("t", Ball::exact(pow2(EXPONENT_MIN + 1)?, 53))?.div(&pt(1024, 53)?, 53),
    )?;
    if !b.rad().is_zero() || b.mid().mul_full_prec(&num(1024)) != pow2(EXPONENT_MIN + 1)? {
        return Err(format!("subnormal quotient: {:?}", b));
    }
    // 2^EMIN / 3 loses bits below the grid; the ball is still a sound enclosure of a positive value.
    let b = ok(
        "under3",
        ok("t", Ball::exact(pow2(EXPONENT_MIN)?, 53))?.div(&pt(3, 53)?, 53),
    )?;
    if b.rad().is_zero() {
        return Err(format!("2^EMIN/3: zero radius {:?}", b));
    }
    for p in [0usize, usize::MAX] {
        want_invalid("div prec", pt(1, 53)?.div(&pt(3, 53)?, p))?;
        want_invalid("inv prec", pt(3, 53)?.inv(p))?;
        want_invalid("hypot prec", pt(3, 53)?.hypot(&pt(4, 53)?, p))?;
    }
    Ok(())
}

/// CONTRACT 8.3 / O3 / O21: zenith-float's `Exponent` is `i64` with range ±(2^62−1); both ends are reachable midpoint
/// exponents and one past either end is refused.
#[test]
fn exponent_min_max_are_i64_62_bit() -> Result<(), String> {
    let max: i64 = EXPONENT_MAX;
    if max != (1i64 << 62) - 1 || EXPONENT_MIN != -max {
        return Err(format!(
            "EXPONENT_MIN {} EXPONENT_MAX {}",
            EXPONENT_MIN, max
        ));
    }
    // 0.1b * 2^e at both ends builds and keeps its exponent; one beyond is refused.
    for e in [EXPONENT_MIN, EXPONENT_MAX] {
        let x =
            ExactNum::from_mag_parts(TOP, e, 64).map_err(|err| format!("exp {}: {:?}", e, err))?;
        if x.exponent() != Some(e) {
            return Err(format!("exp {}: got {:?}", e, x.exponent()));
        }
        ok("ball at the end", Ball::exact(x, 53))?;
    }
    if ExactNum::from_mag_parts(TOP, EXPONENT_MAX + 1, 64).is_ok() {
        return Err("EXPONENT_MAX + 1 accepted".into());
    }
    // A Mag one binade above the zf range has no exact zenith-float image.
    match mag(TOP, EXPONENT_MAX + 1)?.to_exact_num() {
        Err(BoundError::Representation { .. }) => Ok(()),
        other => Err(format!("Mag above EXPONENT_MAX: {:?}", other.map(|_| ()))),
    }
}

/// CONTRACT 8.3 / 10.4: a radius that would fall below the `Mag` range rounds outward to the smallest positive `Mag`,
/// never to zero. `2^EXPONENT_MIN · s` with `s` near the bottom of the `Mag` range underflows `|n| r`.
#[test]
fn radius_below_exponent_min_rounds_out_to_smallest_mag() -> Result<(), String> {
    let smallest = mag(TOP, i64::MIN)?;
    let tiny = ok("tiny", Ball::exact(pow2(EXPONENT_MIN)?, 53))?;
    let low = ball(num(1), mag(TOP, i64::MIN + 1)?, 53)?;
    for (what, r) in [
        ("mul", low.mul(&tiny, 53)),
        (
            "mul_exact",
            ball(num(1), mag(TOP, i64::MIN + 1)?, 53)?.mul_exact(&pow2(EXPONENT_MIN)?, 53),
        ),
    ] {
        let b = ok(what, r)?;
        if b.mid() != &pow2(EXPONENT_MIN)? || b.rad() != &smallest {
            return Err(format!("{}: {:?}", what, b));
        }
    }
    let q = ok(
        "div",
        ball(num(1), mag(TOP, i64::MIN + 1)?, 53)?
            .div(&ok("big", Ball::exact(pow2(EXPONENT_MAX - 1)?, 53))?, 53),
    )?;
    if q.rad().is_zero() || q.rad() < &smallest {
        return Err(format!("div: {:?}", q));
    }
    Ok(())
}

/// 16.S1.8 edges: divisors wider than a 64-bit `Mag` (the `|n|` lower bound and the denominator must round down), and
/// the midpoint of `inv` rounded `ToEven`, including a negative divisor.
#[test]
fn div_inv_wide_divisor_and_mid_rounding() -> Result<(), String> {
    let p = 200;
    let ev = RoundingMode::ToEven;
    for (j, k) in [(1i64, 99i64), (3, 98), (-5, 97), (7, 100), (-1, 96)] {
        for sign in [1i64, -1] {
            // n = ±(2^100 + 2^37 + j), 101 significant bits; s = 2^(k-1). Built afresh for each use.
            let n = || -> Result<ExactNum, String> {
                let base =
                    ExactNum::from_mag_parts(TOP | 1, 101, 64).map_err(|e| format!("{:?}", e))?;
                Ok(base.add(&num(j), p, ev).mul(&num(sign), p, ev))
            };
            let lo = n()?.sub(&pow2(k - 1)?, p, ev);
            let hi = n()?.add(&pow2(k - 1)?, p, ev);
            let inv_b = ok("inv wide", ball(n()?, mag(TOP, k)?, p)?.inv(p))?;
            let div_b = ok("div wide", pt(3, p)?.div(&ball(n()?, mag(TOP, k)?, p)?, p))?;
            for end in [&lo, &hi] {
                for x in brackets(&|rm| num(1).div(end, 512, rm)) {
                    want_contains("inv wide end", &inv_b, &x)?;
                }
                for x in brackets(&|rm| num(3).div(end, 512, rm)) {
                    want_contains("div wide end", &div_b, &x)?;
                }
            }
        }
    }
    for v in [3i64, -3, 10, -7, 11] {
        let b = ok("inv mid", pt(v, 53)?.inv(53))?;
        let want = num(1).div(&num(v), 53, RoundingMode::ToEven);
        if b.mid() != &want {
            return Err(format!(
                "inv {}: mid {:?}, want ToEven {:?}",
                v,
                b.mid(),
                want
            ));
        }
    }
    Ok(())
}

/// 16.S1.8: the denominator `|n| (|n| - s)` is a lower bound. With `n` and `n - s` exact 64-bit values the endpoint
/// `1/(n - s)` sits exactly `s / (n (n - s))` from `1/n`, so an upward-rounded product loses it on some draws.
#[test]
fn div_inv_denominator_rounds_down() -> Result<(), String> {
    let p = 200;
    let ev = RoundingMode::ToEven;
    let mut rng = SplitMix64::new(0x5_1_8_d_e_4);
    for _ in 0..200 {
        let nm = rng.next_u64() | TOP;
        let sm = (rng.next_u64() | TOP) & !3u64;
        for sign in [1i64, -1] {
            let n = || -> Result<ExactNum, String> {
                let x = ExactNum::from_mag_parts(nm, 64, 64).map_err(|e| format!("{:?}", e))?;
                Ok(x.mul(&num(sign), p, ev))
            };
            let sx = ExactNum::from_mag_parts(sm, 62, 64).map_err(|e| format!("{:?}", e))?;
            let near = n()?.sub(&sx.mul(&num(sign), p, ev), p, ev);
            let inv_b = ok("inv", ball(n()?, mag(sm, 62)?, p)?.inv(p))?;
            let div_b = ok("div", pt(1, p)?.div(&ball(n()?, mag(sm, 62)?, p)?, p))?;
            for x in brackets(&|rm| num(1).div(&near, 512, rm)) {
                want_contains("inv near end", &inv_b, &x)?;
                want_contains("div near end", &div_b, &x)?;
            }
        }
    }
    Ok(())
}

/// Harper's oracle rows 1313 and 1343 (arith_oracle_rows.txt @ e634e11, test hypot_wide_rows_1313_1343_tight @ b1c412b):
/// `OP PREC W N AMID ARAD BMID BRAD LO HI RREF`, each value `SIGN EXP MANT`. The first ball touches or straddles 0,
/// so the 1-Lipschitz `r + s` radius is over 4x the reference radius; 16.S1.8 requires rad <= 4 * RREF.
const HYPOT_WIDE_ROWS: [(usize, &str); 2] = [
    (1313, "hypot 64 192 0 0 -28 1 0 -28 1 1 -26 7ffffffffffffff 0 -87 3 0 -26 7fffffffffffffed 0 -25 83f07b357f6837a0c7d62f13ee20be2705fe7f07ebb27d19 0 -31 5e2f274e3236ab25"),
    (1343, "hypot 256 576 0 1 -1058 489706e98a646a056377d 0 -1058 3 1 -1056 ed2d7dd0cc75befbdf7a601e714cae309d5ff5498fc1ca098cc8e4a26cfe994f 0 -1313 1 0 -1056 3b4b5f74331d6fbef7de98079c532b8c2757fd5263f07282633239289b3fa653b 0 -1056 7ddb211a35ecafa48bddea922ae94aa45a81ecab0f05f8151dfd6226c4ef912f6acdb6567175ca04c408b2ff737809a7e9979fa2387738aeacb7da04b1d58e76f25d3b763e9e9823 0 -1060 178be5214e370fb7"),
];

#[test]
fn hypot_wide_rows_1313_1343_tight() -> Result<(), String> {
    for (id, line) in HYPOT_WIDE_ROWS {
        let f: Vec<&str> = line.split_whitespace().collect();
        let v = |i: usize| -> Result<ExactNum, String> {
            match (f.get(i), f.get(i + 1), f.get(i + 2)) {
                (Some(s), Some(e), Some(m)) => row_val(s, e, m),
                _ => Err(format!("row {}: short at field {}", id, i)),
            }
        };
        let prec: usize = f
            .get(1)
            .ok_or("prec")?
            .parse()
            .map_err(|_| format!("row {}: prec", id))?;
        let a = ball(v(4)?, Mag::from_exact_num(&v(7)?), prec)?;
        let b = ball(v(10)?, Mag::from_exact_num(&v(13)?), prec)?;
        let r = ok("hypot", a.hypot(&b, prec))?;
        let what = format!("row {}", id);
        want_contains(&what, &r, &v(16)?)?;
        want_contains(&what, &r, &v(19)?)?;
        let rad = r
            .rad()
            .to_exact_num()
            .map_err(|e| format!("{}: rad {:?}: {:?}", what, r.rad(), e))?;
        let bound = v(22)?.mul_full_prec(&num(4));
        if rad > bound {
            return Err(format!(
                "{}: TIGHT: rad {:?} > 4*RREF {:?}",
                what,
                r.rad(),
                bound
            ));
        }
    }
    Ok(())
}

/// The hypot radius covers both halves of `[lo, hi]`: the `ToEven` midpoint can land on either side of the centre.
/// `x = (1 + k·2^-64) ± 2^-64` with `y = 0` at prec 64 gives `[|x| - r, |x| + r]` exactly; every k parity is checked.
/// A narrow non-straddling pair stays at the monotone width (below the 1-Lipschitz `r + s`).
#[test]
fn hypot_radius_covers_both_endpoints() -> Result<(), String> {
    let ulp = || pow2(-64);
    for k in 1..=8i64 {
        let x = || -> Result<ExactNum, String> {
            Ok(num(1).add_full_prec(&num(k).mul_full_prec(&ulp()?)))
        };
        for neg in [false, true] {
            let xm = if neg { -x()? } else { x()? };
            let r = ok("hypot", ball(xm, mag(1, 0)?, 64)?.hypot(&pt(0, 64)?, 64))?;
            let what = format!("k {} neg {}", k, neg);
            want_contains(&what, &r, &x()?.sub_full_prec(&ulp()?))?;
            want_contains(&what, &r, &x()?.add_full_prec(&ulp()?))?;
        }
    }
    // [2 ± 1, 3 ± 1]: lo = hypot(1, 2) = √5, hi = hypot(3, 4) = 5 exactly; the width (5 - √5)/2 < 1.39 < r + s = 2.
    for (a, b) in [(2i64, 3i64), (-2, 3), (2, -3), (-2, -3)] {
        let r = ok(
            "hypot",
            ball(num(a), mag(1, 64)?, 64)?.hypot(&ball(num(b), mag(1, 64)?, 64)?, 64),
        )?;
        want_contains("hi = 5", &r, &num(5))?;
        if r.rad() > &mag(0x58f5_c28f_5c28_f5c3, 2)? {
            return Err(format!("({}, {}): rad {:?} above 1.39", a, b, r.rad()));
        }
    }
    // [4 ± 1, 5 ± 1]: lo = hypot(3, 4) = 5 exactly.
    for (a, b) in [(4i64, 5i64), (-4, -5)] {
        let r = ok(
            "hypot",
            ball(num(a), mag(1, 64)?, 64)?.hypot(&ball(num(b), mag(1, 64)?, 64)?, 64),
        )?;
        want_contains("lo = 5", &r, &num(5))?;
    }
    Ok(())
}

/// powi carries guard bits through its squaring chain, so an exact input's power is within one ulp at `prec`
/// (after Harper's powi_exact_input_within_one_ulp; kills 246-POWIGUARD).
#[test]
fn powi_exact_input_within_one_ulp() -> Result<(), String> {
    let mut checked = 0u32;
    for (i, prec) in [53usize, 64, 113].into_iter().enumerate() {
        for j in 0..12i64 {
            let x = || -> Result<ExactNum, String> {
                Ok(num(1).add_full_prec(
                    &num(j
                        .checked_mul(37)
                        .and_then(|v| v.checked_add(5))
                        .ok_or("x")?)
                    .mul_full_prec(&pow2(-20)?),
                ))
            };
            let n = 100u64
                .checked_add(
                    u64::try_from(j)
                        .map_err(|_| "j")?
                        .checked_mul(71)
                        .ok_or("n")?,
                )
                .and_then(|v| v.checked_add(u64::try_from(i).ok()?))
                .ok_or("n")?;
            let mut want = num(1);
            for _ in 0..n {
                want = want.mul_full_prec(&x()?);
            }
            let res = ok("powi", ok("exact", Ball::exact(x()?, prec))?.powi(n, prec))?;
            if !res.is_informative() {
                continue;
            }
            want_contains("powi", &res, &want)?;
            let me = res.mid().exponent().ok_or("exp")?;
            let ulp = Mag::half_ulp(me, prec.checked_sub(1).ok_or("p-1")?);
            if res.rad() > &ulp {
                return Err(format!(
                    "prec {} n {}: rad {:?} > one ulp {:?}",
                    prec,
                    n,
                    res.rad(),
                    ulp
                ));
            }
            checked = checked.checked_add(1).ok_or("c")?;
        }
    }
    if checked < 30 {
        return Err(format!("only {} informative powi rows", checked));
    }
    Ok(())
}

// ---- abs, mul_2exp, Tri comparisons (CONTRACT 15.1, 15.2, O14-W, 16.S1.10, 16.S1.11) -------------------------------

fn want_tri(what: &str, got: Tri, want: Tri) -> Result<(), String> {
    if got != want {
        return Err(format!("{}: got {:?}, want {:?}", what, got, want));
    }
    Ok(())
}

#[test]
fn abs_rows() -> Result<(), String> {
    let b = ok("abs -3", pt(-3, 64)?.abs(64))?;
    want_exact("abs -3", &b, &num(3), &Mag::zero())?;
    let b = ok("abs 3", pt(3, 64)?.abs(64))?;
    want_exact("abs 3", &b, &num(3), &Mag::zero())?;
    let b = ok("abs -3±1", ball(num(-3), mag(TOP, 1)?, 64)?.abs(64))?;
    want_exact("abs -3±1", &b, &num(3), &mag(TOP, 1)?)?;
    // Touching 0: [-2, 0] maps to [0, 2] = 1 ± 1, exactly.
    let b = ok("abs -1±1", ball(num(-1), mag(TOP, 1)?, 64)?.abs(64))?;
    want_exact("abs -1±1", &b, &num(1), &mag(TOP, 1)?)?;
    // Straddling: [-3, 1] maps to [0, 3] = 1.5 ± 1.5.
    let b = ok("abs -1±2", ball(num(-1), mag(TOP, 2)?, 64)?.abs(64))?;
    want_exact("abs -1±2", &b, &dy(3, -1)?, &mag(0xC000_0000_0000_0000, 1)?)?;
    let b = ok("abs 1±2", ball(num(1), mag(TOP, 2)?, 64)?.abs(64))?;
    want_exact("abs 1±2", &b, &dy(3, -1)?, &mag(0xC000_0000_0000_0000, 1)?)?;
    let b = ok("abs 0±1", ball(num(0), mag(TOP, 1)?, 64)?.abs(64))?;
    want_exact("abs 0±1", &b, &dy(1, -1)?, &mag(TOP, 0)?)?;
    let b = ok("abs -0", ok("-0", Ball::exact(-num(0), 64))?.abs(64))?;
    if b.mid().is_negative() || !b.is_exact() {
        return Err(format!("abs -0: got {:?}", b));
    }
    want_all_reals("abs all reals", ok("ar", Ball::all_reals(64))?.abs(64))?;
    want_invalid("abs prec 0", pt(-3, 64)?.abs(0))?;
    // The midpoint is rounded to the requested prec and the radius absorbs it.
    let x = || num(-1).sub_full_prec(&pow2(64).unwrap_or_else(|_| num(0)));
    let b = ok("abs round", ok("x", Ball::exact(x(), 128))?.abs(2))?;
    want_contains("abs round", &b, &x().abs())?;
    if b.is_exact() {
        return Err(format!("abs round: want a rounded midpoint, got {:?}", b));
    }
    // Endpoints and midpoint of random balls map inside.
    let mut rng = SplitMix64::new(0x5b1_ab5);
    for &e in EXPS.iter().cycle().take(200) {
        let seed = rng.next_u64();
        let mk = || near_ball(&mut SplitMix64::new(seed), e, 64);
        let got = ok("abs rand", mk()?.abs(64))?;
        let b0 = mk()?;
        want_contains("abs rand mid", &got, &absn(b0.mid()))?;
        if let Ok(rn) = b0.rad().to_exact_num() {
            want_contains("abs rand hi", &got, &b0.mid().add_full_prec(&rn).abs())?;
            want_contains("abs rand lo", &got, &b0.mid().sub_full_prec(&rn).abs())?;
        }
    }
    Ok(())
}

#[test]
fn mul_2exp_rows() -> Result<(), String> {
    let three = || ball(num(3), mag(TOP, 1).unwrap_or_else(|_| Mag::zero()), 64);
    let b = ok("3±1 << 5", three()?.mul_2exp(5, 64))?;
    want_exact("3±1 << 5", &b, &num(96), &mag(TOP, 6)?)?;
    let b = ok("3±1 >> 5", three()?.mul_2exp(-5, 64))?;
    want_exact("3±1 >> 5", &b, &dy(3, -5)?, &mag(TOP, -4)?)?;
    let b = ok("3±1 << 0", three()?.mul_2exp(0, 64))?;
    want_exact("3±1 << 0", &b, &num(3), &mag(TOP, 1)?)?;
    // i64 shifts beyond zenith-float's range: overflow is all reals, underflow is enclosed (never a silent flush).
    want_all_reals("1 << MAX", pt(1, 64)?.mul_2exp(EXPONENT_MAX, 64))?;
    want_all_reals("3±1 << i64::MAX", three()?.mul_2exp(i64::MAX, 64))?;
    let k = EXPONENT_MIN.checked_sub(1000).ok_or("k")?;
    let b = ok("1 >> past MIN", pt(1, 64)?.mul_2exp(k, 64))?;
    if b.is_exact() || b.contains(&num(0)) != Tri::True {
        return Err(format!(
            "1 >> past MIN: want an inexact ball around 0, got {:?}",
            b
        ));
    }
    // The radius alone: exact while the Mag exponent fits, the smallest Mag below that.
    let b = ok(
        "0±1 >> 2^40",
        ball(num(0), mag(TOP, 1)?, 64)?.mul_2exp(-(1 << 40), 64),
    )?;
    want_exact("0±1 >> 2^40", &b, &num(0), &mag(TOP, 1 - (1 << 40))?)?;
    let b = ok(
        "0±1 >> MIN",
        ball(num(0), mag(TOP, 1)?, 64)?.mul_2exp(i64::MIN, 64),
    )?;
    want_exact("0±1 >> MIN", &b, &num(0), &mag(TOP, i64::MIN + 1)?)?;
    let b = ok(
        "0±2^-6 >> MIN",
        ball(num(0), mag(TOP, -5)?, 64)?.mul_2exp(i64::MIN, 64),
    )?;
    want_exact("0±2^-6 >> MIN", &b, &num(0), &mag(TOP, i64::MIN)?)?;
    let b = ok(
        "0±1 << MAX",
        ball(num(0), mag(TOP, 1)?, 64)?.mul_2exp(i64::MAX, 64),
    )?;
    if b.is_informative() {
        return Err(format!("0±1 << MAX: want Infinity radius, got {:?}", b));
    }
    let b = ok(
        "0±max << 1",
        ball(num(0), mag(TOP, i64::MAX)?, 64)?.mul_2exp(1, 64),
    )?;
    if b.is_informative() {
        return Err(format!("0±max << 1: want Infinity radius, got {:?}", b));
    }
    let b = ok("0 >> 1", pt(0, 64)?.mul_2exp(-1, 64))?;
    want_exact("0 >> 1", &b, &num(0), &Mag::zero())?;
    want_all_reals(
        "all reals << 3",
        ok("ar", Ball::all_reals(64))?.mul_2exp(3, 64),
    )?;
    want_invalid("mul_2exp prec 0", pt(3, 64)?.mul_2exp(1, 0))?;
    // The midpoint is rounded to the requested prec.
    let x = || num(1).add_full_prec(&pow2(64).unwrap_or_else(|_| num(0)));
    let b = ok("round", ok("x", Ball::exact(x(), 128))?.mul_2exp(1, 2))?;
    want_contains("round", &b, &x().mul_full_prec(&num(2)))?;
    Ok(())
}

/// `lo ± r` and `hi ± r` style helper: `v ± 2^(k)` with `k` the Mag exponent minus 1.
fn bl(v: i64, rad_exp: i64) -> Result<Ball, String> {
    ball(num(v), mag(TOP, rad_exp)?, 64)
}

/// R > 0 and lo(a) = hi(b): every point of `a` is >= every point of `b`, so `a < b` is certainly False (review nit a).
#[test]
fn lt_touching_from_above_is_false() -> Result<(), String> {
    for (name, a, b) in [
        ("[2,4] < [0,2]", bl(3, 1)?, bl(1, 1)?),
        ("2 < [0,2]", pt(2, 64)?, bl(1, 1)?),
        ("[0,2] < 0", bl(1, 1)?, pt(0, 64)?),
    ] {
        want_tri(&format!("{} lt", name), a.lt(&b), Tri::False)?;
        want_tri(&format!("{} gt", name), b.gt(&a), Tri::False)?;
    }
    Ok(())
}

#[test]
fn tri_compare_rows() -> Result<(), String> {
    use Tri::{False as F, True as T, Unknown as U};
    // (name, a, b, lt, le, gt, ge, eq, ne, overlaps)
    let rows: [(&str, Ball, Ball, [Tri; 7]); 11] = [
        (
            "touching [0,2] [2,4]",
            bl(1, 1)?,
            bl(3, 1)?,
            [U, T, F, U, U, U, T],
        ),
        (
            "touching [2,4] [0,2]",
            bl(3, 1)?,
            bl(1, 1)?,
            [F, U, U, T, U, U, T],
        ),
        (
            "disjoint [0,2] [4,6]",
            bl(1, 1)?,
            bl(5, 1)?,
            [T, T, F, F, F, T, F],
        ),
        (
            "disjoint [4,6] [0,2]",
            bl(5, 1)?,
            bl(1, 1)?,
            [F, F, T, T, F, T, F],
        ),
        (
            "exact 3 vs 3",
            pt(3, 64)?,
            pt(3, 64)?,
            [F, T, F, T, T, F, T],
        ),
        (
            "exact 3 vs 4",
            pt(3, 64)?,
            pt(4, 64)?,
            [T, T, F, F, F, T, F],
        ),
        (
            "[0,2] vs edge point 2",
            bl(1, 1)?,
            pt(2, 64)?,
            [U, T, F, U, U, U, T],
        ),
        (
            "edge point 0 vs [0,2]",
            pt(0, 64)?,
            bl(1, 1)?,
            [U, T, F, U, U, U, T],
        ),
        (
            "interior point 1",
            bl(1, 1)?,
            pt(1, 64)?,
            [U, U, U, U, U, U, T],
        ),
        (
            "nested [0,2] [0.5,1.5]",
            bl(1, 1)?,
            bl(1, 0)?,
            [U, U, U, U, U, U, T],
        ),
        (
            "prec ignored",
            pt(3, 2)?,
            pt(3, 4096)?,
            [F, T, F, T, T, F, T],
        ),
    ];
    for (name, a, b, w) in rows.iter() {
        let [lt, le, gt, ge, eq, ne, ov] = w;
        want_tri(&format!("{} lt", name), a.lt(b), clone_tri(lt))?;
        want_tri(&format!("{} le", name), a.le(b), clone_tri(le))?;
        want_tri(&format!("{} gt", name), a.gt(b), clone_tri(gt))?;
        want_tri(&format!("{} ge", name), a.ge(b), clone_tri(ge))?;
        want_tri(&format!("{} eq", name), a.eq(b), clone_tri(eq))?;
        want_tri(&format!("{} ne", name), a.ne(b), clone_tri(ne))?;
        want_tri(&format!("{} overlaps", name), a.overlaps(b), clone_tri(ov))?;
    }
    // Any comparison with an Infinity radius is Unknown, even a ball against itself; it overlaps everything.
    let ar = ok("ar", Ball::all_reals(64))?;
    let p = pt(3, 64)?;
    for (n, a, b) in [("ar ar", &ar, &ar), ("ar p", &ar, &p), ("p ar", &p, &ar)] {
        for (op, t) in [
            ("lt", a.lt(b)),
            ("le", a.le(b)),
            ("gt", a.gt(b)),
            ("ge", a.ge(b)),
            ("eq", a.eq(b)),
            ("ne", a.ne(b)),
        ] {
            want_tri(&format!("{} {}", n, op), t, U)?;
        }
        want_tri(&format!("{} overlaps", n), a.overlaps(b), T)?;
    }
    // A gap far below every radius and midpoint is still decided exactly.
    let tiny = ok("tiny", Ball::exact(pow2(-(1 << 20))?, 64))?;
    want_tri("0 lt tiny", pt(0, 64)?.lt(&tiny), T)?;
    want_tri("tiny gt 0", tiny.gt(&pt(0, 64)?), T)?;
    let near = ball(pow2(-(1 << 20))?, mag(TOP, 1 - (1 << 20))?, 64)?;
    want_tri("0 le [0, 2^-2^20+1]", pt(0, 64)?.le(&near), T)?;
    want_tri("0 lt [0, 2^-2^20+1]", pt(0, 64)?.lt(&near), U)?;
    Ok(())
}

fn clone_tri(t: &Tri) -> Tri {
    match t {
        Tri::True => Tri::True,
        Tri::False => Tri::False,
        Tri::Unknown => Tri::Unknown,
    }
}

fn absn(x: &ExactNum) -> ExactNum {
    if x.is_negative() {
        num(0).sub_full_prec(x)
    } else {
        num(0).add_full_prec(x)
    }
}

#[test]
fn compare_rad_outside_exact_range_is_unknown() -> Result<(), String> {
    // A finite radius whose exponent is outside the ExactNum range: since S2, cmp_exact_num (16.S2.1) decides that
    // 0 ± 2^huge reaches 1, so overlaps is True; the orders stay Unknown.
    let wide = ball(num(0), mag(TOP, i64::MAX)?, 64)?;
    let p = pt(1, 64)?;
    want_tri("overlaps wide", wide.overlaps(&p), Tri::True)?;
    want_tri("overlaps wide rev", p.overlaps(&wide), Tri::True)?;
    want_tri("lt wide", wide.lt(&p), Tri::Unknown)?;
    want_tri("eq wide", wide.eq(&p), Tri::Unknown)
}

/// 16.S1.7: under `Clip` a ball with a negative midpoint reaching above 0 has image `[0, sqrt(m + r)]`, and `m + r` is
/// `r - |m|`, not `|m| + r`. `-1 ± 3` gives `[0, sqrt 2]`, so 3/2 is outside (kills 286-SQRTNEGHI, which used `|m| + r`).
#[test]
fn sqrt_clip_negative_mid_hull_is_m_plus_r() -> Result<(), String> {
    let b = ok(
        "clip -1±3",
        ball(num(-1), mag(0xC000_0000_0000_0000, 2)?, 53)?.sqrt(53, &DomainPolicy::Clip),
    )?;
    want_contains("clip 0", &b, &num(0))?;
    want_contains("clip sqrt2", &b, &num(2).sqrt(256, RoundingMode::Up))?;
    match b.contains(&dy(3, -1)?) {
        Tri::False => Ok(()),
        t @ (Tri::True | Tri::Unknown) => Err(format!("clip -1±3: {:?} reaches 3/2 ({:?})", b, t)),
    }
}

/// 16.S1.7: for `m - r > 0` the propagated radius is `r / sqrt(m)↓` rounded up, so it is at least `r / sqrt(m)`:
/// `rad^2 · m >= r^2` exactly (at 4096 bits ε is far below the `Mag` rounding). Kills 287-SQRTROOTUP (`sqrt_up`).
#[test]
fn sqrt_radius_is_at_least_r_over_sqrt_m() -> Result<(), String> {
    let r = mag(TOP, -9)?;
    let rn = r.to_exact_num().map_err(|e| format!("r: {:?}", e))?;
    let r2 = rn.mul_full_prec(&rn);
    for v in 2i64..=200 {
        let b = ok(
            "sqrt",
            ball(num(v), mag(TOP, -9)?, 4096)?.sqrt(4096, &DomainPolicy::Strict),
        )?;
        let rad = b
            .rad()
            .to_exact_num()
            .map_err(|e| format!("rad {}: {:?}", v, e))?;
        let lhs = rad.mul_full_prec(&rad).mul_full_prec(&num(v));
        match lhs.partial_cmp(&r2) {
            Some(core::cmp::Ordering::Greater | core::cmp::Ordering::Equal) => {}
            other => {
                return Err(format!(
                    "sqrt({} ± 2^-10): rad {:?} below r/sqrt(m) ({:?})",
                    v,
                    b.rad(),
                    other
                ))
            }
        }
    }
    Ok(())
}

/// CONTRACT 10.4: only an INEXACT result in the top binade counts as overflow. An exact midpoint with exponent
/// `EXPONENT_MAX` stays exact (kills 290-OVFEXACT, which dropped the `inexact()` test).
#[test]
fn exact_top_binade_midpoint_is_not_overflow() -> Result<(), String> {
    let x = ExactNum::from_mag_parts(TOP, EXPONENT_MAX, 64).map_err(|e| format!("top: {:?}", e))?;
    let top = || -> Result<Ball, String> { ok("top", Ball::exact(pow2(EXPONENT_MAX - 1)?, 53)) };
    want_exact(
        "top + 0",
        &ok("add", top()?.add(&pt(0, 53)?, 53))?,
        &x,
        &Mag::zero(),
    )?;
    want_exact(
        "top * 1",
        &ok("mul", top()?.mul(&pt(1, 53)?, 53))?,
        &x,
        &Mag::zero(),
    )?;
    want_exact(
        "2^(EMAX-2) * 2",
        &ok(
            "mul_2exp",
            ok("half", Ball::exact(pow2(EXPONENT_MAX - 2)?, 53))?.mul_2exp(1, 53),
        )?,
        &x,
        &Mag::zero(),
    )
}

// ---- Harper S1 coverage rows (child module under tests/harness/slice_arith/) ----------------------------------------

mod harper_s1cov;

/// The borrowing S1 ops (CONTRACT 15.1, O14-W for abs): one ball feeds several ops and is still usable afterwards,
/// as a complex ball's parts must be (16.S1.5-16.S1.8, 16.S1.10). On the by-value API this does not compile.
#[test]
fn borrowed_ops_reuse_the_same_ball() -> Result<(), String> {
    let p = 64;
    let a = pt(4, p)?;
    let b = pt(3, p)?;
    want_contains("sqr", &ok("sqr", a.sqr(p))?, &num(16))?;
    want_contains("inv", &ok("inv", a.inv(p))?, &dy(1, -2)?)?;
    want_contains(
        "sqrt",
        &ok("sqrt", a.sqrt(p, &DomainPolicy::Strict))?,
        &num(2),
    )?;
    want_contains("hypot", &ok("hypot", a.hypot(&b, p))?, &num(5))?;
    let hh = ok("hypot self", b.hypot(&b, p))?;
    want_contains("hypot self", &ok("sqr", hh.sqr(p))?, &num(18))?;
    want_contains("powi", &ok("powi", a.powi(3, p))?, &num(64))?;
    want_contains("abs", &ok("abs", (-pt(4, p)?).abs(p))?, &num(4))?;
    let n = -pt(4, p)?;
    want_contains("abs twice", &ok("abs", n.abs(p))?, &num(4))?;
    want_contains("abs twice", &ok("abs", n.abs(p))?, &num(4))?;
    // `a` and `b` are still owned here: the by-value ops consume them last.
    want_contains("add after", &ok("add", a.add(&b, p))?, &num(7))
}
