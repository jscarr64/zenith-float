//! Slice CB harness rows (`src/cball.rs`). Owned by slice CB alone.
//! CONTRACT 16.CB.1 to 16.CB.11. Ownership is O9 (16.CB.0): a CB op takes `self` and borrows its second operand, so every
//! use below builds its operands fresh. `mul` and `div` wait on borrowing zf ops (Rusty, after R9);
//! their tests (and sqr's width comparison against mul) are `#[ignore]`d red until then.
//! Oracles: hand-verified exact rows (worked in the comment beside each) and containment of exact images computed with
//! zenith-float's `add_full_prec`/`sub_full_prec`/`mul_full_prec` at sample points (corners and midpoints) of the inputs.
//! The mpmath.mpc golden rows and the RREF tightness factors of 16.CB.10: tests/goldens/slices/cball/cball_rows.txt
//! (scripts/gen_cball_goldens.py).
//! COVERS: CBall::new, CBall::from_real, CBall::all_complex, CBall::re, CBall::im, CBall::prec, CBall::identical, Neg for CBall, CBall::conj, CBall::add, CBall::sub, CBall::mul, CBall::sqr, CBall::abs, CBall::div

use std::fs;

use crate::common::goldens::golden_dir;
use zenith_bound::ball::Ball;
use zenith_bound::cball::CBall;
use zenith_bound::enclosure::Enclosure;
use zenith_bound::error::BoundError;
use zenith_bound::mag::Mag;
use zenith_bound::tri::Tri;
use zenith_float_num::ExactNum;

const TOP: u64 = 0x8000_0000_0000_0000;
const PRECS: [usize; 7] = [1, 2, 53, 64, 113, 256, 1024];

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

fn ok<T>(what: &str, r: Result<T, BoundError>) -> Result<T, String> {
    r.map_err(|e| format!("{}: unexpected {:?}", what, e))
}

fn ball(v: i64, rad: Mag, prec: usize) -> Result<Ball, String> {
    ok("from_mid_rad", Ball::from_mid_rad(num(v), rad, prec))
}

fn pt(v: i64, prec: usize) -> Result<Ball, String> {
    ok("exact", Ball::exact(num(v), prec))
}

/// The exact point `re + i·im`.
fn cpt(re: i64, im: i64, prec: usize) -> Result<CBall, String> {
    ok("CBall::new", CBall::new(pt(re, prec)?, pt(im, prec)?))
}

/// `(re ± 2^-k) + i(im ± 2^-k)`; Mag 2^-k is finite(TOP, 1 - k).
fn cbox(re: i64, im: i64, k: i64, prec: usize) -> Result<CBall, String> {
    let r = 1i64.checked_sub(k).ok_or("cbox k")?;
    ok(
        "CBall::new",
        CBall::new(ball(re, mag(TOP, r)?, prec)?, ball(im, mag(TOP, r)?, prec)?),
    )
}

fn want_point(what: &str, z: &CBall, re: i64, im: i64) -> Result<(), String> {
    if z.re().mid() != &num(re)
        || z.im().mid() != &num(im)
        || !z.re().rad().is_zero()
        || !z.im().rad().is_zero()
    {
        return Err(format!(
            "{}: got {:?}, want exact {} + {}i",
            what, z, re, im
        ));
    }
    Ok(())
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

fn want_all_reals_both(what: &str, z: &CBall) -> Result<(), String> {
    if z.re().is_informative()
        || z.im().is_informative()
        || !z.re().rad().is_infinity()
        || !z.im().rad().is_infinity()
    {
        return Err(format!("{}: want both parts all reals, got {:?}", what, z));
    }
    Ok(())
}

fn want_finite_both(what: &str, z: &CBall) -> Result<(), String> {
    if !z.re().is_informative() || !z.im().is_informative() {
        return Err(format!("{}: want both parts finite, got {:?}", what, z));
    }
    Ok(())
}

/// Sample points (corners and midpoint) of `m ± 2^-k` for integer `m`.
fn samples(m: i64, k: i64) -> Result<[ExactNum; 3], String> {
    let r = pow2(-k)?;
    let c = num(m);
    Ok([c.sub_full_prec(&r), num(m), c.add_full_prec(&r)])
}

// ---- 16.CB.1 construction and invariants ----------------------------------------------------------------------

#[test]
fn cball_construct_invariants() -> Result<(), String> {
    // I2: parts of different prec are InvalidArgument.
    match CBall::new(pt(1, 53)?, pt(2, 64)?) {
        Err(BoundError::InvalidArgument { .. }) => {}
        other => {
            return Err(format!(
                "prec mismatch: want InvalidArgument, got {:?}",
                other
            ))
        }
    }
    let z = ok("new", CBall::new(pt(3, 53)?, pt(-4, 53)?))?;
    if z.prec() != 53 {
        return Err(format!("prec: got {}", z.prec()));
    }
    want_point("new", &z, 3, -4)?;
    want_point(
        "from_real",
        &ok("from_real", CBall::from_real(pt(7, 53)?))?,
        7,
        0,
    )?;
    want_all_reals_both("all_complex", &ok("all_complex", CBall::all_complex(53))?)?;
    // I3: an all-reals part leaves the other informative.
    let h = ok(
        "new",
        CBall::new(ok("all_reals", Ball::all_reals(53))?, pt(1, 53)?),
    )?;
    if h.re().is_informative() || !h.im().is_informative() {
        return Err(format!("I3: got {:?}", h));
    }
    if !z.identical(&ok("new", CBall::new(pt(3, 53)?, pt(-4, 53)?))?) || z.identical(&h) {
        return Err("identical".into());
    }
    // identical compares both parts: same re, different im is not identical (kills 1521-IDENTIM).
    if z.identical(&cpt(3, 4, 53)?) || z.identical(&cpt(-3, -4, 53)?) {
        return Err("identical ignores a part".into());
    }
    Ok(())
}

// ---- 16.CB.3 neg, conj, add, sub -----------------------------------------------------------------------------

#[test]
fn cball_add_sub_neg_conj_rows() -> Result<(), String> {
    let a = || cpt(3, 5, 53);
    let b = cpt(-1, 2, 53)?;
    want_point("add", &ok("add", a()?.add(&b, 53))?, 2, 7)?; // (3+5i)+(-1+2i) = 2+7i
    want_point("sub", &ok("sub", a()?.sub(&b, 53))?, 4, 3)?; // (3+5i)-(-1+2i) = 4+3i
    want_point("neg", &(-a()?), -3, -5)?;
    want_point("conj", &a()?.conj(), 3, -5)?;
    if !a()?.conj().conj().identical(&a()?) {
        return Err("conj∘conj".into());
    }
    // conj keeps radii: (1 ± 2^-4) + (2 ± 2^-4)i → im radius unchanged, mid negated.
    let c = cbox(1, 2, 4, 53)?;
    let cc = cbox(1, 2, 4, 53)?.conj();
    if cc.im().mid() != &num(-2) || cc.im().rad() != c.im().rad() || !cc.re().identical(c.re()) {
        return Err(format!("conj radii: {:?}", cc));
    }
    // containment on samples of a box sum.
    let s = ok("add", c.add(&cbox(-3, 5, 6, 53)?, 53))?;
    for x in samples(1, 4)?.iter() {
        for y in samples(-3, 6)?.iter() {
            want_contains("add re", s.re(), &x.add_full_prec(y))?;
        }
    }
    Ok(())
}

// ---- 16.CB.4 mul ---------------------------------------------------------------------------------------------

#[test]
#[ignore = "red: not implemented: CONTRACT CB.4 (mul): waits zf borrowing ops"]
fn cball_mul_rows() -> Result<(), String> {
    // (1+2i)(3+4i) = (3-8) + (4+6)i = -5+10i.
    want_point(
        "mul",
        &ok("mul", cpt(1, 2, 53)?.mul(&cpt(3, 4, 53)?, 53))?,
        -5,
        10,
    )?;
    // cancellation: (1+i)(1-i) = 2 + 0i.
    want_point(
        "mul cancel",
        &ok("mul", cpt(1, 1, 53)?.mul(&cpt(1, -1, 53)?, 53))?,
        2,
        0,
    )?;
    // i·i = -1.
    want_point(
        "i*i",
        &ok("mul", cpt(0, 1, 53)?.mul(&cpt(0, 1, 53)?, 53))?,
        -1,
        0,
    )?;
    // boxes, parts containing 0 and mixed signs: containment at every sample corner pair.
    for &prec in PRECS.iter() {
        let (ar, ai, br, bi) = (0i64, -3i64, 5i64, 0i64);
        let p = ok(
            "mul",
            cbox(ar, ai, 2, prec)?.mul(&cbox(br, bi, 3, prec)?, prec),
        )?;
        for x in samples(ar, 2)?.iter() {
            for y in samples(ai, 2)?.iter() {
                for u in samples(br, 3)?.iter() {
                    for v in samples(bi, 3)?.iter() {
                        let re = x.mul_full_prec(u).sub_full_prec(&y.mul_full_prec(v));
                        let im = x.mul_full_prec(v).add_full_prec(&y.mul_full_prec(u));
                        want_contains("mul re", p.re(), &re)?;
                        want_contains("mul im", p.im(), &im)?;
                    }
                }
            }
        }
    }
    Ok(())
}

// ---- 16.CB.5 sqr ---------------------------------------------------------------------------------------------

#[test]
fn cball_sqr_rows() -> Result<(), String> {
    // (3+4i)^2 = 9-16 + 24i = -7+24i.
    want_point("sqr", &ok("sqr", cpt(3, 4, 53)?.sqr(53))?, -7, 24)
}

#[test]
fn cball_sqr_contains_images() -> Result<(), String> {
    // Every row and prec: sqr(a) contains (x+iy)^2 at the sample points; im is an exact doubling of re·im (16.CB.5).
    for &prec in PRECS.iter() {
        for &(re, im, k) in SQR_ROWS.iter() {
            let s = ok("sqr", cbox(re, im, k, prec)?.sqr(prec))?;
            want_finite_both("sqr", &s)?;
            if s.prec() != prec {
                return Err(format!("sqr prec {} != {}", s.prec(), prec));
            }
            for x in samples(re, k)?.iter() {
                for y in samples(im, k)?.iter() {
                    want_contains(
                        "sqr re",
                        s.re(),
                        &x.mul_full_prec(x).sub_full_prec(&y.mul_full_prec(y)),
                    )?;
                    let xy = x.mul_full_prec(y);
                    want_contains("sqr im", s.im(), &xy.add_full_prec(&xy))?;
                }
            }
        }
    }
    // An all-reals part spreads to both output parts; a prec mismatch or prec 0 is refused.
    let h = ok(
        "new",
        CBall::new(ok("all_reals", Ball::all_reals(53))?, pt(1, 53)?),
    )?;
    want_all_reals_both("sqr all-reals re", &ok("sqr", h.sqr(53))?)?;
    if !matches!(
        cpt(1, 2, 53)?.sqr(0),
        Err(BoundError::InvalidArgument { .. })
    ) {
        return Err("sqr prec 0 not InvalidArgument".into());
    }
    Ok(())
}

const SQR_ROWS: [(i64, i64, i64); 5] = [(0, 0, 1), (0, 3, 2), (1, -1, 4), (-7, 2, 10), (5, 0, 0)];

/// `a·a` in the 16.CB.4 four-product form, each product `Ball::mul` and each combination `Ball::sub`/`add`: what
/// `CBall::mul(a, a)` computes. Every operand is built fresh (O9).
fn four_product_sq(re: i64, im: i64, k: i64, prec: usize) -> Result<(Ball, Ball), String> {
    let e = 1i64.checked_sub(k).ok_or("k")?;
    let part = |v: i64| -> Result<Ball, String> { ball(v, mag(TOP, e)?, prec) };
    let rr = ok("mul", part(re)?.mul(&part(re)?, prec))?;
    let ii = ok("mul", part(im)?.mul(&part(im)?, prec))?;
    let ri = ok("mul", part(re)?.mul(&part(im)?, prec))?;
    let ir = ok("mul", part(im)?.mul(&part(re)?, prec))?;
    Ok((ok("sub", rr.sub(&ii, prec))?, ok("add", ri.add(&ir, prec))?))
}

#[test]
fn cball_sqr_not_wider_than_four_product() -> Result<(), String> {
    // Every row and prec: each part of sqr(a) is no wider than the four-product a·a of 16.CB.4 (16.CB.5).
    for &prec in PRECS.iter() {
        for &(re, im, k) in SQR_ROWS.iter() {
            let s = ok("sqr", cbox(re, im, k, prec)?.sqr(prec))?;
            let (mre, mim) = four_product_sq(re, im, k, prec)?;
            if s.re().rad() > mre.rad() || s.im().rad() > mim.rad() {
                return Err(format!(
                    "sqr wider than four-product at prec {} for ({}, {}, 2^-{}): {:?} vs {:?} {:?}",
                    prec, re, im, k, s, mre, mim
                ));
            }
        }
    }
    Ok(())
}

#[test]
#[ignore = "red: the width comparison waits CBall::mul (zf borrowing ops)"]
fn cball_sqr_not_wider_than_mul() -> Result<(), String> {
    // Every row: each part of sqr(a) is no wider than the same part of mul(a, a) (16.CB.5).
    for &prec in PRECS.iter() {
        for &(re, im, k) in [
            (0i64, 0i64, 1i64),
            (0, 3, 2),
            (1, -1, 4),
            (-7, 2, 10),
            (5, 0, 0),
        ]
        .iter()
        {
            let a = cbox(re, im, k, prec)?;
            let s = ok("sqr", cbox(re, im, k, prec)?.sqr(prec))?;
            let m = ok("mul", cbox(re, im, k, prec)?.mul(&a, prec))?;
            if s.re().rad() > m.re().rad() || s.im().rad() > m.im().rad() {
                return Err(format!(
                    "sqr wider than mul at prec {} for {:?}: {:?} vs {:?}",
                    prec, a, s, m
                ));
            }
            for x in samples(re, k)?.iter() {
                for y in samples(im, k)?.iter() {
                    want_contains(
                        "sqr re",
                        s.re(),
                        &x.mul_full_prec(x).sub_full_prec(&y.mul_full_prec(y)),
                    )?;
                    let xy = x.mul_full_prec(y);
                    want_contains("sqr im", s.im(), &xy.add_full_prec(&xy))?;
                }
            }
        }
    }
    Ok(())
}

// ---- 16.CB.6 abs ---------------------------------------------------------------------------------------------

#[test]
fn cball_abs_rows() -> Result<(), String> {
    // |3+4i| = 5 exactly (hypot of exact points is exact when the result fits).
    let r = ok("abs", cpt(3, 4, 53)?.abs(53))?;
    want_contains("|3+4i|", &r, &num(5))?;
    // lower endpoint ≥ 0: a box around 0 gives a ball whose mid ≥ rad.
    let z = ok("abs", cbox(0, 0, 3, 53)?.abs(53))?;
    if let Tri::True = z.lt(&ok("exact", Ball::exact(num(0), 53))?) {
        return Err(format!("abs below 0: {:?}", z));
    }
    want_contains("|0|", &z, &num(0))?;
    // axis: |x + 0i| contains |x|.
    want_contains("|-6|", &ok("abs", cpt(-6, 0, 53)?.abs(53))?, &num(6))
}

#[test]
fn cball_abs_i64_extremes() -> Result<(), String> {
    // 16.CB.6: CB uses S1 hypot as is. If this fails, it is a red test for S1's hypot range, not a CB fix.
    // (3·2^K) + (4·2^K)i has |z| = 5·2^K, representable at every K below.
    for &k in [
        1i64 << 40,
        -(1i64 << 40),
        (1i64 << 62) - (1i64 << 10),
        -((1i64 << 62) - (1i64 << 10)),
    ]
    .iter()
    {
        let s = pow2(k)?;
        let re = ok("exact", Ball::exact(num(3).mul_full_prec(&s), 64))?;
        let im = ok("exact", Ball::exact(num(4).mul_full_prec(&s), 64))?;
        let r = ok("abs", ok("new", CBall::new(re, im))?.abs(64))?;
        if !r.is_informative() {
            return Err(format!("abs at 2^{}: not informative: {:?}", k, r));
        }
        want_contains("abs extreme", &r, &num(5).mul_full_prec(&s))?;
    }
    Ok(())
}

// ---- 16.CB.7 div ---------------------------------------------------------------------------------------------

#[test]
#[ignore = "red: not implemented: CONTRACT CB.7 (div): waits zf borrowing ops"]
fn cball_div_rows() -> Result<(), String> {
    // (2+4i)/(1+i) = (2+4i)(1-i)/2 = (2+4 + (4-2)i)/2 = 3+i.
    let q = ok("div", cpt(2, 4, 53)?.div(&cpt(1, 1, 53)?, 53))?;
    want_contains("div re", q.re(), &num(3))?;
    want_contains("div im", q.im(), &num(1))?;
    // by an exact real: (6-9i)/3 = 2-3i.
    let q = ok("div", cpt(6, -9, 53)?.div(&cpt(3, 0, 53)?, 53))?;
    want_contains("div real re", q.re(), &num(2))?;
    want_contains("div real im", q.im(), &num(-3))?;
    // one axis of the divisor touches 0, the other is clear: finite, not all reals.
    let q = ok("div", cpt(1, 0, 53)?.div(&cbox(0, 4, 1, 53)?, 53))?;
    want_finite_both("one axis touches 0", &q)
}

#[test]
#[ignore = "red: not implemented: CONTRACT CB.7 (div all reals R1): waits zf borrowing ops"]
fn cball_div_rectangle_containing_zero_is_all_reals() -> Result<(), String> {
    // R1: both parts reach 0 → Ok with both parts all reals (never Domain).
    let exact_zero = cpt(0, 0, 53)?;
    want_all_reals_both("0/0", &ok("div", cpt(1, 1, 53)?.div(&exact_zero, 53))?)?;
    want_all_reals_both(
        "box at 0",
        &ok("div", cpt(1, 1, 53)?.div(&cbox(0, 0, 3, 53)?, 53))?,
    )?;
    // edge |m| = r on both axes: 1 ± 1 touches 0 exactly.
    let edge = ok(
        "new",
        CBall::new(ball(1, mag(TOP, 1)?, 53)?, ball(-1, mag(TOP, 1)?, 53)?),
    )?;
    want_all_reals_both("edge", &ok("div", cpt(1, 1, 53)?.div(&edge, 53))?)?;
    // just clear on one axis: 2 ± 1 excludes 0, so R1 does not hold.
    let clear = ok(
        "new",
        CBall::new(ball(2, mag(TOP, 1)?, 53)?, ball(-1, mag(TOP, 1)?, 53)?),
    )?;
    want_finite_both("clear on re", &ok("div", cpt(1, 1, 53)?.div(&clear, 53))?)
}

/// `b = (2^-40 ± 2^-41) + 0i` at `prec`: re excludes 0 (|m| > r) and the exact |b|^2 is [2^-82, 9·2^-82], clear of 0.
fn r2_divisor_re(prec: usize) -> Result<Ball, String> {
    ok(
        "from_mid_rad",
        Ball::from_mid_rad(pow2(-40)?, mag(TOP, -40)?, prec),
    )
}

#[test]
fn cball_r2_row_s1_sqr_reaches_zero() -> Result<(), String> {
    // Measured (16.CB.7 R2 note): S1 sqr is mid-rad (16.S1.5), rad = 2|m|r + r^2 + ε = 2^-80 + 2^-82 + ε > m^2 = 2^-80,
    // so the enclosure of (2^-40 ± 2^-41)^2 reaches 0 at every prec. This is the observed S1 behavior the R2 row relies on.
    for &prec in PRECS.iter() {
        let d = ok("sqr", r2_divisor_re(prec)?.sqr(prec))?;
        want_contains("S1 sqr of 2^-40 ± 2^-41 reaches 0", &d, &num(0))?;
    }
    Ok(())
}

#[test]
#[ignore = "red: not implemented: CONTRACT CB.7 (div all reals R2): waits zf borrowing ops"]
fn cball_div_mid_rad_sqr_widening_is_all_reals() -> Result<(), String> {
    // R2 (16.CB.7), kept as a defensive check: R1 does not hold, but d = sqr(b.re) + sqr(b.im) reaches 0 because S1's
    // sqr is mid-rad (see cball_r2_row_s1_sqr_reaches_zero). Observed behavior: both parts all reals, at every prec.
    for &prec in PRECS.iter() {
        let b = ok("new", CBall::new(r2_divisor_re(prec)?, pt(0, prec)?))?;
        want_all_reals_both(
            "R2 mid-rad row",
            &ok("div", cpt(1, 0, prec)?.div(&b, prec))?,
        )?;
    }
    Ok(())
}

#[test]
#[ignore = "red: not implemented: CONTRACT CB.7 (div exact scaling): waits zf borrowing ops"]
fn cball_div_scaling_exact() -> Result<(), String> {
    // The common 2^-k on both operands leaves the quotient unchanged: (2+4i)·2^K / ((1+i)·2^K) contains 3+i at
    // K = ±2^40, where |b|^2 unscaled would leave the midpoint range.
    for &k in [1i64 << 40, -(1i64 << 40)].iter() {
        let s = pow2(k)?;
        let mk = |v: i64| ok("exact", Ball::exact(num(v).mul_full_prec(&s), 64));
        let a = ok("new", CBall::new(mk(2)?, mk(4)?))?;
        let b = ok("new", CBall::new(mk(1)?, mk(1)?))?;
        let q = ok("div", a.div(&b, 64))?;
        want_finite_both("scaled div", &q)?;
        want_contains("scaled re", q.re(), &num(3))?;
        want_contains("scaled im", q.im(), &num(1))?;
    }
    Ok(())
}

// ---- 16.CB.9 all-reals propagation and precision ----------------------------------------------------------------

#[test]
fn cball_all_reals_propagation() -> Result<(), String> {
    let h = || {
        ok(
            "new",
            CBall::new(ok("all_reals", Ball::all_reals(53))?, pt(1, 53)?),
        )
    };
    // add, sub, neg, conj keep I3: im stays finite.
    for (what, s) in [
        ("add", ok("add", h()?.add(&cpt(1, 1, 53)?, 53))?),
        ("sub", ok("sub", h()?.sub(&cpt(1, 1, 53)?, 53))?),
        ("neg", -h()?),
        ("conj", h()?.conj()),
    ] {
        if s.re().is_informative() || !s.im().is_informative() {
            return Err(format!("{} I3: {:?}", what, s));
        }
    }
    // abs of a rectangle with an all-reals part is not informative.
    let r = ok("abs", h()?.abs(53))?;
    if r.is_informative() {
        return Err(format!("abs of all-reals part: {:?}", r));
    }
    Ok(())
}

#[test]
#[ignore = "red: not implemented: CONTRACT CB.9 (all-reals propagation): waits zf borrowing ops"]
fn cball_all_reals_propagation_mul() -> Result<(), String> {
    // mul: both output parts use the all-reals re part.
    let h = ok(
        "new",
        CBall::new(ok("all_reals", Ball::all_reals(53))?, pt(1, 53)?),
    )?;
    want_all_reals_both("mul", &ok("mul", h.mul(&cpt(1, 1, 53)?, 53))?)
}

#[test]
fn cball_prec_sweep() -> Result<(), String> {
    // Both parts carry the op's prec; a 3-bit-or-wider midpoint rounds to prec with the error in rad.
    let a = || cbox(7, -5, 3, 1024);
    for &prec in PRECS.iter() {
        let b = cbox(3, 2, 5, 1024)?;
        for z in [
            ok("add", a()?.add(&b, prec))?,
            ok("sub", a()?.sub(&b, prec))?,
        ]
        .iter()
        {
            if z.prec() != prec || z.re().prec() != prec || z.im().prec() != prec {
                return Err(format!("prec {}: got {:?}", prec, z));
            }
        }
        let r = ok("abs", a()?.abs(prec))?;
        if r.prec() != prec {
            return Err(format!("abs prec {}: got {:?}", prec, r));
        }
    }
    match cpt(1, 1, 53)?.add(&cpt(1, 1, 53)?, 0) {
        Err(BoundError::InvalidArgument { .. }) => Ok(()),
        other => Err(format!("prec 0: want InvalidArgument, got {:?}", other)),
    }
}

#[test]
#[ignore = "red: not implemented: CONTRACT CB.9 (one prec per call): waits zf borrowing ops"]
fn cball_prec_sweep_mul_sqr_div() -> Result<(), String> {
    let a = || cbox(7, -5, 3, 1024);
    for &prec in PRECS.iter() {
        let b = cbox(3, 2, 5, 1024)?;
        for z in [
            ok("mul", a()?.mul(&b, prec))?,
            ok("sqr", a()?.sqr(prec))?,
            ok("div", a()?.div(&b, prec))?,
        ]
        .iter()
        {
            if z.prec() != prec || z.re().prec() != prec || z.im().prec() != prec {
                return Err(format!("prec {}: got {:?}", prec, z));
            }
        }
    }
    Ok(())
}

// ---- 16.CB.10 gold rows (mpmath.mpc point references) and tightness factors --------------------------------------

/// `SIGN EXP MANT` (scripts/gen_cball_goldens.py): value = (-1)^SIGN · MANT · 2^(EXP − bitlen(MANT)), MANT hex.
fn gold_val(sign: &str, exp: &str, hex: &str) -> Result<ExactNum, String> {
    let e: i64 = exp.parse().map_err(|_| format!("bad exponent {:?}", exp))?;
    let first = hex
        .chars()
        .next()
        .ok_or("empty mantissa")?
        .to_digit(16)
        .ok_or("bad digit")?;
    if first == 0 {
        return if hex == "0" && sign == "0" && e == 0 {
            Ok(num(0))
        } else {
            Err(format!("non-canonical zero {} {} {}", sign, exp, hex))
        };
    }
    let len = i64::try_from(hex.len()).map_err(|_| "mantissa length")?;
    let bitlen = 4 * (len - 1) + i64::from(32 - first.leading_zeros());
    let base = e - bitlen;
    let mut acc = num(0);
    let digits = hex.as_bytes();
    let mut shift = 0i64;
    // Least significant 16 hex digits (64 bits) first.
    for bytes in digits.rchunks(16) {
        let chunk = std::str::from_utf8(bytes).map_err(|_| "utf8")?;
        let c = u64::from_str_radix(chunk, 16).map_err(|_| format!("bad hex {:?}", chunk))?;
        if c != 0 {
            let lz = i64::from(c.leading_zeros());
            let part = ExactNum::from_mag_parts(c << lz, base + shift + 64 - lz, 64)
                .map_err(|e| format!("from_mag_parts: {:?}", e))?;
            acc = acc.add_full_prec(&part);
        }
        shift += 64;
    }
    if sign == "0" {
        Ok(acc)
    } else if sign == "1" {
        Ok(num(0).sub_full_prec(&acc))
    } else {
        Err(format!("bad sign {:?}", sign))
    }
}

struct GoldRow {
    op: String,
    prec: usize,
    a: [ExactNum; 4],
    b: Option<[ExactNum; 4]>,
    rref: [ExactNum; 2],
    refs: Vec<(ExactNum, ExactNum)>,
}

fn take_val(f: &[&str], i: &mut usize) -> Result<ExactNum, String> {
    let v = gold_val(
        f.get(*i).ok_or("short row")?,
        f.get(*i + 1).ok_or("short row")?,
        f.get(*i + 2).ok_or("short row")?,
    )?;
    *i += 3;
    Ok(v)
}

fn take4(f: &[&str], i: &mut usize) -> Result<[ExactNum; 4], String> {
    Ok([
        take_val(f, i)?,
        take_val(f, i)?,
        take_val(f, i)?,
        take_val(f, i)?,
    ])
}

fn gold_rows() -> Result<Vec<GoldRow>, String> {
    let path = golden_dir()
        .join("slices")
        .join("cball")
        .join("cball_rows.txt");
    let text = fs::read_to_string(&path).map_err(|e| format!("{}: {}", path.display(), e))?;
    let mut rows = Vec::new();
    for line in text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
    {
        let f: Vec<&str> = line.split(' ').collect();
        let op = String::from(*f.first().ok_or("empty row")?);
        let prec: usize = f
            .get(1)
            .ok_or("short row")?
            .parse()
            .map_err(|_| "bad prec")?;
        let mut i = 3;
        let a = take4(&f, &mut i)?;
        let b = if op == "sqr" || op == "abs" {
            None
        } else {
            Some(take4(&f, &mut i)?)
        };
        let rref = [take_val(&f, &mut i)?, take_val(&f, &mut i)?];
        let mut refs = Vec::new();
        for _ in 0..3 {
            refs.push((take_val(&f, &mut i)?, take_val(&f, &mut i)?));
        }
        if i != f.len() {
            return Err(format!("trailing fields in {:?}", line));
        }
        rows.push(GoldRow {
            op,
            prec,
            a,
            b,
            rref,
            refs,
        });
    }
    if rows.len() != 62 {
        return Err(format!("want 62 gold rows, got {}", rows.len()));
    }
    Ok(rows)
}

fn gold_box(p: &[ExactNum; 4], prec: usize) -> Result<CBall, String> {
    let re = ok(
        "from_mid_rad",
        Ball::from_mid_rad(
            p[0].add_full_prec(&num(0)),
            Mag::from_exact_num(&p[1]),
            prec,
        ),
    )?;
    let im = ok(
        "from_mid_rad",
        Ball::from_mid_rad(
            p[2].add_full_prec(&num(0)),
            Mag::from_exact_num(&p[3]),
            prec,
        ),
    )?;
    ok("CBall::new", CBall::new(re, im))
}

enum GoldOut {
    C(CBall),
    R(Ball),
}

fn gold_eval(r: &GoldRow) -> Result<GoldOut, String> {
    let a = gold_box(&r.a, r.prec)?;
    let b = match &r.b {
        Some(b) => Some(gold_box(b, r.prec)?),
        None => None,
    };
    let need_b = || b.as_ref().ok_or_else(|| format!("{}: no B", r.op));
    Ok(match r.op.as_str() {
        "add" => GoldOut::C(ok("add", a.add(need_b()?, r.prec))?),
        "sub" => GoldOut::C(ok("sub", a.sub(need_b()?, r.prec))?),
        "mul" => GoldOut::C(ok("mul", a.mul(need_b()?, r.prec))?),
        "div" => GoldOut::C(ok("div", a.div(need_b()?, r.prec))?),
        "sqr" => GoldOut::C(ok("sqr", a.sqr(r.prec))?),
        "abs" => GoldOut::R(ok("abs", a.abs(r.prec))?),
        other => return Err(format!("unknown op {:?}", other)),
    })
}

fn gold_contained(ops: &[&str]) -> Result<usize, String> {
    let mut n = 0;
    for r in gold_rows()?.iter().filter(|r| ops.contains(&r.op.as_str())) {
        let what = format!("{} p={}", r.op, r.prec);
        match gold_eval(r)? {
            GoldOut::C(z) => {
                for (re, im) in r.refs.iter() {
                    want_contains(&format!("{} re", what), z.re(), re)?;
                    want_contains(&format!("{} im", what), z.im(), im)?;
                }
            }
            GoldOut::R(b) => {
                for (re, _) in r.refs.iter() {
                    want_contains(&what, &b, re)?;
                }
            }
        }
        n += 1;
    }
    Ok(n)
}

/// 16.CB.10 factors as powers of two: add, sub, mul, sqr 2·RREF; abs 4·RREF; div 8·RREF (gold div rows have d clear by half).
fn factor_log2(op: &str) -> Result<i64, String> {
    match op {
        "add" | "sub" | "mul" | "sqr" => Ok(1),
        "abs" => Ok(2),
        "div" => Ok(3),
        other => Err(format!("no factor for {:?}", other)),
    }
}

fn want_tight(what: &str, b: &Ball, rref: &ExactNum, k: i64) -> Result<(), String> {
    let bound = Mag::from_exact_num(&rref.mul_full_prec(&pow2(k)?));
    if b.rad() > &bound {
        return Err(format!(
            "{}: rad {:?} > 2^{}·RREF {:?}",
            what,
            b.rad(),
            k,
            rref
        ));
    }
    Ok(())
}

fn gold_tight(ops: &[&str]) -> Result<usize, String> {
    let mut n = 0;
    for r in gold_rows()?.iter().filter(|r| ops.contains(&r.op.as_str())) {
        let what = format!("{} p={}", r.op, r.prec);
        let k = factor_log2(&r.op)?;
        match gold_eval(r)? {
            GoldOut::C(z) => {
                want_tight(&format!("{} re", what), z.re(), &r.rref[0], k)?;
                want_tight(&format!("{} im", what), z.im(), &r.rref[1], k)?;
            }
            GoldOut::R(b) => want_tight(&what, &b, &r.rref[0], k)?,
        }
        n += 1;
    }
    Ok(n)
}

#[test]
fn cball_golden_rows_contained() -> Result<(), String> {
    let n = gold_contained(&["add", "sub", "abs", "sqr"])?;
    if n != 42 {
        return Err(format!("want 42 add/sub/abs/sqr gold rows, got {}", n));
    }
    Ok(())
}

#[test]
fn cball_tightness_factors() -> Result<(), String> {
    let n = gold_tight(&["add", "sub", "abs"])?;
    if n != 30 {
        return Err(format!("want 30 add/sub/abs gold rows, got {}", n));
    }
    Ok(())
}

#[test]
#[ignore = "red: not implemented: CONTRACT CB.10 (golds): waits zf borrowing ops"]
fn cball_golden_rows_contained_mul_div() -> Result<(), String> {
    gold_contained(&["mul", "div"]).map(|_| ())
}

/// One part of a gold box as a fresh real ball: `i` 0 is re (`p[0] ± p[1]`), 1 is im (`p[2] ± p[3]`).
fn gold_part(p: &[ExactNum; 4], i: usize, prec: usize) -> Result<Ball, String> {
    let (mid, rad) = if i == 0 {
        (&p[0], &p[1])
    } else {
        (&p[2], &p[3])
    };
    ok(
        "from_mid_rad",
        Ball::from_mid_rad(mid.add_full_prec(&num(0)), Mag::from_exact_num(rad), prec),
    )
}

/// `a·a` for a gold sqr row in the 16.CB.4 four-product form (`Ball::mul` products, `Ball::sub`/`add`): what
/// `CBall::mul(a, a)` computes, while `CBall::mul` waits on zf borrowing ops. Every operand is built fresh (O9).
fn gold_four_product_sq(r: &GoldRow) -> Result<(Ball, Ball), String> {
    let part = |i: usize| gold_part(&r.a, i, r.prec);
    let rr = ok("mul", part(0)?.mul(&part(0)?, r.prec))?;
    let ii = ok("mul", part(1)?.mul(&part(1)?, r.prec))?;
    let ri = ok("mul", part(0)?.mul(&part(1)?, r.prec))?;
    let ir = ok("mul", part(1)?.mul(&part(0)?, r.prec))?;
    Ok((
        ok("sub", rr.sub(&ii, r.prec))?,
        ok("add", ri.add(&ir, r.prec))?,
    ))
}

#[test]
fn cball_sqr_gold_not_wider_than_four_product() -> Result<(), String> {
    // 16.CB.5 on the 12 gold sqr rows: each part of sqr(a) is no wider than mul(a, a), in the four-product form until
    // CBall::mul lands (cball_sqr_not_wider_than_mul keeps the direct comparison, ignored).
    let mut m = 0;
    for r in gold_rows()?.iter().filter(|r| r.op == "sqr") {
        let s = ok("sqr", gold_box(&r.a, r.prec)?.sqr(r.prec))?;
        let (mre, mim) = gold_four_product_sq(r)?;
        if s.re().rad() > mre.rad() || s.im().rad() > mim.rad() {
            return Err(format!(
                "sqr p={} wider than four-product a·a: {:?} vs {:?} {:?}",
                r.prec, s, mre, mim
            ));
        }
        m += 1;
    }
    if m != 12 {
        return Err(format!("want 12 sqr gold rows, got {}", m));
    }
    Ok(())
}

/// 16.CB.10 sqr reference per part, exact: `[RREF_re + rad(re)² + rad(im)², RREF_im + 2·rad(re)·rad(im)]`. RREF is
/// first-order, which cannot cover `re²` when `rad(re)` is not small; the second-order term comes from the gold row's
/// input radii (`a[1]`, `a[3]`) in exact arithmetic (Spence ruling, 2026-10-25).
fn sqr_ref(r: &GoldRow) -> [ExactNum; 2] {
    let (rx, ry) = (&r.a[1], &r.a[3]);
    let so_re = rx.mul_full_prec(rx).add_full_prec(&ry.mul_full_prec(ry));
    let xy = rx.mul_full_prec(ry);
    let so_im = xy.add_full_prec(&xy);
    [
        r.rref[0].add_full_prec(&so_re),
        r.rref[1].add_full_prec(&so_im),
    ]
}

#[test]
fn cball_tightness_factors_sqr() -> Result<(), String> {
    // 16.CB.10: each part of sqr(a) within 2·(RREF + second-order term) on all 12 gold sqr rows.
    let k = factor_log2("sqr")?;
    let mut n = 0;
    for r in gold_rows()?.iter().filter(|r| r.op == "sqr") {
        let what = format!("sqr p={}", r.prec);
        let s = ok("sqr", gold_box(&r.a, r.prec)?.sqr(r.prec))?;
        let [bre, bim] = sqr_ref(r);
        want_tight(&format!("{} re", what), s.re(), &bre, k)?;
        want_tight(&format!("{} im", what), s.im(), &bim, k)?;
        if std::env::var_os("CB10_MARGINS").is_some() {
            eprintln!(
                "CB10 {} row {} re rad {:?} bound {:?} im rad {:?} bound {:?}",
                what,
                n,
                s.re().rad(),
                Mag::from_exact_num(&bre.mul_full_prec(&pow2(k)?)),
                s.im().rad(),
                Mag::from_exact_num(&bim.mul_full_prec(&pow2(k)?))
            );
        }
        n += 1;
    }
    if n != 12 {
        return Err(format!("want 12 sqr gold rows, got {}", n));
    }
    Ok(())
}

#[test]
#[ignore = "red: not implemented: CONTRACT CB.10 (tightness): waits zf borrowing ops"]
fn cball_tightness_factors_mul_div() -> Result<(), String> {
    gold_tight(&["mul", "div"]).map(|_| ())
}
