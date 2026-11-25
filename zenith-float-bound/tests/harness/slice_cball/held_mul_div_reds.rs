//! Slice CB rows for `CBall::mul` (16.CB.4) and `CBall::div` (16.CB.7), and the rows that need them.
//! Taken live from HOLD H11: the ops are implemented, so these tests run green (P12: no ignored reds).
//! COVERS: CBall::mul, CBall::div

use super::*;

#[test]
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

#[test]
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

#[test]
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
    // just clear on one axis: 2 ± 1 excludes 0, so R1 does not hold. Mid-rad
    // d = sqr(2±1)+sqr(−1±1) still reaches 0, but the rectangle itself excludes 0
    // (min |b|² = 1), so the quotient must be a finite enclosure.
    let clear = ok(
        "new",
        CBall::new(ball(2, mag(TOP, 1)?, 53)?, ball(-1, mag(TOP, 1)?, 53)?),
    )?;
    let q = ok("div", cpt(1, 1, 53)?.div(&clear, 53))?;
    want_finite_both("clear on re", &q)?;
    // corners of [1, 3] × [−2, 0]; every corner is nonzero.
    for re in [num(1), num(2), num(3)] {
        for im in [num(-2), num(-1), num(0)] {
            let (tr, ti) = cdiv_refs(&num(1), &num(1), &re, &im)?;
            want_contains("clear re", q.re(), &tr)?;
            want_contains("clear im", q.im(), &ti)?;
        }
    }
    Ok(())
}

#[test]
fn cball_div_excludes_zero_is_finite() -> Result<(), String> {
    // The measured R2 row: b.re = 2^-40 ± 2^-41 excludes 0, b.im = 0, mid-rad sqr(b.re)
    // reaches 0 (see cball_r2_row_s1_sqr_reaches_zero). The rectangle still excludes 0
    // (|b|² ≥ 2^-82), so the quotient must be a finite enclosure of 1/b, not all reals.
    for &prec in PRECS.iter() {
        let b = ok("new", CBall::new(r2_divisor_re(prec)?, pt(0, prec)?))?;
        let q = ok("div", cpt(1, 0, prec)?.div(&b, prec))?;
        want_finite_both("R2 mid-rad row", &q)?;
        want_contains("1/mid = 2^40", q.re(), &pow2(40)?)?;
        want_contains("1/lo = 2^41", q.re(), &pow2(41)?)?;
        want_contains("im of a real reciprocal", q.im(), &num(0))?;
    }
    Ok(())
}

#[test]
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

#[test]
fn cball_all_reals_propagation_mul() -> Result<(), String> {
    // mul: both output parts use the all-reals re part.
    let h = ok(
        "new",
        CBall::new(ok("all_reals", Ball::all_reals(53))?, pt(1, 53)?),
    )?;
    want_all_reals_both("mul", &ok("mul", h.mul(&cpt(1, 1, 53)?, 53))?)
}

#[test]
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

#[test]
fn cball_golden_rows_contained_mul_div() -> Result<(), String> {
    gold_contained(&["mul", "div"]).map(|_| ())
}

/// 16.CB.10 mul reference per part, exact: RREF plus the second-order terms the first-order
/// RREF omits (`re`: rA_re·rB_re + rA_im·rB_im; `im`: rA_re·rB_im + rA_im·rB_re). Same reason
/// `sqr` uses RREF + S2 (Spence 2026-10-25): a first-order RREF cannot cover the product when
/// an input radius is not small. One gold mul row at prec=2 has S2_re larger than RREF_re.
fn four_product_rref(r: &GoldRow) -> Result<[ExactNum; 2], String> {
    let b =
        r.b.as_ref()
            .ok_or_else(|| "mul gold missing b".to_string())?;
    let (ar, ai) = (&r.a[1], &r.a[3]);
    let (br, bi) = (&b[1], &b[3]);
    let so_re = ar.mul_full_prec(br).add_full_prec(&ai.mul_full_prec(bi));
    let so_im = ar.mul_full_prec(bi).add_full_prec(&ai.mul_full_prec(br));
    Ok([
        r.rref[0].add_full_prec(&so_re),
        r.rref[1].add_full_prec(&so_im),
    ])
}

#[test]
fn cball_tightness_factors_mul_div() -> Result<(), String> {
    // 16.CB.10: mul ≤ 2·(RREF + S2); div ≤ 8·RREF when d is clear of 0 by half (gold div
    // rows are drawn in that region).
    let mut n = 0;
    for r in gold_rows()?
        .iter()
        .filter(|r| r.op == "mul" || r.op == "div")
    {
        let GoldOut::C(z) = gold_eval(r)? else {
            return Err(format!("unary gold for {}", r.op));
        };
        let k = factor_log2(&r.op)?;
        if r.op == "mul" {
            let [bre, bim] = four_product_rref(r)?;
            want_tight(&format!("mul p={} re", r.prec), z.re(), &bre, k)?;
            want_tight(&format!("mul p={} im", r.prec), z.im(), &bim, k)?;
        } else {
            want_tight(&format!("div p={} re", r.prec), z.re(), &r.rref[0], k)?;
            want_tight(&format!("div p={} im", r.prec), z.im(), &r.rref[1], k)?;
        }
        n += 1;
    }
    if n == 0 {
        return Err("no mul/div gold rows".into());
    }
    Ok(())
}

#[test]
fn cball_div_thin_strip_away_from_zero() -> Result<(), String> {
    // d0013292-shaped: b.re = 0 ± 2^{-64} (reaches 0), b.im = 2^{-1000} exact.
    // The rectangle excludes 0; 1/|b| ≈ 2^{1000} does not overflow.
    let p = 64;
    let a = ok(
        "new",
        CBall::new(ball(0, mag(TOP, -10)?, p)?, ball(0, mag(TOP, -64)?, p)?),
    )?;
    let b = ok(
        "new",
        CBall::new(
            ball(0, mag(TOP, -64)?, p)?,
            ok("exact", Ball::exact(pow2(-1000)?, p))?,
        ),
    )?;
    let q = ok("div", a.div(&b, p))?;
    want_finite_both("thin strip 2^{-1000}", &q)
}
