//! Harper L0b: adversarial edge table for slice S1 (CONTRACT 2.2–2.4, 3.2, 6.3, 16.S1.1–9). Every row is worked by hand
//! in its comment or checked by containment of an exact image built with zf's exact `*_full_prec` arithmetic.

use zenith_bound::ball::Ball;
use zenith_bound::error::BoundError;
use zenith_bound::mag::Mag;
use zenith_bound::policy::DomainPolicy;
use zenith_float_num::{ExactNum, EXPONENT_MAX, EXPONENT_MIN};

use super::{ball, mag, num, ok, pow2, pt, want_contains, want_exact, TOP};

/// The smallest positive `Mag`, `2^(i64::MIN - 1)`; everything below the Mag range rounds up to it (CONTRACT 3.2).
fn smallest() -> Result<Mag, String> {
    mag(TOP, i64::MIN)
}

fn plus_zero(what: &str, b: &Ball) -> Result<(), String> {
    if !b.mid().is_zero() || b.mid().is_negative() {
        return Err(format!("{}: mid {:?}, want +0", what, b.mid()));
    }
    Ok(())
}

/// Builds a fresh operand Ball from its source inputs.
type Mk<'a> = &'a dyn Fn() -> Result<Ball, String>;

/// Every op applied to a fresh `a` from `mk` (and `b` for the binary ones), each result paired with its name.
fn all_ops(
    mk: &dyn Fn() -> Result<Ball, String>,
    b: &Ball,
    prec: usize,
) -> Result<Vec<(&'static str, Ball)>, String> {
    // `mk` rebuilds the operand from its source inputs, so no Ball is ever copied.
    let a = mk()?;
    Ok(vec![
        ("add", ok("add", mk()?.add(b, prec))?),
        ("sub", ok("sub", mk()?.sub(b, prec))?),
        ("mul", ok("mul", mk()?.mul(b, prec))?),
        (
            "mul_exact",
            ok("mul_exact", mk()?.mul_exact(b.mid(), prec))?,
        ),
        ("sqr", ok("sqr", mk()?.sqr(prec))?),
        ("powi3", ok("powi", mk()?.powi(3, prec))?),
        ("neg", -mk()?),
        // sqrt of |a| (a negative mid would just be the Domain row, covered by Spence's sqrt rows).
        (
            "sqrt",
            ok(
                "sqrt",
                if a.mid().is_negative() { -a } else { a }.sqrt(prec, &DomainPolicy::Clip),
            )?,
        ),
    ])
}

#[test]
fn signed_zero_operands_give_plus_zero() -> Result<(), String> {
    // -0 as an operand (built as -(0) and as 0·(-1)) behaves as +0, and every zero midpoint out is +0 (CONTRACT 6.3).
    let nz = -num(0);
    for prec in [1usize, 53, 4096] {
        let z1 = ok("exact -0", Ball::exact(nz.add_full_prec(&num(0)), prec))?;
        let z2 = ok("mul_exact -1", pt(0, prec)?.mul_exact(&num(-1), prec))?;
        let mk1 = || ok("exact -0", Ball::exact(nz.add_full_prec(&num(0)), prec));
        let mk2 = || ok("mul_exact -1", pt(0, prec)?.mul_exact(&num(-1), prec));
        let pairs: [(&Ball, Mk); 2] = [(&z1, &mk1), (&z2, &mk2)];
        for (z, mk) in pairs {
            plus_zero("operand", z)?;
            for (name, r) in all_ops(mk, z, prec)? {
                plus_zero(name, &r)?;
                if !r.rad().is_zero() {
                    return Err(format!("{} of exact zeros: rad {:?}", name, r.rad()));
                }
            }
        }
        // a zero scalar of either sign on a ball with radius: 0·(3 ± 1) = +0 ± 0.
        let b = ball(num(3), mag(TOP, 1)?, prec)?;
        want_exact(
            "(3±1)·-0",
            &ok("mul_exact", b.mul_exact(&nz, prec))?,
            &num(0),
            &Mag::zero(),
        )?;
        // cancellation to zero keeps the radius sum: (3 ± 1) - (3 ± 1) = +0 ± 2; (-5) + 5 = +0 ± 0.
        let d = ok(
            "sub",
            ball(num(3), mag(TOP, 1)?, prec)?.sub(&ball(num(3), mag(TOP, 1)?, prec)?, prec),
        )?;
        want_exact("x-x", &d, &num(0), &mag(TOP, 2)?)?;
        plus_zero("x-x", &d)?;
        plus_zero("-5+5", &ok("add", pt(-5, prec)?.add(&pt(5, prec)?, prec))?)?;
        // (-3)·0 and powi(-0, 0) = 1 exactly.
        plus_zero("-3·0", &ok("mul", pt(-3, prec)?.mul(&z1, prec))?)?;
        want_exact(
            "(-0)^0",
            &ok("powi", z2.powi(0, prec))?,
            &num(1),
            &Mag::zero(),
        )?;
    }
    Ok(())
}

#[test]
fn radius_zero_vs_nonzero_at_zero_mid() -> Result<(), String> {
    let p = 53;
    let r = mag(TOP, -9)?; // 2^-10
    let z = ball(num(0), cm(&r)?, p)?;
    let zz = pt(0, p)?;
    // (0 ± 2^-10) + (0 ± 0) = 0 ± 2^-10; (0 ± 2^-10) + (0 ± 2^-10) = 0 ± 2^-9.
    want_exact(
        "z+0",
        &ok("add", ball(num(0), cm(&r)?, p)?.add(&zz, p))?,
        &num(0),
        &r,
    )?;
    want_exact(
        "z+z",
        &ok("add", ball(num(0), cm(&r)?, p)?.add(&z, p))?,
        &num(0),
        &mag(TOP, -8)?,
    )?;
    // (0 ± 2^-10)(0 ± 2^-10) = 0 ± 2^-20 (only the r·s term); (0 ± 2^-10)·5 = 0 ± 5·2^-10.
    want_exact(
        "z*z",
        &ok("mul", ball(num(0), cm(&r)?, p)?.mul(&z, p))?,
        &num(0),
        &mag(TOP, -19)?,
    )?;
    want_exact(
        "z*5",
        &ok("mul", ball(num(0), cm(&r)?, p)?.mul(&pt(5, p)?, p))?,
        &num(0),
        &mag(0xA000_0000_0000_0000, -7)?,
    )?;
    // sqr(0 ± 2^-10) = 0 ± 2^-20; powi(0 ± 2^-10, 3) = 0 ± 2^-30 (all three are tight because the mid is 0).
    want_exact(
        "sqr z",
        &ok("sqr", ball(num(0), cm(&r)?, p)?.sqr(p))?,
        &num(0),
        &mag(TOP, -19)?,
    )?;
    want_exact(
        "z^3",
        &ok("powi", ball(num(0), cm(&r)?, p)?.powi(3, p))?,
        &num(0),
        &mag(TOP, -29)?,
    )?;
    // sqrt(0 ± 4) under Clip: endpoint ball over [0, √4] = [0, 2] is 1 ± 1 exactly (16.S1.8); Strict is Domain.
    let four = mag(TOP, 3)?;
    want_exact(
        "sqrt(0±4)",
        &ok(
            "sqrt",
            ball(num(0), cm(&four)?, p)?.sqrt(p, &DomainPolicy::Clip),
        )?,
        &num(1),
        &mag(TOP, 1)?,
    )?;
    // sqrt(4 ± 4) (lo = 0 exactly): [0, √8] — contains 0 and √8's lower bracket, both policies.
    for pol in [DomainPolicy::Clip, DomainPolicy::Strict] {
        let s = ok("sqrt(4±4)", ball(num(4), cm(&four)?, p)?.sqrt(p, &pol))?;
        want_contains("sqrt(4±4) 0", &s, &num(0))?;
        want_contains("sqrt(4±4) 2", &s, &num(2))?;
    }
    // exact zero: every op gives exactly 0 ± 0, P6 (radius zero iff exact) with zero radius in.
    for (name, b) in all_ops(&|| pt(0, p), &zz, p)? {
        want_exact(name, &b, &num(0), &Mag::zero())?;
    }
    Ok(())
}

#[test]
fn mag_exponent_extremes() -> Result<(), String> {
    let p = 64;
    let mids = [num(0), num(1), num(-3)];
    let rads = [
        mag(TOP, i64::MIN)?,
        mag(TOP, i64::MIN.checked_add(1).ok_or("min+1")?)?,
        mag(u64::MAX, i64::MIN)?,
        mag(TOP, i64::MAX)?,
        mag(u64::MAX, i64::MAX)?,
    ];
    for m in &mids {
        for r in &rads {
            let b = ball(num(2), cm(r)?, p)?;
            for (name, res) in all_ops(&|| ball(m.add_full_prec(&num(0)), cm(r)?, p), &b, p)? {
                // Every op answers; the radius never drops below the input radius for the ops that carry it whole.
                if matches!(name, "add" | "sub" | "neg") && res.rad() < r {
                    return Err(format!(
                        "{} mid {:?} rad {:?}: result rad {:?} below input",
                        name,
                        m,
                        r,
                        res.rad()
                    ));
                }
                let x = m.add_full_prec(&num(0));
                {
                    if matches!(name, "neg") {
                        want_contains(name, &res, &-x)?;
                    } else if matches!(name, "add") {
                        want_contains(name, &res, &x.add_full_prec(&num(2)))?;
                    } else if matches!(name, "sub") {
                        want_contains(name, &res, &x.sub_full_prec(&num(2)))?;
                    } else if matches!(name, "mul" | "mul_exact") {
                        want_contains(name, &res, &x.mul_full_prec(&num(2)))?;
                    } else if matches!(name, "sqr") {
                        want_contains(name, &res, &x.mul_full_prec(&x))?;
                    } else if matches!(name, "powi3") {
                        want_contains(name, &res, &x.mul_full_prec(&x).mul_full_prec(&x))?;
                    }
                }
            }
        }
    }
    Ok(())
}

#[test]
fn below_mag_range_rounds_to_smallest() -> Result<(), String> {
    // (0 ± s)(0 ± s) = 0 ± s², s² = 2^(2·MIN-2) is below the Mag range and rounds up to s itself (CONTRACT 3.2).
    let p = 53;
    let s = smallest()?;
    let z = ball(num(0), cm(&s)?, p)?;
    want_exact(
        "s*s",
        &ok("mul", ball(num(0), cm(&s)?, p)?.mul(&z, p))?,
        &num(0),
        &s,
    )?;
    want_exact(
        "sqr s",
        &ok("sqr", ball(num(0), cm(&s)?, p)?.sqr(p))?,
        &num(0),
        &s,
    )?;
    want_exact(
        "s^2",
        &ok("powi", ball(num(0), cm(&s)?, p)?.powi(2, p))?,
        &num(0),
        &s,
    )?;
    want_exact(
        "s^7",
        &ok("powi", ball(num(0), cm(&s)?, p)?.powi(7, p))?,
        &num(0),
        &s,
    )?;
    // a scalar 2^-1000 times s: still s, never 0 (a zero radius would claim exactness).
    want_exact(
        "s*2^-1000",
        &ok(
            "mul_exact",
            ball(num(0), cm(&s)?, p)?.mul_exact(&pow2(-1000)?, p),
        )?,
        &num(0),
        &s,
    )?;
    // the radius stays positive on a nonzero mid too.
    let one = ok(
        "mul",
        ball(num(1), cm(&s)?, p)?.mul(&ball(num(1), cm(&s)?, p)?, p),
    )?;
    if one.rad().is_zero() || one.rad() < &s {
        return Err(format!("(1±s)^2 rad {:?}", one.rad()));
    }
    Ok(())
}

#[test]
fn cancellation_sweep() -> Result<(), String> {
    // ((1 + 2^-k) rounded at p) - 1. For k <= p - 1 the sum is exact: the answer is 2^-k ± 0. For k >= p the sum
    // rounds to 1 ± 2^-p (half_ulp(1, p)) and the difference is 0 ± 2^-p, which must contain 2^-k.
    for p in [2usize, 24, 53, 64, 113, 256] {
        let pi = i64::try_from(p).map_err(|_| "p")?;
        for k in 1..=300i64 {
            let s = ok(
                "add",
                pt(1, p)?.add(&ok("exact", Ball::exact(pow2(-k)?, p))?, p),
            )?;
            let d = ok("sub", s.sub(&pt(1, p)?, p))?;
            let e = pow2(-k)?;
            want_contains("cancel", &d, &e)?;
            if k < pi {
                want_exact("cancel exact", &d, &e, &Mag::zero())?;
            } else {
                let want = mag(TOP, 1i64.checked_sub(pi).ok_or("e")?)?;
                want_exact("cancel lost", &d, &num(0), &want)?;
            }
        }
    }
    Ok(())
}

#[test]
fn ties_round_to_even() -> Result<(), String> {
    // 1 + 2^-p is the tie between 1 and 1 + 2^(1-p): ToEven picks 1. 1 + 3·2^-p ties between 1 + 2^(1-p) and
    // 1 + 2^(2-p): ToEven picks 1 + 2^(2-p). Just above/below the first tie (±2^-(p+40)) go up/down. ε = 2^-p.
    for p in [2usize, 3, 24, 53, 64, 113] {
        let pi = i64::try_from(p).map_err(|_| "p")?;
        let np = pi.checked_neg().ok_or("np")?;
        let eps = mag(TOP, 1i64.checked_sub(pi).ok_or("e")?)?;
        let one = num(1);
        let ulp = pow2(1i64.checked_sub(pi).ok_or("u")?)?;
        let tie = |x: ExactNum| -> Result<Ball, String> {
            ok("add", pt(1, p)?.add(&ok("exact", Ball::exact(x, 4096))?, p))
        };
        let t1 = tie(pow2(np)?)?;
        want_exact("tie down", &t1, &one, &eps)?;
        // at p = 2 the even neighbour is 2 = 0.1b·2^2, one binade up, so there ε = half_ulp(2, 2) = 2^-1.
        let eps2 = if p == 2 {
            mag(TOP, 0)?
        } else {
            mag(TOP, 1i64.checked_sub(pi).ok_or("e")?)?
        };
        let t2 = tie(num(3).mul_full_prec(&pow2(np)?))?;
        want_exact(
            "tie up",
            &t2,
            &one.add_full_prec(&ulp).add_full_prec(&ulp),
            &eps2,
        )?;
        let below = pow2(np)?.sub_full_prec(&pow2(np.checked_sub(40).ok_or("b")?)?);
        want_exact("below tie", &tie(below)?, &one, &eps)?;
        let above = pow2(np)?.add_full_prec(&pow2(np.checked_sub(40).ok_or("a")?)?);
        want_exact("above tie", &tie(above)?, &one.add_full_prec(&ulp), &eps)?;
        // the same tie reached by mul_exact: (1 + 2^-p)·1 at p from an exact 4096-bit ball.
        let b = ok("exact", Ball::exact(one.add_full_prec(&pow2(np)?), 4096))?;
        want_exact(
            "mul tie",
            &ok("mul_exact", b.mul_exact(&one, p))?,
            &one,
            &eps,
        )?;
    }
    Ok(())
}

#[test]
fn overflow_and_error_paths_answer() -> Result<(), String> {
    let p = 53;
    // 2^(EMAX-1) has zf exponent EMAX (the top binade): its double and square overflow.
    let top = pow2(EXPONENT_MAX.checked_sub(1).ok_or("t")?)?;
    let big = ok("exact", Ball::exact(top.add_full_prec(&num(0)), p))?;
    // overflowing mids give the all-reals ball, never an Err or a finite wrong ball.
    for (name, r) in all_ops(
        &|| ok("exact", Ball::exact(top.add_full_prec(&num(0)), p)),
        &big,
        p,
    )? {
        if matches!(name, "add" | "mul" | "mul_exact" | "sqr" | "powi3") && r.is_informative() {
            return Err(format!(
                "{} at 2^(EMAX-1) stayed informative: {:?}",
                name, r
            ));
        }
    }
    // tiny mids: products below EXPONENT_MIN answer (Ok) and contain 0 or the image.
    let tiny = pow2(EXPONENT_MIN.checked_add(2).ok_or("m")?)?;
    let tb = ok("exact", Ball::exact(tiny.add_full_prec(&num(0)), p))?;
    for (name, r) in all_ops(
        &|| ok("exact", Ball::exact(tiny.add_full_prec(&num(0)), p)),
        &tb,
        p,
    )? {
        if matches!(name, "mul" | "mul_exact" | "sqr" | "powi3") {
            want_contains(name, &r, &num(0))?;
        }
    }
    // powi with u64::MAX: 1, -1, 0 stay exact; 2 overflows to all reals; 1 ± 2^-60 answers and contains 1.
    let n = u64::MAX;
    want_exact(
        "1^max",
        &ok("powi", pt(1, p)?.powi(n, p))?,
        &num(1),
        &Mag::zero(),
    )?;
    want_exact(
        "(-1)^max",
        &ok("powi", pt(-1, p)?.powi(n, p))?,
        &num(-1),
        &Mag::zero(),
    )?;
    want_exact(
        "0^max",
        &ok("powi", pt(0, p)?.powi(n, p))?,
        &num(0),
        &Mag::zero(),
    )?;
    if ok("powi", pt(2, p)?.powi(n, p))?.is_informative() {
        return Err("2^(2^64-1) stayed informative".into());
    }
    want_contains(
        "(1±2^-60)^max",
        &ok("powi", ball(num(1), mag(TOP, -59)?, p)?.powi(n, p))?,
        &num(1),
    )?;
    // the all-reals ball in: add/mul/sqr/powi answer all reals; sqrt is Domain under Strict, all reals under Clip.
    for (name, r) in all_ops(&|| ok("all", Ball::all_reals(p)), &pt(1, p)?, p)? {
        if r.is_informative() && name != "mul_exact" {
            return Err(format!("{} of all reals is informative: {:?}", name, r));
        }
    }
    match ok("all", Ball::all_reals(p))?.sqrt(p, &DomainPolicy::Strict) {
        Err(BoundError::Domain { .. }) => {}
        other => return Err(format!("Strict sqrt(all reals): {:?}", other)),
    }
    // bad precisions are InvalidArgument for every op (0 and MAX + 1).
    for bad in [
        0usize,
        zenith_bound::precision::MAX_PRECISION
            .checked_add(1)
            .ok_or("max")?,
    ] {
        let r = pt(1, p)?.add(&pt(1, p)?, bad);
        if !matches!(r, Err(BoundError::InvalidArgument { .. })) {
            return Err(format!("add at prec {}: {:?}", bad, r));
        }
        let r = pt(1, p)?.powi(2, bad);
        if !matches!(r, Err(BoundError::InvalidArgument { .. })) {
            return Err(format!("powi at prec {}: {:?}", bad, r));
        }
    }
    Ok(())
}

#[test]
fn cosmology_scale() -> Result<(), String> {
    // ρ_Λ-style magnitudes: 10^122 (~2^405.3) and 2^-405, at 53 and 113 bits.
    let mut big = num(1);
    for _ in 0..122 {
        big = big.mul_full_prec(&num(10));
    }
    let mut r61 = num(1);
    for _ in 0..61 {
        r61 = r61.mul_full_prec(&num(10));
    }
    let small = pow2(-405)?;
    for p in [53usize, 113] {
        let s = ok("small", Ball::exact(small.add_full_prec(&num(0)), p))?;
        // operands are rebuilt from their exact mids for each use, never copied.
        let cb = || ok("big", Ball::exact(big.add_full_prec(&num(0)), 4096));
        let cs = || ok("small", Ball::exact(small.add_full_prec(&num(0)), p));
        want_contains(
            "big+1",
            &ok("add", cb()?.add(&pt(1, p)?, p))?,
            &big.add_full_prec(&num(1)),
        )?;
        want_contains(
            "big-small",
            &ok("sub", cb()?.sub(&s, p))?,
            &big.sub_full_prec(&small),
        )?;
        want_contains(
            "big*small",
            &ok("mul", cb()?.mul(&s, p))?,
            &big.mul_full_prec(&small),
        )?;
        want_contains(
            "sqr big",
            &ok("sqr", cb()?.sqr(p))?,
            &big.mul_full_prec(&big),
        )?;
        want_contains(
            "sqrt big",
            &ok("sqrt", cb()?.sqrt(p, &DomainPolicy::Strict))?,
            &r61,
        )?;
        let rb = ok("r61", Ball::exact(r61.add_full_prec(&num(0)), 4096))?;
        want_contains("r61^2", &ok("powi", rb.powi(2, p))?, &big)?;
        want_contains(
            "small^3",
            &ok("powi", cs()?.powi(3, p))?,
            &small.mul_full_prec(&small).mul_full_prec(&small),
        )?;
        // a 2^-405 relative radius on 10^122 survives sqr: rad ≥ 2·10^122·(10^122·2^-405).
        let rel = Mag::from_exact_num(&big.mul_full_prec(&small));
        let br = ball(big.add_full_prec(&num(0)), rel, 4096)?;
        let sq = ok("sqr rel", br.sqr(p))?;
        want_contains(
            "sqr rel hi",
            &sq,
            &big.add_full_prec(&big.mul_full_prec(&small))
                .mul_full_prec(&big.add_full_prec(&big.mul_full_prec(&small))),
        )?;
    }
    Ok(())
}

/// A fresh `Mag` equal to `m` (Mag is neither Copy nor Clone by design).
pub(super) fn cm(m: &Mag) -> Result<Mag, String> {
    if m.is_zero() {
        return Ok(Mag::zero());
    }
    if m.is_infinity() {
        return Ok(Mag::infinity());
    }
    match (m.mant(), m.exp()) {
        (Some(mant), Some(exp)) => mag(mant, exp),
        (None, _) | (_, None) => Ok(Mag::zero()),
    }
}

/// div/inv/hypot edges (16.S1.8, Spence's src in 7f5d6d9): exact quotients and roots, the formula's own hand values,
/// divisors reaching 0, the exact zero divisor, all reals in, and bad precisions.
#[test]
fn div_inv_hypot_edges() -> Result<(), String> {
    let p = 53;
    // 1/4 and 0/3 are exact: radius 0 (P6); 0/3 is +0.
    want_exact(
        "1/4",
        &ok("div", pt(1, p)?.div(&pt(4, p)?, p))?,
        &pow2(-2)?,
        &Mag::zero(),
    )?;
    let z = ok("div", pt(0, p)?.div(&pt(3, p)?, p))?;
    want_exact("0/3", &z, &num(0), &Mag::zero())?;
    plus_zero("0/3", &z)?;
    want_exact(
        "1/-2",
        &ok("inv", pt(-2, p)?.inv(p))?,
        &num(0).sub_full_prec(&pow2(-1)?),
        &Mag::zero(),
    )?;
    // (2 ± 1)/(4 ± 0) = 1/2 ± (|2|·0 + |4|·1)/(|4|·(|4| - 0)) = 1/2 ± 1/4, exact; inv(4 ± 2) = 1/4 ± 2/(4·2) = 1/4 ± 1/4.
    want_exact(
        "(2±1)/4",
        &ok("div", ball(num(2), mag(TOP, 1)?, p)?.div(&pt(4, p)?, p))?,
        &pow2(-1)?,
        &mag(TOP, -1)?,
    )?;
    want_exact(
        "1/(4±2)",
        &ok("inv", ball(num(4), mag(TOP, 2)?, p)?.inv(p))?,
        &pow2(-2)?,
        &mag(TOP, -1)?,
    )?;
    // 1/3 is inexact: a nonzero radius (P6); its containment is carried by the oracle rows.
    let third = ok("inv", pt(3, p)?.inv(p))?;
    if third.is_exact() {
        return Err(format!("1/3 claims exact: {:?}", third));
    }
    // hypot(-3, 4) = 5 and hypot(0, 0) = +0, exact; hypot of two balls around 0 contains 0.
    want_exact(
        "hypot(-3,4)",
        &ok("hypot", pt(-3, p)?.hypot(&pt(4, p)?, p))?,
        &num(5),
        &Mag::zero(),
    )?;
    let h0 = ok("hypot", pt(0, p)?.hypot(&pt(0, p)?, p))?;
    want_exact("hypot(0,0)", &h0, &num(0), &Mag::zero())?;
    plus_zero("hypot(0,0)", &h0)?;
    let hz = ok(
        "hypot",
        ball(num(0), mag(TOP, 1)?, p)?.hypot(&ball(num(0), mag(TOP, 1)?, p)?, p),
    )?;
    want_contains("hypot(0±1,0±1) 0", &hz, &num(0))?;
    want_contains("hypot(0±1,0±1) 1", &hz, &num(1))?;
    // the exact zero divisor is Domain; a divisor ball reaching 0 (1 ± 1) is all reals.
    match pt(1, p)?.div(&pt(0, p)?, p) {
        Err(BoundError::Domain { .. }) => {}
        other => return Err(format!("1/(0±0) div: {:?}", other)),
    }
    match pt(0, p)?.inv(p) {
        Err(BoundError::Domain { .. }) => {}
        other => return Err(format!("inv(0±0): {:?}", other)),
    }
    let reach = [
        (
            "div by 1±1",
            ok("div", pt(1, p)?.div(&ball(num(1), mag(TOP, 1)?, p)?, p))?,
        ),
        ("inv 1±1", ok("inv", ball(num(1), mag(TOP, 1)?, p)?.inv(p))?),
        (
            "div by all",
            ok("div", pt(1, p)?.div(&ok("all", Ball::all_reals(p))?, p))?,
        ),
        (
            "all / 2",
            ok("div", ok("all", Ball::all_reals(p))?.div(&pt(2, p)?, p))?,
        ),
        ("inv all", ok("inv", ok("all", Ball::all_reals(p))?.inv(p))?),
        (
            "hypot all",
            ok("hypot", ok("all", Ball::all_reals(p))?.hypot(&pt(1, p)?, p))?,
        ),
    ];
    for (name, r) in reach {
        if r.is_informative() {
            return Err(format!("{} is informative: {:?}", name, r));
        }
    }
    for bad in [
        0usize,
        zenith_bound::precision::MAX_PRECISION
            .checked_add(1)
            .ok_or("max")?,
    ] {
        let rs = [
            ("div", pt(1, p)?.div(&pt(2, p)?, bad)),
            ("inv", pt(2, p)?.inv(bad)),
            ("hypot", pt(1, p)?.hypot(&pt(2, p)?, bad)),
        ];
        for (name, r) in rs {
            if !matches!(r, Err(BoundError::InvalidArgument { .. })) {
                return Err(format!("{} at prec {}: {:?}", name, bad, r));
            }
        }
    }
    Ok(())
}

/// Pins Ball::sqr soundness below the zf range (206-SQRLOST is equivalent at 627fedd, see equivalent.d/arith.txt). Since zf 0cfee92 the full-precision square is exact there (inexact() == false, every bit kept), so the test also asserts the square scaled back up is exactly 1 or 9: a future silent clamp in zf fails here. Midpoints: an exact `2^e` or `3 * 2^e`, `e = EXPONENT_MIN / 2 - drop` with drop ≥ 40, whose square lies below zenith-float's range (16.S1.5). The square is nonzero,
/// so the ball must not be an exact point `v ± 0`; it must carry ε and contain a positive value.
#[test]
fn sqr_below_zf_range_is_not_exact() -> Result<(), String> {
    let half = EXPONENT_MIN.checked_div(2).ok_or("emin/2")?;
    for drop in [40i64, 63, 64, 100, 1000] {
        let e = half.checked_sub(drop).ok_or("e")?;
        let up = pow2(e.checked_neg().ok_or("-e")?)?;
        for (m, want) in [
            (pow2(e)?, 1i64),
            (ExactNum::from(3i64).mul_full_prec(&pow2(e)?), 9),
        ] {
            let sq = m.mul_full_prec(&m);
            if sq.inexact() || sq.mul_full_prec(&up).mul_full_prec(&up) != ExactNum::from(want) {
                return Err(format!(
                    "drop {}: full_prec square not exact: {:?}",
                    drop, sq
                ));
            }
            let b = ok("exact", Ball::exact(m, 53))?;
            let s = ok("sqr", b.sqr(53))?;
            if s.rad().is_zero() {
                return Err(format!(
                    "drop {}: sqr below zf range gave an exact point {:?}",
                    drop, s
                ));
            }
        }
    }
    Ok(())
}
