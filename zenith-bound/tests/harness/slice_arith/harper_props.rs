//! Harper L0b: seeded property tests for slice S1 (TEST_PLAN 6 P4/P6, K8; CONTRACT 16.S1.2).
//! K8 radius monotonicity: same input mid, r1 ≤ r2, f(r2) Ok ⇒ f(r1) Ok and, where the output mid is the same, rad1 ≤ rad2.
//! True error: for exact inputs, rad ≥ |mid − exact| (containment), rad == 0 iff mid == exact (P6), and a single
//! nearest rounding costs exactly half_ulp(e, prec) (16.S1.2).

use zenith_bound::ball::Ball;
use zenith_bound::mag::Mag;
use zenith_bound::policy::DomainPolicy;
use zenith_float_num::ExactNum;

use super::harper_edges::cm;
use super::{ball, mag, num, ok, pow2, want_contains, TOP};
use crate::common::rng::SplitMix64;

const PRECS: [usize; 8] = [1, 2, 24, 53, 64, 65, 113, 1024];

fn below(rng: &mut SplitMix64, n: u64) -> Result<i64, String> {
    i64::try_from(rng.below(n)).map_err(|_| "below".to_string())
}

/// An exact random mid with up to `bits` bits at exponent about `e`, sign random (zero 1 in 16).
fn rmid(rng: &mut SplitMix64, bits: u64, e: i64) -> Result<ExactNum, String> {
    if rng.below(16) == 0 {
        return Ok(num(0));
    }
    let mut x = num(1);
    let mut left = bits;
    while left > 0 {
        let take = left.min(60);
        let chunk = i64::try_from(rng.below(1u64 << take)).map_err(|_| "chunk")?;
        x = x
            .mul_full_prec(&pow2(i64::try_from(take).map_err(|_| "take")?)?)
            .add_full_prec(&num(chunk));
        left = left.checked_sub(take).ok_or("left")?;
    }
    let shift = e
        .checked_sub(i64::try_from(bits).map_err(|_| "bits")?)
        .ok_or("shift")?;
    x = x.mul_full_prec(&pow2(shift)?);
    Ok(if rng.below(2) == 0 { -x } else { x })
}

fn rrad(rng: &mut SplitMix64, e: i64) -> Result<Mag, String> {
    let mant = rng.next_u64() | TOP;
    let off = below(rng, 200)?.checked_sub(150).ok_or("off")?;
    mag(mant, e.checked_add(off).ok_or("rexp")?)
}

/// Unary view of every op: binary ops take `other` as the fixed second operand. sqrt runs Strict (K8 excludes the
/// Clip straddle) and Domain results are `None`.
fn apply(op: usize, a: Ball, other: &Ball, prec: usize) -> Result<Option<Ball>, String> {
    let r = match op {
        0 => a.add(other, prec),
        1 => a.sub(other, prec),
        2 => a.mul(other, prec),
        3 => a.mul_exact(other.mid(), prec),
        4 => a.sqr(prec),
        5 => a.powi(5, prec),
        6 => Ok(-a),
        7 => a.sqrt(prec, &DomainPolicy::Strict),
        8 => a.div(other, prec),
        9 => a.inv(prec),
        10 => a.hypot(other, prec),
        other_op => return Err(format!("op {}", other_op)),
    };
    match r {
        Ok(b) => Ok(Some(b)),
        Err(zenith_bound::error::BoundError::Domain { .. }) if op == 7 || op == 9 => Ok(None),
        Err(e) => Err(format!("op {}: {:?}", op, e)),
    }
}

#[test]
fn k8_radius_monotone_same_mid() -> Result<(), String> {
    let mut rng = SplitMix64::new(0x2180_0000_0000_0001);
    let mut compared = 0u32;
    for case in 0..1200u32 {
        let prec = *PRECS
            .get(usize::try_from(rng.below(8)).map_err(|_| "i")?)
            .ok_or("prec")?;
        let e = below(&mut rng, 400)?.checked_sub(200).ok_or("e")?;
        let m = rmid(
            &mut rng,
            1u64.max(u64::try_from(prec).map_err(|_| "p")?.min(200)),
            e,
        )?;
        let r1 = rrad(&mut rng, e)?;
        let r2 = r1.add_up(&rrad(&mut rng, e)?);
        let oe = below(&mut rng, 40)?.checked_sub(20).ok_or("oe")?;
        let om = rmid(&mut rng, 40, oe)?;
        let other = ball(om, rrad(&mut rng, 0)?, prec)?;
        // operands are rebuilt from (m, r) for each op, never copied.
        let b1 = || ball(m.add_full_prec(&num(0)), cm(&r1)?, prec);
        let b2 = || ball(m.add_full_prec(&num(0)), cm(&r2)?, prec);
        for op in 0..8usize {
            let f2 = apply(op, b2()?, &other, prec)?;
            let f1 = apply(op, b1()?, &other, prec)?;
            match (f1, f2) {
                (None, Some(y)) => {
                    return Err(format!(
                        "case {} op {}: f(r2) = {:?} but f(r1) failed",
                        case, op, y
                    ));
                }
                (Some(x), Some(y)) => {
                    if x.mid() == y.mid() && x.rad() > y.rad() {
                        return Err(format!(
                            "case {} op {}: rad1 {:?} > rad2 {:?} (r1 {:?} <= r2 {:?})",
                            case,
                            op,
                            x.rad(),
                            y.rad(),
                            r1,
                            r2
                        ));
                    }
                    compared = compared.checked_add(1).ok_or("n")?;
                }
                (Some(_), None) | (None, None) => {}
            }
        }
    }
    // the property must actually bite: most pairs compare.
    if compared < 6000 {
        return Err(format!("only {} pairs compared", compared));
    }
    Ok(())
}

#[test]
fn radius_covers_true_error_exact_inputs() -> Result<(), String> {
    let mut rng = SplitMix64::new(0x2180_0000_0000_0002);
    let mut inexact = 0u32;
    for case in 0..1500u32 {
        let prec = *PRECS
            .get(usize::try_from(rng.below(8)).map_err(|_| "i")?)
            .ok_or("prec")?;
        let bits = 1u64.checked_add(rng.below(2 * 128)).ok_or("bits")?;
        let ae = below(&mut rng, 600)?.checked_sub(300).ok_or("e")?;
        let a = rmid(&mut rng, bits, ae)?;
        let be = below(&mut rng, 600)?.checked_sub(300).ok_or("e")?;
        let b = rmid(&mut rng, bits, be)?;
        let ba = || ok("exact", Ball::exact(a.add_full_prec(&num(0)), 4096));
        let bb = ok("exact", Ball::exact(b.add_full_prec(&num(0)), 4096))?;
        let cases: [(&str, Ball, ExactNum, bool); 5] = [
            (
                "add",
                ok("add", ba()?.add(&bb, prec))?,
                a.add_full_prec(&b),
                true,
            ),
            (
                "sub",
                ok("sub", ba()?.sub(&bb, prec))?,
                a.sub_full_prec(&b),
                true,
            ),
            (
                "mul",
                ok("mul", ba()?.mul(&bb, prec))?,
                a.mul_full_prec(&b),
                true,
            ),
            (
                "mul_exact",
                ok("mul_exact", ba()?.mul_exact(&b, prec))?,
                a.mul_full_prec(&b),
                true,
            ),
            (
                "powi3",
                ok("powi", ba()?.powi(3, prec))?,
                a.mul_full_prec(&a).mul_full_prec(&a),
                false,
            ),
        ];
        for (name, res, x, single) in cases.iter() {
            want_contains(name, res, x).map_err(|e| format!("case {}: {}", case, e))?;
            let exact = res.mid() == x;
            if exact != res.rad().is_zero() && *single {
                return Err(format!(
                    "case {} {} P6: mid==exact {} but rad {:?}",
                    case,
                    name,
                    exact,
                    res.rad()
                ));
            }
            if !exact && *single {
                inexact = inexact.checked_add(1).ok_or("n")?;
                let e = res.mid().exponent().ok_or("exp")?;
                let want = Mag::half_ulp(e, prec);
                if res.rad() != &want {
                    return Err(format!(
                        "case {} {} at {}: rad {:?}, want half_ulp {:?}",
                        case,
                        name,
                        prec,
                        res.rad(),
                        want
                    ));
                }
            }
        }
        let sq = ok("sqr", ba()?.sqr(prec))?;
        want_contains("sqr", &sq, &a.mul_full_prec(&a))?;
        if (sq.mid() == &a.mul_full_prec(&a)) != sq.rad().is_zero() {
            return Err(format!("case {} sqr P6: {:?}", case, sq));
        }
    }
    if inexact < 1000 {
        return Err(format!("only {} inexact rows", inexact));
    }
    Ok(())
}

/// powi carries guard bits through its squaring chain, so an exact input's power is within one ulp at `prec`
/// (one final rounding plus a chain error far below it). Without the guard the chain error grows like n·ulp.
#[test]
fn powi_exact_input_within_one_ulp() -> Result<(), String> {
    let mut rng = SplitMix64::new(0x2180_0000_0000_0003);
    let mut checked = 0u32;
    for case in 0..60u32 {
        let prec = *[24usize, 53, 64, 113]
            .get(usize::try_from(rng.below(4)).map_err(|_| "i")?)
            .ok_or("prec")?;
        let n = 100u64.checked_add(rng.below(900)).ok_or("n")?;
        let e = below(&mut rng, 20)?.checked_sub(10).ok_or("e")?;
        let x = rmid(&mut rng, u64::try_from(prec).map_err(|_| "p")?, e)?;
        let mut want = num(1);
        for _ in 0..n {
            want = want.mul_full_prec(&x);
        }
        let res = ok(
            "powi",
            ok("exact", Ball::exact(x.add_full_prec(&num(0)), prec))?.powi(n, prec),
        )?;
        if !res.is_informative() {
            continue;
        }
        want_contains("powi", &res, &want).map_err(|e| format!("case {}: {}", case, e))?;
        let me = res.mid().exponent().ok_or("exp")?;
        let ulp = Mag::half_ulp(me, prec.checked_sub(1).ok_or("p-1")?);
        if res.rad() > &ulp {
            return Err(format!(
                "case {} prec {} n {}: rad {:?} > one ulp {:?}",
                case,
                prec,
                n,
                res.rad(),
                ulp
            ));
        }
        checked = checked.checked_add(1).ok_or("c")?;
    }
    if checked < 50 {
        return Err(format!("only {} informative powi rows", checked));
    }
    Ok(())
}

/// K8 for the 16.S1.8 ops (div, inv, hypot), same generator as `k8_radius_monotone_same_mid`.
#[test]
fn k8_radius_monotone_div_inv_hypot() -> Result<(), String> {
    let mut rng = SplitMix64::new(0x2180_0000_0000_0004);
    let mut compared = 0u32;
    for case in 0..1200u32 {
        let prec = *PRECS
            .get(usize::try_from(rng.below(8)).map_err(|_| "i")?)
            .ok_or("prec")?;
        let e = below(&mut rng, 400)?.checked_sub(200).ok_or("e")?;
        let m = rmid(
            &mut rng,
            1u64.max(u64::try_from(prec).map_err(|_| "p")?.min(200)),
            e,
        )?;
        let r1 = rrad(&mut rng, e)?;
        let r2 = r1.add_up(&rrad(&mut rng, e)?);
        let oe = below(&mut rng, 40)?.checked_sub(20).ok_or("oe")?;
        let om = rmid(&mut rng, 40, oe)?;
        let other = ball(om, rrad(&mut rng, oe.checked_sub(30).ok_or("ore")?)?, prec)?;
        let b1 = || ball(m.add_full_prec(&num(0)), cm(&r1)?, prec);
        let b2 = || ball(m.add_full_prec(&num(0)), cm(&r2)?, prec);
        for op in 8..11usize {
            let f2 = apply(op, b2()?, &other, prec)?;
            let f1 = apply(op, b1()?, &other, prec)?;
            match (f1, f2) {
                (None, Some(y)) => {
                    return Err(format!(
                        "case {} op {}: f(r2) = {:?} but f(r1) failed",
                        case, op, y
                    ));
                }
                (Some(x), Some(y)) => {
                    if x.mid() == y.mid() && x.rad() > y.rad() {
                        return Err(format!(
                            "case {} op {}: rad1 {:?} > rad2 {:?} (r1 {:?} <= r2 {:?})",
                            case,
                            op,
                            x.rad(),
                            y.rad(),
                            r1,
                            r2
                        ));
                    }
                    compared = compared.checked_add(1).ok_or("n")?;
                }
                (Some(_), None) | (None, None) => {}
            }
        }
    }
    if compared < 3000 {
        return Err(format!("only {} pairs compared", compared));
    }
    Ok(())
}

fn exact_abs(x: &ExactNum) -> ExactNum {
    if x < &num(0) {
        num(0).sub_full_prec(x)
    } else {
        x.add_full_prec(&num(0))
    }
}

/// `[mid - rad, mid + rad]` as exact endpoints.
fn ends(what: &str, b: &Ball) -> Result<(ExactNum, ExactNum), String> {
    let r = b
        .rad()
        .to_exact_num()
        .map_err(|e| format!("{}: rad {:?}", what, e))?;
    Ok((b.mid().sub_full_prec(&r), b.mid().add_full_prec(&r)))
}

/// `q ∋ a / b` without dividing: `a` lies between `lo·b` and `hi·b`.
fn want_quot(what: &str, q: &Ball, a: &ExactNum, b: &ExactNum) -> Result<(), String> {
    let (lo, hi) = ends(what, q)?;
    let (x, y) = (lo.mul_full_prec(b), hi.mul_full_prec(b));
    let (x, y) = if x <= y { (x, y) } else { (y, x) };
    if x <= *a && *a <= y {
        Ok(())
    } else {
        Err(format!(
            "{}: {:?} does not contain {:?}/{:?}",
            what, q, a, b
        ))
    }
}

/// Exact inputs, 16.S1.8: containment, P6 (rad 0 iff mid exact), and the single nearest rounding of div/inv costs
/// exactly half_ulp (16.S1.2). hypot is checked by squares: `lo² ≤ a² + b² ≤ hi²`, `lo = max(mid - rad, 0)`.
#[test]
fn div_inv_hypot_true_error_exact_inputs() -> Result<(), String> {
    let mut rng = SplitMix64::new(0x2180_0000_0000_0005);
    let mut inexact = 0u32;
    for case in 0..1500u32 {
        let prec = *PRECS
            .get(usize::try_from(rng.below(8)).map_err(|_| "i")?)
            .ok_or("prec")?;
        let bits = 1u64.checked_add(rng.below(2 * 128)).ok_or("bits")?;
        let ae = below(&mut rng, 600)?.checked_sub(300).ok_or("e")?;
        let a = rmid(&mut rng, bits, ae)?;
        let be = below(&mut rng, 600)?.checked_sub(300).ok_or("e")?;
        let b = rmid(&mut rng, bits, be)?;
        let ba = || ok("exact", Ball::exact(a.add_full_prec(&num(0)), 4096));
        let bb = || ok("exact", Ball::exact(b.add_full_prec(&num(0)), 4096));
        let mut rows: Vec<(&str, Ball, ExactNum)> = Vec::new();
        if !b.is_zero() {
            let q =
                ok("div", ba()?.div(&bb()?, prec)).map_err(|e| format!("case {}: {}", case, e))?;
            want_quot("div", &q, &a, &b).map_err(|e| format!("case {}: {}", case, e))?;
            rows.push(("div", q, b.add_full_prec(&num(0))));
            let i = ok("inv", bb()?.inv(prec)).map_err(|e| format!("case {}: {}", case, e))?;
            want_quot("inv", &i, &num(1), &b).map_err(|e| format!("case {}: {}", case, e))?;
            rows.push(("inv", i, b.add_full_prec(&num(0))));
        }
        for (name, res, den) in rows.iter() {
            let top = if *name == "div" {
                a.add_full_prec(&num(0))
            } else {
                num(1)
            };
            let exact = res.mid().mul_full_prec(den) == top;
            if exact != res.rad().is_zero() {
                return Err(format!(
                    "case {} {} P6: mid exact {} but rad {:?}",
                    case,
                    name,
                    exact,
                    res.rad()
                ));
            }
            if !exact {
                inexact = inexact.checked_add(1).ok_or("n")?;
                let e = res.mid().exponent().ok_or("exp")?;
                let want = Mag::half_ulp(e, prec);
                if res.rad() != &want {
                    return Err(format!(
                        "case {} {} at {}: rad {:?}, want half_ulp {:?}",
                        case,
                        name,
                        prec,
                        res.rad(),
                        want
                    ));
                }
            }
        }
        let h =
            ok("hypot", ba()?.hypot(&bb()?, prec)).map_err(|e| format!("case {}: {}", case, e))?;
        let x2 = a.mul_full_prec(&a).add_full_prec(&b.mul_full_prec(&b));
        let (lo, hi) = ends("hypot", &h)?;
        let lo = if lo < num(0) { num(0) } else { lo };
        if lo.mul_full_prec(&lo) > x2 || hi.mul_full_prec(&hi) < x2 {
            return Err(format!(
                "case {} hypot: {:?} does not contain sqrt({:?})",
                case, h, x2
            ));
        }
        if (h.mid().mul_full_prec(h.mid()) == x2) != h.rad().is_zero() {
            return Err(format!(
                "case {} hypot P6: {:?} for |{:?}|,|{:?}|",
                case,
                h,
                exact_abs(&a),
                exact_abs(&b)
            ));
        }
    }
    if inexact < 1000 {
        return Err(format!("only {} inexact div/inv rows", inexact));
    }
    Ok(())
}
