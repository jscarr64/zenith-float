//! Slice CB rows for the S1 `Ball` API that CB builds on but that no other harness file lists in a COVERS line
//! (`abs`, `mul_2exp` and the Tri comparisons of CONTRACT 15.1, 15.2, 16.S1.10, 16.S1.11). Owned by slice CB.
//! Oracles: hand-verified exact rows (worked in the comment beside each). The S1 mutants 250..=282 are the strong check.
//! COVERS: Ball::abs, Ball::mul_2exp, Ball::lt, Ball::le, Ball::gt, Ball::ge, Ball::eq, Ball::ne, Ball::overlaps

use zenith_bound::ball::Ball;
use zenith_bound::error::BoundError;
use zenith_bound::mag::Mag;
use zenith_bound::tri::Tri;
use zenith_float_num::ExactNum;

const TOP: u64 = 0x8000_0000_0000_0000;

fn ok<T>(what: &str, r: Result<T, BoundError>) -> Result<T, String> {
    r.map_err(|e| format!("{}: unexpected {:?}", what, e))
}

fn pt(v: i64) -> Result<Ball, String> {
    ok("exact", Ball::exact(ExactNum::from(v), 53))
}

/// `v ± 2^-k`; Mag 2^-k is finite(TOP, 1 - k).
fn bl(v: i64, k: i64) -> Result<Ball, String> {
    let e = 1i64.checked_sub(k).ok_or("bl k")?;
    let r = Mag::finite(TOP, e).map_err(|e| format!("Mag::finite: {:?}", e))?;
    ok("from_mid_rad", Ball::from_mid_rad(ExactNum::from(v), r, 53))
}

fn want(what: &str, got: Tri, exp: Tri) -> Result<(), String> {
    if got == exp {
        Ok(())
    } else {
        Err(format!("{}: got {:?}, want {:?}", what, got, exp))
    }
}

#[test]
fn cball_s1_cover_lt_le_gt_ge() -> Result<(), String> {
    want("1 lt 2", pt(1)?.lt(&pt(2)?), Tri::True)?;
    want("2 lt 1", pt(2)?.lt(&pt(1)?), Tri::False)?;
    want("1 lt 1", pt(1)?.lt(&pt(1)?), Tri::False)?;
    // [0, 2] vs 2: some points below, one equal.
    want("1±1 lt 2", bl(1, 0)?.lt(&pt(2)?), Tri::Unknown)?;
    want("1 le 1", pt(1)?.le(&pt(1)?), Tri::True)?;
    want("2 le 1", pt(2)?.le(&pt(1)?), Tri::False)?;
    want("1±1 le 1", bl(1, 0)?.le(&pt(1)?), Tri::Unknown)?;
    want("2 gt 1", pt(2)?.gt(&pt(1)?), Tri::True)?;
    want("1 gt 2", pt(1)?.gt(&pt(2)?), Tri::False)?;
    want("1 ge 1", pt(1)?.ge(&pt(1)?), Tri::True)?;
    want("1 ge 2", pt(1)?.ge(&pt(2)?), Tri::False)
}

#[test]
fn cball_s1_cover_eq_ne_overlaps() -> Result<(), String> {
    want("1 eq 1", pt(1)?.eq(&pt(1)?), Tri::True)?;
    want("1 eq 2", pt(1)?.eq(&pt(2)?), Tri::False)?;
    want("1±1 eq 1", bl(1, 0)?.eq(&pt(1)?), Tri::Unknown)?;
    want("1 ne 2", pt(1)?.ne(&pt(2)?), Tri::True)?;
    want("1 ne 1", pt(1)?.ne(&pt(1)?), Tri::False)?;
    want("1±1 ne 1", bl(1, 0)?.ne(&pt(1)?), Tri::Unknown)?;
    // [0, 2] touches 2; 1 and 2 are disjoint points; an all-reals ball overlaps everything.
    want("1±1 overlaps 2", bl(1, 0)?.overlaps(&pt(2)?), Tri::True)?;
    want("1 overlaps 2", pt(1)?.overlaps(&pt(2)?), Tri::False)?;
    let all = ok("all_reals", Ball::all_reals(53))?;
    want("all overlaps 5", all.overlaps(&pt(5)?), Tri::True)
}

#[test]
fn cball_s1_cover_abs_mul_2exp() -> Result<(), String> {
    // |−3| = 3 exactly.
    let a = ok("abs", pt(-3)?.abs(53))?;
    want("abs(-3) eq 3", a.eq(&pt(3)?), Tri::True)?;
    // |1 ± 2| straddles 0: the endpoint ball [0, 3] holds 0 and 3, not 4 or −1.
    let s = ok("abs", bl(1, -1)?.abs(53))?;
    want("abs(1±2) holds 0", s.overlaps(&pt(0)?), Tri::True)?;
    want("abs(1±2) holds 3", s.overlaps(&pt(3)?), Tri::True)?;
    want("abs(1±2) not 4", s.overlaps(&pt(4)?), Tri::False)?;
    want("abs(1±2) not -1", s.overlaps(&pt(-1)?), Tri::False)?;
    // 3 · 2^2 = 12 exactly.
    let m = ok("mul_2exp", pt(3)?.mul_2exp(2, 53))?;
    want("3·2^2 eq 12", m.eq(&pt(12)?), Tri::True)?;
    // (1 ± 2^-4) · 2^4 = 16 ± 1: holds 15 and 17, not 18.
    let r = ok("mul_2exp", bl(1, 4)?.mul_2exp(4, 53))?;
    if r.mid() != &ExactNum::from(16i64) {
        return Err(format!("mul_2exp mid: {:?}", r));
    }
    want("16±1 holds 17", r.overlaps(&pt(17)?), Tri::True)?;
    want("16±1 holds 15", r.overlaps(&pt(15)?), Tri::True)?;
    want("16±1 not 18", r.overlaps(&pt(18)?), Tri::False)?;
    // prec 0 is rejected.
    let a0 = pt(1)?.abs(0);
    if !matches!(a0, Err(BoundError::InvalidArgument { .. })) {
        return Err(format!("abs prec 0: {:?}", a0));
    }
    let m0 = pt(1)?.mul_2exp(1, 0);
    if !matches!(m0, Err(BoundError::InvalidArgument { .. })) {
        return Err(format!("mul_2exp prec 0: {:?}", m0));
    }
    Ok(())
}
