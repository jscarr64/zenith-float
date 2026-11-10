//! Harper S1/S2: `far_side` (src/arith.rs) with both radii outside zenith-float's range, reached through the public
//! Tri comparisons (`Ball::gt`, `Ball::overlaps` -> `endpoint_signs` -> `far_side`). CONTRACT 16.S1.11, 16.S2.2.
//! Above the range (exponent EXPONENT_MAX + 1): every finite midpoint difference has |D| < 2^E <= r1, so
//! `r1.cmp_exact_num(D)` is Greater and `far_side` is never reached; the overlap is decided. Hand-worked, E = EXPONENT_MAX:
//! a = 0 ± 2^E, b = 2^(E-1) ± 2^E: |D| = 2^(E-1) < R = 2^(E+1), so overlaps is True.
//! Below the range: r1 = r2 = 2^(EXPONENT_MIN-101) (Mag(TOP, EXPONENT_MIN-100), below zf's subnormal floor
//! 2^(EXPONENT_MIN-64), so `to_exact_num` fails for both). a = 1 ± r1, b = 0 ± r2: D = 1, R = 2^(EXPONENT_MIN-100) < 1,
//! so D - R > 0 in exact terms, but `far_side` takes the swap path and `r2.to_exact_num().ok()?` (src/arith.rs:479) returns
//! None, so gt and overlaps are Unknown. That is the contract answer, not a bug: CONTRACT 16.S2.2 says "when neither
//! converts the comparison is None" and "every None (..., both radii outside the range) is Tri::Unknown", and 16.S1.11
//! leaves exactly that case undecided (Unknown, never False). Deciding it (True/False) would be a contract change
//! (Spence/Rusty S2 review 2026-10-31; FIXLIST ZB-32, ZB-33), so the row asserts Unknown and goes red if that ever changes.
//! COVERS: Ball::gt, Ball::overlaps

use zenith_float_bound::tri::Tri;
use zenith_float_num::{EXPONENT_MAX, EXPONENT_MIN};

use super::{ball, mag, num, pow2, TOP};

const P: usize = 53;

fn want(what: &str, got: Tri, w: Tri) -> Result<(), String> {
    if got == w {
        Ok(())
    } else {
        Err(format!("{}: got {:?}, want {:?}", what, got, w))
    }
}

#[test]
fn far_side_edge_radii_above_emax_overlap_decided() -> Result<(), String> {
    let r = || mag(TOP, EXPONENT_MAX.checked_add(1).ok_or("E+1")?); // 2^E
    let a = ball(num(0), r()?, P)?;
    let b = ball(pow2(EXPONENT_MAX.checked_sub(1).ok_or("E-1")?)?, r()?, P)?;
    want("0±2^E overlaps 2^(E-1)±2^E", a.overlaps(&b), Tri::True)?;
    want("2^(E-1)±2^E overlaps 0±2^E", b.overlaps(&a), Tri::True)
}

#[test]
fn far_side_edge_radii_below_emin_far_apart() -> Result<(), String> {
    let r = || mag(TOP, EXPONENT_MIN.checked_sub(100).ok_or("EMIN-100")?); // 2^(EMIN-101)
    let a = ball(num(1), r()?, P)?;
    let b = ball(num(0), r()?, P)?;
    // Both radii below the zf range: undecided by CONTRACT 16.S2.2 and 16.S1.11, so Unknown (never False).
    want("1±tiny gt 0±tiny", a.gt(&b), Tri::Unknown)?;
    want("1±tiny overlaps 0±tiny", a.overlaps(&b), Tri::Unknown)
}
