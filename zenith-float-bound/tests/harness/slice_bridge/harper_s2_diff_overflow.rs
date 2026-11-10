//! Harper S2 red: cmp_exact_diff on opposite-sign operands at zf EXPONENT_MAX, where |a - b| is above the range
//! (CONTRACT 16.S2.2: it must return Some for finite a, b and a finite r). The 8 D rows of bridge_cmp.txt (lines 540,
//! 541, 546, 547, 594, 595, 600, 601) as literals. E = EXPONENT_MAX, x = (2^64 - 1) 2^(E-64) (golden
//! `ffffffffffffffff^4611686018427387839`) and y = (2^65 - 1) 2^(E-65) (`1ffffffffffffffff^4611686018427387838`), both in
//! [2^(E-1), 2^E). a = -x or -y, b = +x or +y, so |a - b| >= 2^E, above every radius here (0 and Mag(TOP, 1) = 1): Less.
//! zb returns None from src/bridge.rs:62 (`cmp_exact_num`: `not_finite(x)`), reached via src/bridge.rs:117: on the
//! gap <= EXACT_GAP path (line 112) `a.sub_full_prec(b)` (line 113) overflows to Infinity (not flagged inexact).
//! COVERS: Mag::cmp_exact_diff

use core::cmp::Ordering;

use zenith_bound::mag::Mag;
use zenith_float_num::{ExactNum, EXPONENT_MAX};

use super::raw;

const TOP: u64 = 0x8000_0000_0000_0000;

/// x (two words false) or y (true), negated or not; `raw` takes little-endian words and the top exponent.
fn operand(is_y: bool, negative: bool) -> Result<ExactNum, String> {
    if is_y {
        raw(&[TOP, u64::MAX], EXPONENT_MAX, negative) // 0.1^65 * 2^E = (2^65 - 1) 2^(E-65)
    } else {
        raw(&[u64::MAX], EXPONENT_MAX, negative) // 0.1^64 * 2^E = (2^64 - 1) 2^(E-64)
    }
}

#[test]
#[ignore = "red: CONTRACT 16.S2.2: cmp_exact_diff returns None when |a-b| exceeds zf EXPONENT_MAX (opposite-sign near-EMAX operands)"]
fn s2_diff_overflow_opposite_sign_near_emax_is_less() -> Result<(), String> {
    let one = Mag::finite(TOP, 1).map_err(|e| format!("Mag 1: {:?}", e))?;
    let zero = Mag::zero();
    let mut bad = Vec::new();
    // (golden line, radius, a is -y, b is +y)
    for (line, r, a_y, b_y) in [
        (540, &zero, false, false),
        (541, &one, false, false),
        (546, &zero, false, true),
        (547, &one, false, true),
        (594, &zero, true, false),
        (595, &one, true, false),
        (600, &zero, true, true),
        (601, &one, true, true),
    ] {
        let (a, b) = (operand(a_y, true)?, operand(b_y, false)?);
        let got = r.cmp_exact_diff(&a, &b);
        if got != Some(Ordering::Less) {
            let d = a.sub_full_prec(&b);
            bad.push(format!(
                "line {}: got {:?}, want Less (a - b: inexact {}, inf {}, nan {})",
                line,
                got,
                d.inexact(),
                d.is_inf(),
                d.is_nan()
            ));
        }
    }
    if bad.is_empty() {
        Ok(())
    } else {
        Err(bad.join("; "))
    }
}
