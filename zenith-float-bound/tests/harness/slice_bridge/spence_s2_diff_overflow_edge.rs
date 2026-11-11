//! Spence S2: cmp_exact_diff when the exact a - b overflows zf EXPONENT_MAX (CONTRACT 16.S2.2, overflow case), with the
//! radius at exponents E, E + 1 and E + 2 (E = EXPONENT_MAX), so `cmp_overflow_diff` decides Less, Equal and Greater
//! on both signs of a. Hand-worked: h = 2^(E-1) (`raw(&[TOP], E)`), x = (2^64 - 1) 2^(E-64) (`raw(&[u64::MAX], E)`),
//! y = (2^65 - 1) 2^(E-65) (`raw(&[TOP, u64::MAX], E)`). |h + h| = 2^E = Mag(TOP, E+1). |x + x| = (2^64 - 1) 2^(E-63) =
//! Mag(u64::MAX, E+1). |x + y| = 2^(E+1) - 3 2^(E-65), strictly above Mag(u64::MAX, E+1) = 2^(E+1) - 2^(E-63) (Less);
//! Mag(u64::MAX, E+1) is above |h + h| = 2^E (Greater).
//! COVERS: Mag::cmp_exact_diff

use core::cmp::Ordering;

use zenith_float_bound::mag::Mag;
use zenith_float_num::{ExactNum, EXPONENT_MAX};

use super::raw;

const TOP: u64 = 0x8000_0000_0000_0000;

/// 0 = h, 1 = x, 2 = y, with the given sign.
fn operand(which: u8, negative: bool) -> Result<ExactNum, String> {
    match which {
        0 => raw(&[TOP], EXPONENT_MAX, negative),
        1 => raw(&[u64::MAX], EXPONENT_MAX, negative),
        _ => raw(&[TOP, u64::MAX], EXPONENT_MAX, negative),
    }
}

#[test]
fn s2_diff_overflow_edge_radius_near_emax() -> Result<(), String> {
    let e1 = EXPONENT_MAX.checked_add(1).ok_or("E+1")?;
    let e2 = EXPONENT_MAX.checked_add(2).ok_or("E+2")?;
    let mut bad = Vec::new();
    // (radius mant, radius exp, a operand, b operand, a negative, want)
    for (i, (mant, exp, wa, wb, a_neg, want)) in [
        (TOP, e1, 0, 0, true, Ordering::Equal),
        (TOP, e1, 0, 0, false, Ordering::Equal),
        (TOP + 1, e1, 0, 0, true, Ordering::Greater),
        (TOP + 1, e1, 0, 0, false, Ordering::Greater),
        (u64::MAX, EXPONENT_MAX, 0, 0, true, Ordering::Less),
        (u64::MAX, e1, 1, 1, true, Ordering::Equal),
        (u64::MAX, e1, 1, 1, false, Ordering::Equal),
        (u64::MAX - 1, e1, 1, 1, true, Ordering::Less),
        (u64::MAX - 1, e1, 1, 1, false, Ordering::Less),
        (u64::MAX, EXPONENT_MAX, 1, 1, false, Ordering::Less),
        (u64::MAX, e1, 1, 2, true, Ordering::Less),
        (u64::MAX, e1, 2, 1, false, Ordering::Less),
        (u64::MAX, e1, 0, 0, true, Ordering::Greater),
        (u64::MAX, e1, 0, 0, false, Ordering::Greater),
        (u64::MAX - 1, e1, 2, 1, true, Ordering::Less),
        (u64::MAX - 1, e1, 1, 2, false, Ordering::Less),
        (TOP, e2, 1, 1, true, Ordering::Greater),
        (TOP, e2, 0, 0, false, Ordering::Greater),
    ]
    .into_iter()
    .enumerate()
    {
        let r = Mag::finite(mant, exp).map_err(|err| format!("row {}: Mag: {:?}", i, err))?;
        let (a, b) = (operand(wa, a_neg)?, operand(wb, !a_neg)?);
        let got = r.cmp_exact_diff(&a, &b);
        if got != Some(want) {
            bad.push(format!("row {}: got {:?}, want {:?}", i, got, want));
        }
        // |a - b| is symmetric: the swapped operands give the same answer.
        let swapped = r.cmp_exact_diff(&b, &a);
        if swapped != Some(want) {
            bad.push(format!(
                "row {} swapped: got {:?}, want {:?}",
                i, swapped, want
            ));
        }
    }
    if bad.is_empty() {
        Ok(())
    } else {
        Err(bad.join("; "))
    }
}
