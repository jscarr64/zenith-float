//! Harper S2: cmp_exact_diff wide-gap path with the radius one binade above zf EXPONENT_MAX (CONTRACT 16.S2.1 to
//! 16.S2.3: a radius outside zf's exponent range still compares). Hand-worked, E = EXPONENT_MAX:
//! a = 2^E (1 - 2^-200) (exponent E), b = -2^(E-130) (exponent E-129), gap 129 > EXACT_GAP so cmp_far_diff runs.
//! |a - b| = 2^E (1 - 2^-200 + 2^-130). Scaled by 2^(200-E) (an exact integer comparison): |a - b| -> 2^200 - 1 + 2^70,
//! r = 2^E = Mag(TOP, E+1) -> 2^200, so |a - b| - r -> 2^70 - 1 > 0 and r vs |a - b| is Less.
//! With b = +2^(E-130): |a - b| -> 2^200 - 1 - 2^70 < 2^200, so Greater. Negating both operands keeps |a - b|.
//! Neighbour: r = Mag(u64::MAX, E) = 2^E (1 - 2^-64) -> 2^200 - 2^136 < 2^200 - 1 + 2^70, so Less.
//! COVERS: Mag::cmp_exact_diff

use core::cmp::Ordering;

use zenith_bound::mag::Mag;
use zenith_float_num::{ExactNum, EXPONENT_MAX};

const TOP: u64 = 0x8000_0000_0000_0000;

fn mag(mant: u64, exp: i64) -> Result<Mag, String> {
    Mag::finite(mant, exp).map_err(|e| format!("Mag::finite({:#x}, {}): {:?}", mant, exp, e))
}

/// `0.mant * 2^exp` exactly at 64 bits.
fn parts(mant: u64, exp: i64) -> Result<ExactNum, String> {
    ExactNum::from_mag_parts(mant, exp, 64)
        .map_err(|e| format!("from_mag_parts({:#x}, {}): {:?}", mant, exp, e))
}

fn neg(x: &ExactNum) -> ExactNum {
    ExactNum::from(0).sub_full_prec(x)
}

fn want(what: &str, got: Option<Ordering>, w: Ordering) -> Result<(), String> {
    if got == Some(w) {
        Ok(())
    } else {
        Err(format!("{}: got {:?}, want {:?}", what, got, w))
    }
}

/// (a, 2^(E-130)) with a = 2^E (1 - 2^-200) = 2^(E-1) + (2^(E-1) - 2^(E-200)), built exactly.
fn operands() -> Result<(ExactNum, ExactNum), String> {
    let e = EXPONENT_MAX;
    let half = parts(TOP, e)?; // 2^(E-1)
    let tiny = parts(TOP, e.checked_sub(199).ok_or("E-199")?)?; // 2^(E-200)
    let rest = half.sub_full_prec(&tiny);
    let a = half.add_full_prec(&rest);
    if a.inexact() || rest.inexact() {
        return Err("a not exact".into());
    }
    let s = parts(TOP, e.checked_sub(129).ok_or("E-129")?)?; // 2^(E-130)
    Ok((a, s))
}

fn r_above() -> Result<Mag, String> {
    mag(TOP, EXPONENT_MAX.checked_add(1).ok_or("E+1")?) // 2^E, exponent E+1
}

#[test]
fn s2_diff_edge_radius_above_emax_less() -> Result<(), String> {
    let (a, s) = operands()?;
    want(
        "2^E vs |a - -s|",
        r_above()?.cmp_exact_diff(&a, &neg(&s)),
        Ordering::Less,
    )
}

#[test]
fn s2_diff_edge_radius_above_emax_mirrored_less() -> Result<(), String> {
    let (a, s) = operands()?;
    want(
        "2^E vs |-a - s|",
        r_above()?.cmp_exact_diff(&neg(&a), &s),
        Ordering::Less,
    )
}

#[test]
fn s2_diff_edge_radius_above_emax_greater() -> Result<(), String> {
    let (a, s) = operands()?;
    want(
        "2^E vs |a - s|",
        r_above()?.cmp_exact_diff(&a, &s),
        Ordering::Greater,
    )
}

#[test]
fn s2_diff_edge_radius_at_emax_less() -> Result<(), String> {
    let (a, s) = operands()?;
    let r = mag(u64::MAX, EXPONENT_MAX)?; // 2^E (1 - 2^-64), exponent E
    want(
        "2^E(1-2^-64) vs |a - -s|",
        r.cmp_exact_diff(&a, &neg(&s)),
        Ordering::Less,
    )?;
    want(
        "2^E(1-2^-64) vs |-a - s|",
        r.cmp_exact_diff(&neg(&a), &s),
        Ordering::Less,
    )
}

// Mantissa-above-TOP branch of cmp_far_diff (me == e + 1, mant != TOP => Greater), E = EXPONENT_MAX, scaled by
// 2^(200-E) as above (|a - s| -> 2^200 - 1 - 2^70, |a + s| -> 2^200 - 1 + 2^70, both < 2^200 + 2^70):
// Mag(TOP+1, E+1) = (2^63 + 1) 2^(E-63) = 2^E + 2^(E-63) -> 2^200 + 2^137 > 2^200 + 2^70, so Greater for both.
// Mag(u64::MAX, E+1) = (2^64 - 1) 2^(E-63) = 2^(E+1) - 2^(E-63) -> 2^201 - 2^137 > 2^200 + 2^70, so Greater for both.
// Mag(TOP, E+1) (the other branch) is the Less / Greater pair above. Kills: `!= TOP` to `== TOP` (TOP+1 |a+s| row goes
// to the 2^E path, Less; the TOP rows above go Greater), Greater to Less or Equal, and the branch deleted (2^E path).

fn r_above_mant(mant: u64) -> Result<Mag, String> {
    mag(mant, EXPONENT_MAX.checked_add(1).ok_or("E+1")?)
}

#[test]
fn s2_diff_edge_radius_above_emax_mant_top_plus_one_greater() -> Result<(), String> {
    let (a, s) = operands()?;
    let r = r_above_mant(TOP.checked_add(1).ok_or("TOP+1")?)?; // 2^E + 2^(E-63)
    want(
        "(2^E+2^(E-63)) vs |a - -s|",
        r.cmp_exact_diff(&a, &neg(&s)),
        Ordering::Greater,
    )?;
    want(
        "(2^E+2^(E-63)) vs |-a - s|",
        r.cmp_exact_diff(&neg(&a), &s),
        Ordering::Greater,
    )?;
    want(
        "(2^E+2^(E-63)) vs |a - s|",
        r.cmp_exact_diff(&a, &s),
        Ordering::Greater,
    )
}

#[test]
fn s2_diff_edge_radius_above_emax_mant_max_greater() -> Result<(), String> {
    let (a, s) = operands()?;
    let r = r_above_mant(u64::MAX)?; // 2^(E+1) - 2^(E-63)
    want(
        "(2^(E+1)-2^(E-63)) vs |a - -s|",
        r.cmp_exact_diff(&a, &neg(&s)),
        Ordering::Greater,
    )?;
    want(
        "(2^(E+1)-2^(E-63)) vs |-a - s|",
        r.cmp_exact_diff(&neg(&a), &s),
        Ordering::Greater,
    )?;
    want(
        "(2^(E+1)-2^(E-63)) vs |a - s|",
        r.cmp_exact_diff(&a, &s),
        Ordering::Greater,
    )
}
