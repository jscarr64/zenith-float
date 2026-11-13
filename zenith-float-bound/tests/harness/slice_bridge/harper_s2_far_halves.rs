//! Harper S2 f3 (FIXLIST ZB-32): the far-halves arm of `cmp_far_diff` (src/bridge.rs:145-165, CONTRACT 16.S2.2/16.S2.3).
//! With the radius r = Mag(TOP, e+1) = 2^e and the larger operand at exponent e, `|big| - r` is formed as two exact
//! subtractions of h = 2^(e-1) (two additions for a negative big). Rows are Rusty's PR2 review edge table
//! (pop-logs/rusty/2026-10-31/zb-pr2-review-0c4e2d5.md) I, L, M, N, and his P3 (operand precision 128..1024 bits,
//! multi-limb operands; Paul 2026-10-31, FIXLIST ZB-36 keeps TEST_PLAN P3 separate). Gap e - es = 129 > EXACT_GAP.
//! Hand-worked, scaled by 2^(200-e) (an exact integer comparison, r -> 2^200):
//! L (a = 2^e (1 - 2^-200), s = 2^(e-130)): |a + s| -> 2^200 - 1 + 2^70 > r, so Less; |a - s| -> 2^200 - 1 - 2^70 < r,
//!   so Greater; negating both operands keeps |a - b|. I is the (-a, -s) row at e = EXPONENT_MAX: Greater.
//! M (big = 2^(e-1), s = 2^(e-130)): |big +- s| -> 2^199 +- 2^70 < 2^200, so Greater for all four sign pairs.
//! N (a = 2^e - 2^(e-130), s = 2^(e-130)): |a + s| = 2^e exactly, so Equal; |a - s| = 2^e - 2^(e-129) < r, so Greater.
//! The answers do not depend on e or on the operand precision; the seeded sweep draws e over
//! [EXPONENT_MIN + 200, EXPONENT_MAX] (both ends included) so every operand stays a normal zf value.
//! COVERS: Mag::cmp_exact_diff

use core::cmp::Ordering;

use zenith_float_bound::mag::Mag;
use zenith_float_num::{ExactNum, RawKind, EXPONENT_MAX, EXPONENT_MIN};

use crate::common::rng::SplitMix64;

const TOP: u64 = 0x8000_0000_0000_0000;

/// Rusty's P3 operand precisions (multi-limb).
const P3: [usize; 4] = [128, 256, 512, 1024];

/// The 64-bit base precision plus P3.
const ALL_P: [usize; 5] = [64, 128, 256, 512, 1024];

/// Seeded exponents drawn per sweep (plus the two range ends).
const SWEEP_DRAWS: usize = 32;

/// Rows per (e, precision) from `shape_rows`: L 4 + M 4 + N 4.
const SHAPE_ROWS: usize = 12;

/// One decided row: label, what cmp_exact_diff returned, what the contract wants.
struct Row {
    what: String,
    got: Option<Ordering>,
    want: Ordering,
}

fn mag(mant: u64, exp: i64) -> Result<Mag, String> {
    Mag::finite(mant, exp).map_err(|e| format!("Mag::finite({:#x}, {}): {:?}", mant, exp, e))
}

/// `0.mant * 2^exp` exactly, stored at `p` bits (`p / 64` words).
fn parts(mant: u64, exp: i64, p: usize) -> Result<ExactNum, String> {
    let x = ExactNum::from_mag_parts(mant, exp, p)
        .map_err(|e| format!("from_mag_parts({:#x}, {}, {}): {:?}", mant, exp, p, e))?;
    let words = p.checked_div(64).ok_or("p / 64")?;
    if x.as_raw_parts().words.len() != words {
        return Err(format!(
            "from_mag_parts({:#x}, {}, {}): {} words, want {}",
            mant,
            exp,
            p,
            x.as_raw_parts().words.len(),
            words
        ));
    }
    Ok(x)
}

fn neg(x: &ExactNum) -> ExactNum {
    ExactNum::from(0).sub_full_prec(x)
}

fn sub(e: i64, k: i64) -> Result<i64, String> {
    e.checked_sub(k).ok_or_else(|| format!("{} - {}", e, k))
}

fn exact(what: &str, x: ExactNum) -> Result<ExactNum, String> {
    if x.inexact() || x.is_nan() || !matches!(x.as_raw_parts().kind, RawKind::Finite) {
        Err(format!("{} not exact and finite", what))
    } else {
        Ok(x)
    }
}

/// r = Mag(TOP, e+1) = 2^e.
fn r_halves(e: i64) -> Result<Mag, String> {
    mag(TOP, e.checked_add(1).ok_or("e+1")?)
}

/// L/I operands: a = 2^e (1 - 2^-200) = 2^(e-1) + (2^(e-1) - 2^(e-200)), s = 2^(e-130), both at `p` bits.
fn l_operands(e: i64, p: usize) -> Result<(ExactNum, ExactNum), String> {
    let half = parts(TOP, e, p)?;
    let tiny = parts(TOP, sub(e, 199)?, p)?;
    let rest = exact("rest", half.sub_full_prec(&tiny))?;
    let a = exact("a", half.add_full_prec(&rest))?;
    let s = parts(TOP, sub(e, 129)?, p)?;
    Ok((a, s))
}

/// M operands: big = 2^(e-1), s = 2^(e-130).
fn m_operands(e: i64, p: usize) -> Result<(ExactNum, ExactNum), String> {
    Ok((parts(TOP, e, p)?, parts(TOP, sub(e, 129)?, p)?))
}

/// N operands: a = 2^e - 2^(e-130), s = delta = 2^(e-130). Rusty's two = half + half is 2^e, which overflows zf at
/// e = EXPONENT_MAX, so a is formed as half + (half - delta): the same exact value, every step in range.
fn n_operands(e: i64, p: usize) -> Result<(ExactNum, ExactNum), String> {
    let half = parts(TOP, e, p)?;
    let delta = parts(TOP, sub(e, 129)?, p)?;
    let rest = exact("half - delta", half.sub_full_prec(&delta))?;
    let a = exact("a", half.add_full_prec(&rest))?;
    Ok((a, delta))
}

fn row(rows: &mut Vec<Row>, what: String, got: Option<Ordering>, want: Ordering) {
    rows.push(Row { what, got, want });
}

/// Row I: r vs |-a - -s| at (e, p), want Greater.
fn i_rows(rows: &mut Vec<Row>, e: i64, p: usize) -> Result<(), String> {
    let (a, s) = l_operands(e, p)?;
    let r = r_halves(e)?;
    row(
        rows,
        format!("I e={} p={} (-a,-s)", e, p),
        r.cmp_exact_diff(&neg(&a), &neg(&s)),
        Ordering::Greater,
    );
    Ok(())
}

fn l_rows(rows: &mut Vec<Row>, e: i64, p: usize) -> Result<(), String> {
    let (a, s) = l_operands(e, p)?;
    let (na, ns) = (neg(&a), neg(&s));
    let r = r_halves(e)?;
    row(
        rows,
        format!("L e={} p={} (a,-s)", e, p),
        r.cmp_exact_diff(&a, &ns),
        Ordering::Less,
    );
    row(
        rows,
        format!("L e={} p={} (-a,s)", e, p),
        r.cmp_exact_diff(&na, &s),
        Ordering::Less,
    );
    row(
        rows,
        format!("L e={} p={} (a,s)", e, p),
        r.cmp_exact_diff(&a, &s),
        Ordering::Greater,
    );
    row(
        rows,
        format!("L e={} p={} (-a,-s)", e, p),
        r.cmp_exact_diff(&na, &ns),
        Ordering::Greater,
    );
    Ok(())
}

fn m_rows(rows: &mut Vec<Row>, e: i64, p: usize) -> Result<(), String> {
    let (big, s) = m_operands(e, p)?;
    let (nb, ns) = (neg(&big), neg(&s));
    let r = r_halves(e)?;
    row(
        rows,
        format!("M e={} p={} (big,s)", e, p),
        r.cmp_exact_diff(&big, &s),
        Ordering::Greater,
    );
    row(
        rows,
        format!("M e={} p={} (big,-s)", e, p),
        r.cmp_exact_diff(&big, &ns),
        Ordering::Greater,
    );
    row(
        rows,
        format!("M e={} p={} (-big,s)", e, p),
        r.cmp_exact_diff(&nb, &s),
        Ordering::Greater,
    );
    row(
        rows,
        format!("M e={} p={} (-big,-s)", e, p),
        r.cmp_exact_diff(&nb, &ns),
        Ordering::Greater,
    );
    Ok(())
}

fn n_rows(rows: &mut Vec<Row>, e: i64, p: usize) -> Result<(), String> {
    let (a, s) = n_operands(e, p)?;
    let (na, ns) = (neg(&a), neg(&s));
    let r = r_halves(e)?;
    row(
        rows,
        format!("N e={} p={} (a,-s)", e, p),
        r.cmp_exact_diff(&a, &ns),
        Ordering::Equal,
    );
    row(
        rows,
        format!("N e={} p={} (a,s)", e, p),
        r.cmp_exact_diff(&a, &s),
        Ordering::Greater,
    );
    row(
        rows,
        format!("N e={} p={} (-a,s)", e, p),
        r.cmp_exact_diff(&na, &s),
        Ordering::Equal,
    );
    row(
        rows,
        format!("N e={} p={} (-a,-s)", e, p),
        r.cmp_exact_diff(&na, &ns),
        Ordering::Greater,
    );
    Ok(())
}

/// L, M and N at (e, p): `SHAPE_ROWS` rows.
fn shape_rows(rows: &mut Vec<Row>, e: i64, p: usize) -> Result<(), String> {
    l_rows(rows, e, p)?;
    m_rows(rows, e, p)?;
    n_rows(rows, e, p)
}

/// Every row must be decided as wanted and there must be exactly `want_rows`; all failures are listed.
fn check(rows: &[Row], want_rows: usize) -> Result<(), String> {
    let bad: Vec<String> = rows
        .iter()
        .filter(|r| r.got != Some(r.want))
        .map(|r| format!("{}: got {:?}, want {:?}", r.what, r.got, r.want))
        .collect();
    if rows.len() != want_rows {
        return Err(format!(
            "{} rows, want {}; failures: {:#?}",
            rows.len(),
            want_rows,
            bad
        ));
    }
    if bad.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "{} of {} far-halves rows fail:\n{:#?}",
            bad.len(),
            rows.len(),
            bad
        ))
    }
}

#[test]
fn s2_far_halves_i_emax_negative_greater() -> Result<(), String> {
    let mut rows = Vec::new();
    for p in ALL_P {
        i_rows(&mut rows, EXPONENT_MAX, p)?;
    }
    check(&rows, ALL_P.len())
}

#[test]
fn s2_far_halves_l_in_range_e() -> Result<(), String> {
    let mut rows = Vec::new();
    for e in [0, 1000] {
        l_rows(&mut rows, e, 64)?;
    }
    check(&rows, 8)
}

#[test]
fn s2_far_halves_m_big_exact_half() -> Result<(), String> {
    let mut rows = Vec::new();
    m_rows(&mut rows, 1000, 64)?;
    check(&rows, 4)
}

#[test]
fn s2_far_halves_n_sum_lands_on_r() -> Result<(), String> {
    let mut rows = Vec::new();
    n_rows(&mut rows, 1000, 64)?;
    check(&rows, 4)
}

/// Rusty's P3: I, L, M, N at operand precision 128, 256, 512 and 1024 bits (multi-limb), at e = 0, 1000, EXPONENT_MAX.
#[test]
fn s2_far_halves_p3_operand_precision() -> Result<(), String> {
    let mut rows = Vec::new();
    for p in P3 {
        i_rows(&mut rows, EXPONENT_MAX, p)?;
        for e in [0, 1000, EXPONENT_MAX] {
            shape_rows(&mut rows, e, p)?;
        }
    }
    let per_p = SHAPE_ROWS
        .checked_mul(3)
        .and_then(|n| n.checked_add(1))
        .ok_or("per_p")?;
    check(&rows, per_p.checked_mul(P3.len()).ok_or("rows")?)
}

/// Seeded sweep: L, M, N at both range ends and `SWEEP_DRAWS` SplitMix64 exponents in [EXPONENT_MIN + 200, EXPONENT_MAX],
/// each at 64, 128, 256, 512 and 1024 bits.
#[test]
fn s2_far_halves_seeded_exponents() -> Result<(), String> {
    let lo = EXPONENT_MIN.checked_add(200).ok_or("lo")?;
    let span = u64::try_from(EXPONENT_MAX.checked_sub(lo).ok_or("span")?)
        .map_err(|_| "span u64")?
        .checked_add(1)
        .ok_or("span+1")?;
    let mut es = vec![lo, EXPONENT_MAX];
    let mut rng = SplitMix64::new(0x5200_00F3);
    for _ in 0..SWEEP_DRAWS {
        let off = rng.next_u64().checked_rem(span).ok_or("rem")?;
        let off = i64::try_from(off).map_err(|_| "off i64")?;
        es.push(lo.checked_add(off).ok_or("lo+off")?);
    }
    let mut rows = Vec::new();
    for &e in &es {
        for p in ALL_P {
            shape_rows(&mut rows, e, p)?;
        }
    }
    let n_e = SWEEP_DRAWS.checked_add(2).ok_or("n_e")?;
    let want = n_e
        .checked_mul(ALL_P.len())
        .and_then(|n| n.checked_mul(SHAPE_ROWS))
        .ok_or("want")?;
    check(&rows, want)
}
