//! Harper S2 rows above 1028 bits (2026-10-29): `x = 2^E +- 2^(E-P)` holds P+1 significant bits, P in 2000..100_000, at
//! E from EXPONENT_MAX-1 down to EXPONENT_MIN+P+2. Hand-exact oracles: the least Mag at or above 2^E + 2^(E-P) is
//! 2^E (1 + 2^-63), the greatest Mag at or below 2^E - 2^(E-P) is 2^E - 2^(E-64), and every difference is a power of two
//! or 2^E - 2^(E-P), so each answer below is exact (16.S2.1, 16.S2.2, 16.S2.3).

use core::cmp::Ordering;

use zenith_float_bound::mag::Mag;
use zenith_float_num::{ExactNum, EXPONENT_MAX, EXPONENT_MIN};

const TOP: u64 = 0x8000_0000_0000_0000;

fn add(a: i64, b: i64) -> Result<i64, String> {
    a.checked_add(b)
        .ok_or_else(|| format!("{a} + {b} overflows"))
}

fn pow2(k: i64) -> Result<ExactNum, String> {
    ExactNum::from_mag_parts(TOP, add(k, 1)?, 64).map_err(|e| format!("pow2({k}): {e:?}"))
}

fn mag(mant: u64, exp: i64) -> Result<Mag, String> {
    Mag::finite(mant, exp).map_err(|e| format!("Mag::finite({mant:#x}, {exp}): {e:?}"))
}

fn mag2(k: i64) -> Result<Mag, String> {
    mag(TOP, add(k, 1)?)
}

fn same<T: PartialEq + core::fmt::Debug>(what: &str, got: T, want: T) -> Result<(), String> {
    if got == want {
        Ok(())
    } else {
        Err(format!("{what}: got {got:?}, want {want:?}"))
    }
}

/// `(E, P)` rows: both edges, the cosmology spread 2^997, and the unit binade.
fn rows() -> Result<Vec<(i64, i64)>, String> {
    let mut out = Vec::new();
    for p in [2000, 10_000, 100_000] {
        for e in [
            EXPONENT_MAX.checked_sub(1).ok_or("emax-1")?,
            997,
            0,
            -997,
            add(add(EXPONENT_MIN, p)?, 2)?,
        ] {
            out.push((e, p));
        }
    }
    Ok(out)
}

/// `cmp_exact_num` and both Mag conversions on P+1-bit points.
#[test]
fn harper_s2_wide_cmp_and_convert() -> Result<(), String> {
    for (e, p) in rows()? {
        let g = add(e, -p)?;
        let x = pow2(e)?.add_full_prec(&pow2(g)?);
        let y = pow2(e)?.sub_full_prec(&pow2(g)?);
        let up = mag(TOP | 1, add(e, 1)?)?;
        let below = mag(u64::MAX, e)?;
        let w = format!("E={e} P={p}");
        same(
            &format!("{w} 2^E vs x"),
            mag2(e)?.cmp_exact_num(&x),
            Some(Ordering::Less),
        )?;
        same(
            &format!("{w} 2^E vs -x"),
            mag2(e)?.cmp_exact_num(&(-x.add_full_prec(&ExactNum::from(0)))),
            Some(Ordering::Less),
        )?;
        same(
            &format!("{w} 2^E(1+2^-63) vs x"),
            up.cmp_exact_num(&x),
            Some(Ordering::Greater),
        )?;
        same(
            &format!("{w} 2^E vs y"),
            mag2(e)?.cmp_exact_num(&y),
            Some(Ordering::Greater),
        )?;
        same(
            &format!("{w} 2^E-2^(E-64) vs y"),
            below.cmp_exact_num(&y),
            Some(Ordering::Less),
        )?;
        same(
            &format!("{w} from_exact_num(x)"),
            Mag::from_exact_num(&x),
            up,
        )?;
        same(
            &format!("{w} from_exact_num_lower(x)"),
            Mag::from_exact_num_lower(&x),
            mag2(e)?,
        )?;
        same(
            &format!("{w} from_exact_num(y)"),
            Mag::from_exact_num(&y),
            mag2(e)?,
        )?;
        same(
            &format!("{w} from_exact_num_lower(y)"),
            Mag::from_exact_num_lower(&y),
            below,
        )?;
    }
    Ok(())
}

/// `cmp_exact_diff` on P+1-bit operands: `x - 2^E = 2^(E-P)`, `x - y = 2^(E-P+1)`, and the far pair
/// `|2^E - 2^(E-P)|` (gap P binades) against 2^E (Greater) and 2^E - 2^(E-64) (Less).
#[test]
fn harper_s2_wide_cmp_diff() -> Result<(), String> {
    for (e, p) in rows()? {
        let g = add(e, -p)?;
        let x = pow2(e)?.add_full_prec(&pow2(g)?);
        let y = pow2(e)?.sub_full_prec(&pow2(g)?);
        let w = format!("E={e} P={p}");
        same(
            &format!("{w} 2^g vs x-2^E"),
            mag2(g)?.cmp_exact_diff(&x, &pow2(e)?),
            Some(Ordering::Equal),
        )?;
        same(
            &format!("{w} 2^(g-1) vs x-2^E"),
            mag2(add(g, -1)?)?.cmp_exact_diff(&x, &pow2(e)?),
            Some(Ordering::Less),
        )?;
        same(
            &format!("{w} 2^(g+1) vs x-y"),
            mag2(add(g, 1)?)?.cmp_exact_diff(&x, &y),
            Some(Ordering::Equal),
        )?;
        same(
            &format!("{w} 2^(g+1) vs y-x"),
            mag2(add(g, 1)?)?.cmp_exact_diff(&y, &x),
            Some(Ordering::Equal),
        )?;
        same(
            &format!("{w} 2^E vs 2^E-2^g"),
            mag2(e)?.cmp_exact_diff(&pow2(e)?, &pow2(g)?),
            Some(Ordering::Greater),
        )?;
        let below = mag(u64::MAX, e)?;
        same(
            &format!("{w} 2^E-2^(E-64) vs 2^E-2^g"),
            below.cmp_exact_diff(&pow2(e)?, &pow2(g)?),
            Some(Ordering::Less),
        )?;
        same(
            &format!("{w} 2^E-2^(E-64) vs 2^g-2^E"),
            below.cmp_exact_diff(&pow2(g)?, &pow2(e)?),
            Some(Ordering::Less),
        )?;
    }
    Ok(())
}
