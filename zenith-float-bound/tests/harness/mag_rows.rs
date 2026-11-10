//! TEST_PLAN 8.1 edge rows (L0a-2): `sqrt` of 1, powers of two and all-ones mantissas; `mul` of all-ones and of powers of two; and the
//! shift distances 0, 63 and 64 in `Mag::finite` normalisation, `add_up` and `sub_lower` (CONTRACT 1.3, 3.3, 14.1, 14.2). A `Mag` is
//! `mant * 2^(exp - 64)` with the top bit set. Oracle: hand-derived exact values; the two irrational roots are `isqrt(2^127)` =
//! `0xb504f333f9de6484` and `isqrt((2^64 - 1) * 2^63)` = `0xb504f333f9de6483` (neither square is exact, so `up` is one more).
//! Mutants this file must kill: 181-S-sqrt-nosticky (sqrt drops the sticky bit) and 182-F-finite-minexp (Mag::finite rejects the
//! smallest exponent `i64::MIN`).
//! COVERS: Mag::sqrt_up, Mag::sqrt_lower, Mag::finite, Mag::mul_up, Mag::mul_lower, Mag::add_up, Mag::sub_lower

use crate::common::expect_eq;
use crate::common::generators::TOP;
use zenith_bound::mag::Mag;

/// `isqrt(2^127)`: the `sqrt(2)` mantissa, rounded down.
const SQRT2_LOWER: u64 = 0xb504_f333_f9de_6484;
/// `isqrt((2^64 - 1) * 2^63)`, rounded down.
const SQRT_ONES_ODD_LOWER: u64 = 0xb504_f333_f9de_6483;

const EVEN: [i64; 5] = [-64, -2, 0, 2, 64];
const ODD: [i64; 5] = [-63, -1, 1, 3, 63];

fn mag(mant: u64, exp: i64) -> Result<Mag, String> {
    Mag::finite(mant, exp).map_err(|e| format!("Mag::finite({:#x}, {}) failed: {:?}", mant, exp, e))
}

/// Checks `sqrt_lower` and `sqrt_up` of `x` against `(lower, up)`.
fn both(label: &str, x: &Mag, lower: &Mag, up: &Mag) -> Result<(), String> {
    expect_eq(&format!("{}: sqrt_lower", label), &x.sqrt_lower(), lower)?;
    expect_eq(&format!("{}: sqrt_up", label), &x.sqrt_up(), up)
}

#[test]
fn sqrt_up_and_sqrt_lower_of_one_and_of_powers_of_two() -> Result<(), String> {
    let one = mag(TOP, 1)?;
    expect_eq("sqrt(1): sqrt_up", &one.sqrt_up(), &one)?;
    expect_eq("sqrt(1): sqrt_lower", &one.sqrt_lower(), &one)?;
    for e in ODD {
        // 2^(e - 1) with e - 1 even: the root 2^((e - 1) / 2) is exact.
        let want = mag(TOP, (e - 1) / 2 + 1)?;
        both(&format!("sqrt(2^63 at {})", e), &mag(TOP, e)?, &want, &want)?;
    }
    for e in EVEN {
        // 2^(e - 1) with e - 1 odd: the root is sqrt(2) * 2^((e - 2) / 2), irrational.
        let x = mag(TOP, e)?;
        both(
            &format!("sqrt(2^63 at {})", e),
            &x,
            &mag(SQRT2_LOWER, e / 2)?,
            &mag(SQRT2_LOWER + 1, e / 2)?,
        )?;
    }
    Ok(())
}

#[test]
fn sqrt_up_of_all_ones_mantissa_carries_and_sqrt_lower_does_not() -> Result<(), String> {
    for e in EVEN {
        // isqrt((2^64 - 1) * 2^64) = 2^64 - 1, not exact: up carries to 2^64, i.e. 2^63 at e / 2 + 1.
        let x = mag(u64::MAX, e)?;
        both(
            &format!("sqrt(all-ones at {})", e),
            &x,
            &mag(u64::MAX, e / 2)?,
            &mag(TOP, e / 2 + 1)?,
        )?;
    }
    for e in ODD {
        let x = mag(u64::MAX, e)?;
        let want_exp = (e + 1) / 2;
        both(
            &format!("sqrt(all-ones at {})", e),
            &x,
            &mag(SQRT_ONES_ODD_LOWER, want_exp)?,
            &mag(SQRT_ONES_ODD_LOWER + 1, want_exp)?,
        )?;
    }
    Ok(())
}

#[test]
fn mul_up_and_mul_lower_of_all_ones_and_powers_of_two() -> Result<(), String> {
    for (a, b) in [(0i64, 0i64), (5, -7), (64, 64), (-64, 3)] {
        // (2^64 - 1)^2 = 2^128 - 2^65 + 1: top word 2^64 - 2, remainder 1.
        let x = mag(u64::MAX, a)?;
        let y = mag(u64::MAX, b)?;
        expect_eq(
            &format!("all-ones {} * {}: mul_lower", a, b),
            &x.mul_lower(&y),
            &mag(u64::MAX - 1, a + b)?,
        )?;
        expect_eq(
            &format!("all-ones {} * {}: mul_up", a, b),
            &x.mul_up(&y),
            &mag(u64::MAX, a + b)?,
        )?;
        // 2^126 exactly: 2^63 at a + b - 1, both directions.
        let p = mag(TOP, a)?;
        let q = mag(TOP, b)?;
        let want = mag(TOP, a + b - 1)?;
        expect_eq(
            &format!("2^63 {} * {}: mul_lower", a, b),
            &p.mul_lower(&q),
            &want,
        )?;
        expect_eq(&format!("2^63 {} * {}: mul_up", a, b), &p.mul_up(&q), &want)?;
    }
    Ok(())
}

#[test]
fn finite_normalises_shifts_of_0_and_63() -> Result<(), String> {
    for e in [0i64, 63, -5, i64::MAX] {
        expect_eq(
            &format!("shift 0 at {}", e),
            &mag(TOP, e)?.mant(),
            &Some(TOP),
        )?;
        expect_eq(
            &format!("shift 0 at {}: exp", e),
            &mag(TOP, e)?.exp(),
            &Some(e),
        )?;
        // mant 1 has 63 leading zeros: exactly 2^63 at e - 63.
        expect_eq(
            &format!("shift 63 at {}", e),
            &mag(1, e)?,
            &mag(TOP, e - 63)?,
        )?;
    }
    expect_eq("zero", &mag(0, 0)?, &Mag::zero())?;
    expect_eq("zero mant, nonzero exp", &Mag::finite(0, 1).is_err(), &true)?;
    // The shift reaching i64::MIN exactly is representable; one further is not.
    expect_eq(
        "shift 63 to i64::MIN",
        &mag(1, i64::MIN + 63)?.exp(),
        &Some(i64::MIN),
    )?;
    expect_eq(
        "shift 63 to i64::MIN: mant",
        &mag(1, i64::MIN + 63)?.mant(),
        &Some(TOP),
    )?;
    expect_eq(
        "shift 63 below i64::MIN",
        &Mag::finite(1, i64::MIN + 62).is_err(),
        &true,
    )?;
    expect_eq(
        "shift 0 at i64::MIN",
        &mag(TOP, i64::MIN)?.exp(),
        &Some(i64::MIN),
    )
}

#[test]
fn add_up_at_gaps_0_63_64_65_200() -> Result<(), String> {
    for e in [0i64, 100, -100] {
        let a = mag(TOP, e)?;
        let rows = [
            ("gap 0", mag(TOP, e)?, mag(TOP, e + 1)?),
            ("gap 63 (exact)", mag(TOP, e - 63)?, mag(TOP + 1, e)?),
            (
                "gap 64 (half ulp, rounds up)",
                mag(TOP, e - 64)?,
                mag(TOP + 1, e)?,
            ),
            (
                "gap 65 (below half ulp, rounds up)",
                mag(TOP, e - 65)?,
                mag(TOP + 1, e)?,
            ),
            ("gap 200 (rounds up)", mag(TOP, e - 200)?, mag(TOP + 1, e)?),
        ];
        for (name, b, want) in rows {
            expect_eq(&format!("{} at {}: a + b", name, e), &a.add_up(&b), &want)?;
            expect_eq(&format!("{} at {}: b + a", name, e), &b.add_up(&a), &want)?;
        }
    }
    Ok(())
}

#[test]
fn sub_lower_at_gaps_63_64_65() -> Result<(), String> {
    for e in [0i64, 100, -100] {
        let a = mag(TOP, e)?;
        let rows = [
            ("gap 0", mag(TOP, e)?, Mag::zero()),
            (
                "gap 63 (exact)",
                mag(TOP, e - 63)?,
                mag(u64::MAX - 1, e - 1)?,
            ),
            ("gap 64 (exact)", mag(TOP, e - 64)?, mag(u64::MAX, e - 1)?),
            (
                "gap 65 (rounds down)",
                mag(TOP, e - 65)?,
                mag(u64::MAX, e - 1)?,
            ),
            (
                "gap 200 (rounds down)",
                mag(TOP, e - 200)?,
                mag(u64::MAX, e - 1)?,
            ),
        ];
        for (name, b, want) in rows {
            expect_eq(
                &format!("{} at {}: a - b", name, e),
                &a.sub_lower(&b),
                &want,
            )?;
        }
    }
    Ok(())
}
