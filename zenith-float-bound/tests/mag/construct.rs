//! CONTRACT 2.1, 15.4; oracle: the clauses (exact expected values); Mag has private fields and is built only through its constructors
//! COVERS: BoundError::Representation, Mag::exp, Mag::finite, Mag::infinity, Mag::is_infinity, Mag::is_zero, Mag::mant, Mag::zero, Eq for Mag, PartialEq for Mag

use zenith_bound::error::BoundError;
use zenith_bound::mag::Mag;

fn made(mant: u64, exp: i64) -> Result<Mag, String> {
    Mag::finite(mant, exp).map_err(|e| format!("Mag::finite({}, {}) failed: {:?}", mant, exp, e))
}

fn parts_are(m: &Mag, exp: i64, mant: u64) -> Result<(), String> {
    match (m.exp(), m.mant()) {
        (Some(e), Some(m)) if e == exp && m == mant => Ok(()),
        _ => Err(format!(
            "expected exp {} mant {}, got {:?} {:?}",
            exp,
            mant,
            m.exp(),
            m.mant()
        )),
    }
}

#[test]
fn finite_keeps_a_canonical_mantissa() -> Result<(), String> {
    let test_cases = [
        (1u64 << 63, i64::MIN),
        ((1u64 << 63) + 1, -65),
        (0xAAAA_AAAA_AAAA_AAAB, 0),
        (u64::MAX, 64),
        (1u64 << 63, i64::MAX),
    ];

    for &(mant, exp) in &test_cases {
        let m = made(mant, exp)?;
        parts_are(&m, exp, mant)?;
    }
    Ok(())
}

#[test]
fn finite_normalises_a_small_mantissa() -> Result<(), String> {
    let test_cases = [
        (1, 64, 1, 1u64 << 63),
        (1u64 << 62, 10, 9, 1u64 << 63),
        (3, 0, -62, 3u64 << 62),
        (0x00FF_FFFF_FFFF_FFFF, 100, 92, 0xFFFF_FFFF_FFFF_FF00),
    ];

    for &(mant, exp, want_exp, want_mant) in &test_cases {
        let m = made(mant, exp)?;
        parts_are(&m, want_exp, want_mant)?;
    }
    Ok(())
}

#[test]
fn finite_rejects_zero_mantissa_with_a_nonzero_exponent() -> Result<(), String> {
    let test_cases = [1, -1, 64, i64::MIN, i64::MAX];

    for exp in &test_cases {
        if Mag::finite(0, *exp).is_ok() {
            return Err(format!("Mag::finite(0, {}) must be rejected", exp));
        }
    }
    Ok(())
}

#[test]
fn finite_zero_is_the_canonical_zero() -> Result<(), String> {
    let m = made(0, 0)?;
    if !m.is_zero() {
        return Err(String::from("made(0, 0) must be zero"));
    }
    parts_are(&m, 0, 0)?;
    if m != Mag::zero() {
        return Err(String::from("made(0, 0) must equal Mag::zero()"));
    }
    if Mag::zero().is_infinity() {
        return Err(String::from("Mag::zero() must not be infinity"));
    }
    Ok(())
}

#[test]
fn finite_rejects_an_exponent_that_cannot_be_represented() -> Result<(), String> {
    let bad_cases = [
        (1, i64::MIN),
        (1, i64::MIN + 62),
        (1u64 << 62, i64::MIN),
        (0x7FFF_FFFF_FFFF_FFFF, i64::MIN),
    ];

    for &(mant, exp) in &bad_cases {
        if Mag::finite(mant, exp).is_ok() {
            return Err(format!("Mag::finite({}, {}) must be rejected", mant, exp));
        }
    }

    let good_cases = [
        (1, i64::MIN + 63, i64::MIN, 1u64 << 63),
        (1u64 << 63, i64::MIN, i64::MIN, 1u64 << 63),
    ];

    for &(mant, exp, want_exp, want_mant) in &good_cases {
        let m = made(mant, exp)?;
        parts_are(&m, want_exp, want_mant)?;
    }
    Ok(())
}

#[test]
fn finite_rejections_are_the_representation_error() -> Result<(), String> {
    let cases = [
        // mant == 0 with exp != 0 (a proven zero is built with Mag::zero(), clause 2.1)
        (0u64, 1i64),
        (0, -1),
        (0, i64::MIN),
        (0, i64::MAX),
        // exp - leading_zeros(mant) below i64::MIN: the value is below the smallest positive Mag
        (1, i64::MIN),
        (1, i64::MIN + 62),
        (1u64 << 62, i64::MIN),
        (0x7FFF_FFFF_FFFF_FFFF, i64::MIN),
    ];

    for &(mant, exp) in &cases {
        let got = Mag::finite(mant, exp);
        if !matches!(got, Err(BoundError::Representation { .. })) {
            return Err(format!(
                "Mag::finite({}, {}) must be Err(BoundError::Representation), got {:?}",
                mant, exp, got
            ));
        }
    }
    Ok(())
}

#[test]
fn eq_treats_all_zeros_as_equal_and_infinity_only_equal_to_itself() -> Result<(), String> {
    let zero = Mag::zero();
    let infinity = Mag::infinity();

    if zero != made(0, 0)? {
        return Err(String::from("zero must equal made(0, 0)"));
    }
    if zero != Mag::zero() {
        return Err(String::from("zero must equal itself"));
    }
    if infinity != Mag::infinity() {
        return Err(String::from("infinity must equal itself"));
    }
    if zero == infinity {
        return Err(String::from("zero must not equal infinity"));
    }

    let nonzero = made(1u64 << 63, 0)?;
    if nonzero == zero {
        return Err(String::from("nonzero must not equal zero"));
    }
    if nonzero == infinity {
        return Err(String::from("nonzero must not equal infinity"));
    }

    let a = made(1, 64)?;
    let b = made(1u64 << 63, 1)?;
    if a != b {
        return Err(String::from("made(1, 64) must equal made(1u64 << 63, 1)"));
    }
    Ok(())
}

#[test]
fn accessors_describe_zero_and_infinity() -> Result<(), String> {
    let infinity = Mag::infinity();
    if infinity.exp().is_some() || infinity.mant().is_some() {
        return Err(String::from("infinity must have None exp and mant"));
    }
    if !infinity.is_infinity() {
        return Err(String::from("infinity must be infinity"));
    }
    if infinity.is_zero() {
        return Err(String::from("infinity must not be zero"));
    }

    let zero = Mag::zero();
    if !zero.is_zero() {
        return Err(String::from("zero must be zero"));
    }
    Ok(())
}
