//! CONTRACT 2.3, 3.3, 3.4, 15.1 (O9), PLAN C5 `from_exact_num`/`to_exact_num`; oracle: exact parts computed by hand from
//! the `0.mant * 2^exp` convention shared by `Mag` (CONTRACT 1.2) and zenith-float `exponent()`; round trips are exact.
//! COVERS: Mag::from_exact_num, Mag::from_exact_num_lower, Mag::to_exact_num

use zenith_bound::error::BoundError;
use zenith_bound::mag::Mag;
use zenith_float_num::{ExactNum, RoundingMode, EXPONENT_MAX, EXPONENT_MIN, INF_NEG, INF_POS, NAN};

const TOP: u64 = 1u64 << 63;

fn parts(m: &Mag) -> String {
    format!("{:?}/{:?}", m.mant(), m.exp())
}

fn is_parts(what: &str, m: &Mag, mant: u64, exp: i64) -> Result<(), String> {
    match (m.mant(), m.exp()) {
        (Some(a), Some(e)) if a == mant && e == exp => Ok(()),
        _ => Err(format!("{}: want {}/{}, got {}", what, mant, exp, parts(m))),
    }
}

fn pow2(e: i64, p: usize) -> Result<ExactNum, String> {
    ExactNum::pow2(e, p).map_err(|err| format!("pow2({}, {}): {:?}", e, p, err))
}

fn mag(mant: u64, exp: i64) -> Result<Mag, String> {
    Mag::finite(mant, exp).map_err(|e| format!("Mag::finite({}, {}): {:?}", mant, exp, e))
}

fn emin() -> i64 {
    EXPONENT_MIN
}

fn emax() -> i64 {
    EXPONENT_MAX
}

#[test]
fn zero_maps_to_zero_both_ways() -> Result<(), String> {
    let z = ExactNum::new(64);
    if !Mag::from_exact_num(&z).is_zero() || !Mag::from_exact_num_lower(&z).is_zero() {
        return Err("zero ExactNum did not give a zero Mag".into());
    }
    let back = Mag::zero().to_exact_num().map_err(|e| format!("{:?}", e))?;
    if !back.is_zero() {
        return Err("Mag::zero().to_exact_num() is not zero".into());
    }
    Ok(())
}

#[test]
fn nan_and_inf_are_infinity_up_and_zero_down() -> Result<(), String> {
    for (name, x) in [("NAN", NAN), ("INF_POS", INF_POS), ("INF_NEG", INF_NEG)] {
        if !Mag::from_exact_num(&x).is_infinity() {
            return Err(format!(
                "from_exact_num({}) = {}",
                name,
                parts(&Mag::from_exact_num(&x))
            ));
        }
        if !Mag::from_exact_num_lower(&x).is_zero() {
            return Err(format!("from_exact_num_lower({}) not zero", name));
        }
    }
    let ovf = ExactNum::max_value(64).mul(&ExactNum::from_word(2, 64), 64, RoundingMode::ToEven);
    if !ovf.is_nan() && !ovf.is_inf() {
        return Err("max_value * 2 did not overflow".into());
    }
    if !Mag::from_exact_num(&ovf).is_infinity() || !Mag::from_exact_num_lower(&ovf).is_zero() {
        return Err("an overflowed zf result must give Infinity up and zero down".into());
    }
    Ok(())
}

#[test]
fn exact_values_keep_their_parts_and_ignore_sign() -> Result<(), String> {
    let three = ExactNum::from_word(3, 64);
    let minus_three = -ExactNum::from_word(3, 64);
    for (name, x) in [("3", &three), ("-3", &minus_three)] {
        is_parts(name, &Mag::from_exact_num(x), 3 << 62, 2)?;
        is_parts(name, &Mag::from_exact_num_lower(x), 3 << 62, 2)?;
    }
    is_parts(
        "1",
        &Mag::from_exact_num(&ExactNum::from_word(1, 64)),
        TOP,
        1,
    )?;
    Ok(())
}

#[test]
fn sticky_bits_round_up_by_one_and_down_by_truncation() -> Result<(), String> {
    // 1 + 2^-100 at 192 bits: top 64 bits are 2^63, the 2^-100 bit is below them
    let x = ExactNum::from_word(1, 192).add(&pow2(-100, 192)?, 192, RoundingMode::ToEven);
    is_parts("up", &Mag::from_exact_num(&x), TOP + 1, 1)?;
    is_parts("down", &Mag::from_exact_num_lower(&x), TOP, 1)?;
    Ok(())
}

#[test]
fn round_up_carry_renormalises() -> Result<(), String> {
    // 1 - 2^-200 at 256 bits: top 64 bits all ones (exp 0) with sticky; up carries to 2^63 at exp 1
    let x = pow2(0, 256)?.sub(&pow2(-200, 256)?, 256, RoundingMode::ToEven);
    is_parts("up", &Mag::from_exact_num(&x), TOP, 1)?;
    is_parts("down", &Mag::from_exact_num_lower(&x), u64::MAX, 0)?;
    Ok(())
}

#[test]
fn huge_exponents_up_and_down() -> Result<(), String> {
    is_parts(
        "max64 up",
        &Mag::from_exact_num(&ExactNum::max_value(64)),
        u64::MAX,
        emax(),
    )?;
    let big = ExactNum::max_value(128);
    // all-ones 128-bit mantissa: up carries past EXPONENT_MAX (fine in i64), down truncates
    is_parts("max128 up", &Mag::from_exact_num(&big), TOP, emax() + 1)?;
    is_parts(
        "max128 down",
        &Mag::from_exact_num_lower(&big),
        u64::MAX,
        emax(),
    )?;
    Ok(())
}

#[test]
fn subnormal_tiny_values_normalise_exactly() -> Result<(), String> {
    for p in [64usize, 128, 256] {
        let tiny = ExactNum::min_positive(p);
        let pp = i64::try_from(p).map_err(|e| format!("{}", e))?;
        let want_exp = emin() - pp + 1;
        is_parts(
            "min_positive up",
            &Mag::from_exact_num(&tiny),
            TOP,
            want_exp,
        )?;
        is_parts(
            "min_positive down",
            &Mag::from_exact_num_lower(&tiny),
            TOP,
            want_exp,
        )?;
    }
    Ok(())
}

#[test]
fn to_exact_num_is_exact_and_round_trips() -> Result<(), String> {
    let three = mag(3 << 62, 2)?
        .to_exact_num()
        .map_err(|e| format!("{:?}", e))?;
    if three != ExactNum::from_word(3, 64) {
        return Err("to_exact_num(3) != 3".into());
    }
    let mants = [TOP, TOP + 1, 3 << 62, 0xAAAA_AAAA_AAAA_AAAB, u64::MAX];
    let exps = [
        emin() - 63,
        emin(),
        emin() + 1,
        -1000,
        -1,
        0,
        1,
        64,
        1000,
        emax() - 1,
        emax(),
    ];
    for &e in &exps {
        for &mt in &mants {
            if e == emin() - 63 && mt != TOP {
                continue; // below the 64-bit subnormal grid (tested in to_exact_num_rejects_unrepresentable)
            }
            let m = mag(mt, e)?;
            let x = m
                .to_exact_num()
                .map_err(|err| format!("to_exact_num {}: {:?}", parts(&m), err))?;
            is_parts("round trip up", &Mag::from_exact_num(&x), mt, e)?;
            is_parts("round trip down", &Mag::from_exact_num_lower(&x), mt, e)?;
        }
    }
    Ok(())
}

#[test]
fn to_exact_num_rejects_unrepresentable() -> Result<(), String> {
    let cases = [
        ("Infinity", Mag::infinity()),
        ("exp EXPONENT_MAX + 1", mag(TOP, emax() + 1)?),
        ("exp i64::MAX", mag(u64::MAX, i64::MAX)?),
        ("exp i64::MIN", mag(TOP, i64::MIN)?),
        ("below the grid", mag(TOP + 1, emin() - 63)?),
        ("far below", mag(TOP, emin() - 65)?),
    ];
    for (name, m) in cases {
        match m.to_exact_num() {
            Err(BoundError::Representation { .. }) => {}
            Ok(_) => return Err(format!("{}: accepted", name)),
            Err(e) => return Err(format!("{}: wrong error {:?}", name, e)),
        }
    }
    Ok(())
}

#[test]
fn wide_exponent_sweep_round_trips_or_rejects() -> Result<(), String> {
    let mut e = emin() - 70;
    let mut checked = 0u32;
    while e <= emax() + 70 {
        let m = mag(TOP, e)?;
        let in_range = e >= emin() - 63 && e <= emax();
        match (m.to_exact_num(), in_range) {
            (Ok(x), true) => is_parts("sweep", &Mag::from_exact_num(&x), TOP, e)?,
            (Err(BoundError::Representation { .. }), false) => {}
            (r, _) => {
                return Err(format!(
                    "exp {}: in_range {} got {:?}",
                    e,
                    in_range,
                    r.map(|_| ())
                ))
            }
        }
        checked = checked.checked_add(1).ok_or("overflow")?;
        e = if (emin() - 70..emin() + 70).contains(&e) || e >= emax() - 70 {
            e.checked_add(1).ok_or("overflow")?
        } else {
            {
                // i64 exponent range (2^63 wide): geometric stride, at least 2^24.
                let stride = (e.checked_abs().ok_or("overflow")? / 64).max(1 << 24);
                e.checked_add(stride).ok_or("overflow")?.min(emax() - 70)
            }
        };
    }
    if checked < 400 {
        return Err(format!("only {} exponents checked", checked));
    }
    Ok(())
}
