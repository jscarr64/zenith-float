//! CONTRACT 14.5, 2.1; oracle: dyadic-rational model (tests/model); zero least, Infinity greatest (decided 2026-10-15)
//! COVERS: Eq for Mag, Ord for Mag, PartialEq for Mag, PartialOrd for Mag

use crate::model::dyadic::MagV;
use crate::model::sets::{extreme_exponents, magv_set, mid_exponents};
use crate::support::{ctx, dcmp, dy_of, from_magv};
use core::cmp::Ordering;
use zenith_bound::mag::Mag;

fn check_pair(a: &MagV, b: &MagV) -> Result<(), String> {
    let got = from_magv(a)?.cmp(&from_magv(b)?);
    let want = dcmp(&dy_of(a)?, &dy_of(b)?)?;
    if got == want {
        Ok(())
    } else {
        Err(format!("expected {:?}, got {:?}", want, got))
    }
}

fn check_set(set: &[MagV]) -> Result<(), String> {
    for a in set.iter() {
        for b in set.iter() {
            ctx("cmp", &[a, b], check_pair(a, b))?;
        }
    }
    Ok(())
}

#[test]
fn cmp_mid_exponent_set_matches_model() -> Result<(), String> {
    check_set(&magv_set(&mid_exponents()))
}

#[test]
fn cmp_extreme_exponent_set_matches_model() -> Result<(), String> {
    check_set(&magv_set(&extreme_exponents()))
}

#[test]
fn cmp_is_antisymmetric_on_mid_exponent_set() -> Result<(), String> {
    let set = magv_set(&mid_exponents());
    for a in set.iter() {
        for b in set.iter() {
            let ord1 = from_magv(a)?.cmp(&from_magv(b)?);
            let ord2 = from_magv(b)?.cmp(&from_magv(a)?);
            if ord1 == ord2.reverse() {
                continue;
            } else {
                return Err(format!("expected {:?}, got {:?}", ord2.reverse(), ord1));
            }
        }
    }
    Ok(())
}

#[test]
fn cmp_zero_is_least_and_infinity_is_greatest_and_equal_to_itself() -> Result<(), String> {
    let zero = Mag::zero();
    if zero.cmp(&Mag::zero()) != Ordering::Equal {
        return Err(String::from("zero must equal zero"));
    }
    if Mag::infinity().cmp(&Mag::infinity()) != Ordering::Equal {
        return Err(String::from("Infinity must equal Infinity"));
    }
    if zero.cmp(&Mag::infinity()) != Ordering::Less
        || Mag::infinity().cmp(&zero) != Ordering::Greater
    {
        return Err(String::from("zero must be below Infinity"));
    }
    for x in magv_set(&extreme_exponents()).iter() {
        let m = from_magv(x)?;
        if matches!(x, MagV::Finite { mant: 0, .. }) {
            continue;
        }
        if zero.cmp(&m) != Ordering::Less || m.cmp(&zero) != Ordering::Greater {
            return Err(format!("zero must be below {:?}", x));
        }
        if m.cmp(&Mag::infinity()) != Ordering::Less || Mag::infinity().cmp(&m) != Ordering::Greater
        {
            return Err(format!("{:?} must be below Infinity", x));
        }
    }
    Ok(())
}

#[test]
fn cmp_agrees_with_eq_and_the_comparison_operators() -> Result<(), String> {
    let mut set: Vec<Mag> = Vec::new();
    for x in magv_set(&mid_exponents()).iter() {
        set.push(from_magv(x)?);
    }
    set.push(Mag::infinity());
    for a in set.iter() {
        for b in set.iter() {
            let ord = a.cmp(b);
            if (a == b) != (ord == Ordering::Equal)
                || (a < b) != (ord == Ordering::Less)
                || (a > b) != (ord == Ordering::Greater)
                || (a <= b) != (ord != Ordering::Greater)
                || (a >= b) != (ord != Ordering::Less)
            {
                return Err(format!(
                    "operators disagree with cmp for {:?} and {:?}",
                    a, b
                ));
            }
        }
    }
    Ok(())
}
