//! CONTRACT 14.1, 14.4, 2.3, 3.2, 3.3; oracle: dyadic-rational model (tests/model)
//! COVERS: Mag::div_up

use crate::model::dyadic::MagV;
use crate::model::oracle::{check_up, cmp_quotient};
use crate::model::sets::{extreme_exponents, magv_set, mid_exponents};
use crate::support::{ctx, dy_of, from_magv, to_magv};
use zenith_bound::mag::Mag;

pub fn verify(a: &MagV, b: &MagV, r: &MagV) -> Result<(), String> {
    let da = dy_of(a)?;
    let db = dy_of(b)?;
    check_up(r, da.is_zero(), &|m| cmp_quotient(m, &da, &db))
}

fn check_pair(a: &MagV, b: &MagV) -> Result<(), String> {
    if dy_of(b)?.is_zero() {
        return Ok(());
    }
    let r = to_magv(&from_magv(a)?.div_up(&from_magv(b)?))?;
    verify(a, b, &r)
}

fn check_set(set: &[MagV], op_name: &str) -> Result<(), String> {
    for a in set.iter() {
        for b in set.iter() {
            let result = check_pair(a, b);
            ctx(op_name, &[a, b], result)?;
        }
    }
    Ok(())
}

#[test]
fn div_up_mid_exponent_set_is_smallest_upper_bound() -> Result<(), String> {
    let set = magv_set(&mid_exponents());
    check_set(&set, "div_up")?;
    Ok(())
}

#[test]
fn div_up_extreme_exponent_set_covers_underflow_and_overflow() -> Result<(), String> {
    let set = magv_set(&extreme_exponents());
    check_set(&set, "div_up")?;
    Ok(())
}

#[test]
fn div_up_infinity_numerator_gives_infinity() -> Result<(), String> {
    let set = magv_set(&mid_exponents());
    for x in set.iter() {
        if dy_of(x)?.is_zero() {
            continue;
        }
        let r = to_magv(&Mag::infinity().div_up(&from_magv(x)?))?;
        if !matches!(r, MagV::Infinity) {
            return Err(String::from("div_up(infinity, nonzero) must be infinity"));
        }
    }
    Ok(())
}
