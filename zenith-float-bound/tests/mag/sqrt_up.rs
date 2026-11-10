//! CONTRACT 14.1, 14.4, 14.6, 2.3; oracle: dyadic-rational model (tests/model)
//! COVERS: Mag::sqrt_up

use crate::model::dyadic::MagV;
use crate::model::oracle::{check_up, cmp_sqrt};
use crate::model::sets::{extreme_exponents, magv_set, mid_exponents};
use crate::support::{ctx, dy_of, from_magv, to_magv};

pub fn verify(a: &MagV, r: &MagV) -> Result<(), String> {
    let da = dy_of(a)?;
    check_up(r, da.is_zero(), &|m| cmp_sqrt(m, &da))
}

fn check_one(a: &MagV) -> Result<(), String> {
    let r = to_magv(&from_magv(a)?.sqrt_up())?;
    verify(a, &r)
}

fn check_set(set: &[MagV]) -> Result<(), String> {
    for a in set.iter() {
        ctx("sqrt_up", &[a], check_one(a))?;
    }
    Ok(())
}

#[test]
fn sqrt_up_mid_exponent_set_is_smallest_upper_bound() -> Result<(), String> {
    let set = magv_set(&mid_exponents());
    check_set(&set)
}

#[test]
fn sqrt_up_extreme_exponent_set_is_smallest_upper_bound() -> Result<(), String> {
    let set = magv_set(&extreme_exponents());
    check_set(&set)
}

#[test]
fn sqrt_up_infinity_gives_infinity() -> Result<(), String> {
    let a = MagV::Infinity;
    let r = to_magv(&from_magv(&a)?.sqrt_up())?;
    if r != MagV::Infinity {
        return Err(String::from("sqrt_up(Infinity) != Infinity"));
    }
    Ok(())
}
