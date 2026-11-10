//! CONTRACT 14.2, 14.4, 2.3; oracle: dyadic-rational model (tests/model)

use crate::model::dyadic::MagV;
use crate::model::oracle::{check_down, cmp_sqrt};
use crate::model::sets::{extreme_exponents, magv_set, mid_exponents};
use crate::support::{ctx, dy_of, from_magv, to_magv};

pub fn verify(a: &MagV, r: &MagV) -> Result<(), String> {
    let da = dy_of(a)?;
    check_down(r, da.is_zero(), &|m| cmp_sqrt(m, &da))
}

fn check_one(a: &MagV) -> Result<(), String> {
    let r = to_magv(&from_magv(a)?.sqrt_lower())?;
    verify(a, &r)
}

fn check_set(set: &[MagV]) -> Result<(), String> {
    for a in set.iter() {
        ctx("sqrt_lower", &[a], check_one(a))?;
    }
    Ok(())
}

#[test]
fn sqrt_lower_mid_exponent_set_is_largest_lower_bound() -> Result<(), String> {
    let set = magv_set(&mid_exponents());
    check_set(&set)
}

#[test]
fn sqrt_lower_extreme_exponent_set_is_largest_lower_bound() -> Result<(), String> {
    let set = magv_set(&extreme_exponents());
    check_set(&set)
}
