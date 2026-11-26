//! CONTRACT 14.2, 14.4, 2.3, 3.2; oracle: dyadic-rational model (tests/model)

use crate::model::dyadic::MagV;
use crate::model::oracle::{check_down, cmp_quotient};
use crate::model::sets::{extreme_exponents, magv_set, mid_exponents};
use crate::support::{ctx, dy_of, from_magv, to_magv};

pub fn verify(a: &MagV, b: &MagV, r: &MagV) -> Result<(), String> {
    let da = dy_of(a)?;
    let db = dy_of(b)?;
    check_down(r, da.is_zero(), &|m| cmp_quotient(m, &da, &db))
}

fn check_pair(a: &MagV, b: &MagV) -> Result<(), String> {
    if dy_of(b)?.is_zero() {
        return Ok(());
    }
    let r = to_magv(&from_magv(a)?.div_lower(&from_magv(b)?))?;
    verify(a, b, &r)
}

fn check_set(set: &[MagV]) -> Result<(), String> {
    for a in set.iter() {
        for b in set.iter() {
            ctx("div_lower", &[a, b], check_pair(a, b))?;
        }
    }
    Ok(())
}

#[test]
fn div_lower_mid_exponent_set_is_largest_lower_bound() -> Result<(), String> {
    let set = magv_set(&mid_exponents());
    check_set(&set)
}

#[test]
fn div_lower_extreme_exponent_set_is_largest_lower_bound_or_zero_on_underflow() -> Result<(), String>
{
    let set = magv_set(&extreme_exponents());
    check_set(&set)
}
