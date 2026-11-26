//! CONTRACT 14.2, 3.2, 2.3; oracle: dyadic-rational model (tests/model)

use crate::model::dyadic::MagV;
use crate::model::oracle::{check_down, cmp_dy};
use crate::model::sets::{extreme_exponents, magv_set, mid_exponents};
use crate::support::{ctx, dy_of, from_magv, to_magv};

pub fn verify(a: &MagV, b: &MagV, r: &MagV) -> Result<(), String> {
    let da = dy_of(a)?;
    let db = dy_of(b)?;
    let exact = da.mul(&db);
    check_down(r, da.is_zero() || db.is_zero(), &|m| cmp_dy(m, &exact))
}

fn check_pair(a: &MagV, b: &MagV) -> Result<(), String> {
    let r = to_magv(&from_magv(a)?.mul_lower(&from_magv(b)?))?;
    verify(a, b, &r)
}

#[test]
fn mul_lower_mid_exponent_set_is_largest_lower_bound() -> Result<(), String> {
    let set = magv_set(&mid_exponents());
    for a in set.iter() {
        for b in set.iter() {
            ctx("mul_lower", &[a, b], check_pair(a, b))?;
        }
    }
    Ok(())
}

#[test]
fn mul_lower_extreme_exponent_set_is_largest_lower_bound_or_zero_on_underflow() -> Result<(), String>
{
    let set = magv_set(&extreme_exponents());
    for a in set.iter() {
        for b in set.iter() {
            ctx("mul_lower", &[a, b], check_pair(a, b))?;
        }
    }
    Ok(())
}
