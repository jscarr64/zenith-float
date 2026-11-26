//! CONTRACT 14.2, 3.2, 2.3; oracle: dyadic-rational model (tests/model)
//! COVERS: Mag::sub_lower

use crate::model::dyadic::MagV;
use crate::model::oracle::{check_down, cmp_minus, pred, zero_mag};
use crate::model::sets::{extreme_exponents, magv_set, mid_exponents};
use crate::support::{ctx, dcmp, dy_of, from_magv, to_magv};
use core::cmp::Ordering;

pub fn verify(a: &MagV, b: &MagV, r: &MagV) -> Result<(), String> {
    let da = dy_of(a)?;
    let db = dy_of(b)?;

    if dcmp(&da, &db)? != Ordering::Greater {
        check_down(r, true, &|_m| Some(Ordering::Less))
    } else if da.far_apart(&db) {
        // the branch is chosen from the INPUT exponent gap, never from the result:
        // the exact difference is just below a, so CONTRACT 14.2 gives pred(a)
        let p = pred(a).ok_or_else(|| String::from("no predecessor"))?;
        if *r == p {
            Ok(())
        } else {
            Err(format!("expected {:?}, got {:?}", p, r))
        }
    } else {
        check_down(r, false, &|m| cmp_minus(m, &da, &db))
    }
}

fn check_pair(a: &MagV, b: &MagV) -> Result<(), String> {
    let r = to_magv(&from_magv(a)?.sub_lower(&from_magv(b)?))?;
    verify(a, b, &r)
}

fn check_set(set: &[MagV]) -> Result<(), String> {
    for a in set.iter() {
        for b in set.iter() {
            ctx("sub_lower", &[a, b], check_pair(a, b))?;
        }
    }
    Ok(())
}

#[test]
fn sub_lower_mid_exponent_set_is_largest_lower_bound() -> Result<(), String> {
    check_set(&magv_set(&mid_exponents()))
}

#[test]
fn sub_lower_extreme_exponent_set_is_largest_lower_bound() -> Result<(), String> {
    check_set(&magv_set(&extreme_exponents()))
}

#[test]
fn sub_lower_equal_operands_are_zero() -> Result<(), String> {
    for x in magv_set(&mid_exponents()).iter() {
        let r = to_magv(&from_magv(x)?.sub_lower(&from_magv(x)?))?;
        if r == zero_mag() {
            continue;
        } else {
            return Err(format!("expected zero, got {:?}", r));
        }
    }
    Ok(())
}
