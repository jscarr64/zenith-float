//! CONTRACT 1.3, 14.1, 2.3, 2.4, 3.2, 3.3; oracle: dyadic-rational model (tests/model)
//! COVERS: Mag::mul_up

use crate::model::dyadic::MagV;
use crate::model::oracle::{check_up, cmp_dy};
use crate::model::sets::{extreme_exponents, magv_set, mid_exponents};
use crate::support::{ctx, dy_of, from_magv, to_magv};
use zenith_bound::mag::Mag;

pub fn verify(a: &MagV, b: &MagV, r: &MagV) -> Result<(), String> {
    let da = dy_of(a)?;
    let db = dy_of(b)?;
    let exact = da.mul(&db);
    check_up(r, da.is_zero() || db.is_zero(), &|m| cmp_dy(m, &exact))
}

fn check_pair(a: &MagV, b: &MagV) -> Result<(), String> {
    let r = to_magv(&from_magv(a)?.mul_up(&from_magv(b)?))?;
    verify(a, b, &r)
}

fn check_set(set: &[MagV]) -> Result<(), String> {
    for a in set.iter() {
        for b in set.iter() {
            ctx("mul_up", &[a, b], check_pair(a, b))?;
        }
    }
    Ok(())
}

#[test]
fn mul_up_mid_exponent_set_is_smallest_upper_bound() -> Result<(), String> {
    check_set(&magv_set(&mid_exponents()))
}

#[test]
fn mul_up_extreme_exponent_set_covers_underflow_and_overflow() -> Result<(), String> {
    check_set(&magv_set(&extreme_exponents()))
}

#[test]
fn mul_up_infinity_operand_gives_infinity_even_for_zero() -> Result<(), String> {
    let set = magv_set(&mid_exponents());
    for x in set.iter() {
        // x * Infinity
        let r1 = to_magv(&from_magv(x)?.mul_up(&Mag::infinity()))?;
        if r1 != MagV::Infinity {
            return Err(String::from("x * Infinity != Infinity"));
        }

        // Infinity * x
        let r2 = to_magv(&Mag::infinity().mul_up(&from_magv(x)?))?;
        if r2 != MagV::Infinity {
            return Err(String::from("Infinity * x != Infinity"));
        }
    }

    // Infinity * Infinity
    let r3 = to_magv(&Mag::infinity().mul_up(&Mag::infinity()))?;
    if r3 != MagV::Infinity {
        return Err(String::from("Infinity * Infinity != Infinity"));
    }

    Ok(())
}

#[test]
fn mul_up_is_commutative() -> Result<(), String> {
    let mut set = magv_set(&mid_exponents());
    set.push(MagV::Infinity);
    for a in set.iter() {
        for b in set.iter() {
            let ab = to_magv(&from_magv(a)?.mul_up(&from_magv(b)?))?;
            let ba = to_magv(&from_magv(b)?.mul_up(&from_magv(a)?))?;
            if ab != ba {
                return Err(format!(
                    "mul_up({:?}, {:?}) = {:?} but swapped gives {:?}",
                    a, b, ab, ba
                ));
            }
        }
    }
    Ok(())
}
