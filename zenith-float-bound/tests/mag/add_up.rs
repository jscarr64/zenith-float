//! CONTRACT 14.1, 14.3, 14.6, 2.3; oracle: dyadic-rational model (tests/model)
//! COVERS: Mag::add_up
use crate::model::dyadic::MagV;
use crate::model::oracle::{check_up, cmp_dy, succ};
use crate::model::sets::{extreme_exponents, magv_set, mid_exponents};
use crate::support::{ctx, dcmp, dy_of, from_magv, to_magv};
use core::cmp::Ordering;
use zenith_bound::mag::Mag;

pub fn verify(a: &MagV, b: &MagV, r: &MagV) -> Result<(), String> {
    let da = dy_of(a)?;
    let db = dy_of(b)?;
    if da.far_apart(&db) {
        // the branch is chosen from the INPUT exponent gap, never from the result;
        // CONTRACT 14.3: the result is the larger operand plus one unit in its last place = succ(larger)
        let larger = match dcmp(&da, &db)? {
            Ordering::Less => b,
            Ordering::Equal | Ordering::Greater => a,
        };
        let want = succ(larger);
        if *r == want {
            Ok(())
        } else {
            Err(format!("expected {:?}, got {:?}", want, r))
        }
    } else {
        let exact = da
            .add(&db)
            .ok_or_else(|| String::from("model cannot align operands that are not far apart"))?;
        check_up(r, da.is_zero() && db.is_zero(), &|m| cmp_dy(m, &exact))
    }
}

fn check_pair(a: &MagV, b: &MagV) -> Result<(), String> {
    let r = to_magv(&from_magv(a)?.add_up(&from_magv(b)?))?;
    verify(a, b, &r)
}

fn check_set(set: &[MagV]) -> Result<(), String> {
    for a in set.iter() {
        for b in set.iter() {
            ctx("add_up", &[a, b], check_pair(a, b))?;
        }
    }
    Ok(())
}

#[test]
fn add_up_mid_exponent_set_is_smallest_upper_bound() -> Result<(), String> {
    check_set(&magv_set(&mid_exponents()))
}

#[test]
fn add_up_extreme_exponent_set_is_smallest_upper_bound() -> Result<(), String> {
    check_set(&magv_set(&extreme_exponents()))
}

#[test]
fn add_up_infinity_operand_gives_infinity() -> Result<(), String> {
    for x in magv_set(&mid_exponents()).iter() {
        let r1 = to_magv(&from_magv(x)?.add_up(&Mag::infinity()))?;
        let r2 = to_magv(&Mag::infinity().add_up(&from_magv(x)?))?;
        if r1 != MagV::Infinity {
            return Err(format!("expected Infinity, got {:?}", r1));
        }
        if r2 != MagV::Infinity {
            return Err(format!("expected Infinity, got {:?}", r2));
        }
    }
    let r3 = to_magv(&Mag::infinity().add_up(&Mag::infinity()))?;
    if r3 != MagV::Infinity {
        return Err(format!("expected Infinity, got {:?}", r3));
    }
    Ok(())
}

#[test]
fn add_up_is_commutative() -> Result<(), String> {
    let mut set = magv_set(&mid_exponents());
    set.push(MagV::Infinity);
    for a in set.iter() {
        for b in set.iter() {
            let ab = to_magv(&from_magv(a)?.add_up(&from_magv(b)?))?;
            let ba = to_magv(&from_magv(b)?.add_up(&from_magv(a)?))?;
            if ab != ba {
                return Err(format!(
                    "add_up({:?}, {:?}) = {:?} but swapped gives {:?}",
                    a, b, ab, ba
                ));
            }
        }
    }
    Ok(())
}
