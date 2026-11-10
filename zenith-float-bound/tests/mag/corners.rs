//! CONTRACT 14.6 (Infinity and zero-divisor rules), 14.2 (round-down overflow), 3.2, 3.3, 3.4, 2.3, 2.4; oracle: the rules stated in the clauses (exact expected values)
//! COVERS: Mag::add_up, Mag::div_lower, Mag::div_up, Mag::mul_lower, Mag::mul_up, Mag::sqrt_lower, Mag::sqrt_up, Mag::sub_lower
use crate::model::dyadic::MagV;
use crate::model::oracle::{max_mag, min_pos_mag, zero_mag};
use crate::model::sets::{magv_set, mid_exponents};
use crate::support::{from_magv, to_magv};
use zenith_bound::mag::Mag;

fn check_is(op: &str, args: &[&MagV], got: &Mag, want: &MagV) -> Result<(), String> {
    if to_magv(got)? == *want {
        Ok(())
    } else {
        Err(format!(
            "{}{:?}: expected {:?}, got {:?}",
            op,
            args,
            want,
            to_magv(got)?
        ))
    }
}

#[test]
fn round_up_ops_with_an_infinity_operand_give_infinity() -> Result<(), String> {
    let inf = Mag::infinity();
    for x in magv_set(&mid_exponents()).iter() {
        let x_mag = &from_magv(x)?;
        check_is(
            "add_up",
            &[x, &MagV::Infinity],
            &x_mag.add_up(&inf),
            &MagV::Infinity,
        )?;
        check_is(
            "add_up",
            &[&MagV::Infinity, x],
            &inf.add_up(x_mag),
            &MagV::Infinity,
        )?;
        check_is(
            "mul_up",
            &[x, &MagV::Infinity],
            &x_mag.mul_up(&inf),
            &MagV::Infinity,
        )?;
        check_is(
            "mul_up",
            &[&MagV::Infinity, x],
            &inf.mul_up(x_mag),
            &MagV::Infinity,
        )?;
        check_is(
            "div_up",
            &[x, &MagV::Infinity],
            &x_mag.div_up(&inf),
            &MagV::Infinity,
        )?;
        check_is(
            "div_up",
            &[&MagV::Infinity, x],
            &inf.div_up(x_mag),
            &MagV::Infinity,
        )?;
    }
    check_is(
        "add_up",
        &[&MagV::Infinity, &MagV::Infinity],
        &inf.add_up(&inf),
        &MagV::Infinity,
    )?;
    check_is(
        "mul_up",
        &[&MagV::Infinity, &MagV::Infinity],
        &inf.mul_up(&inf),
        &MagV::Infinity,
    )?;
    check_is(
        "div_up",
        &[&MagV::Infinity, &MagV::Infinity],
        &inf.div_up(&inf),
        &MagV::Infinity,
    )?;
    check_is(
        "sqrt_up",
        &[&MagV::Infinity],
        &inf.sqrt_up(),
        &MagV::Infinity,
    )?;
    Ok(())
}

#[test]
fn div_up_by_zero_is_infinity_including_zero_over_zero() -> Result<(), String> {
    let zero = Mag::zero();
    for x in magv_set(&mid_exponents()).iter() {
        let x_mag = &from_magv(x)?;
        check_is(
            "div_up",
            &[x, &MagV::Finite { exp: 0, mant: 0 }],
            &x_mag.div_up(&zero),
            &MagV::Infinity,
        )?;
    }
    Ok(())
}

#[test]
fn round_down_ops_with_an_infinity_operand_give_zero() -> Result<(), String> {
    let inf = Mag::infinity();
    for x in magv_set(&mid_exponents()).iter() {
        let x_mag = &from_magv(x)?;
        check_is(
            "mul_lower",
            &[x, &MagV::Infinity],
            &x_mag.mul_lower(&inf),
            &zero_mag(),
        )?;
        check_is(
            "mul_lower",
            &[&MagV::Infinity, x],
            &inf.mul_lower(x_mag),
            &zero_mag(),
        )?;
        check_is(
            "sub_lower",
            &[x, &MagV::Infinity],
            &x_mag.sub_lower(&inf),
            &zero_mag(),
        )?;
        check_is(
            "sub_lower",
            &[&MagV::Infinity, x],
            &inf.sub_lower(x_mag),
            &zero_mag(),
        )?;
        check_is(
            "div_lower",
            &[x, &MagV::Infinity],
            &x_mag.div_lower(&inf),
            &zero_mag(),
        )?;
        check_is(
            "div_lower",
            &[&MagV::Infinity, x],
            &inf.div_lower(x_mag),
            &zero_mag(),
        )?;
    }
    let zero = zero_mag();
    check_is(
        "mul_lower",
        &[&MagV::Infinity, &MagV::Infinity],
        &inf.mul_lower(&inf),
        &zero,
    )?;
    check_is(
        "sub_lower",
        &[&MagV::Infinity, &MagV::Infinity],
        &inf.sub_lower(&inf),
        &zero,
    )?;
    check_is(
        "div_lower",
        &[&MagV::Infinity, &MagV::Infinity],
        &inf.div_lower(&inf),
        &zero,
    )?;
    check_is("sqrt_lower", &[&MagV::Infinity], &inf.sqrt_lower(), &zero)?;
    Ok(())
}

#[test]
fn div_lower_by_zero_is_zero() -> Result<(), String> {
    let zero = Mag::zero();
    for x in magv_set(&mid_exponents()).iter() {
        let x_mag = &from_magv(x)?;
        check_is(
            "div_lower",
            &[x, &MagV::Finite { exp: 0, mant: 0 }],
            &x_mag.div_lower(&zero),
            &zero_mag(),
        )?;
    }
    Ok(())
}

#[test]
fn round_down_overflow_is_the_largest_finite_mag() -> Result<(), String> {
    let big = from_magv(&max_mag())?;
    let tiny = from_magv(&min_pos_mag())?;
    let max = max_mag();
    let min = min_pos_mag();
    check_is("mul_lower", &[&max, &max], &big.mul_lower(&big), &max)?;
    check_is("div_lower", &[&max, &min], &big.div_lower(&tiny), &max)?;
    Ok(())
}

#[test]
fn round_up_overflow_is_infinity_and_round_up_underflow_is_the_smallest_positive_mag(
) -> Result<(), String> {
    let big = from_magv(&max_mag())?;
    let tiny = from_magv(&min_pos_mag())?;
    let max = max_mag();
    let min = min_pos_mag();
    check_is("add_up", &[&max, &max], &big.add_up(&big), &MagV::Infinity)?;
    check_is("mul_up", &[&max, &max], &big.mul_up(&big), &MagV::Infinity)?;
    check_is("div_up", &[&max, &min], &big.div_up(&tiny), &MagV::Infinity)?;
    check_is("mul_up", &[&min, &min], &tiny.mul_up(&tiny), &min_pos_mag())?;
    check_is("div_up", &[&min, &max], &tiny.div_up(&big), &min_pos_mag())?;
    Ok(())
}
