//! Self-tests of the oracle itself: a model failure must surface as a test failure, never as a silent pass; oracle: tests/model

use crate::model::dyadic::{Dy, MagV};
use crate::model::oracle::{check_down, check_up, cmp_minus, max_mag};

#[test]
fn unalignable_exponents_make_dy_add_fail() -> Result<(), String> {
    let a = Dy::from_mag(5000, 1u64 << 63);
    let b = Dy::from_mag(0, 1u64 << 63);
    if a.add(&b).is_some() {
        return Err(String::from("expected None for a gap above MAX_ALIGN"));
    }
    Ok(())
}

#[test]
fn check_up_reports_an_error_when_the_model_cannot_compare() -> Result<(), String> {
    let result = check_up(&max_mag(), false, &|_m| None);
    if result.is_ok() {
        return Err(String::from(
            "check_up accepted a result although the model returned None",
        ));
    }
    Ok(())
}

#[test]
fn check_down_reports_an_error_when_the_model_cannot_compare() -> Result<(), String> {
    let result = check_down(
        &MagV::Finite {
            exp: 1,
            mant: 1u64 << 63,
        },
        false,
        &|_m| None,
    );
    if result.is_ok() {
        return Err(String::from(
            "check_down accepted a result although the model returned None",
        ));
    }
    Ok(())
}

#[test]
fn cmp_minus_is_none_for_unalignable_operands() -> Result<(), String> {
    let a = Dy::from_mag(5000, 1u64 << 63);
    let b = Dy::from_mag(0, 1u64 << 63);
    if cmp_minus(&max_mag(), &a, &b).is_some() {
        return Err(String::from("expected None"));
    }
    Ok(())
}
