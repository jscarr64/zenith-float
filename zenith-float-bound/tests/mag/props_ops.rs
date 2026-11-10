//! CONTRACT 14.1-14.6, 2.3, 2.4, 3.2, 3.3, 3.4; oracles: (1) independent u128 integer oracle (tests/model/direct.rs), (2) the dyadic-rational model (the per-op `verify` functions), (3) metamorphic relations; seeded splitmix64, replay with ZB_PROP_SEED/ZB_PROP_CASES
//! COVERS: Mag::add_up, Mag::div_lower, Mag::div_up, Mag::mul_lower, Mag::mul_up, Mag::sqrt_lower, Mag::sqrt_up, Mag::sub_lower

use crate::model::direct;
use crate::model::draw::{draw_mag, draw_partner};
use crate::model::dyadic::MagV;
use crate::model::oracle::succ;
use crate::model::rng::run_cases;
use crate::support::{from_magv, to_magv};
use zenith_bound::mag::Mag;

const CASES: u64 = 3000;

fn is_finite(m: &MagV) -> bool {
    matches!(m, MagV::Finite { .. })
}
fn is_zero(m: &MagV) -> bool {
    matches!(m, MagV::Finite { mant: 0, .. })
}

fn check_binary(
    name: &str,
    lib: fn(&Mag, &Mag) -> Mag,
    oracle: fn(&MagV, &MagV) -> Option<MagV>,
    verify: fn(&MagV, &MagV, &MagV) -> Result<(), String>,
    needs_nonzero_divisor: bool,
) -> Result<(), String> {
    run_cases(name, CASES, &|rng| {
        let a = draw_mag(rng);
        let b = draw_partner(rng, &a);
        let got = to_magv(&lib(&from_magv(&a)?, &from_magv(&b)?))?;
        let want =
            oracle(&a, &b).ok_or_else(|| format!("direct oracle failed on {:?} {:?}", a, b))?;
        if got != want {
            return Err(format!(
                "{}({:?}, {:?}): library gave {:?}, integer oracle gave {:?}",
                name, a, b, got, want
            ));
        }
        if is_finite(&a) && is_finite(&b) && !(needs_nonzero_divisor && is_zero(&b)) {
            verify(&a, &b, &got)
                .map_err(|e| format!("{}({:?}, {:?}) = {:?}: model: {}", name, a, b, got, e))?;
        }
        Ok(())
    })
}

fn check_unary(
    name: &str,
    lib: fn(&Mag) -> Mag,
    oracle: fn(&MagV) -> Option<MagV>,
    verify: fn(&MagV, &MagV) -> Result<(), String>,
) -> Result<(), String> {
    run_cases(name, CASES, &|rng| {
        let a = draw_mag(rng);
        let got = to_magv(&lib(&from_magv(&a)?))?;
        let want = oracle(&a).ok_or_else(|| format!("direct oracle failed on {:?}", a))?;
        if got != want {
            return Err(format!(
                "{}({:?}): library gave {:?}, integer oracle gave {:?}",
                name, a, got, want
            ));
        }
        if is_finite(&a) {
            verify(&a, &got).map_err(|e| format!("{}({:?}) = {:?}: model: {}", name, a, got, e))?;
        }
        Ok(())
    })
}

fn commutes(name: &str, lib: fn(&Mag, &Mag) -> Mag) -> Result<(), String> {
    run_cases(name, CASES, &|rng| {
        let a = draw_mag(rng);
        let b = draw_partner(rng, &a);
        let ab = to_magv(&lib(&from_magv(&a)?, &from_magv(&b)?))?;
        let ba = to_magv(&lib(&from_magv(&b)?, &from_magv(&a)?))?;
        if ab != ba {
            return Err(format!(
                "{}({:?}, {:?}) = {:?} but swapped gives {:?}",
                name, a, b, ab, ba
            ));
        }
        Ok(())
    })
}

fn up_is_lower_or_successor(
    name: &str,
    up: fn(&Mag, &Mag) -> Mag,
    lower: fn(&Mag, &Mag) -> Mag,
    needs_nonzero_divisor: bool,
) -> Result<(), String> {
    run_cases(name, CASES, &|rng| {
        let a = draw_mag(rng);
        let b = draw_partner(rng, &a);
        if !is_finite(&a) || !is_finite(&b) || (needs_nonzero_divisor && is_zero(&b)) {
            return Ok(());
        }
        let u = to_magv(&up(&from_magv(&a)?, &from_magv(&b)?))?;
        let l = to_magv(&lower(&from_magv(&a)?, &from_magv(&b)?))?;
        if u != l && u != succ(&l) {
            return Err(format!(
                "{}({:?}, {:?}): up {:?} is neither lower {:?} nor its successor",
                name, a, b, u, l
            ));
        }
        Ok(())
    })
}

fn sqrt_up_is_lower_or_successor() -> Result<(), String> {
    run_cases("sqrt_up_is_lower_or_successor", CASES, &|rng| {
        let a = draw_mag(rng);
        if !is_finite(&a) {
            return Ok(());
        }
        let u = to_magv(&Mag::sqrt_up(&from_magv(&a)?))?;
        let l = to_magv(&Mag::sqrt_lower(&from_magv(&a)?))?;
        if u != l && u != succ(&l) {
            return Err(format!(
                "sqrt_up({:?}): up {:?} is neither lower {:?} nor its successor",
                a, u, l
            ));
        }
        Ok(())
    })
}

fn monotone(name: &str, lib: fn(&Mag, &Mag) -> Mag) -> Result<(), String> {
    run_cases(name, CASES, &|rng| {
        let x = draw_mag(rng);
        let y = draw_partner(rng, &x);
        let c = draw_partner(rng, &x);
        // order x <= y by the integer oracle's order
        let (lo, hi) = if direct::order(&x, &y) == core::cmp::Ordering::Greater {
            (y, x)
        } else {
            (x, y)
        };
        let r_lo = to_magv(&lib(&from_magv(&lo)?, &from_magv(&c)?))?;
        let r_hi = to_magv(&lib(&from_magv(&hi)?, &from_magv(&c)?))?;
        if direct::order(&r_lo, &r_hi) == core::cmp::Ordering::Greater {
            return Err(format!("{}({:?}, {:?}) = {:?} exceeds {}({:?}, {:?}) = {:?} although the first operand is not larger", name, lo, c, r_lo, name, hi, c, r_hi));
        }
        Ok(())
    })
}

#[test]
fn prop_add_up_matches_integer_oracle_and_model() -> Result<(), String> {
    check_binary(
        "add_up",
        Mag::add_up,
        direct::add_up,
        crate::add_up::verify,
        false,
    )
}

#[test]
fn prop_mul_up_matches_integer_oracle_and_model() -> Result<(), String> {
    check_binary(
        "mul_up",
        Mag::mul_up,
        direct::mul_up,
        crate::mul_up::verify,
        false,
    )
}

#[test]
fn prop_div_up_matches_integer_oracle_and_model() -> Result<(), String> {
    check_binary(
        "div_up",
        Mag::div_up,
        direct::div_up,
        crate::div_up::verify,
        true,
    )
}

#[test]
fn prop_mul_lower_matches_integer_oracle_and_model() -> Result<(), String> {
    check_binary(
        "mul_lower",
        Mag::mul_lower,
        direct::mul_lower,
        crate::lower_mul::verify,
        false,
    )
}

#[test]
fn prop_div_lower_matches_integer_oracle_and_model() -> Result<(), String> {
    check_binary(
        "div_lower",
        Mag::div_lower,
        direct::div_lower,
        crate::lower_div::verify,
        true,
    )
}

#[test]
fn prop_sub_lower_matches_integer_oracle_and_model() -> Result<(), String> {
    check_binary(
        "sub_lower",
        Mag::sub_lower,
        direct::sub_lower,
        crate::lower_sub::verify,
        false,
    )
}

#[test]
fn prop_sqrt_up_matches_integer_oracle_and_model() -> Result<(), String> {
    check_unary(
        "sqrt_up",
        Mag::sqrt_up,
        direct::sqrt_up,
        crate::sqrt_up::verify,
    )
}

#[test]
fn prop_sqrt_lower_matches_integer_oracle_and_model() -> Result<(), String> {
    check_unary(
        "sqrt_lower",
        Mag::sqrt_lower,
        direct::sqrt_lower,
        crate::lower_sqrt::verify,
    )
}

#[test]
fn prop_add_up_is_commutative() -> Result<(), String> {
    commutes("add_up", Mag::add_up)
}

#[test]
fn prop_mul_up_is_commutative() -> Result<(), String> {
    commutes("mul_up", Mag::mul_up)
}

#[test]
fn prop_mul_up_is_lower_or_successor() -> Result<(), String> {
    up_is_lower_or_successor("mul", Mag::mul_up, Mag::mul_lower, false)
}

#[test]
fn prop_div_up_is_lower_or_successor() -> Result<(), String> {
    up_is_lower_or_successor("div", Mag::div_up, Mag::div_lower, true)
}

#[test]
fn prop_sqrt_up_is_lower_or_successor() -> Result<(), String> {
    sqrt_up_is_lower_or_successor()
}

#[test]
fn prop_add_up_is_monotone() -> Result<(), String> {
    monotone("add_up", Mag::add_up)
}

#[test]
fn prop_mul_up_is_monotone() -> Result<(), String> {
    monotone("mul_up", Mag::mul_up)
}
