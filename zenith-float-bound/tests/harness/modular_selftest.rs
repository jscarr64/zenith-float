//! COVERS: common::modular::add_i64, common::modular::add_u64, common::modular::mul_u64, common::modular::sub_u64
//! Self-tests of tests/common/modular.rs against hand-derived values (mod 2^64) and against `checked_*` where no wrap happens.

use crate::common::expect_eq;
use crate::common::modular::{add_i64, add_u64, mul_u64, sub_u64};

#[test]
fn add_u64_wraps_at_two_to_the_64() -> Result<(), String> {
    expect_eq("no wrap", &add_u64(2, 3), &5)?;
    expect_eq("max + 1", &add_u64(u64::MAX, 1), &0)?;
    expect_eq("max + max", &add_u64(u64::MAX, u64::MAX), &(u64::MAX - 1))?;
    expect_eq(
        "matches checked_add",
        &add_u64(1 << 62, 1 << 62),
        &(1u64 << 63),
    )
}

#[test]
fn sub_u64_wraps_below_zero() -> Result<(), String> {
    expect_eq("no wrap", &sub_u64(5, 3), &2)?;
    expect_eq("0 - 1", &sub_u64(0, 1), &u64::MAX)?;
    expect_eq("0 - max", &sub_u64(0, u64::MAX), &1)?;
    expect_eq("x - x", &sub_u64(u64::MAX, u64::MAX), &0)
}

#[test]
fn mul_u64_keeps_the_low_64_bits() -> Result<(), String> {
    expect_eq("small", &mul_u64(6, 7), &42)?;
    expect_eq("2^32 * 2^32", &mul_u64(1 << 32, 1 << 32), &0)?;
    expect_eq("max * max = 1 mod 2^64", &mul_u64(u64::MAX, u64::MAX), &1)?;
    // the splitmix64 constant squared, derived by hand: (2^64 - k)^2 mod 2^64 = k^2 mod 2^64 with k = 2^64 - 0x9E3779B97F4A7C15
    let k: u128 = (1u128 << 64) - 0x9E37_79B9_7F4A_7C15;
    let want = u64::try_from((k * k) & u128::from(u64::MAX)).map_err(|_| String::from("mask"))?;
    expect_eq(
        "constant squared",
        &mul_u64(0x9E37_79B9_7F4A_7C15, 0x9E37_79B9_7F4A_7C15),
        &want,
    )
}

#[test]
fn add_i64_is_twos_complement() -> Result<(), String> {
    expect_eq("no wrap", &add_i64(-2, 5), &3)?;
    expect_eq("max + 1", &add_i64(i64::MAX, 1), &i64::MIN)?;
    expect_eq("min + -1", &add_i64(i64::MIN, -1), &i64::MAX)?;
    expect_eq("min + min", &add_i64(i64::MIN, i64::MIN), &0)
}
