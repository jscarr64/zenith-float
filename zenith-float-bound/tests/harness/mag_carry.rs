//! TEST_PLAN 8.1 "Carry across the exponent" (L0a-2): `add_up` where the mantissa sum, or its round-up, reaches `2^64` renormalises to
//! `mant = 2^63`, `exp + 1`, including at the zf exponent bound `ZF_EXPONENT_MAX` (whose `+ 1` is still a `Mag` exponent) and at
//! `exp = i64::MAX`, where the carry overflows to `Infinity` (CONTRACT 3.3: a round-up overflow is `Infinity`, never a saturated finite value;
//! also 1.3, 14.3). The saturating `add_up` is mutant 179-C-add-nocarry (tools/mutants), which this file must kill. Oracle: hand-derived exact values. A `Mag` is
//! `mant * 2^(exp - 64)`; `u64::MAX` at `e` plus `2^63` at `e - 63` (one unit of the first) is exactly `2^e` = `2^63` at `e + 1`; with
//! the second operand at `e - 100` instead the exact sum lies strictly between `u64::MAX` at `e` and `2^e`, so rounding up gives `2^e`;
//! `u64::MAX` at `e` doubled is exactly `u64::MAX` at `e + 1` (no carry, the control).
//! COVERS: Mag::add_up

use crate::common::expect_eq;
use crate::common::generators::{TOP, ZF_EXPONENT_MAX};
use zenith_bound::mag::Mag;

/// The exponents the rows run at: an ordinary one, one below the zf bound, the zf bound, and the top of `i64`.
const CARRY_EXPONENTS: [i64; 5] = [
    0,
    ZF_EXPONENT_MAX - 1,
    ZF_EXPONENT_MAX,
    i64::MAX - 1,
    i64::MAX,
];

fn mag(mant: u64, exp: i64) -> Result<Mag, String> {
    Mag::finite(mant, exp).map_err(|e| format!("Mag::finite({:#x}, {}) failed: {:?}", mant, exp, e))
}

/// `mant` at `e + 1`, or `Infinity` when `e + 1` is past `i64::MAX`.
fn at_next(mant: u64, e: i64) -> Result<Mag, String> {
    match e.checked_add(1) {
        Some(n) => mag(mant, n),
        None => Ok(Mag::infinity()),
    }
}

/// Runs every carry row through `add` (both argument orders); returns the number of rows checked.
fn check_carry_rows(add: fn(&Mag, &Mag) -> Mag) -> Result<usize, String> {
    let mut rows = 0usize;
    for e in CARRY_EXPONENTS {
        let a = mag(u64::MAX, e)?;
        let cases = [
            ("exact carry", mag(TOP, e - 63)?, at_next(TOP, e)?),
            ("round-up carry", mag(TOP, e - 100)?, at_next(TOP, e)?),
            (
                "no carry (doubling)",
                mag(u64::MAX, e)?,
                at_next(u64::MAX, e)?,
            ),
        ];
        for (name, b, want) in cases {
            expect_eq(&format!("{} at {}: a + b", name, e), &add(&a, &b), &want)?;
            expect_eq(&format!("{} at {}: b + a", name, e), &add(&b, &a), &want)?;
            rows += 1;
        }
    }
    Ok(rows)
}

#[test]
fn add_up_carries_into_the_exponent_and_overflows_to_infinity_at_i64_max() -> Result<(), String> {
    expect_eq("carry rows", &check_carry_rows(Mag::add_up)?, &15)
}
