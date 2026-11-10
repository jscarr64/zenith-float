//! COVERS: Mag::add_up, Mag::div_lower, Mag::div_up, Mag::half_ulp, Mag::mul_lower, Mag::mul_up, Mag::smaller_than_half_ulp, Mag::sqrt_lower, Mag::sqrt_up, Mag::sub_lower, Mag::ulp_mag
//! L0a-5 wide-range Mag property fuzz (TEST_PLAN 8.1, section 7, checklist L0a-5): one `#[test]` per property family, each run through
//! Harper's property-runner contract in tests/harness/fuzz_mag.rs (`run_property`: per-property stream, `ZB_FUZZ_SEED`/`ZB_FUZZ_CASES`/
//! `ZB_REPLAY`, `ZB_FUZZ_DEEP` deep stream, non-vacuity on outcome-class counts, `<replay line> -- <why>` on failure). The checkers and their
//! self-tests are in fuzz_mag_props.rs and fuzz_mag_props_selftest.rs. Gate count `PROP_CASES` per property (a count, never wall-clock time),
//! then every property's forced pairs (`fuzz_mag_props::run_forced`: the zf exponent bounds, 0 and the 1024 scale).

use crate::common::expect_eq;
use crate::common::rng::SplitMix64;
use crate::fuzz_mag::{run_property, CaseFailure};
use crate::fuzz_mag_props::{self as props, Bad, Mono, Ops, Src};

/// Cases per property at the gate tier (override with `ZB_FUZZ_CASES`).
const PROP_CASES: u64 = 1_000_000;

impl CaseFailure for Bad {
    fn op(&self) -> &str {
        self.op
    }
    fn prec(&self) -> usize {
        self.prec
    }
    fn a(&self) -> &str {
        self.a.as_str()
    }
    fn b(&self) -> &str {
        self.b.as_str()
    }
    fn why(&self) -> &str {
        self.msg.as_str()
    }
}

/// One property over the real library ops at `PROP_CASES` (environment overrides); `need` are the outcome classes that must occur.
fn go(
    prop_index: u64,
    label: &str,
    need: &[(u8, &str)],
    f: &dyn Fn(&Ops, &mut Src) -> Result<u8, Bad>,
) -> Result<(), String> {
    let ops = props::real();
    run_property(
        prop_index,
        label,
        need,
        PROP_CASES,
        &|rng: &mut SplitMix64| f(&ops, &mut Src::random(rng)),
    )?;
    let forced = props::run_forced(&|src| f(&ops, src))
        .map_err(|m| format!("property {} ({}): {}", prop_index, label, m))?;
    expect_eq("forced pairs run", &forced, &props::FORCED_PAIRS)
}

const ALL: [(u8, &str); 4] = [
    (props::FIN, "finite"),
    (props::INF, "Infinity"),
    (props::ZERO, "zero"),
    (props::EXACT, "exact"),
];

#[test]
fn fuzz_add_up_is_sound_and_the_smallest_upper_bound() -> Result<(), String> {
    go(0, "add_up", &[ALL[0], ALL[1]], &props::add_up_case)
}

#[test]
fn fuzz_sub_lower_is_zero_iff_not_above_and_the_largest_lower_bound() -> Result<(), String> {
    go(1, "sub_lower", &[ALL[0], ALL[2]], &props::sub_lower_case)
}

#[test]
fn fuzz_mul_up_and_mul_lower_are_tight_and_bracket_each_other() -> Result<(), String> {
    go(
        2,
        "mul_up mul_lower",
        &[ALL[0], ALL[1], ALL[2], ALL[3]],
        &props::mul_case,
    )
}

#[test]
fn fuzz_div_up_and_div_lower_are_tight_and_bracket_each_other() -> Result<(), String> {
    go(
        3,
        "div_up div_lower",
        &[ALL[0], ALL[1], ALL[2], ALL[3]],
        &props::div_case,
    )
}

#[test]
fn fuzz_sqrt_up_and_sqrt_lower_are_tight_and_bracket_each_other() -> Result<(), String> {
    go(
        4,
        "sqrt_up sqrt_lower",
        &[ALL[0], ALL[3]],
        &props::sqrt_case,
    )
}

#[test]
fn fuzz_add_up_mul_up_and_mul_lower_commute() -> Result<(), String> {
    go(
        5,
        "add_up mul_up mul_lower",
        &[ALL[0]],
        &props::commutative_case,
    )
}

#[test]
fn fuzz_add_up_is_monotone() -> Result<(), String> {
    go(6, "add_up", &[ALL[0]], &|ops, rng| {
        props::mono_case(ops, &Mono::AddUp, rng)
    })
}

#[test]
fn fuzz_mul_up_and_mul_lower_are_monotone() -> Result<(), String> {
    go(7, "mul_up mul_lower", &[ALL[0]], &|ops, rng| {
        props::mono_case(ops, &Mono::MulUp, rng)?;
        props::mono_case(ops, &Mono::MulLower, rng)
    })
}

#[test]
fn fuzz_div_up_and_div_lower_are_nondecreasing_in_the_numerator_and_nonincreasing_in_the_divisor(
) -> Result<(), String> {
    go(8, "div_up div_lower", &[ALL[0]], &|ops, rng| {
        props::mono_case(ops, &Mono::DivUp, rng)?;
        props::mono_case(ops, &Mono::DivLower, rng)
    })
}

#[test]
fn fuzz_sqrt_up_and_sqrt_lower_are_monotone() -> Result<(), String> {
    go(9, "sqrt_up sqrt_lower", &[ALL[0]], &|ops, rng| {
        props::mono_case(ops, &Mono::SqrtUp, rng)?;
        props::mono_case(ops, &Mono::SqrtLower, rng)
    })
}

#[test]
fn fuzz_sub_lower_is_nondecreasing_in_the_minuend_and_nonincreasing_in_the_subtrahend(
) -> Result<(), String> {
    go(10, "sub_lower", &[ALL[0]], &|ops, rng| {
        props::mono_case(ops, &Mono::SubLower, rng)
    })
}

#[test]
fn fuzz_identities_of_add_up_sub_lower_mul_up_mul_lower_div_up_div_lower_and_sqrt(
) -> Result<(), String> {
    go(
        11,
        "add_up sub_lower mul_up mul_lower div_up div_lower sqrt_up sqrt_lower",
        &[ALL[0]],
        &props::identity_case,
    )
}

#[test]
fn fuzz_sqrt_up_and_sqrt_lower_are_exact_on_squares_and_odd_powers_of_two() -> Result<(), String> {
    go(
        12,
        "sqrt_up sqrt_lower",
        &[ALL[0], ALL[3]],
        &props::sqrt_exact_case,
    )
}

#[test]
fn fuzz_every_result_is_canonical_and_zero_and_infinity_operands_give_the_specified_value(
) -> Result<(), String> {
    go(
        13,
        "add_up sub_lower mul_up mul_lower div_up div_lower sqrt_up sqrt_lower",
        &[ALL[0], ALL[1], ALL[2]],
        &props::canonical_case,
    )
}

#[test]
fn fuzz_ulp_mag_half_ulp_and_smaller_than_half_ulp_match_the_exact_values() -> Result<(), String> {
    go(
        14,
        "ulp_mag half_ulp smaller_than_half_ulp",
        &[ALL[0], ALL[1], ALL[2], ALL[3]],
        &props::ulp_case,
    )
}
