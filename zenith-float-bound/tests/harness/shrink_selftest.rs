//! COVERS: common::shrink::is_locally_minimal, common::shrink::shrink
//! Self-tests of tests/common/shrink.rs on a synthetic case type. Oracle: minimal cases derived by hand from the stated predicates
//! (see each test). The evaluation counts are pinned from this implementation and only guard determinism and the step cap.

use crate::common::expect_eq;
use crate::common::generators::{PS, TOP};
use crate::common::shrink::{is_locally_minimal, shrink, ShrinkCase, Shrunk, SHRINK_STEP_CAP};

#[derive(Debug, PartialEq, Eq)]
struct Syn {
    exps: [i64; 2],
    prec: u32,
    mants: [u64; 2],
}

fn syn(e0: i64, e1: i64, prec: u32, m0: u64, m1: u64) -> Syn {
    Syn {
        exps: [e0, e1],
        prec,
        mants: [m0, m1],
    }
}

impl ShrinkCase for Syn {
    fn exponents(&self) -> Vec<i64> {
        Vec::from(self.exps)
    }
    fn with_exponents(&self, exponents: &[i64]) -> Self {
        Syn {
            exps: [
                exponents.first().map_or(0, |&x| x),
                exponents.get(1).map_or(0, |&x| x),
            ],
            prec: self.prec,
            mants: [
                self.mants.first().map_or(0, |&x| x),
                self.mants.get(1).map_or(0, |&x| x),
            ],
        }
    }
    fn prec(&self) -> u32 {
        self.prec
    }
    fn with_prec(&self, prec: u32) -> Self {
        Syn {
            exps: [
                self.exps.first().map_or(0, |&x| x),
                self.exps.get(1).map_or(0, |&x| x),
            ],
            prec,
            mants: [
                self.mants.first().map_or(0, |&x| x),
                self.mants.get(1).map_or(0, |&x| x),
            ],
        }
    }
    fn mantissas(&self) -> Vec<u64> {
        Vec::from(self.mants)
    }
    fn with_mantissas(&self, mantissas: &[u64]) -> Self {
        Syn {
            exps: [
                self.exps.first().map_or(0, |&x| x),
                self.exps.get(1).map_or(0, |&x| x),
            ],
            prec: self.prec,
            mants: [
                mantissas.first().map_or(0, |&x| x),
                mantissas.get(1).map_or(0, |&x| x),
            ],
        }
    }
}

/// Fails iff the exponent gap is at least 70, the precision at least 5, and the first mantissa has its top bit.
fn gap_fails(c: &Syn) -> bool {
    let e0 = c.exps.first().map_or(0, |&x| x);
    let e1 = c.exps.get(1).map_or(0, |&x| x);
    let m0 = c.mants.first().map_or(0, |&x| x);
    e1.abs_diff(e0) >= 70 && c.prec >= 5 && m0 & TOP != 0
}

/// Fails iff the first exponent is at least 1000.
fn high_fails(c: &Syn) -> bool {
    c.exps.first().map_or(0, |&x| x) >= 1000
}

#[test]
fn the_shrinker_reaches_the_known_minimal_case() -> Result<(), String> {
    // minimum by hand: translate to e0 = 0 (gap 123 kept), precision 5 is the smallest PS value >= 5, mantissas TOP
    let start = syn(
        1_000_000,
        1_000_123,
        1000,
        u64::MAX - 1,
        0xAAAA_AAAA_AAAA_AAAB,
    );
    let out = shrink(start, &gap_fails, SHRINK_STEP_CAP)?;
    expect_eq("minimal case", &out.case, &syn(0, 123, 5, TOP, TOP))?;
    if !gap_fails(&out.case) || !is_locally_minimal(&out.case, &gap_fails) {
        return Err(String::from(
            "the result does not fail or is not locally minimal",
        ));
    }
    expect_eq("evaluations (pinned)", &out.evaluations, &EVALS_GAP)?;
    // a negative start translates toward 0 as well
    let negative = shrink(
        syn(-3_000_000_000, -2_999_999_900, 4096, TOP, TOP),
        &gap_fails,
        SHRINK_STEP_CAP,
    )?;
    expect_eq("negative start", &negative.case, &syn(0, 100, 5, TOP, TOP))?;
    // the minimal precision is a schedule value
    if !PS.contains(&out.case.prec) {
        return Err(String::from("the shrunk precision is not on the schedule"));
    }
    Ok(())
}

const EVALS_GAP: usize = 25;
const EVALS_HIGH: usize = 47;

#[test]
fn exponent_translation_halves_and_stops_at_the_boundary() -> Result<(), String> {
    // fails iff e0 >= 1000: halving shifts get as close as the smallest shift (1) allows, so the minimum is exactly 1000
    let out = shrink(
        syn(1_000_000, 5, 53, TOP, TOP),
        &high_fails,
        SHRINK_STEP_CAP,
    )?;
    expect_eq("e0 minimal", &out.case.exps.first(), &Some(&1000))?;
    expect_eq("prec reduced to 1", &out.case.prec, &1)?;
    expect_eq("evaluations (pinned)", &out.evaluations, &EVALS_HIGH)?;
    if !is_locally_minimal(&out.case, &high_fails) {
        return Err(String::from("not locally minimal"));
    }
    Ok(())
}

#[test]
fn the_shrinker_is_deterministic() -> Result<(), String> {
    let a = shrink(
        syn(
            77_777_777,
            77_777_999,
            1024,
            u64::MAX,
            0xC000_0000_0000_0000,
        ),
        &gap_fails,
        SHRINK_STEP_CAP,
    )?;
    let b = shrink(
        syn(
            77_777_777,
            77_777_999,
            1024,
            u64::MAX,
            0xC000_0000_0000_0000,
        ),
        &gap_fails,
        SHRINK_STEP_CAP,
    )?;
    expect_eq("same result", &a.case, &b.case)?;
    expect_eq("same evaluation count", &a.evaluations, &b.evaluations)?;
    // shrinking a shrunk case changes nothing
    let Shrunk { case, .. } = shrink(syn(0, 222, 5, TOP, TOP), &gap_fails, SHRINK_STEP_CAP)?;
    expect_eq("fixed point", &case, &syn(0, 222, 5, TOP, TOP))
}

#[test]
fn shrinking_a_passing_case_and_hitting_the_cap_are_reported() -> Result<(), String> {
    let never = |_: &Syn| false;
    if shrink(syn(1, 2, 3, TOP, TOP), &never, SHRINK_STEP_CAP).is_ok() {
        return Err(String::from("a start that does not fail was accepted"));
    }
    let capped = shrink(
        syn(1_000_000, 1_000_123, 1000, u64::MAX, u64::MAX),
        &gap_fails,
        3,
    )?;
    if capped.evaluations > 3 {
        return Err(format!(
            "{} evaluations with a cap of 3",
            capped.evaluations
        ));
    }
    if is_locally_minimal(&capped.case, &gap_fails) {
        return Err(String::from(
            "a case cut short by the cap must not be locally minimal",
        ));
    }
    Ok(())
}

#[test]
fn a_wrong_expected_minimum_or_a_non_minimal_case_is_detected() -> Result<(), String> {
    // mutation-style: a deliberately wrong expectation must fail the comparison
    let out = shrink(
        syn(1_000_000, 1_000_123, 1000, TOP, TOP),
        &gap_fails,
        SHRINK_STEP_CAP,
    )?;
    if expect_eq("wrong minimum", &out.case, &syn(0, 123, 6, TOP, TOP)).is_ok() {
        return Err(String::from("expect_eq accepted a wrong minimum"));
    }
    // a deliberately wrong predicate (fails for everything) shrinks everything to the floor, which differs from the real minimum
    let always = |_: &Syn| true;
    let floor = shrink(
        syn(1_000_000, 1_000_123, 1000, u64::MAX, u64::MAX),
        &always,
        SHRINK_STEP_CAP,
    )?;
    expect_eq(
        "floor under an always-failing predicate",
        &floor.case.prec,
        &1,
    )?;
    if floor.case == out.case {
        return Err(String::from("a wrong predicate produced the real minimum"));
    }
    // cases that can still be simplified are not locally minimal
    for c in [
        syn(1, 71, 5, TOP, TOP),
        syn(0, 123, 6, TOP, TOP),
        syn(0, 123, 5, u64::MAX, TOP),
    ] {
        if is_locally_minimal(&c, &gap_fails) {
            return Err(format!("{:?} was reported locally minimal", c));
        }
    }
    Ok(())
}
