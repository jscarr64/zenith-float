//! S2/S4 deterministic greedy shrinker (TEST_PLAN section 7). No randomness, a fixed step cap, so a shrunk case reproduces.
//! Steps implemented (those possible without balls): reduce the precision along the schedule `PS`; translate all exponents toward 0 by
//! halving (valid because of W1: exponent translation preserves the result); replace mantissas by simpler patterns (top bit only, then all ones).
//! Plan text for later items, NOT implemented here: drop ball radii to 0 then to 1 ulp (needs L0b balls), and delete steps of a chain (needs chains).
//!
//! Cases are owned values: the shrinker takes the starting case, builds new cases with the `with_*` methods and passes them to the
//! `fails` closure by reference. Nothing is copied or cloned.

use super::generators::{PS, TOP};

/// What the shrinker needs to know about a case.
pub trait ShrinkCase: Sized {
    fn exponents(&self) -> Vec<i64>;
    fn with_exponents(&self, exponents: &[i64]) -> Self;
    fn prec(&self) -> u32;
    fn with_prec(&self, prec: u32) -> Self;
    fn mantissas(&self) -> Vec<u64>;
    fn with_mantissas(&self, mantissas: &[u64]) -> Self;
}

/// The default step cap: the most `fails` evaluations one shrink may spend.
pub const SHRINK_STEP_CAP: usize = 10_000;

/// A shrunk case and the number of `fails` evaluations it took.
pub struct Shrunk<C> {
    pub case: C,
    pub evaluations: usize,
}

/// Simplicity rank of a mantissa: top bit only (0), all ones (1), anything else (2).
fn mantissa_rank(m: u64) -> u8 {
    if m == TOP {
        0
    } else if m == u64::MAX {
        1
    } else {
        2
    }
}

/// Every simpler neighbour of `case`, in the fixed order: precisions along `PS` ascending, exponent translations (shift by `e0`, `e0/2`, `e0/4`, ...
/// where `e0` is the first exponent), mantissa replacements by position (top bit only, then all ones, only when simpler than the current one).
fn candidates<C: ShrinkCase>(case: &C) -> Vec<C> {
    let mut out: Vec<C> = Vec::new();
    let prec = case.prec();
    for p in PS.iter() {
        if *p < prec {
            out.push(case.with_prec(*p));
        }
    }
    let exps = case.exponents();
    if let Some(e0) = exps.first() {
        let mut shift = *e0;
        while shift != 0 {
            let moved: Option<Vec<i64>> = exps.iter().map(|e| e.checked_sub(shift)).collect();
            if let Some(m) = moved {
                out.push(case.with_exponents(&m));
            }
            shift /= 2;
        }
    }
    let mants = case.mantissas();
    for (i, m) in mants.iter().enumerate() {
        for target in [TOP, u64::MAX] {
            if mantissa_rank(target) < mantissa_rank(*m) {
                let replaced: Vec<u64> = mants
                    .iter()
                    .enumerate()
                    .map(|(j, x)| if j == i { target } else { *x })
                    .collect();
                out.push(case.with_mantissas(&replaced));
            }
        }
    }
    out
}

/// Shrinks `start`, a case for which `fails` is true, greedily: take the first candidate (fixed order, see `candidates`) that still fails,
/// restart from it, stop when no candidate fails or `cap` evaluations are spent. `Err` when `start` does not fail.
pub fn shrink<C: ShrinkCase>(
    start: C,
    fails: &dyn Fn(&C) -> bool,
    cap: usize,
) -> Result<Shrunk<C>, String> {
    if !fails(&start) {
        return Err(String::from(
            "the starting case does not fail: nothing to shrink",
        ));
    }
    let mut current = start;
    let mut evaluations = 1usize;
    'outer: loop {
        for cand in candidates(&current) {
            if evaluations >= cap {
                break 'outer;
            }
            evaluations += 1;
            if fails(&cand) {
                current = cand;
                continue 'outer;
            }
        }
        break;
    }
    Ok(Shrunk {
        case: current,
        evaluations,
    })
}

/// True when no simpler neighbour of `case` still fails (a local minimum of the shrinker's moves).
pub fn is_locally_minimal<C: ShrinkCase>(case: &C, fails: &dyn Fn(&C) -> bool) -> bool {
    !candidates(case).iter().any(fails)
}
