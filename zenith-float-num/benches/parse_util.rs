//! Parse helpers for arithmetic and transcendental benches.

use zenith_float_num::{Consts, ExactNum, Radix, RoundingMode};

pub fn parse_list(p: usize, cc: &mut Consts, literals: &[&str]) -> Vec<ExactNum> {
    literals
        .iter()
        .map(|s| ExactNum::parse(s, Radix::Dec, p, RoundingMode::ToEven, cc))
        .collect()
}

pub fn cycle_batch(values: &[ExactNum], batch: usize) -> Vec<ExactNum> {
    let n = values.len();
    (0..batch).map(|i| values[i % n].clone()).collect()
}
