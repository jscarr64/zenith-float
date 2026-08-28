//! Helpers shared by the arithmetic and transcendental Criterion binaries.

use std::hint::black_box;

use zenith_float_num::{Consts, ExactNum, Radix, RoundingMode};

/// Precision tiers: interactive, high-precision scientific, heavy Accumath-style,
/// and a tier large enough to exercise Toom-3 multiplication (512 mantissa words).
pub const PRECISIONS: &[usize] = &[128, 1024, 4096, 32768];

pub fn init_cc() -> Consts {
    Consts::new().expect("constants cache")
}

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

/// Prevent LLVM from eliminating benchmark results.
pub fn sink<T>(v: T) {
    black_box(v);
}
