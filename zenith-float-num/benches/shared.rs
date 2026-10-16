//! Helpers shared by the arithmetic and transcendental Criterion binaries.

use std::hint::black_box;

use zenith_float_num::{Consts, ExactNum, Radix, RoundingMode};

/// Precision tiers: interactive, high-precision scientific, heavy evaluation,
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

/// The literals `literals[(i * mul + add) % n]` for `i` in `0..batch`.
pub fn cycle_literals<'a>(
    literals: &[&'a str],
    batch: usize,
    mul: usize,
    add: usize,
) -> Vec<&'a str> {
    let n = literals.len();
    (0..batch)
        .map(|i| literals[(i.wrapping_mul(mul).wrapping_add(add)) % n])
        .collect()
}

/// `batch` operands parsed from the literals in cyclic order (each is parsed afresh).
pub fn cycle_batch(p: usize, cc: &mut Consts, literals: &[&str], batch: usize) -> Vec<ExactNum> {
    parse_list(p, cc, &cycle_literals(literals, batch, 1, 0))
}

/// `batch` operand pairs parsed from the literals: the first walks the list in order, the second
/// at index `(i * mul + add) % n`.
pub fn cycle_pairs(
    p: usize,
    cc: &mut Consts,
    literals: &[&str],
    batch: usize,
    mul: usize,
    add: usize,
) -> Vec<(ExactNum, ExactNum)> {
    let firsts = parse_list(p, cc, &cycle_literals(literals, batch, 1, 0));
    let seconds = parse_list(p, cc, &cycle_literals(literals, batch, mul, add));
    firsts.into_iter().zip(seconds).collect()
}

/// Prevent LLVM from eliminating benchmark results.
pub fn sink<T>(v: T) {
    black_box(v);
}
