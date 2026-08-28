//! Deterministic fixtures and helpers for zenith-float Criterion benches.
//!
//! All inputs are fixed decimal literals (constants, angles, coefficients) — no
//! random generation and no brute-force sweeps.

#![allow(dead_code)]

use std::hint::black_box;

use zenith_float_num::{Consts, ExactNum, Radix, RoundingMode};

/// Batch size for cheap arithmetic kernels (add / mul / div / rem).
pub const ARITH_BATCH: usize = 256;

/// Batch size for transcendental kernels (each call is substantially more work).
pub const TRANSCENDENTAL_BATCH: usize = 64;

/// Precision tiers: interactive, high-precision scientific, heavy Accumath-style,
/// and a tier large enough to exercise Toom-3 multiplication (512 mantissa words).
pub const PRECISIONS: &[usize] = &[128, 1024, 4096, 32768];

/// Both mantissa operands must exceed this word count to use FFT (`mantissa/mul.rs`).
pub const FFT_MUL_THRESHOLD_WORDS: usize = 5400;

/// Bit precision one word past the Toom-3 → FFT crossover (5401 mantissa words).
pub const FFT_MUL_PRECISION: usize = (FFT_MUL_THRESHOLD_WORDS + 1) * 64;

/// Heavier FFT-regime tier (8192 mantissa words).
pub const FFT_MUL_PRECISION_LARGE: usize = 8192 * 64;

pub const FFT_MUL_PRECISIONS: &[usize] = &[FFT_MUL_PRECISION, FFT_MUL_PRECISION_LARGE];

/// Small batch: each FFT-scale multiply is very expensive.
pub const FFT_MUL_BATCH: usize = 8;

/// Mathematical constants and coefficients that appear in real formula libraries.
pub const DECIMAL_FIXTURES: &[&str] = &[
    "3.141592653589793238462643383279",
    "2.718281828459045235360287471",
    "1.414213562373095048801688724",
    "0.577215664901532860606512090",
    "1.618033988749894848204586834",
    "1.234567890123456789",
    "9.869604401089358618834464944",
    "0.693147180559945309417232121",
    "2.302585092994045684017991455",
    "0.367879441171442321595523770",
    "12.345678901234567890",
    "0.001234567890123456789",
    "987.654321098765432109",
    "0.999999999999999999",
    "1.000000000000000001",
    "42.0",
];

/// Positive operands for logarithms (domain is (0, +Inf)).
pub const LN_FIXTURES: &[&str] = &[
    "1.234567890123456789",
    "2.718281828459045235360287471",
    "10.0",
    "0.001234567890123456789",
    "1234.567890123456789",
    "3.141592653589793238462643383279",
    "1.000000000000000001",
    "0.999999999999999999",
    "987.654321098765432109",
    "42.0",
];

/// Arguments where `exp` stays in range without premature overflow at modest precision.
pub const EXP_FIXTURES: &[&str] = &[
    "0.01",
    "0.1",
    "0.5",
    "1.0",
    "2.3",
    "-0.5",
    "-1.0",
    "0.693147180559945309417232121",
    "2.302585092994045684017991455",
    "3.141592653589793238462643383279",
];

/// Radian arguments for sin / cos / tan (typical after range reduction).
pub const TRIG_FIXTURES: &[&str] = &[
    "0.1",
    "0.5",
    "1.0",
    "1.234567890123456789",
    "2.718281828459045235360287471",
    "3.141592653589793238462643383279",
    "-0.75",
    "-2.3",
    "0.999999999999999999",
    "1.570796326794896619231321691",
];

pub fn init_cc() -> Consts {
    Consts::new().expect("constants cache")
}

pub fn parse_list(p: usize, cc: &mut Consts, literals: &[&str]) -> Vec<ExactNum> {
    literals
        .iter()
        .map(|s| ExactNum::parse(s, Radix::Dec, p, RoundingMode::ToEven, cc))
        .collect()
}

pub fn parse_fixtures(p: usize, cc: &mut Consts) -> Vec<ExactNum> {
    parse_list(p, cc, DECIMAL_FIXTURES)
}

pub fn pair_fixtures(p: usize, cc: &mut Consts, batch: usize) -> Vec<(ExactNum, ExactNum)> {
    let values = parse_fixtures(p, cc);
    let n = values.len();
    (0..batch)
        .map(|i| {
            (
                values[i % n].clone(),
                values[(i.wrapping_mul(7).wrapping_add(3)) % n].clone(),
            )
        })
        .collect()
}

pub fn cycle_batch(values: &[ExactNum], batch: usize) -> Vec<ExactNum> {
    let n = values.len();
    (0..batch)
        .map(|i| values[i % n].clone())
        .collect()
}

/// Prevent LLVM from eliminating benchmark results.
pub fn sink<T>(v: T) {
    black_box(v);
}
