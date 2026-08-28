//! Arithmetic-only fixtures (compiled only from `arithmetic.rs`).

use super::shared::parse_list;
use zenith_float_num::{Consts, ExactNum};

/// Batch size for cheap arithmetic kernels (add / mul / div / rem).
pub const ARITH_BATCH: usize = 256;

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
