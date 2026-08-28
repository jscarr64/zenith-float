//! Helpers used by every Criterion binary in this crate.

use std::hint::black_box;

use zenith_float_num::Consts;

/// Precision tiers: interactive, high-precision scientific, heavy Accumath-style,
/// and a tier large enough to exercise Toom-3 multiplication (512 mantissa words).
pub const PRECISIONS: &[usize] = &[128, 1024, 4096, 32768];

pub fn init_cc() -> Consts {
    Consts::new().expect("constants cache")
}

/// Prevent LLVM from eliminating benchmark results.
pub fn sink<T>(v: T) {
    black_box(v);
}
