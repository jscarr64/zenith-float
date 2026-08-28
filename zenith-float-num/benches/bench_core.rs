//! Helpers used only by the composite Criterion binary.

use std::hint::black_box;

use zenith_float_num::Consts;

/// Precision tiers shared with the other benches.
pub const PRECISIONS: &[usize] = &[128, 1024, 4096, 32768];

pub fn init_cc() -> Consts {
    Consts::new().expect("constants cache")
}

/// Prevent LLVM from eliminating benchmark results.
pub fn sink<T>(v: T) {
    black_box(v);
}
