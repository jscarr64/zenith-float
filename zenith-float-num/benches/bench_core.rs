//! Helpers used only by the composite Criterion binary.

use std::hint::black_box;

use zenith_float_num::Consts;

/// Precision tiers shared with the other benches.
pub const PRECISIONS: &[usize] = &[128, 1024, 4096, 32768];

/// `None` when the constants cache cannot be built; the benchmark that asked for it then does nothing.
pub fn init_cc() -> Option<Consts> {
    Consts::new().ok()
}

/// Prevent LLVM from eliminating benchmark results.
pub fn sink<T>(v: T) {
    black_box(v);
}
