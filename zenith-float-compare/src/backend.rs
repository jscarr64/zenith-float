//! Shared backend trait for cross-library comparison.

use std::time::Duration;

pub trait BenchState: Clone {}

pub trait BenchFloat: Sized + Clone {
    type State: BenchState;

    fn state(precision: usize) -> Self::State;

    fn rand_normal(
        count: usize,
        exp_from: i32,
        exp_to: i32,
        state: &Self::State,
        sign_positive: bool,
    ) -> Vec<Self>;

    fn run_task(state: &Self::State, task: &str, values: &[Self]) -> (Self, Duration);
}

/// Convert base-10 exponent bounds to library exponent units (bigfloat-bench convention).
pub fn convert_exp(exp: i32) -> i32 {
    (exp as i64 * 3_321_928_095 / 1_000_000_000) as i32
}
