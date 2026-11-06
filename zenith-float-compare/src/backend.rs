//! Shared backend trait for cross-library comparison.

use std::time::Duration;

/// Benchmark setup / dispatch failure (replaces panics: every error returns an answer).
#[derive(Debug)]
pub enum BenchError {
    /// The task name is not one of `task::TASKS`.
    UnknownTask(String),
    /// The zenith constants cache could not be built.
    Constants,
}

impl std::fmt::Display for BenchError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BenchError::UnknownTask(t) => write!(f, "unknown task: {t}"),
            BenchError::Constants => write!(f, "constants cache could not be built"),
        }
    }
}

impl std::error::Error for BenchError {}

pub trait BenchState: Clone {}

pub trait BenchFloat: Sized + Clone {
    type State: BenchState;

    fn state(precision: usize) -> Result<Self::State, BenchError>;

    fn rand_normal(
        count: usize,
        exp_from: i32,
        exp_to: i32,
        state: &Self::State,
        sign_positive: bool,
    ) -> Vec<Self>;

    fn run_task(
        state: &Self::State,
        task: &str,
        values: &[Self],
    ) -> Result<(Self, Duration), BenchError>;
}

/// Convert base-10 exponent bounds to library exponent units (bigfloat-bench convention).
pub fn convert_exp(exp: i32) -> i32 {
    (exp as i64 * 3_321_928_095 / 1_000_000_000) as i32
}
