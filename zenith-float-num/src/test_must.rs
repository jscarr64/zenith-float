//! `?` for tests (P2, Jeff: no panic in test code either; tests return `Result` and use `?`).
//!
//! `x.must("what")?` turns a `None` or an `Err(e)` into `Err(String)` naming what was expected, so a test reports the
//! failure instead of panicking. Shared by the library's unit tests and `tests/tier1_regressions.rs` (by `#[path]`).

use alloc::format;
use alloc::string::String;
use core::fmt::Debug;

/// A value a test needs: `Some` / `Ok`, else an `Err(String)` naming it.
pub(crate) trait Must<T> {
    /// `Ok(value)`, or `Err("<what>: None")` / `Err("<what>: <error>")`.
    fn must(self, what: &str) -> Result<T, String>;
}

impl<T> Must<T> for Option<T> {
    fn must(self, what: &str) -> Result<T, String> {
        self.ok_or_else(|| format!("{what}: None"))
    }
}

impl<T, E: Debug> Must<T> for Result<T, E> {
    fn must(self, what: &str) -> Result<T, String> {
        self.map_err(|e| format!("{what}: {e:?}"))
    }
}
