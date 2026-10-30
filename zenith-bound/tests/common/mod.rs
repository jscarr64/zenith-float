//! Shared test support for zenith-bound (checklist item S1): deterministic generator and fixed seeds. Integer-only, no dependencies.
//!
//! RULE: tests/common is included ONLY by the single harness target (tests/harness/main.rs), as
//! `#[path = "../common/mod.rs"] mod common;`. No other test target may include it. The crate denies `unused` and `#[allow]` is banned,
//! so every `pub` item added to common MUST be used by the harness target (its self-tests count), or the build fails with a
//! `never used` error. Add the item and its use in the same commit.
pub mod ballrows;
pub mod digest;
pub mod edges;
pub mod generators;
pub mod goldens;
pub mod modular;
pub mod points;
pub mod replay;
pub mod rng;
pub mod seeds;
pub mod shrink;
pub mod tightness;

use core::fmt::Debug;

/// A failing check returns `Err` (never a panic): `label` names the check, `got` and `want` are printed with `{:?}`.
pub fn expect_eq<T: PartialEq + Debug>(label: &str, got: &T, want: &T) -> Result<(), String> {
    if got == want {
        Ok(())
    } else {
        Err(format!("{}: expected {:?}, got {:?}", label, want, got))
    }
}
