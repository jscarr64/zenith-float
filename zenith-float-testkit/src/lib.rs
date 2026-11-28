//! Test-only helpers shared by the zenith-float unit and integration tests. A dev-dependency only;
//! never published. Its items are `pub`, so a test binary that uses part of them has no dead code.
#![forbid(unsafe_code)]
#![deny(unused)]

pub mod test_time;
