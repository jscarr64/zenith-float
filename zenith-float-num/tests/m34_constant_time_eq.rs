//! ZF-M34 red rows: narrowing casts of size-like values (constant_time_eq, via the crate-root export).

use zenith_float_num::constant_time_eq;

/// hash.rs:63 seeded the length mismatch as `(a.len() ^ b.len()) as u8`: lengths that differ by
/// a multiple of 256 lost it, and all-zero tails then compared equal.
#[test]
fn constant_time_eq_length_mismatch_survives_mod_256() {
    let z = [0u8; 512];
    assert!(!constant_time_eq(&z[..256], &[]));
    assert!(!constant_time_eq(&[], &z[..256]));
    assert!(!constant_time_eq(&z[..300], &z[..44]));
    assert!(!constant_time_eq(&z[..512], &z[..256]));
    assert!(constant_time_eq(&z[..256], &z[..256]));
    assert!(!constant_time_eq(&z[..1], &[]));
}
