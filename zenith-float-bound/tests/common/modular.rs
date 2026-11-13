//! Arithmetic modulo 2^64 without `wrapping_*` (rule: checked arithmetic only, TEST_PLAN section 3 / Spence's review of 7cbf6dc).
//! Where a hash or a generator needs the wrap (splitmix64, FNV-1a, a seed pack), it is stated explicitly: the exact result is computed in `u128`
//! (it cannot overflow: both operands are below 2^64) and the low 64 bits are taken. Integer-only.

fn low64(x: u128) -> u64 {
    let [_, _, _, _, _, _, _, _, a, b, c, d, e, f, g, h] = x.to_be_bytes();
    u64::from_be_bytes([a, b, c, d, e, f, g, h])
}

/// `(a + b) mod 2^64`.
pub fn add_u64(a: u64, b: u64) -> u64 {
    low64(u128::from(a) + u128::from(b))
}

/// `(a - b) mod 2^64`.
pub fn sub_u64(a: u64, b: u64) -> u64 {
    // a + 2^64 - b is at least 1 and below 2^65
    low64(u128::from(a) + (1u128 << 64) - u128::from(b))
}

/// `(a * b) mod 2^64`.
pub fn mul_u64(a: u64, b: u64) -> u64 {
    low64(u128::from(a) * u128::from(b))
}

/// `(a + b)` in two's complement, i.e. modulo 2^64 reinterpreted as `i64`.
pub fn add_i64(a: i64, b: i64) -> i64 {
    i64::from_ne_bytes(
        add_u64(
            u64::from_ne_bytes(a.to_ne_bytes()),
            u64::from_ne_bytes(b.to_ne_bytes()),
        )
        .to_ne_bytes(),
    )
}
