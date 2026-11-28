//! splitmix64 (Steele, Lea, Flood): one `u64` state, integer-only, all arithmetic modulo 2^64 through `modular` (no `wrapping_*`).
//! The draw order is frozen; changing it changes every case of every suite and bumps `seeds::GENERATOR_VERSION`.

use super::modular::{add_u64, mul_u64};

const INCREMENT: u64 = 0x9E37_79B9_7F4A_7C15;
const MUL1: u64 = 0xBF58_476D_1CE4_E5B9;
const MUL2: u64 = 0x94D0_49BB_1331_11EB;

pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    pub const fn new(seed: u64) -> Self {
        SplitMix64 { state: seed }
    }

    pub fn next_u64(&mut self) -> u64 {
        self.state = add_u64(self.state, INCREMENT);
        let mut z = self.state;
        z = mul_u64(z ^ (z >> 30), MUL1);
        z = mul_u64(z ^ (z >> 27), MUL2);
        z ^ (z >> 31)
    }

    /// A value in `0..n` (`0` when `n == 0`): one draw, widening multiply (no modulo bias to speak of, no division).
    pub fn below(&mut self, n: u64) -> u64 {
        let wide = (u128::from(self.next_u64()) * u128::from(n)) >> 64;
        u64::try_from(wide).unwrap_or(0)
    }
}

/// The first output of a generator seeded with `seed` (one-shot mixer, used to derive stream seeds).
pub fn splitmix64(seed: u64) -> u64 {
    SplitMix64::new(seed).next_u64()
}
