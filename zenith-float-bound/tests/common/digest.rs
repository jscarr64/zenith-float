//! Deterministic digests for the cross-target grid job (S13, TEST_PLAN P10 and section 13): FNV-1a 64 over bytes, every integer folded as
//! its LITTLE-ENDIAN bytes, so the digest does not depend on the target's endianness or pointer width. Integer-only, no dependencies.
//! A digest names a drift, not a security property. The result rows `NAME WORDS HEX16` are what tools/digests.sh compares across targets.

use super::modular::mul_u64;

const START: u64 = 0xCBF2_9CE4_8422_2325;
const PRIME: u64 = 0x0000_0100_0000_01B3;

/// A running digest and the number of folded bytes.
pub struct Digest {
    state: u64,
    words: u128,
}

impl Digest {
    pub fn start() -> Digest {
        Digest {
            state: START,
            words: 0,
        }
    }

    /// Folds bytes, FNV-1a: `state = (state ^ byte) * prime mod 2^64` per byte (`modular::mul_u64`).
    pub fn fold_bytes(&mut self, bytes: &[u8]) {
        for b in bytes {
            self.state = mul_u64(self.state ^ u64::from(*b), PRIME);
            self.words += 1;
        }
    }

    /// Folds a `u64` as its 8 little-endian bytes.
    pub fn fold_u64(&mut self, v: u64) {
        self.fold_bytes(&v.to_le_bytes());
    }

    /// Folds an `i64` as the little-endian bytes of its two's complement.
    pub fn fold_i64(&mut self, v: i64) {
        self.fold_bytes(&v.to_le_bytes());
    }

    pub fn value(&self) -> u64 {
        self.state
    }

    pub fn count(&self) -> u128 {
        self.words
    }

    /// The result row `NAME COUNT HEX16` (lowercase hex, 16 digits).
    pub fn row(&self, name: &str) -> String {
        format!("{} {} {:016x}\n", name, self.words, self.state)
    }
}
