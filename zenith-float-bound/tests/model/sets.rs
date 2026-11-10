//! Operand sets shared by all op tests: mantissas, exponent grids and the xorshift generator that seeds the random mantissas. Zero is a member of every set.

use super::dyadic::MagV;
use super::oracle::TWO63;

pub struct Rng {
    state: u64,
}

impl Rng {
    pub fn new(seed: u64) -> Rng {
        Rng { state: seed | 1 }
    }

    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.state = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
}

pub fn mantissas() -> Vec<u64> {
    let mut result = vec![
        TWO63,
        TWO63 + 1,
        TWO63 + 2,
        0xC000_0000_0000_0000,
        0xAAAA_AAAA_AAAA_AAAB,
        u64::MAX - 1,
        u64::MAX,
    ];
    let mut rng = Rng::new(0x5EED);
    for _ in 0..5 {
        result.push(rng.next_u64() | TWO63);
    }
    result
}

pub fn mid_exponents() -> Vec<i64> {
    vec![-200, -130, -65, -64, -63, -1, 0, 1, 63, 64, 65, 130, 200]
}

pub fn extreme_exponents() -> Vec<i64> {
    vec![
        i64::MIN,
        i64::MIN + 1,
        i64::MIN + 64,
        i64::MIN + 65,
        -1,
        0,
        1,
        i64::MAX - 65,
        i64::MAX - 64,
        i64::MAX - 1,
        i64::MAX,
    ]
}

pub fn magv_set(exponents: &[i64]) -> Vec<MagV> {
    let mut result = vec![MagV::Finite { exp: 0, mant: 0 }];
    let mants = mantissas();
    for &exp in exponents.iter() {
        for &mant in mants.iter() {
            result.push(MagV::Finite { exp, mant });
        }
    }
    result
}
