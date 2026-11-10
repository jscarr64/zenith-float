//! `Nat`: a minimal arbitrary-precision natural number (little-endian `u64` limbs), just enough for exact add, mul, shift and compare. Integer-only; every operation is total.

use core::cmp::Ordering;

#[derive(Debug, PartialEq, Eq)]
pub struct Nat {
    limbs: Vec<u64>,
}

fn lo(x: u128) -> u64 {
    let [b0, b1, b2, b3, b4, b5, b6, b7, _, _, _, _, _, _, _, _] = x.to_le_bytes();
    u64::from_le_bytes([b0, b1, b2, b3, b4, b5, b6, b7])
}

fn hi(x: u128) -> u64 {
    let [_, _, _, _, _, _, _, _, b8, b9, b10, b11, b12, b13, b14, b15] = x.to_le_bytes();
    u64::from_le_bytes([b8, b9, b10, b11, b12, b13, b14, b15])
}

fn add_limb_at(limbs: &mut Vec<u64>, pos: usize, x: u64) {
    let mut pos = pos;
    let mut x = x;
    while x != 0 {
        while limbs.len() <= pos {
            limbs.push(0);
        }
        if let Some(slot) = limbs.get_mut(pos) {
            let s = u128::from(*slot) + u128::from(x);
            *slot = lo(s);
            x = hi(s);
        }
        pos += 1;
    }
}

fn add_u128_at(limbs: &mut Vec<u64>, pos: usize, v: u128) {
    add_limb_at(limbs, pos, lo(v));
    add_limb_at(limbs, pos + 1, hi(v));
}

impl Nat {
    pub fn zero() -> Nat {
        Nat { limbs: vec![] }
    }

    pub fn from_u64(x: u64) -> Nat {
        if x == 0 {
            Nat::zero()
        } else {
            Nat::normalize(vec![x])
        }
    }

    pub fn is_zero(&self) -> bool {
        self.limbs.is_empty()
    }

    pub fn bit_len(&self) -> u64 {
        let mut total: u64 = 0;
        let mut top_seen = false;
        for limb in self.limbs.iter().rev() {
            if !top_seen {
                let mut t = *limb;
                while t != 0 {
                    total += 1;
                    t >>= 1;
                }
                top_seen = true;
            } else {
                total += 64;
            }
        }
        total
    }

    pub fn shl(&self, bits: usize) -> Nat {
        if self.is_zero() {
            return Nat::zero();
        }
        let limb_shift = bits / 64;
        let mut limbs: Vec<u64> = vec![0; limb_shift];
        let mut factor: u128 = 1;
        for _ in 0..(bits % 64) {
            factor *= 2;
        }
        let mut carry: u128 = 0;
        for limb in self.limbs.iter() {
            let v = u128::from(*limb) * factor + carry;
            limbs.push(lo(v));
            carry = u128::from(hi(v));
        }
        if carry != 0 {
            limbs.push(lo(carry));
        }
        Nat::normalize(limbs)
    }

    pub fn add(&self, other: &Nat) -> Nat {
        let max_len = self.limbs.len().max(other.limbs.len());
        let mut carry = 0u64;
        let mut result = Vec::with_capacity(max_len);
        for i in 0..max_len {
            let a = self.limbs.get(i).map_or(0, |v| *v);
            let b = other.limbs.get(i).map_or(0, |v| *v);
            let s = u128::from(a) + u128::from(b) + u128::from(carry);
            result.push(lo(s));
            carry = hi(s);
        }
        if carry != 0 {
            result.push(carry);
        }
        Nat::normalize(result)
    }

    pub fn mul(&self, other: &Nat) -> Nat {
        if self.is_zero() || other.is_zero() {
            return Nat::zero();
        }
        let mut limbs: Vec<u64> = Vec::new();
        for (i, a) in self.limbs.iter().enumerate() {
            for (j, b) in other.limbs.iter().enumerate() {
                add_u128_at(&mut limbs, i + j, u128::from(*a) * u128::from(*b));
            }
        }
        Nat::normalize(limbs)
    }

    pub fn cmp(&self, other: &Nat) -> Ordering {
        if self.limbs.len() < other.limbs.len() {
            Ordering::Less
        } else if self.limbs.len() > other.limbs.len() {
            Ordering::Greater
        } else {
            for (a, b) in self.limbs.iter().rev().zip(other.limbs.iter().rev()) {
                match a.cmp(b) {
                    Ordering::Equal => continue,
                    other @ (Ordering::Less | Ordering::Greater) => return other,
                }
            }
            Ordering::Equal
        }
    }

    fn normalize(limbs: Vec<u64>) -> Nat {
        let mut limbs = limbs;
        while limbs.last() == Some(&0) && !limbs.is_empty() {
            limbs.pop();
        }
        Nat { limbs }
    }
}
