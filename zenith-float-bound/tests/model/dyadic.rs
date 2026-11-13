//! `MagV` (the model of a `Mag`) and `Dy`, a nonnegative dyadic rational `n * 2^e` with an `i128` exponent. Comparisons that the model cannot decide return `None`, which tests turn into a failure.

use super::nat::Nat;
use core::cmp::Ordering;

/// Model of the radius type: Finite{exp, mant} denotes mant * 2^(exp - 64); mant == 0 (with exp == 0) is the canonical zero; Infinity is unbounded.
#[derive(Debug, PartialEq, Eq)]
pub enum MagV {
    Finite { exp: i64, mant: u64 },
    Infinity,
}

/// A nonnegative dyadic rational n * 2^e.
pub struct Dy {
    n: Nat,
    e: i128,
}

pub const MAX_ALIGN: i128 = 4096;

impl Dy {
    pub fn from_mag(exp: i64, mant: u64) -> Dy {
        Dy {
            n: Nat::from_u64(mant),
            e: i128::from(exp) - 64,
        }
    }

    pub fn from_magv(m: &MagV) -> Option<Dy> {
        match m {
            MagV::Finite { exp, mant } => Some(Dy::from_mag(*exp, *mant)),
            MagV::Infinity => None,
        }
    }

    pub fn pow2(e: i128) -> Dy {
        Dy {
            n: Nat::from_u64(1),
            e,
        }
    }

    pub fn is_zero(&self) -> bool {
        self.n.is_zero()
    }

    /// True when both values are nonzero and their exponents differ by more than MAX_ALIGN.
    pub fn far_apart(&self, other: &Dy) -> bool {
        !self.is_zero() && !other.is_zero() && (self.e - other.e).abs() > MAX_ALIGN
    }

    pub fn mul(&self, other: &Dy) -> Dy {
        Dy {
            n: self.n.mul(&other.n),
            e: self.e + other.e,
        }
    }

    pub fn add(&self, other: &Dy) -> Option<Dy> {
        if self.is_zero() {
            return Some(Dy {
                n: other.n.add(&Nat::zero()),
                e: other.e,
            });
        }
        if other.is_zero() {
            return Some(Dy {
                n: self.n.add(&Nat::zero()),
                e: self.e,
            });
        }

        if self.far_apart(other) {
            return None;
        }

        let lowest = self.e.min(other.e);
        let d1 = self.e - lowest;
        let d2 = other.e - lowest;

        let shift1 = usize::try_from(d1).ok()?;
        let shift2 = usize::try_from(d2).ok()?;

        let n1_shifted = self.n.shl(shift1);
        let n2_shifted = other.n.shl(shift2);
        let sum = n1_shifted.add(&n2_shifted);

        Some(Dy { n: sum, e: lowest })
    }

    pub fn cmp(&self, other: &Dy) -> Option<Ordering> {
        if self.is_zero() && other.is_zero() {
            return Some(Ordering::Equal);
        }
        if self.is_zero() {
            return Some(Ordering::Less);
        }
        if other.is_zero() {
            return Some(Ordering::Greater);
        }

        let top1 = i128::from(self.n.bit_len()) + self.e;
        let top2 = i128::from(other.n.bit_len()) + other.e;

        match top1.cmp(&top2) {
            Ordering::Equal => {
                let d = self.e - other.e;
                if d >= 0 {
                    let shift = usize::try_from(d).ok()?;
                    Some(self.n.shl(shift).cmp(&other.n))
                } else {
                    let shift = usize::try_from(-d).ok()?;
                    Some(self.n.cmp(&other.n.shl(shift)))
                }
            }
            ord @ (Ordering::Less | Ordering::Greater) => Some(ord),
        }
    }
}
