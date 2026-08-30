//! Dense 1-D arrays of software IEEE values and `ExactNum`.

use super::{Ieee32, Ieee64};
use crate::defs::RoundingMode;
use crate::Consts;
use crate::ExactNum;
use alloc::vec::Vec;

/// Contiguous binary32 lanes (`u32` bits).
#[derive(Clone, Debug)]
pub struct Ieee32Array {
    bits: Vec<u32>,
}

/// Contiguous binary64 lanes (`u64` bits).
#[derive(Clone, Debug)]
pub struct Ieee64Array {
    bits: Vec<u64>,
}

/// Batch of `ExactNum` values at a shared precision `p`.
#[derive(Clone, Debug)]
pub struct ExactNumArray {
    p: usize,
    vals: Vec<ExactNum>,
}

macro_rules! impl_ieee_array {
    ($arr:ident, $scalar:ident, $bits:ty) => {
        impl $arr {
            /// Empty array.
            pub fn new() -> Self {
                Self { bits: Vec::new() }
            }

            /// `n` copies of `fill`.
            pub fn filled(n: usize, fill: $scalar) -> Self {
                Self {
                    bits: alloc::vec![fill.to_bits(); n],
                }
            }

            /// From a slice of IEEE bit patterns.
            pub fn from_bits(bits: &[$bits]) -> Self {
                Self {
                    bits: bits.to_vec(),
                }
            }

            /// From scalars.
            pub fn from_values(vals: &[$scalar]) -> Self {
                Self {
                    bits: vals.iter().map(|v| v.to_bits()).collect(),
                }
            }

            /// Number of lanes.
            pub fn len(&self) -> usize {
                self.bits.len()
            }

            /// True if there are no lanes.
            pub fn is_empty(&self) -> bool {
                self.bits.is_empty()
            }

            /// IEEE bit patterns.
            pub fn as_bits(&self) -> &[$bits] {
                &self.bits
            }

            /// Lane `i`, or `None` if out of range.
            pub fn get(&self, i: usize) -> Option<$scalar> {
                self.bits.get(i).copied().map($scalar::from_bits)
            }

            /// Elementwise add. Lengths must match.
            pub fn add(&self, rhs: &Self) -> Option<Self> {
                self.zip_op(rhs, $scalar::add)
            }

            /// Add a scalar to every lane.
            pub fn add_scalar(&self, s: $scalar) -> Self {
                self.map(|x| x.add(s))
            }

            /// Elementwise sub.
            pub fn sub(&self, rhs: &Self) -> Option<Self> {
                self.zip_op(rhs, $scalar::sub)
            }

            /// Elementwise mul.
            pub fn mul(&self, rhs: &Self) -> Option<Self> {
                self.zip_op(rhs, $scalar::mul)
            }

            /// Multiply every lane by a scalar.
            pub fn mul_scalar(&self, s: $scalar) -> Self {
                self.map(|x| x.mul(s))
            }

            /// Elementwise div.
            pub fn div(&self, rhs: &Self) -> Option<Self> {
                self.zip_op(rhs, $scalar::div)
            }

            /// Elementwise sqrt.
            pub fn sqrt(&self) -> Self {
                self.map($scalar::sqrt)
            }

            /// Sequential IEEE sum (one rounding per add).
            pub fn sum(&self) -> $scalar {
                let mut acc = $scalar::ZERO;
                for &b in &self.bits {
                    acc = acc.add($scalar::from_bits(b));
                }
                acc
            }

            /// Sequential dot product (mul then add, two rounds per term).
            pub fn dot(&self, rhs: &Self) -> Option<$scalar> {
                if self.len() != rhs.len() {
                    return None;
                }
                let mut acc = $scalar::ZERO;
                for (a, b) in self.bits.iter().zip(rhs.bits.iter()) {
                    acc = acc.add($scalar::from_bits(*a).mul($scalar::from_bits(*b)));
                }
                Some(acc)
            }

            fn zip_op(&self, rhs: &Self, op: fn($scalar, $scalar) -> $scalar) -> Option<Self> {
                if self.len() != rhs.len() {
                    return None;
                }
                Some(Self {
                    bits: self
                        .bits
                        .iter()
                        .zip(rhs.bits.iter())
                        .map(|(a, b)| op($scalar::from_bits(*a), $scalar::from_bits(*b)).to_bits())
                        .collect(),
                })
            }

            fn map(&self, op: impl Fn($scalar) -> $scalar) -> Self {
                Self {
                    bits: self
                        .bits
                        .iter()
                        .map(|&b| op($scalar::from_bits(b)).to_bits())
                        .collect(),
                }
            }

            /// Widen each lane, apply `op`, IEEE-round back.
            pub fn map_exact<F>(&self, p: usize, mut op: F) -> Self
            where
                F: FnMut(&ExactNum) -> ExactNum,
            {
                Self {
                    bits: self
                        .bits
                        .iter()
                        .map(|&b| {
                            let x = $scalar::from_bits(b).to_exact(p);
                            $scalar::from_exact(&op(&x)).to_bits()
                        })
                        .collect(),
                }
            }
        }

        impl Default for $arr {
            fn default() -> Self {
                Self::new()
            }
        }
    };
}

impl_ieee_array!(Ieee32Array, Ieee32, u32);
impl_ieee_array!(Ieee64Array, Ieee64, u64);

macro_rules! ieee_array_specials {
    ($arr:ident, $_scalar:ident, $p:expr) => {
        impl $arr {
            /// Elementwise `exp` via `ExactNum`.
            pub fn exp(&self, cc: &mut Consts) -> Self {
                self.map_exact($p, |x| x.exp($p, RoundingMode::ToEven, cc))
            }
            /// Elementwise `ln` via `ExactNum`.
            pub fn ln(&self, cc: &mut Consts) -> Self {
                self.map_exact($p, |x| x.ln($p, RoundingMode::ToEven, cc))
            }
            /// Elementwise `sin` via `ExactNum`.
            pub fn sin(&self, cc: &mut Consts) -> Self {
                self.map_exact($p, |x| x.sin($p, RoundingMode::ToEven, cc))
            }
            /// Elementwise `cos` via `ExactNum`.
            pub fn cos(&self, cc: &mut Consts) -> Self {
                self.map_exact($p, |x| x.cos($p, RoundingMode::ToEven, cc))
            }
            /// Elementwise `tan` via `ExactNum`.
            pub fn tan(&self, cc: &mut Consts) -> Self {
                self.map_exact($p, |x| x.tan($p, RoundingMode::ToEven, cc))
            }
            /// Elementwise `erf` via `ExactNum`.
            pub fn erf(&self, cc: &mut Consts) -> Self {
                self.map_exact($p, |x| x.erf($p, RoundingMode::ToEven, cc))
            }
            /// Elementwise `erfc` via `ExactNum`.
            pub fn erfc(&self, cc: &mut Consts) -> Self {
                self.map_exact($p, |x| x.erfc($p, RoundingMode::ToEven, cc))
            }
            /// Elementwise `gamma` via `ExactNum`.
            pub fn gamma(&self, cc: &mut Consts) -> Self {
                self.map_exact($p, |x| x.gamma($p, RoundingMode::ToEven, cc))
            }
            /// Elementwise `ei` via `ExactNum`.
            pub fn ei(&self, cc: &mut Consts) -> Self {
                self.map_exact($p, |x| x.ei($p, RoundingMode::ToEven, cc))
            }
            /// Elementwise `si` via `ExactNum`.
            pub fn si(&self, cc: &mut Consts) -> Self {
                self.map_exact($p, |x| x.si($p, RoundingMode::ToEven, cc))
            }
        }
    };
}

ieee_array_specials!(Ieee32Array, Ieee32, 64);
ieee_array_specials!(Ieee64Array, Ieee64, 128);

impl ExactNumArray {
    /// Empty array at precision `p`.
    pub fn new(p: usize) -> Self {
        Self {
            p,
            vals: Vec::new(),
        }
    }

    /// `n` copies of `fill` (rounded to `p`).
    pub fn filled(p: usize, n: usize, fill: &ExactNum) -> Self {
        let mut v = fill.clone();
        let _ = v.set_precision(p, RoundingMode::ToEven);
        Self {
            p,
            vals: alloc::vec![v; n],
        }
    }

    /// From values; each is rounded to `p`.
    pub fn from_values(p: usize, vals: &[ExactNum]) -> Self {
        Self {
            p,
            vals: vals
                .iter()
                .map(|x| {
                    let mut y = x.clone();
                    let _ = y.set_precision(p, RoundingMode::ToEven);
                    y
                })
                .collect(),
        }
    }

    /// Shared precision.
    pub fn precision(&self) -> usize {
        self.p
    }

    /// Number of lanes.
    pub fn len(&self) -> usize {
        self.vals.len()
    }

    /// True if there are no lanes.
    pub fn is_empty(&self) -> bool {
        self.vals.is_empty()
    }

    /// Lane `i`.
    pub fn get(&self, i: usize) -> Option<&ExactNum> {
        self.vals.get(i)
    }

    /// All values.
    pub fn as_slice(&self) -> &[ExactNum] {
        &self.vals
    }

    /// Elementwise add at `p`.
    pub fn add(&self, rhs: &Self) -> Option<Self> {
        self.zip(rhs, |a, b| a.add(b, self.p, RoundingMode::ToEven))
    }

    /// Add a scalar to every lane.
    pub fn add_scalar(&self, s: &ExactNum) -> Self {
        Self {
            p: self.p,
            vals: self
                .vals
                .iter()
                .map(|x| x.add(s, self.p, RoundingMode::ToEven))
                .collect(),
        }
    }

    /// Elementwise sub.
    pub fn sub(&self, rhs: &Self) -> Option<Self> {
        self.zip(rhs, |a, b| a.sub(b, self.p, RoundingMode::ToEven))
    }

    /// Elementwise mul.
    pub fn mul(&self, rhs: &Self) -> Option<Self> {
        self.zip(rhs, |a, b| a.mul(b, self.p, RoundingMode::ToEven))
    }

    /// Multiply every lane by a scalar.
    pub fn mul_scalar(&self, s: &ExactNum) -> Self {
        Self {
            p: self.p,
            vals: self
                .vals
                .iter()
                .map(|x| x.mul(s, self.p, RoundingMode::ToEven))
                .collect(),
        }
    }

    /// Elementwise div.
    pub fn div(&self, rhs: &Self) -> Option<Self> {
        self.zip(rhs, |a, b| a.div(b, self.p, RoundingMode::ToEven))
    }

    /// Elementwise sqrt.
    pub fn sqrt(&self) -> Self {
        Self {
            p: self.p,
            vals: self
                .vals
                .iter()
                .map(|x| x.sqrt(self.p, RoundingMode::ToEven))
                .collect(),
        }
    }

    /// Sequential sum at `p`.
    pub fn sum(&self) -> ExactNum {
        let mut acc = ExactNum::from_u8(0, self.p);
        for v in &self.vals {
            acc = acc.add(v, self.p, RoundingMode::ToEven);
        }
        acc
    }

    /// Sequential dot product at `p`.
    pub fn dot(&self, rhs: &Self) -> Option<ExactNum> {
        if self.len() != rhs.len() {
            return None;
        }
        let mut acc = ExactNum::from_u8(0, self.p);
        for (a, b) in self.vals.iter().zip(rhs.vals.iter()) {
            let t = a.mul(b, self.p, RoundingMode::ToEven);
            acc = acc.add(&t, self.p, RoundingMode::ToEven);
        }
        Some(acc)
    }

    /// Elementwise `sin`.
    pub fn sin(&self, cc: &mut Consts) -> Self {
        self.map_unary(|x| x.sin(self.p, RoundingMode::ToEven, cc))
    }

    /// Elementwise `exp`.
    pub fn exp(&self, cc: &mut Consts) -> Self {
        self.map_unary(|x| x.exp(self.p, RoundingMode::ToEven, cc))
    }

    /// Elementwise `ln`.
    pub fn ln(&self, cc: &mut Consts) -> Self {
        self.map_unary(|x| x.ln(self.p, RoundingMode::ToEven, cc))
    }

    /// Elementwise `erf`.
    pub fn erf(&self, cc: &mut Consts) -> Self {
        self.map_unary(|x| x.erf(self.p, RoundingMode::ToEven, cc))
    }

    fn zip(&self, rhs: &Self, op: impl Fn(&ExactNum, &ExactNum) -> ExactNum) -> Option<Self> {
        if self.len() != rhs.len() {
            return None;
        }
        Some(Self {
            p: self.p,
            vals: self
                .vals
                .iter()
                .zip(rhs.vals.iter())
                .map(|(a, b)| op(a, b))
                .collect(),
        })
    }

    fn map_unary(&self, mut op: impl FnMut(&ExactNum) -> ExactNum) -> Self {
        Self {
            p: self.p,
            vals: self.vals.iter().map(|x| op(x)).collect(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bin64_array_add_dot() {
        let one = Ieee64::from_i32(1);
        let two = Ieee64::from_i32(2);
        let a = Ieee64Array::from_values(&[one, two]);
        let b = Ieee64Array::from_values(&[two, one]);
        let s = a.add(&b).unwrap();
        assert_eq!(s.get(0).unwrap().to_bits(), Ieee64::from_i32(3).to_bits());
        assert_eq!(s.get(1).unwrap().to_bits(), Ieee64::from_i32(3).to_bits());
        let d = a.dot(&b).unwrap();
        assert_eq!(d.to_bits(), Ieee64::from_i32(4).to_bits());
        assert_eq!(a.sum().to_bits(), Ieee64::from_i32(3).to_bits());
        let scaled = a.mul_scalar(two);
        assert_eq!(scaled.get(0).unwrap().to_bits(), two.to_bits());
    }

    #[test]
    fn bin32_array_sqrt() {
        let four = Ieee32::from_i32(4);
        let a = Ieee32Array::filled(3, four);
        let r = a.sqrt();
        assert_eq!(r.get(0).unwrap().to_bits(), Ieee32::from_i32(2).to_bits());
        assert_eq!(r.len(), 3);
    }

    #[test]
    fn bin64_array_exp_sin() {
        let mut cc = Consts::new().unwrap();
        let z = Ieee64Array::from_values(&[Ieee64::ZERO]);
        let e = z.exp(&mut cc);
        assert_eq!(e.get(0).unwrap().to_bits(), Ieee64::from_i32(1).to_bits());
        let s = z.sin(&mut cc);
        assert!(s.get(0).unwrap().is_zero());
    }

    #[test]
    fn exact_array_add() {
        let p = 64;
        let one = ExactNum::from_u8(1, p);
        let two = ExactNum::from_u8(2, p);
        let a = ExactNumArray::from_values(p, &[one.clone(), two.clone()]);
        let b = ExactNumArray::from_values(p, &[two.clone(), one.clone()]);
        let s = a.add(&b).unwrap();
        assert_eq!(s.get(0).unwrap().cmp(&ExactNum::from_u8(3, p)), Some(0));
        let d = a.dot(&b).unwrap();
        assert_eq!(d.cmp(&ExactNum::from_u8(4, p)), Some(0));
    }
}
