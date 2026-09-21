//! Software IEEE-754 binary32 / binary64 (`Ieee32` / `Ieee64`).
//! Bits live in `u32` / `u64`. Arithmetic is integer only.

mod arith;
mod array;
mod convert;
mod simd;

use crate::defs::Error;
use crate::RoundingMode;
use crate::Sign;
use arith::{
    add_bits, cmp_bits, div_bits, fma_bits, frexp_bits, from_i32_bits, mul_bits, next_after,
    next_down, next_up, sqrt_bits, sub_bits, unpack, Class, BIN32, BIN64,
};
use core::cmp::Ordering;
use core::num::FpCategory;

pub use array::{ExactNumArray, Ieee32Array, Ieee64Array};
pub use simd::IEEE_SIMD_LANE_WIDTH;

/// Software IEEE-754 binary32 (24-bit significand, 8-bit exponent). Stored as `u32` bits.
#[derive(Clone, Copy, Debug)]
pub struct Ieee32(u32);

/// Software IEEE-754 binary64 (53-bit significand, 11-bit exponent). Stored as `u64` bits.
#[derive(Clone, Copy, Debug)]
pub struct Ieee64(u64);

macro_rules! impl_ieee {
    ($ty:ident, $bits:ty, $fmt:ident, $as_u64:ident, $from_u64:ident) => {
        impl $ty {
            /// \(+0\).
            pub const ZERO: Self = Self(0);

            /// Build from the IEEE bit pattern.
            #[inline]
            pub const fn from_bits(bits: $bits) -> Self {
                Self(bits)
            }

            /// Return the IEEE bit pattern.
            #[inline]
            pub const fn to_bits(self) -> $bits {
                self.0
            }

            fn raw(self) -> u64 {
                $as_u64(self.0)
            }

            /// Integer \(n\) converted with IEEE rounding.
            pub fn from_i32(n: i32) -> Self {
                Self($from_u64(from_i32_bits(n, $fmt)))
            }

            /// Classify the value.
            pub fn classify(self) -> FpCategory {
                match unpack(self.raw(), $fmt).class {
                    Class::Nan => FpCategory::Nan,
                    Class::Inf => FpCategory::Infinite,
                    Class::Zero => FpCategory::Zero,
                    Class::Sub => FpCategory::Subnormal,
                    Class::Norm => FpCategory::Normal,
                }
            }

            /// True if the value is NaN.
            pub fn is_nan(self) -> bool {
                unpack(self.raw(), $fmt).class == Class::Nan
            }

            /// True if the value is \(\pm\infty\).
            pub fn is_infinite(self) -> bool {
                unpack(self.raw(), $fmt).class == Class::Inf
            }

            /// True if the value is finite (including zero and subnormals).
            pub fn is_finite(self) -> bool {
                !self.is_nan() && !self.is_infinite()
            }

            /// True if the value is \(\pm 0\).
            pub fn is_zero(self) -> bool {
                unpack(self.raw(), $fmt).class == Class::Zero
            }

            /// True if the value is subnormal.
            pub fn is_subnormal(self) -> bool {
                unpack(self.raw(), $fmt).class == Class::Sub
            }

            /// True if the sign bit is set.
            pub fn is_sign_negative(self) -> bool {
                unpack(self.raw(), $fmt).sign
            }

            /// Associated error on NaN (`InvalidArgument`); `None` otherwise.
            pub fn err(self) -> Option<Error> {
                if self.is_nan() {
                    Some(Error::InvalidArgument)
                } else {
                    None
                }
            }

            /// IEEE add (to-nearest, ties to even).
            pub fn soft_add(self, rhs: Self) -> Self {
                Self($from_u64(add_bits(self.raw(), rhs.raw(), $fmt)))
            }

            /// IEEE sub.
            pub fn soft_sub(self, rhs: Self) -> Self {
                Self($from_u64(sub_bits(self.raw(), rhs.raw(), $fmt)))
            }

            /// IEEE mul.
            pub fn soft_mul(self, rhs: Self) -> Self {
                Self($from_u64(mul_bits(self.raw(), rhs.raw(), $fmt)))
            }

            /// IEEE div. \(x/0\) is \(\pm\infty\); \(0/0\) is NaN.
            pub fn soft_div(self, rhs: Self) -> Self {
                Self($from_u64(div_bits(self.raw(), rhs.raw(), $fmt)))
            }

            /// IEEE sqrt. Negative finite → NaN.
            pub fn sqrt(self) -> Self {
                Self($from_u64(sqrt_bits(self.raw(), $fmt)))
            }

            /// IEEE fused multiply-add \(a\cdot b + c\).
            pub fn mul_add(self, b: Self, c: Self) -> Self {
                Self($from_u64(fma_bits(self.raw(), b.raw(), c.raw(), $fmt)))
            }

            /// Next representable value toward \(+\infty\).
            pub fn next_up(self) -> Self {
                Self($from_u64(next_up(self.raw(), $fmt)))
            }

            /// Next representable value toward \(-\infty\).
            pub fn next_down(self) -> Self {
                Self($from_u64(next_down(self.raw(), $fmt)))
            }

            /// Next representable value toward `other`.
            pub fn next_after(self, other: Self) -> Self {
                Self($from_u64(next_after(self.raw(), other.raw(), $fmt)))
            }

            /// Split as \(m\cdot 2^e\) with \(m\in[1/2,1)\) (or a special).
            pub fn frexp(self) -> (Self, i32) {
                let (m, e) = frexp_bits(self.raw(), $fmt);
                (Self($from_u64(m)), e)
            }

            /// Negate (flip the sign bit).
            pub fn soft_neg(self) -> Self {
                Self($from_u64(self.raw() ^ $fmt.sign_mask()))
            }

            /// Absolute value.
            pub fn abs(self) -> Self {
                Self($from_u64(self.raw() & !$fmt.sign_mask()))
            }

            /// Copy the sign of `sign` onto `self`.
            pub fn copysign(self, sign: Self) -> Self {
                let mag = self.raw() & !$fmt.sign_mask();
                let s = sign.raw() & $fmt.sign_mask();
                Self($from_u64(mag | s))
            }

            /// Sign, or `None` if NaN.
            pub fn sign(self) -> Option<Sign> {
                if self.is_nan() {
                    None
                } else if unpack(self.raw(), $fmt).sign {
                    Some(Sign::Neg)
                } else {
                    Some(Sign::Pos)
                }
            }

            /// Compare; `None` if either is NaN. \(+0 = -0\).
            pub fn soft_cmp(self, other: Self) -> Option<Ordering> {
                cmp_bits(self.raw(), other.raw(), $fmt).map(|d| {
                    if d < 0 {
                        Ordering::Less
                    } else if d > 0 {
                        Ordering::Greater
                    } else {
                        Ordering::Equal
                    }
                })
            }
        }

        impl PartialEq for $ty {
            fn eq(&self, other: &Self) -> bool {
                self.soft_cmp(*other) == Some(Ordering::Equal)
            }
        }

        impl core::ops::Add for $ty {
            type Output = Self;
            fn add(self, rhs: Self) -> Self {
                $ty::soft_add(self, rhs)
            }
        }
        impl core::ops::Sub for $ty {
            type Output = Self;
            fn sub(self, rhs: Self) -> Self {
                $ty::soft_sub(self, rhs)
            }
        }
        impl core::ops::Mul for $ty {
            type Output = Self;
            fn mul(self, rhs: Self) -> Self {
                $ty::soft_mul(self, rhs)
            }
        }
        impl core::ops::Div for $ty {
            type Output = Self;
            fn div(self, rhs: Self) -> Self {
                $ty::soft_div(self, rhs)
            }
        }
        impl core::ops::Neg for $ty {
            type Output = Self;
            fn neg(self) -> Self {
                $ty::soft_neg(self)
            }
        }
    };
}

const fn u32_as_u64(x: u32) -> u64 {
    x as u64
}
const fn u64_as_u64(x: u64) -> u64 {
    x
}
const fn u64_as_u32(x: u64) -> u32 {
    x as u32
}
const fn u64_id(x: u64) -> u64 {
    x
}

impl_ieee!(Ieee32, u32, BIN32, u32_as_u64, u64_as_u32);
impl_ieee!(Ieee64, u64, BIN64, u64_as_u64, u64_id);

impl Ieee32 {
    /// Canonical quiet NaN.
    pub const NAN: Self = Self(0x7FC0_0000);
    /// \(+\infty\).
    pub const INFINITY: Self = Self(0x7F80_0000);
    /// \(-\infty\).
    pub const NEG_INFINITY: Self = Self(0xFF80_0000);
    /// \(-0\).
    pub const NEG_ZERO: Self = Self(0x8000_0000);
}

impl Ieee64 {
    /// Canonical quiet NaN.
    pub const NAN: Self = Self(0x7FF8_0000_0000_0000);
    /// \(+\infty\).
    pub const INFINITY: Self = Self(0x7FF0_0000_0000_0000);
    /// \(-\infty\).
    pub const NEG_INFINITY: Self = Self(0xFFF0_0000_0000_0000);
    /// \(-0\).
    pub const NEG_ZERO: Self = Self(0x8000_0000_0000_0000);
}

impl Ieee32 {
    /// Widen to `ExactNum` at precision `p` (use at least 32 bits).
    pub fn to_exact(self, p: usize) -> crate::ExactNum {
        convert::to_exact(self.to_bits() as u64, BIN32, p)
    }

    /// Round an `ExactNum` to binary32 (to-nearest, ties to even).
    pub fn from_exact(x: &crate::ExactNum) -> Self {
        Self::from_bits(convert::from_exact(x, BIN32) as u32)
    }

    /// Widen then apply `op`, then IEEE-round back.
    pub fn map_exact<F>(self, p: usize, _rm: RoundingMode, op: F) -> Self
    where
        F: FnOnce(&crate::ExactNum) -> crate::ExactNum,
    {
        Self::from_exact(&op(&self.to_exact(p)))
    }
}

impl Ieee64 {
    /// Widen to `ExactNum` at precision `p` (use at least 64 bits).
    pub fn to_exact(self, p: usize) -> crate::ExactNum {
        convert::to_exact(self.to_bits(), BIN64, p)
    }

    /// Round an `ExactNum` to binary64 (to-nearest, ties to even).
    pub fn from_exact(x: &crate::ExactNum) -> Self {
        Self::from_bits(convert::from_exact(x, BIN64))
    }

    /// Widen then apply `op`, then IEEE-round back.
    pub fn map_exact<F>(self, p: usize, _rm: RoundingMode, op: F) -> Self
    where
        F: FnOnce(&crate::ExactNum) -> crate::ExactNum,
    {
        Self::from_exact(&op(&self.to_exact(p)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const ONE32: u32 = 0x3F80_0000;
    const TWO32: u32 = 0x4000_0000;
    const HALF32: u32 = 0x3F00_0000;
    const THIRD32: u32 = 0x3EAA_AAAB;
    const THREE32: u32 = 0x4040_0000;
    const INF32: u32 = 0x7F80_0000;
    const QNAN32: u32 = 0x7FC0_0000;

    const ONE64: u64 = 0x3FF0_0000_0000_0000;
    const TWO64: u64 = 0x4000_0000_0000_0000;
    const HALF64: u64 = 0x3FE0_0000_0000_0000;
    const THIRD64: u64 = 0x3FD5_5555_5555_5555;
    const THREE64: u64 = 0x4008_0000_0000_0000;
    const INF64: u64 = 0x7FF0_0000_0000_0000;
    const QNAN64: u64 = 0x7FF8_0000_0000_0000;
    const MONE64: u64 = 0xBFF0_0000_0000_0000;

    #[test]
    fn bin32_add_mul_div_bits() {
        let one = Ieee32::from_bits(ONE32);
        let two = Ieee32::from_bits(TWO32);
        let half = Ieee32::from_bits(HALF32);
        let three = Ieee32::from_bits(THREE32);
        assert_eq!((one + one).to_bits(), TWO32);
        assert_eq!((one + two).to_bits(), THREE32);
        assert_eq!((two * half).to_bits(), ONE32);
        assert_eq!((one / Ieee32::from_i32(3)).to_bits(), THIRD32);
        assert_eq!((three - one).to_bits(), TWO32);
        assert!(Ieee32::from_bits(0)
            .soft_add(Ieee32::from_bits(0x8000_0000))
            .is_zero());
        assert_eq!(Ieee32::from_i32(1).to_bits(), ONE32);
        assert_eq!(Ieee32::from_i32(-1).to_bits(), 0xBF80_0000);
        assert_eq!(one.sqrt().to_bits(), ONE32);
        assert!((Ieee32::from_i32(-1).sqrt()).is_nan());
        assert_eq!((one / Ieee32::ZERO).to_bits(), INF32);
        assert!((Ieee32::ZERO / Ieee32::ZERO).is_nan());
        assert_eq!(Ieee32::INFINITY.to_bits(), INF32);
        assert_eq!(Ieee32::NAN.to_bits(), QNAN32);
        let tiny = Ieee32::from_bits(1);
        assert!(tiny.is_subnormal());
        assert_eq!((tiny + tiny).to_bits(), 2);
    }

    #[test]
    fn bin64_add_mul_div_bits() {
        let one = Ieee64::from_bits(ONE64);
        let two = Ieee64::from_bits(TWO64);
        let half = Ieee64::from_bits(HALF64);
        let three = Ieee64::from_bits(THREE64);
        assert_eq!((one + one).to_bits(), TWO64);
        assert_eq!((one + two).to_bits(), THREE64);
        assert_eq!((two * half).to_bits(), ONE64);
        assert_eq!((one / Ieee64::from_i32(3)).to_bits(), THIRD64);
        assert_eq!((three - one).to_bits(), TWO64);
        assert_eq!(Ieee64::from_i32(1).to_bits(), ONE64);
        assert_eq!(Ieee64::from_i32(-1).to_bits(), MONE64);
        assert_eq!(one.sqrt().to_bits(), ONE64);
        assert_eq!(Ieee64::from_i32(4).sqrt().to_bits(), TWO64);
        // IEEE sqrt(2) (to-nearest, ties to even).
        assert_eq!(two.sqrt().to_bits(), 0x3FF6_A09E_667F_3BCD);
        assert_eq!((one / Ieee64::ZERO).to_bits(), INF64);
        assert_eq!(Ieee64::INFINITY.to_bits(), INF64);
        assert_eq!(Ieee64::NAN.to_bits(), QNAN64);
        let tiny = Ieee64::from_bits(1);
        assert!(tiny.is_subnormal());
        assert_eq!((tiny + tiny).to_bits(), 2);
        assert_eq!(one.mul_add(two, one).to_bits(), THREE64);
        assert_eq!(one.soft_neg().to_bits(), MONE64);
        assert_eq!(Ieee64::ZERO.soft_cmp(Ieee64::NEG_ZERO), Some(Ordering::Equal));
        assert_eq!(one.next_up().next_down().to_bits(), ONE64);
        let (m, e) = two.frexp();
        assert_eq!(e, 2);
        assert_eq!(m.to_bits(), HALF64);
    }

    #[test]
    fn widen_roundtrip_one() {
        let x = Ieee64::from_bits(ONE64);
        let e = x.to_exact(64);
        assert_eq!(Ieee64::from_exact(&e).to_bits(), ONE64);
        let y = Ieee32::from_bits(ONE32);
        assert_eq!(Ieee32::from_exact(&y.to_exact(64)).to_bits(), ONE32);
    }
}
