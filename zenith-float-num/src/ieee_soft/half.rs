//! Software IEEE-754 binary16 (`Ieee16`) and bfloat16 (`Bfloat16`): `u16` bit patterns and their exact
//! [`ExactNum`](crate::ExactNum) values, rounded back with any [`RoundingMode`].
//!
//! Integer bit manipulation only: no hardware float type is used in decoding or rounding.

use crate::defs::{Error, Exponent, RoundingMode};
use crate::ExactNum;
use core::num::FpCategory;

/// Layout of a 16-bit interchange format: 1 sign bit, `15 - frac` exponent bits, `frac` fraction bits.
#[derive(Clone, Copy)]
struct Layout {
    frac: u32,
    bias: Exponent,
}

/// IEEE 754 binary16: 5 exponent bits (bias 15), 10 fraction bits.
const BIN16: Layout = Layout { frac: 10, bias: 15 };
/// bfloat16: 8 exponent bits (bias 127), 7 fraction bits.
const BF16: Layout = Layout { frac: 7, bias: 127 };
const SIGN: u16 = 0x8000;

impl Layout {
    const fn exp_all(self) -> u16 {
        (1u16 << (15 - self.frac)) - 1
    }
    const fn frac_mask(self) -> u16 {
        (1u16 << self.frac) - 1
    }
    const fn hidden(self) -> u16 {
        1u16 << self.frac
    }
    const fn inf(self) -> u16 {
        self.exp_all() << self.frac
    }
    const fn qnan(self) -> u16 {
        self.inf() | (1u16 << (self.frac - 1))
    }
    const fn max_finite(self) -> u16 {
        self.inf() - 1
    }
    fn kind(self, bits: u16) -> FpCategory {
        let e = (bits >> self.frac) & self.exp_all();
        let m = bits & self.frac_mask();
        if e == self.exp_all() {
            if m == 0 {
                FpCategory::Infinite
            } else {
                FpCategory::Nan
            }
        } else if e == 0 {
            if m == 0 {
                FpCategory::Zero
            } else {
                FpCategory::Subnormal
            }
        } else {
            FpCategory::Normal
        }
    }
}

/// Red-test stub: +0 (scaled by 2^-bias) for every pattern; the green slice decodes the fields exactly.
fn to_exact(_bits: u16, f: Layout, p: usize) -> ExactNum {
    ExactNum::from_u8(0, p).ldexp(-f.bias, p, RoundingMode::None)
}

/// Red-test stub: the green slice rounds with `rm`.
fn from_exact(_x: &ExactNum, _f: Layout, _rm: RoundingMode) -> Result<u16, Error> {
    Err(Error::InvalidArgument)
}

/// Software IEEE 754 binary16 (FP16: 1 sign, 5 exponent, 10 fraction bits). Stored as `u16` bits.
#[derive(Clone, Copy, Debug)]
pub struct Ieee16(u16);

/// Software bfloat16 (BF16: 1 sign, 8 exponent, 7 fraction bits; binary32's top 16 bits). Stored as
/// `u16` bits.
#[derive(Clone, Copy, Debug)]
pub struct Bfloat16(u16);

macro_rules! impl_half {
    ($ty:ident, $fmt:ident) => {
        impl $ty {
            /// \(+0\).
            pub const ZERO: Self = Self(0);
            /// \(-0\).
            pub const NEG_ZERO: Self = Self(SIGN);
            /// Canonical quiet NaN (positive, top fraction bit set).
            pub const NAN: Self = Self($fmt.qnan());
            /// \(+\infty\).
            pub const INFINITY: Self = Self($fmt.inf());
            /// \(-\infty\).
            pub const NEG_INFINITY: Self = Self(SIGN | $fmt.inf());
            /// Largest finite value.
            pub const MAX: Self = Self($fmt.max_finite());
            /// Smallest positive normal value.
            pub const MIN_POSITIVE: Self = Self($fmt.hidden());
            /// Smallest positive subnormal value.
            pub const MIN_POSITIVE_SUBNORMAL: Self = Self(1);

            /// Build from the bit pattern.
            #[inline]
            pub const fn from_bits(bits: u16) -> Self {
                Self(bits)
            }

            /// Return the bit pattern.
            #[inline]
            pub const fn to_bits(self) -> u16 {
                self.0
            }

            /// Classify the value.
            pub fn classify(self) -> FpCategory {
                $fmt.kind(self.0)
            }

            /// True if the value is NaN.
            pub fn is_nan(self) -> bool {
                self.classify() == FpCategory::Nan
            }

            /// True if the value is \(\pm\infty\).
            pub fn is_infinite(self) -> bool {
                self.classify() == FpCategory::Infinite
            }

            /// True if the value is finite (including zero and subnormals).
            pub fn is_finite(self) -> bool {
                !self.is_nan() && !self.is_infinite()
            }

            /// True if the value is \(\pm 0\).
            pub fn is_zero(self) -> bool {
                self.classify() == FpCategory::Zero
            }

            /// True if the value is subnormal.
            pub fn is_subnormal(self) -> bool {
                self.classify() == FpCategory::Subnormal
            }

            /// True if the sign bit is set (NaN included).
            pub fn is_sign_negative(self) -> bool {
                self.0 & SIGN != 0
            }

            /// Widen to `ExactNum` at precision `p` (use at least 16 bits; precision is rounded up to
            /// the word size). Exact for every finite pattern, subnormals included; \(-0\) keeps its
            /// sign; \(\pm\infty\) maps to `INF_POS` / `INF_NEG`; every NaN maps to `NAN`.
            pub fn to_exact(self, p: usize) -> ExactNum {
                to_exact(self.0, $fmt, p)
            }

            /// Round an `ExactNum` to this format (to nearest, ties to even). Same as
            /// [`Self::from_exact_rm`] with `RoundingMode::ToEven`, which never refuses.
            pub fn from_exact(x: &ExactNum) -> Self {
                match from_exact(x, $fmt, RoundingMode::ToEven) {
                    Ok(b) => Self(b),
                    Err(_) => Self::NAN,
                }
            }

            /// Round an `ExactNum` to this format with `rm`, reading every mantissa word for the
            /// sticky bit. Underflow rounds onto the subnormal grid (or to \(\pm 0\) / the least
            /// subnormal by mode); overflow gives \(\pm\infty\) in `ToEven`, `ToOdd` and `FromZero`,
            /// the largest finite value in `ToZero`, and by sign in `Up` / `Down`. NaN gives the
            /// canonical quiet NaN; \(-0\) keeps its sign. `RoundingMode::None` refuses
            /// (`Error::InvalidArgument`) any input that is not exactly representable; no other
            /// mode refuses.
            pub fn from_exact_rm(x: &ExactNum, rm: RoundingMode) -> Result<Self, Error> {
                from_exact(x, $fmt, rm).map(Self)
            }
        }
    };
}

impl_half!(Ieee16, BIN16);
impl_half!(Bfloat16, BF16);
