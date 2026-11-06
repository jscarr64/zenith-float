//! Software IEEE-754 binary16 (`Ieee16`) and bfloat16 (`Bfloat16`): `u16` bit patterns and their exact
//! [`ExactNum`](crate::ExactNum) values, rounded back with any [`RoundingMode`].
//!
//! Integer bit manipulation only: no hardware float type is used in decoding or rounding.

use crate::defs::{Error, Exponent, RoundingMode, Word, WORD_BIT_SIZE};
use crate::ExactNum;
use core::num::FpCategory;
use core::ops::Neg;

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

/// Exact value of `bits` at precision `p`.
fn to_exact(bits: u16, f: Layout, p: usize) -> ExactNum {
    let neg = bits & SIGN != 0;
    let e = (bits >> f.frac) & f.exp_all();
    let m = bits & f.frac_mask();
    let fr = Exponent::from(f.frac);
    let mag = match f.kind(bits) {
        FpCategory::Nan => return crate::NAN,
        FpCategory::Infinite => return if neg { crate::INF_NEG } else { crate::INF_POS },
        FpCategory::Zero => ExactNum::from_u8(0, p),
        // (-1)^s * m * 2^(1 - bias - frac)
        FpCategory::Subnormal => {
            ExactNum::from_u64(u64::from(m), p).ldexp(1 - f.bias - fr, p, RoundingMode::None)
        }
        // (-1)^s * (2^frac + m) * 2^(e - bias - frac)
        FpCategory::Normal => ExactNum::from_u64(u64::from(m | f.hidden()), p).ldexp(
            Exponent::from(e) - f.bias - fr,
            p,
            RoundingMode::None,
        ),
    };
    if neg {
        mag.neg()
    } else {
        mag
    }
}

/// Bit `t` of the mantissa counted from its most significant bit (`t = 0`). The words are least
/// significant first.
fn bit_from_top(words: &[Word], t: usize) -> bool {
    let wi = t / WORD_BIT_SIZE;
    if wi >= words.len() {
        return false;
    }
    match words.get(words.len() - 1 - wi) {
        Some(w) => (*w >> (WORD_BIT_SIZE - 1 - t % WORD_BIT_SIZE)) & 1 == 1,
        None => false,
    }
}

/// True if any mantissa bit at position `t` or further down (counted from the top) is set. Every
/// lower word is inspected, so a sticky bit in the last limb of a wide operand is seen.
fn any_from_top(words: &[Word], t: usize) -> bool {
    let wi = t / WORD_BIT_SIZE;
    if wi >= words.len() {
        return false;
    }
    let idx = words.len() - 1 - wi;
    let low = WORD_BIT_SIZE - t % WORD_BIT_SIZE;
    let one: Word = 1;
    let mask = if low == WORD_BIT_SIZE {
        Word::MAX
    } else {
        (one << low) - 1
    };
    let partial = words.get(idx).is_some_and(|w| *w & mask != 0);
    partial || words.get(..idx).is_some_and(|lower| lower.iter().any(|w| *w != 0))
}

/// Number of zero bits above the leading set bit, or `None` if every word is zero.
fn leading_zero_bits(words: &[Word]) -> Option<usize> {
    let mut acc = 0usize;
    for w in words.iter().rev() {
        if *w != 0 {
            return usize::try_from(w.leading_zeros())
                .ok()
                .map(|lz| acc + lz);
        }
        acc += WORD_BIT_SIZE;
    }
    None
}

/// The result of a magnitude too large for the format, by rounding direction.
fn overflow(sign: u16, f: Layout, rm: RoundingMode) -> Result<u16, Error> {
    let neg = sign != 0;
    let to_inf = match rm {
        RoundingMode::None => return Err(Error::InvalidArgument),
        RoundingMode::ToEven | RoundingMode::ToOdd | RoundingMode::FromZero => true,
        RoundingMode::ToZero => false,
        RoundingMode::Up => !neg,
        RoundingMode::Down => neg,
    };
    Ok(sign | if to_inf { f.inf() } else { f.max_finite() })
}

/// Round `x` to the format with `rm`. `RoundingMode::None` refuses an inexact result.
fn from_exact(x: &ExactNum, f: Layout, rm: RoundingMode) -> Result<u16, Error> {
    if x.is_nan() {
        return Ok(f.qnan());
    }
    let sign = if x.is_negative() { SIGN } else { 0 };
    if x.is_inf() {
        return Ok(sign | f.inf());
    }
    if x.is_zero() {
        return Ok(sign);
    }
    let (Some(words), Some(e)) = (x.mantissa_digits(), x.exponent()) else {
        return Ok(sign);
    };
    let Some(lead) = leading_zero_bits(words) else {
        return Ok(sign);
    };
    let neg = sign != 0;
    let lead_i = i128::try_from(lead).map_err(|_| Error::InvalidArgument)?;
    let fr = i128::from(f.frac);
    let emin = 1 - i128::from(f.bias);
    let emax = i128::from(f.bias);
    // |x| lies in [2^top, 2^(top + 1)): the mantissa is a fraction in [1/2, 1) times 2^e.
    let top = i128::from(e) - 1 - lead_i;
    if top > emax {
        return overflow(sign, f, rm);
    }
    // Exponent of the last kept bit, and the number of kept bits (frac + 1 for a normal result,
    // fewer, possibly none, below the normal range).
    let mut lsb = top.max(emin) - fr;
    let kept = top - lsb + 1;
    let (mut q, round, sticky) = if kept < 0 {
        (0u16, false, true)
    } else {
        let k = usize::try_from(kept).map_err(|_| Error::InvalidArgument)?;
        let mut q = 0u16;
        for i in 0..k {
            q = (q << 1) | u16::from(bit_from_top(words, lead + i));
        }
        (
            q,
            bit_from_top(words, lead + k),
            any_from_top(words, lead + k + 1),
        )
    };
    let inexact = round || sticky;
    let odd = q & 1 == 1;
    let up = match rm {
        RoundingMode::None => {
            if inexact {
                return Err(Error::InvalidArgument);
            }
            false
        }
        RoundingMode::ToEven => round && (sticky || odd),
        RoundingMode::ToOdd => round && (sticky || !odd),
        RoundingMode::Up => inexact && !neg,
        RoundingMode::Down => inexact && neg,
        RoundingMode::ToZero => false,
        RoundingMode::FromZero => inexact,
    };
    if up {
        q += 1;
        if q == f.hidden() << 1 {
            q >>= 1;
            lsb += 1;
        }
    }
    if q < f.hidden() {
        // Subnormal or zero: lsb is emin - frac here, so q is the fraction field itself.
        return Ok(sign | q);
    }
    let biased = lsb + fr + i128::from(f.bias);
    if biased >= i128::from(f.exp_all()) {
        return overflow(sign, f, rm);
    }
    let ef = u16::try_from(biased).map_err(|_| Error::InvalidArgument)?;
    Ok(sign | (ef << f.frac) | (q - f.hidden()))
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
