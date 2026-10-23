// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! Definitions.

use core::alloc::LayoutError;
use core::fmt::Display;

#[cfg(feature = "std")]
use std::collections::TryReserveError;

#[cfg(not(feature = "std"))]
use alloc::collections::TryReserveError;

/// A word.
#[cfg(not(target_pointer_width = "32"))]
pub type Word = u64;

/// Doubled word.
#[cfg(not(target_pointer_width = "32"))]
pub type DoubleWord = u128;

/// Word with sign.
#[cfg(not(target_pointer_width = "32"))]
pub type SignedWord = i128;

/// A word.
#[cfg(target_pointer_width = "32")]
pub type Word = u32;

/// Doubled word.
#[cfg(target_pointer_width = "32")]
pub type DoubleWord = u64;

/// Word with sign.
#[cfg(target_pointer_width = "32")]
pub type SignedWord = i64;

// Word-sized values are used as indices throughout the implementation.
const _: [(); 1] = [(); (core::mem::size_of::<Word>() <= core::mem::size_of::<usize>()) as usize];

/// An exponent.
pub type Exponent = i32;

/// Maximum exponent value. The same on every target.
pub const EXPONENT_MAX: Exponent = Exponent::MAX;

/// Minimum exponent value. The same on every target.
pub const EXPONENT_MIN: Exponent = Exponent::MIN;

/// Exponent of the smallest positive normal number (leading mantissa bit set).
///
/// Subnormal numbers share this exponent and have the leading bit clear, so
/// `ExactNum::min_positive_normal(p).exponent()` and `ExactNum::min_positive(p).exponent()`
/// are both `Some(MIN_NORMAL_EXPONENT)`. Equal to [`EXPONENT_MIN`].
pub const MIN_NORMAL_EXPONENT: Exponent = EXPONENT_MIN;

/// Smallest accepted precision in bits. Precision `0` is rejected with
/// `Error::InvalidArgument` (a NaN from the `ExactNum` API). Precision is rounded up to whole words.
pub const MIN_PRECISION: usize = 1;

/// Largest accepted precision in bits; anything above it is rejected with `Error::InvalidArgument`.
///
/// Exponent arithmetic is done in `i64` and sizes in `usize` (32 bits on some targets), so the
/// limit is the smaller of `i64::MAX / 2 + EXPONENT_MIN` and `usize::MAX / 4`, minus one.
/// Allocating a mantissa this large still fails with `Error::MemoryAllocation` on real hardware.
pub const MAX_PRECISION: usize = {
    let by_exponent = (i64::MAX / 2 + EXPONENT_MIN as i64) as u64;
    let by_size = (usize::MAX / 4) as u64;
    let lim = if by_exponent < by_size { by_exponent } else { by_size };
    (lim - 1) as usize
};

/// Maximum value of a word.
pub const WORD_MAX: Word = Word::MAX;

/// Base of words.
pub const WORD_BASE: DoubleWord = WORD_MAX as DoubleWord + 1;

/// Size of a word in bits.
pub const WORD_BIT_SIZE: usize = core::mem::size_of::<Word>() * 8;

/// Cases per `proptest` property under `cargo test` (no `mpfr-tests` required). Test support, not API.
#[doc(hidden)]
pub const PROPTEST_CASES: u32 = 1000;

/// Word with the most significant bit set.
pub const WORD_SIGNIFICANT_BIT: Word = WORD_MAX << (WORD_BIT_SIZE - 1);

/// Default precision in bits (two 64-bit words); the precision of the built-in [`crate::Consts`] cache.
pub const DEFAULT_P: usize = 128;

/// The size of exponent type in bits.
pub const EXPONENT_BIT_SIZE: usize = core::mem::size_of::<Exponent>() * 8;

/// Sign.
#[derive(PartialEq, Eq, Copy, Clone, Debug, Hash)]
pub enum Sign {
    /// Negative.
    Neg = -1,

    /// Positive.
    Pos = 1,
}

impl Sign {
    /// Changes the sign to the opposite.
    pub fn invert(&self) -> Self {
        match *self {
            Sign::Pos => Sign::Neg,
            Sign::Neg => Sign::Pos,
        }
    }

    /// Returns true if `self` is positive.
    pub fn is_positive(&self) -> bool {
        *self == Sign::Pos
    }

    /// Returns true if `self` is negative.
    pub fn is_negative(&self) -> bool {
        *self == Sign::Neg
    }

    /// Returns 1 for the positive sign and -1 for the negative sign.
    pub fn to_int(&self) -> i8 {
        *self as i8
    }
}

/// Possible errors.
#[derive(Debug, Clone, Copy)]
pub enum Error {
    /// The exponent value becomes greater than the upper limit of the range of exponent values.
    ExponentOverflow(Sign),

    /// Divizor is zero.
    DivisionByZero,

    /// Invalid argument.
    InvalidArgument,

    /// Correct-rounding retries exhausted (`MAX_PREC_RETRY`). Not a domain error.
    PrecisionRetryExhausted,

    /// Memory allocation error.
    MemoryAllocation,
}

#[cfg(feature = "std")]
impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None
    }
}

impl Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let repr = match self {
            Error::ExponentOverflow(s) => {
                if s.is_positive() {
                    "positive overflow"
                } else {
                    "negative overflow"
                }
            }
            Error::DivisionByZero => "division by zero",
            Error::InvalidArgument => "invalid argument",
            Error::PrecisionRetryExhausted => "precision retry exhausted",
            Error::MemoryAllocation => "memory allocation failure",
        };
        f.write_str(repr)
    }
}

impl PartialEq for Error {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::ExponentOverflow(l0), Self::ExponentOverflow(r0)) => l0 == r0,
            _ => core::mem::discriminant(self) == core::mem::discriminant(other),
        }
    }
}

impl From<TryReserveError> for Error {
    fn from(_: TryReserveError) -> Self {
        Error::MemoryAllocation
    }
}

impl From<LayoutError> for Error {
    fn from(_: LayoutError) -> Self {
        Error::MemoryAllocation
    }
}

/// Radix for parse/format (bases 2 through 36).
///
/// The four common bases are named variants; every other base is [`Radix::Other`], built with
/// [`Radix::try_new`]. `try_new(2)`, `try_new(8)`, `try_new(10)` and `try_new(16)` return the
/// named variants, so equal bases always compare equal.
#[derive(PartialEq, Eq, Copy, Clone, Debug, Hash)]
pub enum Radix {
    /// Binary (base 2).
    Bin,
    /// Octal (base 8).
    Oct,
    /// Decimal (base 10).
    Dec,
    /// Hexadecimal (base 16).
    Hex,
    /// Any other base in 2..=36; only [`Radix::try_new`] can build one.
    Other(OtherRadix),
}

/// A base in 2..=36 other than 2, 8, 10 and 16. Its field is private, so it is always valid.
#[derive(PartialEq, Eq, Copy, Clone, Debug, Hash)]
pub struct OtherRadix(u8);

impl Radix {
    /// Creates a radix in the inclusive range 2..=36.
    ///
    /// ## Errors
    ///
    ///  - InvalidArgument: `base` is outside 2..=36.
    pub fn try_new(base: u8) -> Result<Self, Error> {
        match base {
            2 => Ok(Radix::Bin),
            8 => Ok(Radix::Oct),
            10 => Ok(Radix::Dec),
            16 => Ok(Radix::Hex),
            3..=36 => Ok(Radix::Other(OtherRadix(base))),
            _ => Err(Error::InvalidArgument),
        }
    }

    /// Returns the numeric base.
    pub const fn value(self) -> u8 {
        match self {
            Radix::Bin => 2,
            Radix::Oct => 8,
            Radix::Dec => 10,
            Radix::Hex => 16,
            Radix::Other(OtherRadix(b)) => b,
        }
    }

    /// Returns `log2(base)` when the base is a power of two.
    pub const fn commensurable_shift(self) -> Option<usize> {
        let b = self.value();
        if b.is_power_of_two() {
            Some(b.trailing_zeros() as usize)
        } else {
            None
        }
    }

    /// Whether the scientific exponent must use `_e` (digit `e` appears in mantissa digits).
    pub const fn uses_underscore_exponent(self) -> bool {
        self.value() > 10
    }

    /// Approximate bits per digit (for buffer sizing).
    pub fn bits_per_digit(self) -> usize {
        match self.value() {
            2 => 1,
            8 => 3,
            10 => 3,
            16 => 4,
            b if b.is_power_of_two() => b.trailing_zeros() as usize,
            _ => 4,
        }
    }
}

impl From<Radix> for u8 {
    fn from(r: Radix) -> u8 {
        r.value()
    }
}

/// Rounding modes.
#[derive(Eq, PartialEq, Debug, Copy, Clone)]
pub enum RoundingMode {
    /// Skip rounding operation.
    None = 1,

    /// Round toward positive infinity (directed rounding).
    ///
    /// Every inexact result moves to the next representable value at or above the exact one.
    /// There is no tie case.
    Up = 2,

    /// Round toward negative infinity (directed rounding).
    ///
    /// Every inexact result moves to the next representable value at or below the exact one.
    /// There is no tie case.
    Down = 4,

    /// Round toward zero (truncation; directed rounding).
    ToZero = 8,

    /// Round away from zero (directed rounding).
    FromZero = 16,

    /// Round half to even.
    ToEven = 32,

    /// Round to nearest, ties to odd: the mirror of [`RoundingMode::ToEven`].
    ///
    /// This is **not** the sticky "round-to-odd" (jamming) mode that algorithms use to avoid double rounding
    /// (where any inexact result gets its last bit set to 1): only an exact tie goes to the odd neighbour.
    /// At `p = 4`, 17/16 is an exact tie and rounds to 1.000 under `ToEven` and to 1.001 under `ToOdd`;
    /// 33/32 is below the tie and rounds to 1.000 under both (sticky round-to-odd would give 1.001).
    ToOdd = 64,
}
