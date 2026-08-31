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

/// Maximum exponent value.
#[cfg(not(target_pointer_width = "32"))]
pub const EXPONENT_MAX: Exponent = Exponent::MAX;

/// Maximum exponent value.
#[cfg(target_pointer_width = "32")]
pub const EXPONENT_MAX: Exponent = Exponent::MAX / 4;

/// Minimum exponent value.
#[cfg(not(target_pointer_width = "32"))]
pub const EXPONENT_MIN: Exponent = Exponent::MIN;

/// Minimum exponent value.
#[cfg(target_pointer_width = "32")]
pub const EXPONENT_MIN: Exponent = Exponent::MIN / 4;

/// Maximum value of a word.
pub const WORD_MAX: Word = Word::MAX;

/// Base of words.
pub const WORD_BASE: DoubleWord = WORD_MAX as DoubleWord + 1;

/// Size of a word in bits.
pub const WORD_BIT_SIZE: usize = core::mem::size_of::<Word>() * 8;

/// Word with the most significant bit set.
pub const WORD_SIGNIFICANT_BIT: Word = WORD_MAX << (WORD_BIT_SIZE - 1);

/// Default precision.
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
#[derive(PartialEq, Eq, Copy, Clone, Debug, Hash)]
pub struct Radix(u8);

#[allow(non_upper_case_globals)]
impl Radix {
    /// Binary (base 2).
    pub const Bin: Radix = Radix(2);
    /// Octal (base 8).
    pub const Oct: Radix = Radix(8);
    /// Decimal (base 10).
    pub const Dec: Radix = Radix(10);
    /// Hexadecimal (base 16).
    pub const Hex: Radix = Radix(16);

    /// Creates a radix in the inclusive range 2..=36.
    ///
    /// ## Errors
    ///
    ///  - InvalidArgument: `base` is outside 2..=36.
    pub fn try_new(base: u8) -> Result<Self, Error> {
        if (2..=36).contains(&base) {
            Ok(Radix(base))
        } else {
            Err(Error::InvalidArgument)
        }
    }

    /// Returns the numeric base.
    pub const fn value(self) -> u8 {
        self.0
    }

    /// Returns `log2(base)` when the base is a power of two.
    pub const fn commensurable_shift(self) -> Option<usize> {
        let b = self.0;
        if b.is_power_of_two() {
            Some(b.trailing_zeros() as usize)
        } else {
            None
        }
    }

    /// Whether the scientific exponent must use `_e` (digit `e` appears in mantissa digits).
    pub const fn uses_underscore_exponent(self) -> bool {
        self.0 > 10
    }

    /// Approximate bits per digit (for buffer sizing).
    pub fn bits_per_digit(self) -> usize {
        match self.0 {
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
        r.0
    }
}

/// Rounding modes.
#[derive(Eq, PartialEq, Debug, Copy, Clone)]
pub enum RoundingMode {
    /// Skip rounding operation.
    None = 1,

    /// Round half toward positive infinity.
    Up = 2,

    /// Round half toward negative infinity.
    Down = 4,

    /// Round half toward zero.
    ToZero = 8,

    /// Round half away from zero.
    FromZero = 16,

    /// Round half to even.
    ToEven = 32,

    /// Round half to odd.
    ToOdd = 64,
}
