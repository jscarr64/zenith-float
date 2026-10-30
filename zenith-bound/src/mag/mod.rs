//! Integer-only magnitude, an upper bound on a nonnegative real.
//!
//! A finite nonzero value with exponent `exp` and mantissa `mant` (top bit set) denotes
//! `mant * 2^(exp - 64)`, so the value lies in `[2^(exp - 1), 2^exp)` (CONTRACT 1.2).
//! The canonical zero is `exp == 0`, `mant == 0`, and it is the only finite value whose
//! mantissa is below `2^63` (CONTRACT 2.1). `Infinity` means that no finite bound exists.
//!
//! Fields are private. Values are built only with [`Mag::finite`], [`Mag::zero`] and
//! [`Mag::infinity`]. `Mag` is not `Clone`, `Copy` or `Hash`. Ordering is a total order:
//! zero is least, `Infinity` is greatest, and canonical finite values compare by `(exp, mant)`.
//!
//! Round-up operations return the least representable `Mag` at or above the exact result
//! (CONTRACT 14.1). Round-down operations return the greatest `Mag` at or below it, or zero
//! when that result is not representable (CONTRACT 14.2, 3.2). Exponent overflow is
//! `Infinity` on the way up and the largest finite `Mag` on the way down (CONTRACT 3.3).

mod add;
mod cmp;
mod convert;
mod div;
mod lower;
mod mul;
mod round;
mod sqrt;
mod ulp;

use crate::error::BoundError;
use core::num::NonZeroU64;

/// Least canonical nonzero mantissa, `2^63`.
const TOP: u64 = 0x8000_0000_0000_0000;

#[derive(Debug)]
enum Inner {
    Finite { exp: i64, mant: u64 },
    Infinity,
}

/// Upper bound on a nonnegative real, or `Infinity` when no finite bound exists.
#[derive(Debug)]
pub struct Mag {
    inner: Inner,
}

impl Mag {
    /// The canonical zero (`exp == 0`, `mant == 0`).
    pub fn zero() -> Mag {
        Mag {
            inner: Inner::Finite { exp: 0, mant: 0 },
        }
    }

    /// The infinite bound: every finite `Mag` is strictly smaller.
    pub fn infinity() -> Mag {
        Mag {
            inner: Inner::Infinity,
        }
    }

    /// Builds the finite magnitude denoted by `mant * 2^(exp - 64)`.
    ///
    /// A mantissa in `1..2^63` is normalised by a left shift of `lz` leading zeros, and the
    /// exponent becomes `exp - lz`. That adjustment is exact. `mant == 0` is accepted only
    /// together with `exp == 0`. If `exp - lz` would be below `i64::MIN`, the value is not
    /// representable and this returns [`BoundError::Representation`] without rounding.
    pub fn finite(mant: u64, exp: i64) -> Result<Mag, BoundError> {
        match NonZeroU64::new(mant) {
            None if exp == 0 => Ok(Mag::zero()),
            None => Err(BoundError::Representation {
                op: "Mag::finite",
                detail: "zero mantissa with nonzero exponent",
            }),
            Some(nonzero) => {
                // A nonzero `u64` has at most 63 leading zeros, so the shift below drops only
                // zero bits and loses nothing.
                let shift = nonzero.leading_zeros();
                // One error arm for both steps: the `i128` subtraction cannot fail (it is below
                // 2^65 in magnitude), and `i64::try_from` fails exactly when `exp - shift`
                // is below `i64::MIN`, which is the documented representation error.
                let exp = i128::from(exp)
                    .checked_sub(i128::from(shift))
                    .and_then(|adjusted| i64::try_from(adjusted).ok())
                    .ok_or(BoundError::Representation {
                        op: "Mag::finite",
                        detail: "below the smallest positive Mag",
                    })?;
                Ok(Mag::from_canonical(exp, mant << shift))
            }
        }
    }

    /// Exponent of a finite value. `None` for `Infinity`.
    pub fn exp(&self) -> Option<i64> {
        match &self.inner {
            Inner::Finite { exp, .. } => Some(*exp),
            Inner::Infinity => None,
        }
    }

    /// Mantissa of a finite value. `None` for `Infinity`.
    pub fn mant(&self) -> Option<u64> {
        match &self.inner {
            Inner::Finite { mant, .. } => Some(*mant),
            Inner::Infinity => None,
        }
    }

    /// True for the canonical zero.
    pub fn is_zero(&self) -> bool {
        matches!(&self.inner, Inner::Finite { mant: 0, .. })
    }

    /// True for `Infinity`.
    pub fn is_infinity(&self) -> bool {
        matches!(&self.inner, Inner::Infinity)
    }

    fn from_canonical(exp: i64, mant: u64) -> Mag {
        Mag {
            inner: Inner::Finite { exp, mant },
        }
    }

    fn min_pos() -> Mag {
        Mag::from_canonical(i64::MIN, TOP)
    }

    fn max_finite() -> Mag {
        Mag::from_canonical(i64::MAX, u64::MAX)
    }

    /// `(exp, mant)` when this value is finite and nonzero; the type carries the nonzero mantissa.
    fn nonzero(&self) -> Option<(i64, NonZeroU64)> {
        match &self.inner {
            Inner::Finite { exp, mant } => NonZeroU64::new(*mant).map(|m| (*exp, m)),
            Inner::Infinity => None,
        }
    }
}
