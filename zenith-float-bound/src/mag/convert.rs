//! The `Mag` ⇄ `ExactNum` bridge (CONTRACT 2.3, 3.3, 3.4; PLAN C5).
//!
//! Both directions go through zf's raw parts (R14): `ExactNum::as_raw_parts` gives the value as
//! `(-1)^negative * M * 2^(exponent - precision)` with `M` in little-endian 64-bit words, and `ExactNum::from_raw_parts`
//! builds one back. The reading side takes the top 64 bits of `M` as a `0.mant * 2^exp` `Mag` (round down) plus a
//! sticky flag (any lower bit set), and one ulp more is the round-up result when bits were dropped. The stored precision
//! is never assumed to round-trip: only the value is read.

use super::{Mag, TOP};
use crate::error::BoundError;
use zenith_float_num::{Error, ExactNum, RawKind, RawParts, EXPONENT_MAX, EXPONENT_MIN};

/// zf precision for `to_exact_num`: a `Mag` mantissa has 64 bits, so one 64-bit word is always exact.
const EXACT_PREC: usize = 64;

/// The top 64 bits of finite nonzero `|x|` as `(mant, sticky, exp)`: `0.mant * 2^exp` is `|x|` rounded down and
/// `sticky` says bits below it are set. `None` for zero, NaN and ±Inf (callers test NaN and ±Inf first). A `Mag`
/// exponent below `i64::MIN` (not reachable from zf's range; kept total) reads as `0.1 * 2^i64::MIN` when `upper`
/// (a bound above) and as `None` (zero, a bound below) otherwise.
fn top_of(x: &ExactNum, upper: bool) -> Option<(u64, bool, i64)> {
    let p = x.as_raw_parts();
    match p.kind {
        RawKind::NaN | RawKind::Inf | RawKind::Zero => None,
        RawKind::Finite => match top_of_words(p.words, p.exponent) {
            Ok(top) => top,
            Err(()) => upper.then_some((TOP, false, i64::MIN)),
        },
    }
}

/// `Ok(None)` when every word is zero, `Err(())` on exponent underflow.
fn top_of_words(words: &[u64], exponent: i64) -> Result<Option<(u64, bool, i64)>, ()> {
    // `skipped` zero words above the leading word, `idx` the leading word's index.
    let Some(skipped) = words.iter().rev().position(|w| *w != 0) else {
        return Ok(None);
    };
    let Some(idx) = words
        .len()
        .checked_sub(1)
        .and_then(|n| n.checked_sub(skipped))
    else {
        return Ok(None);
    };
    let Some(&hi) = words.get(idx) else {
        return Ok(None);
    };
    let lo = idx
        .checked_sub(1)
        .and_then(|i| words.get(i))
        .map_or(0, |w| *w);
    let rest_sticky = idx
        .checked_sub(1)
        .and_then(|i| words.get(..i))
        .is_some_and(|ws| ws.iter().any(|w| *w != 0));
    let lz = hi.leading_zeros();
    let window = (u128::from(hi).checked_shl(64).unwrap_or(0) | u128::from(lo))
        .checked_shl(lz)
        .unwrap_or(0);
    let mant = window
        .checked_shr(64)
        .and_then(|t| u64::try_from(t).ok())
        .unwrap_or(0);
    let sticky = rest_sticky || (window & u128::from(u64::MAX)) != 0;
    // Leading bit of M sits at 2^(exponent - 64 * skipped - lz - 1), so the `0.mant` exponent is this.
    let shift = i64::try_from(skipped)
        .ok()
        .and_then(|s| s.checked_mul(64))
        .and_then(|s| s.checked_add(i64::from(lz)));
    match shift.and_then(|s| exponent.checked_sub(s)) {
        Some(exp) => Ok(Some((mant, sticky, exp))),
        None => Err(()),
    }
}

impl Mag {
    /// Least `Mag` at or above `|x|` (CONTRACT 3.3). NaN and ±Inf give `Infinity`; zero gives zero.
    pub fn from_exact_num(x: &ExactNum) -> Mag {
        if x.is_nan() || x.is_inf() {
            return Mag::infinity();
        }
        let Some((mant, sticky, exp)) = top_of(x, true) else {
            return Mag::zero();
        };
        if !sticky {
            return Mag::from_canonical(exp, mant);
        }
        match mant.checked_add(1) {
            Some(m) => Mag::from_canonical(exp, m),
            None => match exp.checked_add(1) {
                Some(e) => Mag::from_canonical(e, TOP),
                None => Mag::infinity(),
            },
        }
    }

    /// Greatest `Mag` at or below `|x|` (CONTRACT 3.4). NaN and ±Inf give zero; zero gives zero.
    pub fn from_exact_num_lower(x: &ExactNum) -> Mag {
        match top_of(x, false) {
            Some((mant, _, exp)) => Mag::from_canonical(exp, mant),
            None => Mag::zero(),
        }
    }

    /// The exact value of `self` as an `ExactNum` (CONTRACT 2.3), built with `from_raw_parts` at 64 bits. Below
    /// `EXPONENT_MIN` the value is stored subnormal (mantissa shifted right, exponent `EXPONENT_MIN`) when no set bit
    /// is shifted out.
    ///
    /// # Errors
    /// `Representation` when `self` is `Infinity`, its exponent is above `EXPONENT_MAX`, or it is below
    /// `EXPONENT_MIN` and not exactly a zf subnormal.
    pub fn to_exact_num(&self) -> Result<ExactNum, BoundError> {
        const OP: &str = "Mag::to_exact_num";
        let Some((exp, mant)) = self.finite_parts() else {
            return Err(BoundError::Representation {
                op: OP,
                detail: "Infinity has no ExactNum value",
            });
        };
        if mant == 0 {
            return Ok(ExactNum::new(EXACT_PREC));
        }
        if exp > EXPONENT_MAX {
            return Err(BoundError::Representation {
                op: OP,
                detail: "exponent above zf EXPONENT_MAX",
            });
        }
        let word = match EXPONENT_MIN.checked_sub(exp).filter(|d| *d > 0) {
            None => mant,
            Some(d) => match u32::try_from(d) {
                Ok(d) if d < 64 && mant.trailing_zeros() >= d => mant.checked_shr(d).unwrap_or(0),
                Ok(_) | Err(_) => {
                    return Err(BoundError::Representation {
                        op: OP,
                        detail: "below zf EXPONENT_MIN and not an exact subnormal",
                    })
                }
            },
        };
        let words = [word];
        let parts = RawParts {
            kind: RawKind::Finite,
            words: &words,
            precision: EXACT_PREC,
            negative: false,
            exponent: exp.max(EXPONENT_MIN),
            inexact: false,
            op_inexact: false,
        };
        ExactNum::from_raw_parts(&parts).map_err(|e| map_zf_error(OP, e))
    }

    /// `(exp, mant)` of a finite `Mag`, including the canonical zero.
    fn finite_parts(&self) -> Option<(i64, u64)> {
        match (self.exp(), self.mant()) {
            (Some(e), Some(m)) => Some((e, m)),
            (None, Some(_)) | (Some(_), None) | (None, None) => None,
        }
    }
}

/// zf errors in zb terms (CONTRACT O2 table; PLAN: `PrecisionRetryExhausted` is `Resource`).
fn map_zf_error(op: &'static str, e: Error) -> BoundError {
    match e {
        Error::ExponentOverflow(_) => BoundError::Representation {
            op,
            detail: "zf exponent overflow",
        },
        Error::InvalidArgument => BoundError::Representation {
            op,
            detail: "zf invalid argument",
        },
        Error::DivisionByZero => BoundError::Domain {
            op,
            detail: "zf division by zero",
        },
        Error::PrecisionRetryExhausted => BoundError::Resource {
            op,
            detail: "zf precision retry exhausted",
        },
        Error::MemoryAllocation => BoundError::Resource {
            op,
            detail: "zf memory allocation",
        },
    }
}
