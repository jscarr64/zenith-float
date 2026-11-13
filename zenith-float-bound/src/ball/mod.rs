//! `Ball`: the closed set `[mid - rad, mid + rad]` (CONTRACT 6.1, 6.5 D15a/e, 8.1, 15.1-15.3).
//!
//! Fields are private (D15a). Constructors are checked: `prec == 0` is `InvalidArgument` (O17-W), a NaN or infinite
//! midpoint is `Representation` (6.2), a −0 midpoint becomes +0 (6.3), and an `Infinity` radius gives the canonical
//! all-reals ball with `mid = +0` (6.1). A constructor never rounds: an exact input wider than `prec` keeps its value
//! and its radius (D15a), so nothing is silently truncated. `Ball` is not `Clone`, `Copy`, `PartialEq` or `PartialOrd`
//! (15.2).

use crate::error::{BoundError, BoundResult};
use crate::mag::Mag;
use crate::precision::{Precision, MAX_PRECISION, MIN_PRECISION};
use zenith_float_num::ExactNum;

/// zf precision of a canonical +0 midpoint (any value works: zero has no mantissa bits; 64 is zf's smallest container).
const ZERO_PREC: usize = 64;

/// A midpoint-radius ball with the requested midpoint precision `prec` (bits, at least 1).
#[derive(Debug)]
pub struct Ball {
    mid: ExactNum,
    rad: Mag,
    prec: Precision,
}

impl Ball {
    /// The exact point `mid ± 0` (D15a: the midpoint is kept exactly even when wider than `prec`).
    ///
    /// # Errors
    /// `InvalidArgument` when `prec == 0`; `Representation` when `mid` is NaN or infinite.
    pub fn exact(mid: ExactNum, prec: Precision) -> BoundResult<Ball> {
        check_prec("Ball::exact", prec)?;
        Ball::from_mid_rad(mid, Mag::zero(), prec)
    }

    /// The ball `mid ± rad`; an `Infinity` radius gives the all-reals ball (`mid = +0`).
    ///
    /// # Errors
    /// `InvalidArgument` when `prec == 0`; `Representation` when `mid` is NaN or infinite.
    pub fn from_mid_rad(mid: ExactNum, rad: Mag, prec: Precision) -> BoundResult<Ball> {
        check_prec("Ball::from_mid_rad", prec)?;
        if mid.is_nan() || mid.is_inf() {
            return Err(BoundError::Representation {
                op: "Ball::from_mid_rad",
                detail: "NaN or infinite midpoint",
            });
        }
        if rad.is_infinity() {
            return Ball::all_reals(prec);
        }
        Ok(Ball {
            mid: positive_zero(mid),
            rad,
            prec,
        })
    }

    /// The all-reals ball `+0 ± Infinity` (CONTRACT 6.1, 6.4).
    ///
    /// # Errors
    /// `InvalidArgument` when `prec == 0`.
    pub fn all_reals(prec: Precision) -> BoundResult<Ball> {
        check_prec("Ball::all_reals", prec)?;
        Ok(Ball {
            mid: ExactNum::new(ZERO_PREC),
            rad: Mag::infinity(),
            prec,
        })
    }

    /// The midpoint.
    pub fn mid(&self) -> &ExactNum {
        &self.mid
    }

    /// The radius.
    pub fn rad(&self) -> &Mag {
        &self.rad
    }

    /// The requested midpoint precision (never read back from zenith-float, CONTRACT 8.1).
    pub fn prec(&self) -> Precision {
        self.prec
    }

    /// False iff the radius is `Infinity` (CONTRACT 6.4, 15.3).
    pub fn is_informative(&self) -> bool {
        !self.rad.is_infinity()
    }

    /// True iff the radius is zero, so the ball is the single point `mid`.
    pub fn is_exact(&self) -> bool {
        self.rad.is_zero()
    }

    /// Moves the parts out: `(mid, rad, prec)`.
    pub fn into_parts(self) -> (ExactNum, Mag, Precision) {
        (self.mid, self.rad, self.prec)
    }

    /// `-self`, exact and infallible (CONTRACT 15.1): the midpoint negated with −0 normalised to +0 (6.3);
    /// radius and prec unchanged. Crate-only: `impl Neg for Ball` lives in `arith` (slice S1).
    pub(crate) fn negated(mut self) -> Ball {
        self.mid.inv_sign();
        self.mid = positive_zero(self.mid);
        self
    }
}

/// A prec outside `MIN_PRECISION..=MAX_PRECISION` is rejected before zenith-float sees it (CONTRACT 8.1, O17-W).
fn check_prec(op: &'static str, prec: Precision) -> Result<(), BoundError> {
    if prec < MIN_PRECISION {
        Err(BoundError::InvalidArgument {
            op,
            detail: "prec must be at least MIN_PRECISION (1)",
        })
    } else if prec > MAX_PRECISION {
        Err(BoundError::InvalidArgument {
            op,
            detail: "prec must be at most MAX_PRECISION",
        })
    } else {
        Ok(())
    }
}

/// A −0 midpoint becomes +0 (CONTRACT 6.3); every other value is returned unchanged.
fn positive_zero(mid: ExactNum) -> ExactNum {
    if mid.is_zero() && mid.is_negative() {
        ExactNum::new(ZERO_PREC)
    } else {
        mid
    }
}
