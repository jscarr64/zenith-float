//! The sealed `Enclosure` trait (PLAN C5, 2026-10-23): what every enclosure type answers.
//!
//! Implemented by [`Ball`]; `ComplexBall` implements it when it lands (L2a). The trait is sealed:
//! no type outside zenith-bound can implement it, so methods may be added in a minor release.

use zenith_float_num::{ExactNum, RoundingMode};

use crate::ball::Ball;
use crate::mag::Mag;
use crate::precision::Precision;
use crate::tri::Tri;

mod sealed {
    /// Only zenith-bound types implement `Enclosure`.
    pub trait Sealed {}
}

/// Extra bits over the wider operand used by [`Ball`]'s `contains` difference.
const GUARD_BITS: usize = 128;

/// zf precision used when an operand reports none.
const MIN_PREC: usize = 64;

/// A set-valued number with a midpoint precision.
pub trait Enclosure: sealed::Sealed {
    /// The point type this enclosure can contain.
    type Point;

    /// The requested midpoint precision in bits.
    fn prec(&self) -> Precision;

    /// Structural identity: same midpoint value, same radius, same `prec` (CONTRACT 6.5 b, O20).
    fn identical(&self, other: &Self) -> bool;

    /// Whether `x` lies in the enclosure: `True`/`False` only when decided, else `Unknown`.
    fn contains(&self, x: &Self::Point) -> Tri;
}

impl sealed::Sealed for Ball {}

impl Enclosure for Ball {
    type Point = ExactNum;

    fn prec(&self) -> Precision {
        Ball::prec(self)
    }

    fn identical(&self, other: &Ball) -> bool {
        self.prec() == other.prec() && self.rad() == other.rad() && self.mid() == other.mid()
    }

    fn contains(&self, x: &ExactNum) -> Tri {
        if x.is_nan() || x.is_inf() {
            return Tri::False;
        }
        if self.rad().is_infinity() {
            return Tri::True;
        }
        let p = match working_prec(self.mid(), x) {
            Some(p) => p,
            None => return Tri::Unknown,
        };
        // d = mid - x, enclosed by [lo, hi]: both steps round in one direction and are monotone.
        let hi = ExactNum::new(p)
            .sub(x, p, RoundingMode::Up)
            .add(self.mid(), p, RoundingMode::Up);
        let lo =
            ExactNum::new(p)
                .sub(x, p, RoundingMode::Down)
                .add(self.mid(), p, RoundingMode::Down);
        if hi.is_nan() || lo.is_nan() {
            return Tri::Unknown;
        }
        let upper = Mag::from_exact_num(&hi).max(Mag::from_exact_num(&lo));
        if upper <= *self.rad() {
            return Tri::True;
        }
        // Tight: a rad zenith-float can hold is exact at p (64-bit mantissa, p >= MIN_PREC + GUARD_BITS), and directed
        // rounding is monotone and fixes representable values, so hi > rad (or lo < -rad) decides |mid - x| > rad.
        match self.rad().to_exact_num() {
            Ok(_) => Tri::False,
            Err(_) => Tri::Unknown,
        }
    }
}

/// Precision for the exact-when-possible difference of `a` and `b`.
fn working_prec(a: &ExactNum, b: &ExactNum) -> Option<usize> {
    let pa = a.precision().map_or(MIN_PREC, |p| p.max(MIN_PREC));
    let pb = b.precision().map_or(MIN_PREC, |p| p.max(MIN_PREC));
    pa.max(pb).checked_add(GUARD_BITS)
}
