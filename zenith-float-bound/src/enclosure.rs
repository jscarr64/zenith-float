//! The sealed `Enclosure` trait (PLAN C5, 2026-10-23): what every enclosure type answers.
//!
//! Implemented by [`Ball`]; `ComplexBall` implements it when it lands (L2a). The trait is sealed:
//! no type outside zenith-float-bound can implement it, so methods may be added in a minor release.

use zenith_float_num::{ExactNum, RoundingMode};

use crate::ball::Ball;
use crate::mag::Mag;
use crate::precision::Precision;
use crate::tri::Tri;

mod sealed {
    /// Only zenith-float-bound types implement `Enclosure`.
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
        both_within(within(&hi, self.rad()), within(&lo, self.rad()), self.rad())
    }
}

/// Whether one endpoint `d` of the enclosure of `mid - x` lies within `rad`: `Unknown` when zf could not compute
/// it (a NaN), otherwise decided by `|d| <= rad`.
fn within(d: &ExactNum, rad: &Mag) -> Tri {
    if d.is_nan() {
        Tri::Unknown
    } else if Mag::from_exact_num(d) <= *rad {
        Tri::True
    } else {
        Tri::False
    }
}

/// `x` is in the ball when both endpoints of `mid - x` are within `rad`. An undecided endpoint leaves the answer
/// undecided whatever the other one says. An endpoint outside decides `False` only when `rad` is exact.
/// Tight: a rad zenith-float can hold is exact at p (64-bit mantissa, p >= MIN_PREC + GUARD_BITS), and directed
/// rounding is monotone and fixes representable values, so hi > rad (or lo < -rad) decides |mid - x| > rad.
fn both_within(hi: Tri, lo: Tri, rad: &Mag) -> Tri {
    match (hi, lo) {
        (Tri::True, Tri::True) => Tri::True,
        (Tri::Unknown, Tri::True | Tri::False | Tri::Unknown)
        | (Tri::True | Tri::False, Tri::Unknown) => Tri::Unknown,
        (Tri::False, Tri::True | Tri::False) | (Tri::True, Tri::False) => {
            match rad.to_exact_num() {
                Ok(_) => Tri::False,
                Err(_) => Tri::Unknown,
            }
        }
    }
}

/// Precision for the exact-when-possible difference of `a` and `b`.
fn working_prec(a: &ExactNum, b: &ExactNum) -> Option<usize> {
    let pa = a.precision().map_or(MIN_PREC, |p| p.max(MIN_PREC));
    let pb = b.precision().map_or(MIN_PREC, |p| p.max(MIN_PREC));
    pa.max(pb).checked_add(GUARD_BITS)
}

/// Rusty (2026-11-10): every endpoint verdict and every pair of them. A one-sided failure (one endpoint
/// undecided, the other decided) cannot be produced through `contains` on demand (zf fails on both endpoints
/// alike), so the combination is checked here directly; the old single `hi.is_nan() || lo.is_nan()` test had a
/// surviving `&&` mutant exactly there.
#[cfg(test)]
mod contains_endpoints {
    use super::{both_within, within};
    use crate::mag::Mag;
    use crate::tri::Tri;
    use zenith_float_num::{ExactNum, RoundingMode};

    const TOP: u64 = 0x8000_0000_0000_0000;

    /// The nine (hi, lo) endpoint pairs with the verdict for an exact rad and for a rad zf cannot hold.
    /// Rows hi = True, False, Unknown; columns lo the same.
    fn table() -> [(Tri, Tri, Tri, Tri); 9] {
        [
            (Tri::True, Tri::True, Tri::True, Tri::True),
            (Tri::True, Tri::False, Tri::False, Tri::Unknown),
            (Tri::True, Tri::Unknown, Tri::Unknown, Tri::Unknown),
            (Tri::False, Tri::True, Tri::False, Tri::Unknown),
            (Tri::False, Tri::False, Tri::False, Tri::Unknown),
            (Tri::False, Tri::Unknown, Tri::Unknown, Tri::Unknown),
            (Tri::Unknown, Tri::True, Tri::Unknown, Tri::Unknown),
            (Tri::Unknown, Tri::False, Tri::Unknown, Tri::Unknown),
            (Tri::Unknown, Tri::Unknown, Tri::Unknown, Tri::Unknown),
        ]
    }

    #[test]
    fn every_pair_of_endpoint_verdicts() -> Result<(), String> {
        // 2^-1 is exact in zf; 2^(2^62 - 1) is above zf's exponent range, so it has no exact zf value.
        let exact = Mag::finite(TOP, 0).map_err(|e| format!("{e:?}"))?;
        let huge = Mag::finite(TOP, (1i64 << 62) + 1).map_err(|e| format!("{e:?}"))?;
        if huge.to_exact_num().is_ok() {
            return Err(String::from("the huge rad should not be an exact zf value"));
        }
        for (rad, rad_exact) in [(&exact, true), (&huge, false), (&Mag::infinity(), false)] {
            for (hi, lo, if_exact, if_not) in table() {
                let what = format!("({hi:?}, {lo:?}) rad {rad:?}");
                let w = if rad_exact { if_exact } else { if_not };
                let got = both_within(hi, lo, rad);
                if got != w {
                    return Err(format!("{what}: got {got:?}, want {w:?}"));
                }
            }
        }
        Ok(())
    }

    fn q(num: i64, den: i64) -> ExactNum {
        ExactNum::from_i64(num, 64).div(&ExactNum::from_i64(den, 64), 64, RoundingMode::ToEven)
    }

    #[test]
    fn one_endpoint_verdict() -> Result<(), String> {
        // rad = 1/2 (Mag 0.5 * 2^0): |d| below, at and above it, both signs; 0; Inf of both signs; NaN.
        let rad = Mag::finite(TOP, 0).map_err(|e| format!("{e:?}"))?;
        let rows = [
            (q(0, 1), Tri::True),
            (q(1, 4), Tri::True),
            (q(-1, 4), Tri::True),
            (q(1, 2), Tri::True),
            (q(-1, 2), Tri::True),
            (q(1, 1), Tri::False),
            (q(-1, 1), Tri::False),
            (q(1, 0), Tri::False),
            (q(-1, 0), Tri::False),
            (ExactNum::nan(None), Tri::Unknown),
        ];
        for (d, w) in rows {
            let got = within(&d, &rad);
            if got != w {
                return Err(format!("within({d:?}, 1/2): got {got:?}, want {w:?}"));
            }
        }
        Ok(())
    }
}
