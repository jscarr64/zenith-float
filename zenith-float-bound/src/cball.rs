//! Slice CB: complex balls `CBall { re, im }`, a closed rectangle of two independent real [`Ball`] parts (CONTRACT 16.CB).
//! Built only on the public `Ball` API and S1; no radius arithmetic of its own. Ownership is O9 (16.CB.0): an op takes
//! `self` and borrows its second operand; a part is consumed where the inner `Ball` op consumes it.
//! `mul` (16.CB.4) and `div` (16.CB.7) are the four-product / scaled-`|b|²` forms, ported from Arb `acb_mul` / `acb_div`.

use zenith_float_num::ExactNum;

use std::ops::Neg;

use crate::ball::Ball;
use crate::enclosure::Enclosure;
use crate::error::{BoundError, BoundResult};
use crate::precision::Precision;

/// A complex ball: the rectangle `re × im` (16.CB.1). Both parts share one prec (I2); the parts are independent (I3).
#[derive(Debug)]
pub struct CBall {
    re: Ball,
    im: Ball,
}

impl CBall {
    /// `re + i·im`. `InvalidArgument` when the parts' precs differ (I2).
    pub fn new(re: Ball, im: Ball) -> BoundResult<CBall> {
        if re.prec() != im.prec() {
            return Err(BoundError::InvalidArgument {
                op: "cball_new",
                detail: "re and im must have the same prec",
            });
        }
        Ok(CBall { re, im })
    }

    /// `re + 0i`, the imaginary part exact zero at `re`'s prec.
    pub fn from_real(re: Ball) -> BoundResult<CBall> {
        let im = Ball::exact(ExactNum::from(0i64), re.prec())?;
        Ok(CBall { re, im })
    }

    /// Both parts all reals.
    pub fn all_complex(prec: Precision) -> BoundResult<CBall> {
        Ok(CBall {
            re: Ball::all_reals(prec)?,
            im: Ball::all_reals(prec)?,
        })
    }

    /// The real part.
    pub fn re(&self) -> &Ball {
        &self.re
    }

    /// The imaginary part.
    pub fn im(&self) -> &Ball {
        &self.im
    }

    /// The shared prec of both parts.
    pub fn prec(&self) -> Precision {
        self.re.prec()
    }

    /// Both parts identical (same mid, rad and prec).
    pub fn identical(&self, other: &CBall) -> bool {
        self.re.identical(&other.re) && self.im.identical(&other.im)
    }

    /// `(re, −im)`, exact.
    pub fn conj(self) -> CBall {
        CBall {
            re: self.re,
            im: self.im.negated(),
        }
    }

    /// Part-wise `Ball::add`.
    pub fn add(self, other: &CBall, prec: Precision) -> BoundResult<CBall> {
        Ok(CBall {
            re: self.re.add(&other.re, prec)?,
            im: self.im.add(&other.im, prec)?,
        })
    }

    /// Part-wise `Ball::sub`.
    pub fn sub(self, other: &CBall, prec: Precision) -> BoundResult<CBall> {
        Ok(CBall {
            re: self.re.sub(&other.re, prec)?,
            im: self.im.sub(&other.im, prec)?,
        })
    }

    /// `|a|` as a real ball: `Ball::hypot(re, im)` as is (16.CB.6).
    pub fn abs(self, prec: Precision) -> BoundResult<Ball> {
        self.re.hypot(&self.im, prec)
    }

    /// `a²` (16.CB.5): `re = sqr(re) − sqr(im)` by `Ball::sqr`, `im = 2·re·im` by one `Ball::mul` and an exact
    /// `mul_2exp(·, 1)`. Each part is no wider than the same part of `mul(a, a)`: `Ball::sqr` is no wider than
    /// `Ball::mul(x, x)`, and the doubling replaces a rounded `add` of two equal products.
    pub fn sqr(self, prec: Precision) -> BoundResult<CBall> {
        let rr = self.re.sqr(prec)?;
        let ii = self.im.sqr(prec)?;
        let re = rr.sub(&ii, prec)?;
        let im = self.re.mul(&self.im, prec)?.mul_2exp(1, prec)?;
        CBall::new(re, im)
    }

    /// Four-product `a·b` (16.CB.4 / Arb `acb_mul`): `re = a.re·b.re − a.im·b.im`,
    /// `im = a.re·b.im + a.im·b.re`, each product `Ball::mul`, each combination `Ball::sub`/`add`.
    pub fn mul(self, other: &CBall, prec: Precision) -> BoundResult<CBall> {
        const OP: &str = "cball_mul";
        let (ar, ai) = (self.re, self.im);
        let ar2 = ball_again(&ar).map_err(|e| remap(e, OP))?;
        let ai2 = ball_again(&ai).map_err(|e| remap(e, OP))?;
        let rr = ar.mul(&other.re, prec).map_err(|e| remap(e, OP))?;
        let ii = ai.mul(&other.im, prec).map_err(|e| remap(e, OP))?;
        let re = rr.sub(&ii, prec).map_err(|e| remap(e, OP))?;
        let ri = ar2.mul(&other.im, prec).map_err(|e| remap(e, OP))?;
        let ir = ai2.mul(&other.re, prec).map_err(|e| remap(e, OP))?;
        let im = ri.add(&ir, prec).map_err(|e| remap(e, OP))?;
        CBall::new(re, im).map_err(|e| remap(e, OP))
    }

    /// `a / b` (16.CB.7 / Arb `acb_div`): scale by the larger part exponent of `b`, then
    /// `d = sqr(b.re) + sqr(b.im)`, `re = (a.re·b.re + a.im·b.im) / d`,
    /// `im = (a.im·b.re − a.re·b.im) / d`. R1 and R2 return both parts all reals, never `Domain`.
    pub fn div(self, other: &CBall, prec: Precision) -> BoundResult<CBall> {
        const OP: &str = "cball_div";
        if reaches_zero(&other.re) && reaches_zero(&other.im) {
            return CBall::all_complex(prec).map_err(|e| remap(e, OP));
        }
        let k = scale_exp(&other.re).max(scale_exp(&other.im));
        let nk = k.checked_neg().unwrap_or(i64::MIN);
        let ar = self.re.mul_2exp(nk, prec).map_err(|e| remap(e, OP))?;
        let ai = self.im.mul_2exp(nk, prec).map_err(|e| remap(e, OP))?;
        let br = ball_again(&other.re)
            .and_then(|b| b.mul_2exp(nk, prec))
            .map_err(|e| remap(e, OP))?;
        let bi = ball_again(&other.im)
            .and_then(|b| b.mul_2exp(nk, prec))
            .map_err(|e| remap(e, OP))?;
        let d = br
            .sqr(prec)
            .and_then(|dr| bi.sqr(prec).and_then(|di| dr.add(&di, prec)))
            .map_err(|e| remap(e, OP))?;
        if reaches_zero(&d) {
            return CBall::all_complex(prec).map_err(|e| remap(e, OP));
        }
        let ar2 = ball_again(&ar).map_err(|e| remap(e, OP))?;
        let ai2 = ball_again(&ai).map_err(|e| remap(e, OP))?;
        let d2 = ball_again(&d).map_err(|e| remap(e, OP))?;
        let num_re = ar
            .mul(&br, prec)
            .and_then(|t| ai.mul(&bi, prec).and_then(|u| t.add(&u, prec)))
            .map_err(|e| remap(e, OP))?;
        let num_im = ai2
            .mul(&br, prec)
            .and_then(|t| ar2.mul(&bi, prec).and_then(|u| t.sub(&u, prec)))
            .map_err(|e| remap(e, OP))?;
        let re = num_re.div(&d, prec).map_err(|e| remap(e, OP))?;
        let im = num_im.div(&d2, prec).map_err(|e| remap(e, OP))?;
        CBall::new(re, im).map_err(|e| remap(e, OP))
    }
}

/// Rebuild a `Ball` from its public mid, rad and prec so a part can be used in two `Ball` ops
/// (`Ball::{add,sub,mul,div,mul_2exp}` still take the first operand by value).
fn ball_again(b: &Ball) -> BoundResult<Ball> {
    Ball::from_mid_rad(b.mid().clone(), mag_again(b.rad()), b.prec())
}

fn mag_again(r: &crate::mag::Mag) -> crate::mag::Mag {
    if r.is_infinity() {
        return crate::mag::Mag::infinity();
    }
    if r.is_zero() {
        return crate::mag::Mag::zero();
    }
    crate::mag::Mag::finite(
        r.mant().expect("finite mag mantissa"),
        r.exp().expect("finite mag exponent"),
    )
    .expect("canonical finite mag is valid")
}

/// `|m| ≤ r`, decided by `Mag::cmp_exact_num`, including an infinite radius (16.CB.7 R1/R2).
fn reaches_zero(b: &Ball) -> bool {
    match b.rad().cmp_exact_num(b.mid()) {
        Some(core::cmp::Ordering::Less) => false,
        Some(core::cmp::Ordering::Equal | core::cmp::Ordering::Greater) | None => true,
    }
}

fn scale_exp(b: &Ball) -> i64 {
    match b.mid().exponent() {
        Some(e) if !b.mid().is_zero() => e,
        _ => 0,
    }
}

fn remap(err: BoundError, op: &'static str) -> BoundError {
    match err {
        BoundError::Representation { detail, .. } => BoundError::Representation { op, detail },
        BoundError::InvalidArgument { detail, .. } => BoundError::InvalidArgument { op, detail },
        BoundError::Domain { detail, .. } => BoundError::Domain { op, detail },
        BoundError::Resource { detail, .. } => BoundError::Resource { op, detail },
        BoundError::NotPositiveDefinite { detail, .. } => {
            BoundError::NotPositiveDefinite { op, detail }
        }
        BoundError::Singular { detail, .. } => BoundError::Singular { op, detail },
        BoundError::NotConverged { detail, .. } => BoundError::NotConverged { op, detail },
        BoundError::DimensionMismatch { detail, .. } => {
            BoundError::DimensionMismatch { op, detail }
        }
        BoundError::NotSymmetric { detail, .. } => BoundError::NotSymmetric { op, detail },
    }
}

impl Neg for CBall {
    type Output = CBall;

    /// `(−re, −im)`: exact and infallible (each part via `Ball`'s exact negation, CONTRACT 16.CB.3).
    fn neg(self) -> CBall {
        CBall {
            re: self.re.negated(),
            im: self.im.negated(),
        }
    }
}
