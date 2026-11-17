//! Slice CB: complex balls `CBall { re, im }`, a closed rectangle of two independent real [`Ball`] parts (CONTRACT 16.CB).
//! Built only on the public `Ball` API and S1; no radius arithmetic of its own. Ownership is O9 (16.CB.0): an op takes
//! `self` and borrows its second operand; a part is consumed where the inner `Ball` op consumes it.
//! `mul` (16.CB.4) and `div` (16.CB.7) are the four-product / scaled-`|b|²` forms, ported from Arb `acb_mul` / `acb_div`.

use zenith_float_num::{ExactNum, RoundingMode};

use std::ops::Neg;

use crate::ball::Ball;
use crate::enclosure::Enclosure;
use crate::error::{BoundError, BoundResult};
use crate::precision::{Precision, MAX_PRECISION};

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

    /// `a / b` (16.CB.7 / Arb `acb_div`): scale by the larger nonzero part exponent of `b`, then
    /// `d` is a positive enclosure of `|b|²` (mid-rad, or endpoint distances when that reaches 0),
    /// `re = (a.re·b.re + a.im·b.im) / d`, `im = (a.im·b.re − a.re·b.im) / d`.
    /// R1 (the divisor rectangle contains 0) returns both parts all reals, never `Domain`.
    /// Mid-rad `sqr` of a part may reach 0 even when the rectangle excludes 0; that is not
    /// all-reals. All reals only if a positive lower bound on `|b|²` cannot be certified
    /// (or the reciprocal overflows).
    pub fn div(self, other: &CBall, prec: Precision) -> BoundResult<CBall> {
        const OP: &str = "cball_div";
        if reaches_zero(&other.re) && reaches_zero(&other.im) {
            return CBall::all_complex(prec).map_err(|e| remap(e, OP));
        }
        let k = scale_k(&other.re, &other.im);
        let nk = k.checked_neg().unwrap_or(i64::MIN);
        let ar = self.re.mul_2exp(nk, prec).map_err(|e| remap(e, OP))?;
        let ai = self.im.mul_2exp(nk, prec).map_err(|e| remap(e, OP))?;
        let br = ball_again(&other.re)
            .and_then(|b| b.mul_2exp(nk, prec))
            .map_err(|e| remap(e, OP))?;
        let bi = ball_again(&other.im)
            .and_then(|b| b.mul_2exp(nk, prec))
            .map_err(|e| remap(e, OP))?;
        let d_mrad = ball_again(&br)
            .and_then(|x| x.sqr(prec))
            .and_then(|dr| {
                ball_again(&bi)
                    .and_then(|y| y.sqr(prec))
                    .and_then(|di| dr.add(&di, prec))
            })
            .map_err(|e| remap(e, OP))?;
        let ar2 = ball_again(&ar).map_err(|e| remap(e, OP))?;
        let ai2 = ball_again(&ai).map_err(|e| remap(e, OP))?;
        let num_re = ar
            .mul(&br, prec)
            .and_then(|t| ai.mul(&bi, prec).and_then(|u| t.add(&u, prec)))
            .map_err(|e| remap(e, OP))?;
        let num_im = ai2
            .mul(&br, prec)
            .and_then(|t| ar2.mul(&bi, prec).and_then(|u| t.sub(&u, prec)))
            .map_err(|e| remap(e, OP))?;
        let (re, im) = if !reaches_zero(&d_mrad) {
            let d2 = ball_again(&d_mrad).map_err(|e| remap(e, OP))?;
            let re = num_re.div(&d_mrad, prec).map_err(|e| remap(e, OP))?;
            let im = num_im.div(&d2, prec).map_err(|e| remap(e, OP))?;
            (re, im)
        } else {
            // Mid-rad sqr can reach 0 even when the rectangle excludes 0 (S1).
            // Divide by the endpoint range of |b|². A mid-rad packing of that
            // range may reach 0 when hi/lo exceeds the Mag mantissa; the
            // quotient range can still be a finite enclosure.
            let wp = match prec.checked_add(64) {
                Some(w) if w <= MAX_PRECISION => w,
                Some(_) | None => prec,
            };
            let Some((lo2, hi2)) = abs2_range(&br, &bi, wp).map_err(|e| remap(e, OP))? else {
                return CBall::all_complex(prec).map_err(|e| remap(e, OP));
            };
            let re = div_by_pos_range(&num_re, &lo2, &hi2, prec).map_err(|e| remap(e, OP))?;
            let im = div_by_pos_range(&num_im, &lo2, &hi2, prec).map_err(|e| remap(e, OP))?;
            (re, im)
        };
        CBall::new(re, im).map_err(|e| remap(e, OP))
    }
}

/// Rebuild a `Ball` from its public mid, rad and prec so a part can be used in two `Ball` ops
/// (`Ball::{add,sub,mul,div,mul_2exp}` still take the first operand by value).
fn exact_from_parts(x: &ExactNum) -> BoundResult<ExactNum> {
    // Midpoints are finite (from_mid_rad). x+0 is a new ExactNum with the same
    // value, so a Ball can be rebuilt without Clone (zero-clone rule).
    Ok(x.add_full_prec(&ExactNum::from(0i64)))
}

fn ball_again(b: &Ball) -> BoundResult<Ball> {
    Ball::from_mid_rad(exact_from_parts(b.mid())?, mag_again(b.rad())?, b.prec())
}

fn mag_again(r: &crate::mag::Mag) -> BoundResult<crate::mag::Mag> {
    if r.is_infinity() {
        return Ok(crate::mag::Mag::infinity());
    }
    if r.is_zero() {
        return Ok(crate::mag::Mag::zero());
    }
    Ok(crate::mag::Mag::from_exact_num(&r.to_exact_num()?))
}

/// `|m| ≤ r`, decided by `Mag::cmp_exact_num`, including an infinite radius (16.CB.7 R1/R2).
fn reaches_zero(b: &Ball) -> bool {
    match b.rad().cmp_exact_num(b.mid()) {
        Some(core::cmp::Ordering::Less) => false,
        Some(core::cmp::Ordering::Equal | core::cmp::Ordering::Greater) | None => true,
    }
}

fn scale_exp(b: &Ball) -> Option<i64> {
    match b.mid().exponent() {
        Some(e) if !b.mid().is_zero() => Some(e),
        _ => None,
    }
}

/// Larger exponent of a nonzero part midpoint. A zero mid no longer contributes
/// exponent 0, which used to block scaling of a tiny other part (e.g. 2^{-40}+0i).
fn scale_k(re: &Ball, im: &Ball) -> i64 {
    match (scale_exp(re), scale_exp(im)) {
        (Some(a), Some(b)) => a.max(b),
        (Some(a), None) => a,
        (None, Some(b)) => b,
        (None, None) => 0,
    }
}

/// Certified range of `|z|²` over the rectangle `re × im`. `None` when a
/// positive lower bound cannot be certified.
fn abs2_range(re: &Ball, im: &Ball, prec: Precision) -> BoundResult<Option<(ExactNum, ExactNum)>> {
    let lo_r = dist_lo(re, prec)?;
    let lo_i = dist_lo(im, prec)?;
    let lo2 = lo_r.sqr(prec, RoundingMode::Down).add(
        &lo_i.sqr(prec, RoundingMode::Down),
        prec,
        RoundingMode::Down,
    );
    if !lo2.is_positive() || lo2.is_nan() {
        return Ok(None);
    }
    let (Some(hi_r), Some(hi_i)) = (dist_hi(re, prec)?, dist_hi(im, prec)?) else {
        return Ok(None);
    };
    let hi2 = hi_r.sqr(prec, RoundingMode::Up).add(
        &hi_i.sqr(prec, RoundingMode::Up),
        prec,
        RoundingMode::Up,
    );
    if hi2.is_nan() || hi2.is_inf() || !hi2.is_positive() {
        return Ok(None);
    }
    Ok(Some((lo2, hi2)))
}

fn dist_lo(b: &Ball, prec: Precision) -> BoundResult<ExactNum> {
    if reaches_zero(b) {
        return Ok(ExactNum::from(0i64));
    }
    let rad = match b.rad().to_exact_num() {
        Ok(r) => r,
        Err(_) => return Ok(ExactNum::from(0i64)),
    };
    Ok(exact_from_parts(b.mid())?
        .abs()
        .sub(&rad, prec, RoundingMode::Down))
}

fn dist_hi(b: &Ball, prec: Precision) -> BoundResult<Option<ExactNum>> {
    if b.rad().is_infinity() {
        return Ok(None);
    }
    let rad = match b.rad().to_exact_num() {
        Ok(r) => r,
        Err(_) => return Ok(None),
    };
    let hi = exact_from_parts(b.mid())?
        .abs()
        .add(&rad, prec, RoundingMode::Up);
    if hi.is_nan() || hi.is_inf() {
        return Ok(None);
    }
    Ok(Some(hi))
}

fn end_lo(b: &Ball, prec: Precision) -> BoundResult<Option<ExactNum>> {
    if b.rad().is_infinity() {
        return Ok(None);
    }
    let rad = match b.rad().to_exact_num() {
        Ok(r) => r,
        Err(_) => return Ok(None),
    };
    let lo = exact_from_parts(b.mid())?.sub(&rad, prec, RoundingMode::Down);
    if lo.is_nan() || lo.is_inf() {
        return Ok(None);
    }
    Ok(Some(lo))
}

fn end_hi(b: &Ball, prec: Precision) -> BoundResult<Option<ExactNum>> {
    if b.rad().is_infinity() {
        return Ok(None);
    }
    let rad = match b.rad().to_exact_num() {
        Ok(r) => r,
        Err(_) => return Ok(None),
    };
    let hi = exact_from_parts(b.mid())?.add(&rad, prec, RoundingMode::Up);
    if hi.is_nan() || hi.is_inf() {
        return Ok(None);
    }
    Ok(Some(hi))
}

/// `num / [dlo, dhi]` with `dlo > 0`, packed as a dest-prec mid-rad ball.
/// The packed ball may reach 0 when the range is very wide; it is still finite
/// unless an endpoint overflows.
fn div_by_pos_range(
    num: &Ball,
    dlo: &ExactNum,
    dhi: &ExactNum,
    prec: Precision,
) -> BoundResult<Ball> {
    let (Some(nlo), Some(nhi)) = (end_lo(num, prec)?, end_hi(num, prec)?) else {
        return Ball::all_reals(prec);
    };
    let q_hi_den_lo = exact_from_parts(&nhi)?.div(dlo, prec, RoundingMode::Up);
    let q_lo_den_lo = exact_from_parts(&nlo)?.div(dlo, prec, RoundingMode::Down);
    let q_hi_den_hi = exact_from_parts(&nhi)?.div(dhi, prec, RoundingMode::Up);
    let q_lo_den_hi = exact_from_parts(&nlo)?.div(dhi, prec, RoundingMode::Down);
    let (qlo, qhi) = if !nlo.is_negative() {
        (q_lo_den_hi, q_hi_den_lo)
    } else if !nhi.is_positive() {
        (q_lo_den_lo, q_hi_den_hi)
    } else {
        (q_lo_den_lo, q_hi_den_lo)
    };
    if qlo.is_nan() || qhi.is_nan() || qlo.is_inf() || qhi.is_inf() {
        return Ball::all_reals(prec);
    }
    pack_interval(qlo, qhi, prec)
}

fn pack_interval(lo: ExactNum, hi: ExactNum, prec: Precision) -> BoundResult<Ball> {
    let mid = lo.add_full_prec(&hi).ldexp(-1, prec, RoundingMode::ToEven);
    if mid.is_nan() || mid.is_inf() {
        return Ball::all_reals(prec);
    }
    let rad = crate::mag::Mag::from_exact_num(&hi.sub_full_prec(&mid))
        .max(crate::mag::Mag::from_exact_num(&mid.sub_full_prec(&lo)));
    Ball::from_mid_rad(mid, rad, prec)
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
