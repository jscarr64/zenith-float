//! Slice S1: Ball core arithmetic: add, sub, neg, mul, mul by an exact scalar, sqr, powi, sqrt, div, inv and hypot on `Ball`, with the
//! radius rules of CONTRACT clauses 2, 3, 9 and 10 and the slice clauses in CONTRACT section 16 "Slice S1".
//!
//! Owner: slice S1 (Spence). Only that slice edits this file.
//! Every midpoint is the zenith-float operation at the target `prec` rounded `ToEven` (O19); the radius is the propagated
//! radius plus the midpoint rounding error ε (9.4, 16.S1.2), all in `Mag` rounding upward. A midpoint overflow gives the
//! all-reals ball (10.4, O2). div, inv and hypot use zenith-float R13 (16.S1.8).

use core::cmp::Ordering;
use core::ops::Neg;

use zenith_float_num::{Error, ExactNum, RoundingMode, EXPONENT_MAX, EXPONENT_MIN};

use crate::ball::Ball;
use crate::error::{BoundError, BoundResult};
use crate::mag::Mag;
use crate::policy::DomainPolicy;
use crate::precision::{Precision, MAX_PRECISION, MIN_PRECISION};
use crate::tri::Tri;

/// Every Ball midpoint rounds to nearest, ties to even (CONTRACT 15.1, O19).
const MID_RM: RoundingMode = RoundingMode::ToEven;
/// Extra bits `powi` carries through its squaring chain before the one final rounding to `prec` (16.S1.6).
const POWI_GUARD_BITS: Precision = 64;
/// Guard bits for the hypot endpoint magnitudes `|m| ± r`.
const HYPOT_GUARD_BITS: Precision = 64;
/// `Mag` mantissa of one half: `0.1b * 2^0`.
const HALF_MANT: u64 = 0x8000_0000_0000_0000;

/// A classified zenith-float midpoint result (CONTRACT 10.4).
enum Mid {
    Value(ExactNum),
    Overflow,
}

impl Neg for Ball {
    type Output = Ball;

    /// `-(m ± r) = (-m) ± r`: exact and infallible; a −0 midpoint becomes +0 (CONTRACT 6.3, 15.1, 16.S1.3).
    fn neg(self) -> Ball {
        self.negated()
    }
}

impl Ball {
    /// `(m ± r) + (n ± s) = round(m + n) ± (r + s + ε)` (CONTRACT 9.3, 16.S1.1).
    ///
    /// # Errors
    /// `InvalidArgument` when `prec` is outside `MIN_PRECISION..=MAX_PRECISION`. A midpoint overflow is `Ok(all reals)`.
    pub fn add(self, other: &Ball, prec: Precision) -> BoundResult<Ball> {
        check_prec("Ball::add", prec)?;
        let (m, r, _) = self.into_parts();
        let rad = r.add_up(other.rad());
        finish(
            "Ball::add",
            m.add(other.mid(), prec, MID_RM),
            &rad,
            prec,
            false,
        )
    }

    /// `(m ± r) - (n ± s) = round(m - n) ± (r + s + ε)` (CONTRACT 9.3, 16.S1.1).
    ///
    /// # Errors
    /// `InvalidArgument` when `prec` is outside `MIN_PRECISION..=MAX_PRECISION`. A midpoint overflow is `Ok(all reals)`.
    pub fn sub(self, other: &Ball, prec: Precision) -> BoundResult<Ball> {
        check_prec("Ball::sub", prec)?;
        let (m, r, _) = self.into_parts();
        let rad = r.add_up(other.rad());
        finish(
            "Ball::sub",
            m.sub(other.mid(), prec, MID_RM),
            &rad,
            prec,
            false,
        )
    }

    /// `(m ± r)(n ± s) = round(m n) ± (|m| s + |n| r + r s + ε)` (CONTRACT 9.3, 16.S1.1).
    ///
    /// # Errors
    /// `InvalidArgument` when `prec` is outside `MIN_PRECISION..=MAX_PRECISION`. A midpoint overflow is `Ok(all reals)`.
    pub fn mul(self, other: &Ball, prec: Precision) -> BoundResult<Ball> {
        check_prec("Ball::mul", prec)?;
        let (m, r, _) = self.into_parts();
        let (n, s) = (other.mid(), other.rad());
        let rad = Mag::from_exact_num(&m)
            .mul_up(s)
            .add_up(&Mag::from_exact_num(n).mul_up(&r))
            .add_up(&r.mul_up(s));
        finish("Ball::mul", m.mul(n, prec, MID_RM), &rad, prec, false)
    }

    /// `(m ± r) · c = round(m c) ± (|c| r + ε)` for an exact scalar `c` (16.S1.4).
    ///
    /// # Errors
    /// `InvalidArgument` for a bad `prec`; `Representation` when `c` is NaN or infinite. A midpoint overflow is
    /// `Ok(all reals)`.
    pub fn mul_exact(self, c: &ExactNum, prec: Precision) -> BoundResult<Ball> {
        check_prec("Ball::mul_exact", prec)?;
        if c.is_nan() || c.is_inf() {
            return Err(BoundError::Representation {
                op: "Ball::mul_exact",
                detail: "NaN or infinite scalar",
            });
        }
        let (m, r, _) = self.into_parts();
        let rad = Mag::from_exact_num(c).mul_up(&r);
        finish("Ball::mul_exact", m.mul(c, prec, MID_RM), &rad, prec, false)
    }

    /// `(m ± r)^2 = round(m^2) ± (2|m| r + r^2 + ε)` (16.S1.5): tighter than `mul` with itself, which would also
    /// bound the same set but needs a second copy of the ball.
    ///
    /// # Errors
    /// `InvalidArgument` when `prec` is outside `MIN_PRECISION..=MAX_PRECISION`. A midpoint overflow is `Ok(all reals)`.
    pub fn sqr(&self, prec: Precision) -> BoundResult<Ball> {
        check_prec("Ball::sqr", prec)?;
        let (m, r) = (self.mid(), self.rad());
        let t = Mag::from_exact_num(m).mul_up(r);
        let rad = t.add_up(&t).add_up(&r.mul_up(r));
        // The exact square, then one rounding to `prec`. `inexact()` of a product is its own per-operation flag, so it
        // reports only a real dropped bit inside the square, which forces ε (16.S1.5). Since zf 0cfee92 a square below the
        // range keeps every bit (up to MAX_PRECISION), so `lost` is false there and ε comes from round_to's own flag.
        // A `mul_full_prec` exponent overflow passes through `finish` to all reals; pinned by
        // tests/harness/slice_arith.rs `overflow_is_all_reals` (line 424).
        let sq = m.mul_full_prec(m);
        let lost = sq.inexact();
        finish("Ball::sqr", round_to(&sq, prec), &rad, prec, lost)
    }

    /// `(m ± r)^n` by left-to-right binary exponentiation with this module's `sqr` and `mul` at `prec + 64` bits and
    /// one final rounding to `prec`; each step is a rigorous enclosure, so the chain is (16.S1.6). `n == 0` is exactly
    /// `1 ± 0` for every ball, including one containing 0 (CONTRACT 10.2e).
    ///
    /// # Errors
    /// `InvalidArgument` when `prec` is outside `MIN_PRECISION..=MAX_PRECISION`. A midpoint overflow is `Ok(all reals)`.
    pub fn powi(&self, n: u64, prec: Precision) -> BoundResult<Ball> {
        check_prec("Ball::powi", prec)?;
        if n == 0 {
            return Ball::exact(ExactNum::from(1i64), prec);
        }
        let wp = match prec.checked_add(POWI_GUARD_BITS) {
            Some(w) if w <= MAX_PRECISION => w,
            Some(_) | None => prec,
        };
        let mut acc = Ball::exact(ExactNum::from(1i64), wp)?;
        let mut started = false;
        for i in (0..u64::BITS).rev() {
            let bit = n.checked_shr(i).map_or(0, |v| v & 1) == 1;
            if started {
                acc = acc.sqr(wp)?;
            }
            if bit {
                acc = acc.mul(self, wp)?;
                started = true;
            }
            if !acc.is_informative() {
                return Ball::all_reals(prec);
            }
        }
        let (m, r, _) = acc.into_parts();
        finish("Ball::powi", round_to(&m, prec), &r, prec, false)
    }

    /// `sqrt(m ± r)` on the closed domain `[0, inf)` (CONTRACT 10.2b, D15d, 16.S1.7). With `lo = m - r` and
    /// `hi = m + r`, decided exactly: `lo > 0` gives `round(sqrt m) ± (r / sqrt(m) + ε)`; `lo == 0` gives the endpoint ball
    /// `[0, sqrt_up(hi)]`; `lo < 0` is `Domain` under `Strict`. Under `Clip`, `hi < 0` is `Domain`, `hi == 0` is exactly
    /// `0`, and `hi > 0` gives the endpoint ball. An all-reals input is `Domain` (`Strict`) or all reals (`Clip`, 10.2a).
    ///
    /// # Errors
    /// `InvalidArgument` for a bad `prec`; `Domain` as above; `Resource` if the sign of `m - r` cannot be decided.
    pub fn sqrt(&self, prec: Precision, policy: &DomainPolicy) -> BoundResult<Ball> {
        const OP: &str = "Ball::sqrt";
        check_prec(OP, prec)?;
        let strict = match policy {
            DomainPolicy::Strict => true,
            DomainPolicy::Clip => false,
        };
        let (m, r) = (self.mid(), self.rad());
        if r.is_infinity() {
            return if strict {
                Err(domain(OP, "the all-reals ball reaches below 0"))
            } else {
                Ball::all_reals(prec)
            };
        }
        let Some(ord) = r.cmp_exact_num(m).map(Ordering::reverse) else {
            return Err(BoundError::Resource {
                op: OP,
                detail: "cannot decide the sign of m - r",
            });
        };
        let neg = m.is_negative() && !m.is_zero();
        let pos = !m.is_negative() && !m.is_zero();
        // (sign of lo, sign of hi) as Ordering against 0.
        let (lo, hi) = match (neg, pos, ord) {
            (false, true, Ordering::Greater) => (Ordering::Greater, Ordering::Greater),
            (false, true, Ordering::Equal) => (Ordering::Equal, Ordering::Greater),
            (true, false, Ordering::Greater) => (Ordering::Less, Ordering::Less),
            (true, false, Ordering::Equal) => (Ordering::Less, Ordering::Equal),
            (false, false, Ordering::Equal | Ordering::Greater) => {
                (Ordering::Equal, Ordering::Equal)
            }
            (_, _, Ordering::Less) | (true, true, _) => (Ordering::Less, Ordering::Greater),
        };
        match (lo, strict, hi) {
            (Ordering::Greater, _, _) => {
                let root_lower = Mag::from_exact_num_lower(m).sqrt_lower();
                let rad = r.div_up(&root_lower);
                finish(OP, m.sqrt(prec, MID_RM), &rad, prec, false)
            }
            (Ordering::Equal, _, _) | (Ordering::Less, false, Ordering::Greater) => {
                let root_up = if neg {
                    // m < 0 < m + r: hi = m + r exactly when r is a zenith-float value (16.S1.7), else r (>= m + r).
                    match r.to_exact_num() {
                        Ok(rn) => Mag::from_exact_num(&rn.add_full_prec(m)).sqrt_up(),
                        Err(_) => r.sqrt_up(),
                    }
                } else {
                    Mag::from_exact_num(m).add_up(r).sqrt_up()
                };
                endpoint_ball(OP, root_up, prec)
            }
            (Ordering::Less, true, _) => Err(domain(OP, "m - r < 0 under Strict")),
            (Ordering::Less, false, Ordering::Less) => {
                Err(domain(OP, "m + r < 0: the ball misses [0, inf)"))
            }
            (Ordering::Less, false, Ordering::Equal) => Ball::exact(ExactNum::from(0i64), prec),
        }
    }
    /// `(m ± r) / (n ± s) = round(m / n) ± ((|m| s + |n| r) / (|n| (|n| - s)) + ε)` when `|n| > s` (16.S1.8). The
    /// midpoint is `div_exact` when the quotient is exact, else zenith-float `div` rounded `ToEven`.
    ///
    /// # Errors
    /// `InvalidArgument` for a bad `prec`; `Domain` when the divisor is exactly `0 ± 0`. A divisor ball reaching 0
    /// (`|n| <= s`), an infinite radius, or a midpoint overflow is `Ok(all reals)`.
    pub fn div(self, other: &Ball, prec: Precision) -> BoundResult<Ball> {
        const OP: &str = "Ball::div";
        check_prec(OP, prec)?;
        let (n, s) = (other.mid(), other.rad());
        if n.is_zero() && s.is_zero() {
            return Err(domain(OP, "division by the exact zero ball"));
        }
        let Some(n_lo) = clear_of_zero(n, s) else {
            return Ball::all_reals(prec);
        };
        let (m, r, _) = self.into_parts();
        let num = Mag::from_exact_num(&m)
            .mul_up(s)
            .add_up(&Mag::from_exact_num(n).mul_up(&r));
        let rad = num.div_up(&n_lo.mul_lower(&n_lo.sub_lower(s)));
        let q = match m.div_exact(n, prec) {
            Ok(Some(q)) => q,
            Ok(None) | Err(_) => m.div(n, prec, MID_RM),
        };
        finish(OP, q, &rad, prec, false)
    }

    /// `1 / (n ± s) = round(1 / n) ± (s / (|n| (|n| - s)) + ε)` when `|n| > s` (16.S1.8).
    ///
    /// # Errors
    /// `InvalidArgument` for a bad `prec`; `Domain` for the exact zero ball. A ball reaching 0, an infinite radius, or a
    /// midpoint overflow is `Ok(all reals)`.
    pub fn inv(&self, prec: Precision) -> BoundResult<Ball> {
        const OP: &str = "Ball::inv";
        check_prec(OP, prec)?;
        let (n, s) = (self.mid(), self.rad());
        if n.is_zero() && s.is_zero() {
            return Err(domain(OP, "inverse of the exact zero ball"));
        }
        let Some(n_lo) = clear_of_zero(n, s) else {
            return Ball::all_reals(prec);
        };
        let rad = s.div_up(&n_lo.mul_lower(&n_lo.sub_lower(s)));
        finish(OP, n.reciprocal(prec, MID_RM), &rad, prec, false)
    }

    /// Monotone enclosure (16.S1.8): hypot is increasing in `|x|` and `|y|`, so the result lies in
    /// `[hypot↓(min|X|, min|Y|), hypot↑(max|X|, max|Y|)]` (endpoints at `prec`), returned as the ball
    /// `round((lo + hi) / 2) ± max(hi - mid, mid - lo)` rounded out. A radius outside zf's exponent range falls back to
    /// the 1-Lipschitz `round(hypot(m, n)) ± (r + s + ε)`.
    ///
    /// # Errors
    /// `InvalidArgument` for a bad `prec`. An infinite radius or an overflowing upper endpoint is `Ok(all reals)`.
    pub fn hypot(&self, other: &Ball, prec: Precision) -> BoundResult<Ball> {
        const OP: &str = "Ball::hypot";
        check_prec(OP, prec)?;
        let (m, r) = (self.mid(), self.rad());
        let (n, s) = (other.mid(), other.rad());
        let (Ok(rn), Ok(sn)) = (r.to_exact_num(), s.to_exact_num()) else {
            let rad = r.add_up(s);
            return finish(OP, m.hypot(n, prec, MID_RM), &rad, prec, false);
        };
        let wp = match prec.checked_add(HYPOT_GUARD_BITS) {
            Some(w) if w <= MAX_PRECISION => w,
            Some(_) | None => prec,
        };
        let lo =
            abs_lower(m, r, &rn, wp).hypot(&abs_lower(n, s, &sn, wp), prec, RoundingMode::Down);
        let hi = abs_upper(m, &rn, wp).hypot(&abs_upper(n, &sn, wp), prec, RoundingMode::Up);
        let (Mid::Value(lo), Mid::Value(hi)) = (classify(OP, lo)?, classify(OP, hi)?) else {
            return Ball::all_reals(prec);
        };
        let Mid::Value(mid) = classify(OP, lo.add_full_prec(&hi).ldexp(-1, prec, MID_RM))? else {
            return Ball::all_reals(prec);
        };
        let rad = Mag::from_exact_num(&hi.sub_full_prec(&mid))
            .max(Mag::from_exact_num(&mid.sub_full_prec(&lo)));
        Ball::from_mid_rad(mid, rad, prec)
    }
}

impl Ball {
    /// `|m ± r|` (CONTRACT 15.1, O14-W, 16.S1.10): `round(|m|) ± (r + ε)` when the ball does not straddle 0 (`|m| ≥ r`,
    /// decided exactly); otherwise the endpoint ball `[0, |m| + r]` with the `Mag` sum rounded up (`add_up`).
    ///
    /// # Errors
    /// `InvalidArgument` when `prec` is outside `MIN_PRECISION..=MAX_PRECISION`.
    pub fn abs(&self, prec: Precision) -> BoundResult<Ball> {
        const OP: &str = "Ball::abs";
        check_prec(OP, prec)?;
        let (m, r) = (self.mid(), self.rad());
        match r.cmp_exact_num(m).map(Ordering::reverse) {
            Some(Ordering::Greater | Ordering::Equal) => {
                // round(|m|) as one product (-1 or 1)·m, so `m` is borrowed, not consumed (O14-W).
                let sign = if m.is_negative() { -1i64 } else { 1i64 };
                let am = ExactNum::from(sign).mul(m, prec, MID_RM);
                finish(OP, am, r, prec, false)
            }
            Some(Ordering::Less) | None => {
                let u = Mag::from_exact_num(m).add_up(r);
                endpoint_ball(OP, u, prec)
            }
        }
    }

    /// `(m ± r) · 2^k = round(m · 2^k) ± (r · 2^k + ε)` (CONTRACT 15.1, O14-W, 16.S1.10). The midpoint is exact unless
    /// it over- or underflows or has more than `prec` bits; the radius exponent is checked (overflow is `Infinity`,
    /// underflow the smallest positive `Mag`).
    ///
    /// # Errors
    /// `InvalidArgument` when `prec` is outside `MIN_PRECISION..=MAX_PRECISION`. A midpoint overflow is `Ok(all reals)`.
    pub fn mul_2exp(self, k: i64, prec: Precision) -> BoundResult<Ball> {
        const OP: &str = "Ball::mul_2exp";
        check_prec(OP, prec)?;
        let (m, r, _) = self.into_parts();
        finish(
            OP,
            m.ldexp(k, prec, MID_RM),
            &mag_mul_2exp_up(&r, k),
            prec,
            false,
        )
    }

    /// `self < other` decided by exact endpoint comparison (CONTRACT 15.2, D14, 16.S1.11).
    pub fn lt(&self, other: &Ball) -> Tri {
        match endpoint_signs(self, other) {
            Some((Ordering::Less, _)) => Tri::True,
            // lo(a) >= hi(b): every point of `a` is >= every point of `b`, so `a < b` is certainly False, also when the
            // balls only touch (R > 0, lo(a) = hi(b)).
            Some((_, Ordering::Greater | Ordering::Equal)) => Tri::False,
            Some((Ordering::Equal | Ordering::Greater, Ordering::Less)) | None => Tri::Unknown,
        }
    }

    /// `self ≤ other` decided by exact endpoint comparison (CONTRACT 15.2, D14, 16.S1.11).
    pub fn le(&self, other: &Ball) -> Tri {
        match endpoint_signs(self, other) {
            Some((Ordering::Less | Ordering::Equal, _)) => Tri::True,
            Some((Ordering::Greater, Ordering::Greater)) => Tri::False,
            Some((Ordering::Greater, Ordering::Equal | Ordering::Less)) | None => Tri::Unknown,
        }
    }

    /// `self > other`, i.e. `other < self` (CONTRACT 15.2, 16.S1.11).
    pub fn gt(&self, other: &Ball) -> Tri {
        other.lt(self)
    }

    /// `self ≥ other`, i.e. `other ≤ self` (CONTRACT 15.2, 16.S1.11).
    pub fn ge(&self, other: &Ball) -> Tri {
        other.le(self)
    }

    /// Value equality: `True` only for two exact equal points, `False` for disjoint balls (CONTRACT 15.2, 16.S1.11).
    pub fn eq(&self, other: &Ball) -> Tri {
        match endpoint_signs(self, other) {
            Some((Ordering::Equal, Ordering::Equal)) => Tri::True,
            Some((Ordering::Less, _) | (_, Ordering::Greater)) => Tri::False,
            Some((Ordering::Equal | Ordering::Greater, Ordering::Equal | Ordering::Less))
            | None => Tri::Unknown,
        }
    }

    /// The negation of [`Self::eq`] (CONTRACT 15.2, 16.S1.11).
    pub fn ne(&self, other: &Ball) -> Tri {
        match self.eq(other) {
            Tri::True => Tri::False,
            Tri::False => Tri::True,
            Tri::Unknown => Tri::Unknown,
        }
    }

    /// Whether the two closed balls share a point (CONTRACT 15.2, O20, 16.S1.11); an `Infinity` radius overlaps
    /// everything.
    pub fn overlaps(&self, other: &Ball) -> Tri {
        if self.rad().is_infinity() || other.rad().is_infinity() {
            return Tri::True;
        }
        match endpoint_signs(self, other) {
            Some((Ordering::Less, _) | (_, Ordering::Greater)) => Tri::False,
            Some((Ordering::Equal | Ordering::Greater, Ordering::Equal | Ordering::Less)) => {
                Tri::True
            }
            None => Tri::Unknown,
        }
    }
}

/// The signs of `D + R` and `D − R` (CONTRACT 15.2, 16.S1.11) with `D = m1 − m2` and `R = r1 + r2`, decided exactly by
/// the S2 comparisons (16.S2): `D` is formed with `sub_full_prec`, and `|D|` against `R` is `r1` against `|D|` and, when
/// `|D| > r1`, `r2` against `||D| − r1|`, so `R` is never rounded. `None` (undecided) for an `Infinity` radius, an
/// inexact difference, or `|D| > r1` with `r1` outside zenith-float's range.
fn endpoint_signs(a: &Ball, b: &Ball) -> Option<(Ordering, Ordering)> {
    let (r1, r2) = (a.rad(), b.rad());
    if r1.is_infinity() || r2.is_infinity() {
        return None;
    }
    let d = a.mid().sub_full_prec(b.mid());
    if d.inexact() || d.is_nan() {
        return None;
    }
    // |D| against R.
    let c = match r1.cmp_exact_num(&d)? {
        Ordering::Greater => Ordering::Less,
        Ordering::Equal if r2.is_zero() => Ordering::Equal,
        Ordering::Equal => Ordering::Less,
        Ordering::Less => far_side(r1, r2, &d)?,
    };
    let sd = if d.is_zero() {
        Ordering::Equal
    } else if d.is_negative() {
        Ordering::Less
    } else {
        Ordering::Greater
    };
    let lo = if sd == Ordering::Less {
        Ordering::Less
    } else {
        c
    };
    let hi = if sd == Ordering::Greater {
        Ordering::Greater
    } else {
        c.reverse()
    };
    Some((hi, lo))
}

/// `|D|` against `r1 + r2` once `|D| > r1`: `r2` against `||D| − r1|` (16.S2.2) with `r1` converted, or, when `r1` lies
/// outside zenith-float's range, the same with the roles swapped. `None` when neither radius converts.
fn far_side(r1: &Mag, r2: &Mag, d: &ExactNum) -> Option<Ordering> {
    if r2.is_zero() {
        return Some(Ordering::Greater);
    }
    if let Ok(mut rn) = r1.to_exact_num() {
        if d.is_negative() {
            rn.inv_sign();
        }
        return Some(r2.cmp_exact_diff(d, &rn)?.reverse());
    }
    // r1 > 0 is outside the range, so |D| against r2 alone first.
    match r2.cmp_exact_num(d)? {
        Ordering::Greater | Ordering::Equal => Some(Ordering::Less),
        Ordering::Less => {
            let mut sn = r2.to_exact_num().ok()?;
            if d.is_negative() {
                sn.inv_sign();
            }
            Some(r1.cmp_exact_diff(d, &sn)?.reverse())
        }
    }
}

/// `r · 2^k` rounded up with a checked exponent: exact while the exponent fits, `Infinity` on overflow and the smallest
/// positive `Mag` on underflow (CONTRACT 15.1, 3.2, 3.3).
fn mag_mul_2exp_up(r: &Mag, k: i64) -> Mag {
    if r.is_zero() {
        return Mag::zero();
    }
    let (Some(mant), Some(e)) = (r.mant(), r.exp()) else {
        return Mag::infinity();
    };
    let shifted = match e.checked_add(k) {
        Some(ne) => Mag::finite(mant, ne),
        None if k < 0 => Mag::finite(HALF_MANT, i64::MIN),
        None => return Mag::infinity(),
    };
    // The Err arm is unreachable: `r` is a canonical Mag (from_canonical keeps the mantissa's top bit set), so
    // Mag::finite's normalising shift is lz = 0 and `ne` itself is the stored exponent, which fits i64 by
    // construction; the only out-of-range case (e + k below i64::MIN) is handled above as finite(HALF_MANT, i64::MIN),
    // which has lz = 0 too. Infinity is kept as the outward-safe answer should that ever change.
    match shifted {
        Ok(m) => m,
        Err(_) => Mag::infinity(),
    }
}

/// A lower bound on `min |m ± r|` at `wp` bits: 0 when the ball reaches 0, else `|m| - r` rounded toward 0.
fn abs_lower(m: &ExactNum, r: &Mag, rn: &ExactNum, wp: Precision) -> ExactNum {
    match r.cmp_exact_num(m).map(Ordering::reverse) {
        Some(Ordering::Greater) => {
            let x = ExactNum::from(1i64).mul(m, wp, RoundingMode::ToZero);
            if m.is_negative() {
                x.add(rn, wp, RoundingMode::ToZero)
            } else {
                x.sub(rn, wp, RoundingMode::ToZero)
            }
        }
        Some(Ordering::Equal | Ordering::Less) | None => ExactNum::from(0i64),
    }
}

/// An upper bound on `max |m ± r|` at `wp` bits: `|m| + r` rounded away from 0.
fn abs_upper(m: &ExactNum, rn: &ExactNum, wp: Precision) -> ExactNum {
    let x = ExactNum::from(1i64).mul(m, wp, RoundingMode::FromZero);
    if m.is_negative() {
        x.sub(rn, wp, RoundingMode::FromZero)
    } else {
        x.add(rn, wp, RoundingMode::FromZero)
    }
}

/// A lower bound on `|n|` when the ball `n ± s` excludes 0 (`|n| > s`, decided exactly), else `None`.
fn clear_of_zero(n: &ExactNum, s: &Mag) -> Option<Mag> {
    if s.is_infinity() {
        return None;
    }
    match s.cmp_exact_num(n).map(Ordering::reverse) {
        Some(Ordering::Greater) => Some(Mag::from_exact_num_lower(n)),
        Some(Ordering::Equal | Ordering::Less) | None => None,
    }
}

/// The ball `[0, u]` as `round(u / 2) ± (u / 2 + ε)`; when `u` is outside zenith-float's range, `0 ± u` (16.S1.7).
fn endpoint_ball(op: &'static str, u: Mag, prec: Precision) -> BoundResult<Ball> {
    let half = Mag::finite(HALF_MANT, 0)?;
    let u_half = u.mul_up(&half);
    match (u.to_exact_num(), ExactNum::from_mag_parts(HALF_MANT, 0, 64)) {
        (Ok(un), Ok(h)) => finish(op, un.mul(&h, prec, MID_RM), &u_half, prec, false),
        (Err(_), _) | (_, Err(_)) => Ball::from_mid_rad(ExactNum::from(0i64), u, prec),
    }
}

/// `x` rounded once to `prec`, `ToEven`; `inexact()` reports this rounding only.
fn round_to(x: &ExactNum, prec: Precision) -> ExactNum {
    ExactNum::from(1i64).mul(x, prec, MID_RM)
}

fn domain(op: &'static str, detail: &'static str) -> BoundError {
    BoundError::Domain { op, detail }
}

fn check_prec(op: &'static str, prec: Precision) -> BoundResult<()> {
    if (MIN_PRECISION..=MAX_PRECISION).contains(&prec) {
        Ok(())
    } else {
        Err(BoundError::InvalidArgument {
            op,
            detail: "prec outside MIN_PRECISION..=MAX_PRECISION",
        })
    }
}

/// CONTRACT 10.4 / O2: overflow detection and the NaN error-kind mapping.
fn classify(op: &'static str, raw: ExactNum) -> BoundResult<Mid> {
    if raw.is_inf() {
        return Ok(Mid::Overflow);
    }
    if raw.is_nan() {
        return match raw.err() {
            Some(Error::ExponentOverflow(_)) => Ok(Mid::Overflow),
            Some(Error::DivisionByZero) => Err(domain(op, "zf division by zero")),
            Some(Error::InvalidArgument) => Err(BoundError::InvalidArgument {
                op,
                detail: "zf invalid argument",
            }),
            Some(Error::MemoryAllocation) => Err(BoundError::Resource {
                op,
                detail: "zf memory allocation",
            }),
            Some(Error::PrecisionRetryExhausted) => Err(BoundError::Resource {
                op,
                detail: "zf precision retry exhausted",
            }),
            None => Err(BoundError::Representation {
                op,
                detail: "zf NaN without an error kind",
            }),
        };
    }
    // Deliberate: an inexact result in the top binade (exponent EXPONENT_MAX) may have been clamped there by zf, so
    // its true value can lie beyond any finite enclosure; it is reported as overflow, i.e. all reals (CONTRACT 15.1).
    if raw.inexact() && raw.exponent() == Some(EXPONENT_MAX) {
        return Ok(Mid::Overflow);
    }
    Ok(Mid::Value(raw))
}

/// ε of a rounded midpoint (CONTRACT 9.4, O19, 16.S1.2): zero when exact; `half_ulp(e, prec)` for a result with exponent
/// `e > EXPONENT_MIN`; and `2^EXPONENT_MIN` for a zero or bottom-binade (subnormal) result, which bounds any nearest
/// rounding there (both the exact value and the result lie in `[-2^EXPONENT_MIN, 2^EXPONENT_MIN]` with one sign).
fn eps(mid: &ExactNum, prec: Precision, lost: bool) -> Mag {
    if !mid.inexact() && !lost {
        return Mag::zero();
    }
    match mid.exponent() {
        Some(e) if e > EXPONENT_MIN && !mid.is_zero() => Mag::half_ulp(e, prec),
        Some(_) | None => Mag::ulp_mag(EXPONENT_MIN, 0),
    }
}

/// Classify the midpoint, add ε to the propagated radius and build the ball.
fn finish(
    op: &'static str,
    raw: ExactNum,
    prop: &Mag,
    prec: Precision,
    lost: bool,
) -> BoundResult<Ball> {
    match classify(op, raw)? {
        Mid::Overflow => Ball::all_reals(prec),
        Mid::Value(m) => {
            let rad = prop.add_up(&eps(&m, prec, lost));
            Ball::from_mid_rad(m, rad, prec)
        }
    }
}
