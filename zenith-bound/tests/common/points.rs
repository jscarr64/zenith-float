//! L0b P1 sample points on exact dyadics (TEST_PLAN 7.1 and S3, Cassie's grid `[Cas-7]`): for a ball described by `generators::BallDesc`
//! the points are `mid`, `mid ± r`, `mid ± r*j/8` (`j = 1..7`) and one fixed-seed point `mid + r*j/2^t` (`j` drawn in `-2^t..=2^t`), each
//! kept EXACT as an unevaluated sum `base + off` of two signed dyadic terms (no rounding, no zenith-float, no Ball type). `cmp_points`
//! compares two points exactly with a small limb integer, refusing (Err) when the exponent spread exceeds `max_shift` bits.
//! Value conventions: a `BallDesc` mid is `(-1)^sign * 0.mid_mant * 2^mid_exp` and a `RandomMag` radius `0.rad_mant * 2^rad_exp` (64-bit
//! mantissas); a `Term` is `(-1)^neg * mant * 2^exp`, normalised (`mant` odd, or the zero term `+0 * 2^0`), so both zeros are the same point.
//! `OneUlp`/`HalfUlp` resolve at the given precision to `2^(mid_exp - prec)` / `2^(mid_exp - prec - 1)`; `Huge` has no decided value and
//! `Infinity` (the entire line) has no finite sample grid: both are an Err. Integer-only; every failure is an Err, nothing panics.

use core::cmp::Ordering;

use super::generators::{BallDesc, RadiusClass, Sign};
use super::modular::sub_u64;
use super::rng::SplitMix64;

/// `(-1)^neg * mant * 2^exp`, normalised. The exponent is `i128` so every grid of an `i64` mid exponent at a `u32` precision is
/// exact (headroom for zenith-float 2.0 `i64` exponents: `mid_exp - prec - 1 - t` stays far inside `i128`).
#[derive(Debug, PartialEq, Eq)]
pub struct Term {
    pub neg: bool,
    pub mant: u128,
    pub exp: i128,
}

impl Term {
    /// The normalised term for `(-1)^neg * mant * 2^exp`.
    pub fn new(neg: bool, mant: u128, exp: i128) -> Result<Term, String> {
        if mant == 0 {
            return Ok(Term {
                neg: false,
                mant: 0,
                exp: 0,
            });
        }
        let tz = mant.trailing_zeros();
        let e = exp
            .checked_add(i128::from(tz))
            .ok_or_else(|| format!("term exponent overflow: {} + {}", exp, tz))?;
        Ok(Term {
            neg,
            mant: mant >> tz,
            exp: e,
        })
    }

    pub fn is_zero(&self) -> bool {
        self.mant == 0
    }
}

/// An exact point `base + off`.
#[derive(Debug, PartialEq, Eq)]
pub struct Point {
    pub base: Term,
    pub off: Term,
}

/// Which grid point a sample is.
#[derive(Debug, PartialEq, Eq)]
pub enum SampleKind {
    Mid,
    /// `mid + r` (`plus`) or `mid - r`.
    Edge {
        plus: bool,
    },
    /// `mid + r*j/8`, `j` in `-7..=-1` and `1..=7`.
    Eighth {
        j: i64,
    },
    /// `mid + r*j/2^t`.
    Seeded {
        j: i64,
        t: u32,
    },
}

/// Largest `t` for the seeded point (`|j| <= 2^t` times a 64-bit radius mantissa stays below 2^128).
pub const MAX_T: u32 = 62;
/// Default exponent spread accepted by `cmp_points` (bits).
pub const DEFAULT_MAX_SHIFT: u64 = 1 << 16;

fn mid_term(b: &BallDesc) -> Result<Term, String> {
    if b.mid_is_zero {
        return Term::new(false, 0, 0);
    }
    if b.mid_mant == 0 {
        return Err(String::from("nonzero mid with a zero mantissa"));
    }
    let e = i128::from(b.mid_exp)
        .checked_sub(64)
        .ok_or_else(|| String::from("mid exponent underflow"))?;
    Term::new(b.sign == Sign::Minus, u128::from(b.mid_mant), e)
}

/// The radius as `(mant, exp)` with value `mant * 2^exp`; `Ok(None)` for a zero radius.
pub fn radius_value(b: &BallDesc, prec: u32) -> Result<Option<(u64, i128)>, String> {
    match b.radius {
        RadiusClass::Zero => Ok(None),
        RadiusClass::RandomMag => {
            if b.rad_mant == 0 {
                return Ok(None);
            }
            let e = i128::from(b.rad_exp)
                .checked_sub(64)
                .ok_or_else(|| String::from("radius exponent underflow"))?;
            Ok(Some((b.rad_mant, e)))
        }
        RadiusClass::OneUlp | RadiusClass::HalfUlp => {
            if b.mid_is_zero {
                return Err(String::from("an ulp radius needs a nonzero mid"));
            }
            let half = i128::from(b.radius == RadiusClass::HalfUlp);
            let e = i128::from(b.mid_exp)
                .checked_sub(i128::from(prec))
                .and_then(|v| v.checked_sub(half))
                .ok_or_else(|| String::from("ulp exponent underflow"))?;
            Ok(Some((1, e)))
        }
        RadiusClass::Huge => Err(String::from("Huge radius has no decided value yet")),
        RadiusClass::Infinity => Err(String::from(
            "Infinity radius (entire line) has no finite sample grid",
        )),
    }
}

fn offset(rad: Option<(u64, i128)>, j: i64, t: u32) -> Result<Term, String> {
    match rad {
        None => Term::new(false, 0, 0),
        Some((m, e)) => {
            let prod = u128::from(m)
                .checked_mul(u128::from(j.unsigned_abs()))
                .ok_or_else(|| String::from("offset mantissa overflow"))?;
            let ex = e
                .checked_sub(i128::from(t))
                .ok_or_else(|| String::from("offset exponent underflow"))?;
            Term::new(j < 0, prod, ex)
        }
    }
}

/// The P1 grid for `b` at precision `prec`: 18 points in a fixed order (Mid, Edge +/-, Eighth -7..=-1, 1..=7, Seeded), the seeded `j` from `rng`.
pub fn sample_points(
    b: &BallDesc,
    prec: u32,
    rng: &mut SplitMix64,
    t: u32,
) -> Result<Vec<(SampleKind, Point)>, String> {
    if prec == 0 {
        return Err(String::from("precision must be >= 1"));
    }
    if t > MAX_T {
        return Err(format!("t = {} is above {}", t, MAX_T));
    }
    let rad = radius_value(b, prec)?;
    let mut out = Vec::with_capacity(18);
    out.push((
        SampleKind::Mid,
        Point {
            base: mid_term(b)?,
            off: offset(rad, 0, 0)?,
        },
    ));
    for plus in [true, false] {
        let j = if plus { 1 } else { -1 };
        out.push((
            SampleKind::Edge { plus },
            Point {
                base: mid_term(b)?,
                off: offset(rad, j, 0)?,
            },
        ));
    }
    for j in (-7i64..=-1).chain(1..=7) {
        out.push((
            SampleKind::Eighth { j },
            Point {
                base: mid_term(b)?,
                off: offset(rad, j, 3)?,
            },
        ));
    }
    let span = 1u64 << t;
    let draw = rng.below(
        span.checked_mul(2)
            .and_then(|v| v.checked_add(1))
            .ok_or_else(|| String::from("span"))?,
    );
    let j = i64::try_from(draw)
        .map_err(|_| String::from("draw"))?
        .checked_sub(i64::try_from(span).map_err(|_| String::from("span"))?)
        .ok_or_else(|| String::from("seeded j"))?;
    out.push((
        SampleKind::Seeded { j, t },
        Point {
            base: mid_term(b)?,
            off: offset(rad, j, t)?,
        },
    ));
    Ok(out)
}

// ---- exact comparison with little-endian u64 limbs ----

fn shifted(mant: u128, shift: u64) -> Result<Vec<u64>, String> {
    let limbs = usize::try_from(shift / 64).map_err(|_| String::from("shift"))?;
    let bits = u32::try_from(shift % 64).map_err(|_| String::from("shift"))?;
    let mut v: Vec<u64> = vec![0; limbs];
    let lo = u64::try_from(mant & u128::from(u64::MAX)).map_err(|_| String::from("limb"))?;
    let hi = u64::try_from(mant >> 64).map_err(|_| String::from("limb"))?;
    let mut carry = 0u64;
    for w in [lo, hi] {
        if bits == 0 {
            v.push(w);
        } else {
            v.push((w << bits) | carry);
            carry = w >> (64 - bits);
        }
    }
    v.push(carry);
    Ok(v)
}

fn mag_cmp(a: &[u64], b: &[u64]) -> Ordering {
    let n = a.len().max(b.len());
    for i in (0..n).rev() {
        let x = a.get(i).map_or(0, |v| *v);
        let y = b.get(i).map_or(0, |v| *v);
        match x.cmp(&y) {
            Ordering::Equal => continue,
            o @ (Ordering::Less | Ordering::Greater) => return o,
        }
    }
    Ordering::Equal
}

fn mag_add(a: &[u64], b: &[u64]) -> Vec<u64> {
    let n = a.len().max(b.len());
    let mut out = Vec::with_capacity(n + 1);
    let mut carry = 0u128;
    for i in 0..n {
        let s =
            u128::from(a.get(i).map_or(0, |v| *v)) + u128::from(b.get(i).map_or(0, |v| *v)) + carry;
        out.push(u64::try_from(s & u128::from(u64::MAX)).unwrap_or(0));
        carry = s >> 64;
    }
    out.push(u64::try_from(carry).unwrap_or(0));
    out
}

/// `a - b` for `a >= b`.
fn mag_sub(a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut out = Vec::with_capacity(a.len());
    let mut borrow = false;
    for i in 0..a.len() {
        let x = a.get(i).map_or(0, |v| *v);
        let y = b.get(i).map_or(0, |v| *v);
        let b1 = x < y;
        let d1 = sub_u64(x, y);
        let b2 = d1 < u64::from(borrow);
        out.push(sub_u64(d1, u64::from(borrow)));
        borrow = b1 || b2;
    }
    out
}

/// Exact order of `a` and `b`; Err when aligning the four terms needs more than `max_shift` bits.
pub fn cmp_points(a: &Point, b: &Point, max_shift: u64) -> Result<Ordering, String> {
    // (term, flip): b's terms enter with their sign flipped, so the sum is a - b.
    let terms = [
        (&a.base, false),
        (&a.off, false),
        (&b.base, true),
        (&b.off, true),
    ];
    let e0 = terms
        .iter()
        .filter(|(t, _)| !t.is_zero())
        .map(|(t, _)| t.exp)
        .min();
    let e0 = match e0 {
        None => return Ok(Ordering::Equal),
        Some(e) => e,
    };
    let (mut neg, mut acc): (bool, Vec<u64>) = (false, vec![0]);
    for (t, flip) in terms.iter().filter(|(t, _)| !t.is_zero()) {
        let t_neg = t.neg != *flip;
        let d = t
            .exp
            .checked_sub(e0)
            .ok_or_else(|| String::from("exponent spread overflow"))?;
        let s = u64::try_from(d).map_err(|_| String::from("negative shift"))?;
        if s > max_shift {
            return Err(format!("exponent spread {} exceeds {} bits", s, max_shift));
        }
        let m = shifted(t.mant, s)?;
        if t_neg == neg {
            acc = mag_add(&acc, &m);
        } else if mag_cmp(&acc, &m) != Ordering::Less {
            acc = mag_sub(&acc, &m);
        } else {
            acc = mag_sub(&m, &acc);
            neg = t_neg;
        }
    }
    if acc.iter().all(|w| *w == 0) {
        Ok(Ordering::Equal)
    } else if neg {
        Ok(Ordering::Less)
    } else {
        Ok(Ordering::Greater)
    }
}
