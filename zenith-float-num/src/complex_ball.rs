//! Midpoint-radius (ball) arithmetic for complex formulas without a closed error bound (ZF-18/68).
//!
//! A [`RB`] is a value `v` and a radius `r` with `|true − v| ≤ r`; a [`CB`] holds one ball per
//! part. Each operation forms its value at the working width `q` with [`RoundingMode::None`] and
//! adds to the radius (1) the propagated input radii, from the operation's derivative bound over
//! the ball, and (2) its own rounding error. The library's real operations at `(q, None)` are
//! truncations of Ziv-certified values, within 2 ulps at `q`; the bound used is 4 ulps.
//! [`ziv_complex_ball`] then claims for each part only the bits its radius leaves, so
//! cancellation and tiny parts lower the claim and force a wider retry instead of being
//! certified. Every bound is an upper bound; an operation that cannot bound its result (a
//! ball that touches a singularity or a branch cut, a non-finite or subnormal value) returns
//! `None` and the loop retries wider or refuses with `PrecisionRetryExhausted`.

use crate::common::util::bump_prec_retry;
use crate::common::util::round_p;
use crate::defs::Exponent;
use crate::Consts;
use crate::Error;
use crate::ExactComplex;
use crate::ExactNum;
use crate::RoundingMode;
use crate::WORD_BIT_SIZE;
use core::ops::Neg;

const N: RoundingMode = RoundingMode::None;

/// Upper bound `m · 2^e` with `m < 2^32` (`m = 0` is zero). Every operation rounds up.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Mag {
    m: u64,
    e: Exponent,
}

impl Mag {
    pub(crate) const ZERO: Mag = Mag { m: 0, e: 0 };

    /// `2^k`.
    pub(crate) fn pow2(k: Exponent) -> Mag {
        Mag {
            m: 1 << 31,
            e: k.saturating_sub(31),
        }
    }

    fn norm(m: u128, e: Exponent) -> Mag {
        if m == 0 {
            return Mag::ZERO;
        }
        let bits = 128 - m.leading_zeros() as i64;
        if bits > 32 {
            let sh = bits - 32;
            let lost = m & ((1u128 << sh) - 1) != 0;
            let mut mm = (m >> sh) + u128::from(lost);
            let mut ee = e.saturating_add(sh);
            if mm >= 1 << 32 {
                mm >>= 1; // mm was exactly 2^32: no bit lost
                ee = ee.saturating_add(1);
            }
            Mag {
                m: mm as u64,
                e: ee,
            }
        } else {
            let sh = 32 - bits;
            Mag {
                m: (m << sh) as u64,
                e: e.saturating_sub(sh),
            }
        }
    }

    pub(crate) fn is_zero(&self) -> bool {
        self.m == 0
    }

    /// `k` with `self < 2^k`; `None` for zero.
    pub(crate) fn ub_exp(&self) -> Option<Exponent> {
        if self.m == 0 {
            None
        } else {
            Some(self.e.saturating_add(32))
        }
    }

    pub(crate) fn add(self, o: Mag) -> Mag {
        if self.is_zero() {
            return o;
        }
        if o.is_zero() {
            return self;
        }
        let (a, b) = if self.e >= o.e { (self, o) } else { (o, self) };
        let d = a.e.saturating_sub(b.e);
        if d > 64 {
            // b < 2^(b.e+32) ≤ 2^(a.e-32): below one unit of a's mantissa.
            return Mag::norm(u128::from(a.m) + 1, a.e);
        }
        Mag::norm((u128::from(a.m) << d) + u128::from(b.m), b.e)
    }

    pub(crate) fn mul(self, o: Mag) -> Mag {
        if self.is_zero() || o.is_zero() {
            return Mag::ZERO;
        }
        Mag::norm(
            u128::from(self.m) * u128::from(o.m),
            self.e.saturating_add(o.e),
        )
    }

    pub(crate) fn mul_pow2(self, k: Exponent) -> Mag {
        if self.is_zero() {
            return self;
        }
        Mag {
            m: self.m,
            e: self.e.saturating_add(k),
        }
    }

    pub(crate) fn mul_u64(self, k: u64) -> Mag {
        self.mul(Mag::norm(u128::from(k), 0))
    }
}

/// Real ball: `|true − v| ≤ r`.
#[derive(Clone, Debug)]
pub(crate) struct RB {
    pub(crate) v: ExactNum,
    pub(crate) r: Mag,
}

/// Complex ball: one real ball per part.
#[derive(Clone, Debug)]
pub(crate) struct CB {
    pub(crate) re: RB,
    pub(crate) im: RB,
}

fn usable(v: &ExactNum) -> bool {
    !v.is_nan() && !v.is_inf() && !v.is_subnormal()
}

/// Exponent `e` with `2^(e−1) ≤ |v| < 2^e`; `None` for zero.
fn ex(v: &ExactNum) -> Option<Exponent> {
    if v.is_zero() {
        None
    } else {
        v.exponent()
    }
}

/// Rounding error bound of a value formed at `q`: 4 ulps; zero and unflagged results are exact.
fn rnd(v: &ExactNum, q: usize) -> Mag {
    // The library's inexact flag reports any rounding by the operation (and the operands'
    // history): an unflagged result is exact (exact sums, products, roots, powers, `ln 1`).
    if !v.inexact() {
        return Mag::ZERO;
    }
    match ex(v) {
        Some(e) => Mag::pow2(e.saturating_sub(q as Exponent).saturating_add(2)),
        None => Mag::ZERO,
    }
}

fn ball(v: ExactNum, r: Mag) -> Option<RB> {
    if usable(&v) {
        Some(RB { v, r })
    } else {
        None
    }
}

impl RB {
    /// An exact value.
    pub(crate) fn exact(v: ExactNum) -> RB {
        RB { v, r: Mag::ZERO }
    }

    /// A value rounded once at `q` (a correctly rounded constant or real function).
    pub(crate) fn rounded(v: ExactNum, q: usize) -> Option<RB> {
        let r = rnd(&v, q);
        ball(v, r)
    }

    pub(crate) fn is_exact(&self) -> bool {
        self.r.is_zero()
    }

    /// Upper bound of `|true|`.
    pub(crate) fn abs_ub(&self) -> Mag {
        match ex(&self.v) {
            Some(e) => Mag::pow2(e).add(self.r),
            None => self.r,
        }
    }

    /// `l` with `|x| ≥ 2^l` for every `x` in the ball, when the ball keeps clear of zero
    /// (radius at most an eighth of `|v|`).
    pub(crate) fn lo_exp(&self) -> Option<Exponent> {
        let e = ex(&self.v)?;
        match self.r.ub_exp() {
            None => Some(e - 1),
            Some(u) if u <= e - 3 => Some(e - 2),
            _ => None,
        }
    }

    /// May the ball contain zero?
    pub(crate) fn may_be_zero(&self) -> bool {
        if self.v.is_zero() {
            return true;
        }
        match (ex(&self.v), self.r.ub_exp()) {
            (_, None) => false,
            (Some(e), Some(u)) => u >= e - 1,
            _ => true,
        }
    }

    /// May the ball contain a negative value (or is it the negative half)?
    fn may_be_neg(&self) -> bool {
        (self.v.is_negative() && !self.v.is_zero()) || (!self.r.is_zero() && self.may_be_zero())
    }

    fn may_be_nonpos(&self) -> bool {
        self.v.is_zero() || self.v.is_negative() || (!self.r.is_zero() && self.may_be_zero())
    }

    pub(crate) fn neg(self) -> RB {
        RB {
            v: self.v.neg(),
            r: self.r,
        }
    }

    pub(crate) fn abs(self) -> RB {
        RB {
            v: self.v.abs(),
            r: self.r,
        }
    }

    /// `self · 2^k` (exact).
    pub(crate) fn scale2(self, k: Exponent) -> Option<RB> {
        let pv = self.v.mantissa_max_bit_len().unwrap_or(WORD_BIT_SIZE);
        let v = self.v.ldexp(k, pv, N);
        ball(v, self.r.mul_pow2(k))
    }

    pub(crate) fn add(&self, o: &RB, q: usize) -> Option<RB> {
        let v = ExactNum::add_ref(&self.v, &o.v, q, N);
        let r = self.r.add(o.r).add(rnd(&v, q));
        ball(v, r)
    }

    pub(crate) fn sub(&self, o: &RB, q: usize) -> Option<RB> {
        let v = ExactNum::sub_ref(&self.v, &o.v, q, N);
        let r = self.r.add(o.r).add(rnd(&v, q));
        ball(v, r)
    }

    pub(crate) fn mul(&self, o: &RB, q: usize) -> Option<RB> {
        let v = ExactNum::mul_ref(&self.v, &o.v, q, N);
        if v.is_zero() && !self.v.is_zero() && !o.v.is_zero() {
            return None; // underflow: a zero here is not exact
        }
        let mut r = rnd(&v, q).add(self.r.mul(o.r));
        if let Some(ea) = ex(&self.v) {
            r = r.add(o.r.mul_pow2(ea));
        }
        if let Some(eb) = ex(&o.v) {
            r = r.add(self.r.mul_pow2(eb));
        }
        ball(v, r)
    }

    pub(crate) fn sqr(&self, q: usize) -> Option<RB> {
        self.mul(self, q)
    }

    /// `self / o`; `None` when the divisor ball may reach zero.
    pub(crate) fn div(&self, o: &RB, q: usize) -> Option<RB> {
        let lb = o.lo_exp()?;
        let eb = ex(&o.v)?;
        let v = ExactNum::div_ref(&self.v, &o.v, q, N);
        if v.is_zero() && !self.v.is_zero() {
            return None; // underflow
        }
        // |a/b − A/B| ≤ (|a|·|B−b| + |b|·|A−a|) / (|b|·|B|), |b| ≥ 2^(eb−1), |B| ≥ 2^lb.
        let mut r = rnd(&v, q).add(self.r.mul_pow2(-lb));
        if let Some(ea) = ex(&self.v) {
            r = r.add(o.r.mul_pow2(ea - (eb - 1) - lb));
        }
        ball(v, r)
    }

    pub(crate) fn sqrt(&self, q: usize) -> Option<RB> {
        if self.v.is_negative() && !self.v.is_zero() {
            return None;
        }
        let v = self.v.sqrt(q, N);
        if self.r.is_zero() {
            return RB::rounded(v, q);
        }
        // |√A − √a| ≤ |A − a| / √min, min ≥ 2^l.
        let l = self.lo_exp()?;
        let r = rnd(&v, q).add(self.r.mul_pow2(-l.div_euclid(2)));
        ball(v, r)
    }

    pub(crate) fn ln(&self, q: usize, cc: &mut Consts) -> Option<RB> {
        if self.v.is_negative() || self.v.is_zero() {
            return None;
        }
        let v = self.v.ln(q, N, cc);
        if self.r.is_zero() {
            return RB::rounded(v, q);
        }
        let l = self.lo_exp()?;
        let r = rnd(&v, q).add(self.r.mul_pow2(-l));
        ball(v, r)
    }

    /// `ln(1 + self)`: the derivative `1/(1 + x)` is bounded over the ball of `1 + self`.
    pub(crate) fn log1p(&self, q: usize, cc: &mut Consts) -> Option<RB> {
        let v = self.v.log1p(q, N, cc);
        if self.r.is_zero() {
            return RB::rounded(v, q);
        }
        let w = RB::exact(ExactNum::from_u8(1, WORD_BIT_SIZE)).add(self, q)?;
        if w.v.is_negative() {
            return None;
        }
        let l = w.lo_exp()?;
        let r = rnd(&v, q).add(self.r.mul_pow2(-l));
        ball(v, r)
    }

    pub(crate) fn exp(&self, q: usize, cc: &mut Consts) -> Option<RB> {
        let v = self.v.exp(q, N, cc);
        if v.is_zero() {
            return None; // underflow: e^x is never zero
        }
        if self.r.is_zero() {
            return RB::rounded(v, q);
        }
        // |e^A − e^a| ≤ e^a |δ| e^|δ| < 2^(ev+1) |δ| for |δ| ≤ 1/2.
        if self.r.ub_exp()? > -1 {
            return None;
        }
        let ev = ex(&v)?;
        let r = rnd(&v, q).add(self.r.mul_pow2(ev + 1));
        ball(v, r)
    }

    /// `(sin, cos)`: both are 1-Lipschitz.
    pub(crate) fn sin_cos(&self, q: usize, cc: &mut Consts) -> Option<(RB, RB)> {
        let (s, c) = self.v.sin_cos(q, N, cc);
        let rs = rnd(&s, q).add(self.r);
        let rc = rnd(&c, q).add(self.r);
        Some((ball(s, rs)?, ball(c, rc)?))
    }

    /// `atan2(self, x)` (`self` is `y`); `None` when the ball pair may touch the cut
    /// `y = 0, x ≤ 0` or the origin.
    pub(crate) fn atan2(&self, x: &RB, q: usize, cc: &mut Consts) -> Option<RB> {
        let v = self.atan2_value(x, q, cc);
        if self.r.is_zero() && x.r.is_zero() {
            return RB::rounded(v, q);
        }
        if self.r.is_zero() && self.v.is_zero() && !x.may_be_nonpos() {
            return Some(RB::exact(self.v.clone())); // atan2(±0, x > 0) = ±0
        }
        if self.may_be_zero() && x.may_be_nonpos() {
            // A y of exact ±0 keeps the side of the cut; x must stay clear of zero.
            if !(self.r.is_zero() && self.v.is_zero() && !x.may_be_zero()) {
                return None;
            }
            return RB::rounded(v, q);
        }
        let m = match (ex(&self.v), ex(&x.v)) {
            (Some(a), Some(b)) => a.max(b),
            (Some(a), None) => a,
            (None, Some(b)) => b,
            (None, None) => return None,
        };
        let rin = self.r.add(x.r);
        if rin.ub_exp()? > m - 3 {
            return None;
        }
        // |∇atan2| = 1/ρ with ρ ≥ max(|x|, |y|) − (rx + ry) ≥ 2^(m−2).
        let r = rnd(&v, q).add(rin.mul_pow2(2 - m));
        ball(v, r)
    }

    fn atan2_value(&self, x: &RB, q: usize, cc: &mut Consts) -> ExactNum {
        self.v.atan2(&x.v, q, N, cc)
    }

    /// `hypot(self, o)` (1-Lipschitz in each argument).
    pub(crate) fn hypot(&self, o: &RB, q: usize) -> Option<RB> {
        let v = self.v.hypot(&o.v, q, N);
        let r = rnd(&v, q).add(self.r).add(o.r);
        ball(v, r)
    }

    fn may_be_neg_pub(&self) -> bool {
        self.may_be_neg()
    }
}

impl CB {
    pub(crate) fn exact(z: &ExactComplex) -> CB {
        CB {
            re: RB::exact(z.re().clone()),
            im: RB::exact(z.im().clone()),
        }
    }

    pub(crate) fn from_parts(re: RB, im: RB) -> CB {
        CB { re, im }
    }

    /// Real ball as `x + 0i` (exact `+0` imaginary part).
    pub(crate) fn from_real(x: RB) -> CB {
        CB {
            re: x,
            im: RB::exact(ExactNum::new(WORD_BIT_SIZE)),
        }
    }

    pub(crate) fn int(k: u32) -> CB {
        CB::from_real(RB::exact(ExactNum::from_u32(k, WORD_BIT_SIZE)))
    }

    pub(crate) fn i() -> CB {
        CB {
            re: RB::exact(ExactNum::new(WORD_BIT_SIZE)),
            im: RB::exact(ExactNum::from_u8(1, WORD_BIT_SIZE)),
        }
    }

    pub(crate) fn add(&self, o: &CB, q: usize) -> Option<CB> {
        Some(CB {
            re: self.re.add(&o.re, q)?,
            im: self.im.add(&o.im, q)?,
        })
    }

    pub(crate) fn sub(&self, o: &CB, q: usize) -> Option<CB> {
        Some(CB {
            re: self.re.sub(&o.re, q)?,
            im: self.im.sub(&o.im, q)?,
        })
    }

    pub(crate) fn neg(self) -> CB {
        CB {
            re: self.re.neg(),
            im: self.im.neg(),
        }
    }

    pub(crate) fn conj(self) -> CB {
        CB {
            re: self.re,
            im: self.im.neg(),
        }
    }

    /// `i · self` (exact).
    pub(crate) fn mul_i(self) -> CB {
        CB {
            re: self.im.neg(),
            im: self.re,
        }
    }

    /// `−i · self` (exact).
    pub(crate) fn mul_neg_i(self) -> CB {
        CB {
            re: self.im,
            im: self.re.neg(),
        }
    }

    pub(crate) fn scale2(self, k: Exponent) -> Option<CB> {
        Some(CB {
            re: self.re.scale2(k)?,
            im: self.im.scale2(k)?,
        })
    }

    pub(crate) fn mul(&self, o: &CB, q: usize) -> Option<CB> {
        let ac = self.re.mul(&o.re, q)?;
        let bd = self.im.mul(&o.im, q)?;
        let ad = self.re.mul(&o.im, q)?;
        let bc = self.im.mul(&o.re, q)?;
        Some(CB {
            re: ac.sub(&bd, q)?,
            im: ad.add(&bc, q)?,
        })
    }

    pub(crate) fn mul_real(&self, x: &RB, q: usize) -> Option<CB> {
        Some(CB {
            re: self.re.mul(x, q)?,
            im: self.im.mul(x, q)?,
        })
    }

    pub(crate) fn div_real(&self, x: &RB, q: usize) -> Option<CB> {
        Some(CB {
            re: self.re.div(x, q)?,
            im: self.im.div(x, q)?,
        })
    }

    pub(crate) fn div(&self, o: &CB, q: usize) -> Option<CB> {
        let den = o.re.sqr(q)?.add(&o.im.sqr(q)?, q)?;
        let nre = self.re.mul(&o.re, q)?.add(&self.im.mul(&o.im, q)?, q)?;
        let nim = self.im.mul(&o.re, q)?.sub(&self.re.mul(&o.im, q)?, q)?;
        Some(CB {
            re: nre.div(&den, q)?,
            im: nim.div(&den, q)?,
        })
    }

    /// Upper bound of `|true|` (`|re| + |im|`).
    pub(crate) fn abs_ub(&self) -> Mag {
        self.re.abs_ub().add(self.im.abs_ub())
    }

    /// Add `r` to both radii (an error of modulus at most `r`).
    pub(crate) fn widen(mut self, r: Mag) -> CB {
        self.re.r = self.re.r.add(r);
        self.im.r = self.im.r.add(r);
        self
    }

    pub(crate) fn exp(&self, q: usize, cc: &mut Consts) -> Option<CB> {
        let e = self.re.exp(q, cc)?;
        let (s, c) = self.im.sin_cos(q, cc)?;
        Some(CB {
            re: e.mul(&c, q)?,
            im: e.mul(&s, q)?,
        })
    }

    /// Principal `ln`: `ln|big| + ½ log1p((small/big)²) + i atan2(im, re)`.
    pub(crate) fn ln(&self, q: usize, cc: &mut Consts) -> Option<CB> {
        let are = self.re.clone().abs();
        let aim = self.im.clone().abs();
        let re_big = match (ex(&are.v), ex(&aim.v)) {
            (Some(a), Some(b)) => a >= b,
            (Some(_), None) => true,
            (None, Some(_)) => false,
            (None, None) => return None,
        };
        let (big, small) = if re_big { (are, aim) } else { (aim, are) };
        let ln_big = big.ln(q, cc)?;
        let re = if small.v.is_zero() && small.r.is_zero() {
            ln_big
        } else {
            let t = small.div(&big, q)?.sqr(q)?;
            let corr = t.log1p(q, cc)?.scale2(-1)?;
            ln_big.add(&corr, q)?
        };
        let im = self.im.atan2(&self.re, q, cc)?;
        Some(CB { re, im })
    }

    /// Principal `√`, Kahan form: `h = (|z| + |x|)/2`, `t = √h`, `u = |y|/(2t)`; `(t, ±u)` for
    /// `x ≥ 0` and `(u, ±t)` for `x < 0`. The two forms agree at `x = 0` up to
    /// `|t − u| ≤ |x|/t`, added when the ball of `x` may change sign.
    pub(crate) fn sqrt(&self, q: usize) -> Option<CB> {
        let (x, y) = (&self.re, &self.im);
        if x.is_exact() && y.is_exact() && x.v.is_zero() && y.v.is_zero() {
            // C99 csqrt(±0 ± 0i) = +0 ± 0i.
            return Some(CB {
                re: RB::exact(x.v.clone().abs()),
                im: RB::exact(y.v.clone()),
            });
        }
        if y.may_be_zero() && !y.r.is_zero() && x.may_be_neg_pub() {
            return None; // the ball may straddle the cut
        }
        let a = x.hypot(y, q)?;
        let ax = x.clone().abs();
        let h = a.add(&ax, q)?.scale2(-1)?;
        let t = h.sqrt(q)?;
        let ay = y.clone().abs();
        let u = if ay.v.is_zero() && ay.r.is_zero() {
            RB::exact(ExactNum::new(WORD_BIT_SIZE))
        } else {
            ay.div(&t.clone().scale2(1)?, q)?
        };
        let swap = if !x.r.is_zero() && x.may_be_zero() {
            let et = t.lo_exp()?;
            x.abs_ub().mul_pow2(-et)
        } else {
            Mag::ZERO
        };
        // A y ball across zero (x ≥ 0 side only, checked above) may take the other sign of
        // ±u: |true − center| ≤ (|Y| + |y|)/(2t) ≤ |y|_ub / t.
        let flip = if !y.r.is_zero() && y.may_be_zero() {
            let et = t.lo_exp()?;
            y.abs_ub().mul_pow2(-et)
        } else {
            Mag::ZERO
        };
        let neg_im = y.v.is_negative();
        let (re, mut im) = if x.v.is_negative() && !x.v.is_zero() {
            (u, if neg_im { t.neg() } else { t })
        } else {
            (t, if neg_im { u.neg() } else { u })
        };
        im.r = im.r.add(flip);
        Some(CB { re, im }.widen(swap))
    }
}

/// Round one part of a ball at `(p, rm)`; `None` if the radius leaves too few bits.
fn round_part(b: &RB, p: usize, rm: RoundingMode) -> Option<ExactNum> {
    let mut v = b.v.clone();
    if b.r.is_zero() {
        v.set_inexact(false);
        v.set_precision(p, rm).ok()?;
        return Some(v);
    }
    let e = ex(&v)?;
    let u = b.r.ub_exp()?;
    // |true − v| < 2^u ≤ 2^(e − s − 2): `s` leading bits of v are certified (MPFR-style).
    let s = e.saturating_sub(u).saturating_sub(2);
    let len = v.mantissa_max_bit_len()?;
    if s <= p as Exponent {
        return None;
    }
    let s = (s as usize).min(len) / WORD_BIT_SIZE * WORD_BIT_SIZE;
    if s <= p {
        return None;
    }
    v.set_inexact(true);
    if v.try_set_precision(p, rm, s) {
        Some(v)
    } else {
        None
    }
}

/// Ziv loop over ball evaluations: `compute(q)` returns a ball at working width `q`; each part
/// is rounded only when its radius certifies the rounding. `None` balls (unbounded at this
/// width) retry wider; retry exhaustion refuses with `PrecisionRetryExhausted`.
pub(crate) fn ziv_complex_ball<F>(p: usize, rm: RoundingMode, mut compute: F) -> ExactComplex
where
    F: FnMut(usize) -> Option<CB>,
{
    if p == 0 {
        return crate::complex_special::nan_pair(Error::InvalidArgument);
    }
    let p = round_p(p);
    if let Err(e) = crate::common::util::reserve_probe(p) {
        return crate::complex_special::nan_pair(e);
    }
    let mut p_inc = WORD_BIT_SIZE;
    let Some(mut p_wrk) = p.checked_add(p_inc) else {
        return crate::complex_special::nan_pair(Error::InvalidArgument);
    };
    p_wrk = round_p(p_wrk);
    loop {
        let Some(q) = p_wrk.checked_add(WORD_BIT_SIZE.saturating_mul(2)) else {
            return crate::complex_special::nan_pair(Error::InvalidArgument);
        };
        if let Some(b) = compute(q) {
            if let (Some(re), Some(im)) = (round_part(&b.re, p, rm), round_part(&b.im, p, rm)) {
                return ExactComplex::new(re, im);
            }
        }
        if bump_prec_retry(&mut p_wrk, &mut p_inc, p).is_err() {
            return crate::complex_special::nan_pair(Error::PrecisionRetryExhausted);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mag_rounds_up() {
        let a = Mag::pow2(0);
        let b = Mag::pow2(-200);
        let s = a.add(b);
        assert!(s.ub_exp() == Some(1));
        assert!(s != a); // the tiny term is not dropped
        let t = Mag::pow2(3).mul(Mag::pow2(-5));
        assert_eq!(t, Mag::pow2(-2));
        assert_eq!(Mag::pow2(4).ub_exp(), Some(5));
    }

    #[test]
    fn cancellation_is_not_certified() {
        // (1 + 2^-300) − 1 at q = 192: the radius of 1 + 2^-300 (rounded) covers the result.
        let q = 192;
        let one = RB::exact(ExactNum::from_u8(1, q));
        let tiny = RB::exact(ExactNum::from_u8(1, q).ldexp(-300, q, N));
        let s = one.add(&tiny, q).unwrap().sub(&one, q).unwrap();
        assert!(round_part(&s, 64, RoundingMode::ToEven).is_none());
    }
}
