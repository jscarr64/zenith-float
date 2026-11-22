//! Slice CB: complex balls `CBall { re, im }`, a closed rectangle of two independent real [`Ball`] parts (CONTRACT 16.CB).
//! `add`, `sub`, `sqr` and `abs` are built on the public `Ball` API and S1. `mul` (16.CB.4) and `div` (16.CB.7) are
//! faithful ports of Arb `acb_mul` / `acb_div` (FLINT 3.6.0), with the `arb_mul`, `arb_addmul`, `arb_div` and
//! `acb_inv` radius bounds they use, so they carry their own `Mag` radius arithmetic (Jeff ruling 2026-11-17).
//! Ownership is O9 (16.CB.0): an op takes `self` and borrows its second operand.

use zenith_float_num::{Error, ExactNum, RoundingMode, EXPONENT_MAX, EXPONENT_MIN};

use std::ops::Neg;

use crate::ball::Ball;
use crate::enclosure::Enclosure;
use crate::error::{BoundError, BoundResult};
use crate::mag::Mag;
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

    /// `a²` (16.CB.5), Arb `acb_sqr` = `acb_mul(a, a)` (FLINT 3.6.0 `src/acb/mul.c`). An exact-zero `im` gives
    /// `re = sqr(re)`, `im = im·re`; an exact-zero `re` gives `re = −sqr(im)`, `im = re·im` (`Ball::sqr` for Arb's
    /// aliased `arb_mul(x, x)`, no wider). Otherwise the midpoints `re² − im²` and `2·re·im` are formed from exact
    /// products and rounded once (Arb `arf_complex_sqr`), with Arb's radii (`_acb_sqr_slow`, `re` grouped per part)
    /// `2(|re| r(re) + |im| r(im)) + r(re)² + r(im)²` and `2(|im| r(re) + |re| r(im) + r(re) r(im))`, plus ε. Each part is
    /// no wider than the same part of `mul(a, a)`: the same midpoints, and radii that drop `mul`'s duplicated cross terms.
    pub fn sqr(self, prec: Precision) -> BoundResult<CBall> {
        const OP: &str = "cball_sqr";
        let (re, im) = acb_sqr(self.re, self.im, prec).map_err(|e| remap(e, OP))?;
        CBall::new(re, im).map_err(|e| remap(e, OP))
    }

    /// `a·b`, Arb `acb_mul` (FLINT 3.6.0 `src/acb/mul.c`, 16.CB.4). An exact-zero part (Arb `arb_is_zero`), tested in
    /// Arb's order `a.im`, `b.im`, `a.re`, `b.re`, reduces the product to two `Ball::mul` (Arb `arb_mul`, the same
    /// radius formula). Otherwise each midpoint, `a.re·b.re − a.im·b.im` and `a.re·b.im + a.im·b.re`, is formed from
    /// exact products and rounded once (Arb `arf_complex_mul`), and the radii are Arb's six-term bounds (`_acb_mul_slow`)
    /// in `Mag` rounding up, plus ε (9.4, O19).
    pub fn mul(self, other: &CBall, prec: Precision) -> BoundResult<CBall> {
        const OP: &str = "cball_mul";
        let (re, im) =
            acb_mul(self.re, self.im, &other.re, &other.im, prec).map_err(|e| remap(e, OP))?;
        CBall::new(re, im).map_err(|e| remap(e, OP))
    }

    /// `a / b`, Arb `acb_div` (FLINT 3.6.0 `src/acb/div.c`, 16.CB.7). An exact-zero part of `b` gives part-wise Arb
    /// `arb_div` (an exact-zero part of `a` gives an exact-zero part); an exact `b` of at most `prec/2` bits per part with
    /// `prec > 256` gives `(a·conj(b)) / (b.re² + b.im²)` by Arb `arb_mul`/`arb_addmul`/`arb_submul`/`arb_div`;
    /// otherwise `a · inv(b)` with Arb `acb_inv` and `acb_mul`. Arb's indeterminate (`nan ± inf`) and `0 ± inf` parts are
    /// all reals here (a zf midpoint is finite, 6.2), so a divisor that reaches 0 gives all reals in each part Arb
    /// cannot bound, never `Domain`.
    pub fn div(self, other: &CBall, prec: Precision) -> BoundResult<CBall> {
        const OP: &str = "cball_div";
        let (re, im) =
            acb_div(self.re, self.im, &other.re, &other.im, prec).map_err(|e| remap(e, OP))?;
        CBall::new(re, im).map_err(|e| remap(e, OP))
    }
}

/// Arb `acb_div` takes its exact-divisor branch above this prec (FLINT 3.6.0 `src/acb/div.c`).
const ACB_DIV_EXACT_PREC: Precision = 256;
/// Arb `MAG_BITS`: `arb_div` keeps this many midpoint bits past the input accuracy; `acb_inv`'s `|b|²` prec floor.
const ARB_MAG_BITS: i64 = 30;
/// `arb_div` switches to the interval quotient `arb_div_wide` at or below this many bits of relative accuracy.
const ARB_DIV_WIDE_ACC: i64 = 20;
/// Arb's lagom exponent bound (`ARF_MAX_LAGOM_EXP` = `MAG_MAX_LAGOM_EXP` = `COEFF_MAX / 4` on a 64-bit limb).
const ARB_LAGOM_EXP: i64 = (1 << 60) - 1;
/// `1/2` as a `Mag` mantissa (Arb `MAG_ONE_HALF`).
const MAG_HALF: u64 = 0x8000_0000_0000_0000;
/// `1 + 2^-16` as a `Mag` mantissa at exponent 1 (Arb `arb_div`: `MAG_ONE_HALF + (MAG_ONE_HALF >> 16)`, exponent 1).
const ONE_PLUS_2M16: u64 = MAG_HALF + (MAG_HALF >> 16);
/// Bits kept for `Mag`-valued sums and differences formed through `ExactNum` (a `Mag` mantissa is 64 bits).
const MAG_SUM_BITS: Precision = 64;
/// Every midpoint rounds to nearest, ties to even (CONTRACT 15.1, O19).
const MID_RM: RoundingMode = RoundingMode::ToEven;

type Parts = (Ball, Ball);

/// Arb `acb_mul` on parts `a + i b` (consumed) and `c + i d`.
fn acb_mul(a: Ball, b: Ball, c: &Ball, d: &Ball, prec: Precision) -> BoundResult<Parts> {
    if exact_zero(&b) {
        let im = ball_again(d)?.mul(&a, prec)?;
        Ok((ball_again(c)?.mul(&a, prec)?, im))
    } else if exact_zero(d) {
        let im = b.mul(c, prec)?;
        Ok((a.mul(c, prec)?, im))
    } else if exact_zero(&a) {
        // i b (c + i d) = −(d b) + i (c b)
        let re = ball_again(d)?.mul(&b, prec)?.negated();
        Ok((re, ball_again(c)?.mul(&b, prec)?))
    } else if exact_zero(c) {
        // (a + i b) i d = −(b d) + i (a d)
        let re = b.mul(d, prec)?.negated();
        Ok((re, a.mul(d, prec)?))
    } else {
        mul_general(&a, &b, c, d, prec)
    }
}

/// Arb `acb_sqr` (`acb_mul(x, x)`) on parts `a + i b` (consumed).
fn acb_sqr(a: Ball, b: Ball, prec: Precision) -> BoundResult<Parts> {
    if exact_zero(&b) {
        let re = a.sqr(prec)?;
        Ok((re, b.mul(&a, prec)?))
    } else if exact_zero(&a) {
        let re = b.sqr(prec)?.negated();
        Ok((re, a.mul(&b, prec)?))
    } else {
        let (am, bm) = (mag_of(a.mid()), mag_of(b.mid()));
        let (ar, br) = (a.rad(), b.rad());
        // Arb's sum 2(|a|ar + |b|br) + ar² + br², grouped as `(2|a|ar + ar²) + (2|b|br + br²)` (each `Ball::sqr`'s
        // radius) so the 16.CB.5 bound against the four-product form holds after `Mag` rounding too.
        let er = mag_double(&am.mul_up(ar))?
            .add_up(&ar.mul_up(ar))
            .add_up(&mag_double(&bm.mul_up(br))?.add_up(&br.mul_up(br)));
        let fr = mag_double(&bm.mul_up(ar).add_up(&am.mul_up(br)).add_up(&ar.mul_up(br)))?;
        let (aa, bb) = (
            a.mid().mul_full_prec(a.mid()),
            b.mid().mul_full_prec(b.mid()),
        );
        let ab = a.mid().mul_full_prec(b.mid());
        let lost_re = aa.inexact() || bb.inexact();
        let lost_im = ab.inexact();
        let re = finish(aa.sub(&bb, prec, MID_RM), &er, prec, prec, lost_re)?;
        let im = finish(ab.ldexp(1, prec, MID_RM), &fr, prec, prec, lost_im)?;
        Ok((re, im))
    }
}

/// Arb `acb_div` on parts `a + i b` (consumed) and `c + i d`.
fn acb_div(a: Ball, b: Ball, c: &Ball, d: &Ball, prec: Precision) -> BoundResult<Parts> {
    if exact_zero(d) {
        if exact_zero(&b) {
            Ok((arb_div(&a, c, prec)?, zero_ball(prec)?))
        } else if exact_zero(&a) {
            Ok((zero_ball(prec)?, arb_div(&b, c, prec)?))
        } else {
            Ok((arb_div(&a, c, prec)?, arb_div(&b, c, prec)?))
        }
    } else if exact_zero(c) {
        if exact_zero(&b) {
            Ok((zero_ball(prec)?, arb_div(&a, d, prec)?.negated()))
        } else if exact_zero(&a) {
            Ok((arb_div(&b, d, prec)?, zero_ball(prec)?))
        } else {
            Ok((arb_div(&b, d, prec)?, arb_div(&a, d, prec)?.negated()))
        }
    } else if prec > ACB_DIV_EXACT_PREC
        && c.is_exact()
        && d.is_exact()
        && fits_bits(c.mid(), prec / 2)
        && fits_bits(d.mid(), prec / 2)
    {
        let t = arb_addmul(ball_again(c)?.mul(c, prec)?, d, d, prec, false)?;
        let u = arb_addmul(ball_again(&a)?.mul(c, prec)?, &b, d, prec, false)?;
        let v = arb_addmul(ball_again(&b)?.mul(c, prec)?, &a, d, prec, true)?;
        Ok((arb_div(&u, &t, prec)?, arb_div(&v, &t, prec)?))
    } else {
        let (ir, ii) = acb_inv(c, d, prec)?;
        acb_mul(a, b, &ir, &ii, prec)
    }
}

fn zero_ball(prec: Precision) -> BoundResult<Ball> {
    Ball::exact(ExactNum::from(0i64), prec)
}

/// Arb `arb_is_zero`: the exact point 0.
fn exact_zero(b: &Ball) -> bool {
    b.mid().is_zero() && b.rad().is_zero()
}

fn lagom_exp(e: Option<i64>) -> bool {
    e.is_some_and(|e| (-ARB_LAGOM_EXP..=ARB_LAGOM_EXP).contains(&e))
}

/// Arb `ARB_IS_LAGOM`: finite radius, and the midpoint and radius exponents within the lagom range.
fn lagom(b: &Ball) -> bool {
    !b.rad().is_infinity()
        && (b.mid().is_zero() || lagom_exp(b.mid().exponent()))
        && (b.rad().is_zero() || lagom_exp(b.rad().exp()))
}

/// `x` has at most `bits` significant bits (Arb `arf_bits(x) <= bits`): rounding it to `bits` drops nothing.
fn fits_bits(x: &ExactNum, bits: Precision) -> bool {
    x.is_zero() || (bits > 0 && !ExactNum::from(1i64).mul(x, bits, MID_RM).inexact())
}

/// A rounded midpoint, or `None` on overflow (CONTRACT 10.4), as S1's `classify`.
fn classify(op: &'static str, raw: ExactNum) -> BoundResult<Option<ExactNum>> {
    if raw.is_inf() {
        return Ok(None);
    }
    if raw.is_nan() {
        return match raw.err() {
            Some(Error::ExponentOverflow(_)) => Ok(None),
            Some(Error::DivisionByZero) => Err(BoundError::Domain {
                op,
                detail: "zf division by zero",
            }),
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
    // An inexact result in the top binade may have been clamped there by zf (S1 `classify`, CONTRACT 15.1).
    if raw.inexact() && raw.exponent() == Some(EXPONENT_MAX) {
        return Ok(None);
    }
    Ok(Some(raw))
}

/// ε of a midpoint rounded to `wp` bits (CONTRACT 9.4, O19; Arb `arf_mag_add_ulp` under nearest rounding): 0 when
/// exact, `half_ulp` above the bottom binade, `2^EXPONENT_MIN` there.
fn eps(mid: &ExactNum, wp: Precision, lost: bool) -> Mag {
    if !mid.inexact() && !lost {
        return Mag::zero();
    }
    match mid.exponent() {
        Some(e) if e > EXPONENT_MIN && !mid.is_zero() => Mag::half_ulp(e, wp),
        Some(_) | None => Mag::ulp_mag(EXPONENT_MIN, 0),
    }
}

/// `mid ± (prop + ε)` for a midpoint `raw` rounded to `wp` bits; the ball carries `prec`. Overflow of the midpoint or
/// an infinite radius is all reals (10.4; Arb's `± inf`).
fn finish(
    raw: ExactNum,
    prop: &Mag,
    wp: Precision,
    prec: Precision,
    lost: bool,
) -> BoundResult<Ball> {
    if prop.is_infinity() {
        return Ball::all_reals(prec);
    }
    match classify("cball", raw)? {
        None => Ball::all_reals(prec),
        Some(m) => {
            let rad = prop.add_up(&eps(&m, wp, lost));
            Ball::from_mid_rad(m, rad, prec)
        }
    }
}

/// Upper `Mag` of `|x|` (Arb `arf_get_mag`).
fn mag_of(x: &ExactNum) -> Mag {
    Mag::from_exact_num(x)
}

/// Lower bound of `min |x|` over the ball (Arb `arb_get_mag_lower`): 0 when the ball reaches 0.
fn mag_lower_ball(b: &Ball) -> Mag {
    if b.rad().is_infinity() {
        return Mag::zero();
    }
    Mag::from_exact_num_lower(b.mid()).sub_lower(b.rad())
}

/// Upper bound of `max |x|` over the ball (Arb `arb_get_mag`).
fn mag_upper_ball(b: &Ball) -> Mag {
    mag_of(b.mid()).add_up(b.rad())
}

/// `2·m`, exact (Arb `mag_mul_2exp_si(·, ·, 1)`).
fn mag_double(m: &Mag) -> BoundResult<Mag> {
    Ok(m.mul_up(&Mag::finite(MAG_HALF, 2)?))
}

/// Lower bound of `x + y` (Arb `mag_add_lower`); `max(x, y)` when either is outside zf's range.
fn mag_add_lower(x: &Mag, y: &Mag) -> BoundResult<Mag> {
    match (x.to_exact_num(), y.to_exact_num()) {
        (Ok(a), Ok(b)) => Ok(Mag::from_exact_num_lower(&a.add(
            &b,
            MAG_SUM_BITS,
            RoundingMode::Down,
        ))),
        (Err(_), _) | (_, Err(_)) => mag_again(if x >= y { x } else { y }),
    }
}

/// Upper bound of `x − |m|`, clamped at 0 (Arb `mag_sub(·, x, arf_get_mag_lower(m))`).
fn mag_sub_upper(x: &Mag, m: &ExactNum) -> BoundResult<Mag> {
    if x.is_infinity() {
        return Ok(Mag::infinity());
    }
    let d = x.to_exact_num()?.sub(
        &ExactNum::from(1i64)
            .mul(m, MAG_SUM_BITS, RoundingMode::ToZero)
            .abs(),
        MAG_SUM_BITS,
        RoundingMode::Up,
    );
    if d.is_negative() || d.is_zero() {
        return Ok(Mag::zero());
    }
    Ok(Mag::from_exact_num(&d))
}

/// The general branch of Arb `acb_mul` (`_acb_mul_slow`): one rounding per midpoint, six-term radii.
fn mul_general(a: &Ball, b: &Ball, c: &Ball, d: &Ball, prec: Precision) -> BoundResult<Parts> {
    let (am, bm, cm, dm) = (
        mag_of(a.mid()),
        mag_of(b.mid()),
        mag_of(c.mid()),
        mag_of(d.mid()),
    );
    let (ar, br, cr, dr) = (a.rad(), b.rad(), c.rad(), d.rad());
    let er = am
        .mul_up(cr)
        .add_up(&bm.mul_up(dr))
        .add_up(&cm.mul_up(ar))
        .add_up(&dm.mul_up(br))
        .add_up(&ar.mul_up(cr))
        .add_up(&br.mul_up(dr));
    let fr = am
        .mul_up(dr)
        .add_up(&bm.mul_up(cr))
        .add_up(&cm.mul_up(br))
        .add_up(&dm.mul_up(ar))
        .add_up(&br.mul_up(cr))
        .add_up(&ar.mul_up(dr));
    let (ac, bd) = (
        a.mid().mul_full_prec(c.mid()),
        b.mid().mul_full_prec(d.mid()),
    );
    let (ad, bc) = (
        a.mid().mul_full_prec(d.mid()),
        b.mid().mul_full_prec(c.mid()),
    );
    let lost_re = ac.inexact() || bd.inexact();
    let lost_im = ad.inexact() || bc.inexact();
    let re = finish(ac.sub(&bd, prec, MID_RM), &er, prec, prec, lost_re)?;
    let im = finish(ad.add(&bc, prec, MID_RM), &fr, prec, prec, lost_im)?;
    Ok((re, im))
}

/// Arb `arb_addmul` (`sub == false`) / `arb_submul`: `z ± x·y`, one midpoint rounding, radius
/// `r(z) + |x| r(y) + |y| r(x) + r(x) r(y)`.
fn arb_addmul(z: Ball, x: &Ball, y: &Ball, prec: Precision, sub: bool) -> BoundResult<Ball> {
    let (zm, zr, _) = z.into_parts();
    let rad = zr
        .add_up(&mag_of(x.mid()).mul_up(y.rad()))
        .add_up(&mag_of(y.mid()).mul_up(x.rad()))
        .add_up(&x.rad().mul_up(y.rad()));
    let p = x.mid().mul_full_prec(y.mid());
    let lost = p.inexact();
    let raw = if sub {
        zm.sub(&p, prec, MID_RM)
    } else {
        zm.add(&p, prec, MID_RM)
    };
    finish(raw, &rad, prec, prec, lost)
}

/// `xm / ym` rounded to `wp` bits (exact when the quotient is), as S1's `div` midpoint.
fn quot(xm: &ExactNum, ym: &ExactNum, wp: Precision) -> BoundResult<ExactNum> {
    match xm.div_exact(ym, wp) {
        Ok(Some(q)) => Ok(q),
        Ok(None) | Err(_) => Ok(exact_from_parts(xm)?.div(ym, wp, MID_RM)),
    }
}

/// Arb `arb_rel_accuracy_bits`: `exp(mid) − exp(rad) − 1`; an exact ball is `+inf`, a zero midpoint `−inf`.
fn rel_accuracy_bits(b: &Ball) -> i64 {
    if b.rad().is_zero() {
        return i64::MAX;
    }
    if b.mid().is_zero() || b.rad().is_infinity() {
        return i64::MIN;
    }
    match (b.mid().exponent(), b.rad().exp()) {
        (Some(em), Some(er)) => exp_diff(em, er, 1),
        (None, _) | (_, None) => i64::MIN,
    }
}

/// `em − er − k` (exponents of a midpoint and a radius), clamped to the `i64` range like Arb's `slong` accuracy.
fn exp_diff(em: i64, er: i64, k: i64) -> i64 {
    // An i128 difference of three i64 values cannot overflow.
    let d = i128::from(em)
        .checked_sub(i128::from(er))
        .and_then(|d| d.checked_sub(i128::from(k)))
        .unwrap_or(0);
    i64::try_from(d).unwrap_or(if d > 0 { i64::MAX } else { i64::MIN })
}

/// Arb `arb_div`'s reduced working precision: `acc = clamp(min(xacc, yacc), 0, prec)` and
/// `wp = max(min(prec, acc + MAG_BITS), 2)`.
fn reduced_prec(xacc: i64, yacc: i64, prec: Precision) -> (i64, Precision) {
    let p = i64::try_from(prec).unwrap_or(i64::MAX);
    let acc = xacc.min(yacc).max(0).min(p);
    let wp = p.min(acc.checked_add(ARB_MAG_BITS).unwrap_or(p)).max(2);
    (acc, Precision::try_from(wp).unwrap_or(prec))
}

/// Arb `arb_div` (FLINT 3.6.0 `src/arb/div.c`). Arb's Newton midpoint (`WANT_NEWTON`) is a speed path with the same
/// enclosure; zf's correctly rounded quotient is used for the midpoint.
fn arb_div(x: &Ball, y: &Ball, prec: Precision) -> BoundResult<Ball> {
    if y.rad().is_zero() {
        return arb_div_arf(x, y.mid(), prec);
    }
    if y.mid().is_zero() {
        return Ball::all_reals(prec);
    }
    if x.mid().is_zero() && !y.rad().is_infinity() {
        return arb_div_wide(x, y, prec, prec);
    }
    if lagom(x) && lagom(y) {
        let yacc = match (y.mid().exponent(), y.rad().exp()) {
            (Some(em), Some(er)) => exp_diff(em, er, 0),
            (None, _) | (_, None) => i64::MIN,
        };
        let (acc, wp) = reduced_prec(rel_accuracy_bits(x), yacc, prec);
        if acc <= ARB_DIV_WIDE_ACC {
            return arb_div_wide(x, y, wp, prec);
        }
        let t = Mag::from_exact_num_lower(y.mid());
        let Some(q) = classify("cball_div", quot(x.mid(), y.mid(), wp)?)? else {
            return Ball::all_reals(prec);
        };
        let v = x.rad().add_up(&y.rad().mul_up(&mag_of(&q)));
        let rad = v.div_up(&t).mul_up(&Mag::finite(ONE_PLUS_2M16, 1)?);
        return finish(q, &rad, wp, prec, false);
    }
    if x.rad().is_infinity() || y.rad().is_infinity() {
        // Arb: finite / [± inf] is indeterminate; [± inf] / y is [± inf] or indeterminate.
        return Ball::all_reals(prec);
    }
    let (acc, wp) = reduced_prec(rel_accuracy_bits(x), rel_accuracy_bits(y), prec);
    if acc <= ARB_DIV_WIDE_ACC {
        return arb_div_wide(x, y, wp, prec);
    }
    let zr = mag_of(x.mid())
        .mul_up(y.rad())
        .add_up(&mag_of(y.mid()).mul_up(x.rad()));
    let yl = Mag::from_exact_num_lower(y.mid()).mul_lower(&mag_lower_ball(y));
    finish(
        quot(x.mid(), y.mid(), wp)?,
        &zr.div_up(&yl),
        wp,
        prec,
        false,
    )
}

/// Arb `arb_div_arf`: division by the exact `ym`.
fn arb_div_arf(x: &Ball, ym: &ExactNum, prec: Precision) -> BoundResult<Ball> {
    if ym.is_zero() || x.rad().is_infinity() {
        return Ball::all_reals(prec);
    }
    if x.rad().is_zero() {
        return finish(quot(x.mid(), ym, prec)?, &Mag::zero(), prec, prec, false);
    }
    let rad = x.rad().div_up(&Mag::from_exact_num_lower(ym));
    if x.mid().is_zero() {
        return finish(ExactNum::from(0i64), &rad, prec, prec, false);
    }
    finish(quot(x.mid(), ym, prec)?, &rad, prec, prec, false)
}

/// Arb `arb_div_wide`: the interval quotient for a finite `y`, midpoint rounded to `wp` bits; the ball carries `prec`.
fn arb_div_wide(x: &Ball, y: &Ball, wp: Precision, prec: Precision) -> BoundResult<Ball> {
    let t = mag_lower_ball(y);
    if t.is_zero() {
        return Ball::all_reals(prec);
    }
    if x.mid().is_zero() {
        return finish(ExactNum::from(0i64), &x.rad().div_up(&t), wp, prec, false);
    }
    let neg = x.mid().is_negative() != y.mid().is_negative();
    if !reaches_zero(x) {
        // [a, b] / [t, u] = [a/u, b/t], negated for opposite signs
        let lo = mag_lower_ball(x).div_lower(&mag_upper_ball(y));
        let hi = mag_upper_ball(x).div_up(&t);
        let z = set_interval(&lo, &hi, false, wp, prec)?;
        return Ok(if neg { z.negated() } else { z });
    }
    // [-a, b] / [t, u] = [-a/t, b/t];  [-a, b] / -[t, u] = [-b/t, a/t]
    let b = mag_upper_ball(x);
    let a = mag_sub_upper(x.rad(), x.mid())?;
    let (a, b) = if neg { (b, a) } else { (a, b) };
    set_interval(&a.div_up(&t), &b.div_up(&t), true, wp, prec)
}

/// `r > |m|`, decided exactly: Arb `arb_div_wide`'s straddling branch (`arf_cmpabs_mag(mid, rad) < 0`); an infinite
/// radius straddles.
fn reaches_zero(b: &Ball) -> bool {
    match b.rad().cmp_exact_num(b.mid()) {
        Some(core::cmp::Ordering::Greater) | None => true,
        Some(core::cmp::Ordering::Less | core::cmp::Ordering::Equal) => false,
    }
}

/// Arb `arb_set_interval_mag` (`[lo, hi]`) or `arb_set_interval_neg_pos_mag` (`[−lo, hi]`): midpoint `(±lo + hi)/2`
/// rounded once to `wp` bits, radius the larger distance to an end, rounded up; the ball carries `prec`. An infinite
/// end is all reals.
fn set_interval(
    lo: &Mag,
    hi: &Mag,
    neg_lo: bool,
    wp: Precision,
    prec: Precision,
) -> BoundResult<Ball> {
    let (Ok(l), Ok(h)) = (lo.to_exact_num(), hi.to_exact_num()) else {
        return Ball::all_reals(prec);
    };
    let l = if neg_lo {
        ExactNum::from(0i64).sub_full_prec(&l)
    } else {
        l
    };
    let Some(mid) = classify(
        "cball_div",
        exact_from_parts(&l)?
            .add(&h, wp, MID_RM)
            .ldexp(-1, wp, MID_RM),
    )?
    else {
        return Ball::all_reals(prec);
    };
    let rad = dist_up(&h, &mid)?.max(dist_up(&mid, &l)?);
    Ball::from_mid_rad(mid, rad, prec)
}

/// Upper `Mag` of `|x − y|`, rounded once away from zero (no full-precision difference of far-apart values).
fn dist_up(x: &ExactNum, y: &ExactNum) -> BoundResult<Mag> {
    Ok(Mag::from_exact_num(
        &exact_from_parts(x)?
            .sub(y, MAG_SUM_BITS, RoundingMode::FromZero)
            .abs(),
    ))
}

/// Arb `acb_inv` (FLINT 3.6.0 `src/acb/inv.c`) of `c + i d`.
fn acb_inv(c: &Ball, d: &Ball, prec: Precision) -> BoundResult<Parts> {
    if exact_zero(d) {
        return Ok((
            arb_div(&Ball::exact(ExactNum::from(1i64), prec)?, c, prec)?,
            zero_ball(prec)?,
        ));
    }
    if exact_zero(c) {
        let im = arb_div(&Ball::exact(ExactNum::from(1i64), prec)?, d, prec)?.negated();
        return Ok((zero_ball(prec)?, im));
    }
    let all = || -> BoundResult<Parts> { Ok((Ball::all_reals(prec)?, Ball::all_reals(prec)?)) };
    if c.rad().is_infinity() || d.rad().is_infinity() {
        return all();
    }
    let hprec = prec
        .checked_add(3)
        .unwrap_or(prec)
        .max(Precision::try_from(ARB_MAG_BITS).unwrap_or(30));
    let (a, b) = (c.mid(), d.mid());
    let (aa, bb) = (a.mul_full_prec(a), b.mul_full_prec(b));
    let lost = aa.inexact() || bb.inexact();
    let raw = aa.add(&bb, hprec, RoundingMode::Down);
    let inexact = lost || raw.inexact();
    let Some(a2b2) = classify("cball_div", raw)? else {
        return all();
    };
    if a2b2.is_zero() {
        return all();
    }
    if c.rad().is_zero() && d.rad().is_zero() {
        let re = div_rounded_den(a, &a2b2, inexact, &Mag::zero(), prec)?;
        let im = div_rounded_den(b, &a2b2, inexact, &Mag::zero(), prec)?.negated();
        return Ok((re, im));
    }
    let (am, bm) = (mag_lower_ball(c), mag_lower_ball(d));
    if am.is_zero() && bm.is_zero() {
        return all();
    }
    let (x, y) = (c.rad(), d.rad());
    // t = ((|a|-x)^2 + (|b|-y)^2)(a^2 + b^2), lower bound
    let t = mag_add_lower(&am.mul_lower(&am), &bm.mul_lower(&bm))?
        .mul_lower(&Mag::from_exact_num_lower(&a2b2));
    let (au, bu) = (mag_of(a), mag_of(b));
    // v = x^2 + y^2
    let v = x.mul_up(x).add_up(&y.mul_up(y));
    // u = |a| ((x^2 + y^2) + 2 |b| y)
    let u = mag_double(&bu)?.mul_up(y).add_up(&v).mul_up(&au);
    // v = |b| ((x^2 + y^2) + 2 |a| x)
    let v = v.add_up(&mag_double(&au)?.mul_up(x)).mul_up(&bu);
    // w = |b^2 - a^2| (upper bound)
    let w = if a.abs_cmp(b).is_some_and(|s| s >= 0) {
        au.mul_up(&au)
    } else {
        bu.mul_up(&bu)
    };
    let u = u.add_up(&w.mul_up(x));
    let v = v.add_up(&w.mul_up(y));
    let re = div_rounded_den(a, &a2b2, inexact, &u.div_up(&t), prec)?;
    let im = div_rounded_den(b, &a2b2, inexact, &v.div_up(&t), prec)?.negated();
    Ok((re, im))
}

/// Arb `_arb_arf_div_rounded_den(_add_err)`: `x / den` rounded to `prec`, plus `prop`. A rounded `den` adds Arb's ulp at
/// `prec − 1`, which covers the quotient rounding and `den`'s relative error below `2^-(prec+2)`; else ε.
fn div_rounded_den(
    x: &ExactNum,
    den: &ExactNum,
    den_inexact: bool,
    prop: &Mag,
    prec: Precision,
) -> BoundResult<Ball> {
    if prop.is_infinity() {
        return Ball::all_reals(prec);
    }
    let Some(q) = classify("cball_div", quot(x, den, prec)?)? else {
        return Ball::all_reals(prec);
    };
    let err = if den_inexact && !q.is_zero() {
        match q.exponent() {
            Some(e) if e > EXPONENT_MIN => prec
                .checked_sub(1)
                .map_or(Mag::ulp_mag(e, 0), |q| Mag::ulp_mag(e, q)),
            Some(_) | None => Mag::ulp_mag(EXPONENT_MIN, 0),
        }
    } else {
        eps(&q, prec, false)
    };
    Ball::from_mid_rad(q, prop.add_up(&err), prec)
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

fn mag_again(r: &Mag) -> BoundResult<Mag> {
    if r.is_infinity() {
        return Ok(Mag::infinity());
    }
    if r.is_zero() {
        return Ok(Mag::zero());
    }
    Ok(Mag::from_exact_num(&r.to_exact_num()?))
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
