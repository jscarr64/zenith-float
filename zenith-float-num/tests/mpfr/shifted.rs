//! Exact MPFR references for operands and results anywhere in the 62-bit exponent range, subnormals included.
//!
//! MPFR's `emin_min` is the zenith `EXPONENT_MIN`, so MPFR cannot hold the zenith subnormal grid
//! (`2^(EXPONENT_MIN - p)`) nor an operand far below it, and an unshifted reference answers 0 or Inf there.
//! Each reference here works in a space scaled by a power of two in which the scaling identity is exact:
//! a sum, product or quotient is rounded at `value * 2^shift`, a `n`-th root of `x` reads `x * 2^(n u)` and
//! gives `root * 2^u`, `hypot` and `atan2` read both operands under one common shift (`hypot` is scaled by it,
//! `atan2` is not). The zenith range is emulated in that space by [`mpfr_reference_shifted`], and the zenith
//! result is compared after the same scaling ([`check`]), so no step rounds or clamps.
//!
//! Where no shift holds both operands (they are further apart than MPFR's range), the smaller one cannot
//! change the rounding and is replaced by a sticky bit of its sign at a documented place (see [`ref_sum`],
//! [`ref_hypot`], [`ref_atan2`]); results far outside the zenith range are replaced by a stand-in far outside
//! it too ([`round_scaled`]). `self_test_*` below checks every reference against plain MPFR where MPFR alone
//! is exact, and against hand-computed results at the subnormal floor and beyond both ends.

use crate::mpfr::common::{
    exact_mpfr_shifted, mpfr_reference_shifted, mpfr_reference_with, show_mp,
};
use gmp_mpfr_sys::mpfr::{self, rnd_t};
use rug::float::Round;
use rug::{Float, Integer};
use std::cmp::Ordering;
use zenith_float_num::{ExactNum, RoundingMode, EXPONENT_MAX, EXPONENT_MIN};

type Res<T> = Result<T, Box<dyn std::error::Error>>;

const MIN: i128 = EXPONENT_MIN as i128;
const MAX: i128 = EXPONENT_MAX as i128;

/// Terms further apart than this (in binary places, top bit to top bit) are not summed exactly: the smaller one
/// is a sticky bit. It exceeds every precision the tests use (20 000) plus the operand widths.
const GAP: i128 = 1 << 16;

/// `m * 2^k` exactly (`m` signed, any size; `k` beyond the zenith range is allowed).
#[derive(Debug)]
pub struct Term {
    pub m: Integer,
    pub k: i128,
}

impl Term {
    pub fn new(m: u64, k: i64, negative: bool) -> Term {
        let m = Integer::from(m);
        Term {
            m: if negative { -m } else { m },
            k: i128::from(k),
        }
    }

    /// The value lies in `[2^(top - 1), 2^top)` in magnitude.
    pub fn top(&self) -> i128 {
        self.k + i128::from(self.m.significant_bits())
    }

    pub fn mul(&self, o: &Term) -> Term {
        Term {
            m: Integer::from(&self.m * &o.m),
            k: self.k + o.k,
        }
    }

    pub fn neg(&self) -> Term {
        Term {
            m: Integer::from(-&self.m),
            k: self.k,
        }
    }

    pub fn ldexp(&self, d: i128) -> Term {
        Term {
            m: Integer::from(&self.m),
            k: self.k + d,
        }
    }

    pub fn abs(&self) -> Term {
        Term {
            m: Integer::from(self.m.abs_ref()),
            k: self.k,
        }
    }

    pub fn is_negative(&self) -> bool {
        self.m < 0
    }

    pub fn is_zero(&self) -> bool {
        self.m == 0
    }
}

/// MPFR's widest exponent range until dropped.
struct Widest(mpfr::exp_t, mpfr::exp_t);

fn widest() -> Widest {
    let w = unsafe { Widest(mpfr::get_emin(), mpfr::get_emax()) };
    unsafe {
        mpfr::set_emin(mpfr::get_emin_min());
        mpfr::set_emax(mpfr::get_emax_max());
    }
    w
}

impl Drop for Widest {
    fn drop(&mut self) {
        unsafe {
            mpfr::set_emin(self.0);
            mpfr::set_emax(self.1);
        }
    }
}

/// `t * 2^r` as an exact MPFR number; its top bit must be well inside MPFR's range.
fn exact(t: &Term, r: i128) -> Res<Float> {
    let top = t.top() + r;
    if !t.is_zero() && top.abs() > 1 << 61 {
        return Err(format!("exact: 2^{top} is outside MPFR's range").into());
    }
    let _w = widest();
    let bits = t.m.significant_bits().max(1);
    let mut f = Float::with_val(bits, &t.m);
    let sh = i64::try_from(t.k + r)?;
    unsafe { mpfr::mul_2si(f.as_raw_mut(), f.as_raw(), sh, rnd_t::RNDN) };
    Ok(f)
}

/// The reference: `v` is the correctly rounded zenith result times `2^shift`.
#[derive(Debug)]
pub struct Want {
    pub v: Float,
    pub inexact: bool,
    pub shift: i64,
}

fn zero(negative: bool) -> Want {
    let z = Float::with_val(2, 0);
    Want {
        v: if negative { -z } else { z },
        inexact: false,
        shift: 0,
    }
}

fn nan() -> Want {
    Want {
        v: Float::with_val(2, rug::float::Special::Nan),
        inexact: false,
        shift: 0,
    }
}

/// Rounds `eval` (the result times `2^r`) to `p` bits in the zenith range scaled by `2^r`. `e_est` is the
/// result's top exponent to within two places: a result more than four places above `EXPONENT_MAX` (or below
/// the subnormal grid `EXPONENT_MIN - p`) overflows (underflows) like any other one that far out, and is
/// replaced by `1.5` of the sign in a space where it is that far out (eight places), since `r` would leave
/// MPFR's range.
pub fn round_scaled(
    p: u32,
    rm: RoundingMode,
    e_est: i128,
    negative: bool,
    r: i128,
    eval: impl Fn(&mut Float, rnd_t) -> i32,
) -> Res<Want> {
    let pi = i128::from(p);
    let stand_in_shift = if e_est > MAX + 4 {
        Some(-(MAX + 8))
    } else if e_est < MIN - pi - 4 {
        Some(-(MIN - pi - 8))
    } else {
        None
    };
    if let Some(s) = stand_in_shift {
        // `ToOdd` is nearest with ties to odd: a stand-in is never a tie, so it rounds as `ToEven`.
        let rm = if rm == RoundingMode::ToOdd { RoundingMode::ToEven } else { rm };
        let v = exact(&Term::new(3, -1, negative), 0)?;
        let shift = i64::try_from(s)?;
        if s > 0 {
            // Far below the grid `g = 2^8` of this space: exact math gives `0` or `g` of the sign, inexact
            // (no MPFR call: the answer is fixed by the rounding direction alone).
            let away = match rm {
                RoundingMode::Up => !negative,
                RoundingMode::Down => negative,
                RoundingMode::FromZero => true,
                _ => false,
            };
            let mut v = Float::new(p);
            if away {
                v = Float::with_val(p, exact(&Term::new(1, 8, negative), 0)?);
            } else if negative {
                v = -v;
            }
            return Ok(Want {
                v,
                inexact: true,
                shift,
            });
        }
        let (v, inexact) = mpfr_reference_shifted(p, rm, shift, |o, rnd| set_checked(o, &v, rnd))?;
        return Ok(Want { v, inexact, shift });
    }
    let shift = i64::try_from(r)?;
    let (v, inexact) = mpfr_reference_shifted(p, rm, shift, eval)?;
    Ok(Want { v, inexact, shift })
}

/// `o = v` rounded in the current range. `mpfr_set` assumes `v` is inside that range (a shifted `v` near `1`
/// can sit above `emax` or below `emin`), so `mpfr_check_range` applies the overflow or underflow.
fn set_checked(o: &mut Float, v: &Float, rnd: rnd_t) -> i32 {
    unsafe {
        let t = mpfr::set(o.as_raw_mut(), v.as_raw(), rnd);
        mpfr::check_range(o.as_raw_mut(), t, rnd)
    }
}

/// `f` plus a sticky bit of the sign `negative`, far below `p` and below `f`'s own last bit.
fn with_sticky(f: &Float, p: u32, negative: bool) -> Res<Float> {
    let _w = widest();
    let prec = f.prec().max(p) + 16;
    let e = unsafe { mpfr::get_exp(f.as_raw()) };
    let tiny = exact(
        &Term::new(1, 0, negative),
        i128::from(e) - i128::from(prec) - 4,
    )?;
    let mut g = Float::with_val(prec + 8, f);
    g += &tiny;
    Ok(g)
}

/// The sum of exact terms. Terms are grouped so that within a group consecutive top bits are at most [`GAP`]
/// apart; the first group with a nonzero sum is the value, rounded at `2^-top` (so it sits at one), and the
/// next nonzero group only decides on which side of it the exact sum lies: it is at least `GAP - 192` places
/// below the value's lowest bit, so a sticky bit of its sign below `p` and below that bit rounds the same.
/// An all-zero sum is `+0` (`-0` in `Down`), as an exact cancellation of nonzero terms.
pub fn ref_sum(terms: &[Term], p: u32, rm: RoundingMode) -> Res<Want> {
    let mut ts: Vec<&Term> = terms.iter().filter(|t| !t.is_zero()).collect();
    ts.sort_by_key(|t| std::cmp::Reverse(t.top()));
    let mut groups: Vec<Vec<&Term>> = Vec::new();
    for t in ts {
        match groups.last_mut() {
            Some(g) if g.last().is_some_and(|l| l.top() - t.top() <= GAP) => g.push(t),
            _ => groups.push(vec![t]),
        }
    }
    let mut sums: Vec<Term> = Vec::new();
    for g in &groups {
        let kmin = g.iter().map(|t| t.k).min().ok_or("empty group")?;
        let mut s = Integer::new();
        for t in g {
            let sh = u32::try_from(t.k - kmin)?;
            s += Integer::from(&t.m << sh);
        }
        if s != 0 {
            sums.push(Term { m: s, k: kmin });
        }
    }
    let Some(main) = sums.first() else {
        return Ok(zero(rm == RoundingMode::Down));
    };
    let e = main.top();
    let f = exact(main, -e)?;
    let f = match sums.get(1) {
        Some(next) => with_sticky(&f, p, next.is_negative())?,
        None => f,
    };
    round_scaled(p, rm, e, main.is_negative(), -e, |o, rnd| {
        set_checked(o, &f, rnd)
    })
}

/// `a / b`: both operands scaled to `[1/2, 1)`, the quotient by the difference of the two scalings.
pub fn ref_div(a: &Term, b: &Term, p: u32, rm: RoundingMode) -> Res<Want> {
    if b.is_zero() {
        return Err("ref_div: zero divisor".into());
    }
    if a.is_zero() {
        return Ok(zero(a.is_negative() != b.is_negative()));
    }
    let (ta, tb) = (a.top(), b.top());
    let (fa, fb) = (exact(a, -ta)?, exact(b, -tb)?);
    round_scaled(
        p,
        rm,
        ta - tb,
        a.is_negative() != b.is_negative(),
        tb - ta,
        |o, rnd| unsafe { mpfr::div(o.as_raw_mut(), fa.as_raw(), fb.as_raw(), rnd) },
    )
}

pub fn ref_recip(x: &Term, p: u32, rm: RoundingMode) -> Res<Want> {
    ref_div(&Term::new(1, 0, false), x, p, rm)
}

/// The `n`-th root of `x`: `x * 2^(n u)` (exact, near one) has the root `root * 2^u`. An even root of a negative
/// number is NaN.
pub fn ref_root(x: &Term, n: u32, p: u32, rm: RoundingMode) -> Res<Want> {
    if x.is_zero() {
        return Ok(zero(x.is_negative()));
    }
    if x.is_negative() && n.is_multiple_of(2) {
        return Ok(nan());
    }
    let ni = i128::from(n);
    let u = (-x.top()).div_euclid(ni);
    let fx = exact(x, ni * u)?;
    round_scaled(p, rm, -u, x.is_negative(), u, |o, rnd| unsafe {
        match n {
            2 => mpfr::sqrt(o.as_raw_mut(), fx.as_raw(), rnd),
            3 => mpfr::cbrt(o.as_raw_mut(), fx.as_raw(), rnd),
            _ => mpfr::rootn_ui(o.as_raw_mut(), fx.as_raw(), n.into(), rnd),
        }
    })
}

/// `hypot(a, b)` under the shift that puts the larger operand at one; operands more than [`GAP`] places apart
/// give `|big| (1 + small^2 / (2 big^2))`, the larger magnitude with a positive sticky bit ([`ref_sum`]).
pub fn ref_hypot(a: &Term, b: &Term, p: u32, rm: RoundingMode) -> Res<Want> {
    let (big, small) = if a.top() >= b.top() { (a, b) } else { (b, a) };
    if small.is_zero() || big.top() - small.top() > GAP {
        return ref_sum(&[big.abs(), small.abs()], p, rm);
    }
    let r = -big.top();
    let (fa, fb) = (exact(a, r)?, exact(b, r)?);
    round_scaled(p, rm, big.top(), false, r, |o, rnd| unsafe {
        mpfr::hypot(o.as_raw_mut(), fa.as_raw(), fb.as_raw(), rnd)
    })
}

/// `f` truncated with an inexact flag `t` made round-to-odd (the last bit set when inexact): rounding the result
/// to fewer than `prec - 1` bits is rounding the exact value.
fn odd(mut f: Float, t: i32) -> Float {
    if t != 0 && unsafe { mpfr::min_prec(f.as_raw()) } < f.prec() as mpfr::prec_t {
        unsafe {
            if f.is_sign_negative() {
                mpfr::nextbelow(f.as_raw_mut());
            } else {
                mpfr::nextabove(f.as_raw_mut());
            }
        }
    }
    f
}

/// `pi * 2^h` of the sign `negative`, round-to-odd at `p + 200` bits (`pi` is not within `2^-(p + 200)` of a
/// `p`-bit boundary at the precisions tested, so a perturbation far smaller rounds the same).
fn pi_odd(p: u32, h: i64, negative: bool) -> Float {
    let _w = widest();
    let mut f = Float::new(p + 200);
    let t = unsafe { mpfr::const_pi(f.as_raw_mut(), rnd_t::RNDZ) };
    let mut f = odd(f, t);
    unsafe { mpfr::mul_2si(f.as_raw_mut(), f.as_raw(), h, rnd_t::RNDN) };
    if negative {
        -f
    } else {
        f
    }
}

/// `atan2(y, x)` (nonzero operands). With `d = top(y) - top(x)` and `T = p + 600`:
/// - `|d| <= T`: MPFR's `atan2` of both operands under the shift that puts the larger at one (the result is not
///   scaled, and is above `2^-(T + 2)`, far inside MPFR's range);
/// - `d < -T`, `x > 0`: the quotient `q = y / x` less `q^3 / 3` (relative `2^-(2T)`): `q` scaled near one,
///   computed round-to-odd at `p + 700` bits and, when exact, moved one place toward zero, then rounded in the
///   zenith range (subnormal and underflowing quotients included);
/// - `d < -T`, `x < 0`: `+-(pi - |q|)`, and `d > T`: `+-(pi/2 -+ |x / y|)`: `pi` round-to-odd ([`pi_odd`]).
pub fn ref_atan2(y: &Term, x: &Term, p: u32, rm: RoundingMode) -> Res<Want> {
    if y.is_zero() || x.is_zero() {
        return Err("ref_atan2: zero operand".into());
    }
    let (ty, tx) = (y.top(), x.top());
    let d = ty - tx;
    let t = i128::from(p) + 600;
    let neg = y.is_negative();
    if d < -t && !x.is_negative() {
        let (fy, fx) = (exact(y, -ty)?, exact(x, -tx)?);
        let q = {
            let _w = widest();
            let mut q = Float::new(p + 700);
            let tq = unsafe { mpfr::div(q.as_raw_mut(), fy.as_raw(), fx.as_raw(), rnd_t::RNDZ) };
            if tq == 0 {
                unsafe {
                    if neg {
                        mpfr::nextabove(q.as_raw_mut());
                    } else {
                        mpfr::nextbelow(q.as_raw_mut());
                    }
                }
                q
            } else {
                odd(q, tq)
            }
        };
        return round_scaled(p, rm, d, neg, tx - ty, |o, rnd| set_checked(o, &q, rnd));
    }
    if d < -t || d > t {
        let h = if d > t { -1 } else { 0 };
        let f = pi_odd(p, h, neg);
        return round_scaled(p, rm, 2, neg, 0, |o, rnd| set_checked(o, &f, rnd));
    }
    let r = -ty.max(tx);
    let (fy, fx) = (exact(y, r)?, exact(x, r)?);
    round_scaled(p, rm, 2, neg, 0, |o, rnd| unsafe {
        mpfr::atan2(o.as_raw_mut(), fy.as_raw(), fx.as_raw(), rnd)
    })
}

/// `ln`, `log2` or `log10` of a positive `x`; `mp_fn` is MPFR's correctly rounded function. An `x` whose top
/// bit is inside [`exact`]'s range is read as is. Further out (the zenith subnormal grid and the top binades),
/// `log(f 2^t) = log(f) + t log(2)` with `f = x 2^-t` in `[1/2, 1)` is enclosed at growing precision until
/// the enclosure fixes one round-to-odd value at `p + 64` bits, which rounds correctly in every mode. An
/// enclosure of width zero is an exact result (`log2` of a power of two).
pub fn ref_log(
    x: &Term,
    p: u32,
    rm: RoundingMode,
    mp_fn: &impl Fn(&mut Float, &Float, Round) -> Ordering,
) -> Res<Want> {
    if x.is_zero() || x.is_negative() {
        return Err("ref_log: operand must be positive".into());
    }
    let t = x.top();
    if t.abs() <= 1 << 61 {
        let fx = exact(x, 0)?;
        return round_scaled(p, rm, 2, false, 0, |o, rnd| {
            mp_fn(o, &fx, round_of_rnd(rnd)) as i32
        });
    }
    let q = {
        let _w = widest();
        let f = exact(x, -t)?;
        let two = Float::with_val(2, 2);
        let ti = Integer::from(i64::try_from(t)?);
        let mut w = p + 128;
        loop {
            let mut lf = [Float::new(w), Float::new(w)];
            let mut l2 = [Float::new(w), Float::new(w)];
            for (i, r) in [Round::Down, Round::Up].into_iter().enumerate() {
                mp_fn(&mut lf[i], &f, r);
                mp_fn(&mut l2[i], &two, r);
            }
            let (l2_for_lo, l2_for_hi) = if t > 0 { (&l2[0], &l2[1]) } else { (&l2[1], &l2[0]) };
            let tl_lo = Float::with_val_round(w, l2_for_lo * &ti, Round::Down).0;
            let tl_hi = Float::with_val_round(w, l2_for_hi * &ti, Round::Up).0;
            let lo = Float::with_val_round(w, &lf[0] + &tl_lo, Round::Down).0;
            let hi = Float::with_val_round(w, &lf[1] + &tl_hi, Round::Up).0;
            if lo == hi {
                break lo;
            }
            let (zl, ol) = Float::with_val_round(p + 64, &lo, Round::Zero);
            let (zh, oh) = Float::with_val_round(p + 64, &hi, Round::Zero);
            let same_side = lo.is_sign_negative() == hi.is_sign_negative();
            if same_side && zl == zh && ol != Ordering::Equal && oh != Ordering::Equal {
                break odd(zl, 1);
            }
            if w > 1 << 20 {
                return Err(format!("ref_log: no enclosure at {w} bits").into());
            }
            w *= 2;
        }
    };
    round_scaled(p, rm, 2, q.is_sign_negative(), 0, |o, rnd| {
        set_checked(o, &q, rnd)
    })
}

fn round_of_rnd(rnd: rnd_t) -> Round {
    match rnd {
        rnd_t::RNDZ => Round::Zero,
        rnd_t::RNDU => Round::Up,
        rnd_t::RNDD => Round::Down,
        rnd_t::RNDA => Round::AwayZero,
        _ => Round::Nearest,
    }
}

fn show_zf(x: &ExactNum) -> String {
    match x.as_raw_parts() {
        Some((w, n, s, e, inx)) => format!("{s:?} {w:x?} n={n} e={e} inexact={inx}"),
        None => format!("special nan={} inf={}", x.is_nan(), x.is_inf()),
    }
}

/// Asserts that `got` is the reference: NaN, the signed infinity, the signed zero, or the value scaled by
/// `2^want.shift`; with `flags`, the inexact flag too (not for NaN and infinities).
pub fn check(got: &ExactNum, want: &Want, flags: bool, info: &str) -> Res<()> {
    if want.v.is_nan() {
        assert!(got.is_nan(), "{info}: MPFR NaN, zenith {}", show_zf(got));
        return Ok(());
    }
    if want.v.is_infinite() {
        assert!(
            got.is_inf(),
            "{info}: MPFR {}, zenith {}",
            show_mp(&want.v),
            show_zf(got)
        );
        assert_eq!(
            got.is_negative(),
            want.v.is_sign_negative(),
            "{info}: sign of Inf"
        );
        return Ok(());
    }
    assert!(
        !got.is_nan() && !got.is_inf(),
        "{info}: MPFR {} (shift {}), zenith {}",
        show_mp(&want.v),
        want.shift,
        show_zf(got)
    );
    if want.v.is_zero() {
        assert!(got.is_zero(), "{info}: MPFR 0, zenith {}", show_zf(got));
        assert_eq!(
            got.is_negative(),
            want.v.is_sign_negative(),
            "{info}: sign of zero"
        );
    } else {
        let _w = widest();
        let g = exact_mpfr_shifted(got, want.shift)?;
        assert!(
            g == want.v,
            "{info}: MPFR {} (shift {}), zenith {} = {}",
            show_mp(&want.v),
            want.shift,
            show_zf(got),
            show_mp(&g)
        );
    }
    if flags {
        assert_eq!(
            got.inexact(),
            want.inexact,
            "{info}: inexact flag, MPFR {} (shift {})",
            show_mp(&want.v),
            want.shift
        );
    }
    Ok(())
}

// ---- self-tests: the references against plain MPFR, and at the ends against hand-computed values ----------------

#[cfg(test)]
mod self_test {
    use super::*;

    const MODES: [RoundingMode; 6] = [
        RoundingMode::ToEven,
        RoundingMode::ToOdd,
        RoundingMode::Up,
        RoundingMode::Down,
        RoundingMode::ToZero,
        RoundingMode::FromZero,
    ];
    const PRECISIONS: [u32; 8] = [1, 2, 3, 24, 64, 65, 128, 200];

    /// `want` (scaled by `2^shift`) equals the plain result `plain` (unscaled), flags included.
    fn agree(want: &Want, plain: &(Float, bool), info: &str) {
        let _w = widest();
        let mut back = Float::with_val(want.v.prec(), &want.v);
        unsafe { mpfr::mul_2si(back.as_raw_mut(), back.as_raw(), -want.shift, rnd_t::RNDN) };
        assert!(
            (back.is_nan() && plain.0.is_nan()) || back == plain.0,
            "{info}: shifted {} (shift {}), plain {}",
            show_mp(&want.v),
            want.shift,
            show_mp(&plain.0)
        );
        if back.is_zero() {
            assert_eq!(
                back.is_sign_negative(),
                plain.0.is_sign_negative(),
                "{info}: zero sign"
            );
        }
        assert_eq!(want.inexact, plain.1, "{info}: inexact");
    }

    /// Deterministic operands with ordinary exponents (MPFR alone is exact there).
    fn operands() -> Vec<(u64, i64, bool)> {
        let mut v = Vec::new();
        let mut s = 0x9e37_79b9_7f4a_7c15u64;
        for i in 0..24 {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            let m = match i % 4 {
                0 => 1,
                1 => 3,
                2 => s,
                _ => s | 1,
            };
            let k = (s % 400) as i64 - 200;
            v.push((m, k, i % 3 == 0));
        }
        v
    }

    fn plain_f(m: u64, k: i64, neg: bool) -> Float {
        let f = Float::with_val(64, m) << (k as i32);
        if neg {
            -f
        } else {
            f
        }
    }

    #[test]
    fn shifted_references_match_plain_mpfr_inside_its_range() -> Res<()> {
        let ops = operands();
        let mut n = 0usize;
        for (i, &(m1, k1, n1)) in ops.iter().enumerate() {
            let (a, fa) = (Term::new(m1, k1, n1), plain_f(m1, k1, n1));
            let (m2, k2, n2) = ops[(i * 7 + 3) % ops.len()];
            let (b, fb) = (Term::new(m2, k2, n2), plain_f(m2, k2, n2));
            let (m3, k3, n3) = ops[(i * 5 + 1) % ops.len()];
            let (c, fc) = (Term::new(m3, k3, n3), plain_f(m3, k3, n3));
            let prod = Float::with_val(128, &fa * &fb);
            for p in PRECISIONS {
                for rm in MODES {
                    let info = format!("{a:?} {b:?} {c:?} p={p} {rm:?}");
                    let plain =
                        |f: &dyn Fn(&mut Float, rnd_t) -> i32| mpfr_reference_with(p, rm, f);
                    agree(
                        &ref_sum(&[a.mul(&b), c.ldexp(0)], p, rm)?,
                        &plain(&|o, rnd| unsafe {
                            mpfr::fma(o.as_raw_mut(), fa.as_raw(), fb.as_raw(), fc.as_raw(), rnd)
                        })?,
                        &format!("fma {info}"),
                    );
                    agree(
                        &ref_sum(&[a.ldexp(0), b.ldexp(0)], p, rm)?,
                        &plain(&|o, rnd| unsafe {
                            mpfr::add(o.as_raw_mut(), fa.as_raw(), fb.as_raw(), rnd)
                        })?,
                        &format!("add {info}"),
                    );
                    agree(
                        &ref_sum(&[a.mul(&b), a.mul(&b).neg()], p, rm)?,
                        &plain(&|o, rnd| unsafe {
                            mpfr::sub(o.as_raw_mut(), prod.as_raw(), prod.as_raw(), rnd)
                        })?,
                        &format!("cancel {info}"),
                    );
                    agree(
                        &ref_div(&a, &b, p, rm)?,
                        &plain(&|o, rnd| unsafe {
                            mpfr::div(o.as_raw_mut(), fa.as_raw(), fb.as_raw(), rnd)
                        })?,
                        &format!("div {info}"),
                    );
                    agree(
                        &ref_recip(&a, p, rm)?,
                        &plain(&|o, rnd| unsafe {
                            mpfr::ui_div(o.as_raw_mut(), 1, fa.as_raw(), rnd)
                        })?,
                        &format!("recip {info}"),
                    );
                    for nr in [2u32, 3, 5, 7, 33] {
                        agree(
                            &ref_root(&a, nr, p, rm)?,
                            &plain(&|o, rnd| unsafe {
                                mpfr::rootn_ui(o.as_raw_mut(), fa.as_raw(), nr.into(), rnd)
                            })?,
                            &format!("root {nr} {info}"),
                        );
                    }
                    agree(
                        &ref_hypot(&a, &b, p, rm)?,
                        &plain(&|o, rnd| unsafe {
                            mpfr::hypot(o.as_raw_mut(), fa.as_raw(), fb.as_raw(), rnd)
                        })?,
                        &format!("hypot {info}"),
                    );
                    agree(
                        &ref_atan2(&a, &b, p, rm)?,
                        &plain(&|o, rnd| unsafe {
                            mpfr::atan2(o.as_raw_mut(), fa.as_raw(), fb.as_raw(), rnd)
                        })?,
                        &format!("atan2 {info}"),
                    );
                    n += 1;
                }
            }
        }
        assert!(n >= 1000, "{n}");
        Ok(())
    }

    /// The sticky and far-apart branches against plain MPFR where MPFR can still hold both terms (a gap above
    /// `GAP` and `p + 600`, inside MPFR's range).
    #[test]
    fn far_apart_branches_match_plain_mpfr() -> Res<()> {
        for gap in [GAP as i64 + 1, 100_000, 1_000_000] {
            for (m1, m2) in [(1u64, 1u64), (3, 0xdead_beef_1234_5677), (u64::MAX, 1)] {
                for (n1, n2) in [(false, false), (false, true), (true, false), (true, true)] {
                    let (a, b) = (Term::new(m1, 0, n1), Term::new(m2, -gap, n2));
                    let fa = plain_f(m1, 0, n1);
                    let fb = {
                        let _w = widest();
                        let f = Float::with_val(64, m2) >> (gap as i32);
                        if n2 {
                            -f
                        } else {
                            f
                        }
                    };
                    for p in PRECISIONS {
                        for rm in MODES {
                            let info = format!("{a:?} {b:?} p={p} {rm:?}");
                            let plain = |f: &dyn Fn(&mut Float, rnd_t) -> i32| {
                                mpfr_reference_with(p, rm, f)
                            };
                            agree(
                                &ref_sum(&[a.ldexp(0), b.ldexp(0)], p, rm)?,
                                &plain(&|o, rnd| unsafe {
                                    mpfr::add(o.as_raw_mut(), fa.as_raw(), fb.as_raw(), rnd)
                                })?,
                                &format!("add {info}"),
                            );
                            agree(
                                &ref_hypot(&a, &b, p, rm)?,
                                &plain(&|o, rnd| unsafe {
                                    mpfr::hypot(o.as_raw_mut(), fa.as_raw(), fb.as_raw(), rnd)
                                })?,
                                &format!("hypot {info}"),
                            );
                            for (y, x, fy, fx) in [(&b, &a, &fb, &fa), (&a, &b, &fa, &fb)] {
                                agree(
                                    &ref_atan2(y, x, p, rm)?,
                                    &plain(&|o, rnd| unsafe {
                                        mpfr::atan2(o.as_raw_mut(), fy.as_raw(), fx.as_raw(), rnd)
                                    })?,
                                    &format!("atan2 {y:?} {x:?} p={p} {rm:?}"),
                                );
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }

    fn is(want: &Want, t: &Term, inexact: bool, info: &str) -> Res<()> {
        let e = exact(t, i128::from(want.shift))?;
        let _w = widest();
        assert!(
            want.v == e,
            "{info}: {} (shift {}), expected {t:?}",
            show_mp(&want.v),
            want.shift
        );
        assert_eq!(want.inexact, inexact, "{info}: inexact");
        Ok(())
    }

    /// Hand-computed results at the subnormal grid `2^(EXPONENT_MIN - p)` and beyond both ends of the range.
    #[test]
    fn ends_of_the_range_are_hand_computed() -> Res<()> {
        let min = EXPONENT_MIN;
        let max = EXPONENT_MAX;
        for p in [1u32, 24, 64, 128] {
            let pi = i64::from(p);
            let g = Term::new(1, min - pi, false);
            // 2^(MIN - p - 36): Up gives the grid, ToEven and ToZero zero
            let tiny = [Term::new(1, min - pi - 36, false)];
            is(&ref_sum(&tiny, p, RoundingMode::Up)?, &g, true, "up tiny")?;
            let w = ref_sum(&tiny, p, RoundingMode::ToEven)?;
            assert!(
                w.v.is_zero() && !w.v.is_sign_negative() && w.inexact,
                "even tiny"
            );
            // 1.5 g is a tie between g (odd) and 2 g (even)
            let tie = [Term::new(3, min - pi - 1, false)];
            is(
                &ref_sum(&tie, p, RoundingMode::ToEven)?,
                &Term::new(2, min - pi, false),
                true,
                "tie even",
            )?;
            is(&ref_sum(&tie, p, RoundingMode::ToOdd)?, &g, true, "tie odd")?;
            // far beyond the bottom (stand-in): FromZero gives -g
            let deep = [Term::new(1, min - 5_000_000, true)];
            is(
                &ref_sum(&deep, p, RoundingMode::FromZero)?,
                &g.neg(),
                true,
                "deep",
            )?;
            // the quotient 2^(MAX - 10) / 2^(MIN + 10) is far beyond the top: Inf, or the largest finite in ToZero
            let (a, b) = (Term::new(1, max - 10, false), Term::new(1, min + 10, false));
            let w = ref_div(&a, &b, p, RoundingMode::ToEven)?;
            assert!(w.v.is_infinite() && !w.v.is_sign_negative(), "overflow");
            let big = Term {
                m: (Integer::from(1) << p) - 1u32,
                k: MAX - i128::from(p),
            };
            is(
                &ref_div(&a, &b, p, RoundingMode::ToZero)?,
                &big,
                true,
                "largest",
            )?;
            // sqrt of the grid at an even exponent is exact
            let sq = Term::new(1, 2 * ((min - pi) / 2), false);
            is(
                &ref_root(&sq, 2, p, RoundingMode::ToEven)?,
                &Term::new(1, (min - pi) / 2, false),
                false,
                "sqrt",
            )?;
            // 1 / 2^(MIN - 2): 2^(-(MIN - 2)) = 2^(MAX + 3) overflows
            let w = ref_recip(&Term::new(1, min - 2, true), p, RoundingMode::Up)?;
            is(&w, &big.neg(), true, "recip")?;
            // atan2(2^(MIN - p - 1), 1) = 2^(MIN - p - 1) (1 - eps): Down 0, Up the grid, ToEven 0
            let (y, one) = (Term::new(1, min - pi - 1, false), Term::new(1, 0, false));
            let w = ref_atan2(&y, &one, p, RoundingMode::Down)?;
            assert!(w.v.is_zero() && w.inexact, "atan2 down");
            is(
                &ref_atan2(&y, &one, p, RoundingMode::Up)?,
                &g,
                true,
                "atan2 up",
            )?;
            let w = ref_atan2(&y, &one, p, RoundingMode::ToEven)?;
            assert!(w.v.is_zero(), "atan2 even (below the midpoint)");
        }
        Ok(())
    }
}
