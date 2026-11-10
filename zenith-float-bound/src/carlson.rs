//! Slice S4b: Carlson symmetric elliptic integrals (RF, RD, RJ, RC), the polylogarithm Li_s and the Lambert W function (principal and -1 branches).
//!
//! Owner: slice S4b (Cassie). Only that slice edits this file.
//! Part A (CONTRACT 16.S4b.1 to 16.S4b.3): `Ball::lambert_w0` and `Ball::lambert_wm1`. Each endpoint value of W is a
//! Newton estimate turned into a certified bracket by a directed-rounding sign test (16.S4b.2); the Carlson integrals
//! and Li_s are later S4b parts. Contract text: `docs/CONTRACT.md` section "Slice S4b"; harness rows:
//! `tests/harness/slice_carlson.rs`; goldens: `tests/goldens/slices/carlson/`; mutants: `tools/mutants/manifest.d/carlson.txt`.

use core::cmp::Ordering;

use zenith_float_num::{Consts, ExactNum, RoundingMode, EXPONENT_MIN};

use crate::ball::Ball;
use crate::error::{BoundError, BoundResult};
use crate::mag::Mag;
use crate::policy::DomainPolicy;
use crate::precision::{Precision, MAX_PRECISION, MIN_PRECISION};

/// Every Ball midpoint rounds to nearest, ties to even (CONTRACT 15.1, O19).
const MID_RM: RoundingMode = RoundingMode::ToEven;
/// Guard bits of the first working precision over `prec` (16.S4b.2).
const GUARD_BITS: Precision = 64;
/// The certified bracket half-width is `2^(e(w) - prec - BRACKET_SHIFT)` (16.S4b.2).
const BRACKET_SHIFT: i64 = 8;
/// The working-precision cap is `CAP_FACTOR * prec + CAP_EXTRA` (16.S4b.2).
const CAP_FACTOR: Precision = 8;
const CAP_EXTRA: Precision = 4096;
/// Newton steps per working precision (the iteration then is certified; an unconverged estimate only fails the test).
const NEWTON_STEPS: usize = 200;
/// Endpoints `m ± r` are computed exactly when the exponents of `m` and `r` differ by at most this.
const EXACT_GAP: i64 = 1 << 16;
/// `0.1b` as a `Mag`/`from_mag_parts` mantissa.
const TOP: u64 = 0x8000_0000_0000_0000;

/// The two real branches of W.
#[derive(Debug, PartialEq, Eq)]
enum Branch {
    /// W0, values in `[-1, inf)`.
    Principal,
    /// W-1, values in `(-inf, -1]`.
    Lower,
}

/// One certification attempt: the branch, the test form (log or product), `prec` and the working precision.
struct At<'a> {
    branch: &'a Branch,
    log_form: bool,
    prec: Precision,
    wp: Precision,
}

/// Which side of `W(x)` a certified value bounds.
#[derive(Debug, PartialEq, Eq)]
enum Side {
    Below,
    Above,
}

impl Ball {
    /// Principal-branch Lambert W, `W0`, on the ball (CONTRACT 16.S4b.1, 16.S4b.2): domain `[-1/e, inf)`.
    ///
    /// # Errors
    /// `InvalidArgument` when `prec` is outside `MIN_PRECISION..=MAX_PRECISION`; `Domain` per `policy` (16.S4b.2);
    /// `Resource` when a bound cannot be certified within the working-precision cap.
    pub fn lambert_w0(&self, prec: Precision, policy: &DomainPolicy) -> BoundResult<Ball> {
        lambert("Ball::lambert_w0", self, prec, policy, &Branch::Principal)
    }

    /// Branch -1 of Lambert W, `W-1`, on the ball (CONTRACT 16.S4b.1, 16.S4b.2): domain `[-1/e, 0)`.
    ///
    /// # Errors
    /// `InvalidArgument` when `prec` is outside `MIN_PRECISION..=MAX_PRECISION`; `Domain` per `policy` (16.S4b.2);
    /// `Resource` when a bound cannot be certified within the working-precision cap.
    pub fn lambert_wm1(&self, prec: Precision, policy: &DomainPolicy) -> BoundResult<Ball> {
        lambert("Ball::lambert_wm1", self, prec, policy, &Branch::Lower)
    }
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

fn domain(op: &'static str) -> BoundError {
    BoundError::Domain {
        op,
        detail: "input outside the Lambert W domain",
    }
}

fn resource(op: &'static str, detail: &'static str) -> BoundError {
    BoundError::Resource { op, detail }
}

/// `Strict`: `Domain`; `Clip`: all reals (16.S4b.2, an unbounded image or input).
fn unbounded(op: &'static str, prec: Precision, strict: bool) -> BoundResult<Ball> {
    if strict {
        Err(domain(op))
    } else {
        Ball::all_reals(prec)
    }
}

/// `a` compared with `b`; `None` when either is NaN.
fn cmp(a: &ExactNum, b: &ExactNum) -> Option<Ordering> {
    a.compare(b).map(|s| s.cmp(&0))
}

fn num(v: i64) -> ExactNum {
    ExactNum::from(v)
}

/// `x` rounded once to `p` bits with `rm` (`x` itself when it fits).
fn round(x: &ExactNum, p: Precision, rm: RoundingMode) -> ExactNum {
    num(1).mul(x, p, rm)
}

/// `2^k`, with `k` raised to `EXPONENT_MIN` when below it.
fn pow2(k: i64) -> Option<ExactNum> {
    let e = k.max(EXPONENT_MIN).checked_add(1)?;
    ExactNum::from_mag_parts(TOP, e, 64).ok()
}

fn finite(x: &ExactNum) -> bool {
    !x.is_nan() && !x.is_inf()
}

/// The working-precision cap `CAP_FACTOR * prec + CAP_EXTRA` (`Resource` on overflow: unreachable for a valid `prec`).
fn wp_cap(op: &'static str, prec: Precision) -> BoundResult<Precision> {
    prec.checked_mul(CAP_FACTOR)
        .and_then(|p| p.checked_add(CAP_EXTRA))
        .ok_or_else(|| resource(op, "precision cap"))
}

/// The first working precision `prec + GUARD_BITS` (`Resource` on overflow: unreachable for a valid `prec`).
fn wp_start(op: &'static str, prec: Precision) -> BoundResult<Precision> {
    prec.checked_add(GUARD_BITS)
        .ok_or_else(|| resource(op, "working precision"))
}

fn lambert(
    op: &'static str,
    x: &Ball,
    prec: Precision,
    policy: &DomainPolicy,
    branch: &Branch,
) -> BoundResult<Ball> {
    check_prec(op, prec)?;
    let strict = match policy {
        DomainPolicy::Strict => true,
        DomainPolicy::Clip => false,
    };
    let mut cc = Consts::new().map_err(|_| resource(op, "zf constants cache"))?;
    let m = x.mid();
    let r = x.rad();
    if r.is_infinity() {
        return unbounded(op, prec, strict);
    }
    let (lo_own, hi_own) = if r.is_zero() {
        (None, None)
    } else {
        let rn = match r.to_exact_num() {
            Ok(v) => v,
            Err(_) => match r.exp() {
                Some(e) if e > 0 => return unbounded(op, prec, strict),
                Some(_) | None => pow2(EXPONENT_MIN).ok_or_else(|| resource(op, "tiny radius"))?,
            },
        };
        endpoints(op, m, rn, prec)?
    };
    let lo = lo_own.as_ref().unwrap_or(m);
    let hi = hi_own.as_ref().unwrap_or(m);
    if !finite(hi) {
        return unbounded(op, prec, strict);
    }
    let cap = wp_cap(op, prec)?;
    let lo_in = finite(lo) && inside(lo, cap, &mut cc);
    let (l, h) = match branch {
        Branch::Principal => {
            if lo_in {
                let l = w_bound(op, branch, lo, &Side::Below, prec, &mut cc)?;
                (l, w_bound(op, branch, hi, &Side::Above, prec, &mut cc)?)
            } else if strict || !inside(hi, cap, &mut cc) {
                return Err(domain(op));
            } else {
                (
                    num(-1),
                    w_bound(op, branch, hi, &Side::Above, prec, &mut cc)?,
                )
            }
        }
        Branch::Lower => {
            if hi.is_zero() || hi.is_positive() {
                return unbounded(op, prec, strict);
            }
            if lo_in {
                let l = w_bound(op, branch, hi, &Side::Below, prec, &mut cc)?;
                (l, w_bound(op, branch, lo, &Side::Above, prec, &mut cc)?)
            } else if strict || !inside(hi, cap, &mut cc) {
                return Err(domain(op));
            } else {
                (
                    w_bound(op, branch, hi, &Side::Below, prec, &mut cc)?,
                    num(-1),
                )
            }
        }
    };
    result(op, &l, &h, prec)
}

/// `(m - r, m + r)`: exact when the exponents are close, else rounded outward at `prec + 2 * GUARD_BITS` bits.
fn endpoints(
    op: &'static str,
    m: &ExactNum,
    rn: ExactNum,
    prec: Precision,
) -> BoundResult<(Option<ExactNum>, Option<ExactNum>)> {
    let gap = match (m.exponent(), rn.exponent()) {
        (Some(a), Some(b)) => a.checked_sub(b).map(i64::abs),
        (Some(_), None) | (None, Some(_)) | (None, None) => None,
    };
    let close = m.is_zero() || gap.is_some_and(|g| g <= EXACT_GAP);
    if close {
        Ok((Some(m.sub_full_prec(&rn)), Some(m.add_full_prec(&rn))))
    } else {
        let p = GUARD_BITS
            .checked_mul(2)
            .and_then(|g| prec.checked_add(g))
            .ok_or_else(|| resource(op, "endpoint precision"))?;
        let neg = num(0).sub(&rn, p, RoundingMode::Up);
        Ok((
            Some(neg.add(m, p, RoundingMode::Down)),
            Some(rn.add(m, p, RoundingMode::Up)),
        ))
    }
}

/// `x >= -1/e`, decided with `e` rounded `Down` and `Up` at doubling precision (16.S4b.1); undecided at `cap`: false.
fn inside(x: &ExactNum, cap: Precision, cc: &mut Consts) -> bool {
    if !x.is_negative() {
        return true;
    }
    let minus_one = num(-1);
    if cmp(x, &minus_one) != Some(Ordering::Greater) {
        return false;
    }
    let mut wp = GUARD_BITS;
    while wp <= cap {
        let e_hi = cc.e(wp, RoundingMode::Up);
        if cmp(&x.mul_full_prec(&e_hi), &minus_one).is_some_and(|o| o != Ordering::Less) {
            return true;
        }
        let e_lo = cc.e(wp, RoundingMode::Down);
        if cmp(&x.mul_full_prec(&e_lo), &minus_one).is_some_and(|o| o != Ordering::Greater) {
            return false;
        }
        match wp.checked_mul(2) {
            Some(next) => wp = next,
            None => return false,
        }
    }
    false
}

/// A certified bound of `W(x)` on `side` (16.S4b.2), `x` inside the branch's domain.
fn w_bound(
    op: &'static str,
    branch: &Branch,
    x: &ExactNum,
    side: &Side,
    prec: Precision,
    cc: &mut Consts,
) -> BoundResult<ExactNum> {
    if x.is_zero() && *branch == Branch::Principal {
        return Ok(num(0));
    }
    let log_form = *branch == Branch::Lower || cmp(x, &num(4)) == Some(Ordering::Greater);
    let cap = wp_cap(op, prec)?;
    let mut wp = wp_start(op, prec)?;
    while wp <= cap {
        let at = At {
            branch,
            log_form,
            prec,
            wp,
        };
        let w = newton(&at, x, cc);
        let doubled = wp.checked_mul(2);
        // Near the branch point `w*e^w` is flat, so `w` carries only about `wp - 2k` good bits when `|1 + w| ~ 2^-k`:
        // a certified bound from it would be valid but loose. Raise `wp` first while the cap allows.
        if !flat_enough(&w, prec, wp) {
            if let Some(next) = doubled.filter(|n| *n <= cap) {
                wp = next;
                continue;
            }
        }
        if let Some(c) = candidate(&at, &w, x, side, cc) {
            return Ok(c);
        }
        match doubled {
            Some(next) => wp = next,
            None => break,
        }
    }
    Err(resource(
        op,
        "Lambert W bound not certified within the precision cap",
    ))
}

/// `wp >= prec + GUARD_BITS + 2k` where `|1 + w| ~ 2^-k` (true away from `-1`).
fn flat_enough(w: &ExactNum, prec: Precision, wp: Precision) -> bool {
    let k = match w.add_full_prec(&num(1)).exponent() {
        Some(e) if e < 0 => usize::try_from(e.unsigned_abs()).unwrap_or(usize::MAX),
        Some(_) | None => return true,
    };
    k.checked_mul(2)
        .and_then(|d| prec.checked_add(GUARD_BITS)?.checked_add(d))
        .is_some_and(|need| wp >= need)
}

/// The candidate `w -/+ 2^(e(w) - prec - BRACKET_SHIFT)` (rounded outward), clamped at `-1`, when the sign test proves it.
fn candidate(
    at: &At<'_>,
    w: &ExactNum,
    x: &ExactNum,
    side: &Side,
    cc: &mut Consts,
) -> Option<ExactNum> {
    let (branch, prec, wp) = (at.branch, at.prec, at.wp);
    if !finite(w) {
        return None;
    }
    let shift = i64::try_from(prec).ok()?.checked_add(BRACKET_SHIFT)?;
    let delta = pow2(w.exponent()?.checked_sub(shift)?)?;
    let c = match side {
        Side::Below => round(w, wp, RoundingMode::Down).sub(&delta, wp, RoundingMode::Down),
        Side::Above => round(w, wp, RoundingMode::Up).add(&delta, wp, RoundingMode::Up),
    };
    let minus_one = num(-1);
    let vs_one = cmp(&c, &minus_one)?;
    let c = match (branch, side, vs_one) {
        // W0 >= -1 and W-1 <= -1 always: these clamps need no test.
        (Branch::Principal, Side::Below, Ordering::Less)
        | (Branch::Lower, Side::Above, Ordering::Greater | Ordering::Equal) => {
            return Some(minus_one);
        }
        (Branch::Principal, Side::Above, Ordering::Less)
        | (Branch::Lower, Side::Below, Ordering::Greater) => minus_one,
        (Branch::Principal, Side::Below | Side::Above, Ordering::Equal | Ordering::Greater)
        | (Branch::Lower, Side::Below, Ordering::Less | Ordering::Equal)
        | (Branch::Lower, Side::Above, Ordering::Less) => c,
    };
    let ok = if at.log_form {
        log_test(&c, x, side, wp, cc)
    } else {
        product_test(&c, x, side, wp, cc)
    };
    ok.then_some(c)
}

/// `w·e^w <= x` (`Below`) or `>= x` (`Above`) proved with `e^w` rounded the safe way and an exact product; `w >= -1`.
fn product_test(w: &ExactNum, x: &ExactNum, side: &Side, wp: Precision, cc: &mut Consts) -> bool {
    let up = match side {
        Side::Below => !w.is_negative(),
        Side::Above => w.is_negative(),
    };
    let rm = if up {
        RoundingMode::Up
    } else {
        RoundingMode::Down
    };
    let ew = w.exp(wp, rm, cc);
    if !finite(&ew) {
        return false;
    }
    let f = w.mul_full_prec(&ew);
    match (side, cmp(&f, x)) {
        (Side::Below, Some(Ordering::Less | Ordering::Equal))
        | (Side::Above, Some(Ordering::Greater | Ordering::Equal)) => true,
        (Side::Below, Some(Ordering::Greater) | None)
        | (Side::Above, Some(Ordering::Less) | None) => false,
    }
}

/// `ln|w| + w <= ln|x|` (`Below`) or `>=` (`Above`), the logs rounded outward and the sums exact; `w > 0` or `w <= -1`.
fn log_test(w: &ExactNum, x: &ExactNum, side: &Side, wp: Precision, cc: &mut Consts) -> bool {
    if w.is_zero() || (w.is_negative() != x.is_negative()) {
        return false;
    }
    let (rw, rx) = match side {
        Side::Below => (RoundingMode::Up, RoundingMode::Down),
        Side::Above => (RoundingMode::Down, RoundingMode::Up),
    };
    let lw = abs_of(w).ln(wp, rw, cc);
    let lx = abs_of(x).ln(wp, rx, cc);
    if !finite(&lw) || !finite(&lx) {
        return false;
    }
    let lhs = lw.add_full_prec(w);
    match (side, cmp(&lhs, &lx)) {
        (Side::Below, Some(Ordering::Less | Ordering::Equal))
        | (Side::Above, Some(Ordering::Greater | Ordering::Equal)) => true,
        (Side::Below, Some(Ordering::Greater) | None)
        | (Side::Above, Some(Ordering::Less) | None) => false,
    }
}

/// `|x|`, exact.
fn abs_of(x: &ExactNum) -> ExactNum {
    if x.is_negative() {
        num(0).sub_full_prec(x)
    } else {
        num(0).add_full_prec(x)
    }
}

/// `p = sqrt(2 (1 + e x))` at `wp` (the branch-point series variable), zero when rounding makes `1 + e x` negative.
fn branch_p(x: &ExactNum, wp: Precision, cc: &mut Consts) -> ExactNum {
    let ex = cc.e(wp, MID_RM).mul(x, wp, MID_RM);
    let t = ex.add(&num(1), wp, MID_RM);
    if t.is_negative() || !finite(&t) {
        return num(0);
    }
    t.ldexp(1, wp, MID_RM).sqrt(wp, MID_RM)
}

/// A Newton estimate of `W(x)` at `wp` (not trusted: `candidate` certifies it).
fn newton(at: &At<'_>, x: &ExactNum, cc: &mut Consts) -> ExactNum {
    let (branch, log_form, wp) = (at.branch, at.log_form, at.wp);
    // x < -1/4: seed from the branch-point series in p = sqrt(2 (1 + e x)).
    let near_branch =
        pow2(-2).is_some_and(|q| cmp(x, &num(0).sub_full_prec(&q)) == Some(Ordering::Less));
    let minus_one = num(-1);
    let mut w = if near_branch {
        let p = branch_p(x, wp, cc);
        let p2 = p.mul_full_prec(&p).mul(&third(wp), wp, MID_RM);
        match branch {
            Branch::Principal => minus_one.add_full_prec(&p).sub(&p2, wp, MID_RM),
            Branch::Lower => minus_one.sub_full_prec(&p).sub(&p2, wp, MID_RM),
        }
    } else if log_form {
        let big_l = abs_of(x).ln(wp, MID_RM, cc);
        let inner = abs_of(&big_l).ln(wp, MID_RM, cc);
        big_l.sub(&inner, wp, MID_RM)
    } else if x.exponent().is_some_and(|e| e < -4) {
        round(x, wp, MID_RM)
    } else {
        x.add_full_prec(&num(1)).ln(wp, MID_RM, cc)
    };
    let lx = if log_form {
        abs_of(x).ln(wp, MID_RM, cc)
    } else {
        num(0)
    };
    let tol = i64::try_from(wp).ok().and_then(|p| p.checked_sub(4));
    for _ in 0..NEWTON_STEPS {
        if !finite(&w) {
            break;
        }
        let w1 = w.add_full_prec(&num(1));
        let step = if log_form {
            // g = ln|w| + w - ln|x|, w -= g·w/(w+1)
            let g = abs_of(&w)
                .ln(wp, MID_RM, cc)
                .add(&w, wp, MID_RM)
                .sub(&lx, wp, MID_RM);
            g.mul(&w, wp, MID_RM).div(&w1, wp, MID_RM)
        } else {
            // f = w - x e^-w, w -= f/(w+1)
            let e = num(0).sub_full_prec(&w).exp(wp, MID_RM, cc);
            let f = round(&w, wp, MID_RM).sub(&e.mul(x, wp, MID_RM), wp, MID_RM);
            f.div(&w1, wp, MID_RM)
        };
        if !finite(&step) {
            break;
        }
        let next = round(&w, wp, MID_RM).sub(&step, wp, MID_RM);
        let next = keep_branch(branch, log_form, &w, next, wp);
        let small = step.is_zero()
            || match (step.exponent(), w.exponent(), tol) {
                (Some(s), Some(e), Some(t)) => e.checked_sub(s).is_some_and(|d| d > t),
                (Some(_) | None, Some(_) | None, Some(_) | None) => false,
            };
        w = next;
        if small {
            break;
        }
    }
    w
}

/// Keeps an iterate on its branch: halfway back to the old iterate (or towards `-1`) when it jumps across.
fn keep_branch(
    branch: &Branch,
    log_form: bool,
    old: &ExactNum,
    next: ExactNum,
    wp: Precision,
) -> ExactNum {
    let minus_one = num(-1);
    let bad = match branch {
        Branch::Principal if log_form => !next.is_positive(),
        Branch::Principal => cmp(&next, &minus_one) != Some(Ordering::Greater),
        Branch::Lower => cmp(&next, &minus_one) != Some(Ordering::Less),
    };
    if !bad {
        return next;
    }
    let bound = if *branch == Branch::Principal && log_form {
        num(0)
    } else {
        minus_one
    };
    old.add_full_prec(&bound).ldexp(-1, wp, MID_RM)
}

/// `1/3` at `wp`.
fn third(wp: Precision) -> ExactNum {
    num(1).div(&num(3), wp, MID_RM)
}

/// The ball `[l, h]`: `mid = round((l + h) / 2, prec)`, `rad = max(h - mid, mid - l)` rounded up.
fn result(op: &'static str, l: &ExactNum, h: &ExactNum, prec: Precision) -> BoundResult<Ball> {
    let wide = prec
        .checked_add(2)
        .ok_or_else(|| resource(op, "midpoint precision"))?;
    let mid = round(l, wide, MID_RM)
        .add(h, wide, MID_RM)
        .ldexp(-1, prec, MID_RM);
    if !finite(&mid) {
        return Ball::all_reals(prec);
    }
    let above = round(&h.sub_full_prec(&mid), 64, RoundingMode::Up);
    let below = round(&mid.sub_full_prec(l), 64, RoundingMode::Up);
    if above.is_nan() || below.is_nan() {
        return Err(resource(op, "radius"));
    }
    let rad = Mag::from_exact_num(&above).max(Mag::from_exact_num(&below));
    Ball::from_mid_rad(mid, rad, prec)
}

/// In-module checks of private helpers (mutate.sh runs only integration suites, so these back the harness, they do not
/// replace it). Mutant-killers: they pass on the code at 41df22b.
#[cfg(test)]
mod slice_carlson_unit {
    use super::{log_test, num, product_test, result, Consts, RoundingMode, Side};
    use crate::enclosure::Enclosure;
    use crate::tri::Tri;

    /// 614-PRODBELOWGT: `1·e^1 > 1`, so `Below` refuses `w = 1` at `x = 1`; `0·e^0 <= 1` is accepted.
    #[test]
    fn product_test_below_refuses_greater() -> Result<(), String> {
        let mut cc = Consts::new().map_err(|e| format!("Consts::new: {:?}", e))?;
        // no-source: hand-derived: 1·e^1 = e > 1 (refused by Below), 0·e^0 = 0 <= 1 (accepted).
        if product_test(&num(1), &num(1), &Side::Below, 128, &mut cc) {
            return Err("product_test(1, 1, Below) accepted 1·e > 1".into());
        }
        if !product_test(&num(0), &num(1), &Side::Below, 128, &mut cc) {
            return Err("product_test(0, 1, Below) refused 0 <= 1".into());
        }
        Ok(())
    }

    /// 615-LOGABOVELT: `ln 1 + 1 = 1 < ln 100`, so `Above` refuses `w = 1` at `x = 100`; at `x = 1` (`1 >= 0`) it accepts.
    #[test]
    fn log_test_above_refuses_less() -> Result<(), String> {
        let mut cc = Consts::new().map_err(|e| format!("Consts::new: {:?}", e))?;
        // no-source: hand-derived: ln 1 + 1 = 1 < ln 100 ≈ 4.6 (refused by Above), 1 >= ln 1 = 0 (accepted).
        if log_test(&num(1), &num(100), &Side::Above, 128, &mut cc) {
            return Err("log_test(1, 100, Above) accepted 1 >= ln 100".into());
        }
        if !log_test(&num(1), &num(1), &Side::Above, 128, &mut cc) {
            return Err("log_test(1, 1, Above) refused 1 >= 0".into());
        }
        Ok(())
    }

    /// 626-RADBELOW64: the ball from `result` contains its certified `[l, h]`, here `l = -1`, `h = -1 + k/7` at 400 bits,
    /// where `mid - l` needs more than 64 bits.
    #[test]
    fn result_contains_endpoints() -> Result<(), String> {
        let mut bad = Vec::new();
        // no-source: hand-chosen inputs: l = -1, h = k/7 - 1 (non-dyadic, so mid - l needs > 64 bits); the expected
        // outcome is the enclosure property (result contains l and h), not a computed number.
        for prec in [100, 200, 300] {
            for k in 1..17 {
                let l = num(-1);
                let h = num(k)
                    .div(&num(7), 400, RoundingMode::Down)
                    .sub_full_prec(&num(1));
                match result("lambert_w0", &l, &h, prec) {
                    Ok(b) if b.contains(&l) == Tri::True && b.contains(&h) == Tri::True => {}
                    got => bad.push(format!("prec {} k {}: {:?}", prec, k, got)),
                }
            }
        }
        if bad.is_empty() {
            Ok(())
        } else {
            Err(format!(
                "{} of 48 lose an endpoint:\n{}",
                bad.len(),
                bad.join("\n")
            ))
        }
    }

    /// 626-RADBELOW64 at low `prec`: with a 2-to-40-bit midpoint and `mid - l` needing about 90 bits, rounding that half-width
    /// down can drop `l` when the midpoint rounds towards `h`.
    #[test]
    fn result_contains_endpoints_low_prec() -> Result<(), String> {
        // no-source: hand-chosen inputs: l = 1 + 2^-90, h = j/7 (j = 8..63) at 400 bits; the expected outcome is the
        // enclosure property (result contains l and h), not a computed number.
        let tiny = num(1).ldexp(-90, 8, RoundingMode::Down);
        let mut bad = Vec::new();
        for prec in 2..41 {
            for j in 8..64 {
                let l = num(1).add_full_prec(&tiny);
                let h = num(j).div(&num(7), 400, RoundingMode::Down);
                match result("lambert_w0", &l, &h, prec) {
                    Ok(b) if b.contains(&l) == Tri::True && b.contains(&h) == Tri::True => {}
                    got => bad.push(format!("prec {} j {}: {:?}", prec, j, got)),
                }
            }
        }
        if bad.is_empty() {
            Ok(())
        } else {
            Err(format!(
                "{} of {} lose an endpoint:\n{}",
                bad.len(),
                39 * 56,
                bad.join("\n")
            ))
        }
    }
}
