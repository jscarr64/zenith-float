//! Paul's ruling 2026-10-23: `contains` is tight. Oracle: the exact integer order (tests/model/nat.rs) on
//! `|mid - x|` vs `rad`, all scaled to one common exponent. True and False must never be wrong, and with `rad`
//! inside zenith-float's exponent range the answer is never Unknown. Seeded splitmix64, replay with ZB_PROP_SEED/ZB_PROP_CASES.

use crate::model::nat::Nat;
use crate::model::rng::{run_cases, SplitMix64};
use core::cmp::Ordering;
use zenith_bound::ball::Ball;
use zenith_bound::enclosure::Enclosure;
use zenith_bound::mag::Mag;
use zenith_bound::tri::Tri;
use zenith_float_num::{ExactNum, RoundingMode};

const CASES: u64 = 3000;

/// A small signed mantissa (sometimes zero) and an exponent in [-400, 400].
fn draw_dyadic(rng: &mut SplitMix64) -> (i64, i64) {
    let a = match rng.next_u64() % 8 {
        0 => 0,
        _ => rng.next_i64() % (1i64 << 20),
    };
    let e = rng.next_i64() % 401;
    (a, e)
}

fn exact(a: i64, e: i64) -> Result<ExactNum, String> {
    Ok(ExactNum::from(a).ldexp(e, 64, RoundingMode::None))
}

/// `|v| << (e - base)`.
fn scaled(v: u64, e: i64, base: i64) -> Result<Nat, String> {
    let shift = usize::try_from(e - base).map_err(|_| "shift")?;
    Ok(Nat::from_u64(v).shl(shift))
}

/// Exact `|sum(mid) - sum(x)| <= m 2^er`, every value a sum of signed dyadic terms `a 2^e`.
fn reference(mid: &[(i64, i64)], x: &[(i64, i64)], m: u64, er: i64) -> Result<bool, String> {
    let base = mid.iter().chain(x.iter()).map(|t| t.1).fold(er, i64::min);
    let mut pos = Nat::zero();
    let mut neg = Nat::zero();
    for (a, e) in mid {
        let v = scaled(a.unsigned_abs(), *e, base)?;
        if *a < 0 {
            neg = neg.add(&v)
        } else {
            pos = pos.add(&v)
        }
    }
    for (a, e) in x {
        let v = scaled(a.unsigned_abs(), *e, base)?;
        if *a < 0 {
            pos = pos.add(&v)
        } else {
            neg = neg.add(&v)
        }
    }
    let big_r = scaled(m, er, base)?;
    let (hi, lo) = match pos.cmp(&neg) {
        Ordering::Less => (neg, pos),
        Ordering::Equal | Ordering::Greater => (pos, neg),
    };
    Ok(hi.cmp(&lo.add(&big_r)) != Ordering::Greater)
}

fn sum(terms: &[(i64, i64)]) -> Result<ExactNum, String> {
    let mut acc = ExactNum::new(4096);
    for (a, e) in terms {
        acc = acc.add(&exact(*a, *e)?, 4096, RoundingMode::None);
    }
    Ok(acc)
}

fn check(mid: (i64, i64), rad: Mag, x: &[(i64, i64)]) -> Result<(), String> {
    let m = rad.mant().ok_or("mant")?;
    let er = rad.exp().ok_or("exp")? - 64;
    let ball = Ball::from_mid_rad(exact(mid.0, mid.1)?, rad, 64).map_err(|e| format!("{:?}", e))?;
    let got = ball.contains(&sum(x)?);
    let want = reference(&[mid], x, m, er)?;
    let ok = match got {
        Tri::True => want,
        Tri::False => !want,
        Tri::Unknown => false,
    };
    if ok {
        Ok(())
    } else {
        Err(format!(
            "{:?} ± {}*2^{} vs {:?}: got {:?}, exact {}",
            mid, m, er, x, got, want
        ))
    }
}

#[test]
fn prop_contains_is_tight_against_the_exact_order() -> Result<(), String> {
    run_cases(
        "prop_contains_is_tight_against_the_exact_order",
        CASES,
        &|rng| {
            let (a, ea) = draw_dyadic(rng);
            let (b, eb) = match rng.next_u64() % 4 {
                // x a hair away from mid: the 2^-300 shape of the ruling's red case.
                0 => (rng.next_i64() % 3, ea - 300 + rng.next_i64() % 8),
                _ => draw_dyadic(rng),
            };
            let top = ea.max(eb);
            let rad = Mag::finite(rng.next_u64() | 1, top + 44 + rng.next_i64() % 48)
                .map_err(|e| format!("{:?}", e))?;
            check((a, ea), rad, &[(b, eb)])
        },
    )
}

#[test]
fn prop_contains_is_tight_at_the_edges() -> Result<(), String> {
    run_cases("prop_contains_is_tight_at_the_edges", CASES, &|rng| {
        // x = mid ± rad + t 2^et with 2^et far below rad: just inside, on, or just outside either edge.
        let (a, ea) = draw_dyadic(rng);
        let r = i64::try_from(rng.next_u64() % (1u64 << 20)).map_err(|_| "r")? + 1;
        let er = ea + rng.next_i64() % 24;
        let s = if rng.next_u64() % 2 == 0 { r } else { -r };
        let t = rng.next_i64() % 2;
        let et = er - 300 + rng.next_i64() % 8;
        let rad = Mag::finite(r.unsigned_abs(), er + 64).map_err(|e| format!("{:?}", e))?;
        check((a, ea), rad, &[(a, ea), (s, er), (t, et)])
    })
}
