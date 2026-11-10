//! Harper L0b: abs, mul_2exp and the Tri comparisons of slice S1 (CONTRACT 15.1, 15.2, 16.S1.10, 16.S1.11).
//! Hand-verified rows: zero and signed zero, straddling 0, touching endpoints, undecidable (Unknown) orderings,
//! all-reals and out-of-range radii, mul_2exp near EXPONENT_MIN / EXPONENT_MAX, and bad precision.
//! COVERS: Ball::abs, Ball::mul_2exp, Ball::lt, Ball::le, Ball::gt, Ball::ge, Ball::eq, Ball::ne, Ball::overlaps

use zenith_bound::ball::Ball;
use zenith_bound::enclosure::Enclosure;
use zenith_bound::error::BoundError;
use zenith_bound::mag::Mag;
use zenith_bound::tri::Tri;
use zenith_float_num::{ExactNum, EXPONENT_MAX, EXPONENT_MIN};

use super::{ball, mag, num, ok, pow2, pt, want_contains, want_exact, TOP};

const P: usize = 53;
const T: Tri = Tri::True;
const F: Tri = Tri::False;
const U: Tri = Tri::Unknown;

/// `v ± 2^k` (Mag `finite(TOP, k + 1)` is `2^k`).
fn bk(v: i64, k: i64) -> Result<Ball, String> {
    ball(num(v), mag(TOP, k.checked_add(1).ok_or("k")?)?, P)
}

fn all() -> Result<Ball, String> {
    ok("all_reals", Ball::all_reals(P))
}

/// One row: want `[lt, le, gt, ge, eq, ne, overlaps]` of `a` against `b`.
fn want_cmp(what: &str, a: &Ball, b: &Ball, want: [Tri; 7]) -> Result<(), String> {
    let got = [
        a.lt(b),
        a.le(b),
        a.gt(b),
        a.ge(b),
        a.eq(b),
        a.ne(b),
        a.overlaps(b),
    ];
    if got != want {
        return Err(format!(
            "{}: [lt le gt ge eq ne overlaps] got {:?}, want {:?} ({:?} vs {:?})",
            what, got, want, a, b
        ));
    }
    Ok(())
}

#[test]
fn compare_exact_points() -> Result<(), String> {
    // 1 < 2: every relation decided (kills 274 eq-False, 275 ne-True, 282 swapped D).
    want_cmp("1 vs 2", &pt(1, P)?, &pt(2, P)?, [T, T, F, F, F, T, F])?;
    want_cmp("2 vs 1", &pt(2, P)?, &pt(1, P)?, [F, F, T, T, F, T, F])?;
    // equal exact points: eq True (273), overlaps True (278).
    want_cmp("2 vs 2", &pt(2, P)?, &pt(2, P)?, [F, T, F, T, T, F, T])?;
    // -0 and +0 are the same value.
    let nz = -pt(0, P)?;
    want_cmp("-0 vs +0", &nz, &pt(0, P)?, [F, T, F, T, T, F, T])?;
    want_cmp("-0 vs -0", &nz, &-pt(0, P)?, [F, T, F, T, T, F, T])?;
    want_cmp("-1 vs 0", &pt(-1, P)?, &pt(0, P)?, [T, T, F, F, F, T, F])
}

#[test]
fn compare_balls_disjoint_touching_straddling() -> Result<(), String> {
    // [-1, 1] vs 5 and its mirror: disjoint, decided.
    want_cmp("0±1 vs 5", &bk(0, 0)?, &pt(5, P)?, [T, T, F, F, F, T, F])?;
    want_cmp("5 vs 0±1", &pt(5, P)?, &bk(0, 0)?, [F, F, T, T, F, T, F])?;
    // [0, 2] vs [2, 4] touch at 2 (CONTRACT 15.2 touching edges): le/ge on the right side True; the upper ball is
    // certainly not below the lower one (R > 0, lo(upper) = hi(lower)), so its lt (the lower one's gt) is False (S1 nit a,
    // spence/zb-s1 026e979); the other strict order and eq stay Unknown; they overlap (279/280 flip the D ± R signs).
    want_cmp("1±1 vs 3±1", &bk(1, 0)?, &bk(3, 0)?, [U, T, F, U, U, U, T])?;
    want_cmp("3±1 vs 1±1", &bk(3, 0)?, &bk(1, 0)?, [F, U, U, T, U, U, T])?;
    // [-2, 2] vs 1: straddles, nothing decided but overlap (281 uses r1 - r2).
    want_cmp("0±2 vs 1", &bk(0, 1)?, &pt(1, P)?, [U, U, U, U, U, U, T])?;
    // R = r1 + r2 must be the sum: (1 ± 1/2) vs (2 ± 1/4) is disjoint by 1/4; (1 ± 1/2) vs (2 ± 1/2) touches at 3/2.
    want_cmp(
        "1±1/2 vs 2±1/4",
        &bk(1, -1)?,
        &bk(2, -2)?,
        [T, T, F, F, F, T, F],
    )?;
    want_cmp(
        "1±1/2 vs 2±1/2",
        &bk(1, -1)?,
        &bk(2, -1)?,
        [U, T, F, U, U, U, T],
    )?;
    // overlapping but neither contains the other: (0 ± 1) vs (1 ± 1).
    want_cmp("0±1 vs 1±1", &bk(0, 0)?, &bk(1, 0)?, [U, U, U, U, U, U, T])?;
    // a ball against itself (rebuilt): not a point, so = is Unknown but it overlaps.
    want_cmp("3±1 vs 3±1", &bk(3, 0)?, &bk(3, 0)?, [U, U, U, U, U, U, T])
}

#[test]
fn compare_infinite_and_out_of_range_radius() -> Result<(), String> {
    // all reals: every order Unknown, overlaps True (276), in both argument orders.
    want_cmp("all vs 1", &all()?, &pt(1, P)?, [U, U, U, U, U, U, T])?;
    want_cmp("1 vs all", &pt(1, P)?, &all()?, [U, U, U, U, U, U, T])?;
    want_cmp("all vs all", &all()?, &all()?, [U, U, U, U, U, U, T])?;
    // a finite radius outside zf's range: orders undecidable; since S2, cmp_exact_num (16.S2.1) decides overlaps True, never
    // False (277).
    let huge = ball(
        num(1),
        mag(TOP, EXPONENT_MAX.checked_add(8).ok_or("emax+8")?)?,
        P,
    )?;
    want_cmp("huge vs 1", &huge, &pt(1, P)?, [U, U, U, U, U, U, T])?;
    want_cmp("1 vs huge", &pt(1, P)?, &huge, [U, U, U, U, U, U, T])
}

#[test]
fn abs_rows() -> Result<(), String> {
    // exact points: |-3| = 3, |3| = 3, rad 0.
    want_exact(
        "|-3|",
        &ok("abs", pt(-3, P)?.abs(P))?,
        &num(3),
        &Mag::zero(),
    )?;
    want_exact("|3|", &ok("abs", pt(3, P)?.abs(P))?, &num(3), &Mag::zero())?;
    // clear of 0: |-3 ± 1| = 3 ± 1 keeps the radius (283 drops it).
    want_exact(
        "|-3±1|",
        &ok("abs", bk(-3, 0)?.abs(P))?,
        &num(3),
        &mag(TOP, 1)?,
    )?;
    // touching 0: |-1 ± 1| covers [0, 2].
    let t = ok("abs", bk(-1, 0)?.abs(P))?;
    for x in [num(0), num(1), num(2)] {
        want_contains("|-1±1|", &t, &x)?;
    }
    if t.mid().is_negative() {
        return Err(format!("|-1±1|: negative mid {:?}", t));
    }
    // straddling: |-1 ± 2| = [0, 3]; contains 0 and 3, nothing below 0 beyond the endpoint enclosure.
    let s = ok("abs", bk(-1, 1)?.abs(P))?;
    for x in [num(0), num(3), num(1)] {
        want_contains("|-1±2|", &s, &x)?;
    }
    if s.contains(&num(-1)) == Tri::True {
        return Err(format!("|-1±2| contains -1: {:?}", s));
    }
    // signed zero: |-0| = +0 exactly.
    let z = ok("abs", (-pt(0, P)?).abs(P))?;
    if !z.mid().is_zero() || z.mid().is_negative() || !z.rad().is_zero() {
        return Err(format!("|-0|: {:?}, want +0 ± 0", z));
    }
    // rounding to prec: |-(2^60 + 1)| at 53 bits is inexact, so the radius is positive and contains the value.
    let big = num(1 << 60).add_full_prec(&num(1));
    let r = ok(
        "abs",
        ok("exact", Ball::exact(-big.add_full_prec(&num(0)), 128))?.abs(P),
    )?;
    want_contains("|-(2^60+1)|@53", &r, &big)?;
    if r.rad().is_zero() {
        return Err(format!("|-(2^60+1)|@53: rad 0 {:?}", r));
    }
    // all reals stays all reals.
    if ok("abs", all()?.abs(P))?.is_informative() {
        return Err("|all reals| must be all reals".into());
    }
    Ok(())
}

#[test]
fn mul_2exp_rows() -> Result<(), String> {
    // (3 ± 1) · 2^4 = 48 ± 16 and (3 ± 1) · 2^-2 = 3/4 ± 1/4, exact.
    want_exact(
        "(3±1)·2^4",
        &ok("m2e", bk(3, 0)?.mul_2exp(4, P))?,
        &num(48),
        &mag(TOP, 5)?,
    )?;
    want_exact(
        "(3±1)·2^-2",
        &ok("m2e", bk(3, 0)?.mul_2exp(-2, P))?,
        &num(3).mul_full_prec(&pow2(-2)?),
        &mag(TOP, -1)?,
    )?;
    want_exact(
        "1·2^0",
        &ok("m2e", pt(1, P)?.mul_2exp(0, P))?,
        &num(1),
        &Mag::zero(),
    )?;
    // -0 · 2^5 stays a zero with rad 0.
    let z = ok("m2e", (-pt(0, P)?).mul_2exp(5, P))?;
    if !z.mid().is_zero() || !z.rad().is_zero() {
        return Err(format!("-0·2^5: {:?}", z));
    }
    // more bits than prec: (2^60 + 1) · 2 at 53 bits rounds, so it gets a radius that contains the exact value.
    let big = num(1 << 60).add_full_prec(&num(1));
    let r = ok(
        "m2e",
        ok("exact", Ball::exact(big.add_full_prec(&num(0)), 128))?.mul_2exp(1, P),
    )?;
    want_contains("(2^60+1)·2@53", &r, &big.mul_full_prec(&num(2)))?;
    if r.rad().is_zero() {
        return Err(format!("(2^60+1)·2@53: rad 0 {:?}", r));
    }
    // midpoint overflow: 1 · 2^EXPONENT_MAX, and 2^(EMAX-1) · 2, are all reals.
    for (what, b) in [
        ("1·2^EMAX", pt(1, P)?.mul_2exp(EXPONENT_MAX, P)),
        (
            "2^(EMAX-1)·2",
            ok("big", Ball::exact(pow2(EXPONENT_MAX - 1)?, P))?.mul_2exp(1, P),
        ),
        ("1·2^i64::MAX", pt(1, P)?.mul_2exp(i64::MAX, P)),
    ] {
        if ok(what, b)?.is_informative() {
            return Err(format!("{}: want all reals", what));
        }
    }
    // just inside: 2^(EMAX-2) · 2 = 2^(EMAX-1) exactly.
    want_exact(
        "2^(EMAX-2)·2",
        &ok(
            "m2e",
            ok("big", Ball::exact(pow2(EXPONENT_MAX - 2)?, P))?.mul_2exp(1, P),
        )?,
        &pow2(EXPONENT_MAX - 1)?,
        &Mag::zero(),
    )?;
    // underflow: 1 ± 2^-10 scaled far below EXPONENT_MIN is enclosed, never a silent flush: it contains 0 and has a
    // positive radius; so does an exact point.
    for k in [EXPONENT_MIN - 200, EXPONENT_MIN.saturating_mul(2), i64::MIN] {
        for (what, b) in [("1±2^-10", bk(1, -10)?), ("1", pt(1, P)?)] {
            let r = ok("m2e", b.mul_2exp(k, P))?;
            if r.rad().is_zero() {
                return Err(format!("{}·2^{}: rad 0 {:?}", what, k, r));
            }
            want_contains(&format!("{}·2^{}", what, k), &r, &num(0))?;
        }
    }
    // a radius pushed past Mag's exponent range is Infinity (all reals), not a wrapped small Mag.
    let r = ok(
        "m2e",
        ball(num(0), mag(TOP, 1 << 40)?, P)?.mul_2exp(i64::MAX, P),
    )?;
    if r.is_informative() {
        return Err(format!("rad·2^i64::MAX: {:?}", r));
    }
    // all reals stays all reals.
    if ok("m2e", all()?.mul_2exp(-3, P))?.is_informative() {
        return Err("all reals · 2^-3 must be all reals".into());
    }
    Ok(())
}

#[test]
fn abs_mul_2exp_bad_prec() -> Result<(), String> {
    let max1 = zenith_bound::MAX_PRECISION.checked_add(1).ok_or("max+1")?;
    let checks: [Result<Ball, BoundError>; 4] = [
        pt(1, P)?.abs(0),
        pt(1, P)?.abs(max1),
        pt(1, P)?.mul_2exp(1, 0),
        pt(1, P)?.mul_2exp(1, max1),
    ];
    for (i, r) in checks.into_iter().enumerate() {
        match r {
            Err(BoundError::InvalidArgument { .. }) => {}
            other => return Err(format!("case {}: want InvalidArgument, got {:?}", i, other)),
        }
    }
    Ok(())
}

/// Seeded cross-check: the Tri answers agree with exact endpoint arithmetic and with each other.
#[test]
fn compare_consistency_seeded() -> Result<(), String> {
    let mut rng = crate::common::rng::SplitMix64::new(0x2180_0000_0000_0006);
    for case in 0..2000u32 {
        let v = |rng: &mut crate::common::rng::SplitMix64| -> Result<i64, String> {
            i64::try_from(rng.below(9))
                .map_err(|_| "v".to_string())?
                .checked_sub(4)
                .ok_or_else(|| "v".to_string())
        };
        let (ma, mb) = (v(&mut rng)?, v(&mut rng)?);
        let ka = i64::try_from(rng.below(4))
            .map_err(|_| "k")?
            .checked_sub(2)
            .ok_or("k")?;
        let kb = i64::try_from(rng.below(4))
            .map_err(|_| "k")?
            .checked_sub(2)
            .ok_or("k")?;
        let a = if rng.below(4) == 0 {
            pt(ma, P)?
        } else {
            bk(ma, ka)?
        };
        let b = if rng.below(4) == 0 {
            pt(mb, P)?
        } else {
            bk(mb, kb)?
        };
        let ends = |x: &Ball| -> Result<(ExactNum, ExactNum), String> {
            let r = x.rad().to_exact_num().map_err(|e| format!("{:?}", e))?;
            Ok((x.mid().sub_full_prec(&r), x.mid().add_full_prec(&r)))
        };
        let ((alo, ahi), (blo, bhi)) = (ends(&a)?, ends(&b)?);
        let point = a.rad().is_zero() && b.rad().is_zero();
        let tri = |t: bool, f: bool| {
            if t {
                T
            } else if f {
                F
            } else {
                U
            }
        };
        // 16.S1.11 from exact endpoints: touching edges decide le/ge, and (R > 0) the upper ball's lt / lower ball's gt
        // is False (S1 nit a).
        let peq = point && alo == blo;
        let a_on_top = !point && alo == bhi;
        let b_on_top = !point && blo == ahi;
        let disjoint = ahi < blo || bhi < alo;
        let want = [
            tri(ahi < blo, alo > bhi || peq || a_on_top),
            tri(ahi <= blo, alo > bhi),
            tri(bhi < alo, blo > ahi || peq || b_on_top),
            tri(bhi <= alo, blo > ahi),
            tri(peq, disjoint),
            tri(disjoint, peq),
            tri(!disjoint, disjoint),
        ];
        want_cmp(&format!("case {}", case), &a, &b, want)?;
    }
    Ok(())
}

/// 16.S1.11 with 16.S2.1/16.S2.2: a finite radius below zenith-float's range is still compared exactly when the other
/// radius is zero or inside the range (the in-range radius is converted, the tiny one only compared), so `0 ± 2^(EXPONENT_MIN
/// − 101)` is disjoint from 1 and from `1 ± 1/2`.
#[test]
fn compare_radius_below_zf_range_decided() -> Result<(), String> {
    let e = EXPONENT_MIN.checked_sub(100).ok_or("emin-100")?;
    let tiny = ball(num(0), mag(TOP, e)?, P)?;
    let one = pt(1, P)?;
    let one_half = ball(num(1), mag(TOP, 0)?, P)?;
    want_cmp("tiny vs 1", &tiny, &one, [T, T, F, F, F, T, F])?;
    want_cmp("1 vs tiny", &one, &tiny, [F, F, T, T, F, T, F])?;
    want_cmp("tiny vs 1±1/2", &tiny, &one_half, [T, T, F, F, F, T, F])?;
    want_cmp("1±1/2 vs tiny", &one_half, &tiny, [F, F, T, T, F, T, F])
}

/// Kills 277-OVFALSE: both radii `2^(EXPONENT_MIN − 65)` lie below zenith-float's range, so `0 ± r` against `2^(EXPONENT_MIN
/// − 64) ± r` (touching: `|D| = R`) is undecided; overlaps must be Unknown, never False.
#[test]
fn compare_both_radii_below_zf_range_touching_not_false() -> Result<(), String> {
    let e = EXPONENT_MIN.checked_sub(64).ok_or("emin-64")?;
    let a = ball(num(0), mag(TOP, e)?, P)?;
    let b = ball(pow2(e)?, mag(TOP, e)?, P)?;
    for (what, t) in [("a b", a.overlaps(&b)), ("b a", b.overlaps(&a))] {
        if t == Tri::False {
            return Err(format!(
                "{}: touching balls below the zf range reported disjoint",
                what
            ));
        }
    }
    Ok(())
}

/// Kills 343-FARSWAPSIGN: a = 0 ± 2^(EXPONENT_MIN − 150) (radius below zf's range), b = (2^(EXPONENT_MIN − 1) + 2^(EXPONENT_MIN −
/// 200)) ± 2^(EXPONENT_MIN − 1) at 256 bits. lo(b) = 2^(EXPONENT_MIN − 200) < hi(a), so the balls overlap: True both ways.
#[test]
fn compare_tiny_radius_against_fine_subnormal_gap() -> Result<(), String> {
    let e1 = EXPONENT_MIN.checked_sub(1).ok_or("e1")?;
    // 2^(EXPONENT_MIN − 200) is below pow2's range; an exact product of two in-range powers gives it as a wide subnormal.
    let h = EXPONENT_MIN.checked_div(2).ok_or("h")?;
    let fine = pow2(h.checked_sub(100).ok_or("h1")?)?.mul_full_prec(&pow2(
        EXPONENT_MIN
            .checked_sub(h)
            .and_then(|v| v.checked_sub(100))
            .ok_or("h2")?,
    )?);
    let m2 = pow2(e1)?.add_full_prec(&fine);
    let a = ball(
        num(0),
        mag(TOP, EXPONENT_MIN.checked_sub(149).ok_or("e149")?)?,
        256,
    )?;
    let b = ball(m2, mag(TOP, EXPONENT_MIN)?, 256)?;
    for (what, t) in [("a b", a.overlaps(&b)), ("b a", b.overlaps(&a))] {
        if t != Tri::True {
            return Err(format!(
                "{}: overlapping balls with a fine subnormal gap gave {:?}",
                what, t
            ));
        }
    }
    Ok(())
}
