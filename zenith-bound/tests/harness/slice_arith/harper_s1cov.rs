//! Harper S1 coverage rows for abs, mul_2exp and the Tri comparisons (CONTRACT 15.1, 15.2, 16.S1.10, 16.S1.11).
//! Hand-verified rows at the edges the parent rows leave out: i64::MIN and i64::MIN + 1 midpoints and shift counts,
//! signed zero on both sides, round-to-even ties of the midpoint, |m| = r exactly, touching endpoints, and every
//! precision error path.

use zenith_bound::ball::Ball;
use zenith_bound::enclosure::Enclosure;
use zenith_bound::error::BoundError;
use zenith_bound::mag::Mag;
use zenith_bound::tri::Tri;
use zenith_float_num::ExactNum;

use super::{
    ball, dy, mag, num, ok, pow2, pt, want_all_reals, want_contains, want_exact, want_invalid, TOP,
};

const P: usize = 64;
const T: Tri = Tri::True;
const F: Tri = Tri::False;
const U: Tri = Tri::Unknown;

fn min1() -> Result<i64, String> {
    i64::MIN
        .checked_add(1)
        .ok_or_else(|| "i64::MIN + 1".to_string())
}

/// `2^e + v` exactly.
fn p2v(e: i64, v: i64) -> Result<ExactNum, String> {
    Ok(pow2(e)?.add_full_prec(&num(v)))
}

fn exact_at(x: ExactNum, prec: usize) -> Result<Ball, String> {
    ok("exact", Ball::exact(x, prec))
}

fn max1() -> Result<usize, String> {
    zenith_bound::MAX_PRECISION
        .checked_add(1)
        .ok_or_else(|| "MAX_PRECISION + 1".to_string())
}

/// Want a rounded (inexact) ball whose midpoint is exactly `mid` and which still contains `x`.
fn want_rounded(what: &str, b: &Ball, mid: &ExactNum, x: &ExactNum) -> Result<(), String> {
    if b.mid() != mid || b.is_exact() {
        return Err(format!("{}: got {:?}, want inexact mid {:?}", what, b, mid));
    }
    want_contains(what, b, x)
}

/// Want `[lt, le, gt, ge, eq, ne, overlaps]` of `a` against `b`.
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
fn s1cov_abs_extremes_ties_and_zero() -> Result<(), String> {
    // |i64::MIN| = 2^63, one bit, exact at 64 bits.
    want_exact(
        "|MIN|",
        &ok("abs", pt(i64::MIN, P)?.abs(P))?,
        &pow2(63)?,
        &Mag::zero(),
    )?;
    // |i64::MIN + 1| = 2^63 - 1 needs 63 bits: exact at 64, rounded at 53 (2^63 - 1 -> 2^63, ToEven nearest).
    let m1 = min1()?;
    want_exact(
        "|MIN+1|@64",
        &ok("abs", pt(m1, P)?.abs(P))?,
        &num(i64::MAX),
        &Mag::zero(),
    )?;
    want_rounded(
        "|MIN+1|@53",
        &ok("abs", pt(m1, P)?.abs(53))?,
        &pow2(63)?,
        &num(i64::MAX),
    )?;
    // Ties at 53 bits round to even: |-(2^53 + 1)| -> 2^53, |-(2^53 + 3)| -> 2^53 + 4.
    let b = ok("abs", exact_at(-p2v(53, 1)?, 128)?.abs(53))?;
    want_rounded("|-(2^53+1)|@53", &b, &pow2(53)?, &p2v(53, 1)?)?;
    let b = ok("abs", exact_at(-p2v(53, 3)?, 128)?.abs(53))?;
    want_rounded("|-(2^53+3)|@53", &b, &p2v(53, 4)?, &p2v(53, 3)?)?;
    // |m| = r exactly: |-2 ± 2| = 2 ± 2 (the ball touches 0 from the left; the non-straddling branch keeps r).
    want_exact(
        "|-2±2|",
        &ok("abs", ball(num(-2), mag(TOP, 2)?, P)?.abs(P))?,
        &num(2),
        &mag(TOP, 2)?,
    )?;
    // +0 is +0; -0 ± 1 is the endpoint ball [0, 1] = 1/2 ± 1/2.
    let z = ok("abs", pt(0, P)?.abs(P))?;
    if !z.mid().is_zero() || z.mid().is_negative() || !z.is_exact() {
        return Err(format!("|+0|: got {:?}", z));
    }
    let b = ok("abs", ball(-num(0), mag(TOP, 1)?, P)?.abs(P))?;
    want_exact("|-0±1|", &b, &dy(1, -1)?, &mag(TOP, 0)?)?;
    // A ball around i64::MIN straddling nothing: |MIN ± 1| = 2^63 ± 1.
    want_exact(
        "|MIN±1|",
        &ok("abs", ball(num(i64::MIN), mag(TOP, 1)?, P)?.abs(P))?,
        &pow2(63)?,
        &mag(TOP, 1)?,
    )?;
    Ok(())
}

#[test]
fn s1cov_mul_2exp_extremes_ties_and_zero() -> Result<(), String> {
    let m1 = min1()?;
    // i64::MIN · 2^-63 = -1 and i64::MIN · 2 = -2^64, both exact.
    want_exact(
        "MIN·2^-63",
        &ok("m2e", pt(i64::MIN, P)?.mul_2exp(-63, P))?,
        &num(-1),
        &Mag::zero(),
    )?;
    want_exact(
        "MIN·2",
        &ok("m2e", pt(i64::MIN, P)?.mul_2exp(1, P))?,
        &-pow2(64)?,
        &Mag::zero(),
    )?;
    // Radius exponent 1 + (MIN + 1) = MIN + 2 still fits: 0 ± 1 scaled by 2^(MIN+1) is exact.
    let b = ok("m2e", ball(num(0), mag(TOP, 1)?, P)?.mul_2exp(m1, P))?;
    want_exact(
        "0±1·2^(MIN+1)",
        &b,
        &num(0),
        &mag(TOP, i64::MIN.checked_add(2).ok_or("MIN+2")?)?,
    )?;
    // Radius exponent -5 + (MIN + 1) underflows: the smallest positive Mag, never zero.
    let b = ok("m2e", ball(num(0), mag(TOP, -5)?, P)?.mul_2exp(m1, P))?;
    want_exact("0±2^-6·2^(MIN+1)", &b, &num(0), &mag(TOP, i64::MIN)?)?;
    // A nonzero exact midpoint shifted by MIN + 1 underflows to an enclosure of 0, not a silent exact 0.
    let b = ok("m2e", pt(1, P)?.mul_2exp(m1, P))?;
    if b.is_exact() || b.contains(&num(0)) != Tri::True {
        return Err(format!(
            "1·2^(MIN+1): want an inexact ball around 0, got {:?}",
            b
        ));
    }
    // Signed zero stays an exact zero under any shift.
    for k in [i64::MIN, m1, -1, 0, 1, i64::MAX] {
        let z = ok("m2e", (-pt(0, P)?).mul_2exp(k, P))?;
        if !z.mid().is_zero() || !z.is_exact() {
            return Err(format!("-0·2^{}: got {:?}", k, z));
        }
    }
    // Ties at 53 bits round to even: (2^53 + 1)·2 -> 2^54, (2^53 + 3)·2^-1 -> 2^52 + 2.
    let b = ok("m2e", exact_at(p2v(53, 1)?, 128)?.mul_2exp(1, 53))?;
    want_rounded("(2^53+1)·2@53", &b, &pow2(54)?, &p2v(54, 2)?)?;
    let b = ok("m2e", exact_at(p2v(53, 3)?, 128)?.mul_2exp(-1, 53))?;
    want_rounded(
        "(2^53+3)/2@53",
        &b,
        &p2v(52, 2)?,
        &p2v(53, 3)?.mul_full_prec(&pow2(-1)?),
    )?;
    // Midpoint overflow at the top shift is all reals.
    want_all_reals("MIN·2^MAX", pt(i64::MIN, P)?.mul_2exp(i64::MAX, P))?;
    Ok(())
}

#[test]
fn s1cov_abs_mul_2exp_prec_errors() -> Result<(), String> {
    let hi = max1()?;
    want_invalid("abs prec 0", pt(i64::MIN, P)?.abs(0))?;
    want_invalid("abs prec max+1", pt(i64::MIN, P)?.abs(hi))?;
    want_invalid("abs -0 prec max+1", (-pt(0, P)?).abs(hi))?;
    want_invalid("m2e prec 0 k MIN", pt(1, P)?.mul_2exp(i64::MIN, 0))?;
    want_invalid("m2e prec max+1 k MAX", pt(1, P)?.mul_2exp(i64::MAX, hi))?;
    // The error is checked before any overflow: all reals with a bad prec is still InvalidArgument.
    want_invalid(
        "m2e all prec 0",
        ok("ar", Ball::all_reals(P))?.mul_2exp(1, 0),
    )?;
    match pt(1, P)?.abs(0) {
        Err(BoundError::InvalidArgument { .. }) => Ok(()),
        other => Err(format!("abs prec 0: got {:?}", other)),
    }
}

#[test]
fn s1cov_compare_extremes_signed_zero_and_ties() -> Result<(), String> {
    let m1 = min1()?;
    // -0 and +0 are one value in every argument order.
    let nz = -pt(0, P)?;
    want_cmp("-0 vs +0", &nz, &pt(0, P)?, [F, T, F, T, T, F, T])?;
    want_cmp("+0 vs -0", &pt(0, P)?, &nz, [F, T, F, T, T, F, T])?;
    // i64::MIN against i64::MIN + 1, both orders, and MIN against itself.
    want_cmp(
        "MIN vs MIN+1",
        &pt(i64::MIN, P)?,
        &pt(m1, P)?,
        [T, T, F, F, F, T, F],
    )?;
    want_cmp(
        "MIN+1 vs MIN",
        &pt(m1, P)?,
        &pt(i64::MIN, P)?,
        [F, F, T, T, F, T, F],
    )?;
    want_cmp(
        "MIN vs MIN",
        &pt(i64::MIN, P)?,
        &pt(i64::MIN, P)?,
        [F, T, F, T, T, F, T],
    )?;
    want_cmp(
        "MIN vs MAX",
        &pt(i64::MIN, P)?,
        &pt(i64::MAX, P)?,
        [T, T, F, F, F, T, F],
    )?;
    // Ties (touching endpoints): MIN ± 1 touches MIN + 1; 0 ± 1 touches 2 ± 1 at 1; -0 ± 1 touches 1. The upper ball is
    // certainly not below the lower one, so its lt (the lower one's gt) is False (review nit a).
    let mb = ball(num(i64::MIN), mag(TOP, 1)?, P)?;
    want_cmp("MIN±1 vs MIN+1", &mb, &pt(m1, P)?, [U, T, F, U, U, U, T])?;
    want_cmp("MIN+1 vs MIN±1", &pt(m1, P)?, &mb, [F, U, U, T, U, U, T])?;
    let a = ball(num(0), mag(TOP, 1)?, P)?;
    want_cmp(
        "0±1 vs 2±1",
        &a,
        &ball(num(2), mag(TOP, 1)?, P)?,
        [U, T, F, U, U, U, T],
    )?;
    let nzb = ball(-num(0), mag(TOP, 1)?, P)?;
    want_cmp("-0±1 vs 1", &nzb, &pt(1, P)?, [U, T, F, U, U, U, T])?;
    want_cmp("-0±1 vs +0", &nzb, &pt(0, P)?, [U, U, U, U, U, U, T])?;
    // One step past the tie is decided: 0 ± 1 vs 2 ± 1/2 is disjoint by 1/2.
    want_cmp(
        "0±1 vs 2±1/2",
        &a,
        &ball(num(2), mag(TOP, 0)?, P)?,
        [T, T, F, F, F, T, F],
    )?;
    // All reals overlaps everything and orders nothing, including signed zero and i64::MIN.
    let all = ok("ar", Ball::all_reals(P))?;
    want_cmp("all vs -0", &all, &nz, [U, U, U, U, U, U, T])?;
    want_cmp("MIN vs all", &pt(i64::MIN, P)?, &all, [U, U, U, U, U, U, T])?;
    Ok(())
}

fn tri_not(t: Tri) -> Tri {
    match t {
        Tri::True => Tri::False,
        Tri::False => Tri::True,
        Tri::Unknown => Tri::Unknown,
    }
}

/// Duality over a pool with the extremes, signed zero, touching balls and all reals: `gt` is the swapped `lt`, `ge`
/// the swapped `le`, `ne` the negated `eq`, `overlaps` is symmetric, and a decided `lt` implies `le`.
#[test]
fn s1cov_compare_duality_pool() -> Result<(), String> {
    let m1 = min1()?;
    let pool = [
        pt(i64::MIN, P)?,
        pt(m1, P)?,
        pt(i64::MAX, P)?,
        pt(0, P)?,
        -pt(0, P)?,
        ball(-num(0), mag(TOP, 1)?, P)?,
        ball(num(i64::MIN), mag(TOP, 1)?, P)?,
        ball(num(2), mag(TOP, 1)?, P)?,
        pt(1, P)?,
        ok("ar", Ball::all_reals(P))?,
    ];
    for (i, a) in pool.iter().enumerate() {
        for (j, b) in pool.iter().enumerate() {
            let w = format!("pool {} vs {}", i, j);
            if a.gt(b) != b.lt(a) {
                return Err(format!(
                    "{}: gt {:?} != swapped lt {:?}",
                    w,
                    a.gt(b),
                    b.lt(a)
                ));
            }
            if a.ge(b) != b.le(a) {
                return Err(format!(
                    "{}: ge {:?} != swapped le {:?}",
                    w,
                    a.ge(b),
                    b.le(a)
                ));
            }
            if a.ne(b) != tri_not(a.eq(b)) {
                return Err(format!("{}: ne {:?} vs eq {:?}", w, a.ne(b), a.eq(b)));
            }
            if a.overlaps(b) != b.overlaps(a) {
                return Err(format!("{}: overlaps not symmetric", w));
            }
            if a.lt(b) == Tri::True && a.le(b) != Tri::True {
                return Err(format!("{}: lt True but le {:?}", w, a.le(b)));
            }
            if a.eq(b) == Tri::True && a.overlaps(b) != Tri::True {
                return Err(format!("{}: eq True but overlaps {:?}", w, a.overlaps(b)));
            }
        }
    }
    Ok(())
}
