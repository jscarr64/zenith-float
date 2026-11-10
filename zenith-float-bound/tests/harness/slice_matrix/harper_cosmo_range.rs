//! Harper cosmology-range rows for slice S6 (`src/matrix.rs`), 2026-10-29. Hand-worked oracles: every input is a power
//! of two or a small integer, so each exact answer is written by exponent bookkeeping (value beside each row).
//! Scales: 2^203 ~ 1.3e61 (Planck-to-Hubble), 2^399 ~ 1e120, 2^997 ~ 1e300, 2^20, 2^40, 2^61 and the zf limits
//! EXPONENT_MAX = 2^62 - 1 and EXPONENT_MIN = -EXPONENT_MAX. Precisions 1 to 1_000_000 bits.
//! CONTRACT 16.S6.3 (no error from spread or overflow), 16.S6.4 (rank-N add/sub), 16.S6.6 (norms), 16.S6.7 to 16.S6.9.

use zenith_float_bound::array::{BallArray, BallMat, BallVec};
use zenith_float_bound::ball::Ball;
use zenith_float_bound::enclosure::Enclosure;
use zenith_float_bound::error::BoundError;
use zenith_float_bound::mag::Mag;
use zenith_float_bound::tri::Tri;
use zenith_float_bound::MAX_PRECISION;
use zenith_float_num::{ExactNum, EXPONENT_MAX, EXPONENT_MIN};

const TOP: u64 = 0x8000_0000_0000_0000;
/// Exponent spreads: Planck-to-Hubble, 1e120, 1e300, 2^20, 2^40, 2^61.
const SPREADS: [i64; 6] = [203, 399, 997, 1 << 20, 1 << 40, 1 << 61];
/// Precisions from 1 bit to a million bits.
const PRECS: [usize; 8] = [1, 2, 53, 64, 1000, 10_000, 100_000, 1_000_000];

fn ok<T>(what: &str, r: Result<T, BoundError>) -> Result<T, String> {
    r.map_err(|e| format!("{what}: unexpected {e:?}"))
}

/// `2^k` exactly.
fn pow2(k: i64) -> Result<ExactNum, String> {
    let e = k.checked_add(1).ok_or("pow2 exponent")?;
    ExactNum::from_mag_parts(TOP, e, 64).map_err(|e| format!("pow2({k}): {e:?}"))
}

/// `-2^k` exactly.
fn npow2(k: i64) -> Result<ExactNum, String> {
    Ok(-pow2(k)?)
}

/// The Mag `mant * 2^(exp - 64)`.
fn mag(mant: u64, exp: i64) -> Result<Mag, String> {
    ok("Mag::finite", Mag::finite(mant, exp))
}

fn pt(x: ExactNum, prec: usize) -> Result<Ball, String> {
    ok("Ball::exact", Ball::exact(x, prec))
}

fn vec_of(xs: Vec<ExactNum>, prec: usize) -> Result<BallVec, String> {
    let mut balls = Vec::new();
    for x in xs {
        balls.push(pt(x, prec)?);
    }
    Ok(BallVec::from_balls(balls))
}

fn mat_of(nrows: usize, ncols: usize, xs: Vec<ExactNum>, prec: usize) -> Result<BallMat, String> {
    let mut balls = Vec::new();
    for x in xs {
        balls.push(pt(x, prec)?);
    }
    ok(
        "BallMat::from_balls",
        BallMat::from_balls(nrows, ncols, balls),
    )
}

fn want_contains(what: &str, b: &Ball, x: &ExactNum) -> Result<(), String> {
    match b.contains(x) {
        Tri::True => Ok(()),
        t @ (Tri::False | Tri::Unknown) => {
            Err(format!("{what}: {b:?} does not contain {x:?} ({t:?})"))
        }
    }
}

fn want_not_false(what: &str, t: Tri) -> Result<(), String> {
    match t {
        Tri::True | Tri::Unknown => Ok(()),
        Tri::False => Err(format!(
            "{what}: decided False, but the exact value says otherwise"
        )),
    }
}

fn want_le(what: &str, x: &Mag, bound: &Mag) -> Result<(), String> {
    if x <= bound {
        Ok(())
    } else {
        Err(format!("{what}: {x:?} exceeds {bound:?}"))
    }
}

fn want_ge(what: &str, x: &Mag, bound: &Mag) -> Result<(), String> {
    if x >= bound {
        Ok(())
    } else {
        Err(format!(
            "{what}: {x:?} is below {bound:?} (not an upper bound)"
        ))
    }
}

/// `x = (2^E, 1, -2^E)`, `y = (1, 1, 1)`: the exact dot is 1 (huge cancellation around a unit term), for every spread and
/// precision. Tightness (16.S6.9): TREF = 2^(1-prec)(2^(E+1) + 1), and 2 TREF <= 2^(E+3-prec)(1 + 2^-63) since E >= 203.
#[test]
fn cosmo_dot_cancel_spread_sweep() -> Result<(), String> {
    for e in SPREADS {
        for prec in PRECS {
            let what = format!("dot cancel E={e} prec={prec}");
            let x = vec_of(vec![pow2(e)?, ExactNum::from(1), npow2(e)?], prec)?;
            let y = vec_of(
                vec![ExactNum::from(1), ExactNum::from(1), ExactNum::from(1)],
                prec,
            )?;
            let d = ok(&what, x.dot(&y, prec))?;
            want_contains(&what, &d, &ExactNum::from(1))?;
            let p = i64::try_from(prec).map_err(|_| "prec as i64")?;
            let k = e
                .checked_add(4)
                .and_then(|v| v.checked_sub(p))
                .ok_or("bound exponent")?;
            want_le(&format!("{what} rad <= 2 TREF"), d.rad(), &mag(TOP | 1, k)?)?;
        }
    }
    Ok(())
}

/// `x = (2^E, 2^-E)`, `y = (2^-E, 2^E)`: both products are exactly 1, the dot is exactly 2 with radius 0 (16.S6.9 last
/// sentence: an all-exact entry has radius R1 = 0), at prec 2 to 10^6 and every spread up to 2^61.
#[test]
fn cosmo_dot_mixed_scale_exact() -> Result<(), String> {
    for e in SPREADS {
        for prec in [2, 53, 10_000, 1_000_000] {
            let what = format!("dot mixed E={e} prec={prec}");
            let ne = e.checked_neg().ok_or("neg")?;
            let x = vec_of(vec![pow2(e)?, pow2(ne)?], prec)?;
            let y = vec_of(vec![pow2(ne)?, pow2(e)?], prec)?;
            let d = ok(&what, x.dot(&y, prec))?;
            want_contains(&what, &d, &ExactNum::from(2))?;
            if !d.rad().is_zero() {
                return Err(format!(
                    "{what}: all-exact entry has radius {:?}, want 0",
                    d.rad()
                ));
            }
        }
    }
    Ok(())
}

/// Products that cross overflow: `2^(EMAX-2) * 2^4 = 2^(EMAX+2)` is past zf's range. 16.S6.3/16.S6.7: `Ok`, and that entry
/// alone is all reals; a mixed row `(2^(EMAX-2), 1) . (2^4, 1)` is all reals too (it holds the overflowing term).
#[test]
fn cosmo_dot_overflow_crossing() -> Result<(), String> {
    let top = EXPONENT_MAX.checked_sub(2).ok_or("emax-2")?;
    for prec in [1, 64, 100_000] {
        let what = format!("dot overflow prec={prec}");
        let x = vec_of(vec![pow2(top)?, ExactNum::from(1)], prec)?;
        let y = vec_of(vec![pow2(4)?, ExactNum::from(1)], prec)?;
        let d = ok(&what, x.dot(&y, prec))?;
        if !d.rad().is_infinity() {
            return Err(format!(
                "{what}: overflowing term gave a finite ball {d:?}, want all reals"
            ));
        }
        // matvec: row 0 overflows, row 1 = 2^(EMAX-2) * 2^-4 = 2^(EMAX-6) stays finite and exact.
        let a = mat_of(2, 1, vec![pow2(top)?, pow2(top)?], prec)?;
        let v = vec_of(vec![pow2(4)?], prec)?;
        let a2 = mat_of(2, 1, vec![pow2(-4)?, pow2(4)?], prec)?;
        let r = ok(&what, a.matvec(&v, prec))?;
        let r0 = r.get(0).ok_or("row 0")?;
        let r1 = r.get(1).ok_or("row 1")?;
        if !(r0.rad().is_infinity() && r1.rad().is_infinity()) {
            return Err(format!(
                "{what}: matvec rows {r0:?} {r1:?}, want both all reals"
            ));
        }
        let w = vec_of(vec![pow2(top)?], prec)?;
        let s = ok(&what, a2.matvec(&w, prec))?;
        let s0 = s.get(0).ok_or("row 0")?;
        let s1 = s.get(1).ok_or("row 1")?;
        let e6 = EXPONENT_MAX.checked_sub(6).ok_or("emax-6")?;
        want_contains(&format!("{what} finite row"), s0, &pow2(e6)?)?;
        if !s1.rad().is_infinity() {
            return Err(format!(
                "{what}: overflowing row 1 is {s1:?}, want all reals (row 0 stays finite)"
            ));
        }
    }
    Ok(())
}

/// Products that cross underflow: `(+-2^a)^2` with `a = EMIN/2 - 8` is `2^(2a)`, about `2^(EMIN-17)`, below zf's range.
/// The exact value is positive and below 2^EMIN, so `dot > 0` and `dot < 2^EMIN` must never be decided False.
#[test]
fn cosmo_dot_underflow_crossing() -> Result<(), String> {
    let a = EXPONENT_MIN
        .checked_div(2)
        .and_then(|v| v.checked_sub(8))
        .ok_or("a")?;
    for prec in [1, 64, 100_000] {
        let what = format!("dot underflow prec={prec}");
        let zero = pt(ExactNum::from(0), prec)?;
        let tiny = pt(pow2(EXPONENT_MIN)?, prec)?;
        let x = vec_of(vec![pow2(a)?], prec)?;
        let d = ok(&what, x.dot(&x, prec))?;
        want_not_false(&format!("{what} > 0"), d.gt(&zero))?;
        want_not_false(&format!("{what} < 2^EMIN"), d.lt(&tiny))?;
        let xn = vec_of(vec![npow2(a)?], prec)?;
        let dn = ok(&what, xn.dot(&x, prec))?;
        want_not_false(&format!("{what} -x.x < 0"), dn.lt(&zero))?;
    }
    Ok(())
}

/// Signed zero at huge scale: `(-0, 2^-E) . (2^E, -2^E) = -1`, and `(-0) . (2^(EMAX-1)) = 0` with a +0 midpoint (CONTRACT 6.3).
#[test]
fn cosmo_dot_signed_zero() -> Result<(), String> {
    let nz = -ExactNum::from(0);
    let big = EXPONENT_MAX.checked_sub(1).ok_or("emax-1")?;
    for prec in [1, 53, 100_000] {
        for e in SPREADS {
            let what = format!("dot -0 E={e} prec={prec}");
            let ne = e.checked_neg().ok_or("neg")?;
            let nz1 = -ExactNum::from(0);
            let x = vec_of(vec![nz1, pow2(ne)?], prec)?;
            let y = vec_of(vec![pow2(e)?, npow2(e)?], prec)?;
            let d = ok(&what, x.dot(&y, prec))?;
            want_contains(&what, &d, &ExactNum::from(-1))?;
        }
        let what = format!("dot -0 * 2^(EMAX-1) prec={prec}");
        let nzp = ExactNum::from(0).add_full_prec(&nz);
        let x = vec_of(vec![nzp], prec)?;
        let y = vec_of(vec![pow2(big)?], prec)?;
        let d = ok(&what, x.dot(&y, prec))?;
        want_contains(&what, &d, &ExactNum::from(0))?;
        if !d.mid().is_zero() || d.mid().is_negative() {
            return Err(format!("{what}: midpoint {:?}, want +0", d.mid()));
        }
    }
    Ok(())
}

/// `diag(2^E, 2^-E, 1) * diag(2^-E, 2^E, 1) = I` exactly: 1 on the diagonal, 0 elsewhere, radius 0 everywhere.
#[test]
fn cosmo_matmul_diag_spread_identity() -> Result<(), String> {
    for e in SPREADS {
        for prec in [1, 53, 10_000] {
            let what = format!("matmul diag E={e} prec={prec}");
            let ne = e.checked_neg().ok_or("neg")?;
            let z = || ExactNum::from(0);
            let a = mat_of(
                3,
                3,
                vec![
                    pow2(e)?,
                    z(),
                    z(),
                    z(),
                    pow2(ne)?,
                    z(),
                    z(),
                    z(),
                    ExactNum::from(1),
                ],
                prec,
            )?;
            let b = mat_of(
                3,
                3,
                vec![
                    pow2(ne)?,
                    z(),
                    z(),
                    z(),
                    pow2(e)?,
                    z(),
                    z(),
                    z(),
                    ExactNum::from(1),
                ],
                prec,
            )?;
            let c = ok(&what, a.matmul(&b, prec))?;
            for i in 0..3 {
                for j in 0..3 {
                    let cij = c.get(i, j).ok_or("entry")?;
                    let want = ExactNum::from(i64::from(i == j));
                    want_contains(&format!("{what} [{i}][{j}]"), cij, &want)?;
                    if !cij.rad().is_zero() {
                        return Err(format!("{what} [{i}][{j}]: radius {:?}, want 0", cij.rad()));
                    }
                }
            }
        }
    }
    Ok(())
}

/// Planck-to-Hubble matvec: rows `(2^203, -2^203, 2^-203)` and `(2^-399, 2^399, -2^399)` times `(1, 1, 1)`:
/// exact rows `2^-203` and `2^-399` (cancellation of the huge pair leaves the tiny term).
#[test]
fn cosmo_matvec_planck_hubble() -> Result<(), String> {
    for prec in PRECS {
        let what = format!("matvec planck/hubble prec={prec}");
        let a = mat_of(
            2,
            3,
            vec![
                pow2(203)?,
                npow2(203)?,
                pow2(-203)?,
                pow2(-399)?,
                pow2(399)?,
                npow2(399)?,
            ],
            prec,
        )?;
        let x = vec_of(
            vec![ExactNum::from(1), ExactNum::from(1), ExactNum::from(1)],
            prec,
        )?;
        let r = ok(&what, a.matvec(&x, prec))?;
        want_contains(
            &format!("{what} row 0"),
            r.get(0).ok_or("row 0")?,
            &pow2(-203)?,
        )?;
        want_contains(
            &format!("{what} row 1"),
            r.get(1).ok_or("row 1")?,
            &pow2(-399)?,
        )?;
    }
    Ok(())
}

/// Rank-3 `[2, 2, 2]` BallArray with entries `+-2^e`, `e` from EMIN+2000 to EMAX-1 (the whole zf range):
/// `a + (-a)` contains 0 entrywise, and `a + 2^(e-1000)` contains the exact sum `2^e + 2^(e-1000)` (16.S6.4 = 16.S1.1).
#[test]
fn cosmo_rank3_add_sub_full_range() -> Result<(), String> {
    let lo = EXPONENT_MIN.checked_add(2000).ok_or("emin+2000")?;
    let hi = EXPONENT_MAX.checked_sub(1).ok_or("emax-1")?;
    let es: [i64; 8] = [lo, -(1 << 61), -997, -203, 0, 399, 1 << 40, hi];
    for prec in [1, 64, 1001, 100_000] {
        let what = format!("rank3 prec={prec}");
        let shape = || -> Box<[usize]> { Box::new([2, 2, 2]) };
        let build = |f: &dyn Fn(i64) -> Result<ExactNum, String>| -> Result<BallArray, String> {
            let mut balls = Vec::new();
            for e in es {
                balls.push(pt(f(e)?, prec)?);
            }
            ok(
                "BallArray::from_balls",
                BallArray::from_balls(shape(), balls),
            )
        };
        let a = build(&|e| pow2(e))?;
        let s = ok(&what, a.add(&build(&|e| npow2(e))?, prec))?;
        for (b, e) in s.as_slice().iter().zip(es) {
            want_contains(&format!("{what} 2^{e} - 2^{e}"), b, &ExactNum::from(0))?;
        }
        let a = build(&|e| pow2(e))?;
        let t = build(&|e| pow2(e.checked_sub(1000).ok_or("e-1000")?))?;
        let u = ok(&what, a.add(&t, prec))?;
        for (b, e) in u.as_slice().iter().zip(es) {
            let exact = pow2(e)?.add_full_prec(&pow2(e.checked_sub(1000).ok_or("e-1000")?)?);
            want_contains(&format!("{what} 2^{e} + 2^({e}-1000)"), b, &exact)?;
        }
        let d = ok(&what, build(&|e| pow2(e))?.sub(&build(&|e| pow2(e))?, prec))?;
        for (b, e) in d.as_slice().iter().zip(es) {
            want_contains(&format!("{what} sub 2^{e}"), b, &ExactNum::from(0))?;
        }
    }
    Ok(())
}

/// Norms near EMAX (16.S6.6): `v = (2^(EMAX-1), 2^(EMAX-1), -2^(EMAX-1))`. norm_inf = 2^(EMAX-1), norm_1 = 3 * 2^(EMAX-1);
/// each must be an upper bound and at most V (1 + 2^-56). Mag holds these (its exponent is all of i64).
#[test]
fn cosmo_norms_near_emax() -> Result<(), String> {
    let b = EXPONENT_MAX.checked_sub(1).ok_or("emax-1")?;
    let v = vec_of(vec![pow2(b)?, pow2(b)?, npow2(b)?], 64)?;
    let ninf = v.norm_inf();
    want_ge("norm_inf", &ninf, &mag(TOP, EXPONENT_MAX)?)?;
    want_le("norm_inf tight", &ninf, &mag(TOP | 128, EXPONENT_MAX)?)?;
    let n1 = v.norm_1();
    want_ge(
        "norm_1",
        &n1,
        &mag(0xC000_0000_0000_0000, b.checked_add(2).ok_or("e")?)?,
    )?;
    want_le(
        "norm_1 tight",
        &n1,
        &mag(0xC000_0000_0000_00C0, b.checked_add(2).ok_or("e")?)?,
    )?;
    let m = mat_of(1, 3, vec![pow2(b)?, pow2(b)?, npow2(b)?], 64)?;
    want_ge(
        "mat norm_inf",
        &m.norm_inf(),
        &mag(0xC000_0000_0000_0000, b.checked_add(2).ok_or("e")?)?,
    )?;
    want_ge("mat norm_1", &m.norm_1(), &mag(TOP, EXPONENT_MAX)?)?;
    Ok(())
}

/// norm_2_up / norm_frob_up near EMAX: `(2^(EMAX-1), 2^(EMAX-1))` has 2-norm `sqrt(2) * 2^(EMAX-1)`; the squares
/// (2^(2 EMAX - 2), about 2^(2^63)) leave even the i64 Mag exponent range, but the norm itself does not.
/// 16.S6.6: an upper bound, at most V (1 + 2^-56). floor(sqrt(2) * 2^63) = 0xB504_F333_F9DE_6484.
#[test]
fn cosmo_norm2_near_emax() -> Result<(), String> {
    let b = EXPONENT_MAX.checked_sub(1).ok_or("emax-1")?;
    let lo = mag(0xB504_F333_F9DE_6484, EXPONENT_MAX)?;
    let hi = mag(0xB504_F333_F9DE_6585, EXPONENT_MAX)?;
    let v = vec_of(vec![pow2(b)?, npow2(b)?], 64)?;
    let n2 = v.norm_2_up();
    want_ge("norm_2_up", &n2, &lo)?;
    want_le("norm_2_up tight", &n2, &hi)?;
    let m = mat_of(1, 2, vec![pow2(b)?, npow2(b)?], 64)?;
    let nf = m.norm_frob_up();
    want_ge("norm_frob_up", &nf, &lo)?;
    want_le("norm_frob_up tight", &nf, &hi)
}

/// norm_2_up / norm_frob_up near EMIN: `(2^(EMIN+1), 2^(EMIN+1))` has 2-norm `sqrt(2) * 2^(EMIN+1)`; the squares fall below
/// the i64 Mag range, and an upper bound must still be at least the true norm (never 0).
#[test]
fn cosmo_norm2_near_emin() -> Result<(), String> {
    let s = EXPONENT_MIN.checked_add(1).ok_or("emin+1")?;
    let e = s.checked_add(1).ok_or("e")?;
    let lo = mag(0xB504_F333_F9DE_6484, e)?;
    let hi = mag(0xB504_F333_F9DE_6585, e)?;
    let v = vec_of(vec![pow2(s)?, npow2(s)?], 64)?;
    let n2 = v.norm_2_up();
    want_ge("norm_2_up", &n2, &lo)?;
    want_le("norm_2_up tight", &n2, &hi)?;
    let m = mat_of(2, 1, vec![pow2(s)?, pow2(s)?], 64)?;
    let nf = m.norm_frob_up();
    want_ge("norm_frob_up", &nf, &lo)?;
    want_le("norm_frob_up tight", &nf, &hi)
}

/// Precision edges on a huge-spread dot (16.S6.3): `MAX_PRECISION + 1` is InvalidArgument before any work.
#[test]
fn cosmo_dot_prec_edge() -> Result<(), String> {
    let x = vec_of(vec![pow2(1 << 61)?, ExactNum::from(1)], 64)?;
    let over = MAX_PRECISION.checked_add(1).ok_or("max+1")?;
    match x.dot(&x, over) {
        Err(BoundError::InvalidArgument { .. }) => Ok(()),
        other => Err(format!(
            "dot at MAX_PRECISION + 1: want InvalidArgument, got {other:?}"
        )),
    }
}
