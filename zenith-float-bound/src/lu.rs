//! Slice LA2: LU factorisation with partial pivoting (`BallLu`), the determinant, the general log-determinant and the
//! linear solve; a certified singular matrix is `BoundError::Singular` (CONTRACT 16.LA2, 16.LA1.4).
//!
//! Owner: slice LA2 (Cassie). Only that slice edits this file. Contract text: `docs/CONTRACT.md` section "Slice LA2" and
//! 16.LA1.4/16.LA1.6 (general `logdet`); harness rows: `tests/harness/slice_lu.rs`; goldens: `tests/goldens/slices/lu/`;
//! mutants: `tools/mutants/manifest.d/lu.txt`.
//!
//! Every op works at `w = prec + 64 + 2·⌈log2(n + 1)⌉` (capped at `MAX_PRECISION`) and rounds outward to `prec` once at the
//! end, with no retry (16.LA.6).

use crate::array::{BallMat, BallVec};
use crate::ball::Ball;
use crate::enclosure::Enclosure;
use crate::error::{BoundError, BoundResult};
use crate::matrix::{check_prec, mismatch, reserve, resource};
use crate::precision::{Precision, MAX_PRECISION};
use crate::tri::Tri;
use zenith_float_num::ExactNum;

/// Fixed guard bits of 16.LA.6 (`64`), before `2·⌈log2(n + 1)⌉`.
const GUARD_FIXED: Precision = 64;

/// Sign, determinant and `ln|det|` of a square ball matrix (CONTRACT 16.LA1.1, 16.LA1.4).
#[derive(Debug)]
pub struct LogDet {
    sign: i8,
    det: Ball,
    ln_abs_det: Ball,
}

/// `P·A = L·U` with unit lower `L` (CONTRACT 16.LA2.3).
#[derive(Debug)]
pub struct BallLu {
    perm: Box<[usize]>,
    l: BallMat,
    u: BallMat,
}

/// `w = prec + 64 + 2·⌈log2(n + 1)⌉`, capped at `MAX_PRECISION` (16.LA.6). `⌈log2(n + 1)⌉` is the bit length of `n`.
fn work_prec(prec: Precision, n: usize) -> Precision {
    let bits = usize::BITS.checked_sub(n.leading_zeros());
    let guard = bits
        .and_then(|b| Precision::try_from(b).ok())
        .and_then(|b| b.checked_mul(2))
        .and_then(|b| b.checked_add(GUARD_FIXED));
    match guard.and_then(|g| prec.checked_add(g)) {
        Some(w) => w.min(MAX_PRECISION),
        None => MAX_PRECISION,
    }
}

/// The exact ball `0 ± 0` at `prec`.
fn zero(prec: Precision) -> BoundResult<Ball> {
    Ball::exact(ExactNum::from(0i64), prec)
}

/// `x` widened to `prec` (`0 + x` by S1 `add`; exact when the midpoint fits in `prec` bits).
fn widen(x: &Ball, prec: Precision) -> BoundResult<Ball> {
    zero(prec)?.add(x, prec)
}

/// The one outward rounding of a result to `prec` (16.LA.6).
fn round(x: Ball, prec: Precision) -> BoundResult<Ball> {
    x.add(&zero(prec)?, prec)
}

fn is_exact_zero(x: &Ball) -> bool {
    x.mid().is_zero() && x.rad().is_zero()
}

/// Whether the ball contains 0 (an undecided `Tri` counts as containing it).
fn has_zero(x: &Ball) -> bool {
    match x.contains(&ExactNum::from(0i64)) {
        Tri::True | Tri::Unknown => true,
        Tri::False => false,
    }
}

/// Row-major index `i·n + j`, checked.
fn at(op: &'static str, i: usize, n: usize, j: usize) -> BoundResult<usize> {
    i.checked_mul(n)
        .and_then(|r| r.checked_add(j))
        .ok_or_else(|| resource(op, "index overflow"))
}

fn entry<'v>(op: &'static str, v: &'v [Ball], k: usize) -> BoundResult<&'v Ball> {
    v.get(k).ok_or_else(|| resource(op, "index out of range"))
}

fn slot<'v>(op: &'static str, v: &'v mut [Ball], k: usize) -> BoundResult<&'v mut Ball> {
    v.get_mut(k)
        .ok_or_else(|| resource(op, "index out of range"))
}

fn diag<'m>(op: &'static str, m: &'m BallMat, k: usize) -> BoundResult<&'m Ball> {
    m.get(k, k)
        .ok_or_else(|| resource(op, "index out of range"))
}

/// Whether candidate `c` beats the current pivot `best`: a strictly larger `|mid|` (so ties keep the lowest index,
/// 16.LA2.1); among equal zero midpoints a ball that is not exactly `0 ± 0` beats an exact zero, so an exact-zero pivot is
/// chosen only when the whole column is exactly zero (the skip of 16.LA2.3).
fn beats(c: &Ball, best: &Ball) -> bool {
    match (c.mid().is_zero(), best.mid().is_zero()) {
        (true, true) => is_exact_zero(best) && !is_exact_zero(c),
        (true, false) => false,
        (false, true) => true,
        (false, false) => matches!(c.mid().abs_cmp(best.mid()), Some(s) if s > 0),
    }
}

/// Classifies every `u_kk` in index order: exactly `0 ± 0` → `Singular "singular"`, containing 0 → `Resource undecided`.
fn classify_pivots(op: &'static str, u: &BallMat, undecided: &'static str) -> BoundResult<()> {
    for k in 0..u.nrows() {
        let ukk = diag(op, u, k)?;
        if is_exact_zero(ukk) {
            return Err(BoundError::Singular {
                op,
                detail: "singular",
            });
        }
        if has_zero(ukk) {
            return Err(resource(op, undecided));
        }
    }
    Ok(())
}

/// `Σ_k ln|u_kk|` (16.LA1.4). Before S3 there is no `ln`: `Resource "needs_ln"` (16.LA1.6).
fn sum_ln_abs(op: &'static str, _lu: &BallLu, _prec: Precision) -> BoundResult<Ball> {
    Err(resource(op, "needs_ln"))
}

impl BallMat {
    /// LU factorisation with partial pivoting on midpoints (CONTRACT 16.LA2.3): `P·A = L·U`, `L` unit lower.
    /// Never refuses a singular matrix: an exactly zero column is skipped (`u_kk = 0 ± 0`), a pivot containing 0 gives
    /// all-reals multipliers.
    ///
    /// # Errors
    /// `InvalidArgument` for `prec`, `DimensionMismatch "rectangular"`, `Resource` for storage (16.LA.3).
    pub fn lu(&self, prec: Precision) -> BoundResult<BallLu> {
        const OP: &str = "BallMat::lu";
        check_prec(OP, prec)?;
        let n = self.nrows();
        if self.ncols() != n {
            return Err(mismatch(OP, "rectangular"));
        }
        let w = work_prec(prec, n);
        let nn = n
            .checked_mul(n)
            .ok_or_else(|| resource(OP, "index overflow"))?;
        let mut a: Vec<Ball> = reserve(OP, nn)?;
        for x in self.as_array().as_slice() {
            a.push(widen(x, w)?);
        }
        let mut l: Vec<Ball> = reserve(OP, nn)?;
        for _ in 0..nn {
            l.push(zero(w)?);
        }
        let mut perm: Vec<usize> = reserve(OP, n)?;
        perm.extend(0..n);
        for k in 0..n {
            let mut p = k;
            for i in (k..n).skip(1) {
                if beats(
                    entry(OP, &a, at(OP, i, n, k)?)?,
                    entry(OP, &a, at(OP, p, n, k)?)?,
                ) {
                    p = i;
                }
            }
            let kk = at(OP, k, n, k)?;
            if is_exact_zero(entry(OP, &a, at(OP, p, n, k)?)?) {
                continue; // the whole column is exactly zero: no swap, no elimination (16.LA2.3)
            }
            if p != k {
                for j in 0..n {
                    a.swap(at(OP, k, n, j)?, at(OP, p, n, j)?);
                }
                for j in 0..k {
                    l.swap(at(OP, k, n, j)?, at(OP, p, n, j)?);
                }
                perm.swap(k, p);
            }
            let pivot = widen(entry(OP, &a, kk)?, w)?;
            for i in (k..n).skip(1) {
                let ik = at(OP, i, n, k)?;
                let aik = core::mem::replace(slot(OP, &mut a, ik)?, zero(w)?);
                let lik = aik.div(&pivot, w)?;
                for j in (k..n).skip(1) {
                    let prod = widen(entry(OP, &a, at(OP, k, n, j)?)?, w)?.mul(&lik, w)?;
                    let s = slot(OP, &mut a, at(OP, i, n, j)?)?;
                    let aij = core::mem::replace(s, zero(w)?);
                    *s = aij.sub(&prod, w)?;
                }
                *slot(OP, &mut l, ik)? = lik;
            }
        }
        let mut lo: Vec<Ball> = reserve(OP, nn)?;
        for (idx, x) in l.into_iter().enumerate() {
            if idx.checked_rem(
                n.checked_add(1)
                    .ok_or_else(|| resource(OP, "index overflow"))?,
            ) == Some(0)
            {
                lo.push(Ball::exact(ExactNum::from(1i64), prec)?);
            } else {
                lo.push(round(x, prec)?);
            }
        }
        let mut up: Vec<Ball> = reserve(OP, nn)?;
        for x in a {
            up.push(round(x, prec)?);
        }
        Ok(BallLu {
            perm: perm.into_boxed_slice(),
            l: BallMat::from_balls(n, n, lo)?,
            u: BallMat::from_balls(n, n, up)?,
        })
    }

    /// General log-determinant backed by `lu` (CONTRACT 16.LA1.1, 16.LA1.4): exact `sign`, `det` and `ln|det|`.
    ///
    /// # Errors
    /// `InvalidArgument` for `prec`, `DimensionMismatch "rectangular"`, `Resource` for storage; `Singular "singular"` for a
    /// skipped zero column, `Resource "undecided_det"` for a `u_kk` containing 0; `Resource "needs_ln"` until S3 (16.LA1.6).
    pub fn logdet(&self, prec: Precision) -> BoundResult<LogDet> {
        const OP: &str = "BallMat::logdet";
        check_prec(OP, prec)?;
        if self.ncols() != self.nrows() {
            return Err(mismatch(OP, "rectangular"));
        }
        let lu = self.lu(prec)?;
        classify_pivots(OP, &lu.u, "undecided_det")?;
        let mut sign = lu.perm_sign();
        for k in 0..lu.u.nrows() {
            if diag(OP, &lu.u, k)?.mid().is_negative() {
                sign = sign.checked_neg().ok_or_else(|| resource(OP, "sign"))?;
            }
        }
        let det = lu.det(prec)?;
        let ln_abs_det = sum_ln_abs(OP, &lu, prec)?;
        Ok(LogDet {
            sign,
            det,
            ln_abs_det,
        })
    }

    /// `A⁻¹b` by `lu` then `BallLu::solve` (CONTRACT 16.LA2.3).
    ///
    /// # Errors
    /// As `lu`, then `DimensionMismatch "b_len"`, then as `BallLu::solve`.
    pub fn solve(&self, b: &BallVec, prec: Precision) -> BoundResult<BallVec> {
        const OP: &str = "BallMat::solve";
        check_prec(OP, prec)?;
        if self.ncols() != self.nrows() {
            return Err(mismatch(OP, "rectangular"));
        }
        if b.len() != self.nrows() {
            return Err(mismatch(OP, "b_len"));
        }
        self.lu(prec)?.solve(b, prec)
    }

    /// Many-RHS form of `solve`, column by column (CONTRACT 16.LA2.3).
    ///
    /// # Errors
    /// As `solve`; `DimensionMismatch "b_len"` when `b` has the wrong row count.
    pub fn solve_mat(&self, b: &BallMat, prec: Precision) -> BoundResult<BallMat> {
        const OP: &str = "BallMat::solve_mat";
        check_prec(OP, prec)?;
        let n = self.nrows();
        if self.ncols() != n {
            return Err(mismatch(OP, "rectangular"));
        }
        if b.nrows() != n {
            return Err(mismatch(OP, "b_len"));
        }
        let m = b.ncols();
        let total = n
            .checked_mul(m)
            .ok_or_else(|| resource(OP, "index overflow"))?;
        let mut cols: Vec<Ball> = reserve(OP, total)?;
        let lu = self.lu(prec)?;
        for j in 0..m {
            let mut col: Vec<Ball> = reserve(OP, n)?;
            for i in 0..n {
                let x = b
                    .get(i, j)
                    .ok_or_else(|| resource(OP, "index out of range"))?;
                col.push(widen(x, work_prec(prec, n))?);
            }
            let (_shape, xs) = lu
                .solve(&BallVec::from_balls(col), prec)?
                .into_array()
                .into_parts();
            cols.extend(xs);
        }
        let mut out: Vec<Ball> = reserve(OP, total)?;
        for i in 0..n {
            for j in 0..m {
                let s = slot(OP, &mut cols, at(OP, j, n, i)?)?;
                out.push(core::mem::replace(s, zero(prec)?));
            }
        }
        BallMat::from_balls(n, m, out)
    }
}

impl BallLu {
    /// The row permutation `P`: row `i` of `P·A` is row `perm[i]` of `A` (CONTRACT 16.LA2.3).
    pub fn perm(&self) -> &[usize] {
        &self.perm
    }

    /// The unit lower factor `L` (CONTRACT 16.LA2.3).
    pub fn l(&self) -> &BallMat {
        &self.l
    }

    /// The upper factor `U` (CONTRACT 16.LA2.3).
    pub fn u(&self) -> &BallMat {
        &self.u
    }

    /// Parity of `perm`: `1` even, `-1` odd, by inversion count (CONTRACT 16.LA2.3).
    pub fn perm_sign(&self) -> i8 {
        let mut odd = false;
        for (i, a) in self.perm.iter().enumerate() {
            for b in self.perm.iter().skip(i) {
                if a > b {
                    odd = !odd;
                }
            }
        }
        if odd {
            -1
        } else {
            1
        }
    }

    /// `perm_sign·Π u_kk` in index order (CONTRACT 16.LA2.3, 16.LA1.4): exactly `0 ± 0` if a column was skipped; an
    /// overflowing product is an all-reals `Ok`.
    ///
    /// # Errors
    /// `InvalidArgument` for `prec`.
    pub fn det(&self, prec: Precision) -> BoundResult<Ball> {
        const OP: &str = "BallLu::det";
        check_prec(OP, prec)?;
        let n = self.u.nrows();
        for k in 0..n {
            if is_exact_zero(diag(OP, &self.u, k)?) {
                return zero(prec);
            }
        }
        let w = work_prec(prec, n);
        let mut acc = Ball::exact(ExactNum::from(1i64), w)?;
        for k in 0..n {
            acc = acc.mul(diag(OP, &self.u, k)?, w)?;
        }
        if self.perm_sign() < 0 {
            acc = zero(w)?.sub(&acc, w)?;
        }
        round(acc, prec)
    }

    /// `A⁻¹b`: classify every `u_kk` (16.LA2.4), apply `P`, then `solve_lower` (unit), then `solve_upper` (16.S6.12).
    ///
    /// # Errors
    /// `InvalidArgument` for `prec`, `DimensionMismatch "b_len"`, `Resource` for storage; `Singular "singular"` for an
    /// exactly zero `u_kk`, `Resource "undecided_pivot"` for a `u_kk` containing 0; then the errors of 16.S6.12.
    pub fn solve(&self, b: &BallVec, prec: Precision) -> BoundResult<BallVec> {
        const OP: &str = "BallLu::solve";
        check_prec(OP, prec)?;
        let n = self.u.nrows();
        if b.len() != n {
            return Err(mismatch(OP, "b_len"));
        }
        classify_pivots(OP, &self.u, "undecided_pivot")?;
        let w = work_prec(prec, n);
        let mut pb: Vec<Ball> = reserve(OP, n)?;
        for src in self.perm.iter() {
            let x = b
                .get(*src)
                .ok_or_else(|| resource(OP, "index out of range"))?;
            pb.push(widen(x, w)?);
        }
        let y = self.l.solve_lower(&BallVec::from_balls(pb), true, w)?;
        let x = self.u.solve_upper(&y, false, w)?;
        let (_shape, xs) = x.into_array().into_parts();
        let mut out: Vec<Ball> = reserve(OP, n)?;
        for v in xs {
            out.push(round(v, prec)?);
        }
        Ok(BallVec::from_balls(out))
    }
}

impl LogDet {
    /// The exact sign of the determinant, `1` or `-1` (CONTRACT 16.LA1.1).
    pub fn sign(&self) -> i8 {
        self.sign
    }

    /// A ball enclosing the determinant (CONTRACT 16.LA1.1).
    pub fn det(&self) -> &Ball {
        &self.det
    }

    /// A ball enclosing `ln|det|` (CONTRACT 16.LA1.1).
    pub fn ln_abs_det(&self) -> &Ball {
        &self.ln_abs_det
    }
}

#[cfg(test)]
mod slice_lu_work_prec {
    use super::work_prec;
    use crate::precision::MAX_PRECISION;

    /// 16.LA.6: `w = prec + 64 + 2·⌈log2(n + 1)⌉`, capped.
    #[test]
    fn slice_lu_work_prec_formula_and_cap() -> Result<(), String> {
        // no-source: hand-derived: prec + 64 + 2·ceil(log2(n + 1)) per 16.LA.6, e.g. n = 3: 53 + 64 + 2·2 = 121.
        for (prec, n, want) in [
            (53, 0, 117),
            (53, 1, 119),
            (53, 3, 121),
            (53, 4, 123),
            (113, 7, 183),
            (192, 8, 264),
        ] {
            let got = work_prec(prec, n);
            if got != want {
                return Err(format!("work_prec({prec}, {n}) = {got}, want {want}"));
            }
        }
        if work_prec(MAX_PRECISION, 5) != MAX_PRECISION {
            return Err("work_prec must cap at MAX_PRECISION".to_string());
        }
        Ok(())
    }
}
