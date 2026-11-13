//! Slice LA3: QR factorisation (`BallQr`) of a `BallMat` by Householder reflections, and least-squares solve.
//!
//! Owner: slice LA3 (Rusty). Only that slice edits this file. Contract: `docs/CONTRACT.md` section "Slice LA3"
//! (16.LA3, with 16.LA.2 to 16.LA.6). Harness rows: `tests/harness/slice_qr.rs`; goldens: `tests/goldens/slices/qr/`;
//! mutants: `tools/mutants/manifest.d/qr.txt` (1200 to 1212).
//!
//! Every step works at `w = la_work_prec(prec, m)` (16.LA.6) and reads its operands by borrow: each new entry is one
//! 16.S6.7 product entry (`dot_pairs`, `sub_pairs`), so no ball is copied. `qr` rounds `Q` and `R` outward to `prec` once
//! at the end; `lstsq` keeps the unrounded factors and hands `R` and `Qᵀb` to 16.S6.12 `solve_upper`.

use crate::array::{BallMat, BallVec};
use crate::ball::Ball;
use crate::error::{BoundError, BoundResult};
use crate::matrix::{
    check_prec, dot_pairs, excludes_zero, is_exact_zero, la_work_prec, mismatch, reserve, resource,
    round_out, sub_pairs,
};
use crate::policy::DomainPolicy;
use crate::precision::Precision;
use zenith_float_num::ExactNum;

/// Thin QR: `Q` is `m×p`, `R` is `p×p` upper with `diag(R) > 0` (CONTRACT 16.LA3.1).
#[derive(Debug)]
pub struct BallQr {
    q: BallMat,
    r: BallMat,
}

/// One Householder step (16.LA3.2): `v = x + s·σ·e_1` with `v_1` owned here and `v_i = x_i` (`i ≥ 2`) read back from
/// the step's input level, `vᵀv`, and whether `s = +1` (so the 16.LA3.2 sign fix negates row `k` of `R` and column `k`
/// of `Q`).
struct Reflector {
    v1: Ball,
    vtv: Ball,
    s_pos: bool,
}

/// Unrounded factors at the working precision, `Q` as `m×p` and `R` as `p×p`, both row-major, signs fixed.
struct Factors {
    q: Vec<Ball>,
    r: Vec<Ball>,
}

fn singular(op: &'static str, detail: &'static str) -> BoundError {
    BoundError::Singular { op, detail }
}

fn at<'a>(
    op: &'static str,
    data: &'a [Ball],
    i: usize,
    j: usize,
    cols: usize,
) -> BoundResult<&'a Ball> {
    i.checked_mul(cols)
        .and_then(|o| o.checked_add(j))
        .and_then(|o| data.get(o))
        .ok_or_else(|| resource(op, "index out of range"))
}

fn exact(v: i64, prec: Precision) -> BoundResult<Ball> {
    Ball::exact(ExactNum::from(v), prec)
}

fn size(op: &'static str, a: usize, b: usize) -> BoundResult<usize> {
    a.checked_mul(b)
        .ok_or_else(|| resource(op, "result storage cannot be reserved"))
}

/// `v_i` of a reflector: `v_1` for `i = 0`, else entry `(i, 0)` of the step's input level.
fn v_at<'a>(
    op: &'static str,
    v1: &'a Ball,
    level: &'a [Ball],
    i: usize,
    cols: usize,
) -> BoundResult<&'a Ball> {
    if i == 0 {
        Ok(v1)
    } else {
        at(op, level, i, 0, cols)
    }
}

/// Householder QR of the `m×p` row-major `a` at working precision `w` (16.LA3.2, 16.LA3.3), before the final rounding.
fn factor(op: &'static str, a: &[Ball], m: usize, p: usize, w: Precision) -> BoundResult<Factors> {
    let one = exact(1, w)?;
    let minus_one = exact(-1, w)?;
    let mut r = reserve::<Ball>(op, size(op, p, p)?)?;
    let mut hs = reserve::<Reflector>(op, p)?;
    // levels[k - 1] is the trailing (m-k)×(p-k) block after k steps; level 0 is `a` itself.
    let mut levels = reserve::<Vec<Ball>>(op, p)?;
    let mut r_rows = reserve::<Vec<Ball>>(op, p)?;
    for k in 0..p {
        let level: &[Ball] = match k.checked_sub(1) {
            None => a,
            Some(j) => levels
                .get(j)
                .map(Vec::as_slice)
                .ok_or_else(|| resource(op, "level out of range"))?,
        };
        let rows = m
            .checked_sub(k)
            .ok_or_else(|| resource(op, "step out of range"))?;
        let cols = p
            .checked_sub(k)
            .ok_or_else(|| resource(op, "step out of range"))?;
        let x = (0..rows)
            .map(|i| at(op, level, i, 0, cols))
            .collect::<BoundResult<Vec<&Ball>>>()?;
        if x.iter().all(|b| is_exact_zero(b)) {
            return Err(singular(op, "rank_deficient"));
        }
        let sigma2 = dot_pairs(op, x.iter().map(|b| (*b, *b)), rows, w)?;
        if !excludes_zero(&sigma2) {
            return Err(resource(op, "undecided_rank"));
        }
        let sigma = sigma2.sqrt(w, &DomainPolicy::Strict)?;
        if !excludes_zero(&sigma) {
            return Err(resource(op, "undecided_rank"));
        }
        let x1 = x.first().ok_or_else(|| resource(op, "empty column"))?;
        // s = sign(mid(x_1)), +1 at 0, chosen on the midpoint (16.LA3.2).
        let s_pos = !x1.mid().is_negative() || x1.mid().is_zero();
        let sign = if s_pos { &minus_one } else { &one };
        // v_1 = x_1 + s·σ = x_1 − σ·(−s).
        let v1 = sub_pairs(op, x1, std::iter::once((&sigma, sign)), 1, w)?;
        let v = (0..rows)
            .map(|i| v_at(op, &v1, level, i, cols))
            .collect::<BoundResult<Vec<&Ball>>>()?;
        let vtv = dot_pairs(op, v.iter().map(|b| (*b, *b)), rows, w)?;
        // f_j = 2·(vᵀa_j)/(vᵀv) for the trailing columns j ≥ 1.
        let mut f = reserve::<Ball>(op, cols)?;
        for j in 1..cols {
            let col = (0..rows)
                .map(|i| at(op, level, i, j, cols))
                .collect::<BoundResult<Vec<&Ball>>>()?;
            let tau = dot_pairs(op, v.iter().zip(col.iter()).map(|(a, c)| (*a, *c)), rows, w)?;
            f.push(tau.mul_2exp(1, w)?.div(&vtv, w)?);
        }
        // Row k of R: σ on the diagonal (−s·σ after the sign fix), then a_kj − f_j·v_1, negated when s = +1.
        let mut row = reserve::<Ball>(op, cols)?;
        for (j, fj) in (1..cols).zip(f.iter()) {
            let e = sub_pairs(
                op,
                at(op, level, 0, j, cols)?,
                std::iter::once((fj, &v1)),
                1,
                w,
            )?;
            row.push(if s_pos { e.negated() } else { e });
        }
        r_rows.push(row);
        let mut next = reserve::<Ball>(
            op,
            size(
                op,
                rows.checked_sub(1)
                    .ok_or_else(|| resource(op, "step out of range"))?,
                cols.checked_sub(1)
                    .ok_or_else(|| resource(op, "step out of range"))?,
            )?,
        )?;
        for i in 1..rows {
            let vi = at(op, level, i, 0, cols)?;
            for (j, fj) in (1..cols).zip(f.iter()) {
                next.push(sub_pairs(
                    op,
                    at(op, level, i, j, cols)?,
                    std::iter::once((fj, vi)),
                    1,
                    w,
                )?);
            }
        }
        drop(v);
        drop(x);
        hs.push(Reflector { v1, vtv, s_pos });
        levels.push(next);
        r.push(sigma);
    }
    // Assemble R row-major: zeros below the diagonal, the step diagonals, then the step rows.
    let mut diag = r.into_iter();
    let mut r_full = reserve::<Ball>(op, size(op, p, p)?)?;
    for (k, row) in r_rows.into_iter().enumerate() {
        for _ in 0..k {
            r_full.push(exact(0, w)?);
        }
        r_full.push(
            diag.next()
                .ok_or_else(|| resource(op, "diagonal out of range"))?,
        );
        r_full.extend(row);
    }
    let q = accumulate_q(op, a, &hs, &levels, m, p, w)?;
    Ok(Factors { q, r: r_full })
}

/// `Q = H_1·…·H_p·[e_1 … e_p]`, applying `H_p` first (16.LA3.2), then the sign fix on the columns of `Q`.
fn accumulate_q(
    op: &'static str,
    a: &[Ball],
    hs: &[Reflector],
    levels: &[Vec<Ball>],
    m: usize,
    p: usize,
    w: Precision,
) -> BoundResult<Vec<Ball>> {
    let mut y = reserve::<Ball>(op, size(op, m, p)?)?;
    for i in 0..m {
        for c in 0..p {
            y.push(exact(i64::from(i == c), w)?);
        }
    }
    for k in (0..p).rev() {
        let h = hs
            .get(k)
            .ok_or_else(|| resource(op, "reflector out of range"))?;
        let level: &[Ball] = match k.checked_sub(1) {
            None => a,
            Some(j) => levels
                .get(j)
                .map(Vec::as_slice)
                .ok_or_else(|| resource(op, "level out of range"))?,
        };
        let rows = m
            .checked_sub(k)
            .ok_or_else(|| resource(op, "step out of range"))?;
        let cols = p
            .checked_sub(k)
            .ok_or_else(|| resource(op, "step out of range"))?;
        let v = (0..rows)
            .map(|i| v_at(op, &h.v1, level, i, cols))
            .collect::<BoundResult<Vec<&Ball>>>()?;
        let mut f = reserve::<Ball>(op, p)?;
        for c in 0..p {
            let col = (k..m)
                .map(|i| at(op, &y, i, c, p))
                .collect::<BoundResult<Vec<&Ball>>>()?;
            let tau = dot_pairs(op, v.iter().zip(col.iter()).map(|(a, c)| (*a, *c)), rows, w)?;
            f.push(tau.mul_2exp(1, w)?.div(&h.vtv, w)?);
        }
        let mut upd = reserve::<Ball>(op, size(op, rows, p)?)?;
        for (i, vi) in (k..m).zip(v.iter()) {
            for (c, fc) in f.iter().enumerate() {
                upd.push(sub_pairs(
                    op,
                    at(op, &y, i, c, p)?,
                    std::iter::once((fc, *vi)),
                    1,
                    w,
                )?);
            }
        }
        drop(v);
        let keep = size(op, k, p)?;
        let mut fresh = upd.into_iter();
        let mut out = reserve::<Ball>(op, size(op, m, p)?)?;
        for (idx, old) in y.into_iter().enumerate() {
            if idx < keep {
                out.push(old);
            } else {
                out.push(
                    fresh
                        .next()
                        .ok_or_else(|| resource(op, "update out of range"))?,
                );
            }
        }
        y = out;
    }
    let mut out = reserve::<Ball>(op, size(op, m, p)?)?;
    for (idx, b) in y.into_iter().enumerate() {
        let c = idx
            .checked_rem(p)
            .ok_or_else(|| resource(op, "column out of range"))?;
        let flip = hs.get(c).is_some_and(|h| h.s_pos);
        out.push(if flip { b.negated() } else { b });
    }
    Ok(out)
}

/// 16.LA.2 computed triangles: `R`'s diagonal is classified before 16.S6.12 sees it.
fn check_r_diag(op: &'static str, r: &[Ball], p: usize) -> BoundResult<()> {
    for k in 0..p {
        let d = at(op, r, k, k, p)?;
        if is_exact_zero(d) {
            return Err(singular(op, "rank_deficient"));
        }
        if !excludes_zero(d) {
            return Err(resource(op, "undecided_pivot"));
        }
    }
    Ok(())
}

/// `solve_upper(R, Qᵀb)` with `Q` (`m×p`) and `R` (`p×p`) row-major at working precision `w` (16.LA3.2).
fn solve_from(
    op: &'static str,
    q: &[Ball],
    r: &BallMat,
    b: &BallVec,
    w: Precision,
    prec: Precision,
) -> BoundResult<BallVec> {
    let (m, p) = (b.len(), r.nrows());
    check_r_diag(op, r.as_array().as_slice(), p)?;
    let bs = b.as_array().as_slice();
    let mut qtb = reserve::<Ball>(op, p)?;
    for c in 0..p {
        let col = (0..m)
            .map(|i| at(op, q, i, c, p))
            .collect::<BoundResult<Vec<&Ball>>>()?;
        qtb.push(dot_pairs(
            op,
            col.iter().zip(bs.iter()).map(|(c, b)| (*c, b)),
            m,
            w,
        )?);
    }
    r.solve_upper(&BallVec::from_balls(qtb), false, prec)
}

fn into_mat(rows: usize, cols: usize, v: Vec<Ball>) -> BoundResult<BallMat> {
    BallMat::from_balls(rows, cols, v)
}

impl BallMat {
    /// Householder QR (CONTRACT 16.LA3.1 to 16.LA3.4): thin `Q` (`m×p`) and upper `R` (`p×p`) with `diag(R) > 0`,
    /// computed at `w = prec + 64 + 2·⌈log2(m + 1)⌉` and rounded outward to `prec` once.
    ///
    /// # Errors
    /// `InvalidArgument` for a bad `prec`; `DimensionMismatch "wide"` when `p > m`; `Resource` for storage;
    /// `Singular "rank_deficient"` for an exactly zero trailing column; `Resource "undecided_rank"` when `σ` reaches 0.
    pub fn qr(&self, prec: Precision) -> BoundResult<BallQr> {
        const OP: &str = "BallMat::qr";
        check_prec(OP, prec)?;
        let (m, p) = (self.nrows(), self.ncols());
        if p > m {
            return Err(mismatch(OP, "wide"));
        }
        let w = la_work_prec(prec, m);
        let f = factor(OP, self.as_array().as_slice(), m, p, w)?;
        let q =
            f.q.into_iter()
                .map(|b| round_out(OP, b, prec))
                .collect::<BoundResult<Vec<Ball>>>()?;
        let r =
            f.r.into_iter()
                .map(|b| round_out(OP, b, prec))
                .collect::<BoundResult<Vec<Ball>>>()?;
        Ok(BallQr {
            q: into_mat(m, p, q)?,
            r: into_mat(p, p, r)?,
        })
    }

    /// `argmin ‖A·x − b‖₂` = `solve_upper(R, Qᵀb)` (CONTRACT 16.LA3.1, 16.LA3.2); an `m×0` system gives the empty `x`.
    ///
    /// # Errors
    /// As `qr`, then `DimensionMismatch "b_len"` when `b.len() ≠ m`, `Resource "undecided_pivot"` for an `R` diagonal
    /// reaching 0 (16.LA.2), and the errors of 16.S6.12 `solve_upper`.
    pub fn lstsq(&self, b: &BallVec, prec: Precision) -> BoundResult<BallVec> {
        const OP: &str = "BallMat::lstsq";
        check_prec(OP, prec)?;
        let (m, p) = (self.nrows(), self.ncols());
        if p > m {
            return Err(mismatch(OP, "wide"));
        }
        if b.len() != m {
            return Err(mismatch(OP, "b_len"));
        }
        if p == 0 {
            return Ok(BallVec::from_balls(Vec::new()));
        }
        let w = la_work_prec(prec, m);
        let f = factor(OP, self.as_array().as_slice(), m, p, w)?;
        solve_from(OP, &f.q, &into_mat(p, p, f.r)?, b, w, prec)
    }
}

impl BallQr {
    /// The factor `Q` (`m×p`, CONTRACT 16.LA3.1).
    pub fn q(&self) -> &BallMat {
        &self.q
    }

    /// The factor `R` (`p×p` upper, `diag(R) > 0`, CONTRACT 16.LA3.1).
    pub fn r(&self) -> &BallMat {
        &self.r
    }

    /// `(Q, R)` by value (CONTRACT 16.LA3.1).
    pub fn into_parts(self) -> (BallMat, BallMat) {
        (self.q, self.r)
    }

    /// `solve_upper(R, Qᵀb)` from the stored factors (CONTRACT 16.LA3.1); an `m×0` factorisation gives the empty `x`.
    ///
    /// # Errors
    /// `InvalidArgument` for a bad `prec`; `DimensionMismatch "b_len"` when `b.len() ≠ m`; `Singular "rank_deficient"`
    /// or `Resource "undecided_pivot"` from the `R` diagonal (16.LA.2); the errors of 16.S6.12 `solve_upper`.
    pub fn lstsq(&self, b: &BallVec, prec: Precision) -> BoundResult<BallVec> {
        const OP: &str = "BallQr::lstsq";
        check_prec(OP, prec)?;
        let (m, p) = (self.q.nrows(), self.q.ncols());
        if b.len() != m {
            return Err(mismatch(OP, "b_len"));
        }
        if p == 0 {
            return Ok(BallVec::from_balls(Vec::new()));
        }
        solve_from(
            OP,
            self.q.as_array().as_slice(),
            &self.r,
            b,
            la_work_prec(prec, m),
            prec,
        )
    }
}
