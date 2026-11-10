//! Slice LA2: LU factorisation with partial pivoting (`BallLu`), the determinant and linear solve; a certified singular matrix is `BoundError::Singular`.
//!
//! Owner: slice LA2 (owner to be assigned). Only that slice edits this file.
//! S0 stub (2026-10-24): no items yet. The slice adds its code here, its contract text in
//! `docs/CONTRACT.md` section "Slice LA2", its harness rows in `tests/harness/slice_lu.rs`, its
//! goldens under `tests/goldens/slices/lu/` and its mutants in `tools/mutants/manifest.d/lu.txt`.

use crate::array::{BallMat, BallVec};
use crate::ball::Ball;
use crate::error::{BoundError, BoundResult};
use crate::precision::Precision;

/// The skeleton refusal of an op whose algorithm has not landed yet.
fn not_implemented(op: &'static str) -> BoundError {
    BoundError::Resource {
        op,
        detail: "not_implemented",
    }
}

/// Sign, determinant and `ln|det|` of a square ball matrix (CONTRACT 16.LA1.1, 16.LA1.4).
/// Not implemented yet: `BallMat::logdet` returns Resource "not_implemented", so no value exists.
#[derive(Debug)]
pub struct LogDet {
    sign: i8,
    det: Ball,
    ln_abs_det: Ball,
}

/// `P·A = L·U` with unit lower `L` (CONTRACT 16.LA2.3).
/// Not implemented yet: `BallMat::lu` returns Resource "not_implemented", so no value exists.
#[derive(Debug)]
pub struct BallLu {
    perm: Box<[usize]>,
    l: BallMat,
    u: BallMat,
}

impl BallMat {
    /// LU factorisation with partial pivoting on midpoints (CONTRACT 16.LA2.3).
    /// Not implemented yet: returns Resource "not_implemented".
    ///
    /// # Errors
    /// `Resource` "not_implemented".
    pub fn lu(&self, _prec: Precision) -> BoundResult<BallLu> {
        Err(not_implemented("BallMat::lu"))
    }

    /// General log-determinant, backed by LU (CONTRACT 16.LA1.1, 16.LA1.4).
    /// Not implemented yet: returns Resource "not_implemented".
    ///
    /// # Errors
    /// `Resource` "not_implemented".
    pub fn logdet(&self, _prec: Precision) -> BoundResult<LogDet> {
        Err(not_implemented("BallMat::logdet"))
    }

    /// `A⁻¹b` by `lu` then `BallLu::solve` (CONTRACT 16.LA2.3).
    /// Not implemented yet: returns Resource "not_implemented".
    ///
    /// # Errors
    /// `Resource` "not_implemented".
    pub fn solve(&self, _b: &BallVec, _prec: Precision) -> BoundResult<BallVec> {
        Err(not_implemented("BallMat::solve"))
    }

    /// Many-RHS form of `solve` (CONTRACT 16.LA2.3).
    /// Not implemented yet: returns Resource "not_implemented".
    ///
    /// # Errors
    /// `Resource` "not_implemented".
    pub fn solve_mat(&self, _b: &BallMat, _prec: Precision) -> BoundResult<BallMat> {
        Err(not_implemented("BallMat::solve_mat"))
    }
}

impl BallLu {
    /// The row permutation `P` (CONTRACT 16.LA2.3).
    /// Not implemented yet: no `BallLu` can be built (`lu` returns Resource "not_implemented").
    pub fn perm(&self) -> &[usize] {
        &self.perm
    }

    /// The unit lower factor `L` (CONTRACT 16.LA2.3).
    /// Not implemented yet: no `BallLu` can be built (`lu` returns Resource "not_implemented").
    pub fn l(&self) -> &BallMat {
        &self.l
    }

    /// The upper factor `U` (CONTRACT 16.LA2.3).
    /// Not implemented yet: no `BallLu` can be built (`lu` returns Resource "not_implemented").
    pub fn u(&self) -> &BallMat {
        &self.u
    }

    /// Parity of `perm`: `1` even, `-1` odd, by inversion count (CONTRACT 16.LA2.3).
    /// Not implemented yet: no `BallLu` can be built (`lu` returns Resource "not_implemented").
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

    /// `perm_sign·Π u_kk` (CONTRACT 16.LA2.3).
    /// Not implemented yet: returns Resource "not_implemented".
    ///
    /// # Errors
    /// `Resource` "not_implemented".
    pub fn det(&self, _prec: Precision) -> BoundResult<Ball> {
        Err(not_implemented("BallLu::det"))
    }

    /// `A⁻¹b`: apply `P`, then `solve_lower` (unit), then `solve_upper` (CONTRACT 16.LA2.3, 16.LA2.4).
    /// Not implemented yet: returns Resource "not_implemented".
    ///
    /// # Errors
    /// `Resource` "not_implemented".
    pub fn solve(&self, _b: &BallVec, _prec: Precision) -> BoundResult<BallVec> {
        Err(not_implemented("BallLu::solve"))
    }
}

impl LogDet {
    /// The exact sign of the determinant, `1` or `-1` (CONTRACT 16.LA1.1).
    /// Not implemented yet: no `LogDet` can be built (`logdet` returns Resource "not_implemented").
    pub fn sign(&self) -> i8 {
        self.sign
    }

    /// A ball enclosing the determinant (CONTRACT 16.LA1.1).
    /// Not implemented yet: no `LogDet` can be built (`logdet` returns Resource "not_implemented").
    pub fn det(&self) -> &Ball {
        &self.det
    }

    /// A ball enclosing `ln|det|` (CONTRACT 16.LA1.1).
    /// Not implemented yet: no `LogDet` can be built (`logdet` returns Resource "not_implemented").
    pub fn ln_abs_det(&self) -> &Ball {
        &self.ln_abs_det
    }
}
