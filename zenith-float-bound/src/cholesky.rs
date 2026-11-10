//! Slice LA1: Cholesky factorisation of a symmetric positive definite `BallMat` (`BallCholesky`) and the log-determinant; a certified non-SPD pivot is `BoundError::NotPositiveDefinite`.
//!
//! Owner: slice LA1 (owner to be assigned). Only that slice edits this file.
//! S0 stub (2026-10-24): no items yet. The slice adds its code here, its contract text in
//! `docs/CONTRACT.md` section "Slice LA1", its harness rows in `tests/harness/slice_cholesky.rs`, its
//! goldens under `tests/goldens/slices/cholesky/` and its mutants in `tools/mutants/manifest.d/cholesky.txt`.

use crate::array::{BallMat, BallVec};
use crate::ball::Ball;
use crate::error::{BoundError, BoundResult};
use crate::precision::Precision;
use crate::tri::Tri;

/// The skeleton refusal of an op whose algorithm has not landed yet.
fn not_implemented(op: &'static str) -> BoundError {
    BoundError::Resource {
        op,
        detail: "not_implemented",
    }
}

/// Cholesky factor `G` (lower, `A = G·Gᵀ`) of an SPD ball matrix (CONTRACT 16.LA1.1).
/// Not implemented yet: `BallMat::cholesky` returns Resource "not_implemented", so no value exists.
#[derive(Debug)]
pub struct BallCholesky {
    g: BallMat,
}

/// `A = L·diag(d)·Lᵀ` with unit lower `L` (CONTRACT 16.LA1.1).
/// Not implemented yet: `BallMat::ldl` returns Resource "not_implemented", so no value exists.
#[derive(Debug)]
pub struct BallLdl {
    l: BallMat,
    d: BallVec,
}

impl BallMat {
    /// Cholesky factorisation (CONTRACT 16.LA1.1).
    /// Not implemented yet: returns Resource "not_implemented".
    ///
    /// # Errors
    /// `Resource` "not_implemented".
    pub fn cholesky(&self, _prec: Precision) -> BoundResult<BallCholesky> {
        Err(not_implemented("BallMat::cholesky"))
    }

    /// `LDLᵀ` factorisation (CONTRACT 16.LA1.1).
    /// Not implemented yet: returns Resource "not_implemented".
    ///
    /// # Errors
    /// `Resource` "not_implemented".
    pub fn ldl(&self, _prec: Precision) -> BoundResult<BallLdl> {
        Err(not_implemented("BallMat::ldl"))
    }

    /// Whether every symmetric matrix in the ball is SPD (CONTRACT 16.LA1.1).
    /// Not implemented yet: returns Resource "not_implemented".
    ///
    /// # Errors
    /// `Resource` "not_implemented".
    pub fn is_spd(&self, _prec: Precision) -> BoundResult<Tri> {
        Err(not_implemented("BallMat::is_spd"))
    }
}

impl BallCholesky {
    /// The lower factor `G` (CONTRACT 16.LA1.1).
    /// Not implemented yet: no `BallCholesky` can be built (`cholesky` returns Resource "not_implemented").
    pub fn g(&self) -> &BallMat {
        &self.g
    }

    /// The lower factor `G`, by value (CONTRACT 16.LA1.1).
    /// Not implemented yet: no `BallCholesky` can be built (`cholesky` returns Resource "not_implemented").
    pub fn into_g(self) -> BallMat {
        self.g
    }

    /// `A⁻¹b` through `G` (CONTRACT 16.LA1.1).
    /// Not implemented yet: returns Resource "not_implemented".
    ///
    /// # Errors
    /// `Resource` "not_implemented".
    pub fn solve(&self, _b: &BallVec, _prec: Precision) -> BoundResult<BallVec> {
        Err(not_implemented("BallCholesky::solve"))
    }

    /// `rᵀA⁻¹r = ‖G⁻¹r‖²` (CONTRACT 16.LA1.1).
    /// Not implemented yet: returns Resource "not_implemented".
    ///
    /// # Errors
    /// `Resource` "not_implemented".
    pub fn quad_form(&self, _r: &BallVec, _prec: Precision) -> BoundResult<Ball> {
        Err(not_implemented("BallCholesky::quad_form"))
    }

    /// `ln det A = 2·Σ ln g_kk` (CONTRACT 16.LA1.1, 16.LA1.4).
    /// Not implemented yet: returns Resource "not_implemented".
    ///
    /// # Errors
    /// `Resource` "not_implemented".
    pub fn logdet_spd(&self, _prec: Precision) -> BoundResult<Ball> {
        Err(not_implemented("BallCholesky::logdet_spd"))
    }
}

impl BallLdl {
    /// The unit lower factor `L` (CONTRACT 16.LA1.1).
    /// Not implemented yet: no `BallLdl` can be built (`ldl` returns Resource "not_implemented").
    pub fn l(&self) -> &BallMat {
        &self.l
    }

    /// The diagonal `d` (CONTRACT 16.LA1.1).
    /// Not implemented yet: no `BallLdl` can be built (`ldl` returns Resource "not_implemented").
    pub fn d(&self) -> &BallVec {
        &self.d
    }
}
