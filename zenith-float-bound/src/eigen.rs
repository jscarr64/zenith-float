//! Slice LA5: Eigenvalue enclosures of a symmetric `BallMat` (Jacobi sweeps under `EigOpts::max_sweeps`); running out of sweeps is `BoundError::NotConverged`, never a loose answer.
//!
//! Owner: slice LA5 (owner to be assigned). Only that slice edits this file.
//! S0 stub (2026-10-24): no items yet. The slice adds its code here, its contract text in
//! `docs/CONTRACT.md` section "Slice LA5", its harness rows in `tests/harness/slice_eigen.rs`, its
//! goldens under `tests/goldens/slices/eigen/` and its mutants in `tools/mutants/manifest.d/eigen.txt`.

use crate::array::{BallMat, BallVec};
use crate::cball::CBall;
use crate::error::{BoundError, BoundResult};
use crate::policy::EigOpts;
use crate::precision::Precision;

/// The skeleton refusal of an op whose algorithm has not landed yet.
fn not_implemented(op: &'static str) -> BoundError {
    BoundError::Resource {
        op,
        detail: "not_implemented",
    }
}

impl BallMat {
    /// Ascending eigenvalue enclosures of a symmetric ball matrix (CONTRACT 16.LA5a.1).
    /// Not implemented yet: returns Resource "not_implemented".
    ///
    /// # Errors
    /// `Resource` "not_implemented".
    pub fn eig_sym(&self, _opts: &EigOpts, _prec: Precision) -> BoundResult<BallVec> {
        Err(not_implemented("BallMat::eig_sym"))
    }

    /// Eigenvalue enclosures of a general ball matrix, with multiplicity (CONTRACT 16.LA5b.1).
    /// Not implemented yet: returns Resource "not_implemented".
    ///
    /// # Errors
    /// `Resource` "not_implemented".
    pub fn eig(&self, _opts: &EigOpts, _prec: Precision) -> BoundResult<Box<[CBall]>> {
        Err(not_implemented("BallMat::eig"))
    }
}
