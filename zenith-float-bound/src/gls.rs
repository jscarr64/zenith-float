//! Slice LA4: Generalised least squares (`GlsResult`): estimates and their covariance from a design matrix and an SPD covariance, built on LA1 and LA3.
//!
//! Owner: slice LA4 (owner to be assigned). Only that slice edits this file.
//! S0 stub (2026-10-24): no items yet. The slice adds its code here, its contract text in
//! `docs/CONTRACT.md` section "Slice LA4", its harness rows in `tests/harness/slice_gls.rs`, its
//! goldens under `tests/goldens/slices/gls/` and its mutants in `tools/mutants/manifest.d/gls.txt`.

use crate::array::{BallMat, BallVec};
use crate::ball::Ball;
use crate::error::{BoundError, BoundResult};
use crate::precision::Precision;

/// Generalised least-squares result (CONTRACT 16.LA4.2).
/// Not implemented yet: `gls` returns Resource "not_implemented", so no value exists.
#[derive(Debug)]
pub struct GlsResult {
    beta: BallVec,
    cov: BallMat,
    chi2: Ball,
    logdet_c: Ball,
}

impl GlsResult {
    /// The estimates `β` (CONTRACT 16.LA4.2).
    /// Not implemented yet: no `GlsResult` can be built (`gls` returns Resource "not_implemented").
    pub fn beta(&self) -> &BallVec {
        &self.beta
    }

    /// The covariance of `β` (CONTRACT 16.LA4.2).
    /// Not implemented yet: no `GlsResult` can be built (`gls` returns Resource "not_implemented").
    pub fn cov(&self) -> &BallMat {
        &self.cov
    }

    /// The whitened residual sum of squares (the gold's `rss`) (CONTRACT 16.LA4.2).
    /// Not implemented yet: no `GlsResult` can be built (`gls` returns Resource "not_implemented").
    pub fn chi2(&self) -> &Ball {
        &self.chi2
    }

    /// `ln det C` (CONTRACT 16.LA4.2).
    /// Not implemented yet: no `GlsResult` can be built (`gls` returns Resource "not_implemented").
    pub fn logdet_c(&self) -> &Ball {
        &self.logdet_c
    }
}

/// GLS of `y` on the `m×p` design `x` with SPD `m×m` covariance `c` (CONTRACT 16.LA4.2, 16.LA4.4).
/// Not implemented yet: returns Resource "not_implemented".
///
/// # Errors
/// `Resource` "not_implemented".
pub fn gls(_x: &BallMat, _c: &BallMat, _y: &BallVec, _prec: Precision) -> BoundResult<GlsResult> {
    Err(BoundError::Resource {
        op: "gls",
        detail: "not_implemented",
    })
}
