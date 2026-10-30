//! Error type of zenith-bound (variant set: CONTRACT O11-W; struct form and the two
//! linear-algebra variants: 2026-10-23 rulings a and b; DimensionMismatch: S0, 2026-10-24; NotSymmetric: LA plan, 2026-10-24).

use core::fmt;

/// Why an operation returned no value (CONTRACT 6.4, O11-W).
///
/// Every variant carries `op`, the operation that refused, and `detail`, a short reason.
/// The enum is `#[non_exhaustive]`: new variants may be added in a minor release.
#[derive(Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum BoundError {
    /// The value is not representable (a `Mag::finite` rejection; a NaN or infinite midpoint).
    Representation {
        /// The operation that refused.
        op: &'static str,
        /// A short reason.
        detail: &'static str,
    },
    /// An argument is invalid (`prec == 0`; a container shape that does not fit its data). Two
    /// operands whose shapes disagree are [`BoundError::DimensionMismatch`].
    InvalidArgument {
        /// The operation that refused.
        op: &'static str,
        /// A short reason.
        detail: &'static str,
    },
    /// The input lies outside the operation's domain (CONTRACT 4.2, 10.2).
    Domain {
        /// The operation that refused.
        op: &'static str,
        /// A short reason.
        detail: &'static str,
    },
    /// The engine gave up (budget, iteration, precision or memory cap; "cannot certify at this prec").
    Resource {
        /// The operation that refused.
        op: &'static str,
        /// A short reason.
        detail: &'static str,
    },
    /// The matrix is certified not positive definite (ruling b).
    NotPositiveDefinite {
        /// The operation that refused.
        op: &'static str,
        /// A short reason.
        detail: &'static str,
    },
    /// The matrix is certified singular (ruling b).
    Singular {
        /// The operation that refused.
        op: &'static str,
        /// A short reason.
        detail: &'static str,
    },
    /// An iteration ran out of its budget (e.g. `EigOpts::max_sweeps`) before its stopping test held
    /// (Paul, 2026-10-23); never a silently loose answer.
    NotConverged {
        /// The operation that refused.
        op: &'static str,
        /// A short reason.
        detail: &'static str,
    },
    /// Operand shapes disagree (matmul inner dimensions, a non-square matrix where a square one is
    /// required, a right-hand side of the wrong length; S0, 2026-10-24).
    DimensionMismatch {
        /// The operation that refused.
        op: &'static str,
        /// A short reason.
        detail: &'static str,
    },
    /// A matrix required to be symmetric is not (symmetric eigenvalue enclosures, LA5a; LA plan, 2026-10-24).
    NotSymmetric {
        /// The operation that refused.
        op: &'static str,
        /// A short reason.
        detail: &'static str,
    },
}

impl BoundError {
    /// The operation that refused.
    pub fn op(&self) -> &'static str {
        match self {
            BoundError::Representation { op, .. } => op,
            BoundError::InvalidArgument { op, .. } => op,
            BoundError::Domain { op, .. } => op,
            BoundError::Resource { op, .. } => op,
            BoundError::NotPositiveDefinite { op, .. } => op,
            BoundError::Singular { op, .. } => op,
            BoundError::NotConverged { op, .. } => op,
            BoundError::DimensionMismatch { op, .. } => op,
            BoundError::NotSymmetric { op, .. } => op,
        }
    }

    /// The short reason given by the operation.
    pub fn detail(&self) -> &'static str {
        match self {
            BoundError::Representation { detail, .. } => detail,
            BoundError::InvalidArgument { detail, .. } => detail,
            BoundError::Domain { detail, .. } => detail,
            BoundError::Resource { detail, .. } => detail,
            BoundError::NotPositiveDefinite { detail, .. } => detail,
            BoundError::Singular { detail, .. } => detail,
            BoundError::NotConverged { detail, .. } => detail,
            BoundError::DimensionMismatch { detail, .. } => detail,
            BoundError::NotSymmetric { detail, .. } => detail,
        }
    }

    fn kind_text(&self) -> &'static str {
        match self {
            BoundError::Representation { .. } => "value not representable",
            BoundError::InvalidArgument { .. } => "invalid argument",
            BoundError::Domain { .. } => "input outside the domain",
            BoundError::Resource { .. } => "resource limit reached",
            BoundError::NotPositiveDefinite { .. } => "matrix not positive definite",
            BoundError::Singular { .. } => "matrix singular",
            BoundError::NotConverged { .. } => "iteration did not converge",
            BoundError::DimensionMismatch { .. } => "dimensions do not match",
            BoundError::NotSymmetric { .. } => "matrix not symmetric",
        }
    }
}

impl fmt::Display for BoundError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} in {}: {}",
            self.kind_text(),
            self.op(),
            self.detail()
        )
    }
}

impl core::error::Error for BoundError {}

/// Result of a zenith-bound operation (CONTRACT 15.3).
pub type BoundResult<T> = Result<T, BoundError>;
