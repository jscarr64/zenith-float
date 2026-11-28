//! Domain policy and solver options (CONTRACT 10.1, 15.1, O6; options unit, Paul's rulings 2026-10-23).

/// How a restricted-domain operation treats input outside its domain (CONTRACT 10.2).
#[non_exhaustive]
#[derive(Debug, PartialEq, Eq, Default)]
pub enum DomainPolicy {
    /// `Domain` if any point of the input lies outside the domain.
    #[default]
    Strict,
    /// Enclose the function over the intersection of the input with the domain.
    Clip,
}

/// Options of the linear solvers (CONTRACT 10.1). Start from `SolveOpts::default()` and set fields.
#[non_exhaustive]
#[derive(Debug, PartialEq, Eq)]
pub struct SolveOpts {
    /// Iterative-refinement steps after the first solve (default 2).
    pub refine_steps: usize,
}

impl Default for SolveOpts {
    fn default() -> Self {
        SolveOpts { refine_steps: 2 }
    }
}

/// Options of the eigensolvers (CONTRACT 10.1). Start from `EigOpts::default()` and set fields.
#[non_exhaustive]
#[derive(Debug, PartialEq, Eq)]
pub struct EigOpts {
    /// Maximum number of Jacobi sweeps before giving up (default 64).
    pub max_sweeps: usize,
    /// Enclose clustered eigenvalues together instead of one by one (default true).
    pub cluster: bool,
}

impl Default for EigOpts {
    fn default() -> Self {
        EigOpts {
            max_sweeps: 64,
            cluster: true,
        }
    }
}
