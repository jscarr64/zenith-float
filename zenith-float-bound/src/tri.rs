//! Three-valued comparison result (CONTRACT 15.2, O15-W, D14).

/// Result of a ball comparison; `True`/`False` only when decided exactly.
///
/// The enum is deliberately **not** `#[non_exhaustive]` (API freeze F2): three values are the
/// complete set, so callers may match it exhaustively without a wildcard arm. Adding a
/// variant would be a breaking change.
#[derive(Debug, PartialEq, Eq)]
pub enum Tri {
    /// Decided true.
    True,
    /// Decided false.
    False,
    /// Not decided.
    Unknown,
}
