//! The equivalent-mutant ledger holds only real equivalents (Paul's no-dead-code ruling).
//!
//! A mutant of a branch that can never run survives as "equivalent"; the fix is to delete the
//! branch, not to document it. The `None => Some(Mag::zero())` arm of `sub_wide` (mutant 025-L12)
//! was such a branch and is removed, so neither the manifest nor the ledger may name it.

const EQUIVALENT: &str = include_str!("../../tools/mutants/EQUIVALENT.txt");
const MANIFEST: &str = include_str!("../../tools/mutants/manifest.txt");

/// Mutant ids of branches removed as dead code.
const REMOVED: [&str; 1] = ["025-L12"];

/// Ids of the non-comment rows of a ledger file (first TAB-separated column for EQUIVALENT.txt,
/// second for manifest.txt).
fn ids(text: &str, column: usize) -> Vec<&str> {
    text.lines()
        .map(str::trim_start)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .filter_map(|l| l.split('\t').nth(column))
        .collect()
}

#[test]
fn equivalent_ledger_names_no_removed_dead_branch() -> Result<(), String> {
    let left: Vec<&str> = ids(EQUIVALENT, 0)
        .into_iter()
        .filter(|id| REMOVED.contains(id))
        .collect();
    if left.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "EQUIVALENT.txt still lists removed dead branches: {:?}",
            left
        ))
    }
}

#[test]
fn manifest_has_no_mutant_of_a_removed_dead_branch() -> Result<(), String> {
    let left: Vec<&str> = ids(MANIFEST, 1)
        .into_iter()
        .filter(|id| REMOVED.contains(id))
        .collect();
    if left.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "manifest still mutates removed dead branches: {:?}",
            left
        ))
    }
}
