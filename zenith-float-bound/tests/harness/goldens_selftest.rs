//! COVERS: common::goldens::golden_dir, common::goldens::repo_root
//! Every golden file exists, is non-empty and starts with its expected header line; every *.txt in tests/goldens has a table row and a
//! SHA256SUMS line and vice versa. (The sha256 values themselves are checked by tools/goldens_regen.sh, not here: no hashing crate.)

use crate::common::expect_eq;
use crate::common::goldens::{golden_dir, repo_root, GOLDEN_FILES};

#[test]
fn golden_dir_exists() -> Result<(), String> {
    if golden_dir().is_dir() {
        Ok(())
    } else {
        Err(format!("{:?} is not a directory", golden_dir()))
    }
}

/// The goldens are read from the checkout the test runs in (cargo test runs it from the package root), never from a path baked in at
/// compile time: a test binary built in one worktree and run in another must read the second worktree's files.
#[test]
fn goldens_come_from_the_checkout_the_test_runs_in() -> Result<(), String> {
    let canon =
        |p: std::path::PathBuf| std::fs::canonicalize(&p).map_err(|e| format!("{:?}: {}", p, e));
    let here = canon(std::env::current_dir().map_err(|e| format!("current_dir: {}", e))?)?;
    expect_eq(
        "repo_root() is the directory the test runs in",
        &canon(repo_root())?,
        &here,
    )?;
    expect_eq(
        "golden_dir() is <checkout>/tests/goldens",
        &canon(golden_dir())?,
        &here.join("tests").join("goldens"),
    )
}

#[test]
fn each_golden_file_exists_is_non_empty_and_starts_with_its_header() -> Result<(), String> {
    for g in GOLDEN_FILES.iter() {
        let path = golden_dir().join(g.name);
        let text = std::fs::read_to_string(&path).map_err(|e| format!("{:?}: {}", path, e))?;
        if text.is_empty() {
            return Err(format!("{} is empty", g.name));
        }
        let first = text.lines().next().unwrap_or("");
        expect_eq(&format!("{} header line", g.name), &first, &g.header)?;
        if !text.ends_with('\n') {
            return Err(format!("{} does not end with a newline", g.name));
        }
        if text.contains('\r') {
            return Err(format!("{} contains CR (LF line endings only)", g.name));
        }
    }
    Ok(())
}

#[test]
fn golden_table_matches_the_directory_and_sha256sums_list() -> Result<(), String> {
    let mut on_disk: Vec<String> = Vec::new();
    let entries = std::fs::read_dir(golden_dir()).map_err(|e| format!("read_dir: {}", e))?;
    for entry in entries {
        let entry = entry.map_err(|e| format!("dir entry: {}", e))?;
        let name = entry.file_name().to_string_lossy().into_owned();
        if name.ends_with(".txt") {
            on_disk.push(name);
        }
    }
    on_disk.sort();
    let mut table: Vec<String> = GOLDEN_FILES.iter().map(|g| String::from(g.name)).collect();
    table.sort();
    expect_eq("*.txt files vs GOLDEN_FILES", &on_disk, &table)?;

    let sums_path = golden_dir().join("SHA256SUMS");
    let sums =
        std::fs::read_to_string(&sums_path).map_err(|e| format!("{:?}: {}", sums_path, e))?;
    let mut listed: Vec<String> = Vec::new();
    for line in sums.lines() {
        let (hash, name) = line
            .split_once("  ")
            .ok_or_else(|| format!("SHA256SUMS line without two spaces: {:?}", line))?;
        let ok = hash.len() == 64
            && hash
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b));
        if !ok {
            return Err(format!(
                "SHA256SUMS: {:?} is not a lowercase 64-digit hex hash",
                hash
            ));
        }
        listed.push(String::from(name));
    }
    listed.sort();
    expect_eq("SHA256SUMS names vs GOLDEN_FILES", &listed, &table)
}
