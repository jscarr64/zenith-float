//! CRLF guards (ZB-TD41). No new function is added here, so there is no
//! A-Z edge-case list: these tests guard the line-ending behaviour of the
//! parsers and golden files that already exist.
//!
//! The guards feed CRLF text through the REAL loaders, never through a
//! `lines()` / `split_whitespace()` call in the test itself (Spence review
//! of PR #13: the first version of the first two guards pre-stripped `\r`
//! in the test, so they only tested `std`):
//! - Ball rows go through `ball_goldens_selftest::rows()` and then
//!   `ballrows::parse_row`. The CR-safety lives in `rows()`' `lines()`:
//!   a bare `\r`-terminated row is rejected by `parse_row` (its last
//!   field, the inexact flag, would be `0\r`/`1\r`, not `0`/`1`), so if
//!   `rows()` ever switched to `split('\n')` this guard would fail.
//! - Replay lines go through `ReplayLine::parse_line`, whose
//!   `split_whitespace()` treats `\r` as whitespace; if it ever switched
//!   to a single-space split, the trailing `\r` would stick to the `b=`
//!   value and the round-trip equality below would fail.
//!
//! The on-disk guard (third test) is unchanged. The genuinely failing red
//! check for this PR remains `tools/gitattributes_selftest.sh`, which fails
//! while the repo has no `.gitattributes`. One honest limit, stated plainly:
//! `common::ballrows::parse_row` itself does not trim a raw line — every
//! caller feeds it `lines()` output, so no change is made to it here.

use crate::ball_goldens_selftest::rows;
use crate::common::ballrows::{format_row, parse_row, GoldRow};
use crate::common::expect_eq;
use crate::common::goldens::{golden_dir, GOLDEN_FILES};
use crate::common::replay::ReplayLine;

#[test]
fn crlf_ball_rows_parse_identically_through_real_loader() -> Result<(), String> {
    let path = golden_dir().join("ball_exact_ops.txt");
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{:?}: {}", path, e))?;
    let crlf_text = text.replace('\n', "\r\n");
    let lf_rows = rows(&text);
    let crlf_rows = rows(&crlf_text);
    expect_eq(
        "CRLF loader row count equals LF loader row count",
        &crlf_rows.len(),
        &lf_rows.len(),
    )?;
    let back = |row: &GoldRow| match format_row(row) {
        Ok(s) => s,
        Err(_) => String::from("wide"),
    };
    for (crlf_line, lf_line) in crlf_rows.iter().zip(lf_rows.iter()) {
        if crlf_line.contains('\r') {
            return Err(format!(
                "real loader yielded a CR-terminated row: {:?}",
                crlf_line
            ));
        }
        let crlf_row = parse_row(crlf_line)
            .map_err(|e| format!("CRLF row rejected: {:?}: {}", crlf_line, e))?;
        let lf_row =
            parse_row(lf_line).map_err(|e| format!("LF row rejected: {:?}: {}", lf_line, e))?;
        expect_eq(
            "CRLF row formats identically to LF row",
            &back(&crlf_row),
            &back(&lf_row),
        )?;
    }
    Ok(())
}

#[test]
fn crlf_replay_line_parses_identically_through_real_loader() -> Result<(), String> {
    let line = ReplayLine {
        suite: String::from("mag_add"),
        seed: 0x0123_4567_89AB_CDEF,
        case: 42,
        op: String::from("add"),
        prec: 53,
        rm: String::from("ToEven"),
        a: String::from("+8000000000000000e-5~r=ulp"),
        b: String::from("-0~r=zero"),
    };
    let lf_text = line.format_line()?;
    let crlf_text = format!("{}\r", lf_text);
    let lf_parsed = ReplayLine::parse_line(&lf_text)?;
    let crlf_parsed = ReplayLine::parse_line(&crlf_text)?;
    expect_eq(
        "CRLF replay line parses identically to LF replay line",
        &crlf_parsed,
        &lf_parsed,
    )?;
    expect_eq(
        "CRLF replay line round-trips to the original line",
        &crlf_parsed,
        &line,
    )
}

fn check_no_cr_in_dir(dir: &std::path::Path) -> Result<(), String> {
    let entries = std::fs::read_dir(dir).map_err(|e| format!("{:?}: {}", dir, e))?;
    for entry in entries {
        let entry = entry.map_err(|e| format!("dir entry: {}", e))?;
        let path = entry.path();
        let kind = entry
            .file_type()
            .map_err(|e| format!("{:?}: {}", path, e))?;
        if kind.is_dir() {
            check_no_cr_in_dir(&path)?;
        } else if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
            if name.ends_with(".txt") || name == "SHA256SUMS" {
                let text =
                    std::fs::read_to_string(&path).map_err(|e| format!("{:?}: {}", path, e))?;
                if text.contains('\r') {
                    return Err(format!("{:?} contains CR", path));
                }
            }
        }
    }
    Ok(())
}

#[test]
fn golden_files_and_slice_sums_contain_no_cr() -> Result<(), String> {
    for g in GOLDEN_FILES.iter() {
        let path = golden_dir().join(g.name);
        let text = std::fs::read_to_string(&path).map_err(|e| format!("{:?}: {}", path, e))?;
        if text.contains('\r') {
            return Err(format!("{:?} contains CR", path));
        }
    }
    check_no_cr_in_dir(&golden_dir())
}
