//! Harper S2 goldens: every row of tests/goldens/slices/bridge/bridge_cmp.txt (scripts/gen_bridge_cmp_goldens.py, an
//! exact-integer oracle cross-checked with mpmath) is asserted exactly (CONTRACT 16.S2.1, 16.S2.2, 16.S2.4). N rows:
//! `Mag::cmp_exact_num`; D rows: `Mag::cmp_exact_diff`, including the few 100000-bit rows and the radius-above-EXPONENT_MAX
//! rows of harper_s2_diff_edge.rs (Mag(TOP, E+1), Mag(TOP+1, E+1), Mag(u64::MAX, E+1)); U rows: `Mag::from_exact_num` and
//! `Mag::from_exact_num_lower`; T rows: `Mag::to_exact_num`, values and error details. Each kind's row count is pinned.
//! COVERS: Mag::cmp_exact_num, Mag::cmp_exact_diff, Mag::from_exact_num, Mag::from_exact_num_lower, Mag::to_exact_num

use core::cmp::Ordering;

use zenith_float_bound::error::BoundError;
use zenith_float_bound::mag::Mag;
use zenith_float_num::{ExactNum, RawParts, RoundingMode, EXPONENT_MIN};

use super::{golden_value, raw};
use crate::common::goldens::golden_dir;

const HEADER: &str = "# zenith-float-bound S2 cmp/conversion rows, scripts/gen_bridge_cmp_goldens.py v2 (see its docstring)";

/// The bridge_cmp.txt rows of one kind; see [`rows_in`].
fn rows(kind: &str, want: usize) -> Result<Vec<(usize, String)>, String> {
    rows_in("bridge_cmp.txt", HEADER, kind, want)
}

/// The rows of one kind (`N`, `D`, `U` or `T`) of slices/bridge/`file` with their 1-based line numbers; fails unless the
/// first line is `header` and there are exactly `want` rows.
pub(super) fn rows_in(
    file: &str,
    header: &str,
    kind: &str,
    want: usize,
) -> Result<Vec<(usize, String)>, String> {
    let path = golden_dir().join("slices").join("bridge").join(file);
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{}: {}", path.display(), e))?;
    if text.lines().next() != Some(header) {
        return Err(format!("{}: header is not {:?}", path.display(), header));
    }
    let mut out = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if line.starts_with('#') {
            continue;
        }
        let (k, rest) = line
            .split_once(' ')
            .ok_or_else(|| format!("bad row {}", line))?;
        if k == kind {
            out.push((i.checked_add(1).ok_or("line")?, rest.to_string()));
        }
    }
    if out.len() != want {
        return Err(format!("{} {} rows, want {}", out.len(), kind, want));
    }
    Ok(out)
}

/// `x` with its sign flipped through raw parts (no arithmetic, so a subnormal is not flushed).
fn flip(x: &ExactNum) -> Result<ExactNum, String> {
    let p = x.as_raw_parts();
    ExactNum::from_raw_parts(&RawParts {
        negative: !p.negative,
        ..p
    })
    .map_err(|e| format!("flip: {:?}", e))
}

/// An ExactNum token: `nan`, `+inf`, `-inf`, `+0`, `-0` or `[+-]HEX^E` = +-HEX * 2^E.
pub(super) fn exact(tok: &str) -> Result<ExactNum, String> {
    match tok {
        "nan" => Ok(ExactNum::from(0i64).div(&ExactNum::from(0i64), 64, RoundingMode::ToEven)),
        "+inf" => Ok(ExactNum::from(1i64).div(&ExactNum::from(0i64), 64, RoundingMode::ToEven)),
        "-inf" => flip(&exact("+inf")?),
        "+0" => raw(&[0], 0, false),
        "-0" => raw(&[0], 0, true),
        other => {
            let (negative, body) = match (other.strip_prefix('-'), other.strip_prefix('+')) {
                (Some(b), None) => (true, b),
                (None, Some(b)) => (false, b),
                (Some(_), Some(_)) | (None, None) => return Err(format!("bad ExactNum {}", other)),
            };
            let (mant, exp) = body
                .split_once('^')
                .ok_or_else(|| format!("bad ExactNum {}", other))?;
            let v = golden_value(exp, mant)?;
            if negative {
                flip(&v)
            } else {
                Ok(v)
            }
        }
    }
}

/// A T-row value token `+HEX^E` (HEX at most 64 bits). Below zf's range it is built the way zb's `Mag::to_exact_num`
/// and zf store a subnormal: one word, mantissa shifted right (leading zeros), exponent `EXPONENT_MIN`; else as [`exact`].
fn t_value(tok: &str) -> Result<ExactNum, String> {
    let bad = || format!("bad T value {}", tok);
    // Zero and multi-word tokens are never subnormal-shaped here; `exact` handles them.
    let Some((mant, exp)) = tok.strip_prefix('+').and_then(|b| b.split_once('^')) else {
        return exact(tok);
    };
    let Ok(m) = u64::from_str_radix(mant, 16) else {
        return exact(tok);
    };
    let e: i64 = exp.parse().map_err(|_| bad())?;
    let lz = m.leading_zeros();
    let bits = 64u32.checked_sub(lz).ok_or_else(bad)?;
    let top = e.checked_add(i64::from(bits)).ok_or_else(bad)?;
    match EXPONENT_MIN.checked_sub(top).filter(|d| *d > 0) {
        None => exact(tok),
        Some(d) => {
            let d = u32::try_from(d).ok().filter(|d| *d < 64).ok_or_else(bad)?;
            let norm = m.checked_shl(lz).ok_or_else(bad)?;
            if norm.trailing_zeros() < d {
                return Err(format!("{}: not an exact subnormal", tok));
            }
            raw(&[norm.checked_shr(d).ok_or_else(bad)?], EXPONENT_MIN, false)
        }
    }
}

/// A Mag token: `0`, `inf` or `HEX16:E` = 0.HEX * 2^E.
pub(super) fn mag(tok: &str) -> Result<Mag, String> {
    match tok {
        "inf" => Ok(Mag::infinity()),
        "0" => Ok(Mag::zero()),
        other => {
            let (m, e) = other
                .split_once(':')
                .ok_or_else(|| format!("bad Mag {}", other))?;
            let mant = u64::from_str_radix(m, 16).map_err(|_| format!("bad Mag {}", other))?;
            let exp: i64 = e.parse().map_err(|_| format!("bad Mag {}", other))?;
            Mag::finite(mant, exp).map_err(|err| format!("Mag::finite {}: {:?}", other, err))
        }
    }
}

pub(super) fn ord(tok: &str) -> Result<Option<Ordering>, String> {
    match tok {
        "lt" => Ok(Some(Ordering::Less)),
        "eq" => Ok(Some(Ordering::Equal)),
        "gt" => Ok(Some(Ordering::Greater)),
        "none" => Ok(None),
        other => Err(format!("bad ORD {}", other)),
    }
}

/// Collects mismatches (every row is checked) and fails naming the first few.
pub(super) fn verdict(kind: &str, bad: &[String]) -> Result<(), String> {
    if bad.is_empty() {
        Ok(())
    } else {
        let shown: Vec<&String> = bad.iter().take(10).collect();
        Err(format!(
            "{} {} row(s) disagree with the oracle: {:?}",
            bad.len(),
            kind,
            shown
        ))
    }
}

#[test]
fn s2_goldens_cmp_exact_num_rows() -> Result<(), String> {
    let mut bad = Vec::new();
    for (line, row) in rows("N", 795)? {
        let f: Vec<&str> = row.split(' ').collect();
        let [r, x, o] = f.as_slice() else {
            return Err(format!("line {}: bad N row", line));
        };
        let got = mag(r)?.cmp_exact_num(&exact(x)?);
        if got != ord(o)? {
            bad.push(format!("line {}: N {} got {:?}", line, row, got));
        }
    }
    verdict("N", &bad)
}

// Expected-red exception set (CONTRACT 16.S2.2, FIXLIST ZB-34; the bug is FIXLIST ZB-35, harper_s2_diff_overflow.rs): these 8
// opposite-sign near-EMAX D rows return None instead of Less because `a.sub_full_prec(b)` overflows (src/bridge.rs:113/117).
// Each is matched by its row content and must still return None at its line, so the test goes red when the bug is fixed and
// the exception is not removed. Every other D row must match the oracle.
const D_EXPECTED_NONE: [(usize, &str); 8] = [
    (540, "0 -ffffffffffffffff^4611686018427387839 +ffffffffffffffff^4611686018427387839 lt"),
    (541, "8000000000000000:1 -ffffffffffffffff^4611686018427387839 +ffffffffffffffff^4611686018427387839 lt"),
    (546, "0 -ffffffffffffffff^4611686018427387839 +1ffffffffffffffff^4611686018427387838 lt"),
    (547, "8000000000000000:1 -ffffffffffffffff^4611686018427387839 +1ffffffffffffffff^4611686018427387838 lt"),
    (594, "0 -1ffffffffffffffff^4611686018427387838 +ffffffffffffffff^4611686018427387839 lt"),
    (595, "8000000000000000:1 -1ffffffffffffffff^4611686018427387838 +ffffffffffffffff^4611686018427387839 lt"),
    (600, "0 -1ffffffffffffffff^4611686018427387838 +1ffffffffffffffff^4611686018427387838 lt"),
    (601, "8000000000000000:1 -1ffffffffffffffff^4611686018427387838 +1ffffffffffffffff^4611686018427387838 lt"),
];

#[test]
fn s2_goldens_cmp_exact_diff_rows() -> Result<(), String> {
    let mut bad = Vec::new();
    let mut hit = Vec::new();
    let mut passed: usize = 0;
    for (line, row) in rows("D", 1028)? {
        let f: Vec<&str> = row.split(' ').collect();
        let [r, a, b, o] = f.as_slice() else {
            return Err(format!("line {}: bad D row", line));
        };
        let got = mag(r)?.cmp_exact_diff(&exact(a)?, &exact(b)?);
        match D_EXPECTED_NONE
            .iter()
            .find(|(_, text)| *text == row.as_str())
        {
            Some((at, _)) => {
                hit.push(line);
                if *at != line || got.is_some() {
                    bad.push(format!("line {}: D {} is a CONTRACT 16.S2.2 / FIXLIST ZB-34 exception (listed at line {}) but got {:?}, want None: remove the exception", line, row, at, got));
                }
            }
            None => {
                if got == ord(o)? {
                    passed = passed.checked_add(1).ok_or("D pass count overflow")?;
                } else {
                    bad.push(format!("line {}: D {} got {:?}", line, row, got));
                }
            }
        }
    }
    let listed: Vec<usize> = D_EXPECTED_NONE.iter().map(|(at, _)| *at).collect();
    if hit != listed {
        bad.push(format!(
            "D exceptions seen at lines {:?}, listed {:?}",
            hit, listed
        ));
    }
    if passed != 1020 {
        bad.push(format!(
            "D rows matching the oracle: {}, want 1020 (1028 - 8 exceptions)",
            passed
        ));
    }
    verdict("D", &bad)
}

#[test]
fn s2_goldens_from_exact_num_rows() -> Result<(), String> {
    let mut bad = Vec::new();
    for (line, row) in rows("U", 615)? {
        let f: Vec<&str> = row.split(' ').collect();
        let [x, up, lo] = f.as_slice() else {
            return Err(format!("line {}: bad U row", line));
        };
        let v = exact(x)?;
        let (gu, gl) = (Mag::from_exact_num(&v), Mag::from_exact_num_lower(&v));
        if gu != mag(up)? || gl != mag(lo)? {
            bad.push(format!("line {}: U {} got {:?} {:?}", line, row, gu, gl));
        }
    }
    verdict("U", &bad)
}

#[test]
fn s2_goldens_to_exact_num_rows() -> Result<(), String> {
    let mut bad = Vec::new();
    for (line, row) in rows("T", 25)? {
        let f: Vec<&str> = row.split(' ').collect();
        let [r, res] = f.as_slice() else {
            return Err(format!("line {}: bad T row", line));
        };
        let want_err = match res.strip_prefix("err:") {
            Some("inf") => Some("Infinity has no ExactNum value"),
            Some("max") => Some("exponent above zf EXPONENT_MAX"),
            Some("min") => Some("below zf EXPONENT_MIN and not an exact subnormal"),
            Some(other) => return Err(format!("line {}: bad error token {}", line, other)),
            None => None,
        };
        let ok = match (mag(r)?.to_exact_num(), want_err) {
            (Err(BoundError::Representation { op, detail }), Some(w)) => {
                op == "Mag::to_exact_num" && detail == w
            }
            (Ok(got), None) => {
                let want = t_value(res)?;
                got.as_raw_parts().negative == want.as_raw_parts().negative
                    && got.abs_cmp(&want) == Some(0)
            }
            (Err(e), want) => {
                bad.push(format!(
                    "line {}: T {} got {:?}, want {:?}",
                    line, row, e, want
                ));
                continue;
            }
            (Ok(got), Some(w)) => {
                bad.push(format!(
                    "line {}: T {} got {:?}, want {}",
                    line,
                    row,
                    got.as_raw_parts(),
                    w
                ));
                continue;
            }
        };
        if !ok {
            bad.push(format!("line {}: T {}", line, row));
        }
    }
    verdict("T", &bad)
}
