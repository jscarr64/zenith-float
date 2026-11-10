//! Slice LA1 harness rows (`src/cholesky.rs`). Owned by slice LA1 alone; gate-row shape tests green; the per-prec replay tests (16.LA1.7: every row at 53, 113, 192) are ignored (red) until LA2 lands BallMat::logdet.
//! Exact golds: Cassie (owns the zb LA exact golds, 2026-10-24).
//! Goldens: tests/goldens/slices/cholesky/ (its own SHA256SUMS and GENERATORS). Mutants: tools/mutants/manifest.d/cholesky.txt, ids 1000..=1099.
//! Golden file headers for this slice are declared here, not in tests/common/goldens.rs.

use crate::common::goldens::golden_dir;

/// The LA1 logdet gate file, read at run time from the checkout under test.
fn logdet_gate_text() -> Result<String, String> {
    let path = golden_dir()
        .join("slices")
        .join("cholesky")
        .join("la_logdet.txt");
    std::fs::read_to_string(&path).map_err(|e| format!("cannot read {}: {e}", path.display()))
}

/// The `# rows: N` header value.
fn header_rows(text: &str) -> Result<usize, String> {
    let line = text
        .lines()
        .find_map(|l| l.strip_prefix("# rows: "))
        .ok_or("no '# rows:' header")?;
    line.trim()
        .parse::<usize>()
        .map_err(|e| format!("bad '# rows:' header {line:?}: {e}"))
}

/// Checks every data row against `width(dims, refusal)` and the `# rows:` header; returns the number of rows whose id starts with `cosmo`.
fn check_rows(
    text: &str,
    cosmo_prefix: &str,
    width: &dyn Fn(&str, bool) -> Result<usize, String>,
) -> Result<usize, String> {
    let want = header_rows(text)?;
    let mut seen = 0_usize;
    let mut cosmo = 0_usize;
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let toks: Vec<&str> = line.split_whitespace().collect();
        let mut it = toks.iter();
        let id = it.next().ok_or("empty row")?;
        let dims = it.next().ok_or_else(|| format!("{id}: no dims"))?;
        let refusal = toks.contains(&"ERR");
        let w = width(dims, refusal)?;
        if toks.len() != w {
            return Err(format!(
                "{id}: {} fields, want {w} for dims {dims} (refusal {refusal})",
                toks.len()
            ));
        }
        if id.starts_with(cosmo_prefix) {
            cosmo = cosmo.checked_add(1).ok_or("count overflow")?;
        }
        seen = seen.checked_add(1).ok_or("count overflow")?;
    }
    if seen != want {
        return Err(format!("{seen} data rows, header says {want}"));
    }
    Ok(cosmo)
}

/// Gold format of a logdet file: `# format: N` (2 adds the six cond columns of scripts/la_cond_cols.py); no line means format 1.
fn gold_format(text: &str) -> Result<u32, String> {
    match text.lines().find_map(|l| l.strip_prefix("# format: ")) {
        None => Ok(1),
        Some(rest) => rest
            .split_whitespace()
            .next()
            .ok_or("empty '# format:' header")?
            .parse::<u32>()
            .map_err(|e| format!("bad '# format:' header {rest:?}: {e}")),
    }
}

/// Cond columns per non-refusal row: none in format 1, six in format 2 (cond_row k_row scale_row containment_only reason rref).
fn cond_cols_for(format: u32) -> Result<usize, String> {
    match format {
        1 => Ok(0),
        2 => Ok(6),
        other => Err(format!("unknown gold format {other}")),
    }
}

/// Token count of a logdet row of order n: case_id n A[n*n], then sign det lo256 hi256 and `cond` cond columns, or ERR Singular singular.
fn logdet_row_width(dims: &str, refusal: bool, cond: usize) -> Result<usize, String> {
    let n = dims
        .parse::<usize>()
        .map_err(|e| format!("dims {dims:?}: {e}"))?;
    let tail = if refusal {
        3_usize
    } else {
        4_usize.checked_add(cond).ok_or("row width overflows")?
    };
    n.checked_mul(n)
        .and_then(|nn| nn.checked_add(2))
        .and_then(|v| v.checked_add(tail))
        .ok_or_else(|| "row width overflows".to_string())
}

/// The replay precisions of CONTRACT 16.LA1.7 (the file has no prec column).
const REPLAY_PRECS: [u32; 3] = [53, 113, 192];

/// Largest replay precision; `containment_only` is 1 (k_row_too_large) exactly when `k_row + 8 >= MAX_REPLAY_PREC`; 2 is lstsq_x_zero (any k_row) (scripts/la_cond_cols.py).
const MAX_REPLAY_PREC: u32 = 192;

/// The format-2 cond columns the replay uses (cond_row, scale_row and rref are exact p/q, checked for shape only).
struct CondCols {
    k_row: u32,
    containment_only: bool,
}

/// The closed containment_only reason set (scripts/la_cond_cols.py REASONS): tight (flag 0), k_row_too_large (flag 1,
/// k_row + 8 >= MAX_REPLAY_PREC), lstsq_x_zero (flag 2, LA3 lstsq rows with x = 0: scale_row undefined, any k_row).
const COND_REASONS: [&str; 3] = ["tight", "k_row_too_large", "lstsq_x_zero"];

/// Parses and cross-checks the six cond columns (in order cond_row k_row scale_row containment_only reason rref).
fn parse_cond_cols(id: &str, c: &[&str]) -> Result<CondCols, String> {
    let [cond_row, k_row, scale_row, flag, reason, rref] = c else {
        return Err(format!("{id}: want 6 cond columns, got {}", c.len()));
    };
    for (name, v) in [
        ("cond_row", cond_row),
        ("scale_row", scale_row),
        ("rref", rref),
    ] {
        let (p, q) = v
            .split_once('/')
            .ok_or_else(|| format!("{id}: {name} {v:?} is not p/q"))?;
        let digits = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
        if !(digits(p.strip_prefix('-').unwrap_or(p)) && digits(q)) {
            return Err(format!("{id}: {name} {v:?} is not an exact rational p/q"));
        }
    }
    let k_row = k_row
        .parse::<u32>()
        .map_err(|e| format!("{id}: k_row {k_row:?}: {e}"))?;
    let too_large = k_row.checked_add(8).ok_or("k_row overflows")? >= MAX_REPLAY_PREC;
    if !COND_REASONS.contains(reason) {
        return Err(format!(
            "{id}: unknown containment_only reason {reason:?} (closed set {COND_REASONS:?})"
        ));
    }
    // (flag, reason) -> (containment_only, the k_row + 8 >= MAX_REPLAY_PREC value it implies; None: independent of k_row)
    let (containment_only, implied) = match (*flag, *reason) {
        ("0", "tight") => (false, Some(false)),
        ("1", "k_row_too_large") => (true, Some(true)),
        ("2", "lstsq_x_zero") => (true, None),
        (f, r) => return Err(format!("{id}: bad containment_only/reason {f} {r}")),
    };
    if implied.is_some_and(|t| t != too_large) {
        return Err(format!(
            "{id}: containment_only {flag} disagrees with k_row {k_row} (k_row + 8 >= {MAX_REPLAY_PREC} is {too_large})"
        ));
    }
    Ok(CondCols {
        k_row,
        containment_only,
    })
}

/// Replay tightness of CONTRACT 16.LA.5: `rad <= 2^(k_row - prec) * max(|v|, scale_row)` when `prec > k_row + 8`.
#[derive(Debug, PartialEq, Eq)]
enum Tightness {
    /// Enforced bound; the factor is `2^log2_factor` with `log2_factor = k_row - prec`.
    Bound { log2_factor: i64 },
    /// Containment only (prec too low for the row, a format-1 row, or a `containment_only` row).
    ContainmentOnly,
}

fn tightness(prec: u32, cond: Option<&CondCols>) -> Result<Tightness, String> {
    let Some(c) = cond else {
        return Ok(Tightness::ContainmentOnly);
    };
    let floor = c.k_row.checked_add(8).ok_or("k_row overflows")?;
    if c.containment_only || prec <= floor {
        return Ok(Tightness::ContainmentOnly);
    }
    let log2_factor = i64::from(c.k_row)
        .checked_sub(i64::from(prec))
        .ok_or("exponent overflows")?;
    Ok(Tightness::Bound { log2_factor })
}

/// One replay row: its id and, for a non-refusal format-2 row, its cond columns.
struct LogdetRow {
    id: String,
    refusal: bool,
    cond: Option<CondCols>,
}

fn logdet_rows(text: &str) -> Result<Vec<LogdetRow>, String> {
    let ncond = cond_cols_for(gold_format(text)?)?;
    check_rows(text, "logdet.cosmo_", &|d, r| logdet_row_width(d, r, ncond))?;
    let mut rows = Vec::new();
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let toks: Vec<&str> = line.split_whitespace().collect();
        let id = toks.first().ok_or("empty row")?;
        let refusal = toks.contains(&"ERR");
        let cond = if refusal || ncond == 0 {
            None
        } else {
            let skip = toks.len().checked_sub(ncond).ok_or("short row")?;
            Some(parse_cond_cols(id, toks.get(skip..).ok_or("short row")?)?)
        };
        rows.push(LogdetRow {
            id: (*id).to_string(),
            refusal,
            cond,
        });
    }
    Ok(rows)
}

/// The red replay body for one prec: plans every row's rule (refusal, enforced bound, containment) and stops at the missing API.
fn replay_red(prec: u32, prefix: &str) -> Result<(), String> {
    let all = logdet_rows(&logdet_gate_text()?)?;
    let rows: Vec<&LogdetRow> = all.iter().filter(|r| r.id.starts_with(prefix)).collect();
    if rows.is_empty() {
        return Err(format!("no la_logdet rows match {prefix:?}"));
    }
    let (mut refusals, mut bound, mut contain) = (0_usize, 0_usize, 0_usize);
    let mut widest = None;
    for r in rows.iter() {
        if r.refusal {
            refusals = refusals.checked_add(1).ok_or("count overflow")?;
            continue;
        }
        match tightness(prec, r.cond.as_ref())? {
            Tightness::Bound { log2_factor } => {
                bound = bound.checked_add(1).ok_or("count overflow")?;
                widest = Some((r.id.as_str(), log2_factor));
            }
            Tightness::ContainmentOnly => {
                contain = contain.checked_add(1).ok_or("count overflow")?
            }
        }
    }
    Err(format!(
        "not implemented: BallMat::logdet (CONTRACT 16.LA1.4, LA2 src/lu.rs) is not on main; prec {prec}: {} {prefix}* rows planned \
         ({refusals} refusal, {bound} k_row bound, {contain} containment only; last bound row {widest:?})",
        rows.len()
    ))
}

#[test]
fn la1_logdet_gate_rows_match_header_count_and_field_widths() -> Result<(), String> {
    let text = logdet_gate_text()?;
    let cond = cond_cols_for(gold_format(&text)?)?;
    let cosmo = check_rows(&text, "logdet.cosmo_", &|d, r| logdet_row_width(d, r, cond))?;
    if cosmo == 0 {
        return Err("no logdet.cosmo_* rows in the gate file".to_string());
    }
    Ok(())
}

#[test]
fn la1_logdet_row_width_rejects_a_long_row() -> Result<(), String> {
    let text = "# rows: 1\nlogdet.t 1 2 1 2 0 1 9\n";
    match check_rows(text, "logdet.cosmo_", &|d, r| logdet_row_width(d, r, 0)) {
        Err(e) if e.contains("want 7") => Ok(()),
        Err(e) => Err(format!("wrong error: {e}")),
        Ok(n) => Err(format!("long row accepted ({n} cosmo)")),
    }
}

/// A format-2 row (n = 1, cond 1, k_row = 4 + 0 + 0) parses; a format-1 row of the same file shape is rejected by width.
const V2_ROW: &str = "# format: 2 (x)\n# rows: 2\nlogdet.t 1 2/1 1 2/1 1*2^-1 3*2^-2 1/1 4 1/1 0 tight 1/1\nlogdet.z 1 0/1 ERR Singular singular\n";

#[test]
fn la1_logdet_format2_rows_parse_cond_columns() -> Result<(), String> {
    if gold_format(V2_ROW)? != 2 || gold_format("# rows: 0\n")? != 1 {
        return Err("format header not read".to_string());
    }
    let rows = logdet_rows(V2_ROW)?;
    let ks: Vec<Option<u32>> = rows
        .iter()
        .map(|r| r.cond.as_ref().map(|c| c.k_row))
        .collect();
    if ks != [Some(4), None] {
        return Err(format!("k_row per row {ks:?}, want [Some(4), None]"));
    }
    match logdet_rows(&V2_ROW.replace("# format: 2 (x)\n", "")) {
        Err(e) if e.contains("want 7") => Ok(()),
        Err(e) => Err(format!("wrong error: {e}")),
        Ok(_) => Err("format-2 row accepted as format 1".to_string()),
    }
}

#[test]
fn la1_logdet_format2_rejects_bad_cond_columns() -> Result<(), String> {
    let good = " 1/1 4 1/1 0 tight 1/1\n";
    for (bad, want) in [
        (" 1/1 190 1/1 0 tight 1/1\n", "disagrees"),
        (" 1/1 4 1/1 1 tight 1/1\n", "bad containment_only"),
        (" 1.0/1 4 1/1 0 tight 1/1\n", "not an exact rational"),
    ] {
        match logdet_rows(&V2_ROW.replace(good, bad)) {
            Err(e) if e.contains(want) => {}
            Err(e) => return Err(format!("{bad:?}: wrong error {e}")),
            Ok(_) => return Err(format!("{bad:?}: bad cond columns accepted")),
        }
    }
    Ok(())
}

/// The reason token is from the closed set: an unknown one is rejected, lstsq_x_zero (code 2) is accepted at any k_row.
#[test]
fn la1_logdet_format2_reason_set_is_closed() -> Result<(), String> {
    let good = " 1/1 4 1/1 0 tight 1/1\n";
    for (bad, want) in [
        (
            " 1/1 4 1/1 1 x_zero 1/1\n",
            "unknown containment_only reason",
        ),
        (
            " 1/1 4 1/1 0 Tight 1/1\n",
            "unknown containment_only reason",
        ),
        (" 1/1 4 1/1 0 lstsq_x_zero 1/1\n", "bad containment_only"),
        (" 1/1 4 1/1 1 lstsq_x_zero 1/1\n", "bad containment_only"),
    ] {
        match logdet_rows(&V2_ROW.replace(good, bad)) {
            Err(e) if e.contains(want) => {}
            Err(e) => return Err(format!("{bad:?}: wrong error {e}")),
            Ok(_) => return Err(format!("{bad:?}: bad reason accepted")),
        }
    }
    let rows = logdet_rows(&V2_ROW.replace(good, " 1/1 4 1/1 2 lstsq_x_zero 1/1\n"))?;
    let flags: Vec<Option<bool>> = rows
        .iter()
        .map(|r| r.cond.as_ref().map(|c| c.containment_only))
        .collect();
    if flags != [Some(true), None] {
        return Err(format!(
            "lstsq_x_zero row flags {flags:?}, want [Some(true), None]"
        ));
    }
    Ok(())
}

#[test]
fn la1_tightness_rule_is_enforced_only_above_k_row_plus_8() -> Result<(), String> {
    let at = |k_row: u32, containment_only: bool, prec: u32| {
        tightness(
            prec,
            Some(&CondCols {
                k_row,
                containment_only,
            }),
        )
    };
    let got = [
        at(44, false, 53)?,
        at(45, false, 53)?,
        at(4, false, 192)?,
        at(44, true, 192)?,
        tightness(192, None)?,
    ];
    let want = [
        Tightness::Bound { log2_factor: -9 },
        Tightness::ContainmentOnly,
        Tightness::Bound { log2_factor: -188 },
        Tightness::ContainmentOnly,
        Tightness::ContainmentOnly,
    ];
    if got != want {
        return Err(format!("tightness {got:?}, want {want:?}"));
    }
    if REPLAY_PRECS.iter().max() != Some(&MAX_REPLAY_PREC) {
        return Err("MAX_REPLAY_PREC is not the largest replay prec".to_string());
    }
    Ok(())
}

/// Red until LA2 lands `BallMat::logdet` (stub `Resource "not_implemented"` on spence/zb-la-api; absent on main). All 51 rows, prec 53.
#[test]
#[ignore = "red: CONTRACT 16.LA1.4: BallMat::logdet (LA2, src/lu.rs) not implemented; 16.LA1.7 replay of all la_logdet rows at prec 53"]
fn la1_logdet_rows_replay_through_zb_prec53() -> Result<(), String> {
    replay_red(53, "logdet.")
}

/// As above at prec 113.
#[test]
#[ignore = "red: CONTRACT 16.LA1.4: BallMat::logdet (LA2, src/lu.rs) not implemented; 16.LA1.7 replay of all la_logdet rows at prec 113"]
fn la1_logdet_rows_replay_through_zb_prec113() -> Result<(), String> {
    replay_red(113, "logdet.")
}

/// As above at prec 192.
#[test]
#[ignore = "red: CONTRACT 16.LA1.4: BallMat::logdet (LA2, src/lu.rs) not implemented; 16.LA1.7 replay of all la_logdet rows at prec 192"]
fn la1_logdet_rows_replay_through_zb_prec192() -> Result<(), String> {
    replay_red(192, "logdet.")
}

/// Red until slice LA1/LA2 lands zb logdet (the LA4 cosmo rows are tracked under mutants 1300..=1399): the `cosmo_*` rows at every
/// replay prec, kept as their own CONTRACT LA4.2 marker so tools/la4_red_markers.sh keeps LA4 from being marked done meanwhile.
#[test]
#[ignore = "red: not implemented: zb logdet (CONTRACT LA4.2: sign/det/ln enclosure on the cosmo_* la_logdet gate rows; API from slice LA1); un-ignore when src/cholesky.rs exports it"]
fn la4_logdet_cosmo_rows_replay_through_zb() -> Result<(), String> {
    for prec in REPLAY_PRECS {
        replay_red(prec, "logdet.cosmo_")?;
    }
    Ok(())
}

/// Every cosmo family keeps at least one gate row (Spence 2026-10-24) even when its big rows live in the heavy tier.
fn families_present(text: &str, prefix: &str, families: &[&str]) -> Result<(), String> {
    for fam in families {
        let want = format!("{prefix}{fam}.");
        if !text.lines().any(|l| l.starts_with(&want)) {
            return Err(format!("no gate row for cosmo family {want}*"));
        }
    }
    Ok(())
}

#[test]
fn la1_logdet_every_cosmo_family_has_a_gate_row() -> Result<(), String> {
    let fams = [
        "spread",
        "cond",
        "toeplitz",
        "blockdiag",
        "lowrank",
        "fisher",
        "refuse",
    ];
    families_present(&logdet_gate_text()?, "logdet.cosmo_", &fams)
}
