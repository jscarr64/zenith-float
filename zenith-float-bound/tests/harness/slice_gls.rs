//! Slice LA4 harness rows (`src/gls.rs`). Owned by slice LA4 alone; gate-row shape tests green; the replay test is ignored (red) until the API lands.
//! Exact golds: Cassie (owns the zb LA exact golds, 2026-10-24).
//! Goldens: tests/goldens/slices/gls/ (its own SHA256SUMS and GENERATORS). Mutants: tools/mutants/manifest.d/gls.txt, ids 1300..=1399.
//! Golden file headers for this slice are declared here, not in tests/common/goldens.rs.

use crate::common::goldens::golden_dir;

/// The LA4 gate file, read at run time from the checkout under test.
fn gls_gate_text() -> Result<String, String> {
    let path = golden_dir().join("slices").join("gls").join("la_gls.txt");
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
    width: fn(&str, bool) -> Result<usize, String>,
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

/// Token count of a GLS row with design m x p: case_id dims X[m*p] S[m*m] y[m], then beta[p] cov[p*p] rss or ERR variant detail.
fn gls_row_width(dims: &str, refusal: bool) -> Result<usize, String> {
    let (m_s, p_s) = dims
        .split_once('x')
        .ok_or_else(|| format!("dims {dims:?} is not mxp"))?;
    let m = m_s
        .parse::<usize>()
        .map_err(|e| format!("dims {dims:?}: {e}"))?;
    let p = p_s
        .parse::<usize>()
        .map_err(|e| format!("dims {dims:?}: {e}"))?;
    let mp = m.checked_mul(p).ok_or("m*p overflows")?;
    let mm = m.checked_mul(m).ok_or("m*m overflows")?;
    let pp = p.checked_mul(p).ok_or("p*p overflows")?;
    let tail = if refusal {
        Some(3_usize)
    } else {
        p.checked_add(pp).and_then(|v| v.checked_add(1))
    };
    let parts = [mp, mm, m, tail.ok_or("tail overflows")?];
    parts
        .iter()
        .try_fold(2_usize, |acc, v| acc.checked_add(*v))
        .ok_or_else(|| "row width overflows".to_string())
}

#[test]
fn la4_gls_gate_rows_match_header_count_and_field_widths() -> Result<(), String> {
    let cosmo = check_rows(&gls_gate_text()?, "gls.cosmo_", gls_row_width)?;
    if cosmo == 0 {
        return Err("no gls.cosmo_* rows in the gate file".to_string());
    }
    Ok(())
}

#[test]
fn la4_gls_row_width_rejects_a_short_row() -> Result<(), String> {
    let text = "# rows: 1\ngls.t 1x1 1 1 1 1 1\n";
    match check_rows(text, "gls.cosmo_", gls_row_width) {
        Err(e) if e.contains("want 8") => Ok(()),
        Err(e) => Err(format!("wrong error: {e}")),
        Ok(n) => Err(format!("short row accepted ({n} cosmo)")),
    }
}

/// Red until slice LA4 lands zb::gls (mutant block 1300..=1399): replays every gate row through the implementation.
#[test]
#[ignore = "red: not implemented: zb::gls (CONTRACT LA4.1: GLS beta/cov/rss and refusals on every la_gls gate row); un-ignore when src/gls.rs exports it"]
fn la4_gls_gate_rows_replay_through_zb() -> Result<(), String> {
    let rows = header_rows(&gls_gate_text()?)?;
    Err(format!(
        "not implemented: zb::gls API missing (src/gls.rs is an S0 stub); {rows} gate rows waiting"
    ))
}

/// The LA slices whose gate files share the 1 MB budget (Spence 2026-10-24); heavy files live outside the repo.
const LA_SLICES: [&str; 5] = ["cholesky", "eigen", "gls", "lu", "qr"];
const LA_GATE_BUDGET_BYTES: u64 = 1_000_000;

#[test]
fn la_gate_set_totals_under_1_mb() -> Result<(), String> {
    let mut total = 0_u64;
    for slice in LA_SLICES {
        let dir = golden_dir().join("slices").join(slice);
        for entry in
            std::fs::read_dir(&dir).map_err(|e| format!("cannot list {}: {e}", dir.display()))?
        {
            let entry = entry.map_err(|e| format!("{}: {e}", dir.display()))?;
            let name = entry.file_name();
            let name = name.to_string_lossy();
            if name.starts_with("la_") && name.ends_with(".txt") && !name.ends_with("_heavy.txt") {
                let len = entry.metadata().map_err(|e| format!("{name}: {e}"))?.len();
                total = total.checked_add(len).ok_or("size overflow")?;
            }
        }
    }
    if total >= LA_GATE_BUDGET_BYTES {
        return Err(format!("LA gate set is {total} bytes, budget {LA_GATE_BUDGET_BYTES}: move rows to the heavy tier"));
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
fn la4_gls_every_cosmo_family_has_a_gate_row() -> Result<(), String> {
    let fams = [
        "bao",
        "cmb_toeplitz",
        "blockdiag",
        "lowrank",
        "cond",
        "spread",
        "refuse",
    ];
    families_present(&gls_gate_text()?, "gls.cosmo_", &fams)
}
