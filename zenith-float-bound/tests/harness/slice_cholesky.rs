//! Slice LA1 harness rows (`src/cholesky.rs`). Owned by slice LA1 alone; gate-row shape tests green; the replay test is ignored (red) until the API lands.
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

/// Token count of a logdet row of order n: case_id n A[n*n], then sign det lo256 hi256 or ERR Singular singular.
fn logdet_row_width(dims: &str, refusal: bool) -> Result<usize, String> {
    let n = dims
        .parse::<usize>()
        .map_err(|e| format!("dims {dims:?}: {e}"))?;
    let tail = if refusal { 3_usize } else { 4_usize };
    n.checked_mul(n)
        .and_then(|nn| nn.checked_add(2))
        .and_then(|v| v.checked_add(tail))
        .ok_or_else(|| "row width overflows".to_string())
}

#[test]
fn la1_logdet_gate_rows_match_header_count_and_field_widths() -> Result<(), String> {
    let cosmo = check_rows(&logdet_gate_text()?, "logdet.cosmo_", logdet_row_width)?;
    if cosmo == 0 {
        return Err("no logdet.cosmo_* rows in the gate file".to_string());
    }
    Ok(())
}

#[test]
fn la1_logdet_row_width_rejects_a_long_row() -> Result<(), String> {
    let text = "# rows: 1\nlogdet.t 1 2 1 2 0 1 9\n";
    match check_rows(text, "logdet.cosmo_", logdet_row_width) {
        Err(e) if e.contains("want 7") => Ok(()),
        Err(e) => Err(format!("wrong error: {e}")),
        Ok(n) => Err(format!("long row accepted ({n} cosmo)")),
    }
}

/// Red until slice LA1 lands zb logdet (LA1 mutants 1000..=1099; the LA4 cosmo rows are tracked under 1300..=1399).
#[test]
#[ignore = "red: not implemented: zb logdet (CONTRACT LA4.2: sign/det/ln enclosure on the cosmo_* la_logdet gate rows; API from slice LA1); un-ignore when src/cholesky.rs exports it"]
fn la1_logdet_gate_rows_replay_through_zb() -> Result<(), String> {
    let rows = header_rows(&logdet_gate_text()?)?;
    Err(format!("not implemented: zb logdet API missing (src/cholesky.rs is an S0 stub); {rows} gate rows waiting"))
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
