//! COVERS: common::TightnessLog::empty, common::TightnessLog::record, common::TightnessLog::render, common::TightnessLog::stats, common::TightnessLog::write_results, common::tightness::check_tight, common::tightness::p99_rank, common::tightness::parse_results, common::tightness::slack_q16
//! S12 self-tests of tests/common/tightness.rs: the exact q16 ratio, the tight check, the per-cell statistics, the results text that
//! tools/tightness_report.py reads, and the strict parser. The pinned text below is the same fixture as RESULTS_TEXT in
//! tools/tightness_selftest.sh, so the Rust writer and the Python reader agree on the format. Oracle: hand arithmetic. Every checker
//! gets a wrong input it must reject (a ratio one unit over the tight bound, malformed rows, an unsorted file).

use crate::common::expect_eq;
use crate::common::tightness::{
    check_tight, p99_rank, parse_results, slack_q16, CellStat, TightnessLog, Q16_ONE,
};

const PINNED: &str = "# zenith-bound tightness results (S12), written by the tightness sets; ratio = rad / R_ref in units of 1/65536 (q16), rounded up\n# columns: OP PREC COUNT MIN MEDIAN P99 MAX\nadd 24 3 65536 65536 98304 98304\nadd 53 2 0 0 131072 131072\nmul 24 4 10 20 40 40\n";

fn fixture_log() -> Result<TightnessLog, String> {
    let mut log = TightnessLog::empty();
    // recorded out of order on purpose: the text is sorted by (op, prec)
    log.record("mul", 24, 40, Q16_ONE.into())?;
    log.record("add", 53, 2, 1)?;
    log.record("add", 24, 3, 2)?;
    log.record("mul", 24, 10, Q16_ONE.into())?;
    log.record("add", 53, 0, 7)?;
    log.record("add", 24, 5, 5)?;
    log.record("mul", 24, 30, Q16_ONE.into())?;
    log.record("add", 24, 9, 9)?;
    log.record("mul", 24, 20, Q16_ONE.into())?;
    Ok(log)
}

#[test]
fn slack_q16_hand_values() -> Result<(), String> {
    expect_eq("rad == R_ref", &slack_q16(7, 7)?, &Q16_ONE)?;
    expect_eq("rad 0", &slack_q16(0, 5)?, &0)?;
    expect_eq("3/2", &slack_q16(3, 2)?, &98_304)?;
    // 65536 / 3 = 21845.33..., rounded UP
    expect_eq("1/3 rounds up", &slack_q16(1, 3)?, &21_846)?;
    expect_eq("exact quotient is not bumped", &slack_q16(1, 4)?, &16_384)?;
    expect_eq(
        "exact quotient 65535/65536",
        &slack_q16(65_535, 65_536)?,
        &65_535,
    )?;
    Ok(())
}

#[test]
fn slack_q16_rejects_a_zero_reference_and_overflow() -> Result<(), String> {
    for (label, r) in [
        ("R_ref 0", slack_q16(1, 0)),
        ("rad * 65536 overflows u128", slack_q16(u128::MAX, 1)),
        ("ratio overflows u64", slack_q16(1 << 100, 1)),
    ] {
        if r.is_ok() {
            return Err(format!("slack_q16 accepted: {}", label));
        }
    }
    Ok(())
}

#[test]
fn check_tight_is_exact_at_the_bound() -> Result<(), String> {
    check_tight(8, 2, 4)?;
    check_tight(0, 2, 4)?;
    if check_tight(9, 2, 4).is_ok() {
        return Err(String::from("one unit over tau * R_ref was accepted"));
    }
    if check_tight(1, u128::MAX, 2).is_ok() {
        return Err(String::from("an overflowing bound was accepted"));
    }
    Ok(())
}

#[test]
fn log_renders_the_pinned_text_and_stats_follow_the_definition() -> Result<(), String> {
    let log = fixture_log()?;
    expect_eq("rendered text", &log.render()?.as_str(), &PINNED)?;
    let stats = log.stats()?;
    let want = [
        CellStat {
            op: String::from("add"),
            prec: 24,
            count: 3,
            min: 65_536,
            median: 65_536,
            p99: 98_304,
            max: 98_304,
        },
        CellStat {
            op: String::from("add"),
            prec: 53,
            count: 2,
            min: 0,
            median: 0,
            p99: 131_072,
            max: 131_072,
        },
        CellStat {
            op: String::from("mul"),
            prec: 24,
            count: 4,
            min: 10,
            median: 20,
            p99: 40,
            max: 40,
        },
    ];
    expect_eq("cell count", &stats.len(), &want.len())?;
    for (got, w) in stats.iter().zip(want.iter()) {
        expect_eq(&format!("{} {}", w.op, w.prec), got, w)?;
    }
    // 100 samples 1..=100 (q16 = i): lower median is the 50th, p99 the 99th, max the 100th
    let mut big = TightnessLog::empty();
    for i in 1u128..=100 {
        big.record("sqrt", 64, i, Q16_ONE.into())?;
    }
    let s = big.stats()?;
    expect_eq("100-sample cell", &s.len(), &1)?;
    for c in s.iter() {
        expect_eq("min", &c.min, &1)?;
        expect_eq("median", &c.median, &50)?;
        expect_eq("p99", &c.p99, &99)?;
        expect_eq("max", &c.max, &100)?;
        expect_eq("count", &c.count, &100)?;
    }
    Ok(())
}

#[test]
fn parse_results_inverts_render_and_rejects_malformed_files() -> Result<(), String> {
    let parsed = parse_results(PINNED)?;
    expect_eq("rows", &parsed.len(), &3)?;
    expect_eq("render then parse", &parsed, &fixture_log()?.stats()?)?;
    let head = "# h\n";
    for (label, body) in [
        ("six fields", "add 24 3 1 1 1\n"),
        ("eight fields", "add 24 3 1 1 1 1 1\n"),
        ("MIN above MAX", "add 24 3 5 5 5 4\n"),
        ("MEDIAN above P99", "add 24 3 1 3 2 4\n"),
        ("count 0", "add 24 0 0 0 0 0\n"),
        ("prec 0", "add 0 1 0 0 0 0\n"),
        ("leading zero", "add 024 3 1 1 1 1\n"),
        ("not a number", "add 24 3 1 1 1 x\n"),
        ("plus sign", "add 24 3 +1 1 1 1\n"),
        ("duplicate", "add 24 3 1 1 1 1\nadd 24 3 1 1 1 1\n"),
        ("unsorted by op", "mul 24 3 1 1 1 1\nadd 24 3 1 1 1 1\n"),
        ("unsorted by prec", "add 53 3 1 1 1 1\nadd 24 3 1 1 1 1\n"),
        ("double space", "add  24 3 1 1 1 1\n"),
        ("CR", "add 24 3 1 1 1 1\r\n"),
        ("no final newline", "add 24 3 1 1 1 1"),
    ] {
        if parse_results(&format!("{}{}", head, body)).is_ok() {
            return Err(format!("parse_results accepted: {}", label));
        }
    }
    // an empty file (comments only) is the empty result set
    expect_eq("comments only", &parse_results(head)?.len(), &0)
}

#[test]
fn record_rejects_bad_names_precisions_and_references() -> Result<(), String> {
    let mut log = TightnessLog::empty();
    for (label, r) in [
        ("empty op", log.record("", 24, 1, 1)),
        ("upper case op", log.record("Add", 24, 1, 1)),
        ("space in op", log.record("a b", 24, 1, 1)),
        ("prec 0", log.record("add", 0, 1, 1)),
        ("R_ref 0", log.record("add", 24, 1, 0)),
    ] {
        if r.is_ok() {
            return Err(format!("record accepted: {}", label));
        }
    }
    expect_eq(
        "nothing recorded by the rejected calls",
        &log.stats()?.len(),
        &0,
    )
}

#[test]
fn write_results_writes_the_rendered_text() -> Result<(), String> {
    let path =
        std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("tightness_selftest_results.txt");
    fixture_log()?.write_results(&path)?;
    let back = std::fs::read_to_string(&path).map_err(|e| format!("{:?}: {}", path, e))?;
    expect_eq("file text", &back.as_str(), &PINNED)?;
    let unwritable = std::path::Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("no-such-dir")
        .join("x.txt");
    if fixture_log()?.write_results(&unwritable).is_ok() {
        return Err(String::from("writing into a missing directory succeeded"));
    }
    Ok(())
}

#[test]
fn p99_rank_follows_ceil_99_n_over_100_minus_1_and_refuses_the_unrepresentable(
) -> Result<(), String> {
    for (n, want) in [
        (1usize, 0usize),
        (2, 1),
        (4, 3),
        (50, 49),
        (100, 98),
        (101, 99),
        (200, 197),
    ] {
        expect_eq(&format!("n = {}", n), &p99_rank(n)?, &want)?;
    }
    if p99_rank(0).is_ok() {
        return Err(String::from(
            "p99_rank(0) accepted: there is no p99 of zero samples",
        ));
    }
    if p99_rank(usize::MAX).is_ok() {
        return Err(String::from(
            "p99_rank(usize::MAX) accepted: 99 n overflows",
        ));
    }
    Ok(())
}
