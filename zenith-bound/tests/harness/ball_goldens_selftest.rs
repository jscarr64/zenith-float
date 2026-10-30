//! COVERS: common::BallOp::Add, common::BallOp::Div, common::BallOp::Mul, common::BallOp::Sqrt, common::BallOp::Sub, common::BallOp::name, common::BallOp::parse, common::Checked::Recomputed, common::Checked::Structural, common::ballrows::check_row, common::ballrows::format_row, common::ballrows::parse_row, common::ballrows::round_int
//! L0b-7: the committed exact-input ball goldens (tests/goldens/ball_exact_ops.txt) are parsed row by row (a malformed row is an `Err`, never a
//! panic) and every row is re-verified by the independent exact-integer checker in tests/common/ballrows.rs. SCOPE, stated plainly: rows with
//! `prec <= 64` and `iprec <= 64` are fully RECOMPUTED (result, sign of zero and inexact flag); the rows with `prec > 64` (65, 127, 128, 129,
//! 256, 1024) only get the STRUCTURAL checks (syntax, canonical form, mantissa widths, zero and sign rules, exponent window), because the
//! checker is `u128` based. Their values are guaranteed only by the generator's own exact-Fraction and mpmath cross-check.
//! The checker's negative self-tests (a perturbed result, flag, sign or mode must be rejected on EVERY recomputed row) and hand rows live here too.

use std::collections::BTreeMap;

use crate::common::ballrows::{
    check_row, format_row, parse_row, round_int, BallOp, Checked, Rounded, MAX_EXACT_PREC,
};
use crate::common::expect_eq;
use crate::common::generators::RoundingMode;
use crate::common::goldens::golden_dir;

const FILE: &str = "ball_exact_ops.txt";
const PREC_LIST: [u32; 18] = [
    1, 2, 3, 4, 8, 24, 31, 32, 33, 53, 63, 64, 65, 127, 128, 129, 256, 1024,
];
const OPS: [BallOp; 5] = [
    BallOp::Add,
    BallOp::Sub,
    BallOp::Mul,
    BallOp::Div,
    BallOp::Sqrt,
];

fn read_file() -> Result<String, String> {
    let path = golden_dir().join(FILE);
    std::fs::read_to_string(&path).map_err(|e| format!("{:?}: {}", path, e))
}

fn op_index(op: &BallOp) -> usize {
    match op {
        BallOp::Add => 0,
        BallOp::Sub => 1,
        BallOp::Mul => 2,
        BallOp::Div => 3,
        BallOp::Sqrt => 4,
    }
}

/// The data lines (every line that is not a `#` header line).
fn rows(text: &str) -> Vec<&str> {
    text.lines().filter(|l| !l.starts_with('#')).collect()
}

fn rm_named(name: &str) -> Result<RoundingMode, String> {
    RoundingMode::parse(name).ok_or_else(|| format!("no mode {}", name))
}

fn verdict(line: &str) -> Result<Checked, String> {
    check_row(&parse_row(line)?)
}

fn must_reject(label: &str, line: &str) -> Result<(), String> {
    match verdict(line) {
        Err(_) => Ok(()),
        Ok(c) => Err(format!("{}: {:?} was accepted as {:?}", label, line, c)),
    }
}

// ---- the committed file -------------------------------------------------------------------------------------------------------

#[test]
fn header_block_is_ten_comment_lines_and_data_follows() -> Result<(), String> {
    let text = read_file()?;
    let mut header = 0usize;
    for l in text.lines() {
        if l.starts_with('#') {
            if header == 0 && !l.contains("generator version 1") {
                return Err(format!("first header line {:?}", l));
            }
            header += 1;
        } else {
            break;
        }
    }
    expect_eq("header lines before the first row", &header, &10)?;
    expect_eq(
        "comment lines in the whole file",
        &text.lines().filter(|l| l.starts_with('#')).count(),
        &10,
    )?;
    let last = text.lines().nth(9).unwrap_or("");
    // provenance: the oracle venv pins (scripts/goldlib.py require_oracle); a 1.2.1 header is stale
    if last != "# generated in the oracle venv (mpmath 1.3.0, sympy 1.14.0, gmpy2 2.3.2, python-flint 0.9.0); see scripts/README.md" {
        return Err(format!("last header line {:?}", last));
    }
    Ok(())
}

#[test]
fn every_row_parses_is_sorted_and_verifies() -> Result<(), String> {
    let text = read_file()?;
    let mut prev: Option<(usize, u32, usize, &str)> = None;
    let mut recomputed = 0usize;
    let mut structural = 0usize;
    let mut per_op: BTreeMap<usize, usize> = BTreeMap::new();
    let mut per_cell: BTreeMap<(usize, u32, usize), usize> = BTreeMap::new();
    for line in rows(&text) {
        let row = parse_row(line).map_err(|e| format!("{:?}: {}", line, e))?;
        let key = (op_index(&row.op), row.prec, row.rm.index(), line);
        if let Some(p) = &prev {
            if *p >= key {
                return Err(format!("not strictly sorted at {:?}", line));
            }
        }
        prev = Some(key);
        match check_row(&row).map_err(|e| format!("{:?}: {}", line, e))? {
            Checked::Recomputed => recomputed += 1,
            Checked::Structural => structural += 1,
        }
        if let Ok(back) = format_row(&row) {
            expect_eq("format_row round trip", &back.as_str(), &line)?;
        }
        *per_op.entry(op_index(&row.op)).or_insert(0) += 1;
        *per_cell
            .entry((op_index(&row.op), row.prec, row.rm.index()))
            .or_insert(0) += 1;
    }
    expect_eq("recomputed rows", &recomputed, &RECOMPUTED)?;
    expect_eq("structural-only rows", &structural, &STRUCTURAL)?;
    let counts: Vec<&usize> = per_op.values().collect();
    let want_counts: Vec<&usize> = PER_OP.iter().collect();
    expect_eq("rows per op (add sub mul div sqrt)", &counts, &want_counts)?;
    // every (op, prec) cell carries the same operand sets under all six modes
    for op in OPS.iter() {
        for p in PREC_LIST.iter() {
            let mut n: Option<usize> = None;
            for rm in RoundingMode::all() {
                let c = match per_cell.get(&(op_index(op), *p, rm.index())) {
                    Some(&c) => c,
                    None => {
                        return Err(format!(
                            "no rows for {} prec {} {}",
                            op.name(),
                            p,
                            rm.name()
                        ))
                    }
                };
                if c == 0 {
                    return Err(String::from("empty cell"));
                }
                match n {
                    None => n = Some(c),
                    Some(first) => expect_eq(
                        &format!("rows of {} prec {} under {}", op.name(), p, rm.name()),
                        &c,
                        &first,
                    )?,
                }
            }
        }
    }
    expect_eq("number of cells", &per_cell.len(), &(5 * 18 * 6))?;
    Ok(())
}

#[test]
fn file_holds_the_directed_hand_rows() -> Result<(), String> {
    let text = read_file()?;
    let data = rows(&text);
    // 9/8 * 9/8 = 81/64 at p = 3 lies between 5/4 and 3/2 (not a tie: 81/64 < 11/8): ToEven, ToOdd, Down and ToZero give 5/4, Up and FromZero 3/2
    for (rm, r) in [
        ("ToEven", "0 1 5"),
        ("ToOdd", "0 1 5"),
        ("Down", "0 1 5"),
        ("ToZero", "0 1 5"),
        ("Up", "0 1 3"),
        ("FromZero", "0 1 3"),
    ] {
        let want = format!("mul 3 {} 4 0 1 9 0 1 9 {} 1", rm, r);
        if !data.contains(&want.as_str()) {
            return Err(format!("missing directed row {:?}", want));
        }
        expect_eq(&want, &verdict(&want)?, &Checked::Recomputed)?;
    }
    // sqrt(81/64) = 9/8 is an exact tie at p = 3 between 1 and 5/4
    for (rm, r) in [("ToEven", "0 1 1"), ("ToOdd", "0 1 5")] {
        let want = format!("sqrt 3 {} 7 0 1 51 {} 1", rm, r);
        if !data.contains(&want.as_str()) {
            return Err(format!("missing tie row {:?}", want));
        }
    }
    Ok(())
}

// ---- negative self-tests of the checker over the whole file -------------------------------------------------------------------

#[test]
fn every_recomputed_row_rejects_perturbed_results_flags_signs_and_modes() -> Result<(), String> {
    let text = read_file()?;
    let mut tried = [0usize; 5];
    for line in rows(&text) {
        let base = parse_row(line)?;
        if check_row(&base)? != Checked::Recomputed {
            continue;
        }
        // 1. the inexact flag flipped
        let mut r = parse_row(line)?;
        r.inexact = !r.inexact;
        if check_row(&r).is_ok() {
            return Err(format!("flipped flag accepted: {:?}", line));
        }
        tried[0] += 1;
        // 2. the result sign flipped (also flips the sign of a zero result)
        let mut r = parse_row(line)?;
        r.r.neg = !r.r.neg;
        if check_row(&r).is_ok() {
            return Err(format!("flipped result sign accepted: {:?}", line));
        }
        tried[1] += 1;
        // 3. the result exponent moved by one (expected_mid perturbed)
        if base.r.bits != 0 {
            let mut r = parse_row(line)?;
            r.r.exp += 1;
            if check_row(&r).is_ok() {
                return Err(format!("moved result exponent accepted: {:?}", line));
            }
            tried[2] += 1;
        }
        // 4. Up <-> Down and ToZero <-> FromZero on an inexact row give different results, so the wrong mode must be rejected
        if base.inexact {
            for (from, to) in [
                ("Up", "Down"),
                ("Down", "Up"),
                ("ToZero", "FromZero"),
                ("FromZero", "ToZero"),
            ] {
                if base.rm.name() == from {
                    let mut r = parse_row(line)?;
                    r.rm = rm_named(to)?;
                    if check_row(&r).is_ok() {
                        return Err(format!("wrong mode {} accepted: {:?}", to, line));
                    }
                    tried[3] += 1;
                }
            }
        }
        // 5. the result mantissa changed to a different canonical odd value (kept within prec)
        if let Some(m) = base.r.mant {
            if m > 1 {
                let mut r = parse_row(line)?;
                r.r.mant = Some(m - 2);
                if r.r.bits == base.r.bits
                    && (m - 2) >= (1u128 << (base.r.bits - 1))
                    && check_row(&r).is_ok()
                {
                    return Err(format!("changed result mantissa accepted: {:?}", line));
                }
                tried[4] += 1;
            }
        }
    }
    for (i, t) in tried.iter().enumerate() {
        if *t < 1000 {
            return Err(format!("negative self-test {} ran only {} times", i + 1, t));
        }
    }
    Ok(())
}

// ---- hand rows: values, zero signs, malformed text ----------------------------------------------------------------------------

#[test]
fn hand_rows_with_known_values_are_recomputed_ok() -> Result<(), String> {
    for line in [
        // 1 + 1/8 at p = 3 is a tie between 1 and 5/4
        "add 3 ToEven 3 0 1 1 0 -2 1 0 1 1 1",
        "add 3 ToOdd 3 0 1 1 0 -2 1 0 1 5 1",
        "add 3 Up 3 0 1 1 0 -2 1 0 1 5 1",
        "add 3 Down 3 0 1 1 0 -2 1 0 1 1 1",
        // 7/8 + 1/8 = 1 exactly (a carry)
        "add 3 ToEven 3 0 0 7 0 -2 1 0 1 1 0",
        // 1/3 at p = 2
        "div 2 ToEven 2 0 1 1 0 2 3 0 -1 3 1",
        "div 2 Down 2 0 1 1 0 2 3 0 -1 1 1",
        // sqrt 2 at p = 3
        "sqrt 3 ToEven 3 0 2 1 0 1 3 1",
        "sqrt 3 Down 3 0 2 1 0 1 5 1",
        // 3 - 1 = 2 exactly at p = 1
        "sub 1 ToEven 2 0 2 3 0 1 1 0 2 1 0",
        // 1 - 1/8 = 7/8 exactly at p = 3, and with the 3-bit grid at p = 2: 7/8 is a tie between 3/4 and 1
        "sub 3 ToEven 3 0 1 1 0 -2 1 0 0 7 0",
        "sub 2 ToEven 3 0 1 1 0 -2 1 0 1 1 1",
        "sub 2 ToOdd 3 0 1 1 0 -2 1 0 0 3 1",
        // a huge gap: 1 + 2^-30000 rounds back to 1 or one tick up
        "add 8 ToEven 1 0 1 1 0 -29999 1 0 1 1 1",
        "add 8 Up 1 0 1 1 0 -29999 1 0 1 81 1",
        "sub 8 Down 1 0 1 1 0 -29999 1 0 0 ff 1",
        // exact products and quotients
        "mul 4 ToEven 4 0 1 9 1 1 1 1 1 9 0",
        "div 4 ToEven 4 0 1 9 1 1 1 1 1 9 0",
    ] {
        expect_eq(line, &verdict(line)?, &Checked::Recomputed)?;
    }
    Ok(())
}

#[test]
fn zero_results_follow_the_ieee_sign_rules() -> Result<(), String> {
    for line in [
        "add 4 ToEven 4 0 0 0 1 0 0 0 0 0 0",
        "add 4 Down 4 0 0 0 1 0 0 1 0 0 0",
        "add 4 Up 4 1 0 0 0 0 0 0 0 0 0",
        "add 4 ToEven 4 1 0 0 1 0 0 1 0 0 0",
        "add 4 Down 4 0 0 0 0 0 0 0 0 0 0",
        "add 4 ToEven 4 0 1 1 1 1 1 0 0 0 0",
        "add 4 Down 4 0 1 1 1 1 1 1 0 0 0",
        "sub 4 ToEven 4 0 1 1 0 1 1 0 0 0 0",
        "sub 4 Down 4 0 1 1 0 1 1 1 0 0 0",
        "sub 4 ToEven 4 0 0 0 0 0 0 0 0 0 0",
        "sub 4 ToEven 4 1 0 0 0 0 0 1 0 0 0",
        "mul 4 ToEven 4 1 0 0 0 1 1 1 0 0 0",
        "mul 4 Down 4 1 0 0 1 0 0 0 0 0 0",
        "div 4 ToEven 4 0 0 0 1 1 1 1 0 0 0",
        "div 4 ToEven 4 1 0 0 1 1 1 0 0 0 0",
        "sqrt 4 ToEven 4 1 0 0 1 0 0 0",
        "sqrt 4 ToEven 4 0 0 0 0 0 0 0",
    ] {
        expect_eq(line, &verdict(line)?, &Checked::Recomputed)?;
    }
    for (label, line) in [
        (
            "+0 + -0 under Down must be -0",
            "add 4 Down 4 0 0 0 1 0 0 0 0 0 0",
        ),
        (
            "+0 + -0 under ToEven must be +0",
            "add 4 ToEven 4 0 0 0 1 0 0 1 0 0 0",
        ),
        (
            "x + -x under Down must be -0",
            "add 4 Down 4 0 1 1 1 1 1 0 0 0 0",
        ),
        (
            "x + -x under ToEven must be +0",
            "add 4 ToEven 4 0 1 1 1 1 1 1 0 0 0",
        ),
        ("-0 + -0 must be -0", "add 4 ToEven 4 1 0 0 1 0 0 0 0 0 0"),
        ("+0 * -1 must be -0", "mul 4 ToEven 4 0 0 0 1 1 1 0 0 0 0"),
        ("-0 / -1 must be +0", "div 4 ToEven 4 1 0 0 1 1 1 1 0 0 0"),
        ("sqrt(-0) must be -0", "sqrt 4 ToEven 4 1 0 0 0 0 0 0"),
        ("an inexact zero", "add 4 ToEven 4 0 0 0 0 0 0 0 0 0 1"),
    ] {
        must_reject(label, line)?;
    }
    Ok(())
}

#[test]
fn malformed_and_wrong_rows_are_errors_not_panics() -> Result<(), String> {
    for (label, line) in [
        ("empty row", ""),
        ("unknown op", "pow 3 ToEven 3 0 1 1 0 1 1 0 1 1 0"),
        ("too few fields", "add 3 ToEven 3 0 1 1 0 1 1 0 1 1"),
        ("too many fields", "add 3 ToEven 3 0 1 1 0 1 1 0 1 1 0 0"),
        (
            "sqrt with a b operand",
            "sqrt 3 ToEven 3 0 1 1 0 1 1 0 1 1 0",
        ),
        ("repeated space", "add  3 ToEven 3 0 1 1 0 1 1 0 1 1 0"),
        ("tab separator", "add 3\tToEven 3 0 1 1 0 1 1 0 1 1 0"),
        ("unknown mode", "add 3 Nearest 3 0 1 1 0 1 1 0 1 1 0"),
        ("None mode", "add 3 None 3 0 1 1 0 1 1 0 1 1 0"),
        ("prec 0", "add 0 ToEven 3 0 1 1 0 1 1 0 1 1 0"),
        ("plus sign in prec", "add +3 ToEven 3 0 1 1 0 1 1 0 1 1 0"),
        (
            "leading zero in prec",
            "add 03 ToEven 3 0 1 1 0 1 1 0 1 1 0",
        ),
        ("bad sign", "add 3 ToEven 3 2 1 1 0 1 1 0 1 1 0"),
        (
            "negative zero exponent",
            "add 3 ToEven 3 0 -0 1 0 1 1 0 1 1 0",
        ),
        (
            "exponent out of bound",
            "add 3 ToEven 3 0 9999999 1 0 1 1 0 1 1 0",
        ),
        (
            "exponent not a number",
            "add 3 ToEven 3 0 x 1 0 1 1 0 1 1 0",
        ),
        (
            "exponent overflowing i64",
            "add 3 ToEven 3 0 99999999999999999999 1 0 1 1 0 1 1 0",
        ),
        ("uppercase hex", "add 3 ToEven 8 0 1 B 0 1 1 0 1 1 0"),
        ("even mantissa", "add 3 ToEven 3 0 1 2 0 1 1 0 1 1 0"),
        ("leading zero hex", "add 3 ToEven 3 0 1 01 0 1 1 0 1 1 0"),
        ("non-hex digit", "add 3 ToEven 3 0 1 g 0 1 1 0 1 1 0"),
        (
            "zero with an exponent",
            "add 3 ToEven 3 0 5 0 0 1 1 0 1 1 0",
        ),
        ("bad inexact flag", "add 3 ToEven 3 0 1 1 0 1 1 0 1 1 2"),
        (
            "operand wider than iprec",
            "add 3 ToEven 2 0 1 7 0 1 1 0 1 1 0",
        ),
        (
            "result wider than prec",
            "add 2 ToEven 3 0 1 1 0 1 1 0 1 7 0",
        ),
        ("inexact zero result", "mul 3 ToEven 3 0 0 0 0 1 1 0 0 0 1"),
        ("division by zero", "div 3 ToEven 3 0 1 1 0 0 0 0 1 1 0"),
        ("sqrt of a negative", "sqrt 3 ToEven 3 1 1 1 0 1 1 0"),
        (
            "mul zero result from nonzero operands",
            "mul 3 ToEven 3 0 1 1 0 1 1 0 0 0 0",
        ),
        ("mul result sign", "mul 3 ToEven 3 0 1 1 0 1 1 1 1 1 0"),
        (
            "mul result exponent window",
            "mul 3 ToEven 3 0 1 1 0 1 1 0 5 1 0",
        ),
        ("div result sign", "div 3 ToEven 3 0 1 1 1 1 1 0 1 1 0"),
        (
            "div result exponent window",
            "div 3 ToEven 3 0 1 1 0 1 1 0 9 1 0",
        ),
        (
            "sqrt result exponent window",
            "sqrt 3 ToEven 3 0 2 1 0 9 1 0",
        ),
        (
            "add result exponent above both operands",
            "add 3 ToEven 3 0 1 1 0 1 1 0 9 1 0",
        ),
        (
            "wrong value: 1 + 1 is not 3",
            "add 3 ToEven 3 0 1 1 0 1 1 0 2 3 0",
        ),
        ("wrong tie direction", "add 3 ToEven 3 0 1 1 0 -2 1 0 1 5 1"),
        (
            "exact result flagged inexact",
            "add 3 ToEven 3 0 0 7 0 -2 1 0 1 1 1",
        ),
        (
            "inexact result flagged exact",
            "add 3 ToEven 3 0 1 1 0 -2 1 0 1 1 0",
        ),
    ] {
        must_reject(label, line)?;
    }
    Ok(())
}

#[test]
fn rows_wider_than_64_bits_get_only_the_structural_checks() -> Result<(), String> {
    // a p = 65 row is not recomputed (u128 limit): a consistent one is accepted as Structural, and so is a WRONG value (that is the stated limit)
    let ok = "mul 65 ToEven 4 0 1 9 0 1 9 0 1 51 0";
    expect_eq(ok, &verdict(ok)?, &Checked::Structural)?;
    let wrong_value = "mul 65 ToEven 4 0 1 9 0 1 9 0 1 53 1";
    expect_eq(wrong_value, &verdict(wrong_value)?, &Checked::Structural)?;
    // but the structural rules still bite
    must_reject("p=128 result sign", "mul 128 ToEven 4 0 1 9 0 1 9 1 1 51 0")?;
    must_reject(
        "p=128 result wider than prec",
        &format!("mul 8 ToEven 4 0 1 9 0 1 9 0 1 {}1 0", "f".repeat(3)),
    )?;
    must_reject(
        "p=1024 sqrt of a negative",
        "sqrt 1024 ToEven 4 1 1 9 0 1 9 0",
    )?;
    // a mantissa wider than 128 bits parses (as a value without u128 form) and is accepted structurally only under a wide prec
    let wide = format!("mul 256 ToEven 4 0 1 9 0 1 9 0 1 {}1 1", "f".repeat(40));
    expect_eq("wide mantissa", &verdict(&wide)?, &Checked::Structural)?;
    if format_row(&parse_row(&wide)?).is_ok() {
        return Err(String::from(
            "format_row printed a mantissa wider than 128 bits",
        ));
    }
    must_reject(
        "wide mantissa over prec",
        &format!("mul 64 ToEven 4 0 1 9 0 1 9 0 1 {}1 1", "f".repeat(40)),
    )?;
    Ok(())
}

#[test]
fn round_int_hand_cases() -> Result<(), String> {
    let r = |m: u128, t: i64, inexact: bool| Rounded { m, t, inexact };
    // 0b1011 at p = 3 is the tie 5.5 units: ToEven takes 6 (5 is odd), ToOdd keeps 5
    for (neg, rm, want) in [
        (false, "ToEven", r(3, 2, true)),
        (false, "ToOdd", r(5, 1, true)),
        (false, "Up", r(3, 2, true)),
        (false, "Down", r(5, 1, true)),
        (false, "ToZero", r(5, 1, true)),
        (false, "FromZero", r(3, 2, true)),
        (true, "Up", r(5, 1, true)),
        (true, "Down", r(3, 2, true)),
        (true, "ToZero", r(5, 1, true)),
        (true, "FromZero", r(3, 2, true)),
    ] {
        let got = round_int(neg, 0b1011, false, 0, 3, &rm_named(rm)?)?;
        expect_eq(&format!("1011 neg={} {}", neg, rm), &got, &want)?;
    }
    // exact: nothing to round
    expect_eq(
        "exact",
        &round_int(false, 0b1010, false, 0, 3, &rm_named("Up")?)?,
        &r(5, 1, false),
    )?;
    // a sticky fraction below half, and above half
    expect_eq(
        "sticky below",
        &round_int(false, 0b1010, true, 0, 3, &rm_named("ToEven")?)?,
        &r(5, 1, true),
    )?;
    expect_eq(
        "sticky below up",
        &round_int(false, 0b1010, true, 0, 3, &rm_named("Up")?)?,
        &r(3, 2, true),
    )?;
    expect_eq(
        "sticky above ToOdd",
        &round_int(false, 0b1011, true, 0, 3, &rm_named("ToOdd")?)?,
        &r(3, 2, true),
    )?;
    expect_eq(
        "sticky above ToEven",
        &round_int(false, 0b1011, true, 0, 3, &rm_named("ToEven")?)?,
        &r(3, 2, true),
    )?;
    expect_eq(
        "sticky above Down",
        &round_int(false, 0b1011, true, 0, 3, &rm_named("Down")?)?,
        &r(5, 1, true),
    )?;
    // a carry: 0b111 at p = 2 is the tie 3.5 units
    expect_eq(
        "carry ToEven",
        &round_int(false, 0b111, false, 0, 2, &rm_named("ToEven")?)?,
        &r(1, 3, true),
    )?;
    expect_eq(
        "tie ToOdd",
        &round_int(false, 0b111, false, 0, 2, &rm_named("ToOdd")?)?,
        &r(3, 1, true),
    )?;
    // p = 1: 3 is the tie between 2 and 4, the lower neighbour has an odd unit count
    expect_eq(
        "p1 ToEven",
        &round_int(false, 3, false, 0, 1, &rm_named("ToEven")?)?,
        &r(1, 2, true),
    )?;
    expect_eq(
        "p1 ToOdd",
        &round_int(false, 3, false, 0, 1, &rm_named("ToOdd")?)?,
        &r(1, 1, true),
    )?;
    // zero and canonicalisation
    expect_eq(
        "zero",
        &round_int(true, 0, false, 77, 5, &rm_named("Down")?)?,
        &r(0, 0, false),
    )?;
    expect_eq(
        "trailing zeros",
        &round_int(false, 0b1100, false, 3, 8, &rm_named("Down")?)?,
        &r(3, 5, false),
    )?;
    // errors
    for (label, res) in [
        ("prec 0", round_int(false, 5, false, 0, 0, &rm_named("Up")?)),
        (
            "prec above the limit",
            round_int(false, 5, false, 0, MAX_EXACT_PREC + 1, &rm_named("Up")?),
        ),
        (
            "fraction without a unit",
            round_int(false, 0, true, 0, 3, &rm_named("Up")?),
        ),
        (
            "fraction with too few bits",
            round_int(false, 0b101, true, 0, 3, &rm_named("Up")?),
        ),
    ] {
        if res.is_ok() {
            return Err(format!("round_int accepted {}", label));
        }
    }
    Ok(())
}

#[test]
fn op_variants_have_the_file_order_and_names() -> Result<(), String> {
    for (op, i, name) in [
        (BallOp::Add, 0, "add"),
        (BallOp::Sub, 1, "sub"),
        (BallOp::Mul, 2, "mul"),
        (BallOp::Div, 3, "div"),
        (BallOp::Sqrt, 4, "sqrt"),
    ] {
        expect_eq(name, &op_index(&op), &i)?;
        expect_eq(name, &op.name(), &name)?;
    }
    Ok(())
}

#[test]
fn op_names_round_trip() -> Result<(), String> {
    for op in OPS.iter() {
        expect_eq("parse(name)", &BallOp::parse(op.name()).as_ref(), &Some(op))?;
    }
    expect_eq("unknown op", &BallOp::parse("pow"), &None)?;
    expect_eq("upper case", &BallOp::parse("ADD"), &None)
}

const RECOMPUTED: usize = 25_044;
const STRUCTURAL: usize = 12_756;
const PER_OP: [usize; 5] = [9_222, 9_216, 6_564, 6_330, 6_468];
