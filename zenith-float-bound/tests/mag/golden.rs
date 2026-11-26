//! CONTRACT 1.3, 14.1, 14.2, 14.3, 14.4; oracle: mpmath.libmp directed rounding (scripts/gen_goldens.py -> tests/goldens/mag_ops.txt); the model checkers are validated against the goldens (accepts them, rejects perturbed ones)
use crate::model::direct;
use crate::model::dyadic::MagV;
use crate::model::oracle::{max_mag, min_pos_mag, pred, succ};
use crate::support::{from_magv, to_magv};
use core::str::SplitWhitespace;

const GOLDEN: &str = include_str!("../goldens/mag_ops.txt");

const OPS: [&str; 8] = [
    "add_up",
    "mul_up",
    "mul_lower",
    "sub_lower",
    "div_up",
    "div_lower",
    "sqrt_up",
    "sqrt_lower",
];

pub(crate) struct Row {
    pub(crate) op: String,
    pub(crate) args: Vec<MagV>,
    pub(crate) want: MagV,
}

pub(crate) fn parse_row(line: &str) -> Result<Option<Row>, String> {
    let line = line.trim();
    if line.is_empty() || line.starts_with('#') {
        return Ok(None);
    }

    let mut it = line.split_whitespace();
    let op = it
        .next()
        .ok_or_else(|| String::from("unexpected end of line"))?;
    if !OPS.contains(&op) {
        return Err(format!("unknown op {}", op));
    }

    let a = next_magv(&mut it)?;
    let args = match op {
        "add_up" | "mul_up" | "mul_lower" | "sub_lower" | "div_up" | "div_lower" => {
            let b = next_magv(&mut it)?;
            vec![a, b]
        }
        "sqrt_up" | "sqrt_lower" => {
            vec![a]
        }
        _ => return Err(format!("unknown op {}", op)),
    };

    let want = next_result(&mut it)?;
    if let Some(extra) = it.next() {
        return Err(format!(
            "unexpected trailing token {:?} after the result of a {} row",
            extra, op
        ));
    }
    Ok(Some(Row {
        op: String::from(op),
        args,
        want,
    }))
}

fn next_i64(it: &mut SplitWhitespace) -> Result<i64, String> {
    it.next()
        .ok_or_else(|| String::from("unexpected end of line"))?
        .parse::<i64>()
        .map_err(|e| format!("failed to parse i64: {}", e))
}

fn next_u64(it: &mut SplitWhitespace) -> Result<u64, String> {
    it.next()
        .ok_or_else(|| String::from("unexpected end of line"))?
        .parse::<u64>()
        .map_err(|e| format!("failed to parse u64: {}", e))
}

fn next_magv(it: &mut SplitWhitespace) -> Result<MagV, String> {
    let exp = next_i64(it)?;
    let mant = next_u64(it)?;
    Ok(MagV::Finite { exp, mant })
}

fn next_result(it: &mut SplitWhitespace) -> Result<MagV, String> {
    let peek = it
        .next()
        .ok_or_else(|| String::from("unexpected end of line"))?;
    if peek == "inf" {
        Ok(MagV::Infinity)
    } else {
        let exp = peek
            .parse::<i64>()
            .map_err(|e| format!("failed to parse i64: {}", e))?;
        let mant = next_u64(it)?;
        Ok(MagV::Finite { exp, mant })
    }
}

pub(crate) fn run_line(line: &str) -> Result<(), String> {
    match parse_row(line)? {
        None => Ok(()),
        Some(row) => {
            let a = row
                .args
                .first()
                .ok_or_else(|| String::from("missing argument a"))?;
            let got_and_want: (MagV, MagV) = match row.op.as_str() {
                "add_up" | "mul_up" | "mul_lower" | "sub_lower" | "div_up" | "div_lower" => {
                    let b = row
                        .args
                        .get(1)
                        .ok_or_else(|| String::from("missing argument b"))?;
                    let ma = from_magv(a)?;
                    let mb = from_magv(b)?;
                    let r = match row.op.as_str() {
                        "add_up" => ma.add_up(&mb),
                        "mul_up" => ma.mul_up(&mb),
                        "mul_lower" => ma.mul_lower(&mb),
                        "sub_lower" => ma.sub_lower(&mb),
                        "div_up" => ma.div_up(&mb),
                        "div_lower" => ma.div_lower(&mb),
                        _ => return Err(format!("unknown op {}", row.op)),
                    };
                    (to_magv(&r)?, row.want)
                }
                "sqrt_up" | "sqrt_lower" => {
                    let ma = from_magv(a)?;
                    let r = if row.op == "sqrt_up" {
                        ma.sqrt_up()
                    } else {
                        ma.sqrt_lower()
                    };
                    (to_magv(&r)?, row.want)
                }
                _ => return Err(format!("unknown op {}", row.op)),
            };

            if got_and_want.0 != got_and_want.1 {
                Err(format!(
                    "{}: expected {:?}, got {:?}",
                    line, got_and_want.1, got_and_want.0
                ))
            } else {
                Ok(())
            }
        }
    }
}

pub(crate) fn verify_row(op: &str, args: &[MagV], r: &MagV) -> Result<(), String> {
    let a = args
        .first()
        .ok_or_else(|| String::from("missing argument a"))?;
    match op {
        "sqrt_up" => crate::sqrt_up::verify(a, r),
        "sqrt_lower" => crate::lower_sqrt::verify(a, r),
        _ => {
            let b = args
                .get(1)
                .ok_or_else(|| String::from("missing argument b"))?;
            match op {
                "add_up" => crate::add_up::verify(a, b, r),
                "mul_up" => crate::mul_up::verify(a, b, r),
                "div_up" => crate::div_up::verify(a, b, r),
                "mul_lower" => crate::lower_mul::verify(a, b, r),
                "sub_lower" => crate::lower_sub::verify(a, b, r),
                "div_lower" => crate::lower_div::verify(a, b, r),
                _ => Err(format!("unknown op {}", op)),
            }
        }
    }
}

fn is_up(op: &str) -> bool {
    matches!(op, "add_up" | "mul_up" | "div_up" | "sqrt_up")
}

pub(crate) fn perturb(op: &str, want: &MagV) -> Option<MagV> {
    if is_up(op) {
        match want {
            MagV::Infinity => Some(max_mag()),
            MagV::Finite { .. } => Some(succ(want)),
        }
    } else {
        match want {
            MagV::Infinity => None,
            MagV::Finite { mant: 0, .. } => Some(min_pos_mag()),
            MagV::Finite { .. } => pred(want),
        }
    }
}

#[test]
fn model_accepts_every_mpmath_golden() -> Result<(), String> {
    for line in GOLDEN.lines() {
        if let Some(row) = parse_row(line)? {
            verify_row(&row.op, &row.args, &row.want).map_err(|e| format!("{}: {}", line, e))?;
        }
    }
    Ok(())
}

#[test]
fn model_rejects_every_perturbed_mpmath_golden() -> Result<(), String> {
    for line in GOLDEN.lines() {
        if let Some(row) = parse_row(line)? {
            if let Some(p) = perturb(&row.op, &row.want) {
                if verify_row(&row.op, &row.args, &p).is_ok() {
                    return Err(format!(
                        "{}: model accepted the perturbed result {:?}",
                        line, p
                    ));
                }
            }
        }
    }
    Ok(())
}

#[test]
fn mag_ops_match_mpmath_goldens() -> Result<(), String> {
    for line in GOLDEN.lines() {
        run_line(line)?;
    }
    Ok(())
}

#[test]
fn golden_file_has_all_eight_ops() -> Result<(), String> {
    for op in OPS.iter() {
        let prefix = format!("{} ", op);
        if !GOLDEN.lines().any(|l| l.starts_with(prefix.as_str())) {
            return Err(format!("missing {}", op));
        }
    }
    Ok(())
}

fn direct_result(op: &str, args: &[MagV]) -> Result<MagV, String> {
    let a = args
        .first()
        .ok_or_else(|| String::from("missing argument a"))?;
    let r = match op {
        "sqrt_up" => direct::sqrt_up(a),
        "sqrt_lower" => direct::sqrt_lower(a),
        _ => {
            let b = args
                .get(1)
                .ok_or_else(|| String::from("missing argument b"))?;
            match op {
                "add_up" => direct::add_up(a, b),
                "mul_up" => direct::mul_up(a, b),
                "mul_lower" => direct::mul_lower(a, b),
                "sub_lower" => direct::sub_lower(a, b),
                "div_up" => direct::div_up(a, b),
                "div_lower" => direct::div_lower(a, b),
                _ => return Err(format!("unknown op {}", op)),
            }
        }
    };
    r.ok_or_else(|| format!("direct oracle failed on {} {:?}", op, args))
}

#[test]
fn direct_integer_oracle_matches_every_mpmath_golden() -> Result<(), String> {
    for line in GOLDEN.lines() {
        if let Some(row) = parse_row(line)? {
            let got = direct_result(&row.op, &row.args)?;
            if got != row.want {
                return Err(format!(
                    "{}: direct oracle gave {:?}, mpmath golden is {:?}",
                    line, got, row.want
                ));
            }
        }
    }
    Ok(())
}

// ---- L0a-6 trap regressions (TEST_PLAN section 12) -------------------------------------------------------------------------
// Trap "golden runner second operand" (commit 0eec0fe): a row's arity is fixed by its op (binary: two operands then the result;
// unary: one operand then the result), and a row with a missing or an extra column is an error, never a silent misparse.
// Trap "round-down overflow dropped from the goldens": the generator once dropped the rows whose exact result is above the largest
// Mag and whose correct round-down result is that largest Mag (`i64::MAX`, `u64::MAX`).

const TOP: &str = "9223372036854775808";

/// The first and the last of the specific rows that must stay in the golden file (operands at the extremes).
const NAMED_OVERFLOW_ROWS: [&str; 2] = [
    "mul_lower 9223372036854775807 18446744073709551615 9223372036854775807 18446744073709551615 9223372036854775807 18446744073709551615",
    "div_lower 9223372036854775807 18446744073709551615 -9223372036854775808 9223372036854775808 9223372036854775807 18446744073709551615",
];

/// Floors on the number of round-down-overflow rows (rows whose result is the largest Mag) at the committed generator output:
/// 512 `mul_lower` and 617 `div_lower` (1,129 together). Regenerating with a wider operand set may only raise them.
const MUL_LOWER_MIN_ROWS: usize = 512;
const DIV_LOWER_MIN_ROWS: usize = 617;

fn is_largest_mag_row(line: &str, op: &str) -> bool {
    matches!(parse_row(line), Ok(Some(r)) if r.op == op && r.want == max_mag())
}

fn count_largest_mag_rows(text: &str, op: &str) -> usize {
    text.lines().filter(|l| is_largest_mag_row(l, op)).count()
}

fn require_round_down_overflow_rows(text: &str) -> Result<(), String> {
    for named in NAMED_OVERFLOW_ROWS.iter() {
        if !text.lines().any(|l| l.trim() == *named) {
            return Err(format!("golden file lacks the row: {}", named));
        }
    }
    let mul = count_largest_mag_rows(text, "mul_lower");
    if mul < MUL_LOWER_MIN_ROWS {
        return Err(format!(
            "only {} mul_lower rows with the largest Mag as result, need at least {}",
            mul, MUL_LOWER_MIN_ROWS
        ));
    }
    let div = count_largest_mag_rows(text, "div_lower");
    if div < DIV_LOWER_MIN_ROWS {
        return Err(format!(
            "only {} div_lower rows with the largest Mag as result, need at least {}",
            div, DIV_LOWER_MIN_ROWS
        ));
    }
    Ok(())
}

fn without(text: &str, drop: impl Fn(&str) -> bool) -> String {
    text.lines()
        .filter(|l| !drop(l))
        .collect::<Vec<&str>>()
        .join("\n")
}

#[test]
fn golden_file_has_the_round_down_overflow_rows() -> Result<(), String> {
    require_round_down_overflow_rows(GOLDEN)
}

#[test]
fn round_down_overflow_check_rejects_a_golden_file_without_those_rows() -> Result<(), String> {
    let no_mul = without(GOLDEN, |l| is_largest_mag_row(l, "mul_lower"));
    let no_div = without(GOLDEN, |l| is_largest_mag_row(l, "div_lower"));
    let no_named = without(GOLDEN, |l| NAMED_OVERFLOW_ROWS.contains(&l.trim()));
    let no_both = without(GOLDEN, |l| {
        is_largest_mag_row(l, "mul_lower") || is_largest_mag_row(l, "div_lower")
    });
    for (what, text) in [
        ("without the mul_lower rows", &no_mul),
        ("without the div_lower rows", &no_div),
        ("without the two named rows", &no_named),
        ("without any round-down-overflow row", &no_both),
    ] {
        if require_round_down_overflow_rows(text).is_ok() {
            return Err(format!("the check accepted a golden file {}", what));
        }
    }
    Ok(())
}

#[test]
fn round_down_overflow_rows_give_the_largest_mag_and_the_library_agrees() -> Result<(), String> {
    // Every such row is also run against the library by mag_ops_match_mpmath_goldens; here: the rows are real overflows, i.e. the
    // exact result is above the largest Mag (the named rows are max*max and max/min-positive).
    for named in NAMED_OVERFLOW_ROWS.iter() {
        run_line(named)?;
    }
    Ok(())
}

fn expect_rejected(line: &str) -> Result<(), String> {
    match parse_row(line) {
        Err(_) => Ok(()),
        Ok(_) => Err(format!("parse_row accepted the malformed row: {}", line)),
    }
}

#[test]
fn parse_row_accepts_well_formed_rows_of_every_arity() -> Result<(), String> {
    let binary = format!("add_up 0 {} 0 {} 1 {}", TOP, TOP, TOP);
    let binary_inf = format!("mul_up 0 {} 0 {} inf", TOP, TOP);
    let unary = format!("sqrt_up 0 {} 0 {}", TOP, TOP);
    let unary_inf = format!("sqrt_lower 0 {} inf", TOP);
    for line in [&binary, &binary_inf, &unary, &unary_inf] {
        match parse_row(line)? {
            Some(_) => {}
            None => return Err(format!("parse_row skipped a data row: {}", line)),
        }
    }
    for line in ["", "   ", "# a comment 1 2 3"] {
        if parse_row(line)?.is_some() {
            return Err(format!("parse_row returned a row for: {:?}", line));
        }
    }
    Ok(())
}

#[test]
fn parse_row_rejects_a_stray_second_operand_on_a_unary_row() -> Result<(), String> {
    // `sqrt_up A B R`: the generator once wrote the binary column layout for a unary op.
    expect_rejected(&format!("sqrt_up 0 {} 0 {} 0 {}", TOP, TOP, TOP))?;
    expect_rejected(&format!("sqrt_lower 0 {} 0 {} inf", TOP, TOP))?;
    expect_rejected(&format!("sqrt_up 0 {} 0 {} 7", TOP, TOP))
}

#[test]
fn parse_row_rejects_a_binary_row_missing_its_second_operand() -> Result<(), String> {
    expect_rejected(&format!("add_up 0 {} 0 {}", TOP, TOP))?;
    expect_rejected(&format!("mul_lower 0 {} inf", TOP))?;
    expect_rejected(&format!("div_up 0 {} 1 {}", TOP, TOP))?;
    expect_rejected(&format!("sub_lower 0 {}", TOP))
}

#[test]
fn parse_row_rejects_trailing_tokens() -> Result<(), String> {
    expect_rejected(&format!("add_up 0 {} 0 {} 1 {} 7", TOP, TOP, TOP))?;
    expect_rejected(&format!("mul_up 0 {} 0 {} inf 7", TOP, TOP))?;
    expect_rejected(&format!("sqrt_up 0 {} 1 {} 7", TOP, TOP))?;
    expect_rejected(&format!("sqrt_lower 0 {} inf inf", TOP))
}

#[test]
fn parse_row_rejects_non_numeric_and_unknown_fields() -> Result<(), String> {
    expect_rejected(&format!("add_up 0 {} 0 {} x {}", TOP, TOP, TOP))?;
    expect_rejected(&format!("add_up 0 {} 0 {} 1 -{}", TOP, TOP, TOP))?;
    expect_rejected(&format!("pow_up 0 {} 0 {} 1 {}", TOP, TOP, TOP))
}
