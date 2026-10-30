//! CONTRACT 14.1 to 14.4, 14.7, 3.2 to 3.4; oracle: an independent big-integer script (scripts/gen_review_goldens.py -> tests/goldens/mag_ops_review.txt), not mpmath and not the Rust model
//! COVERS: Mag::add_up, Mag::div_lower, Mag::div_up, Mag::half_ulp, Mag::mul_lower, Mag::mul_up, Mag::sub_lower, Mag::ulp_mag

use crate::model::dyadic::MagV;
use crate::support::{from_magv, to_magv};
use core::str::SplitWhitespace;
use zenith_bound::mag::Mag;

const ROWS: &str = include_str!("../goldens/mag_ops_review.txt");

fn word<'a>(it: &mut SplitWhitespace<'a>) -> Result<&'a str, String> {
    it.next()
        .ok_or_else(|| String::from("unexpected end of line"))
}

fn i64_of(it: &mut SplitWhitespace) -> Result<i64, String> {
    word(it)?
        .parse::<i64>()
        .map_err(|e| format!("bad i64: {}", e))
}

fn u64_of(it: &mut SplitWhitespace) -> Result<u64, String> {
    word(it)?
        .parse::<u64>()
        .map_err(|e| format!("bad u64: {}", e))
}

fn finite(it: &mut SplitWhitespace) -> Result<MagV, String> {
    let exp = i64_of(it)?;
    let mant = u64_of(it)?;
    Ok(MagV::Finite { exp, mant })
}

fn result(it: &mut SplitWhitespace) -> Result<MagV, String> {
    let first = word(it)?;
    if first == "inf" {
        return Ok(MagV::Infinity);
    }
    let exp = first
        .parse::<i64>()
        .map_err(|e| format!("bad i64: {}", e))?;
    let mant = u64_of(it)?;
    Ok(MagV::Finite { exp, mant })
}

/// Runs one row if its op is `op`; returns whether it ran.
fn run_line(line: &str, op: &str) -> Result<bool, String> {
    let line = line.trim();
    if line.is_empty() || line.starts_with('#') {
        return Ok(false);
    }
    let mut it = line.split_whitespace();
    if word(&mut it)? != op {
        return Ok(false);
    }
    let (got, want) = match op {
        "ulp_mag" | "half_ulp" => {
            let exp = i64_of(&mut it)?;
            let prec = usize::try_from(u64_of(&mut it)?).map_err(|e| format!("bad prec: {}", e))?;
            let m = if op == "ulp_mag" {
                Mag::ulp_mag(exp, prec)
            } else {
                Mag::half_ulp(exp, prec)
            };
            (to_magv(&m)?, result(&mut it)?)
        }
        _ => {
            let a = from_magv(&finite(&mut it)?)?;
            let b = from_magv(&finite(&mut it)?)?;
            let m = match op {
                "add_up" => a.add_up(&b),
                "mul_up" => a.mul_up(&b),
                "mul_lower" => a.mul_lower(&b),
                "sub_lower" => a.sub_lower(&b),
                "div_up" => a.div_up(&b),
                "div_lower" => a.div_lower(&b),
                _ => return Err(format!("unknown op {}", op)),
            };
            (to_magv(&m)?, result(&mut it)?)
        }
    };
    if let Some(extra) = it.next() {
        return Err(format!(
            "unexpected trailing token {:?} after the result of a {} row",
            extra, op
        ));
    }
    if got != want {
        return Err(format!("{}: expected {:?}, got {:?}", line, want, got));
    }
    Ok(true)
}

/// Runs every row whose op is `op`; returns how many rows ran.
fn run(op: &str) -> Result<usize, String> {
    let mut ran = 0usize;
    for line in ROWS.lines() {
        if run_line(line, op)? {
            ran += 1;
        }
    }
    if ran == 0 {
        Err(format!("no rows for {}", op))
    } else {
        Ok(ran)
    }
}

#[test]
fn review_rows_sub_lower() -> Result<(), String> {
    run("sub_lower").map(|_| ())
}

#[test]
fn review_rows_add_up() -> Result<(), String> {
    run("add_up").map(|_| ())
}

#[test]
fn review_rows_div_up() -> Result<(), String> {
    run("div_up").map(|_| ())
}

#[test]
fn review_rows_div_lower() -> Result<(), String> {
    run("div_lower").map(|_| ())
}

#[test]
fn review_rows_mul_up() -> Result<(), String> {
    run("mul_up").map(|_| ())
}

#[test]
fn review_rows_mul_lower() -> Result<(), String> {
    run("mul_lower").map(|_| ())
}

#[test]
fn review_rows_ulp_mag() -> Result<(), String> {
    run("ulp_mag").map(|_| ())
}

#[test]
fn review_rows_half_ulp() -> Result<(), String> {
    run("half_ulp").map(|_| ())
}

const REVIEW_OPS: [&str; 8] = [
    "sub_lower",
    "add_up",
    "div_up",
    "div_lower",
    "mul_up",
    "mul_lower",
    "ulp_mag",
    "half_ulp",
];

fn first_row_of(op: &str) -> Result<&'static str, String> {
    ROWS.lines()
        .find(|l| l.split_whitespace().next() == Some(op))
        .ok_or_else(|| format!("no {} row in the review goldens", op))
}

fn expect_row_rejected(line: &str, op: &str) -> Result<(), String> {
    match run_line(line, op) {
        Err(_) => Ok(()),
        Ok(_) => Err(format!("the review-row runner accepted: {}", line)),
    }
}

#[test]
fn review_row_runner_accepts_the_first_row_of_every_op() -> Result<(), String> {
    for op in REVIEW_OPS {
        let line = first_row_of(op)?;
        if !run_line(line, op)? {
            return Err(format!("{} row was not run: {}", op, line));
        }
    }
    Ok(())
}

#[test]
fn review_row_runner_rejects_trailing_tokens_on_binary_and_ulp_rows() -> Result<(), String> {
    for op in REVIEW_OPS {
        let line = first_row_of(op)?;
        expect_row_rejected(&format!("{} 7", line), op)?;
        expect_row_rejected(&format!("{} inf", line), op)?;
    }
    Ok(())
}

#[test]
fn review_row_runner_rejects_rows_missing_a_column() -> Result<(), String> {
    for op in REVIEW_OPS {
        let line = first_row_of(op)?;
        let cut = line
            .rsplit_once(' ')
            .map(|(head, _)| head)
            .ok_or_else(|| format!("row has no space: {}", line))?;
        expect_row_rejected(cut, op)?;
    }
    Ok(())
}
