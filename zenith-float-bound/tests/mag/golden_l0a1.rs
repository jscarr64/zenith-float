//! TEST_PLAN 21 L0a-1 remainder; CONTRACT 1.3, 14.1 to 14.4; oracle: mpmath.libmp directed rounding (scripts/gen_goldens.py ->
//! tests/goldens/mag_ops_l0a1.txt): i64 +-3 exponent neighbours, add_up/sub_lower exponent gaps, sqrt radicand exponent parity
//! COVERS: Mag::add_up, Mag::sub_lower, Mag::mul_up, Mag::mul_lower, Mag::div_up, Mag::div_lower, Mag::sqrt_up, Mag::sqrt_lower
use crate::golden::{parse_row, perturb, run_line, verify_row, Row};
use crate::model::dyadic::MagV;

const ROWS: &str = include_str!("../goldens/mag_ops_l0a1.txt");

fn rows(text: &str) -> Result<Vec<Row>, String> {
    let mut out = Vec::new();
    for line in text.lines() {
        if let Some(r) = parse_row(line)? {
            out.push(r);
        }
    }
    Ok(out)
}

fn exp_of(m: Option<&MagV>) -> Option<i64> {
    match m {
        Some(MagV::Finite { exp, mant }) if *mant != 0 => Some(*exp),
        _ => None,
    }
}

#[test]
fn l0a1_rows_match_the_library() -> Result<(), String> {
    for line in ROWS.lines() {
        run_line(line)?;
    }
    Ok(())
}

#[test]
fn model_accepts_every_l0a1_row() -> Result<(), String> {
    for r in rows(ROWS)?.iter() {
        verify_row(&r.op, &r.args, &r.want).map_err(|e| format!("{} {:?}: {}", r.op, r.args, e))?;
    }
    Ok(())
}

#[test]
fn model_rejects_every_perturbed_l0a1_row() -> Result<(), String> {
    for r in rows(ROWS)?.iter() {
        if let Some(p) = perturb(&r.op, &r.want) {
            if verify_row(&r.op, &r.args, &p).is_ok() {
                return Err(format!(
                    "{} {:?}: model accepted the perturbed result {:?}",
                    r.op, r.args, p
                ));
            }
        }
    }
    Ok(())
}

/// Every op named in COVERS has golden rows here, so l0a1_rows_match_the_library really runs it (an op whose rows vanished from the
/// regenerated file would otherwise still read as covered).
#[test]
fn l0a1_rows_exercise_every_covered_op() -> Result<(), String> {
    let all = rows(ROWS)?;
    for op in [
        "add_up",
        "sub_lower",
        "mul_up",
        "mul_lower",
        "div_up",
        "div_lower",
        "sqrt_up",
        "sqrt_lower",
    ] {
        if !all.iter().any(|r| r.op == op) {
            return Err(format!("mag_ops_l0a1.txt has no {} row", op));
        }
    }
    Ok(())
}

const NB3: [i64; 4] = [i64::MIN + 2, i64::MIN + 3, i64::MAX - 3, i64::MAX - 2];

/// The coverage the L0a-1 remainder asks for; returns the first missing item.
fn missing(text: &str) -> Result<Option<String>, String> {
    let all = rows(text)?;
    for op in [
        "add_up",
        "mul_up",
        "mul_lower",
        "sub_lower",
        "div_up",
        "div_lower",
        "sqrt_up",
        "sqrt_lower",
    ] {
        for e in NB3 {
            let hit = all
                .iter()
                .any(|r| r.op == op && r.args.iter().any(|a| exp_of(Some(a)) == Some(e)));
            if !hit {
                return Ok(Some(format!("{}: no row with exponent {}", op, e)));
            }
        }
    }
    for op in ["add_up", "sub_lower"] {
        for gap in [1_i64, 2, 62, 63, 64, 65, 66, 127, 128, 129, 200] {
            for top in [100_i64, -100, i64::MIN + 300] {
                let hit = all.iter().any(|r| {
                    r.op == op
                        && exp_of(r.args.first()) == Some(top)
                        && exp_of(r.args.get(1)).and_then(|b| top.checked_sub(b)) == Some(gap)
                });
                if !hit {
                    return Ok(Some(format!(
                        "{}: no row with a = 2^{} scale and an exponent gap of {}",
                        op, top, gap
                    )));
                }
            }
        }
    }
    for op in ["sqrt_up", "sqrt_lower"] {
        for e in [
            -3_i64,
            -2,
            -1,
            1,
            2,
            3,
            i64::MIN,
            i64::MIN + 1,
            i64::MIN + 2,
            i64::MIN + 3,
            i64::MAX - 3,
            i64::MAX - 2,
            i64::MAX - 1,
            i64::MAX,
        ] {
            let hit = all
                .iter()
                .any(|r| r.op == op && exp_of(r.args.first()) == Some(e));
            if !hit {
                return Ok(Some(format!(
                    "{}: no radicand with exponent {} (parity)",
                    op, e
                )));
            }
        }
    }
    Ok(None)
}

#[test]
fn l0a1_file_has_neighbour_gap_and_parity_rows() -> Result<(), String> {
    match missing(ROWS)? {
        None => Ok(()),
        Some(m) => Err(m),
    }
}

#[test]
fn l0a1_coverage_check_rejects_a_file_without_the_gap_rows() -> Result<(), String> {
    let cut: String = ROWS
        .lines()
        .filter(|l| {
            !(l.starts_with("add_up ")
                && l.split_whitespace().nth(1) == Some("100")
                && l.split_whitespace().nth(3) == Some("35"))
        })
        .map(|l| format!("{}\n", l))
        .collect();
    match missing(&cut)? {
        Some(_) => Ok(()),
        None => Err(String::from(
            "coverage check accepted a file without the add_up gap-65 rows at exponent 100",
        )),
    }
}
