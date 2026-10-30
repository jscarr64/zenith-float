//! Harper L0b: the slice oracle rows (tests/goldens/slices/arith/arith_oracle_rows.txt, scripts/gen_arith_oracle.py).
//! Each row is an input ball pair at PREC and the exact image rounded outward to W = 2·PREC + 64 bits ([LO, HI],
//! cross-checked with mpmath.iv). Checks per row: [LO, HI] ⊆ result (containment); rad ≤ τ·RREF (tight); rad ≤ 2^10·RREF
//! (loose), RREF = hull spread + U (TEST_PLAN 9, placeholder τ). Radii compare as exact numbers; a radius that does not
//! convert, or a row that does not parse, fails the test (never skipped).

use core::cmp::Ordering;

use zenith_bound::ball::Ball;
use zenith_bound::mag::Mag;
use zenith_bound::policy::DomainPolicy;
use zenith_float_num::{ExactNum, EXPONENT_MIN};

use super::{num, ok, pow2, row_val, want_contains};
use crate::common::goldens::golden_dir;

/// Rows the generator writes (v3): 13 random ops × 124 (add sub mul mul_exact sqr powi sqrt neg div inv hypot abs mul_2exp),
/// then the fixed sections: 36 powi rows, 32 mul_2exp exponent-edge rows, 22 hypot edge rows, 2 mul_2exp_unf_rad rows
/// (mid 0, rad 2^-53, k = i64::MIN; appended last). A short file fails.
const WANT_RANDOM: usize = 1612;
const WANT_POWI_FIXED: usize = 36;
const WANT_ROWS: usize = 1704;
/// Exact row count per OP tag (every row is checked; a tag count off by one fails).
const WANT_PER_OP: [(&str, usize); 16] = [
    ("add", 124),
    ("sub", 124),
    ("mul", 124),
    ("mul_exact", 124),
    ("sqr", 124),
    ("powi", 160),
    ("sqrt", 124),
    ("neg", 124),
    ("div", 124),
    ("inv", 124),
    ("hypot", 146),
    ("abs", 124),
    ("mul_2exp", 140),
    ("mul_2exp_ovf", 8),
    ("mul_2exp_unf", 8),
    ("mul_2exp_unf_rad", 2),
];
const LAMBDA: i64 = 1024;

fn field<'a>(v: &[&'a str], i: usize) -> Result<&'a str, String> {
    match v.get(i) {
        Some(&s) => Ok(s),
        None => Err(format!("missing field {}", i)),
    }
}

fn value(v: &[&str], at: usize) -> Result<ExactNum, String> {
    let s = field(v, at)?;
    let e = field(v, at.checked_add(1).ok_or("at")?)?;
    let m = field(v, at.checked_add(2).ok_or("at")?)?;
    row_val(s, e, m)
}

/// Whether the finite-or-infinite `rad` is at least `|field at| · 2^n`, comparing exponents in i128 (the scaled truth
/// lies below the ExactNum and Mag ranges). A zero truth endpoint is an error: these rows carry the nonzero mpmath image.
fn rad_covers_scaled(rad: &Mag, v: &[&str], at: usize, n: i64) -> Result<bool, String> {
    let e = field(v, at.checked_add(1).ok_or("at")?)?;
    let m = field(v, at.checked_add(2).ok_or("at")?)?;
    if m == "0" {
        return Err("truth endpoint is 0; want the nonzero mpmath image".into());
    }
    let (Some(rexp), Some(rmant)) = (rad.exp(), rad.mant()) else {
        return Ok(rad.is_infinity());
    };
    let te: i64 = e.parse().map_err(|_| format!("exp {}", e))?;
    let et = i128::from(te)
        .checked_add(i128::from(n))
        .ok_or("scaled exp")?;
    let er = i128::from(rexp);
    if er != et {
        return Ok(er > et);
    }
    // Same exponent: compare the fractions 0.MANT in [1/2, 1).
    let fr = ExactNum::from(rmant).mul_full_prec(&pow2(-64)?);
    let fh = row_val("0", "0", m)?;
    Ok(fr.partial_cmp(&fh) != Some(Ordering::Less))
}

/// Placeholder τ (TEST_PLAN 9): add/sub 2, mul/mul_exact/sqr/sqrt/div/inv/hypot 4, powi 2n, neg/abs/mul_2exp 1.
fn tau(op: &str, n: i64) -> Result<i64, String> {
    match op {
        "add" | "sub" => Ok(2),
        "mul" | "mul_exact" | "sqr" | "sqrt" | "div" | "inv" | "hypot" => Ok(4),
        "neg" | "abs" | "mul_2exp" => Ok(1),
        "powi" => n
            .checked_mul(2)
            .map(|k| k.max(1))
            .ok_or_else(|| format!("powi n {}", n)),
        other => Err(format!("unknown op {}", other)),
    }
}

/// `rad <= k·rref`, exactly. An infinite radius or one that does not convert is never within a bound.
fn within(rad: &Mag, k: i64, rref: &ExactNum) -> Result<bool, String> {
    if rad.is_infinity() {
        return Ok(false);
    }
    let r = rad
        .to_exact_num()
        .map_err(|e| format!("rad to exact: {:?}", e))?;
    match r.partial_cmp(&rref.mul_full_prec(&num(k))) {
        Some(Ordering::Less | Ordering::Equal) => Ok(true),
        Some(Ordering::Greater) => Ok(false),
        None => Err("rad vs bound unordered".into()),
    }
}

fn run(op: &str, prec: usize, k: i64, a: Ball, b: &Ball, bmid: &ExactNum) -> Result<Ball, String> {
    let n = || u64::try_from(k).map_err(|_| format!("{}: N {} is negative", op, k));
    match op {
        "add" => ok("add", a.add(b, prec)),
        "sub" => ok("sub", a.sub(b, prec)),
        "mul" => ok("mul", a.mul(b, prec)),
        "mul_exact" => ok("mul_exact", a.mul_exact(bmid, prec)),
        "sqr" => ok("sqr", a.sqr(prec)),
        "powi" => ok("powi", a.powi(n()?, prec)),
        "neg" => Ok(-a),
        "sqrt" => ok("sqrt", a.sqrt(prec, &DomainPolicy::Clip)),
        "div" => ok("div", a.div(b, prec)),
        "inv" => ok("inv", a.inv(prec)),
        "hypot" => ok("hypot", a.hypot(b, prec)),
        "abs" => ok("abs", a.abs(prec)),
        "mul_2exp" | "mul_2exp_ovf" | "mul_2exp_unf" | "mul_2exp_unf_rad" => {
            ok("mul_2exp", a.mul_2exp(k, prec))
        }
        other => Err(format!("unknown op {}", other)),
    }
}

fn check_row(line: &str) -> Result<(), String> {
    let v: Vec<&str> = line.split_whitespace().collect();
    if v.len() != 25 {
        return Err(format!("{} fields", v.len()));
    }
    let op = field(&v, 0)?;
    let prec: usize = field(&v, 1)?.parse().map_err(|_| "prec")?;
    let n: i64 = field(&v, 3)?.parse().map_err(|_| "n")?;
    let amid = value(&v, 4)?;
    let arad = Mag::from_exact_num(&value(&v, 7)?);
    let bmid = value(&v, 10)?;
    let brad = Mag::from_exact_num(&value(&v, 13)?);
    let lo = value(&v, 16)?;
    let hi = value(&v, 19)?;
    let rref = value(&v, 22)?;
    let a = ok(
        "a",
        Ball::from_mid_rad(amid.add_full_prec(&num(0)), arad, prec),
    )?;
    let b = ok(
        "b",
        Ball::from_mid_rad(bmid.add_full_prec(&num(0)), brad, prec),
    )?;
    // Strict must agree with Clip wherever the input ball lies in [0, inf).
    let strict_twin = if op == "sqrt" {
        let rn = a
            .rad()
            .to_exact_num()
            .map_err(|e| format!("rad: {:?}", e))?;
        match amid.sub_full_prec(&rn).partial_cmp(&num(0)) {
            Some(Ordering::Less) => None,
            Some(Ordering::Equal | Ordering::Greater) => Some(ok(
                "a2",
                Ball::from_mid_rad(amid.add_full_prec(&num(0)), Mag::from_exact_num(&rn), prec),
            )?),
            None => return Err("sqrt row: cannot order m - r against 0".into()),
        }
    } else {
        None
    };
    let sign_neg = amid.partial_cmp(&num(0)) == Some(Ordering::Less);
    let res = run(op, prec, n, a, &b, &bmid)?;
    // Out-of-range mul_2exp outcomes (decided by mpmath on the exact exponent): overflow is all reals; underflow
    // must contain [0, 2^(EXPONENT_MIN - 1)] on the input's side (the exact image lies strictly inside it).
    if op == "mul_2exp_ovf" {
        return if res.is_informative() {
            Err(format!("overflow row gave an informative ball {:?}", res))
        } else {
            Ok(())
        };
    }
    if op == "mul_2exp_unf" || op == "mul_2exp_unf_rad" {
        if op == "mul_2exp_unf_rad" {
            // A positive radius shifted below the Mag range clamps up to min_pos = finite(HALF_MANT, i64::MIN), never 0.
            let min_pos = Mag::finite(0x8000_0000_0000_0000, i64::MIN)
                .map_err(|e| format!("min_pos: {:?}", e))?;
            if res.rad() != &min_pos {
                return Err(format!(
                    "underflow rad {:?}, want min_pos {:?}",
                    res.rad(),
                    min_pos
                ));
            }
            // [LO, HI] of this tag is the mpmath truth at scale 2^-N (its true exponent is below every number range);
            // a mid-0 result contains it iff rad >= |LO|, |HI| scaled by 2^N.
            for (end, at) in [("LO", 16usize), ("HI", 19usize)] {
                if !rad_covers_scaled(res.rad(), &v, at, n)? {
                    return Err(format!(
                        "underflow rad {:?} does not cover truth {} * 2^{}",
                        res.rad(),
                        end,
                        n
                    ));
                }
            }
            return if res.mid() == &num(0) {
                Ok(())
            } else {
                Err(format!("mid 0 shifted to {:?}", res.mid()))
            };
        }
        let bmax = pow2(EXPONENT_MIN.checked_sub(1).ok_or("EXPONENT_MIN - 1")?)?;
        want_contains("underflow 0", &res, &num(0))?;
        return want_contains(
            "underflow bound",
            &res,
            &if sign_neg { -bmax } else { bmax },
        );
    }
    want_contains("LO", &res, &lo)?;
    want_contains("HI", &res, &hi)?;
    let t = tau(op, n)?;
    if !within(res.rad(), t.checked_mul(LAMBDA).ok_or("loose k")?, &rref)? {
        return Err(format!(
            "LOOSE: rad {:?} > 2^10·{}·RREF {:?}",
            res.rad(),
            t,
            rref
        ));
    }
    if !within(res.rad(), t, &rref)? {
        return Err(format!(
            "TIGHT: rad {:?} > {}·RREF {:?}",
            res.rad(),
            t,
            rref
        ));
    }
    if let Some(a2) = strict_twin {
        let s = ok("sqrt Strict", a2.sqrt(prec, &DomainPolicy::Strict))?;
        if s.mid() != res.mid() || s.rad() != res.rad() {
            return Err(format!(
                "Strict {:?} != Clip {:?} on a nonnegative ball",
                s, res
            ));
        }
    }
    Ok(())
}

#[test]
fn oracle_rows_contained_tight_loose() -> Result<(), String> {
    let path = golden_dir()
        .join("slices")
        .join("arith")
        .join("arith_oracle_rows.txt");
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{:?}: {}", path, e))?;
    let mut rows = 0usize;
    let mut per_op = WANT_PER_OP.map(|(op, _)| (op, 0usize));
    let mut bad: Vec<String> = Vec::new();
    for line in text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
    {
        rows = rows.checked_add(1).ok_or("rows")?;
        let tag = line.split_whitespace().next().ok_or("empty row")?;
        let slot = per_op
            .iter_mut()
            .find(|(op, _)| *op == tag)
            .ok_or_else(|| format!("row {}: unknown OP {}", rows, tag))?;
        slot.1 = slot.1.checked_add(1).ok_or("op count")?;
        if let Err(e) = check_row(line) {
            bad.push(format!(
                "row {}: {}\n    row: {}",
                rows,
                e,
                line.get(..160).unwrap_or(line)
            ));
        }
    }
    if rows != WANT_ROWS {
        return Err(format!("{} rows, want {}", rows, WANT_ROWS));
    }
    for ((op, got), (_, want)) in per_op.iter().zip(WANT_PER_OP.iter()) {
        if got != want {
            return Err(format!("{} {} rows, want {}", got, op, want));
        }
    }
    if !bad.is_empty() {
        let shown: Vec<&String> = bad.iter().take(12).collect();
        return Err(format!(
            "{} of {} oracle rows fail:\n{:#?}",
            bad.len(),
            rows,
            shown
        ));
    }
    Ok(())
}

/// The two hypot rows (data rows 1313 and 1343, 1-based, comments skipped) whose radius exceeds 4·RREF under the
/// 1-Lipschitz `r + s` radius: the second argument's ball straddles 0 widely, so the true spread is far below `s`.
/// Spence's monotone enclosure (hypot(min|X|, min|Y|) .. hypot(max|X|, max|Y|)) is to make them pass; RED until then.
const HYPOT_WIDE_ROWS: [usize; 2] = [1313, 1343];

#[test]
fn hypot_wide_rows_1313_1343_tight() -> Result<(), String> {
    let path = golden_dir()
        .join("slices")
        .join("arith")
        .join("arith_oracle_rows.txt");
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{:?}: {}", path, e))?;
    let data: Vec<&str> = text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .collect();
    let mut bad: Vec<String> = Vec::new();
    for id in HYPOT_WIDE_ROWS {
        let line = data
            .get(id.checked_sub(1).ok_or("row id 0")?)
            .ok_or_else(|| format!("row {} missing", id))?;
        if !line.starts_with("hypot ") {
            return Err(format!("row {} is not a hypot row: {}", id, line));
        }
        if let Err(e) = check_row(line) {
            bad.push(format!("row {}: {}", id, e));
        }
    }
    if bad.is_empty() {
        Ok(())
    } else {
        Err(format!("{:#?}", bad))
    }
}

/// Fixed powi section (Spence): data rows WANT_RANDOM+1 .. +36 are (1 + (37j+5)·2^-20)^n, n = 100 + 71j + i, at precision
/// [53, 64, 113][i], j in 0..12 (the loop of powi_exact_input_within_one_ulp). The exact power is LO = HI (mpmath.iv checked
/// at W + 64); RREF = U = one ulp of the power at PREC. Each row must pass check_row, and an informative result must have
/// rad <= RREF (one ulp); at least 30 of the 36 are informative, as in the loop.
#[test]
fn powi_fixed_rows_within_one_ulp() -> Result<(), String> {
    let path = golden_dir()
        .join("slices")
        .join("arith")
        .join("arith_oracle_rows.txt");
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{:?}: {}", path, e))?;
    let rows: Vec<&str> = text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .skip(WANT_RANDOM)
        .take(WANT_POWI_FIXED)
        .collect();
    if rows.len() != WANT_POWI_FIXED {
        return Err(format!(
            "{} fixed powi rows, want {}",
            rows.len(),
            WANT_POWI_FIXED
        ));
    }
    let mut informative = 0usize;
    for (at, (line, (i, j))) in rows
        .iter()
        .zip(
            [53usize, 64, 113]
                .into_iter()
                .enumerate()
                .flat_map(|(i, p)| (0..12i64).map(move |j| (i, (p, j)))),
        )
        .enumerate()
    {
        let (prec, jj) = j;
        let ii = i64::try_from(i).map_err(|_| "i")?;
        let n = jj
            .checked_mul(71)
            .and_then(|v| v.checked_add(100))
            .and_then(|v| v.checked_add(ii))
            .ok_or("n")?;
        let v: Vec<&str> = line.split_whitespace().collect();
        let want_head = format!("powi {}", prec);
        if !line.starts_with(&want_head) || field(&v, 3)? != n.to_string() {
            return Err(format!(
                "fixed powi row {}: want {} n {}: {}",
                at, want_head, n, line
            ));
        }
        check_row(line)
            .map_err(|e| format!("fixed powi row {} (p {} n {}): {}", at, prec, n, e))?;
        let a = ok(
            "a",
            Ball::from_mid_rad(value(&v, 4)?, Mag::from_exact_num(&value(&v, 7)?), prec),
        )?;
        let res = ok("powi", a.powi(u64::try_from(n).map_err(|_| "n")?, prec))?;
        if !res.is_informative() {
            continue;
        }
        if !within(res.rad(), 1, &value(&v, 22)?)? {
            return Err(format!(
                "fixed powi p {} n {}: rad {:?} > one ulp RREF",
                prec,
                n,
                res.rad()
            ));
        }
        informative = informative.checked_add(1).ok_or("c")?;
    }
    if informative < 30 {
        return Err(format!("only {} informative fixed powi rows", informative));
    }
    Ok(())
}
