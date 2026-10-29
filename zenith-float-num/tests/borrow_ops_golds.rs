//! The borrowing operations `ExactNum::{add_ref, sub_ref, mul_ref, div_ref}` (both operands borrowed,
//! only the result allocated, rounded once at exactly `p` bits).
//!
//! `tests/oracle/goldens/borrow_ops.jsonl` comes from `scripts/gen_borrow_goldens.py`: mpmath at p+256
//! and p+512 bits (the two must agree, or be the same exact value) and an exact integer model. Rows cover
//! every rounding mode, both signs, signed zeros, subnormal operands and results, operands next to
//! `EXPONENT_MAX` (directional overflow) and aliasing (`op(&x, &x)`). Every row is also checked
//! bit-identical (limbs, exponent, sign, `inexact`, storage length) to the by-value operation whose first
//! operand is moved in.

#[path = "support/gold_limbs.rs"]
mod gold_limbs;

use serde_json::Value;
use zenith_float_num::{Error, ExactNum, Exponent, RoundingMode, Sign, Word, WORD_BIT_SIZE};

const GOLDS: &str = include_str!("oracle/goldens/borrow_ops.jsonl");
const SHOW: usize = 8;

fn mode_of(name: &str) -> Result<RoundingMode, String> {
    match name {
        "ToEven" => Ok(RoundingMode::ToEven),
        "Up" => Ok(RoundingMode::Up),
        "Down" => Ok(RoundingMode::Down),
        "ToZero" => Ok(RoundingMode::ToZero),
        "FromZero" => Ok(RoundingMode::FromZero),
        other => Err(format!("unknown mode {other}")),
    }
}

fn sign_str(n: &ExactNum) -> &'static str {
    if n.sign() == Some(Sign::Neg) {
        "-"
    } else {
        "+"
    }
}

/// Canonical text of a result; a zero is `zero:<sign>` (its stored exponent is not part of the value).
fn canon_of(n: &ExactNum) -> String {
    if n.is_nan() {
        return format!("nan:{:?}", n.err());
    }
    if n.is_inf() {
        return format!("inf:{}", sign_str(n));
    }
    if n.is_zero() {
        return format!("zero:{}", sign_str(n));
    }
    let limbs = gold_limbs::hex_of(n.mantissa_digits().unwrap_or(&[]));
    format!("f:{}:{}:{}", sign_str(n), n.exponent().unwrap_or(0), limbs)
}

/// Everything observable about a result, for the bit-identity check against the by-value operation.
fn full_of(n: &ExactNum) -> String {
    format!(
        "{} e={:?} m={:?} inexact={} len={:?}",
        canon_of(n),
        n.exponent(),
        n.mantissa_digits(),
        n.inexact(),
        n.mantissa_max_bit_len()
    )
}

fn field<'a>(v: &'a Value, key: &str) -> Result<&'a Value, String> {
    v.get(key).ok_or_else(|| format!("missing {key}"))
}

fn text<'a>(v: &'a Value, key: &str) -> Result<&'a str, String> {
    field(v, key)?
        .as_str()
        .ok_or_else(|| format!("{key} is not a string"))
}

fn limbs_of(v: &Value) -> Result<Vec<Word>, String> {
    let arr = field(v, "m")?.as_array().ok_or("m is not an array")?;
    let mut out = Vec::with_capacity(arr.len());
    for w in arr {
        let s = w.as_str().ok_or("limb is not a string")?;
        out.extend(gold_limbs::words_of(s)?);
    }
    Ok(out)
}

fn canon_expect(v: &Value) -> Result<String, String> {
    match text(v, "k")? {
        "nan" => Ok(format!("nan:Some({})", text(v, "err")?)),
        "inf" => Ok(format!("inf:{}", text(v, "s")?)),
        "f" => {
            let limbs = limbs_of(v)?;
            if limbs.iter().all(|w| *w == 0) {
                return Ok(format!("zero:{}", text(v, "s")?));
            }
            let shown = gold_limbs::hex_of(&limbs);
            let e: Exponent = field(v, "e")?.as_i64().ok_or("e is not an i64")?;
            Ok(format!("f:{}:{}:{}", text(v, "s")?, e, shown))
        }
        other => Err(format!("unknown kind {other}")),
    }
}

fn arg_of(v: &Value) -> Result<ExactNum, String> {
    let limbs = limbs_of(v)?;
    let sign = match text(v, "s")? {
        "-" => Sign::Neg,
        "+" => Sign::Pos,
        other => return Err(format!("unknown sign {other}")),
    };
    let e: Exponent = field(v, "e")?.as_i64().ok_or("e is not an i64")?;
    let x = ExactNum::from_words(&limbs, sign, e);
    if x.is_nan() {
        return Err(format!("argument out of range: e={e}"));
    }
    Ok(x)
}

type RefOp = fn(&ExactNum, &ExactNum, usize, RoundingMode) -> ExactNum;
type ValOp = fn(ExactNum, &ExactNum, usize, RoundingMode) -> ExactNum;

fn ops(func: &str) -> Result<(RefOp, ValOp), String> {
    match func {
        "add" => Ok((ExactNum::add_ref, ExactNum::add)),
        "sub" => Ok((ExactNum::sub_ref, ExactNum::sub)),
        "mul" => Ok((ExactNum::mul_ref, ExactNum::mul)),
        "div" => Ok((ExactNum::div_ref, ExactNum::div)),
        other => Err(format!("unknown function {other}")),
    }
}

/// A fresh number with `x`'s value, to be moved into the by-value operation.
fn words_of(x: &ExactNum) -> ExactNum {
    match (x.mantissa_digits(), x.sign(), x.exponent()) {
        (Some(m), Some(s), Some(e)) => ExactNum::from_words(m, s, e),
        _ => ExactNum::nan(Some(Error::InvalidArgument)),
    }
}

fn run_row(row: &Value) -> Result<(), String> {
    let (by_ref, by_val) = ops(text(row, "fn")?)?;
    let p = usize::try_from(field(row, "p")?.as_u64().ok_or("p is not a u64")?)
        .map_err(|e| format!("p: {e}"))?;
    let rm = mode_of(text(row, "rm")?)?;
    let alias = field(row, "alias")?
        .as_bool()
        .ok_or("alias is not a bool")?;
    let raw = field(row, "args")?
        .as_array()
        .ok_or("args is not an array")?;
    let x = arg_of(raw.first().ok_or("no x")?)?;
    let y = arg_of(raw.get(1).ok_or("no y")?)?;
    let want = canon_expect(field(row, "expect")?)?;
    let before = (full_of(&x), full_of(&y));
    let got = if alias { by_ref(&x, &x, p, rm) } else { by_ref(&x, &y, p, rm) };
    if (full_of(&x), full_of(&y)) != before {
        return Err("an operand changed".to_string());
    }
    if canon_of(&got) != want {
        return Err(format!("want {want} got {}", canon_of(&got)));
    }
    let moved = by_val(words_of(&x), if alias { &x } else { &y }, p, rm);
    if full_of(&got) != full_of(&moved) {
        return Err(format!(
            "borrowed {} differs from by-value {}",
            full_of(&got),
            full_of(&moved)
        ));
    }
    // Only the output is allocated, at the output's size: a finite result is stored in p rounded up to whole
    // words (the size of the by-value result, checked above), never the operands' or a working width.
    if let Some(len) = got.mantissa_max_bit_len() {
        let words = p.div_ceil(WORD_BIT_SIZE);
        if len != words * WORD_BIT_SIZE {
            return Err(format!("result storage {len} bits for p={p}"));
        }
    }
    Ok(())
}

fn check(func: &str) -> Result<(), String> {
    let mut total = 0usize;
    let mut fails: Vec<String> = Vec::new();
    for line in GOLDS.lines().filter(|l| !l.trim().is_empty()) {
        let row: Value = serde_json::from_str(line).map_err(|e| format!("json: {e}"))?;
        if row.get("fn").and_then(Value::as_str) != Some(func) {
            continue;
        }
        total += 1;
        if let Err(why) = run_row(&row) {
            let id = row.get("id").and_then(Value::as_str).unwrap_or("?");
            fails.push(format!("{id}: {}", why.get(..300).unwrap_or(&why)));
        }
    }
    eprintln!(
        "borrow {func}: rows={total} pass={} fail={}",
        total - fails.len(),
        fails.len()
    );
    for f in fails.iter().take(SHOW) {
        eprintln!("  {f}");
    }
    if total == 0 {
        return Err(format!("{func}: no golden rows"));
    }
    if fails.is_empty() {
        Ok(())
    } else {
        Err(format!("{func}: {}/{total} rows differ", fails.len()))
    }
}

#[test]
fn borrow_add_ref() -> Result<(), String> {
    check("add")
}

#[test]
fn borrow_sub_ref() -> Result<(), String> {
    check("sub")
}

#[test]
fn borrow_mul_ref() -> Result<(), String> {
    check("mul")
}

#[test]
fn borrow_div_ref() -> Result<(), String> {
    check("div")
}

/// The golds cover every mode, both signs, aliasing, and both ends of the exponent range.
#[test]
fn golds_cover_the_edges() -> Result<(), String> {
    use zenith_float_num::{EXPONENT_MAX, EXPONENT_MIN};
    for func in ["add", "sub", "mul", "div"] {
        let (mut alias, mut top, mut bottom, mut neg, mut modes) = (0, 0, 0, 0, Vec::new());
        for line in GOLDS.lines().filter(|l| !l.trim().is_empty()) {
            let row: Value = serde_json::from_str(line).map_err(|e| e.to_string())?;
            if text(&row, "fn")? != func {
                continue;
            }
            alias += usize::from(field(&row, "alias")?.as_bool() == Some(true));
            let args = field(&row, "args")?.as_array().ok_or("args")?;
            let es: Vec<Option<i64>> = args
                .iter()
                .chain([field(&row, "expect")?])
                .map(|v| v.get("e").and_then(Value::as_i64))
                .collect();
            top += usize::from(es.contains(&Some(EXPONENT_MAX)));
            bottom += usize::from(es.contains(&Some(EXPONENT_MIN)));
            neg += usize::from(
                args.iter()
                    .any(|a| a.get("s").and_then(Value::as_str) == Some("-")),
            );
            let rm = text(&row, "rm")?.to_string();
            if !modes.contains(&rm) {
                modes.push(rm);
            }
        }
        if alias == 0 || top == 0 || bottom == 0 || neg == 0 || modes.len() != 5 {
            return Err(format!(
                "{func}: alias={alias} top={top} bottom={bottom} neg={neg} modes={modes:?}"
            ));
        }
    }
    Ok(())
}
