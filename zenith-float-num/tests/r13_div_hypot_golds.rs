//! R13 (zenith-bound S1): div, reciprocal, hypot in every directed mode, and div_exact.
//!
//! `tests/oracle/goldens/r13_div_hypot.jsonl` comes from `scripts/gen_r13_goldens.py`:
//! mpmath at p+256 and p+512 bits (the two must agree, or both be the same exact value)
//! and an exact integer model; every row is limb-identical at exactly `p` bits over the
//! i64 exponent range, subnormal grid included. `div_exact` rows expect the quotient or
//! `none`.

#[path = "support/gold_limbs.rs"]
mod gold_limbs;

use serde_json::Value;
use zenith_float_num::{Error, ExactNum, Exponent, RoundingMode, Sign, Word};

const GOLDS: &str = include_str!("oracle/goldens/r13_div_hypot.jsonl");
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

fn canon_of(n: &ExactNum) -> String {
    if n.is_nan() {
        return format!("nan:{:?}", n.err());
    }
    if n.is_inf() {
        return format!("inf:{}", sign_str(n));
    }
    let limbs = gold_limbs::hex_of(n.mantissa_digits().unwrap_or(&[]));
    format!("f:{}:{}:{}", sign_str(n), n.exponent().unwrap_or(0), limbs)
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
        "none" => Ok("none".to_string()),
        "inf" => Ok(format!("inf:{}", text(v, "s")?)),
        "f" => {
            let limbs = gold_limbs::hex_of(&limbs_of(v)?);
            let e: Exponent = field(v, "e")?.as_i64().ok_or("e is not an i64")?;
            Ok(format!("f:{}:{}:{}", text(v, "s")?, e, limbs))
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

/// The call under test, in its canonical text.
fn eval(func: &str, args: &[ExactNum], p: usize, rm: &str) -> Result<String, String> {
    let x = args.first().ok_or("no argument")?;
    let y = args.get(1);
    let two = || y.ok_or_else(|| format!("{func} needs two arguments"));
    if func == "div_exact" {
        return match x.div_exact(two()?, p) {
            Ok(Some(q)) if q.inexact() => {
                Err(format!("div_exact returned an inexact {}", canon_of(&q)))
            }
            Ok(Some(q)) => Ok(canon_of(&q)),
            Ok(None) => Ok("none".to_string()),
            Err(e) => Err(format!("div_exact error {e:?}")),
        };
    }
    let rm = mode_of(rm)?;
    let r = match func {
        "div" => words_of(x).div(two()?, p, rm),
        "reciprocal" => x.reciprocal(p, rm),
        "hypot" => x.hypot(two()?, p, rm),
        other => return Err(format!("unknown function {other}")),
    };
    Ok(canon_of(&r))
}

/// `ExactNum::div` takes `self` by value: a new number built from the argument's words.
fn words_of(x: &ExactNum) -> ExactNum {
    match (x.mantissa_digits(), x.sign(), x.exponent()) {
        (Some(m), Some(s), Some(e)) => ExactNum::from_words(m, s, e),
        _ => ExactNum::nan(Some(Error::InvalidArgument)),
    }
}

fn run_row(row: &Value) -> Result<(), String> {
    let func = text(row, "fn")?;
    let p = usize::try_from(field(row, "p")?.as_u64().ok_or("p is not a u64")?)
        .map_err(|e| format!("p: {e}"))?;
    let raw = field(row, "args")?
        .as_array()
        .ok_or("args is not an array")?;
    let mut args = Vec::with_capacity(raw.len());
    for a in raw {
        args.push(arg_of(a)?);
    }
    let want = canon_expect(field(row, "expect")?)?;
    let got = eval(func, &args, p, text(row, "rm")?)?;
    if got == want {
        Ok(())
    } else {
        Err(format!("want {want} got {got}"))
    }
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
            fails.push(format!("{id}: {}", why.get(..200).unwrap_or(&why)));
        }
    }
    eprintln!(
        "R13 {func}: rows={total} pass={} fail={}",
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
fn r13_div() -> Result<(), String> {
    check("div")
}

#[test]
fn r13_reciprocal() -> Result<(), String> {
    check("reciprocal")
}

#[test]
fn r13_hypot() -> Result<(), String> {
    check("hypot")
}

#[test]
fn r13_div_exact() -> Result<(), String> {
    check("div_exact")
}

/// The golds reach both ends of the exponent range `±(2^62 − 1)` for div, reciprocal and hypot
/// (in an argument or in the expected result).
#[test]
fn golds_reach_exponent_extremes() -> Result<(), String> {
    use zenith_float_num::{EXPONENT_MAX, EXPONENT_MIN};
    assert_eq!(EXPONENT_MAX, (1i64 << 62) - 1);
    for func in ["div", "reciprocal", "hypot"] {
        for extreme in [EXPONENT_MAX, EXPONENT_MIN] {
            let mut hits = 0usize;
            for line in GOLDS.lines().filter(|l| !l.trim().is_empty()) {
                let row: Value = serde_json::from_str(line).map_err(|e| e.to_string())?;
                if text(&row, "fn")? != func {
                    continue;
                }
                let args = field(&row, "args")?
                    .as_array()
                    .ok_or("args is not an array")?;
                let hit = args
                    .iter()
                    .chain([field(&row, "expect")?])
                    .any(|v| v.get("e").and_then(Value::as_i64) == Some(extreme));
                hits += usize::from(hit);
            }
            if hits == 0 {
                return Err(format!("no {func} gold at exponent {extreme}"));
            }
        }
    }
    Ok(())
}
