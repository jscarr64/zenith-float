//! S3: the Euler–Mascheroni constant γ from `Consts::euler_gamma(p, rm)` is correctly rounded
//! at exactly `p` bits in all five modes, like `Consts::pi` and `Consts::ln_2`, from a cold
//! cache and from a warm one, and agrees with `Consts::get_exact_prec(ConstantId::EulerGamma)`.
//!
//! `tests/oracle/goldens/s3_euler_gamma.jsonl` comes from `scripts/gen_s3_euler_goldens.py`:
//! mpmath `mpf_euler` at p+256 and p+512 bits, each error interval settled to one p-bit
//! result and the two agreeing, for p in {24, 53, 64, 113, 256, 1024, 4096, 65536, 131072}.

use serde_json::Value;
use std::time::Instant;
use zenith_float_num::{ConstantId, Consts, ExactNum, RoundingMode, Sign, Word};

const GOLDS: &str = include_str!("oracle/goldens/s3_euler_gamma.jsonl");

/// Rows at or above this precision run in `euler_gamma_golds_high`.
const HIGH_P: usize = 65536;

fn value_of(n: &ExactNum) -> String {
    if n.is_nan() || n.is_inf() {
        return format!("special:{:?}", n.err());
    }
    let limbs: Vec<String> = n
        .mantissa_digits()
        .unwrap_or(&[])
        .iter()
        .map(|w| format!("{w:016x}"))
        .collect();
    let s = if n.sign() == Some(Sign::Neg) { "-" } else { "+" };
    format!("f:{s}:{}:{}", n.exponent().unwrap_or(0), limbs.join(","))
}

fn field<'a>(v: &'a Value, key: &str) -> Result<&'a Value, String> {
    v.get(key).ok_or_else(|| format!("missing {key}"))
}

fn text<'a>(v: &'a Value, key: &str) -> Result<&'a str, String> {
    field(v, key)?
        .as_str()
        .ok_or_else(|| format!("{key} is not a string"))
}

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

fn expect_of(v: &Value) -> Result<ExactNum, String> {
    let arr = field(v, "m")?.as_array().ok_or("m is not an array")?;
    let mut limbs = Vec::with_capacity(arr.len());
    for w in arr {
        let s = w.as_str().ok_or("limb is not a string")?;
        limbs.push(Word::from_str_radix(s, 16).map_err(|e| format!("limb {s}: {e}"))?);
    }
    let e = field(v, "e")?.as_i64().ok_or("e is not an i64")?;
    let sign = if text(v, "s")? == "-" { Sign::Neg } else { Sign::Pos };
    Ok(ExactNum::from_words(&limbs, sign, e))
}

/// (id, p, mode, expected) for every gold row with `keep(p)`.
fn rows(keep: fn(usize) -> bool) -> Result<Vec<(String, usize, RoundingMode, ExactNum)>, String> {
    let mut out = Vec::new();
    for line in GOLDS.lines().filter(|l| !l.trim().is_empty()) {
        let row: Value = serde_json::from_str(line).map_err(|e| format!("json: {e}"))?;
        let p = field(&row, "p")?.as_u64().ok_or("p is not a u64")?;
        let p = usize::try_from(p).map_err(|e| e.to_string())?;
        if keep(p) {
            out.push((
                text(&row, "id")?.to_string(),
                p,
                mode_of(text(&row, "rm")?)?,
                expect_of(field(&row, "expect")?)?,
            ));
        }
    }
    Ok(out)
}

/// Every row from a cold cache, then again from the warm one; both getters must match the gold.
fn check(keep: fn(usize) -> bool) -> Result<(), String> {
    let rows = rows(keep)?;
    if rows.is_empty() {
        return Err("no gold rows".to_string());
    }
    let mut bad = Vec::new();
    let mut warm = Consts::new().map_err(|e| format!("{e:?}"))?;
    for (id, p, rm, want) in &rows {
        let mut cold = Consts::new().map_err(|e| format!("{e:?}"))?;
        let t = Instant::now();
        let got = cold.euler_gamma(*p, *rm);
        let ms = t.elapsed().as_millis();
        if *rm == RoundingMode::ToEven {
            println!("s3 timing: euler_gamma p={p} cold {ms} ms");
        }
        let warm_got = warm.euler_gamma(*p, *rm);
        let exact = cold
            .get_exact_prec(ConstantId::EulerGamma, *p, *rm)
            .map_err(|e| format!("{id}: get_exact_prec {e:?}"))?;
        let w = value_of(want);
        for (what, g) in [("cold", &got), ("warm", &warm_got), ("get_exact_prec", &exact)] {
            if value_of(g) != w || !g.inexact() {
                bad.push(format!(
                    "{id} {what}: got {} inexact={} want {w}",
                    value_of(g),
                    g.inexact()
                ));
            }
        }
    }
    if bad.is_empty() {
        return Ok(());
    }
    let shown: Vec<&str> = bad.iter().take(8).map(String::as_str).collect();
    Err(format!(
        "{} of {} checks failed:\n{}",
        bad.len(),
        3 * rows.len(),
        shown.join("\n")
    ))
}

#[test]
fn euler_gamma_golds_low() -> Result<(), String> {
    check(|p| p < HIGH_P)
}

#[test]
fn euler_gamma_golds_high() -> Result<(), String> {
    check(|p| p >= HIGH_P)
}
