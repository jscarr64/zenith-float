//! R15b: `exp`, `exp2` and `pow` (bases 1.5 and 4) at p in {8, 24, 32} of integer-valued
//! arguments with exponent 33..61. The argument is held at exactly `p` bits (one 32-bit word on
//! a 32-bit target), so its integer part is wider than its mantissa there. Golds from
//! `scripts/gen_r15b_int_goldens.py` (exact powers of two, or mpmath at p+256 and p+1024 bits),
//! compared by value on any word size, inexact flag included.

use serde_json::Value;
use zenith_float_num::{Consts, ExactNum, Exponent, RoundingMode, Sign, Word};

const GOLDS: &str = include_str!("oracle/goldens/r15b_int_exp.jsonl");
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
    format!(
        "f:{}:{}:{}",
        sign_str(n),
        n.exponent().unwrap_or(0),
        hex_of(n.mantissa_digits().unwrap_or(&[]))
    )
}

/// The mantissa as hex, most significant word first, trailing zeros dropped: the same text for
/// 64- and 32-bit words (the value, not the layout).
fn hex_of(m: &[Word]) -> String {
    let width = usize::try_from(Word::BITS / 4).unwrap_or(16);
    let s: String = m.iter().rev().map(|w| format!("{w:0width$x}")).collect();
    s.trim_end_matches('0').to_string()
}

fn field<'a>(v: &'a Value, key: &str) -> Result<&'a Value, String> {
    v.get(key).ok_or_else(|| format!("missing {key}"))
}

fn text<'a>(v: &'a Value, key: &str) -> Result<&'a str, String> {
    field(v, key)?
        .as_str()
        .ok_or_else(|| format!("{key} is not a string"))
}

/// The gold limbs (64-bit, little-endian) as this target's words.
fn limbs_of(v: &Value) -> Result<Vec<Word>, String> {
    let arr = field(v, "m")?.as_array().ok_or("m is not an array")?;
    let mut out = Vec::with_capacity(2 * arr.len());
    for w in arr {
        let s = w.as_str().ok_or("limb is not a string")?;
        let x = u64::from_str_radix(s, 16).map_err(|e| format!("limb {s}: {e}"))?;
        let halves = if Word::BITS == 64 { vec![x] } else { vec![x & 0xffff_ffff, x >> 32] };
        for h in halves {
            out.push(Word::try_from(h).map_err(|e| format!("limb {s}: {e}"))?);
        }
    }
    Ok(out)
}

fn canon_expect(v: &Value) -> Result<String, String> {
    match text(v, "k")? {
        "inf" => Ok(format!("inf:{}", text(v, "s")?)),
        "f" => {
            let limbs = hex_of(&limbs_of(v)?);
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

fn run_row(row: &Value, cc: &mut Consts) -> Result<(), String> {
    let func = text(row, "fn")?;
    let p = usize::try_from(field(row, "p")?.as_u64().ok_or("p is not a u64")?)
        .map_err(|e| format!("p: {e}"))?;
    let rm = mode_of(text(row, "rm")?)?;
    let mut args = Vec::new();
    for a in field(row, "args")?
        .as_array()
        .ok_or("args is not an array")?
    {
        let mut x = arg_of(a)?;
        // Exact: every argument fits p bits. On a 32-bit target this leaves one word.
        x.set_precision(p, RoundingMode::ToZero)
            .map_err(|e| format!("set_precision: {e:?}"))?;
        args.push(x);
    }
    let r = match (func, args.as_slice()) {
        ("exp", [x]) => x.exp(p, rm, cc),
        ("exp2", [x]) => x.exp2(p, rm, cc),
        ("pow15" | "pow4", [b, y]) => b.pow(y, p, rm, cc),
        _ => return Err(format!("bad row {func}/{}", args.len())),
    };
    let want = canon_expect(field(row, "expect")?)?;
    let got = canon_of(&r);
    if got != want {
        return Err(format!("want {want} got {got}"));
    }
    let flag = field(row, "inexact")?
        .as_bool()
        .ok_or("inexact is not a bool")?
        && !r.is_inf();
    if r.inexact() != flag {
        return Err(format!("inexact flag {} want {flag}", r.inexact()));
    }
    Ok(())
}

fn check(func: &str) -> Result<(), String> {
    let mut cc = Consts::new().map_err(|e| format!("consts: {e:?}"))?;
    let mut total = 0usize;
    let mut fails: Vec<String> = Vec::new();
    for line in GOLDS.lines().filter(|l| !l.trim().is_empty()) {
        let row: Value = serde_json::from_str(line).map_err(|e| format!("json: {e}"))?;
        if row.get("fn").and_then(Value::as_str) != Some(func) {
            continue;
        }
        total += 1;
        if let Err(why) = run_row(&row, &mut cc) {
            let id = row.get("id").and_then(Value::as_str).unwrap_or("?");
            fails.push(format!("{id}: {}", why.get(..200).unwrap_or(&why)));
        }
    }
    eprintln!(
        "R15b {func}: rows={total} pass={} fail={}",
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
fn r15b_exp_int_past_32_bits() -> Result<(), String> {
    check("exp")
}

#[test]
fn r15b_exp2_int_past_32_bits() -> Result<(), String> {
    check("exp2")
}

#[test]
fn r15b_pow_int_past_32_bits() -> Result<(), String> {
    check("pow15")?;
    check("pow4")
}
