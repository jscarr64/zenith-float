//! R9: overflow and underflow at the i64 exponent limits, every rounding mode.
//!
//! `tests/oracle/goldens/r9_overflow.jsonl` comes from
//! `scripts/gen_r9_overflow_goldens.py`: an exact integer model (add, sub, mul,
//! sqr, fma, ldexp, powi, powsi, exact exp2/pow, short decimal parse) checked
//! against mpmath at p+256 and p+512, and mpmath Ziv pairs for exp, exp2, exp10,
//! sinh, cosh, pow and huge-exponent parse. Rows are limb-identical. Overflow is
//! +-Inf or +-MAX by direction, underflow lands on the subnormal grid (or 0 / the
//! min subnormal by direction) with the sign kept. A zero's exponent is not
//! compared (only its sign).

use std::sync::mpsc;
use std::thread;
use std::time::Duration;

use serde_json::Value;
use zenith_float_num::{Consts, ExactNum, Exponent, Radix, RoundingMode, Sign, Word};

const GOLDS: &str = include_str!("oracle/goldens/r9_overflow.jsonl");
const ROW_BUDGET: Duration = Duration::from_secs(30);
const SHOW: usize = 400;

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
    if n.is_negative() {
        "-"
    } else {
        "+"
    }
}

fn canon_of(n: &ExactNum) -> String {
    if n.is_nan() {
        return "nan".to_string();
    }
    if n.is_inf() {
        return format!("inf:{}", sign_str(n));
    }
    let digits = n.mantissa_digits().unwrap_or(&[]);
    let is_zero = digits.iter().all(|w| *w == 0);
    let limbs: Vec<String> = digits.iter().map(|w| format!("{w:016x}")).collect();
    let e = if is_zero { 0 } else { n.exponent().unwrap_or(0) };
    format!("f:{}:{}:{}", sign_str(n), e, limbs.join(","))
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
        out.push(Word::from_str_radix(s, 16).map_err(|e| format!("limb {s}: {e}"))?);
    }
    Ok(out)
}

fn canon_expect(v: &Value) -> Result<String, String> {
    match text(v, "k")? {
        "nan" => Ok("nan".to_string()),
        "inf" => Ok(format!("inf:{}", text(v, "s")?)),
        "f" => {
            let limbs: Vec<String> = limbs_of(v)?
                .iter()
                .map(|w| format!("{w:016x}"))
                .collect();
            let e: Exponent = field(v, "e")?.as_i64().ok_or("e is not an i64")?;
            Ok(format!("f:{}:{}:{}", text(v, "s")?, e, limbs.join(",")))
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

fn arg_n(raw: &[Value], i: usize) -> Result<ExactNum, String> {
    arg_of(raw.get(i).ok_or_else(|| format!("no argument {i}"))?)
}

fn int_of(row: &Value) -> Result<i64, String> {
    field(row, "n")?.as_i64().ok_or_else(|| "n is not an i64".to_string())
}

fn eval(row: &Value) -> Result<ExactNum, String> {
    let func = text(row, "fn")?;
    let p = usize::try_from(field(row, "p")?.as_u64().ok_or("p is not a u64")?)
        .map_err(|e| format!("p: {e}"))?;
    let rm = mode_of(text(row, "rm")?)?;
    let raw = field(row, "args")?.as_array().ok_or("args is not an array")?;
    let mut cc = Consts::new().map_err(|e| format!("consts: {e:?}"))?;
    let r = match func {
        "add" => arg_n(raw, 0)?.add(&arg_n(raw, 1)?, p, rm),
        "sub" => arg_n(raw, 0)?.sub(&arg_n(raw, 1)?, p, rm),
        "mul" => arg_n(raw, 0)?.mul(&arg_n(raw, 1)?, p, rm),
        "sqr" => arg_n(raw, 0)?.sqr(p, rm),
        "fma" => arg_n(raw, 0)?.fma(&arg_n(raw, 1)?, &arg_n(raw, 2)?, p, rm),
        "ldexp" => arg_n(raw, 0)?.ldexp(int_of(row)?, p, rm),
        "powi" => {
            let n = usize::try_from(int_of(row)?).map_err(|e| format!("n: {e}"))?;
            arg_n(raw, 0)?.powi(n, p, rm)
        }
        "powsi" => arg_n(raw, 0)?.powsi(int_of(row)?, p, rm),
        "pow" => arg_n(raw, 0)?.pow(&arg_n(raw, 1)?, p, rm, &mut cc),
        "exp" => arg_n(raw, 0)?.exp(p, rm, &mut cc),
        "exp2" => arg_n(raw, 0)?.exp2(p, rm, &mut cc),
        "exp10" => arg_n(raw, 0)?.exp10(p, rm, &mut cc),
        "sinh" => arg_n(raw, 0)?.sinh(p, rm, &mut cc),
        "cosh" => arg_n(raw, 0)?.cosh(p, rm, &mut cc),
        "parse" => ExactNum::parse(text(row, "s")?, Radix::Dec, p, rm, &mut cc),
        other => return Err(format!("unknown function {other}")),
    };
    Ok(r)
}

fn run_line(line: &'static str) -> Result<(), String> {
    let row: Value = serde_json::from_str(line).map_err(|e| format!("json: {e}"))?;
    let want = canon_expect(field(&row, "expect")?)?;
    let got = canon_of(&eval(&row)?);
    if got == want {
        Ok(())
    } else {
        Err(format!("want {want} got {got}"))
    }
}

fn label(line: &str) -> String {
    let row: Value = serde_json::from_str(line).unwrap_or(Value::Null);
    let g = |k: &str| row.get(k).map(Value::to_string).unwrap_or_default();
    format!("{} p={} {} [{}]", g("fn"), g("p"), g("rm"), g("note"))
}

fn check(func: &str) -> Result<(), String> {
    let tag = format!("\"fn\": \"{func}\"");
    let mut rows = 0usize;
    let mut bad = Vec::new();
    for line in GOLDS.lines().filter(|l| l.contains(&tag)) {
        rows += 1;
        let (tx, rx) = mpsc::channel();
        thread::spawn(move || {
            let _ = tx.send(run_line(line));
        });
        let outcome = match rx.recv_timeout(ROW_BUDGET) {
            Ok(r) => r,
            Err(_) => Err("timed out".to_string()),
        };
        if let Err(why) = outcome {
            bad.push(format!("{}: {why}", label(line)));
        }
    }
    if rows == 0 {
        return Err(format!("no rows for {func}"));
    }
    eprintln!("r9 {func}: {rows} rows, {} failed", bad.len());
    if bad.is_empty() {
        return Ok(());
    }
    let shown: Vec<&String> = bad.iter().take(SHOW).collect();
    Err(format!(
        "{func}: {}/{rows} rows wrong:\n{}",
        bad.len(),
        shown.iter().map(|s| s.as_str()).collect::<Vec<_>>().join("\n")
    ))
}

macro_rules! golds {
    ($($name:ident => $f:literal),* $(,)?) => {
        $(#[test]
        fn $name() -> Result<(), String> {
            check($f)
        })*
    };
}

golds! {
    r9_add => "add", r9_sub => "sub", r9_mul => "mul", r9_sqr => "sqr", r9_fma => "fma",
    r9_ldexp => "ldexp", r9_powi => "powi", r9_powsi => "powsi", r9_pow => "pow",
    r9_exp => "exp", r9_exp2 => "exp2", r9_exp10 => "exp10", r9_sinh => "sinh",
    r9_cosh => "cosh", r9_parse => "parse",
}
