//! S9 directed special functions (R7, R7b): red goldens.
//!
//! `tests/oracle/goldens/s9_specials.jsonl` comes from
//! `scripts/gen_s9_specials_goldens.py` (mpmath, Ziv pair at p+256 and p+512
//! bits). Every row is limb-identical: one correctly rounded result per
//! rounding mode (`ToEven`, `Up`, `Down`, `ToZero`, `FromZero`) at exactly `p`
//! bits over the i64 exponent range.
//!
//! One test per function. A test fails when any of its rows differs, times out,
//! or has no directed implementation yet (`Unsupported`). `beta` has no public
//! real API today, so it goes through a local stub that returns `Unsupported`.

use std::sync::mpsc;
use std::thread;
use std::time::{Duration, Instant};

use serde_json::Value;
use zenith_float_num::{Consts, Error, ExactNum, Exponent, RoundingMode, Sign, Word};

const GOLDS: &str = include_str!("oracle/goldens/s9_specials.jsonl");

/// Wall budget for one row, and for one function's whole table.
const ROW_BUDGET: Duration = Duration::from_secs(10);
const FN_BUDGET: Duration = Duration::from_secs(120);
/// Mismatch lines printed per function.
const SHOW: usize = 8;

/// The typed "no directed implementation yet" outcome.
#[derive(Debug)]
struct Unsupported;

/// Stub for the directed real beta `B(a, b)`; no public real API exists yet.
fn beta_directed(
    _a: &ExactNum,
    _b: &ExactNum,
    _p: usize,
    _rm: RoundingMode,
    _cc: &mut Consts,
) -> Result<ExactNum, Unsupported> {
    Err(Unsupported)
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

fn err_name(err: Option<Error>) -> String {
    match err {
        None => "null".to_string(),
        Some(Error::InvalidArgument) => "InvalidArgument".to_string(),
        Some(Error::DivisionByZero) => "DivisionByZero".to_string(),
        Some(Error::MemoryAllocation) => "MemoryAllocation".to_string(),
        Some(Error::PrecisionRetryExhausted) => "PrecisionRetryExhausted".to_string(),
        Some(Error::ExponentOverflow(_)) => "ExponentOverflow".to_string(),
    }
}

fn sign_str(n: &ExactNum) -> &'static str {
    if n.sign() == Some(Sign::Neg) {
        "-"
    } else {
        "+"
    }
}

/// Canonical text of a result, in the golden's shape.
fn canon_of(n: &ExactNum) -> String {
    if n.is_nan() {
        return format!("nan:{}", err_name(n.err()));
    }
    if n.is_inf() {
        return format!("inf:{}", sign_str(n));
    }
    let limbs: Vec<String> = n
        .mantissa_digits()
        .unwrap_or(&[])
        .iter()
        .map(|word| format!("{word:016x}"))
        .collect();
    format!(
        "f:{}:{}:{}",
        sign_str(n),
        n.exponent().unwrap_or(0),
        limbs.join(",")
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
        out.push(Word::from_str_radix(s, 16).map_err(|e| format!("limb {s}: {e}"))?);
    }
    Ok(out)
}

/// Canonical text of a golden `expect`.
fn canon_expect(v: &Value) -> Result<String, String> {
    match text(v, "k")? {
        "nan" => {
            let err = match field(v, "err")? {
                Value::String(s) => s.as_str(),
                Value::Null => "null",
                Value::Bool(_) | Value::Number(_) | Value::Array(_) | Value::Object(_) => {
                    return Err("err is not a string".to_string())
                }
            };
            Ok(format!("nan:{err}"))
        }
        "inf" => Ok(format!("inf:{}", text(v, "s")?)),
        "f" => {
            let limbs: Vec<String> = limbs_of(v)?
                .iter()
                .map(|word| format!("{word:016x}"))
                .collect();
            let e: Exponent = field(v, "e")?.as_i64().ok_or("e is not an i64")?;
            Ok(format!("f:{}:{}:{}", text(v, "s")?, e, limbs.join(",")))
        }
        other => Err(format!("unknown kind {other}")),
    }
}

/// The exact argument a golden encodes.
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

/// The directed call under test.
fn eval(
    func: &str,
    args: &[ExactNum],
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Result<ExactNum, String> {
    let x = args.first().ok_or("no argument")?;
    let r = match func {
        "gamma" => x.gamma(p, rm, cc),
        "ln_gamma" => x.ln_gamma(p, rm, cc),
        "erf" => x.erf(p, rm, cc),
        "erfc" => x.erfc(p, rm, cc),
        "expm1" => x.expm1(p, rm, cc),
        "log1p" => x.log1p(p, rm, cc),
        "digamma" => x.digamma(p, rm, cc),
        "bessel_j0" => x.bessel_j(0, p, rm, cc),
        "bessel_j1" => x.bessel_j(1, p, rm, cc),
        "zeta" => x.riemann_zeta(p, rm, cc),
        "beta" => {
            let y = args.get(1).ok_or("beta needs two arguments")?;
            beta_directed(x, y, p, rm, cc).map_err(|u| format!("{u:?}"))?
        }
        other => return Err(format!("unknown function {other}")),
    };
    Ok(r)
}

/// Outcome of one row: `Ok(())` when limb-identical, else the reason.
fn run_row(row: &Value, cc: &mut Consts) -> Result<(), String> {
    let func = text(row, "fn")?;
    let p = field(row, "p")?.as_u64().ok_or("p is not a u64")?;
    let p = usize::try_from(p).map_err(|e| format!("p: {e}"))?;
    let rm = mode_of(text(row, "rm")?)?;
    let raw = field(row, "args")?
        .as_array()
        .ok_or("args is not an array")?;
    let mut args = Vec::with_capacity(raw.len());
    for a in raw {
        args.push(arg_of(a)?);
    }
    let want = canon_expect(field(row, "expect")?)?;
    let got = canon_of(&eval(func, &args, p, rm, cc)?);
    if got == want {
        Ok(())
    } else {
        Err(format!("want {want} got {got}"))
    }
}

fn rows_for(func: &str) -> Result<Vec<Value>, String> {
    let mut out = Vec::new();
    for line in GOLDS.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let row: Value = serde_json::from_str(line).map_err(|e| format!("json: {e}"))?;
        if row.get("fn").and_then(Value::as_str) == Some(func) {
            out.push(row);
        }
    }
    Ok(out)
}

fn short(s: &str) -> &str {
    s.get(..160).unwrap_or(s)
}

/// Runs every row of `func` on a worker thread with per-row and per-table budgets.
fn check(func: &'static str) -> Result<(), String> {
    let rows = rows_for(func)?;
    let total = rows.len();
    if total == 0 {
        return Err(format!("{func}: no golden rows"));
    }
    let ids: Vec<String> = rows
        .iter()
        .map(|r| {
            r.get("id")
                .and_then(Value::as_str)
                .unwrap_or("?")
                .to_string()
        })
        .collect();
    let (tx, rx) = mpsc::channel::<Result<(), String>>();
    thread::spawn(move || {
        let mut cc = match Consts::new() {
            Ok(cc) => cc,
            Err(e) => {
                let _sent = tx.send(Err(format!("constants: {e:?}")));
                return;
            }
        };
        for row in &rows {
            if tx.send(run_row(row, &mut cc)).is_err() {
                return;
            }
        }
    });
    let start = Instant::now();
    let mut pass = 0usize;
    let mut fails: Vec<String> = Vec::new();
    let mut unsupported = 0usize;
    let mut timed_out = 0usize;
    for id in &ids {
        let left = FN_BUDGET.saturating_sub(start.elapsed());
        if left.is_zero() {
            timed_out = total - pass - fails.len();
            break;
        }
        match rx.recv_timeout(ROW_BUDGET.min(left)) {
            Ok(Ok(())) => pass += 1,
            Ok(Err(why)) => {
                if why.contains("Unsupported") {
                    unsupported += 1;
                }
                fails.push(format!("{id}: {}", short(&why)));
            }
            Err(mpsc::RecvTimeoutError::Timeout) => {
                timed_out = total - pass - fails.len();
                fails.push(format!("{id}: timeout; {timed_out} rows not run"));
                break;
            }
            Err(mpsc::RecvTimeoutError::Disconnected) => {
                timed_out = total - pass - fails.len();
                fails.push(format!("{id}: worker stopped; {timed_out} rows not run"));
                break;
            }
        }
    }
    let failed = total - pass;
    eprintln!(
        "S9 {func}: rows={total} pass={pass} fail={failed} unsupported={unsupported} not_run={timed_out}"
    );
    for line in fails.iter().take(SHOW) {
        eprintln!("  {line}");
    }
    if failed == 0 {
        Ok(())
    } else {
        Err(format!(
            "{func}: {failed}/{total} rows are not the directed correctly rounded result"
        ))
    }
}

#[test]
fn s9_gamma() -> Result<(), String> {
    check("gamma")
}

#[test]
fn s9_ln_gamma() -> Result<(), String> {
    check("ln_gamma")
}

#[test]
fn s9_erf() -> Result<(), String> {
    check("erf")
}

#[test]
fn s9_erfc() -> Result<(), String> {
    check("erfc")
}

#[test]
fn s9_expm1() -> Result<(), String> {
    check("expm1")
}

#[test]
fn s9_log1p() -> Result<(), String> {
    check("log1p")
}

#[test]
fn s9_digamma() -> Result<(), String> {
    check("digamma")
}

#[test]
fn s9_bessel_j0() -> Result<(), String> {
    check("bessel_j0")
}

#[test]
fn s9_bessel_j1() -> Result<(), String> {
    check("bessel_j1")
}

#[test]
fn s9_zeta() -> Result<(), String> {
    check("zeta")
}

#[test]
fn s9_beta() -> Result<(), String> {
    check("beta")
}
