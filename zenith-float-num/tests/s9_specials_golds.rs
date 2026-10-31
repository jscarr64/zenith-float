//! S9 directed special functions (R7, R7b): red goldens.
//!
//! `tests/oracle/goldens/s9_specials.jsonl` comes from
//! `scripts/gen_s9_specials_goldens.py` (mpmath, Ziv pair at p+256 and p+512
//! bits). Every row is limb-identical: one correctly rounded result per
//! rounding mode (`ToEven`, `Up`, `Down`, `ToZero`, `FromZero`) at exactly `p`
//! bits over the i64 exponent range.
//!
//! Digamma rows at p >= 16384, and in debug builds the rows in `DEBUG_FULL_TIER_IDS`, run only
//! with `MPFR_TIER=full` (see `full_tier_only`). `EXCLUDED_IDS` (empty since S3) never run.
//!
//! One test per function. A test fails when any of its rows differs, times out,
//! or has no directed implementation yet (`Unsupported`).

use std::sync::mpsc;
use std::thread;
use std::time::{Duration, Instant};

use serde_json::Value;
use zenith_float_num::{Consts, Error, ExactNum, Exponent, RoundingMode, Sign, Word};

const GOLDS: &str = include_str!("oracle/goldens/s9_specials.jsonl");

/// Wall budget for one row, and for one function's whole table; the full tier gets more.
/// Release fits these: erf(0.5) at p=65536 is 5109 series terms, 2.07 s per row (debug is slower).
const ROW_BUDGET: Duration = Duration::from_secs(10);
const FN_BUDGET: Duration = Duration::from_secs(120);
const FULL_ROW_BUDGET: Duration = Duration::from_secs(120);
const FULL_FN_BUDGET: Duration = Duration::from_secs(900);
/// ZF-T02 slowdown factor for emulated runners, in thousandths (no IEEE float): unset is
/// native (1000). Accepts ASCII digits with at most 3 fraction digits and a value >= 1
/// ("4.27" -> 4270); anything else (non-numeric, inf, NaN, exponent, sign, < 1) is Err,
/// so the run FAILs instead of panicking.
fn time_factor_milli(raw: Option<&str>) -> Result<u64, String> {
    let Some(s) = raw else { return Ok(1000) };
    let bad = || {
        format!(
            "ZF_TIME_FACTOR={s:?}: want a decimal >= 1 with at most 3 fraction digits (e.g. 4.27)"
        )
    };
    let (int, frac, dot) = match s.split_once('.') {
        Some((i, f)) => (i, f, true),
        None => (s, "", false),
    };
    let digits = |t: &str| !t.is_empty() && t.bytes().all(|b| b.is_ascii_digit());
    if !digits(int) || (dot && !digits(frac)) || frac.len() > 3 {
        return Err(bad());
    }
    let whole: u64 = int.parse().map_err(|_| bad())?;
    let mut part: u64 = if dot { frac.parse().map_err(|_| bad())? } else { 0 };
    for _ in frac.len()..3 {
        part = part.checked_mul(10).ok_or_else(bad)?;
    }
    let milli = whole
        .checked_mul(1000)
        .and_then(|m| m.checked_add(part))
        .ok_or_else(bad)?;
    if milli < 1000 {
        return Err(bad());
    }
    Ok(milli)
}

/// Exact budget scaling: native ms times the milli-factor is the scaled budget in µs.
fn scaled(native: Duration, milli: u64) -> Result<Duration, String> {
    let ms = u64::try_from(native.as_millis()).map_err(|e| format!("budget {native:?}: {e}"))?;
    ms.checked_mul(milli)
        .map(Duration::from_micros)
        .ok_or_else(|| format!("budget {native:?} x {milli}/1000 overflows"))
}
/// Mismatch lines printed per function.
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
            zenith_float_num::beta(x, y, p, rm, cc)
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

/// `MPFR_TIER=full` (as in `tools/full_gate.sh`) also runs the full-tier rows.
fn full_tier() -> bool {
    std::env::var("MPFR_TIER").is_ok_and(|t| t == "full")
}

/// Full-tier rows: digamma at p >= 16384. ψ(½) = −γ − 2 ln 2 needs Euler γ at the working
/// precision; since S3 (binary-splitting B1) that is 0.49 s release at p=65536 (was 44 s).
fn full_tier_only(row: &Value) -> bool {
    let high_p_digamma = row.get("fn").and_then(Value::as_str) == Some("digamma")
        && row
            .get("p")
            .and_then(Value::as_u64)
            .is_some_and(|p| p >= 16384);
    high_p_digamma || (cfg!(debug_assertions) && row_id_in(row, DEBUG_FULL_TIER_IDS))
}

/// Debug builds only: rows that exceed the default row budget unoptimized (measured in the
/// debug gate) run in the full tier. Chosen by id, never by timing; release keeps them.
const DEBUG_FULL_TIER_IDS: &[&str] = &[
    "s9-erf-65536-ToEven-1695",
    "s9-erf-65536-Up-1696",
    "s9-erf-65536-Down-1697",
    "s9-erf-65536-ToZero-1698",
    "s9-erf-65536-FromZero-1699",
    "s9-digamma-4096-ToEven-3170",
    "s9-digamma-4096-Up-3171",
    "s9-digamma-4096-Down-3172",
    "s9-digamma-4096-ToZero-3173",
    "s9-digamma-4096-FromZero-3174",
];

/// Rows excluded from every tier by id. Empty: ψ(½)@65536 (s9-digamma-65536-*) runs in the
/// full tier again since S3 made Euler γ 0.49 s (release) at p=65536.
const EXCLUDED_IDS: &[&str] = &[];

fn row_id_in(row: &Value, ids: &[&str]) -> bool {
    row.get("id")
        .and_then(Value::as_str)
        .is_some_and(|id| ids.contains(&id))
}

fn rows_for(func: &str) -> Result<Vec<Value>, String> {
    let mut out = Vec::new();
    for line in GOLDS.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let row: Value = serde_json::from_str(line).map_err(|e| format!("json: {e}"))?;
        if row.get("fn").and_then(Value::as_str) == Some(func)
            && (full_tier() || !full_tier_only(&row))
            && !row_id_in(&row, EXCLUDED_IDS)
        {
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
    let milli = match std::env::var("ZF_TIME_FACTOR") {
        Err(std::env::VarError::NotPresent) => time_factor_milli(None)?,
        Err(e) => return Err(format!("ZF_TIME_FACTOR: {e}")),
        Ok(s) => time_factor_milli(Some(&s))?,
    };
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
    let (native_row, native_fn) = if full_tier() {
        (FULL_ROW_BUDGET, FULL_FN_BUDGET)
    } else {
        (ROW_BUDGET, FN_BUDGET)
    };
    // ZF-T02: an emulated runner (qemu, wine) gets the native budget times its measured
    // slowdown factor (ZF_TIME_FACTOR; aarch64 under qemu measured 4.27 on Pop). Same work
    // allowance, not a looser one; native runs leave it unset (factor 1, budget unchanged).
    let (row_budget, fn_budget) = (scaled(native_row, milli)?, scaled(native_fn, milli)?);
    eprintln!(
        "S9 {func}: budget native row={native_row:?} fn={native_fn:?} factor={milli}/1000 scaled row={row_budget:?} fn={fn_budget:?}"
    );
    let start = Instant::now();
    let mut pass = 0usize;
    let mut fails: Vec<String> = Vec::new();
    let mut unsupported = 0usize;
    let mut timed_out = 0usize;
    for id in &ids {
        let left = fn_budget.saturating_sub(start.elapsed());
        if left.is_zero() {
            timed_out = total - pass - fails.len();
            break;
        }
        match rx.recv_timeout(row_budget.min(left)) {
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

#[test]
fn s9_time_factor_rejects_bad() -> Result<(), String> {
    for bad in [
        "abc", "", "0", "0.999", "0.5", "-1", "+2", "inf", "NaN", "1e3", "4.2701", "1.", ".5",
        "1.2.3", " 2", "8.001", "1000",
    ] {
        if let Ok(v) = time_factor_milli(Some(bad)) {
            return Err(format!("ZF_TIME_FACTOR={bad:?} accepted as {v}"));
        }
    }
    for (raw, want) in [
        (None, 1000u64),
        (Some("1"), 1000),
        (Some("1.34"), 1340),
        (Some("4.27"), 4270),
        (Some("2.5"), 2500),
        (Some("8.000"), 8000),
        (Some("8"), 8000),
        (Some("1.000"), 1000),
        (Some("01.5"), 1500),
    ] {
        let got = time_factor_milli(raw)?;
        if got != want {
            return Err(format!("ZF_TIME_FACTOR={raw:?}: got {got}, want {want}"));
        }
    }
    let row = scaled(ROW_BUDGET, 4270)?;
    if row != Duration::from_millis(42_700) {
        return Err(format!("ROW_BUDGET x 4.27 = {row:?}, want 42.7s"));
    }
    let full = scaled(FULL_FN_BUDGET, 1340)?;
    if full != Duration::from_secs(1206) {
        return Err(format!("FULL_FN_BUDGET x 1.34 = {full:?}, want 1206s"));
    }
    if scaled(Duration::from_secs(u64::MAX), 2000).is_ok() {
        return Err("overflowing budget accepted".to_string());
    }
    // from_secs(u64::MAX) already fails at u64::try_from(as_millis); this one reaches checked_mul.
    if scaled(Duration::from_millis(u64::MAX), 2000).is_ok() {
        return Err("budget ms x milli overflow accepted".to_string());
    }
    Ok(())
}
