//! ln, log2 and log10 at the edges of the i64 exponent range, every directed mode.
//!
//! `tests/oracle/goldens/log_extreme.tsv` comes from `scripts/gen_log_extreme_goldens.py`:
//! 33 operands (range ends, subnormals, values next to one, powers of 2 and 10) at
//! p = 24, 53, 113, 256, 1024. Each expected value is mpmath at p+256 and p+512 bits
//! (both must round to it), or an exact integer (ln 1, log2 2^k, log10 10^k). A result
//! must match limb for limb at exactly `p` bits.

use zenith_float_num::{Consts, ExactNum, Exponent, RoundingMode, Sign, Word};

const GOLDS: &str = include_str!("oracle/goldens/log_extreme.tsv");
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

/// `f:sign:exponent:limb,limb,...` (low limb first) as a number.
fn arg_of(text: &str) -> Result<ExactNum, String> {
    let mut parts = text.split(':');
    if parts.next() != Some("f") {
        return Err(format!("argument {text} is not finite"));
    }
    let sign = match parts.next() {
        Some("+") => Sign::Pos,
        Some("-") => Sign::Neg,
        other => return Err(format!("bad sign {other:?}")),
    };
    let e: Exponent = parts
        .next()
        .ok_or("no exponent")?
        .parse()
        .map_err(|e| format!("exponent: {e}"))?;
    let mut limbs = Vec::new();
    for w in parts.next().ok_or("no limbs")?.split(',') {
        limbs.push(Word::from_str_radix(w, 16).map_err(|e| format!("limb {w}: {e}"))?);
    }
    let x = ExactNum::from_words(&limbs, sign, e);
    if x.is_nan() {
        return Err(format!("argument out of range: e={e}"));
    }
    Ok(x)
}

fn run_row(line: &str, cc: &mut Consts) -> Result<(), String> {
    let cols: Vec<&str> = line.split('\t').collect();
    let [_id, func, p, rm, arg, want, _how] = cols.as_slice() else {
        return Err(format!("{} columns", cols.len()));
    };
    let p: usize = p.parse().map_err(|e| format!("p: {e}"))?;
    let rm = mode_of(rm)?;
    let x = arg_of(arg)?;
    let r = match *func {
        "ln" => x.ln(p, rm, cc),
        "log2" => x.log2(p, rm, cc),
        "log10" => x.log10(p, rm, cc),
        other => return Err(format!("unknown function {other}")),
    };
    let got = canon_of(&r);
    if got == *want {
        Ok(())
    } else {
        Err(format!("want {want} got {got}"))
    }
}

fn check(func: &str) -> Result<(), String> {
    let mut cc = Consts::new().map_err(|e| format!("consts: {e:?}"))?;
    let mut total = 0usize;
    let mut fails: Vec<String> = Vec::new();
    for line in GOLDS.lines().filter(|l| !l.trim().is_empty()) {
        if line.split('\t').nth(1) != Some(func) {
            continue;
        }
        total += 1;
        if let Err(why) = run_row(line, &mut cc) {
            let id = line.split('\t').next().unwrap_or("?");
            fails.push(format!("{id}: {}", why.get(..240).unwrap_or(&why)));
        }
    }
    eprintln!(
        "log_extreme {func}: rows={total} pass={} fail={}",
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
fn log_extreme_ln() -> Result<(), String> {
    check("ln")
}

#[test]
fn log_extreme_log2() -> Result<(), String> {
    check("log2")
}

#[test]
fn log_extreme_log10() -> Result<(), String> {
    check("log10")
}
