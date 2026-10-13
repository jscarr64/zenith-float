//! Step-1 oracle harness.
//!
//! Checked-in JSONL rows were produced by `scripts/gen_oracle_goldens.py` and
//! cross-checked with `scripts/mpfr_oracle.c`. `cargo test --test oracle` runs
//! the fast tier. `ORACLE_TIER=deep` also runs the deep tier. Each row has two
//! seconds; a panic or a timeout is a failure.

use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::mpsc;
use std::thread;
use std::time::Duration;

use zenith_float_num::{
    Consts, Error, ExactNum, Exponent, Radix, RoundingMode, Sign, EXPONENT_MIN, INF_NEG, INF_POS,
};

const CASE_BUDGET: Duration = Duration::from_secs(2);

#[derive(Clone, Debug, PartialEq)]
enum Canon {
    Finite {
        sign: String,
        exp: Exponent,
        limbs: Vec<String>,
    },
    Inf(String),
    Nan(Option<String>),
}

fn tier_selected(tier: &str) -> bool {
    match std::env::var("ORACLE_TIER") {
        Ok(value) if value == "deep" => true,
        Ok(value) if value == "fast" => tier == "fast",
        _ => tier == "fast",
    }
}

fn mode_of(name: &str) -> RoundingMode {
    match name {
        "ToEven" => RoundingMode::ToEven,
        "ToOdd" => RoundingMode::ToOdd,
        "Up" => RoundingMode::Up,
        "Down" => RoundingMode::Down,
        "ToZero" => RoundingMode::ToZero,
        "FromZero" => RoundingMode::FromZero,
        "None" => RoundingMode::None,
        _ => RoundingMode::ToEven,
    }
}

fn err_name(err: Option<Error>) -> Option<String> {
    match err {
        None => None,
        Some(Error::InvalidArgument) => Some("InvalidArgument".to_string()),
        Some(Error::DivisionByZero) => Some("DivisionByZero".to_string()),
        Some(Error::MemoryAllocation) => Some("MemoryAllocation".to_string()),
        Some(Error::PrecisionRetryExhausted) => Some("PrecisionRetryExhausted".to_string()),
        Some(Error::ExponentOverflow(_)) => Some("ExponentOverflow".to_string()),
    }
}

fn storage_words(p: usize) -> usize {
    let rounded = p.div_ceil(64) * 64;
    rounded / 64
}

fn canon_of(n: &ExactNum) -> Canon {
    if n.is_nan() {
        return Canon::Nan(err_name(n.err()));
    }
    if n.is_inf() {
        let sign = if n.sign() == Some(Sign::Neg) { "-" } else { "+" };
        return Canon::Inf(sign.to_string());
    }
    let sign = if n.sign() == Some(Sign::Neg) { "-" } else { "+" };
    let exp = n.exponent().unwrap_or(0);
    let limbs = n
        .mantissa_digits()
        .unwrap_or(&[])
        .iter()
        .map(|word| format!("{word:016x}"))
        .collect();
    Canon::Finite {
        sign: sign.to_string(),
        exp,
        limbs,
    }
}

fn expect_of(v: &serde_json::Value) -> Canon {
    let kind = v.get("k").and_then(|k| k.as_str()).unwrap_or("");
    if kind == "nan" {
        let err = match v.get("err") {
            Some(serde_json::Value::Null) | None => None,
            Some(serde_json::Value::String(s)) => Some(s.clone()),
            Some(_) => Some("?".to_string()),
        };
        return Canon::Nan(err);
    }
    if kind == "inf" {
        return Canon::Inf(
            v.get("s")
                .and_then(|s| s.as_str())
                .unwrap_or("+")
                .to_string(),
        );
    }
    let limbs = v
        .get("m")
        .and_then(|m| m.as_array())
        .map(|arr| {
            arr.iter()
                .filter_map(|w| w.as_str().map(|s| s.to_string()))
                .collect()
        })
        .unwrap_or_default();
    Canon::Finite {
        sign: v
            .get("s")
            .and_then(|s| s.as_str())
            .unwrap_or("+")
            .to_string(),
        exp: v.get("e").and_then(|e| e.as_i64()).unwrap_or(0) as Exponent,
        limbs,
    }
}

fn arg_value(tok: &str, p: usize, cc: &mut Consts) -> ExactNum {
    match tok {
        "+0" | "0" => ExactNum::new(p),
        "-0" => {
            let mut z = ExactNum::new(p);
            z.set_sign(Sign::Neg);
            z
        }
        "+inf" | "inf" => INF_POS,
        "-inf" => INF_NEG,
        "min" | "+min" => ExactNum::min_positive(p),
        "-min" => {
            let mut n = ExactNum::min_positive(p);
            n.set_sign(Sign::Neg);
            n
        }
        "3min" | "+3min" | "-3min" => {
            let mut words = vec![0u64; storage_words(p).max(1)];
            words[0] = 3;
            let mut n = ExactNum::from_words(&words, Sign::Pos, EXPONENT_MIN);
            if tok.starts_with('-') {
                n.set_sign(Sign::Neg);
            }
            n
        }
        "max" | "+max" => ExactNum::max_value(p),
        "-max" => {
            let mut n = ExactNum::max_value(p);
            n.set_sign(Sign::Neg);
            n
        }
        other => ExactNum::parse(other, Radix::Dec, p, RoundingMode::ToEven, cc),
    }
}

fn eval_row(row: &serde_json::Value) -> Result<Canon, String> {
    let op = row.get("op").and_then(|v| v.as_str()).ok_or("missing op")?;
    let p = row.get("p").and_then(|v| v.as_u64()).ok_or("missing p")? as usize;
    let rm = mode_of(row.get("rm").and_then(|v| v.as_str()).unwrap_or("ToEven"));
    let args: Vec<String> = row
        .get("args")
        .and_then(|v| v.as_array())
        .map(|arr| {
            arr.iter()
                .filter_map(|v| v.as_str().map(|s| s.to_string()))
                .collect()
        })
        .unwrap_or_default();
    let mut cc = Consts::new().map_err(|_| "constants".to_string())?;
    let got = match op {
        "parse" => ExactNum::parse(&args[0], Radix::Dec, p, rm, &mut cc),
        "chop" => {
            let src_p = args[1]
                .parse::<usize>()
                .map_err(|_| "src precision".to_string())?;
            let mut n = ExactNum::parse(&args[0], Radix::Dec, src_p, RoundingMode::ToEven, &mut cc);
            match n.set_precision(p, rm) {
                Ok(()) => n,
                Err(Error::ExponentOverflow(sign)) => {
                    if sign.is_negative() {
                        INF_NEG
                    } else {
                        INF_POS
                    }
                }
                Err(err) => return Err(format!("set_precision {err:?}")),
            }
        }
        "add" | "sub" | "mul" | "div" | "rem" | "hypot" | "fma" => {
            let a = arg_value(&args[0], p, &mut cc);
            let b = arg_value(args.get(1).map(String::as_str).unwrap_or("0"), p, &mut cc);
            match op {
                "add" => a.add(&b, p, rm),
                "sub" => a.sub(&b, p, rm),
                "mul" => a.mul(&b, p, rm),
                "div" => a.div(&b, p, rm),
                "rem" => a.rem(&b),
                "hypot" => a.hypot(&b, p, rm),
                "fma" => {
                    let c = arg_value(args.get(2).map(String::as_str).unwrap_or("0"), p, &mut cc);
                    a.fma(&b, &c, p, rm)
                }
                _ => unreachable!(),
            }
        }
        "sqrt" => arg_value(&args[0], p, &mut cc).sqrt(p, rm),
        "cbrt" => arg_value(&args[0], p, &mut cc).cbrt(p, rm),
        other => return Err(format!("unknown op {other}")),
    };
    if row.get("roundtrip").and_then(|v| v.as_bool()) == Some(true)
        && !got.is_nan()
        && !got.is_inf()
    {
        let text = got
            .format(Radix::Dec, rm, &mut cc)
            .map_err(|err| format!("format {err:?}"))?;
        let back = ExactNum::parse(&text, Radix::Dec, p, rm, &mut cc);
        if canon_of(&back) != canon_of(&got) {
            return Err(format!("round-trip {text}"));
        }
    }
    Ok(canon_of(&got))
}

fn run_file(name: &str) {
    let path = format!(
        "{}/tests/oracle/goldens/{name}.jsonl",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(&path).unwrap_or_else(|err| panic!("read {path}: {err}"));
    let mut ran = 0usize;
    let mut failed = 0usize;
    let mut samples: Vec<String> = Vec::new();
    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        let row: serde_json::Value = serde_json::from_str(line).unwrap_or_else(|err| {
            panic!("json {name}: {err}");
        });
        let tier = row.get("tier").and_then(|v| v.as_str()).unwrap_or("fast");
        if !tier_selected(tier) {
            continue;
        }
        ran += 1;
        let id = row
            .get("id")
            .and_then(|v| v.as_str())
            .unwrap_or("?")
            .to_string();
        let expect = expect_of(row.get("expect").unwrap_or(&serde_json::Value::Null));
        let (tx, rx) = mpsc::channel();
        let row_for_thread = row.clone();
        thread::spawn(move || {
            let result = catch_unwind(AssertUnwindSafe(|| eval_row(&row_for_thread)));
            let _ = tx.send(result);
        });
        let outcome = match rx.recv_timeout(CASE_BUDGET) {
            Ok(Ok(Ok(got))) => {
                if got == expect {
                    None
                } else {
                    Some(format!("{id}: got {got:?} expect {expect:?}"))
                }
            }
            Ok(Ok(Err(err))) => Some(format!("{id}: {err}")),
            Ok(Err(_)) => Some(format!("{id}: panic")),
            Err(_) => Some(format!("{id}: timeout")),
        };
        if let Some(msg) = outcome {
            failed += 1;
            if samples.len() < 12 {
                samples.push(msg);
            }
        }
    }
    if failed != 0 {
        panic!(
            "{name}.jsonl fast-or-selected: {failed} failed / {ran}\n{}",
            samples.join("\n")
        );
    }
    eprintln!("{name}.jsonl passed {ran}");
}

#[test]
fn rounding() {
    run_file("rounding");
}

#[test]
fn arith() {
    run_file("arith");
}

#[test]
fn edges() {
    run_file("edges");
}

#[test]
fn parse() {
    run_file("parse");
}
