//! R15: directed real `exp2` and `cbrt` in every mode, p 24..65536.
//!
//! `tests/oracle/goldens/r15_exp2_cbrt.jsonl` comes from `scripts/gen_r15_goldens.py`: mpmath
//! at p+256 and p+512 bits (both error intervals round to one p-bit result), or a known exact
//! value (`exp2` of an integer, `cbrt` of a perfect cube of either sign) rounded at p bits.
//! Every row is value-identical at exactly `p` bits (sign, exponent and mantissa bits, compared
//! whatever the word size, so the same golds run on 64- and 32-bit targets) over the i64
//! exponent range (overflow to Inf or MAX by direction, the subnormal grid), and carries the
//! inexact flag the result must have (false on Inf, as `ExactNum::inexact` documents). Fast
//! and deep tiers always run; the full tier (p = 65536) with `MPFR_TIER=full`.

use serde_json::Value;
use zenith_float_num::{
    Consts, ExactNum, Exponent, RoundingMode, Sign, Word, EXPONENT_MAX, EXPONENT_MIN,
};

const GOLDS: &str = include_str!("oracle/goldens/r15_exp2_cbrt.jsonl");
const SHOW: usize = 8;

fn full_tier() -> bool {
    std::env::var("MPFR_TIER").is_ok_and(|t| t == "full")
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
    let raw = field(row, "args")?
        .as_array()
        .ok_or("args is not an array")?;
    let x = arg_of(raw.first().ok_or("no argument")?)?;
    let rm = mode_of(text(row, "rm")?)?;
    let r = match func {
        "exp2" => x.exp2(p, rm, cc),
        "cbrt" => x.cbrt(p, rm),
        other => return Err(format!("unknown function {other}")),
    };
    let want = canon_expect(field(row, "expect")?)?;
    let got = canon_of(&r);
    if got != want {
        return Err(format!("want {want} got {got}"));
    }
    // `ExactNum::inexact` is documented false for Inf (an overflow to Inf carries no flag);
    // overflow to MAX and every underflow keep it.
    let flag = field(row, "inexact")?
        .as_bool()
        .ok_or("inexact is not a bool")?
        && !r.is_inf();
    if r.inexact() != flag {
        return Err(format!("inexact flag {} want {flag}", r.inexact()));
    }
    Ok(())
}

fn check(func: &str, tiers: &[&str]) -> Result<(), String> {
    let mut cc = Consts::new().map_err(|e| format!("consts: {e:?}"))?;
    let mut total = 0usize;
    let mut fails: Vec<String> = Vec::new();
    for line in GOLDS.lines().filter(|l| !l.trim().is_empty()) {
        let row: Value = serde_json::from_str(line).map_err(|e| format!("json: {e}"))?;
        if row.get("fn").and_then(Value::as_str) != Some(func) {
            continue;
        }
        if !tiers.contains(&text(&row, "tier")?) {
            continue;
        }
        total += 1;
        if let Err(why) = run_row(&row, &mut cc) {
            let id = row.get("id").and_then(Value::as_str).unwrap_or("?");
            fails.push(format!("{id}: {}", why.get(..200).unwrap_or(&why)));
        }
    }
    eprintln!(
        "R15 {func} {tiers:?}: rows={total} pass={} fail={}",
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
fn r15_exp2_fast_deep() -> Result<(), String> {
    check("exp2", &["fast", "deep"])
}

#[test]
fn r15_cbrt_fast_deep() -> Result<(), String> {
    check("cbrt", &["fast", "deep"])
}

#[test]
fn r15_full_tier_p65536() -> Result<(), String> {
    if !full_tier() {
        eprintln!("r15 full tier skipped: set MPFR_TIER=full");
        return Ok(());
    }
    check("exp2", &["full"])?;
    check("cbrt", &["full"])
}

/// The golds cover every mode, p from 24 to 65536, exact and Ziv rows, negative cubes, and both
/// ends of the exponent range `±(2^62 − 1)` for each function (in an argument or a result).
#[test]
fn golds_cover_modes_precisions_and_extremes() -> Result<(), String> {
    assert_eq!(EXPONENT_MAX, (1i64 << 62) - 1);
    for func in ["exp2", "cbrt"] {
        let mut modes = std::collections::BTreeSet::new();
        let mut ps = std::collections::BTreeSet::new();
        let (mut top, mut bottom, mut exact, mut ziv, mut neg_exact) = (0, 0, 0, 0, 0);
        for line in GOLDS.lines().filter(|l| !l.trim().is_empty()) {
            let row: Value = serde_json::from_str(line).map_err(|e| e.to_string())?;
            if text(&row, "fn")? != func {
                continue;
            }
            modes.insert(text(&row, "rm")?.to_string());
            ps.insert(field(&row, "p")?.as_u64().ok_or("p")?);
            let arg = field(&row, "args")?.get(0).ok_or("args")?;
            let vals = [arg, field(&row, "expect")?];
            let es = vals
                .iter()
                .filter_map(|v| v.get("e").and_then(Value::as_i64));
            for e in es {
                top += usize::from(e == EXPONENT_MAX);
                bottom += usize::from(e == EXPONENT_MIN);
            }
            match text(&row, "how")? {
                "exact" => {
                    exact += 1;
                    neg_exact += usize::from(arg.get("s").and_then(Value::as_str) == Some("-"));
                }
                _ => ziv += 1,
            }
        }
        assert_eq!(modes.len(), 5, "{func} modes {modes:?}");
        assert_eq!(ps.first(), Some(&24), "{func} p {ps:?}");
        assert_eq!(ps.last(), Some(&65536), "{func} p {ps:?}");
        assert!(
            top > 0 && bottom > 0,
            "{func} extremes top={top} bottom={bottom}"
        );
        assert!(exact > 0 && ziv > 0, "{func} exact={exact} ziv={ziv}");
        if func == "cbrt" {
            assert!(neg_exact > 0, "no negative perfect cube");
        }
    }
    Ok(())
}
