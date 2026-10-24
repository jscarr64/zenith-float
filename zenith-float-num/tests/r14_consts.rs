//! R14 part 1: `Consts::get_exact_prec` returns π, e, ln 2, ln 10, √2, φ and γ correctly
//! rounded at exactly `p` bits in every mode, bit-identical whatever the cache holds and
//! whatever order the requests come in; `Down ≤ c ≤ Up` for every `p ≥ 1`.
//!
//! `tests/oracle/goldens/r14_consts.jsonl` comes from `scripts/gen_r14_const_goldens.py`:
//! mpmath at p+256 and p+512 bits, each error interval settled to one p-bit result and the
//! two agreeing, for p in 1..=65536 (1, 2, 63, 64, 65 and 65536 included) and six modes.

use serde_json::Value;
use zenith_float_num::{ConstantId, Consts, Error, ExactNum, RoundingMode, Sign, Word};

const GOLDS: &str = include_str!("oracle/goldens/r14_consts.jsonl");
const NAMES: [&str; 7] = ["Pi", "E", "Ln2", "Ln10", "Sqrt2", "Phi", "EulerGamma"];
const MODES: [RoundingMode; 6] = [
    RoundingMode::ToEven,
    RoundingMode::ToOdd,
    RoundingMode::Up,
    RoundingMode::Down,
    RoundingMode::ToZero,
    RoundingMode::FromZero,
];
const EDGE_P: [usize; 6] = [1, 2, 63, 64, 65, 65536];
const SHOW: usize = 8;

fn id_of(name: &str) -> Result<ConstantId, String> {
    match name {
        "Pi" => Ok(ConstantId::Pi),
        "E" => Ok(ConstantId::E),
        "Ln2" => Ok(ConstantId::Ln2),
        "Ln10" => Ok(ConstantId::Ln10),
        "Sqrt2" => Ok(ConstantId::Sqrt2),
        "Phi" => Ok(ConstantId::Phi),
        "EulerGamma" => Ok(ConstantId::EulerGamma),
        other => Err(format!("unknown constant {other}")),
    }
}

fn mode_of(name: &str) -> Result<RoundingMode, String> {
    match name {
        "ToEven" => Ok(RoundingMode::ToEven),
        "ToOdd" => Ok(RoundingMode::ToOdd),
        "Up" => Ok(RoundingMode::Up),
        "Down" => Ok(RoundingMode::Down),
        "ToZero" => Ok(RoundingMode::ToZero),
        "FromZero" => Ok(RoundingMode::FromZero),
        other => Err(format!("unknown mode {other}")),
    }
}

/// Sign, exponent and limbs.
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

/// Value, representation and both inexact flags.
fn canon_of(n: &ExactNum) -> String {
    format!(
        "{} p={:?} max={:?} acc={} op={}",
        value_of(n),
        n.precision(),
        n.mantissa_max_bit_len(),
        n.inexact_accumulated(),
        n.inexact()
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

/// The expected value as an `ExactNum` built from its limbs.
fn expect_of(v: &Value) -> Result<ExactNum, String> {
    let arr = field(v, "m")?.as_array().ok_or("m is not an array")?;
    let mut limbs = Vec::with_capacity(arr.len());
    for w in arr {
        let s = w.as_str().ok_or("limb is not a string")?;
        limbs.push(Word::from_str_radix(s, 16).map_err(|e| format!("limb {s}: {e}"))?);
    }
    let e = field(v, "e")?.as_i64().ok_or("e is not an i64")?;
    let sign = if text(v, "s")? == "-" { Sign::Neg } else { Sign::Pos };
    let x = ExactNum::from_words(&limbs, sign, e);
    if x.is_nan() {
        return Err(format!("expect out of range: e={e}"));
    }
    Ok(x)
}

/// (constant, p, mode, expected) for every gold row.
fn rows() -> Result<Vec<(String, usize, String, ExactNum)>, String> {
    let mut out = Vec::new();
    for line in GOLDS.lines().filter(|l| !l.trim().is_empty()) {
        let row: Value = serde_json::from_str(line).map_err(|e| format!("json: {e}"))?;
        let p = field(&row, "p")?.as_u64().ok_or("p is not a u64")?;
        out.push((
            text(&row, "c")?.to_string(),
            usize::try_from(p).map_err(|e| e.to_string())?,
            text(&row, "rm")?.to_string(),
            expect_of(field(&row, "expect")?)?,
        ));
    }
    Ok(out)
}

fn get(cc: &mut Consts, name: &str, p: usize, rm: RoundingMode) -> Result<ExactNum, String> {
    cc.get_exact_prec(id_of(name)?, p, rm)
        .map_err(|e| format!("{name} p={p} {rm:?}: {e:?}"))
}

fn report(bad: &[String], total: usize, what: &str) -> Result<(), String> {
    if bad.is_empty() {
        return Ok(());
    }
    let shown: Vec<&str> = bad.iter().take(SHOW).map(String::as_str).collect();
    Err(format!(
        "{} of {total} {what} failed:\n{}",
        bad.len(),
        shown.join("\n")
    ))
}

#[test]
fn golds_every_constant_every_mode() -> Result<(), String> {
    let rows = rows()?;
    let mut cc = Consts::new().map_err(|e| format!("{e:?}"))?;
    let mut bad = Vec::new();
    for (name, p, rm, want) in &rows {
        let got = get(&mut cc, name, *p, mode_of(rm)?)?;
        // Expect the gold value at exactly p bits, flagged inexact (no constant is exact).
        // The gold value, flagged inexact (no constant is exact).
        let flags = got.inexact() && got.inexact_accumulated();
        if value_of(&got) != value_of(want) || !flags {
            bad.push(format!(
                "{name} p={p} {rm}: got {} want {}",
                canon_of(&got),
                value_of(want)
            ));
        }
    }
    assert_eq!(rows.len(), NAMES.len() * 25 * MODES.len(), "gold row count");
    for name in NAMES {
        for p in EDGE_P {
            assert!(
                rows.iter().any(|(n, q, _, _)| n == name && *q == p),
                "{name} p={p} has golds"
            );
        }
    }
    report(&bad, rows.len(), "gold rows")
}

#[test]
fn down_and_up_bracket_the_constant() -> Result<(), String> {
    // The p=65536 golds bracket each constant to 2^-65535 relative: a reference far finer
    // than any p below it, and exactly the constant's own bracket at p=65536.
    let rows = rows()?;
    let mut cc = Consts::new().map_err(|e| format!("{e:?}"))?;
    let mut bad = Vec::new();
    let mut checked = 0;
    for name in NAMES {
        let reference = |rm: &str| {
            rows.iter()
                .find(|(n, q, m, _)| n == name && *q == 65536 && m == rm)
                .map(|(_, _, _, v)| v)
                .ok_or(format!("{name} has no p=65536 {rm} gold"))
        };
        let (lo, hi) = (reference("Down")?, reference("Up")?);
        let mut ps: Vec<usize> = rows
            .iter()
            .filter(|(n, _, _, _)| n == name)
            .map(|r| r.1)
            .collect();
        ps.dedup();
        for p in ps {
            let down = get(&mut cc, name, p, RoundingMode::Down)?;
            let up = get(&mut cc, name, p, RoundingMode::Up)?;
            checked += 1;
            if !(down < up && down <= *lo && *hi <= up) {
                bad.push(format!(
                    "{name} p={p}: Down {} Up {}",
                    canon_of(&down),
                    canon_of(&up)
                ));
            }
        }
    }
    report(&bad, checked, "brackets")
}

#[test]
fn bit_identical_whatever_the_cache_and_order() -> Result<(), String> {
    let mut bad = Vec::new();
    let mut checked = 0;
    for name in NAMES {
        for rm in MODES {
            // 1000 then 64 against 64 then 1000, each on a fresh cache, and a cache warmed by
            // the word-precision entry points first.
            let mut a = Consts::new().map_err(|e| format!("{e:?}"))?;
            let a1000 = get(&mut a, name, 1000, rm)?;
            let a64 = get(&mut a, name, 64, rm)?;
            let mut b = Consts::new().map_err(|e| format!("{e:?}"))?;
            let b64 = get(&mut b, name, 64, rm)?;
            let b1000 = get(&mut b, name, 1000, rm)?;
            let mut w = Consts::new().map_err(|e| format!("{e:?}"))?;
            for q in [4096, 65, 1] {
                let _ = (w.pi(q, rm), w.e(q, rm), w.ln_2(q, rm), w.ln_10(q, rm));
                let _ = (w.sqrt2(q, rm), w.phi(q, rm), w.euler_gamma(q, rm));
            }
            let w64 = get(&mut w, name, 64, rm)?;
            let w1000 = get(&mut w, name, 1000, rm)?;
            for p in EDGE_P.into_iter().filter(|&p| p < 65536) {
                let mut f = Consts::new().map_err(|e| format!("{e:?}"))?;
                let fresh = get(&mut f, name, p, rm)?;
                let warm = get(&mut w, name, p, rm)?;
                let after = get(&mut a, name, p, rm)?;
                checked += 1;
                if canon_of(&fresh) != canon_of(&warm) || canon_of(&fresh) != canon_of(&after) {
                    bad.push(format!("{name} p={p} {rm:?}: fresh/warm/after differ"));
                }
            }
            checked += 1;
            let k64 = [canon_of(&a64), canon_of(&b64), canon_of(&w64)];
            let k1000 = [canon_of(&a1000), canon_of(&b1000), canon_of(&w1000)];
            if k64[0] != k64[1] || k64[0] != k64[2] || k1000[0] != k1000[1] || k1000[0] != k1000[2]
            {
                bad.push(format!(
                    "{name} {rm:?}: order changed the result\n{k64:?}\n{k1000:?}"
                ));
            }
        }
    }
    report(&bad, checked, "order checks")
}

#[test]
fn no_constant_is_exact() -> Result<(), String> {
    for name in NAMES {
        assert!(!Consts::is_exact(id_of(name)?), "{name}");
    }
    Ok(())
}

#[test]
fn bad_precision_and_mode_are_errors() -> Result<(), String> {
    let mut cc = Consts::new().map_err(|e| format!("{e:?}"))?;
    for name in NAMES {
        assert!(
            matches!(
                cc.get_exact_prec(id_of(name)?, 0, RoundingMode::ToEven),
                Err(Error::InvalidArgument)
            ),
            "{name} p=0"
        );
        assert!(
            matches!(
                cc.get_exact_prec(id_of(name)?, 64, RoundingMode::None),
                Err(Error::InvalidArgument)
            ),
            "{name} None"
        );
    }
    Ok(())
}
