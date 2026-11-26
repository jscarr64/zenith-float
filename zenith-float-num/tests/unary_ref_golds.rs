//! The borrowed unary kernels `ExactNum::{abs_ref, neg_ref, floor_ref, ldexp_ref}`: the operand is only
//! read and the result's limbs are written once (ldexp and floor round straight into the result).
//!
//! `tests/oracle/goldens/unary_ref.jsonl` comes from `scripts/gen_unary_ref_goldens.py`: the exact result
//! from FLINT (arf/arb), the same exact result from mpmath, then an exact integer model rounds it at `p`
//! bits in all 7 rounding modes (mpmath's own rounding cross-checks the in-range IEEE modes). Rows cover
//! signed zeros, Inf, NaN, subnormal operands and results, exact and inexact ldexp, ties, overflow at
//! `EXPONENT_MAX`, underflow below `EXPONENT_MIN` and the `Exponent` extremes.
//!
//! Every row also runs with the operand's flags in each legal state, and the borrowed result must be
//! bit-identical (raw parts, both inexact flags, error kind) to the by-value operation on a fresh copy.

use serde_json::Value;
use zenith_float_num::{
    Error, ExactNum, Exponent, RawKind, RawParts, RoundingMode, Sign, Word, INF_NEG, INF_POS,
};

const GOLDS: &str = include_str!("oracle/goldens/unary_ref.jsonl");
const SHOW: usize = 8;

fn mode_of(name: &str) -> Result<RoundingMode, String> {
    match name {
        "ToEven" => Ok(RoundingMode::ToEven),
        "ToOdd" => Ok(RoundingMode::ToOdd),
        "Up" => Ok(RoundingMode::Up),
        "Down" => Ok(RoundingMode::Down),
        "ToZero" => Ok(RoundingMode::ToZero),
        "FromZero" => Ok(RoundingMode::FromZero),
        "None" => Ok(RoundingMode::None),
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

/// Canonical text of a value; a zero is `zero:<sign>` (its stored exponent is not part of the value).
fn canon_of(n: &ExactNum) -> String {
    if n.is_nan() {
        return format!("nan:{:?}", n.err());
    }
    if n.is_inf() {
        return format!("inf:{}", sign_str(n));
    }
    if n.is_zero() {
        return format!("zero:{}", sign_str(n));
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

/// Everything observable about a number: raw parts with both flags, and the NaN error kind.
fn full_of(n: &ExactNum) -> String {
    format!("{:?} err={:?}", n.as_raw_parts(), n.err())
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
        "nan" => Ok(format!("nan:Some({})", text(v, "err")?)),
        "inf" => Ok(format!("inf:{}", text(v, "s")?)),
        "f" => {
            let limbs = limbs_of(v)?;
            if limbs.iter().all(|w| *w == 0) {
                return Ok(format!("zero:{}", text(v, "s")?));
            }
            let shown: Vec<String> = limbs.iter().map(|w| format!("{w:016x}")).collect();
            let e: Exponent = field(v, "e")?.as_i64().ok_or("e is not an i64")?;
            Ok(format!("f:{}:{}:{}", text(v, "s")?, e, shown.join(",")))
        }
        other => Err(format!("unknown kind {other}")),
    }
}

/// The operand of a row with the flags `(inexact, op_inexact)`; `None` when that flag state is not
/// legal for it (Inf and NaN carry no flags).
fn arg_of(v: &Value, flags: (bool, bool)) -> Result<Option<ExactNum>, String> {
    match text(v, "k")? {
        "nan" => {
            return Ok(
                (flags == (false, false)).then(|| ExactNum::nan(Some(Error::InvalidArgument)))
            )
        }
        "inf" => {
            let inf = if text(v, "s")? == "-" { INF_NEG } else { INF_POS };
            return Ok((flags == (false, false)).then_some(inf));
        }
        _ => {}
    }
    let limbs = limbs_of(v)?;
    let negative = text(v, "s")? == "-";
    let zero = limbs.iter().all(|w| *w == 0);
    let e: Exponent = field(v, "e")?.as_i64().ok_or("e is not an i64")?;
    let parts = RawParts {
        kind: if zero { RawKind::Zero } else { RawKind::Finite },
        words: &limbs,
        precision: limbs.len() * zenith_float_num::WORD_BIT_SIZE,
        negative,
        exponent: if zero { 0 } else { e },
        inexact: flags.0,
        op_inexact: flags.1,
    };
    ExactNum::from_raw_parts(&parts)
        .map(Some)
        .map_err(|err| format!("argument {parts:?}: {err:?}"))
}

/// A fresh number equal to `x` bit for bit, to be moved into the by-value operation.
fn fresh(x: &ExactNum) -> Result<ExactNum, String> {
    if x.is_nan() {
        return Ok(ExactNum::nan(x.err()));
    }
    if x.is_inf() {
        return Ok(if x.sign() == Some(Sign::Neg) { INF_NEG } else { INF_POS });
    }
    ExactNum::from_raw_parts(&x.as_raw_parts()).map_err(|e| format!("fresh: {e:?}"))
}

struct Row<'a> {
    func: String,
    p: usize,
    rm: RoundingMode,
    k: Exponent,
    arg: &'a Value,
    want: String,
}

fn row_of(row: &Value) -> Result<Row<'_>, String> {
    let raw = field(row, "args")?
        .as_array()
        .ok_or("args is not an array")?;
    Ok(Row {
        func: text(row, "fn")?.to_string(),
        p: usize::try_from(field(row, "p")?.as_u64().ok_or("p is not a u64")?)
            .map_err(|e| format!("p: {e}"))?,
        rm: mode_of(text(row, "rm")?)?,
        k: field(row, "k")?.as_i64().ok_or("k is not an i64")?,
        arg: raw.first().ok_or("no argument")?,
        want: canon_expect(field(row, "expect")?)?,
    })
}

fn apply(r: &Row<'_>, x: &ExactNum) -> Result<(ExactNum, ExactNum), String> {
    let moved = fresh(x)?;
    Ok(match r.func.as_str() {
        "abs" => (x.abs_ref(), moved.abs()),
        "neg" => (x.neg_ref(), -moved),
        "floor" => (x.floor_ref(), moved.floor()),
        "ldexp" => (x.ldexp_ref(r.k, r.p, r.rm), moved.ldexp(r.k, r.p, r.rm)),
        other => return Err(format!("unknown function {other}")),
    })
}

fn run_row(row: &Value) -> Result<usize, String> {
    let r = row_of(row)?;
    let mut ran = 0;
    for flags in [(false, false), (true, false), (true, true)] {
        let Some(x) = arg_of(r.arg, flags)? else {
            continue;
        };
        let before = full_of(&x);
        let (got, by_val) = apply(&r, &x)?;
        if full_of(&x) != before {
            return Err(format!("{flags:?}: the operand changed"));
        }
        if canon_of(&got) != r.want {
            return Err(format!("{flags:?}: want {} got {}", r.want, canon_of(&got)));
        }
        if full_of(&got) != full_of(&by_val) {
            return Err(format!(
                "{flags:?}: borrowed {} differs from by-value {}",
                full_of(&got),
                full_of(&by_val)
            ));
        }
        ran += 1;
    }
    Ok(ran)
}

fn check(func: &str) -> Result<(), String> {
    let (mut total, mut runs) = (0usize, 0usize);
    let mut fails: Vec<String> = Vec::new();
    for line in GOLDS.lines().filter(|l| !l.trim().is_empty()) {
        let row: Value = serde_json::from_str(line).map_err(|e| format!("json: {e}"))?;
        if row.get("fn").and_then(Value::as_str) != Some(func) {
            continue;
        }
        total += 1;
        match run_row(&row) {
            Ok(n) => runs += n,
            Err(why) => {
                let id = row.get("id").and_then(Value::as_str).unwrap_or("?");
                fails.push(format!("{id}: {}", why.get(..400).unwrap_or(&why)));
            }
        }
    }
    eprintln!(
        "unary {func}: rows={total} runs={runs} pass={} fail={}",
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
fn unary_abs_ref() -> Result<(), String> {
    check("abs")
}

#[test]
fn unary_neg_ref() -> Result<(), String> {
    check("neg")
}

#[test]
fn unary_floor_ref() -> Result<(), String> {
    check("floor")
}

#[test]
fn unary_ldexp_ref() -> Result<(), String> {
    check("ldexp")
}

/// The ldexp golds cover all 7 modes, signed zeros, Inf, NaN, overflow and underflow.
#[test]
fn golds_cover_the_edges() -> Result<(), String> {
    let (mut modes, mut kinds) = (Vec::new(), Vec::new());
    let (mut neg_zero, mut over, mut under, mut sub_out) = (0, 0, 0, 0);
    for line in GOLDS.lines().filter(|l| !l.trim().is_empty()) {
        let row: Value = serde_json::from_str(line).map_err(|e| e.to_string())?;
        if text(&row, "fn")? != "ldexp" {
            continue;
        }
        let rm = text(&row, "rm")?.to_string();
        if !modes.contains(&rm) {
            modes.push(rm);
        }
        let arg = field(&row, "args")?.as_array().ok_or("args")?;
        let a = arg.first().ok_or("arg")?;
        let k = text(a, "k")?.to_string();
        if !kinds.contains(&k) {
            kinds.push(k);
        }
        let want = canon_expect(field(&row, "expect")?)?;
        neg_zero += usize::from(canon_expect(a)? == "zero:-");
        over += usize::from(want.starts_with("inf:") && text(a, "k")? == "f");
        let e = field(&row, "expect")?.get("e").and_then(Value::as_i64);
        sub_out += usize::from(e == Some(zenith_float_num::EXPONENT_MIN));
        under += usize::from(want.starts_with("zero:") && canon_expect(a)?.starts_with("f:"));
    }
    if modes.len() != 7
        || kinds.len() != 3
        || neg_zero == 0
        || over == 0
        || under == 0
        || sub_out == 0
    {
        return Err(format!(
            "modes={modes:?} kinds={kinds:?} neg_zero={neg_zero} over={over} under={under} sub_out={sub_out}"
        ));
    }
    Ok(())
}
