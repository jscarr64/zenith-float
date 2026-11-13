//! COVERS: common::BallOp::parse, common::BallOp::name, common::goldens::golden_dir
//! L0b-7 (TEST_PLAN section 9): the committed ball enclosure rows (tests/goldens/ball_enclosure_rows.txt, scripts/gen_ball_enclosure.py) at the
//! full precision schedule PS. SCOPE, stated plainly: these are STRUCTURAL checks only (row shape, W = 2*PREC + 64, canonical values, non-negative
//! radii, LO <= HI, a positive RREF, every op at every PS precision). The interval values themselves are guaranteed by the generator's exact
//! Fraction computation and its mpmath.iv cross-check; there is no Ball type in src yet to test against. Consumers assert CONTAINMENT only: the
//! tight/loose checks against RREF need tau and lambda, which are undecided (D2), so nothing here asserts them.

use std::cmp::Ordering;
use std::collections::BTreeSet;

use crate::common::ballrows::BallOp;
use crate::common::goldens::golden_dir;

const FILE: &str = "ball_enclosure_rows.txt";
const PS: [u32; 41] = [
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 20, 24, 31, 32, 33, 48, 53, 63, 64, 65,
    96, 113, 127, 128, 129, 192, 256, 257, 384, 512, 1000, 1024, 2048, 4095, 4096,
];

/// A canonical value: sign, exponent, mantissa bits from the top (no leading zero; empty for zero).
#[derive(Debug, PartialEq, Eq)]
struct Val {
    neg: bool,
    exp: i64,
    bits: String,
}

fn parse_val(s: &str, e: &str, m: &str) -> Result<Val, String> {
    let neg = match s {
        "0" => false,
        "1" => true,
        _ => return Err(format!("bad sign {:?}", s)),
    };
    let exp: i64 = e.parse().map_err(|_| format!("bad exponent {:?}", e))?;
    if m.is_empty() || m.chars().any(|c| !matches!(c, '0'..='9' | 'a'..='f')) {
        return Err(format!("bad mantissa {:?}", m));
    }
    let mut bits = String::new();
    for c in m.chars() {
        let d = c.to_digit(16).ok_or_else(|| format!("bad digit {:?}", c))?;
        bits.push_str(&format!("{:04b}", d));
    }
    let bits = String::from(bits.trim_start_matches('0'));
    if bits.is_empty() {
        if neg || exp != 0 || m != "0" {
            return Err(format!("non-canonical zero {} {} {}", s, e, m));
        }
    } else if !bits.ends_with('1') || m.starts_with('0') {
        return Err(format!("non-canonical mantissa {:?}", m));
    }
    Ok(Val { neg, exp, bits })
}

fn cmp_mag(a: &Val, b: &Val) -> Ordering {
    match (a.bits.is_empty(), b.bits.is_empty()) {
        (true, true) => return Ordering::Equal,
        (true, false) => return Ordering::Less,
        (false, true) => return Ordering::Greater,
        _ => {}
    }
    let by_exp = a.exp.cmp(&b.exp);
    if by_exp != Ordering::Equal {
        return by_exp;
    }
    let n = a.bits.len().max(b.bits.len());
    let pa = format!("{:0<n$}", a.bits, n = n);
    let pb = format!("{:0<n$}", b.bits, n = n);
    pa.cmp(&pb)
}

fn cmp_val(a: &Val, b: &Val) -> Ordering {
    let za = a.bits.is_empty();
    let zb = b.bits.is_empty();
    match (a.neg && !za, b.neg && !zb) {
        (false, true) => Ordering::Greater,
        (true, false) => Ordering::Less,
        (false, false) => cmp_mag(a, b),
        (true, true) => cmp_mag(b, a),
    }
}

struct Row {
    op: BallOp,
    prec: u32,
}

fn take_val<'a, I: Iterator<Item = &'a str>>(it: &mut I) -> Result<Val, String> {
    let s = it.next().ok_or("missing sign")?;
    let e = it.next().ok_or("missing exponent")?;
    let m = it.next().ok_or("missing mantissa")?;
    parse_val(s, e, m)
}

fn check_line(line: &str) -> Result<Row, String> {
    let mut it = line.split(' ');
    let op = BallOp::parse(it.next().ok_or("empty row")?)
        .ok_or_else(|| format!("bad op in {:?}", line))?;
    let prec: u32 = it
        .next()
        .ok_or("missing prec")?
        .parse()
        .map_err(|_| "bad prec")?;
    let w: u64 = it.next().ok_or("missing W")?.parse().map_err(|_| "bad W")?;
    if !PS.contains(&prec) {
        return Err(format!("prec {} is not in PS", prec));
    }
    if w != u64::from(prec) * 2 + 64 {
        return Err(format!("W {} is not 2*{}+64", w, prec));
    }
    let _amid = take_val(&mut it)?;
    let arad = take_val(&mut it)?;
    let mut rads = vec![arad];
    if !matches!(op, BallOp::Sqrt) {
        let _bmid = take_val(&mut it)?;
        rads.push(take_val(&mut it)?);
    }
    let lo = take_val(&mut it)?;
    let hi = take_val(&mut it)?;
    let rref = take_val(&mut it)?;
    if it.next().is_some() {
        return Err(format!("extra fields in {:?}", line));
    }
    if rads.iter().any(|r| r.neg) {
        return Err(format!("negative radius in {:?}", line));
    }
    if cmp_val(&lo, &hi) == Ordering::Greater {
        return Err(format!("LO > HI in {:?}", line));
    }
    if rref.neg || rref.bits.is_empty() || rref.bits.len() > 64 {
        return Err(format!("RREF is not a positive 64-bit value in {:?}", line));
    }
    if lo.bits.len() > 2 * 4096 + 64 || hi.bits.len() > 2 * 4096 + 64 {
        return Err(format!("endpoint wider than W in {:?}", line));
    }
    let wl = usize::try_from(w).map_err(|_| "W too big")?;
    if lo.bits.len() > wl || hi.bits.len() > wl {
        return Err(format!("endpoint wider than W in {:?}", line));
    }
    Ok(Row { op, prec })
}

#[test]
fn every_enclosure_row_is_well_formed() -> Result<(), String> {
    let path = golden_dir().join(FILE);
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{:?}: {}", path, e))?;
    let mut seen = BTreeSet::new();
    let mut n = 0usize;
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let r = check_line(line)?;
        seen.insert((r.op.name(), r.prec));
        n += 1;
    }
    if n < 1000 {
        return Err(format!("only {} enclosure rows", n));
    }
    for op in ["add", "sub", "mul", "div", "sqrt"] {
        for p in PS {
            if !seen.contains(&(op, p)) {
                return Err(format!("no {} row at prec {}", op, p));
            }
        }
    }
    Ok(())
}

#[test]
fn checker_rejects_bad_rows() -> Result<(), String> {
    if BallOp::parse("pow").is_some() || BallOp::parse("sqrt").map(|o| o.name()) != Some("sqrt") {
        return Err(String::from(
            "BallOp::parse accepts the enclosure op names only",
        ));
    }
    let good = "add 8 80 0 1 1 0 0 1 0 1 1 0 0 1 0 1 1 0 2 3 0 1 1";
    check_line(good)?;
    let bad = [
        ("W", "add 8 81 0 1 1 0 0 1 0 1 1 0 0 1 0 1 1 0 2 3 0 1 1"),
        (
            "prec not in PS",
            "add 17 98 0 1 1 0 0 1 0 1 1 0 0 1 0 1 1 0 2 3 0 1 1",
        ),
        (
            "LO > HI",
            "add 8 80 0 1 1 0 0 1 0 1 1 0 0 1 0 2 3 0 1 1 0 1 1",
        ),
        (
            "negative radius",
            "add 8 80 0 1 1 1 0 1 0 1 1 0 0 1 0 1 1 0 2 3 0 1 1",
        ),
        (
            "even mantissa",
            "add 8 80 0 1 2 0 0 1 0 1 1 0 0 1 0 1 1 0 2 3 0 1 1",
        ),
        (
            "zero RREF",
            "add 8 80 0 1 1 0 0 1 0 1 1 0 0 1 0 1 1 0 2 3 0 0 0",
        ),
        (
            "sqrt with B",
            "sqrt 8 80 0 1 1 0 0 1 0 1 1 0 0 1 0 1 1 0 2 3 0 1 1",
        ),
        ("short", "mul 8 80 0 1 1 0 0 1"),
        (
            "negative zero",
            "add 8 80 1 0 0 0 0 1 0 1 1 0 0 1 0 1 1 0 2 3 0 1 1",
        ),
    ];
    for (label, line) in bad {
        if check_line(line).is_ok() {
            return Err(format!("{}: {:?} was accepted", label, line));
        }
    }
    Ok(())
}

#[test]
fn value_order_is_exact() -> Result<(), String> {
    let v = |s: &str, e: &str, m: &str| parse_val(s, e, m);
    // 0.11b * 2^1 = 1.5 > 0.1b * 2^1 = 1 ; -1 < 0 < 1 ; 0.101b*2^2 = 2.5 vs 0.11b*2^2 = 3
    if cmp_val(&v("0", "1", "3")?, &v("0", "1", "1")?) != Ordering::Greater {
        return Err(String::from("1.5 vs 1"));
    }
    if cmp_val(&v("1", "1", "1")?, &v("0", "0", "0")?) != Ordering::Less {
        return Err(String::from("-1 vs 0"));
    }
    if cmp_val(&v("0", "2", "5")?, &v("0", "2", "3")?) != Ordering::Less {
        return Err(String::from("2.5 vs 3"));
    }
    if cmp_val(&v("1", "2", "5")?, &v("1", "2", "3")?) != Ordering::Greater {
        return Err(String::from("-2.5 vs -3"));
    }
    Ok(())
}
