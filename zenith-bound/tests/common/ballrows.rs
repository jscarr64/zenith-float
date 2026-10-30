//! L0b-7 row checker: parses the rows of tests/goldens/ball_exact_ops.txt and re-verifies each one with exact integer arithmetic only.
//!
//! INDEPENDENT of the Python generator (scripts/gen_ball_goldens.py): no Fraction, no mpmath, no goldlib; it shares nothing with it but the
//! row format. SCOPE: `u128` arithmetic limits the full recomputation to rows with `prec <= 64` and `iprec <= 64` (`Checked::Recomputed`);
//! every other row gets only the structural checks (`Checked::Structural`): field syntax and canonical form, mantissa widths, the zero
//! rules, the sign rules of mul, div and sqrt, and a window for the result exponent. Nothing here claims a value check for `prec > 64`.
//! Value model: `(-1)^neg * 0.MANT * 2^exp = mant * 2^(exp - bits)`, `mant` odd, zero is `neg 0 0`. Unbounded exponent range; |exp| <= EXP_BOUND.
//! Every failure is an `Err(String)`; nothing here panics.

use crate::common::generators::RoundingMode;

/// The largest exponent magnitude accepted in a row (the generator stays within 61000).
const EXP_BOUND: i64 = 1_000_000;
/// `round_int` and the recomputation work for precisions up to this (the intermediate integers need `prec + 2 <= 66` bits of room).
pub const MAX_EXACT_PREC: u32 = 64;

#[derive(Debug, PartialEq, Eq)]
pub enum BallOp {
    Add,
    Sub,
    Mul,
    Div,
    Sqrt,
}

impl BallOp {
    pub fn parse(name: &str) -> Option<BallOp> {
        match name {
            "add" => Some(BallOp::Add),
            "sub" => Some(BallOp::Sub),
            "mul" => Some(BallOp::Mul),
            "div" => Some(BallOp::Div),
            "sqrt" => Some(BallOp::Sqrt),
            _ => None,
        }
    }

    pub fn name(&self) -> &'static str {
        match self {
            BallOp::Add => "add",
            BallOp::Sub => "sub",
            BallOp::Mul => "mul",
            BallOp::Div => "div",
            BallOp::Sqrt => "sqrt",
        }
    }
}

/// A parsed value. `bits` is the bit length of the mantissa (0 for zero); `mant` is `None` when it does not fit `u128`.
#[derive(Debug, PartialEq, Eq)]
pub struct Num {
    pub neg: bool,
    pub exp: i64,
    pub bits: u32,
    pub mant: Option<u128>,
}

/// One golden row. `b` is `None` exactly for sqrt.
#[derive(Debug)]
pub struct GoldRow {
    pub op: BallOp,
    pub prec: u32,
    pub rm: RoundingMode,
    pub iprec: u32,
    pub a: Num,
    pub b: Option<Num>,
    pub r: Num,
    pub inexact: bool,
}

/// How far a row was verified.
#[derive(Debug, PartialEq, Eq)]
pub enum Checked {
    /// The result and the inexact flag were recomputed from the operands with exact integers and agree.
    Recomputed,
    /// Only the structural checks ran (`prec > 64` or `iprec > 64`).
    Structural,
}

/// `round_int` result: `m * 2^t` with `m` odd, or `m == 0` (an exact zero, `t == 0`).
#[derive(Debug, PartialEq, Eq)]
pub struct Rounded {
    pub m: u128,
    pub t: i64,
    pub inexact: bool,
}

// ---- parsing ------------------------------------------------------------------------------------------------------------------

fn canon_dec<T: core::str::FromStr + core::fmt::Display>(
    what: &str,
    tok: &str,
) -> Result<T, String> {
    let v: T = tok
        .parse()
        .map_err(|_| format!("{}: {:?} is not a number", what, tok))?;
    if v.to_string() != tok {
        return Err(format!("{}: {:?} is not canonical decimal", what, tok));
    }
    Ok(v)
}

fn parse_num(sign: &str, exp: &str, hex: &str, what: &str) -> Result<Num, String> {
    let neg = match sign {
        "0" => false,
        "1" => true,
        _ => return Err(format!("{}: bad sign {:?}", what, sign)),
    };
    let exp: i64 = canon_dec(what, exp)?;
    if exp.abs() > EXP_BOUND {
        return Err(format!("{}: exponent {} out of bound", what, exp));
    }
    let mut chars = hex.chars();
    let first = chars
        .next()
        .ok_or_else(|| format!("{}: empty mantissa", what))?;
    if hex
        .chars()
        .any(|c| !(c.is_ascii_digit() || ('a'..='f').contains(&c)))
    {
        return Err(format!("{}: mantissa {:?} is not lowercase hex", what, hex));
    }
    let d = first
        .to_digit(16)
        .ok_or_else(|| format!("{}: bad digit", what))?;
    if hex == "0" {
        if exp != 0 {
            return Err(format!("{}: zero with exponent {}", what, exp));
        }
        return Ok(Num {
            neg,
            exp: 0,
            bits: 0,
            mant: Some(0),
        });
    }
    if d == 0 {
        return Err(format!("{}: leading zero in {:?}", what, hex));
    }
    let digits = u32::try_from(hex.len()).map_err(|_| format!("{}: mantissa too long", what))?;
    let bits = 4 * (digits - 1) + (32 - d.leading_zeros());
    let mant = if bits <= 128 {
        Some(
            u128::from_str_radix(hex, 16)
                .map_err(|_| format!("{}: bad mantissa {:?}", what, hex))?,
        )
    } else {
        None
    };
    let last = hex
        .chars()
        .last()
        .and_then(|c| c.to_digit(16))
        .ok_or_else(|| format!("{}: bad digit", what))?;
    if last % 2 == 0 {
        return Err(format!(
            "{}: mantissa {:?} is even (not canonical)",
            what, hex
        ));
    }
    Ok(Num {
        neg,
        exp,
        bits,
        mant,
    })
}

/// Parses one non-comment row: 14 fields (add, sub, mul, div) or 11 (sqrt), canonical single-space text.
pub fn parse_row(line: &str) -> Result<GoldRow, String> {
    if line.split(' ').any(|t| t.is_empty()) {
        return Err(format!("empty field (repeated space) in {:?}", line));
    }
    let f: Vec<&str> = line.split(' ').collect();
    let op_name = f.first().ok_or_else(|| String::from("empty row"))?;
    let op = BallOp::parse(op_name).ok_or_else(|| format!("unknown op {:?}", op_name))?;
    let want = if op == BallOp::Sqrt { 11 } else { 14 };
    if f.len() != want {
        return Err(format!(
            "{} row has {} fields, want {}",
            op.name(),
            f.len(),
            want
        ));
    }
    let get = |i: usize| match f.get(i) {
        Some(&s) => Ok(s),
        None => Err(format!("field {} missing", i)),
    };
    let prec: u32 = canon_dec("prec", get(1)?)?;
    let rm = RoundingMode::parse(get(2)?)
        .ok_or_else(|| format!("unknown rounding mode {:?}", get(2).unwrap_or("")))?;
    let iprec: u32 = canon_dec("iprec", get(3)?)?;
    if prec == 0 || iprec == 0 {
        return Err(String::from("precision 0"));
    }
    let a = parse_num(get(4)?, get(5)?, get(6)?, "a")?;
    let (b, rest) = if op == BallOp::Sqrt {
        (None, 7)
    } else {
        (Some(parse_num(get(7)?, get(8)?, get(9)?, "b")?), 10)
    };
    let r = parse_num(get(rest)?, get(rest + 1)?, get(rest + 2)?, "r")?;
    let inexact = match get(rest + 3)? {
        "0" => false,
        "1" => true,
        other => return Err(format!("bad inexact flag {:?}", other)),
    };
    Ok(GoldRow {
        op,
        prec,
        rm,
        iprec,
        a,
        b,
        r,
        inexact,
    })
}

fn num_text(n: &Num) -> Result<String, String> {
    let sign = if n.neg { 1 } else { 0 };
    match n.mant {
        Some(m) => Ok(format!("{} {} {:x}", sign, n.exp, m)),
        None => Err(String::from(
            "a mantissa wider than 128 bits has no text here",
        )),
    }
}

/// The canonical text of a row with every mantissa at most 128 bits (the generator emits wider ones for `prec > 64`; those rows are
/// compared by text in the harness test only through `parse_row` field counts, so this returns `Err` for them).
pub fn format_row(row: &GoldRow) -> Result<String, String> {
    let mut s = format!(
        "{} {} {} {} {}",
        row.op.name(),
        row.prec,
        row.rm.name(),
        row.iprec,
        num_text(&row.a)?
    );
    if let Some(b) = &row.b {
        s.push(' ');
        s.push_str(&num_text(b)?);
    }
    s.push(' ');
    s.push_str(&num_text(&row.r)?);
    s.push(' ');
    s.push(if row.inexact { '1' } else { '0' });
    Ok(s)
}

// ---- exact rounding -----------------------------------------------------------------------------------------------------------

enum Frac {
    Zero,
    Below,
    Tie,
    Above,
}

fn canonical(m: u128, t: i64, inexact: bool) -> Rounded {
    if m == 0 {
        return Rounded {
            m: 0,
            t: 0,
            inexact,
        };
    }
    let tz = m.trailing_zeros();
    Rounded {
        m: m >> tz,
        t: t + i64::from(tz),
        inexact,
    }
}

/// Rounds `(-1)^neg * (n + f) * 2^t` to `prec` significant bits under `rm`, where `f` is a fraction in `(0, 1)` when `sticky` and 0
/// otherwise. A nonzero `f` needs `n` to have more than `prec` bits (so the discarded part always contains the leading unit of `f`'s
/// influence). `prec` is 1..=64. A nearest tie is resolved by the unit count of the LOWER neighbour: ToEven takes the upper one iff it is
/// odd, ToOdd iff it is even (the upper neighbour of an all-ones mantissa is a carry to a power of two, which counts as even).
pub fn round_int(
    neg: bool,
    n: u128,
    sticky: bool,
    t: i64,
    prec: u32,
    rm: &RoundingMode,
) -> Result<Rounded, String> {
    if prec == 0 || prec > MAX_EXACT_PREC {
        return Err(format!(
            "round_int: prec {} outside 1..={}",
            prec, MAX_EXACT_PREC
        ));
    }
    if n == 0 {
        if sticky {
            return Err(String::from("round_int: a fraction without a leading unit"));
        }
        return Ok(canonical(0, 0, false));
    }
    let len = 128 - n.leading_zeros();
    if len <= prec {
        if sticky {
            return Err(String::from(
                "round_int: sticky fraction needs more than prec bits",
            ));
        }
        return Ok(canonical(n, t, false));
    }
    let drop = len - prec;
    let q = n >> drop;
    let rem = n & ((1u128 << drop) - 1);
    let half = 1u128 << (drop - 1);
    let frac = if rem > half || (rem == half && sticky) {
        Frac::Above
    } else if rem == half {
        Frac::Tie
    } else if rem == 0 && !sticky {
        Frac::Zero
    } else {
        Frac::Below
    };
    let inexact = !matches!(frac, Frac::Zero);
    let tie = matches!(frac, Frac::Tie);
    let above = matches!(frac, Frac::Above);
    let q_odd = q % 2 == 1;
    let up = match rm {
        RoundingMode::ToEven => above || (tie && q_odd),
        RoundingMode::ToOdd => above || (tie && !q_odd),
        RoundingMode::Up => inexact && !neg,
        RoundingMode::Down => inexact && neg,
        RoundingMode::ToZero => false,
        RoundingMode::FromZero => inexact,
    };
    let kept = if up { q + 1 } else { q };
    Ok(canonical(kept, t + i64::from(drop), inexact))
}

// ---- recomputation ------------------------------------------------------------------------------------------------------------

#[derive(Debug, PartialEq, Eq)]
enum Out {
    Zero { neg: bool },
    Fin { neg: bool, r: Rounded },
}

fn from_rounded(neg: bool, r: Rounded) -> Out {
    if r.m == 0 {
        Out::Zero { neg }
    } else {
        Out::Fin { neg, r }
    }
}

/// An operand of the recomputation: sign, mantissa and the exponent of `0.m * 2^exp`; the bit length is derived from `m`.
struct Operand {
    neg: bool,
    m: u128,
    exp: i64,
}

impl Operand {
    fn bits(&self) -> u32 {
        128 - self.m.leading_zeros()
    }
}

fn operand_of(n: &Num, neg: bool) -> Result<Operand, String> {
    let m = n
        .mant
        .ok_or_else(|| String::from("mantissa too wide to recompute"))?;
    if n.bits > 64 {
        return Err(String::from("operand wider than 64 bits"));
    }
    Ok(Operand { neg, m, exp: n.exp })
}

/// The mantissa shifted so its top bit is bit 63.
fn top64(x: &Operand) -> u128 {
    x.m << (64 - x.bits())
}

fn add_mags(x: &Operand, y: &Operand, prec: u32, rm: &RoundingMode) -> Result<Out, String> {
    if x.m == 0 && y.m == 0 {
        let neg = if x.neg == y.neg {
            x.neg
        } else {
            matches!(rm, RoundingMode::Down)
        };
        return Ok(Out::Zero { neg });
    }
    if x.m == 0 {
        return Ok(from_rounded(
            y.neg,
            round_int(y.neg, y.m, false, y.exp - i64::from(y.bits()), prec, rm)?,
        ));
    }
    if y.m == 0 {
        return Ok(from_rounded(
            x.neg,
            round_int(x.neg, x.m, false, x.exp - i64::from(x.bits()), prec, rm)?,
        ));
    }
    let (xa, ya) = (top64(x), top64(y));
    let (hi, lo, ha, la) = if x.exp > y.exp || (x.exp == y.exp && xa >= ya) {
        (x, y, xa, ya)
    } else {
        (y, x, ya, xa)
    };
    let gap = u32::try_from(hi.exp - lo.exp).map_err(|_| String::from("exponent order"))?;
    let aw = ha << 62;
    let (bw, sticky) = if gap >= 126 {
        (0u128, true)
    } else {
        let l = la << 62;
        (l >> gap, l & ((1u128 << gap) - 1) != 0)
    };
    let t = hi.exp - 126;
    if hi.neg == lo.neg {
        return Ok(from_rounded(
            hi.neg,
            round_int(hi.neg, aw + bw, sticky, t, prec, rm)?,
        ));
    }
    if gap == 0 && ha == la {
        return Ok(Out::Zero {
            neg: matches!(rm, RoundingMode::Down),
        });
    }
    Ok(from_rounded(
        hi.neg,
        round_int(hi.neg, aw - bw - u128::from(sticky), sticky, t, prec, rm)?,
    ))
}

fn mul_mags(x: &Operand, y: &Operand, prec: u32, rm: &RoundingMode) -> Result<Out, String> {
    let neg = x.neg != y.neg;
    if x.m == 0 || y.m == 0 {
        return Ok(Out::Zero { neg });
    }
    let p =
        x.m.checked_mul(y.m)
            .ok_or_else(|| String::from("product overflow"))?;
    let t = (x.exp - i64::from(x.bits())) + (y.exp - i64::from(y.bits()));
    Ok(from_rounded(neg, round_int(neg, p, false, t, prec, rm)?))
}

fn div_mags(x: &Operand, y: &Operand, prec: u32, rm: &RoundingMode) -> Result<Out, String> {
    let neg = x.neg != y.neg;
    if y.m == 0 {
        return Err(String::from("zero divisor"));
    }
    if x.m == 0 {
        return Ok(Out::Zero { neg });
    }
    let (a, b) = (top64(x), top64(y));
    let num = a << 64;
    let (q1, r1) = (num / b, num % b);
    let (q2, r2) = ((r1 << 2) / b, (r1 << 2) % b);
    let n = (q1 << 2) + q2;
    let t = x.exp - y.exp - 66;
    Ok(from_rounded(neg, round_int(neg, n, r2 != 0, t, prec, rm)?))
}

/// floor(sqrt(n)) and n - floor(sqrt(n))^2, digit by digit.
fn isqrt_rem(n: u128) -> (u128, u128) {
    let mut x = n;
    let mut res = 0u128;
    let mut bit = 1u128 << 126;
    while bit > x {
        bit >>= 2;
    }
    while bit != 0 {
        if x >= res + bit {
            x -= res + bit;
            res = (res >> 1) + bit;
        } else {
            res >>= 1;
        }
        bit >>= 2;
    }
    (res, x)
}

fn sqrt_mag(x: &Operand, prec: u32, rm: &RoundingMode) -> Result<Out, String> {
    if x.m == 0 {
        return Ok(Out::Zero { neg: x.neg });
    }
    if x.neg {
        return Err(String::from("sqrt of a negative nonzero value"));
    }
    let mut m = x.m;
    let mut t = x.exp - i64::from(x.bits());
    if t.rem_euclid(2) == 1 {
        m <<= 1;
        t -= 1;
    }
    let len = 128 - m.leading_zeros();
    let k = (127 - len) / 2;
    let scaled = m << (2 * k);
    let t = t - 2 * i64::from(k);
    let (mut r, mut rem) = isqrt_rem(scaled);
    for _ in 0..3 {
        let four_rem = rem * 4;
        let thresh = r * 4 + 1;
        if four_rem >= thresh {
            rem = four_rem - thresh;
            r = 2 * r + 1;
        } else {
            rem = four_rem;
            r *= 2;
        }
    }
    Ok(from_rounded(
        false,
        round_int(false, r, rem != 0, t / 2 - 3, prec, rm)?,
    ))
}

fn recompute(row: &GoldRow) -> Result<Out, String> {
    let a = operand_of(&row.a, row.a.neg)?;
    let two = |b: &Option<Num>| -> Result<Operand, String> {
        let n = b
            .as_ref()
            .ok_or_else(|| String::from("missing operand b"))?;
        operand_of(n, n.neg)
    };
    match row.op {
        BallOp::Add => add_mags(&a, &two(&row.b)?, row.prec, &row.rm),
        BallOp::Sub => {
            let n = row
                .b
                .as_ref()
                .ok_or_else(|| String::from("missing operand b"))?;
            add_mags(&a, &operand_of(n, !n.neg)?, row.prec, &row.rm)
        }
        BallOp::Mul => mul_mags(&a, &two(&row.b)?, row.prec, &row.rm),
        BallOp::Div => div_mags(&a, &two(&row.b)?, row.prec, &row.rm),
        BallOp::Sqrt => sqrt_mag(&a, row.prec, &row.rm),
    }
}

// ---- the checks ---------------------------------------------------------------------------------------------------------------

fn structural(row: &GoldRow) -> Result<(), String> {
    if row.a.bits > row.iprec {
        return Err(format!(
            "a has {} bits, more than iprec {}",
            row.a.bits, row.iprec
        ));
    }
    if let Some(b) = &row.b {
        if b.bits > row.iprec {
            return Err(format!(
                "b has {} bits, more than iprec {}",
                b.bits, row.iprec
            ));
        }
    }
    if row.r.bits > row.prec {
        return Err(format!(
            "r has {} bits, more than prec {}",
            row.r.bits, row.prec
        ));
    }
    if row.r.bits == 0 && row.inexact {
        return Err(String::from("an inexact row with a zero result"));
    }
    let r_zero = row.r.bits == 0;
    let a_zero = row.a.bits == 0;
    let b_zero = row.b.as_ref().map(|b| b.bits == 0).unwrap_or(false);
    let b_neg = row.b.as_ref().map(|b| b.neg).unwrap_or(false);
    let (ea, eb) = (row.a.exp, row.b.as_ref().map(|b| b.exp).unwrap_or(0));
    match row.op {
        BallOp::Add | BallOp::Sub => {
            if !r_zero && !a_zero && !b_zero && row.r.exp > ea.max(eb) + 2 {
                return Err(String::from("add result exponent above max(ea, eb) + 2 (a sum below 2^(e+1) can round up to it)"));
            }
        }
        BallOp::Mul => {
            if r_zero != (a_zero || b_zero) {
                return Err(String::from("mul result is zero iff an operand is zero"));
            }
            if row.r.neg != (row.a.neg != b_neg) {
                return Err(String::from(
                    "mul result sign is not the xor of the operand signs",
                ));
            }
            if !r_zero && (row.r.exp < ea + eb - 1 || row.r.exp > ea + eb + 1) {
                return Err(String::from("mul result exponent outside ea + eb -1..=+1"));
            }
        }
        BallOp::Div => {
            if b_zero {
                return Err(String::from("zero divisor"));
            }
            if r_zero != a_zero {
                return Err(String::from("div result is zero iff the dividend is zero"));
            }
            if row.r.neg != (row.a.neg != b_neg) {
                return Err(String::from(
                    "div result sign is not the xor of the operand signs",
                ));
            }
            if !r_zero && (row.r.exp < ea - eb || row.r.exp > ea - eb + 2) {
                return Err(String::from("div result exponent outside ea - eb..=+2"));
            }
        }
        BallOp::Sqrt => {
            if row.a.neg && !a_zero {
                return Err(String::from("sqrt of a negative nonzero value"));
            }
            if r_zero != a_zero || row.r.neg != row.a.neg {
                return Err(String::from("sqrt zero or sign rule violated"));
            }
            if !r_zero && (2 * row.r.exp < ea || 2 * row.r.exp > ea + 3) {
                return Err(String::from(
                    "sqrt result exponent outside ea..=ea+3 (doubled)",
                ));
            }
        }
    }
    Ok(())
}

fn expected_out(row: &GoldRow) -> Result<Out, String> {
    if row.r.bits == 0 {
        return Ok(Out::Zero { neg: row.r.neg });
    }
    let m = row
        .r
        .mant
        .ok_or_else(|| String::from("result mantissa too wide"))?;
    Ok(Out::Fin {
        neg: row.r.neg,
        r: Rounded {
            m,
            t: row.r.exp - i64::from(row.r.bits),
            inexact: row.inexact,
        },
    })
}

/// Verifies one parsed row: the structural checks always, the exact recomputation when `prec <= 64` and `iprec <= 64`.
pub fn check_row(row: &GoldRow) -> Result<Checked, String> {
    structural(row)?;
    if row.prec > MAX_EXACT_PREC || row.iprec > MAX_EXACT_PREC {
        return Ok(Checked::Structural);
    }
    let got = recompute(row)?;
    let want = expected_out(row)?;
    if got != want {
        return Err(format!("recomputed {:?}, row says {:?}", got, want));
    }
    Ok(Checked::Recomputed)
}
