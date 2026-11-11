//! ACCEPTANCE-pr39 I10 (T3): exhaustive small-precision enumeration for the binary ops whose public precision check
//! changed in #39 (8da300e: `add`, `sub`, `mul`, `div`). For every p in 1..=P_MAX, every operand representable at p
//! (both signs, every p-bit mantissa, exponents in a window, and signed zero) is paired with every other, and each op
//! is checked bit for bit in all 7 rounding modes against an exact integer oracle written here (no zenith-float
//! arithmetic in the oracle). A second pass shifts the operands to the top of the exponent range so results overflow,
//! and checks the overflow result: the largest finite value at exactly p bits when the mode rounds toward zero for the
//! result's sign (`ToZero`, `Down` positive, `Up` negative), else infinity (`None` included; CHANGELOG 2.0 overflow
//! entry, tests/overflow_modes.rs). The point count and the overflow count are printed.
//!
//! P_MAX is 8 in release (the claim: p in 1..=8) and 5 in debug (the debug build is about 30 times slower).
//! Underflow is not enumerated here: the bottom of the range has denormal results (exponent_range.rs covers them).

use std::ops::Neg;
use zenith_float_num::{ExactNum, RoundingMode, EXPONENT_MAX};

const MODES: [RoundingMode; 7] = [
    RoundingMode::None,
    RoundingMode::Up,
    RoundingMode::Down,
    RoundingMode::ToZero,
    RoundingMode::FromZero,
    RoundingMode::ToEven,
    RoundingMode::ToOdd,
];

/// Exponent window of the operands: value = m * 2^(e - p) with m a p-bit mantissa, e in E_LO..=E_HI.
const E_LO: i64 = -2;
const E_HI: i64 = 2;
/// Storage precision of every p <= 64 (one word): `None` keeps it (RNDZ at 64 bits).
const STORAGE: u32 = 64;

fn p_max() -> u32 {
    if cfg!(debug_assertions) {
        5
    } else {
        8
    }
}

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
enum Op {
    Add,
    Sub,
    Mul,
    Div,
}

/// An exact operand: (-1)^neg * n * 2^e (n = 0 is a signed zero).
#[derive(Copy, Clone, Debug)]
struct Val {
    neg: bool,
    n: u128,
    e: i64,
}

fn bitlen(n: u128) -> u32 {
    u128::BITS - n.leading_zeros()
}

fn operands(p: u32) -> Vec<Val> {
    let mut v = vec![
        Val {
            neg: false,
            n: 0,
            e: 0,
        },
        Val {
            neg: true,
            n: 0,
            e: 0,
        },
    ];
    let lo = 1u128 << (p - 1);
    let hi = 1u128 << p;
    for neg in [false, true] {
        for m in lo..hi {
            for e in E_LO..=E_HI {
                v.push(Val {
                    neg,
                    n: m,
                    e: e - i64::from(p),
                });
            }
        }
    }
    v
}

fn to_num(x: Val, shift: i64) -> Result<ExactNum, String> {
    let lo = u64::try_from(x.n).map_err(|e| format!("mantissa {e}"))?;
    let r = ExactNum::from_u64(lo, 64).ldexp(x.e + shift, 64, RoundingMode::ToZero);
    Ok(if x.neg { r.neg() } else { r })
}

/// A rounded oracle result.
#[derive(Debug, PartialEq)]
enum Want {
    Zero { neg: bool },
    Finite { neg: bool, n: u128, e: i64 },
    Overflow { neg: bool },
}

/// Round (-1)^neg * (n * 2^e + sticky tail below bit 0) to `bits` significant bits in `rm`.
fn round(
    neg: bool,
    n: u128,
    e: i64,
    sticky: bool,
    bits: u32,
    rm: RoundingMode,
) -> Result<(u128, i64), String> {
    let l = bitlen(n);
    if l <= bits {
        if sticky {
            return Err("oracle needs more quotient bits".to_string());
        }
        return Ok((n, e));
    }
    let s = l - bits;
    let kept = n >> s;
    let rem = n & ((1u128 << s) - 1);
    let half = 1u128 << (s - 1);
    let inexact = rem != 0 || sticky;
    let above_half = rem > half || (rem == half && sticky);
    let tie = rem == half && !sticky;
    let up = match rm {
        RoundingMode::ToZero | RoundingMode::None => false,
        RoundingMode::FromZero => inexact,
        RoundingMode::Up => inexact && !neg,
        RoundingMode::Down => inexact && neg,
        RoundingMode::ToEven => above_half || (tie && kept & 1 == 1),
        RoundingMode::ToOdd => above_half || (tie && kept & 1 == 0),
    };
    let k = if up { kept + 1 } else { kept };
    let e2 = e + i64::from(s);
    if bitlen(k) > bits {
        Ok((k >> 1, e2 + 1))
    } else {
        Ok((k, e2))
    }
}

/// Exact result of `op` as (neg, n, e, sticky); `None` for x/0.
fn exact(op: Op, a: Val, b: Val) -> Option<(bool, u128, i64, bool)> {
    match op {
        Op::Add | Op::Sub => {
            let bneg = if op == Op::Sub { !b.neg } else { b.neg };
            if a.n == 0 && b.n == 0 {
                return Some((a.neg && bneg, 0, 0, false));
            }
            if a.n == 0 {
                return Some((bneg, b.n, b.e, false));
            }
            if b.n == 0 {
                return Some((a.neg, a.n, a.e, false));
            }
            let e = a.e.min(b.e);
            let an = a.n << (a.e - e);
            let bn = b.n << (b.e - e);
            if a.neg == bneg {
                Some((a.neg, an + bn, e, false))
            } else if an >= bn {
                Some((a.neg, an - bn, e, false))
            } else {
                Some((bneg, bn - an, e, false))
            }
        }
        Op::Mul => Some((a.neg != b.neg, a.n * b.n, a.e + b.e, false)),
        Op::Div => {
            if b.n == 0 {
                return None;
            }
            let k = 80u32;
            let num = a.n << k;
            let q = num / b.n;
            let r = num % b.n;
            Some((a.neg != b.neg, q, a.e - b.e - i64::from(k), r != 0))
        }
    }
}

/// The oracle: exact result rounded at p (`None`: at the 64-bit storage, toward zero), `shift` added to the exponent.
fn want(
    op: Op,
    a: Val,
    b: Val,
    p: u32,
    rm: RoundingMode,
    shift: i64,
) -> Result<Option<Want>, String> {
    let Some((neg, n, e, sticky)) = exact(op, a, b) else {
        return Ok(None);
    };
    if n == 0 {
        // Exact zero: a difference of equal magnitudes is +0, -0 only in Down; otherwise the sign of the operands.
        let zneg = match op {
            Op::Add | Op::Sub if a.n != 0 || b.n != 0 => rm == RoundingMode::Down,
            Op::Add | Op::Sub => {
                if neg {
                    true
                } else {
                    let bneg = if op == Op::Sub { !b.neg } else { b.neg };
                    a.neg != bneg && rm == RoundingMode::Down
                }
            }
            Op::Mul | Op::Div => neg,
        };
        return Ok(Some(Want::Zero { neg: zneg }));
    }
    let bits = if rm == RoundingMode::None { STORAGE } else { p };
    let (k, e2) = round(neg, n, e, sticky, bits, rm)?;
    let top = i64::from(bitlen(k)) + e2 + shift;
    if top > EXPONENT_MAX {
        return Ok(Some(Want::Overflow { neg }));
    }
    Ok(Some(Want::Finite { neg, n: k, e: e2 }))
}

fn toward_zero(rm: RoundingMode, neg: bool) -> bool {
    match rm {
        RoundingMode::ToZero => true,
        RoundingMode::None => false,
        RoundingMode::Down => !neg,
        RoundingMode::Up => neg,
        RoundingMode::FromZero | RoundingMode::ToEven | RoundingMode::ToOdd => false,
    }
}

fn run(op: Op, a: &ExactNum, b: &ExactNum, p: usize, rm: RoundingMode) -> ExactNum {
    match op {
        Op::Add => ExactNum::add_ref(a, b, p, rm),
        Op::Sub => ExactNum::sub_ref(a, b, p, rm),
        Op::Mul => ExactNum::mul_ref(a, b, p, rm),
        Op::Div => ExactNum::div_ref(a, b, p, rm),
    }
}

/// Ok(None): x/0 (not counted); Ok(Some(overflowed)) on a match.
fn check(
    op: Op,
    a: Val,
    b: Val,
    p: u32,
    rm: RoundingMode,
    shift: (i64, i64),
) -> Result<Option<bool>, String> {
    let total = match op {
        Op::Add | Op::Sub => shift.0,
        Op::Mul => shift.0 + shift.1,
        Op::Div => shift.0 - shift.1,
    };
    let Some(w) = want(op, a, b, p, rm, total)? else {
        return Ok(None);
    };
    let overflowed = matches!(w, Want::Overflow { .. });
    let pu = usize::try_from(p).map_err(|e| format!("{e}"))?;
    let got = run(op, &to_num(a, shift.0)?, &to_num(b, shift.1)?, pu, rm);
    let info =
        || format!("{op:?} p={p} {rm:?} a={a:?} b={b:?} shift={shift:?}: want {w:?}, got {got:?}");
    let ok = match w {
        Want::Zero { neg } => got.is_zero() && got.is_negative() == neg,
        Want::Overflow { neg } => {
            if toward_zero(rm, neg) {
                // The largest finite value at exactly p bits: (2^p - 1) * 2^(EXPONENT_MAX - p).
                let n = (1u128 << p) - 1;
                let max = to_num(
                    Val {
                        neg,
                        n,
                        e: EXPONENT_MAX - i64::from(p),
                    },
                    0,
                )?;
                !got.is_inf() && got == max && got.is_negative() == neg
            } else {
                got.is_inf() && got.is_negative() == neg
            }
        }
        Want::Finite { neg, n, e } => {
            let exp = to_num(Val { neg, n, e }, total)?;
            !got.is_inf() && !got.is_nan() && got == exp && got.is_negative() == neg
        }
    };
    if ok {
        Ok(Some(overflowed))
    } else {
        Err(info())
    }
}

fn sweep(shift: (i64, i64), ops: &[Op], label: &str, want_overflow: bool) -> Result<(), String> {
    let mut points: u64 = 0;
    let mut overflows: u64 = 0;
    for p in 1..=p_max() {
        let v = operands(p);
        for &op in ops {
            for &a in &v {
                for &b in &v {
                    for rm in MODES {
                        if let Some(o) = check(op, a, b, p, rm, shift)? {
                            points += 1;
                            if o {
                                overflows += 1;
                            }
                        }
                    }
                }
            }
        }
    }
    println!(
        "small_p_exhaustive {label}: p 1..={} x exponents {E_LO}..={E_HI} x 7 modes: {points} points checked, {overflows} overflowed",
        p_max()
    );
    if points == 0 {
        return Err("no points".to_string());
    }
    if want_overflow && overflows == 0 {
        return Err(format!("{label}: no point overflowed"));
    }
    Ok(())
}

#[test]
fn add_sub_mul_div_every_small_p_operand_pair_in_7_modes() -> Result<(), String> {
    sweep(
        (0, 0),
        &[Op::Add, Op::Sub, Op::Mul, Op::Div],
        "core window",
        false,
    )
}

#[test]
fn add_sub_at_the_top_of_the_range_overflow_in_7_modes() -> Result<(), String> {
    // Operand exponents reach E_HI + s = EXPONENT_MAX; sums one binade up overflow.
    let s = EXPONENT_MAX - E_HI;
    sweep((s, s), &[Op::Add, Op::Sub], "add/sub at EXPONENT_MAX", true)
}

#[test]
fn mul_div_at_the_top_of_the_range_overflow_in_7_modes() -> Result<(), String> {
    // a at the top (exponent up to EXPONENT_MAX), b in the window: products and quotients above 1 overflow.
    sweep(
        (EXPONENT_MAX - E_HI, 0),
        &[Op::Mul, Op::Div],
        "mul/div at EXPONENT_MAX",
        true,
    )
}
