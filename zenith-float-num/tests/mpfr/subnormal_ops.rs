//! `add`, `sub`, `mul`, `div` and `sqr` with subnormal operands (either side, either sign) are correctly
//! rounded in every directed and nearest mode: checked against MPFR on a grid of subnormal and tiny normal
//! operands, through the borrowed forms (`*_ref`) and the owned forms. The kernels take a subnormal operand
//! as stored instead of normalizing a copy of it (`num.rs`), and this test is the outside oracle for that.

use crate::mpfr::common::{assert_exact_eq, mpfr_reference, BinOp};
use crate::mpfr::retry_budget::use_full_exponent_range;
use rug::Float;
use std::ops::Neg;
use zenith_float_num::{Consts, ExactNum, Exponent, RoundingMode, EXPONENT_MIN};

const MODES: [RoundingMode; 5] = [
    RoundingMode::ToEven,
    RoundingMode::Up,
    RoundingMode::Down,
    RoundingMode::ToZero,
    RoundingMode::FromZero,
];

/// Result precisions (whole words: the library rounds `p` up to a word) and operand precisions.
const RESULT_PRECISIONS: [usize; 3] = [64, 128, 192];
const OPERAND_PRECISION: usize = 192;

/// `(significand, exponent)` pairs: the operand is `significand * 2^exponent`. The exponents straddle the
/// subnormal boundary (`2^(EXPONENT_MIN - 1)` is the smallest normal number); the last three are ordinary.
fn operands() -> Vec<(u64, i64)> {
    let min = EXPONENT_MIN as i64;
    let mut v = Vec::new();
    for m in [1u64, 3, 0xdead_beef_1234_5677, u64::MAX] {
        for k in [min - 120, min - 70, min - 21, min - 3, min, min + 5] {
            v.push((m, k));
        }
    }
    v.push((1, 1));
    v.push((5, -7));
    v.push((0xffff_ffff_ffff_fff1, 3));
    v
}

/// `k` as two `i32` shifts (an exponent below `EXPONENT_MIN` does not fit one).
fn split_exponent(k: i64) -> (Exponent, Exponent) {
    let hi = k.max(EXPONENT_MIN as i64);
    ((hi as Exponent), (k - hi) as Exponent)
}

/// MPFR's exponent of `m * 2^k`: it must lie inside the emulated range (`emin = EXPONENT_MIN - p + 1`
/// in `mpfr_reference`) for MPFR's answer to be the reference.
fn in_mpfr_range(m: u64, k: i64, p: usize) -> bool {
    k + 64 - m.leading_zeros() as i64 > EXPONENT_MIN as i64 - p as i64 + 1
}

fn build(m: u64, k: i64, negative: bool) -> ExactNum {
    let p = OPERAND_PRECISION;
    let (hi, lo) = split_exponent(k);
    let x = ExactNum::from_word(m, p)
        .ldexp(hi, p, RoundingMode::ToEven)
        .ldexp(lo, p, RoundingMode::ToEven);
    if negative {
        x.neg()
    } else {
        x
    }
}

fn mpfr_operand(m: u64, k: i64, negative: bool) -> Float {
    let (hi, lo) = split_exponent(k);
    let f = (Float::with_val(OPERAND_PRECISION as u32, m) << hi) << lo;
    if negative {
        -f
    } else {
        f
    }
}

#[test]
fn subnormal_operands_match_mpfr_in_every_mode() -> Result<(), Box<dyn std::error::Error>> {
    use_full_exponent_range(); // the operands are below MPFR's default `emin`
    let mut cc = Consts::new()?;
    let ops = operands();
    let mut checked = 0usize;
    for &(m1, k1) in &ops {
        for &(m2, k2) in &ops {
            for (n1, n2) in [(false, false), (true, false), (false, true)] {
                let f1 = mpfr_operand(m1, k1, n1);
                let f2 = mpfr_operand(m2, k2, n2);
                let a = build(m1, k1, n1);
                let b = build(m2, k2, n2);
                for p in RESULT_PRECISIONS {
                    if !in_mpfr_range(m1, k1, p) || !in_mpfr_range(m2, k2, p) {
                        continue;
                    }
                    for rm in MODES {
                        for (op, name) in [
                            (BinOp::Add, "add"),
                            (BinOp::Sub, "sub"),
                            (BinOp::Mul, "mul"),
                            (BinOp::Div, "div"),
                        ] {
                            let want = mpfr_reference(&f1, &f2, op, p as u32, rm)?;
                            let label = format!(
                                "{name} ({m1:#x}, {k1}, neg {n1}) ({m2:#x}, {k2}, neg {n2}) p={p} {rm:?}"
                            );
                            let by_ref = match op {
                                BinOp::Add => ExactNum::add_ref(&a, &b, p, rm),
                                BinOp::Sub => ExactNum::sub_ref(&a, &b, p, rm),
                                BinOp::Mul => ExactNum::mul_ref(&a, &b, p, rm),
                                BinOp::Div => ExactNum::div_ref(&a, &b, p, rm),
                            };
                            assert_exact_eq(&by_ref, &want, p, &format!("{label} (ref)"), &mut cc)?;
                            // The owned forms consume the first operand: rebuild it.
                            let first = build(m1, k1, n1);
                            let owned = match op {
                                BinOp::Add => first.add(&b, p, rm),
                                BinOp::Sub => first.sub(&b, p, rm),
                                BinOp::Mul => first.mul(&b, p, rm),
                                BinOp::Div => first.div(&b, p, rm),
                            };
                            assert_exact_eq(
                                &owned,
                                &want,
                                p,
                                &format!("{label} (owned)"),
                                &mut cc,
                            )?;
                            checked += 2;
                        }
                    }
                }
            }
        }
    }
    assert!(checked > 50_000, "grid too small: {checked}");
    Ok(())
}

#[test]
fn subnormal_square_matches_mpfr_in_every_mode() -> Result<(), Box<dyn std::error::Error>> {
    use_full_exponent_range();
    let mut cc = Consts::new()?;
    for &(m, k) in &operands() {
        let f = mpfr_operand(m, k, false);
        let x = build(m, k, false);
        for p in RESULT_PRECISIONS {
            if !in_mpfr_range(m, k, p) {
                continue;
            }
            for rm in MODES {
                let want = mpfr_reference(&f, &f, BinOp::Mul, p as u32, rm)?;
                let got = x.sqr(p, rm);
                assert_exact_eq(
                    &got,
                    &want,
                    p,
                    &format!("sqr ({m:#x}, {k}) p={p} {rm:?}"),
                    &mut cc,
                )?;
            }
        }
    }
    Ok(())
}
