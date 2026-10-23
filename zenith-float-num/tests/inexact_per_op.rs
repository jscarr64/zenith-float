//! The inexact flag reports whether this operation rounded (2.0.0 audit, defect 3).

use std::ops::Neg;
use zenith_float_num::{ExactNum, RoundingMode};

/// `-x` for an operand that stays borrowed: fresh `-1` times `x` at `x`'s own precision.
fn neg_of(x: &ExactNum) -> ExactNum {
    let p = x
        .mantissa_max_bit_len()
        .unwrap_or(zenith_float_num::WORD_BIT_SIZE);
    ExactNum::from_i8(-1, p).mul(x, p, RoundingMode::None)
}

fn with_exp(mut x: ExactNum, e: zenith_float_num::Exponent) -> ExactNum {
    x.set_exponent(e);
    x
}

const MODES: [RoundingMode; 6] = [
    RoundingMode::ToEven,
    RoundingMode::ToOdd,
    RoundingMode::Up,
    RoundingMode::Down,
    RoundingMode::ToZero,
    RoundingMode::FromZero,
];

const PRECS: [usize; 4] = [64, 128, 256, 1024];

fn exact(x: &ExactNum, what: &str, p: usize, rm: RoundingMode, negative: bool) {
    assert!(
        !x.inexact(),
        "Expected {} to be exact for p={}, rm={:?}, negative={}",
        what,
        p,
        rm,
        negative
    );
}

#[test]
fn inexact_is_set_when_this_op_rounds() {
    for &p in &PRECS {
        for &rm in &MODES {
            for &negative in &[false, true] {
                let one = ExactNum::from_i8(1, p);
                let one = if negative { one.neg() } else { one };
                let three = ExactNum::from_i8(3, p);

                let third = ExactNum::div_ref(&one, &three, p, rm);
                assert_eq!(
                    third.is_negative(),
                    negative,
                    "the quotient of the signed one keeps its sign: p={p}, rm={rm:?}"
                );
                assert!(
                    third.inexact(),
                    "Expected division to be inexact for p={}, rm={:?}, negative={}",
                    p,
                    rm,
                    negative
                );

                let tiny = with_exp(ExactNum::from_i8(1, p), -(p as i64) - 8);
                let s = one.add(&tiny, p, rm);
                assert!(
                    s.inexact(),
                    "Expected addition to be inexact for p={}, rm={:?}, negative={}",
                    p,
                    rm,
                    negative
                );
            }
        }
    }
}

#[test]
fn inexact_is_not_inherited_from_operands() {
    for &p in &PRECS {
        for &rm in &MODES {
            for &negative in &[false, true] {
                let one = ExactNum::from_i8(1, p);
                let one = if negative { one.neg() } else { one };
                let three = ExactNum::from_i8(3, p);
                let two = ExactNum::from_i8(2, p);
                let half = with_exp(ExactNum::from_i8(1, p), 0);

                let third = one.div(&three, p, rm);
                assert!(
                    third.inexact(),
                    "Expected division to be inexact for p={}, rm={:?}, negative={}",
                    p,
                    rm,
                    negative
                );

                // These operations should be exact and not inexact
                exact(
                    &ExactNum::mul_ref(&third, &two, p, rm),
                    "third * two",
                    p,
                    rm,
                    negative,
                );
                exact(
                    &ExactNum::div_ref(&third, &two, p, rm),
                    "third / two",
                    p,
                    rm,
                    negative,
                );
                exact(
                    &ExactNum::mul_ref(&third, &half, p, rm),
                    "third * half",
                    p,
                    rm,
                    negative,
                );
                exact(
                    &ExactNum::add_ref(&third, &third, p, rm),
                    "third + third",
                    p,
                    rm,
                    negative,
                );
                exact(
                    &ExactNum::sub_ref(&third, &third, p, rm),
                    "third - third",
                    p,
                    rm,
                    negative,
                );
                exact(
                    &neg_of(&third).neg().add(&third, p, rm),
                    "third - (-third)",
                    p,
                    rm,
                    negative,
                );
                exact(
                    &ExactNum::from_i8(0, p).add(&third, p, rm),
                    "0 + third",
                    p,
                    rm,
                    negative,
                );
                exact(
                    &ExactNum::from_i8(0, p).add(&third, p, rm),
                    "third + 0",
                    p,
                    rm,
                    negative,
                );
                exact(
                    &ExactNum::from_i8(1, p).mul(&third, p, rm),
                    "third * 1",
                    p,
                    rm,
                    negative,
                );
                exact(
                    &third.div(&ExactNum::from_i8(-1, p), p, rm),
                    "third / (-1)",
                    p,
                    rm,
                    negative,
                );

                // Exact operations on exact operands
                exact(
                    &ExactNum::from_i8(2, p).mul(&three, p, rm),
                    "two * three",
                    p,
                    rm,
                    negative,
                );
                exact(&three.sub(&two, p, rm), "three - two", p, rm, negative);
            }
        }
    }
}

#[test]
fn accumulated_flag_still_inherits() {
    for &p in &PRECS {
        for &rm in &MODES {
            let third = ExactNum::from_i8(1, p).div(&ExactNum::from_i8(3, p), p, rm);
            let two = ExactNum::from_i8(2, p);
            let doubled = ExactNum::mul_ref(&third, &two, p, rm);
            assert!(!doubled.inexact(), "per-op p={p} rm={rm:?}");
            assert!(doubled.inexact_accumulated(), "accumulated p={p} rm={rm:?}");
            let back = doubled.sub(&third, p, rm);
            assert!(!back.inexact(), "per-op sub p={p} rm={rm:?}");
            assert!(
                back.inexact_accumulated(),
                "accumulated sub p={p} rm={rm:?}"
            );
            let exact = ExactNum::from_i8(2, p).add(&two, p, rm);
            assert!(!exact.inexact_accumulated(), "exact chain p={p} rm={rm:?}");
        }
    }
}

// ExactComplex: `inexact()` reports this operation, `inexact_accumulated()` the history
// (2.0.0 audit, round 2: item c).

fn cx(re: i8, im: i8, p: usize) -> zenith_float_num::ExactComplex {
    let mut z =
        zenith_float_num::ExactComplex::new(ExactNum::from_i8(re, p), ExactNum::from_i8(im, p));
    // operands that were rounded earlier: the flag must not leak into this operation
    z.set_inexact(true);
    z
}

#[test]
fn complex_inexact_is_per_operation() -> Result<(), Box<dyn std::error::Error>> {
    use zenith_float_num::{Consts, ExactComplex};
    let mut cc = Consts::new()?;
    for &p in &[64usize, 128, 256] {
        for &rm in &MODES {
            let z34 = cx(3, 4, p);
            let z12 = cx(1, 2, p);
            let checks: [(&str, ExactComplex); 10] = [
                ("add", ExactComplex::add_ref(&z34, &z12, p, rm)),
                ("sub", ExactComplex::sub_ref(&z34, &z12, p, rm)),
                ("mul", z34.mul(&z12, p, rm)),
                ("div", cx(-5, 10, p).div(&z12, p, rm)),
                ("sqrt", z34.sqrt(p, rm, &mut cc)),
                ("exp(0)", cx(0, 0, p).exp(p, rm, &mut cc)),
                ("ln(1)", cx(1, 0, p).ln(p, rm, &mut cc)),
                ("conj", z34.conj()),
                ("ldexp", z34.ldexp(3, p, rm)),
                ("fma", z34.fma(&z12, &z12, p, rm)),
            ];
            for (what, r) in checks.iter() {
                assert!(
                    !r.inexact(),
                    "{what} p={p} {rm:?}: exact result reported inexact {r:?}"
                );
            }
            let abs = z34.abs(p, rm);
            assert!(!abs.inexact(), "abs p={p} {rm:?}");
            let r = cx(2, 0, p).sqrt(p, rm, &mut cc);
            assert!(r.inexact(), "sqrt(2) p={p} {rm:?} must be inexact");
            let r = z34.div(&cx(3, 0, p), p, rm);
            assert!(r.inexact(), "(3+4i)/3 p={p} {rm:?} must be inexact");
            assert!(
                z34.add(&z12, p, rm).inexact_accumulated(),
                "accumulated flag inherits"
            );
        }
    }
    Ok(())
}
