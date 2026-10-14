//! The inexact flag reports whether this operation rounded (2.0.0 audit, defect 3).

use zenith_float_num::{ExactNum, RoundingMode};

fn with_exp(mut x: ExactNum, e: i32) -> ExactNum {
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

                let third = one.div(&three, p, rm);
                assert!(
                    third.inexact(),
                    "Expected division to be inexact for p={}, rm={:?}, negative={}",
                    p,
                    rm,
                    negative
                );

                let tiny = with_exp(ExactNum::from_i8(1, p), -(p as i32) - 8);
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
                exact(&third.mul(&two, p, rm), "third * two", p, rm, negative);
                exact(&third.div(&two, p, rm), "third / two", p, rm, negative);
                exact(&third.mul(&half, p, rm), "third * half", p, rm, negative);
                exact(&third.add(&third, p, rm), "third + third", p, rm, negative);
                exact(&third.sub(&third, p, rm), "third - third", p, rm, negative);
                exact(
                    &third.sub(&third.neg(), p, rm),
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
                    &third.add(&ExactNum::from_i8(0, p), p, rm),
                    "third + 0",
                    p,
                    rm,
                    negative,
                );
                exact(
                    &third.mul(&ExactNum::from_i8(1, p), p, rm),
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
                exact(&two.mul(&three, p, rm), "two * three", p, rm, negative);
                exact(&three.sub(&two, p, rm), "three - two", p, rm, negative);
            }
        }
    }
}
