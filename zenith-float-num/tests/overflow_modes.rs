//! Exponent overflow follows the rounding mode (2.0.0 audit, defect 1).

use std::ops::Neg;
use zenith_float_num::{ExactNum, RoundingMode, EXPONENT_MAX};

/// `-x` for an operand that stays borrowed: fresh `-1` times `x` at `x`'s own precision.
fn neg_of(x: &ExactNum) -> ExactNum {
    let p = x
        .mantissa_max_bit_len()
        .unwrap_or(zenith_float_num::WORD_BIT_SIZE);
    ExactNum::from_i8(-1, p).mul(x, p, RoundingMode::None)
}

fn with_exp(mut x: ExactNum, e: i32) -> ExactNum {
    x.set_exponent(e);
    x
}

const MODES: [RoundingMode; 7] = [
    RoundingMode::ToEven,
    RoundingMode::ToOdd,
    RoundingMode::Up,
    RoundingMode::Down,
    RoundingMode::ToZero,
    RoundingMode::FromZero,
    RoundingMode::None,
];

fn expect_overflow(r: &ExactNum, negative: bool, p: usize, rm: RoundingMode, what: &str) {
    let toward_zero = rm == RoundingMode::ToZero
        || (rm == RoundingMode::Down && !negative)
        || (rm == RoundingMode::Up && negative);

    if toward_zero {
        assert!(
            !r.is_inf(),
            "{what}: p={p} rm={rm:?} want max finite, got {r:?}"
        );
        let max = ExactNum::max_value(p);
        let expected = if negative { max.neg() } else { max };
        assert_eq!(*r, expected, "{what}: p={p} rm={rm:?} negative={negative}");
    } else {
        assert!(r.is_inf(), "{what}: p={p} rm={rm:?} want Inf, got {r:?}");
        assert_eq!(r.is_negative(), negative, "{what}: p={p} rm={rm:?} sign");
    }
}

#[test]
fn overflow_follows_rounding_mode() {
    for &p in &[64, 128, 256, 1024] {
        let max = ExactNum::max_value(p);
        let half = with_exp(ExactNum::from_i8(1, p), 0);
        let two = ExactNum::from_i8(2, p);
        let onehalf = with_exp(ExactNum::from_i8(3, p), 1);

        for &rm in &MODES {
            // max.add(&max) -> pos
            expect_overflow(
                &ExactNum::add_ref(&ExactNum::max_value(p), &max, p, rm),
                false,
                p,
                rm,
                "add pos",
            );

            // max.neg().add(&max.neg()) -> neg
            expect_overflow(
                &neg_of(&max).add(&neg_of(&max), p, rm),
                true,
                p,
                rm,
                "neg add neg",
            );

            // max.sub(&max.neg()) -> pos
            expect_overflow(
                &ExactNum::sub_ref(&max, &neg_of(&max), p, rm),
                false,
                p,
                rm,
                "sub pos",
            );

            // max.neg().sub(&max) -> neg
            expect_overflow(&neg_of(&max).sub(&max, p, rm), true, p, rm, "neg sub pos");

            // max.mul(&two) -> pos
            expect_overflow(&max.mul(&two, p, rm), false, p, rm, "mul pos by 2");

            // max.mul(&onehalf) -> pos
            expect_overflow(&max.mul(&onehalf, p, rm), false, p, rm, "mul pos by 1.5");

            // max.neg().mul(&two) -> neg
            expect_overflow(&neg_of(&max).mul(&two, p, rm), true, p, rm, "neg mul by 2");

            // max.mul(&two.neg()) -> neg
            expect_overflow(&max.mul(&neg_of(&two), p, rm), true, p, rm, "pos mul by -2");

            // max.neg().mul(&two.neg()) -> pos
            expect_overflow(
                &neg_of(&max).mul(&neg_of(&two), p, rm),
                false,
                p,
                rm,
                "neg mul by -2",
            );

            // max.div(&half) -> pos
            expect_overflow(&max.div(&half, p, rm), false, p, rm, "div pos by 0.5");

            // max.neg().div(&half) -> neg
            expect_overflow(
                &neg_of(&max).div(&half, p, rm),
                true,
                p,
                rm,
                "neg div by 0.5",
            );

            // max.div(&half.neg()) -> neg
            expect_overflow(
                &max.div(&neg_of(&half), p, rm),
                true,
                p,
                rm,
                "pos div by -0.5",
            );

            // big = from_i8(3,p) with set_exponent(EXPONENT_MAX / 2 + 2); big.mul(&big) -> pos
            let big = with_exp(ExactNum::from_i8(3, p), EXPONENT_MAX / 2 + 2);
            expect_overflow(&big.mul(&big, p, rm), false, p, rm, "large mul pos");

            // big.neg().mul(&big) -> neg
            expect_overflow(
                &neg_of(&big).mul(&big, p, rm),
                true,
                p,
                rm,
                "large neg mul pos",
            );
        }
    }
}

#[test]
fn finite_max_is_not_overflow() {
    for &p in &[64, 128, 256, 1024] {
        let max = ExactNum::max_value(p);
        let one = ExactNum::from_i8(1, p);
        let zero = ExactNum::from_i8(0, p);

        for &rm in &MODES {
            assert_eq!(max.mul(&one, p, rm), max, "max*1 p={p} rm={rm:?}");
            assert_eq!(
                neg_of(&max).add(&zero, p, rm),
                neg_of(&max),
                "-max+0 p={p} rm={rm:?}"
            );
        }
    }
}
