//! The exponent range is `i32::MIN..=i32::MAX` on every target (2.0.0 audit, defect 4).
//! Run on a 32-bit target (`scripts/ci_hex_32bit.sh`, `i686-unknown-linux-musl`) as well as
//! on the host: results near both limits must be the same on both.

use std::ops::Neg;
use zenith_float_num::{ExactNum, RoundingMode, EXPONENT_MAX, EXPONENT_MIN, WORD_BIT_SIZE};

const MODES: [RoundingMode; 6] = [
    RoundingMode::ToEven,
    RoundingMode::ToOdd,
    RoundingMode::Up,
    RoundingMode::Down,
    RoundingMode::ToZero,
    RoundingMode::FromZero,
];

fn precisions() -> [usize; 4] {
    [WORD_BIT_SIZE, 2 * WORD_BIT_SIZE, 3 * WORD_BIT_SIZE, 64 * WORD_BIT_SIZE]
}

fn signed(x: ExactNum, negative: bool) -> ExactNum {
    if negative {
        x.neg()
    } else {
        x
    }
}

#[test]
fn limits_are_i32_bounds() {
    assert_eq!(EXPONENT_MAX, i32::MAX);
    assert_eq!(EXPONENT_MIN, i32::MIN);
    for p in precisions() {
        assert_eq!(
            ExactNum::max_value(p).exponent(),
            Some(i32::MAX),
            "max_value({p})"
        );
        assert_eq!(
            ExactNum::min_positive(p).exponent(),
            Some(i32::MIN),
            "min_positive({p})"
        );
        assert_eq!(
            ExactNum::min_positive_normal(p).exponent(),
            Some(i32::MIN),
            "min_positive_normal({p})"
        );
    }
}

#[test]
fn ldexp_reaches_both_limits() {
    for p in precisions() {
        for rm in MODES {
            for negative in [false, true] {
                let one = || signed(ExactNum::from_u8(1, p), negative);
                // 1 = 0.5 * 2^1, so 1 * 2^(k - 1) has exponent k.
                for k in [EXPONENT_MAX - 2, EXPONENT_MAX - 1, EXPONENT_MAX] {
                    let r = one().ldexp(k - 1, p, rm);
                    assert_eq!(
                        r.exponent(),
                        Some(k),
                        "ldexp to {k}: p={p} {rm:?} neg={negative}"
                    );
                    assert_eq!(r.is_negative(), negative);
                }
                for k in [EXPONENT_MIN + 1, EXPONENT_MIN + 2, EXPONENT_MIN + 3] {
                    let r = one().ldexp(k - 1, p, rm);
                    assert_eq!(
                        r.exponent(),
                        Some(k),
                        "ldexp to {k}: p={p} {rm:?} neg={negative}"
                    );
                }
            }
        }
    }
}

#[test]
fn arithmetic_at_the_top_of_the_range() {
    for p in precisions() {
        for rm in MODES {
            for negative in [false, true] {
                let two = ExactNum::from_u8(2, p);
                let half = ExactNum::from_u8(1, p).ldexp(-1, p, rm);
                let big = signed(
                    ExactNum::from_u8(1, p).ldexp(EXPONENT_MAX - 2, p, rm),
                    negative,
                );
                let info = format!("p={p} {rm:?} neg={negative}");
                assert_eq!(
                    ExactNum::mul_ref(&big, &two, p, rm).exponent(),
                    Some(EXPONENT_MAX),
                    "{info} mul"
                );
                assert_eq!(
                    ExactNum::add_ref(&big, &big, p, rm).exponent(),
                    Some(EXPONENT_MAX),
                    "{info} add"
                );
                assert_eq!(
                    big.div(&half, p, rm).exponent(),
                    Some(EXPONENT_MAX),
                    "{info} div"
                );
                let max = signed(ExactNum::max_value(p), negative);
                assert_eq!(
                    ExactNum::mul_ref(&max, &half, p, rm).exponent(),
                    Some(EXPONENT_MAX - 1),
                    "{info} max/2"
                );
                assert!(
                    ExactNum::sub_ref(&max, &max, p, rm).is_zero(),
                    "{info} max-max"
                );
                let toward_zero = rm == RoundingMode::ToZero
                    || (rm == RoundingMode::Down && !negative)
                    || (rm == RoundingMode::Up && negative);
                for r in [
                    two.mul(&max, p, rm),
                    ExactNum::add_ref(&max, &max, p, rm),
                    ExactNum::div_ref(&max, &half, p, rm),
                ] {
                    if toward_zero {
                        assert_eq!(r, max, "{info} overflow toward zero");
                    } else {
                        assert!(
                            r.is_inf() && r.is_negative() == negative,
                            "{info} overflow: {r:?}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn arithmetic_at_the_bottom_of_the_range() {
    for p in precisions() {
        for rm in MODES {
            for negative in [false, true] {
                let two = ExactNum::from_u8(2, p);
                let info = format!("p={p} {rm:?} neg={negative}");
                let low = signed(
                    ExactNum::from_u8(1, p).ldexp(EXPONENT_MIN + 1, p, rm),
                    negative,
                );
                // low = 0.5 * 2^(EXPONENT_MIN + 2), so each halving lowers the exponent by one.
                let r = low.div(&two, p, rm);
                assert_eq!(r.exponent(), Some(EXPONENT_MIN + 1), "{info} div");
                let r = r.div(&two, p, rm);
                assert_eq!(r.exponent(), Some(EXPONENT_MIN), "{info} div to the limit");
                let tiny = signed(ExactNum::min_positive_normal(p), negative);
                let r = ExactNum::div_ref(&tiny, &two, p, rm);
                assert_eq!(r.exponent(), Some(EXPONENT_MIN), "{info} subnormal");
                assert!(!r.is_zero(), "{info} subnormal is not zero");
                assert_eq!(r.mul(&two, p, rm), tiny, "{info} subnormal * 2");
                let m = signed(ExactNum::min_positive(p), negative);
                assert_eq!(
                    ExactNum::add_ref(&m, &m, p, rm).exponent(),
                    Some(EXPONENT_MIN),
                    "{info} min+min"
                );
                assert_eq!(
                    ExactNum::mul_ref(&tiny, &tiny, p, rm).is_zero(),
                    !matches!(rm, RoundingMode::FromZero | RoundingMode::Up),
                    "{info} underflow of min_normal^2"
                );
            }
        }
    }
}
