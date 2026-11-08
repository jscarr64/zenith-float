//! The exponent range is `-(2^62 - 1)..=2^62 - 1` on every target (2.0.0 audit, defect 4; widened from `i32`).
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
fn limits_are_62_bit_bounds() {
    assert_eq!(EXPONENT_MAX, (1i64 << 62) - 1);
    assert_eq!(EXPONENT_MIN, -((1i64 << 62) - 1));
    for p in precisions() {
        assert_eq!(
            ExactNum::max_value(p).exponent(),
            Some(EXPONENT_MAX),
            "max_value({p})"
        );
        assert_eq!(
            ExactNum::min_positive(p).exponent(),
            Some(EXPONENT_MIN),
            "min_positive({p})"
        );
        assert_eq!(
            ExactNum::min_positive_normal(p).exponent(),
            Some(EXPONENT_MIN),
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
                    tiny.sqr(p, rm).is_zero(),
                    !matches!(rm, RoundingMode::FromZero | RoundingMode::Up),
                    "{info} underflow of min_normal^2"
                );
            }
        }
    }
}

/// `rem` of the largest finite number by the smallest subnormal is exactly zero.
///
/// At precision `p` the mantissa represents `0.m * 2^e`, with the last bit worth
/// `2^(e - p)`. `max_value(p) = (2^p - 1) * 2^(EXPONENT_MAX - p)` and
/// `min_positive(p) = 2^(EXPONENT_MIN - p)` (only the last bit set, exponent
/// `EXPONENT_MIN`). Their ratio is `(2^p - 1) * 2^(EXPONENT_MAX - EXPONENT_MIN)`,
/// an integer, so `|x| mod |y| = 0` exactly. The result carries the sign of `self`
/// (`+`), so it is `+0`, not NaN, and exact. Negating the dividend yields `-0`.
///
/// Before the fix a debug build panics at `zenith-float-num/src/num.rs:670`
/// ("attempt to subtract with overflow", `e1eff - e2eff` in `i64`), and a release
/// build returned a nonzero value with exponent `EXPONENT_MIN`.
#[test]
fn rem_of_max_value_by_min_positive_is_exact_zero() {
    for p in [WORD_BIT_SIZE, 2 * WORD_BIT_SIZE, 64 * WORD_BIT_SIZE] {
        let x = ExactNum::max_value(p);
        let y = ExactNum::min_positive(p);
        assert_eq!(x.exponent(), Some(EXPONENT_MAX), "p = {p}");
        assert_eq!(y.exponent(), Some(EXPONENT_MIN), "p = {p}");

        let r = x.rem(&y);
        assert!(!r.is_nan(), "p = {p}: rem returned NaN");
        assert!(
            r.is_zero(),
            "p = {p}: rem returned {:?}, expected exact +0",
            r.exponent()
        );
        assert!(!r.is_negative(), "p = {p}: rem of +max must be +0");
        assert!(!r.inexact(), "p = {p}: rem must be exact");

        let r_neg = x.neg().rem(&y);
        assert!(!r_neg.is_nan(), "p = {p}: rem of -max returned NaN");
        assert!(
            r_neg.is_zero(),
            "p = {p}: rem of -max returned {:?}, expected exact -0",
            r_neg.exponent()
        );
        assert!(r_neg.is_negative(), "p = {p}: rem of -max must be -0");
        assert!(!r_neg.inexact(), "p = {p}: rem of -max must be exact");
    }
}

/// Independent u128 oracle for rem (used by the edge tests below).
///
/// 2^d mod m by square-and-multiply. The products stay below 2^128
/// because every modulus here fits in u64.
fn powmod_u128(mut base: u128, mut exp: i128, modulus: u128) -> u128 {
    if modulus == 1 {
        return 0;
    }
    let mut result: u128 = 1;
    base %= modulus;
    while exp > 0 {
        if exp & 1 == 1 {
            result = result * base % modulus;
        }
        exp >>= 1;
        if exp > 0 {
            base = base * base % modulus;
        }
    }
    result
}

/// Absolute remainder for x = mx * 2^ex, y = my * 2^ey, as
/// (mantissa, effective exponent) of the exact remainder. No 2^spread
/// integer is ever built: the wide case reduces 2^(ex - ey) modulo my,
/// and the narrow case reduces mx modulo my * 2^(ey - ex). Callers keep
/// the narrow-case shift small (the exhaustive grid below uses at most 6);
/// the wide case is the one that needs the exact exponent gap.
fn rem_oracle_u128(mx: u128, ex: i128, my: u128, ey: i128) -> (u128, i128) {
    if ex >= ey {
        let r = ((mx % my) * powmod_u128(2, ex - ey, my)) % my;
        (r, ey)
    } else {
        let mut scaled = my;
        let mut s = ey - ex;
        while s > 0 {
            scaled *= 2;
            s -= 1;
        }
        (mx % scaled, ex)
    }
}

/// Edge cases for rem at extreme exponent spreads (PR #33).
///
/// A. rem_extreme_spread_odd_divisor_is_exact: the dividend is the largest
///    finite number and the divisor is 7 at the bottom of the range, so the
///    remainder is nonzero and depends on 2^d mod 7 with
///    d = EXPONENT_MAX - EXPONENT_MIN. A saturated or wrapped exponent gap
///    gives a different remainder, not just a panic.
/// B. rem_boundary_cases: equal operands, dividend below divisor, negative
///    divisor, and the smallest normal by the smallest subnormal, each with
///    the sign the result must carry.
/// C. rem_small_precision_exhaustive_matches_u128_oracle: every mantissa
///    pair at precisions 1..=5 over a small exponent grid, against the
///    independent u128 remainder above.
#[test]
fn rem_extreme_spread_odd_divisor_is_exact() {
    for p in [WORD_BIT_SIZE, 2 * WORD_BIT_SIZE] {
        let x = ExactNum::max_value(p);
        let y = ExactNum::from_u8(7, p).mul(&ExactNum::min_positive(p), p, RoundingMode::None);
        // p is 64 or 128 here, so this conversion cannot fail; the default
        // only keeps the binding total without an unwrap.
        let mut p_sub_1 = u64::MAX;
        if let Ok(v) = u64::try_from(p - 1) {
            p_sub_1 = v;
        }
        // Values are mantissa * 2^(e - max_bit_len), so max_value(p) is
        // (2^p - 1) * 2^(EXPONENT_MAX - p) and the divisor is
        // 7 * 2^(EXPONENT_MIN - p). The gap is EXPONENT_MAX - EXPONENT_MIN.
        let p128 = i128::from(p_sub_1) + 1;
        let ex = i128::from(EXPONENT_MAX) - p128;
        let ey = i128::from(EXPONENT_MIN) - p128;
        let (expected_mant, expected_exp) = rem_oracle_u128(u128::MAX >> (128 - p), ex, 7, ey);
        assert_eq!(expected_exp, ey, "p = {p}");
        let expected = match u8::try_from(expected_mant) {
            Ok(v) => ExactNum::from_u8(v, p).mul(&ExactNum::min_positive(p), p, RoundingMode::None),
            Err(_) => ExactNum::from_u8(0, p),
        };
        let r = x.rem(&y);
        assert_eq!(r, expected, "p = {p}: rem of max by 7 * min_positive");
        assert!(!r.inexact(), "p = {p}: rem must be exact");

        let r_neg = x.neg().rem(&y);
        assert_eq!(r_neg, expected.neg(), "p = {p}: rem of -max");
        assert!(!r_neg.inexact(), "p = {p}: rem of -max must be exact");
    }
}

#[test]
fn rem_boundary_cases() {
    for p in [WORD_BIT_SIZE, 2 * WORD_BIT_SIZE] {
        // Equal operands: exact +0.
        let x = ExactNum::max_value(p);
        let r = x.rem(&x);
        assert!(r.is_zero(), "p = {p}: max rem max must be +0");
        assert!(!r.is_negative(), "p = {p}: max rem max must be +0");
        assert!(!r.inexact(), "p = {p}: rem must be exact");

        // Dividend below divisor: the dividend itself.
        let small = ExactNum::min_positive(p);
        let big = ExactNum::min_positive_normal(p);
        let r_below = small.rem(&big);
        assert_eq!(r_below, small, "p = {p}: min_positive rem min_normal");

        // Negative divisor: the sign follows the dividend, and the
        // extreme spread still has to survive (red before the fix).
        let y = ExactNum::min_positive(p);
        let r_negdiv = x.rem(&y.neg());
        assert!(r_negdiv.is_zero(), "p = {p}: max rem -min must be +0");
        assert!(!r_negdiv.is_negative(), "p = {p}: sign follows dividend");

        // Smallest normal by smallest subnormal: an exact power of two.
        let r_norm = big.rem(&small);
        assert!(r_norm.is_zero(), "p = {p}: min_normal rem min_positive");
        assert!(!r_norm.is_negative(), "p = {p}: must be +0");
    }
}

#[test]
fn rem_small_precision_exhaustive_matches_u128_oracle() {
    // Operands are built at word precision: the integer constructors reject
    // p below the integer width, so the small grid is over mantissa widths
    // (the values are exact at any precision that holds them).
    for p in 1..=5 {
        for m1 in 1..(1u64 << p) {
            for m2 in 1..(1u64 << p) {
                for e1 in -3i64..=3 {
                    for e2 in -3i64..=3 {
                        let x = ExactNum::from_u64(m1, WORD_BIT_SIZE).ldexp(
                            e1,
                            WORD_BIT_SIZE,
                            RoundingMode::None,
                        );
                        let y = ExactNum::from_u64(m2, WORD_BIT_SIZE).ldexp(
                            e2,
                            WORD_BIT_SIZE,
                            RoundingMode::None,
                        );
                        let (rm, re) = rem_oracle_u128(
                            u128::from(m1),
                            i128::from(e1),
                            u128::from(m2),
                            i128::from(e2),
                        );
                        let expected = match u64::try_from(rm) {
                            Ok(m) => match i64::try_from(re) {
                                Ok(e) => ExactNum::from_u64(m, WORD_BIT_SIZE).ldexp(
                                    e,
                                    WORD_BIT_SIZE,
                                    RoundingMode::None,
                                ),
                                Err(_) => ExactNum::from_u64(0, WORD_BIT_SIZE),
                            },
                            Err(_) => ExactNum::from_u64(0, WORD_BIT_SIZE),
                        };
                        assert_eq!(x.rem(&y), expected, "p={p} m1={m1} m2={m2} e1={e1} e2={e2}");
                        if m1 == 1 {
                            assert_eq!(
                                x.neg().rem(&y),
                                expected.neg(),
                                "p={p} m1={m1} m2={m2} e1={e1} e2={e2} neg"
                            );
                        }
                    }
                }
            }
        }
    }
}
