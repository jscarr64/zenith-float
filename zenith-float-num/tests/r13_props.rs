//! R13 bracket properties: `Down <= true value <= Up` for div, reciprocal and hypot, checked
//! with exact products (no oracle), plus hypot's special values and directional overflow.

use proptest::prelude::*;
use std::ops::Neg;
use zenith_float_num::{
    Error, ExactNum, Exponent, RoundingMode, Sign, Word, EXPONENT_MAX, EXPONENT_MIN, INF_NEG,
    INF_POS, PROPTEST_CASES, WORD_BIT_SIZE,
};

const MODES: [RoundingMode; 5] = [
    RoundingMode::ToEven,
    RoundingMode::Up,
    RoundingMode::Down,
    RoundingMode::ToZero,
    RoundingMode::FromZero,
];

fn cfg() -> ProptestConfig {
    ProptestConfig::with_cases(PROPTEST_CASES)
}

/// `n * 2^shift` at `p` bits (exact: `n` has at most 63 bits and `p >= 64`).
fn num(n: i64, shift: Exponent, p: usize) -> ExactNum {
    let mut x = ExactNum::from_i64(n, p);
    if let Some(e) = x.exponent() {
        x.set_exponent(e + shift);
    }
    x
}

fn le(a: &ExactNum, b: &ExactNum) -> bool {
    matches!(a.compare(b), Some(c) if c <= 0)
}

/// `a * b` exactly (`a`, `b` at most `p` bits each).
fn prod(a: &ExactNum, b: &ExactNum, p: usize) -> Result<ExactNum, TestCaseError> {
    let r = ExactNum::from_i64(1, p)
        .mul(a, 2 * p, RoundingMode::None)
        .mul(b, 4 * p, RoundingMode::None);
    prop_assert!(!r.is_nan() && !r.inexact(), "product not exact");
    Ok(r)
}

/// Each mode's result lies between Down and Up; Down <= Up.
fn ordered(f: &dyn Fn(RoundingMode) -> ExactNum) -> Result<(ExactNum, ExactNum), TestCaseError> {
    let lo = f(RoundingMode::Down);
    let hi = f(RoundingMode::Up);
    for rm in MODES {
        let r = f(rm);
        prop_assert!(
            le(&lo, &r) && le(&r, &hi),
            "mode {:?} outside [Down, Up]",
            rm
        );
    }
    Ok((lo, hi))
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn hypot_brackets(a in any::<i32>(), b in any::<i32>(), ea in -300i64..300, eb in -300i64..300, w in 1usize..4) {
        let p = 64 * w;
        let x = num(i64::from(a), ea, p);
        let y = num(i64::from(b), eb, p);
        let (lo, hi) = ordered(&|rm| x.hypot(&y, p, rm))?;
        // x^2 + y^2 exactly: the exponents differ by under 700 bits.
        let wide = 4 * p + 2048;
        let s = prod(&x, &x, p)?.add(&prod(&y, &y, p)?, wide, RoundingMode::None);
        prop_assert!(!s.inexact());
        prop_assert!(le(&prod(&lo, &lo, p)?, &s) && le(&s, &prod(&hi, &hi, p)?));
    }

    #[test]
    fn div_brackets(a in any::<i64>(), b in 1i64..i64::MAX, ea in -300i64..300, eb in -300i64..300, w in 1usize..4) {
        let p = 64 * w;
        let x = num(a, ea, p);
        let d = num(b, eb, p);
        let (lo, hi) = ordered(&|rm| num(a, ea, p).div(&d, p, rm))?;
        prop_assert!(le(&prod(&lo, &d, p)?, &x) && le(&x, &prod(&hi, &d, p)?));
    }

    #[test]
    fn reciprocal_brackets(a in 1i64..i64::MAX, ea in -300i64..300, w in 1usize..4) {
        let p = 64 * w;
        let x = num(a, ea, p);
        let one = ExactNum::from_i64(1, p);
        let (lo, hi) = ordered(&|rm| x.reciprocal(p, rm))?;
        prop_assert!(le(&prod(&lo, &x, p)?, &one) && le(&one, &prod(&hi, &x, p)?));
    }
}

/// Exponents over the whole `i64` range, weighted toward both ends of the representable range
/// `±(2^62 − 1)` and just past them.
fn full_exp() -> impl Strategy<Value = Exponent> {
    prop_oneof![
        any::<i64>(),
        (EXPONENT_MIN - 70)..(EXPONENT_MIN + 300),
        (EXPONENT_MAX - 300)..(EXPONENT_MAX + 70),
        -300i64..300,
    ]
}

fn in_range(e: Exponent) -> bool {
    (EXPONENT_MIN..=EXPONENT_MAX).contains(&e)
}

/// Outside `±(2^62 − 1)` no number exists: building one gives NaN.
fn unrepresentable(e: Exponent) -> Result<(), TestCaseError> {
    let top: [Word; 1] = [1 << (WORD_BIT_SIZE - 1)];
    prop_assert!(
        ExactNum::from_words(&top, Sign::Pos, e).is_nan(),
        "e={} built a number",
        e
    );
    Ok(())
}

/// `n` (non-zero) with exponent set to `e` (`0.m · 2^e`); `e` must be in range.
fn at_exp(n: i64, e: Exponent, p: usize) -> ExactNum {
    let mut x = ExactNum::from_i64(n, p);
    x.set_exponent(e);
    x
}

/// `r` must equal `r0` scaled by `2^shift` whenever that is a normal in-range value (scaling is
/// exact), and lie in [Down, Up] in every case (overflow / underflow included).
fn scaled_matches(r: &ExactNum, mut r0: ExactNum, shift: i64) -> Result<(), TestCaseError> {
    let Some(e0) = r0.exponent() else {
        return Ok(());
    };
    let Some(e) = e0.checked_add(shift) else {
        return Ok(());
    };
    if r0.is_zero() || !in_range(e) {
        return Ok(());
    }
    r0.set_exponent(e);
    prop_assert!(
        r.compare(&r0) == Some(0),
        "scaled result differs at e={}",
        e
    );
    Ok(())
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn div_full_exponent_range(a in 1i64..i64::MAX, na in any::<bool>(), b in 1i64..i64::MAX, ex in full_exp(), ed in full_exp(), w in 1usize..4) {
        let p = 64 * w;
        for e in [ex, ed] {
            if !in_range(e) {
                return unrepresentable(e);
            }
        }
        let a = if na { -a } else { a };
        let d = at_exp(b, ed, p);
        ordered(&|rm| at_exp(a, ex, p).div(&d, p, rm))?;
        for rm in MODES {
            let r = at_exp(a, ex, p).div(&d, p, rm);
            let r0 = at_exp(a, 0, p).div(&at_exp(b, 0, p), p, rm);
            scaled_matches(&r, r0, ex - ed)?;
        }
    }

    #[test]
    fn reciprocal_full_exponent_range(a in 1i64..i64::MAX, ex in full_exp(), w in 1usize..4) {
        let p = 64 * w;
        if !in_range(ex) {
            return unrepresentable(ex);
        }
        let x = at_exp(a, ex, p);
        ordered(&|rm| x.reciprocal(p, rm))?;
        for rm in MODES {
            let r0 = at_exp(a, 0, p).reciprocal(p, rm);
            scaled_matches(&x.reciprocal(p, rm), r0, -ex)?;
        }
    }

    #[test]
    fn hypot_full_exponent_range(a in 1i32..i32::MAX, b in 1i32..i32::MAX, e in full_exp(), da in -64i64..64, db in -64i64..64, w in 1usize..4) {
        let p = 64 * w;
        let (Some(ea), Some(eb)) = (e.checked_add(da), e.checked_add(db)) else { return Ok(()) };
        for e in [ea, eb] {
            if !in_range(e) {
                return unrepresentable(e);
            }
        }
        let x = at_exp(i64::from(a), ea, p);
        let y = at_exp(i64::from(b), eb, p);
        ordered(&|rm| x.hypot(&y, p, rm))?;
        for rm in MODES {
            let r0 = at_exp(i64::from(a), da, p).hypot(&at_exp(i64::from(b), db, p), p, rm);
            scaled_matches(&x.hypot(&y, p, rm), r0, e)?;
        }
    }
}

fn pos_inf(r: &ExactNum) -> bool {
    r.is_inf() && r.sign() == Some(Sign::Pos)
}

#[test]
fn hypot_specials() {
    let p = 128;
    let nan = ExactNum::nan(Some(Error::InvalidArgument));
    let one = ExactNum::from_i64(1, p);
    for rm in MODES {
        for inf in [INF_POS, INF_NEG] {
            assert!(pos_inf(&inf.hypot(&nan, p, rm)), "hypot(Inf, NaN)");
            assert!(pos_inf(&nan.hypot(&inf, p, rm)), "hypot(NaN, Inf)");
            assert!(pos_inf(&inf.hypot(&one, p, rm)), "hypot(Inf, 1)");
        }
        assert!(nan.hypot(&one, p, rm).is_nan(), "hypot(NaN, 1)");
        let nz = ExactNum::from_i64(0, p).neg();
        let z = nz.hypot(&nz, p, rm);
        assert!(
            z.is_zero() && z.sign() == Some(Sign::Pos),
            "hypot(-0, -0) = +0"
        );
        let five = ExactNum::from_i64(3, p)
            .neg()
            .hypot(&ExactNum::from_i64(4, p), p, rm);
        assert_eq!(five.compare(&ExactNum::from_i64(5, p)), Some(0));
        assert!(
            !five.inexact() && !five.inexact_accumulated(),
            "3,4,5 is exact"
        );
    }
}

#[test]
fn hypot_overflow_is_directional() {
    for p in [64usize, 192, 1024] {
        let max = ExactNum::max_value(p);
        for rm in MODES {
            let r = max.hypot(&max, p, rm);
            match rm {
                RoundingMode::Down | RoundingMode::ToZero => {
                    assert_eq!(
                        r.compare(&ExactNum::max_value(p)),
                        Some(0),
                        "{rm:?} gives MAX"
                    )
                }
                _ => assert!(pos_inf(&r), "{rm:?} gives +Inf"),
            }
        }
    }
}
