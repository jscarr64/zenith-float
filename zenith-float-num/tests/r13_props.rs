//! R13 bracket properties: `Down <= true value <= Up` for div, reciprocal and hypot, checked
//! with exact products (no oracle), plus hypot's special values and directional overflow.

use proptest::prelude::*;
use std::ops::Neg;
use zenith_float_num::{
    Error, ExactNum, Exponent, RoundingMode, Sign, INF_NEG, INF_POS, PROPTEST_CASES,
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
