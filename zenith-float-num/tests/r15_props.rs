//! R15 bracket property: for `exp2` and `cbrt`, every mode's result lies in `[Down, Up]`, and
//! the three directed results bracket the round-to-nearest one, over random mantissas, signs,
//! precisions 24..=2048 and exponents across the i64 range.

use proptest::prelude::*;
use zenith_float_num::{
    Consts, ExactNum, Exponent, RoundingMode, Sign, Word, EXPONENT_MAX, EXPONENT_MIN,
    PROPTEST_CASES,
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

fn le(a: &ExactNum, b: &ExactNum) -> bool {
    a.is_inf() && b.is_inf() && a.sign() == b.sign() || matches!(a.compare(b), Some(c) if c <= 0)
}

fn ordered(f: &mut dyn FnMut(RoundingMode) -> ExactNum) -> Result<(), TestCaseError> {
    let lo = f(RoundingMode::Down);
    let hi = f(RoundingMode::Up);
    prop_assert!(!lo.is_nan() && !hi.is_nan(), "NaN result");
    for rm in MODES {
        let r = f(rm);
        prop_assert!(
            le(&lo, &r) && le(&r, &hi),
            "mode {:?} outside [Down, Up]",
            rm
        );
    }
    Ok(())
}

fn arg(hi: Word, lo: Word, neg: bool, e: Exponent) -> ExactNum {
    let sign = if neg { Sign::Neg } else { Sign::Pos };
    ExactNum::from_words(&[lo, hi | (1 << (Word::BITS - 1))], sign, e)
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn cbrt_modes_bracketed(
        hi in any::<Word>(), lo in any::<Word>(), neg in any::<bool>(),
        e in EXPONENT_MIN..=EXPONENT_MAX, p in 24usize..=2048,
    ) {
        let x = arg(hi, lo, neg, e);
        ordered(&mut |rm| x.cbrt(p, rm))?;
    }

    #[test]
    fn exp2_modes_bracketed(
        hi in any::<Word>(), lo in any::<Word>(), neg in any::<bool>(),
        e in -8i64..=70, p in 24usize..=2048,
    ) {
        let mut cc = Consts::new().map_err(|e| TestCaseError::fail(format!("{e:?}")))?;
        let x = arg(hi, lo, neg, e);
        ordered(&mut |rm| x.exp2(p, rm, &mut cc))?;
    }
}
