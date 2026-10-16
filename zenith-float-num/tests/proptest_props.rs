//! Property tests (`PROPTEST_CASES`). No `mpfr-tests` feature.

use proptest::prelude::*;
use std::ops::Neg;
use zenith_float_num::{
    Consts, ExactNum, ExactNumArray, ExactRational, RoundingMode, PROPTEST_CASES, WORD_BIT_SIZE,
};

fn cfg() -> ProptestConfig {
    ProptestConfig::with_cases(PROPTEST_CASES)
}

fn eq(a: &ExactNum, b: &ExactNum) -> bool {
    a.compare(b) == Some(0)
}

fn array_eq(a: &ExactNumArray, b: &ExactNumArray) -> bool {
    a.shape() == b.shape() && a.as_slice().iter().zip(b.as_slice()).all(|(x, y)| eq(x, y))
}

fn rat_eq(a: &ExactRational, b: &ExactRational) -> bool {
    !a.is_nan() && !b.is_nan() && eq(a.num(), b.num()) && eq(a.den(), b.den())
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn add_commutes(n1 in -10_000i64..10_000, n2 in -10_000i64..10_000, w in 1usize..4) {
        // `from_i64` needs p >= 64 (documented: smaller precisions give NaN(InvalidArgument)),
        // so step in 64-bit units rather than words (a 32-bit word would give p = 32 at w = 1).
        let p = w * 64;
        let rm = RoundingMode::ToEven;
        let a = ExactNum::from_i64(n1, p);
        let b = ExactNum::from_i64(n2, p);
        prop_assert!(eq(&ExactNum::from_i64(n1, p).add(&b, p, rm), &b.add(&a, p, rm)));
    }

    /// Directed rounding only. `ToEven` double-rounds and is not this identity.
    #[test]
    fn round_then_coarser_eq_direct(
        n in -1_000_000i64..1_000_000,
        mode in 0u8..3,
    ) {
        let rm = match mode {
            0 => RoundingMode::ToZero,
            1 => RoundingMode::Down,
            _ => RoundingMode::Up,
        };
        let q = WORD_BIT_SIZE;
        let p = 2 * WORD_BIT_SIZE;
        let p_hi = 4 * WORD_BIT_SIZE;
        let mut via = ExactNum::from_i64(n, p_hi);
        via.set_precision(p, rm)?;
        via.set_precision(q, rm)?;
        let mut direct = ExactNum::from_i64(n, p_hi);
        direct.set_precision(q, rm)?;
        prop_assert!(eq(&via, &direct));
    }

    #[test]
    fn erf_is_odd(n in -8i32..9) {
        let p = 2 * WORD_BIT_SIZE;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let x = ExactNum::from_i32(n, p);
        let e = x.erf(p, rm, &mut cc);
        let en = x.neg().erf(p, rm, &mut cc);
        prop_assert!(eq(&e, &en.neg()));
    }

    #[test]
    fn matmul_associative_small_int(
        ents in prop::collection::vec(-2i64..=2, 12),
    ) {
        let p = 2 * WORD_BIT_SIZE;
        let cell = |v: i64| ExactNum::from_i64(v, p);
        let a = ExactNumArray::from_shape(p, 2, 2, &[cell(ents[0]), cell(ents[1]), cell(ents[2]), cell(ents[3])])
            .ok_or_else(|| proptest::test_runner::TestCaseError::fail("shape"))?;
        let b = ExactNumArray::from_shape(p, 2, 2, &[cell(ents[4]), cell(ents[5]), cell(ents[6]), cell(ents[7])])
            .ok_or_else(|| proptest::test_runner::TestCaseError::fail("shape"))?;
        let c = ExactNumArray::from_shape(p, 2, 2, &[cell(ents[8]), cell(ents[9]), cell(ents[10]), cell(ents[11])])
            .ok_or_else(|| proptest::test_runner::TestCaseError::fail("shape"))?;
        let ab = a.matmul(&b)
            .ok_or_else(|| proptest::test_runner::TestCaseError::fail("matmul"))?;
        let left = ab.matmul(&c)
            .ok_or_else(|| proptest::test_runner::TestCaseError::fail("matmul"))?;
        let bc = b.matmul(&c)
            .ok_or_else(|| proptest::test_runner::TestCaseError::fail("matmul"))?;
        let right = a.matmul(&bc)
            .ok_or_else(|| proptest::test_runner::TestCaseError::fail("matmul"))?;
        prop_assert!(array_eq(&left, &right));
    }

    #[test]
    fn rational_add_sub_inverse(
        n1 in -10_000i64..10_000,
        d1 in 1i64..10_000,
        n2 in -10_000i64..10_000,
        d2 in 1i64..10_000,
    ) {
        let a = ExactRational::from_i64(n1, d1);
        let b = ExactRational::from_i64(n2, d2);
        prop_assert!(rat_eq(&a.add(&b).sub(&b), &a));
    }
}
