//! 2.0 public limits: the exponent and precision constants and what they mean.

use zenith_float_num::{
    Error, ExactNum, DEFAULT_P, EXPONENT_MIN, MAX_DEC_SCALE, MAX_PRECISION, MIN_NORMAL_EXPONENT,
    MIN_PRECISION, WORD_BIT_SIZE,
};

#[test]
fn min_normal_exponent_is_the_exponent_of_the_smallest_normal() {
    assert_eq!(MIN_NORMAL_EXPONENT, EXPONENT_MIN);
    for p in [1, WORD_BIT_SIZE, 3 * WORD_BIT_SIZE] {
        let n = ExactNum::min_positive_normal(p);
        assert!(!n.is_subnormal(), "p={p}");
        assert_eq!(n.exponent(), Some(MIN_NORMAL_EXPONENT), "p={p}");
        let s = ExactNum::min_positive(p);
        assert!(s.is_subnormal(), "p={p}");
        assert_eq!(s.exponent(), Some(MIN_NORMAL_EXPONENT), "p={p}");
    }
}

#[test]
fn precision_zero_is_rejected() {
    assert_eq!(MIN_PRECISION, 1);
    let z = ExactNum::new(0);
    assert!(z.is_nan(), "ExactNum::new(0) must be NaN");
    assert_eq!(z.err(), Some(Error::InvalidArgument));
    assert!(ExactNum::from_word(1, 0).is_nan(), "from_word(1, 0) must be NaN");
    assert!(!ExactNum::new(MIN_PRECISION).is_nan());
}

#[test]
fn max_precision_is_the_largest_accepted_precision() {
    let over = ExactNum::new(MAX_PRECISION + 1);
    assert_eq!(over.err(), Some(Error::InvalidArgument));
    assert!(MAX_PRECISION < usize::MAX / 4);
}

#[test]
fn default_precision_and_decimal_scale_are_public() {
    assert_eq!(DEFAULT_P, 128);
    assert_eq!(MAX_DEC_SCALE, 1_000_000);
}
