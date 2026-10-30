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
    assert!(
        ExactNum::from_word(1, 0).is_nan(),
        "from_word(1, 0) must be NaN"
    );
    assert!(!ExactNum::new(MIN_PRECISION).is_nan());
}

const _: () = assert!(MAX_PRECISION == (1usize << 32) - 1);

#[test]
fn max_precision_is_the_largest_accepted_precision() {
    let over = ExactNum::new(MAX_PRECISION + 1);
    assert_eq!(over.err(), Some(Error::InvalidArgument));
}

#[test]
fn default_precision_and_decimal_scale_are_public() {
    assert_eq!(DEFAULT_P, 128);
    assert_eq!(MAX_DEC_SCALE, 1_000_000);
}

// `Precision` is the public name of the precision type: it holds the whole accepted range, and an exponent
// plus twice the largest precision still fits the 64-bit exponent arithmetic.
const _: zenith_float_num::Precision = zenith_float_num::MAX_PRECISION;
const _: zenith_float_num::Precision = zenith_float_num::MIN_PRECISION;
const _: zenith_float_num::Precision = zenith_float_num::DEFAULT_P;
const _: () = assert!(
    zenith_float_num::EXPONENT_MAX as i128 + 2 * zenith_float_num::MAX_PRECISION as i128
        <= i64::MAX as i128
);

#[test]
fn precision_alias_reports_precisions() {
    let x = zenith_float_num::ExactNum::from_u8(3, zenith_float_num::DEFAULT_P);
    let p: Option<zenith_float_num::Precision> = x.precision();
    let m: Option<zenith_float_num::Precision> = x.mantissa_max_bit_len();
    assert!(
        p == Some(zenith_float_num::DEFAULT_P) && m == p,
        "{p:?} {m:?}"
    );
}
