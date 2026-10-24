//! R14 (zenith-bound bridge): `ExactNum::as_raw_parts` / `ExactNum::from_raw_parts`.
//! Round-trips over the whole representable set and named rejections of malformed parts.

use std::ops::Neg;
use zenith_float_num::{
    Error, ExactNum, RawKind, RawParts, RoundingMode, Sign, Word, EXPONENT_MAX, EXPONENT_MIN,
    INF_NEG, INF_POS, NAN, WORD_BIT_SIZE,
};

fn same(a: &RawParts<'_>, b: &RawParts<'_>) -> bool {
    a.kind == b.kind
        && a.words == b.words
        && a.precision == b.precision
        && a.negative == b.negative
        && a.exponent == b.exponent
        && a.inexact == b.inexact
}

/// `x -> parts -> y -> parts` gives the same parts, and `y` compares equal to `x`.
fn round_trip(x: &ExactNum) -> Result<(), String> {
    let a = x.as_raw_parts();
    let y = ExactNum::from_raw_parts(&a).map_err(|e| format!("rejected own parts {a:?}: {e:?}"))?;
    let b = y.as_raw_parts();
    if !same(&a, &b) {
        return Err(format!("parts differ: {a:?} vs {b:?}"));
    }
    if !x.is_nan() && x.compare(&y) != Some(0) {
        return Err("value differs".to_string());
    }
    Ok(())
}

fn finite<'a>(words: &'a [Word], precision: usize, exponent: i64) -> RawParts<'a> {
    RawParts {
        kind: RawKind::Finite,
        words,
        precision,
        negative: false,
        exponent,
        inexact: false,
    }
}

fn rejected(parts: &RawParts<'_>) -> bool {
    matches!(ExactNum::from_raw_parts(parts), Err(Error::InvalidArgument))
}

const TOP: Word = 1 << (WORD_BIT_SIZE - 1);

#[test]
fn round_trip_exponent_extremes() -> Result<(), String> {
    assert_eq!(EXPONENT_MAX, (1i64 << 62) - 1);
    assert_eq!(EXPONENT_MIN, -((1i64 << 62) - 1));
    for e in [EXPONENT_MAX, EXPONENT_MAX - 1, EXPONENT_MIN, EXPONENT_MIN + 1, 0, 1, -1] {
        for s in [Sign::Pos, Sign::Neg] {
            for m in [[TOP, 0], [TOP | 5, Word::MAX], [Word::MAX, Word::MAX]] {
                let x = ExactNum::from_words(&m, s, e);
                round_trip(&x)?;
                let parts = x.as_raw_parts();
                assert!(parts.kind == RawKind::Finite && parts.exponent == e);
                assert_eq!(parts.negative, s == Sign::Neg);
                assert_eq!(parts.words, &m);
                assert_eq!(parts.precision, 2 * WORD_BIT_SIZE);
            }
        }
    }
    round_trip(&ExactNum::max_value(256))?;
    round_trip(&ExactNum::max_value(256).neg())
}

#[test]
fn round_trip_subnormals() -> Result<(), String> {
    for p in [1usize, 64, 65, 128, 1000] {
        let x = ExactNum::min_positive(p);
        let parts = x.as_raw_parts();
        assert_eq!(parts.exponent, EXPONENT_MIN);
        assert_eq!(parts.words.first(), Some(&1));
        round_trip(&x)?;
        round_trip(&x.neg())?;
    }
    // Leading zero words at EXPONENT_MIN are a valid subnormal form.
    let m: [Word; 3] = [7, 0, 0];
    let parts = finite(&m, 3 * WORD_BIT_SIZE, EXPONENT_MIN);
    let x = ExactNum::from_raw_parts(&parts).map_err(|e| format!("{e:?}"))?;
    assert!(same(&x.as_raw_parts(), &parts));
    Ok(())
}

fn zero(p: usize) -> ExactNum {
    ExactNum::min_positive(p).sub(&ExactNum::min_positive(p), p, RoundingMode::ToEven)
}

#[test]
fn round_trip_zero_signs() -> Result<(), String> {
    for p in [1usize, 64, 192] {
        for (x, negative) in [(zero(p), false), (zero(p).neg(), true)] {
            let parts = x.as_raw_parts();
            assert!(parts.kind == RawKind::Zero && parts.negative == negative);
            assert!(parts.words.iter().all(|w| *w == 0));
            round_trip(&x)?;
        }
    }
    Ok(())
}

#[test]
fn round_trip_inf_nan() -> Result<(), String> {
    for (x, kind, negative) in [
        (INF_POS, RawKind::Inf, false),
        (INF_NEG, RawKind::Inf, true),
        (NAN, RawKind::NaN, false),
    ] {
        let parts = x.as_raw_parts();
        assert!(parts.kind == kind && parts.negative == negative);
        assert!(parts.words.is_empty() && parts.precision == 0 && parts.exponent == 0);
        round_trip(&x)?;
    }
    let y = ExactNum::from_raw_parts(&NAN.as_raw_parts()).map_err(|e| format!("{e:?}"))?;
    assert!(y.is_nan());
    Ok(())
}

#[test]
fn round_trip_precisions_1_to_65536() -> Result<(), String> {
    let mut ps: Vec<usize> = (1..=200).collect();
    ps.extend([255, 256, 257, 1000, 4095, 4096, 4097, 65535, 65536]);
    for p in ps {
        for x in [
            ExactNum::max_value(p),
            ExactNum::min_positive(p),
            ExactNum::min_positive_normal(p).neg(),
        ] {
            let parts = x.as_raw_parts();
            assert_eq!(
                parts.precision,
                p.div_ceil(WORD_BIT_SIZE) * WORD_BIT_SIZE,
                "p={p}"
            );
            assert_eq!(parts.words.len() * WORD_BIT_SIZE, parts.precision);
            round_trip(&x).map_err(|e| format!("p={p}: {e}"))?;
        }
    }
    Ok(())
}

#[test]
fn round_trip_keeps_inexact() -> Result<(), String> {
    let third =
        ExactNum::from_i64(1, 128).div(&ExactNum::from_i64(3, 128), 128, RoundingMode::ToEven);
    assert!(third.as_raw_parts().inexact);
    round_trip(&third)
}

#[test]
fn reject_non_normalized_mantissa() {
    let m: [Word; 2] = [Word::MAX, TOP >> 1];
    assert!(rejected(&finite(&m, 2 * WORD_BIT_SIZE, 5)));
    assert!(rejected(&finite(&m, 2 * WORD_BIT_SIZE, EXPONENT_MIN + 1)));
    // All-zero words are not a finite non-zero value.
    let z: [Word; 1] = [0];
    assert!(rejected(&finite(&z, WORD_BIT_SIZE, 0)));
}

#[test]
fn reject_precision_zero() {
    let empty: [Word; 0] = [];
    assert!(rejected(&finite(&empty, 0, 1)));
    let one: [Word; 1] = [TOP];
    assert!(rejected(&finite(&one, 0, 1)));
}

#[test]
fn reject_word_count_not_matching_precision() {
    let one: [Word; 1] = [TOP];
    let two: [Word; 2] = [0, TOP];
    assert!(rejected(&finite(&one, 2 * WORD_BIT_SIZE, 1)));
    assert!(rejected(&finite(&two, WORD_BIT_SIZE, 1)));
    // Precision is the stored width, a whole number of words.
    assert!(rejected(&finite(&one, WORD_BIT_SIZE - 1, 1)));
    assert!(rejected(&finite(&two, WORD_BIT_SIZE + 1, 1)));
}

#[test]
fn reject_exponent_out_of_range() {
    let one: [Word; 1] = [TOP];
    let two62 = 1i64 << 62;
    assert!(rejected(&finite(&one, WORD_BIT_SIZE, two62)));
    assert!(rejected(&finite(&one, WORD_BIT_SIZE, -two62)));
    assert!(rejected(&finite(&one, WORD_BIT_SIZE, i64::MAX)));
    assert!(rejected(&finite(&one, WORD_BIT_SIZE, i64::MIN)));
    assert!(!rejected(&finite(&one, WORD_BIT_SIZE, two62 - 1)));
    assert!(!rejected(&finite(&one, WORD_BIT_SIZE, -(two62 - 1))));
}

#[test]
fn reject_malformed_special_parts() {
    let one: [Word; 1] = [TOP];
    let zero: [Word; 1] = [0];
    let empty: [Word; 0] = [];
    let with = |kind, words, precision, exponent| RawParts {
        kind,
        words,
        precision,
        negative: false,
        exponent,
        inexact: false,
    };
    assert!(rejected(&with(RawKind::Zero, &one, WORD_BIT_SIZE, 0)));
    assert!(rejected(&with(RawKind::Zero, &zero, WORD_BIT_SIZE, 3)));
    assert!(rejected(&with(RawKind::Zero, &empty, 0, 0)));
    assert!(rejected(&with(RawKind::Inf, &one, WORD_BIT_SIZE, 0)));
    assert!(rejected(&with(RawKind::NaN, &empty, 0, 7)));
    let inexact_inf = RawParts {
        inexact: true,
        ..with(RawKind::Inf, &empty, 0, 0)
    };
    assert!(rejected(&inexact_inf));
}
