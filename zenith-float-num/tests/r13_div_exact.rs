//! R13: `ExactNum::div_exact` returns the quotient only when it is exact at `p` bits.
//! The operands are 64-bit words.

use zenith_float_num::{Error, ExactNum, RoundingMode, INF_NEG, INF_POS, NAN};

const P: usize = 64;

fn n(v: u64) -> ExactNum {
    ExactNum::from_word(v, P)
}

fn same(a: &ExactNum, b: &ExactNum) -> bool {
    a.mantissa_digits() == b.mantissa_digits()
        && a.exponent() == b.exponent()
        && a.sign() == b.sign()
}

#[test]
fn div_exact_returns_exact_quotients() -> Result<(), Error> {
    for (a, b, q) in [(6u64, 3u64, 2u64), (15, 5, 3), (1 << 40, 1 << 20, 1 << 20), (21, 7, 3)] {
        let got = n(a).div_exact(&n(b), P)?.ok_or(Error::InvalidArgument)?;
        assert!(same(&got, &n(q)), "{a}/{b}");
        assert!(!got.inexact(), "{a}/{b} must not be inexact");
    }
    // 3/4 = 0.75 and 1/1024 are dyadic: exact.
    let q = n(3).div_exact(&n(4), P)?.ok_or(Error::InvalidArgument)?;
    assert!(same(&q, &n(3).div(&n(4), P, RoundingMode::ToEven)));
    assert!(n(1).div_exact(&n(1024), P)?.is_some());
    // sign
    let q = (-n(6)).div_exact(&n(3), P)?.ok_or(Error::InvalidArgument)?;
    assert!(same(&q, &-n(2)));
    Ok(())
}

#[test]
fn div_exact_reports_inexact_quotients() -> Result<(), Error> {
    assert!(n(1).div_exact(&n(3), P)?.is_none());
    assert!(n(22).div_exact(&n(7), P)?.is_none());
    // 2^64-1 over 2 needs 64 bits: exact at 64, not at 24 bits (rounded to a word: 64), so
    // use 2^65-ish: (2^64-1)/2^-1 is exact; (2^64-1) at 128 bits / 3 is exact (divisible).
    assert!(n(u64::MAX).div_exact(&n(3), P)?.is_some());
    Ok(())
}

#[test]
fn div_exact_range_ends() -> Result<(), Error> {
    let big = ExactNum::max_value(P);
    let tiny = ExactNum::min_positive(P);
    // overflow and underflow are not exact
    assert!(big.div_exact(&tiny, P)?.is_none());
    assert!(tiny.div_exact(&big, P)?.is_none());
    // the subnormal grid: min_positive / 1 is exact; min_positive / 2 is not
    assert!(tiny.div_exact(&n(1), P)?.is_some());
    assert!(tiny.div_exact(&n(2), P)?.is_none());
    Ok(())
}

#[test]
fn div_exact_specials() -> Result<(), Error> {
    assert!(matches!(
        n(1).div_exact(&n(0), P),
        Err(Error::DivisionByZero)
    ));
    assert!(matches!(
        n(0).div_exact(&n(0), P),
        Err(Error::InvalidArgument)
    ));
    assert!(NAN.div_exact(&n(1), P).is_err());
    assert!(n(1).div_exact(&NAN, P).is_err());
    for (a, b) in [
        (&INF_POS, &INF_NEG),
        (&INF_POS, &INF_POS),
        (&INF_NEG, &INF_NEG),
        (&INF_NEG, &INF_POS),
    ] {
        assert!(matches!(a.div_exact(b, P), Err(Error::InvalidArgument)));
    }
    let q = INF_POS
        .div_exact(&-n(2), P)?
        .ok_or(Error::InvalidArgument)?;
    assert!(q.is_inf_neg());
    let q = n(2).div_exact(&INF_NEG, P)?.ok_or(Error::InvalidArgument)?;
    assert!(q.is_zero() && q.sign() == Some(zenith_float_num::Sign::Neg));
    let q = n(0).div_exact(&n(5), P)?.ok_or(Error::InvalidArgument)?;
    assert!(q.is_zero());
    assert!(matches!(
        n(1).div_exact(&n(3), 0),
        Err(Error::InvalidArgument)
    ));
    Ok(())
}

/// The precision is validated before the operands: `x / 0` and `0 / 0` at `p = 0` report
/// InvalidArgument, not DivisionByZero.
#[test]
fn div_exact_precision_zero_wins_over_division_by_zero() {
    assert!(matches!(
        n(1).div_exact(&n(0), 0),
        Err(Error::InvalidArgument)
    ));
    assert!(matches!(
        (-n(5)).div_exact(&n(0), 0),
        Err(Error::InvalidArgument)
    ));
    assert!(matches!(
        n(0).div_exact(&n(0), 0),
        Err(Error::InvalidArgument)
    ));
    assert!(matches!(
        INF_POS.div_exact(&n(0), 0),
        Err(Error::InvalidArgument)
    ));
    assert!(matches!(
        n(1).div_exact(&n(3), 0),
        Err(Error::InvalidArgument)
    ));
}
