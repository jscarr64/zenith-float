//! MP18 integer getters and fits predicates (MPFR `get_sj/uj/si/ui`, rounded with `rm`).
//!
//! Edge-case list A to Z (every row below gets a test; "no answer" = Err, never a
//! saturated value — MPFR saturates, our rule forbids it, so the divergence is
//! documented on the functions):
//! A. Normal values: to_i64(5)=5; to_i64(-5)=-5; to_u64(5)=5; to_i32(5)=5; to_u32(5)=5.
//! B. Boundaries, exact: i64::MAX; i64::MIN; u64::MAX; i32::MAX; i32::MIN; u32::MAX;
//!    each fits in its own type and not in any smaller signed/unsigned type it exceeds.
//! C. Signed zeros: -0 -> 0 in every getter (to_i64, to_u64, to_i32, to_u32).
//! D. Infinities: +Inf, -Inf -> Err(InvalidArgument) in every getter; fits_* false.
//! E. NaN with error: NaN input -> Err; a NaN carrying an error propagates that error, a
//!    bare NaN gives Err(InvalidArgument); fits_* false. Never a saturated value.
//! F. Negatives: to_u64(-1) Err in every mode (a negative value rounds to a negative
//!    integer or to -0, and -0 is the only negative-side value that converts, to 0);
//!    to_u32(-1) Err; to_i64(-1)=-1; to_i32(-1)=-1.
//! G. Huge values: 2^63 -> to_i64 Err, to_u64 Ok(2^63); 2^64 -> to_u64 Err; 2^100 ->
//!    Err in every getter; fits_* false for all of them.
//! H. Tiny values: 2^-1000 -> 0 in every getter in ToZero/Down/ToEven (rounds to +0),
//!    and 1 in Up/FromZero (rounds to +1); -2^-1000 -> 0 in ToZero/Up/ToEven, -1 in
//!    Down/FromZero, and to_u64 gives 0 in ToZero/Up/ToEven but Err in Down/FromZero
//!    (the rounded value is -1, which does not fit an unsigned type). -0.5 to u64:
//!    0 in ToEven/ToZero/Up, Err in Down/ToOdd/FromZero.
//! I. Out-of-domain: 0.5 in each mode (fits in down, not in up — the pack's
//!    i64::MAX+0.5 row generalised): to_i64(i64::MAX+0.5) at ToOdd -> i64::MAX (tie
//!    to odd: i64::MAX is odd), at Down/ToZero -> i64::MAX, at Up/FromZero/ToEven ->
//!    Err (the tie rounds to the even neighbour 2^63, which does not fit i64);
//!    fits_i64 mirrors this (true in Down/ToZero/ToOdd, false in Up/FromZero/ToEven).
//! J. Every RoundingMode: 2.5 -> 2 ToEven, 3 ToOdd, 3 Up, 2 Down, 2 ToZero, 3 FromZero;
//!    -2.5 -> -2 ToEven, -3 ToOdd, -2 Up, -3 Down, -2 ToZero, -3 FromZero; 5.5 with
//!    RoundingMode::None -> Err (None does not round to an integer, so a non-integer
//!    input is refused; an already-integer input converts), fits false for 5.5 None.
//! K. Precisions: getters take no precision argument (rounding is to an integer, not
//!    to p bits), so the rows use operands at p=128 and at p=64 with identical results.
//! L. Golds: exact integers and the i64/u64/i32/u32 boundary values (identities from
//!    the integer types themselves); no decimal float golds are used.

use zenith_float_num::{Consts, Error, ExactNum, Radix, RoundingMode};

const MODES: [RoundingMode; 7] = [
    RoundingMode::ToEven,
    RoundingMode::ToOdd,
    RoundingMode::Up,
    RoundingMode::Down,
    RoundingMode::ToZero,
    RoundingMode::FromZero,
    RoundingMode::None,
];

fn parse(s: &str, p: usize, cc: &mut Consts) -> ExactNum {
    ExactNum::parse(s, Radix::Dec, p, RoundingMode::ToEven, cc)
}

fn expect_i64_err(got: Result<i64, Error>, want: Error, what: &str) {
    assert!(got.is_err(), "{what}: want Err({want:?}), got Ok");
    if let Err(e) = got {
        assert_eq!(e, want, "{what}: wrong error");
    }
}

fn expect_u64_err(got: Result<u64, Error>, want: Error, what: &str) {
    assert!(got.is_err(), "{what}: want Err({want:?}), got Ok");
    if let Err(e) = got {
        assert_eq!(e, want, "{what}: wrong error");
    }
}

fn expect_i32_err(got: Result<i32, Error>, want: Error, what: &str) {
    assert!(got.is_err(), "{what}: want Err({want:?}), got Ok");
    if let Err(e) = got {
        assert_eq!(e, want, "{what}: wrong error");
    }
}

fn expect_u32_err(got: Result<u32, Error>, want: Error, what: &str) {
    assert!(got.is_err(), "{what}: want Err({want:?}), got Ok");
    if let Err(e) = got {
        assert_eq!(e, want, "{what}: wrong error");
    }
}

#[test]
fn int_getters_normal_and_zeros() -> Result<(), Error> {
    let mut cc = Consts::new()?;
    let p = 128;
    assert_eq!(parse("5", p, &mut cc).to_i64(RoundingMode::ToEven)?, 5, "to_i64(5)");
    assert_eq!(parse("-5", p, &mut cc).to_i64(RoundingMode::ToEven)?, -5, "to_i64(-5)");
    assert_eq!(parse("5", p, &mut cc).to_u64(RoundingMode::ToEven)?, 5, "to_u64(5)");
    assert_eq!(parse("5", p, &mut cc).to_i32(RoundingMode::ToEven)?, 5, "to_i32(5)");
    assert_eq!(parse("5", p, &mut cc).to_u32(RoundingMode::ToEven)?, 5, "to_u32(5)");
    assert_eq!(parse("0", p, &mut cc).to_i64(RoundingMode::ToEven)?, 0, "to_i64(0)");
    let mut neg_zero = ExactNum::new(p);
    neg_zero.set_sign(zenith_float_num::Sign::Neg);
    assert_eq!(neg_zero.to_i64(RoundingMode::ToEven)?, 0, "to_i64(-0)");
    assert_eq!(neg_zero.to_u64(RoundingMode::ToEven)?, 0, "to_u64(-0)");
    assert_eq!(neg_zero.to_i32(RoundingMode::ToEven)?, 0, "to_i32(-0)");
    assert_eq!(neg_zero.to_u32(RoundingMode::ToEven)?, 0, "to_u32(-0)");
    // Same results at p=64.
    assert_eq!(parse("5", 64, &mut cc).to_i64(RoundingMode::ToEven)?, 5, "to_i64(5) p=64");
    assert!(parse("-5", 64, &mut cc).to_u32(RoundingMode::ToEven).is_err(), "to_u32(-5) p=64 Err");
    Ok(())
}

#[test]
fn int_getters_boundaries() -> Result<(), Error> {
    let mut cc = Consts::new()?;
    let p = 128;
    assert_eq!(
        parse("9223372036854775807", p, &mut cc).to_i64(RoundingMode::ToEven)?,
        i64::MAX,
        "to_i64(i64::MAX)"
    );
    assert_eq!(
        parse("-9223372036854775808", p, &mut cc).to_i64(RoundingMode::ToEven)?,
        i64::MIN,
        "to_i64(i64::MIN)"
    );
    assert_eq!(
        parse("18446744073709551615", p, &mut cc).to_u64(RoundingMode::ToEven)?,
        u64::MAX,
        "to_u64(u64::MAX)"
    );
    assert_eq!(
        parse("2147483647", p, &mut cc).to_i32(RoundingMode::ToEven)?,
        i32::MAX,
        "to_i32(i32::MAX)"
    );
    assert_eq!(
        parse("-2147483648", p, &mut cc).to_i32(RoundingMode::ToEven)?,
        i32::MIN,
        "to_i32(i32::MIN)"
    );
    assert_eq!(
        parse("4294967295", p, &mut cc).to_u32(RoundingMode::ToEven)?,
        u32::MAX,
        "to_u32(u32::MAX)"
    );
    // Cross-type: i64::MAX does not fit i32/u32; u64::MAX does not fit i64.
    expect_i32_err(
        parse("9223372036854775807", p, &mut cc).to_i32(RoundingMode::ToEven),
        Error::InvalidArgument,
        "to_i32(i64::MAX)",
    );
    expect_u32_err(
        parse("9223372036854775807", p, &mut cc).to_u32(RoundingMode::ToEven),
        Error::InvalidArgument,
        "to_u32(i64::MAX)",
    );
    expect_i64_err(
        parse("18446744073709551615", p, &mut cc).to_i64(RoundingMode::ToEven),
        Error::InvalidArgument,
        "to_i64(u64::MAX)",
    );
    // fits mirrors the getters at the boundaries.
    assert!(parse("9223372036854775807", p, &mut cc).fits_i64(RoundingMode::ToEven));
    assert!(!parse("9223372036854775807", p, &mut cc).fits_i32(RoundingMode::ToEven));
    assert!(!parse("9223372036854775807", p, &mut cc).fits_u32(RoundingMode::ToEven));
    assert!(parse("18446744073709551615", p, &mut cc).fits_u64(RoundingMode::ToEven));
    assert!(!parse("18446744073709551615", p, &mut cc).fits_i64(RoundingMode::ToEven));
    Ok(())
}

#[test]
fn int_getters_inf_nan_and_negatives() -> Result<(), Error> {
    let mut cc = Consts::new()?;
    let p = 128;
    let pos_inf = ExactNum::from_i8(1, p).div(&ExactNum::new(p), p, RoundingMode::ToEven);
    let neg_inf = ExactNum::from_i8(-1, p).div(&ExactNum::new(p), p, RoundingMode::ToEven);
    assert!(pos_inf.is_inf_pos() && neg_inf.is_inf_neg(), "setup infinities");
    for inf in [&pos_inf, &neg_inf] {
        expect_i64_err(inf.to_i64(RoundingMode::ToEven), Error::InvalidArgument, "to_i64(Inf)");
        expect_u64_err(inf.to_u64(RoundingMode::ToEven), Error::InvalidArgument, "to_u64(Inf)");
        assert!(!inf.fits_i64(RoundingMode::ToEven), "fits_i64(Inf) false");
        assert!(!inf.fits_u64(RoundingMode::ToEven), "fits_u64(Inf) false");
    }
    for nan in [
        ExactNum::nan(Some(Error::InvalidArgument)),
        ExactNum::nan(None),
    ] {
        expect_i64_err(nan.to_i64(RoundingMode::ToEven), Error::InvalidArgument, "to_i64(NaN)");
        assert!(!nan.fits_i64(RoundingMode::ToEven), "fits_i64(NaN) false");
    }
    expect_i64_err(
        ExactNum::nan(Some(Error::DivisionByZero)).to_i64(RoundingMode::ToEven),
        Error::DivisionByZero,
        "to_i64(NaN(Div)) propagates",
    );
    // Negatives into unsigned types: -1 Err in every mode, -0 -> 0 (covered in the
    // zeros test); -1 fits_i64 true, fits_u64 false.
    for rm in MODES {
        expect_u64_err(
            parse("-1", p, &mut cc).to_u64(rm),
            Error::InvalidArgument,
            "to_u64(-1) in every mode",
        );
        expect_u32_err(
            parse("-1", p, &mut cc).to_u32(rm),
            Error::InvalidArgument,
            "to_u32(-1) in every mode",
        );
        assert_eq!(parse("-1", p, &mut cc).to_i64(rm)?, -1, "to_i64(-1) in every mode");
        assert!(parse("-1", p, &mut cc).fits_i64(rm), "fits_i64(-1)");
        assert!(!parse("-1", p, &mut cc).fits_u64(rm), "fits_u64(-1) false");
    }
    Ok(())
}

#[test]
fn int_getters_huge_and_tiny() -> Result<(), Error> {
    let p = 128;
    let two63 = ExactNum::pow2(63, p)?;
    expect_i64_err(two63.to_i64(RoundingMode::ToEven), Error::InvalidArgument, "to_i64(2^63)");
    assert_eq!(two63.to_u64(RoundingMode::ToEven)?, 1u64 << 63, "to_u64(2^63)");
    assert!(two63.fits_u64(RoundingMode::ToEven));
    assert!(!two63.fits_i64(RoundingMode::ToEven));
    let two64 = ExactNum::pow2(64, p)?;
    expect_u64_err(two64.to_u64(RoundingMode::ToEven), Error::InvalidArgument, "to_u64(2^64)");
    assert!(!two64.fits_u64(RoundingMode::ToEven));
    let two100 = ExactNum::pow2(100, p)?;
    expect_i64_err(two100.to_i64(RoundingMode::ToEven), Error::InvalidArgument, "to_i64(2^100)");
    expect_u64_err(two100.to_u64(RoundingMode::ToEven), Error::InvalidArgument, "to_u64(2^100)");
    expect_i32_err(two100.to_i32(RoundingMode::ToEven), Error::InvalidArgument, "to_i32(2^100)");
    expect_u32_err(two100.to_u32(RoundingMode::ToEven), Error::InvalidArgument, "to_u32(2^100)");
    // Tiny: 2^-1000 rounds to +0 in ToZero/Down/ToEven, to +1 in Up/FromZero.
    let tiny = ExactNum::pow2(-1000, p)?;
    for rm in [RoundingMode::ToZero, RoundingMode::Down, RoundingMode::ToEven] {
        assert_eq!(tiny.to_i64(rm)?, 0, "to_i64(2^-1000) down modes");
        assert_eq!(tiny.to_u64(rm)?, 0, "to_u64(2^-1000) down modes");
    }
    for rm in [RoundingMode::Up, RoundingMode::FromZero] {
        assert_eq!(tiny.to_i64(rm)?, 1, "to_i64(2^-1000) up modes");
        assert_eq!(tiny.to_u64(rm)?, 1, "to_u64(2^-1000) up modes");
    }
    let neg_tiny = ExactNum::pow2(-1000, p)?.mul(&ExactNum::from_i8(-1, p), p, RoundingMode::ToEven);
    for rm in [RoundingMode::ToZero, RoundingMode::Up, RoundingMode::ToEven] {
        assert_eq!(neg_tiny.to_i64(rm)?, 0, "to_i64(-2^-1000) toward-zero modes");
        assert_eq!(neg_tiny.to_u64(rm)?, 0, "to_u64(-2^-1000) toward-zero modes");
    }
    for rm in [RoundingMode::Down, RoundingMode::FromZero] {
        assert_eq!(neg_tiny.to_i64(rm)?, -1, "to_i64(-2^-1000) away modes");
        expect_u64_err(neg_tiny.to_u64(rm), Error::InvalidArgument, "to_u64(-2^-1000) away modes");
    }
    Ok(())
}

#[test]
fn int_getters_max_plus_half_and_rounding_modes() -> Result<(), Error> {
    let mut cc = Consts::new()?;
    let p = 128;
    // i64::MAX + 0.5: tie between i64::MAX (odd) and 2^63 (even). ToOdd and
    // Down/ToZero give i64::MAX; ToEven, Up and FromZero give 2^63, which does
    // not fit i64 (probe: round0(MAX+0.5, ToEven) = 9223372036854775808).
    let x = parse("9223372036854775807.5", p, &mut cc);
    assert_eq!(x.to_i64(RoundingMode::ToOdd)?, i64::MAX, "to_i64(MAX+0.5) ToOdd");
    assert_eq!(x.to_i64(RoundingMode::Down)?, i64::MAX, "to_i64(MAX+0.5) Down");
    assert_eq!(x.to_i64(RoundingMode::ToZero)?, i64::MAX, "to_i64(MAX+0.5) ToZero");
    expect_i64_err(x.to_i64(RoundingMode::Up), Error::InvalidArgument, "to_i64(MAX+0.5) Up");
    expect_i64_err(x.to_i64(RoundingMode::FromZero), Error::InvalidArgument, "to_i64(MAX+0.5) FromZero");
    expect_i64_err(x.to_i64(RoundingMode::ToEven), Error::InvalidArgument, "to_i64(MAX+0.5) ToEven");
    assert!(x.fits_i64(RoundingMode::ToOdd));
    assert!(x.fits_i64(RoundingMode::Down));
    assert!(x.fits_i64(RoundingMode::ToZero));
    assert!(!x.fits_i64(RoundingMode::Up));
    assert!(!x.fits_i64(RoundingMode::FromZero));
    assert!(!x.fits_i64(RoundingMode::ToEven));
    // -0.5 to u64 by mode: 0 in ToEven/ToZero/Up (tie to even is 0, toward zero is 0,
    // toward +inf is 0); Err in Down/ToOdd/FromZero (all round to -1).
    let neg_half = parse("-0.5", p, &mut cc);
    for rm in [RoundingMode::ToEven, RoundingMode::ToZero, RoundingMode::Up] {
        assert_eq!(neg_half.to_u64(rm)?, 0, "to_u64(-0.5) zero modes");
        assert!(neg_half.fits_u64(rm), "fits_u64(-0.5) zero modes");
    }
    for rm in [RoundingMode::Down, RoundingMode::ToOdd, RoundingMode::FromZero] {
        expect_u64_err(neg_half.to_u64(rm), Error::InvalidArgument, "to_u64(-0.5) minus-one modes");
        assert!(!neg_half.fits_u64(rm), "fits_u64(-0.5) minus-one modes false");
    }
    // 2.5 / -2.5 in every mode.
    let pos = parse("2.5", p, &mut cc);
    assert_eq!(pos.to_i64(RoundingMode::ToEven)?, 2);
    assert_eq!(pos.to_i64(RoundingMode::ToOdd)?, 3);
    assert_eq!(pos.to_i64(RoundingMode::Up)?, 3);
    assert_eq!(pos.to_i64(RoundingMode::Down)?, 2);
    assert_eq!(pos.to_i64(RoundingMode::ToZero)?, 2);
    assert_eq!(pos.to_i64(RoundingMode::FromZero)?, 3);
    let neg = parse("-2.5", p, &mut cc);
    assert_eq!(neg.to_i64(RoundingMode::ToEven)?, -2);
    assert_eq!(neg.to_i64(RoundingMode::ToOdd)?, -3);
    assert_eq!(neg.to_i64(RoundingMode::Up)?, -2);
    assert_eq!(neg.to_i64(RoundingMode::Down)?, -3);
    assert_eq!(neg.to_i64(RoundingMode::ToZero)?, -2);
    assert_eq!(neg.to_i64(RoundingMode::FromZero)?, -3);
    // None: a non-integer input is refused, an integer input converts.
    expect_i64_err(parse("5.5", p, &mut cc).to_i64(RoundingMode::None), Error::InvalidArgument, "to_i64(5.5) None");
    assert_eq!(parse("5", p, &mut cc).to_i64(RoundingMode::None)?, 5, "to_i64(5) None");
    assert!(!parse("5.5", p, &mut cc).fits_i64(RoundingMode::None), "fits_i64(5.5) None false");
    assert!(parse("5", p, &mut cc).fits_i64(RoundingMode::None), "fits_i64(5) None true");
    Ok(())
}
