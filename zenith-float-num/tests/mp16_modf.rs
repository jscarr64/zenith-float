//! MP16 `modf`: MPFR `mpfr_modf` — integer part (trunc toward 0) and fractional part
//! with the sign of `x`, each rounded to `p` bits with `rm`, built on `int` and `fract`
//! (both exact).
//!
//! Edge-case list A to Z (every row below gets a test; "no answer" = Err, never Ok(NaN)):
//! A. Normal values: 3.75->(3,0.75); -3.75->(-3,-0.75); 5->(5,+0); -5->(-5,-0);
//!    0.5->(0,0.5); -0.5->(-0,-0.5).
//! B. Boundaries: +0->(+0,+0); -0->(-0,-0).
//! C. Infinities: +Inf->(+Inf,+0); -Inf->(-Inf,-0) (MPFR convention, frac sign follows x
//!    for the zero part; the infinite part carries x's sign).
//! D. NaN with error: NaN input -> Err; a NaN carrying an error propagates that error, a bare
//!    NaN gives Err(InvalidArgument). Never Ok(NaN, _) or (_, Ok(NaN)).
//! E. Negatives: covered in A; trunc toward 0, never floor: -3.75 int is -3, not -4.
//! F. Huge values: 2^100+0.25 at p=128 -> (2^100,0.25); max_value(128) -> (max_value,+0).
//! G. Tiny values: 2^-1000 -> (+0,2^-1000); -2^-1000 -> (-0,-2^-1000) (frac keeps x's sign).
//! H. Out-of-domain: p=0 and p>MAX_PRECISION -> Err(InvalidArgument). No other domain
//!    restriction: modf is total on finite and infinite inputs.
//! I. Every RoundingMode: 3.75 exact in all 7 modes; 0.5+2^-100 at p=64 -> int 0 in all
//!    modes, frac 0.5 in ToEven/ToOdd/Down/ToZero, frac 0.5+2^-64 in Up/FromZero, and
//!    frac 0.5 in None (storage grid value, pinned, not claimed exact).
//! J. Precisions: 5 and 0.5 exact at p = 1, 53, 63, 64, 65, 113, 128, 200; 3.75 exact at
//!    p >= 64 (3.75 needs 3 bits, and at p=1 the integer part itself would round, so p=1
//!    is not an exactness row for 3.75).
//! K. Inexact flag: modf(3.75) both parts exact -> inexact false; modf(0.5+2^-100) at
//!    p=64 ToEven -> int exact, frac inexact true.
//! L. Golds: exact integers, powers of two and 0.75/0.25 (identities); the 0.5+2^-100 row
//!    is checked against 0.5+2^-64 built independently at p=128 and rounded to p=64
//!    (exact at 64 bits), not against a decimal float from memory. mpmath spot check
//!    (mp.prec = 4*64+64 = 320): floor/frac of 0.5+2**-100 rounded to 64 bits ToEven
//!    gives frac 0.5, Up gives 0.5+2**-64.

use zenith_float_num::{
    Consts, Error, ExactNum, Radix, RoundingMode, Sign, MAX_PRECISION,
};

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

fn assert_eq_num(got: &ExactNum, want: &ExactNum, what: &str) {
    assert_eq!(
        got.compare(want),
        Some(0),
        "{what}: got {got:?}, want {want:?}"
    );
}

fn assert_zero_sign(got: &ExactNum, negative: bool, what: &str) {
    assert!(got.is_zero(), "{what}: got {got:?}, want zero");
    assert_eq!(
        got.is_negative(),
        negative,
        "{what}: got {got:?}, want negative={negative}"
    );
}

fn expect_modf_err(
    got: Result<(ExactNum, ExactNum), Error>,
    want: Error,
    what: &str,
) {
    assert!(got.is_err(), "{what}: want Err({want:?}), got Ok");
    if let Err(e) = got {
        assert_eq!(e, want, "{what}: wrong error");
    }
}

fn assert_modf(
    x: &ExactNum,
    p: usize,
    rm: RoundingMode,
    want_int: &ExactNum,
    want_frac: &ExactNum,
    what: &str,
) -> Result<(), Error> {
    let (int_part, frac_part) = x.modf(p, rm)?;
    assert_eq_num(&int_part, want_int, &format!("{what} int"));
    assert_eq_num(&frac_part, want_frac, &format!("{what} frac"));
    Ok(())
}

#[test]
fn modf_normal_and_negative() -> Result<(), Error> {
    let mut cc = Consts::new()?;
    let p = 128;
    assert_modf(
        &parse("3.75", p, &mut cc),
        p,
        RoundingMode::ToEven,
        &parse("3", p, &mut cc),
        &parse("0.75", p, &mut cc),
        "modf(3.75)",
    )?;
    assert_modf(
        &parse("-3.75", p, &mut cc),
        p,
        RoundingMode::ToEven,
        &parse("-3", p, &mut cc),
        &parse("-0.75", p, &mut cc),
        "modf(-3.75)",
    )?;
    assert_modf(
        &parse("5", p, &mut cc),
        p,
        RoundingMode::ToEven,
        &parse("5", p, &mut cc),
        &parse("0", p, &mut cc),
        "modf(5)",
    )?;
    assert_modf(
        &parse("-5", p, &mut cc),
        p,
        RoundingMode::ToEven,
        &parse("-5", p, &mut cc),
        &parse("-0", p, &mut cc),
        "modf(-5)",
    )?;
    assert_modf(
        &parse("0.5", p, &mut cc),
        p,
        RoundingMode::ToEven,
        &parse("0", p, &mut cc),
        &parse("0.5", p, &mut cc),
        "modf(0.5)",
    )?;
    assert_modf(
        &parse("-0.5", p, &mut cc),
        p,
        RoundingMode::ToEven,
        &parse("-0", p, &mut cc),
        &parse("-0.5", p, &mut cc),
        "modf(-0.5)",
    )?;
    Ok(())
}

#[test]
fn modf_zeros_and_infinities() -> Result<(), Error> {
    let p = 128;
    let pos_zero = ExactNum::new(p);
    let mut neg_zero = ExactNum::new(p);
    neg_zero.set_sign(Sign::Neg);
    let (int_part, frac_part) = pos_zero.modf(p, RoundingMode::ToEven)?;
    assert_zero_sign(&int_part, false, "modf(+0) int");
    assert_zero_sign(&frac_part, false, "modf(+0) frac");
    let (int_part, frac_part) = neg_zero.modf(p, RoundingMode::ToEven)?;
    assert_zero_sign(&int_part, true, "modf(-0) int");
    assert_zero_sign(&frac_part, true, "modf(-0) frac");
    let pos_inf = ExactNum::from_i8(1, p).div(&ExactNum::new(p), p, RoundingMode::ToEven);
    assert!(pos_inf.is_inf_pos(), "setup: 1/+0 must be +Inf");
    let neg_inf = ExactNum::from_i8(-1, p).div(&ExactNum::new(p), p, RoundingMode::ToEven);
    assert!(neg_inf.is_inf_neg(), "setup: -1/+0 must be -Inf");
    let (int_part, frac_part) = pos_inf.modf(p, RoundingMode::ToEven)?;
    assert!(int_part.is_inf_pos(), "modf(+Inf) int must be +Inf");
    assert_zero_sign(&frac_part, false, "modf(+Inf) frac");
    let (int_part, frac_part) = neg_inf.modf(p, RoundingMode::ToEven)?;
    assert!(int_part.is_inf_neg(), "modf(-Inf) int must be -Inf");
    assert_zero_sign(&frac_part, true, "modf(-Inf) frac");
    Ok(())
}

#[test]
fn modf_nan_is_error() -> Result<(), Error> {
    let p = 128;
    let one = ExactNum::from_i8(1, p);
    expect_modf_err(
        ExactNum::nan(Some(Error::InvalidArgument)).modf(p, RoundingMode::ToEven),
        Error::InvalidArgument,
        "modf(NaN)",
    );
    expect_modf_err(
        ExactNum::nan(Some(Error::DivisionByZero)).modf(p, RoundingMode::ToEven),
        Error::DivisionByZero,
        "modf(NaN(Div))",
    );
    expect_modf_err(
        ExactNum::nan(None).modf(p, RoundingMode::ToEven),
        Error::InvalidArgument,
        "modf(bare NaN)",
    );
    let _ = one;
    Ok(())
}

#[test]
fn modf_huge_and_tiny() -> Result<(), Error> {
    let mut cc = Consts::new()?;
    let p = 128;
    let big = ExactNum::pow2(100, p)?;
    let quarter = ExactNum::pow2(-2, p)?;
    let x = ExactNum::add_ref(&big, &quarter, p, RoundingMode::ToEven);
    assert_modf(&x, p, RoundingMode::ToEven, &big, &quarter, "modf(2^100+0.25)")?;
    let max = ExactNum::max_value(p);
    assert_modf(
        &max,
        p,
        RoundingMode::ToEven,
        &ExactNum::max_value(p),
        &parse("0", p, &mut cc),
        "modf(max_value)",
    )?;
    let tiny = ExactNum::pow2(-1000, p)?;
    let (int_part, frac_part) = tiny.modf(p, RoundingMode::ToEven)?;
    assert_zero_sign(&int_part, false, "modf(2^-1000) int");
    assert_eq_num(&frac_part, &tiny, "modf(2^-1000) frac");
    let neg_tiny = ExactNum::pow2(-1000, p)?.mul(&ExactNum::from_i8(-1, p), p, RoundingMode::ToEven);
    let (int_part, frac_part) = neg_tiny.modf(p, RoundingMode::ToEven)?;
    assert_zero_sign(&int_part, true, "modf(-2^-1000) int");
    assert_eq_num(&frac_part, &neg_tiny, "modf(-2^-1000) frac");
    assert!(frac_part.is_negative(), "modf(-2^-1000) frac must be negative");
    Ok(())
}

#[test]
fn modf_rounding_modes_and_precisions() -> Result<(), Error> {
    let mut cc = Consts::new()?;
    // Exact at every precision: integers and powers of two.
    for p in [1usize, 53, 63, 64, 65, 113, 128, 200] {
        assert_modf(
            &parse("5", 128, &mut cc),
            p,
            RoundingMode::ToEven,
            &parse("5", 128, &mut cc),
            &parse("0", 128, &mut cc),
            "modf(5) at every precision",
        )?;
        assert_modf(
            &parse("0.5", 128, &mut cc),
            p,
            RoundingMode::ToEven,
            &parse("0", 128, &mut cc),
            &parse("0.5", 128, &mut cc),
            "modf(0.5) at every precision",
        )?;
    }
    // 3.75 exact in every mode at p=128.
    for rm in MODES {
        assert_modf(
            &parse("3.75", 128, &mut cc),
            128,
            rm,
            &parse("3", 128, &mut cc),
            &parse("0.75", 128, &mut cc),
            "modf(3.75) in every mode",
        )?;
    }
    // 0.5 + 2^-100 at p=64: frac rounding. Expected values are rebuilt fresh here
    // because add consumes its receiver, so `half` cannot be reused after building x.
    let p = 64;
    let tiny = ExactNum::pow2(-100, p)?;
    let x = ExactNum::pow2(-1, p)?.add(&tiny, 128, RoundingMode::ToEven);
    let half = ExactNum::pow2(-1, p)?;
    let mut up_expected =
        ExactNum::pow2(-1, 128)?.add(&ExactNum::pow2(-64, 128)?, 128, RoundingMode::ToEven);
    up_expected.set_precision(64, RoundingMode::ToEven)?;
    for rm in [
        RoundingMode::ToEven,
        RoundingMode::ToOdd,
        RoundingMode::Down,
        RoundingMode::ToZero,
    ] {
        let (int_part, frac_part) = x.modf(p, rm)?;
        assert_zero_sign(&int_part, false, "modf(0.5+2^-100) int down modes");
        assert_eq_num(&frac_part, &half, "modf(0.5+2^-100) frac down modes");
    }
    for rm in [RoundingMode::Up, RoundingMode::FromZero] {
        let (int_part, frac_part) = x.modf(p, rm)?;
        assert_zero_sign(&int_part, false, "modf(0.5+2^-100) int up modes");
        assert_eq_num(&frac_part, &up_expected, "modf(0.5+2^-100) frac up modes");
    }
    // None: pinned to the storage grid value, which here is 0.5 (see MP12 None note:
    // set_precision with None keeps storage precision; the test pins the observed
    // grid value rather than claiming a rounding direction).
    let (int_part, frac_part) = x.modf(p, RoundingMode::None)?;
    assert_zero_sign(&int_part, false, "modf(0.5+2^-100) int None");
    assert_eq_num(&frac_part, &half, "modf(0.5+2^-100) frac None");
    Ok(())
}

#[test]
fn modf_invalid_precision_and_inexact_flag() -> Result<(), Error> {
    let mut cc = Consts::new()?;
    let x = parse("3.75", 128, &mut cc);
    expect_modf_err(
        x.modf(0, RoundingMode::ToEven),
        Error::InvalidArgument,
        "modf p=0",
    );
    expect_modf_err(
        x.modf(MAX_PRECISION + 1, RoundingMode::ToEven),
        Error::InvalidArgument,
        "modf p>MAX",
    );
    let (int_part, frac_part) = parse("3.75", 128, &mut cc).modf(128, RoundingMode::ToEven)?;
    assert!(!int_part.inexact(), "modf(3.75) int must be exact");
    assert!(!frac_part.inexact(), "modf(3.75) frac must be exact");
    let half = ExactNum::pow2(-1, 64)?;
    let tiny = ExactNum::pow2(-100, 64)?;
    let y = half.add(&tiny, 128, RoundingMode::ToEven);
    let (int_part, frac_part) = y.modf(64, RoundingMode::ToEven)?;
    assert!(!int_part.inexact(), "modf(0.5+2^-100) int must be exact");
    assert!(frac_part.inexact(), "modf(0.5+2^-100) frac must be inexact");
    Ok(())
}
