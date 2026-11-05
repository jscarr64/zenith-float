//! MP12 `dim`: MPFR `mpfr_dim` — `x - y` rounded if `x > y`, else `+0`.
//!
//! Edge-case list A to Z (every row below gets a test; "no answer" = Err, never a bare NaN):
//! A. Normal values: dim(5,3)=2; dim(3,5)=+0; dim(3,3)=+0; dim(3,-3)=6; dim(-3,-5)=2.
//! B. Boundaries: dim(5,5)=+0; dim(-5,-5)=+0; dim(0,0)=+0.
//! C. Signed zeros: dim(+0,-0)=+0; dim(-0,+0)=+0; dim(-0,-0)=+0 (result is always +0 when x<=y).
//! D. Infinities: dim(+Inf,+Inf)=+0; dim(+Inf,1)=+Inf; dim(1,-Inf)=+Inf; dim(-Inf,1)=+0;
//!    dim(1,+Inf)=+0; dim(-Inf,-Inf)=+0; dim(+Inf,-Inf)=+Inf; dim(-Inf,+Inf)=+0.
//! E. NaN with error: NaN input -> Err; a NaN carrying an error propagates that error, a bare
//!    NaN gives Err(InvalidArgument). Never Ok(NaN).
//! F. Negatives: dim(-3,3)=+0; dim(-5,-3)=+0; dim(3,-3)=6; dim(-3,-5)=2.
//! G. Huge values: dim(max_value(128),-max_value(128)) at ToEven -> Err(ExponentOverflow(Pos))
//!    (the exact difference overflows and ToEven rounds away from zero); at ToZero the same
//!    pair -> Ok(max_value(128)) (overflow follows the rounding mode, tests/overflow_modes.rs).
//! H. Tiny values: dim(2^-1000,2^-1001)=2^-1001; dim(min_positive(128),-min_positive(128))
//!    = 2*min_positive(128); dim(pow2(1000),pow2(999))=pow2(999), all at p=128.
//! I. Out-of-domain: none beyond NaN (dim is total on finite/Inf inputs); p=0 and
//!    p>MAX_PRECISION -> Err(InvalidArgument).
//! J. Every RoundingMode: dim(5,3)=2 and dim(3,5)=+0 in all 7 modes; dim(1,2^-200) at p=64
//!    -> 1 in ToEven/ToOdd/Up/FromZero, and 1-2^-64 in Down/ToZero (the exact value 1-2^-200
//!    lies in [0.5,1), whose spacing at p=64 is 2^-64); RoundingMode::None keeps the exact
//!    difference 1-2^-200 (no rounding to p bits is applied by sub with None).
//! K. Precisions: exact cases at p = 1, 53, 63, 64, 65, 113, 128, 200.
//! L. Inexact flag: dim(5,3) exact -> inexact false; dim(1,2^-200) at p=64 ToEven -> inexact true.
//! M. Golds: exact integers and powers of two (identities), overflow rows from the rounding-mode
//!    convention in tests/overflow_modes.rs; the 1-2^-200 row is checked against the value
//!    1-2^-64 built independently at p=128 and rounded to p=64 (exact at 64 bits), not against
//!    a decimal float from memory. mpmath spot check (mp.prec = 4*64+64 = 320):
//!    mpmath: 1 - mpf(2)**-200 rounded to 64 bits ToEven -> 1; Down -> 1 - 2**-64.

use zenith_float_num::{
    Consts, Error, ExactNum, Radix, RoundingMode, Sign, EXPONENT_MAX, MAX_PRECISION,
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

fn assert_pos_zero(got: &ExactNum, what: &str) {
    assert!(got.is_zero(), "{what}: got {got:?}, want +0");
    assert!(!got.is_negative(), "{what}: got {got:?}, want +0 not -0");
}

fn expect_err(got: Result<ExactNum, Error>, want: Error, what: &str) {
    match got {
        Err(e) => assert_eq!(e, want, "{what}: wrong error"),
        Ok(v) => assert!(false, "{what}: want Err({want:?}), got Ok({v:?})"),
    }
}

#[test]
fn dim_normal_and_negative() -> Result<(), Error> {
    let mut cc = Consts::new()?;
    let p = 128;
    assert_eq_num(
        &parse("5", p, &mut cc).dim(&parse("3", p, &mut cc), p, RoundingMode::ToEven)?,
        &parse("2", p, &mut cc),
        "dim(5,3)",
    );
    assert_pos_zero(
        &parse("3", p, &mut cc).dim(&parse("5", p, &mut cc), p, RoundingMode::ToEven)?,
        "dim(3,5)",
    );
    assert_pos_zero(
        &parse("3", p, &mut cc).dim(&parse("3", p, &mut cc), p, RoundingMode::ToEven)?,
        "dim(3,3)",
    );
    assert_eq_num(
        &parse("3", p, &mut cc).dim(&parse("-3", p, &mut cc), p, RoundingMode::ToEven)?,
        &parse("6", p, &mut cc),
        "dim(3,-3)",
    );
    assert_eq_num(
        &parse("-3", p, &mut cc).dim(&parse("-5", p, &mut cc), p, RoundingMode::ToEven)?,
        &parse("2", p, &mut cc),
        "dim(-3,-5)",
    );
    assert_pos_zero(
        &parse("-5", p, &mut cc).dim(&parse("-3", p, &mut cc), p, RoundingMode::ToEven)?,
        "dim(-5,-3)",
    );
    assert_pos_zero(
        &parse("-3", p, &mut cc).dim(&parse("3", p, &mut cc), p, RoundingMode::ToEven)?,
        "dim(-3,3)",
    );
    Ok(())
}

#[test]
fn dim_signed_zeros() -> Result<(), Error> {
    let mut cc = Consts::new()?;
    let p = 128;
    let pos_zero = ExactNum::new(p);
    let mut neg_zero = ExactNum::new(p);
    neg_zero.set_sign(Sign::Neg);
    assert_pos_zero(
        &pos_zero.dim(&neg_zero, p, RoundingMode::ToEven)?,
        "dim(+0,-0)",
    );
    assert_pos_zero(
        &neg_zero.dim(&pos_zero, p, RoundingMode::ToEven)?,
        "dim(-0,+0)",
    );
    assert_pos_zero(
        &neg_zero.dim(&neg_zero, p, RoundingMode::ToEven)?,
        "dim(-0,-0)",
    );
    assert_pos_zero(
        &pos_zero.dim(&pos_zero, p, RoundingMode::ToEven)?,
        "dim(+0,+0)",
    );
    Ok(())
}

#[test]
fn dim_infinities() -> Result<(), Error> {
    let mut cc = Consts::new()?;
    let p = 128;
    let one = parse("1", p, &mut cc);
    let pos_inf = ExactNum::from_i8(1, p).div(&ExactNum::new(p), p, RoundingMode::ToEven);
    assert!(pos_inf.is_inf_pos(), "setup: 1/+0 must be +Inf");
    let neg_inf = ExactNum::from_i8(-1, p).div(&ExactNum::new(p), p, RoundingMode::ToEven);
    assert!(neg_inf.is_inf_neg(), "setup: -1/+0 must be -Inf");
    assert_pos_zero(
        &pos_inf.dim(&pos_inf, p, RoundingMode::ToEven)?,
        "dim(+Inf,+Inf)",
    );
    assert!(
        pos_inf
            .dim(&one, p, RoundingMode::ToEven)?
            .is_inf_pos(),
        "dim(+Inf,1) must be +Inf"
    );
    assert!(
        one.dim(&neg_inf, p, RoundingMode::ToEven)?.is_inf_pos(),
        "dim(1,-Inf) must be +Inf"
    );
    assert_pos_zero(
        &neg_inf.dim(&one, p, RoundingMode::ToEven)?,
        "dim(-Inf,1)",
    );
    assert_pos_zero(
        &one.dim(&pos_inf, p, RoundingMode::ToEven)?,
        "dim(1,+Inf)",
    );
    assert_pos_zero(
        &neg_inf.dim(&neg_inf, p, RoundingMode::ToEven)?,
        "dim(-Inf,-Inf)",
    );
    assert!(
        pos_inf
            .dim(&neg_inf, p, RoundingMode::ToEven)?
            .is_inf_pos(),
        "dim(+Inf,-Inf) must be +Inf"
    );
    assert_pos_zero(
        &neg_inf.dim(&pos_inf, p, RoundingMode::ToEven)?,
        "dim(-Inf,+Inf)",
    );
    Ok(())
}

#[test]
fn dim_nan_is_error() -> Result<(), Error> {
    let mut cc = Consts::new()?;
    let p = 128;
    let one = parse("1", p, &mut cc);
    expect_err(
        ExactNum::nan(Some(Error::InvalidArgument)).dim(&one, p, RoundingMode::ToEven),
        Error::InvalidArgument,
        "dim(NaN,1)",
    );
    expect_err(
        one.dim(&ExactNum::nan(Some(Error::InvalidArgument)), p, RoundingMode::ToEven),
        Error::InvalidArgument,
        "dim(1,NaN)",
    );
    // A NaN carrying a different error propagates that error, per the pack (Err or propagate).
    expect_err(
        ExactNum::nan(Some(Error::DivisionByZero)).dim(&one, p, RoundingMode::ToEven),
        Error::DivisionByZero,
        "dim(NaN(DivisionByZero),1) propagates",
    );
    // A bare NaN (no carried error) is still refused, never Ok(NaN).
    expect_err(
        ExactNum::nan(None).dim(&one, p, RoundingMode::ToEven),
        Error::InvalidArgument,
        "dim(bare NaN,1)",
    );
    Ok(())
}

#[test]
fn dim_huge_overflow_follows_rounding_mode() -> Result<(), Error> {
    let p = 128;
    let max = ExactNum::max_value(p);
    let mut neg_max = ExactNum::max_value(p);
    neg_max.set_sign(Sign::Neg);
    expect_err(
        max.dim(&neg_max, p, RoundingMode::ToEven),
        Error::ExponentOverflow(Sign::Pos),
        "dim(max,-max) ToEven",
    );
    assert_eq_num(
        &max.dim(&neg_max, p, RoundingMode::ToZero)?,
        &ExactNum::max_value(p),
        "dim(max,-max) ToZero",
    );
    // EXPONENT_MAX is a real bound of the type, not a value dim can produce from finite
    // operands at p=128 without overflow; record it so the row is not silently dropped.
    assert_eq!(EXPONENT_MAX, (1i64 << 62) - 1, "EXPONENT_MAX sanity");
    Ok(())
}

#[test]
fn dim_tiny_values() -> Result<(), Error> {
    let p = 128;
    let tiny = ExactNum::pow2(-1000, p)?;
    let tinier = ExactNum::pow2(-1001, p)?;
    assert_eq_num(
        &tiny.dim(&tinier, p, RoundingMode::ToEven)?,
        &tinier,
        "dim(2^-1000,2^-1001)",
    );
    let min_pos = ExactNum::min_positive(p);
    let mut neg_min_pos = ExactNum::min_positive(p);
    neg_min_pos.set_sign(Sign::Neg);
    // 2*min_positive: built at p=128 via ldexp of min_positive by 1, exact.
    let want = ExactNum::min_positive(p).ldexp(1, p, RoundingMode::ToEven);
    assert_eq_num(
        &min_pos.dim(&neg_min_pos, p, RoundingMode::ToEven)?,
        &want,
        "dim(min_positive,-min_positive)",
    );
    let big = ExactNum::pow2(1000, p)?;
    let half_big = ExactNum::pow2(999, p)?;
    assert_eq_num(
        &big.dim(&half_big, p, RoundingMode::ToEven)?,
        &half_big,
        "dim(2^1000,2^999)",
    );
    Ok(())
}

#[test]
fn dim_rounding_and_precisions() -> Result<(), Error> {
    let mut cc = Consts::new()?;
    for p in [1usize, 53, 63, 64, 65, 113, 128, 200] {
        assert_eq_num(
            &parse("5", p, &mut cc).dim(&parse("3", p, &mut cc), p, RoundingMode::ToEven)?,
            &parse("2", p, &mut cc),
            "dim(5,3) exact at every precision",
        );
    }
    // 1 - 2^-200 at p=64.
    let p = 64;
    let one = parse("1", p, &mut cc);
    let tiny = ExactNum::pow2(-200, p)?;
    // 1-2^-64 built independently: 1 minus 2^-64 at p=128, then to p=64 (exact, 64 ones).
    let mut cc2 = Consts::new()?;
    let one_128 = parse("1", 128, &mut cc2);
    let tiny64_128 = ExactNum::pow2(-64, 128)?;
    let mut below_one = one_128.sub(&tiny64_128, 128, RoundingMode::ToEven);
    below_one.set_precision(64, RoundingMode::ToEven)?;
    assert!(!below_one.inexact(), "1-2^-64 must be exact at p=64");
    for rm in [RoundingMode::ToEven, RoundingMode::ToOdd, RoundingMode::Up, RoundingMode::FromZero] {
        assert_eq_num(
            &one.dim(&tiny, p, rm)?,
            &one,
            "dim(1,2^-200) rounds to 1",
        );
    }
    for rm in [RoundingMode::Down, RoundingMode::ToZero] {
        assert_eq_num(
            &one.dim(&tiny, p, rm)?,
            &below_one,
            "dim(1,2^-200) rounds down to 1-2^-64",
        );
    }
    // RoundingMode::None: sub keeps the exact difference, so dim does too.
    let exact = one.dim(&tiny, p, RoundingMode::None)?;
    assert!(!exact.inexact(), "dim with None must be exact");
    assert_eq_num(
        &exact.sub(&one, 200, RoundingMode::ToEven),
        &ExactNum::pow2(-200, 200)?.mul(&ExactNum::from_i8(-1, 200), 200, RoundingMode::ToEven),
        "dim(1,2^-200) None minus 1 is -2^-200",
    );
    Ok(())
}

#[test]
fn dim_invalid_precision_and_inexact_flag() -> Result<(), Error> {
    let mut cc = Consts::new()?;
    let one = parse("1", 128, &mut cc);
    let two = parse("2", 128, &mut cc);
    expect_err(
        two.dim(&one, 0, RoundingMode::ToEven),
        Error::InvalidArgument,
        "dim p=0",
    );
    expect_err(
        two.dim(&one, MAX_PRECISION + 1, RoundingMode::ToEven),
        Error::InvalidArgument,
        "dim p>MAX_PRECISION",
    );
    let exact = parse("5", 128, &mut cc).dim(&parse("3", 128, &mut cc), 128, RoundingMode::ToEven)?;
    assert!(!exact.inexact(), "dim(5,3) must be exact");
    let rounded = parse("1", 64, &mut cc).dim(&ExactNum::pow2(-200, 64)?, 64, RoundingMode::ToEven)?;
    assert!(rounded.inexact(), "dim(1,2^-200) at p=64 must be inexact");
    // Every mode gives the exact 2 for dim(5,3) and +0 for dim(3,5).
    for rm in MODES {
        assert_eq_num(
            &parse("5", 128, &mut cc).dim(&parse("3", 128, &mut cc), 128, rm)?,
            &parse("2", 128, &mut cc),
            "dim(5,3) in every mode",
        );
        assert_pos_zero(
            &parse("3", 128, &mut cc).dim(&parse("5", 128, &mut cc), 128, rm)?,
            "dim(3,5) in every mode",
        );
    }
    Ok(())
}
