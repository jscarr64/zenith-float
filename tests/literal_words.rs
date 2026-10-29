//! `exact!` literals carry host-independent mantissa words: the proc macro runs on the host, so the
//! expansion must compile and round-trip exactly on a 32-bit target too.

use zenith_float::exact;
use zenith_float_num::{Consts, ExactNum, Radix, RoundingMode};

/// The literal's words, precision, sign and exponent equal a runtime parse at the same precision.
fn check(lit: &ExactNum, s: &str) -> Result<(), Box<dyn std::error::Error>> {
    let mut cc = Consts::new()?;
    let lp = lit.as_raw_parts();
    let run = ExactNum::parse(s, Radix::Dec, lp.precision, RoundingMode::ToEven, &mut cc);
    assert!(!run.inexact(), "{s}");
    let rp = run.as_raw_parts();
    assert_eq!(lp.words, rp.words, "{s}");
    assert_eq!(
        (lp.precision, lp.negative, lp.exponent, lp.inexact),
        (rp.precision, rp.negative, rp.exponent, rp.inexact),
        "{s}"
    );
    assert_eq!(lit.compare(&run), Some(0), "{s}");
    Ok(())
}

#[test]
fn single_word_literal_round_trips() -> Result<(), Box<dyn std::error::Error>> {
    check(&exact!("2.5"), "2.5")
}

/// 2^64 + 1 and a three-64-bit-word value: every 32-bit half, including the inner zero halves,
/// must survive.
#[test]
fn multi_word_literal_round_trips() -> Result<(), Box<dyn std::error::Error>> {
    check(&exact!("18446744073709551617"), "18446744073709551617")?;
    check(
        &exact!("-3138550867693340381917894711603833208051177722232017256449.0625"),
        "-3138550867693340381917894711603833208051177722232017256449.0625",
    )
}
