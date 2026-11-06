//! ZF-TD41: CRLF input handling for the CSV parsers and the goldens.
//!
//! Guard list (no new function, so no A-Z edge-case list):
//! - `Ieee64Array`: a CRLF text parses to the same array as the LF text
//!   (shape and bit patterns), including the last field of each row.
//! - `ExactNumArray`: the same, via a `to_csv_string` round-trip and via a
//!   literal CRLF text, including the last field of each row.
//! - Every golden file the tests load with `include_str!` contains no `\r`
//!   on this checkout.
//!
//! Honesty note, checked against current main (49a358a): the parser part of
//! this file is a GUARD, not a bug reproduction. `data_lines` in
//! `zenith-float-num/src/csvfmt.rs` already does
//! `.split('\n').map(|l| l.strip_suffix('\r').unwrap_or(l))`, and it has
//! done so since the original CSV commit 42eac59, and both parsers also
//! trim cells. So the parser tests below PASS on main. The part of ZF-TD41
//! that actually fails on main is the `.gitattributes` protection, covered
//! by `tools/gitattributes_selftest.sh` (there is no `.gitattributes` on
//! main), and the golden-sums CRLF hazard, covered by the CRLF case in
//! `tools/golden_sums_selftest.sh`. This PR therefore changes no parser
//! code and no golden bytes; it adds the protection and the guards.

use zenith_float_num::{Error, ExactNum, ExactNumArray, Ieee64, Ieee64Array};

const BORROW_OPS: &str = include_str!("oracle/goldens/borrow_ops.jsonl");
const LOG_EXTREME: &str = include_str!("oracle/goldens/log_extreme.tsv");
const R13_DIV_HYPOT: &str = include_str!("oracle/goldens/r13_div_hypot.jsonl");
const R14_CONSTS: &str = include_str!("oracle/goldens/r14_consts.jsonl");
const R9_OVERFLOW: &str = include_str!("oracle/goldens/r9_overflow.jsonl");
const S3_EULER_GAMMA: &str = include_str!("oracle/goldens/s3_euler_gamma.jsonl");
const S9_SPECIALS: &str = include_str!("oracle/goldens/s9_specials.jsonl");

#[test]
fn crlf_ieee64_array_parses_like_lf() -> Result<(), Error> {
    let a = Ieee64Array::from_shape(
        2,
        2,
        &[
            Ieee64::from_i32(0),
            Ieee64::from_i32(1),
            Ieee64::from_i32(2),
            Ieee64::from_i32(3),
        ],
    )
    .ok_or(Error::InvalidArgument)?;
    let lf = a.to_csv_string()?;
    let crlf = lf.replace('\n', "\r\n");
    let b = Ieee64Array::from_csv_str(&crlf)?;
    assert_eq!(b.shape(), (2, 2));
    assert_eq!(a.as_bits(), b.as_bits());
    // The last field of each row, directly: bit patterns 2 and 4, with no
    // stray '\r' surviving into the parse.
    let c = Ieee64Array::from_csv_str("1,2\r\n3,4\r\n")?;
    assert_eq!(c.shape(), (2, 2));
    assert_eq!(c.get2(0, 1).ok_or(Error::InvalidArgument)?.to_bits(), 2);
    assert_eq!(c.get2(1, 1).ok_or(Error::InvalidArgument)?.to_bits(), 4);
    Ok(())
}

#[test]
fn crlf_exact_num_array_parses_like_lf() -> Result<(), Error> {
    let p = 64;
    let a = ExactNumArray::from_shape(
        p,
        2,
        2,
        &[
            ExactNum::from_u8(1, p),
            ExactNum::from_u8(2, p),
            ExactNum::from_u8(3, p),
            ExactNum::from_u8(4, p),
        ],
    )
    .ok_or(Error::InvalidArgument)?;
    let lf = a.to_csv_string()?;
    let crlf = lf.replace('\n', "\r\n");
    let b = ExactNumArray::from_csv_str(&crlf)?;
    assert_eq!(b.shape(), (2, 2));
    assert_eq!(b.precision(), p);
    for i in 0..4 {
        assert_eq!(
            a.get(i)
                .ok_or(Error::InvalidArgument)?
                .compare(b.get(i).ok_or(Error::InvalidArgument)?),
            Some(0)
        );
    }
    // The last field of each row, directly.
    let c = ExactNumArray::from_csv_str("1@p=64,2@p=64\r\n3@p=64,4@p=64\r\n")?;
    assert_eq!(c.shape(), (2, 2));
    assert_eq!(
        c.get(1)
            .ok_or(Error::InvalidArgument)?
            .compare(&ExactNum::from_u8(2, p)),
        Some(0)
    );
    assert_eq!(
        c.get(3)
            .ok_or(Error::InvalidArgument)?
            .compare(&ExactNum::from_u8(4, p)),
        Some(0)
    );
    Ok(())
}

#[test]
fn goldens_have_no_carriage_returns() {
    assert!(!BORROW_OPS.contains('\r'), "borrow_ops.jsonl contains CR");
    assert!(!LOG_EXTREME.contains('\r'), "log_extreme.tsv contains CR");
    assert!(
        !R13_DIV_HYPOT.contains('\r'),
        "r13_div_hypot.jsonl contains CR"
    );
    assert!(!R14_CONSTS.contains('\r'), "r14_consts.jsonl contains CR");
    assert!(!R9_OVERFLOW.contains('\r'), "r9_overflow.jsonl contains CR");
    assert!(
        !S3_EULER_GAMMA.contains('\r'),
        "s3_euler_gamma.jsonl contains CR"
    );
    assert!(!S9_SPECIALS.contains('\r'), "s9_specials.jsonl contains CR");
}
