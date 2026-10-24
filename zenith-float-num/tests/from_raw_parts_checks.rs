//! `ExactNum::from_raw_parts` enforces its documented checks (2.0.0 audit, defect 5; migrated
//! to the `RawParts` API in R14).

use std::ops::Neg;
use zenith_float_num::{
    Error, ExactNum, RawKind, RawParts, Word, EXPONENT_MAX, EXPONENT_MIN, WORD_BIT_SIZE,
};

fn parts(words: &[Word], negative: bool, exponent: i64) -> RawParts<'_> {
    RawParts {
        kind: RawKind::Finite,
        words,
        precision: words.len() * WORD_BIT_SIZE,
        negative,
        exponent,
        inexact: false,
        op_inexact: false,
    }
}

fn assert_invalid(p: &RawParts<'_>, what: &str) {
    assert!(
        matches!(ExactNum::from_raw_parts(p), Err(Error::InvalidArgument)),
        "Expected InvalidArgument for {what}"
    );
}

fn round_trip(x: &ExactNum) -> Result<ExactNum, Error> {
    ExactNum::from_raw_parts(&x.as_raw_parts())
}

#[test]
fn from_raw_parts_round_trips_valid_parts() -> Result<(), Error> {
    let ps = [64, 128, 256, 1024];
    let vs = [1u64, 3, 0xFFFF_FFFF_FFFF_FFFF, 0x8000_0000_0000_0001];
    let es = [EXPONENT_MIN, EXPONENT_MIN + 1, -1, 0, 1, EXPONENT_MAX - 1, EXPONENT_MAX];
    for &p in &ps {
        for &v in &vs {
            for &e in &es {
                for negative in [false, true] {
                    let mut x = ExactNum::from_u64(v, p);
                    x.set_exponent(e);
                    if negative {
                        x = x.neg();
                    }
                    let y = round_trip(&x)?;
                    assert_eq!(x, y, "p={p}, v={v}, e={e}, negative={negative}");
                }
            }
        }
        assert_eq!(
            round_trip(&ExactNum::min_positive(p))?,
            ExactNum::min_positive(p)
        );
        assert_eq!(
            round_trip(&ExactNum::min_positive(p).neg())?,
            ExactNum::min_positive(p).neg()
        );
    }
    Ok(())
}

#[test]
fn from_raw_parts_rejects_inconsistent_parts() -> Result<(), Error> {
    for p in [64usize, 128, 256, 1024] {
        let x = ExactNum::from_u64(3, p);
        let good = x.as_raw_parts();
        // (a) precision one more / one less than the stored width.
        for precision in [p + 1, p - 1] {
            let bad = RawParts {
                precision,
                ..parts(good.words, false, good.exponent)
            };
            assert_invalid(&bad, &format!("precision {precision} for p={p}"));
        }
        // (b)-(c) a subnormal mantissa (lowest bit only) needs e == EXPONENT_MIN.
        let mut m: Vec<Word> = vec![0; p / WORD_BIT_SIZE];
        if let Some(w) = m.first_mut() {
            *w = 1;
        }
        for e in [EXPONENT_MIN + 1, -1, 0, 1, EXPONENT_MAX] {
            for negative in [false, true] {
                assert_invalid(&parts(&m, negative, e), &format!("subnormal e={e} p={p}"));
            }
        }
        for negative in [false, true] {
            let y = ExactNum::from_raw_parts(&parts(&m, negative, EXPONENT_MIN))?;
            let want = if negative {
                ExactNum::min_positive(p).neg()
            } else {
                ExactNum::min_positive(p)
            };
            assert_eq!(y, want, "subnormal p={p} negative={negative}");
        }
        // (d) all-zero words are not a finite non-zero value.
        let z: Vec<Word> = vec![0; p / WORD_BIT_SIZE];
        for negative in [false, true] {
            assert_invalid(
                &parts(&z, negative, EXPONENT_MIN),
                &format!("all-zero p={p}"),
            );
        }
    }
    Ok(())
}
