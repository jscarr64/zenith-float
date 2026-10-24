//! R14 RawParts follow-ups (Spence): the precision a value was built at survives the round
//! trip, and a proptest over every kind checks `as_raw_parts -> from_raw_parts` is
//! bit-identical, both inexact flags included.

use proptest::prelude::*;
use std::ops::Neg;
use zenith_float_num::{
    ExactNum, Exponent, RawKind, RoundingMode, Sign, Word, EXPONENT_MAX, EXPONENT_MIN, INF_NEG,
    INF_POS, NAN, PROPTEST_CASES, WORD_BIT_SIZE,
};

const MODES: [RoundingMode; 5] = [
    RoundingMode::ToEven,
    RoundingMode::Up,
    RoundingMode::Down,
    RoundingMode::ToZero,
    RoundingMode::FromZero,
];

/// Everything observable about a value's representation and flags.
fn fingerprint(x: &ExactNum) -> String {
    let r = x.as_raw_parts();
    format!(
        "{:?} neg={} {:x?} prec={} e={} acc={} op={} precision()={:?} max_bits={:?} nan={}",
        r.kind,
        r.negative,
        r.words,
        r.precision,
        r.exponent,
        x.inexact_accumulated(),
        x.inexact(),
        x.precision(),
        x.mantissa_max_bit_len(),
        x.is_nan()
    )
}

fn rebuild(x: &ExactNum) -> Result<ExactNum, String> {
    ExactNum::from_raw_parts(&x.as_raw_parts())
        .map_err(|e| format!("own parts rejected: {e:?} for {}", fingerprint(x)))
}

fn check_identical(x: &ExactNum) -> Result<(), String> {
    let y = rebuild(x)?;
    let (a, b) = (fingerprint(x), fingerprint(&y));
    if a != b {
        return Err(format!("round trip changed the value:\n  {a}\n  {b}"));
    }
    Ok(())
}

/// `1/3` rounded at `p` bits (64-bit operands): inexact, with the low `round_p(p) - p` bits clear.
fn third(p: usize, rm: RoundingMode) -> ExactNum {
    ExactNum::from_i64(1, 64).div(&ExactNum::from_i64(3, 64), p, rm)
}

/// The representation keeps no requested precision apart from the word-rounded mantissa:
/// a value built at p = 100 stores 2 words (128 bits) and reports `precision == 128`. The
/// round trip gives the same representation, and every later operation behaves the same.
#[test]
fn precision_round_trip_p_1_63_65_100() -> Result<(), String> {
    for p in [1usize, 63, 64, 65, 100, 128, 129] {
        let words = p.div_ceil(WORD_BIT_SIZE);
        for rm in MODES {
            let x = third(p, rm);
            let parts = x.as_raw_parts();
            if parts.words.len() != words || parts.precision != words * WORD_BIT_SIZE {
                return Err(format!("p={p}: {}", fingerprint(&x)));
            }
            check_identical(&x)?;
            let y = rebuild(&x)?;
            // Same behavior afterwards: arithmetic at p and p+1 in every mode.
            for q in [p, p + 1] {
                for rm2 in MODES {
                    let seven = ExactNum::from_i64(7, 64);
                    let ops: [(
                        &str,
                        fn(ExactNum, &ExactNum, usize, RoundingMode) -> ExactNum,
                    ); 3] = [
                        ("add", |a, b, q, r| a.add(b, q, r)),
                        ("mul", |a, b, q, r| a.mul(b, q, r)),
                        ("div", |a, b, q, r| a.div(b, q, r)),
                    ];
                    for (name, f) in ops {
                        let (ra, rb) = (
                            f(rebuild(&x)?, &seven, q, rm2),
                            f(rebuild(&y)?, &seven, q, rm2),
                        );
                        if fingerprint(&ra) != fingerprint(&rb) {
                            return Err(format!("p={p} q={q} {name} {rm2:?} differs"));
                        }
                    }
                    let (sa, sb) = (rebuild(&x)?.sqrt(q, rm2), y.sqrt(q, rm2));
                    if fingerprint(&sa) != fingerprint(&sb) {
                        return Err(format!("p={p} q={q} sqrt {rm2:?} differs"));
                    }
                }
            }
        }
    }
    Ok(())
}

/// A value whose producing operation was exact but whose inputs were rounded: the
/// per-operation flag is false while the accumulated flag is true.
#[test]
fn round_trip_keeps_split_inexact_flags() -> Result<(), String> {
    for p in [1usize, 64, 100] {
        let x =
            third(p, RoundingMode::ToEven).mul(&ExactNum::from_i64(1, 64), p, RoundingMode::ToEven);
        if x.inexact() || !x.inexact_accumulated() {
            return Err(format!(
                "setup: expected op=false acc=true: {}",
                fingerprint(&x)
            ));
        }
        check_identical(&x)?;
    }
    Ok(())
}

/// Kind selector: 0..=3 finite (rounded quotient, optionally through an exact multiply),
/// 4 zero, 5 Inf, 6 NaN, 7..=8 subnormal at `EXPONENT_MIN`.
fn exponent() -> impl Strategy<Value = Exponent> {
    prop_oneof![
        Just(EXPONENT_MAX),
        Just(EXPONENT_MIN),
        EXPONENT_MIN..EXPONENT_MIN + 70,
        EXPONENT_MAX - 70..=EXPONENT_MAX,
        -300i64..300,
    ]
}

/// `(kind, n, d, e, mode, flag, w, neg)`; see [`exponent`] for the kind selector.
type Shape = (u8, i64, i64, Exponent, usize, bool, Word, bool);

fn build(s: &Shape, p: usize) -> ExactNum {
    let (kind, n, d, e, rm, flag, w, neg) = *s;
    let rm = MODES[rm];
    match kind {
        0..=3 => {
            let n = if n == 0 { 1 } else { n };
            let mut x = ExactNum::from_i64(n, 64).div(&ExactNum::from_i64(d, 64), p, rm);
            if flag {
                x = x.mul(&ExactNum::from_i64(1, 64), p, rm);
            }
            if x.exponent().is_some() {
                x.set_exponent(e);
            }
            x
        }
        4 => {
            let z = ExactNum::new(p);
            if neg {
                z.neg()
            } else {
                z
            }
        }
        5 => {
            if neg {
                INF_NEG
            } else {
                INF_POS
            }
        }
        6 => NAN,
        _ => {
            let words = p.div_ceil(WORD_BIT_SIZE);
            let mut m = vec![0 as Word; words];
            if let Some(top) = m.last_mut() {
                *top = w;
            }
            let x = ExactNum::from_words(&m, if neg { Sign::Neg } else { Sign::Pos }, EXPONENT_MIN);
            if flag {
                // Halving at p rounds a subnormal whose low bit is set.
                x.div(&ExactNum::from_i64(2, 64), p, rm)
            } else {
                x
            }
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(PROPTEST_CASES))]

    #[test]
    fn round_trip_every_kind_bit_identical(
        s in (0u8..=8, any::<i64>(), 1i64..1_000_000, exponent(), 0usize..5, any::<bool>(), 1..Word::MAX >> 1, any::<bool>()),
        p in 1usize..=4096,
    ) {
        let x = build(&s, p);
        if s.0 >= 7 && !x.is_zero() {
            prop_assert!(x.as_raw_parts().kind == RawKind::Finite);
        }
        if let Err(e) = check_identical(&x) {
            return Err(TestCaseError::fail(format!("{s:?} p={p}: {e}")));
        }
    }
}
