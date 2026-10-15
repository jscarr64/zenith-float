//! Properties beside the checked-in JSONL rows: signed zeros, strict parsing,
//! decimal round-trip, elementary symmetries, and constant-cache order.

use proptest::prelude::*;
use std::ops::Neg;
use zenith_float_num::{Consts, Error, ExactNum, Radix, RoundingMode, Sign};

/// `-x` for an operand that stays borrowed: fresh `-1` times `x` at `x`'s own precision.
fn neg_of(x: &ExactNum) -> ExactNum {
    let p = x
        .mantissa_max_bit_len()
        .unwrap_or(zenith_float_num::WORD_BIT_SIZE);
    ExactNum::from_i8(-1, p).mul(x, p, RoundingMode::None)
}

fn cc() -> Consts {
    match Consts::new() {
        Ok(c) => c,
        Err(_) => panic!("constants"),
    }
}

fn is_invalid(n: &ExactNum) -> bool {
    n.is_nan() && n.err() == Some(Error::InvalidArgument)
}

#[test]
fn one_tenth_round_trips_directed() {
    for rm in [RoundingMode::Up, RoundingMode::ToZero, RoundingMode::ToEven] {
        let mut consts = cc();
        let n = ExactNum::parse("0.1", Radix::Dec, 53, rm, &mut consts);
        let spelled = match n.format(Radix::Dec, rm, &mut consts) {
            Ok(s) => s,
            Err(err) => panic!("format {err:?}"),
        };
        let back = ExactNum::parse(&spelled, Radix::Dec, 53, rm, &mut consts);
        assert_eq!(n.compare(&back), Some(0), "{rm:?} {spelled}");
        assert_eq!(n.sign(), back.sign());
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(48))]

    #[test]
    fn parser_does_not_panic(bytes in proptest::collection::vec(any::<u8>(), 0..48)) {
        let mut text = String::new();
        for b in bytes {
            if b.is_ascii() && b != 0 {
                text.push(b as char);
            }
        }
        let mut consts = cc();
        let n = ExactNum::parse(&text, Radix::Dec, 53, RoundingMode::ToEven, &mut consts);
        let _ = (n.is_nan(), n.is_inf(), n.is_zero(), n.err());
    }

    #[test]
    fn trailing_mark_is_invalid(prefix in "[0-9]{1,6}", junk in "[A-DF-Za-df-z!]{1,4}") {
        let mut consts = cc();
        let text = format!("{prefix}{junk}");
        let n = ExactNum::parse(&text, Radix::Dec, 53, RoundingMode::ToEven, &mut consts);
        prop_assert!(is_invalid(&n), "{text} -> {:?}", n.err());
    }

    #[test]
    fn zero_sum_sign(down in proptest::bool::ANY) {
        let rm = if down { RoundingMode::Down } else { RoundingMode::ToEven };
        let pos = ExactNum::new(53);
        let mut neg = ExactNum::new(53);
        neg.set_sign(Sign::Neg);
        let sum = pos.add(&neg, 53, rm);
        prop_assert!(sum.is_zero());
        let want = if down { Sign::Neg } else { Sign::Pos };
        prop_assert_eq!(sum.sign(), Some(want));
    }

    #[test]
    fn decimal_round_trip(bits in 0u64..10_000, mode in 0u8..7) {
        let rm = match mode {
            0 => RoundingMode::ToEven,
            1 => RoundingMode::ToOdd,
            2 => RoundingMode::Up,
            3 => RoundingMode::Down,
            4 => RoundingMode::ToZero,
            5 => RoundingMode::FromZero,
            _ => RoundingMode::None,
        };
        let mut consts = cc();
        let text = format!("{}.5", bits);
        let n = ExactNum::parse(&text, Radix::Dec, 53, rm, &mut consts);
        prop_assume!(!n.is_nan() && !n.is_inf());
        let spelled = match n.format(Radix::Dec, rm, &mut consts) {
            Ok(s) => s,
            Err(_) => return Err(TestCaseError::fail("format failed")),
        };
        let back = ExactNum::parse(&spelled, Radix::Dec, 53, rm, &mut consts);
        prop_assert_eq!(n.compare(&back), Some(0));
        prop_assert_eq!(n.sign(), back.sign());
    }
}

fn same_bits(a: &ExactNum, b: &ExactNum) -> bool {
    a.sign() == b.sign()
        && a.exponent() == b.exponent()
        && a.mantissa_digits() == b.mantissa_digits()
}

#[test]
fn elementary_symmetries() {
    let mut consts = cc();
    let p = 53;
    let rm = RoundingMode::ToEven;
    for text in ["0.5", "1", "1.25", "2"] {
        let x = ExactNum::parse(text, Radix::Dec, p, rm, &mut consts);
        let neg_x = neg_of(&x);
        let erf_x = x.erf(p, rm, &mut consts);
        let erf_neg = neg_x.erf(p, rm, &mut consts);
        assert!(
            same_bits(&erf_x.neg(), &erf_neg),
            "erf(-{text}) is not -erf({text})"
        );
        let cos_x = x.cos(p, rm, &mut consts);
        let cos_neg = neg_x.cos(p, rm, &mut consts);
        assert!(
            same_bits(&cos_x, &cos_neg),
            "cos(-{text}) is not cos({text})"
        );
        let sinh_x = x.sinh(p, rm, &mut consts);
        let sinh_neg = neg_x.sinh(p, rm, &mut consts);
        assert!(
            same_bits(&sinh_x.neg(), &sinh_neg),
            "sinh(-{text}) is not -sinh({text})"
        );
    }
}

#[test]
fn sqrt2_and_phi_match_a_fresh_cache() {
    let directed = [RoundingMode::Up, RoundingMode::FromZero];
    for warm in [24usize, 53, 64, 256] {
        for rm in directed {
            let mut fresh = cc();
            let mut warmed = cc();
            let _ = warmed.sqrt2(warm, RoundingMode::ToEven);
            let cached = warmed.sqrt2(113, rm);
            let direct = fresh.sqrt2(113, rm);
            assert!(
                same_bits(&cached, &direct),
                "sqrt2(113, {rm:?}) after {warm} differs from a fresh cache"
            );
            let mut fresh = cc();
            let mut warmed = cc();
            let _ = warmed.phi(warm, RoundingMode::ToEven);
            let cached = warmed.phi(113, rm);
            let direct = fresh.phi(113, rm);
            assert!(
                same_bits(&cached, &direct),
                "phi(113, {rm:?}) after {warm} differs from a fresh cache"
            );
        }
    }
}
