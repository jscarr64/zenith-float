//! Step-1 properties: signed zeros, strict parsing, and decimal round-trip.
//! These sit beside the checked-in JSONL rows and do not replace them.

use proptest::prelude::*;
use zenith_float_num::{Consts, Error, ExactNum, Radix, RoundingMode, Sign};

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
        assert_eq!(n.cmp(&back), Some(0), "{rm:?} {spelled}");
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
        let mut pos = ExactNum::new(53);
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
        prop_assert_eq!(n.cmp(&back), Some(0));
        prop_assert_eq!(n.sign(), back.sign());
    }
}
