//! Structured stress tests for parse/format round-trip across radices 2–36.

use zenith_float_num::{
    seeded_random as random, Consts, Error, ExactNum, Exponent, Radix, RoundingMode, Sign,
    DEFAULT_RANDOM_SEED, WORD_BIT_SIZE, WORD_SIGNIFICANT_BIT,
};

const LOOPS_PER_BASE: usize = 24;
const TEST_EXP_BOUND: Exponent = 64;

fn roundtrip_exact(rdx: Radix, n: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) {
    let s = n.format(rdx, rm, cc).expect("format");
    let g = ExactNum::parse(&s, rdx, p, rm, cc);
    assert_eq!(n.cmp(&g), Some(0), "rdx={:?} s={s}", rdx.value());
}

/// A significand character is legal when it is `.` or a digit of `base`.
fn literal_in_radix(lit: &str, base: u8) -> bool {
    lit.chars().all(|c| c == '.' || c.is_digit(u32::from(base)))
}

/// In-radix text must round-trip as a finite value. A digit outside the radix is
/// `NaN(InvalidArgument)`; formatting it prints `NaN`, which parses back as NaN.
fn roundtrip_literal(rdx: Radix, base: u8, lit: &str, p: usize, cc: &mut Consts) {
    let rm = RoundingMode::ToEven;
    let n = ExactNum::parse(lit, rdx, p, rm, cc);
    if literal_in_radix(lit, base) {
        roundtrip_exact(rdx, &n, p, rm, cc);
        return;
    }
    assert!(n.is_nan(), "rdx={} lit={lit}", base);
    assert_eq!(
        n.err(),
        Some(Error::InvalidArgument),
        "rdx={} lit={lit}",
        base
    );
    let s = n.format(rdx, rm, cc).expect("format");
    assert_eq!(s, "NaN", "rdx={base} lit={lit}");
    let g = ExactNum::parse(&s, rdx, p, rm, cc);
    assert!(g.is_nan(), "rdx={base} s={s}");
}

#[test]
fn radix_roundtrip_commensurable() {
    zenith_float_num::reseed_random(DEFAULT_RANDOM_SEED);
    let mut cc = Consts::new().unwrap();

    for rdx in [Radix::Bin, Radix::Oct, Radix::Hex] {
        for _ in 0..LOOPS_PER_BASE {
            let p = (random::<usize>() % 6 + 2) * WORD_BIT_SIZE;
            let n = ExactNum::random_normal(p, -TEST_EXP_BOUND, TEST_EXP_BOUND);
            roundtrip_exact(rdx, &n, p, RoundingMode::None, &mut cc);
        }
    }
}

#[test]
fn radix_roundtrip_general_bases() {
    let mut cc = Consts::new().unwrap();
    let p = 192;

    for base in [3u8, 5, 7, 9, 11, 12, 15, 20, 36] {
        let rdx = Radix::try_new(base).unwrap();
        for lit in ["0.1", "1.0", "2.5"] {
            roundtrip_literal(rdx, base, lit, p, &mut cc);
        }
        if base >= 11 {
            roundtrip_literal(rdx, base, "A.5", p, &mut cc);
        }
    }
}

#[test]
fn radix_roundtrip_all_bases_smoke() {
    let mut cc = Consts::new().unwrap();
    let p = 128;
    let n = ExactNum::parse("1.25", Radix::Dec, p, RoundingMode::ToEven, &mut cc);

    for base in 2u8..=36 {
        roundtrip_exact(
            Radix::try_new(base).unwrap(),
            &n,
            p,
            RoundingMode::ToEven,
            &mut cc,
        );
    }
}

#[test]
fn radix_roundtrip_extremes() {
    let mut cc = Consts::new().unwrap();
    let p = 4 * WORD_BIT_SIZE;
    let rm = RoundingMode::None;

    for rdx in [Radix::Bin, Radix::Oct, Radix::Hex] {
        for n in [
            ExactNum::max_value(p),
            ExactNum::min_value(p),
            ExactNum::min_positive_normal(p),
            ExactNum::from_words(&[WORD_SIGNIFICANT_BIT], Sign::Pos, 0),
        ] {
            roundtrip_exact(rdx, &n, p, rm, &mut cc);
        }
    }
}
