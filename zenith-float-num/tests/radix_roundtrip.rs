//! Structured stress tests for parse/format round-trip across radices 2–36.

use zenith_float_num::{
    seeded_random as random, Consts, ExactNum, Exponent, Radix, RoundingMode, Sign,
    DEFAULT_RANDOM_SEED, WORD_BIT_SIZE, WORD_SIGNIFICANT_BIT,
};

const LOOPS_PER_BASE: usize = 24;
const TEST_EXP_BOUND: Exponent = 64;

fn roundtrip_exact(rdx: Radix, n: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) {
    let s = n.format(rdx, rm, cc).expect("format");
    let g = ExactNum::parse(&s, rdx, p, rm, cc);
    assert_eq!(n.cmp(&g), Some(0), "rdx={:?} s={s}", rdx.value());
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
            if base < 11 && lit.chars().any(|c| c.is_ascii_alphabetic()) {
                continue;
            }
            let n = ExactNum::parse(lit, rdx, p, RoundingMode::ToEven, &mut cc);
            roundtrip_exact(rdx, &n, p, RoundingMode::ToEven, &mut cc);
        }
        if base >= 11 {
            let n = ExactNum::parse("A.5", rdx, p, RoundingMode::ToEven, &mut cc);
            roundtrip_exact(rdx, &n, p, RoundingMode::ToEven, &mut cc);
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
