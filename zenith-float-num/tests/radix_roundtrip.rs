//! Structured stress tests for parse/format round-trip across radices 2–36.

use zenith_float_num::{
    Consts, DEFAULT_RANDOM_SEED, Error, ExactNum, Exponent, Radix, RoundingMode, Sign,
    WORD_BIT_SIZE, WORD_SIGNIFICANT_BIT, Word, seeded_random as random,
};

const LOOPS_PER_BASE: usize = 24;
const TEST_EXP_BOUND: Exponent = 64;

fn roundtrip_exact(
    rdx: Radix,
    n: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Result<(), Box<dyn std::error::Error>> {
    let s = n.format(rdx, rm, cc)?;
    let g = ExactNum::parse(&s, rdx, p, rm, cc);
    assert_eq!(n.compare(&g), Some(0), "rdx={:?} s={s}", rdx.value());
    Ok(())
}

/// A significand character is legal when it is `.` or a digit of `base`.
fn literal_in_radix(lit: &str, base: u8) -> bool {
    lit.chars().all(|c| c == '.' || c.is_digit(u32::from(base)))
}

/// In-radix text must round-trip as a finite value. A digit outside the radix is
/// `NaN(InvalidArgument)`; formatting it prints `NaN`, which parses back as NaN.
fn roundtrip_literal(
    rdx: Radix,
    base: u8,
    lit: &str,
    p: usize,
    cc: &mut Consts,
) -> Result<(), Box<dyn std::error::Error>> {
    let rm = RoundingMode::ToEven;
    let n = ExactNum::parse(lit, rdx, p, rm, cc);
    if literal_in_radix(lit, base) {
        roundtrip_exact(rdx, &n, p, rm, cc)?;
        return Ok(());
    }
    assert!(n.is_nan(), "rdx={} lit={lit}", base);
    assert_eq!(
        n.err(),
        Some(Error::InvalidArgument),
        "rdx={} lit={lit}",
        base
    );
    let s = n.format(rdx, rm, cc)?;
    assert_eq!(s, "NaN", "rdx={base} lit={lit}");
    let g = ExactNum::parse(&s, rdx, p, rm, cc);
    assert!(g.is_nan(), "rdx={base} s={s}");
    Ok(())
}

#[test]
fn radix_roundtrip_commensurable() -> Result<(), Box<dyn std::error::Error>> {
    zenith_float_num::reseed_random(DEFAULT_RANDOM_SEED);
    let mut cc = Consts::new()?;

    for rdx in [Radix::Bin, Radix::Oct, Radix::Hex] {
        for _ in 0..LOOPS_PER_BASE {
            let p = (random::<usize>() % 6 + 2) * WORD_BIT_SIZE;
            let n = ExactNum::random_normal(p, -TEST_EXP_BOUND, TEST_EXP_BOUND);
            roundtrip_exact(rdx, &n, p, RoundingMode::None, &mut cc)?;
        }
    }
    Ok(())
}

#[test]
fn radix_roundtrip_general_bases() -> Result<(), Box<dyn std::error::Error>> {
    let mut cc = Consts::new()?;
    let p = 192;
    for base in [3u8, 5, 7, 9, 11, 12, 15, 20, 36] {
        let rdx = Radix::try_new(base)?;
        for lit in ["0.1", "1.0", "2.5"] {
            roundtrip_literal(rdx, base, lit, p, &mut cc)?;
        }
        if base >= 11 {
            roundtrip_literal(rdx, base, "A.5", p, &mut cc)?;
        }
    }
    Ok(())
}

#[test]
fn radix_roundtrip_all_bases_smoke() -> Result<(), Box<dyn std::error::Error>> {
    let mut cc = Consts::new()?;
    let p = 128;
    let n = ExactNum::parse("1.25", Radix::Dec, p, RoundingMode::ToEven, &mut cc);

    for base in 2u8..=36 {
        roundtrip_exact(Radix::try_new(base)?, &n, p, RoundingMode::ToEven, &mut cc)?;
    }
    Ok(())
}

#[test]
fn radix_roundtrip_extremes() -> Result<(), Box<dyn std::error::Error>> {
    let mut cc = Consts::new()?;
    let p = 4 * WORD_BIT_SIZE;
    let rm = RoundingMode::None;

    for rdx in [Radix::Bin, Radix::Oct, Radix::Hex] {
        for n in [
            ExactNum::max_value(p),
            ExactNum::min_value(p),
            ExactNum::min_positive_normal(p),
            ExactNum::from_words(&[WORD_SIGNIFICANT_BIT], Sign::Pos, 0),
        ] {
            roundtrip_exact(rdx, &n, p, rm, &mut cc)?;
        }
    }
    Ok(())
}

/// Official W1 rows: format must parse back in the same mode. 16 extra digits
/// were not enough (and for general radix were added to the bit width, not the
/// digit count). The digit count is MPFR's `1 + ceil(p·log(2)/log(b))`.
#[test]
fn zf_format_minpos_maxv_and_dense_round_trip() -> Result<(), Box<dyn std::error::Error>> {
    let mut cc = Consts::new()?;
    let modes = [
        RoundingMode::Up,
        RoundingMode::Down,
        RoundingMode::ToZero,
        RoundingMode::FromZero,
        RoundingMode::ToEven,
        RoundingMode::ToOdd,
        RoundingMode::None,
    ];

    // minpos:128 radix 10 — official WRONG in Up, FromZero, ToEven, ToOdd.
    let minpos = ExactNum::min_positive(128);
    let dec = Radix::Dec;
    for rm in [
        RoundingMode::Up,
        RoundingMode::FromZero,
        RoundingMode::ToEven,
        RoundingMode::ToOdd,
    ] {
        let s = minpos.format(dec, rm, &mut cc)?;
        let back = ExactNum::parse(&s, dec, 128, rm, &mut cc);
        assert!(
            !back.is_nan() && minpos.compare(&back) == Some(0) && minpos.sign() == back.sign(),
            "min_positive(128) radix 10 {rm:?}: format {s} parsed back to {back:?}"
        );
    }

    // maxv:1024 radix 3 — official WRONG in every mode.
    let maxv = ExactNum::max_value(1024);
    let r3 = Radix::try_new(3)?;
    for rm in modes {
        let s = maxv.format(r3, rm, &mut cc)?;
        let back = ExactNum::parse(&s, r3, 1024, rm, &mut cc);
        assert!(
            !back.is_nan() && maxv.compare(&back) == Some(0) && maxv.sign() == back.sign(),
            "max_value(1024) radix 3 {rm:?}: format {s:?} ({} chars) parsed back to {back:?}",
            s.len()
        );
    }

    // Dense 1024-bit mantissa, radix 3 (official c0576495): hex a..ab times 2^{-1025}
    // (harness token `...abp-1025`), i.e. that significand at exponent -1.
    let mut words = [0xaaaaaaaaaaaaaaaa as Word; 16];
    words[0] = 0xaaaaaaaaaaaaaaab as Word;
    let n = ExactNum::from_words(&words, Sign::Pos, -1);
    assert!(!n.is_nan(), "dense a..abp-1025 construction");
    for rm in [RoundingMode::Down, RoundingMode::ToZero, RoundingMode::None] {
        let s = n.format(r3, rm, &mut cc)?;
        let back = ExactNum::parse(&s, r3, 1024, rm, &mut cc);
        assert!(
            !back.is_nan() && n.compare(&back) == Some(0) && n.sign() == back.sign(),
            "dense a..abp-1025 radix 3 {rm:?}: format ({} chars) did not round-trip",
            s.len()
        );
    }
    Ok(())
}

/// Official leftovers after a09e84b: `r.format` of minpos in radices 3 and 36
/// must return a spelling that parses back, never `Err`. That includes p=256
/// Up/FromZero (official 1 ulp) and Down/ToZero/None (CORRECT at fb66cec, then
/// `S:Err`). `convert_to_radix` radix 10 of minpos must round-trip too.
#[test]
fn zf_format_minpos_radix3_36_and_radix_conv_round_trip() -> Result<(), Box<dyn std::error::Error>>
{
    let mut cc = Consts::new()?;
    let modes = [
        RoundingMode::Up,
        RoundingMode::Down,
        RoundingMode::ToZero,
        RoundingMode::FromZero,
        RoundingMode::ToEven,
        RoundingMode::ToOdd,
        RoundingMode::None,
    ];
    for p in [256usize, 1024] {
        let minpos = ExactNum::min_positive(p);
        for &base in &[3u8, 36] {
            let rdx = Radix::try_new(base)?;
            for rm in modes {
                let s = minpos.format(rdx, rm, &mut cc)?;
                assert!(
                    s != "Err" && !s.ends_with("Err") && !s.contains("Err"),
                    "min_positive({p}) radix {base} {rm:?} gave {s:?}"
                );
                let back = ExactNum::parse(&s, rdx, p, rm, &mut cc);
                assert!(
                    !back.is_nan()
                        && minpos.compare(&back) == Some(0)
                        && minpos.sign() == back.sign(),
                    "min_positive({p}) radix {base} {rm:?}: format {s:?} parsed back to {back:?}"
                );
            }
        }
    }
    for p in [128usize, 2048] {
        let minpos = ExactNum::min_positive(p);
        let (sgn, digits, exp) =
            minpos.convert_to_radix(Radix::Dec, RoundingMode::ToEven, &mut cc)?;
        let back = ExactNum::convert_from_radix(
            sgn,
            &digits,
            exp,
            Radix::Dec,
            p,
            RoundingMode::ToEven,
            &mut cc,
        );
        assert!(
            !back.is_nan() && minpos.compare(&back) == Some(0) && minpos.sign() == back.sign(),
            "convert_to_radix min_positive({p}) radix 10 ToEven did not round-trip"
        );
    }
    Ok(())
}
