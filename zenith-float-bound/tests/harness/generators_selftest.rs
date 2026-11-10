//! COVERS: common::BallDesc::encode, common::ExpBand::Edge, common::ExpBand::Gap, common::ExpBand::Near0, common::ExpBand::PowerOfTwo, common::ExpBand::Uniform, common::ExpBand::from_index, common::ExpBand::index, common::ExpDomain::edge_values, common::ExpDomain::in_raw_range, common::GapClass::One, common::GapClass::P, common::GapClass::PGuard, common::GapClass::PGuardMinus1, common::GapClass::PGuardPlus1, common::GapClass::PMinus1, common::GapClass::PPlus1, common::GapClass::Random, common::GapClass::Two, common::GapClass::Zero, common::GapClass::from_index, common::GapClass::index, common::GapClass::magnitude, common::MantPattern::Aaab, common::MantPattern::AllOnes, common::MantPattern::Alt5555, common::MantPattern::C000, common::MantPattern::MaxMinus1, common::MantPattern::SingleBit, common::MantPattern::TieAtP, common::MantPattern::TieMinus1, common::MantPattern::TiePlus1, common::MantPattern::Top, common::MantPattern::TopPlus1, common::MantPattern::Uniform, common::MantPattern::from_index, common::MantPattern::index, common::ModeOrNone::None, common::ModeOrNone::Valid, common::ModeOrNone::all_with_none, common::ModeOrNone::is_valid, common::ModeOrNone::name, common::PrecSource::Fast, common::PrecSource::Operand, common::PrecSource::Schedule, common::PrecSource::Uniform, common::RadiusClass::HalfUlp, common::RadiusClass::Huge, common::RadiusClass::Infinity, common::RadiusClass::OneUlp, common::RadiusClass::RandomMag, common::RadiusClass::Zero, common::RadiusClass::from_index, common::RadiusClass::index, common::RadiusClass::tag, common::RoundingMode::Down, common::RoundingMode::FromZero, common::RoundingMode::ToEven, common::RoundingMode::ToOdd, common::RoundingMode::ToZero, common::RoundingMode::Up, common::RoundingMode::all, common::RoundingMode::index, common::RoundingMode::name, common::RoundingMode::parse, common::Sign::Minus, common::Sign::Plus, common::Sign::symbol, common::generators::check_exp_mix, common::generators::check_mant_mix, common::generators::draw_ball, common::generators::draw_exponent, common::generators::draw_exponent_pair, common::generators::draw_exponent_with, common::generators::draw_gap, common::generators::draw_mantissa, common::generators::draw_mantissa_with, common::generators::draw_mode, common::generators::draw_precision, common::generators::forced_exponents, common::generators::forced_prop_exponents, common::generators::gap_pair, common::generators::tie_mantissa, common::generators::uniform_in
//! Self-tests of tests/common/generators.rs. Oracle: the plan (TEST_PLAN section 7, 9), hand-derived expectations, and pinned counts/digests
//! for the fixed seed `seed_for("common", CASE_CLASS_GATE)` at GENERATOR_VERSION 5 (a drift detector, produced by this generator; the band
//! counts were also recomputed independently by scripts/goldlib.py, see the comment at EXP_BAND_COUNTS).

use crate::common::expect_eq;
use crate::common::generators::{
    check_exp_mix, check_mant_mix, draw_ball, draw_exponent, draw_exponent_pair,
    draw_exponent_with, draw_gap, draw_mantissa, draw_mantissa_with, draw_mode, draw_precision,
    forced_exponents, forced_prop_exponents, gap_pair, tie_mantissa, uniform_in, BallDesc, ExpBand,
    ExpDomain, GapClass, MantPattern, ModeOrNone, PrecSource, RadiusClass, RoundingMode, Sign,
    EXP_MIX, EXP_MIX_SPEC_PERCENT, GUARD, MANT_MIX, MAX_PREC, PS, PS_FAST, RADIUS_MIX, TOP,
    ZF_EXPONENT_MAX, ZF_EXPONENT_MIN, ZF_MIN_NORMAL_EXPONENT,
};
use crate::common::rng::SplitMix64;
use crate::common::seeds::{seed_for, CASE_CLASS_GATE, GENERATOR_VERSION};

const N: usize = 10_000;

fn pinned_rng() -> Result<SplitMix64, String> {
    seed_for("common", CASE_CLASS_GATE)
        .map(SplitMix64::new)
        .ok_or_else(|| String::from("suite common is not in the seeds table"))
}

fn fold(h: u64, v: u64) -> u64 {
    crate::common::modular::mul_u64(h ^ v, 0x0000_0100_0000_01B3)
}

const FNV_START: u64 = 0xCBF2_9CE4_8422_2325;

fn bits(v: i64) -> u64 {
    u64::from_ne_bytes(v.to_ne_bytes())
}

fn near_power_of_two(v: i64) -> bool {
    let a = v.unsigned_abs();
    a.checked_sub(1).is_some_and(u64::is_power_of_two)
        || a.is_power_of_two()
        || a.checked_add(1).is_some_and(u64::is_power_of_two)
}

/// Band counts (Edge, Near0, Uniform, PowerOfTwo, Gap) and the FNV-1a digest of the values of `n` draws.
fn exponent_run(
    mix: &[(ExpBand, u64)],
    dom: &ExpDomain,
    n: usize,
) -> Result<([usize; 5], u64), String> {
    let mut rng = pinned_rng()?;
    let mut counts = [0usize; 5];
    let mut digest = FNV_START;
    for _ in 0..n {
        let d = draw_exponent_with(&mut rng, dom, 0, 53, mix)?;
        if let Some(c) = counts.get_mut(d.band.index()) {
            *c += 1;
        }
        digest = fold(digest, bits(d.value));
    }
    Ok((counts, digest))
}

// Pinned for GENERATOR_VERSION 5, unchanged since 4 (MAG domain; band counts unchanged since 1; the digest moved at 3 with the log-scale Uniform band), seed_for("common", CASE_CLASS_GATE), 10,000 draws, MAG domain, anchor 0, prec 53.
// The digests are pinned from this generator (drift detectors only). The counts were recomputed independently in Python (scripts/goldlib.py splitmix64 + the band rule `below(100)` against the cumulative
// 30/50/80/90/100 thresholds, with the fixed consumption per band: Edge 1, Near0 1, Uniform 1, PowerOfTwo 3, Gap 3 extra outputs).
const EXP_BAND_COUNTS: [usize; 5] = [2976, 2040, 3008, 994, 982];
const EXP_DIGEST_MAG: u64 = 3_632_968_743_847_705_314;

#[test]
fn schedules_constants_and_zf_range_are_as_documented() -> Result<(), String> {
    expect_eq("PS length", &PS.len(), &41)?;
    expect_eq("PS_FAST length", &PS_FAST.len(), &18)?;
    for (i, p) in PS.iter().take(16).enumerate() {
        expect_eq(
            "PS starts with 1..=16",
            p,
            &(u32::try_from(i).map_err(|_| String::from("index"))? + 1),
        )?;
    }
    if !PS.windows(2).all(|w| w.first() < w.get(1)) {
        return Err(String::from("PS is not strictly increasing"));
    }
    if !PS_FAST.windows(2).all(|w| w.first() < w.get(1)) {
        return Err(String::from("PS_FAST is not strictly increasing"));
    }
    if !PS_FAST.iter().all(|p| PS.contains(p)) {
        return Err(String::from("PS_FAST is not a subset of PS"));
    }
    expect_eq("last PS", &PS.last(), &Some(&MAX_PREC))?;
    expect_eq("MAX_PREC", &MAX_PREC, &4096)?;
    expect_eq("GUARD", &GUARD, &64)?;
    expect_eq("TOP", &TOP, &(1u64 << 63))?;
    // zenith-float 2.0 moves `Exponent` from i32 to i64 (Jeff, 2026-10-23 12:36 PM CT) with range +-(2^62 - 1) (Rusty)
    expect_eq("ZF_EXPONENT_MIN", &ZF_EXPONENT_MIN, &(-(1i64 << 62) + 1))?;
    expect_eq("ZF_EXPONENT_MAX", &ZF_EXPONENT_MAX, &((1i64 << 62) - 1))?;
    expect_eq(
        "ZF_MIN_NORMAL_EXPONENT",
        &ZF_MIN_NORMAL_EXPONENT,
        &(ZF_EXPONENT_MIN + 4096),
    )?;
    expect_eq(
        "EXP_MIX_SPEC_PERCENT total",
        &EXP_MIX_SPEC_PERCENT.iter().sum::<u64>(),
        &100,
    )?;
    expect_eq("RADIUS_MIX total", &RADIUS_MIX.iter().sum::<u64>(), &100)
}

#[test]
fn exp_mix_table_matches_the_plan_and_wrong_tables_are_rejected() -> Result<(), String> {
    check_exp_mix(&EXP_MIX)?;
    for i in 0..5 {
        let band = ExpBand::from_index(i).ok_or_else(|| String::from("from_index"))?;
        expect_eq("index round trip", &band.index(), &i)?;
    }
    if ExpBand::from_index(5).is_some() {
        return Err(String::from("from_index(5) is not None"));
    }
    // deliberately wrong tables must all be rejected
    let swapped_weights = vec![
        (ExpBand::Edge, 20),
        (ExpBand::Near0, 30),
        (ExpBand::Uniform, 30),
        (ExpBand::PowerOfTwo, 10),
        (ExpBand::Gap, 10),
    ];
    let swapped_bands = vec![
        (ExpBand::Near0, 30),
        (ExpBand::Edge, 20),
        (ExpBand::Uniform, 30),
        (ExpBand::PowerOfTwo, 10),
        (ExpBand::Gap, 10),
    ];
    let short_total = vec![
        (ExpBand::Edge, 30),
        (ExpBand::Near0, 20),
        (ExpBand::Uniform, 30),
        (ExpBand::PowerOfTwo, 10),
        (ExpBand::Gap, 0),
    ];
    let four_bands = vec![
        (ExpBand::Edge, 30),
        (ExpBand::Near0, 20),
        (ExpBand::Uniform, 40),
        (ExpBand::PowerOfTwo, 10),
    ];
    for (name, table) in [
        ("swapped weights", &swapped_weights),
        ("swapped bands", &swapped_bands),
        ("short total", &short_total),
        ("four bands", &four_bands),
    ] {
        if check_exp_mix(table).is_ok() {
            return Err(format!("check_exp_mix accepted the wrong table: {}", name));
        }
    }
    // draw_exponent_with honors the table it is given: a table of one band draws only that band
    let only_near0 = vec![(ExpBand::Near0, 100)];
    let mut rng = pinned_rng()?;
    for _ in 0..200 {
        let d = draw_exponent_with(&mut rng, &ExpDomain::MAG, 0, 53, &only_near0)?;
        if d.band != ExpBand::Near0 || d.value.abs() > 64 {
            return Err(format!("a Near0-only table gave {:?} {}", d.band, d.value));
        }
    }
    let empty: Vec<(ExpBand, u64)> = Vec::new();
    if draw_exponent_with(&mut rng, &ExpDomain::MAG, 0, 53, &empty).is_ok() {
        return Err(String::from("an empty table was accepted"));
    }
    // and a generator driven by a wrong table produces different counts than the pinned ones
    let (right, _) = exponent_run(&EXP_MIX, &ExpDomain::MAG, N)?;
    let (wrong, _) = exponent_run(&swapped_weights, &ExpDomain::MAG, N)?;
    if right == wrong {
        return Err(String::from(
            "a swapped-weights table gave the same band counts as the plan's table",
        ));
    }
    Ok(())
}

#[test]
fn exponent_band_counts_are_pinned_for_the_fixed_seed() -> Result<(), String> {
    expect_eq(
        "GENERATOR_VERSION (re-pin the counts below when it changes)",
        &GENERATOR_VERSION,
        &5,
    )?;
    let (counts, digest) = exponent_run(&EXP_MIX, &ExpDomain::MAG, N)?;
    expect_eq(
        "band counts Edge, Near0, Uniform, PowerOfTwo, Gap",
        &counts,
        &EXP_BAND_COUNTS,
    )?;
    expect_eq("digest of the values", &digest, &EXP_DIGEST_MAG)?;
    // statistical sanity, independent of the pin: each band within 2.5 percentage points of its weight
    for (i, (band, pct)) in EXP_MIX.iter().enumerate() {
        let got =
            u64::try_from(counts.get(i).map_or(0, |&c| c)).map_err(|_| String::from("count"))?;
        let want = u64::try_from(N).map_err(|_| String::from("N"))? * pct / 100;
        if got.abs_diff(want) * 40 > u64::try_from(N).map_err(|_| String::from("N"))? {
            return Err(format!(
                "band {:?}: {} draws, expected about {}",
                band, got, want
            ));
        }
    }
    // the band choice does not depend on the domain (same consumption per band), so the MID domain has the same counts
    let (mid_counts, _) = exponent_run(&EXP_MIX, &ExpDomain::MID, N)?;
    expect_eq("MID band counts", &mid_counts, &counts)
}

#[test]
fn every_exponent_draw_obeys_its_band() -> Result<(), String> {
    for dom in [&ExpDomain::MAG, &ExpDomain::MID] {
        let mut rng = pinned_rng()?;
        let edge = dom.edge_values();
        let mut seen_edge: Vec<i64> = Vec::new();
        let mut big_uniform = false;
        for _ in 0..N {
            let d = draw_exponent(&mut rng, dom, 0, 24)?;
            if !dom.in_raw_range(d.value) {
                return Err(format!(
                    "{:?} draw {} is outside the raw range",
                    d.band, d.value
                ));
            }
            match d.band {
                ExpBand::Edge => {
                    if !edge.contains(&d.value) {
                        return Err(format!("Edge draw {} is not in the edge table", d.value));
                    }
                    if !seen_edge.contains(&d.value) {
                        seen_edge.push(d.value);
                    }
                }
                ExpBand::Near0 => {
                    if d.value.abs() > 64 {
                        return Err(format!("Near0 draw {} is outside +-64", d.value));
                    }
                }
                ExpBand::Uniform => {
                    if d.value.unsigned_abs() > dom.max.unsigned_abs() / 2 {
                        big_uniform = true;
                    }
                }
                ExpBand::PowerOfTwo => {
                    if !near_power_of_two(d.value) {
                        return Err(format!("PowerOfTwo draw {} is not 2^k +-1", d.value));
                    }
                }
                ExpBand::Gap => {}
            }
        }
        if !big_uniform {
            return Err(String::from(
                "no Uniform draw above half the domain maximum",
            ));
        }
        let mut want = dom.edge_values();
        want.sort_unstable();
        seen_edge.sort_unstable();
        expect_eq("every edge value is drawn", &seen_edge, &want)?;
    }
    Ok(())
}

#[test]
fn edge_tables_hold_the_extremes_the_subnormal_boundary_and_the_raw_off_by_ones(
) -> Result<(), String> {
    let mag = ExpDomain::MAG.edge_values();
    expect_eq("MAG edge length", &mag.len(), &8)?;
    for v in [i64::MIN, i64::MIN + 3, i64::MAX, i64::MAX - 3] {
        if !mag.contains(&v) {
            return Err(format!("MAG edge table lacks {}", v));
        }
    }
    let mid = ExpDomain::MID.edge_values();
    // zf 2.0 range +-(2^62 - 1): the raw off-by-ones min - 1 and max + 1 are representable, so all 13 entries are present
    expect_eq("MID edge length", &mid.len(), &13)?;
    for v in [
        ZF_EXPONENT_MIN - 1,
        ZF_EXPONENT_MIN,
        ZF_EXPONENT_MIN + 3,
        ZF_EXPONENT_MAX + 1,
        ZF_EXPONENT_MAX,
        ZF_EXPONENT_MAX - 3,
        ZF_MIN_NORMAL_EXPONENT - 1,
        ZF_MIN_NORMAL_EXPONENT,
        ZF_MIN_NORMAL_EXPONENT + 1,
    ] {
        if !mid.contains(&v) {
            return Err(format!("MID edge table lacks {}", v));
        }
    }
    // MID raw range: one tick outside +-(2^62 - 1) and no further
    if !ExpDomain::MID.in_raw_range(ZF_EXPONENT_MIN - 1)
        || !ExpDomain::MID.in_raw_range(ZF_EXPONENT_MAX + 1)
    {
        return Err(String::from("MID raw range excludes a raw off-by-one"));
    }
    if ExpDomain::MID.in_raw_range(ZF_EXPONENT_MIN - 2)
        || ExpDomain::MID.in_raw_range(ZF_EXPONENT_MAX + 2)
        || ExpDomain::MID.in_raw_range(i64::MIN)
        || ExpDomain::MID.in_raw_range(i64::MAX)
    {
        return Err(String::from(
            "MID raw range reaches past one tick outside the zf range",
        ));
    }
    if !ExpDomain::MAG.in_raw_range(i64::MIN) || !ExpDomain::MAG.in_raw_range(i64::MAX) {
        return Err(String::from("MAG raw range excludes an i64 extreme"));
    }
    Ok(())
}

#[test]
fn forced_exponents_are_the_zf_bounds_then_cosmology_scale_with_both_signs() -> Result<(), String> {
    // hand-derived: +-(2^62 - 1) and one inside each, then 2^8 .. 2^20 magnitudes (Jeff's cosmology range) each with its negation
    let m = (1i64 << 62) - 1;
    let want: Vec<i64> = vec![
        -m,
        -m + 1,
        m - 1,
        m,
        256,
        -256,
        405,
        -405,
        1000,
        -1000,
        1024,
        -1024,
        4097,
        -4097,
        65_536,
        -65_536,
        1_048_575,
        -1_048_575,
        1_048_576,
        -1_048_576,
    ];
    expect_eq("forced exponents", &forced_exponents(), &want)?;
    for e in forced_exponents() {
        if !ExpDomain::MID.in_raw_range(e) || !(ZF_EXPONENT_MIN..=ZF_EXPONENT_MAX).contains(&e) {
            return Err(format!("forced exponent {} is outside the zf range", e));
        }
    }
    Ok(())
}

#[test]
fn draws_at_the_extremes_never_wrap() -> Result<(), String> {
    let anchors = [
        i64::MIN,
        i64::MIN + 1,
        -1,
        0,
        1,
        i64::MAX - 1,
        i64::MAX,
        ZF_EXPONENT_MIN,
        ZF_EXPONENT_MAX,
    ];
    for dom in [&ExpDomain::MAG, &ExpDomain::MID] {
        for anchor in anchors.iter() {
            for prec in [1u32, 53, 4096] {
                let mut rng = pinned_rng()?;
                for _ in 0..500 {
                    match draw_exponent(&mut rng, dom, *anchor, prec) {
                        Ok(d) => {
                            if !dom.in_raw_range(d.value) {
                                return Err(format!(
                                    "anchor {} prec {}: {:?} gave {}",
                                    anchor, prec, d.band, d.value
                                ));
                            }
                        }
                        // an anchor outside the domain may be refused (Gap band), never wrapped into range silently
                        Err(_) if !dom.in_raw_range(*anchor) => {}
                        Err(e) => return Err(format!("anchor {} prec {}: {}", anchor, prec, e)),
                    }
                    // overflow is an error, never a wrap
                    if let Ok(((a, b), g)) = draw_exponent_pair(&mut rng, *anchor, prec, dom) {
                        expect_eq("pair base", &a, anchor)?;
                        expect_eq("pair gap", &(i128::from(b) - i128::from(a)), &g.gap)?;
                    }
                }
            }
        }
    }
    Ok(())
}

#[test]
fn gap_pair_reports_overflow_as_an_error() -> Result<(), String> {
    expect_eq(
        "in range",
        &gap_pair(100, -72, &ExpDomain::MAG),
        &Ok((100, 28)),
    )?;
    expect_eq(
        "zero gap",
        &gap_pair(i64::MAX, 0, &ExpDomain::MAG),
        &Ok((i64::MAX, i64::MAX)),
    )?;
    for (base, gap) in [
        (i64::MAX, 1i128),
        (i64::MIN, -1),
        (i64::MAX, i128::from(u64::MAX)),
        (0, i128::MAX),
        (0, i128::MIN),
        (i64::MIN, -(1i128 << 64)),
    ] {
        if gap_pair(base, gap, &ExpDomain::MAG).is_ok() {
            return Err(format!(
                "gap_pair({}, {}) did not report overflow",
                base, gap
            ));
        }
    }
    // MID is +-(2^62 - 1): one tick past either bound is a raw value, two ticks past is an error, never a wrap
    expect_eq(
        "MID one tick past max",
        &gap_pair(ZF_EXPONENT_MAX, 1, &ExpDomain::MID),
        &Ok((ZF_EXPONENT_MAX, ZF_EXPONENT_MAX + 1)),
    )?;
    expect_eq(
        "MID one tick past min",
        &gap_pair(ZF_EXPONENT_MIN, -1, &ExpDomain::MID),
        &Ok((ZF_EXPONENT_MIN, ZF_EXPONENT_MIN - 1)),
    )?;
    if gap_pair(ZF_EXPONENT_MAX, 2, &ExpDomain::MID).is_ok()
        || gap_pair(ZF_EXPONENT_MIN, -2, &ExpDomain::MID).is_ok()
    {
        return Err(String::from(
            "MID gap two ticks past a zf bound was accepted",
        ));
    }
    Ok(())
}

#[test]
fn gap_draws_hit_every_listed_value() -> Result<(), String> {
    let p = 7u32;
    let g = i128::from(GUARD);
    let pi = i128::from(p);
    let mut want: Vec<i128> = vec![0];
    for m in [1, 2, pi - 1, pi, pi + 1, pi + g - 1, pi + g, pi + g + 1] {
        want.push(m);
        want.push(-m);
    }
    want.sort_unstable();
    want.dedup();
    let mut rng = pinned_rng()?;
    let mut seen: Vec<i128> = Vec::new();
    let mut class_seen = [false; 10];
    let mut random_large = false;
    for _ in 0..4000 {
        let d = draw_gap(&mut rng, p, &ExpDomain::MAG)?;
        if let Some(f) = class_seen.get_mut(d.class.index()) {
            *f = true;
        }
        if d.class == GapClass::Random {
            if d.gap.unsigned_abs() > (1u128 << 40) {
                random_large = true;
            }
        } else if !seen.contains(&d.gap) {
            seen.push(d.gap);
        }
    }
    seen.sort_unstable();
    expect_eq("fixed gap values seen", &seen, &want)?;
    expect_eq("every class seen", &class_seen, &[true; 10])?;
    if !random_large {
        return Err(String::from("Random gaps never exceeded 2^40"));
    }
    // class table: magnitudes at p = 7 and the Random class passes its argument through
    let expected: [(GapClass, i128); 10] = [
        (GapClass::Zero, 0),
        (GapClass::One, 1),
        (GapClass::Two, 2),
        (GapClass::PMinus1, 6),
        (GapClass::P, 7),
        (GapClass::PPlus1, 8),
        (GapClass::PGuardMinus1, 70),
        (GapClass::PGuard, 71),
        (GapClass::PGuardPlus1, 72),
        (GapClass::Random, 12345),
    ];
    for (i, (class, mag)) in expected.iter().enumerate() {
        expect_eq("class index", &class.index(), &i)?;
        expect_eq(
            "from_index",
            &GapClass::from_index(i).as_ref(),
            &Some(class),
        )?;
        expect_eq("magnitude", &class.magnitude(p, 12345), &Ok(*mag))?;
    }
    if GapClass::from_index(10).is_some() {
        return Err(String::from("GapClass::from_index(10) is not None"));
    }
    // precision 1 gives p - 1 = 0; precision 0 is an error
    expect_eq(
        "PMinus1 at p = 1",
        &GapClass::PMinus1.magnitude(1, 0),
        &Ok(0),
    )?;
    if draw_gap(&mut rng, 0, &ExpDomain::MAG).is_ok() {
        return Err(String::from("draw_gap accepted precision 0"));
    }
    if GapClass::PMinus1.magnitude(0, 0).is_ok() {
        return Err(String::from(
            "p - 1 at p = 0 is negative and must be an error",
        ));
    }
    Ok(())
}

#[test]
fn uniform_in_stays_in_range_and_rejects_an_empty_range() -> Result<(), String> {
    let mut rng = pinned_rng()?;
    expect_eq("single value", &uniform_in(&mut rng, 5, 5), &Ok(5))?;
    for _ in 0..1000 {
        let v = uniform_in(&mut rng, -3, 4)?;
        if !(-3..=4).contains(&v) {
            return Err(format!("{} outside -3..=4", v));
        }
        uniform_in(&mut rng, i64::MIN, i64::MAX)?;
    }
    if uniform_in(&mut rng, 1, 0).is_ok() {
        return Err(String::from("empty range accepted"));
    }
    Ok(())
}

#[test]
fn mantissa_patterns_all_appear_have_the_top_bit_and_match_their_definition() -> Result<(), String>
{
    let p = 24u32;
    let tie = TOP | (1u64 << 39);
    expect_eq("tie_mantissa(24)", &tie_mantissa(p), &Ok(tie))?;
    expect_eq("tie_mantissa(63)", &tie_mantissa(63), &Ok(TOP | 1))?;
    expect_eq("tie_mantissa(64) degenerates", &tie_mantissa(64), &Ok(TOP))?;
    expect_eq(
        "tie_mantissa(4096) degenerates",
        &tie_mantissa(4096),
        &Ok(TOP),
    )?;
    if tie_mantissa(0).is_ok() {
        return Err(String::from("tie_mantissa(0) accepted"));
    }
    let mut rng = pinned_rng()?;
    let mut counts = [0usize; 12];
    let mut digest = FNV_START;
    for _ in 0..N {
        let d = draw_mantissa(&mut rng, p)?;
        if d.value & TOP == 0 {
            return Err(format!(
                "{:?} mantissa {:#x} has no top bit",
                d.pattern, d.value
            ));
        }
        if let Some(c) = counts.get_mut(d.pattern.index()) {
            *c += 1;
        }
        digest = fold(digest, d.value);
        let ok = match d.pattern {
            MantPattern::Top => d.value == TOP,
            MantPattern::TopPlus1 => d.value == TOP + 1,
            MantPattern::AllOnes => d.value == u64::MAX,
            MantPattern::MaxMinus1 => d.value == u64::MAX - 1,
            MantPattern::C000 => d.value == 0xC000_0000_0000_0000,
            MantPattern::Aaab => d.value == 0xAAAA_AAAA_AAAA_AAAB,
            MantPattern::Alt5555 => d.value == 0xD555_5555_5555_5555,
            MantPattern::SingleBit => (d.value & !TOP).count_ones() == 1,
            MantPattern::TieAtP => d.value == tie,
            MantPattern::TieMinus1 => d.value == tie - 1,
            MantPattern::TiePlus1 => d.value == tie + 1,
            MantPattern::Uniform => true,
        };
        if !ok {
            return Err(format!("{:?} gave {:#x}", d.pattern, d.value));
        }
    }
    if counts.contains(&0) {
        return Err(format!("a mantissa pattern never appeared: {:?}", counts));
    }
    expect_eq("pattern counts", &counts, &MANT_PATTERN_COUNTS)?;
    expect_eq("mantissa digest", &digest, &MANT_DIGEST)?;
    for i in 0..12 {
        let pat = MantPattern::from_index(i).ok_or_else(|| String::from("from_index"))?;
        expect_eq("MantPattern index", &pat.index(), &i)?;
    }
    if MantPattern::from_index(12).is_some() {
        return Err(String::from("MantPattern::from_index(12) is not None"));
    }
    Ok(())
}

const MANT_PATTERN_COUNTS: [usize; 12] =
    [829, 577, 749, 623, 619, 613, 651, 766, 815, 596, 580, 2582];
const MANT_DIGEST: u64 = 10_607_778_329_162_618_863;

#[test]
fn mantissa_mix_table_is_checked_and_a_wrong_table_is_rejected() -> Result<(), String> {
    check_mant_mix(&MANT_MIX)?;
    let mut wrong: Vec<(MantPattern, u64)> = Vec::new();
    for (i, (_, pct)) in MANT_MIX.iter().enumerate() {
        let pat = MantPattern::from_index(i).ok_or_else(|| String::from("from_index"))?;
        wrong.push((pat, if i == 0 { pct + 1 } else { *pct }));
    }
    if check_mant_mix(&wrong).is_ok() {
        return Err(String::from("a mix totalling 101 was accepted"));
    }
    let mut zero_weight: Vec<(MantPattern, u64)> = Vec::new();
    for (i, (_, pct)) in MANT_MIX.iter().enumerate() {
        let pat = MantPattern::from_index(i).ok_or_else(|| String::from("from_index"))?;
        zero_weight.push((
            pat,
            if i == 1 {
                0
            } else if i == 0 {
                pct + 6
            } else {
                *pct
            },
        ));
    }
    if check_mant_mix(&zero_weight).is_ok() {
        return Err(String::from("a zero-weight pattern was accepted"));
    }
    let wrong_order = vec![(MantPattern::TopPlus1, 50), (MantPattern::Top, 50)];
    if check_mant_mix(&wrong_order).is_ok() {
        return Err(String::from("a short, reordered table was accepted"));
    }
    // a wrong mix changes what is drawn
    let only_top = vec![(MantPattern::Top, 100)];
    let mut rng = pinned_rng()?;
    expect_eq(
        "only Top",
        &draw_mantissa_with(&mut rng, 10, &only_top)?.value,
        &TOP,
    )?;
    let short = vec![(MantPattern::Top, 50)];
    let mut found_error = false;
    for _ in 0..100 {
        if draw_mantissa_with(&mut rng, 10, &short).is_err() {
            found_error = true;
        }
    }
    if !found_error {
        return Err(String::from(
            "a table that does not cover 100 never produced an error",
        ));
    }
    if draw_mantissa(&mut rng, 0).is_ok() {
        return Err(String::from("draw_mantissa accepted precision 0"));
    }
    Ok(())
}

#[test]
fn precision_draws_come_from_their_source() -> Result<(), String> {
    let mut rng = pinned_rng()?;
    let mut seen_ps = [false; 41];
    let mut seen_fast = [false; 18];
    let mut seen_op: Vec<u32> = Vec::new();
    let mut min_u = u32::MAX;
    let mut max_u = 0u32;
    let mut above_max_prec = 0usize;
    let mut lengths = [false; 32];
    for _ in 0..N {
        let s = draw_precision(&mut rng, &PrecSource::Schedule)?;
        let i = PS
            .iter()
            .position(|p| *p == s)
            .ok_or_else(|| format!("{} is not in PS", s))?;
        if let Some(f) = seen_ps.get_mut(i) {
            *f = true;
        }
        let f = draw_precision(&mut rng, &PrecSource::Fast)?;
        let j = PS_FAST
            .iter()
            .position(|p| *p == f)
            .ok_or_else(|| format!("{} is not in PS_FAST", f))?;
        if let Some(x) = seen_fast.get_mut(j) {
            *x = true;
        }
        let u = draw_precision(&mut rng, &PrecSource::Uniform)?;
        if u == 0 {
            return Err(String::from("uniform precision 0"));
        }
        min_u = min_u.min(u);
        max_u = max_u.max(u);
        if u > MAX_PREC {
            above_max_prec += 1;
        }
        if let Some(slot) = usize::try_from(31 - u.leading_zeros())
            .ok()
            .and_then(|b| lengths.get_mut(b))
        {
            *slot = true;
        }
        let o = draw_precision(&mut rng, &PrecSource::Operand(53))?;
        if !seen_op.contains(&o) {
            seen_op.push(o);
        }
    }
    expect_eq("all 41 schedule values drawn", &seen_ps, &[true; 41])?;
    expect_eq("all 18 fast values drawn", &seen_fast, &[true; 18])?;
    // headroom (Jeff, 2026-10-23 12:33 PM CT): uniform precisions span all of u32 (CONTRACT 8.1), every bit length, both ends exactly
    if min_u != 1 || max_u != u32::MAX {
        return Err(format!(
            "uniform precisions only span {}..={}",
            min_u, max_u
        ));
    }
    if above_max_prec < N / 4 {
        return Err(format!(
            "only {} of {} uniform precisions above {}",
            above_max_prec, N, MAX_PREC
        ));
    }
    expect_eq(
        "uniform precision bit lengths 1..=32",
        &lengths,
        &[true; 32],
    )?;
    seen_op.sort_unstable();
    expect_eq("operand +-1 around 53", &seen_op, &vec![52, 53, 54])?;
    for (op, want) in [
        (1u32, vec![1, 2]),
        (MAX_PREC, vec![MAX_PREC - 1, MAX_PREC, MAX_PREC + 1]),
    ] {
        let mut got: Vec<u32> = Vec::new();
        for _ in 0..200 {
            let v = draw_precision(&mut rng, &PrecSource::Operand(op))?;
            if !got.contains(&v) {
                got.push(v);
            }
        }
        got.sort_unstable();
        expect_eq("operand at the limit", &got, &want)?;
    }
    if draw_precision(&mut rng, &PrecSource::Operand(0)).is_ok() {
        return Err(String::from("operand precision 0 was accepted"));
    }
    Ok(())
}

// Headroom (Jeff, 2026-10-23 12:33 PM CT): an operand precision at the top of u32 keeps its -1 neighbour and itself; +1 would wrap, so it is dropped.
#[test]
fn operand_precision_at_u32_max_keeps_its_lower_neighbour() -> Result<(), String> {
    let mut rng = pinned_rng()?;
    let mut got: Vec<u32> = Vec::new();
    for _ in 0..200 {
        let v = draw_precision(&mut rng, &PrecSource::Operand(u32::MAX))?;
        if !got.contains(&v) {
            got.push(v);
        }
    }
    got.sort_unstable();
    expect_eq("operand at u32::MAX", &got, &vec![u32::MAX - 1, u32::MAX])
}

#[test]
fn ball_descriptors_cover_signs_zeros_radius_classes_and_independent_exponents(
) -> Result<(), String> {
    let mut rng = pinned_rng()?;
    let mut sign_pairs = [false; 4];
    let mut plus_zero = false;
    let mut minus_zero = false;
    let mut classes = [false; 6];
    let mut rad_above = false;
    let mut rad_below = false;
    let mut digest = FNV_START;
    for _ in 0..4000 {
        let a = draw_ball(&mut rng, 53, &ExpDomain::MID)?;
        let b = draw_ball(&mut rng, 53, &ExpDomain::MID)?;
        let k = usize::from(a.sign == Sign::Minus) * 2 + usize::from(b.sign == Sign::Minus);
        if let Some(f) = sign_pairs.get_mut(k) {
            *f = true;
        }
        for ball in [&a, &b] {
            if ball.mid_is_zero {
                match ball.sign {
                    Sign::Plus => plus_zero = true,
                    Sign::Minus => minus_zero = true,
                }
            }
            if let Some(f) = classes.get_mut(ball.radius.index()) {
                *f = true;
            }
            if ball.mid_mant & TOP == 0 || ball.rad_mant & TOP == 0 {
                return Err(format!("mantissa without the top bit in {}", ball.encode()));
            }
            if ball.rad_exp > ball.mid_exp {
                rad_above = true;
            }
            if ball.rad_exp < ball.mid_exp {
                rad_below = true;
            }
            let text = ball.encode();
            if text.is_empty() || text.chars().any(char::is_whitespace) {
                return Err(format!(
                    "encoding {:?} is empty or contains whitespace",
                    text
                ));
            }
            for byte in text.bytes() {
                digest = fold(digest, u64::from(byte));
            }
        }
    }
    expect_eq("all four sign pairs", &sign_pairs, &[true; 4])?;
    if !plus_zero || !minus_zero {
        return Err(String::from("both zero mids (+0 and -0) must appear"));
    }
    expect_eq("all radius classes", &classes, &[true; 6])?;
    if !rad_above || !rad_below {
        return Err(String::from(
            "radius exponents are not independent of the mid exponent",
        ));
    }
    expect_eq("ball digest", &digest, &BALL_DIGEST)?;
    for i in 0..6 {
        let c = RadiusClass::from_index(i).ok_or_else(|| String::from("from_index"))?;
        expect_eq("RadiusClass index", &c.index(), &i)?;
    }
    if RadiusClass::from_index(6).is_some() {
        return Err(String::from("RadiusClass::from_index(6) is not None"));
    }
    // encodings of hand-built descriptors
    let zero = BallDesc {
        sign: Sign::Minus,
        mid_is_zero: true,
        mid_exp: 7,
        mid_mant: TOP,
        radius: RadiusClass::Zero,
        rad_exp: 0,
        rad_mant: TOP,
    };
    expect_eq("signed zero", &zero.encode(), &String::from("-0~r=zero"))?;
    let finite = BallDesc {
        sign: Sign::Plus,
        mid_is_zero: false,
        mid_exp: -5,
        mid_mant: 0xC000_0000_0000_0000,
        radius: RadiusClass::RandomMag,
        rad_exp: 9,
        rad_mant: TOP,
    };
    expect_eq(
        "finite",
        &finite.encode(),
        &String::from("+c000000000000000e-5~r=mag(8000000000000000,e9)"),
    )?;
    for (class, tag) in [
        (RadiusClass::OneUlp, "ulp"),
        (RadiusClass::HalfUlp, "half"),
        (RadiusClass::Huge, "huge"),
        (RadiusClass::Infinity, "inf"),
    ] {
        expect_eq("tag", &class.tag(), &tag)?;
    }
    expect_eq(
        "symbols",
        &(Sign::Plus.symbol(), Sign::Minus.symbol()),
        &('+', '-'),
    )
}

// Re-pinned at GENERATOR_VERSION 5: the MID domain moved to +-(2^62 - 1) (Rusty), so mid exponents drawn by ball descriptors changed.
const BALL_DIGEST: u64 = 6_896_748_639_994_193_576;

#[test]
fn rounding_modes_are_six_named_distinct_and_none_is_separate() -> Result<(), String> {
    let all: Vec<RoundingMode> = RoundingMode::all().collect();
    expect_eq("six modes", &all.len(), &6)?;
    for (i, m) in all.iter().enumerate() {
        expect_eq("index in order", &m.index(), &i)?;
        expect_eq(
            "parse(name) round trip",
            &RoundingMode::parse(m.name()).as_ref(),
            &Some(m),
        )?;
    }
    let mut names: Vec<&str> = all.iter().map(|m| m.name()).collect();
    names.sort_unstable();
    names.dedup();
    expect_eq("distinct names", &names.len(), &6)?;
    expect_eq("ToEven", &RoundingMode::ToEven.name(), &"ToEven")?;
    expect_eq("ToOdd", &RoundingMode::ToOdd.name(), &"ToOdd")?;
    expect_eq(
        "Up Down",
        &(RoundingMode::Up.name(), RoundingMode::Down.name()),
        &("Up", "Down"),
    )?;
    expect_eq(
        "ToZero FromZero",
        &(RoundingMode::ToZero.name(), RoundingMode::FromZero.name()),
        &("ToZero", "FromZero"),
    )?;
    // None is not a RoundingMode
    if RoundingMode::parse("None").is_some()
        || RoundingMode::parse("").is_some()
        || RoundingMode::parse("toeven").is_some()
    {
        return Err(String::from(
            "RoundingMode::parse accepted a name that is not one of the six",
        ));
    }
    let with_none = ModeOrNone::all_with_none();
    expect_eq("seven entries", &with_none.len(), &7)?;
    let valid = with_none.iter().filter(|m| m.is_valid()).count();
    expect_eq("six valid", &valid, &6)?;
    expect_eq("last is None", &with_none.last(), &Some(&ModeOrNone::None))?;
    expect_eq("None name", &ModeOrNone::None.name(), &"None")?;
    expect_eq(
        "Valid name",
        &ModeOrNone::Valid(RoundingMode::Up).name(),
        &"Up",
    )?;
    if ModeOrNone::None.is_valid() || !ModeOrNone::Valid(RoundingMode::Down).is_valid() {
        return Err(String::from("is_valid is wrong"));
    }
    let mut rng = pinned_rng()?;
    let mut seen = [false; 6];
    for _ in 0..600 {
        let m = draw_mode(&mut rng)?;
        if let Some(f) = seen.get_mut(m.index()) {
            *f = true;
        }
    }
    expect_eq("every mode drawn", &seen, &[true; 6])
}

/// Cosmology-scale exponents (Jeff, 2026-10-23 12:20 PM CT): 1e±122 is about 2^±405, so a fixed-seed run must reach magnitudes
/// `2^8 <= |e| < 2^20` with exponents that are not near a power of two (the PowerOfTwo band only gives `±2^k ±1`), on both domains,
/// and every bit length up to the domain's must still be reachable (headroom well past cosmology).
#[test]
fn uniform_band_reaches_cosmology_scale_exponents_and_every_bit_length() -> Result<(), String> {
    for dom in [&ExpDomain::MAG, &ExpDomain::MID] {
        let mut rng = pinned_rng()?;
        let mut cosmology = 0u32;
        let mut lengths = [false; 64];
        for _ in 0..N {
            let d = draw_exponent(&mut rng, dom, 0, 53)?;
            if d.band != ExpBand::Uniform {
                continue;
            }
            let m = d.value.unsigned_abs();
            let near_pow2 = m.is_power_of_two()
                || m.checked_sub(1).is_some_and(|v| v.is_power_of_two())
                || m.checked_add(1).is_some_and(|v| v.is_power_of_two());
            if (1u64 << 8..1u64 << 20).contains(&m) && !near_pow2 {
                cosmology += 1;
            }
            if let Some(slot) = usize::try_from(64 - m.leading_zeros())
                .ok()
                .and_then(|b| b.checked_sub(1))
                .and_then(|b| lengths.get_mut(b))
            {
                *slot = true;
            }
        }
        if cosmology < 100 {
            return Err(format!(
                "only {} cosmology-scale Uniform draws in {}",
                cosmology, N
            ));
        }
        // MAG reaches bit length 63; MID stops at 62 (|e| <= 2^62 - 1)
        let top = if dom.max == ZF_EXPONENT_MAX { 62 } else { 63 };
        let missing = lengths.iter().take(top).filter(|hit| !**hit).count();
        expect_eq("bit lengths 1..=top missing", &missing, &0)?;
        let past = lengths.iter().skip(top).filter(|hit| **hit).count();
        expect_eq("bit lengths past the domain", &past, &0)?;
    }
    Ok(())
}

#[test]
fn forced_prop_exponents_are_the_zf_bounds_zero_then_the_1024_scale_with_both_signs(
) -> Result<(), String> {
    let max = (1i64 << 62) - 1;
    let want: Vec<i64> = vec![
        -max,
        -max + 1,
        max - 1,
        max,
        0,
        256,
        -256,
        405,
        -405,
        1000,
        -1000,
        1024,
        -1024,
    ];
    expect_eq("forced prop exponents", &forced_prop_exponents(), &want)
}
