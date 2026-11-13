//! COVERS: common::FuzzConfig::case_range, common::FuzzConfig::seed_of, common::ReplayLine::format_line, common::ReplayLine::parse_line, common::replay::case_seed, common::replay::config_from_env, common::replay::parse_config, common::replay::parse_hex_u64
//! Self-tests of tests/common/replay.rs. Oracle: hand-written lines and variable sets; `case_seed` values recomputed with an independent Python
//! integer-mask splitmix64 (2026-10-15). The process environment is never mutated: `parse_config` takes an injected lookup.

use crate::common::expect_eq;
use crate::common::generators::{draw_ball, ExpDomain};
use crate::common::replay::{
    case_seed, config_from_env, parse_config, parse_hex_u64, FuzzConfig, ReplayLine, ReplayTarget,
};
use crate::common::rng::SplitMix64;

fn sample() -> ReplayLine {
    ReplayLine {
        suite: String::from("mag_add"),
        seed: 0x0123_4567_89AB_CDEF,
        case: 42,
        op: String::from("add"),
        prec: 53,
        rm: String::from("ToEven"),
        a: String::from("+8000000000000000e-5~r=ulp"),
        b: String::from("-0~r=zero"),
    }
}

/// A lookup over fixed pairs; a name mapped to `None` simulates a variable that is not valid unicode.
type Pairs = Vec<(&'static str, Option<&'static str>)>;

fn env(pairs: Pairs) -> impl Fn(&str) -> Result<Option<String>, String> {
    move |name: &str| match pairs.iter().find(|(k, _)| *k == name) {
        Some((_, Some(v))) => Ok(Some(String::from(*v))),
        Some((_, None)) => Err(format!("{} is not valid unicode", name)),
        None => Ok(None),
    }
}

#[test]
fn replay_line_has_the_documented_shape_and_round_trips() -> Result<(), String> {
    let line = sample();
    let text = line.format_line()?;
    expect_eq(
        "line text",
        &text,
        &String::from("suite=mag_add seed=0123456789abcdef case=42 op=add prec=53 rm=ToEven a=+8000000000000000e-5~r=ulp b=-0~r=zero"),
    )?;
    expect_eq("round trip", &ReplayLine::parse_line(&text)?, &line)?;
    // generated ball encodings survive the round trip too
    let mut rng = SplitMix64::new(77);
    for case in 0..200u64 {
        let a = draw_ball(&mut rng, 53, &ExpDomain::MID)?;
        let b = draw_ball(&mut rng, 53, &ExpDomain::MID)?;
        let l = ReplayLine {
            suite: String::from("s"),
            seed: rng.next_u64(),
            case,
            op: String::from("fma"),
            prec: 4096,
            rm: String::from("FromZero"),
            a: a.encode(),
            b: b.encode(),
        };
        expect_eq(
            "generated round trip",
            &ReplayLine::parse_line(&l.format_line()?)?,
            &l,
        )?;
    }
    Ok(())
}

#[test]
fn malformed_replay_lines_are_errors_not_panics() -> Result<(), String> {
    let good = sample().format_line()?;
    let bad: [(&str, String); 13] = [
        ("empty", String::new()),
        ("missing key", good.replace(" b=-0~r=zero", "")),
        ("extra token", format!("{} extra=1", good)),
        (
            "swapped keys",
            good.replace(
                "seed=0123456789abcdef case=42",
                "case=42 seed=0123456789abcdef",
            ),
        ),
        ("no equals", good.replace("op=add", "opadd")),
        ("empty value", good.replace("rm=ToEven", "rm=")),
        (
            "non-hex seed",
            good.replace("0123456789abcdef", "0123456789abcdeg"),
        ),
        ("non-numeric case", good.replace("case=42", "case=x")),
        ("negative prec", good.replace("prec=53", "prec=-1")),
        (
            "prec too big",
            good.replace("prec=53", "prec=99999999999999999999"),
        ),
        (
            "prec above MAX_PRECISION",
            good.replace("prec=53", "prec=1152921504606846976"),
        ),
        (
            "seed overflow",
            good.replace("0123456789abcdef", "10123456789abcdef"),
        ),
        ("wrong key name", good.replace("suite=", "suit=")),
    ];
    for (name, text) in bad.iter() {
        if ReplayLine::parse_line(text).is_ok() {
            return Err(format!(
                "parse_line accepted a malformed line ({}): {}",
                name, text
            ));
        }
    }
    // fields that would not parse back are refused when formatting
    for field in 0..5 {
        let mut l = sample();
        match field {
            0 => l.suite = String::new(),
            1 => l.op = String::from("a b"),
            2 => l.rm = String::new(),
            3 => l.a = String::from("x y"),
            _ => l.b = String::from("tab\there"),
        }
        if l.format_line().is_ok() {
            return Err(format!("format_line accepted a bad field (case {})", field));
        }
    }
    Ok(())
}

#[test]
fn hex_numbers_parse_with_or_without_the_prefix_and_reject_junk() -> Result<(), String> {
    expect_eq("plain", &parse_hex_u64("x", "ff"), &Ok(255))?;
    expect_eq("prefixed", &parse_hex_u64("x", "0xFF"), &Ok(255))?;
    expect_eq(
        "max",
        &parse_hex_u64("x", "ffffffffffffffff"),
        &Ok(u64::MAX),
    )?;
    for text in ["", "0x", "g", "-1", "1ffffffffffffffff", " 1", "1 "] {
        if parse_hex_u64("x", text).is_ok() {
            return Err(format!("parse_hex_u64 accepted {:?}", text));
        }
    }
    Ok(())
}

#[test]
fn case_seeds_match_the_builders_rule_and_independent_values() -> Result<(), String> {
    expect_eq("case 0 is the base", &case_seed(0xABCD, 0), &0xABCD)?;
    expect_eq("0,1", &case_seed(0, 1), &0x6e78_9e6a_a1b9_65f4)?;
    expect_eq("0x1234,7", &case_seed(0x1234, 7), &0xf6b1_0aa2_f186_18d0)?;
    expect_eq("max,3", &case_seed(u64::MAX, 3), &0x738b_10af_1713_67ff)?;
    expect_eq(
        "0xdeadbeef,1000",
        &case_seed(0xDEAD_BEEF, 1000),
        &0xeb56_b74c_68d4_f7cb,
    )?;
    // distinct cases get distinct seeds
    let mut seen: Vec<u64> = (0..1000).map(|c| case_seed(5, c)).collect();
    seen.sort_unstable();
    seen.dedup();
    expect_eq("1000 distinct seeds", &seen.len(), &1000)
}

#[test]
fn config_defaults_and_the_preferred_names() -> Result<(), String> {
    let none = env(vec![]);
    let c = parse_config(&none, 0x55, 300)?;
    expect_eq(
        "defaults",
        &(c.seed, c.cases, c.replay.is_some()),
        &(0x55, 300, false),
    )?;
    expect_eq("range", &c.case_range(), &Ok(0..300))?;
    expect_eq("seed_of(0)", &c.seed_of(0), &0x55)?;
    expect_eq("seed_of(5)", &c.seed_of(5), &case_seed(0x55, 5))?;
    let fuzz = env(vec![
        ("ZB_FUZZ_SEED", Some("0xBEEF")),
        ("ZB_FUZZ_CASES", Some("12")),
    ]);
    let c = parse_config(&fuzz, 1, 2)?;
    expect_eq("ZB_FUZZ_*", &(c.seed, c.cases), &(0xBEEF, 12))?;
    let prop = env(vec![
        ("ZB_PROP_SEED", Some("0x10")),
        ("ZB_PROP_CASES", Some("7")),
    ]);
    let c = parse_config(&prop, 1, 2)?;
    expect_eq("ZB_PROP_* hex", &(c.seed, c.cases), &(0x10, 7))?;
    let prop_dec = env(vec![("ZB_PROP_SEED", Some("16"))]);
    expect_eq(
        "ZB_PROP_SEED decimal",
        &parse_config(&prop_dec, 1, 2)?.seed,
        &16,
    )?;
    let both = env(vec![
        ("ZB_FUZZ_SEED", Some("aa")),
        ("ZB_PROP_SEED", Some("0xbb")),
        ("ZB_FUZZ_CASES", Some("3")),
        ("ZB_PROP_CASES", Some("9")),
    ]);
    let c = parse_config(&both, 1, 2)?;
    expect_eq("ZB_FUZZ_* wins", &(c.seed, c.cases), &(0xAA, 3))
}

#[test]
fn replay_selects_exactly_one_case_from_its_own_seed() -> Result<(), String> {
    let replay = env(vec![
        ("ZB_REPLAY", Some("0xdeadbeef:1000")),
        ("ZB_FUZZ_SEED", Some("1")),
        ("ZB_FUZZ_CASES", Some("50")),
    ]);
    let c = parse_config(&replay, 9, 9)?;
    expect_eq(
        "replay target",
        &c.replay,
        &Some(ReplayTarget {
            seed: 0xDEAD_BEEF,
            case: 1000,
        }),
    )?;
    expect_eq("only that case", &c.case_range(), &Ok(1000..1001))?;
    expect_eq("seed of the case", &c.seed_of(1000), &0xeb56_b74c_68d4_f7cb)?;
    // the failure line carries the seed and case that ZB_REPLAY takes back
    let line = ReplayLine {
        seed: 0xDEAD_BEEF,
        case: 1000,
        ..sample()
    };
    let text = line.format_line()?;
    let parsed = ReplayLine::parse_line(&text)?;
    let spec = format!("{:x}:{}", parsed.seed, parsed.case);
    let again = |name: &str| -> Result<Option<String>, String> {
        if name == "ZB_REPLAY" {
            Ok(Some(String::from(spec.as_str())))
        } else {
            Ok(None)
        }
    };
    expect_eq(
        "replay from a line",
        &parse_config(&again, 0, 1)?.replay,
        &c.replay,
    )?;
    let last = FuzzConfig {
        seed: 0,
        cases: 1,
        replay: Some(ReplayTarget {
            seed: 0,
            case: u64::MAX,
        }),
    };
    if last.case_range().is_ok() {
        return Err(String::from(
            "a replay of case u64::MAX must not wrap its range end",
        ));
    }
    Ok(())
}

#[test]
fn malformed_environment_values_are_errors_not_panics() -> Result<(), String> {
    let cases: [(&str, Pairs); 15] = [
        ("seed not hex", vec![("ZB_FUZZ_SEED", Some("xyz"))]),
        ("seed empty", vec![("ZB_FUZZ_SEED", Some(""))]),
        (
            "seed too long",
            vec![("ZB_FUZZ_SEED", Some("10000000000000000"))],
        ),
        ("prop seed junk", vec![("ZB_PROP_SEED", Some("0xzz"))]),
        ("prop seed negative", vec![("ZB_PROP_SEED", Some("-5"))]),
        ("cases zero", vec![("ZB_FUZZ_CASES", Some("0"))]),
        ("cases junk", vec![("ZB_FUZZ_CASES", Some("many"))]),
        ("cases negative", vec![("ZB_FUZZ_CASES", Some("-1"))]),
        ("prop cases zero", vec![("ZB_PROP_CASES", Some("0"))]),
        ("replay without colon", vec![("ZB_REPLAY", Some("abc"))]),
        ("replay bad seed", vec![("ZB_REPLAY", Some("zz:1"))]),
        ("replay bad case", vec![("ZB_REPLAY", Some("ab:x"))]),
        ("replay empty case", vec![("ZB_REPLAY", Some("ab:"))]),
        ("not unicode", vec![("ZB_FUZZ_SEED", None)]),
        ("replay not unicode", vec![("ZB_REPLAY", None)]),
    ];
    for (name, pairs) in cases {
        let lookup = env(pairs);
        if parse_config(&lookup, 1, 1).is_ok() {
            return Err(format!("parse_config accepted: {}", name));
        }
    }
    Ok(())
}

#[test]
fn config_from_env_reads_the_process_environment() -> Result<(), String> {
    const NAMES: [&str; 5] = [
        "ZB_FUZZ_SEED",
        "ZB_FUZZ_CASES",
        "ZB_REPLAY",
        "ZB_PROP_SEED",
        "ZB_PROP_CASES",
    ];
    if NAMES.iter().any(|n| std::env::var_os(n).is_some()) {
        // the caller set the knobs: the result must be a valid config or a clean Err, and replaying must honor ZB_REPLAY
        return match config_from_env(1, 1) {
            Ok(c) => {
                if std::env::var_os("ZB_REPLAY").is_some() && c.replay.is_none() {
                    return Err(String::from(
                        "ZB_REPLAY is set but the config has no replay target",
                    ));
                }
                Ok(())
            }
            Err(e) => {
                if e.is_empty() {
                    Err(String::from("empty error"))
                } else {
                    Ok(())
                }
            }
        };
    }
    let c = config_from_env(0x77, 5)?;
    expect_eq(
        "unset environment gives the defaults",
        &(c.seed, c.cases, c.replay.is_some()),
        &(0x77, 5, false),
    )
}
