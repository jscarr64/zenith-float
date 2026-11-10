//! COVERS: common::seeds::stream_seed, common::seeds::seed_for
//! Also pins (not function coverage): BASE, GENERATOR_VERSION, ZF_PIN, SUITES, CASE_CLASS_GATE, CASE_CLASS_DEEP.
//! Oracle: stream_seed(suite, class) = splitmix64(BASE + suite * 2^32 + class), pinned values recomputed with an independent Python integer-mask implementation (2026-10-15).

use crate::common::expect_eq;
use crate::common::seeds::{
    seed_for, stream_seed, BASE, CASE_CLASS_DEEP, CASE_CLASS_GATE, GENERATOR_VERSION, SUITES,
    ZF_PIN,
};

#[test]
fn stream_seed_matches_independent_python_values() -> Result<(), String> {
    let cases: [(u32, u32, u64); 7] = [
        (0, 0, 0xA6A3_ED83_6C68_6CD0),
        (0, 1, 0xF226_8EFD_823F_4F33),
        (1, 0, 0x7596_AB86_7D07_4637),
        (1, 1, 0x8BC7_7992_C364_0E5B),
        (6, 1, 0xD243_F5FC_52AD_2C5A),
        (3, u32::MAX, 0x66AE_6CF5_7B32_420B),
        // BASE + (2^64 - 1) wraps: pinned so a checked/panicking variant would show up.
        (u32::MAX, u32::MAX, 0x86CF_DEB4_D158_6C75),
    ];
    for (suite, class, want) in cases.iter() {
        expect_eq(
            &format!("stream_seed({}, {})", suite, class),
            &stream_seed(*suite, *class),
            want,
        )?;
    }
    Ok(())
}

#[test]
fn stream_seed_is_deterministic() -> Result<(), String> {
    for suite in 0u32..16 {
        for class in 0u32..16 {
            expect_eq(
                "repeat",
                &stream_seed(suite, class),
                &stream_seed(suite, class),
            )?;
        }
    }
    Ok(())
}

#[test]
fn stream_seed_is_distinct_over_a_grid() -> Result<(), String> {
    let mut v: Vec<u64> = Vec::new();
    let suites = (0u32..64).chain([1 << 16, u32::MAX - 1, u32::MAX]);
    for suite in suites {
        let classes = (0u32..64).chain([1 << 16, u32::MAX - 1, u32::MAX]);
        for class in classes {
            v.push(stream_seed(suite, class));
        }
    }
    let n = v.len();
    v.sort_unstable();
    v.dedup();
    expect_eq("distinct stream seeds on a 67 x 67 grid", &v.len(), &n)
}

#[test]
fn suite_table_is_consistent() -> Result<(), String> {
    let mut ids: Vec<u32> = SUITES.iter().map(|s| s.suite_id).collect();
    let mut names: Vec<&str> = SUITES.iter().map(|s| s.name).collect();
    let n = SUITES.len();
    ids.sort_unstable();
    ids.dedup();
    names.sort_unstable();
    names.dedup();
    expect_eq("distinct suite ids", &ids.len(), &n)?;
    expect_eq("distinct suite names", &names.len(), &n)?;
    for s in SUITES.iter() {
        for class in [CASE_CLASS_GATE, CASE_CLASS_DEEP] {
            expect_eq(
                s.name,
                &seed_for(s.name, class),
                &Some(stream_seed(s.suite_id, class)),
            )?;
        }
    }
    expect_eq(
        "unknown suite",
        &seed_for("no_such_suite", CASE_CLASS_GATE),
        &None,
    )?;
    if CASE_CLASS_GATE == CASE_CLASS_DEEP {
        return Err(String::from("gate and deep case classes must differ"));
    }
    Ok(())
}

#[test]
fn constants_are_recorded_and_zf_pin_is_a_placeholder_or_a_sha() -> Result<(), String> {
    expect_eq("BASE", &BASE, &0x5A45_4E49_5442_4F55_u64)?;
    expect_eq("GENERATOR_VERSION", &GENERATOR_VERSION, &5u32)?;
    let is_sha = ZF_PIN.len() == 40
        && ZF_PIN
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b));
    if ZF_PIN == "pending" || is_sha {
        Ok(())
    } else {
        Err(format!(
            "ZF_PIN must be \"pending\" or a 40-digit lowercase hex SHA, got {:?}",
            ZF_PIN
        ))
    }
}
