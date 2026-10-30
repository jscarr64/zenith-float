//! COVERS: common::SplitMix64::new, common::SplitMix64::next_u64, common::rng::splitmix64, common::SplitMix64::below, common::expect_eq
//! Oracle: splitmix64 reference outputs (seed 0: the published values; the other seeds recomputed with an independent Python integer-mask implementation, 2026-10-15).

use crate::common::expect_eq;
use crate::common::rng::{splitmix64, SplitMix64};

fn first_outputs(seed: u64, n: usize) -> Vec<u64> {
    let mut g = SplitMix64::new(seed);
    (0..n).map(|_| g.next_u64()).collect()
}

#[test]
fn seed_zero_matches_reference_outputs() -> Result<(), String> {
    let want = vec![
        0xE220_A839_7B1D_CDAF_u64,
        0x6E78_9E6A_A1B9_65F4,
        0x06C4_5D18_8009_454F,
    ];
    expect_eq("seed 0", &first_outputs(0, 3), &want)
}

#[test]
fn other_seeds_match_independent_python_outputs() -> Result<(), String> {
    expect_eq(
        "seed 1",
        &first_outputs(1, 3),
        &vec![
            0x910A_2DEC_8902_5CC1_u64,
            0xBEEB_8DA1_658E_EC67,
            0xF893_A2EE_FB32_555E,
        ],
    )?;
    expect_eq(
        "seed 0xDEADBEEF",
        &first_outputs(0xDEAD_BEEF, 5),
        &vec![
            0x4ADF_B90F_68C9_EB9B_u64,
            0xDE58_6A31_41A1_0922,
            0x021F_BC2F_8E1C_FC1D,
            0x7466_CE73_7BE1_6790,
            0x3BFA_8764_F685_BD1C,
        ],
    )?;
    // The state wraps on the first increment for u64::MAX: wrapping arithmetic is required.
    expect_eq(
        "seed u64::MAX",
        &first_outputs(u64::MAX, 3),
        &vec![
            0xE4D9_7177_1B65_2C20_u64,
            0xE99F_F867_DBF6_82C9,
            0x382F_F84C_B272_81E9,
        ],
    )
}

#[test]
fn one_shot_equals_first_draw() -> Result<(), String> {
    for seed in [0u64, 1, 2, 0x5EED, 0xDEAD_BEEF, u64::MAX - 1, u64::MAX] {
        expect_eq(
            "one-shot",
            &splitmix64(seed),
            &SplitMix64::new(seed).next_u64(),
        )?;
    }
    Ok(())
}

#[test]
fn same_seed_same_sequence_different_seed_different_sequence() -> Result<(), String> {
    expect_eq(
        "same seed",
        &first_outputs(42, 1000),
        &first_outputs(42, 1000),
    )?;
    if first_outputs(42, 8) == first_outputs(43, 8) {
        return Err(String::from(
            "seeds 42 and 43 gave the same first 8 outputs",
        ));
    }
    Ok(())
}

#[test]
fn a_long_run_has_no_repeats() -> Result<(), String> {
    let mut v = first_outputs(7, 100_000);
    v.sort_unstable();
    let n = v.len();
    v.dedup();
    expect_eq("distinct outputs of 100000 draws", &v.len(), &n)
}

#[test]
fn below_is_in_range_deterministic_and_covers_the_range() -> Result<(), String> {
    let mut g = SplitMix64::new(9);
    let mut seen = [false; 10];
    for _ in 0..1000 {
        let v = g.below(10);
        if v >= 10 {
            return Err(format!("below(10) gave {}", v));
        }
        if let Some(f) = seen.get_mut(usize::try_from(v).map_err(|_| String::from("index"))?) {
            *f = true;
        }
    }
    expect_eq("every value of below(10) seen", &seen, &[true; 10])?;
    expect_eq("below(0)", &g.below(0), &0)?;
    expect_eq("below(1)", &g.below(1), &0)?;
    // below is a pure function of the draw: one next_u64 output scaled by n (widening multiply)
    let mut a = SplitMix64::new(0);
    let mut b = SplitMix64::new(0);
    let wide = (u128::from(b.next_u64()) * 100) >> 64;
    expect_eq("below(100) from seed 0", &u128::from(a.below(100)), &wide)
}
