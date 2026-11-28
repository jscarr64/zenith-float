//! The record of fixed seeds. One `BASE`, one table of suites, one derivation (`stream_seed`).
//! Frozen draw order: a change to the generator or to how suites consume draws bumps `GENERATOR_VERSION`.

use super::rng::splitmix64;

/// Base of every stream seed (ASCII "ZENITBOU"). Never changed without bumping `GENERATOR_VERSION`.
pub const BASE: u64 = 0x5A45_4E49_5442_4F55;

/// Generator version; goldens are unaffected (they come from Python), grid digests are not.
pub const GENERATOR_VERSION: u32 = 5;

/// PLACEHOLDER: the zenith-float pin SHA is not decided yet (TEST_PLAN section 15). "pending" until Jeff/Rusty set it; then a 40-digit lowercase hex SHA.
pub const ZF_PIN: &str = "pending";

/// Case class of the fixed-seed gate run.
pub const CASE_CLASS_GATE: u32 = 0;
/// Case class of the deep / nightly run (extra seeds).
pub const CASE_CLASS_DEEP: u32 = 1;

pub struct SuiteSeed {
    pub name: &'static str,
    pub suite_id: u32,
}

/// The single table of per-suite seeds (suite ids are never reused or renumbered). Names follow TEST_PLAN section 1 levels.
pub const SUITES: [SuiteSeed; 8] = [
    SuiteSeed {
        name: "common",
        suite_id: 0,
    },
    SuiteSeed {
        name: "l0a_mag",
        suite_id: 1,
    },
    SuiteSeed {
        name: "l0b_ball",
        suite_id: 2,
    },
    SuiteSeed {
        name: "l1_elementary",
        suite_id: 3,
    },
    SuiteSeed {
        name: "l2",
        suite_id: 4,
    },
    SuiteSeed {
        name: "l3_quadrature",
        suite_id: 5,
    },
    SuiteSeed {
        name: "l4_workflows",
        suite_id: 6,
    },
    SuiteSeed {
        name: "s6_matrix",
        suite_id: 17,
    },
];

/// `splitmix64(BASE + suite_id * 2^32 + case_class)`. `suite_id * 2^32 + case_class` fits `u64` for every pair of `u32`
/// inputs (it is the pair packed high/low); adding `BASE` wraps, deterministically.
pub fn stream_seed(suite_id: u32, case_class: u32) -> u64 {
    let packed = (u64::from(suite_id) << 32) | u64::from(case_class);
    splitmix64(super::modular::add_u64(BASE, packed))
}

/// The stream seed of a suite from the table, `None` for a name that is not in the table.
pub fn seed_for(name: &str, case_class: u32) -> Option<u64> {
    SUITES
        .iter()
        .find(|s| s.name == name)
        .map(|s| stream_seed(s.suite_id, case_class))
}
