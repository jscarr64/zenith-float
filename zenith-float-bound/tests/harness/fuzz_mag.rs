//! COVERS: Mag::add_up, Mag::div_lower, Mag::div_up, Mag::mul_lower, Mag::mul_up, Mag::sqrt_lower, Mag::sqrt_up, Mag::sub_lower
//! S3/S4 fuzz scaffold for the real `Mag` (TEST_PLAN section 7, 13; Paul's ruling: Harper owns this module, Spence writes the property checkers
//! in `fuzz_mag_props.rs` and adds them to `PROPS` below). What is here: the seeded, deterministic case draw (S2 generators: edge-weighted
//! exponents, mantissa patterns, controlled exponent gaps, zero and `Infinity` operands), the runner (`run`) with `ZB_FUZZ_SEED`, `ZB_FUZZ_CASES`
//! and `ZB_REPLAY` from tests/common/replay.rs, the shrinker hookup (tests/common/shrink.rs), and the failure report. A failing property returns
//! an `Err` that is ONE replay line (`suite=l0a_mag seed=<hex> case=<n> op=<property> prec=<p> rm=none a=<enc> b=<enc>`, the case
//! already shrunk), then ` -- `, then why. Spence's rng-driven checkers (fuzz_mag_props.rs) plug in through the property-runner contract:
//! `run_cfg` / `run_property` / `config_for` over any failure type implementing `CaseFailure` (same report format, per-property streams,
//! `ZB_FUZZ_DEEP` deep stream, non-vacuity check on outcome-class counts). Re-run it with `ZB_REPLAY=<seed>:<case> cargo test --test harness fuzz_mag_gate`.
//! The gate test `fuzz_mag_gate` runs `GATE_CASES` cases (1,000,000) by default; `tools/fuzz_deep.sh` runs the deep counts outside the gate.
//! Only properties that follow from CONTRACT 14.1/14.2 (lower <= exact <= upper) are built in, as a smoke set; nothing proposed-only (CONTRACT
//! 10.2a-d, O13-C) is asserted. The properties take an `Ops` table, so a self-test can hand them a deliberately wrong op and require rejection.

use crate::common::expect_eq;
use crate::common::generators::{
    draw_exponent, draw_exponent_pair, draw_mantissa, draw_precision, forced_exponents, ExpDomain,
    PrecSource, ZF_EXPONENT_MAX, ZF_EXPONENT_MIN,
};
use crate::common::modular;
use crate::common::replay::{
    config_from_env, parse_config, parse_hex_u64, FuzzConfig, ReplayLine, ReplayTarget,
};
use crate::common::rng::splitmix64;
use crate::common::rng::SplitMix64;
use crate::common::seeds::{seed_for, CASE_CLASS_DEEP, CASE_CLASS_GATE};
use crate::common::shrink::{shrink, ShrinkCase, SHRINK_STEP_CAP};
use zenith_bound::mag::Mag;

/// Cases of the gate run (override with `ZB_FUZZ_CASES`).
const GATE_CASES: u64 = 1_000_000;

// ---- cases ----------------------------------------------------------------------------------------------------------------

/// One operand as plain data.
#[derive(Debug, PartialEq, Eq)]
enum MagArg {
    Zero,
    Infinity,
    Finite { exp: i64, mant: u64 },
}

/// A finite `Mag` from plain integers (the operand is data, not another `Mag`).
fn finite_from_parts(mant: u64, exp: i64) -> Result<Mag, String> {
    Mag::finite(mant, exp).map_err(|e| format!("Mag::finite({:#x}, {}) failed: {:?}", mant, exp, e))
}

/// The replay target a failure line names.
fn target_of(seed: u64, case: u64) -> ReplayTarget {
    ReplayTarget { seed, case }
}

impl MagArg {
    fn to_mag(&self) -> Result<Mag, String> {
        match self {
            MagArg::Zero => Ok(Mag::zero()),
            MagArg::Infinity => Ok(Mag::infinity()),
            MagArg::Finite { exp, mant } => finite_from_parts(*mant, *exp),
        }
    }

    /// `z`, `inf`, or `e<exp>m<mant as 16 hex digits>`; no whitespace.
    fn encode(&self) -> String {
        match self {
            MagArg::Zero => String::from("z"),
            MagArg::Infinity => String::from("inf"),
            MagArg::Finite { exp, mant } => format!("e{}m{:016x}", exp, mant),
        }
    }

    fn decode(text: &str) -> Result<MagArg, String> {
        match text {
            "z" => Ok(MagArg::Zero),
            "inf" => Ok(MagArg::Infinity),
            _ => {
                let body = text
                    .strip_prefix('e')
                    .ok_or_else(|| format!("operand {:?} is not z, inf or e<exp>m<hex>", text))?;
                let (exp_text, mant_text) = body
                    .split_once('m')
                    .ok_or_else(|| format!("operand {:?} has no 'm'", text))?;
                let exp = exp_text
                    .parse::<i64>()
                    .map_err(|_| format!("operand {:?}: bad exponent", text))?;
                if mant_text.len() != 16 {
                    return Err(format!(
                        "operand {:?}: the mantissa needs 16 hex digits",
                        text
                    ));
                }
                Ok(MagArg::Finite {
                    exp,
                    mant: parse_hex_u64("mantissa", mant_text)?,
                })
            }
        }
    }

    fn with_exp(&self, replacement: Option<&i64>) -> MagArg {
        match (self, replacement) {
            (MagArg::Finite { mant, .. }, Some(e)) => MagArg::Finite {
                exp: *e,
                mant: *mant,
            },
            (MagArg::Finite { exp, mant }, None) => MagArg::Finite {
                exp: *exp,
                mant: *mant,
            },
            (MagArg::Zero, _) => MagArg::Zero,
            (MagArg::Infinity, _) => MagArg::Infinity,
        }
    }

    fn with_mant(&self, replacement: Option<&u64>) -> MagArg {
        match (self, replacement) {
            (MagArg::Finite { exp, .. }, Some(m)) => MagArg::Finite {
                exp: *exp,
                mant: *m,
            },
            (MagArg::Finite { exp, mant }, None) => MagArg::Finite {
                exp: *exp,
                mant: *mant,
            },
            (MagArg::Zero, _) => MagArg::Zero,
            (MagArg::Infinity, _) => MagArg::Infinity,
        }
    }
}

/// A drawn case: two operands and the precision that steered the exponent-gap draw (a `Mag` has no precision; it keeps the replay line uniform).
#[derive(Debug, PartialEq, Eq)]
struct MagCase {
    a: MagArg,
    b: MagArg,
    prec: u32,
}

impl MagCase {
    fn mags(&self) -> Result<(Mag, Mag), String> {
        Ok((self.a.to_mag()?, self.b.to_mag()?))
    }
}

impl ShrinkCase for MagCase {
    fn exponents(&self) -> Vec<i64> {
        let mut out = Vec::new();
        for arg in [&self.a, &self.b] {
            if let MagArg::Finite { exp, .. } = arg {
                out.push(*exp);
            }
        }
        out
    }

    fn with_exponents(&self, exponents: &[i64]) -> MagCase {
        // finite operands take the replacements in order; a zero or Infinity operand consumes none
        let mut it = exponents.iter();
        let a = match &self.a {
            MagArg::Finite { .. } => self.a.with_exp(it.next()),
            MagArg::Zero => MagArg::Zero,
            MagArg::Infinity => MagArg::Infinity,
        };
        let b = match &self.b {
            MagArg::Finite { .. } => self.b.with_exp(it.next()),
            MagArg::Zero => MagArg::Zero,
            MagArg::Infinity => MagArg::Infinity,
        };
        MagCase {
            a,
            b,
            prec: self.prec,
        }
    }

    fn prec(&self) -> u32 {
        self.prec
    }

    fn with_prec(&self, prec: u32) -> MagCase {
        MagCase {
            a: self.a.with_exp(None),
            b: self.b.with_exp(None),
            prec,
        }
    }

    fn mantissas(&self) -> Vec<u64> {
        let mut out = Vec::new();
        for arg in [&self.a, &self.b] {
            if let MagArg::Finite { mant, .. } = arg {
                out.push(*mant);
            }
        }
        out
    }

    fn with_mantissas(&self, mantissas: &[u64]) -> MagCase {
        let mut it = mantissas.iter();
        let a = match &self.a {
            MagArg::Finite { .. } => self.a.with_mant(it.next()),
            MagArg::Zero => MagArg::Zero,
            MagArg::Infinity => MagArg::Infinity,
        };
        let b = match &self.b {
            MagArg::Finite { .. } => self.b.with_mant(it.next()),
            MagArg::Zero => MagArg::Zero,
            MagArg::Infinity => MagArg::Infinity,
        };
        MagCase {
            a,
            b,
            prec: self.prec,
        }
    }
}

/// Draws one case from `rng`. Everything is drawn in a fixed order whatever the operand kinds are, so the stream position never depends on a kind:
/// precision, kind of `a` (1/32 zero, 1/32 Infinity), the exponent of `a`, its mantissa, kind of `b`, the gap pair for `b`'s exponent, its mantissa.
fn draw_case(rng: &mut SplitMix64) -> Result<MagCase, String> {
    let prec = draw_precision(rng, &PrecSource::Fast)?;
    let kind_a = rng.below(32);
    let exp_a = draw_exponent(rng, &ExpDomain::MAG, 0, prec)?.value;
    let mant_a = draw_mantissa(rng, 64)?.value;
    let kind_b = rng.below(32);
    // the partner exponent keeps a controlled gap to `a`'s; where the gap would leave the Mag range, `b` shares `a`'s exponent
    let exp_b = match draw_exponent_pair(rng, exp_a, prec, &ExpDomain::MAG) {
        Ok(((_, partner), _)) => partner,
        Err(_) => exp_a,
    };
    let mant_b = draw_mantissa(rng, 64)?.value;
    let arg = |kind: u64, exp: i64, mant: u64| match kind {
        0 => MagArg::Zero,
        1 => MagArg::Infinity,
        _ => MagArg::Finite { exp, mant },
    };
    Ok(MagCase {
        a: arg(kind_a, exp_a, mant_a),
        b: arg(kind_b, exp_b, mant_b),
        prec,
    })
}

// ---- operations and properties ----------------------------------------------------------------------------------------------

/// The `Mag` operations under test. `Ops::lib()` is the library; a self-test builds a deliberately wrong table.
struct Ops {
    add_up: fn(&Mag, &Mag) -> Mag,
    sub_lower: fn(&Mag, &Mag) -> Mag,
    mul_up: fn(&Mag, &Mag) -> Mag,
    mul_lower: fn(&Mag, &Mag) -> Mag,
    div_up: fn(&Mag, &Mag) -> Mag,
    div_lower: fn(&Mag, &Mag) -> Mag,
    sqrt_up: fn(&Mag) -> Mag,
    sqrt_lower: fn(&Mag) -> Mag,
}

impl Ops {
    fn lib() -> Ops {
        Ops {
            add_up: Mag::add_up,
            sub_lower: Mag::sub_lower,
            mul_up: Mag::mul_up,
            mul_lower: Mag::mul_lower,
            div_up: Mag::div_up,
            div_lower: Mag::div_lower,
            sqrt_up: Mag::sqrt_up,
            sqrt_lower: Mag::sqrt_lower,
        }
    }
}

/// A named property: `Err(why)` on a violation. Spence's checkers join `PROPS` with the same shape.
struct FuzzMagProp {
    name: &'static str,
    check: fn(&Ops, &Mag, &Mag) -> Result<(), String>,
}

/// CONTRACT 14.1/14.2: a round-down result never exceeds the round-up result of the same exact value (zero is least, `Infinity` greatest).
fn lower_not_above_upper(ops: &Ops, a: &Mag, b: &Mag) -> Result<(), String> {
    let pairs = [
        ("mul", (ops.mul_lower)(a, b), (ops.mul_up)(a, b)),
        ("div", (ops.div_lower)(a, b), (ops.div_up)(a, b)),
        (
            "sub_lower_vs_add_up",
            (ops.sub_lower)(a, b),
            (ops.add_up)(a, b),
        ),
        ("sqrt", (ops.sqrt_lower)(a), (ops.sqrt_up)(a)),
    ];
    for (name, lower, upper) in pairs.iter() {
        if lower > upper {
            return Err(format!(
                "{}: lower {:?} is above upper {:?}",
                name, lower, upper
            ));
        }
    }
    Ok(())
}

/// CONTRACT 14.1: the round-up sum of two nonnegative values is at least each of them.
fn add_up_dominates_operands(ops: &Ops, a: &Mag, b: &Mag) -> Result<(), String> {
    let sum = (ops.add_up)(a, b);
    if sum < *a {
        return Err(format!("add_up {:?} is below a = {:?}", sum, a));
    }
    if sum < *b {
        return Err(format!("add_up {:?} is below b = {:?}", sum, b));
    }
    Ok(())
}

/// The built-in smoke set. `fuzz_mag_props.rs` (Spence, L0a-5) appends the contract properties: enclosure against the independent model,
/// monotonicity, W1 translation.
const PROPS: [FuzzMagProp; 2] = [
    FuzzMagProp {
        name: "lower_not_above_upper",
        check: lower_not_above_upper,
    },
    FuzzMagProp {
        name: "add_up_dominates_operands",
        check: add_up_dominates_operands,
    },
];

// ---- runner -----------------------------------------------------------------------------------------------------------------

/// What a passing run did: `cases` drawn ones, then `forced` cases from `forced_cases`; `checks` counts both.
#[derive(Debug, PartialEq, Eq)]
struct RunStats {
    forced: u64,
    cases: u64,
    checks: u64,
}

/// Precision of the forced cases.
const FORCED_PREC: u32 = 53;

/// The forced cases: every ordered pair `(a, b)` of `forced_exponents` (the zf bounds and the cosmology-scale exponents), each with the
/// mantissas 2^63 and all ones on both sides (20 * 20 * 4 = 1,600 cases). They are run after the random cases of every non-replay run (so a
/// random failure still reports its replay line first),
/// so these exponents are always exercised, not only when a band happens to draw them.
fn forced_cases() -> Vec<MagCase> {
    let exps = forced_exponents();
    let mut out: Vec<MagCase> = Vec::new();
    for ea in exps.iter() {
        for eb in exps.iter() {
            for ma in [TOP_MANT, u64::MAX] {
                for mb in [TOP_MANT, u64::MAX] {
                    out.push(MagCase {
                        a: MagArg::Finite { exp: *ea, mant: ma },
                        b: MagArg::Finite { exp: *eb, mant: mb },
                        prec: FORCED_PREC,
                    });
                }
            }
        }
    }
    out
}

const TOP_MANT: u64 = 1 << 63;

fn fails(ops: &Ops, prop: &FuzzMagProp, case: &MagCase) -> bool {
    match case.mags() {
        Ok((a, b)) => (prop.check)(ops, &a, &b).is_err(),
        Err(_) => false,
    }
}

/// The report of a failing case: the replay line of the shrunk case, then the reason.
fn report(
    ops: &Ops,
    prop: &FuzzMagProp,
    base_seed: u64,
    case_index: u64,
    case: MagCase,
) -> Result<String, String> {
    let shrunk = shrink(case, &|c| fails(ops, prop, c), SHRINK_STEP_CAP)?;
    let (a, b) = shrunk.case.mags()?;
    let why = match (prop.check)(ops, &a, &b) {
        Err(why) => why,
        Ok(()) => String::from("(the shrunk case no longer fails: shrinker bug)"),
    };
    let line = ReplayLine {
        suite: String::from(SUITE),
        seed: base_seed,
        case: case_index,
        op: String::from(prop.name),
        prec: usize::try_from(shrunk.case.prec).unwrap_or(usize::MAX),
        rm: String::from("none"),
        a: shrunk.case.a.encode(),
        b: shrunk.case.b.encode(),
    };
    Ok(format!("{} -- {}", line.format_line()?, why))
}

/// Runs every property on every case of `cfg` (`0..cases`, or just the replayed case), then, outside a replay, on every forced case.
/// `Err` is the report of the first failure: a replay line for a drawn case, `forced case <i> ...` for a forced one.
fn run(ops: &Ops, props: &[FuzzMagProp], cfg: &FuzzConfig) -> Result<RunStats, String> {
    let base_seed = match &cfg.replay {
        Some(t) => t.seed,
        None => cfg.seed,
    };
    let mut forced = 0u64;
    let mut cases = 0u64;
    let mut checks = 0u64;
    for index in cfg.case_range()? {
        let mut rng = SplitMix64::new(cfg.seed_of(index));
        let case = draw_case(&mut rng)?;
        let (a, b) = case.mags()?;
        for prop in props.iter() {
            checks += 1;
            if (prop.check)(ops, &a, &b).is_err() {
                return Err(report(ops, prop, base_seed, index, case)?);
            }
        }
        cases += 1;
    }
    if cfg.replay.is_none() {
        for (i, case) in forced_cases().into_iter().enumerate() {
            let (a, b) = case.mags()?;
            for prop in props.iter() {
                checks += 1;
                if let Err(why) = (prop.check)(ops, &a, &b) {
                    return Err(format!(
                        "forced case {} op={} prec={} a={} b={} -- {}",
                        i,
                        prop.name,
                        case.prec,
                        case.a.encode(),
                        case.b.encode(),
                        why
                    ));
                }
            }
            forced += 1;
        }
    }
    Ok(RunStats {
        forced,
        cases,
        checks,
    })
}

fn gate_seed() -> Result<u64, String> {
    seed_for("l0a_mag", CASE_CLASS_GATE)
        .ok_or_else(|| String::from("suite l0a_mag is not in the seeds table"))
}

// ---- property-runner contract (Spence's rng-driven checkers in fuzz_mag_props.rs plug in here) ---------------------------------

/// Outcome-class bits a property returns for a passing case; `Counts` is indexed by bit position.
pub const OUT_FINITE: u8 = 1;
pub const OUT_INFINITY: u8 = 2;
pub const OUT_ZERO: u8 = 4;
pub const OUT_EXACT: u8 = 8;

/// Per-class pass counts `[finite, Infinity, zero, exact]`.
pub type Counts = [u64; 4];

/// The suite name in replay lines and the seeds table.
const SUITE: &str = "l0a_mag";

/// Offset of the per-property streams from the suite seed.
const PROP_BASE: u64 = 1000;

/// Below this many cases (or in a replay) the non-vacuity check is skipped.
pub const NON_VACUITY_CASES: u64 = 100_000;

/// What a failed case reports: the replay-line fields and the reason. Spence's `Bad` implements this in the post-merge glue.
pub trait CaseFailure {
    fn op(&self) -> &str;
    fn prec(&self) -> usize;
    fn a(&self) -> &str;
    fn b(&self) -> &str;
    fn why(&self) -> &str;
}

/// Base seed of property `prop_index`, gate or deep stream (no wrapping arithmetic: tests/common/modular.rs).
fn property_seed(prop_index: u64, deep: bool) -> Result<u64, String> {
    let class = if deep {
        CASE_CLASS_DEEP
    } else {
        CASE_CLASS_GATE
    };
    let s = seed_for(SUITE, class)
        .ok_or_else(|| String::from("suite l0a_mag is not in the seeds table"))?;
    Ok(splitmix64(modular::add_u64(
        modular::add_u64(s, PROP_BASE),
        prop_index,
    )))
}

/// The configuration of property `prop_index` at `cases` cases: `ZB_FUZZ_DEEP` picks the deep stream; the environment overrides.
pub fn config_for(prop_index: u64, cases: u64) -> Result<FuzzConfig, String> {
    let deep = std::env::var_os("ZB_FUZZ_DEEP").is_some();
    config_from_env(property_seed(prop_index, deep)?, cases)
}

/// Runs the cases of `cfg` for one rng-driven property and counts outcome classes. `Err` is `<replay line> -- <why>` of the first failure.
pub fn run_cfg<E: CaseFailure>(
    cfg: &FuzzConfig,
    f: &dyn Fn(&mut SplitMix64) -> Result<u8, E>,
) -> Result<Counts, String> {
    let mut counts: Counts = [0; 4];
    for case in cfg.case_range()? {
        let mut rng = SplitMix64::new(cfg.seed_of(case));
        match f(&mut rng) {
            Ok(flags) => {
                for (bit, slot) in counts.iter_mut().enumerate() {
                    if (flags >> bit) & 1 == 1 {
                        *slot += 1;
                    }
                }
            }
            Err(bad) => {
                let line = ReplayLine {
                    suite: String::from(SUITE),
                    seed: match &cfg.replay {
                        Some(t) => t.seed,
                        None => cfg.seed,
                    },
                    case,
                    op: String::from(bad.op()),
                    prec: bad.prec(),
                    rm: String::from("none"),
                    a: String::from(bad.a()),
                    b: String::from(bad.b()),
                };
                return Err(format!("{} -- {}", line.format_line()?, bad.why()));
            }
        }
    }
    Ok(counts)
}

/// One property at `cases` (environment overrides); `need` lists outcome classes that must occur, so a property cannot pass vacuously.
pub fn run_property<E: CaseFailure>(
    prop_index: u64,
    label: &str,
    need: &[(u8, &str)],
    cases: u64,
    f: &dyn Fn(&mut SplitMix64) -> Result<u8, E>,
) -> Result<(), String> {
    let cfg = config_for(prop_index, cases)?;
    let counts = run_cfg(&cfg, f)?;
    eprintln!(
        "fuzz_mag property {} ({}): cases={} [finite, Infinity, zero, exact] = {:?}",
        prop_index, label, cfg.cases, counts
    );
    if cfg.replay.is_none() && cfg.cases >= NON_VACUITY_CASES {
        for (bit, name) in need {
            let idx =
                usize::try_from(bit.trailing_zeros()).map_err(|_| String::from("bit index"))?;
            if counts.get(idx).map_or(0, |n| *n) == 0 {
                return Err(format!(
                    "property {} ({}): no case had the outcome class `{}` in {} cases ({:?}), the property is vacuous",
                    prop_index, label, name, cfg.cases, counts
                ));
            }
        }
    }
    Ok(())
}

// ---- tests ------------------------------------------------------------------------------------------------------------------

/// The gate run: `GATE_CASES` seeded cases (`ZB_FUZZ_SEED`, `ZB_FUZZ_CASES`, `ZB_REPLAY` override).
#[test]
fn fuzz_mag_gate() -> Result<(), String> {
    let cfg = config_from_env(gate_seed()?, GATE_CASES)?;
    let stats = run(&Ops::lib(), &PROPS, &cfg)?;
    let want_cases = u64::try_from(cfg.case_range()?.count()).map_err(|_| String::from("count"))?;
    let want_forced = if cfg.replay.is_none() { 1_600 } else { 0 };
    expect_eq("forced cases run", &stats.forced, &want_forced)?;
    expect_eq("cases run", &stats.cases, &want_cases)?;
    expect_eq(
        "checks run",
        &stats.checks,
        &((want_cases + want_forced) * 2),
    )
}

#[test]
fn the_ops_table_is_the_library() -> Result<(), String> {
    let ops = Ops::lib();
    let a = Mag::finite(0xC000_0000_0000_0000, 3).map_err(|e| format!("{:?}", e))?;
    let b = Mag::finite(0x8000_0000_0000_0001, -70).map_err(|e| format!("{:?}", e))?;
    expect_eq("add_up", &(ops.add_up)(&a, &b), &a.add_up(&b))?;
    expect_eq("sub_lower", &(ops.sub_lower)(&a, &b), &a.sub_lower(&b))?;
    expect_eq("mul_up", &(ops.mul_up)(&a, &b), &a.mul_up(&b))?;
    expect_eq("mul_lower", &(ops.mul_lower)(&a, &b), &a.mul_lower(&b))?;
    expect_eq("div_up", &(ops.div_up)(&a, &b), &a.div_up(&b))?;
    expect_eq("div_lower", &(ops.div_lower)(&a, &b), &a.div_lower(&b))?;
    expect_eq("sqrt_up", &(ops.sqrt_up)(&a), &a.sqrt_up())?;
    expect_eq("sqrt_lower", &(ops.sqrt_lower)(&a), &a.sqrt_lower())
}

#[test]
fn draws_are_deterministic_and_cover_every_operand_kind() -> Result<(), String> {
    let mut x = SplitMix64::new(gate_seed()?);
    let mut y = SplitMix64::new(gate_seed()?);
    let (mut zero, mut inf, mut finite, mut cross_zero_inf) = (0u32, 0u32, 0u32, 0u32);
    for _ in 0..20_000 {
        let (cx, cy) = (draw_case(&mut x)?, draw_case(&mut y)?);
        expect_eq("same stream, same case", &cx, &cy)?;
        for arg in [&cx.a, &cx.b] {
            match arg {
                MagArg::Zero => zero += 1,
                MagArg::Infinity => inf += 1,
                MagArg::Finite { mant, .. } => {
                    finite += 1;
                    if mant & 0x8000_0000_0000_0000 == 0 {
                        return Err(format!("non-canonical mantissa {:#x}", mant));
                    }
                }
            }
        }
        if cx.a == MagArg::Zero && cx.b == MagArg::Infinity {
            cross_zero_inf += 1;
        }
    }
    // 40,000 operands at 1/32 each: about 1,250 zeros and 1,250 Infinity; the bounds are loose, a wrong draw table is not
    if !(800..=1700).contains(&zero) || !(800..=1700).contains(&inf) || finite < 36_000 {
        return Err(format!(
            "operand kinds off: zero {}, inf {}, finite {}",
            zero, inf, finite
        ));
    }
    if cross_zero_inf == 0 {
        return Err(String::from("no (zero, Infinity) pair in 20,000 draws"));
    }
    Ok(())
}

#[test]
fn operands_encode_and_decode_round_trip_and_bad_text_is_refused() -> Result<(), String> {
    for arg in [
        MagArg::Zero,
        MagArg::Infinity,
        MagArg::Finite {
            exp: i64::MIN,
            mant: 0x8000_0000_0000_0000,
        },
        MagArg::Finite {
            exp: -5,
            mant: u64::MAX,
        },
        MagArg::Finite {
            exp: i64::MAX,
            mant: 0xAAAA_AAAA_AAAA_AAAB,
        },
    ] {
        let text = arg.encode();
        if text.chars().any(char::is_whitespace) {
            return Err(format!("encoding has whitespace: {:?}", text));
        }
        expect_eq(&text, &MagArg::decode(&text)?, &arg)?;
    }
    for bad in [
        "",
        "Z",
        "e5",
        "e5m",
        "e5m8000",
        "exm8000000000000000",
        "5m8000000000000000",
        "e5mg000000000000000",
    ] {
        if MagArg::decode(bad).is_ok() {
            return Err(format!("decode accepted {:?}", bad));
        }
    }
    Ok(())
}

/// A deliberately wrong table: `mul_lower` is the round-up product plus one more `Mag` unit above, so `lower > upper` whenever the product is finite.
fn broken_mul_lower(a: &Mag, b: &Mag) -> Mag {
    a.mul_up(b)
        .add_up(&Mag::finite(0x8000_0000_0000_0000, 1).unwrap_or_else(|_| Mag::infinity()))
}

fn broken_ops() -> Ops {
    Ops {
        mul_lower: broken_mul_lower,
        ..Ops::lib()
    }
}

#[test]
fn a_wrong_op_is_caught_with_a_replay_line_that_reproduces_it() -> Result<(), String> {
    let cfg = FuzzConfig {
        seed: gate_seed()?,
        cases: 5_000,
        replay: None,
    };
    let msg = match run(&broken_ops(), &PROPS, &cfg) {
        Ok(s) => return Err(format!("the wrong mul_lower passed {:?}", s)),
        Err(m) => m,
    };
    let (first, rest) = msg
        .split_once(" -- ")
        .ok_or_else(|| format!("report has no reason: {:?}", msg))?;
    if !rest.starts_with("mul: lower ") {
        return Err(format!("unexpected reason: {:?}", rest));
    }
    let line = ReplayLine::parse_line(first)?;
    expect_eq("suite", &line.suite.as_str(), &SUITE)?;
    expect_eq("op", &line.op.as_str(), &"lower_not_above_upper")?;
    expect_eq("seed", &line.seed, &cfg.seed)?;
    // the replay of (seed, case) re-draws the original case and fails the same way, with the same shrunk report
    let replay = FuzzConfig {
        seed: 0,
        cases: 1,
        replay: Some(target_of(line.seed, line.case)),
    };
    let again = match run(&broken_ops(), &PROPS, &replay) {
        Ok(s) => return Err(format!("the replay of the failing case passed {:?}", s)),
        Err(m) => m,
    };
    expect_eq("replay report", &again, &msg)?;
    // the shrunk operands decode, still fail, and are minimal in the shrinker's sense: a finite one has the mantissa 2^63 or all ones
    let (a, b) = (MagArg::decode(&line.a)?, MagArg::decode(&line.b)?);
    let shrunk = MagCase {
        a,
        b,
        prec: u32::try_from(line.prec)
            .map_err(|_| format!("replay prec {} is not a fuzz precision", line.prec))?,
    };
    if !fails(&broken_ops(), &PROPS[0], &shrunk) {
        return Err(format!("the shrunk case does not fail: {}", first));
    }
    // and the real library passes the very same case
    if fails(&Ops::lib(), &PROPS[0], &shrunk) {
        return Err(format!("the library fails the shrunk case: {}", first));
    }
    Ok(())
}

/// `mul_lower` that is wrong (above `mul_up`) only for one pair of operand exponents.
fn wrong_only_at(a: &Mag, b: &Mag, ea: i64, eb: i64) -> Mag {
    if a.exp() == Some(ea) && b.exp() == Some(eb) {
        broken_mul_lower(a, b)
    } else {
        a.mul_lower(b)
    }
}

fn wrong_at_the_zf_bounds(a: &Mag, b: &Mag) -> Mag {
    wrong_only_at(a, b, ZF_EXPONENT_MAX - 1, ZF_EXPONENT_MIN + 1)
}

fn wrong_at_cosmology_scale(a: &Mag, b: &Mag) -> Mag {
    wrong_only_at(a, b, 405, -1000)
}

#[test]
fn a_wrong_op_only_at_the_zf_bounds_or_cosmology_scale_is_caught_by_the_forced_cases(
) -> Result<(), String> {
    // 5,000 random cases essentially never draw these exact exponent pairs (MAG edges are the i64 extremes; 405 and -1000 are not near a
    // power of two), so only the forced cases catch them, and they do in every run
    let cfg = FuzzConfig {
        seed: gate_seed()?,
        cases: 5_000,
        replay: None,
    };
    let zf_bounds: fn(&Mag, &Mag) -> Mag = wrong_at_the_zf_bounds;
    for (label, f) in [
        ("zf bounds", zf_bounds),
        ("cosmology scale", wrong_at_cosmology_scale),
    ] {
        let ops = Ops {
            mul_lower: f,
            ..Ops::lib()
        };
        match run(&ops, &PROPS, &cfg) {
            Ok(s) => return Err(format!("{}: the wrong mul_lower passed {:?}", label, s)),
            Err(m) if m.starts_with("forced case ") && m.contains("op=lower_not_above_upper") => {}
            Err(m) => return Err(format!("{}: unexpected report {:?}", label, m)),
        }
    }
    Ok(())
}

#[test]
fn forced_cases_cover_every_pair_of_forced_exponents() -> Result<(), String> {
    let cases = forced_cases();
    expect_eq("forced case count", &cases.len(), &1_600)?;
    let exps = forced_exponents();
    for ea in exps.iter() {
        for eb in exps.iter() {
            let hit = cases.iter().any(|c| {
                c.a == MagArg::Finite {
                    exp: *ea,
                    mant: u64::MAX,
                } && c.b
                    == MagArg::Finite {
                        exp: *eb,
                        mant: TOP_MANT,
                    }
            });
            if !hit {
                return Err(format!("no forced case at exponents ({}, {})", ea, eb));
            }
        }
    }
    Ok(())
}

#[test]
fn a_wrong_add_up_is_caught_by_the_dominance_property() -> Result<(), String> {
    fn forgetful(a: &Mag, b: &Mag) -> Mag {
        // forgets the smaller operand when the two exponents differ by more than 64: drops below `a` when b > a
        if a < b {
            a.add_up(&Mag::zero())
        } else {
            a.add_up(b)
        }
    }
    let ops = Ops {
        add_up: forgetful,
        ..Ops::lib()
    };
    let cfg = FuzzConfig {
        seed: gate_seed()?,
        cases: 5_000,
        replay: None,
    };
    match run(&ops, &PROPS[1..], &cfg) {
        Ok(s) => Err(format!("the forgetful add_up passed {:?}", s)),
        Err(m) if m.contains("op=add_up_dominates_operands") && m.contains(" -- add_up ") => Ok(()),
        Err(m) => Err(format!("unexpected report {:?}", m)),
    }
}

#[test]
fn the_run_is_env_driven_and_replay_runs_exactly_one_case() -> Result<(), String> {
    let lookup = |vars: &'static [(&'static str, &'static str)]| {
        move |name: &str| -> Result<Option<String>, String> {
            Ok(vars
                .iter()
                .find(|(k, _)| *k == name)
                .map(|(_, v)| String::from(*v)))
        }
    };
    let cfg = parse_config(
        &lookup(&[("ZB_FUZZ_SEED", "0x2a"), ("ZB_FUZZ_CASES", "1234")]),
        7,
        GATE_CASES,
    )?;
    expect_eq("seed", &cfg.seed, &0x2a)?;
    expect_eq("cases", &cfg.cases, &1234)?;
    let dflt = parse_config(&lookup(&[]), 7, GATE_CASES)?;
    expect_eq("default cases", &dflt.cases, &1_000_000)?;
    let one = parse_config(&lookup(&[("ZB_REPLAY", "2a:17")]), 7, GATE_CASES)?;
    let stats = run(&Ops::lib(), &PROPS, &one)?;
    expect_eq(
        "replay runs one case",
        &stats,
        &RunStats {
            forced: 0,
            cases: 1,
            checks: 2,
        },
    )?;
    expect_eq(
        "replay case seed",
        &one.seed_of(17),
        &crate::common::replay::case_seed(0x2a, 17),
    )
}

#[test]
fn the_shrinker_reaches_a_minimal_case_for_a_known_failure() -> Result<(), String> {
    // fails exactly when a's exponent is above b's: the minimum keeps that order with the simplest mantissas
    let prop = FuzzMagProp {
        name: "a_not_above_b",
        check: |_, a, b| {
            if a > b {
                Err(String::from("a above b"))
            } else {
                Ok(())
            }
        },
    };
    let start = MagCase {
        a: MagArg::Finite {
            exp: 900,
            mant: 0xAAAA_AAAA_AAAA_AAAB,
        },
        b: MagArg::Finite {
            exp: 5,
            mant: 0xC123_4567_89AB_CDEF,
        },
        prec: 129,
    };
    let s = shrink(start, &|c| fails(&Ops::lib(), &prop, c), SHRINK_STEP_CAP)?;
    expect_eq("prec", &s.case.prec, &1)?;
    for m in s.case.mantissas() {
        if m != 0x8000_0000_0000_0000 && m != u64::MAX {
            return Err(format!("mantissa not simplified: {:#x}", m));
        }
    }
    if !fails(&Ops::lib(), &prop, &s.case) {
        return Err(String::from("the shrunk case does not fail"));
    }
    Ok(())
}

// ---- self-tests of the property-runner contract (run_cfg, run_property, config_for) -------------------------------------------

/// A planted failure for the runner self-tests (Spence's `Bad` plays this role for the real properties).
struct Planted {
    why: &'static str,
}

impl CaseFailure for Planted {
    fn op(&self) -> &str {
        "planted"
    }
    fn prec(&self) -> usize {
        5
    }
    fn a(&self) -> &str {
        "00000000000000a1e1"
    }
    fn b(&self) -> &str {
        "-"
    }
    fn why(&self) -> &str {
        self.why
    }
}

fn planted_on_three(rng: &mut SplitMix64) -> Result<u8, Planted> {
    if rng.next_u64() % 7 == 3 {
        Err(Planted {
            why: "planted failure",
        })
    } else {
        Ok(OUT_FINITE)
    }
}

#[test]
fn run_cfg_counts_outcome_classes_and_stops_at_the_first_failure_with_a_replay_line(
) -> Result<(), String> {
    let cfg = FuzzConfig {
        seed: 0xABCD,
        cases: 1000,
        replay: None,
    };
    let counts = run_cfg(&cfg, &|rng| -> Result<u8, Planted> {
        Ok(if rng.next_u64().is_multiple_of(2) {
            OUT_FINITE | OUT_EXACT
        } else {
            OUT_ZERO
        })
    })?;
    let [finite, infinity, zero, exact] = counts;
    expect_eq(
        "finite + zero = cases",
        &modular::add_u64(finite, zero),
        &1000,
    )?;
    expect_eq("exact = finite", &exact, &finite)?;
    expect_eq("no Infinity", &infinity, &0)?;
    let text = match run_cfg(&cfg, &planted_on_three) {
        Err(t) => t,
        Ok(c) => return Err(format!("the planted failure was not found: {:?}", c)),
    };
    let (line, why) = text
        .split_once(" -- ")
        .ok_or_else(|| format!("no separator in {:?}", text))?;
    expect_eq("reason", &why, &"planted failure")?;
    let parsed = ReplayLine::parse_line(line)?;
    expect_eq("seed", &parsed.seed, &0xABCD)?;
    expect_eq("op", &parsed.op.as_str(), &"planted")?;
    expect_eq("suite", &parsed.suite.as_str(), &SUITE)?;
    let replay = parse_config(
        &|name: &str| {
            Ok(if name == "ZB_REPLAY" {
                Some(format!("{:x}:{}", parsed.seed, parsed.case))
            } else {
                None
            })
        },
        0,
        1000,
    )?;
    match run_cfg(&replay, &planted_on_three) {
        Err(t) if t.contains(&format!("case={}", parsed.case)) => Ok(()),
        other => Err(format!("replay of the failing case: {:?}", other)),
    }
}

#[test]
fn run_property_refuses_a_vacuous_property_and_accepts_a_covering_one() -> Result<(), String> {
    let need = [(OUT_FINITE, "finite"), (OUT_INFINITY, "Infinity")];
    let only_finite = |_: &mut SplitMix64| -> Result<u8, Planted> { Ok(OUT_FINITE) };
    match run_property(3, "vacuous", &need, NON_VACUITY_CASES, &only_finite) {
        Err(m) if m.contains("`Infinity`") && m.contains("vacuous") => {}
        other => {
            return Err(format!(
                "a property that never reaches Infinity passed: {:?}",
                other
            ))
        }
    }
    let both = |rng: &mut SplitMix64| -> Result<u8, Planted> {
        Ok(if rng.next_u64().is_multiple_of(3) {
            OUT_INFINITY
        } else {
            OUT_FINITE
        })
    };
    run_property(3, "covering", &need, NON_VACUITY_CASES, &both)?;
    match run_property(3, "failing", &need, NON_VACUITY_CASES, &planted_on_three) {
        Err(m) if m.contains(" -- planted failure") => Ok(()),
        other => Err(format!("a failing property passed: {:?}", other)),
    }
}

#[test]
fn property_streams_are_distinct_per_index_and_class() -> Result<(), String> {
    let g0 = property_seed(0, false)?;
    let g1 = property_seed(1, false)?;
    let d0 = property_seed(0, true)?;
    if g0 == g1 || g0 == d0 || g1 == d0 {
        return Err(format!(
            "property seeds collide: {:x} {:x} {:x}",
            g0, g1, d0
        ));
    }
    expect_eq("deterministic", &property_seed(1, false)?, &g1)?;
    // known answers computed outside Rust (Python, same SplitMix64 + seeds formula): the streams keep L0a-5's numbering
    // (`splitmix64(seed_for(l0a_mag, class) + 1000 + index)`), so replay lines written against 93c8b49 stay valid
    expect_eq("gate stream 0", &g0, &0x3979_fa53_2cd0_15a8)?;
    expect_eq("gate stream 1", &g1, &0x5ec5_e017_70d3_6fab)?;
    expect_eq("deep stream 0", &d0, &0x8485_a395_1d0a_4533)?;
    let cfg = config_for(1, 77)?;
    if std::env::var_os("ZB_FUZZ_SEED").is_none() && std::env::var_os("ZB_FUZZ_DEEP").is_none() {
        expect_eq("config seed", &cfg.seed, &g1)?;
    }
    Ok(())
}
