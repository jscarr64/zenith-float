//! Self-tests of the L0a-5 Mag property checkers (tests/harness/fuzz_mag_props.rs): a checker must reject a wrong answer.
//! (1) the exact-comparison oracle accepts every row of the mpmath goldens and rejects every perturbed result;
//! (2) every property, run with a deliberately wrong operation in place of the real one, fails within a few thousand cases and reports a replay line;
//! (3) the runner stops at the first failure, replays exactly one case on request, and the case draws reach the range boundaries.
//! No test here depends on the correctness of the library except where it says so (the planted operations are real `Mag` operations used in the
//! wrong place, e.g. `mul_lower` as `mul_up`).

use crate::common::expect_eq;
use crate::common::generators::TOP;
use crate::common::replay::{parse_config, FuzzConfig, ReplayLine};
use crate::common::rng::SplitMix64;
use crate::fuzz_mag::run_cfg;
use crate::fuzz_mag_props::{
    self as props, check_bracket, check_down, check_up, cmp_add, cmp_div, cmp_mul, cmp_scaled,
    cmp_sqrt, cmp_sub, draw_case, next_down, next_up, Bad, Mono, Ops, Out, Rel, Src,
};
use core::cmp::Ordering;
use zenith_bound::mag::Mag;

const GOLDEN: &str = include_str!("../goldens/mag_ops.txt");

// ---- (1) the oracle against the mpmath goldens --------------------------------------------------------------------------------

fn num_i64(it: &mut core::str::SplitWhitespace) -> Result<i64, String> {
    it.next()
        .ok_or_else(|| String::from("row ends early"))?
        .parse::<i64>()
        .map_err(|e| format!("bad i64: {}", e))
}

fn num_u64(it: &mut core::str::SplitWhitespace) -> Result<u64, String> {
    it.next()
        .ok_or_else(|| String::from("row ends early"))?
        .parse::<u64>()
        .map_err(|e| format!("bad u64: {}", e))
}

fn result_of(it: &mut core::str::SplitWhitespace) -> Result<Out, String> {
    let first = it.next().ok_or_else(|| String::from("row ends early"))?;
    if first == "inf" {
        return Ok(Out::Infinity);
    }
    let exp = first
        .parse::<i64>()
        .map_err(|e| format!("bad i64: {}", e))?;
    let mant = num_u64(it)?;
    if mant == 0 {
        Ok(Out::Zero)
    } else {
        Ok(Out::Finite { exp, mant })
    }
}

/// Runs `visit(op, a, b, result)` on every golden row whose operands are finite and nonzero (`b` is `(0, 0)` for the unary ops).
fn each_golden_row(
    mut visit: impl FnMut(&str, (i64, u64), (i64, u64), Out) -> Result<(), String>,
) -> Result<usize, String> {
    let mut n = 0usize;
    for line in GOLDEN.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut it = line.split_whitespace();
        let op = it.next().ok_or_else(|| String::from("empty row"))?;
        let a = (num_i64(&mut it)?, num_u64(&mut it)?);
        let unary = op.starts_with("sqrt");
        let b = if unary {
            (0, 0)
        } else {
            (num_i64(&mut it)?, num_u64(&mut it)?)
        };
        let r = result_of(&mut it)?;
        if a.1 == 0 || (!unary && b.1 == 0) {
            continue;
        }
        visit(op, a, b, r).map_err(|e| format!("{}: {}", line, e))?;
        n += 1;
    }
    Ok(n)
}

fn check_row(op: &str, a: (i64, u64), b: (i64, u64), r: &Out) -> Result<(), String> {
    let cmp = |er: i64, mr: u64| -> Result<Ordering, String> {
        match op {
            "add_up" => cmp_add(a, b, (er, mr)),
            "mul_up" | "mul_lower" => cmp_mul(a, b, (er, mr)),
            "div_up" | "div_lower" => cmp_div(a, b, (er, mr)),
            "sqrt_up" | "sqrt_lower" => cmp_sqrt(a, (er, mr)),
            "sub_lower" => cmp_sub(a, b, (er, mr)),
            _ => Err(format!("unknown op {}", op)),
        }
    };
    match op {
        "add_up" | "mul_up" | "div_up" | "sqrt_up" => check_up(r, &cmp),
        "sub_lower" if a <= b => {
            if *r == Out::Zero {
                Ok(())
            } else {
                Err(String::from("a <= b but the golden result is not zero"))
            }
        }
        _ => check_down(r, &cmp),
    }
}

#[test]
fn the_exact_oracle_accepts_every_mpmath_golden_row() -> Result<(), String> {
    let n = each_golden_row(|op, a, b, r| check_row(op, a, b, &r))?;
    if n < 15_000 {
        return Err(format!("only {} golden rows were checked", n));
    }
    Ok(())
}

#[test]
fn the_exact_oracle_rejects_every_perturbed_golden_result() -> Result<(), String> {
    let mut perturbed = 0usize;
    each_golden_row(|op, a, b, r| {
        let mut wrong: Vec<Out> = Vec::new();
        match &r {
            Out::Finite { exp, mant } => {
                wrong.push(match next_up(*exp, *mant) {
                    Some((e, m)) => Out::Finite { exp: e, mant: m },
                    None => Out::Infinity,
                });
                wrong.push(match next_down(*exp, *mant) {
                    Some((e, m)) => Out::Finite { exp: e, mant: m },
                    None => Out::Zero,
                });
            }
            Out::Zero => wrong.push(Out::Finite {
                exp: i64::MIN,
                mant: TOP,
            }),
            Out::Infinity => wrong.push(Out::Finite {
                exp: i64::MAX,
                mant: u64::MAX,
            }),
        }
        wrong.push(Out::Zero);
        wrong.push(Out::Infinity);
        for w in wrong.iter() {
            if *w == r {
                continue;
            }
            if check_row(op, a, b, w).is_ok() {
                return Err(format!("the oracle accepted {:?} instead of {:?}", w, r));
            }
            perturbed += 1;
        }
        Ok(())
    })?;
    if perturbed < 40_000 {
        return Err(format!("only {} perturbed results were tried", perturbed));
    }
    Ok(())
}

#[test]
fn cmp_scaled_orders_values_with_different_scales() -> Result<(), String> {
    let cases: [(u128, i128, u128, i128, Ordering); 9] = [
        (1, 0, 1, 0, Ordering::Equal),
        (1, 5, 32, 0, Ordering::Equal),
        (3, 0, 1, 2, Ordering::Less),
        (5, 0, 1, 2, Ordering::Greater),
        (u128::MAX, 0, 1, 128, Ordering::Less),
        (1, 128, u128::MAX, 0, Ordering::Greater),
        (0, 100, 1, -100, Ordering::Less),
        (0, 1, 0, 2, Ordering::Equal),
        (1, -70, 1u128 << 64, -134, Ordering::Equal),
    ];
    for (x, sx, y, sy, want) in cases {
        let got = cmp_scaled(x, sx, y, sy)?;
        if got != want {
            return Err(format!(
                "{}*2^{} vs {}*2^{}: expected {:?}, got {:?}",
                x, sx, y, sy, want, got
            ));
        }
    }
    Ok(())
}

#[test]
fn next_up_and_next_down_are_inverse_neighbours_with_the_range_ends() -> Result<(), String> {
    for (e, m) in [
        (0i64, TOP),
        (7, TOP + 1),
        (-3, u64::MAX),
        (i64::MIN + 1, TOP),
        (i64::MAX - 1, u64::MAX),
    ] {
        let (ue, um) = next_up(e, m).ok_or_else(|| format!("no successor of {}:{}", e, m))?;
        let back = next_down(ue, um).ok_or_else(|| format!("no predecessor of {}:{}", ue, um))?;
        if back != (e, m) {
            return Err(format!("next_down(next_up({}:{})) = {:?}", e, m, back));
        }
    }
    if next_up(i64::MAX, u64::MAX).is_some() || next_down(i64::MIN, TOP).is_some() {
        return Err(String::from("the ends of the range have a neighbour"));
    }
    if next_up(i64::MIN, TOP) != Some((i64::MIN, TOP + 1)) {
        return Err(String::from("successor of the smallest positive Mag"));
    }
    Ok(())
}

#[test]
fn check_up_check_down_and_check_bracket_reject_wrong_results() -> Result<(), String> {
    // exact value 1.5 = (1, 0xC000..): the candidate order is against (1, 0xC000..)
    let cmp = |e: i64, m: u64| cmp_scaled(3, -1, u128::from(m), i128::from(e) - 64);
    let exact = Out::Finite {
        exp: 1,
        mant: 0xC000_0000_0000_0000,
    };
    let above = Out::Finite {
        exp: 1,
        mant: 0xC000_0000_0000_0001,
    };
    let below = Out::Finite {
        exp: 1,
        mant: 0xBFFF_FFFF_FFFF_FFFF,
    };
    for ok in [check_up(&exact, &cmp), check_down(&exact, &cmp)] {
        ok?;
    }
    for (what, r) in [
        ("up accepted a result below", check_up(&below, &cmp)),
        ("up accepted a non-smallest result", check_up(&above, &cmp)),
        ("up accepted zero", check_up(&Out::Zero, &cmp)),
        ("up accepted Infinity", check_up(&Out::Infinity, &cmp)),
        ("down accepted a result above", check_down(&above, &cmp)),
        (
            "down accepted a non-largest result",
            check_down(&below, &cmp),
        ),
        ("down accepted zero", check_down(&Out::Zero, &cmp)),
        ("down accepted Infinity", check_down(&Out::Infinity, &cmp)),
        ("bracket accepted two apart", check_bracket(&above, &below)),
    ] {
        if r.is_ok() {
            return Err(String::from(what));
        }
    }
    check_bracket(&exact, &exact)?;
    check_bracket(&above, &exact)
}

// ---- (2) every property fails with a wrong operation ---------------------------------------------------------------------------

/// Runs `f` for `cases` cases of the stream of property `idx` (no environment); `Err` is the expected result of a planted wrong op.
fn run_planted(
    idx: u64,
    cases: u64,
    f: &dyn Fn(&mut SplitMix64) -> Result<u8, Bad>,
) -> Result<(), String> {
    let cfg = FuzzConfig {
        seed: 0x5EED_0000 + idx,
        cases,
        replay: None,
    };
    match run_cfg(&cfg, f) {
        Ok(counts) => Err(format!(
            "the wrong operation survived {} cases (outcome counts {:?})",
            cases, counts
        )),
        Err(text) => {
            let (line, _) = text
                .split_once(" -- ")
                .ok_or_else(|| format!("no replay line in {:?}", text))?;
            ReplayLine::parse_line(line)?;
            Ok(())
        }
    }
}

fn wrong_sqrt_flip(a: &Mag) -> Mag {
    match a.exp() {
        Some(e) if e % 2 == 0 => a.sqrt_up(),
        _ => a.sqrt_up().sqrt_up(),
    }
}

fn wrong_smaller(x: &Mag, mid: i64, prec: usize) -> bool {
    x.smaller_than_half_ulp(mid, prec.checked_add(1).unwrap_or(prec))
}

fn wrong_smaller_inverted(x: &Mag, mid: i64, prec: usize) -> bool {
    !x.smaller_than_half_ulp(mid, prec)
}

macro_rules! planted {
    ($name:ident, $idx:expr, $cases:expr, $ops:expr, $case:path) => {
        #[test]
        fn $name() -> Result<(), String> {
            let ops = $ops;
            run_planted($idx, $cases, &|rng| $case(&ops, &mut Src::random(rng)))
        }
    };
}

planted!(
    add_up_case_rejects_sub_lower_as_add_up,
    0,
    5000,
    Ops {
        add_up: Mag::sub_lower,
        ..props::real()
    },
    props::add_up_case
);
planted!(
    add_up_case_rejects_mul_up_as_add_up,
    1,
    5000,
    Ops {
        add_up: Mag::mul_up,
        ..props::real()
    },
    props::add_up_case
);
planted!(
    sub_lower_case_rejects_add_up_as_sub_lower,
    2,
    5000,
    Ops {
        sub_lower: Mag::add_up,
        ..props::real()
    },
    props::sub_lower_case
);
planted!(
    mul_case_rejects_mul_lower_as_mul_up,
    3,
    5000,
    Ops {
        mul_up: Mag::mul_lower,
        ..props::real()
    },
    props::mul_case
);
planted!(
    mul_case_rejects_mul_up_as_mul_lower,
    4,
    5000,
    Ops {
        mul_lower: Mag::mul_up,
        ..props::real()
    },
    props::mul_case
);
planted!(
    div_case_rejects_div_lower_as_div_up,
    5,
    5000,
    Ops {
        div_up: Mag::div_lower,
        ..props::real()
    },
    props::div_case
);
planted!(
    div_case_rejects_div_up_as_div_lower,
    6,
    5000,
    Ops {
        div_lower: Mag::div_up,
        ..props::real()
    },
    props::div_case
);
planted!(
    sqrt_case_rejects_sqrt_lower_as_sqrt_up,
    7,
    5000,
    Ops {
        sqrt_up: Mag::sqrt_lower,
        ..props::real()
    },
    props::sqrt_case
);
planted!(
    sqrt_case_rejects_sqrt_up_as_sqrt_lower,
    8,
    5000,
    Ops {
        sqrt_lower: Mag::sqrt_up,
        ..props::real()
    },
    props::sqrt_case
);
planted!(
    commutative_case_rejects_sub_lower_as_add_up,
    9,
    5000,
    Ops {
        add_up: Mag::sub_lower,
        ..props::real()
    },
    props::commutative_case
);
planted!(
    commutative_case_rejects_div_up_as_mul_up,
    10,
    5000,
    Ops {
        mul_up: Mag::div_up,
        ..props::real()
    },
    props::commutative_case
);
planted!(
    identity_case_rejects_mul_up_as_add_up,
    11,
    5000,
    Ops {
        add_up: Mag::mul_up,
        ..props::real()
    },
    props::identity_case
);
planted!(
    identity_case_rejects_mul_lower_as_mul_up,
    12,
    5000,
    Ops {
        mul_up: Mag::mul_lower,
        ..props::real()
    },
    props::identity_case
);
planted!(
    identity_case_rejects_div_up_as_div_lower,
    13,
    5000,
    Ops {
        div_lower: Mag::div_up,
        ..props::real()
    },
    props::identity_case
);
planted!(
    sqrt_exact_case_rejects_sqrt_lower_as_sqrt_up,
    14,
    5000,
    Ops {
        sqrt_up: Mag::sqrt_lower,
        ..props::real()
    },
    props::sqrt_exact_case
);
planted!(
    sqrt_exact_case_rejects_sqrt_up_as_sqrt_lower,
    15,
    5000,
    Ops {
        sqrt_lower: Mag::sqrt_up,
        ..props::real()
    },
    props::sqrt_exact_case
);
planted!(
    ulp_case_rejects_half_ulp_as_ulp_mag,
    16,
    5000,
    Ops {
        ulp_mag: Mag::half_ulp,
        ..props::real()
    },
    props::ulp_case
);
planted!(
    ulp_case_rejects_ulp_mag_as_half_ulp,
    17,
    5000,
    Ops {
        half_ulp: Mag::ulp_mag,
        ..props::real()
    },
    props::ulp_case
);
planted!(
    ulp_case_rejects_a_shifted_threshold,
    18,
    5000,
    Ops {
        smaller_than_half_ulp: wrong_smaller,
        ..props::real()
    },
    props::ulp_case
);
planted!(
    ulp_case_rejects_an_inverted_predicate,
    19,
    5000,
    Ops {
        smaller_than_half_ulp: wrong_smaller_inverted,
        ..props::real()
    },
    props::ulp_case
);

fn wrong_add_zero_zero(a: &Mag, b: &Mag) -> Mag {
    if a.is_zero() && b.is_zero() {
        Mag::infinity()
    } else {
        a.add_up(b)
    }
}

fn wrong_mul_lower_infinity(a: &Mag, b: &Mag) -> Mag {
    if a.is_infinity() || b.is_infinity() {
        Mag::infinity()
    } else {
        a.mul_lower(b)
    }
}

fn wrong_div_up_zero_by_zero(a: &Mag, b: &Mag) -> Mag {
    if a.is_zero() && b.is_zero() {
        Mag::zero()
    } else {
        a.div_up(b)
    }
}

fn wrong_sub_infinity_subtrahend(a: &Mag, b: &Mag) -> Mag {
    if b.is_infinity() {
        a.add_up(&Mag::zero())
    } else {
        a.sub_lower(b)
    }
}

fn wrong_sqrt_lower_infinity(a: &Mag) -> Mag {
    if a.is_infinity() {
        Mag::infinity()
    } else {
        a.sqrt_lower()
    }
}

planted!(
    canonical_case_rejects_zero_plus_zero_as_infinity,
    40,
    5000,
    Ops {
        add_up: wrong_add_zero_zero,
        ..props::real()
    },
    props::canonical_case
);
planted!(
    canonical_case_rejects_mul_lower_of_infinity_as_infinity,
    41,
    5000,
    Ops {
        mul_lower: wrong_mul_lower_infinity,
        ..props::real()
    },
    props::canonical_case
);
planted!(
    canonical_case_rejects_zero_over_zero_as_zero,
    42,
    5000,
    Ops {
        div_up: wrong_div_up_zero_by_zero,
        ..props::real()
    },
    props::canonical_case
);
planted!(
    canonical_case_rejects_subtracting_infinity_as_a_no_op,
    43,
    5000,
    Ops {
        sub_lower: wrong_sub_infinity_subtrahend,
        ..props::real()
    },
    props::canonical_case
);
planted!(
    canonical_case_rejects_sqrt_lower_of_infinity_as_infinity,
    44,
    5000,
    Ops {
        sqrt_lower: wrong_sqrt_lower_infinity,
        ..props::real()
    },
    props::canonical_case
);

macro_rules! planted_mono {
    ($name:ident, $idx:expr, $which:expr, $ops:expr) => {
        #[test]
        fn $name() -> Result<(), String> {
            let ops = $ops;
            run_planted($idx, 20_000, &|rng| {
                props::mono_case(&ops, &$which, &mut Src::random(rng))
            })
        }
    };
}

planted_mono!(
    mono_add_up_rejects_sub_lower_as_add_up,
    20,
    Mono::AddUp,
    Ops {
        add_up: Mag::sub_lower,
        ..props::real()
    }
);
planted_mono!(
    mono_mul_up_rejects_div_up_as_mul_up,
    21,
    Mono::MulUp,
    Ops {
        mul_up: Mag::div_up,
        ..props::real()
    }
);
planted_mono!(
    mono_mul_lower_rejects_div_lower_as_mul_lower,
    22,
    Mono::MulLower,
    Ops {
        mul_lower: Mag::div_lower,
        ..props::real()
    }
);
planted_mono!(
    mono_div_up_rejects_mul_up_as_div_up,
    23,
    Mono::DivUp,
    Ops {
        div_up: Mag::mul_up,
        ..props::real()
    }
);
planted_mono!(
    mono_div_lower_rejects_mul_lower_as_div_lower,
    24,
    Mono::DivLower,
    Ops {
        div_lower: Mag::mul_lower,
        ..props::real()
    }
);
planted_mono!(
    mono_sub_lower_rejects_add_up_as_sub_lower,
    25,
    Mono::SubLower,
    Ops {
        sub_lower: Mag::add_up,
        ..props::real()
    }
);
planted_mono!(
    mono_sqrt_up_rejects_a_parity_dependent_root,
    26,
    Mono::SqrtUp,
    Ops {
        sqrt_up: wrong_sqrt_flip,
        ..props::real()
    }
);
planted_mono!(
    mono_sqrt_lower_rejects_a_parity_dependent_root,
    27,
    Mono::SqrtLower,
    Ops {
        sqrt_lower: wrong_sqrt_flip,
        ..props::real()
    }
);

#[test]
fn the_real_operations_pass_every_property_for_a_short_run() -> Result<(), String> {
    // the positive control of the planted tests above: same seeds scheme, the real `Mag`
    let ops = props::real();
    let cfg = FuzzConfig {
        seed: 0x5EED_0000,
        cases: 5000,
        replay: None,
    };
    let first: &dyn Fn(&Ops, &mut Src) -> Result<u8, Bad> = &props::add_up_case;
    for f in [
        first,
        &props::sub_lower_case,
        &props::mul_case,
        &props::div_case,
        &props::sqrt_case,
        &props::commutative_case,
        &props::identity_case,
        &props::sqrt_exact_case,
        &props::canonical_case,
        &props::ulp_case,
    ] {
        run_cfg(&cfg, &|rng| f(&ops, &mut Src::random(rng)))?;
    }
    Ok(())
}

// ---- (3) the runner and the draws ----------------------------------------------------------------------------------------------

/// The configuration of `ZB_REPLAY=<seed>:<case>`, parsed by the real `parse_config` over an injected lookup (the process environment is untouched).
fn replay_cfg(seed: u64, case: u64) -> Result<FuzzConfig, String> {
    parse_config(
        &|name: &str| {
            Ok(if name == "ZB_REPLAY" {
                Some(format!("{:x}:{}", seed, case))
            } else {
                None
            })
        },
        0,
        1000,
    )
}

fn fail_on_case_three(rng: &mut SplitMix64) -> Result<u8, Bad> {
    let first = rng.next_u64();
    if first % 7 == 3 {
        Err(Bad {
            op: "planted",
            prec: 5,
            a: String::from("00000000000000a1e1"),
            b: String::from("-"),
            msg: String::from("planted failure"),
        })
    } else {
        Ok(props::FIN)
    }
}

#[test]
fn the_runner_stops_at_the_first_failure_and_reports_a_replayable_line() -> Result<(), String> {
    let cfg = FuzzConfig {
        seed: 0xABCD,
        cases: 1000,
        replay: None,
    };
    let text = match run_cfg(&cfg, &fail_on_case_three) {
        Err(t) => t,
        Ok(c) => return Err(format!("the planted failure was not found: {:?}", c)),
    };
    let (line, msg) = text
        .split_once(" -- ")
        .ok_or_else(|| String::from("no separator"))?;
    if msg != "planted failure" {
        return Err(format!("reason {:?}", msg));
    }
    let parsed = ReplayLine::parse_line(line)?;
    if parsed.seed != 0xABCD || parsed.op != "planted" || parsed.prec != 5 {
        return Err(format!("replay line fields: {:?}", parsed));
    }
    // replay exactly that case: the same failure; replay the case before it: a pass with one case counted
    let replay = replay_cfg(parsed.seed, parsed.case)?;
    match run_cfg(&replay, &fail_on_case_three) {
        Err(t) if t.contains("planted failure") && t.contains(&format!("case={}", parsed.case)) => {
        }
        other => return Err(format!("replay of the failing case: {:?}", other)),
    }
    let earlier = match parsed.case.checked_sub(1) {
        Some(c) => c,
        None => return Ok(()),
    };
    let ok = replay_cfg(parsed.seed, earlier)?;
    let counts = run_cfg(&ok, &fail_on_case_three)?;
    if counts != [1, 0, 0, 0] {
        return Err(format!("replay of a passing case ran {:?}", counts));
    }
    Ok(())
}

#[test]
fn the_case_draws_reach_the_range_boundaries_and_every_shape() -> Result<(), String> {
    let (mut sum_hi, mut sum_lo, mut diff_hi, mut diff_lo, mut equal_exp, mut wide_gap) =
        (0u32, 0u32, 0u32, 0u32, 0u32, 0u32);
    let cfg = FuzzConfig {
        seed: 0x77,
        cases: 40_000,
        replay: None,
    };
    for case in cfg.case_range()? {
        let mut rng = SplitMix64::new(cfg.seed_of(case));
        let s = draw_case(&mut Src::random(&mut rng), &Rel::Sum)?;
        let sum = i128::from(s.ea) + i128::from(s.eb);
        if (i128::from(i64::MAX) - 3..=i128::from(i64::MAX) + 3).contains(&sum) {
            sum_hi += 1;
        }
        if (i128::from(i64::MIN) - 3..=i128::from(i64::MIN) + 3).contains(&sum) {
            sum_lo += 1;
        }
        let mut rng = SplitMix64::new(cfg.seed_of(case));
        let d = draw_case(&mut Src::random(&mut rng), &Rel::Diff)?;
        let diff = i128::from(d.ea) - i128::from(d.eb);
        if (i128::from(i64::MAX) - 3..=i128::from(i64::MAX) + 3).contains(&diff) {
            diff_hi += 1;
        }
        if (i128::from(i64::MIN) - 3..=i128::from(i64::MIN) + 3).contains(&diff) {
            diff_lo += 1;
        }
        let mut rng = SplitMix64::new(cfg.seed_of(case));
        let a = draw_case(&mut Src::random(&mut rng), &Rel::Any)?;
        if a.ea == a.eb {
            equal_exp += 1;
        }
        if i128::from(a.ea).abs_diff(i128::from(a.eb)) >= 128 {
            wide_gap += 1;
        }
        if a.ma < TOP || a.mb < TOP {
            return Err(format!(
                "a non-canonical mantissa was drawn: {:#x} {:#x}",
                a.ma, a.mb
            ));
        }
    }
    for (name, n) in [
        ("sum near i64::MAX", sum_hi),
        ("sum near i64::MIN", sum_lo),
        ("difference near i64::MAX", diff_hi),
        ("difference near i64::MIN", diff_lo),
        ("equal exponents", equal_exp),
        ("gap of 128 or more", wide_gap),
    ] {
        if n == 0 {
            return Err(format!("no draw had the shape `{}`", name));
        }
    }
    Ok(())
}

// ---- (4) the forced pairs --------------------------------------------------------------------------------------------------------

/// `mul_up` that is wrong (zero) only on the operand exponent pairs (zf MAX - 1, zf MIN + 1) and (405, -1000).
fn wrong_mul_up_only_at_forced_pairs(a: &Mag, b: &Mag) -> Mag {
    let max = (1i64 << 62) - 1;
    match (a.exp(), b.exp()) {
        (Some(x), Some(y)) if (x, y) == (max - 1, 1 - max) || (x, y) == (405, -1000) => Mag::zero(),
        _ => Mag::mul_up(a, b),
    }
}

#[test]
fn a_wrong_op_only_at_the_zf_bounds_or_the_1024_scale_survives_the_random_band_and_is_caught_by_the_forced_pairs(
) -> Result<(), String> {
    let ops = Ops {
        mul_up: wrong_mul_up_only_at_forced_pairs,
        ..props::real()
    };
    let cfg = FuzzConfig {
        seed: 0x5EED_F0CE,
        cases: 20_000,
        replay: None,
    };
    run_cfg(&cfg, &|rng| props::mul_case(&ops, &mut Src::random(rng)))?;
    match props::run_forced(&|src| props::mul_case(&ops, src)) {
        Ok(n) => Err(format!("the wrong mul_up survived {} forced pairs", n)),
        Err(m) if m.starts_with("forced pair ") && m.contains("op=mul_up") => Ok(()),
        Err(m) => Err(format!("unexpected forced report: {}", m)),
    }
}

#[test]
fn every_property_runs_on_every_ordered_forced_pair() -> Result<(), String> {
    let seen = core::cell::RefCell::new(Vec::new());
    let n = props::run_forced(&|src| {
        let c = props::draw_case(src, &Rel::Any).map_err(|m| Bad {
            op: "draw",
            prec: 0,
            a: String::from("-"),
            b: String::from("-"),
            msg: m,
        })?;
        seen.borrow_mut().push((c.ea, c.eb));
        Ok(0)
    })?;
    expect_eq("forced pairs", &n, &props::FORCED_PAIRS)?;
    let es = crate::common::generators::forced_prop_exponents();
    let want: Vec<(i64, i64)> = es
        .iter()
        .flat_map(|a| es.iter().map(move |b| (*a, *b)))
        .collect();
    expect_eq("forced pairs in order", &seen.into_inner(), &want)
}
