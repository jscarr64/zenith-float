//! COVERS: Mag::add_up, Mag::div_lower, Mag::div_up, Mag::mul_lower, Mag::mul_up, Mag::sqrt_lower, Mag::sqrt_up, Mag::sub_lower
//! L0a-3 (TEST_PLAN section 21), the W1 property sets: EXPONENT TRANSLATION. A value `0.m * 2^e` shifted by `k` is `0.m * 2^(e+k)`, so,
//! as long as no range limit is involved, rounding commutes with the shift and the result mantissa is bit-for-bit the same:
//!   mul_up, mul_lower:  op(a<<ka, b<<kb) = op(a, b) << (ka + kb)     div_up, div_lower: op(a<<ka, b<<kb) = op(a, b) << (ka - kb)
//!   add_up, sub_lower:  op(a<<k, b<<k)   = op(a, b) << k              sqrt_up, sqrt_lower: op(a<<k) = op(a) << (k / 2), for EVEN k only
//! (an odd k changes the parity of the exponent, which changes the mantissa: CONTRACT 14.4 fixes the parity first). The sqrt sets are run
//! twice, once with even base exponents and once with odd ones (the separate odd-parity set).
//! Oracle: the relation itself (a metamorphic oracle, no model needed); the value model is covered by tests/mag. Scope: finite nonzero operands
//! near exponent 0 and shifts over the whole i64 range; zero and Infinity operands are corner cases (CONTRACT 14.6) and belong to tests/mag/corners.rs.
//! RANGE BANDS (documented predicate): a case is EXCLUDED when a translated operand exponent or the expected translated result exponent is not
//! in `i64::MIN..=i64::MAX` (checked i128 arithmetic, never wrapping). That is exactly the set where CONTRACT 3.2 to 3.4 (underflow to
//! min_pos or zero, overflow to Infinity or the largest Mag, a carry past `i64::MAX`) replaces the translated value. Inside the predicate the
//! relation is exact, including results at `i64::MIN` and `i64::MAX`. The excluded bands are enumerated by the edge set: every shift in
//! `i64::MIN..i64::MIN+3` and `i64::MAX-3..i64::MAX` (and cancelling pairs) against small base exponents, with both outcomes counted.
//! Non-vacuity: for finite nonzero operands near 0 the base result must be finite and nonzero (`sub_lower` zero iff a <= b), otherwise the case
//! fails, so a constant-zero or constant-Infinity operation cannot pass by being translation invariant.
//! Draws (order frozen, GENERATOR_VERSION governs): prec (PS-fast), a.exp (near domain), a.mant, b.exp (anchored at a.exp: controlled gap), b.mant,
//! ka and kb (full i64 domain: edge band, near 0, uniform, powers of two, gaps), selector below(3), delta uniform in -70..=70.
//! Failures return `Err` whose text starts with one S2 replay line (`suite=l0a_mag seed=<hex> case=<n> op=<op>[ka=..,kb=..] prec=.. rm=None a=.. b=..`).
//! Cases: 20,000 per operation at the fixed seed of suite `l0a_mag` (gate tier); ZB_FUZZ_CASES widens it, ZB_FUZZ_SEED and ZB_REPLAY work as in common/replay.rs.
//! The checker's own self-tests use operations that translate exactly by construction (and deliberately wrong rules and operations): they do not depend on `Mag`.

use crate::common::expect_eq;
use crate::common::generators::{
    draw_exponent, draw_mantissa, draw_precision, uniform_in, ExpDomain, ModeOrNone, PrecSource,
    TOP,
};
use crate::common::replay::{config_from_env, FuzzConfig, ReplayLine};
use crate::common::rng::{splitmix64, SplitMix64};
use crate::common::seeds::{seed_for, CASE_CLASS_GATE};
use zenith_bound::mag::Mag;

/// Cases per operation at the gate tier.
const GATE_CASES: u64 = 20_000;

/// The base exponents are drawn in this window around 0 (raw range one tick wider).
const NEAR: ExpDomain = ExpDomain {
    min: -1000,
    max: 1000,
    subnormal: None,
};

/// A finite nonzero `Mag` operand: `0.mant * 2^exp`, top bit of `mant` set.
#[derive(Debug, PartialEq, Eq)]
struct Val {
    exp: i64,
    mant: u64,
}

/// The result of an operation, as plain data.
#[derive(Debug, PartialEq, Eq)]
enum Out {
    Zero,
    Infinity,
    Finite { exp: i64, mant: u64 },
}

/// The shifts of a case: `ka` for the first operand, `kb` for the second (the rule decides which of them are used).
struct Shifts {
    ka: i64,
    kb: i64,
}

/// How a rule translates: which shifts are applied to the operands and which shift the result must show.
struct Rule {
    name: &'static str,
    unary: bool,
    uses_kb: bool,
    even_shifts: bool,
    applied: fn(&Shifts) -> (i128, i128),
    result: fn(&Shifts) -> i128,
    /// The `kb` that makes the result shift small (a huge `ka` cancelled by `kb`), for `delta` in -70..=70.
    cancel: fn(i64, i64) -> i128,
    /// The class check of the base result (non-vacuity).
    base: fn(&Val, &Val, &Out) -> Result<(), String>,
}

fn base_finite(_a: &Val, _b: &Val, out: &Out) -> Result<(), String> {
    match out {
        Out::Finite { .. } => Ok(()),
        other @ (Out::Zero | Out::Infinity) => Err(format!(
            "finite nonzero operands near exponent 0 must give a finite nonzero result, the operation gave {:?}",
            other
        )),
    }
}

fn base_sub(a: &Val, b: &Val, out: &Out) -> Result<(), String> {
    let a_le_b = (a.exp, a.mant) <= (b.exp, b.mant);
    match (a_le_b, out) {
        (true, Out::Zero) | (false, Out::Finite { .. }) => Ok(()),
        (true, other) => Err(format!("sub_lower with a <= b must be the canonical zero (CONTRACT 14.2), the operation gave {:?}", other)),
        (false, other) => Err(format!("sub_lower with a > b must be finite and nonzero near exponent 0, the operation gave {:?}", other)),
    }
}

fn applied_both(s: &Shifts) -> (i128, i128) {
    (i128::from(s.ka), i128::from(s.kb))
}
fn applied_same(s: &Shifts) -> (i128, i128) {
    (i128::from(s.ka), i128::from(s.ka))
}
fn applied_first(s: &Shifts) -> (i128, i128) {
    (i128::from(s.ka), 0)
}
fn result_sum(s: &Shifts) -> i128 {
    i128::from(s.ka) + i128::from(s.kb)
}
fn result_diff(s: &Shifts) -> i128 {
    i128::from(s.ka) - i128::from(s.kb)
}
fn result_first(s: &Shifts) -> i128 {
    i128::from(s.ka)
}
fn result_half(s: &Shifts) -> i128 {
    i128::from(s.ka).div_euclid(2)
}
fn cancel_mul(ka: i64, delta: i64) -> i128 {
    -i128::from(ka) + i128::from(delta)
}
fn cancel_div(ka: i64, delta: i64) -> i128 {
    i128::from(ka) + i128::from(delta)
}
fn cancel_none(_ka: i64, _delta: i64) -> i128 {
    0
}

const MUL: Rule = Rule {
    name: "mul",
    unary: false,
    uses_kb: true,
    even_shifts: false,
    applied: applied_both,
    result: result_sum,
    cancel: cancel_mul,
    base: base_finite,
};
const DIV: Rule = Rule {
    name: "div",
    unary: false,
    uses_kb: true,
    even_shifts: false,
    applied: applied_both,
    result: result_diff,
    cancel: cancel_div,
    base: base_finite,
};
const ADD: Rule = Rule {
    name: "add",
    unary: false,
    uses_kb: false,
    even_shifts: false,
    applied: applied_same,
    result: result_first,
    cancel: cancel_none,
    base: base_finite,
};
const SUB: Rule = Rule {
    name: "sub",
    unary: false,
    uses_kb: false,
    even_shifts: false,
    applied: applied_same,
    result: result_first,
    cancel: cancel_none,
    base: base_sub,
};
const SQRT: Rule = Rule {
    name: "sqrt",
    unary: true,
    uses_kb: false,
    even_shifts: true,
    applied: applied_first,
    result: result_half,
    cancel: cancel_none,
    base: base_finite,
};

// Deliberately wrong rules (negative self-tests): each must be rejected by the checker.
const MUL_RESULT_OFF_BY_ONE: Rule = Rule {
    result: |s| result_sum(s) + 1,
    ..MUL
};
const MUL_FIRST_OPERAND_ONLY: Rule = Rule {
    applied: applied_first,
    ..MUL
};
const DIV_SUM_INSTEAD_OF_DIFFERENCE: Rule = Rule {
    result: result_sum,
    ..DIV
};
const ADD_FIRST_OPERAND_ONLY: Rule = Rule {
    applied: applied_first,
    ..ADD
};
const SQRT_ANY_PARITY: Rule = Rule {
    even_shifts: false,
    ..SQRT
};
const SQRT_FULL_SHIFT: Rule = Rule {
    result: result_first,
    ..SQRT
};

type FakeOp = fn(&Val, &Val) -> Result<Out, String>;
type FakeCase<'a> = (&'a str, u64, &'a Rule, &'a Parity, FakeOp);
type Op<'a> = dyn Fn(&Val, &Val) -> Result<Out, String> + 'a;

/// The parity forced on the first exponent of the base operands.
enum Parity {
    Any,
    Even,
    Odd,
}

struct Case {
    prec: u32,
    a: Val,
    b: Val,
    shifts: Shifts,
}

/// Draws one case; always consumes the same draws (see the header).
fn draw_case(rng: &mut SplitMix64, rule: &Rule, parity: &Parity) -> Result<Case, String> {
    let prec = draw_precision(rng, &PrecSource::Fast)?;
    let ea = draw_exponent(rng, &NEAR, 0, prec)?.value;
    let ma = draw_mantissa(rng, prec)?.value;
    let eb = draw_exponent(rng, &NEAR, ea, prec)?.value;
    let mb = draw_mantissa(rng, prec)?.value;
    let k1 = draw_exponent(rng, &ExpDomain::MAG, 0, prec)?.value;
    let k2 = draw_exponent(rng, &ExpDomain::MAG, 0, prec)?.value;
    let selector = rng.below(3);
    let delta = uniform_in(rng, -70, 70)?;
    let ka = if rule.even_shifts {
        k1 - k1.rem_euclid(2)
    } else {
        k1
    };
    let kb = match selector {
        0 => k2,
        1 => i64::try_from((rule.cancel)(k1, delta)).unwrap_or(k2),
        _ => 0,
    };
    let ea = match (parity, ea.rem_euclid(2)) {
        (Parity::Even, 1) | (Parity::Odd, 0) => ea
            .checked_add(1)
            .ok_or_else(|| String::from("exponent overflow"))?,
        _ => ea,
    };
    Ok(Case {
        prec,
        a: Val { exp: ea, mant: ma },
        b: Val { exp: eb, mant: mb },
        shifts: Shifts { ka, kb },
    })
}

fn translate(v: &Val, k: i128) -> Option<Val> {
    let e = i128::from(v.exp).checked_add(k)?;
    Some(Val {
        exp: i64::try_from(e).ok()?,
        mant: v.mant,
    })
}

fn translate_out(o: &Out, k: i128) -> Option<Out> {
    match o {
        Out::Zero => Some(Out::Zero),
        Out::Infinity => Some(Out::Infinity),
        Out::Finite { exp, mant } => {
            let e = i128::from(*exp).checked_add(k)?;
            Some(Out::Finite {
                exp: i64::try_from(e).ok()?,
                mant: *mant,
            })
        }
    }
}

/// What happened to one case.
enum Verdict {
    Checked { near_min: bool, near_max: bool },
    OperandOut,
    ResultOut,
}

#[derive(Debug, PartialEq, Eq, Default)]
struct Tally {
    checked: u64,
    operand_out: u64,
    result_out: u64,
    near_min: u64,
    near_max: u64,
}

impl Tally {
    fn add(&mut self, v: &Verdict) {
        match v {
            Verdict::Checked { near_min, near_max } => {
                self.checked += 1;
                self.near_min += u64::from(*near_min);
                self.near_max += u64::from(*near_max);
            }
            Verdict::OperandOut => self.operand_out += 1,
            Verdict::ResultOut => self.result_out += 1,
        }
    }
}

fn near(e: i64) -> (bool, bool) {
    (e <= i64::MIN + 4, e >= i64::MAX - 4)
}

/// Checks the translation relation for one case.
fn check_case(rule: &Rule, op: &Op, a: &Val, b: &Val, s: &Shifts) -> Result<Verdict, String> {
    let base = op(a, b)?;
    (rule.base)(a, b, &base)?;
    let (da, db) = (rule.applied)(s);
    let Some(ta) = translate(a, da) else {
        return Ok(Verdict::OperandOut);
    };
    let Some(tb) = translate(b, db) else {
        return Ok(Verdict::OperandOut);
    };
    let Some(want) = translate_out(&base, (rule.result)(s)) else {
        return Ok(Verdict::ResultOut);
    };
    let got = op(&ta, &tb)?;
    if got != want {
        return Err(format!(
            "translated result {:?} differs from the base result {:?} shifted to {:?} (operands moved to {:?} and {:?})",
            got, base, want, ta, tb
        ));
    }
    let mut flags = vec![near(ta.exp), near(tb.exp)];
    if let Out::Finite { exp, .. } = &want {
        flags.push(near(*exp));
    }
    Ok(Verdict::Checked {
        near_min: flags.iter().any(|f| f.0),
        near_max: flags.iter().any(|f| f.1),
    })
}

fn enc(v: &Val) -> String {
    format!("{:016x}e{}", v.mant, v.exp)
}

fn replay_text(label: &str, seed: u64, case: u64, prec: u32, rule: &Rule, c: &Case) -> String {
    let line = ReplayLine {
        suite: String::from("l0a_mag"),
        seed,
        case,
        op: format!("{}[ka={},kb={}]", label, c.shifts.ka, c.shifts.kb),
        prec: usize::try_from(prec).unwrap_or(usize::MAX),
        rm: String::from(ModeOrNone::None.name()),
        a: enc(&c.a),
        b: if rule.unary {
            String::from("-")
        } else {
            enc(&c.b)
        },
    };
    match line.format_line() {
        Ok(t) => t,
        Err(e) => format!("(replay line not formattable: {})", e),
    }
}

/// The base seed of one operation's stream: the `l0a_mag` gate seed of seeds.rs, mixed with the operation index.
fn base_seed(op_index: u64) -> Result<u64, String> {
    let s = seed_for("l0a_mag", CASE_CLASS_GATE)
        .ok_or_else(|| String::from("suite l0a_mag is not in the seeds table"))?;
    Ok(splitmix64(crate::common::modular::add_u64(s, op_index)))
}

/// Runs the random set. `honor_env`: ZB_FUZZ_* / ZB_PROP_* / ZB_REPLAY apply (the real-Mag tests); the self-tests pass `false`.
fn run_random(
    label: &str,
    op_index: u64,
    rule: &Rule,
    parity: &Parity,
    op: &Op,
    honor_env: bool,
    cases: u64,
) -> Result<Tally, String> {
    let base = base_seed(op_index)?;
    let cfg = if honor_env {
        config_from_env(base, cases)?
    } else {
        FuzzConfig {
            seed: base,
            cases,
            replay: None,
        }
    };
    let seed = match &cfg.replay {
        Some(t) => t.seed,
        None => cfg.seed,
    };
    let mut tally = Tally::default();
    for case in cfg.case_range()? {
        let mut rng = SplitMix64::new(cfg.seed_of(case));
        let c = draw_case(&mut rng, rule, parity)?;
        match check_case(rule, op, &c.a, &c.b, &c.shifts) {
            Ok(v) => tally.add(&v),
            Err(e) => {
                return Err(format!(
                    "{} -- {}",
                    replay_text(label, seed, case, c.prec, rule, &c),
                    e
                ))
            }
        }
    }
    Ok(tally)
}

const EDGE_EXPS: [i64; 5] = [-70, -1, 0, 1, 70];
const EDGE_MANTS: [u64; 4] = [TOP, TOP + 1, u64::MAX, 0xAAAA_AAAA_AAAA_AAAB];

/// The shifts of the edge set: i64::MIN..MIN+3 and MAX-3..MAX.
fn edge_shifts() -> Vec<i64> {
    let mut v: Vec<i64> = Vec::new();
    for d in 0..4 {
        v.push(i64::MIN + d);
        v.push(i64::MAX - d);
    }
    v
}

/// The (ka, kb) pairs of the edge set for `rule`.
fn edge_pairs(rule: &Rule) -> Vec<Shifts> {
    let mut out: Vec<Shifts> = Vec::new();
    for ka in edge_shifts() {
        if rule.even_shifts && ka.rem_euclid(2) != 0 {
            continue;
        }
        if !rule.uses_kb {
            out.push(Shifts { ka, kb: 0 });
            continue;
        }
        for kb in edge_shifts() {
            out.push(Shifts { ka, kb });
        }
        out.push(Shifts { ka, kb: 0 });
        out.push(Shifts { ka: 0, kb: ka });
        if let Some(neg) = ka.checked_neg() {
            out.push(Shifts { ka, kb: neg });
        }
    }
    out
}

/// The edge set: small base exponents and mantissa patterns against the shifts at the ends of the i64 range.
fn run_edge(label: &str, rule: &Rule, parity: &Parity, op: &Op) -> Result<Tally, String> {
    let mut tally = Tally::default();
    let mut n = 0u64;
    for ea in EDGE_EXPS {
        if matches!(
            (parity, ea.rem_euclid(2)),
            (Parity::Even, 1) | (Parity::Odd, 0)
        ) {
            continue;
        }
        for ma in EDGE_MANTS {
            for eb in EDGE_EXPS {
                for mb in EDGE_MANTS {
                    if rule.unary && (eb != 0 || mb != TOP) {
                        continue;
                    }
                    let a = Val { exp: ea, mant: ma };
                    let b = Val { exp: eb, mant: mb };
                    for s in edge_pairs(rule) {
                        match check_case(rule, op, &a, &b, &s) {
                            Ok(v) => tally.add(&v),
                            Err(e) => {
                                let c = Case {
                                    prec: 64,
                                    a,
                                    b,
                                    shifts: s,
                                };
                                let text =
                                    replay_text(&format!("{}-edge", label), 0, n, 64, rule, &c);
                                return Err(format!("{} -- {}", text, e));
                            }
                        }
                        n += 1;
                    }
                }
            }
        }
    }
    Ok(tally)
}

// ---- the real Mag ------------------------------------------------------------------------------------------------------------

fn mag_from(mant: u64, exp: i64) -> Result<Mag, String> {
    Mag::finite(mant, exp).map_err(|e| format!("Mag::finite({:#x}, {}) failed: {:?}", mant, exp, e))
}

fn out_of(m: &Mag) -> Result<Out, String> {
    match (m.exp(), m.mant()) {
        (Some(_), Some(0)) => Ok(Out::Zero),
        (Some(exp), Some(mant)) => Ok(Out::Finite { exp, mant }),
        (None, None) => Ok(Out::Infinity),
        (e, m) => Err(format!(
            "inconsistent accessors: exp() {:?}, mant() {:?}",
            e, m
        )),
    }
}

fn binary(f: fn(&Mag, &Mag) -> Mag) -> impl Fn(&Val, &Val) -> Result<Out, String> {
    move |a: &Val, b: &Val| out_of(&f(&mag_from(a.mant, a.exp)?, &mag_from(b.mant, b.exp)?))
}

fn unary(f: fn(&Mag) -> Mag) -> impl Fn(&Val, &Val) -> Result<Out, String> {
    move |a: &Val, _b: &Val| out_of(&f(&mag_from(a.mant, a.exp)?))
}

fn w1_random(
    label: &str,
    op_index: u64,
    rule: &Rule,
    parity: &Parity,
    op: &Op,
) -> Result<(), String> {
    let t = run_random(label, op_index, rule, parity, op, true, GATE_CASES)?;
    if t.checked == 0 {
        return Err(format!(
            "{}: no case was inside the predicate ({:?})",
            label, t
        ));
    }
    Ok(())
}

fn w1_edge(label: &str, rule: &Rule, parity: &Parity, op: &Op) -> Result<(), String> {
    let t = run_edge(label, rule, parity, op)?;
    if t.near_min == 0 || t.near_max == 0 || t.operand_out + t.result_out == 0 {
        return Err(format!(
            "{}: the edge set did not reach both ends of the range and the excluded bands ({:?})",
            label, t
        ));
    }
    Ok(())
}

#[test]
fn mag_adapter_round_trips_finite_values() -> Result<(), String> {
    for (exp, mant) in [
        (0i64, TOP),
        (5, u64::MAX),
        (-7, 0xAAAA_AAAA_AAAA_AAAB),
        (i64::MIN, TOP),
        (i64::MAX, u64::MAX),
    ] {
        let v = Val { exp, mant };
        let m = mag_from(v.mant, v.exp)?;
        expect_eq(
            "Mag::finite then exp and mant",
            &out_of(&m)?,
            &Out::Finite { exp, mant },
        )?;
    }
    expect_eq("zero", &out_of(&Mag::zero())?, &Out::Zero)?;
    expect_eq("infinity", &out_of(&Mag::infinity())?, &Out::Infinity)
}

#[test]
fn w1_add_up_translates_exactly() -> Result<(), String> {
    w1_random("add_up", 0, &ADD, &Parity::Any, &binary(Mag::add_up))
}
#[test]
fn w1_sub_lower_translates_exactly() -> Result<(), String> {
    w1_random("sub_lower", 1, &SUB, &Parity::Any, &binary(Mag::sub_lower))
}
#[test]
fn w1_mul_up_translates_by_the_sum() -> Result<(), String> {
    w1_random("mul_up", 2, &MUL, &Parity::Any, &binary(Mag::mul_up))
}
#[test]
fn w1_mul_lower_translates_by_the_sum() -> Result<(), String> {
    w1_random("mul_lower", 3, &MUL, &Parity::Any, &binary(Mag::mul_lower))
}
#[test]
fn w1_div_up_translates_by_the_difference() -> Result<(), String> {
    w1_random("div_up", 4, &DIV, &Parity::Any, &binary(Mag::div_up))
}
#[test]
fn w1_div_lower_translates_by_the_difference() -> Result<(), String> {
    w1_random("div_lower", 5, &DIV, &Parity::Any, &binary(Mag::div_lower))
}
#[test]
fn w1_sqrt_up_translates_by_half_the_even_shift() -> Result<(), String> {
    w1_random("sqrt_up", 6, &SQRT, &Parity::Even, &unary(Mag::sqrt_up))
}
#[test]
fn w1_sqrt_lower_translates_by_half_the_even_shift() -> Result<(), String> {
    w1_random(
        "sqrt_lower",
        7,
        &SQRT,
        &Parity::Even,
        &unary(Mag::sqrt_lower),
    )
}
#[test]
fn w1_sqrt_up_odd_exponent_parity_set() -> Result<(), String> {
    w1_random("sqrt_up-odd", 16, &SQRT, &Parity::Odd, &unary(Mag::sqrt_up))
}
#[test]
fn w1_sqrt_lower_odd_exponent_parity_set() -> Result<(), String> {
    w1_random(
        "sqrt_lower-odd",
        17,
        &SQRT,
        &Parity::Odd,
        &unary(Mag::sqrt_lower),
    )
}

#[test]
fn w1_edge_add_up() -> Result<(), String> {
    w1_edge("add_up", &ADD, &Parity::Any, &binary(Mag::add_up))
}
#[test]
fn w1_edge_sub_lower() -> Result<(), String> {
    w1_edge("sub_lower", &SUB, &Parity::Any, &binary(Mag::sub_lower))
}
#[test]
fn w1_edge_mul_up() -> Result<(), String> {
    w1_edge("mul_up", &MUL, &Parity::Any, &binary(Mag::mul_up))
}
#[test]
fn w1_edge_mul_lower() -> Result<(), String> {
    w1_edge("mul_lower", &MUL, &Parity::Any, &binary(Mag::mul_lower))
}
#[test]
fn w1_edge_div_up() -> Result<(), String> {
    w1_edge("div_up", &DIV, &Parity::Any, &binary(Mag::div_up))
}
#[test]
fn w1_edge_div_lower() -> Result<(), String> {
    w1_edge("div_lower", &DIV, &Parity::Any, &binary(Mag::div_lower))
}
#[test]
fn w1_edge_sqrt_up() -> Result<(), String> {
    // both parities of the base exponent, even shifts
    w1_edge("sqrt_up", &SQRT, &Parity::Any, &unary(Mag::sqrt_up))
}
#[test]
fn w1_edge_sqrt_lower() -> Result<(), String> {
    w1_edge("sqrt_lower", &SQRT, &Parity::Any, &unary(Mag::sqrt_lower))
}

// ---- self-tests of the checker: operations that translate exactly by construction, and deliberately wrong ones -------------------
// None of them calls Mag.

fn fin(exp: i128, mant: u64) -> Result<Out, String> {
    let e = i64::try_from(exp).map_err(|_| String::from("fake result exponent out of range"))?;
    Ok(Out::Finite {
        exp: e,
        mant: mant | TOP,
    })
}

fn gap_rot(a: &Val, b: &Val) -> Result<u32, String> {
    let gap = i128::from(b.exp) - i128::from(a.exp);
    u32::try_from(gap.rem_euclid(64)).map_err(|_| String::from("rotation"))
}

fn fake_mul(a: &Val, b: &Val) -> Result<Out, String> {
    fin(
        i128::from(a.exp) + i128::from(b.exp),
        a.mant ^ b.mant.rotate_left(7),
    )
}
fn fake_div(a: &Val, b: &Val) -> Result<Out, String> {
    fin(
        i128::from(a.exp) - i128::from(b.exp),
        a.mant.rotate_left(3) ^ b.mant,
    )
}
fn fake_add(a: &Val, b: &Val) -> Result<Out, String> {
    fin(
        i128::from(a.exp.max(b.exp)),
        a.mant ^ b.mant.rotate_left(gap_rot(a, b)?),
    )
}
fn fake_sub(a: &Val, b: &Val) -> Result<Out, String> {
    if (a.exp, a.mant) <= (b.exp, b.mant) {
        Ok(Out::Zero)
    } else {
        fin(
            i128::from(a.exp),
            a.mant ^ b.mant.rotate_left(gap_rot(a, b)?),
        )
    }
}
fn fake_sqrt(a: &Val, _b: &Val) -> Result<Out, String> {
    let odd = u32::from(a.exp.rem_euclid(2) == 1);
    fin(
        (i128::from(a.exp) + 1).div_euclid(2),
        a.mant.rotate_left(odd),
    )
}

// not translation covariant
fn bad_mul(a: &Val, b: &Val) -> Result<Out, String> {
    fin(i128::from(a.exp) + 2 * i128::from(b.exp), a.mant ^ b.mant)
}
fn bad_div(a: &Val, b: &Val) -> Result<Out, String> {
    fin(
        i128::from(a.exp) - i128::from(b.exp),
        a.mant ^ u64::try_from(a.exp.rem_euclid(64)).map_err(|_| String::from("rem"))?,
    )
}
fn bad_add(a: &Val, b: &Val) -> Result<Out, String> {
    fin(
        i128::from(a.exp.max(b.exp)),
        a.mant ^ u64::try_from(a.exp.rem_euclid(64)).map_err(|_| String::from("rem"))?,
    )
}
fn bad_sub(a: &Val, b: &Val) -> Result<Out, String> {
    if (a.exp, a.mant) <= (b.exp, b.mant) {
        Ok(Out::Zero)
    } else {
        fin(
            i128::from(a.exp),
            a.mant ^ u64::try_from(a.exp.rem_euclid(8)).map_err(|_| String::from("rem"))?,
        )
    }
}
fn bad_sqrt(a: &Val, _b: &Val) -> Result<Out, String> {
    fin(
        (i128::from(a.exp) + 1).div_euclid(2),
        a.mant ^ u64::try_from(a.exp.rem_euclid(4)).map_err(|_| String::from("rem"))?,
    )
}

fn always_zero(_a: &Val, _b: &Val) -> Result<Out, String> {
    Ok(Out::Zero)
}
fn always_infinity(_a: &Val, _b: &Val) -> Result<Out, String> {
    Ok(Out::Infinity)
}

fn rejected(what: &str, needle: &str, r: Result<Tally, String>) -> Result<(), String> {
    match r {
        // it must be the relation that failed, not a generator or adapter error
        // (a wrong rule may also drive the exact fake out of its range, which is a rejection too)
        Err(e) if e.contains(needle) || e.contains("fake result exponent out of range") => Ok(()),
        Err(e) => Err(format!("{} was rejected for the wrong reason: {}", what, e)),
        Ok(t) => Err(format!("the checker accepted {} ({:?})", what, t)),
    }
}

#[test]
fn the_checker_accepts_exactly_translating_operations_and_pins_its_counts() -> Result<(), String> {
    let cases: [FakeCase; 6] = [
        ("fake_add", 100, &ADD, &Parity::Any, fake_add),
        ("fake_sub", 101, &SUB, &Parity::Any, fake_sub),
        ("fake_mul", 102, &MUL, &Parity::Any, fake_mul),
        ("fake_div", 103, &DIV, &Parity::Any, fake_div),
        ("fake_sqrt_even", 104, &SQRT, &Parity::Even, fake_sqrt),
        ("fake_sqrt_odd", 105, &SQRT, &Parity::Odd, fake_sqrt),
    ];
    let mut got: Vec<(u64, u64, u64)> = Vec::new();
    for (label, index, rule, parity, op) in cases {
        let random = run_random(label, index, rule, parity, &op, false, GATE_CASES)?;
        if random.checked * 10 < GATE_CASES * 4 {
            return Err(format!(
                "{}: fewer than 40% of the cases are inside the predicate ({:?})",
                label, random
            ));
        }
        got.push((random.checked, random.operand_out, random.result_out));
        let edge = run_edge(label, rule, parity, &op)?;
        if edge.checked == 0
            || edge.near_min == 0
            || edge.near_max == 0
            || edge.operand_out + edge.result_out == 0
        {
            return Err(format!(
                "{}: the edge set misses a band ({:?})",
                label, edge
            ));
        }
    }
    expect_eq(
        "pinned (checked, operand out, result out) per fake operation",
        &got,
        &Vec::from(PINNED_RANDOM_TALLIES),
    )
}

// Pinned from this checker at GENERATOR_VERSION 3, unchanged at 4 (a drift detector for the draw order and the predicate).
const PINNED_RANDOM_TALLIES: [(u64, u64, u64); 6] = [
    (15712, 4288, 0),
    (15739, 4261, 0),
    (14950, 3939, 1111),
    (14932, 3920, 1148),
    (17219, 2781, 0),
    (17154, 2846, 0),
];

#[test]
fn deliberately_wrong_rules_are_rejected() -> Result<(), String> {
    rejected(
        "a result shift off by one",
        "translated result",
        run_random(
            "mul_up",
            200,
            &MUL_RESULT_OFF_BY_ONE,
            &Parity::Any,
            &fake_mul,
            false,
            2000,
        ),
    )?;
    rejected(
        "a translation of only the first operand (mul)",
        "translated result",
        run_random(
            "mul_up",
            201,
            &MUL_FIRST_OPERAND_ONLY,
            &Parity::Any,
            &fake_mul,
            false,
            2000,
        ),
    )?;
    rejected(
        "the sum as the division shift",
        "translated result",
        run_random(
            "div_up",
            202,
            &DIV_SUM_INSTEAD_OF_DIFFERENCE,
            &Parity::Any,
            &fake_div,
            false,
            2000,
        ),
    )?;
    rejected(
        "a translation of only the first operand (add)",
        "translated result",
        run_random(
            "add_up",
            203,
            &ADD_FIRST_OPERAND_ONLY,
            &Parity::Any,
            &fake_add,
            false,
            2000,
        ),
    )?;
    rejected(
        "the full shift as the sqrt shift",
        "translated result",
        run_random(
            "sqrt_up",
            204,
            &SQRT_FULL_SHIFT,
            &Parity::Even,
            &fake_sqrt,
            false,
            2000,
        ),
    )?;
    rejected(
        "odd shifts for sqrt",
        "translated result",
        run_random(
            "sqrt_up",
            205,
            &SQRT_ANY_PARITY,
            &Parity::Even,
            &fake_sqrt,
            false,
            2000,
        ),
    )?;
    // the same wrong rules are rejected by the edge set too
    rejected(
        "an off-by-one result shift at the edge",
        "translated result",
        run_edge("mul_up", &MUL_RESULT_OFF_BY_ONE, &Parity::Any, &fake_mul),
    )?;
    rejected(
        "odd sqrt shifts at the edge",
        "translated result",
        run_edge("sqrt_up", &SQRT_ANY_PARITY, &Parity::Any, &fake_sqrt),
    )
}

#[test]
fn operations_that_are_not_translation_covariant_are_rejected() -> Result<(), String> {
    rejected(
        "a multiplication that weighs b.exp twice",
        "translated result",
        run_random("mul_up", 210, &MUL, &Parity::Any, &bad_mul, false, 2000),
    )?;
    rejected(
        "a division whose mantissa depends on the exponent",
        "translated result",
        run_random("div_up", 211, &DIV, &Parity::Any, &bad_div, false, 2000),
    )?;
    rejected(
        "an addition whose mantissa depends on the exponent",
        "translated result",
        run_random("add_up", 212, &ADD, &Parity::Any, &bad_add, false, 2000),
    )?;
    rejected(
        "a subtraction whose mantissa depends on the exponent mod 8",
        "translated result",
        run_random("sub_lower", 213, &SUB, &Parity::Any, &bad_sub, false, 2000),
    )?;
    rejected(
        "a square root whose mantissa depends on the exponent mod 4",
        "translated result",
        run_random("sqrt_up", 214, &SQRT, &Parity::Even, &bad_sqrt, false, 2000),
    )
}

#[test]
fn the_sqrt_parity_condition_is_necessary_and_the_odd_set_differs_from_the_even_set(
) -> Result<(), String> {
    // the exact fake passes with even shifts for even and for odd base exponents
    run_random(
        "fake_sqrt",
        220,
        &SQRT,
        &Parity::Even,
        &fake_sqrt,
        false,
        3000,
    )?;
    run_random(
        "fake_sqrt",
        221,
        &SQRT,
        &Parity::Odd,
        &fake_sqrt,
        false,
        3000,
    )?;
    // the two sets draw different operands from the same seed (the parity is forced on the first exponent)
    let mut r1 = SplitMix64::new(5);
    let mut r2 = SplitMix64::new(5);
    let even = draw_case(&mut r1, &SQRT, &Parity::Even)?;
    let odd = draw_case(&mut r2, &SQRT, &Parity::Odd)?;
    expect_eq("even set", &even.a.exp.rem_euclid(2), &0)?;
    expect_eq("odd set", &odd.a.exp.rem_euclid(2), &1)?;
    expect_eq(
        "shifts are even for sqrt",
        &(even.shifts.ka.rem_euclid(2), odd.shifts.ka.rem_euclid(2)),
        &(0, 0),
    )?;
    // an operation whose mantissa rule depends on the parity (the fake) is not invariant under an odd shift
    let a = Val {
        exp: 4,
        mant: TOP + 12345,
    };
    let moved = Val {
        exp: 5,
        mant: TOP + 12345,
    };
    let base = fake_sqrt(&a, &a)?;
    let odd_moved = fake_sqrt(&moved, &moved)?;
    if let (Out::Finite { mant: m0, .. }, Out::Finite { mant: m1, .. }) = (&base, &odd_moved) {
        if m0 == m1 {
            return Err(String::from(
                "the parity-sensitive fake did not change its mantissa under an odd shift",
            ));
        }
    }
    Ok(())
}

#[test]
fn the_range_predicate_excludes_exactly_the_out_of_range_bands() -> Result<(), String> {
    // translate: operand exponents (checked, in i128)
    let v = Val { exp: 3, mant: TOP };
    expect_eq(
        "to MAX",
        &translate(&v, i128::from(i64::MAX) - 3).map(|x| x.exp),
        &Some(i64::MAX),
    )?;
    expect_eq("past MAX", &translate(&v, i128::from(i64::MAX) - 2), &None)?;
    expect_eq(
        "to MIN",
        &translate(&v, i128::from(i64::MIN) - 3).map(|x| x.exp),
        &Some(i64::MIN),
    )?;
    expect_eq("past MIN", &translate(&v, i128::from(i64::MIN) - 4), &None)?;
    expect_eq("huge shift", &translate(&v, i128::MAX), &None)?;
    // translate_out: results; zero and Infinity never move and never leave the range
    let r = Out::Finite { exp: -5, mant: TOP };
    expect_eq(
        "result to MIN",
        &translate_out(&r, i128::from(i64::MIN) + 5),
        &Some(Out::Finite {
            exp: i64::MIN,
            mant: TOP,
        }),
    )?;
    expect_eq(
        "result below MIN",
        &translate_out(&r, i128::from(i64::MIN) + 4),
        &None,
    )?;
    expect_eq(
        "result to MAX",
        &translate_out(&Out::Finite { exp: 7, mant: TOP }, i128::from(i64::MAX) - 7),
        &Some(Out::Finite {
            exp: i64::MAX,
            mant: TOP,
        }),
    )?;
    expect_eq(
        "result above MAX",
        &translate_out(&Out::Finite { exp: 7, mant: TOP }, i128::from(i64::MAX) - 6),
        &None,
    )?;
    expect_eq(
        "zero stays",
        &translate_out(&Out::Zero, i128::MAX),
        &Some(Out::Zero),
    )?;
    expect_eq(
        "infinity stays",
        &translate_out(&Out::Infinity, i128::MIN),
        &Some(Out::Infinity),
    )?;
    // check_case classifies: operand out, result out, checked at the very ends
    let one = Val { exp: 0, mant: TOP };
    let s = Shifts {
        ka: i64::MAX,
        kb: 0,
    };
    let kind = |r: Result<Verdict, String>| -> Result<&'static str, String> {
        Ok(match r? {
            Verdict::Checked { .. } => "checked",
            Verdict::OperandOut => "operand",
            Verdict::ResultOut => "result",
        })
    };
    // fake_mul: result exp = 0 + 0 at base; shifting a by MAX gives result MAX (in), b by 1 more pushes the result out
    expect_eq(
        "result exactly MAX",
        &kind(check_case(&MUL, &fake_mul, &one, &one, &s))?,
        &"checked",
    )?;
    expect_eq(
        "result past MAX",
        &kind(check_case(
            &MUL,
            &fake_mul,
            &one,
            &one,
            &Shifts {
                ka: i64::MAX,
                kb: 1,
            },
        ))?,
        &"result",
    )?;
    expect_eq(
        "operand past MAX",
        &kind(check_case(
            &MUL,
            &fake_mul,
            &Val { exp: 1, mant: TOP },
            &one,
            &s,
        ))?,
        &"operand",
    )?;
    expect_eq(
        "result exactly MIN",
        &kind(check_case(
            &MUL,
            &fake_mul,
            &one,
            &one,
            &Shifts {
                ka: i64::MIN,
                kb: 0,
            },
        ))?,
        &"checked",
    )?;
    expect_eq(
        "result below MIN",
        &kind(check_case(
            &MUL,
            &fake_mul,
            &one,
            &one,
            &Shifts {
                ka: i64::MIN,
                kb: -1,
            },
        ))?,
        &"result",
    )?;
    // a wrapping operation passes inside the predicate and is not asked about the excluded bands
    let wrapping = |a: &Val, b: &Val| -> Result<Out, String> {
        Ok(Out::Finite {
            exp: crate::common::modular::add_i64(a.exp, b.exp),
            mant: (a.mant ^ b.mant) | TOP,
        })
    };
    run_random(
        "modular_add_op",
        230,
        &MUL,
        &Parity::Any,
        &wrapping,
        false,
        3000,
    )?;
    Ok(())
}

#[test]
fn a_constant_zero_or_infinity_operation_cannot_pass_by_being_translation_invariant(
) -> Result<(), String> {
    for rule in [&ADD, &MUL, &DIV, &SQRT] {
        rejected(
            "a constant zero operation",
            "must",
            run_random(rule.name, 240, rule, &Parity::Any, &always_zero, false, 200),
        )?;
        rejected(
            "a constant Infinity operation",
            "must",
            run_random(
                rule.name,
                241,
                rule,
                &Parity::Any,
                &always_infinity,
                false,
                200,
            ),
        )?;
        rejected(
            "a constant zero operation at the edge",
            "must",
            run_edge(rule.name, rule, &Parity::Any, &always_zero),
        )?;
    }
    // sub_lower: zero is right exactly when a <= b, so constant zero and constant Infinity are both caught
    rejected(
        "constant zero for sub_lower",
        "must",
        run_random(
            "sub_lower",
            242,
            &SUB,
            &Parity::Any,
            &always_zero,
            false,
            200,
        ),
    )?;
    rejected(
        "constant Infinity for sub_lower",
        "must",
        run_random(
            "sub_lower",
            243,
            &SUB,
            &Parity::Any,
            &always_infinity,
            false,
            200,
        ),
    )
}

#[test]
fn a_failure_carries_a_replay_line_that_parses_back() -> Result<(), String> {
    let err = match run_random(
        "mul_up",
        250,
        &MUL_RESULT_OFF_BY_ONE,
        &Parity::Any,
        &fake_mul,
        false,
        2000,
    ) {
        Err(e) => e,
        Ok(t) => return Err(format!("the wrong rule passed ({:?})", t)),
    };
    let (line, reason) = err
        .split_once(" -- ")
        .ok_or_else(|| format!("no replay line in {:?}", err))?;
    let parsed = ReplayLine::parse_line(line)?;
    expect_eq("suite", &parsed.suite.as_str(), &"l0a_mag")?;
    expect_eq("rm", &parsed.rm.as_str(), &"None")?;
    expect_eq("seed", &parsed.seed, &base_seed(250)?)?;
    if !parsed.op.starts_with("mul_up[ka=") || reason.is_empty() {
        return Err(format!("unexpected failure text: {}", err));
    }
    // the line reproduces the failure: redraw the case from (seed, case) and check it again
    let cfg = FuzzConfig {
        seed: parsed.seed,
        cases: 1,
        replay: None,
    };
    let mut rng = SplitMix64::new(cfg.seed_of(parsed.case));
    let c = draw_case(&mut rng, &MUL_RESULT_OFF_BY_ONE, &Parity::Any)?;
    expect_eq("a", &enc(&c.a), &parsed.a)?;
    expect_eq("b", &enc(&c.b), &parsed.b)?;
    if check_case(&MUL_RESULT_OFF_BY_ONE, &fake_mul, &c.a, &c.b, &c.shifts).is_ok() {
        return Err(String::from("the replayed case does not fail"));
    }
    Ok(())
}

#[test]
fn the_gate_tier_runs_twenty_thousand_cases_per_operation_from_distinct_seeds() -> Result<(), String>
{
    expect_eq("cases", &GATE_CASES, &20_000)?;
    let mut seeds: Vec<u64> = Vec::new();
    for i in [0u64, 1, 2, 3, 4, 5, 6, 7, 16, 17] {
        seeds.push(base_seed(i)?);
    }
    seeds.sort_unstable();
    seeds.dedup();
    expect_eq("distinct per-operation seeds", &seeds.len(), &10)?;
    // the draw is deterministic and consumes the same draws for every rule
    let mut r1 = SplitMix64::new(9);
    let mut r2 = SplitMix64::new(9);
    let c1 = draw_case(&mut r1, &MUL, &Parity::Any)?;
    let _ = draw_case(&mut r2, &SQRT, &Parity::Any)?;
    let after1 = r1.next_u64();
    let after2 = r2.next_u64();
    expect_eq("same consumption for every rule", &after1, &after2)?;
    let mut r3 = SplitMix64::new(9);
    let c3 = draw_case(&mut r3, &MUL, &Parity::Any)?;
    expect_eq(
        "same case",
        &(c1.a, c1.b, c1.shifts.ka, c1.shifts.kb),
        &(c3.a, c3.b, c3.shifts.ka, c3.shifts.kb),
    )
}
