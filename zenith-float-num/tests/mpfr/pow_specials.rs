//! W6: `pow` with an infinite operand, and `powi`/`powsi` of a signed zero, against MPFR (`mpfr_pow`, `mpfr_pow_ui`,
//! `mpfr_pow_si`, which follow IEEE 754 / C99 Annex F). Results zenith-float still gets wrong are listed in
//! `KNOWN_WRONG`: the test fails if one of them starts to agree (remove it from the list with its fix) or if any
//! other case disagrees.

use crate::mpfr::common::{assert_exact_eq, mpfr_reference_with};
use crate::mpfr::retry_budget::{round_of, ternary};
use rug::float::Special;
use rug::ops::{AssignRound, Pow};
use rug::Float;
use std::collections::BTreeSet;
use std::ops::Neg;
use std::panic::{catch_unwind, AssertUnwindSafe};
use zenith_float_num::{Consts, ExactNum, RoundingMode, INF_NEG, INF_POS, NAN};

type Res = Result<(), Box<dyn std::error::Error>>;

const P: usize = 64;
const RM: RoundingMode = RoundingMode::ToEven;

/// Results that differ from MPFR today; each fix removes its entries.
const KNOWN_WRONG: &[&str] = &[
    // (-Inf)^(negative odd integer) is -0
    "pow(-Inf, -1)",
    "pow(-Inf, -3)",
    // powi(-0, even) is +0
    "powi(-0, 2)",
    "powi(-0, 4)",
    // powsi(-0, n): -Inf for negative odd n, +0 for positive even n
    "powsi(-0, -1)",
    "powsi(-0, -3)",
    "powsi(-0, 2)",
    "powsi(-0, 4)",
];

/// A base or exponent: `(label, zenith value, MPFR value)`.
fn val(label: &str, sig: u64, exp2: i32, neg: bool) -> (String, ExactNum, Float) {
    let mut z = ExactNum::from_word(sig, P).ldexp(exp2 as _, P, RM);
    let mut f = Float::with_val(P as u32, sig) << exp2;
    if neg {
        z = z.neg();
        f = -f;
    }
    (label.to_string(), z, f)
}

fn inf(neg: bool) -> (String, ExactNum, Float) {
    let f = Float::with_val(
        P as u32,
        if neg { Special::NegInfinity } else { Special::Infinity },
    );
    (
        if neg { "-Inf" } else { "+Inf" }.to_string(),
        if neg { INF_NEG } else { INF_POS },
        f,
    )
}

struct Cases {
    cc: Consts,
    wrong: BTreeSet<String>,
    checked: usize,
}

impl Cases {
    fn check(
        &mut self,
        label: String,
        got: &ExactNum,
        f: impl Fn(&mut Float, rug::float::Round) -> std::cmp::Ordering,
    ) {
        self.checked += 1;
        let cc = &mut self.cc;
        let ok = catch_unwind(AssertUnwindSafe(|| -> Res {
            let (want, _) =
                mpfr_reference_with(P as u32, RM, |r, rnd| ternary(f(r, round_of(rnd))))?;
            assert_exact_eq(got, &want, P, &label, cc)
        }));
        if !matches!(ok, Ok(Ok(()))) {
            self.wrong.insert(label);
        }
    }
}

#[test]
fn pow_inf_and_signed_zero_powers_match_mpfr() -> Res {
    let mut c = Cases {
        cc: Consts::new()?,
        wrong: BTreeSet::new(),
        checked: 0,
    };
    let finite_bases = [
        val("0.5", 1, -1, false),
        val("1", 1, 0, false),
        val("2", 1, 1, false),
        val("-0.5", 1, -1, true),
        val("-1", 1, 0, true),
        val("-2", 1, 1, true),
        val("-3", 3, 0, true),
        val("+0", 0, 0, false),
        val("-0", 0, 0, true),
    ];
    // x^±Inf
    for (bl, b, fb) in &finite_bases {
        for (el, e, fe) in [inf(false), inf(true)] {
            let mut cc = Consts::new()?;
            let got = b.pow(&e, P, RM, &mut cc);
            c.check(format!("pow({bl}, {el})"), &got, |r, rd| {
                r.assign_round(fb.pow(&fe), rd)
            });
        }
    }
    // ±Inf^y, y odd / even / non-integer, both signs
    let exps = [
        val("-1", 1, 0, true),
        val("-2", 1, 1, true),
        val("-3", 3, 0, true),
        val("-0.5", 1, -1, true),
        val("1", 1, 0, false),
        val("2", 1, 1, false),
        val("3", 3, 0, false),
        val("0.5", 1, -1, false),
    ];
    for (bl, b, fb) in [inf(false), inf(true)] {
        for (el, e, fe) in &exps {
            let mut cc = Consts::new()?;
            let got = b.pow(e, P, RM, &mut cc);
            c.check(format!("pow({bl}, {el})"), &got, |r, rd| {
                r.assign_round((&fb).pow(fe), rd)
            });
        }
    }
    // powi / powsi of a signed zero
    for (bl, b, fb) in [val("+0", 0, 0, false), val("-0", 0, 0, true)] {
        for n in 0..=4usize {
            let got = b.powi(n, P, RM);
            c.check(format!("powi({bl}, {n})"), &got, |r, rd| {
                r.assign_round((&fb).pow(n as u32), rd)
            });
        }
        for n in -4i64..=4 {
            let got = b.powsi(n, P, RM);
            c.check(format!("powsi({bl}, {n})"), &got, |r, rd| {
                r.assign_round((&fb).pow(n as i32), rd)
            });
        }
    }
    assert!(c.checked >= 62, "too few cases: {}", c.checked);
    let known: BTreeSet<String> = KNOWN_WRONG.iter().map(|s| s.to_string()).collect();
    assert_eq!(
        c.wrong, known,
        "cases that differ from MPFR (left) vs KNOWN_WRONG (right)"
    );
    Ok(())
}

/// zenith-float's documented NaN-operand rule (ext.rs): an operation with a NaN operand returns NaN, whatever
/// IEEE 754 says for `pow(1, NaN)` and `pow(NaN, 0)` (both 1 there). Not a bug; this pins the rule.
#[test]
fn pow_with_a_nan_operand_is_nan() -> Res {
    let mut cc = Consts::new()?;
    let one = ExactNum::from_word(1, P);
    let zero = ExactNum::from_word(0, P);
    assert!(one.pow(&NAN, P, RM, &mut cc).is_nan(), "pow(1, NaN)");
    assert!(NAN.pow(&zero, P, RM, &mut cc).is_nan(), "pow(NaN, 0)");
    assert!(
        NAN.pow(&zero.neg(), P, RM, &mut cc).is_nan(),
        "pow(NaN, -0)"
    );
    assert!(NAN.powi(0, P, RM).is_nan(), "powi(NaN, 0)");
    assert!(NAN.powsi(0, P, RM).is_nan(), "powsi(NaN, 0)");
    Ok(())
}
