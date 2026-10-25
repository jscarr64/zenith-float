//! sinh/cosh overflow by mode near the i64 exponent limit (replaces the internal
//! `test_sinh_cosh` ExponentOverflow asserts removed with `ExactNumNumber::sinh_cosh` in b73337a):
//! a positive overflowing result is `+max` under `ToZero` and `Down`, `+inf` otherwise; a negative
//! one is `-max` under `ToZero` and `Up`, `-inf` otherwise.

use std::ops::Neg;
use zenith_float_num::{Consts, ExactNum, RoundingMode, INF_NEG, INF_POS};

const MODES: [RoundingMode; 6] = [
    RoundingMode::ToEven,
    RoundingMode::ToOdd,
    RoundingMode::Up,
    RoundingMode::Down,
    RoundingMode::ToZero,
    RoundingMode::FromZero,
];

fn expected(neg: bool, rm: RoundingMode, p: usize) -> ExactNum {
    let to_max = match rm {
        RoundingMode::ToZero => true,
        RoundingMode::Down => !neg,
        RoundingMode::Up => neg,
        _ => false,
    };
    match (to_max, neg) {
        (true, false) => ExactNum::max_value(p),
        (true, true) => ExactNum::max_value(p).neg(),
        (false, false) => INF_POS,
        (false, true) => INF_NEG,
    }
}

#[test]
fn sinh_cosh_overflow_is_directional_in_every_mode_and_sign() {
    let mut cc = Consts::new().unwrap();
    let mut bad = Vec::new();
    let mut checked = 0;
    for p in [64usize, 128, 320] {
        // 3 * 2^60 = 0.75 * 2^62: e^x has exponent about 1.08 * EXPONENT_MAX, just past the limit.
        let xs = [
            (
                "3*2^60",
                ExactNum::from_u8(3, p).ldexp(60, p, RoundingMode::None),
            ),
            (
                "2^63",
                ExactNum::from_u8(1, p).ldexp(63, p, RoundingMode::None),
            ),
            ("max", ExactNum::max_value(p)),
        ];
        for (label, x) in xs {
            for neg in [false, true] {
                let x = if neg { x.clone().neg() } else { x.clone() };
                for rm in MODES {
                    for (name, got, rneg) in [
                        ("sinh", x.sinh(p, rm, &mut cc), neg),
                        ("cosh", x.cosh(p, rm, &mut cc), false),
                    ] {
                        let want = expected(rneg, rm, p);
                        let ok = !got.is_nan()
                            && got.is_inf() == want.is_inf()
                            && got.compare(&want) == Some(0);
                        checked += 1;
                        if !ok {
                            bad.push(format!(
                                "{name}({}{label}) p={p} {rm:?}: got {got:?}, want {want:?}",
                                if neg { "-" } else { "+" }
                            ));
                        }
                    }
                }
            }
        }
    }
    assert_eq!(checked, 216);
    assert!(
        bad.is_empty(),
        "{} of {checked} wrong:\n{}",
        bad.len(),
        bad.join("\n")
    );
}

/// The rest of the old `test_sinh_cosh`: sinh(0) = 0 and cosh(0) = 1 exactly in every mode; for the
/// least positive subnormal x, sinh(x) rounds by mode next to x and cosh(x) next to 1.
#[test]
fn sinh_cosh_at_zero_and_at_the_least_subnormal() {
    let mut cc = Consts::new().unwrap();
    for p in [64usize, 320] {
        let one = ExactNum::from_u8(1, p);
        let zero = ExactNum::from_u8(0, p);
        let tiny = ExactNum::min_positive(p);
        for rm in MODES {
            let s = zero.sinh(p, rm, &mut cc);
            assert!(s.is_zero() && !s.inexact(), "sinh(0) p={p} {rm:?}: {s:?}");
            let c = zero.cosh(p, rm, &mut cc);
            assert!(
                c.compare(&one) == Some(0) && !c.inexact(),
                "cosh(0) p={p} {rm:?}: {c:?}"
            );
            // sinh(x) = x + x^3/6: x itself, or its successor when rounding up.
            let s = tiny.sinh(p, rm, &mut cc);
            let up = matches!(rm, RoundingMode::Up | RoundingMode::FromZero);
            let c = s.compare(&tiny);
            assert!(
                s.inexact() && c == Some(if up { 1 } else { 0 }),
                "sinh(tiny) p={p} {rm:?}: {s:?}"
            );
            // cosh(x) = 1 + x^2/2: 1 itself, or the successor of 1 when rounding up (ToOdd is
            // round-half-to-odd, so nearest).
            let c = tiny.cosh(p, rm, &mut cc);
            let k = c.compare(&one);
            assert!(
                c.inexact() && k == Some(if up { 1 } else { 0 }),
                "cosh(tiny) p={p} {rm:?}: {c:?}"
            );
        }
    }
}
