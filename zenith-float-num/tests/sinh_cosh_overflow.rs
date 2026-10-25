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
