//! Correctly rounded cos, sin, tan, atan, ln, and pow against MPFR.
//!
//! Arguments are built from integers: ±1/3, ±3/4, ±5·2^−20, ±7, ±(2^20+1).
//! Rounding modes are the six modes `exact_bit_functions` sweeps (ToEven, ToOdd, Up, Down,
//! ToZero, FromZero). `RoundingMode::None` has no MPFR counterpart.
//!
//! Pow uses the exponent 1/3. That value is not an integer and, at a full word of precision,
//! six exponent increments do not land on an integer, so the dyadic-root shortcut does not
//! apply. The general path is `exp(y · ln x)`. A non-integer power of `e` is split into an
//! integer part and a fractional part, and the fractional part is `expf`, which calls
//! `sinh_series`. A negative base to this exponent is `InvalidArgument` (NaN), as in MPFR.

use crate::mpfr::common::assert_exact_eq;
use crate::mpfr::common::conv_to_mpfr;
use crate::mpfr::common::mpfr_reference_with;
use gmp_mpfr_sys::mpfr::{self, rnd_t};
use rug::Float;
use zenith_float_num::{Consts, ExactNum, RoundingMode, WORD_BIT_SIZE};

const MODES: [RoundingMode; 6] = [
    RoundingMode::ToEven,
    RoundingMode::ToOdd,
    RoundingMode::Up,
    RoundingMode::Down,
    RoundingMode::ToZero,
    RoundingMode::FromZero,
];

enum Spec {
    Ratio(i64, i64),
    Shift(i64, i64),
    Int(i64),
}

struct Case {
    label: &'static str,
    spec: Spec,
}

fn cases() -> [Case; 5] {
    [
        Case {
            label: "1/3",
            spec: Spec::Ratio(1, 3),
        },
        Case {
            label: "3/4",
            spec: Spec::Ratio(3, 4),
        },
        Case {
            label: "5*2^-20",
            spec: Spec::Shift(5, -20),
        },
        Case {
            label: "7",
            spec: Spec::Int(7),
        },
        Case {
            label: "2^20+1",
            spec: Spec::Int(1_048_577),
        },
    ]
}

fn storage(p: usize) -> usize {
    p.div_ceil(WORD_BIT_SIZE) * WORD_BIT_SIZE
}

fn signed(n: i64, neg: bool) -> i64 {
    if !neg {
        n
    } else {
        match n.checked_neg() {
            Some(v) => v,
            None => {
                assert!(n > 0, "argument {n} has no negative");
                0
            }
        }
    }
}

fn build(spec: &Spec, neg: bool, p: usize) -> ExactNum {
    let s = storage(p);
    let rm = RoundingMode::ToEven;
    match spec {
        Spec::Ratio(n, d) => {
            let num = ExactNum::from_i64(signed(*n, neg), s);
            let den = ExactNum::from_i64(*d, s);
            num.div(&den, s, rm)
        }
        Spec::Shift(n, exp) => ExactNum::from_i64(signed(*n, neg), s).ldexp(*exp, s, rm),
        Spec::Int(n) => ExactNum::from_i64(signed(*n, neg), s),
    }
}

fn each_precision(
    mut body: impl FnMut(usize) -> Result<(), Box<dyn std::error::Error>>,
) -> Result<(), Box<dyn std::error::Error>> {
    for p in 2..=256usize {
        body(p)?;
    }
    for p in [511usize, 512, 513, 1024, 1625, 4096] {
        body(p)?;
    }
    Ok(())
}

fn sweep_unary(
    name: &str,
    eval: impl Fn(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
    mp: impl Fn(&mut Float, &Float, rnd_t) -> i32,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut cc = Consts::new()?;
    let args = cases();
    each_precision(|p| {
        for case in &args {
            for neg in [false, true] {
                let x = build(&case.spec, neg, p);
                let fx = conv_to_mpfr(storage(p), &x, &mut cc)?;
                let p32 = match u32::try_from(p) {
                    Ok(v) => v,
                    Err(_) => return Err("precision does not fit u32".into()),
                };
                for rm in MODES {
                    let got = eval(&x, p, rm, &mut cc);
                    let (want, _) = mpfr_reference_with(p32, rm, |r, rnd| mp(r, &fx, rnd))?;
                    let sign = if neg { "-" } else { "+" };
                    let info = format!("{name}({sign}{}) p={p} {rm:?}", case.label);
                    assert_exact_eq(&got, &want, p, &info, &mut cc)?;
                }
            }
        }
        Ok(())
    })
}

#[test]
fn mpfr_grid_cos() -> Result<(), Box<dyn std::error::Error>> {
    sweep_unary(
        "cos",
        |x, p, rm, cc| x.cos(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::cos(r.as_raw_mut(), a.as_raw(), rnd) },
    )
}

#[test]
fn mpfr_grid_sin() -> Result<(), Box<dyn std::error::Error>> {
    sweep_unary(
        "sin",
        |x, p, rm, cc| x.sin(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::sin(r.as_raw_mut(), a.as_raw(), rnd) },
    )
}

#[test]
fn mpfr_grid_tan() -> Result<(), Box<dyn std::error::Error>> {
    sweep_unary(
        "tan",
        |x, p, rm, cc| x.tan(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::tan(r.as_raw_mut(), a.as_raw(), rnd) },
    )
}

#[test]
fn mpfr_grid_atan() -> Result<(), Box<dyn std::error::Error>> {
    sweep_unary(
        "atan",
        |x, p, rm, cc| x.atan(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::atan(r.as_raw_mut(), a.as_raw(), rnd) },
    )
}

#[test]
fn mpfr_grid_ln() -> Result<(), Box<dyn std::error::Error>> {
    sweep_unary(
        "ln",
        |x, p, rm, cc| x.ln(p, rm, cc),
        |r, a, rnd| unsafe { mpfr::log(r.as_raw_mut(), a.as_raw(), rnd) },
    )
}

#[test]
fn mpfr_grid_pow_sinh_path() -> Result<(), Box<dyn std::error::Error>> {
    let mut cc = Consts::new()?;
    let args = cases();
    each_precision(|p| {
        let exponent = build(&Spec::Ratio(1, 3), false, p);
        let fy = conv_to_mpfr(storage(p), &exponent, &mut cc)?;
        for case in &args {
            for neg in [false, true] {
                let x = build(&case.spec, neg, p);
                let fx = conv_to_mpfr(storage(p), &x, &mut cc)?;
                for rm in MODES {
                    let got = x.pow(&exponent, p, rm, &mut cc);
                    let p32 = match u32::try_from(p) {
                        Ok(v) => v,
                        Err(_) => return Err("precision does not fit u32".into()),
                    };
                    let (want, _) = mpfr_reference_with(p32, rm, |r, rnd| unsafe {
                        mpfr::pow(r.as_raw_mut(), fx.as_raw(), fy.as_raw(), rnd)
                    })?;
                    let sign = if neg { "-" } else { "+" };
                    let info = format!("pow({sign}{}, 1/3) p={p} {rm:?}", case.label);
                    assert_exact_eq(&got, &want, p, &info, &mut cc)?;
                }
            }
        }
        Ok(())
    })
}
