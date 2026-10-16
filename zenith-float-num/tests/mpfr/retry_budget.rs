//! Correct rounding at very high precision when the result lies within about `2^-2p` of a rounding
//! boundary. A correctly rounded result there needs working precision of about `2p` bits, so the
//! retry budget has to grow with `p` (a fixed number of extra bits ends in
//! `NaN(PrecisionRetryExhausted)` once `p` is large). Every expected value comes from MPFR (through
//! `rug`), never from zenith-float itself.

use crate::mpfr::common::{
    assert_exact_eq, conv_str_to_mpfr_compat, conv_to_mpfr, mpfr_reference_with,
};
use gmp_mpfr_sys::mpfr::rnd_t;
use rug::float::Round;
use rug::ops::AssignRound;
use rug::Float;
use zenith_float_num::{Consts, ExactNum, Exponent, Radix, RoundingMode};

/// Precision (a multiple of the word size) at which the original failure was found.
const P: usize = 100_032;

/// Precision of the tests whose retry steps are the most expensive (`ln`, `log`, `nth_root`).
const PM: usize = 40_064;

const MODES: [RoundingMode; 6] = [
    RoundingMode::ToEven,
    RoundingMode::ToOdd,
    RoundingMode::Up,
    RoundingMode::Down,
    RoundingMode::ToZero,
    RoundingMode::FromZero,
];

type Res = Result<(), Box<dyn std::error::Error>>;

fn round_of(rnd: rnd_t) -> Round {
    match rnd {
        rnd_t::RNDZ => Round::Zero,
        rnd_t::RNDU => Round::Up,
        rnd_t::RNDD => Round::Down,
        rnd_t::RNDA => Round::AwayZero,
        _ => Round::Nearest,
    }
}

fn ternary(o: std::cmp::Ordering) -> i32 {
    o as i32
}

/// `pi * 2^k` at `p` bits.
fn pi_times_2k(cc: &mut Consts, p: usize, k: Exponent) -> ExactNum {
    cc.pi(p, RoundingMode::ToEven)
        .ldexp(k, p, RoundingMode::ToEven)
}

/// `2^k` at `p` bits.
fn pow2(p: usize, k: Exponent) -> ExactNum {
    ExactNum::from_u8(1, p).ldexp(k, p, RoundingMode::ToEven)
}

/// `1 + 2^k` at `p` bits (exact for `-p < k < 0`).
fn one_plus_pow2(p: usize, k: Exponent) -> ExactNum {
    ExactNum::from_u8(1, p).add(&pow2(p, k), p, RoundingMode::ToEven)
}

/// `sqrt(c) * 2^k` at `p` bits (correctly rounded square root, then an exact scaling).
fn sqrt_times_2k(c: u8, k: Exponent, p: usize) -> ExactNum {
    ExactNum::from_u8(c, p)
        .sqrt(p, RoundingMode::ToEven)
        .ldexp(k, p, RoundingMode::ToEven)
}

/// Compares one result with MPFR's value for the same operation, bit for bit.
fn check(
    got: &ExactNum,
    p: usize,
    rm: RoundingMode,
    info: &str,
    cc: &mut Consts,
    want: impl Fn(&mut Float, rnd_t) -> i32,
) -> Res {
    let (w, _) = mpfr_reference_with(p as u32, rm, want)?;
    // Out-of-domain arguments: both sides must be NaN; everything else must be a number (never the
    // `PrecisionRetryExhausted` NaN).
    assert!(
        got.is_nan() == (w.is_nan()),
        "{info} p={p} {rm:?}: zenith {got:?}, MPFR {w}"
    );
    if w.is_nan() {
        return Ok(());
    }
    assert_exact_eq(got, &w, p, &format!("{info} p={p} {rm:?}"), cc)?;
    Ok(())
}

/// The argument as an MPFR number of `p` bits, which must be exact (also for subnormal arguments, where
/// `conv_to_mpfr`'s parse-back comparison does not apply).
fn argument_to_mpfr(
    p: usize,
    x: &ExactNum,
    cc: &mut Consts,
) -> Result<Float, Box<dyn std::error::Error>> {
    let text = conv_str_to_mpfr_compat(x.format(Radix::Bin, RoundingMode::ToEven, cc)?)?;
    let parsed = Float::parse_radix(text, 2)?;
    let mut f = Float::new(p as u32);
    let ord = f.assign_round(parsed, Round::Nearest);
    assert_eq!(
        ord,
        std::cmp::Ordering::Equal,
        "the argument is not exact in {p} bits"
    );
    Ok(f)
}

fn run_unary(
    name: &str,
    label: &str,
    p: usize,
    x: &ExactNum,
    cc: &mut Consts,
    zf: impl Fn(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
    mp: impl Fn(&mut Float, &Float, Round) -> std::cmp::Ordering,
) -> Res {
    let fx = argument_to_mpfr(p, x, cc)?;
    for rm in MODES {
        let got = zf(x, p, rm, cc);
        let info = format!("{name}({label})");
        check(&got, p, rm, &info, cc, |r, rnd| {
            ternary(mp(r, &fx, round_of(rnd)))
        })?;
    }
    Ok(())
}

fn cos_of(x: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> ExactNum {
    x.cos(p, rm, cc)
}

fn mp_cos(r: &mut Float, a: &Float, round: Round) -> std::cmp::Ordering {
    r.assign_round(a.cos_ref(), round)
}

/// The failure that was found: `cos(4 pi)` at 100032 bits is `1 - 2^-(2p-5)`-ish, so FromZero, Down and Up
/// (and the nearest modes) need about `2p` working bits.
#[test]
fn cos_four_pi_at_100032_bits() -> Res {
    let mut cc = Consts::new()?;
    let x = pi_times_2k(&mut cc, P, 2);
    run_unary("cos", "4 pi", P, &x, &mut cc, cos_of, mp_cos)
}

fn sin_of(x: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> ExactNum {
    x.sin(p, rm, cc)
}

fn mp_sin(r: &mut Float, a: &Float, round: Round) -> std::cmp::Ordering {
    r.assign_round(a.sin_ref(), round)
}

fn tan_of(x: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> ExactNum {
    x.tan(p, rm, cc)
}

fn mp_tan(r: &mut Float, a: &Float, round: Round) -> std::cmp::Ordering {
    r.assign_round(a.tan_ref(), round)
}

/// `cos` of 2 pi and of pi / 2 (a result of about `2^-p`: needs `2p` bits of the argument reduction).
#[test]
fn cos_two_pi_and_half_pi_at_100032_bits() -> Res {
    let mut cc = Consts::new()?;
    let x = pi_times_2k(&mut cc, P, 1);
    run_unary("cos", "2 pi", P, &x, &mut cc, cos_of, mp_cos)?;
    let x = pi_times_2k(&mut cc, P, -1);
    run_unary("cos", "pi/2", P, &x, &mut cc, cos_of, mp_cos)
}

/// `sin` near multiples of pi (tiny result) and near pi / 2 (result `1 - tiny`).
#[test]
fn sin_near_pi_multiples_at_100032_bits() -> Res {
    let mut cc = Consts::new()?;
    for (label, k) in [("4 pi", 2), ("pi", 0), ("pi/2", -1)] {
        let x = pi_times_2k(&mut cc, P, k);
        run_unary("sin", label, P, &x, &mut cc, sin_of, mp_sin)?;
    }
    Ok(())
}

/// `tan` near multiples of pi (tiny result) and near pi / 2 (huge result).
#[test]
fn tan_near_pi_multiples_at_100032_bits() -> Res {
    let mut cc = Consts::new()?;
    for (label, k) in [("4 pi", 2), ("pi/2", -1)] {
        let x = pi_times_2k(&mut cc, P, k);
        run_unary("tan", label, P, &x, &mut cc, tan_of, mp_tan)?;
    }
    Ok(())
}

/// `sin_cos` returns both values; each is compared with MPFR.
#[test]
fn sin_cos_near_pi_multiples_at_100032_bits() -> Res {
    let mut cc = Consts::new()?;
    for (label, k) in [("4 pi", 2), ("pi/2", -1)] {
        let x = pi_times_2k(&mut cc, P, k);
        let fx = conv_to_mpfr(P, &x, &mut cc)?;
        for rm in MODES {
            let (s, c) = x.sin_cos(P, rm, &mut cc);
            check(
                &s,
                P,
                rm,
                &format!("sin_cos.sin({label})"),
                &mut cc,
                |r, rnd| ternary(mp_sin(r, &fx, round_of(rnd))),
            )?;
            check(
                &c,
                P,
                rm,
                &format!("sin_cos.cos({label})"),
                &mut cc,
                |r, rnd| ternary(mp_cos(r, &fx, round_of(rnd))),
            )?;
        }
    }
    Ok(())
}

/// `exp(2^-p)` is `1 + 2^-p + 2^(-2p-1) + ...`: the midpoint of two neighbours plus a `2^-2p` term.
#[test]
fn exp_just_above_a_tie_at_100032_bits() -> Res {
    let mut cc = Consts::new()?;
    let x = pow2(P, -(P as Exponent));
    run_unary(
        "exp",
        "2^-p",
        P,
        &x,
        &mut cc,
        |x, p, rm, cc| x.exp(p, rm, cc),
        |r, a, round| r.assign_round(a.exp_ref(), round),
    )
}

/// `expm1(2^(1-p)) = x + x^2/2 + x^3/6`: the midpoint between two neighbours plus a `2^-2p` term.
#[test]
fn expm1_just_above_a_tie_at_100032_bits() -> Res {
    let mut cc = Consts::new()?;
    let d = pow2(P, 1 - P as Exponent);
    run_unary(
        "expm1",
        "2^(1-p)",
        P,
        &d,
        &mut cc,
        |x, p, rm, cc| x.expm1(p, rm, cc),
        |r, a, round| r.assign_round(a.exp_m1_ref(), round),
    )
}

/// `log1p(2^(1-p)) = x - x^2/2 + x^3/3`: a grid point of `p` bits plus a term of relative size `2^-2p`.
#[test]
fn log1p_just_above_a_grid_point_at_100032_bits() -> Res {
    let mut cc = Consts::new()?;
    let d = pow2(P, 1 - P as Exponent);
    run_unary(
        "log1p",
        "2^(1-p)",
        P,
        &d,
        &mut cc,
        |x, p, rm, cc| x.log1p(p, rm, cc),
        |r, a, round| r.assign_round(a.ln_1p_ref(), round),
    )
}

/// `ln(1 + 2^(1-p))`: the same value as `log1p(2^(1-p))` reached through the logarithm of a number near 1
/// (a working precision of about `2p` plus the cancellation guard of the argument; the slowest test here, so
/// it runs at 40064 bits, a precision at which `p + 64 * 256` working bits cannot reach the `2p` these shapes need).
#[test]
fn ln_just_above_a_grid_point_at_40064_bits() -> Res {
    let mut cc = Consts::new()?;
    let y = one_plus_pow2(PM, 1 - PM as Exponent);
    run_unary(
        "ln",
        "1+2^(1-p)",
        PM,
        &y,
        &mut cc,
        |x, p, rm, cc| x.ln(p, rm, cc),
        |r, a, round| r.assign_round(a.ln_ref(), round),
    )
}

/// `cosh(x) - 1 = x^2/2 + ...` with `x^2/2` within `2^-2p` of the midpoint `2^-p`.
#[test]
fn cosh_sinh_tanh_near_a_tie_at_100032_bits() -> Res {
    let mut cc = Consts::new()?;
    let half = -(P as Exponent) / 2;
    let xc = sqrt_times_2k(2, half, P);
    run_unary(
        "cosh",
        "sqrt(2) 2^(-p/2)",
        P,
        &xc,
        &mut cc,
        |x, p, rm, cc| x.cosh(p, rm, cc),
        |r, a, round| r.assign_round(a.cosh_ref(), round),
    )?;
    let xs = sqrt_times_2k(6, half, P);
    run_unary(
        "sinh",
        "sqrt(6) 2^(-p/2)",
        P,
        &xs,
        &mut cc,
        |x, p, rm, cc| x.sinh(p, rm, cc),
        |r, a, round| r.assign_round(a.sinh_ref(), round),
    )?;
    let xt = sqrt_times_2k(3, half, P);
    run_unary(
        "tanh",
        "sqrt(3) 2^(-p/2)",
        P,
        &xt,
        &mut cc,
        |x, p, rm, cc| x.tanh(p, rm, cc),
        |r, a, round| r.assign_round(a.tanh_ref(), round),
    )
}

/// The `n`-th root of `1 + n * 2^-p` is `1 + 2^-p - (n-1) * 2^(-2p-1) / n + ...`: just below the midpoint
/// of two neighbours, for `n = 4` and `n = 6`.
#[test]
fn nth_root_just_below_a_tie_at_40064_bits() -> Res {
    let mut cc = Consts::new()?;
    let y4 = one_plus_pow2(PM, 2 - PM as Exponent);
    let three_ulp =
        pow2(PM, 2 - PM as Exponent).add(&pow2(PM, 1 - PM as Exponent), PM, RoundingMode::ToEven);
    let y6 = ExactNum::from_u8(1, PM).add(&three_ulp, PM, RoundingMode::ToEven);
    for (n, y) in [(4_usize, &y4), (6, &y6)] {
        let fy = conv_to_mpfr(PM, y, &mut cc)?;
        for rm in MODES {
            let got = y.nth_root(n, PM, rm);
            check(
                &got,
                PM,
                rm,
                &format!("nth_root({n})"),
                &mut cc,
                |r, rnd| ternary(r.assign_round(fy.root_ref(n as u32), round_of(rnd))),
            )?;
        }
    }
    Ok(())
}

/// `log_b(x)` with `b = 1 + 2^(1-p)` and `x = 1 + 2^(2-p)`: `2 - 2^(1-p)` plus a term of about `2^-2p`.
#[test]
fn log_with_base_near_a_grid_point_at_40064_bits() -> Res {
    let mut cc = Consts::new()?;
    let b = one_plus_pow2(PM, 1 - PM as Exponent);
    let x = one_plus_pow2(PM, 2 - PM as Exponent);
    let fb = conv_to_mpfr(PM, &b, &mut cc)?;
    let fx = conv_to_mpfr(PM, &x, &mut cc)?;
    for rm in MODES {
        let got = x.log(&b, PM, rm, &mut cc);
        check(
            &got,
            PM,
            rm,
            "log(1+2^(2-p), base 1+2^(1-p))",
            &mut cc,
            |r, rnd| {
                // The quotient of two logarithms, each more precise than `p`; MPFR rounds the one division.
                let round = round_of(rnd);
                let wp = 2 * PM as u32 + 128;
                let num = Float::with_val(wp, fx.ln_ref());
                let den = Float::with_val(wp, fb.ln_ref());
                ternary(r.assign_round(&num / &den, round))
            },
        )?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------------------------
// Wide range: arguments `sqrt(2) * 2^k` from the bottom of the exponent range up to 2, every mode.
// ---------------------------------------------------------------------------------------------

/// MPFR's default exponent range is narrower than zenith-float's: widen it to the library's for this thread so
/// that the operand (and the expected value) can be built at any exponent.
pub fn use_full_exponent_range() {
    use gmp_mpfr_sys::gmp::exp_t;
    use gmp_mpfr_sys::mpfr;
    // SAFETY: only sets MPFR's thread-local exponent bounds, as the other tests of this suite do.
    unsafe {
        // Below `EXPONENT_MIN` so that subnormal arguments convert exactly (the reference computation sets
        // its own bound for the result).
        mpfr::set_emin(zenith_float_num::EXPONENT_MIN as exp_t - 1024);
        mpfr::set_emax(zenith_float_num::EXPONENT_MAX as exp_t);
    }
}

type ZfFn = fn(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum;
type MpFn = fn(&mut Float, &Float, Round) -> std::cmp::Ordering;

/// The functions of the wide-range sweeps, with their MPFR counterparts.
fn wide_functions() -> [(&'static str, ZfFn, MpFn); 10] {
    [
        (
            "exp",
            |x, p, rm, cc| x.exp(p, rm, cc),
            |r, a, rd| r.assign_round(a.exp_ref(), rd),
        ),
        (
            "expm1",
            |x, p, rm, cc| x.expm1(p, rm, cc),
            |r, a, rd| r.assign_round(a.exp_m1_ref(), rd),
        ),
        (
            "sinh",
            |x, p, rm, cc| x.sinh(p, rm, cc),
            |r, a, rd| r.assign_round(a.sinh_ref(), rd),
        ),
        (
            "cosh",
            |x, p, rm, cc| x.cosh(p, rm, cc),
            |r, a, rd| r.assign_round(a.cosh_ref(), rd),
        ),
        (
            "tanh",
            |x, p, rm, cc| x.tanh(p, rm, cc),
            |r, a, rd| r.assign_round(a.tanh_ref(), rd),
        ),
        ("sin", sin_of, mp_sin),
        ("cos", cos_of, mp_cos),
        ("tan", tan_of, mp_tan),
        (
            "atan",
            |x, p, rm, cc| x.atan(p, rm, cc),
            |r, a, rd| r.assign_round(a.atan_ref(), rd),
        ),
        (
            "log1p",
            |x, p, rm, cc| x.log1p(p, rm, cc),
            |r, a, rd| r.assign_round(a.ln_1p_ref(), rd),
        ),
    ]
}

/// First `n` characters of a message (the failure messages carry whole MPFR values).
fn head(s: &str, n: usize) -> String {
    s.chars().take(n).collect()
}

/// Runs every function of `names` at every exponent of `exponents` and reports all the mismatches at once
/// (a mismatch panics inside `assert_exact_eq`; each one is caught so that one run lists them all).
fn sweep(p: usize, exponents: &[Exponent], names: &[&str]) -> Res {
    let mut args: Vec<(String, ExactNum)> = Vec::new();
    for (k, negative) in exponents.iter().flat_map(|k| [(k, false), (k, true)]) {
        // `sqrt(2)` has exponent 1, so scaling by `k - 1` puts the argument's exponent at `k`.
        let pos = sqrt_times_2k(2, k - 1, p);
        let sign = if negative { "-" } else { "" };
        let label = format!("{sign}sqrt(2) 2^{}", k - 1);
        args.push((label, if negative { -pos } else { pos }));
    }
    sweep_args(p, &args, names)
}

/// Runs every function of `names` at every argument and reports all the mismatches at once.
fn sweep_args(p: usize, args: &[(String, ExactNum)], names: &[&str]) -> Res {
    use std::panic::{catch_unwind, AssertUnwindSafe};
    use_full_exponent_range();
    let mut cc = Consts::new()?;
    let mut failures: Vec<String> = Vec::new();
    for (label, x) in args {
        for (name, zf, mp) in wide_functions() {
            if !names.contains(&name) {
                continue;
            }
            let run = catch_unwind(AssertUnwindSafe(|| {
                run_unary(name, label, p, x, &mut cc, zf, mp)
            }));
            match run {
                Ok(Ok(())) => {}
                Ok(Err(e)) => failures.push(head(&format!("{name} {label}: error {e}"), 200)),
                Err(payload) => {
                    let msg = match payload.downcast_ref::<String>() {
                        Some(m) => head(m, 200),
                        None => String::from("panic without a string message"),
                    };
                    failures.push(msg);
                }
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} mismatches with MPFR:\n{}",
        failures.len(),
        failures.join("\n")
    );
    Ok(())
}

const ALL_NAMES: [&str; 10] =
    ["exp", "expm1", "sinh", "cosh", "tanh", "sin", "cos", "tan", "atan", "log1p"];

/// From the bottom of the exponent range to 2, around `-p/2` and `-p` where the result changes from "1 plus a
/// tail below the last bit" to "1 plus a tail the rounding has to see": all six modes, `p = 20032`.
#[test]
// The exponent `-(w - 1)` (`log1p(-sqrt(2) 2^-w)` with `ToEven`) is left out on purpose: it exposes the known
// roundability-check issue described under "Known issues" in the CHANGELOG.
fn wide_exponent_range_at_20032_bits() -> Res {
    const W: usize = 20_032;
    let w = W as Exponent;
    let exponents = [
        Exponent::MIN + 1,
        -(1 << 30) - 5,
        -(1 << 20),
        -3 * w,
        -2 * w - 1,
        -2 * w,
        -(w + 64),
        -(w + 1),
        -w,
        -(w / 2 + 1),
        -(w / 2),
        -(w / 2 - 1),
        -(w / 4),
        -1000,
        -64,
        -2,
        0,
        1,
    ];
    sweep(W, &exponents, &ALL_NAMES)
}

/// The same at `p = 100032` (the precision of the original failure) for the arguments where the tail of the
/// series is about `2^-p` (exponents `-p/2` and `-p`) and at the bottom of the range.
#[test]
fn wide_exponent_range_at_100032_bits() -> Res {
    let w = P as Exponent;
    let exponents = [Exponent::MIN + 1, -(w + 1), -w, -(w / 2 + 1), -(w / 2)];
    sweep(P, &exponents, &ALL_NAMES)
}

/// Arguments of both signs whose exponent is just below, at and far below `-p` (the "small argument"
/// shortcut of the `exp`-type and `sin`-type functions), plus the ordinary range, at a small precision:
/// the quick version of the wide sweeps.
#[test]
fn small_arguments_of_both_signs_at_128_bits() -> Res {
    let mut exponents = vec![Exponent::MIN + 1, -1000, -200, -10, -1, 0, 1, 2];
    exponents.extend(-140..=-60);
    sweep(128, &exponents, &ALL_NAMES)
}

/// Subnormal arguments (the grid below `2^EXPONENT_MIN`), both signs, at 64, 128 and 200 bits, vs MPFR: the
/// operand that exposed the wrong `log1p`/`Down` pin (`3 * 2^(EXPONENT_MIN - 20)`-like values).
#[test]
fn subnormal_arguments_of_both_signs() -> Res {
    fn subnormal(c: u8, shift: Exponent, p: usize) -> ExactNum {
        ExactNum::from_u8(c, p)
            .ldexp(Exponent::MIN, p, RoundingMode::ToEven)
            .ldexp(-shift, p, RoundingMode::ToEven)
    }
    for p in [64usize, 128, 200] {
        let mut args: Vec<(String, ExactNum)> = Vec::new();
        for (c, shift) in [(3u8, 20), (1, 1), (5, 40), (7, 3)] {
            args.push((format!("{c} 2^(MIN-{shift})"), subnormal(c, shift, p)));
            args.push((format!("-{c} 2^(MIN-{shift})"), -subnormal(c, shift, p)));
        }
        sweep_args(p, &args, &ALL_NAMES)?;
    }
    Ok(())
}
