//! Correct rounding at very high precision when the result lies within about `2^-2p` of a rounding
//! boundary. A correctly rounded result there needs working precision of about `2p` bits, so the
//! retry budget has to grow with `p` (a fixed number of extra bits ends in
//! `NaN(PrecisionRetryExhausted)` once `p` is large). Every expected value comes from MPFR (through
//! `rug`), never from zenith-float itself.

use crate::mpfr::common::{assert_exact_eq, conv_to_mpfr, mpfr_reference_with};
use gmp_mpfr_sys::mpfr::rnd_t;
use rug::float::Round;
use rug::ops::AssignRound;
use rug::Float;
use zenith_float_num::{Consts, ExactNum, Exponent, RoundingMode};

/// Precision (a multiple of the word size) at which the original failure was found.
const P: usize = 100_032;

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
    assert!(
        !got.is_nan(),
        "{info} p={p} {rm:?}: zenith {got:?}, MPFR {w}"
    );
    assert_exact_eq(got, &w, p, &format!("{info} p={p} {rm:?}"), cc)?;
    Ok(())
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
    let fx = conv_to_mpfr(p, x, cc)?;
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

/// `expm1(2^(1-p))` and `log1p(2^(1-p))`, `ln(1 + 2^(1-p))`: a neighbour of the grid plus a `2^-2p` term.
#[test]
fn expm1_log1p_ln_near_a_grid_point_at_100032_bits() -> Res {
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
    )?;
    run_unary(
        "log1p",
        "2^(1-p)",
        P,
        &d,
        &mut cc,
        |x, p, rm, cc| x.log1p(p, rm, cc),
        |r, a, round| r.assign_round(a.ln_1p_ref(), round),
    )?;
    let y = one_plus_pow2(P, 1 - P as Exponent);
    run_unary(
        "ln",
        "1+2^(1-p)",
        P,
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
fn nth_root_just_below_a_tie_at_100032_bits() -> Res {
    let mut cc = Consts::new()?;
    let y4 = one_plus_pow2(P, 2 - P as Exponent);
    let three_ulp =
        pow2(P, 2 - P as Exponent).add(&pow2(P, 1 - P as Exponent), P, RoundingMode::ToEven);
    let y6 = ExactNum::from_u8(1, P).add(&three_ulp, P, RoundingMode::ToEven);
    for (n, y) in [(4_usize, &y4), (6, &y6)] {
        let fy = conv_to_mpfr(P, y, &mut cc)?;
        for rm in MODES {
            let got = y.nth_root(n, P, rm);
            check(&got, P, rm, &format!("nth_root({n})"), &mut cc, |r, rnd| {
                ternary(r.assign_round(fy.root_ref(n as u32), round_of(rnd)))
            })?;
        }
    }
    Ok(())
}

/// `log_b(x)` with `b = 1 + 2^(1-p)` and `x = 1 + 2^(2-p)`: `2 - 2^(1-p)` plus a term of about `2^-2p`.
#[test]
fn log_with_base_near_a_grid_point_at_100032_bits() -> Res {
    let mut cc = Consts::new()?;
    let b = one_plus_pow2(P, 1 - P as Exponent);
    let x = one_plus_pow2(P, 2 - P as Exponent);
    let fb = conv_to_mpfr(P, &b, &mut cc)?;
    let fx = conv_to_mpfr(P, &x, &mut cc)?;
    for rm in MODES {
        let got = x.log(&b, P, rm, &mut cc);
        check(
            &got,
            P,
            rm,
            "log(1+2^(2-p), base 1+2^(1-p))",
            &mut cc,
            |r, rnd| {
                // The quotient of two logarithms, each more precise than `p`; MPFR rounds the one division.
                let round = round_of(rnd);
                let wp = 2 * P as u32 + 128;
                let num = Float::with_val(wp, fx.ln_ref());
                let den = Float::with_val(wp, fb.ln_ref());
                ternary(r.assign_round(&num / &den, round))
            },
        )?;
    }
    Ok(())
}
