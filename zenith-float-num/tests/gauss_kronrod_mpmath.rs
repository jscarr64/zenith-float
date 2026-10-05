//! (7, 15) Gauss–Kronrod and Wynn ε against mpmath 1.4.1 closed forms.
//!
//! mpmath `quad` on the singular integrands stopped short of the antiderivative
//! at 60 decimals, so those golds are the closed forms evaluated by mpmath
//! (`e - 1`, `2`, `2/3`). SymPy is not an oracle.

use zenith_float_num::{
    gauss_kronrod_interval, gauss_legendre, integrate_adaptive_gk, Consts, ExactNum,
    QuadratureError, Radix, RoundingMode,
};

const RM: RoundingMode = RoundingMode::ToEven;
const P: usize = 128;

fn parse(s: &str, p: usize, cc: &mut Consts) -> ExactNum {
    ExactNum::parse(s, Radix::Dec, p, RM, cc)
}

fn bits(got: &ExactNum, want: &ExactNum, p: usize) -> i32 {
    if got.is_nan() {
        return i32::MIN;
    }
    let d = got.sub(want, 2 * p, RoundingMode::None);
    if d.is_zero() {
        return i32::MAX;
    }
    want.exponent().unwrap_or(0) - d.exponent().unwrap_or(0)
}

fn tol(bits_rel: i32) -> ExactNum {
    ExactNum::from_u8(1, P).ldexp(-bits_rel, P, RM)
}

fn check(label: &str, got: &ExactNum, want: &str, min_bits: i32, cc: &mut Consts) {
    let w = parse(want, 4 * P, cc);
    let bts = bits(got, &w, P);
    assert!(
        bts >= min_bits,
        "{label}: {bts} correct bits, wanted {min_bits}, got {got}"
    );
}

#[test]
fn gk_smooth_exp_matches_e_minus_one() {
    // mpmath 1.4.1, dps=60: exp(1) - 1. One (7, 15) panel is enough on a
    // short analytic interval; the Kronrod sum is the value under test.
    let mut cc = Consts::new().unwrap();
    let a = ExactNum::from_u8(0, P);
    let b = ExactNum::from_u8(1, P);
    let got = gauss_kronrod_interval(|x, p, rm, cc| x.exp(p, rm, cc), &a, &b, P, RM, &mut cc)
        .expect("exp");
    check(
        "int_0^1 exp",
        &got.estimate_kronrod,
        "1.7182818284590452353602874713526624977572470936999595749669676277",
        (P as i32) - 8,
        &mut cc,
    );
}

#[test]
fn gk_oscillatory_sin_matches_two() {
    // The antiderivative of sin on [0, pi] is 2. mpmath nstr(quad(sin, [0, pi]))
    // at 50 decimals is 2; the adaptive driver may bisect until |K-G| agrees.
    let mut cc = Consts::new().unwrap();
    let a = ExactNum::from_u8(0, P);
    let b = cc.pi(P, RM);
    let rel = tol(48);
    let got = integrate_adaptive_gk(
        |x, p, rm, cc| x.sin(p, rm, cc),
        &a,
        &b,
        &rel,
        16,
        P,
        RM,
        &mut cc,
    )
    .expect("sin");
    check("int_0^pi sin", &got.estimate_kronrod, "2", 96, &mut cc);
}

#[test]
fn gk_mild_endpoint_sqrt_matches_two_thirds() {
    // Hand antiderivative of sqrt(x) on [0, 1] is 2/3. mpmath at 50 decimals
    // prints the same repeating expansion. The derivative blows up at 0, so
    // the driver bisects that panel instead of trusting one Gauss rule.
    let mut cc = Consts::new().unwrap();
    let a = ExactNum::from_u8(0, P);
    let b = ExactNum::from_u8(1, P);
    let rel = tol(40);
    let got = integrate_adaptive_gk(
        |x, p, rm, _cc| x.sqrt(p, rm),
        &a,
        &b,
        &rel,
        16,
        P,
        RM,
        &mut cc,
    )
    .expect("sqrt");
    check(
        "int_0^1 sqrt(x)",
        &got.estimate_kronrod,
        "0.6666666666666666666666666666666666666666666666666666666666666667",
        56,
        &mut cc,
    );
}

#[test]
fn gk_epsilon_beats_plain_gauss_on_invsqrt() {
    // Hand antiderivative of x^(-1/2) on (0, 1] is 2. A single 7-point Gauss
    // panel misses by several hundredths because the endpoint is not a
    // polynomial. Wynn epsilon on the adaptive Kronrod sums cancels that
    // nearly geometric panel error. mpmath's closed form at 50 decimals is 2.
    let mut cc = Consts::new().unwrap();
    let a = ExactNum::from_u8(0, P);
    let b = ExactNum::from_u8(1, P);
    let plain = gauss_legendre(
        |x, p, rm, _cc| x.sqrt(p, rm).reciprocal(p, rm),
        &a,
        &b,
        7,
        P,
        RM,
        &mut cc,
    )
    .expect("gauss7");
    let rel = tol(40);
    let adaptive = integrate_adaptive_gk(
        |x, p, rm, _cc| x.sqrt(p, rm).reciprocal(p, rm),
        &a,
        &b,
        &rel,
        12,
        P,
        RM,
        &mut cc,
    )
    .expect("adaptive invsqrt");
    let two = parse("2", 4 * P, &mut cc);
    let plain_bits = bits(&plain, &two, P);
    let adapt_bits = bits(&adaptive.estimate_kronrod, &two, P);
    assert!(
        adapt_bits >= plain_bits + 32,
        "epsilon/adaptive bits {adapt_bits} should beat plain Gauss bits {plain_bits}"
    );
    assert!(
        adapt_bits >= 48,
        "adaptive invsqrt bits {adapt_bits}, value {}",
        adaptive.estimate_kronrod
    );
}

#[test]
fn gk_nonfinite_and_reversed_bounds_are_errors() {
    // A reversed interval and a NaN sample must not come back as a finite
    // integral. The typed error is the whole assertion.
    let mut cc = Consts::new().unwrap();
    let a = ExactNum::from_u8(0, P);
    let b = ExactNum::from_u8(1, P);
    let reversed = gauss_kronrod_interval(|x, _p, _rm, _cc| x.clone(), &b, &a, P, RM, &mut cc);
    assert_eq!(reversed.unwrap_err(), QuadratureError::InvalidArgument);
    let nan = gauss_kronrod_interval(
        |_x, _p, _rm, _cc| ExactNum::nan(None),
        &a,
        &b,
        P,
        RM,
        &mut cc,
    );
    assert_eq!(nan.unwrap_err(), QuadratureError::InvalidArgument);
}
