//! Regressions for 1.0.5: `ExactNum::cmp` / `abs_cmp` return exactly -1/0/1, which the
//! quadrature, root, ODE and Chebyshev routines rely on (`== Some(-1)`). In 1.0.4 equal-exponent
//! comparisons returned the mantissa word difference, so e.g. intervals like `[2, 3]` were
//! rejected and `rk45_adaptive` with a non-zero `rtol` never accepted a step. Root finders now
//! also stop (instead of returning `None`) when `tol` is below the working resolution.

use zenith_float_num::{
    bisect, brent, chebyshev_coeffs, chebyshev_eval, gauss_legendre, illinois, rk4, rk45_adaptive,
    root_default_tol, tanh_sinh, Consts, ExactNum, Radix, RoundingMode,
};

const RM: RoundingMode = RoundingMode::ToEven;
const P: usize = 128;

fn n(v: u8) -> ExactNum {
    ExactNum::from_u8(v, P)
}

/// `|got − want| < 2^-bits · max(1, |want|)`.
fn close(got: &ExactNum, want: &ExactNum, bits: i32) -> bool {
    let d = got.sub(want, 4 * P, RoundingMode::None);
    if d.is_zero() {
        return true;
    }
    let scale = want.exponent().unwrap_or(0).max(1);
    d.exponent().unwrap_or(i32::MAX) <= scale - bits
}

#[test]
fn cmp_is_normalized() {
    let (two, three) = (n(2), n(3));
    assert_eq!(three.cmp(&two), Some(1));
    assert_eq!(two.cmp(&three), Some(-1));
    assert_eq!(three.cmp(&three), Some(0));
    assert_eq!(three.neg().abs_cmp(&two), Some(1));
    assert_eq!(two.abs_cmp(&three.neg()), Some(-1));
}

#[test]
fn quadrature_equal_exponent_interval() {
    let mut cc = Consts::new().unwrap();
    let want = ExactNum::parse("1.5", Radix::Dec, 4 * P, RM, &mut cc).ln(4 * P, RM, &mut cc);
    let recip = |x: &ExactNum, p: usize, rm: RoundingMode, _: &mut Consts| {
        ExactNum::from_u8(1, p).div(x, p, rm)
    };
    let g = gauss_legendre(recip, &n(2), &n(3), 30, P, RM, &mut cc).expect("gauss_legendre [2,3]");
    assert!(close(&g, &want, 120));
    let t = tanh_sinh(recip, &n(2), &n(3), P, RM, &mut cc).expect("tanh_sinh [2,3]");
    assert!(close(&t, &want, 120));
}

#[test]
fn roots_default_tol_and_equal_exponent_bracket() {
    let mut cc = Consts::new().unwrap();
    let tol = root_default_tol(P, RM);
    // Brent at 128 bits with the default 2^-256 tolerance (1.0.4: None).
    let dottie = ExactNum::parse(
        "0.739085133215160641655312087673873404013411758900757465",
        Radix::Dec,
        4 * P,
        RM,
        &mut cc,
    );
    let r = brent(
        |x, p, rm, cc| x.cos(p, rm, cc).sub(x, p, rm),
        &n(0),
        &n(1),
        &tol,
        P,
        RM,
        &mut cc,
    )
    .expect("brent cos x = x");
    assert!(close(&r, &dottie, 125));
    // Bracket [2, 3] (1.0.4: rejected as not a < b).
    let six = n(6);
    let r = bisect(
        |x, p, rm, _| x.mul(x, p, rm).sub(&ExactNum::from_u8(6, p), p, rm),
        &n(2),
        &n(3),
        &tol,
        P,
        RM,
        &mut cc,
    )
    .expect("bisect x² = 6");
    assert!(close(&r, &six.sqrt(4 * P, RM), 125));
    let r = illinois(
        |x, p, rm, _| x.mul(x, p, rm).sub(&ExactNum::from_u8(6, p), p, rm),
        &n(2),
        &n(3),
        &tol,
        P,
        RM,
        &mut cc,
    )
    .expect("illinois x² = 6");
    assert!(close(&r, &six.sqrt(4 * P, RM), 125));
    let r = brent(
        |x, p, rm, _| x.mul(x, p, rm).sub(&ExactNum::from_u8(6, p), p, rm),
        &n(2),
        &n(3),
        &tol,
        P,
        RM,
        &mut cc,
    )
    .expect("brent x² = 6");
    assert!(close(&r, &six.sqrt(4 * P, RM), 125));
}

#[test]
fn ode_rtol_and_equal_exponent_span() {
    let mut cc = Consts::new().unwrap();
    let e = n(1).exp(4 * P, RM, &mut cc);
    let tol = ExactNum::parse("1e-12", Radix::Dec, P, RM, &mut cc);
    // Non-zero rtol (1.0.4: None).
    let (_, y) = rk45_adaptive(
        |_, y, _, _, _| y.clone(),
        &ExactNum::new(P),
        &n(1),
        &n(1),
        &tol,
        &tol,
        P,
        RM,
        &mut cc,
    )
    .expect("rk45 rtol");
    assert!(close(y.get(y.len() - 1).unwrap(), &e, 33));
    // Span [2, 3] (1.0.4: rejected as not t0 < t1): y' = y, y(2) = 1 → y(3) = e.
    let (_, y) = rk4(
        |_, y, _, _, _| y.clone(),
        &n(2),
        &n(1),
        &n(3),
        1000,
        P,
        RM,
        &mut cc,
    )
    .expect("rk4 [2,3]");
    assert!(close(y.get(y.len() - 1).unwrap(), &e, 40));
}

#[test]
fn chebyshev_equal_exponent_interval() {
    let mut cc = Consts::new().unwrap();
    let c = chebyshev_coeffs(
        |x, p, rm, cc| x.ln(p, rm, cc),
        60,
        &n(2),
        &n(3),
        P,
        RM,
        &mut cc,
    )
    .expect("chebyshev [2,3]");
    let x = ExactNum::parse("2.5", Radix::Dec, P, RM, &mut cc);
    let v = chebyshev_eval(&c, &x, &n(2), &n(3), P, RM);
    assert!(close(&v, &x.ln(4 * P, RM, &mut cc), 120));
}
