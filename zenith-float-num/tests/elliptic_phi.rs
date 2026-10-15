//! Angle-form elliptic integrals: sine-form consistency, period π, and oddness.
//!
//! `elliptic_f` takes `sin φ`. `elliptic_f_phi` takes `φ`. The same split holds
//! for `E` and `Π`.

use std::ops::Neg;
use zenith_float_num::{Consts, Error, ExactNum, Radix, RoundingMode, INF_POS};

/// `-x` for an operand that stays borrowed: fresh `-1` times `x` at `x`'s own precision.
fn neg_of(x: &ExactNum) -> ExactNum {
    let p = x
        .mantissa_max_bit_len()
        .unwrap_or(zenith_float_num::WORD_BIT_SIZE);
    ExactNum::from_i8(-1, p).mul(x, p, RoundingMode::None)
}

const HALF_PI: &str = "1.57079632679489661923132169163975144209858469968755291048747229615390820314310449931401741267105853399107404325664115332354692230477529111586267970406424055872514205135096926055278";
const M_CUBIC: &str = "0.933012701892219323381861585376468091735701313452595157013951744862983254227200009270286546689312143918906535353851675757492486273749738119702913878023593412132023307975576395516994";
const PHI0: &str = "1.29953257306887582080339773503650160302735297255117556883708654212600457346289335365431960644782250061318071518211065710940897438485509079286935599598508141086119376863166755133786";
const PHI1: &str = "1.64265492275483293872763014546957250537481163728018102972221743087594415906527765906601436244541554006160551657304897484811247961936506448745328548521739496845797918755628708095511";

fn cc() -> Consts {
    match Consts::new() {
        Ok(c) => c,
        Err(_) => panic!("constants"),
    }
}

fn parse(text: &str, p: usize, rm: RoundingMode, consts: &mut Consts) -> ExactNum {
    ExactNum::parse(text, Radix::Dec, p, rm, consts)
}

fn same_bits(a: &ExactNum, b: &ExactNum) -> bool {
    a.sign() == b.sign()
        && a.exponent() == b.exponent()
        && a.mantissa_digits() == b.mantissa_digits()
}

fn is_invalid(n: &ExactNum) -> bool {
    n.is_nan() && n.err() == Some(Error::InvalidArgument)
}

/// Ulps of `a - b` at the stored width of `p`. Zero when the values match.
fn ulp_gap(a: &ExactNum, b: &ExactNum, p: usize) -> i32 {
    let width = p.div_ceil(64) * 64;
    let diff = ExactNum::sub_ref(a, b, width, RoundingMode::None);
    if diff.is_zero() {
        return 0;
    }
    let scale = a.exponent().unwrap_or(0).max(b.exponent().unwrap_or(0));
    let gap = scale - diff.exponent().unwrap_or(i32::MIN);
    width as i32 - gap
}

fn round_down(mut n: ExactNum, p: usize, rm: RoundingMode) -> ExactNum {
    match n.set_precision(p, rm) {
        Ok(()) => n,
        Err(err) => panic!("set_precision {err:?}"),
    }
}

#[test]
fn sine_form_matches_through_asin() {
    let mut consts = cc();
    let p = 64;
    let rm = RoundingMode::ToEven;
    let m = parse("0.5", p, rm, &mut consts);
    for text in ["0", "0.5", "-0.5", "1", "-1", "0.3"] {
        let x = parse(text, p, rm, &mut consts);
        let angle = x.asin(p, rm, &mut consts);
        let by_angle = angle.elliptic_f_phi(&m, p, rm, &mut consts);
        let by_sine = x.elliptic_f(&m, p, rm, &mut consts);
        let gap = ulp_gap(&by_angle, &by_sine, p);
        assert!(
            gap <= 1,
            "f_phi(asin({text})) and elliptic_f({text}) differ by {gap} ulps"
        );
        if text == "0" || text == "1" || text == "-1" {
            assert!(
                same_bits(&by_angle, &by_sine),
                "endpoint {text} must be limb-identical"
            );
        }
    }
}

#[test]
fn angle_two_is_not_the_sine_form() {
    let mut consts = cc();
    let p = 64;
    let rm = RoundingMode::ToEven;
    let two = parse("2", p, rm, &mut consts);
    let m = parse("0.5", p, rm, &mut consts);
    assert!(is_invalid(&two.elliptic_f(&m, p, rm, &mut consts)));
    let angle = two.elliptic_f_phi(&m, p, rm, &mut consts);
    assert!(!angle.is_nan() && !angle.is_inf());
    let via_sin = two
        .sin(p, rm, &mut consts)
        .elliptic_f(&m, p, rm, &mut consts);
    assert!(
        !same_bits(&angle, &via_sin),
        "F(2|m) collapsed to F(arcsin(sin 2)|m)"
    );
}

#[test]
fn half_pi_matches_complete() {
    let mut consts = cc();
    let p = 64;
    let rm = RoundingMode::ToEven;
    let phi = parse(HALF_PI, p, rm, &mut consts);
    let m = parse("0.5", p, rm, &mut consts);
    let n = parse("0.2", p, rm, &mut consts);
    assert!(same_bits(
        &phi.elliptic_f_phi(&m, p, rm, &mut consts),
        &m.elliptic_k(p, rm, &mut consts)
    ));
    assert!(same_bits(
        &phi.elliptic_e_phi(&m, p, rm, &mut consts),
        &m.elliptic_e_complete(p, rm, &mut consts)
    ));
    assert!(same_bits(
        &n.elliptic_pi_phi(&phi, &m, p, rm, &mut consts),
        &n.elliptic_pi_complete(&m, p, rm, &mut consts)
    ));
}

#[test]
fn periodicity_adds_two_complete_periods() {
    let mut consts = cc();
    let p = 64;
    let hi = 192;
    let rm = RoundingMode::ToEven;
    let none = RoundingMode::None;
    let m_text = "0.5";
    let phi_text = "0.4";
    let n_text = "0.2";
    let phi = parse(phi_text, hi, rm, &mut consts);
    let m = parse(m_text, hi, rm, &mut consts);
    let n = parse(n_text, hi, rm, &mut consts);
    let pi = consts.pi(hi, rm);
    let two = parse("2", hi, rm, &mut consts);
    let shifted = pi.add(&phi, hi, none);

    let f_shift = shifted.elliptic_f_phi(&m, p, rm, &mut consts);
    let f_sum = phi.elliptic_f_phi(&m, hi, rm, &mut consts).add(
        &two.mul(&m.elliptic_k(hi, rm, &mut consts), hi, none),
        hi,
        none,
    );
    assert!(same_bits(&f_shift, &round_down(f_sum, p, rm)));

    let e_shift = shifted.elliptic_e_phi(&m, p, rm, &mut consts);
    let e_sum = phi.elliptic_e_phi(&m, hi, rm, &mut consts).add(
        &two.mul(&m.elliptic_e_complete(hi, rm, &mut consts), hi, none),
        hi,
        none,
    );
    assert!(same_bits(&e_shift, &round_down(e_sum, p, rm)));

    let p_shift = n.elliptic_pi_phi(&shifted, &m, p, rm, &mut consts);
    let p_sum = n.elliptic_pi_phi(&phi, &m, hi, rm, &mut consts).add(
        &two.mul(&n.elliptic_pi_complete(&m, hi, rm, &mut consts), hi, none),
        hi,
        none,
    );
    assert!(same_bits(&p_shift, &round_down(p_sum, p, rm)));
}

#[test]
fn odd_in_phi_under_to_even() {
    let mut consts = cc();
    let p = 64;
    let rm = RoundingMode::ToEven;
    let m = parse("0.3", p, rm, &mut consts);
    let n = parse("0.2", p, rm, &mut consts);
    for text in ["0.4", "2", "5"] {
        let phi = parse(text, p, rm, &mut consts);
        let neg = neg_of(&phi);
        let f = phi.elliptic_f_phi(&m, p, rm, &mut consts);
        let fn_ = neg.elliptic_f_phi(&m, p, rm, &mut consts);
        assert!(same_bits(&f.neg(), &fn_), "F odd at {text}");
        let e = phi.elliptic_e_phi(&m, p, rm, &mut consts);
        let en = neg.elliptic_e_phi(&m, p, rm, &mut consts);
        assert!(same_bits(&e.neg(), &en), "E odd at {text}");
        let p_pos = n.elliptic_pi_phi(&phi, &m, p, rm, &mut consts);
        let p_neg = n.elliptic_pi_phi(&neg, &m, p, rm, &mut consts);
        assert!(same_bits(&p_pos.neg(), &p_neg), "Pi odd at {text}");
    }
}

#[test]
fn infinity_is_invalid_argument() {
    let mut consts = cc();
    let p = 64;
    let rm = RoundingMode::ToEven;
    let m = parse("0.5", p, rm, &mut consts);
    let n = parse("0.2", p, rm, &mut consts);
    let phi = parse("0.4", p, rm, &mut consts);
    assert!(is_invalid(&INF_POS.elliptic_f_phi(&m, p, rm, &mut consts)));
    assert!(is_invalid(&INF_POS.elliptic_e_phi(&m, p, rm, &mut consts)));
    assert!(is_invalid(&n.elliptic_pi_phi(
        &INF_POS,
        &m,
        p,
        rm,
        &mut consts
    )));
    assert!(is_invalid(&INF_POS.elliptic_pi_phi(
        &phi,
        &m,
        p,
        rm,
        &mut consts
    )));
    assert!(is_invalid(&phi.elliptic_f_phi(
        &INF_POS,
        p,
        rm,
        &mut consts
    )));
}

#[test]
fn nan_keeps_its_tag() {
    let mut consts = cc();
    let p = 64;
    let rm = RoundingMode::ToEven;
    let tagged = ExactNum::nan(Some(Error::DivisionByZero));
    let m = parse("0.5", p, rm, &mut consts);
    let y = tagged.elliptic_f_phi(&m, p, rm, &mut consts);
    assert!(y.is_nan());
    assert_eq!(y.err(), Some(Error::DivisionByZero));
}

#[test]
fn cubic_integral_identity() {
    // ∫_0^1 dx/√(x³+1) = 3^{-1/4} (F(φ1|m) - F(φ0|m)), m = (2+√3)/4.
    // φ1 > π/2, so the sine-form API cannot express F(φ1|m).
    // The printed 20 digits stop at a truncation whose residual is about 2.5e-21.
    let mut consts = cc();
    let p = 128;
    let rm = RoundingMode::ToEven;
    let m = parse(M_CUBIC, p, rm, &mut consts);
    let phi0 = parse(PHI0, p, rm, &mut consts);
    let phi1 = parse(PHI1, p, rm, &mut consts);
    let half = consts.pi(p, rm).div(&parse("2", p, rm, &mut consts), p, rm);
    assert_eq!(phi1.compare(&half), Some(1));
    let diff = phi1.elliptic_f_phi(&m, p, rm, &mut consts).sub(
        &phi0.elliptic_f_phi(&m, p, rm, &mut consts),
        p,
        rm,
    );
    let fourth = parse("3", p, rm, &mut consts).sqrt(p, rm).sqrt(p, rm);
    let scale = parse("1", p, rm, &mut consts).div(&fourth, p, rm);
    let got = scale.mul(&diff, p, rm);
    let printed = parse("0.90960424263889577391", p, rm, &mut consts);
    let tol = parse("5e-21", p, rm, &mut consts);
    let err = got.sub(&printed, p, RoundingMode::None).abs();
    assert!(
        matches!(err.compare(&tol), Some(c) if c <= 0),
        "cubic integral missed 0.90960424263889577391"
    );
}
