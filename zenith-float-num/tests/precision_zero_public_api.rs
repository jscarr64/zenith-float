//! 2.0: every public entry point that takes a precision rejects `p = 0` (`MIN_PRECISION` is 1).
//! Generated from the `-> Self` public methods of `ExactNum` (ext.rs, catalog.rs) and `ExactComplex`
//! (complex.rs) that take `p: usize`; Ball/ComplexBall are out of scope.

use zenith_float_num::{Consts, Error, ExactComplex, ExactNum, RoundingMode};

fn bad(v: &ExactNum) -> bool {
    v.err() == Some(Error::InvalidArgument)
}

fn bad_c(v: &ExactComplex) -> bool {
    bad(v.re()) || bad(v.im())
}

#[test]
fn exact_num_entry_points_reject_precision_zero() -> Result<(), Error> {
    let mut cc = Consts::new()?;
    let rm = RoundingMode::ToEven;
    let x =
        || ExactNum::from_word(3, 64).div(&ExactNum::from_word(7, 64), 64, RoundingMode::ToEven);
    let y =
        || ExactNum::from_word(5, 64).div(&ExactNum::from_word(4, 64), 64, RoundingMode::ToEven);
    let mut accepted = Vec::new();
    if !bad(&ExactNum::new(0)) {
        accepted.push("new");
    }
    if !bad(&x().add(&y(), 0, rm)) {
        accepted.push("add");
    }
    if !bad(&ExactNum::add_ref(&y(), &y(), 0, rm)) {
        accepted.push("add_ref");
    }
    if !bad(&ExactNum::sub_ref(&y(), &y(), 0, rm)) {
        accepted.push("sub_ref");
    }
    if !bad(&x().sub(&y(), 0, rm)) {
        accepted.push("sub");
    }
    if !bad(&x().mul(&y(), 0, rm)) {
        accepted.push("mul");
    }
    if !bad(&ExactNum::mul_ref(&y(), &y(), 0, rm)) {
        accepted.push("mul_ref");
    }
    if !bad(&x().sqr(0, rm)) {
        accepted.push("sqr");
    }
    if !bad(&ExactNum::div_ref(&y(), &y(), 0, rm)) {
        accepted.push("div_ref");
    }
    if !bad(&x().fma(&y(), &y(), 0, rm)) {
        accepted.push("fma");
    }
    if !bad(&x().mul_add(&y(), &y(), 0, rm)) {
        accepted.push("mul_add");
    }
    if !bad(&x().div(&y(), 0, rm)) {
        accepted.push("div");
    }
    if !bad(&x().pow(&y(), 0, rm, &mut cc)) {
        accepted.push("pow");
    }
    if !bad(&x().powi(3, 0, rm)) {
        accepted.push("powi");
    }
    if !bad(&x().log(&y(), 0, rm, &mut cc)) {
        accepted.push("log");
    }
    if !bad(&x().atan(0, rm, &mut cc)) {
        accepted.push("atan");
    }
    if !bad(&x().atan2(&y(), 0, rm, &mut cc)) {
        accepted.push("atan2");
    }
    if !bad(&x().hypot(&y(), 0, rm)) {
        accepted.push("hypot");
    }
    if !bad(&x().log1p(0, rm, &mut cc)) {
        accepted.push("log1p");
    }
    if !bad(&x().expm1(0, rm, &mut cc)) {
        accepted.push("expm1");
    }
    if !bad(&x().tanh(0, rm, &mut cc)) {
        accepted.push("tanh");
    }
    if !bad(&ExactNum::max_value(0)) {
        accepted.push("max_value");
    }
    if !bad(&ExactNum::min_value(0)) {
        accepted.push("min_value");
    }
    if !bad(&ExactNum::min_positive(0)) {
        accepted.push("min_positive");
    }
    if !bad(&ExactNum::min_positive_normal(0)) {
        accepted.push("min_positive_normal");
    }
    if !bad(&ExactNum::from_word(5, 0)) {
        accepted.push("from_word");
    }
    if !bad(&x().reciprocal(0, rm)) {
        accepted.push("reciprocal");
    }
    if !bad(&x().ldexp(2, 0, rm)) {
        accepted.push("ldexp");
    }
    if !bad(&x().scalb(2, 0, rm)) {
        accepted.push("scalb");
    }
    if !bad(&x().logb(0, rm)) {
        accepted.push("logb");
    }
    if !bad(&x().copysign(&y(), 0, rm)) {
        accepted.push("copysign");
    }
    if !bad(&x().next_after(&y(), 0, rm)) {
        accepted.push("next_after");
    }
    if !bad(&x().sqrt(0, rm)) {
        accepted.push("sqrt");
    }
    if !bad(&x().nth_root(3, 0, rm)) {
        accepted.push("nth_root");
    }
    if !bad(&x().gammainc(&y(), 0, rm, &mut cc)) {
        accepted.push("gammainc");
    }
    if !bad(&x().gammainc_upper(&y(), 0, rm, &mut cc)) {
        accepted.push("gammainc_upper");
    }
    if !bad(&x().si(0, rm, &mut cc)) {
        accepted.push("si");
    }
    if !bad(&x().fresnel_s(0, rm, &mut cc)) {
        accepted.push("fresnel_s");
    }
    if !bad(&x().fresnel_c(0, rm, &mut cc)) {
        accepted.push("fresnel_c");
    }
    if !bad(&x().ai(0, rm, &mut cc)) {
        accepted.push("ai");
    }
    if !bad(&x().bi(0, rm, &mut cc)) {
        accepted.push("bi");
    }
    if !bad(&x().ai_prime(0, rm, &mut cc)) {
        accepted.push("ai_prime");
    }
    if !bad(&x().bi_prime(0, rm, &mut cc)) {
        accepted.push("bi_prime");
    }
    if !bad(&x().bessel_j(3, 0, rm, &mut cc)) {
        accepted.push("bessel_j");
    }
    if !bad(&x().bessel_j_nu(&y(), 0, rm, &mut cc)) {
        accepted.push("bessel_j_nu");
    }
    if !bad(&x().bessel_y(&y(), 0, rm, &mut cc)) {
        accepted.push("bessel_y");
    }
    if !bad(&x().bessel_i(&y(), 0, rm, &mut cc)) {
        accepted.push("bessel_i");
    }
    if !bad(&x().bessel_k(&y(), 0, rm, &mut cc)) {
        accepted.push("bessel_k");
    }
    if !bad(&x().elliptic_f(&y(), 0, rm, &mut cc)) {
        accepted.push("elliptic_f");
    }
    if !bad(&x().elliptic_e(&y(), 0, rm, &mut cc)) {
        accepted.push("elliptic_e");
    }
    if !bad(&x().elliptic_f_phi(&y(), 0, rm, &mut cc)) {
        accepted.push("elliptic_f_phi");
    }
    if !bad(&x().elliptic_e_phi(&y(), 0, rm, &mut cc)) {
        accepted.push("elliptic_e_phi");
    }
    if !bad(&x().jacobi_am(&y(), 0, rm, &mut cc)) {
        accepted.push("jacobi_am");
    }
    if !bad(&x().jacobi_sn(&y(), 0, rm, &mut cc)) {
        accepted.push("jacobi_sn");
    }
    if !bad(&x().jacobi_cn(&y(), 0, rm, &mut cc)) {
        accepted.push("jacobi_cn");
    }
    if !bad(&x().jacobi_dn(&y(), 0, rm, &mut cc)) {
        accepted.push("jacobi_dn");
    }
    if !bad(&x().jacobi_cd(&y(), 0, rm, &mut cc)) {
        accepted.push("jacobi_cd");
    }
    if !bad(&x().jacobi_ns(&y(), 0, rm, &mut cc)) {
        accepted.push("jacobi_ns");
    }
    if !bad(&x().jacobi_nc(&y(), 0, rm, &mut cc)) {
        accepted.push("jacobi_nc");
    }
    if !bad(&x().jacobi_nd(&y(), 0, rm, &mut cc)) {
        accepted.push("jacobi_nd");
    }
    if !bad(&x().jacobi_sc(&y(), 0, rm, &mut cc)) {
        accepted.push("jacobi_sc");
    }
    if !bad(&x().jacobi_sd(&y(), 0, rm, &mut cc)) {
        accepted.push("jacobi_sd");
    }
    if !bad(&x().jacobi_cs(&y(), 0, rm, &mut cc)) {
        accepted.push("jacobi_cs");
    }
    if !bad(&x().jacobi_ds(&y(), 0, rm, &mut cc)) {
        accepted.push("jacobi_ds");
    }
    if !bad(&x().jacobi_dc(&y(), 0, rm, &mut cc)) {
        accepted.push("jacobi_dc");
    }
    if !bad(&x().betainc(&y(), &y(), 0, rm, &mut cc)) {
        accepted.push("betainc");
    }
    if !bad(&x().scorer_gi(0, rm, &mut cc)) {
        accepted.push("scorer_gi");
    }
    if !bad(&x().scorer_hi(0, rm, &mut cc)) {
        accepted.push("scorer_hi");
    }
    if !bad(&x().kelvin_ber_nu(&y(), 0, rm, &mut cc)) {
        accepted.push("kelvin_ber_nu");
    }
    if !bad(&x().kelvin_bei_nu(&y(), 0, rm, &mut cc)) {
        accepted.push("kelvin_bei_nu");
    }
    if !bad(&x().kelvin_ker_nu(&y(), 0, rm, &mut cc)) {
        accepted.push("kelvin_ker_nu");
    }
    if !bad(&x().kelvin_kei_nu(&y(), 0, rm, &mut cc)) {
        accepted.push("kelvin_kei_nu");
    }
    if !bad(&x().kelvin_ber(0, rm, &mut cc)) {
        accepted.push("kelvin_ber");
    }
    if !bad(&x().kelvin_bei(0, rm, &mut cc)) {
        accepted.push("kelvin_bei");
    }
    if !bad(&x().kelvin_ker(0, rm, &mut cc)) {
        accepted.push("kelvin_ker");
    }
    if !bad(&x().kelvin_kei(0, rm, &mut cc)) {
        accepted.push("kelvin_kei");
    }
    if !bad(&x().struve_h(&y(), 0, rm, &mut cc)) {
        accepted.push("struve_h");
    }
    if !bad(&x().anger_j(&y(), 0, rm, &mut cc)) {
        accepted.push("anger_j");
    }
    if !bad(&x().weber_e(&y(), 0, rm, &mut cc)) {
        accepted.push("weber_e");
    }
    if !bad(&x().clausen_cl2(0, rm, &mut cc)) {
        accepted.push("clausen_cl2");
    }
    if !bad(&x().clausen_cl3(0, rm, &mut cc)) {
        accepted.push("clausen_cl3");
    }
    if !bad(&x().hurwitz_zeta(&y(), 0, rm, &mut cc)) {
        accepted.push("hurwitz_zeta");
    }
    if !bad(&x().riemann_zeta(0, rm, &mut cc)) {
        accepted.push("riemann_zeta");
    }
    if !bad(&x().polygamma(3, 0, rm, &mut cc)) {
        accepted.push("polygamma");
    }
    if !bad(&x().barnes_g(0, rm, &mut cc)) {
        accepted.push("barnes_g");
    }
    if !bad(&x().jacobi_arcsn(&y(), 0, rm, &mut cc)) {
        accepted.push("jacobi_arcsn");
    }
    if !bad(&x().jacobi_arccn(&y(), 0, rm, &mut cc)) {
        accepted.push("jacobi_arccn");
    }
    if !bad(&x().jacobi_arcdn(&y(), 0, rm, &mut cc)) {
        accepted.push("jacobi_arcdn");
    }
    if !bad(&x().hypergeom_0f1(&y(), 0, rm, &mut cc)) {
        accepted.push("hypergeom_0f1");
    }
    if !bad(&x().lambert_w0(0, rm, &mut cc)) {
        accepted.push("lambert_w0");
    }
    if !bad(&x().lambert_wm1(0, rm, &mut cc)) {
        accepted.push("lambert_wm1");
    }
    if !bad(&x().polylog(3, 0, rm, &mut cc)) {
        accepted.push("polylog");
    }
    assert!(accepted.is_empty(), "accepted p = 0: {accepted:?}");
    Ok(())
}

#[test]
fn exact_complex_entry_points_reject_precision_zero() -> Result<(), Error> {
    let mut cc = Consts::new()?;
    let rm = RoundingMode::ToEven;
    let x = || {
        ExactComplex::new(
            ExactNum::from_word(3, 64).div(&ExactNum::from_word(7, 64), 64, RoundingMode::ToEven),
            ExactNum::from_word(1, 64).div(&ExactNum::from_word(5, 64), 64, RoundingMode::ToEven),
        )
    };
    let y = || {
        ExactComplex::new(
            ExactNum::from_word(5, 64).div(&ExactNum::from_word(4, 64), 64, RoundingMode::ToEven),
            ExactNum::from_word(1, 64).div(&ExactNum::from_word(3, 64), 64, RoundingMode::ToEven),
        )
    };
    let mut accepted = Vec::new();
    if !bad_c(&ExactComplex::zero(0)) {
        accepted.push("zero");
    }
    if !bad_c(&ExactComplex::one(0)) {
        accepted.push("one");
    }
    if !bad_c(&ExactComplex::i(0)) {
        accepted.push("i");
    }
    if !bad_c(&x().reciprocal(0, rm)) {
        accepted.push("reciprocal");
    }
    if !bad_c(&ExactComplex::add_ref(&y(), &y(), 0, rm)) {
        accepted.push("add_ref");
    }
    if !bad_c(&ExactComplex::sub_ref(&y(), &y(), 0, rm)) {
        accepted.push("sub_ref");
    }
    if !bad_c(&x().add(&y(), 0, rm)) {
        accepted.push("add");
    }
    if !bad_c(&x().sub(&y(), 0, rm)) {
        accepted.push("sub");
    }
    if !bad_c(&x().mul(&y(), 0, rm)) {
        accepted.push("mul");
    }
    if !bad_c(&x().div(&y(), 0, rm)) {
        accepted.push("div");
    }
    if !bad_c(&x().exp(0, rm, &mut cc)) {
        accepted.push("exp");
    }
    if !bad_c(&x().ln(0, rm, &mut cc)) {
        accepted.push("ln");
    }
    if !bad_c(&x().sin(0, rm, &mut cc)) {
        accepted.push("sin");
    }
    if !bad_c(&x().cos(0, rm, &mut cc)) {
        accepted.push("cos");
    }
    if !bad_c(&x().tan(0, rm, &mut cc)) {
        accepted.push("tan");
    }
    if !bad_c(&x().sinh(0, rm, &mut cc)) {
        accepted.push("sinh");
    }
    if !bad_c(&x().cosh(0, rm, &mut cc)) {
        accepted.push("cosh");
    }
    if !bad_c(&x().tanh(0, rm, &mut cc)) {
        accepted.push("tanh");
    }
    if !bad_c(&x().sqrt(0, rm, &mut cc)) {
        accepted.push("sqrt");
    }
    if !bad_c(&x().log2(0, rm, &mut cc)) {
        accepted.push("log2");
    }
    if !bad_c(&x().log10(0, rm, &mut cc)) {
        accepted.push("log10");
    }
    if !bad_c(&x().log(&y(), 0, rm, &mut cc)) {
        accepted.push("log");
    }
    if !bad_c(&x().log1p(0, rm, &mut cc)) {
        accepted.push("log1p");
    }
    if !bad_c(&x().exp2(0, rm, &mut cc)) {
        accepted.push("exp2");
    }
    if !bad_c(&x().exp10(0, rm, &mut cc)) {
        accepted.push("exp10");
    }
    if !bad_c(&x().expm1(0, rm, &mut cc)) {
        accepted.push("expm1");
    }
    if !bad_c(&x().ldexp(2, 0, rm)) {
        accepted.push("ldexp");
    }
    if !bad_c(&x().scalb(2, 0, rm)) {
        accepted.push("scalb");
    }
    if !bad_c(&x().nth_root(3, 0, rm, &mut cc)) {
        accepted.push("nth_root");
    }
    if !bad_c(&x().cbrt(0, rm, &mut cc)) {
        accepted.push("cbrt");
    }
    if !bad_c(&x().hypot(&y(), 0, rm, &mut cc)) {
        accepted.push("hypot");
    }
    if !bad_c(&x().fma(&y(), &y(), 0, rm)) {
        accepted.push("fma");
    }
    if !bad_c(&x().mul_add(&y(), &y(), 0, rm)) {
        accepted.push("mul_add");
    }
    if !bad_c(&x().pow(&y(), 0, rm, &mut cc)) {
        accepted.push("pow");
    }
    if !bad_c(&x().asin(0, rm, &mut cc)) {
        accepted.push("asin");
    }
    if !bad_c(&x().acos(0, rm, &mut cc)) {
        accepted.push("acos");
    }
    if !bad_c(&x().atan(0, rm, &mut cc)) {
        accepted.push("atan");
    }
    if !bad_c(&x().asinh(0, rm, &mut cc)) {
        accepted.push("asinh");
    }
    if !bad_c(&x().acosh(0, rm, &mut cc)) {
        accepted.push("acosh");
    }
    if !bad_c(&x().atanh(0, rm, &mut cc)) {
        accepted.push("atanh");
    }
    assert!(accepted.is_empty(), "accepted p = 0: {accepted:?}");
    Ok(())
}
