// Additional tests of the library.

use zenith_float::{cexpr, exact, fbig, ExactComplex, RadixFloat};
use zenith_float_macro::expr;
use zenith_float_num::{
    ctx::Context, Consts, ExactNum, Radix, RoundingMode, Sign, EXPONENT_MAX, EXPONENT_MIN,
    WORD_BIT_SIZE, WORD_MAX, WORD_SIGNIFICANT_BIT,
};

#[test]
fn native_operators() {
    let a = ExactNum::from(6);
    let b = ExactNum::from(2);
    assert_eq!((&a + &b).cmp(&ExactNum::from(8)), Some(0));
    assert_eq!((&a - &b).cmp(&ExactNum::from(4)), Some(0));
    assert_eq!((&a * &b).cmp(&ExactNum::from(12)), Some(0));
    assert_eq!((&a / &b).cmp(&ExactNum::from(3)), Some(0));
}

#[test]
fn compile_time_literals() {
    let a = exact!("2.5");
    let b = fbig!("2.5");
    assert_eq!(a.cmp(&b), Some(0));
}

fn cplx_eq(name: &str, a: &ExactComplex, b: &ExactComplex) {
    if a.is_nan() && b.is_nan() {
        return;
    }
    assert_eq!(
        a.re().cmp(b.re()),
        Some(0),
        "{name} re {:?} vs {:?}",
        a.re(),
        b.re()
    );
    assert_eq!(
        a.im().cmp(b.im()),
        Some(0),
        "{name} im {:?} vs {:?}",
        a.im(),
        b.im()
    );
}

fn tiny_part(x: &ExactNum, p: usize) -> bool {
    x.is_zero() || x.exponent().unwrap_or(0) < -((p as i32) / 4)
}

#[test]
fn cexpr_i_squared() {
    let mut ctx = Context::new(
        256,
        RoundingMode::ToEven,
        Consts::new().unwrap(),
        -10000,
        10000,
    );
    let z = cexpr!(I * I, &mut ctx);
    assert_eq!(z.re().cmp(&ExactNum::from_i8(-1, 256)), Some(0));
    assert!(z.im().is_zero());

    let z0 = ExactComplex::zero(256);
    let s = cexpr!(sin(z0), &mut ctx);
    assert!(s.re().is_zero() || s.re().exponent().unwrap_or(0) < -40);
    assert!(s.im().is_zero() || s.im().exponent().unwrap_or(0) < -40);
}

#[test]
fn macro_run_cexpr_tests() {
    let p = 320;
    let rm = RoundingMode::None;
    let mut cc = Consts::new().unwrap();
    let mut ctx = Context::new(p, rm, Consts::new().unwrap(), -1000000, 1000000);

    let x = ExactComplex::from_real(ExactNum::from(2), p);
    let y = ExactComplex::from_real(ExactNum::from(5), p);
    let i = ExactComplex::i(p);
    let z0 = ExactComplex::zero(p);
    let z1 = ExactComplex::one(p);

    cplx_eq(
        "neg1",
        &cexpr!(-1, &mut ctx),
        &ExactComplex::from_real(ExactNum::from(-1), p),
    );
    cplx_eq(
        "add",
        &cexpr!(2 + 3, &mut ctx),
        &ExactComplex::from_real(ExactNum::from(5), p),
    );
    cplx_eq(
        "sub",
        &cexpr!(3 - 4, &mut ctx),
        &ExactComplex::from_real(ExactNum::from(-1), p),
    );
    cplx_eq(
        "mul",
        &cexpr!(4 * 5, &mut ctx),
        &ExactComplex::from_real(ExactNum::from(20), p),
    );
    cplx_eq(
        "div",
        &cexpr!(5 / 6, &mut ctx),
        &ExactComplex::from_real(ExactNum::from(5), p).div(
            &ExactComplex::from_real(ExactNum::from(6), p),
            p,
            rm,
        ),
    );
    cplx_eq("i2", &cexpr!(I * I, &mut ctx), &i.mul(&i, p, rm));

    cplx_eq("recip", &cexpr!(recip(x), &mut ctx), &x.reciprocal(p, rm));
    cplx_eq("sqrt", &cexpr!(sqrt(x), &mut ctx), &x.sqrt(p, rm, &mut cc));
    cplx_eq("ln", &cexpr!(ln(x), &mut ctx), &x.ln(p, rm, &mut cc));
    cplx_eq("log2", &cexpr!(log2(x), &mut ctx), &x.log2(p, rm, &mut cc));
    cplx_eq(
        "log10",
        &cexpr!(log10(x), &mut ctx),
        &x.log10(p, rm, &mut cc),
    );
    cplx_eq(
        "log",
        &cexpr!(log(x, y), &mut ctx),
        &x.log(&y, p, rm, &mut cc),
    );
    cplx_eq(
        "log1p",
        &cexpr!(log1p(z1), &mut ctx),
        &z1.log1p(p, rm, &mut cc),
    );
    cplx_eq("exp", &cexpr!(exp(x), &mut ctx), &x.exp(p, rm, &mut cc));
    cplx_eq("exp2", &cexpr!(exp2(x), &mut ctx), &x.exp2(p, rm, &mut cc));
    cplx_eq(
        "exp10",
        &cexpr!(exp10(x), &mut ctx),
        &x.exp10(p, rm, &mut cc),
    );
    cplx_eq(
        "expm1",
        &cexpr!(expm1(x), &mut ctx),
        &x.expm1(p, rm, &mut cc),
    );
    cplx_eq(
        "pow",
        &cexpr!(pow(x, y), &mut ctx),
        &x.pow(&y, p, rm, &mut cc),
    );

    cplx_eq("sin", &cexpr!(sin(x), &mut ctx), &x.sin(p, rm, &mut cc));
    cplx_eq("cos", &cexpr!(cos(x), &mut ctx), &x.cos(p, rm, &mut cc));
    cplx_eq("tan0", &cexpr!(tan(z0), &mut ctx), &z0.tan(p, rm, &mut cc));
    let u = ExactComplex::from_real(ExactNum::from(1), p);
    let half = ExactComplex::from_real(ExactNum::from(1).div(&ExactNum::from(2), p, rm), p);
    cplx_eq("asin", &cexpr!(asin(u), &mut ctx), &u.asin(p, rm, &mut cc));
    cplx_eq("acos", &cexpr!(acos(u), &mut ctx), &u.acos(p, rm, &mut cc));
    cplx_eq("atan", &cexpr!(atan(u), &mut ctx), &u.atan(p, rm, &mut cc));
    cplx_eq("sinh", &cexpr!(sinh(x), &mut ctx), &x.sinh(p, rm, &mut cc));
    cplx_eq("cosh", &cexpr!(cosh(x), &mut ctx), &x.cosh(p, rm, &mut cc));
    cplx_eq("tanh", &cexpr!(tanh(x), &mut ctx), &x.tanh(p, rm, &mut cc));
    cplx_eq(
        "asinh",
        &cexpr!(asinh(x), &mut ctx),
        &x.asinh(p, rm, &mut cc),
    );
    cplx_eq(
        "acosh",
        &cexpr!(acosh(x), &mut ctx),
        &x.acosh(p, rm, &mut cc),
    );
    cplx_eq(
        "atanh",
        &cexpr!(atanh(half), &mut ctx),
        &half.atanh(p, rm, &mut cc),
    );

    let a = ExactComplex::new(ExactNum::from_u8(3, p), ExactNum::from_u8(4, p));
    let mag = cexpr!(abs(a), &mut ctx);
    assert_eq!(mag.re().cmp(&a.abs(p, rm)), Some(0));
    assert!(tiny_part(mag.im(), p));
    let th = cexpr!(arg(a), &mut ctx);
    assert_eq!(th.re().cmp(&a.arg(p, rm, &mut cc)), Some(0));
    assert!(tiny_part(th.im(), p));
    cplx_eq("conj", &cexpr!(conj(a), &mut ctx), &a.conj());
    cplx_eq("ldexp", &cexpr!(ldexp(x, 3), &mut ctx), &x.ldexp(3, p, rm));
    cplx_eq("scalb", &cexpr!(scalb(x, 3), &mut ctx), &x.scalb(3, p, rm));
    cplx_eq("cbrt", &cexpr!(cbrt(x), &mut ctx), &x.cbrt(p, rm, &mut cc));
    cplx_eq(
        "root",
        &cexpr!(root(x, 3), &mut ctx),
        &x.nth_root(3, p, rm, &mut cc),
    );
    cplx_eq(
        "hypot",
        &cexpr!(hypot(x, y), &mut ctx),
        &x.hypot(&y, p, rm, &mut cc),
    );
    cplx_eq(
        "fma",
        &cexpr!(fma(x, y, z1), &mut ctx),
        &x.fma(&y, &z1, p, rm),
    );
    let lb = cexpr!(logb(x), &mut ctx);
    assert_eq!(lb.re().cmp(&x.logb(p, rm)), Some(0));
    assert!(tiny_part(lb.im(), p));

    let m1 = ExactComplex::from_real(ExactNum::from_i8(-1, p), p);
    let sqm = cexpr!(sqrt(m1), &mut ctx);
    assert!(tiny_part(sqm.re(), p));
    assert!(sqm.im().is_positive());
    let lnm = cexpr!(ln(m1), &mut ctx);
    assert!(tiny_part(lnm.re(), p));
    assert!(lnm.im().is_positive());

    cplx_eq(
        "pi",
        &cexpr!(pi, &mut ctx),
        &ExactComplex::from_real(ctx.const_pi(), p),
    );
    cplx_eq(
        "e",
        &cexpr!(e, &mut ctx),
        &ExactComplex::from_real(ctx.const_e(), p),
    );
    cplx_eq(
        "ln2",
        &cexpr!(ln_2, &mut ctx),
        &ExactComplex::from_real(ctx.const_ln2(), p),
    );
    cplx_eq(
        "ln10",
        &cexpr!(ln_10, &mut ctx),
        &ExactComplex::from_real(ctx.const_ln10(), p),
    );
    cplx_eq(
        "sqrt2",
        &cexpr!(sqrt2, &mut ctx),
        &ExactComplex::from_real(ctx.const_sqrt2(), p),
    );
    cplx_eq(
        "phi",
        &cexpr!(phi, &mut ctx),
        &ExactComplex::from_real(ctx.const_phi(), p),
    );
    cplx_eq(
        "gamma",
        &cexpr!(euler_gamma, &mut ctx),
        &ExactComplex::from_real(ctx.const_euler_gamma(), p),
    );
    cplx_eq("cerf0", &cexpr!(erf(z0), &mut ctx), &z0.erf(p, rm, &mut cc));
    cplx_eq(
        "cerfc0",
        &cexpr!(erfc(z0), &mut ctx),
        &z0.erfc(p, rm, &mut cc),
    );
    cplx_eq(
        "cgamma1",
        &cexpr!(gamma(z1), &mut ctx),
        &z1.gamma(p, rm, &mut cc),
    );
    cplx_eq(
        "clngamma1",
        &cexpr!(ln_gamma(z1), &mut ctx),
        &z1.ln_gamma(p, rm, &mut cc),
    );
    cplx_eq(
        "cdigamma1",
        &cexpr!(digamma(z1), &mut ctx),
        &z1.digamma(p, rm, &mut cc),
    );
    cplx_eq("cei1", &cexpr!(ei(z1), &mut ctx), &z1.ei(p, rm, &mut cc));
    cplx_eq("csi0", &cexpr!(si(z0), &mut ctx), &z0.si(p, rm, &mut cc));
    cplx_eq("cci1", &cexpr!(ci(z1), &mut ctx), &z1.ci(p, rm, &mut cc));
    cplx_eq("cai0", &cexpr!(ai(z0), &mut ctx), &z0.ai(p, rm, &mut cc));
    cplx_eq("cbi0", &cexpr!(bi(z0), &mut ctx), &z0.bi(p, rm, &mut cc));
    cplx_eq("cli2", &cexpr!(li(x), &mut ctx), &x.li(p, rm, &mut cc));
    cplx_eq(
        "cfs0",
        &cexpr!(fresnel_s(z0), &mut ctx),
        &z0.fresnel_s(p, rm, &mut cc),
    );
    cplx_eq(
        "cfc0",
        &cexpr!(fresnel_c(z0), &mut ctx),
        &z0.fresnel_c(p, rm, &mut cc),
    );
    cplx_eq(
        "cj0",
        &cexpr!(bessel_j_nu(z1, z0), &mut ctx),
        &z1.bessel_j_nu(&z0, p, rm, &mut cc),
    );
    cplx_eq(
        "cy0",
        &cexpr!(bessel_y(z1, z0), &mut ctx),
        &z1.bessel_y(&z0, p, rm, &mut cc),
    );
    cplx_eq(
        "ci0",
        &cexpr!(bessel_i(z1, z0), &mut ctx),
        &z1.bessel_i(&z0, p, rm, &mut cc),
    );
    let ck = cexpr!(bessel_k(z1, z0), &mut ctx);
    let k0 = z1.bessel_k(&z0, p, rm, &mut cc);
    assert_eq!(ck.re().cmp(k0.re()), Some(0), "ck0 re");
    assert!(tiny_part(ck.im(), p) && tiny_part(k0.im(), p));

    let mh = ExactComplex::from_real(
        ExactNum::from(1).div(&ExactNum::from(2), p, rm),
        p,
    );
    cplx_eq(
        "cellk",
        &cexpr!(elliptic_k(mh), &mut ctx),
        &mh.elliptic_k(p, rm, &mut cc),
    );
    cplx_eq(
        "celle",
        &cexpr!(elliptic_e(z0), &mut ctx),
        &z0.elliptic_e_complete(p, rm, &mut cc),
    );
    let cf = cexpr!(elliptic_f(mh, z0), &mut ctx);
    let f0 = mh.elliptic_f(&z0, p, rm, &mut cc);
    assert_eq!(cf.re().cmp(f0.re()), Some(0), "cellf re");
    assert!(tiny_part(cf.im(), p) && tiny_part(f0.im(), p));
    let cei = cexpr!(elliptic_e_inc(mh, z0), &mut ctx);
    let ei0 = mh.elliptic_e(&z0, p, rm, &mut cc);
    assert_eq!(cei.re().cmp(ei0.re()), Some(0), "cellei re");
    assert!(tiny_part(cei.im(), p) && tiny_part(ei0.im(), p));
    cplx_eq(
        "cellpi",
        &cexpr!(elliptic_pi(z0, mh), &mut ctx),
        &z0.elliptic_pi_complete(&mh, p, rm, &mut cc),
    );
    let cpi = cexpr!(elliptic_pi_inc(z0, mh, z0), &mut ctx);
    let pi0 = z0.elliptic_pi(&mh, &z0, p, rm, &mut cc);
    assert_eq!(cpi.re().cmp(pi0.re()), Some(0), "cellpii re");
    assert!(tiny_part(cpi.im(), p) && tiny_part(pi0.im(), p));

    let back = cexpr!(ln(exp(x)), &mut ctx);
    let d = back.re().sub(x.re(), p, RoundingMode::None).abs();
    assert!(d.exponent().unwrap_or(0) < -((p as i32) / 8));
    assert!(tiny_part(back.im(), p));

    let two = ExactComplex::from_real(ExactNum::from(2), p);
    cplx_eq(
        "c2f1",
        &cexpr!(hypergeom_2f1(z1, z1, two, mh), &mut ctx),
        &z1.hypergeom_2f1(&z1, &two, &mh, p, rm, &mut cc),
    );
    cplx_eq(
        "exp2lit",
        &cexpr!(exp2(2), &mut ctx),
        &two.exp2(p, rm, &mut cc),
    );

    let mut tight = Context::new(p, RoundingMode::ToEven, Consts::new().unwrap(), -10, 10);
    let huge = cexpr!(exp(40), &mut tight);
    assert!(huge.re().is_inf() || huge.im().is_inf() || huge.is_nan() || huge.re().is_zero());
}

#[test]
fn cexpr_cancellation_matches_high_prec() {
    let p = 192;
    let rm = RoundingMode::ToEven;
    let mut cc = Consts::new().unwrap();
    let mut ctx = Context::new(p, rm, Consts::new().unwrap(), EXPONENT_MIN, EXPONENT_MAX);

    let x = ExactNum::parse(
        "0.00000000000000000000000000000000000001",
        Radix::Dec,
        p,
        RoundingMode::None,
        &mut cc,
    );
    let y = ExactNum::parse(
        "1.57079632679489661923132169163975144209",
        Radix::Dec,
        p,
        RoundingMode::None,
        &mut cc,
    );
    let zx = ExactComplex::from_real(x.clone(), p);
    let zy = ExactComplex::from_real(y.clone(), p);

    let z = cexpr!(cos(zx) - sin(zy), &mut ctx);
    let cx = zx.cos(p + 256, RoundingMode::None, &mut cc);
    let sy = zy.sin(p + 256, RoundingMode::None, &mut cc);
    let r = cx.sub(&sy, p, rm);
    cplx_eq("cancel", &z, &r);
}

#[test]
fn radix_base12_roundtrip() {
    let mut cc = Consts::new().unwrap();
    let rdx = Radix::try_new(12).unwrap();
    let p = 192;
    let rm = RoundingMode::ToEven;
    let n = ExactNum::parse("10.A", rdx, p, rm, &mut cc);
    let s = n.format(rdx, rm, &mut cc).unwrap();
    let g = ExactNum::parse(&s, rdx, p, rm, &mut cc);
    assert_eq!(n.cmp(&g), Some(0));

    let rf = n.clone().with_radix(rdx);
    assert_eq!(rf.radix(), rdx);
    let s2 = rf.format(rm, &mut cc).unwrap();
    let g2 = RadixFloat::parse(&s2, rdx, p, rm, &mut cc).unwrap();
    assert_eq!(n.cmp(g2.value()), Some(0));
}

#[test]
fn macro_compile_tests() {
    let t = trybuild::TestCases::new();
    t.pass("./tests/tests/expr.rs");
    t.pass("./tests/tests/cexpr.rs");
}

#[test]
fn macro_run_basic_tests() {
    let p = 320;
    let rm = RoundingMode::None;
    let mut cc = Consts::new().unwrap();

    let mut ctx = Context::new(p, rm, Consts::new().unwrap(), -1000000, 1000000);

    let x = ExactNum::from(2);
    let y = ExactNum::from(5); // Change 4.56 to 5

    let res: ExactNum = expr!(-1, &mut ctx);
    debug_assert_eq!(res, ExactNum::from(-1));

    let res: ExactNum = expr!(2 + 3, &mut ctx);
    debug_assert_eq!(res, ExactNum::from(5));

    let res: ExactNum = expr!(3 - 4, &mut ctx);
    debug_assert_eq!(res, ExactNum::from(-1));

    let res: ExactNum = expr!(4 * 5, &mut ctx);
    debug_assert_eq!(res, ExactNum::from(20));

    let res: ExactNum = expr!(5 / 6, &mut ctx);
    debug_assert_eq!(res, ExactNum::from(5).div(&ExactNum::from(6), p, rm));

    let res: ExactNum = expr!(6 % 7, &mut ctx);
    debug_assert_eq!(res, ExactNum::from(6));

    let res: ExactNum = expr!(ln(x), &mut ctx);
    debug_assert_eq!(res, x.ln(p, rm, &mut cc));

    let res: ExactNum = expr!(log2(x), &mut ctx);
    debug_assert_eq!(res, x.log2(p, rm, &mut cc));

    let res: ExactNum = expr!(log10(x), &mut ctx);
    debug_assert_eq!(res, x.log10(p, rm, &mut cc));

    let res: ExactNum = expr!(log(x, y), &mut ctx);
    debug_assert_eq!(res, x.log(&y, p, rm, &mut cc));

    let res: ExactNum = expr!(exp(x), &mut ctx);
    debug_assert_eq!(res, x.exp(p, rm, &mut cc));

    let res: ExactNum = expr!(exp2(x), &mut ctx);
    debug_assert_eq!(res, x.exp2(p, rm, &mut cc));

    let res: ExactNum = expr!(exp10(x), &mut ctx);
    debug_assert_eq!(res, x.exp10(p, rm, &mut cc));

    let res: ExactNum = expr!(rem_pi(x), &mut ctx);
    debug_assert_eq!(res, x.rem_pi(p, rm, &mut cc));

    let res: ExactNum = expr!(pow(x, y), &mut ctx);
    debug_assert_eq!(res, x.pow(&y, p, rm, &mut cc));

    let res: ExactNum = expr!(sin(x), &mut ctx);
    debug_assert_eq!(res, x.sin(p, rm, &mut cc));

    let res: ExactNum = expr!(cos(x), &mut ctx);
    debug_assert_eq!(res, x.cos(p, rm, &mut cc));

    let (ps, pc) = x.sin_cos(p, rm, &mut cc);
    debug_assert_eq!(ps, x.sin(p, rm, &mut cc));
    debug_assert_eq!(pc, x.cos(p, rm, &mut cc));

    let res: ExactNum = expr!(tan(x), &mut ctx);
    debug_assert_eq!(res, x.tan(p, rm, &mut cc));

    let x = ExactNum::from(1);

    let res: ExactNum = expr!(asin(x), &mut ctx);
    debug_assert_eq!(res, x.asin(p, rm, &mut cc));

    let res: ExactNum = expr!(acos(x), &mut ctx);
    debug_assert_eq!(res, x.acos(p, rm, &mut cc));

    let res: ExactNum = expr!(atan(x), &mut ctx);
    debug_assert_eq!(res, x.atan(p, rm, &mut cc));

    let x = ExactNum::from(2);

    let res: ExactNum = expr!(sinh(x), &mut ctx);
    debug_assert_eq!(res, x.sinh(p, rm, &mut cc));

    let res: ExactNum = expr!(cosh(x), &mut ctx);
    debug_assert_eq!(res, x.cosh(p, rm, &mut cc));

    let (psh, pch) = x.sinh_cosh(p, rm, &mut cc);
    debug_assert_eq!(psh, x.sinh(p, rm, &mut cc));
    debug_assert_eq!(pch, x.cosh(p, rm, &mut cc));

    let res: ExactNum = expr!(tanh(x), &mut ctx);
    debug_assert_eq!(res, x.tanh(p, rm, &mut cc));

    let res: ExactNum = expr!(asinh(x), &mut ctx);
    debug_assert_eq!(res, x.asinh(p, rm, &mut cc));

    let res: ExactNum = expr!(acosh(x), &mut ctx);
    debug_assert_eq!(res, x.acosh(p, rm, &mut cc));

    let x = ExactNum::from(1);

    let res: ExactNum = expr!(atanh(x), &mut ctx);
    debug_assert_eq!(res, x.atanh(p, rm, &mut cc));

    let y = ExactNum::from(3);
    let x = ExactNum::from(4);
    let res: ExactNum = expr!(hypot(x, y), &mut ctx);
    debug_assert_eq!(res, x.hypot(&y, p, rm));

    let z = ExactNum::from(2);
    let res: ExactNum = expr!(fma(x, y, z), &mut ctx);
    debug_assert_eq!(res, x.fma(&y, &z, p, rm));

    let res: ExactNum = expr!(mul_add(x, y, z), &mut ctx);
    debug_assert_eq!(res, x.mul_add(&y, &z, p, rm));

    let res: ExactNum = expr!(root(x, 5), &mut ctx);
    debug_assert_eq!(res, x.nth_root(5, p, rm));

    let res: ExactNum = expr!(atan2(y, x), &mut ctx);
    debug_assert_eq!(res, y.atan2(&x, p, rm, &mut cc));

    let x = ExactNum::from(1);
    let res: ExactNum = expr!(log1p(x), &mut ctx);
    debug_assert_eq!(res, x.log1p(p, rm, &mut cc));

    let res: ExactNum = expr!(expm1(x), &mut ctx);
    debug_assert_eq!(res, x.expm1(p, rm, &mut cc));

    let res: ExactNum = expr!(erf(x), &mut ctx);
    debug_assert_eq!(res, x.erf(p, rm, &mut cc));

    let res: ExactNum = expr!(erfc(x), &mut ctx);
    debug_assert_eq!(res, x.erfc(p, rm, &mut cc));

    let res: ExactNum = expr!(gamma(x), &mut ctx);
    debug_assert_eq!(res, x.gamma(p, rm, &mut cc));

    let res: ExactNum = expr!(ln_gamma(x), &mut ctx);
    debug_assert_eq!(res, x.ln_gamma(p, rm, &mut cc));

    let res: ExactNum = expr!(digamma(x), &mut ctx);
    debug_assert_eq!(res, x.digamma(p, rm, &mut cc));

    let res: ExactNum = expr!(gammainc(x, x), &mut ctx);
    debug_assert_eq!(res, x.gammainc(&x, p, rm, &mut cc));

    let res: ExactNum = expr!(bessel_j(x, 0), &mut ctx);
    debug_assert_eq!(res, x.bessel_j(0, p, rm, &mut cc));

    let nu0 = ExactNum::from(0);
    let res: ExactNum = expr!(bessel_j_nu(x, nu0), &mut ctx);
    debug_assert_eq!(res, x.bessel_j_nu(&nu0, p, rm, &mut cc));
    let res: ExactNum = expr!(bessel_i(x, nu0), &mut ctx);
    debug_assert_eq!(res, x.bessel_i(&nu0, p, rm, &mut cc));

    let mh = ExactNum::from(1).div(&ExactNum::from(2), p, rm);
    let res: ExactNum = expr!(elliptic_k(mh), &mut ctx);
    debug_assert_eq!(res, mh.elliptic_k(p, rm, &mut cc));
    let res: ExactNum = expr!(elliptic_e(mh), &mut ctx);
    debug_assert_eq!(res, mh.elliptic_e_complete(p, rm, &mut cc));
    let res: ExactNum = expr!(elliptic_f(mh, mh), &mut ctx);
    debug_assert_eq!(res, mh.elliptic_f(&mh, p, rm, &mut cc));
    let res: ExactNum = expr!(elliptic_e_inc(mh, mh), &mut ctx);
    debug_assert_eq!(res, mh.elliptic_e(&mh, p, rm, &mut cc));
    let res: ExactNum = expr!(elliptic_pi(mh, mh), &mut ctx);
    debug_assert_eq!(res, mh.elliptic_pi_complete(&mh, p, rm, &mut cc));
    let res: ExactNum = expr!(elliptic_pi_inc(mh, mh, mh), &mut ctx);
    debug_assert_eq!(res, mh.elliptic_pi(&mh, &mh, p, rm, &mut cc));
    let res: ExactNum = expr!(jacobi_sn(mh, mh), &mut ctx);
    debug_assert_eq!(res, mh.jacobi_sn(&mh, p, rm, &mut cc));
    let res: ExactNum = expr!(jacobi_cn(mh, mh), &mut ctx);
    debug_assert_eq!(res, mh.jacobi_cn(&mh, p, rm, &mut cc));
    let res: ExactNum = expr!(jacobi_dn(mh, mh), &mut ctx);
    debug_assert_eq!(res, mh.jacobi_dn(&mh, p, rm, &mut cc));
    let res: ExactNum = expr!(jacobi_am(mh, mh), &mut ctx);
    debug_assert_eq!(res, mh.jacobi_am(&mh, p, rm, &mut cc));
    let res: ExactNum = expr!(jacobi_cd(mh, mh), &mut ctx);
    debug_assert_eq!(res, mh.jacobi_cd(&mh, p, rm, &mut cc));

    let res: ExactNum = expr!(legendre_p(x, 2), &mut ctx);
    debug_assert_eq!(res, x.legendre_p(2, p, rm));
    let res: ExactNum = expr!(legendre_p_assoc(mh, 1, 1), &mut ctx);
    debug_assert_eq!(res, mh.assoc_legendre_p(1, 1, p, rm));
    let res: ExactNum = expr!(hypergeom_2f1(mh, mh, x, mh), &mut ctx);
    debug_assert_eq!(res, mh.hypergeom_2f1(&mh, &x, &mh, p, rm, &mut cc));
    let res: ExactNum = expr!(betainc(x, x, mh), &mut ctx);
    debug_assert_eq!(res, x.betainc(&x, &mh, p, rm, &mut cc));
    let res: ExactNum = expr!(normal_pdf(x, x, x), &mut ctx);
    debug_assert_eq!(res, x.normal_pdf(&x, &x, p, rm, &mut cc));
    let z0 = ExactNum::from(0);
    let res: ExactNum = expr!(poisson_pmf(z0, x), &mut ctx);
    debug_assert_eq!(res, z0.poisson_pmf(&x, p, rm, &mut cc));

    let res: ExactNum = expr!(ei(x), &mut ctx);
    debug_assert_eq!(res, x.ei(p, rm, &mut cc));

    let res: ExactNum = expr!(si(x), &mut ctx);
    debug_assert_eq!(res, x.si(p, rm, &mut cc));

    let res: ExactNum = expr!(ci(x), &mut ctx);
    debug_assert_eq!(res, x.ci(p, rm, &mut cc));

    let ee = ExactNum::from(1).exp(p, rm, &mut cc);
    let res: ExactNum = expr!(li(ee), &mut ctx);
    debug_assert_eq!(res, ee.li(p, rm, &mut cc));

    let res: ExactNum = expr!(fresnel_s(x), &mut ctx);
    debug_assert_eq!(res, x.fresnel_s(p, rm, &mut cc));

    let res: ExactNum = expr!(fresnel_c(x), &mut ctx);
    debug_assert_eq!(res, x.fresnel_c(p, rm, &mut cc));

    let res: ExactNum = expr!(ai(x), &mut ctx);
    debug_assert_eq!(res, x.ai(p, rm, &mut cc));

    let res: ExactNum = expr!(bi(x), &mut ctx);
    debug_assert_eq!(res, x.bi(p, rm, &mut cc));

    let res: ExactNum = expr!(ldexp(x, 3), &mut ctx);
    debug_assert_eq!(res, x.ldexp(3, p, rm));

    let res: ExactNum = expr!(scalb(x, -1), &mut ctx);
    debug_assert_eq!(res, x.scalb(-1, p, rm));

    let res: ExactNum = expr!(logb(x), &mut ctx);
    debug_assert_eq!(res, x.logb(p, rm));
}

#[test]
fn macro_run_err_test() {
    let p = 192;
    let rm = RoundingMode::ToEven;
    let mut cc = Consts::new().unwrap();

    // Extra-precision cases below produce exponents far outside ±10k; Inf would
    // mask the compensation the test is actually checking.
    let mut ctx = Context::new(p, rm, Consts::new().unwrap(), EXPONENT_MIN, EXPONENT_MAX);

    let two = ExactNum::from(2);
    let ten = ExactNum::from(10);

    // sub cancellation of 256 bits
    let x = ExactNum::parse(
        "0.00000000000000000000000000000000000001",
        Radix::Dec,
        p,
        RoundingMode::None,
        &mut cc,
    );
    let y = ExactNum::parse(
        "1.57079632679489661923132169163975144209",
        Radix::Dec,
        p,
        RoundingMode::None,
        &mut cc,
    );

    let z = expr!(cos(x) - sin(y), &mut ctx);
    let cx = x.cos(p + 256, RoundingMode::None, &mut cc);
    let sy = y.sin(p + 256, RoundingMode::None, &mut cc);
    let r = cx.sub(&sy, p, rm);

    assert_eq!(r, z);

    // constants: pi, e, ln_2, ln_10
    let x = expr!(pi, &mut ctx);
    assert_eq!(x, ctx.const_pi());

    let x = expr!(e, &mut ctx);
    assert_eq!(x, ctx.const_e());

    let x = expr!(ln_2, &mut ctx);
    assert_eq!(x, ctx.const_ln2());

    let x = expr!(ln_10, &mut ctx);
    assert_eq!(x, ctx.const_ln10());

    let x = expr!(sqrt2, &mut ctx);
    assert_eq!(x, ctx.const_sqrt2());

    let x = expr!(phi, &mut ctx);
    assert_eq!(x, ctx.const_phi());

    let x = expr!(euler_gamma, &mut ctx);
    assert_eq!(x, ctx.const_euler_gamma());

    // ln
    for x in [
        ExactNum::from_words(&[234, 0, WORD_SIGNIFICANT_BIT], Sign::Pos, -123),
        ExactNum::from_words(&[234, 0, WORD_SIGNIFICANT_BIT], Sign::Neg, -123),
    ] {
        let y1 = x.exp(p * 3, RoundingMode::None, &mut cc);
        let mut z1 = y1.ln(p * 3, RoundingMode::None, &mut cc);
        z1.set_precision(p, rm).unwrap();

        let y2 = x.exp(p + 1, RoundingMode::None, &mut cc);
        let z2 = y2.ln(p, rm, &mut cc);

        let y = expr!(ln(exp(x)), &mut ctx);

        assert_eq!(z1, y);
        assert_ne!(z2, y);
    }

    // exp
    for x in [
        ExactNum::from_words(&[234, 0, WORD_SIGNIFICANT_BIT], Sign::Pos, 1000000),
        ExactNum::from_words(&[234, 0, WORD_SIGNIFICANT_BIT], Sign::Pos, -1000000),
    ] {
        let y1 = x.ln(p * 3, RoundingMode::None, &mut cc);
        let mut z1 = y1.exp(p * 3, RoundingMode::None, &mut cc);
        z1.set_precision(p, rm).unwrap();

        let y2 = x.ln(p, RoundingMode::None, &mut cc);
        let z2 = y2.exp(p, rm, &mut cc);

        let y = expr!(exp(ln(x)), &mut ctx);

        assert_eq!(z1, y);
        assert_ne!(z2, y);
    }

    // log2
    for x in [
        ExactNum::from_words(&[234, WORD_SIGNIFICANT_BIT], Sign::Pos, -123),
        ExactNum::from_words(&[234, WORD_SIGNIFICANT_BIT], Sign::Neg, -123),
    ] {
        let y1 = two.pow(&x, p * 2, RoundingMode::None, &mut cc);
        let mut z1 = y1.log2(p * 2, RoundingMode::None, &mut cc);
        z1.set_precision(p, rm).unwrap();

        let y2 = two.pow(&x, p + 1, RoundingMode::None, &mut cc);
        let z2 = y2.log2(p, rm, &mut cc);

        let y = expr!(log2(pow(2, x)), &mut ctx);

        assert_eq!(z1, y);
        assert_ne!(z2, y);
    }

    // log10
    for x in [
        ExactNum::from_words(&[234, WORD_SIGNIFICANT_BIT], Sign::Pos, -123),
        ExactNum::from_words(&[234, WORD_SIGNIFICANT_BIT], Sign::Neg, -123),
    ] {
        let y1 = ten.pow(&x, p * 2, RoundingMode::None, &mut cc);
        let mut z1 = y1.log10(p * 2, RoundingMode::None, &mut cc);
        z1.set_precision(p, rm).unwrap();

        let y2 = ten.pow(&x, p + 1, RoundingMode::None, &mut cc);
        let z2 = y2.log10(p, rm, &mut cc);

        let y = expr!(log10(pow(10, x)), &mut ctx);

        assert_eq!(z1, y);
        assert_ne!(z2, y);
    }

    // log
    for x in [
        ExactNum::from_words(&[234, WORD_SIGNIFICANT_BIT], Sign::Pos, -123),
        ExactNum::from_words(&[234, WORD_SIGNIFICANT_BIT], Sign::Neg, -123),
    ] {
        for b in [
            ExactNum::from_words(&[123, WORD_MAX], Sign::Pos, 0),
            ExactNum::from_words(&[123, WORD_SIGNIFICANT_BIT], Sign::Pos, 1),
        ] {
            let y1 = b.pow(&x, p * 4, RoundingMode::None, &mut cc);
            let mut z1 = y1.log(&b, p * 4, RoundingMode::None, &mut cc);
            z1.set_precision(p, rm).unwrap();

            let y2 = b.pow(&x, p + 1, RoundingMode::None, &mut cc);
            let z2 = y2.log(&b, p, rm, &mut cc);

            let y = expr!(log(pow(b, x), b), &mut ctx);

            assert_eq!(z1, y);
            assert_ne!(z2, y);
        }
    }

    let s1 = "1.0000000000000000000234";
    let s2 = "1.234567890123456789e+20";
    let b = ExactNum::parse(
        s1,
        zenith_float_num::Radix::Dec,
        p,
        RoundingMode::None,
        &mut cc,
    );
    let n = ExactNum::parse(
        s2,
        zenith_float_num::Radix::Dec,
        p,
        RoundingMode::None,
        &mut cc,
    );
    let y1 = b.pow(&n, p, rm, &mut cc);
    let b = ExactNum::parse(
        s1,
        zenith_float_num::Radix::Dec,
        p + 128,
        RoundingMode::None,
        &mut cc,
    );
    let n = ExactNum::parse(
        s2,
        zenith_float_num::Radix::Dec,
        p + 128,
        RoundingMode::None,
        &mut cc,
    );
    let mut y2 = b.pow(&n, p + 128, RoundingMode::None, &mut cc);
    y2.set_precision(p, rm).unwrap();

    let z = expr!(pow(s1, s2), &mut ctx);

    assert_eq!(y2, z);
    assert_ne!(y1, z);

    let s1 = "0.9999999999999999923456";
    let s2 = "-1.234567890123456789e+25";
    let b = ExactNum::parse(
        s1,
        zenith_float_num::Radix::Dec,
        p,
        RoundingMode::None,
        &mut cc,
    );
    let n = ExactNum::parse(
        s2,
        zenith_float_num::Radix::Dec,
        p,
        RoundingMode::None,
        &mut cc,
    );
    let y1 = b.pow(&n, p, rm, &mut cc);
    let b = ExactNum::parse(
        s1,
        zenith_float_num::Radix::Dec,
        p + 128,
        RoundingMode::None,
        &mut cc,
    );
    let n = ExactNum::parse(
        s2,
        zenith_float_num::Radix::Dec,
        p + 128,
        RoundingMode::None,
        &mut cc,
    );
    let mut y2 = b.pow(&n, p + 128, RoundingMode::None, &mut cc);
    y2.set_precision(p, rm).unwrap();

    let z = expr!(pow(s1, s2), &mut ctx);

    assert_eq!(y2, z);
    assert_ne!(y1, z);

    // sin
    let s1 = "1.2345678901234e+77";
    let n = ExactNum::parse(
        s1,
        zenith_float_num::Radix::Dec,
        128,
        RoundingMode::None,
        &mut cc,
    );
    let y1 = n.sin(128, rm, &mut cc);
    let n = ExactNum::parse(
        s1,
        zenith_float_num::Radix::Dec,
        320,
        RoundingMode::ToEven,
        &mut cc,
    );
    let mut y2 = n.sin(320, RoundingMode::None, &mut cc);
    y2.set_precision(128, rm).unwrap();

    let z = expr!(sin(s1), (128, rm, &mut cc));

    assert_ne!(y1, z);
    assert_eq!(y2, z);

    // cos
    let n = ExactNum::parse(
        s1,
        zenith_float_num::Radix::Dec,
        128,
        RoundingMode::None,
        &mut cc,
    );
    let y1 = n.cos(128, rm, &mut cc);
    let n = ExactNum::parse(
        s1,
        zenith_float_num::Radix::Dec,
        320,
        RoundingMode::None,
        &mut cc,
    );
    let mut y2 = n.cos(320, RoundingMode::None, &mut cc);
    y2.set_precision(128, rm).unwrap();

    let z = expr!(cos(s1), (128, rm, &mut cc));

    assert_ne!(y1, z);
    assert_eq!(y2, z);

    // tan
    let s1 = "1.1001001000011111101101010100010001000010110100011000010001101001100010011000110011000101000101110000000110111000001110011010001e+0";
    let n = ExactNum::parse(
        s1,
        zenith_float_num::Radix::Bin,
        p,
        RoundingMode::None,
        &mut cc,
    );
    let y1 = n.tan(p, rm, &mut cc);

    let n = ExactNum::parse(
        s1,
        zenith_float_num::Radix::Dec,
        y1.exponent().unwrap() as usize * 2 + p,
        RoundingMode::None,
        &mut cc,
    );
    let mut y2 = n.tan(
        y1.exponent().unwrap() as usize * 2 + p,
        RoundingMode::None,
        &mut cc,
    );
    y2.set_precision(p, rm).unwrap();

    let z = expr!(tan(s1), &mut ctx);

    assert_ne!(y1, z);
    assert_eq!(y2, z);

    // asin
    let mut x = cc.pi(p - WORD_BIT_SIZE, RoundingMode::None);
    x.set_exponent(1);

    let z = x.sin(p + 1, RoundingMode::None, &mut cc);
    let y1 = z.asin(p, rm, &mut cc);

    let z = x.sin(p * 2, RoundingMode::None, &mut cc);
    let mut y2 = z.asin(p * 2, RoundingMode::None, &mut cc);
    y2.set_precision(p, rm).unwrap();

    let z = expr!(asin(sin(x)), &mut ctx);

    assert_ne!(y1, z);
    assert_eq!(y2, z);

    // acos
    let x = ExactNum::from_words(&[234, WORD_MAX], Sign::Pos, -100);

    let z = x.cos(p + 1, RoundingMode::None, &mut cc);
    let y1 = z.acos(p, rm, &mut cc);

    let z = x.cos(p * 3, RoundingMode::None, &mut cc);
    let mut y2 = z.acos(p * 3, RoundingMode::None, &mut cc);
    y2.set_precision(p, rm).unwrap();

    let z = expr!(acos(cos(x)), &mut ctx);

    assert_ne!(y1, z);
    assert_eq!(y2, z);

    // atan
    let x = ExactNum::from_words(&[234, WORD_MAX], Sign::Neg, 128);
    let z = x.atan(p + 1, RoundingMode::None, &mut cc);
    let y1 = z.tan(p, rm, &mut cc);

    let z = x.atan(p + 256, RoundingMode::None, &mut cc);
    let mut y2 = z.tan(p + 256, RoundingMode::None, &mut cc);
    y2.set_precision(p, rm).unwrap();

    let z = expr!(tan(atan(x)), &mut ctx);

    assert_ne!(y1, z);
    assert_eq!(y2, z);

    // sinh
    let x = ExactNum::from_words(&[234, 0, WORD_SIGNIFICANT_BIT], Sign::Pos, 1000000);
    let y1 = x.asinh(p * 3, RoundingMode::None, &mut cc);
    let mut z1 = y1.sinh(p * 3, RoundingMode::None, &mut cc);
    z1.set_precision(p, rm).unwrap();

    let y2 = x.asinh(p, RoundingMode::None, &mut cc);
    let z2 = y2.sinh(p, rm, &mut cc);

    let y = expr!(sinh(asinh(x)), &mut ctx);

    assert_eq!(z1, y);
    assert_ne!(z2, y);

    // cosh
    let x = ExactNum::from_words(&[234, 0, WORD_SIGNIFICANT_BIT], Sign::Pos, 1000000);
    let y1 = x.acosh(p * 3, RoundingMode::None, &mut cc);
    let mut z1 = y1.cosh(p * 3, RoundingMode::None, &mut cc);
    z1.set_precision(p, rm).unwrap();

    let y2 = x.acosh(p, RoundingMode::None, &mut cc);
    let z2 = y2.cosh(p, rm, &mut cc);

    let y = expr!(cosh(acosh(x)), &mut ctx);

    assert_eq!(z1, y);
    assert_ne!(z2, y);

    // asinh
    let s1 = "6.1705892816164049e-1";

    let x = ExactNum::parse(
        s1,
        zenith_float_num::Radix::Dec,
        p,
        RoundingMode::None,
        &mut cc,
    );
    let y1 = x.asinh(p, rm, &mut cc);

    let x = ExactNum::parse(
        s1,
        zenith_float_num::Radix::Dec,
        p + 1,
        RoundingMode::None,
        &mut cc,
    );
    let mut y2 = x.asinh(p + 1, RoundingMode::None, &mut cc);
    y2.set_precision(p, rm).unwrap();

    let z = expr!(asinh(x), &mut ctx);

    assert_ne!(y1, z);
    assert_eq!(y2, z);

    // acosh
    let x = ExactNum::from_words(&[123, 123, WORD_SIGNIFICANT_BIT], Sign::Pos, -100);

    let z = x.cosh(p + 1, RoundingMode::None, &mut cc);
    let y1 = z.acosh(p, rm, &mut cc);

    let z = x.cosh(p + 256, RoundingMode::None, &mut cc);
    let mut y2 = z.acosh(p + 256, RoundingMode::None, &mut cc);
    y2.set_precision(p, rm).unwrap();

    let z = expr!(acosh(cosh(x)), &mut ctx);

    assert_ne!(y1, z);
    assert_eq!(y2, z);

    // atanh
    let x = ExactNum::from_words(&[123, 123, WORD_SIGNIFICANT_BIT], Sign::Pos, 7);

    let z = x.tanh(p + 1, RoundingMode::None, &mut cc);
    let y1 = z.atanh(p, rm, &mut cc);

    let z = x.tanh(p + 256, RoundingMode::None, &mut cc);
    let mut y2 = z.atanh(p + 256, RoundingMode::None, &mut cc);
    y2.set_precision(p, rm).unwrap();

    let z = expr!(atanh(tanh(x)), &mut ctx);

    assert_ne!(y1, z);
    assert_eq!(y2, z);
}

// test precision range for error compensation
#[test]
fn macro_run_misc_test() {
    let p = 256;
    let rm = RoundingMode::None;
    let mut cc = Consts::new().unwrap();
    let pi1 = cc.pi(p, RoundingMode::None);
    let mut pi2 = pi1.clone();
    pi2.set_exponent(1);
    let mut ctx = Context::new(p, rm, cc, -1000, 1000);
    let one = ExactNum::from_u32(1, p);

    // literals: float
    let z = expr!(2e-301, &mut ctx);
    assert!(!z.is_zero());

    let z = expr!(2e-302, &mut ctx);
    assert!(z.is_zero());

    // literals: strings
    let z = expr!("2e-301", &mut ctx);
    assert!(!z.is_zero());

    let z = expr!("2e-302", &mut ctx);
    assert!(z.is_zero());

    // exceed output exponent
    let z = expr!(2e+151 / 2e-151, &mut ctx);
    assert!(z.is_inf_pos());

    let z = expr!(2e-151 / 2e+151, &mut ctx);
    assert!(z.is_zero());

    // ln, log2, lg
    let z = expr!(1 + 2e-301, &mut ctx);
    assert!(z == one);

    let z = expr!(ln(1 + 2e-301), &mut ctx);
    assert!(!z.is_zero());

    let z = expr!(ln(1 + 2e-302), &mut ctx);
    assert!(z.is_zero());

    let z = expr!(log10(1 + 2e-301), &mut ctx);
    assert!(!z.is_zero());

    let z = expr!(log10(1 + 2e-302), &mut ctx);
    assert!(z.is_zero());

    let z = expr!(log2(1 + 2e-301), &mut ctx);
    assert!(!z.is_zero());

    let z = expr!(log2(1 + 2e-302), &mut ctx);
    assert!(z.is_zero());

    // log
    let z = expr!(log(1 + 2e-301, 1 + 2e-301), &mut ctx);
    assert!(z == one);

    let z = expr!(log(2e+300, 2e+300), &mut ctx);
    assert!(z == one);

    let z = expr!(log(2e+300, 1 + 2e-301), &mut ctx);
    assert!(z.is_inf());

    let z = expr!(log(1 + 2e-301, 2e+300), &mut ctx);
    assert!(z.is_zero());

    let z = expr!(log(1 + 2e-301, 1 + 2e-302), &mut ctx);
    assert!(z.is_inf());

    let z = expr!(log(1 + 2e-302, 1 + 2e-301), &mut ctx);
    assert!(z.is_zero());

    let z = expr!(log(1 + 2e-302, 1 + 2e-302), &mut ctx);
    assert!(z.is_nan());

    // pow
    let z = expr!(pow(1 + 2e-301, 2e+300), &mut ctx);
    assert!(!z.is_inf());
    assert!(z > one);

    let z = expr!(pow(1 + 2e-280, 2e+300), &mut ctx);
    assert!(z.is_inf());

    let z = expr!(pow(1 + 2e-300, 2e+300 / 2e-2), &mut ctx);
    assert!(!z.is_inf());
    assert!(z > one);

    // sin, cos, tan
    let z = expr!(sin(pi1 + 2e-300), &mut ctx);
    assert!(z.is_zero() || z.exponent().unwrap() < 0);

    let z = expr!(cos(pi2 + 2e-300), &mut ctx);
    assert!(z.is_zero() || z.exponent().unwrap() < 0);

    let z = expr!(tan(pi1 + 2e-300), &mut ctx);
    assert!(z.is_zero() || z.exponent().unwrap() < 0);

    let z = expr!(tan(pi2 + 2e-300), &mut ctx);
    assert!(z.is_zero() || z.exponent().unwrap() > 1);

    // asin
    let x = expr!(asin(1), &mut ctx);
    let z = expr!(asin(1 - 2e-301), &mut ctx);
    assert!(z == x);
    let z = expr!(asin(1 - 2e-256), &mut ctx);
    assert!(z == x);
    let z = expr!(asin(1 - 2e-150), &mut ctx);
    assert!(z != x);

    // acos
    let x = expr!(acos(1), &mut ctx);
    let z = expr!(acos(1 - 2e-302), &mut ctx);
    assert!(z == x);

    let z = expr!(acos(1 - 2e-301), &mut ctx);
    assert!(z.exponent().unwrap() < 0);

    // acos
    let x = expr!(acos(-1), &mut ctx);
    let z = expr!(acos(-1 + 2e-302), &mut ctx);
    assert!(z == x);

    let z = expr!(acos(-1 + 2e-151), &mut ctx);
    assert!(z != x);

    // acosh
    let x = expr!(acosh(1), &mut ctx);
    let z = expr!(acosh(1 + 2e-302), &mut ctx);
    assert!(z == x);

    let z = expr!(acosh(1 + 2e-301), &mut ctx);
    assert!(z != x);

    // atanh
    let x = expr!(atanh(1), &mut ctx);
    let z = expr!(atanh(1 - 2e-302), &mut ctx);
    assert!(z == x);

    let z = expr!(atanh(1 - 2e-301), &mut ctx);
    assert!(z != x);

    let x = expr!(atanh(-1), &mut ctx);
    let z = expr!(atanh(-1 + 2e-302), &mut ctx);
    assert!(z == x);

    let z = expr!(atanh(-1 + 2e-301), &mut ctx);
    assert!(z != x);

    // infinitely close to 1: 0.99999999(9)
    let z = expr!(5 * (1 / 5), &mut ctx);
    assert!(z < one);

    let z = expr!(ln(5 * (1 / 5)), &mut ctx);
    assert!(z.is_zero());
}

#[test]
fn expr_erf_plus_erfc_is_one_at_256() {
    let p = 256;
    let rm = RoundingMode::ToEven;
    let mut ctx = Context::new(p, rm, Consts::new().unwrap(), -100000, 100000);
    let one = ExactNum::from_u8(1, p);
    let x = ExactNum::from_u8(1, p);
    let res: ExactNum = expr!(erf(x) + erfc(x), &mut ctx);
    assert_eq!(res.cmp(&one), Some(0));
    assert!(res.err().is_none());
}

#[test]
fn expr_bessel_j0_y0_sq_working_prec() {
    let p = 256;
    let rm = RoundingMode::ToEven;
    let mut ctx = Context::new(p, rm, Consts::new().unwrap(), -100000, 100000);
    let mut cc = Consts::new().unwrap();
    let x = ExactNum::from_u8(1, p);
    let nu0 = ExactNum::from_u8(0, p);
    let got: ExactNum = expr!(
        bessel_j(x, 0) * bessel_j(x, 0) + bessel_y(x, nu0) * bessel_y(x, nu0),
        &mut ctx
    );
    let p_wrk = p + WORD_BIT_SIZE + 16;
    let j = x.bessel_j(0, p_wrk, RoundingMode::None, &mut cc);
    let y = x.bessel_y(&nu0, p_wrk, RoundingMode::None, &mut cc);
    let mut expect = j
        .mul(&j, p_wrk, RoundingMode::None)
        .add(&y.mul(&y, p_wrk, RoundingMode::None), p_wrk, RoundingMode::None);
    expect.set_precision(p, rm).unwrap();
    assert_eq!(got.cmp(&expect), Some(0));
    assert!(got.err().is_none());
}

#[test]
fn cexpr_erf_uses_p_wrk() {
    let p = 256;
    let rm = RoundingMode::ToEven;
    let mut ctx = Context::new(p, rm, Consts::new().unwrap(), -100000, 100000);
    let mut cc = Consts::new().unwrap();
    let z = ExactComplex::new(ExactNum::from_u8(1, p), ExactNum::from_u8(1, p));
    let got = cexpr!(erf(z), &mut ctx);
    let p_wrk = p + WORD_BIT_SIZE + 2;
    let mut expect = z.erf(p_wrk, RoundingMode::None, &mut cc);
    expect.set_precision(p, rm).unwrap();
    cplx_eq("cexpr erf p_wrk", &got, &expect);
    assert!(got.re().err().is_none());
    assert!(got.im().err().is_none());
}
