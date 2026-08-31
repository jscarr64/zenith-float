//! Identity golds for functions GNU MPFR / MPC do not implement.
//!
//! No invented decimal tables. Each check is an identity already used as a
//! unit gold: odd/even, complementary definitions, `K_{1/2}`, Wronskian,
//! `_2F1(...,0)=1`, `K(0)=π/2`.

use crate::mpfr::common::reset_test_rng;
use zenith_float_num::{Consts, ExactComplex, ExactNum, RoundingMode, WORD_BIT_SIZE};

fn bits_agree(a: &ExactNum, b: &ExactNum, p: usize, min_bits: i32, label: &str) {
    let d = a.sub(b, p, RoundingMode::None).abs();
    if d.is_zero() {
        return;
    }
    let rel = d.exponent().unwrap_or(0) - a.exponent().unwrap_or(0);
    assert!(
        rel < -min_bits,
        "{label}: relative exponent {rel} (a.exp={:?} b.exp={:?})",
        a.exponent(),
        b.exponent()
    );
}

fn tiny(x: &ExactNum, p: usize) -> bool {
    x.is_zero() || x.exponent().is_some_and(|e| e < -((p as i32) / 4))
}

#[test]
fn mpfr_identity_no_oracle_specials() {
    reset_test_rng();
    let p = 2 * WORD_BIT_SIZE;
    let rm = RoundingMode::ToEven;
    let mut cc = Consts::new().unwrap();
    let one = ExactNum::from_u8(1, p);
    let zero = ExactNum::new(p);
    let two = ExactNum::from_u8(2, p);
    let half = one.div(&two, p, rm);

    let s1 = one.si(p, rm, &mut cc);
    let sn = one.neg().si(p, rm, &mut cc);
    bits_agree(&s1, &sn.neg(), p, (p as i32) / 4, "Si odd");
    assert!(zero.si(p, rm, &mut cc).is_zero());
    let hp = cc.pi(p, rm).div(&two, p, rm);
    bits_agree(&zenith_float_num::INF_POS.si(p, rm, &mut cc), &hp, p, (p as i32) / 4, "Si(+∞)");

    let e = one.exp(p, rm, &mut cc);
    bits_agree(&e.li(p, rm, &mut cc), &one.ei(p, rm, &mut cc), p, (p as i32) / 4, "li(e)=Ei(1)");

    assert!(zero.fresnel_s(p, rm, &mut cc).is_zero());
    assert!(zero.fresnel_c(p, rm, &mut cc).is_zero());
    let fs = one.fresnel_s(p, rm, &mut cc);
    bits_agree(&fs, &one.neg().fresnel_s(p, rm, &mut cc).neg(), p, (p as i32) / 4, "S odd");
    bits_agree(
        &zenith_float_num::INF_POS.fresnel_s(p, rm, &mut cc),
        &half,
        p,
        (p as i32) / 4,
        "S(+∞)=1/2",
    );

    let k_half = one.bessel_k(&half, p, rm, &mut cc);
    let want_k = hp
        .sqrt(p, rm)
        .mul(&one.neg().exp(p, rm, &mut cc), p, rm);
    bits_agree(&k_half, &want_k, p, (p as i32) / 4, "K_{1/2}(1)");

    let ai = one.ai(p, rm, &mut cc);
    let bi = one.bi(p, rm, &mut cc);
    let aip = one.ai_prime(p, rm, &mut cc);
    let bip = one.bi_prime(p, rm, &mut cc);
    let wr = ai.mul(&bip, p, rm).sub(&aip.mul(&bi, p, rm), p, rm);
    bits_agree(&wr, &one.div(&cc.pi(p, rm), p, rm), p, 40, "Airy Wronskian");

    bits_agree(&zero.elliptic_k(p, rm, &mut cc), &hp, p, (p as i32) / 4, "K(0)=π/2");

    let a = ExactNum::from_u8(2, p);
    let b = ExactNum::from_u8(3, p);
    let c = ExactNum::from_u8(4, p);
    bits_agree(
        &a.hypergeom_2f1(&b, &c, &zero, p, rm, &mut cc),
        &one,
        p,
        (p as i32) / 4,
        "_2F1(a,b;c;0)=1",
    );

    let z = ExactComplex::new(one.clone(), one);
    let ez = z.erf(p, rm, &mut cc);
    let em = ExactComplex::new(z.re().neg(), z.im().neg()).erf(p, rm, &mut cc);
    assert!(tiny(&ez.re().add(&em.re(), p, rm), p) && tiny(&ez.im().add(&em.im(), p, rm), p));
    let erfc_z = z.erfc(p, rm, &mut cc);
    let one_c = ExactComplex::one(p);
    let id = one_c.sub(&ez, p, rm);
    assert!(
        tiny(&erfc_z.re().sub(id.re(), p, rm), p) && tiny(&erfc_z.im().sub(id.im(), p, rm), p),
        "erfc=1-erf"
    );

    let five = ExactComplex::from_real(ExactNum::from_u8(5, p), p);
    let g5 = five.gamma(p, rm, &mut cc);
    let tf = ExactNum::from_u8(24, p);
    assert!(tiny(&g5.re().sub(&tf, p, rm), p) && tiny(g5.im(), p), "Γ(5)=24");
}
