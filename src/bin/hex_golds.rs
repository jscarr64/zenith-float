//! Print canonical `to_bytes()` hex for one representative of each special.
//!
//! Output is independent of `WORD_BIT_SIZE` (binary interchange uses big-endian `u32` limbs).
//! Host and cross scripts diff this against `golds/hex/reference.txt`.

use std::fmt::Write;

use zenith_float::{Consts, ExactComplex, ExactNum, RoundingMode};

fn hex_bytes(bytes: &[u8]) -> Result<String, std::fmt::Error> {
    let mut s = String::with_capacity(bytes.len() * 2);
    for b in bytes {
        write!(&mut s, "{b:02x}")?;
    }
    Ok(s)
}

fn emit(name: &str, x: &ExactNum) -> Result<(), Box<dyn std::error::Error>> {
    let bytes = x.to_bytes()?;
    println!("{name}\t{}", hex_bytes(&bytes)?);
    Ok(())
}

fn emit_c(name: &str, z: &ExactComplex) -> Result<(), Box<dyn std::error::Error>> {
    let re = z.re().to_bytes()?;
    let im = z.im().to_bytes()?;
    println!("{name}\t{}\t{}", hex_bytes(&re)?, hex_bytes(&im)?);
    Ok(())
}

fn n(v: i32, p: usize) -> ExactNum {
    ExactNum::from_i32(v, p)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let p = 64usize;
    let rm = RoundingMode::ToEven;
    let mut cc = Consts::new()?;
    println!("# zenith-float hex limb golds p={p}");

    let one = n(1, p);
    let two = n(2, p);
    let half = one.div(&two, p, rm);

    emit("add_1_2", &one.add(&two, p, rm))?;
    emit("sqrt_2", &two.sqrt(p, rm))?;
    emit("exp_1", &one.exp(p, rm, &mut cc))?;
    emit("sin_1", &one.sin(p, rm, &mut cc))?;
    emit("erf_1", &one.erf(p, rm, &mut cc))?;
    emit("erfc_1", &one.erfc(p, rm, &mut cc))?;
    emit("gamma_5", &n(5, p).gamma(p, rm, &mut cc))?;
    emit("ln_gamma_5", &n(5, p).ln_gamma(p, rm, &mut cc))?;
    emit("digamma_1", &one.digamma(p, rm, &mut cc))?;
    emit("gammainc_1_1", &one.gammainc(&one, p, rm, &mut cc))?;
    emit(
        "gammainc_upper_1_1",
        &one.gammainc_upper(&one, p, rm, &mut cc),
    )?;
    emit("ei_1", &one.ei(p, rm, &mut cc))?;
    emit("si_1", &one.si(p, rm, &mut cc))?;
    emit("ci_1", &one.ci(p, rm, &mut cc))?;
    emit("li_2", &two.li(p, rm, &mut cc))?;
    emit("fresnel_s_1", &one.fresnel_s(p, rm, &mut cc))?;
    emit("fresnel_c_1", &one.fresnel_c(p, rm, &mut cc))?;
    emit("ai_0", &n(0, p).ai(p, rm, &mut cc))?;
    emit("bi_0", &n(0, p).bi(p, rm, &mut cc))?;
    emit("bessel_j_0_1", &one.bessel_j(0, p, rm, &mut cc))?;
    emit("bessel_j_nu_half_1", &one.bessel_j_nu(&half, p, rm, &mut cc))?;
    emit("bessel_y_0_1", &one.bessel_y(&n(0, p), p, rm, &mut cc))?;
    emit("bessel_i_0_1", &one.bessel_i(&n(0, p), p, rm, &mut cc))?;
    emit("bessel_k_0_1", &one.bessel_k(&n(0, p), p, rm, &mut cc))?;
    emit("elliptic_k_half", &half.elliptic_k(p, rm, &mut cc))?;
    emit(
        "elliptic_e_complete_half",
        &half.elliptic_e_complete(p, rm, &mut cc),
    )?;
    emit(
        "elliptic_f_half_half",
        &half.elliptic_f(&half, p, rm, &mut cc),
    )?;
    emit("jacobi_sn_1_half", &one.jacobi_sn(&half, p, rm, &mut cc))?;
    emit("legendre_p2_0", &n(0, p).legendre_p(2, p, rm))?;
    emit(
        "hypergeom_2f1_gauss",
        &one.hypergeom_2f1(&one, &two, &half, p, rm, &mut cc),
    )?;
    emit("betainc_1_1_half", &one.betainc(&one, &half, p, rm, &mut cc))?;
    emit(
        "normal_pdf_0_0_1",
        &n(0, p).normal_pdf(&n(0, p), &one, p, rm, &mut cc),
    )?;

    let z = ExactComplex::new(n(1, p), n(1, p));
    emit_c("complex_erf_1_i", &z.erf(p, rm, &mut cc))?;
    emit_c("complex_gamma_1_i", &z.gamma(p, rm, &mut cc))?;
    emit_c("complex_ai_1_i", &z.ai(p, rm, &mut cc))?;
    emit_c(
        "complex_bessel_j_nu_0",
        &z.bessel_j_nu(&ExactComplex::from_real(n(0, p), p), p, rm, &mut cc),
    )?;
    Ok(())
}
