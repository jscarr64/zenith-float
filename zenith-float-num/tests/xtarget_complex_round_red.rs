//! RED: complex special functions must return p-bit results (rm honoured; None truncates at
//! rp(p, 64)) identically on every target. Today the complex paths round to `round_p(p)`, i.e. to
//! native words (p = 20 gives 64 bits on x86_64 and 32 on i686), and the None path of asin and
//! bessel_i is not a correct truncation on i686.
//!
//! Value rows: scripts/gen_complex_xtarget_ref.py, env mpmath-1.3.0 (mpmath 1.3.0, Python 3.12.14),
//! mp.prec = 1024. asin: principal branch, DLMF 4.23.19. bessel_i: I_nu(x), nu = 1/2 + i/4,
//! x = 3/2 + i/4, DLMF 10.25.2. Generator output SHA256:
//! 25dd70bd96c58429c6829fd4c8dce6c3c69cda39a44b3f3bae8538ff7fcee1d1
//! Width rows need no reference value: a p-bit rounding is unchanged by rounding to p bits again,
//! and is stored in at most rp(p, 64) bits.

use zenith_float_num::{Consts, ExactComplex, ExactNum, Radix, RoundingMode};

const REF: &[(&str, usize, &str, &str, &str)] = &[
    ("asin", 20, "ToEven", "0.10000000011010010111", "0.01001000000010011001"),
    ("asin", 20, "Up", "0.10000000011010011", "0.010010000000100110011"),
    ("asin", 20, "None", "0.10000000011010010111000000010100001100111100110001110000101", "0.01001000000010011001001001101010110011000000001000111010011000111"),
    ("asin", 65, "ToEven", "0.10000000011010010111000000010100001100111100110001110000101000001", "0.010010000000100110010010011010101100110000000010001110100110001111"),
    ("asin", 65, "Up", "0.10000000011010010111000000010100001100111100110001110000101000001", "0.01001000000010011001001001101010110011000000001000111010011001"),
    ("asin", 65, "None", "0.10000000011010010111000000010100001100111100110001110000101000000101111010100010010110000001101001011110111011011101101110110101", "0.010010000000100110010010011010101100110000000010001110100110001111000001111011001101101011000101010010111010110000011100001001101"),
    ("bessel_i", 20, "ToEven", "1.011000001001001101", "0.00010101000111110011"),
    ("bessel_i", 20, "Up", "1.011000001001001101", "0.00010101000111110011001"),
    ("bessel_i", 20, "None", "1.0110000010010011001110100110001100100010001100011100111001111", "0.00010101000111110011000000010101000001000100110001101001100110011"),
    ("bessel_i", 65, "ToEven", "1.011000001001001100111010011000110010001000110001110011100111101", "0.0001010100011111001100000001010100000100010011000110100110011001101"),
    ("bessel_i", 65, "Up", "1.011000001001001100111010011000110010001000110001110011100111101", "0.0001010100011111001100000001010100000100010011000110100110011001101"),
    ("bessel_i", 65, "None", "1.0110000010010011001110100110001100100010001100011100111001111001100000011100010110011001111101011110000111110000001010111001111", "0.0001010100011111001100000001010100000100010011000110100110011001100110010011110000010100110100110001111100111100111101101110111"),
];

fn rm_of(name: &str) -> Option<RoundingMode> {
    match name {
        "ToEven" => Some(RoundingMode::ToEven),
        "Up" => Some(RoundingMode::Up),
        "None" => Some(RoundingMode::None),
        _ => None,
    }
}

fn num(s: &str, cc: &mut Consts) -> ExactNum {
    ExactNum::parse(s, Radix::Bin, 512, RoundingMode::None, cc)
}

fn cnum(re: &str, im: &str, cc: &mut Consts) -> ExactComplex {
    ExactComplex::new(num(re, cc), num(im, cc))
}

/// z = 1/2 + i/4, w = 1/4 + 3i/4, v = 3/2 + i/4 (exact binary values).
fn inputs(cc: &mut Consts) -> (ExactComplex, ExactComplex, ExactComplex) {
    (
        cnum("0.1", "0.01", cc),
        cnum("0.01", "0.11", cc),
        cnum("1.1", "0.01", cc),
    )
}

const FNS: [&str; 26] = [
    "ei",
    "si",
    "ci",
    "li",
    "fresnel_s",
    "fresnel_c",
    "ln",
    "asin",
    "elliptic_k",
    "elliptic_e_complete",
    "elliptic_f",
    "elliptic_e",
    "elliptic_pi_complete",
    "elliptic_pi",
    "ai",
    "bi",
    "erf",
    "erfc",
    "gamma",
    "ln_gamma",
    "digamma",
    "bessel_j_nu",
    "bessel_y",
    "bessel_i",
    "bessel_k",
    "hypergeom_2f1",
];

fn eval(f: &str, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<ExactComplex, String> {
    let (z, w, v) = inputs(cc);
    Ok(match f {
        "ei" => z.ei(p, rm, cc),
        "si" => z.si(p, rm, cc),
        "ci" => z.ci(p, rm, cc),
        "li" => v.li(p, rm, cc),
        "fresnel_s" => z.fresnel_s(p, rm, cc),
        "fresnel_c" => z.fresnel_c(p, rm, cc),
        "ln" => z.ln(p, rm, cc),
        "asin" => z.asin(p, rm, cc),
        "elliptic_k" => z.elliptic_k(p, rm, cc),
        "elliptic_e_complete" => z.elliptic_e_complete(p, rm, cc),
        "elliptic_f" => z.elliptic_f(&w, p, rm, cc),
        "elliptic_e" => z.elliptic_e(&w, p, rm, cc),
        "elliptic_pi_complete" => z.elliptic_pi_complete(&w, p, rm, cc),
        "elliptic_pi" => z.elliptic_pi(&w, &w, p, rm, cc),
        "ai" => z.ai(p, rm, cc),
        "bi" => z.bi(p, rm, cc),
        "erf" => z.erf(p, rm, cc),
        "erfc" => z.erfc(p, rm, cc),
        "gamma" => v.gamma(p, rm, cc),
        "ln_gamma" => v.ln_gamma(p, rm, cc),
        "digamma" => v.digamma(p, rm, cc),
        "bessel_j_nu" => v.bessel_j_nu(&z, p, rm, cc),
        "bessel_y" => v.bessel_y(&z, p, rm, cc),
        "bessel_i" => v.bessel_i(&z, p, rm, cc),
        "bessel_k" => v.bessel_k(&z, p, rm, cc),
        "hypergeom_2f1" => z.hypergeom_2f1(&w, &v, &z, p, rm, cc),
        _ => return Err(format!("bad fn {f}")),
    })
}

/// `a` holds a p-bit value (re-rounding `b`, an identical evaluation, to p bits changes nothing)
/// and is stored in at most `cap` bits.
fn p_bit_part(a: &ExactNum, mut b: ExactNum, p: usize, cap: usize) -> Result<bool, String> {
    b.set_precision(p, RoundingMode::ToEven)
        .map_err(|e| format!("{e:?}"))?;
    let sig = a.mantissa_max_bit_len().ok_or("no mantissa")?;
    Ok(*a == b && sig <= cap)
}

#[test]
fn complex_values_are_target_independent() -> Result<(), String> {
    let mut cc = Consts::new().map_err(|e| format!("{e:?}"))?;
    let mut bad = Vec::new();
    for &(f, p, name, re, im) in REF {
        let rm = rm_of(name).ok_or_else(|| format!("bad rm {name}"))?;
        let got = eval(f, p, rm, &mut cc)?;
        let want = cnum(re, im, &mut cc);
        if *got.re() != *want.re() || *got.im() != *want.im() {
            bad.push(format!("{f} p={p} {name}: got {got:?}, want {re} + i {im}"));
        }
    }
    if bad.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "{} of {} rows differ:\n{}",
            bad.len(),
            REF.len(),
            bad.join("\n")
        ))
    }
}

#[test]
fn complex_results_have_p_bits() -> Result<(), String> {
    let mut cc = Consts::new().map_err(|e| format!("{e:?}"))?;
    let mut bad = Vec::new();
    let mut total = 0usize;
    for p in [20usize, 65] {
        let cap = p.div_ceil(64).checked_mul(64).ok_or("cap overflow")?;
        for rm in [RoundingMode::ToEven, RoundingMode::Up] {
            for f in FNS {
                total = total.checked_add(1).ok_or("count overflow")?;
                let a = eval(f, p, rm, &mut cc)?;
                let (bre, bim) = eval(f, p, rm, &mut cc)?.into_parts();
                if !p_bit_part(a.re(), bre, p, cap)? || !p_bit_part(a.im(), bim, p, cap)? {
                    bad.push(format!("{f} p={p} {rm:?}: {a:?} is not a p-bit value"));
                }
            }
        }
    }
    if bad.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "{} of {total} rows are not p-bit values:\n{}",
            bad.len(),
            bad.join("\n")
        ))
    }
}
