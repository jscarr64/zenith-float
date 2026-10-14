//! mpmath goldens for `tests/oracle/goldens/specials`.
//!
//! `cargo test --test specials_oracle` runs the fast tier (2 s per row).
//! `ORACLE_TIER=deep` also runs the deep tier (10 s per row). A panic or a
//! timeout is a failure. `CorrectlyRounded` rows and rows marked `exact` are
//! limb-identical. `lost_max` is the bit slack a folded fixture already
//! allowed. `WorkingPrecision` uses the class budget. `Unspecified` checks
//! the domain, the sign, and two decimal digits.

use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::mpsc;
use std::thread;
use std::time::Duration;

use zenith_float_num::{
    accuracy_of, Accuracy, Consts, Error, ExactComplex, ExactNum, Exponent, Radix, RoundingMode,
    Sign, Word, INF_NEG, INF_POS,
};

#[derive(Clone, Debug, PartialEq)]
enum Canon {
    Finite {
        sign: String,
        exp: Exponent,
        limbs: Vec<String>,
    },
    Inf(String),
    Nan(Option<String>),
    Complex {
        re: Box<Canon>,
        im: Box<Canon>,
    },
}

enum Tok {
    N(ExactNum),
    I(usize),
    Z(ExactComplex),
}

fn tier_selected(tier: &str) -> bool {
    match std::env::var("ORACLE_TIER") {
        Ok(value) if value == "deep" => true,
        Ok(value) if value == "fast" => tier == "fast",
        _ => tier == "fast",
    }
}

fn mode_of(name: &str) -> RoundingMode {
    match name {
        "ToEven" => RoundingMode::ToEven,
        "ToOdd" => RoundingMode::ToOdd,
        "Up" => RoundingMode::Up,
        "Down" => RoundingMode::Down,
        "ToZero" => RoundingMode::ToZero,
        "FromZero" => RoundingMode::FromZero,
        "None" => RoundingMode::None,
        _ => RoundingMode::None,
    }
}

fn err_name(err: Option<Error>) -> Option<String> {
    match err {
        None => None,
        Some(Error::InvalidArgument) => Some("InvalidArgument".to_string()),
        Some(Error::DivisionByZero) => Some("DivisionByZero".to_string()),
        Some(Error::MemoryAllocation) => Some("MemoryAllocation".to_string()),
        Some(Error::PrecisionRetryExhausted) => Some("PrecisionRetryExhausted".to_string()),
        Some(Error::ExponentOverflow(_)) => Some("ExponentOverflow".to_string()),
    }
}

fn canon_of(n: &ExactNum) -> Canon {
    if n.is_nan() {
        return Canon::Nan(err_name(n.err()));
    }
    if n.is_inf() {
        let sign = if n.sign() == Some(Sign::Neg) { "-" } else { "+" };
        return Canon::Inf(sign.to_string());
    }
    let sign = if n.sign() == Some(Sign::Neg) { "-" } else { "+" };
    let exp = n.exponent().unwrap_or(0);
    let limbs = n
        .mantissa_digits()
        .unwrap_or(&[])
        .iter()
        .map(|word| format!("{word:016x}"))
        .collect();
    Canon::Finite {
        sign: sign.to_string(),
        exp,
        limbs,
    }
}

fn complex_of(z: &ExactComplex) -> Canon {
    Canon::Complex {
        re: Box::new(canon_of(z.re())),
        im: Box::new(canon_of(z.im())),
    }
}

fn expect_of(v: &serde_json::Value) -> Canon {
    let kind = v.get("k").and_then(|k| k.as_str()).unwrap_or("");
    if kind == "c" {
        return Canon::Complex {
            re: Box::new(expect_of(v.get("re").unwrap_or(&serde_json::Value::Null))),
            im: Box::new(expect_of(v.get("im").unwrap_or(&serde_json::Value::Null))),
        };
    }
    if kind == "nan" {
        let err = match v.get("err") {
            Some(serde_json::Value::Null) | None => None,
            Some(serde_json::Value::String(s)) => Some(s.clone()),
            Some(_) => Some("?".to_string()),
        };
        return Canon::Nan(err);
    }
    if kind == "inf" {
        return Canon::Inf(
            v.get("s")
                .and_then(|s| s.as_str())
                .unwrap_or("+")
                .to_string(),
        );
    }
    let limbs = v
        .get("m")
        .and_then(|m| m.as_array())
        .map(|arr| {
            arr.iter()
                .filter_map(|w| w.as_str().map(|s| s.to_string()))
                .collect()
        })
        .unwrap_or_default();
    Canon::Finite {
        sign: v
            .get("s")
            .and_then(|s| s.as_str())
            .unwrap_or("+")
            .to_string(),
        exp: v.get("e").and_then(|e| e.as_i64()).unwrap_or(0) as Exponent,
        limbs,
    }
}

fn arg_value(tok: &str, p: usize, cc: &mut Consts) -> ExactNum {
    ExactNum::parse(tok, Radix::Dec, p, RoundingMode::ToEven, cc)
}

fn parse_tok(tok: &str, p: usize, cc: &mut Consts) -> Result<Tok, String> {
    if let Some(rest) = tok.strip_prefix("n:") {
        let n = rest.parse::<usize>().map_err(|_| format!("int {tok}"))?;
        return Ok(Tok::I(n));
    }
    if let Some(rest) = tok.strip_prefix("c:") {
        let (re, im) = rest
            .split_once(',')
            .ok_or_else(|| format!("complex {tok}"))?;
        return Ok(Tok::Z(ExactComplex::new(
            arg_value(re, p, cc),
            arg_value(im, p, cc),
        )));
    }
    Ok(Tok::N(arg_value(tok, p, cc)))
}

fn num(toks: &[Tok], i: usize) -> Result<&ExactNum, String> {
    match toks.get(i) {
        Some(Tok::N(n)) => Ok(n),
        _ => Err(format!("arg {i} is not real")),
    }
}

fn idx(toks: &[Tok], i: usize) -> Result<usize, String> {
    match toks.get(i) {
        Some(Tok::I(n)) => Ok(*n),
        _ => Err(format!("arg {i} is not an integer")),
    }
}

fn cpx(toks: &[Tok], i: usize) -> Result<&ExactComplex, String> {
    match toks.get(i) {
        Some(Tok::Z(z)) => Ok(z),
        _ => Err(format!("arg {i} is not complex")),
    }
}

fn json_nums(
    v: Option<&serde_json::Value>,
    p: usize,
    cc: &mut Consts,
) -> Result<Vec<ExactNum>, String> {
    let Some(arr) = v.and_then(|v| v.as_array()) else {
        return Ok(Vec::new());
    };
    let mut out = Vec::with_capacity(arr.len());
    for item in arr {
        let s = item.as_str().ok_or("parameter is not a string")?;
        out.push(arg_value(s, p, cc));
    }
    Ok(out)
}

fn json_pairs(
    v: Option<&serde_json::Value>,
    p: usize,
    cc: &mut Consts,
) -> Result<Vec<(ExactNum, ExactNum)>, String> {
    let Some(arr) = v.and_then(|v| v.as_array()) else {
        return Ok(Vec::new());
    };
    let mut out = Vec::with_capacity(arr.len());
    for item in arr {
        let pair = item.as_array().ok_or("pair")?;
        if pair.len() != 2 {
            return Err("pair length".to_string());
        }
        let a = pair[0].as_str().ok_or("pair a")?;
        let b = pair[1].as_str().ok_or("pair b")?;
        out.push((arg_value(a, p, cc), arg_value(b, p, cc)));
    }
    Ok(out)
}

fn eval_unary(
    op: &str,
    x: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNum> {
    Some(match op {
        "erf" => x.erf(p, rm, cc),
        "erfc" => x.erfc(p, rm, cc),
        "gamma" => x.gamma(p, rm, cc),
        "ln_gamma" => x.ln_gamma(p, rm, cc),
        "digamma" => x.digamma(p, rm, cc),
        "ai" => x.ai(p, rm, cc),
        "bi" => x.bi(p, rm, cc),
        "ai_prime" => x.ai_prime(p, rm, cc),
        "bi_prime" => x.bi_prime(p, rm, cc),
        "ei" => x.ei(p, rm, cc),
        "si" => x.si(p, rm, cc),
        "ci" => x.ci(p, rm, cc),
        "li" => x.li(p, rm, cc),
        "fresnel_s" => x.fresnel_s(p, rm, cc),
        "fresnel_c" => x.fresnel_c(p, rm, cc),
        "elliptic_k" => x.elliptic_k(p, rm, cc),
        "elliptic_e_complete" => x.elliptic_e_complete(p, rm, cc),
        "lambert_w0" => x.lambert_w0(p, rm, cc),
        "lambert_wm1" => x.lambert_wm1(p, rm, cc),
        "kelvin_ber" => x.kelvin_ber(p, rm, cc),
        "kelvin_bei" => x.kelvin_bei(p, rm, cc),
        "kelvin_ker" => x.kelvin_ker(p, rm, cc),
        "kelvin_kei" => x.kelvin_kei(p, rm, cc),
        "scorer_gi" => x.scorer_gi(p, rm, cc),
        "scorer_hi" => x.scorer_hi(p, rm, cc),
        "clausen_cl2" => x.clausen_cl2(p, rm, cc),
        "clausen_cl3" => x.clausen_cl3(p, rm, cc),
        "riemann_zeta" => x.riemann_zeta(p, rm, cc),
        "barnes_g" => x.barnes_g(p, rm, cc),
        _ => return None,
    })
}

fn eval_real(
    op: &str,
    toks: &[Tok],
    row: &serde_json::Value,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Result<Canon, String> {
    if let Some(y) = eval_unary(op, num(toks, 0)?, p, rm, cc) {
        return Ok(canon_of(&y));
    }
    let y = match op {
        "gammainc" => num(toks, 0)?.gammainc(num(toks, 1)?, p, rm, cc),
        "gammainc_upper" => num(toks, 0)?.gammainc_upper(num(toks, 1)?, p, rm, cc),
        "bessel_j" => num(toks, 0)?.bessel_j(idx(toks, 1)?, p, rm, cc),
        "bessel_j_nu" => num(toks, 0)?.bessel_j_nu(num(toks, 1)?, p, rm, cc),
        "bessel_y" => num(toks, 0)?.bessel_y(num(toks, 1)?, p, rm, cc),
        "bessel_i" => num(toks, 0)?.bessel_i(num(toks, 1)?, p, rm, cc),
        "bessel_k" => num(toks, 0)?.bessel_k(num(toks, 1)?, p, rm, cc),
        "elliptic_f" => num(toks, 0)?.elliptic_f(num(toks, 1)?, p, rm, cc),
        "elliptic_e" => num(toks, 0)?.elliptic_e(num(toks, 1)?, p, rm, cc),
        "elliptic_f_phi" => num(toks, 0)?.elliptic_f_phi(num(toks, 1)?, p, rm, cc),
        "elliptic_e_phi" => num(toks, 0)?.elliptic_e_phi(num(toks, 1)?, p, rm, cc),
        "elliptic_pi" => num(toks, 0)?.elliptic_pi(num(toks, 1)?, num(toks, 2)?, p, rm, cc),
        "elliptic_pi_phi" => num(toks, 0)?.elliptic_pi_phi(num(toks, 1)?, num(toks, 2)?, p, rm, cc),
        "elliptic_pi_complete" => num(toks, 0)?.elliptic_pi_complete(num(toks, 1)?, p, rm, cc),
        "jacobi_am" => num(toks, 0)?.jacobi_am(num(toks, 1)?, p, rm, cc),
        "jacobi_sn" => num(toks, 0)?.jacobi_sn(num(toks, 1)?, p, rm, cc),
        "jacobi_cn" => num(toks, 0)?.jacobi_cn(num(toks, 1)?, p, rm, cc),
        "jacobi_dn" => num(toks, 0)?.jacobi_dn(num(toks, 1)?, p, rm, cc),
        "jacobi_cd" => num(toks, 0)?.jacobi_cd(num(toks, 1)?, p, rm, cc),
        "jacobi_dc" => num(toks, 0)?.jacobi_dc(num(toks, 1)?, p, rm, cc),
        "jacobi_arcsn" => num(toks, 0)?.jacobi_arcsn(num(toks, 1)?, p, rm, cc),
        "hypergeom_2f1" => {
            num(toks, 1)?.hypergeom_2f1(num(toks, 2)?, num(toks, 3)?, num(toks, 0)?, p, rm, cc)
        }
        "hypergeom_1f1" => num(toks, 0)?.hypergeom_1f1(num(toks, 1)?, num(toks, 2)?, p, rm, cc),
        "hypergeom_0f1" => num(toks, 0)?.hypergeom_0f1(num(toks, 1)?, p, rm, cc),
        "hurwitz_zeta" => num(toks, 0)?.hurwitz_zeta(num(toks, 1)?, p, rm, cc),
        "anger_j" => num(toks, 0)?.anger_j(num(toks, 1)?, p, rm, cc),
        "weber_e" => num(toks, 0)?.weber_e(num(toks, 1)?, p, rm, cc),
        "struve_h" => num(toks, 0)?.struve_h(num(toks, 1)?, p, rm, cc),
        "hermite_he" => num(toks, 0)?.hermite_he(idx(toks, 1)?, p, rm),
        "hermite_h" => num(toks, 0)?.hermite_h(idx(toks, 1)?, p, rm),
        "laguerre" => num(toks, 0)?.laguerre(idx(toks, 1)?, p, rm),
        "gen_laguerre" => num(toks, 0)?.gen_laguerre(idx(toks, 1)?, num(toks, 2)?, p, rm),
        "chebyshev_t" => num(toks, 0)?.chebyshev_t(idx(toks, 1)?, p, rm),
        "chebyshev_u" => num(toks, 0)?.chebyshev_u(idx(toks, 1)?, p, rm),
        "gegenbauer" => num(toks, 0)?.gegenbauer(idx(toks, 1)?, num(toks, 2)?, p, rm),
        "jacobi_p" => num(toks, 0)?.jacobi_p(idx(toks, 1)?, num(toks, 2)?, num(toks, 3)?, p, rm),
        "legendre_p" => {
            num(toks, 0)?.legendre_p(u32::try_from(idx(toks, 1)?).map_err(|_| "n")?, p, rm)
        }
        "polylog" => num(toks, 0)?.polylog(idx(toks, 1)?, p, rm, cc),
        "polygamma" => num(toks, 0)?.polygamma(idx(toks, 1)?, p, rm, cc),
        "betainc" => num(toks, 0)?.betainc(num(toks, 1)?, num(toks, 2)?, p, rm, cc),
        "normal_cdf" => num(toks, 0)?.normal_cdf(num(toks, 1)?, num(toks, 2)?, p, rm, cc),
        "log" => num(toks, 0)?.log(num(toks, 1)?, p, rm, cc),
        "nth_root" => num(toks, 0)?.nth_root(idx(toks, 1)?, p, rm),
        "meijer_g" => {
            let an = json_nums(row.get("an"), p, cc)?;
            let ap = json_nums(row.get("ap"), p, cc)?;
            let bm = json_nums(row.get("bm"), p, cc)?;
            let bq = json_nums(row.get("bq"), p, cc)?;
            num(toks, 0)?.meijer_g(&an, &ap, &bm, &bq, p, rm, cc)
        }
        "fox_h" => {
            let an = json_pairs(row.get("an"), p, cc)?;
            let ap = json_pairs(row.get("ap"), p, cc)?;
            let bm = json_pairs(row.get("bm"), p, cc)?;
            let bq = json_pairs(row.get("bq"), p, cc)?;
            num(toks, 0)?.fox_h(&an, &ap, &bm, &bq, p, rm, cc)
        }
        "heun_g" => num(toks, 0)?.heun_g(
            num(toks, 1)?,
            num(toks, 2)?,
            num(toks, 3)?,
            num(toks, 4)?,
            num(toks, 5)?,
            num(toks, 6)?,
            p,
            rm,
            cc,
        ),
        other => return Err(format!("unknown op {other}")),
    };
    Ok(canon_of(&y))
}

fn eval_complex(
    op: &str,
    toks: &[Tok],
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Result<Canon, String> {
    let z = cpx(toks, 0)?;
    let y = match op {
        "c.erf" => z.erf(p, rm, cc),
        "c.gamma" => z.gamma(p, rm, cc),
        "c.ln_gamma" => z.ln_gamma(p, rm, cc),
        "c.ei" => z.ei(p, rm, cc),
        "c.si" => z.si(p, rm, cc),
        "c.ci" => z.ci(p, rm, cc),
        "c.acosh" => z.acosh(p, rm, cc),
        "c.ai" => z.ai(p, rm, cc),
        "c.bi" => z.bi(p, rm, cc),
        "c.bessel_j" => z.bessel_j_nu(cpx(toks, 1)?, p, rm, cc),
        "c.bessel_y" => z.bessel_y(cpx(toks, 1)?, p, rm, cc),
        "c.bessel_i" => z.bessel_i(cpx(toks, 1)?, p, rm, cc),
        "c.bessel_k" => z.bessel_k(cpx(toks, 1)?, p, rm, cc),
        "c.hyp2f1" => z.hypergeom_2f1(cpx(toks, 1)?, cpx(toks, 2)?, cpx(toks, 3)?, p, rm, cc),
        other => return Err(format!("unknown op {other}")),
    };
    Ok(complex_of(&y))
}

fn eval_row(row: &serde_json::Value) -> Result<Canon, String> {
    let op = row.get("op").and_then(|v| v.as_str()).ok_or("missing op")?;
    let p = usize::try_from(row.get("p").and_then(|v| v.as_u64()).ok_or("missing p")?)
        .map_err(|_| "p")?;
    let rm = mode_of(row.get("rm").and_then(|v| v.as_str()).unwrap_or("ToEven"));
    let args: Vec<String> = row
        .get("args")
        .and_then(|v| v.as_array())
        .map(|arr| {
            arr.iter()
                .filter_map(|v| v.as_str().map(|s| s.to_string()))
                .collect()
        })
        .unwrap_or_default();
    let mut cc = Consts::new().map_err(|_| "constants".to_string())?;
    let mut toks = Vec::with_capacity(args.len());
    for arg in &args {
        toks.push(parse_tok(arg, p, &mut cc)?);
    }
    if op.starts_with("c.") {
        eval_complex(op, &toks, p, rm, &mut cc)
    } else {
        eval_real(op, &toks, row, p, rm, &mut cc)
    }
}

fn realize(c: &Canon) -> Result<ExactNum, String> {
    match c {
        Canon::Finite { sign, exp, limbs } => {
            let mut words = Vec::with_capacity(limbs.len().max(1));
            for limb in limbs {
                words.push(Word::from_str_radix(limb, 16).map_err(|_| format!("limb {limb}"))?);
            }
            if words.is_empty() {
                words.push(0);
            }
            let s = if sign == "-" { Sign::Neg } else { Sign::Pos };
            Ok(ExactNum::from_words(&words, s, *exp))
        }
        Canon::Inf(sign) => Ok(if sign == "-" { INF_NEG } else { INF_POS }),
        Canon::Nan(_) => Err("nan has no magnitude".to_string()),
        Canon::Complex { .. } => Err("complex part".to_string()),
    }
}

fn round_p(p: usize) -> usize {
    p.div_ceil(64) * 64
}

fn mag(x: &ExactNum) -> i32 {
    if x.is_zero() {
        i32::MIN
    } else {
        x.exponent().unwrap_or(i32::MIN)
    }
}

fn lost_bits(got: &ExactNum, want: &ExactNum, scale: &ExactNum, p: usize) -> Result<u32, String> {
    if got.is_nan() || want.is_nan() || scale.is_nan() {
        return Err("nan".to_string());
    }
    if got.is_inf() || want.is_inf() {
        return if got.is_inf() && want.is_inf() && got.sign() == want.sign() {
            Ok(0)
        } else {
            Err("inf".to_string())
        };
    }
    let width = round_p(p).max(64).saturating_mul(2);
    let diff = got.sub(want, width, RoundingMode::None);
    if diff.is_zero() {
        return Ok(0);
    }
    let correct =
        i64::from(scale.exponent().unwrap_or(0)) - i64::from(diff.exponent().unwrap_or(0));
    let word = i64::try_from(round_p(p).max(64)).unwrap_or(i64::MAX);
    if correct >= word {
        Ok(0)
    } else if correct <= 0 {
        u32::try_from(word).map_err(|_| "lost".to_string())
    } else {
        u32::try_from(word - correct).map_err(|_| "lost".to_string())
    }
}

fn lost_ok(got: &Canon, expect: &Canon, p: usize, max: u32) -> Result<(), String> {
    match (got, expect) {
        (Canon::Nan(g), Canon::Nan(e)) => {
            if g == e {
                Ok(())
            } else {
                Err("nan tag".to_string())
            }
        }
        (Canon::Complex { re: gr, im: gi }, Canon::Complex { re: er, im: ei }) => {
            let grn = realize(gr)?;
            let gin = realize(gi)?;
            let ern = realize(er)?;
            let ein = realize(ei)?;
            let scale = if mag(&ern) >= mag(&ein) { &ern } else { &ein };
            let lost = lost_bits(&grn, &ern, scale, p)?.max(lost_bits(&gin, &ein, scale, p)?);
            if lost <= max {
                Ok(())
            } else {
                Err(format!("lost {lost}>{max}"))
            }
        }
        _ => {
            let gn = realize(got)?;
            let en = realize(expect)?;
            let lost = lost_bits(&gn, &en, &en, p)?;
            if lost <= max {
                Ok(())
            } else {
                Err(format!("lost {lost}>{max}"))
            }
        }
    }
}

fn near_two_digits(got: &ExactNum, want: &ExactNum, scale: &ExactNum, cc: &mut Consts) -> bool {
    if want.is_inf() {
        return got.is_inf() && got.sign() == want.sign();
    }
    if got.is_nan() || got.is_inf() || want.is_nan() {
        return false;
    }
    if !want.is_zero() && got.sign() != want.sign() {
        return false;
    }
    let width = 256usize;
    let diff = got.sub(want, width, RoundingMode::None).abs();
    if diff.is_zero() {
        return true;
    }
    if scale.is_zero() || scale.is_nan() || scale.is_inf() {
        return false;
    }
    let tol_unit = ExactNum::parse("0.005", Radix::Dec, width, RoundingMode::ToEven, cc);
    let tol = scale.abs().mul(&tol_unit, width, RoundingMode::None);
    matches!(diff.compare(&tol), Some(c) if c <= 0)
}

fn smoke_ok(got: &Canon, expect: &Canon) -> Result<(), String> {
    let mut cc = Consts::new().map_err(|_| "constants".to_string())?;
    match (got, expect) {
        (Canon::Nan(g), Canon::Nan(e)) => {
            if e.is_none() || g == e {
                Ok(())
            } else {
                Err("nan tag".to_string())
            }
        }
        (Canon::Inf(g), Canon::Inf(e)) => {
            if g == e {
                Ok(())
            } else {
                Err("inf sign".to_string())
            }
        }
        (Canon::Finite { .. }, Canon::Finite { .. }) => {
            let gn = realize(got)?;
            let en = realize(expect)?;
            if near_two_digits(&gn, &en, &en, &mut cc) {
                Ok(())
            } else {
                Err("two digits".to_string())
            }
        }
        (Canon::Complex { re: gr, im: gi }, Canon::Complex { re: er, im: ei }) => {
            if matches!(er.as_ref(), Canon::Nan(_)) || matches!(ei.as_ref(), Canon::Nan(_)) {
                return if matches!(gr.as_ref(), Canon::Nan(_))
                    || matches!(gi.as_ref(), Canon::Nan(_))
                {
                    Ok(())
                } else {
                    Err("complex nan".to_string())
                };
            }
            let grn = realize(gr)?;
            let gin = realize(gi)?;
            let ern = realize(er)?;
            let ein = realize(ei)?;
            let scale = if mag(&ern) >= mag(&ein) { &ern } else { &ein };
            if near_two_digits(&grn, &ern, scale, &mut cc)
                && near_two_digits(&gin, &ein, scale, &mut cc)
            {
                Ok(())
            } else {
                Err("two digits".to_string())
            }
        }
        _ => Err("kind".to_string()),
    }
}

fn accuracy_path(op: &str) -> String {
    if let Some(rest) = op.strip_prefix("c.") {
        let name = match rest {
            "bessel_j" => "bessel_j_nu",
            "hyp2f1" => "hypergeom_2f1",
            other => other,
        };
        format!("ExactComplex::{name}")
    } else {
        format!("ExactNum::{op}")
    }
}

fn verdict(got: &Canon, expect: &Canon, row: &serde_json::Value) -> Result<(), String> {
    let p = usize::try_from(row.get("p").and_then(|v| v.as_u64()).unwrap_or(53)).unwrap_or(53);
    if let Some(max) = row.get("lost_max").and_then(|v| v.as_u64()) {
        let max = u32::try_from(max).unwrap_or(u32::MAX);
        return lost_ok(got, expect, p, max);
    }
    let op = row.get("op").and_then(|v| v.as_str()).unwrap_or("");
    let exact = row.get("exact").and_then(|v| v.as_bool()) == Some(true);
    let class = accuracy_of(&accuracy_path(op));
    if exact || matches!(class, Some(Accuracy::CorrectlyRounded)) {
        return if got == expect { Ok(()) } else { Err("limbs".to_string()) };
    }
    if let Some(Accuracy::WorkingPrecision { lost_bits }) = class {
        return lost_ok(got, expect, p, lost_bits);
    }
    // Unspecified, and functions outside the accuracy scan, use the two-digit smoke.
    if !matches!(
        class,
        Some(Accuracy::Unspecified) | None | Some(Accuracy::Faithful) | Some(Accuracy::Enclosure)
    ) {
        return Err(format!("unhandled class for {op}"));
    }
    smoke_ok(got, expect)
}

fn compress(lines: &[usize]) -> String {
    if lines.is_empty() {
        return String::new();
    }
    let mut out = String::new();
    let mut start = lines[0];
    let mut prev = lines[0];
    for &n in &lines[1..] {
        if n == prev + 1 {
            prev = n;
            continue;
        }
        push_range(&mut out, start, prev);
        start = n;
        prev = n;
    }
    push_range(&mut out, start, prev);
    out
}

fn push_range(out: &mut String, start: usize, end: usize) {
    if !out.is_empty() {
        out.push(',');
    }
    if start == end {
        out.push_str(&start.to_string());
    } else {
        out.push_str(&format!("{start}-{end}"));
    }
}

fn run_file(name: &str) {
    let path = format!(
        "{}/tests/oracle/goldens/specials/{name}.jsonl",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(&path).unwrap_or_else(|err| panic!("read {path}: {err}"));
    let mut ran = 0usize;
    let mut failed = 0usize;
    let mut samples: Vec<String> = Vec::new();
    let mut lines: Vec<usize> = Vec::new();
    let sample_cap = std::env::var("ORACLE_FAILS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(12usize);
    for (idx, line) in text.lines().enumerate() {
        if line.is_empty() {
            continue;
        }
        let row: serde_json::Value =
            serde_json::from_str(line).unwrap_or_else(|err| panic!("json {name}: {err}"));
        let tier = row.get("tier").and_then(|v| v.as_str()).unwrap_or("fast");
        if !tier_selected(tier) {
            continue;
        }
        ran += 1;
        let id = row
            .get("id")
            .and_then(|v| v.as_str())
            .unwrap_or("?")
            .to_string();
        let line_no = idx + 1;
        let budget = row
            .get("budget_s")
            .and_then(|v| v.as_u64())
            .map_or(Duration::from_secs(2), Duration::from_secs)
            // Budgets are for optimised builds; an unoptimised build that shares the machine
            // with the rest of the workspace tests runs several times slower.
            .saturating_mul(if cfg!(debug_assertions) { 10 } else { 1 });
        let expect = expect_of(row.get("expect").unwrap_or(&serde_json::Value::Null));
        let (tx, rx) = mpsc::channel();
        let row_for_thread = row.clone();
        thread::spawn(move || {
            let result = catch_unwind(AssertUnwindSafe(|| eval_row(&row_for_thread)));
            let _ = tx.send(result);
        });
        let outcome = match rx.recv_timeout(budget) {
            Ok(Ok(Ok(got))) => match verdict(&got, &expect, &row) {
                Ok(()) => None,
                Err(why) => Some(format!("{name}:{line_no} {id}: {why}")),
            },
            Ok(Ok(Err(err))) => Some(format!("{name}:{line_no} {id}: {err}")),
            Ok(Err(_)) => Some(format!("{name}:{line_no} {id}: panic")),
            Err(_) => Some(format!("{name}:{line_no} {id}: timeout")),
        };
        if let Some(msg) = outcome {
            failed += 1;
            lines.push(line_no);
            if samples.len() < sample_cap {
                samples.push(msg);
            }
        }
    }
    if failed != 0 {
        panic!(
            "{name}.jsonl: {failed} failed / {ran}\nlines {}\n{}",
            compress(&lines),
            samples.join("\n")
        );
    }
    eprintln!("{name}.jsonl passed {ran}");
}

#[test]
fn erf() {
    run_file("erf");
}

#[test]
fn gamma() {
    run_file("gamma");
}

#[test]
fn bessel() {
    run_file("bessel");
}

#[test]
fn elliptic() {
    run_file("elliptic");
}

#[test]
fn jacobi() {
    run_file("jacobi");
}

#[test]
fn airy() {
    run_file("airy");
}

#[test]
fn expint() {
    run_file("expint");
}

#[test]
fn hypergeom_2f1() {
    run_file("hypergeom_2f1");
}

#[test]
fn hypergeom_1f1() {
    run_file("hypergeom_1f1");
}

#[test]
fn hypergeom_0f1() {
    run_file("hypergeom_0f1");
}

#[test]
fn lambert() {
    run_file("lambert");
}

#[test]
fn kelvin() {
    run_file("kelvin");
}

#[test]
fn scorer() {
    run_file("scorer");
}

#[test]
fn clausen() {
    run_file("clausen");
}

#[test]
fn zeta() {
    run_file("zeta");
}

#[test]
fn barnes() {
    run_file("barnes");
}

#[test]
fn anger_weber_struve() {
    run_file("anger_weber_struve");
}

#[test]
fn orthopoly() {
    run_file("orthopoly");
}

#[test]
fn polylog() {
    run_file("polylog");
}

#[test]
fn meijer_fox_heun() {
    run_file("meijer_fox_heun");
}

#[test]
fn complex() {
    run_file("complex");
}

#[test]
fn other() {
    run_file("other");
}
