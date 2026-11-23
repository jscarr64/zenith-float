//! Tier 1 regression rows (Spence B2 on 11377ab, U18/T1/T2): one test per Tier 1 fix family, each row a case from
//! Harper's 5227f77 sweep (`~/pop-logs/harper/2026-11-13/zf-full-sweep`) that failed there and passes after the fix.
//!
//! `tier1/rows.tsv` columns: family, case id, function, precision, mode, judge, expected value, arguments. The expected
//! value is the correctly rounded result as judged by FLINT (python-flint 0.9.0, Arb) where Arb has the function, else
//! mpmath (correct rounding, no time limit), else exact arithmetic (next_after, logb, log base 0, IEEE, rationals);
//! the judge column names which. Inexact flags are not compared (the judges treat them as notes). The input and output
//! encodings follow Harper's sweep harness (`zf-sweep-harness/h2/src/main.rs`): `<hex>p<exp>` is hex·2^exp exactly,
//! `maxv:p` / `minpos:p` / `minnorm:p` are the extreme values at p bits, `re,im` is complex, IEEE operands are bits.
//!
//! Every test here fails on 5227f77 (the Tier 1 base) and passes on the land39 head; the red run is logged in
//! `~/pop-logs/paul/2026-11-15/land39/b2/`.

use zenith_float::{
    Consts, Error, ExactComplex, ExactInt, ExactNum, ExactRational, Ieee32, Ieee64, Radix, RawKind,
    RawParts, RoundingMode, Sign, INF_NEG, INF_POS, NAN, WORD_BIT_SIZE,
};

const ROWS: &str = include_str!("tier1/rows.tsv");

const MODES: [RoundingMode; 7] = [
    RoundingMode::Up,
    RoundingMode::Down,
    RoundingMode::ToZero,
    RoundingMode::FromZero,
    RoundingMode::ToEven,
    RoundingMode::ToOdd,
    RoundingMode::None,
];

fn hex_to_words(h: &str) -> Option<Vec<u64>> {
    let h = h.trim_start_matches('0');
    let mut words = Vec::new();
    let mut end = h.len();
    while end > 0 {
        let start = end.saturating_sub(16);
        words.push(u64::from_str_radix(h.get(start..end)?, 16).ok()?);
        end = start;
    }
    Some(words)
}

fn signed_zero(neg: bool) -> ExactNum {
    let mut z = ExactNum::new(64);
    if neg {
        z.set_sign(Sign::Neg);
    }
    z
}

fn parse_real(t: &str) -> Result<ExactNum, String> {
    let (neg, body) = match t.strip_prefix('-') {
        Some(b) => (true, b),
        None => (false, t.strip_prefix('+').unwrap_or(t)),
    };
    let sgn = |mut v: ExactNum| {
        if neg {
            v.set_sign(Sign::Neg);
        }
        v
    };
    match body {
        "0" => return Ok(signed_zero(neg)),
        "inf" => return Ok(if neg { INF_NEG } else { INF_POS }),
        "nan" => return Ok(NAN),
        _ => {}
    }
    for (pre, f) in [
        ("maxv:", ExactNum::max_value as fn(usize) -> ExactNum),
        ("minpos:", ExactNum::min_positive),
        ("minnorm:", ExactNum::min_positive_normal),
    ] {
        if let Some(w) = body.strip_prefix(pre) {
            return Ok(sgn(f(w.parse().map_err(|_| "bad p")?)));
        }
    }
    let (m, e) = body.split_once('p').ok_or("no p")?;
    let exp2: i64 = e.parse().map_err(|_| "bad exp")?;
    let mut words = hex_to_words(m).ok_or("bad hex")?;
    let Some(&top) = words.last() else {
        return Ok(signed_zero(neg));
    };
    let n = (words.len() - 1) * 64 + (64 - top.leading_zeros() as usize);
    let width = words.len() * 64;
    let s = u32::try_from(width - n).map_err(|_| "shift")?;
    if s > 0 {
        for i in (0..words.len()).rev() {
            let lo = if i > 0 { words[i - 1] >> (64 - s) } else { 0 };
            words[i] = (words[i] << s) | lo;
        }
    }
    let exponent = exp2
        .checked_add(i64::try_from(n).map_err(|_| "n")?)
        .ok_or("exp overflow")?;
    let parts = RawParts {
        kind: RawKind::Finite,
        words: &words,
        precision: width,
        negative: neg,
        exponent,
        inexact: false,
        op_inexact: false,
    };
    ExactNum::from_raw_parts(&parts).map_err(|e| format!("from_raw_parts: {e:?}"))
}

fn parse_complex(t: &str) -> Result<ExactComplex, String> {
    let (a, b) = t.split_once(',').ok_or("no comma")?;
    Ok(ExactComplex::new(parse_real(a)?, parse_real(b)?))
}

/// The value only (kind, sign, exponent, mantissa words): the harness format without its inexact-flag suffix.
fn fmt_real(x: &ExactNum) -> String {
    let rp = x.as_raw_parts();
    let s = if rp.negative { '-' } else { '+' };
    match rp.kind {
        RawKind::NaN => format!("N:{:?}", x.err()),
        RawKind::Inf => format!("I:{s}"),
        RawKind::Zero => format!("Z:{s}"),
        RawKind::Finite => {
            let h: String = rp.words.iter().rev().map(|w| format!("{w:016x}")).collect();
            format!("F:{s}:{}:{h}", rp.exponent)
        }
    }
}

fn fmt_c(z: &ExactComplex) -> String {
    format!("{};{}", fmt_real(z.re()), fmt_real(z.im()))
}

fn mode(s: &str) -> Result<RoundingMode, String> {
    Ok(match s {
        "none" => RoundingMode::None,
        "up" => RoundingMode::Up,
        "down" => RoundingMode::Down,
        "tozero" => RoundingMode::ToZero,
        "fromzero" => RoundingMode::FromZero,
        "toeven" => RoundingMode::ToEven,
        "toodd" => RoundingMode::ToOdd,
        _ => return Err(format!("mode {s}")),
    })
}

fn arg<'a>(a: &[&'a str], i: usize) -> Result<&'a str, String> {
    a.get(i).copied().ok_or_else(|| format!("missing arg {i}"))
}

fn real(a: &[&str], i: usize) -> Result<ExactNum, String> {
    parse_real(arg(a, i)?)
}

fn int(s: &str) -> Result<ExactInt, String> {
    ExactInt::from_dec_string(s).ok_or_else(|| format!("bad int {s}"))
}

/// `r.meijer_g` / `r.fox_h` rows: z, |an|, |ap|, |bm|, |bq|, then the values
/// (fox_h: (value, scale) pairs).
fn meijer_fox(
    func: &str,
    p: usize,
    rm: RoundingMode,
    a: &[&str],
    cc: &mut Consts,
) -> Result<ExactNum, String> {
    let fox = func == "r.fox_h";
    let z = real(a, 0)?;
    let mut k = 5;
    let mut groups: Vec<Vec<(ExactNum, ExactNum)>> = Vec::new();
    for i in 1..5 {
        let c: usize = arg(a, i)?.parse().map_err(|_| "count")?;
        let mut g = Vec::new();
        for _ in 0..c {
            let s = if fox { real(a, k + 1)? } else { ExactNum::from_u8(1, 64) };
            g.push((real(a, k)?, s));
            k += if fox { 2 } else { 1 };
        }
        groups.push(g);
    }
    Ok(if fox {
        z.fox_h(&groups[0], &groups[1], &groups[2], &groups[3], p, rm, cc)
    } else {
        let v: Vec<Vec<ExactNum>> = groups
            .into_iter()
            .map(|g| g.into_iter().map(|x| x.0).collect())
            .collect();
        z.meijer_g(&v[0], &v[1], &v[2], &v[3], p, rm, cc)
    })
}

fn run(
    func: &str,
    p: usize,
    rm: RoundingMode,
    a: &[&str],
    cc: &mut Consts,
) -> Result<String, String> {
    let pair = |(x, y): (ExactNum, ExactNum)| format!("{};{}", fmt_real(&x), fmt_real(&y));
    let ieee64 = |i: usize| -> Result<Ieee64, String> {
        Ok(Ieee64::from_bits(
            u64::from_str_radix(arg(a, i)?, 16).map_err(|_| "bits")?,
        ))
    };
    let ieee32 = |i: usize| -> Result<Ieee32, String> {
        Ok(Ieee32::from_bits(
            u32::from_str_radix(arg(a, i)?, 16).map_err(|_| "bits")?,
        ))
    };
    let out = match func {
        "r.barnes_g" => real(a, 0)?.barnes_g(p, rm, cc),
        "r.bessel_y" => real(a, 0)?.bessel_y(&real(a, 1)?, p, rm, cc),
        "r.bi" => real(a, 0)?.bi(p, rm, cc),
        "r.bi_prime" => real(a, 0)?.bi_prime(p, rm, cc),
        "r.ai" => real(a, 0)?.ai(p, rm, cc),
        "r.ai_prime" => real(a, 0)?.ai_prime(p, rm, cc),
        "r.digamma" => real(a, 0)?.digamma(p, rm, cc),
        "r.gamma" => real(a, 0)?.gamma(p, rm, cc),
        "r.ln_gamma" => real(a, 0)?.ln_gamma(p, rm, cc),
        "r.riemann_zeta" => real(a, 0)?.riemann_zeta(p, rm, cc),
        "r.scorer_gi" => real(a, 0)?.scorer_gi(p, rm, cc),
        "r.scorer_hi" => real(a, 0)?.scorer_hi(p, rm, cc),
        "r.two_sum" => return Ok(pair(real(a, 0)?.two_sum(&real(a, 1)?, p, rm))),
        "r.two_product" => return Ok(pair(real(a, 0)?.two_product(&real(a, 1)?, p, rm))),
        "r.sin_cos" => return Ok(pair(real(a, 0)?.sin_cos(p, rm, cc))),
        "r.elliptic_pi" => real(a, 0)?.elliptic_pi(&real(a, 1)?, &real(a, 2)?, p, rm, cc),
        "r.gammainc_upper" => real(a, 0)?.gammainc_upper(&real(a, 1)?, p, rm, cc),
        "r.gammainc" => real(a, 0)?.gammainc(&real(a, 1)?, p, rm, cc),
        "r.betainc" => real(a, 0)?.betainc(&real(a, 1)?, &real(a, 2)?, p, rm, cc),
        "r.bessel_i" => real(a, 0)?.bessel_i(&real(a, 1)?, p, rm, cc),
        "r.bessel_j_nu" => real(a, 0)?.bessel_j_nu(&real(a, 1)?, p, rm, cc),
        "r.bessel_k" => real(a, 0)?.bessel_k(&real(a, 1)?, p, rm, cc),
        "r.beta" => zenith_float::beta(&real(a, 0)?, &real(a, 1)?, p, rm, cc),
        "r.hypergeom_2f1" => {
            real(a, 0)?.hypergeom_2f1(&real(a, 1)?, &real(a, 2)?, &real(a, 3)?, p, rm, cc)
        }
        "r.jacobi_am" => real(a, 0)?.jacobi_am(&real(a, 1)?, p, rm, cc),
        "r.jacobi_sn" => real(a, 0)?.jacobi_sn(&real(a, 1)?, p, rm, cc),
        "r.jacobi_ns" => real(a, 0)?.jacobi_ns(&real(a, 1)?, p, rm, cc),
        "r.jacobi_nc" => real(a, 0)?.jacobi_nc(&real(a, 1)?, p, rm, cc),
        "r.jacobi_nd" => real(a, 0)?.jacobi_nd(&real(a, 1)?, p, rm, cc),
        "r.jacobi_sc" => real(a, 0)?.jacobi_sc(&real(a, 1)?, p, rm, cc),
        "r.jacobi_sd" => real(a, 0)?.jacobi_sd(&real(a, 1)?, p, rm, cc),
        "r.jacobi_arccn" => real(a, 0)?.jacobi_arccn(&real(a, 1)?, p, rm, cc),
        "r.jacobi_arcdn" => real(a, 0)?.jacobi_arcdn(&real(a, 1)?, p, rm, cc),
        "r.elliptic_f_phi" => real(a, 0)?.elliptic_f_phi(&real(a, 1)?, p, rm, cc),
        "r.elliptic_e_phi" => real(a, 0)?.elliptic_e_phi(&real(a, 1)?, p, rm, cc),
        "r.elliptic_pi_phi" => real(a, 0)?.elliptic_pi_phi(&real(a, 1)?, &real(a, 2)?, p, rm, cc),
        "r.meijer_g" | "r.fox_h" => meijer_fox(func, p, rm, a, cc)?,
        "r.log" => real(a, 0)?.log(&real(a, 1)?, p, rm, cc),
        "r.logb" => real(a, 0)?.logb(p, rm),
        "r.rem" => real(a, 0)?.rem(&real(a, 1)?),
        "r.next_after" => real(a, 0)?.next_after(&real(a, 1)?, p, rm),
        "r.atan" => real(a, 0)?.atan(p, rm, cc),
        "r.atan2" => real(a, 0)?.atan2(&real(a, 1)?, p, rm, cc),
        "r.sin" => real(a, 0)?.sin(p, rm, cc),
        "r.cos" => real(a, 0)?.cos(p, rm, cc),
        "r.tan" => real(a, 0)?.tan(p, rm, cc),
        "r.clausen_cl2" => real(a, 0)?.clausen_cl2(p, rm, cc),
        "r.clausen_cl3" => real(a, 0)?.clausen_cl3(p, rm, cc),
        "r.normal_cdf" => real(a, 0)?.normal_cdf(&real(a, 1)?, &real(a, 2)?, p, rm, cc),
        "r.normal_pdf" => real(a, 0)?.normal_pdf(&real(a, 1)?, &real(a, 2)?, p, rm, cc),
        "r.chi_squared_cdf" => real(a, 0)?.chi_squared_cdf(&real(a, 1)?, p, rm, cc),
        "r.poisson_pmf" => real(a, 0)?.poisson_pmf(&real(a, 1)?, p, rm, cc),
        "r.student_t_pdf" => real(a, 0)?.student_t_pdf(&real(a, 1)?, p, rm, cc),
        "c.exp" => return Ok(fmt_c(&parse_complex(arg(a, 0)?)?.exp(p, rm, cc))),
        "c.sin" => return Ok(fmt_c(&parse_complex(arg(a, 0)?)?.sin(p, rm, cc))),
        "c.cos" => return Ok(fmt_c(&parse_complex(arg(a, 0)?)?.cos(p, rm, cc))),
        "c.atan" => return Ok(fmt_c(&parse_complex(arg(a, 0)?)?.atan(p, rm, cc))),
        "c.ai" => return Ok(fmt_c(&parse_complex(arg(a, 0)?)?.ai(p, rm, cc))),
        "c.bi" => return Ok(fmt_c(&parse_complex(arg(a, 0)?)?.bi(p, rm, cc))),
        "c.mul" => {
            return Ok(fmt_c(&parse_complex(arg(a, 0)?)?.mul(
                &parse_complex(arg(a, 1)?)?,
                p,
                rm,
            )))
        }
        "c.arg" => parse_complex(arg(a, 0)?)?.arg(p, rm, cc),
        "misc.from_i128" => ExactNum::from_i128(arg(a, 0)?.parse::<i128>().map_err(|_| "i128")?, p),
        "x.rat" => {
            let x = ExactRational::from_ints(int(arg(a, 1)?)?, int(arg(a, 2)?)?);
            let r = match arg(a, 0)? {
                "round" => x.round(),
                "floor" => x.floor(),
                "ceil" => x.ceil(),
                op => return Err(format!("rat op {op}")),
            };
            return Ok(match r.format_exact(Radix::Dec) {
                Some(s) => format!("S:{s}"),
                None => String::from(if r.is_nan() { "NAN" } else { "NOFMT" }),
            });
        }
        _ => {
            let (kind, f) = func.split_once('.').ok_or("func")?;
            return match kind {
                "ieee64" => {
                    let r = match f {
                        "add" => ieee64(0)?.soft_add(ieee64(1)?),
                        "sub" => ieee64(0)?.soft_sub(ieee64(1)?),
                        "mul" => ieee64(0)?.soft_mul(ieee64(1)?),
                        "div" => ieee64(0)?.soft_div(ieee64(1)?),
                        "fma" => ieee64(0)?.mul_add(ieee64(1)?, ieee64(2)?),
                        "from_exact" => Ieee64::from_exact(&real(a, 0)?),
                        _ => return Err(format!("ieee64 {f}")),
                    };
                    Ok(format!("B:{:016x}", r.to_bits()))
                }
                "ieee32" => {
                    let r = match f {
                        "add" => ieee32(0)?.soft_add(ieee32(1)?),
                        "sub" => ieee32(0)?.soft_sub(ieee32(1)?),
                        "mul" => ieee32(0)?.soft_mul(ieee32(1)?),
                        "div" => ieee32(0)?.soft_div(ieee32(1)?),
                        "fma" => ieee32(0)?.mul_add(ieee32(1)?, ieee32(2)?),
                        "from_exact" => Ieee32::from_exact(&real(a, 0)?),
                        _ => return Err(format!("ieee32 {f}")),
                    };
                    Ok(format!("B:{:08x}", r.to_bits()))
                }
                _ => Err(format!("no dispatch for {func}")),
            };
        }
    };
    Ok(fmt_real(&out))
}

/// Runs every row of `family`; fails listing each row whose value differs from the judged expectation.
fn family(name: &str) {
    let mut cc = Consts::new().expect("Consts");
    let (mut n, mut bad) = (0usize, Vec::new());
    for line in ROWS.lines().filter(|l| l.split('\t').next() == Some(name)) {
        let c: Vec<&str> = line.split('\t').collect();
        let (id, func, judge, want) = (c[1], c[2], c[5], c[6]);
        let p: usize = c[3].parse().expect("precision");
        let rm = mode(c[4]).expect("mode");
        let got = run(func, p, rm, &c[7..], &mut cc).unwrap_or_else(|e| format!("HARNESS-ERR:{e}"));
        n += 1;
        if got != want {
            bad.push(format!(
                "{id} {func} p={p} {} args={:?} [{judge}]\n    want {want}\n    got  {got}",
                c[4],
                &c[7..]
            ));
        }
    }
    assert!(n > 0, "no rows for {name}");
    assert!(
        bad.is_empty(),
        "{name}: {} of {n} rows wrong:\n{}",
        bad.len(),
        bad.join("\n")
    );
}

macro_rules! families {
    ($($test:ident => $fam:literal),* $(,)?) => {
        $(#[test] fn $test() { family($fam); })*
    };
}

families! {
    barnes_g_negint => "T1-barnes_g-negint",
    barnes_g_tiny => "T1-barnes_g-tiny",
    bessel_y_x0 => "T1-bessel_y-x0",
    c_atan_branchcut => "T1-c.atan-branchcut",
    c_mul_double_round => "T1-c.mul-double-round",
    cexp_csin_ccos_no_ziv => "T1-cexp-csin-ccos-no-ziv",
    directed_overflow => "T1-directed-overflow",
    elliptic_pi_x1 => "T1-elliptic_pi-x1",
    from_i128 => "T1-from_i128",
    gamma_stirling_cap => "T1-gamma-stirling-cap",
    gammainc_upper_s0 => "T1-gammainc_upper-s0",
    ieee_soft_rounding => "T1-ieee_soft-rounding",
    jacobi_m1_huge => "T1-jacobi-m1-huge",
    jacobi_period_bits => "T1-jacobi-period-bits",
    jacobi_am_quasiperiod => "T1-jacobi_am-quasiperiod",
    jacobi_arccn_sign => "T1-jacobi_arccn-sign",
    log_base0 => "T1-log-base0",
    logb_inexact => "T1-logb-inexact",
    next_after => "T1-next_after",
    normal_sigma0 => "T1-normal-sigma0",
    pi_directed_sign => "T1-pi-directed-sign",
    rational_round => "T1-rational-round",
    rem_overflow => "T1-rem-overflow",
    trig_huge_arg => "T1-trig-huge-arg",
}

/// Exact fmod, no rounding allowance: each finite binary float is an integer multiple of minpos at its precision.
/// Across the whole exponent range MAX rem minpos must be signed zero. Odd multiples also exercise every exponent
/// bit of the modular power (expected residues from Python integer pow, independent of the floating library).
#[test]
fn rem_overflow_exact_top_bit_and_odd_divisors() {
    for p in [64usize, 128, 256, 1024, 2048] {
        for neg in [false, true] {
            let mut x = ExactNum::max_value(p);
            x.set_sign(if neg { Sign::Neg } else { Sign::Pos });
            let r = x.rem(&ExactNum::min_positive(p));
            assert!(
                r.is_zero(),
                "MAX rem minpos p={p} neg={neg}: {}",
                fmt_real(&r)
            );
            assert_eq!(r.sign(), Some(if neg { Sign::Neg } else { Sign::Pos }));
        }
    }
    for (p, modulus, residue) in
        [(64, 7, 1), (64, 11, 3), (64, 13, 11), (128, 7, 3), (128, 11, 7), (128, 13, 5)]
    {
        let tiny = ExactNum::min_positive(p);
        let y = ExactNum::mul_ref(
            &tiny,
            &ExactNum::from_i64(modulus, p),
            p,
            RoundingMode::ToEven,
        );
        let want = ExactNum::mul_ref(
            &tiny,
            &ExactNum::from_i64(residue, p),
            p,
            RoundingMode::ToEven,
        );
        let got = ExactNum::max_value(p).rem(&y);
        // Value equality: rem keeps the dividend's width, so the mantissa words differ in length.
        assert_eq!(
            got.compare(&want),
            Some(0),
            "MAX{p} rem {modulus} minpos: {}",
            fmt_real(&got)
        );
        assert_eq!(got.sign(), Some(Sign::Pos));
    }
}

/// FLINT/Arb Ziv references (`b2/trig_refs.py`) for sin(2^k) at p = 64. exponent(2^16) = 2^16 + 1 is one bit past
/// main's silent 2^16 cap (which returned a value with no correct bits); 2^18 is further out. Tier 1's 2^22 refusal
/// boundary (exponent 2^22 + 1) is covered by `reduction_width` in the num unit tests: one 4-million-bit π takes
/// about 15 minutes in this library, a speed issue (Tier 2), not a limit.
#[test]
#[ignore = "65536- and 262144-bit pi; run by the full gate's ignored (release) step"]
fn trig_no_cap_past_old_limit_values() -> Result<(), String> {
    let mut cc = Consts::new().map_err(|e| format!("Consts {e:?}"))?;
    let refs = [
        (
            65536i64,
            "F:-:0:e7e8101e7aad1954",
            "F:-:0:e7e8101e7aad1953",
            "F:-:0:e7e8101e7aad1954",
        ),
        (
            262144,
            "F:+:0:c7ef2df578b3285c",
            "F:+:0:c7ef2df578b3285c",
            "F:+:0:c7ef2df578b3285b",
        ),
    ];
    for (k, even, up, down) in refs {
        let x = ExactNum::pow2(k, 64).map_err(|e| format!("pow2 {e:?}"))?;
        for (rm, want) in
            [(RoundingMode::ToEven, even), (RoundingMode::Up, up), (RoundingMode::Down, down)]
        {
            assert_eq!(fmt_real(&x.sin(64, rm, &mut cc)), want, "sin(2^{k}) {rm:?}");
        }
    }
    Ok(())
}

/// Every public precision up to MAX_PRECISION is valid. When there is no memory for its working width the result
/// is MemoryAllocation / NotConverged / None, never a panic or InvalidArgument just for being near MAX_PRECISION.
#[test]
fn precision_limit_panics_valid_widths() {
    use zenith_float::{
        gauss_kronrod_interval, wynn_epsilon, Error, QuadratureError, MAX_PRECISION, WORD_BIT_SIZE,
    };
    let mut cc = Consts::new().expect("Consts");
    let one = ExactNum::from_i64(1, 64);
    let zero = ExactNum::new(64);
    let z = ExactComplex::new(ExactNum::from_i64(1, 64), ExactNum::from_i64(1, 64));
    for p in [
        MAX_PRECISION,
        MAX_PRECISION - 1,
        MAX_PRECISION - WORD_BIT_SIZE,
        MAX_PRECISION - 2 * WORD_BIT_SIZE,
    ] {
        let v = cc.ln_10(p, RoundingMode::ToEven);
        assert_eq!(v.err(), Some(Error::MemoryAllocation), "ln10 valid p={p}");
        // The complex functions can report MemoryAllocation but must not panic.
        for v in [
            z.log10(p, RoundingMode::ToEven, &mut cc),
            z.exp10(p, RoundingMode::ToEven, &mut cc),
        ] {
            assert!(v.re().is_nan() && v.im().is_nan(), "p={p}");
        }
        let q = gauss_kronrod_interval(
            |x, p, rm, _| x.sqr(p, rm),
            &zero,
            &one,
            p,
            RoundingMode::ToEven,
            &mut cc,
        );
        assert!(
            !matches!(q, Err(QuadratureError::InvalidArgument)),
            "valid GK precision {p}"
        );
        // A two-term sequence needs no table arithmetic at the working width: a finite value or None, no panic.
        let w = wynn_epsilon(
            &[ExactNum::from_i64(1, 64), ExactNum::from_i64(2, 64)],
            p,
            RoundingMode::ToEven,
        );
        assert!(w.as_ref().is_none_or(|v| !v.is_nan()), "wynn p={p}");
    }
}

/// T2 boundary matrix: all seven modes and public precisions (32, 64, 128, 256, 1024), signed zeros, infinities,
/// NaN and the subnormal floor. Exact next_after, from_i128 and normal(sigma=0) expectations need no transcendental
/// reference; log base ±0 is ln(x)/ln(0), per IEEE division. Broader family-specific values live in rows.tsv.
#[test]
fn tier1_boundary_matrix() {
    let modes = [
        RoundingMode::None,
        RoundingMode::Up,
        RoundingMode::Down,
        RoundingMode::ToZero,
        RoundingMode::FromZero,
        RoundingMode::ToEven,
        RoundingMode::ToOdd,
    ];
    let mut cc = Consts::new().expect("Consts");
    for p in [32, 64, 128, 256, 1024] {
        let tiny = ExactNum::min_positive(p);
        // Inputs are 64-bit values (as in the sweep harness); p is the result precision.
        let one = parse_real("1p0").expect("one");
        for rm in modes {
            for neg in [false, true] {
                let mut zero = ExactNum::new(p);
                zero.set_sign(if neg { Sign::Neg } else { Sign::Pos });
                let next = zero.next_after(&INF_POS, p, rm);
                assert_eq!(
                    next.compare(&tiny),
                    Some(0),
                    "next_after zero -> inf p={p} {rm:?}"
                );
                assert_eq!(next.sign(), Some(Sign::Pos));
                let log = one.log(&zero, p, rm, &mut cc);
                assert!(log.is_zero(), "log_0(1) p={p} {rm:?}: {}", fmt_real(&log));
                assert_eq!(log.sign(), Some(Sign::Neg));
                for v in [
                    one.normal_pdf(&zero, &zero, p, rm, &mut cc),
                    one.normal_cdf(&zero, &zero, p, rm, &mut cc),
                ] {
                    assert!(v.is_nan(), "sigma=0 p={p} {rm:?}");
                }
            }
            assert!(NAN.next_after(&one, p, rm).is_nan());
            assert!(
                INF_POS
                    .next_after(&INF_NEG, p, rm)
                    .compare(&ExactNum::max_value(p))
                    == Some(0)
            );
            let n = ExactNum::from_i128(17, p);
            assert_eq!(
                n.compare(&parse_real("11p0").expect("17")),
                Some(0),
                "from_i128(17) p={p}"
            );
        }
    }
}

// ---- Harper ZF bug regressions (Jeff 2026-11-15 11:14 CT: every ZF fixed for every input in its rows, or BROKEN).

/// ZF-05: a decimal numeral below the least subnormal rounds by mode (Up away from zero for +, Down for -): the
/// least subnormal or a signed zero, never a silent +0 in a direction that must not reach zero.
#[test]
fn zf_05_parse_below_least_subnormal_rounds_by_mode() {
    use zenith_float::Radix;
    let mut cc = Consts::new().expect("Consts");
    for p in [64usize, 128, 256, 1024, 2048] {
        for s in ["1e-4611686018427387904", "1e-4611686018427387903", "7e-9000000000000000000"] {
            for neg in [false, true] {
                let txt = if neg { format!("-{s}") } else { s.to_string() };
                for rm in MODES {
                    let v = ExactNum::parse(&txt, Radix::Dec, p, rm, &mut cc);
                    let away = match rm {
                        RoundingMode::Up => !neg,
                        RoundingMode::Down => neg,
                        RoundingMode::FromZero => true,
                        _ => false,
                    };
                    if away {
                        let mut want = ExactNum::min_positive(p);
                        if neg {
                            want.set_sign(Sign::Neg);
                        }
                        assert_eq!(
                            v.compare(&want),
                            Some(0),
                            "{txt} p={p} {rm:?}: {}",
                            fmt_real(&v)
                        );
                    } else {
                        assert!(v.is_zero(), "{txt} p={p} {rm:?}: {}", fmt_real(&v));
                        assert_eq!(v.is_negative(), neg, "{txt} p={p} {rm:?} sign");
                    }
                }
            }
        }
    }
}

/// ZF-71: x^-k with k > 2^62 underflows (it is not an overflow): 2^-(2^62 + 1000) is below every subnormal, so it
/// rounds by mode to a signed zero or the least subnormal; i64::MIN and powers of a non-power-of-two base too.
#[test]
fn zf_71_powsi_huge_negative_exponent_underflows_by_mode() {
    let n = -((1i64 << 62) + 1000);
    for p in [64usize, 256, 1024] {
        for base in [2i64, 3, -2, -3] {
            for k in [n, n - 1, i64::MIN, i64::MIN + 1] {
                if base.abs() == 2 && p == 1024 && k > i64::MIN + 1 {
                    continue; // 2^-(2^62 + 1000) is a 1024-bit subnormal: checked exactly below
                }
                let neg = base < 0 && k & 1 == 1;
                for rm in MODES {
                    let v = ExactNum::from_i64(base, p).powsi(k, p, rm);
                    let away = match rm {
                        RoundingMode::Up => !neg,
                        RoundingMode::Down => neg,
                        RoundingMode::FromZero => true,
                        _ => false,
                    };
                    let tag = format!("{base}^{k} p={p} {rm:?}: {}", fmt_real(&v));
                    if away {
                        let mut want = ExactNum::min_positive(p);
                        if neg {
                            want.set_sign(Sign::Neg);
                        }
                        assert_eq!(v.compare(&want), Some(0), "{tag}");
                    } else {
                        assert!(v.is_zero(), "{tag}");
                        assert_eq!(v.is_negative(), neg, "{tag} sign");
                    }
                }
            }
        }
    }
    // 2^-(2^62 + 1000) = 2^(EXPONENT_MIN - 1001) and its half are exact 1024-bit subnormals: minpos * 2^23, 2^22.
    for (k, up) in [(n, 23i64), (n - 1, 22)] {
        for rm in MODES {
            let v = ExactNum::from_i64(2, 1024).powsi(k, 1024, rm);
            let want = ExactNum::min_positive(1024).ldexp(up, 1024, RoundingMode::None);
            assert_eq!(
                v.compare(&want),
                Some(0),
                "2^{k} p=1024 {rm:?}: {}",
                fmt_real(&v)
            );
        }
    }
    // Still in range: 2^-(2^61) is a normal value, exactly.
    let v = ExactNum::from_i64(2, 64).powsi(-(1i64 << 61), 64, RoundingMode::ToEven);
    assert_eq!(v.exponent(), Some(1 - (1i64 << 61)), "{}", fmt_real(&v));
}

/// ZF-06: format(±MAX) is a finite numeral in every mode and radix, never "Inf": a directed mode that rounds the
/// decimal past MAX gets the other spelling, and the text parses back to ±MAX in the same mode and in ToEven.
#[test]
fn zf_06_format_max_is_finite_and_round_trips() {
    let mut cc = Consts::new().expect("Consts");
    for p in [64usize, 128, 256, 1024, 2048] {
        for neg in [false, true] {
            let mut x = ExactNum::max_value(p);
            if neg {
                x.set_sign(Sign::Neg);
            }
            for rdx in [Radix::Bin, Radix::Oct, Radix::Dec, Radix::Hex] {
                for rm in MODES {
                    let s = x.format(rdx, rm, &mut cc).expect("format");
                    assert!(
                        !s.contains("Inf") && s != "Err",
                        "p={p} neg={neg} {rdx:?} {rm:?}: {s}"
                    );
                    for prm in [rm, RoundingMode::ToEven] {
                        let y = ExactNum::parse(&s, rdx, p, prm, &mut cc);
                        assert_eq!(
                            y.compare(&x),
                            Some(0),
                            "p={p} {rdx:?} {rm:?} parse {prm:?} {s}: {}",
                            fmt_real(&y)
                        );
                    }
                }
            }
        }
    }
}

/// ZF-76: on the imaginary axis erf(±0 + iy) = ±0 + i·erfi(y) and erfc(±0 + iy) = 1 − i·erfi(y): the real parts
/// are exact (a signed zero carrying the sign of Re z, and exactly 1), in every mode; the imaginary part is
/// erfi(y) rounded by mode (FLINT acb erf, mpmath erfi agree to 130 digits).
#[test]
fn zf_76_c_erf_imaginary_axis_real_part_exact() {
    let mut cc = Consts::new().expect("Consts");
    let erfi: [(&str, &str); 5] = [
        ("5p1", "1.524307422708669699360546614726544062463812074698857310557201311229700714271763759939487330097258173455013441042776677333802311927e+42"),
        ("3p0", "1.629994622601565651061647952076274162778889910915259988647982150695656826775680269243153186081103257149596077135803271842564204772e+3"),
        ("1p-1", "6.149520946965109808396811856236413930513456178954035723266386934995887895613341957110356633853898031990244274675974500382667930578e-1"),
        ("fp1", "1.379047205635398320143951493011934592102794303811758758732511651393839606187768538637067194569236345548044800918610209213444574807e+389"),
        ("1p-4", "7.061563338595019687242633333808628245633743337566627759223830419468093689378912767943556359155055522339060979687785301946058337896e-2"),
    ];
    for p in [64usize, 256] {
        for (y, val) in erfi {
            for yneg in [false, true] {
                for reneg in [false, true] {
                    let ys = if yneg { format!("-{y}") } else { y.to_string() };
                    let rs = if reneg { "-0" } else { "0" };
                    let z = parse_complex(&format!("{rs},{ys}")).expect("z");
                    for rm in MODES {
                        if rm == RoundingMode::None {
                            continue;
                        }
                        let w = z.erf(p, rm, &mut cc);
                        let txt = if yneg { format!("-{val}") } else { val.to_string() };
                        let want = ExactNum::parse(&txt, Radix::Dec, p, rm, &mut cc);
                        let tag = format!("erf({rs}+i{ys}) p={p} {rm:?}: {}", fmt_c(&w));
                        assert!(w.re().is_zero() && w.re().is_negative() == reneg, "{tag}");
                        assert_eq!(w.im().compare(&want), Some(0), "{tag}");
                        let c = z.erfc(p, rm, &mut cc);
                        let ntxt = if yneg { val.to_string() } else { format!("-{val}") };
                        let nwant = ExactNum::parse(&ntxt, Radix::Dec, p, rm, &mut cc);
                        let tag = format!("erfc({rs}+i{ys}) p={p} {rm:?}: {}", fmt_c(&c));
                        assert_eq!(c.re().compare(&ExactNum::from_u8(1, p)), Some(0), "{tag}");
                        assert_eq!(c.im().compare(&nwant), Some(0), "{tag}");
                    }
                }
            }
        }
    }
}

/// ZF-10: gcd over integer polynomials must not use rounded remainders. FLINT fmpz_poly gives x−1 for
/// (x−1)(x−2)(x−3)(x−4) and x⁴−1, and for Harper's x⁴−10x³+23x²−38x+24 and x⁴−1: division by a nonunit leading coefficient creates rational
/// remainders, whose rounding used to turn the exact common factor into a nonzero constant.
#[test]
fn zf_10_poly_gcd_integer_remainders_exact() {
    use zenith_float::ExactNumPoly;
    for p in [64usize, 128, 256] {
        for rm in MODES {
            let poly = |v: &[i64]| {
                ExactNumPoly::from_coeffs(
                    p,
                    rm,
                    &v.iter()
                        .map(|&x| ExactNum::from_i64(x, p))
                        .collect::<Vec<_>>(),
                )
            };
            for scale in [1i64, -7, 13] {
                let a = poly(&[24 * scale, -50 * scale, 35 * scale, -10 * scale, scale]);
                let h = poly(&[24 * scale, -38 * scale, 23 * scale, -10 * scale, scale]); // Harper d0228098
                let b = poly(&[-1, 0, 0, 0, 1]);
                for (a, b) in [(&a, &b), (&b, &a), (&h, &b), (&b, &h)] {
                    let g = a.gcd(b, p, rm);
                    assert_eq!(g.degree(), Some(1), "p={p} {rm:?} scale={scale}");
                    assert_eq!(g.coeff(0).compare(&ExactNum::from_i8(-1, p)), Some(0));
                    assert_eq!(g.coeff(1).compare(&ExactNum::from_u8(1, p)), Some(0));
                }
            }
        }
    }
}

/// ZF-47: Ball::exp encloses exp(m−r)..exp(m+r) even when `rm` is None: the radius cannot be based
/// on a rounded-down exp(m) without accounting for the *radius* error that this introduces. FLINT
/// arb exp, mpmath exp endpoint references at 2600 bits (720 decimal digits).
#[test]
fn zf_47_ball_exp_encloses_endpoints() {
    use zenith_float::Ball;
    let mut cc = Consts::new().expect("Consts");
    let refs: [(&str, &str, &str, &str); 7] = [
        ("-2", "5", "9.11881965554516208003136084409282626473724527436053840816133421889479889310306529326279006464900442087462481380818271013199660566979106506670204549524834385945835040211210998958128337025691930041681590202201853633082331489618774325874647273446850107217308472599704921431557988258328268814343533200818800770476398207361203467580236085562374300907106235122459320211396640772040483823215609549342831306161721616176933765564983461333503936234487075670534515769060007517458722752470902753316960785929958251393252109465270602716043150298757409338389713831458597564245530307431176833024951122102699669657628081674762760664604054201574988928402849278910852420887832576144904977176190096644087082399109193867856569562616399692377e-4", "2.00855369231876677409285296545817178969879078385541501443789342296988458780919737312044971602530177021536076158519490028818110124793535066902326217844772505039456771000660778518122290478843839402581525347093526229814655384245556977335151081501184047549338384978431776760709137728624917873493960378227937176871312540605975534266408260309486639202162586673911818102143442597648610144743968207859211681419683748223689094322882074779471238071775123479594835563129187071772398032557072240696132676430886975618359507223192869047808426732636908842014209942605405569985196862091776780376305707135430357271227505944350140953683366699214534148855026835050422511531422504252485965345769230848886056984445015631635513964479308659536e+1"),
        ("0", "1", "3.67879441171442321595523770161460867445811131031767834507836801697461495744899803357147274345919643746627325276843995208246975792790129008626653589494098783092194367377338115048638991125145616344987719978684475957939747302549892495453239366207964810514647520612294223089164926566600365074577283705532853738388106804787611956829893454497350739318599216617433003569937208207102277518021584994233781690715667671762336608230376122915623757209470007040509733425677576252528030376886165157093653799542740637071787844541946749093130698056016370211138977422821401738023283246528729138900466098665951244409769985145916428780372020251022457873211105953777680743711220624000516796528097544478028648600683856420043368466248434938692e-1", "2.71828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193200305992181741359662904357290033429526059563073813232862794349076323382988075319525101901157383418793070215408914993488416750924476146066808226480016847741185374234544243710753907774499206955170276183860626133138458300075204493382656029760673711320070932870912744374704723069697720931014169283681902551510865746377211125238978442505695369677078544996996794686445490598793163688923009879312773617821542499922957635148220826989519366803318252886939849646510582093923982948879332036250944311730123819706841614039701983767932068328237646480429531180232878250981945581530175671736133206981125099618188159304169035"),
        ("3", "0.25", "1.56426318841881716102126980461566588450380350341076228557306530116729257751859973695336920080009251501442042818073858963351411413278076229410081576392272592530797134999328076786375554536119890934960685418177807379769464919104853662165338890630044364798688843388065483479816611342383473031587761586030859465037576919019762313699995251441111557313354983707985267081519067636843791542024584131098731904534200460889249799795804514979820890012554894069001354142043884540988787476588766933562932587460066988792793852247670137562296042696161875077695873100425507972928287738747809353780901414593301840741567420675442700630525915427066970222092943736613609647618168326209643336514839656809372794827762772670536198670149990658013e+1", "2.57903399171930620890801076693772218766552688489539736340424977271388720447207274888484147322397605262439439122840412454152164674793952271566093381069299817951987603654957505688725715830662486037892062474725936881525655049663636583678236975809435830044230209822776914902621386703981743612156350722183125642160424023438529351359584877609039507090340549487526420954549894164582212262226980286278281965344204036083363088010405771130138494572390132364021949441836493176972675095033693961721612571235617136295399231153979857230936924976028717661330985218050301255535748701361897719446996544876826504924075042109789425272898286935256628016960570672396746317067279468683518886820150397712996246378365475911024177817040017699199e+1"),
        ("-50", "10", "8.75651076269652033848873280073916603655710748178175890605671542387940202673337676377376373533697467547587557020541937343568555903118409328375744699672238633545462035286934488513505108365550042941298970036685237720741399595074249554661429995537856477801873096676674350133052257141929348965929841685715901044546711277824694187686118436913704980776759067269975668649354653195290332720067162166869929709808167958079806891220176624120264622379923735717676112153658485315405656190081330939565796267056426925906308493599647219902610592771593275457097325237472052865206023754801450679641241403412873499930577123204791561805235517681467937609862505511546386926425517952411391744172826145273897101577471084504886032633040289443064e-27", "4.24835425529158899532923478285865801787956555416644628805081891892603306392691465410438922859472778091036076743217617879470936916039988403905692470827008696314677574665161032219496143471467436667484897089329405982809650254450776062180522708447793080352358557497604028136593819449866314952916829363764401225963004625825819165906249232180700276804224168703857715087159049661988215275588772382460269329268924768404756459421683904186174686586474656880667818198827003408779336434409957087397570846221240584025763086949066540289036296429377616434670180403996303539028122974175803948692853429500649789053056534895698547683464311387911887360834411355944842206090319805854480554286755842906040667924207500942988303714391313169601e-18"),
        ("10", "9.31322574615478515625e-10", "2.2026465774292971692811193938167589148242647044422084954026055604103002291547498718270504488108298685291639083210160234667064070779243560069500139475290392202802068525179776050054037370081704407174561599065737890398939256048260057419964947216591777168810773611601112836588743977994252962355464165872893034157934905600091054243510841470530384873335952739485832326189478291628888135063040490468158387606931146896758886934929067349799675599295316330236194309540245461458631092930979144643794000686983789312400553814191810526623122319829720439037301734493871666928336318630452469518367151980860361853182868109672148633428523779719870889302328793047247281296995230240427845851652449511458536618798381842421977628349982515141e+4", "2.20264658153204613602095210059265757272582318126411805966845368167727862819152644358595231198815991047631492563965707190436744722554869045822910837534201376793304553639798243714759355782768968078627224788903575750577787960821196820115493124117651874380073291251438906602597812963985198018540920618374185959258846430454447087969706192159865639736374452228049423904909460724073867432355932672775035220259878210167533105534540796389132094836671694110656424565276329433397162786983357906695502551181011027072226803200473867937644000463853583555391392338133362744148709377340966969736817961736883890938725778723925347892745779035257871059962325592958396619981769806082643775977720227776650389356423717785733140594339917828893e+4"),
        ("1", "0", "2.71828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193200305992181741359662904357290033429526059563073813232862794349076323382988075319525101901157383418793070215408914993488416750924476146066808226480016847741185374234544243710753907774499206955170276183860626133138458300075204493382656029760673711320070932870912744374704723069697720931014169283681902551510865746377211125238978442505695369677078544996996794686445490598793163688923009879312773617821542499922957635148220826989519366803318252886939849646510582093923982948879332036250944311730123819706841614039701983767932068328237646480429531180232878250981945581530175671736133206981125099618188159304169035", "2.71828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193200305992181741359662904357290033429526059563073813232862794349076323382988075319525101901157383418793070215408914993488416750924476146066808226480016847741185374234544243710753907774499206955170276183860626133138458300075204493382656029760673711320070932870912744374704723069697720931014169283681902551510865746377211125238978442505695369677078544996996794686445490598793163688923009879312773617821542499922957635148220826989519366803318252886939849646510582093923982948879332036250944311730123819706841614039701983767932068328237646480429531180232878250981945581530175671736133206981125099618188159304169035"),
        ("-700", "3", "4.90884390169192177195612911178243419830913088739799883384794203038742453453491117855609963333051245725440211825323762992353408844411571066917879586170649044228244020649696128386393842173012986503695091295623723755383423130309143773301760584980849859045143225033571940200082809326240589537894178343199408448712145712616651646042044453169250422203833678906317455153175073938900351302967555974519922835310572345943902605988213167896676246972949876168846256396695926430057017171882573660664829912875301666974793402900034042932768118367653080548663555492018876853850322292559503446013090361177767752540404168420074374364335969496069960635895595400596576005304965550576992566232211136293323135389308610220377572927781441128691e-306", "1.98036897270374246007628519831576272198292371701921179439676014495165161788645957376301307353141960806149918624403577122135569407601310626014226758757565498726676448273539695546903543523974025884285302173390082913377396624749566200951159605838152469938430456726679672551236920344688215222234796619524103701073041961657815915473890321423781989604739998909195053165022335554099614639384584512533364472663115251085016691829628965333154212356276742825505646738064072497568260334471134063957303723532841825432286582421170167238349370835208806848469281916286066543738373308224185808611327436318417748179163078392786964483896370151174038981208575538430153348371933621390742206510408279244666910455784607581856368604150517817171e-303"),
    ];
    for p in [64usize, 256, 2048] {
        for (m, r, low, high) in refs {
            let mm = ExactNum::parse(m, Radix::Dec, p, RoundingMode::ToEven, &mut cc);
            let rr = ExactNum::parse(r, Radix::Dec, p, RoundingMode::ToEven, &mut cc);
            let b = Ball::new(mm, rr);
            let lo = ExactNum::parse(low, Radix::Dec, 2304, RoundingMode::Down, &mut cc);
            let hi = ExactNum::parse(high, Radix::Dec, 2304, RoundingMode::Up, &mut cc);
            for rm in MODES {
                let w = b.exp(p, rm, &mut cc);
                let lb = ExactNum::sub_ref(w.mid(), w.rad(), 2304, RoundingMode::Down);
                let ub = ExactNum::add_ref(w.mid(), w.rad(), 2304, RoundingMode::Up);
                assert!(
                    matches!(lb.compare(&lo), Some(o) if o <= 0),
                    "lower m={m} r={r} p={p} {rm:?}"
                );
                assert!(
                    matches!(ub.compare(&hi), Some(o) if o >= 0),
                    "upper m={m} r={r} p={p} {rm:?}"
                );
                assert!(
                    !w.rad().is_inf(),
                    "unnecessarily infinite radius m={m} r={r}"
                );
            }
        }
    }
}

/// ZF-46: a disk that meets the cut (−∞, 0] of the principal log has an image on both sides of the jump
/// (imaginary parts near +π and −π), so its log disk must cover both. Points: the disk's lowest and highest
/// points straight above/below the center, and the center. Principal values agree with mpmath/FLINT log.
#[test]
fn zf_46_cball_ln_disk_across_cut_encloses_both_sides() {
    use zenith_float::ComplexBall;
    let mut cc = Consts::new().expect("Consts");
    for p in [64usize, 128, 256] {
        for (re, im, r) in [
            ("-1p0", "0", -60i64),
            ("-1p0", "-0", -8),
            ("-1p0", "0", -8),
            ("-1p0", "1p-10", -8),
            ("-3p0", "-1p-20", -4),
        ] {
            let mid = parse_complex(&format!("{re},{im}")).expect("mid");
            let rad = ExactNum::from_u8(1, p).ldexp(r, p, RoundingMode::None);
            let b = ComplexBall::new(mid.clone(), rad.clone());
            let w = b.ln(p, RoundingMode::ToEven, &mut cc);
            for s in [-1i8, 1] {
                // a point of the disk on the far side of the cut, inside by a quarter of the radius
                let dy = ExactNum::from_i8(s, p)
                    .mul(&rad, p, RoundingMode::None)
                    .ldexp(-1, p, RoundingMode::None);
                let pt = ExactComplex::new(
                    mid.re().clone(),
                    ExactNum::add_ref(mid.im(), &dy, p, RoundingMode::None),
                );
                assert!(b.contains(&pt, p), "point not in source disk");
                let v = pt.ln(p + 64, RoundingMode::ToEven, &mut cc);
                assert!(
                    w.contains(&v, p),
                    "p={p} {re},{im} r=2^{r} s={s}: ln point {} not in {} ± {}",
                    fmt_c(&v),
                    fmt_c(w.mid()),
                    fmt_real(w.rad())
                );
            }
        }
        // A disk away from the cut keeps a small radius.
        let mid = parse_complex("-1p0,1p-4").expect("mid");
        let b = ComplexBall::new(
            mid,
            ExactNum::from_u8(1, p).ldexp(-8, p, RoundingMode::None),
        );
        let w = b.ln(p, RoundingMode::ToEven, &mut cc);
        assert!(
            matches!(w.rad().compare(&ExactNum::from_u8(1, p)), Some(o) if o < 0),
            "{}",
            fmt_real(w.rad())
        );
    }
}

/// Every original ZF-46 source disk, sampled at the cardinal and diagonal boundary points.
#[test]
fn zf_46_cball_ln_all_harper_disks() {
    use zenith_float::ComplexBall;
    let mut cc = Consts::new().expect("Consts");
    for row in include_str!("tier1/zf46.tsv").lines() {
        let t: Vec<_> = row.split('\t').collect();
        let p: usize = t[1].parse().expect("p");
        let mid = parse_complex(t[2]).expect("mid");
        let rad = parse_real(t[3]).expect("rad");
        let b = ComplexBall::new(mid.clone(), rad.clone());
        let w = b.ln(p, RoundingMode::ToEven, &mut cc);
        for (dx, dy) in [(0i8, 1i8), (0, -1), (1, 0), (-1, 0), (1, 1), (-1, -1), (1, -1), (-1, 1)] {
            let offset = |d: i8| {
                ExactNum::from_i8(d, p)
                    .mul(&rad, p, RoundingMode::None)
                    .ldexp(-1, p, RoundingMode::None)
            };
            let pt = ExactComplex::new(
                ExactNum::add_ref(mid.re(), &offset(dx), p, RoundingMode::None),
                ExactNum::add_ref(mid.im(), &offset(dy), p, RoundingMode::None),
            );
            assert!(b.contains(&pt, p), "{} point not in source", t[0]);
            let v = pt.ln(p + 64, RoundingMode::ToEven, &mut cc);
            assert!(w.contains(&v, p), "{} dx={dx} dy={dy}", t[0]);
        }
    }
}

/// ZF-64: digamma of a large negative non-integer. ψ(x) = ψ(1−x) − π·cot(πx); cot has period 1, so the
/// reduction must use the exact fractional part of x, not π·x rounded at a precision below the bits of x.
/// References: FLINT arb_digamma at 9000 bits, cross-checked with mpmath to 700 digits.
#[test]
fn zf_64_digamma_large_negative_reduces_exactly() {
    let mut cc = Consts::new().expect("Consts");
    let refs: [(&str, &str); 4] = [
        ("-9d4370565e678de07a9d9ccdca9a47f476b9f5a11e09273263078f65dcffb4e79a2d3076b707e26064cddd0ad78f592e9ca7fa5a53af47e6dc044010181cc4f8e76bbf994252eda83c921eee0922f530cd33b1242e08197331dc8110f334f743a86e0b1844d42de873468e878501362553790ed5a2025ddda4a5b797ef1d79a2db7e143440078a4317fcb6f17087d1af014aaff0ffc124384270bb522d0d4b2d479c7f023d6e9a279867257b1ca5c760cd68b9e1a4f7c0b53d5f450b25025333366efa767d6e7c84e6bee61d333f27a6d343215c51f2a8d2c3f5ccd5efb5d84a827ae2bd92c4bbe00db351ba441160d480bfc83fe6447a75849a67c740d15d77p-622", "977.797103447905361069154798629579685310265718656030793567766385352886652514901110357702712611656578724138560564501203677483483276887069965702143548550803024432057589331070787547776106796032733927227308549365504464775025700266595102196319693637555506877457063040282476045304427987683901303707143904203656511839285326102341452453150321484714477480814179317984953443123205346988945836335965169624118093303563939959791405012128911958488124718645167742715505377226721379715182790300079757708533432137131010526434963903665924533479735225461562834383887673736069461679989698908004005862049242504680512257284302351715171054708445052786552927398790081797374819657147857895936877993544949841319109514338478513832824478034562404598"),
        ("-1000000000000000000000000000000000000000000000000003p-2", "136.874137819519159263818245151272163867053857741396455511402417519826922664752801902008358908563340333947531352178469120967257554538948954051027888159041744010353472560184268620887359860951060451339345197796828924246015488442113578237410593945867546789189378199766509238353714695712834357824293705099168947979897074218955945364537665105790242847458203457737830253690531649807933050546275730568143604223647526281889320269910543306577056850821134186187933427798118344351356569566677187714053372457367442916477419540074168643415096180783543138029343251049026144342869258976333697786054411422718789346644229077837958031602793158890617113791591164047913066754385885632835435558940449013349441128211159488884000211958239193419"),
        ("-c0000000000000000000000000000000000000000000000000000000000000001p-1", "179.237437692574054211623900451673903700051025088408349760749456773439655139430152316941916641967439018140994451419219285299607444196000910262360879806282326697026919787640207419123720743366703360778908941734226035821369078782330384721557705265661900598352066915068865188164191694554927557574720908028379809620537924874355832459388678663284985050060299502657099130406550345380494830454035857501051454593308560780338563960508728844346756222813698818750732044343471155115335995395798050880592915402320326882130482248889468752638064394759262556643706758656182388446984565559438128846889717783343605642915889216955582093449725247579386701585619680841077024232082602966223739671744217846291700387059502800642552795460313594341"),
        ("-1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000003p-4", "351.275315949282704662428307585503308150471089379757804686930311253399943802214502270594957641968814475143982218154347025282709676419834472117563019982460225542272206021886065600432841521442548785804767007137010533322499548950128581932539614283917939511972427821743057073973463712236529542978166399193524987869862319637973721133823833658662937226724403684476322060646508170445561554374268527684566871071187471648716572547902577722935758149070351571721315087050877359582840505449335144097233704454835171527100656384334516278510389733914547716084417171778244043735182981645587361847294163862003130940356818398891379998609979301248375449459307958085141654187367667617178368160814703225002096426549142123862597140229927021170"),
    ];
    for p in [64usize, 256, 2048] {
        for (x, want) in refs {
            let x = parse_real(x).expect("x");
            for rm in MODES {
                if rm == RoundingMode::None {
                    continue;
                }
                let got = x.digamma(p, rm, &mut cc);
                let w = ExactNum::parse(want, Radix::Dec, p, rm, &mut cc);
                assert_eq!(
                    got.compare(&w),
                    Some(0),
                    "x={} p={p} {rm:?}: {} vs {}",
                    fmt_real(&x),
                    fmt_real(&got),
                    fmt_real(&w)
                );
            }
        }
    }
}

/// ZF-52: default tolerance below the working grid must not make Newton reject a converged point;
/// bisection must be allowed enough halvings for the requested tolerance. References: mpmath roots
/// at 1200 bits; sqrt(2) and ln(3) cross-checked with FLINT arb.
#[test]
fn zf_52_roots_default_tolerance_converge() {
    use zenith_float::{bisect, newton, root_default_tol};
    let mut cc = Consts::new().expect("Consts");
    let refs = [
        "1.41421356237309504880168872420969807856967187537694807317667973799073247846210703885038753432764157273501384623091229702492483605585073721264412149709993583141322266592750559275579995050115278206057147010955997160597027453459686201472851741864088919860955232923048430871432145083976260362799525140798968725339654633180882964062062",
        "0.739085133215160641655312087673873404013411758900757464965680635773284654883547594599376106931766531849801246643987163027714903691308420315780440574620778688524903891539289438845095234801335631276772231580956353776572451204373419936433512538409780034340646700479402143478080271801883771136138204206631633503727799169673122323006139",
        "1.32471795724474602596090885447809734073440405690173336453401505030282785124554759405469934798178728032991092099474220742510890263904589779559431475709672347175416683903886741875173693158425354990824662235453372735045898799095681506277455098024862130121698941575245745486250756265246103689389048399322699520749759628288685569081507",
        "1.09861228866810969139524523692252570464749055782274945173469433363749429321860896687361575481373208878797002906595786574236800422593051982105280187076727741060316276918338136717937369884436095990374257031679591152114559191775067134705494016677558022220317025294689756069010652150564286813803631737329857778236699165479213181814902",
    ];
    for p in [64usize, 128, 256, 1024] {
        let rm = RoundingMode::ToEven;
        let tol = root_default_tol(p, rm);
        for (id, txt) in refs.iter().enumerate() {
            let f = |x: &ExactNum, q: usize, r: RoundingMode, c: &mut Consts| match id {
                0 => x.sqr(q, r).sub(&ExactNum::from_u8(2, q), q, r),
                1 => x.cos(q, r, c).sub(x, q, r),
                2 => x
                    .powi(3, q, r)
                    .sub(x, q, r)
                    .sub(&ExactNum::from_u8(1, q), q, r),
                _ => x.exp(q, r, c).sub(&ExactNum::from_u8(3, q), q, r),
            };
            let df = |x: &ExactNum, q: usize, r: RoundingMode, c: &mut Consts| match id {
                0 => ExactNum::from_u8(2, q).mul(x, q, r),
                1 => -x.sin(q, r, c).add(&ExactNum::from_u8(1, q), q, r),
                2 => ExactNum::from_u8(3, q).mul(&x.sqr(q, r), q, r).sub(
                    &ExactNum::from_u8(1, q),
                    q,
                    r,
                ),
                _ => x.exp(q, r, c),
            };
            let x0 = ExactNum::from_u8(1, p);
            let want = ExactNum::parse(txt, Radix::Dec, 1152, rm, &mut cc);
            for got in [
                newton(f, df, &x0, &tol, 200, p, rm, &mut cc),
                bisect(
                    f,
                    &ExactNum::new(p),
                    &ExactNum::from_u8(2, p),
                    &tol,
                    p,
                    rm,
                    &mut cc,
                ),
            ] {
                let got = got.unwrap_or_else(|| panic!("root {id} p={p} refused"));
                let diff = ExactNum::sub_ref(&got, &want, 1152, rm).abs();
                let allowed = ExactNum::from_u8(1, p).ldexp(-((p.min(256) - 4) as i64), p, rm);
                assert!(
                    matches!(diff.compare(&allowed),Some(o) if o < 0),
                    "root {id} p={p}: {}",
                    fmt_real(&got)
                );
            }
        }
    }
}

/// `|got - want| < 2^(e(want) - p + slack)`; `want` parsed at 1152 bits.
fn quad_near(got: &ExactNum, want: &str, p: usize, slack: i64, cc: &mut Consts) -> bool {
    let w = ExactNum::parse(want, Radix::Dec, 1152, RoundingMode::ToEven, cc);
    let d = ExactNum::sub_ref(got, &w, 1152, RoundingMode::ToEven).abs();
    let lim = w.abs().ldexp(slack - p as i64, 1152, RoundingMode::ToEven);
    matches!(d.compare(&lim), Some(o) if o < 0)
}

const SQRT_PI_HALF: &str = "0.886226925452758013649083741670572591398774728061193564106903894926455642295516090687475328369272332708113411812141285333118076432862211301262546854801393534231018849326552561424962586514475413114466047689633981400087319507675739860258350095092617009292723487247456320156960887762953108202709666250453199203806866738737576716833995";
const SQRT_PI_3_4: &str = "1.32934038817913702047362561250585888709816209209179034616035584238968346344327413603121299255390849906217011771821192799967711464929331695189382028220209030134652827398982884213744387977171311967169907153445097210013097926151360979038752514263892551393908523087118448023544133164442966230406449937567979880571030010810636507525099";

/// ZF-49: odd-n Gauss–Hermite counted the origin twice (H_n(0)=0 opened a bracket at
/// the left end of the positive search), dropping the outermost pair. The rule is
/// exact for degree <= 2n-1: x^2 -> sqrt(pi)/2 (n >= 2), x^4 -> 3 sqrt(pi)/4 (n >= 3).
/// References: FLINT arb at 1400 bits.
#[test]
fn zf_49_gauss_hermite_odd_n_exact_on_polynomials() {
    use zenith_float::gauss_hermite;
    let mut cc = Consts::new().expect("Consts");
    for p in [64usize, 128, 256, 1024] {
        for n in 2usize..=13 {
            let rm = RoundingMode::ToEven;
            let x2 = gauss_hermite(|x, q, r, _| x.sqr(q, r), n, p, rm, &mut cc)
                .unwrap_or_else(|| panic!("n={n} p={p} refused"));
            assert!(
                quad_near(&x2, SQRT_PI_HALF, p, 4, &mut cc),
                "x^2 n={n} p={p}: {}",
                fmt_real(&x2)
            );
            if n >= 3 {
                let x4 = gauss_hermite(|x, q, r, _| x.sqr(q, r).sqr(q, r), n, p, rm, &mut cc)
                    .unwrap_or_else(|| panic!("n={n} p={p} refused"));
                assert!(
                    quad_near(&x4, SQRT_PI_3_4, p, 4, &mut cc),
                    "x^4 n={n} p={p}: {}",
                    fmt_real(&x4)
                );
            }
        }
    }
}

/// ZF-48: Gauss–Laguerre refused n >= 15. The uniform isolation grid on [0, 4n+4]
/// was far coarser than the clustered small zeros of L_n (x_1 ~ 1.4/(4n+2)), so
/// sign changes were missed and the root count fell short. The rule is exact for
/// degree <= 2n-1: 1 -> 1, x^2 -> 2, x^5 -> 120 (integer values, exact).
#[test]
fn zf_48_gauss_laguerre_large_n_exact_on_polynomials() {
    use zenith_float::gauss_laguerre;
    let mut cc = Consts::new().expect("Consts");
    let cases: &[(usize, &[usize])] =
        &[(64, &[5, 14, 15, 16, 20, 30, 32, 48, 64]), (256, &[15, 30, 64]), (1024, &[30])];
    for &(p, ns) in cases {
        for &n in ns {
            let rm = RoundingMode::ToEven;
            for (k, want) in [(0u32, "1"), (2, "2"), (5, "120")] {
                let g = gauss_laguerre(|x, q, r, _| x.powi(k as usize, q, r), n, p, rm, &mut cc)
                    .unwrap_or_else(|| panic!("n={n} p={p} k={k} refused"));
                assert!(
                    quad_near(&g, want, p, 6, &mut cc),
                    "x^{k} n={n} p={p}: {}",
                    fmt_real(&g)
                );
            }
        }
    }
}

/// Harper's ten integrands (h2 aux.rs `integrand`), shared by the quadrature fixtures.
fn harper_integrand(
    id: u32,
    x: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> ExactNum {
    let one = ExactNum::from_u8(1, p);
    match id {
        0 => x.exp(p, rm, cc),
        1 => x.sin(p, rm, cc),
        2 => x.sqr(p, rm).add(&one, p, rm).reciprocal(p, rm),
        3 => x.sqrt(p, rm),
        4 => x.ln(p, rm, cc),
        5 => x.sqr(p, rm),
        6 => (-x.sqr(p, rm)).exp(p, rm, cc),
        8 => ExactNum::from_u8(50, p).mul(x, p, rm).cos(p, rm, cc),
        9 => x.cos(p, rm, cc),
        _ => one,
    }
}

/// ZF-50: tanh–sinh stopped every sum at the hard TANH_SINH_K_MAX = 512 abscissae.
/// With h = 2π/(p ln 2), the tails needed for p >= 256 run past that, so the result
/// was a truncated sum, scaled down by roughly 512/k_needed. It was also returned
/// unconverged after the last level. References: FLINT arb at 1400 bits.
#[test]
fn zf_50_tanh_sinh_high_precision() {
    use zenith_float::tanh_sinh;
    let mut cc = Consts::new().expect("Consts");
    let refs: &[(u32, &str, &str, &str)] = &[
        (0, "0", "1", "1.71828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193200305992181741359662904357290033429526059563073813232862794349076323382988075319525101901157383418793070215408914993488416750924476146066808226480016847741185374234544243710753907774499206955170276183860626133138458300"),
        (1, "0", "1", "0.459697694131860282599063392557023396267689579382077772329902744618899605225528235482048143912816910656428268839969910902139366239978336543593487734582681415282028835525520505766882075448606745664056482243297107403624263845672450358245508224884868777269899368642921767763228598482531004063321269325772379754922362559324125018382157"),
        (2, "-1", "1", "1.57079632679489661923132169163975144209858469968755291048747229615390820314310449931401741267105853399107404325664115332354692230477529111586267970406424055872514205135096926055277982231147447746519098221440548783296672306423782411689339158263560095457282428346173017430522716332410669680363012457063686229350330315779408744076046"),
        (2, "0", "1", "0.785398163397448309615660845819875721049292349843776455243736148076954101571552249657008706335529266995537021628320576661773461152387645557931339852032120279362571025675484630276389911155737238732595491107202743916483361532118912058446695791317800477286412141730865087152613581662053348401815062285318431146751651578897043720380230"),
        (3, "0", "1", "0.666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666667"),
        (4, "1", "2", "0.386294361119890618834464242916353136151000268720510508241360018986787243939389431211726653992837375084002962041141371467371040471516261114065341503270151923861455141656742870380614077247783346942246700230728995910478240950345363149864130311049468279051765900906014190652733285308208478315629904087480860771001603888341283343037289"),
        (5, "-2", "3", "11.6666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666667"),
        (5, "0", "1", "0.333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333"),
        (6, "-3", "3", "1.77241469651904246778896915582369115913928386949056681168932665253569068978410511638758646210551290988862782979587731223697725520635433992700384918158463645183260796910421404838098631032866808383917414280986457976121565354387638714549482138025787014843094526931583544382012711208533257681772610460568994468579911262686447639885280"),
        (6, "0", "1", "0.746824132812427025399467436131853005354499686812606329027654498958605327561772831497848429822901919730618954017089471349498097131491048662627352146982179779594358148099229331600834712708984912805053638950832959892321512680543298847843764313959957959978504210515364118213547939212780905389905695880666452198774298219871339902211690"),
        (8, "0", "1", "-0.00524749707407857571828787293825245091773331929942477045887871184828917905179901503577847382570704739317198027405306417222213647304310017250348792583103888631908764942515594820929406369254832813980963530228047177295338223389756335547195140623093431833416637568940168458545776033945301316421969951103154288667798249064536349136985842"),
        (9, "0", "1", "0.841470984807896506652502321630298999622563060798371065672751709991910404391239668948639743543052695854349037907920674293259118920991898881193410327729212409480791955826766606999907764011978408782732566347484802870298656157017962455394893572924670127086486281053382030561377218203868449667761674266239013382753397956764255565477964"),
        (10, "0", "1", "1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"),
    ];
    for p in [64usize, 256, 512] {
        for &(id, a, b, want) in refs {
            let rm = RoundingMode::ToEven;
            let a = ExactNum::parse(a, Radix::Dec, p, rm, &mut cc);
            let b = ExactNum::parse(b, Radix::Dec, p, rm, &mut cc);
            let g = tanh_sinh(
                |x, q, r, c| harper_integrand(id, x, q, r, c),
                &a,
                &b,
                p,
                rm,
                &mut cc,
            )
            .unwrap_or_else(|| panic!("f{id} p={p} refused"));
            assert!(
                quad_near(&g, want, p, 12, &mut cc),
                "f{id} p={p}: {}",
                fmt_real(&g)
            );
        }
    }
    let one = ExactNum::from_u8(1, 1024);
    let g = tanh_sinh(
        |x, q, r, c| harper_integrand(0, x, q, r, c),
        &ExactNum::new(1024),
        &one,
        1024,
        RoundingMode::ToEven,
        &mut cc,
    )
    .expect("f0 p=1024 refused");
    assert!(
        quad_near(&g, refs[0].3, 1024, 12, &mut cc),
        "f0 p=1024: {}",
        fmt_real(&g)
    );
}

/// ZF-21: C99 Annex G: Re catan(±0 + iy) = copysign(π/2, Re z) for both cuts
/// |y| > 1, and Im = atanh(1/y). Harper's remaining 60 rows (y = -2, -3/2 with
/// Re = +0) expect -π/2: mpmath's raw-real side, not C99. References: FLINT arb.
#[test]
fn zf_21_c_atan_imaginary_cuts_follow_c99() {
    let mut cc = Consts::new().expect("Consts");
    let half_pi = "1.57079632679489661923132169163975144209858469968755291048747229615390820314310449931401741267105853399107404325664115332354692230477529111586267970406424055872514205135096926055277982231147447746519098221440548783296672306423782411689339158263560095457282428346173017430522716332410669680363012457063686229350330315779408744076046";
    let inv = [("2", "0.549306144334054845697622618461262852323745278911374725867347166818747146609304483436807877406866044393985014532978932871184002112965259910526400935383638705301581384591690683589686849422180479951871285158397955760572795958875335673527470083387790111101585126473448780345053260752821434069018158686649288891183495827396065909074510"), ("1.5", "0.804718956217050187300379666613093819762800677134258860956323945737089493853828882315066939046589805399983151510857781449862002614662338099816808308731852863776089818748591622826746428101170762528635077596800439888694862844096770356383077365610904763974260646410679029861283836142643620230794724089182335664336999212318879796594712"), ("5", "0.202732554054082190989006557732174568285995211731247098807007162072050335624457125633876213908656700622984274022693590004341241995086194632010065559566100724336217599177504965994533330110234643216309610100715706782953235721288994886061437505625423041338643651246945232681869939425769314490110682664779073698433094855310745073315188")];
    for p in [64usize, 128, 1024] {
        for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
            for (y, at) in inv {
                for (re_neg, y_neg) in [(false, false), (false, true), (true, false), (true, true)]
                {
                    let mut yv = ExactNum::parse(y, Radix::Dec, p, rm, &mut cc);
                    if y_neg {
                        yv = -yv;
                    }
                    let z = ExactComplex::new(signed_zero(re_neg), yv);
                    let w = z.atan(p, rm, &mut cc);
                    let tag = format!(
                        "atan({}0 {}{y}i) p={p} {rm:?}",
                        if re_neg { "-" } else { "+" },
                        if y_neg { "-" } else { "+" }
                    );
                    assert_eq!(
                        w.re().is_negative(),
                        re_neg,
                        "{tag}: Re sign {}",
                        fmt_real(w.re())
                    );
                    assert!(
                        quad_near(&w.re().clone().abs(), half_pi, p, 2, &mut cc),
                        "{tag}: |Re| {}",
                        fmt_real(w.re())
                    );
                    assert_eq!(
                        w.im().is_negative(),
                        y_neg,
                        "{tag}: Im sign {}",
                        fmt_real(w.im())
                    );
                    assert!(
                        quad_near(&w.im().clone().abs(), at, p, 2, &mut cc),
                        "{tag}: |Im| {}",
                        fmt_real(w.im())
                    );
                }
            }
        }
    }
}

/// Correctly rounded `want` (a decimal reference good to its printed digits, parsed at 8960 bits) at `p` bits in `rm`.
fn cr_ref(want: &str, p: usize, rm: RoundingMode, cc: &mut Consts) -> ExactNum {
    let mut w = ExactNum::parse(want, Radix::Dec, 8960, RoundingMode::ToEven, cc);
    let _ = w.set_precision(p, rm);
    w
}

/// ZF-20: c.expm1 was exp(z) - 1 at one fixed working width. For small z the
/// subtraction cancels about log2(1/|z|) bits of the real part, and for z = iy the
/// real part -2 sin²(y/2) disappears below the grid, so directed modes came back
/// 1 ulp off. Now: Ziv on (expm1(x) cos y - 2 sin²(y/2)) + i e^x sin y.
/// References: FLINT acb expm1 at 1400 bits, rounded once per mode.
#[test]
fn zf_20_c_expm1_small_arguments_correctly_rounded() {
    let mut cc = Consts::new().expect("Consts");
    let refs: &[(&str, &str, &str, &str)] = &[
        ("0", "1p-60", "-3.7615819226313200254999569191111861687839048056595677675907888118907323193057665362845823745350500340693842487841299284149949698316892228266996644532460714933674816471634135245010951272787307441215099368126394208979248267302881966825731332361082285399936145427703179256939639231656745515858075940852395063066841216285894e-37", "8.6736173798840354720596224069595336903186992553335792455482843151239087407573817358659866299826920237132650895911463153195694128700658996152528572403483008172903295898706172258532732006558207095470277585160849183364043353298282444794502978058085541516495296391384851280147005982487503128541727515721075427069291537514599e-19"),
        ("1p-60", "0", "8.6736173798840354758212043295908537179937576637776069237072704425578836293896585968467802871185525450575019930342518361422796890094594259626224871759093042035782510950622932259303154693349842950025823472325803002756275505605974387484565953260009967809136071440485787697276016488556274989928500084888115984656552582721853e-19", "0"),
        ("1p-60", "1p-30", "4.3368086899420177363432763636990435142138374675537734857259846480165036270525105799921820635383966806135340791716474314325236716583577822178381052089597491914657577690158727297339946252353058005789483660739349593343167687945437935820781139099669791392711520318772866034003556618056545328233418892343445826627545787131619e-19", "9.3132257461547851629816130578859674085739686336346716565952671390880972822246406916511256507188306349577849165554087875086581270779652070527259388484581890682583665994629381113807580549402889650964820759751260824123643781653128667307752954613035361027622992178367605024671922455773657109682799635337643965511797296557384e-10"),
        ("1p-41", "1p-20", "-6.8931717712769024381162421334863835905379877893253540134834699680380971000866761689796226623506312221196131852728873959536013560846069847863880949051317948495095115504863546022126209803318154894421685743423882161152798180893113086991487310365896641000157701887087185162131461557416248129979416296115662256499967427616819e-26", "9.5367431640653912057932950729211391513366933763699227814091851005468559082430865268502400209103917721815738057627897780759531816433698392308985675605504827059454423532023780467833505655161912978054582095033604357720388382348538611436980485259867124375513077894719271184143369133353487155740679235446159952428694421476338e-7"),
        ("-1p-40", "3p-2", "-0.26831113112684458083774791130233651084079640935874388786267461872457884352801180395260542187292002999304989218228155461657415123685002568740535610799208811519921927580717001549017253910351905294140479834699167674654627551530015373994108866933406469240395420655794223828841763772156542680809763312689711620658023055032151", "0.68163876002271421989247794379240190136191422315472499885734812817076599385111159343973468420530921652514650400148975362080628500434098886617613340745948905071666208874199962025853902977931255992124021775813390623606041508695109197065598443491422254772452050305993614210598890699434411769418568429876547079441672307177345"),
        ("0", "-1p-100", "-3.1115076389305708535720320268900621202951260843605835665505581460368529649627636865073626912981964930526836260779418457389054794963876028000538817632884502756923941155947703768580612761462312969919646226633598371580865904460833376517747173163313812234488503646704572772380229701095314565510622366083615167391426640166068e-61", "-7.8886090522101180541172856528278622967320643510902300477027884884583807837122411507046741689540595041307479167584058147418392599521533764197076898119627739255061988255761752002193612542137083185029762058325304258948677852605808340784904354102492225824714249798896155834447566852292791342645301742720107313691862869435072e-31"),
    ];
    for p in [64usize, 128, 256] {
        for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
            for &(x, y, wr, wi) in refs {
                let z = ExactComplex::new(parse_real(x).expect("x"), parse_real(y).expect("y"));
                let w = z.expm1(p, rm, &mut cc);
                let (er, ei) = (cr_ref(wr, p, rm, &mut cc), cr_ref(wi, p, rm, &mut cc));
                assert_eq!(
                    fmt_real(w.re()),
                    fmt_real(&er),
                    "Re expm1({x}+{y}i) p={p} {rm:?}"
                );
                if ei.is_zero() {
                    assert!(
                        w.im().is_zero(),
                        "Im expm1({x}+{y}i) p={p} {rm:?}: {}",
                        fmt_real(w.im())
                    );
                } else {
                    assert_eq!(
                        fmt_real(w.im()),
                        fmt_real(&ei),
                        "Im expm1({x}+{y}i) p={p} {rm:?}"
                    );
                }
            }
        }
    }
}

/// ZF-20 at Harper's high widths (p up to 8192). References: FLINT acb at 9400 bits.
#[test]
fn zf_20_c_expm1_high_precision_correctly_rounded() {
    let mut cc = Consts::new().expect("Consts");
    let refs: &[(&str, &str, &str, &str)] = &[
        ("0", "1p-60", "-3.76158192263132002549995691911118616878390480565956776759078881189073231930576653628458237453505003406938424878412992841499496983168922282669966445324607149336748164716341352450109512727873074412150993681263942089792482673028819668257313323610822853999361454277031792569396392316567455158580759408523950630668412162858939904149389178160478770340861234071999792368620854807170122872927274448555922077988688650009621668936391151525363358743586949110625125194217919477574938197223597301785695225382906993750682391365201383144000258489939050835987508742996133796298320690178714756845626703252921812623973568329845503413966189762017958874648952598933309091740915001155455042952298286281124191226885200698389561840893883171868761331015715276976220780688792908048754507506481554866990501538206388291949085382055193635596188800517894576100385334563273356099360533194670410580773381939723753172521762416399587188283668079290380248306771871823996911027286719412799845272287679366418559397316406622053519609276812910906563778275984804373616750955432102479980010669785019576785182833672478224982386600632574025104999445694303861744745702142983987641245566833503555502599596286820100704727450326302737375641095668255968387966287690428163270308755084814313911242519189036758790157134702855960461613248354392924513961410665805128724451657434740360188743256137679626878296648452713111973059896662181346360057058364111835597878073706096087035399126043769294521508842352340554697811881360274288181513392797038452078330652771325967027764963426938484758326232745601070252175675470712655974469858620320044438068513396822843479750139189491166643577664128624753474975729907053837995178370150782000062700945062910143242689686781290232757052296843307104572909213164349280940747906729154331769263211527949595581727609597723109697771551876241859889885531923432007575494394861837421664960177652586674170626707012025722509140415207983636962901669664158299150556596280570006904036452084620548573319892508649098232232020801095844961375282398415660493157590831803780521260669856391883471234930288951301183771916480058026777733367958035978877907532932947562582587580605830935197196847661123437014843005399861325996236449110915965380753089642592906038541409626559528387695722902814591582395989587468479386231280801976955839220472497675359341697474567430238791446924466989286631430391748632507595710768749085580310516105835072286829629659023644276361530970849820525402499533393443992015293794331046613286015493586872493296729399873110868317114507421421121884441890572246715248791306249087435698790887523735591125063903356357469392755419879755345090661379728559639889992326608319537601296570097371385829546503709294969127635621828200142e-37", "8.67361737988403547205962240695953369031869925533357924554828431512390874075738173586598662998269202371326508959114631531956941287006589961525285724034830081729032958987061722585327320065582070954702775851608491833640433532982824447945029780580855415164952963913848512801470059824875031285417275157210754270692915375145994801345464892768807697491910558621368882660646036132103549194997399302120945424924603158994639053330811817164849443198521051899730068356200676463582368509172211520474599619121903712681080290267620292159185707118627315369347627420251253369381905374492102521083063089242194910598455246650798023168506010127407520476506891872370026916264701583024170576533537945079472859608595756524815710701180191864745993598711186311786259341799912446129104043242948644088887969263583042931790368711286638156295469037954335424924293985638249217537751710932282163461379043969359758882756121207382649235869118023905262611417586957211861381421315564231848813255012289142799843850297221891936243542882654455226749592757349793291577632005514632323022251454823016193756859821777116786234035743367958574356885884014377485689686344279139290955817649206615095640942271910999374685554691071409584757919224491972000392960903172570274284920882195890809692748843074870990529325821069934551888185221503741499406465614260202738940165937721454023784222349985198578174097003902578016841336835628954768296137984782762676857742118568180941412484506651926139338784341481997682655821721805316190756075599201256855400955626835437784763484372866687751259678396135577414668109008656838684916123660642853311349314907753179811082123522666320558442050976703452244627038055465387826217166708813265928393802754854566406029859152189332686489553055050370705899807859039256928995459979909592518744888231502850603759325693960868256172965698439208793702217074103022876060171622214987673145040550160584083206742043932410506847462135345276545632290858543478431130294775567932453146715636418931124948648201443680770529519786639234191234979605447406424122433240804800346684653299907938564987429241468895533560315848779259834824679658459167645109130061720505230768292298087426324703530564561965291002997034399006439339909819769321616823889271343889186525748381584408435046392673763897920424389863418564252618140374863836917141488272042135915058437836265906071586307079462729834573600225476171503267846090823935030510329858781619297622362795359609980625588560255800263898741149540570226857155806756528145820885317052664416288762056762459735621680234744852096567538979854857175421428501741902026246550105934680699816797842957048743672093409202998518519515599606994367903611933184493666304182093277531550968838870886247801481112358049284923e-19"),
        ("1p-60", "1p-60", "8.67361737988403547205962240695953368923114851066715849015326872620336827632044220566627528848949844987105616895933723145365141710126043537835665103639778881515906056116617003815487417901557288955989830149690809741837949174003888136704033229443705975250115737064878903842728242231699103786170479955027212702784377602036579350791088311947382738124503418430704740797107470745427993643922751185875034291335730169381251541281155740556438699176208553143382290595026909350557186809951077689763830885662348890905821276844407339563904567529804301142985976578045272433050857643405276640389255354092588254366896481718117161037786041136276885681410124913912887366466613418014139766093246977254002058210219714891282216381740682715636152268512109061964823388546185481834854770226638716222362323262592243205680352887914059534458649149484979055220306886117252359203733380857009151096966577763598026847926525519739338647116840852477444723615578677757867560628175778834696462226886394810184562351861959640033239929383641799029366234685084078256499313935093105052200288626893858042949788215243424081641945987460446191414647223224181808630884323550282040934430800045415312901687124752233353368525762200102231905850340936711040391494924460630252986034784139067275881149319845984990657287556977619061467892966630043928626926265268184077131549723252510542346026537085349514153907181774821211426282714678866611545342312327125249420766771100423285543617193519720774463241016150441278628150253144720749867055561932350590251213553541406183634624356006270235822487646087569333250940044635796766515904576903727703784100578077912403557395336189599749511747508450569018993654162125505546157050735962316464338996351345818505621802674768269852242981743343092292340462042510133238711050339479731645944705353978905832029421220345563536278173253715969962665916355153807207999364301590280931394552195075628382323719016403116225136933601566510682757982574238521290904825649652806535382707404914281663948262887427595097948900388361726627850046438158623521994264003726132814332366068282332838835731865261308746030682065544076460751207369319612049250669383797371402605762334230022047011828290409340034480660632773609602243880111684691607108304395120677800163030610315977259752765928808937118263056707367512381868027934565778027816433540875646240529093168495049038724489166910789921526853896530453157956262052181818209772223890313003580579911843393380152943279911393402535910948684396058567817626976085438183838722449831170391370436158936917641415279286044144210987154215907517047269542814590008733443932238585696318462001089017784073520406926389517017489629793996175825271749558839057528060613575039119569688504003895469206860992473405424741e-19", "8.67361737988403547958278625222217374458126532755521384694289093309354265605942922402826842461262966135138577431734764205313379289196834640141842345598459283822507929335143135834549962090586657087866638662074073898589046158519191238021732682592572583073366401515723125475276237236940243303449160685176416614145086997554151248637725917748781488572731129282940913251452318746298471256503577116188428393513486708252419316617107182452201082299563117551969864480950400013652783686268177847002986229350580064303703671220312997406581453820332615077218334792913241330330196506675790800169640575447017783193318564858730121515752809683812074160489633546402125556142198730621487231325138596104857450951807480116161894703007242367254431976487875866436029722767182464066038802359176923519666835563760248835334383288982123646279473319835100078089647854753872397417292704460608671929149546948153410956776476983296357771673467658776881139228476039789745747115926281548328099764253259164327951660623977790944867428735960130022725614714625952585101671265799351925614066657465912900861020765021788497540208153035200136521578675748788625618122089014961697182768357339539647256696627659992838565005206972123173727006953727518916605779928295582214501156000818917478646819904780478291222684944835148129713396731211091579196663068735273293898038752257315142675927476752141084397555970387501959338127627787783954940211491971428035815254696155417529341350267403556585327804774597423450983914292033774969941005609255708154089386737027875848080002011783912635261548633266295954765116459955994749141105686829972996329771453508792037571180890994630817401085615477565524741182031638365115655068270210053520915116463914422086521683732532420793542253742538088960295206447840952614351514994441041872732325243299057162663854431794685612623891768240741656887394676823135024662227240951411645413197201759236032030777555787665958800490107300312975854641440922541251357479282375307159717886487249188942081862271560972858712521504245518673882751494714068256687712611208811018491137186173490065633387432874860514753715432903427676098151223971605020705200489122415451661580179721829771153703222540385486595556680095470772800517222518918825092233226453028976207563534189657633778552306348256562482880128710915814296960148761880271384496326813498097474659938478394627412967098321518778007379953556006260412092780129193872809545483717890895258701074793409319684626764325189833145309410477009729945329131774912840888392083447730665004047160512805828105060353545333464495319635712607460904573978127156037490632838346709410146342851047658241057782373514860640991350649880801578664393644247886797961807433231794314582663075075625860246972061863857237e-19"),
    ];
    let cr = |s: &str, p: usize, rm: RoundingMode, cc: &mut Consts| {
        let mut w = ExactNum::parse(s, Radix::Dec, 8960, RoundingMode::ToEven, cc);
        let _ = w.set_precision(p, rm);
        w
    };
    for p in [1024usize, 2048, 4096, 8192] {
        for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
            for &(x, y, wr, wi) in refs {
                let z = ExactComplex::new(parse_real(x).expect("x"), parse_real(y).expect("y"));
                let w = z.expm1(p, rm, &mut cc);
                assert_eq!(
                    fmt_real(w.re()),
                    fmt_real(&cr(wr, p, rm, &mut cc)),
                    "Re expm1({x}+{y}i) p={p} {rm:?}"
                );
                assert_eq!(
                    fmt_real(w.im()),
                    fmt_real(&cr(wi, p, rm, &mut cc)),
                    "Im expm1({x}+{y}i) p={p} {rm:?}"
                );
            }
        }
    }
}

/// ZF-69: c.tan / c.tanh as sin/cos and sinh/cosh quotients at a fixed working width.
/// The complex division cancels the small part entirely: Im tanh(0x100001 + i/2)
/// ~ 2 sin(1) e^{-2097154} came back +0. Now tanh(x+iy) =
/// (sinh 2x + i sin 2y) / (cosh 2x + cos 2y), with an e^{-2|x|} form for large |x|,
/// tan(z) = -i tanh(iz), and Ziv rounding. References: FLINT acb at 1400 bits.
#[test]
fn zf_69_c_tan_tanh_small_parts_not_lost() {
    let mut cc = Consts::new().expect("Consts");
    let refs: &[(&str, &str, &str, &str, &str)] = &[
        ("tanh", "100001p0", "1p-1", "1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000", "6.54889170640323520743151517912337410109558326556177870459202766641107156296447075869000302842985255259865767050336185103701198712495050421119154197673439214999752754619284766341186505680292800413851303245682901185927552135769749572071227899588425835088428260638123610203080888855968057886255929894199961615584252597481898608611621274022552343605790271266206556126509457486224105801246260903967069799973006421458209621931906736755546967346851829451683928780240428881736191928859040480545602417455432846200425189176442247791802841725597852899627171104014696993781326197065054710754466864561956029696702735555035813003351874643267877751128867364488980130586781480192483575980420713640481503491337006003597061123997584350269539194638503821934932941593114095808020822562339727376668131781021586984520047777699224366794861542424848707846976636861781803308679687607435095424645749774197409548547178440486470936809856849177677679036067885670303144137680496300963501726177261606610733703111265697534036605310963490910632315331658967956581613893335497923224967778473368634408394641948959883435424884354819932702275316517649799793479064338055924179907491714529328499906666353149323583481711491843355404282905994986711361021011969546249840760776826861815687551099988814886347489913783680260383777616525586461550410883777309034238706978932758227117613758818337000321836112675083945544061140164916858712653411321695478893967491212823748061995122403943064963539503595572911508624204883385222018464822987211554939597205906697912974482260132319233441590882301384855409470553678624844820797934959570787635573038407644743219901695962866085499793412138044635220979692243905489834566707454499945189557158799191252648837819750524602750104640562855700475331157429414947003154080178688745019666611466542653212100382803895989443787625594238036513942827287428958865850741809588418389847046807041443657428223137197364911044267663174864969009738747777258011692933146691690633279445940270584627765508832006162835561836114967033138059177508828320614918066249000248020209826778940873003267299973614872164386732182864102014032027280521919359564048997397945144769657135142369210983512953763236479900719986821190068031587237229478483322163805134786713663911503691797964223636954793712844025911604776072410134240922857628904169208023780692063174048167651319833307768542422205209165903796401336046787637962019647821529983630032828646613431354953570116069669996276237289698353603592063132955226020584363233237142229528830878403705754816640634181171972506947747882718845899667982920838166819691809814890896386157931239012258990263038479109929863149160365509423770090096443858669893999339241412686759379122490116209587983104224256881808746e-910783"),
        ("tanh", "-b3db40f09dff5b13p-58", "aad5eef4569d62adp-61", "-1.00000000000000000000000000000000000000055110314768364587742003752274440797163013504488047393662461342755857635252887679664679377399448657949898991297482820966338393408384297483861120327123644547325718902343056291642948774496279978617826740912066521931322257658436539903792440604926510942536023065215635717635395691011031208911384315947462727802433371413539889224843168248842137570971203172639405316711350734450346666473364217409461735111785445603647665788747687872719206316399935252376257731043054705884026486643657415235903999099590673460066423994224459054417608653249925201844196085414018425241779358472903406689211716134937061155829107554789846197371253759151925712444629791475836181954364494232634610068129689414171280999309664564919006704243444030171280384123240277047758665026002508037509712656734223205144517127713199002873787757049830066661615229510844083579098030047923104651563926556666447267224775259016634848357672408672337962275389868555762473091741174082962916117318530454794116336037807856349488445941887820819376933660929306873167762020789688903700232024838135416479714771453128462125907811931996027894938100181569338242310413318118121297564706470210996903569244877559046835190437553870329777774985653022662604393035827230548337059313788770668093751957493691768364833917990633494297782725307702608185713938363991710804873685513800338734052386724725013290728548473469415035319859339226472181712182427089238218777287158100870974569650495144977815616559090402664143788137181934198618074394005313593051582086400774222033587676817485213041910547892968846168990078750858680581762646720030714862133831328263480767405840771836884061287544988735384010371275611952772576571414297163639470081949840999437234748850014171263930435341363600576180214640911572692923217603590708683804329815243493924121518985212474558825301106724738456270525342172418224955037585190970123408144002770684386994717301881637702315388540144335454053287235141604307055565410424516400096380366628010868169457256670622383302921951123095241453250123347953763561462682421327952908060038092384418689802841406845021870691772986136403112251821641841523829395954770330600639811454467907245769303075411766998755180991325465072456981005803222378206140153685958221495923808275407964514195013801688176914006762785162913420692347964647429055319004718436572224727034141006672715527038853956975595562435478631486205471624413917232698457581386593779877352569633051734355594704038395094173509545098042514513555029383958119903393867605096121850453767635504391336146001410481301940909037638894443295358862478372721779761737178705233944506178987546579221862588976884264687631535994105792412650044690065235033031876180887746683520483340905641", "-1.67227445694201778763379600579414763374446573471962639400201226984821891037227604292150038917725374925177402581207081637891839692192980156657085545673890704393253858176888972246138399979079023927475158567504497760514825492418497974193082859328870647041441305316235870931121413375451818506691320003709706641157694548288903404330140799403725572665949127462208120923146822259315705533815992508617069246135320724979236538815238837391029631590861359400270492627447952765615081370522815012783285759770474379815718886829678996402039217853759360783687456498714040870284547489811734199237212960456149959810514298274425897538849244777675664276837555297905351415863367055348481949611657510111893412102238929631391838367762421485283482402334338593568189965832713011048291285845689562608289779799868591778009410295740202894809690905081909741133632093322795342730843683736896285970391690571636943315230245527817013061222616085289380751843762484341828938678428557387873994880823117249806924664568582067723898320099363292807056972232983587418078546237330809911812710621614312426794941685612274178119819252489617238723781873947839631350920566598709205729843234747783994354298893946471308580545738315667354129923966467761181998443575222191776172737351398146997167480433343485146165170294472480102340073255517014099563109874895053870020199525964595015582882658001120356847662566211163371598302956108313767778820329365215346510347699342474048507914375556781806899818539459287898848416617033681979827574768955888631592686090562568985129325319077837317801444815993104372646508731697624771150941289932217586950882669009449298709005938858892942345174616723183138841693246229983569750276637389959194769278888114713340509066189169287625417553290674339793747889143872516758121133189526980760500100667842861417400152854355824828317930513788864259203616737901498546545717768909596995203105605745204296296293093603856700121390575987646935955974181414923515506847728701099915403480700745637660322906552095549479336952869675103157045737122473925610658940483259361595196142314156986640855054925675190067446394099279669079615310119226447311974800790401045383424980381929994819865353059688387720004440743976844960887058262096707272717735492777358680849139057474169004889360787372002980480096925586398565829500090824801600648664443222692860779909201529610028051643077119671535004567097478214287015464136952387391233335368260372571393692867068392048754709963953387319134686026335367452750916564459637866402295221178106760753081145915112661169730107600340794315900940604581063162768512489653741613975278459524747975421062619316356466796717946723517954972549067370544150436508076957363637385027252812676224360392334538542869262940323175230e-39"),
        ("tan", "f6e67cba37a77ae9p-66", "-f526110ca96bbb49p-58", "5.41943412376321424663580278654670761441833353885545358177360761414148472488754881810452778008432828937005847061055762561391927184634918740432638458951252121584556449620671233943390401237573160249798000183531601276581299356568063907165911616478892221548669707022275130521105918552889723320348041671065463787115695347022514161836663040655279154387576158860290774027381224271952377380296310879220330348836948565346271480812807740726233477119221496207123258375816624301145988940652968116248583317465472529288012808350259435012462158897416736656438496913013529220582153950860308741534888907968486808686174393943330747530354647084908810039283654074323153430334385851648883510650599926536992454316509272363045779418885080440952847991690950185913892509715645979599049369747831638149169726495005905275876557494592352508852459113673624842050398158358014586579305720957194580525068666238740919619373827455239680426765768662111662936626868246541991091342471421040260455973127527017153147101370688791340456337977927602441639565110840208446645841992844385284776708359590280369126805728461676129957598688878448294638990668569063650554502636719943037915225454267574096398639138578267926026231389593858732586476596723423836469060089023434826336979222768889891713328526451716514117276558763712087127857021517301626608009551700284281178487669032623344279546452775020762858291252632465083228546936398628080522990120327786933964208805596283784219176898293092643288959977882824801797368904925435743517321622861914229944484322681576336775419467481258775833891936737697240117937470960790264411556671458745791263610532784542959530040855245987925710342161811768536189119172901549998260066861212577815145194269694133600366169747227251353578590109222105432194238727375537128925385896582749844973832034437198229816555463611266992610681552838443554467144148472196503872209308695499424157731914060478861092366858270867135826003856748814368465060838074701911317457950836378034035086370213004356289789136223164744874922810297616779711444422593383654283135613404447230867433877891825495838150119286617450402670503410150922128046313992004234911751887653306914616583161645124988780246101918511421569710995496545399161694719694728309027277713773999593027297144987601086203153198451552962152358192235172854628661251633866675949226162187609780638800568297402752332301644683700948989401872510014127731541808965313207283928225594243925260845232257869210138997803120931019906591311109590912337190676819181009349439699479298497268525631778047021768554927437206416495617071881076028904737306580981686456114822364225605569407070259724505145674113915679471665181708819033981006005176214785765180485298355295870189582750220017632056433622170439399e-54", "-0.999999999999999999999999999999999999999999999999999989646602978516332025331885442528843718991709783008386092962107886241622933752156549261263864990774225554027375868847225804877374207058778141343557267148313458045582720881361678561202298745501332473283940676657816899651010695335661894814796333659837768817684968149599913027312429938795747597292312199718493545557458159605707056619779081734023350133031289239108085677537091034210381310871407212498185900097315022319904351731771191699944260012334425957925281721106186125965634099964704498604773732082805240602479458051700435779100222663734869339596419762220188092303919504935499297589165867606842310328534124512845054378573726673979471627457718882676280271482269468083080366473178329219572468784528677518009299959462867524451241665982118456651789265598256188680267952614895340030570363160005214206545136459920883129697027901941978723725869108696454204519718394944295578362847738469007060822240289142767599518891897096142600228018213799274756214882412272520512428752061664873059033578004297277000545080192414798388566221709021751485556388141474860183989262844519022156515828208697586094534564188476849781937008070832742759094495976351242153657783395705393605991937322970518395493335394909426171704887136731331178034781946500999569671253499212949310524739269373847605329966033215936847077234767117921364984488306235210833756653686073951093875392560383888034140402206974272319039461503219972475576293326373068301943349383521150858806343541161975825573459579076054106167078977583964403075153624616064693319913628767675696401121161657910909164048939368891396000759194546695100778315971641328481337132905649102561146737672507946066044930812065476971318138116603591689062905448825362172339091304822601038159027329630813259290011560517512269893324697662759201569237819592479042086497321744090082752347738164139849649369407546744761592574523763560058778569254057013528046371240083824432631413333498087071575211118810436453687738172749148753358564450258564129107714399522718288749874065492363569682625944177955366211247949666637802269436617866117340313859969176382074676514631869930703078311412606109111473458130689205403152274211915682022913087522306320157310111285744244040600498257859918447586124138242101760395138074053119252401973645031853073203756633958191230958552910361326633404130860123974306112999174827143312006630783816140535974525845150662284707390955024110253982686857526813357637864775154363269098051558656318441533833379103817849875908900081953602200145321069863406010125155723859439684816050661757691168686538896325035952148584001139654090542760351151553617438846839129275111383604268893220593514912499311958980471167875843735472300835953712150"),
        ("tanh", "3p-2", "-5p-3", "0.798160937059225037342931830177067868620290755634583274292498906708912595863664162528438442516083344899337539152184906263421506608677940001658214667921876335354883592596940465917105370925466219140312802368020938113575667301066437260580745615418661406982332671413167156555566056615670040898520463893294768407623168089302616541413867153163540712347792017117156788991604209069824661194740428309038290663759296839862719959347841361342193472976616075221913776236511614506095912121755185781248725391418985910012719603013346268981474215104463751805952131898573351586173415652518872953168764870804147997426386783794942927746865810086168781720952791852324219396544984983357128472218526027868381523977959141434774063515984004250298806981109495312652067717070907418244532015790478248520840019497700089711186918778305846161708167724779792580149889975194058012826286996945313532025405279922820075082861270039336422950300709924073298637074842210575736418561062143806194540249814852325060268528293124533091302710983639674256815898460925635681716144016468801230674309775696687498834132491314006182075686087202712609186619497975990054633733392069150088495710884251407547599370955502562730803269440138769043151557561963355736451748126894689689153256568504461669287937292519584182484947154764333976171773879636181433992590188923011031283870985902673060153906293339322523734874927515101203907304209413951100048516987760277608433830663359144873672782316894310152611637909719053048840281914766427171383229626781536730480270907309168788376095706437392387738289028397353363628311164966682376522024816995679592237261515211977053996285929249311696048436321698063300213179948239554095394917050399000893532385349363327760875100715665607218141821623067137920788907335869665291862486128967057519836988268658852440679708567456155221742448602088736597624524024815493546957000214523752361763729798028543684685946744444203888419026177164203677521910215852452352159131399518749293972800811502176262862949489316015787708116726142917762138726344504179869330733060856330804426971611690123043051076668458432091268881038567995432642563295371066963511154582352377860030405548257547745557267727011700169340715051635463408953609675323970913111618401529941368655939626473525422814002591001674485774702765867298667765850382412651084373831233174416042387633288737751414701355425848815988721110811934152016630698489863225891164807622524147096713108623718924020535326813779002297917369619596691995324224214455649398005428328333271875026943302792306015758946915424307990876179061022182039665334394626304505338928603836626060826529464008120493897410150924594118207004721425074558776005176879027331564644879564830023303710566603083741766357935084953073", "-0.355727122255973452418186839126939510179679505466256045745574258468485864731349042687793078626774340406098057209678298214534823774124578228724352854129514912414845352300798064308047195926167801665432411119358019114599942178553093332087063564457289490462230290981595764225200838209920561586241767653122322449049902270716291165266898591536751627205169996611403797739237339009500557147361169759449747454588543293600412661478840653456139730389030045355388679975952349719165347561043284576942228426148855086443181747282260856488763722918873642808773373391465903568108382867285507081310754050196367984791004685841333732455079303922638569152841398853670013992624447153145122287461457668967373723832558731512087183464338761463820500366203656023379439642466513437275988658555190209609555586308445628373473366618611348206193871726338893150382066848339651724202174801022601267845571890289876499638909948504185361508186952961136847324966916679685681710729201849443587227148289352541085293422210951908344975935975081359511580183112403548969653477541993173928711025165668234591111372614221807819496249305167289499779469732819225777036869613369409185946244040597625655133845168480495346267849128035210623698256911262168189934797260981631324510836064515209042848727473157402313415955334387388296590229099316485640480052799657638852112540215477930721937676146011335656012065697553118109607204017150296226431854955545467061069936698571400488042605962589563378509332150071555438784189284668854392911706277436176864338808010394927754638624196586141774217220192908926717050652512036582858843098726527186447397426904348487650157384175481802242199422369339049609827504326849662838306021160044850789509551165499406859752447445144524865789521171614833934452772303364397357872266285938125898808262328440418961193859404751874913717029236922078632958428650967726796961352956196101112946371717513903509423519759580448628810337953521213023388833553476129634617278118958952520449797649459671640669133428018568791924258156816978606297970949765120888558637447565241336201973922705955062748981120492326320110294183529332649068043987754285258342787615209782063548991626556983186602051303379920862314139586063093152204434877565370473008351748528455658128431756700636816817874311322029273988592440325692510700282655238999903911144441678514493936454571820100551962035353195531177162642527498722450891030470513028813709826376022094178738616814557216853482160908660674916288684568644955541585180747660220951414441085855698237075525801368261302309099238975967297703164283081981084751659398900522902652108713345334397646808508908562368676495203847861309963919282664585453031849775082116507878730745749740144932691564194698211899404753032270394"),
        ("tan", "-1p0", "7p-1", "-0.00165960203060817949068787429456295705453871319334516147070057231979331780494241797121686541955930276644171486464165852651485544376590090175387006051987266407328951442621268348328859778591618028624532809599111340148398335377277735487819683485206046148772736599047728728106999235755765459597466978478562878718115802122940608916094808062109853708364408852051872139512693703648936870528218766428531422787064151097855906449748357222892750664782431497246433655115129909456952513958646387832162808691205366232913152002034978034355426022622428778817551234820587087923301650306620165042307979293863593913951121161340421742238198283543017464573802038571757980874778191594289775060306426673719502165683083130737365691562422814838155065994312641704913445099778372706031823989659640539171730161209373679926019336962619807777692520552861848496983729449998382632958601283671079927678258812771917482421497364557948085999703527028720651452895689172301288287490156120914350887923469473990790736801281940899389395275768638796724550821378933305559823968447313335585908100630970752246214327244090851684408875291315696541285754117542972239680970286113336375464488265329209594433118916198541959791315108887594913797104122578107798007541415688701874543158738501444290118797767999957733631534703128857293856042020825091957033540864908265128695009020508329832565774810236129149957832357819196807545744799927840554484043107503209744107709022414522646115081939303685546517538504345799969958848625242870827677706509270295058672370460686633748027958889245810025355099755658864732951361568524009824170591672420283161298241377709439447157385650704472246715876237486940836550842040699178467742033380787706702268813524362973248745147860334582377839203370308598704029352698115975254811940507078894034134845174387743288632259089783230781485896421966416533221163100568054063478870878859136644577852337124617385087618845829663028739810119429666206957338002704340715879907205722303731277011519814212648063448944515732960450275508510657479375172042773437570006296691690463379866398032189184902151070461821513332597910971819911948115753421809958013123170150158780592934991117273482892515148198961286196617219050958207780525123660884485047384688160135716136715015106555624517573751889701279376554143951830945043082861206791946821609650738169940172083837879111071265691638357555590777893440973813227562635284280452586322985419058441969063996687330813767133886248755795025258816227345607165386195503411154243575245224440117744689576266772379253519129797881024933594315341487871837291865856481808243121739687191493282171284159786756263204909143977117944087737213385477167307131318811251849951735892077751221279633378782558553140267437166952753112", "1.00075786508734419653353243274821407981319851002406338049184234934903578052966663497951926389345287848809911853883973492472077379904396302660085118051183274169688674769682663719129518169724830085878409338620871663155396658689927086523839062555470954527641227895013309083679933827903903157519078816982045226644059402463204547339057363215980106162012370598728388555135577303583617856961024420171823009307976166947493623184852664491323107532049036472555913637208789266324461874737711698820306832489119838775725552896397473602317711377174748423279347772660182743283713312993851710527634302086846427180732326284031545223431457985305883372592368239549701800510555638877863379480568948502015734397582676714357746135258947062608708071713526161391333105450400988761200988759711789148763605187970492163980614238576959383073870728965732877389990036113364668011371276492201492109856711287736824481335336877858787144252926483890033921003667724393818420930340914480086517192086949690816192815289374180041302271251078458490845661884226480890491720794480793665669072820453505100523485798973980413312803812413670239730044665175071854970590195707329609258730304661023449051278884821723369139384851054178903109720828188955008410903068875076465779676931864584524934438288729798222330970013837684725884585151563268975107188924247401921146883381560858650346857021122019449004893601779216379734423013274601589574792152724924847860751267451201856377241823538722872427731750557392325816541232046784512865759459917570233715038311417490525654314754501185128308487884206204386756519281126996387123693079162979401510335151051260953899182875867470279552879133422757932617164098207210436104188509293481958277380553599699150936260791864221039148172499550450827022086165802750305433281329809027811399907555491562256278089405327511661623495437405506507013812564030402332837869030624635752563402727403602416848658043606080767740766358932374320255552587250441831465056298618977377344602154529393152124777816247683795979162290419577002915888822814229921954980432315478707042684483454991341667498733806177076045688298687382527205073189261435055283573104480213514195096299486334911744413913132529657756338711215832896963059024087850342972360745212093877683987748303783933991431185935365667948349055751799412807148265914803235319175585103785470958440713313437412712967195715071727449125083652000146306281222818868862008754876987861101354239331914379418490795306755957797620007072931470670347533805841829355759631725090388577364837104829547526841056297974920621313647009121956707024878713191148455479365553968642815729715611065139542118014215372085645558494680182788618369547508568616646831722181440419983407741143418049406362359332197838649"),
        ("tanh", "0", "3p0", "0", "-0.142546543074277805295635410533913493226092284901804647633238976688858595221538538059106058347766911365259878245507888772472019076920087846369343999408979649330759312837267324176849873375274245333747976187752189532112586806299307693271291571383776659064099698457847360557362490157466292617484083735829816725441222845483599610030924358943120617137741246743485064206120323431719185038291162004075049365113737263939948169167841734994399163071061240284060244786081950663714823651981748237102807926554503747104781888953826902503679115208997074574159679600329301545300392810947536216950830874703569362340527236915279133247195957508214459887379900750807529590997101565587843429615328348638398486209968717149420636233727707192311456407731935794918153194381138001102900654657581918973027392154024648317891925556280182466803807766154813933311342380274325685544859482188436757012451017957011782382362378289963065753221839773063271916003206144904359597643902934844377665696205488880014640565785164109545166875193648978817500967529731945283890423670935029542746274516610658757756657390577429214749443176228170915802924602556636892478428358904602660080591054407546977758101020376853508868753515833918991258278049912433449889092399107825147821382712831533302787948872236789624413075503134292853928439710361253717266757345162441359242223025420871707083665833742313891414433867137093690277237561991929412699488192607829809993896637041226820717518481991801923517479209525377088908456168580119305890966565237637028934255200830973036564554274007438505001001432481545728780950125324124058874540935373402934229293039248003655139248145015078122551296956411050549167852026161648720387256805516502492002334416877049644168351452980595560813020207676624076577427314778778942114202786778184721360437485357697269741145249513620505851829681682414191618650069596471070260315936470211422276533573654476748018944521610568658481253988798395344998493212925642770515500234444508624153161375612681917157017795134604973265560201659229006472558701798505880918110971163802972935906338123425262971558472229811117391555205581079325640421699626495420276059446873702398467270956562543249857803016856602333111410817910125555979128003781622411893042283238619281349228831908133949785151813793574451303927717579302196173339762999393606251028967071299680118637257821723829742496183409906097575995973637366336032682815277221670995212819408987427865108577748811896521668233322261391284432016244572323882361191923473751022348545751339502195819988293936848384830700714070210338776707089716616677899918186694209698636683496707447361425724696787513957308655955510302924715040395499135447815353516527566203567454937661536435850509698897836758806926943715713"),
        ("tan", "-0", "2p0", "0", "0.964027580075816883946413724100923150255029976240934776048263217413107946317610202559474850045207689149461612218924596668273067021913020334547485956472914624336751594506316812007245950658562057196950550739835054131848435723071198131843339062586263829722006258367481879218305929297275083923803297700207172355848002085069538548200137017285152482271698397724274243225820094581477045763389460081245553195088494466588107653125282885083584614166508866377219424512881516524242694642518573784048608452648156177106053134205210471114883342823262944045332966852370212324308322909318660573723962507323206215636805872310682321745819884493948112773945133947509343658001366991702682860540736937439779733085285610020708699389336369816136102710298320299050304441942727946762826851534920285754601525760697971440898155981605495458845118219474687093290550710774635760199807987695543647627582934980179194919649325139630985930581182004567267318605267923853148938139124128543698177230449145457068846642019367802941658723005092030840741636850677743970934748808007225199106618388407601044449810648665587951585482832891494744316810622497789311676411437561313245586933425936487683794125962166826718246919971932988560064523455951840011854157486455917681057988783644248630787775931070978130112136801305616010217765632373523618433199540012251422259709883708458656882956271211989900847298705124241768206429450960775522903913796361587134328889830485892328523409337458835254601265978724500797547995535169795014872736557332791071660157796873811851150575386270257649837008861499979708726385247009227614833150792127986917066649867945587486120391692194562984212545780221377945336902812531445678799827521846024096467382562343996630566971131569492724379028827593506848650002916622420359298326034647429510836843883305804327889153555121930373163583749787956549107480120181302581814329007662928161101470946120780534452150624130158218871015768476087727034132851355458094546389225528180828218117089931172395109886639967785470456187455013983565596656727612532419468409175721824570992458255293898700349907236776175505563726216065570857132824367902853368206926643606701901908916481214094702479601656168857656794786108211945142943946148547473615237835152364824126538548873467229201545132322668598849627305092219630583070569386687766597807568639620752918903375945354136212453339852919929790981188063064182263068196678300140813247536443667936441865367447788970347312475208040513529887753146823542523413987024530721058729136215804456902284841753047687912322741373515956947935122429223391928533224236332146055091632626445522636419528644214521931575313449651051473858846786056107695259539713854653616544830411463746933594574645646485027976091109137223440"),
    ];
    for p in [64usize, 128, 256, 1024, 2048] {
        for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
            for (k, &(f, x, y, wr, wi)) in refs.iter().enumerate() {
                let z = ExactComplex::new(parse_real(x).expect("x"), parse_real(y).expect("y"));
                let w = if f == "tan" { z.tan(p, rm, &mut cc) } else { z.tanh(p, rm, &mut cc) };
                let tag = format!("{f}({x}+{y}i) p={p} {rm:?}");
                // Row 0: Re = 1 - 2cos(1)e^{-2097154}, just below 1.
                let er = if k == 0 {
                    let one = ExactNum::from_u8(1, p);
                    match rm {
                        RoundingMode::Down | RoundingMode::ToZero => {
                            one.sub(&ExactNum::from_u8(1, p).ldexp(-(p as i64), p, rm), p, rm)
                        }
                        _ => one,
                    }
                } else if wr == "0" {
                    signed_zero(x.starts_with('-'))
                } else {
                    cr_ref(wr, p, rm, &mut cc)
                };
                assert_eq!(fmt_real(w.re()), fmt_real(&er), "Re {tag}");
                assert_eq!(
                    fmt_real(w.im()),
                    fmt_real(&cr_ref(wi, p, rm, &mut cc)),
                    "Im {tag}"
                );
            }
        }
    }
}

/// ZF-40: pow_real(0, +0) returned 0 because the positive-sign check preceded
/// is_zero. At z=0, the b=0 residue is multiplied by 0 rather than 1.
/// Exact identities: G( ; ; 0 ; | 0) = 1;
/// G( ; ; 0, 1/2 ; | 0) = sqrt(pi), H((0,1/2) | 0) = 2.
/// mpmath meijerg for the nontrivial identity; sqrt(pi) from FLINT.
#[test]
fn zf_40_meijer_fox_at_zero() {
    let mut cc = Consts::new().expect("Consts");
    for p in [64usize, 128, 256, 1024, 2048] {
        {
            let rm = RoundingMode::ToEven;
            let z = ExactNum::new(p);
            let one = ExactNum::from_u8(1, p);
            let half = ExactNum::from_u8(1, p).ldexp(-1, p, rm);
            let two = ExactNum::from_u8(2, p);
            let sp = cc.pi(p + 128, RoundingMode::ToEven).sqrt(p, rm);
            let g = z.meijer_g(&[], &[], std::slice::from_ref(&z), &[], p, rm, &mut cc);
            assert_eq!(fmt_real(&g), fmt_real(&one), "G exp(-z) at 0 p={p}");
            let f = z.fox_h(&[], &[], &[(z.clone(), one.clone())], &[], p, rm, &mut cc);
            assert_eq!(fmt_real(&f), fmt_real(&one), "H exp(-z) at 0 p={p}");
            let g = z.meijer_g(&[], &[], &[z.clone(), half.clone()], &[], p, rm, &mut cc);
            let d = ExactNum::sub_ref(&g, &sp, p + 128, rm).abs();
            assert!(
                matches!(d.compare(&ExactNum::from_u8(1,p).ldexp(4-p as i64,p,rm)),Some(o) if o<0),
                "G sqrt pi at 0 p={p}: {}",
                fmt_real(&g)
            );
            let f = z.fox_h(
                &[],
                &[],
                &[(z.clone(), one.clone()), (half.clone(), one.clone())],
                &[],
                p,
                rm,
                &mut cc,
            );
            assert_eq!(fmt_real(&f), fmt_real(&g), "H sqrt pi at 0 p={p}");
            let f = z.fox_h(&[], &[], &[(z.clone(), half)], &[], p, rm, &mut cc);
            assert_eq!(fmt_real(&f), fmt_real(&two), "H at 0 scale 1/2 p={p}");
        }
    }
}

/// ZF-39: capped direct pfq series at |z| < 1 returned its last truncated sum.
/// Harper's four shapes at z=3/4; references from FLINT acb at 4000 bits.
#[test]
fn zf_39_pfq_inside_disk_converges_past_old_cap() {
    let mut cc = Consts::new().expect("Consts");
    let refs: &[(&str, &[&str], &[&str], &str)] = &[
        ("3p-2", &["1p0"], &[], "4.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"),
        ("3p-2", &["1p-1", "1p0"], &["3p-1"], "1.520691992601892695062188509760811648403255461894352854041141340520110245309821506056916729596974693430518671353748608720025605317633562192764714045310686747343910060176758784032202332210619240278502297079489791974306504093342872653238506615028626506080057441843799774749555134524212996496395626165787567959314251295658193026918417046641117111543724982470764424398305334336694894234688401590531242787625205050924795804471218249842141094818721045335900548168238109695830354938939257836155711558749673567047886682430451244170950920142261573407166721329596769036845706002341768040460801764408133841080847142244341502796346642897478713019673722602427907349985765221275866918005474234428475673802485133372477522368989764900307979997845520681649571521677114746368677851177953242412568046642667856544835363726179829777615338324345651467947093162116204687550084117325472510293751289995269230236680173312871574270759270381017745867243400989377578215819130139486805237400198062105748629284738521603603815536182713257170572727949252332576266169041931254998980527822353277565888545251584262350516725426073114902531166012831990340227267126348760100622742458846371"),
        ("3p-2", &["1p0", "1p0", "1p0"], &["1p1", "1p1"], "1.304625857240408138324088888699415330348645699492649669473105599334189018390153076252195293468939089267470754613660220441903569239408998562355562837288542656564358956605558739071695909650558292965557983032144172424175925626881086846307082445738852410219024637178486077633066457780224933116724966005055615609569622455462136794728763404665651689879667729621699707475479264330820328917222307344677750227488702330924746701962947484311097047386385228198532939481298047732350590704084326109871820578276236088449160890846980816693201729480735633679070859397752937715842078154463194107648856275781210333439737176197778137987385777473600661014542592810722655099905874096775403996911675214386957340412256371785863595684371321032312520792432325322125876395125833315431984930635074644600953084529236692045836133985977922483488029233914525131924210587730210336697700583420614582507144070717561215373782724375309801048988634730938268253602981269027784175848738350617118529253262746084532114546561107640625706648343308141955159721129757564605041305201038963512508912844828236531300531578279967055755939442172935677666843528679214764850618630969385856212275119403521"),
        ("3p-2", &["1p-1", "1p-1", "1p-1"], &["1p0", "1p0"], "1.151719620469138368609328506960781742789965277542611238701709933668888216914227253637277667534584081748303232629587828538457214063154522238154435407103788348243409575052341875958704821491998122215884817129457848392001864626365430698056181429339018858389353580690730646858481334417112330818878887236635835258703817566885497452188287500497298676640309542685287872203304804666807101454254972019273873923195674892481529621675979744680982790102906159761030621797910047520371463589653154637083136754073001312526915425191426169207705930068977832883738892689285898542387504271355675899911677379031870475537029594893930197656523513772661597870712809556674096761708582864578529798917705008093637607855120155488753777830680545776321203754183528379933296901700851249353241869426481418369168484760629487100024949851745974083476267000243184186818725416141044991205903584198503943443593253612532810110770999275440801978942515650994317344932810622574819559597605042228919786986903110851927385886474509277679306491742187990077446182141454326702455570830542650073974675694358188370023860132822834308821701282649275669157369124036792813430405063690775390538027506921145"),
    ];
    for p in [64usize, 128, 256, 1024] {
        for &(z, aa, bb, want) in refs {
            let z = parse_real(z).expect("z");
            let a: Vec<_> = aa.iter().map(|s| parse_real(s).expect("a")).collect();
            let b: Vec<_> = bb.iter().map(|s| parse_real(s).expect("b")).collect();
            for rm in [RoundingMode::ToEven, RoundingMode::Up] {
                let g = z.hypergeom_pfq(&a, &b, p, rm, &mut cc);
                assert!(
                    quad_near(&g, want, p, 8, &mut cc),
                    "p={p} {rm:?} a={aa:?} b={bb:?}: {}",
                    fmt_real(&g)
                );
            }
        }
    }
}

/// ZF-38: p = q+1 with |z| > 1 summed the divergent series and returned its last
/// truncated sum (1F0(1;;-3) gave 5.6e106). Now continued through DLMF 16.18.1,
/// qFq-1 = Π Γ(b)/Π Γ(a) · G^{1,p}_{p,q+1}(-z | 1-a; 0, 1-b), whose type-2 residue
/// series is in 1/z. A complex value on the cut z > 1 ("C") must be refused (NaN).
/// References: FLINT acb hypgeom at 4000 bits; the 3F2 rows (FLINT returns no
/// value outside the disk) from mpmath hyper at 3600 bits.
#[test]
fn zf_38_pfq_outside_disk_continued() {
    let mut cc = Consts::new().expect("Consts");
    let refs: &[(&str, &[&str], &[&str], &str)] = &[
        ("1p1", &["1p0"], &[], "-1.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"),
        ("1p1", &["1p-1", "1p0"], &["3p-1"], "C"),
        ("1p1", &["1p0", "1p0", "1p0"], &["1p1", "1p1"], "C"),
        ("1p1", &["1p-1", "1p-1", "1p-1"], &["1p0", "1p0"], "C"),
        ("-3p0", &["1p0"], &[], "0.2500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"),
        ("-3p0", &["1p-1", "1p0"], &["3p-1"], "0.6045997880780726168646927525473852440946887493642468585232949784627077270421179612280416627373533896187408048270221751902653508313448027165994173438210206238727144690015912459983646071251381122280044319220091497469332127356585458688283611427866798297286974414999287275229223056868597398702087317959411112567429801190248167282198343147141804614390873928977991210702349886957688178551241231896660173652790621673460689983331298260956878970075982667878400702599783549060227893594090669452885393803212154420897241331446374714851461338008455225356264791889782311375232958758284896716470219739074565714984375850118133417003093235782157786195519323018639978218531770363284780013889055723765486785185268299577786493847714743672196825678189961116552375786128819039931786304379530161150467729610575258034098769650530635380629428728721655202259996336007828552762463765857735236635146619485591750930258763056296572625479975018068757679752680943746982502111862893551782927652545383541463692924141550952115989104584610315360424205129155898390306290955332495432707906862503242564326635213053480367002031675287723423250497816240086131639044395156737936436126158937660"),
        ("-3p0", &["1p0", "1p0", "1p0"], &["1p1", "1p1"], "0.64645847358890298435909057305929714707419672593619280861285248891599084277076125052968099899651943293867808709432098676464597791875776197604059423924269852470054618235065395432475892652840201196997698301538913584791541676127135013688409509865272945659425767033670698624131312412129161969960846087158589435287450547059990665142118824043918087787714940323244389643008611936801658700899192735025183527259751585830109650862697199569369815077211741696902561847682515626457536801459929468380538495257245463178996496845913989423162178560737518356382929104855272310792468329382752946745839164138382414285125048072562155314272909742639404506224827367401353594723875674728821327233671795938994569623037665191092906862179839077844799954891957184590322169721774728939171989249851675528353955293817220599021845816472736494987858302180836534204971198456899004669582788448834538201175527330327484521727567737095525918330309696623316336883523830801873673046461771428450340050786351244591239707749499031416437117806752971869383158922255628281869144492564268435847061312646543446686843419834396498273809988078520289845"),
        ("-3p0", &["1p-1", "1p-1", "1p-1"], &["1p0", "1p0"], "0.81231699341484982485444599008022356349830125247211513049834599409782375098049396493328973648958547582245440606836297071979745908280145374413812192293859710106355619598627687993558434302072555674806084219216861984908463804236618247906814720881826445876598122450803572293139513632716139901185210184589389647922328384645610546005692510871480072435306328125513016885164237492301223468680367563236798533356533280116847370932641124539606514097443112178139914723589610254773823684176282837088789551916233508120042532506092086706974451293698154927303683405270669943672552068004364472801523571401565253480410285007167648701136478694574262141173694738050972998733665341550511775396089476463710414577989333216498100625286805920040769612554453269944711419913553961433979039879946410269899397696836429190695156954573342078633649257636240791469175627611841165137648109811954173486070104648384711315075213291514699684431385517597316784648517934227753352757018004242501860777847420212411860805460063757846389999107725210117299391621308097194600239643229614841184389250454850409608103203832506974586225103300396559417"),
        ("-1p4", &["1p-2", "3p-2"], &["5p-1"], "0.7009889082827711916836145772286121661516897886053581483272092792411825389494568706646750857530739939818165795748903328930400295234865651007716832249071617426286563497475524359729054873470589568040779240930648741765191561479061348048892943381617179954836360651323207190243199833304072278953949343978506752385139673634786631523934082656875012038435913914264512455855901675442091737865685137087072133922488843889894377811176531627556683398499298726948647795179710838442931643227622295995035591154939217219912858509062585028743511620806490989898913512065988777997668870990575076878645687756381958522946148275066076132311723720355544703508241514325052041436081201313077403576505370619898113292118329783392304388357300678040098857646527805014238559395782485637624529279803633803952240292558165852946185218015769940291852184969173910912620451129772608754603256906370601949620412994805889285313692550896878904668889965097654268778926120658919146846722819523020629109198896511016030945115031358682448930808798086933224350980522237406409587707190592459969659085985095953368921170197770562222933094372199663332723786484670209132920689824939196840069662834294539"),
        ("5p0", &["1p0", "2p0"], &["7p-1"], "C"),
    ];
    for p in [64usize, 128, 256, 1024] {
        for &(z, aa, bb, want) in refs {
            let zz = parse_real(z).expect("z");
            let a: Vec<_> = aa.iter().map(|s| parse_real(s).expect("a")).collect();
            let b: Vec<_> = bb.iter().map(|s| parse_real(s).expect("b")).collect();
            let g = zz.hypergeom_pfq(&a, &b, p, RoundingMode::ToEven, &mut cc);
            if want == "C" {
                assert!(
                    g.is_nan(),
                    "z={z} p={p} a={aa:?} b={bb:?}: complex value must be refused, got {}",
                    fmt_real(&g)
                );
            } else {
                assert!(
                    quad_near(&g, want, p, 8, &mut cc),
                    "z={z} p={p} a={aa:?} b={bb:?}: {}",
                    fmt_real(&g)
                );
            }
        }
    }
}

/// ZF-73: 0F1 at 2^20 and 1F1 at ±1000 / -400 (huge growing terms, or e^{|z|}-sized
/// cancellation) were off by orders of magnitude from the p+32 term cap and fixed
/// width. Since bb3e4e0 the series has no cap, a proven tail, and cancellation-driven
/// precision. Harper's 39 shapes, references from FLINT acb hypgeom at 3000 bits.
#[test]
fn zf_73_hypergeom_0f1_1f1_large_arguments() {
    let mut cc = Consts::new().expect("Consts");
    let refs: &[(&str, &[&str], &[&str], &str)] = &[
        ("1p20", &[], &["19p2"], "1956309906482646147095375499632301679088436267994987800942366714675280268240935497342562189197194099259062168853090711562938250104493846582456180291380134448264897790839477268296764543158157033258060690396932212221279940047179705531962699859229682043225981865594420811631854760309330542179564811443662151433259583851184281903798708739397772827615042804072143666853773468232918613702888593619234114460288980094303462343643226780582067779382168874975381192010779104189381895881967572831255592462658431689179881470865621895624573634390601379653805129899957654186868993741645615418653009299063839707575543206259691150001071044303639834504847100102793823536859930690988949006624359203490765643999162625029995786966279694820300316859935957120281547856.37929305810739146662111590573476230000369232722266633019972538586483625631856801477050590935786719362090896388045523191852446610654393442924433862145048566"),
        ("1p20", &[], &["1p0"], "240088048317073938876169624847079874459857769219129691808124000367058402167497016800284020229323980480462140706885182531789148508024184272734821284052717486446627533403978618792031147454160161529522516536784342685332242085579177513877195973740474805804202383036429055128408316307841532736732062089759773103377113604369647034927602695199932471591395954663324326439675797732841371492721166291950032598852803326707182004755934009929777637808012934184652152501710009489607769768964588455793642760800863780601742721017021709970108838434051752736874723862125200414987182457281589389311087324218714352194841211042069348227177605275150003745131745325327675543540293558177101327090951974472012084667620849184660064615620949734472294451043544992563452050461360813380898366437837216060907491287072711111102306162505643856919091210546263577839570896716535465660545463548689446045666927490861154907790.341427478900"),
        ("1p20", &[], &["-1p-1"], "-27873197037145286525027388670447242686571331842145122211320986904563993276908738802613880472677208631116155468356995863031916502660509280176876178048537614151918838346543877537312466214703494441997004128173747803644770140522358082039470901187817521569648323198271753336023712944924707563142813829388657256650418859456476646169802975682971590460725483825979109265927308987904716456079194938389909571861091049245295017694033142799108222738410375433242922108490463099816847957830508811096176580478933070286155974770696979721269893120925149243036203216090140198475065915330676143340612872613825364051158856529404115043608640136643458283946316444766700540500993293348472313818848267161523620837469372478470167377506821520073937346390152786932544443796640088425170564993959361449914860620757005337166182841568965110548998955949120212522887191220340916976296769876349322561028312591055540054624959126.2419831"),
        ("1p20", &[], &["1p-1"], "13616608225278596250624029638713845963151603244819307382179280363734241952568997949493835111224821021551614786691253474856822912877630327394663496848332981998983311356396618240015860388228380284317051357192842112186013747201933601387137714307678320258743685001598316236455160207584126801730734650409700662750571010970433144196288703313615823380911325757683981077639134825551888840292718582506062321378158792987442607569141740497854529916175073489615496877621134880223179266160483053784160518064940434922401550938298475682105468061028407055708941483190102686113857310860125131089698521061956699585324307049049396699369145157129193104028488737062384240596479381215667959852881420205922628645563933795051376344654040801208567340688887536361770612504465114032814149972623039301375115105401565870623440567449421158060087423521797856630623933180430345371908534380239043752334300239890346875732759.7099374612"),
        ("1p20", &[], &["1p1"], "234403736242407179412824898232848247658740488137317173743264739949200506015035106328862557794350452828352013635294624644755600543266865563758062543643348859749745865524332500484904886912516699825823392857168700885629252811373220173847145828286381998056531455928176570448131593943515456032185323776204477074867501249063011714570087308024694081033031742693841152477311939550831819536785935918005080341404337458737597251388535842398029388950638680725097626250688415041862278502062244760659846914711809895294583855677494092325930165324198060078182696377001827223349448283864848587706389102455259419284483720143990979440376106149036717281868118314332856711566104185265790543832808402657657043355253680083305182194202330228822221544267868048712846197316160053639375515335421955707602584386417872056901838696619682959066817692717764052040734901273449584327352132271488330015699921009881788482.855107758961402"),
        ("1p20", &[], &["3e9p-1"], "13825905832883644284313406003255377940103901922894401819585065176841862991442471322023961129454283775983342441736920336328972280837179793608801623221671979917649104301197219297526466421727707841712086860052024922981297582668751980029774396867776561837736434318425769195487934248257043941100513929985299097358638912960826596437951513585869094913288626401864364768861159531370957531565855342597034552457685215811808735266477049580481963352057424541523723233189874936573479959517382913363762086.4395209683696240829942316647747551297629584666686415625769852028028539398363242675933961823646352702706666284520032536699333808738821306057867456172974787236366960557611884043601450659527493263624842876952167294129397256901630296595630947111242193316281396731976543473806451556479077170180957292392037590118842430545248461255739559572309265670574538040083070621749110959923116516186488362488650470285700215922"),
        ("1p20", &[], &["3p-1"], "6648734484999314575500014472028245099195118771884427432704726740104610328402831030026286675402744639429499407564088610769933062928530433298175535570475088866691069998240536250007744330189638810701685233004317437590827025000944141302313337064296054813839689942186677849831621195109436914907585278520361651733677251450406808689594093414851476260210608280119131385565983801539008222799178995364288242860429098138399710727119990477468032185632360102351316834775944765733974251054923366105547127961396696739453882294091052579153060576674026882670381583588917327204031890068420474164891855987283544719396634301293650732113840408754488820326410516143742304978749697859212871021914755959923158518341764548364929855788105859965120771820745867364145806886945856461335034166319843408874567922559358335265351839574912674834027062266502859682921842373257004576127214052851095582194482539008958435416.38657711790099"),
        ("1p20", &[], &["5p1"], "68999117683922831424193262651469680621169184248427174502107666036071715827584872343260828954728692759744604644468526056188039123236289921735429889944344555770719228125493257302229213006387646533122183078062802559403259968370797335950516153926164444608789968485483261140719512136389092412495135948311924546620541869064566970010727562801958071626266922768772764386476154987645229787185629866522770628000272996201680543134672532424609039097977476017071803078010571138052328915794064327035336804671473027433474705598954600040035634922386972001708800485119581451361682898367883551016984468474568179016782984341345189982213589357994300056294375019084477820442897621807146813905283600787645065570330363373728521327038777563878335504051587537622335591983857130732329884872055289500938089609672589427756846966848598410533277693205229575899440830488476229811192851947735882138.4052298989768506738686231022993373"),
        ("1p20", &[], &["aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabp-257"], "65336037648205360149077459619841721749717299881046706550221897364406510746924651097379628241495603745597420297228830902177984356888427252635672634204845894100656375417882889806756582850873291899702180443025152528492759718516084242544017228282758280621873701407619459059500856063135075556974585840900483330384765202394807802978487454042847038864159497087089600845978868381829267815637824462311573814418685703937077557750110258292109502582922689178662314957244416214500076338474624301725569999637844423277494033817954055144580301023970486127202053806885557287769449589136377367507592289464105458056148932302949325252434143498092056556093611029509909563370813560910622480044907516644680328080300377317554436356437986499954195355927277448542264480938850051348312391274055912038225344489600947108348352438286582277338129466725844044009887678198471141757969919801595011750880606112956935009954924.239339101"),
        ("1p20", &[], &["aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabp-129"], "65336037648205360149077459619841721749395269135263389035620073705948124151023778303555110605306784805810987895410412112327149246693196072945674283337066305802206480682387775672368904564252684717712080467957441651704048784967323601958918057129219359893725747419013409419684992516332977877021200670247585325547474321758909347681142470374634337549307408088852402844202902857667663435390553588668733198515555356356893617907785018300462599782594410296308108685871134619379881213565120458847679991514932671827521907952117174078437601069371886456934240805366616022822110829774047975828299397707538042269090204666309691740225729057213820299090770559777980343375727597577985066153287386419374337386095829139983181654739109608278650567758487691648255209236544197750694195798661295404285224984486387459984741698963375545413107230951375156822015388273268971287173744279796660767526481647343641768450989.9189006613"),
        ("1p20", &[], &["aaaaaaaaaaaaaaabp-65"], "65336037648205360143137040868511042209183755713526026170740142103227626184041765243913600113417581874475320840183779819976514725603941469743941539418330760886759330836292648387352596973448049946948892050628760803491613368110873054211640008060175242668053618909214343835756164403512412182087889581447666034993868682385091800375550011462289316738184347042053983556371700910763908701778251307449046612409257374992700418667086206347943626829047342931896290034110628122573608307264983449929864562705557567521328386220541479621772412033607920342251047484020524403892414042574214187223158028798248592667206408640169139709537203071868955472788423172406053673698034981448531751076910276797451603363428983521134811090910241907422439739103893290197022045101207285980535092222517820031070266739078879794104269777385039374213310876471992860266942612421323447896327079056046513951371017607291436444962984.8726820974"),
        ("-19p2", &["1p0"], &["1p0"], "3.72007597602083596295969580386311833735889229237678196712061387666329047589581571815711877864228149660193561764231106980024798564205253560026618568828390755743881911602284486914975858551028166117416087723707013450821757552574968763804789272795294006197962264770505210979350924055716149816993739806507943850173926661166690848203558527673492647359653340123511241557018885560613215556672375974832799430815805582000240564091832624500149510415866202796436562644679688947009387929162532550106322217632016527933498038909405956058413877000844162154012414142500778256654159836909760099849613602567411940341383948768030536756039333030247443746537738391943980381211393591259336532973137296481439761663365206733489804617788594010106836121728411922670484004250047645085883371641091044714170300911306279943414442310667065160562248063073664771895620932001059174369711374912492100332939004308975375654231729998332149e-44"),
        ("-19p2", &["1p0"], &["1p-1"], "-0.00507694375197056065496839721430971775813500566420420126367386426823732564963126009685435003569582978095907168627634084489092470470150297460963458359190552816729709464714505944936603656325808126317232261467261115018091058653598916299106660774483216791690672727175719081510884005646501478968470797073444927856951133056977321542831374844663175607231996692241590285420664823636655399760144938164157929816140752321187406212629256824239854793554834545780969647266479029163695771960441807077269895510293466856253073885910844119463062170297993761418679779872349318951793937600092345749986235297394388531093878261489620746510210025149934905864884314667012678475993436262817943444581672847719065438553998529514607784340113779420719814433445311886434424041590755928791870772533420773264940701547345565505062851858035198579548699007947046874100677310663865934608276459214373670432047643399944560153440369084696630"),
        ("-19p2", &["1p0"], &["5p1"], "0.0832752132984748799999999999999999999999999999999999999865005882982155904576118558669415161773920516494231335977127163643642515210692639219514473760628889051308960306995815899088601097021197588137540653743553562555660131917763005388293560449002899079312050086821198958965800619321595334990518068688043711030883453360079069039813136678661723554409311699014397335648885493195831225003892681477842988126128995935980240663788986807764676338794728206252673742545560470403752704102355777221385745660290472329229100147498634474909233308265500188217417793665693842343492231640554766665522772313933670437551975155969317586225338478382186254966572215900337555088918592671057078821968446829983807613256385492331368399266009493603811958914707937852815139287598033780551219700296875005612431308147193881501334764537742710551234642098880881694121881205305177134133767174305135394535171422871828527574517076195455647"),
        ("-19p2", &["1p-1"], &["1p0"], "0.0565616266474541925299391880157635346285608435047457839961044034794656034403366519033385290146154249331512103304181670933706421476703832381536629281151243176619761886290864057203098556957862620796004395455791684282151599203775158513551686935954932924457192296160032034685873453166297336077318694121577634475096683854187482403335058901529082691740052793697037167607678183953741240078890265637910356706653189386078695605133395343280052090505755011792470164650395111732136202913174871484769680259149928436511358889707559922704639853683045513379022650043035326087386690216014496823674181584354363412342785909807478606010819934177984510633130803408156233129847438659521972937029119699046253312074801224996994117923298920064464823126167252238262691706106651659408829086764261067620951923173342370801827573154559876812720439079180438397642357147242576947822658024023611212651254009897569048747907649121129608"),
        ("-19p2", &["1p-1"], &["1p-1"], "3.72007597602083596295969580386311833735889229237678196712061387666329047589581571815711877864228149660193561764231106980024798564205253560026618568828390755743881911602284486914975858551028166117416087723707013450821757552574968763804789272795294006197962264770505210979350924055716149816993739806507943850173926661166690848203558527673492647359653340123511241557018885560613215556672375974832799430815805582000240564091832624500149510415866202796436562644679688947009387929162532550106322217632016527933498038909405956058413877000844162154012414142500778256654159836909760099849613602567411940341383948768030536756039333030247443746537738391943980381211393591259336532973137296481439761663365206733489804617788594010106836121728411922670484004250047645085883371641091044714170300911306279943414442310667065160562248063073664771895620932001059174369711374912492100332939004308975375654231729998332149e-44"),
        ("-19p2", &["5p1"], &["3p-1"], "-9.80608357023662418354482673842597106254241539953138784977993490761720852483580432236628611549346122394279234682276139942217460061469142809856112643557902760687638219414267235580538014343021849525252283780216608949372778888783196462735386792769089040632723083556283198628165336522561191561428719933585092013680770682280411517585554909727605442211774863435257321655128269851843352388166909620835867097254200274737403769645506166253020121621386779744342925378545215678102959794504199319763499961082059202244856694464381866307669897742495402568917430224337561213269587586665983932698917000018907125551350440853192470545844317925522903657830302595543151076800355432319118365825098799535291127799922016044242232710866942325324167990624129446277537287942369095046472626645662366752908299473997032611239854339794263818358746524896046544367551813577786805525416297548400101617507630425703875260419710084982577e-16"),
        ("-19p2", &["5p1"], &["5p1"], "3.72007597602083596295969580386311833735889229237678196712061387666329047589581571815711877864228149660193561764231106980024798564205253560026618568828390755743881911602284486914975858551028166117416087723707013450821757552574968763804789272795294006197962264770505210979350924055716149816993739806507943850173926661166690848203558527673492647359653340123511241557018885560613215556672375974832799430815805582000240564091832624500149510415866202796436562644679688947009387929162532550106322217632016527933498038909405956058413877000844162154012414142500778256654159836909760099849613602567411940341383948768030536756039333030247443746537738391943980381211393591259336532973137296481439761663365206733489804617788594010106836121728411922670484004250047645085883371641091044714170300911306279943414442310667065160562248063073664771895620932001059174369711374912492100332939004308975375654231729998332149e-44"),
        ("-19p2", &["aaaaaaaaaaaaaaabp-65"], &["5p1"], "0.441200984848807549436141185221275552456659680422546732705192977513555833582130092751166585446256647746388292288318873061075842366953100296653729610162044800375134751155470146492790418499265051231123385745630404619838515284830147090610836332064110309912526577580265453237373544027323918826898781566353727091177232101624598580239994973038912038329546580588079767167651895943550308947180715638564325717561793234218513207167653182319270455562935459992553249433401071963526098078628986437308890133845944591656450159537485709628432052132725117434674074063915340351489772859034593173997670670302856871372969156977207058807985730032820727431391357959732607595988516064241005752724214640614570412163400453896454236261837073650598458053038871507738985949768200045994080490463683606982451876795177007360679367963360631596196216665432615300912190494267441845189611578868452908778094989327991657402208379275569545"),
        ("7dp3", &["1p0"], &["1p0"], "197007111401704699388887935224332312531693798532384578995280299138506385078244119347497807656302688993096381798752022693598298173054461289923262783660152825232320535169584566756192271567602788071422466826314006855168508653497941660316045367817938092905299728580132869945856470286534375900456564355589156220422320260518826112288638358372248724725214506150418881937494100871264232248436315760560377439930623959705844189509050047074217568.226757808330810207066881891196853644591820658492943388594373441606683399590492828162770613598773090497956651224670222796547028060016974015433216920112279419476911933498024014771208957692397594254436621593942610178129942185855427185229801528630341105804209568586616823953605342858090073518818427307513671712518312938822368831025594914114667498754443872668606582490770720339578911220032562819555103422010728982107295731574962192206277209720805104756889364954963599063"),
        ("7dp3", &["1p0"], &["1p-1"], "11042231289744830710957907154773505059804913259935022241423924900504912424121110679841414027258305814701015050577721010426139561662902111667860267814103891185057153990576071771123582340863628707270693268402753344571715171971513056587564793988284746328706552006385523141546298334945402476780981987219463540776889235074487373410901932943706090493377559767413621425129784676388404397083565284472006583630269680975488454054539080408376669624.8900855163843993785004331907664573566305248030291392198754969474358356787128867788906085682784294014212259266479363426965799900653942564161140203818764087251869261885181885802059242435309756641499247680408329892805947627334784268831617987209112419424221328339617389328182851904454411993755783248673297859753923456784973593623876533850497756857210582040839788066492801590884589008045521571640765304093566646504846628079781731497512999881898946389821115442569609214"),
        ("7dp3", &["1p0"], &["3p-1"], "5521115644872415355478953577386752529902456629967511120711962450252456212060555339920707013629152907350507525288860505213069780831451055833930133907051945592528576995288035885561791170431814353635346634201376672285857585985756528293782396994142373164353276003192761570773149167472701238390490993609731770388444617537243686705450966471853045246688779883706810712564892338194202198541782642236003291815134840487744227027269540204188334.81194504275819219968925021659538322867831526240151456960993774847371791783935644338944530428413921470071061296332396817134828999503269712820805701019093820436259346309425909429010296212176548783207496238402041649464029738136673921344158089936045562097121106641698086946640914259522272059968778916243366489298769617283924867968119382669252488784286052910204198940332464007954422945040227607858203826520467833232524233140398908657487564999409494731949105577212848046071"),
        ("7dp3", &["1p0"], &["5p1"], "71489940585450601314239653934205709571501045611431716025807314951381197017193226028820004442319119781814815027131133995052950441038002912887353598934596257220304475802338847584487051506451699735357784761932826807603548420181333069695486543073773375153475165507158615845952395937577594326757678073356193009266851576137071619627301087486121617228285839991864003877477859324164364598312570263192149765402024822498056.7304165370419016302679152458734870844079404301006775142505494798805539192968531262344902332720433980374837062002627227790798985095964083304484109855424189595307204057559703439557397818024277629544801763065674172390030499612440098943814397934204032174169761903787093781804742315682487115130762843068163437258785088309013505611910386453992398612014025678824339305419480125925139839566542508789968263953035254163959601559297792533330270954750739222803078138738634857564141800127548571908279"),
        ("7dp3", &["1p-1"], &["1p0"], "3515731088517800699796201397105864408023954354635858567003961941776273981191737835374593016913564143831906635982235023918617724812196767576200249885395171652444921234443858317973084920155925505886514243351775052436192435956427731436841106171311422528788276860502911065049469237245728831848784380592000615254091274666582455062744725174592882837975978136666736232834700606433888310111850388067293887570773930981269128949666573869990813.87139581883007252713664486960384891392359618794838347669066367368523199976212286162566618605660525523873574031723768441242567878891344346623227121974409936515651116613934431382448905351777621569133066494074600052865679269598732245258194290034836728072427521012867485740710647727337006009946910800214408866841983917844696042557511936159688638456336186880194333495815199563743980518309384601876843538590993564372481813036201336493760261589299321679775718199968471487081"),
        ("7dp3", &["1p-1"], &["1p-1"], "197007111401704699388887935224332312531693798532384578995280299138506385078244119347497807656302688993096381798752022693598298173054461289923262783660152825232320535169584566756192271567602788071422466826314006855168508653497941660316045367817938092905299728580132869945856470286534375900456564355589156220422320260518826112288638358372248724725214506150418881937494100871264232248436315760560377439930623959705844189509050047074217568.226757808330810207066881891196853644591820658492943388594373441606683399590492828162770613598773090497956651224670222796547028060016974015433216920112279419476911933498024014771208957692397594254436621593942610178129942185855427185229801528630341105804209568586616823953605342858090073518818427307513671712518312938822368831025594914114667498754443872668606582490770720339578911220032562819555103422010728982107295731574962192206277209720805104756889364954963599063"),
        ("7dp3", &["1p-1"], &["3p-1"], "98552881541713073068980971904668767089051270887771895644398883239255596308047443767135103509237742668166138632942450962103197010686163793136368390768615420720664259985262074968117132027397326827151206157934635686474527205383961591391590356445164060947980322002346147401360078771687579483210572513574954480968453569750446821536399576946044028116216074080233156356612779722684952604375277046351349584559404254388460669134109596271310.9442944546699092951797575480663059395685891924137657472174414478859877731511796147785439458129748904480188413630568254124627015717216464432983770899035025410243871000388168497849054213893376965111653433916250152003974140301308978872710675354369878542703905150160638213353446712540079936270595914886001628347403170133269187790202916660671515108207049418317198801906835091613819519096559359853654822370568945549373024622120332450608499561653863860482165882331787267313805"),
        ("7dp3", &["1p-1"], &["5p1"], "1281575922454432690721311238954662981307181811421987841159192348739138110452845508758694929819637499902137293706224416100671487383076197399978635801844607652285323166743707331320503752450187563581795129829839630913600556582839953334316345442152518723262721657353588730542888368196625616533913623335512047712704670433818200524407991508061255274623377040518125376297685192757882471989148146432781322190078161722949.05090301539284415189637580682867268001886833767703623429508354726338183390902502619615180748176963790685567868028028433000587987850280973931009476334189112177553821052255476068300514025816036364855895210003499311088506925164774836124017280782179976312851823410047032268805474298884531433187026606810900599737111767070867404882067614699464717242665500093242151688382571888925476390387939488919124592816287924589282914539908166998501365457540027436153136356867486694008972453930929804495844"),
        ("7dp3", &["5p1"], &["1p0"], "588303973191761052187257556369844148210213663891632530671214758784591280754827023273124658200041872428849607145592395154125181510868128598573494892308437427829544353334354743145174085870655134321721069465330607695005673431531986116040229347130213799282896390943332427720434352343423399510837034487319019414733386066509002618329495112163887544874362668847543320771317517039194425569648821308872192426988299205059999964279858758391124314988819523674234802997.127098170596862368990282182087197613164372690880729790013277691796527283011476192559925185849611757848171251597120218891044985374510377475619256646772016577709740797228641166430477659761318802626100896720468644909659198122117542009515513089914449780883427193698273823465815916171923311644568129034184385373790610396369465606341908547283788139627941425292747043833373055394999898443869683696770460892861375880490761888385029366285621859957501214"),
        ("7dp3", &["5p1"], &["1p-1"], "33121167009825164143504441127293654748681559202196399945532607508412548089918312167138330524774397732631734180760038414601670846102255694621314745454036089029024263329479961181684472695229200993358329238134710492576693339252660636864248807396889502898169190190867743697139486703664994507118483785366849294071629565755594905499445431582240282906839991780621312949882183374465340466786251738199612193369391241864475559146931639373839955283505647944322230272228.6659893701188504428500984674174159421019206142384654121356773222074620299451409120215764745805226188666432589946154722977499784572977006642006247604550167355472397426704586412273111474113310003333816452148124421449035457606268010319312542160021184632274128033577112883772510762920790085678773929744489385590158586284526633416167256518535935591789494033159256265402953870504552679204273700909326926762904504977853212046558147430260687527101161"),
        ("7dp3", &["5p1"], &["3p-1"], "16414092044058847678712804673884506640085646208075115776943351604612399269806085265827673184857278008070102345597474754238932824339411199723692139566148681956537617273911332182520232841746801054794596415419106825319207774784064389242786032635787351143664402874717170719983731764774765875018179139765750671520253019573131538024222290633695979905263797757039664531787341015921320130499602636333083771402853192480413919340128844248775418522871119694909306086.7965745318543922127742307004935953048360364961929513947445138720539819201901738604559623216694827554258666804877586616785320150092293752471914789300847965411061060716287166718863580866858683695143356281607112153859234860445342270735415882907055316139104978293772581378516431433518086111055303856360817933774948065781861741260715147825357450230373358080827532679659716153319655865539829486006982198344007820235294313463874331295087316951382715643"),
        ("7dp3", &["5p1"], &["5p1"], "197007111401704699388887935224332312531693798532384578995280299138506385078244119347497807656302688993096381798752022693598298173054461289923262783660152825232320535169584566756192271567602788071422466826314006855168508653497941660316045367817938092905299728580132869945856470286534375900456564355589156220422320260518826112288638358372248724725214506150418881937494100871264232248436315760560377439930623959705844189509050047074217568.226757808330810207066881891196853644591820658492943388594373441606683399590492828162770613598773090497956651224670222796547028060016974015433216920112279419476911933498024014771208957692397594254436621593942610178129942185855427185229801528630341105804209568586616823953605342858090073518818427307513671712518312938822368831025594914114667498754443872668606582490770720339578911220032562819555103422010728982107295731574962192206277209720805104756889364954963599063"),
        ("7dp3", &["aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabp-257"], &["1p0"], "735719724297239634996639984313073572051277214434754469878734295292475205796779334551925165811722975091899606433315323350941163028761726249084758471541536343790438707958964419856991770132220635337440352538639906645149861854929486130523964463254752737972991857437228281714548796415353207280468059561372603443164191456383591902934967422767668832844808732262146806782891050251186095521768685670274245644769060790310320428006437029986663.174741017447196735287315029924556121057998936549491522901427843029121924620787640233167403859413306409882141140307053305024400540083505143780173791314688105961485809674026989295225163528060957819414527485629600034889315893171487457567142468948739368217053701770210764025566957111627116487847608816548239183394942733986980914273390653717064254416672847080230127224322992858762710389883550328799156703430486268590498358631628548940959439579303485783869959210827695744384"),
        ("7dp3", &["aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabp-257"], &["1p-1"], "41223265320922977597662944903981814914472701775013361635433891337359793273732363764160904205475631809728025671132263541224445338362590547217914118969111720691257249476947136330348776660218320337544063462497313801124198997633568002896787837692093669137069223420047549374154309355999611365027812420028541076867438673959393078868436924198581293358087786781155767219816798149271048888140650252895156288608680468355483408080540742575724896.8011371583572228298407619731754680898460852640525404815990978005913434273135174108331359156345214049574406535649264375473670073099758377699247940585130902005203033546083202038795627411436183304427614021623588677622002539740699819671532411295461425148851343390726399408110268749595824965147353931207201835814760524925133997884485273318355959876745301011953046734604673162153046377092120300842503026307754539958036465800450429244271876056392338301946083720162490510714"),
        ("7dp3", &["aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabp-257"], &["5p1"], "268594410802596292320525710401797478364948912747053185910711397594378944728564415034085240623147086423248682508965972419536858052633330033013005334892928287348912615761632179073444223031691038820417822900298443618041525357230152081529570614380101492374543633633808657841016554977120992754516289380726232309502494191694898479065937286127689096094626545248934207479159027859426012343825229222907530011158561799843.491499790732270602922804723134734847632692374367342594379481243615573890544530166963434729827736822428139958935938301168103927006549100753623995116455556528078271571476279354709355244130426007444060541262326527641366927746341513023961433135378310475477997320508391243907032684739835566769318008327248058418169726447052236512407572040859501198775217348629134443632916723423413283193877086649455826592929833960050305402582223850645366552245771551562941710337797893176290120600541067053704303"),
        ("7dp3", &["aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabp-129"], &["3p-1"], "20625399014087849892084975534564615873995007924331971217245816849250962797668007618657177442053390517493246572797146154508492260628763503950857686517826949114504817321926961154331250349291907666941610798745004821096177796907291049112413388470263218962859110433747438671216190296960787554037516979974080520004080407002049730228073912349261317357569121445359032224284268948914450876241147398033881848901255270522487556762363990456876.4280430775624485826079055682896528072190609136169982077777394381143422313810288192787077021807615884889478620498149437219521926595837029821949455087541388352676107964394527710619450638051675417275528721823468404325750162268893270262093899126078821429455910539356137188063775717412605804693627018250223994762246760432046978623709917305642075468646129582057372017736865286654967078400005256567403078658068286074338223875979593832484796818536675625551679784935579705005190"),
        ("7dp3", &["aaaaaaaaaaaaaaabp-65"], &["1p0"], "735719724297239635063368005744428552601964344766313557529686575123613183042611730818150157819292300261731195920720028924132337514934291391118359166674574248467417614418399987270671096468222223813039007574194942207542874848317839692823716500901340463128486189760209143581402660815571004072279906434427379985547072123178910936745156013449265706807179856510738918752424002724722878927175201118497402165564306831690194059357551891653046.439303095345938064624459753615104263919916053423821151138605676150647365227601917524163928988477654825823867512754287760716879909253029966948398458274890824076861197148828384105635850599498736842861453041661022922181369764193017669315301213044800020927885270872572914468877588675325417013194271102365268989792754115854538387606982903972423140913198374501897409412121898412725307018363553372453117955172893979764781888440601321044010210614040182472347545803862639608863"),
        ("7dp3", &["aaaaaaaaaaaaaaabp-65"], &["1p-1"], "41223265320922977601401983237104336928351494697253902892313237201201028413283976422845390448958780897032905033668786982444029506987852221460064760819689444935006773377968697719378649905599571145764831445035664816425360439666523973918305426017161425287103140480100671945718067094718246325516435825558188815805529728836169512099879209010007917797059798098509817623340346522665807990042097089502012331333026450106966721565109207815068282.5594742022036571935612888522971664872402276304643026164626519220977458953470990089003166385532392801178746970649619427816481122454728928479567299857395959536149655361904294132149392251993724742423011525364186394735155359051053137938833285243103160945884507683044203296268866772605434134568772152758230795676371941570469956691124554480758884808689634688398279543365424369773082947725670814288475752481482473853918563331739807060300242093871530508847219556686714206412"),
        ("7dp3", &["aaaaaaaaaaaaaaabp-65"], &["3p-1"], "20625399014087849893955556520391172006647947979083515426162391228826600788483183346622125838469154455575399672655254537925866793822867053049127582528704592398181746131074951220679697653080371491748960726479816367825529282019514405604366606119399986629279458178472484208180122254996150583317188221236776559426275620905169858914955054011596941730876798850470371469122483066718958827722135843144101070143032233920319622342864462827561.3785825867226538682536016664904140320042219304650785478646281648388783490123629891851613854831419625651117128790733431624522758982974001254685731986631093616899538453527126206415488876166755545241618815530208033701934432279209735228612038944326550442768405571735235398803806264214415804947488492435635582465034152828923892341160946671910915851200155207870439955532781084781368404804173155268959653590490145152048378278402635249206274601051552953904287797823892818137997"),
        ("7dp3", &["aaaaaaaaaaaaaaabp-65"], &["5p1"], "268594410802596292344864586604019609794645749372686967696924733160175763900645039434004737647745940196986559148640310150310063202356942126990015223437009480348641183533443678618284475749590197391879608555867207840332127549123688299091906624562243086042873301777497797076696013988920691247190657381158049273910827457768520406718571865110439497220352474299546031521539139086459608931675719594930621952548920774466.135446403323989799245414081258847036317368950051285384026695406718906312006090180687972839914455198479199635650500806911110831317891454220289255501839456611593462841542735562330008659598910093702132003754669702642452380179360899599453003860028428087139992633532434244412230791948333647803727048358487773140973621625611940835286192170384599912961655133739988766334155476454542395018836138080149088280798495966134487066886369284586713393593631713303495819700118146981275248797630437872841132"),
    ];
    for p in [64usize, 128, 256] {
        for &(z, aa, bb, want) in refs {
            let zz = parse_real(z).expect("z");
            let b = parse_real(bb[0]).expect("b");
            let g = if aa.is_empty() {
                zz.hypergeom_0f1(&b, p, RoundingMode::ToEven, &mut cc)
            } else {
                zz.hypergeom_1f1(
                    &parse_real(aa[0]).expect("a"),
                    &b,
                    p,
                    RoundingMode::ToEven,
                    &mut cc,
                )
            };
            assert!(
                quad_near(&g, want, p, 8, &mut cc),
                "z={z} p={p} a={aa:?} b={bb:?}: {}",
                fmt_real(&g)
            );
        }
    }
}

/// ZF-51: integrate_adaptive_gk reported Ok when two successive Wynn ε values agreed,
/// although the estimate missed the requested tolerance by up to 10^100. Harper's
/// rows: an Ok result must now lie within its own abs_err_est of the truth, and that
/// bound must meet the requested tolerance. Otherwise the call must be NotConverged.
/// References: FLINT arb at 1400 bits.
#[test]
fn zf_51_gk_adapt_success_meets_tolerance() {
    use zenith_float::{integrate_adaptive_gk, QuadratureError};
    let mut cc = Consts::new().expect("Consts");
    let e1 = "1.71828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193200305992181741359662904357290033429526059563073813232862794349076323382988075319525101901157383418793070215408914993488416750924476146066808226480016847741185374234544243710753907774499206955170276183860626133138458300";
    let s50 = "-0.00524749707407857571828787293825245091773331929942477045887871184828917905179901503577847382570704739317198027405306417222213647304310017250348792583103888631908764942515594820929406369254832813980963530228047177295338223389756335547195140623093431833416637568940168458545776033945301316421969951103154288667798249064536349136985842";
    let rows: &[(u32, usize, i64, usize, &str)] = &[
        (3, 256, -128, 8, "0.6666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666667"),
        (3, 256, -128, 32, ""),
        (8, 256, -128, 32, ""),
        (0, 1024, -512, 32, ""),
        (3, 1024, -512, 8, ""),
        (3, 1024, -512, 32, ""),
        (8, 1024, -512, 32, ""),
        (0, 64, -50, 8, ""),
    ];
    let mut ok = 0;
    for &(id, p, tl, maxsub, _) in rows {
        let want = match id {
            0 => e1,
            8 => s50,
            _ => rows[0].4,
        };
        let rm = RoundingMode::ToEven;
        let a = ExactNum::new(p);
        let b = ExactNum::from_u8(1, p);
        let tol = ExactNum::from_u8(1, p).ldexp(tl, p, rm);
        match integrate_adaptive_gk(
            |x, q, r, c| harper_integrand(id, x, q, r, c),
            &a,
            &b,
            &tol,
            maxsub,
            p,
            rm,
            &mut cc,
        ) {
            Ok(g) => {
                ok += 1;
                let w = ExactNum::parse(want, Radix::Dec, 1152, rm, &mut cc);
                let d = ExactNum::sub_ref(&g.estimate_kronrod, &w, 1152, rm).abs();
                assert!(
                    matches!(d.compare(&g.abs_err_est), Some(o) if o <= 0),
                    "f{id} p={p} n={maxsub}: error above abs_err_est"
                );
                let lim = tol.clone().mul(&w.abs(), 1152, RoundingMode::Up);
                let floor = ExactNum::from_u8(1, p).ldexp(1 - p as i64, p, rm);
                assert!(
                    matches!(d.compare(&lim), Some(o) if o <= 0)
                        || matches!(d.compare(&floor), Some(o) if o <= 0),
                    "f{id} p={p} n={maxsub}: Ok but error misses tolerance 2^{tl}"
                );
            }
            Err(QuadratureError::NotConverged) => {}
            Err(e) => panic!("f{id} p={p}: {e:?}"),
        }
    }
    assert!(ok >= 1, "smooth p=64 control must still converge");
    // The _estimate form keeps the extrapolated value: within its own bound, flagged unmet.
    use zenith_float::integrate_adaptive_gk_estimate;
    let p = 256;
    let rm = RoundingMode::ToEven;
    let tol = ExactNum::from_u8(1, p).ldexp(-128, p, rm);
    let (g, met) = integrate_adaptive_gk_estimate(
        |x, q, r, c| harper_integrand(3, x, q, r, c),
        &ExactNum::new(p),
        &ExactNum::from_u8(1, p),
        &tol,
        8,
        p,
        rm,
        &mut cc,
    )
    .expect("estimate");
    assert!(!met, "f3 8 panels cannot meet 2^-128");
    let w = ExactNum::parse(rows[0].4, Radix::Dec, 1152, rm, &mut cc);
    let d = ExactNum::sub_ref(&g.estimate_kronrod, &w, 1152, rm).abs();
    assert!(
        matches!(d.compare(&g.abs_err_est), Some(o) if o <= 0),
        "estimate outside its bound"
    );
}

/// ZF-11 (Jeff's 11:54 ruling): `1 + MAX` has an exact value too wide to store (~2^62 bits).
/// `add_full_prec` / `sub_full_prec` refuse with `MemoryAllocation` (no substituted value);
/// `two_sum` still returns `hi` = the correctly rounded sum (MAX, or +Inf in Up/FromZero, as
/// Harper's 2x-checked refs) and `lo` = `a + b − hi` rounded to `p` in `rm`, inexact-flagged
/// when that rounding drops bits. Exact values: lo = 1 for 1 + MAX; for MAX + (1 + 2^-200)
/// at p = 64, lo = 1 + 2^-200 rounded (1 + 2^-63 Up, else 1); for -MAX + 1 in Up,
/// hi = -(MAX − ulp) and lo = 1 − ulp(MAX) rounded up.
#[test]
fn zf_11_two_sum_wide_exact_sum() {
    for p in [64usize, 128, 256, 1024, 2048, 4096, 8192] {
        let one = ExactNum::from_u8(1, p);
        let max = ExactNum::max_value(p);
        let s = one.add_full_prec(&max);
        assert_eq!(
            s.err(),
            Some(zenith_float::Error::MemoryAllocation),
            "add_full_prec p={p}"
        );
        let s = one.sub_full_prec(&max);
        assert_eq!(
            s.err(),
            Some(zenith_float::Error::MemoryAllocation),
            "sub_full_prec p={p}"
        );
        for rm in MODES {
            let (hi, lo) = one.clone().two_sum(&max, p, rm);
            let up = matches!(rm, RoundingMode::Up | RoundingMode::FromZero);
            if up {
                assert!(hi.is_inf() && hi.is_positive(), "hi p={p} {rm:?}: {hi}");
                assert!(lo.is_zero(), "lo p={p} {rm:?}: {lo}");
            } else {
                assert_eq!(hi.compare(&max), Some(0), "hi p={p} {rm:?}: {hi}");
                assert_eq!(lo.compare(&one), Some(0), "lo p={p} {rm:?}: {lo}");
                assert!(!lo.inexact(), "lo exact p={p} {rm:?}");
            }
        }
    }
    let p = 64usize;
    let max = ExactNum::max_value(p);
    let one = ExactNum::from_u8(1, p);
    let tiny = ExactNum::from_u8(1, 256).ldexp(-200, 256, RoundingMode::None);
    let b = ExactNum::add_ref(&ExactNum::from_u8(1, 256), &tiny, 256, RoundingMode::None);
    for rm in [RoundingMode::Down, RoundingMode::ToZero, RoundingMode::ToEven] {
        let (hi, lo) = max.clone().two_sum(&b, p, rm);
        assert_eq!(hi.compare(&max), Some(0), "hi {rm:?}");
        assert_eq!(lo.compare(&one), Some(0), "lo {rm:?}: {lo}");
        assert!(lo.inexact(), "lo inexact {rm:?}");
    }
    let nmax = -max.clone();
    let (hi, lo) = nmax.clone().two_sum(&one, p, RoundingMode::Up);
    let want_hi = nmax.next_after(&INF_POS, p, RoundingMode::None);
    assert_eq!(hi.compare(&want_hi), Some(0), "hi -MAX+1 Up: {hi}");
    let ulp = ExactNum::sub_ref(&want_hi, &nmax, p, RoundingMode::None);
    let want_lo = ExactNum::sub_ref(&one, &ulp, p, RoundingMode::Up);
    assert_eq!(
        lo.compare(&want_lo),
        Some(0),
        "lo -MAX+1 Up: {lo} want {want_lo}"
    );
    assert!(lo.inexact(), "lo inexact -MAX+1 Up");
}

/// ZF-11 bessel_j rows: `J_1` at the minimum exponent (`minpos`, the other subnormals, the
/// least normal) returned NaN(MemoryAllocation) because the series widened by the result's
/// ~2^62-bit exponent; `J_0` there skipped the tiny-argument branch. Exact semantics (C99
/// correctly rounded, the library's subnormal grid): `J_1(x) = x/2 − δ`, `J_0(x) = 1 − δ`
/// with `δ` far below an ulp. Harper's mpmath refs are not rounded to the subnormal grid,
/// so the expected values are derived from exact `ExactNum` arithmetic instead.
#[test]
fn zf_11_bessel_j_at_min_exponent() {
    let mut cc = Consts::new().expect("Consts");
    for p in [64usize, 128, 256, 1024, 2048] {
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);
        let mp = ExactNum::min_positive(p);
        let mn = ExactNum::min_positive_normal(p);
        let two_mp = ExactNum::add_ref(&mp, &mp, p, RoundingMode::None);
        let half_mn = mn.clone().ldexp(-1, p, RoundingMode::None);
        assert!(
            !half_mn.is_zero() && !half_mn.inexact(),
            "minnorm/2 on grid p={p}"
        );
        let below = |v: &ExactNum| v.next_after(&zero, p, RoundingMode::None);
        for rm in MODES {
            let away = matches!(rm, RoundingMode::Up | RoundingMode::FromZero);
            // `None` truncates the excess bits, as the library's tiny-argument branches do.
            let toward = matches!(
                rm,
                RoundingMode::Down | RoundingMode::ToZero | RoundingMode::None
            );
            // minpos/2 is the midpoint between 0 and minpos.
            let want = if away { mp.clone() } else { zero.clone() };
            let j = mp.bessel_j(1, p, rm, &mut cc);
            assert_eq!(
                j.compare(&want),
                Some(0),
                "J1(minpos) p={p} {rm:?}: {j} {:?}",
                j.err()
            );
            assert!(
                j.is_zero() == want.is_zero() && !j.is_negative(),
                "J1(minpos) sign p={p} {rm:?}"
            );
            // -minpos: odd function; Down is away from zero for a negative result.
            let away_n = matches!(rm, RoundingMode::Down | RoundingMode::FromZero);
            let want = if away_n { -mp.clone() } else { -zero.clone() };
            let j = (-mp.clone()).bessel_j(1, p, rm, &mut cc);
            assert_eq!(j.compare(&want), Some(0), "J1(-minpos) p={p} {rm:?}: {j}");
            assert!(j.is_negative(), "J1(-minpos) sign p={p} {rm:?}: {j}");
            // 2·minpos and minnorm: x/2 is on the grid.
            for (x, h) in [(&two_mp, &mp), (&mn, &half_mn)] {
                let want = if toward { below(h) } else { h.clone() };
                let j = x.bessel_j(1, p, rm, &mut cc);
                assert_eq!(
                    j.compare(&want),
                    Some(0),
                    "J1({x}) p={p} {rm:?}: {j} want {want}"
                );
            }
            // J0 = 1 − δ.
            let want = if toward { below(&one) } else { one.clone() };
            for x in [&mp, &mn] {
                let j = x.bessel_j(0, p, rm, &mut cc);
                assert_eq!(j.compare(&want), Some(0), "J0({x}) p={p} {rm:?}: {j}");
            }
        }
    }
}

/// ZF-11 riemann_zeta row: ζ(s) for s = -0xc4eea7fd3d243a17·2^-39 ≈ -2.58·10^7 reflects to
/// ζ(1-s), whose Euler–Maclaurin cut grew as 2(1-s): ~5·10^7 trial-divided terms, hours,
/// then NaN(MemoryAllocation). Reference: FLINT arb zeta at 900 bits (mpmath agrees to 60
/// digits), rounded once per mode.
#[test]
fn zf_11_riemann_zeta_large_negative() {
    let mut cc = Consts::new().expect("Consts");
    let want = "1.621596734986527778777610221795721975338643829734992759188335692762124964237794899915069460092625623907760252390890221657635590894413672477676035983145507405973e+159503325";
    let s = parse_real("-c4eea7fd3d243a17p-39").expect("s");
    for p in [64usize, 128] {
        for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
            let z = s.riemann_zeta(p, rm, &mut cc);
            let w = cr_ref(want, p, rm, &mut cc);
            assert_eq!(
                z.compare(&w),
                Some(0),
                "p={p} {rm:?}: {z} {:?} want {w}",
                z.err()
            );
        }
    }
}

/// ZF-53 (Jeff's 11:54 ruling): primality is proven, like FLINT's `fmpz_is_prime`, and
/// `pollard_rho` has no iteration cap. Before, `pollard_rho(M89)` spent 2^20 · 32 steps
/// (144 s in release) and its `None` meant only "cap hit"; a composite whose factors are
/// past the cap also came back `None`. Primality and factors from FLINT (`fmpz.is_prime`,
/// `fmpz.factor`). 3317044064679887385961981 is the least strong pseudoprime to the first
/// 13 prime bases (Sorenson–Webster), so it exercises the Pocklington path.
#[test]
fn zf_53_pollard_rho_proven_primality_uncapped() {
    use zenith_float::{is_prime, pollard_rho};
    let int = ExactInt::from_u128;
    let primes: [u128; 6] = [
        (1u128 << 89) - 1,
        (1u128 << 127) - 1,
        18446744073709551629,
        4835703278458516698824713,
        1_000_000_007,
        2,
    ];
    for &n in &primes {
        assert!(is_prime(&int(n)), "{n} is prime");
        if n > 2 {
            assert!(pollard_rho(&int(n)).is_none(), "rho({n}) on a prime");
        }
    }
    let composites: [(u128, &[u128]); 4] = [
        (3317044064679887385961981, &[1287836182261, 2575672364521]),
        ((1u128 << 67) - 1, &[193707721, 761838257287]),
        (561, &[3, 11, 17]),
        (1_000_000_007u128 * 1_000_000_007, &[1_000_000_007]),
    ];
    for &(n, fs) in &composites {
        assert!(!is_prime(&int(n)), "{n} is composite");
        let f = pollard_rho(&int(n)).unwrap_or_else(|| panic!("rho({n}) found no factor"));
        let (q, r) = int(n).div_rem(&f).expect("div");
        assert!(
            r.is_zero() && f != int(1) && f != int(n),
            "rho({n}) = {f:?}"
        );
        // The factor is a product of the listed primes.
        let mut g = f.clone();
        for &p in fs {
            while let Some((d, r)) = g.div_rem(&int(p)) {
                if !r.is_zero() {
                    break;
                }
                g = d;
            }
        }
        assert!(g == int(1), "rho({n}) = {f:?}, cofactor {q:?}");
    }
    for n in [0u128, 1, 4, 9, 15] {
        assert!(!is_prime(&int(n)), "{n}");
    }
}

/// ZF-51 (standard reading): the reported bound includes the final rounding ulp, and `Ok`
/// requires that whole bound to meet the tolerance. ∫_0^1 1 dx at p = 64 has Σ|K−G| = 0,
/// but the result's rounding ulp 2^-63 exceeds a requested 2^-70, so no p = 64 result can
/// carry that guarantee. At 2^-60 it is met.
#[test]
fn zf_51_reported_bound_includes_final_rounding() {
    use zenith_float::{integrate_adaptive_gk, integrate_adaptive_gk_estimate, QuadratureError};
    let mut cc = Consts::new().expect("Consts");
    let p = 64;
    let rm = RoundingMode::ToEven;
    let (a, b) = (ExactNum::new(p), ExactNum::from_u8(1, p));
    let one = |_x: &ExactNum, q: usize, _r: RoundingMode, _c: &mut Consts| ExactNum::from_u8(1, q);
    let tight = ExactNum::from_u8(1, p).ldexp(-70, p, rm);
    let r = integrate_adaptive_gk(one, &a, &b, &tight, 4, p, rm, &mut cc);
    assert!(
        matches!(r, Err(QuadratureError::NotConverged)),
        "2^-70 at p=64: {r:?}"
    );
    let (g, met) =
        integrate_adaptive_gk_estimate(one, &a, &b, &tight, 4, p, rm, &mut cc).expect("estimate");
    assert!(
        !met && g.estimate_kronrod.compare(&b) == Some(0),
        "estimate {g:?}"
    );
    assert!(
        g.abs_err_est.compare(&tight) == Some(1),
        "bound above 2^-70"
    );
    let loose = ExactNum::from_u8(1, p).ldexp(-60, p, rm);
    let g = integrate_adaptive_gk(one, &a, &b, &loose, 4, p, rm, &mut cc).expect("2^-60 met");
    assert!(
        g.abs_err_est.compare(&loose) != Some(1),
        "bound meets 2^-60"
    );
}

/// ZF-11 follow-up: ζ(s) by reflection for large |s|. The factors were formed at the Ziv
/// precision from rounded arguments (π·s, 1 − s), so an error ε grew by |s| ln|s|, which is
/// silent for |s| ≥ 2^100 when the margin is 128 guard bits. Now the arguments are exact
/// (s mod 4, exact 1 − s), the precision gains log2|s| + log2 log2|s| bits, large |s|
/// uses the log-domain product, and |s| ≥ 2^64 overflows with the sign of sin(πs/2).
/// References: FLINT arb zeta at 1200 bits, rounded once per mode.
#[test]
fn zf_11_riemann_zeta_reflection_large_s() {
    let mut cc = Consts::new().expect("Consts");
    let rows: &[(&str, &str)] = &[
        ("-10001p-1", "-2.5807773263234762242686736608678069832777872362024525369031599519101624469532782025158025749453974044076352582365374318350339094757241552377e+107580"),
        ("-100000001p-2", "-7.3149691613544282365887632717126612596241972256522213841379695897846438755602988516562540871499607097550476731729235225593151515339289066546e+8373495647"),
        ("-40000000003p-2", "-3.3449479113386980706833361499789638222658442499664767958603353568195491386477465027405858952472906904914675684466536870028970604744879109652e+11884319342642"),
        ("-20000000000003p-3", "-5.3110969919362335732524370696865304216932764417203059129646495388257533981727632915420429572230056423795400428108870223719625509067451511160e+15558839447600412"),
    ];
    for p in [64usize, 128, 256] {
        for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
            for &(s, want) in rows {
                let s = parse_real(s).expect("s");
                let z = s.riemann_zeta(p, rm, &mut cc);
                let w = cr_ref(want, p, rm, &mut cc);
                assert_eq!(
                    z.compare(&w),
                    Some(0),
                    "zeta({s}) p={p} {rm:?}: {z} {:?} want {w}",
                    z.err()
                );
            }
        }
    }
    // Overflow: s = -(2^k + 1/2) has sin(πs/2) < 0, s = -(2^k + 5/2) has sin > 0.
    for (s, neg) in [
        ("-8000000000000001p-1", true),
        ("-8000000000000005p-1", false),
        ("-20000000000000000000000001p-1", true),
    ] {
        let s = parse_real(s).expect("s");
        for p in [128usize, 256] {
            for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
                let z = s.riemann_zeta(p, rm, &mut cc);
                let toward = match rm {
                    RoundingMode::ToZero => true,
                    RoundingMode::Up => neg,
                    RoundingMode::Down => !neg,
                    _ => false,
                };
                let m = ExactNum::max_value(p);
                let want = match (toward, neg) {
                    (true, true) => -m,
                    (true, false) => m,
                    (false, true) => INF_NEG,
                    (false, false) => INF_POS,
                };
                assert_eq!(
                    z.compare(&want),
                    Some(0),
                    "zeta({s}) p={p} {rm:?}: {z} {:?}",
                    z.err()
                );
                assert_eq!(z.is_negative(), neg, "sign zeta({s}) p={p} {rm:?}");
            }
        }
    }
}

/// ZF-38 / ZF-41 / Jeff's 11:54 2F1 ruling: real ₂F₁ at z = ±1 returned NaN ("for now")
/// although the value is defined. z = −1 by Pfaff, z = 1 by Gauss for s = c−a−b > 0, and
/// ±Inf when the series diverges (mpmath's ±inf). References: FLINT acb hypgeom_2f1 at
/// 1400 bits; the exact-zero Gauss row (c−a = −1) from mpmath, where FLINT gives nan.
#[test]
fn zf_38_hypergeom_2f1_unit_argument() {
    let mut cc = Consts::new().expect("Consts");
    let rows: &[(&str, &str, &str, &str, &str)] = &[
        ("-1p0", "1p-1", "1p0", "3p-1", "0.78539816339744830961566084581987572104929234984377645524373614807695410157155224965700870633552926699553702162832057666177346115238764555793133985203212027936257102567548463027638991115573723873259549110720274391648336153211891205844669579131780047728641214173086508715261358166205334840181506228531843114675165157889704372038023024070731352292884109197314759000283263263720511663034603673798537790235826431759143989799"),
        ("-1p0", "1p0", "1p0", "1p1", "0.69314718055994530941723212145817656807550013436025525412068000949339362196969471560586332699641868754200148102057068573368552023575813055703267075163507596193072757082837143519030703862389167347112335011536449795523912047517268157493206515552473413952588295045300709532636664265410423915781495204374043038550080194417064167151864471283996817178454695702627163106454615025720740248163777338963855069526066834113727387372"),
        ("-1p0", "3p-1", "5p-1", "1p-1", "-0.265165042944955321650316635789318389731813476633177763720627450873262339711645069784447662686432794887815096168296055692173406760472013227370772780706237968389979249861407298641712490718966146636357150645542494676119426475236911627761597015995166724739291061730715807883935272032455488180249109638998066360011852437214155557616365360941074072764828178954992824344170353813287223189628945238139093533842572543744506063403"),
        ("1p0", "1p-1", "1p-2", "1p0", "1.669253683348146372562859465598093617987986026980694004899654740207363985419052823739382320702550648258135700704824020173449578015269500785318121053485425120641371997947503042923148214381335429524622489233288103742767935690280567738400902762671781517106040996960105608938632871618737957913781378348241474582188133116350374050954802149356759433052713506663813090905035385922061120108136989325097928437228589441504297307316322152404292075040697550719241713542336584400650718219179799786139666399393954671767478209918171934805836493913749383297149038213455148157222876285522116036864069085230987509220191417823504213684027160662495078438762780464781819046868177484453836877772598844451852670432789560220441981504298138792861672076515435882389741193479285453116701696643387942673901190232477656303094940457284155515507421203776094464978092370506037791940513796360277576097914718586930807892662643166775710657090449868832835303464046993654067494418829529929929212474758879914103817333965096251924931552992"),
        ("1p0", "5p-2", "-1p-2", "9p-2", "0.7725311155422382993474376446646916859239995098007362106608033521066757403492166576488632756876549571990716129322749176552046915602730307989724457143199940872798489213805402716628433031841806111331109634801322053149773002023416800849925711282795805881879260780917341163567208136143771148686863071780150974223544287209563865156467950819765437433484922160777402980147816535255103175129717618451319704173"),
        ("1p0", "7p-1", "-3p-1", "5p-1", "0"),
        ("1p0", "1p0", "1p0", "1p1", "+Inf"),
        ("1p0", "1p-1", "1p0", "1p0", "+Inf"),
        ("1p0", "-1p-1", "1p0", "1p-2", "-Inf"),
    ];
    for p in [64usize, 128, 256, 1024] {
        for &(z, a0, a1, c, want) in rows {
            let zz = parse_real(z).expect("z");
            let a = [parse_real(a0).expect("a0"), parse_real(a1).expect("a1")];
            let b = [parse_real(c).expect("c")];
            let g = zz.hypergeom_pfq(&a, &b, p, RoundingMode::ToEven, &mut cc);
            match want {
                "+Inf" => assert!(
                    g.is_inf() && g.is_positive(),
                    "{z} {a0} {a1} {c} p={p}: {}",
                    fmt_real(&g)
                ),
                "-Inf" => assert!(
                    g.is_inf() && g.is_negative(),
                    "{z} {a0} {a1} {c} p={p}: {}",
                    fmt_real(&g)
                ),
                "0" => assert!(g.is_zero(), "{z} {a0} {a1} {c} p={p}: {}", fmt_real(&g)),
                w => assert!(
                    quad_near(&g, w, p, 8, &mut cc),
                    "{z} {a0} {a1} {c} p={p}: {}",
                    fmt_real(&g)
                ),
            }
        }
    }
}

/// ZF-41: G^{1,2}_{2,2}(z | 1,1; 1,0) = ln(1+z) was ~10 bits at z = 1 (the residue series
/// was a perturbed limit around the refused ₂F₁ at |z| = 1) and finite garbage at z = −1.
/// With ₂F₁ at ±1 computed (5a29121), z = 1 is ln 2 and z = −1 is the pole −Inf (a
/// divergent ₂F₁(1,1;2;1) residue term is now a signed infinity, not a limit). Fox H with
/// unit scales is the same function. ln 2 from FLINT arb at 3000 bits.
#[test]
fn zf_41_meijer_fox_ln1pz_at_unit() {
    let mut cc = Consts::new().expect("Consts");
    let ln2 = "0.6931471805599453094172321214581765680755001343602552541206800094933936219696947156058633269964186875420014810205706857336855202357581305570326707516350759619307275708283714351903070386238916734711233501153644979552391204751726815749320651555247341395258829504530070953263666426541042391578149520437404303855008019441706416715186447128399681717845469570262716310645461502572074024816377733896385506952606683411372738737229289564935470257626520988596932019650585547647033067936544325476327449512504060694381471046899465062201677204245245296126879465461931651746813926725041038025462596568691441928716082938031727143677826548775664850856740776484514644399404614226031930967354025744460703080960850474866385231381816767514386674766478908814371419854942315199735488037516586127535291661000710535582498794147295092931138971559982056543928717000721808576103";
    for p in [64usize, 128, 256, 1024, 2048] {
        let one = ExactNum::from_u8(1, p);
        let zero = ExactNum::new(p);
        for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
            for z in [one.clone(), -one.clone()] {
                let g = z.meijer_g(
                    &[one.clone(), one.clone()],
                    &[],
                    std::slice::from_ref(&one),
                    std::slice::from_ref(&zero),
                    p,
                    rm,
                    &mut cc,
                );
                let h = z.fox_h(
                    &[(one.clone(), one.clone()), (one.clone(), one.clone())],
                    &[],
                    &[(one.clone(), one.clone())],
                    &[(zero.clone(), one.clone())],
                    p,
                    rm,
                    &mut cc,
                );
                for (name, v) in [("G", &g), ("H", &h)] {
                    if z.is_negative() {
                        assert!(
                            v.is_inf() && v.is_negative(),
                            "{name}(-1) p={p} {rm:?}: {v}"
                        );
                    } else {
                        let w = cr_ref(ln2, p, rm, &mut cc);
                        assert_eq!(
                            v.compare(&w),
                            Some(0),
                            "{name}(1) p={p} {rm:?}: {v} want {w}"
                        );
                    }
                }
            }
        }
    }
}

/// ZF-44 / ZF-67: complex 2F1 past the unit disc where the 1/z transform is degenerate (integer a − b), where a
/// Γ in a denominator has a pole (c = a or c = b), and on the cut z > 1. At a real z > 1, both
/// signs of a zero imaginary part take the lower side, like Arb (the same sheet as z − 0i).
/// References: FLINT acb_hypgeom_2f1 at exact real z > 1 and at z − i·2^-3000, 1400 bits,
/// cross-checked against mpmath to 85 digits. Parts that FLINT shows as ~1e-904 off the axis
/// are exactly 0 on it (e.g. (1−z)^{−1/2} = +i at z = 2 − 0i) and are expected as 0.
#[test]
fn zf_44_67_c_hypergeom_2f1_outside_unit_disc() {
    let mut cc = Consts::new().expect("Consts");
    let refs: [(&str, &str, &str, &str, &str, &str); 37] = [
        ("1p-1,0", "1p-1,0", "1p0,0", "1p1,0", "0.834626841674073186281429732799046808993993013490347002449827370103681992709526411869691160", "-0.834626841674073186281429732799046808993993013490347002449827370103681992709526411869691160"),
        ("1p-1,0", "1p-1,0", "1p0,0", "1p1,-0", "0.834626841674073186281429732799046808993993013490347002449827370103681992709526411869691160", "-0.834626841674073186281429732799046808993993013490347002449827370103681992709526411869691160"),
        ("1p-1,0", "1p-1,0", "1p0,0", "1p0,1p0", "0.960810086136863131467986827938751397011613243714171883506450459923558444882247604964663070", "0.397980568542750383168216809736942112184308404780933296819710861042718924233724217979965421"),
        ("1p-1,0", "1p-1,0", "3p-1,0", "1p1,0", "1.11072073453959156175397024751517342465365542234392255577134890173910869827486847764383173", "-0.623225240140230513394020080250568002650695312346567252898714776096170004547014180467669073"),
        ("1p-1,0", "1p-1,0", "3p-1,0", "1p1,-0", "1.11072073453959156175397024751517342465365542234392255577134890173910869827486847764383173", "-0.623225240140230513394020080250568002650695312346567252898714776096170004547014180467669073"),
        ("1p-1,0", "1p-1,0", "3p-1,0", "1p0,1p0", "1.02200637560365811103709779776167746474404098608460339057958663185211904353716266841054303", "0.272308224999276444944346988699785356524619518337697181666521127021501580272957140466250113"),
        ("1p-1,0", "1p-1,0", "1p1,1p0", "1p1,0", "1.04377236877632735236564482042483199210638012525100503050471394648506459617136795109733123", "-0.304776271660972671835830552406891883887652910752649433351367835072316485862962505258439927"),
        ("1p-1,0", "1p-1,0", "1p1,1p0", "1p1,-0", "1.04377236877632735236564482042483199210638012525100503050471394648506459617136795109733123", "-0.304776271660972671835830552406891883887652910752649433351367835072316485862962505258439927"),
        ("1p-1,0", "1p-1,0", "1p1,1p0", "1p0,1p0", "1.14869874162767852959863889234171079455340197870563393711730636844440811511480974183405513", "0.170554643540511340364224157388069729495197407572957021666435162638932948490607149827656379"),
        ("1p-1,0", "1p0,0", "1p0,0", "1p1,0", "0", "-1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"),
        ("1p-1,0", "1p0,0", "1p0,0", "1p1,-0", "0", "-1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"),
        ("1p-1,0", "1p0,0", "1p0,0", "1p0,1p0", "0.707106781186547524400844362104849039284835937688474036588339868995366239231053519425193767", "0.707106781186547524400844362104849039284835937688474036588339868995366239231053519425193767"),
        ("1p-1,0", "1p0,0", "1p0,0", "7ffff79c842fa509p-63,431bde82d7b634dbp-82", "776.886987015015888201030570501776566026117178899495215317047486182467141933622487180402317", "321.797126452784635958171535212363696516941843379893194105340944356674923876204882459912078"),
        ("1p-1,0", "-1p-1,0", "1p0,0", "1p1,0", "0.381379881750906594031162990481807896646314618674036457608535799132211346478128096097733073", "0.381379881750906594031162990481807896646314618674036457608535799132211346478128096097733073"),
        ("1p-1,0", "-1p-1,0", "1p0,0", "1p1,-0", "0.381379881750906594031162990481807896646314618674036457608535799132211346478128096097733073", "0.381379881750906594031162990481807896646314618674036457608535799132211346478128096097733073"),
        ("1p-1,0", "-1p-1,0", "1p0,0", "1p0,1p0", "0.817318538373358049311725309802364445263380613423935612833323966743180785420660834874594908", "-0.338544423373199823467714308569407991542338901003064991702865079474996198615359453467183420"),
        ("1p-1,0", "-1p-1,0", "3p-1,0", "1p1,0", "0.555360367269795780876985123757586712326827711171961277885674450869554349137434238821915867", "0.188387379929884743302989959874715998674652343826716373550642611951914997726492909766165463"),
        ("1p-1,0", "-1p-1,0", "3p-1,0", "1p1,-0", "0.555360367269795780876985123757586712326827711171961277885674450869554349137434238821915867", "0.188387379929884743302989959874715998674652343826716373550642611951914997726492909766165463"),
        ("1p-1,0", "-1p-1,0", "3p-1,0", "1p0,1p0", "0.864556578395102817718971079933263252014438461886538713583963250423742641384108093917868399", "-0.217399278093635539728248686702531841380108209675388427460909370986932329479048189479471827"),
        ("1p-1,0", "-1p-1,0", "1p1,1p0", "1p1,0", "0.801873841601499683336523883921601584523067055193576790602683897224535445694255147047908199", "0.186750841232378239189960738778897099048608182293316328247452053361303268276685251341882122"),
        ("1p-1,0", "-1p-1,0", "1p1,1p0", "1p1,-0", "0.801873841601499683336523883921601584523067055193576790602683897224535445694255147047908199", "0.186750841232378239189960738778897099048608182293316328247452053361303268276685251341882122"),
        ("1p-1,0", "-1p-1,0", "1p1,1p0", "1p0,1p0", "0.837842118445094408758976240447433747350760425242786601243337857725828477694224828898596513", "-0.0838106481652454371912254234827097022047987681531621166655938091299576426402452174617410014"),
        ("1p0,0", "1p-1,0", "1p0,0", "1p1,0", "0", "-1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"),
        ("1p0,0", "1p-1,0", "1p0,0", "-5p1,3p0", "0.293531247524081783963613014637468034562682954473426466804824785394911719837502384568263164", "0.0393091421638675936610573316350128101109253557801601641145307943057197383279450089802544610"),
        ("1p0,0", "1p-1,0", "1p0,0", "1p0,1p0", "0.707106781186547524400844362104849039284835937688474036588339868995366239231053519425193767", "0.707106781186547524400844362104849039284835937688474036588339868995366239231053519425193767"),
        ("1p0,0", "1p0,0", "1p0,0", "1p1,0", "-1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000", "-0"),
        ("1p0,0", "1p0,0", "1p0,0", "1p1,-0", "-1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000", "-0"),
        ("1p0,0", "1p0,0", "1p0,0", "1p0,1p0", "0", "1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"),
        ("1p0,0", "1p0,0", "1p0,0", "1p-1,1p-1", "1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000", "1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"),
        ("1p0,0", "1p0,0", "3p-1,0", "1p1,0", "-0.623225240140230513394020080250568002650695312346567252898714776096170004547014180467669073", "-1.11072073453959156175397024751517342465365542234392255577134890173910869827486847764383173"),
        ("1p0,0", "1p0,0", "3p-1,0", "1p1,-0", "-0.623225240140230513394020080250568002650695312346567252898714776096170004547014180467669073", "-1.11072073453959156175397024751517342465365542234392255577134890173910869827486847764383173"),
        ("1p0,0", "1p0,0", "3p-1,0", "1p0,1p0", "0.530116646135371858449400652341695990385764884757438943579284905471062392552807905295240466", "0.915218631075092897572048777678051030781244817158149494917156307249996856128559958038351380"),
        ("1p0,0", "1p0,0", "1p1,1p0", "1p1,0", "0.451904712350074751721003696903098579780122706548984104159078756776980440244041579168493171", "-0.857110025516451727255841125814462926516354572554424162767032236372018540460463274415067880"),
        ("1p0,0", "1p0,0", "1p1,1p0", "1p1,-0", "0.451904712350074751721003696903098579780122706548984104159078756776980440244041579168493171", "-0.857110025516451727255841125814462926516354572554424162767032236372018540460463274415067880"),
        ("1p0,0", "1p0,0", "1p1,1p0", "1p0,1p0", "1.23197561004635516930748692988399069924925306955296508592360443707565757852162304195378142", "1.01293219295213684827213984775517212769105482544314407133129331082312809072693621001532325"),
        ("1p-1,0", "3p-1,0", "3p0,0", "-3p0,1p1", "0.627197658547468993247607607425446180001995307821144066705790867065189853764081241126896145", "0.105837111913788062894585271318001148881668701901892821045542814639791892883392376459135474"),
        ("1p0,0", "3p0,0", "1p1,0", "-5p0,1p-1", "0.0963614744351961950059453032104637336504161712247324613555291319857312722948870392390011891", "0.00917954815695600475624256837098692033293697978596908442330558858501783590963139120095124851"),
    ];
    for p in [64usize, 256] {
        for (a, b, c, z, re, im) in refs {
            let (a, b, c, zz) = (
                parse_complex(a).expect("a"),
                parse_complex(b).expect("b"),
                parse_complex(c).expect("c"),
                parse_complex(z).expect("z"),
            );
            for rm in [RoundingMode::ToEven, RoundingMode::Up, RoundingMode::Down] {
                let w = a.hypergeom_2f1(&b, &c, &zz, p, rm, &mut cc);
                let wr = ExactNum::parse(re, Radix::Dec, p, rm, &mut cc);
                let wi = ExactNum::parse(im, Radix::Dec, p, rm, &mut cc);
                let tag = format!("2f1 z={z} p={p} {rm:?}: {}", fmt_c(&w));
                assert_eq!(w.re().compare(&wr), Some(0), "re {tag}");
                assert_eq!(w.im().compare(&wi), Some(0), "im {tag}");
                // Exact zeros keep the signed zero of the gold (lower-side values
                // at z = 2+0i may have a −0 part).
                for (got, want) in [(w.re(), re), (w.im(), im)] {
                    if want == "0" || want == "-0" {
                        assert_eq!(got.is_negative(), want == "-0", "zero sign {tag}");
                    }
                }
            }
        }
    }
    // z = 1 with real parameters and c − a − b ≤ 0: the sum diverges to +Inf (imaginary part +0).
    for (a, b, c) in [
        ("1p-1,0", "1p-1,0", "1p0,0"),
        ("1p-1,0", "1p0,0", "3p-1,0"),
        ("1p0,0", "1p0,0", "3p-1,0"),
    ] {
        let w = parse_complex(a).expect("a").hypergeom_2f1(
            &parse_complex(b).expect("b"),
            &parse_complex(c).expect("c"),
            &parse_complex("1p0,0").expect("z"),
            64,
            RoundingMode::ToEven,
            &mut cc,
        );
        assert!(
            w.re().is_inf() && w.re().is_positive(),
            "2f1({a},{b},{c};1): {}",
            fmt_c(&w)
        );
        assert!(
            w.im().is_zero() && !w.im().is_negative(),
            "2f1({a},{b},{c};1): {}",
            fmt_c(&w)
        );
    }
}

/// ZF-11 follow-up (large-|s| spot check, 2026-11-15): with a wide-mantissa argument the
/// working precision follows the input width, and `ln_gamma` sent z ≈ 2^62 to Γ(z), which
/// overflows although ln Γ(z) ≈ 2^67.4 is finite. The ζ reflection then returned an
/// ExponentOverflow with the wrong sign: ζ(−(2^62+1)) came back +Inf / +MAX (FLINT: negative).
/// Now large z always takes the Stirling series with enough terms. References: FLINT arb
/// lgamma at 1200 bits; ζ signs from FLINT arb zeta (the value overflows EXPONENT_MAX).
#[test]
fn zf_11_wide_argument_ln_gamma_and_zeta_sign() {
    let mut cc = Consts::new().expect("Consts");
    let lg: &[(&str, &str)] = &[
        ("4611686018427387905", "193576097982213774845.13027697657098014944313953265391747471882753939795149052818771495315937924339432188369712088602790649758523644844215971"),
        ("216172782113783808", "8412332341206235586.0881246231845805591580605550737661624681615854095037107394831826491636613235279472084371792054678902369170253322419789864"),
        ("281474976710663", "9083497175537599.9941806189084721115639498679731370425980111174612984549371271696075125087720783022347650839764140646196196702565384594329106"),
    ];
    for &(x, want) in lg {
        let x = ExactNum::parse(x, Radix::Dec, 4096, RoundingMode::ToEven, &mut cc);
        for p in [64usize, 128, 256] {
            for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
                let g = x.ln_gamma(p, rm, &mut cc);
                let w = cr_ref(want, p, rm, &mut cc);
                assert_eq!(
                    g.compare(&w),
                    Some(0),
                    "ln_gamma({x}) p={p} {rm:?}: {g} {:?} want {w}",
                    g.err()
                );
            }
        }
    }
    // (numerator, halve, negative result)
    for (n, half, neg) in [
        ("-4611686018427387905", false, true),
        ("-9223372036854775809", true, true),
        ("-9223372036854775813", true, false),
    ] {
        let mut s = ExactNum::parse(n, Radix::Dec, 4096, RoundingMode::ToEven, &mut cc);
        if half {
            s = s.ldexp(-1, 4096, RoundingMode::ToEven);
        }
        for p in [64usize, 128] {
            for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
                let z = s.riemann_zeta(p, rm, &mut cc);
                let toward = match rm {
                    RoundingMode::ToZero => true,
                    RoundingMode::Up => neg,
                    RoundingMode::Down => !neg,
                    _ => false,
                };
                let m = ExactNum::max_value(p);
                let want = match (toward, neg) {
                    (true, true) => -m,
                    (true, false) => m,
                    (false, true) => INF_NEG,
                    (false, false) => INF_POS,
                };
                assert_eq!(
                    z.compare(&want),
                    Some(0),
                    "zeta({s}) p={p} {rm:?}: {z} {:?}",
                    z.err()
                );
                assert_eq!(z.is_negative(), neg, "sign zeta({s}) p={p} {rm:?}");
            }
        }
    }
}

/// ZF-08: to_bytes / to_inline_bytes refused finite values with |exponent| ≥ 2^31 (the
/// exponent field was an i32; EXPONENT_MAX is 2^62 − 1). New records in format v1: a
/// wide-exponent heap record (flag 0x12/0x13, 8-byte exponent after the header) encodes every
/// finite value, and a wide-exponent inline record (0x0C/0x0D) holds the top 40 mantissa
/// bits. A full 64-bit mantissa beside a 63-bit exponent cannot fit 16 bytes, so the inline
/// form refuses it with MemoryAllocation, as it does a wide mantissa. Rows: Harper's 5227f77
/// T2-bytes-format-v1 values (exact round trip judged by the raw parts).
#[test]
fn zf_08_bytes_wide_exponent_round_trip() {
    let vals = [
        "1p1099511627776",
        "1p-1099511627776",
        "-1p1099511627776",
        "3p-4611686018427387000",
        "maxv:P",
        "-maxv:P",
        "minnorm:P",
        "-minnorm:P",
        "minpos:P",
    ];
    for p in [64usize, 128, 256, 1024, 2048] {
        let ps = p.to_string();
        let mut xs = Vec::new();
        for v in vals {
            let mut x = parse_real(&v.replace('P', &ps)).expect("parse");
            if !v.contains(':') {
                let _ = x.set_precision(p, RoundingMode::ToEven);
            }
            let b = x
                .to_bytes()
                .unwrap_or_else(|e| panic!("to_bytes {v} p={p}: {e:?}"));
            let y =
                ExactNum::from_bytes(&b).unwrap_or_else(|e| panic!("from_bytes {v} p={p}: {e:?}"));
            assert_eq!(fmt_real(&y), fmt_real(&x), "bytes {v} p={p}");
            assert_eq!(y.inexact(), x.inexact(), "bytes inexact {v} p={p}");
            let small = x.as_raw_parts().words.len() == 1;
            let top40 = small && x.as_raw_parts().words[0] & ((1u64 << 24) - 1) == 0;
            match x.to_inline_bytes() {
                Ok(buf) => {
                    assert!(
                        top40 || i32::try_from(x.as_raw_parts().exponent).is_ok(),
                        "inline accepted {v} p={p}"
                    );
                    let y = ExactNum::from_inline_bytes(buf.as_bytes()).expect("from_inline");
                    assert_eq!(fmt_real(&y), fmt_real(&x), "inline {v} p={p}");
                }
                Err(e) => {
                    assert_eq!(e, zenith_float::Error::InvalidArgument, "inline {v} p={p}");
                    assert!(!top40, "inline refused a value that fits: {v} p={p}");
                }
            }
            xs.push(x);
        }
        let a = zenith_float::ExactNumArray::from_values(p, &xs);
        let b = a.to_bytes().expect("array to_bytes");
        let c = zenith_float::ExactNumArray::from_bytes(&b).expect("array from_bytes");
        assert_eq!(c.len(), xs.len());
        for (x, y) in xs.iter().zip(c.as_slice()) {
            assert_eq!(fmt_real(y), fmt_real(x), "array p={p}");
        }
    }
}

/// ZF-70: exact parts must stay exact (FLINT arb/acb at 4000 bits; short dyadics exact).
#[test]
fn zf_70_complex_exact_parts() {
    let mut cc = Consts::new().expect("Consts");
    let rows: &[(&str, usize, &str, &str, &str, &str)] = &[
        ("cbrt", 3, "0", "1p0", "0.866025403784438646763723170752936183471402626905190314027903489725966508454400018540573093378624287837813070707703351514984972547499476239405827756047186824264046615951152791033987410050542337461632507656171633451661443325336127334460918985613523565830183930794009524993268689929694733825173753288025378309174064803050473801093595162541572914761979916498894912254144357231916458673612081992293927698833979031909176833055421586890447189158051044152762450835011760355572144347995478182898543584249036449746648242141510393204301994369348768791158658915697996491503919351438526956684781656051853632009624553384115599644187820570711008371376051186497135415529949229737993832144448898073918979195114427426451788016926404032190986172330529844861436432632076911332349210010597742077639220590643267253517595825008344647207740423035638571999881463417314788719180947555063574319373488272991225894275487689506940332480955981111478555277621461861596098869131280815734421016426858341469324805958524869418197747969072878835926686816562955449827712312417393598802617998884596161785110152651420192957707485536214779603353101254760087981592936383179987641831715540075332926855323664269312961130291110255201840135148752399364039730829050208526340970009547866731087971946835124660211345517184906231860055592630542134455149860155601050031753588187291202601923777598639966898807453053947492772111663002009425651815780572443423640794644081622593632533322692438799581288321436055620421034008381758550051471590357757595480821230453519704064608451758746480682002069835215520626816163515461288663980146188298872768547734557871070212115390996163808700953212256227438758431348052666847108106802697302122827070064266503906816724928368203351988671114905980521462767250800702970239773577267274203635868828131183274583320116503006632872035053918429234226156580233877240250201120319955985181109301460119433575355508584700148434379831750020447658107116712613978417033506735092951014180358381073871746724781797904041065221293234734261305454131676504378306301730336097702027992064456487997405000386028720115023836629400045757185744737722439578595647329541785436980757768898820131034185424023036484691359792844879879813052079576328788895391167490283920011450766026069467686887768283213523413437144981720697871833036872222791543238946606492651074098697670739085258307476275881645996849782872261319556259546770693494683408715469113212368463101036495483915577065973242188957729957961964143857347574637013204606822827020822290745100972874709652634501306986323040532535719801603038755297093899142396993097624982082606569857646799710948708323537593617894314733054280085214434802899197026453203715405916694338940781317933578004198380122674614971969433529877157721471715478629235494107731555630383387033932285789030323723749875177272279656643274594924668286373813148707369117843457418915681680641813951912420081903335803589924364277921465674613046620282978776825568377322193917141656733322277090195410949491647313172508085561008464809730084660311", "0.5"),
        ("cbrt", 3, "0", "-1p0", "0.866025403784438646763723170752936183471402626905190314027903489725966508454400018540573093378624287837813070707703351514984972547499476239405827756047186824264046615951152791033987410050542337461632507656171633451661443325336127334460918985613523565830183930794009524993268689929694733825173753288025378309174064803050473801093595162541572914761979916498894912254144357231916458673612081992293927698833979031909176833055421586890447189158051044152762450835011760355572144347995478182898543584249036449746648242141510393204301994369348768791158658915697996491503919351438526956684781656051853632009624553384115599644187820570711008371376051186497135415529949229737993832144448898073918979195114427426451788016926404032190986172330529844861436432632076911332349210010597742077639220590643267253517595825008344647207740423035638571999881463417314788719180947555063574319373488272991225894275487689506940332480955981111478555277621461861596098869131280815734421016426858341469324805958524869418197747969072878835926686816562955449827712312417393598802617998884596161785110152651420192957707485536214779603353101254760087981592936383179987641831715540075332926855323664269312961130291110255201840135148752399364039730829050208526340970009547866731087971946835124660211345517184906231860055592630542134455149860155601050031753588187291202601923777598639966898807453053947492772111663002009425651815780572443423640794644081622593632533322692438799581288321436055620421034008381758550051471590357757595480821230453519704064608451758746480682002069835215520626816163515461288663980146188298872768547734557871070212115390996163808700953212256227438758431348052666847108106802697302122827070064266503906816724928368203351988671114905980521462767250800702970239773577267274203635868828131183274583320116503006632872035053918429234226156580233877240250201120319955985181109301460119433575355508584700148434379831750020447658107116712613978417033506735092951014180358381073871746724781797904041065221293234734261305454131676504378306301730336097702027992064456487997405000386028720115023836629400045757185744737722439578595647329541785436980757768898820131034185424023036484691359792844879879813052079576328788895391167490283920011450766026069467686887768283213523413437144981720697871833036872222791543238946606492651074098697670739085258307476275881645996849782872261319556259546770693494683408715469113212368463101036495483915577065973242188957729957961964143857347574637013204606822827020822290745100972874709652634501306986323040532535719801603038755297093899142396993097624982082606569857646799710948708323537593617894314733054280085214434802899197026453203715405916694338940781317933578004198380122674614971969433529877157721471715478629235494107731555630383387033932285789030323723749875177272279656643274594924668286373813148707369117843457418915681680641813951912420081903335803589924364277921465674613046620282978776825568377322193917141656733322277090195410949491647313172508085561008464809730084660311", "-0.5"),
        ("cbrt", 3, "-1p0", "0", "0.5", "0.866025403784438646763723170752936183471402626905190314027903489725966508454400018540573093378624287837813070707703351514984972547499476239405827756047186824264046615951152791033987410050542337461632507656171633451661443325336127334460918985613523565830183930794009524993268689929694733825173753288025378309174064803050473801093595162541572914761979916498894912254144357231916458673612081992293927698833979031909176833055421586890447189158051044152762450835011760355572144347995478182898543584249036449746648242141510393204301994369348768791158658915697996491503919351438526956684781656051853632009624553384115599644187820570711008371376051186497135415529949229737993832144448898073918979195114427426451788016926404032190986172330529844861436432632076911332349210010597742077639220590643267253517595825008344647207740423035638571999881463417314788719180947555063574319373488272991225894275487689506940332480955981111478555277621461861596098869131280815734421016426858341469324805958524869418197747969072878835926686816562955449827712312417393598802617998884596161785110152651420192957707485536214779603353101254760087981592936383179987641831715540075332926855323664269312961130291110255201840135148752399364039730829050208526340970009547866731087971946835124660211345517184906231860055592630542134455149860155601050031753588187291202601923777598639966898807453053947492772111663002009425651815780572443423640794644081622593632533322692438799581288321436055620421034008381758550051471590357757595480821230453519704064608451758746480682002069835215520626816163515461288663980146188298872768547734557871070212115390996163808700953212256227438758431348052666847108106802697302122827070064266503906816724928368203351988671114905980521462767250800702970239773577267274203635868828131183274583320116503006632872035053918429234226156580233877240250201120319955985181109301460119433575355508584700148434379831750020447658107116712613978417033506735092951014180358381073871746724781797904041065221293234734261305454131676504378306301730336097702027992064456487997405000386028720115023836629400045757185744737722439578595647329541785436980757768898820131034185424023036484691359792844879879813052079576328788895391167490283920011450766026069467686887768283213523413437144981720697871833036872222791543238946606492651074098697670739085258307476275881645996849782872261319556259546770693494683408715469113212368463101036495483915577065973242188957729957961964143857347574637013204606822827020822290745100972874709652634501306986323040532535719801603038755297093899142396993097624982082606569857646799710948708323537593617894314733054280085214434802899197026453203715405916694338940781317933578004198380122674614971969433529877157721471715478629235494107731555630383387033932285789030323723749875177272279656643274594924668286373813148707369117843457418915681680641813951912420081903335803589924364277921465674613046620282978776825568377322193917141656733322277090195410949491647313172508085561008464809730084660311"),
        ("cbrt", 3, "-1p0", "-0", "0.5", "-0.866025403784438646763723170752936183471402626905190314027903489725966508454400018540573093378624287837813070707703351514984972547499476239405827756047186824264046615951152791033987410050542337461632507656171633451661443325336127334460918985613523565830183930794009524993268689929694733825173753288025378309174064803050473801093595162541572914761979916498894912254144357231916458673612081992293927698833979031909176833055421586890447189158051044152762450835011760355572144347995478182898543584249036449746648242141510393204301994369348768791158658915697996491503919351438526956684781656051853632009624553384115599644187820570711008371376051186497135415529949229737993832144448898073918979195114427426451788016926404032190986172330529844861436432632076911332349210010597742077639220590643267253517595825008344647207740423035638571999881463417314788719180947555063574319373488272991225894275487689506940332480955981111478555277621461861596098869131280815734421016426858341469324805958524869418197747969072878835926686816562955449827712312417393598802617998884596161785110152651420192957707485536214779603353101254760087981592936383179987641831715540075332926855323664269312961130291110255201840135148752399364039730829050208526340970009547866731087971946835124660211345517184906231860055592630542134455149860155601050031753588187291202601923777598639966898807453053947492772111663002009425651815780572443423640794644081622593632533322692438799581288321436055620421034008381758550051471590357757595480821230453519704064608451758746480682002069835215520626816163515461288663980146188298872768547734557871070212115390996163808700953212256227438758431348052666847108106802697302122827070064266503906816724928368203351988671114905980521462767250800702970239773577267274203635868828131183274583320116503006632872035053918429234226156580233877240250201120319955985181109301460119433575355508584700148434379831750020447658107116712613978417033506735092951014180358381073871746724781797904041065221293234734261305454131676504378306301730336097702027992064456487997405000386028720115023836629400045757185744737722439578595647329541785436980757768898820131034185424023036484691359792844879879813052079576328788895391167490283920011450766026069467686887768283213523413437144981720697871833036872222791543238946606492651074098697670739085258307476275881645996849782872261319556259546770693494683408715469113212368463101036495483915577065973242188957729957961964143857347574637013204606822827020822290745100972874709652634501306986323040532535719801603038755297093899142396993097624982082606569857646799710948708323537593617894314733054280085214434802899197026453203715405916694338940781317933578004198380122674614971969433529877157721471715478629235494107731555630383387033932285789030323723749875177272279656643274594924668286373813148707369117843457418915681680641813951912420081903335803589924364277921465674613046620282978776825568377322193917141656733322277090195410949491647313172508085561008464809730084660311"),
        ("nth", 3, "-1p3", "0", "1.0", "1.73205080756887729352744634150587236694280525381038062805580697945193301690880003708114618675724857567562614141540670302996994509499895247881165551209437364852809323190230558206797482010108467492326501531234326690332288665067225466892183797122704713166036786158801904998653737985938946765034750657605075661834812960610094760218719032508314582952395983299778982450828871446383291734722416398458785539766795806381835366611084317378089437831610208830552490167002352071114428869599095636579708716849807289949329648428302078640860398873869753758231731783139599298300783870287705391336956331210370726401924910676823119928837564114142201674275210237299427083105989845947598766428889779614783795839022885485290357603385280806438197234466105968972287286526415382266469842002119548415527844118128653450703519165001668929441548084607127714399976292683462957743836189511012714863874697654598245178855097537901388066496191196222295711055524292372319219773826256163146884203285371668293864961191704973883639549593814575767185337363312591089965542462483478719760523599776919232357022030530284038591541497107242955920670620250952017596318587276635997528366343108015066585371064732853862592226058222051040368027029750479872807946165810041705268194001909573346217594389367024932042269103436981246372011118526108426891029972031120210006350717637458240520384755519727993379761490610789498554422332600401885130363156114488684728158928816324518726506664538487759916257664287211124084206801676351710010294318071551519096164246090703940812921690351749296136400413967043104125363232703092257732796029237659774553709546911574214042423078199232761740190642451245487751686269610533369421621360539460424565414012853300781363344985673640670397734222981196104292553450160140594047954715453454840727173765626236654916664023300601326574407010783685846845231316046775448050040224063991197036221860292023886715071101716940029686875966350004089531621423342522795683406701347018590202836071676214774349344956359580808213044258646946852261090826335300875661260346067219540405598412891297599481000077205744023004767325880009151437148947544487915719129465908357087396151553779764026206837084804607296938271958568975975962610415915265757779078233498056784002290153205213893537377553656642704682687428996344139574366607374444558308647789321298530214819739534147817051661495255176329199369956574452263911251909354138698936681743093822642473692620207299096783115413194648437791545991592392828771469514927402640921364565404164458149020194574941930526900261397264608106507143960320607751059418779828479398619524996416521313971529359942189741664707518723578862946610856017042886960579839405290640743081183338867788156263586715600839676024534922994393886705975431544294343095725847098821546311126076677406786457157806064744749975035454455931328654918984933657274762629741473823568691483783136336128362790382484016380667160717984872855584293134922609324056595755365113675464438783428331346664455418039082189898329462634501617112201692961946016932062"),
        ("nth", 3, "-1p3", "-0", "1.0", "-1.73205080756887729352744634150587236694280525381038062805580697945193301690880003708114618675724857567562614141540670302996994509499895247881165551209437364852809323190230558206797482010108467492326501531234326690332288665067225466892183797122704713166036786158801904998653737985938946765034750657605075661834812960610094760218719032508314582952395983299778982450828871446383291734722416398458785539766795806381835366611084317378089437831610208830552490167002352071114428869599095636579708716849807289949329648428302078640860398873869753758231731783139599298300783870287705391336956331210370726401924910676823119928837564114142201674275210237299427083105989845947598766428889779614783795839022885485290357603385280806438197234466105968972287286526415382266469842002119548415527844118128653450703519165001668929441548084607127714399976292683462957743836189511012714863874697654598245178855097537901388066496191196222295711055524292372319219773826256163146884203285371668293864961191704973883639549593814575767185337363312591089965542462483478719760523599776919232357022030530284038591541497107242955920670620250952017596318587276635997528366343108015066585371064732853862592226058222051040368027029750479872807946165810041705268194001909573346217594389367024932042269103436981246372011118526108426891029972031120210006350717637458240520384755519727993379761490610789498554422332600401885130363156114488684728158928816324518726506664538487759916257664287211124084206801676351710010294318071551519096164246090703940812921690351749296136400413967043104125363232703092257732796029237659774553709546911574214042423078199232761740190642451245487751686269610533369421621360539460424565414012853300781363344985673640670397734222981196104292553450160140594047954715453454840727173765626236654916664023300601326574407010783685846845231316046775448050040224063991197036221860292023886715071101716940029686875966350004089531621423342522795683406701347018590202836071676214774349344956359580808213044258646946852261090826335300875661260346067219540405598412891297599481000077205744023004767325880009151437148947544487915719129465908357087396151553779764026206837084804607296938271958568975975962610415915265757779078233498056784002290153205213893537377553656642704682687428996344139574366607374444558308647789321298530214819739534147817051661495255176329199369956574452263911251909354138698936681743093822642473692620207299096783115413194648437791545991592392828771469514927402640921364565404164458149020194574941930526900261397264608106507143960320607751059418779828479398619524996416521313971529359942189741664707518723578862946610856017042886960579839405290640743081183338867788156263586715600839676024534922994393886705975431544294343095725847098821546311126076677406786457157806064744749975035454455931328654918984933657274762629741473823568691483783136336128362790382484016380667160717984872855584293134922609324056595755365113675464438783428331346664455418039082189898329462634501617112201692961946016932062"),
        ("nth", 3, "0", "1p3", "1.73205080756887729352744634150587236694280525381038062805580697945193301690880003708114618675724857567562614141540670302996994509499895247881165551209437364852809323190230558206797482010108467492326501531234326690332288665067225466892183797122704713166036786158801904998653737985938946765034750657605075661834812960610094760218719032508314582952395983299778982450828871446383291734722416398458785539766795806381835366611084317378089437831610208830552490167002352071114428869599095636579708716849807289949329648428302078640860398873869753758231731783139599298300783870287705391336956331210370726401924910676823119928837564114142201674275210237299427083105989845947598766428889779614783795839022885485290357603385280806438197234466105968972287286526415382266469842002119548415527844118128653450703519165001668929441548084607127714399976292683462957743836189511012714863874697654598245178855097537901388066496191196222295711055524292372319219773826256163146884203285371668293864961191704973883639549593814575767185337363312591089965542462483478719760523599776919232357022030530284038591541497107242955920670620250952017596318587276635997528366343108015066585371064732853862592226058222051040368027029750479872807946165810041705268194001909573346217594389367024932042269103436981246372011118526108426891029972031120210006350717637458240520384755519727993379761490610789498554422332600401885130363156114488684728158928816324518726506664538487759916257664287211124084206801676351710010294318071551519096164246090703940812921690351749296136400413967043104125363232703092257732796029237659774553709546911574214042423078199232761740190642451245487751686269610533369421621360539460424565414012853300781363344985673640670397734222981196104292553450160140594047954715453454840727173765626236654916664023300601326574407010783685846845231316046775448050040224063991197036221860292023886715071101716940029686875966350004089531621423342522795683406701347018590202836071676214774349344956359580808213044258646946852261090826335300875661260346067219540405598412891297599481000077205744023004767325880009151437148947544487915719129465908357087396151553779764026206837084804607296938271958568975975962610415915265757779078233498056784002290153205213893537377553656642704682687428996344139574366607374444558308647789321298530214819739534147817051661495255176329199369956574452263911251909354138698936681743093822642473692620207299096783115413194648437791545991592392828771469514927402640921364565404164458149020194574941930526900261397264608106507143960320607751059418779828479398619524996416521313971529359942189741664707518723578862946610856017042886960579839405290640743081183338867788156263586715600839676024534922994393886705975431544294343095725847098821546311126076677406786457157806064744749975035454455931328654918984933657274762629741473823568691483783136336128362790382484016380667160717984872855584293134922609324056595755365113675464438783428331346664455418039082189898329462634501617112201692961946016932062", "1.0"),
        ("nth", 2, "-1p0", "0", "0.0", "1.0"),
        ("nth", 4, "1p4", "0", "2.0", "0.0"),
        ("nth", 3, "-2p0", "2p0", "1.0", "1.0"),
        ("nth", 5, "-4p0", "-4p0", "1.26007351067010099621873275056220120226585322402863165351502995578278340149464472657239641306550695041827676158909645008505038312488146839221180230069654555600247198126837082417103626225166148219773794603281608262991958788098193625572758872880583370944793140302323952873769794328137509589514277718471660762388249199858345296245365949807662481865886142358959320454155557805761925727681421541374788601010558976479574483078597494010755420772275581585135921459707030178492245087616807217612310455900569419064377559832930697353837143216090319132037712308729594437196726932629383981219929878422034161442978149731967312876482903364050961143319535896855171597013623799783261550670572039699450975195983766281995799356658857237501491376471552802317203738896643247226031445647527937196265105290830130293234014555165890144513054475512109303943976146914996597187388197513009047311794686586955697496791259309258888987125247741659194559153569042647319375160312601879973494454309576438555440658869840171304667488634463232230459516913181701071806524318226258825638432111211403706001884728076997276351839674618619253920814856902114418898070522291674960405880998643671320973774770201052697040055255411479767514307666405926051523924176568357640395640659417816447037161118238435127318308422490913605118474247430216116226723504952155743031282826109755941351876545609420373772842605717066735969712252444981456689680388645853590966884340489068326055169802735163633314331010922295304685342230108813046815586468230895833479069896987170207802353863442444518852362401028641665305530834952839310938587881638417309718835092358312353947190297356247513128914724187000770443352548328070804597905015195430786123779116473450818036704696700067043482490849977479446734278791524885344806399767901201723628663595302971606569452387020026531880626072395590821591123418298831467772017936626712011671934370572630979103937567048643254651150744376009151405115149380792331147141920118916007532851532847084559141253711301447256941648254333273148788834715315465398491055907435446845406914877774313487989191432908921029356427810746149251620442186221213414449560194439427605334354932584161171755369323050118908709426200783543119929041008177465863410865790196739165239098923612678173239183581206453050946017275174079211828982134364896376324360876017104796079292808448939487895075649288659110601666497158122263682771403457889617518692770193674772196664942415228128575708251780397569280448227685551394061413657108170997092583855317098458344917004000116895764143988443381801173332582044050727401434308109530068604780196915544615556273228201332966188483664317438318766342955192417264590638142857069511917230173487264616719933997748156999709949029230381350615776027545078196906090690811524230322145031566715850511727813846881579050875634204346625569976233007970378915745101249134899396432148205529229955518258889663185442278663603672929251661840052595729262634382477714046278723806286273378626647362648837232787994701836331716331631605891607776834074390953", "-0.642039521920206148014145916196563084545544044222868791379581333077522938675742276865189208876115812933522673513709558332923749262657931460418622239929819201569082894672431241265397939590341553907470065357607315737417876184774904034684372459257207413134316964873480827517357137393420641145896591489351743174641447677812108015404093651288114944319419211044716138063463989596869268564413450196690368031271423508546337241716270939826343780100984698073306061425660184738256459409188198999987503888130887089325823229386554778694840864330880312750398825303817485589678158350043809542583681781715698231992016636291290001931791790713756980267856111651840594811317599484516411666653498739081642999013288594706786394223352598880029823355239314800191020127895836355804015358312373819877403149384241101204342673798231125517344288560800265534308543312864453747064487262788638770880588808141695374099074297862129625083450688424392918513845300894276373890088747987299993837617007090089327649028175228935526321272408421554915218713525756418552170575586743584250825231609818727636894023728647848609078332042780535191231508845071644723134916492520518700308470047874610786047496315206381190443571329657325894801345399527056351266920425239081530933252417409226949340923803580909572243540128629577064397510379177892226694138004425600586295273062278341392878236400904854261802418708647846348815176079007591377983610939562315886654405702191453182520845415544091140255605088126110553226286375900522689309206320806569199857367746766362845887295441787156081818590017463764017493218051541968381788153161298018221752602658864982088201896949672313053556902247516918705338508679999586102774230508077549228089805445620530026388972205492948317316639041135474010821798216693175306994313950171429775252854898585806628858688056388953628867704954201066724465892673858773977218231611405489362492926411779506117841317052147249260520203483039946954546653701687162063614960039580789480340655258768951492548529468652763235243773351899081350730687038350839678648178560033294972559092876397524761884009265274995022633419298173922171153641047666252153861988961596134878809622013760682744761471900694759618311834520433686143883631324691421323132742784204413209106241327950263634604567395946898005598961460303957741587470972532633037373725501460711638554818275477284679958988077969078774676512237737385402326779248529103229195345500678643123735723895431039915320933415749257137843301400428937429553590689606105539141039968631253297147544591240020437424375480322157244782728062607936878903464850268194207108780853464406878724278418680114003588238912140022864437972115552131401842292601833110779635505803640432729878344227592100040534754285007813667028949955708150717263765476102931023669555537494467089226100146512334534708571791387136727101171329171081986179887223231357344017893289032998747208507829569680200208765544673711097793910510403870005264636476127142306531765251489253914609520877587557889474619534581938916624571635952586297525043174280"),
        ("exp2", 0, "1p0", "0", "2.0", "0.0"),
        ("exp2", 0, "1p0", "-0", "2.0", "-0"),
        ("exp2", 0, "1p1", "0", "4.0", "0.0"),
        ("exp2", 0, "5p1", "-0", "1024.0", "-0"),
        ("exp2", 0, "-3p0", "0", "0.125000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000", "0.0"),
        ("exp10", 0, "1p0", "0", "10.0", "0.0"),
        ("exp10", 0, "1p1", "-0", "100.0", "-0"),
        ("exp10", 0, "5p1", "0", "10000000000.0", "0.0"),
        ("log2", 0, "0", "1p-60", "-60.0", "2.26618007091359690481384147285833340508593073386189779209300827397030047686065255112954193980200804482653323135388951403811127118581148075087270222783191740557844345953986336849528862553327474290916734597508770443847390594092789601370188719805158356206956902779669440667091453542952316343583518022939293354520619181132165906329432211645494501531778719242358305109452152486005915069524964461728353894927196283747385212197464417609420327488535092269701983704033067342012299083556121215109874319915287748608301752123737607556360395468155482592075117296048594809781855413508703672317681314324135421535931761756746451084125873139939902711395153459628199290864212875468256359105801088324235324782476792239840033519944202311670990932307146232746778531995255012144964033979089525968629525724053346811492227109980784555235490792067912170217465196682384585564834654983122052403066033707908483486295094314637412306279707099837961938729332328341169225477138648341562170933987229805152379947869393427569265178181035536993095358400916414895738045440679544725329629748383446957543586022114331395293915471010467715441078763887856207769447731893035194039911488884275787992270256934488545001600736567547926349219775391712103827985095468445227854089399243647012689702844850389934588240268798127738623321009276970883765637627268605950245579511352233617756425109269080912895197353425403591588584139972500133086487213626889304428692446412220742459910397613125271089975855019689336661588852887977345451504265561387318279342052715372433169676112777118622451544796309456035961381112320886851634969490240887414504677171996921291040845307299970633965381976574612528563378947016334834968041511570501383927524153505201980554628170093174565618023210087196956245751290431325022196231077700785299726537104494168396337431841450096106365687921280294829711222946270302755036660452675302607754280495006334451083432466477927231345756161999972745670177919623747945517517034030139510902657110958331025216494174411352736455430434545810138107544753803802679098705677969326547289372016534352083870316384434213367709945685722096968408144487312123610083614230780096941306264740896420979226207375455028161587980578144893272438811137359218331195273964713455742953944432113844196761903165370494130741705485368175327181772380879612274334237303669464121487731676488832649920369266038966688685874742066680247002021725997752685038006217602662000376770678301263384826405799599021344872708148442482256663729363298889329714740140511271542152861536964269701421499907418993908616176258710900958397435672299230114316624140469060231854554398341826466492932996027078080516816077978071838505399778664736320949178490987419172506646820087186965934526700813715774829173701232141523100107415526380534423992822931144501454002862918617485164653790368555537644480337540337907529358714518477326748824095315359170658107484672026295756297129211357139817011481265465892267318886818508865571898608562537784931937362437004479909896816476662527470589847073278879918565524654"),
        ("log2", 0, "1p-60", "0", "-60.0", "0.0"),
        ("log2", 0, "1p-60", "-0", "-60.0", "-0"),
        ("log2", 0, "1p-60", "1p-60", "-59.5", "1.13309003545679845240692073642916670254296536693094889604650413698515023843032627556477096990100402241326661567694475701905563559290574037543635111391595870278922172976993168424764431276663737145458367298754385221923695297046394800685094359902579178103478451389834720333545726771476158171791759011469646677260309590566082953164716105822747250765889359621179152554726076243002957534762482230864176947463598141873692606098732208804710163744267546134850991852016533671006149541778060607554937159957643874304150876061868803778180197734077741296037558648024297404890927706754351836158840657162067710767965880878373225542062936569969951355697576729814099645432106437734128179552900544162117662391238396119920016759972101155835495466153573116373389265997627506072482016989544762984314762862026673405746113554990392277617745396033956085108732598341192292782417327491561026201533016853954241743147547157318706153139853549918980969364666164170584612738569324170781085466993614902576189973934696713784632589090517768496547679200458207447869022720339772362664814874191723478771793011057165697646957735505233857720539381943928103884723865946517597019955744442137893996135128467244272500800368283773963174609887695856051913992547734222613927044699621823506344851422425194967294120134399063869311660504638485441882818813634302975122789755676116808878212554634540456447598676712701795794292069986250066543243606813444652214346223206110371229955198806562635544987927509844668330794426443988672725752132780693659139671026357686216584838056388559311225772398154728017980690556160443425817484745120443707252338585998460645520422653649985316982690988287306264281689473508167417484020755785250691963762076752600990277314085046587282809011605043598478122875645215662511098115538850392649863268552247084198168715920725048053182843960640147414855611473135151377518330226337651303877140247503167225541716233238963615672878080999986372835088959811873972758758517015069755451328555479165512608247087205676368227715217272905069053772376901901339549352838984663273644686008267176041935158192217106683854972842861048484204072243656061805041807115390048470653132370448210489613103687727514080793990289072446636219405568679609165597636982356727871476972216056922098380951582685247065370852742684087663590886190439806137167118651834732060743865838244416324960184633019483344342937371033340123501010862998876342519003108801331000188385339150631692413202899799510672436354074221241128331864681649444664857370070255635771076430768482134850710749953709496954308088129355450479198717836149615057158312070234530115927277199170913233246466498013539040258408038989035919252699889332368160474589245493709586253323410043593482967263350406857887414586850616070761550053707763190267211996411465572250727001431459308742582326895184277768822240168770168953764679357259238663374412047657679585329053742336013147878148564605678569908505740632732946133659443409254432785949304281268892465968681218502239954948408238331263735294923536639439959282762327"),
        ("log2", 0, "-1p3", "0", "3.0", "4.53236014182719380962768294571666681017186146772379558418601654794060095372130510225908387960401608965306646270777902807622254237162296150174540445566383481115688691907972673699057725106654948581833469195017540887694781188185579202740377439610316712413913805559338881334182907085904632687167036045878586709041238362264331812658864423290989003063557438484716610218904304972011830139049928923456707789854392567494770424394928835218840654977070184539403967408066134684024598167112242430219748639830575497216603504247475215112720790936310965184150234592097189619563710827017407344635362628648270843071863523513492902168251746279879805422790306919256398581728425750936512718211602176648470649564953584479680067039888404623341981864614292465493557063990510024289928067958179051937259051448106693622984454219961569110470981584135824340434930393364769171129669309966244104806132067415816966972590188629274824612559414199675923877458664656682338450954277296683124341867974459610304759895738786855138530356362071073986190716801832829791476090881359089450659259496766893915087172044228662790587830942020935430882157527775712415538895463786070388079822977768551575984540513868977090003201473135095852698439550783424207655970190936890455708178798487294025379405689700779869176480537596255477246642018553941767531275254537211900491159022704467235512850218538161825790394706850807183177168279945000266172974427253778608857384892824441484919820795226250542179951710039378673323177705775954690903008531122774636558684105430744866339352225554237244903089592618912071922762224641773703269938980481774829009354343993842582081690614599941267930763953149225057126757894032669669936083023141002767855048307010403961109256340186349131236046420174393912491502580862650044392462155401570599453074208988336792674863682900192212731375842560589659422445892540605510073320905350605215508560990012668902166864932955854462691512323999945491340355839247495891035034068060279021805314221916662050432988348822705472910860869091620276215089507607605358197411355938653094578744033068704167740632768868426735419891371444193936816288974624247220167228461560193882612529481792841958452414750910056323175961156289786544877622274718436662390547929426911485907888864227688393523806330740988261483410970736350654363544761759224548668474607338928242975463352977665299840738532077933377371749484133360494004043451995505370076012435205324000753541356602526769652811599198042689745416296884964513327458726597778659429480281022543084305723073928539402842999814837987817232352517421801916794871344598460228633248280938120463709108796683652932985865992054156161033632155956143677010799557329472641898356981974838345013293640174373931869053401627431549658347402464283046200214831052761068847985645862289002908005725837234970329307580737111075288960675080675815058717429036954653497648190630718341316214969344052591512594258422714279634022962530931784534637773637017731143797217125075569863874724874008959819793632953325054941179694146557759837131049307"),
        ("log10", 0, "6p0", "8p0", "1.0", "0.402719196273373142077900460387457894508681695348164852149756118784031761699616919347122335881865726958199911934782354774379982920586968502890606488334121302667417665560476730579085435283587107781441890644260236278841484747576609053753560890966483627076359822578704042876383138419095926232512833407851951629412006925245263112528016650654305163497081897309917749591380854433921025391705829875349684727988180980654816162063050301883405327565924994640751793274973520500821338628887310803990024211874587948164172845263191709966357517124530182057495268321176478575317110272325470896564524187912489678624935910560887366829859773085862948664250673715248947780670450042343478227343249253036165385974058749867471078944990685825686817239897380231637127966181077421527081178462113887238298764815496232011181062548536466791291693097562772524547145705107573368738216561729705219370297628388622618221175968545576774887328929302420574753590109822234389955736054277716765642370181499081614655686175132312498052451122017818186814097952051632401762944116651402833259800647364171991374639610007234828269706979968219599907880803342410517047552800085090311975870238939418044266541307416017178933259514390094373079874334037069852739621952960974134548350575291795892459558593709345663835118909126500074175590087517531224624173046478110943567861370013199375974140017586509389351344750925310758390704253461686848195552816664963638657077026300551274960484970072585780728752997521327498607748796877180923638779846070537419766484096059342285405000327919180914578575393992630640502838568164671602296142089715541611321555906313373089399401972697414330629397669288470027932878748125080428881171608978949872426072043743810660189689704490202957231341597505862950653768933037355263463896130694452396618962171099677117056714869437698835917261679986230716003271316065062826968621909985204184277824701024042275335241831514406215137087860933456327016030227931570262933236743193723486242390663848712257142803870401538045689498589198749211273815666324462521384056543705022647048602420377452573105174133054453270881317400460464236444876339800314677160056529430682997821671849860388509350466524294620138515810291637103380008076102006280406592690658532538892373912167349320919736472778964509782651583628053900725490218986933477283327113413265679487709548816211132876009000035055571583286815664988298186330106028478965518175822471670118228049180548582081292059924228071544816897620445829932177382060370089434016388771064872334867467610761000637944824073265283386244444972934581939833463797064450793737952473525891204830775313429149924476087460191347564518175820406984426697084249837999782216332341273529051850336042297206917686192086615696856900525839086843701229565229886095889358080915053311178577718035858904117670640977433195446699440555336726089859186184300774789581654441910696562908327429358110851255328696681990991131673411307037488947215984210346268944821438568128588807655141245140779492079560134890392742977845942996003921007756152331"),
        ("log10", 0, "1p0", "-0", "0.0", "-0"),
    ];
    for p in [64usize, 128, 256, 1024, 2048, 4096, 8192] {
        for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
            for &(f, n, re, im, wr, wi) in rows {
                let z = ExactComplex::new(parse_real(re).expect("re"), parse_real(im).expect("im"));
                let w = match f {
                    "exp2" => z.exp2(p, rm, &mut cc),
                    "exp10" => z.exp10(p, rm, &mut cc),
                    "log2" => z.log2(p, rm, &mut cc),
                    "log10" => z.log10(p, rm, &mut cc),
                    "cbrt" => z.cbrt(p, rm, &mut cc),
                    _ => z.nth_root(n, p, rm, &mut cc),
                };
                for (name, got, want) in [("re", w.re(), wr), ("im", w.im(), wi)] {
                    let r = cr_ref(want, p, rm, &mut cc);
                    assert_eq!(
                        got.compare(&r),
                        Some(0),
                        "{name} {f}({re},{im}), n={n} p={p} {rm:?}: {} want {r}",
                        fmt_c(&w)
                    );
                    if matches!(want, "0" | "-0" | "0.0") {
                        assert!(got.is_zero(), "{name} {f}({re},{im}) zero p={p} {rm:?}");
                        assert_eq!(
                            got.is_negative(),
                            want == "-0",
                            "{name} {f}({re},{im}) zero sign p={p} {rm:?}"
                        );
                    }
                }
            }
        }
    }
}

/// ZF-70: terminating real 2F1 with an exactly representable value came back 1 ulp off in
/// directed modes (the Ziv loop cannot certify an exact value). The polynomial is now summed
/// exactly and rounded once. References: exact rational sums (Python fractions).
#[test]
fn zf_70_hypergeom_2f1_terminating_exact() {
    let mut cc = Consts::new().expect("Consts");
    let rows: &[(&str, &str, &str, &str, &str)] = &[
        ("-1p-1", "-3p0", "1p1", "1p1", "2.125"),
        ("-3p0", "-1p-1", "1p1", "1p1", "2.125"),
        ("5p1", "-3p0", "1p0", "-1p-1", "84.75"),
        ("1p-1", "-3p0", "1p1", "1p1", "0.375"),
        ("3p0", "-3p0", "1p0", "1p-1", "-0.25"),
        ("-5p1", "-1p-1", "1p1", "-5p1", "-29052098.70703125"),
        ("c9p-1", "-3p0", "1p0", "1p0", "-159262.1875"),
        ("1p0", "-3p0", "1p0", "-7dp3", "1003003001"),
        ("-1p0", "1p0", "1p0", "1p0", "0"),
        ("-3p0", "3p0", "1p0", "1p0", "0"),
        ("-2p0", "1p0", "3p0", "1p-1", "0.708333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333"),
    ];
    for p in [64usize, 128, 256, 1024, 8192] {
        for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
            for &(a, b, c, z, want) in rows {
                let (a, b, c, z) = (
                    parse_real(a).expect("a"),
                    parse_real(b).expect("b"),
                    parse_real(c).expect("c"),
                    parse_real(z).expect("z"),
                );
                let v = a.hypergeom_2f1(&b, &c, &z, p, rm, &mut cc);
                let w = cr_ref(want, p, rm, &mut cc);
                assert_eq!(
                    v.compare(&w),
                    Some(0),
                    "2F1({a},{b};{c};{z}) p={p} {rm:?}: {v} {:?} want {w}",
                    v.err()
                );
            }
        }
    }
}

/// ZF-30: ln_gamma(x) for small positive x was ln(Γ(x)). For x near EXPONENT_MIN, Γ(x) ≈ 1/x
/// overflows although ln Γ(x) ≈ −ln x is finite (≈ 3.2e18), so the result was Inf/MAX. Near
/// x = 1 and x = 2, ln(Γ(x)) at p bits has only absolute accuracy 2^-p while ln Γ(x) is tiny,
/// so the last ~40 bits were wrong (Ziv accepted them). Now: −ln x − γx for tiny x,
/// ln Γ(1+x) − ln x for x < 1/2, −γ(x−1) / (1−γ)(x−2) at the poles of cancellation, and
/// otherwise extra working bits log2(1/|x−1|), log2(1/|x−2|). ln Γ(1) = ln Γ(2) = +0 (C99).
/// References: FLINT arb lgamma at 4000+ bits (zf30/refs.py, refs2.py), rounded once per mode.
#[test]
fn zf_30_ln_gamma_small_and_near_one_two() {
    let mut cc = Consts::new().expect("Consts");
    let rows: &[(&str, &[usize], &str)] = &[
        ("minnorm:64", &[64], "3196577161300663914.9471576749356129486404926450780118465153"),
        ("minpos:64", &[64], "3196577161300663958.6154300502121674419261162969431356352718"),
        ("minnorm:128", &[128], "3196577161300663914.94715767493561294864049264507801184651532475839641542212010"),
        ("minpos:128", &[128], "3196577161300664002.97684960604866724462897207026643599210384182214883269544646"),
        ("minnorm:256", &[256], "3196577161300663914.94715767493561294864049264507801184651532475839641542212010034162864247995101201484682501613286513"),
        ("minpos:256", &[256], "3196577161300664091.69968871772166685003468361691303670576785902026150522289350276244401608222316449434197340021963045"),
        ("8000000000000000000000000000000000000000000000000000000000000001p-255", &[256], "-9.96986354946110035861583444182881795741157466946105185177113219267729426630858711200752995788853696427347082446349553e-78"),
        ("ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffp-256", &[256], "4.98493177473055017930791722091440897870578733473052592588556609633864713315447758318502867290804906080341321815994684e-78"),
        ("minnorm:1024", &[1024], "3196577161300663914.947157674935612948640492645078011846515324758396415422120100341628642479951012014846825016132865126716987071539864391008396372544241576064067256311206662675069081350568244357647313187128220955325029632719073305273367960593968046153629925695757930255514844745638198965378077917542631524803951137116236064535868898210839116350125173"),
        ("minpos:1024", &[1024], "3196577161300664624.036723387759664482468952896792640987751962208937540387575750053370317754948706079645008533469182482184502155583675896568683573724809135908489435233889371730203386307992222557331413699369402916284216800736954713482988206695621297309132579797560954990493103059064457484251153352691268183248647077862696348903189287097405546313698714"),
        ("8000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001p-1023", &[1024], "-6.421737433466267479214569199268797732325729828269747790409945209499784131942992628680923650241863666850653572001572449424990263883993607352215937949781908015336293305190364728930784888548682243740132501275628379423858582851563685540958110047039834453305961942725912695591949684960311599991606712807156822243448698601621250260016713871195791790539771e-309"),
        ("ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffp-1024", &[1024], "3.210868716733133739607284599634398866162864914134873895204972604749892065971496314340461825120931833425326786000786224712495131941996803676107968974890954007668146652595182364465392444274341121870066250637814189711929291425781842770479055023519917226652980971362956347795974842480155799995803356403578411121800699229229976518926847640136824850112675e-309"),
        ("3p-3", &[64, 256, 1024], "0.86307398227064746240508909413401549533247062934842713501214158155473398936618356558114368929134299687589643426085011770515992375830666255075703713274325370506086803686340963785071788805972925729520097677008461880093595991208892001353987919836314143040186149824626596671208396608051567115031753346224702579107454771944193818560214093395001627311916757656042103325117111368272012792754102997897661825460663575388203044836396398666878873135829397245557317019147842876708529467687291214699138798593134790154590863982132549318627801132134353964429271658988759005728070840287396773512110891736183424432782270984964757130444389233676667809900067381151621765537519384475655246902704079384384870018818962513896757528602871228724439017551591886187210396245199961950977887830027827948905370974318513407391373130167538109221755909678710227177149688354556819451199560010082748806931988157284754972353090410914685360042654465265489225911057835495839421810163009713051296685906131198851568033146042213970497103748378888005607283090466450566702633109826254720091838574352714868226076081840831837334188256858948317685"),
        ("1p-130", &[64, 256, 1024], "90.109133472792890224240175789562953849814593395737371859179254552081188100121122978844220730527749425165172800112996958348851836082213313610485441707138546996084560103630180543063539116197084509927963474963986434511426598303950961459162344146430197278383118754271456045161593740278587747203817166540423056583682295778531694494810999957998431529378930415410893537564619341595428833213893095209841990602050676539543580987010918971286667678457176808440452263309161644590991905104685960966400497560960354762786027312701571143696242556332709866858417031177271328881170641999256205807399668906921522955603203754777336957638851196441075971251394495930791451057785876188407746938888444879727568198221664340796759642928454963733919645751959665879015422133258340058964680120132229269737534249222509606137825851435762363843921316222074169906556021135922435123356263635327695918613122042228941869632055475470288993538990244391474465157515615385950935339225380559774133461917744485051112458614279023534663301597138050137189063269564257082854809024145601948016350601197458791198061139973443781246105433417175893430"),
        ("1p-600", &[64, 256, 1024], "415.88830833596718565033927287490594084530008061615315247240800569603617318181682936351799619785121252520088861234241144021131214145487833421960245098104557715843654249702286111418422303523065435319664665367429293567841438835935958933175210082755431318772381445820999063005343366603953456457565668094461556455504998146477092475249308514152219705609584459640227822506974465510798941348457741546976173184453331979963032847179483496726610678901408482310764188479184480288137825947677089579945416419030111085683314092694306664282791760946604490488045706473930160134662131129277533776188619780534358422728557939627629216761049533338521036987409115351618038731668942316790814078363182853289303037823213168313710736124907735795572578369888855970515383734902799825655058841898275128184270085057939848736704576010302068592970123630010309747423163427533803804340910496093628497239126381313216071934920240942968127090518859658106265646312549205474377412965924971291164314279882893500985223067416104572583064235833351987962086694126603135937335243240219618316860535289652958960416251566841780641626276397305141514"),
        ("1p-1", &[64, 256, 1024], "0.57236494292470008707171367567652935582364740645765578575681153573606888494241303989181163513774485385100490611434899457952410276396172829363954053940514341263819695713317295145124238667943496889460155981541237839700595801410861368994406328158902491184865665534750180003220274363194011163504821675247975590753311862623417169563494828987570238851928899769991292128301142425074068108957962625283538193430140381728444875256167180390719957072132147983564488905682632261725205295035804092854124905940915934488364229640551288284375862111691685946365216441086743255213807661875805141961106700718483587928082212368593902530233460281416886550668108137257949376007564982727328698457641261958479262268967973007001899782598330182690000563299292500648828495303723337277363355225424753342793716953871256707962063261588858924588997940478839402551482223754507544557016390403841686689694747440375764450459376818830433537179166725540696162677870338421637155990248169998809015231106431807979299182023790049308999691323146388231379742424482070537415662967310268175365230275153841077472220773894422797676142200239254591086"),
        ("fp-1", &[64, 256, 1024], "7.5343642367587329551583676324366857670272790219521205641257857314005935525077508561640373587281857253979065413766952156360649549290312113095323961106624562312279353616444121763869929871734132747895352074023178498757809684096190050457833378835111357261323875243553699745160556821968603417171471794747948490297276731829779409328248985756889509668178791533303798855001074850655743484179339811603673511269165520615395574651441547017102181822335217973677018052818091014630861511479935120668881965929590268899461628298306938981455015612950367818531685960338719826820423488517430416825510173369397105222329964675743996867278753670998967666861999741030579945362152130260106247707139874317765223424805256970561505861538960695102662969062236293392671785579676661323251885204091306550796981694053902957272759775812435029520208111316841508807112853406178528610709870266332294779009647982454380393779127633890723374086275964117784231929839706078679117643639021738992158373087581680988601093672660222414033781736842289857592321068626702506938919166757143532942716995771571188219149652021584995009787683896878436156"),
        ("800000000000000000000000000000000000000000000000000000000000000000000000001p-298", &[64, 256, 1024], "8.3019512666485171532222799257533338104854306714190794517392927760509898462442170447962630548484560145787318458447279192438393581576920655785680948173325357153660787191232894643520437163915339074769529484458253784326207686251230641184218901016380167307662849381036896166435912402321970682147781249325614179143283096867833767141050713679124093604232826696073445185356014233445586513358615484848955711181983046644945440463041613252684968582585014585086641862178930188271033841888410282150853407183131343701799380543488809696427134423728642953975556718163453840453398405041502483121511238005332418985663354995805498647530419276315006507422910086332472175994608944808384558784547431665274762063069415774042385574349065174836243402995036196154313277162917043947660653096819350627399992249656477066221480767010697183996640567039795604039181114896657742062699216383743114386047516479407731060056413808524178778855350256572617238777523582147551676323566571873672551319355893503635466494302952282789420750702761124402013485127194432691512313699474083772980407484377338409789359037543589278689232277102208314781e-91"),
        ("ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffp-304", &[64, 256, 1024], "1.7710035303972482904938762334674638689711611529029248138604537019839736207254134970381450772234854781442676416850690507015419447502752417844699517835217718599702621772941202134835037899299687509645004455824377101656457087081965245246939812359074585411719303905878478505371902855190973209266464788116178146051593132530677853501321666090406036548710840309964611145017661899791686825270606099479049789061038708587293615954853043696286187538289875632798381482475520639760191704656565861832929979576240727236677820893325505908967348101959368831557205244506862852071515161931720888152135057305528879995309910236588859793750082561172141317462178455557665645776237577076772197327650771645135812650351226967673691924428478777237808391855067767995269150968228321917980696678652890869553596145048453169207351889585802720825749578944940790452659971609172003234350175518771729310378461511885234433143368784618925921164945092487075014562231749107273712544420558590612533099527343167268159868409286184797198784605218554159158892653976391981925056078322632351076099412452131207402473513981650486575487887667403044444e-92"),
        ("1fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffp-300", &[64, 256, 1024], "-2.0754878166621292883055699814383334526213576678547698629348231940127474615610542611990657598265154551968522071070131873980452201351525185256547989484412801856683340833397203870017025525404261251553220693890363987580498658167591124919999044515146521324754269736232612368023126781277782774056647808037160792133014012973174344252555846539417624808804692687299511780837435958061072839843917039755589962597889239298991214934525046601026915946968622118654233394541974363698872800040766459564013609003282959522701746630587048000772688329531865423333108370188808311146846434742587569258232202184953335509923003036753111819742046923720466380075553046198475586342407935283890241298685454801603684864271647483245393719998563327438357687959598886643692985357027677120798740515497956952838646960819162546932813289672940869705714803757973222840316389206309638130977306476731880712409603125282511527697219523541677764597685190023611892398716733202351292952608900989825635769196228590650598759976659113371213617501654067634236888265740279300090293103720228697257482059438083521586155495436514864212807849865957917887e-91"),
        ("1p0", &[64, 256, 1024], "0"),
        ("1p1", &[64, 256, 1024], "0"),
    ];
    for &(xs, ps, want) in rows {
        let x = parse_real(xs).expect("x");
        for &p in ps {
            for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
                let g = x.ln_gamma(p, rm, &mut cc);
                let w = cr_ref(want, p, rm, &mut cc);
                assert_eq!(
                    g.compare(&w),
                    Some(0),
                    "ln_gamma({xs}) p={p} {rm:?}: {g} {:?} want {w}",
                    g.err()
                );
                if want == "0" {
                    assert!(
                        !g.is_negative() && !g.inexact(),
                        "ln_gamma({xs}) p={p} {rm:?} must be exact +0"
                    );
                }
            }
        }
    }
}

/// ZF-70 (follow-up): the exact terminating-2F1 path stopped at 2^20 numerator bits and sent
/// longer polynomials to the floating path, where an exactly representable value could not be
/// certified in directed modes. It is now exact at any degree (binary splitting). Rows use
/// 2F1(−N, b; b; z) = (1 − z)^N: z = 1/2 gives 2^-N, z = 2 gives (−1)^N, z = 3 gives (−2)^N,
/// all exact, so every mode must return the value itself with the inexact flag clear.
#[test]
fn zf_70_hypergeom_2f1_long_terminating_exact() {
    let mut cc = Consts::new().expect("Consts");
    let b = parse_real("3p-2").expect("b");
    for n in [2000i64, 3001, 20000] {
        let a = parse_real(&format!("-{n:x}p0")).expect("a");
        for (z, e, neg) in [("1p-1", -n, false), ("1p1", 0, n % 2 == 1), ("3p0", n, n % 2 == 1)] {
            let zz = parse_real(z).expect("z");
            let mut want = parse_real(&format!("1p{e}")).expect("want");
            if neg {
                want = -want;
            }
            for p in [64usize, 256] {
                for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
                    let v = a.hypergeom_2f1(&b, &b, &zz, p, rm, &mut cc);
                    assert_eq!(
                        v.compare(&want),
                        Some(0),
                        "2F1(-{n},3/4;3/4;{z}) p={p} {rm:?}: {v} {:?}",
                        v.err()
                    );
                    assert!(
                        !v.inexact(),
                        "2F1(-{n},3/4;3/4;{z}) p={p} {rm:?} flagged inexact"
                    );
                }
            }
        }
    }
}

/// ZF-70: no size-threshold refusal. The old 2^46-bit estimate (Paul's addition,
/// mislabelled as a Jeff ruling) refused when `n * 3 * (w0 + 2*WORD_BIT_SIZE) > 2^46`.
/// A degree just past that estimate, at high input precision, with small integer
/// inputs: ₂F₁(−n, 1; 1; 1/2) = 2^{−n}, a tiny exact dyadic. It must compute, not
/// return MemoryAllocation or NaN.
#[test]
fn zf_70_terminating_2f1_no_size_threshold() {
    let mut cc = Consts::new().expect("Consts");
    let p_in = 1 << 30;
    let dest = 64usize;
    let rm = RoundingMode::ToEven;
    let one = ExactNum::from_u8(1, p_in);
    let half = ExactNum::from_u8(1, p_in).ldexp(-1, p_in, rm);
    let w0 = one
        .mantissa_max_bit_len()
        .expect("one")
        .max(half.mantissa_max_bit_len().expect("half"));
    let per = 3 * (w0 + 2 * WORD_BIT_SIZE);
    let n = (1usize << 46) / per + 1;
    let a = ExactNum::from_i64(-(n as i64), p_in);
    let v = a.hypergeom_2f1(&one, &one, &half, dest, rm, &mut cc);
    let want = ExactNum::from_u8(1, dest).ldexp(-(n as i64), dest, rm);
    assert!(
        !v.is_nan(),
        "2F1(-{n},1;1;1/2) must compute, not error: {v} {:?}",
        v.err()
    );
    assert_eq!(
        v.compare(&want),
        Some(0),
        "2F1(-{n},1;1;1/2) p={dest}: {v} want {want} {:?}",
        v.err()
    );
    assert!(!v.inexact(), "2F1(-{n},1;1;1/2) flagged inexact");
}

/// ZF-08 (follow-up): values the 16-byte inline record cannot hold were refused with
/// MemoryAllocation although nothing is allocated; it is a format limit, now InvalidArgument.
/// Both cases: a mantissa wider than BINARY_INLINE_MANT_BITS, and |exponent| >= 2^31 with
/// bits below the top BINARY_INLINE_WIDE_MANT_BITS. to_bytes still encodes both.
#[test]
fn zf_08_inline_refusal_is_format_limit() {
    let mut buf = [0u8; 32];
    for v in ["10000000000000001p0", "maxv:64", "-ffffffffffffffffp-3000000000"] {
        let x = parse_real(v).expect("x");
        assert_eq!(
            x.to_inline_bytes().err(),
            Some(zenith_float::Error::InvalidArgument),
            "to_inline_bytes {v}"
        );
        assert_eq!(
            x.write_inline_bytes(&mut buf).err(),
            Some(zenith_float::Error::InvalidArgument),
            "write_inline_bytes {v}"
        );
        let y = ExactNum::from_bytes(&x.to_bytes().expect("to_bytes")).expect("from_bytes");
        assert_eq!(fmt_real(&y), fmt_real(&x), "to_bytes round trip {v}");
    }
}
/// ZF-25: gammainc_upper was Γ(s) − γ(s, x) at p_wrk + 128 bits. For large x the two cancel
/// about x / ln 2 bits (Γ(1, 1000) = e^-1000 came back wrong at 8192 bits and NaN at 256),
/// and for tiny x the γ tail vanished under Γ(s), so Γ(10, 2^-30) rounded as 9! in directed
/// modes. Now: DLMF 8.11.2 asymptotic series for x > s (remainder bounded by the first
/// omitted term, 8.11(ii)); otherwise the cancellation is measured and added to the working
/// precision, and the subtraction keeps γ's tail. References: FLINT arb gamma_upper at
/// 12000 bits (zf25/refs.py), rounded once per mode.
#[test]
fn zf_25_gammainc_upper_cancellation() {
    let mut cc = Consts::new().expect("Consts");
    let rows: &[(&str, &str, &[usize], &str)] = &[
        ("5p1", "1p-30", &[64, 256, 1024], "362879.99999999999999999999999999999999999999999999999999999999999999999999999999999999999999999995090906538858591231611078047515518694482987889134827987930995061707535518851467746173122729771633649677224519338748567780239312598821020843944191151931286816553446167684699804312434751919910533431956894065928959250863575251551210124210979225104936565368265493959689418488254723827779873666138725782483040552645620348991563294456047254853363741743495694936937242414830913856117200458057762195503845841821346210683729939822120327027365873270107809453480939686990947971592339424515495823378630842366926169064885300327474140405767607751861557594284855184546447985155570524807760126924589475383783880649842015375454475013122128704872895813450184594635657661035251240217207872593022813717819658087001095355287754695253066416897279951371660393088713041909199344992584169374132776061699659689837304030290860368426033432835860743127095950142728639126669879036231104809594858957028056237055375118186498811684960604607412998125734463758515663237263912426259965959421203933022778562650989927344702636873818009852672267818691464953955214505299550414077611779985123439936852516036620312341068867413111085259476595041947508766484399517313128542627301893800556595739949659591042518158745797311538330680390906212556223629443453595231934480023336620852018935363015399324074501871342854678856538219587097321692439655487500139046822141214656562200232448988232174417876694888649541334593834235928669053594354117272566719765167285564224221574781167345710223602309753905115678973611843102437125822630412034750618691384633166913368802686306099872979470144784897078405917223572665570847408296164715572345571850154227149344024992991989361006597909873958992289119246406324203376378795399598770025834560962768954502733554725755228532941987134302321279955593685332006342833444306940807511214360555372421619325188996777276931573608363255978520648083705558462835745388752300713216396564271764333540062845931964360915337460908661847740887294562398493926190987386845163740622446737037039831932616630718476704149460292192705495737888360049457105155018871282826881601179030572635387611877273366121025551445579547293621793737622226557463770475825816871086121143007963080944017917172488875728894638923667295431447233634155653401457272266537519159897990435533928481071411260159225751754869681377473862197861839076244005750641797047048118800373721678742740622338509586928122529665310110089205938622723659394009888874502268526497759882248593473630386395835541708269924149868763958338196820606405424553332132001036588432747084943599936638737523"),
        ("1p-30", "19p2", &[64, 256, 1024], "3.6835977775142537994474906603098708291391015800291810557384600819145819716725476840139182086049558748501137327157524898882235745592820597778241181969317236504166205436850835541673353046879775728931491493713072961572248444878545042140142240225366261880870893638537760748484618713930845576656227078964201473017024234896818291996131972002304535523299723543616897110566268209118273841108386668324177700400443731305043085514068096655523659809147872247512357408315751678016146982508945580129354278993918952160336644639762930813509616762232498978243639116176245976187359634986551531323051756120779993380725350432792730194630848164784135818387734085227104092098468626539749956691974060705342819178070681832037932312479332140097635051170152198367955001721201031698056094627175860272425718805557176911679920114685879045706943374750576931405054395528137069351425854127029412375411083414615270952381889954291837139432927731762254252785635791309069372271524039888176021512159068647434318903669986441644408783344863260142980771014795341990626617455117450802874137427209934311129169358847158685238294477258038563224983987489523810234341080862411400906033632636828728473058448608601456269625146709947714112351015431937287932705537890829018440969131893305743137289172735158999901161256039369701428402197556630907469315785718530737949386617000506786085300414750625415381255944627297509476512836291345056312533532892814205794595817031628309588087079379264091564432479183719309493318177783970786575501693530162388777963786635405768027791181781013224800283223111054392451941494610470513152294019971167544034306400979416846017959873203438933738851817122651750949335511690013826190295080743969120207186680380048116029574670354650570601207079941436327335682876871642691405293270900966143198222117481957043602648779013230273883177158080741635654233544319818090903763740700794578882661595783819987297559364705086878186304724551803925007313532807616434192393592898497199163646326085013033470066695112103944587231052989622069482380837383976132378141637281927463372787118420362455075796974055961280346060126895481147900300049656475481908997794872373310125508532093416019752021236418542830981979023612750576539452955504610977812076852706673917828186450163691940749714517977739903944958034027850773794069949941344773631210060355341024650657877172131831219115144114397402297493597641180754902324987303755073393997852791133397870793438994422066154303685383347108539887053284461952345000808429907412237633520303343844361996802346995276752941254497205025314906358660002297506392851298862294343315211043973802776779131933667675033063205e-46"),
        ("984af71cd55895c171e0628e33521e934c54871dc4906bfd7451b3f3773dd14dp-258", "907edd236e118db3cbc04f037333fb9365774f56cc015d81bc2f0c95391763c5p-248", &[64, 256, 1024], "2.5421875207686896190901313037441254300292725539448046454872191951638246782025466467579900594267466196829331528062442071184024977973712288696329505610956799663504557090965058927375104280677612326883130970682099388738546410949087223719872193157595335541939051656713033004153977566849933281505753549962634188714728478504282501637569131968488091455244027898360827280683557183445868402587641461824025833046731523694587677123516811007992777823332917715213892249634579142188288672745391817721322509008534225114918781611092685364488493738909390292921970473433991556955209415363044717107732249785369391667536493426651434639016607371070457576273004168023882364017546210081822192372258316825369565261120413885687056329957564765037608226660495590651839727344748963944844413958739492148501591484662905653975241509321018262609155846564503864261583154116201107999429135774302870280794711044207940478999703542423870792598958146040218760004760295674708804113362612759352117134913579378196434245863217540630566596677747623861305740592150875211925689477300242370653047563647097566201747497295935685511037478693609542098136021231592441808073562243386859506906001440125319895025496321968730988465912684931368910641689134329853099213883612025989376608619653833253911707750205286414588413567159527676130406435199275794650414136989414395482463559179454434434691488773663398733923617259304685906480655658271686216353937358250894847833052909333146315482459646591097878941106219132770842390305597989467253414144159910230855291421533944779632428394323222051233656128954541392503438022122108553398405074865311877013782401334407092049768054881175373182337361035501923221772565000848767230951497721902246595387248729140491540368641102433778097699004394467609773069755968978866625280436385258219769637871414451005025411351261620927500240258081596225993166110178245313493367256665001239570539233975778875404247796775830343427156779065261534911950755567766789253668359331184184832801475657964663072672590049224537172921517533191185891147133323815655367361449340123320706661315917562426358204700256492893597654111192970045472663521618841629817052303749456158750700511248451354671367481484450439859562129589237586505456504359287543215006153070341093658653568378413087828335990125905564219073334843607947832595461023980490112628879263953446086932579746207357850279966902362610888943987199777607975142776903701244287009766203512731721658033978791359156237909944494815926322010443447784699931021943382436408123552143600438897417138291658855440612191724197056209064116718915692042531216165544509633490328581008814469458543914582658992772257e-65"),
        ("b8a485fcec7076dd9113e25ab0c53a5c5e663ea0f043076b684bed20554e6fc7p-252", "e0fb1710738037213961b9edfdba6fa0885f6240342c95b1f82b9c10fa0e080dp-248", &[64, 256, 1024], "1.2729660415978606412290262044481256235517504846507036172820644472952024285496748723707921058009874594200200749459687878038375488376168449442670864413233360427291592408153050543457242124463949603746790109242863334387189080754086114355218849429723713431606995066294907328228094856197352638033071452971632108412160031727280669706960142748755442561836369046782927133896356774322302812532181846938349907563422534784386958155173454667547659071501964345610397929149242607521168326102892919398880631107115179325836791917529640153639281527039222128990234955690749606771041689958686911754140080244659932065994291199360038525207403011064025246262772067124070773131744517939282565261719278336217420734878555674263603549104078281866925595450069414347519528231085974178227401994223316625161937245651697240882948682835545994147593890504126507414814097765030388396455919090147381163577842857495800611400855346831152974492468549217803524548690348939350032758458381332359898908901804460020401306985062725378126279444820198724164341838010735208003697424336288182262308473467070885494586651530544930380526345926656910641103155273542460113255338146859883494857448771858892791782855881298173871343018081695185229160725722913702340470537344575936067640273963045036349443982120374066504671924123401278404703590391667763886488517558659121042086744521892081622052063614508033745423058944584494276076294569199912571524463771543478886712117674304804125852025543197204901814366170943554054214802500714686547655336892311375546198677148876856125085087994244693919703106037565099621773965078512866064138900048490259519986864175681419328072379122192563288407762440016097051835268235746374525128540864401222017701092531555069946515155084533378146483712908531937016555084684374258914653068530095343188549868268735518083101505798839461367726905156077266898431578545818340183026560980671739292582654717085397012089367982252522866659090722964612275499678721512660298547601556687581526200414317720422451257381827087680381253830436997888266038473675259338055635294338560228123381935948790092255115045935009942508110182618302805497850744116036946088946793159040104303498149246600290925598397637681344362554307521512245620420285387221434585203358207095890215920966104213540348425872461274015048369077198171294687873382169803520244004295398931494578525272321122740467588007273934272771557281644859666809982202272683194446660177407241028427201104927188518681839919305382383130675707803821550004596756374350635777678286184343668450026095023096341390814333803887521440314595610189628169637365901341245455293897577714768226174612318389820583171346e-73"),
        ("8437800ebff961a7b75b49a00d002cedf7e32056a6cde8cadf058e7893c172b9p-258", "d90119b2ec65629925f3cab074dee5af1dae506c63bc3d7d22c19c04a6bdd8b1p-248", &[64, 256, 1024], "5.2445815034168372519673221221128407871606479684387776585030863998351916545213820919366591750516149377002137983909147065416411836917041346501652328244919327611229795445648719313823390227083263842182388327328953448367871985029788894823789867916562241042222394009863282355229157798145789102621806224233038568056554488812179752236141923262697378940474480544712965251378700268554594476916308139169038732700031409924511596992284570993488905518083398381090738949788856976090415959618969511614011423157540642975497881052182921504867899896926284484845746463832296415242963387094190377116356047621460581717984017871337733447518161972814908980215829884008753549965611681055060259871531908152227818649665412262769011901059259499863618448866898207327056273865112490591997271489710116648118997257509030147418134879948635542277154323974680899162457329348816422279014787872663747591137543277595189278313607713322141010254153489208225182472750214823386222833780199131897420467209066120525506948456993397607116085652363036664514338792005750480727151232036057961628428480518406175743341280363464728789172708527438898655378498285693198560719630512753979922377366987032021329410358073639656228280965021877833145440753470358491403966808859040628184253862630068996186315622100444196266665897544612092194058236616815679965139784846551039181893736812813402104982827748087620097601507550794090135845354541731419051547752438190821392919523713262566220903995955594611286291166206744639059724585387078193444935620310134896708854456013696589337607126687717400825441826963826685761920187395400879187894761815490247113150805967151605209923664218819234762005237640706205895053149400326777914638788003528474454808716109698996278528360183710072702751106368116316167009797535513640660739313400745733657009091550241619686755647930139913686537924633357569017968513474308604821583103140434768707832449097049851122015728961864501819161233064293959551804571054729492844137062168799735469178711235642231646393189460072040305568779076274274663511457850100055420872079796086996294111055271955715675000499196429682620256991855635250448052414521915448885906453214611278208136988007278263308033384583471643336809403908537300311540158445034983879507230784452494935127846645950824952879519852000583605292482696231129237989008179094522309756585843202977269386238020764256235327015240621678189696026989582638662930098153515196564574852751627828387627227789853556993344228061215635537607241459816155522704081921289596813051119298386413725416523749447716864643156845670819566031782191411407325075825054843334604761925800647513500901592235329671621638058e-97"),
        ("1p-1", "7dp3", &[64, 256, 2048], "1.6043577638642438113772375618601748734459186946199968089228844478245098891209212029490092485168108886560967084820632052152597184133126014658524885626641566100669701644009857491949667891377562879446044711673944164177291216111632772564487467431681089544732146562006820592011473910128610555149464681829806992316898615625660730426974193216638755520992806511445135344553566713827727906513341797413216560845408028495704413956286454332480420429600739770439720884068187985704097319023430872791108080168158901880714570041416739100612883488474134051598784131633407636938798724183432537461244818802730551488148174031877326831158142410021038676370824076527083899033057129419234805882770127324193557214199907364057616262147514547680014232300512339299421412923471791616804202619189639062763513578592315799541887320675114386700090065178466379288005597204981935237627549260617777579487912962734463036476253640976285373346448676148321349731874077165409654206697344486029954137076876651283557137020279785851034224005330346898260785470518837578850268103452242136907378412705734924704568992532657116518510063142056512096721616902860146522271303459238933253033643052771197785337223582271384455330469614531082611932398826862198569803648132134094255697943118258639712032956371058843673962661359674520296429757843216226371360194266980261479284608475649870652914554453559704038213644609310691157650478859148349504568879555142123579933660302776177275697258125506833217134393028319938101651088565687243330341576951867239148357092415006870420675228155377653877536738724450631469990249947570825086114301398502139343986061466679797533746692569781797310431549925250910469253805973684524717684804895029059907287170929446131268565496280628358270980038374860559807258119220582951404612039913631727688588618578387381181353635199667649461131893574471328589358552666954191849154692084433968209959272569435470482209667336736344226072691143464981199230505097327501363264366531043899446326019743510249510632150779001606080941987847014015118427671928317635886053664075633493562717143783502642365481694316545135189260652336351036012862065912497078324662233392513938219512431959806896199158195237452143047381306012963621473815731123656861610539877736294635835907136919799938390831624176468290790523995206101854905873109720112583210771457822106332563083806686134406299785896735688915528578193866034173892327304406333244436571258738979079694865141155847080000828683344534495908840927727173806536949512851830227030702957282487818055792359439779247627545392045966988567066407654149762433453290541059081833086282959698180580822978649855126949755656e-436"),
        ("1p0", "7dp3", &[64, 256, 2048, 8192], "5.0759588975494567652918094795743369193055992828928373618323938454105405429748191756796621690465428678636671068310652851135787934480190632251259072300213915638091771495398351108574919194309548129952421441572726108465407163812260104924530270737073247546217081943180823516857873407345613076984468096760005536701904004361380296144254899617340297251706669750000570059682610377163358310506430299251123729911069692020702293072960714387795991606507061281606125451899980021585749973908979115477589463559076631105369581825728355661271532526542763867507763553707531903008307478183606488322891882266274479887514197253430346787661973689973287567728837484932921259627218043915405648417551234455887469753663861713150289624405538369220597632574027132836904060897683899611864713873289750189631568407214092711576946752759134734941753791122188574460653736449732714603603391674234175451495161445510453833937601031893071722128023689749004145108809616839713923899063996531729092040957023295770228890632386537336856779368491350040944508733858834307744320891367716615473224487092211285230401570273158036432389262561609431049525947493603822742518180342416479647608788319009644515435918586632415033154675317637001388598734563555166743379207523651803253831502425243786364160427340619439842059064589941613609433628799760431856870867286359745929670601163974827607249720038272692134287889668702688646515670694713634033192190713908885436378150421648976478211476073745223029805246307659253903388413453702697274464813154043347106379046257801619255743856774047065432980471288369014612858345597942172442318488198755111624983720000619234727175528450222737814594755762810312409738817292117920037345974730550941545930409544544901283795220817378832289072305758150252302895105283576862167579478520987016483151903093952788881599722031632601144744825408941122761271061091031726888113182659535376558156733021064664346795960835662218826497426525238061618560695138216043231697088137582206181366929130833278180819001019872211244402818487628942926029378423954660496267279194963584321339049070294393498028871859360435529919249539043852789939758166402995106110659892372679074768188333469135158365234934198104193212215931359022540551892819541474512444933997283147648229908651737467542607358173088993248900057435067112619999683865756281646693778091690105158108448814238425634139569331148633806903280448711845025638089128288904186047456783347482633406246864246464624085661580189455905279014573403282491143708975907720325993518978228113910525629685472989415387676795184973489736241980343510538418118825909048074093473416299502519486880264909711158913049e-435"),
        ("1p1", "7dp3", &[64, 256, 2048], "5.0810348564470062220571012890539112562249048821757301991942262392559510835177939948553418312155894107315307739378963503986923722414670822883510331372514129553729863266893749459683494113503857678082373863014298834573872570976072365029454801007810320793763299025124004340374731280752958690061452564856765542238605908365741676440399154516957637548958376419750570629742292987540521668816936729550374853640980761712722995366033675102183787598113568342887731577351880001607335723882888094593067053022635707736474951407554084016932804059069306631375271317261239434911315785661790094811214774148540754367401711450683777134449635663663260855296566322417854180886845261959321054065968785690343357223417525574863439914029943907589818230206601159969740964958581583511476578587163039939821199975621306804288523699511893869676695544913310763035114390186182447318206995065908409626946656606955964287771538632924964793850151713438753149253918426456553637822963060528260821132997980319065999119523018923874193636147859841390985453242592693142052065212259084332088697711579303496515631971843431194468821651824171040480575473441097426565260698522758896127256397107328654159951354505219047448187829992954638389987333298118721910122586731175455057085333927669030150524587767960059281901123654531555223043062428560192288727738153646105675600271765138802434856969758310964826422177558371391335162186365408347667225382904622794321814528572070625454689687549818968252835051553966913157291801867156399971739277967197390453485425304059420874999600630821112498413451759657383627471203943540114614760806686953866736608703720619853961902703978672960552409350518573122722148556109410037957383320705281492487476339954089446185079016038196211121361378063908402555198000388860439029747057999508003499635054997046741670481321753664233745889570234350063884032332152122758615001295842194911934714889754085729011142756796497881045323923951763299680179255833354259274928785225719788387548296059964111458999820020892083455647221306116571868955407802378615156763546474158547905660388119364687891526900731219795965449168788582896642729697924569398101216770552265051753842956521802604293523600169132302297405428147290381563092444712361015986957378931280430795878138560389205010149965531262082242148957492502179732619683549622037928340471869781795263266557263052664059773708900479782440710183729160556870663727217417193090233504240130830116039653111110711088709747241769645361184293587976685773634852684883628046319512497206342024436155315158462404803064471980158463225978222323854048956536944734957122167566889715802022006367145174620870071962e-432"),
        ("5p-1", "7dp3", &[64, 256, 2048], "1.6075680845475271767933617526850629148606660908185346794438181898896061248281308884210104269183792500032067028920091762089635784623531877272708885529907227822028884110618975698127849617369290133098770539366628460263684438893824256656456749438284488619022027082857502690268977968865698340694234336175042147434572137957626817109627592207001447077774237828098490757497212079152119330508384171129895176678028434357022296182555889000358352034882477132022192579952466601592240272513725920833738241976890975983742533827525206565214340666709154156649088846316165269540475245480006720941396625783002386313205258551703018107799614504617051209163406133822144429656943291303223231854964195998717476055969206298649960071341182320601685188829512356559674698025146156183866585819888732718453958932967727490388620274626653581929815209454784609569870430442805380167980328003214601470812246579858250292680246997423749552338924937726294450327691301134680790679225538711384388229161420740560703972280929268725515711139123914094200096756335922489915542993993733214733875326115116802487226775829732032692628672686378381937689190268675709708964323617356637914531295210562790458298048072594649080005653312947921773082593655744873881190792781799059847004953465180547460299690143299601405465726778509637004911294576396188333341961167166715707885432995182101080888901485073522019819406622496888830722992225048747987944842894016447881792317550349343531242332505297597931051306760101879882075160413423151478742607900205971480141999342363994673501315504320578304667114123856345290327902093400725176140125950300325028265797346364679886637535489397871175741590059996686702935735394097094608065908199020814804054591114868555937292724619480060982390012132760611714946622579533448253770145403730981619604752548988732829990828666491213358723517570012758324005280317851067969179659802035467180599985733509516879289267502490355682446944997752875838988850762585328629584525410202736069072938021725145280163772114980586951237591036501954285994991045006520400619892620433262515552525536859084365482530121424297416443572919480789108894321080461354079482349878593342656077674711891899346393041633895990327617584629387668941278269044216902731660058318285828611443705577378262250492653159340980454579170379006463259825249969864868946415665688082372987767790684170518299752160258006817275140536812929191418516016227763599910206120139235856750276268290498834454941066961511969100900194556267520773272262273285746803381750083365333290530051869411062786700723663018240550642013101037450913542691390184655493941339920114587682249725969331195815223907e-430"),
        ("5p1", "7dp3", &[64, 256, 2048], "5.1220105703780659235331636180047469959106233734457662111057042995186513658901129279976398235873146215502553891231859550860989003037838005978973177931202987422707161075583475715864461003734082259285798018515675579126531879768057578539524467538111276854762233600933756590481662194990267168651471384657336114340805629955542824376494562025978242589565100204255807484983748688102957540085959590542158988386760145999483608572328506806040522982204109476637013708952614534659606998405483522505100233892523127199561036217638653123107178285635964022392070592832645052913918188715105584424197081406495360816107886590218726373528229807242690455858185029581906332677590673071977769518585976324057320784745207892658219053654656561132343325248286419720298515054934636752493550487959672347251709737921375631662950050494696986507438520980394208995842266244170384144063860650560221503972917827350412292740265603998568287441465264112564890375036691609479190711710594525008515224802097353870889878718324519001395258768431855195303451350816485897502418721301719958356372572045222378203357591137729543765737872460541411539552476272338863070409771477870126825169962224783736415724267546222537326691271827536477953927665863141039012450367245720996802679497876560139203100919045460597008863116589418820622894385207036586132707264755679486295815313557239138081029080637363546599500304866543989206993598242959791157801523777592104945158212133716438343836962960052925027725581487848614222346672972455482514272591650516428474431280038424878374550743041916772327194288371263541675622105574265208867684866885360150172659627163820607278809056754651970508333349589316065850348824489818802844271994194238647443002225089567529292527534563429627603443973711728534089274751674489244805812109122840173044441570751159535705560979170804718711962498804565999541358591048121330438442601528870303895124050500339145731359059813997170790261420673223610452705292006076925455673907743835158402661430008210876709708606484085794167143624279190113094659206550605059760631434281472388322813667537323914188181286319038050574531516896695327584353803779422784080458878313743091275670358720642059035383793704656499132748583740540610828225441558349440321678020114220096493188298589995587272424047732709772741289492203907144251505352355747892540051636967522650017187313542569008119302693518708855610838672462218665598787367427076986134628187251425075357615506638360514148990502805588447821108661332049500472771484210672767143312682310178250765744265991678147453146430325936737656331911167932761815247688176507652278274003811852355327521280281317779063697181e-408"),
        ("3p-1", "1p-200", &[64, 256, 1024], "0.88622692545275801364908374167057259139877472806119356410690389492645564229551609068747532804199943502159830827242315500069970023451451046796590918069824125254520851873277754326083219540033683741185083641406449171810554619734842339356440174149313924062256934307603757279354214319804292511944689199806605378443508979217286245430421384143374803370983998246259485007007763830485415667148508645746216060687778990245862429707626249694852358548300326180968330085426816401036168966494395374968021871716499353773151847560121747945635436871985298301568719709407233703503666691295184409245172429831035487461066109563143821284744099548256165872332172634185756383048085112941613062658312661726895370098060691436375293448745336517608460044605414286172284005243051439460178202283732321225491613014967510797246311816172607135346517020505427038834149179019383697324187305250960247105079724667929203183404020843644313170509038987985705224216461478467450700085265301737496840723374700341907380010935085298642104777789631659986360006019743931822288394744550307713179793380362550713616821360837281165514817778844393032533589598146630723832213457754956351181058401035415028558979015486964149074086003694240718331060224714105201119452732595973703534918390059821525121737297537209685886113717722533709597792132307009581274490524287161863606888728490611517375214527454158700039111423767233749726870167391096372843670567392710777925122493777336071104737779592476961360226321606994887685929267701750509619141344883249642914183779984122286653061940571717933747591606973921646877378946233755913634179183301135806625905959824535365498719994627718266686178411513154754656972809334995986723350064387921455635736341862338541390993859436602432973341573441633372396223320563218393237185469477444263292686842881175867349340764652811404026690122240557447383734374323870609068200359258912292651133587011320273953873792848456191496082217971044128826242009418021322649919013378891111896959732233965096978615934450129697337915078207085668461100220320263207453962794387452598861204423174440597083185844664117660948406819165519839803551212003390573267889938172173037280513865361506221710199812372732607245437260861688471690529406659405049512341400035600730559654259446269955639941115390061747850079158800046611407733733570649495652975564798659948526397534376736605102582973461681198446564536092909898177009428254961118100502177128707981507397384404653595894706376918859240433917124740048714273512769993851313058357047728147753395620454376182133001434732775864904978070382372878386160158222349618854963464539597202965146133352452885375165108780"),
        ("1p2", "3p0", &[64, 256, 1024], "3.8833913326933875523887084207048185772725681906970108142749627532727320829555649120162310504574008730593036347697071921674421841654215951756320969648740303388774000529615695548763290324528857847609809263824729023279974414382250401904070600628356425168202066865781211759649437742071537418739421553292894388511688840233639308116361989275972725014247684632371668191586673551869497118140604800214387081867257763643917263576722572629057573964762536543806742562710829073670298607945686300538518242314460980023787095026619693877290954498301340350639109386361109726684845949265414459352554393841060753034594139575901400702651560626534973055775871952118611178338374775799295740085943769161226988659028191145041163252774565702673840060791020000440724441266398715807854566065013267471431416807411453219029166641211151940676490970606075333248420449989977972081036499283613740651929552133595905165572684209997872228624363053222107574898980788950010132354003706371454550795863587542714670826280184624488374640380608367856333272036200143597809828200578870350409834212135039429437852503165478055568410745027698789585125446577335428959741079707690978181413410383186081541859531689237082091233366167411438915762044412896667459355417287346602926757023588870834910395885073096028531084729301838486261828329231526731219857250289093192387954570565714233569338181378049946554990776081722586347660480008247499737071561756710387949099413578428691411694651139415127030872657498020587824011230951086052815432857714438250037610449998565857429058223060797219916135125696865106851054163786205054242861700077806381227690899694792423138237254462168158777686535373642233767450150181935336522226665072452691595704348201041354061556778316697915990823405621866221647893758320017917883497792693416363430566665791010975538271955253305877021705022817112647762960065923701121912101035918360145624797337361947522181631347758051097530590395779532305536541869772142174854180783223142584492454075155155167094773341340256581685639605557418786831190566848047192081435437768785063275155314939950102528989754307797613851182495904117871966782349338638241236573151172931180272571405869802289695373499892046782004042574716546727472734658118450094792407521748828316010074537253116910481189541546591986691687231297303599795906491975733234711684470529260445257575744259104746396171113714477945750089394381448110076668031975942783588003964386484832615321885750376931260314008996736061506360455106080682917784048485180685848640284647695506051846035391350610312858429102059425754200315738131895328328218289352489574994857860284329277032585512219836209608663"),
        ("7p-2", "5p0", &[64, 256, 1024], "0.025770386159647409877938681158143868769059170589517451630307753281635538584340554817469302368841375585299115090698159426228735000741008345078133218156261702864131853017324780014304617227332949417532631944803075988332452307065508532407601395255427487312344512391363688268961471250877149021299185522800149574699810367928054816735023330484205081119806671937486834433019702485580445049719159410205921297024248564173129330101453663117313238136171510471401844355593000203299213475322893777010408858305285980514599373617186072141868651629031491433093679743710064709065444370654576869944391596880827255186560171085160166837547691171451214014014562098013479181291145609371562999127919401309673538033897747469779350224244542287345497013350140779355774005823778828830333531535887695032898736592214444165801970646929689577866299366447935894377749951391095650506854393982903846734809242867848026848749986691309926971594128846988459954613946283305786276479294644807746455298183734194050815919580321154249267655777323234938208781989408757826944118249341130936817768011787624564409073387203199010614327509736506324291707944774926797921788830929752962925407588490259713978215760470624319972749572518768170812496333850488688920744435052425740226768126880057562929141796633305491246383181111853414534331438908687713000652241131405675754511002370018573129318509354408856708385757861626714502179304061536709202237982208999417438818816747595685576371661370440191572005250653913368687028963185091494703849199976421573645062357828434932905831431125353411669659464833602355527410127632780842797483700288907328841402246206422671939193890328621531717824631452504274170165028535276554384454922873582629059041915455252950567703567389408684935313682499744963851484712116802705505252516538872420019837650855551464252318332735597788868499560057079531562717282706712320906383678917750356985867097190850918492371598126003056330010789570456158849613595311167705939075215825184194063833058468594631215805757026815869234758084964241921158378517939304815761028513815457609836559283073291918805686985201576087916382103019687370495497415163113798256388009394269115364936851920930735484750530137805416402566584419047908695486762392607440520652870981029430445280939651245338911121779846333652874553022934440459486249233773874204924197350988899900587086375130888625565122261562176938983294360545305942596287453311568957854874641259371211060763673504050520889883849439006467860236028202609269151711793583011352356159982304340467249931406026325668420511056186269270613399520007178681252026664159469782479238439350206876877870226676984382385290674"),
        ("1p6", "1p6", &[64, 256, 1024], "958345301841954115900523840912859299205575359540812255045779173831750733742683039141825.12942877382029299716772735032342250268910872007757418283103825204316887544113561082218125811304024513020308615947702142900535505495372053454724576241002456219495301865343273035977202740145355013845176678142479809295539214256896982193690715904287250488999759134491128007267704763230273946959950596940734420873496409883339331658995006730433486827202022833421645036120738380699481844348009103023312561320575664262581691975230068304084650851474288422718318176897662067341009868351625722669801271934657974109786101783273894251285776050933575282703748799210493242499248241223531932704849966102463752739714626420447150463591206554152119795654474770058704547390964781819213480267532248053805235479342498315991082855455526832471158115063431997117610478631425634285176922511260367442071155964019438217943621867872357448126201151146552482488865866087342101705954304486195323825696491493790016941029726718185979295941240875645522804857778531472530297503380013980305954438517401963137405379448399092323461909588737353433000340230212430959444610288858403517299628483599113025047021287027679268140042330708044878567220804524680503970272229577173022912399462372374305570320564294728733534015978697407822339132412126479313785246074325934135027722250055891015152041962844418027151578756557640132698122845250926828700824274460373493470075511142698907022674635544798490169747426814528283063658139975769566195622437337094094697809469405936916434560032761287671168426183482650829042289819808425630595362008891526508771627930037010274697347146315559239569844837615887582494159693809274653872752859423746836537832653169849398524842319748665773253014244985897215249179206360265475364307144505318427290932723759258953724919485007385713381262192883758196657954772943627489969616847652902607010333523948478982113266558102350078915713849842085981700282129227323389038559923853026461134944434020429512433109865292065824556023403605190537757470141557291914131086309337179327778671778385388389875250386430874085212450304954803088084150720734421899551778234835737769839246301381604749191476607842322212598437142771092730182605785388277066760865223450277799968458377640940093561331343138927308648891233985821996363196407965331329950268045080897531792462133229604348835461713578458585037337159637805035232452914414297721952429833796939714962785843733349453085830454982960778481485416072093300093006915901941678326086535839452802330326299427500574632158057468070805316084591811348339262423001473136069430964129824715035187003074134814186242015016736"),
        ("3p4", "5p3", &[64, 256, 1024], "227696373622022100770879535610832541596504880359268138671883.56959298738493687869854082060524683787055618171724703793002471941809669648296671543118944363038627296178718936891268341683831427881626259795950539485078486083832222333890465572782428646543496939483838041011299969737680249290476343778778860125053076245943551547285970057446802806728492092009700263161441612208609723285877608304565867209601108366283331222262004806696354644722833045682495147302310312912565654004974437794260875085662858956479436470330095242832251822917730473042215912534451017785421753209801216907359836608579326458765823665953522236663223649223449839663017567354002981883309791678578350413935000483798932941085766412120841959518927976842998014020695533733565827825503020783651435998067915964937551417917471094488361935939309630560747184497082185569339456650939041638772047110549672083499304697862650531735709154108182979557698305505022726776092383234215491020429013241884513468126026145196153814952046069607725432176660840484628021570077109663500865416541361249866932237910276900734788504064588871529258238073184507598853161325587547094216317161151333380002473797362913304992601323493146836486569535459952808378466543956193404196996857560428197525335606538070081531485501219446452970695730389658906494708320212066225527131276127960425299085504426818612480091026115437457722891486137259308409475334471813331493718727672877500552795310453160197233514413857437278581626493260237332545416096359439439449621202065653071609421338810524707487842729016718681923440917009061999554523713697221678375822136911974789192178179019403765109197897687900858672449697130573837949396556243548305391258704066718415444885622097442561812735890629460793899250310834501930966221597818103080763453993390904589863396884653943396134440358452614666119597725657692462888747242097818634778984635551971599468953475096011896178396326972160731009617279907161434862909858232771362597318019405934232563130496365751820936990867403473681413109097776928864375054171501432168380268348820398255571898328676863853894042781488847632640037411186986501553980122849933972435909796835315439210312093233410083700645887632440246751091334537073822718307207639457935300055071714180100752972814939989480450902893419863549307790477783341836629034371350662871486356514026365991477063411531310498946243112000427108017854436525062390930838129864754711431548275554901132836993637552324050644828329713570465389959134672403877036174936656041989740318959968956908472475128950662812705176167699727368553065047468730797826495510150542546887257135813069802367969594446070477735382767473"),
        ("1p0", "1p40", &[64, 256, 1024], "2.6363825696815466856704062519983232570352956115998952335329623103316346705123875977287785721662818794227723547845573747178692481956897982957258574938126676629789970872907103794408580725312240311533514969745677005164527174343907255412437406848318163332733698570030367265303997607164199867835067884329676326006064000289339165429640727825732734102774744855342076553567656663512642815359735061379726592353794920150283215754552704567043385353987662817935491447879484307975878823076551917813250060622980737737315670846968567914700958809658680216109131484376483910816134372766984940121292888750954695283153650723019037442446648104545077945411711075004000700718079859523177044526014916634623770967076192176535801782118616762028630783172876081893175516572809876413628992928560907547850794553258023310623900094187568784370492189464829225535066774516483603961086351193241281686598157807259499028007747209062161257641491920907211954585847454501028771545985943792512723940762750387689537378180306980069076408154684203009351390497351233271717807151531662470841649926317783428474217339734613704401653749229388162489292422586858234397300164942312313891023608188711625902859267158875602026394526091907650682448063130363558249899475346222433261666924231728137468779573646899999937690460052500681518730610005709340293738654617805166328149857383034245609604954984803387063868885703884813901472218272490022134592013779425717709637740958967277215378559714648093721742855431277333971608686824032981452294442927711296682555482245371344571295685735225810746944369908162817599812380875446409902159244949663638973840235370951533965185283905376294680188708350079276738108509303101963366155436597110749680077483624175129350702119359350053963381174342772311104747283097073976182993358422072058722395045416121220815048882742842265544840801644562031869095240162006863544827317639066082565777984787604749039381597710862181253750857073485126590083411590374096551361827455612499129214892995216408710859509330018943298405752528388664637687262145063399615325923196157977915099978102745174898658217085412626007727210300256216049189877020261391911412999702126828102392461170578734880045249362200328388173824217682680736233460894770135332578878106432662386729128220527312030545508967362100243447085235077559302184305912559145242354597911469632820056330424518552013241577401024341036838614771379006884734160893072585979882232117674078076125085107590170641127641771230265089653945362834004956746582873572195813673585313265168970376423643542068624020458318268245858937655637014581747999633518141151539816554204901646349258390449124408861091232e-477511832732"),
        ("3p-1", "1p20", &[64, 256, 1024], "1.7363768200289580487839649034974771823799599009785300349089651212144937596693775371016806760720384481648061499836242031474786196963932589589089307211846662685598406055942968962185154220173321303242047772274233036807964668094569627996553089842978225583044157675662382281214572152697941323531457516053933573492579061917383903413883332421035542531424416361917811089754329593620803669564548878562402117441644425281846507593062302769159995417377315391102825749946702442893944491538207435265629009373952013785759838870216826693837837931320874492641988010367170390374403753651994056344639116464980805128973356526176058752769944122416465520846980518670997197702720103681267379951170975307856823302826051926239321642944305823020160620844724251576162822436395014518064856894903274188249790231822653087169508463904855506766275770149349169603604030500852222750788628035393687373488053223349108333681533143961808205272903285932055604862587147683305105050872888611768006166399999189128769858440802426376620069378485749735544785933586703445119792353402096832899514959926430778213713041460743549476280181575068414573789842293409575602231328608304301222115901523641501227638825579403340647167064341028750766565446675988958799013830958432329406689380308398942955771377631827142485394239580561231723723089331872448872895455488832603549884503492704318440753985735954533136208876946406510270317856896566350823416345521959335444612664822120050016874591201115249320287800722315364441533664828125056969600977981053746741881217364659818706107610902239365330849451541298832330718498342827609208559598056579322635778936902359343879116454030803575487927103216347039038436794216542544712452141654498779848565269915709832837184460516701144378764601001386497947096063619644104270928606419118862081935274362646094389808494693520534144370401049510108521453037431149253869982580561481884512436888213469775778179115785966110482043686577387871171601517306373707583562974114628125063255532469825695046744514823590597847219387760938546762369412813640359403952643174235629129405833212031534718734172890336393519585495243102124904240086045051055616435327975628881925096919660452035741722907498429575340001082853811589033618555796419734653417815136865700383095562334740343417182400769091975323978064570579109107694745522157556841202969572518948404474451708235673194083084392017595423368790723576224720877609628386551106340486722899948408701978432717487522107375996778060330208276579817137948450845372360013906987999766024034040656212431349376188617586207667049665060121921097799446131726949156278599067800801864267146950139643325120752470115e-455388"),
    ];
    for &(ss, xs, ps, want) in rows {
        let s = parse_real(ss).expect("s");
        let x = parse_real(xs).expect("x");
        for &p in ps {
            for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
                let g = s.gammainc_upper(&x, p, rm, &mut cc);
                let w = cr_ref(want, p, rm, &mut cc);
                assert_eq!(
                    g.compare(&w),
                    Some(0),
                    "gammainc_upper({ss}, {xs}) p={p} {rm:?}: {g} {:?} want {w}",
                    g.err()
                );
            }
        }
    }
}

/// ZF-24: The alternating Si/Ci series loses about 1.44x bits to cancellation.
/// References: FLINT arb si/ci at 12000 bits (zf24/refs.py).
#[test]
fn zf_24_large_arguments() {
    let mut cc = Consts::new().expect("Consts");
    let rows: &[(&str, &str, &str)] = &[
        ("si", "19p2", "1.5622254668890562933523451388045026772278249805410834563843116927091608709090404901824073524559347599047987121847137155269263300803229352386782035600882314026615094015126096075749662470221680410806698996367752493045787896757238025307932036460114400067197089476118776887283479132400264236909559625340295182783984900186081824209816640158649544156388536683390468652136408215610293606820070200379026305505422991512108303976491413283031016002078972869128072010482352146284013057857406479121029811987683525713459187848994880790382193523857452562998166628290146768126995055944457901378112581574255280139034633924200577378515427030304501000085530910796546955145345181864363628631530927593345335610278328889754724960535926211236543606399342124662879830134798201931584235596173481190313831295240621783444739326688507105900081934701104859671716338230948033837331860999797306695804611466445199798782858073150283623231854555787288940303089712544403868024629070127608943410075231750654308356122535186668630503054966803627745194764108938057366993382411551997917972193547730010937883378694818017022716972450242348535454927152756451678217466112769814110058982426070949671225158329693227236393045547547257655577523741638241919455813025550586427647947895369690372042438352334966914892967353433061137936720419679709492642870393554570239619885074211171383292254398502308693309568671013722946156384954965739269143333702160912975632421637054723806861313454925036960500668745256782165136715997702967920873471056173419437645123275060844528753932013445179958532740881672290821624889954114902345315689974515244239396268275183753211657649256579405865249395994684727836506822439445161206486290981494304687520830249285230840045949428507407427435283217288740103618508096574460342503488179635342881967760131815859783431169727890127043761387686098744501164750542414910038945074819964373282438531173796465870717710035290797945988314704748045232063275478407898718184636688799279772784113257773278517646933807534747568941920826132557979336253466225210403756965628590272468387234794320546758057028582677281512287242209953399687538312447987838400476192771559413641163404858821966813506831673265801721031967395485045884207503716537120434751414692957892138337958528888929674896056905534614613529774531525866439788800058767261690087932944875596893346229677291126810404526627401223309023150662250594230744982685683331818523497103161582303760563434740828736434190494709597094535443984449783233396187593642810354727098608781494595467284479685504891896694182725348912860298439976113069934172748971025279703007949528617290442950157632159846239830"),
        ("si", "c9p-1", "1.5608558784061177878892908807431710731394812631057846315500829098922183898908830941597962382426775202868333996951804718938212962206840708740355238066246420114408235116642254041244883253298455970630761537753369744668284858610814689455367163529418705281137393540622565100524678829395817662448928132431193434637575790833913850244132427809727726937621944585524120246186927979973549429568905226301945950733504695116172953972495735752727148404060406314452505423502706751248347459235443031387746056745005919127773573279447271126330778896835715912021635909356086865945547494005035824172492348298320781261157637180070330232804264521422659539106622930609671347275252774658113081662617558988497874709882361240139096813665105609837274737185041529914050801063602614116222744054481589220001438574331648329872582657585008407625564870930630907636933195473111139523923264614948809187858140797219431109510232149338851207594467512266010258724751790788501252400778493342493786294420594570020880796126955118256864669143639518949119565322204213024912568820510449104518342621452582503587874976596450974662827518334658243169336023908702306387071345117901554298975348745967474659293428979888856810659313258227476793518800226475428927689075256507739154489943762203017320274538953168274928317945384910569707311864248428598224381400629150856510215724706587572331521507472645151282040830070756266739177702225412956667396222414689149697745297648795771016829927917791754111540299843866835835616664505061709435719789227669297716462095618436773555945409301580228256027173535478687489381559917825054350482960687388118321138662479156448773853284438413642365000338583004390279525624703009727590943710914979580389172671475257755118772468253108303721845054059285251765196377838780274219380502728156401406580867829760937763874688006037748495000752738501413967562917363568028675078263038263828184037921476849526739601318377057141356143068919374642932149676220435497941618088986857706801282591552813350995969187373864686102980003001354975520687229496035344202471345792861231070320401571693470703740674991878846540602958497740198470372469198605034182231976228627154347476278651262446999683763374530866947088377350048384207383568893512883505500658700690134258893902375934803837258838612713550164770459613871806316611256905146476094556153892582095007856347214384709298828175388290555561699626188652000011564739758397045031630691484462892197490325151159663979168046568147225200407151667792931005736480263432189452742360411775210267008001701010434546055927301072215207507210089726405048239830823678867249783116786606394856776803251716918279268136"),
        ("si", "1p9", "1.5727429488260625275912304359551161714027668697937757230580000922943694514500253292114401348973046116886010576099870672751232309443163360458425278408646407070042584614547669447538573154942998947254351464206869947634077147420865087497532263223034294810197992633666392498703203344905030803624938050140836042299403835585576580760060762592209594857153027667244554352359676331682127613149932736237736776407876803518392298684036327553448225659610939842914706968468071779498280239376325493124648121397610618384668264460099964791688203541816202405415866589064622177632459052128293516036776936601275214025239032490654969509654816955822810840883321629430666162868697487265462173394782106908658271906982691313767813248327884502015660084981770784837557416498400512927105544586192297561184634657271370829597067011625545874856842414235839147663518529318088326837402919137909973720216945802592312910574716267133684707032782303894573818400339727706372849215932187234574644457792441967796839996816907649088148671763990784429828391866466082144279431415078606061426268681477526856138017277012096901982046379724867317300134048238910675211242870252305521830571140813795948314095313246164177154602943118374599388826257794875024620399521785273522558317642548702522324603324326562904708348589451895399205305995828023003159437097387167910697534035244115802065976198755504635967847296420701587842605207427766148663684850794507816213837671685229785377196299512369468404290685375139614198461837135476022259641798823658542866315556686472700257049128358634728120027276799643141932589648279349943863819645407610914504040735833892304088650259327204227254722784722251868030619996299248548367359562704009139149236254935068599710661293977409712701294376969650849115935314635403205286257060744632734704380643662872723413920392057299926116586108376572786929999816048744846849836946996962918947140694467965688645579415566223960604370353188092007251223007183525299205136580143096888676418515883133882982141486101860945309343085444014112754464165554624588977902779991502558839453643028849462367278953092232540431234517436401952575117526263349363474192497511779455918743296574948741418026219043528944529004475085966986610487481993467891937963219765349870815990477817306899845405586902759927769025830284739149787404288613679628213805955459816244349841536515289110594259510393712570720432987803493577923526885499735239458488959681283208862685713391787402039853967580707123840079948509948090429882429821563799615666173427383815173546080642334415746658793559212096642999229826393409630584529766984709560670404054184819620867324983046664842135240"),
        ("si", "3p6", "1.5756749066802059484874551074969763380125661891822897526503220841905202814374907886061199634375692134371860381855336907916178011434947991100760698889612093463772861512328899313948773234140593894043754103096652507164500431988068180638089289646388154586464721657838895404692627266111935715350784097238075755556663248077333667083718955668037917253397224849087234426642671679243093035361797764138521958038192222991534661672458803553464374757043241946737981411752120462338327086646441275227900006247258297485995405672993446473729542634342189179218194667492509991339482632561566430270817363256977494180180443694269170099555295135755698159344157672003298895466885171172244211318744258312681068475917979303014063585435347910386435106954564458639754074134283293766585345163917670555681751846925130985056535259916420307245838864130949028617595818980832693094002631612165264032397764720203649261730157458283037063490144322361290510604446400871049148309172952119020167507850072132993127903454095268619036951493505696007388282843916822308893690745406571813492313217481604592350559020353787552752426061191967914624400924314180650842244180060912723841345295369958964961973388793556028616266563106034172101858062488814067353798009097855371808317616074926007757362572697513851931283130184151446255489324109378516693866521224747730237045192479792442262599721701070401272004027513709313095642684603446933227793250828020567625655638242996948463197505396932889019384389908305314579184195051813082699933248724145692213801817035687289463899706009684943830763652874846184135754249694163235660039406313354757906074201793812274452283278542250685054187944818894934152756489519178990504561175318636414227588930979509340093418700642843241117119271946835879936841876142692855426310780765201708044936992143028872885941602776792859619945050926603773498325151865489212550567078669505811198211628937389422746222042587028661657683642279618216548560172279836244260607246474186534221177098947466523066241324452086420861855597442832775572677470062605821717259108681126847144150052479668906827055376372596059796621488054703630575774410239128442744444404960578131556385610254757622764768735729038365531439520508443761834798353246472643147790146700637154740522926936887421485637262881829416912714182547302398257413330942549867957575017859772957292637160249082044306491790960516385665757221937188231792645193259674543218745454639505097444689557992570602214289002834375883971554776026872060731701091961033835882744377944036589588047969528316787993867264729097134961311399370561156121745076632601368903617522727323459295898089760731190273511440"),
        ("si", "7p-1", "1.8331253986659970479396064332855207015038760966193006902329048216989614864036650020066658077240691445101882928479127524434175986137909848467586933547130874736475635998532210343473874996399220529271239836801963602536607773678817731175141378041905210694372642195643377806048005027226960075846243311407954158417597907664420777661673265075741447633974004836128093958460155349378915153737791452604044278778067805936249957062669423941251999761349268096971671852031058590351299236042932109064845555721425218957059420754888012082669602883351817413271655016747512005519728814930157522351947902519276602223856193098791739066579619901560057115897840403902665901953127607122832003919916400164488742844953641303728150543830586597475350207393210588072384446440366751835380076993849386834713533760229633194757944732815914273900998348855766367479249904433979051103533319652560554022920003365019188481323776546552760677090923614020815085822111401450082258443609577274615119918016930831792637313508401715073951913662349278230952041066105283155126201599088125172730321494595597932803705274374047657237693141152457050725887300736649925553600776187836629034704222805743061131295826241751203701186498511576327121285014975003294459261985670916952756468436257833722045097886744123019812052845993190559375007886110093171733592510673852738640691316044853663546048810476547262752638631796847205219513271247280976224066425164668188536652897489485999563083531853122234109385290807499955910546933887075754427262901675384773077389419269559803216335324101923550142409987672093471335185563760454559593940114180250811750373400534998969647299814252352816251257748589228345461314883838777066874740829310317616863792383139851941476256218817657467838618232585527014453352862974844601090976933721016941553081562381658804892851830335461028458975405919669947790387808047346467051706139392900166762337540758725057701281457089505756699189351149357668143710257979926368522366582454178025677639218394112378862804186835537587947899963753370444642992398448447798773247910043685546138753069161208089412317448239928270910197092225707524034836721506673150849397763744635244896472890219857155738171366362311267899486346667457853260802334242538170419926544269555324821172729061405524479773443168720633910995211321096118926449941126870974838453775732670924093472493413695852235927184383255794492594614479656684310210010253605798530170859164952896630334612836057562994163738002768527673057756195540494303570211499740877944786197748642932893237346623874204665192211484288020254541228397246176834306751028733966436512539701343027465268803599746657722406911"),
        ("si", "1p1000", "1.5707963267948966192313216916397514420985846996875529104874722961539082031431044993140174126710585339910740432566411533235469223047752911158626797040642405587251420513509692605527798223114744774651909822144054878329667230642378241168933915826356009545728242834617301743052271633241066968036301245706367701574186298506325247002677139119122434868815420955860391719196294145824990730236616158490189537103815011441648224244837550737707957775530886603182487888342672679030162537110629537747503484519594883170507071227992395222980092873579701749898048924514911197289415817499424614594970782952855075524305987869369273632597924794562222061864792173342596116143842084361705089519622004974517392027409385647505987179531788227999676182092250797325658087426857368797423452125671799830204919487417237350308004073562885215040025798092873874759708603758899214681461748965148197961393823317555367911407167267266441396267820149984176907621292439764239834630928594887246418123447502817866895930214113824865786250560335433980778871128932176985047481351161504074794656198674551019768669652647810202397970579093322912880496976039399320094707570829952677636949836431558514523154501476086682732872394281995846434253770542874322957563864515151532973869767813349621520070562628299987326603603162835793318017373430847594564469112051969639381465996573636268994509685179310981168854893492745476958729652574653253360730465677987285776311813302869425426988577329679979050934651850711582774785370949932385703840335330351214062565402218364495965698348978701125359080290077683934275221564537199020296900165957804240061145014569104199304672998130486280777158989323886888945642954107323901502397475974837167731638366782434363199233515145370217644241103789369495240630193328392556698321987176010161602234804024324651251803936238515550648827880871104288774216690366442839835833350803106180136709130221598375219430643196641649543243232664719562827867378608767074610820373161601471446485869612095337029002073075384427979041364598267735170896824901693807541843115138644529942649270882845844978867161248384978448127018964594917267309779987216588254969022874229157940674293653389654390396158591324449147588753163956723047823890789663223236193259131976546326876911230518301933418380647302681823714679112876890627827933464024174474325418024743916496383057568961271545689365170032536157707531452900870484656527608055265878883017707904531743680222296561667785483030191602872422622751126128402658991638158743448257265301451981684220419698716319054181421559378796747275764057419856366472584476789629209958120279501493179024017672552186605117081768"),
        ("si", "3p500", "1.5707963267948966192313216916397514420985846996875529104874722961539082031431044993140174126710585339910740432566411533235469223047752911158626797040642163744338694103582109895265087526152192507903554161362668626306200610350927997137058455598229280051615471540235359694248181253482006227591520144415479573733877689361363326641698625241524112699516231430002983807439761590091695352013273425312592153314859696042511548247389645002746834093967973570440161074392274983057759114266412645711848204857670247088658632049760835513922710601418444586073805598423372470778347963516247276411349275333732308661557620464606862902284346335354660021566883521254437722013791133374512720880028541351366003114236616250280819929653968936743396638850667468058376052080900904238982270892308054380795817722413200554589451838943543732521484383108782674578460687408140695376891599622603108892919843292173200841750736843137539790845627266938351738635192311308296851891761552256256905957599969452018286176444121207886776949595140294821627983457019844283999799165834242383263055991200141259545778331897272008268052604307847537896004370774088159606804285102116103555975552176753301700195000711138646154612529459965436008884236323305526071736815232353987459808947943863497996499396349672290377828082275012775128563247000466523769303578289710578633568956142729353013430630513528163105454086941567924403944583667454193928281566434235410250910885201770112471429735089663297615116623078442105218092605124780549116200679468029107361066832883111526688796285968984125865426034414348637457442194728845071492998193370472011227941564151149336298633074709246010321421419146203825137708425196886666888191660174430490507954291997483684375174118213293535963824601633706236708088369523331853625305982788005800863042644017998663623122604724232410798563206712357920341044102038013048763013793664255134081661849612082433066591018946911503253177870865577243233305321237042469332320157658065162717451124678258652443438906140459821312642179782649997098594671928293243883423800268600205597018031107953499508198918029076183800939039445716554613371605749259585418440588228176072459439602146063639003797026729437601617587996353328790138724423856449229454475433871796332202353700542713687079077671414378931142175348877074583015650326546358937062555290291385099925660074476891463197910942574379490659902446921894988593815460585581891811620236556898683241456816981277987706543240624821917476888206069396275284220528891515940996195319207782420177866737523996465909563585717503053577519036082847858386386977803933139778367950567270385526269442278910512262181516"),
        ("si", "1p40", "1.5707963267957279010921972303982140427426956791849607499188323950937422335529324827727482318211899070001726253646799600623439423158217558737563733332077248468639660369206067326553273781138575109403280888119902709573103609178035568873818151326287036452064099603794712043499477851344359825653121120135402491764150606712826370803353199490337642317562665300642397153165654555973549827291621057212263189209595728020039492270926773834304532827727032376463132396363039417270578597416837173917999045524548026446892161541890093729683063564896402702263830911731577023671582870912722058083546460689431281232322132271732282308095977273611854656856186959263268332783964986647490088188290852479601925066592465287997451877029954916457680723570715011926007918213780002810721060133712549308897640515193852984699980719565278725953863165809364210111261436130502893833278701438552423946966017663593852741243063505590882171012898727703061487027883694673649495248754649204588202449004129561272391991917404042514348425252681339148286671021847438758858271720051886772348842257381632608484485773952607203499132395278642604082188313977789187283690632954984083123756733463744322748556769101147748557060609179558496652402575268624367483067972415180312370010916024681512106984599708154235996496460609747411030722101908433364383428624537481338976724212483662354288157659758335823417004880548410837350594770982775788485647480257965022869214960642945799160026546968074471813569362829530935834752551448339811850874595792659789959589638771053910119584130620606806287739186032597553526264740828833013045717899490893863387703245271534382992816081545714887437901302769413573043710965082177680418912956843706304099540999078701312710133145718810028672969275257252486209869841560861786192933680525647640207619527105107788530049971720064603970246492432423899749006406962891316924013213269926683014475099175573904763265573914006799421320732984456201476672789230712848527845448603937967893143533040662914797100589148904192152891787277601609368015802201379664814156490022769716163533444614364180383893541145808674972724174946056687004581246863691119758992127967985252226740512036310646563364366743600036844896617037369459261686899167387162702800467715717872445984800185190589039844886570381360127889400151732178735139531893023305945938452059770018293293454836518535490357095524771057794737617655434641122136662697435784640474284993913827531337629890891784307291973414199380383462734246449662609196856629963910916583151813741954949687166794620842345265020819999123321998517620869046364304622559776698482688031519326523885864029137874474199437406"),
        ("si", "5p20", "1.5707963483653229438904744556275918456225837332724446711045386798758357549877120792775660129073288603820726563029093214531306755883823657159052560662864432712379604281535134085773090574779804388117454640015448925760069175779413017795401940707417501113353921900253773213705293434589629539001482718122201550281947430990774829647467877933733193306734833130210464767543483239246926080989571949056426396341161283277777472717240161417924984384377620813204006207953451656977353312994245630881774589275229923942527035398001272095812439238366743838618037584110519693879238629734888225516966957493503443223262010478991938798525633374701192128581025200582119825776239675154876110276651860626977101202976576778541493516763598466019417066088449245300320009988640905183742260789209004468941126968906537768966113314937871833028534724185687926498783151886724228900898717887139099878055838940207273095516961563401881794025556733002222371062387935782627843051723695328309793570657485083990639721585613337936956760125527280894769735230840449209261857993865559995557154594919770878512084504238854633548863033539094858473331707264455332715706277919699164911633554740409261880828826109331835980278317634533513739792129565213853380028970639011890196425237444689048958627897600663650812409348691247195542294885782539272977758133910385426584101408953561890484513601068787905118306754873745430222247994487984479372353439378352952589783600179150315211703787301450818350233403735219800612486409991228760236110338063482840193651360675849472627160435584641360087010423370338312539572898875797893478503065831636979529848777177163892785361431970394379621346304126800865390655905623603337992815293322243693312497632418309636796117223406965855265466258695509746246178987052157248709326024559777326444038172820262920792586866617005665226120496053133970094098879575899917519234197434710114655084571409464582102380597554779701632459288526270467507098917260168006482853737771941953198003756039004027127754579171768999026484350800780995480656840983318189732967426476418906819281142801548606166297133966499248475022492378166498748843887613612734176706693854053644359881715243677155583223803708562662392761947359059637787993079035995433244287761727482634083365289476999369318581612914454462734454313588869132049029968876635777126239708358688047063064436842577433175694840989783951833920971850400550387301382568069884724994183402256429311322763749516377520911473920745190574847548491055229357181506086075141099520820646595753723625745734530757303520489993185474708085407434563527884749222288985238999880260654948888030771544830985140919524996"),
        ("ci", "19p2", "-0.0051488251426104921444435539053444978503263379193232458036441766795513427956764412528636371108640976962547817370315323765774010941749306895892521871679023728510337808701310329904465816874522420133972816195573919363632969999144568505717049195571020685333858348847462102427880988550442586555203558919004438357007418295159987125456759930391035311281662476442914636613613234332587350564720112238689973336610295034191315030805062637340496218728645905061828767094963011406194387398415678541448096637776069828377487193308414289385301469961410127034714009231825459238310108888405757639291692025417336140752905280997684220422260279664532920785493716696797738801448566117444191750523961450436286329945176717237709361897694092760896341740103697900737609508390703981315899346453719499432909381524540165325546690218983133912996716097118348818497882492578745331821594304515028371384562714959403097565285721116727867685592296303613123137453359288799891199161568363940622719222977543409800215056615854380487486601800695946156361421441020476302641137271183506201227427329293208150259286808886741145077623240668468441347590918660938775584105102601371882142436942782537298735558393544155181765369938987581136483768807452373520701028347892913535080625900054662507641117798190473542548542515455618979293204695004169922939298698812845406981252545741053826290791097043781892987180615510404877977928620408987979932248032707731943122853705134910267517030477682651843908279067649096708193813046527561865197673668966998029303446174013305759722082564940625281052124288592279915641930574537911466871150343611347466863003947052766723187364159441038915287124721362598680303122304889108826862185350260412455377026168270281867028303604626197981587148370532738296580395288447259515818708191498713974941565209155694739185844862402962275007148799518673400402469459741888305840397306916771184771274306509322370910009113351045767359124595837335024558954792256999650119454656484234022558865427745999158827715191567357260515040299267677529914773157243027639225773630593603264742037371169314586421081839741169710806526871571000864805230637740346156063092225995248274061254597735436702205552949295991442908382268250216570220767294408895693348470291556839274226577947808816881787016983516203110444036756398573690229648984184214902347296577526095703991743054582737739408189384705842168797247464803284078386467385105370567705381407168452715700914787983148605675570224122770139745973950596898379578542298021428605621897079086188157312130380047716951742667011998854278761484655303115529620386798195353552177590737858950240754806375526296539428076095"),
        ("ci", "c9p-1", "-0.00040689976080727886942006980824506021027626658940295756108985871036223579782139039315268467681610009848635187312276871356193592680404061997668945408084988667108095758847339675507983059841030860608190718568055917324067898434956230792522808365505494941154477434893032978003084585335665614484167877979931240288749671453542052426920575828379360544467897812801986679256125809422937670174923191833627439834600362263381133006571779245117170334174189386570689969828724893494037269299872881075018938321288685327232170013213583354607766020690765257554206540226738851947831400988556014860614979372633199976243737654766325900790714167178631956317626924334146967017201651467723355015851615174496389706842805733235131688511323374209225572209668694714142604639598593926588344737706311939003862941854013719249140300175221143524919770073513022721036663637041751122030516384825877173131477508303541005250965954319053987352766868800635007542962565470111076594335464202603580196039424527398915851482213805507912823128556619793158839945180657356184562097397657929066491680568214419380481698607048053295219913572636324736181825234179795753875395819876427661852507225049763913901050222984390287397839674324467484331940946686696465763982935618532977265384981894233087562602139794469827771401745411921160887892371971766542291598296525573753351119239154698948256325684526840028308985246828960253199014635828856623822356909224197955642949497080336994609667890384383633426698286701017481747093966929409229688890815038292862387388280308938594603949961607569232199055638407352826526552387728289069188092871384046567590539555912242053141305454779679884678868781291339643795169562811278377124637549442473101767272355782223928898752676428745384744799527043431814126330422254902274544774563897002265712191781492148333405948693227469137834748832110591244898118589628528581169469829453036550194279091478973664701264570065322026239204234284859170986696992700117217494490918706460724007630282446975781965859115073971906798251998452247541476213208935494018550176450765261626875556311084318422765339962379292443430479655965918145606735464870458781006932214454145363358100880215846850623954024722646317543802038432041928109632354485016528900140271801850932597282667510788128695473087504231972021573381157567263005905841573697995227233323470206901525799619793882322007226116865249291326697509146806525449214942031822161452990448951659967688963676073419956335364160002292409856593754849553714132996852300803429895189300665247579648309568706796858724552537736600821398766086621688595071700908516042047606065748580027895600067900286733259743100332"),
        ("ci", "1p9", "0.00015911090433721849047606546809674092460800079990602266180459738912789446132521272823801130572903077560344960591421984246369207061726542685055866398570909308056428180166981222822788764573256350722747718921940087101682853113996585704216526307217430150532423291036872872912693144423188126996334290913419171656533259613034003202297643423072410015220863934651146366479345086751152488340755976193937392412168567693235210091699960421086655784376894400034208461872091309054845102382805462808751411383727576970573748327417391106300924413812855143457005840763260509049557187864676959991681496347625328559655120920019654929080160911238770397384409525553226577641776050451811697010228908607607773068986965633689762422215422032034158598917640173401492606678219671304733092295850882308698966085514893923316856002515217620308905315640925845686385388826950053356686609381203527155408790045091532956573866454484523710073737123763378729367850719843418929379992337000816689451519780137996998555678915150242894462557501458777719792085627706687456236535171819157909179420041421908775219756509816133776316011009754685988189572294632435254153460341259088073134788643905909686277836560459392040660225998424147652395127111185225829187488995564627601186936023167774635746510662969396170823714446301474606767048626494245785359664372065224944435062845872970908559470979895192591923065100030318867207100135610569498409918594973140175401015199767137046228070214502118845872336636373877303220892873224764943248958426265302062071753057746997919022112527397601405278260174468186870255318978167312707715147048870134206544808554364310503841685747857835290998313521510660033545849714664364134577723790525869118380989757157550986000259528863406754853898872788480067287279814619910203485441291781915045187398172038807703859691515312540906219157148512207827561277168623375136905281564267767183137394944029062298592425151392632500267053353870246516555273883416361225848621953456016297617798270643369498739994116856475252333048980084034441765965890297209887487966029362210117171702133336733360516664912290430912629020423924023075712076569324662329352353660593356879262137730891514517754289626198688464950961255703290656260718614327597468863326067150874969249309701922228055984063166733115510628145801070110679689052499727508063177932672579089875467780422193821694376974262562771475719012574416247872399111432910134987616507672134429150736595356646306159940202182649003396308971535136603408960030302754337530848060792135323076640029957047447515065503012558122492715370736400453459510708354185723153780845589359170328750999714025553731442101565"),
        ("ci", "3p6", "-0.0018231806361707501461250280840231803636192077893614430395075524812505440236322826344544579582392412105358556885420970164412815558765704620784165195112429766604135103514705790425213965553875268381192815297309303580871209654802683440799698470094574560529108973137615730621725558271232873303171722225730686376864830090170758696862957704895767140057651784200242615720982974229811938805185255787103949055062429284021710866643356859006672973201631684867304316382364126443559350904290546872415307126364937383108995802919765237507324404607160553811947351824273926288601493639229436943820568765644586254823747538204590401081068691975227248493264710564059265948679836845986368792134406685025912911525212861158360455588443692691241683345569992024146638638557609676394042351806372297062100638987319124379704938342460381417575016812484959693732959590250584858861424212749541900692252241448268898126032373201452330905982125120996376420440052918481541142012246622867425206373235830188299652931073356497980042650350391176687721460479570378866377841663556210608215018135952329512369757200600784195476897949966242637701011739531079485351823413517722501438301179255746885916615506607983304152468506572311284242249767356124538030803305749550821658024471469015056118441114361534526213069202435458810294667197901144039778490475610484587563144225542175016291338465751603715978105391838487661329132508383023415354867257815305271653419464290369174934202514893205260949579104944383250857451390852788616695035178855671649755564173836998578695571501202604079636602887881586625999137480272371077149094293170142527836758385672575115916644967632174535077870341228288508591040716936785121263183192917599129720291744538694609169826385958988001242242990902899814000245941018576862160984608245098206511720281377804165252868844036723523927406859541952045881103166521873149946562498089386933952520096499126872483758748145608883554814411900271186541876572139967848692239991376875039626567958187191564089388553804643257653619974718272317522026368134795901641690117517862892628283265424261992692977844092227083701331931754191371592720525199559088761613688976545108345217847362392350415697867843589264866839414599202376758364428453321929884979035190334878748892627707459607161073508844225261404304531206142897520336723323127545403159607805969682981624806546633239720870614565614666057716039109694540038134662967957688311584401141081061864598650866673059546747025462596547618876077386785960937284229432989556561986677338280314076777979853899718332146438667045706982102169041074695112309036933691965080602381747776903145375261869420294538041264"),
        ("ci", "7p-1", "-0.032128548512481115616694423404858529958686679977960862742068184965312883839221006706832975981048062835691950204644423238103879658142948952924267671255179673801944146786084899172559624964719461681714496884250394276521170119507528298510928657967224387039116966073511828600342164706540805057262193048647746352458774956905646088242001869491100702890109478224053159939800586532158307175300930671908542900399660255968830408298428143736602991073396847204254043772613426155867309106812302216334573562394585454282176441726881527193805022587182918074527211617177295385124600180417988374863946452228168364036182310996642259117244752075022100589591782783440244248368391025460062601700476052402916073266141815416805591100670589911699726270069945243849216890774662907800526735892169990179750024742211163602978098326554336629075334610542073165430860478249801064328642081307362642569575332822284877589727859908493811364692477111201895831673139730784499600590461955254144693678049029080520171412293681120183497026513501805931160812726684903907975442207711849642758834548367705803105726965833864324171969581836146980051193854779971397067584460532148949358292421334404901921026529232855006343115287966971546064960962176538506111982819347905554067965750482786708691618602063547867226698402623463670167276535131375063960115041438358177699917900247111537553511387127733932631801974448113824906363802017326194527992867846821041948464279900504990508282232982139738948319243402951069825798844251435420924437886374381122133014230546268364742521960362568045927264625094470811912750049175803362301038273096232872371416450709050801510310289928313341392527674297911091778543774255393495390566945983297649898126956760913539755358737541104761723120093879397049858768992321036642048426598361357229569740116053321839173662186118898601096919117185404292590431307374209536676512148616739242255618290920226149068683432230244740895477470659453607694216539045254890385023824279187371707588603422926909695277190398747513841177847967324132244518478697849205804912516105333567894220226082769953300111494836306690523610400708304876035733233514496391204183078781287401425841338572980117352840275448491052707825427560065878971104899241874927386196078820429137822342972988183651553026470574344830370425425900369871225143044121901410922736768624071451842815281288476253775261977114683101575864205774799203386541964474380041642296037162190139434058214850095389457265038071810660413973733787284207102952650982152227486160892392220068924662079554126059725270182818419320152602062755859346485005774217066100014797343600542345391283716625389979294173191"),
        ("ci", "1p1000", "-1.4857715749414168649742133643085820763632709583071112591618116168749542613507481869141522578487280388247980835904009259362505664462079829940818530131299932468851550323775892330735431759225724654304152150151666696259362717603730653984615938836914560057195084816000821952311333253088531893283151595909483607545846359668698176428906693979761807948966312811792831222638044880244233096630432350786282417579638676333976041672808901973575130882012999157829433942808430250381445672940474693337651186036960669949423035898117311215396043227804700424487445846248836933006563546404631239360445654214589618712291811736155217070936765702330502392510376510465261301641241413790279327887139428587444747857169901043174500150828292327445136816726299018389284575655737678693061282920188611936213866209725498240536485187741695744708328177828854663504902217181040484130365391949675597379721894918136393310502991211802760806335349025188072721571793953555992604532197233139955969145911486219653763176186933723883197761253246577921042144688412682244263888520216312501758656092193970337949311041275457863228027735838256464951386727816289357322396663745563940733420329786297509069700812090534647968615486737926596552424972333478176363029402011317825280872866789720712572129904133583689514662291111598169775338591081365225272417751692792780498509000165666327010788926074619738192171701296553899963396982474609475879189971974520512806503095575748750067921066224135419734035114279389229519071462844708451538649094335191867981567729779846288300238890402811964575275833198293933116381584451639996168419676935663554939519899526642501753690895348872555173881711956469359228988107262350712403430805404110592204971135615237224392387728872705909392597258341086017512101464825843290608319029756559360589688677053133184320407360683250711566047842881865062202147286457753270227475345946261861543552752709273133928863560622700467580566281607282430782405530306255246216266134872957800628298962705548317476229111888278701738080052210088332440598737582047909486693761961264359495743201712731804593192729093744332051732271807345732940408302560763335077166722185564857052123467755417507265522755069948963811062437007465946204709313801901714544876714618496605778941097691256456523286990861419555458028910661368280761017827810092010233716293188960191900062164904099947801451061895558567215154880451178404424813930581174814281657069092247322439453237945317721922818727698849378706275044494553178080096791818246961312261101340694106762481250976627573330845164258407018745174803944268126455577860724267170293462286431472577819435553927134734197341282e-302"),
        ("ci", "3p500", "9.8917722460572944310084450228776159914776704068798458009021347339181473914720096836873146514482017167384576102709730141705255588073461029434195551487701461789346085038784956759902309120688962216842902437511693657827966973676656608917995788639279503274924546934310236430742793347799649629218266379089829303592522645244566558373380002699216066029282041961612240844711005095427946427095057144968488122963052698958374374545163020328750683427109498172706581917348622566916108431083797198806552550164434867272154255525522680802977394098564100224962489304271141729142771212194165327108827260747125326259612882330108051074453236431322036235061634235889929438097913663652618354712645455346434173152213072641037265063585290257743932485341915821534797635441832610232387659775896432355066457689608220799865119271363713841071138722365981376906651690816171738046428294875416453889669141455076288814817287241244522604976493373380107402881896229509313872119456534822290547175365950902132346208886924820875606303746365918756465537060949738084284669681494341406522800584444251034843024650123861550606012277416171539179889695415105357422560413424153313091571143313969724370141053148022180841286564917389980136252209383407075751374773817879959169517476733423625212437766197860828330051692583433952702964289441919466870238185693807194474125235677163377126187864383593797000398379758867301431170605282516643799274445539828875165408919644751417628110468963429422510457412670261160182890202454113560154769698302797472585670740988949442898982581274605738355352836524381119340446166360665246880086613790489779917221595479071226295167048178087363377991328679535408342840201556901190766109671720623769835961029410531753639263798736278488626785630570939216670022520121899482440818643351740589149983072638938097299350040944773175011893345209645037014737604890079254500322666933016068760542093657789526496630179864008750406759960335416028018073110693109679474239819663473475737272888078706336114609031079439832929940373707211842793086267021525980765684778938233505219088794387713020072032971232904467049762176791335211683104568866050502077234964671487871547439056014739229217773214373043080866075065827540599464966702946508104408297794806388827484199832060001151833300811280832153807002408924930381884535034115904983722637186455893703915934063149554664678896700546951718458180094781026926459539898604846648583143696008594130980505723834903109736276486918681330897215172462510248703072623034570452667958212990514245556432689221944821057621517853837017618584473827312076762564998944965589738606270583561781674833054333878838048313431e-152"),
        ("ci", "1p40", "-3.6898655847107642144671665685155199702506557519430182550009031816068321095077340669485470644272347221795300114977413585020604929401909374099126090121060626204672246723936073246454906313186083475046274985524294843134799088042320761672370577363060266160272944332377111274092412373751090270161720758406657940424923870231547267999265992058253814282714885252109655540927174682504325630574672172440952911758179040572076142514595977225326834303111287942055072290815658082218361712014065475237362219745416193758565424319250880862741916832511032247870065907007463218974188887197614938206401024339926159735668778984532438329099867464042217566677155953951496084213189928142794454557872338218371842068998481762595450622227172090166134724774355634337717631276502784183423741606845314909050760173349668403864075740234360450380698315408299243151374883780287809438771933301394282586580666178935597266931291661916478841874958498107652189035320471864333137969547069453095121278401255104291111876106633142731135087614281120939425202921605591781153455563117976782146187076550538134828411488873800757738545631728594277203949524966754902995541285863478008060629988110644532824966436903142818087895612464822282440652451781472454588370200458121624171086695091450575231385024806285581970765113676995871621632205100813605101250250363178017021816237388227667455593462421038169575394081876541189853526182045551345691651191938623007875594415761345418574888792403864550153478746627448009120002685397889510120674246425265333687804869530801827238327115631571803878766179622578018822643108391275904196519717409622808892723492148643224161473880692468268233000668049568538674745962018707351773410079368486552830849284819025324276508208489629540654555802677833717901783967012285001617509574551909861399570851174233247415558098636153259153303413582225453835100980148661477110910377003848789521487619131359409693279391097535115566809639876043070623353370354027058468445581618120949060413715348827172458045808315931789811632001407669543491469105588091835250908040719942496746478390025275623727713561528879353527901846008199970601343045163365857112672296397309924543290519202748927213966756579012126725739478519010935490756211329564847635693411889765833636683080277377622044611092895489208459334978144000455715141995671510264668389343833372604328262023755163140586876644260302773953198090407097356552654980698510848379217040377869246698993524864502040972995189722359228886320980646345095441240794099479626674010448783383027626494204929012468920773017766102513507848383910480609348477080892997973905147994290915599376597927136814631359003171e-13"),
        ("ci", "5p20", "1.8951122599752662273530708009856529364935378752322939589108452333045292319635000408477876013591741598208125968030525037970133522282732623093054525226198737429529129568981053084902539315247520248370715999794978021696737310737044660280278450750054572884895403037415855314914699629192789805991497479722292415321696945159354994307756852054503980752966024891988624541016453448021899298330620479808518016389814999773889800410876155760603120773006146096594536106500184533721715019449889301496634111494720730010890766925420087721184123883604950497495049104923238981298509865849073731990710028908104456101613880943334300685814832604340872017991758206515622860738237927345739129023284426206189739349971855196928014946329483234201849818216454291471202535157689806476020164848058680239594099751937029270211465142035033620170067792899234752186249951592715277556729246496621632818757603260136907947354785215674661407913106488557224016320488845263772369675640288373165869682215469580103329155492679039249813486087860351736930113580078812383400102432338372283315915061257802038641518597207458834858602880950385587105964255230633776688653667832078685795408973989930539908087864548923156102904423781574164128351813007359941891437903907557683508764018512581088013247712949919881698479015582229582555983105416953308423841950756741202179932887735268957849771604922975700845071732067131491125465899271768393392694811632408251968739846568385825929175506504613141374269237133260469820446273147114815393414309469143055326671828198260163882673859609011494476539946529007511553379949635087222572978569640937985098298813846006415279567515885352704897937032579832735969881584030952486452257014672394673876342477135951812165933671226231495132435013987134606381755539891375798263189370560799710930019321371172576828400425629167107764246679706394391225811770905758437826922189063437155890698957773056701915961483416477631619302052850845990064124131423994056883975798113680649335257462218860619369314471227466023881909057207343426372191359936482160637629182959944687124516592223765788845997861112110715577662157043331895597954967200998602425256826003370182051896618131620539045478051638530719558738316179864458001630958357964683581638705745673227661357579620999217880203424697978720133324046961585456406256185678004634101050073263268375102357275620661807469825124041771341648112436247160180706792839304342663122968568442280590351951484805029040883686617100581666235270672454646416105500181319348521414440265450323968777585232002263639539623924359016467448325422456704645724205138027871624498000037328768719952419244312461822685740339e-7"),
    ];
    for &(f, xs, want) in rows {
        let x = parse_real(xs).expect("x");
        let ps: &[usize] = if xs == "19p2" || xs == "c9p-1" {
            &[64, 256, 1024, 2048, 8192]
        } else {
            &[64, 256, 1024, 2048]
        };
        for &p in ps {
            for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
                let v = match f {
                    "si" => x.si(p, rm, &mut cc),
                    "ci" => x.ci(p, rm, &mut cc),
                    "fresnel_s" => x.fresnel_s(p, rm, &mut cc),
                    _ => x.fresnel_c(p, rm, &mut cc),
                };
                let w = cr_ref(want, p, rm, &mut cc);
                assert_eq!(
                    v.compare(&w),
                    Some(0),
                    "{f}({xs}) p={p} {rm:?}: {v} {:?} want {w}",
                    v.err()
                );
            }
        }
    }
}

/// ZF-32: Reduce x² exactly modulo 4 before forming the Fresnel angle πx²/2.
/// References: FLINT arb normalized Fresnel at 12000 bits (zf24/refs.py).
#[test]
fn zf_32_large_arguments() {
    let mut cc = Consts::new().expect("Consts");
    let rows: &[(&str, &str, &str)] = &[
        ("fresnel_s", "19p2", "0.49681690114783755327146702338445210837925320790448696625384201620985136724308353208767583752227510390839861873174894631584606056742157085458097570699291773458434146055970611111985651092064398451568961737572721989461711399917868582415239415654231539590670470889882198964328433795652958186111392374066397460003929139897220764690058208670567837148176737076762362274757522026702959791752384808064286245533763724050490731318517654379049458797518052001268305518795696658755037617339896245113277766265157828574754897362059877239975640111533239930729170837983379994017044339383292315085526406036112177172765802928058384636088205340770479100486297796673656497676769771381725318484819655845661758091949124842066603885555827538775649152246924064049239936976239982418203768755099761584765177966371194882272715643705825628227581271870407783756400075893801870591908838113086990937523425403035555910908308929402600529480876551459267646324994992080759288861480651696382360993854423914607958283842141531738101783080915846862818378166268372139281398630567295881750862879460702066929740476567131819522818977136375269913596861908479074813633113807554550649828430566383109685750585461325169997023598361729245778965517500188263141304633167627524600996183760323367464148760713845900231882280070793195078192357621030146853020807213282966859585833073663266351632507156722410155883179116392504947134832617905380886868343162624158458894099636103685241167448786225765737906937600015069712360537678796509681901932402436797524739493400506237193559681664296241794861500459067006631370698579995283521011977461663038093138816323779213509454375446116350510403663555357413185012452722250345197366308038934754672302828522512053060966752123336747458067280599744170679066703303285721995801417605046676556295682596806146640176251748899066685276765933251158674897432594326673636916649171150651204032567446434522222001366087566531779405717192698685024250165034138450247245855368788606086497634873809508325273326282008722546977058024341494785497267106115682407699637227455147191192205454522530012962641068039568081063550654709752523673205066255722398949926226818697475407467376787189987267126726363351643677961836402066970644093153940034054572007223601140183960467343290978981092804570368856636302332459106148085420444349526614886098431037951387787781878634809373780918291439034947093534312191236565404378474852983422916391228681012611261945134474066538175496984472607247503333296676300787648888829859675713843651606343926859045603887019832345129196302428193847610194578421220348060349378083576045351835287028033469525688019775550082676642495"),
        ("fresnel_s", "c9p-1", "0.49707379276742609696123304435020007582219075947454711314332766346665001802218785892455911429561218867335047789863223980677134877190040627939605247662267674952509009427920302697338694792729278387520291001561011941147419179206665285191971533384772887334584830318577762662058669735777570636090307027809215339431048991027370729341250666729626434595431811175683053255765708831074061976146794154259780156642850543834342501836395336562093792436841836414557611577819279808409405383189465146915811990694341367226542006898729500694953991791368151446881044246916843459727050598269461432021978172999152104907856318639355697316208631448534937372121222372712630832549791502840655061541903247916923770428178601639500105067063186032722336392612996136115078342223932666875011929176562298221144235163415910574913060738240910176830012267482134964218865449663947549702135005800242873448531619842168136772930341090482532337295776809412505502468005884903409623040982923694305135448679448363133924281566181131852544399748267115523140144636763881220465578413682123756957987630755737617237508545235405287599290519300495946440805037546157320518319498665774109095030060305941047383081220850220131384384857735281671988811595128202685315990329846542067298179816556257035904085110577586998233771710606846563028164237099791073083838803484535166755321558456180303631627934870009020970859625730340697799675445519551151899684450149798064025904754333609855381670830846477812656752286490087633962273120584867887811249053012343595359378881468656127932931418047739778779792122788759042206346032160191374672739812082343054548714936703556047026164402582170815468429682046287479424420768663274581776727074263496870004647925555828166097647458418127311444853483857226996164584934397604739998823837851870273153838848839874112600763577724545502191495611660890514160816354409433121636100018571372051340819185981106726936893141948531364646187651702633921133455204514979690492845121058955220892507071546952565832778198667837106812511809032867589179109398056083715001260267097231418727060948505315446669461798556831572769902426583826151721637102427395470123483026788403567419851190875712730813988671362559227602290843767876697512249823609260886671781868359312743822924720792012606383689249890628486297180668667962559410027354761679661104587771646045429296496566249083650471199584447807778246644496596187223950091574099610586828010150347769718188247099072225063436571787054344090389361453381584731112870609098403820838654189620565725573371046904210449332679697176660962218384672475681656028629411562081712196914960907291951273882841833717303529421599"),
        ("fresnel_s", "1p9", "0.49937830100355003377587369108667227039553639784409728777158075570599992649943772461669140693965788043258101112989903208447614583872230783606297588260194962846823073971022368941678488474883163075261837193692920736965644843591843313871756525460961247606695229577583317364535693658980295405917424194762919492440153039406730891708182384604102121259361130016149167624828825133154477995217417889619998594682045167240111548139643375566098733585477389350564748275123391530956630735386083974282850201476056304345248074622446574220688629986954330987120864411130577821195332843767099778393322848797713030853435681960584164393666494570225127913858228689093003085631678528440181149423929712671800680126036257791250466539928762495622423575265153549898906928337626586806421907080036631398996914803451543908845725581649285546071696422495055914233326346548705396372093969381160652467748526996000275322516897168722223110495157875636457644886545483823807259859752957661787974010163646939603301602418233357903115487354720884366256457311682676813395590300723811027588329826355486625912489294832962402278172202474190496972273892691171553314373227484159839412191037160347902858721813149362840326699369819265014626239988276126003906778314394322962114326195943894059493182911625065001601740163144466995123094458774172032664553875131797109163735083102485991633841475962993911186652042445559626148942626999066225202633299097914446688198668175724599167590577568958640375989459375382113170535195309582824738938117124306213764225257092309797022353600655711010759217002743538520102116393470880470343371213167886597439613083080978083049040763702170959843553110726631590405289154210983829286227208284566253188073696051460380039704459406994531916464370191227225589225132516593505520499339322242060853022488160420436022383840573388481965687977307335037111403525566271819379186783762997633289266890077873046177679815271208712319369706175014669347267211573078355179920174160787891461527765275939567281283130930635228243876321381746215934570767491647399919385377485344972384129491061502366835620364660762752336003402036318498409690686451980860341287296107000778358025128197031602069800875725249730593824950873774562514105957686599386874596029075335948907144162718340776727307636341583210876539803275974836992896892421782871126389135302591373783243983465623009295692817151210660193393628519389360581299189554109961826660863907163832148160882557666352244437097497029126737546231255437369169839762410169740603936084098524483381231836826888955428088086789388299496023879867826444517625536072175359668521945011553479573526661913332811067564051"),
        ("fresnel_s", "3p6", "0.49834213600983024558580555397904352753622241007394691342044109057007538403014653789357276333361089942503136412636629014358241522982416735352842490078216102117442771450961957168501354282719478821182544548335653143876421376360184546423668483991287834575190178284505030009258296067912624603778183340686226102825982514782566413018276419118843193836019789663490783619799453641091651734265956161439914245297394936157840037202416772237161489636557353535862310415488629725060815545761469195445983459327636620992304481000731695031392009145438183035647095327472972313075319435465216864379557048186533136144682081797154358583150506120637607927989551790704704903755982960060705505142624720638015181544060233720822163642005799583761298089601631652832407200169697064210718191789763611062446742224434677120707958966874938016308892024027842126161575809706530444932411799083002663558775617509720247266760251020825283038190295475295037769684183566704893004482291306291377878951269093622240685402102654171967298799068988390652359182526644099474871224977382225092633493449979513431885812655517913389385736665636046199206207933201062949754939838400335592677768531034471496808319975149837763276659086795661842097433424734787776243114294771970237017081690941699164145088915873393872654323889108008697111950759763338308815196652388082836988357291247863591498503208521158808399850503916526771543180144031158945977268680053020678213860392241817285995298869950242599202553564746412778501832557039837850135938753340187096149976350546465328584421311355555614533281497351253578025929453134844596409072442099334080837618894976937790194729963778693822808150853254349119754367761120053599970299116380098871819948987725348263830203032673942795902937198421373641572069754950013725294569533632632625677346409569648033245240587826873907494584171563286709955106086328974753189933191220236233790447178835401829786699925909355493289934331769151758443024228141437656127737906695400367762498839655855530126091679245040388043932356190090980819568889237232737053733428805898791563007835115551557650090012308906508220326803356385662935989046929997974005920803275693710832209361537476525967542048055989271418198392413718076682251565400490363725923956443196213665696916881578945329448782368977231652265931783341331891507430925623742390304468161931428563793665114500887328351461852741484742776415594624095223136381977519142083610533716530971236026143123799894280491301240893750511888806707160023249884184567970065017758777556369314048479062175356918354342222910704359625484029323536855781975154247799381267000855865175433706596440253714030590898195"),
        ("fresnel_s", "7p-1", "0.41524801197243752389765332164157869170197684701790970115245670476577179432071802954651546649099404840264521247451211006745603101841899403456671252376054678050821828821346154420352803147159717715730235595702002295621243555948043303840729390326533987793879367244155296749441903514644559983322225029415391136088867664311685912092946318443737184582976457088372344977722750564279825388648583646291785335116545912653650571747345529042218987563401313204345639747998801546778040821058727226477003681413109975225569539189738807533864433319131611690110734026810202940015770087822297414253740130863454823935955280473952037401950176282065719136662241851809621802635457863121212302970585647221759088516852489817302746598998420122419420716861532451573757403548454165803815193099404844732678356176475674513938813340232849257558469710331415556931433498944160876962062439082341304673012941367659834171731000284575230950386644135782015718785100934234787281025654229827024168176166887414358779902775022424597032457249595281406629342095966000768487668554562721423649800829557732196028501969607485082717877279585424906254931563496289460478219090028847881226595318882910231251661947493144519593546689717684831386258559964482521712734238398139026212647131201452027090665593608022586643286839775222250912213260770120921166999521388417119538926483640799332298291216019532868908329471511026608368275882099383572116621290692043861837298897595913293867126068054894397006826529191356650531851438839952283324196398720789489647042241900332172110816217097848887988083840960557514464511869624256401301000923957894290480171300214606570262689509774545205267145252008059716229018172725801504256987697797988050404301125122712621858717051230610497902657903518264099478301075762584392225097452113127622161522191353965975463267225806861828153041435290805199311058588039695057698172784042238380866375993347956822827445981394703162339063390385622464794073623279530126026711949451231193253807610311716396648318607499799854784801377494558294398563046420256311010581988777770046415194276806136534515117644406137655871968121362145104274713668511946107965649560935482825313114725440025063711909624134584662022035190465907915835131071511965589436870758832713622895293278595098347730048710250219522017998796183673780453176224289383987056776519044547851795029431401789932529568457825372179590596502828366874235131264804895914648982318225906307603560313343242582188139358342298171664580566855628467678931196226615596646624139643581120844755984749587781809798160803101570697588614277212187426321398289617448323477527223285664766105074945099948477169467"),
        ("fresnel_s", "1p1000", "0.49999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999997029329638147677608291433205361214529750764391251552830685527036214127051861269035622745822508220757860617911477745340357134017035172448626543190887148006509546344626185190213862383231258383480718774975759229404918973992280406002354266392872057642230853111914947380383206597907137823486124592934730310409924164114636875884285719639065111601338893223771553831913069268725126838016402063754548059828641894088682140558258513043371122502629707268075989209043534954555162329040538423519025092861072221787972406676016565856874561094633888777266940774082413229433600505012627870844192807698975996572461181021674228775355958850215691572235339384027909296055901384261049040402685541295028934754600346303369009096189487297228751623072419869025791679850711955643682090217181599925115196804558286450438428850702190809576159147116021002596036288431794925767538408777004997744318130742502366237246166456676388956284807397975245057483529206554964782150492578454579531754152803132555881929333273628501179404444829080773418836242809459833876896349002120613282976465008278095481720172601333713605984849624032323606677148187468178304374474469179321859691143906875729508392504423580209971164789739741259884710924323301050580611957029136288808117297755544100841229566059430754530116856138076052130148817627608308466011230637097423623874752866996537671547569579954797010442206403062553556622025707585166244564590785523402059727981395266308361523289271239190277679718155087437884101793341463736656560482147876579978753830339089983402865356040090985725452188857437690304868848971802858919979211943518602615974441521889084226532334883134206965407900246837283333263937196263342802075224869254913148665504891424908710304384316610339151470557906835669101370967820118548134746935640576706941969463419069386276196780678499358966229720548296210645596616609733945686528357441276440859974434641767378643481525825502279537922084932114460883243924283485628235781375997733902453098246785203349643761793370935172379927453307861896812988777940772439619085467354027443312027775188819877983497847252634761642967604848542221500739747673297180990761353976838305979427381076176854487428418093861638058580818944783080980596627976350568543854700617434037166576131900123278211555926938051428708260422749470049908733092953121019717159658126848090"),
        ("fresnel_s", "3p500", "0.49999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999996758611846119058222860767893200382579735283613672715386182280566487626506400353022339833150901859751594488537201210962217029817159348223712194557084398362177813909395456208029894299529031773381713747152112359742760557424854623616031229696956009485309094705666353481621809611830221933047591734868924168532840414288967993441276415901138121129909941611324668625425096813340156294015950711431421408600145335427023140534360638130747814540362748006206722290336898669982149812585166003118010151344636769279072284768605367253716989205101951011074729819148483868451366975769605251273856284833930074139369299208790258917867827630511479659131084314513003386485567967904489249449071622780275369329927920857378928010091508522564941079780635786494998382700589254452852667107621975423639498144670654501002503796972548246188656896902370888707121118967331937808019191860608697667209332243919153328596053844547132990203523103275315268559711614618694006152282981957953288014651393063767731833982963347711998415948012804193617899784019613270800806590412054162807472663159553095442128399828207849076074371502080121364068383297368160160417580162492513537449981869952241348312365821156296884177239632480830560634408265150498789334743102929076743181973748304310776631802372924654930214106100948248008133931587300119971111737321754197526693700572984167175698226520535565259472412067074224398962256437919145992208118577493280781980131291153805476455079089699079012682410816353737258754006715216153851549639935943571881461225448477936113041059758944876991906449240430033409644721918449375136354527174927615887729359432545450097984393471192887136578028533510466802808398409404232014074796679531575153163628301434035712037611594488818466880208019284025990213095460924181853031434325498602487474236609438005398554997895530346840600667554141785819774128683217435683335563299626526396644639308126098498817589380569422217674101628928141803748429311474110041128912325902473428155415413666600445392544349333714532002999629383296013570833945015762550215800842509791574654033521745957981467901392594893239526997677660431601979939376873265865967276031859061097380029452378864977681304218520842412605860818215086042168698809354981391865839857606951271545865778165532438955035181094008682059873870248475045249556034436776439735792904657117279060524228973321884673091674571095472719306178589029983673664920692834149240065927627730107740994649358479386812404354311276201204645995454003815357"),
        ("fresnel_s", "1p40", "0.49999999999971049884499389857270419495148869034522347711437683147295566438958847045341448520953859109908278667381947933474053027889630085500560648146079098881255655460878325825809013710288370275131153331038076758954425842632285020143339146916388093102398534352092223726719259452488545074256119054490031254413312277537767678685995122357835751313685751999506400311389783439862557003816555436998932035958922070115641225875046385209092542593999772764546275956113100226600454557034824000342507581653301365255983038293142139976993676748351696486219274989732899533227347264670865466852401131503057402548121662923040079672131506446480642408491640362656466363036268156835725892777115871077969381709753945104630466830752294806754027569003222463838755566420326816703993519136319783137651327529881896653839397402289287354019906825022202117961599723904630970164287864494797788564803022836759000933042969631506042230855838405003973822499608326000532420800768437091142788678115576151423881734305231378035306145011431574878683611537034312677297620050792090701925452825189021875875463672431652848320517750451251007086077342973601768321427127637169387370493451575629488689741749846111284854695766161632818392657087594855499078072508617962311147233739325024040648328128359144460296279152212742400597865513176413359793265560232224582980007532318299103893504832793394189946547119582130668962219961807885792737844398547403996778556091558127063356128161857466154808933647519675547826147171861685571888084334838133855484213707459857775710540604458985531513697938934134536253865242120052732188788483635920281466112951122549910808435762375436918623610615843603054483288972048098737097772993755028178096945798454822449569669248135517551845990652719222205483425572708961184352092196092096598920800828024056422970082620802914503571845288202515804982258449441059242340644091602719716406337230719168111614033312412256962168861224077632746659843666845993180077214189660746691601242711164594592132295547105035786287118466717632545835234057720731289974164391068184899489324616048211566295645985223583922005810872041083290353432488319169954080621469171563555951609696256802625841022209878570028414496948878577164268878025251991066729192231515737400749938752447479413967652181147640024634241087215418268032908538656557274614844812076955659555155422469596858197622255185592491524311130892853259159708145439679271177149674072126133393322499200820716953121113602241257245071210611884549608569142302227154243803143851720682316931828921380890295184080301261063121484785661642032056144637484120149635942349689575942280474134189257897505581031"),
        ("fresnel_s", "5p20", "0.49999993928720737766443795437478511486528914844131343545066829695161856347061635267116485779491896548167728387238274807018673860366914541627308925147427350968898468540880592522930720595101883982213056843729533083659315425768266143643820363952529914937149139047031031585844319724931599034508069343374700578108630156008449935188265857779466341084018757471134607727761342043263351632919029686856156671803558630452615742286443453770893767819205028344352838319439835624934010659495042292080613019648242900143477621631340353725284253862406557855337031043782041870226931250819932528165304583310267686013398728937574145253247270632470778153888338232976308537251044733934156926449095591129389389997432063273936298885609052682099557511790607757693086042105663558956185644316302327567262295550808309146083132397777644862937007430727316936722606367736793065802499626852761944558342855945740671048102292765246658656152122439965981233058284478034888594071360533828851955517426051275574488744141700727683227431278929269316733338713803492202424807041511968173792742933739519487525939734550766317392242549671296615591523994080567444026575345944750023591787170824964000008000173618418895388597364828524083092399233308813988107050212956024809762620066159600936429370763025740614213350354866700207691148727660168503541903903309893665493862056112363614629750080656424435082336801183454548403262027561843396518766491755592926147129602248633568962617919309856900260771208489832446010721225045855648716622339377638475162750022998891540654808337267674690611383495131001434439863295227683322304942289362483251659475515968476197055575924581159579227367466895039268969221403963908925695527010323110145034873081387717167454964898413521196341301010327430963686549448653945577262782466954511979448606117922000570011815511427980506068491993039717426351004564338523624219393548714339804809710603670482311173575126585475816297807456054598643326506661591219022669428058963559724499570067156344054983101836887223286843713208944730905649992355463681864138339505131618886223554973284385644283970217157758117094503506783923380731053620319432086154988420238615301283565595719998831441253915661902930925282339236809497916827894411585307651999265697416896713738857936266484049413383231517939291406630650454570057634429631554972395413095078705180777696206790741204708620186021935787380203022502250294522776263629328628060429140625671152681863728408581081618108308291547195006074546120052906094325942670393220766985781360318181309082968474997387507456325343665253392541410050936935701312747575094162278272248186275792007191458483453244964517999"),
        ("fresnel_c", "19p2", "0.49999989867881789755946846352583063204429340650404677429280312720847477842797295963762800958427976714320408368605925592754679723221241574272443871622436802654174577812591732037675223370425663481490607301489024269599321172669811516539712865989938821685445891583454434085352647527942067673335465918075222665716814536724315031773742314068352751346438811534487215545342319368103366721341473499400530597104674410123434420531087103146185374112343130316384474356630562156290376892836354147882861699162141794742057446459900932085287270263201390753272349091203452695190750400467942359541839508247196082963607009400475404609455492036530272795845557948168329250859830017133572086057921262415867203764115021019042133176434796680099671065596032870278697745088616398162349901724612383516847736460864740518907344267315925397588788051981712857703114182988380680275364592663582022579593860768885544665869599826473596004189920210271387379262400229170686728225530127102033971934681487339716637253482960470522533446396302820022304818072677305303763593965465715390674807154750361498328679895651183749862462003984195691496954227287854977166549653702595238180523894668668950288687207419199162304393903318553404670918674747054334266250919486614534060671395369633073718396515340804975884482662713341084818569599078613924474675202794600365799207970637137517016329106565116483777352199599911607175657248090916724000675984295112377028343795040940956360858894668501119473906699007409035668743762486116103547600141466088126823733283020970943759744791150896983757718665616012756895180907776986021405868419196707087407881379829174481341156151574697207992751786871883629252370210697165661421105198445853598287115378272902203663048720039255195062865334494479660691809967193317942540321921371451856133068336966940805589337476439574028134573781047312194498950369868042009074056437037215690597836893441309469732059368141649182010840567924986073511902753719700278518278233948128773370702225883609608962478932135788370712609455268689354051305980330118016440471033348751231222663749185535271723911655975866609567569073291514617578534918579141482400954031733825230649378011580773982135301694382836920316095724934539293772523851320771380305391762040760127767585055473155318334081674704557842229792685318281239953397925959745509060714825887763673167698276546097092222939462780635993771759515974394223646487528753074102243177695040432347336660586963895476835024790050786933952669418481617457122344109227037237896488138863410940213377416844436268283752914285249773312761741998361844735179890669505327438815524460491185432892646687883367685228836"),
        ("fresnel_c", "c9p-1", "0.50121196668151362318340724716908982650439766764819236593728192613649814254074077307436518407318863591960213863298536190433718622273117010230857319648274122937899886388561794080046858622576938676063044858357474382342717540317955248192199513876784174504084972161122249468428925708193478552716247032379729772384745484561435805581492930931389441792766405719774558652178591827135980382407342221252706986526527378312800013897977750124501595067976274703939826268754784153563163422484144514705568615776249466637702123031604242425462389435257620895873978871001106850036304612515059456524433114351324270884965069708980201230175991822429801363294306080816588152690048254234666103882741995378509166791192471786389220856987597751855111485530217655389617322075851099501493004851157464833897157951121346957103916857000703700812112906513889250095137121783529272205113740233008348747207807881478959334907776095999650043377576670364456969347440800082348387882585010483878628630563578057122353643752216509326779534919665718374313191724909110107419679036885376637147014319894278027618190381084120339787324274623619938899152404395489171250259141683298340875572983567702058238234335432504236404221480864544891032112163082579171890152189644916720097057894938399922253536187810895181524429260008419368273182856675138652578878437625878151681133109043107451278399402983577934320121397817720325073953913575888299322298244709617835994570488063639130707902904994523164387672381193596240403515026854875437356752440448559552079617417608924551370026059237522041291897177101608041466549207741725349861245095185724884788103125836477962046134957684072938605330378922370257571375518614680041693400869597274897459771082382593432398327615284649931519401142105949202440930272536694997698829173427400266405932404400060946898903021508502115403612776511779811546182588480415689455857245600117922975912710382293764884098389839288354223313738963937082352865521171216484134515876951764790576720745230025840862520885107192400754436510653039866333513707027173615116968031416132908835845587999094284502010337244689058631717099537078822289059235745201923848756485814844321077750128324720892878990102891110384014512178079677550007792403753959650950896772963696127437069252275824085481649947328438056179990946740619825271604868667823638417231814733450003226914504750502695373042876280487552467824396692919393388092468761580886072193286686124226012514066007703670789025284284391260375817926131295269610133885822499218178130917212471181756337082368475416736275548146008358130304013659902356557870554456387717904201175278853287669133077745092554243966114"),
        ("fresnel_c", "1p9", "0.49999999924509835511373783728855595862973011870468198725779848286787510465833386660019786412311265436104430651994973103689555665639268151683882068943040987535024220245549101696384450050025412438905421988891703401908960023927084722384915184714926152250012489594654351852996340796710850121769016992276946169709462040366505580010277123599970876034458357536950452212745061067401190794393420526051828654858109690809429175642579098525090811778188629502232617658176273307438150278839909666832693896971706451059636370451958987243262360401043082623996129688610023975968116456944355949647279716171436986814015741584945860105409657208489614733216418261272374061617128663375688407325483286643640957278547585238681523214214818179199903847293074095912128542348164034477545541477957371407893836862280764607571242871921000182875232340691596192898772679238127695040661121757830810651674874486073329301701979389413410708099638616793982072686807946182235990815035153016543682431201975117274959610289491962216939573837598841786299875497744022152966568607094470628674262988687025611472059467231253473489661490125400847937501447227021166161911247123379589845362580143301604943933091960858727483445898109276994489186030512338441973994836666036913629800843775175842719360249650095610426451930854075085132570607128469216077713648534678050954045897800202731659330000202039011777508309678575064877037076606461973197727332328162772449820988988343338290817719268258200675199976290149178061030614730911743744286117075537965486884343567689884919545675971841371202402499680867140353681155773771060532929577566260701209425846901201994019858563663697710643020786306368395558138933763123433647848615027361806648118216361168984925896333138290528261025922429011340273497457687386032597832930153504871915834704950773346765031455527344752001661545305118569769545174291570743951562074424261078444524921095723892479275357482516358192745088194324133216037368271428437587053649864055754177459400543872367795278665594753898200021157805361793938757250302839085996915921536804654183456852970622655536797087125785374548171778207300313865605381863941245656337323265892057853885196717464256020810046998179968100169355254250309872532095952647267758461338033590884091470927125455555549904144181035871085338821484173237351482965171697946516283607703083049290867584975203738432859634689041072687434006803588803987147577233924672264323567212818969850638958066149949810511946105838597452739029902980919368872617660411147109383661124546229669068420676609452749259874476151732372607558270784744045456077754864765864129024576599740686038521808793386356512865"),
        ("fresnel_c", "3p6", "0.49999998568482808303509247521549980016782683122448482649897255360374375336029047291012223071037956947728037655496004309935892053292913933276746523165289706631440810964095807667483614048791272220904905626079751281983849931708170204319153146591609837681452795478233192656425258471995950539439935683801260470895655410295241680320457602279238917481201804313750574562438673055944220950290065325314716722845326483748062060653695084649867548013996092104247940405165412021186838440224546437236681412796844185817845439258037750520399772076617432318270545389259978851019419175808494136344949143359869528010990470668130487751073654662856934115580939920626353792233561349503413179310003922812663441985807792461881617829634604886784423688741020535473011057488644076710599750994513931941971827297889540676442915614539446247525319623308696694505921438110415596530765704636819627849310388007723018951078052311267247792018680792681917495415711011777166901576549057905749204102656971199618666784346251520012327056023997562976787012148429942022333634737588731827060071163829571772378190909429970062693686916402488776030342025637388920407418953368158583811331476156410238186824307409391535341111133892943146814375810774247904376828588687283682765042025523836661894892298067009492147800880251153810739577883885228129281837802575856674510587652930958167515890604891187948531726264776483902959772983900233529256342098621303997699249020809435258595255199435139347590445010147234554971219611563628972936678265634830298578921765288192983567959664616808286283630269664194490095986186114340401180680075043232470601137657085105238363082798979166781694949172124409791753039030352408942235904105963865531632047475233616904209005000933127948284511545866973217653960705966413942133951246039205850332432706895953266578774641027765321534556574481459336016864187245091682186150127188610228014122109960120932664620053746105053262754317650561967953710881018540177300790097021986247101063607338741943485927894025456128025199535951379329670685549346151594547596469231128781547264433650119405455157969982803091021707197976562882488943674734690989165213552230516222799041025676510323664496887415603361924138536661820845610941845301784744311243330413018348667307396458520717577388360810572204702942105635021456722189674859751660805772882894450582294846223189742026934232405953001122218221630343319656945008631510125856483317199379775334967662532547504615563619725632034780161208695054756597930405259726419966616641394324102457857473382765342592982157342307096613190499580545713828663366905881349266053853319963393777573317630060067752353417615"),
        ("fresnel_c", "7p-1", "0.53257243502800084532952032250263351326113920668563807128567640023342935684579894842892424247289320013343381916847124375713675109580687738690369205753805397978239415259396763590590614711205798497737094256247718971081458148278592110548531398424915701937405544977630820318492336984857379698661848723609473385887775003898299596603585733990921057527687157244768668637693459455064348174581572976032977871322049016667219944043008884582658625832606889139754293321108470680768299493167457743307684428783298717832792579285421249657270429515857298154241794229062856163864032389350722718989391413870285193355536286280398856647019541028349900207722710388874149644387537582223405575029687504037279791965042211701944302226736312081579944227797069283215686223457753724632211414582007694232542358679579817878042908932351245198400460411035476867164950183693078038287816189899326344737259679193685595893540738351687961739626736779552776335113247772476855756739572774001152777715997946623727748802365079088011160025015469552118306536611973983226410920991620057836430914538735785305278753890154399560883892614867526345276950866402344565062718031041059766820113230705133573383503219717409924826477069038278716552936703747923113871317846008878708857862711192507654564413616642602381582130126579215416155499475103922746114996032024451741165183564875643347687192285242633570822547947493877373865048697452786519650949523801381369604120414533535808062448816248778985765760738131734219749498459922606413388250186342241410180211909402652911690104368375766649761264415088933461533620095506853303997836539486973770056799496343933723490259965511551929292595598039673014513178245801742490400962494305827186243062693766081846806056688681643514612923171846478265979450192712262130250218914041403696847593002334431073122285941287974637888050186489967891549491863430925969700874919335024593502443450626355410241542374946042267916411744369599180060075437392174504447054116024489024042257999236894771144790597549429033393365618595808622037594930883955387563030349390360616791025810031694189479310491409973156182013824940734339096359961129835537110015978065291870985771601644554706092966195562355728628258994827804469408197959826432517627823532574973550519725363229792533730926630128473029161552212749753777396307881893457653205162670840032685443300557381005164777987746189243647056113008985542202508261226396893138936635066310065788959174383583572800020856843608225441142134465829486753124219972074728618336792693395444136018185206595637403910837234716225991897082073588526041125387218156679554733229714450665816054626040522774174219519849"),
        ("fresnel_c", "1p1000", "0.49999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999991764058319641924011563980217584021142650648962420126269490146691416321668034591572607972240679963319151273272954781703533257104019512708154546273258094982990699690350120158863273718153025233464484865938320670592673375014012656855239479406871476610883334906563640750677336236704521491211119752227657530973567961381327722978715960539449612181114410261882043953659651860045481920806669387044862739343684940140422911434021831045981650184227090726879183338066061518798420002825300984720702470622179804693787830779152299751492369884269870432171592952108638958468325944801813916319686713362402291900165566179261235693140095974915418698431090306205497589430887948916075069483583184628539558044708196404496669186692933400127917600870713265783131923151522141453512210972112158158146462276758305151586003093803682519937026993907726483758525163065093429181366810026052573400589244611759255226218790679600523691796556167450099029298932276688748887276278362699865052288266847106645165396148323216776352184196191164061925998388063751828909039049713500514819218326316186267596567501111355705319633916135931881663770316380417498491065433228081922015133725000656873666578528166738496051162333389158051951920049377295356211738423595878516328646200739209265218111923926816883999028345687452182253201013766260394221391300002211147587261551270930865127160462719457073026068738557632167513849400753730380215117498101414184770134337012847388372432422618295320682422237622374779553099336821402475690915902107561612420936664811313521413053497979267481624236491728519273385014509717202772968907991485401968132480191224450610273068728079815611551067782262758007497328349695056007678189408785319168699813952424619807815216867"),
        ("fresnel_c", "3p500", "0.49999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999989301004755562968262349106189813490317159513030143692226938734633761008944169962802139355423080658328086550494658560294673271811462302192773664598168635523959195669882612104342425980565744023626904935755826342969588561859032430963992360858572002683027196860191713644979313376733702806677906284469042920369817294568948594439151695816704835634384871755735174262632130129858920976010119944778097431044421545517987779620575588291782310289571969046285068284520449272299947632724036640915680118807150535711645217478370970997941401685673834442182964753492841019613816863022666532013791336478124112096774283493244233744016165325479222868871789731662400431375755977188929217591761605428330825056076921869728591846744485988984880761250718549109160903902762643077452651619060082584794911270781935802262621708880191221748379931888331368142607474031767981661514234476566569200109242837855587080351259137136829711590909578427081772480976038019416961466834191835042399061264532838863822629753513196492522138509875928616789566952819762629058556164819886144925008496709621275313762128655499348899775103959789491668821508322727774211889311279670412868484958290068604401772092958622139823143935398735987351013646616265338207368897156911046727100054278893665381291758000630346878506558868283665421963695276240498607705421946124047819994326449997336674212319055933216229751974501858023768457344178726601118952729817261066315623126831629743275438285251719338733425042902640671902678824433190832150856484907059019831677669864567121847332762249141942934850117845962763170766783751844214958648162279797642767946429617109793379198206181977702381567997186517403060152564168942431132055365859016128959554815835211593026469771784798736141878313531922297376400160062947358345068567598609484848632135849013812460076466573046715874583273182710446031350846419608139925468892476009267354670814895571531000444529956254229581913466278630118189678789727329302197856487572412021051613045168455444916618528835876819910264411888106811674432037757051702626662303069401226656206834402241102425789645439754728706018081258118611011287172367698642448028871229084963358557300088608627284321761"),
        ("fresnel_c", "1p40", "0.49999999999999999999999999999999999992377441344627480657563373478046164914498482145569628916930013809654269111142050125530836060466945543600263867872410495229120940231025032694362823370400295188372971398084239431791604038602135743541272555567904736021616543979598656826442210340401157526608472617919284903533569973123219543306873786024463597633185956407716215973365993513895019840274534961340771085384287976578260318642365522035088462877578301497533873657339309024875824653351209113865608804960594624716821763743971948252168511290873161878099895200119349137572478188655476889140145753288007524942748411479155066228760365223290712162234146759765421970430344335676285618855686680113007748598254740987721488310258096983921233808724524040071261604679899006963071375380410307903155226170356619160086643820674839971442430887021268877784410661190661865418562668205896610559234133039317991284997037367140539800362379849056740384756592512463751726531672254406674098234760198256833142173299382234372515778435460608352956785627701936231845299124131382933068426811308914473650186350148053367498343419080308066061548929716319631735274466473310378750897054314236750912285172100830028980267837177059957192227047646952397859261551543867081211998419484545954621358205033598522793364811886330176384663587317100791090708798563934818793808637085047605581420138224189457719424801532446579543224359752378371253225555169652570017448896988189813291588769179752235852797473563497899170334957320164752006439262402035790538228962866126984282029504430234181311167886125194304223809751919318850200727869033375794956646112945552100567287034484730807865615136372006617930942461293624720119437281208767603173279668392664222333194062725825698122072045934751321161981725877314356694059706107326961017770248462717323191468424701711365635064427494029538619722172757145637268925345308468620990766633437385190324002606805194613801536264910334416371149764320719073912115485800152295120615539572995260594069543261246961803387653024732947868004749107175738363756319645776092921820215114414977114412047210805270786967391663808434043920216194915483847381774282798375051455781422224624327706926510114896862758912603737711206678576416946164537412819620718056985936607676897028497238474032960633754924769195264353006053220824091744488575584873572780924293303314344859482507612896208956852961412766668408482025021525094265426033105962342031532421523740140131691989242466793398833280448071180201589451539673063999575306044022094523644627888865225982763385628573652189063497149226181352541211512841940117380896354637909534654483235843084482956122609"),
        ("fresnel_c", "5p20", "0.49999999999999999999929694305648751765200427472214725397351542695031570125119617006746603846665164939456114369907716588485164600823234725132671654702665912407128382982478578554251418132757138759240099054954624401024842740356743889569544367748698315804072390495218436490056474562018914719056893363477535142071721613649341730093813880000838905569231743794248076945138031985908926669547045627111206029859535870706831184707625309756382972324070132608629263714046105272697055342010588500337878943286507349363728731292643614683677250534870932913964592004933898553449317014740423943045799994065082663657279176429934084707196639303525122954472950795791708917246397392040904490198891843749472312318860070637949523233021075819940326677062323440772864929788997470667643837717550852908125245836456871082146562209587277705611142007619693313141507040677678552771492627383183374469179987510684625723486474430790761532174791503824296625097265498752334868029015589353930319998808341803504706988130868420077640445995479687791214418965974621766447255948285945338227332715036745651886450623907629614656574884709364923295384550445648592803117572668133451890385018554743957633386003436760406122497063287851340906657621154180305386347675724676496573363454942168262140397663052976692245133743645549728922006181857101715454771063061994603327052284132626531043790674408758600168976950407477601635396320291422978353454557747761564638180373968742460321688741411149750654458193450237139902248812894742847153526179498038065363874500833565714396686190095587076653582913597690119923420251710337700210529135549032783182786476114137295433604614529506367035817201131706016154598085872761673548747894430245879598741191007273850582356275704581187107744380793009424440773743368801787913259164737748327901774263455533806481210914385046158706971505670907845251926017317599412324425755045477598631826567210706259157919541705833489715486300377141831828411395730330743772121335451512454415257576167732828479363618392135318675719135915610843116551836078902197895196284423046694155777963301620082359099592686395416826842356650084968478195851601546067992781138716410091522081110978958678565005874906275586197940439069788903153139797496467320733853876561635140412649606810073869257444381314996337392347393819298483289232511566721578023207707305077401987805906306237740132861226081625754875142472717093384459559055490026529818713253201863353763746746005155539147504025827627397631037568198730518876362284095504326999272320104195493765108052490444301109414981571133761644930850688026704527217031448081094428133681826290847663998991868522329799374191"),
    ];
    for &(f, xs, want) in rows {
        let x = parse_real(xs).expect("x");
        let ps: &[usize] = if xs == "19p2" || xs == "c9p-1" {
            &[64, 256, 1024, 2048, 8192]
        } else {
            &[64, 256, 1024, 2048]
        };
        for &p in ps {
            for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
                let v = match f {
                    "si" => x.si(p, rm, &mut cc),
                    "ci" => x.ci(p, rm, &mut cc),
                    "fresnel_s" => x.fresnel_s(p, rm, &mut cc),
                    _ => x.fresnel_c(p, rm, &mut cc),
                };
                let w = cr_ref(want, p, rm, &mut cc);
                assert_eq!(
                    v.compare(&w),
                    Some(0),
                    "{f}({xs}) p={p} {rm:?}: {v} {:?} want {w}",
                    v.err()
                );
            }
        }
    }
}

/// ZF-31: lambert_w0 / lambert_wm1 ran a fixed number of Halley steps at a fixed width and
/// returned whatever they reached (W0(1 + 2^-1023) at p=1024 came back as 565.3, W-1(-minnorm)
/// with wrong digits). Now Newton on w e^w = x or, away from the branch point, on
/// w + ln|w| = ln|x| (no e^w underflow at |w| ~ 2^62), with a certified sign-change bracket
/// and Ziv rounding. Domain: x < -1/e is NaN. References: FLINT arb lambertw at 12000 bits
/// (zf31/refs.py), Harper's inputs plus branch-point, large and tiny arguments.
#[test]
fn zf_31_lambert_w_correctly_rounded() {
    let mut cc = Consts::new().expect("Consts");
    let rows: &[(bool, &str, &str)] = &[
        (false, "8000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001p-1023", "0.56714329040978387299996866221035554975381578718651250813513107922304579308668456669321944696175229455763802497286678978545235846594007299560851643928999461431157149295980359437669847463560613422684613569895704539776248557078658773370635663330123843045563542978608509015429081920856055752374819658465950807273491673107401135829275145251337631376134005476935791418935544554347458745675020103267896288915054076515129277387382563226194000274603792935151282853760427856580899140888750935776569099504280560319445860923973670113242094222763346981308395042505553693595498589171136894826320132290482894200010093532131414446491206635212092255780885809620421700993965418704378399093547794413996726326776966985217810284138501695671362843197511272629189359937851235877726884545245685013835052647838942536672975613501824051070968391977921071052269899406233099095198632839291854733459761743556489332235974900461809859632577244694823155799635407673027170995056770568075566644544768766740420399488984871770877493798243761073164896088287964733623685394003282530578007973591561928800494680466153579977739174343002480486867962349369662610623915243966257690474812402170859616356546484040069390416971799433852823746492023531244069372750762432289393725626545064169283186105987854281976795196719269470945521836395847827079567731737059502859296936418834200419496985191002500539375434328348127049608042474116607674054887739628014991156481105894085515907083780674940745643609591724507554731120561394224924025996708415203023031540122593308271305832155483274843890729390094416212626756439699727260333330772011244136615336237859816954564100440801312554375891875524244391889391154148450745821914145986959756881392986954557711822440030432987372899602046947962605107786472150564124534156821071103978941983008681248366810361501986448492932648303248516219715733907583315265662854935871539421941723137455050355311314939264405624554593165601494995893925166729623891924002273201068695261439517391371309597817033500861749531569432764100581446220841416846513690638473177543715772828591571366965691492469633487578314075511837249118510935743370183515673163588829511853557693740807026114918386259193250791163591901738330244978747318807198536434893899830012331732862923280115323416384885107989838007208172446848014347144756755176359287733786648703227561012502858262326851324147316481889375174665431899343599842502133531187971114417537805573029614747968776115929806922727054377090113639133819105922765289972132032362249875599392463098400216538231987784924887733737620686939304721899619035534805550831500895639612551946516279747224077024162564102"),
        (false, "80000000000000000000000000000001p-127", "0.56714329040978387299996866221035554975594282221280273345380271004084211859501131049893175603237409297715198466456034776701485475036171632823990406962852309606765455016086358862917158699877822435215346375284910096855521413995518421793871365452815409631803442360638303365231429737776625796976512258855331702056121697869231895153626333309094713591324935867162277648265883950745054501452266418986743842883748765465404760086564191085646387301676456670434617260935002667702250939228186943537725790481176944071424905362770132947610080760445124340637139517762757003030774036419992670985510294536292349826480916352173317852164843060926805170017055397948291476698628233768268884877665018880213084596851499235972057786731730129599397281890257470634205346940593179761205331818404266342760415706379517817503483172499659163162948821826168788187080314181217750929554482274827153561119590875500290535501996940197961792314815412376695702499767935218795939470672238947754204435321289331143946317290913410768318695511993862092930545083227090532798644628709287635067981577046121752299244695707105705371937847976400269054555176648030148964581997231560534487769826846661025535200057462188564764586076320914851177613592906515695100281878385094647492497041489394498082650194360446071706687848367770658265550662801200252112315327366420712776397262746585612137936374613765198963629156981534947144504244836164109503505183920636220255681150690758271769480679830089857040725328169546829842360450354946217530255105936261601681170410006443275983322972373566403809916000826947020210240456007473577054609248955852020499789964207374464275534107416131297453816130240384997636893219863886224690264898526445134863510486247607037625820804812333405712655434809901082910445186426000756335903345950310920179128372830093274436181341439592496740749969456433301895529416774098038097902653399477410743123343083929353804908359498968735650424100136458608639721606419467449082194065382738416450372489393930444717099733662311657043471704319074989174453976748941234101261149987638870757190739687423921754420824652526081891470093655174844867044138710078270128893853780761622713270602196478718903372709229462292190397066002511994219550978363506263533185863239636281389293491921037047630975670846946124767810234792900363641667603701464857047907349410761421819823692965212849996997084233519339931895802913813566307584390509682323266574916947894506224628241650616442095576747843213683564786137642324252022613627779755980616857104429469127506644123574892145384437125909784167276264494960802763625502400528576239261058667974860283446388402581136334186950478"),
        (false, "80000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001p-2047", "0.56714329040978387299996866221035554975381578718651250813513107922304579308668456669321944696175229455763802497286678978545235846594007299560851643928999461431157149295980359437669847463560613422684613569895704539776248557078658773370635663330123843045563542978608509015429081920856055752374819658465950807273089050157336183159607066710803928391836014949964634934844831746591593363689336809714908569837175100935467921667475528897314755889250305728224604865124854109688318448770433467727016574464765200627013360494805788387577491463598303480868698562734209915119830613025027022372928387272165424265726984306938906858750536086180795696303059483952700913335139521528955267173975273488596948880117740676434178669094013112884382343781242420385441615419639759208655839555284999764768459844462100298341434660635607825033266960209477398558223930632653108452661499473519766025880693020451390468870924311740477296845417909315878564467475853505228397021295584218675199653102345532968435066823936754136902543027660760430897475160799435113097998178729556074649831837270384518263216932134516607803802713837874081574862658075959448428264392732080241596571607976590109980150299684423003671890547381582913014312543512866261665385793351405061250410222623332443077431838506126389354389267411460651858391854366503283637333277952037935541638664028357521564089539980593453223916175449730618196291048913653742217614007704052396871978423532945146320143910359403152666349648946603633664935871268925515383300671862484297058273902750444220452227299182862863004305340518059045076439262197472941702144641268361730818922388475020167000487434949049726602870728770888475673828956849450304818478007035670452760577649415316224143784503462221018823413421278867350277199801771496324352256522738538417048100782971656113791531058252905321179487036654662534155771957379870294747368948863569290797633261367267962640709870549301693230307738550154562864105649888679984247670106577862342105631744353235712046123019751412449303118631592059512947407631160874961077178083943383368022683729218325716957069186836979313168855445840400665330920982944380088346962437549053138453899854068157103851120275775184843965950047140043759372945621233152875392167010195263897313569151703623747809567024612452592375354857346125812002259698575182464732481510379764652808276195071444321465528082100334000885956908234630896187377917124528451194194978642191196321989086480988487153521950895288571569938779646443805097310656955499542364171897851161492492572603735539844504571105726429484274975255784834540221485805854660871076106684328163323474344200580306310494394837"),
        (false, "8000000000000000000000000000000000000000000000000000000000000001p-255", "0.56714329040978387299996866221035554975381578718651250813513107922304579308669081748736296020088593419992123261565657734616247767412453059124613100016290431540458914251214413487609267937388560744587208924515562192694032369016329792451593879170011358295274630327176056478508497319190072499525093679086722424252507281276274546613213073813610978924536318520100284775165218317288908869014685257521763951245620084845809227582673401588868972004066681611116843908132123747505796015069100771389173165950790532420844008394652437440527450469886312484529778172033455694650841148668469122493638444671535722765723481191290398525371327074152798053055699539303780719889012983368803518945962729036480561676066578247744496341635194187854662084304701683626703236295361282242519219401074479558892530199533317470683943826656233447728004026351161698019586943920489991973816964757939791732707943596878855119446625455768367622462668796426366298828326666143504764506998643189909423510991846258231047897136487462058191663586046296284850783749388962802872447658296935456832023142271955346475084227574212872514651051273102071027375215070792373601935166755415061589216858516706811657429007210465491637461081346115122255816955687033323440360754700542074628483572695687951995166371151658001406486677037613771992820647123303213068200189050915492156636502052483583244556236760337179988765032483533634104222912178762378541684619981724866689136076950046062507365392426501285155309322848959691231159663701283926742317808011778800985895925555393901690233728938612290028841650522537001462319932184786361212707847920831402664143009962859392459260085453724602666509828961193631655420808780786939119443219444752787529215802989639355353978440821010434361244828716108718857848655303752710904720377462001558129604619385331021888759520194350954583598139159064716272478874874576628452045323614581784508349822762320105218841233353320590399441979624856796005199821197667463253237034243666769563791445395212361007999793314411874545437258599767927045218234292482129769698595285474562975100949095178949735181389021436247476773873975862635265672470961837532044073838643186188837381120601158763978761676932306525738860534623626343336968243836942953522861879317661456285768005602231717064875925025985800310713764594151730332526773256762598477637913766451730308363078763278700180356923530831816313203867879989888836344046444354108462546186055486527016405206205700922925568263571582998290376107992781440175991402556872188267763543734839460830862284790721846496203962465419632660899388793991503509994705410867014788986886867126513279603390476698447528889048"),
        (false, "8000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001p-4095", "0.56714329040978387299996866221035554975381578718651250813513107922304579308668456669321944696175229455763802497286678978545235846594007299560851643928999461431157149295980359437669847463560613422684613569895704539776248557078658773370635663330123843045563542978608509015429081920856055752374819658465950807273089050157336183159607066710803928391836014949964634934844831746591593363689336809714908569837175100935467921667475528897314755889250305728224604865124854109688318448770433467727016574464765200627013360494805788387577491463598303480868698562734209915119830613025027022372928387272165424265726984306938906858748296421678234255042003071966795220895590936913343950051154949542716768789494702443830337878400260663760909856364582878781879533830423747555569697542866561354805400901104771237324730168088420093342591933743019354662353076567270975761218841385994442828002635250684447525596225061138784812897842769388047292026888923851648475342384495390278960997110605478422120613619831119733762279760964910117710881374070497324611062417462642371228844836331612237767232395672044941712482078124546406074239932977793693621417342171703795002833882744532085064816371022285818135905329253088215908316499482150654308765887848700093471483056995335342309725618741850072849647633657168786823359704003842175324788340540840078085163212011118484036391277940331803756388128964477481103416897999080258084948831702013776155886748229740129563028166727978856884278039470309882476223349829199167264039162857192317839390367382047052276521546236550425903445810119348848043130419046747599458626778005337011116687395474138280349091563337124683212099123257109328623501176256846187564067662282083228063832689665671723935108070720191011916609312692668662866913299818170392069762444719530732883772309776724923432002239640853733202610837676218838853278050859805575727041311106618298213707609457136235519040350889999051660006639249223109032376618008247248533910644834391144725172054179690050993846858110052030727770984836749910348642622373124418174549359690086713464258024171199284247390659072053872420519135320693925041876944485419858933242963738898653729424084111216012939896787033532871041022301336030881194671213286480331660915044448417543078580374489476935416325583525550040194016183353217181115178804738405357374064031756042146567896511182725504319420934456731637817008513731320675022884570970280676722606726694292234089792194564581917986257293733587857874671744765764929302778294026787699592468513485285331390754346946020643238664225034776447504913735592192463194989620115437134611256574516629410750823962929564588940581896"),
        (false, "8000000000000001p-63", "0.56714329040978387303920553297614740396778810447917647661979322306171582112162323587884962778675532372780204399835547421125387333481477943386323777109929848458276581349302400579205748674485917163869116375314453584158778168306332963931796275233822612467302449263703873044280482736216827757618430435488651386709307389781659991113619171983062798859112647213148842064428761354406722753697361606286917466493948523172542257249613680832551627984719153450611980213353382779670141653529911863986168285902621095907611601606366699250907985645597891558529228961516095205290486173906898990799749236613458134945583483467002066907223847442166611247532696300705982585500700651081970624790179730945517804250637251878821093973282904757465105763632986454981754781530759470227232443021880781508882121692596366429914178757925544219571165431571620132616257955789299822801863274003730746881575945379797177881208467737071922032395533609938412695105689794431028999715214427023826457294106703523078500295894631728126792313036714814826934674760118926034167176097747094695293291423545616930588660825395798752452634779388350169379183589466837291469834438547135631120575580854677261291583200346790968208280926127184037115924853439167966280561914361306226265546464901016701987696939806123211820064156109411659608229332365724326031226650066545047285581175573480733771925827866271925774113880757263858533675010987471116478067728504136966016019384367239089285429303100702305039620907179811394542629104479614855396242671167266416711380832492177442786882690238572800940434420430559620148583677099449038841312570643494653483872720693888941737018437711703455681801661298622260705277190001109423610061400020242506224368694721541092424603866238496386551210209764103521425049269403157040264941533252593738517169479523478275968287701777795951012840859905126393453305743012945053363917033040012252004031384826610340865127132692261727871183400163316360752019071024372609245209643335849980571007670135571876577119574715595394740327538733132390577066865843648081104981719506155256567694507924161110214409800799310118562385341097966546412052701897308382388218231015254254811960603282944529871434339189649033959334150227079881665391911776750930974528095869048914200361121090647983404109841066616321841095005094148560710891112013833703597745460256943233205068155903109093049017765162551695697712517345246018989265206363578540713936453711095646292140914644231662562409813014819282414348119626574364925144582560596801004434740354028163430946457066679897415231723889011294169624188142383659920401857359780112982942702625801877092141692082991937257817660"),
        (true, "-d385f359df8592a9b05c8ae937ec20d4a1dadcb75bbdec2b33c2c28b61742af6c4b323b1e954d5141bbb3126e8c281fd21f4dc5d8e3a488378a50a01cf2a36da2d1f172a0773006930c5e008b39cb5f39aeecfcf092b18bfe96c314d09a39344bf61a861903e72ff1711e99ca553684052fe2d1dbce06f787a160213b8a369c3p-1962", "-656.85035413527969071366505447847441484157595423067446791400290354825604875373877052318475930441804934961263368095717184010454761049052789174403128229404105675448965221041460723610392877284154940386761087363099142866802731408102612295431066326895552204979745540105316287778798210144751790865673153026969242703258541816850894173352018245235240142219417760066618464076225721371652308610856922238344726198610745164531022074128448951289119166515256437551787423287717324422043498542371802275177006117388881260899721769713013472357676575103263900925572409733437426541077392370161777300924718322301817758640764156891573004453315581115441925373374046841952029870735883375768587689865111970268494120151995788850107990132743969625325208627449059666511151732718366154822822255781655165216806752998532271784755489135255794057815822884956382516381622103355532006894060227313686563338180325058776719854383675480448268056822245356187803768569577693555803359781944689780273360661573243517939280863725520032570576190113135517432406010481451476551875206831125387282707947067355871367697293090028938067812688018934836365144765375412498116605061999326319595036566907450989341676584514786674571202794216676628416453598761895197139572858243278992797186655930037857571571211074024872683275537865032803623049282837884486073352946535715647915295976588254996195910528766611283409781617264014728581498728811335436763324090876118367182944759683144695476007828360172759373367181047609175689242881839758502949936098601961809400969597104037014119702693185313204900043921333354831446535132300834327064938554829347078024112414204905224982630888646365385438332188329288000580017830452085754416929935913958897217563026507611129646567375496643653986831913649523462305831713621086017303309648583827051337325754965637946812201222015961173816485214087593683613574974410328893449635086161790876956071921687989999089938454658136584647700323056005613234647757832159717072523860190817120938130513169460758662715230916862059355231097126825852441813530462014538840573814304415179810099781293509010675079591568215398131819868977513235138903119376423043456631162207155572646483480755403544126493924585774867336546422182990794557259005222709438332673858296860264870370069103491246239804223500283094279620684943681879243151012247561518483820498357006051510952040247473138856165874673598396419598178398717332802956203280836327810556740344014745835909872266221940219506014357241951023249164997125637411106165314135304998752477869261416910904097160703396683644856127896033044941839666301070342189738488252085280695101827316183385745685004671765857932487"),
        (true, "-minnorm:1024", "-3196577161300663957.5557699490705578188258940260157168051065837503912724646154612146184526718068209301974001429672739325674156349839781916355198926775934635148876107462778819236244174501524145001429092582594505782620500131271075563850510998948221280374603279142621900733169659329401935990285754396226443032298765303548195884616688982039975196253039306436382699834768416426133570764091442810149516949571596266710320393414227754678991260710257269320069268161745176518167391287551921864262222014246129334376846066214849807832818277358004100356798902195428159545599296940955610231529960301681916044881002704471398877277208374621827144349283625389823007727274179441369261284936035144081118216410572151452985863879962362832353498884145324319115666637216742426591044275584909692199372933123992256876114783859155693094610888215850972975171984831683023101028620172049436665240465276249776153454419353545691289630367273670132371979198677290186592968539777600839338879868678747796146038936264301704409374345885620955809146207567035535497397681427671856622162343109047074483384150272937085716551072432406205505418199321258560188105420105386618582141863069824002939567342721867383135025811130966721696685580407138228120613873599178787108921143116150061114582200881150496037080305178532324802422449957465610521060150462081386239889736455440148620446843037825602897175019267146187796510232011660092706065678362331493623834260838779167642640142825215472906911618850475209860199082810761825360445381208017390606943367047913941719071970629436890372924192544836655509571249388642779826647250137734571232138816027315861801784829666199323990341042234338076126887341253172842787123650829463069603599065353739904872810547234189021466299371812406953616280129404033244141606199620163906937230536888152921002185379144219988745926973434873374131052287447385560642504882599112874788088905116280415007668929834147013519824281300855811198855061783291194222298440649025175511557073637258513239508994428466605532446054630634496208912281247320947114743539224836045198962266597857086348048026136860114274099358118122865687111126820721019672339218507255029409952747129753933757566493022435298368688757074863552571633132566042874053621464596496856177086056534976073519622076284491232581539401214620140926947093418641356127051294035241438100782298271593231397220964129771984184657746086040344296836078637663708645566544569496019975776205297087460700346539480834415849643791772889781030611962589739658543219209837721703314904192203706414841277714698336404619854867654887432943785358122536271758488176174685634713511235271826575843545280887"),
        (true, "-1p-60", "-45.404440747686248841183027384535827451763543917472449854274465521683695332231817048339725879104112289594522547796300830188961101282778384267814462130191199740413229504769988876294655181554723121533948557503095538830149500741368760869603970897285617507841559777590432351072332502769622201467898102208242162510040366181874395419207329350104743794942001084045013543849085182655981305712873079557733591652650949038353350448492436995203744488476762689044389853690197044464862263938250376327313478727338853682980451906538790008214767718191587873660467894469699842415832197300404451144684782225050210978575459037808459334288898911404005212510018525368864990916016891820285448135898799847396172108180552524693354467695751022368610394819913808197321709170965943824358684971023986795310718861875002670763177705929135350808927378701482619081589539773163666540778335162877880335190380352658232770734020251294348661096380255828720341543016575365434426748904760439999030158686034367924415866120039899087630075069579118613739898925864921414528873481111937121840829979197377121248957291371578328373625605119453668870520734493661262727676877773795210838861137190943095379389835480731246293757233920192087905361713523787093176323114558089128513721775513070084790644447370392842845327424264643823195395906851116043530844595693768699390907957553221530885472630103717380682299062807987058062598228670722533921126443919599869603658966909609300567693389754148713339762628039756460759883128993932964997166608464260795366525279324573579731497038004514157519688149727276837402695379685515769942827399805551364385467901534762006830381438906231332684865409934037112307146080611726535580848093125816392564905491143554180680502482288585160973793674238569237823114595053496232003430909148264646656324705022717305355919980571252979736376826253566124662046104888431913031293104681037768431787359900090865384077008487114746844581712779352603514854599237582248906932622097154564586922017075955679844846922109779716316627882115693510590410735836548527540673430774374329562030571445755909498343978095222055880078195421790512678417649976657642287006366515247103543849447898057261591454051817570775125346525605759092925934303331915087485442284887775517661700155149536639305578312059913294018624263334987252052375456641233966555166856107682862343291376278486747290165414362372293771781481893506078961876672216868423163425039218994915784638458131156859994432573363523098607464794533179082441927270551773691818139344150945531703817207444700506882035598206202124760060556113360014398787876562622786253668568260034004192562581083338569678105558"),
        (true, "-8e394f1cf32fe5f532943f667d1860bfp-182", "-41.749408150699907593154312188761938139192710438571812414334892520328574145733076397152000639779821201057015207705572345565855427875388917756957280113091151483800985922551193803790475297843409638652505206284457624380099295157607933523730987297245587002931501510481650794696226264669676857843170379883281990657102176851705004182427949104192653249406847756493685904355567081263793448182088558463446316391575251468510333640735032778807418094060773671716143942679522742142005494744557223664982197691270287748440924175501381141531637117711344531556967266786890053148239601222996100151218283212325883504337585262064928390269214660457878409744020440013366275367252586864746986316690015212366726032974915689825082552927767640094807003203579945974172059740363316345721759423894308269475962495063541451000573123439062288832287709892431776969301135071222230810198546330780808305018370013467036044884621045092165087197310551304384161399254883786047725740899783875434855271010476187996510573706395846521821270746194094242043840419718342301023948184544628205425157320761554516572075295897278087691564583314002497836150362529350783382199695720019215068454646273777451143906309737339689722841547739693097619875242661183892362086001876325141097471062702773463756962726310694158152272282418677675777284485954463843756223369415969019684439816517251830596220926133981445949788727797160594220914608587717543330063798655494427613751079615921526889350019477137214105458509348528950583690811920349690155314167326626666644941224730989425775177637635111564233933079711650747484722349984741985856324890327295814553651204640890103948774782978919518365672769454976971304487797096757708083631336552628081836598267688533638649400195663706789900844895561223041382738400583928509367864761629267345729333341486266499611975118676111160896692686677492042830113498905369725736909499447560455007784382509146101061232155239647564820406136421281637909776623862889930003860696067565730304681821884797284477180067545226978755200781345975654641279429649548536560358856280918331520457016476968715034176373119703660084357601715624517880275936414830544988956267170512255298321380698546548238676166949340070152396281587112567441743261345009110142441143774446206132615912363927323675481946609379086544247524317295465740530359399476901909500999503603762226739998086867856111143180602818935141305037506747428445955110543072912555079329420093585561467399755310702219163635174286597653240250092278725775155146536458819928423825821975086670351853418704771040611309735836212767012886899991361261127997217436362376445573988569754655132511314240830225910376"),
        (true, "-cf8d0ff327b106d1355b7ddf6394a63bp-3063", "-2042.2185672055165218125880847038528122158342945757651984118049561985409193263151505243950450617490093639955598139541134609478262833620707392205739861836050357264191159985816622121362860965565931444996223340947582794043893825933292420540180289019036215085660244850174354522497833065901657015868224549059280584375369965731869858095669853181928316727380256039250132248181584147419059598219377176752518975699325196910631683588714100873190101475655324152726192982336945391896138574631887274542339147374812756588219278345988318507737610021617554490831200246208225803938710574997078670541593143079726654309585239360650074566954840796298136165379108544307578636080796967070338376780778788901838353249632539620194075789134326522193567103098165544824367618346967859416470252619713509540726693987791054624528655301584276612089856216840817851428845570186962522377097131497190103745006771403575505441926320449018171659553840047744225062417804094318281229608669525813516995021410090144648610457490501607116186719819983511068932032945386320376257960595083092589679949059111992078124355658845246891320867713453229069724690825212369284767960310594475628880198907923191007692871248864949288759435610824080000736813730974485699475935564330111823836937198279411124160912720768815874681890644954863113323911609876184663064465594058305124294859001066156015857830143124698872070259541661543425991902268713062656114043464981959909292388241859446541798694314600979231802887170146931616695857817744941136380995453334284384807104943183439797990763144173278044590196969503571115831267500987816507843104688867707970351771735520595091700220674279689226824354750188088185575596318516406194136429852632647874199424545593162653667523224944416149154513900802951388659952398889863803418861572047335881578246200171093329117368064111441068535038977482564825697811110855063812528239968595001650063003638738233796420381090300109166163799328809513022570670479219623330352234212379920698122955561349819125436303688033614863298935783502396845729067086453160306323675055923895778963109014218344890779090832574517513832355351651319649981887009012067289835483010172977862939420755093820405210992410958195724690644034905561361014513983160049285979275797490826488528332400027670319453633242764764998008101938651523487139765306782017025858699466577225425654856073757922666179152562916744288974511172383869602693275396562527429464952130837091414711838356132781591470632623157215407269739641444735454961146190651386273077197234225748297097934542559670236886877285253934376830614435729999519313974112303605997447992482598920792591451270707075664792335"),
        (true, "-minnorm:128", "-3196577161300663957.5557699490705578188258940260157168051065837503912724646154612146184526718068209301974001429672739325674156349839781916355198926775934635148876107462778819236244174501524145001429092582594505782620500131271075563850510998948221280374603279142621900733169659329401935990285754396226443032298765303548195884616688982039975196253039306436382699834768416426133570764091442810149516949571596266710320393414227754678991260710257269320069268161745176518167391287551921864262222014246129334376846066214849807832818277358004100356798902195428159545599296940955610231529960301681916044881002704471398877277208374621827144349283625389823007727274179441369261284936035144081118216410572151452985863879962362832353498884145324319115666637216742426591044275584909692199372933123992256876114783859155693094610888215850972975171984831683023101028620172049436665240465276249776153454419353545691289630367273670132371979198677290186592968539777600839338879868678747796146038936264301704409374345885620955809146207567035535497397681427671856622162343109047074483384150272937085716551072432406205505418199321258560188105420105386618582141863069824002939567342721867383135025811130966721696685580407138228120613873599178787108921143116150061114582200881150496037080305178532324802422449957465610521060150462081386239889736455440148620446843037825602897175019267146187796510232011660092706065678362331493623834260838779167642640142825215472906911618850475209860199082810761825360445381208017390606943367047913941719071970629436890372924192544836655509571249388642779826647250137734571232138816027315861801784829666199323990341042234338076126887341253172842787123650829463069603599065353739904872810547234189021466299371812406953616280129404033244141606199620163906937230536888152921002185379144219988745926973434873374131052287447385560642504882599112874788088905116280415007668929834147013519824281300855811198855061783291194222298440649025175511557073637258513239508994428466605532446054630634496208912281247320947114743539224836045198962266597857086348048026136860114274099358118122865687111126820721019672339218507255029409952747129753933757566493022435298368688757074863552571633132566042874053621464596496856177086056534976073519622076284491232581539401214620140926947093418641356127051294035241438100782298271593231397220964129771984184657746086040344296836078637663708645566544569496019975776205297087460700346539480834415849643791772889781030611962589739658543219209837721703314904192203706414841277714698336404619854867654887432943785358122536271758488176174685634713511235271826575843545280887"),
        (true, "-8e760751e769398ae1b708d712d8c203d0f083db590983fbcea296275b01ea09065ce6126894fd44d8237994443304de8cefa1fc1aa39585444a5b422ef600e74ccc416f4bcd35ff42ec9fd06a3b165d38162a00e8054f1356d917d7d06bc345da8c878bdb99fe418dcf4179d080cc247559078b8fd760dffa78a48faaf55508f8f0fa67f17ec2e1a82ff3d318681e04239ec34efc9384d981444f5cb1ab245acd82c18995c338e049bca2aa909418ed444dcb35d67b2e63bf57ea8c290d11f3a53d0fe1162821c9fbaf793b1dd182345818408ce0db63c74d6f6a9b889ffc3ed234a0cb2bf5455a115a4aa4dc204b120ee661fd0520468142d1ec13be0f4ec7p-4800", "-1915.6849802517937807020428910932420688389662343484903381780865681305901440409892608659065344619097905515109011009806667337110359186588611948288956908194538873663058355422263874471701293364109432482557584927112343228717452993755568920877263171105589744899512719510691742844914295913579781759686293001139064240888173505926369474099654265213575955349988167548896794021225745341691312586414945926831355578300381014815328297944287440787607198717777209155235795007145547924347435385359429285986502976832418953896050085762700748781824889533511354974325535030318312028227420574706191435029360731024326850381827315848275614463065419403995171632296999100907218438069294859395861986716490330361241371685195986727746659437151617108454204874167741233677813656558049618777439790811375355501509185904476471869529095663113566541481080353037807221017277032042069019618672502075595933099683072093489080693822079882940218223766330888063906885621535838260541030979964147582872580200810791054719791623705520514070864261559930816657914852511391087275006023445480102678384197973671396179958073543483738053440876202627835032048935862366075817838247191734334444131506274090824196723984397266934812459324104508817306218936480190389308408238419309660458101050300027833738226041837164437553838068767495710108455384181436828047787081016928039129728977819100955810276763454091969230969848256324972932865226882263318302900056754833888208418374453148565891475831778478409131605906863255186785249780988339680193709072000000353194317333299376674823829478104141519955841754419760621250347262047152302569371416290148754199450383184771298159532147128839375309069207675671459561343426559209704876369615072797253068031392037548245851825435114722825992866064201518431920962405551955147125171976507985674106090132915561946325964040836505776479373483694427915415396037062337362944078414612113766590152298559845041793386193703797050547726024115036356768283645106947364545799349347124557367219910679349750331593623865015490042335758628992684047408685901712539078891730031207200410618629546210643291074310880916477694771504534042135154154153013451936985241581333020757201018008950005276882090443157705871208960382646803733655920365202810799918887033250922095700207482358316961646667326421751023919537409788568068653292983713164103802695173834500550408646737498635486592369175690495661484475947731326817621263703873727771068353290076758211734893549723704633428205756352427742300508306708541643105370155872414198317631649809475191196384602848352259191156546349118901561910772907675294985020492844842393051884952764535851654170870486565289453714666"),
        (true, "-minnorm:2048", "-3196577161300663957.5557699490705578188258940260157168051065837503912724646154612146184526718068209301974001429672739325674156349839781916355198926775934635148876107462778819236244174501524145001429092582594505782620500131271075563850510998948221280374603279142621900733169659329401935990285754396226443032298765303548195884616688982039975196253039306436382699834768416426133570764091442810149516949571596266710320393414227754678991260710257269320069268161745176518167391287551921864262222014246129334376846066214849807832818277358004100356798902195428159545599296940955610231529960301681916044881002704471398877277208374621827144349283625389823007727274179441369261284936035144081118216410572151452985863879962362832353498884145324319115666637216742426591044275584909692199372933123992256876114783859155693094610888215850972975171984831683023101028620172049436665240465276249776153454419353545691289630367273670132371979198677290186592968539777600839338879868678747796146038936264301704409374345885620955809146207567035535497397681427671856622162343109047074483384150272937085716551072432406205505418199321258560188105420105386618582141863069824002939567342721867383135025811130966721696685580407138228120613873599178787108921143116150061114582200881150496037080305178532324802422449957465610521060150462081386239889736455440148620446843037825602897175019267146187796510232011660092706065678362331493623834260838779167642640142825215472906911618850475209860199082810761825360445381208017390606943367047913941719071970629436890372924192544836655509571249388642779826647250137734571232138816027315861801784829666199323990341042234338076126887341253172842787123650829463069603599065353739904872810547234189021466299371812406953616280129404033244141606199620163906937230536888152921002185379144219988745926973434873374131052287447385560642504882599112874788088905116280415007668929834147013519824281300855811198855061783291194222298440649025175511557073637258513239508994428466605532446054630634496208912281247320947114743539224836045198962266597857086348048026136860114274099358118122865687111126820721019672339218507255029409952747129753933757566493022435298368688757074863552571633132566042874053621464596496856177086056534976073519622076284491232581539401214620140926947093418641356127051294035241438100782298271593231397220964129771984184657746086040344296836078637663708645566544569496019975776205297087460700346539480834415849643791772889781030611962589739658543219209837721703314904192203706414841277714698336404619854867654887432943785358122536271758488176174685634713511235271826575843545280887"),
        (true, "-a5b61642ad87d98ba0db1f89d61c6345877aea1cf9b6b104b74d3750ef9015dbp-1516", "-880.58095984143926310317679208505892792879847470488041017075514470471694708662544044855520364545680678591533952026548553369970990289057715374719060799083525011801142217926629382038726477138152065491192257585768546383081027598294485086648762103836907021313035579720854929063202999436801744153597485906741160912263852630668424730184554402511308627854216447326314703107959947208321435519027613154923692205310943346187718404152057929969920773302757401022264506239011259428686584797952581887276749655851991697404473892249345884127032715407761169140002247308278844537126440918331643217078221855718071089015792243273815285152937765473856648959571526050302140793508288712890585460995629811211656578805213812535011650794436187021922494247215560744835478644649428398358679812887644290566221542583073338909878260459666828797257686600799754361756621144080477739525655796615185724370139504785939548314572725603178644287464242554040679271112728157622007317079102854571588345024877487404194334685024976405874224727894581775512131254506715254526424711406385179550017963719300520806392931918804592761072452418131001234896995952018395918020567674696237746147798685257689752968693368762304805375867701117690135703525083881770982342186765944032164180578763103414031564882153135027356019091421006305456490206774159016441922375098094265392152653982085496311912283726153751826937176462370988567652909247955630932069044946204449058693632934635614871196124826395526791648770860729112970960571703927957617942356755372668520728302125947496043497835357076520394625796281322271800743688746112147224194484676827968356769395020577614260748613627803785133047020026792425227616163601160653233025469079987283927355882266687107196086362150214835230776532049851068849094626355007210933063922541707974550282111314879721205664486996512117469729155353958652680246440467349730287112040251257132800712170399197343455916366415528458144670029554033126892419002931778793285404083033008424034642886034905209445397804989149865694954457652629915134685764292129400749054650669844130909642175431866760307881539596345912904901072714703639613232154590132857020257753244105849972963630831853059141006721688400737994719247401072378897902166858700023935918849769541331732577369586181301462764383385572164626047878478312127447985013564403364974742193577305812593368389282858175165721612924520855082954862850400153305254341632662646141525407872738806244036822590705722391680448332225303879550396043822008983967633306202545272205706130392775517784488995316964527572224186295458716984923195237440767027059251555160682529267403076110119373301962873722248399635"),
        (true, "-minnorm:256", "-3196577161300663957.5557699490705578188258940260157168051065837503912724646154612146184526718068209301974001429672739325674156349839781916355198926775934635148876107462778819236244174501524145001429092582594505782620500131271075563850510998948221280374603279142621900733169659329401935990285754396226443032298765303548195884616688982039975196253039306436382699834768416426133570764091442810149516949571596266710320393414227754678991260710257269320069268161745176518167391287551921864262222014246129334376846066214849807832818277358004100356798902195428159545599296940955610231529960301681916044881002704471398877277208374621827144349283625389823007727274179441369261284936035144081118216410572151452985863879962362832353498884145324319115666637216742426591044275584909692199372933123992256876114783859155693094610888215850972975171984831683023101028620172049436665240465276249776153454419353545691289630367273670132371979198677290186592968539777600839338879868678747796146038936264301704409374345885620955809146207567035535497397681427671856622162343109047074483384150272937085716551072432406205505418199321258560188105420105386618582141863069824002939567342721867383135025811130966721696685580407138228120613873599178787108921143116150061114582200881150496037080305178532324802422449957465610521060150462081386239889736455440148620446843037825602897175019267146187796510232011660092706065678362331493623834260838779167642640142825215472906911618850475209860199082810761825360445381208017390606943367047913941719071970629436890372924192544836655509571249388642779826647250137734571232138816027315861801784829666199323990341042234338076126887341253172842787123650829463069603599065353739904872810547234189021466299371812406953616280129404033244141606199620163906937230536888152921002185379144219988745926973434873374131052287447385560642504882599112874788088905116280415007668929834147013519824281300855811198855061783291194222298440649025175511557073637258513239508994428466605532446054630634496208912281247320947114743539224836045198962266597857086348048026136860114274099358118122865687111126820721019672339218507255029409952747129753933757566493022435298368688757074863552571633132566042874053621464596496856177086056534976073519622076284491232581539401214620140926947093418641356127051294035241438100782298271593231397220964129771984184657746086040344296836078637663708645566544569496019975776205297087460700346539480834415849643791772889781030611962589739658543219209837721703314904192203706414841277714698336404619854867654887432943785358122536271758488176174685634713511235271826575843545280887"),
        (true, "-c1c52405b81e054bp-117", "-40.722076062702048906270430584912549606308390772806828512798660398176376395930387113022907921409511987438104515401671247225398200588824453301465219347691044372508174265469950430533987025997775902358857466291665639849542036882125558357717877148665217486045739026037502471517221609001509362145602663633364676495554533605093661255572700664389779389588906289601756106754352404504591064896745347489451301864264120436345610496517482401138799115569641856207665799848342505555442429946591934488144405183561972328760643307389524173059015730349175051182491585344775652704920198988615115207420579072213459781004348186165643275216751216736662734837839769520375420692273496016696795917628665816638228489397089410575887403903639044335601355896101738116557095205508323810205298205827683215715961843832397204069606784660792538193809475791681084780377415630074728533095685815750411557392569575273028333671244918286117369786601895413339271464734884584025300353948713240034991974550039501292793048730047529778390574579860147127540250809156959796756900488688142720204816575761313448121378409256561893643413318528045756656795517568440262102679194810842996808560959020949852700751451536674775739761294663806150095999649524618955316922249655285858748542696255967789101766403441789353418158629988455602599131990238733580237544558054250384187804996310061533258846109300779249541354797264868423468672688991988315385264841535020366627594754943794606686662477989063065875118773063935314516530691242300912665894748515665336778832955641889482511789372755002410742186495246923805057563303829137463314230706795566251866951575397163478013903550433254592664189951878144779904167282916761905970701011020417084750575137866623734496181136189256976752237798173416893883131479046358973527248683099077158526663192429028134660083697439956971846798068724625549779383412737684685341353966098643169975257849928184143418602767880585791147292760898342440649313606707715220412101624171152109980019632511384458346105090078549447781188778173430402437907472406659635103710892255916330132776888215357814137367897641175658790679156570384894047393720292954526231938784583955554761611866545172065296524493727913011082691616953730623194271468904905519871358502852906668881407733158387764417796811912344996824412809478478937047089338251720174386812064048015397687291695832355165893844604334419433605122413183243593459720881020150038167291854296334164449640209589646288331009765045277709550392299895091328171006452584782489572312230995820359558802232354255483845111772940526038349404830108511530454042130465560601313887683551313342597855576018926269967187050"),
        (true, "-edb64f730c9c3175p-104", "-31.241759835962802121511060092652923350081084425881997172845931461561566332959879899127291283765459213363843526882198386662458338199364851737198763920612301755547274534439645676118910867329840768836115941502945067993976560807136249020701321378231357914670811748924365451201859060423091633384084973739472926689451252442834041136701178690853159717783647074721073846826527995898863832149683182506603724984147530059675768887925185169123631995611494216968455456099274291824230719202410013300374975012297006789173790084310806182943197561100649790780875714345977458821630612737334444971995213022954583497071558340014833250399145821711987622134443214472466456344520904870829959362844099575342476572622853633423133805624594552832154478614730739861827070652147580492884075276168755534935061553975339870289731194285295770543519584318200319264877444939882334924600786433668452133635145954004524168454038251630453486544413180956726588688137450830297727707344147680858107223717484792933567806441894566071831933850350847285920540671235950947945218604383151744961141365092182775738711083453792737046542747525732301547305184516528578965917088899001353263520688519327311948838069344478716706275537052018172818251328061160450384300754187522018720647747600322312325683271191719711951812567123233927723670705776010550388638938534277930815883843794985069531583003916083353274696252732129299459291253256873919150952729210051441722910986096266265090868945882805060102948548546904867314899143389020141192655360841485133559549067276504666201645272922860219187017057958626314573409384215594613238536487096791604780985816595102903023012833964595682455352987503610905070990133133014266883539418848528684729947419340631107639196148299308387710354603811811657025161290244633474319463620550382075121481938557518108447051316494786523128425974640346030832145408489327260776313262599445658912064436031076183622766443920727075036770867455188582077046835148275765338941472397861935544513772532307977738429045268870451062504659913409387789188000161071529096896158026536469722529769067846491283036017333884013189172813623823696763048332802757954586645179098887047131564353413303358269482451388290436194694327097615896845703209521041984939413414070000973670923905868606123281574882408417344929708563741474123080641817148961449587873767701963802158330465321852645053636405357804049615335195316792453062454656149455162095217422283395198475234593742431933872652145217494210605788040638365901868189101896586780159871918226262623926259269434607626224520024851004414595978389335557765105640813292765271561310628047215525821680948417498188078756559"),
        (true, "-minnorm:64", "-3196577161300663957.5557699490705578188258940260157168051065837503912724646154612146184526718068209301974001429672739325674156349839781916355198926775934635148876107462778819236244174501524145001429092582594505782620500131271075563850510998948221280374603279142621900733169659329401935990285754396226443032298765303548195884616688982039975196253039306436382699834768416426133570764091442810149516949571596266710320393414227754678991260710257269320069268161745176518167391287551921864262222014246129334376846066214849807832818277358004100356798902195428159545599296940955610231529960301681916044881002704471398877277208374621827144349283625389823007727274179441369261284936035144081118216410572151452985863879962362832353498884145324319115666637216742426591044275584909692199372933123992256876114783859155693094610888215850972975171984831683023101028620172049436665240465276249776153454419353545691289630367273670132371979198677290186592968539777600839338879868678747796146038936264301704409374345885620955809146207567035535497397681427671856622162343109047074483384150272937085716551072432406205505418199321258560188105420105386618582141863069824002939567342721867383135025811130966721696685580407138228120613873599178787108921143116150061114582200881150496037080305178532324802422449957465610521060150462081386239889736455440148620446843037825602897175019267146187796510232011660092706065678362331493623834260838779167642640142825215472906911618850475209860199082810761825360445381208017390606943367047913941719071970629436890372924192544836655509571249388642779826647250137734571232138816027315861801784829666199323990341042234338076126887341253172842787123650829463069603599065353739904872810547234189021466299371812406953616280129404033244141606199620163906937230536888152921002185379144219988745926973434873374131052287447385560642504882599112874788088905116280415007668929834147013519824281300855811198855061783291194222298440649025175511557073637258513239508994428466605532446054630634496208912281247320947114743539224836045198962266597857086348048026136860114274099358118122865687111126820721019672339218507255029409952747129753933757566493022435298368688757074863552571633132566042874053621464596496856177086056534976073519622076284491232581539401214620140926947093418641356127051294035241438100782298271593231397220964129771984184657746086040344296836078637663708645566544569496019975776205297087460700346539480834415849643791772889781030611962589739658543219209837721703314904192203706414841277714698336404619854867654887432943785358122536271758488176174685634713511235271826575843545280887"),
        (false, "1p0", "0.56714329040978387299996866221035554975381578718651250813513107922304579308668456669321944696175229455763802497286678978545235846594007299560851643928999461431157149295980359437669847463560613422684613569895704539776248557078658773370635663330123843045563542978608509015429081920856055752374819658465950807273089050157336183159607066710803928391836014949964634934844831746591593363689336809714908569837175100935467921667475528897314755889250305728224604865124854109688318448770433467727016574464765200627013360494805788387577491463598303480868698562734209915119830613025027022372928387272165424265726984306938906858748296421678234255042003071966795220895590936913343950051154949542716768789494702443830337878400260663760909856364582878781879533830423747555569697542866561354805400901104771237324730168088420093342591933743019354662353076567270975761218841385994442828002635250684447525596225061138784812897842769388047292026888923851648475342384495390278960997110605478422120613619831119733762279760964910117710881374070497324611062417462642371228844836331612237767232395672044941712482078124546406074239932977793693621417342171703795002833882744532085064816371022285818135905329253088215908316499482150654308765887848700093471483056926032372300084977503682648847604581227881174562546539470850477742619395054067469971432668150929603817669944369040799756951991590881882389593267075480745092942258708849833121628839553386481416637732768587338186519256094520476862847500709881960647767107724962878553636425319095807530050220726129673878630280023108131157320704583489987906679499474812472495381535295734596510849207375000428388927276165557885654309089540475364294256863051875056204927868135657180674533222413030170592923361501207889939844744448491377309663853803684013156977894881774023752332813105734765079119208750917929528078423226023202524922869567726053909014023702582698170077949887036355265587825209818299166866095332488483081765618884488773188457964962984491514098462564296246550150981947859601820197030939131637411708897700764940997739065072852681620192991789302852463799476672439204442427788817076826123536997010591684815488234400397880469438253411277676841477798315925506677563604914955255155743418581648067164571470078088616519231327726516026000550066256047009430468316665242618199824802593215082576630966945787425939719080728935849017970662441561515731510546872419921765800914112915395127731084724210731722651454509383087229392934382549884861325165475557739020757629517898705234902927829714925391507942324953167802811061027750240958128036664960665366674091684878637029239538072727726973866183"),
        (false, "-1p-2", "-0.35740295618138890306881110405590475331659055507601204362762044858967140259614579628961685134444118514972510003129829014363079137615153915679879519089416183922405311669076159806886582644984474335477565998189259225194441270708746581164752453594256723176777251814270218652834918878204005374522495385663211531867815424702132760136406363751632006987382406334148285805739423422387293283962506841656665793792012231794078512034316839731212925158279340379909231773775836488227613360022442438166405550417247044781042845651024048276011299898266677289948006712824426550776598051838210939593939327015962718596691035671149368217654139787436343609367318034685439372961704417955416072999448475676222632711788823840882795748534665612440472556784740202140467062398138311147489047771979822091118320470347004669181205568388210980344075681275855946167257124797939194638469606760831619911349760344107316001979654817647606441099736028531141012895824808305951804729506502207685524003074702627355991995949819807131853359895597651204334455021149150327487220262216313168113922721214682196176903090392702686555645925670416705600217092648765859623042290889245733015067574288702583449564265991792314456824243123978354104717638381909518932428401682322611107221648380394630491172478139398169415313067682498853295413082310417906912110511670985361927799359562574958363179508042099163398000028437640287415371164557876847882946892862292031032441370767595751580136461562158684568796758618596787479298666560805281260953064598573166532156628557263787751224232081397348718258249627304427566869096165848129948548824063485853379106768923760522506767407669220157251252143258451719918400093778150549493291985218008519058039870740843022037404770350562547962907424856501743861903510585482006970489760881511626483948431877703763855928601129275026899915201237217091817442618906203451221107336555618986463961965533055866890762603150158332219888941742520096053198413543933686323669706215782185941053059898036566044741606894660453804080405850587439800340855669708732974475962839458566457239546562192664008797016282958381938572224470866302284563131393743981495715369439473482376632341280546384481628041712622931745522039559250957624667786452682552656033360790131143902472635494981143139268156104251562090347748801625174275229116869339832065878573135782284809146589875120426573364381225623936343527834735502288894899954573112678931505415202640065972574455736235973088452599327673101484362458743126553169147053271698654594687859641540717246437709018229301913136632185205571513791526333549008631961308602858997725882501250360509743543448387583375445843077"),
        (false, "-bc5ab1b16779be35p-65", "-0.99999999973966637381919200439320102934123588046947929369807731288211641302036366054383237735483547370452685176612430093023310964760897051956036984357292103181234136868714156909820767750438691774051573992201246097255322686112432467289273446634356204427496180221222279588021122959809320189950476421067222057961091011840885873749660938321156495012986548098018253261065277962129240948668851990772934972323466137226394181884162266940399189080051923956116785816718421236223512718895971221859116615909843171855096268860948657329407691129407826844825497951577034152475289551353773081924434986833605454671128197415758982181684048509946134297234839073518381931300916182617126183819852379590266100698691794429804801416347208850303832530474128263360946263193680062469392763279329011484561900210719360018511400965102440751726584455587150531199504685895179356119786517609098137743272388076769055205602926542658241823607124962569923338371067168865127941939163001588376016142572209445366434131760967394239502570384550760133201527638811399735308769817342449335584407333246349524302848257412221433065284795259938915557436884247473227137943397142526420698072163491172949266952680357771702888095658434501558494937830979444687785978356644699930830817305457171573384678295934789396853847426483909724131907365682511671364102026971105553901652254633493197094849552113171166242995804368617375654706845366229292012602801749810618414725967603901949198784625341977440992358386515801242403598752530744965951958500123648768590470903665241289753604613839243045838576267349331763636209591614908420567209873762034891165703024219243186604038814895078590431841222555578894267164701830112489493215639233796398468416401534374933912056537148645457504855149765405280719535501300783530207472407251503245360376195113920714400897088125924634586503895324029689665456377078764014673672506620294820976839851506877290405563704084292915796664070411549968522854721000078268911688229988242787348790851901113478766695388796177624245312748458741198285245905022847869724938283712981703924931196734787643612993774222951361257716386196108313658078650141034351493381080215813301202986497196509510411939979675849864208647046323904733305021519869655479934937169443848900403660029553023369364420427911871459666788863057392820650541922533130500523792892714399470619760388127824368993344912627639559852431666230211641336267088096288222560477898088409625443593257725135624499094040423440683896341538991909606308213369305100503725695011298314066608071055921987783384259225247503952646288430711372141655771590733055763151616152934138501809764690172716955295713224"),
        (false, "1p1000", "686.61540624072206937737503082164570364260495703600219568992626193258385411753375539279254455032462341032303678966407554786799129443957577721780141100283855203650616162406556245852979751933745423431770856407981559986247310497713428730417730315604262717408520504459941952305507935136955354019883019846131259036984302269220168625385543326279504724920706975065078868623134695503691919941436025538091575421256701670029885697838162616867695716217132853872813892726554246712596218230955426726256312164415760834639395941578052391471255777557914990707597397433262598369373807713203575259616735319718885662238778157804966074525483902121881221069407131752367505940132044345783325499228282133347046432416889381916298985427200353167301188681200825727078169432095986807944319941096009954812268311400810063045530473115508826925966170888525924216680206006898315485668041873128091786924188171206912317840144833441249399226596363558646068694892022114790752076927273010911837844694248126374511582574150335352920742755820585584395838597990160809786360784856249085606936528359608023533685599432848700128443781122999149909353500375313297110718468394945786861166252478452112740117958664566692429974440289551912153744860081217100649916348450932462484486229742192432173912205583717116192898955455907070035971834765828702826781150192593614295376540688441900061485183441827020685941401856506259831746815138963113518482663393597704083998554247478201049379032916762985191774319549326644032737472399178948146536443342647270906253049471129983802674446817100904491460325924572031184904366252238014984827838289462462474142912838555204010809566260558439498489225199061723711365338973215754864162428463565874220491683295669420914344926874218384855757306596162478840350570396903834288418241678388180359946519034664133407134545454176766984378952308808513907718518996321917574953347618409112229473512971756436644560941017866160841207536765963424539043004565238244600855173119481408286976408991787456242961882622752025125739687375043997586685101148761883665132179653706088935855969616059537783533258615998781807020530693507211786875169145913360247648728814467503267368863039238055406099230849335106997098371979446848362812643988737642199571185222297859934649742196647589920114330939569084294587780260088386991770821803050832735205994235716509710328106617524004476776429973908409644860357055248396207750640624935977551405849294720621618786043310346900900559889646897516103848462203436829094616915764118653372227687379160039156735432372342900368201865792733162235058205274904688924042442947274063888087658826957671077224664680580534869202939"),
        (false, "3p0", "1.0499088949640399599886970705528979045894669437063414529328715833166490504444429578856786668224346744856585059655764803453462982505088160401420535232851596183341340977340407662783137190600537497575529411984255920994948852604584468169394889641563527397176479807016676800294741603119656300120792169904016763283903967057181906918883720328491432169562483452295398616594215659564098129445117098606252480836305783306436057131773377302807102071504086592760852095260458153581728778507034111268965604357378815930091204391451070194772866267028861708846942892050717580419523397505927284876633805180790075411407143582064590949958987662059289561662282581316606808196438247860834617236678612099660555955215939970375805474935408109857051608732260044194806940853713875884490461138995812952229046548206844895074100563851066657624900654240462329377730062333734334206936109711716143030950556026425766945062638685487779744802004092393680183022437537965499048329020025601557461094677123013204239400893520595410011533885804469936022651900730911288731672150656913068191168065314669137461671754668935771279355079857901903989008571854277291809753626952491874580870822376721985757718812945410304768513643813220446998621282947376187419916346438257317908352692146984242695079938285541034003585318763149236675868316739431831070043109025284261290610849780234474090867972294360929175940627273783662345016886098173790151043387461804733428618709034186467710353586733150264679386190633589875065682451597218987855696180532007938877019819213364096805820459299862243427731141744280272167882305695545773615672789620972362667934769857192042890701606015034041296750363293158953356683706418023765872889654158402743138415626797454771740263909653226435245038221740272659261280668380279472678107165056710234388694544613569434719677289940070396186648259673833181682162143538078994171151381809633594404056817624935733243095936915463884012434955566367056150986920056542621690186484456585917375101144861373819112921632853781488838849949712567805997090768753190773304982932460607714552849346625876756631092839698789867135894944457609073005893198114829645034297283769837491705396401620364228351631306269715789996733834981801260056289443208273222322546921182946960803243351852988553916307580008271149610741754920473267656527692352807712349979586476103410334373208555321064224268889167965834513903983060121990814208041777811095259088689265208482723764821762010944426088989410759748544918461284088910421119040994556755472194061533251819628927454456155836574037596935128474001213360049439390045091695374310996319790037994999636480164808194980332138841008"),
        (false, "1p4611686018427387000", "3196577161300663245.7334941746101085879599042688547673645042527530135232915400353570598373918888869241102727906668420937691169553982847730567162063781381217197885532621994019450574909202604852496894468614317509340286109377378852976141497726715592566761585955624912500209505883306728941840980524951027968813427363281803495758418411035078483488499298082088098969817480798716916305010897693882089933073959997075755702726689264630173288207331742418880911597057896187035608889950186916436793994045119965729579312450674154064137140453738615140825690215213656571151836715834238129556756212221348115266963652334091445924139118712023713933555205657579767008475715016291233650443799903946806312235768834933507204381881834885666520879107450990182392438869193980659694516065033042153102981473061718487410846982914152026380083140465607381063213377153150211523772394246796496541714425095665613303235649216638818194107915988996452439720745818143140322650408540633103044403062751406948017811499869351178080905600534764358081535696735957336220341713846915186973766776963639148923995263663260070522071464926820482686151189757883137649643397143099363767277351609431322144438334283197754398874500127312030756769107603706209177806595089983022612418524223112654128486160391161466953096268835969765110886786090839971174284658387336760023631828137340249153314910530367445458433594109940611789817830289500350046758985298051628514881434055942910151483733366809272598569703070964129114057539729642499735418966324213085712914794389130261857876425220454532043079320532297788472125075315111208871435704355710229438853889633485577820942191900913014168896350438953620960990679911418120516422548294697715208139416126041616150930893278893860161708392575583677087311893253278161815909400406056292860682386222094911987343945124528566187338656988357358669096302154365881432911325826103733587532923066428620367109042434318618265362068867299565195423377287344784673965647959982446220433593190225642689560266663790007528973592404579583677651846192542014074465236066370813704191586088292659996162476536078802689720839215828656736694013226336557135924132597150801436730554424696896140695249063025326044764383609086676096925752642692000466349930162903191965760570488705277778625195737004689662772307453747116180973284228750804128265578452700101873399204507922150650453124623837960552201015820887413017163774208385086074864487604179787154597100800939978609023503337617748942064227347398993312504205249749506437904545828970498418836761648447082452780692745938984270009829461761714953418921444663255157751272243637001730249565932676095690818380154"),
        (true, "-1p-2", "-2.1532923641103496491690991500929813755362064853194776958845115077213625846506493789457279095380274391182277200936280737293428580696253913399592325083620716748379232980650523958025453954364074192466006979985698892036250762094152659979731045605023374576258952883008101806719822665311062608263579704179589301936198472864880108436929669116493604260992237104683225561490419970131605031462920271177808845699549823423385873774502105702150972312666684661338947438359143176760728553540428917464014931605086730709530329900004008118101632275292268633789358710301672604257809202925666702627716115406584989004622643956569789407947347887325570979345015112779899732684833246107348888354660073907675075119346977736390583132821978505731926772609516699895464944587549152409887233434834844342222838291242920885304564631169924905007611894386054972042881377934676940379127006303937932590430762642585379305409567941436234141585505383455064969974332763866334861817860562480810597465527691186793749679520619336791943096274938550662439620498846748236380968581913080627943474977135482644799892869996678571522230293217657639465661717170609590436795634799120101574591201561648099056237719639629684842545049556141421036740876243798254803567386048008135575617895758202720451145519676926852640805775733679692455413442018518154331125167626339093620890896141250297907305696351133537819181961451576378047384693714259128687771735924757332945723669449648270650414492727501827534705482762669835149722375952478483160704139941412112405073496336710844016298451110347881253459895953753486024344535478054943735840525894955780621342654932246704262845144059769985064327986794322688822883402149992819018304319713772730935097120550843206202380000561558325212423798686439017817352184753765622487569284476757599171352719306402321748529134793908379093789160611523470159063096207405144003141126343710087058856706766218698499819358179345347914362367459342490889558978097539657108360064707534029780277793960936986784016986020873977954670779003968688265606485056482465866417449257920118358030212635168942701571355957124930907410633852544780449514925374903598592548995137936057554629007528689122301523650909915073992931488638967431136719302585516693632154872409664833174731398508235178876496324648497305422371545124683429612162188737115281768836060131575019149452398388174891506408752040495654542930975955013188231301346488904988633674477894315007157957052717391528097952871981792265376382440116150658668193134767229494095680049064375833012034861947914743992335436325296005904742224815803680596184008418370213719730501245981961982788637301420954507544562"),
        (true, "-bc5ab1b16779be35p-65", "-1.0000000002603336262259903935616064211549208630528991669461407641786577259015288706227095257850309301409679728279162186299601348705418445703331665654906172428719910142949755647869498244139256471202621298523681452169844054476408555811179129782338947400153994304113642154301068266967720366084673976413076843141711955472255958488595217079938841875327156457809536187339330257079947901014488857495346752222488407327393094912927860117938856298098876407699874038704877589942070232088403600036058927076677353377443866611736792224895655202882423720519196292524151928726073323892769884202802724823189761204127988405282053396731408360501553925934508946359407120390360715080178011506854466832263080923016640411976121932697733765235335192188837095726056314203469540704335598908129326287569884716484502930019201030560707466355677930244615136232162420618894208727724293268567787751546184417913629976877218642507429147493787448043277538094354300061360781467345642886066361816928472300431499784534243483574982498890340011885765988539842111250284371993874846044070184834428427708270077459973994098716510416437244301641114893839707636892679131461689852882892234470019268732602034599092302310473162389147308856836189542247678282542158773619738001371330829302806700615192253957701760852062091254258792156380981319823314762715303042286470376015969377061105040247454003250331346917741863257083979258043913191187003476327331764874537801356511578570779901918501358563652482682838852148795313410404214239765901221814169869157431317956942686854905734396783255726985777343078926487318056677740694116460885084908872998060293847197416185091166436775429895766435995144433866026195154121637067272232504251766161374592987750081856502176034404308204060180136118698389031383473801585504443212275797737215879644698150799309588346122642234241708490630351344586726529810760323632598248854641957532221071696441923989425660034847299933884227160510669003626938291654874959555652247227288507569815992067516293631028758808051438065012397461743036770698104538501764887630798553844385373440382032841505157501549444308409839291365261044291642566124194359939796185616016615852018716917822807155125775474052909909207084448509214108942548423912812341541953130314541914714246692273904502842589159933274975769771093376224778482710042179889871028822967449612129754529490277479115804538574556471989784435322721265962750060708572167225297421047787142869490103988013113160002724112787423165684066920589998826356611756729395982823520714041898746861111995831960184907606935929946843618007270995824580751559268464253548490053450074009221570143685338629682576"),
    ];
    for &(m1, xs, want) in rows {
        let x = parse_real(xs).expect("x");
        for p in [64usize, 128, 256, 1024, 2048] {
            for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
                let v =
                    if m1 { x.lambert_wm1(p, rm, &mut cc) } else { x.lambert_w0(p, rm, &mut cc) };
                let w = cr_ref(want, p, rm, &mut cc);
                assert_eq!(
                    v.compare(&w),
                    Some(0),
                    "W{}({xs}) p={p} {rm:?}: {v} {:?} want {w}",
                    if m1 { "-1" } else { "0" },
                    v.err()
                );
            }
        }
    }
    for (m1, xs) in [(false, "-3p-3"), (true, "-3p-3"), (true, "1p0"), (true, "0")] {
        let x = parse_real(xs).expect("x");
        let v = if m1 {
            x.lambert_wm1(64, RoundingMode::ToEven, &mut cc)
        } else {
            x.lambert_w0(64, RoundingMode::ToEven, &mut cc)
        };
        assert!(v.is_nan(), "W({xs}) branch -1={m1} outside the domain: {v}");
    }
    let z = parse_real("-0").expect("-0");
    let v = z.lambert_w0(64, RoundingMode::ToEven, &mut cc);
    assert!(
        v.is_zero() && v.is_negative() && !v.inexact(),
        "W0(-0) = -0 exact: {v}"
    );
}

/// ZF-72 (Harper's abort rows): c.si on the imaginary axis recursed forever when the real part was -0 (the half-plane
/// test read the sign of -1·(+0) = -0 and flipped back), Barnes G at large negative non-integers recursed G(z) =
/// G(z+1)/Γ(z) once per unit and overflowed the stack, and riemann_zeta at large negative s aborted. Si(±0 + iy) =
/// ±0 + i Shi(y) (FLINT acb.si), G(w) for w < -64 by the reflection formula (FLINT acb.barnes_g), zeta (FLINT arb.zeta).
#[test]
fn zf_72_no_abort_si_barnes_g_zeta() {
    let mut cc = Consts::new().expect("Consts");
    let si_rows: &[(&str, &str)] = &[
        ("1p0", "1.0572508753757285145718423548958779590240539375698078996121035224780412665464728491355554805806381034632846649744093816663822638423500803490224455852078030412818429029454807105756432660001479985587471543013256491402636053452675965441454649907632465787981404398184026973765652062134265825257518069738990515302189435177273079877256476149807313544558162046966522547540076731870100328835087700871808952748382066389690047510909364878207485010833135529437737596160298707748958963951177862455985403517859377837704295622769568368576430761695547573324862563890098181021033968446096528606079122499990255312259741603128443197927783961705510635464467412762853542760810792243862531709425526507133554949508835199566709730524693845630536825899008917743295895186125117869521851508942946256043524057242530009199186343609840541279074706886606160929624564841853660183029106316634899189498761986265676870004977656419558158356459743705775349533738840163325424867683241875516323563516555521971733228834019779236916406732837325098169099612368069974900454388036934717777696099211538671624921093999123579786196817162988621743336718805259182363959962324719953368033931930831675760752862911145855894335193529249523229272780537599548049801325660149655217665701976879542811251836245512308203741389977971974300205520616500352227799468228430668019590335933248463605500417563236433307559772749586756743243092977145576872278118198961287184946887539212950224525742033310798500134376967106953572370253122171882063633337933120370613024375233499462347848811788142495251105643697857276589678585734819807635972337612910516381899146620248046624612652624605777007849219813116311279860112170904988781833983090731974746452668612373946202689150940328061387950036062113637181361301665640711007967790471313785856308022459554468699335011980190153003856121934231370257034271719676091066974336254433591715648230964429478762149200313032387492068310141136089430691359939468909544226571959008415593574501847931612695548126353443365269301192571013302512440773175045672363398033467790116374367073792479186373931662315132721317588021501600474498860817415150338871327922097154254619473606413546129155301921468527822394095851747609048648402879568229125806436521883636000328057622000119362744188591007160311426147251265814740577850901786431095370855068784610057869855123350106772655544402090836942013842312919573552778717621191857385184127124640914028778702710904473365421199987949630052683825527346132899576236644017469276791653706643932429581025810792279537907902350293344054121517144805822469686481060357895406679685611802859424501036302238547716665907942317022217576558126508594226388530280688833105584022636697574634011684323431366706674160289607624205864127560686306933312538389953721905631405222749722"),
        ("1p1", "2.5015674333549756414733724827275423989162728736915325854272090641412004210062796910811829664301509024544269652992904949357110015675280567835475738071172859372186785181447800318898173253127528325962026738048195446829391031650602318562259372649462713095022823832599984258340634443797595944218417136559074551267437833065784412220907458790912701502833197812008106596071572194085927015340819735807381183879676585367822915811038036168509700750372357029838812172768465511388371742704323318755459763232129830448829865111604942528786316727929354770025348126530212230759613342284806199107441620819947428764752454803725739422634187273586679519110827882307940391378783117116852574449254530624934281655350587208011418606275339772779911676569053133774684739491620839378669985028352504220374645703991272755788090667764295170031133698331823090160490092778587733826604689797211196014712795634284197215111868347468855922680756198599668580130403797205852382837238609688765460236670675245955338811412566423820938583203880304978503985670875508691537481406082880342185282966665165501069779777078900647920271456178728635689865279128358439123009410188810727837841048656574314143368894308192009575288001873525591806766952251252984908098934530428455367752826852652634011736049029696537904709230750631081984355180189917131943300364870180770396129235714015505146477071445202922623183172691279765157663317089091124113787643649534826715004612796686916154667837869440016860291068849626627108021910516614224692521251021623333316643649381791826297696738752779655036041956582718661162475492298914515535176875845728076424154449922506187810804359920967693176351562635295676109809086006988062742044420720326638277278565980987927484383827986309882025379505571730940422072238651054680467808612618839120409170298618645408346184189529182359095361264967356535923592776018658130041433875055177373411817444221658935013284963020885955499708548307860272278585683394903905429426631433079826202589511327823448296352477300243067701492142528830138865417153317745590020391783622313433467146256153964944129276397722732992564059935463897627561137726689210604875882573750396621450191396113632450239343108144933219914452544099290408954347959085839026507269449248877743493475689876606648703903234237735150232809238844924030855664792424414897318090597865078895047370994012390478046097040766101420537675452303839217868788235473462442957276439626577815709038285456850488664043727338536491777443678963010362410537486521465834809574308933773306537541429646194206907367845001085748982611618860651768670260933846079468808776528069447435702977480605066096191300814148993223717860043822554940651448005009475751594535099859208845185098970821530220716870514774020881959044887180966091354571769390998344655643595162042"),
        ("3p-1", "1.7006525157682152449296603280232048829591684957968863429067405152550412996628855737697183502588850306372880029492154071472661726677452490767720774595289890492826611128697401954575562743639480759298202263030312965953710389863631789623957614610587203873989332812748978972596525180366101384411190790698607555402725228096231408985562410164605921250244104159186922473459516053779415766327009620967023388128060322525551101265799378705306362642086295329133272759465306058920147727384484325487591648413670999333033572351748324050857007339711176551904384567545098115922235964867491513090189024750733480053356645992145111517371472311666488737106390755976854144438352117819587205182430952418210295270190532539377108122359892039544691519054519235961314157192183274004469155882221607966881206018369513181727303373631797176402865400093247966754242858996327521649881390901229385755607403606456374226966472954854679931250648241998651896421682999086482471177825685979837925137855153494818386664648557554761865855258139515286373052258001595218881057623957832868945448599913497730973786222143437556036908232719101350911072279404956912405404792282459054879949607241861047931310724434407225864024820190974523725527230742788014324696652205721567223121710746347887785697253330207596868882123882422335174980416031473390427447293320226333831184218461135611049759126721165620321922413331971909632230260437473127004245646179948587482585612814001991228366869944772055224992767256144100253584335640241313981105653116008532326673422984724513145412114157637461534512282117370384893176320336292700371091136169142603070209799128240508154893513169436099529191614664398687651934656983989493699714826486581750449147184180104769871751641025542650462591889529855014793950067765617987089310033680722128056476604208535940156167415599523804523321296090445699187220065220381835165033528175354136984492130815108234724445888853681374080618172053813380053022534396974002273351542968961831946982510493649929549036929193319115880347893879480590206053728650175092397733903218108167216035572996200459412273769647646822859455643858462757905612216044972155814503908982195027217818772003086540400784329694499813062249730113089611545294403745003321997412254924988893799326750504850673019253187578167793985488108433506990363262420204336661597171742063609503963879709718804859218535812199942802084630033689842770134595400291464374452530433060996129541867791689273348874333568893204082377350706833115106325341112824315885993023648016282066159910600617210848310325380557558698445736659624991031519863865915708033834921853435676462927732165896276781027175465165524979617778010280811380756590394120777555389906970239110131498689399162312678624766524773886522324394662985682778913753550769769712330666449604646"),
        ("1p-1", "0.50699674981966719583365987598894380025412622233449778020555550289774639787846826249826070552238121325086843289232479895917707312620664740818204302778963225207461717607704366444571366637814264872175254197065703702661914617952410447054144789327349551491239014308034719437236281854348918281035494245898279339095004044057718775352448278739508341346815291296684703793169064080635224108298951236419581422260424163949111845453875616324217770250126186419171137191504900867488534772880404218232799144703052989692552727305347280964387625974596084083192684273873942100458532170044037184334996909599901003419828528732137135477268835658470290614909346929947462303656690251953044529178200366121330178007595647792357628538989875475391473105611803156088526774651727508634966032834694147994620071579177530588126310127494183307708140268242939130108547482823313469394685467451345818002608748311459883207781933773667925828730160786955325731258522120380793048874789583175645815410923013006275906520133904966287420237488191829681804104900113725988470068811745412467883220108907873612904798743693649604907145820336863004626188272045739693784740034482568242857236415828603242246194285958870717942301799918072466525953186051724791630074628980829670792314720212661702383091595106435830046619142630862630110609896674480253548696287833215869322878635034748443864631274218315379098068714275123034292922775295946106863359581636999826803287769449243158833309932381889380345968488905139480143768257838922200558426869191445363304716946401245263361642967368945108045170743090928456569116188634015048595694658685224210721692467995507076337638951542472448732430101909383213009827512122368865338630023118572214716523292549118808502706069206598642932644350177743548205215826340618964858520692334753624748148993171256509688311623676375598149603791074056241081085442452951619216440932354529075105167109660317665857765373277895038589370977763796521663667713344426286773751041770242509560685731148638018208207787680409120257936725306089815938114045957882998304288341686004863831800918643776312474063940844160909452504237580191275153058568754744090522889868628760554168064062323646069135916291508034336296368781355318368561885675478192073638699467430412930593759504091611060827433793146610201909660925727632957360825406416779940110789792142532712608407433657041731759710152355310496005024921858751388925687375000366641851957318316180791704269658659322529139670853231972856112331934265050484865987156421990199208215911840196414477872867951774771043024784375161243168386944009696909975027628183012252024395686609374353079070323409878848555421369935489033389072864220072960299246348219160384177220394940693641472967169369985885437065622643058193482208686895608272021571545641978826221350708495450"),
        ("-1p0", "-1.0572508753757285145718423548958779590240539375698078996121035224780412665464728491355554805806381034632846649744093816663822638423500803490224455852078030412818429029454807105756432660001479985587471543013256491402636053452675965441454649907632465787981404398184026973765652062134265825257518069738990515302189435177273079877256476149807313544558162046966522547540076731870100328835087700871808952748382066389690047510909364878207485010833135529437737596160298707748958963951177862455985403517859377837704295622769568368576430761695547573324862563890098181021033968446096528606079122499990255312259741603128443197927783961705510635464467412762853542760810792243862531709425526507133554949508835199566709730524693845630536825899008917743295895186125117869521851508942946256043524057242530009199186343609840541279074706886606160929624564841853660183029106316634899189498761986265676870004977656419558158356459743705775349533738840163325424867683241875516323563516555521971733228834019779236916406732837325098169099612368069974900454388036934717777696099211538671624921093999123579786196817162988621743336718805259182363959962324719953368033931930831675760752862911145855894335193529249523229272780537599548049801325660149655217665701976879542811251836245512308203741389977971974300205520616500352227799468228430668019590335933248463605500417563236433307559772749586756743243092977145576872278118198961287184946887539212950224525742033310798500134376967106953572370253122171882063633337933120370613024375233499462347848811788142495251105643697857276589678585734819807635972337612910516381899146620248046624612652624605777007849219813116311279860112170904988781833983090731974746452668612373946202689150940328061387950036062113637181361301665640711007967790471313785856308022459554468699335011980190153003856121934231370257034271719676091066974336254433591715648230964429478762149200313032387492068310141136089430691359939468909544226571959008415593574501847931612695548126353443365269301192571013302512440773175045672363398033467790116374367073792479186373931662315132721317588021501600474498860817415150338871327922097154254619473606413546129155301921468527822394095851747609048648402879568229125806436521883636000328057622000119362744188591007160311426147251265814740577850901786431095370855068784610057869855123350106772655544402090836942013842312919573552778717621191857385184127124640914028778702710904473365421199987949630052683825527346132899576236644017469276791653706643932429581025810792279537907902350293344054121517144805822469686481060357895406679685611802859424501036302238547716665907942317022217576558126508594226388530280688833105584022636697574634011684323431366706674160289607624205864127560686306933312538389953721905631405222749722"),
        ("-1p1", "-2.5015674333549756414733724827275423989162728736915325854272090641412004210062796910811829664301509024544269652992904949357110015675280567835475738071172859372186785181447800318898173253127528325962026738048195446829391031650602318562259372649462713095022823832599984258340634443797595944218417136559074551267437833065784412220907458790912701502833197812008106596071572194085927015340819735807381183879676585367822915811038036168509700750372357029838812172768465511388371742704323318755459763232129830448829865111604942528786316727929354770025348126530212230759613342284806199107441620819947428764752454803725739422634187273586679519110827882307940391378783117116852574449254530624934281655350587208011418606275339772779911676569053133774684739491620839378669985028352504220374645703991272755788090667764295170031133698331823090160490092778587733826604689797211196014712795634284197215111868347468855922680756198599668580130403797205852382837238609688765460236670675245955338811412566423820938583203880304978503985670875508691537481406082880342185282966665165501069779777078900647920271456178728635689865279128358439123009410188810727837841048656574314143368894308192009575288001873525591806766952251252984908098934530428455367752826852652634011736049029696537904709230750631081984355180189917131943300364870180770396129235714015505146477071445202922623183172691279765157663317089091124113787643649534826715004612796686916154667837869440016860291068849626627108021910516614224692521251021623333316643649381791826297696738752779655036041956582718661162475492298914515535176875845728076424154449922506187810804359920967693176351562635295676109809086006988062742044420720326638277278565980987927484383827986309882025379505571730940422072238651054680467808612618839120409170298618645408346184189529182359095361264967356535923592776018658130041433875055177373411817444221658935013284963020885955499708548307860272278585683394903905429426631433079826202589511327823448296352477300243067701492142528830138865417153317745590020391783622313433467146256153964944129276397722732992564059935463897627561137726689210604875882573750396621450191396113632450239343108144933219914452544099290408954347959085839026507269449248877743493475689876606648703903234237735150232809238844924030855664792424414897318090597865078895047370994012390478046097040766101420537675452303839217868788235473462442957276439626577815709038285456850488664043727338536491777443678963010362410537486521465834809574308933773306537541429646194206907367845001085748982611618860651768670260933846079468808776528069447435702977480605066096191300814148993223717860043822554940651448005009475751594535099859208845185098970821530220716870514774020881959044887180966091354571769390998344655643595162042"),
        ("-3p-1", "-1.7006525157682152449296603280232048829591684957968863429067405152550412996628855737697183502588850306372880029492154071472661726677452490767720774595289890492826611128697401954575562743639480759298202263030312965953710389863631789623957614610587203873989332812748978972596525180366101384411190790698607555402725228096231408985562410164605921250244104159186922473459516053779415766327009620967023388128060322525551101265799378705306362642086295329133272759465306058920147727384484325487591648413670999333033572351748324050857007339711176551904384567545098115922235964867491513090189024750733480053356645992145111517371472311666488737106390755976854144438352117819587205182430952418210295270190532539377108122359892039544691519054519235961314157192183274004469155882221607966881206018369513181727303373631797176402865400093247966754242858996327521649881390901229385755607403606456374226966472954854679931250648241998651896421682999086482471177825685979837925137855153494818386664648557554761865855258139515286373052258001595218881057623957832868945448599913497730973786222143437556036908232719101350911072279404956912405404792282459054879949607241861047931310724434407225864024820190974523725527230742788014324696652205721567223121710746347887785697253330207596868882123882422335174980416031473390427447293320226333831184218461135611049759126721165620321922413331971909632230260437473127004245646179948587482585612814001991228366869944772055224992767256144100253584335640241313981105653116008532326673422984724513145412114157637461534512282117370384893176320336292700371091136169142603070209799128240508154893513169436099529191614664398687651934656983989493699714826486581750449147184180104769871751641025542650462591889529855014793950067765617987089310033680722128056476604208535940156167415599523804523321296090445699187220065220381835165033528175354136984492130815108234724445888853681374080618172053813380053022534396974002273351542968961831946982510493649929549036929193319115880347893879480590206053728650175092397733903218108167216035572996200459412273769647646822859455643858462757905612216044972155814503908982195027217818772003086540400784329694499813062249730113089611545294403745003321997412254924988893799326750504850673019253187578167793985488108433506990363262420204336661597171742063609503963879709718804859218535812199942802084630033689842770134595400291464374452530433060996129541867791689273348874333568893204082377350706833115106325341112824315885993023648016282066159910600617210848310325380557558698445736659624991031519863865915708033834921853435676462927732165896276781027175465165524979617778010280811380756590394120777555389906970239110131498689399162312678624766524773886522324394662985682778913753550769769712330666449604646"),
        ("-1p-1", "-0.50699674981966719583365987598894380025412622233449778020555550289774639787846826249826070552238121325086843289232479895917707312620664740818204302778963225207461717607704366444571366637814264872175254197065703702661914617952410447054144789327349551491239014308034719437236281854348918281035494245898279339095004044057718775352448278739508341346815291296684703793169064080635224108298951236419581422260424163949111845453875616324217770250126186419171137191504900867488534772880404218232799144703052989692552727305347280964387625974596084083192684273873942100458532170044037184334996909599901003419828528732137135477268835658470290614909346929947462303656690251953044529178200366121330178007595647792357628538989875475391473105611803156088526774651727508634966032834694147994620071579177530588126310127494183307708140268242939130108547482823313469394685467451345818002608748311459883207781933773667925828730160786955325731258522120380793048874789583175645815410923013006275906520133904966287420237488191829681804104900113725988470068811745412467883220108907873612904798743693649604907145820336863004626188272045739693784740034482568242857236415828603242246194285958870717942301799918072466525953186051724791630074628980829670792314720212661702383091595106435830046619142630862630110609896674480253548696287833215869322878635034748443864631274218315379098068714275123034292922775295946106863359581636999826803287769449243158833309932381889380345968488905139480143768257838922200558426869191445363304716946401245263361642967368945108045170743090928456569116188634015048595694658685224210721692467995507076337638951542472448732430101909383213009827512122368865338630023118572214716523292549118808502706069206598642932644350177743548205215826340618964858520692334753624748148993171256509688311623676375598149603791074056241081085442452951619216440932354529075105167109660317665857765373277895038589370977763796521663667713344426286773751041770242509560685731148638018208207787680409120257936725306089815938114045957882998304288341686004863831800918643776312474063940844160909452504237580191275153058568754744090522889868628760554168064062323646069135916291508034336296368781355318368561885675478192073638699467430412930593759504091611060827433793146610201909660925727632957360825406416779940110789792142532712608407433657041731759710152355310496005024921858751388925687375000366641851957318316180791704269658659322529139670853231972856112331934265050484865987156421990199208215911840196414477872867951774771043024784375161243168386944009696909975027628183012252024395686609374353079070323409878848555421369935489033389072864220072960299246348219160384177220394940693641472967169369985885437065622643058193482208686895608272021571545641978826221350708495450"),
    ];
    for &(ys, want) in si_rows {
        for re_neg in [false, true] {
            let z = ExactComplex::new(signed_zero(re_neg), parse_real(ys).expect("y"));
            for p in [64usize, 256, 1024, 4096, 8192] {
                for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
                    // Harper's rows go to 8192 bits; one row and two modes there keep the test quick.
                    if p > 1024
                        && (ys != "1p1" || !matches!(rm, RoundingMode::Up | RoundingMode::Down))
                    {
                        continue;
                    }
                    let v = z.si(p, rm, &mut cc);
                    let w = cr_ref(want, p, rm, &mut cc);
                    assert!(
                        v.re().is_zero() && v.re().is_negative() == re_neg,
                        "Si({}0+i{ys}) p={p} {rm:?}: real part {}",
                        if re_neg { "-" } else { "+" },
                        fmt_c(&v)
                    );
                    assert_eq!(
                        v.im().compare(&w),
                        Some(0),
                        "Si(i{ys}) p={p} {rm:?}: {} want {w}",
                        fmt_c(&v)
                    );
                }
            }
        }
    }
    let bg_rows: &[(&str, &str)] = &[
        ("-5a3p-4", "6.066050422392870809991731319895049612461408383231845543231871945905528436280873872067089762172207998068532964429007503707395193217157547429530890924100246247618462235258501415628496456898370690929885755389008083867867463695542999482123498458206149687972815983783183936591107151849058554759047847119083617563867467885448927408634424377413298851473535884702605772894334017711263366960660965741484708972259735075522481245163764002073564076908510518805235908883041032718974910538868153360244407848068955866184279035212443919959539856368842339861826878586098229657880973740553864385596199215449983213403442512872268372399263504432798222754072563687177187366404876916903232276319051840563565950261047569871707598985841197941314514246447052594418527866674834255557084978158298944969260791921076108253277060438030421143173662660700895018746861134015442221316138278678804653007582020182776570439643706789164502682255050724826441465395499134700909421158469524171132242879401534951711907064195094485696951369287e+5407"),
        ("-1d3b5p-9", "-2.382256879930204464171291710756232435928457564178274241073692534885574269944101866029313921548324649619786496455687423919390107330407418119710076260539709965891025199130332504037409632793440960989249193468977584114956412840922293250297229950450689174635600782414745753103382276544486038141264815394937140449482931893705261348557523240582427025850271538683817156537093745827093650797049143167940583726657813772051544612486727470747307515559134762487706266467785707993048130968707691676836267459683159123781338882620552442482472860574878297069133320466367460780463827860920677578524867872723235302221205513076101148448023855330490226902559012548676592013788920429256965494277342442811690772039159439683923583935314637493815106461942462467171537296778941344795513645985195890000209171467614290944712074484737260321008904972834973256001709798800925130515443356096248321932107587544959997319085109883121023632427879926468399490603766411534262216373854958794419330982449616551128883904345863107746359174809e+47310"),
        ("-829p-5", "-6.853478159684898090414221652978812344467726613792040725712277571421705186075003424240942746000553282842742267201815984278333682171714584256114544547329912748138856549518429786080795102539444766011119266422291103881153312939800290372367643128182801146552061704621279791427558269058382669393642891986020007523490878209732388366329859643344740924905679481726596171089348115071942194836921297297422530564141933816877081881140003510102166817784896243183442031982913487264352514796219795611052640585953792760361461716770878105180226729312185913341795283191686536054108404705139547007999775526394805316419169612564728468495929576995874736109717245093229101953861859792435191377890123059733234056298121610262776524435053867582694180312392817239477739297303597510706883838724684128712935851727196849833113293267512249843027229462124310466126662575750603918318992160725744951764898987179206110870401907967499280175601852241152606696914456329606413888961292377620533060031086991789403492875838222419783167692381e+2555"),
        ("-fa3p-2", "-2.099945786305473375328108502229946573984939299451185358807628505263856531114475646677272439569152695784986583206949427826749845237046476245543095069289111399910708776557266819208289063294634848519443129578743695889445056305169523675954440092698606530588745233292982393577349877921456662061500577533463968215294728498937371491002897471571911017587517463365529201096827449757714722262766505604301127319307711685427792523832679801266189716046141675911923084028756744097398563098665806910922028783046674478223457345328099715325284982130798348162626760227126462813193194280833006291063635982098767103714111621147477347487610250671450421870951327202291042117745013596959793019150352352845478833260922424340637809092851704930105805981512504396528764335361160969303926300593826414995723122978677039430620076464568049004817152647292628171960701746595745332976883394394971767511247564813514580332344090107285290001650851560892962010092369785844837754089169816969219039816272671797880466807164814521196468174075e+1178524"),
    ];
    for &(ws, want) in bg_rows {
        let x = parse_real(ws).expect("w");
        for p in [64usize, 256] {
            for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
                let v = x.barnes_g(p, rm, &mut cc);
                let w = cr_ref(want, p, rm, &mut cc);
                assert_eq!(
                    v.compare(&w),
                    Some(0),
                    "G({ws}) p={p} {rm:?}: {v} {:?} want {w}",
                    v.err()
                );
            }
        }
    }
    // G(w) for w = -0x80f7671c63e7f577p-33 ≈ -2.16e9: ln G ≈ w² ln|w| / 2 ≈ 5e19 nats, past EXPONENT_MAX; the sign is
    // + ((1 - w) mod 4 ≈ 3.195, outside [1, 3), checked against FLINT's G on 400 points in (-100, -64)).
    let x = parse_real("-80f7671c63e7f577p-33").expect("w");
    for rm in MODES {
        let v = x.barnes_g(64, rm, &mut cc);
        let max = matches!(rm, RoundingMode::Down | RoundingMode::ToZero);
        assert!(
            !v.is_negative()
                && if max { v.compare(&ExactNum::max_value(64)) == Some(0) } else { v.is_inf() },
            "G(-2.16e9) {rm:?}: {v}"
        );
    }
    for ws in ["-1p1000", "-maxv:64"] {
        let v = parse_real(ws)
            .expect("w")
            .barnes_g(64, RoundingMode::ToEven, &mut cc);
        assert!(v.is_zero(), "G({ws}) at a negative integer: {v}");
    }
    let z_rows: &[(&str, &str)] = &[
        ("-9652ccec4aaa65be78d5f7a848ef30a3p-79", "-3.0715589505054535490687899860240189159939890354139872094847903984457530336047913284066730302960892893413515728089153069619660628026435225586630770770874053646900301773890188225589647571434290743236421598140615943644706218374494853349996621162650979847157817020270067502588285696745182836513824289286397017725771631147393185966662887637654457985431446312922468541650510171195377068189637075424965683371940147661403739361480619535083557518139818192509333120693250924430123493487424191692919368147092775415295351550928559443037535293955597835529039737031069065701326716011140235937635957011959804202778962038663327105005005029426056939621003072489001412719178739761741968214984796713760203874753902727471298271907470563138915001846753939873858861578037333560578849956923840581705553417237403901931349926831023849072872111404513614962478535680331670695496603344887641684964641165427759499487878107557699438649621482705491703059482044207753721104050454186941738212411036343369000717704571022000044417399606783019947921410330015088329435257445012755070431878924876197745364511560557741382102719015078935968055404723008626532333502909200363222162511071608471820394677963934395620127957129985405378200964038949820545873620817298741670670305077813824842971899196052452340730762209923223776775726358211312852858356358690199972542490097570648012845194545283719177885828661621794485553683635653065192451950473777323788852547677459331533834357831698648033976016257717943408417336326011179285259375013080871828302128604812034669686591086794793143672328939091545000514589543198704170370012827253469303832869704768453372767626488377154058268131683832419792696183620226416807749699653338845029629011658583570414913543053109107268571891755228368273654405512733643419633502768216911575233623729693539805284844516329141095844363833345349993364635723630525892613802123550962286268973371470432404777859989410658342697761600409382881721311950182206627282888987698237703101442660630866045318598454990361343331092868820248753144414289323281042679593403856785781303607212422516213131393668037558484234594222377135302810140093798155865202486146990504914309411962459529023658940661801862372741053225942032111773040449513499671851953181651568853094027370763891203666954008447993216361576029884386028570814934702531587458445783762492405072874591142202716797443891502186261586476630866514990839695180114458378039795803121462255436636743307212375375499586323078843699878592772267021136269281180119188102961218619340674436371222835781630398627813459730379773356352749536762452862529823310668420487357036489522461329250884814814281403191758292581739434391599231001357494232569702850728529875566078295631858893877812947531152073086601698433175366112189994866034466283081699492603268930e+4392141682673324"),
        ("-9c9e5ad5c72d06b5p-9", "1.5087127509705042221809258572617172437929691294736689798520239811548305800666826901613842044037830497268979471912264946733491219522696363973632290346802278499094506622382347731748489643057247992869149197696303929060220554331034703786287171007715686058612370012335851437330484794104016976184189268252519089787067079154788828109763462033812907845525001477845004332166018795822126517076885205392875172205273989920135128266875290779554162913737418097419511621935100445683086231398745160177819116143145579150487110235812752529322087435061724131936113735488783234110445242847722627686883783400695515458661307164587396860366227010769039925623420096341879242332749693977015526964013407594591725025416902220701191192020767196737262080776614713198056035527370525989113226711071793579998803615016026641018669450218897143073676110866338498886957782782019257763147900054465185598419192880690315432865229555588720965198401790586741388323554786742799691089638203857299278593535325942010388915282404966280987808680465244692829212760371411653286795711649630542232901319120253950774387607683079608913720238726365707081944431153917261412789187683978186921084150122464458219725376470893474696741155201643332895643532343640425766114564877281210057771774237408535206242440643133025368502400227741783450932660715959731372435643636676605063235016451659973349053008302406643864183416035544951183182096449424136946932662698103809376041914148396156510927457209769616012552440131354433899319734141905332247554178379000172400326500181882032266312494778835306792070951096995659557286674308909008574825234554767154913098364311646671433532412054028599842312682528533906520945821836825410166199814150085582393064153092555853218648573622739908069805139892210612948128911689258544216333819379457053309349586865201827211828194031855769684154702966216318121521248624370236014726612095490935145829884635725724976250920327852394354922982482053094535980764779569585030140618085799214656016798843942028511637572931323122161511832070650682592292225538238536365152242259320285528233404500636879233944329029131068827809936048767217757754238669582345280636858363852480521821548389802119443392262759241983145459422726792149082188276431218618555075741344836189664197775525041885972720661299638482334687773525985939724595022190673533440682305011498044578950895388107027723860500896723889913909520324188185240728441924977872877525351944942687431287249719460164738117273292991998724993757535915290978576750719334596012448412236912613057672025981867221365747071928781159014101342894246148312659721095397877450058721553407586112513693373872968614222963278774483932638732130269434252850460188791498445128267569689418321174483489280096839449900464577524816605780350035106590604694508421669248685997077762e+333073365955539963"),
        ("-aca0228d212f0c8d58512243adad3583p-93", "-1.5270011905073977573940457189568979692463904426520284971544704420935611342104774496407909913638441247704668747595425152132807134005317198876943511697730830834245111046306982483945068306085666572535305278780811393231438539056241330533190136948672668631992234476784837450326900344972265825863764151147940737428111230485131749046694946950366742443564026852958919484675692411509705085891687818275402024039726145011604137307925210833047872552726234207693366129856745411773861043692162303679299689725244056173954636877639841292201201905079503514359414334338981727599982094870945609516968510043615994327710164869950470834197657734597281911602407215624055797258437920371762246671516834484125170115077809052796154183068639926478578786628850838287306669805529059142845176627904583399167412121211291905009531923901039649410882048507569227334706379757652799834227252816485182317505303778932524728033890166015953229306395750341676016901599631642553672668901943501846197942698444119745521412754049647516029740794975341530826801471373644242498633037317419277452287906741336511125811439312006524307623316869938912205818110125658806099899154080251521016619968609547899126980975701532863907515159311450718604066261764917175753127361226807233720993763618847141452095112309930680691462246452198101975927389148414223468367712233987435149189044756123962925664675424660495538603588249000175713612144331742873136290241316630661194122650721271534451785340354696064086749203114223912409263935597875609241363512077014950075141840986076665370072360869116556375677005600782896371988384815223085550028434022747432998742063694231246130883757371138935352909359083080930757245605762491103719786138567871893739945062838411816956546301565802239481499492152824856254108570244270229625712328286484885808104720859304568444014532024726757201964414004613044848591559788411069173876076815548336697479670077755548983140998797942151491876968841340511071444412650530426006015703026407753219483421415449679085466731428339556368782562116461266065044423571004943734840484512961629223497261515523092920517177604712099701716331677501719292516093565139798337403340655583644962424951235601753458671379496245785081802622358697135126818685758210762164340854807057103632916467784686533541357598889515527442716891216227295024812460063787815100565300336739513831647157333595694695083148875051435783027606346791347984993268854164631619428797973575744737592879491547044951529090009003447953827531668610784808921224710851525180215834898092157447456392803933714997552793658236197318676851566730617181323090604869110824482105066569777795358996304245225763565749238141486394997471613285990087482524663177782974241263942691551811216706727114977174345707805426543652770255446107612297302133953514664663209931871481e+211593223414"),
        ("-c2c116580addcec354c6d3e886527ed36127123fd8be6654c16e144c82893bafp-227", "2.4929240134458664463332011691217802745147866710240572592432068203515433856317237794023612276772584539534501883550923145178870064172202428494679481039379629727735832978805157631015295773401548132573626432274815065509151779286553312862666096504667292270013863779049874210036265700478196798472641102075014367484285980735037934491668204800744223987246458878836476410309849027226101725490453825758318961772117037868464589064491373442134486111797531676292253127418844852655930870885150458565331629494624409367766818478638287981589165906141988146653132692102255662020262001814839725288513924780703196232579147351944577279278133573681762334727315286495498458677372347001673050900787658127710767711591382032137202336668865046382959669613123194914471476104539599135716380204322095804774988691291635551017347474549502565834742703199598503754611190948819303362021142974081918999576061791903371698408449683830931292295911146558182162680338712713177171005092608199205854148563957767105744126098578619096737221050600311210596660997440012971947963457809257166577448698099947264149993045097082759741577408481311208294384894008050609975858371528750552272318609113082999292127659386725542790510257460421419238277535814024196555662612306211353900874001924608014174522191954777002329997317129520979064115686241137082509790737803397829800961411262835987272799574907698301254788230992037436517554089813153731901053714740630067777632290263279044109132039971939418523631003901398729470133194564294739008098257736541618075899941591828471145093666305309354391035623195446166276970608181731015967624790233874366370270326864479760932900524203231830184278436814191003562432635983384281929421670375918611634407317797316652799295165911519525099501797883351588258016567922024918789788659384889579872279220173816130888835184694246374281827983237560745155075753258949227283053562501910628500151297704866094121969437077370160708734441118466080368704626746002489527720382232311900425002165641953188108360337855686670501513564687457256052331781690375766450321426273353633403142477513909384751571833511924782572723809017907023689833674706996370599120179698967425171538231959024153779574896662095855096003607316456790921677031404065249959572962292175100914700359815035888593984257853157505093625729232270135123535407728018018922505191786002181969664288859156956641790303840315765794404253062946473476459653664545480735284284134059335669012252627957917480121109594020708677634550193224362687444393161470977968751293564128510123615638346130180849773029144938312812688716193600394909160031150308199697509157166074922266173045402596199223635297610636166034729531091063293777944812539360482486391733998426813276692807973699892383126704097286777485509842526741133136179305661341101039581542183049e+3013653491"),
        ("-c4eea7fd3d243a17p-39", "1.6215967349865277787776102217957219753386438297349927591883356927621249642377948999150694600926256239077602523908902216576355908944136724776760359831455074059732687516574443913692637885779067123679536679225806571313309526434877058021244628107498650062076639399143867154249764632747842220233624106426979091073864688833452342365836878932821758972338347796873886928286135066487752284414074002956721595301982923764833887350506524033539659919870965867141659976362078645431526708931354295193907756236355831341438256650434112348763026043459529764589930622151732253102170678974383131145030216621198045701417293147654358987313409412548558056857161873845467555206127441297837892566229336270020394789811565681634448290050196510245240917572980320419890071551833132104222406285242812272818379806940327405738354230760653022015846186859535100052864082781645273163107201021408017564489420604515587591143310262580489396533458247343140033466399965357365590165369214282193896962383231612359540067146668111858530958809239012154339732495118732312132594048145284637294466979304005260267838001810948400280001478777636762192646289486415401038732622325517474816737663870766096314143995443093779163655923679741848881269962288434520365955499403246058906068128001292507470877140057134684452407661392351158130995500462098751787241224273373265764552594334100414013811184381644156226917543628049800746178574922008398524840445163820437678235523582423594330815713669128648441082300124039276462433044024758525240685071815748777548010036943550939455950239462439506365530182304211387006138914246056833377090892291041456062125522329036093847891808645844844352637355545200421639906484548477151210757413982866497429666154891398833452702784299899025338533322239524586592193769615177974891491544684305621190800233273690375463845065717271735822101362995022473039950599219482902265280545375930866645513199139668096579222143742953560459976261098568049306795225329300242969907347000129545074260083186523281911158518101484270298647903043554849316408311295208055135260297731963215477920762318479197119314275835976305208238445966382522292278537873260063795286669317841054862013587181924683949429381378477464748013802714876634541808109594768296762166648301284070896546767666104439072143261044619487806188654769424971635020305563946800945430258973923352055274659844440888876839264656703351084805035823498658567694357069599072121594841563245974769079558291296297293099129979021793090032702513675431341362945754153582729905673150599148649958734274777904887002447999141095494593211611347178883590698816824477286120930500365101137935135123880229306119151646600307527924158305827137788993202311564230677326643890830981712126276558691215790416593471866439854866276550715523189591806934431921034736213779215e+159503325"),
        ("-cf3a4dee1dd93d1394d76afacb649087p-102", "1.0135845124423200089746634347468035740847326911199961175147352551721809537947827698337935832703362490063169946647228465204128353795492490191877882479159540860891582051633544763195564624751129352764961317776162281664924853365059310998173975050775619298654346590595650352608566161422896538324939614487830331721104373308880179030949063339237278179786185670973699882918726570603369354610227572365052767198390515873765293924889404313594662759602176340340823486094489227255876807545001896169351897197108359194558704822452329577543334270168374963227170548586948984723194333624305964365496958797040772129796438424109118667082290190995996667043507890200862026796236037639862302327777021611660732978175001590972879957816832125755192809830018160484603773936223170013602999409325315451771480712081431746226631598730154880156701302907978698684814228957268482893221109334716034533482025275438441950160319914996875052741003459902364988024560222736864160310734992904418546387459717957942521603139868687911453707401365769124564593829245296766608516497440597159248679334021987782129669037759748407042140791517608917403032636220693166065335108203650422120557556033820073832724134894159002049539809544868062942867346238591737613358648035278453168211077191364505933979912951744466760637059531079595719950246238849365172357040803755274286419160066434410537552619241544587357578789496417251237047246951308566462425823268486047458016883261005967864977795380582083562361937506270620909443264917468933118108786468398897940452657221013367310463779302805327492448933905419059967361295501062957819506541721636594193258634197022810055122788096231919428102740999970328298917423047824873670150555570521166280509640698218317242647833987890044042078737749343458692683293226362568931119541040114684330622514303857199601807435403388083212920376908738938931107756369036222423273059517185422744557789409695990961760994270383046204715025509484049941696670164069078495532741968982386840621786596576600269386086186690504567851146907693601036041965950642080530578034597431079303533742331511771081989025579302203661337817866295844748608812638840708066385932519615820082787582434172479203014556604313749685418716350309124301425618627788839665082952851652605259294422181014110014697645984642788114437904994843549067716020417196262872934493803984971969687188776648187539437746793645394587581445695196929765562643083003313699426222410012721789980637982981746617667650767783974910635501059614834668748131834356768135106357600216573330505372113739885286264683621751718111497387510659720748879741887672834382826007513188464885146438742618077050990617269548168042919277227586959310711839615900155169517607463358221129814951198108516778776686370260623356362059371781520437071716817780303728875874156825e+353239371"),
        ("-dae5c6cc8145cf498b76386c093258c1dbb5f1380e914f256d260588e9d73bf5p-209", "4.6578333898502516235016740535518300122194407894176245283755932619792632299648119596141964517042595446203848479133493646182356386875578525070606477792322856286319954904373292926001762615003461298198220898443963143292606533594625406713452415250754737618870396837300485251179827942754722675149816928170712867872158574145149597338554048669138727253635567043872334715696479239439053456917329341801850295868605380366248478188598243055914348107716943651405525596694302369390107107604835494299535910847188547528430749747662916239887992490555590395024491865753234885459682016971178646169821560778042495595746260669710318636756386156900299855229984506883991117294405690784156537647461033386393639673175662945327807487007441861441771983089506861247956663677655097225378794025128148065924841638577577456703947283603781940424573992191915055640821712825736697426659554972878176421498786706545247906779064431758109339908660297499054820385073521229206112116002639186912815544940444334371618187918816511241688770523812875816372580460282305003902221440357655020515222929880746914904931114043850745918719370331150010112429760749956792459006647946710132184522036292784161738201431579429265188062975974474949916398288101701489189948649712707284237270937997436592704111662692104892913853574508057446970350304928316211466261883965957754435476698881275561000952191923954480033706983533555595042513339708733669471911511679709712210515555374889482732591405919640879270100432025020646240258197032335366878351970932078989356477397074163246963571694323146934894254286949987005394557304529770945635511710529037275336032396177049739137305737047545359333596075478270580846767181290012543395480771082161647119744445352363894006460496047787432441705816241280590105018197508953156567865398778749698012187767156157608555658716718879552183073972163907866245536333262040913176734004015777821182534541176502182491125542628452938982061668608062176578419455586040631215584630249763591070805863079409506146513739552645092420953712482111014736298476475544042448260296219359938393062779190692805849120246175150528416839014012793950441795491561716053102865706285201872811477266527155374511021774512245716413075121212846074821178874269360432191797182368925817545180936416946704808726222057902401487390534478395012633086558243654382310390629803214773698457207221837372121185409243843394417160093001231357892329646544562387717246208259069256578280586791976958797080444661721573154171843758720046343238529504173304603226105712623838056721137515973888719918260754049879688855524081248722177403655415052426508225978745438584629738298914023803429837081634582390709122597109335629142741842072100651271178682858647510850098255560946443450402441403595993409577382159774717290334103913369746857773813898601e+1546123350659986"),
        ("-e9d67343faafcd05p-37", "2.1855058594798886781599558282282614306193656231738566329257838750701785993087691412051533784159743776502811717514953779259556196436792789545040659861075478717717721441794924186048058060179155333897831985115426919931656897067991265965275406231927888661361315033698413986236873444822247873844863170472009122193473279917818079651483519714916438606226717033465988048705354623959333527074139494471269965747398889275892752793353488724725662331756649342590733198078421551403591437734423352729910794143915102183441672251877423647594546163229195810030160325334403460265916250981757515988383237997847744885973231993292679439919931743174833164272044684350790190422262117167199571470858380168782561862392524416724663506348926468635729958989020426579421717377762153255658467136042197783035116575309071032860471787393720026626689959435759983424237822359926756701867042417967149179850545554976652894341670863589427964583544131934173045997925866701366230499898126630909918665095426991822790372011349614157307554923888362052844245461322495275334659518322209124971538189598059209666443746726300838707011398512876413086641508984195339715933148512467832935624968432158069654027326115803724587702686840016661373338259091447516519138564114600205753974838379999610862256693720217751654546052511416278133487015008756855945735742548177873262948465371525595348717846341029243069073556777567456140344602907063859296076518310761290560522582085624343131211824512856947279882925769868980680007781219800575666261785523582664661100481212073019328164888413385178128410293979497392915479266314189703289399124622709139349068592379883559265931448467601919253235496564725501090306911047736367219213137154585388789957451989915141770546012943773147214877300338236533040707284498090377716863527398828537697163580915299946610181077922093245573349445068221014208822694636605859422574077427023849983914336258776491468621990639732023884736451054626097867693484674446800327982072556288418420814896112698050588709973044397832163466713186918866507401576829657216014332993042730290787228872365536886620762488933388335989082729894180510546996807910999973360590685705656525103950650568694377126466397114552589577238962016414697472964914991073719073745240292365837170496089883133787091189928638066821520929540387038767818158629024325899328141153669160140784245510534331193971085528434771487057147670593908359464397649136789020178426649197483656554114879070200495525114184512550845305217250883086297499578512557747884413055834234465466957590436114568214478545664448067833830166830067537983404421286078094380283785537325405648355435399393012294469675542714509311740070586814751984168273860404716612912270985210278526133708912969192519834978047429649662176888233595860295592088148004122060878310168129404e+840535172"),
    ];
    for &(ss, want) in z_rows {
        let x = parse_real(ss).expect("s");
        for p in [64usize, 128, 256] {
            for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
                let v = x.riemann_zeta(p, rm, &mut cc);
                let w = cr_ref(want, p, rm, &mut cc);
                assert_eq!(
                    v.compare(&w),
                    Some(0),
                    "zeta({ss}) p={p} {rm:?}: {v} {:?} want {w}",
                    v.err()
                );
            }
        }
    }
}

/// lambert_wm1(±0) is NaN(InvalidArgument). Jeff wrote "Lambert W staying NaN at
/// zero" (2026-11-15); 6227d48's −∞ was a mislabelled ruling and is reverted.
#[test]
fn zf_31_lambert_wm1_zero_is_nan() {
    let mut cc = Consts::new().expect("Consts");
    for xs in ["0", "-0"] {
        let x = parse_real(xs).expect("x");
        for p in [64usize, 128, 1024] {
            for rm in MODES {
                let v = x.lambert_wm1(p, rm, &mut cc);
                assert!(
                    v.is_nan() && v.err() == Some(Error::InvalidArgument),
                    "W-1({xs}) p={p} {rm:?}: {v} {:?}",
                    v.err()
                );
            }
        }
    }
}

/// ZF-35: the Clausen Bernoulli series stopped at 80 terms (~173 bits at π); near ±π its absolute error swamped the
/// tiny Cl2 value. The same error corrupted Barnes G at negative half-integers through the reflection formula.
/// FLINT polylog(exp(iθ)), acb.barnes_g; all reference balls certify >8700 relative bits. Logs: land39/zf35.
#[test]
fn zf_35_clausen_relative_accuracy_and_barnes_half_integers() {
    let mut cc = Consts::new().expect("Consts");
    for l in include_str!("tier1/zf35.tsv").lines() {
        let t: Vec<&str> = l.split('\t').collect();
        let x = parse_real(t[1]).expect("x");
        for p in t[2].split(',').map(|p| p.parse::<usize>().expect("p")) {
            for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
                let v = match t[0] {
                    "cl2" => x.clausen_cl2(p, rm, &mut cc),
                    "cl3" => x.clausen_cl3(p, rm, &mut cc),
                    "bg" => x.barnes_g(p, rm, &mut cc),
                    _ => panic!("function"),
                };
                let w = cr_ref(t[3], p, rm, &mut cc);
                assert_eq!(
                    v.compare(&w),
                    Some(0),
                    "{}({}) p={p} {rm:?}: {v} {:?} want {w}",
                    t[0],
                    t[1],
                    v.err()
                );
            }
        }
    }
}

/// ZF-19: complex asin/acos/asinh/atanh on the real and imaginary axes (off the cuts) went through the log
/// formulas and returned rounding noise such as 2^-191 where the exact part is a signed zero; elliptic_pi(0, x, 0)
/// inherits it through asin. z^k for a diagonal base |re z| = |im z| fell to exp(k ln z) whenever the exact
/// power path's width limit tripped (at p = 64 already for (1+i)^-2), with the same noise in the zero part.
/// Now: the axis cases are the real function of the nonzero part plus a zero signed by the C99 Annex G
/// symmetries, and diagonal powers are r^k 2^⌊k/2⌋ placed by k·arg mod 2π. `*0` is a zero of either sign
/// (C99 leaves cpow's unspecified). c.sqrt off the axes took cos(θ/2) at θ ≈ π and lost the small part
/// (sqrt(−4 + 2^-200 i) re = 2^-202(1 + 2^-117) at p = 256); now the Kahan form under Ziv.
/// FLINT, >3000 relative bits. Logs: land39/zf19.
#[test]
fn zf_19_exact_zero_parts_are_signed_zeros() {
    let mut cc = Consts::new().expect("Consts");
    let cx = |t: &str| {
        let (a, b) = t.split_once(',').expect("re,im");
        ExactComplex::new(parse_real(a).expect("re"), parse_real(b).expect("im"))
    };
    for l in include_str!("tier1/zf19.tsv").lines() {
        let t: Vec<&str> = l.split('\t').collect();
        let a: Vec<ExactComplex> = t[2].split(';').map(cx).collect();
        for p in t[1].split(',').map(|p| p.parse::<usize>().expect("p")) {
            for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
                let v = match t[0] {
                    "c.sqrt" => a[0].sqrt(p, rm, &mut cc),
                    "c.asin" => a[0].asin(p, rm, &mut cc),
                    "c.acos" => a[0].acos(p, rm, &mut cc),
                    "c.asinh" => a[0].asinh(p, rm, &mut cc),
                    "c.atanh" => a[0].atanh(p, rm, &mut cc),
                    "c.pow" => a[0].pow(&a[1], p, rm, &mut cc),
                    "c.elliptic_pi" => a[0].elliptic_pi(&a[1], &a[2], p, rm, &mut cc),
                    _ => panic!("function"),
                };
                for (part, want, name) in [(v.re(), t[3], "re"), (v.im(), t[4], "im")] {
                    let ok = match want {
                        "+0" => part.is_zero() && !part.is_negative(),
                        "-0" => part.is_zero() && part.is_negative(),
                        "*0" => part.is_zero(),
                        _ => part.compare(&cr_ref(want, p, rm, &mut cc)) == Some(0),
                    };
                    assert!(
                        ok,
                        "{}({}) p={p} {rm:?} {name}: {part} want {:.40}",
                        t[0], t[2], want
                    );
                }
            }
        }
    }
}

/// ZF-16 (rem_pi): the reduction was exact enough but returned the full working width (about p + exponent bits)
/// and ignored `rm`. It is now x − 2πk (k toward zero) correctly rounded to `p`, with the width grown when x is
/// close to a multiple of 2π. FLINT, >8700 relative bits; the near-multiples are the 128-bit dyadics just above
/// 355·2π, 710·2π and (10^12+39)·2π. (Trig and Clausen rows of ZF-16 pass; see zf_35.)
#[test]
fn zf_16_rem_pi_correctly_rounded() {
    let mut cc = Consts::new().expect("Consts");
    let rows: &[(&str, &str)] = &[
        ("1p70000", "1.3669972433977399031113762695610780923659018123264177754142735792976289472451179720971505969858167329171977151151796933865912634426391283316018763494252408487044505699728884926075090372772709503804421160409196588739240672904173258221171922346172530104441704933948050199185739310163937735598605840491325029574931647010922178749717200883577486848142489654483932162004747810125460239525738763801123936678364431219302095049208850111403464894920721023453248179316348974787755392664716145520225930671864724094932705819974683198601380948356084328569780977373101438753234840001551812368652523665947553663702457615845135684919365901149727685086868869536645801045070266188154500514549625647931531948920477843902207940811191418317881131322690067652818570971787396285761800972312912543796759480418127628932176860803376550371760401097475954798650137637828094603561942769866270667935308259974887706604327921348052238874210701657016259013936265037026019508239251141578780783609799364685744553687590358636159117967929107412434941987661056865680943747807003109918908718927006455591750486734388046279879806678107003221513636482642676334111227540635447766043448965449473715061260107900029863728408302549890821289837598664108518328020269186858409021493412376653092060004896690768534806141454933457837246540419292559449236284934100848483680040096364744331461996132391151451952510734537946845010743134431577818265769714492189901766223242671280248137311572560378767373780758111871315817174257560450356227416639459634439117866928997418461355393686902224983309571231714717071158474626530385064145864162167660851806786818451321359762840620131830517578476479388645290907018700366448548103932084455128146107321252247018399994929691444389939988165650044280233846020470270332033626083858939186479314723357531380531655321732741086320370959895642001870534386512742401211625448862176787093827853186611197423124798984446640609206295799919124794542391553025066806752622469689695196989247602638912076082328147984701597925940688099829091371776931084532139738150712274801496469024148066384909374126967890440013297217047353618496281435248461412993987399752359243195557814495940971913800932364397516535365491751393901820826474669508650357101871885757434174514208113474130872081163478840687033975642377465535448190033785691365043245236419285747509965649927236005007056030748107024572032311843748644471552224425101972370855910423077077862548407142592677342689586544102208499639545683081021464181296541885396038316341560793239386352527446420851370217529489705776637443284091946860118004836891465419393439199118534704190167171236797970826430878313293638789024424796481789227405880285587699546381322534936972160137045198409962982257526241284482409949879438621268827163874364542839521221166346318"),
        ("1p100000", "3.5500407910754740683254611430587799296401550995038442093781511890081852023791757609104464353969822404658545048672361443946865273796693387625220939392675429147541840394781454731830262618751560549434708942058780625320366060590922370110299346596694661330089538569524337725612044417854924791302647801753790291469218313804603184806799553282673615993182495659736553489250246491716695153815527596410821267939479056489003735879161028276064083651147142827406014396492656477487266080767538584143449496232850622326579585486688716347458993080080422928035900786536059261616630201121154256678553072157951247363876601185774734265928860750293037020220381784858670660314869628819244909421442933327101303610250506288682001805886054938470153687695359995402973108543952809325443953328934633698823234868615643609705395656959593319125566401630299194705353622105062114472426490043198056049343341239972137571096334964881861074503246937007707714829579887191748253820989547987509787047745895801214331578349572981454406126645545223538075003605974691675025657309966750814030536072126623411533071992643241252117819900901149877889595329925536436024074913097940205646012838435279280610118624024245499787425397829290594850698503576011117783158579528378767153815166521262229184286420205934776266186329890580599811041332494806420452579545159406551254796550834169603393081461417873000282550322481674941398134865939815693678082780027805277309095091711438955185324637863480886060796959489397981444284229812067202661140341434949053289749011150495832669680307348420448783218145355165666825810445870998958941165128070217851796157592789813461253221319574272052190125375633808293843047531547639171196158186076029190183497564873661299417911189722386755506321287787867533631158592622482540808189821258156698859216062972216982766464410961126965878066352590486269114713068005527394230123950378085514507877626757408657571786945258112548099459961356283603367227962321311222257574338401201595588899019819974047269396565084167122513525597386550445931774956502609172877153704100924506565034177299822776186159041924928297105180407591713949350906530594144185360512129327966433469651883468971933171844864779148964786150063964566103263502286679105204691063506800748823622433216270666213971828175750172554677759900684982567670764174340055388198227391148941632502908332446109259803334540183293806274401991413361677537334393223289982944524947836587248020939503989856104753338324430250835099393869190251830575865047412378165635630401264719254298764542424592904181027837007479427496933816858750558354734758827227565103736276259533868324958052049032196346601707356835762688607821998224239045693940205789190649383315535207688432860265386368858031394718602460292046230349471336317229503434454261523816308260847002"),
        ("8b687e176a2a22a7966f8d53b9d89dfcp-116", "1.1758047903021125405882129245252032043304443265672116187899368569585106378303937280563370327219086615474994820228778406610298287081623650015052652317706241999428805255544207277187253248782289754011383952657446644506589517484343152486577428079768490538845223109695655543225309515932395834120146116391229594882091298796260844391955323310337468769817255664121526757302529338040310689713304547856658827845626895365830867131710226447807180161115565327228163890061516117818091912667133783408924809588059390436594601452298729971032917155804831844505806880878887954722356279596629726862871812272102346641969377486717603107757406960130283064595949794386553937481420626164285165831160388394407786535563106532322859912369928031717879941500001155185584534747590469700571613797532633786703695105123809161140125273624964028856398582897772134485450522851422049191361803854905414445240021795583617200876710366643528810808436109974880975822158428967304768631207606050904456341658753952367387243505037072005185331098902142556485614054935684392013000727954416137212697539548729816623941851357857553005874911031487028173222318756680148954207399893083328907115341331813491303007334950870563261469767561325655770820501526698747794835410891738680359943938556126250401192286854176060661181747289309417702005625964939554122918792913362607010098474674732775424277503982620918516450968354814284528105038433367377634733853543716328496810088065448207614661098842103281833059441680333768394901721972605059086043268232098727495624883897293594967877436727154805217588208923521982247405781316401942030587388536873073557015777511631399576198766330805433578279687985232223286877408648440025750158232471088684027611521665665250909554352249166489319477609753028748832826038568458074664717831879093087482142736556167588322742511683725017154400161104366650980723780268773074868620988204647634230430557722211562178352171589016542163211257720981295938725897887093122486686137485258821342972492760600411029379307503944598677437169657046023665939429032520514029170429377352951336459125503159949406733401575517660759257926247483684835474631931834572595512323846000544938475500968561187395759585154710787621829234453050630349245782293075372688803954519507476599349035712265947455775917281935469621528283616842082779616663049668425633621859080441632907972576496396643325347631489070328719661773059807496219760728041654988890097947984370981648989532186623244969615816757388030776290005886394728519456596229880114295990669589721416773242087495405265866574912824267588252202327932650881270636546178919748344293601803102414130050464430459654855877852776605070746406912248551129549867960935856010454569418140173601174570504195361233847053689638137624938589466712067483461629275675503735694597151796507e-35"),
        ("-8b687e176a2a22a7966f8d53b9d89dfcp-116", "-1.1758047903021125405882129245252032043304443265672116187899368569585106378303937280563370327219086615474994820228778406610298287081623650015052652317706241999428805255544207277187253248782289754011383952657446644506589517484343152486577428079768490538845223109695655543225309515932395834120146116391229594882091298796260844391955323310337468769817255664121526757302529338040310689713304547856658827845626895365830867131710226447807180161115565327228163890061516117818091912667133783408924809588059390436594601452298729971032917155804831844505806880878887954722356279596629726862871812272102346641969377486717603107757406960130283064595949794386553937481420626164285165831160388394407786535563106532322859912369928031717879941500001155185584534747590469700571613797532633786703695105123809161140125273624964028856398582897772134485450522851422049191361803854905414445240021795583617200876710366643528810808436109974880975822158428967304768631207606050904456341658753952367387243505037072005185331098902142556485614054935684392013000727954416137212697539548729816623941851357857553005874911031487028173222318756680148954207399893083328907115341331813491303007334950870563261469767561325655770820501526698747794835410891738680359943938556126250401192286854176060661181747289309417702005625964939554122918792913362607010098474674732775424277503982620918516450968354814284528105038433367377634733853543716328496810088065448207614661098842103281833059441680333768394901721972605059086043268232098727495624883897293594967877436727154805217588208923521982247405781316401942030587388536873073557015777511631399576198766330805433578279687985232223286877408648440025750158232471088684027611521665665250909554352249166489319477609753028748832826038568458074664717831879093087482142736556167588322742511683725017154400161104366650980723780268773074868620988204647634230430557722211562178352171589016542163211257720981295938725897887093122486686137485258821342972492760600411029379307503944598677437169657046023665939429032520514029170429377352951336459125503159949406733401575517660759257926247483684835474631931834572595512323846000544938475500968561187395759585154710787621829234453050630349245782293075372688803954519507476599349035712265947455775917281935469621528283616842082779616663049668425633621859080441632907972576496396643325347631489070328719661773059807496219760728041654988890097947984370981648989532186623244969615816757388030776290005886394728519456596229880114295990669589721416773242087495405265866574912824267588252202327932650881270636546178919748344293601803102414130050464430459654855877852776605070746406912248551129549867960935856010454569418140173601174570504195361233847053689638137624938589466712067483461629275675503735694597151796507e-35"),
        ("8b687e176a2a22a7966f8d53b9d89dfcp-115", "2.3516095806042250811764258490504064086608886531344232375798737139170212756607874561126740654438173230949989640457556813220596574163247300030105304635412483998857610511088414554374506497564579508022767905314893289013179034968686304973154856159536981077690446219391311086450619031864791668240292232782459189764182597592521688783910646620674937539634511328243053514605058676080621379426609095713317655691253790731661734263420452895614360322231130654456327780123032235636183825334267566817849619176118780873189202904597459942065834311609663689011613761757775909444712559193259453725743624544204693283938754973435206215514813920260566129191899588773107874962841252328570331662320776788815573071126213064645719824739856063435759883000002310371169069495180939401143227595065267573407390210247618322280250547249928057712797165795544268970901045702844098382723607709810828890480043591167234401753420733287057621616872219949761951644316857934609537262415212101808912683317507904734774487010074144010370662197804285112971228109871368784026001455908832274425395079097459633247883702715715106011749822062974056346444637513360297908414799786166657814230682663626982606014669901741126522939535122651311541641003053397495589670821783477360719887877112252500802384573708352121322363494578618835404011251929879108245837585826725214020196949349465550848555007965241837032901936709628569056210076866734755269467707087432656993620176130896415229322197684206563666118883360667536789803443945210118172086536464197454991249767794587189935754873454309610435176417847043964494811562632803884061174777073746147114031555023262799152397532661610867156559375970464446573754817296880051500316464942177368055223043331330501819108704498332978638955219506057497665652077136916149329435663758186174964285473112335176645485023367450034308800322208733301961447560537546149737241976409295268460861115444423124356704343178033084326422515441962591877451795774186244973372274970517642685944985521200822058758615007889197354874339314092047331878858065041028058340858754705902672918251006319898813466803151035321518515852494967369670949263863669145191024647692001089876951001937122374791519170309421575243658468906101260698491564586150745377607909039014953198698071424531894911551834563870939243056567233684165559233326099336851267243718160883265815945152992793286650695262978140657439323546119614992439521456083309977780195895968741963297979064373246489939231633514776061552580011772789457038913192459760228591981339179442833546484174990810531733149825648535176504404655865301762541273092357839496688587203606204828260100928860919309711755705553210141492813824497102259099735921871712020909138836280347202349141008390722467694107379276275249877178933424134966923258551351007471389194303593013e-35"),
        ("-8b687e176a2a22a7966f8d53b9d89dfcp-115", "-2.3516095806042250811764258490504064086608886531344232375798737139170212756607874561126740654438173230949989640457556813220596574163247300030105304635412483998857610511088414554374506497564579508022767905314893289013179034968686304973154856159536981077690446219391311086450619031864791668240292232782459189764182597592521688783910646620674937539634511328243053514605058676080621379426609095713317655691253790731661734263420452895614360322231130654456327780123032235636183825334267566817849619176118780873189202904597459942065834311609663689011613761757775909444712559193259453725743624544204693283938754973435206215514813920260566129191899588773107874962841252328570331662320776788815573071126213064645719824739856063435759883000002310371169069495180939401143227595065267573407390210247618322280250547249928057712797165795544268970901045702844098382723607709810828890480043591167234401753420733287057621616872219949761951644316857934609537262415212101808912683317507904734774487010074144010370662197804285112971228109871368784026001455908832274425395079097459633247883702715715106011749822062974056346444637513360297908414799786166657814230682663626982606014669901741126522939535122651311541641003053397495589670821783477360719887877112252500802384573708352121322363494578618835404011251929879108245837585826725214020196949349465550848555007965241837032901936709628569056210076866734755269467707087432656993620176130896415229322197684206563666118883360667536789803443945210118172086536464197454991249767794587189935754873454309610435176417847043964494811562632803884061174777073746147114031555023262799152397532661610867156559375970464446573754817296880051500316464942177368055223043331330501819108704498332978638955219506057497665652077136916149329435663758186174964285473112335176645485023367450034308800322208733301961447560537546149737241976409295268460861115444423124356704343178033084326422515441962591877451795774186244973372274970517642685944985521200822058758615007889197354874339314092047331878858065041028058340858754705902672918251006319898813466803151035321518515852494967369670949263863669145191024647692001089876951001937122374791519170309421575243658468906101260698491564586150745377607909039014953198698071424531894911551834563870939243056567233684165559233326099336851267243718160883265815945152992793286650695262978140657439323546119614992439521456083309977780195895968741963297979064373246489939231633514776061552580011772789457038913192459760228591981339179442833546484174990810531733149825648535176504404655865301762541273092357839496688587203606204828260100928860919309711755705553210141492813824497102259099735921871712020909138836280347202349141008390722467694107379276275249877178933424134966923258551351007471389194303593013e-35"),
        ("b6dd60e7e4142eb9f4d5d94e3c327c2fp-85", "1.0568293248677303508001613999654076848129124611146940977466647659018695670827835957776290862425515657004335769668459731265410267435000797778129506107868437858646052657808255442685851641105894734405333107222227094817549117239718714559061439511957775306918857016086251392967436851037134944342280433582859981468589445061020859281716890821475285574423709707449013318882447553732702349005286379592894944259257999019183988875993324692799659134028296770455411092460322562057721242991660191557585991946499881719799768110777237824429767956831311587555479782576378656553895520650145413428564427097899569917130616917618053041519418257275554749411354461706467547781679520355100197891739464399501525867450732978055833991916956530423400379540840768743900224788248638185148287650295568238894275557197896535887193889885715252190907549361046983713982073510858083331175390910897252945299529461825290284237947711204846521429672489003971786228672639505264640453330963457807191486321204928710080697598211759160230357484801043575613514407432741803254908223219501874937802165915563438496738913256876859743186917179482111840705106635030477142867396344499767613931242701946449619572003131515068519869632474712610438625921413410555927821583255597919393582386169476236633468198929077633376041165225621143228982234522021897283243480349736450517520046051060889962231217532866876540206152611143302236458578957643589108831496539354785309901478526264228307643884335422717430708365272741089071647151943679510543591170525223704103314744421473142084146160368738099295626045991571612381473629926903570427963561341030934394236133465788774363809551871203339338538309562775017234652015089093647006894066127731986180236564345752395083141681787715800161112965736586973558293407811093775123627314904273473029372864698510295911079007403695704535395596544755654264208675090200851546769403167672824705385395654079870223995556377386695366712402437595873210828643461136429589528927499266424520838272656271375343782319814339131465045095790045276852293588607020145806869660350677741076551804521920719158751097471493699776628493005900687208807029973008685185503865287959226953435273862427413897277384403226511106656211747329494346238751765340727318968715646115016410168527604621295348040251359748536325419061780511791100488744445277395753687120225782775808760000506500207911850987545799679852776101070068465828431606744297626151330964932633244001070423861785455563082037496346613041769405877943771124016875048278652250580264921873388414764007745244245856044276600180309443235954576019076225380086052041052521891470390224970705271714314173585827227547368969878886463408472284504584680154896622143467796626895498052062076829975125305663866979600486531928103757984992880741971570508471204073075485883369e-26"),
        ("-b6dd60e7e4142eb9f4d5d94e3c327c2fp-85", "-1.0568293248677303508001613999654076848129124611146940977466647659018695670827835957776290862425515657004335769668459731265410267435000797778129506107868437858646052657808255442685851641105894734405333107222227094817549117239718714559061439511957775306918857016086251392967436851037134944342280433582859981468589445061020859281716890821475285574423709707449013318882447553732702349005286379592894944259257999019183988875993324692799659134028296770455411092460322562057721242991660191557585991946499881719799768110777237824429767956831311587555479782576378656553895520650145413428564427097899569917130616917618053041519418257275554749411354461706467547781679520355100197891739464399501525867450732978055833991916956530423400379540840768743900224788248638185148287650295568238894275557197896535887193889885715252190907549361046983713982073510858083331175390910897252945299529461825290284237947711204846521429672489003971786228672639505264640453330963457807191486321204928710080697598211759160230357484801043575613514407432741803254908223219501874937802165915563438496738913256876859743186917179482111840705106635030477142867396344499767613931242701946449619572003131515068519869632474712610438625921413410555927821583255597919393582386169476236633468198929077633376041165225621143228982234522021897283243480349736450517520046051060889962231217532866876540206152611143302236458578957643589108831496539354785309901478526264228307643884335422717430708365272741089071647151943679510543591170525223704103314744421473142084146160368738099295626045991571612381473629926903570427963561341030934394236133465788774363809551871203339338538309562775017234652015089093647006894066127731986180236564345752395083141681787715800161112965736586973558293407811093775123627314904273473029372864698510295911079007403695704535395596544755654264208675090200851546769403167672824705385395654079870223995556377386695366712402437595873210828643461136429589528927499266424520838272656271375343782319814339131465045095790045276852293588607020145806869660350677741076551804521920719158751097471493699776628493005900687208807029973008685185503865287959226953435273862427413897277384403226511106656211747329494346238751765340727318968715646115016410168527604621295348040251359748536325419061780511791100488744445277395753687120225782775808760000506500207911850987545799679852776101070068465828431606744297626151330964932633244001070423861785455563082037496346613041769405877943771124016875048278652250580264921873388414764007745244245856044276600180309443235954576019076225380086052041052521891470390224970705271714314173585827227547368969878886463408472284504584680154896622143467796626895498052062076829975125305663866979600486531928103757984992880741971570508471204073075485883369e-26"),
    ];
    for &(xs, want) in rows {
        let x = parse_real(xs).expect("x");
        for p in [64usize, 256, 1024] {
            for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
                let v = x.rem_pi(p, rm, &mut cc);
                let w = cr_ref(want, p, rm, &mut cc);
                assert_eq!(
                    v.compare(&w),
                    Some(0),
                    "rem_pi({xs}) p={p} {rm:?}: {v} want {w}"
                );
                assert!(
                    v.precision().unwrap_or(0) <= p,
                    "rem_pi({xs}) p={p}: width {:?}",
                    v.precision()
                );
            }
        }
    }
}

/// ZF-37 (PrecisionRetryExhausted on valid inputs): each root cause was a value Ziv can never
/// certify or a series that gave up early, not a precision shortage. Exact or ±1/x-plus-tiny
/// results are now rounded once (expm1(x≪0) = −1 + tiny, si/C(x) = x ∓ tiny, si(maxv) = π/2,
/// erf(2^20+1 + i/2) = 1 − tiny, I_x(a, b) = 1 − tiny near x = 1); exact closed forms are summed
/// exactly (nth_root of a subnormal, Γ(n) for any n, B(a, n), I_x(a, n), Gauss and Euler ₂F₁);
/// ei/li and γ(s, x) series lost their term caps and the wrong Γ(s) shortcut; K_ν(x) at ν ≈ x
/// falls back from the Hankel series; Ci(x ± 0i) pins the zero imaginary part. Expected values
/// are FLINT-certified (rjudge/cjudge/erfbig, 2p+64 bits) except I_{1−2^-256}(½, 100), whose
/// complement I_{2^-256}(100, ½) < 2^-25000 is far below an ulp (land39/zf37).
#[test]
fn zf_37_valid_inputs_do_not_exhaust_precision() {
    let mut cc = Consts::new().expect("Consts");
    let mode = |s: &str| match s {
        "up" => RoundingMode::Up,
        "down" => RoundingMode::Down,
        "tozero" => RoundingMode::ToZero,
        "fromzero" => RoundingMode::FromZero,
        "toeven" => RoundingMode::ToEven,
        "toodd" => RoundingMode::ToOdd,
        _ => RoundingMode::None,
    };
    let cx = |t: &str| {
        let (a, b) = t.split_once(',').expect("re,im");
        ExactComplex::new(parse_real(a).expect("re"), parse_real(b).expect("im"))
    };
    // Runner format: F:sign:exponent:hex (value 0.hex · 2^exponent) or Z:sign.
    let check = |v: &ExactNum, want: &str, what: &str| {
        let k: Vec<&str> = want.split(':').collect();
        match k[0] {
            "Z" => assert!(
                v.is_zero() && v.is_negative() == (k[1] == "-"),
                "{what}: {v} want {want}"
            ),
            "F" => {
                let e: i64 = k[2].parse().expect("exp");
                let w = parse_real(&format!("{}{}p{}", k[1], k[3], e - 4 * k[3].len() as i64))
                    .expect("want");
                assert_eq!(v.compare(&w), Some(0), "{what}: {v} want {want}");
            }
            _ => panic!("expected kind {want}"),
        }
    };
    for l in include_str!("tier1/zf37.tsv").lines() {
        let t: Vec<&str> = l.split('\t').collect();
        let (f, p, rm) = (t[0], t[1].parse::<usize>().expect("p"), mode(t[2]));
        let a: Vec<&str> = t[3].split('|').collect();
        let r = |i: usize| parse_real(a[i]).expect("arg");
        let what = format!("{f}({}) p={p} {rm:?}", t[3]);
        if let Some(cf) = f.strip_prefix("c.") {
            let z = cx(a[0]);
            let v = match cf {
                "erf" => z.erf(p, rm, &mut cc),
                "ci" => z.ci(p, rm, &mut cc),
                _ => panic!("function {f}"),
            };
            let (wr, wi) = t[4].split_once(';').expect("re;im");
            check(v.re(), wr, &format!("{what} re"));
            check(v.im(), wi, &format!("{what} im"));
            continue;
        }
        let v = match f {
            "r.expm1" => r(0).expm1(p, rm, &mut cc),
            "r.nth_root" => r(0).nth_root(a[1].parse().expect("n"), p, rm),
            "r.si" => r(0).si(p, rm, &mut cc),
            "r.fresnel_c" => r(0).fresnel_c(p, rm, &mut cc),
            "r.ei" => r(0).ei(p, rm, &mut cc),
            "r.li" => r(0).li(p, rm, &mut cc),
            "r.gamma" => r(0).gamma(p, rm, &mut cc),
            "r.gammainc" => r(0).gammainc(&r(1), p, rm, &mut cc),
            "r.bessel_k" => r(0).bessel_k(&r(1), p, rm, &mut cc),
            "r.hypergeom_2f1" => r(0).hypergeom_2f1(&r(1), &r(2), &r(3), p, rm, &mut cc),
            "r.betainc" => r(0).betainc(&r(1), &r(2), p, rm, &mut cc),
            "r.beta" => zenith_float::beta(&r(0), &r(1), p, rm, &mut cc),
            _ => panic!("function {f}"),
        };
        check(&v, t[4], &what);
    }
}

/// ZF-65: complex J_ν / I_ν at a real negative integer order use J_{−n} = (−1)^n J_n and
/// I_{−n} = I_n (DLMF 10.4.1, 10.27.1); K_ν with Re ν < 0 uses K_{−ν} = K_ν (DLMF 10.27.3);
/// J_ν(0) = I_ν(0) = 0 for Re ν > 0 (DLMF 10.2.2). Expected values are FLINT-certified
/// (cjudge, 2p+64 bits; land39/zf65). J/I_{−1/2}(0) is a pole and stays InvalidArgument.
#[test]
fn zf_65_complex_bessel_negative_order_and_zero_argument() {
    let mut cc = Consts::new().expect("Consts");
    let mode = |s: &str| match s {
        "up" => RoundingMode::Up,
        "down" => RoundingMode::Down,
        "tozero" => RoundingMode::ToZero,
        "fromzero" => RoundingMode::FromZero,
        "toeven" => RoundingMode::ToEven,
        "toodd" => RoundingMode::ToOdd,
        _ => RoundingMode::None,
    };
    let cx = |t: &str| {
        let (a, b) = t.split_once(',').expect("re,im");
        ExactComplex::new(parse_real(a).expect("re"), parse_real(b).expect("im"))
    };
    let check = |v: &ExactNum, want: &str, what: &str| {
        let k: Vec<&str> = want.split(':').collect();
        match k[0] {
            "Z" => assert!(
                v.is_zero() && v.is_negative() == (k[1] == "-"),
                "{what}: {v} want {want}"
            ),
            "F" => {
                let e: i64 = k[2].parse().expect("exp");
                let w = parse_real(&format!("{}{}p{}", k[1], k[3], e - 4 * k[3].len() as i64))
                    .expect("want");
                assert_eq!(v.compare(&w), Some(0), "{what}: {v} want {want}");
            }
            _ => panic!("expected kind {want}"),
        }
    };
    for l in include_str!("tier1/zf65.tsv").lines() {
        let t: Vec<&str> = l.split('\t').collect();
        let (f, p, rm) = (t[0], t[1].parse::<usize>().expect("p"), mode(t[2]));
        let a: Vec<&str> = t[3].split('|').collect();
        let (z, nu) = (cx(a[0]), cx(a[1]));
        let what = format!("{f}({}) p={p} {rm:?}", t[3]);
        let v = match f {
            "c.bessel_j_nu" => z.bessel_j_nu(&nu, p, rm, &mut cc),
            "c.bessel_i" => z.bessel_i(&nu, p, rm, &mut cc),
            "c.bessel_k" => z.bessel_k(&nu, p, rm, &mut cc),
            _ => panic!("function {f}"),
        };
        let (wr, wi) = t[4].split_once(';').expect("re;im");
        check(v.re(), wr, &format!("{what} re"));
        check(v.im(), wi, &format!("{what} im"));
    }
    let half = cx("-1p-1,0");
    let zero = cx("0,0");
    // Official judge / mpmath: J_{-1/2}(0)=I_{-1/2}(0)=+∞ (DLMF 10.7.3 / 10.30.1
    // one-sided z→0+). 10.2.2 / 10.25.2 do not define a finite value at exact 0.
    // c3dfd98's InvalidArgument refusal made the 98 official-judge rows WRONG
    // (NAN_UNEXPECTED, expected Inf). FLINT 3.6 at exact 0 is indeterminate.
    for v in [
        zero.bessel_j_nu(&half, 64, RoundingMode::ToEven, &mut cc),
        zero.bessel_i(&half, 64, RoundingMode::ToEven, &mut cc),
    ] {
        assert!(
            v.re().is_inf_pos() && v.im().is_zero() && v.im().is_negative(),
            "J/I_-1/2(0+0i) got {v:?}, want +Inf-0i (imaginary zero -sgnGamma(1/2)*sgn(Im z), FLINT direction)"
        );
    }
}

/// J_{-1/2}(0) and I_{-1/2}(0) are compact +∞: official judge / mpmath
/// `besselj(-0.5,0)=besseli(-0.5,0)=+inf`, DLMF 10.7.3 / 10.30.1 one-sided
/// z→0+ (arg(0±0i)=0). This compact +Inf (a pole), the way it worked
/// before c3dfd98, is a Paul call per those DLMF limits. It is not a
/// Jeff ruling. DLMF 10.2.2 / 10.25.2 and the closed forms 10.16.1 /
/// 10.39.2 are singular at exact 0; FLINT 3.6 `acb_bessel_j` / `acb_bessel_i`
/// return an indeterminate ball, so the exact-zero gold is not a FLINT
/// midpoint. Neighbors at z=2^{-16} and 2^{-8} are FLINT 3.6 midpoints
/// uniquely correctly rounded to dest p (ctx.prec = p+128; ball radius ≪ ½ ulp).
/// Covers the 98 official WRONG rows (7 dest p × 7 modes × {J,I} at z=+0+0i). The imaginary zero is
/// −sgn Γ(1/2)·sgn(Im z)·0, i.e. +∞−0i at +0+0i and +∞+0i at +0−0i (Jeff ruling 2026-11-22 ~11:58 PM CT:
/// FLINT's direction; land39 REF-CORRECTIONS RC-14). Those 98 rows now change their imaginary zero from +0 to −0.
#[test]
fn zf_j_i_neg_half_at_zero_is_pos_inf_and_small_z_matches_flint() {
    let mut cc = Consts::new().expect("Consts");
    let half = parse_complex("-1p-1,0").expect("nu");
    let z0 = parse_complex("0,0").expect("0+0i");
    let z0m = ExactComplex::new(parse_real("0").expect("0"), signed_zero(true));
    let modes = [
        RoundingMode::Up,
        RoundingMode::Down,
        RoundingMode::ToZero,
        RoundingMode::FromZero,
        RoundingMode::ToEven,
        RoundingMode::ToOdd,
        RoundingMode::None,
    ];
    for rm in modes {
        for p in [64usize, 128, 256, 1024, 2048, 4096, 8192] {
            let j = z0.bessel_j_nu(&half, p, rm, &mut cc);
            assert!(
                j.re().is_inf_pos() && j.im().is_zero() && j.im().is_negative(),
                "J_-1/2(0+0i) p={p} {rm:?} got {j:?}, want +Inf-0i"
            );
            let jm = z0m.bessel_j_nu(&half, p, rm, &mut cc);
            assert!(
                jm.re().is_inf_pos() && jm.im().is_zero() && !jm.im().is_negative(),
                "J_-1/2(0-0i) p={p} {rm:?} got {jm:?}, want +Inf+0i"
            );
            let i = z0.bessel_i(&half, p, rm, &mut cc);
            assert!(
                i.re().is_inf_pos() && i.im().is_zero() && i.im().is_negative(),
                "I_-1/2(0+0i) p={p} {rm:?} got {i:?}, want +Inf-0i"
            );
        }
    }

    let want = |s: &str| -> ExactNum {
        let k: Vec<&str> = s.split(':').collect();
        let e: i64 = k[2].parse().expect("exp");
        parse_real(&format!("{}{}p{}", k[1], k[3], e - 4 * k[3].len() as i64)).expect("want")
    };
    let same = |v: &ExactNum, w: &ExactNum| v.compare(w) == Some(0);
    let j_golds = [
        ("1p-16", 64usize, "F:+:8:cc42299e3b916f99"),
        ("1p-16", 128, "F:+:8:cc42299e3b916f99360362080f4e5060"),
        ("1p-8", 64, "F:+:4:cc41c37d9565f503"),
        ("1p-8", 128, "F:+:4:cc41c37d9565f5028d7f933271dc9e45"),
    ];
    let i_golds = [
        ("1p-16", 64usize, "F:+:8:cc42299f07d39938"),
        ("1p-16", 128, "F:+:8:cc42299f07d39937d7b5e6708e3972fe"),
        ("1p-8", 64, "F:+:4:cc428fbfbf049746"),
        ("1p-8", 128, "F:+:4:cc428fbfbf0497465205aa2cdaafc498"),
    ];
    let rm = RoundingMode::ToEven;
    for &(zt, p, tok) in &j_golds {
        let z = parse_complex(&format!("{zt},0")).expect("z");
        let v = z.bessel_j_nu(&half, p, rm, &mut cc);
        let w = want(tok);
        assert!(
            same(v.re(), &w),
            "J_-1/2({zt}) p={p} re {} want {tok}",
            v.re()
        );
        assert!(
            v.im().is_zero() && !v.im().is_negative(),
            "J_-1/2({zt}) p={p} im {}",
            v.im()
        );
    }
    for &(zt, p, tok) in &i_golds {
        let z = parse_complex(&format!("{zt},0")).expect("z");
        let v = z.bessel_i(&half, p, rm, &mut cc);
        let w = want(tok);
        assert!(
            same(v.re(), &w),
            "I_-1/2({zt}) p={p} re {} want {tok}",
            v.re()
        );
        assert!(
            v.im().is_zero() && !v.im().is_negative(),
            "I_-1/2({zt}) p={p} im {}",
            v.im()
        );
    }
}

/// ZF-60: Jacobi functions for large u and m near 1. The AGM / descending Landen evaluation
/// claimed the full working precision, but cn and dn near their zeros (and asin near ±1 in the
/// recurrence) lose hundreds of bits, and dn = sqrt(1 - m sn^2) cancelled. The recurrence now carries
/// an error bound that the Ziv loop uses, and dn = sqrt(cn^2 + m' sn^2). Expected values are
/// FLINT-certified (am as the root of acb.elliptic_f, Newton; 2p+64 bits; land39/zf60; ToEven and Up
/// rows of the 590, all 590 judged), except
/// sn(1000 | 1) below 1 in Down/ToZero (tanh 1000 = 1 - 2e^-2000, verified by hand).
#[test]
fn zf_60_jacobi_large_argument_m_near_one() {
    let mut cc = Consts::new().expect("Consts");
    let mode = |s: &str| match s {
        "up" => RoundingMode::Up,
        "down" => RoundingMode::Down,
        "tozero" => RoundingMode::ToZero,
        "fromzero" => RoundingMode::FromZero,
        "toeven" => RoundingMode::ToEven,
        "toodd" => RoundingMode::ToOdd,
        _ => RoundingMode::None,
    };
    for l in include_str!("tier1/zf60.tsv").lines() {
        let t: Vec<&str> = l.split('\t').collect();
        let (f, p, rm) = (t[0], t[1].parse::<usize>().expect("p"), mode(t[2]));
        let a: Vec<&str> = t[3].split('|').collect();
        let (u, m) = (parse_real(a[0]).expect("u"), parse_real(a[1]).expect("m"));
        let what = format!("{f}({}) p={p} {rm:?}", t[3]);
        let v = match f {
            "r.jacobi_am" => u.jacobi_am(&m, p, rm, &mut cc),
            "r.jacobi_sn" => u.jacobi_sn(&m, p, rm, &mut cc),
            "r.jacobi_cn" => u.jacobi_cn(&m, p, rm, &mut cc),
            "r.jacobi_dn" => u.jacobi_dn(&m, p, rm, &mut cc),
            "r.jacobi_cd" => u.jacobi_cd(&m, p, rm, &mut cc),
            "r.jacobi_ns" => u.jacobi_ns(&m, p, rm, &mut cc),
            "r.jacobi_nc" => u.jacobi_nc(&m, p, rm, &mut cc),
            "r.jacobi_nd" => u.jacobi_nd(&m, p, rm, &mut cc),
            "r.jacobi_sc" => u.jacobi_sc(&m, p, rm, &mut cc),
            "r.jacobi_sd" => u.jacobi_sd(&m, p, rm, &mut cc),
            "r.jacobi_cs" => u.jacobi_cs(&m, p, rm, &mut cc),
            "r.jacobi_ds" => u.jacobi_ds(&m, p, rm, &mut cc),
            "r.jacobi_dc" => u.jacobi_dc(&m, p, rm, &mut cc),
            _ => panic!("function {f}"),
        };
        let k: Vec<&str> = t[4].split(':').collect();
        assert_eq!(k[0], "F", "{what}: expected kind {}", t[4]);
        let e: i64 = k[2].parse().expect("exp");
        let w =
            parse_real(&format!("{}{}p{}", k[1], k[3], e - 4 * k[3].len() as i64)).expect("want");
        assert_eq!(v.compare(&w), Some(0), "{what}: {v} want {}", t[4]);
    }
}

/// ZF-28: Ai, Bi, Ai', Bi' at large negative x. The oscillatory expansion formed the phase
/// χ = ζ + π/4, ζ = (2/3)|x|^{3/2}, at the working precision only, so for x = -2^1000
/// (ζ ~ 2^1500) the phase was off by far more than 2π; and near a zero of the combination
/// (sin χ P - cos χ Q, or c1 f - c2 g in the series) the cancelled bits were still certified.
/// ζ and π now carry exponent(ζ) extra bits, and the Ziv claim drops the cancelled bits.
/// Expected values are FLINT-certified (2p+64 bits; land39/zf28: the 156 ZF-28 rows and the 91
/// sample rows that changed, all wrong before).
#[test]
fn zf_28_airy_large_negative_argument() {
    let mut cc = Consts::new().expect("Consts");
    let mode = |s: &str| match s {
        "up" => RoundingMode::Up,
        "down" => RoundingMode::Down,
        "tozero" => RoundingMode::ToZero,
        "fromzero" => RoundingMode::FromZero,
        "toeven" => RoundingMode::ToEven,
        "toodd" => RoundingMode::ToOdd,
        _ => RoundingMode::None,
    };
    for l in include_str!("tier1/zf28.tsv").lines() {
        let t: Vec<&str> = l.split('\t').collect();
        let (f, p, rm) = (t[0], t[1].parse::<usize>().expect("p"), mode(t[2]));
        let x = parse_real(t[3]).expect("x");
        let what = format!("{f}({}) p={p} {rm:?}", t[3]);
        let v = match f {
            "r.ai" => x.ai(p, rm, &mut cc),
            "r.bi" => x.bi(p, rm, &mut cc),
            "r.ai_prime" => x.ai_prime(p, rm, &mut cc),
            "r.bi_prime" => x.bi_prime(p, rm, &mut cc),
            _ => panic!("function {f}"),
        };
        let k: Vec<&str> = t[4].split(':').collect();
        assert_eq!(k[0], "F", "{what}: expected kind {}", t[4]);
        let e: i64 = k[2].parse().expect("exp");
        let w =
            parse_real(&format!("{}{}p{}", k[1], k[3], e - 4 * k[3].len() as i64)).expect("want");
        assert_eq!(v.compare(&w), Some(0), "{what}: {v} want {}", t[4]);
    }
}

/// Dest-128 ToEven \(I_{100}(1000)=J_{100}(1000i)\). Gold is the FLINT 3.6
/// (python-flint 0.9.0) midpoint at prec 1024 (rel_accuracy_bits 1012),
/// cross-checked with mpmath at prec 1024, then rounded ToEven at 128 bits:
/// `F:+:1430:908ee542bcfc5f9d74534a3955e4f135`. The fractional part of the
/// high-prec mid is about 0.79, so ToEven rounds up. Analytically
/// \(J_{100}(1000i)=i^{100}I_{100}(1000)=I_{100}(1000)\) is real, so the
/// imaginary parts are +0. The 0F1 bsplit / prefix path
/// (`sum_bessel_i_series`, `hypgeom_0f1`) already produces this value; the
/// harness token `64p0` is 0x64 = 100, not decimal 64.
#[test]
fn zf_i100_z1000_toeven_128_matches_flint() {
    let mut cc = Consts::new().expect("Consts");
    let p = 128;
    let rm = RoundingMode::ToEven;
    let z = ExactComplex::new(
        parse_real("3e8p0").expect("1000"),
        parse_real("0").expect("0"),
    );
    let zi = ExactComplex::new(
        parse_real("0").expect("0"),
        parse_real("3e8p0").expect("1000i"),
    );
    let nu = ExactComplex::new(
        parse_real("64p0").expect("100"),
        parse_real("0").expect("0"),
    );
    let gold_re = "F:+:1430:908ee542bcfc5f9d74534a3955e4f135";
    let gold_im = "Z:+";
    let want = |s: &str| -> ExactNum {
        let k: Vec<&str> = s.split(':').collect();
        if k[0] == "Z" {
            let z = ExactNum::new(64);
            return if k[1] == "-" { std::ops::Neg::neg(z) } else { z };
        }
        let e: i64 = k[2].parse().expect("exp");
        parse_real(&format!("{}{}p{}", k[1], k[3], e - 4 * k[3].len() as i64)).expect("want")
    };
    let same = |v: &ExactNum, w: &ExactNum| {
        v.compare(w) == Some(0) && (!v.is_zero() || v.is_negative() == w.is_negative())
    };
    let i = z.bessel_i(&nu, p, rm, &mut cc);
    assert!(
        same(i.re(), &want(gold_re)),
        "I_100(1000) re {} want {gold_re}",
        i.re()
    );
    assert!(
        same(i.im(), &want(gold_im)),
        "I_100(1000) im {} want {gold_im}",
        i.im()
    );
    let j = zi.bessel_j_nu(&nu, p, rm, &mut cc);
    assert!(
        same(j.re(), &want(gold_re)),
        "J_100(1000i) re {} want {gold_re}",
        j.re()
    );
    assert!(
        same(j.im(), &want(gold_im)),
        "J_100(1000i) im {} want {gold_im}",
        j.im()
    );

    // J_5(1000i) = i I_5(1000): real +0, imag = I_5(1000).
    let nu5 = ExactComplex::new(parse_real("5p0").expect("5"), parse_real("0").expect("0"));
    let i5 = z.bessel_i(&nu5, p, rm, &mut cc);
    let j5 = zi.bessel_j_nu(&nu5, p, rm, &mut cc);
    assert!(
        same(i5.im(), &want(gold_im)),
        "I_5(1000) im {} want {gold_im}",
        i5.im()
    );
    assert!(
        same(j5.re(), &want(gold_im)),
        "J_5(1000i) re {} want {gold_im}",
        j5.re()
    );
    assert!(
        same(j5.im(), i5.re()),
        "J_5(1000i) im {} want I_5(1000) {}",
        j5.im(),
        i5.re()
    );

    // I_100(1000i) = i^{100} J_100(1000) = J_100(1000), imag +0.
    let i_on_i = zi.bessel_i(&nu, p, rm, &mut cc);
    let j_on_r = z.bessel_j_nu(&nu, p, rm, &mut cc);
    assert!(
        same(i_on_i.im(), &want(gold_im)),
        "I_100(1000i) im {} want {gold_im}",
        i_on_i.im()
    );
    assert!(
        same(j_on_r.im(), &want(gold_im)),
        "J_100(1000) im {} want {gold_im}",
        j_on_r.im()
    );
    assert!(
        same(i_on_i.re(), j_on_r.re()),
        "I_100(1000i) re {} want J_100(1000) {}",
        i_on_i.re(),
        j_on_r.re()
    );
}

/// ZF-58: Kelvin functions at large order (ν = 100, 498.5 at x = 1000) and the complex
/// \(J_\nu, I_\nu, K_\nu\) under them. The Hankel expansion diverges from its first term when
/// \(|\nu|^2\gtrsim|z|\), and `hankel_pq` used to stop at the first growing term and return
/// that truncation (errors up to \(2^{148}\) relative). \(K\) past the switch also went through
/// a \({}_2F_0\) Lentz step good to about 13 bits. Expected values are FLINT-certified
/// Expected values are FLINT 3.6 correctly-rounded (python-flint 0.9.0, dest-p
/// directed rounding; mpmath cross-check of complex midpoints). 7ca563a copied
/// implementation output (`land39/zf58/cx-new.tsv`, `out-new.tsv`) and accepted
/// it with a 2-ULP mpmath relative check (`land39/run/mpcx.py`); 17 of those
/// golds were 1–2 ULP stale or held an algebraically-zero leaked part.
#[test]
fn zf_58_kelvin_and_complex_bessel_large_order() {
    let mut cc = Consts::new().expect("Consts");
    let mode = |s: &str| match s {
        "up" => RoundingMode::Up,
        "down" => RoundingMode::Down,
        "tozero" => RoundingMode::ToZero,
        "fromzero" => RoundingMode::FromZero,
        "toeven" => RoundingMode::ToEven,
        "toodd" => RoundingMode::ToOdd,
        _ => RoundingMode::None,
    };
    let want = |s: &str| -> ExactNum {
        let k: Vec<&str> = s.split(':').collect();
        if k[0] == "Z" {
            let z = ExactNum::new(64);
            return if k[1] == "-" { std::ops::Neg::neg(z) } else { z };
        }
        assert_eq!(k[0], "F", "expected kind {s}");
        let e: i64 = k[2].parse().expect("exp");
        parse_real(&format!("{}{}p{}", k[1], k[3], e - 4 * k[3].len() as i64)).expect("want")
    };
    let same = |v: &ExactNum, w: &ExactNum| {
        v.compare(w) == Some(0) && (!v.is_zero() || v.is_negative() == w.is_negative())
    };
    let cpx = |s: &str| {
        let (a, b) = s.split_once(',').expect("re,im");
        ExactComplex::new(parse_real(a).expect("re"), parse_real(b).expect("im"))
    };
    for l in include_str!("tier1/zf58.tsv").lines() {
        let t: Vec<&str> = l.split('\t').collect();
        let (f, p, rm) = (t[0], t[1].parse::<usize>().expect("p"), mode(t[2]));
        let what = format!("{f}({}) p={p} {rm:?}", t[3..t.len() - 1].join(", "));
        let exp = t[t.len() - 1];
        if let Some(cf) = f.strip_prefix("c.") {
            let (z, nu) = (cpx(t[3]), cpx(t[4]));
            let v = match cf {
                "bessel_j_nu" => z.bessel_j_nu(&nu, p, rm, &mut cc),
                "bessel_i" => z.bessel_i(&nu, p, rm, &mut cc),
                "bessel_k" => z.bessel_k(&nu, p, rm, &mut cc),
                _ => panic!("function {f}"),
            };
            let (a, b) = exp.split_once(';').expect("re;im");
            assert!(same(v.re(), &want(a)), "{what}: re {} want {a}", v.re());
            assert!(same(v.im(), &want(b)), "{what}: im {} want {b}", v.im());
            continue;
        }
        let x = parse_real(t[3]).expect("x");
        let v = match f {
            "r.kelvin_ker" => x.kelvin_ker(p, rm, &mut cc),
            _ => {
                let nu = parse_real(t[4]).expect("nu");
                match f {
                    "r.kelvin_ber_nu" => x.kelvin_ber_nu(&nu, p, rm, &mut cc),
                    "r.kelvin_bei_nu" => x.kelvin_bei_nu(&nu, p, rm, &mut cc),
                    "r.kelvin_ker_nu" => x.kelvin_ker_nu(&nu, p, rm, &mut cc),
                    "r.kelvin_kei_nu" => x.kelvin_kei_nu(&nu, p, rm, &mut cc),
                    _ => panic!("function {f}"),
                }
            }
        };
        assert!(same(&v, &want(exp)), "{what}: {v} want {exp}");
    }
}

/// ZF-18 / ZF-68: complex functions whose parts were certified by two approximations agreeing
/// (`ziv_complex_agree`), which is not a bound: tiny parts beside large ones, cancellation in
/// ln/atan/acos formulas, the e^{-z²} phase in erf/erfc/Fresnel. Each part now carries a
/// rigorous error bound (midpoint-radius balls, analytic series-tail and DLMF 7.12(i)
/// remainder bounds) and is rounded only when the bound certifies it. Also: the C99 sign of
/// the imaginary part of acos on the cut, log(z, z) = 1 + 0i exactly, Fresnel zero parts on
/// the axes signed by the limit rule, and erf/erfc on the real axis through the real
/// functions (erfc(−76) took about 60 s). Expected values are FLINT-certified for nonzero
/// parts (2p+64 bits; land39/zf18, zf18r, zf68, zf18t); zero parts follow C99 (acos) or were
/// checked against mpmath at z ± 2^-300 (Fresnel); log(z, z) and ±1 − δ / 2 − δ rows are
/// exact by construction.
#[test]
fn zf_18_68_complex_parts_carry_error_bounds() {
    let mut cc = Consts::new().expect("Consts");
    let mut bad = Vec::new();
    for l in include_str!("tier1/zf18.tsv").lines() {
        let t: Vec<&str> = l.split('\t').collect();
        let (f, p, rm) = (
            t[0],
            t[1].parse::<usize>().expect("p"),
            mode(t[2]).expect("mode"),
        );
        let want = t[t.len() - 1];
        let z = parse_complex(t[3]).expect("z");
        let c = |i: usize| parse_complex(t[i]).expect("arg");
        let v = match f {
            "c.acos" => z.acos(p, rm, &mut cc),
            "c.asin" => z.asin(p, rm, &mut cc),
            "c.atan" => z.atan(p, rm, &mut cc),
            "c.acosh" => z.acosh(p, rm, &mut cc),
            "c.asinh" => z.asinh(p, rm, &mut cc),
            "c.atanh" => z.atanh(p, rm, &mut cc),
            "c.cosh" => z.cosh(p, rm, &mut cc),
            "c.log1p" => z.log1p(p, rm, &mut cc),
            "c.log" => z.log(&c(4), p, rm, &mut cc),
            "c.pow" => z.pow(&c(4), p, rm, &mut cc),
            "c.erf" => z.erf(p, rm, &mut cc),
            "c.erfc" => z.erfc(p, rm, &mut cc),
            "c.fresnel_s" => z.fresnel_s(p, rm, &mut cc),
            "c.fresnel_c" => z.fresnel_c(p, rm, &mut cc),
            "c.elliptic_k" => z.elliptic_k(p, rm, &mut cc),
            "c.elliptic_e_complete" => z.elliptic_e_complete(p, rm, &mut cc),
            "c.elliptic_pi_complete" => z.elliptic_pi_complete(&c(4), p, rm, &mut cc),
            _ => panic!("function {f}"),
        };
        let got = fmt_c(&v);
        // log(z, z) is exactly 1 + 0i: no inexact flag on either part.
        if f == "c.log" && (v.re().inexact() || v.im().inexact()) {
            bad.push(format!("{f}({}) p={p} {rm:?}: flagged inexact", t[3]));
        }
        if got != want {
            bad.push(format!(
                "{f}({}) p={p} {rm:?}\n    want {want}\n    got  {got}",
                t[3..t.len() - 1].join(", ")
            ));
        }
    }
    assert!(
        bad.is_empty(),
        "{} rows wrong:\n{}",
        bad.len(),
        bad.join("\n")
    );
}

/// ZF-36: J_{−n}(x) = (−1)^n J_n(x), I_n(−x) = (−1)^n I_n(x). A negation after rounding
/// must mirror Up/Down, not reuse the same direction. Expected values: FLINT >2p+64 bits.
#[test]
fn zf_36_bessel_reflection_directed_rounding() {
    let mut cc = Consts::new().expect("Consts");
    let mut bad = Vec::new();
    for l in include_str!("tier1/zf36-reflect.tsv").lines() {
        let t: Vec<&str> = l.split('\t').collect();
        let p = t[1].parse::<usize>().expect("p");
        let got = run(t[0], p, mode(t[2]).expect("rm"), &t[3..5], &mut cc).expect("dispatch");
        if got != t[5] {
            bad.push(format!(
                "{} p={p} {} args={:?} want {} got {got}",
                t[0],
                t[2],
                &t[3..5],
                t[5]
            ));
        }
    }
    assert!(
        bad.is_empty(),
        "{} wrong rows:\n{}",
        bad.len(),
        bad.join("\n")
    );
}

/// ZF-36: the J_ν power series near a zero of J_ν (J_{1/2}(π) at 256 bits) certified p_wrk
/// bits although the sum was ~2^64 below its largest term: 64 wrong bits returned as correctly
/// rounded. The claim now subtracts the cancellation (peak term vs sum) and the accumulated
/// per-term rounding. Expected value: FLINT (2p+64 bits).
#[test]
fn zf_36_bessel_series_cancellation_not_certified() {
    let mut cc = Consts::new().expect("Consts");
    for l in include_str!("tier1/zf36-cancel.tsv").lines() {
        let t: Vec<&str> = l.split('\t').collect();
        let p = t[1].parse::<usize>().expect("p");
        let got = run(t[0], p, mode(t[2]).expect("rm"), &t[3..5], &mut cc).expect("dispatch");
        assert_eq!(got, t[5], "{} p={p} {} args={:?}", t[0], t[2], &t[3..5]);
    }
}

/// ZF-33 / ZF-75 / ZF-36 (elliptic): amplitude-form F, E, Pi and inverse Jacobi
/// cn, dn. A rounded π/2 is not π/2; `1 - sin²` from a rounded sine drops the
/// cosine near π/2; claims ignored cancellation and the conditioning of a
/// reduced angle; arccn/arcdn formed `sqrt(1 - x²)` once at a fixed width.
/// Every expected value is FLINT-checked (2p+64 bits).
#[test]
fn zf_33_75_elliptic_amplitude_and_inverse_jacobi() {
    let mut cc = Consts::new().expect("Consts");
    let mut bad = Vec::new();
    for l in include_str!("tier1/zf33-ell.tsv").lines() {
        let t: Vec<&str> = l.split('\t').collect();
        let n = t.len();
        let p = t[1].parse::<usize>().expect("p");
        let got = run(t[0], p, mode(t[2]).expect("rm"), &t[3..n - 1], &mut cc).expect("dispatch");
        if got != t[n - 1] {
            bad.push(format!(
                "{} p={p} {} args={:?} got {got} want {}",
                t[0],
                t[2],
                &t[3..n - 1],
                t[n - 1]
            ));
        }
    }
    assert!(bad.is_empty(), "{} wrong:\n{}", bad.len(), bad.join("\n"));
}

/// CAPS: `hurwitz_em` capped its Bernoulli tail at 48 terms (`clamp(4, 48)`), so the
/// Euler–Maclaurin remainder stayed near 2^-500 and `hurwitz_zeta` and `polygamma`
/// (n ≥ 1) were silently wrong past about 500 bits (e.g. ψ'(1) at p = 1344 had 610
/// correct bits). The tail length now comes from the remainder bound.
/// References: FLINT acb zeta(s, a) / polygamma at 9400 bits.
#[test]
fn caps_hurwitz_zeta_and_polygamma_above_500_bits() {
    let mut cc = Consts::new().expect("Consts");
    let mut bad = Vec::new();
    for l in include_str!("tier1/caps-hurwitz.tsv").lines() {
        let t: Vec<&str> = l.split('\t').collect();
        let k = t[1].parse::<usize>().expect("order");
        let x = parse_real(t[2]).expect("x");
        for p in t[3].split(',').map(|p| p.parse::<usize>().expect("p")) {
            for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
                let v = match t[0] {
                    "pg" => x.polygamma(k, p, rm, &mut cc),
                    "hz" => ExactNum::from_u32(k as u32, 64).hurwitz_zeta(&x, p, rm, &mut cc),
                    _ => panic!("function"),
                };
                let w = cr_ref(t[4], p, rm, &mut cc);
                if v.compare(&w) != Some(0) {
                    bad.push(format!(
                        "{}({k}, {}) p={p} {rm:?}: {v} want {w}",
                        t[0], t[2]
                    ));
                }
            }
        }
    }
    assert!(bad.is_empty(), "{} wrong:\n{}", bad.len(), bad.join("\n"));
}

/// ZF-42: coincident poles of the Meijer G / Fox H residue sum went through a
/// hypercomb limit (three parameter shifts at most, two had to agree, 12,000-bit
/// cap), and at z = 0 every such row came back +0. Now: a non-removable pole
/// (a_k - b_h a positive integer, k <= n, h <= m; DLMF 16.17) is NaN; z = 0 is the
/// limit from the leftmost Mellin-Barnes pole (G^{2,0}_{1,2}(z|;1/2|0,0) ~ -ln z/sqrt(pi)
/// -> +inf); z != 0 is the constant term of the Laurent expansion in the parameter
/// shift (polygamma jets, rigorous pFq tail). See RC-3.
/// References: mpmath meijerg at 4p+512 bits, agreeing with 2p+256 to 2^-(2p+64);
/// z = 0 rows analytic. Every mode must be within 1 ulp.
#[test]
fn zf_42_meijer_fox_coincident_poles() {
    let mut cc = Consts::new().expect("Consts");
    let mut bad = Vec::new();
    for l in include_str!("tier1/zf42.tsv").lines() {
        let t: Vec<&str> = l.split('\t').collect();
        let n = t.len();
        let p = t[1].parse::<usize>().expect("p");
        let want = t[n - 1];
        for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
            let g = match meijer_fox(t[0], p, rm, &t[2..n - 1], &mut cc) {
                Ok(g) => g,
                Err(e) => panic!("dispatch {e}"),
            };
            let v = fmt_real(&g);
            let ok = match want {
                "nan" => v == "N:Some(InvalidArgument)",
                "+inf" => v == "I:+",
                "+0" => v == "Z:+",
                _ => {
                    let w = ExactNum::parse(want, Radix::Dec, 8960, RoundingMode::ToEven, &mut cc);
                    if g.is_nan() || g.is_inf() || g.is_zero() {
                        false
                    } else {
                        let d = ExactNum::sub_ref(&g, &w, p + 512, RoundingMode::None).abs();
                        let ulp = ExactNum::from_u8(1, 64).ldexp(
                            g.exponent().expect("exp") - p as i64,
                            64,
                            RoundingMode::None,
                        );
                        matches!(d.compare(&ulp), Some(c) if c <= 0)
                    }
                }
            };
            if !ok {
                bad.push(format!(
                    "{} p={p} {rm:?} args={:?} got {v} want {want}",
                    t[0],
                    &t[2..n - 1]
                ));
            }
        }
    }
    assert!(bad.is_empty(), "{} wrong:\n{}", bad.len(), bad.join("\n"));
}

/// ZF-43: the inside-disk Heun series (heun_g, heun_c) stopped after `series_cap`
/// = p + 160 terms, or when one term fell below 2^-p of the sum. At z = 3/4 with
/// radius 1 the terms shrink like (3/4)^j, so the sum kept about p/2 bits. Now the
/// three-term recurrence's tail is bounded by a majorant (characteristic roots 1 and
/// 1/a, or 1 and 0) and summing stops only below that bound; no term cap.
/// References: the same series in mpmath at 3p+256 and 4p+512 bits with N and 2N
/// terms, all four agreeing to 2^-(2p+64); cross-checked by mpmath odefun.
/// Every mode within 1 ulp.
#[test]
fn zf_43_heun_series_rigorous_tail() {
    let mut cc = Consts::new().expect("Consts");
    let mut bad = Vec::new();
    for l in include_str!("tier1/zf43.tsv").lines() {
        let t: Vec<&str> = l.split('\t').collect();
        let n = t.len();
        let p = t[1].parse::<usize>().expect("p");
        let x: Vec<ExactNum> = t[2..n - 1]
            .iter()
            .map(|s| parse_real(s).expect("arg"))
            .collect();
        let w = ExactNum::parse(t[n - 1], Radix::Dec, 8960, RoundingMode::ToEven, &mut cc);
        for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
            let v = match t[0] {
                "r.heun_g" => x[0].heun_g(&x[1], &x[2], &x[3], &x[4], &x[5], &x[6], p, rm, &mut cc),
                "r.heun_c" => x[0].heun_c(&x[1], &x[2], &x[3], &x[4], &x[5], p, rm, &mut cc),
                _ => panic!("function"),
            };
            let ok = !v.is_nan() && !v.is_inf() && !v.is_zero() && {
                let d = ExactNum::sub_ref(&v, &w, p + 512, RoundingMode::None).abs();
                let ulp = ExactNum::from_u8(1, 64).ldexp(
                    v.exponent().expect("exp") - p as i64,
                    64,
                    RoundingMode::None,
                );
                matches!(d.compare(&ulp), Some(c) if c <= 0)
            };
            if !ok {
                bad.push(format!(
                    "{} p={p} {rm:?} args={:?} got {}",
                    t[0],
                    &t[2..n - 1],
                    fmt_real(&v)
                ));
            }
        }
    }
    assert!(bad.is_empty(), "{} wrong:\n{}", bad.len(), bad.join("\n"));
}

/// ZF-63: integer-order J_n for n > 48 used Miller's backward recurrence started at
/// n + p/2 + 8 whatever x was. For x past that start (x > 4n + 2p + 256, e.g.
/// J_100(1000) at 64 bits) the recurrence had not settled and the normalized value
/// was wrong from about bit 60 on. Now the power series covers every n (its guard
/// bits cover the e^{|x|} cancellation). References: FLINT arb bessel_j at 9400 bits.
/// Every mode correctly rounded.
#[test]
fn zf_63_bessel_j_integer_order_large_x() {
    let mut cc = Consts::new().expect("Consts");
    let mut bad = Vec::new();
    for l in include_str!("tier1/zf63.tsv").lines() {
        let t: Vec<&str> = l.split('\t').collect();
        let x = parse_real(t[0]).expect("x");
        let nu = parse_real(t[1]).expect("nu");
        for p in t[2].split(',').map(|p| p.parse::<usize>().expect("p")) {
            for rm in MODES.iter().copied().filter(|m| *m != RoundingMode::None) {
                let v = x.bessel_j_nu(&nu, p, rm, &mut cc);
                let w = cr_ref(t[3], p, rm, &mut cc);
                if v.compare(&w) != Some(0) {
                    bad.push(format!("J_{}({}) p={p} {rm:?}: {v} want {w}", t[1], t[0]));
                }
            }
        }
    }
    assert!(bad.is_empty(), "{} wrong:\n{}", bad.len(), bad.join("\n"));
}

/// ZF-57: Scorer's large-|x| asymptotics had the wrong coefficients and added terms
/// past the minimum; at x=+-10, p>=1024 they returned numbers near 2^453 / 2^697.
/// At +-maxv, x^3 overflowed and the representable 1/(pi x) was lost. Now: Taylor
/// with a geometric bound and cancellation-aware Ziv, Watson with a remainder bound
/// stopped before the minimum, and 1/(pi x) directly for huge x. References: mpmath
/// scorergi/scorerhi at 6000 and 9000 bits (agree to 2^-5000 relative); maxv on its
/// subnormal grid (q = 2^p/(pi (1-2^-p))), overflow C99 7.12.1. All modes CR.
/// Follow-up: at |x| = minnorm / minpos, x^3 in the tail test wrapped and the loop
/// never stopped; those rows are Gi(0) / Hi(0) rounded (|x| far below one ulp).
#[test]
fn zf_57_scorer_bounded_series_and_watson() {
    let mut cc = Consts::new().expect("Consts");
    let mut bad = Vec::new();
    for l in include_str!("tier1/zf57.tsv").lines() {
        let t: Vec<&str> = l.split('\t').collect();
        let x = parse_real(t[1]).expect("x");
        for p in t[2].split(',').map(|p| p.parse::<usize>().expect("p")) {
            for rm in MODES.iter().copied().filter(|m| {
                if t[3] == "*" {
                    *m != RoundingMode::None
                } else {
                    *m == mode(t[3]).expect("mode")
                }
            }) {
                let v = match t[0] {
                    "r.scorer_gi" => x.scorer_gi(p, rm, &mut cc),
                    "r.scorer_hi" => x.scorer_hi(p, rm, &mut cc),
                    _ => panic!("function"),
                };
                let ok = if let Some(w) = t[4].strip_prefix('=') {
                    fmt_real(&v) == w
                } else {
                    let w = cr_ref(t[4], p, rm, &mut cc);
                    v.compare(&w) == Some(0)
                };
                if !ok {
                    bad.push(format!("{}({}) p={p} {rm:?}: {}", t[0], t[1], fmt_real(&v)));
                }
            }
        }
    }
    assert!(bad.is_empty(), "{} wrong:\n{}", bad.len(), bad.join("\n"));
}

/// ZF-37 (re-check): J_ν(0) for non-integer ν > 0 went through the series, whose sum is an
/// exact zero; the cancellation claim (0 bits for a zero sum, 19cff72) never certified it
/// and the Ziv loop ran to PrecisionRetryExhausted (8-39 s at p >= 2048). Real I_ν(0) and
/// K_ν(0) refused poles and I_{-n}(0) = 0. Now x = ±0 is decided before any series:
/// 1 at ν = 0; a signed zero for integer orders (J_n(-x) = (-1)^n J_n(x), J_{-n} = (-1)^n J_n,
/// I_{-n} = I_n); +0 for non-integer ν > 0; ±Inf for non-integer ν < 0 signed by Γ(ν+1), and
/// +Inf for K (C99 7.12.1 poles, as Y_ν(0)). Expected strings from the limits, pole and zero
/// signs cross-checked with mpmath at x = 2^-200 (land39/zf37z/zgen.py).
#[test]
fn zf_37_bessel_j_i_k_at_zero() {
    let mut cc = Consts::new().expect("Consts");
    let mut bad = Vec::new();
    for l in include_str!("tier1/zf37-zero.tsv").lines() {
        let t: Vec<&str> = l.split('\t').collect();
        let x = parse_real(t[1]).expect("x");
        let nu = parse_real(t[2]).expect("nu");
        let p = t[3].parse::<usize>().expect("p");
        let rm = mode(t[4]).expect("mode");
        let v = match t[0] {
            "r.bessel_j_nu" => x.bessel_j_nu(&nu, p, rm, &mut cc),
            "r.bessel_i" => x.bessel_i(&nu, p, rm, &mut cc),
            "r.bessel_k" => x.bessel_k(&nu, p, rm, &mut cc),
            _ => panic!("function"),
        };
        if fmt_real(&v) != t[5] {
            bad.push(format!(
                "{}({}, ν={}) p={p} {rm:?}: {} want {}",
                t[0],
                t[1],
                t[2],
                fmt_real(&v),
                t[5]
            ));
        }
    }
    assert!(bad.is_empty(), "{} wrong:\n{}", bad.len(), bad.join("\n"));
}

/// ZF-37 (tiny-x): J_ν and I_ν at minpos, minnorm, ±2^{-100000}, and other tiny x took the
/// power series. For |ν| ≥ 1 the leading term (x/2)^ν / Γ(ν+1) underflows or the Ziv loop
/// hits PrecisionRetryExhausted (I_2(minpos), J_{5/2}(minpos), I_n(2^{-100000})). I_{±1}(±minpos)
/// ToEven flushed to a signed zero instead of minpos; J_2(minpos) up flushed to +0; I_0(minpos)
/// up missed the ulp above 1; J_{-3/2}(minpos) overflowed with the wrong sign. Now the leading
/// term is placed as a significand times a power of two (`bessel_ji_scaled`) so those rows
/// finish in one rounding. Expected strings from FLINT (python-flint 0.9.0); mpmath agrees
/// (`land39/zf37t/tgen.py`). On 0bca621, before the scaled path, 8 I_{±1}(±minpos) ToEven rows
/// and 18 of a 36-row minpos sample were wrong (signed zero, exhaustion, wrong-sign Inf).
#[test]
fn zf_37_bessel_j_i_tiny_and_extreme_x() {
    let mut cc = Consts::new().expect("Consts");
    let mut bad = Vec::new();
    for l in include_str!("tier1/zf37-tiny.tsv").lines() {
        let t: Vec<&str> = l.split('\t').collect();
        let x = parse_real(t[1]).expect("x");
        let nu = parse_real(t[2]).expect("nu");
        let p = t[3].parse::<usize>().expect("p");
        let rm = mode(t[4]).expect("mode");
        let v = match t[0] {
            "r.bessel_j_nu" => x.bessel_j_nu(&nu, p, rm, &mut cc),
            "r.bessel_i" => x.bessel_i(&nu, p, rm, &mut cc),
            _ => panic!("function"),
        };
        if fmt_real(&v) != t[5] {
            bad.push(format!(
                "{}({}, ν={}) p={p} {rm:?}: {} want {}",
                t[0],
                t[1],
                t[2],
                fmt_real(&v),
                t[5]
            ));
        }
    }
    assert!(bad.is_empty(), "{} wrong:\n{}", bad.len(), bad.join("\n"));
}

/// ZF-37: boundary rows for the 4-word `bessel_ji_scaled` skip. `x = 2^{-30}`
/// is the point that was 1.36 s at p = 8192; `ν = -1 + 2^{-300}` is mapped
/// onto a pole of `Γ(ν+1)` by 256-bit rounding and must not leak that pole.
/// Golds from FLINT/Arb unique-rounding enclosures (`tgen_ji_edge.py`).
#[test]
fn zf_37_bessel_j_i_estimate_boundaries() {
    let mut cc = Consts::new().expect("Consts");
    let mut bad = Vec::new();
    for l in include_str!("tier1/zf37-ji-edge.tsv").lines() {
        let t: Vec<&str> = l.split('\t').collect();
        let x = parse_real(t[1]).expect("x");
        let nu = parse_real(t[2]).expect("nu");
        let p = t[3].parse::<usize>().expect("p");
        let rm = mode(t[4]).expect("mode");
        let v = match t[0] {
            "r.bessel_j_nu" => x.bessel_j_nu(&nu, p, rm, &mut cc),
            "r.bessel_i" => x.bessel_i(&nu, p, rm, &mut cc),
            _ => panic!("function"),
        };
        if fmt_real(&v) != t[5] {
            bad.push(format!(
                "{}({}, ν={}) p={p} {rm:?}: {} want {}",
                t[0],
                t[1],
                t[2],
                fmt_real(&v),
                t[5]
            ));
        }
    }
    assert!(bad.is_empty(), "{} wrong:\n{}", bad.len(), bad.join("\n"));
}

/// ZF-37 (tiny-x Y_ν, K_ν): Y_ν(minpos) for ν ≠ 0 overflowed as +Inf (the
/// dominant term is −(1/π) Γ(ν) (x/2)^{-ν} for ν > 0, or a signed multiple for
/// ν < 0; Y_{-(n+1/2)}(minpos) underflows with the sign of (−1)^n J_{n+1/2}).
/// K_ν(minpos) overflowed as +Inf in every mode, including down / tozero, which
/// must return +MAX (C99 7.12.1). K_n(2^{-100000}) hit PrecisionRetryExhausted
/// in the directed modes (the value is a tiny correction to an exact dyadic
/// leading term, which Ziv could not certify). Now `bessel_yk_scaled` places
/// Γ(|ν|) (x/2)^{-|ν|} as a significand times a power of two and rounds once.
/// Expected strings from the small-x leading term (FLINT when it forms a tight
/// ball; L(1+t_1) certified when t_2 is below 2^{-(p+200)}; `tgen_yk.py`).
/// On d3bc81f, Y_1 / Y_{5/2} / Y_{-3/2}(minpos) and K_1(minpos) down / K_n(2^{-100000})
/// directed were wrong.
#[test]
fn zf_37_bessel_y_k_tiny_and_extreme_x() {
    let mut cc = Consts::new().expect("Consts");
    let mut bad = Vec::new();
    for l in include_str!("tier1/zf37-yk-tiny.tsv").lines() {
        let t: Vec<&str> = l.split('\t').collect();
        let x = parse_real(t[1]).expect("x");
        let nu = parse_real(t[2]).expect("nu");
        let p = t[3].parse::<usize>().expect("p");
        let rm = mode(t[4]).expect("mode");
        let v = match t[0] {
            "r.bessel_y" => x.bessel_y(&nu, p, rm, &mut cc),
            "r.bessel_k" => x.bessel_k(&nu, p, rm, &mut cc),
            _ => panic!("function"),
        };
        if fmt_real(&v) != t[5] {
            bad.push(format!(
                "{}({}, ν={}) p={p} {rm:?}: {} want {}",
                t[0],
                t[1],
                t[2],
                fmt_real(&v),
                t[5]
            ));
        }
    }
    assert!(bad.is_empty(), "{} wrong:\n{}", bad.len(), bad.join("\n"));
}

/// ZF-56: Anger \(\mathbf{J}_\nu\), Weber \(\mathbf{E}_\nu\), and Struve
/// \(\mathbf{H}_\nu\). The 11.6.1 path for \(\mathbf{H}_0(100)\) had a relative
/// error of 2.7e-5 because the term recurrence divided by \((\nu-k-1/2)\) instead
/// of multiplying; at this argument that is exactly the factor-of-four error in
/// \(t_1\). \(A_\nu\) at large \(z\) uses the Laplace / Watson expansion
/// (DLMF 11.11.4) with a Cauchy remainder (panels otherwise). Expected values
/// are dual-precision mpmath that agree after rounding in every mode
/// (land39/zf56/gen.py and gen_hangs.py), except \(\mathbf{H}_{-3/2}(0)\) and
/// \(\mathbf{H}_{-5/2}(0)\), which are 0 by the closed form in RC-4 (mpmath
/// reports ±Inf). The table includes the 16 s56 hangs at \(x=2^{20}\), the
/// five InvalidArgument refusals, and the high-p extras.
#[test]
fn zf_56_anger_weber_struve_correctly_rounded() {
    let mut cc = Consts::new().expect("Consts");
    let mut bad = Vec::new();
    for l in include_str!("tier1/zf56.tsv").lines() {
        let t: Vec<&str> = l.split('\t').collect();
        let p = t[3].parse::<usize>().expect("p");
        // Debug smoke: p=64 ToEven on modest z, plus the H_0(100) regression
        // (r.struve_h 64p0 0). Anger J at z=100 is ~3.8s per mode in debug;
        // the full 462-row table did not finish in 15 min and official
        // build-debug at 6e2165e hit the 20m step timeout (2026-11-17).
        // build-release and verify.sh still run every row.
        if cfg!(debug_assertions) {
            let h0_100 = t[0] == "r.struve_h" && t[1] == "64p0" && t[2] == "0";
            let large_z = matches!(t[1], "1p20" | "1p10" | "7dp3" | "3e8p0" | "64p0" | "19p2");
            if p != 64 || t[4] != "toeven" || (large_z && !h0_100) {
                continue;
            }
        }
        let z = parse_real(t[1]).expect("z");
        let nu = parse_real(t[2]).expect("nu");
        let rm = mode(t[4]).expect("mode");
        let v = match t[0] {
            "r.struve_h" => z.struve_h(&nu, p, rm, &mut cc),
            "r.anger_j" => z.anger_j(&nu, p, rm, &mut cc),
            "r.weber_e" => z.weber_e(&nu, p, rm, &mut cc),
            _ => panic!("function"),
        };
        if fmt_real(&v) != t[5] {
            bad.push(format!(
                "{}(z={}, ν={}) p={p} {rm:?}: {} want {}",
                t[0],
                t[1],
                t[2],
                fmt_real(&v),
                t[5]
            ));
        }
    }
    assert!(bad.is_empty(), "{} wrong:\n{}", bad.len(), bad.join("\n"));
}

/// ZF-66: complex complete Π(n,m) = R_F(0,1-m,1)+(n/3)R_J(0,1-m,1,1-n).
/// The R_J duplication follows Carlson 1995 (22)–(32) (d_m form, series in
/// E_2..E_5). Complete Π at m=2+0i uses the lower side of the cut
/// (Im(m)→0−), like Arb, the same sheet as m=2−0i. O(ε) leakage is
/// flushed. Golds for m=2+0i are the 0− values.
#[test]
fn zf_66_complex_elliptic_pi_complete() {
    let mut cc = Consts::new().expect("Consts");
    let mut bad = Vec::new();
    for l in include_str!("tier1/zf66.tsv").lines() {
        if l.is_empty() {
            continue;
        }
        let t: Vec<&str> = l.split('\t').collect();
        let (f, p, rm) = (
            t[0],
            t[1].parse::<usize>().expect("p"),
            mode(t[2]).expect("mode"),
        );
        let want = t[t.len() - 1];
        let n = parse_complex(t[3]).expect("n");
        let m = parse_complex(t[4]).expect("m");
        assert_eq!(f, "c.elliptic_pi_complete");
        let got = fmt_c(&n.elliptic_pi_complete(&m, p, rm, &mut cc));
        if got != want {
            bad.push(format!(
                "{f}({}) p={p} {rm:?}\n    want {want}\n    got  {got}",
                t[3..t.len() - 1].join(", ")
            ));
        }
    }
    assert!(
        bad.is_empty(),
        "{} rows wrong:\n{}",
        bad.len(),
        bad.join("\n")
    );
}

/// ZF-78: `r.format` of `max_value` in a non-power-of-two radix (3 and 36) hung in
/// `conv_to_general`, which divided by the radix once per unit of the radix exponent
/// (~10^18 iterations at `EXPONENT_MAX`). The child-process budget makes a hang fail
/// this test (and verify.sh) instead of stalling it. The same scale is used on the way
/// back, so the spelling must parse to the same value.
#[test]
fn zf_78_format_max_value_general_radix() {
    const BUDGET: std::time::Duration = if cfg!(debug_assertions) {
        std::time::Duration::from_secs(300)
    } else {
        std::time::Duration::from_secs(20)
    };
    if std::env::var_os("ZF78_FORMAT_CHILD").is_some() {
        zf_78_format_max_value_general_radix_body();
        return;
    }
    let exe = std::env::current_exe().expect("current_exe");
    let mut child = std::process::Command::new(exe)
        .arg("--exact")
        .arg("zf_78_format_max_value_general_radix")
        .env("ZF78_FORMAT_CHILD", "1")
        .env("RUST_TEST_THREADS", "1")
        .spawn()
        .expect("spawn ZF-78 format child");
    let start = std::time::Instant::now();
    loop {
        match child.try_wait().expect("try_wait") {
            Some(status) => {
                assert!(
                    status.success(),
                    "ZF-78 r.format child exited {status}: format or parse of max_value in radix 3/36 failed"
                );
                return;
            }
            None if start.elapsed() > BUDGET => {
                let _ = child.kill();
                let _ = child.wait();
                panic!(
                    "ZF-78 r.format of max_value in radix 3/36 timed out after {}s (hung in conv_to_general or convert_from)",
                    BUDGET.as_secs()
                );
            }
            None => std::thread::sleep(std::time::Duration::from_millis(20)),
        }
    }
}

fn zf_78_format_max_value_general_radix_body() {
    let mut cc = Consts::new().expect("Consts");
    let cases: [(usize, RoundingMode, u8); 4] = [
        (64, RoundingMode::ToEven, 3),
        (64, RoundingMode::ToEven, 36),
        (128, RoundingMode::ToZero, 3),
        (128, RoundingMode::ToEven, 36),
    ];
    for (p, rm, base) in cases {
        let rdx = Radix::try_new(base).expect("radix");
        let v = ExactNum::max_value(p);
        let s = v
            .format(rdx, rm, &mut cc)
            .unwrap_or_else(|e| panic!("format max_value p={p} radix={base} {rm:?}: {e:?}"));
        assert!(
            !s.is_empty() && s.contains('e'),
            "format max_value p={p} radix={base}: {s}"
        );
        let back = ExactNum::parse(&s, rdx, p, rm, &mut cc);
        assert!(
            !back.is_nan() && back.compare(&v) == Some(0) && back.sign() == v.sign(),
            "round-trip max_value p={p} radix={base} {rm:?}: {s} -> {back}"
        );
        let n = ExactNum::min_value(p);
        let sn = n
            .format(rdx, rm, &mut cc)
            .unwrap_or_else(|e| panic!("format min_value p={p} radix={base} {rm:?}: {e:?}"));
        let backn = ExactNum::parse(&sn, rdx, p, rm, &mut cc);
        assert!(
            !backn.is_nan() && backn.compare(&n) == Some(0) && backn.sign() == n.sign(),
            "round-trip min_value p={p} radix={base} {rm:?}: {sn} -> {backn}"
        );
    }
    // Modest values still format and parse (the new scale must not break the old path).
    let rdx = Radix::try_new(3).expect("radix 3");
    let three = ExactNum::from_u8(3, 64);
    let s = three
        .format(rdx, RoundingMode::ToEven, &mut cc)
        .expect("format 3");
    let back = ExactNum::parse(&s, rdx, 64, RoundingMode::ToEven, &mut cc);
    assert_eq!(
        back.compare(&three),
        Some(0),
        "format/parse 3 in radix 3: {s} -> {back}"
    );
}

/// ZF-78 Bessel/Kelvin at p=2048: Kelvin is Re/Im of complex J/K, and those
/// call complex gamma. A child process enforces the budget so a hang fails
/// verify instead of stalling it.
#[test]
fn zf_78_bessel_kelvin_high_prec() {
    const BUDGET: std::time::Duration = if cfg!(debug_assertions) {
        std::time::Duration::from_secs(300)
    } else {
        std::time::Duration::from_secs(20)
    };
    if std::env::var_os("ZF78_BESSEL_CHILD").is_some() {
        zf_78_bessel_kelvin_high_prec_body();
        return;
    }
    let exe = std::env::current_exe().expect("current_exe");
    let mut child = std::process::Command::new(exe)
        .arg("--exact")
        .arg("zf_78_bessel_kelvin_high_prec")
        .env("ZF78_BESSEL_CHILD", "1")
        .env("RUST_TEST_THREADS", "1")
        .spawn()
        .expect("spawn ZF-78 Bessel/Kelvin child");
    let start = std::time::Instant::now();
    loop {
        match child.try_wait().expect("try_wait") {
            Some(status) => {
                assert!(
                    status.success(),
                    "ZF-78 Bessel/Kelvin child exited {status}: a p=2048 case returned a NaN or failed"
                );
                return;
            }
            None if start.elapsed() > BUDGET => {
                let _ = child.kill();
                let _ = child.wait();
                panic!(
                    "ZF-78 Bessel/Kelvin or complex gamma at p=2048 timed out after {}s",
                    BUDGET.as_secs()
                );
            }
            None => std::thread::sleep(std::time::Duration::from_millis(20)),
        }
    }
}

fn zf_78_bessel_kelvin_high_prec_body() {
    let mut cc = Consts::new().expect("Consts");
    let p = 2048usize;
    let rm = RoundingMode::ToEven;
    let check = |name: &str, v: ExactNum| {
        assert!(!v.is_nan() && !v.is_inf(), "ZF-78 {name} at p=2048: {v:?}");
    };
    let checkc = |name: &str, v: ExactComplex| {
        assert!(
            !v.is_nan() && !v.re().is_inf() && !v.im().is_inf(),
            "ZF-78 {name} at p=2048: {v:?}"
        );
    };

    let z = ExactComplex::new(ExactNum::from_i64(3, p), ExactNum::from_i64(1, p));
    checkc("c.gamma(3+i)", z.gamma(p, rm, &mut cc));
    checkc("c.digamma(3+i)", z.digamma(p, rm, &mut cc));

    let z = ExactComplex::new(ExactNum::from_i64(1, p), ExactNum::from_i64(1, p));
    let nu = ExactComplex::new(ExactNum::from_i64(2, p), ExactNum::from_i64(1, p));
    checkc(
        "c.bessel_j_nu(1+i, 2+i)",
        z.bessel_j_nu(&nu, p, rm, &mut cc),
    );

    let x = ExactNum::from_i64(100, p);
    let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, RoundingMode::None);
    check(
        "r.kelvin_ber_nu(100, 1/2)",
        x.kelvin_ber_nu(&half, p, rm, &mut cc),
    );

    let x = ExactNum::from_i64(700, p);
    check("r.kelvin_ker(700)", x.kelvin_ker(p, rm, &mut cc));
    check("r.kelvin_kei(700)", x.kelvin_kei(p, rm, &mut cc));

    let x = ExactNum::from_i64(1000, p);
    let nu = ExactNum::from_i64(100, p);
    check("r.bessel_y(1000, 100)", x.bessel_y(&nu, p, rm, &mut cc));

    let nu = ExactComplex::new(ExactNum::from_i64(2, p), ExactNum::from_i64(1, p));
    for below in [false, true] {
        let z = zf_78_neg_real(-2, p, below);
        let tag = if below { "-0" } else { "+0" };
        checkc(
            &format!("c.bessel_i(-2{tag}i, 2+i)"),
            z.bessel_i(&nu, p, rm, &mut cc),
        );
        checkc(
            &format!("c.bessel_j_nu(-2{tag}i, 2+i)"),
            z.bessel_j_nu(&nu, p, rm, &mut cc),
        );
        checkc(
            &format!("c.bessel_k(-2{tag}i, 2+i)"),
            z.bessel_k(&nu, p, rm, &mut cc),
        );
        if below {
            checkc(
                &format!("c.bessel_y(-2{tag}i, 2+i)"),
                z.bessel_y(&nu, p, rm, &mut cc),
            );
        }
    }

    let x = ExactNum::from_u8(1, p).ldexp(-30, p, RoundingMode::None);
    let n = ExactNum::from_i64(100, p);
    check(
        "r.kelvin_kei_nu(2^-30, 100)",
        x.kelvin_kei_nu(&n, p, rm, &mut cc),
    );

    // x < 0 is outside the real ker/kei domain (DLMF 10.61(iii) / DOM-07).
    let x = ExactNum::from_i64(-700, p);
    for (name, v) in [
        ("r.kelvin_ker(-700)", x.kelvin_ker(p, rm, &mut cc)),
        ("r.kelvin_kei(-700)", x.kelvin_kei(p, rm, &mut cc)),
    ] {
        assert!(
            v.is_nan() && v.err() == Some(Error::InvalidArgument),
            "ZF-78 {name} at p={p}: {v:?} {:?}",
            v.err()
        );
    }
}

/// FLINT gold for a Bessel value on the cut: mid of a 1024-bit Arb ball for the
/// DLMF continuation (or the principal value, which agrees with it from above),
/// parsed at dest p and `rm`.
fn zf_78_assert_flint_gold(
    name: &str,
    below: bool,
    got: &ExactComplex,
    re: &str,
    im: &str,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) {
    let want_re = ExactNum::parse(re, Radix::Dec, p, rm, cc);
    let want_im = ExactNum::parse(im, Radix::Dec, p, rm, cc);
    let tag = if below { "-0" } else { "+0" };
    assert_eq!(
        got.re().compare(&want_re),
        Some(0),
        "ZF-78 {name}_{{2+i}}(-2{tag}i) re got {:?} want {:?}",
        got.re(),
        want_re
    );
    assert_eq!(
        got.im().compare(&want_im),
        Some(0),
        "ZF-78 {name}_{{2+i}}(-2{tag}i) im got {:?} want {:?}",
        got.im(),
        want_im
    );
}

/// `-re + 0i` or `-re - 0i`.
fn zf_78_neg_real(re: i64, p: usize, below: bool) -> ExactComplex {
    let mut im = ExactNum::new(p);
    if below {
        im.set_sign(Sign::Neg);
    }
    ExactComplex::new(ExactNum::from_i64(re, p), im)
}

/// I_{2+i}(-2±0i) at p=64 was NaN(PrecisionRetryExhausted): complex pow of an axis
/// base with a non-real exponent (the J series evaluates (z/2)^ν = (-1)^{2+i} or
/// i^{2+i}) could not certify an exact-axis result. kei_100(2^{-30}) was
/// NaN(InvalidArgument) because BESSEL_INTEGER_MAX=64 refused a computable order.
#[test]
fn zf_78_bessel_cut_and_large_integer_order() {
    let mut cc = Consts::new().expect("Consts");
    let p = 64usize;
    let rm = RoundingMode::ToEven;
    let nu = ExactComplex::new(ExactNum::from_i64(2, p), ExactNum::from_i64(1, p));

    let m1 = ExactComplex::from_real(-ExactNum::from_u8(1, p), p);
    let pr = m1.pow(&nu, p, rm, &mut cc);
    assert!(
        !pr.is_nan() && pr.im().is_zero(),
        "(-1)^{{2+i}} must be finite and exactly real, got {pr:?}"
    );
    let i = ExactComplex::i(p);
    let pr = i.pow(&nu, p, rm, &mut cc);
    assert!(
        !pr.is_nan() && pr.im().is_zero(),
        "i^{{2+i}} must be finite and exactly real, got {pr:?}"
    );

    for below in [false, true] {
        let z = zf_78_neg_real(-2, p, below);
        for (name, v) in [
            ("bessel_i", z.bessel_i(&nu, p, rm, &mut cc)),
            ("bessel_j_nu", z.bessel_j_nu(&nu, p, rm, &mut cc)),
            ("bessel_k", z.bessel_k(&nu, p, rm, &mut cc)),
        ] {
            assert!(
                !v.is_nan() && !v.re().is_inf() && !v.im().is_inf(),
                "ZF-78 {name}(-2{}i, 2+i) at p=64: {v:?}",
                if below { "-0" } else { "+0" }
            );
        }
        if below {
            let v = z.bessel_y(&nu, p, rm, &mut cc);
            assert!(
                !v.is_nan() && !v.re().is_inf() && !v.im().is_inf(),
                "ZF-78 bessel_y(-2-0i, 2+i) at p=64: {v:?}"
            );
        }
    }

    let x = ExactNum::from_u8(1, p).ldexp(-30, p, RoundingMode::None);
    let n = ExactNum::from_i64(100, p);
    let v = x.kelvin_kei_nu(&n, p, rm, &mut cc);
    assert!(
        !v.is_nan() && !v.is_inf(),
        "ZF-78 kelvin_kei_nu(2^-30, 100) at p=64: {v:?}"
    );

    // Jeff 2026-11-18 10:37 CT: both signs of a zero imag take Arb's
    // principal value (Arg = +π, the from-above / FLINT-at-the-real-point
    // sheet). Both must be finite, they must agree, and each part must
    // match the principal FLINT gold at p=64, ToEven.
    let z_up = zf_78_neg_real(-2, p, false);
    let z_dn = zf_78_neg_real(-2, p, true);
    let cases: [(&str, ExactComplex, ExactComplex, [&str; 4]); 4] = [
        (
            "J",
            z_up.bessel_j_nu(&nu, p, rm, &mut cc),
            z_dn.bessel_j_nu(&nu, p, rm, &mut cc),
            [
                "0.01289969292761315243851446812712015859682302002749544162683870",
                "-0.01430184864595774337140776114452022311739960572981015521554012",
                "6.907677921568666159196563207501494876492990721647332611063157",
                "-7.658520608488526894222929735809437183625240906142614899680370",
            ],
        ),
        (
            "Y",
            z_up.bessel_y(&nu, p, rm, &mut cc),
            z_dn.bessel_y(&nu, p, rm, &mut cc),
            [
                "-7.749291772708986615068442909890537724203549729903309683446119",
                "0.06337691413685817313114306027523968820442896394198164109011553",
                "-7.701622372260576175622524266260842865391746735056883430167623",
                "-6.933383044104028443121702084924531502095861823516815494963442",
            ],
        ),
        (
            "I",
            z_up.bessel_i(&nu, p, rm, &mut cc),
            z_dn.bessel_i(&nu, p, rm, &mut cc),
            [
                "0.01780419567714751260293720585820168400046382551450928538466223",
                "-0.03025478149847086275393487467602418024589508504774901062635675",
                "9.533998218442581256047246609482398876981936003884589403597782",
                "-16.20118303215618470638483804891152037276873886144499266690323",
            ],
        ),
        (
            "K",
            z_up.bessel_k(&nu, p, rm, &mut cc),
            z_dn.bessel_k(&nu, p, rm, &mut cc),
            [
                "1.300877352591998256865287720150978803539513694735036658853111",
                "2.206017884177394529560202667673656632668947708363570896773199",
                "2.206017884177394529560202667673656632668947708363570896773199",
                "1.300877352591998256865287720150978803539513694735036658853111",
            ],
        ),
    ];
    for (name, up, dn, golds) in cases {
        assert!(
            !up.is_nan() && !dn.is_nan(),
            "ZF-78 {name}_{{2+i}}(-2±0i) must be finite, got {up:?} / {dn:?}"
        );
        assert!(
            up.re().compare(dn.re()) == Some(0) && up.im().compare(dn.im()) == Some(0),
            "ZF-78 {name}_{{2+i}}(-2+0i) and (-2-0i) must agree, got {up:?} / {dn:?}"
        );
        zf_78_assert_flint_gold(name, false, &up, golds[0], golds[1], p, rm, &mut cc);
        zf_78_assert_flint_gold(name, true, &dn, golds[0], golds[1], p, rm, &mut cc);
    }
}

/// BESSEL_K_RECUR_MAX=1<<22 refused a computable K_{2^{22}+1}(1) (P9).
/// The forward recurrence scales with ν; this case is ~3.5s at p=64.
#[test]
fn zf_78_k_recur_scales_with_nu() {
    const BUDGET: std::time::Duration = if cfg!(debug_assertions) {
        std::time::Duration::from_secs(300)
    } else {
        std::time::Duration::from_secs(20)
    };
    if std::env::var_os("ZF78_K_RECUR_CHILD").is_some() {
        zf_78_k_recur_scales_with_nu_body();
        return;
    }
    let exe = std::env::current_exe().expect("current_exe");
    let mut child = std::process::Command::new(exe)
        .arg("--exact")
        .arg("zf_78_k_recur_scales_with_nu")
        .env("ZF78_K_RECUR_CHILD", "1")
        .env("RUST_TEST_THREADS", "1")
        .spawn()
        .expect("spawn ZF-78 K-recur child");
    let start = std::time::Instant::now();
    loop {
        match child.try_wait().expect("try_wait") {
            Some(status) => {
                assert!(
                    status.success(),
                    "ZF-78 K_{{2^22+1}}(1) child exited {status}: NaN or failed"
                );
                return;
            }
            None if start.elapsed() > BUDGET => {
                let _ = child.kill();
                let _ = child.wait();
                panic!(
                    "ZF-78 K_{{2^22+1}}(1) at p=64 timed out after {}s",
                    BUDGET.as_secs()
                );
            }
            None => std::thread::sleep(std::time::Duration::from_millis(20)),
        }
    }
}

fn zf_78_k_recur_scales_with_nu_body() {
    let mut cc = Consts::new().expect("Consts");
    let p = 64usize;
    let z = ExactComplex::from_real(ExactNum::from_u8(1, p), p);
    let n = (1u32 << 22) + 1;
    let nu = ExactComplex::from_real(ExactNum::from_u32(n, p), p);
    let v = z.bessel_k(&nu, p, RoundingMode::ToEven, &mut cc);
    assert!(
        !v.is_nan() && !v.re().is_inf() && !v.im().is_inf(),
        "ZF-78 K_{{2^22+1}}(1) at p=64: {v:?}"
    );
}

/// The six unique ZF-78 \( {}_2F_1 \) rows, plus a series point whose term count
/// at p=2048 exceeded the old `HYPERGEOM_SERIES_MAX_TERMS = 10_000` cap (P9).
/// p=64 ToEven matches the FLINT mid. At a real \( z>1 \), both signs of a zero
/// imaginary part take the lower side, like Arb (the below-cut principal value).
/// p=2048 must finish in 20s, be finite, and round to the p=64 value.
#[test]
fn zf_78_hypergeom_2f1_cluster() {
    const BUDGET: std::time::Duration = if cfg!(debug_assertions) {
        std::time::Duration::from_secs(300)
    } else {
        std::time::Duration::from_secs(20)
    };
    if std::env::var_os("ZF78_2F1_CHILD").is_some() {
        zf_78_hypergeom_2f1_cluster_body();
        return;
    }
    let exe = std::env::current_exe().expect("current_exe");
    let mut child = std::process::Command::new(exe)
        .arg("--exact")
        .arg("zf_78_hypergeom_2f1_cluster")
        .env("ZF78_2F1_CHILD", "1")
        .env("RUST_TEST_THREADS", "1")
        .spawn()
        .expect("spawn ZF-78 2F1 child");
    let start = std::time::Instant::now();
    loop {
        match child.try_wait().expect("try_wait") {
            Some(status) => {
                assert!(
                    status.success(),
                    "ZF-78 2F1 child exited {status}: a case returned a NaN or missed the FLINT gold"
                );
                return;
            }
            None if start.elapsed() > BUDGET => {
                let _ = child.kill();
                let _ = child.wait();
                panic!(
                    "ZF-78 2F1 cluster at p=2048 timed out after {}s",
                    BUDGET.as_secs()
                );
            }
            None => std::thread::sleep(std::time::Duration::from_millis(20)),
        }
    }
}

fn zf_78_assert_2f1_gold(
    name: &str,
    got: &ExactComplex,
    re: &str,
    im: &str,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) {
    let want_re = ExactNum::parse(re, Radix::Dec, p, rm, cc);
    let want_im = ExactNum::parse(im, Radix::Dec, p, rm, cc);
    assert_eq!(
        got.re().compare(&want_re),
        Some(0),
        "ZF-78 2F1 {name} re got {:?} want {:?}",
        got.re(),
        want_re
    );
    assert_eq!(
        got.im().compare(&want_im),
        Some(0),
        "ZF-78 2F1 {name} im got {:?} want {:?}",
        got.im(),
        want_im
    );
}

fn zf_78_hypergeom_2f1_cluster_body() {
    let mut cc = Consts::new().expect("Consts");
    let wide = include_str!("tier1/zf78-2f1-wide-z.txt").trim();
    let rm = RoundingMode::ToEven;
    // FLINT mid at 256 bits. At z = 2+0i and z = 2−0i both take the lower side.
    let cases: [(&str, &str, &str, &str, &str, &str, &str); 8] = [
        (
            "c0761826 F(1/2,1/2;2+i;-1)",
            "1p-1,0",
            "1p-1,0",
            "1p1,1p0",
            "-1p0,0",
            "0.92072262708275837786712360805728628602799069219411548973238180029657339613",
            "0.03235611082650910346536440749915671160351782362927500643352641950839688426",
        ),
        (
            "c0761861 F(1/2,1/2;2+i;-10+3i)",
            "1p-1,0",
            "1p-1,0",
            "1p1,1p0",
            "-5p1,3p0",
            "0.66681839428499363368129640978596526311968096669941884583066293640544022011",
            "0.11984026591530631213159843424749347544401185206401276434360092219006987221",
        ),
        (
            "c0761843 F(1/2,1/2;2+i;1)",
            "1p-1,0",
            "1p-1,0",
            "1p1,1p0",
            "1p0,0",
            "1.127140861772219510746400137255872644156689401450758418015484299124369638697",
            "-0.138385675471678180466555538286756565464739926118396650212084059792167987619",
        ),
        (
            "c0761854 F(1/2,1/2;2+i;z~1)",
            "1p-1,0",
            "1p-1,0",
            "1p1,1p0",
            wide,
            "1.127140911414630255855879582461655609150710095565160575502083986395663267282",
            "-0.138385396885359943499287365770071517115573303013655895387859957146745433783",
        ),
        (
            "c0762058 F(1/2,1;2+i;1)",
            "1p-1,0",
            "1p0,0",
            "1p1,1p0",
            "1p0,0",
            "1.2",
            "-0.4",
        ),
        (
            "c0761974 F(1/2,1;3/2;2+0)",
            "1p-1,0",
            "1p0,0",
            "3p-1,0",
            "1p1,0",
            "0.623225240140230513394020309942232395602743196240974504536684781385770235240",
            "-1.110720734539591561753969911840061311970824862924995003680889242810802882592",
        ),
        (
            "c0761974 F(1/2,1;3/2;2-0)",
            "1p-1,0",
            "1p0,0",
            "3p-1,0",
            "1p1,-0",
            "0.623225240140230513394020309942232395602743196240974504536684781385770235240",
            "-1.110720734539591561753969911840061311970824862924995003680889242810802882592",
        ),
        (
            "P9 F(1/2,1/2;1;1/2+3/4i)",
            "1p-1,0",
            "1p-1,0",
            "1p0,0",
            "1p-1,3p-2",
            "1.01134360641644933657737508920502311133426764512561677104958460340215240502",
            "0.25014127863591117960917551424098165809387102641562108318368164308968223889",
        ),
    ];

    let p64 = 64usize;
    let mut v64: Vec<ExactComplex> = Vec::with_capacity(cases.len());
    for (name, a, b, c, z, re, im) in cases {
        let a = parse_complex(a).expect("a");
        let b = parse_complex(b).expect("b");
        let c = parse_complex(c).expect("c");
        let z = parse_complex(z).expect("z");
        let v = a.hypergeom_2f1(&b, &c, &z, p64, rm, &mut cc);
        assert!(
            !v.is_nan() && !v.re().is_inf() && !v.im().is_inf(),
            "ZF-78 2F1 {name} at p=64: {v:?}"
        );
        zf_78_assert_2f1_gold(name, &v, re, im, p64, rm, &mut cc);
        v64.push(v);
    }
    assert!(
        v64[5].re().compare(v64[6].re()) == Some(0) && v64[5].im().compare(v64[6].im()) == Some(0),
        "ZF-78 2F1 F(1/2,1;3/2;2±0i) must match on the lower side, got {:?} / {:?}",
        v64[5],
        v64[6]
    );

    let p = 2048usize;
    for ((name, a, b, c, z, _, _), want64) in cases.iter().zip(v64.iter()) {
        let a = parse_complex(a).expect("a");
        let b = parse_complex(b).expect("b");
        let c = parse_complex(c).expect("c");
        let z = parse_complex(z).expect("z");
        let v = a.hypergeom_2f1(&b, &c, &z, p, rm, &mut cc);
        assert!(
            !v.is_nan() && !v.re().is_inf() && !v.im().is_inf(),
            "ZF-78 2F1 {name} at p=2048: {v:?}"
        );
        let mut rounded = v.clone();
        rounded
            .set_precision(p64, rm)
            .expect("set_precision 2048 -> 64");
        assert_eq!(
            rounded.re().compare(want64.re()),
            Some(0),
            "ZF-78 2F1 {name} p=2048 rounded to 64 re got {:?} want {:?}",
            rounded.re(),
            want64.re()
        );
        assert_eq!(
            rounded.im().compare(want64.im()),
            Some(0),
            "ZF-78 2F1 {name} p=2048 rounded to 64 im got {:?} want {:?}",
            rounded.im(),
            want64.im()
        );
    }
}

/// Incomplete F, E, Π and complete K, E inherit the Carlson cut (−∞,0] of
/// 1−x², 1−m x², 1−n x², or 1−m. IEEE subtraction collapses the ±0 imag of
/// those arguments to +0, so E(−2−0i | 2) (c0880554) used to return the +0
/// (above-cut) value. C99: the sign of a zero imag selects the side. Each
/// side must match the one-sided FLINT limit (ε=2^{−80}) at p=64 ToEven.
/// FLINT at the exact real argument is NaN; the stock judge expected for
/// c0880512 / c0880554 is not a one-sided limit (RC-8). There is no public
/// elliptic_d; D = (F−E)/m or (x³/3) R_D and is covered by the same RF/RD
/// arguments.
#[test]
fn zf_78_elliptic_cut_signed_zeros() {
    const BUDGET: std::time::Duration = if cfg!(debug_assertions) {
        std::time::Duration::from_secs(300)
    } else {
        std::time::Duration::from_secs(20)
    };
    if std::env::var_os("ZF78_ELLIPTIC_CHILD").is_some() {
        zf_78_elliptic_cut_signed_zeros_body();
        return;
    }
    let exe = std::env::current_exe().expect("current_exe");
    let mut child = std::process::Command::new(exe)
        .arg("--exact")
        .arg("zf_78_elliptic_cut_signed_zeros")
        .env("ZF78_ELLIPTIC_CHILD", "1")
        .env("RUST_TEST_THREADS", "1")
        .spawn()
        .expect("spawn ZF-78 elliptic-cut child");
    let start = std::time::Instant::now();
    loop {
        match child.try_wait().expect("try_wait") {
            Some(status) => {
                assert!(
                    status.success(),
                    "ZF-78 elliptic cut child exited {status}: a case returned a NaN or missed the FLINT gold"
                );
                return;
            }
            None if start.elapsed() > BUDGET => {
                let _ = child.kill();
                let _ = child.wait();
                panic!(
                    "ZF-78 elliptic cut signed-zero cases at p=64 timed out after {}s",
                    BUDGET.as_secs()
                );
            }
            None => std::thread::sleep(std::time::Duration::from_millis(20)),
        }
    }
}

fn zf_78_assert_elliptic_gold(
    name: &str,
    got: &ExactComplex,
    re: &str,
    im: &str,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) {
    let want_re = ExactNum::parse(re, Radix::Dec, p, rm, cc);
    let want_im = ExactNum::parse(im, Radix::Dec, p, rm, cc);
    assert_eq!(
        got.re().compare(&want_re),
        Some(0),
        "ZF-78 elliptic {name} re got {:?} want {:?}",
        got.re(),
        want_re
    );
    assert_eq!(
        got.im().compare(&want_im),
        Some(0),
        "ZF-78 elliptic {name} im got {:?} want {:?}",
        got.im(),
        want_im
    );
}

fn zf_78_neg_dec(s: &str) -> String {
    if let Some(rest) = s.strip_prefix('-') {
        rest.to_string()
    } else {
        format!("-{s}")
    }
}

fn zf_78_elliptic_cut_signed_zeros_body() {
    let mut cc = Consts::new().expect("Consts");
    let p = 64usize;
    let rm = RoundingMode::ToEven;
    let m2 = parse_complex("1p1,0").expect("m=2");
    let n_half = parse_complex("1p-1,0").expect("n=1/2");
    let m_half = parse_complex("1p-1,0").expect("m=1/2");

    // (name, up, dn, re, im_up). The below-cut imag is −im_up in every case.
    let cases: [(&str, ExactComplex, ExactComplex, &str, &str); 7] = [
        (
            "E(-2±0i | 2) c0880512/c0880554",
            parse_complex("-1p1,0")
                .expect("x")
                .elliptic_e(&m2, p, rm, &mut cc),
            parse_complex("-1p1,-0")
                .expect("x")
                .elliptic_e(&m2, p, rm, &mut cc),
            "-2.63669247722547559233620709360741588042513526117868563359698996500171402407",
            "-0.5990701173677961037199599826009051310508626225268478439922876726089812428",
        ),
        (
            "F(-2±0i | 2)",
            parse_complex("-1p1,0")
                .expect("x")
                .elliptic_f(&m2, p, rm, &mut cc),
            parse_complex("-1p1,-0")
                .expect("x")
                .elliptic_f(&m2, p, rm, &mut cc),
            "-0.37874246652437093562046128447701577951284771373622504753823190852372026605",
            "1.31102877714605990523241961443995159140384266013281473353133230849122081297",
        ),
        (
            "Pi(1/2; -2±0i | 2)",
            n_half.elliptic_pi(&parse_complex("-1p1,0").expect("x"), &m2, p, rm, &mut cc),
            n_half.elliptic_pi(&parse_complex("-1p1,-0").expect("x"), &m2, p, rm, &mut cc),
            "0.0995915243930391820610697880792957204629788145473009342595312260467085250",
            "0.8215835656063735014928681931220844436377086709437080027147165182805491517",
        ),
        (
            "E(2±0i | 2)",
            parse_complex("1p1,0")
                .expect("x")
                .elliptic_e(&m2, p, rm, &mut cc),
            parse_complex("1p1,-0")
                .expect("x")
                .elliptic_e(&m2, p, rm, &mut cc),
            "2.63669247722547559233620709360741588042513526117868563359698996500171402407",
            "-0.5990701173677961037199599826009051310508626225268478439922876726089812428",
        ),
        (
            "F(3/2±0i | 1/2)",
            parse_complex("3p-1,0")
                .expect("x")
                .elliptic_f(&m_half, p, rm, &mut cc),
            parse_complex("3p-1,-0")
                .expect("x")
                .elliptic_f(&m_half, p, rm, &mut cc),
            "1.382085160391020806804427985499581068510236276318718233307052585290649564",
            "1.854074677301371918433848254575442445549518350151002472692264583523955947",
        ),
        (
            "K(2±0i)",
            m2.elliptic_k(p, rm, &mut cc),
            parse_complex("1p1,-0")
                .expect("m")
                .elliptic_k(p, rm, &mut cc),
            "1.311028777146059905232419942175159783795596939280507209689025452359131603269",
            "1.311028777146059905232419399946366306633275882577287003212506857059021323791",
        ),
        (
            "E_complete(2±0i)",
            m2.elliptic_e_complete(p, rm, &mut cc),
            parse_complex("1p1,-0")
                .expect("m")
                .elliptic_e_complete(p, rm, &mut cc),
            "0.59907011736779610371996164113935533932170792474625901822465827710326058245",
            "-0.59907011736779610371996139336976201606782579517276645916851812891986794387",
        ),
    ];
    for (name, up, dn, re, im_up) in cases {
        assert!(
            !up.is_nan() && !dn.is_nan() && !up.re().is_inf() && !dn.re().is_inf(),
            "ZF-78 elliptic {name} must be finite, got {up:?} / {dn:?}"
        );
        assert!(
            up.re().compare(dn.re()) != Some(0) || up.im().compare(dn.im()) != Some(0),
            "ZF-78 elliptic {name} +0 and −0 must differ, both {up:?}"
        );
        zf_78_assert_elliptic_gold(name, &up, re, im_up, p, rm, &mut cc);
        let im_dn = zf_78_neg_dec(im_up);
        zf_78_assert_elliptic_gold(name, &dn, re, &im_dn, p, rm, &mut cc);
    }

    // The on-cut RF/RD path must stay near the off-cut cost at dest p=8192
    // (the survey cap is 16s). Extra work of 2p+128, needed only for R_J,
    // made E(-2±0i | 2) time out. Check finite, sides differ, and the
    // p=8192 value rounds to the p=64 FLINT gold. Debug is too slow for
    // this dest (the 20s child budget is for release / verify.sh).
    if cfg!(debug_assertions) {
        return;
    }
    let p8 = 8192usize;
    let m2 = parse_complex("1p1,0").expect("m=2");
    let up8 = parse_complex("-1p1,0")
        .expect("x")
        .elliptic_e(&m2, p8, rm, &mut cc);
    let dn8 = parse_complex("-1p1,-0")
        .expect("x")
        .elliptic_e(&m2, p8, rm, &mut cc);
    assert!(
        !up8.is_nan() && !dn8.is_nan() && !up8.re().is_inf() && !dn8.re().is_inf(),
        "ZF-78 elliptic E(-2±0i | 2) at p=8192 must be finite, got {up8:?} / {dn8:?}"
    );
    assert!(
        up8.re().compare(dn8.re()) != Some(0) || up8.im().compare(dn8.im()) != Some(0),
        "ZF-78 elliptic E(-2+0i | 2) and E(-2-0i | 2) at p=8192 must differ, both {up8:?}"
    );
    for (name, v, re, im) in [
        (
            "E(-2+0i | 2) p=8192",
            up8,
            "-2.63669247722547559233620709360741588042513526117868563359698996500171402407",
            "-0.5990701173677961037199599826009051310508626225268478439922876726089812428",
        ),
        (
            "E(-2-0i | 2) p=8192",
            dn8,
            "-2.63669247722547559233620709360741588042513526117868563359698996500171402407",
            "0.5990701173677961037199599826009051310508626225268478439922876726089812428",
        ),
    ] {
        let mut rounded = v.clone();
        rounded
            .set_precision(p, rm)
            .expect("set_precision 8192 -> 64");
        zf_78_assert_elliptic_gold(name, &rounded, re, im, p, rm, &mut cc);
    }
}

/// Ai(x) for huge positive x. DLMF 9.7.15: 0 < Ai(x) < exp(-(2/3) x^{3/2}).
/// For x > 2^{64} that upper bound underflows the exponent range, so the correctly
/// rounded value is +0 in every mode that rounds toward zero (Down, ToZero, ToEven,
/// ToOdd, None) and min_positive(p) in Up / FromZero. The large-positive path used
/// to evaluate exp(+ξ) for Bi in the same call; that overflow was mapped to +Inf
/// and poisoned Ai.
#[test]
fn zf_airy_large_positive_underflow() {
    let mut cc = Consts::new().expect("Consts");
    for p in [4096usize, 8192] {
        let xs = [parse_real("1p100").expect("2^100"), ExactNum::max_value(p)];
        for x in xs {
            for rm in MODES {
                let v = x.ai(p, rm, &mut cc);
                let away = matches!(rm, RoundingMode::Up | RoundingMode::FromZero);
                if away {
                    let want = ExactNum::min_positive(p);
                    assert_eq!(
                        v.compare(&want),
                        Some(0),
                        "Ai(x) p={p} {rm:?}: {} want min_positive",
                        fmt_real(&v)
                    );
                } else {
                    assert!(
                        v.is_zero() && !v.is_negative(),
                        "Ai(x) p={p} {rm:?}: {} want +0",
                        fmt_real(&v)
                    );
                }
            }
        }
    }
}

/// Jeff 2026-11-17 12:59 PM CT: arccn(x|1) for x < 0 is NaN(InvalidArgument).
/// cn(u|1) = sech(u) never takes a non-positive value (DLMF 22.5.4). The guard
/// landed in 6bc4bcb; this pins x in {-1, -3/4} and is not a library change.
#[test]
fn zf_arccn_m1_negative_is_invalid() {
    let mut cc = Consts::new().expect("Consts");
    let m = parse_real("1p0").expect("m=1");
    for xs in ["-1p0", "-3p-2"] {
        let x = parse_real(xs).expect("x");
        for p in [64usize, 128, 1024, 4096] {
            for rm in MODES {
                let v = x.jacobi_arccn(&m, p, rm, &mut cc);
                assert!(
                    v.is_nan() && v.err() == Some(Error::InvalidArgument),
                    "arccn({xs}|1) p={p} {rm:?}: {v} {:?}",
                    v.err()
                );
            }
        }
    }
}

/// Real Kelvin domain at 0, x<0, and +Inf, from the recommended inventory in
/// `q-rulings.md` C-kelvin-domain (P20 / DOM-07 / DLMF 10.61; Jeff 2026-11-17).
/// `kei(0) = −π/4`, `ker(0) = +Inf`, `ker_1(0) = kei_1(0) = −∞`, integer
/// `ber_n(−x) = (−1)^n ber_n(x)` (and the same for `bei`); otherwise a negative
/// argument is `NaN(InvalidArgument)`. `ker_ν(+Inf) = kei_ν(+Inf) = +0`; `ber`/`bei`
/// at `+Inf` have no limit.
#[test]
fn zf_kelvin_domain_inventory() {
    let mut cc = Consts::new().expect("Consts");
    let invalid = |v: &ExactNum, what: &str| {
        assert!(
            v.is_nan() && v.err() == Some(Error::InvalidArgument),
            "{what}: {v} {:?}",
            v.err()
        );
    };
    let plus_zero = |v: &ExactNum, what: &str| {
        assert!(v.is_zero() && !v.is_negative(), "{what}: {v}");
    };
    let plus_inf = |v: &ExactNum, what: &str| {
        assert!(v.is_inf() && v.is_positive(), "{what}: {v}");
    };
    let minus_inf = |v: &ExactNum, what: &str| {
        assert!(v.is_inf() && v.is_negative(), "{what}: {v}");
    };
    let zeros = [signed_zero(false), signed_zero(true)];
    for p in [64usize, 128, 1024] {
        for rm in MODES {
            let one = ExactNum::from_u8(1, p);
            let half = parse_real("1p-1").expect("1/2");
            let nu1 = ExactNum::from_u8(1, p);
            let nu_m3 = parse_real("-3p0").expect("-3");
            let nu_mhalf = parse_real("-1p-1").expect("-1/2");
            // −π/4 correctly rounded: not −(π rounded)/4, which is 1 ulp off in directed modes.
            let wp = p + 128;
            let mut minus_pi4 = -cc
                .pi(wp, RoundingMode::None)
                .ldexp(-2, wp, RoundingMode::None);
            minus_pi4.set_precision(p, rm).expect("round -pi/4");

            for (zx, ztag) in [(&zeros[0], "+0"), (&zeros[1], "-0")] {
                let ber = zx.kelvin_ber(p, rm, &mut cc);
                assert!(
                    ber.compare(&one) == Some(0),
                    "ber({ztag}) p={p} {rm:?}: {ber}"
                );
                plus_zero(
                    &zx.kelvin_bei(p, rm, &mut cc),
                    &format!("bei({ztag}) p={p} {rm:?}"),
                );
                plus_inf(
                    &zx.kelvin_ker(p, rm, &mut cc),
                    &format!("ker({ztag}) p={p} {rm:?}"),
                );
                let kei = zx.kelvin_kei(p, rm, &mut cc);
                assert!(
                    kei.compare(&minus_pi4) == Some(0),
                    "kei({ztag}) p={p} {rm:?}: {kei} want {minus_pi4}"
                );
                minus_inf(
                    &zx.kelvin_ker_nu(&nu1, p, rm, &mut cc),
                    &format!("ker_1({ztag}) p={p} {rm:?}"),
                );
                minus_inf(
                    &zx.kelvin_kei_nu(&nu1, p, rm, &mut cc),
                    &format!("kei_1({ztag}) p={p} {rm:?}"),
                );
                plus_zero(
                    &zx.kelvin_ber_nu(&nu_m3, p, rm, &mut cc),
                    &format!("ber_{{-3}}({ztag}) p={p} {rm:?}"),
                );
                plus_zero(
                    &zx.kelvin_bei_nu(&nu_m3, p, rm, &mut cc),
                    &format!("bei_{{-3}}({ztag}) p={p} {rm:?}"),
                );
                plus_inf(
                    &zx.kelvin_ber_nu(&nu_mhalf, p, rm, &mut cc),
                    &format!("ber_{{-1/2}}({ztag}) p={p} {rm:?}"),
                );
                minus_inf(
                    &zx.kelvin_bei_nu(&nu_mhalf, p, rm, &mut cc),
                    &format!("bei_{{-1/2}}({ztag}) p={p} {rm:?}"),
                );
                plus_zero(
                    &zx.kelvin_ber_nu(&half, p, rm, &mut cc),
                    &format!("ber_{{1/2}}({ztag}) p={p} {rm:?}"),
                );
            }

            let xpos = ExactNum::from_u8(1, p);
            let xneg = -xpos.clone();
            invalid(
                &xneg.kelvin_ker(p, rm, &mut cc),
                &format!("ker(-1) p={p} {rm:?}"),
            );
            invalid(
                &xneg.kelvin_kei(p, rm, &mut cc),
                &format!("kei(-1) p={p} {rm:?}"),
            );
            invalid(
                &xneg.kelvin_ber_nu(&half, p, rm, &mut cc),
                &format!("ber_{{1/2}}(-1) p={p} {rm:?}"),
            );
            invalid(
                &xneg.kelvin_kei_nu(&nu_m3, p, rm, &mut cc),
                &format!("kei_{{-3}}(-1) p={p} {rm:?}"),
            );

            // DLMF 10.61.5: integer-order ber/bei at x<0.
            let ber_pos = xpos.kelvin_ber(p, rm, &mut cc);
            let ber_neg = xneg.kelvin_ber(p, rm, &mut cc);
            assert!(
                ber_pos.compare(&ber_neg) == Some(0),
                "ber(-1)=ber(1) p={p} {rm:?}: {ber_neg} vs {ber_pos}"
            );
            let ber1_pos = xpos.kelvin_ber_nu(&nu1, p, rm, &mut cc);
            let ber1_neg = xneg.kelvin_ber_nu(&nu1, p, rm, &mut cc);
            assert!(
                (-ber1_pos.clone()).compare(&ber1_neg) == Some(0),
                "ber_1(-1)=-ber_1(1) p={p} {rm:?}: {ber1_neg} vs -{ber1_pos}"
            );

            let inf = INF_POS;
            plus_zero(
                &inf.kelvin_ker(p, rm, &mut cc),
                &format!("ker(+Inf) p={p} {rm:?}"),
            );
            plus_zero(
                &inf.kelvin_kei(p, rm, &mut cc),
                &format!("kei(+Inf) p={p} {rm:?}"),
            );
            plus_zero(
                &inf.kelvin_ker_nu(&nu1, p, rm, &mut cc),
                &format!("ker_1(+Inf) p={p} {rm:?}"),
            );
            invalid(
                &inf.kelvin_ber(p, rm, &mut cc),
                &format!("ber(+Inf) p={p} {rm:?}"),
            );
            invalid(
                &inf.kelvin_bei(p, rm, &mut cc),
                &format!("bei(+Inf) p={p} {rm:?}"),
            );
            invalid(
                &inf.kelvin_ber_nu(&half, p, rm, &mut cc),
                &format!("ber_{{1/2}}(+Inf) p={p} {rm:?}"),
            );
        }
    }
}

/// In-domain Kelvin at x = 2^100 is past the exponent range (DLMF 10.67 leading
/// term). ber/bei overflow with the sign of the leading trigonometric factor;
/// ker/kei underflow with that sign. Today's single 10.61 evaluation returns
/// NaN(InvalidArgument) for ber/bei and +0 for ker/kei.
#[test]
fn zf_kelvin_in_domain_large_x_exp_range() {
    let mut cc = Consts::new().expect("Consts");
    let x = parse_real("1p100").expect("2^100");
    let half = parse_real("1p-1").expect("1/2");
    let underflow_neg = |v: &ExactNum, p: usize, rm: RoundingMode, what: &str| {
        let away = matches!(rm, RoundingMode::Down | RoundingMode::FromZero);
        if away {
            let want = -ExactNum::min_positive(p);
            assert!(
                v.compare(&want) == Some(0) && v.is_negative(),
                "{what} p={p} {rm:?}: {v} want -min_positive"
            );
        } else {
            assert!(
                v.is_zero() && v.is_negative(),
                "{what} p={p} {rm:?}: {v} want -0"
            );
        }
    };
    for p in [64usize, 128, 4096, 8192] {
        for rm in MODES {
            let ber = x.kelvin_ber(p, rm, &mut cc);
            assert!(
                ber.is_inf() && ber.is_negative(),
                "ber(2^100) p={p} {rm:?}: {ber}"
            );
            let bei = x.kelvin_bei(p, rm, &mut cc);
            assert!(
                bei.is_inf() && bei.is_positive(),
                "bei(2^100) p={p} {rm:?}: {bei}"
            );
            underflow_neg(&x.kelvin_ker(p, rm, &mut cc), p, rm, "ker(2^100)");
            underflow_neg(&x.kelvin_kei(p, rm, &mut cc), p, rm, "kei(2^100)");
            let ber_h = x.kelvin_ber_nu(&half, p, rm, &mut cc);
            assert!(
                ber_h.is_inf() && ber_h.is_negative(),
                "ber_{{1/2}}(2^100) p={p} {rm:?}: {ber_h}"
            );
            let bei_h = x.kelvin_bei_nu(&half, p, rm, &mut cc);
            assert!(
                bei_h.is_inf() && bei_h.is_positive(),
                "bei_{{1/2}}(2^100) p={p} {rm:?}: {bei_h}"
            );
            underflow_neg(
                &x.kelvin_ker_nu(&half, p, rm, &mut cc),
                p,
                rm,
                "ker_{1/2}(2^100)",
            );
            let kei_h = x.kelvin_kei_nu(&half, p, rm, &mut cc);
            let away = matches!(rm, RoundingMode::Up | RoundingMode::FromZero);
            if away {
                let want = ExactNum::min_positive(p);
                assert!(
                    kei_h.compare(&want) == Some(0) && !kei_h.is_negative(),
                    "kei_{{1/2}}(2^100) p={p} {rm:?}: {kei_h} want +min_positive"
                );
            } else {
                assert!(
                    kei_h.is_zero() && !kei_h.is_negative(),
                    "kei_{{1/2}}(2^100) p={p} {rm:?}: {kei_h} want +0"
                );
            }
        }
    }
}

/// The 10.31.2 / 10.61 small-z leading term for ker(x): −ln(x/2) − γ.
/// The remainder is O(x² ln x), which for x = min_positive is far below any
/// destination ulp at a magnitude of ~2^62.
fn ker0_small_z_leading(x: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> ExactNum {
    let wp = p + 256;
    let half = ExactNum::from_u8(1, wp).ldexp(-1, wp, RoundingMode::None);
    let arg = x.clone().mul(&half, wp, RoundingMode::None);
    let ln = arg.ln(wp, RoundingMode::None, cc);
    let g = cc.euler_gamma(wp, RoundingMode::None);
    let mut v = (-ln).sub(&g, wp, RoundingMode::None);
    let _ = v.set_precision(p, rm);
    v
}

/// maxv ker/kei: |e^{−x/√2}| underflows the exponent range (DLMF 10.67).
/// The result is a signed zero, never InvalidArgument, and never a min_positive
/// whose sign was not certified (the argument of the leading trig factor at
/// x = max_value needs ~2^62 bits of 1/π, which is not feasible).
#[test]
fn zf_kelvin_ker_kei_maxv_underflow() {
    let mut cc = Consts::new().expect("Consts");
    for p in [64usize, 2048, 4096] {
        let x = ExactNum::max_value(p);
        for rm in MODES {
            for (name, v) in
                [("ker", x.kelvin_ker(p, rm, &mut cc)), ("kei", x.kelvin_kei(p, rm, &mut cc))]
            {
                assert!(
                    !v.is_nan(),
                    "{name}(maxv) p={p} {rm:?}: {v} {:?} want a signed underflow, not NaN",
                    v.err()
                );
                let nearest = matches!(
                    rm,
                    RoundingMode::ToZero
                        | RoundingMode::ToEven
                        | RoundingMode::ToOdd
                        | RoundingMode::None
                );
                if nearest {
                    assert!(v.is_zero(), "{name}(maxv) p={p} {rm:?}: {v} want 0");
                } else {
                    // Sign is uncertified: do not invent a directed min_positive.
                    assert!(
                        v.is_zero(),
                        "{name}(maxv) p={p} {rm:?}: {v} want 0 (sign of the 10.67 trig factor is uncertified)"
                    );
                }
            }
        }
    }
}

/// ker(min_positive) is correctly rounded from the 10.31.2 leading term
/// (leftover: only ~253 bits agreed at p=2048 ToZero and at p=4096 every mode).
#[test]
fn zf_kelvin_ker_minpos_cr() {
    let mut cc = Consts::new().expect("Consts");
    let cases: &[(usize, RoundingMode)] = &[
        (2048, RoundingMode::ToZero),
        (4096, RoundingMode::ToEven),
        (4096, RoundingMode::Up),
        (4096, RoundingMode::Down),
    ];
    for &(p, rm) in cases {
        let x = ExactNum::min_positive(p);
        let got = x.kelvin_ker(p, rm, &mut cc);
        let want = ker0_small_z_leading(&x, p, rm, &mut cc);
        assert!(
            !got.is_nan() && got.compare(&want) == Some(0),
            "ker(minpos) p={p} {rm:?}: got {} want {}",
            fmt_real(&got),
            fmt_real(&want)
        );
    }
}

/// kei(x) for x in [2^62, 2^63): e^{−x/√2} already underflows, even though
/// x/√2 has not passed 2^63. The 10.67 path must give a signed underflow,
/// not an unsigned flush to +0 (leftover: six p=2048 modes of a 2048-bit x).
#[test]
fn zf_kelvin_kei_just_below_old_1067_threshold() {
    let mut cc = Consts::new().expect("Consts");
    let x = parse_real("1p62").expect("2^62");
    let p = 64usize;
    // Sign of −sin(x/√2 + π/8), certified at a modest working precision.
    let wrk = 512usize;
    let sqrt2 = ExactNum::from_u8(2, wrk).sqrt(wrk, RoundingMode::None);
    let xi = ExactNum::div_ref(&x, &sqrt2, wrk, RoundingMode::None);
    let eighth = cc
        .pi(wrk, RoundingMode::None)
        .ldexp(-3, wrk, RoundingMode::None);
    let t = -xi
        .add(&eighth, wrk, RoundingMode::None)
        .sin(wrk, RoundingMode::None, &mut cc);
    assert!(
        !t.is_nan() && !t.is_zero(),
        "sign of kei 10.67 factor at 2^62 must be certifiable, got {t}"
    );
    let neg = t.is_negative();
    for rm in MODES {
        let v = x.kelvin_kei(p, rm, &mut cc);
        let away = if neg {
            matches!(rm, RoundingMode::Down | RoundingMode::FromZero)
        } else {
            matches!(rm, RoundingMode::Up | RoundingMode::FromZero)
        };
        if away {
            let want = if neg { -ExactNum::min_positive(p) } else { ExactNum::min_positive(p) };
            assert!(
                v.compare(&want) == Some(0) && v.is_negative() == neg,
                "kei(2^62) p={p} {rm:?}: {v} want signed min_positive (neg={neg})"
            );
        } else {
            assert!(
                v.is_zero() && v.is_negative() == neg,
                "kei(2^62) p={p} {rm:?}: {v} want signed 0 (neg={neg})"
            );
        }
    }
}

/// a050bc5 stopped Carlson \(R_F\)/\(R_D\)/\(R_J\) duplication when
/// \((\mathrm{dmax}/\max(|A|,1))^{\mathrm{order}}<2^{-p}\). FLINT
/// `acb_elliptic_rf` uses \(\mathrm{rel}=\mathrm{dmax}/|A|\). When \(|A|<1\)
/// the floor of 1 stops too early. At \(m=1+2^{-60}i\), dest \(p=8192\),
/// that left about 492 trailing bits wrong (first error near bit 7700).
/// Golds are FLINT 3.6 midpoints (python-flint 0.9.0, `ctx.prec=8448`,
/// `rel_accuracy_bits` 8442 for \(K\) and 8435 for \(E\), radius
/// \(\sim 10^{-2541}\)), uniquely correctly rounded to dest \(p\).
#[test]
fn zf_elliptic_k_e_complete_near_one_p8192_matches_flint() {
    let mut cc = Consts::new().expect("Consts");
    let m = parse_complex("1p0,1p-60").expect("m");
    let p = 8192usize;
    let want = |s: &str| -> ExactNum {
        let k: Vec<&str> = s.split(':').collect();
        let e: i64 = k[2].parse().expect("exp");
        parse_real(&format!("{}{}p{}", k[1], k[3], e - 4 * k[3].len() as i64)).expect("want")
    };
    let same = |v: &ExactNum, w: &ExactNum| v.compare(w) == Some(0);
    for line in include_str!("tier1/zf_ellip_k_e_m1_p8192.tsv").lines() {
        if line.is_empty() || line.starts_with('#') || line.starts_with("func") {
            continue;
        }
        let t: Vec<&str> = line.split('\t').collect();
        let (func, rm) = (t[0], mode(t[1]).expect("mode"));
        let wr = want(t[2]);
        let wi = want(t[3]);
        let v = match func {
            "c.elliptic_k" => m.elliptic_k(p, rm, &mut cc),
            "c.elliptic_e_complete" => m.elliptic_e_complete(p, rm, &mut cc),
            other => panic!("unexpected {other}"),
        };
        assert!(
            same(v.re(), &wr),
            "{func} {rm:?} re differs from the FLINT dest-p gold"
        );
        assert!(
            same(v.im(), &wi),
            "{func} {rm:?} im differs from the FLINT dest-p gold"
        );
    }
}

/// B5 (Spence PR #40 review, Jeff ruling 2026-11-22 8:33 PM CT): J_ν / I_ν at the four signed zeros for real non-integer
/// ν < 0 return the signed one-sided limit. z = +0±0i (arg ±0): sgn Γ(ν+1)·∞ (DLMF 10.7.3 / 10.30.1), so −∞ on
/// (−2,−1), (−4,−3), …; the imaginary zero is −sgn Γ(ν+1)·sgn(Im z)·0 (Jeff ruling 2026-11-22 ~11:58 PM CT: match
/// FLINT, whose Im along z = r(1 ± i·r²) has that sign: arg (z/2)^ν = ±ν·r², ν < 0). z = −0±0i (arg ±π, C99 Annex G): sgn Γ(ν+1)·e^{±νπi}·∞
/// (DLMF 10.11.1 / 10.34.1): Re −sgn cos(π(ν−⌊ν⌋))·∞ (+0 at half-odd ν), Im ∓∞.
/// Expected signs: FLINT 3.6.0 `acb_bessel_j` / `acb_bessel_i` (python-flint 0.9.0) along each direction at |z| = 2^-200
/// (2^-3000 for ν = −1 ± 2^-2048), mpmath 1.3.0 for ν = −(2^20+1/2) where FLINT gives [±inf], cross-checked against
/// FLINT's sgn Γ(ν+1) (`land39/final/pr40-fixes/b5/b5_flint_oracle.{py,out}`). mpmath returns +inf at exact 0 for every
/// ν (land39 REF-CORRECTIONS RC-14). Every row fails on e6f439f where ν ∈ (−2,−1) etc. gave +∞ and −0±0i gave +∞±0i.
#[test]
fn b5_j_i_negative_order_at_signed_zeros_is_signed_limit() -> Result<(), String> {
    let mut cc = Consts::new().map_err(|e| format!("Consts: {e:?}"))?;
    let f2048 = format!("-{}p-2048", "f".repeat(512));
    let g2048 = format!("-1{}1p-2048", "0".repeat(511));
    // (ν, sign of Re at +0±0i, Re at −0±0i: '+' +∞, '-' −∞, '0' +0)
    let rows: [(&str, char, char); 14] = [
        ("-1p-2", '+', '+'),
        ("-1p-1", '+', '0'),
        ("-ffffffffffffffffp-64", '+', '-'),
        ("-10000000000000001p-64", '-', '+'),
        (f2048.as_str(), '+', '-'),
        (g2048.as_str(), '-', '+'),
        ("-5p-2", '-', '+'),
        ("-3p-1", '-', '0'),
        ("-5p-1", '+', '0'),
        ("-7p-1", '-', '0'),
        ("-200001p-1", '+', '0'),
        ("-3p-2", '+', '-'),
        ("-7p-2", '-', '-'),
        ("-bp-2", '+', '-'),
    ];
    let zeros = [(false, false), (false, true), (true, false), (true, true)];
    for (nu_t, plus_re, minus_re) in rows {
        let nu = ExactComplex::new(parse_real(nu_t)?, signed_zero(false));
        for rm in MODES {
            for p in [64usize, 128, 256, 1024, 2048] {
                for (re_neg, im_neg) in zeros {
                    let z = ExactComplex::new(signed_zero(re_neg), signed_zero(im_neg));
                    for (name, v) in [
                        ("J", z.bessel_j_nu(&nu, p, rm, &mut cc)),
                        ("I", z.bessel_i(&nu, p, rm, &mut cc)),
                    ] {
                        let what = format!(
                            "{name}_{nu_t}({}0{}0i) p={p} {rm:?}: got {v:?}",
                            if re_neg { "-" } else { "+" },
                            if im_neg { "-" } else { "+" }
                        );
                        if !re_neg {
                            let re_ok = match plus_re {
                                '+' => v.re().is_inf_pos(),
                                _ => v.re().is_inf_neg(),
                            };
                            assert!(re_ok, "{what}, want {plus_re}Inf re");
                            let want_im_neg = (plus_re == '+') != im_neg;
                            assert!(
                                v.im().is_zero() && v.im().is_negative() == want_im_neg,
                                "{what}, want imaginary zero -sgnGamma(nu+1)*sgn(Im z) (FLINT direction)"
                            );
                        } else {
                            let re_ok = match minus_re {
                                '+' => v.re().is_inf_pos(),
                                '-' => v.re().is_inf_neg(),
                                _ => v.re().is_zero() && !v.re().is_negative(),
                            };
                            assert!(re_ok, "{what}, want re {minus_re}");
                            let im_ok =
                                if im_neg { v.im().is_inf_pos() } else { v.im().is_inf_neg() };
                            assert!(
                                im_ok,
                                "{what}, want im {}Inf",
                                if im_neg { "+" } else { "-" }
                            );
                        }
                    }
                }
            }
        }
    }
    Ok(())
}
