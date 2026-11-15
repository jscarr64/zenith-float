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
    Consts, ExactComplex, ExactInt, ExactNum, ExactRational, Ieee32, Ieee64, Radix, RawKind,
    RawParts, RoundingMode, Sign, INF_NEG, INF_POS, NAN,
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
