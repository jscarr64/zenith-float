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
