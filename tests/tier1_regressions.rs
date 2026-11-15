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
        for rm in [RoundingMode::ToEven] {
            let z = ExactNum::new(p);
            let one = ExactNum::from_u8(1, p);
            let half = ExactNum::from_u8(1, p).ldexp(-1, p, rm);
            let two = ExactNum::from_u8(2, p);
            let sp = cc.pi(p + 128, RoundingMode::ToEven).sqrt(p, rm);
            let g = z.meijer_g(&[], &[], &[z.clone()], &[], p, rm, &mut cc);
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
