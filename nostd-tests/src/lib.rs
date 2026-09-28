//! `no_std` reference-value tests for zenith-float-num.
//!
//! This crate is `#![no_std]` (core + alloc only) and depends on `zenith-float-num` with
//! `default-features = false`, so neither the library nor the harness links `std`. The cases in
//! [`cases`] carry mpmath references; [`run`] evaluates each one and reports the number of
//! correct bits through a caller-supplied sink. Two runners use it:
//!
//! - `nostd-host`: a `#![no_std]` / `#![no_main]` binary for hosted Linux targets that links only
//!   libc (allocator via `malloc`, output via `write`).
//! - `nostd-qemu`: a bare-metal `thumbv7em-none-eabihf` binary (Cortex-M4F) run under
//!   `qemu-system-arm -machine mps2-an386` with semihosting.
//!
//! See `scripts/ci_nostd.sh`.
#![no_std]

extern crate alloc;

pub mod cases;

use alloc::format;
use alloc::string::String;
use alloc::vec::Vec;
use cases::{Case, CASES};
use zenith_float_num::{
    brent, gauss_legendre, rk45_adaptive, root_default_tol, tanh_sinh, Consts, ExactComplex,
    ExactNum, Radix, RoundingMode,
};

const RM: RoundingMode = RoundingMode::ToEven;

/// Outcome of one case.
pub struct Outcome {
    /// Correct bits (`i32::MAX` if exact, `i32::MIN` if NaN / `None`).
    pub bits: i32,
    /// Whether `bits >= min_bits`.
    pub pass: bool,
}

fn parse(s: &str, p: usize, cc: &mut Consts) -> ExactNum {
    ExactNum::parse(s, Radix::Dec, p, RM, cc)
}

fn parse_c(s: &str, p: usize, cc: &mut Consts) -> ExactComplex {
    let (a, b) = s.split_once(',').unwrap_or((s, "0"));
    ExactComplex::new(parse(a, p, cc), parse(b, p, cc))
}

fn mag(x: &ExactNum) -> i32 {
    if x.is_zero() {
        i32::MIN
    } else {
        x.exponent().unwrap_or(i32::MIN)
    }
}

/// Correct bits of `got` against `want`, relative to the magnitude of `scale`.
fn bits(got: &ExactNum, want: &ExactNum, scale: &ExactNum, p: usize) -> i32 {
    if got.is_nan() || got.is_inf() {
        return i32::MIN;
    }
    let d = got.sub(want, 2 * p, RoundingMode::None);
    if d.is_zero() {
        return i32::MAX;
    }
    mag(scale).saturating_sub(d.exponent().unwrap_or(0))
}

fn real(op: &str, p: usize, a: &[ExactNum], cc: &mut Consts) -> ExactNum {
    let x = &a[0];
    match op {
        "sqrt" => x.sqrt(p, RM),
        "exp" => x.exp(p, RM, cc),
        "ln" => x.ln(p, RM, cc),
        "sin" => x.sin(p, RM, cc),
        "atan" => x.atan(p, RM, cc),
        "pow" => x.pow(&a[1], p, RM, cc),
        "erf" => x.erf(p, RM, cc),
        "erfc" => x.erfc(p, RM, cc),
        "gamma" => x.gamma(p, RM, cc),
        "ln_gamma" => x.ln_gamma(p, RM, cc),
        "digamma" => x.digamma(p, RM, cc),
        "hyp2f1" => x.hypergeom_2f1(&a[1], &a[2], &a[3], p, RM, cc),
        "betainc" => x.betainc(&a[1], &a[2], p, RM, cc),
        "gammainc" => x.gammainc(&a[1], p, RM, cc),
        "ei" => x.ei(p, RM, cc),
        "si" => x.si(p, RM, cc),
        "ci" => x.ci(p, RM, cc),
        "fresnel_s" => x.fresnel_s(p, RM, cc),
        "ai" => x.ai(p, RM, cc),
        "bi" => x.bi(p, RM, cc),
        "bessel_j" => x.bessel_j_nu(&a[1], p, RM, cc),
        "bessel_y" => x.bessel_y(&a[1], p, RM, cc),
        "bessel_k" => x.bessel_k(&a[1], p, RM, cc),
        "bessel_i" => x.bessel_i(&a[1], p, RM, cc),
        "elliptic_k" => x.elliptic_k(p, RM, cc),
        "elliptic_pi" => x.elliptic_pi(&a[1], &a[2], p, RM, cc),
        "jacobi_sn" => x.jacobi_sn(&a[1], p, RM, cc),
        "log" => x.log(&a[1], p, RM, cc),
        "normal_cdf" => x.normal_cdf(&a[1], &a[2], p, RM, cc),
        _ => ExactNum::nan(None),
    }
}

fn cplx(op: &str, p: usize, a: &[ExactComplex], cc: &mut Consts) -> ExactComplex {
    let z = &a[0];
    match op {
        "c.gamma" => z.gamma(p, RM, cc),
        "c.erf" => z.erf(p, RM, cc),
        "c.ci" => z.ci(p, RM, cc),
        "c.ai" => z.ai(p, RM, cc),
        "c.bessel_j" => z.bessel_j_nu(&a[1], p, RM, cc),
        "c.bessel_i" => z.bessel_i(&a[1], p, RM, cc),
        "c.hyp2f1" => z.hypergeom_2f1(&a[1], &a[2], &a[3], p, RM, cc),
        _ => ExactComplex::new(ExactNum::nan(None), ExactNum::nan(None)),
    }
}

fn method(op: &str, p: usize, cc: &mut Consts) -> ExactNum {
    let n = |v: u8| ExactNum::from_u8(v, p);
    let nan = || ExactNum::nan(None);
    match op {
        "m.gauss_legendre_recip_2_3" => gauss_legendre(
            |x, p, rm, _| ExactNum::from_u8(1, p).div(x, p, rm),
            &n(2),
            &n(3),
            30,
            p,
            RM,
            cc,
        )
        .unwrap_or_else(nan),
        "m.tanh_sinh_ln_0_1" => {
            tanh_sinh(|x, p, rm, cc| x.ln(p, rm, cc), &n(0), &n(1), p, RM, cc).unwrap_or_else(nan)
        }
        "m.brent_cos_x_eq_x" => brent(
            |x, p, rm, cc| x.cos(p, rm, cc).sub(x, p, rm),
            &n(0),
            &n(1),
            &root_default_tol(p, RM),
            p,
            RM,
            cc,
        )
        .unwrap_or_else(nan),
        "m.rk45_exp_0_1" => {
            let tol = parse("1e-12", p, cc);
            match rk45_adaptive(
                |_, y, _, _, _| y.clone(),
                &ExactNum::new(p),
                &n(1),
                &n(1),
                &tol,
                &tol,
                p,
                RM,
                cc,
            ) {
                Some((_, y)) => y.get(y.len() - 1).cloned().unwrap_or_else(nan),
                None => nan(),
            }
        }
        _ => nan(),
    }
}

/// Evaluates one case.
pub fn eval(case: &Case, cc: &mut Consts) -> Outcome {
    let p = case.p;
    let want_re = parse(case.re, 4 * p, cc);
    let bits = if case.op.starts_with("c.") {
        let zs: Vec<ExactComplex> = case.args.iter().map(|s| parse_c(s, p, cc)).collect();
        let want_im = parse(case.im, 4 * p, cc);
        let g = cplx(case.op, p, &zs, cc);
        let scale = if mag(&want_re) >= mag(&want_im) { &want_re } else { &want_im };
        bits(g.re(), &want_re, scale, p).min(bits(g.im(), &want_im, scale, p))
    } else if case.op.starts_with("m.") {
        bits(&method(case.op, p, cc), &want_re, &want_re, p)
    } else if case.op == "nth_root" {
        let x = parse(case.args[0], p, cc);
        let k: usize = case.args[1].parse().unwrap_or(0);
        bits(&x.nth_root(k, p, RM), &want_re, &want_re, p)
    } else {
        let xs: Vec<ExactNum> = case.args.iter().map(|s| parse(s, p, cc)).collect();
        bits(&real(case.op, p, &xs, cc), &want_re, &want_re, p)
    };
    Outcome {
        bits,
        pass: bits >= case.min_bits,
    }
}

/// Runs every case, sending one line per case to `sink`. Returns `(passed, failed)`.
pub fn run(sink: &mut dyn FnMut(&str)) -> (usize, usize) {
    let mut cc = match Consts::new() {
        Ok(c) => c,
        Err(_) => {
            sink("FAIL Consts::new");
            return (0, 1);
        }
    };
    let (mut ok, mut bad) = (0, 0);
    for case in CASES {
        let o = eval(case, &mut cc);
        let shown: String = if o.bits == i32::MAX {
            String::from("exact")
        } else if o.bits == i32::MIN {
            String::from("NaN")
        } else {
            format!("{}", o.bits)
        };
        let line = format!(
            "{} {} p={} {:?}: {} bits (need {})",
            if o.pass { "ok  " } else { "FAIL" },
            case.op,
            case.p,
            case.args,
            shown,
            case.min_bits
        );
        sink(&line);
        if o.pass {
            ok += 1;
        } else {
            bad += 1;
        }
    }
    (ok, bad)
}
