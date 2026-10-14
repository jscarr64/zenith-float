//! Probe and regression for ε-extrapolation GK on endpoint singularities.
//!
//! For each candidate in `CASES` at 256 and 512 bits this prints the
//! estimate, `abs_err_est`, and `|true − estimate|`. Names in `CHOSEN` must
//! also be covered: `|true − estimate| ≤ abs_err_est` and the estimate within
//! `1e-6 · max(|true|, 1)`. Settings match the `gk-eps-*` rows in quad.jsonl.

use zenith_float_num::{integrate_adaptive_gk, Consts, ExactNum, Radix, RoundingMode};

const WIDTH: usize = 1088;

const PI_TEXT: &str = "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038";

const CASES: [(&str, &str); 7] = [
    ("xpow_m0p9", "10"),
    ("log_invsqrt", "-4"),
    ("xpow_m0p99", "100"),
    ("log2_invsqrt", "16"),
    ("log3_invsqrt", "-96"),
    ("arcsine", PI_TEXT),
    ("xpow_m0p75", "4"),
];

const CHOSEN: [&str; 6] = [
    "xpow_m0p9",
    "log_invsqrt",
    "xpow_m0p99",
    "log2_invsqrt",
    "log3_invsqrt",
    "xpow_m0p75",
];

fn power(x: &ExactNum, s: &str, p: usize, rm: RoundingMode, cc: &mut Consts) -> ExactNum {
    let s = ExactNum::parse(s, Radix::Dec, p + 128, RoundingMode::ToEven, cc);
    x.pow(&s, p, rm, cc)
}

fn inv_sqrt(x: &ExactNum, p: usize, rm: RoundingMode) -> ExactNum {
    x.reciprocal(p, rm).sqrt(p, rm)
}

fn eval(name: &str, x: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> ExactNum {
    if x.is_zero() {
        return ExactNum::nan(None);
    }
    match name {
        "xpow_m0p9" => power(x, "-0.9", p, rm, cc),
        "xpow_m0p99" => power(x, "-0.99", p, rm, cc),
        "xpow_m0p75" => power(x, "-0.75", p, rm, cc),
        "log_invsqrt" => x.ln(p, rm, cc).mul(&inv_sqrt(x, p, rm), p, rm),
        "log2_invsqrt" => {
            let l = x.ln(p, rm, cc);
            l.mul(&l, p, rm).mul(&inv_sqrt(x, p, rm), p, rm)
        }
        "log3_invsqrt" => {
            let l = x.ln(p, rm, cc);
            let l2 = l.mul(&l, p, rm);
            l2.mul(&l, p, rm).mul(&inv_sqrt(x, p, rm), p, rm)
        }
        "arcsine" => {
            let one = ExactNum::from_u8(1, p);
            let q = x.mul(&one.sub(x, p, rm), p, rm);
            if q.is_zero() {
                ExactNum::nan(None)
            } else {
                inv_sqrt(&q, p, rm)
            }
        }
        _ => ExactNum::nan(None),
    }
}

#[test]
fn gk_eps_candidates() {
    let Ok(mut cc) = Consts::new() else {
        panic!("constants");
    };
    let rm = RoundingMode::ToEven;
    let mut failures = Vec::new();
    for p in [256usize, 512] {
        let a = ExactNum::from_u8(0, p);
        let b = ExactNum::from_u8(1, p);
        let tol = ExactNum::parse("1e-8", Radix::Dec, p, rm, &mut cc);
        let rel = ExactNum::parse("1e-6", Radix::Dec, WIDTH, RoundingMode::Up, &mut cc);
        for (name, truth_text) in CASES {
            let truth =
                ExactNum::parse(truth_text, Radix::Dec, WIDTH, RoundingMode::ToEven, &mut cc);
            let got = integrate_adaptive_gk(
                |x, wp, wrm, cc| eval(name, x, wp, wrm, cc),
                &a,
                &b,
                &tol,
                32,
                p,
                rm,
                &mut cc,
            );
            let chosen = CHOSEN.contains(&name);
            match got {
                Ok(r) => {
                    let gap = r
                        .estimate_kronrod
                        .sub(&truth, WIDTH, RoundingMode::None)
                        .abs();
                    let covered = matches!(gap.cmp(&r.abs_err_est), Some(c) if c <= 0);
                    let one = ExactNum::from_u8(1, WIDTH);
                    let mag = truth.abs();
                    let scale = if mag.cmp(&one) == Some(1) { mag } else { one };
                    let bar = rel.mul(&scale, WIDTH, RoundingMode::Up);
                    let est_ok = matches!(gap.cmp(&bar), Some(c) if c <= 0);
                    println!(
                        "p={p} {name} est={} abs_err_est={} gap={} covered={covered} est_ok={est_ok}",
                        r.estimate_kronrod, r.abs_err_est, gap
                    );
                    if chosen && !(covered && est_ok) {
                        failures.push(format!("p={p} {name} covered={covered} est_ok={est_ok}"));
                    }
                }
                Err(err) => {
                    println!("p={p} {name} error={err:?}");
                    if chosen {
                        failures.push(format!("p={p} {name} error={err:?}"));
                    }
                }
            }
        }
    }
    assert!(
        failures.is_empty(),
        "gk-eps candidates:\n{}",
        failures.join("\n")
    );
}
