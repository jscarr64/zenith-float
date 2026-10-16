//! Quadrature, ODE, and ball goldens.
//!
//! `quad.jsonl` rows that return `GkResult` fail unless
//! `|true − estimate| ≤ abs_err_est`. Rows with `est_tol` also fail unless
//! the estimate itself is inside that relative tolerance. ODE rows fail unless the final value
//! stays inside the tolerance the solver was given. Ball rows fail unless
//! every listed true value lies in the returned ball. Value quadratures
//! (Legendre, Laguerre, Hermite, tanh–sinh) and Wynn have no error field;
//! those rows fail unless the true error is at most 64 ulps of the requested
//! precision (4 ulps for Wynn). The unit is `2^{e-p}`, the same unit
//! `set_precision(p)` rounds to.
//!
//! Fast tier is the default. `ORACLE_TIER=deep` also runs the deep rows.

use std::ops::Neg;
use zenith_float_num::{
    euler, gauss_hermite, gauss_kronrod_interval, gauss_laguerre, gauss_legendre,
    integrate_adaptive_gk, integrate_adaptive_gk_infinite, rk4, rk45_adaptive, tanh_sinh,
    wynn_epsilon, Ball, Consts, ExactNum, QuadratureError, Radix, RoundingMode,
};

/// `|x|` for an operand that stays borrowed: fresh `±1` times `x` at `x`'s own precision.
fn abs_of(x: &ExactNum) -> ExactNum {
    let p = x
        .mantissa_max_bit_len()
        .unwrap_or(zenith_float_num::WORD_BIT_SIZE);
    let sign = if x.is_negative() { -1 } else { 1 };
    ExactNum::from_i8(sign, p).mul(x, p, RoundingMode::None)
}

/// `-x` for an operand that stays borrowed: fresh `-1` times `x` at `x`'s own precision.
fn neg_of(x: &ExactNum) -> ExactNum {
    let p = x
        .mantissa_max_bit_len()
        .unwrap_or(zenith_float_num::WORD_BIT_SIZE);
    ExactNum::from_i8(-1, p).mul(x, p, RoundingMode::None)
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
        "Up" => RoundingMode::Up,
        "Down" => RoundingMode::Down,
        "ToZero" => RoundingMode::ToZero,
        "FromZero" => RoundingMode::FromZero,
        "ToOdd" => RoundingMode::ToOdd,
        "None" => RoundingMode::None,
        _ => RoundingMode::ToEven,
    }
}

fn parse_at(text: &str, p: usize, rm: RoundingMode, cc: &mut Consts) -> ExactNum {
    ExactNum::parse(text, Radix::Dec, p, rm, cc)
}

fn stored(p: usize) -> usize {
    p.div_ceil(64) * 64
}

fn sample(name: &str, x: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> ExactNum {
    match name {
        "exp" => x.exp(p, rm, cc),
        "sin" => x.sin(p, rm, cc),
        "one" => ExactNum::from_u8(1, p),
        "sin40" => {
            let forty = ExactNum::from_u8(40, p);
            forty.mul(x, p, rm).sin(p, rm, cc)
        }
        "invsqrt" => {
            if x.is_zero() {
                ExactNum::nan(None)
            } else {
                x.reciprocal(p, rm).sqrt(p, rm)
            }
        }
        "sq" => x.sqr(p, rm),
        "inv" => {
            if x.is_zero() {
                ExactNum::nan(None)
            } else {
                x.reciprocal(p, rm)
            }
        }
        "exp_neg" => neg_of(x).exp(p, rm, cc),
        "gauss" => x.sqr(p, rm).neg().exp(p, rm, cc),
        "invsq" => x.sqr(p, rm).reciprocal(p, rm),
        "xpow_m0p9" => nonzero(x, |x| power(x, "-0.9", p, rm, cc)),
        "xpow_m0p99" => nonzero(x, |x| power(x, "-0.99", p, rm, cc)),
        "xpow_m0p75" => nonzero(x, |x| power(x, "-0.75", p, rm, cc)),
        "log_invsqrt" => nonzero(x, |x| x.ln(p, rm, cc).mul(&inv_sqrt(x, p, rm), p, rm)),
        "log2_invsqrt" => nonzero(x, |x| {
            let l = x.ln(p, rm, cc);
            l.sqr(p, rm).mul(&inv_sqrt(x, p, rm), p, rm)
        }),
        "log3_invsqrt" => nonzero(x, |x| {
            let l = x.ln(p, rm, cc);
            let l2 = l.sqr(p, rm);
            l2.mul(&l, p, rm).mul(&inv_sqrt(x, p, rm), p, rm)
        }),
        _ => ExactNum::nan(None),
    }
}

fn nonzero(x: &ExactNum, f: impl FnOnce(&ExactNum) -> ExactNum) -> ExactNum {
    if x.is_zero() {
        ExactNum::nan(None)
    } else {
        f(x)
    }
}

fn power(x: &ExactNum, s: &str, p: usize, rm: RoundingMode, cc: &mut Consts) -> ExactNum {
    let s = ExactNum::parse(s, Radix::Dec, p + 128, RoundingMode::ToEven, cc);
    x.pow(&s, p, rm, cc)
}

fn inv_sqrt(x: &ExactNum, p: usize, rm: RoundingMode) -> ExactNum {
    x.reciprocal(p, rm).sqrt(p, rm)
}

fn covers(est: &ExactNum, truth: &ExactNum, err: &ExactNum) -> bool {
    let width = 1088usize;
    let gap = ExactNum::sub_ref(est, truth, width, RoundingMode::None).abs();
    matches!(gap.compare(err), Some(c) if c <= 0)
}

/// `2^{extra_bits}` units in the last place of a `p`-bit rounding of `est`.
fn ulps(est: &ExactNum, p: usize, extra_bits: i32) -> ExactNum {
    let width = est.mantissa_max_bit_len().unwrap_or(64).max(32);
    let e = est.exponent().unwrap_or(0);
    ExactNum::from_u8(1, width).ldexp(
        e.saturating_sub(p as i32).saturating_add(extra_bits),
        width,
        RoundingMode::Up,
    )
}

fn truth_of(row: &serde_json::Value, cc: &mut Consts) -> ExactNum {
    let text = row.get("true").and_then(|v| v.as_str()).unwrap_or("0");
    parse_at(text, 1088, RoundingMode::ToEven, cc)
}

fn run_quad(row: &serde_json::Value) -> Result<(), String> {
    let op = row.get("op").and_then(|v| v.as_str()).unwrap_or("");
    let p = row.get("p").and_then(|v| v.as_u64()).unwrap_or(53) as usize;
    let rm = mode_of(row.get("rm").and_then(|v| v.as_str()).unwrap_or("ToEven"));
    let mut cc = Consts::new().map_err(|_| "constants".to_string())?;
    let integ = row
        .get("integrand")
        .and_then(|v| v.as_str())
        .unwrap_or("one");
    if let Some(status) = row.get("status").and_then(|v| v.as_str()) {
        let a = parse_at(row["a"].as_str().unwrap_or("0"), p, rm, &mut cc);
        let b = parse_at(row["b"].as_str().unwrap_or("1"), p, rm, &mut cc);
        let tol = parse_at(row["tol"].as_str().unwrap_or("1e-6"), p, rm, &mut cc);
        let panels = row.get("panels").and_then(|v| v.as_u64()).unwrap_or(8) as usize;
        let got = integrate_adaptive_gk(
            |x, p, rm, cc| sample(integ, x, p, rm, cc),
            &a,
            &b,
            &tol,
            panels,
            p,
            rm,
            &mut cc,
        );
        let err = got.err().ok_or("expected NotConverged, got a value")?;
        if status == "NotConverged" && err == QuadratureError::NotConverged {
            return Ok(());
        }
        return Err(format!("status {err:?}"));
    }
    let truth = truth_of(row, &mut cc);
    let (est, err) = match op {
        "gk_adaptive" => {
            let a = parse_at(row["a"].as_str().unwrap_or("0"), p, rm, &mut cc);
            let b = parse_at(row["b"].as_str().unwrap_or("1"), p, rm, &mut cc);
            let tol = parse_at(row["tol"].as_str().unwrap_or("1e-6"), p, rm, &mut cc);
            let panels = row.get("panels").and_then(|v| v.as_u64()).unwrap_or(8) as usize;
            let got = integrate_adaptive_gk(
                |x, p, rm, cc| sample(integ, x, p, rm, cc),
                &a,
                &b,
                &tol,
                panels,
                p,
                rm,
                &mut cc,
            )
            .map_err(|e| format!("{e:?}"))?;
            (got.estimate_kronrod, got.abs_err_est)
        }
        "gk_infinite" => {
            let a = parse_at(row["a"].as_str().unwrap_or("0"), p, rm, &mut cc);
            let b = parse_at(row["b"].as_str().unwrap_or("inf"), p, rm, &mut cc);
            let tol = parse_at(row["tol"].as_str().unwrap_or("1e-6"), p, rm, &mut cc);
            let panels = row.get("panels").and_then(|v| v.as_u64()).unwrap_or(8) as usize;
            let got = integrate_adaptive_gk_infinite(
                |x, p, rm, cc| sample(integ, x, p, rm, cc),
                &a,
                &b,
                &tol,
                panels,
                p,
                rm,
                &mut cc,
            )
            .map_err(|e| format!("{e:?}"))?;
            (got.estimate_kronrod, got.abs_err_est)
        }
        "gk_interval" => {
            let a = parse_at(row["a"].as_str().unwrap_or("0"), p, rm, &mut cc);
            let b = parse_at(row["b"].as_str().unwrap_or("1"), p, rm, &mut cc);
            let got = gauss_kronrod_interval(
                |x, p, rm, cc| sample(integ, x, p, rm, cc),
                &a,
                &b,
                p,
                rm,
                &mut cc,
            )
            .map_err(|e| format!("{e:?}"))?;
            (got.estimate_kronrod, got.abs_err_est)
        }
        "gauss_legendre" | "tanh_sinh" | "gauss_laguerre" | "gauss_hermite" => {
            let n = row.get("n").and_then(|v| v.as_u64()).unwrap_or(8) as usize;
            let got = match op {
                "gauss_legendre" => {
                    let a = parse_at(row["a"].as_str().unwrap_or("-1"), p, rm, &mut cc);
                    let b = parse_at(row["b"].as_str().unwrap_or("1"), p, rm, &mut cc);
                    gauss_legendre(
                        |x, p, rm, cc| sample(integ, x, p, rm, cc),
                        &a,
                        &b,
                        n,
                        p,
                        rm,
                        &mut cc,
                    )
                }
                "tanh_sinh" => {
                    let a = parse_at(row["a"].as_str().unwrap_or("0"), p, rm, &mut cc);
                    let b = parse_at(row["b"].as_str().unwrap_or("1"), p, rm, &mut cc);
                    tanh_sinh(
                        |x, p, rm, cc| sample(integ, x, p, rm, cc),
                        &a,
                        &b,
                        p,
                        rm,
                        &mut cc,
                    )
                }
                "gauss_laguerre" => gauss_laguerre(
                    |x, p, rm, cc| sample(integ, x, p, rm, cc),
                    n,
                    p,
                    rm,
                    &mut cc,
                ),
                _ => gauss_hermite(
                    |x, p, rm, cc| sample(integ, x, p, rm, cc),
                    n,
                    p,
                    rm,
                    &mut cc,
                ),
            }
            .ok_or("quadrature returned None")?;
            let err = ulps(&got, p, 6);
            (got, err)
        }
        "wynn" => {
            let half = ExactNum::from_u8(1, p).ldexp(-1, p, rm);
            let partial = |n: usize| {
                let mut acc = ExactNum::new(p);
                let mut term = ExactNum::from_u8(1, p);
                for _ in 0..n {
                    acc = acc.add(&term, p, rm);
                    term = term.mul(&half, p, rm);
                }
                acc
            };
            let seq: Vec<ExactNum> = (1..=4).map(partial).collect();
            let got = wynn_epsilon(&seq, p, rm).ok_or("wynn returned None")?;
            let err = ulps(&got, p, 2);
            (got, err)
        }
        other => return Err(format!("unknown op {other}")),
    };
    if let Some(est_tol) = row.get("est_tol").and_then(|v| v.as_str()) {
        let rel = parse_at(est_tol, 384, RoundingMode::Up, &mut cc);
        let one = ExactNum::from_u8(1, 384);
        let mag = abs_of(&truth);
        let scale = if mag.compare(&one) == Some(1) { mag } else { one };
        let bar = rel.mul(&scale, 384, RoundingMode::Up);
        if !covers(&est, &truth, &bar) {
            return Err(format!(
                "estimate {est} misses accuracy {est_tol} around {truth}"
            ));
        }
    }
    if covers(&est, &truth, &err) {
        Ok(())
    } else {
        Err(format!("estimate {est} error {err} does not cover {truth}"))
    }
}

fn ode_rhs(
    problem: &str,
    t: &ExactNum,
    y: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> ExactNum {
    let _ = (p, rm, cc);
    if problem == "exp" {
        y.clone()
    } else {
        t.clone()
    }
}

fn run_ode(row: &serde_json::Value) -> Result<(), String> {
    let op = row.get("op").and_then(|v| v.as_str()).unwrap_or("");
    let problem = row.get("problem").and_then(|v| v.as_str()).unwrap_or("exp");
    let p = row.get("p").and_then(|v| v.as_u64()).unwrap_or(53) as usize;
    let rm = mode_of(row.get("rm").and_then(|v| v.as_str()).unwrap_or("ToEven"));
    let mut cc = Consts::new().map_err(|_| "constants".to_string())?;
    let t0 = parse_at(row["t0"].as_str().unwrap_or("0"), p, rm, &mut cc);
    let t1 = parse_at(row["t1"].as_str().unwrap_or("1"), p, rm, &mut cc);
    let y0 = parse_at(row["y0"].as_str().unwrap_or("1"), p, rm, &mut cc);
    let truth = truth_of(row, &mut cc);
    let ys = match op {
        "euler" => {
            let steps = row.get("steps").and_then(|v| v.as_u64()).unwrap_or(100) as usize;
            euler(
                |t, y, p, rm| ode_rhs(problem, t, y, p, rm, &mut cc),
                &t0,
                &y0,
                &t1,
                steps,
                p,
                rm,
            )
        }
        "rk4" => {
            let steps = row.get("steps").and_then(|v| v.as_u64()).unwrap_or(40) as usize;
            rk4(
                |t, y, p, rm, cc| ode_rhs(problem, t, y, p, rm, cc),
                &t0,
                &y0,
                &t1,
                steps,
                p,
                rm,
                &mut cc,
            )
        }
        "rk45" => {
            let atol = parse_at(row["atol"].as_str().unwrap_or("1e-8"), p, rm, &mut cc);
            let rtol = parse_at(row["rtol"].as_str().unwrap_or("1e-6"), p, rm, &mut cc);
            rk45_adaptive(
                |t, y, p, rm, cc| ode_rhs(problem, t, y, p, rm, cc),
                &t0,
                &y0,
                &t1,
                &atol,
                &rtol,
                p,
                rm,
                &mut cc,
            )
        }
        other => return Err(format!("unknown ode {other}")),
    }
    .ok_or("solver returned None")?;
    let last = ys.1.get(ys.1.len().saturating_sub(1)).ok_or("empty y")?;
    let err = if let Some(tol) = row.get("tol").and_then(|v| v.as_str()) {
        parse_at(tol, stored(p).max(64), RoundingMode::Up, &mut cc)
    } else {
        let atol = parse_at(
            row["atol"].as_str().unwrap_or("0"),
            384,
            RoundingMode::Up,
            &mut cc,
        );
        let rtol = parse_at(
            row["rtol"].as_str().unwrap_or("0"),
            384,
            RoundingMode::Up,
            &mut cc,
        );
        let scale = abs_of(&truth);
        let one = ExactNum::from_u8(1, 384);
        let mag = if scale.compare(&one) == Some(1) { scale } else { one };
        atol.add(
            &rtol.mul(&mag, 384, RoundingMode::Up),
            384,
            RoundingMode::Up,
        )
    };
    if covers(last, &truth, &err) {
        Ok(())
    } else {
        Err(format!("y {last} exceeds reported tolerance {err}"))
    }
}

fn run_ball(row: &serde_json::Value) -> Result<(), String> {
    let op = row.get("op").and_then(|v| v.as_str()).unwrap_or("");
    let p = row.get("p").and_then(|v| v.as_u64()).unwrap_or(53) as usize;
    let rm = mode_of(row.get("rm").and_then(|v| v.as_str()).unwrap_or("ToEven"));
    let mut cc = Consts::new().map_err(|_| "constants".to_string())?;
    let mid = parse_at(row["mid"].as_str().unwrap_or("0"), p, rm, &mut cc);
    let rad = parse_at(row["rad"].as_str().unwrap_or("0"), p, rm, &mut cc);
    let ball = Ball::new(mid, rad);
    let out = match op {
        "ball_exp" => ball.exp(p, rm, &mut cc),
        "ball_sin" => ball.sin(p, rm, &mut cc),
        "ball_add" | "ball_mul" => {
            let other = Ball::new(
                parse_at(row["other"].as_str().unwrap_or("0"), p, rm, &mut cc),
                parse_at(row["other_rad"].as_str().unwrap_or("0"), p, rm, &mut cc),
            );
            if op == "ball_add" {
                ball.add(&other, p, rm)
            } else {
                ball.mul(&other, p, rm)
            }
        }
        other => return Err(format!("unknown ball {other}")),
    };
    let points = row
        .get("true")
        .and_then(|v| v.as_array())
        .ok_or("missing true points")?;
    for point in points {
        let text = point.as_str().ok_or("true point")?;
        let truth = parse_at(text, 384, RoundingMode::ToEven, &mut cc);
        if !out.contains(&truth, 384) {
            return Err(format!(
                "{text} is outside mid {} rad {}",
                out.mid(),
                out.rad()
            ));
        }
    }
    Ok(())
}

fn run_file(name: &str, eval: fn(&serde_json::Value) -> Result<(), String>) {
    let path = format!(
        "{}/tests/oracle/goldens/{name}.jsonl",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(&path).unwrap_or_else(|err| panic!("read {path}: {err}"));
    let mut ran = 0usize;
    let mut failed = 0usize;
    let mut samples = Vec::new();
    for (idx, line) in text.lines().enumerate() {
        if line.is_empty() {
            continue;
        }
        let row: serde_json::Value = serde_json::from_str(line).unwrap_or_else(|err| {
            panic!("json {name}: {err}");
        });
        let tier = row.get("tier").and_then(|v| v.as_str()).unwrap_or("fast");
        if !tier_selected(tier) {
            continue;
        }
        ran += 1;
        let id = row.get("id").and_then(|v| v.as_str()).unwrap_or("?");
        if let Err(err) = eval(&row) {
            failed += 1;
            if samples.len() < 12 {
                samples.push(format!("{name}:{} {id}: {err}", idx + 1));
            }
        }
    }
    if ran == 0 {
        panic!("{name}.jsonl ran no rows");
    }
    if failed != 0 {
        panic!(
            "{name}.jsonl: {failed} failed / {ran}\n{}",
            samples.join("\n")
        );
    }
}

#[test]
fn quad() {
    run_file("quad", run_quad);
}

#[test]
fn ode() {
    run_file("ode", run_ode);
}

#[test]
fn balls() {
    run_file("balls", run_ball);
}
