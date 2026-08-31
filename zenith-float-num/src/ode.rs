//! Fixed-step and adaptive ODE solvers on [`ExactNum`].

use crate::defs::WORD_BIT_SIZE;
use crate::Consts;
use crate::ExactNum;
use crate::ExactNumArray;
use crate::RoundingMode;
use alloc::vec::Vec;

/// Maximum accepted steps for any solver in this module.
pub const ODE_MAX_STEPS: usize = 65536;

/// Default minimum step exponent: `h_min = 2^{ODE_MIN_STEP}`.
pub const ODE_MIN_STEP: i32 = -256;

/// Safety factor numerator for Dormand–Prince step updates (`9/10`).
const ODE_FAC_NUM: i64 = 9;
/// Safety factor denominator for Dormand–Prince step updates.
const ODE_FAC_DEN: i64 = 10;
/// Maximum step growth.
const ODE_H_GROW: i64 = 5;
/// Maximum step shrink is `1/ODE_H_SHRINK_DEN`.
const ODE_H_SHRINK_DEN: i64 = 5;

fn work_p(p: usize) -> usize {
    p.saturating_add(WORD_BIT_SIZE)
}

fn finite(x: &ExactNum) -> bool {
    !x.is_nan() && !x.is_inf()
}

fn frac(n: i64, d: i64, p: usize, rm: RoundingMode) -> ExactNum {
    ExactNum::from_i64(n, p).div(&ExactNum::from_i64(d, p), p, rm)
}

/// Minimum step `2^{ODE_MIN_STEP}` at precision `p`.
pub fn ode_min_step(p: usize, rm: RoundingMode) -> ExactNum {
    ExactNum::from_u8(1, p).ldexp(ODE_MIN_STEP, p, rm)
}

fn to_row(p: usize, rm: RoundingMode, vals: &[ExactNum]) -> Option<ExactNumArray> {
    let rounded: Vec<ExactNum> = vals
        .iter()
        .map(|x| {
            let mut y = x.clone();
            let _ = y.set_precision(p, rm);
            y
        })
        .collect();
    ExactNumArray::from_shape(p, 1, rounded.len(), &rounded)
}

fn push_pair(
    ts: &mut Vec<ExactNum>,
    ys: &mut Vec<ExactNum>,
    t: ExactNum,
    y: ExactNum,
    cap: usize,
) -> bool {
    if ts.len() >= cap {
        return false;
    }
    ts.push(t);
    ys.push(y);
    true
}

/// Classical RK4 for `y' = f(t, y)` on `[t0, t1]` with `n_steps` equal steps.
///
/// Returns row arrays `(t, y)` of length `n_steps+1`. `None` if `n_steps` is 0
/// or greater than [`ODE_MAX_STEPS`], `t1 ≤ t0`, or a value is non-finite.
pub fn rk4<F>(
    mut f: F,
    t0: &ExactNum,
    y0: &ExactNum,
    t1: &ExactNum,
    n_steps: usize,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<(ExactNumArray, ExactNumArray)>
where
    F: FnMut(&ExactNum, &ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    if n_steps == 0 || n_steps > ODE_MAX_STEPS {
        return None;
    }
    if !finite(t0) || !finite(y0) || !finite(t1) || t0.cmp(t1) != Some(-1) {
        return None;
    }
    let wrk = work_p(p);
    let n_f = ExactNum::from_u32(n_steps as u32, wrk);
    let h = t1
        .sub(t0, wrk, RoundingMode::None)
        .div(&n_f, wrk, RoundingMode::None);
    let two = ExactNum::from_u8(2, wrk);
    let six = ExactNum::from_u8(6, wrk);
    let half_h = h.div(&two, wrk, RoundingMode::None);
    let mut t = t0.clone();
    let mut y = y0.clone();
    let _ = t.set_precision(wrk, RoundingMode::None);
    let _ = y.set_precision(wrk, RoundingMode::None);
    let mut ts = Vec::with_capacity(n_steps + 1);
    let mut ys = Vec::with_capacity(n_steps + 1);
    ts.push(t.clone());
    ys.push(y.clone());
    for k in 0..n_steps {
        let k1 = f(&t, &y, wrk, RoundingMode::None, cc);
        let y2 = y.add(
            &half_h.mul(&k1, wrk, RoundingMode::None),
            wrk,
            RoundingMode::None,
        );
        let t2 = t.add(&half_h, wrk, RoundingMode::None);
        let k2 = f(&t2, &y2, wrk, RoundingMode::None, cc);
        let y3 = y.add(
            &half_h.mul(&k2, wrk, RoundingMode::None),
            wrk,
            RoundingMode::None,
        );
        let k3 = f(&t2, &y3, wrk, RoundingMode::None, cc);
        let y4 = y.add(
            &h.mul(&k3, wrk, RoundingMode::None),
            wrk,
            RoundingMode::None,
        );
        let t4 = t.add(&h, wrk, RoundingMode::None);
        let k4 = f(&t4, &y4, wrk, RoundingMode::None, cc);
        if !finite(&k1) || !finite(&k2) || !finite(&k3) || !finite(&k4) {
            return None;
        }
        let sum = k1
            .add(
                &two.mul(&k2, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            )
            .add(
                &two.mul(&k3, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            )
            .add(&k4, wrk, RoundingMode::None);
        y = y.add(
            &h.div(&six, wrk, RoundingMode::None)
                .mul(&sum, wrk, RoundingMode::None),
            wrk,
            RoundingMode::None,
        );
        t = t0.add(
            &h.mul(
                &ExactNum::from_u32((k + 1) as u32, wrk),
                wrk,
                RoundingMode::None,
            ),
            wrk,
            RoundingMode::None,
        );
        if !finite(&t) || !finite(&y) {
            return None;
        }
        ts.push(t.clone());
        ys.push(y.clone());
    }
    Some((to_row(p, rm, &ts)?, to_row(p, rm, &ys)?))
}

/// Explicit Euler for `y' = f(t, y)` on `[t0, t1]` with `n_steps` equal steps.
///
/// Returns row arrays `(t, y)` of length `n_steps+1`. Same rejection as [`rk4`].
pub fn euler<F>(
    mut f: F,
    t0: &ExactNum,
    y0: &ExactNum,
    t1: &ExactNum,
    n_steps: usize,
    p: usize,
    rm: RoundingMode,
) -> Option<(ExactNumArray, ExactNumArray)>
where
    F: FnMut(&ExactNum, &ExactNum, usize, RoundingMode) -> ExactNum,
{
    if n_steps == 0 || n_steps > ODE_MAX_STEPS {
        return None;
    }
    if !finite(t0) || !finite(y0) || !finite(t1) || t0.cmp(t1) != Some(-1) {
        return None;
    }
    let wrk = work_p(p);
    let n_f = ExactNum::from_u32(n_steps as u32, wrk);
    let h = t1
        .sub(t0, wrk, RoundingMode::None)
        .div(&n_f, wrk, RoundingMode::None);
    let mut t = t0.clone();
    let mut y = y0.clone();
    let _ = t.set_precision(wrk, RoundingMode::None);
    let _ = y.set_precision(wrk, RoundingMode::None);
    let mut ts = Vec::with_capacity(n_steps + 1);
    let mut ys = Vec::with_capacity(n_steps + 1);
    ts.push(t.clone());
    ys.push(y.clone());
    for k in 0..n_steps {
        let yp = f(&t, &y, wrk, RoundingMode::None);
        if !finite(&yp) {
            return None;
        }
        y = y.add(
            &h.mul(&yp, wrk, RoundingMode::None),
            wrk,
            RoundingMode::None,
        );
        t = t0.add(
            &h.mul(
                &ExactNum::from_u32((k + 1) as u32, wrk),
                wrk,
                RoundingMode::None,
            ),
            wrk,
            RoundingMode::None,
        );
        if !finite(&t) || !finite(&y) {
            return None;
        }
        ts.push(t.clone());
        ys.push(y.clone());
    }
    Some((to_row(p, rm, &ts)?, to_row(p, rm, &ys)?))
}

struct Dp45 {
    c2: ExactNum,
    c3: ExactNum,
    c4: ExactNum,
    c5: ExactNum,
    a21: ExactNum,
    a31: ExactNum,
    a32: ExactNum,
    a41: ExactNum,
    a42: ExactNum,
    a43: ExactNum,
    a51: ExactNum,
    a52: ExactNum,
    a53: ExactNum,
    a54: ExactNum,
    a61: ExactNum,
    a62: ExactNum,
    a63: ExactNum,
    a64: ExactNum,
    a65: ExactNum,
    b5_1: ExactNum,
    b5_3: ExactNum,
    b5_4: ExactNum,
    b5_5: ExactNum,
    b5_6: ExactNum,
    b4_1: ExactNum,
    b4_3: ExactNum,
    b4_4: ExactNum,
    b4_5: ExactNum,
    b4_6: ExactNum,
    b4_7: ExactNum,
}

impl Dp45 {
    fn new(p: usize, rm: RoundingMode) -> Self {
        Self {
            c2: frac(1, 5, p, rm),
            c3: frac(3, 10, p, rm),
            c4: frac(4, 5, p, rm),
            c5: frac(8, 9, p, rm),
            a21: frac(1, 5, p, rm),
            a31: frac(3, 40, p, rm),
            a32: frac(9, 40, p, rm),
            a41: frac(44, 45, p, rm),
            a42: frac(-56, 15, p, rm),
            a43: frac(32, 9, p, rm),
            a51: frac(19372, 6561, p, rm),
            a52: frac(-25360, 2187, p, rm),
            a53: frac(64448, 6561, p, rm),
            a54: frac(-212, 729, p, rm),
            a61: frac(9017, 3168, p, rm),
            a62: frac(-355, 33, p, rm),
            a63: frac(46732, 5247, p, rm),
            a64: frac(49, 176, p, rm),
            a65: frac(-5103, 18656, p, rm),
            b5_1: frac(35, 384, p, rm),
            b5_3: frac(500, 1113, p, rm),
            b5_4: frac(125, 192, p, rm),
            b5_5: frac(-2187, 6784, p, rm),
            b5_6: frac(11, 84, p, rm),
            b4_1: frac(5179, 57600, p, rm),
            b4_3: frac(7571, 16695, p, rm),
            b4_4: frac(393, 640, p, rm),
            b4_5: frac(-92097, 339200, p, rm),
            b4_6: frac(187, 2100, p, rm),
            b4_7: frac(1, 40, p, rm),
        }
    }
}

fn axpy(
    y: &ExactNum,
    h: &ExactNum,
    a: &ExactNum,
    k: &ExactNum,
    p: usize,
    rm: RoundingMode,
) -> ExactNum {
    y.add(&h.mul(a, p, rm).mul(k, p, rm), p, rm)
}

fn dp45_step<F>(
    tab: &Dp45,
    f: &mut F,
    t: &ExactNum,
    y: &ExactNum,
    h: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<(ExactNum, ExactNum)>
where
    F: FnMut(&ExactNum, &ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    let k1 = f(t, y, p, rm, cc);
    let t2 = t.add(&h.mul(&tab.c2, p, rm), p, rm);
    let y2 = axpy(y, h, &tab.a21, &k1, p, rm);
    let k2 = f(&t2, &y2, p, rm, cc);

    let t3 = t.add(&h.mul(&tab.c3, p, rm), p, rm);
    let y3 = axpy(&axpy(y, h, &tab.a31, &k1, p, rm), h, &tab.a32, &k2, p, rm);
    let k3 = f(&t3, &y3, p, rm, cc);

    let t4 = t.add(&h.mul(&tab.c4, p, rm), p, rm);
    let y4s = axpy(
        &axpy(&axpy(y, h, &tab.a41, &k1, p, rm), h, &tab.a42, &k2, p, rm),
        h,
        &tab.a43,
        &k3,
        p,
        rm,
    );
    let k4 = f(&t4, &y4s, p, rm, cc);

    let t5 = t.add(&h.mul(&tab.c5, p, rm), p, rm);
    let y5s = axpy(
        &axpy(
            &axpy(&axpy(y, h, &tab.a51, &k1, p, rm), h, &tab.a52, &k2, p, rm),
            h,
            &tab.a53,
            &k3,
            p,
            rm,
        ),
        h,
        &tab.a54,
        &k4,
        p,
        rm,
    );
    let k5 = f(&t5, &y5s, p, rm, cc);

    let t6 = t.add(h, p, rm);
    let y6 = axpy(
        &axpy(
            &axpy(
                &axpy(&axpy(y, h, &tab.a61, &k1, p, rm), h, &tab.a62, &k2, p, rm),
                h,
                &tab.a63,
                &k3,
                p,
                rm,
            ),
            h,
            &tab.a64,
            &k4,
            p,
            rm,
        ),
        h,
        &tab.a65,
        &k5,
        p,
        rm,
    );
    let k6 = f(&t6, &y6, p, rm, cc);

    let y5 = axpy(
        &axpy(
            &axpy(
                &axpy(&axpy(y, h, &tab.b5_1, &k1, p, rm), h, &tab.b5_3, &k3, p, rm),
                h,
                &tab.b5_4,
                &k4,
                p,
                rm,
            ),
            h,
            &tab.b5_5,
            &k5,
            p,
            rm,
        ),
        h,
        &tab.b5_6,
        &k6,
        p,
        rm,
    );
    let k7 = f(&t6, &y5, p, rm, cc);

    let y4 = axpy(
        &axpy(
            &axpy(
                &axpy(
                    &axpy(&axpy(y, h, &tab.b4_1, &k1, p, rm), h, &tab.b4_3, &k3, p, rm),
                    h,
                    &tab.b4_4,
                    &k4,
                    p,
                    rm,
                ),
                h,
                &tab.b4_5,
                &k5,
                p,
                rm,
            ),
            h,
            &tab.b4_6,
            &k6,
            p,
            rm,
        ),
        h,
        &tab.b4_7,
        &k7,
        p,
        rm,
    );

    if !finite(&y5) || !finite(&y4) {
        return None;
    }
    Some((y5, y4))
}

/// Dormand–Prince RK5(4) with absolute / relative step control.
///
/// Accepts a step when `|y₅−y₄| ≤ atol + rtol max(|y|,|y₅|)`. Step size is
/// updated by `(atol_scale / err)^{1/5}` with safety `9/10`, growth cap 5,
/// shrink cap `1/5`. `None` if `t1 ≤ t0`, a value is non-finite, the step
/// falls below [`ode_min_step`], or [`ODE_MAX_STEPS`] accepted steps are
/// exhausted before `t1`.
pub fn rk45_adaptive<F>(
    mut f: F,
    t0: &ExactNum,
    y0: &ExactNum,
    t1: &ExactNum,
    atol: &ExactNum,
    rtol: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<(ExactNumArray, ExactNumArray)>
where
    F: FnMut(&ExactNum, &ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    if !finite(t0) || !finite(y0) || !finite(t1) || !finite(atol) || !finite(rtol) {
        return None;
    }
    if t0.cmp(t1) != Some(-1) || atol.is_negative() || rtol.is_negative() {
        return None;
    }
    let wrk = work_p(p);
    let tab = Dp45::new(wrk, RoundingMode::None);
    let hmin = ode_min_step(wrk, RoundingMode::None);
    let fac = frac(ODE_FAC_NUM, ODE_FAC_DEN, wrk, RoundingMode::None);
    let grow = ExactNum::from_i64(ODE_H_GROW, wrk);
    let shrink = frac(1, ODE_H_SHRINK_DEN, wrk, RoundingMode::None);
    let mut t = t0.clone();
    let mut y = y0.clone();
    let _ = t.set_precision(wrk, RoundingMode::None);
    let _ = y.set_precision(wrk, RoundingMode::None);
    let mut h = t1
        .sub(t0, wrk, RoundingMode::None)
        .ldexp(-4, wrk, RoundingMode::None);
    let mut ts = Vec::new();
    let mut ys = Vec::new();
    if !push_pair(&mut ts, &mut ys, t.clone(), y.clone(), ODE_MAX_STEPS + 1) {
        return None;
    }
    let mut accepted = 0usize;
    let mut attempts = 0usize;
    while t.cmp(t1) == Some(-1) {
        if accepted >= ODE_MAX_STEPS || attempts >= ODE_MAX_STEPS.saturating_mul(4) {
            return None;
        }
        attempts += 1;
        let remain = t1.sub(&t, wrk, RoundingMode::None);
        if h.cmp(&remain) == Some(1) {
            h = remain;
        }
        if h.cmp(&hmin) == Some(-1) || h.is_zero() || h.is_negative() {
            return None;
        }
        let (y5, y4) = dp45_step(&tab, &mut f, &t, &y, &h, wrk, RoundingMode::None, cc)?;
        let err = y5.sub(&y4, wrk, RoundingMode::None).abs();
        let ymax = if y.abs().cmp(&y5.abs()) == Some(1) { y.abs() } else { y5.abs() };
        let scale = atol.add(
            &rtol.mul(&ymax, wrk, RoundingMode::None),
            wrk,
            RoundingMode::None,
        );
        let ok = scale.is_zero()
            || err.is_zero()
            || err.cmp(&scale) == Some(-1)
            || err.cmp(&scale) == Some(0);
        if ok {
            t = t.add(&h, wrk, RoundingMode::None);
            y = y5;
            if t.cmp(t1) == Some(1) {
                t = t1.clone();
            }
            if !push_pair(&mut ts, &mut ys, t.clone(), y.clone(), ODE_MAX_STEPS + 1) {
                return None;
            }
            accepted += 1;
        }
        let ratio = if err.is_zero() {
            grow.clone()
        } else {
            let q = scale
                .div(&err, wrk, RoundingMode::None)
                .nth_root(5, wrk, RoundingMode::None);
            let raw = fac.mul(&q, wrk, RoundingMode::None);
            if raw.cmp(&grow) == Some(1) {
                grow.clone()
            } else if raw.cmp(&shrink) == Some(-1) {
                shrink.clone()
            } else {
                raw
            }
        };
        h = h.mul(&ratio, wrk, RoundingMode::None);
        if !ok && !finite(&h) {
            return None;
        }
    }
    Some((to_row(p, rm, &ts)?, to_row(p, rm, &ys)?))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Consts;

    fn gold_p() -> (usize, RoundingMode) {
        (256, RoundingMode::ToEven)
    }

    /// Global RK4 error on `y'=-y` with `h=10^{-3}` is `O(h^4) ≈ 10^{-12}`.
    const RK4_EXP_ERR_DIGITS: isize = 12;
    /// Adaptive gold the step cap can actually meet (plan `1e-50` needs ~10¹² steps).
    const RK45_ATOL_DIGITS: isize = 12;
    const EULER_OH_STEPS: usize = 1000;

    #[test]
    fn ode_rk4_rk45_euler() {
        let (p, rm) = gold_p();
        let mut cc = Consts::new().expect("consts");
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);
        let ten = ExactNum::from_u8(10, p);

        let (_t, y) = rk4(
            |_t, y, _p, _rm, _cc| y.neg(),
            &zero,
            &one,
            &one,
            EULER_OH_STEPS,
            p,
            rm,
            &mut cc,
        )
        .expect("rk4");
        let yend = y.get(y.len() - 1).expect("yend").clone();
        let einv = one.neg().exp(p, rm, &mut cc);
        let err = yend.sub(&einv, p, rm).abs();
        let tol12 = one.div(&ten.powsi(RK4_EXP_ERR_DIGITS, p, rm), p, rm);
        assert!(err.is_zero() || err.cmp(&tol12) == Some(-1));

        let atol = one.div(&ten.powsi(RK45_ATOL_DIGITS, p, rm), p, rm);
        let rtol = ExactNum::new(p);
        let (_t45, y45) = rk45_adaptive(
            |_t, y, _p, _rm, _cc| y.neg(),
            &zero,
            &one,
            &one,
            &atol,
            &rtol,
            p,
            rm,
            &mut cc,
        )
        .expect("rk45");
        let y45e = y45.get(y45.len() - 1).expect("y45e").clone();
        let err45 = y45e.sub(&einv, p, rm).abs();
        assert!(err45.is_zero() || err45.cmp(&atol) == Some(-1) || err45.cmp(&atol) == Some(0));

        let (_te, ye) = euler(
            |_t, y, _p, _rm| y.clone(),
            &zero,
            &one,
            &one,
            EULER_OH_STEPS,
            p,
            rm,
        )
        .expect("euler 1000");
        let (_te2, ye2) = euler(
            |_t, y, _p, _rm| y.clone(),
            &zero,
            &one,
            &one,
            EULER_OH_STEPS * 2,
            p,
            rm,
        )
        .expect("euler 2000");
        assert_eq!(ye.len(), EULER_OH_STEPS + 1);
        assert_eq!(ye2.len(), EULER_OH_STEPS * 2 + 1);
        let two = ExactNum::from_u8(2, p);
        let ee = one.exp(p, rm, &mut cc);
        let e1 = ye.get(ye.len() - 1).expect("e1").sub(&ee, p, rm).abs();
        let e2 = ye2.get(ye2.len() - 1).expect("e2").sub(&ee, p, rm).abs();
        let two_h_bound = two.div(&ExactNum::from_u32(EULER_OH_STEPS as u32, p), p, rm);
        assert!(e1.cmp(&two_h_bound) == Some(-1));
        assert!(e2.cmp(&e1) == Some(-1));

        assert!(rk4(
            |_t, y, _p, _rm, _cc| y.neg(),
            &one,
            &one,
            &zero,
            4,
            p,
            rm,
            &mut cc
        )
        .is_none());
    }
}
