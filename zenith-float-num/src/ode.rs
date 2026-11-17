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
    ExactNum::from_u8(1, p).ldexp(crate::defs::Exponent::from(ODE_MIN_STEP), p, rm)
}

/// Why an ODE call stopped.
///
/// Mirrors [`crate::QuadratureError`]: a distinct kind, not a silent `None`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OdeError {
    /// Non-finite values, `t1 ≤ t0`, a negative tolerance, a non-finite
    /// derivative or step, or a failed allocation of the output arrays.
    InvalidArgument,
    /// A probe step of size `(t1-t0)/`[`ODE_MAX_STEPS`] already misses the
    /// tolerance, or the accepted-step (or attempt) cap is hit before `t1`.
    /// Jeff ruling 2026-11-16 8:12 PM CT: this must be a clear error, not
    /// a silent `None` (ZF-54).
    ToleranceUnreachable,
}

impl core::fmt::Display for OdeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let repr = match self {
            Self::InvalidArgument => "invalid argument",
            Self::ToleranceUnreachable => "tolerance unreachable",
        };
        f.write_str(repr)
    }
}

#[cfg(feature = "std")]
impl std::error::Error for OdeError {}

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
    if !finite(t0) || !finite(y0) || !finite(t1) || t0.compare(t1) != Some(-1) {
        return None;
    }
    let wrk = work_p(p);
    let n_f = ExactNum::from_u32(n_steps as u32, wrk);
    let h = ExactNum::sub_ref(t1, t0, wrk, RoundingMode::None).div(&n_f, wrk, RoundingMode::None);
    let two = ExactNum::from_u8(2, wrk);
    let six = ExactNum::from_u8(6, wrk);
    let half_h = ExactNum::div_ref(&h, &two, wrk, RoundingMode::None);
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
        let y2 = ExactNum::mul_ref(&half_h, &k1, wrk, RoundingMode::None).add(
            &y,
            wrk,
            RoundingMode::None,
        );
        let t2 = ExactNum::add_ref(&t, &half_h, wrk, RoundingMode::None);
        let k2 = f(&t2, &y2, wrk, RoundingMode::None, cc);
        let y3 = ExactNum::mul_ref(&half_h, &k2, wrk, RoundingMode::None).add(
            &y,
            wrk,
            RoundingMode::None,
        );
        let k3 = f(&t2, &y3, wrk, RoundingMode::None, cc);
        let y4 =
            ExactNum::mul_ref(&h, &k3, wrk, RoundingMode::None).add(&y, wrk, RoundingMode::None);
        let t4 = t.add(&h, wrk, RoundingMode::None);
        let k4 = f(&t4, &y4, wrk, RoundingMode::None, cc);
        if !finite(&k1) || !finite(&k2) || !finite(&k3) || !finite(&k4) {
            return None;
        }
        let sum = k1
            .add(
                &ExactNum::from_u8(2, wrk).mul(&k2, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            )
            .add(
                &ExactNum::from_u8(2, wrk).mul(&k3, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            )
            .add(&k4, wrk, RoundingMode::None);
        y = y.add(
            &ExactNum::div_ref(&h, &six, wrk, RoundingMode::None).mul(
                &sum,
                wrk,
                RoundingMode::None,
            ),
            wrk,
            RoundingMode::None,
        );
        t = ExactNum::from_u32((k + 1) as u32, wrk)
            .mul(&h, wrk, RoundingMode::None)
            .add(t0, wrk, RoundingMode::None);
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
    if !finite(t0) || !finite(y0) || !finite(t1) || t0.compare(t1) != Some(-1) {
        return None;
    }
    let wrk = work_p(p);
    let n_f = ExactNum::from_u32(n_steps as u32, wrk);
    let h = ExactNum::sub_ref(t1, t0, wrk, RoundingMode::None).div(&n_f, wrk, RoundingMode::None);
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
            &yp.mul(&h, wrk, RoundingMode::None),
            wrk,
            RoundingMode::None,
        );
        t = ExactNum::from_u32((k + 1) as u32, wrk)
            .mul(&h, wrk, RoundingMode::None)
            .add(t0, wrk, RoundingMode::None);
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
    ExactNum::mul_ref(h, a, p, rm).mul(k, p, rm).add(y, p, rm)
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
    let t2 = ExactNum::mul_ref(h, &tab.c2, p, rm).add(t, p, rm);
    let y2 = axpy(y, h, &tab.a21, &k1, p, rm);
    let k2 = f(&t2, &y2, p, rm, cc);

    let t3 = ExactNum::mul_ref(h, &tab.c3, p, rm).add(t, p, rm);
    let y3 = axpy(&axpy(y, h, &tab.a31, &k1, p, rm), h, &tab.a32, &k2, p, rm);
    let k3 = f(&t3, &y3, p, rm, cc);

    let t4 = ExactNum::mul_ref(h, &tab.c4, p, rm).add(t, p, rm);
    let y4s = axpy(
        &axpy(&axpy(y, h, &tab.a41, &k1, p, rm), h, &tab.a42, &k2, p, rm),
        h,
        &tab.a43,
        &k3,
        p,
        rm,
    );
    let k4 = f(&t4, &y4s, p, rm, cc);

    let t5 = ExactNum::mul_ref(h, &tab.c5, p, rm).add(t, p, rm);
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

    let t6 = ExactNum::add_ref(t, h, p, rm);
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

fn step_scale(
    atol: &ExactNum,
    rtol: &ExactNum,
    y: &ExactNum,
    y5: &ExactNum,
    wrk: usize,
) -> ExactNum {
    let big = if y.abs_cmp(y5) == Some(1) { y } else { y5 };
    let negated;
    let ymax: &ExactNum = if big.is_negative() {
        negated = big.neg_via_mul();
        &negated
    } else {
        big
    };
    ExactNum::mul_ref(rtol, ymax, wrk, RoundingMode::None).add(atol, wrk, RoundingMode::None)
}

fn step_ok(err: &ExactNum, scale: &ExactNum) -> bool {
    scale.is_zero()
        || err.is_zero()
        || err.compare(scale) == Some(-1)
        || err.compare(scale) == Some(0)
}

/// Dormand–Prince RK5(4) with absolute / relative step control.
///
/// Accepts a step when `|y₅−y₄| ≤ atol + rtol max(|y|,|y₅|)`. Step size is
/// updated by `(atol_scale / err)^{1/5}` with safety `9/10`, growth cap 5,
/// shrink cap `1/5`.
///
/// Returns [`OdeError::InvalidArgument`] if `t1 ≤ t0`, a value is non-finite,
/// a tolerance is negative, or the output arrays cannot be allocated.
/// Returns [`OdeError::ToleranceUnreachable`] when a probe step of size
/// `(t1-t0)/`[`ODE_MAX_STEPS`] already misses the tolerance, when the
/// accepted-step or attempt cap is hit before `t1`, or when the step falls
/// below [`ode_min_step`] (ZF-54; Jeff ruling 2026-11-16 8:12 PM CT).
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
) -> Result<(ExactNumArray, ExactNumArray), OdeError>
where
    F: FnMut(&ExactNum, &ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    if !finite(t0) || !finite(y0) || !finite(t1) || !finite(atol) || !finite(rtol) {
        return Err(OdeError::InvalidArgument);
    }
    if t0.compare(t1) != Some(-1) || atol.is_negative() || rtol.is_negative() {
        return Err(OdeError::InvalidArgument);
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
    let span = ExactNum::sub_ref(t1, t0, wrk, RoundingMode::None);
    // Coarsest step that can still reach t1 in ODE_MAX_STEPS accepted steps.
    // If a Dormand–Prince step of that size already misses the tolerance,
    // smaller accepted steps cannot finish in the cap either (local error
    // shrinks like h^5). ZF-54: atol=2^{-512} on [0,1] at p=1024.
    let cap_n = ExactNum::from_u32(ODE_MAX_STEPS as u32, wrk);
    let h_budget = span.clone().div(&cap_n, wrk, RoundingMode::None);
    if h_budget.compare(&hmin) == Some(-1) || h_budget.is_zero() || h_budget.is_negative() {
        return Err(OdeError::ToleranceUnreachable);
    }
    let mut h = span.ldexp(-4, wrk, RoundingMode::None);
    let mut ts = Vec::new();
    let mut ys = Vec::new();
    if !push_pair(&mut ts, &mut ys, t.clone(), y.clone(), ODE_MAX_STEPS + 1) {
        return Err(OdeError::InvalidArgument);
    }
    // Probe at h_budget: if this step fails the error test, the cap will be
    // exhausted before t1.
    {
        let (y5, y4) = dp45_step(&tab, &mut f, &t, &y, &h_budget, wrk, RoundingMode::None, cc)
            .ok_or(OdeError::InvalidArgument)?;
        let err = ExactNum::sub_ref(&y5, &y4, wrk, RoundingMode::None).abs();
        let scale = step_scale(atol, rtol, &y, &y5, wrk);
        if !step_ok(&err, &scale) {
            return Err(OdeError::ToleranceUnreachable);
        }
    }
    let mut accepted = 0usize;
    let mut attempts = 0usize;
    while t.compare(t1) == Some(-1) {
        if accepted >= ODE_MAX_STEPS || attempts >= ODE_MAX_STEPS.saturating_mul(4) {
            return Err(OdeError::ToleranceUnreachable);
        }
        attempts += 1;
        let remain = ExactNum::sub_ref(t1, &t, wrk, RoundingMode::None);
        if h.compare(&remain) == Some(1) {
            h = remain;
        }
        if h.compare(&hmin) == Some(-1) || h.is_zero() || h.is_negative() {
            return Err(OdeError::ToleranceUnreachable);
        }
        let (y5, y4) = dp45_step(&tab, &mut f, &t, &y, &h, wrk, RoundingMode::None, cc)
            .ok_or(OdeError::InvalidArgument)?;
        let err = ExactNum::sub_ref(&y5, &y4, wrk, RoundingMode::None).abs();
        let scale = step_scale(atol, rtol, &y, &y5, wrk);
        let ok = step_ok(&err, &scale);
        if ok {
            t = t.add(&h, wrk, RoundingMode::None);
            y = y5;
            if t.compare(t1) == Some(1) {
                t = t1.clone();
            }
            if !push_pair(&mut ts, &mut ys, t.clone(), y.clone(), ODE_MAX_STEPS + 1) {
                return Err(OdeError::InvalidArgument);
            }
            accepted += 1;
        }
        let ratio = if err.is_zero() {
            grow.clone()
        } else {
            let q = scale
                .div(&err, wrk, RoundingMode::None)
                .nth_root(5, wrk, RoundingMode::None);
            let raw = q.mul(&fac, wrk, RoundingMode::None);
            if raw.compare(&grow) == Some(1) {
                grow.clone()
            } else if raw.compare(&shrink) == Some(-1) {
                shrink.clone()
            } else {
                raw
            }
        };
        h = h.mul(&ratio, wrk, RoundingMode::None);
        if !ok && !finite(&h) {
            return Err(OdeError::InvalidArgument);
        }
        let leftover = ODE_MAX_STEPS.saturating_sub(accepted);
        if leftover == 0 {
            return Err(OdeError::ToleranceUnreachable);
        }
        let remain = ExactNum::sub_ref(t1, &t, wrk, RoundingMode::None);
        if remain.is_positive() {
            let reach = h.clone().mul(
                &ExactNum::from_u32(leftover as u32, wrk),
                wrk,
                RoundingMode::None,
            );
            if remain.compare(&reach) == Some(1) {
                return Err(OdeError::ToleranceUnreachable);
            }
        }
    }
    let ts = to_row(p, rm, &ts).ok_or(OdeError::InvalidArgument)?;
    let ys = to_row(p, rm, &ys).ok_or(OdeError::InvalidArgument)?;
    Ok((ts, ys))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Consts;

    fn gold_p() -> (usize, RoundingMode) {
        (256, RoundingMode::ToEven)
    }

    /// Global RK4 error on `y'=-y` with `h=10^{-3}` is `O(h^4) ≈ 10^{-12}`.
    const RK4_EXP_ERR_DIGITS: i64 = 12;
    /// Adaptive gold the step cap can actually meet (plan `1e-50` needs ~10¹² steps).
    const RK45_ATOL_DIGITS: i64 = 12;
    const EULER_OH_STEPS: usize = 1000;

    #[test]
    fn ode_rk4_rk45_euler() -> Result<(), crate::Error> {
        let (p, rm) = gold_p();
        let mut cc = Consts::new()?;
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);
        let ten = ExactNum::from_u8(10, p);

        let (_t, y) = rk4(
            |_t, y, _p, _rm, _cc| y.neg_via_mul(),
            &zero,
            &one,
            &one,
            EULER_OH_STEPS,
            p,
            rm,
            &mut cc,
        )
        .ok_or(crate::Error::InvalidArgument)?;
        let yend = y.get(y.len() - 1).ok_or(crate::Error::InvalidArgument)?;
        let einv = one.neg_via_mul().exp(p, rm, &mut cc);
        let err = ExactNum::sub_ref(yend, &einv, p, rm).abs();
        let tol12 = ExactNum::from_u8(1, p).div(&ten.powsi(RK4_EXP_ERR_DIGITS, p, rm), p, rm);
        assert!(err.is_zero() || err.compare(&tol12) == Some(-1));

        let atol = ExactNum::from_u8(1, p).div(&ten.powsi(RK45_ATOL_DIGITS, p, rm), p, rm);
        let rtol = ExactNum::new(p);
        let (_t45, y45) = rk45_adaptive(
            |_t, y, _p, _rm, _cc| y.neg_via_mul(),
            &zero,
            &one,
            &one,
            &atol,
            &rtol,
            p,
            rm,
            &mut cc,
        )
        .map_err(|_| crate::Error::InvalidArgument)?;
        let y45e = y45
            .get(y45.len() - 1)
            .ok_or(crate::Error::InvalidArgument)?;
        let err45 = ExactNum::sub_ref(y45e, &einv, p, rm).abs();
        assert!(
            err45.is_zero() || err45.compare(&atol) == Some(-1) || err45.compare(&atol) == Some(0)
        );

        let (_te, ye) = euler(
            |_t, y, _p, _rm| y.clone(),
            &zero,
            &one,
            &one,
            EULER_OH_STEPS,
            p,
            rm,
        )
        .ok_or(crate::Error::InvalidArgument)?;
        let (_te2, ye2) = euler(
            |_t, y, _p, _rm| y.clone(),
            &zero,
            &one,
            &one,
            EULER_OH_STEPS * 2,
            p,
            rm,
        )
        .ok_or(crate::Error::InvalidArgument)?;
        assert_eq!(ye.len(), EULER_OH_STEPS + 1);
        assert_eq!(ye2.len(), EULER_OH_STEPS * 2 + 1);
        let two = ExactNum::from_u8(2, p);
        let ee = one.exp(p, rm, &mut cc);
        let e1 = ExactNum::sub_ref(
            ye.get(ye.len() - 1).ok_or(crate::Error::InvalidArgument)?,
            &ee,
            p,
            rm,
        )
        .abs();
        let e2 = ExactNum::sub_ref(
            ye2.get(ye2.len() - 1)
                .ok_or(crate::Error::InvalidArgument)?,
            &ee,
            p,
            rm,
        )
        .abs();
        let two_h_bound = two.div(&ExactNum::from_u32(EULER_OH_STEPS as u32, p), p, rm);
        assert!(e1.compare(&two_h_bound) == Some(-1));
        assert!(e2.compare(&e1) == Some(-1));

        assert!(rk4(
            |_t, y, _p, _rm, _cc| y.neg_via_mul(),
            &one,
            &one,
            &zero,
            4,
            p,
            rm,
            &mut cc
        )
        .is_none());
        Ok(())
    }

    /// ZF-54: atol=rtol=2^{-512} on [0,1] at p=1024 cannot finish inside
    /// [`ODE_MAX_STEPS`] (a 5th-order step of size 1/`ODE_MAX_STEPS` still
    /// misses the tolerance). The solver must return
    /// [`OdeError::ToleranceUnreachable`] without running the cap (65536
    /// Dormand–Prince steps, historically ~67 s). Jeff ruling 2026-11-16
    /// 8:12 PM CT: a clear error, not a silent `None`.
    #[test]
    fn zf_54_rk45_impossible_tol_refuses_early() -> Result<(), crate::Error> {
        let p = 1024;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);
        let tol = ExactNum::from_u8(1, p).ldexp(-512, p, rm);
        let mut evals = 0u32;
        let got = rk45_adaptive(
            |t, _y, q, r, cc| {
                evals += 1;
                t.cos(q, r, cc)
            },
            &zero,
            &one,
            &one,
            &tol,
            &tol,
            p,
            rm,
            &mut cc,
        );
        assert_eq!(got.err(), Some(OdeError::ToleranceUnreachable));
        assert!(
            evals <= 16,
            "refused only after {evals} f-evals; cap would be ~7*65536"
        );
        Ok(())
    }

    /// ZF-54: after the opening probe succeeds, a later stiff region forces
    /// the step below what the remaining budget can cover. That cap exit is
    /// the same [`OdeError::ToleranceUnreachable`] as the up-front probe.
    #[test]
    fn zf_54_rk45_step_cap_is_tolerance_unreachable() -> Result<(), crate::Error> {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);
        let half = ExactNum::from_u8(1, p).ldexp(-1, p, rm);
        let tol = ExactNum::from_u8(1, p).ldexp(-20, p, rm);
        let mut evals = 0u32;
        let got = rk45_adaptive(
            |t, y, q, r, _cc| {
                evals += 1;
                if t.compare(&half) == Some(-1) {
                    ExactNum::new(q)
                } else {
                    y.clone().ldexp(64, q, r)
                }
            },
            &zero,
            &one,
            &one,
            &tol,
            &tol,
            p,
            rm,
            &mut cc,
        );
        assert_eq!(got.err(), Some(OdeError::ToleranceUnreachable));
        assert!(
            evals <= 10_000,
            "cap exit only after {evals} f-evals; must not run the full 65536 steps"
        );
        Ok(())
    }
}
