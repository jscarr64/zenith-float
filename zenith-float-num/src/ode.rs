//! Fixed-step and adaptive ODE solvers on [`ExactNum`].

use crate::defs::WORD_BIT_SIZE;
use crate::Consts;
use crate::ExactNum;
use crate::ExactNumArray;
use crate::RoundingMode;
use alloc::vec::Vec;

/// Minimum step at working precision `p`: `h_min = 2^{-p}`, the unit in the
/// last place of a unit-scale value. A step below this cannot be represented.
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

/// Minimum representable step `2^{-p}` at precision `p`.
pub fn ode_min_step(p: usize, rm: RoundingMode) -> ExactNum {
    let exp = i64::try_from(p).ok().map(|n| -n).unwrap_or(i64::MIN / 4);
    ExactNum::from_u8(1, p).ldexp(exp, p, rm)
}

/// Why an ODE call stopped.
///
/// Mirrors [`crate::QuadratureError`]: a distinct kind, not a silent `None`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OdeError {
    /// Non-finite values, `t1 ≤ t0`, a negative tolerance, a non-finite
    /// derivative or step, or a failed allocation of the output arrays.
    InvalidArgument,
    /// The step size has fallen below what the working precision can
    /// represent, the error estimate cannot decrease at the current
    /// precision (for example a zero tolerance with a method whose
    /// truncation term is nonzero), or an a priori order-5 bound from
    /// the opening step shows that the trajectory needed to meet the
    /// tolerance cannot be allocated or represented. A clear error,
    /// not a silent `None` (ZF-54). This is never a step counter.
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

fn push_pair(ts: &mut Vec<ExactNum>, ys: &mut Vec<ExactNum>, t: ExactNum, y: ExactNum) {
    ts.push(t);
    ys.push(y);
}

/// Classical RK4 for `y' = f(t, y)` on `[t0, t1]` with `n_steps` equal steps.
///
/// Returns row arrays `(t, y)` of length `n_steps+1`. `None` if `n_steps` is 0,
/// `t1 ≤ t0`, or a value is non-finite. A huge `n_steps` may fail to allocate.
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
    if n_steps == 0 {
        return None;
    }
    if !finite(t0) || !finite(y0) || !finite(t1) || t0.compare(t1) != Some(-1) {
        return None;
    }
    let wrk = work_p(p);
    let n_f = ExactNum::from_u64(n_steps as u64, wrk);
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
        t = ExactNum::from_u64((k + 1) as u64, wrk)
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
    if n_steps == 0 {
        return None;
    }
    if !finite(t0) || !finite(y0) || !finite(t1) || t0.compare(t1) != Some(-1) {
        return None;
    }
    let wrk = work_p(p);
    let n_f = ExactNum::from_u64(n_steps as u64, wrk);
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
        t = ExactNum::from_u64((k + 1) as u64, wrk)
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

fn wrk_supports_tol(atol: &ExactNum, rtol: &ExactNum, y: &ExactNum, wrk: usize) -> bool {
    let scale = step_scale(atol, rtol, y, y, wrk);
    if scale.is_zero() || !finite(&scale) {
        return false;
    }
    let floor = ode_min_step(wrk, RoundingMode::None);
    let padded = floor.ldexp(8, wrk, RoundingMode::None);
    matches!(padded.compare(&scale), Some(c) if c < 0)
}

/// Mantissa limbs of one [`ExactNum`] at precision `p`.
#[cfg(test)]
fn mantissa_words(p: usize) -> usize {
    p.div_ceil(WORD_BIT_SIZE).max(1)
}

/// Witness for the unsound 4e6c7e7 opening-probe rule, kept only so
/// the counterexample can show that rule would refuse a finishable
/// run. Not used in production.
#[cfg(test)]
fn max_representable_trajectory_len(p: usize) -> usize {
    let per = mantissa_words(p).saturating_mul(2);
    if per == 0 {
        usize::MAX
    } else {
        usize::MAX / per
    }
}

/// The 4e6c7e7 opening-probe rule, kept only as a test witness.
/// One DP5 step at `h0` estimates the step size needed at `t0`. That
/// is not a lower bound on the step count over the whole span, so
/// this must not be used as a production refusal.
#[cfg(test)]
fn a_priori_trajectory_unrepresentable<F>(
    tab: &Dp45,
    f: &mut F,
    t0: &ExactNum,
    y0: &ExactNum,
    h0: &ExactNum,
    span: &ExactNum,
    atol: &ExactNum,
    rtol: &ExactNum,
    wrk: usize,
    cc: &mut Consts,
) -> Result<bool, OdeError>
where
    F: FnMut(&ExactNum, &ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    let (y5, y4) = dp45_step(tab, f, t0, y0, h0, wrk, RoundingMode::None, cc)
        .ok_or(OdeError::InvalidArgument)?;
    let err = ExactNum::sub_ref(&y5, &y4, wrk, RoundingMode::None).abs();
    if err.is_zero() || !finite(&err) {
        return Ok(false);
    }
    let scale = step_scale(atol, rtol, y0, &y5, wrk);
    if scale.is_zero() || !finite(&scale) {
        return Ok(false);
    }
    if step_ok(&err, &scale) {
        return Ok(false);
    }
    let ratio = scale.div(&err, wrk, RoundingMode::None);
    if !finite(&ratio) || ratio.is_zero() || ratio.is_negative() {
        return Ok(false);
    }
    let root = ratio.nth_root(5, wrk, RoundingMode::None);
    if !finite(&root) || root.is_zero() || root.is_negative() {
        return Ok(false);
    }
    let h_req = h0.clone().mul(&root, wrk, RoundingMode::None);
    if !finite(&h_req) || h_req.is_zero() || h_req.is_negative() {
        return Ok(false);
    }
    let max_n = max_representable_trajectory_len(wrk);
    if max_n <= 1 {
        return Ok(true);
    }
    let cap = match i64::try_from(max_n - 1) {
        Ok(n) => ExactNum::from_i64(n, wrk),
        // `max_n` too large to name as `i64`: the representable
        // length is huge; do not refuse.
        Err(_) => return Ok(false),
    };
    let max_span = h_req.mul(&cap, wrk, RoundingMode::None);
    if !finite(&max_span) {
        return Ok(false);
    }
    Ok(span.clone().abs().compare(&max_span) == Some(1))
}

fn rk45_at_prec<F>(
    mut f: F,
    t0: &ExactNum,
    y0: &ExactNum,
    t1: &ExactNum,
    atol: &ExactNum,
    rtol: &ExactNum,
    wrk: usize,
    cc: &mut Consts,
) -> Result<(Vec<ExactNum>, Vec<ExactNum>), OdeError>
where
    F: FnMut(&ExactNum, &ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
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
    let mut h = span.clone().ldexp(-4, wrk, RoundingMode::None);
    let mut ts = Vec::new();
    let mut ys = Vec::new();
    push_pair(&mut ts, &mut ys, t.clone(), y.clone());
    while t.compare(t1) == Some(-1) {
        let remain = ExactNum::sub_ref(t1, &t, wrk, RoundingMode::None);
        if h.compare(&remain) == Some(1) {
            h = remain.clone();
        }
        let t_try = ExactNum::add_ref(&t, &h, wrk, RoundingMode::None);
        let no_progress = t_try.compare(&t) == Some(0)
            || h.is_zero()
            || h.is_negative()
            || !finite(&h)
            || h.compare(&hmin) == Some(-1);
        if no_progress {
            let yp = f(&t, &y, wrk, RoundingMode::None, cc);
            if !finite(&yp) {
                return Err(OdeError::InvalidArgument);
            }
            let leftover = remain.abs().mul(&yp.abs(), wrk, RoundingMode::None);
            let scale = step_scale(atol, rtol, &y, &y, wrk);
            if leftover.is_zero()
                || leftover.compare(&scale) == Some(-1)
                || leftover.compare(&scale) == Some(0)
            {
                t = t1.clone();
                push_pair(&mut ts, &mut ys, t.clone(), y.clone());
                break;
            }
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
            if !finite(&t) || !finite(&y) {
                return Err(OdeError::InvalidArgument);
            }
            push_pair(&mut ts, &mut ys, t.clone(), y.clone());
        } else {
            // Error stuck at the rounding floor cannot be driven down at this wrk.
            let floor = ode_min_step(wrk, RoundingMode::None).mul(
                &ExactNum::from_u8(1, wrk).add(&y.clone().abs(), wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            );
            if !err.is_zero()
                && (err.compare(&floor) == Some(-1) || err.compare(&floor) == Some(0))
                && err.compare(&scale) == Some(1)
            {
                return Err(OdeError::ToleranceUnreachable);
            }
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
    }
    Ok((ts, ys))
}

/// Dormand–Prince RK5(4) with absolute / relative step control.
///
/// Accepts a step when `|y₅−y₄| ≤ atol + rtol max(|y|,|y₅|)`. Step size is
/// updated by `(atol_scale / err)^{1/5}` with safety `9/10`, growth cap 5,
/// shrink cap `1/5`.
///
/// Returns [`OdeError::InvalidArgument`] if `t1 ≤ t0`, a value is non-finite,
/// a tolerance is negative, or the output arrays cannot be allocated.
/// Returns [`OdeError::ToleranceUnreachable`] when the step falls below
/// [`ode_min_step`] (what the working precision can represent), when the
/// error estimate cannot decrease at the current precision, or when an
/// a priori order-5 bound from the opening step shows that the needed
/// trajectory cannot be allocated or represented (ZF-54). The bound is
/// computed up front from the inputs and is not a step counter; it does
/// not fire when a run would finish. When the tolerance is finer than
/// the current working precision can support, the working precision is
/// raised and the integration is retried.
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
    if atol.is_zero() && rtol.is_zero() {
        return Err(OdeError::ToleranceUnreachable);
    }
    let mut wrk = work_p(p);
    loop {
        match rk45_at_prec(&mut f, t0, y0, t1, atol, rtol, wrk, cc) {
            Ok((ts, ys)) => {
                let ts = to_row(p, rm, &ts).ok_or(OdeError::InvalidArgument)?;
                let ys = to_row(p, rm, &ys).ok_or(OdeError::InvalidArgument)?;
                return Ok((ts, ys));
            }
            Err(OdeError::InvalidArgument) => return Err(OdeError::InvalidArgument),
            Err(OdeError::ToleranceUnreachable) => {
                if wrk_supports_tol(atol, rtol, y0, wrk) {
                    return Err(OdeError::ToleranceUnreachable);
                }
                let bump = crate::defs::WORD_BIT_SIZE.max(32);
                let next = wrk.saturating_add(bump);
                if next <= wrk {
                    return Err(OdeError::InvalidArgument);
                }
                wrk = next;
            }
        }
    }
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
    /// Adaptive gold: 10^{-12} is well above the representable step at p=256.
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

    /// ZF-54: a zero tolerance cannot be met by Dormand–Prince (nonzero
    /// truncation term). The solver returns
    /// [`OdeError::ToleranceUnreachable`] immediately, a clear error, not a
    /// silent `None`.
    #[test]
    fn zf_54_rk45_impossible_tol_refuses_early() -> Result<(), crate::Error> {
        let p = 64;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);
        let mut evals = 0u32;
        let got = rk45_adaptive(
            |t, _y, q, r, cc| {
                evals += 1;
                t.cos(q, r, cc)
            },
            &zero,
            &one,
            &one,
            &zero,
            &zero,
            p,
            rm,
            &mut cc,
        );
        assert_eq!(got.err(), Some(OdeError::ToleranceUnreachable));
        assert_eq!(evals, 0, "zero tolerance refuses before any step");
        Ok(())
    }

    /// The 4e6c7e7 opening-probe rule is not a whole-span proof: a
    /// sharp transient at t0 plus a long smooth tail can make one DP5
    /// step at h0 = span/16 predict more than max_n steps, while the
    /// true adaptive run shrinks through the layer and then grows h
    /// across the tail and finishes. The right-hand side is the
    /// non-stiff quadrature y' = λ exp(−λ t) (Lipschitz constant in y
    /// is 0). A stiff linear sink cannot be used here: explicit RK45
    /// would stay limited by the stability region.
    #[test]
    fn zf_54_probe_must_not_refuse_a_finishable_transient() -> Result<(), crate::Error> {
        let p = 64;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);
        let t1 = one.clone().ldexp(400, p, RoundingMode::None);
        let lam = ExactNum::from_u32(256, p);
        let atol = one.clone().ldexp(-16, p, RoundingMode::None);
        let rtol = ExactNum::new(p); // zero: scale stays at atol
        let wrk = super::work_p(p);
        let tab = super::Dp45::new(wrk, RoundingMode::None);
        let span = ExactNum::sub_ref(&t1, &zero, wrk, RoundingMode::None);
        let h0 = span.clone().ldexp(-4, wrk, RoundingMode::None);
        let mut y0 = zero.clone();
        let _ = y0.set_precision(wrk, RoundingMode::None);
        let mut t0 = zero.clone();
        let _ = t0.set_precision(wrk, RoundingMode::None);
        let mut lam_w = lam.clone();
        let _ = lam_w.set_precision(wrk, RoundingMode::None);
        let mut pulse = |tt: &ExactNum, _y: &ExactNum, q: usize, r: RoundingMode, c: &mut Consts| {
            let nt = lam_w.clone().mul(tt, q, r).neg_via_mul();
            lam_w.clone().mul(&nt.exp(q, r, c), q, r)
        };
        let decided = super::a_priori_trajectory_unrepresentable(
            &tab,
            &mut pulse,
            &t0,
            &y0,
            &h0,
            &span,
            &atol,
            &rtol,
            wrk,
            &mut cc,
        );
        assert_eq!(
            decided,
            Ok(true),
            "the 4e6c7e7 probe must predict N > max_n on this problem"
        );

        let mut evals = 0u32;
        let wall = std::time::Instant::now();
        let got = rk45_adaptive(
            |tt, _y, q, r, c| {
                evals += 1;
                let nt = lam.clone().mul(tt, q, r).neg_via_mul();
                lam.clone().mul(&nt.exp(q, r, c), q, r)
            },
            &zero,
            &zero,
            &t1,
            &atol,
            &rtol,
            p,
            rm,
            &mut cc,
        );
        let ms = wall.elapsed().as_millis();
        let (ts, ys) = match got {
            Ok(v) => v,
            Err(e) => panic!("adaptive run must finish, not {e:?}; evals={evals} in {ms} ms"),
        };
        let yend = ys.get(ys.len() - 1).ok_or(crate::Error::InvalidArgument)?;
        let tend = ts.get(ts.len() - 1).ok_or(crate::Error::InvalidArgument)?;
        assert_eq!(tend.compare(&t1), Some(0), "must land on t1; evals={evals}");
        // Exact y(t) = 1 - exp(-λ t) → 1 at t = 2^{400}.
        let err = ExactNum::sub_ref(yend, &one, p, rm).abs();
        let band = one.clone().ldexp(-12, p, RoundingMode::None);
        assert!(
            err.is_zero() || err.compare(&band) != Some(1),
            "y(t1) must be 1 within a generous multiple of atol; evals={evals}"
        );
        assert!(evals > 16, "the run must actually integrate; evals={evals}");
        assert!(evals < 100_000, "the tail must grow, not crawl; evals={evals}");
        assert!(ms < 2_000, "the finishable run must be quick; took {ms} ms, evals={evals}");
        let _ = ts;
        Ok(())
    }

    /// Tight but reachable tolerances on y' = cos(t) and y' = -y still
    /// return the correctly controlled values. The Harper ZF-54
    /// tolerance is not refused here: a single opening probe is not a
    /// whole-span proof, and no general substitute is claimed.
    #[test]
    fn zf_54_rk45_reachable_tight_tol_still_integrates() -> Result<(), crate::Error> {
        let rm = RoundingMode::ToEven;

        // Reachable tight tolerance: y' = cos(t), y(1) = 1 + sin(1).
        let p = 256;
        let mut cc = Consts::new()?;
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);
        let tol = ExactNum::from_u8(1, p).ldexp(-48, p, RoundingMode::None);
        let mut evals_r = 0u32;
        let (ts, ys) = rk45_adaptive(
            |t, _y, q, r, cc| {
                evals_r += 1;
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
        )
        .map_err(|_| crate::Error::InvalidArgument)?;
        let yend = ys.get(ys.len() - 1).ok_or(crate::Error::InvalidArgument)?;
        let tend = ts.get(ts.len() - 1).ok_or(crate::Error::InvalidArgument)?;
        assert_eq!(tend.compare(&one), Some(0));
        let expect = one.clone().add(&one.sin(p, rm, &mut cc), p, rm);
        let err = ExactNum::sub_ref(yend, &expect, p, rm).abs();
        // The accepted local scale is atol + rtol |y| ≈ 2^{-47}. Allow a
        // generous global multiple of that scale, still far below 2^{-40}.
        let band = ExactNum::from_u8(1, p).ldexp(-40, p, RoundingMode::None);
        assert!(
            err.is_zero() || err.compare(&band) != Some(1),
            "reachable tight tol must return 1+sin(1); evals={evals_r}"
        );
        assert!(evals_r > 16, "a finishable run actually integrates; evals={evals_r}");

        // Same reachable bar on a different f, so the bound is class-wide.
        let mut evals_e = 0u32;
        let ( _, ye) = rk45_adaptive(
            |_t, y, _q, _r, _cc| {
                evals_e += 1;
                y.neg_via_mul()
            },
            &zero,
            &one,
            &one,
            &tol,
            &tol,
            p,
            rm,
            &mut cc,
        )
        .map_err(|_| crate::Error::InvalidArgument)?;
        let yende = ye.get(ye.len() - 1).ok_or(crate::Error::InvalidArgument)?;
        let einv = one.neg_via_mul().exp(p, rm, &mut cc);
        let erre = ExactNum::sub_ref(yende, &einv, p, rm).abs();
        assert!(
            erre.is_zero() || erre.compare(&band) != Some(1),
            "y'=-y at the same tight tol must return e^{{-1}}; evals={evals_e}"
        );
        Ok(())
    }

    /// A caller-requested RK4 / Euler length above the old 65536 cap is a
    /// computable grid, not a refusal.
    #[test]
    fn caps_ode_fixed_step_no_length_cap() -> Result<(), crate::Error> {
        let p = 64;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);
        let n = 65_537usize;
        let (t, y) = rk4(
            |_t, y, _p, _rm, _cc| y.neg_via_mul(),
            &zero,
            &one,
            &one,
            n,
            p,
            rm,
            &mut cc,
        )
        .ok_or(crate::Error::InvalidArgument)?;
        assert_eq!(t.len(), n + 1);
        assert_eq!(y.len(), n + 1);
        let (te, ye) = euler(
            |_t, y, _p, _rm| y.neg_via_mul(),
            &zero,
            &one,
            &one,
            n,
            p,
            rm,
        )
        .ok_or(crate::Error::InvalidArgument)?;
        assert_eq!(te.len(), n + 1);
        assert_eq!(ye.len(), n + 1);
        Ok(())
    }
}
