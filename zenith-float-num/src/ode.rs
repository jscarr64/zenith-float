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
/// Maximum accepted Dormand–Prince steps under a *stable* order-5
/// error constant before `rk45_adaptive` switches to Bulirsch–Stoer.
///
/// After each trial the controller proposes `h ← h (9/10) (scale/err)^{1/5}`
/// (growth cap 5, shrink cap 1/5). After each accepted step the
/// remaining span would be covered in `N ≈ remain / h` accepted
/// steps of the proposed size. That projection is trusted only
/// when `h` is stable: over [`ODE_DP_STABILITY_WINDOW`] successive
/// accepted steps it has grown by less than
/// `2^{ODE_DP_STABILITY_GROW_BITS}`. Rejected steps reset the
/// window (the controller is still hunting). A finishable
/// transient grows `h` through the tail faster than that bound, so
/// it never meets the test with a large `N`. When `N` exceeds
/// this bound at a stable `h`, a 5th-order explicit method is the
/// wrong algorithm (more than `7 · 2^{14} ≈ 10^5` right-hand-side
/// evaluations on a constant-C problem). Measured in release,
/// rustc 1.99.1, 2026-11-19 16:54–16:55 CST America/Chicago, on
/// `y' = cos(t)` over a unit span (log
/// `$L/final/timeouts-11/rk45/bound-2e14.log`): `2^{14}` accepted
/// Dormand–Prince steps cost 969 ms / 1,386 ms / 4,127 ms at
/// wrk = 128 / 320 / 1088 (114,688 evaluations). One Bulirsch–Stoer
/// sweep at the same wrk and the tols used by the probe, the p=256
/// closed-form test, and Harper ZF-54 accepted in one `bs_try_step`
/// and cost <1 ms / 1 ms / 46 ms (15 / 120 / 1,295 evaluations).
/// Dormand–Prince at the bound is already ~90× (wrk = 1088) to
/// ~1,000× (wrk = 128) more expensive than one BS sweep; a higher
/// threshold would only widen that gap. The bound is a
/// method-selection threshold, not a wall-clock and not a
/// working-precision carve-out. Harper ZF-54 (`h ≲ 2^{-102}` on
/// `[0, 1]`) holds `h` steady with `N ≈ 2^{102}` and switches
/// after one window.
const ODE_DP_MAX_CONST_C_STEPS: u32 = 1 << 14;

/// Length of the accepted-step window used to decide whether the
/// order-5 constant is stable. Eight accepted steps is a short
/// method-selection prefix (`8 · 7` right-hand-side evaluations),
/// far below [`ODE_DP_MAX_CONST_C_STEPS`].
const ODE_DP_STABILITY_WINDOW: u32 = 8;

/// If `h` grows by at least `2^{ODE_DP_STABILITY_GROW_BITS}` over
/// one [`ODE_DP_STABILITY_WINDOW`], the controller is still
/// expanding and `N` will fall (mean ratio `8^{1/8} = 2^{3/8} ≈
/// 1.30` per step, well below the grow cap of 5 and above
/// equilibrium chatter). A finishable transient's tail grows
/// faster than this (the P9 probe grows `16×` in the same
/// window). A crawl at nearly-constant `h` grows by `1×`.
const ODE_DP_STABILITY_GROW_BITS: i64 = 3;

/// `true` when a stable order-5 step of size `h` would need more than
/// [`ODE_DP_MAX_CONST_C_STEPS`] accepted steps to cover `remain`.
fn const_c_step_count_exceeds(remain: &ExactNum, h: &ExactNum, wrk: usize) -> bool {
    if !finite(remain) || !finite(h) || h.is_zero() || h.is_negative() || remain.is_zero() {
        return false;
    }
    let n = remain.abs_ref().div(&h.abs_ref(), wrk, RoundingMode::None);
    if !finite(&n) {
        return true;
    }
    let cap = ExactNum::from_u32(ODE_DP_MAX_CONST_C_STEPS, wrk);
    matches!(n.compare(&cap), Some(1))
}

/// Maximum extrapolation rows in one Bulirsch–Stoer step (even sequence
/// `n = 2,4,...,2 K_MAX`, order `2 K_MAX`). Tighter remaining error
/// shrinks `H` and retries.
const ODE_BS_K_MAX: usize = 64;

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

/// The interval and initial value of [`rk4`] and [`rk45_adaptive`]: \(y(t_0) = y_0\), integrated to `t1`
/// (zenith-float 2.0: grouped so no public function takes more than 7 arguments).
#[derive(Clone, Copy, Debug)]
pub struct OdeSpan<'a> {
    /// Start abscissa.
    pub t0: &'a ExactNum,
    /// Initial value \(y(t_0)\).
    pub y0: &'a ExactNum,
    /// End abscissa (`t0 < t1`).
    pub t1: &'a ExactNum,
}

/// Absolute and relative step tolerances of [`rk45_adaptive`] (non-negative, not both zero).
#[derive(Clone, Copy, Debug)]
pub struct OdeTol<'a> {
    /// Absolute tolerance.
    pub atol: &'a ExactNum,
    /// Relative tolerance.
    pub rtol: &'a ExactNum,
}

fn to_row(p: usize, rm: RoundingMode, vals: &[ExactNum]) -> Option<ExactNumArray> {
    let rounded: Vec<ExactNum> = vals.iter().map(|x| x.rounded_to(p, rm)).collect();
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
    OdeSpan { t0, y0, t1 }: OdeSpan<'_>,
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
    let mut t = t0.rounded_to(wrk, RoundingMode::None);
    let mut y = y0.rounded_to(wrk, RoundingMode::None);
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
    let mut t = t0.rounded_to(wrk, RoundingMode::None);
    let mut y = y0.rounded_to(wrk, RoundingMode::None);
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
    (p, rm): (usize, RoundingMode),
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
    usize::MAX.checked_div(per).unwrap_or(usize::MAX)
}

/// The 4e6c7e7 opening-probe rule, kept only as a test witness.
/// One DP5 step at `h0` estimates the step size needed at `t0`. That
/// is not a lower bound on the step count over the whole span, so
/// this must not be used as a production refusal.
#[cfg(test)]
fn a_priori_trajectory_unrepresentable<F>(
    tab: &Dp45,
    f: &mut F,
    [t0, y0]: [&ExactNum; 2],
    [h0, span]: [&ExactNum; 2],
    [atol, rtol]: [&ExactNum; 2],
    wrk: usize,
    cc: &mut Consts,
) -> Result<bool, OdeError>
where
    F: FnMut(&ExactNum, &ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    let (y5, y4) = dp45_step(tab, f, t0, y0, h0, (wrk, RoundingMode::None), cc)
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
    let h_req = ExactNum::mul_ref(h0, &root, wrk, RoundingMode::None);
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
    Ok(span.abs_ref().compare(&max_span) == Some(1))
}

fn rk45_at_prec<F>(
    mut f: F,
    t0: &ExactNum,
    y0: &ExactNum,
    t1: &ExactNum,
    [atol, rtol]: [&ExactNum; 2],
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
    let mut t = t0.rounded_to(wrk, RoundingMode::None);
    let mut y = y0.rounded_to(wrk, RoundingMode::None);
    let span = ExactNum::sub_ref(t1, t0, wrk, RoundingMode::None);
    let mut h = span.ldexp_ref(-4, wrk, RoundingMode::None);
    let mut ts = Vec::new();
    let mut ys = Vec::new();
    let mut stable_n: u32 = 0;
    let mut stable_h0: Option<ExactNum> = None;
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
                push_pair(&mut ts, &mut ys, t, y);
                break;
            }
            return Err(OdeError::ToleranceUnreachable);
        }
        let (y5, y4) = dp45_step(&tab, &mut f, &t, &y, &h, (wrk, RoundingMode::None), cc)
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
                &ExactNum::from_u8(1, wrk).add(&y.abs_ref(), wrk, RoundingMode::None),
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
        let raw;
        let ratio: &ExactNum = if err.is_zero() {
            &grow
        } else {
            let q = scale
                .div(&err, wrk, RoundingMode::None)
                .nth_root(5, wrk, RoundingMode::None);
            raw = q.mul(&fac, wrk, RoundingMode::None);
            if raw.compare(&grow) == Some(1) {
                &grow
            } else if raw.compare(&shrink) == Some(-1) {
                &shrink
            } else {
                &raw
            }
        };
        let h_next = ExactNum::mul_ref(&h, ratio, wrk, RoundingMode::None);
        let h_used = core::mem::replace(&mut h, h_next);
        if !ok && !finite(&h) {
            return Err(OdeError::InvalidArgument);
        }
        // Rejected steps reset the stability window (still hunting).
        // After a full window of accepts, a small net growth in `h`
        // plus a large constant-C projection means a crawl; a
        // finishable tail grows `h` faster than the documented bound.
        if ok {
            if stable_n == 0 {
                stable_h0 = Some(h_used.clone());
            }
            stable_n = stable_n.saturating_add(1);
            if stable_n >= ODE_DP_STABILITY_WINDOW {
                let grew = match stable_h0.as_ref() {
                    Some(h0) => {
                        let mark =
                            h0.ldexp_ref(ODE_DP_STABILITY_GROW_BITS, wrk, RoundingMode::None);
                        matches!(h_used.compare(&mark), Some(c) if c >= 0)
                    }
                    None => true,
                };
                if !grew {
                    let rem = ExactNum::sub_ref(t1, &t, wrk, RoundingMode::None);
                    if const_c_step_count_exceeds(&rem, &h, wrk) {
                        return Err(OdeError::ToleranceUnreachable);
                    }
                }
                stable_n = 0;
                stable_h0 = None;
            }
        } else {
            stable_n = 0;
            stable_h0 = None;
        }
    }
    Ok((ts, ys))
}

/// Modified midpoint (Gragg) over one span `H` with `n` substeps.
fn modified_midpoint<F>(
    f: &mut F,
    t: &ExactNum,
    y: &ExactNum,
    span: &ExactNum,
    n: u64,
    wrk: usize,
    cc: &mut Consts,
) -> Result<ExactNum, OdeError>
where
    F: FnMut(&ExactNum, &ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    if n == 0 {
        return Err(OdeError::InvalidArgument);
    }
    let rm = RoundingMode::None;
    let h = ExactNum::div_ref(span, &ExactNum::from_u64(n, wrk), wrk, rm);
    if !finite(&h) || h.is_zero() {
        return Err(OdeError::InvalidArgument);
    }
    let two_h = h.ldexp_ref(1, wrk, rm);
    let yp = f(t, y, wrk, rm, cc);
    if !finite(&yp) {
        return Err(OdeError::InvalidArgument);
    }
    // z_{m-1}; `None` while it is the start value y (read in place).
    let mut z_prev_slot: Option<ExactNum> = None;
    let mut z = ExactNum::add_ref(y, &ExactNum::mul_ref(&h, &yp, wrk, rm), wrk, rm);
    if !finite(&z) {
        return Err(OdeError::InvalidArgument);
    }
    for i in 1..n {
        let ti = ExactNum::add_ref(
            t,
            &ExactNum::mul_ref(&h, &ExactNum::from_u64(i, wrk), wrk, rm),
            wrk,
            rm,
        );
        let fi = f(&ti, &z, wrk, rm, cc);
        if !finite(&fi) {
            return Err(OdeError::InvalidArgument);
        }
        let z_prev = z_prev_slot.as_ref().unwrap_or(y);
        let z_next = ExactNum::add_ref(z_prev, &ExactNum::mul_ref(&two_h, &fi, wrk, rm), wrk, rm);
        if !finite(&z_next) {
            return Err(OdeError::InvalidArgument);
        }
        z_prev_slot = Some(core::mem::replace(&mut z, z_next));
    }
    let t_end = ExactNum::add_ref(t, span, wrk, rm);
    let f_end = f(&t_end, &z, wrk, rm, cc);
    if !finite(&f_end) {
        return Err(OdeError::InvalidArgument);
    }
    // (z_prev + z + h f(t+H, z)) / 2
    let z_prev = z_prev_slot.as_ref().unwrap_or(y);
    let sum = ExactNum::add_ref(z_prev, &z, wrk, rm).add(
        &ExactNum::mul_ref(&h, &f_end, wrk, rm),
        wrk,
        rm,
    );
    Ok(sum.ldexp(-1, wrk, rm))
}

/// One Bulirsch–Stoer attempt over `H`: modified midpoint with
/// `n = 2,4,...,2 K` and a Neville tableau in `H²`. Returns `Ok(Some(y))`
/// when the last-column increment is within `atol + rtol max(|y0|,|y|)`.
/// Returns `Ok(None)` if the tableau is exhausted without meeting the
/// tolerance (caller shrinks `H`).
fn bs_try_step<F>(
    f: &mut F,
    t: &ExactNum,
    y: &ExactNum,
    span: &ExactNum,
    [atol, rtol]: [&ExactNum; 2],
    wrk: usize,
    cc: &mut Consts,
) -> Result<Option<ExactNum>, OdeError>
where
    F: FnMut(&ExactNum, &ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    let rm = RoundingMode::None;
    let mut row: Vec<ExactNum> = Vec::with_capacity(ODE_BS_K_MAX);
    for k in 0..ODE_BS_K_MAX {
        let n = 2 * (k as u64 + 1);
        let t0k = modified_midpoint(f, t, y, span, n, wrk, cc)?;
        let mut cur = Vec::with_capacity(k + 1);
        cur.push(t0k);
        for j in 1..=k {
            let n_k = 2 * (k as u64 + 1);
            let n_kj = 2 * ((k - j) as u64 + 1);
            // den = (n_k / n_{k-j})² - 1
            let ratio = ExactNum::div_ref(
                &ExactNum::from_u64(n_k, wrk),
                &ExactNum::from_u64(n_kj, wrk),
                wrk,
                rm,
            );
            let den =
                ExactNum::mul_ref(&ratio, &ratio, wrk, rm).sub(&ExactNum::from_u8(1, wrk), wrk, rm);
            if den.is_zero() || !finite(&den) {
                return Err(OdeError::InvalidArgument);
            }
            let diff = ExactNum::sub_ref(&cur[j - 1], &row[j - 1], wrk, rm);
            let nxt = ExactNum::add_ref(
                &cur[j - 1],
                &ExactNum::div_ref(&diff, &den, wrk, rm),
                wrk,
                rm,
            );
            if !finite(&nxt) {
                return Err(OdeError::InvalidArgument);
            }
            cur.push(nxt);
        }
        if k >= 1 {
            let err = ExactNum::sub_ref(&cur[k], &cur[k - 1], wrk, rm).abs();
            let scale = step_scale(atol, rtol, y, &cur[k], wrk);
            if step_ok(&err, &scale) {
                // cur[k] is the last entry: move it out as the result.
                return Ok(Some(cur.swap_remove(k)));
            }
            // If the column error has grown twice in a row, this H will
            // not meet the tolerance: shrink it rather than filling the
            // remaining rows.
            if k >= 3 {
                let err_prev = ExactNum::sub_ref(&cur[k - 1], &cur[k - 2], wrk, rm).abs();
                let err_prev2 = ExactNum::sub_ref(&row[k - 2], &row[k - 3], wrk, rm).abs();
                if err.compare(&err_prev) != Some(-1) && err_prev.compare(&err_prev2) != Some(-1) {
                    return Ok(None);
                }
            }
        }
        row = cur;
    }
    Ok(None)
}

/// Adaptive Bulirsch–Stoer from `t0` to `t1` at working precision `wrk`.
fn bs_at_prec<F>(
    f: &mut F,
    t0: &ExactNum,
    y0: &ExactNum,
    t1: &ExactNum,
    [atol, rtol]: [&ExactNum; 2],
    wrk: usize,
    cc: &mut Consts,
) -> Result<(Vec<ExactNum>, Vec<ExactNum>), OdeError>
where
    F: FnMut(&ExactNum, &ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    let rm = RoundingMode::None;
    let hmin = ode_min_step(wrk, rm);
    let mut t = t0.rounded_to(wrk, rm);
    let mut y = y0.rounded_to(wrk, rm);
    let mut h = ExactNum::sub_ref(t1, t0, wrk, rm);
    let mut ts = Vec::new();
    let mut ys = Vec::new();
    push_pair(&mut ts, &mut ys, t.clone(), y.clone());
    while t.compare(t1) == Some(-1) {
        let remain = ExactNum::sub_ref(t1, &t, wrk, rm);
        if h.compare(&remain) == Some(1) {
            h = remain;
        }
        if !finite(&h) || h.is_zero() || h.is_negative() || h.compare(&hmin) == Some(-1) {
            return Err(OdeError::ToleranceUnreachable);
        }
        match bs_try_step(f, &t, &y, &h, [atol, rtol], wrk, cc)? {
            Some(yn) => {
                t = ExactNum::add_ref(&t, &h, wrk, rm);
                if t.compare(t1) == Some(1) {
                    t = t1.clone();
                }
                if !finite(&t) || !finite(&yn) {
                    return Err(OdeError::InvalidArgument);
                }
                y = yn;
                push_pair(&mut ts, &mut ys, t.clone(), y.clone());
                // Grow cautiously after an accepted step.
                h = h.ldexp(1, wrk, rm);
            }
            None => {
                h = h.ldexp(-1, wrk, rm);
            }
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
/// [`ode_min_step`] (what the working precision can represent) or when
/// the error estimate cannot decrease at the current precision, and
/// Bulirsch–Stoer extrapolation cannot reach the tolerance either; or
/// when both tolerances are zero (Dormand–Prince's truncation term is
/// nonzero). When the tolerance is finer than the current working
/// precision can support, the working precision is raised and the
/// integration is retried.
///
/// Switch (Jeff, 2026-11-19 11:50 AM CT): "rk45 switches to
/// Bulirsch–Stoer extrapolation for tolerances Dormand–Prince can't
/// reach, documented". If Dormand–Prince returns
/// [`OdeError::ToleranceUnreachable`] at a working precision that
/// already supports the tolerance — including when a stable
/// order-5 step (net growth of `h` below
/// `2^{ODE_DP_STABILITY_GROW_BITS}` over
/// [`ODE_DP_STABILITY_WINDOW`] accepted steps) would need more
/// than [`ODE_DP_MAX_CONST_C_STEPS`] accepted steps on the
/// remaining span — the same call continues with a
/// modified-midpoint
/// / polynomial extrapolation (Gragg–Bulirsch–Stoer, even sequence
/// `n = 2, 4, …, 2 K`) rather than refusing. The switch is
/// deterministic and derived from Dormand–Prince's error model
/// (`err ≈ C h^5`, `h ← (9/10) (scale/err)^{1/5}`): it fires only
/// when `h` is stable over [`ODE_DP_STABILITY_WINDOW`] accepted
/// steps and the implied step would still crawl, and it does not
/// use a fixed attempt count or a working-precision carve-out. A
/// finishable transient grows `h` through the tail faster than
/// `2^{ODE_DP_STABILITY_GROW_BITS}` per window, so it stays on
/// Dormand–Prince
/// (`zf_54_probe_must_not_refuse_a_finishable_transient`). A smooth
/// problem whose equilibrium step is too small (Harper ZF-54,
/// `h ≲ 2^{-102}` on `[0, 1]`) reaches equilibrium with
/// `N ≈ 2^{102}` and switches.
pub fn rk45_adaptive<F>(
    mut f: F,
    OdeSpan { t0, y0, t1 }: OdeSpan<'_>,
    OdeTol { atol, rtol }: OdeTol<'_>,
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
        match rk45_at_prec(&mut f, t0, y0, t1, [atol, rtol], wrk, cc) {
            Ok((ts, ys)) => {
                let ts = to_row(p, rm, &ts).ok_or(OdeError::InvalidArgument)?;
                let ys = to_row(p, rm, &ys).ok_or(OdeError::InvalidArgument)?;
                return Ok((ts, ys));
            }
            Err(OdeError::InvalidArgument) => return Err(OdeError::InvalidArgument),
            Err(OdeError::ToleranceUnreachable) => {
                if wrk_supports_tol(atol, rtol, y0, wrk) {
                    match bs_at_prec(&mut f, t0, y0, t1, [atol, rtol], wrk, cc) {
                        Ok((ts, ys)) => {
                            let ts = to_row(p, rm, &ts).ok_or(OdeError::InvalidArgument)?;
                            let ys = to_row(p, rm, &ys).ok_or(OdeError::InvalidArgument)?;
                            return Ok((ts, ys));
                        }
                        Err(e) => return Err(e),
                    }
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
    use crate::test_must::Must;
    use crate::Consts;
    use core::ops::Neg;

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
            OdeSpan {
                t0: &zero,
                y0: &one,
                t1: &one,
            },
            EULER_OH_STEPS,
            p,
            rm,
            &mut cc,
        )
        .ok_or(crate::Error::InvalidArgument)?;
        let yend = y.get(y.len() - 1).ok_or(crate::Error::InvalidArgument)?;
        let einv = ExactNum::from_i8(-1, p).exp(p, rm, &mut cc);
        let err = ExactNum::sub_ref(yend, &einv, p, rm).abs();
        let tol12 = ExactNum::from_u8(1, p).div(&ten.powsi(RK4_EXP_ERR_DIGITS, p, rm), p, rm);
        assert!(err.is_zero() || err.compare(&tol12) == Some(-1));

        let atol = ExactNum::from_u8(1, p).div(&ten.powsi(RK45_ATOL_DIGITS, p, rm), p, rm);
        let rtol = ExactNum::new(p);
        let (_t45, y45) = rk45_adaptive(
            |_t, y, _p, _rm, _cc| y.neg_via_mul(),
            OdeSpan {
                t0: &zero,
                y0: &one,
                t1: &one,
            },
            OdeTol {
                atol: &atol,
                rtol: &rtol,
            },
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
            OdeSpan {
                t0: &one,
                y0: &one,
                t1: &zero
            },
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
            OdeSpan {
                t0: &zero,
                y0: &one,
                t1: &one,
            },
            OdeTol {
                atol: &zero,
                rtol: &zero,
            },
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
    fn zf_54_probe_must_not_refuse_a_finishable_transient() -> Result<(), String> {
        let p = 64;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);
        let t1 = ExactNum::from_u8(1, p).ldexp(400, p, RoundingMode::None);
        let lam = ExactNum::from_u32(256, p);
        let atol = ExactNum::from_u8(1, p).ldexp(-16, p, RoundingMode::None);
        let rtol = ExactNum::new(p); // zero: scale stays at atol
        let wrk = super::work_p(p);
        let tab = super::Dp45::new(wrk, RoundingMode::None);
        let span = ExactNum::sub_ref(&t1, &zero, wrk, RoundingMode::None);
        let h0 = span.ldexp_ref(-4, wrk, RoundingMode::None);
        let mut y0 = ExactNum::new(p);
        let _ = y0.set_precision(wrk, RoundingMode::None);
        let mut t0 = ExactNum::new(p);
        let _ = t0.set_precision(wrk, RoundingMode::None);
        let lam_w = lam.rounded_to(wrk, RoundingMode::None);
        let mut pulse =
            |tt: &ExactNum, _y: &ExactNum, q: usize, r: RoundingMode, c: &mut Consts| {
                let nt = ExactNum::mul_ref(&lam_w, tt, q, r).neg();
                ExactNum::mul_ref(&lam_w, &nt.exp(q, r, c), q, r)
            };
        let decided = super::a_priori_trajectory_unrepresentable(
            &tab,
            &mut pulse,
            [&t0, &y0],
            [&h0, &span],
            [&atol, &rtol],
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
                let nt = ExactNum::mul_ref(&lam, tt, q, r).neg();
                ExactNum::mul_ref(&lam, &nt.exp(q, r, c), q, r)
            },
            OdeSpan {
                t0: &zero,
                y0: &zero,
                t1: &t1,
            },
            OdeTol {
                atol: &atol,
                rtol: &rtol,
            },
            p,
            rm,
            &mut cc,
        );
        let ms = wall.elapsed().as_millis();
        let (ts, ys) = match got {
            Ok(v) => v,
            Err(e) => {
                return Err(format!(
                    "adaptive run must finish, not {e:?}; evals={evals} in {ms} ms"
                ))
            }
        };
        let yend = ys.get(ys.len() - 1).ok_or(crate::Error::InvalidArgument)?;
        let tend = ts.get(ts.len() - 1).ok_or(crate::Error::InvalidArgument)?;
        assert_eq!(tend.compare(&t1), Some(0), "must land on t1; evals={evals}");
        // Exact y(t) = 1 - exp(-λ t) → 1 at t = 2^{400}.
        let err = ExactNum::sub_ref(yend, &one, p, rm).abs();
        let band = one.ldexp(-12, p, RoundingMode::None);
        assert!(
            err.is_zero() || err.compare(&band) != Some(1),
            "y(t1) must be 1 within a generous multiple of atol; evals={evals}"
        );
        assert!(evals > 16, "the run must actually integrate; evals={evals}");
        assert!(
            evals < 100_000,
            "the tail must grow, not crawl; evals={evals}"
        );
        assert!(
            crate::test_time::within_ms(ms, 2_000),
            "the finishable run must be quick; took {ms} ms, evals={evals}"
        );
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
            OdeSpan {
                t0: &zero,
                y0: &one,
                t1: &one,
            },
            OdeTol {
                atol: &tol,
                rtol: &tol,
            },
            p,
            rm,
            &mut cc,
        )
        .map_err(|_| crate::Error::InvalidArgument)?;
        let yend = ys.get(ys.len() - 1).ok_or(crate::Error::InvalidArgument)?;
        let tend = ts.get(ts.len() - 1).ok_or(crate::Error::InvalidArgument)?;
        assert_eq!(tend.compare(&one), Some(0));
        let expect = ExactNum::add_ref(&one, &one.sin(p, rm, &mut cc), p, rm);
        let err = ExactNum::sub_ref(yend, &expect, p, rm).abs();
        // The accepted local scale is atol + rtol |y| ≈ 2^{-47}. Allow a
        // generous global multiple of that scale, still far below 2^{-40}.
        let band = ExactNum::from_u8(1, p).ldexp(-40, p, RoundingMode::None);
        assert!(
            err.is_zero() || err.compare(&band) != Some(1),
            "reachable tight tol must return 1+sin(1); evals={evals_r}"
        );
        assert!(
            evals_r > 16,
            "a finishable run actually integrates; evals={evals_r}"
        );

        // Same reachable bar on a different f, so the bound is class-wide.
        let mut evals_e = 0u32;
        let (_, ye) = rk45_adaptive(
            |_t, y, _q, _r, _cc| {
                evals_e += 1;
                y.neg_via_mul()
            },
            OdeSpan {
                t0: &zero,
                y0: &one,
                t1: &one,
            },
            OdeTol {
                atol: &tol,
                rtol: &tol,
            },
            p,
            rm,
            &mut cc,
        )
        .map_err(|_| crate::Error::InvalidArgument)?;
        let yende = ye.get(ye.len() - 1).ok_or(crate::Error::InvalidArgument)?;
        let einv = one.neg().exp(p, rm, &mut cc);
        let erre = ExactNum::sub_ref(yende, &einv, p, rm).abs();
        assert!(
            erre.is_zero() || erre.compare(&band) != Some(1),
            "y'=-y at the same tight tol must return e^{{-1}}; evals={evals_e}"
        );
        Ok(())
    }

    fn cos_closed_form_bits(
        p: usize,
        tol_exp: i64,
        min_bits: i32,
        max_ms: u128,
        label: &str,
    ) -> Result<(), String> {
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);
        let tol = ExactNum::from_u8(1, p).ldexp(tol_exp, p, RoundingMode::None);
        let t0 = std::time::Instant::now();
        let mut evals = 0u32;
        let got = rk45_adaptive(
            |t, _y, q, r, cc| {
                evals = evals.saturating_add(1);
                t.cos(q, r, cc)
            },
            OdeSpan {
                t0: &zero,
                y0: &one,
                t1: &one,
            },
            OdeTol {
                atol: &tol,
                rtol: &tol,
            },
            p,
            rm,
            &mut cc,
        );
        let ms = t0.elapsed().as_millis();
        let (ts, ys) = match got {
            Ok(v) => v,
            Err(e) => return Err(format!("{label} must integrate, not {e:?} after {ms} ms")),
        };
        let yend = ys.get(ys.len() - 1).ok_or(crate::Error::InvalidArgument)?;
        let tend = ts.get(ts.len() - 1).ok_or(crate::Error::InvalidArgument)?;
        assert_eq!(
            tend.compare(&one),
            Some(0),
            "{label} must land on t=1; {ms} ms"
        );
        let expect = ExactNum::add_ref(&one, &one.sin(p, rm, &mut cc), p, rm);
        let err = ExactNum::sub_ref(yend, &expect, p, rm).abs();
        let bits = if err.is_zero() {
            10_000
        } else {
            (expect.exponent().unwrap_or(0) - err.exponent().unwrap_or(0)) as i32
        };
        eprintln!("{label}: {ms} ms, {bits} bits, {evals} f-evals");
        let _ = std::io::Write::flush(&mut std::io::stderr());
        assert!(
            bits >= min_bits || err.is_zero(),
            "{label}: y(1) vs 1+sin(1) only {bits} bits after {ms} ms"
        );
        assert!(
            crate::test_time::limit_holds(ms < max_ms),
            "{label} must finish under {max_ms} ms; took {ms} ms"
        );
        let _ = ts;
        Ok(())
    }

    /// Same closed form as Harper ZF-54, at a width where a BS tableau
    /// is cheap. atol = 2^{-96} is already out of reach for Dormand–Prince
    /// on [0, 1] (order-5 h ≲ 2^{-19}, tens of millions of steps).
    #[test]
    fn zf_54_bs_cos_matches_closed_form_p256() -> Result<(), String> {
        cos_closed_form_bits(256, -96, 80, 10_000, "y'=cos(t) p=256 atol=2^{-96}")
    }

    /// Harper ZF-54 / official d0229592: y' = cos(t), y(0) = 1, [0, 1],
    /// atol = rtol = 2^{-512}, dest p = 1024. Closed form y(1) = 1 + sin(1),
    /// never a synthetic reference. Dormand–Prince cannot reach this
    /// tolerance in a feasible number of order-5 steps; the documented
    /// switch takes Bulirsch–Stoer and must return the closed form.
    #[test]
    fn zf_54_d0229592_matches_closed_form() -> Result<(), String> {
        cos_closed_form_bits(
            1024,
            -512,
            500,
            90_000,
            "d0229592 y'=cos(t) p=1024 atol=2^{-512}",
        )
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
            OdeSpan {
                t0: &zero,
                y0: &one,
                t1: &one,
            },
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

    /// A zero span (`t1 = t0`) is not an empty successful trajectory:
    /// it is `InvalidArgument`, the same as a reversed span.
    #[test]
    fn zf_54_zero_span_is_invalid() -> Result<(), crate::Error> {
        let p = 64;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);
        let tol = ExactNum::from_u8(1, p).ldexp(-16, p, RoundingMode::None);
        let mut evals = 0u32;
        let got = rk45_adaptive(
            |_t, _y, _q, _r, _cc| {
                evals += 1;
                ExactNum::new(p)
            },
            OdeSpan {
                t0: &zero,
                y0: &one,
                t1: &zero,
            },
            OdeTol {
                atol: &tol,
                rtol: &tol,
            },
            p,
            rm,
            &mut cc,
        );
        assert_eq!(got.err(), Some(OdeError::InvalidArgument));
        assert_eq!(evals, 0, "zero span refuses before any evaluation");
        Ok(())
    }

    /// A reversed span (`t1 < t0`) is `InvalidArgument`, not a backward
    /// integration and not `ToleranceUnreachable`.
    #[test]
    fn zf_54_reversed_span_is_invalid() -> Result<(), crate::Error> {
        let p = 64;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);
        let two = ExactNum::from_u8(2, p);
        let tol = ExactNum::from_u8(1, p).ldexp(-16, p, RoundingMode::None);
        let mut evals = 0u32;
        let got = rk45_adaptive(
            |_t, _y, _q, _r, _cc| {
                evals += 1;
                ExactNum::new(p)
            },
            OdeSpan {
                t0: &two,
                y0: &one,
                t1: &zero,
            },
            OdeTol {
                atol: &tol,
                rtol: &tol,
            },
            p,
            rm,
            &mut cc,
        );
        assert_eq!(got.err(), Some(OdeError::InvalidArgument));
        assert_eq!(evals, 0, "reversed span refuses before any evaluation");
        Ok(())
    }

    /// A tolerance at the working-precision floor (`atol = 2^{-(wrk-9)}`,
    /// just inside [`wrk_supports_tol`]) on `y' = cos(t)` must finish,
    /// either on Dormand–Prince or by the documented Bulirsch–Stoer
    /// switch, and match the closed form.
    #[test]
    fn zf_54_tol_at_precision_floor_matches_closed_form() -> Result<(), String> {
        let p = 64;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);
        let wrk = super::work_p(p);
        // wrk_supports_tol requires scale > 2^{-(wrk-8)}. Use 2^{-(wrk-9)}.
        let exp = -i64::try_from(wrk.saturating_sub(9)).must("i64::try_from")?;
        let atol = ExactNum::from_u8(1, p).ldexp(exp, p, RoundingMode::None);
        assert!(
            super::wrk_supports_tol(&atol, &ExactNum::new(p), &one, wrk),
            "this atol must sit just inside the wrk floor"
        );
        let wall = std::time::Instant::now();
        let got = rk45_adaptive(
            |t, _y, q, r, cc| t.cos(q, r, cc),
            OdeSpan {
                t0: &zero,
                y0: &one,
                t1: &one,
            },
            OdeTol {
                atol: &atol,
                rtol: &ExactNum::new(p),
            },
            p,
            rm,
            &mut cc,
        );
        let ms = wall.elapsed().as_millis();
        let (ts, ys) = match got {
            Ok(v) => v,
            Err(e) => {
                return Err(format!(
                    "tol at wrk floor must finish, not {e:?} in {ms} ms"
                ))
            }
        };
        let yend = ys.get(ys.len() - 1).ok_or(crate::Error::InvalidArgument)?;
        let tend = ts.get(ts.len() - 1).ok_or(crate::Error::InvalidArgument)?;
        assert_eq!(tend.compare(&one), Some(0));
        let expect = ExactNum::add_ref(&one, &one.sin(p, rm, &mut cc), p, rm);
        let err = ExactNum::sub_ref(yend, &expect, p, rm).abs();
        let band = one.ldexp(-20, p, RoundingMode::None);
        assert!(
            err.is_zero() || err.compare(&band) != Some(1),
            "y(1) vs 1+sin(1) at the wrk floor; {ms} ms"
        );
        assert!(
            crate::test_time::within_ms(ms, 5_000),
            "floor-tol run must be quick; took {ms} ms"
        );
        let _ = ts;
        Ok(())
    }

    /// A moderately stiff linear sink `y' = −λ y` with `λ = 40` on
    /// `[0, 1]` (stability bound `h ≲ 3/λ ≈ 0.075` for explicit RK)
    /// is finishable on Dormand–Prince and must match `e^{{-40}}`.
    /// This is stiffer than `y' = −y` but not a stability-limited
    /// crawl; the error-model switch must not fire.
    #[test]
    fn zf_54_stiffish_transient_matches_closed_form() -> Result<(), String> {
        let p = 64;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);
        let lam = ExactNum::from_u32(40, p);
        let atol = ExactNum::from_u8(1, p).ldexp(-24, p, RoundingMode::None);
        let mut evals = 0u32;
        let wall = std::time::Instant::now();
        let got = rk45_adaptive(
            |_t, y, q, r, _cc| {
                evals += 1;
                ExactNum::mul_ref(y, &lam, q, r).neg()
            },
            OdeSpan {
                t0: &zero,
                y0: &one,
                t1: &one,
            },
            OdeTol {
                atol: &atol,
                rtol: &atol,
            },
            p,
            rm,
            &mut cc,
        );
        let ms = wall.elapsed().as_millis();
        let (ts, ys) = match got {
            Ok(v) => v,
            Err(e) => {
                return Err(format!(
                    "stiff-ish y'=-40y must finish, not {e:?}; evals={evals} in {ms} ms"
                ))
            }
        };
        let yend = ys.get(ys.len() - 1).ok_or(crate::Error::InvalidArgument)?;
        let tend = ts.get(ts.len() - 1).ok_or(crate::Error::InvalidArgument)?;
        assert_eq!(tend.compare(&one), Some(0));
        let expect = lam.neg().exp(p, rm, &mut cc);
        let err = ExactNum::sub_ref(yend, &expect, p, rm).abs();
        let band = one.ldexp(-16, p, RoundingMode::None);
        assert!(
            err.is_zero() || err.compare(&band) != Some(1),
            "y(1) vs e^{{-40}}; evals={evals} in {ms} ms"
        );
        assert!(evals > 16, "must actually integrate; evals={evals}");
        assert!(evals < 10_000, "must not crawl; evals={evals}");
        assert!(
            crate::test_time::within_ms(ms, 2_000),
            "stiff-ish run must be quick; took {ms} ms, evals={evals}"
        );
        let _ = ts;
        Ok(())
    }

    /// One-off cost of [`ODE_DP_MAX_CONST_C_STEPS`] Dormand–Prince
    /// accepted steps versus one Bulirsch–Stoer sweep at the same
    /// working precision. Run with `--ignored --release`.
    #[test]
    #[ignore = "measures the 2^14 switch bound; run with --ignored --release"]
    fn zf_54_measure_const_c_bound_vs_bs() -> Result<(), crate::Error> {
        let rm = RoundingMode::None;
        let n_dp = super::ODE_DP_MAX_CONST_C_STEPS;
        // wrk = dest p + 64, the values used by the probe, the p=256
        // closed-form test, and Harper ZF-54.
        let cases: [(usize, i64); 3] = [(128, -16), (320, -96), (1088, -512)];
        for (wrk, atol_exp) in cases {
            let mut cc = Consts::new()?;
            let zero = ExactNum::new(wrk);
            let one = ExactNum::from_u8(1, wrk);
            let atol = ExactNum::from_u8(1, wrk).ldexp(atol_exp, wrk, rm);
            let rtol = ExactNum::new(wrk);
            let tab = super::Dp45::new(wrk, rm);
            let h = ExactNum::div_ref(&one, &ExactNum::from_u32(n_dp, wrk), wrk, rm);
            let mut t = ExactNum::new(wrk);
            let mut y = ExactNum::from_u8(1, wrk);
            let wall = std::time::Instant::now();
            let mut evals = 0u32;
            for i in 0..n_dp {
                if i == 0 || (i + 1) % 2048 == 0 {
                    let done = i;
                    // Duration arithmetic (no hardware float tokens: scripts/ci_no_intrinsics.sh).
                    let spent = wall.elapsed();
                    let eta = if done == 0 {
                        std::time::Duration::ZERO
                    } else {
                        spent * (n_dp - done) / done
                    };
                    let (elapsed, eta) = (spent.as_secs_f64(), eta.as_secs_f64());
                    eprintln!(
                        "BOUND wrk={wrk} phase=dp45 step={}/{n_dp} evals={evals} elapsed={elapsed:.3}s eta={eta:.3}s",
                        done
                    );
                }
                let mut f =
                    |tt: &ExactNum, _y: &ExactNum, q: usize, r: RoundingMode, c: &mut Consts| {
                        evals += 1;
                        tt.cos(q, r, c)
                    };
                let (y5, _y4) = super::dp45_step(&tab, &mut f, &t, &y, &h, (wrk, rm), &mut cc)
                    .ok_or(crate::Error::InvalidArgument)?;
                t = ExactNum::add_ref(&t, &h, wrk, rm);
                y = y5;
            }
            let dp_ms = wall.elapsed().as_millis();
            eprintln!("BOUND wrk={wrk} phase=dp45_done steps={n_dp} evals={evals} ms={dp_ms}");

            // One BS sweep over the same unit span at the same atol.
            let mut evals_bs = 0u32;
            let wall_bs = std::time::Instant::now();
            let mut f =
                |tt: &ExactNum, _y: &ExactNum, q: usize, r: RoundingMode, c: &mut Consts| {
                    evals_bs += 1;
                    tt.cos(q, r, c)
                };
            let got = super::bs_try_step(&mut f, &zero, &one, &one, [&atol, &rtol], wrk, &mut cc);
            let try_ms = wall_bs.elapsed().as_millis();
            let try_ok = matches!(got, Ok(Some(_)));
            eprintln!("BOUND wrk={wrk} phase=bs_try_step ok={try_ok} evals={evals_bs} ms={try_ms}");

            // Full adaptive BS on [0, 1] (may shrink H and retry).
            let mut evals_full = 0u32;
            let wall_full = std::time::Instant::now();
            let mut f =
                |tt: &ExactNum, _y: &ExactNum, q: usize, r: RoundingMode, c: &mut Consts| {
                    evals_full += 1;
                    tt.cos(q, r, c)
                };
            let full = super::bs_at_prec(&mut f, &zero, &one, &one, [&atol, &rtol], wrk, &mut cc);
            let full_ms = wall_full.elapsed().as_millis();
            let full_ok = full.is_ok();
            eprintln!(
                "BOUND wrk={wrk} phase=bs_at_prec ok={full_ok} evals={evals_full} ms={full_ms}"
            );
            let _ = (dp_ms, try_ms, full_ms, y);
        }
        Ok(())
    }
}
