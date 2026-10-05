//! Gauss–Kronrod quadrature and Wynn’s ε table on [`ExactNum`].
//!
//! The (7, 15) pair is the fixed rule QUADPACK uses inside `QAG` / `QAGS`:
//! a 7-point Gauss–Legendre estimate embedded in a 15-point Kronrod estimate,
//! with absolute error `|K − G|`. Endpoint singularities are accelerated by
//! Wynn’s ε-algorithm on the sequence of partial sums, the same device
//! Piessens, de Doncker-Kapenga, Überhuber, and Kahaner describe for `QAGS`
//! (*QUADPACK*, Springer, 1983; the SLATEC port of that design is public
//! domain). This module recomputes the rule in SoftFloat. It does not paste
//! the Fortran.
//!
//! Nodes and weights come from Laurie’s Jacobi–Kronrod recurrence (D. P. Laurie,
//! *Math. Comp.* 66 (1997), 1133–1145) for the Legendre weight on `[-1, 1]`,
//! followed by Golub–Welsch weights (the squared first component of the
//! normalized eigenvector). The matrix is symmetric, so the diagonal is zero
//! and only the squared sub-diagonals `β_j = j² / (4j² − 1)` are stored.
//! Everything is evaluated at `p + WORD_BIT_SIZE` bits and then rounded to `p`.
//!
//! # Not in this slice
//!
//! Infinite intervals (`QAGI`), oscillatory weights (`QAWO` / `QAWF`),
//! algebraic endpoint weights (`QAWS`), Cauchy principal values (`QAWC`),
//! Kronrod pairs other than (7, 15), and the full `QAGS` heap with its
//! four-way split and roundoff flags.

use crate::defs::WORD_BIT_SIZE;
use crate::Consts;
use crate::ExactNum;
use crate::RoundingMode;
use alloc::vec;
use alloc::vec::Vec;
use core::cmp::Ordering;

/// Gauss order of the supported Kronrod pair. Only this order is built.
pub const GK_GAUSS_ORDER: usize = 7;

/// Number of Kronrod nodes in the [`GK_GAUSS_ORDER`] pair (`2n + 1`).
pub const GK_KRONROD_ORDER: usize = 2 * GK_GAUSS_ORDER + 1;

/// Largest `max_subintervals` accepted by [`integrate_adaptive_gk`].
pub const GK_SUBINTERVAL_MAX: usize = 32;

/// Destination precision below this cannot host a guarded (7, 15) rule.
const GK_MIN_PREC: usize = 32;

/// Newton corrections after the sign-change bracket is bisected.
const GK_NEWTON_MAX: usize = 32;

/// Bisection steps used to separate a characteristic-polynomial root.
const GK_BISECT_STEPS: usize = 48;

/// Why a Kronrod call stopped.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum QuadratureError {
    /// Non-finite bounds, an empty interval, a bad tolerance or precision,
    /// an unsupported setting, or a non-finite integrand sample.
    InvalidArgument,
    /// The panel budget or the ε table ran out before the tolerance was met,
    /// or the (7, 15) rule could not be built at this precision.
    NotConverged,
}

impl core::fmt::Display for QuadratureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let repr = match self {
            Self::InvalidArgument => "invalid argument",
            Self::NotConverged => "quadrature did not converge",
        };
        f.write_str(repr)
    }
}

#[cfg(feature = "std")]
impl std::error::Error for QuadratureError {}

/// One (7, 15) panel, or the value accepted by [`integrate_adaptive_gk`].
#[derive(Clone, Debug)]
pub struct GkResult {
    /// Kronrod sum `K`. When [`integrate_adaptive_gk`] accepts an ε value,
    /// this field is that extrapolation rather than the raw panel sum.
    pub estimate_kronrod: ExactNum,
    /// Embedded Gauss sum `G` on the same nodes (not extrapolated).
    pub estimate_gauss: ExactNum,
    /// Absolute error estimate. A single panel reports `|K − G|` after the
    /// factor `(b − a) / 2`. An accepted ε value reports the difference of
    /// the last two even-order ε entries.
    pub abs_err_est: ExactNum,
}

struct Panel {
    a: ExactNum,
    b: ExactNum,
    kronrod: ExactNum,
    gauss: ExactNum,
    err: ExactNum,
}

struct WynnValue {
    value: ExactNum,
    /// Even column index. `0` means the table never left the raw sequence.
    order: usize,
}

fn work_p(p: usize) -> usize {
    p.saturating_add(WORD_BIT_SIZE)
}

fn finite(x: &ExactNum) -> bool {
    !x.is_nan() && !x.is_inf()
}

fn tiny(p: usize, extra: i32) -> ExactNum {
    let e = (p as i32).saturating_sub(extra).saturating_neg();
    ExactNum::from_u8(1, p).ldexp(e, p, RoundingMode::ToEven)
}

fn leq(x: &ExactNum, bound: &ExactNum) -> bool {
    matches!(x.cmp(bound), Some(-1) | Some(0))
}

fn round_out(mut x: ExactNum, p: usize, rm: RoundingMode) -> ExactNum {
    let _ = x.set_precision(p, rm);
    x
}

fn ord_num(a: &ExactNum, b: &ExactNum) -> Ordering {
    match a.cmp(b) {
        Some(c) if c < 0 => Ordering::Less,
        Some(c) if c > 0 => Ordering::Greater,
        _ => Ordering::Equal,
    }
}

/// `β_j = j² / (4j² − 1)` for `j = 1..=count` (Laurie’s squared sub-diagonals).
fn legendre_beta_sq(count: usize, p: usize) -> Option<Vec<ExactNum>> {
    let rm = RoundingMode::None;
    let one = ExactNum::from_u8(1, p);
    let four = ExactNum::from_u8(4, p);
    let mut b = Vec::with_capacity(count);
    for j in 1..=count {
        let jf = ExactNum::from_u32(j as u32, p);
        let js = jf.mul(&jf, p, rm);
        let den = four.mul(&js, p, rm).sub(&one, p, rm);
        let beta = js.div(&den, p, rm);
        if !finite(&beta) || !beta.is_positive() {
            return None;
        }
        b.push(beta);
    }
    Some(b)
}

/// Laurie’s eastward / southward phases for a zero diagonal.
///
/// `b[j]` is `β_{j+1}`. Entries `j >= floor((3n+1)/2)` start at zero and are
/// filled by the southward phase. The algebra is Laurie’s appendix specialized
/// to `a_k = 0` (symmetric Legendre weight): the diagonal corrections drop out
/// and only the `β` updates remain.
fn laurie_legendre_beta(n: usize, p: usize) -> Option<Vec<ExactNum>> {
    let known = (3 * n + 1) / 2;
    let mut b = legendre_beta_sq(known, p)?;
    b.resize(2 * n, ExactNum::new(p));
    let slen = n / 2 + 2;
    let mut s = vec![ExactNum::new(p); slen];
    let mut t = vec![ExactNum::new(p); slen];
    let rm = RoundingMode::None;
    t[1] = b[n].clone();

    for m in 0..n.saturating_sub(1) {
        let mut u = ExactNum::new(p);
        let k_hi = (m + 1) / 2;
        for k in (0..=k_hi).rev() {
            let mut term = b[k + n].mul(&s[k], p, rm);
            if m > k {
                let bl = b[m - k - 1].mul(&s[k + 1], p, rm);
                term = term.sub(&bl, p, rm);
            }
            if !finite(&term) {
                return None;
            }
            u = u.add(&term, p, rm);
            s[k + 1] = u.clone();
        }
        core::mem::swap(&mut s, &mut t);
    }

    for j in (0..=n / 2).rev() {
        s[j + 1] = s[j].clone();
    }

    let m_stop = 2 * n - 2;
    for m in (n - 1)..m_stop {
        let mut u = ExactNum::new(p);
        let k0 = m + 1 - n;
        let k1 = (m - 1) / 2;
        if k0 <= k1 {
            for k in k0..=k1 {
                let j = n as i32 - (m as i32 - k as i32) - 1;
                let sj = s_at(&s, j + 1)?;
                let sj1 = s_at(&s, j + 2)?;
                let bk = b[k + n].mul(&sj, p, rm);
                let bl = b[m - k - 1].mul(&sj1, p, rm);
                let piece = bk.sub(&bl, p, rm);
                if !finite(&piece) {
                    return None;
                }
                u = u.sub(&piece, p, rm);
                s_set(&mut s, j + 1, u.clone())?;
            }
        }
        let k = (m + 1) / 2;
        let j = n as i32 - (m as i32 - k as i32 + 2);
        if 2 * k != m {
            let num = s_at(&s, j + 1)?;
            let den = s_at(&s, j + 2)?;
            if den.is_zero() || !finite(&den) || !finite(&num) {
                return None;
            }
            let beta = num.div(&den, p, rm);
            if !finite(&beta) {
                return None;
            }
            b[k + n] = beta;
        }
        core::mem::swap(&mut s, &mut t);
    }

    if b.iter().any(|beta| !finite(beta) || !beta.is_positive()) {
        None
    } else {
        Some(b)
    }
}

fn s_at(s: &[ExactNum], idx: i32) -> Option<ExactNum> {
    if idx < 0 {
        return None;
    }
    s.get(idx as usize).cloned()
}

fn s_set(s: &mut [ExactNum], idx: i32, value: ExactNum) -> Option<()> {
    if idx < 0 {
        return None;
    }
    let slot = s.get_mut(idx as usize)?;
    *slot = value;
    Some(())
}

/// Characteristic polynomial of the hollow Jacobi matrix, and its derivative.
///
/// `p_0 = 1`, `p_1 = λ`, `p_{k+1} = λ p_k − β_k p_{k−1}`.
fn char_poly(lam: &ExactNum, beta: &[ExactNum], p: usize) -> Option<(ExactNum, ExactNum)> {
    let rm = RoundingMode::None;
    let mut p0 = ExactNum::from_u8(1, p);
    let mut d0 = ExactNum::new(p);
    let mut p1 = lam.clone();
    let mut d1 = ExactNum::from_u8(1, p);
    if !finite(&p1) {
        return None;
    }
    for b in beta {
        let p_new = lam.mul(&p1, p, rm).sub(&b.mul(&p0, p, rm), p, rm);
        let d_new = p1
            .add(&lam.mul(&d1, p, rm), p, rm)
            .sub(&b.mul(&d0, p, rm), p, rm);
        if !finite(&p_new) || !finite(&d_new) {
            return None;
        }
        p0 = p1;
        d0 = d1;
        p1 = p_new;
        d1 = d_new;
    }
    Some((p1, d1))
}

fn polish_root(
    beta: &[ExactNum],
    mut lo: ExactNum,
    mut hi: ExactNum,
    p: usize,
) -> Option<ExactNum> {
    let rm = RoundingMode::None;
    let two = ExactNum::from_u8(2, p);
    let (mut flo, _) = char_poly(&lo, beta, p)?;
    if flo.is_zero() {
        return Some(lo);
    }
    for _ in 0..GK_BISECT_STEPS {
        let mid = lo.add(&hi, p, rm).div(&two, p, rm);
        let (fm, _) = char_poly(&mid, beta, p)?;
        if fm.is_zero() {
            return Some(mid);
        }
        if flo.is_positive() == fm.is_positive() {
            lo = mid;
            flo = fm;
        } else {
            hi = mid;
        }
    }
    let mut x = lo.add(&hi, p, rm).div(&two, p, rm);
    let guard = tiny(p, 8);
    for _ in 0..GK_NEWTON_MAX {
        let (y, d) = char_poly(&x, beta, p)?;
        if d.is_zero() {
            break;
        }
        let step = y.div(&d, p, rm);
        if !finite(&step) {
            break;
        }
        let next = x.sub(&step, p, rm);
        if !finite(&next) {
            break;
        }
        x = next;
        if leq(&step.abs(), &guard) {
            break;
        }
    }
    if finite(&x) {
        Some(x)
    } else {
        None
    }
}

fn positive_roots(beta: &[ExactNum], want: usize, p: usize) -> Option<Vec<ExactNum>> {
    let mut steps: usize = (want.saturating_mul(8)).max(64);
    while steps <= 4096 {
        if let Some(roots) = isolate_roots(beta, want, steps, p) {
            return Some(roots);
        }
        steps = steps.saturating_mul(2);
    }
    None
}

fn isolate_roots(beta: &[ExactNum], want: usize, steps: usize, p: usize) -> Option<Vec<ExactNum>> {
    let rm = RoundingMode::None;
    let denom = ExactNum::from_u32(steps as u32, p);
    let mut prev_x = ExactNum::from_u32(1, p).div(&denom, p, rm);
    let (mut prev_y, _) = char_poly(&prev_x, beta, p)?;
    if !finite(&prev_y) {
        return None;
    }
    let mut brackets = Vec::new();
    for i in 2..=steps {
        let x = ExactNum::from_u32(i as u32, p).div(&denom, p, rm);
        let (y, _) = char_poly(&x, beta, p)?;
        if !finite(&y) {
            return None;
        }
        let sign_change =
            prev_y.is_zero() || y.is_zero() || prev_y.is_positive() != y.is_positive();
        if sign_change {
            brackets.push((prev_x.clone(), x.clone()));
        }
        if brackets.len() > want {
            return None;
        }
        prev_x = x;
        prev_y = y;
    }
    if brackets.len() != want {
        return None;
    }
    let one = ExactNum::from_u8(1, p);
    let mut roots = Vec::with_capacity(want);
    for (lo, hi) in brackets {
        let root = polish_root(beta, lo, hi, p)?;
        if !root.is_positive() || root.cmp(&one) != Some(-1) {
            return None;
        }
        roots.push(root);
    }
    roots.sort_by(ord_num);
    for w in roots.windows(2) {
        if w[0].cmp(&w[1]) != Some(-1) {
            return None;
        }
    }
    Some(roots)
}

/// Golub–Welsch weight `μ0 * v_0²` with `v` the normalized eigenvector.
fn jacobi_weight(
    lam: &ExactNum,
    beta_sq: &[ExactNum],
    mu0: &ExactNum,
    p: usize,
) -> Option<ExactNum> {
    let rm = RoundingMode::None;
    let m = beta_sq.len() + 1;
    if m == 1 {
        return Some(mu0.clone());
    }
    let mut off = Vec::with_capacity(beta_sq.len());
    for b in beta_sq {
        if !b.is_positive() {
            return None;
        }
        let s = b.sqrt(p, rm);
        if !finite(&s) || s.is_zero() {
            return None;
        }
        off.push(s);
    }
    let mut v_prev2 = ExactNum::from_u8(1, p);
    let mut sumsq = ExactNum::from_u8(1, p);
    let mut v_prev = lam.div(&off[0], p, rm);
    if !finite(&v_prev) {
        return None;
    }
    sumsq = sumsq.add(&v_prev.mul(&v_prev, p, rm), p, rm);
    for i in 2..m {
        let numer = lam
            .mul(&v_prev, p, rm)
            .sub(&off[i - 2].mul(&v_prev2, p, rm), p, rm);
        let v = numer.div(&off[i - 1], p, rm);
        if !finite(&v) {
            return None;
        }
        sumsq = sumsq.add(&v.mul(&v, p, rm), p, rm);
        v_prev2 = v_prev;
        v_prev = v;
    }
    if !finite(&sumsq) || !sumsq.is_positive() {
        return None;
    }
    let w = mu0.div(&sumsq, p, rm);
    if finite(&w) && w.is_positive() {
        Some(w)
    } else {
        None
    }
}

fn weight_sum_ok(weights: &[ExactNum], target: &ExactNum, p: usize) -> bool {
    let rm = RoundingMode::None;
    let mut acc = ExactNum::new(p);
    for w in weights {
        acc = acc.add(w, p, rm);
    }
    let err = acc.sub(target, p, rm).abs();
    // Working precision is `p`. Losing a quarter of those bits still leaves
    // the caller's destination (`p - WORD_BIT_SIZE`) intact.
    let slack_bits = (p as i32) / 4;
    leq(&err, &tiny(p, slack_bits))
}

/// (7, 15) nodes and both weight rows on `[-1, 1]`, still at precision `p`.
fn kronrod_legendre_at(
    n: usize,
    p: usize,
    _cc: &mut Consts,
) -> Option<(Vec<ExactNum>, Vec<ExactNum>, Vec<ExactNum>)> {
    if n != GK_GAUSS_ORDER || p < GK_MIN_PREC {
        return None;
    }
    let beta = laurie_legendre_beta(n, p)?;
    let positive = positive_roots(&beta, n, p)?;
    let two = ExactNum::from_u8(2, p);
    let gbeta = legendre_beta_sq(n - 1, p)?;
    let mut nodes = Vec::with_capacity(2 * n + 1);
    let mut w_kronrod = Vec::with_capacity(2 * n + 1);
    let mut w_gauss = vec![ExactNum::new(p); 2 * n + 1];
    let zero = ExactNum::new(p);
    let w0 = jacobi_weight(&zero, &beta, &two, p)?;
    for r in positive.iter().rev() {
        let w = jacobi_weight(r, &beta, &two, p)?;
        nodes.push(r.neg());
        w_kronrod.push(w);
    }
    nodes.push(zero);
    w_kronrod.push(w0);
    for r in &positive {
        let w = jacobi_weight(r, &beta, &two, p)?;
        nodes.push(r.clone());
        w_kronrod.push(w);
    }
    if nodes.len() != GK_KRONROD_ORDER {
        return None;
    }
    for i in (1..nodes.len()).step_by(2) {
        w_gauss[i] = jacobi_weight(&nodes[i], &gbeta, &two, p)?;
    }
    if !weight_sum_ok(&w_kronrod, &two, p) || !weight_sum_ok(&w_gauss, &two, p) {
        return None;
    }
    Some((nodes, w_gauss, w_kronrod))
}

/// (7, 15) Gauss–Kronrod nodes and weights on `[-1, 1]`.
///
/// `n_gauss` must be [`GK_GAUSS_ORDER`]. The three vectors have length
/// [`GK_KRONROD_ORDER`] and are ordered from the left endpoint toward the
/// right. `w_gauss[i]` is zero when `nodes[i]` is a Kronrod-only abscissa
/// (every even 0-based index). The embedded Gauss nodes are the odd indices.
///
/// Nodes and weights are computed at `p + WORD_BIT_SIZE` and rounded to `p`
/// with `rm`. Returns `None` when `n_gauss` is not 7, `p` is below 32, or the
/// Jacobi–Kronrod matrix does not yield a positive-weight rule.
pub fn kronrod_pair(
    n_gauss: usize,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<(Vec<ExactNum>, Vec<ExactNum>, Vec<ExactNum>)> {
    if n_gauss != GK_GAUSS_ORDER || p < GK_MIN_PREC {
        return None;
    }
    let (nodes, w_gauss, w_kronrod) = kronrod_legendre_at(n_gauss, work_p(p), cc)?;
    let round_row = |row: Vec<ExactNum>| -> Vec<ExactNum> {
        row.into_iter().map(|x| round_out(x, p, rm)).collect()
    };
    Some((round_row(nodes), round_row(w_gauss), round_row(w_kronrod)))
}

fn eval_panel<F>(
    f: &mut F,
    a: &ExactNum,
    b: &ExactNum,
    nodes: &[ExactNum],
    w_gauss: &[ExactNum],
    w_kronrod: &[ExactNum],
    p: usize,
    cc: &mut Consts,
) -> Result<(ExactNum, ExactNum), QuadratureError>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    let rm = RoundingMode::None;
    let two = ExactNum::from_u8(2, p);
    let half = b.sub(a, p, rm).div(&two, p, rm);
    let mid = a.add(b, p, rm).div(&two, p, rm);
    if !finite(&half) || !finite(&mid) || !half.is_positive() {
        return Err(QuadratureError::InvalidArgument);
    }
    let mut sk = ExactNum::new(p);
    let mut sg = ExactNum::new(p);
    for i in 0..nodes.len() {
        let x = mid.add(&half.mul(&nodes[i], p, rm), p, rm);
        let y = f(&x, p, rm, cc);
        if !finite(&y) {
            return Err(QuadratureError::InvalidArgument);
        }
        sk = sk.add(&w_kronrod[i].mul(&y, p, rm), p, rm);
        if !w_gauss[i].is_zero() {
            sg = sg.add(&w_gauss[i].mul(&y, p, rm), p, rm);
        }
        if !finite(&sk) || !finite(&sg) {
            return Err(QuadratureError::NotConverged);
        }
    }
    Ok((half.mul(&sk, p, rm), half.mul(&sg, p, rm)))
}

fn interval_ok(a: &ExactNum, b: &ExactNum) -> bool {
    finite(a) && finite(b) && a.cmp(b) == Some(-1)
}

/// Integrate `f` on the finite interval `[a, b]` with one (7, 15) panel.
///
/// # Precision
///
/// `p` is the destination precision in bits. The rule and the weighted sum
/// run at `p + WORD_BIT_SIZE`, then `estimate_kronrod`, `estimate_gauss`, and
/// `abs_err_est = |K − G|` are rounded with `rm`. There is no tolerance
/// argument: a single panel always returns if every sample is finite.
///
/// # Errors
///
/// [`QuadratureError::InvalidArgument`] if `p < 32`, `a` or `b` is non-finite,
/// `a >= b`, or `f` returns a non-finite value.
/// [`QuadratureError::NotConverged`] if the (7, 15) rule cannot be built.
pub fn gauss_kronrod_interval<F>(
    mut f: F,
    a: &ExactNum,
    b: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Result<GkResult, QuadratureError>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    if p < GK_MIN_PREC || !interval_ok(a, b) {
        return Err(QuadratureError::InvalidArgument);
    }
    let wrk = work_p(p);
    let (nodes, w_gauss, w_kronrod) =
        kronrod_legendre_at(GK_GAUSS_ORDER, wrk, cc).ok_or(QuadratureError::NotConverged)?;
    let mut aa = a.clone();
    let mut bb = b.clone();
    let _ = aa.set_precision(wrk, RoundingMode::None);
    let _ = bb.set_precision(wrk, RoundingMode::None);
    let (k, g) = eval_panel(&mut f, &aa, &bb, &nodes, &w_gauss, &w_kronrod, wrk, cc)?;
    let err = k.sub(&g, wrk, RoundingMode::None).abs();
    if !finite(&k) || !finite(&g) || !finite(&err) {
        return Err(QuadratureError::NotConverged);
    }
    Ok(GkResult {
        estimate_kronrod: round_out(k, p, rm),
        estimate_gauss: round_out(g, p, rm),
        abs_err_est: round_out(err, p, rm),
    })
}

fn wynn_table(seq: &[ExactNum], p: usize) -> Option<WynnValue> {
    if seq.is_empty() || seq.iter().any(|s| !finite(s)) {
        return None;
    }
    let rm = RoundingMode::None;
    let n = seq.len();
    let mut cols: Vec<Vec<ExactNum>> = Vec::with_capacity(n);
    cols.push(seq.to_vec());
    for k in 0..n.saturating_sub(1) {
        let mut next = Vec::with_capacity(cols[k].len() - 1);
        let mut broke = false;
        for i in 0..cols[k].len() - 1 {
            let diff = cols[k][i + 1].sub(&cols[k][i], p, rm);
            if !finite(&diff) || diff.is_zero() {
                broke = true;
                break;
            }
            let inv = ExactNum::from_u8(1, p).div(&diff, p, rm);
            if !finite(&inv) {
                broke = true;
                break;
            }
            let base = if k == 0 { ExactNum::new(p) } else { cols[k - 1][i + 1].clone() };
            let val = base.add(&inv, p, rm);
            if !finite(&val) {
                broke = true;
                break;
            }
            next.push(val);
        }
        if broke || next.len() != cols[k].len() - 1 {
            break;
        }
        cols.push(next);
    }
    let mut order = cols.len() - 1;
    if order % 2 == 1 {
        order -= 1;
    }
    let col = cols.get(order)?;
    let value = col.last()?.clone();
    Some(WynnValue { value, order })
}

/// Wynn ε-algorithm on `seq`.
///
/// The tableau is `ε_{-1}^{(n)} = 0`, `ε_0^{(n)} = seq[n]`, and
/// `ε_{k+1}^{(n)} = ε_{k-1}^{(n+1)} + 1 / (ε_k^{(n+1)} − ε_k^{(n)})`.
/// The returned value is the last entry of the highest completed even
/// column (the Shanks entry that consumes the newest samples).
///
/// A zero divisor stops the table. The last completed even column is kept,
/// so a pure geometric series still returns its limit after the column that
/// found it. `None` if `seq` is empty, any entry is non-finite, or `p < 32`.
///
/// The table is built at `p + WORD_BIT_SIZE` and rounded to `p` with `rm`.
pub fn wynn_epsilon(seq: &[ExactNum], p: usize, rm: RoundingMode) -> Option<ExactNum> {
    if p < GK_MIN_PREC {
        return None;
    }
    let wrk = work_p(p);
    let mut lifted = Vec::with_capacity(seq.len());
    for s in seq {
        if !finite(s) {
            return None;
        }
        let mut y = s.clone();
        let _ = y.set_precision(wrk, RoundingMode::None);
        lifted.push(y);
    }
    let value = wynn_table(&lifted, wrk)?.value;
    Some(round_out(value, p, rm))
}

fn abs_tolerance(tol_rel: &ExactNum, estimate: &ExactNum, p: usize) -> ExactNum {
    let rm = RoundingMode::None;
    let floor = tiny(p, 1);
    let one = ExactNum::from_u8(1, p);
    let mut scale = estimate.abs();
    if leq(&scale, &one) {
        scale = one;
    }
    let rel = tol_rel.mul(&scale, p, rm);
    if !finite(&rel) || leq(&rel, &floor) {
        floor
    } else {
        rel
    }
}

fn panel_totals(panels: &[Panel], p: usize) -> (ExactNum, ExactNum, ExactNum) {
    let rm = RoundingMode::None;
    let mut k = ExactNum::new(p);
    let mut g = ExactNum::new(p);
    let mut e = ExactNum::new(p);
    for panel in panels {
        k = k.add(&panel.kronrod, p, rm);
        g = g.add(&panel.gauss, p, rm);
        e = e.add(&panel.err, p, rm);
    }
    (k, g, e)
}

/// Adaptive (7, 15) integration of `f` on `[a, b]`, with Wynn ε.
///
/// The panel with the largest `|K − G|` is bisected until the sum of those
/// estimates meets the tolerance, or until two successive even-order ε
/// values on the sequence of global Kronrod sums agree to the tolerance.
/// This is the QUADPACK error idea (local `|K − G|`, global ε) without the
/// `QAGS` heap, four-way split, or roundoff flags.
///
/// # Arguments
///
/// * `tol_rel` — relative tolerance. The absolute stop is
///   `max(tol_rel * max(1, |estimate|), 2^(1−p))`, where `2^(1−p)` is the
///   unit in the last place of a `p`-bit significand. `tol_rel = 0` asks
///   for that floor only.
/// * `max_subintervals` — panel cap, from 1 through [`GK_SUBINTERVAL_MAX`].
///   Each panel evaluates `f` at [`GK_KRONROD_ORDER`] (15) nodes. A call
///   performs at most `15 * (2 * max_subintervals − 1)` evaluations, because
///   each bisection discards one panel and samples two new ones.
/// * `p` — destination precision in bits (at least 32). Nodes, weights, and
///   panel sums use `p + WORD_BIT_SIZE` guard bits.
/// * `rm` — rounding mode for the three returned fields.
///
/// # Errors
///
/// [`QuadratureError::InvalidArgument`] if bounds are non-finite or not
/// `a < b`, `tol_rel` is non-finite or negative, `p < 32`,
/// `max_subintervals` is outside `1..=GK_SUBINTERVAL_MAX`, or any integrand
/// sample is non-finite.
/// [`QuadratureError::NotConverged`] if the budget is spent or a panel can
/// no longer be bisected at the working precision.
pub fn integrate_adaptive_gk<F>(
    mut f: F,
    a: &ExactNum,
    b: &ExactNum,
    tol_rel: &ExactNum,
    max_subintervals: usize,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Result<GkResult, QuadratureError>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    if p < GK_MIN_PREC
        || max_subintervals == 0
        || max_subintervals > GK_SUBINTERVAL_MAX
        || !interval_ok(a, b)
        || !finite(tol_rel)
        || tol_rel.is_negative()
    {
        return Err(QuadratureError::InvalidArgument);
    }
    let wrk = work_p(p);
    let (nodes, w_gauss, w_kronrod) =
        kronrod_legendre_at(GK_GAUSS_ORDER, wrk, cc).ok_or(QuadratureError::NotConverged)?;
    let mut aa = a.clone();
    let mut bb = b.clone();
    let _ = aa.set_precision(wrk, RoundingMode::None);
    let _ = bb.set_precision(wrk, RoundingMode::None);
    let mut tol_work = tol_rel.clone();
    let _ = tol_work.set_precision(wrk, RoundingMode::None);
    let (k0, g0) = eval_panel(&mut f, &aa, &bb, &nodes, &w_gauss, &w_kronrod, wrk, cc)?;
    let err0 = k0.sub(&g0, wrk, RoundingMode::None).abs();
    let mut panels = vec![Panel {
        a: aa,
        b: bb,
        kronrod: k0,
        gauss: g0,
        err: err0,
    }];
    let mut seq: Vec<ExactNum> = Vec::new();
    let mut prev_eps: Option<ExactNum> = None;

    loop {
        let (total_k, total_g, total_err) = panel_totals(&panels, wrk);
        if !finite(&total_k) || !finite(&total_g) || !finite(&total_err) {
            return Err(QuadratureError::NotConverged);
        }
        seq.push(total_k.clone());
        let tol = abs_tolerance(&tol_work, &total_k, wrk);
        if leq(&total_err, &tol) {
            return Ok(GkResult {
                estimate_kronrod: round_out(total_k, p, rm),
                estimate_gauss: round_out(total_g, p, rm),
                abs_err_est: round_out(total_err, p, rm),
            });
        }
        if seq.len() >= 3 {
            if let Some(ex) = wynn_table(&seq, wrk) {
                if ex.order >= 2 {
                    if let Some(prev) = &prev_eps {
                        let de = ex.value.sub(prev, wrk, RoundingMode::None).abs();
                        if finite(&de) && leq(&de, &tol) {
                            return Ok(GkResult {
                                estimate_kronrod: round_out(ex.value, p, rm),
                                estimate_gauss: round_out(total_g, p, rm),
                                abs_err_est: round_out(de, p, rm),
                            });
                        }
                    }
                    prev_eps = Some(ex.value);
                }
            }
        }
        if panels.len() >= max_subintervals {
            return Err(QuadratureError::NotConverged);
        }
        let mut worst = 0usize;
        for i in 1..panels.len() {
            if panels[i].err.cmp(&panels[worst].err) == Some(1) {
                worst = i;
            }
        }
        let left_a = panels[worst].a.clone();
        let right_b = panels[worst].b.clone();
        let two = ExactNum::from_u8(2, wrk);
        let mid = left_a
            .add(&right_b, wrk, RoundingMode::None)
            .div(&two, wrk, RoundingMode::None);
        if mid.cmp(&left_a) != Some(1) || mid.cmp(&right_b) != Some(-1) {
            return Err(QuadratureError::NotConverged);
        }
        let (k1, g1) = eval_panel(&mut f, &left_a, &mid, &nodes, &w_gauss, &w_kronrod, wrk, cc)?;
        let (k2, g2) = eval_panel(
            &mut f, &mid, &right_b, &nodes, &w_gauss, &w_kronrod, wrk, cc,
        )?;
        let e1 = k1.sub(&g1, wrk, RoundingMode::None).abs();
        let e2 = k2.sub(&g2, wrk, RoundingMode::None).abs();
        panels[worst] = Panel {
            a: left_a,
            b: mid.clone(),
            kronrod: k1,
            gauss: g1,
            err: e1,
        };
        panels.push(Panel {
            a: mid,
            b: right_b,
            kronrod: k2,
            gauss: g2,
            err: e2,
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gold_p() -> (usize, RoundingMode) {
        (128, RoundingMode::ToEven)
    }

    fn bits(got: &ExactNum, want: &ExactNum, p: usize) -> i32 {
        if got.is_nan() {
            return i32::MIN;
        }
        let d = got.sub(want, 2 * p, RoundingMode::None);
        if d.is_zero() {
            return i32::MAX;
        }
        want.exponent().unwrap_or(0) - d.exponent().unwrap_or(0)
    }

    #[test]
    fn gk715_x22_exact_and_weights_symmetric() {
        let (p, rm) = gold_p();
        let mut cc = Consts::new().expect("consts");
        let (nodes, w_g, w_k) = kronrod_pair(GK_GAUSS_ORDER, p, rm, &mut cc).expect("pair");
        assert_eq!(nodes.len(), GK_KRONROD_ORDER);
        assert_eq!(w_g.len(), GK_KRONROD_ORDER);
        assert_eq!(w_k.len(), GK_KRONROD_ORDER);
        assert!(nodes[GK_GAUSS_ORDER].is_zero());
        for i in 0..GK_GAUSS_ORDER {
            let mirror = nodes.len() - 1 - i;
            assert_eq!(nodes[i].cmp(&nodes[mirror].neg()), Some(0));
            assert_eq!(w_k[i].cmp(&w_k[mirror]), Some(0));
        }
        for (i, w) in w_g.iter().enumerate() {
            if i % 2 == 0 {
                assert!(w.is_zero());
            } else {
                assert!(w.is_positive());
            }
        }
        let a = ExactNum::from_i8(-1, p);
        let b = ExactNum::from_u8(1, p);
        // ∫_{-1}^{1} x^22 dx = 2/23. The (7, 15) rule is exact through degree
        // 3*7+1 = 22, so this rational is the hand check of the whole pair.
        let got = gauss_kronrod_interval(
            |x, p, rm, _cc| {
                let mut y = ExactNum::from_u8(1, p);
                for _ in 0..22 {
                    y = y.mul(x, p, rm);
                }
                y
            },
            &a,
            &b,
            p,
            rm,
            &mut cc,
        )
        .expect("x^22");
        let two = ExactNum::from_u8(2, p);
        let twenty_three = ExactNum::from_u8(23, p);
        let want = two.div(&twenty_three, p, rm);
        let correct = bits(&got.estimate_kronrod, &want, p);
        assert!(
            correct > (p as i32) - 8,
            "x^22 bits {correct} estimate {}",
            got.estimate_kronrod
        );
    }

    #[test]
    fn wynn_geometric_series_hits_two() {
        let (p, rm) = gold_p();
        // Partial sums of Σ (1/2)^k from k = 0. One geometric component, so
        // the even ε column of order 2 is the limit 2 in exact arithmetic.
        let mut seq = Vec::new();
        let mut acc = ExactNum::new(p);
        let half = ExactNum::from_u8(1, p).ldexp(-1, p, rm);
        let mut term = ExactNum::from_u8(1, p);
        for _ in 0..4 {
            acc = acc.add(&term, p, rm);
            seq.push(acc.clone());
            term = term.mul(&half, p, rm);
        }
        let got = wynn_epsilon(&seq, p, rm).expect("wynn");
        let two = ExactNum::from_u8(2, p);
        let correct = bits(&got, &two, p);
        assert!(
            correct > (p as i32) - 8,
            "geometric bits {correct} got {got}"
        );
        assert!(wynn_epsilon(&[], p, rm).is_none());
        let nan = ExactNum::nan(None);
        assert!(wynn_epsilon(&[nan], p, rm).is_none());
    }

    #[test]
    fn gk_rejects_bad_interval_and_nonfinite_sample() {
        let (p, rm) = gold_p();
        let mut cc = Consts::new().expect("consts");
        let a = ExactNum::from_u8(0, p);
        let b = ExactNum::from_u8(1, p);
        let swapped = gauss_kronrod_interval(|x, _p, _rm, _cc| x.clone(), &b, &a, p, rm, &mut cc);
        assert_eq!(swapped.unwrap_err(), QuadratureError::InvalidArgument);
        let nan_f = gauss_kronrod_interval(
            |_x, _p, _rm, _cc| ExactNum::nan(None),
            &a,
            &b,
            p,
            rm,
            &mut cc,
        );
        assert_eq!(nan_f.unwrap_err(), QuadratureError::InvalidArgument);
        assert!(kronrod_pair(10, p, rm, &mut cc).is_none());
        let tol = ExactNum::from_u8(1, p).ldexp(-8, p, rm);
        let tight = integrate_adaptive_gk(
            |x, p, rm, _cc| x.reciprocal(p, rm).sqrt(p, rm),
            &a,
            &b,
            &tol,
            1,
            p,
            rm,
            &mut cc,
        );
        assert_eq!(tight.unwrap_err(), QuadratureError::NotConverged);
    }
}
