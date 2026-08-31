//! Chebyshev interpolation and Clenshaw evaluation on [`ExactNum`].

use crate::defs::WORD_BIT_SIZE;
use crate::Consts;
use crate::ExactNum;
use crate::RoundingMode;
use crate::NAN;
use alloc::vec::Vec;

/// Maximum number of Chebyshev coefficients [`chebyshev_coeffs`] will compute.
pub const CHEBYSHEV_MAX_DEGREE: usize = 256;

fn work_p(p: usize) -> usize {
    p.saturating_add(WORD_BIT_SIZE)
}

fn finite(x: &ExactNum) -> bool {
    !x.is_nan() && !x.is_inf()
}

/// Map `x ∈ [a, b]` to `t ∈ [-1, 1]`. `None` if `b ≤ a` or a part is non-finite.
fn map_to_unit(
    x: &ExactNum,
    a: &ExactNum,
    b: &ExactNum,
    p: usize,
    rm: RoundingMode,
) -> Option<ExactNum> {
    if !finite(x) || !finite(a) || !finite(b) {
        return None;
    }
    if a.cmp(b) != Some(-1) {
        return None;
    }
    let two = ExactNum::from_u8(2, p);
    let num = two.mul(x, p, rm).sub(a, p, rm).sub(b, p, rm);
    let den = b.sub(a, p, rm);
    if den.is_zero() {
        return None;
    }
    Some(num.div(&den, p, rm))
}

fn map_from_unit(t: &ExactNum, a: &ExactNum, b: &ExactNum, p: usize, rm: RoundingMode) -> ExactNum {
    let two = ExactNum::from_u8(2, p);
    let mid = a.add(b, p, rm).div(&two, p, rm);
    let half = b.sub(a, p, rm).div(&two, p, rm);
    mid.add(&half.mul(t, p, rm), p, rm)
}

/// Chebyshev–Gauss node `t_k = cos(π (2k+1) / (2n))` on `[-1, 1]`.
fn unit_node(k: usize, n: usize, p: usize, rm: RoundingMode, cc: &mut Consts) -> ExactNum {
    let pi = cc.pi(p, rm);
    let num = ExactNum::from_u32((2 * k + 1) as u32, p);
    let den = ExactNum::from_u32((2 * n) as u32, p);
    let theta = pi.mul(&num, p, rm).div(&den, p, rm);
    theta.cos(p, rm, cc)
}

/// Discrete cosine of order `j` at node `k`: `cos(π j (2k+1) / (2n))`.
fn node_cos_j(
    j: usize,
    k: usize,
    n: usize,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> ExactNum {
    if j == 0 {
        return ExactNum::from_u8(1, p);
    }
    let pi = cc.pi(p, rm);
    let num = ExactNum::from_u32((j * (2 * k + 1)) as u32, p);
    let den = ExactNum::from_u32((2 * n) as u32, p);
    let theta = pi.mul(&num, p, rm).div(&den, p, rm);
    theta.cos(p, rm, cc)
}

/// Interpolation coefficients of `f` on `[a, b]` at the `n` Chebyshev–Gauss nodes.
///
/// `c_0 = (1/n) Σ f(x_k)`, `c_j = (2/n) Σ f(x_k) T_j(t_k)` for `j ≥ 1`,
/// so `f(x) ≈ Σ_{j=0}^{n-1} c_j T_j(t(x))` with `t` the affine map to `[-1, 1]`.
///
/// `None` if `n` is 0 or greater than [`CHEBYSHEV_MAX_DEGREE`], or if the
/// interval is not a finite `a < b`.
pub fn chebyshev_coeffs<F>(
    mut f: F,
    n: usize,
    a: &ExactNum,
    b: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<Vec<ExactNum>>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    if n == 0 || n > CHEBYSHEV_MAX_DEGREE {
        return None;
    }
    if !finite(a) || !finite(b) || a.cmp(b) != Some(-1) {
        return None;
    }
    let wrk = work_p(p);
    let mut fx = Vec::with_capacity(n);
    for k in 0..n {
        let t = unit_node(k, n, wrk, RoundingMode::None, cc);
        let x = map_from_unit(&t, a, b, wrk, RoundingMode::None);
        let y = f(&x, wrk, RoundingMode::None, cc);
        if !finite(&y) {
            return None;
        }
        fx.push(y);
    }
    let n_f = ExactNum::from_u32(n as u32, wrk);
    let two = ExactNum::from_u8(2, wrk);
    let mut coeffs = Vec::with_capacity(n);
    for j in 0..n {
        let mut s = ExactNum::new(wrk);
        for k in 0..n {
            let w = node_cos_j(j, k, n, wrk, RoundingMode::None, cc);
            s = s.add(
                &fx[k].mul(&w, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            );
        }
        let scale = if j == 0 {
            ExactNum::from_u8(1, wrk).div(&n_f, wrk, RoundingMode::None)
        } else {
            two.div(&n_f, wrk, RoundingMode::None)
        };
        let mut c = s.mul(&scale, wrk, RoundingMode::None);
        let _ = c.set_precision(p, rm);
        coeffs.push(c);
    }
    Some(coeffs)
}

/// Clenshaw recurrence for `Σ_{k=0}^{n-1} c_k T_k(x)` at precision `p`.
///
/// `x` is the Chebyshev variable on `[-1, 1]` (no interval map). Empty
/// `coeffs` is zero. A non-finite coefficient or `x` yields `NaN`.
pub fn clenshaw(coeffs: &[ExactNum], x: &ExactNum, p: usize, rm: RoundingMode) -> ExactNum {
    if coeffs.is_empty() {
        return ExactNum::new(p);
    }
    if !finite(x) || coeffs.iter().any(|c| !finite(c)) {
        return NAN;
    }
    let wrk = work_p(p);
    let two = ExactNum::from_u8(2, wrk);
    let mut b1 = ExactNum::new(wrk);
    let mut b2 = ExactNum::new(wrk);
    for c in coeffs.iter().skip(1).rev() {
        let t = two
            .mul(x, wrk, RoundingMode::None)
            .mul(&b1, wrk, RoundingMode::None)
            .sub(&b2, wrk, RoundingMode::None)
            .add(c, wrk, RoundingMode::None);
        b2 = b1;
        b1 = t;
    }
    let mut y = x
        .mul(&b1, wrk, RoundingMode::None)
        .sub(&b2, wrk, RoundingMode::None)
        .add(&coeffs[0], wrk, RoundingMode::None);
    let _ = y.set_precision(p, rm);
    y
}

/// Evaluate the Chebyshev expansion of `coeffs` at `x ∈ [a, b]`.
///
/// Maps `x` to `[-1, 1]` and calls [`clenshaw`]. A bad interval or
/// non-finite input is `NaN`.
pub fn chebyshev_eval(
    coeffs: &[ExactNum],
    x: &ExactNum,
    a: &ExactNum,
    b: &ExactNum,
    p: usize,
    rm: RoundingMode,
) -> ExactNum {
    match map_to_unit(x, a, b, p, rm) {
        Some(t) => clenshaw(coeffs, &t, p, rm),
        None => NAN,
    }
}

/// ℓ¹ tail `Σ_{k≥1} |c_k|` at precision `p`.
///
/// On `[-1, 1]`, `|T_k| ≤ 1`, so this bounds `|Σ c_k T_k(x) - c_0|`.
/// Empty or a single coefficient is zero.
pub fn chebyshev_error_bound(coeffs: &[ExactNum], p: usize) -> ExactNum {
    let mut acc = ExactNum::new(p);
    for c in coeffs.iter().skip(1) {
        acc = acc.add(&c.abs(), p, RoundingMode::ToEven);
    }
    acc
}

/// Direct `Σ c_k T_k(x)` via the three-term recurrence (test / comparison).
#[cfg(test)]
fn chebyshev_sum_direct(coeffs: &[ExactNum], x: &ExactNum, p: usize, rm: RoundingMode) -> ExactNum {
    if coeffs.is_empty() {
        return ExactNum::new(p);
    }
    let wrk = work_p(p);
    let one = ExactNum::from_u8(1, wrk);
    let two = ExactNum::from_u8(2, wrk);
    let mut t_prev = one;
    let mut acc = coeffs[0].clone();
    let _ = acc.set_precision(wrk, RoundingMode::None);
    if coeffs.len() == 1 {
        let _ = acc.set_precision(p, rm);
        return acc;
    }
    let mut t_cur = x.clone();
    let _ = t_cur.set_precision(wrk, RoundingMode::None);
    acc = acc.add(
        &coeffs[1].mul(&t_cur, wrk, RoundingMode::None),
        wrk,
        RoundingMode::None,
    );
    for c in coeffs.iter().skip(2) {
        let t_next = two
            .mul(x, wrk, RoundingMode::None)
            .mul(&t_cur, wrk, RoundingMode::None)
            .sub(&t_prev, wrk, RoundingMode::None);
        acc = acc.add(
            &c.mul(&t_next, wrk, RoundingMode::None),
            wrk,
            RoundingMode::None,
        );
        t_prev = t_cur;
        t_cur = t_next;
    }
    let _ = acc.set_precision(p, rm);
    acc
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Consts;

    /// Plan gold: 20-term `exp` on `[-1, 1]`.
    const CHEBYSHEV_EXP_TERMS: usize = 20;
    /// Plan gold: error smaller than `10^{-15}`.
    const CHEBYSHEV_EXP_ERR_DIGITS: isize = 15;

    fn gold_p() -> (usize, RoundingMode) {
        (256, RoundingMode::ToEven)
    }

    #[test]
    fn chebyshev_exp_nodes_clenshaw() {
        let (p, rm) = gold_p();
        let mut cc = Consts::new().expect("consts");
        let a = ExactNum::from_i64(-1, p);
        let b = ExactNum::from_u8(1, p);
        let coeffs = chebyshev_coeffs(
            |x, p, rm, cc| x.exp(p, rm, cc),
            CHEBYSHEV_EXP_TERMS,
            &a,
            &b,
            p,
            rm,
            &mut cc,
        )
        .expect("coeffs");
        assert_eq!(coeffs.len(), CHEBYSHEV_EXP_TERMS);

        let ten = ExactNum::from_u8(10, p);
        let tol = ExactNum::from_u8(1, p).div(&ten.powsi(CHEBYSHEV_EXP_ERR_DIGITS, p, rm), p, rm);
        for &xi in &[-1i64, 0, 1] {
            let x = ExactNum::from_i64(xi, p);
            let approx = chebyshev_eval(&coeffs, &x, &a, &b, p, rm);
            let exact = x.exp(p, rm, &mut cc);
            let err = approx.sub(&exact, p, rm).abs();
            assert!(
                err.is_zero() || err.cmp(&tol) == Some(-1),
                "exp({}) error not < 10^{{-15}}",
                xi
            );
        }

        let node_slack =
            ExactNum::from_u8(1, p).ldexp(-((p as i32) - (crate::WORD_BIT_SIZE as i32)), p, rm);
        for k in 0..CHEBYSHEV_EXP_TERMS {
            let t = unit_node(k, CHEBYSHEV_EXP_TERMS, p, rm, &mut cc);
            let x = map_from_unit(&t, &a, &b, p, rm);
            let approx = chebyshev_eval(&coeffs, &x, &a, &b, p, rm);
            let exact = x.exp(p, rm, &mut cc);
            let err = approx.sub(&exact, p, rm).abs();
            assert!(
                err.is_zero() || err.cmp(&node_slack) == Some(-1),
                "node {k} interpolant farther than a working word from exp"
            );
        }

        const IDENTITY_NODES: usize = 4;
        let id_c = chebyshev_coeffs(
            |x, _p, _rm, _cc| x.clone(),
            IDENTITY_NODES,
            &a,
            &b,
            p,
            rm,
            &mut cc,
        )
        .expect("identity coeffs");
        let one = ExactNum::from_u8(1, p);
        assert_eq!(id_c[1].cmp(&one), Some(0));
        for k in 0..IDENTITY_NODES {
            let t = unit_node(k, IDENTITY_NODES, p, rm, &mut cc);
            let x = map_from_unit(&t, &a, &b, p, rm);
            let approx = chebyshev_eval(&id_c, &x, &a, &b, p, rm);
            let err = approx.sub(&x, p, rm).abs();
            assert!(
                err.is_zero() || err.cmp(&node_slack) == Some(-1),
                "identity node {k} farther than a working word from x"
            );
        }

        let c = [ExactNum::from_u8(1, p), ExactNum::from_u8(2, p), ExactNum::from_u8(3, p)];
        let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, rm);
        let via_clenshaw = clenshaw(&c, &half, p, rm);
        let via_direct = chebyshev_sum_direct(&c, &half, p, rm);
        assert_eq!(via_clenshaw.cmp(&half), Some(0));
        assert_eq!(via_direct.cmp(&half), Some(0));
        assert_eq!(via_clenshaw.cmp(&via_direct), Some(0));

        assert!(
            chebyshev_coeffs(|x, p, rm, cc| x.exp(p, rm, cc), 0, &a, &b, p, rm, &mut cc).is_none()
        );
        assert!(chebyshev_error_bound(&c, p).cmp(&ExactNum::from_u8(5, p)) == Some(0));
    }
}
