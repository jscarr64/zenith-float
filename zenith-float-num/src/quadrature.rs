//! Gauss and tanh–sinh quadrature on [`ExactNum`].

use crate::defs::WORD_BIT_SIZE;
use crate::Consts;
use crate::ExactNum;
use crate::RoundingMode;
use alloc::vec::Vec;

/// Maximum number of Gauss nodes. Larger `n` is `None`.
pub const QUADRATURE_MAX_NODES: usize = 64;

/// Maximum tanh–sinh step halvings after the precision-based `h`.
pub const TANH_SINH_LEVELS_MAX: usize = 8;

/// Newton sweeps per Gauss root.
const QUADRATURE_NEWTON_MAX: usize = 64;

/// Grid points per expected root when isolating Laguerre / Hermite zeros.
const QUADRATURE_BRACKET_MUL: usize = 8;

/// Bisection steps inside a sign-change bracket before Newton.
const QUADRATURE_BISECT_STEPS: usize = 16;

/// Hard cap on tanh–sinh sample index `|k|`.
const TANH_SINH_K_MAX: usize = 512;

fn work_p(p: usize) -> usize {
    p.saturating_add(WORD_BIT_SIZE)
}

/// Extra word so `1-x²` stays nonzero at tanh–sinh nodes until the Jacobian
/// is smaller than a unit in the last place of a `p`-bit result.
fn work_p_de(p: usize) -> usize {
    p.saturating_mul(2).saturating_add(WORD_BIT_SIZE)
}

fn finite(x: &ExactNum) -> bool {
    !x.is_nan() && !x.is_inf()
}

fn tiny(p: usize, extra: i32) -> ExactNum {
    let e = (p as i32).saturating_sub(extra).saturating_neg();
    ExactNum::from_u8(1, p).ldexp(e, p, RoundingMode::ToEven)
}

fn below(x: &ExactNum, bound: &ExactNum) -> bool {
    x.is_zero() || x.cmp(bound) == Some(-1)
}

fn factorial(n: usize, p: usize, rm: RoundingMode) -> ExactNum {
    let mut acc = ExactNum::from_u8(1, p);
    for k in 2..=n {
        acc = acc.mul(&ExactNum::from_u32(k as u32, p), p, rm);
    }
    acc
}

/// Newton for `f(x)=0` with analytic `df`, starting at `x0`.
fn newton<F, D>(
    mut f: F,
    mut df: D,
    x0: ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNum>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
    D: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    let guard = tiny(p, 8);
    let mut x = x0;
    for _ in 0..QUADRATURE_NEWTON_MAX {
        let y = f(&x, p, rm, cc);
        let d = df(&x, p, rm, cc);
        if !finite(&y) || !finite(&d) || d.is_zero() {
            return None;
        }
        let step = y.div(&d, p, rm);
        x = x.sub(&step, p, rm);
        if !finite(&x) {
            return None;
        }
        if below(&step.abs(), &guard) {
            return Some(x);
        }
    }
    None
}

fn bisect_zero<F>(
    mut f: F,
    mut a: ExactNum,
    mut b: ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNum>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    let fa0 = f(&a, p, rm, cc);
    let fb0 = f(&b, p, rm, cc);
    if !finite(&fa0) || !finite(&fb0) {
        return None;
    }
    if fa0.is_zero() {
        return Some(a);
    }
    if fb0.is_zero() {
        return Some(b);
    }
    if fa0.is_positive() == fb0.is_positive() {
        return None;
    }
    let two = ExactNum::from_u8(2, p);
    for _ in 0..QUADRATURE_BISECT_STEPS {
        let m = a.add(&b, p, rm).div(&two, p, rm);
        let fm = f(&m, p, rm, cc);
        if !finite(&fm) {
            return None;
        }
        if fm.is_zero() {
            return Some(m);
        }
        if fm.is_positive() == fa0.is_positive() {
            a = m;
        } else {
            b = m;
        }
    }
    Some(a.add(&b, p, rm).div(&two, p, rm))
}

fn isolate_positive<F>(
    mut f: F,
    lo: &ExactNum,
    hi: &ExactNum,
    want: usize,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<Vec<(ExactNum, ExactNum)>>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    if want == 0 {
        return Some(Vec::new());
    }
    let n_grid = want
        .saturating_mul(QUADRATURE_BRACKET_MUL)
        .saturating_add(4);
    let step = hi
        .sub(lo, p, rm)
        .div(&ExactNum::from_u32(n_grid as u32, p), p, rm);
    let mut prev_x = lo.clone();
    let mut prev_f = f(lo, p, rm, cc);
    if !finite(&prev_f) {
        return None;
    }
    let mut out = Vec::new();
    for i in 1..=n_grid {
        let x = if i == n_grid {
            hi.clone()
        } else {
            lo.add(&step.mul(&ExactNum::from_u32(i as u32, p), p, rm), p, rm)
        };
        let y = f(&x, p, rm, cc);
        if !finite(&y) {
            return None;
        }
        if prev_f.is_zero() {
            out.push((prev_x.clone(), x.clone()));
        } else if y.is_zero() || prev_f.is_positive() != y.is_positive() {
            out.push((prev_x, x.clone()));
        }
        if out.len() == want {
            return Some(out);
        }
        prev_x = x;
        prev_f = y;
    }
    if out.len() == want {
        Some(out)
    } else {
        None
    }
}

fn map_ab(xi: &ExactNum, a: &ExactNum, b: &ExactNum, p: usize, rm: RoundingMode) -> ExactNum {
    let two = ExactNum::from_u8(2, p);
    let mid = a.add(b, p, rm).div(&two, p, rm);
    let half = b.sub(a, p, rm).div(&two, p, rm);
    mid.add(&half.mul(xi, p, rm), p, rm)
}

fn gauss_legendre_nodes(
    n: usize,
    p: usize,
    _rm: RoundingMode,
    cc: &mut Consts,
) -> Option<Vec<(ExactNum, ExactNum)>> {
    if n == 0 || n > QUADRATURE_MAX_NODES {
        return None;
    }
    let wrk = work_p(p);
    let nu = n as u32;
    let pi = cc.pi(wrk, RoundingMode::None);
    let den = ExactNum::from_u32((4 * n + 2) as u32, wrk);
    let n_f = ExactNum::from_u32(nu, wrk);
    let one = ExactNum::from_u8(1, wrk);
    let two = ExactNum::from_u8(2, wrk);
    let mut nodes = Vec::with_capacity(n);
    for k in 1..=n {
        let num = ExactNum::from_u32((4 * k - 1) as u32, wrk);
        let theta = pi
            .mul(&num, wrk, RoundingMode::None)
            .div(&den, wrk, RoundingMode::None);
        let x0 = theta.cos(wrk, RoundingMode::None, cc);
        let x = newton(
            |x, p, rm, _cc| x.legendre_p(nu, p, rm),
            |x, p, rm, _cc| {
                let pn = x.legendre_p(nu, p, rm);
                let pnm = if nu == 0 { ExactNum::new(p) } else { x.legendre_p(nu - 1, p, rm) };
                let d = one.sub(&x.mul(x, p, rm), p, rm);
                n_f.mul(&pnm.sub(&x.mul(&pn, p, rm), p, rm), p, rm)
                    .div(&d, p, rm)
            },
            x0,
            wrk,
            RoundingMode::None,
            cc,
        )?;
        let pn = x.legendre_p(nu, wrk, RoundingMode::None);
        let pnm = if nu == 0 {
            ExactNum::new(wrk)
        } else {
            x.legendre_p(nu - 1, wrk, RoundingMode::None)
        };
        let d = one.sub(&x.mul(&x, wrk, RoundingMode::None), wrk, RoundingMode::None);
        let dp = n_f
            .mul(
                &pnm.sub(
                    &x.mul(&pn, wrk, RoundingMode::None),
                    wrk,
                    RoundingMode::None,
                ),
                wrk,
                RoundingMode::None,
            )
            .div(&d, wrk, RoundingMode::None);
        let w = two.div(
            &d.mul(
                &dp.mul(&dp, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            ),
            wrk,
            RoundingMode::None,
        );
        if !finite(&x) || !finite(&w) {
            return None;
        }
        nodes.push((x, w));
    }
    Some(nodes)
}

/// Gauss–Legendre quadrature of `f` on `[a, b]` with `n` nodes.
///
/// Nodes are roots of `P_n`; weights `w_i = 2 / ((1-x_i²) [P_n'(x_i)]²)`.
/// Exact for polynomials of degree `≤ 2n−1`. `None` if `n` is 0 or greater
/// than [`QUADRATURE_MAX_NODES`], the interval is not a finite `a < b`, or
/// a root fails to converge.
pub fn gauss_legendre<F>(
    mut f: F,
    a: &ExactNum,
    b: &ExactNum,
    n: usize,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNum>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    if !finite(a) || !finite(b) || a.cmp(b) != Some(-1) {
        return None;
    }
    let wrk = work_p(p);
    let nodes = gauss_legendre_nodes(n, wrk, RoundingMode::None, cc)?;
    let two = ExactNum::from_u8(2, wrk);
    let half = b
        .sub(a, wrk, RoundingMode::None)
        .div(&two, wrk, RoundingMode::None);
    let mut acc = ExactNum::new(wrk);
    for (xi, wi) in &nodes {
        let x = map_ab(xi, a, b, wrk, RoundingMode::None);
        let y = f(&x, wrk, RoundingMode::None, cc);
        if !finite(&y) {
            return None;
        }
        acc = acc.add(
            &wi.mul(&y, wrk, RoundingMode::None),
            wrk,
            RoundingMode::None,
        );
    }
    let mut out = half.mul(&acc, wrk, RoundingMode::None);
    let _ = out.set_precision(p, rm);
    Some(out)
}

fn tanh_sinh_sum<F>(
    mut f: F,
    a: &ExactNum,
    b: &ExactNum,
    h: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNum>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    let wrk = work_p_de(p);
    let two = ExactNum::from_u8(2, wrk);
    let half = b
        .sub(a, wrk, RoundingMode::None)
        .div(&two, wrk, RoundingMode::None);
    let mid = a
        .add(b, wrk, RoundingMode::None)
        .div(&two, wrk, RoundingMode::None);
    let pi = cc.pi(wrk, RoundingMode::None);
    let half_pi = pi.div(&two, wrk, RoundingMode::None);
    let wmin = tiny(wrk, 8);
    let mut acc = ExactNum::new(wrk);
    for k in 0..TANH_SINH_K_MAX {
        let t = h.mul(&ExactNum::from_u32(k as u32, wrk), wrk, RoundingMode::None);
        let (sh, ch) = t.sinh_cosh(wrk, RoundingMode::None, cc);
        let u = half_pi.mul(&sh, wrk, RoundingMode::None);
        let (su, cu) = u.sinh_cosh(wrk, RoundingMode::None, cc);
        if !finite(&cu) || cu.is_zero() {
            break;
        }
        let xi = su.div(&cu, wrk, RoundingMode::None);
        let w = half_pi.mul(&ch, wrk, RoundingMode::None).div(
            &cu.mul(&cu, wrk, RoundingMode::None),
            wrk,
            RoundingMode::None,
        );
        if !finite(&w) {
            break;
        }
        if k > 0 && below(&w.abs(), &wmin) {
            break;
        }
        let mut eval = |s: &ExactNum| {
            let x = mid.add(
                &half.mul(s, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            );
            f(&x, wrk, RoundingMode::None, cc)
        };
        let y0 = eval(&xi);
        if !finite(&y0) {
            if k == 0 {
                return None;
            }
            break;
        }
        acc = acc.add(
            &w.mul(&y0, wrk, RoundingMode::None),
            wrk,
            RoundingMode::None,
        );
        if k > 0 {
            let y1 = eval(&xi.neg());
            if finite(&y1) {
                acc = acc.add(
                    &w.mul(&y1, wrk, RoundingMode::None),
                    wrk,
                    RoundingMode::None,
                );
            }
        }
    }
    let mut out = half
        .mul(h, wrk, RoundingMode::None)
        .mul(&acc, wrk, RoundingMode::None);
    let _ = out.set_precision(p, rm);
    Some(out)
}

/// Tanh–sinh (double-exponential) quadrature of `f` on `[a, b]`.
///
/// `x = tanh((π/2) sinh t)`, step `h = 2π / (p ln 2)` then up to
/// [`TANH_SINH_LEVELS_MAX`] successive halvings until successive sums
/// agree to working precision. `None` if the interval is not a finite `a < b`.
pub fn tanh_sinh<F>(
    mut f: F,
    a: &ExactNum,
    b: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNum>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    if !finite(a) || !finite(b) || a.cmp(b) != Some(-1) {
        return None;
    }
    let wrk = work_p(p);
    let two = ExactNum::from_u8(2, wrk);
    let pi = cc.pi(wrk, RoundingMode::None);
    let ln2 = cc.ln_2(wrk, RoundingMode::None);
    let p_f = ExactNum::from_u32(p as u32, wrk);
    let h0 = two.mul(&pi, wrk, RoundingMode::None).div(
        &p_f.mul(&ln2, wrk, RoundingMode::None),
        wrk,
        RoundingMode::None,
    );
    let guard = tiny(p, 8);
    let mut prev: Option<ExactNum> = None;
    let mut best = None;
    let two_wrk = ExactNum::from_u8(2, wrk);
    let mut h = h0;
    for level in 0..TANH_SINH_LEVELS_MAX {
        if level > 0 {
            h = h.div(&two_wrk, wrk, RoundingMode::None);
        }
        let s = tanh_sinh_sum(&mut f, a, b, &h, p, rm, cc)?;
        if let Some(pr) = &prev {
            let err = s.sub(pr, p, rm).abs();
            if below(&err, &guard) {
                return Some(s);
            }
        }
        prev = Some(s.clone());
        best = Some(s);
    }
    best
}

fn gauss_laguerre_nodes(
    n: usize,
    p: usize,
    _rm: RoundingMode,
    cc: &mut Consts,
) -> Option<Vec<(ExactNum, ExactNum)>> {
    if n == 0 || n > QUADRATURE_MAX_NODES {
        return None;
    }
    let wrk = work_p(p);
    let nu = n as u32;
    let zero = ExactNum::new(wrk);
    let hi = ExactNum::from_u32((4 * n + 4) as u32, wrk);
    let brackets = isolate_positive(
        |x, p, rm, _cc| x.laguerre(n, p, rm),
        &zero,
        &hi,
        n,
        wrk,
        RoundingMode::None,
        cc,
    )?;
    let n_f = ExactNum::from_u32(nu, wrk);
    let np1 = ExactNum::from_u32((n + 1) as u32, wrk);
    let np1sq = np1.mul(&np1, wrk, RoundingMode::None);
    let mut nodes = Vec::with_capacity(n);
    for (a, b) in brackets {
        let x0 = bisect_zero(
            |x, p, rm, _cc| x.laguerre(n, p, rm),
            a,
            b,
            wrk,
            RoundingMode::None,
            cc,
        )?;
        let x = newton(
            |x, p, rm, _cc| x.laguerre(n, p, rm),
            |x, p, rm, _cc| {
                if x.is_zero() {
                    return ExactNum::nan(None);
                }
                let ln = x.laguerre(n, p, rm);
                let lnm = if n == 0 { ExactNum::new(p) } else { x.laguerre(n - 1, p, rm) };
                n_f.mul(&ln.sub(&lnm, p, rm), p, rm).div(x, p, rm)
            },
            x0,
            wrk,
            RoundingMode::None,
            cc,
        )?;
        let ln1 = x.laguerre(n + 1, wrk, RoundingMode::None);
        let w = x.div(
            &np1sq.mul(
                &ln1.mul(&ln1, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            ),
            wrk,
            RoundingMode::None,
        );
        if !finite(&x) || !finite(&w) {
            return None;
        }
        nodes.push((x, w));
    }
    Some(nodes)
}

/// Gauss–Laguerre quadrature: `∫₀^∞ f(x) e^{-x} dx` with `n` nodes.
///
/// Nodes are roots of `L_n`; weights `w_i = x_i / ((n+1)² [L_{n+1}(x_i)]²)`.
/// Exact for `f` a polynomial of degree `≤ 2n−1`.
pub fn gauss_laguerre<F>(
    mut f: F,
    n: usize,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNum>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    let wrk = work_p(p);
    let nodes = gauss_laguerre_nodes(n, wrk, RoundingMode::None, cc)?;
    let mut acc = ExactNum::new(wrk);
    for (xi, wi) in &nodes {
        let y = f(xi, wrk, RoundingMode::None, cc);
        if !finite(&y) {
            return None;
        }
        acc = acc.add(
            &wi.mul(&y, wrk, RoundingMode::None),
            wrk,
            RoundingMode::None,
        );
    }
    let _ = acc.set_precision(p, rm);
    Some(acc)
}

fn gauss_hermite_nodes(
    n: usize,
    p: usize,
    _rm: RoundingMode,
    cc: &mut Consts,
) -> Option<Vec<(ExactNum, ExactNum)>> {
    if n == 0 || n > QUADRATURE_MAX_NODES {
        return None;
    }
    let wrk = work_p(p);
    let nu = n as u32;
    let two = ExactNum::from_u8(2, wrk);
    let hi = ExactNum::from_u32((4 * n + 2) as u32, wrk).sqrt(wrk, RoundingMode::None);
    let zero = ExactNum::new(wrk);
    let want_pos = n / 2;
    let brackets = isolate_positive(
        |x, p, rm, _cc| x.hermite_h(n, p, rm),
        &zero,
        &hi,
        want_pos,
        wrk,
        RoundingMode::None,
        cc,
    )?;
    let n_f = ExactNum::from_u32(nu, wrk);
    let nsq = n_f.mul(&n_f, wrk, RoundingMode::None);
    let fact = factorial(n, wrk, RoundingMode::None);
    let pi = cc.pi(wrk, RoundingMode::None);
    let sqrt_pi = pi.sqrt(wrk, RoundingMode::None);
    let two_pow = ExactNum::from_u8(1, wrk).ldexp((n as i32) - 1, wrk, RoundingMode::None);
    let scale = two_pow
        .mul(&fact, wrk, RoundingMode::None)
        .mul(&sqrt_pi, wrk, RoundingMode::None);
    let weight = |x: &ExactNum| {
        let hm = if n == 0 {
            ExactNum::from_u8(1, wrk)
        } else {
            x.hermite_h(n - 1, wrk, RoundingMode::None)
        };
        scale.div(
            &nsq.mul(
                &hm.mul(&hm, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            ),
            wrk,
            RoundingMode::None,
        )
    };
    let mut nodes = Vec::with_capacity(n);
    if n % 2 == 1 {
        let w0 = weight(&zero);
        if !finite(&w0) {
            return None;
        }
        nodes.push((zero.clone(), w0));
    }
    for (a, b) in brackets {
        let x0 = bisect_zero(
            |x, p, rm, _cc| x.hermite_h(n, p, rm),
            a,
            b,
            wrk,
            RoundingMode::None,
            cc,
        )?;
        let x = newton(
            |x, p, rm, _cc| x.hermite_h(n, p, rm),
            |x, p, rm, _cc| {
                if n == 0 {
                    return ExactNum::new(p);
                }
                two.mul(&n_f, p, rm).mul(&x.hermite_h(n - 1, p, rm), p, rm)
            },
            x0,
            wrk,
            RoundingMode::None,
            cc,
        )?;
        let w = weight(&x);
        if !finite(&x) || !finite(&w) {
            return None;
        }
        nodes.push((x.clone(), w.clone()));
        nodes.push((x.neg(), w));
    }
    if nodes.len() != n {
        return None;
    }
    Some(nodes)
}

/// Gauss–Hermite quadrature: `∫_{-∞}^{∞} f(x) e^{-x²} dx` with `n` nodes.
///
/// Nodes are roots of the physicist's `H_n`; weights
/// `w_i = 2^{n-1} n! √π / (n² [H_{n-1}(x_i)]²)`.
/// Exact for `f` a polynomial of degree `≤ 2n−1`.
pub fn gauss_hermite<F>(
    mut f: F,
    n: usize,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNum>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    let wrk = work_p(p);
    let nodes = gauss_hermite_nodes(n, wrk, RoundingMode::None, cc)?;
    let mut acc = ExactNum::new(wrk);
    for (xi, wi) in &nodes {
        let y = f(xi, wrk, RoundingMode::None, cc);
        if !finite(&y) {
            return None;
        }
        acc = acc.add(
            &wi.mul(&y, wrk, RoundingMode::None),
            wrk,
            RoundingMode::None,
        );
    }
    let _ = acc.set_precision(p, rm);
    Some(acc)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Consts;

    fn gold_p() -> (usize, RoundingMode) {
        (256, RoundingMode::ToEven)
    }

    #[test]
    fn quadrature_gl_tanh_laguerre_hermite() {
        let (p, rm) = gold_p();
        let mut cc = Consts::new().expect("consts");
        let a = ExactNum::from_i64(-1, p);
        let b = ExactNum::from_u8(1, p);
        let two = ExactNum::from_u8(2, p);
        let three = ExactNum::from_u8(3, p);

        let gl = gauss_legendre(|x, p, rm, _cc| x.mul(x, p, rm), &a, &b, 10, p, rm, &mut cc)
            .expect("gl x^2");
        let two_thirds = two.div(&three, p, rm);
        assert_eq!(gl.cmp(&two_thirds), Some(0));

        const GL_EXACT_DEG: u32 = 38;
        const GL_EXACT_NODES: usize = 20;
        let gl38 = gauss_legendre(
            |x, p, rm, _cc| {
                let mut y = ExactNum::from_u8(1, p);
                for _ in 0..GL_EXACT_DEG {
                    y = y.mul(x, p, rm);
                }
                y
            },
            &a,
            &b,
            GL_EXACT_NODES,
            p,
            rm,
            &mut cc,
        )
        .expect("gl x^38");
        let thirty_nine = ExactNum::from_u32(GL_EXACT_DEG + 1, p);
        let two_over = two.div(&thirty_nine, p, rm);
        assert_eq!(gl38.cmp(&two_over), Some(0));

        let ts = tanh_sinh(
            |x, p, rm, _cc| {
                let one = ExactNum::from_u8(1, p);
                one.sub(&x.mul(x, p, rm), p, rm)
                    .sqrt(p, rm)
                    .reciprocal(p, rm)
            },
            &a,
            &b,
            p,
            rm,
            &mut cc,
        )
        .expect("tanh-sinh arcsine");
        let pi = cc.pi(p, rm);
        assert_eq!(ts.cmp(&pi), Some(0));

        let lag = gauss_laguerre(|x, p, rm, _cc| x.mul(x, p, rm), 10, p, rm, &mut cc)
            .expect("laguerre x^2");
        assert_eq!(lag.cmp(&two), Some(0));

        let gh = gauss_hermite(|_x, p, _rm, _cc| ExactNum::from_u8(1, p), 1, p, rm, &mut cc)
            .expect("hermite 1");
        let sqrt_pi = pi.sqrt(p, rm);
        assert_eq!(gh.cmp(&sqrt_pi), Some(0));

        assert!(gauss_legendre(|x, _p, _rm, _cc| x.clone(), &b, &a, 4, p, rm, &mut cc).is_none());
        assert!(gauss_legendre(|x, _p, _rm, _cc| x.clone(), &a, &b, 0, p, rm, &mut cc).is_none());
    }
}
