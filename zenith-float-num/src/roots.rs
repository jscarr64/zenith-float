//! Scalar real root finders on [`ExactNum`].

use crate::defs::WORD_BIT_SIZE;
use crate::Consts;
use crate::ExactNum;
use crate::RoundingMode;

/// Maximum iterations for every method in this module.
pub const ROOT_MAX_ITER: usize = 256;

/// Default absolute tolerance exponent: `tol = 2^{ROOT_DEFAULT_TOL}`.
pub const ROOT_DEFAULT_TOL: i32 = -256;

fn work_p(p: usize) -> usize {
    p.saturating_add(WORD_BIT_SIZE)
}

fn finite(x: &ExactNum) -> bool {
    !x.is_nan() && !x.is_inf()
}

fn opposite_signs(a: &ExactNum, b: &ExactNum) -> bool {
    if a.is_zero() || b.is_zero() || !finite(a) || !finite(b) {
        return false;
    }
    a.is_positive() != b.is_positive()
}

fn below_tol(x: &ExactNum, tol: &ExactNum) -> bool {
    x.is_zero() || x.cmp(tol) == Some(-1) || x.cmp(tol) == Some(0)
}

/// Suggested absolute tolerance `2^{ROOT_DEFAULT_TOL}` at precision `p`.
pub fn root_default_tol(p: usize, rm: RoundingMode) -> ExactNum {
    ExactNum::from_u8(1, p).ldexp(ROOT_DEFAULT_TOL, p, rm)
}

fn bisect_counted<F>(
    mut f: F,
    mut a: ExactNum,
    mut b: ExactNum,
    tol: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<(ExactNum, usize)>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    if !finite(&a) || !finite(&b) || !finite(tol) || a.cmp(&b) != Some(-1) {
        return None;
    }
    let wrk = work_p(p);
    let mut fa = f(&a, wrk, RoundingMode::None, cc);
    let mut fb = f(&b, wrk, RoundingMode::None, cc);
    if !opposite_signs(&fa, &fb) {
        return None;
    }
    let two = ExactNum::from_u8(2, wrk);
    for k in 1..=ROOT_MAX_ITER {
        let mid = a
            .add(&b, wrk, RoundingMode::None)
            .div(&two, wrk, RoundingMode::None);
        if !strictly_inside(&mid, &a, &b) {
            // The bracket is two adjacent working-precision values: `tol` is below resolution.
            let mut out = mid;
            let _ = out.set_precision(p, rm);
            return Some((out, k));
        }
        let fm = f(&mid, wrk, RoundingMode::None, cc);
        if !finite(&fm) {
            return None;
        }
        if opposite_signs(&fa, &fm) {
            b = mid.clone();
            fb = fm.clone();
        } else {
            a = mid.clone();
            fa = fm.clone();
        }
        if fm.is_zero() || below_tol(&b.sub(&a, wrk, RoundingMode::None).abs(), tol) {
            let mut out = mid;
            let _ = out.set_precision(p, rm);
            return Some((out, k));
        }
        let _ = fb;
    }
    None
}

/// Bisection on `[a, b]`. `None` if `f(a)` and `f(b)` do not have opposite signs.
pub fn bisect<F>(
    f: F,
    a: &ExactNum,
    b: &ExactNum,
    tol: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNum>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    bisect_counted(f, a.clone(), b.clone(), tol, p, rm, cc).map(|(x, _)| x)
}

fn newton_counted<F, D>(
    mut f: F,
    mut df: D,
    x0: ExactNum,
    tol: &ExactNum,
    max_iter: usize,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<(ExactNum, usize)>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
    D: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    if !finite(&x0) || !finite(tol) || max_iter == 0 || max_iter > ROOT_MAX_ITER {
        return None;
    }
    let wrk = work_p(p);
    let mut x = x0;
    for k in 1..=max_iter {
        let y = f(&x, wrk, RoundingMode::None, cc);
        let d = df(&x, wrk, RoundingMode::None, cc);
        if !finite(&y) || !finite(&d) || d.is_zero() {
            return None;
        }
        let step = y.div(&d, wrk, RoundingMode::None);
        x = x.sub(&step, wrk, RoundingMode::None);
        if !finite(&x) {
            return None;
        }
        if y.is_zero() || below_tol(&step.abs(), tol) {
            let _ = x.set_precision(p, rm);
            return Some((x, k));
        }
    }
    None
}

/// Newton–Raphson from `x0` with exact derivative `df`.
pub fn newton<F, D>(
    f: F,
    df: D,
    x0: &ExactNum,
    tol: &ExactNum,
    max_iter: usize,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNum>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
    D: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    newton_counted(f, df, x0.clone(), tol, max_iter, p, rm, cc).map(|(x, _)| x)
}

fn illinois_counted<F>(
    mut f: F,
    mut a: ExactNum,
    mut b: ExactNum,
    tol: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<(ExactNum, usize)>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    if !finite(&a) || !finite(&b) || !finite(tol) || a.cmp(&b) != Some(-1) {
        return None;
    }
    let wrk = work_p(p);
    let two = ExactNum::from_u8(2, wrk);
    let mut fa = f(&a, wrk, RoundingMode::None, cc);
    let mut fb = f(&b, wrk, RoundingMode::None, cc);
    if !opposite_signs(&fa, &fb) {
        return None;
    }
    for k in 1..=ROOT_MAX_ITER {
        let den = fb.sub(&fa, wrk, RoundingMode::None);
        if den.is_zero() {
            return None;
        }
        let c = a
            .mul(&fb, wrk, RoundingMode::None)
            .sub(
                &b.mul(&fa, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            )
            .div(&den, wrk, RoundingMode::None);
        if finite(&c) && !strictly_inside(&c, &a, &b) {
            // The secant point hit an endpoint: the bracket is at working resolution.
            let mut out = c;
            let _ = out.set_precision(p, rm);
            return Some((out, k));
        }
        let fc = f(&c, wrk, RoundingMode::None, cc);
        if !finite(&c) || !finite(&fc) {
            return None;
        }
        if fc.is_zero() || below_tol(&b.sub(&a, wrk, RoundingMode::None).abs(), tol) {
            let mut out = c;
            let _ = out.set_precision(p, rm);
            return Some((out, k));
        }
        if opposite_signs(&fa, &fc) {
            b = c;
            fb = fc;
            fa = fa.div(&two, wrk, RoundingMode::None);
        } else {
            a = c;
            fa = fc;
            fb = fb.div(&two, wrk, RoundingMode::None);
        }
    }
    None
}

/// Illinois (modified regula falsi) on `[a, b]`.
pub fn illinois<F>(
    f: F,
    a: &ExactNum,
    b: &ExactNum,
    tol: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNum>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    illinois_counted(f, a.clone(), b.clone(), tol, p, rm, cc).map(|(x, _)| x)
}

fn strictly_inside(x: &ExactNum, lo: &ExactNum, hi: &ExactNum) -> bool {
    x.cmp(lo) == Some(1) && x.cmp(hi) == Some(-1)
}

fn lt(x: &ExactNum, y: &ExactNum) -> bool {
    matches!(x.cmp(y), Some(c) if c < 0)
}

fn same_sign(x: &ExactNum, y: &ExactNum) -> bool {
    !x.is_zero() && !y.is_zero() && x.is_positive() == y.is_positive()
}

/// Brent's zero finder (R. P. Brent, *Algorithms for Minimization without Derivatives*, 1973,
/// ch. 4): keeps a sign-changing bracket `[b, c]`, tries inverse quadratic / secant steps and
/// falls back to bisection whenever they do not shrink the bracket fast enough. Stops when the
/// half-bracket is at most `2·2^{-wrk}|b| + tol/2`, so it also terminates when `tol` is below
/// the working resolution.
fn brent_counted<F>(
    mut f: F,
    mut a: ExactNum,
    mut b: ExactNum,
    tol: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<(ExactNum, usize)>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    if !finite(&a) || !finite(&b) || !finite(tol) || a.cmp(&b) != Some(-1) {
        return None;
    }
    let wrk = work_p(p);
    let none = RoundingMode::None;
    let one = ExactNum::from_u8(1, wrk);
    let two = ExactNum::from_u8(2, wrk);
    let three = ExactNum::from_u8(3, wrk);
    let half_tol = tol.abs().div(&two, wrk, none);
    let eps2 = one.ldexp(1 - wrk as i32, wrk, none);
    let mut fa = f(&a, wrk, none, cc);
    let mut fb = f(&b, wrk, none, cc);
    if !opposite_signs(&fa, &fb) {
        return None;
    }
    let mut c = a.clone();
    let mut fc = fa.clone();
    let mut d = b.sub(&a, wrk, none);
    let mut e = d.clone();
    for k in 1..=ROOT_MAX_ITER {
        if same_sign(&fb, &fc) {
            c = a.clone();
            fc = fa.clone();
            d = b.sub(&a, wrk, none);
            e = d.clone();
        }
        if lt(&fc.abs(), &fb.abs()) {
            a = b;
            b = c;
            c = a.clone();
            fa = fb;
            fb = fc;
            fc = fa.clone();
        }
        let tol1 = eps2.mul(&b.abs(), wrk, none).add(&half_tol, wrk, none);
        let m = c.sub(&b, wrk, none).div(&two, wrk, none);
        if fb.is_zero() || !lt(&tol1, &m.abs()) {
            let mut out = b;
            let _ = out.set_precision(p, rm);
            return Some((out, k));
        }
        if !lt(&e.abs(), &tol1) && lt(&fb.abs(), &fa.abs()) {
            // Secant (a == c) or inverse quadratic interpolation, as the correction p/q.
            let s = fb.div(&fa, wrk, none);
            let two_m = two.mul(&m, wrk, none);
            let (mut pn, mut q) = if a.cmp(&c) == Some(0) {
                (two_m.mul(&s, wrk, none), one.sub(&s, wrk, none))
            } else {
                let qa = fa.div(&fc, wrk, none);
                let r = fb.div(&fc, wrk, none);
                let inner = two_m
                    .mul(&qa, wrk, none)
                    .mul(&qa.sub(&r, wrk, none), wrk, none)
                    .sub(
                        &b.sub(&a, wrk, none).mul(&r.sub(&one, wrk, none), wrk, none),
                        wrk,
                        none,
                    );
                let den = qa
                    .sub(&one, wrk, none)
                    .mul(&r.sub(&one, wrk, none), wrk, none)
                    .mul(&s.sub(&one, wrk, none), wrk, none);
                (s.mul(&inner, wrk, none), den)
            };
            if pn.is_positive() {
                q = q.neg();
            } else {
                pn = pn.neg();
            }
            let lim1 = three.mul(&m, wrk, none).mul(&q, wrk, none).sub(
                &tol1.mul(&q, wrk, none).abs(),
                wrk,
                none,
            );
            let lim2 = e.mul(&q, wrk, none).abs();
            let lim = if lt(&lim1, &lim2) { lim1 } else { lim2 };
            if lt(&two.mul(&pn, wrk, none), &lim) {
                e = d;
                d = pn.div(&q, wrk, none);
            } else {
                d = m.clone();
                e = m.clone();
            }
        } else {
            d = m.clone();
            e = m.clone();
        }
        a = b.clone();
        fa = fb.clone();
        b = if lt(&tol1, &d.abs()) {
            b.add(&d, wrk, none)
        } else if m.is_positive() {
            b.add(&tol1, wrk, none)
        } else {
            b.sub(&tol1, wrk, none)
        };
        fb = f(&b, wrk, none, cc);
        if !finite(&fb) {
            return None;
        }
    }
    None
}

/// Brent's method (bisection + secant + inverse quadratic) on `[a, b]`.
pub fn brent<F>(
    f: F,
    a: &ExactNum,
    b: &ExactNum,
    tol: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNum>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    brent_counted(f, a.clone(), b.clone(), tol, p, rm, cc).map(|(x, _)| x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Consts;

    fn gold_p() -> (usize, RoundingMode) {
        (256, RoundingMode::ToEven)
    }

    #[test]
    fn roots_bisect_newton_brent_illinois() {
        let (p, rm) = gold_p();
        let mut cc = Consts::new().expect("consts");
        let tol = root_default_tol(p, rm);
        let three = ExactNum::from_u8(3, p);
        let four = ExactNum::from_u8(4, p);
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);
        let two = ExactNum::from_u8(2, p);
        let pi = cc.pi(p, rm);

        let (broot, biters) = bisect_counted(
            |x, p, rm, cc| x.sin(p, rm, cc),
            three.clone(),
            four.clone(),
            &tol,
            p,
            rm,
            &mut cc,
        )
        .expect("bisect sin");
        assert_eq!(broot.cmp(&pi), Some(0));

        let nroot = newton(
            |x, p, rm, _cc| x.mul(x, p, rm).sub(&two, p, rm),
            |x, p, rm, _cc| two.mul(x, p, rm),
            &one,
            &tol,
            ROOT_MAX_ITER,
            p,
            rm,
            &mut cc,
        )
        .expect("newton sqrt2");
        let s2 = two.sqrt(p, rm);
        assert_eq!(nroot.cmp(&s2), Some(0));

        let (rbrent, briters) = brent_counted(
            |x, p, rm, cc| x.sin(p, rm, cc),
            three.clone(),
            four.clone(),
            &tol,
            p,
            rm,
            &mut cc,
        )
        .expect("brent sin");
        assert_eq!(rbrent.cmp(&pi), Some(0));
        assert!(briters < biters);

        let iroot = illinois(
            |x, p, rm, cc| x.sin(p, rm, cc),
            &three,
            &four,
            &tol,
            p,
            rm,
            &mut cc,
        )
        .expect("illinois sin");
        assert_eq!(iroot.cmp(&pi), Some(0));

        assert!(bisect(
            |x, p, rm, cc| x.sin(p, rm, cc),
            &zero,
            &one,
            &tol,
            p,
            rm,
            &mut cc
        )
        .is_none());
    }
}
