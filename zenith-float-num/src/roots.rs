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
    let two = ExactNum::from_u8(2, wrk);
    let mut fa = f(&a, wrk, RoundingMode::None, cc);
    let mut fb = f(&b, wrk, RoundingMode::None, cc);
    if !opposite_signs(&fa, &fb) {
        return None;
    }
    let mut c = a.clone();
    let mut fc = fa.clone();
    for k in 1..=ROOT_MAX_ITER {
        if below_tol(&b.sub(&a, wrk, RoundingMode::None).abs(), tol) || fb.is_zero() {
            let mut out = b;
            let _ = out.set_precision(p, rm);
            return Some((out, k));
        }
        let den_ba = fb.sub(&fa, wrk, RoundingMode::None);
        let mut s = if !den_ba.is_zero() {
            a.mul(&fb, wrk, RoundingMode::None)
                .sub(
                    &b.mul(&fa, wrk, RoundingMode::None),
                    wrk,
                    RoundingMode::None,
                )
                .div(&den_ba, wrk, RoundingMode::None)
        } else {
            a.add(&b, wrk, RoundingMode::None)
                .div(&two, wrk, RoundingMode::None)
        };
        if a.cmp(&c) != Some(0) && b.cmp(&c) != Some(0) {
            let d1 = fa.sub(&fb, wrk, RoundingMode::None);
            let d2 = fa.sub(&fc, wrk, RoundingMode::None);
            let d3 = fb.sub(&fc, wrk, RoundingMode::None);
            if !d1.is_zero() && !d2.is_zero() && !d3.is_zero() {
                let t1 = a
                    .mul(&fb, wrk, RoundingMode::None)
                    .mul(&fc, wrk, RoundingMode::None)
                    .div(
                        &d1.mul(&d2, wrk, RoundingMode::None),
                        wrk,
                        RoundingMode::None,
                    );
                let t2 = b
                    .mul(&fa, wrk, RoundingMode::None)
                    .mul(&fc, wrk, RoundingMode::None)
                    .div(
                        &fb.sub(&fa, wrk, RoundingMode::None)
                            .mul(&d3, wrk, RoundingMode::None),
                        wrk,
                        RoundingMode::None,
                    );
                let t3 = c
                    .mul(&fa, wrk, RoundingMode::None)
                    .mul(&fb, wrk, RoundingMode::None)
                    .div(
                        &fc.sub(&fa, wrk, RoundingMode::None).mul(
                            &fc.sub(&fb, wrk, RoundingMode::None),
                            wrk,
                            RoundingMode::None,
                        ),
                        wrk,
                        RoundingMode::None,
                    );
                let iqi = t1
                    .add(&t2, wrk, RoundingMode::None)
                    .add(&t3, wrk, RoundingMode::None);
                if strictly_inside(&iqi, &a, &b) {
                    s = iqi;
                }
            }
        }
        if !strictly_inside(&s, &a, &b) {
            s = a
                .add(&b, wrk, RoundingMode::None)
                .div(&two, wrk, RoundingMode::None);
        }
        let fs = f(&s, wrk, RoundingMode::None, cc);
        if !finite(&fs) {
            return None;
        }
        c = b.clone();
        fc = fb.clone();
        if opposite_signs(&fa, &fs) {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }
        if a.cmp(&b) == Some(1) {
            core::mem::swap(&mut a, &mut b);
            core::mem::swap(&mut fa, &mut fb);
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
