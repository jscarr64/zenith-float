//! Scalar real root finders on [`ExactNum`].

use crate::defs::Exponent;
use crate::defs::WORD_BIT_SIZE;
use crate::Consts;
use crate::ExactNum;
use crate::RoundingMode;
use core::ops::Neg;

/// Former fixed iteration limit, kept for API compatibility. The bracketing methods
/// stop on their tolerance or at working resolution, and Newton honours only the caller's
/// `max_iter`; no method uses this as a hidden cap.
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
    x.is_zero() || x.compare(tol) == Some(-1) || x.compare(tol) == Some(0)
}

/// Suggested absolute tolerance `2^{ROOT_DEFAULT_TOL}` at precision `p`.
pub fn root_default_tol(p: usize, rm: RoundingMode) -> ExactNum {
    ExactNum::from_u8(1, p).ldexp(Exponent::from(ROOT_DEFAULT_TOL), p, rm)
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
    if !finite(&a) || !finite(&b) || !finite(tol) || a.compare(&b) != Some(-1) {
        return None;
    }
    let wrk = work_p(p);
    let mut fa = f(&a, wrk, RoundingMode::None, cc);
    let mut fb = f(&b, wrk, RoundingMode::None, cc);
    if !opposite_signs(&fa, &fb) {
        return None;
    }
    let two = ExactNum::from_u8(2, wrk);
    let mut k = 0usize;
    loop {
        k += 1;
        let mid =
            ExactNum::add_ref(&a, &b, wrk, RoundingMode::None).div(&two, wrk, RoundingMode::None);
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
        // The new bracket is [a, mid] or [mid, b]; its width is checked before mid moves into it.
        let left = opposite_signs(&fa, &fm);
        if fm.is_zero()
            || below_tol(
                &if left {
                    ExactNum::sub_ref(&mid, &a, wrk, RoundingMode::None)
                } else {
                    ExactNum::sub_ref(&b, &mid, wrk, RoundingMode::None)
                }
                .abs(),
                tol,
            )
        {
            let mut out = mid;
            let _ = out.set_precision(p, rm);
            return Some((out, k));
        }
        if left {
            b = mid;
            fb = fm;
        } else {
            a = mid;
            fa = fm;
        }
        let _ = fb;
    }
}

/// Bisection on `[a, b]`. `None` if `f(a)` and `f(b)` do not have opposite signs. The bracket ends
/// are taken by value: they become the working bracket (no copy).
pub fn bisect<F>(
    f: F,
    a: ExactNum,
    b: ExactNum,
    tol: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNum>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    bisect_counted(f, a, b, tol, p, rm, cc).map(|(x, _)| x)
}

/// Stop tolerance and iteration cap of [`newton`] (zenith-float 2.0: grouped so no public function
/// takes more than 7 arguments).
#[derive(Clone, Copy, Debug)]
pub struct NewtonOptions<'a> {
    /// Absolute step tolerance (finite); e.g. [`root_default_tol`].
    pub tol: &'a ExactNum,
    /// Most Newton steps (at least 1).
    pub max_iter: usize,
}

fn newton_counted<F, D>(
    mut f: F,
    mut df: D,
    x0: ExactNum,
    NewtonOptions { tol, max_iter }: NewtonOptions<'_>,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<(ExactNum, usize)>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
    D: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    if !finite(&x0) || !finite(tol) || max_iter == 0 {
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
        let y_zero = y.is_zero();
        let step = y.div(&d, wrk, RoundingMode::None);
        x = x.sub(&step, wrk, RoundingMode::None);
        if !finite(&x) {
            return None;
        }
        // A step below 2^-(p+8) relative to x is beneath the output grid: further
        // iterations cannot change the rounded result, even if `tol` is finer.
        let at_resolution = !x.is_zero()
            && step.exponent().zip(x.exponent()).is_some_and(|(es, ex)| {
                i128::from(es) + i128::try_from(p).unwrap_or(i128::MAX) + 8 <= i128::from(ex)
            });
        if y_zero || at_resolution || below_tol(&step.abs(), tol) {
            let _ = x.set_precision(p, rm);
            return Some((x, k));
        }
    }
    None
}

/// Newton–Raphson from `x0` with exact derivative `df`. `x0` is taken by value: it becomes the
/// first iterate (no copy).
pub fn newton<F, D>(
    f: F,
    df: D,
    x0: ExactNum,
    NewtonOptions { tol, max_iter }: NewtonOptions<'_>,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNum>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
    D: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    newton_counted(f, df, x0, NewtonOptions { tol, max_iter }, p, rm, cc).map(|(x, _)| x)
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
    if !finite(&a) || !finite(&b) || !finite(tol) || a.compare(&b) != Some(-1) {
        return None;
    }
    let wrk = work_p(p);
    let two = ExactNum::from_u8(2, wrk);
    let mut fa = f(&a, wrk, RoundingMode::None, cc);
    let mut fb = f(&b, wrk, RoundingMode::None, cc);
    if !opposite_signs(&fa, &fb) {
        return None;
    }
    let mut k = 0usize;
    loop {
        k += 1;
        let den = ExactNum::sub_ref(&fb, &fa, wrk, RoundingMode::None);
        if den.is_zero() {
            return None;
        }
        let c = ExactNum::mul_ref(&a, &fb, wrk, RoundingMode::None)
            .sub(
                &ExactNum::mul_ref(&b, &fa, wrk, RoundingMode::None),
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
        if fc.is_zero()
            || below_tol(
                &ExactNum::sub_ref(&b, &a, wrk, RoundingMode::None).abs(),
                tol,
            )
        {
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
}

/// Illinois (modified regula falsi) on `[a, b]`. The bracket ends are taken by value: they become
/// the working bracket (no copy).
pub fn illinois<F>(
    f: F,
    a: ExactNum,
    b: ExactNum,
    tol: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNum>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    illinois_counted(f, a, b, tol, p, rm, cc).map(|(x, _)| x)
}

fn strictly_inside(x: &ExactNum, lo: &ExactNum, hi: &ExactNum) -> bool {
    x.compare(lo) == Some(1) && x.compare(hi) == Some(-1)
}

fn lt(x: &ExactNum, y: &ExactNum) -> bool {
    matches!(x.compare(y), Some(c) if c < 0)
}

/// `lt(|x|, |y|)` without forming either magnitude.
fn abs_lt(x: &ExactNum, y: &ExactNum) -> bool {
    matches!(x.abs_cmp(y), Some(c) if c < 0)
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
    if !finite(&a) || !finite(&b) || !finite(tol) || a.compare(&b) != Some(-1) {
        return None;
    }
    let wrk = work_p(p);
    let none = RoundingMode::None;
    let one = ExactNum::from_u8(1, wrk);
    let two = ExactNum::from_u8(2, wrk);
    let half_tol = if tol.is_negative() {
        tol.neg_via_mul().div(&two, wrk, none)
    } else {
        ExactNum::div_ref(tol, &two, wrk, none)
    };
    let eps2 = ExactNum::from_u8(1, wrk).ldexp(1 - wrk as Exponent, wrk, none);
    let mut fa = f(&a, wrk, none, cc);
    let mut fb = f(&b, wrk, none, cc);
    if !opposite_signs(&fa, &fb) {
        return None;
    }
    // Brent's c equals a after every `c = a`, and e equals d after every `e = d`: `None` stands for
    // that equal value, so it is read in place instead of copied. c and fc are always set together.
    let mut c: Option<ExactNum> = None;
    let mut fc: Option<ExactNum> = None;
    let mut d = ExactNum::sub_ref(&b, &a, wrk, none);
    let mut e: Option<ExactNum> = None;
    let mut k = 0usize;
    loop {
        k += 1;
        if same_sign(&fb, fc.as_ref().unwrap_or(&fa)) {
            c = None;
            fc = None;
            d = ExactNum::sub_ref(&b, &a, wrk, none);
            e = None;
        }
        if abs_lt(fc.as_ref().unwrap_or(&fa), &fb) {
            // a = b; b = c; c = a
            let old_a = core::mem::replace(&mut a, b);
            b = c.take().unwrap_or(old_a);
            let old_fa = core::mem::replace(&mut fa, fb);
            fb = fc.take().unwrap_or(old_fa);
        }
        let tol1 = if b.is_negative() {
            b.neg_via_mul().mul(&eps2, wrk, none)
        } else {
            ExactNum::mul_ref(&eps2, &b, wrk, none)
        }
        .add(&half_tol, wrk, none);
        let m = ExactNum::sub_ref(c.as_ref().unwrap_or(&a), &b, wrk, none).div(&two, wrk, none);
        let m_positive = m.is_positive();
        if fb.is_zero() || !abs_lt(&tol1, &m) {
            let mut out = b;
            let _ = out.set_precision(p, rm);
            return Some((out, k));
        }
        if !abs_lt(e.as_ref().unwrap_or(&d), &tol1) && abs_lt(&fb, &fa) {
            // Secant (a == c) or inverse quadratic interpolation, as the correction p/q.
            let s = ExactNum::div_ref(&fb, &fa, wrk, none);
            let two_m = ExactNum::from_u8(2, wrk).mul(&m, wrk, none);
            let (mut pn, mut q) = if c.as_ref().is_none_or(|cv| a.compare(cv) == Some(0)) {
                (two_m.mul(&s, wrk, none), s.neg().add(&one, wrk, none))
            } else {
                // c is Some here (c == a takes the secant branch); a is set again below, so its
                // value is moved into -a + b, and the stand-in is never read
                let fcr = fc.as_ref().unwrap_or(&fa);
                let qa = ExactNum::div_ref(&fa, fcr, wrk, none);
                let r = ExactNum::div_ref(&fb, fcr, wrk, none);
                let inner = two_m
                    .mul(&qa, wrk, none)
                    .mul(&ExactNum::sub_ref(&qa, &r, wrk, none), wrk, none)
                    .sub(
                        &core::mem::replace(&mut a, ExactNum::nan(None))
                            .neg()
                            .add(&b, wrk, none)
                            .mul(
                                &ExactNum::from_u8(1, wrk).neg().add(&r, wrk, none),
                                wrk,
                                none,
                            ),
                        wrk,
                        none,
                    );
                let den = qa
                    .sub(&one, wrk, none)
                    .mul(&r.sub(&one, wrk, none), wrk, none)
                    .mul(
                        &ExactNum::from_u8(1, wrk).neg().add(&s, wrk, none),
                        wrk,
                        none,
                    );
                (s.mul(&inner, wrk, none), den)
            };
            if pn.is_positive() {
                q = q.neg();
            } else {
                pn = pn.neg();
            }
            let lim1 = ExactNum::from_u8(3, wrk)
                .mul(&m, wrk, none)
                .mul(&q, wrk, none)
                .sub(&ExactNum::mul_ref(&tol1, &q, wrk, none).abs(), wrk, none);
            let lim2 = match e.take() {
                Some(ev) => ev.mul(&q, wrk, none),
                None => ExactNum::mul_ref(&d, &q, wrk, none),
            }
            .abs();
            let lim = if lt(&lim1, &lim2) { lim1 } else { lim2 };
            if lt(&ExactNum::from_u8(2, wrk).mul(&pn, wrk, none), &lim) {
                e = Some(core::mem::replace(&mut d, pn.div(&q, wrk, none)));
            } else {
                d = m;
                e = None;
            }
        } else {
            d = m;
            e = None;
        }
        let next_b = if abs_lt(&tol1, &d) {
            ExactNum::add_ref(&b, &d, wrk, none)
        } else if m_positive {
            ExactNum::add_ref(&b, &tol1, wrk, none)
        } else {
            ExactNum::sub_ref(&b, &tol1, wrk, none)
        };
        // a = b: the old a keeps living as c where c stood for it
        let old_a = core::mem::replace(&mut a, core::mem::replace(&mut b, next_b));
        if c.is_none() {
            c = Some(old_a);
        }
        let old_fa = core::mem::replace(&mut fa, fb);
        if fc.is_none() {
            fc = Some(old_fa);
        }
        fb = f(&b, wrk, none, cc);
        if !finite(&fb) {
            return None;
        }
    }
}

/// Brent's method (bisection + secant + inverse quadratic) on `[a, b]`. The bracket ends are taken
/// by value: they become the working bracket (no copy).
pub fn brent<F>(
    f: F,
    a: ExactNum,
    b: ExactNum,
    tol: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNum>
where
    F: FnMut(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
{
    brent_counted(f, a, b, tol, p, rm, cc).map(|(x, _)| x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Consts;

    fn gold_p() -> (usize, RoundingMode) {
        (256, RoundingMode::ToEven)
    }

    #[test]
    fn roots_bisect_newton_brent_illinois() -> Result<(), crate::Error> {
        let (p, rm) = gold_p();
        let mut cc = Consts::new()?;
        let tol = root_default_tol(p, rm);
        let three = ExactNum::from_u8(3, p);
        let four = ExactNum::from_u8(4, p);
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);
        let two = ExactNum::from_u8(2, p);
        let pi = cc.pi(p, rm);

        let (broot, biters) = bisect_counted(
            |x, p, rm, cc| x.sin(p, rm, cc),
            ExactNum::from_u8(3, p),
            ExactNum::from_u8(4, p),
            &tol,
            p,
            rm,
            &mut cc,
        )
        .ok_or(crate::Error::InvalidArgument)?;
        assert_eq!(broot.compare(&pi), Some(0));

        let nroot = newton(
            |x, p, rm, _cc| x.sqr(p, rm).sub(&two, p, rm),
            |x, p, rm, _cc| ExactNum::mul_ref(&two, x, p, rm),
            one,
            NewtonOptions {
                tol: &tol,
                max_iter: ROOT_MAX_ITER,
            },
            p,
            rm,
            &mut cc,
        )
        .ok_or(crate::Error::InvalidArgument)?;
        let s2 = two.sqrt(p, rm);
        assert_eq!(nroot.compare(&s2), Some(0));

        let (rbrent, briters) = brent_counted(
            |x, p, rm, cc| x.sin(p, rm, cc),
            ExactNum::from_u8(3, p),
            ExactNum::from_u8(4, p),
            &tol,
            p,
            rm,
            &mut cc,
        )
        .ok_or(crate::Error::InvalidArgument)?;
        assert_eq!(rbrent.compare(&pi), Some(0));
        assert!(briters < biters);

        let iroot = illinois(
            |x, p, rm, cc| x.sin(p, rm, cc),
            three,
            four,
            &tol,
            p,
            rm,
            &mut cc,
        )
        .ok_or(crate::Error::InvalidArgument)?;
        assert_eq!(iroot.compare(&pi), Some(0));

        assert!(bisect(
            |x, p, rm, cc| x.sin(p, rm, cc),
            zero,
            ExactNum::from_u8(1, p),
            &tol,
            p,
            rm,
            &mut cc
        )
        .is_none());
        Ok(())
    }
}
