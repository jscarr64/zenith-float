//! Complex Gaussian \({}_2F_1(a,b;c;z)\).
//!
//! Series, Euler, Pfaff, and Kummer / linear transformations in \(\mathbb{C}\).
//! Not the real series at \(\lvert z\rvert\).

use crate::common::util::round_p;
use crate::complex_special::is_nonpos_integer;
use crate::complex_special::nan_pair;
use crate::complex_special::neg_c;
use crate::complex_special::term_negligible;
use crate::complex_special::ziv_complex;
use crate::Consts;
use crate::Error;
use crate::ExactComplex;
use crate::ExactNum;
use crate::RoundingMode;

/// Series term cap (named; same order as the real kernel).
const HYPERGEOM_SERIES_MAX_TERMS: u32 = 10_000;

/// Linear-transform attempts before `NaN`.
const HYPERGEOM_TRANSFORM_MAX: u32 = 8;

fn rm() -> RoundingMode {
    RoundingMode::None
}

fn is_c_zero(z: &ExactComplex) -> bool {
    z.re().is_zero() && z.im().is_zero()
}

fn is_c_one(z: &ExactComplex, p: usize) -> bool {
    z.im().is_zero() && z.re().cmp(&ExactNum::from_u8(1, p)) == Some(0)
}

fn c_eq(a: &ExactComplex, b: &ExactComplex) -> bool {
    a.re().cmp(b.re()) == Some(0) && a.im().cmp(b.im()) == Some(0)
}

fn abs_lt_one(z: &ExactComplex, p: usize) -> bool {
    let az = z.abs(p, rm());
    matches!(az.cmp(&ExactNum::from_u8(1, p)), Some(c) if c < 0)
}

fn re_lt_half(z: &ExactComplex, p: usize) -> bool {
    let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, rm());
    matches!(z.re().cmp(&half), Some(c) if c < 0)
}

fn re_positive_strict(z: &ExactComplex) -> bool {
    z.re().is_positive() && !z.re().is_zero()
}

fn any_nan(args: &[&ExactComplex]) -> bool {
    args.iter().any(|z| z.is_nan())
}

fn as_i32_real_int(z: &ExactComplex, p: usize) -> Option<i32> {
    if !z.im().is_zero() || !z.re().is_int() {
        return None;
    }
    if z.re().is_zero() {
        return Some(0);
    }
    for n in 1i32..=10_000 {
        let w = ExactNum::from_u32(n as u32, p);
        if z.re().abs().cmp(&w) == Some(0) {
            return Some(if z.re().is_negative() { -n } else { n });
        }
    }
    None
}

fn terminating_neg_int(v: &ExactComplex, p: usize) -> bool {
    as_i32_real_int(v, p).is_some_and(|n| n <= 0)
}

fn terminating_after(v: &ExactComplex, next_n: u32, p: usize) -> bool {
    as_i32_real_int(v, p).is_some_and(|m| m <= 0 && next_n as i32 > -m)
}

fn pole_c(c: &ExactComplex, a: &ExactComplex, b: &ExactComplex, p: usize) -> bool {
    let Some(cn) = as_i32_real_int(c, p) else {
        return false;
    };
    if cn > 0 {
        return false;
    }
    let pole_k = -cn;
    let stop_a = as_i32_real_int(a, p).filter(|&n| n <= 0);
    let stop_b = as_i32_real_int(b, p).filter(|&n| n <= 0);
    match (stop_a, stop_b) {
        (Some(sa), _) if -sa < pole_k => false,
        (_, Some(sb)) if -sb < pole_k => false,
        _ => true,
    }
}

fn gamma_fixed(z: &ExactComplex, p: usize, cc: &mut Consts) -> ExactComplex {
    if is_nonpos_integer(z) {
        return nan_pair(Error::InvalidArgument);
    }
    z.gamma_at(p, cc)
}

fn hypergeom_series(
    a: &ExactComplex,
    b: &ExactComplex,
    c: &ExactComplex,
    z: &ExactComplex,
    p: usize,
) -> ExactComplex {
    let one = ExactComplex::one(p);
    let mut term = ExactComplex::one(p);
    let mut sum = ExactComplex::one(p);
    let tiny = ExactNum::from_u8(1, p).ldexp(-((p as i32) - 8), p, rm());
    // Stop when the term is below 2^-p relative to the sum, past the largest parameter (so the
    // terms are decreasing). Hitting the cap means the series has not converged: NaN rather than
    // a silently truncated value.
    let big = max_abs_param(&[a, b, c], p);
    let mut converged = false;
    for n in 0..HYPERGEOM_SERIES_MAX_TERMS as usize {
        if n > 0 && (n as u64) > big && term_small_rel(&term, &sum, p) {
            converged = true;
            break;
        }
        let nn = ExactComplex::from_real(ExactNum::from_u32(n as u32, p), p);
        let den = c.add(&nn, p, rm()).mul(&one.add(&nn, p, rm()), p, rm());
        if matches!(den.abs(p, rm()).cmp(&tiny), Some(c) if c < 0) {
            return nan_pair(Error::InvalidArgument);
        }
        term = term
            .mul(&a.add(&nn, p, rm()), p, rm())
            .mul(&b.add(&nn, p, rm()), p, rm())
            .div(&den, p, rm())
            .mul(z, p, rm());
        sum = sum.add(&term, p, rm());
        if terminating_after(a, (n + 1) as u32, p) || terminating_after(b, (n + 1) as u32, p) {
            sum.set_inexact(false);
            return sum;
        }
    }
    if !converged {
        return nan_pair(Error::InvalidArgument);
    }
    sum
}

/// Upper bound on `max(|a|, |b|, |c|)` as an integer (saturating).
fn max_abs_param(v: &[&ExactComplex], p: usize) -> u64 {
    v.iter()
        .map(|z| match z.abs(p.min(128), rm()).exponent() {
            Some(e) if e > 0 => 1u64 << (e as u32).min(62),
            _ => 1,
        })
        .max()
        .unwrap_or(1)
}

/// `|t| < 2^-p |s|` (or `t = 0`).
fn term_small_rel(t: &ExactComplex, s: &ExactComplex, p: usize) -> bool {
    let at = t.abs(64, rm());
    if at.is_zero() {
        return true;
    }
    let as_ = s.abs(64, rm());
    match (at.exponent(), as_.exponent()) {
        (Some(et), Some(es)) => (et as i64) + (p as i64) < es as i64,
        _ => term_negligible(t, p),
    }
}

/// Kummer: \({}_2F_1(a,b;c;1)=\Gamma(c)\Gamma(c-a-b)/(\Gamma(c-a)\Gamma(c-b))\)
/// when \(\mathrm{Re}(c-a-b)>0\).
fn kummer_z_one(
    a: &ExactComplex,
    b: &ExactComplex,
    c: &ExactComplex,
    p: usize,
    cc: &mut Consts,
) -> ExactComplex {
    let cab = c.sub(a, p, rm()).sub(b, p, rm());
    if !re_positive_strict(&cab) {
        return nan_pair(Error::InvalidArgument);
    }
    let gc = gamma_fixed(c, p, cc);
    let gcab = gamma_fixed(&cab, p, cc);
    let gca = gamma_fixed(&c.sub(a, p, rm()), p, cc);
    let gcb = gamma_fixed(&c.sub(b, p, rm()), p, cc);
    if gc.is_nan() || gcab.is_nan() || gca.is_nan() || gcb.is_nan() {
        return nan_pair(Error::InvalidArgument);
    }
    gc.mul(&gcab, p, rm()).div(&gca.mul(&gcb, p, rm()), p, rm())
}

/// \({}_2F_1(1,1;2;z)=-\ln(1-z)/z\).
fn is_112(a: &ExactComplex, b: &ExactComplex, c: &ExactComplex, p: usize) -> bool {
    let one = ExactComplex::one(p);
    let two = ExactComplex::from_real(ExactNum::from_u8(2, p), p);
    c_eq(a, &one) && c_eq(b, &one) && c_eq(c, &two)
}

fn f112(z: &ExactComplex, p: usize, cc: &mut Consts) -> ExactComplex {
    let one = ExactComplex::one(p);
    let omz = one.sub(z, p, rm());
    neg_c(&omz.ln(p, rm(), cc)).div(z, p, rm())
}

fn pfaff(
    a: &ExactComplex,
    b: &ExactComplex,
    c: &ExactComplex,
    z: &ExactComplex,
    p: usize,
    cc: &mut Consts,
    left: u32,
) -> ExactComplex {
    let one = ExactComplex::one(p);
    let zm1 = z.sub(&one, p, rm());
    if is_c_zero(&zm1) {
        return nan_pair(Error::InvalidArgument);
    }
    let w = z.div(&zm1, p, rm());
    let pref = one.sub(z, p, rm()).pow(&neg_c(a), p, rm(), cc);
    let cb = c.sub(b, p, rm());
    pref.mul(&hypergeom_at(a, &cb, c, &w, p, cc, left), p, rm())
}

/// A&S 15.3.6: argument \(1-z\).
fn transform_1mz(
    a: &ExactComplex,
    b: &ExactComplex,
    c: &ExactComplex,
    z: &ExactComplex,
    p: usize,
    cc: &mut Consts,
    left: u32,
) -> ExactComplex {
    let one = ExactComplex::one(p);
    let omz = one.sub(z, p, rm());
    let cab = c.sub(a, p, rm()).sub(b, p, rm());
    if is_nonpos_integer(&cab) || is_nonpos_integer(&neg_c(&cab)) {
        return nan_pair(Error::InvalidArgument);
    }
    let gc = gamma_fixed(c, p, cc);
    let gcab = gamma_fixed(&cab, p, cc);
    let gca = gamma_fixed(&c.sub(a, p, rm()), p, cc);
    let gcb = gamma_fixed(&c.sub(b, p, rm()), p, cc);
    let abc = a.add(b, p, rm()).sub(c, p, rm());
    let gabc = gamma_fixed(&abc, p, cc);
    let ga = gamma_fixed(a, p, cc);
    let gb = gamma_fixed(b, p, cc);
    if gc.is_nan()
        || gcab.is_nan()
        || gca.is_nan()
        || gcb.is_nan()
        || gabc.is_nan()
        || ga.is_nan()
        || gb.is_nan()
    {
        return nan_pair(Error::InvalidArgument);
    }
    let a_pref = gc.mul(&gcab, p, rm()).div(&gca.mul(&gcb, p, rm()), p, rm());
    let b_pref = gc.mul(&gabc, p, rm()).div(&ga.mul(&gb, p, rm()), p, rm());
    let cap1 = a.add(b, p, rm()).sub(c, p, rm()).add(&one, p, rm());
    let t1 = a_pref.mul(&hypergeom_at(a, b, &cap1, &omz, p, cc, left), p, rm());
    let f2 = hypergeom_at(
        &c.sub(a, p, rm()),
        &c.sub(b, p, rm()),
        &cab.add(&one, p, rm()),
        &omz,
        p,
        cc,
        left,
    );
    let t2 = b_pref
        .mul(&omz.pow(&cab, p, rm(), cc), p, rm())
        .mul(&f2, p, rm());
    t1.add(&t2, p, rm())
}

/// `1/Γ(z)`, exactly zero at the poles.
fn rgamma_c(z: &ExactComplex, p: usize, cc: &mut Consts) -> ExactComplex {
    if is_nonpos_integer(z) {
        return ExactComplex::zero(p);
    }
    ExactComplex::one(p).div(&z.gamma_at(p, cc), p, rm())
}

/// `z^n` for `n ≥ 0` by repeated multiplication (no branch cut).
fn powi_c(z: &ExactComplex, n: u32, p: usize) -> ExactComplex {
    let mut acc = ExactComplex::one(p);
    for _ in 0..n {
        acc = acc.mul(z, p, rm());
    }
    acc
}

/// Binary exponent of `|z|` (`i32::MIN` for zero / non-finite).
fn cexp(z: &ExactComplex) -> i32 {
    match z.abs(64, rm()).exponent() {
        Some(e) if !z.is_nan() => e,
        _ => i32::MIN,
    }
}

/// Integer \(m=c-a-b\) (real, \(\lvert m\rvert\le 10\,000\)), if it is one.
fn int_cab(a: &ExactComplex, b: &ExactComplex, c: &ExactComplex, p: usize) -> Option<i32> {
    let cab = c.sub(a, p, rm()).sub(b, p, rm());
    as_i32_real_int(&cab, p)
}

/// Degenerate \(1-z\) connection, \(c=a+b+m\) with integer \(m\ge 0\) and \(\lvert 1-z\rvert<1\)
/// (A&S 15.3.10–15.3.11, DLMF 15.8.10). With \(w=1-z\):
///
/// \[F=\frac{\Gamma(m)\Gamma(c)}{\Gamma(a+m)\Gamma(b+m)}\sum_{n<m}\frac{(a)_n(b)_n}{n!(1-m)_n}w^n
///   -(-w)^m\frac{\Gamma(c)}{\Gamma(a)\Gamma(b)}\sum_{n\ge0}\frac{(a+m)_n(b+m)_n}{n!(n+m)!}w^n
///   \bigl[\ln w-\psi(n+1)-\psi(n+m+1)+\psi(a+m+n)+\psi(b+m+n)\bigr].\]
///
/// Returns the value and the largest binary exponent among the summed parts (for the loss
/// estimate). `None` if the series does not converge within the term cap.
fn log_case_sum(
    a: &ExactComplex,
    b: &ExactComplex,
    c: &ExactComplex,
    z: &ExactComplex,
    m: u32,
    p: usize,
    cc: &mut Consts,
) -> Option<(ExactComplex, i32)> {
    let one = ExactComplex::one(p);
    let w = one.sub(z, p, rm());
    let re = |v: u32| ExactComplex::from_real(ExactNum::from_u32(v, p), p);
    let mf = re(m);
    let gc = c.gamma_at(p, cc);
    let mut max_e = i32::MIN;
    let mut total = ExactComplex::zero(p);
    if m >= 1 {
        // Finite part. (1-m)_n ≠ 0 for n < m; Γ(a+m) = Γ(c-b), Γ(b+m) = Γ(c-a).
        let pref = factorial_c(m - 1, p)
            .mul(&gc, p, rm())
            .mul(&rgamma_c(&a.add(&mf, p, rm()), p, cc), p, rm())
            .mul(&rgamma_c(&b.add(&mf, p, rm()), p, cc), p, rm());
        let mut t = ExactComplex::one(p);
        let mut s = ExactComplex::one(p);
        let one_m = ExactComplex::from_real(ExactNum::from_i64(1 - m as i64, p), p);
        for n in 0..m - 1 {
            let nn = re(n);
            let num = a.add(&nn, p, rm()).mul(&b.add(&nn, p, rm()), p, rm());
            let den = re(n + 1).mul(&one_m.add(&nn, p, rm()), p, rm());
            t = t.mul(&num, p, rm()).div(&den, p, rm()).mul(&w, p, rm());
            max_e = max_e.max(cexp(&t.mul(&pref, p, rm())));
            s = s.add(&t, p, rm());
        }
        let part = pref.mul(&s, p, rm());
        max_e = max_e.max(cexp(&pref));
        total = part;
    }
    // Logarithmic part; skipped when 1/(Γ(a)Γ(b)) = 0.
    let rg = rgamma_c(a, p, cc).mul(&rgamma_c(b, p, cc), p, rm());
    if is_c_zero(&rg) {
        return Some((total, max_e));
    }
    let mut pref2 = powi_c(&neg_c(&w), m, p).mul(&gc, p, rm()).mul(&rg, p, rm());
    pref2 = neg_c(&pref2);
    let am = a.add(&mf, p, rm());
    let bm = b.add(&mf, p, rm());
    let lnw = w.ln(p, rm(), cc);
    let mut psi_n1 = one.digamma_at(p, cc);
    let mut psi_nm1 = re(m + 1).digamma_at(p, cc);
    let mut psi_a = am.digamma_at(p, cc);
    let mut psi_b = bm.digamma_at(p, cc);
    let mut coef = ExactComplex::one(p).div(&factorial_c(m, p), p, rm());
    let mut wn = ExactComplex::one(p);
    let mut sum = ExactComplex::zero(p);
    let big = max_abs_param(&[a, b, c], p);
    let mut converged = false;
    for n in 0..HYPERGEOM_SERIES_MAX_TERMS {
        let bracket = lnw
            .sub(&psi_n1, p, rm())
            .sub(&psi_nm1, p, rm())
            .add(&psi_a, p, rm())
            .add(&psi_b, p, rm());
        let term = coef.mul(&wn, p, rm()).mul(&bracket, p, rm());
        sum = sum.add(&term, p, rm());
        max_e = max_e.max(cexp(&term.mul(&pref2, p, rm())));
        if (n as u64) > big && n > 0 && term_small_rel(&term, &sum, p) {
            converged = true;
            break;
        }
        // Advance n → n+1.
        let nn = re(n);
        let a_n = am.add(&nn, p, rm());
        let b_n = bm.add(&nn, p, rm());
        coef = coef.mul(&a_n.mul(&b_n, p, rm()), p, rm()).div(
            &re(n + 1).mul(&re(n + m + 1), p, rm()),
            p,
            rm(),
        );
        wn = wn.mul(&w, p, rm());
        psi_n1 = psi_n1.add(&one.div(&re(n + 1), p, rm()), p, rm());
        psi_nm1 = psi_nm1.add(&one.div(&re(n + m + 1), p, rm()), p, rm());
        psi_a = psi_a.add(&one.div(&a_n, p, rm()), p, rm());
        psi_b = psi_b.add(&one.div(&b_n, p, rm()), p, rm());
    }
    if !converged {
        return None;
    }
    let v = total.add(&pref2.mul(&sum, p, rm()), p, rm());
    Some((v, max_e))
}

fn factorial_c(n: u32, p: usize) -> ExactComplex {
    let mut acc = ExactNum::from_u8(1, p);
    for k in 2..=n {
        acc = acc.mul(&ExactNum::from_u32(k, p), p, rm());
    }
    ExactComplex::from_real(acc, p)
}

/// \(1-z\) transform for integer \(m=c-a-b\), \(\lvert 1-z\rvert<1\). Negative \(m\) goes through
/// Euler's transformation \(F(a,b;c;z)=(1-z)^{m}F(c-a,c-b;c;z)\) (then \(c-a'-b'=-m>0\)).
/// Cancellation between the parts is measured and the sum re-evaluated with that many extra bits.
fn transform_1mz_log(
    a: &ExactComplex,
    b: &ExactComplex,
    c: &ExactComplex,
    z: &ExactComplex,
    m: i32,
    p: usize,
    cc: &mut Consts,
    left: u32,
) -> ExactComplex {
    if m < 0 {
        let one = ExactComplex::one(p);
        let w = one.sub(z, p, rm());
        let wm = one.div(&powi_c(&w, m.unsigned_abs(), p), p, rm());
        let f = hypergeom_at(&c.sub(a, p, rm()), &c.sub(b, p, rm()), c, z, p, cc, left);
        return wm.mul(&f, p, rm());
    }
    let m = m as u32;
    let mut pw = p + 32;
    for _ in 0..3 {
        let Some((v, max_e)) = log_case_sum(a, b, c, z, m, pw, cc) else {
            return nan_pair(Error::InvalidArgument);
        };
        let ve = cexp(&v);
        if ve == i32::MIN {
            return v;
        }
        let lost = max_e.saturating_sub(ve).max(0) as usize;
        if lost + 16 <= pw - p {
            return v;
        }
        pw = p + lost + 48;
    }
    nan_pair(Error::PrecisionRetryExhausted)
}

/// Bits lost to cancellation in the generic \(1-z\) transform: \(\Gamma(\pm(c-a-b))\) grow like
/// \(1/\delta\), \(\delta\) = distance of \(c-a-b\) to the nearest integer, and the two terms cancel
/// to \(O(1)\). Returns \(\lceil\log_2(1/\delta)\rceil\) (0 when \(\delta\ge 1/2\)).
fn near_int_loss(a: &ExactComplex, b: &ExactComplex, c: &ExactComplex, p: usize) -> usize {
    let cab = c.sub(a, p, rm()).sub(b, p, rm());
    let e = |x: &ExactNum| -> i64 {
        if x.is_zero() {
            i64::MIN / 2
        } else {
            x.exponent().map_or(0, i64::from)
        }
    };
    let f = cab.re().fract().abs();
    let g = ExactNum::from_u8(1, p).sub(&f, p, rm());
    let re_e = e(&f).min(e(&g));
    let dist_e = re_e.max(e(cab.im()));
    usize::try_from(-dist_e).unwrap_or(0)
}

/// Generic \(1-z\) transform with guard bits for \(c-a-b\) close to an integer.
fn transform_1mz_guarded(
    a: &ExactComplex,
    b: &ExactComplex,
    c: &ExactComplex,
    z: &ExactComplex,
    p: usize,
    cc: &mut Consts,
    left: u32,
) -> ExactComplex {
    let loss = near_int_loss(a, b, c, p);
    if loss > 8 * p + 1024 {
        return nan_pair(Error::PrecisionRetryExhausted);
    }
    let pw = round_p(p + loss + 32);
    transform_1mz(a, b, c, z, pw, cc, left)
}

/// A&S 15.3.7: argument \(1/z\). Degenerate when \(a=b\).
fn transform_inv(
    a: &ExactComplex,
    b: &ExactComplex,
    c: &ExactComplex,
    z: &ExactComplex,
    p: usize,
    cc: &mut Consts,
    left: u32,
) -> ExactComplex {
    if c_eq(a, b) {
        return nan_pair(Error::InvalidArgument);
    }
    let one = ExactComplex::one(p);
    let inv = one.div(z, p, rm());
    let mz = neg_c(z);
    let t1 = inv_term(a, b, c, &mz, &inv, p, cc, left);
    let t2 = inv_term(b, a, c, &mz, &inv, p, cc, left);
    t1.add(&t2, p, rm())
}

fn inv_term(
    a: &ExactComplex,
    b: &ExactComplex,
    c: &ExactComplex,
    mz: &ExactComplex,
    inv: &ExactComplex,
    p: usize,
    cc: &mut Consts,
    left: u32,
) -> ExactComplex {
    let one = ExactComplex::one(p);
    let gc = gamma_fixed(c, p, cc);
    let gbma = gamma_fixed(&b.sub(a, p, rm()), p, cc);
    let gb = gamma_fixed(b, p, cc);
    let gcma = gamma_fixed(&c.sub(a, p, rm()), p, cc);
    if gc.is_nan() || gbma.is_nan() || gb.is_nan() || gcma.is_nan() {
        return nan_pair(Error::InvalidArgument);
    }
    let pref = gc.mul(&gbma, p, rm()).div(&gb.mul(&gcma, p, rm()), p, rm());
    let pow = mz.pow(&neg_c(a), p, rm(), cc);
    let ac1 = a.sub(c, p, rm()).add(&one, p, rm());
    let ab1 = a.sub(b, p, rm()).add(&one, p, rm());
    pref.mul(&pow, p, rm())
        .mul(&hypergeom_at(a, &ac1, &ab1, inv, p, cc, left), p, rm())
}

fn hypergeom_at(
    a: &ExactComplex,
    b: &ExactComplex,
    c: &ExactComplex,
    z: &ExactComplex,
    p: usize,
    cc: &mut Consts,
    left: u32,
) -> ExactComplex {
    if any_nan(&[a, b, c, z]) {
        return nan_pair(Error::InvalidArgument);
    }
    if is_c_zero(z) || is_c_zero(a) || is_c_zero(b) {
        return ExactComplex::one(p);
    }
    if pole_c(c, a, b, p) {
        return nan_pair(Error::InvalidArgument);
    }
    if is_c_one(z, p) {
        return kummer_z_one(a, b, c, p, cc);
    }
    if is_112(a, b, c, p) {
        return f112(z, p, cc);
    }
    if terminating_neg_int(a, p) || terminating_neg_int(b, p) {
        return hypergeom_series(a, b, c, z, p);
    }
    // Integer c − a − b near z = 1: the direct series converges too slowly (or not at all for
    // |z| ≥ 1) and the generic 1 − z transform is singular, so use the logarithmic connection.
    if let Some(m) = int_cab(a, b, c, p) {
        let omz = ExactComplex::one(p).sub(z, p, rm());
        let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, rm());
        let near = matches!(omz.abs(p, rm()).cmp(&half), Some(o) if o <= 0);
        if left > 0 && (near || (abs_lt_one(&omz, p) && !abs_lt_one(z, p))) {
            return transform_1mz_log(a, b, c, z, m, p, cc, left - 1);
        }
    } else if left > 0 {
        // Non-integer c − a − b with z close to 1 inside the unit disc: the series would need
        // ≈ p / log2(1/|z|) terms, so use the 1 − z transform (guarded for near-integer c − a − b).
        let omz = ExactComplex::one(p).sub(z, p, rm());
        let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, rm());
        let nine_tenths = ExactNum::from_u8(9, p).div(&ExactNum::from_u8(10, p), p, rm());
        let near = matches!(omz.abs(p, rm()).cmp(&half), Some(o) if o <= 0);
        let slow = matches!(z.abs(p, rm()).cmp(&nine_tenths), Some(o) if o > 0);
        if near && slow && abs_lt_one(z, p) {
            return transform_1mz_guarded(a, b, c, z, p, cc, left - 1);
        }
    }
    if abs_lt_one(z, p) {
        return hypergeom_series(a, b, c, z, p);
    }
    if left == 0 {
        return nan_pair(Error::InvalidArgument);
    }
    let next = left - 1;
    let one = ExactComplex::one(p);
    let omz = one.sub(z, p, rm());
    let zm1 = z.sub(&one, p, rm());
    if re_lt_half(z, p) && !is_c_zero(&zm1) {
        return pfaff(a, b, c, z, p, cc, next);
    }
    if abs_lt_one(&omz, p) && !is_nonpos_integer(&c.sub(a, p, rm()).sub(b, p, rm())) {
        return transform_1mz_guarded(a, b, c, z, p, cc, next);
    }
    let inv = one.div(z, p, rm());
    if abs_lt_one(&inv, p) && !c_eq(a, b) {
        return transform_inv(a, b, c, z, p, cc, next);
    }
    if !is_c_zero(&zm1) {
        let w = z.div(&zm1, p, rm());
        if abs_lt_one(&w, p) {
            return pfaff(a, b, c, z, p, cc, next);
        }
    }
    nan_pair(Error::InvalidArgument)
}

impl ExactComplex {
    /// Gaussian \({}_2F_1(a=\mathrm{self},b;c;z)\) in \(\mathbb{C}\).
    ///
    /// All-real arguments with \(z<1\) use the real [`ExactNum::hypergeom_2f1`].
    /// Otherwise: series when \(\lvert z\rvert<1\) (converged to \(2^{-p}\) relative, NaN if not
    /// converged within the term cap); Pfaff when \(\mathrm{Re}(z)<1/2\);
    /// Euler / \(1-z\) and \(1/z\) linear transforms otherwise. Kummer at \(z=1\)
    /// when \(\mathrm{Re}(c-a-b)>0\). Cut on \([1,+\infty)\) in \(z\) (principal
    /// value from above). Non-positive integer \(c\) (uncanceled) → NaN.
    ///
    /// Integer \(m=c-a-b\) with \(z\) near 1 (\(\lvert 1-z\rvert\le 1/2\), or
    /// \(\lvert 1-z\rvert<1\) with \(\lvert z\rvert\ge 1\)): logarithmic \(1-z\) connection
    /// (DLMF 15.8.10, A&S 15.3.10–15.3.11) for \(m\ge 0\), Euler's transformation
    /// \((1-z)^m F(c-a,c-b;c;z)\) first for \(m<0\). Near-integer \(c-a-b\) uses the generic
    /// \(1-z\) transform with \(\log_2(1/\delta)\) guard bits (\(\delta\) = distance to the
    /// nearest integer).
    ///
    /// # Precision
    ///
    /// - Algorithm: series / Euler / Pfaff / Kummer / logarithmic \(1-z\) connection. Caps `HYPERGEOM_SERIES_MAX_TERMS = 10_000`, `HYPERGEOM_TRANSFORM_MAX = 8`.
    /// - Bound: Ziv on each part (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    pub fn hypergeom_2f1(
        &self,
        b: &Self,
        c: &Self,
        z: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        if any_nan(&[self, b, c, z]) {
            return nan_pair(Error::InvalidArgument);
        }
        let dest = round_p(p);
        if is_c_zero(z) || is_c_zero(self) || is_c_zero(b) {
            let mut one = ExactComplex::one(dest);
            one.set_inexact(self.inexact() | b.inexact() | c.inexact() | z.inexact());
            return one;
        }
        if pole_c(c, self, b, dest) {
            return nan_pair(Error::InvalidArgument);
        }
        // All-real arguments with z < 1: the real kernel handles |z| → 1 (including the
        // logarithmic integer-(c − a − b) case) and certifies its rounding.
        if self.im().is_zero()
            && b.im().is_zero()
            && c.im().is_zero()
            && z.im().is_zero()
            && matches!(z.re().cmp(&ExactNum::from_u8(1, dest)), Some(o) if o < 0)
        {
            let r = self
                .re()
                .hypergeom_2f1(b.re(), c.re(), z.re(), dest, rm, cc);
            if !r.is_nan() {
                return ExactComplex::new(r, ExactNum::new(dest));
            }
        }
        ziv_complex(dest, rm, |pw| {
            hypergeom_at(self, b, c, z, pw, cc, HYPERGEOM_TRANSFORM_MAX)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(a: &ExactNum, b: &ExactNum, p: usize, slack: i32) -> bool {
        let d = a.sub(b, p, RoundingMode::None).abs();
        d.is_zero() || d.exponent().unwrap_or(0) < -((p as i32) - slack)
    }

    fn tiny(x: &ExactNum, p: usize) -> bool {
        x.is_zero() || x.exponent().unwrap_or(0) < -((p as i32) / 4)
    }

    #[test]
    fn complex_hypergeom_2f1_golds() {
        let p = 256;
        let r = RoundingMode::ToEven;
        let mut cc = Consts::new().unwrap();
        let one = ExactComplex::one(p);
        let two = ExactComplex::from_real(ExactNum::from_u8(2, p), p);
        let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, r);
        let zh = ExactComplex::from_real(half.clone(), p);
        let a_h = zh.clone();
        let z0 = ExactComplex::zero(p);

        let f1 = one.hypergeom_2f1(&one, &two, &zh, p, r, &mut cc);
        let ln2 = cc.ln_2(p, r);
        let want = ln2.mul(&ExactNum::from_u8(2, p), p, r);
        assert!(near(f1.re(), &want, p, 40), "2F1(1,1;2;1/2)=2ln2");
        assert!(tiny(f1.im(), p));

        let fk = a_h.hypergeom_2f1(&a_h, &one, &zh, p, r, &mut cc);
        let k = half.elliptic_k(p, r, &mut cc);
        let pi = cc.pi(p, r);
        let want_k = k.mul(&ExactNum::from_u8(2, p), p, r).div(&pi, p, r);
        assert!(near(fk.re(), &want_k, p, 40), "2F1(1/2,1/2;1;1/2)=2K/π");
        assert!(tiny(fk.im(), p));

        let a = ExactComplex::new(
            ExactNum::from_u8(2, p).div(&ExactNum::from_u8(3, p), p, r),
            ExactNum::from_u8(1, p).div(&ExactNum::from_u8(4, p), p, r),
        );
        let b = ExactComplex::new(
            ExactNum::from_u8(1, p).div(&ExactNum::from_u8(5, p), p, r),
            ExactNum::from_u8(1, p).div(&ExactNum::from_u8(7, p), p, r),
        );
        let c = ExactComplex::new(
            ExactNum::from_u8(3, p).div(&ExactNum::from_u8(2, p), p, r),
            ExactNum::from_u8(1, p).div(&ExactNum::from_u8(9, p), p, r),
        );
        let f0 = a.hypergeom_2f1(&b, &c, &z0, p, r, &mut cc);
        assert!(
            near(f0.re(), &ExactNum::from_u8(1, p), p, 8),
            "2F1(*,*,*;0)=1"
        );
        assert!(tiny(f0.im(), p));

        let z = ExactComplex::new(
            ExactNum::from_u8(3, p).div(&ExactNum::from_u8(10, p), p, r),
            ExactNum::from_u8(1, p).div(&ExactNum::from_u8(10, p), p, r),
        );
        let lhs = a.hypergeom_2f1(&b, &c, &z, p, r, &mut cc);
        let cab = c.sub(&a, p, r).sub(&b, p, r);
        let pref = ExactComplex::one(p).sub(&z, p, r).pow(&cab, p, r, &mut cc);
        let rhs = pref.mul(
            &c.sub(&a, p, r)
                .hypergeom_2f1(&c.sub(&b, p, r), &c, &z, p, r, &mut cc),
            p,
            r,
        );
        assert!(near(lhs.re(), rhs.re(), p, 30), "Euler re");
        assert!(near(lhs.im(), rhs.im(), p, 30), "Euler im");

        let c0 = ExactComplex::zero(p);
        let bad = one.hypergeom_2f1(&one, &c0, &zh, p, r, &mut cc);
        assert!(bad.is_nan(), "c=0 → NaN");

        let two_r = ExactNum::from_u8(2, p);
        let eps = ExactNum::from_u8(1, p).ldexp(-40, p, RoundingMode::None);
        let above = ExactComplex::new(two_r.clone(), eps.clone());
        let below = ExactComplex::new(two_r, eps.neg());
        let fa = one.hypergeom_2f1(&one, &two, &above, p, r, &mut cc);
        let fb = one.hypergeom_2f1(&one, &two, &below, p, r, &mut cc);
        assert!(!fa.is_nan() && !fb.is_nan(), "cut defined");
        assert!(near(fa.re(), fb.re(), p, 20), "cut Re");
        assert_ne!(fa.im().cmp(fb.im()), Some(0), "cut Im differs");
        assert!(near(&fa.im().abs(), &fb.im().abs(), p, 20), "cut |Im|");
    }
}
