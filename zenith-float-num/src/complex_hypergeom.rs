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
use crate::defs::Exponent;
use crate::Consts;
use crate::Error;
use crate::ExactComplex;
use crate::ExactNum;
use crate::RoundingMode;
use crate::WORD_BIT_SIZE;

/// Linear-transform attempts before `NaN`.
const HYPERGEOM_TRANSFORM_MAX: u32 = 8;

fn rm() -> RoundingMode {
    RoundingMode::None
}

/// `1 − z` computed as `−(z − 1)`, so a signed zero imaginary part of `z` keeps its side of the
/// cut: `z = 2 + 0i` (the limit from above) gives `−1 − 0i`, as the limit `1 − (2 + iε)` does.
/// `ExactComplex::one(p).sub(z, ..)` gives `0 − 0 = +0` for both signs and loses the side.
fn one_minus(z: &ExactComplex, p: usize) -> ExactComplex {
    neg_c(ExactComplex::sub_ref(z, &ExactComplex::one(p), p, rm()))
}

fn is_c_zero(z: &ExactComplex) -> bool {
    z.re().is_zero() && z.im().is_zero()
}

fn is_c_one(z: &ExactComplex, p: usize) -> bool {
    z.im().is_zero() && z.re().compare(&ExactNum::from_u8(1, p)) == Some(0)
}

fn c_eq(a: &ExactComplex, b: &ExactComplex) -> bool {
    a.re().compare(b.re()) == Some(0) && a.im().compare(b.im()) == Some(0)
}

fn abs_lt_one(z: &ExactComplex, p: usize) -> bool {
    let az = z.abs(p, rm());
    matches!(az.compare(&ExactNum::from_u8(1, p)), Some(c) if c < 0)
}

fn mag64(z: &ExactComplex) -> ExactNum {
    z.abs(64, rm())
}

fn mag_lt(a: &ExactNum, b: &ExactNum) -> bool {
    matches!(a.compare(b), Some(c) if c < 0)
}

/// Terms needed for a geometric series at `|z|` and precision `p`, plus `big` (the
/// largest parameter, so the terms are already decreasing) and a 128-bit guard.
/// For `|z| < 1`, `n ≲ (p + guard) / −log2(|z|) ≤ 2(p + guard) / (1 − |z|)`.
fn series_term_budget(z: &ExactComplex, p: usize, big: u64) -> usize {
    let base = p.saturating_add(128).saturating_add(big as usize).max(32);
    let one = ExactNum::from_u8(1, 64);
    let az = mag64(z);
    if !mag_lt(&az, &one) {
        return base;
    }
    let mut d = ExactNum::sub_ref(&one, &az, 64, RoundingMode::ToZero);
    if d.is_zero() {
        return base;
    }
    let mut shift = 0u32;
    while shift < 62 && mag_lt(&d, &one) {
        d = d.ldexp(1, 64, RoundingMode::ToZero);
        shift += 1;
    }
    let scale = 1usize.checked_shl(shift).unwrap_or(usize::MAX);
    base.saturating_mul(2)
        .saturating_mul(scale)
        .min(u32::MAX as usize)
}

fn re_positive_strict(z: &ExactComplex) -> bool {
    z.re().is_positive() && !z.re().is_zero()
}

fn any_nan(args: &[&ExactComplex]) -> bool {
    args.iter().any(|z| z.is_nan())
}

fn as_i32_real_int(z: &ExactComplex, p: usize) -> Option<i32> {
    let _ = p;
    if !z.im().is_zero() || !z.re().is_int() {
        return None;
    }
    z.re().to_i32(RoundingMode::ToZero).ok()
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
    let mut term = ExactComplex::one(p);
    let mut sum = ExactComplex::one(p);
    let tiny = ExactNum::from_u8(1, p).ldexp(-((p as Exponent) - 8), p, rm());
    // Stop when the term is below 2^-p relative to the sum, past the largest parameter (so the
    // terms are decreasing). Hitting the cap means the series has not converged: NaN rather than
    // a silently truncated value.
    let big = max_abs_param(&[a, b, c], p);
    let budget = series_term_budget(z, p, big);
    let mut converged = false;
    for n in 0..budget {
        if n > 0 && (n as u64) > big && term_small_rel(&term, &sum, p) {
            converged = true;
            break;
        }
        let nn = ExactComplex::from_real(ExactNum::from_u32(n as u32, p), p);
        let n_plus_c = int_plus(n as u32, c, p);
        let n_plus_1 = trim_c(ExactComplex::one(p).add(&nn, p, rm()));
        let den = trim_c(n_plus_c.mul(&n_plus_1, p, rm()));
        if matches!(den.abs(p, rm()).compare(&tiny), Some(c) if c < 0) {
            return nan_pair(Error::InvalidArgument);
        }
        term = term
            .mul(&int_plus(n as u32, a, p), p, rm())
            .mul(&int_plus(n as u32, b, p), p, rm())
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

/// `x` stored without the all-zero low words of its mantissa: the same value, sign and accumulated
/// flag (the per-operation flag falls back to it), so every later operation gives the same result, but a product or quotient with it costs by its
/// significant words instead of by the working precision. Zero, Inf and NaN come back unchanged.
fn trim_zero_words(mut x: ExactNum) -> ExactNum {
    let Some((words, ..)) = x.value_parts() else {
        return x;
    };
    let zeros = words.iter().take_while(|w| **w == 0).count();
    if zeros == 0 || zeros >= words.len() {
        return x;
    }
    let q = (words.len() - zeros) * WORD_BIT_SIZE;
    match x.set_precision(q, RoundingMode::ToZero) {
        Ok(()) => x,
        Err(e) => ExactNum::nan(Some(e)),
    }
}

/// Both parts of `z` through [`trim_zero_words`].
fn trim_c(z: ExactComplex) -> ExactComplex {
    let (re, im) = z.into_parts();
    ExactComplex::new(trim_zero_words(re), trim_zero_words(im))
}

/// `n + a` at `p` bits, trimmed ([`trim_c`]). For the short parameters of a series (and small `n`)
/// the sum is exact in a few words, which keeps the per-term products and quotient of the
/// `p`-bit term cheap at a high working precision.
fn int_plus(n: u32, a: &ExactComplex, p: usize) -> ExactComplex {
    trim_c(ExactComplex::from_real(ExactNum::from_u32(n, p), p).add(a, p, rm()))
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

/// Bits of each part [`abs_lead`] keeps.
const ABS_LEAD_BITS: usize = 128;

/// `|z|` at 64 bits from the parts cut (toward zero) to their leading [`ABS_LEAD_BITS`] bits:
/// within a relative `2^-126` below `|z|`. Only its exponent is used, and that matches the one of
/// `|z|` unless `|z|` lies within `2^-126` above a power of two. `hypot` reads every input bit, so
/// on the full-width term and sum this check cost more than the series step it guards.
fn abs_lead(z: &ExactComplex) -> ExactNum {
    // `x + 0` at ABS_LEAD_BITS toward zero is `x` cut to its leading bits, as a new value.
    let lead = |x: &ExactNum| {
        ExactNum::add_ref(
            x,
            &ExactNum::new(ABS_LEAD_BITS),
            ABS_LEAD_BITS,
            RoundingMode::ToZero,
        )
    };
    ExactComplex::new(lead(z.re()), lead(z.im())).abs(64, rm())
}

/// `|t| < 2^-p |s|` (or `t = 0`).
fn term_small_rel(t: &ExactComplex, s: &ExactComplex, p: usize) -> bool {
    let at = abs_lead(t);
    if at.is_zero() {
        return true;
    }
    let as_ = abs_lead(s);
    match (at.exponent(), as_.exponent()) {
        (Some(et), Some(es)) => et + (p as i64) < es,
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
    let cab = ExactComplex::sub_ref(c, a, p, rm()).sub(b, p, rm());
    if !re_positive_strict(&cab) {
        return nan_pair(Error::InvalidArgument);
    }
    let gc = gamma_fixed(c, p, cc);
    let gcab = gamma_fixed(&cab, p, cc);
    if gc.is_nan() || gcab.is_nan() {
        return nan_pair(Error::InvalidArgument);
    }
    let rgca = rgamma_c(&ExactComplex::sub_ref(c, a, p, rm()), p, cc);
    let rgcb = rgamma_c(&ExactComplex::sub_ref(c, b, p, rm()), p, cc);
    gc.mul(&gcab, p, rm())
        .mul(&rgca.mul(&rgcb, p, rm()), p, rm())
}

/// \({}_2F_1(1,1;2;z)=-\ln(1-z)/z\).
fn is_112(a: &ExactComplex, b: &ExactComplex, c: &ExactComplex, p: usize) -> bool {
    let one = ExactComplex::one(p);
    let two = ExactComplex::from_real(ExactNum::from_u8(2, p), p);
    c_eq(a, &one) && c_eq(b, &one) && c_eq(c, &two)
}

fn f112(z: &ExactComplex, p: usize, cc: &mut Consts) -> ExactComplex {
    let omz = one_minus(z, p);
    neg_c(omz.ln(p, rm(), cc)).div(z, p, rm())
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
    let zm1 = ExactComplex::sub_ref(z, &one, p, rm());
    if is_c_zero(&zm1) {
        return nan_pair(Error::InvalidArgument);
    }
    let w = z.div(&zm1, p, rm());
    let pref = one_minus(z, p).pow(&a.neg_via_mul(), p, rm(), cc);
    let cb = ExactComplex::sub_ref(c, b, p, rm());
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
    let omz = one_minus(z, p);
    let cab = ExactComplex::sub_ref(c, a, p, rm()).sub(b, p, rm());
    if is_nonpos_integer(&cab) || is_nonpos_integer(&cab.neg_via_mul()) {
        return nan_pair(Error::InvalidArgument);
    }
    let gc = gamma_fixed(c, p, cc);
    let gcab = gamma_fixed(&cab, p, cc);
    let abc = ExactComplex::add_ref(a, b, p, rm()).sub(c, p, rm());
    let gabc = gamma_fixed(&abc, p, cc);
    if gc.is_nan() || gcab.is_nan() || gabc.is_nan() {
        return nan_pair(Error::InvalidArgument);
    }
    // Denominator Γs as 1/Γ: a pole there (c − a, c − b, a or b a non-positive integer) makes
    // that term exactly zero, not NaN.
    let rgca = rgamma_c(&ExactComplex::sub_ref(c, a, p, rm()), p, cc);
    let rgcb = rgamma_c(&ExactComplex::sub_ref(c, b, p, rm()), p, cc);
    let rga = rgamma_c(a, p, cc);
    let rgb = rgamma_c(b, p, cc);
    let a_pref = gc
        .mul(&gcab, p, rm())
        .mul(&rgca.mul(&rgcb, p, rm()), p, rm());
    let b_pref = gc.mul(&gabc, p, rm()).mul(&rga.mul(&rgb, p, rm()), p, rm());
    let cap1 = ExactComplex::add_ref(a, b, p, rm())
        .sub(c, p, rm())
        .add(&one, p, rm());
    let t1 = a_pref.mul(&hypergeom_at(a, b, &cap1, &omz, p, cc, left), p, rm());
    let f2 = hypergeom_at(
        &ExactComplex::sub_ref(c, a, p, rm()),
        &ExactComplex::sub_ref(c, b, p, rm()),
        &one.add(&cab, p, rm()),
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

/// Binary exponent of `|z|` (`Exponent::MIN` for zero / non-finite).
fn cexp(z: &ExactComplex) -> Exponent {
    match z.abs(64, rm()).exponent() {
        Some(e) if !z.is_nan() => e,
        _ => Exponent::MIN,
    }
}

/// Integer \(m=c-a-b\) (real, \(\lvert m\rvert\le 10\,000\)), if it is one.
fn int_cab(a: &ExactComplex, b: &ExactComplex, c: &ExactComplex, p: usize) -> Option<i32> {
    let cab = ExactComplex::sub_ref(c, a, p, rm()).sub(b, p, rm());
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
) -> Option<(ExactComplex, Exponent)> {
    let one = ExactComplex::one(p);
    let w = one_minus(z, p);
    let re = |v: u32| ExactComplex::from_real(ExactNum::from_u32(v, p), p);
    let mf = re(m);
    let gc = c.gamma_at(p, cc);
    let mut max_e = Exponent::MIN;
    let mut total = ExactComplex::zero(p);
    if m >= 1 {
        // Finite part. (1-m)_n ≠ 0 for n < m; Γ(a+m) = Γ(c-b), Γ(b+m) = Γ(c-a).
        let pref = factorial_c(m - 1, p)
            .mul(&gc, p, rm())
            .mul(&rgamma_c(&re(m).add(a, p, rm()), p, cc), p, rm())
            .mul(&rgamma_c(&re(m).add(b, p, rm()), p, cc), p, rm());
        let mut t = ExactComplex::one(p);
        let mut s = ExactComplex::one(p);
        let one_m = ExactComplex::from_real(ExactNum::from_i64(1 - m as i64, p), p);
        for n in 0..m - 1 {
            let nn = re(n);
            let num = re(n).add(a, p, rm()).mul(&re(n).add(b, p, rm()), p, rm());
            let den = re(n + 1).mul(&nn.add(&one_m, p, rm()), p, rm());
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
    let mut pref2 = powi_c(&w.neg_via_mul(), m, p)
        .mul(&gc, p, rm())
        .mul(&rg, p, rm());
    pref2 = neg_c(pref2);
    let am = re(m).add(a, p, rm());
    let bm = mf.add(b, p, rm());
    let lnw = w.ln(p, rm(), cc);
    let mut psi_n1 = one.digamma_at(p, cc);
    let mut psi_nm1 = re(m + 1).digamma_at(p, cc);
    let mut psi_a = am.digamma_at(p, cc);
    let mut psi_b = bm.digamma_at(p, cc);
    let mut coef = ExactComplex::one(p).div(&factorial_c(m, p), p, rm());
    let mut wn = ExactComplex::one(p);
    let mut sum = ExactComplex::zero(p);
    let big = max_abs_param(&[a, b, c], p);
    let budget = series_term_budget(&w, p, big);
    let mut converged = false;
    for n in 0..budget as u32 {
        let bracket = ExactComplex::sub_ref(&lnw, &psi_n1, p, rm())
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
        let a_n = re(n).add(&am, p, rm());
        let b_n = nn.add(&bm, p, rm());
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
        let w = one_minus(z, p);
        let wm = one.div(&powi_c(&w, m.unsigned_abs(), p), p, rm());
        let f = hypergeom_at(
            &ExactComplex::sub_ref(c, a, p, rm()),
            &ExactComplex::sub_ref(c, b, p, rm()),
            c,
            z,
            p,
            cc,
            left,
        );
        return wm.mul(&f, p, rm());
    }
    let m = m as u32;
    let mut pw = p + 32;
    for _ in 0..3 {
        let Some((v, max_e)) = log_case_sum(a, b, c, z, m, pw, cc) else {
            return nan_pair(Error::InvalidArgument);
        };
        let ve = cexp(&v);
        if ve == Exponent::MIN {
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
    let cab = ExactComplex::sub_ref(c, a, p, rm()).sub(b, p, rm());
    let e = |x: &ExactNum| -> i64 {
        if x.is_zero() {
            i64::MIN / 2
        } else {
            x.exponent().map_or(0, i64::from)
        }
    };
    let (re, im) = cab.into_parts();
    let f = re.fract().abs();
    let g = ExactNum::from_u8(1, p).sub(&f, p, rm());
    let re_e = e(&f).min(e(&g));
    let dist_e = re_e.max(e(&im));
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
    let mz = z.neg_via_mul();
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
    let gbma = gamma_fixed(&ExactComplex::sub_ref(b, a, p, rm()), p, cc);
    if gc.is_nan() || gbma.is_nan() {
        return nan_pair(Error::InvalidArgument);
    }
    // 1/Γ(b) and 1/Γ(c − a): zero at a pole (the term vanishes), never NaN.
    let rgb = rgamma_c(b, p, cc);
    let rgcma = rgamma_c(&ExactComplex::sub_ref(c, a, p, rm()), p, cc);
    if is_c_zero(&rgb) || is_c_zero(&rgcma) {
        return ExactComplex::zero(p);
    }
    let pref = gc
        .mul(&gbma, p, rm())
        .mul(&rgb.mul(&rgcma, p, rm()), p, rm());
    let pow = mz.pow(&a.neg_via_mul(), p, rm(), cc);
    let ac1 = ExactComplex::sub_ref(a, c, p, rm()).add(&one, p, rm());
    let ab1 = ExactComplex::sub_ref(a, b, p, rm()).add(&one, p, rm());
    pref.mul(&pow, p, rm())
        .mul(&hypergeom_at(a, &ac1, &ab1, inv, p, cc, left), p, rm())
}

/// `(lo, hi, m)` with `hi = lo + m`, `m` a non-negative integer, if `a − b` is an integer.
fn int_ab<'x>(
    a: &'x ExactComplex,
    b: &'x ExactComplex,
    p: usize,
) -> Option<(&'x ExactComplex, &'x ExactComplex, u32)> {
    let d = as_i32_real_int(&ExactComplex::sub_ref(b, a, p, rm()), p)?;
    if d >= 0 {
        Some((a, b, d.unsigned_abs()))
    } else {
        Some((b, a, d.unsigned_abs()))
    }
}

/// Degenerate \(1/z\) connection, \(b=a+m\) with integer \(m\ge 0\), \(\lvert z\rvert>1\)
/// (DLMF 15.8.8):
///
/// \[\frac{F}{\Gamma(c)}=\frac{(-z)^{-a}}{\Gamma(a+m)}\sum_{k<m}\frac{(a)_k(m-k-1)!}{k!\,\Gamma(c-a-k)}z^{-k}
///   +\frac{(-z)^{-a}}{\Gamma(a)}\sum_{k\ge0}\frac{(a+m)_k(-1)^k z^{-k-m}}{k!(k+m)!}
///   \Bigl[\frac{\ln(-z)+\psi(1+m+k)+\psi(1+k)-\psi(a+m+k)}{\Gamma(c-a-m-k)}-\frac{\psi(c-a-m-k)}{\Gamma(c-a-m-k)}\Bigr].\]
///
/// At \(x=c-a-m-k=-n\) (a non-positive integer) \(1/\Gamma(x)=0\) and \(\psi(x)/\Gamma(x)=(-1)^{n+1}n!\).
/// Returns the value and the largest binary exponent among the summed parts; `None` if the
/// series does not converge within the term cap.
fn inv_log_sum(
    a: &ExactComplex,
    c: &ExactComplex,
    z: &ExactComplex,
    m: u32,
    p: usize,
    cc: &mut Consts,
) -> Option<(ExactComplex, Exponent)> {
    let re = |v: u32| ExactComplex::from_real(ExactNum::from_u32(v, p), p);
    let one = ExactComplex::one(p);
    let mz = z.neg_via_mul();
    let lnmz = mz.ln(p, rm(), cc);
    let pw = mz.pow(&a.neg_via_mul(), p, rm(), cc);
    let zinv = one.div(z, p, rm());
    let am = re(m).add(a, p, rm());
    let ca = ExactComplex::sub_ref(c, a, p, rm());
    let mut max_e = Exponent::MIN;
    let mut total = ExactComplex::zero(p);
    if m >= 1 {
        let pref1 = pw.mul(&rgamma_c(&am, p, cc), p, rm());
        let mut poch = ExactComplex::one(p);
        let mut rg = rgamma_c(&ca, p, cc);
        let mut zk = ExactComplex::one(p);
        let mut s1 = ExactComplex::zero(p);
        for k in 0..m {
            let t = poch
                .mul(&factorial_c(m - k - 1, p), p, rm())
                .div(&factorial_c(k, p), p, rm())
                .mul(&rg, p, rm())
                .mul(&zk, p, rm());
            max_e = max_e.max(cexp(&t.mul(&pref1, p, rm())));
            s1 = s1.add(&t, p, rm());
            poch = poch.mul(&re(k).add(a, p, rm()), p, rm());
            // 1/Γ(x − 1) = (x − 1)/Γ(x): exact zero once x passes a pole.
            rg = rg.mul(&ExactComplex::sub_ref(&ca, &re(k + 1), p, rm()), p, rm());
            zk = zk.mul(&zinv, p, rm());
        }
        total = pref1.mul(&s1, p, rm());
    }
    let rga = rgamma_c(a, p, cc);
    if is_c_zero(&rga) {
        return Some((total, max_e));
    }
    let pref2 = pw.mul(&rga, p, rm());
    // x = c − a − m − k, with r = 1/Γ(x) and q = ψ(x)/Γ(x).
    let mut x = ExactComplex::sub_ref(&ca, &re(m), p, rm());
    let mut pole_n: Option<u32> = None;
    let mut r = rgamma_c(&x, p, cc);
    let mut psi_x = ExactComplex::zero(p);
    let mut q;
    if let Some(n0) = as_i32_real_int(&x, p).filter(|&n| n <= 0) {
        pole_n = Some(n0.unsigned_abs());
        q = factorial_c(n0.unsigned_abs(), p);
        if n0 % 2 == 0 {
            q = neg_c(q);
        }
    } else {
        psi_x = x.digamma_at(p, cc);
        q = psi_x.mul(&r, p, rm());
    }
    let mut coef = powi_c(&zinv, m, p).div(&factorial_c(m, p), p, rm());
    let mut psi1m = re(m + 1).digamma_at(p, cc);
    let mut psi1 = one.digamma_at(p, cc);
    let mut psiam = am.digamma_at(p, cc);
    let mzinv = zinv.neg_via_mul();
    let big = max_abs_param(&[a, c], p) + u64::from(m);
    let mut sum = ExactComplex::zero(p);
    let budget = series_term_budget(&zinv, p, big);
    for k in 0..budget as u32 {
        let bracket = r
            .mul(
                &ExactComplex::add_ref(&lnmz, &psi1m, p, rm())
                    .add(&psi1, p, rm())
                    .sub(&psiam, p, rm()),
                p,
                rm(),
            )
            .sub(&q, p, rm());
        let term = coef.mul(&bracket, p, rm());
        sum = sum.add(&term, p, rm());
        max_e = max_e.max(cexp(&term.mul(&pref2, p, rm())));
        if u64::from(k) > big && term_small_rel(&term, &sum, p) {
            let v = total.add(&pref2.mul(&sum, p, rm()), p, rm());
            return Some((v, max_e));
        }
        let amk = re(k).add(&am, p, rm());
        coef = coef
            .mul(&amk, p, rm())
            .div(&re(k + 1).mul(&re(k + m + 1), p, rm()), p, rm())
            .mul(&mzinv, p, rm());
        psi1m = psi1m.add(&one.div(&re(m + k + 1), p, rm()), p, rm());
        psi1 = psi1.add(&one.div(&re(k + 1), p, rm()), p, rm());
        psiam = psiam.add(&one.div(&amk, p, rm()), p, rm());
        let xm1 = ExactComplex::sub_ref(&x, &one, p, rm());
        match pole_n {
            Some(n) => {
                // ψ(−n−1)/Γ(−n−1) = −(n+1)·ψ(−n)/Γ(−n).
                q = neg_c(q.mul(&re(n + 1), p, rm()));
                pole_n = Some(n + 1);
            }
            None if is_c_zero(&xm1) => {
                // x = 1 → 0: 1/Γ(0) = 0, ψ(0)/Γ(0) = −1.
                r = ExactComplex::zero(p);
                q = neg_c(ExactComplex::one(p));
                pole_n = Some(0);
            }
            None => {
                psi_x = psi_x.sub(&one.div(&xm1, p, rm()), p, rm());
                r = r.mul(&xm1, p, rm());
                q = psi_x.mul(&r, p, rm());
            }
        }
        x = xm1;
    }
    None
}

/// [`inv_log_sum`] times \(\Gamma(c)\), re-evaluated with the bits lost to cancellation.
fn transform_inv_log(
    a: &ExactComplex,
    c: &ExactComplex,
    z: &ExactComplex,
    m: u32,
    p: usize,
    cc: &mut Consts,
) -> ExactComplex {
    let mut pw = p + 32;
    for _ in 0..3 {
        let Some((v, max_e)) = inv_log_sum(a, c, z, m, pw, cc) else {
            return nan_pair(Error::InvalidArgument);
        };
        let ve = cexp(&v);
        let lost = if ve == Exponent::MIN { 0 } else { max_e.saturating_sub(ve).max(0) as usize };
        if ve == Exponent::MIN || lost + 16 <= pw - p {
            let gc = gamma_fixed(c, pw, cc);
            return gc.mul(&v, pw, rm());
        }
        pw = p + lost + 48;
    }
    nan_pair(Error::PrecisionRetryExhausted)
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
    // c = b: F(a, b; b; z) = (1 − z)^(−a) on the whole cut plane (and c = a likewise).
    if c_eq(c, b) {
        return one_minus(z, p).pow(&a.neg_via_mul(), p, rm(), cc);
    }
    if c_eq(c, a) {
        return one_minus(z, p).pow(&b.neg_via_mul(), p, rm(), cc);
    }
    // Integer c − a − b near z = 1: the direct series converges too slowly (or not at all for
    // |z| ≥ 1) and the generic 1 − z transform is singular, so use the logarithmic connection.
    if let Some(m) = int_cab(a, b, c, p) {
        let omz = ExactComplex::one(p).sub(z, p, rm());
        let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, rm());
        let near = matches!(omz.abs(p, rm()).compare(&half), Some(o) if o <= 0);
        if left > 0 && (near || (abs_lt_one(&omz, p) && !abs_lt_one(z, p))) {
            return transform_1mz_log(a, b, c, z, m, p, cc, left - 1);
        }
    }

    // Among the four maps with |w| < 1, use the smallest |w|. A fixed 10_000-term
    // cap used to refuse a computable series (P9); the term budget now scales with
    // p and 1 − |w|, and this picker prefers a cheaper map when one exists
    // (e.g. 1/z over Pfaff at z = −10 + 3i).
    let one = ExactComplex::one(p);
    let mut best_m: Option<ExactNum> = None;
    let mut best: u8 = 0; // 1 series, 2 pfaff, 3 1-z, 4 1/z, 5 1/z log
    let consider = |m: ExactNum, tag: u8, best_m: &mut Option<ExactNum>, best: &mut u8| {
        if !mag_lt(&m, &ExactNum::from_u8(1, 64)) {
            return;
        }
        if best_m.as_ref().is_none_or(|b| mag_lt(&m, b)) {
            *best_m = Some(m);
            *best = tag;
        }
    };
    consider(mag64(z), 1, &mut best_m, &mut best);
    let zm1 = ExactComplex::sub_ref(z, &one, p, rm());
    if !is_c_zero(&zm1) {
        consider(mag64(&z.div(&zm1, p, rm())), 2, &mut best_m, &mut best);
    }
    let cab = ExactComplex::sub_ref(c, a, p, rm()).sub(b, p, rm());
    if int_cab(a, b, c, p).is_none() && !is_nonpos_integer(&cab) {
        consider(mag64(&one_minus(z, p)), 3, &mut best_m, &mut best);
    }
    let inv = one.div(z, p, rm());
    if let Some((lo, _hi, m)) = int_ab(a, b, p) {
        consider(mag64(&inv), 5, &mut best_m, &mut best);
        if best == 5 {
            if left == 0 {
                return nan_pair(Error::InvalidArgument);
            }
            return transform_inv_log(lo, c, z, m, p, cc);
        }
    } else {
        consider(mag64(&inv), 4, &mut best_m, &mut best);
    }
    match best {
        1 => hypergeom_series(a, b, c, z, p),
        0 => nan_pair(Error::InvalidArgument),
        _ if left == 0 => nan_pair(Error::InvalidArgument),
        2 => pfaff(a, b, c, z, p, cc, left - 1),
        3 => transform_1mz_guarded(a, b, c, z, p, cc, left - 1),
        4 => transform_inv(a, b, c, z, p, cc, left - 1),
        _ => nan_pair(Error::InvalidArgument),
    }
}

/// `z = 1` with real `a, b, c`, `c − a − b ≤ 0` and a non-terminating series: the sum diverges,
/// `F → ±∞` as `z → 1⁻` (DLMF 15.4(ii)): the sign of `Γ(c)Γ(a+b−c)/(Γ(a)Γ(b))` for `c − a − b < 0`,
/// of `Γ(c)/(Γ(a)Γ(b))` for `c − a − b = 0` (`F ~ −Γ(c)/(Γ(a)Γ(b))·ln(1 − z)`). Imaginary part +0.
fn divergent_at_one(
    a: &ExactComplex,
    b: &ExactComplex,
    c: &ExactComplex,
    z: &ExactComplex,
    p: usize,
    cc: &mut Consts,
) -> Option<ExactComplex> {
    if !is_c_one(z, p)
        || !a.im().is_zero()
        || !b.im().is_zero()
        || !c.im().is_zero()
        || terminating_neg_int(a, p)
        || terminating_neg_int(b, p)
    {
        return None;
    }
    let cab = ExactNum::sub_ref(c.re(), a.re(), p, rm()).sub(b.re(), p, rm());
    if cab.is_positive() && !cab.is_zero() {
        return None;
    }
    let neg = |x: &ExactNum, cc: &mut Consts| x.gamma(64, rm(), cc).is_negative();
    let flip = neg(c.re(), cc) ^ neg(a.re(), cc) ^ neg(b.re(), cc);
    let inf = if flip { crate::INF_NEG } else { crate::INF_POS };
    Some(ExactComplex::new(inf, ExactNum::new(p)))
}

impl ExactComplex {
    /// Gaussian \({}_2F_1(a=\mathrm{self},b;c;z)\) in \(\mathbb{C}\).
    ///
    /// All-real arguments with \(z<1\) use the real [`ExactNum::hypergeom_2f1`].
    /// Otherwise: among the series (\(\lvert z\rvert<1\)), Pfaff, Euler / \(1-z\), and \(1/z\)
    /// maps, the argument of smallest modulus is used (NaN if a series does not converge
    /// within the \(p\)-scaled term bound). Kummer at \(z=1\)
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
    /// - Algorithm: series / Euler / Pfaff / Kummer / logarithmic \(1-z\) connection. The series term bound scales with `p` and `1 − |z|` (P9); transform recursion is limited by `HYPERGEOM_TRANSFORM_MAX = 8`.
    /// - Bound: Ziv on each part (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
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
            && matches!(z.re().compare(&ExactNum::from_u8(1, dest)), Some(o) if o < 0)
        {
            let r = self
                .re()
                .hypergeom_2f1(b.re(), c.re(), z.re(), dest, rm, cc);
            if !r.is_nan() {
                return ExactComplex::new(r, ExactNum::new(dest));
            }
        }
        if let Some(inf) = divergent_at_one(self, b, c, z, dest, cc) {
            return inf;
        }
        ziv_complex(dest, rm, |pw| {
            hypergeom_at(self, b, c, z, pw, cc, HYPERGEOM_TRANSFORM_MAX)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::ops::Neg;

    #[test]
    fn test_as_i32_real_int_beyond_10000() {
        let p = 64;
        let z = ExactComplex::new(ExactNum::from_u32(10_001, p), ExactNum::new(p));
        assert_eq!(as_i32_real_int(&z, p), Some(10_001));
        let zn = ExactComplex::new(ExactNum::from_i64(-10_001, p), ExactNum::new(p));
        assert_eq!(as_i32_real_int(&zn, p), Some(-10_001));
    }

    fn near(a: &ExactNum, b: &ExactNum, p: usize, slack: Exponent) -> bool {
        let d = ExactNum::sub_ref(a, b, p, RoundingMode::None).abs();
        d.is_zero() || d.exponent().unwrap_or(0) < -((p as Exponent) - slack)
    }

    fn tiny(x: &ExactNum, p: usize) -> bool {
        x.is_zero() || x.exponent().unwrap_or(0) < -((p as Exponent) / 4)
    }

    /// `trim_zero_words` drops only all-zero low words: 5 at 1024 bits keeps one word, 1/3 at
    /// 256 bits plus 7 at 1024 bits keeps the five words of the exact sum, a value with a non-zero
    /// lowest word keeps all of them, and zero, NaN and the values themselves are unchanged.
    #[test]
    fn trim_zero_words_keeps_the_value() {
        let r = RoundingMode::None;
        let five = trim_zero_words(ExactNum::from_u8(5, 1024));
        assert_eq!(five.mantissa_max_bit_len(), Some(64));
        assert_eq!(five.compare(&ExactNum::from_u8(5, 64)), Some(0));

        let third = ExactNum::from_u8(1, 256).div(&ExactNum::from_u8(3, 256), 256, r);
        let sum = ExactNum::from_u8(7, 1024).add(&third, 1024, r);
        let want = ExactNum::from_u8(7, 1024).add(&third, 1024, r);
        let t = trim_zero_words(sum);
        // 7 + 1/3 spans 2^2 down to 2^-257: 260 bits, five words.
        assert_eq!(t.mantissa_max_bit_len(), Some(320));
        assert_eq!(t.compare(&want), Some(0));
        assert_eq!(t.is_negative(), want.is_negative());

        let dense = ExactNum::from_u8(1, 1024).div(&ExactNum::from_u8(3, 1024), 1024, r);
        let d = trim_zero_words(dense);
        assert_eq!(d.mantissa_max_bit_len(), Some(1024));

        let neg = trim_zero_words(ExactNum::from_i8(-3, 512));
        assert_eq!(neg.mantissa_max_bit_len(), Some(64));
        assert_eq!(neg.compare(&ExactNum::from_i8(-3, 64)), Some(0));

        let z = trim_zero_words(ExactNum::new(512));
        assert!(z.is_zero());
        assert!(trim_zero_words(crate::NAN).is_nan());
    }

    /// `int_plus(n, a, p)` is `n + a` at `p` bits, both parts trimmed to their exact words.
    #[test]
    fn int_plus_is_the_trimmed_sum() {
        let p = 2048;
        let r = RoundingMode::None;
        let a = ExactComplex::new(
            ExactNum::from_u8(1, 256).div(&ExactNum::from_u8(3, 256), 256, r),
            ExactNum::from_u8(1, 128).div(&ExactNum::from_u8(7, 128), 128, r),
        );
        let s = int_plus(9, &a, p);
        let want = ExactComplex::from_real(ExactNum::from_u32(9, p), p).add(&a, p, r);
        assert_eq!(s.re().compare(want.re()), Some(0));
        assert_eq!(s.im().compare(want.im()), Some(0));
        assert_eq!(s.re().mantissa_max_bit_len(), Some(320));
        assert_eq!(s.im().mantissa_max_bit_len(), Some(128));
    }

    /// `abs_lead` is `|z|` from the leading 128 bits of each part: |3 + 4i| = 5 exactly, and on a
    /// 4096-bit 1/3 + i/7 its exponent is that of the full-width `abs`. `term_small_rel` against 1
    /// at p = 256: 2^-256 is not small, 2^-257 is.
    #[test]
    fn abs_lead_matches_abs() {
        let r = RoundingMode::None;
        let z = ExactComplex::new(ExactNum::from_u8(3, 64), ExactNum::from_u8(4, 64));
        assert_eq!(abs_lead(&z).compare(&ExactNum::from_u8(5, 64)), Some(0));
        let p = 4096;
        let w = ExactComplex::new(
            ExactNum::from_u8(1, p).div(&ExactNum::from_u8(3, p), p, r),
            ExactNum::from_u8(1, p).div(&ExactNum::from_u8(7, p), p, r),
        );
        let lead = abs_lead(&w);
        assert_eq!(lead.exponent(), w.abs(64, r).exponent());
        assert_eq!(lead.mantissa_max_bit_len(), Some(64));
        let tiny = ExactComplex::new(ExactNum::from_u8(1, p).ldexp(-300, p, r), ExactNum::new(p));
        let one = ExactComplex::one(p);
        assert!(term_small_rel(&tiny, &one, 256));
        assert!(!term_small_rel(&one, &tiny, 256));
        let edge = ExactComplex::new(ExactNum::from_u8(1, p).ldexp(-256, p, r), ExactNum::new(p));
        assert!(!term_small_rel(&edge, &one, 256));
        let below = ExactComplex::new(ExactNum::from_u8(1, p).ldexp(-257, p, r), ExactNum::new(p));
        assert!(term_small_rel(&below, &one, 256));
    }

    /// Exponent of `|z|` from exact arithmetic: `|z|^2 = re^2 + im^2` exactly, and
    /// `2^(e-1) <= |z| < 2^e` iff `2^(2e-2) <= |z|^2 < 2^(2e)`, so `e = ceil(e2 / 2)` for the
    /// exponent `e2` of `|z|^2`.
    fn exact_abs_exponent(z: &ExactComplex) -> Option<Exponent> {
        let (r, w) = (RoundingMode::None, 8192);
        let re2 = ExactNum::mul_ref(z.re(), z.re(), w, r);
        let im2 = ExactNum::mul_ref(z.im(), z.im(), w, r);
        let e2 = re2.add(&im2, w, r).exponent()?;
        Some(e2.div_euclid(2) + e2.rem_euclid(2))
    }

    /// `abs_lead` keeps the leading 128 bits of each part: its exponent equals the exact one of
    /// `|z|` just below, at and just above a power of two `2^m` (real and complex, both signs,
    /// m = 0, 1, 63, 64, 65, -64, 200, p = 256 and 1024), except when `|z|` lies within `2^-126`
    /// above `2^m` with a part that loses bits, where it may be one below (the documented
    /// deviation). Oracle: exact arithmetic (`exact_abs_exponent`).
    #[test]
    fn abs_lead_exponent_near_a_power_of_two() {
        let r = RoundingMode::None;
        for p in [256usize, 1024] {
            for m in [0i64, 1, 63, 64, 65, -64, 200] {
                let pow = |e: i64| ExactNum::from_u8(1, p).ldexp(e, p, r);
                for neg in [false, true] {
                    let s = |x: ExactNum| if neg { x.neg() } else { x };
                    let zero = || ExactNum::new(p);
                    // (value, may the lead exponent be one below the exact one?)
                    let cases = [
                        // just below 2^m: 2^m - 2^(m-200)
                        (
                            ExactComplex::new(s(pow(m).sub(&pow(m - 200), p, r)), zero()),
                            false,
                        ),
                        // at 2^m
                        (ExactComplex::new(s(pow(m)), zero()), false),
                        // just above 2^m, real: the cut never goes below 2^m
                        (
                            ExactComplex::new(s(pow(m).add(&pow(m - 200), p, r)), zero()),
                            false,
                        ),
                        // just above 2^m from exact parts: |z|^2 = 4^m + 4^(m-100)
                        (ExactComplex::new(s(pow(m)), s(pow(m - 100))), false),
                        // just above 2^m with a re part that loses bits in the cut:
                        // re = 2^m - 2^(m-150), im = 2^(m-70): |z|^2 = 4^m (1 - 2^-149 + 2^-140 + ...)
                        (
                            ExactComplex::new(s(pow(m).sub(&pow(m - 150), p, r)), s(pow(m - 70))),
                            true,
                        ),
                        // on the diagonal: re = im = 2^(m-1), |z| = 2^(m-1) sqrt 2
                        (ExactComplex::new(s(pow(m - 1)), s(pow(m - 1))), false),
                    ];
                    for (z, may_dip) in cases {
                        let lead = abs_lead(&z).exponent();
                        let exact = exact_abs_exponent(&z);
                        if may_dip {
                            assert!(
                                lead == exact
                                    || lead.zip(exact).map(|(l, e)| l + 1 == e) == Some(true),
                                "p={p} m={m} neg={neg}: lead {lead:?} exact {exact:?}"
                            );
                        } else {
                            assert_eq!(lead, exact, "p={p} m={m} neg={neg}: {z:?}");
                        }
                    }
                }
            }
        }
    }

    /// `term_small_rel(t, s, p)` is `|t| < 2^-p |s|` (exact for powers of two): with t = ±2^a and
    /// s = ±2^b it holds for b = a + p + 1 and fails for b = a + p and a + p - 1; t = 0 always holds.
    /// p = 64, 128, 256; a = 0, -64, 100. Oracle: exact (powers of two).
    #[test]
    fn term_small_rel_at_the_exponent_boundary() {
        let r = RoundingMode::None;
        let w = 512;
        for p in [64usize, 128, 256] {
            for a in [0i64, -64, 100] {
                let pt = |e: i64, neg: bool| {
                    let x = ExactNum::from_u8(1, w).ldexp(e, w, r);
                    ExactComplex::new(if neg { x.neg() } else { x }, ExactNum::new(w))
                };
                let pi = p as i64;
                for (nt, ns) in [(false, false), (true, false), (false, true), (true, true)] {
                    let t = pt(a, nt);
                    assert!(
                        term_small_rel(&t, &pt(a + pi + 1, ns), p),
                        "p={p} a={a} b=a+p+1"
                    );
                    assert!(!term_small_rel(&t, &pt(a + pi, ns), p), "p={p} a={a} b=a+p");
                    assert!(
                        !term_small_rel(&t, &pt(a + pi - 1, ns), p),
                        "p={p} a={a} b=a+p-1"
                    );
                }
                assert!(term_small_rel(&ExactComplex::zero(w), &pt(a, false), p));
            }
        }
    }

    #[test]
    fn complex_hypergeom_2f1_golds() -> Result<(), crate::Error> {
        let p = 256;
        let r = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let one = ExactComplex::one(p);
        let two = ExactComplex::from_real(ExactNum::from_u8(2, p), p);
        let mk_half = || ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, r);
        let half = mk_half();
        let zh = ExactComplex::from_real(mk_half(), p);
        let a_h = ExactComplex::from_real(mk_half(), p);
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
        assert!(!lhs.is_nan(), "Euler lhs defined");
        let cab = ExactComplex::sub_ref(&c, &a, p, r).sub(&b, p, r);
        let pref = ExactComplex::one(p).sub(&z, p, r).pow(&cab, p, r, &mut cc);
        let rhs = pref.mul(
            &a.neg()
                .add(&c, p, r)
                .hypergeom_2f1(&b.neg().add(&c, p, r), &c, &z, p, r, &mut cc),
            p,
            r,
        );
        assert!(near(lhs.re(), rhs.re(), p, 30), "Euler re");
        assert!(near(lhs.im(), rhs.im(), p, 30), "Euler im");

        let c0 = ExactComplex::zero(p);
        let bad = one.hypergeom_2f1(&one, &c0, &zh, p, r, &mut cc);
        assert!(bad.is_nan(), "c=0 → NaN");

        let two_r = ExactNum::from_u8(2, p);
        let mk_eps = || ExactNum::from_u8(1, p).ldexp(-40, p, RoundingMode::None);
        let above = ExactComplex::new(ExactNum::from_u8(2, p), mk_eps());
        let below = ExactComplex::new(two_r, mk_eps().neg());
        let fa = one.hypergeom_2f1(&one, &two, &above, p, r, &mut cc);
        let fb = one.hypergeom_2f1(&one, &two, &below, p, r, &mut cc);
        assert!(!fa.is_nan() && !fb.is_nan(), "cut defined");
        assert!(near(fa.re(), fb.re(), p, 20), "cut Re");
        assert_ne!(fa.im().compare(fb.im()), Some(0), "cut Im differs");
        assert!(
            near(&fa.into_parts().1.abs(), &fb.into_parts().1.abs(), p, 20),
            "cut |Im|"
        );
        Ok(())
    }
}
