//! Additive catalog specials on [`ExactNum`]: Scorer, Kelvin, Struve, Anger–Weber,
//! Clausen, Barnes \(G\), polygamma / Hurwitz, inverse Jacobi, \({}_pF_q\), Lambert \(W\),
//! and the polylogarithm. SoftFloat only (integer limbs; no IEEE arithmetic).
//!
//! Meijer \(G\), Fox \(H\), and Heun live in `heun_fox`.

use crate::common::util::bump_prec_retry;
use crate::defs::{MAX_PRECISION, WORD_BIT_SIZE};
use crate::Consts;
use crate::Error;
use crate::ExactComplex;
use crate::ExactNum;
use crate::RoundingMode;
use alloc::vec::Vec;
use core::ops::Neg;

#[path = "catalog_aws.rs"]
mod catalog_aws;

pub(crate) fn work_p(p: usize) -> usize {
    p.saturating_add(WORD_BIT_SIZE.saturating_mul(2))
}

/// Destination width: the requested precision rounded up to a whole word.
pub(crate) fn word_p(p: usize) -> usize {
    let w = WORD_BIT_SIZE.max(1);
    p.div_ceil(w).saturating_mul(w)
}

pub(crate) fn cat_nan() -> ExactNum {
    ExactNum::nan(Some(Error::InvalidArgument))
}

pub(crate) fn finite(x: &ExactNum) -> bool {
    !x.is_nan() && !x.is_inf()
}

pub(crate) fn finish(mut y: ExactNum, p: usize, rm: RoundingMode) -> ExactNum {
    let _ = y.set_precision(word_p(p), rm);
    y
}

/// Absolute residual below `2^{-(word width + 4)}`.
fn under_word(err: &ExactNum, p: usize) -> bool {
    let guard = (word_p(p) as i64).saturating_add(4);
    err.is_zero() || err.exponent().is_some_and(|e| e + guard < 0)
}

/// Residual below `2^{-(word width + 4)}` relative to `sum`.
fn under_word_rel(term: &ExactNum, sum: &ExactNum, p: usize) -> bool {
    let guard = (word_p(p) as i64).saturating_add(4);
    let rel = sum.exponent().unwrap_or(0);
    term.is_zero() || term.exponent().is_some_and(|e| e + guard < rel)
}

fn even_bernoulli_ext(n: usize, p: usize) -> Vec<ExactNum> {
    // The shared exact cache (O(K²) additions once) rather than a fresh O(K²) multiply triangle per call.
    let got = crate::ops::special::with_even_bernoulli(n, p, |v| {
        v.iter()
            .map(|b| b.clone().map(ExactNum::from))
            .collect::<Result<Vec<_>, _>>()
    });
    match got {
        Ok(v) => v,
        Err(e) => (0..n).map(|_| ExactNum::nan(Some(e))).collect(),
    }
}

fn factorial_n(n: usize, p: usize) -> ExactNum {
    let mut acc = ExactNum::from_u8(1, p);
    for k in 2..=n {
        acc = acc.mul(&ExactNum::from_u64(k as u64, p), p, RoundingMode::None);
    }
    acc
}

/// Positive integer value of `x` when it is an exact integer in `1..=u64::MAX`.
fn positive_int_u64(x: &ExactNum) -> Option<u64> {
    if !x.is_int() || !x.is_positive() {
        return None;
    }
    x.to_u64(RoundingMode::None).ok()
}

pub(crate) fn small_int(x: &ExactNum, p: usize, lo: i32, hi: i32) -> Option<i32> {
    if !x.is_int() {
        return None;
    }
    (lo..=hi).find(|&n| x.compare(&ExactNum::from_i32(n, p)) == Some(0))
}

/// `x` is 0, -1, -2, ... (any size).
fn nonpositive_int(x: &ExactNum) -> bool {
    x.is_int() && (x.is_zero() || x.is_negative())
}

/// Upper bound (64-bit, rounded up) on `|x|`.
fn abs_up(x: &ExactNum) -> ExactNum {
    let mut y = x.clone().abs();
    let _ = y.set_precision(64, RoundingMode::Up);
    y
}

/// DLMF 16.2.5: `{}_1F_0(a;;z) = (1-z)^{-a}` on the principal branch.
///
/// The real API refuses when `1-z < 0` and `a` is not an integer (the value is not real).
fn pfq_1f0(z: &ExactNum, a: &ExactNum, p: usize, cc: &mut Consts) -> ExactNum {
    let wrk = work_p(p).saturating_add(WORD_BIT_SIZE);
    let n = RoundingMode::None;
    let one = ExactNum::from_u8(1, wrk);
    let w = wrk
        .max(z.precision().unwrap_or(0))
        .max(a.precision().unwrap_or(0))
        .saturating_add(WORD_BIT_SIZE);
    let omz = ExactNum::sub_ref(&one, z, w, n);
    if !finite(&omz) || omz.is_zero() {
        return cat_nan();
    }
    if omz.is_negative() && !a.is_int() {
        return cat_nan();
    }
    let v = omz.pow(&a.clone().neg_via_mul(), wrk, n, cc);
    if finite(&v) {
        v
    } else {
        cat_nan()
    }
}

/// DLMF 16.18.1: F(a;b;z) = ΠΓ(b)/ΠΓ(a) · G(-z | 1-a; ; 0; 1-b).
/// For p = q+1, |z| > 1, the type-2 G residue series is in 1/z, inside the disk.
/// Negative G arguments (the cut z > 1) are refused unless every power is integer.
/// Real \({}_2F_1(a_0, a_1; c; z)\) at `z = ±1`, nonterminating, `c` not a pole.
///
/// - `z = −1`: Pfaff (DLMF 15.8.1), `2^{−a_0} ₂F₁(a_0, c−a_1; c; 1/2)`, a convergent
///   series, so the analytic value is computed whether or not the series at −1 converges.
/// - `z = 1`, `s = c − a_0 − a_1 > 0`: Gauss (DLMF 15.4.20),
///   `Γ(c)Γ(s)/(Γ(c−a_0)Γ(c−a_1))`, with a pole of a denominator Γ giving exactly 0.
/// - `z = 1`, `s ≤ 0`: the series diverges to ±∞ with the sign of
///   `Γ(c)/(Γ(a_0)Γ(a_1))`, the sign of its terms for large `k` (Jeff's ruling: mpmath's
///   ±Inf). An exact infinity, as for a pole.
fn hyp2f1_unit(z: &ExactNum, a: &[ExactNum], c: &ExactNum, p: usize, cc: &mut Consts) -> ExactNum {
    let wrk = work_p(p).saturating_add(WORD_BIT_SIZE);
    let n = RoundingMode::None;
    let wide = |x: &ExactNum| x.precision().unwrap_or(0).max(wrk) + WORD_BIT_SIZE;
    if z.is_negative() {
        let cma1 = ExactNum::sub_ref(c, &a[1], wide(c).max(wide(&a[1])) + 64, n);
        let half = ExactNum::from_u8(1, wrk).ldexp(-1, wrk, n);
        let f = pfq_series(&half, &[a[0].clone(), cma1], core::slice::from_ref(c), wrk);
        if !finite(&f) {
            return cat_nan();
        }
        let ln2 = ExactNum::from_u8(2, wrk).ln(wrk, n, cc);
        let scale = a[0].clone().neg_via_mul().mul(&ln2, wrk, n).exp(wrk, n, cc);
        return f.mul(&scale, wrk, n);
    }
    let w = wide(c).max(wide(&a[0])).max(wide(&a[1])) + 64;
    let cma0 = ExactNum::sub_ref(c, &a[0], w, n);
    let cma1 = ExactNum::sub_ref(c, &a[1], w, n);
    let s = ExactNum::sub_ref(&cma0, &a[1], w, n);
    if s.is_zero() || s.is_negative() {
        let g = c
            .gamma(wrk, n, cc)
            .div(&a[0].gamma(wrk, n, cc), wrk, n)
            .div(&a[1].gamma(wrk, n, cc), wrk, n);
        if !finite(&g) || g.is_zero() {
            return cat_nan();
        }
        return if g.is_negative() { crate::ext::INF_NEG } else { crate::ext::INF_POS };
    }
    if nonpositive_int(&cma0) || nonpositive_int(&cma1) {
        return ExactNum::new(wrk);
    }
    let v = c
        .gamma(wrk, n, cc)
        .mul(&s.gamma(wrk, n, cc), wrk, n)
        .div(&cma0.gamma(wrk, n, cc), wrk, n)
        .div(&cma1.gamma(wrk, n, cc), wrk, n);
    if finite(&v) {
        v
    } else {
        cat_nan()
    }
}

fn pfq_continued(
    z: &ExactNum,
    a: &[ExactNum],
    b: &[ExactNum],
    p: usize,
    cc: &mut Consts,
) -> ExactNum {
    let wrk = work_p(p).saturating_add(WORD_BIT_SIZE);
    let n = RoundingMode::None;
    let one = ExactNum::from_u8(1, wrk);
    let an: Vec<ExactNum> = a
        .iter()
        .map(|x| ExactNum::sub_ref(&one, x, wrk, n))
        .collect();
    let bq: Vec<ExactNum> = b
        .iter()
        .map(|x| ExactNum::sub_ref(&one, x, wrk, n))
        .collect();
    let mz = z.clone().neg();
    let g = mz.meijer_g(&an, &[], &[ExactNum::new(wrk)], &bq, wrk, n, cc);
    if !finite(&g) {
        return cat_nan();
    }
    let mut c = one;
    for x in b {
        c = c.mul(&x.gamma(wrk, n, cc), wrk, n);
    }
    for x in a {
        c = c.div(&x.gamma(wrk, n, cc), wrk, n);
    }
    let v = c.mul(&g, wrk, n);
    if finite(&v) {
        v
    } else {
        cat_nan()
    }
}

/// Direct \({}_pF_q\) series with a proven truncation bound and working precision
/// raised until rounding error is below the result (Ziv on the cancellation).
///
/// No term cap (ZF-39, ZF-73): after the index `K` exceeds every `|b_j|`, the term ratio
/// `|z| Π|K+a_j| / ((K+1) Π|K+b_j|)` is bounded for all later `k` by
/// `B(K) = |z| Π pair(K)`, where each `(K+|a_j|)/(K-|b_j|)` (and `(K+|a_0|)/(K+1)`,
/// taken as at least 1) does not increase. Once `B(K) < 1` the tail is at most
/// `|t_K| B/(1-B)`. Accumulated rounding is at most `(n+1) 2^{e_max - wrk + 2}`. When
/// that plus the tail is not `2^{-(p+8)}` relative to the sum, the sum is redone wider.
/// Callers make the series convergent: `p <= q`, or `p = q+1` with `|z| < 1`, or a
/// terminating parameter.
fn pfq_series(z: &ExactNum, a: &[ExactNum], b: &[ExactNum], p: usize) -> ExactNum {
    let lo = RoundingMode::Down;
    let up = RoundingMode::Up;
    let n = RoundingMode::None;
    let az = abs_up(z);
    let aa: Vec<ExactNum> = a.iter().map(abs_up).collect();
    let ab: Vec<ExactNum> = b.iter().map(abs_up).collect();
    let bmax = ab.iter().fold(ExactNum::new(64), |m, x| {
        if x.compare(&m) == Some(1) {
            x.clone()
        } else {
            m
        }
    });
    // For p = q+1 one a pairs with (K+1) (the largest |a|); the rest pair with b's.
    let pair_first = a.len() == b.len() + 1;
    let mut wrk = work_p(p);
    let one_lp = ExactNum::from_u8(1, 64);
    loop {
        let mut term = ExactNum::from_u8(1, wrk);
        let mut sum = term.clone();
        let mut emax: i128 = 1;
        let mut nterms: u64 = 1;
        let mut km1 = ExactNum::new(128); // k - 1, exact
        let tail: ExactNum;
        loop {
            let kf = km1.clone().add(&ExactNum::from_u8(1, 128), 128, n);
            let mut num = ExactNum::from_u8(1, wrk);
            for aj in a {
                num = num.mul(&ExactNum::add_ref(&km1, aj, wrk, n), wrk, n);
            }
            let mut den = kf.clone();
            let _ = den.set_precision(wrk, n);
            for bj in b {
                den = den.mul(&ExactNum::add_ref(&km1, bj, wrk, n), wrk, n);
            }
            term = term.mul(&num, wrk, n).div(&den, wrk, n).mul(z, wrk, n);
            km1 = kf;
            if term.is_zero() {
                tail = ExactNum::new(64);
                break;
            }
            if !finite(&term) {
                return cat_nan();
            }
            sum = sum.add(&term, wrk, n);
            nterms = nterms.saturating_add(1);
            if let Some(e) = term.exponent() {
                emax = emax.max(i128::from(e));
            }
            // Ratio bound for every step after this index K = km1 (= k).
            let kk = &km1;
            if kk.compare(&bmax.clone().add(&one_lp, 64, up)) != Some(1) {
                continue;
            }
            let mut bound = az.clone();
            if pair_first {
                // Largest |a| pairs with K+1.
                let mut idx = 0usize;
                for (i, x) in aa.iter().enumerate() {
                    if x.compare(&aa[idx]) == Some(1) {
                        idx = i;
                    }
                }
                let f = ExactNum::add_ref(kk, &aa[idx], 64, up).div(
                    &ExactNum::add_ref(kk, &one_lp, 64, lo),
                    64,
                    up,
                );
                bound = bound.mul(
                    if f.compare(&one_lp) == Some(1) { &f } else { &one_lp },
                    64,
                    up,
                );
                let rest = aa
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i != idx)
                    .map(|(_, x)| x);
                for (x, y) in rest.zip(ab.iter()) {
                    let f = ExactNum::add_ref(kk, x, 64, up).div(
                        &ExactNum::sub_ref(kk, y, 64, lo),
                        64,
                        up,
                    );
                    bound = bound.mul(&f, 64, up);
                }
            } else {
                // p <= q: pair a_j with b_j, then (K+1) and the remaining b's divide.
                for (j, y) in ab.iter().enumerate() {
                    let d = ExactNum::sub_ref(kk, y, 64, lo);
                    let f = match aa.get(j) {
                        Some(x) => ExactNum::add_ref(kk, x, 64, up).div(&d, 64, up),
                        None => one_lp.clone().div(&d, 64, up),
                    };
                    bound = bound.mul(&f, 64, up);
                }
                bound = bound.div(&ExactNum::add_ref(kk, &one_lp, 64, lo), 64, up);
            }
            if bound.compare(&one_lp) == Some(-1) {
                // tail <= |t| B/(1-B)
                let gap = ExactNum::sub_ref(&one_lp, &bound, 64, lo);
                let t = abs_up(&term).mul(&bound, 64, up).div(&gap, 64, up);
                let lim = sum.exponent().map(|e| i128::from(e) - p as i128 - 16);
                let small = t.is_zero()
                    || match (t.exponent(), lim) {
                        (Some(te), Some(l)) => i128::from(te) < l,
                        _ => false,
                    };
                if small {
                    tail = t;
                    break;
                }
            }
        }
        if sum.is_zero() {
            return sum;
        }
        let es = i128::from(sum.exponent().unwrap_or(0));
        // log2(nterms + 1) <= 64
        let ln = 64 - i128::from((nterms.saturating_add(1)).leading_zeros());
        let round_e = emax - wrk as i128 + 2 + ln;
        let tail_e = tail.exponent().map(i128::from).unwrap_or(i128::MIN / 2);
        let need = es - p as i128 - 8;
        if round_e < need && (tail.is_zero() || tail_e < need) {
            return sum;
        }
        let deficit = (round_e - need).max(0) as usize;
        let next = wrk.saturating_add(deficit).saturating_add(WORD_BIT_SIZE);
        wrk = next.max(wrk.saturating_add(wrk / 2));
    }
}

/// \(\zeta(2k)=(-1)^{k-1}B_{2k}(2\pi)^{2k}/(2\,(2k)!)\).
fn zeta_even(k: usize, p: usize, cc: &mut Consts) -> ExactNum {
    if k == 0 {
        return cat_nan();
    }
    let bs = even_bernoulli_ext(k, p);
    let b = &bs[k - 1];
    let pi = cc.pi(p, RoundingMode::None);
    let two = ExactNum::from_u8(2, p);
    let two_pi = ExactNum::from_u8(2, p).mul(&pi, p, RoundingMode::None);
    let mut pow = ExactNum::from_u8(1, p);
    for _ in 0..(2 * k) {
        pow = pow.mul(&two_pi, p, RoundingMode::None);
    }
    let den = two.mul(&factorial_n(2 * k, p), p, RoundingMode::None);
    let mut z = pow
        .mul(b, p, RoundingMode::None)
        .div(&den, p, RoundingMode::None);
    if k.is_multiple_of(2) {
        z = z.neg();
    }
    z
}

/// Front length \(N\) and Bernoulli count \(M\) for [`hurwitz_em`] (any \(a\ge0\)),
/// from the remainder ratio \(r_M\) evaluated with upward rounding, \(2\pi\ge 6\)
/// and \(N+a\ge N\). \(N\) doubles when the terms of \(r_M\) start growing first.
fn hurwitz_em_lengths(s: u64, p: usize) -> Option<(usize, usize)> {
    let up = RoundingMode::Up;
    let target = i64::try_from(p).ok()?.checked_add(8)?;
    let mut n_front = (16 + p / 4).max(24);
    while n_front <= (1 << 28) {
        let den = ExactNum::from_u64(6 * n_front as u64, 64);
        let mut prod = ExactNum::from_u64(s - 1, 64).mul(&ExactNum::from_u8(4, 64), 64, up);
        let mut m = 0usize;
        loop {
            let i = 2 * m as u64;
            let f0 = ExactNum::from_u64(s.checked_add(i)?, 64);
            let f1 = ExactNum::from_u64(s.checked_add(i + 1)?, 64);
            if !matches!(f1.compare(&den), Some(c) if c < 0) {
                break;
            }
            prod = prod
                .mul(&f0, 64, up)
                .div(&den, 64, up)
                .mul(&f1, 64, up)
                .div(&den, 64, up);
            m += 1;
            let r = prod.clone().div(
                &ExactNum::from_u64(s.checked_add(2 * m as u64 - 1)?, 64),
                64,
                up,
            );
            if r.is_zero() || r.exponent().is_some_and(|e| e < -target) {
                return Some((n_front, m.max(1)));
            }
        }
        n_front = n_front.saturating_mul(2);
    }
    None
}

/// Euler–Maclaurin for \(\zeta(s,a)=\sum_{k=0}^\infty(k+a)^{-s}\), integer \(s\ge 2\), \(a>0\).
fn hurwitz_em(s: u64, a: &ExactNum, p: usize, cc: &mut Consts) -> ExactNum {
    if s < 2 || !finite(a) || a.is_negative() || a.is_zero() {
        return cat_nan();
    }
    // The former `m_bern = clamp(4, 48)` stopped the Bernoulli tail near 2^-500, so
    // ψ^{(n)} and ζ(s, a) were silently wrong past about 500 bits. The tail length now
    // comes from the remainder bound (Johansson 2015, Thm 1, real s):
    // |R_M| ≤ 4 (s)_{2M} / ((2π)^{2M} (s+2M−1) (N+a)^{s+2M−1}). Against the lower bound
    // ζ(s,a) ≥ (N+a)^{1−s}/(s−1) that is
    // r_M = 4 (s−1)/(s+2M−1) · Π_{i<2M} (s+i)/(2π(N+a)), required below 2^{−p−8}.
    let Some((n_front, m_bern)) = hurwitz_em_lengths(s, p) else {
        return cat_nan();
    };
    let Some(s_i) = i64::try_from(s).ok() else {
        return cat_nan();
    };
    let mut sum = ExactNum::new(p);
    for k in 0..n_front {
        let t = ExactNum::from_u64(k as u64, p).add(a, p, RoundingMode::None);
        let term = ExactNum::from_u8(1, p).div(
            &t.powsi(s_i, p, RoundingMode::None),
            p,
            RoundingMode::None,
        );
        sum = sum.add(&term, p, RoundingMode::None);
    }
    let n_a = ExactNum::from_u64(n_front as u64, p).add(a, p, RoundingMode::None);
    let sm1 = ExactNum::from_u64(s - 1, p);
    let tail_int = ExactNum::from_u8(1, p).div(
        &n_a.powsi(s_i - 1, p, RoundingMode::None)
            .mul(&sm1, p, RoundingMode::None),
        p,
        RoundingMode::None,
    );
    let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, RoundingMode::None);
    let tail_half = half.div(
        &n_a.powsi(s_i, p, RoundingMode::None),
        p,
        RoundingMode::None,
    );
    sum = sum
        .add(&tail_int, p, RoundingMode::None)
        .add(&tail_half, p, RoundingMode::None);
    let bs = even_bernoulli_ext(m_bern, p);
    // Term j: B_2j (s)_{2j-1} / ((2j)! (N+a)^{s+2j-1}), built incrementally (the former
    // per-j rising product, factorial and power were O(M^2) multiplies).
    let nr = RoundingMode::None;
    let inv_na = ExactNum::from_u8(1, p).div(&n_a, p, nr);
    let inv_na2 = ExactNum::mul_ref(&inv_na, &inv_na, p, nr);
    let mut coef = ExactNum::from_u64(s, p).div(&ExactNum::from_u8(2, p), p, nr);
    let Some(sp1) = s_i.checked_add(1) else {
        return cat_nan();
    };
    let mut inv_pow = inv_na.powsi(sp1, p, nr);
    for j in 1..=m_bern {
        if j > 1 {
            let j2 = 2 * j as u64;
            // (s+2j-3)(s+2j-2) as ExactNum, so large s cannot overflow u64.
            let t0 = ExactNum::from_u64(s, p).add(&ExactNum::from_u64(j2 - 3, p), p, nr);
            let t1 = ExactNum::from_u64(s, p).add(&ExactNum::from_u64(j2 - 2, p), p, nr);
            coef = coef.mul(&t0.mul(&t1, p, nr), p, nr).div(
                &ExactNum::from_u64((j2 - 1) * j2, p),
                p,
                nr,
            );
            inv_pow = inv_pow.mul(&inv_na2, p, nr);
        }
        let term = ExactNum::mul_ref(&bs[j - 1], &coef, p, nr).mul(&inv_pow, p, nr);
        sum = sum.add(&term, p, nr);
    }
    let _ = cc;
    sum
}

/// \(m\)-th derivative of \(\ln x/x^2\) at `x` (Leibniz).
#[cfg(test)]
fn zeta_prime_2(p: usize, cc: &mut Consts) -> ExactNum {
    // Euler–Maclaurin for Σ ln(n)/n² = −ζ'(2). The former 36-term / N≤256 caps left an
    // error around 2^-440, corrupting ln A and Barnes G at high precision.
    // d^m(ln x/x²) = (−1)^m (m+1)! x^(−m−2) [ln x − H_{m+1} + 1].
    let nr = RoundingMode::None;
    let n_front = (p + 64).div_ceil(2).max(48);
    let mut sum = ExactNum::new(p);
    for n in 2..=n_front {
        let nf = ExactNum::from_u64(n as u64, p);
        let term = nf.ln(p, nr, cc).div(&nf.sqr(p, nr), p, nr);
        sum = sum.add(&term, p, nr);
    }
    let n = ExactNum::from_u64((n_front + 1) as u64, p);
    let ln = n.ln(p, nr, cc);
    let one = ExactNum::from_u8(1, p);
    let invn = ExactNum::from_u8(1, p).div(&n, p, nr);
    let invn2 = invn.sqr(p, nr);
    let integral = ExactNum::add_ref(&one, &ln, p, nr).mul(&invn, p, nr);
    let mid = ExactNum::mul_ref(&ln, &invn2, p, nr).ldexp(-1, p, nr);
    sum = sum.add(&integral, p, nr).add(&mid, p, nr);
    let mut bs: Vec<ExactNum> = Vec::new();
    let mut invpower = invn;
    let mut h = one.clone();
    let mut j = 1usize;
    let mut last: Option<crate::Exponent> = None;
    loop {
        if j > bs.len() {
            bs = even_bernoulli_ext((2 * j).max(32), p);
        }
        // H_{2j}: previous h is H_{2j−2} (H_1 on the first step).
        if j > 1 {
            h = h.add(
                &ExactNum::from_u8(1, p).div(&ExactNum::from_u64((2 * j - 1) as u64, p), p, nr),
                p,
                nr,
            );
        }
        h = h.add(
            &ExactNum::from_u8(1, p).div(&ExactNum::from_u64((2 * j) as u64, p), p, nr),
            p,
            nr,
        );
        invpower = invpower.mul(&invn2, p, nr);
        let term = ExactNum::sub_ref(&ln, &h, p, nr)
            .add(&one, p, nr)
            .mul(&bs[j - 1], p, nr)
            .mul(&invpower, p, nr);
        if term.is_nan() {
            return term;
        }
        let et = term.exponent().unwrap_or(0);
        sum = sum.add(&term, p, nr);
        if term.is_zero() || et < -(p as crate::Exponent) - 16 {
            return sum.neg();
        }
        if matches!(last, Some(l) if et > l) {
            return cat_nan();
        }
        last = Some(et);
        j += 1;
    }
}

/// \(\ln A = (\gamma+\ln(2\pi))/12 - \zeta'(2)/(2\pi^2)\).
/// A constant at `p` bits from a shared cache holding the widest value computed so far.
fn cached_const(
    #[cfg(feature = "std")] slot: &std::sync::RwLock<Option<ExactNum>>,
    p: usize,
    cc: &mut Consts,
    f: fn(usize, &mut Consts) -> ExactNum,
) -> ExactNum {
    let round = |mut r: ExactNum| {
        let _ = r.set_precision(p, RoundingMode::ToEven);
        r
    };
    #[cfg(feature = "std")]
    {
        if let Ok(g) = slot.read() {
            if let Some(v) = g
                .as_ref()
                .filter(|v| v.precision().unwrap_or(0) >= p + WORD_BIT_SIZE)
            {
                return round(v.clone());
            }
        }
        let v = f(p + WORD_BIT_SIZE, cc);
        if !v.is_nan() {
            if let Ok(mut g) = slot.write() {
                *g = Some(v.clone());
            }
        }
        round(v)
    }
    #[cfg(not(feature = "std"))]
    round(f(p + WORD_BIT_SIZE, cc))
}

#[cfg(feature = "std")]
#[cfg(test)]
static LN_GLAISHER: std::sync::RwLock<Option<ExactNum>> = std::sync::RwLock::new(None);
#[cfg(feature = "std")]
static ZETA3: std::sync::RwLock<Option<ExactNum>> = std::sync::RwLock::new(None);

/// ln A at `p` bits (ζ'(2) costs about p/2 logarithms; Barnes G would pay it on every call).
#[cfg(test)]
fn ln_glaisher(p: usize, cc: &mut Consts) -> ExactNum {
    cached_const(
        #[cfg(feature = "std")]
        &LN_GLAISHER,
        p,
        cc,
        ln_glaisher_uncached,
    )
}

/// ζ(3) at `p` bits (Clausen Cl3 needs it on every call).
fn zeta3(p: usize, cc: &mut Consts) -> ExactNum {
    cached_const(
        #[cfg(feature = "std")]
        &ZETA3,
        p,
        cc,
        |q, c| ExactNum::from_u8(3, q).riemann_zeta(q, RoundingMode::None, c),
    )
}

#[cfg(test)]
fn ln_glaisher_uncached(p: usize, cc: &mut Consts) -> ExactNum {
    let g = cc.euler_gamma(p, RoundingMode::None);
    let two = ExactNum::from_u8(2, p);
    let pi = cc.pi(p, RoundingMode::None);
    let twopi = ExactNum::from_u8(2, p).mul(&pi, p, RoundingMode::None);
    let twelve = ExactNum::from_u32(12, p);
    let a = g
        .add(&twopi.ln(p, RoundingMode::None, cc), p, RoundingMode::None)
        .div(&twelve, p, RoundingMode::None);
    let zp = zeta_prime_2(p, cc);
    let pi2 = pi.sqr(p, RoundingMode::None);
    let b = zp.div(&two.mul(&pi2, p, RoundingMode::None), p, RoundingMode::None);
    a.sub(&b, p, RoundingMode::None)
}

fn reduce_clausen_angle(theta: &ExactNum, p: usize, cc: &mut Consts) -> ExactNum {
    // Use the same exponent-aware π reduction as sin/cos; a p-bit 2π loses all bits of the
    // remainder when the argument's exponent is much larger than p.
    let mut t = theta.rem_pi(p, RoundingMode::None, cc);
    let pi = cc.pi(p, RoundingMode::None);
    let two_pi = ExactNum::from_u8(2, p).mul(&pi, p, RoundingMode::None);
    if matches!(t.compare(&pi), Some(c) if c > 0) {
        t = t.sub(&two_pi, p, RoundingMode::None);
    } else if matches!(t.compare(&pi.neg()), Some(c) if c < 0) {
        t = t.add(&two_pi, p, RoundingMode::None);
    }
    t
}

/// `θ ≡ x (mod 2π)` in `[−π, π]` and `u = π − |θ|`, each to `wrk + 64` bits relative to its own size: the
/// reduction is widened until neither a near-multiple of 2π (θ ≈ 0) nor of π (u ≈ 0) cancels the bits.
/// `None` only when the reduction itself fails (NaN).
fn clausen_reduce(x: &ExactNum, wrk: usize, cc: &mut Consts) -> Option<(ExactNum, ExactNum)> {
    let reduced = x.exponent().is_some_and(|e| e > 1);
    let mut pr = word_p(wrk.saturating_add(128));
    loop {
        let th = if reduced { reduce_clausen_angle(x, pr, cc) } else { x.clone() };
        if th.is_nan() {
            return None;
        }
        let pi = cc.pi(pr, RoundingMode::None);
        let u = pi.sub(&ExactNum::clone(&th).abs(), pr, RoundingMode::None);
        let small = |v: &ExactNum| -> Option<usize> {
            if v.is_zero() {
                None
            } else {
                Some(v.exponent().unwrap_or(0).min(0).unsigned_abs() as usize)
            }
        };
        let lost = match (if reduced { small(&th) } else { Some(0) }, small(&u)) {
            (Some(a), Some(b)) => Some(a.max(b)),
            _ => None,
        };
        match lost {
            Some(l) if pr >= wrk.saturating_add(128).saturating_add(l) => {
                let mut th = th;
                let mut u = u;
                let _ = th.set_precision(word_p(wrk.saturating_add(64)), RoundingMode::None);
                let _ = u.set_precision(word_p(wrk.saturating_add(64)), RoundingMode::None);
                return Some((th, u));
            }
            Some(l) => pr = word_p(wrk.saturating_add(128).saturating_add(l)),
            None => pr = word_p(pr.saturating_mul(2)),
        }
    }
}

/// Largest integer `n ≤ cap` with `n ≤ x`. `None` when `x` is not a positive finite value.
fn floor_small(x: &ExactNum, p: usize, cap: u32) -> Option<u32> {
    if !finite(x) || !x.is_positive() {
        return None;
    }
    let mut n = 0u32;
    while n < cap {
        if matches!(x.compare(&ExactNum::from_u32(n + 1, p)), Some(c) if c < 0) {
            return Some(n);
        }
        n += 1;
    }
    Some(cap)
}

/// Terms of the Bernoulli series so the tail is below `2^{-bits}` on \((-\pi,\pi)\).
fn clausen_nterms(theta: &ExactNum, bits: usize, cc: &mut Consts) -> usize {
    let p = 64;
    let negated;
    let abs: &ExactNum = if theta.is_negative() {
        negated = theta.neg_via_mul();
        &negated
    } else {
        theta
    };
    if abs.is_zero() {
        return 8;
    }
    let two_pi = ExactNum::from_u8(2, p).mul(&cc.pi(p, RoundingMode::None), p, RoundingMode::None);
    let ratio = two_pi.div(abs, p, RoundingMode::None);
    let ln2 = ExactNum::from_u8(2, p).ln(p, RoundingMode::None, cc);
    let log2_ratio = ratio
        .ln(p, RoundingMode::None, cc)
        .div(&ln2, p, RoundingMode::None);
    let per_term = ExactNum::from_u8(2, p).mul(&log2_ratio, p, RoundingMode::None);
    let per = floor_small(&per_term, p, 32).unwrap_or(2).clamp(2, 32) as usize;
    // No upper bound: |θ| ≤ 2π/3 here (callers reflect about ±π), so a term gains >= 3 bits.
    bits.div_ceil(per).saturating_add(6).max(8)
}

/// \(\mathrm{Cl}_2(\theta)=\theta-\theta\ln|\theta|-\sum_{k\ge1}(-1)^k B_{2k}\theta^{2k+1}/(2k\cdot(2k+1)!)\)
/// on \((-\pi,\pi)\setminus\{0\}\).
fn clausen_cl2_series(theta: &ExactNum, p: usize, bits: usize, cc: &mut Consts) -> ExactNum {
    let nterms = clausen_nterms(theta, bits, cc);
    let bs = even_bernoulli_ext(nterms, p);
    let th2 = theta.sqr(p, RoundingMode::None);
    let mut u = ExactNum::mul_ref(theta, &th2, p, RoundingMode::None).div(
        &ExactNum::from_u8(6, p),
        p,
        RoundingMode::None,
    );
    let ln_abs = if theta.is_negative() {
        theta.neg_via_mul().ln(p, RoundingMode::None, cc)
    } else {
        theta.ln(p, RoundingMode::None, cc)
    };
    let mut sum = ln_abs
        .mul(theta, p, RoundingMode::None)
        .neg()
        .add(theta, p, RoundingMode::None);
    for (k, b) in (1..=nterms).zip(bs) {
        if b.is_nan() {
            return cat_nan();
        }
        let raw = b.mul(&u, p, RoundingMode::None).div(
            &ExactNum::from_u32((2 * k) as u32, p),
            p,
            RoundingMode::None,
        );
        let signed = if k % 2 == 1 { raw.neg() } else { raw };
        sum = sum.sub(&signed, p, RoundingMode::None);
        if under_word_rel(&signed, &sum, bits) {
            break;
        }
        let step = ExactNum::from_u32((2 * k + 2) as u32, p).mul(
            &ExactNum::from_u32((2 * k + 3) as u32, p),
            p,
            RoundingMode::None,
        );
        u = u
            .mul(&th2, p, RoundingMode::None)
            .div(&step, p, RoundingMode::None);
    }
    sum
}

/// \(\mathrm{Cl}_3(|\theta|)=\zeta(3)-\int_0^{|\theta|}\mathrm{Cl}_2\), integrated term by term.
fn clausen_cl3_series(theta_abs: &ExactNum, p: usize, bits: usize, cc: &mut Consts) -> ExactNum {
    let nterms = clausen_nterms(theta_abs, bits, cc);
    let bs = even_bernoulli_ext(nterms, p);
    let th2 = theta_abs.sqr(p, RoundingMode::None);
    let mut u = ExactNum::mul_ref(theta_abs, &th2, p, RoundingMode::None).div(
        &ExactNum::from_u8(6, p),
        p,
        RoundingMode::None,
    );
    let ln_abs = theta_abs.ln(p, RoundingMode::None, cc);
    let mut integ = ExactNum::from_u8(3, p)
        .mul(&th2, p, RoundingMode::None)
        .div(&ExactNum::from_u8(4, p), p, RoundingMode::None)
        .sub(
            &ExactNum::div_ref(&th2, &ExactNum::from_u8(2, p), p, RoundingMode::None).mul(
                &ln_abs,
                p,
                RoundingMode::None,
            ),
            p,
            RoundingMode::None,
        );
    for (k, b) in (1..=nterms).zip(bs) {
        if b.is_nan() {
            return cat_nan();
        }
        let raw = b.mul(&u, p, RoundingMode::None).div(
            &ExactNum::from_u32((2 * k) as u32, p),
            p,
            RoundingMode::None,
        );
        let signed = if k % 2 == 1 { raw.neg() } else { raw };
        let piece = signed.mul(theta_abs, p, RoundingMode::None).div(
            &ExactNum::from_u32((2 * k + 2) as u32, p),
            p,
            RoundingMode::None,
        );
        integ = integ.sub(&piece, p, RoundingMode::None);
        if under_word(&piece, bits) {
            break;
        }
        let step = ExactNum::from_u32((2 * k + 2) as u32, p).mul(
            &ExactNum::from_u32((2 * k + 3) as u32, p),
            p,
            RoundingMode::None,
        );
        u = u
            .mul(&th2, p, RoundingMode::None)
            .div(&step, p, RoundingMode::None);
    }
    let z3 = zeta3(p, cc);
    z3.sub(&integ, p, RoundingMode::None)
}

/// Ziv loop for a catalog kernel. `f(wp)` returns a value and its correct bits `s`
/// (`|error| < 2^(e(v) - s - 1)`); `f` is retried with more precision until the value
/// rounds to `p` bits with certainty. Past `8p + 16384` bits the answer is NaN, never a guess.
fn ziv_round<F>(p: usize, rm: RoundingMode, mut f: F) -> ExactNum
where
    F: FnMut(usize) -> Option<(ExactNum, i64)>,
{
    if let Err(e) = crate::common::util::reserve_probe(crate::common::util::round_p(p)) {
        return ExactNum::nan(Some(e));
    }
    let w = WORD_BIT_SIZE as i64;
    let limit = p.saturating_mul(8).saturating_add(16_384);
    let mut wp = work_p(p);
    loop {
        let Some((mut v, good)) = f(wp) else {
            return cat_nan();
        };
        if !finite(&v) {
            return cat_nan();
        }
        let s = good.div_euclid(w).saturating_mul(w);
        if !v.is_zero() && s > p as i64 && v.try_set_precision(p, rm, s as usize) {
            return v;
        }
        let short = (p as i64).saturating_add(2 * w).saturating_sub(good).max(w);
        let next = (wp as i64).saturating_add(short.max(wp as i64 / 2));
        wp = word_p(usize::try_from(next).unwrap_or(usize::MAX));
        if wp > limit {
            return cat_nan();
        }
    }
}

/// Bits in `n` (`0` for `0`).
fn bitlen(n: usize) -> i64 {
    (usize::BITS - n.leading_zeros()) as i64
}

/// Taylor series of Scorer \(Gi\) (`hi = false`) or \(Hi\) at `x`: \(w''=xw\mp1/\pi\) gives
/// \(a_n=a_{n-3}/(n(n-1))\) with \(a_0, a_1 = (1\text{ or }2)\,\mathrm{Bi}(0)/3,\ (1\text{ or }2)\,\mathrm{Bi}'(0)/3\)
/// and \(a_2=\mp1/(2\pi)\) (DLMF 9.12.1, 9.12.4-7). Once \(q=|x|^3/(n(n+1))\le1/2\) each residue class
/// mod 3 is dominated by a geometric series of ratio \(q\), so the tail is below
/// \(3\cdot2^{e_{\max}}\) for the last three terms. Summing stops when that is under
/// \(2^{\mathrm{peak}-wp}\); every term carries at most \(2n+6\) ulps and each addition one ulp
/// of at most \(n\) times the peak, so the error is below \(2^{\mathrm{peak}-wp+2\,\mathrm{bitlen}(n)+3}\).
fn scorer_taylor(x: &ExactNum, hi: bool, wp: usize, cc: &mut Consts) -> Option<(ExactNum, i64)> {
    let nr = RoundingMode::None;
    let three = ExactNum::from_u8(3, wp);
    let scale = ExactNum::from_u8(if hi { 2 } else { 1 }, wp);
    let a0 = ExactNum::new(wp)
        .bi(wp, nr, cc)
        .mul(&scale, wp, nr)
        .div(&three, wp, nr);
    // |x| < 2^-(wp+8): the rest is below |a1 x| + |a2| x^2 < 2^(ex-1), under a0's own error.
    // (Forming x^2 or a1 x would also drag the arithmetic through the subnormal range.)
    if x.is_zero() || x.exponent().is_some_and(|e| e < -(wp as i64) - 8) {
        return Some((a0, wp as i64 - 8));
    }
    let a1 = ExactNum::new(wp)
        .bi_prime(wp, nr, cc)
        .mul(&scale, wp, nr)
        .div(&three, wp, nr);
    let pi = cc.pi(wp, nr);
    let mut a2 = ExactNum::from_u8(1, wp).div(&ExactNum::from_u8(2, wp).mul(&pi, wp, nr), wp, nr);
    if !hi {
        a2 = a2.neg();
    }
    let ex = x.exponent()?;
    let t1 = ExactNum::mul_ref(&a1, x, wp, nr);
    let mut xn = ExactNum::mul_ref(x, x, wp, nr);
    let t2 = ExactNum::mul_ref(&a2, &xn, wp, nr);
    let exp_of = |t: &ExactNum| t.exponent().unwrap_or(i64::MIN / 4);
    let mut last = [exp_of(&a0), exp_of(&t1), exp_of(&t2)];
    let mut peak = last[0].max(last[1]).max(last[2]);
    let mut sum = ExactNum::add_ref(&a0, &t1, wp, nr).add(&t2, wp, nr);
    let mut a = [a0, a1, a2];
    let mut n = 2usize;
    loop {
        // q = |x|^3 / (n(n+1)) < 2^(3 ex - bitlen(n(n+1)) + 1).
        let nn = n.checked_mul(n + 1)?;
        let qexp = ex
            .saturating_mul(3)
            .saturating_sub(bitlen(nn))
            .saturating_add(1);
        let emax = last[0].max(last[1]).max(last[2]);
        if qexp <= -1 && emax + 2 < peak - wp as i64 {
            break;
        }
        n += 1;
        let den = ExactNum::from_u32(u32::try_from(n).ok()?, wp).mul(
            &ExactNum::from_u32(u32::try_from(n - 1).ok()?, wp),
            wp,
            nr,
        );
        let an = ExactNum::div_ref(&a[n % 3], &den, wp, nr);
        xn = xn.mul(x, wp, nr);
        let t = ExactNum::mul_ref(&an, &xn, wp, nr);
        a[n % 3] = an;
        last[n % 3] = exp_of(&t);
        peak = peak.max(last[n % 3]);
        sum = sum.add(&t, wp, nr);
    }
    if sum.is_zero() {
        return Some((sum, i64::MIN / 4));
    }
    let good = exp_of(&sum) - peak + wp as i64 - 2 * bitlen(n) - 6;
    Some((sum, good))
}

/// Watson's lemma for \(y>0\): \(\pi\,\mathrm{Hi}(-y)=\int_0^\infty e^{-yt-t^3/3}dt\) (`alt`) and
/// \(\pi\,\mathrm{Gi}(y)=\operatorname{Im}\int_0^\infty e^{i(yt+t^3/3)}dt\) both equal
/// \(\sum_{k<K}(\pm1)^k c_k+R_K\), \(c_k=(3k)!/(3^kk!\,y^{3k+1})\) (DLMF 9.12.20, 9.12.27, 9.12.29).
/// For \(\mathrm{Hi}\) the alternating Taylor remainder of \(e^{-u}\), \(u\ge0\), gives \(|R_K|\le c_K\).
/// For \(\mathrm{Gi}\) the contour is turned to \(\arg t=\pi/3\), where \(e^{-is^3/3}\) has
/// \(|R|\le(s^3/3)^K/K!\) and \(e^{-\lambda s}\) has \(\operatorname{Re}\lambda=y\sqrt3/2\), so
/// \(|R_K|\le c_K(2/\sqrt3)^{3K+1}<c_K2^{0.21(3K+1)}\). Summing stops at the first \(K\) with that
/// bound under \(2^{e(\mathrm{sum})-wp-2}\); `None` when the bound stops shrinking first.
fn scorer_watson(y: &ExactNum, alt: bool, wp: usize, cc: &mut Consts) -> Option<(ExactNum, i64)> {
    let nr = RoundingMode::None;
    let ey = y.exponent()?;
    if ey < 2 {
        // y < 2: c_1/c_0 = 2/y^3 > 1/4, the bound cannot reach working precision.
        return None;
    }
    let inv_pi = ExactNum::from_u8(1, wp).div(&cc.pi(wp, nr), wp, nr);
    if ey > wp as i64 + 8 {
        // K = 1: |R_1| <= c_1 (2/sqrt 3)^4 < 3.6 c_0 / y^3 < 2^(2 - 3(ey - 1)) c_0, far below 2^-wp.
        // 1/pi first: pi * y would overflow for y near the largest value.
        return Some((inv_pi.div(y, wp, nr), wp as i64 - 4));
    }
    let y3 = ExactNum::mul_ref(y, y, wp, nr).mul(y, wp, nr);
    let mut c = ExactNum::from_u8(1, wp).div(y, wp, nr);
    let e0 = c.exponent()?;
    let mut sum = c.clone();
    let mut k = 0usize;
    loop {
        let f = ExactNum::from_word(((3 * k + 1) as u64) * ((3 * k + 2) as u64), wp);
        let cn = c.mul(&f, wp, nr).div(&y3, wp, nr);
        let kk = 3 * (k + 1) + 1;
        let gb = if alt { 0 } else { (kk * 21).div_ceil(100) as i64 + 1 };
        let bexp = cn.exponent().unwrap_or(i64::MIN / 4).saturating_add(gb);
        if bexp < sum.exponent().unwrap_or(0) - wp as i64 - 2 {
            k += 1;
            break;
        }
        // The bound's next ratio, (3k+4)(3k+5)/y^3 (times 1.54 for Gi), must stay below 1.
        let m = ((3 * k + 4) as u64) * ((3 * k + 5) as u64);
        let m = if alt { m } else { (m * 154).div_ceil(100) };
        if !matches!(ExactNum::from_word(m, 64).compare(&y3), Some(o) if o < 0) {
            return None;
        }
        sum = if alt && k.is_multiple_of(2) {
            sum.sub(&cn, wp, nr)
        } else {
            sum.add(&cn, wp, nr)
        };
        c = cn;
        k += 1;
    }
    let es = sum.exponent()?;
    let good = es - e0 + wp as i64 - 2 * bitlen(k) - 8;
    Some((sum.mul(&inv_pi, wp, nr), good))
}

/// `a - b` with absolute errors below `2^ea` and `2^eb`: the difference and its correct bits.
fn diff_good(a: &ExactNum, ea: i64, b: &ExactNum, eb: i64, wp: usize) -> (ExactNum, i64) {
    let v = ExactNum::sub_ref(a, b, wp, RoundingMode::None);
    let err = ea.max(eb) + 1;
    let e = v.exponent().unwrap_or(i64::MIN / 4);
    // The subtraction itself adds one ulp at wp.
    (v, e - err.max(e - wp as i64) - 2)
}

fn kelvin_j_arg(x: &ExactNum, p: usize, _cc: &mut Consts) -> ExactComplex {
    // x exp(3πi/4) = x (−√2/2 + i √2/2)
    let two = ExactNum::from_u8(2, p);
    let s = two.sqrt(p, RoundingMode::None);
    // `-h` and `h` are two separate numbers; each is computed (no copy of `h`).
    let re = ExactNum::div_ref(x, &s, p, RoundingMode::None).neg();
    let im = ExactNum::div_ref(x, &s, p, RoundingMode::None);
    ExactComplex::new(re, im)
}

fn kelvin_k_arg(x: &ExactNum, p: usize, _cc: &mut Consts) -> ExactComplex {
    // x exp(πi/4) = x (√2/2 + i √2/2)
    let two = ExactNum::from_u8(2, p);
    let s = two.sqrt(p, RoundingMode::None);
    let h = ExactNum::div_ref(x, &s, p, RoundingMode::None);
    ExactComplex::new(h.clone(), h)
}

#[derive(Clone, Copy)]
enum KelvinWhich {
    Ber,
    Bei,
    Ker,
    Kei,
}

fn signed_inf(neg: bool) -> ExactNum {
    if neg {
        crate::ext::INF_NEG
    } else {
        crate::ext::INF_POS
    }
}

fn plus_zero(p: usize) -> ExactNum {
    ExactNum::new(p)
}

fn is_odd_integer(n: &ExactNum) -> bool {
    if !n.is_int() {
        return false;
    }
    let two = ExactNum::from_u8(2, 64);
    !n.clone().abs().rem(&two).is_zero()
}

/// `|n| mod 8` for an integer `n`.
fn int_abs_mod8(n: &ExactNum) -> Option<u8> {
    if !n.is_int() {
        return None;
    }
    let eight = ExactNum::from_u8(8, 64);
    n.clone()
        .abs()
        .rem(&eight)
        .to_u64(RoundingMode::ToZero)
        .ok()
        .map(|u| u as u8)
}

/// Sign of a clearly nonzero value: `None` if it is zero.
fn nz_sign_neg(x: &ExactNum) -> Option<bool> {
    if x.is_zero() {
        None
    } else {
        Some(x.is_negative())
    }
}

/// Origin, infinity, and out-of-domain edges. `None` means evaluate 10.61 at a
/// positive argument (`x > 0`, or integer-order `ber`/`bei` at `|x|`).
fn kelvin_edge(
    x: &ExactNum,
    nu: &ExactNum,
    which: KelvinWhich,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNum> {
    if !finite(nu) || x.is_nan() {
        return Some(cat_nan());
    }
    if x.is_inf() {
        if x.is_negative() {
            return Some(cat_nan());
        }
        return Some(match which {
            KelvinWhich::Ker | KelvinWhich::Kei => plus_zero(p),
            KelvinWhich::Ber | KelvinWhich::Bei => cat_nan(),
        });
    }
    if x.is_zero() {
        return Some(kelvin_at_zero(nu, which, p, rm, cc));
    }
    if x.is_negative() {
        match which {
            KelvinWhich::Ker | KelvinWhich::Kei => return Some(cat_nan()),
            KelvinWhich::Ber | KelvinWhich::Bei => {
                if !nu.is_int() {
                    return Some(cat_nan());
                }
                return None;
            }
        }
    }
    None
}

fn kelvin_at_zero(
    nu: &ExactNum,
    which: KelvinWhich,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> ExactNum {
    match which {
        KelvinWhich::Ber | KelvinWhich::Bei => kelvin_jb_at_zero(nu, which, p, rm, cc),
        KelvinWhich::Ker | KelvinWhich::Kei => kelvin_k_at_zero(nu, which, p, rm, cc),
    }
}

fn kelvin_jb_at_zero(
    nu: &ExactNum,
    which: KelvinWhich,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> ExactNum {
    let _ = rm;
    if nu.is_zero() {
        return if matches!(which, KelvinWhich::Ber) {
            ExactNum::from_u8(1, p)
        } else {
            plus_zero(p)
        };
    }
    if nu.is_int() || nu.is_positive() {
        return plus_zero(p);
    }
    // ν < 0, non-integer: (z/2)^ν / Γ(ν+1) with z = x e^{3πi/4} → ±Inf.
    let wrk = work_p(p);
    let theta = nu
        .clone()
        .mul(&ExactNum::from_u8(3, wrk), wrk, RoundingMode::None)
        .mul(&cc.pi(wrk, RoundingMode::None), wrk, RoundingMode::None)
        .ldexp(-2, wrk, RoundingMode::None);
    let mut c = theta.cos(wrk, RoundingMode::None, cc);
    let mut s = theta.sin(wrk, RoundingMode::None, cc);
    let g = nu
        .clone()
        .add(&ExactNum::from_u8(1, wrk), wrk, RoundingMode::None)
        .gamma(wrk, RoundingMode::None, cc);
    if g.is_nan() {
        return cat_nan();
    }
    if g.is_negative() {
        c = c.neg();
        s = s.neg();
    }
    let part = if matches!(which, KelvinWhich::Ber) { c } else { s };
    match nz_sign_neg(&part) {
        None => plus_zero(p),
        Some(neg) => signed_inf(neg),
    }
}

fn kelvin_k_at_zero(
    nu: &ExactNum,
    which: KelvinWhich,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> ExactNum {
    if nu.is_zero() {
        return if matches!(which, KelvinWhich::Ker) {
            crate::ext::INF_POS
        } else {
            finish(
                cc.pi(work_p(p), RoundingMode::None)
                    .ldexp(-2, work_p(p), RoundingMode::None)
                    .neg(),
                p,
                rm,
            )
        };
    }
    if let Some(m) = int_abs_mod8(nu) {
        return kelvin_k_at_zero_int(nu.is_negative(), m, which);
    }
    let wrk = work_p(p);
    let pi = cc.pi(wrk, RoundingMode::None);
    let (c, s) = if nu.is_positive() {
        // ker ~ cos(3νπ/4), kei ~ −sin(3νπ/4)
        let theta = nu
            .clone()
            .mul(&ExactNum::from_u8(3, wrk), wrk, RoundingMode::None)
            .mul(&pi, wrk, RoundingMode::None)
            .ldexp(-2, wrk, RoundingMode::None);
        let c = theta.cos(wrk, RoundingMode::None, cc);
        let s = theta.sin(wrk, RoundingMode::None, cc).neg();
        (c, s)
    } else {
        // ν = −μ: ker ~ cos(μπ/4), kei ~ sin(μπ/4)
        let mu = nu.clone().abs();
        let theta = mu
            .mul(&pi, wrk, RoundingMode::None)
            .ldexp(-2, wrk, RoundingMode::None);
        let c = theta.cos(wrk, RoundingMode::None, cc);
        let s = theta.sin(wrk, RoundingMode::None, cc);
        (c, s)
    };
    let part = if matches!(which, KelvinWhich::Ker) { c } else { s };
    match nz_sign_neg(&part) {
        None => plus_zero(p),
        Some(neg) => signed_inf(neg),
    }
}

/// Leading (and, when that part vanishes, k=1) signs of ker_n(0), kei_n(0)
/// for integer n ≠ 0. `m = |n| mod 8`. A `true` sign is negative infinity.
fn kelvin_k_at_zero_int(neg_nu: bool, m: u8, which: KelvinWhich) -> ExactNum {
    let (lead_ker, lead_kei) = if !neg_nu {
        // ker ~ cos(3nπ/4), kei ~ −sin(3nπ/4)
        match m {
            0 => (Some(false), None),
            1 => (Some(true), Some(true)),
            2 => (None, Some(false)),
            3 => (Some(false), Some(true)),
            4 => (Some(true), None),
            5 => (Some(false), Some(false)),
            6 => (None, Some(true)),
            7 => (Some(true), Some(false)),
            _ => return cat_nan(),
        }
    } else {
        // ν = −μ: ker ~ cos(μπ/4), kei ~ sin(μπ/4)
        match m {
            0 => (Some(false), None),
            1 => (Some(false), Some(false)),
            2 => (None, Some(false)),
            3 => (Some(true), Some(false)),
            4 => (Some(true), None),
            5 => (Some(true), Some(true)),
            6 => (None, Some(true)),
            7 => (Some(false), Some(true)),
            _ => return cat_nan(),
        }
    };
    let lead = if matches!(which, KelvinWhich::Ker) { lead_ker } else { lead_kei };
    if let Some(neg) = lead {
        return signed_inf(neg);
    }
    // The k=1 term of K_|n| is rotated by π/2 and carries a minus from (−z²/4).
    let next_neg = if !neg_nu {
        match (which, m) {
            (KelvinWhich::Kei, 0) => true,
            (KelvinWhich::Kei, 4) => false,
            (KelvinWhich::Ker, 2) => true,
            (KelvinWhich::Ker, 6) => false,
            _ => return plus_zero(64),
        }
    } else {
        match (which, m) {
            (KelvinWhich::Kei, 0) => true,
            (KelvinWhich::Kei, 4) => false,
            (KelvinWhich::Ker, 2) => false,
            (KelvinWhich::Ker, 6) => true,
            _ => return plus_zero(64),
        }
    };
    signed_inf(next_neg)
}

/// True when e^{±x/√2} is past the exponent range (DLMF 10.67): e^{-ξ} is below
/// `min_positive` and e^{ξ} is above `EXPONENT_MAX`. That happens once
/// ξ > ((1<<62)−1) ln 2 ≈ 2^{61.47}. `x ≥ 2^{62}` (exponent field ≥ 63) is
/// already past; in the binade `[2^{61}, 2^{62})` the comparison is made
/// against the exact threshold, rounded so the short-circuit fires only when
/// underflow/overflow is certain.
fn kelvin_xi_past_exp_range(x: &ExactNum, cc: &mut Consts) -> bool {
    match x.exponent() {
        Some(e) if e >= 63 => true,
        Some(e) if e < 62 => false,
        Some(_) => {
            let wrk = 128;
            let ln2 = cc.ln_2(wrk, RoundingMode::None);
            let sqrt2 = ExactNum::from_u8(2, wrk).sqrt(wrk, RoundingMode::None);
            // (2^62 − 1) via 2^62 − 1, both exact at this width.
            let two62 = ExactNum::from_u8(1, wrk).ldexp(62, wrk, RoundingMode::None);
            let two62m1 = two62.sub(&ExactNum::from_u8(1, wrk), wrk, RoundingMode::None);
            let thresh =
                sqrt2
                    .mul(&ln2, wrk, RoundingMode::Up)
                    .mul(&two62m1, wrk, RoundingMode::Up);
            matches!(x.compare(&thresh), Some(ord) if ord > 0)
        }
        None => false,
    }
}

/// Far-underflow of a signed value: signed `min_positive` when `rm` rounds
/// away from zero, otherwise a signed zero.
fn kelvin_underflow(neg: bool, p: usize, rm: RoundingMode) -> ExactNum {
    let away = if neg {
        matches!(rm, RoundingMode::Down | RoundingMode::FromZero)
    } else {
        matches!(rm, RoundingMode::Up | RoundingMode::FromZero)
    };
    if away {
        let m = ExactNum::min_positive(p);
        if neg {
            m.neg()
        } else {
            m
        }
    } else {
        let z = ExactNum::new(p);
        if neg {
            z.neg()
        } else {
            z
        }
    }
}

fn kelvin_overflow(neg: bool) -> ExactNum {
    let mut v = signed_inf(neg);
    v.set_inexact(true);
    v
}

/// Sign of the DLMF 10.67 leading trigonometric factor, when that factor is
/// clearly bounded away from zero. `None` if the sign cannot be certified.
fn kelvin_asymp_trig_neg(
    x: &ExactNum,
    nu: &ExactNum,
    which: KelvinWhich,
    cc: &mut Consts,
) -> Option<bool> {
    let xe = x.exponent().unwrap_or(0).max(0) as usize;
    // Sign of cos(x/√2 + ·) at |x| ~ 2^{EXPONENT_MAX} needs ~2^62 bits of 1/π,
    // which is not feasible (and `word_p(xe+256)` would exceed MAX_PRECISION).
    if xe > 200_000 || xe.saturating_add(256) > MAX_PRECISION {
        return None;
    }
    let wrk = word_p(xe.saturating_add(256).max(256));
    if wrk > MAX_PRECISION {
        return None;
    }
    let sqrt2 = ExactNum::from_u8(2, wrk).sqrt(wrk, RoundingMode::None);
    let xi = ExactNum::div_ref(x, &sqrt2, wrk, RoundingMode::None);
    let pi = cc.pi(wrk, RoundingMode::None);
    let half_pi_nu =
        nu.clone()
            .mul(&pi, wrk, RoundingMode::None)
            .ldexp(-1, wrk, RoundingMode::None);
    let eighth = pi.ldexp(-3, wrk, RoundingMode::None);
    // ber/bei: x/√2 + (ν/2 − 1/8)π; ker/kei: x/√2 + (ν/2 + 1/8)π.
    let theta =
        match which {
            KelvinWhich::Ber | KelvinWhich::Bei => xi
                .add(&half_pi_nu, wrk, RoundingMode::None)
                .sub(&eighth, wrk, RoundingMode::None),
            KelvinWhich::Ker | KelvinWhich::Kei => xi
                .add(&half_pi_nu, wrk, RoundingMode::None)
                .add(&eighth, wrk, RoundingMode::None),
        };
    let c = theta.cos(wrk, RoundingMode::None, cc);
    let s = theta.sin(wrk, RoundingMode::None, cc);
    if c.is_nan() || s.is_nan() {
        return None;
    }
    // Leading factor: ber ~ cos, bei ~ sin, ker ~ cos, kei ~ −sin.
    let t = match which {
        KelvinWhich::Ber | KelvinWhich::Ker => c,
        KelvinWhich::Bei => s,
        KelvinWhich::Kei => s.neg(),
    };
    // |t| > 2^{-16}: a reduction error cannot flip a clearly nonzero sign.
    let bound = ExactNum::from_u8(1, wrk).ldexp(-16, wrk, RoundingMode::None);
    let at = t.clone().abs();
    if !matches!(at.compare(&bound), Some(ord) if ord > 0) {
        return None;
    }
    Some(t.is_negative())
}

fn kelvin_past_exp_range(
    x: &ExactNum,
    nu: &ExactNum,
    which: KelvinWhich,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNum> {
    if !kelvin_xi_past_exp_range(x, cc) {
        return None;
    }
    match kelvin_asymp_trig_neg(x, nu, which, cc) {
        Some(neg) => Some(match which {
            KelvinWhich::Ber | KelvinWhich::Bei => kelvin_overflow(neg),
            KelvinWhich::Ker | KelvinWhich::Kei => kelvin_underflow(neg, p, rm),
        }),
        None => match which {
            // Magnitude underflow is certain. Do not return InvalidArgument and
            // do not invent a sign for a directed min_positive.
            KelvinWhich::Ker | KelvinWhich::Kei => {
                let mut z = plus_zero(p);
                z.set_inexact(true);
                Some(z)
            }
            KelvinWhich::Ber | KelvinWhich::Bei => None,
        },
    }
}

/// Order-zero small-`x` leading terms (DLMF 10.31.2 / 10.61). The remainder is
/// `O(x² |ln x|)` for ker/kei and `O(x²)` for ber/bei, which is below a quarter
/// ulp of the leading term once `x.exponent() < −(p/2) − 32`.
fn kelvin_small_x(
    x: &ExactNum,
    nu: &ExactNum,
    which: KelvinWhich,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNum> {
    if !nu.is_zero() {
        return None;
    }
    let e = x.exponent()?;
    let thresh = -((p as i64) / 2) - 32;
    if e >= 0 || e > thresh {
        return None;
    }
    let wp = word_p(p.saturating_add(256));
    let mut v = match which {
        KelvinWhich::Ker => {
            let half = ExactNum::from_u8(1, wp).ldexp(-1, wp, RoundingMode::None);
            let arg = x.clone().mul(&half, wp, RoundingMode::None);
            let ln = arg.ln(wp, RoundingMode::None, cc);
            let g = cc.euler_gamma(wp, RoundingMode::None);
            (-ln).sub(&g, wp, RoundingMode::None)
        }
        KelvinWhich::Kei => {
            let q = cc
                .pi(wp, RoundingMode::None)
                .ldexp(-2, wp, RoundingMode::None);
            -q
        }
        KelvinWhich::Ber => ExactNum::from_u8(1, wp),
        KelvinWhich::Bei => ExactNum::new(wp),
    };
    if v.is_nan() {
        return None;
    }
    let _ = v.set_precision(p, rm);
    Some(v)
}

fn kelvin_eval_1061(
    x: &ExactNum,
    nu: &ExactNum,
    which: KelvinWhich,
    wp: usize,
    cc: &mut Consts,
) -> ExactNum {
    let mut xc = x.clone();
    let _ = xc.set_precision(wp, RoundingMode::None);
    let mut nc = nu.clone();
    let _ = nc.set_precision(wp, RoundingMode::None);
    match which {
        KelvinWhich::Ber | KelvinWhich::Bei => {
            let z = kelvin_j_arg(&xc, wp, cc);
            let n = ExactComplex::from_real(nc, wp);
            let j = z.bessel_j_nu(&n, wp, RoundingMode::None, cc);
            if matches!(which, KelvinWhich::Ber) {
                j.re().clone()
            } else {
                j.im().clone()
            }
        }
        KelvinWhich::Ker | KelvinWhich::Kei => {
            let z = kelvin_k_arg(&xc, wp, cc);
            let n = ExactComplex::from_real(nc, wp);
            let k = z.bessel_k(&n, wp, RoundingMode::None, cc);
            let half_pi = ExactComplex::new(
                ExactNum::new(wp),
                cc.pi(wp, RoundingMode::None)
                    .ldexp(-1, wp, RoundingMode::None)
                    .neg(),
            );
            let phase = n
                .mul(&half_pi, wp, RoundingMode::None)
                .exp(wp, RoundingMode::None, cc);
            let w = phase.mul(&k, wp, RoundingMode::None);
            if matches!(which, KelvinWhich::Ker) {
                w.re().clone()
            } else {
                w.im().clone()
            }
        }
    }
}

/// In-domain (`x > 0`) correctly rounded 10.61 value. No precision or `|x|` cap:
/// past the exponent range the 10.67 leading term is used; otherwise a Ziv
/// loop around complex `J_ν` / `K_ν`, with the usual `2p` retry budget.
fn kelvin_positive(
    x: &ExactNum,
    nu: &ExactNum,
    which: KelvinWhich,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> ExactNum {
    if let Some(v) = kelvin_past_exp_range(x, nu, which, p, rm, cc) {
        return v;
    }
    if let Some(v) = kelvin_small_x(x, nu, which, p, rm, cc) {
        return v;
    }
    let mut p_inc = WORD_BIT_SIZE;
    let Some(mut p_wrk) = p.checked_add(WORD_BIT_SIZE) else {
        return cat_nan();
    };
    p_wrk = word_p(p_wrk);
    loop {
        let Some(p_x) = p_wrk.checked_add(WORD_BIT_SIZE.saturating_mul(2)) else {
            return ExactNum::nan(Some(Error::InvalidArgument));
        };
        let mut v = kelvin_eval_1061(x, nu, which, p_x, cc);
        if v.is_nan() {
            return v;
        }
        if !finite(&v) {
            return v;
        }
        if v.try_set_precision(p, rm, p_wrk) {
            return v;
        }
        if bump_prec_retry(&mut p_wrk, &mut p_inc, p).is_err() {
            return ExactNum::nan(Some(Error::PrecisionRetryExhausted));
        }
    }
}

impl ExactNum {
    /// Scorer \(\mathrm{Gi}(\mathrm{self})\). Entire; NaN in → NaN out.
    /// \(\mathrm{Gi}(+\infty)=+0\) (DLMF 9.12.25). \(\mathrm{Gi}(-\infty)=+0\)
    /// (9.12.11, 9.12.27, 9.7.9–11: the \(x^{-1/4}\) oscillation tends to 0).
    ///
    /// # Precision
    ///
    /// - Algorithm: Taylor series of \(w''-xw=-1/\pi\) with a geometric tail bound and guard
    ///   bits for its cancellation; Watson's lemma (DLMF 9.12.27, remainder below
    ///   \(c_K(2/\sqrt3)^{3K+1}\)) for large `+x`, stopped where that bound is smallest;
    ///   \(\mathrm{Bi}-\mathrm{Hi}\) with Watson's \(\mathrm{Hi}(-y)\) (DLMF 9.12.29) for large `-x`;
    ///   \(1/(\pi x)\) directly once \(x^3\) passes the working precision.
    /// - Bound: Ziv correct-rounding; NaN past `8p + 16384` working bits.
    ///
    /// Accuracy: correctly-rounded
    pub fn scorer_gi(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.scorer_gi_word(q, r, cc)
        })
    }

    /// Word-precision implementation of [`Self::scorer_gi`].
    fn scorer_gi_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return self.clone();
        }
        // DLMF 9.12.25: Gi(+∞)=+0. DLMF 9.12.11 + 9.12.27 + 9.7.9–11:
        // Gi(−∞)=+0 (the x^{-1/4} oscillation tends to 0; no eventual sign).
        if self.is_inf() {
            return plus_zero(p);
        }
        let x = self.clone();
        ziv_round(p, rm, |wp| {
            if !x.is_negative() {
                if !x.is_zero() {
                    if let Some(r) = scorer_watson(&x, false, wp, cc) {
                        return Some(r);
                    }
                }
                return scorer_taylor(&x, false, wp, cc);
            }
            let y = x.clone().neg();
            match scorer_watson(&y, true, wp, cc) {
                Some((h, gh)) => {
                    let bi = x.bi(wp, RoundingMode::None, cc);
                    let eb = bi.exponent()? - wp as i64 + 1;
                    let eh = h.exponent()? - gh - 1;
                    Some(diff_good(&bi, eb, &h, eh, wp))
                }
                None => scorer_taylor(&x, false, wp, cc),
            }
        })
    }

    /// Scorer \(\mathrm{Hi}(\mathrm{self})=\mathrm{Bi}(\mathrm{self})-\mathrm{Gi}(\mathrm{self})\).
    /// \(\mathrm{Hi}(+\infty)=+\infty\) (DLMF 9.12.29). \(\mathrm{Hi}(-\infty)=+0\)
    /// (9.12.27, leading term \(1/(\pi x)\)).
    ///
    /// # Precision
    ///
    /// - Algorithm: Taylor series of \(w''-xw=1/\pi\) (tail bound, cancellation guard);
    ///   Watson's lemma (DLMF 9.12.29, alternating remainder) for large `-x`;
    ///   \(\mathrm{Bi}-\mathrm{Gi}\) with Watson's \(\mathrm{Gi}\) for large `+x`; overflow
    ///   follows \(\mathrm{Bi}\).
    /// - Bound: Ziv correct-rounding; NaN past `8p + 16384` working bits.
    ///
    /// Accuracy: correctly-rounded
    pub fn scorer_hi(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.scorer_hi_word(q, r, cc)
        })
    }

    /// Word-precision implementation of [`Self::scorer_hi`].
    fn scorer_hi_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return self.clone();
        }
        // DLMF 9.12.29: Hi(+∞)=+∞. DLMF 9.12.27: Hi(−∞)=+0 (∼1/(πx)).
        if self.is_inf_pos() {
            return signed_inf(false);
        }
        if self.is_inf_neg() {
            return plus_zero(p);
        }
        if self.is_positive() {
            // Hi = Bi - Gi with 0 < Gi < 1/(pi x): when Bi overflows, so does Hi.
            let b = self.bi(p, rm, cc);
            if b.is_inf() || b.exponent() == Some(crate::defs::EXPONENT_MAX) {
                return b;
            }
        }
        let x = self.clone();
        ziv_round(p, rm, |wp| {
            if !x.is_positive() {
                if !x.is_zero() {
                    if let Some(r) = scorer_watson(&x.clone().neg(), true, wp, cc) {
                        return Some(r);
                    }
                }
                return scorer_taylor(&x, true, wp, cc);
            }
            match scorer_watson(&x, false, wp, cc) {
                Some((g, gg)) => {
                    let bi = x.bi(wp, RoundingMode::None, cc);
                    if !finite(&bi) {
                        return None;
                    }
                    let eb = bi.exponent()? - wp as i64 + 1;
                    let eg = g.exponent()? - gg - 1;
                    Some(diff_good(&bi, eb, &g, eg, wp))
                }
                None => scorer_taylor(&x, true, wp, cc),
            }
        })
    }

    /// Kelvin \(\mathrm{ber}_\nu(\mathrm{self})=\mathrm{Re}\,J_\nu(x e^{3\pi i/4})\) (DLMF 10.61.1).
    ///
    /// The real API follows DLMF 10.61 (`x ≥ 0`) and P20 / DOM-07. At `x = ±0`:
    /// `ber(0) = 1`, `ber_n(0) = 0` for integer `n ≠ 0`, `ber_ν(0) = 0` for `ν > 0`,
    /// and a signed infinity for non-integer `ν < 0` from the small-`z` limit of `J_ν`.
    /// For `x < 0`, the integer reflection `ber_n(−x) = (−1)^n ber_n(x)` (10.61.5);
    /// otherwise `NaN(InvalidArgument)`. `ber_ν(+∞)` has no limit → `NaN(InvalidArgument)`.
    /// In-domain values use complex `J_ν` at the 10.61 argument (Arb `acb_hypgeom_bessel_j`).
    ///
    /// # Precision
    ///
    /// - Algorithm: DLMF 10.61 via complex `J_ν`; DLMF 10.67 leading term when `x/√2`
    ///   is past the exponent range.
    /// - Bound: Ziv on the real part (`prec_retry_budget`; no domain or precision cap).
    ///
    /// Accuracy: correctly-rounded
    pub fn kelvin_ber_nu(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let nu_c = nu.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || nu.inexact_accumulated(),
            p,
            rm,
            |q, r| this.kelvin_ber_nu_word(&nu_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::kelvin_ber_nu`].
    fn kelvin_ber_nu_word(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if let Some(v) = kelvin_edge(self, nu, KelvinWhich::Ber, p, rm, cc) {
            return v;
        }
        let (arg, flip) = if self.is_negative() {
            (self.clone().abs(), is_odd_integer(nu))
        } else {
            (self.clone(), false)
        };
        let y = kelvin_positive(&arg, nu, KelvinWhich::Ber, p, rm, cc);
        if flip {
            y.neg()
        } else {
            y
        }
    }

    /// Kelvin \(\mathrm{bei}_\nu(\mathrm{self})=\mathrm{Im}\,J_\nu(x e^{3\pi i/4})\) (DLMF 10.61.1).
    ///
    /// Domain as [`Self::kelvin_ber_nu`]: `bei(0) = 0`, integer `bei_n(0) = 0`, signed
    /// infinities for non-integer `ν < 0`, the reflection `bei_n(−x) = (−1)^n bei_n(x)`,
    /// and `NaN(InvalidArgument)` for a non-integer negative argument or at `+∞`.
    ///
    /// # Precision
    ///
    /// - Algorithm: DLMF 10.61 via complex `J_ν`; DLMF 10.67 leading term when `x/√2`
    ///   is past the exponent range.
    /// - Bound: Ziv on the imaginary part (`prec_retry_budget`; no domain or precision cap).
    ///
    /// Accuracy: correctly-rounded
    pub fn kelvin_bei_nu(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let nu_c = nu.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || nu.inexact_accumulated(),
            p,
            rm,
            |q, r| this.kelvin_bei_nu_word(&nu_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::kelvin_bei_nu`].
    fn kelvin_bei_nu_word(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if let Some(v) = kelvin_edge(self, nu, KelvinWhich::Bei, p, rm, cc) {
            return v;
        }
        let (arg, flip) = if self.is_negative() {
            (self.clone().abs(), is_odd_integer(nu))
        } else {
            (self.clone(), false)
        };
        let y = kelvin_positive(&arg, nu, KelvinWhich::Bei, p, rm, cc);
        if flip {
            y.neg()
        } else {
            y
        }
    }

    /// Kelvin \(\mathrm{ker}_\nu(\mathrm{self})=\mathrm{Re}\,e^{-\nu\pi i/2}K_\nu(x e^{\pi i/4})\)
    /// (DLMF 10.61.2).
    ///
    /// The real API follows DLMF 10.61 and P20 / DOM-07. `ker(0) = +∞` (a simple pole);
    /// the sign of `ker_ν(0)` is the small-`z` limit of 10.61.2 (`ker_1(0) = −∞`).
    /// `x < 0` is outside the real domain → `NaN(InvalidArgument)`. `ker_ν(+∞) = +0`.
    /// In-domain values use complex `K_ν` at the 10.61 argument (Arb `acb_hypgeom_bessel_k`).
    ///
    /// # Precision
    ///
    /// - Algorithm: DLMF 10.61 via complex `K_ν`; DLMF 10.67 leading term when `x/√2`
    ///   is past the exponent range.
    /// - Bound: Ziv on the real part (`prec_retry_budget`; no domain or precision cap).
    ///
    /// Accuracy: correctly-rounded
    pub fn kelvin_ker_nu(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let nu_c = nu.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || nu.inexact_accumulated(),
            p,
            rm,
            |q, r| this.kelvin_ker_nu_word(&nu_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::kelvin_ker_nu`].
    fn kelvin_ker_nu_word(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if let Some(v) = kelvin_edge(self, nu, KelvinWhich::Ker, p, rm, cc) {
            return v;
        }
        kelvin_positive(self, nu, KelvinWhich::Ker, p, rm, cc)
    }

    /// Kelvin \(\mathrm{kei}_\nu(\mathrm{self})=\mathrm{Im}\,e^{-\nu\pi i/2}K_\nu(x e^{\pi i/4})\)
    /// (DLMF 10.61.2).
    ///
    /// Domain as [`Self::kelvin_ker_nu`], except `kei(0) = −π/4` (finite, not a pole) and
    /// `kei_1(0) = −∞`. `x < 0` → `NaN(InvalidArgument)`. `kei_ν(+∞) = +0`.
    ///
    /// # Precision
    ///
    /// - Algorithm: DLMF 10.61 via complex `K_ν`; DLMF 10.67 leading term when `x/√2`
    ///   is past the exponent range.
    /// - Bound: Ziv on the imaginary part (`prec_retry_budget`; no domain or precision cap).
    ///
    /// Accuracy: correctly-rounded
    pub fn kelvin_kei_nu(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let nu_c = nu.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || nu.inexact_accumulated(),
            p,
            rm,
            |q, r| this.kelvin_kei_nu_word(&nu_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::kelvin_kei_nu`].
    fn kelvin_kei_nu_word(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if let Some(v) = kelvin_edge(self, nu, KelvinWhich::Kei, p, rm, cc) {
            return v;
        }
        kelvin_positive(self, nu, KelvinWhich::Kei, p, rm, cc)
    }

    /// Order-0 Kelvin \(\mathrm{ber}(\mathrm{self})\).
    ///
    /// Accuracy: correctly-rounded
    pub fn kelvin_ber(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.kelvin_ber_word(q, r, cc)
        })
    }

    /// Word-precision implementation of [`Self::kelvin_ber`].
    fn kelvin_ber_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.kelvin_ber_nu(&ExactNum::new(p), p, rm, cc)
    }

    /// Order-0 Kelvin \(\mathrm{bei}(\mathrm{self})\).
    ///
    /// Accuracy: correctly-rounded
    pub fn kelvin_bei(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.kelvin_bei_word(q, r, cc)
        })
    }

    /// Word-precision implementation of [`Self::kelvin_bei`].
    fn kelvin_bei_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.kelvin_bei_nu(&ExactNum::new(p), p, rm, cc)
    }

    /// Order-0 Kelvin \(\mathrm{ker}(\mathrm{self})\).
    ///
    /// Accuracy: correctly-rounded
    pub fn kelvin_ker(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.kelvin_ker_word(q, r, cc)
        })
    }

    /// Word-precision implementation of [`Self::kelvin_ker`].
    fn kelvin_ker_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.kelvin_ker_nu(&ExactNum::new(p), p, rm, cc)
    }

    /// Order-0 Kelvin \(\mathrm{kei}(\mathrm{self})\).
    ///
    /// Accuracy: correctly-rounded
    pub fn kelvin_kei(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.kelvin_kei_word(q, r, cc)
        })
    }

    /// Word-precision implementation of [`Self::kelvin_kei`].
    fn kelvin_kei_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.kelvin_kei_nu(&ExactNum::new(p), p, rm, cc)
    }

    /// Struve \(\mathbf{H}_\nu(\mathrm{self})\).
    ///
    /// # Precision
    ///
    /// - Algorithm: power series (DLMF 11.2.1) with a tail bound and cancellation guard;
    ///   DLMF 11.6.1 (\(Y_\nu+K_\nu\)) only when \(z>0\) and \(m\ge\nu-1/2\), so the remainder
    ///   is of the same sign as and no larger than the first omitted term.
    ///
    /// Accuracy: unspecified
    pub fn struve_h(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let nu_c = nu.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || nu.inexact_accumulated(),
            p,
            rm,
            |q, r| this.struve_h_word(&nu_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::struve_h`].
    fn struve_h_word(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        catalog_aws::struve_h_word(self, nu, p, rm, cc)
    }

    /// Anger \(\mathbf{J}_\nu(\mathrm{self})\). Integer \(\nu\) reduces to Bessel \(J\).
    ///
    /// # Precision
    ///
    /// - Algorithm: Maclaurin series (DLMF 11.10.8–11) with a tail bound and cancellation
    ///   guard; for large \(z>0\), DLMF 11.10.15 with \(A_\nu\) from the Laplace / Watson
    ///   expansion (DLMF 11.11.4) and a Cauchy remainder, else Taylor panels.
    ///   Negative \(z\) uses DLMF 11.10.12. \(z=0\) uses the closed form \(\sin(\nu\pi)/(\nu\pi)\).
    ///
    /// Accuracy: unspecified
    pub fn anger_j(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let nu_c = nu.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || nu.inexact_accumulated(),
            p,
            rm,
            |q, r| this.anger_j_word(&nu_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::anger_j`].
    fn anger_j_word(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        catalog_aws::anger_j_word(self, nu, p, rm, cc)
    }

    /// Weber \(\mathbf{E}_\nu(\mathrm{self})\).
    ///
    /// # Precision
    ///
    /// - Algorithm: Maclaurin series (DLMF 11.10.8–11) with a tail bound and cancellation
    ///   guard; for large \(z>0\), DLMF 11.10.16 with \(A_\nu\) from the Laplace / Watson
    ///   expansion (DLMF 11.11.4) and a Cauchy remainder, else Taylor panels.
    ///   Negative \(z\) uses DLMF 11.10.12. \(z=0\) uses the closed form
    ///   \((1-\cos(\nu\pi))/(\nu\pi)\).
    ///
    /// Accuracy: unspecified
    pub fn weber_e(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let nu_c = nu.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || nu.inexact_accumulated(),
            p,
            rm,
            |q, r| this.weber_e_word(&nu_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::weber_e`].
    fn weber_e_word(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        catalog_aws::weber_e_word(self, nu, p, rm, cc)
    }

    /// Clausen \(\mathrm{Cl}_2(\mathrm{self})=\sum_{k\ge1}\sin(k\theta)/k^2=-\int_0^\theta\ln\lvert 2\sin(t/2)\rvert\,dt\).
    ///
    /// Accuracy: unspecified
    pub fn clausen_cl2(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.clausen_cl2_word(q, r, cc)
        })
    }

    /// Word-precision implementation of [`Self::clausen_cl2`].
    fn clausen_cl2_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) {
            return cat_nan();
        }
        if self.is_zero() {
            return finish(self.clone(), p, rm);
        }
        let wrk = work_p(p);
        let Some((th, u)) = clausen_reduce(self, wrk, cc) else {
            return cat_nan();
        };
        let bits = word_p(p).saturating_add(32);
        let third =
            cc.pi(wrk, RoundingMode::None)
                .div(&ExactNum::from_u8(3, wrk), wrk, RoundingMode::None);
        if matches!(u.compare(&third), Some(c) if c >= 0) {
            // |θ| ≤ 2π/3: the Bernoulli series directly (θ carries wrk relative bits).
            return finish(clausen_cl2_series(&th, wrk, bits, cc), p, rm);
        }
        // |θ| = π − u, u < π/3: Cl2(π − u) = Cl2(u) − Cl2(2u)/2 (odd in θ). The u ln u parts cancel to
        // u ln 2, losing log2|ln u| bits; the guard covers them so the zero at ±π keeps relative accuracy.
        let eu = u.exponent().unwrap_or(0).min(0).unsigned_abs() as usize;
        let pw = word_p(
            wrk.saturating_add(64)
                .saturating_add(usize::BITS as usize - eu.leading_zeros() as usize),
        );
        let two_u = ExactNum::from_u8(2, pw).mul(&u, pw, RoundingMode::None);
        let a = clausen_cl2_series(&u, pw, bits.saturating_add(pw - wrk), cc);
        let b = clausen_cl2_series(&two_u, pw, bits.saturating_add(pw - wrk), cc);
        let v = a.sub(
            &b.div(&ExactNum::from_u8(2, pw), pw, RoundingMode::None),
            pw,
            RoundingMode::None,
        );
        finish(if th.is_negative() { v.neg() } else { v }, p, rm)
    }

    /// Clausen \(\mathrm{Cl}_3(\mathrm{self})=\sum_{k\ge1}\cos(k\theta)/k^3=\zeta(3)-\int_0^\theta\mathrm{Cl}_2(t)\,dt\).
    ///
    /// Accuracy: unspecified
    pub fn clausen_cl3(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.clausen_cl3_word(q, r, cc)
        })
    }

    /// Word-precision implementation of [`Self::clausen_cl3`].
    fn clausen_cl3_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) {
            return cat_nan();
        }
        let wrk = work_p(p);
        if self.is_zero() {
            let z3 = zeta3(wrk, cc);
            return finish(z3, p, rm);
        }
        // Cl3 has zeros inside (0, π): the series' absolute error is widened until it sits
        // below the value's own last bits.
        let mut extra = 0usize;
        loop {
            let w = word_p(wrk.saturating_add(extra));
            let Some((th, u)) = clausen_reduce(self, w, cc) else {
                return cat_nan();
            };
            let bits = word_p(p).saturating_add(32).saturating_add(extra);
            let third =
                cc.pi(w, RoundingMode::None)
                    .div(&ExactNum::from_u8(3, w), w, RoundingMode::None);
            let v = if matches!(u.compare(&third), Some(c) if c >= 0) {
                clausen_cl3_series(&th.abs(), w, bits, cc)
            } else {
                // Cl3(π − u) = Cl3(2u)/4 − Cl3(u), both arguments ≤ 2π/3.
                let two_u = ExactNum::from_u8(2, w).mul(&u, w, RoundingMode::None);
                clausen_cl3_series(&two_u, w, bits, cc)
                    .div(&ExactNum::from_u8(4, w), w, RoundingMode::None)
                    .sub(&clausen_cl3_series(&u, w, bits, cc), w, RoundingMode::None)
            };
            if v.is_nan() {
                return v;
            }
            let need = match v.exponent() {
                Some(e) if !v.is_zero() => (-e).max(0) as usize + 16,
                _ => extra.saturating_mul(2).max(64),
            };
            if extra >= need {
                return finish(v, p, rm);
            }
            extra = need;
        }
    }

    /// Hurwitz \(\zeta(s,a)\) for integer \(s\ge 2\) and \(a>0\). `self` is \(s\).
    ///
    /// Accuracy: unspecified
    pub fn hurwitz_zeta(&self, a: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let a_c = a.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || a.inexact_accumulated(),
            p,
            rm,
            |q, r| this.hurwitz_zeta_word(&a_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::hurwitz_zeta`].
    fn hurwitz_zeta_word(&self, a: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) || !finite(a) {
            return cat_nan();
        }
        let wrk = work_p(p);
        let Some(s) = positive_int_u64(self).filter(|&s| s >= 2) else {
            return cat_nan();
        };
        let Some(s_i) = i64::try_from(s).ok() else {
            return cat_nan();
        };
        let mut aa = a.clone();
        let mut extra = ExactNum::new(wrk);
        while matches!(aa.compare(&ExactNum::from_u8(1, wrk)), Some(c) if c < 0) {
            if aa.is_zero() || (aa.is_int() && !aa.is_positive()) {
                return cat_nan();
            }
            extra = extra.add(
                &ExactNum::from_u8(1, wrk).div(
                    &aa.powsi(s_i, wrk, RoundingMode::None),
                    wrk,
                    RoundingMode::None,
                ),
                wrk,
                RoundingMode::None,
            );
            aa = aa.add(&ExactNum::from_u8(1, wrk), wrk, RoundingMode::None);
        }
        finish(
            hurwitz_em(s, &aa, wrk, cc).add(&extra, wrk, RoundingMode::None),
            p,
            rm,
        )
    }

    /// Riemann \(\zeta(s)\) for real \(s\). \(\zeta(-2k)=+0\); the pole \(s=1\) is
    /// `NaN(InvalidArgument)`. `self` is \(s\).
    ///
    /// Accuracy: unspecified
    pub fn riemann_zeta(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.riemann_zeta_word(q, r, cc)
        })
    }

    /// Word-precision implementation of [`Self::riemann_zeta`].
    fn riemann_zeta_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let wrk = work_p(p);
        if let Some(s) = small_int(self, wrk, -1, 64) {
            if s == 0 {
                let half = ExactNum::from_u8(1, wrk)
                    .div(&ExactNum::from_u8(2, wrk), wrk, RoundingMode::None)
                    .neg();
                return finish(half, p, rm);
            }
            if s == -1 {
                let v = ExactNum::from_u8(1, wrk)
                    .div(&ExactNum::from_u32(12, wrk), wrk, RoundingMode::None)
                    .neg();
                return finish(v, p, rm);
            }
            if s >= 2 && s % 2 == 0 {
                return finish(zeta_even((s / 2) as usize, wrk, cc), p, rm);
            }
        }
        self.zeta_real_word(p, rm, cc)
    }

    /// Polygamma \(\psi^{(n)}(\mathrm{self})\). \(n=0\) is digamma; \(n\ge1\) uses
    /// \(\psi^{(n)}(z)=(-1)^{n+1}n!\,\zeta(n+1,z)\).
    ///
    /// Accuracy: unspecified
    pub fn polygamma(&self, n: usize, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.polygamma_word(n, q, r, cc)
        })
    }

    /// Word-precision implementation of [`Self::polygamma`].
    fn polygamma_word(&self, n: usize, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) {
            return cat_nan();
        }
        if n == 0 {
            return self.digamma(p, rm, cc);
        }
        let wrk = work_p(p);
        let Some(s) = u64::try_from(n).ok().and_then(|n| n.checked_add(1)) else {
            return cat_nan();
        };
        let z = ExactNum::from_u64(s, wrk).hurwitz_zeta(self, wrk, RoundingMode::None, cc);
        let fact = factorial_n(n, wrk);
        let mut y = fact.mul(&z, wrk, RoundingMode::None);
        // ψ^{(n)} = (-1)^{n+1} n! ζ(n+1, z)
        if n.is_multiple_of(2) {
            y = y.neg();
        }
        finish(y, p, rm)
    }

    /// Barnes \(G(\mathrm{self})\). Entire; zeros at non-positive integers.
    ///
    /// Accuracy: unspecified
    pub fn barnes_g(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.barnes_g_word(q, r, cc)
        })
    }

    /// Word-precision implementation of [`Self::barnes_g`].
    fn barnes_g_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        // G(+∞)=+∞ (the limit is exact, every mode). G(−∞) has no limit
        // (zeros at every non-positive integer).
        if self.is_inf_pos() {
            return crate::ext::INF_POS;
        }
        if !finite(self) {
            return cat_nan();
        }
        // Large positive z: G overflows the exponent range. C99 7.12.1 maps
        // that by rounding mode (+MAX toward zero, +Inf otherwise), not a
        // refusal and not +Inf in every mode.
        if self.is_positive() && barnes_g_pos_overflows(self) {
            return catalog_aws::overflow_signed(false, p, rm);
        }
        let wrk = work_p(p);
        // Zeros at every non-positive integer, however large.
        if self.is_int() && (self.is_zero() || self.is_negative()) {
            return finish(ExactNum::new(wrk), p, rm);
        }
        // Tiny |z|: G(z) = z + c z² + O(z³), c = γ + (ln 2π − 1)/2 > 0; exp(ln G) would underflow.
        // The z² term is formed rounding up so that, when it is below z's last bit, it still
        // moves the sum above z for the directed modes (G(z) > z for every real z ≠ 0 here).
        if matches!(self.exponent(), Some(e) if e < -((p / 2) as crate::defs::Exponent) - 4) {
            let up = RoundingMode::Up;
            let two_pi = ExactNum::from_u8(2, wrk).mul(
                &cc.pi(wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            );
            let half_ln = two_pi
                .ln(wrk, RoundingMode::None, cc)
                .sub(&ExactNum::from_u8(1, wrk), wrk, RoundingMode::None)
                .ldexp(-1, wrk, RoundingMode::None);
            let c = cc
                .euler_gamma(wrk, RoundingMode::None)
                .add(&half_ln, wrk, RoundingMode::None);
            let z2c = self.sqr(wrk, up).mul(&c, wrk, up);
            return finish(ExactNum::add_ref(self, &z2c, wrk, up), p, rm);
        }
        if let Some(n) = small_int(self, wrk, -64, 64) {
            if n <= 0 {
                return finish(ExactNum::new(wrk), p, rm);
            }
            let mut g = ExactNum::from_u8(1, wrk);
            for k in 2..n {
                g = g.mul(&factorial_n((k - 1) as usize, wrk), wrk, RoundingMode::None);
            }
            return finish(g, p, rm);
        }
        if !self.is_positive() && self.compare(&ExactNum::from_i32(-64, wrk)) == Some(-1) {
            return barnes_g_reflect(self, p, rm, cc);
        }
        if !self.is_positive() {
            // Recur up to a positive argument: G(z) = G(z+1)/Γ(z) (at most 64 levels; below
            // −64 the reflection formula is used, so the depth never grows with |z|).
            let zp = ExactNum::from_u8(1, wrk).add(self, wrk, RoundingMode::None);
            let g1 = zp.barnes_g(wrk, RoundingMode::None, cc);
            let gz = self.gamma(wrk, RoundingMode::None, cc);
            return finish(g1.div(&gz, wrk, RoundingMode::None), p, rm);
        }
        let g = barnes_ln_g(self, wrk, cc).exp(wrk, RoundingMode::None, cc);
        // The lb≥30 bound is sufficient but not tight: z ∈ [2^{29}, 2^{30}) can
        // still overflow. If exp(ln G) overflows, map by C99 7.12.1.
        if g.is_inf() {
            return catalog_aws::overflow_signed(g.is_negative(), p, rm);
        }
        finish(g, p, rm)
    }

    /// Inverse Jacobi \(\mathrm{arcsn}(\mathrm{self}|m)=F(\mathrm{self}|m)\).
    ///
    /// Accuracy: correctly-rounded
    pub fn jacobi_arcsn(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let m_c = m.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || m.inexact_accumulated(),
            p,
            rm,
            |q, r| this.jacobi_arcsn_word(&m_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::jacobi_arcsn`].
    fn jacobi_arcsn_word(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.elliptic_f(m, p, rm, cc)
    }

    /// Inverse Jacobi \(\mathrm{arccn}(x|m)=F(\sqrt{1-x^2}|m)\).
    ///
    /// Accuracy: correctly-rounded
    pub fn jacobi_arccn(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let m_c = m.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || m.inexact_accumulated(),
            p,
            rm,
            |q, r| this.jacobi_arccn_word(&m_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::jacobi_arccn`].
    fn jacobi_arccn_word(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) || !finite(m) {
            return cat_nan();
        }
        self.jacobi_arc_cd_word(m, false, p, rm, cc)
    }

    /// Inverse Jacobi \(\mathrm{arcdn}(x|m)=F(\sqrt{(1-x^2)/m}|m)\) for \(m>0\).
    ///
    /// Accuracy: correctly-rounded
    pub fn jacobi_arcdn(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let m_c = m.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || m.inexact_accumulated(),
            p,
            rm,
            |q, r| this.jacobi_arcdn_word(&m_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::jacobi_arcdn`].
    fn jacobi_arcdn_word(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) || !finite(m) || !m.is_positive() {
            return cat_nan();
        }
        self.jacobi_arc_cd_word(m, true, p, rm, cc)
    }

    /// \({}_0F_1(;b;z)\) with `self = z`.
    ///
    /// Accuracy: unspecified
    pub fn hypergeom_0f1(&self, b: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        let b_c = b.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || b.inexact_accumulated(),
            p,
            rm,
            |q, r| this.hypergeom_0f1_word(&b_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::hypergeom_0f1`].
    fn hypergeom_0f1_word(&self, b: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.hypergeom_pfq(&[], core::slice::from_ref(b), p, rm, cc)
    }

    /// Kummer \({}_1F_1(a;b;z)\) with `self = z`.
    ///
    /// Accuracy: unspecified
    pub fn hypergeom_1f1(
        &self,
        a: &Self,
        b: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        let this = self.exact_copy_ext();
        let a_c = a.exact_copy_ext();
        let b_c = b.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated() || a.inexact_accumulated() || b.inexact_accumulated(),
            p,
            rm,
            |q, r| this.hypergeom_1f1_word(&a_c, &b_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::hypergeom_1f1`].
    fn hypergeom_1f1_word(
        &self,
        a: &Self,
        b: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        self.hypergeom_pfq(
            core::slice::from_ref(a),
            core::slice::from_ref(b),
            p,
            rm,
            cc,
        )
    }

    /// Generalized \({}_pF_q(a_1,\ldots,a_p;b_1,\ldots,b_q;z)\) with `self = z`.
    ///
    /// Direct series with a geometric tail bound; no fixed term cap. Nonterminating
    /// `p > q+1` is `NaN(InvalidArgument)`. `p = q+1, |z| > 1` is continued through Meijer
    /// G. On the cut `z > 1` the value is not real (the limits from above and below are
    /// complex conjugates), so the real function is undefined there and returns NaN, as
    /// FLINT's real `arb_hypgeom_pfq` and C99 domain errors do. `₂F₁` at `z = ±1` is
    /// computed ([`hyp2f1_unit`]); `q ≥ 2` at `|z| = 1` is still `NaN`.
    /// Poles in any `b_j = 0, -1, ...`.
    ///
    /// Accuracy: unspecified
    pub fn hypergeom_pfq(
        &self,
        a: &[Self],
        b: &[Self],
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        let this = self.exact_copy_ext();
        let a_c: Vec<Self> = a.iter().map(|v| v.exact_copy_ext()).collect();
        let b_c: Vec<Self> = b.iter().map(|v| v.exact_copy_ext()).collect();
        Self::round_exact(
            self.inexact_accumulated()
                || a.iter().any(|v| v.inexact_accumulated())
                || b.iter().any(|v| v.inexact_accumulated()),
            p,
            rm,
            |q, r| this.hypergeom_pfq_word(&a_c, &b_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::hypergeom_pfq`].
    fn hypergeom_pfq_word(
        &self,
        a: &[Self],
        b: &[Self],
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        let _ = cc;
        if !finite(self) || a.iter().any(|x| !finite(x)) || b.iter().any(|x| !finite(x)) {
            return cat_nan();
        }
        // DLMF 16.2.5: 1F0(a;;z) = (1-z)^{-a}. On the cut z>1 this is real precisely
        // when a is an integer; otherwise the real API refuses, matching ZF-38.
        if a.len() == 1 && b.is_empty() {
            return finish(pfq_1f0(self, &a[0], p, cc), p, rm);
        }
        let terminating = a.iter().any(nonpositive_int);
        if !terminating && !self.is_zero() {
            if a.len() > b.len() + 1 {
                return cat_nan();
            }
            if a.len() == b.len() + 1 {
                match self.abs_cmp(&ExactNum::from_u8(1, WORD_BIT_SIZE)) {
                    Some(c) if c > 0 => return finish(pfq_continued(self, a, b, p, cc), p, rm),
                    Some(0) if a.len() == 2 => {
                        if nonpositive_int(&b[0]) {
                            return cat_nan();
                        }
                        return finish(hyp2f1_unit(self, a, &b[0], p, cc), p, rm);
                    }
                    Some(0) => return cat_nan(),
                    _ => {}
                }
            }
        }
        for bj in b {
            if nonpositive_int(bj) {
                return cat_nan();
            }
        }
        finish(pfq_series(self, a, b, p), p, rm)
    }

    /// Principal Lambert \(W_0(\mathrm{self})\). Domain \([-1/e,\infty)\).
    ///
    /// Accuracy: correctly-rounded
    pub fn lambert_w0(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.lambert_w0_word(q, r, cc)
        })
    }

    /// Word-precision implementation of [`Self::lambert_w0`].
    fn lambert_w0_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.lambert_w_branch(false, p, rm, cc)
    }

    /// Lambert \(W_{-1}(\mathrm{self})\). Domain \([-1/e,0)\).
    /// \(\pm 0\) and \(x>0\) and \(x<-1/e\) give NaN(InvalidArgument).
    ///
    /// Accuracy: correctly-rounded
    pub fn lambert_wm1(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.lambert_wm1_word(q, r, cc)
        })
    }

    /// Word-precision implementation of [`Self::lambert_wm1`].
    fn lambert_wm1_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.lambert_w_branch(true, p, rm, cc)
    }

    fn lambert_w_branch(&self, minus1: bool, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.lambert_w_word(minus1, p, rm, cc)
    }

    /// Polylogarithm \(\mathrm{Li}_n(\mathrm{self})\) for integer \(n\ge 2\) and \(\lvert x\rvert\le 1\).
    ///
    /// Accuracy: unspecified
    pub fn polylog(&self, n: usize, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy_ext();
        Self::round_exact(self.inexact_accumulated(), p, rm, |q, r| {
            this.polylog_word(n, q, r, cc)
        })
    }

    /// Defining series Li_n(z) = ∑_{k≥1} z^k / k^n, no fixed term cap (P9).
    /// Tail after term t_K is at most |t_K| |z| / (1-|z|) when |z|<1.
    fn polylog_series(z: &ExactNum, n: usize, p: usize) -> ExactNum {
        let wrk = work_p(p);
        let nr = RoundingMode::None;
        let mut xk = z.clone();
        let mut sum = ExactNum::new(wrk);
        let mut last_e = None;
        let mut k: u64 = 1;
        loop {
            let den = ExactNum::from_u64(k, wrk).powsi(n as i64, wrk, nr);
            let term = ExactNum::div_ref(&xk, &den, wrk, nr);
            sum = sum.add(&term, wrk, nr);
            if under_word(&term, p) {
                break;
            }
            let e = term.exponent();
            if last_e.is_some_and(|le| e.is_some_and(|ee| ee >= le)) && k > 4 {
                // Term magnitude stalled (should not happen for |z|<1); stop rather than loop forever.
                break;
            }
            last_e = e;
            xk = xk.mul(z, wrk, nr);
            k = k.saturating_add(1);
        }
        sum
    }

    /// DLMF 25.12.2 at working precision. `omz` is 1-z, already |1-z| ≤ |z|.
    fn polylog_dilog_identity(z: &ExactNum, omz: &ExactNum, p: usize, cc: &mut Consts) -> ExactNum {
        let wrk = work_p(p).saturating_add(WORD_BIT_SIZE);
        let nr = RoundingMode::None;
        let li2 = Self::polylog_series(omz, 2, wrk);
        let ln_z = z.ln(wrk, nr, cc);
        let ln_omz = omz.ln(wrk, nr, cc);
        let zeta2 = ExactNum::from_u8(2, wrk).riemann_zeta(wrk, nr, cc);
        zeta2
            .sub(&li2, wrk, nr)
            .sub(&ln_z.mul(&ln_omz, wrk, nr), wrk, nr)
    }

    /// Word-precision implementation of [`Self::polylog`].
    fn polylog_word(&self, n: usize, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) || !(2..=64).contains(&n) {
            return cat_nan();
        }
        let wrk = work_p(p);
        let one = ExactNum::from_u8(1, wrk);
        if self.compare(&one) == Some(0) {
            return ExactNum::from_u32(n as u32, wrk).riemann_zeta(p, rm, cc);
        }
        if matches!(self.abs_cmp(&one), Some(c) if c > 0) {
            return cat_nan();
        }
        if self.compare(&ExactNum::from_i8(-1, wrk)) == Some(0) {
            let z = ExactNum::from_u32(n as u32, wrk).riemann_zeta(wrk, RoundingMode::None, cc);
            let two = ExactNum::from_u8(2, wrk);
            let f = two.powsi(1 - (n as i64), wrk, RoundingMode::None).sub(
                &one,
                wrk,
                RoundingMode::None,
            );
            return finish(f.mul(&z, wrk, RoundingMode::None), p, rm);
        }
        // DLMF 25.12.2: Li_2(z) + Li_2(1-z) = ζ(2) − ln(z) ln(1-z).
        // Use it when 0 < z < 1 and |1-z| ≤ |z|, so the companion series is the
        // faster one and we do not recurse. The defining series at z = 1-2^{-100}
        // needs ~2^{100} terms; the identity evaluates Li_2(2^{-100}) in a few.
        if n == 2 && self.is_positive() {
            let omz = ExactNum::sub_ref(&one, self, wrk, RoundingMode::None);
            if finite(&omz) && !omz.is_zero() && matches!(omz.abs_cmp(self), Some(c) if c <= 0) {
                return finish(Self::polylog_dilog_identity(self, &omz, p, cc), p, rm);
            }
        }
        finish(Self::polylog_series(self, n, p), p, rm)
    }
}

/// Proven overflow of Barnes G at a positive `z`.
///
/// DLMF 5.17.5: `ln G(w+1) = (w²/2) ln w − 3w²/4 + (w/2) ln(2π) − (ln w)/12
/// + ζ'(−1) + R(w)`. Hence `ln G(w+1) − ((w²/2) ln w − w²)` equals
/// `w²/4 + (w/2) ln(2π) − (ln w)/12 + ζ'(−1) + R(w)`. For `w ≥ 2^{29}`,
/// `|R(w)| ≤ 1/(120 w²)` and that difference is positive, so
/// `ln G(w+1) > (w²/2) ln w − w²`. If `z ∈ [2^{lb}, 2^{lb+1})` with `lb ≥ 30`,
/// then `w = z−1 ≥ 2^{30}−1 > 2^{29}` and `(w²/2) ln w − w² > EMAX ln 2`, so
/// `G(z)` overflows the exponent range.
fn barnes_g_pos_overflows(z: &ExactNum) -> bool {
    let Some(e) = z.exponent() else {
        return false;
    };
    // `exponent()` e means the value lies in [2^{e-1}, 2^e); lb = floor(log2 z) = e-1.
    e.checked_sub(1).is_some_and(|lb| lb >= 30)
}

/// Barnes G at a non-integer `w < −64` by reflection (ZF-72), with `z = 1 − w > 65`:
/// `ln|G(1−z)| = ln G(1+z) − z ln 2π + z ln|2 sin πz| + Cl₂(2πz)/(2π)`, from
/// `ln G(1−z) − ln G(1+z) = −z ln 2π + ∫₀^z πt cot πt dt` (DLMF 5.17) and
/// `∫₀^z πt cot πt dt = z ln|2 sin πz| + Cl₂(2πz)/(2π)` (both sides have derivative πz cot πz
/// and agree on every interval; checked against mpmath barnesg). sin and Cl₂ take the
/// exactly reduced `z mod 2` and `2π (z mod 1)`. Sign: with `n = ⌊z⌋`, `G(w)` has the sign
/// `(−1)^(n(n+1)/2)` (product of the Γ signs on the way up), negative iff `z mod 4 ∈ [1, 3)`.
/// The former one-step-per-unit recursion `G(z) = G(z+1)/Γ(z)` overflowed the stack and grew
/// the working precision by a word per level.
fn barnes_g_reflect(w: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> ExactNum {
    let nr = RoundingMode::None;
    let wrk = work_p(p);
    let exact_w = word_p(
        wrk.max(w.precision().unwrap_or(0))
            .saturating_add(2 * WORD_BIT_SIZE),
    );
    let z = ExactNum::from_u8(1, exact_w).sub(w, exact_w, nr);
    let ez = z.exponent().unwrap_or(0).max(0) as usize;
    // ln G(1+z) ≈ z² ln z / 2: carry its magnitude in extra bits.
    let pw = word_p(wrk.saturating_add(2 * ez).saturating_add(2 * WORD_BIT_SIZE));
    let t2 = z.rem(&ExactNum::from_u8(2, exact_w));
    let f1 = z.rem(&ExactNum::from_u8(1, exact_w));
    let t4 = z.rem(&ExactNum::from_u8(4, exact_w));
    let neg = matches!(t4.compare(&ExactNum::from_u8(1, exact_w)), Some(c) if c >= 0)
        && matches!(t4.compare(&ExactNum::from_u8(3, exact_w)), Some(c) if c < 0);
    let pi = cc.pi(pw, nr);
    let two_pi = ExactNum::from_u8(2, pw).mul(&pi, pw, nr);
    let s = pi.mul(&t2, pw, nr).sin(pw, nr, cc).abs();
    let ln2s = ExactNum::from_u8(2, pw).mul(&s, pw, nr).ln(pw, nr, cc);
    let cl = ExactNum::mul_ref(&two_pi, &f1, pw, nr)
        .clausen_cl2(pw, nr, cc)
        .div(&two_pi, pw, nr);
    let lng1 = barnes_ln_g(
        &ExactNum::add_ref(&z, &ExactNum::from_u8(1, pw), pw, nr),
        pw,
        cc,
    );
    let ln2pi = ExactNum::clone(&two_pi).ln(pw, nr, cc);
    let a = lng1
        .sub(&ExactNum::mul_ref(&z, &ln2pi, pw, nr), pw, nr)
        .add(&ExactNum::mul_ref(&z, &ln2s, pw, nr), pw, nr)
        .add(&cl, pw, nr);
    let mag_rm = if neg { crate::common::util::invert_rm_for_sign(rm) } else { rm };
    let mut g = a.exp(p, mag_rm, cc);
    if neg {
        g = g.neg();
    }
    g
}

/// ln G(z0) for z0 > 0: shift up until Stirling (DLMF 5.17.5) is past the working width, then subtract the Γ
/// product. NaN if the asymptotic series did not reach the width (never a truncated value).
fn barnes_ln_g(z0: &ExactNum, wrk: usize, cc: &mut Consts) -> ExactNum {
    // FLINT `_acb_log_barnes_g_zeta`: LogG(z)=(z−1)LogΓ(z) − ζ'(−1,z) + ζ'(−1).
    // The former shift-to-height 32+wrk/5 plus a 13,000-term log product cannot
    // finish at dest p=2^16 inside ZFS_WALL. The EM jet is bound-driven.
    crate::hurwitz_jet::log_barnes_g_zeta(z0, wrk, cc)
}

#[cfg(test)]
fn barnes_ln_g_plus1(z: &ExactNum, p: usize, cc: &mut Consts) -> ExactNum {
    // DLMF 5.17.5: ln G(z+1) ∼ (z²/2) ln z − (3/4) z² + (z/2) ln(2π) − (1/12) ln z
    //              + ζ'(−1) + Σ_{k≥1} B_{2k+2} / (4k(k+1) z^{2k}),  ζ'(−1)=1/12−ln A.
    let two = ExactNum::from_u8(2, p);
    let three = ExactNum::from_u8(3, p);
    let four = ExactNum::from_u8(4, p);
    let twelve = ExactNum::from_u32(12, p);
    let z2 = z.sqr(p, RoundingMode::None);
    let lnz = z.ln(p, RoundingMode::None, cc);
    let pi = cc.pi(p, RoundingMode::None);
    let twopi = ExactNum::from_u8(2, p).mul(&pi, p, RoundingMode::None);
    let t1 = ExactNum::div_ref(&z2, &two, p, RoundingMode::None).mul(&lnz, p, RoundingMode::None);
    let t2 = three
        .mul(&z2, p, RoundingMode::None)
        .div(&four, p, RoundingMode::None);
    let t3 = ExactNum::div_ref(z, &two, p, RoundingMode::None).mul(
        &twopi.ln(p, RoundingMode::None, cc),
        p,
        RoundingMode::None,
    );
    let t4 = lnz.div(&twelve, p, RoundingMode::None);
    let lna = ln_glaisher(p, cc);
    let zeta_pm1 = ExactNum::from_u8(1, p)
        .div(&twelve, p, RoundingMode::None)
        .sub(&lna, p, RoundingMode::None);
    let mut s = t1
        .sub(&t2, p, RoundingMode::None)
        .add(&t3, p, RoundingMode::None)
        .sub(&t4, p, RoundingMode::None)
        .add(&zeta_pm1, p, RoundingMode::None);
    // Sum until a term falls below the last working bit of s (no fixed term count: a cap of 48
    // terms stopped near 2^-400 at z ≈ 100, a silent error at p ≥ 512). The series is asymptotic:
    // if the terms turn upward first, z was too small for this width and the result is NaN.
    let mut bs: Vec<ExactNum> = Vec::new();
    let mut z2k = ExactNum::from_u8(1, p);
    let mut last: Option<crate::Exponent> = None;
    let mut k = 1usize;
    loop {
        if k + 1 > bs.len() {
            bs = even_bernoulli_ext((2 * (k + 1)).max(32), p);
        }
        let b = &bs[k];
        if b.is_nan() {
            return cat_nan();
        }
        z2k = z2k.mul(&z2, p, RoundingMode::None);
        let den = ExactNum::from_u8(4, p)
            .mul(&ExactNum::from_u32(k as u32, p), p, RoundingMode::None)
            .mul(
                &ExactNum::from_u32((k + 1) as u32, p),
                p,
                RoundingMode::None,
            )
            .mul(&z2k, p, RoundingMode::None);
        let term = ExactNum::div_ref(b, &den, p, RoundingMode::None);
        let et = term.exponent().unwrap_or(crate::Exponent::MIN / 2);
        s = s.add(&term, p, RoundingMode::None);
        if term.is_zero() || et + (p as crate::Exponent) + 16 < s.exponent().unwrap_or(0) {
            return s;
        }
        if matches!(last, Some(l) if et > l) {
            return cat_nan();
        }
        last = Some(et);
        k += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(a: &ExactNum, b: &ExactNum, p: usize) -> bool {
        let d = ExactNum::sub_ref(a, b, p, RoundingMode::None).abs();
        d.is_zero()
            || d.exponent()
                .is_some_and(|e| e < -((p as crate::Exponent) / 4))
    }

    #[test]
    fn test_catalog_identities() -> Result<(), crate::Error> {
        let p = 128;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);
        let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, rm);

        let gi0 = zero.scorer_gi(p, rm, &mut cc);
        let bi0 = zero.bi(p, rm, &mut cc);
        assert!(near(
            &gi0,
            &ExactNum::div_ref(&bi0, &ExactNum::from_u8(3, p), p, rm),
            p
        ));
        let hi0 = zero.scorer_hi(p, rm, &mut cc);
        assert!(near(
            &hi0,
            &bi0.mul(&ExactNum::from_u8(2, p), p, rm)
                .div(&ExactNum::from_u8(3, p), p, rm),
            p
        ));
        let x2 = ExactNum::from_u8(2, p);
        let gi2 = x2.scorer_gi(p, rm, &mut cc);
        let hi2 = x2.scorer_hi(p, rm, &mut cc);
        let bi2 = x2.bi(p, rm, &mut cc);
        assert!(near(&hi2.add(&gi2, p, rm), &bi2, p), "Hi+Gi=Bi at 2");

        let ber0 = zero.kelvin_ber(p, rm, &mut cc);
        assert!(near(&ber0, &one, p));
        let bei0 = zero.kelvin_bei(p, rm, &mut cc);
        assert!(bei0.is_zero() || near(&bei0, &zero, p));

        let cl2pi = cc.pi(p, rm).clausen_cl2(p, rm, &mut cc);
        assert!(cl2pi.is_zero() || near(&cl2pi, &zero, p), "Cl2(pi)");
        let z2 = ExactNum::from_u8(2, p).riemann_zeta(p, rm, &mut cc);
        let pi2_6 = cc
            .pi(p, rm)
            .mul(&cc.pi(p, rm), p, rm)
            .div(&ExactNum::from_u32(6, p), p, rm);
        assert!(near(&z2, &pi2_6, p), "zeta(2)");

        let z3 = ExactNum::from_u8(3, p).riemann_zeta(p, rm, &mut cc);
        let cl30 = zero.clausen_cl3(p, rm, &mut cc);
        assert!(near(&cl30, &z3, p));

        assert!(near(
            &ExactNum::from_u8(1, p).barnes_g(p, rm, &mut cc),
            &one,
            p
        ));
        assert!(near(
            &ExactNum::from_u8(2, p).barnes_g(p, rm, &mut cc),
            &one,
            p
        ));
        assert!(near(
            &ExactNum::from_u8(4, p).barnes_g(p, rm, &mut cc),
            &ExactNum::from_u8(2, p),
            p
        ));
        assert!(near(
            &ExactNum::from_u8(5, p).barnes_g(p, rm, &mut cc),
            &ExactNum::from_u32(12, p),
            p
        ));

        let psi0 = one.polygamma(0, p, rm, &mut cc);
        let g = cc.euler_gamma(p, rm).neg();
        assert!(near(&psi0, &g, p));

        let j = one.bessel_j(1, p, rm, &mut cc);
        let aj = one.anger_j(&one, p, rm, &mut cc);
        assert!(near(&aj, &j, p));

        let w0 = zero.lambert_w0(p, rm, &mut cc);
        assert!(w0.is_zero() || near(&w0, &zero, p));
        // −1/e rounded toward zero stays inside the domain (a value below −1/e is NaN).
        let inv_e = ExactNum::from_u8(1, p)
            .exp(p, RoundingMode::Up, &mut cc)
            .reciprocal(p, RoundingMode::ToZero)
            .neg();
        let wm = inv_e.lambert_w0(p, rm, &mut cc);
        assert!(near(&wm, &ExactNum::from_i8(-1, p), p));

        let e = one.exp(p, rm, &mut cc);
        let m = one.hypergeom_1f1(&half, &half, p, rm, &mut cc);
        assert!(near(&m, &e, p));

        let li2h = half.polylog(2, p, rm, &mut cc);
        // Li_2(1/2) = π²/12 − (ln 2)²/2
        let pi = cc.pi(p, rm);
        let ln2 = cc.ln_2(p, rm);
        let want = pi.sqr(p, rm).div(&ExactNum::from_u32(12, p), p, rm).sub(
            &ln2.sqr(p, rm).div(&ExactNum::from_u8(2, p), p, rm),
            p,
            rm,
        );
        assert!(near(&li2h, &want, p));

        let lna = ln_glaisher(p, &mut cc);
        let lna_ref = ExactNum::parse(
            "0.2487544770337842625472529935761139760974",
            crate::Radix::Dec,
            4 * p,
            rm,
            &mut cc,
        );
        let dl = lna.sub(&lna_ref, 2 * p, RoundingMode::None).abs();
        let lb = if dl.is_zero() {
            i32::MAX
        } else {
            i32::try_from(lna_ref.exponent().unwrap_or(0) - dl.exponent().unwrap_or(0))
                .unwrap_or(i32::MAX)
        };
        assert!(lb >= 80, "ln A only {lb} bits");

        let mpar = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(4, p), p, rm);
        let u = half.elliptic_f(&mpar, p, rm, &mut cc);
        let sn = u.jacobi_sn(&mpar, p, rm, &mut cc);
        let back = sn.jacobi_arcsn(&mpar, p, rm, &mut cc);
        assert!(near(&back, &u, p));
        Ok(())
    }

    /// Compact limits of the real Scorer functions at exact ±∞.
    /// Gi(+∞)=+0 (DLMF 9.12.25); Hi(+∞)=+∞ (9.12.29);
    /// Hi(−∞)=+0 (9.12.27); Gi(−∞)=+0 (9.12.11 + 9.12.27 + 9.7.9–11).
    /// The via-`finite` path used to return InvalidArgument at every infinity.
    #[test]
    fn zf_scorer_gi_hi_at_infinity() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let modes = [
            RoundingMode::Up,
            RoundingMode::Down,
            RoundingMode::ToZero,
            RoundingMode::FromZero,
            RoundingMode::ToEven,
            RoundingMode::ToOdd,
            RoundingMode::None,
        ];
        for p in [64usize, 256, 2048, 4096] {
            let pos = crate::INF_POS.clone();
            let neg = crate::INF_NEG.clone();
            for rm in modes {
                let gi_pos = pos.scorer_gi(p, rm, &mut cc);
                assert!(
                    gi_pos.is_zero() && !gi_pos.is_negative(),
                    "Gi(+∞) p={p} {rm:?} {:?}",
                    gi_pos
                );
                let gi_neg = neg.scorer_gi(p, rm, &mut cc);
                assert!(
                    gi_neg.is_zero() && !gi_neg.is_negative(),
                    "Gi(−∞) p={p} {rm:?} {:?}",
                    gi_neg
                );
                let hi_pos = pos.scorer_hi(p, rm, &mut cc);
                assert!(
                    hi_pos.is_inf_pos(),
                    "Hi(+∞) p={p} {rm:?} {:?}",
                    hi_pos
                );
                let hi_neg = neg.scorer_hi(p, rm, &mut cc);
                assert!(
                    hi_neg.is_zero() && !hi_neg.is_negative(),
                    "Hi(−∞) p={p} {rm:?} {:?}",
                    hi_neg
                );
            }
            // neighbors stay finite and are not the compact limit (modest p only)
            if p > 256 {
                continue;
            }
            let big = ExactNum::from_u32(100, p);
            let gi_big = big.scorer_gi(p, RoundingMode::ToEven, &mut cc);
            assert!(
                !gi_big.is_nan() && gi_big.is_positive() && !gi_big.is_zero(),
                "Gi(100) should be a tiny positive, got {:?}",
                gi_big
            );
            let hi_big = big.scorer_hi(p, RoundingMode::ToEven, &mut cc);
            assert!(
                !hi_big.is_nan() && (hi_big.is_positive() || hi_big.is_inf_pos()),
                "Hi(100) should be a huge positive or +∞, got {:?}",
                hi_big
            );
            let gi_mbig = big.clone().neg().scorer_gi(p, RoundingMode::ToEven, &mut cc);
            assert!(
                !gi_mbig.is_nan() && !gi_mbig.is_inf() && !gi_mbig.is_zero(),
                "Gi(−100) should be a finite oscillation, got {:?}",
                gi_mbig
            );
            let hi_mbig = big.clone().neg().scorer_hi(p, RoundingMode::ToEven, &mut cc);
            assert!(
                !hi_mbig.is_nan() && hi_mbig.is_positive() && !hi_mbig.is_zero(),
                "Hi(−100) should be a tiny positive, got {:?}",
                hi_mbig
            );
            let nan = ExactNum::nan(None);
            assert!(
                nan.scorer_gi(p, RoundingMode::ToEven, &mut cc).is_nan(),
                "Gi(NaN)"
            );
            assert!(
                nan.scorer_hi(p, RoundingMode::ToEven, &mut cc).is_nan(),
                "Hi(NaN)"
            );
        }
        Ok(())
    }

    /// LogG via the EM jet agrees with the Stirling formula at a height
    /// where the old path is already accurate (z=32, p=256).
    #[test]
    fn zf_barnes_logg_em_matches_stirling_at_height_32() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let p = 256usize;
        let z = ExactNum::from_u32(32, p);
        let em = crate::hurwitz_jet::log_barnes_g_zeta(&z, p, &mut cc);
        let st = barnes_ln_g_plus1(&z.sub(&ExactNum::from_u8(1, p), p, RoundingMode::None), p, &mut cc);
        assert!(near(&em, &st, p), "EM {:?} vs Stirling {:?}", em, st);
        Ok(())
    }

    /// FLINT 3.6 golds for G(±1/2) at p=64 and p=256. The p=65536 landing is
    /// the six remaining-28 TIMEOUT rows; that width is tested separately.
    #[test]
    fn zf_barnes_g_half_matches_flint_golds() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let cases: &[(usize, &str, &str, &str)] = &[
            (
                64,
                "0.6032442812094462061914292245347020799",
                "-0.1701720698965615191650119108353393216",
                "-0.5054330544896953827976849898083449517",
            ),
            (
                256,
                "0.60324428120944620619142922453470207988300342038945976538776920411943278567933529074829868836987",
                "-0.17017206989656151916501191083533932159543634434543587506565817193148051677954844820299361047756",
                "-0.50543305448969538279768498980834495172139910146661993278982756034184925772233499926828249610945",
            ),
            (
                1024,
                "0.60324428120944620619142922453470207988300342038945976538776920411943278567933529074829868836987374145499098282178654258945350301298541945816732560358510209692692633243411756714351602421054351840512029944936374185787771519482152901122412567242217387753020259455670990756936630586008503618662704901631303620538430522485783118283",
                "-0.17017206989656151916501191083533932159543634434543587506565817193148051677954844820299361047756152558220446726300626182221833720769976988706338711073742052629396032776735761700020545059602742370195421981524096266998385349195698524738205413652794130600284864016812012148785749644645775865137635398164260924315154612356632101720",
                "-0.50543305448969538279768498980834495172139910146661993278982756034184925772233499926828249610945144825215140496314472505251613924092530457949088976177276258115793711047920896359726607916443054161767942128102038575551812482701832682374028795290197380809339690939009861812749702849212581763201032051583311442539622205135106285824",
            ),
        ];
        for &(p, gpos, gneg, lng) in cases {
            let half = ExactNum::from_u8(1, p).ldexp(-1, p, RoundingMode::None);
            let g = half.barnes_g(p, RoundingMode::ToEven, &mut cc);
            let want = ExactNum::parse(gpos, crate::Radix::Dec, 4 * p, RoundingMode::ToEven, &mut cc);
            let d = ExactNum::sub_ref(&g, &want, p, RoundingMode::None).abs();
            let bits = if d.is_zero() { 10_000 } else { (want.exponent().unwrap_or(0) - d.exponent().unwrap_or(0)) as i32 };
            assert!(
                bits >= (p as i32) - 16 || d.is_zero(),
                "G(1/2) p={p} only {bits} bits"
            );
            let gm = half.clone().neg().barnes_g(p, RoundingMode::ToEven, &mut cc);
            let wantn = ExactNum::parse(gneg, crate::Radix::Dec, 4 * p, RoundingMode::ToEven, &mut cc);
            assert!(near(&gm, &wantn, p), "G(-1/2) p={p} got {:?}", gm);
            let lg = crate::hurwitz_jet::log_barnes_g_zeta(&half, p, &mut cc);
            let wantl = ExactNum::parse(lng, crate::Radix::Dec, 4 * p, RoundingMode::ToEven, &mut cc);
            assert!(near(&lg, &wantl, p), "LogG(1/2) p={p} got {:?}", lg);
        }
        Ok(())
    }

    /// FLINT 3.6 golds for G(±1/2) at p=4096 and 16384 (dest widths below
    /// the six remaining-28 TIMEOUT rows at p=65536, which live in the
    /// ignored `zf_barnes_g_half_matches_flint_golds_p65536` test).
    #[test]
    fn zf_barnes_g_half_matches_flint_golds_high_p() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let raw = include_str!("testdata/barnes_g_half_flint_golds.txt");
        for ln in raw.lines() {
            let ln = ln.trim();
            if ln.is_empty() || ln.starts_with('#') {
                continue;
            }
            let mut parts = ln.splitn(3, ' ');
            let p: usize = parts.next().unwrap().parse().unwrap();
            let sign = parts.next().unwrap();
            let dec = parts.next().unwrap();
            if p >= 65_536 {
                continue;
            }
            let t0 = std::time::Instant::now();
            let mut z = ExactNum::from_u8(1, p).ldexp(-1, p, RoundingMode::None);
            if sign == "-" {
                z = z.neg();
            }
            let g = z.barnes_g(p, RoundingMode::ToEven, &mut cc);
            let ms = t0.elapsed().as_millis();
            let want = ExactNum::parse(dec, crate::Radix::Dec, 4 * p, RoundingMode::ToEven, &mut cc);
            let d = ExactNum::sub_ref(&g, &want, p, RoundingMode::None).abs();
            let bits = if d.is_zero() {
                10_000
            } else {
                (want.exponent().unwrap_or(0) - d.exponent().unwrap_or(0)) as i32
            };
            eprintln!("G({sign}1/2) p={p}: {ms} ms, {bits} bits");
            assert!(
                bits >= (p as i32) - 16 || d.is_zero(),
                "G({sign}1/2) p={p} only {bits} bits after {ms} ms"
            );
            assert!(
                ms < 90_000,
                "G({sign}1/2) p={p} must finish under 90 s; took {ms} ms"
            );
        }
        Ok(())
    }


    /// FLINT golds for G(±1/2) at dest p=2^16. Ignored in the default suite:
    /// the sieved EM jet is correct (~330 s, full width) but over `ZFS_WALL`.
    #[test]
    #[ignore = "G(±1/2) at p=65536 is ~330 s; over ZFS_WALL because of ln cost"]
    fn zf_barnes_g_half_matches_flint_golds_p65536() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let raw = include_str!("testdata/barnes_g_half_flint_golds.txt");
        for ln in raw.lines() {
            let ln = ln.trim();
            if ln.is_empty() || ln.starts_with('#') {
                continue;
            }
            let mut parts = ln.splitn(3, ' ');
            let p: usize = parts.next().unwrap().parse().unwrap();
            let sign = parts.next().unwrap();
            let dec = parts.next().unwrap();
            if p < 65_536 {
                continue;
            }
            let t0 = std::time::Instant::now();
            let mut z = ExactNum::from_u8(1, p).ldexp(-1, p, RoundingMode::None);
            if sign == "-" {
                z = z.neg();
            }
            let g = z.barnes_g(p, RoundingMode::ToEven, &mut cc);
            let ms = t0.elapsed().as_millis();
            let want = ExactNum::parse(dec, crate::Radix::Dec, 4 * p, RoundingMode::ToEven, &mut cc);
            let d = ExactNum::sub_ref(&g, &want, p, RoundingMode::None).abs();
            let bits = if d.is_zero() {
                10_000
            } else {
                (want.exponent().unwrap_or(0) - d.exponent().unwrap_or(0)) as i32
            };
            eprintln!("G({sign}1/2) p={p}: {ms} ms, {bits} bits");
            assert!(
                bits >= (p as i32) - 16 || d.is_zero(),
                "G({sign}1/2) p={p} only {bits} bits after {ms} ms"
            );
        }
        Ok(())
    }

    /// G(2^{30}) overflows the exponent range (DLMF 5.17.5: log₂ G(z) ∼ (z²/2) log₂ z).
    /// C99 7.12.1: Down/ToZero must return +MAX, not +Inf; Up returns +Inf.
    /// G(+Inf) = +Inf in every mode; G(+MAX) overflows by the same rule (not InvalidArgument).
    #[test]
    fn zf_barnes_g_large_positive_overflows_by_mode() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let p = 128;
        let two_30 = ExactNum::from_u32(1u32 << 30, p);
        let down = two_30.barnes_g(p, RoundingMode::Down, &mut cc);
        let tozero = two_30.barnes_g(p, RoundingMode::ToZero, &mut cc);
        let up = two_30.barnes_g(p, RoundingMode::Up, &mut cc);
        let max = ExactNum::max_value(p);
        // 3·2^{28} ∈ [2^{29}, 2^{30}): the lb≥30 bound does not fire, but
        // (z²/2) log₂ z still exceeds EMAX.
        let mid = ExactNum::from_u32(3u32 << 28, p);
        let mid_down = mid.barnes_g(p, RoundingMode::Down, &mut cc);
        let mid_up = mid.barnes_g(p, RoundingMode::Up, &mut cc);
        assert!(
            finite(&mid_down) && !mid_down.is_negative() && mid_down.compare(&max) == Some(0),
            "G(3·2^28) Down should be +MAX, got inf={} nan={}",
            mid_down.is_inf(),
            mid_down.is_nan()
        );
        assert!(
            mid_up.is_inf_pos(),
            "G(3·2^28) Up should be +Inf, got nan={}",
            mid_up.is_nan()
        );
        assert!(
            finite(&down) && !down.is_negative() && down.compare(&max) == Some(0),
            "G(2^30) Down should be +MAX, got inf={} nan={}",
            down.is_inf(),
            down.is_nan()
        );
        assert!(
            finite(&tozero) && !tozero.is_negative() && tozero.compare(&max) == Some(0),
            "G(2^30) ToZero should be +MAX, got inf={} nan={}",
            tozero.is_inf(),
            tozero.is_nan()
        );
        assert!(
            up.is_inf_pos(),
            "G(2^30) Up should be +Inf, got nan={}",
            up.is_nan()
        );

        let inf = crate::ext::INF_POS;
        for rm in [RoundingMode::Down, RoundingMode::ToZero, RoundingMode::Up, RoundingMode::ToEven]
        {
            let g = inf.barnes_g(p, rm, &mut cc);
            assert!(
                g.is_inf_pos(),
                "G(+Inf) should be +Inf in {rm:?}, nan={}",
                g.is_nan()
            );
        }

        let gmax_down = max.barnes_g(p, RoundingMode::Down, &mut cc);
        let gmax_up = max.barnes_g(p, RoundingMode::Up, &mut cc);
        assert!(
            finite(&gmax_down) && !gmax_down.is_negative() && gmax_down.compare(&max) == Some(0),
            "G(+MAX) Down should be +MAX, got inf={} nan={}",
            gmax_down.is_inf(),
            gmax_down.is_nan()
        );
        assert!(
            gmax_up.is_inf_pos(),
            "G(+MAX) Up should be +Inf, got nan={}",
            gmax_up.is_nan()
        );
        Ok(())
    }

    /// 1F0(1;;2) = (1-2)^{-1} = -1 exactly (DLMF 16.2.5). The Meijer continuation
    /// at f7cdde7/82bd851 missed the last ulp in Up (d0008803 and the four
    /// higher-precision Up twins).
    #[test]
    fn zf_hypergeom_1f0_at_two_is_exact_minus_one() -> Result<(), crate::Error> {
        let p = 64;
        let mut cc = Consts::new()?;
        let z = ExactNum::from_u8(2, p);
        let a = ExactNum::from_u8(1, p);
        let want = ExactNum::from_i8(-1, p);
        for rm in [RoundingMode::Down, RoundingMode::ToZero, RoundingMode::Up, RoundingMode::ToEven]
        {
            let g = z.hypergeom_pfq(core::slice::from_ref(&a), &[], p, rm, &mut cc);
            assert!(
                finite(&g) && g.compare(&want) == Some(0),
                "1F0(1;;2) in {rm:?} should be exactly -1, got inf={} nan={} exp={:?}",
                g.is_inf(),
                g.is_nan(),
                g.exponent()
            );
        }
        Ok(())
    }

    /// Li_2(1 − 2^{-100}) is within a few ulps of ζ(2) (DLMF 25.12.2). The defining
    /// series at 82bd851/291594b only recovered 7 bits (t-zf-156-polylog).
    #[test]
    fn zf_polylog_dilog_near_one() -> Result<(), crate::Error> {
        let p = 256;
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;
        let one = ExactNum::from_u8(1, p);
        let eps = ExactNum::from_u8(1, p).ldexp(-100, p, RoundingMode::None);
        let z = one.sub(&eps, p, RoundingMode::None);
        let got = z.polylog(2, p, rm, &mut cc);
        // Gold from the same identity, with Li_2(2^{-100}) from the series (tiny z).
        let wrk = work_p(p).saturating_add(WORD_BIT_SIZE);
        let nr = RoundingMode::None;
        let li2_eps = eps.polylog(2, wrk, nr, &mut cc);
        let ln_z = z.ln(wrk, nr, &mut cc);
        let ln_eps = eps.ln(wrk, nr, &mut cc);
        let zeta2 = ExactNum::from_u8(2, wrk).riemann_zeta(wrk, nr, &mut cc);
        let want = zeta2
            .sub(&li2_eps, wrk, nr)
            .sub(&ln_z.mul(&ln_eps, wrk, nr), wrk, nr);
        let want = finish(want, p, rm);
        assert!(
            finite(&got) && got.compare(&want) == Some(0),
            "Li_2(1-2^{{-100}}) ToEven should match DLMF 25.12.2, got inf={} nan={}",
            got.is_inf(),
            got.is_nan()
        );
        Ok(())
    }

    /// The former `n > 62` refusal (and `hurwitz_zeta`'s `small_int(..., 2, 64)`
    /// range that made `n = 64` impossible) is gone. `n = 63` uses
    /// `ζ(64, z)`; `n = 64` needs `ζ(65, z)`.
    ///
    /// Golds: mpmath `polygamma` / `zeta` at 260 dps (no limit), agreeing with
    /// FLINT `acb` (`z.polygamma(63)` at z=1, `acb.zeta(65, 3/2)`).
    #[test]
    fn caps_polygamma_order_above_62() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let modes = [
            RoundingMode::ToEven,
            RoundingMode::ToOdd,
            RoundingMode::Up,
            RoundingMode::Down,
            RoundingMode::ToZero,
            RoundingMode::FromZero,
            RoundingMode::None,
        ];
        fn dec(s: &str, wp: usize, cc: &mut Consts) -> ExactNum {
            ExactNum::parse(s, crate::Radix::Dec, wp, RoundingMode::ToEven, cc)
        }
        let half3 = |wp: usize, rm: RoundingMode| {
            ExactNum::from_u8(3, wp).div(&ExactNum::from_u8(2, wp), wp, rm)
        };
        // (kind, n, z-tag, gold)
        let rows: &[(&str, u64, &str, &str)] = &[
            (
                "pg",
                63,
                "1",
                "1.982608315404440064223624120499941542539576696683473298878843684608867521539038300839454224920286282265801298832647262360747899047118665117470397350554297883462598454622584527035715463837940718445621638691808802548235194e+87",
            ),
            (
                "pg",
                63,
                "3/2",
                "1.065114614704354774083764694179462660512011471463142372110324241059047586063749768088287670860477432803596786010891101003768290121765596839008910161155155849918185055881767138602262546860427719705762614570805628259383470e+76",
            ),
            (
                "pg",
                64,
                "3/2",
                "-4.544489022738568855454557905100303452821656239017902378244380322112111447056476691735435395961717998511182663752882031613632858590970875789831817123311242573234578925986576826574896117772306508619996501460057888470552830e+77",
            ),
            (
                "hz",
                65,
                "3/2",
                "3.581526438105602847926440864063557031949546716502448624191871882983324210396937084017395540706111037533681084763596594639820758387295454430614404409367941426177117392975198825398343054876253866376732416947990863959777414e-12",
            ),
        ];
        let mut bad = Vec::new();
        for &(kind, n, ztag, gold_s) in rows {
            for &p in &[64usize, 256] {
                for &rm in &modes {
                    let wp = 1024;
                    let want = finish(dec(gold_s, wp, &mut cc), p, rm);
                    let z = if ztag == "1" {
                        ExactNum::from_u8(1, wp)
                    } else {
                        half3(wp, RoundingMode::None)
                    };
                    let got = if kind == "pg" {
                        z.polygamma(n as usize, p, rm, &mut cc)
                    } else {
                        ExactNum::from_u64(n, wp).hurwitz_zeta(&z, p, rm, &mut cc)
                    };
                    if got.compare(&want) != Some(0) {
                        bad.push(format!(
                            "{kind}({n}, {ztag}) p={p} {rm:?} got {got} want {want}"
                        ));
                    }
                }
            }
        }
        assert!(bad.is_empty(), "{} wrong:\n{}", bad.len(), bad.join("\n"));
        Ok(())
    }

    /// arcdn(0|1) = K(1) = +∞, because dn(u|1) = sech(u) → 0 as u → +∞.
    /// Both signs of zero; every rounding mode returns +Inf (a pole, not a
    /// finite overflow that C99 7.12.1 would map to +MAX).
    #[test]
    fn zf_jacobi_arcdn_zero_at_m_one_is_pos_inf() {
        let mut cc = Consts::new().expect("Consts");
        let p = 64;
        let modes = [
            RoundingMode::ToEven,
            RoundingMode::ToOdd,
            RoundingMode::Up,
            RoundingMode::Down,
            RoundingMode::ToZero,
            RoundingMode::FromZero,
            RoundingMode::None,
        ];
        let one = ExactNum::from_u8(1, p);
        let pos0 = ExactNum::new(p);
        let mut neg0 = ExactNum::new(p);
        neg0.set_sign(crate::Sign::Neg);
        for (tag, x) in [("+0", &pos0), ("-0", &neg0)] {
            for &rm in &modes {
                let g = x.jacobi_arcdn(&one, p, rm, &mut cc);
                assert!(
                    g.is_inf_pos(),
                    "arcdn({tag}|1) p={p} {rm:?} should be +Inf, got inf={} nan={} err={:?}",
                    g.is_inf(),
                    g.is_nan(),
                    g.err()
                );
            }
        }
        let p2 = 2048;
        let one2 = ExactNum::from_u8(1, p2);
        let z2 = ExactNum::new(p2);
        let g2 = z2.jacobi_arcdn(&one2, p2, RoundingMode::ToEven, &mut cc);
        assert!(
            g2.is_inf_pos(),
            "arcdn(+0|1) p=2048 ToEven should be +Inf, got inf={} nan={} err={:?}",
            g2.is_inf(),
            g2.is_nan(),
            g2.err()
        );
    }
}
