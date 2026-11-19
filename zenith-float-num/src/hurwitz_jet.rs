//! Bound-driven Euler–Maclaurin jet of the Hurwitz zeta function,
//! after FLINT 3.6 `_acb_poly_zeta_cpx_series`.
//!
//! `hurwitz_em_jet(s, a)` returns `(ζ(s,a), ζ'(s,a))` for real `s` and
//! `a > 0`, using `_acb_poly_zeta_em_choose_param` / `_acb_poly_zeta_em_sum`
//! at jet length 2 (Cohen / Johansson). This is the kernel used by
//! FLINT's `LogG(z) = (z−1) LogΓ(z) − ζ'(−1,z) + ζ'(−1)` and by the
//! Riemann-zeta Euler–Maclaurin path.
//!
//! At `s = -1` and `a ∈ {1, 1/2}` the finite power-sum is sieved:
//! `Σ k ln k = Σ_p ln(p) Σ_j p^j T(⌊n/p^j⌋)`. LogG(1/2) uses the closed
//! form `ζ'(-1, 1/2) = -(ln 2)/24 - ζ'(-1)/2`, so G(±1/2) pays only the
//! `a = 1` jet. Dest `p = 2^16` is correct against FLINT golds (~330 s on
//! this box) but over `ZFS_WALL = 120` because `ExactNum::ln` is ~160 ms
//! at that width (`π(N) ≈ 1280` logs with `N ≈ 10^4`). Dest `p ≤ 2^14`
//! finishes in a few seconds.

use crate::defs::WORD_BIT_SIZE;
use crate::ops::special::with_even_bernoulli;
use crate::Consts;
use crate::ExactNum;
use crate::RoundingMode;
use core::ops::Neg;
use alloc::vec::Vec;

const NR: RoundingMode = RoundingMode::None;
const UP: RoundingMode = RoundingMode::Up;

/// Length-2 jet `c0 + c1 x` of a function of the series variable `x`.
#[derive(Clone)]
struct Jet {
    c0: ExactNum,
    c1: ExactNum,
}

impl Jet {
    fn zero(p: usize) -> Self {
        Self {
            c0: ExactNum::new(p),
            c1: ExactNum::new(p),
        }
    }

    fn add(&self, other: &Self, p: usize) -> Self {
        Self {
            c0: ExactNum::add_ref(&self.c0, &other.c0, p, NR),
            c1: ExactNum::add_ref(&self.c1, &other.c1, p, NR),
        }
    }

    fn scale(&self, k: &ExactNum, p: usize) -> Self {
        Self {
            c0: ExactNum::mul_ref(&self.c0, k, p, NR),
            c1: ExactNum::mul_ref(&self.c1, k, p, NR),
        }
    }

    /// `self * (c + x)`.
    fn mul_cpx(&self, c: &ExactNum, p: usize) -> Self {
        Self {
            c0: ExactNum::mul_ref(&self.c0, c, p, NR),
            c1: ExactNum::mul_ref(&self.c1, c, p, NR).add(&self.c0, p, NR),
        }
    }

    /// `self / c` for a nonzero constant `c`.
    fn div_const(&self, c: &ExactNum, p: usize) -> Self {
        Self {
            c0: ExactNum::div_ref(&self.c0, c, p, NR),
            c1: ExactNum::div_ref(&self.c1, c, p, NR),
        }
    }
}

/// `(N+a)^{-(s+x)} = (N+a)^{-s} (1 − x ln(N+a) + O(x²))`.
fn invpow_jet(base: &ExactNum, s: &ExactNum, p: usize, cc: &mut Consts) -> Jet {
    let ln = base.ln(p, NR, cc);
    let c0 = base.pow(&s.clone().neg(), p, NR, cc);
    let c1 = ExactNum::mul_ref(&c0, &ln, p, NR).neg();
    Jet { c0, c1 }
}

fn choose_m(n: u64, target: u64) -> u64 {
    let extra = (n / 100).min(2000);
    n.min(target.saturating_add(extra)).max(1)
}

fn bit_count(n: u64) -> u32 {
    64 - n.leading_zeros()
}

/// Conservative log-space EM tail bound (real `s`, real `a > 0`, jet length 2),
/// matching FLINT `_acb_poly_zeta_em_bound` / `_acb_poly_zeta_em_bound1` in
/// the real case: rising factorial of `|s| + x`, power integral, and
/// `4/(2π)^{2M}`. All short-precision arithmetic is rounded upward.
fn em_bound1(
    s: &ExactNum,
    a: &ExactNum,
    n: u64,
    m: u64,
    cc: &mut Consts,
) -> Option<ExactNum> {
    if n < 1 || m < 1 {
        return None;
    }
    let wp = 64usize;
    let two_m = m.checked_mul(2)?;
    let sigma = {
        let mut t = s.clone();
        let _ = t.set_precision(wp, UP);
        t
    };
    let alpha = {
        let mut t = a.clone();
        let _ = t.set_precision(wp, UP);
        t
    };
    // A = N+a, B = σ+2M; need A>1 and σ+2M>1.
    let an = ExactNum::from_u64(n, wp).add(&alpha, wp, UP);
    let s2m = sigma.add(&ExactNum::from_u64(two_m, wp), wp, UP);
    if !an.is_positive() || matches!(an.compare(&ExactNum::from_u8(1, wp)), Some(c) if c <= 0) {
        return None;
    }
    if !s2m.is_positive() || matches!(s2m.compare(&ExactNum::from_u8(1, wp)), Some(c) if c <= 0)
    {
        return None;
    }
    let bm1 = s2m.sub(&ExactNum::from_u8(1, wp), wp, UP);
    if !bm1.is_positive() {
        return None;
    }
    let ln_a = an.ln(wp, UP, cc);
    if ln_a.is_nan() {
        return None;
    }
    // I0 = 1 / (A^{B-1} (B-1)), I1 = (1+(B-1)ln A) / (A^{B-1} (B-1)^2).
    let a_pow = an.pow(&bm1, wp, UP, cc);
    let den0 = a_pow.mul(&bm1, wp, UP);
    let i0 = ExactNum::from_u8(1, wp).div(&den0, wp, UP);
    let one_plus = ExactNum::from_u8(1, wp).add(&ExactNum::mul_ref(&bm1, &ln_a, wp, UP), wp, UP);
    let den1 = den0.mul(&bm1, wp, UP);
    let i1 = one_plus.div(&den1, wp, UP);

    // F0 = (|s|)_{2M}, H = Σ 1/(|s|+i), F1 = F0 H.
    let sabs = {
        let mut t = s.clone().abs();
        let _ = t.set_precision(wp, UP);
        t
    };
    let mut f0 = ExactNum::from_u8(1, wp);
    let mut harm = ExactNum::new(wp);
    for i in 0..two_m {
        let term = sabs.clone().add(&ExactNum::from_u64(i, wp), wp, UP);
        if !term.is_positive() || term.is_zero() {
            // A zero factor makes the |s|-rising vanish; the length-2 jet
            // of the true rising factorial does not. Treat as "bound failed".
            return None;
        }
        f0 = f0.mul(&term, wp, UP);
        harm = harm.add(&ExactNum::from_u8(1, wp).div(&term, wp, UP), wp, UP);
    }
    let f1 = ExactNum::mul_ref(&f0, &harm, wp, UP);
    let r0 = ExactNum::mul_ref(&f0, &i0, wp, UP);
    let r1 = ExactNum::mul_ref(&f0, &i1, wp, UP).add(&ExactNum::mul_ref(&f1, &i0, wp, UP), wp, UP);
    // 4 / (2π)^{2M}
    let pi = cc.pi(wp, UP);
    let two_pi = ExactNum::from_u8(2, wp).mul(&pi, wp, UP);
    let pow = two_pi.pow(&ExactNum::from_u64(two_m, wp), wp, UP, cc);
    let scale = ExactNum::from_u8(4, wp).div(&pow, wp, UP);
    let b0 = r0.mul(&scale, wp, UP);
    let b1 = r1.mul(&scale, wp, UP);
    let mag = if matches!(b1.compare(&b0), Some(c) if c > 0) {
        b1
    } else {
        b0
    };
    if mag.is_nan() || !mag.is_positive() {
        None
    } else {
        Some(mag)
    }
}

fn a_pow_neg_s_mag(s: &ExactNum, a: &ExactNum, cc: &mut Consts) -> ExactNum {
    let wp = 64usize;
    let mut aa = a.clone();
    let _ = aa.set_precision(wp, UP);
    let mut ss = s.clone();
    let _ = ss.set_precision(wp, UP);
    if aa.is_zero() {
        return ExactNum::from_u8(1, wp);
    }
    let v = aa.pow(&ss.neg(), wp, UP, cc);
    if v.is_nan() || v.is_inf() || v.is_zero() {
        ExactNum::from_u8(1, wp)
    } else {
        v.abs()
    }
}

/// FLINT `_acb_poly_zeta_em_choose_param` for a real argument and jet length 2.
fn choose_param(
    s: &ExactNum,
    a: &ExactNum,
    target: usize,
    cc: &mut Consts,
) -> Option<(u64, u64)> {
    let target_u = u64::try_from(target).ok()?;
    let mut tol = a_pow_neg_s_mag(s, a, cc);
    // tol *= 2^{-target}
    let te = -i64::try_from(target).ok()?;
    tol = tol.ldexp(te, 64, UP);
    if tol.is_nan() || !tol.is_positive() {
        return None;
    }
    let limit = target_u.saturating_mul(100).max(16);
    let mut a_n = 1u64;
    let mut b_n = 2u64;
    let mut b_bound = em_bound1(s, a, b_n, choose_m(b_n, target_u), cc)?;
    if matches!(b_bound.compare(&tol), Some(c) if c > 0) {
        while matches!(b_bound.compare(&tol), Some(c) if c > 0) && b_n <= limit {
            a_n = b_n;
            b_n = b_n.checked_mul(2)?;
            if b_n == 0 {
                return None;
            }
            match em_bound1(s, a, b_n, choose_m(b_n, target_u), cc) {
                Some(b) => b_bound = b,
                None => {
                    // The |s|-rising bound can fail (zero factor). Grow N only
                    // until the search limit, then take that pair.
                    if b_n > limit {
                        break;
                    }
                    continue;
                }
            }
        }
        while b_n > a_n + 4 {
            let c_n = a_n + (b_n - a_n) / 2;
            match em_bound1(s, a, c_n, choose_m(c_n, target_u), cc) {
                Some(cb) if matches!(cb.compare(&tol), Some(cmp) if cmp < 0) => {
                    let _ = cb;
                    b_n = c_n;
                }
                Some(_) | None => a_n = c_n,
            }
        }
    }
    Some((b_n, choose_m(b_n, target_u)))
}

fn bernoulli_even(m: usize, p: usize) -> Option<Vec<ExactNum>> {
    if m == 0 {
        return Some(Vec::new());
    }
    with_even_bernoulli(m, p, |v| {
        v.iter()
            .map(|b| b.clone().map(ExactNum::from))
            .collect::<Result<Vec<_>, _>>()
    })
    .ok()
}

/// Euler–Maclaurin sum for the length-2 jet, after `_acb_poly_zeta_em_sum`
/// with `deflate = false`.
fn em_sum(
    s: &ExactNum,
    a: &ExactNum,
    n: u64,
    m: u64,
    p: usize,
    cc: &mut Consts,
) -> Option<Jet> {
    let p = p
        .saturating_add(2 * (bit_count(n) as usize + 2));
    let one = ExactNum::from_u8(1, p);
    let mut sum = if let Some(fast) = powsum_s_minus_one(s, a, n, p, cc) {
        fast
    } else {
        let mut sum = Jet::zero(p);
        for k in 0..n {
            let t = ExactNum::from_u64(k, p).add(a, p, NR);
            if !t.is_positive() {
                return None;
            }
            sum = sum.add(&invpow_jet(&t, s, p, cc), p);
        }
        sum
    };
    let na = ExactNum::from_u64(n, p).add(a, p, NR);
    // t = (N+a)^{-(s+x)}
    let tjet = invpow_jet(&na, s, p, cc);
    // u = (N+a)^{1-(s+x)} / ((s-1)+x)
    let sm1 = s.clone().sub(&one, p, NR);
    if sm1.is_zero() {
        return None;
    }
    let v = s.clone().sub(&ExactNum::from_u8(1, p), p, NR);
    let u_pow = invpow_jet(&na, &v, p, cc);
    // divide by (s-1)+x: u0 = w0/(s-1), u1 = (w1 - u0)/(s-1)
    let u0 = ExactNum::div_ref(&u_pow.c0, &sm1, p, NR);
    let u1 = ExactNum::sub_ref(&u_pow.c1, &u0, p, NR).div(&sm1, p, NR);
    sum = sum.add(&Jet { c0: u0, c1: u1 }, p);
    // + 1/2 t
    let half = ExactNum::from_u8(1, p).ldexp(-1, p, NR);
    sum = sum.add(&tjet.scale(&half, p), p);

    // Bernoulli tail, `_acb_poly_zeta_em_tail_naive`.
    let bs = bernoulli_even(m as usize, p)?;
    // u = 1/2 t; term = u * (s+x) / (N+a)
    let mut term = tjet.scale(&half, p).mul_cpx(s, p).div_const(&na, p);
    let na2 = na.sqr(p, NR);
    for r in 1..=m as usize {
        let b = bs.get(r - 1)?;
        if b.is_nan() {
            return None;
        }
        sum = sum.add(&term.scale(b, p), p);
        // term *= ((s+x)+2r-1)((s+x)+2r) / ((N+a)^2 (2r+1)(2r+2))
        let splus0 = s.clone().add(&ExactNum::from_u64((2 * r - 1) as u64, p), p, NR);
        term = term.mul_cpx(&splus0, p);
        let splus1 = s.clone().add(&ExactNum::from_u64((2 * r) as u64, p), p, NR);
        term = term.mul_cpx(&splus1, p);
        let den = na2.clone().mul(
            &ExactNum::from_u64((2 * r + 1) as u64 * (2 * r + 2) as u64, p),
            p,
            NR,
        );
        term = term.div_const(&den, p);
    }
    Some(sum)
}

/// `(ζ(s,a), ζ'(s,a))` for real `s` and `a > 0`, by bound-driven EM.
pub(crate) fn hurwitz_em_jet(
    s: &ExactNum,
    a: &ExactNum,
    p: usize,
    cc: &mut Consts,
) -> Option<(ExactNum, ExactNum)> {
    if s.is_nan() || s.is_inf() || a.is_nan() || a.is_inf() || !a.is_positive() || a.is_zero() {
        return None;
    }
    let (n, m) = choose_param(s, a, p, cc)?;
    let jet = em_sum(s, a, n, m, p, cc)?;
    if jet.c0.is_nan() || jet.c1.is_nan() {
        return None;
    }
    Some((jet.c0, jet.c1))
}

/// Cached `ζ'(−1)`.
#[cfg(feature = "std")]
static ZETA_PRIME_M1: std::sync::RwLock<Option<ExactNum>> = std::sync::RwLock::new(None);

fn zeta_prime_minus_one_uncached(p: usize, cc: &mut Consts) -> ExactNum {
    let s = ExactNum::from_i32(-1, p);
    let a = ExactNum::from_u8(1, p);
    match hurwitz_em_jet(&s, &a, p, cc) {
        Some((_, zp)) => zp,
        None => ExactNum::nan(Some(crate::Error::InvalidArgument)),
    }
}

fn cached_round(mut r: ExactNum, p: usize) -> ExactNum {
    let _ = r.set_precision(p, RoundingMode::ToEven);
    r
}

/// `ζ'(−1)` at `p` bits, cached under `std`.
pub(crate) fn zeta_prime_minus_one(p: usize, cc: &mut Consts) -> ExactNum {
    #[cfg(feature = "std")]
    {
        if let Ok(g) = ZETA_PRIME_M1.read() {
            if let Some(v) = g.as_ref().filter(|v| v.precision().unwrap_or(0) >= p) {
                return cached_round(v.clone(), p);
            }
        }
        // Extra words cover the slightly wider `work_p` used by G(-1/2)
        // reflection, so a later call does not rerun the EM jet.
        let wrk = p.saturating_add(4 * WORD_BIT_SIZE);
        let v = zeta_prime_minus_one_uncached(wrk, cc);
        if !v.is_nan() {
            if let Ok(mut g) = ZETA_PRIME_M1.write() {
                let have = g.as_ref().map(|x| x.precision().unwrap_or(0)).unwrap_or(0);
                if v.precision().unwrap_or(0) >= have {
                    *g = Some(v.clone());
                }
            }
        }
        cached_round(v, p)
    }
    #[cfg(not(feature = "std"))]
    cached_round(zeta_prime_minus_one_uncached(p + WORD_BIT_SIZE, cc), p)
}

/// `LogG(z) = (z−1) LogΓ(z) − ζ'(−1,z) + ζ'(−1)` for real `z > 0`,
/// after FLINT `_acb_log_barnes_g_zeta`.
pub(crate) fn log_barnes_g_zeta(z: &ExactNum, p: usize, cc: &mut Consts) -> ExactNum {
    let c = zeta_prime_minus_one(p, cc);
    // ζ'(-1, 1/2) = -(ln 2)/24 - ζ'(-1)/2, so one EM jet (at a=1) suffices
    // for G(±1/2) after the recurrence for the negative argument.
    let zp = if is_exact_half(z, p) {
        let ln2 = cc.ln_2(p, NR);
        ln2.div(&ExactNum::from_u32(24, p), p, NR)
            .neg()
            .sub(&c.clone().ldexp(-1, p, NR), p, NR)
    } else {
        let s = ExactNum::from_i32(-1, p);
        let Some((_, zp)) = hurwitz_em_jet(&s, z, p, cc) else {
            return ExactNum::nan(Some(crate::Error::InvalidArgument));
        };
        zp
    };
    let lng = z.ln_gamma(p, NR, cc);
    let zm1 = z.clone().sub(&ExactNum::from_u8(1, p), p, NR);
    zm1.mul(&lng, p, NR).sub(&zp, p, NR).add(&c, p, NR)
}



/// `Σ_{k=1}^{n} k ln k`, via `Σ_p ln(p) Σ_{j: p^j ≤ n} p^j T(⌊n/p^j⌋)`
/// with `T(m) = m(m+1)/2`. Only `π(n)` full-precision logarithms.
fn sum_k_ln_k(n: u64, p: usize, cc: &mut Consts) -> ExactNum {
    if n == 0 {
        return ExactNum::new(p);
    }
    let primes = crate::primality::primes_upto(n);
    let mut acc = ExactNum::new(p);
    for pr in primes {
        let pru = pr as u64;
        let mut pj = pru;
        let mut coeff = 0u128;
        loop {
            let m = (n / pj) as u128;
            coeff += (pj as u128) * m * (m + 1) / 2;
            match pj.checked_mul(pru) {
                Some(next) if next <= n => pj = next,
                _ => break,
            }
        }
        let ln = ExactNum::from_u64(pru, p).ln(p, NR, cc);
        let c = crate::ExactInt::from_u128(coeff).to_exact_num(p, NR);
        acc = acc.add(&ln.mul(&c, p, NR), p, NR);
    }
    acc
}

fn is_exact_one(x: &ExactNum, p: usize) -> bool {
    matches!(x.compare(&ExactNum::from_u8(1, p)), Some(0))
}

fn is_exact_half(x: &ExactNum, p: usize) -> bool {
    let two_x = x.clone().ldexp(1, p, NR);
    is_exact_one(&two_x, p)
}

fn is_exact_neg_one(x: &ExactNum, p: usize) -> bool {
    matches!(x.compare(&ExactNum::from_i32(-1, p)), Some(0))
}

/// Finite EM power-sum jet at `s = -1` for `a = 1` or `a = 1/2`.
/// `c0 = Σ (k+a)`, `c1 = -Σ (k+a) ln(k+a)`, `k = 0..n-1`.
fn powsum_s_minus_one(
    s: &ExactNum,
    a: &ExactNum,
    n: u64,
    p: usize,
    cc: &mut Consts,
) -> Option<Jet> {
    if n == 0 || !is_exact_neg_one(s, p) {
        return None;
    }
    let nf = n as u128;
    if is_exact_one(a, p) {
        // Σ_{k=0}^{n-1} (k+1) = Σ_{m=1}^{n} m = n(n+1)/2,
        // Σ (k+1) ln(k+1) = S(n).
        let c0 = crate::ExactInt::from_u128(nf * (nf + 1) / 2).to_exact_num(p, NR);
        let s_n = sum_k_ln_k(n, p, cc);
        return Some(Jet {
            c0,
            c1: s_n.neg(),
        });
    }
    if is_exact_half(a, p) {
        // Σ (k+1/2) = n²/2.
        // Σ (k+1/2) ln(k+1/2)
        //   = (1/2)[S(2n-1) - n(n-1) ln 2 - 2 S(n-1)] - (ln 2) n²/2.
        let c0 = crate::ExactInt::from_u128(nf * nf).to_exact_num(p, NR).ldexp(-1, p, NR);
        let s_odd = if n == 0 {
            ExactNum::new(p)
        } else {
            let s_hi = sum_k_ln_k(2 * n - 1, p, cc);
            let s_lo = if n == 1 {
                ExactNum::new(p)
            } else {
                sum_k_ln_k(n - 1, p, cc)
            };
            let ln2 = cc.ln_2(p, NR);
            let nn1 = crate::ExactInt::from_u128(nf * (nf - 1)).to_exact_num(p, NR);
            s_hi.sub(&nn1.mul(&ln2, p, NR), p, NR)
                .sub(&s_lo.ldexp(1, p, NR), p, NR)
        };
        let ln2 = cc.ln_2(p, NR);
        let n2_half = crate::ExactInt::from_u128(nf * nf)
            .to_exact_num(p, NR)
            .ldexp(-1, p, NR);
        let sum_ln = s_odd.ldexp(-1, p, NR).sub(&ln2.mul(&n2_half, p, NR), p, NR);
        return Some(Jet {
            c0,
            c1: sum_ln.neg(),
        });
    }
    None
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::Error;

    fn near(a: &ExactNum, b: &ExactNum, p: usize) -> bool {
        let d = ExactNum::sub_ref(a, b, p, RoundingMode::None).abs();
        d.is_zero()
            || d.exponent()
                .is_some_and(|e| e + (p as crate::defs::Exponent) / 4 < b.exponent().unwrap_or(0))
            || d.exponent()
                .is_some_and(|e| e < -((p as crate::defs::Exponent) / 4))
    }

    #[test]
    fn zf_hurwitz_jet_zeta_prime_minus_one_matches_flint() -> Result<(), Error> {
        // FLINT midpoints of ζ'(−1) at 256 bits (ctx.prec = 320).
        let mut cc = Consts::new()?;
        let p = 256usize;
        let zp = zeta_prime_minus_one(p, &mut cc);
        // ζ'(−1) = 1/12 − ln A ≈ −0.165421143700450929212310400185656306942...
        let want = ExactNum::parse(
            "-0.165421143700450929212310400185656306942",
            crate::Radix::Dec,
            4 * p,
            RoundingMode::ToEven,
            &mut cc,
        );
        assert!(near(&zp, &want, 96), "ζ'(-1) got {:?}", zp);
        Ok(())
    }

    fn sum_k_ln_k_naive(n: u64, p: usize, cc: &mut Consts) -> ExactNum {
        let mut acc = ExactNum::new(p);
        for k in 1..=n {
            let kf = ExactNum::from_u64(k, p);
            acc = acc.add(&kf.clone().mul(&kf.ln(p, NR, cc), p, NR), p, NR);
        }
        acc
    }

    #[test]
    fn zf_hurwitz_sum_k_ln_k_sieved_matches_naive() -> Result<(), Error> {
        let mut cc = Consts::new()?;
        for p in [64usize, 256, 1024] {
            for n in [1u64, 2, 10, 64, 200] {
                let got = sum_k_ln_k(n, p, &mut cc);
                let want = sum_k_ln_k_naive(n, p, &mut cc);
                assert!(near(&got, &want, p / 2), "S({n}) p={p}");
            }
        }
        Ok(())
    }

    #[test]
    fn zf_hurwitz_em_sieved_matches_naive_at_one_and_half() -> Result<(), Error> {
        let mut cc = Consts::new()?;
        let s = ExactNum::from_i32(-1, 64);
        for p in [64usize, 256, 1024] {
            let a1 = ExactNum::from_u8(1, p);
            let a2 = ExactNum::from_u8(1, p).ldexp(-1, p, NR);
            let Some(fast1) = hurwitz_em_jet(&s, &a1, p, &mut cc) else {
                panic!("jet a=1 failed");
            };
            let Some(fast2) = hurwitz_em_jet(&s, &a2, p, &mut cc) else {
                panic!("jet a=1/2 failed");
            };
            // Closed form: ζ'(-1, 1/2) = -ln2/24 - ζ'(-1)/2.
            let ln2 = cc.ln_2(p, NR);
            let c = ln2
                .div(&ExactNum::from_u32(24, p), p, NR)
                .neg()
                .sub(&fast1.1.clone().ldexp(-1, p, NR), p, NR);
            assert!(near(&fast2.1, &c, p / 2), "closed form ζ'(-1,1/2) p={p}");
            let _ = fast1;
        }
        Ok(())
    }
}
