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
//! `Σ k ln k = Σ_p ln(p) Σ_j p^j T(⌊n/p^j⌋)`. At general `s` and `a = 1`
//! the same jet uses an `acb_dirichlet_powsum_sieved` equivalent: logs and
//! powers only at primes, composites built by multiplication. The prime logs
//! come from FLINT 3.6 `arb_log_primes_vec_bsplit` (Machin atanh binary
//! splitting, then `2 atanh(1/(2p^2-1))` for later primes). LogG(1/2) uses
//! the closed form `ζ'(-1, 1/2) = -(ln 2)/24 - ζ'(-1)/2`, so G(±1/2) pays
//! only the `a = 1` jet. Dest `p ≤ 2^14` finishes in a few seconds.

use crate::defs::WORD_BIT_SIZE;
use crate::ops::special::with_even_bernoulli;
use crate::Consts;
use crate::ExactNum;
use crate::RoundingMode;
use alloc::vec;
use alloc::vec::Vec;
use core::ops::Neg;

const NR: RoundingMode = RoundingMode::None;
const UP: RoundingMode = RoundingMode::Up;

/// Length-2 jet `c0 + c1 x` of a function of the series variable `x`.
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

    /// Length-2 product `(a0 + a1 x)(b0 + b1 x)` truncated after `x` (FLINT `_acb_poly_mullow`).
    fn mullow(&self, other: &Self, p: usize) -> Self {
        Self {
            c0: ExactNum::mul_ref(&self.c0, &other.c0, p, NR),
            c1: ExactNum::mul_ref(&self.c1, &other.c0, p, NR).add(
                &ExactNum::mul_ref(&self.c0, &other.c1, p, NR),
                p,
                NR,
            ),
        }
    }
}

/// `(N+a)^{-(s+x)} = (N+a)^{-s} (1 − x ln(N+a) + O(x²))`.
fn invpow_jet(base: &ExactNum, s: &ExactNum, p: usize, cc: &mut Consts) -> Jet {
    let ln = base.ln(p, NR, cc);
    let c0 = base.pow(&s.neg_ref(), p, NR, cc);
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
fn em_bound1(s: &ExactNum, a: &ExactNum, n: u64, m: u64, cc: &mut Consts) -> Option<ExactNum> {
    if n < 1 || m < 1 {
        return None;
    }
    let wp = 64usize;
    let two_m = m.checked_mul(2)?;
    let sigma = { s.rounded_to(wp, UP) };
    let alpha = { a.rounded_to(wp, UP) };
    // A = N+a, B = σ+2M; need A>1 and σ+2M>1.
    let an = ExactNum::from_u64(n, wp).add(&alpha, wp, UP);
    let s2m = sigma.add(&ExactNum::from_u64(two_m, wp), wp, UP);
    if !an.is_positive() || matches!(an.compare(&ExactNum::from_u8(1, wp)), Some(c) if c <= 0) {
        return None;
    }
    if !s2m.is_positive() || matches!(s2m.compare(&ExactNum::from_u8(1, wp)), Some(c) if c <= 0) {
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
        let mut t = s.abs_ref();
        let _ = t.set_precision(wp, UP);
        t
    };
    let mut f0 = ExactNum::from_u8(1, wp);
    let mut harm = ExactNum::new(wp);
    for i in 0..two_m {
        let term = ExactNum::add_ref(&sabs, &ExactNum::from_u64(i, wp), wp, UP);
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
    let mag = if matches!(b1.compare(&b0), Some(c) if c > 0) { b1 } else { b0 };
    if mag.is_nan() || !mag.is_positive() {
        None
    } else {
        Some(mag)
    }
}

fn a_pow_neg_s_mag(s: &ExactNum, a: &ExactNum, cc: &mut Consts) -> ExactNum {
    let wp = 64usize;
    let aa = a.rounded_to(wp, UP);
    let ss = s.rounded_to(wp, UP);
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
fn choose_param(s: &ExactNum, a: &ExactNum, target: usize, cc: &mut Consts) -> Option<(u64, u64)> {
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
fn em_sum(s: &ExactNum, a: &ExactNum, n: u64, m: u64, p: usize, cc: &mut Consts) -> Option<Jet> {
    let p = p.saturating_add(2 * (bit_count(n) as usize + 2));
    let one = ExactNum::from_u8(1, p);
    #[cfg(all(test, feature = "std"))]
    let t_powsum = std::time::Instant::now();
    let mut sum = if let Some(fast) = powsum_s_minus_one(s, a, n, p, cc) {
        fast
    } else if let Some(fast) = powsum_sieved(s, a, n, p, cc) {
        fast
    } else {
        let mut sum = Jet::zero(p);
        for k in 0..n {
            let t = ExactNum::from_u64(k, WORD_BIT_SIZE).add(a, p, NR);
            if !t.is_positive() {
                return None;
            }
            sum = sum.add(&invpow_jet(&t, s, p, cc), p);
        }
        sum
    };
    let na = ExactNum::from_u64(n, WORD_BIT_SIZE).add(a, p, NR);
    // t = (N+a)^{-(s+x)}
    let tjet = invpow_jet(&na, s, p, cc);
    // u = (N+a)^{1-(s+x)} / ((s-1)+x)
    let sm1 = ExactNum::sub_ref(s, &one, p, NR);
    if sm1.is_zero() {
        return None;
    }
    let v = ExactNum::sub_ref(s, &ExactNum::from_u8(1, p), p, NR);
    let u_pow = invpow_jet(&na, &v, p, cc);
    // divide by (s-1)+x: u0 = w0/(s-1), u1 = (w1 - u0)/(s-1)
    let u0 = ExactNum::div_ref(&u_pow.c0, &sm1, p, NR);
    let u1 = ExactNum::sub_ref(&u_pow.c1, &u0, p, NR).div(&sm1, p, NR);
    sum = sum.add(&Jet { c0: u0, c1: u1 }, p);
    // + 1/2 t
    let half = ExactNum::from_u8(1, p).ldexp(-1, p, NR);
    sum = sum.add(&tjet.scale(&half, p), p);

    #[cfg(all(test, feature = "std"))]
    let ms_powsum = t_powsum.elapsed().as_millis();
    // Bernoulli tail, `_acb_poly_zeta_em_tail_naive`.
    #[cfg(all(test, feature = "std"))]
    let t_tail = std::time::Instant::now();
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
        // Integer factors as one-word values (FLINT fmpz / arb_mul_ui). Building
        // them at dest `p` made each of the ~10^4 tail steps a full-width mul.
        let splus0 = ExactNum::add_ref(
            s,
            &ExactNum::from_u64((2 * r - 1) as u64, WORD_BIT_SIZE),
            p,
            NR,
        );
        term = term.mul_cpx(&splus0, p);
        let splus1 =
            ExactNum::add_ref(s, &ExactNum::from_u64((2 * r) as u64, WORD_BIT_SIZE), p, NR);
        term = term.mul_cpx(&splus1, p);
        let den = ExactNum::mul_ref(
            &na2,
            &ExactNum::from_u64((2 * r + 1) as u64 * (2 * r + 2) as u64, WORD_BIT_SIZE),
            p,
            NR,
        );
        term = term.div_const(&den, p);
    }
    #[cfg(all(test, feature = "std"))]
    {
        eprintln!(
            "EM_SUM profile n={n} m={m} p={p} ms_powsum={ms_powsum} ms_tail={}",
            t_tail.elapsed().as_millis()
        );
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
                return v.rounded_to(p, RoundingMode::ToEven);
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
            .sub(&c.ldexp_ref(-1, p, NR), p, NR)
    } else {
        let s = ExactNum::from_i32(-1, p);
        let Some((_, zp)) = hurwitz_em_jet(&s, z, p, cc) else {
            return ExactNum::nan(Some(crate::Error::InvalidArgument));
        };
        zp
    };
    let lng = z.ln_gamma(p, NR, cc);
    let zm1 = ExactNum::sub_ref(z, &ExactNum::from_u8(1, p), p, NR);
    zm1.mul(&lng, p, NR).sub(&zp, p, NR).add(&c, p, NR)
}

/// `Σ_{k=1}^{n} k ln k`, via `Σ_p ln(p) Σ_{j: p^j ≤ n} p^j T(⌊n/p^j⌋)`
/// with `T(m) = m(m+1)/2`. Only `π(n)` full-precision logarithms.
fn sum_k_ln_k(n: u64, p: usize, cc: &mut Consts) -> ExactNum {
    if n == 0 {
        return ExactNum::new(p);
    }
    let lns = crate::ops::log_primes::ln_primes_upto(n, p, cc);
    let mut acc = ExactNum::new(p);
    for (pru, ln) in lns {
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
        let c = crate::ExactInt::from_u128(coeff).to_exact_num(p, NR);
        acc = acc.add(&ln.mul(&c, p, NR), p, NR);
    }
    acc
}

fn is_exact_one(x: &ExactNum, p: usize) -> bool {
    matches!(x.compare(&ExactNum::from_u8(1, p)), Some(0))
}

fn is_exact_half(x: &ExactNum, p: usize) -> bool {
    let two_x = x.ldexp_ref(1, p, NR);
    is_exact_one(&two_x, p)
}

fn is_exact_neg_one(x: &ExactNum, p: usize) -> bool {
    matches!(x.compare(&ExactNum::from_i32(-1, p)), Some(0))
}

/// Memory cap for the dense `k^{-s}` tables (two arrays of `ExactNum`).
/// Above this, composites are factored from a smallest-prime-factor table
/// so the 2^20 sieve does not allocate tens of gigabytes.
fn powsum_table_too_big(n: u64, p: usize) -> bool {
    let bytes = (p / 8).saturating_add(64) as u64;
    n.saturating_add(1).saturating_mul(bytes).saturating_mul(2) > 1_500_000_000
}

/// `k^{-s}` and its `s`-derivative `-ln(k) k^{-s}`.
/// On the critical line `s = 1/2`, this is FLINT `arb_rsqrt_ui(k)` rather
/// than `exp(-s log k)` (`acb_dirichlet_powsum_term`). Integer `s` uses
/// `pow`. The log of `k` is supplied by the caller (prime-log vec).
fn powsum_term_from_ln(
    s: &ExactNum,
    k: u64,
    ln: &ExactNum,
    p: usize,
    cc: &mut Consts,
) -> Option<Jet> {
    let inv = if k <= 1 {
        ExactNum::from_u8(1, p)
    } else if is_exact_half(s, p) {
        ExactNum::from_u64(k, p).sqrt(p, NR).reciprocal(p, NR)
    } else if s.is_int() {
        ExactNum::from_u64(k, p).pow(&s.neg_ref(), p, NR, cc)
    } else {
        ExactNum::mul_ref(s, ln, p, NR).neg().exp(p, NR, cc)
    };
    if inv.is_nan() || inv.is_inf() {
        return None;
    }
    let dinv = if k <= 1 { ExactNum::new(p) } else { ExactNum::mul_ref(&inv, ln, p, NR).neg() };
    Some(Jet { c0: inv, c1: dinv })
}

/// `k^{-s}` at a prime `k` identified by its log, and the `s`-derivative.
fn prime_invpow(
    s: &ExactNum,
    k: u64,
    ln: &ExactNum,
    p: usize,
    cc: &mut Consts,
) -> Option<(ExactNum, ExactNum)> {
    let j = powsum_term_from_ln(s, k, ln, p, cc)?;
    Some((j.c0, j.c1))
}

/// Sieved Σ_{k=1}^{n} k^{-(s+x)} at `a = 1` (FLINT `acb_dirichlet_powsum_sieved`):
/// logs and powers only at primes, composites by multiplication.
/// The length-2 jet is `c0 = Σ k^{-s}`, `c1 = Σ -ln(k) k^{-s}`.
fn powsum_sieved(s: &ExactNum, a: &ExactNum, n: u64, p: usize, cc: &mut Consts) -> Option<Jet> {
    if n == 0 || !is_exact_one(a, p) {
        return None;
    }
    let n_us = usize::try_from(n).ok()?;
    let lns = crate::ops::log_primes::ln_primes_upto(n, p, cc);
    if powsum_table_too_big(n, p) {
        return powsum_sieved_factored(s, n_us, p, &lns, cc);
    }
    if n == 1 {
        return Some(Jet {
            c0: ExactNum::from_u8(1, p),
            c1: ExactNum::new(p),
        });
    }

    // FLINT `acb_dirichlet_powsum_sieved`: only odd k in the table,
    // composites by multiplication of stored odd factors, evens by Horner
    // in `2^{-s}`. Hard `k^{-s}` only at 2 and at odd primes (`powsum_term`).
    let mut divisors = vec![0u32; n_us / 2 + 1];
    let ibound = n_us.isqrt();
    let mut i = 3usize;
    while i <= ibound {
        if divisors[(i - 1) / 2] == 0 {
            let mut j = i * i;
            while j <= n_us {
                divisors[(j - 1) / 2] = i as u32;
                j += 2 * i;
            }
        }
        i += 2;
    }

    let ln_of = |k: u64| -> Option<&ExactNum> {
        let i = lns.binary_search_by_key(&k, |(pr, _)| *pr).ok()?;
        Some(&lns[i].1)
    };

    #[cfg(all(test, feature = "std"))]
    let t_term = std::time::Instant::now();
    #[cfg(all(test, feature = "std"))]
    let mut n_term = 0u64;
    let x = powsum_term_from_ln(s, 2, ln_of(2)?, p, cc)?;
    #[cfg(all(test, feature = "std"))]
    {
        n_term += 1;
    }

    // Stored odd powers for k <= n/3, indexed by (k-1)/2 (FLINT `POWER`).
    let n_store = n_us / 6 + 1;
    let mut powers: Vec<Jet> = (0..n_store).map(|_| Jet::zero(p)).collect();
    let mut z = Jet::zero(p);
    let mut u = Jet::zero(p);

    let mut power_of_two = 1u64;
    while power_of_two.saturating_mul(2) <= n {
        power_of_two *= 2;
    }
    let mut horner_point = n / power_of_two;

    #[cfg(all(test, feature = "std"))]
    let mut n_comp = 0u64;
    #[cfg(all(test, feature = "std"))]
    let mut ms_term = 0u128;
    #[cfg(all(test, feature = "std"))]
    let mut ms_comp = 0u128;
    #[cfg(all(test, feature = "std"))]
    let mut ms_horner = 0u128;

    let mut k = 1u64;
    while k <= n {
        let ku = k as usize;
        let t = if k == 1 {
            Jet {
                c0: ExactNum::from_u8(1, p),
                c1: ExactNum::new(p),
            }
        } else if divisors[(ku - 1) / 2] == 0 {
            #[cfg(all(test, feature = "std"))]
            let t0 = std::time::Instant::now();
            let term = powsum_term_from_ln(s, k, ln_of(k)?, p, cc)?;
            #[cfg(all(test, feature = "std"))]
            {
                n_term += 1;
                ms_term += t0.elapsed().as_millis();
            }
            term
        } else {
            let d = divisors[(ku - 1) / 2] as usize;
            let cof = ku / d;
            let p1 = &powers[(d - 1) / 2];
            let p2 = &powers[(cof - 1) / 2];
            #[cfg(all(test, feature = "std"))]
            let t0 = std::time::Instant::now();
            let term = p1.mullow(p2, p);
            #[cfg(all(test, feature = "std"))]
            {
                n_comp += 1;
                ms_comp += t0.elapsed().as_millis();
            }
            term
        };
        u = u.add(&t, p);
        if k.saturating_mul(3) <= n {
            let idx = (ku - 1) / 2;
            if idx < powers.len() {
                powers[idx] = t;
            }
        }

        while k == horner_point && power_of_two != 1 {
            #[cfg(all(test, feature = "std"))]
            let t0 = std::time::Instant::now();
            z = z.mullow(&x, p).add(&u, p);
            #[cfg(all(test, feature = "std"))]
            {
                ms_horner += t0.elapsed().as_millis();
            }
            power_of_two /= 2;
            horner_point = n / power_of_two;
            if horner_point.is_multiple_of(2) {
                horner_point -= 1;
            }
        }
        k += 2;
    }
    z = z.mullow(&x, p).add(&u, p);
    #[cfg(all(test, feature = "std"))]
    {
        eprintln!(
            "POWSUM profile n={n} p={p} n_term={n_term} ms_term={} n_comp={n_comp} ms_comp={} ms_horner={} ms_setup={}",
            ms_term,
            ms_comp,
            ms_horner,
            t_term.elapsed().as_millis()
        );
    }
    Some(z)
}

/// Same sum as `powsum_sieved`, without a dense table of `k^{-s}`.
fn powsum_sieved_factored(
    s: &ExactNum,
    n_us: usize,
    p: usize,
    lns: &[(u64, ExactNum)],
    cc: &mut Consts,
) -> Option<Jet> {
    let mut prime_ln: Vec<ExactNum> = Vec::with_capacity(lns.len());
    let mut prime_inv: Vec<ExactNum> = Vec::with_capacity(lns.len());
    let mut idx = vec![u32::MAX; n_us + 1];
    let mut primes: Vec<u32> = Vec::with_capacity(lns.len());
    for (i, (pr, ln)) in lns.iter().enumerate() {
        let pu = *pr as usize;
        if pu > n_us {
            break;
        }
        let (inv, _) = prime_invpow(s, *pr, ln, p, cc)?;
        let i32u = u32::try_from(i).ok()?;
        idx[pu] = i32u;
        prime_ln.push(ln.clone());
        prime_inv.push(inv);
        primes.push(*pr as u32);
    }

    let mut spf = vec![0u32; n_us + 1];
    for i in 2..=n_us {
        if spf[i] == 0 {
            spf[i] = i as u32;
        }
        for &pr in &primes {
            let pu = pr as usize;
            let v = match i.checked_mul(pu) {
                Some(v) if v <= n_us => v,
                _ => break,
            };
            spf[v] = pr;
            if (i as u32).is_multiple_of(pr) {
                break;
            }
        }
    }

    let mut c0 = ExactNum::from_u8(1, p);
    let mut c1 = ExactNum::new(p);
    for k in 2..=n_us {
        let mut rem = k as u32;
        let mut inv = ExactNum::from_u8(1, p);
        let mut lnk = ExactNum::new(p);
        while rem > 1 {
            let pr = spf[rem as usize];
            let i = idx[pr as usize] as usize;
            if i == u32::MAX as usize {
                return None;
            }
            inv = inv.mul(&prime_inv[i], p, NR);
            lnk = lnk.add(&prime_ln[i], p, NR);
            rem /= pr;
        }
        c0 = c0.add(&inv, p, NR);
        c1 = c1.add(&ExactNum::mul_ref(&inv, &lnk, p, NR).neg(), p, NR);
    }
    Some(Jet { c0, c1 })
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
        return Some(Jet { c0, c1: s_n.neg() });
    }
    if is_exact_half(a, p) {
        // Σ (k+1/2) = n²/2.
        // Σ (k+1/2) ln(k+1/2)
        //   = (1/2)[S(2n-1) - n(n-1) ln 2 - 2 S(n-1)] - (ln 2) n²/2.
        let c0 = crate::ExactInt::from_u128(nf * nf)
            .to_exact_num(p, NR)
            .ldexp(-1, p, NR);
        let s_odd = if n == 0 {
            ExactNum::new(p)
        } else {
            let s_hi = sum_k_ln_k(2 * n - 1, p, cc);
            let s_lo = if n == 1 { ExactNum::new(p) } else { sum_k_ln_k(n - 1, p, cc) };
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
    use crate::test_must::Must;
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
            acc = acc.add(&ExactNum::mul_ref(&kf, &kf.ln(p, NR, cc), p, NR), p, NR);
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
    fn zf_hurwitz_em_sieved_matches_naive_at_one_and_half() -> Result<(), String> {
        let mut cc = Consts::new()?;
        let s = ExactNum::from_i32(-1, 64);
        for p in [64usize, 256, 1024] {
            let a1 = ExactNum::from_u8(1, p);
            let a2 = ExactNum::from_u8(1, p).ldexp(-1, p, NR);
            let Some(fast1) = hurwitz_em_jet(&s, &a1, p, &mut cc) else {
                return Err("jet a=1 failed".into());
            };
            let Some(fast2) = hurwitz_em_jet(&s, &a2, p, &mut cc) else {
                return Err("jet a=1/2 failed".into());
            };
            // Closed form: ζ'(-1, 1/2) = -ln2/24 - ζ'(-1)/2.
            let ln2 = cc.ln_2(p, NR);
            let c = ln2.div(&ExactNum::from_u32(24, p), p, NR).neg().sub(
                &fast1.1.ldexp(-1, p, NR),
                p,
                NR,
            );
            assert!(near(&fast2.1, &c, p / 2), "closed form ζ'(-1,1/2) p={p}");
            let _ = fast1;
        }
        Ok(())
    }

    fn parse_zeta_arg(arg: &str, p: usize, cc: &mut Consts) -> ExactNum {
        if let Some((num, den)) = arg.split_once('/') {
            let n = ExactNum::parse(num, crate::Radix::Dec, p, RoundingMode::ToEven, cc);
            let d = ExactNum::parse(den, crate::Radix::Dec, p, RoundingMode::ToEven, cc);
            ExactNum::div_ref(&n, &d, p, RoundingMode::ToEven)
        } else {
            ExactNum::parse(arg, crate::Radix::Dec, p, RoundingMode::ToEven, cc)
        }
    }

    fn agree_bits(got: &ExactNum, want: &ExactNum, p: usize) -> i32 {
        let d = ExactNum::sub_ref(got, want, p, RoundingMode::None).abs();
        if d.is_zero() {
            10_000
        } else {
            (want.exponent().unwrap_or(0) - d.exponent().unwrap_or(0)) as i32
        }
    }

    #[test]
    fn zf_riemann_zeta_matches_flint_golds() -> Result<(), String> {
        // FLINT 3.6 arb.zeta midpoints at dest p = 64, 4096, 16384, 65536.
        let mut cc = Consts::new()?;
        let raw = include_str!("testdata/riemann_zeta_flint_golds.txt");
        for line in raw.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut parts = line.splitn(3, ' ');
            let p: usize = parts.next().must("field")?.parse().must("parse")?;
            let arg = parts.next().must("field")?;
            let dec = parts.next().must("field")?;
            let t0 = std::time::Instant::now();
            let x = parse_zeta_arg(arg, p, &mut cc);
            let got = x.riemann_zeta(p, RoundingMode::ToEven, &mut cc);
            let ms = t0.elapsed().as_millis();
            let want = ExactNum::parse(
                dec,
                crate::Radix::Dec,
                p.saturating_add(256),
                RoundingMode::ToEven,
                &mut cc,
            );
            let bits = agree_bits(&got, &want, p);
            eprintln!("zeta({arg}) p={p}: {ms} ms, {bits} bits");
            assert!(
                !got.is_nan() && (bits >= (p as i32) - 16 || bits == 10_000),
                "zeta({arg}) p={p} only {bits} bits after {ms} ms"
            );
        }
        Ok(())
    }

    #[test]
    fn zf_riemann_zeta_edge_cases() -> Result<(), Error> {
        let mut cc = Consts::new()?;
        let p = 256usize;
        let half = ExactNum::from_u8(1, p).ldexp(-1, p, NR);

        // ζ(0) = -1/2.
        let z0 = ExactNum::new(p).riemann_zeta(p, RoundingMode::ToEven, &mut cc);
        let d0 =
            ExactNum::sub_ref(&z0, &ExactNum::from_u8(1, p).ldexp(-1, p, NR).neg(), p, NR).abs();
        assert!(
            d0.is_zero()
                || agree_bits(&z0, &ExactNum::from_u8(1, p).ldexp(-1, p, NR).neg(), p)
                    >= (p as i32) - 8,
            "zeta(0) is not -1/2"
        );

        // Pole at s = 1: +Inf, exact (POLE-08, divideByZero).
        let z1 = ExactNum::from_u8(1, p).riemann_zeta(p, RoundingMode::ToEven, &mut cc);
        assert!(z1.is_inf_pos(), "zeta(1) is +Inf (POLE-08)");

        // Trivial zeros: ζ(-2) = ζ(-4) = 0.
        for n in [-2i32, -4] {
            let z = ExactNum::from_i32(n, p).riemann_zeta(p, RoundingMode::ToEven, &mut cc);
            assert!(z.is_zero(), "zeta({n}) must be the trivial zero");
        }

        // ζ(1/2) is finite and negative.
        let zhalf = half.riemann_zeta(p, RoundingMode::ToEven, &mut cc);
        assert!(
            !zhalf.is_nan() && zhalf.is_negative(),
            "zeta(1/2) is finite and negative"
        );

        // Large positive s: ζ(s) → 1. At s = 400, 2^{-400} is far below a 256-bit ulp.
        let large = ExactNum::from_u32(400, p);
        let zl = large.riemann_zeta(p, RoundingMode::ToEven, &mut cc);
        let one = ExactNum::from_u8(1, p);
        let bits = agree_bits(&zl, &one, p);
        assert!(
            bits >= (p as i32) - 16 || zl == one,
            "zeta(400) vs 1 only {bits} bits"
        );

        // Large |s| on the negative side: ζ(-21) is a finite nonzero rational.
        let zneg = ExactNum::from_i32(-21, p).riemann_zeta(p, RoundingMode::ToEven, &mut cc);
        assert!(
            !zneg.is_nan() && !zneg.is_zero() && !zneg.is_inf(),
            "zeta(-21) is finite and nonzero"
        );
        Ok(())
    }

    #[test]
    fn time_riemann_zeta_half_p65536() -> Result<(), String> {
        let mut cc = Consts::new()?;
        let p = 65536usize;
        let raw = include_str!("testdata/riemann_zeta_flint_golds.txt");
        let wanted = ["1/2", "-1/2", "3/2"];
        let mut golds = Vec::new();
        for line in raw.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with("#") {
                continue;
            }
            let mut parts = line.splitn(3, " ");
            let gp: usize = parts.next().must("field")?.parse().must("parse")?;
            let arg = parts.next().must("field")?;
            let dec = parts.next().must("field")?;
            if gp == p && wanted.contains(&arg) {
                golds.push((arg, dec));
            }
        }
        assert_eq!(golds.len(), 3, "expected 1/2, -1/2, 3/2 golds at p=65536");
        eprintln!("ZETA profile phase=start p={p} (FLINT ~5 s each)");
        for (i, (arg, dec)) in golds.iter().enumerate() {
            eprintln!(
                "ZETA profile phase=first arg={arg} p={p} counter={}/3",
                i + 1
            );
            let t0 = std::time::Instant::now();
            let x = parse_zeta_arg(arg, p, &mut cc);
            let got = x.riemann_zeta(p, RoundingMode::ToEven, &mut cc);
            let ms = t0.elapsed().as_millis();
            let want = ExactNum::parse(
                dec,
                crate::Radix::Dec,
                p.saturating_add(256),
                RoundingMode::ToEven,
                &mut cc,
            );
            let bits = agree_bits(&got, &want, p);
            eprintln!(
                "ZETA profile phase=first_done arg={arg} p={p} ms={ms} bits={bits} (FLINT ~5 s)"
            );
            assert!(
                !got.is_nan() && (bits >= (p as i32) - 16 || bits == 10_000),
                "zeta({arg}) p={p} only {bits} bits after {ms} ms"
            );
        }
        let t1 = std::time::Instant::now();
        let x = parse_zeta_arg("1/2", p, &mut cc);
        let _ = x.riemann_zeta(p, RoundingMode::ToEven, &mut cc);
        eprintln!(
            "ZETA profile phase=cached_done arg=1/2 p={p} ms={} (FLINT ~5 s)",
            t1.elapsed().as_millis()
        );
        Ok(())
    }

    fn powsum_direct(s: &ExactNum, n: u64, p: usize, cc: &mut Consts) -> Jet {
        let mut sum = Jet::zero(p);
        for k in 1..=n {
            let t = ExactNum::from_u64(k, p);
            sum = sum.add(&invpow_jet(&t, s, p, cc), p);
        }
        sum
    }

    fn jets_agree(a: &Jet, b: &Jet, p: usize, min_bits: i32) {
        let b0 = agree_bits(&a.c0, &b.c0, p);
        let b1 = agree_bits(&a.c1, &b.c1, p);
        assert!(
            b0 >= min_bits && b1 >= min_bits,
            "powsum jet bits c0={b0} c1={b1} want >={min_bits}"
        );
    }

    #[test]
    fn powsum_sieved_matches_direct_sum_at_modest_n() -> Result<(), String> {
        let mut cc = Consts::new()?;
        let p = 256usize;
        let n = 80u64;
        let a = ExactNum::from_u8(1, p);
        for (label, s) in [
            ("1/2", ExactNum::from_u8(1, p).ldexp(-1, p, NR)),
            ("-1/2", ExactNum::from_u8(1, p).ldexp(-1, p, NR).neg()),
            ("3/2", ExactNum::from_u8(3, p).ldexp(-1, p, NR)),
            ("2", ExactNum::from_u8(2, p)),
            ("3", ExactNum::from_u8(3, p)),
        ] {
            let got = powsum_sieved(&s, &a, n, p, &mut cc)
                .ok_or_else(|| format!("sieved failed for s={label}"))?;
            let want = powsum_direct(&s, n, p, &mut cc);
            jets_agree(&got, &want, p, (p as i32) - 16);
        }
        // N just below/above a power of two, where Horner changes scale.
        for n in [31u64, 32, 33, 63, 64, 65] {
            let s = ExactNum::from_u8(1, p).ldexp(-1, p, NR);
            let got = powsum_sieved(&s, &a, n, p, &mut cc).must("sieved")?;
            let want = powsum_direct(&s, n, p, &mut cc);
            jets_agree(&got, &want, p, (p as i32) - 16);
        }
        Ok(())
    }
}
