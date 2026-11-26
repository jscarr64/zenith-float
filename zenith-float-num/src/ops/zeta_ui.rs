//! Integer zeta for odd positive `n`, after FLINT 3.6 `arb_zeta_ui`
//! (Fredrik Johansson, LGPL-3.0-or-later): Euler product, Borwein
//! binary splitting, and the asymptotic product for `n > 0.7 prec`.
//! Even `n` stay on the existing Bernoulli path in `catalog`. `n = 3`
//! uses the same Borwein splitter FLINT uses past the Apery special
//! case; at dest p ≤ 2^20 that is still the Borwein / high-precision
//! branch (`n < prec * 0.0006` is false only at modest p, where the
//! vec-Borwein fallback is the same series).

use crate::defs::RoundingMode;
use crate::ext::ExactNum;
use crate::Consts;
use alloc::vec;
use alloc::vec::Vec;
use core::ops::Neg;

const NR: RoundingMode = RoundingMode::None;

/// Reciprocal of FLINT `ERROR_B` = `log2(3+sqrt(8))`, as 10000/25431
/// (1/2.543106606327224 ≈ 0.39322).
const BORWEIN_N_NUM: u64 = 10_000;
const BORWEIN_N_DEN: u64 = 25_431;

/// `arb_zeta_ui` for a positive integer `n ≠ 1`, at precision `p`.
pub(crate) fn zeta_ui(n: u64, p: usize, _cc: &mut Consts) -> ExactNum {
    if n == 0 {
        return ExactNum::from_u8(1, p).ldexp(-1, p, NR).neg();
    }
    if n == 1 {
        return ExactNum::nan(Some(crate::Error::InvalidArgument));
    }
    let prec = p as u64;
    // FLINT: n > 0.7 * prec → Euler product (called zeta_ui_asymp).
    if n.saturating_mul(10) > prec.saturating_mul(7) {
        return zeta_ui_euler_product(n, p);
    }
    if n.is_multiple_of(2) {
        // Even n is handled by the caller (Bernoulli). This is a
        // fallback if someone asks for a large even n here.
        return zeta_ui_euler_product(n, p);
    }
    // Odd n. FLINT: n == 3 → Apery; n < prec*0.0006 → Borwein bsplit
    // (small odd n, extremely high precision); n > euler_product_cutoff
    // → Euler product; else vec-Borwein. A single-value Borwein bsplit
    // is the right method for n = 3 and n = 5 at dest p ≥ 64, and for
    // any odd n below the Euler-product cutoff.
    // n < prec * 0.0006, integer: 10000 n < 6 prec.
    if n.saturating_mul(10_000) < prec.saturating_mul(6) || n == 3 {
        return zeta_ui_borwein_bsplit(n, p);
    }
    if n > euler_product_cutoff(prec) {
        return zeta_ui_euler_product(n, p);
    }
    zeta_ui_borwein_bsplit(n, p)
}

/// Integer form of FLINT `euler_product_cutoff` in `zeta_ui.c`:
/// `7 + 0.535 * prec / ln(prec)`, with `ln(prec) ≈ ln(2) · log2(prec)`
/// and `ln(2) ≈ 545/786`. The mid-precision "bulge" `0.39 prec^{0.8}`
/// is within a few percent of this expression and is not needed for
/// a correct method choice.
fn euler_product_cutoff(prec: u64) -> u64 {
    let log2 = (u64::BITS - prec.leading_zeros()) as u64;
    if log2 == 0 {
        return 7;
    }
    // 0.535 * prec / (ln2 * log2) ≈ 420 * prec / (545 * log2)
    7 + (420 * prec) / (545 * log2)
}

/// `ζ(n) = ∏_q 1/(1 − q^{-n})` over primes, after
/// `arb_zeta_ui_euler_product` / `_acb_dirichlet_euler_product_real_ui`.
fn zeta_ui_euler_product(n: u64, p: usize) -> ExactNum {
    let wp = p.saturating_add(32);
    let one = ExactNum::from_u8(1, wp);
    // Stop when q^n ≥ 2^{wp+8}.
    let pmax = euler_prime_limit(n, wp);
    let primes = primes_upto(pmax);
    let mut z = ExactNum::from_u8(1, wp);
    for q in primes {
        let qn = pow_ui(&ExactNum::from_u64(q, wp), n, wp);
        let den = ExactNum::sub_ref(&qn, &one, wp, NR);
        if den.is_zero() || den.is_nan() || den.is_inf() {
            return ExactNum::nan(Some(crate::Error::InvalidArgument));
        }
        let term = qn.div(&den, wp, NR);
        z = z.mul(&term, wp, NR);
        if z.is_nan() || z.is_inf() {
            return ExactNum::nan(Some(crate::Error::InvalidArgument));
        }
    }
    let mut out = z;
    let _ = out.set_precision(p, RoundingMode::None);
    out
}

/// Largest prime `q` with `q^n < 2^{wp+8}`.
fn euler_prime_limit(n: u64, wp: usize) -> u64 {
    if n == 0 {
        return 2;
    }
    // 2^{ceil((wp+8)/n)} is an upper bound on 2^{(wp+8)/n}.
    let exp = wp.saturating_add(8).div_ceil(n as usize);
    if exp >= 63 {
        return u64::MAX / 4;
    }
    (1u64 << exp).saturating_add(2).max(3)
}

fn primes_upto(max: u64) -> Vec<u64> {
    if max < 2 {
        return Vec::new();
    }
    let m = usize::try_from(max).unwrap_or(usize::MAX / 4);
    let mut is_p = vec![true; m + 1];
    is_p[0] = false;
    is_p[1] = false;
    let mut p = 2usize;
    while p.saturating_mul(p) <= m {
        if is_p[p] {
            let mut k = p.saturating_mul(p);
            while k <= m {
                is_p[k] = false;
                k = k.saturating_add(p);
            }
        }
        p += 1;
    }
    (2..=m).filter(|&i| is_p[i]).map(|i| i as u64).collect()
}

fn pow_ui(base: &ExactNum, exp: u64, p: usize) -> ExactNum {
    if exp == 0 {
        return ExactNum::from_u8(1, p);
    }
    let mut acc = ExactNum::from_u8(1, p);
    // base^(2^j); `None` while j = 0 (base read in place).
    let mut b_slot: Option<ExactNum> = None;
    let mut e = exp;
    loop {
        let b = b_slot.as_ref().unwrap_or(base);
        if e & 1 == 1 {
            acc = acc.mul(b, p, NR);
        }
        e >>= 1;
        if e == 0 {
            return acc;
        }
        let sq = ExactNum::mul_ref(b, b, p, NR);
        b_slot = Some(sq);
    }
}

struct ZetaSplit {
    a: ExactNum,
    b: ExactNum,
    c: ExactNum,
    d: ExactNum,
    q1: ExactNum,
    q2: ExactNum,
    q3: ExactNum,
}

/// Leaf of FLINT `zeta_coeff_k` at working precision `wp`.
fn zeta_coeff_k(k: u64, n: u64, s: u64, wp: usize) -> ZetaSplit {
    // D = 2 (n+k) (n-k)
    let d = ExactNum::from_u64(2 * (n + k) * (n - k), wp);
    // Q1 = (k+1) (2k+1)
    let q1 = ExactNum::from_u64((k + 1) * (2 * k + 1), wp);
    let (a, q2) = if k == 0 {
        (ExactNum::new(wp), ExactNum::from_u8(1, wp))
    } else {
        let mut a = q1.clone();
        if k.is_multiple_of(2) {
            a = a.neg();
        }
        (a, pow_ui(&ExactNum::from_u64(k, wp), s, wp))
    };
    let q3 = ExactNum::mul_ref(&q1, &q2, wp, NR);
    ZetaSplit {
        a,
        b: ExactNum::new(wp),
        c: q1.clone(),
        d,
        q1,
        q2,
        q3,
    }
}

fn zeta_bsplit(a: u64, b: u64, n: u64, s: u64, cont: bool, wp: usize) -> ZetaSplit {
    if a + 1 == b {
        return zeta_coeff_k(a, n, s, wp);
    }
    let m = (a + b) / 2;
    let mut l = zeta_bsplit(m, b, n, s, true, wp);
    let r = zeta_bsplit(a, m, n, s, true, wp);
    // B = (B*D + A*C) * Q2  +  R.B * Q3   [L and R as in FLINT]
    l.b =
        l.b.mul(&r.d, wp, NR)
            .add(&ExactNum::mul_ref(&l.a, &r.c, wp, NR), wp, NR);
    l.b = l.b.mul(&r.q2, wp, NR).add(&r.b.mul(&l.q3, wp, NR), wp, NR);
    l.a = l.a.mul(&r.q3, wp, NR).add(&r.a.mul(&l.q3, wp, NR), wp, NR);
    l.c = l.c.mul(&r.d, wp, NR).add(&r.c.mul(&l.q1, wp, NR), wp, NR);
    if cont {
        l.d = l.d.mul(&r.d, wp, NR);
        l.q2 = l.q2.mul(&r.q2, wp, NR);
    }
    l.q1 = l.q1.mul(&r.q1, wp, NR);
    l.q3 = l.q3.mul(&r.q3, wp, NR);
    l
}

/// `arb_zeta_ui_borwein_bsplit`.
fn zeta_ui_borwein_bsplit(s: u64, p: usize) -> ExactNum {
    if s == 0 {
        return ExactNum::from_u8(1, p).ldexp(-1, p, NR).neg();
    }
    if s == 1 {
        return ExactNum::nan(Some(crate::Error::InvalidArgument));
    }
    let n = (p as u64).saturating_mul(BORWEIN_N_NUM) / BORWEIN_N_DEN + 2;
    let wp = p.saturating_add(32);
    let sum = zeta_bsplit(0, n + 1, n, s, false, wp);
    // eta = (A*C - B*Q1) / (Q3*C)
    let num = sum
        .a
        .mul(&sum.c, wp, NR)
        .sub(&sum.b.mul(&sum.q1, wp, NR), wp, NR);
    let den = sum.q3.mul(&sum.c, wp, NR);
    if den.is_zero() || den.is_nan() || den.is_inf() {
        return ExactNum::nan(Some(crate::Error::InvalidArgument));
    }
    let eta = num.div(&den, wp, NR);
    // zeta = eta * 2^{s-1} / (2^{s-1} - 1)
    let two_sm1 = ExactNum::from_u8(1, wp).ldexp((s - 1) as i64, wp, NR);
    let den_z = ExactNum::sub_ref(&two_sm1, &ExactNum::from_u8(1, wp), wp, NR);
    if den_z.is_zero() || den_z.is_nan() || den_z.is_inf() {
        return ExactNum::nan(Some(crate::Error::InvalidArgument));
    }
    let mut z = eta.mul(&two_sm1, wp, NR).div(&den_z, wp, NR);
    let _ = z.set_precision(p, RoundingMode::None);
    z
}
