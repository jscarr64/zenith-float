//! Arb/Flint variable-order Carlson Taylor sums for complex \(R_F\) and \(R_J\).
//!
//! Faithful port of `acb_elliptic_rf_taylor_sum` and
//! `acb_elliptic_rj_taylor_sum` (FLINT 3.6.0), except that the unused
//! trailing `d` update in the R_F `x`-loop is skipped: that step is not
//! exact, and FLINT's `fmpz_divexact_ui` ignores the remainder. Used so the duplication
//! loops can stop when \((\mathrm{rel})^{order} < 2^{-p}\) with
//! \(order \sim p^{1/3}\) to \(p^{2/5}\), instead of waiting for a
//! deviation of \(2^{-(p/3+16)}\) that forces \(\sim p/6\) iterations.

use super::*;
use crate::complex_special::nan_pair;
use crate::Error;
use crate::ExactInt;
use alloc::vec;
use alloc::vec::Vec;

const RF_DEN_RATIO: [u16; 512] = [
    1, 1, 10, 7, 12, 11, 26, 1, 136, 19, 2, 23, 20, 1, 58, 31, 16, 1, 74, 1, 164, 43, 2, 47, 56, 1,
    106, 1, 4, 59, 122, 1, 32, 67, 2, 71, 292, 1, 2, 79, 24, 83, 2, 1, 356, 1, 2, 1, 1552, 1, 202,
    103, 4, 107, 218, 1, 904, 1, 2, 1, 44, 1, 10, 127, 64, 131, 2, 1, 548, 139, 2, 1, 8, 1, 298,
    151, 4, 1, 314, 1, 16, 163, 2, 167, 52, 1, 346, 1, 8, 179, 362, 1, 4, 1, 2, 191, 6176, 1, 394,
    199, 4, 1, 2, 1, 8, 211, 2, 1, 4, 1, 2, 223, 16, 227, 458, 1, 932, 1, 2, 239, 1928, 1, 2, 1, 4,
    251, 2, 1, 32896, 1, 2, 263, 4, 1, 538, 271, 8, 1, 554, 1, 1124, 283, 2, 1, 272, 1, 586, 1, 4,
    1, 2, 1, 8, 307, 2, 311, 1252, 1, 634, 1, 32, 1, 2, 1, 4, 331, 2, 1, 2696, 1, 2, 7, 4, 347,
    698, 1, 5648, 1, 2, 359, 76, 1, 2, 367, 8, 1, 746, 1, 4, 379, 2, 383, 64, 1, 778, 1, 4, 1, 794,
    1, 3208, 1, 2, 1, 1636, 1, 2, 1, 16, 419, 842, 1, 4, 1, 2, 431, 3464, 1, 2, 439, 4, 443, 2, 1,
    14368, 1, 2, 1, 1828, 1, 922, 463, 8, 467, 2, 1, 4, 1, 2, 479, 16, 1, 2, 487, 4, 491, 2, 1, 8,
    499, 2, 503, 4, 1, 1018, 1, 256, 1, 2, 1, 2084, 523, 2, 1, 184, 1, 2, 1, 4, 1, 1082, 1, 16,
    547, 2, 1, 4, 1, 1114, 1, 8, 563, 2, 1, 2276, 571, 2, 1, 18464, 1, 2, 1, 4, 587, 2, 1, 4744, 1,
    2, 599, 2404, 1, 2, 607, 16, 1, 1226, 1, 2468, 619, 2, 1, 40, 1, 2, 631, 4, 1, 2, 1, 41024,
    643, 2, 647, 4, 1, 1306, 1, 8, 659, 1322, 1, 4, 1, 2, 1, 10768, 1, 1354, 1, 4, 683, 2, 1, 8,
    691, 2, 1, 4, 1, 1402, 1, 32, 1, 1418, 1, 4, 1, 2, 719, 8, 1, 2, 727, 12, 1, 1466, 1, 16, 739,
    2, 743, 4, 1, 2, 751, 8, 1, 1514, 1, 3044, 1, 2, 2, 49216, 1, 1546, 1, 4, 1, 2, 1, 8, 787, 2,
    1, 4, 1, 1594, 1, 16, 1, 2, 1, 3236, 811, 2, 1, 8, 1, 1642, 823, 4, 827, 1658, 1, 32, 1, 2,
    839, 116, 1, 2, 1, 8, 1, 1706, 1, 3428, 859, 2, 863, 16, 1, 2, 1, 4, 1, 1754, 1, 7048, 883, 2,
    887, 4, 1, 2, 1, 64, 1, 2, 1, 4, 907, 2, 911, 8, 1, 2, 919, 4, 1, 2, 1, 14864, 1, 2, 1, 3748,
    1, 1882, 1, 8, 947, 2, 1, 3812, 1, 2, 1, 992, 1, 2, 967, 4, 971, 2, 1, 7816, 1, 2, 983, 4, 1,
    2, 991, 16, 1, 1994, 1, 4, 1, 2, 1, 8072, 1, 2026, 1, 4, 1019, 2042, 1,
];

const RJ_DEN_RATIO: [u16; 512] = [
    1, 1, 14, 3, 44, 13, 10, 17, 152, 1, 46, 1, 12, 29, 62, 1, 16, 37, 2, 41, 172, 1, 94, 7, 8, 53,
    2, 1, 236, 61, 2, 1, 2144, 1, 142, 73, 20, 1, 158, 3, 664, 1, 2, 89, 4, 1, 2, 97, 16, 101, 206,
    1, 428, 109, 2, 113, 8, 1, 2, 11, 4, 1, 254, 1, 8384, 1, 2, 137, 556, 1, 2, 1, 8, 149, 302, 1,
    4, 157, 2, 1, 2608, 1, 334, 13, 4, 173, 2, 1, 1432, 181, 2, 1, 4, 1, 382, 193, 32, 197, 398, 1,
    4, 1, 2, 1, 1688, 1, 2, 1, 4, 1, 446, 1, 3632, 229, 2, 233, 4, 1, 478, 241, 24, 1, 2, 1, 1004,
    1, 2, 257, 128, 1, 526, 1, 4, 269, 542, 1, 8, 277, 2, 281, 1132, 1, 2, 17, 16, 293, 2, 1, 4, 1,
    2, 1, 2456, 1, 622, 313, 4, 317, 2, 1, 32, 1, 2, 1, 1324, 1, 2, 337, 8, 1, 14, 1, 1388, 349, 2,
    353, 16, 1, 718, 19, 4, 1, 734, 1, 8, 373, 10, 1, 1516, 1, 766, 1, 64, 389, 2, 1, 4, 397, 2,
    401, 8, 1, 2, 409, 4, 1, 2, 1, 6704, 421, 2, 1, 4, 1, 862, 433, 8, 1, 878, 1, 1772, 1, 2, 449,
    32, 1, 2, 457, 4, 461, 926, 1, 3736, 1, 2, 1, 4, 1, 958, 1, 16, 1, 974, 1, 1964, 1, 2, 1, 3992,
    1, 1006, 1, 4, 509, 2, 1, 256, 1, 2, 521, 2092, 1, 2, 23, 8, 1, 2, 1, 4, 541, 2, 1, 8752, 1, 2,
    1, 4, 557, 2, 1, 4504, 1, 2, 569, 2284, 1, 2, 577, 32, 1, 2, 1, 2348, 1, 2, 593, 8, 1, 1198,
    601, 4, 1, 1214, 1, 16, 613, 2, 617, 2476, 1, 2, 1, 8, 1, 1262, 1, 4, 1, 2, 641, 41152, 1,
    1294, 1, 4, 653, 2, 1, 5272, 661, 2, 1, 4, 1, 2, 673, 16, 677, 2, 1, 2732, 1, 2, 1, 5528, 1, 2,
    1, 4, 701, 2, 1, 32, 709, 2, 1, 4, 1, 1438, 1, 8, 1, 1454, 3, 4, 733, 2, 1, 11824, 1, 1486, 1,
    4, 1, 1502, 1, 8, 757, 2, 761, 4, 1, 2, 769, 128, 773, 2, 1, 4, 1, 2, 1, 6296, 1, 2, 1, 4, 797,
    2, 1, 16, 1, 2, 809, 3244, 1, 2, 1, 8, 821, 1646, 1, 3308, 829, 2, 1, 32, 1, 1678, 29, 4, 1, 2,
    1, 8, 853, 2, 857, 3436, 1, 1726, 1, 16, 1, 2, 1, 4, 877, 2, 881, 7064, 1, 1774, 1, 4, 1, 2, 1,
    64, 1, 2, 1, 3628, 1, 1822, 1, 8, 1, 1838, 1, 4, 1, 2, 929, 16, 1, 2, 937, 4, 941, 2, 1, 7576,
    1, 2, 953, 4, 1, 2, 31, 32, 1, 1934, 1, 3884, 1, 2, 977, 8, 1, 1966, 1, 4, 1, 1982, 1, 16, 997,
    2, 1, 4, 1, 2, 1009, 8, 1013, 2, 1, 4076, 1021, 2, 1,
];

fn int_from_u(n: u64) -> ExactInt {
    ExactInt::from_u64(n)
}

/// Exact `a / d`. `None` if `d` does not divide `a`.
/// Never returns a truncated (wrong) quotient. FLINT's
/// `fmpz_divexact_ui` ignores a nonzero remainder (truncating `idiv`
/// for small fmpz; `mpz_divexact_ui` is undefined). A failed exact
/// division is a checked error: the Taylor sum returns a typed NaN
/// (`InvalidArgument`) instead of a wrong coefficient.
fn try_divexact_u(a: &ExactInt, d: u64) -> Option<ExactInt> {
    // `div_rem` is `None` only for a zero divisor; that is not an exact division either.
    let (q, r) = a.div_rem(&int_from_u(d))?;
    if r.is_zero() {
        Some(q)
    } else {
        None
    }
}

fn mul_u(a: &ExactInt, d: u64) -> ExactInt {
    a.mul(&int_from_u(d))
}

fn universal_den(tab: &[u16; 512], nmax: usize) -> ExactInt {
    let mut den = ExactInt::one();
    for t in tab.iter().take(nmax + 1).skip(1) {
        den = mul_u(&den, u64::from(*t));
    }
    den
}

fn addmul_int(s: ExactComplex, term: &ExactComplex, c: &ExactInt, p: usize) -> ExactComplex {
    if c.is_zero() {
        return s;
    }
    let cf = ExactComplex::from_real(c.to_exact_num(p, rm()), p);
    s.add(&term.mul(&cf, p, rm()), p, rm())
}

// FLINT computes each order as a `double` expression (`2.1 * pow(prec, 0.4)` and so on) and truncates.
// Here it is the exact floor of the same real value, from integers only (`fixlog::floor_scaled_root`: the
// largest `n` with `(cd n)^b <= cn^b prec^a` for the factor `cn/cd` and the power `a/b`), so it is the same on
// every target and in `no_std`. The order only sets how many Taylor terms are summed; the truncation error is
// bounded separately.

/// `p^(a/b)` scaled by `cn/cd`, floored, capped at 500 (the top of every clamp below).
fn order_floor(cn: u64, cd: u64, prec: usize, a: u64, b: u64) -> usize {
    let p = u64::try_from(prec.max(2)).unwrap_or(u64::MAX);
    usize::try_from(crate::common::fixlog::floor_scaled_root(
        cn, cd, p, 1, a, b, 500,
    ))
    .unwrap_or(500)
}

/// Arb `order = 2.1 p^{0.4}` (real) or `2.5 p^{0.4}` (complex), clamped to `[2, 500]`.
pub(super) fn rf_order(prec: usize, real: bool) -> usize {
    let (cn, cd) = if real { (21, 10) } else { (25, 10) };
    order_floor(cn, cd, prec, 2, 5).clamp(2, 500)
}

/// Arb `order = 2.3 p^{0.34}` (real R_J) or `2.5 p^{0.35}` (complex), clamped to `[5, 500]`.
pub(super) fn rj_order(prec: usize, real: bool) -> usize {
    let o = if real {
        order_floor(23, 10, prec, 17, 50)
    } else {
        order_floor(25, 10, prec, 7, 20)
    };
    o.clamp(5, 500)
}

/// Arb RD-as-RJ: `2.0 p^{0.33}` (real) or `2.2 p^{0.33}` (complex), clamped to `[5, 500]`.
pub(super) fn rd_order(prec: usize, real: bool) -> usize {
    let (cn, cd) = if real { (20, 10) } else { (22, 10) };
    order_floor(cn, cd, prec, 33, 100).clamp(5, 500)
}

/// FLINT 3.6.0's truncation bound for the Carlson Taylor sums (`acb_elliptic/rf.c:285-295`, `rj.c:466-477`).
/// With `eps = max(|X|, |Y|, |Z|)` (and `|P|` for R_J) at the stop, Arb bounds every coefficient of total degree
/// `n` by `2 (9/8)^n`, so the terms of degree `>= order` that the sum leaves out total at most
/// `2 sum_{n >= order} (9 eps / 8)^n = 2 q^order / (1 - q)` with `q = 9 eps / 8` (`mag_geom_series`; infinite
/// for `q >= 1`). FLINT adds that to the result's radius. zf has no radius here, so it keeps duplicating until
/// the bound is below `2^-prec` (the series is `1 + O(eps^2)`, so this is relative), and only then sums.
/// Every rounding is upward: `q` at 64 bits, then `Q32` `log2` bounds (`fixlog`).
pub(super) fn trunc_bound_below_prec(eps: &ExactNum, order: usize, prec: usize) -> bool {
    use crate::common::fixlog;
    if eps.is_nan() || eps.is_inf() {
        return false;
    }
    if eps.is_zero() {
        return true;
    }
    let q = ExactNum::mul_ref(eps, &ExactNum::from_u8(9, 64), 64, RoundingMode::Up).ldexp(
        -3,
        64,
        RoundingMode::Up,
    );
    let one = ExactNum::from_u8(1, 64);
    if !matches!(q.compare(&one), Some(c) if c < 0) {
        return false;
    }
    let om = one.sub(&q, 64, RoundingMode::Down);
    let (Some(lq), Some(lom)) = (fixlog::log2_abs_up(&q), fixlog::log2_abs_down(&om)) else {
        return false;
    };
    // log2(2 q^order / (1 - q)) <= 1 + order lq - lom
    let lb = fixlog::ONE
        .saturating_add((order as i128).saturating_mul(lq))
        .saturating_sub(lom);
    lb < (prec as i128).saturating_mul(fixlog::ONE).saturating_neg()
}

/// `rel^{order} < 2^{-prec}`, i.e. `dmax < scale * 2^{-ceil(prec/order)}`.
pub(super) fn rel_pow_below_prec(
    dmax: &ExactNum,
    scale: &ExactNum,
    order: usize,
    prec: usize,
    p: usize,
) -> bool {
    if scale.is_zero() {
        return dmax.is_zero();
    }
    let order = order.max(1);
    let shift = prec.div_ceil(order) as Exponent;
    let thresh = scale.ldexp_ref(-shift, p, rm());
    matches!(dmax.compare(&thresh), Some(c) if c < 0)
}

/// Arb `acb_elliptic_rf_taylor_sum`. Returns the series without the final \(A^{-1/2}\).
pub(super) fn rf_taylor_sum(
    e2: &ExactComplex,
    e3: &ExactComplex,
    nterms: usize,
    p: usize,
) -> ExactComplex {
    rf_taylor_sum_checked(e2, e3, nterms, p).unwrap_or_else(|| nan_pair(Error::InvalidArgument))
}

fn rf_taylor_sum_checked(
    e2: &ExactComplex,
    e3: &ExactComplex,
    nterms: usize,
    p: usize,
) -> Option<ExactComplex> {
    let nterms = nterms.clamp(1, 500);
    let nmax = nterms - 1;
    let ymax = nmax / 3;
    let xmax = nmax / 2;
    let den = universal_den(&RF_DEN_RATIO, nmax);
    // powers of E2
    let mut e2pow: Vec<ExactComplex> = Vec::with_capacity(xmax + 1);
    e2pow.push(c_u32(1, p));
    if xmax >= 1 {
        e2pow.push(e2.clone());
        for i in 2..=xmax {
            e2pow.push(e2pow[i - 1].mul(e2, p, rm()));
        }
    }
    let mut c = den.clone();
    for y in 0..ymax {
        c = mul_u(&c, (2 * y + 1) as u64);
        c = try_divexact_u(&c, (2 * y + 2) as u64)?;
    }
    let mut res = ExactComplex::zero(p);
    for y in (0..=ymax).rev() {
        if y != ymax {
            c = mul_u(&c, (2 * y + 2) as u64);
            c = try_divexact_u(&c, (2 * y + 1) as u64)?;
        }
        let mut d = c.clone();
        let mut s = ExactComplex::zero(p);
        for (x, e2x) in e2pow.iter().enumerate() {
            let n = 2 * x + 3 * y;
            if n <= nmax {
                let mut e = try_divexact_u(&d, (2 * n + 1) as u64)?;
                if x % 2 == 1 {
                    e = e.neg();
                }
                if x != 0 || y != 0 {
                    s = addmul_int(s, e2x, &e, p);
                }
                // FLINT always updates `d` here, including after the last
                // used `x` for this `y`. That trailing step is unused and
                // is not exact (nmax=21: 317246823152966183226750 / 20
                // remainder 10). FLINT's `fmpz_divexact_ui` ignores the
                // remainder. The R_J port already skips the matching
                // trailing `c2`/`c3` updates (`m2 < m2start`). Skip the
                // unused step so every remaining `divexact` is exact for
                // nmax in 0..=499 (the clamp).
                let next_n = 2 * (x + 1) + 3 * y;
                if next_n <= nmax {
                    d = mul_u(&d, (2 * x + 2 * y + 1) as u64);
                    d = try_divexact_u(&d, (2 * x + 2) as u64)?;
                }
            }
        }
        res = res.mul(e3, p, rm()).add(&s, p, rm());
    }
    let den_c = ExactComplex::from_real(den.to_exact_num(p, rm()), p);
    Some(res.div(&den_c, p, rm()).add(&c_u32(1, p), p, rm()))
}

fn rj_fdiv(x: i32, y: i32) -> i32 {
    if x < 0 {
        -1
    } else {
        x / y
    }
}

/// Arb `acb_elliptic_rj_taylor_sum`. Returns the series without \(4^{-k} A^{-3/2}\).
pub(super) fn rj_taylor_sum(
    e2: &ExactComplex,
    e3: &ExactComplex,
    e4: &ExactComplex,
    e5: &ExactComplex,
    nterms: usize,
    p: usize,
) -> ExactComplex {
    rj_taylor_sum_checked(e2, e3, e4, e5, nterms, p)
        .unwrap_or_else(|| nan_pair(Error::InvalidArgument))
}

fn rj_taylor_sum_checked(
    e2: &ExactComplex,
    e3: &ExactComplex,
    e4: &ExactComplex,
    e5: &ExactComplex,
    nterms: usize,
    p: usize,
) -> Option<ExactComplex> {
    let nterms = nterms.clamp(1, 500);
    let nmax = nterms as i32 - 1;
    let m2dim = (nmax / 2 + 1) as usize;
    let m3dim = (nmax / 3 + 1) as usize;
    let den = universal_den(&RJ_DEN_RATIO, nterms - 1);
    // powtab[m3 * m2dim + m2] = E2^{m2} E3^{m3}
    let mut powtab: Vec<ExactComplex> = vec![c_u32(1, p); m2dim * m3dim];
    let m2max = (nmax / 2) as usize;
    for m2 in 0..=m2max {
        let m3cap = rj_fdiv(nmax - 2 * m2 as i32, 3);
        if m3cap < 0 {
            continue;
        }
        for m3 in 0..=m3cap as usize {
            let i = m3 * m2dim + m2;
            if m2 <= 1 && m3 <= 1 {
                powtab[i] = match (m2, m3) {
                    (0, 0) => c_u32(1, p),
                    (0, 1) => e3.clone(),
                    (1, 0) => e2.clone(),
                    _ => e2.mul(e3, p, rm()),
                };
            } else {
                let j = (m3 / 2) * m2dim + (m2 / 2);
                let k = (m3 - m3 / 2) * m2dim + (m2 - m2 / 2);
                powtab[i] = powtab[j].mul(&powtab[k], p, rm());
            }
        }
    }
    let mut c5 = mul_u(&den, 3);
    let m5start = (nmax / 5) as usize;
    for m5 in 0..m5start {
        c5 = mul_u(&c5, (2 * m5 + 1) as u64);
        c5 = try_divexact_u(&c5, (2 * m5 + 2) as u64)?;
    }
    let mut s5 = ExactComplex::zero(p);
    for m5 in (0..=m5start).rev() {
        let m4start = rj_fdiv(nmax - 5 * m5 as i32, 4);
        if m5 != m5start {
            c5 = mul_u(&c5, (2 * m5 + 2) as u64);
            c5 = try_divexact_u(&c5, (2 * m5 + 1) as u64)?;
        }
        let mut c4 = c5.clone();
        if m4start > 0 {
            for m4 in 0..m4start as usize {
                c4 = mul_u(&c4, (2 * m5 + 2 * m4 + 1) as u64);
                c4 = try_divexact_u(&c4, (2 * m4 + 2) as u64)?;
            }
        }
        let mut s4 = ExactComplex::zero(p);
        let m4_hi = if m4start < 0 { 0 } else { m4start as usize };
        // If m4start is -1, the inner loop should not run. Arb uses m4start = rj_fdiv which can be -1.
        if m4start >= 0 {
            for m4 in (0..=m4_hi).rev() {
                let m3start = rj_fdiv(nmax - 5 * m5 as i32 - 4 * m4 as i32, 3);
                if m4 != m4_hi {
                    c4 = mul_u(&c4, (2 * m4 + 2) as u64);
                    c4 = try_divexact_u(&c4, (2 * m5 + 2 * m4 + 1) as u64)?;
                }
                let mut c3 = c4.clone();
                let mut s3 = ExactComplex::zero(p);
                if m3start >= 0 {
                    for m3 in 0..=m3start as usize {
                        let m2start =
                            rj_fdiv(nmax - 5 * m5 as i32 - 4 * m4 as i32 - 3 * m3 as i32, 2);
                        let mut c2 = c3.clone();
                        if m2start >= 0 {
                            for m2 in 0..=m2start as usize {
                                let mm = (m5 + m4 + m3 + m2) as i32;
                                let nn =
                                    5 * m5 as i32 + 4 * m4 as i32 + 3 * m3 as i32 + 2 * m2 as i32;
                                let mut t = try_divexact_u(&c2, (2 * nn + 3) as u64)?;
                                if (mm + nn) % 2 != 0 {
                                    t = t.neg();
                                }
                                s3 = addmul_int(s3, &powtab[m3 * m2dim + m2], &t, p);
                                if (m2 as i32) < m2start {
                                    c2 = mul_u(&c2, (2 * m5 + 2 * m4 + 2 * m3 + 2 * m2 + 1) as u64);
                                    c2 = try_divexact_u(&c2, (2 * m2 + 2) as u64)?;
                                }
                            }
                        }
                        if (m3 as i32) < m3start {
                            c3 = mul_u(&c3, (2 * m5 + 2 * m4 + 2 * m3 + 1) as u64);
                            c3 = try_divexact_u(&c3, (2 * m3 + 2) as u64)?;
                        }
                    }
                }
                s4 = s4.mul(e4, p, rm()).add(&s3, p, rm());
            }
        }
        s5 = s5.mul(e5, p, rm()).add(&s4, p, rm());
    }
    let den_c = ExactComplex::from_real(den.to_exact_num(p, rm()), p);
    Some(s5.div(&den_c, p, rm()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_must::Must;

    /// nterms=22 (nmax=21) is the first order that traps in ZF-78 / ZF-66:
    /// after the last used x at y=1, FLINT still does
    /// `d *= (2x+2y+1); d /= (2x+2)` with dividend
    /// 317246823152966183226750 and divisor 20 (remainder 10). That
    /// trailing update is unused. The coefficient divisions are exact.
    #[test]
    fn rf_taylor_sum_at_nterms_22_does_not_trap() {
        let p = 64;
        let e2 = c_u32(1, p);
        let e3 = c_u32(1, p);
        let _ = rf_taylor_sum(&e2, &e3, 22, p);
    }

    /// Clamped RF orders must keep every used `divexact` exact.
    /// E2=E3=0: exactness of the integer recurrence does not depend on
    /// the deviations. 1..=40 covers nmax=21; larger values are sampled
    /// through the 500 clamp.
    #[test]
    fn rf_taylor_sum_all_clamped_orders_do_not_trap() {
        let p = 64;
        let e2 = ExactComplex::zero(p);
        let e3 = ExactComplex::zero(p);
        let orders = (1..=40).chain([50, 80, 100, 150, 200, 300, 400, 499, 500]);
        for nterms in orders {
            let _ = rf_taylor_sum(&e2, &e3, nterms, p);
        }
    }

    /// An inexact dividend must not become a truncated quotient.
    /// 15/4 is the trailing unused RF update at nmax=2.
    #[test]
    fn divexact_u_rejects_inexact_dividend() -> Result<(), String> {
        let a = int_from_u(15);
        assert!(
            try_divexact_u(&a, 4).is_none(),
            "15/4 must be refused, not truncated to 3"
        );
        assert_eq!(
            try_divexact_u(&int_from_u(20), 4).must("try_divexact_u")?,
            int_from_u(5)
        );
        Ok(())
    }
}
