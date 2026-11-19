//! Arb `{}_0F_1` via `acb_hypgeom_pfq_sum_bs` / `pfq_choose_n` (FLINT 3.6.0).
//!
//! `{}_0F_1(;a;z)=\sum_{k\ge0} z^k/((a)_k k!)` is `acb_hypgeom_pfq_direct` with
//! no upper parameters and lower parameters `(a,1)`. Binary splitting is the
//! method Arb uses for this series at dest width once `n>4` and `prec>=128`
//! (and, for exact low-bit `z`, in preference to rectangular splitting).

use crate::ExactComplex;
use crate::ExactNum;
use crate::RoundingMode;

fn one(p: usize) -> ExactComplex {
    ExactComplex::one(p)
}

fn ui(k: usize, p: usize) -> ExactComplex {
    ExactComplex::from_real(ExactNum::from_u32(k as u32, p), p)
}

fn add_ui(x: &ExactComplex, k: usize, p: usize) -> ExactComplex {
    x.clone().add(&ui(k, p), p, RoundingMode::None)
}

/// `log2(x)` for `x>0`, used only by the Arb term-count heuristic.
fn log2_f64(x: f64) -> f64 {
    if !(x > 0.0) {
        return f64::NEG_INFINITY;
    }
    let bits = x.to_bits();
    let exp_bits = (bits >> 52) & 0x7ff;
    if exp_bits == 0 {
        // subnormal: scale and retry
        return log2_f64(x * (1u64 << 52) as f64) - 52.0;
    }
    if exp_bits == 0x7ff {
        return x;
    }
    let e = exp_bits as i32 - 1023;
    let m = f64::from_bits((bits & ((1u64 << 52) - 1)) | (0x3ffu64 << 52));
    let y = m - 1.0;
    // log(1+y) = y * (1 - y/2 + y^2/3 - y^3/4 + y^4/5 - y^5/6 + y^6/7)
    let y2 = y * y;
    let lny = y * (1.0 - y * 0.5 + y2 * (1.0 / 3.0 - y * 0.25 + y2 * (0.2 - y / 6.0 + y2 / 7.0)));
    lny * 1.4426950408889634074 + f64::from(e)
}

/// Upper bound on `log2(|z|)` from the leading bit and the `0.m · 2^e` exponent.
fn log2_mag(z: &ExactComplex) -> f64 {
    let a = z.abs(64, RoundingMode::None);
    if a.is_zero() {
        return f64::NEG_INFINITY;
    }
    // value is in [2^(e-1), 2^e); e = exponent(). Using e overestimates |z|
    // by at most 1 bit, which only adds a few terms.
    match a.exponent() {
        Some(e) => e as f64,
        None => 0.0,
    }
}

/// Arb `acb_hypgeom_pfq_choose_n` specialised to `{}_0F_1(;a;z)`
/// (upper none, lower `(a,1)`). `n_max = 50 + 10 prec`.
fn choose_n_0f1(a: &ExactComplex, z: &ExactComplex, prec: usize) -> usize {
    if z.re().is_zero() && z.im().is_zero() {
        return 1;
    }
    let log2_z = log2_mag(z);
    // midpoints of a at 64 bits, as f64 (Arb uses arf_get_d of the mid)
    let a64 = a.abs(64, RoundingMode::None);
    let are = signed_f64_part(a.re());
    let aim = signed_f64_part(a.im());
    let _ = a64;
    let n_max = (50 + 10 * prec).min(usize::MAX / 2);
    let n_skip = 1usize;
    let n_min = 1usize;
    let required_decrease = 0.01;
    let mut term = 0.0;
    let mut term_max = 0.0;
    let mut accuracy_best = 0.0;
    let mut n_best = n_skip;
    let mut success = false;
    let mut n = n_skip;
    while n < n_max {
        // t = 1 / (|a+n-1|^2 |n|^2); increase = 0.5 log2(t) + log2|z|
        let br = are + (n as f64) - 1.0;
        let b1 = n as f64; // (1+n-1)
        let den = (br * br + aim * aim) * (b1 * b1);
        let increase = if den > 1e-100 {
            -0.5 * log2_f64(den) + log2_z
        } else {
            1.0e6
        };
        term += increase;
        if term > term_max {
            term_max = term;
        }
        let accuracy = term_max - term;
        if accuracy > accuracy_best && n >= n_min && increase < -required_decrease {
            n_best = n;
            accuracy_best = accuracy;
        }
        if accuracy_best > (prec as f64) + 4.0 {
            success = true;
            break;
        }
        n += 1;
    }
    if success {
        n_best
    } else {
        n_best.max(n_min).min(n_max)
    }
}

fn signed_f64_part(x: &ExactNum) -> f64 {
    if x.is_zero() {
        return 0.0;
    }
    let Some(e) = x.leading_bit_exponent() else {
        return 0.0;
    };
    // floor(log2(|x|)) = e, so |x| is in [2^e, 2^{e+1}). Use 2^e with the sign;
    // choose_n only needs a rough mid (Arb uses a double of the mid).
    let mag = if (0..63).contains(&e) {
        (1u64 << e) as f64
    } else if e == 63 {
        8.0 * ((1u64 << 60) as f64)
    } else if (-1022..0).contains(&e) {
        f64::from_bits(((e + 1023) as u64) << 52)
    } else {
        1.0
    };
    if x.is_negative() { -mag } else { mag }
}

/// Leaf / product factor: `z` (no upper parameters) or 1 when `z` is omitted.
fn factor_upper(z: Option<&ExactComplex>, p: usize) -> ExactComplex {
    match z {
        Some(z) => z.clone(),
        None => one(p),
    }
}

/// `C(k) = (a+k)(1+k)`.
fn factor_lower(a: &ExactComplex, k: usize, p: usize) -> ExactComplex {
    add_ui(a, k, p).mul(&ui(k + 1, p), p, RoundingMode::None)
}

/// Arb `bsplit` for `{}_0F_1`, `invz = 0`.
///
/// Returns `(A, B, C)` over `k ∈ [aa, bb)`. For a length-1 interval
/// `B` is unused and treated as `C` by the caller, as in `pfq_sum_bs.c`.
fn bsplit(
    a: &ExactComplex,
    z: &ExactComplex,
    aa: usize,
    bb: usize,
    p: usize,
) -> (ExactComplex, ExactComplex, ExactComplex) {
    if bb - aa == 1 {
        let a1 = factor_upper(Some(z), p);
        let c1 = factor_lower(a, aa, p);
        // B unused at a leaf (Arb skips `acb_set(B1, C1)`).
        (a1, one(p), c1)
    } else {
        let m = aa + (bb - aa) / 2;
        let (a1, b1, c1) = bsplit(a, z, aa, m, p);
        let (a2, b2, c2) = bsplit(a, z, m, bb, p);
        let b_out = if bb - m == 1 {
            // B2 = C2
            let sum = if m - aa == 1 {
                a1.clone().add(&c1, p, RoundingMode::None)
            } else {
                a1.clone().add(&b1, p, RoundingMode::None)
            };
            sum.mul(&c2, p, RoundingMode::None)
        } else {
            let t = if m - aa == 1 {
                c1.mul(&c2, p, RoundingMode::None)
            } else {
                b1.mul(&c2, p, RoundingMode::None)
            };
            t.add(&a1.clone().mul(&b2, p, RoundingMode::None), p, RoundingMode::None)
        };
        let a_out = a1.mul(&a2, p, RoundingMode::None);
        let c_out = c1.mul(&c2, p, RoundingMode::None);
        (a_out, b_out, c_out)
    }
}

/// Arb `acb_hypgeom_pfq_sum_forward` for `{}_0F_1`.
/// Returns `(s, last_term)`.
fn sum_forward(a: &ExactComplex, z: &ExactComplex, n: usize, p: usize) -> (ExactComplex, ExactComplex) {
    let mut term = one(p);
    let mut s = one(p);
    let mut k = 0usize;
    while k + 1 < n {
        // term *= z / ((a+k)(k+1))
        let den = factor_lower(a, k, p);
        term = term.mul(z, p, RoundingMode::None).div(&den, p, RoundingMode::None);
        s = s.add(&term, p, RoundingMode::None);
        k += 1;
    }
    (s, term)
}

/// `{}_0F_1(;a;z)` at working precision `p`.
pub(crate) fn hypgeom_0f1(a: &ExactComplex, z: &ExactComplex, p: usize) -> ExactComplex {
    if z.re().is_zero() && z.im().is_zero() {
        return one(p);
    }
    let n = choose_n_0f1(a, z, p);
    if n < 4 || p < 128 {
        return sum_forward(a, z, n.max(1), p).0;
    }
    // bsplit to n-1, then split off the last factor (Arb pfq_sum_bs.c).
    // Arb computes the last term after the n-1 bsplit only to attach a tail
    // bound in pfq_direct. Ziv on the Bessel entry already covers the remainder;
    // the sum itself is (A+B)/C over 0..n-1, same as pfq_sum_bs.
    let (u, v, w) = bsplit(a, z, 0, n - 1, p);
    u.add(&v, p, RoundingMode::None).div(&w, p, RoundingMode::None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Consts;
    use crate::Error;
    use core::ops::Neg;

    fn near(a: &ExactComplex, b: &ExactComplex, p: usize, slack: i64) -> bool {
        let d = a.clone().sub(b, p, RoundingMode::None);
        if d.is_nan() {
            return false;
        }
        let ad = d.abs(64, RoundingMode::None);
        if ad.is_zero() {
            return true;
        }
        let ab = b.abs(64, RoundingMode::None);
        match (ad.exponent(), ab.exponent()) {
            (_, Some(_)) if b.re().is_zero() && b.im().is_zero() => ad.exponent().map_or(false, |e| e + slack < 0),
            (Some(ed), Some(eb)) => ed + (p as i64) < eb + slack,
            _ => false,
        }
    }

    /// `{}_0F_1(;2;-i/2)` is `J_1(1+i) / ((1+i)/2)`, checked against the
    /// forward sum at dest 64/256/1024.
    #[test]
    fn zf_0f1_bsplit_matches_forward() -> Result<(), Error> {
        for p in [64usize, 256, 1024, 8192] {
            let two = ExactComplex::from_real(ExactNum::from_u8(2, p), p);
            let half = ExactComplex::from_real(ExactNum::from_u8(1, p), p)
                .div(&two, p, RoundingMode::ToEven);
            let w = ExactComplex::new(ExactNum::new(p), half.re().clone().neg()); // -i/2
            let n = choose_n_0f1(&two, &w, p);
            let fwd = sum_forward(&two, &w, n.max(1), p).0;
            // force the bsplit path even at p=64 by calling bsplit+close for p>=4
            let bs = if n >= 4 {
                let (u, v, ww) = bsplit(&two, &w, 0, n - 1, p);
                u.add(&v, p, RoundingMode::None).div(&ww, p, RoundingMode::None)
            } else {
                fwd.clone()
            };
            assert!(!bs.is_nan() && !fwd.is_nan(), "p={p} nan");
            assert!(near(&bs, &fwd, p, 8), "p={p} bsplit vs forward");
        }
        let _ = Consts::new()?;
        Ok(())
    }
}
