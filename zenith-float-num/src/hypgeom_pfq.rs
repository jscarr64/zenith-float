//! Arb `{}_0F_1` via `acb_hypgeom_pfq_sum_bs` / `pfq_choose_n` (FLINT 3.6.0).
//!
//! `{}_0F_1(;a;z)=\sum_{k\ge0} z^k/((a)_k k!)` is `acb_hypgeom_pfq_direct` with
//! no upper parameters and lower parameters `(a,1)`. Binary splitting is the
//! method Arb uses for this series at dest width once `n>4` and `prec>=128`
//! (and, for exact low-bit `z`, in preference to rectangular splitting).

use crate::common::fixlog::ONE;
use crate::ExactComplex;
use crate::ExactNum;
use crate::RoundingMode;
use core::ops::Neg;

fn one(p: usize) -> ExactComplex {
    ExactComplex::one(p)
}

fn ui(k: usize, p: usize) -> ExactComplex {
    ExactComplex::from_real(ExactNum::from_u32(k as u32, p), p)
}

fn add_ui(x: &ExactComplex, k: usize, p: usize) -> ExactComplex {
    x.clone().add(&ui(k, p), p, RoundingMode::None)
}

/// `log2(|z|)` rounded up to an integer, in `Q32`: `|z|` in `[2^(e-1), 2^e)` at 64 bits gives `e`. Overestimating by at
/// most one bit only adds a few terms. `0` for a non-finite `z`; the caller handles `z = 0`.
fn log2_mag(z: &ExactComplex) -> i128 {
    match z.abs(64, RoundingMode::None).exponent() {
        Some(e) => i128::from(e).saturating_mul(ONE),
        None => 0,
    }
}

/// Arb `required_decrease = 0.01`, in `Q32` (floor).
const REQUIRED_DECREASE: i128 = 42_949_672;

/// Arb `acb_hypgeom_pfq_choose_n` specialised to `{}_0F_1(;a;z)`
/// (upper none, lower `(a,1)`), without Arb's `n_max = 50 + 10 prec`: the walk runs until the terms have
/// fallen 2^-(prec+4) below their peak. Arb's Bessel I/J switch to an asymptotic expansion for large |z| (and past
/// `n_max` its `pfq_direct` tail bound is infinite, never a silent truncation); here the `I_ν` series is used
/// for every Re z > 0, and a capped walk returned `n_best = 1`, i.e. only the first term (I_ν(20000), 7b972c3).
/// `None` only if the index would pass `u32::MAX` (the terms' index type).
///
/// The term-size walk is in `Q32` fixed point (`crate::common::fixlog`), so the choice is the same on every
/// target. The parts of `a` enter as signed powers of two ([`rough_mid`]), as the earlier `double` walk had them;
/// Arb uses the `double` midpoint itself. It only sets the number of terms.
fn choose_n_0f1(a: &ExactComplex, z: &ExactComplex, prec: usize) -> Option<usize> {
    if z.re().is_zero() && z.im().is_zero() {
        return Some(1);
    }
    let log2_z = log2_mag(z);
    let rm = RoundingMode::None;
    let are = rough_mid(a.re());
    let aim2 = rough_mid(a.im()).sqr(64, rm);
    let n_max = u32::MAX as usize;
    let n_skip = 1usize;
    let n_min = 1usize;
    let target = (prec as i128).saturating_add(4).saturating_mul(ONE);
    let mut term: i128 = 0;
    let mut term_max: i128 = 0;
    let mut accuracy_best: i128 = 0;
    let mut n_best = n_skip;
    let mut n = n_skip;
    while n < n_max {
        // t = 1 / (|a+n-1|^2 |n|^2); increase = 0.5 log2(t) + log2|z|
        let br = ExactNum::from_word((n - 1) as crate::Word, 64).add(&are, 64, rm);
        let nn = ExactNum::from_word(n as crate::Word, 64);
        let den = br
            .sqr(64, rm)
            .add(&aim2, 64, rm)
            .mul(&nn.sqr(64, rm), 64, rm);
        let increase = match crate::common::fixlog::log2_abs_down(&den) {
            // Arb: den > 1e-100
            Some(l) if l > crate::common::fixlog::LOG2_1E_M100_DOWN => {
                log2_z.saturating_sub(l >> 1)
            }
            _ => 1_000_000 * ONE,
        };
        term = term.saturating_add(increase);
        if term > term_max {
            term_max = term;
        }
        let accuracy = term_max - term;
        if accuracy > accuracy_best && n >= n_min && increase < -REQUIRED_DECREASE {
            n_best = n;
            accuracy_best = accuracy;
        }
        if accuracy_best > target {
            return Some(n_best);
        }
        n += 1;
    }
    None
}

/// A rough midpoint for [`choose_n_0f1`]: `sign(x) 2^floor(log2|x|)` for `floor(log2|x|)` in `[-1022, 63]`,
/// `0` for zero, and `1` otherwise (a huge, tiny or non-finite part), the values the earlier `double` walk used.
fn rough_mid(x: &ExactNum) -> ExactNum {
    let one = ExactNum::from_word(1, 64);
    if x.is_zero() {
        return ExactNum::new(64);
    }
    let Some(e) = x.leading_bit_exponent() else {
        return one;
    };
    if !(-1022..=63).contains(&e) {
        return one;
    }
    let m = one.ldexp(e, 64, RoundingMode::None);
    if x.is_negative() {
        m.neg()
    } else {
        m
    }
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
            t.add(&a1.mul(&b2, p, RoundingMode::None), p, RoundingMode::None)
        };
        let a_out = a1.mul(&a2, p, RoundingMode::None);
        let c_out = c1.mul(&c2, p, RoundingMode::None);
        (a_out, b_out, c_out)
    }
}

/// Arb `acb_hypgeom_pfq_sum_forward` for `{}_0F_1`.
/// Returns `(s, last_term)`.
fn sum_forward(
    a: &ExactComplex,
    z: &ExactComplex,
    n: usize,
    p: usize,
) -> (ExactComplex, ExactComplex) {
    let mut term = one(p);
    let mut s = one(p);
    let mut k = 0usize;
    while k + 1 < n {
        // term *= z / ((a+k)(k+1))
        let den = factor_lower(a, k, p);
        term = term
            .mul(z, p, RoundingMode::None)
            .div(&den, p, RoundingMode::None);
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
    let Some(n) = choose_n_0f1(a, z, p) else {
        return crate::complex_special::nan_pair(crate::Error::PrecisionRetryExhausted);
    };
    if n < 4 || p < 128 {
        return sum_forward(a, z, n.max(1), p).0;
    }
    // bsplit to n-1, then split off the last factor (Arb pfq_sum_bs.c).
    // Arb computes the last term after the n-1 bsplit only to attach a tail
    // bound in pfq_direct. Ziv on the Bessel entry already covers the remainder;
    // the sum itself is (A+B)/C over 0..n-1, same as pfq_sum_bs.
    let (u, v, w) = bsplit(a, z, 0, n - 1, p);
    u.add(&v, p, RoundingMode::None)
        .div(&w, p, RoundingMode::None)
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
            (_, Some(_)) if b.re().is_zero() && b.im().is_zero() => {
                ad.exponent().is_some_and(|e| e + slack < 0)
            }
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
            let half = ExactComplex::from_real(ExactNum::from_u8(1, p), p).div(
                &two,
                p,
                RoundingMode::ToEven,
            );
            let w = ExactComplex::new(ExactNum::new(p), half.re().clone().neg()); // -i/2
            let n = choose_n_0f1(&two, &w, p).ok_or(Error::PrecisionRetryExhausted)?;
            let fwd = sum_forward(&two, &w, n.max(1), p).0;
            // force the bsplit path even at p=64 by calling bsplit+close for p>=4
            let bs = if n >= 4 {
                let (u, v, ww) = bsplit(&two, &w, 0, n - 1, p);
                u.add(&v, p, RoundingMode::None)
                    .div(&ww, p, RoundingMode::None)
            } else {
                fwd.clone()
            };
            assert!(!bs.is_nan() && !fwd.is_nan(), "p={p} nan");
            assert!(near(&bs, &fwd, p, 8), "p={p} bsplit vs forward");
        }
        let _ = Consts::new()?;
        Ok(())
    }

    /// `{}_0F_1(;3/2;x^2/4) = \sinh(x)/x` at x = 20000 (I_{1/2}(20000), complex_bessel_mpmath): the terms
    /// grow until n ~ x/2 = 10^4, past Arb's `n_max = 50 + 10 prec`. Capping the walk there summed only
    /// the first term (7b972c3); the walk now runs until the terms have decreased by 2^-(prec+4).
    #[test]
    fn zf_0f1_large_z_runs_to_convergence() -> Result<(), Error> {
        let p = 128;
        let mut cc = Consts::new()?;
        let x = ExactNum::from_u32(20_000, p);
        let a = ExactComplex::from_real(
            ExactNum::from_u8(3, p).div(&ExactNum::from_u8(2, p), p, RoundingMode::None),
            p,
        );
        let w = ExactComplex::from_real(ExactNum::from_u32(100_000_000, p), p);
        let got = hypgeom_0f1(&a, &w, p + 64);
        let want = ExactComplex::from_real(
            x.sinh(p + 64, RoundingMode::None, &mut cc)
                .div(&x, p + 64, RoundingMode::None),
            p + 64,
        );
        assert!(
            near(&got, &want, p, 8),
            "0F1(;3/2;1e8): got {:?}, want {:?}",
            got.re(),
            want.re()
        );
        Ok(())
    }
}
