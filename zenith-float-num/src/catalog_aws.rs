//! Anger \(\mathbf{J}_\nu\), Weber \(\mathbf{E}_\nu\), and Struve \(\mathbf{H}_\nu\)
//! kernels with rigorous tails (DLMF 11.2, 11.6, 11.10).
//!
//! Method selection is made from \((\nu,z)\) by a cost estimate before any
//! kernel is started. The 11.2.1 / 11.10.8–11 Maclaurin series is used only
//! when its term count (about \(|z|/2\)) plus cancellation bits is cheap.
//! Otherwise Anger/Weber use the 11.10.15–16 connection with \(A_\nu\). The
//! 11.6.1 expansion for \(\mathbf{H}_\nu-Y_\nu\) is used only when \(\nu\) is
//! real, \(z>0\), \(m\ge\nu-1/2\), and the first-omitted remainder (which does
//! not shrink past \(\sim e^{-z}\)) is small enough for the target precision.
//! If that bound is not small enough and Maclaurin is not cheap, \(K_\nu=
//! \mathbf{H}_\nu-Y_\nu\) is computed from the Laplace integral DLMF 11.5.2
//! with a proven tail. The amplitude \((1+t^2)^{\nu-1/2}\) is expanded on
//! panels inside its disk of analyticity (singularities at \(\pm i\)), and
//! \(e^{-z t}\) is integrated exactly; a Cauchy remainder on a disk of
//! radius \(9/10\) bounds the truncation.
//! Where \(|\nu|\gtrsim|z|\), 11.2.1 is used with cancellation guard bits;
//! the series is entire and has little cancellation when \(\nu\ge z\), so it
//! must not return `None` merely because the first working precision is
//! tight (the local `ziv_round` treats `None` as a hard NaN).
//!
//! \(A_\nu(z)=(1/\pi)\int_0^\infty\exp(-\nu t-z\sinh t)\,dt\) (DLMF 11.10.4,
//! \(\Re z>0\), \(\nu>-z\)) is computed, for large \(z\), from the Laplace /
//! Watson expansion (DLMF 11.11.4) stopped at the smallest term, with a
//! Cauchy remainder on the disk of analyticity of the integrand in the
//! \(u=\sinh t\) coordinate. When that bound is not tight enough, a monotone
//! integrand, a convexity tail bound, and Taylor panels sized to the decay
//! scale \(\sim 1/(\nu+z)\) are used instead. The catalog's uncontrolled
//! tanh-sinh path is not used.

use super::cat_nan;
use super::finish;
use super::finite;
use super::work_p;
use super::ziv_round;
use crate::defs::WORD_BIT_SIZE;
use crate::Consts;
use crate::ExactNum;
use crate::RoundingMode;
use crate::Sign;
use alloc::vec;
use alloc::vec::Vec;
use core::ops::Neg;

const N: RoundingMode = RoundingMode::None;
const UP: RoundingMode = RoundingMode::Up;
const DN: RoundingMode = RoundingMode::Down;

/// Series, jets, and Taylor panels stop on a remainder bound or a stall
/// (the bound stops improving, a term-index overflows, or a panel width
/// is already below \(2^{-(wp+16)}\)). There is no pre-set term or panel
/// refusal; a genuinely huge request may fail in the allocator.

fn bitlen(n: usize) -> i64 {
    (usize::BITS - n.leading_zeros()) as i64
}

fn exp_of(x: &ExactNum) -> i64 {
    x.exponent().unwrap_or(i64::MIN / 4)
}

fn one(p: usize) -> ExactNum {
    ExactNum::from_u8(1, p)
}

fn two(p: usize) -> ExactNum {
    ExactNum::from_u8(2, p)
}

fn half(p: usize) -> ExactNum {
    one(p).ldexp(-1, p, N)
}

/// C99 7.12.1 overflow: MAX when the mode rounds toward zero, else signed Inf.
fn overflow_signed(neg: bool, p: usize, rm: RoundingMode) -> ExactNum {
    let toward_zero = matches!(
        (!neg, rm),
        (_, RoundingMode::ToZero) | (true, RoundingMode::Down) | (false, RoundingMode::Up)
    );
    if toward_zero {
        let mut m = ExactNum::max_value(p);
        if neg {
            m.set_sign(Sign::Neg);
        }
        m
    } else if neg {
        crate::ext::INF_NEG
    } else {
        crate::ext::INF_POS
    }
}

fn signed_zero(neg: bool, p: usize) -> ExactNum {
    let mut z = ExactNum::new(p);
    if neg {
        z.set_sign(Sign::Neg);
    }
    z
}

fn nonpos_int(x: &ExactNum) -> bool {
    x.is_int() && (x.is_zero() || x.is_negative())
}

fn is_odd_int(n: &ExactNum, p: usize) -> bool {
    if !n.is_int() {
        return false;
    }
    !ExactNum::div_ref(n, &two(p), p, N).is_int()
}

/// \(1/\Gamma(z)\) is zero at a non-positive integer; otherwise \(\Gamma(z)\).
fn gamma_or_pole(z: &ExactNum, p: usize, cc: &mut Consts) -> Option<ExactNum> {
    if nonpos_int(z) {
        return None;
    }
    let g = z.gamma(p, N, cc);
    if !finite(&g) || g.is_zero() {
        None
    } else {
        Some(g)
    }
}

fn bits_good(wp: usize, peak: i64, result_exp: i64, lost: i64, tail_exp: i64) -> i64 {
    let cancel = (peak - result_exp).max(0);
    let arith = (wp as i64).saturating_sub(lost).saturating_sub(cancel);
    if tail_exp <= i64::MIN / 8 {
        arith
    } else {
        (result_exp - tail_exp - 1).min(arith)
    }
}

/// Negating a p-bit value swaps Up and Down; the other modes are symmetric.
fn mode_for_neg(rm: RoundingMode) -> RoundingMode {
    match rm {
        RoundingMode::Up => RoundingMode::Down,
        RoundingMode::Down => RoundingMode::Up,
        other => other,
    }
}

/// A conservative integer estimate of `|x|` when `|x| < 2^32`; `None` if huge.
fn abs_u32(x: &ExactNum) -> Option<u32> {
    if x.is_zero() {
        return Some(0);
    }
    let (_mant, _st, exp) = x.top_u64_mantissa_and_sticky()?;
    // |x| < 2^exp. If exp <= 0 then |x| < 1.
    if exp <= 0 {
        return Some(0);
    }
    if exp > 32 {
        return None;
    }
    // |x| >= 2^{exp-1}, so the integer ceiling of |x| is at most 2^exp
    // and at least 1. Use the mantissa to get a tight ceiling.
    let (mant, sticky, e) = x.top_u64_mantissa_and_sticky()?;
    // |x| = mant * 2^{e-64} plus a fraction if sticky.
    // For 1 <= e <= 32, shift mant down by (64-e).
    let shift = 64 - e;
    if !(1..=63).contains(&shift) {
        return None;
    }
    let mut v = mant >> shift;
    if sticky || (mant & ((1u64 << shift) - 1)) != 0 {
        v = v.saturating_add(1);
    }
    u32::try_from(v).ok()
}

/// Lower bound on the number of bits in `e^{-|z|}`: `|z| / ln 2 = |z| log2(e)`.
/// `log2(e) > 11/8 = 1.375`, so `|z| * 11/8` is a lower bound once `|z|` is known.
fn rem_bits_neg_exp(z: &ExactNum) -> i64 {
    let Some(a) = abs_u32(z) else {
        // |z| >= 2^32: e^{-z} has more than 2^32 bits. Saturate.
        return i64::MAX / 4;
    };
    if a == 0 {
        return 0;
    }
    // bits >= a * 11/8
    i64::from(a).saturating_mul(11) / 8
}

/// Term index at which the Maclaurin ratio `|(z/2)^2| / k^2` drops to 1/2:
/// `k ≈ |z| / sqrt(2)`.
fn maclaurin_k_est(z: &ExactNum) -> Option<u32> {
    let a = abs_u32(z)?;
    // ceil(a / sqrt(2)) <= ceil(a * 3 / 4) + 1, and a/sqrt(2) ≈ 0.707 a,
    // so ceil(a * 3 / 4) is a modest over-estimate (0.75 a).
    Some(a.saturating_mul(3).div_ceil(4).saturating_add(1).max(1))
}

/// Maclaurin is cheap only when the term count (about `|z|/2`) plus a
/// cancellation budget is a small amount of work. When `|ν| ≥ |z|` the
/// 11.2.1 series has little cancellation; otherwise cancellation can reach
/// about `|z|/ln 2` bits (Bessel-like).
fn series_cheap(z: &ExactNum, nu: &ExactNum, p: usize) -> bool {
    let Some(k) = maclaurin_k_est(z) else {
        return false;
    };
    // A hard cap: past a few dozen terms the connection / 11.5.2 paths
    // are faster than Maclaurin at large |z|.
    const K_MAX: u32 = 80;
    if k > K_MAX {
        return false;
    }
    let cancel = if nu.abs_cmp(&z.clone().abs()) != Some(-1) {
        8i64
    } else {
        rem_bits_neg_exp(z).min(1 << 20)
    };
    let work = i64::from(k).saturating_mul(1 + cancel.saturating_add(31) / 32);
    // Compare to a panel budget of about `p` panels: Maclaurin wins only
    // when it is clearly the cheaper option.
    work <= 256 && work <= (p as i64).saturating_add(64)
}

/// The 11.6.1 / Laplace remainder floor `~ e^{-z}` is small enough to
/// correctly round at `p` bits (the local Ziv loop needs about `p+64` good
/// bits, rounded down to a whole word).
fn exp_neg_z_enough(z: &ExactNum, p: usize) -> bool {
    rem_bits_neg_exp(z) >= (p as i64).saturating_add(64)
}

/// Overflow / exact-zero / finite closed form of \(\mathbf{H}_\nu(0)\).
fn struve_at_zero(
    z: &ExactNum,
    nu: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> ExactNum {
    let wrk = work_p(p).saturating_add(WORD_BIT_SIZE);
    let three_half = ExactNum::from_u8(3, wrk).div(&two(wrk), wrk, N);
    let nu_p_th = ExactNum::add_ref(nu, &three_half, wrk, N);
    // ν = -3/2, -5/2, ... : leading 11.2.1 coefficient vanishes, value is 0.
    if nonpos_int(&nu_p_th) {
        return finish(
            signed_zero(z.is_negative() && !is_odd_int(nu, wrk) && nu.is_int(), p),
            p,
            rm,
        );
    }
    // ν > -1 => (z/2)^{ν+1} → 0.
    let minus_one = one(wrk).neg();
    if nu.compare(&minus_one) == Some(1) {
        return finish(signed_zero(false, p), p, rm);
    }
    // ν = -1: H_{-1}(0) = 2/π.
    if nu.compare(&minus_one) == Some(0) {
        let pi = cc.pi(wrk, N);
        return finish(two(wrk).div(&pi, wrk, N), p, rm);
    }
    // Other ν < -1 with finite Γ(ν+3/2): signed Inf, sign of 1/(Γ(3/2)Γ(ν+3/2)), z→0+.
    // For z→0- the function is real only for integer ν; that sign follows the parity.
    if !nu.is_int() && z.is_negative() {
        return cat_nan();
    }
    let g1 = three_half.gamma(wrk, N, cc);
    let g2 = nu_p_th.gamma(wrk, N, cc);
    if !finite(&g1) || !finite(&g2) {
        return cat_nan();
    }
    let mut neg = g1.is_negative() != g2.is_negative();
    // (z/2)^{ν+1}: for z→0+ this is positive. For z→0- and integer ν,
    // H_n(-0) = (-1)^{n+1} H_n(+0).
    if z.is_negative() && nu.is_int() {
        // even n: H is odd, so flip; odd n: even function, keep.
        if !is_odd_int(nu, wrk) {
            neg = !neg;
        }
    }
    overflow_signed(neg, p, rm)
}

/// DLMF 11.2.1 Maclaurin series with a geometric tail and a cancellation guard.
fn struve_series(
    z: &ExactNum,
    nu: &ExactNum,
    wp: usize,
    cc: &mut Consts,
) -> Option<(ExactNum, i64)> {
    if z.is_zero() {
        return None;
    }
    let halfz = ExactNum::div_ref(z, &two(wp), wp, N);
    let z2 = halfz.sqr(wp, N);
    let three_half = ExactNum::from_u8(3, wp).div(&two(wp), wp, N);
    let mut sum = ExactNum::new(wp);
    let mut k = 0usize;
    // Fold the prefactor (z/2)^{ν+1} into the first term in log space so
    // (z/2)^{ν+1} / (Γ(3/2)Γ(ν+3/2)) does not overflow or underflow
    // intermediately when |ν| is large (H_{500.5}(100) is ~2^{-948}).
    let mut term = {
        let nu_p = ExactNum::add_ref(nu, &three_half, wp, N);
        if nonpos_int(&nu_p) {
            ExactNum::new(wp)
        } else if !halfz.is_positive() {
            return None;
        } else {
            let ln_h = halfz.ln(wp, N, cc);
            let ln_g1 = three_half.ln_gamma(wp, N, cc);
            let ln_g2 = nu_p.ln_gamma(wp, N, cc);
            if finite(&ln_h) && finite(&ln_g1) && finite(&ln_g2) {
                let ln = ExactNum::add_ref(nu, &one(wp), wp, N)
                    .mul(&ln_h, wp, N)
                    .sub(&ln_g1, wp, N)
                    .sub(&ln_g2, wp, N);
                let t = ln.exp(wp, N, cc);
                if finite(&t) {
                    t
                } else {
                    // Fall back to the factored form when the log path
                    // overflows; the factored form may still be representable
                    // after the division, or may correctly overflow.
                    let pref = halfz.pow(&ExactNum::add_ref(nu, &one(wp), wp, N), wp, N, cc);
                    let g1 = gamma_or_pole(&three_half, wp, cc)?;
                    let g2 = gamma_or_pole(&nu_p, wp, cc)?;
                    if !finite(&pref) || !finite(&g1) || !finite(&g2) || g2.is_zero() {
                        return None;
                    }
                    pref.div(&g1.mul(&g2, wp, N), wp, N)
                }
            } else {
                let pref = halfz.pow(&ExactNum::add_ref(nu, &one(wp), wp, N), wp, N, cc);
                let g1 = gamma_or_pole(&three_half, wp, cc)?;
                let g2 = gamma_or_pole(&nu_p, wp, cc)?;
                if !finite(&pref) || !finite(&g1) || !finite(&g2) || g2.is_zero() {
                    return None;
                }
                pref.div(&g1.mul(&g2, wp, N), wp, N)
            }
        }
    };
    let mut peak = exp_of(&term);
    // If the k=0 gamma in (ν+3/2) is a pole, term is 0; recurrence would divide
    // by zero, so each term is formed from the closed factor when needed.
    let mut formed = 0usize;
    loop {
        if k > 0 {
            // t_k / t_{k-1} = -(z/2)^2 / ((k+1/2)(k+ν+1/2))
            let half = half(wp);
            let kp = ExactNum::from_u64(k as u64, wp).add(&half, wp, N);
            let kpn = ExactNum::add_ref(&kp, nu, wp, N);
            if nonpos_int(&kp) || nonpos_int(&kpn) {
                term = ExactNum::new(wp);
            } else if term.is_zero() {
                // Rebuild the full 11.2.1 term, including (z/2)^{ν+1},
                // after a skipped pole: (z/2)^{2k+ν+1} (-1)^k /
                // (Γ(k+3/2) Γ(k+ν+3/2)).
                let a = ExactNum::from_u64(k as u64, wp).add(&three_half, wp, N);
                let b = ExactNum::add_ref(&a, nu, wp, N);
                let ga = match gamma_or_pole(&a, wp, cc) {
                    Some(g) => g,
                    None => {
                        k = k.checked_add(1)?;
                        continue;
                    }
                };
                let gb = match gamma_or_pole(&b, wp, cc) {
                    Some(g) => g,
                    None => {
                        k = k.checked_add(1)?;
                        continue;
                    }
                };
                let ln_h = halfz.ln(wp, N, cc);
                let ln_ga = a.ln_gamma(wp, N, cc);
                let ln_gb = b.ln_gamma(wp, N, cc);
                term = if finite(&ln_h) && finite(&ln_ga) && finite(&ln_gb) {
                    let two_k_nu1 = ExactNum::from_u64((k as u64).saturating_mul(2), wp)
                        .add(nu, wp, N)
                        .add(&one(wp), wp, N);
                    let ln = two_k_nu1
                        .mul(&ln_h, wp, N)
                        .sub(&ln_ga, wp, N)
                        .sub(&ln_gb, wp, N);
                    let t = ln.exp(wp, N, cc);
                    if !k.is_multiple_of(2) {
                        t.neg()
                    } else {
                        t
                    }
                } else if finite(&ga) && finite(&gb) && !gb.is_zero() {
                    let pow = halfz.pow(
                        &ExactNum::from_u64((k as u64).saturating_mul(2), wp)
                            .add(nu, wp, N)
                            .add(&one(wp), wp, N),
                        wp,
                        N,
                        cc,
                    );
                    let t = pow.div(&ga.mul(&gb, wp, N), wp, N);
                    if !k.is_multiple_of(2) {
                        t.neg()
                    } else {
                        t
                    }
                } else {
                    ExactNum::new(wp)
                };
            } else {
                let den = ExactNum::mul_ref(&kp, &kpn, wp, N);
                if den.is_zero() || !finite(&den) {
                    term = ExactNum::new(wp);
                } else {
                    term = term.mul(&z2, wp, N).div(&den, wp, N).neg();
                }
            }
        }
        if finite(&term) && !term.is_zero() {
            peak = peak.max(exp_of(&term));
            sum = sum.add(&term, wp, N);
            formed += 1;
            // Geometric tail once |t_{k+1}/t_k| ≤ 1/2, i.e.
            // |(z/2)^2| / |(k+3/2)(k+ν+3/2)| ≤ 1/2.
            let kp = ExactNum::from_u64(k as u64, wp).add(&three_half, wp, N);
            let kpn = ExactNum::add_ref(&kp, nu, wp, N);
            let den = ExactNum::mul_ref(&kp, &kpn, wp, N);
            if finite(&den) && !den.is_zero() {
                let ratio_up = ExactNum::div_ref(&z2, &den, wp, UP).abs();
                if ratio_up.compare(&half(wp)) != Some(1) {
                    // r ≤ 1/2: the remaining tail after this term is at most
                    // |term| r / (1-r) ≤ |term|. Stop once that is 2^{-(wp-lost)}
                    // below the sum; otherwise keep adding terms at this wp.
                    let one_m = ExactNum::sub_ref(&one(wp), &ratio_up, wp, DN);
                    let tail = if one_m.is_positive() {
                        ExactNum::mul_ref(&term, &ratio_up, wp, UP).div(&one_m, wp, UP)
                    } else {
                        term.clone()
                    };
                    let result = sum.clone();
                    if !finite(&result) {
                        return None;
                    }
                    let lost = 2 * bitlen(k.saturating_add(2)) + 6;
                    let tail_exp = exp_of(&tail);
                    let good = bits_good(wp, peak, exp_of(&result), lost, tail_exp);
                    let cancel = (peak - exp_of(&result)).max(0);
                    let budget = (wp as i64).saturating_sub(lost).saturating_sub(cancel);
                    // Stop only once the tail is down to the arithmetic budget,
                    // so Ziv sees more than a word of good bits.
                    if good >= budget.saturating_sub(4) && good > 64 {
                        return Some((result, good));
                    }
                    // Cancellation already ate this working precision; a later
                    // term will not restore those bits. Return the partial sum
                    // with the bits we actually have so the Ziv loop can raise
                    // wp (a local ziv_round treats None as a hard NaN).
                    if budget <= 64 {
                        return Some((result, good.max(1)));
                    }
                }
            }
        }
        k = k.checked_add(1)?;
        if formed == 0 && k > 64 {
            // All early terms vanished; keep scanning a little, then give up.
            if k > 1024 {
                return None;
            }
        }
    }
}

/// Whether 11.6.1 may be used: z>0 and m_min = ceil(ν-1/2)_+ is defined, and
/// |ν| is not comparable to or larger than z.
fn asym_gate(z: &ExactNum, nu: &ExactNum, wp: usize) -> bool {
    if !z.is_positive() {
        return false;
    }
    // |ν| ≥ z/2 → series / other methods only.
    let bound = ExactNum::div_ref(z, &two(wp), wp, DN);
    if nu.abs_cmp(&bound) != Some(-1) {
        return false;
    }
    true
}

fn ceil_m_min(nu: &ExactNum, wp: usize) -> Option<usize> {
    // m ≥ ν - 1/2, m ≥ 0. Compute ceil(ν-1/2) if positive.
    // If that integer does not fit in usize, this path is unavailable
    // (the caller falls back); that is an index-type limit, not a cap.
    let t = ExactNum::sub_ref(nu, &half(wp), wp, UP);
    if t.is_negative() || t.is_zero() {
        return Some(0);
    }
    let u = t.to_u64(RoundingMode::Up).ok()?;
    usize::try_from(u).ok()
}

/// DLMF 11.6.1: \(\mathbf{H}_\nu=Y_\nu+(1/\pi)\sum_{k=0}^{m-1}
/// \Gamma(k+1/2)(z/2)^{\nu-2k-1}/\Gamma(\nu+1/2-k)+R_m\), with
/// \(|R_m|\) less than the first omitted term when \(m\ge\nu-1/2\).
fn struve_asym(z: &ExactNum, nu: &ExactNum, wp: usize, cc: &mut Consts) -> Option<(ExactNum, i64)> {
    if !asym_gate(z, nu, wp) {
        return None;
    }
    let m_min = ceil_m_min(nu, wp)?;
    let halfz = ExactNum::div_ref(z, &two(wp), wp, N);
    let pi = cc.pi(wp, N);
    let halfv = half(wp);
    // t_0 = Γ(1/2) (z/2)^{ν-1} / (π Γ(ν+1/2))
    let g_half = halfv.gamma(wp, N, cc);
    let nu_ph = ExactNum::add_ref(nu, &halfv, wp, N);
    let g_nup = gamma_or_pole(&nu_ph, wp, cc)?;
    let pow0 = halfz.pow(&ExactNum::sub_ref(nu, &one(wp), wp, N), wp, N, cc);
    if !finite(&g_half) || !finite(&pow0) {
        return None;
    }
    let mut term = g_half.mul(&pow0, wp, N).div(&pi.mul(&g_nup, wp, N), wp, N);
    let mut sum = ExactNum::new(wp);
    let mut k = 0usize;
    let mut first_omitted;
    let mut terminated = false;
    loop {
        if k >= m_min {
            // Candidate truncation: first omitted is the current `term` (not yet added).
            first_omitted = term.clone();
            // Accept if |term| is small enough vs the running sum, or we terminated.
            break;
        }
        if !finite(&term) {
            return None;
        }
        // Pole in the next factor means the series of 11.6.1 terminates exactly
        // after the terms already added (1/Γ = 0 from here on).
        let den_shift = ExactNum::sub_ref(
            nu,
            &ExactNum::from_u64(k as u64, wp).add(&halfv, wp, N),
            wp,
            N,
        );
        // den_shift = ν - k - 1/2, used for the recurrence to t_{k+1}.
        if nonpos_int(&den_shift) && den_shift.is_zero() {
            // t_{k+1} has Γ(ν+1/2-(k+1)) at a non-positive integer starting at 0.
            terminated = true;
            first_omitted = ExactNum::new(wp);
            break;
        }
        if !term.is_zero() {
            sum = sum.add(&term, wp, N);
        }
        // t_{k+1} = t_k * (k+1/2) * 4 / z^2 * (ν-k-1/2)
        let num = ExactNum::from_u64(k as u64, wp).add(&halfv, wp, N).mul(
            &ExactNum::from_u8(4, wp),
            wp,
            N,
        );
        let z2 = z.sqr(wp, N);
        if den_shift.is_zero() || !finite(&den_shift) || !finite(&z2) || z2.is_zero() {
            terminated = true;
            first_omitted = ExactNum::new(wp);
            break;
        }
        let next = term
            .clone()
            .mul(&num, wp, N)
            .mul(&den_shift, wp, N)
            .div(&z2, wp, N);
        // If |next| > |term| we have passed the useful region before m_min.
        if k + 1 < m_min && finite(&next) && !term.is_zero() && next.abs_cmp(&term) == Some(1) {
            return None;
        }
        term = next;
        k = k.checked_add(1)?;
    }
    // Add more terms past m_min while they shrink and improve the bound, not
    // past the upturn and never without m ≥ m_min.
    if !terminated {
        let mut prev = term.clone();
        // `term` is the first omitted at m = k = m_min (nothing added at this k yet).
        // Add terms k = m_min, m_min+1, ... as long as they decrease, and stop
        // when the omitted term is small or starts growing.
        loop {
            if !finite(&term) {
                return None;
            }
            let den_shift = ExactNum::sub_ref(
                nu,
                &ExactNum::from_u64(k as u64, wp).add(&halfv, wp, N),
                wp,
                N,
            );
            if nonpos_int(&den_shift) {
                // Next term vanishes identically; current omitted can be added and
                // the remainder is then zero.
                if !term.is_zero() {
                    sum = sum.add(&term, wp, N);
                }
                first_omitted = ExactNum::new(wp);
                terminated = true;
                break;
            }
            // Stop if this omitted term already underflows the working grid
            // relative to the sum, or if the sequence has started to grow.
            let grow = k > m_min && !prev.is_zero() && term.abs_cmp(&prev) == Some(1);
            if grow {
                first_omitted = term.clone();
                break;
            }
            // If |term| is already far below the running (Y-less) sum, stop and
            // use it as the remainder bound. We still need Y to judge the scale;
            // keep going a few terms and leave the decision until after Y is
            // added, but cap the number of extra terms.
            if !term.is_zero() {
                sum = sum.add(&term, wp, N);
            }
            let num = ExactNum::from_u64(k as u64, wp).add(&halfv, wp, N).mul(
                &ExactNum::from_u8(4, wp),
                wp,
                N,
            );
            let z2 = z.sqr(wp, N);
            if den_shift.is_zero() || z2.is_zero() {
                first_omitted = ExactNum::new(wp);
                terminated = true;
                break;
            }
            prev = term;
            term = prev
                .clone()
                .mul(&num, wp, N)
                .mul(&den_shift, wp, N)
                .div(&z2, wp, N);
            first_omitted = term.clone();
            k = k.checked_add(1)?;
            if k.saturating_sub(m_min) > 65_536 {
                break;
            }
        }
    }
    let y = z.bessel_y(nu, wp, N, cc);
    if !finite(&y) {
        return None;
    }
    let h = ExactNum::add_ref(&y, &sum, wp, N);
    if !finite(&h) {
        return None;
    }
    let rem = if terminated { ExactNum::new(wp) } else { first_omitted };
    let rem = rem.abs();
    let peak = exp_of(&y).max(exp_of(&sum)).max(exp_of(&h));
    let lost = 2 * bitlen(k.saturating_add(2)) + 8;
    // A zero remainder (the 11.6.1 series terminated) has no tail; a zero's
    // exponent is not a bound and must not shrink bits_good.
    let tail = if rem.is_zero() { i64::MIN / 4 } else { exp_of(&rem) };
    let good = bits_good(wp, peak, exp_of(&h), lost, tail);
    // The remainder bound is only valid if we actually truncated at m ≥ m_min
    // (guaranteed) and |R| ≤ |first omitted|. If that is not small enough
    // relative to H after cancellation, refuse rather than guess.
    if good <= 0 {
        // Cancellation ate the working precision; Ziv will retry at a wider wp
        // only if we return a value with a (tiny) bit count. Zero or negative
        // would look like a certified empty result, so refuse this wp only by
        // handing a 1-bit claim: ziv_round will increase wp. A structural
        // failure (no bound) already returned None above.
        return Some((h, 1));
    }
    Some((h, good))
}

/// Tail of DLMF 11.5.2: \(\int_T^\infty e^{-z t}(1+t^2)^\mu\,dt\), \(z>0\), \(T\ge0\).
/// \(f'/f = -z + 2\mu t/(1+t^2)\). On \([T,\infty)\) the second term is at most
/// \(2|\mu| T/(1+T^2)\) in magnitude when \(T\ge1\), or \(2|\mu|T\) when \(T<1\)
/// (since \(1+t^2\ge1\)). If the resulting lower bound on \(-f'/f\) is positive,
/// the integral is at most \(f(T)\) over that bound.
fn k_tail(
    t: &ExactNum,
    z: &ExactNum,
    mu: &ExactNum,
    wp: usize,
    cc: &mut Consts,
) -> Option<ExactNum> {
    let tt = t.sqr(wp, UP);
    let one_tt_dn = one(wp).add(&t.sqr(wp, DN), wp, DN);
    if !one_tt_dn.is_positive() {
        return None;
    }
    // Upper bound of |2 μ t / (1+t^2)| on [T, ∞) for the decaying estimate.
    let two_mu_abs = two(wp).mul(&mu.clone().abs(), wp, UP);
    let frac_up = if ExactNum::from_u8(1, wp).abs_cmp(t) != Some(1) {
        // T >= 1: t/(1+t^2) <= T/(1+T^2)
        two_mu_abs
            .mul(t, wp, UP)
            .div(&one(wp).add(&tt, wp, DN), wp, UP)
    } else {
        // T < 1: t/(1+t^2) <= T
        two_mu_abs.mul(t, wp, UP)
    };
    let den_dn = ExactNum::sub_ref(z, &frac_up, wp, DN);
    if !den_dn.is_positive() {
        return None;
    }
    let one_tt = one(wp).add(&tt, wp, UP);
    let ln_u = one_tt.ln(wp, UP, cc);
    if !finite(&ln_u) {
        return None;
    }
    // f(T) = exp(-z T + μ ln(1+T^2)), upper bound.
    let expo_up = mu
        .clone()
        .mul(&ln_u, wp, UP)
        .sub(&ExactNum::mul_ref(z, t, wp, DN), wp, UP);
    let f_up = expo_up.exp(wp, UP, cc);
    if !finite(&f_up) {
        return None;
    }
    let tb = ExactNum::div_ref(&f_up, &den_dn, wp, UP);
    if finite(&tb) {
        Some(tb)
    } else {
        None
    }
}

fn pick_k_t(z: &ExactNum, mu: &ExactNum, wp: usize, cc: &mut Consts) -> Option<ExactNum> {
    let target = one(wp).ldexp(-((wp as i64) + 8), wp, N);
    // T ≈ (wp+16) ln 2 / z, then bump until the tail is small.
    let ln2 = two(wp).ln(wp, N, cc);
    let num = ExactNum::from_u64((wp as u64).saturating_add(16), wp).mul(&ln2, wp, N);
    let mut t = num.div(z, wp, UP);
    if !finite(&t) || t.is_negative() || t.is_zero() {
        t = one(wp);
    }
    for _ in 0..64 {
        if let Some(tb) = k_tail(&t, z, mu, wp, cc) {
            if tb.abs_cmp(&target) != Some(1) {
                return Some(t);
            }
        }
        t = t.mul(&ExactNum::from_u8(3, wp).div(&two(wp), wp, N), wp, UP);
        if t.exponent()? > 16 {
            return None;
        }
    }
    None
}

/// Derivatives of \(g=(1+t^2)^\mu\) at `t`, via
/// \((1+t^2)g^{(n)}=2t(\mu-(n-1))g^{(n-1)}+(n-1)(2\mu-(n-2))g^{(n-2)}\).
fn g_jet(
    t: &ExactNum,
    mu: &ExactNum,
    n_jet: usize,
    wp: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<Vec<ExactNum>> {
    let one_tt = one(wp).add(&t.sqr(wp, rm), wp, rm);
    if !one_tt.is_positive() {
        return None;
    }
    let ln = one_tt.ln(wp, rm, cc);
    if !finite(&ln) {
        return None;
    }
    let g0 = mu.clone().mul(&ln, wp, rm).exp(wp, rm, cc);
    if !finite(&g0) {
        return None;
    }
    let mut g: Vec<ExactNum> = vec![g0];
    for n in 1..=n_jet {
        // (1+t^2) g^{(n)} = 2t (μ-(n-1)) g^{(n-1)} + (n-1)(2μ-(n-2)) g^{(n-2)}
        let nm1 = ExactNum::from_u64(n as u64 - 1, wp);
        let coef1 = two(wp)
            .mul(t, wp, rm)
            .mul(&ExactNum::sub_ref(mu, &nm1, wp, rm), wp, rm);
        let mut rhs = ExactNum::mul_ref(&coef1, &g[n - 1], wp, rm);
        if n >= 2 {
            let nm2 = ExactNum::from_u64(n as u64 - 2, wp);
            let coef2 = nm1.mul(&two(wp).mul(mu, wp, rm).sub(&nm2, wp, rm), wp, rm);
            rhs = rhs.add(&coef2.mul(&g[n - 2], wp, rm), wp, rm);
        }
        let gn = rhs.div(&one_tt, wp, rm);
        if !finite(&gn) {
            return None;
        }
        g.push(gn);
    }
    Some(g)
}

/// Jet length for an 11.5.2 amplitude expansion at working precision `wp`.
fn k_n_jet(wp: usize) -> usize {
    let need = wp.saturating_add(16);
    // n ≈ (wp+16)/log2(e) at the near-optimal width h=ρ/e.
    ((need as u64 * 3) / 4).max(16) as usize
}

/// Panel width so a `k_n_jet`-term Cauchy remainder on the disk of
/// radius \(ρ=9/10\) (distance 1 from the real axis to \(\pm i\)) meets
/// \(2^{-(wp+16)}\). Also capped at \(ρ/2\).
fn k_panel_width(wp: usize) -> ExactNum {
    let n = k_n_jet(wp);
    let need = (wp as i64).saturating_add(16);
    let bits = (need + n as i64) / (n as i64 + 1);
    let rho = ExactNum::from_u8(9, wp).div(&ExactNum::from_u8(10, wp), wp, DN);
    let h = rho.clone().ldexp(-bits, wp, DN);
    let half_rho = rho.ldexp(-1, wp, DN);
    if finite(&h) && h.is_positive() && h.abs_cmp(&half_rho) == Some(-1) {
        h
    } else {
        half_rho
    }
}

/// \(\int_0^h e^{-z s}s^k\,ds\) for \(k=0,\ldots,n\).
///
/// For \(k\) at which the upper tail is at most half of \(k!/z^{k+1}\),
/// this is the complete moment minus the forward recurrence for
/// \(\int_h^\infty e^{-zs}s^k\,ds\). Otherwise the lower-end series
/// \(e^{-zh}h^{k+1}\sum (zh)^j/[(k+1)\cdots(k+1+j)]\) (all positive).
fn inc_exp_moments(
    z: &ExactNum,
    h: &ExactNum,
    n: usize,
    wp: usize,
    cc: &mut Consts,
) -> Option<Vec<ExactNum>> {
    if !z.is_positive() || !h.is_positive() {
        return None;
    }
    let zh = ExactNum::mul_ref(z, h, wp, N);
    if !finite(&zh) {
        return None;
    }
    let ln_h = h.ln(wp, N, cc);
    if !finite(&ln_h) {
        return None;
    }
    // I_n from the lower-end series (all positive, first term dominates
    // when n > zh), then the stable backward recurrence
    // I_{k-1} = (z I_k + e^{-zh} h^k)/k down to I_0.
    let lower = |k: usize, zh: &ExactNum, ln_h: &ExactNum, cc: &mut Consts| -> Option<ExactNum> {
        if k as u64 + 1 > u32::MAX as u64 {
            return None;
        }
        let kp1 = ExactNum::from_u64((k as u64).saturating_add(1), wp);
        let expo = kp1.clone().mul(ln_h, wp, N).sub(zh, wp, N);
        let mut term = expo.exp(wp, N, cc).div(&kp1, wp, N);
        if !finite(&term) {
            return None;
        }
        let mut sum = term.clone();
        let mut j = 0usize;
        loop {
            let den = ExactNum::from_u64((k as u64).saturating_add(2).saturating_add(j as u64), wp);
            term = term.mul(zh, wp, N).div(&den, wp, N);
            if !finite(&term) {
                return None;
            }
            sum = sum.add(&term, wp, N);
            let small = term.is_zero()
                || term
                    .clone()
                    .abs()
                    .abs_cmp(&one(wp).ldexp(-((wp as i64) + 16), wp, N).mul(
                        &sum.clone().abs(),
                        wp,
                        UP,
                    ))
                    != Some(1);
            if small && j >= 2 {
                break;
            }
            j = j.checked_add(1)?;
        }
        if finite(&sum) {
            Some(sum)
        } else {
            None
        }
    };
    let mut out: Vec<ExactNum> = (0..=n).map(|_| ExactNum::new(wp)).collect();
    out[n] = lower(n, &zh, &ln_h, cc)?;
    for k in (1..=n).rev() {
        let expo = ExactNum::from_u64(k as u64, wp)
            .mul(&ln_h, wp, N)
            .sub(&zh, wp, N);
        let pow = expo.exp(wp, N, cc);
        if !finite(&pow) {
            return None;
        }
        let num = ExactNum::mul_ref(z, &out[k], wp, N).add(&pow, wp, N);
        let ik = num.div(&ExactNum::from_u64(k as u64, wp), wp, N);
        if !finite(&ik) {
            return None;
        }
        out[k - 1] = ik;
    }
    Some(out)
}

/// Integral of \(f(t)=e^{-z t}(1+t^2)^\mu\) on `[a,b]`, by Taylor-expanding
/// only the amplitude \(g=(1+t^2)^\mu\) about \(a\) and integrating
/// \(e^{-z(a+s)}s^k\) exactly. The remainder uses a Cauchy bound on a disk
/// of radius \(9/10\), which stays a positive distance from the singularities
/// of \(g\) at \(\pm i\).
fn k_panel(
    a: &ExactNum,
    b: &ExactNum,
    z: &ExactNum,
    mu: &ExactNum,
    wp: usize,
    cc: &mut Consts,
) -> Option<(ExactNum, ExactNum)> {
    if a.is_negative() || ExactNum::sub_ref(b, a, wp, DN).is_negative() {
        return None;
    }
    let h = ExactNum::sub_ref(b, a, wp, N);
    if !h.is_positive() || !finite(&h) {
        return None;
    }
    let n = k_n_jet(wp);
    let g = g_jet(a, mu, n, wp, N, cc)?;
    if g.len() < n + 1 {
        return None;
    }
    let moments = inc_exp_moments(z, &h, n, wp, cc)?;
    let emza = ExactNum::mul_ref(z, a, wp, N).neg().exp(wp, N, cc);
    if !finite(&emza) {
        return None;
    }
    let mut integ = ExactNum::new(wp);
    let mut fac = one(wp);
    for k in 0..n {
        if k > 0 {
            fac = fac.mul(&ExactNum::from_u64(k as u64, wp), wp, N);
        }
        if fac.is_zero() || !finite(&fac) || !finite(&moments[k]) || !finite(&g[k]) {
            return None;
        }
        let term = g[k].clone().div(&fac, wp, N).mul(&moments[k], wp, N);
        integ = integ.add(&term, wp, N);
    }
    integ = integ.mul(&emza, wp, N);
    if !finite(&integ) {
        return None;
    }
    // |g^{(n)}(ξ)| ≤ n! / ρ^n * M,  |R| ≤ e^{-za} M I_n / ρ^n
    // ≤ e^{-za} M h^{n+1} / ((n+1) ρ^n).
    let rho = ExactNum::from_u8(9, wp).div(&ExactNum::from_u8(10, wp), wp, DN);
    if !rho.is_positive() {
        return None;
    }
    // |1+ζ^2| ≥ 1-ρ^2 on the ρ-neighborhood of the real axis;
    // |1+ζ^2| ≤ 1+(|x|+ρ)^2 with |x| ≤ max(a,b).
    let one_minus_rho2 = one(wp).sub(&rho.sqr(wp, UP), wp, DN);
    if !one_minus_rho2.is_positive() {
        return None;
    }
    let xmax = if b.abs_cmp(a) == Some(-1) { a.clone().abs() } else { b.clone().abs() };
    let max_abs = one(wp).add(&xmax.add(&rho, wp, UP).sqr(wp, UP), wp, UP);
    let m = if mu.is_negative() {
        let ln = one_minus_rho2.ln(wp, DN, cc);
        if !finite(&ln) {
            return None;
        }
        // (min)^μ with μ<0: smaller min gives a larger value.
        mu.clone().mul(&ln, wp, UP).exp(wp, UP, cc)
    } else {
        let ln = max_abs.ln(wp, UP, cc);
        if !finite(&ln) {
            return None;
        }
        mu.clone().mul(&ln, wp, UP).exp(wp, UP, cc)
    };
    if !finite(&m) {
        return None;
    }
    let e_up = ExactNum::mul_ref(z, a, wp, DN).neg().exp(wp, UP, cc);
    if !finite(&e_up) {
        return None;
    }
    let h_up = ExactNum::sub_ref(b, a, wp, UP);
    let h_over_rho = h_up.clone().div(&rho, wp, UP);
    let pow = h_over_rho.pow(&ExactNum::from_u64(n as u64, wp), wp, UP, cc);
    let np1 = ExactNum::from_u64((n as u64).saturating_add(1), wp);
    let rem = e_up
        .mul(&m, wp, UP)
        .mul(&h_up, wp, UP)
        .mul(&pow, wp, UP)
        .div(&np1, wp, UP);
    if !finite(&rem) {
        return None;
    }
    Some((integ, rem))
}

/// DLMF 11.5.2: \(K_\nu(z)=\dfrac{2(z/2)^\nu}{\sqrt{\pi}\,\Gamma(\nu+1/2)}
/// \int_0^\infty e^{-z t}(1+t^2)^{\nu-1/2}\,dt\), \(\Re z>0\).
/// Returns \(K_\nu\) and a good-bit count, or `None`.
fn struve_k_integral(
    z: &ExactNum,
    nu: &ExactNum,
    wp: usize,
    cc: &mut Consts,
) -> Option<(ExactNum, i64)> {
    if !z.is_positive() {
        return None;
    }
    let mu = ExactNum::sub_ref(nu, &half(wp), wp, N);
    let t_end = pick_k_t(z, &mu, wp, cc)?;
    let tb = k_tail(&t_end, z, &mu, wp, cc)?;
    let target = one(wp).ldexp(-((wp as i64) + 8), wp, N);
    // Panels sized so the Cauchy remainder on g, not a Taylor jet of the
    // exponential, meets the working precision (width ~ ρ/e, independent
    // of 1/z except for a last short panel).
    let mut h = k_panel_width(wp);
    if !finite(&h) || !h.is_positive() {
        h = t_end.clone();
    }
    if let Some(e) = ExactNum::div_ref(&t_end, &h, wp, UP).exponent() {
        if e > 16 {
            h = t_end.clone().div(&ExactNum::from_u32(65_536, wp), wp, UP);
        }
    }
    let mut stack: Vec<(ExactNum, ExactNum)> = Vec::new();
    let mut left = ExactNum::new(wp);
    loop {
        if ExactNum::sub_ref(&t_end, &left, wp, DN).is_negative()
            || ExactNum::sub_ref(&t_end, &left, wp, DN).is_zero()
        {
            break;
        }
        let mut right = left.clone().add(&h, wp, UP);
        if right.abs_cmp(&t_end) != Some(-1) {
            right = t_end.clone();
        }
        stack.push((left.clone(), right.clone()));
        if right.abs_cmp(&t_end) != Some(-1) {
            break;
        }
        left = right;
    }
    let mut acc = ExactNum::new(wp);
    let mut rem_acc = tb.clone();
    let mut panels = 0usize;
    while let Some((a, b)) = stack.pop() {
        panels = panels.saturating_add(1);
        match k_panel(&a, &b, z, &mu, wp, cc) {
            Some((val, rem)) if rem.abs_cmp(&target) != Some(1) => {
                acc = acc.add(&val, wp, N);
                rem_acc = rem_acc.add(&rem, wp, UP);
            }
            Some(_) | None => {
                let mid = ExactNum::add_ref(&a, &b, wp, N).ldexp(-1, wp, N);
                if ExactNum::sub_ref(&b, &a, wp, N).exponent()? < -((wp as i64) + 16) {
                    return None;
                }
                stack.push((mid.clone(), b));
                stack.push((a, mid));
            }
        }
    }
    if !finite(&acc) {
        return None;
    }
    // C = 2 (z/2)^ν / (√π Γ(ν+1/2))
    let nu_ph = ExactNum::add_ref(nu, &half(wp), wp, N);
    if nonpos_int(&nu_ph) {
        return None;
    }
    let g = gamma_or_pole(&nu_ph, wp, cc)?;
    let pi = cc.pi(wp, N);
    let sqrt_pi = pi.sqrt(wp, N);
    let halfz = ExactNum::div_ref(z, &two(wp), wp, N);
    let c = if halfz.is_positive() {
        let ln_h = halfz.ln(wp, N, cc);
        let ln_g = nu_ph.ln_gamma(wp, N, cc);
        let ln_spi = sqrt_pi.ln(wp, N, cc);
        if finite(&ln_h) && finite(&ln_g) && finite(&ln_spi) && !g.is_negative() {
            let ln = ExactNum::mul_ref(nu, &ln_h, wp, N)
                .add(&two(wp).ln(wp, N, cc), wp, N)
                .sub(&ln_spi, wp, N)
                .sub(&ln_g, wp, N);
            ln.exp(wp, N, cc)
        } else {
            let pow = halfz.pow(nu, wp, N, cc);
            two(wp).mul(&pow, wp, N).div(&sqrt_pi.mul(&g, wp, N), wp, N)
        }
    } else {
        return None;
    };
    if !finite(&c) {
        return None;
    }
    let k = c.clone().mul(&acc, wp, N);
    if !finite(&k) {
        return None;
    }
    let rem = c.abs().mul(&rem_acc, wp, UP);
    let lost = 2 * bitlen(panels) + 8;
    let tail = if rem.is_zero() || !finite(&rem) { i64::MIN / 4 } else { exp_of(&rem) };
    let good = bits_good(wp, exp_of(&k), exp_of(&k), lost, tail);
    if good <= 8 {
        return None;
    }
    Some((k, good))
}

/// \(\mathbf{H}_\nu = Y_\nu + K_\nu\) from 11.5.2 and \(Y_\nu\).
fn struve_via_k(
    z: &ExactNum,
    nu: &ExactNum,
    wp: usize,
    cc: &mut Consts,
) -> Option<(ExactNum, i64)> {
    let (k, gk) = struve_k_integral(z, nu, wp, cc)?;
    let y = z.bessel_y(nu, wp, N, cc);
    if !finite(&y) {
        return None;
    }
    let h = ExactNum::add_ref(&y, &k, wp, N);
    if !finite(&h) {
        return None;
    }
    let peak = exp_of(&y).max(exp_of(&k)).max(exp_of(&h));
    let good = bits_good(wp, peak, exp_of(&h), 8, i64::MIN / 4).min(gk);
    if good <= 0 {
        return Some((h, 1));
    }
    Some((h, good))
}

fn struve_h_positive(
    z: &ExactNum,
    nu: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> ExactNum {
    let zc = z.clone();
    let nuc = nu.clone();
    // Pick the method from (ν, z) before starting any kernel. 11.6.1 is
    // used only when its e^{-z} remainder floor is enough for p; Maclaurin
    // only when its term count is cheap, or when |ν| ≳ |z| (entire, little
    // cancellation). Otherwise 11.5.2 + Y_ν.
    let use_asym = zc.is_positive() && exp_neg_z_enough(&zc, p);
    // 11.2.1 when it is cheap, or when |ν| is comparable to |z| (11.6.1
    // is gated at |ν| >= z/2) and the term count is not a hang.
    let nu_ge_half_z = {
        let halfz = zc.clone().ldexp(-1, 64, N);
        nuc.clone().abs().abs_cmp(&halfz) != Some(-1)
    };
    let k_safe = maclaurin_k_est(&zc).is_some_and(|k| k <= 16_384);
    let use_series = series_cheap(&zc, &nuc, p) || (nu_ge_half_z && k_safe);
    ziv_round(p, rm, |wp| {
        if use_asym {
            if let Some(r) = struve_asym(&zc, &nuc, wp, cc) {
                return Some(r);
            }
        }
        if use_series {
            if let Some(r) = struve_series(&zc, &nuc, wp, cc) {
                return Some(r);
            }
        }
        // 11.5.2 + Y is the bounded large-z path when 11.6.1 is not tight
        // enough. Never fall through into Maclaurin at large |z|.
        struve_via_k(&zc, &nuc, wp, cc)
    })
}

/// Struve \(\mathbf{H}_\nu(z)\).
pub(super) fn struve_h_word(
    z: &ExactNum,
    nu: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> ExactNum {
    if !finite(z) || !finite(nu) {
        return cat_nan();
    }
    if z.is_zero() {
        return struve_at_zero(z, nu, p, rm, cc);
    }
    if z.is_negative() {
        if !nu.is_int() {
            return cat_nan();
        }
        let x = z.clone().abs();
        // H_n(-x) = (-1)^{n+1} H_n(x): flip the sign if n is even. The inner
        // call must use the negated mode so the final rounding direction is
        // the one the caller asked for (Up/Down swap under a sign flip).
        let odd = is_odd_int(nu, work_p(p));
        let inner_rm = if odd { rm } else { mode_for_neg(rm) };
        let h = struve_h_positive(&x, nu, p, inner_rm, cc);
        if !finite(&h) && !h.is_inf() {
            return h;
        }
        if odd {
            h
        } else {
            h.neg()
        }
    } else {
        struve_h_positive(z, nu, p, rm, cc)
    }
}

// ---------------- Anger–Weber ----------------

/// \(\mathbf{J}_\nu(0)=\sin(\nu\pi)/(\nu\pi)\) (1 at ν=0),
/// \(\mathbf{E}_\nu(0)=(1-\cos(\nu\pi))/(\nu\pi)\) (0 at ν=0).
fn anger_weber_at_zero(
    weber: bool,
    nu: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> ExactNum {
    let wrk = work_p(p).saturating_add(WORD_BIT_SIZE);
    if nu.is_zero() {
        return finish(if weber { signed_zero(false, p) } else { one(p) }, p, rm);
    }
    if nu.is_int() {
        if !weber {
            return finish(signed_zero(false, p), p, rm);
        }
        // E_n(0) = 0 for even n, 2/(nπ) for odd n (sign of n).
        if !is_odd_int(nu, wrk) {
            return finish(signed_zero(false, p), p, rm);
        }
        let pi = cc.pi(wrk, N);
        let v = two(wrk).div(&ExactNum::mul_ref(nu, &pi, wrk, N), wrk, N);
        return finish(v, p, rm);
    }
    let pi = cc.pi(wrk, N);
    let a = ExactNum::mul_ref(nu, &pi, wrk, N);
    if weber {
        let c = a.cos(wrk, N, cc);
        let num = ExactNum::sub_ref(&one(wrk), &c, wrk, N);
        finish(num.div(&a, wrk, N), p, rm)
    } else {
        let s = a.sin(wrk, N, cc);
        finish(s.div(&a, wrk, N), p, rm)
    }
}

/// DLMF 11.10.10–11 Maclaurin sums \(S_1,S_2\).
fn maclaurin_s(
    z: &ExactNum,
    nu: &ExactNum,
    wp: usize,
    cc: &mut Consts,
) -> Option<((ExactNum, ExactNum), i64)> {
    let halfz = ExactNum::div_ref(z, &two(wp), wp, N);
    let z2 = halfz.sqr(wp, N);
    let half_nu = ExactNum::div_ref(nu, &two(wp), wp, N);
    let mut s1 = ExactNum::new(wp);
    let mut s2 = ExactNum::new(wp);
    let mut peak = i64::MIN / 4;
    // S1 term k=0: 1 / (Γ(1+ν/2) Γ(1-ν/2))
    let form_s1 = |k: usize, cc: &mut Consts| -> Option<ExactNum> {
        let kk = ExactNum::from_u64(k as u64, wp);
        let a = ExactNum::add_ref(&kk, &half_nu, wp, N).add(&one(wp), wp, N);
        let b = ExactNum::sub_ref(&kk, &half_nu, wp, N).add(&one(wp), wp, N);
        let ga = gamma_or_pole(&a, wp, cc)?;
        let gb = gamma_or_pole(&b, wp, cc)?;
        let mut zk = one(wp);
        for _ in 0..k {
            zk = zk.mul(&z2, wp, N);
        }
        if !k.is_multiple_of(2) {
            zk = zk.neg();
        }
        Some(zk.div(&ga.mul(&gb, wp, N), wp, N))
    };
    let form_s2 = |k: usize, cc: &mut Consts| -> Option<ExactNum> {
        let kk = ExactNum::from_u64(k as u64, wp);
        let three_half = ExactNum::from_u8(3, wp).div(&two(wp), wp, N);
        let a = ExactNum::add_ref(&kk, &half_nu, wp, N).add(&three_half, wp, N);
        let b = ExactNum::sub_ref(&kk, &half_nu, wp, N).add(&three_half, wp, N);
        let ga = gamma_or_pole(&a, wp, cc)?;
        let gb = gamma_or_pole(&b, wp, cc)?;
        let mut zk = halfz.clone();
        for _ in 0..k {
            zk = zk.mul(&z2, wp, N);
        }
        if !k.is_multiple_of(2) {
            zk = zk.neg();
        }
        Some(zk.div(&ga.mul(&gb, wp, N), wp, N))
    };
    let mut t1 = form_s1(0, cc).unwrap_or_else(|| ExactNum::new(wp));
    let mut t2 = form_s2(0, cc).unwrap_or_else(|| ExactNum::new(wp));
    let mut k = 0usize;
    loop {
        if finite(&t1) && !t1.is_zero() {
            peak = peak.max(exp_of(&t1));
            s1 = s1.add(&t1, wp, N);
        }
        if finite(&t2) && !t2.is_zero() {
            peak = peak.max(exp_of(&t2));
            s2 = s2.add(&t2, wp, N);
        }
        // Ratio for S1: |t_{k+1}/t_k| = (z/2)^2 / |((k+1+ν/2)(k+1-ν/2))|
        let kp1 = ExactNum::from_u64((k as u64).saturating_add(1), wp);
        let d1 = ExactNum::add_ref(&kp1, &half_nu, wp, N).mul(
            &ExactNum::sub_ref(&kp1, &half_nu, wp, N),
            wp,
            N,
        );
        let ratio_ok = if finite(&d1) && !d1.is_zero() {
            let r = ExactNum::div_ref(&z2, &d1, wp, UP).abs();
            r.compare(&half(wp)) != Some(1)
        } else {
            false
        };
        if ratio_ok && k >= 1 {
            let tail = exp_of(&t1).max(exp_of(&t2));
            let res_e = exp_of(&s1).max(exp_of(&s2));
            let lost = 2 * bitlen(k + 2) + 6;
            let good = bits_good(wp, peak, res_e, lost, tail);
            let cancel = (peak - res_e).max(0);
            let budget = (wp as i64).saturating_sub(lost).saturating_sub(cancel);
            if good >= budget.saturating_sub(4) && good > 64 {
                return Some(((s1, s2), good));
            }
            // Tight first wp: return a low good-bit count so Ziv can raise
            // it. None is a hard NaN (E_100(100), E_{500.5}(100)).
            if budget <= 64 {
                return Some(((s1, s2), good.max(1)));
            }
        }
        k = k.checked_add(1)?;
        t1 = form_s1(k, cc).unwrap_or_else(|| ExactNum::new(wp));
        t2 = form_s2(k, cc).unwrap_or_else(|| ExactNum::new(wp));
    }
}

fn anger_weber_maclaurin(
    weber: bool,
    z: &ExactNum,
    nu: &ExactNum,
    wp: usize,
    cc: &mut Consts,
) -> Option<(ExactNum, i64)> {
    let ((s1, s2), good) = maclaurin_s(z, nu, wp, cc)?;
    let half_pi_nu = ExactNum::mul_ref(nu, &cc.pi(wp, N), wp, N).ldexp(-1, wp, N);
    let (c, s) = (half_pi_nu.cos(wp, N, cc), half_pi_nu.sin(wp, N, cc));
    // J = cos(πν/2) S1 + sin(πν/2) S2
    // E = sin(πν/2) S1 - cos(πν/2) S2
    let v = if weber {
        ExactNum::mul_ref(&s, &s1, wp, N).sub(&ExactNum::mul_ref(&c, &s2, wp, N), wp, N)
    } else {
        ExactNum::mul_ref(&c, &s1, wp, N).add(&ExactNum::mul_ref(&s, &s2, wp, N), wp, N)
    };
    if !finite(&v) {
        return None;
    }
    Some((v, good))
}

/// Factorial \(n!\) as an ExactNum, rounded up when `up` is set (for majorants).
fn factorial_rm(n: usize, p: usize, rm: RoundingMode) -> ExactNum {
    let mut acc = one(p);
    for k in 2..=n {
        acc = acc.mul(&ExactNum::from_u64(k as u64, p), p, rm);
    }
    acc
}

fn binom(n: usize, j: usize, p: usize, rm: RoundingMode) -> ExactNum {
    if j > n {
        return ExactNum::new(p);
    }
    let j = j.min(n - j);
    let mut acc = one(p);
    for i in 0..j {
        acc = acc.mul(&ExactNum::from_u64((n - i) as u64, p), p, rm).div(
            &ExactNum::from_u64((i as u64).saturating_add(1), p),
            p,
            rm,
        );
    }
    acc
}

/// Integrand \(f(t)=\exp(-\nu t-z\sinh t)\) and a majorant of its jet on `[a,b]`.
fn f_and_gprime(
    t: &ExactNum,
    nu: &ExactNum,
    z: &ExactNum,
    wp: usize,
    cc: &mut Consts,
) -> Option<(ExactNum, ExactNum, ExactNum, ExactNum)> {
    // returns (f, sinh t, cosh t, g' = -ν - z cosh t)
    let (sh, ch) = t.sinh_cosh(wp, N, cc);
    if !finite(&sh) || !finite(&ch) {
        return None;
    }
    let g = ExactNum::mul_ref(nu, t, wp, N)
        .add(&ExactNum::mul_ref(z, &sh, wp, N), wp, N)
        .neg();
    let f = g.exp(wp, N, cc);
    if !finite(&f) {
        return None;
    }
    let gp = ExactNum::mul_ref(z, &ch, wp, N).add(nu, wp, N).neg();
    Some((f, sh, ch, gp))
}

/// Upper bound of \(\int_T^\infty f\) for \(z>0\), \(\nu>-z\):
/// \(f(T)/(\nu+z\cosh T)\), computed with directed rounding so it is an
/// upper bound (sinh is convex, \(f\) is decreasing).
fn tail_bound(
    t: &ExactNum,
    nu: &ExactNum,
    z: &ExactNum,
    wp: usize,
    cc: &mut Consts,
) -> Option<ExactNum> {
    let (sh, ch) = t.sinh_cosh(wp, DN, cc); // DN on t>0: sinh, cosh are increasing, DN is a lower bound
                                            // Wait: we need an UPPER bound of f(T)/(ν+z cosh T).
                                            // f=exp(-νT-z sinh T). Larger f means smaller (νT+z sinh T).
                                            // So we want sinh T and T rounded down, νT+z sinh rounded down if all positive.
                                            // Denominator ν+z cosh T rounded down (positive). Then quotient rounded up.
    let (sh_dn, ch_dn) = t.sinh_cosh(wp, DN, cc);
    if !finite(&sh_dn) || !finite(&ch_dn) {
        return None;
    }
    let _ = sh;
    let _ = ch;
    let s_dn = ExactNum::mul_ref(nu, t, wp, DN).add(&ExactNum::mul_ref(z, &sh_dn, wp, DN), wp, DN);
    if !s_dn.is_positive() {
        return None;
    }
    let f_up = s_dn.neg().exp(wp, UP, cc);
    let den_dn = ExactNum::mul_ref(z, &ch_dn, wp, DN).add(nu, wp, DN);
    if !den_dn.is_positive() {
        return None;
    }
    let tb = ExactNum::div_ref(&f_up, &den_dn, wp, UP);
    if finite(&tb) {
        Some(tb)
    } else {
        None
    }
}

/// Pick \(T>0\) so the tail bound is below \(2^{-(wp+8)}\).
fn pick_t(nu: &ExactNum, z: &ExactNum, wp: usize, cc: &mut Consts) -> Option<ExactNum> {
    let target = one(wp).ldexp(-((wp as i64) + 8), wp, N);
    // T ≈ asinh((wp+16) ln2 / z) is a good start when z is not tiny.
    let ln2 = two(wp).ln(wp, N, cc);
    let num = ExactNum::from_u64((wp as u64).saturating_add(16), wp).mul(&ln2, wp, N);
    let mut t = num.div(z, wp, UP).asinh(wp, UP, cc);
    if !finite(&t) || t.is_negative() || t.is_zero() {
        t = one(wp);
    }
    for _ in 0..64 {
        let tb = tail_bound(&t, nu, z, wp, cc)?;
        if tb.abs_cmp(&target) != Some(1) {
            return Some(t);
        }
        t = t.mul(&ExactNum::from_u8(3, wp).div(&two(wp), wp, N), wp, UP);
        if t.exponent()? > 16 {
            return None;
        }
    }
    None
}

/// Taylor-panel integral of \(f\) on `[a,b]` plus a majorant remainder.
#[allow(clippy::needless_range_loop)]
fn panel_integral(
    a: &ExactNum,
    b: &ExactNum,
    nu: &ExactNum,
    z: &ExactNum,
    wp: usize,
    cc: &mut Consts,
) -> Option<(ExactNum, ExactNum)> {
    let mid = ExactNum::add_ref(a, b, wp, N).ldexp(-1, wp, N);
    let h = ExactNum::sub_ref(b, a, wp, N);
    let half_h = h.clone().ldexp(-1, wp, N);
    let (f0, sh0, ch0, _) = f_and_gprime(&mid, nu, z, wp, cc)?;
    // Jet u_n = f^{(n)}(mid) via Leibniz: u_n = Σ_{j=0}^{n-1} C(n-1,j) a_{n-j} u_j
    // a_1 = -ν - z cosh(mid); a_k = -z sinh(mid) (k even ≥2); -z cosh(mid) (k odd ≥3).
    let n_jet = (16 + wp / 16).max(12);
    let mut jet: Vec<ExactNum> = vec![f0.clone()];
    let a1 = ExactNum::mul_ref(z, &ch0, wp, N).add(nu, wp, N).neg();
    for n in 1..=n_jet {
        let ak = |k: usize| -> ExactNum {
            if k == 1 {
                a1.clone()
            } else if k.is_multiple_of(2) {
                ExactNum::mul_ref(z, &sh0, wp, N).neg()
            } else {
                ExactNum::mul_ref(z, &ch0, wp, N).neg()
            }
        };
        let mut un = ExactNum::new(wp);
        for j in 0..n {
            let c = binom(n - 1, j, wp, N);
            let term = c.mul(&ak(n - j), wp, N).mul(&jet[j], wp, N);
            un = un.add(&term, wp, N);
        }
        jet.push(un);
    }
    // Integrate Σ u_k (t-mid)^k / k! over [a,b] = Σ u_k / k! * (h/2)^{k+1} / (k+1) * (1-(-1)^{k+1})
    // i.e. only even k survive: 2 (h/2)^{k+1} /(k+1)! * u_k
    let mut integ = ExactNum::new(wp);
    for k in 0..=n_jet {
        if !k.is_multiple_of(2) {
            continue;
        }
        let pow = half_h.pow(
            &ExactNum::from_u64((k as u64).saturating_add(1), wp),
            wp,
            N,
            cc,
        );
        let fac = factorial_rm(k + 1, wp, N);
        let term = jet[k]
            .clone()
            .mul(&pow, wp, N)
            .div(&fac, wp, N)
            .mul(&two(wp), wp, N);
        integ = integ.add(&term, wp, N);
    }
    // Majorant of |f^{(n_jet+1)}| on [a,b]:
    // A_1 = |ν| + z cosh(t_max), A_k = z cosh(t_max) (k≥2), t_max = b (all positive here).
    let t_max = if b.abs_cmp(a) == Some(-1) { a.clone() } else { b.clone() };
    let ch_max = t_max.cosh(wp, UP, cc);
    if !finite(&ch_max) {
        return None;
    }
    let a1u = ExactNum::mul_ref(z, &ch_max, wp, UP).add(&nu.clone().abs(), wp, UP);
    let aku = ExactNum::mul_ref(z, &ch_max, wp, UP);
    // U_0 = max f on the panel. f is decreasing, so max is f(a), upper rounded.
    let (fa, _, _, _) = f_and_gprime(a, nu, z, wp, cc)?;
    // Slightly inflate fa to cover rounding in the sample.
    let mut u0 = fa.clone().mul(
        &one(wp).add(&one(wp).ldexp(-((wp as i64) / 2), wp, N), wp, UP),
        wp,
        UP,
    );
    if u0.abs_cmp(&fa) != Some(1) {
        u0 = fa.clone();
        // one ulp up
        let _ = u0.set_precision(wp, UP);
    }
    let mut major: Vec<ExactNum> = vec![u0];
    for n in 1..=n_jet + 1 {
        let mut un = ExactNum::new(wp);
        for j in 0..n {
            let c = binom(n - 1, j, wp, UP);
            let ak = if n - j == 1 { a1u.clone() } else { aku.clone() };
            let term = c.mul(&ak, wp, UP).mul(&major[j], wp, UP);
            un = un.add(&term, wp, UP);
        }
        major.push(un);
    }
    // |R| ≤ (h/2)^{n_jet+1} / (n_jet+1)! * U_{n_jet+1} * h  (bound the integral of the remainder)
    let pow = half_h.pow(
        &ExactNum::from_u64((n_jet as u64).saturating_add(1), wp),
        wp,
        UP,
        cc,
    );
    let fac = factorial_rm(n_jet + 1, wp, DN);
    if fac.is_zero() || !finite(&fac) || !finite(&pow) {
        return None;
    }
    let rem = pow
        .div(&fac, wp, UP)
        .mul(&major[n_jet + 1], wp, UP)
        .mul(&h, wp, UP);
    if !finite(&rem) || !finite(&integ) {
        return None;
    }
    Some((integ, rem))
}

/// \(A_\nu(z)=(1/\pi)\int_0^\infty\exp(-\nu t-z\sinh t)\,dt\), \(z>0\), \(\nu>-z\).
/// Returns the integral (not yet divided by π) and good bits, or `None`.
fn anger_a_integral(
    nu: &ExactNum,
    z: &ExactNum,
    wp: usize,
    cc: &mut Consts,
) -> Option<(ExactNum, i64)> {
    if !z.is_positive() {
        return None;
    }
    // Monotone decreasing iff ν + z > 0.
    if ExactNum::add_ref(nu, z, wp, DN).is_negative() || ExactNum::add_ref(nu, z, wp, DN).is_zero()
    {
        return None;
    }
    let t_end = pick_t(nu, z, wp, cc)?;
    let tb = tail_bound(&t_end, nu, z, wp, cc)?;
    // Panels sized to the decay scale ~ 1/(ν+z), then bisect any whose
    // remainder is still above the target. A single [0, T] panel at large
    // z produces many failed jets, so the first grid follows the decay.
    let target = one(wp).ldexp(-((wp as i64) + 8), wp, N);
    let scale = ExactNum::add_ref(nu, z, wp, UP);
    let mut h = if scale.is_positive() && finite(&scale) {
        one(wp).div(&scale, wp, UP)
    } else {
        t_end.clone()
    };
    if !finite(&h) || !h.is_positive() {
        h = t_end.clone();
    }
    if let Some(e) = ExactNum::div_ref(&t_end, &h, wp, UP).exponent() {
        if e > 16 {
            h = t_end.clone().div(&ExactNum::from_u32(65_536, wp), wp, UP);
        }
    }
    let mut stack: Vec<(ExactNum, ExactNum)> = Vec::new();
    let mut left = ExactNum::new(wp);
    loop {
        if ExactNum::sub_ref(&t_end, &left, wp, DN).is_negative()
            || ExactNum::sub_ref(&t_end, &left, wp, DN).is_zero()
        {
            break;
        }
        let mut right = left.clone().add(&h, wp, UP);
        if right.abs_cmp(&t_end) != Some(-1) {
            right = t_end.clone();
        }
        stack.push((left.clone(), right.clone()));
        if right.abs_cmp(&t_end) != Some(-1) {
            break;
        }
        left = right;
    }
    let mut acc = ExactNum::new(wp);
    let mut rem_acc = tb.clone();
    let mut panels = 0usize;
    while let Some((a, b)) = stack.pop() {
        panels = panels.saturating_add(1);
        match panel_integral(&a, &b, nu, z, wp, cc) {
            Some((val, rem)) if rem.abs_cmp(&target) != Some(1) => {
                acc = acc.add(&val, wp, N);
                rem_acc = rem_acc.add(&rem, wp, UP);
            }
            Some(_) | None => {
                let mid = ExactNum::add_ref(&a, &b, wp, N).ldexp(-1, wp, N);
                if ExactNum::sub_ref(&b, &a, wp, N).exponent()? < -((wp as i64) + 16) {
                    return None;
                }
                stack.push((mid.clone(), b));
                stack.push((a, mid));
            }
        }
    }
    if !finite(&acc) {
        return None;
    }
    let lost = 2 * bitlen(panels) + 8;
    let good = bits_good(wp, exp_of(&acc), exp_of(&acc), lost, exp_of(&rem_acc));
    if good <= 8 {
        return None;
    }
    Some((acc, good))
}

/// Majorant of \(|h(u)|\) on the circle \(|u|=1/2\), where
/// \(h(u)=(\sqrt{1+u^2}-u)^\nu/\sqrt{1+u^2}\) is the Watson integrand of
/// \(A_\nu\) after \(u=\sinh t\). The nearest singularities are at \(u=\pm i\).
fn h_circle_majorant(nu: &ExactNum, wp: usize, cc: &mut Consts) -> Option<ExactNum> {
    let _ = cc;
    // |1+u^2| lies in [3/4, 5/4], so |sqrt(1+u^2)| is in [sqrt(3/4), sqrt(5/4)]
    // and |sqrt(1+u^2)-u| is in [sqrt(3/4)-1/2, sqrt(5/4)+1/2].
    let three_q = ExactNum::from_u8(3, wp).ldexp(-2, wp, DN);
    let five_q = ExactNum::from_u8(5, wp).ldexp(-2, wp, UP);
    let lo = three_q.sqrt(wp, DN);
    if !finite(&lo) {
        return None;
    }
    let hi = five_q.sqrt(wp, UP);
    if !finite(&hi) {
        return None;
    }
    let halfv = half(wp);
    let diff_lo = ExactNum::sub_ref(&lo, &halfv, wp, DN);
    if !diff_lo.is_positive() {
        return None;
    }
    let sum_hi = ExactNum::add_ref(&hi, &halfv, wp, UP);
    // |(sqrt-u)^ν| = |sqrt-u|^ν for real ν: use the larger endpoint.
    let base = if nu.is_negative() { diff_lo } else { sum_hi };
    let pow = if nu.is_zero() { one(wp) } else { base.pow(nu, wp, UP, cc) };
    if !finite(&pow) {
        return None;
    }
    // |1/sqrt(1+u^2)| ≤ 1/sqrt(3/4)
    let inv = one(wp).div(&lo, wp, UP);
    let m = pow.mul(&inv, wp, UP);
    if finite(&m) {
        Some(m)
    } else {
        None
    }
}

/// Upper bound of \(\int_\tau^\infty e^{-z u} u^k\,du\) for \(k=0,\ldots,k_{\max}\).
fn incgamma_tail_uk(
    z: &ExactNum,
    tau: &ExactNum,
    kmax: usize,
    wp: usize,
    cc: &mut Consts,
) -> Option<Vec<ExactNum>> {
    let mut z_dn = z.clone();
    let _ = z_dn.set_precision(wp, DN);
    if !z_dn.is_positive() {
        return None;
    }
    let ztau_dn = ExactNum::mul_ref(&z_dn, tau, wp, DN);
    if !ztau_dn.is_positive() {
        return None;
    }
    let exp_up = ztau_dn.neg().exp(wp, UP, cc);
    if !finite(&exp_up) {
        return None;
    }
    let e_over_z = ExactNum::div_ref(&exp_up, &z_dn, wp, UP);
    let mut out = Vec::with_capacity(kmax + 1);
    out.push(e_over_z.clone());
    let mut tau_pow = one(wp);
    for k in 1..=kmax {
        tau_pow = tau_pow.mul(tau, wp, UP);
        let term1 = tau_pow.clone().mul(&e_over_z, wp, UP);
        let k_over_z = ExactNum::from_u64(k as u64, wp).div(&z_dn, wp, UP);
        let nxt = term1.add(&k_over_z.mul(&out[k - 1], wp, UP), wp, UP);
        if !finite(&nxt) {
            return None;
        }
        out.push(nxt);
    }
    Some(out)
}

/// Upper bound of \(\int_\tau^\infty e^{-z u} |h(u)|\,du\).
fn h_real_tail(
    nu: &ExactNum,
    z: &ExactNum,
    tau: &ExactNum,
    wp: usize,
    cc: &mut Consts,
) -> Option<ExactNum> {
    let mut z_dn = z.clone();
    let _ = z_dn.set_precision(wp, DN);
    if !z_dn.is_positive() {
        return None;
    }
    if !nu.is_negative() {
        let ztau_dn = ExactNum::mul_ref(&z_dn, tau, wp, DN);
        if !ztau_dn.is_positive() {
            return None;
        }
        let exp_up = ztau_dn.neg().exp(wp, UP, cc);
        if !finite(&exp_up) {
            return None;
        }
        let tb = ExactNum::div_ref(&exp_up, &z_dn, wp, UP);
        return if finite(&tb) { Some(tb) } else { None };
    }
    // ν<0: |h(u)| ≤ (1+2u)^{-ν} ≤ exp(-2ν u), so the tail is at most
    // e^{-(z+2ν)τ}/(z+2ν) provided z+2ν>0 (true when |ν|<z/2).
    let two_nu = nu.clone().ldexp(1, wp, DN); // more negative if ν<0
    let den = ExactNum::add_ref(z, &two_nu, wp, DN);
    if !den.is_positive() {
        return None;
    }
    let ztau_dn = ExactNum::mul_ref(&den, tau, wp, DN);
    if !ztau_dn.is_positive() {
        return None;
    }
    let exp_up = ztau_dn.neg().exp(wp, UP, cc);
    if !finite(&exp_up) {
        return None;
    }
    let tb = ExactNum::div_ref(&exp_up, &den, wp, UP);
    if finite(&tb) {
        Some(tb)
    } else {
        None
    }
}

/// DLMF 11.11.4 via Watson's lemma: \(A_\nu(z)\sim(1/(\pi z))\sum F_k(\nu)/z^{2k}
/// -(\nu/(\pi z^2))\sum G_k(\nu)/z^{2k}\), truncated at the smallest term, with
/// a Cauchy remainder on \(|u|<1/2\) in the \(u=\sinh t\) coordinate.
///
/// \(F_0=G_0=1\), \(F_k=(\nu^2-1^2)\cdots(\nu^2-(2k-1)^2)\),
/// \(G_k=(\nu^2-2^2)\cdots(\nu^2-(2k)^2)\).
fn anger_a_laplace(
    nu: &ExactNum,
    z: &ExactNum,
    wp: usize,
    cc: &mut Consts,
) -> Option<(ExactNum, i64)> {
    if !z.is_positive() {
        return None;
    }
    if ExactNum::add_ref(nu, z, wp, DN).is_negative() || ExactNum::add_ref(nu, z, wp, DN).is_zero()
    {
        return None;
    }
    // The Cauchy split is at τ=1/4, radius ρ=1/2. If z is too small the bound
    // cannot beat a few bits, so refuse and let the panel integrator run.
    if z.exponent()? < 3 {
        return None;
    }
    let m_maj = h_circle_majorant(nu, wp, cc)?;
    let rho = half(wp);
    let tau = one(wp).ldexp(-2, wp, N);
    let pi = cc.pi(wp, N);
    let mut z_dn = z.clone();
    let _ = z_dn.set_precision(wp, DN);
    if !z_dn.is_positive() {
        return None;
    }
    let h_tail = h_real_tail(nu, z, &tau, wp, cc)?;

    // Build terms t_n = a_n n! / (π z^{n+1}) and a running Cauchy bound.
    // a_{2p} = F_p / (2p)!, a_{2p+1} = -ν G_p / (2p+1)!.
    let nu2 = nu.sqr(wp, N);
    let mut f = one(wp);
    let mut g = one(wp);
    let mut sum = ExactNum::new(wp);
    let mut prev_abs = ExactNum::new(wp);
    let mut best: Option<(ExactNum, ExactNum, usize)> = None;
    let mut used = 0usize;

    // Precompute I_k up to the cap we might use, lazily extended.
    let mut ik: Vec<ExactNum> = Vec::new();
    let ensure_ik = |ik: &mut Vec<ExactNum>, n: usize, cc: &mut Consts| -> Option<()> {
        if ik.len() > n {
            return Some(());
        }
        let got = incgamma_tail_uk(z, &tau, n.max(ik.len().saturating_add(8)), wp, cc)?;
        *ik = got;
        Some(())
    };

    // Running factors for the Cauchy rem, updated in O(1) per term.
    // After step n (n_terms = n+1): fac = (n+1)!, rho_pow = ρ^{n+1},
    // zpow_rem = z^{n+2}, r_poly = (M/π) Σ_{k=0}^{n} ρ^{-k} I_k.
    let mut fac = one(wp);
    let mut rho_pow = one(wp);
    let mut zpow_rem = z_dn.clone();
    let mut zpow_term = one(wp); // becomes z^{n+1} at the start of step n
    let mut r_poly = ExactNum::new(wp);
    let mut rho_k = one(wp); // ρ^k for the newest poly term

    let mut n = 0usize;
    loop {
        let p = n / 2;
        zpow_term = zpow_term.mul(z, wp, N);
        if !finite(&zpow_term) || zpow_term.is_zero() {
            return None;
        }
        let term = if n % 2 == 0 {
            // F_p / (π z^{2p+1})
            f.clone()
                .div(&ExactNum::mul_ref(&pi, &zpow_term, wp, N), wp, N)
        } else {
            // -ν G_p / (π z^{2p+2})
            nu.clone()
                .neg()
                .mul(&g, wp, N)
                .div(&ExactNum::mul_ref(&pi, &zpow_term, wp, N), wp, N)
        };
        if !finite(&term) {
            return None;
        }
        let grow = n > 0 && !prev_abs.is_zero() && term.clone().abs().abs_cmp(&prev_abs) == Some(1);
        if grow {
            break;
        }
        sum = sum.add(&term, wp, N);
        used = n + 1;
        prev_abs = term.abs();

        // Update F or G for the next pair.
        if n % 2 == 0 {
            // After even term 2p, prepare G already present; F_{p+1} is needed
            // after the following odd term. Nothing to do to F/G here except
            // after the odd term we step both.
        } else {
            // Finished the pair (F_p, G_p). Advance to p+1.
            let odd = ExactNum::from_u64((p as u64).saturating_mul(2).saturating_add(1), wp);
            let even = ExactNum::from_u64((p as u64).saturating_mul(2).saturating_add(2), wp);
            f = f.mul(&ExactNum::sub_ref(&nu2, &odd.sqr(wp, N), wp, N), wp, N);
            g = g.mul(&ExactNum::sub_ref(&nu2, &even.sqr(wp, N), wp, N), wp, N);
            if !finite(&f) || !finite(&g) {
                break;
            }
        }

        // Cauchy remainder after N = n+1 terms (indices 0..n inclusive).
        // |R_series| ≤ 2 M ρ^{-N} N! / (π z^{N+1})
        // |R_poly tail| ≤ (M/π) Σ_{k=0}^{N-1} ρ^{-k} I_k
        // |R_h tail| ≤ h_tail / π
        let n_terms = n + 1;
        fac = fac.mul(&ExactNum::from_u64(n_terms as u64, wp), wp, UP);
        rho_pow = rho_pow.mul(&rho, wp, DN);
        zpow_rem = zpow_rem.mul(&z_dn, wp, DN);
        if !finite(&fac)
            || !finite(&rho_pow)
            || !finite(&zpow_rem)
            || rho_pow.is_zero()
            || zpow_rem.is_zero()
        {
            break;
        }
        let r_series = two(wp)
            .mul(&m_maj, wp, UP)
            .mul(&fac, wp, UP)
            .div(&rho_pow, wp, UP)
            .div(&zpow_rem, wp, UP)
            .div(&pi, wp, UP);
        ensure_ik(&mut ik, n, cc)?;
        // Add the new k=n term of the poly tail (k=0 uses rho_k=1).
        if n > 0 {
            rho_k = rho_k.mul(&rho, wp, DN);
        }
        if rho_k.is_zero() || !finite(&rho_k) {
            return None;
        }
        if n < ik.len() {
            let piece = m_maj
                .clone()
                .mul(&ik[n], wp, UP)
                .div(&rho_k, wp, UP)
                .div(&pi, wp, UP);
            r_poly = r_poly.add(&piece, wp, UP);
        }
        let r_h = ExactNum::div_ref(&h_tail, &pi, wp, UP);
        let rem = r_series.add(&r_poly, wp, UP).add(&r_h, wp, UP);
        if !finite(&rem) {
            break;
        }
        let improved = match &best {
            None => {
                best = Some((sum.clone(), rem.clone(), n_terms));
                true
            }
            Some((_, brem, _)) if rem.abs_cmp(brem) == Some(-1) => {
                best = Some((sum.clone(), rem.clone(), n_terms));
                true
            }
            Some(_) => false,
        };
        // Stop if the remainder is far below the working precision relative
        // to the sum, or if the bound has stalled (asymptotic terms will
        // start growing, or a tighter rem is not coming).
        if let Some((bs, brem, _)) = &best {
            if !bs.is_zero()
                && brem.abs_cmp(&one(wp).ldexp(-((wp as i64) + 8), wp, N).mul(bs, wp, UP))
                    != Some(1)
            {
                break;
            }
        }
        if !improved {
            break;
        }
        match n.checked_add(1) {
            Some(n1) => n = n1,
            None => break,
        }
    }

    let (val, rem, n_terms) = best?;
    if !finite(&val) {
        return None;
    }
    // Rounding/term-count budget. +32 was so conservative that a rem well
    // below 2^{-wp} still missed Ziv's p+64 bar at 8192 bits (g=8268 at
    // wp=8320, need 8256) and the large_z fallback fed the short count
    // to Ziv, which raised wp forever.
    let lost = 2 * bitlen(n_terms) + 16;
    let good = bits_good(wp, exp_of(&val), exp_of(&val), lost, exp_of(&rem));
    if good <= 8 {
        return None;
    }
    let _ = used;
    Some((val, good))
}

fn anger_a(nu: &ExactNum, z: &ExactNum, wp: usize, cc: &mut Consts) -> Option<(ExactNum, i64)> {
    if let Some((v, g)) = anger_a_laplace(nu, z, wp, cc) {
        // A loose remainder must not be fed to the Ziv loop (it would raise
        // wp forever and never switch to panels). Use Laplace only when the
        // rem meets this working precision.
        let tight = g >= (wp as i64).saturating_sub(64);
        if tight {
            return Some((v, g));
        }
    }
    let (i, good) = anger_a_integral(nu, z, wp, cc)?;
    let pi = cc.pi(wp, N);
    let a = i.div(&pi, wp, N);
    if !finite(&a) {
        return None;
    }
    Some((a, good))
}

/// DLMF 11.10.15–16, only for \(z>0\).
fn anger_weber_connection(
    weber: bool,
    z: &ExactNum,
    nu: &ExactNum,
    wp: usize,
    cc: &mut Consts,
) -> Option<(ExactNum, i64)> {
    if !z.is_positive() {
        return None;
    }
    // 11.10.4 needs ν > -z for A_ν. Weber also needs A_{-ν}, so ν < z.
    // The old |ν| < z/2 cut refused J_{500.5}(1000) and E_{500.5}(1000),
    // which are finite.
    if ExactNum::add_ref(nu, z, wp, DN).is_negative() || ExactNum::add_ref(nu, z, wp, DN).is_zero()
    {
        return None;
    }
    if weber
        && (ExactNum::sub_ref(z, nu, wp, DN).is_negative()
            || ExactNum::sub_ref(z, nu, wp, DN).is_zero())
    {
        return None;
    }
    let pi_nu = ExactNum::mul_ref(nu, &cc.pi(wp, N), wp, N);
    let (s_pi, c_pi) = (pi_nu.sin(wp, N, cc), pi_nu.cos(wp, N, cc));
    if weber {
        let y = z.bessel_y(nu, wp, N, cc);
        if !finite(&y) {
            return None;
        }
        let (ap, gp) = anger_a(nu, z, wp, cc)?;
        let nneg = nu.clone().neg();
        let (am, gm) = anger_a(&nneg, z, wp, cc)?;
        // E = -Y - cos(πν) A_ν - A_{-ν}
        let v = y
            .clone()
            .neg()
            .sub(&ExactNum::mul_ref(&c_pi, &ap, wp, N), wp, N)
            .sub(&am, wp, N);
        if !finite(&v) {
            return None;
        }
        let peak = exp_of(&y).max(exp_of(&ap)).max(exp_of(&am));
        let good = bits_good(wp, peak, exp_of(&v), 8, i64::MIN / 4)
            .min(gp)
            .min(gm);
        if good <= 0 {
            return Some((v, 1));
        }
        Some((v, good))
    } else {
        let j = z.bessel_j_nu(nu, wp, N, cc);
        if !finite(&j) {
            return None;
        }
        // For integer ν, sin(πν)=0 and J_Anger = J_Bessel; A_ν is not needed.
        if nu.is_int() {
            return Some((j, wp as i64 - 8));
        }
        let (ap, gp) = anger_a(nu, z, wp, cc)?;
        let v = ExactNum::add_ref(&j, &ExactNum::mul_ref(&s_pi, &ap, wp, N), wp, N);
        if !finite(&v) {
            return None;
        }
        let peak = exp_of(&j).max(exp_of(&ap));
        let good = bits_good(wp, peak, exp_of(&v), 8, i64::MIN / 4).min(gp);
        if good <= 0 {
            return Some((v, 1));
        }
        Some((v, good))
    }
}

fn anger_weber_positive(
    weber: bool,
    z: &ExactNum,
    nu: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> ExactNum {
    let zc = z.clone();
    let nuc = nu.clone();
    let cheap = series_cheap(&zc, &nuc, p);
    ziv_round(p, rm, |wp| {
        if !weber && nuc.is_int() {
            let j = zc.bessel_j_nu(&nuc, wp, N, cc);
            if finite(&j) {
                return Some((j, wp as i64 - 8));
            }
            return None;
        }
        // Pick from (ν, z) before starting a kernel. Maclaurin only when
        // ~|z|/2 terms plus cancellation bits are cheap; otherwise the
        // 11.10.15–16 connection (Laplace if the bound is tight, else
        // A_ν panels sized to 1/(ν+z)).
        if cheap {
            if let Some(r) = anger_weber_maclaurin(weber, &zc, &nuc, wp, cc) {
                return Some(r);
            }
            return anger_weber_connection(weber, &zc, &nuc, wp, cc);
        }
        if let Some(r) = anger_weber_connection(weber, &zc, &nuc, wp, cc) {
            return Some(r);
        }
        // Connection refuses A_{-ν} when ν >= z (Weber). Maclaurin of
        // 11.10.8–11 is entire; use it when the term count is not a hang.
        let k_safe = maclaurin_k_est(&zc).is_some_and(|k| k <= 16_384);
        if k_safe {
            if let Some(r) = anger_weber_maclaurin(weber, &zc, &nuc, wp, cc) {
                return Some(r);
            }
        }
        None
    })
}

fn anger_weber_word(
    weber: bool,
    z: &ExactNum,
    nu: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> ExactNum {
    if !finite(z) || !finite(nu) {
        return cat_nan();
    }
    if z.is_zero() {
        return anger_weber_at_zero(weber, nu, p, rm, cc);
    }
    if z.is_negative() {
        // 11.10.12: J_ν(-z) = J_{-ν}(z),  E_ν(-z) = -E_{-ν}(z).
        let x = z.clone().abs();
        let nneg = nu.clone().neg();
        // 11.10.12: J_ν(-z) = J_{-ν}(z) (no sign flip), E_ν(-z) = -E_{-ν}(z).
        let inner_rm = if weber { mode_for_neg(rm) } else { rm };
        let v = anger_weber_positive(weber, &x, &nneg, p, inner_rm, cc);
        if weber {
            v.neg()
        } else {
            v
        }
    } else {
        anger_weber_positive(weber, z, nu, p, rm, cc)
    }
}

/// Anger \(\mathbf{J}_\nu(z)\).
pub(super) fn anger_j_word(
    z: &ExactNum,
    nu: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> ExactNum {
    anger_weber_word(false, z, nu, p, rm, cc)
}

/// Weber \(\mathbf{E}_\nu(z)\).
pub(super) fn weber_e_word(
    z: &ExactNum,
    nu: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> ExactNum {
    anger_weber_word(true, z, nu, p, rm, cc)
}

// Silence a helper kept for structure; bits_good is the real one.
#[allow(dead_code)]
#[cfg(test)]
mod aws_bounds {
    use super::*;
    use crate::Consts;
    use crate::Error;
    use crate::ExactNum;
    use crate::RoundingMode;

    fn rel_closer_than(a: &ExactNum, b: &ExactNum, bits: i64, wp: usize) -> bool {
        if !finite(a) || !finite(b) {
            return false;
        }
        if a.is_zero() && b.is_zero() {
            return true;
        }
        let diff = ExactNum::sub_ref(a, b, wp, N).abs();
        let scale = if b.clone().abs().abs_cmp(&a.clone().abs()) == Some(1) {
            b.clone().abs()
        } else {
            a.clone().abs()
        };
        if scale.is_zero() {
            return diff.is_zero();
        }
        let tol = scale.ldexp(-bits, wp, N);
        diff.abs_cmp(&tol) != Some(1)
    }

    /// The 11.6.1 recurrence must multiply by (ν−k−1/2). Dividing produced a
    /// relative error of 2.7e-5 for H_0(100), exactly the size of the t_1
    /// discrepancy (mpmath: K_0≈0.00636556, t_0=1/(50π), t_1=-t_0/z^2).
    #[test]
    fn zf56_struve_h0_100_recurrence_and_remainder() -> Result<(), Error> {
        let mut cc = Consts::new()?;
        let wp = 192usize;
        let z = ExactNum::from_u32(100, wp);
        let nu = ExactNum::new(wp);
        let (h, good) = struve_asym(&z, &nu, wp, &mut cc).expect("11.6.1 must apply at z=100, ν=0");
        assert!(finite(&h), "H_0(100) not finite");
        // Y_0(100) + 1/(50π) is the one-term approximation; the true K_0 is
        // smaller than t_0 by about t_0/z^2 ≈ 6.4e-7. A correctly summed
        // 11.6.1 series plus a correct Y_0 must beat the old 2.7e-5 error
        // by many orders of magnitude.
        let y = z.bessel_y(&nu, wp, N, &mut cc);
        let pi = cc.pi(wp, N);
        let t0 = one(wp).div(&pi.mul(&ExactNum::from_u32(50, wp), wp, N), wp, N);
        let one_term = ExactNum::add_ref(&y, &t0, wp, N);
        let diff = ExactNum::sub_ref(&h, &one_term, wp, N).abs();
        // |t_1| = t_0/z^2 ≈ 6.37e-7; |H - (Y+t_0)| should be about that, not
        // the 4×-too-large 2.55e-6 of the inverted recurrence.
        let t1_abs = t0.div(&z.sqr(wp, N), wp, N);
        // Allow a factor of two for t_2 and rounding, but not the factor of four.
        let cap = t1_abs.mul(&ExactNum::from_u8(2, wp), wp, N);
        assert!(
            diff.abs_cmp(&cap) != Some(1),
            "H_0(100) is farther from Y+t_0 than 2|t_1|; recurrence is wrong"
        );
        // The first omitted term after a valid m is far below 1e-20, so good
        // bits must be well above 40 (the old 2.7e-5 error is ~15 bits).
        assert!(
            good >= 40,
            "bits_good={good} is too small for a correct 11.6.1 remainder at z=100"
        );
        // Word-level value at 64 bits must match a 192-bit evaluation rounded
        // to 64 (Ziv consistency), and beat the old relative error.
        let word = struve_h_word(&z, &nu, 64, RoundingMode::ToEven, &mut cc);
        let mut rounded = h.clone();
        let _ = rounded.set_precision(64, RoundingMode::ToEven);
        assert_eq!(
            word.compare(&rounded),
            Some(0),
            "64-bit ToEven H_0(100) disagrees with the 192-bit kernel rounded down"
        );
        Ok(())
    }

    #[test]
    fn zf56_struve_series_and_asym_agree_in_overlap() -> Result<(), Error> {
        let mut cc = Consts::new()?;
        let wp = 256usize;
        let z = ExactNum::from_u8(20, wp);
        let nu = ExactNum::new(wp);
        let (s, sg) = struve_series(&z, &nu, wp, &mut cc).expect("series at z=20");
        let (a, ag) = struve_asym(&z, &nu, wp, &mut cc).expect("asym at z=20");
        let bits = sg.min(ag).saturating_sub(8).max(16);
        assert!(
            rel_closer_than(&s, &a, bits, wp),
            "series and 11.6.1 disagree at z=20, ν=0 (good {sg} vs {ag})"
        );
        Ok(())
    }

    #[test]
    fn zf56_struve_half_integer_terminates() -> Result<(), Error> {
        let mut cc = Consts::new()?;
        let wp = 192usize;
        let z = ExactNum::from_u32(100, wp);
        let nu = half(wp);
        let (h, good) =
            struve_asym(&z, &nu, wp, &mut cc).expect("11.6.1 terminates after t_0 for ν=1/2");
        assert!(good >= 40, "bits_good={good} for terminated H_{{1/2}}(100)");
        // K_{1/2} = t_0 exactly, so H = Y + t_0.
        let y = z.bessel_y(&nu, wp, N, &mut cc);
        let pi = cc.pi(wp, N);
        // t_0 = Γ(1/2)(z/2)^{ν-1} / (π Γ(ν+1/2)) = sqrt(2/(π z)) for ν=1/2.
        let t0 = two(wp).div(&pi.mul(&z, wp, N), wp, N).sqrt(wp, N);
        assert!(finite(&t0), "sqrt(2/(πz))");
        let expect = ExactNum::add_ref(&y, &t0, wp, N);
        assert!(
            rel_closer_than(&h, &expect, 80, wp),
            "H_{{1/2}}(100) must equal Y_{{1/2}}(100)+sqrt(2/(πz))"
        );
        Ok(())
    }

    #[test]
    fn zf56_anger_a_laplace_large_z() -> Result<(), Error> {
        let mut cc = Consts::new()?;
        let wp = 192usize;
        let z = ExactNum::from_u32(100, wp);
        let nu = ExactNum::new(wp);
        let (a, good) =
            anger_a_laplace(&nu, &z, wp, &mut cc).expect("Laplace/Watson must apply at z=100, ν=0");
        assert!(finite(&a) && good >= 16, "A_0(100) good={good}");
        // Leading term is 1/(πz); the next is O(1/z^3) since ν=0 kills the 1/z^2 term.
        let pi = cc.pi(wp, N);
        let lead = one(wp).div(&pi.mul(&z, wp, N), wp, N);
        let diff = ExactNum::sub_ref(&a, &lead, wp, N).abs();
        // F_1(0) = -1, so the next term is -1/(π z^3) ≈ -3.2e-8. Stay within 2× of that.
        let t2 = lead.div(&z.sqr(wp, N), wp, N);
        assert!(
            diff.abs_cmp(&t2.mul(&two(wp), wp, N)) != Some(1),
            "A_0(100) is farther from 1/(πz) than 2/ (π z^3)"
        );
        // Huge z must return quickly with a tight bound (this is the path that
        // hung in the panel integrator).
        let zbig = one(wp).ldexp(20, wp, N);
        let (ab, gb) =
            anger_a_laplace(&nu, &zbig, wp, &mut cc).expect("Laplace must apply at z=2^20");
        assert!(finite(&ab) && gb >= 40, "A_0(2^20) good={gb}");
        Ok(())
    }

    #[test]
    fn zf56_h_neg_half_integers_at_zero_are_zero() -> Result<(), Error> {
        let mut cc = Consts::new()?;
        let p = 64usize;
        // RC-4: H_{-3/2}(0)=H_{-5/2}(0)=0, contrary to mpmath's ±Inf.
        for times in [-3i32, -5] {
            let nu = ExactNum::from_i128(times as i128, 64).div(&two(64), 64, N);
            let z = ExactNum::new(64);
            let h = struve_h_word(&z, &nu, p, RoundingMode::ToEven, &mut cc);
            assert!(h.is_zero(), "H_{{{times}/2}}(0) must be 0, got a non-zero");
        }
        Ok(())
    }

    #[test]
    fn zf56_laplace_good_bits_at_high_wp() -> Result<(), Error> {
        let mut cc = Consts::new()?;
        for wp in [192usize, 2176, 8320] {
            let z = one(wp).ldexp(20, wp, N);
            let nu = half(wp);
            match anger_a_laplace(&nu, &z, wp, &mut cc) {
                Some((a, g)) => {
                    let tight = g >= (wp as i64).saturating_sub(64);
                    assert!(finite(&a), "A_{{1/2}}(2^20) not finite at wp={wp}");
                    assert!(
                        tight,
                        "A_{{1/2}}(2^20) Laplace rem not tight at wp={wp}, good={g}"
                    );
                }
                None => panic!("anger_a_laplace returned None at wp={wp}"),
            }
        }
        Ok(())
    }
}
