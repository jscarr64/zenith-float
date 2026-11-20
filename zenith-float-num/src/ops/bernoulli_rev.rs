//! Exact even Bernoulli numbers by Arb/FLINT's reverse ζ walk (`bernoulli_rev`).
//!
//! For even `n ≥ 32`, `|B_n| = 2 n! / (2π)^n · ζ(n)`, the denominator is the
//! von Staudt–Clausen product, and the numerator is the unique nearby integer.
//! The walk starts at `nmax` and steps down by 2, updating `2 n! / (2π)^n` and
//! the odd inverse powers `j^{-n}` rather than rebuilding ζ at each index
//! (Johansson, FLINT `bernoulli_rev_{init,next}`).

use super::factorial;
use super::from_usize_p;
use super::pow_ui;
use super::primes_upto;
use crate::common::util::round_p;
use crate::defs::Error;
use crate::defs::Exponent;
use crate::defs::RoundingMode;
use crate::defs::Sign;
use crate::defs::WORD_BIT_SIZE;
use crate::integer::ExactInt;
use crate::num::ExactNumNumber;
use crate::ops::consts::Consts;
use alloc::vec::Vec;

/// Avoid the reverse walk for very small `n` (FLINT `BERNOULLI_REV_MIN`).
const REV_MIN: usize = 32;

fn bitlen_usize(n: usize) -> usize {
    (usize::BITS - n.leading_zeros()) as usize
}

/// `log2(n)` in Q32.32 for `n >= 1`, via the artanh series for `ln(1+x)`.
/// Accurate to better than `2^{-30}`; used only to size the reverse walk
/// the way FLINT does (`arith_bernoulli_number_size`, `bernoulli_zeta_terms`).
fn log2_q32(n: u64) -> u64 {
    if n <= 1 {
        return 0;
    }
    let k = n.ilog2();
    let m_q32 = if k >= 32 { n >> (k - 32) } else { n << (32 - k) };
    let x_q32 = m_q32 - (1u64 << 32);
    // z = x / (2 + x) in Q48, so z ∈ [0, 1/3] and the series converges quickly.
    let z_q48 = (u128::from(x_q32) << 48) / ((2u128 << 32) + u128::from(x_q32));
    let z2 = (z_q48 * z_q48) >> 48;
    let mut zpow = z_q48;
    let mut sum = 0u128;
    for kterm in 0..12u128 {
        sum += zpow / (2 * kterm + 1);
        zpow = (zpow * z2) >> 48;
    }
    let ln_q48 = sum * 2;
    // ln 2 * 2^48, rounded nearest.
    const LN2_Q48: u128 = 195_103_586_505_167;
    let log2_frac_q32 = (ln_q48 << 32) / LN2_Q48;
    (u64::from(k) << 32).saturating_add(log2_frac_q32 as u64)
}

/// Upper bound on the bits of `|B_n|` (FLINT `arith_bernoulli_number_size`).
fn bernoulli_number_size_bits(n: usize) -> usize {
    if n <= 14 {
        return 8;
    }
    // 2 + (n+1) log2(n+1) - n log2(2πe) + 2, with log2(2πe) = 4.0941911703612822.
    const LOG2_2PIE_Q32: u128 = 17_584_417_180;
    let log2n1 = u128::from(log2_q32((n as u64).saturating_add(1)));
    let x_q32 = (4u128 << 32)
        .saturating_add((n as u128 + 1).saturating_mul(log2n1))
        .saturating_sub((n as u128).saturating_mul(LOG2_2PIE_Q32));
    ((x_q32 >> 32) as usize).max(16)
}

/// FLINT `bernoulli_denom_size`: `0.5 log2(n) · n^{1.066 / ln ln n}`.
/// Integer fixed-point, padded by 8 bits so the bound stays above FLINT's.
fn bernoulli_denom_size(n: usize) -> usize {
    if n < 32 {
        return 16;
    }
    const LN2_Q32: u128 = 2_977_044_472;
    const C1066_Q32: u128 = 4_578_038_080; // 1.066 * 2^32
    let log2n = u128::from(log2_q32(n as u64));
    let ln_n = (log2n * LN2_Q32) >> 32;
    if ln_n < 2 {
        return 16;
    }
    // ln(ln n) = ln(ln_n / 2^32) = log2(ln_n)*ln2 - 32*ln2.
    let log2_ln = u128::from(log2_q32(ln_n as u64));
    let log2_ln_n = log2_ln.saturating_sub(32u128 << 32);
    let ln_ln = (log2_ln_n * LN2_Q32) >> 32;
    if ln_ln == 0 {
        return 16;
    }
    let alpha = (C1066_Q32 << 32) / ln_ln;
    let exp_q32 = (alpha * log2n) >> 32;
    let npow = exp2_q32_ceil(exp_q32).max(1);
    let bits = ((log2n >> 1) * npow) >> 32;
    (bits as usize).saturating_add(8)
}

fn exp2_q32_ceil(x_q32: u128) -> u128 {
    let ip = (x_q32 >> 32) as u32;
    let frac = x_q32 & 0xffff_ffff;
    const LN2_Q32: u128 = 2_977_044_472;
    let u = (frac * LN2_Q32) >> 32;
    let u2 = (u * u) >> 32;
    let u3 = (u2 * u) >> 32;
    let u4 = (u2 * u2) >> 32;
    // 2^frac ≈ 1 + u + u^2/2 + u^3/6 + u^4/24, plus a few ulps so we ceil.
    let two_frac = (1u128 << 32) + u + (u2 / 2) + (u3 / 6) + (u4 / 24) + (1 << 16);
    if ip >= 96 {
        return u128::MAX;
    }
    (two_frac << ip) >> 32
}

/// `denom(B_n)` by von Staudt–Clausen: `∏_{p-1 | n} p` for even `n`.
/// `primes` must contain every prime up to `n+1`.
fn bernoulli_denom_int(n: usize, primes: &[u32]) -> Result<ExactInt, Error> {
    if n % 2 == 1 {
        return Ok(if n == 1 {
            ExactInt::from_u64(2)
        } else {
            ExactInt::one()
        });
    }
    let mut den = ExactInt::from_u64(6);
    for &p in primes.iter().skip(2) {
        let p = p as usize;
        if p.saturating_sub(1) > n {
            break;
        }
        if n % (p - 1) == 0 {
            den.mul_u64_assign(p as u64);
        }
    }
    Ok(den)
}

fn bernoulli_denom_ui(n: usize, primes: &[u32]) -> Result<ExactNumNumber, Error> {
    int_to_num(&bernoulli_denom_int(n, primes)?)
}

fn bernoulli_global_prec(nmax: usize) -> Result<usize, Error> {
    // FLINT `bernoulli_global_prec`: size + denom_size, no extra 128-bit pad.
    let wp = bernoulli_number_size_bits(nmax).saturating_add(bernoulli_denom_size(nmax));
    Ok(round_p(wp.max(16)))
}

fn bernoulli_zeta_terms(s: usize, prec: usize) -> usize {
    // FLINT: N = 2^{(prec+1)/(s-1)}, truncated to slong, then made odd.
    let den = s.saturating_sub(1).max(1);
    let target = prec.saturating_add(1);
    let q = (target / den).min((usize::BITS - 2) as usize);
    let r = target % den;
    let base = 1usize << q;
    if r == 0 {
        let mut n = base;
        if n.is_multiple_of(2) {
            n = n.saturating_add(1);
        }
        return n.max(3);
    }
    let target_q32 = u128::from(target as u64) << 32;
    let den_u = u128::from(den as u64);
    let mut lo = base;
    let mut hi = base.saturating_mul(2).saturating_sub(1).max(lo + 1);
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        let v = u128::from(log2_q32(mid as u64)).saturating_mul(den_u);
        if v >= target_q32 {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    let mut n = lo.max(3);
    if n.is_multiple_of(2) {
        n = n.saturating_add(1);
    }
    n
}

fn bernoulli_power_prec(j: usize, s: usize, wp: usize) -> usize {
    // FLINT: wp - s log2(j), at least 10; we keep a 2-bit guard and min 16.
    let log2j_q32 = u128::from(log2_q32(j.max(2) as u64));
    let sub = ((u128::from(s as u64).saturating_mul(log2j_q32)) >> 32) as usize;
    let p = wp.saturating_sub(sub).saturating_add(2);
    round_p(p.max(16))
}

fn nearest_unique_int(x: &ExactNumNumber) -> Result<ExactNumNumber, Error> {
    let p = round_p(x.mantissa_max_bit_len().max(64));
    let mut quarter = ExactNumNumber::from_word(1, p)?;
    quarter = quarter.ldexp(-2, p, RoundingMode::None)?;
    let mut half = ExactNumNumber::from_word(1, p)?;
    half = half.ldexp(-1, p, RoundingMode::None)?;
    let mut n = if x.is_zero() {
        x.clone()?
    } else if x.is_positive() {
        x.add(&half, p, RoundingMode::None)?.floor()?
    } else {
        x.sub(&half, p, RoundingMode::None)?.ceil()?
    };
    let err = x.sub(&n, p, RoundingMode::None)?.abs()?;
    if err.cmp(&quarter) >= 0 {
        return Err(Error::PrecisionRetryExhausted);
    }
    // Unique nearby integer: exact by construction (the source may carry an
    // inexact flag from the ζ approximation).
    n.set_inexact(false);
    Ok(n)
}

fn pair_to_ints(num: ExactNumNumber, den: ExactNumNumber) -> Result<(ExactInt, ExactInt), Error> {
    Ok((num_to_int(&num)?, num_to_int(&den)?))
}

fn exp_as_exponent(n: usize) -> Result<Exponent, Error> {
    Exponent::try_from(n).map_err(|_| Error::InvalidArgument)
}

/// Truncate a finite number toward zero as an integer (mantissa is little-endian).
pub(super) fn num_to_int(x: &ExactNumNumber) -> Result<ExactInt, Error> {
    if x.is_zero() {
        return Ok(ExactInt::zero());
    }
    let (words, _n, sign, e, _) = x.as_raw_parts();
    let p = words.len().saturating_mul(WORD_BIT_SIZE);
    if p == 0 {
        return Ok(ExactInt::zero());
    }
    let mut z = ExactInt::from_le_words(sign, words);
    let e = i128::from(e);
    let p = i128::from(p as u64);
    if e >= p {
        let sh = usize::try_from(e - p).map_err(|_| Error::InvalidArgument)?;
        z = z.shl(sh);
    } else {
        let sh = usize::try_from(p - e).map_err(|_| Error::InvalidArgument)?;
        z = z.shr_towards_zero(sh);
    }
    Ok(z)
}

fn int_to_num(x: &ExactInt) -> Result<ExactNumNumber, Error> {
    if x.is_zero() {
        return ExactNumNumber::from_word(0, WORD_BIT_SIZE);
    }
    let limbs = x.as_limbs();
    let e = Exponent::try_from(limbs.len().saturating_mul(WORD_BIT_SIZE))
        .map_err(|_| Error::InvalidArgument)?;
    let sign = if x.is_negative() { Sign::Neg } else { Sign::Pos };
    ExactNumNumber::from_words(limbs, sign, e)
}


/// FLINT `arb_div_2expm1_ui`: y = z / (2^n - 1).
fn div_2expm1(z: &ExactNumNumber, n: usize, p: usize) -> Result<ExactNumNumber, Error> {
    if n == 0 {
        return Err(Error::InvalidArgument);
    }
    // Small n: divide by the word 2^n - 1 (FLINT `n < FLINT_BITS`).
    if n < WORD_BIT_SIZE {
        let den = from_usize_p((1usize << n) - 1, WORD_BIT_SIZE)?;
        return z.div(&den, p, RoundingMode::None);
    }
    // Mid n: a general division is cheaper than a long geometric sum
    // (FLINT: n < 1024 + prec/32).
    if n < 1024 + p / 32 {
        let mut two_n = ExactNumNumber::from_word(1, p)?;
        two_n = two_n.ldexp(exp_as_exponent(n)?, p, RoundingMode::None)?;
        let one = ExactNumNumber::from_word(1, p)?;
        let den = two_n.sub(&one, p, RoundingMode::None)?;
        return z.div(&den, p, RoundingMode::None);
    }
    // Large n: z/(2^n-1) = sum_{k>=1} z 2^{-k n}.
    let nexp = exp_as_exponent(n)?;
    let mut term = z.clone()?.ldexp(-nexp, p, RoundingMode::ToEven)?;
    let mut h = term.clone()?;
    let last = p / n + 1;
    for _ in 2..=last {
        term = term.ldexp(-nexp, p, RoundingMode::ToEven)?;
        if term.is_zero() {
            break;
        }
        h = h.add(&term, p, RoundingMode::None)?;
    }
    Ok(h)
}

/// One even `B_n` from the zeta formula, used when the reverse walk cannot
/// isolate a unique numerator (same fallback as FLINT `_bernoulli_fmpq_ui`).
fn exact_even_one(n: usize) -> Result<(ExactInt, ExactInt), Error> {
    if n < 2 || n % 2 == 1 {
        return Err(Error::InvalidArgument);
    }
    let wp = round_p(bernoulli_global_prec(n)?.saturating_add(256));
    let terms = bernoulli_zeta_terms(n, wp);
    let mut cc = Consts::new()?;
    let mut z = ExactNumNumber::from_word(1, wp)?;
    for j in 2..=terms {
        let pp = bernoulli_power_prec(j, n, wp);
        let jn = pow_ui(&from_usize_p(j, pp)?, n, pp)?;
        let inv = ExactNumNumber::from_word(1, pp)?.div(&jn, pp, RoundingMode::None)?;
        z = z.add(&inv, wp, RoundingMode::None)?;
    }
    let fac = factorial(n, wp)?;
    let two = ExactNumNumber::from_word(2, wp)?;
    let pi = cc.pi_num(wp, RoundingMode::None)?;
    let two_pi = two.mul(&pi, wp, RoundingMode::None)?;
    let pow = pow_ui(&two_pi, n, wp)?;
    let pre = two
        .mul(&fac, wp, RoundingMode::None)?
        .div(&pow, wp, RoundingMode::None)?;
    z = z.mul(&pre, wp, RoundingMode::None)?;
    let primes = primes_upto(n.saturating_add(1));
    let mut denom = bernoulli_denom_ui(n, &primes)?;
    z = z.mul(&denom, wp, RoundingMode::None)?;
    if n.is_multiple_of(4) {
        z = z.neg()?;
    }
    let numer = nearest_unique_int(&z)?;
    denom.set_inexact(false);
    pair_to_ints(numer, denom)
}

/// Exact `(num, den)` for `B_2, …, B_{2n}`. `n = 0` is empty.
pub(super) fn exact_even_bernoulli_rev(n: usize) -> Result<Vec<(ExactInt, ExactInt)>, Error> {
    exact_even_bernoulli_rev_range(1, n)
}

/// Exact `(num, den)` for `B_{2 k_lo}, …, B_{2 k_hi}` in index order.
/// The reverse walk starts at `2 k_hi` and stops once `k < k_lo`, so a
/// short high-end extension does not recompute the already-cached prefix.
pub(super) fn exact_even_bernoulli_rev_range(
    k_lo: usize,
    k_hi: usize,
) -> Result<Vec<(ExactInt, ExactInt)>, Error> {
    if k_hi == 0 || k_lo > k_hi {
        return Ok(Vec::new());
    }
    let n = k_hi;
    let nmax = n.checked_mul(2).ok_or(Error::InvalidArgument)?;
    let slot = k_hi - k_lo + 1;
    let mut out: Vec<Option<(ExactInt, ExactInt)>> = (0..slot).map(|_| None).collect();

    // B_2 .. B_30 from the existing exact triangle (instant at this size).
    if k_lo <= 15 {
        let small_k = n.min(15);
        let small = super::exact_even_bernoulli_boustrophedon(small_k)?;
        for (i, pair) in small.into_iter().enumerate() {
            let k = i + 1;
            if k >= k_lo {
                out[k - k_lo] = Some(pair);
            }
        }
        if n <= 15 {
            return out
                .into_iter()
                .map(|b| b.ok_or(Error::InvalidArgument))
                .collect();
        }
    }

    let mut wp = bernoulli_global_prec(nmax)?;
    let primes = primes_upto(nmax.saturating_add(1));

    let mut max_power = bernoulli_zeta_terms(nmax, wp);
    let mut cc = Consts::new()?;
    let piwp = round_p(wp.saturating_add(2 * bitlen_usize(nmax)).saturating_add(64));

    // Odd inverse powers j^{-nmax} as wp-bit fixed-point integers, plus a running
    // absolute error in the same unit (FLINT `powers` / `pow_error`).
    #[cfg(all(test, feature = "std"))]
    let t_init0 = std::time::Instant::now();
    let one_wp = ExactNumNumber::from_word(1, wp)?;
    let mut powers: Vec<ExactInt> = Vec::with_capacity(max_power + 1);
    for _ in 0..=max_power {
        powers.push(ExactInt::zero());
    }
    let mut pow_error = ExactInt::zero();
    let mut wp_exp = exp_as_exponent(wp)?;
    for j in (3..=max_power).step_by(2) {
        let pp = bernoulli_power_prec(j, nmax, wp);
        let jn = pow_ui(&from_usize_p(j, pp)?, nmax, pp)?;
        let inv = ExactNumNumber::from_word(1, pp)?.div(&jn, pp, RoundingMode::None)?;
        let scaled = inv.ldexp(wp_exp, pp.max(wp), RoundingMode::ToEven)?;
        let tint = nearest_unique_int(&scaled).or_else(|_| {
            // Fall back to a rounded integer; the walk accounts for 1 ulp later.
            let mut half = ExactNumNumber::from_word(1, pp.max(wp))?;
            half = half.ldexp(-1, pp.max(wp), RoundingMode::None)?;
            if scaled.is_positive() || scaled.is_zero() {
                scaled
                    .clone()?
                    .add(&half, pp.max(wp), RoundingMode::None)?
                    .floor()
            } else {
                scaled
                    .clone()?
                    .sub(&half, pp.max(wp), RoundingMode::None)?
                    .ceil()
            }
        })?;
        if !tint.sub(&scaled, wp, RoundingMode::None)?.abs()?.is_zero()
            && pow_error.soft_cmp(&ExactInt::one()) == core::cmp::Ordering::Less
        {
            pow_error = ExactInt::one();
        }
        powers[j] = num_to_int(&tint)?;
    }

    #[cfg(all(test, feature = "std"))]
    let init_pow_ms = t_init0.elapsed().as_millis();
    #[cfg(all(test, feature = "std"))]
    let t_init1 = std::time::Instant::now();
    let fac = factorial(nmax, piwp)?;
    let two = ExactNumNumber::from_word(2, piwp)?;
    let pi = cc.pi_num(piwp, RoundingMode::None)?;
    let two_pi = two.mul(&pi, piwp, RoundingMode::None)?;
    let mut two_pi_squared = two_pi.mul(&two_pi, piwp, RoundingMode::None)?;
    let pow = pow_ui(&two_pi, nmax, piwp)?;
    let mut prefactor =
        two.mul(&fac, piwp, RoundingMode::None)?
            .div(&pow, piwp, RoundingMode::None)?;
    let _ = prefactor.set_precision(wp, RoundingMode::ToEven);
    let _ = two_pi_squared.set_precision(wp, RoundingMode::ToEven);
    #[cfg(all(test, feature = "std"))]
    let init_pref_ms = t_init1.elapsed().as_millis();
    #[cfg(all(test, feature = "std"))]
    if nmax >= 4096 {
        eprintln!(
            "bernoulli_rev init nmax={nmax} wp={wp} terms={max_power} ms_pow={init_pow_ms} ms_pref={init_pref_ms}"
        );
    }

    let mut ncur = nmax;
    #[cfg(all(test, feature = "std"))]
    let mut t_intnum = 0u128;
    #[cfg(all(test, feature = "std"))]
    let mut t_near = 0u128;
    #[cfg(all(test, feature = "std"))]
    let mut t_sum = 0u128;
    #[cfg(all(test, feature = "std"))]
    let mut t_conv = 0u128;
    #[cfg(all(test, feature = "std"))]
    let mut t_denom = 0u128;
    #[cfg(all(test, feature = "std"))]
    let mut t_pref = 0u128;
    #[cfg(all(test, feature = "std"))]
    let mut t_pow = 0u128;
    #[cfg(all(test, feature = "std"))]
    let mut t_reduce = 0u128;
    #[cfg(all(test, feature = "std"))]
    let mut steps = 0u64;

    while ncur >= REV_MIN {
        #[cfg(all(test, feature = "std"))]
        let t0 = std::time::Instant::now();
        let mut sum = ExactInt::zero();
        let mut j = max_power;
        while j >= 3 {
            sum.add_assign(&powers[j]);
            j -= 2;
        }
        #[cfg(all(test, feature = "std"))]
        {
            t_sum += t0.elapsed().as_nanos();
        }
        // z = 1 + sum * 2^{-wp}.
        #[cfg(all(test, feature = "std"))]
        let t01 = std::time::Instant::now();
        let conv = round_p(wp.saturating_add(32));
        let mut z = int_to_num(&sum)?;
        let p_sum = conv.max(z.mantissa_max_bit_len());
        z = z.ldexp(-wp_exp, p_sum, RoundingMode::ToEven)?;
        z = z.add(&one_wp, conv, RoundingMode::None)?;
        #[cfg(all(test, feature = "std"))]
        {
            t_intnum += t01.elapsed().as_nanos();
        }

        #[cfg(all(test, feature = "std"))]
        let t1 = std::time::Instant::now();
        // Convert the odd-zeta sum to ζ(n): z * 2^n / (2^n - 1) (FLINT arb_div_2expm1_ui).
        let h = div_2expm1(&z, ncur, conv)?;
        z = z.add(&h, conv, RoundingMode::None)?;
        z = z.mul(&prefactor, conv, RoundingMode::None)?;
        #[cfg(all(test, feature = "std"))]
        {
            t_conv += t1.elapsed().as_nanos();
        }
        #[cfg(all(test, feature = "std"))]
        let t2 = std::time::Instant::now();
        let denom = bernoulli_denom_ui(ncur, &primes)?;
        z = z.mul(&denom, conv, RoundingMode::None)?;
        #[cfg(all(test, feature = "std"))]
        {
            t_denom += t2.elapsed().as_nanos();
        }
        if ncur.is_multiple_of(4) {
            z = z.neg()?;
        }
        #[cfg(all(test, feature = "std"))]
        let t_n = std::time::Instant::now();
        let pair = match nearest_unique_int(&z) {
            Ok(n) => {
                let mut denom = denom;
                denom.set_inexact(false);
                pair_to_ints(n, denom)?
            }
            Err(_) => exact_even_one(ncur)?,
        };
        #[cfg(all(test, feature = "std"))]
        {
            t_near += t_n.elapsed().as_nanos();
        }

        if ncur % 2 == 0 && ncur >= 2 {
            let k = ncur / 2;
            if k >= k_lo && k <= k_hi {
                out[k - k_lo] = Some(pair);
            }
            if k <= k_lo {
                break;
            }
        }

        if ncur <= 2 {
            break;
        }
        // prefactor_{n-2} = prefactor_n * (2π)^2 / (n(n-1))
        #[cfg(all(test, feature = "std"))]
        let t3 = std::time::Instant::now();
        prefactor = prefactor.mul(&two_pi_squared, wp, RoundingMode::None)?;
        // Divide by a one-word integer (FLINT arb_div_ui). Building n at wp
        // bits made this a 3400-limb / 3400-limb division every step.
        let nn = from_usize_p(ncur, WORD_BIT_SIZE)?;
        let nm1 = from_usize_p(ncur - 1, WORD_BIT_SIZE)?;
        prefactor = prefactor.div(&nn, wp, RoundingMode::None)?;
        prefactor = prefactor.div(&nm1, wp, RoundingMode::None)?;
        #[cfg(all(test, feature = "std"))]
        {
            t_pref += t3.elapsed().as_nanos();
        }

        #[cfg(all(test, feature = "std"))]
        let t4 = std::time::Instant::now();
        let n2 = max_power.saturating_mul(max_power) as u64;
        for j in (3..=max_power).step_by(2) {
            powers[j].mul_u64_assign(j.saturating_mul(j) as u64);
        }
        pow_error.mul_u64_assign(n2);
        #[cfg(all(test, feature = "std"))]
        {
            t_pow += t4.elapsed().as_nanos();
        }

        #[cfg(test)]
        if nmax >= 4096 && ncur.is_multiple_of(512) {
            let done = nmax - ncur;
            let pct = (done as u128) * 100 / (nmax as u128);
            eprintln!(
                "bernoulli_rev phase=walk ncur={ncur} nmax={nmax} counter={done}/{nmax} pct={pct}"
            );
        }

        if ncur.is_multiple_of(64) && ncur > REV_MIN {
            #[cfg(all(test, feature = "std"))]
            let t5 = std::time::Instant::now();
            let new_prec = bernoulli_global_prec(ncur)?;
            let new_max = bernoulli_zeta_terms(ncur, new_prec);
            if new_prec < wp && new_max <= max_power {
                let shift_bits = wp - new_prec;
                for j in (3..=new_max).step_by(2) {
                    powers[j] = powers[j].shr_towards_zero(shift_bits);
                }
                pow_error = pow_error.shr_towards_zero(shift_bits).add(&ExactInt::one());
                let _ = two_pi_squared.set_precision(new_prec, RoundingMode::ToEven);
                let _ = prefactor.set_precision(new_prec, RoundingMode::ToEven);
                max_power = new_max;
                wp = new_prec;
                wp_exp = exp_as_exponent(wp)?;
            }
            #[cfg(all(test, feature = "std"))]
            {
                t_reduce += t5.elapsed().as_nanos();
            }
        }

        #[cfg(all(test, feature = "std"))]
        {
            steps += 1;
        }
        ncur -= 2;
    }

    #[cfg(all(test, feature = "std"))]
    if nmax >= 4096 {
        eprintln!(
            "bernoulli_rev profile nmax={nmax} steps={steps} ms_sum={} ms_intnum={} ms_conv={} ms_denom={} ms_near={} ms_pref={} ms_pow={} ms_reduce={} wp0={} terms0={}",
            t_sum / 1_000_000,
            t_intnum / 1_000_000,
            t_conv / 1_000_000,
            t_denom / 1_000_000,
            t_near / 1_000_000,
            t_pref / 1_000_000,
            t_pow / 1_000_000,
            t_reduce / 1_000_000,
            bernoulli_global_prec(nmax).unwrap_or(0),
            bernoulli_zeta_terms(nmax, bernoulli_global_prec(nmax).unwrap_or(0)),
        );
    }

    out.into_iter()
        .map(|b| b.ok_or(Error::InvalidArgument))
        .collect()
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::defs::Exponent;

    #[test]
    fn flint_sizing_matches_published_formulas() {
        // FLINT 3.6.0 arith_bernoulli_number_size + bernoulli_denom_size
        // + bernoulli_zeta_terms for nmax = 21192 (the EM fill for ζ(1/2) at p=65536).
        let size = bernoulli_number_size_bits(21192);
        assert!(
            (217800..=217820).contains(&size),
            "size bits {size}, expected ~217811"
        );
        let den = bernoulli_denom_size(21192);
        assert!(
            (728..=740).contains(&den),
            "denom size {den}, expected ~729 plus an 8-bit pad"
        );
        let wp = bernoulli_global_prec(21192).unwrap();
        // round_p of ~218539.
        assert!(
            (218500..=218600).contains(&wp),
            "wp {wp}, expected ~218539..218560"
        );
        let terms = bernoulli_zeta_terms(21192, wp);
        assert!(
            (1269..=1275).contains(&terms),
            "zeta terms {terms}, FLINT has 1271"
        );
        // The old integer-ceil power-of-two formula gave 2049 terms here.
        assert!(terms < 1500);
    }

    #[test]
    fn div_2expm1_matches_direct_div() -> Result<(), Error> {
        let p = 256;
        let z = ExactNumNumber::from_word(3, p)?;
        for n in [16usize, 64, 200, 1024] {
            let got = div_2expm1(&z, n, p)?;
            let mut two_n = ExactNumNumber::from_word(1, p)?;
            two_n = two_n.ldexp(exp_as_exponent(n)?, p, RoundingMode::None)?;
            let one = ExactNumNumber::from_word(1, p)?;
            let den = two_n.sub(&one, p, RoundingMode::None)?;
            let expect = z.div(&den, p, RoundingMode::None)?;
            let d = got.sub(&expect, p, RoundingMode::None)?.abs()?;
            if !d.is_zero() {
                assert!(
                    d.exponent() < -((p as Exponent) / 2),
                    "n={n}: exponent {}",
                    d.exponent()
                );
            }
        }
        Ok(())
    }
}
