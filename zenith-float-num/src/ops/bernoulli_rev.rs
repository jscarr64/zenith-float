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
use crate::defs::WORD_BIT_SIZE;
use crate::defs::Sign;
use crate::integer::ExactInt;
use crate::num::ExactNumNumber;
use crate::ops::consts::Consts;
use alloc::vec::Vec;

/// Avoid the reverse walk for very small `n` (FLINT `BERNOULLI_REV_MIN`).
const REV_MIN: usize = 32;

fn bitlen_usize(n: usize) -> usize {
    (usize::BITS - n.leading_zeros()) as usize
}

/// Upper bound on the bits of `|B_n|` (FLINT `arith_bernoulli_number_size`).
fn bernoulli_number_size_bits(n: usize) -> usize {
    if n <= 14 {
        return 8;
    }
    // 2 + (n+1) log2(n+1) - n log2(2πe) + 2, with log2(2πe) > 4.
    let log2n1 = bitlen_usize(n + 1);
    4usize
        .saturating_add((n + 1).saturating_mul(log2n1))
        .saturating_sub(n.saturating_mul(4))
        .max(16)
}

/// `denom(B_n)` by von Staudt–Clausen: `∏_{p-1 | n} p` for even `n`.
fn bernoulli_denom_ui(n: usize) -> Result<ExactNumNumber, Error> {
    if n % 2 == 1 {
        return ExactNumNumber::from_word(u64::from(n == 1) + 1, WORD_BIT_SIZE);
    }
    let primes = primes_upto(n.saturating_add(1));
    let pd = round_p(n.saturating_add(2 * WORD_BIT_SIZE).max(64));
    let mut den = ExactNumNumber::from_word(6, pd)?;
    for p in primes.into_iter().skip(2) {
        let p = p as usize;
        if p.saturating_sub(1) > n {
            break;
        }
        if n % (p - 1) == 0 {
            den = den.mul(&from_usize_p(p, pd)?, pd, RoundingMode::None)?;
        }
    }
    Ok(den)
}

fn bernoulli_global_prec(nmax: usize) -> Result<usize, Error> {
    let size = bernoulli_number_size_bits(nmax);
    let den = bernoulli_denom_ui(nmax)?;
    let dbits = den.mantissa_max_bit_len().max(8);
    Ok(round_p(size.saturating_add(dbits).saturating_add(128)))
}

fn bernoulli_zeta_terms(s: usize, prec: usize) -> usize {
    let den = s.saturating_sub(1).max(1);
    let exp = (prec.saturating_add(1) / den).min((usize::BITS - 2) as usize);
    let mut n = 1usize << exp;
    if n.is_multiple_of(2) {
        n = n.saturating_add(1);
    }
    n.max(3)
}

fn bernoulli_power_prec(j: usize, s: usize, wp: usize) -> usize {
    let log2j = bitlen_usize(j).saturating_sub(1);
    let p = wp.saturating_sub(s.saturating_mul(log2j));
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

fn exp_as_exponent(n: usize) -> Result<Exponent, Error> {
    Exponent::try_from(n).map_err(|_| Error::InvalidArgument)
}

/// Truncate a finite number toward zero as an integer (mantissa is little-endian).
fn num_to_int(x: &ExactNumNumber) -> Result<ExactInt, Error> {
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

/// One even `B_n` from the zeta formula, used when the reverse walk cannot
/// isolate a unique numerator (same fallback as FLINT `_bernoulli_fmpq_ui`).
fn exact_even_one(n: usize) -> Result<(ExactNumNumber, ExactNumNumber), Error> {
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
    let mut denom = bernoulli_denom_ui(n)?;
    z = z.mul(&denom, wp, RoundingMode::None)?;
    if n.is_multiple_of(4) {
        z = z.neg()?;
    }
    let numer = nearest_unique_int(&z)?;
    denom.set_inexact(false);
    Ok((numer, denom))
}

/// Exact `(num, den)` for `B_2, …, B_{2n}`. `n = 0` is empty.
pub(super) fn exact_even_bernoulli_rev(
    n: usize,
) -> Result<Vec<(ExactNumNumber, ExactNumNumber)>, Error> {
    if n == 0 {
        return Ok(Vec::new());
    }
    let nmax = n.checked_mul(2).ok_or(Error::InvalidArgument)?;
    let mut out: Vec<Option<(ExactNumNumber, ExactNumNumber)>> = (0..n).map(|_| None).collect();

    // B_2 .. B_30 from the existing exact triangle (instant at this size).
    let small_k = n.min(15);
    let small = super::exact_even_bernoulli_boustrophedon(small_k)?;
    for (i, pair) in small.into_iter().enumerate() {
        out[i] = Some(pair);
    }
    if n <= 15 {
        return out
            .into_iter()
            .map(|b| b.ok_or(Error::InvalidArgument))
            .collect();
    }

    let mut wp = bernoulli_global_prec(nmax)?;

    let mut max_power = bernoulli_zeta_terms(nmax, wp);
    let mut cc = Consts::new()?;
    let piwp = round_p(wp.saturating_add(2 * bitlen_usize(nmax)).saturating_add(64));

    // Odd inverse powers j^{-nmax} as wp-bit fixed-point integers, plus a running
    // absolute error in the same unit (FLINT `powers` / `pow_error`).
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
                scaled.clone()?.add(&half, pp.max(wp), RoundingMode::None)?.floor()
            } else {
                scaled.clone()?.sub(&half, pp.max(wp), RoundingMode::None)?.ceil()
            }
        })?;
        if !tint.sub(&scaled, wp, RoundingMode::None)?.abs()?.is_zero()
            && pow_error.soft_cmp(&ExactInt::one()) == core::cmp::Ordering::Less
        {
            pow_error = ExactInt::one();
        }
        powers[j] = num_to_int(&tint)?;
    }

    let fac = factorial(nmax, piwp)?;
    let two = ExactNumNumber::from_word(2, piwp)?;
    let pi = cc.pi_num(piwp, RoundingMode::None)?;
    let two_pi = two.mul(&pi, piwp, RoundingMode::None)?;
    let mut two_pi_squared = two_pi.mul(&two_pi, piwp, RoundingMode::None)?;
    let pow = pow_ui(&two_pi, nmax, piwp)?;
    let mut prefactor = two
        .mul(&fac, piwp, RoundingMode::None)?
        .div(&pow, piwp, RoundingMode::None)?;
    let _ = prefactor.set_precision(wp, RoundingMode::ToEven);
    let _ = two_pi_squared.set_precision(wp, RoundingMode::ToEven);

    let mut ncur = nmax;
    while ncur >= REV_MIN {

        let mut sum = ExactInt::zero();
        let mut j = max_power;
        while j >= 3 {
            sum = sum.add(&powers[j]);
            j -= 2;
        }
        // z = 1 + sum * 2^{-wp}.
        let conv = round_p(wp.saturating_add(32));
        let mut z = int_to_num(&sum)?;
        let p_sum = conv.max(z.mantissa_max_bit_len());
        z = z.ldexp(-wp_exp, p_sum, RoundingMode::ToEven)?;
        z = z.add(&one_wp, conv, RoundingMode::None)?;

        // Convert the odd-zeta sum to ζ(n): z * 2^n / (2^n - 1).
        let mut two_n = ExactNumNumber::from_word(1, conv)?;
        two_n = two_n.ldexp(exp_as_exponent(ncur)?, conv, RoundingMode::None)?;
        let den_pow = two_n.sub(&one_wp, conv, RoundingMode::None)?;
        let h = z.div(&den_pow, conv, RoundingMode::None)?;
        z = z.add(&h, conv, RoundingMode::None)?;
        z = z.mul(&prefactor, conv, RoundingMode::None)?;
        let denom = bernoulli_denom_ui(ncur)?;
        z = z.mul(&denom, conv, RoundingMode::None)?;
        if ncur.is_multiple_of(4) {
            z = z.neg()?;
        }
        let (numer, denom) = match nearest_unique_int(&z) {
            Ok(n) => (n, denom),
            Err(_) => exact_even_one(ncur)?,
        };

        if ncur % 2 == 0 && ncur >= 2 {
            let k = ncur / 2;
            if k >= 1 && k <= n {
                let mut denom = denom;
                denom.set_inexact(false);
                out[k - 1] = Some((numer, denom));
            }
        }

        if ncur <= 2 {
            break;
        }
        // prefactor_{n-2} = prefactor_n * (2π)^2 / (n(n-1))
        prefactor = prefactor.mul(&two_pi_squared, wp, RoundingMode::None)?;
        let nn = from_usize_p(ncur, wp)?;
        let nm1 = from_usize_p(ncur - 1, wp)?;
        prefactor = prefactor.div(&nn, wp, RoundingMode::None)?;
        prefactor = prefactor.div(&nm1, wp, RoundingMode::None)?;

        let n2 = ExactInt::from_u64(max_power.saturating_mul(max_power) as u64);
        for j in (3..=max_power).step_by(2) {
            let j2 = ExactInt::from_u64(j.saturating_mul(j) as u64);
            powers[j] = powers[j].mul(&j2);
        }
        pow_error = pow_error.mul(&n2);

        if ncur.is_multiple_of(64) && ncur > REV_MIN {
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
        }

        ncur -= 2;
    }

    out.into_iter()
        .map(|b| b.ok_or(Error::InvalidArgument))
        .collect()
}
