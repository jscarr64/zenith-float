//! Modular exponentiation, inverses, Miller–Rabin, and Pollard–Brent on [`ExactInt`].

use crate::ExactInt;
use alloc::borrow::Cow;
use alloc::vec::Vec;
use core::cmp::Ordering;

/// Former cap on `f` evaluations in one [`pollard_rho`] run. No longer used (ZF-53): a
/// composite always yields a factor, and a prime is recognised by [`is_prime`] first.
/// Kept for API compatibility.
pub const POLLARD_RHO_ITER_MAX: usize = 1_048_576;

/// Brent inner-loop product batch.
const POLLARD_BRENT_M: usize = 128;

pub(crate) fn two() -> ExactInt {
    ExactInt::from_u64(2)
}

/// `|a|`: `a` itself (borrowed) when it is not negative, its negation otherwise. An owned value is
/// [`ExactInt::into_abs`].
pub(crate) fn abs_int(a: &ExactInt) -> Cow<'_, ExactInt> {
    if a.is_negative() {
        Cow::Owned(a.neg())
    } else {
        Cow::Borrowed(a)
    }
}

pub(crate) fn is_even(a: &ExactInt) -> bool {
    a.low_word() & 1 == 0
}

pub(crate) fn rem_nonneg(a: &ExactInt, m: &ExactInt) -> Option<ExactInt> {
    let mabs = abs_int(m);
    if mabs.is_zero() {
        return None;
    }
    let (_, r) = a.div_rem(&mabs)?;
    if r.is_negative() {
        Some(r.add(&mabs))
    } else {
        Some(r)
    }
}

pub(crate) fn mul_mod(a: &ExactInt, b: &ExactInt, m: &ExactInt) -> Option<ExactInt> {
    rem_nonneg(&a.mul(b), m)
}

pub(crate) fn add_mod(a: &ExactInt, b: &ExactInt, m: &ExactInt) -> Option<ExactInt> {
    rem_nonneg(&a.add(b), m)
}

pub(crate) fn sub_mod(a: &ExactInt, b: &ExactInt, m: &ExactInt) -> Option<ExactInt> {
    rem_nonneg(&a.sub(b), m)
}

/// `base^exp mod modulus`. `exp < 0` or a zero modulus is `None`.
/// `|modulus| = 1` is `0`.
pub fn mod_pow(base: &ExactInt, exp: &ExactInt, modulus: &ExactInt) -> Option<ExactInt> {
    if exp.is_negative() {
        return None;
    }
    let m = abs_int(modulus);
    if m.is_zero() {
        return None;
    }
    if m.is_one() {
        return Some(ExactInt::zero());
    }
    let mut acc = ExactInt::one();
    let mut b = rem_nonneg(base, &m)?;
    // `None` is the caller's exponent itself, read through the borrow; every halving is owned.
    let mut e: Option<ExactInt> = None;
    let two = two();
    while !e.as_ref().unwrap_or(exp).is_zero() {
        let er = e.as_ref().unwrap_or(exp);
        if !is_even(er) {
            acc = mul_mod(&acc, &b, &m)?;
        }
        let half = er.div_rem(&two)?.0;
        let more = !half.is_zero();
        e = Some(half);
        if more {
            b = mul_mod(&b, &b, &m)?;
        }
    }
    Some(acc)
}

/// Bezout coefficients `(x, y)` with `a x + b y = 1` when `gcd(a, b) = 1`, `None` otherwise. Only the
/// coefficients are built: the gcd itself is only tested, so neither input is copied.
fn bezout_if_coprime(a: &ExactInt, b: &ExactInt) -> Option<(ExactInt, ExactInt)> {
    if b.is_zero() {
        return a.is_one().then(|| (ExactInt::one(), ExactInt::zero()));
    }
    let (q, r) = a.div_rem(b).unwrap_or((ExactInt::zero(), ExactInt::zero()));
    let (x, y) = bezout_if_coprime(b, &r)?;
    // x1 = y, y1 = x - q y
    let y1 = x.sub(&q.mul(&y));
    Some((y, y1))
}

/// Modular inverse of `a` modulo `modulus`. `None` if not invertible.
pub fn mod_inv(a: &ExactInt, modulus: &ExactInt) -> Option<ExactInt> {
    let m = abs_int(modulus);
    if m.is_zero() || m.is_one() {
        return None;
    }
    let aa = rem_nonneg(a, &m)?;
    if aa.is_zero() {
        return None;
    }
    let (x, _) = bezout_if_coprime(&aa, &m)?;
    rem_nonneg(&x, &m)
}

/// Miller–Rabin on `|n|` with the given bases. `n < 2` is `false`.
///
/// Deterministic for `n < 3·10^{18}` when the witness list is a complete
/// Jaeschke set; this function only uses the callers' bases.
pub fn miller_rabin(n: &ExactInt, witnesses: &[ExactInt]) -> bool {
    let n = abs_int(n);
    if n.soft_cmp(&two()) == Ordering::Less {
        return false;
    }
    if *n == two() || *n == ExactInt::from_u64(3) {
        return true;
    }
    if is_even(&n) {
        return false;
    }
    let one = ExactInt::one();
    let n_minus = n.sub(&one);
    // `d` starts as `n - 1` itself (`None`, read through the borrow); each halving is owned.
    let mut d: Option<ExactInt> = None;
    let mut s = 0usize;
    let two = two();
    loop {
        let dr = d.as_ref().unwrap_or(&n_minus);
        if dr.is_zero() || !is_even(dr) {
            break;
        }
        let half = dr
            .div_rem(&two)
            .map(|(q, _)| q)
            .unwrap_or_else(ExactInt::zero);
        d = Some(half);
        s += 1;
    }
    if s == 0 {
        return false;
    }
    let d = d.as_ref().unwrap_or(&n_minus);
    'wit: for a in witnesses {
        let a = rem_nonneg(a, &n).unwrap_or_else(ExactInt::zero);
        if a.is_zero() || a.is_one() {
            continue;
        }
        let mut x = match mod_pow(&a, d, &n) {
            Some(v) => v,
            None => return false,
        };
        if x.is_one() || x == n_minus {
            continue;
        }
        for _ in 1..s {
            x = match mul_mod(&x, &x, &n) {
                Some(v) => v,
                None => return false,
            };
            if x == n_minus {
                continue 'wit;
            }
            if x.is_one() {
                return false;
            }
        }
        return false;
    }
    !witnesses.is_empty()
}

fn pollard_f(x: &ExactInt, c: &ExactInt, n: &ExactInt) -> Option<ExactInt> {
    add_mod(&mul_mod(x, x, n)?, c, n)
}

fn brent_once(n: &ExactInt, c: &ExactInt) -> Option<ExactInt> {
    let mut y = ExactInt::zero();
    let mut g = ExactInt::one();
    let mut r = 1usize;
    let mut q = ExactInt::one();
    let mut x = ExactInt::zero();
    let mut ys = ExactInt::zero();
    while g.is_one() {
        // `x` takes the current `y` and `y` steps on from it (`r >= 1`): a move, not a copy.
        let y1 = pollard_f(&y, c, n)?;
        x = core::mem::replace(&mut y, y1);
        for _ in 1..r {
            y = pollard_f(&y, c, n)?;
        }
        let mut k = 0usize;
        while k < r && g.is_one() {
            let steps = POLLARD_BRENT_M.min(r - k);
            for step in 0..steps {
                // The first step's start is saved in `ys` (moved out as `y` advances).
                let y_next = pollard_f(&y, c, n)?;
                let y_prev = core::mem::replace(&mut y, y_next);
                if step == 0 {
                    ys = y_prev;
                }
                let diff = x.sub(&y).into_abs();
                q = mul_mod(&q, &diff, n)?;
            }
            g = q.gcd(n);
            k += steps;
        }
        r = r.saturating_mul(2);
    }
    if g == *n {
        loop {
            ys = pollard_f(&ys, c, n)?;
            // gcd reads magnitudes only: no |.| needed
            g = x.sub(&ys).gcd(n);
            if !g.is_one() {
                break;
            }
        }
    }
    if g.is_one() || g == *n {
        None
    } else {
        Some(g)
    }
}

/// Brent Pollard ρ. A proper factor of `|n|`, or `None` when `|n| ≤ 2` or `|n|` is prime.
///
/// Primality is proven first ([`is_prime`]), so `None` never means "gave up": there is no
/// iteration cap, and a composite always yields a factor (ZF-53).
pub fn pollard_rho(n: &ExactInt) -> Option<ExactInt> {
    let n = abs_int(n);
    if n.soft_cmp(&two()) != Ordering::Greater {
        return None;
    }
    if is_even(&n) {
        return Some(two());
    }
    if is_prime(&n) {
        return None;
    }
    Some(split_composite(&n))
}

/// A proper factor of the odd composite `n` (Brent ρ over `c = 1, 2, …`, uncapped). A
/// perfect power or a cycle that closes on `n` moves on to the next `c`; every odd
/// composite has a `c` for which the cycle separates a factor.
fn split_composite(n: &ExactInt) -> ExactInt {
    if let Some(d) = crate::primality::trial_prime_factor(n) {
        if d != *n {
            return d;
        }
    }
    // ECM with growing B1 first. Brent ρ is uncapped, but a single brent_once
    // on a ~55-bit factor is O(2^{27}) ExactInt operations, which is the
    // hang we are replacing. Defer ρ until ECM has tried a large B1; then
    // ρ still runs without a cap if ECM has not found a factor.
    let bits = n.bit_length();
    let mut b1 = if bits > 140 {
        11_000u64
    } else if bits > 90 {
        2_000u64
    } else {
        200u64
    };
    let mut sigma = 6u64;
    let mut c = 1u64;
    loop {
        let curves: u32 = if b1 < 2_000 {
            20
        } else if b1 < 20_000 {
            50
        } else {
            80
        };
        // B2 = 40 B1: stage 2 still helps, without a 100× sieve/walk cost.
        let b2 = b1.saturating_mul(40);
        for _ in 0..curves {
            if let Some(f) = crate::primality::ecm_curve(n, b1, b2, sigma) {
                if !f.is_one() && f != *n {
                    return f;
                }
            }
            sigma = sigma.saturating_add(1).max(6);
        }
        if b1 >= 50_000 {
            if let Some(f) = brent_once(n, &ExactInt::from_u64(c)) {
                if !f.is_one() && f != *n {
                    return f;
                }
            }
            c = c.checked_add(1).unwrap_or(1);
        }
        match b1.checked_mul(3) {
            Some(v) if v > b1 => b1 = v,
            _ => b1 = u64::MAX / 4,
        }
    }
}

/// Small primes for trial division and the deterministic Miller–Rabin set.
const SMALL_PRIMES: [u64; 13] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41];

/// Miller–Rabin with the first 13 primes is deterministic for `n < 3 317 044 064 679 887
/// 385 961 981` (Sorenson & Webster, Math. Comp. 86 (2017)), which covers `2^64`.
const MR13_BOUND: u128 = 3_317_044_064_679_887_385_961_981;

/// Proven primality of `|n|`, like FLINT 3.6 `fmpz_is_prime` (ZF-53).
///
/// - `|n| < 3.3·10^24` (past `2^64`): Miller–Rabin with the first 13 prime bases, which is
///   deterministic there (Sorenson & Webster).
/// - Larger: Baillie–PSW (Miller–Rabin base 2 plus a strong Lucas–Selfridge test)
///   rejects composites. A probable-prime result is never returned as a proof. A
///   Brillhart–Lehmer–Selfridge `n-1` / `n+1` test with only a partial
///   factorization of `n±1` then proves primality when the factored part is
///   large enough. If that front end does not decide, the Cohen–Lenstra
///   Jacobi-sum APRCL test does. Neither step ECM-factors `n-1`.
///
/// No step is capped. The answer is a proof either way.
pub fn is_prime(n: &ExactInt) -> bool {
    let n = abs_int(n);
    if n.soft_cmp(&two()) == Ordering::Less {
        return false;
    }
    if let Some(d) = crate::primality::trial_prime_factor(&n) {
        return *n == d;
    }
    if n.soft_cmp(&ExactInt::from_u128(MR13_BOUND)) == Ordering::Less {
        let wits: Vec<ExactInt> = SMALL_PRIMES
            .iter()
            .map(|&q| ExactInt::from_u64(q))
            .collect();
        return miller_rabin(&n, &wits);
    }
    // Filter only. A BPSW probable prime is not a proof.
    if !crate::primality::bpsw(&n) {
        return false;
    }
    crate::aprcl::prove_prime(&n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_must::Must;

    #[test]
    fn modular_pow_inv_miller_rho() -> Result<(), crate::Error> {
        let two = ExactInt::from_u64(2);
        let exp = ExactInt::from_u64(100);
        let m = ExactInt::from_u64(1_000_000_007);
        let got = mod_pow(&two, &exp, &m).ok_or(crate::Error::InvalidArgument)?;
        assert_eq!(got, ExactInt::from_u64(976_371_285));

        let inv = mod_inv(&ExactInt::from_u64(3), &ExactInt::from_u64(7))
            .ok_or(crate::Error::InvalidArgument)?;
        assert_eq!(inv, ExactInt::from_u64(5));
        assert!(mod_inv(&ExactInt::from_u64(2), &ExactInt::from_u64(4)).is_none());

        let mersenne = two.pow(31).sub(&ExactInt::one());
        let wits: Vec<ExactInt> = [2u64, 3, 5, 7]
            .into_iter()
            .map(ExactInt::from_u64)
            .collect();
        assert!(miller_rabin(&mersenne, &wits));
        assert!(!miller_rabin(&ExactInt::from_u64(9), &wits));

        let n = ExactInt::from_u64(8051);
        let f = pollard_rho(&n).ok_or(crate::Error::InvalidArgument)?;
        let other = n.div_rem(&f).ok_or(crate::Error::InvalidArgument)?.0;
        let a = ExactInt::from_u64(83);
        let b = ExactInt::from_u64(97);
        assert!(
            (f == a && other == b) || (f == b && other == a),
            "factor {f:?}"
        );
        assert!(pollard_rho(&ExactInt::from_u64(31)).is_none());
        assert!(mod_pow(&two, &ExactInt::from_i64(-1), &m).is_none());
        Ok(())
    }

    /// ZF-53 / d0225257: a FLINT 3.6-style proof path (BPSW filter,
    /// BLS n−1/n+1 with a partial factorization, then Jacobi-sum APRCL).
    /// Probable-prime filters are never returned as a proof. Covers
    /// 2^255−19, the ZF-53 Mersennes, Carmichael numbers, strong
    /// pseudoprimes to several bases, and semiprimes with large
    /// balanced factors. 2^255−19 must finish under 2 s native: 1.5 × the
    /// slowest run measured, 1.32 s (i9-14900KF E-core beside 3 other test
    /// slices; 0.76 s on a P-core; FLINT ~14 ms). Measurements:
    /// `~/pop-logs/paul/2026-11-15/land39/final/pr40-fixes/r2nd/pollard/limit/LIMIT.md`.
    #[test]
    fn zf_53_pollard_rho_proof_path() -> Result<(), String> {
        fn dec(s: &str) -> Result<ExactInt, String> {
            ExactInt::from_dec_string(s).must("test integer")
        }
        fn assert_proved_prime(n: &ExactInt) {
            // pollard_rho proves primality via is_prime and returns None.
            // Do not also call is_prime first: the certificate is not cached,
            // and a second proof would double the wall time.
            assert!(
                pollard_rho(n).is_none(),
                "must be a proved prime (pollard_rho yields None) {}",
                n.to_dec_string()
            );
        }
        fn assert_proved_composite(n: &ExactInt) -> Result<(), String> {
            assert!(
                !is_prime(n),
                "must be proved composite {}",
                n.to_dec_string()
            );
            let f = pollard_rho(n).must("a composite must yield a proper factor")?;
            assert!(
                f.soft_cmp(&ExactInt::one()) == Ordering::Greater
                    && f.soft_cmp(n) == Ordering::Less,
                "factor {} is not proper for {}",
                f.to_dec_string(),
                n.to_dec_string()
            );
            let (q, r) = n.div_rem(&f).must("div")?;
            assert!(r.is_zero(), "factor does not divide");
            assert!(!q.is_one() && q != *n);
            Ok(())
        }

        // Curve25519 prime, d0225257: APRCL, not an ECM/ρ search of n−1.
        let p255 = ExactInt::from_u64(2).pow(255).sub(&ExactInt::from_u64(19));
        // 2 s native (1.5 × the 1.32 s worst measured run); times ZF_TIME_FACTOR on an emulated
        // runner (crate::test_time).
        let limit = crate::test_time::budget(std::time::Duration::from_secs(2));
        assert!(limit.is_ok(), "{limit:?}");
        let limit = limit.unwrap_or_default();
        let t0 = std::time::Instant::now();
        assert_proved_prime(&p255);
        let took = t0.elapsed();
        assert!(
            crate::test_time::limit_holds(took < limit),
            "2^255-19 must finish under 2 s native ({limit:?} here); took {took:?}"
        );

        // ZF-53 Mersennes (d0225248, d0225251).
        assert_proved_prime(&ExactInt::from_u64(2).pow(89).sub(&ExactInt::one()));
        assert_proved_prime(&ExactInt::from_u64(2).pow(127).sub(&ExactInt::one()));

        // Carmichael numbers (Fermat PSP to every coprime base).
        for s in [
            "561",
            "1105",
            "1729",
            "2465",
            "2821",
            "6601",
            "8911",
            "10585",
            "15841",
            "29341",
            "172081",
            "2152302898747",
            "105098405405550221251",
            "3211971524423007740460469",
            // Chernick U_k, k=13679106, above the MR13 deterministic bound.
            "3317249643051242788534009",
        ] {
            assert_proved_composite(&dec(s)?)?;
        }

        // Strong pseudoprimes to the first several prime bases.
        for s in [
            "1373653", // strong PSP to 2 and 3
            "25326001",
            "161304001",
            "960946321",
            "1156916633",
            "3215031751", // strong PSP to 2, 3, 5, and 7
            "75130373920513",
        ] {
            assert_proved_composite(&dec(s)?)?;
        }

        // Semiprimes with large balanced factors.
        let p40a = dec("1099511627791")?;
        let p40b = dec("1099511627803")?;
        assert_proved_composite(&p40a.mul(&p40b))?;
        let p48a = dec("281474976710677")?;
        let p48b = dec("281474976710717")?;
        assert_proved_composite(&p48a.mul(&p48b))?;
        Ok(())
    }
}
