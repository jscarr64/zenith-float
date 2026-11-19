//! Modular exponentiation, inverses, Miller–Rabin, and Pollard–Brent on [`ExactInt`].

use crate::ExactInt;
use alloc::vec;
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

pub(crate) fn abs_int(a: &ExactInt) -> ExactInt {
    if a.is_negative() {
        a.neg()
    } else {
        a.clone()
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
    let mut e = exp.clone();
    let two = two();
    while !e.is_zero() {
        if !is_even(&e) {
            acc = mul_mod(&acc, &b, &m)?;
        }
        e = e.div_rem(&two)?.0;
        if !e.is_zero() {
            b = mul_mod(&b, &b, &m)?;
        }
    }
    Some(acc)
}

fn egcd(a: &ExactInt, b: &ExactInt) -> (ExactInt, ExactInt, ExactInt) {
    if b.is_zero() {
        return (a.clone(), ExactInt::one(), ExactInt::zero());
    }
    let (q, r) = a.div_rem(b).unwrap_or((ExactInt::zero(), ExactInt::zero()));
    let (g, x, y) = egcd(b, &r);
    // x1 = y, y1 = x - q y
    (g, y.clone(), x.sub(&q.mul(&y)))
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
    let (g, x, _) = egcd(&aa, &m);
    if !g.is_one() {
        return None;
    }
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
    if n == two() || n == ExactInt::from_u64(3) {
        return true;
    }
    if is_even(&n) {
        return false;
    }
    let one = ExactInt::one();
    let n_minus = n.sub(&one);
    let mut d = n_minus.clone();
    let mut s = 0usize;
    let two = two();
    while !d.is_zero() && is_even(&d) {
        d = d
            .div_rem(&two)
            .map(|(q, _)| q)
            .unwrap_or_else(ExactInt::zero);
        s += 1;
    }
    if s == 0 {
        return false;
    }
    'wit: for a in witnesses {
        let a = rem_nonneg(a, &n).unwrap_or_else(ExactInt::zero);
        if a.is_zero() || a.is_one() {
            continue;
        }
        let mut x = match mod_pow(&a, &d, &n) {
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
        x = y.clone();
        for _ in 0..r {
            y = pollard_f(&y, c, n)?;
        }
        let mut k = 0usize;
        while k < r && g.is_one() {
            ys = y.clone();
            let steps = POLLARD_BRENT_M.min(r - k);
            for _ in 0..steps {
                y = pollard_f(&y, c, n)?;
                let diff = abs_int(&x.sub(&y));
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
            g = abs_int(&x.sub(&ys)).gcd(n);
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
            c = match c.checked_add(1) {
                Some(v) => v,
                None => 1,
            };
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

/// Proven primality of `|n|`, like FLINT's `fmpz_is_prime` (ZF-53).
///
/// - `|n| < 3.3·10^24` (past `2^64`): Miller–Rabin with the first 13 prime bases, which is
///   deterministic there (Sorenson & Webster).
/// - Larger: Baillie–PSW (Miller–Rabin base 2 plus a strong Lucas–Selfridge test)
///   rejects composites. A probable-prime result is never returned as a proof. A
///   Pocklington–Lehmer `n − 1` certificate then proves primality. `n − 1` is factored
///   (trial division, then [`is_prime`] and ECM or Brent ρ on the cofactors,
///   recursively) until the proven part `F` satisfies `F² > n`. Each prime `q | F`
///   then needs a witness `a` with `a^(n−1) ≡ 1` and `gcd(a^((n−1)/q) − 1, n) = 1`.
///   Witnesses are tried from `a = 2` up, and an `a` sharing a factor with `n` proves
///   it composite, so the search terminates.
///
/// No step is capped. The answer is a proof either way.
pub fn is_prime(n: &ExactInt) -> bool {
    let n = abs_int(n);
    if n.soft_cmp(&two()) == Ordering::Less {
        return false;
    }
    if let Some(d) = crate::primality::trial_prime_factor(&n) {
        return n == d;
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
    pocklington(&n)
}

/// Distinct proven prime factors `q` of `n − 1` whose product with multiplicity, `F`,
/// divides `n − 1` and satisfies `F² > n` (the Pocklington bound). Factoring stops as soon
/// as the bound holds, so a hard cofactor of `n − 1` left over is never split.
fn certified_part(nm1: &ExactInt, n: &ExactInt) -> Vec<ExactInt> {
    let mut f = ExactInt::one();
    let mut qs: Vec<ExactInt> = Vec::new();
    let done = |f: &ExactInt| f.mul(f).soft_cmp(n) == Ordering::Greater;
    let mut m = nm1.clone();
    for &q in SMALL_PRIMES.iter() {
        let qi = ExactInt::from_u64(q);
        while let Some((d, r)) = m.div_rem(&qi) {
            if !r.is_zero() {
                break;
            }
            if !qs.contains(&qi) {
                qs.push(qi.clone());
            }
            f = f.mul(&qi);
            m = d;
        }
    }
    crate::primality::divide_out_primes_upto(&mut m, &mut f, &mut qs, crate::primality::TRIAL_LIMIT);
    let mut stack = vec![m];
    while let Some(x) = stack.pop() {
        if done(&f) {
            break;
        }
        if x.is_one() || x.is_zero() {
            continue;
        }
        if is_prime(&x) {
            f = f.mul(&x);
            if !qs.contains(&x) {
                qs.push(x);
            }
            continue;
        }
        let d = if is_even(&x) { two() } else { split_composite(&x) };
        let g = x.div_rem(&d).map(|(q, _)| q).unwrap_or_else(ExactInt::one);
        // The smaller part first: it is cheaper to certify.
        if d.soft_cmp(&g) == Ordering::Less {
            stack.push(g);
            stack.push(d);
        } else {
            stack.push(d);
            stack.push(g);
        }
    }
    qs
}

/// Pocklington–Lehmer proof for an odd strong probable prime `n` (see [`is_prime`]):
/// `F | n − 1`, `F² > n`, and for each prime `q | F` some `a` with `a^(n−1) ≡ 1 (mod n)`
/// and `gcd(a^((n−1)/q) − 1, n) = 1` (Crandall & Pomerance, Thm 4.1.3).
fn pocklington(n: &ExactInt) -> bool {
    let one = ExactInt::one();
    let nm1 = n.sub(&one);
    let qs = certified_part(&nm1, n);
    let mut f = ExactInt::one();
    for q in &qs {
        // Reconstruct F from the distinct primes with the multiplicity in n−1.
        let mut tmp = nm1.clone();
        while let Some((d, r)) = tmp.div_rem(q) {
            if !r.is_zero() {
                break;
            }
            f = f.mul(q);
            tmp = d;
        }
    }
    if f.mul(&f).soft_cmp(n) != Ordering::Greater {
        // Incomplete n−1 factorization is not a proof. certified_part keeps
        // splitting until the bound holds, so this path is defensive.
        return false;
    }
    for q in &qs {
        let e = match nm1.div_rem(q) {
            Some((e, _)) => e,
            None => return false,
        };
        let mut a = two();
        loop {
            // An `a` sharing a factor with n proves n composite; it also bounds the search
            // by n's least prime factor.
            if !a.gcd(n).is_one() {
                return false;
            }
            match mod_pow(&a, &nm1, n) {
                Some(x) if x.is_one() => {}
                _ => return false,
            }
            let t = match mod_pow(&a, &e, n) {
                Some(t) => t,
                None => return false,
            };
            let h = abs_int(&t.sub(&one)).gcd(n);
            if h.is_one() {
                break;
            }
            if h != *n {
                return false;
            }
            a = a.add(&one);
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

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

    /// ZF-53 / d0225257: a FLINT-style proof path (BPSW, then a
    /// Pocklington–Lehmer n−1 certificate whose factoring uses ECM).
    /// Probable-prime filters are never returned as a proof. Covers
    /// 2^255−19, the ZF-53 Mersennes, Carmichael numbers, strong
    /// pseudoprimes to several bases, and semiprimes with large
    /// balanced factors.
    #[test]
    fn zf_53_pollard_rho_proof_path() -> Result<(), crate::Error> {
        fn dec(s: &str) -> ExactInt {
            ExactInt::from_dec_string(s).expect("test integer")
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
        fn assert_proved_composite(n: &ExactInt) {
            assert!(
                !is_prime(n),
                "must be proved composite {}",
                n.to_dec_string()
            );
            let f = pollard_rho(n).expect("a composite must yield a proper factor");
            assert!(
                f.soft_cmp(&ExactInt::one()) == Ordering::Greater
                    && f.soft_cmp(n) == Ordering::Less,
                "factor {} is not proper for {}",
                f.to_dec_string(),
                n.to_dec_string()
            );
            let (q, r) = n.div_rem(&f).expect("div");
            assert!(r.is_zero(), "factor does not divide");
            assert!(!q.is_one() && q != *n);
        }

        // Curve25519 prime, d0225257: needs ECM in the n−1 proof, not ρ.
        let p255 = ExactInt::from_u64(2).pow(255).sub(&ExactInt::from_u64(19));
        let t0 = std::time::Instant::now();
        assert_proved_prime(&p255);
        assert!(
            t0.elapsed().as_secs() < 110,
            "2^255-19 must finish under the 120 s wall, not a ρ search; took {:?}"
            , t0.elapsed()
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
            assert_proved_composite(&dec(s));
        }

        // Strong pseudoprimes to the first several prime bases.
        for s in [
            "1373653",    // strong PSP to 2 and 3
            "25326001",
            "161304001",
            "960946321",
            "1156916633",
            "3215031751", // strong PSP to 2, 3, 5, and 7
            "75130373920513",
        ] {
            assert_proved_composite(&dec(s));
        }

        // Semiprimes with large balanced factors.
        let p40a = dec("1099511627791");
        let p40b = dec("1099511627803");
        assert_proved_composite(&p40a.mul(&p40b));
        let p48a = dec("281474976710677");
        let p48b = dec("281474976710717");
        assert_proved_composite(&p48a.mul(&p48b));
        Ok(())
    }
}
