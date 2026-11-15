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

fn two() -> ExactInt {
    ExactInt::from_u64(2)
}

fn abs_int(a: &ExactInt) -> ExactInt {
    if a.is_negative() {
        a.neg()
    } else {
        a.clone()
    }
}

fn is_even(a: &ExactInt) -> bool {
    a.low_word() & 1 == 0
}

fn rem_nonneg(a: &ExactInt, m: &ExactInt) -> Option<ExactInt> {
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

fn mul_mod(a: &ExactInt, b: &ExactInt, m: &ExactInt) -> Option<ExactInt> {
    rem_nonneg(&a.mul(b), m)
}

fn add_mod(a: &ExactInt, b: &ExactInt, m: &ExactInt) -> Option<ExactInt> {
    rem_nonneg(&a.add(b), m)
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
    let mut c = 1u64;
    loop {
        if let Some(f) = brent_once(n, &ExactInt::from_u64(c)) {
            if !f.is_one() && f != *n {
                return f;
            }
        }
        // A cycle that never separates for every c below 2^64 would need n's prime
        // factors to share every orbit; trial division settles it before c wraps.
        c = match c.checked_add(1) {
            Some(v) => v,
            None => return trial_factor(n),
        };
    }
}

/// Smallest prime factor of the composite `n` by trial division (the last resort of
/// [`split_composite`]; never reached in practice).
fn trial_factor(n: &ExactInt) -> ExactInt {
    let mut d = ExactInt::from_u64(3);
    let two = two();
    loop {
        if rem_nonneg(n, &d).is_some_and(|r| r.is_zero()) {
            return d;
        }
        d = d.add(&two);
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
///   deterministic there.
/// - Larger: a strong probable-prime test with those bases rejects composites, then a
///   Pocklington–Lehmer `n − 1` certificate proves primality. `n − 1` is factored (trial
///   division, then [`is_prime`] and Brent ρ on the cofactors, recursively) until the
///   proven part `F` satisfies `F² > n`. Each prime `q | F` then needs a witness `a` with
///   `a^(n−1) ≡ 1` and `gcd(a^((n−1)/q) − 1, n) = 1`. Witnesses are tried from `a = 2` up,
///   and an `a` sharing a factor with `n` proves it composite, so the search terminates.
///
/// No step is capped. The answer is a proof either way.
pub fn is_prime(n: &ExactInt) -> bool {
    let n = abs_int(n);
    if n.soft_cmp(&two()) == Ordering::Less {
        return false;
    }
    for &q in SMALL_PRIMES.iter() {
        let qi = ExactInt::from_u64(q);
        if n == qi {
            return true;
        }
        if rem_nonneg(&n, &qi).is_some_and(|r| r.is_zero()) {
            return false;
        }
    }
    let wits: Vec<ExactInt> = SMALL_PRIMES
        .iter()
        .map(|&q| ExactInt::from_u64(q))
        .collect();
    if !miller_rabin(&n, &wits) {
        return false;
    }
    if n.soft_cmp(&ExactInt::from_u128(MR13_BOUND)) == Ordering::Less {
        return true;
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
    for q in &certified_part(&nm1, n) {
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
}
