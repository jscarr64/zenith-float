//! Modular exponentiation, inverses, Miller–Rabin, and Pollard–Brent on [`ExactInt`].

use crate::ExactInt;
use core::cmp::Ordering;

/// Maximum `f` evaluations in one [`pollard_rho`] run before `None`.
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
    if n.cmp(&two()) == Ordering::Less {
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

fn brent_once(n: &ExactInt, c: &ExactInt, max_f: usize) -> Option<ExactInt> {
    let mut y = ExactInt::zero();
    let mut g = ExactInt::one();
    let mut r = 1usize;
    let mut q = ExactInt::one();
    let mut f_used = 0usize;
    let mut x = ExactInt::zero();
    let mut ys = ExactInt::zero();
    while g.is_one() {
        x = y.clone();
        for _ in 0..r {
            y = pollard_f(&y, c, n)?;
            f_used += 1;
            if f_used >= max_f {
                return None;
            }
        }
        let mut k = 0usize;
        while k < r && g.is_one() {
            ys = y.clone();
            let steps = POLLARD_BRENT_M.min(r - k);
            for _ in 0..steps {
                y = pollard_f(&y, c, n)?;
                f_used += 1;
                if f_used >= max_f {
                    return None;
                }
                let diff = abs_int(&x.sub(&y));
                q = mul_mod(&q, &diff, n)?;
            }
            g = q.gcd(n);
            k += steps;
        }
        if r > usize::MAX / 2 {
            return None;
        }
        r = r.saturating_mul(2);
    }
    if g == *n {
        loop {
            ys = pollard_f(&ys, c, n)?;
            f_used += 1;
            if f_used >= max_f {
                return None;
            }
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

/// Brent Pollard ρ. A proper factor of `|n|`, or `None` if prime / cap hit.
pub fn pollard_rho(n: &ExactInt) -> Option<ExactInt> {
    let n = abs_int(n);
    if n.cmp(&two()) != Ordering::Greater {
        return None;
    }
    if is_even(&n) {
        return Some(two());
    }
    for c in 1u64..=32 {
        if let Some(f) = brent_once(&n, &ExactInt::from_u64(c), POLLARD_RHO_ITER_MAX) {
            if !f.is_one() && f != n {
                return Some(f);
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn modular_pow_inv_miller_rho() {
        let two = ExactInt::from_u64(2);
        let exp = ExactInt::from_u64(100);
        let m = ExactInt::from_u64(1_000_000_007);
        let got = mod_pow(&two, &exp, &m).expect("mod_pow");
        assert_eq!(got, ExactInt::from_u64(976_371_285));

        let inv = mod_inv(&ExactInt::from_u64(3), &ExactInt::from_u64(7)).expect("inv");
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
        let f = pollard_rho(&n).expect("rho");
        let other = n.div_rem(&f).expect("div").0;
        let a = ExactInt::from_u64(83);
        let b = ExactInt::from_u64(97);
        assert!(
            (f == a && other == b) || (f == b && other == a),
            "factor {f:?}"
        );
        assert!(pollard_rho(&ExactInt::from_u64(31)).is_none());
        assert!(mod_pow(&two, &ExactInt::from_i64(-1), &m).is_none());
    }
}
