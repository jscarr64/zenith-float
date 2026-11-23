//! Baillie–PSW (as a filter only) and Montgomery ECM for the Pocklington
//! proof path in [`crate::modular::is_prime`]. A probable-prime result is
//! never returned as a proof.

use crate::modular::{abs_int, add_mod, is_even, mod_inv, mul_mod, rem_nonneg, sub_mod, two};
use crate::ExactInt;
use alloc::vec;
use alloc::vec::Vec;
use core::cmp::Ordering;

/// Trial-division limit used before BPSW, ECM, or Brent ρ. This is a cheap
/// first pass, not a cap: a larger prime factor is found by ECM or ρ.
pub(crate) const TRIAL_LIMIT: u64 = 2_000_003;

pub(crate) fn trial_prime_factor(n: &ExactInt) -> Option<ExactInt> {
    if n.soft_cmp(&two()) == Ordering::Less {
        return None;
    }
    if is_even(n) {
        return Some(two());
    }
    for p in primes_upto(TRIAL_LIMIT).into_iter() {
        let q = ExactInt::from_u64(u64::from(p));
        if n.soft_cmp(&q) != Ordering::Greater {
            return Some(n.clone());
        }
        if rem_nonneg(n, &q).is_some_and(|r| r.is_zero()) {
            return Some(q);
        }
    }
    None
}

/// Baillie–PSW: Miller–Rabin with base 2, then a strong Lucas–Selfridge test.
/// This is a composite filter only, not a primality proof.
pub(crate) fn bpsw(n: &ExactInt) -> bool {
    let n = abs_int(n);
    if n.soft_cmp(&two()) == Ordering::Less {
        return false;
    }
    if n == two() {
        return true;
    }
    if is_even(&n) {
        return false;
    }
    let two_w = two();
    if !crate::modular::miller_rabin(&n, &[two_w]) {
        return false;
    }
    strong_lucas_selfridge(&n)
}

pub(crate) fn jacobi(mut a: ExactInt, mut n: ExactInt) -> i8 {
    if n.is_zero() || is_even(&n) || n.is_negative() {
        return 0;
    }
    a = rem_nonneg(&a, &n).unwrap_or_else(ExactInt::zero);
    let mut t: i8 = 1;
    while !a.is_zero() {
        while is_even(&a) {
            a = a
                .div_rem(&two())
                .map(|(q, _)| q)
                .unwrap_or_else(ExactInt::zero);
            let n8 = n.low_word() & 7;
            if n8 == 3 || n8 == 5 {
                t = -t;
            }
        }
        // quadratic reciprocity
        if a.low_word() & 3 == 3 && n.low_word() & 3 == 3 {
            t = -t;
        }
        let r = rem_nonneg(&n, &a).unwrap_or_else(ExactInt::zero);
        n = a;
        a = r;
    }
    if n.is_one() {
        t
    } else {
        0
    }
}

/// Selfridge's method: the first `D` in 5, −7, 9, −11, … with `(D/n) = −1`.
/// `P = 1`, `Q = (1−D)/4`. Returns `None` if a proper factor of `n` is found.
fn selfridge_params(n: &ExactInt) -> Option<(ExactInt, ExactInt, ExactInt)> {
    let mut d = ExactInt::from_i64(5);
    let mut sign = 1i8;
    loop {
        let dd = if sign < 0 { d.neg() } else { d.clone() };
        let j = jacobi(dd.clone(), n.clone());
        if j == 0 {
            let g = abs_int(&dd).gcd(n);
            if g.is_one() || g == *n {
                // |D| = n, try the next D
            } else {
                return None;
            }
        } else if j == -1 {
            // Q = (1-D)/4
            let one = ExactInt::one();
            let four = ExactInt::from_u64(4);
            let num = one.sub(&dd);
            let q = num.div_rem(&four).map(|(q, r)| {
                debug_assert!(r.is_zero());
                q
            })?;
            return Some((dd, ExactInt::one(), q));
        }
        d = d.add(&ExactInt::from_u64(2));
        sign = -sign;
    }
}

/// Strong Lucas–Selfridge probable-prime test (the Lucas half of BPSW).
fn strong_lucas_selfridge(n: &ExactInt) -> bool {
    let Some((d, p, q)) = selfridge_params(n) else {
        return false;
    };
    let two = two();
    // n+1 = 2^s * d_odd
    let mut delta = n.add(&ExactInt::one());
    let mut s = 0usize;
    while is_even(&delta) {
        delta = delta
            .div_rem(&two)
            .map(|(q, _)| q)
            .unwrap_or_else(ExactInt::zero);
        s += 1;
    }
    if s == 0 {
        return false;
    }
    let Some((u, v0, qk0)) = lucas_uv(n, &p, &q, &d, &delta) else {
        return false;
    };
    if u.is_zero() || v0.is_zero() {
        return true;
    }
    let mut vk = v0;
    let mut qk = qk0;
    for _ in 1..s {
        vk = match mul_mod(&vk, &vk, n) {
            Some(x) => x,
            None => return false,
        };
        let two_qk = match mul_mod(&two, &qk, n) {
            Some(x) => x,
            None => return false,
        };
        vk = match sub_mod(&vk, &two_qk, n) {
            Some(x) => x,
            None => return false,
        };
        if vk.is_zero() {
            return true;
        }
        qk = match mul_mod(&qk, &qk, n) {
            Some(x) => x,
            None => return false,
        };
    }
    false
}

pub(crate) fn lucas_uv(
    n: &ExactInt,
    p: &ExactInt,
    q: &ExactInt,
    d: &ExactInt,
    k: &ExactInt,
) -> Option<(ExactInt, ExactInt, ExactInt)> {
    if k.is_zero() {
        return Some((ExactInt::zero(), ExactInt::from_u64(2), ExactInt::one()));
    }
    let one = ExactInt::one();
    let two = two();
    let inv2 = mod_inv(&two, n)?;
    // Start at k's highest bit: U=1, V=P, Qk=Q after processing the MSB (which is 1).
    let bits = k.bit_length();
    let mut u = one.clone();
    let mut v = rem_nonneg(p, n)?;
    let mut qk = rem_nonneg(q, n)?;
    let dm = rem_nonneg(d, n)?;
    // Process bits from bits-2 down to 0 (MSB already consumed as the start).
    for i in (0..bits.saturating_sub(1)).rev() {
        // Double: U2=U V, V2=V^2-2 Qk, Qk := Qk^2
        let u2 = mul_mod(&u, &v, n)?;
        let vsq = mul_mod(&v, &v, n)?;
        let two_qk = mul_mod(&two, &qk, n)?;
        let v2 = sub_mod(&vsq, &two_qk, n)?;
        qk = mul_mod(&qk, &qk, n)?;
        u = u2;
        v = v2;
        if bit_set(k, i) {
            // Add one: U_{k+1} = (P U + V)/2, V_{k+1} = (D U + P V)/2, Qk *= Q
            let pu = mul_mod(p, &u, n)?;
            let un = mul_mod(&add_mod(&pu, &v, n)?, &inv2, n)?;
            let du = mul_mod(&dm, &u, n)?;
            let pv = mul_mod(p, &v, n)?;
            let vn = mul_mod(&add_mod(&du, &pv, n)?, &inv2, n)?;
            u = un;
            v = vn;
            qk = mul_mod(&qk, q, n)?;
        }
    }
    Some((u, v, qk))
}

fn bit_set(k: &ExactInt, i: usize) -> bool {
    k.shr_towards_zero(i).low_word() & 1 == 1
}

#[derive(Clone)]
struct MPt {
    x: ExactInt,
    z: ExactInt,
}

fn mont_dbl(p: &MPt, a24: &ExactInt, n: &ExactInt) -> Option<MPt> {
    let xpz = add_mod(&p.x, &p.z, n)?;
    let xmz = sub_mod(&p.x, &p.z, n)?;
    let s = mul_mod(&xpz, &xpz, n)?;
    let t = mul_mod(&xmz, &xmz, n)?;
    let m = sub_mod(&s, &t, n)?;
    let x2 = mul_mod(&s, &t, n)?;
    let a24m = mul_mod(a24, &m, n)?;
    let z2 = mul_mod(&m, &add_mod(&t, &a24m, n)?, n)?;
    Some(MPt { x: x2, z: z2 })
}

fn mont_add(p: &MPt, q: &MPt, diff: &MPt, n: &ExactInt) -> Option<MPt> {
    let u = mul_mod(&sub_mod(&p.x, &p.z, n)?, &add_mod(&q.x, &q.z, n)?, n)?;
    let v = mul_mod(&add_mod(&p.x, &p.z, n)?, &sub_mod(&q.x, &q.z, n)?, n)?;
    let upv = add_mod(&u, &v, n)?;
    let umv = sub_mod(&u, &v, n)?;
    let x = mul_mod(&diff.z, &mul_mod(&upv, &upv, n)?, n)?;
    let z = mul_mod(&diff.x, &mul_mod(&umv, &umv, n)?, n)?;
    Some(MPt { x, z })
}

fn mont_mul(p: &MPt, k: &ExactInt, a24: &ExactInt, n: &ExactInt) -> Option<MPt> {
    if k.is_zero() {
        return Some(MPt {
            x: ExactInt::one(),
            z: ExactInt::zero(),
        });
    }
    if k.is_one() {
        return Some(p.clone());
    }
    let bits = k.bit_length();
    let mut r0 = p.clone();
    let mut r1 = mont_dbl(p, a24, n)?;
    for i in (0..bits.saturating_sub(1)).rev() {
        if bit_set(k, i) {
            r0 = mont_add(&r0, &r1, p, n)?;
            r1 = mont_dbl(&r1, a24, n)?;
        } else {
            r1 = mont_add(&r0, &r1, p, n)?;
            r0 = mont_dbl(&r0, a24, n)?;
        }
    }
    Some(r0)
}

fn as_factor(g: ExactInt, n: &ExactInt) -> Option<ExactInt> {
    let g = abs_int(&g);
    if g.is_one() || g == *n || g.is_zero() {
        None
    } else {
        Some(g)
    }
}

/// One Suyama curve, Montgomery stage 1 and stage 2. `None` if this curve
/// finds no proper factor.
pub(crate) fn ecm_curve(n: &ExactInt, b1: u64, b2: u64, sigma: u64) -> Option<ExactInt> {
    if b1 < 2 || n.soft_cmp(&ExactInt::from_u64(3)) != Ordering::Greater {
        return None;
    }
    let sig = ExactInt::from_u64(sigma.max(6));
    if let Some(g) = as_factor(sig.gcd(n), n) {
        return Some(g);
    }
    let five = ExactInt::from_u64(5);
    let four = ExactInt::from_u64(4);
    let u = sub_mod(&mul_mod(&sig, &sig, n)?, &five, n)?;
    let v = mul_mod(&four, &sig, n)?;
    if let Some(g) = as_factor(u.gcd(n), n) {
        return Some(g);
    }
    if let Some(g) = as_factor(v.gcd(n), n) {
        return Some(g);
    }
    let u2 = mul_mod(&u, &u, n)?;
    let u3 = mul_mod(&u2, &u, n)?;
    let v2 = mul_mod(&v, &v, n)?;
    let v3 = mul_mod(&v2, &v, n)?;
    let vmu = sub_mod(&v, &u, n)?;
    let vmu2 = mul_mod(&vmu, &vmu, n)?;
    let vmu3 = mul_mod(&vmu2, &vmu, n)?;
    let three_u = mul_mod(&ExactInt::from_u64(3), &u, n)?;
    let three_u_v = add_mod(&three_u, &v, n)?;
    let num = mul_mod(&vmu3, &three_u_v, n)?;
    let den = mul_mod(&ExactInt::from_u64(16), &mul_mod(&u3, &v, n)?, n)?;
    if let Some(g) = as_factor(den.gcd(n), n) {
        return Some(g);
    }
    let inv = match mod_inv(&den, n) {
        Some(i) => i,
        None => return as_factor(den.gcd(n), n),
    };
    let a24 = mul_mod(&num, &inv, n)?;
    let p0 = MPt { x: u3, z: v3 };
    if let Some(g) = as_factor(p0.z.gcd(n), n) {
        return Some(g);
    }

    // Stage 1.
    let primes = primes_upto(b1.max(b2));
    let mut qpt = p0;
    let mut acc = ExactInt::one();
    let mut since = 0u32;
    for &p in &primes {
        if u64::from(p) > b1 {
            break;
        }
        let mut pk = u64::from(p);
        while let Some(next) = pk.checked_mul(u64::from(p)) {
            if next > b1 {
                break;
            }
            pk = next;
        }
        qpt = mont_mul(&qpt, &ExactInt::from_u64(pk), &a24, n)?;
        acc = mul_mod(&acc, &qpt.z, n)?;
        since += 1;
        if since == 32 {
            if let Some(g) = as_factor(acc.gcd(n), n) {
                return Some(g);
            }
            acc = ExactInt::one();
            since = 0;
        }
    }
    if let Some(g) = as_factor(acc.gcd(n), n).or_else(|| as_factor(qpt.z.gcd(n), n)) {
        return Some(g);
    }
    if b2 <= b1 {
        return None;
    }
    ecm_stage2(&qpt, b1, b2, &a24, n, &primes)
}

fn ecm_stage2(
    q: &MPt,
    b1: u64,
    b2: u64,
    a24: &ExactInt,
    n: &ExactInt,
    primes: &[u32],
) -> Option<ExactInt> {
    let d = 210u64;
    let two_q = mont_dbl(q, a24, n)?;
    // Odd multiples [1]Q, [3]Q, …, [d-1]Q.
    let mut odds: Vec<(u64, MPt)> = Vec::new();
    let mut cur = q.clone();
    odds.push((1, cur.clone()));
    // [3]Q = [1]Q + [2]Q, difference [1]Q (X is even).
    let mut nxt = mont_add(&cur, &two_q, q, n)?;
    let mut prev = cur;
    cur = nxt;
    let mut j = 3u64;
    while j < d {
        odds.push((j, cur.clone()));
        // [j+2] = [j] + [2], difference [j-2]
        nxt = mont_add(&cur, &two_q, &prev, n)?;
        prev = cur;
        cur = nxt;
        j += 2;
    }
    let d_q = mont_mul(q, &ExactInt::from_u64(d), a24, n)?;

    let mut is_prime_gap = vec![false; (b2 as usize).saturating_add(1)];
    for &p in primes {
        let pu = u64::from(p);
        if pu > b1 && pu <= b2 {
            is_prime_gap[pu as usize] = true;
        }
    }

    let mut k = (b1 / d) * d;
    if k <= b1 {
        k = k.saturating_add(d);
    }
    if k < d {
        k = d;
    }
    let mut g_pt = mont_mul(q, &ExactInt::from_u64(k), a24, n)?;
    let mut g_prev = if k >= d {
        mont_mul(q, &ExactInt::from_u64(k - d), a24, n)?
    } else {
        MPt {
            x: ExactInt::one(),
            z: ExactInt::zero(),
        }
    };
    let mut acc = ExactInt::one();
    while k <= b2.saturating_add(d) {
        for (j, sj) in &odds {
            let up = k.saturating_add(*j);
            let down = k.saturating_sub(*j);
            if (up <= b2
                && up > b1
                && (up as usize) < is_prime_gap.len()
                && is_prime_gap[up as usize])
                || (down > b1
                    && down <= b2
                    && (down as usize) < is_prime_gap.len()
                    && is_prime_gap[down as usize])
            {
                let cross = sub_mod(
                    &mul_mod(&g_pt.x, &sj.z, n)?,
                    &mul_mod(&sj.x, &g_pt.z, n)?,
                    n,
                )?;
                acc = mul_mod(&acc, &cross, n)?;
            }
        }
        if let Some(g) = as_factor(acc.gcd(n), n) {
            return Some(g);
        }
        let nxt = mont_add(&g_pt, &d_q, &g_prev, n)?;
        g_prev = g_pt;
        g_pt = nxt;
        k = k.saturating_add(d);
    }
    as_factor(acc.gcd(n), n)
}

pub(crate) fn primes_upto(limit: u64) -> Vec<u32> {
    let lim = limit.min(u32::MAX as u64) as usize;
    if lim < 2 {
        return Vec::new();
    }
    let mut mark = vec![true; lim + 1];
    mark[0] = false;
    mark[1] = false;
    let mut p = 2usize;
    while p * p <= lim {
        if mark[p] {
            let mut m = p * p;
            while m <= lim {
                mark[m] = false;
                m += p;
            }
        }
        p += 1;
    }
    let mut out = Vec::new();
    for (i, &ok) in mark.iter().enumerate() {
        if ok {
            out.push(i as u32);
        }
    }
    out
}
