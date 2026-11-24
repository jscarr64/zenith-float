//! FLINT 3.6-style primality proofs: a Brillhart-Lehmer-Selfridge `n-1` /
//! `n+1` front end with only a partial factorization, then the Cohen-Lenstra
//! Jacobi-sum APRCL test when that front end does not decide (Cohen, *A Course
//! in Computational Algebraic Number Theory*, algorithm 9.1.28).
//!
//! This is the proof path used by FLINT 3.6.0 `fmpz_is_prime` once trial
//! division, a square check, and a base-2 SPRP (plus BPSW as a composite
//! filter) have run. It does not ECM-factor `n-1`. A probable-prime result
//! is never returned as a proof.

use crate::modular::{is_even, miller_rabin, mod_pow, rem_nonneg, two};
use crate::primality::{bpsw, jacobi, lucas_uv, primes_upto};
use crate::ExactInt;
use alloc::vec;
use alloc::vec::Vec;
use core::cmp::Ordering;

fn one() -> ExactInt {
    ExactInt::one()
}

fn rem(a: &ExactInt, n: &ExactInt) -> ExactInt {
    rem_nonneg(a, n).unwrap_or_else(ExactInt::zero)
}

fn rem_u64(a: &ExactInt, m: u64) -> u64 {
    if m == 0 {
        return 0;
    }
    rem(a, &ExactInt::from_u64(m)).low_word()
}

fn bit_set(k: &ExactInt, i: usize) -> bool {
    k.shr_towards_zero(i).low_word() & 1 == 1
}

fn is_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    if n == 2 || n == 3 || n == 5 || n == 7 {
        return true;
    }
    if n.is_multiple_of(2) || n.is_multiple_of(3) || n.is_multiple_of(5) || n.is_multiple_of(7) {
        return false;
    }
    let v = ExactInt::from_u64(n);
    let wits: [u64; 13] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41];
    let ws: Vec<ExactInt> = wits
        .iter()
        .filter(|&&q| q < n)
        .map(|&q| ExactInt::from_u64(q))
        .collect();
    miller_rabin(&v, &ws)
}

fn factor_u64(mut m: u64) -> Vec<(u64, u32)> {
    let mut out = Vec::new();
    if m == 0 {
        return out;
    }
    let mut p = 2u64;
    while p * p <= m {
        if m.is_multiple_of(p) {
            let mut e = 0u32;
            while m.is_multiple_of(p) {
                m /= p;
                e += 1;
            }
            out.push((p, e));
        }
        p = if p == 2 { 3 } else { p + 2 };
    }
    if m > 1 {
        out.push((m, 1));
    }
    out
}

fn p_power_in_q(mut q: u64, p: u64) -> u32 {
    if p == 0 {
        return 0;
    }
    let mut k = 0u32;
    while q != 0 && q.is_multiple_of(p) {
        k += 1;
        q /= p;
    }
    k
}

fn pow_u64(base: u64, exp: u32) -> Option<u64> {
    let mut acc = 1u64;
    let mut b = base;
    let mut e = exp;
    while e > 0 {
        if e & 1 == 1 {
            acc = acc.checked_mul(b)?;
        }
        e >>= 1;
        if e > 0 {
            b = b.checked_mul(b)?;
        }
    }
    Some(acc)
}

fn primitive_root_prime(q: u64) -> u64 {
    let factors = factor_u64(q - 1);
    let mut g = 2u64;
    loop {
        let mut ok = true;
        for &(p, _) in &factors {
            let e = (q - 1) / p;
            let mut acc = 1u64;
            let mut b = g % q;
            let mut k = e;
            while k > 0 {
                if k & 1 == 1 {
                    acc = acc.wrapping_mul(b) % q;
                }
                b = b.wrapping_mul(b) % q;
                k >>= 1;
            }
            if acc == 1 {
                ok = false;
                break;
            }
        }
        if ok {
            return g;
        }
        g += 1;
        if g >= q {
            return 2;
        }
    }
}

fn f_table(q: u64) -> Vec<u64> {
    let g = primitive_root_prime(q);
    let mut g_table = vec![0u64; q as usize];
    let mut g_pow = g % q;
    for i in 1..q {
        g_table[g_pow as usize] = i;
        g_pow = g_pow.wrapping_mul(g) % q;
    }
    let mut table = vec![0u64; q as usize];
    g_pow = g % q;
    for i in 1..q {
        let one_minus = if g_pow == 0 { 1 } else { (q + 1 - g_pow) % q };
        table[i as usize] = g_table[one_minus as usize];
        g_pow = g_pow.wrapping_mul(g) % q;
    }
    table
}

fn isqrt(n: &ExactInt) -> ExactInt {
    if n.is_zero() || n.is_one() || n.is_negative() {
        // n is 0, 1 or negative.
        return if n.is_one() { ExactInt::one() } else { ExactInt::zero() };
    }
    let bits = n.bit_length();
    let mut hi = ExactInt::from_u64(2).pow((bits.div_ceil(2) + 1) as u64);
    let mut lo = ExactInt::zero();
    let one = one();
    while hi.sub(&lo).soft_cmp(&one) == Ordering::Greater {
        let mid = lo.add(&hi).shr_towards_zero(1);
        if mid.is_zero() {
            lo = ExactInt::one();
            continue;
        }
        let sq = mid.mul(&mid);
        match sq.soft_cmp(n) {
            Ordering::Greater => hi = mid,
            Ordering::Equal => return mid,
            Ordering::Less => lo = mid,
        }
    }
    lo
}

fn is_square(n: &ExactInt) -> bool {
    if n.is_negative() {
        return false;
    }
    let r = isqrt(n);
    r.mul(&r) == *n
}

fn nth_root(n: &ExactInt, k: u32) -> ExactInt {
    if k == 0 || n.is_zero() {
        return ExactInt::zero();
    }
    if k == 1 || n.is_one() {
        return n.clone();
    }
    if k == 2 {
        return isqrt(n);
    }
    let bits = n.bit_length();
    let mut hi = ExactInt::from_u64(2).pow((bits / k as usize + 2) as u64);
    if hi.is_zero() {
        hi = ExactInt::from_u64(2);
    }
    let mut lo = ExactInt::zero();
    let one = one();
    while hi.sub(&lo).soft_cmp(&one) == Ordering::Greater {
        let mid = lo.add(&hi).shr_towards_zero(1);
        if mid.is_zero() {
            lo = ExactInt::one();
            continue;
        }
        let pow = mid.pow(k as u64);
        match pow.soft_cmp(n) {
            Ordering::Greater => hi = mid,
            Ordering::Equal => return mid,
            Ordering::Less => lo = mid,
        }
    }
    lo
}

fn is_perfect_power(n: &ExactInt) -> bool {
    if n.soft_cmp(&two()) != Ordering::Greater {
        return false;
    }
    let bits = n.bit_length();
    for k in 2..=bits as u32 {
        let r = nth_root(n, k);
        if r.pow(k as u64) == *n {
            return true;
        }
    }
    false
}

/// Integer trial-limit used by FLINT: `floor((ln n)^3 / 100) + 20`.
fn dlog_limit(n: &ExactInt) -> u64 {
    let bits = n.bit_length() as u128;
    // ln(2) ~= 709/1024, so ln(n) ~= bits * 709/1024.
    let logd = bits.saturating_mul(709) / 1024;
    let cube = logd.saturating_mul(logd).saturating_mul(logd);
    u64::try_from(cube / 100 + 20).unwrap_or(u64::MAX / 4)
}

fn r_value(bits: usize) -> u64 {
    if bits <= 17 {
        6
    } else if bits <= 31 {
        12
    } else if bits <= 54 {
        36
    } else if bits <= 68 {
        72
    } else if bits <= 101 {
        180
    } else if bits <= 127 {
        360
    } else if bits <= 152 {
        720
    } else if bits <= 204 {
        1260
    } else if bits <= 268 {
        2520
    } else if bits <= 344 {
        5040
    } else if bits <= 525 {
        27720
    } else if bits <= 650 {
        55440
    } else if bits <= 774 {
        110880
    } else if bits <= 1566 {
        720720
    } else if bits <= 1999 {
        1_441_440
    } else if bits <= 2096 {
        1_663_200
    } else if bits <= 2165 {
        1_965_600
    } else if bits <= 2321 {
        2_162_160
    } else if bits <= 2377 {
        2_827_440
    } else if bits <= 2514 {
        3_326_400
    } else if bits <= 2588 {
        3_341_520
    } else if bits <= 2636 {
        3_603_600
    } else if bits <= 3028 {
        4_324_320
    } else if bits <= 3045 {
        5_654_880
    } else if bits <= 3080 {
        6_652_800
    } else if bits <= 3121 {
        6_683_040
    } else if bits <= 3283 {
        7_207_200
    } else if bits <= 3491 {
        8_648_640
    } else if bits <= 3726 {
        10_810_800
    } else if bits <= 3818 {
        12_972_960
    } else if bits <= 3977 {
        14_414_400
    } else if bits <= 4762 {
        21_621_600
    } else if bits <= 5068 {
        36_756_720
    } else if bits <= 5658 {
        43_243_200
    } else if bits <= 5960 {
        64_864_800
    } else if bits <= 6423 {
        73_513_440
    } else if bits <= 6900 {
        122_522_400
    } else if bits <= 9977 {
        367_567_200
    } else if bits <= 12713 {
        1_396_755_360
    } else {
        6_983_776_800
    }
}

struct Config {
    r: u64,
    s: ExactInt,
    /// Distinct prime factors of `R`.
    rs: Vec<u64>,
    /// Prime powers `q^e` that currently contribute to `s`.
    qs: Vec<(u64, u32)>,
    qs_used: Vec<bool>,
}

impl Config {
    fn new(n: &ExactInt) -> Self {
        let r = r_value(n.bit_length());
        let mut qs: Vec<(u64, u32)> = Vec::new();
        let mut s = ExactInt::one();
        let e2 = p_power_in_q(r, 2) + 2;
        if let Some(pow) = pow_u64(2, e2) {
            qs.push((2, e2));
            s = s.mul(&ExactInt::from_u64(pow));
        }
        let mut prime = 3u64;
        while 2u64.saturating_mul(prime.saturating_sub(1)) <= r {
            if is_prime_u64(prime) && r.is_multiple_of(prime - 1) {
                let e = p_power_in_q(r, prime) + 1;
                if let Some(pow) = pow_u64(prime, e) {
                    qs.push((prime, e));
                    s = s.mul(&ExactInt::from_u64(pow));
                }
            }
            prime = prime.saturating_add(if prime == 2 { 1 } else { 2 });
        }
        if is_prime_u64(r.saturating_add(1)) {
            qs.push((r + 1, 1));
            s = s.mul(&ExactInt::from_u64(r + 1));
        }
        let rs: Vec<u64> = factor_u64(r).into_iter().map(|(p, _)| p).collect();
        let mut conf = Self {
            r,
            s,
            rs,
            qs,
            qs_used: Vec::new(),
        };
        conf.qs_used = vec![true; conf.qs.len()];
        conf.reduce_s2(n);
        conf
    }

    fn p_ind(&self, p: u64) -> Option<usize> {
        self.rs.iter().position(|&x| x == p)
    }

    fn reduce_s2(&mut self, n: &ExactInt) {
        // Drop a prime power of s while s^2 stays > n, preferring the
        // largest work-per-bit (FLINT's w heuristic, in integers).
        let mut w_num = vec![0u128; self.qs.len()];
        let mut w_den = vec![1u128; self.qs.len()];
        for (i, &(q, exp)) in self.qs.iter().enumerate() {
            let mut num = 0u128;
            for (p, e) in factor_u64(q - 1) {
                // phi(p^e)^2 = (p^{e-1} (p-1))^2
                let phi = pow_u64(p, e.saturating_sub(1))
                    .unwrap_or(1)
                    .saturating_mul(p - 1);
                let phi = phi as u128;
                num = num.saturating_add(phi.saturating_mul(phi));
            }
            w_num[i] = num;
            let bits = 64u32.saturating_sub(q.leading_zeros());
            w_den[i] = (exp as u128).saturating_mul(bits.max(1) as u128).max(1);
        }
        loop {
            let mut best: Option<usize> = None;
            for i in 0..self.qs.len() {
                if !self.qs_used[i] {
                    continue;
                }
                let (q, exp) = self.qs[i];
                let Some(pow) = pow_u64(q, exp) else {
                    continue;
                };
                let new_s = match self.s.div_rem(&ExactInt::from_u64(pow)) {
                    Some((d, r)) if r.is_zero() => d,
                    _ => continue,
                };
                if new_s.mul(&new_s).soft_cmp(n) == Ordering::Greater {
                    let take = match best {
                        None => true,
                        Some(j) => {
                            // w_i >= w_j iff num_i * den_j >= num_j * den_i
                            w_num[i].saturating_mul(w_den[j]) >= w_num[j].saturating_mul(w_den[i])
                        }
                    };
                    if take {
                        best = Some(i);
                    }
                }
            }
            let Some(i) = best else {
                break;
            };
            let (q, exp) = self.qs[i];
            if let Some(pow) = pow_u64(q, exp) {
                if let Some((d, r)) = self.s.div_rem(&ExactInt::from_u64(pow)) {
                    if r.is_zero() {
                        self.s = d;
                    }
                }
            }
            self.qs_used[i] = false;
        }
    }
}

/// Element of `Z[zeta_{p^k}] / n`: coefficients of `1, zeta, ..., zeta^{p^k-1}`.
/// Little-endian unsigned limbs of `x` (a single `0` limb when `x` is zero).
fn limbs_of(x: &ExactInt) -> &[u64] {
    let l = x.as_limbs();
    if l.is_empty() {
        &[0]
    } else {
        l
    }
}

#[cfg(test)]
fn exact_from_limbs(limbs: &[u64]) -> ExactInt {
    ExactInt::from_le_words(crate::defs::Sign::Pos, limbs)
}

fn limb_len(a: &[u64]) -> usize {
    a.iter().rposition(|&w| w != 0).map_or(0, |i| i + 1)
}

fn limb_cmp(a: &[u64], b: &[u64]) -> Ordering {
    let la = limb_len(a);
    let lb = limb_len(b);
    if la != lb {
        return la.cmp(&lb);
    }
    for i in (0..la).rev() {
        if a[i] != b[i] {
            return a[i].cmp(&b[i]);
        }
    }
    Ordering::Equal
}

fn limbs_is_zero(a: &[u64]) -> bool {
    a.iter().all(|&w| w == 0)
}

fn shl_limbs(a: &mut [u64], shift: u32) -> u64 {
    if shift == 0 {
        return 0;
    }
    let spare = 64 - shift;
    let mut carry = 0u64;
    for w in a.iter_mut() {
        let new_carry = *w >> spare;
        *w = (*w << shift) | carry;
        carry = new_carry;
    }
    carry
}

fn shr_limbs(a: &mut [u64], shift: u32) {
    if shift == 0 {
        return;
    }
    let spare = 64 - shift;
    let mut carry = 0u64;
    for w in a.iter_mut().rev() {
        let new_carry = *w << spare;
        *w = (*w >> shift) | carry;
        carry = new_carry;
    }
}

fn add_raw(dest: &mut [u64], src: &[u64]) {
    let mut carry = 0u128;
    let mut i = 0;
    while i < src.len() {
        let t = dest[i] as u128 + src[i] as u128 + carry;
        dest[i] = t as u64;
        carry = t >> 64;
        i += 1;
    }
    // No overflow past `dest` (Jeff 2026-11-24 4:40 PM CT item 5, was a debug_assert): the only caller,
    // `UnityZp::mul`, adds fewer than 2^64 products a*b with a, b < n < 2^(64*nl) into a slot of 2*nl + 2 limbs, so
    // the sum is below 2^(64*(2*nl + 1)) and the carry stops inside `dest`.
    while carry > 0 {
        let t = dest[i] as u128 + carry;
        dest[i] = t as u64;
        carry = t >> 64;
        i += 1;
    }
}

fn add_mod_limbs(dest: &mut [u64], src: &[u64], n: &[u64]) {
    let mut carry = 0u128;
    for i in 0..dest.len() {
        let s = if i < src.len() { src[i] as u128 } else { 0 };
        let t = dest[i] as u128 + s + carry;
        dest[i] = t as u64;
        carry = t >> 64;
    }
    if carry != 0 || limb_cmp(dest, n) != Ordering::Less {
        let _ = sub_raw(dest, n);
    }
}

fn sub_raw(dest: &mut [u64], src: &[u64]) -> bool {
    let mut borrow = 0u64;
    let mut i = 0;
    while i < src.len() {
        let (d1, b1) = dest[i].overflowing_sub(src[i]);
        let (d2, b2) = d1.overflowing_sub(borrow);
        dest[i] = d2;
        borrow = u64::from(b1 || b2);
        i += 1;
    }
    while borrow != 0 && i < dest.len() {
        let (d2, b2) = dest[i].overflowing_sub(borrow);
        dest[i] = d2;
        borrow = u64::from(b2);
        i += 1;
    }
    borrow != 0
}

fn sub_mod_limbs(dest: &mut [u64], src: &[u64], n: &[u64]) {
    if sub_raw(dest, src) {
        let mut carry = 0u128;
        for i in 0..dest.len() {
            let s = if i < n.len() { n[i] as u128 } else { 0 };
            let t = dest[i] as u128 + s + carry;
            dest[i] = t as u64;
            carry = t >> 64;
        }
    }
}

fn add_one_mod(dest: &mut [u64], n: &[u64]) {
    add_mod_limbs(dest, &[1], n);
}

fn sub_one_mod(dest: &mut [u64], n: &[u64]) {
    sub_mod_limbs(dest, &[1], n);
}

fn mul_raw(a: &[u64], b: &[u64], out: &mut [u64]) {
    out.fill(0);
    for i in 0..a.len() {
        if a[i] == 0 {
            continue;
        }
        let mut carry = 0u128;
        for j in 0..b.len() {
            let t = out[i + j] as u128 + a[i] as u128 * b[j] as u128 + carry;
            out[i + j] = t as u64;
            carry = t >> 64;
        }
        let mut k = i + b.len();
        while carry > 0 {
            let t = out[k] as u128 + carry;
            out[k] = t as u64;
            carry = t >> 64;
            k += 1;
        }
    }
}

/// `dest -= q * v`, with `dest` one limb longer than `v`. Returns true if a
/// leftover borrow remains in the high limb (the estimate `q` was too large).
fn submul_q(dest: &mut [u64], v: &[u64], q: u64) -> bool {
    let mut hold = 0u128;
    for i in 0..v.len() {
        hold += v[i] as u128 * q as u128;
        let sub = hold as u64;
        hold >>= 64;
        let (d, br) = dest[i].overflowing_sub(sub);
        dest[i] = d;
        if br {
            hold += 1;
        }
    }
    let i = v.len();
    let (d, br) = dest[i].overflowing_sub(hold as u64);
    dest[i] = d;
    br || (hold >> 64) != 0
}

fn add_assign_n(dest: &mut [u64], v: &[u64]) {
    let mut carry = 0u128;
    for i in 0..v.len() {
        let t = dest[i] as u128 + v[i] as u128 + carry;
        dest[i] = t as u64;
        carry = t >> 64;
    }
    let i = v.len();
    if carry > 0 && i < dest.len() {
        let t = dest[i] as u128 + carry;
        dest[i] = t as u64;
    }
}

fn rem_by_u64(u: &[u64], v: u64) -> u64 {
    if v == 0 {
        return 0;
    }
    let mut r = 0u128;
    let len = limb_len(u);
    for i in (0..len).rev() {
        r = ((r << 64) | u[i] as u128) % v as u128;
    }
    r as u64
}

/// Remainder `u mod v` for unsigned little-endian limbs (`v` nonzero). The
/// result is padded or truncated to `v.len()`.
fn rem_limbs(u: &[u64], v: &[u64]) -> Vec<u64> {
    // `v` is never zero: every caller passes a `UnityZp::n`, which `UnityZp::new` refuses when zero.
    let n = limb_len(v);
    let mut r = vec![0u64; v.len()];
    let ulen = limb_len(u);
    if ulen == 0 {
        return r;
    }
    if n == 1 {
        r[0] = rem_by_u64(u, v[0]);
        return r;
    }
    if limb_cmp(&u[..ulen], &v[..n]) == Ordering::Less {
        r[..ulen].copy_from_slice(&u[..ulen]);
        return r;
    }

    let shift = v[n - 1].leading_zeros();
    let mut vn = vec![0u64; n];
    vn.copy_from_slice(&v[..n]);
    let _ = shl_limbs(&mut vn, shift);

    // `un` holds the normalized numerator plus one extra high limb.
    let mut un = vec![0u64; ulen + 2];
    un[..ulen].copy_from_slice(&u[..ulen]);
    let carry = shl_limbs(&mut un[..ulen], shift);
    un[ulen] = carry;

    let un_len = if un[ulen] != 0 { ulen + 1 } else { ulen };
    let m = un_len - n;
    for j in (0..=m).rev() {
        let uj_n = un[j + n] as u128;
        let uj_n1 = un[j + n - 1] as u128;
        let vnm1 = vn[n - 1] as u128;
        let dividend = (uj_n << 64) | uj_n1;
        let mut qhat = if uj_n == vnm1 { u64::MAX as u128 } else { dividend / vnm1 };
        let mut rhat = dividend - qhat * vnm1;
        let vn_2 = vn[n - 2] as u128;
        let uj_n2 = un[j + n - 2] as u128;
        // `rhat` can start at or above 2^64 when `qhat` was capped at `u64::MAX`.
        // Shifting that in `u128` would wrap and pull `qhat` down too far.
        loop {
            if rhat >= (1u128 << 64) {
                break;
            }
            if qhat * vn_2 <= (rhat << 64) | uj_n2 {
                break;
            }
            qhat -= 1;
            rhat += vnm1;
        }
        let too_much = submul_q(&mut un[j..j + n + 1], &vn, qhat as u64);
        if too_much {
            add_assign_n(&mut un[j..j + n + 1], &vn);
        }
    }

    shr_limbs(&mut un[..n], shift);
    r[..n].copy_from_slice(&un[..n]);
    r
}

/// The modulus `n` is borrowed from the caller's `ExactInt`: every element of one test shares it.
struct UnityZp<'a> {
    p: u64,
    exp: u32,
    n: &'a [u64],
    /// Packed little-endian coefficients: slot `i` is `c[i * nlen .. (i + 1) * nlen]`.
    c: Vec<u64>,
    slots: usize,
}

impl<'a> UnityZp<'a> {
    /// `None` when `n` is zero: every reduction below is mod `n`, so `n != 0` is checked once here and
    /// [`rem_limbs`] can rely on it (Jeff 2026-11-24 4:40 PM CT item 5, was a debug_assert there).
    fn new(p: u64, exp: u32, n: &'a ExactInt) -> Option<Self> {
        let n_limbs = limbs_of(n);
        if limb_len(n_limbs) == 0 {
            return None;
        }
        let slots = pow_u64(p, exp).unwrap_or(1) as usize;
        let nlen = n_limbs.len();
        Some(Self {
            p,
            exp,
            n: n_limbs,
            c: vec![0; slots.max(1) * nlen],
            slots: slots.max(1),
        })
    }

    fn zeroed(&self) -> Self {
        Self {
            p: self.p,
            exp: self.exp,
            n: self.n,
            c: vec![0; self.slots.max(1) * self.n.len()],
            slots: self.slots.max(1),
        }
    }

    fn nlen(&self) -> usize {
        self.n.len()
    }

    fn p_pow(&self) -> usize {
        pow_u64(self.p, self.exp).unwrap_or(1) as usize
    }

    fn p_pow_prev(&self) -> usize {
        if self.exp == 0 {
            1
        } else {
            pow_u64(self.p, self.exp - 1).unwrap_or(1) as usize
        }
    }

    fn set_zero(&mut self) {
        self.c.fill(0);
    }

    fn slot(&self, i: usize) -> Option<&[u64]> {
        let l = self.nlen();
        if i >= self.slots {
            return None;
        }
        Some(&self.c[i * l..(i + 1) * l])
    }

    fn slot_mut(&mut self, i: usize) -> &mut [u64] {
        let l = self.nlen();
        &mut self.c[i * l..(i + 1) * l]
    }

    fn ensure(&mut self, i: usize) {
        if i < self.slots {
            return;
        }
        let l = self.nlen();
        self.c.resize((i + 1) * l, 0);
        self.slots = i + 1;
    }

    fn set_coeff_limbs(&mut self, i: usize, x: &[u64]) {
        self.ensure(i);
        let l = self.nlen();
        let r = rem_limbs(x, self.n);
        let dest = self.slot_mut(i);
        dest.fill(0);
        let take = r.len().min(l);
        dest[..take].copy_from_slice(&r[..take]);
    }

    fn set_ui(&mut self, i: usize, x: u64) {
        self.ensure(i);
        let n = self.n;
        let dest = self.slot_mut(i);
        dest.fill(0);
        if x != 0 {
            let r = rem_limbs(&[x], n);
            let take = r.len().min(dest.len());
            dest[..take].copy_from_slice(&r[..take]);
        }
    }

    fn inc(&mut self, i: usize) {
        self.ensure(i);
        let n = self.n;
        add_one_mod(self.slot_mut(i), n);
    }

    fn dec(&mut self, i: usize) {
        self.ensure(i);
        let n = self.n;
        sub_one_mod(self.slot_mut(i), n);
    }

    fn mul_scalar_ui(&self, s: u64) -> Self {
        let mut out = self.zeroed();
        if s == 0 {
            return out;
        }
        let sv = [s];
        let mut prod = vec![0u64; self.nlen() + 1];
        for i in 0..self.slots {
            if let Some(a) = self.slot(i) {
                if limbs_is_zero(a) {
                    continue;
                }
                mul_raw(a, &sv, &mut prod);
                out.set_coeff_limbs(i, &prod);
            }
        }
        out
    }

    /// Adds (or subtracts) slot `from` mod `n` onto each slot of `to`, all below `from`, then zeroes
    /// it. The source limbs are read in place.
    fn fold_slot(&mut self, from: usize, to: impl Iterator<Item = usize>, add: bool) {
        let l = self.n.len();
        let (lo, hi) = self.c.split_at_mut(from * l);
        let src = &mut hi[..l];
        for t in to {
            let dest = &mut lo[t * l..(t + 1) * l];
            if add {
                add_mod_limbs(dest, src, self.n);
            } else {
                sub_mod_limbs(dest, src, self.n);
            }
        }
        src.fill(0);
    }

    fn reduce_xpk(&mut self) {
        let pk = self.p_pow();
        let mut i = self.slots;
        while i > pk {
            i -= 1;
            let Some(slot) = self.slot(i) else {
                continue;
            };
            if limbs_is_zero(slot) {
                continue;
            }
            self.fold_slot(i, core::iter::once(i - pk), true);
        }
    }

    fn reduce_cyclotomic(&mut self) {
        if self.slots == 0 {
            return;
        }
        let ppow = self.p_pow_prev();
        let cycl = (self.p as usize - 1) * ppow;
        let p = self.p as usize;
        let mut i = self.slots;
        while i > cycl {
            i -= 1;
            if i >= self.slots {
                continue;
            }
            let Some(slot) = self.slot(i) else {
                continue;
            };
            if limbs_is_zero(slot) {
                continue;
            }
            self.fold_slot(i, (0..p - 1).map(|j| i - cycl + j * ppow), false);
        }
        self.trim();
    }

    fn trim(&mut self) {
        let l = self.nlen();
        while self.slots > 0 && limbs_is_zero(&self.c[(self.slots - 1) * l..self.slots * l]) {
            self.slots -= 1;
        }
        self.c.truncate(self.slots * l);
    }

    fn reduce_divmod(&mut self) {
        self.reduce_xpk();
        self.reduce_cyclotomic();
        self.trim();
    }

    fn mul(&self, other: &Self) -> Self {
        let mut out = Self {
            p: self.p,
            exp: self.exp,
            n: self.n,
            c: Vec::new(),
            slots: 0,
        };
        if self.slots == 0 || other.slots == 0 {
            return out;
        }
        let nl = self.nlen();
        let wide = 2 * nl + 2;
        let out_slots = self.slots + other.slots;
        let mut acc = vec![0u64; out_slots * wide];
        let mut prod = vec![0u64; 2 * nl + 1];
        for i in 0..self.slots {
            let a = match self.slot(i) {
                Some(a) if !limbs_is_zero(a) => a,
                _ => continue,
            };
            for j in 0..other.slots {
                let b = match other.slot(j) {
                    Some(b) if !limbs_is_zero(b) => b,
                    _ => continue,
                };
                mul_raw(a, b, &mut prod);
                let off = (i + j) * wide;
                add_raw(&mut acc[off..off + wide], &prod);
            }
        }
        out.c = vec![0u64; out_slots * nl];
        out.slots = out_slots;
        for k in 0..out_slots {
            let r = rem_limbs(&acc[k * wide..(k + 1) * wide], self.n);
            let dest = &mut out.c[k * nl..(k + 1) * nl];
            let take = r.len().min(nl);
            dest[..take].copy_from_slice(&r[..take]);
        }
        out.reduce_divmod();
        out
    }

    fn sqr(&self) -> Self {
        self.mul(self)
    }

    fn pow_int(&self, e: &ExactInt) -> Self {
        let mut out = self.zeroed();
        out.set_ui(0, 1);
        if e.is_zero() {
            return out;
        }
        if e.is_negative() {
            return out;
        }
        let bits = e.bit_length();
        let mut i = bits;
        while i > 0 {
            i -= 1;
            out = out.sqr();
            if bit_set(e, i) {
                out = out.mul(self);
            }
        }
        out
    }

    fn pow_ui(&self, e: u64) -> Self {
        self.pow_int(&ExactInt::from_u64(e))
    }

    fn aut_inv(&self, x: u64) -> Self {
        let p = self.p as usize;
        let p_pow1 = self.p_pow_prev();
        let p_pow2 = self.p_pow();
        let m = (p - 1) * p_pow1;
        let mut out = self.zeroed();
        let n = self.n;
        for i in 0..m {
            let g_ind = ((x as u128 * i as u128) % p_pow2 as u128) as usize;
            out.ensure(i);
            if let Some(src) = self.slot(g_ind) {
                out.slot_mut(i).copy_from_slice(src);
            }
        }
        for i in m..p_pow2 {
            let g_ind = ((x as u128 * i as u128) % p_pow2 as u128) as usize;
            // A missing slot is zero, and subtracting zero leaves the slot as it is.
            let g_coeff = self.slot(g_ind);
            for j in 1..p {
                let f_ind = i - j * p_pow1;
                out.ensure(f_ind);
                if let Some(g) = g_coeff {
                    sub_mod_limbs(out.slot_mut(f_ind), g, n);
                }
            }
        }
        out
    }

    /// `self == other` in `Z[zeta]/n`: the difference, written into a new element, reduces to zero
    /// (the reduction is linear and every coefficient is a residue in `[0, n)`).
    fn equal_to(&self, other: &Self) -> bool {
        let slots = self.slots.max(other.slots);
        let mut d = Self {
            p: self.p,
            exp: self.exp,
            n: self.n,
            c: vec![0; slots.max(1) * self.nlen()],
            slots: slots.max(1),
        };
        for i in 0..slots {
            let dest = d.slot_mut(i);
            if let Some(a) = self.slot(i) {
                dest.copy_from_slice(a);
            }
            if let Some(b) = other.slot(i) {
                sub_mod_limbs(dest, b, self.n);
            }
        }
        d.reduce_cyclotomic();
        (0..d.slots).all(|i| d.slot(i).is_none_or(limbs_is_zero))
    }

    /// `h` such that `self = zeta_{p^k}^h`, or `-1`.
    fn is_unity(&self) -> i64 {
        let p_pow = self.p_pow();
        for i in 0..p_pow {
            let mut u = self.zeroed();
            u.set_ui(i, 1);
            if u.equal_to(self) {
                return i as i64;
            }
        }
        -1
    }
}

fn jacobi_sum_general(f: &mut UnityZp<'_>, table: &[u64], p: u64, q: u64, k: u32, a: u64, b: u64) {
    f.set_zero();
    let pow_dec = pow_u64(p, k.saturating_sub(1)).unwrap_or(1);
    let size = (p - 1) * pow_dec;
    let pow = pow_dec.saturating_mul(p);
    if pow == 0 {
        return;
    }
    let q_i = q as usize;
    for i in 1..q.saturating_sub(1) {
        let t = if (i as usize) < table.len() { table[i as usize] } else { 0 };
        let l = (a as u128 * i as u128 + b as u128 * t as u128) % pow as u128;
        if l < size as u128 {
            f.inc(l as usize);
        } else {
            let mut l = l;
            for _ in 0..p - 1 {
                l = l.saturating_sub(pow_dec as u128);
                f.dec(l as usize);
            }
        }
        let _ = q_i;
    }
}

fn jacobi_sum_pq(n: &ExactInt, q: u64, p: u64) -> Option<UnityZp<'_>> {
    let table = f_table(q);
    let k = p_power_in_q(q - 1, p);
    let mut f = UnityZp::new(p, k, n)?;
    jacobi_sum_general(&mut f, &table, p, q, k, 1, 1);
    Some(f)
}

fn jacobi_sum_2q_one(n: &ExactInt, q: u64) -> Option<UnityZp<'_>> {
    let table = f_table(q);
    let k = p_power_in_q(q - 1, 2);
    let mut f = UnityZp::new(2, k, n)?;
    jacobi_sum_general(&mut f, &table, 2, q, k, 2, 1);
    Some(f)
}

fn jacobi_sum_2q_two(n: &ExactInt, q: u64) -> Option<UnityZp<'_>> {
    let table = f_table(q);
    let k = p_power_in_q(q - 1, 2);
    let b = pow_u64(2, k.saturating_sub(3)).unwrap_or(1);
    let a = 3 * b;
    let mut f = UnityZp::new(2, k, n)?;
    jacobi_sum_general(&mut f, &table, 2, q, k, a, b);
    Some(f)
}

fn check_pk(j: &UnityZp<'_>, u: &ExactInt, v: u64) -> i64 {
    let r = j.p_pow() as u64;
    let mut j0 = j.zeroed();
    j0.set_ui(0, 1);
    let mut jv = j.zeroed();
    jv.set_ui(0, 1);
    for i in 1..=r {
        if i % j.p == 0 {
            continue;
        }
        let mut temp = j.pow_ui(i);
        temp.reduce_cyclotomic();
        let aut = temp.aut_inv(i);
        j0 = j0.mul(&aut);

        let e = (v as u128 * i as u128) / r as u128;
        let mut temp = j.pow_ui(e as u64);
        temp.reduce_cyclotomic();
        let aut = temp.aut_inv(i);
        jv = jv.mul(&aut);
    }
    let temp = j0.pow_int(u);
    let j0 = jv.mul(&temp);
    j0.is_unity()
}

fn check_21(q: u64, n: &ExactInt) -> i64 {
    let ndec = n.sub(&one());
    let qpow = n.sub(&ExactInt::from_u64(q));
    let exp = ndec.shr_towards_zero(1);
    let Some(qpow) = mod_pow(&qpow, &exp, n) else {
        return -1;
    };
    if qpow.is_one() {
        0
    } else if qpow == ndec {
        1
    } else {
        -1
    }
}

fn check_22(j: &UnityZp<'_>, u: &ExactInt, v: u64, q: u64) -> i64 {
    let j_pow = j.mul(j);
    let j0 = j_pow.mul_scalar_ui(q);
    let jv = if v == 1 {
        let mut t = j.zeroed();
        t.set_ui(0, 1);
        t
    } else if v == 3 {
        j_pow
    } else {
        j.zeroed()
    };
    let tmp = j0.pow_int(u);
    let j0 = jv.mul(&tmp);
    j0.is_unity()
}

fn check_2k(j: &UnityZp<'_>, j2_1: &UnityZp<'_>, j2_2: &UnityZp<'_>, u: &ExactInt, v: u64) -> i64 {
    let r = j.p_pow() as u64;
    let mut j0 = j.zeroed();
    j0.set_ui(0, 1);
    let mut jv = j.zeroed();
    jv.set_ui(0, 1);
    let j_j0 = j.mul(j2_1);
    let mut i = 1u64;
    while i < r {
        // i == 1 (mod 8), then i == 3 (mod 8)
        for _pass in 0..2 {
            let mut temp = j_j0.pow_ui(i);
            temp.reduce_cyclotomic();
            let aut = temp.aut_inv(i);
            j0 = j0.mul(&aut);

            let e = (v as u128 * i as u128) / r as u128;
            let mut temp = j_j0.pow_ui(e as u64);
            temp.reduce_cyclotomic();
            let aut = temp.aut_inv(i);
            jv = jv.mul(&aut);

            if _pass == 0 {
                i += 2;
            }
        }
        i += 6;
    }
    if v % 8 != 1 && v % 8 != 3 {
        let tmp = j2_2.mul(j2_2);
        jv = jv.mul(&tmp);
    }
    let temp = j0.pow_int(u);
    let j0 = jv.mul(&temp);
    j0.is_unity()
}

/// Additional L_p test. Returns 2 if composite, 1 if this p is now satisfied,
/// 0 if no suitable q was found (cannot prove).
fn additional_test(n: &ExactInt, p: u64) -> u8 {
    let mut result = 0u8;
    let mut p_counter = 50i32;
    let mut m = 9u64;
    let mut found_q = 0u64;
    while p_counter > 0 {
        let q = 2u64.saturating_mul(m).saturating_mul(p).saturating_add(1);
        if is_prime_u64(q) && rem_u64(n, q) != 0 {
            let qmod = ExactInt::from_u64(q);
            let e = (q - 1) / p;
            if let Some(npow) = mod_pow(n, &ExactInt::from_u64(e), &qmod) {
                if !npow.is_one() {
                    found_q = q;
                    break;
                }
            }
            p_counter -= 1;
        }
        m = m.saturating_add(2);
        if m > (1u64 << 20) {
            p_counter = 0;
            break;
        }
    }
    if p_counter != 0 && found_q != 0 {
        let q = found_q;
        if rem_u64(n, q) == 0 && *n != ExactInt::from_u64(q) {
            result = 2;
        } else {
            let k = p_power_in_q(q - 1, p);
            let Some(jacobi_sum) = jacobi_sum_pq(n, q, p) else {
                return 0;
            };
            let r = pow_u64(p, k).unwrap_or(1);
            let u = n
                .div_rem(&ExactInt::from_u64(r))
                .map(|(q, _)| q)
                .unwrap_or_else(ExactInt::zero);
            let v = rem_u64(n, r);
            if p == 2 {
                let h = check_22(&jacobi_sum, &u, v, q);
                if h < 0 || h % 2 == 0 {
                    result = 2;
                } else {
                    let ndec = n.sub(&one());
                    let exp = ndec.shr_towards_zero(1);
                    match mod_pow(&ExactInt::from_u64(q), &exp, n) {
                        Some(qpow) if qpow == ndec => result = 1,
                        _ => result = 2,
                    }
                }
            } else {
                let h = check_pk(&jacobi_sum, &u, v);
                if h < 0 || h % p as i64 == 0 {
                    result = 2;
                } else {
                    result = 1;
                }
            }
        }
    }
    if p_counter == 0 && (rem_u64(n, p) == 0 || is_perfect_power(n)) {
        result = 2;
    }
    result
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Status {
    Unknown,
    Prime,
    Composite,
    Probabprime,
}

fn is_mul_coprime(r: u64, s: &ExactInt, n: &ExactInt) -> bool {
    let g_r = if r == 0 { 1 } else { gcd_u64(r, rem_u64(n, r)) };
    if g_r != 1 {
        return false;
    }
    rem(n, s).gcd(s).is_one()
}

fn gcd_u64(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

fn is_prime_jacobi(n: &ExactInt) -> Status {
    if n == &two() || *n == ExactInt::from_u64(3) {
        return Status::Prime;
    }
    let conf = Config::new(n);
    let ndec = n.sub(&one());
    let ndecdiv = ndec.shr_towards_zero(1);
    let mut result = Status::Probabprime;
    let mut lambdas = vec![0i8; conf.rs.len()];
    let nmod4 = rem_u64(n, 4);
    for (i, &p) in conf.rs.iter().enumerate() {
        if p > 2 {
            let p2 = ExactInt::from_u64(p.saturating_mul(p));
            if let Some(t) = mod_pow(n, &ExactInt::from_u64(p - 1), &p2) {
                lambdas[i] = if t.is_one() { 0 } else { 1 };
            }
        }
    }
    if !is_mul_coprime(conf.r, &conf.s, n) {
        result = Status::Composite;
    }
    'qs: for (qi, &(q, _)) in conf.qs.iter().enumerate() {
        if !conf.qs_used[qi] {
            continue;
        }
        if result == Status::Composite {
            break;
        }
        if *n == ExactInt::from_u64(q) {
            result = Status::Prime;
            break;
        }
        let q_factors = factor_u64(q - 1);
        for &(p, k) in &q_factors {
            if result == Status::Composite {
                break 'qs;
            }
            let r = pow_u64(p, k).unwrap_or(1);
            let Some(pind) = conf.p_ind(p) else {
                continue;
            };
            let mut q_pow = ExactInt::from_u64(q);
            if lambdas[pind] == 0 && p == 2 {
                q_pow = mod_pow(&q_pow, &ndecdiv, n).unwrap_or(q_pow);
            }
            let u = n
                .div_rem(&ExactInt::from_u64(r))
                .map(|(d, _)| d)
                .unwrap_or_else(ExactInt::zero);
            let v = rem_u64(n, r);
            let Some(jacobi_sum) = jacobi_sum_pq(n, q, p) else {
                return Status::Unknown;
            };
            if p == 2 && k == 1 {
                let h = check_21(q, n);
                if h < 0 {
                    result = Status::Composite;
                }
                if lambdas[pind] == 0 && h == 1 && nmod4 == 1 {
                    lambdas[pind] = 1;
                }
            }
            if p == 2 && k == 2 {
                let h = check_22(&jacobi_sum, &u, v, q);
                if h < 0 {
                    result = Status::Composite;
                }
                if h % 2 != 0 && lambdas[pind] == 0 && q_pow == ndec {
                    lambdas[pind] = 1;
                }
            }
            if p == 2 && k >= 3 {
                let (Some(j31), Some(j2)) = (jacobi_sum_2q_one(n, q), jacobi_sum_2q_two(n, q))
                else {
                    return Status::Unknown;
                };
                let h = check_2k(&jacobi_sum, &j31, &j2, &u, v);
                if h < 0 {
                    result = Status::Composite;
                }
                if h % 2 != 0 && lambdas[pind] == 0 && q_pow == ndec {
                    lambdas[pind] = 1;
                }
            }
            if p != 2 {
                let h = check_pk(&jacobi_sum, &u, v);
                if h < 0 {
                    result = Status::Composite;
                }
                if h % p as i64 != 0 && lambdas[pind] == 0 {
                    lambdas[pind] = 1;
                }
            }
        }
    }
    if result == Status::Probabprime {
        for (lam, &r) in lambdas.iter_mut().zip(conf.rs.iter()) {
            if *lam == 0 {
                match additional_test(n, r) {
                    2 => result = Status::Composite,
                    1 => {}
                    _ => result = Status::Unknown,
                }
                if result == Status::Probabprime {
                    *lam = 1;
                }
            }
        }
    }
    if result == Status::Probabprime {
        if final_division(n, &conf.s, conf.r) {
            result = Status::Prime;
        } else {
            result = Status::Composite;
        }
    }
    result
}

fn final_division(n: &ExactInt, s: &ExactInt, r: u64) -> bool {
    let nmul = rem(n, s);
    // nmul^i mod s; `None` while i = 1 (nmul read in place).
    let mut npow_owned: Option<ExactInt> = None;
    for _ in 1..=r {
        let npow = npow_owned.as_ref().unwrap_or(&nmul);
        if npow.is_one() {
            break;
        }
        if let Some((_, remn)) = n.div_rem(npow) {
            if remn.is_zero() && n != npow && !npow.is_one() {
                return false;
            }
        }
        let next = rem(&npow.mul(&nmul), s);
        npow_owned = Some(next);
    }
    true
}

/// Cohen-Lenstra Jacobi-sum APRCL proof for an odd strong probable prime `n`.
pub(crate) fn aprcl_is_prime(n: &ExactInt) -> bool {
    matches!(is_prime_jacobi(n), Status::Prime)
}

fn trial_pm1_pp1(n: &ExactInt, limit: u64) -> (Vec<u64>, Vec<u64>) {
    let mut pm1 = Vec::new();
    let mut pp1 = Vec::new();
    for p in primes_upto(limit) {
        let p = p as u64;
        if p >= limit {
            break;
        }
        let r = rem_u64(n, p);
        if r == 1 {
            pm1.push(p);
        }
        if r + 1 == p {
            pp1.push(p);
        }
    }
    (pm1, pp1)
}

fn remove_primes(m: &mut ExactInt, primes: &[u64], fac: &mut Vec<(ExactInt, u32)>) {
    for &p in primes {
        let pv = ExactInt::from_u64(p);
        let mut e = 0u32;
        loop {
            match m.div_rem(&pv) {
                Some((d, r)) if r.is_zero() => {
                    *m = d;
                    e += 1;
                }
                _ => break,
            }
        }
        if e > 0 {
            fac.push((pv, e));
        }
    }
}

fn product_fac(fac: &[(ExactInt, u32)]) -> ExactInt {
    let mut f = ExactInt::one();
    for (p, e) in fac {
        f = f.mul(&p.pow(*e as u64));
    }
    f
}

/// Try to prove a BPSW cofactor prime (recursive FLINT step).
fn prove_if_prime(r: &ExactInt) -> bool {
    if r.is_one() || r.is_zero() {
        return false;
    }
    if r == &two() {
        return true;
    }
    if is_even(r) {
        return false;
    }
    if !bpsw(r) {
        return false;
    }
    // Recurse through the same proof path; `r` is strictly smaller than the
    // original `n` because it divides `n-1` or `n+1`.
    prove_prime(r)
}

fn pocklington(n: &ExactInt, pm1: &[u64]) -> Option<(ExactInt, ExactInt, bool)> {
    let mut r = n.sub(&one());
    let mut fac: Vec<(ExactInt, u32)> = Vec::new();
    remove_primes(&mut r, pm1, &mut fac);
    if bpsw(&r) && prove_if_prime(&r) {
        fac.push((r, 1));
        r = ExactInt::one();
    }
    let f = product_fac(&fac);
    // Witness search.
    let mut a = two();
    loop {
        let Some(pow) = mod_pow(&a, &r, n) else {
            return Some((f, r, false));
        };
        let Some(pow2) = mod_pow(&pow, &f, n) else {
            return Some((f, r, false));
        };
        if !pow2.is_one() {
            return Some((f, r, false));
        }
        let mut c = ExactInt::one();
        let mut ok = true;
        for (q, _) in &fac {
            let ex = match f.div_rem(q) {
                Some((e, _)) => e,
                None => return Some((f, r, false)),
            };
            let Some(mut t) = mod_pow(&pow, &ex, n) else {
                return Some((f, r, false));
            };
            t = t.sub(&one());
            if t.is_negative() {
                t = t.add(n);
            }
            if t.is_zero() {
                ok = false;
                break;
            }
            c = rem(&c.mul(&t), n);
        }
        if ok {
            let g = c.gcd(n);
            return Some((f, r, g.is_one()));
        }
        a = a.add(&one());
        if a.soft_cmp(n) != Ordering::Less {
            return Some((f, r, false));
        }
    }
}

fn morrison(n: &ExactInt, pp1: &[u64]) -> Option<(ExactInt, ExactInt, bool)> {
    let mut r = n.add(&one());
    let mut fac: Vec<(ExactInt, u32)> = Vec::new();
    remove_primes(&mut r, pp1, &mut fac);
    if bpsw(&r) && prove_if_prime(&r) {
        fac.push((r, 1));
        r = ExactInt::one();
    }
    let f = product_fac(&fac);
    let mut b = 1u64;
    loop {
        let bv = ExactInt::from_u64(b);
        let g = bv.gcd(n);
        if g == *n {
            b = b.saturating_add(1);
            if b == 0 {
                return Some((f, r, false));
            }
            continue;
        }
        if !g.is_one() {
            return Some((f, r, false));
        }
        let mut a = 3u64;
        let mut d;
        loop {
            let av = ExactInt::from_u64(a);
            d = av.mul(&av).sub(&ExactInt::from_u64(4 * b));
            if jacobi(&d, n) == -1 {
                break;
            }
            a = a.saturating_add(1);
            if a > 10_000 {
                return Some((f, r, false));
            }
        }
        let av = ExactInt::from_u64(a);
        // U((n+1)/F) via the existing Lucas chain (k = R).
        let Some((um, _, _)) = lucas_uv(n, &av, &bv, &d, &r) else {
            return Some((f, r, false));
        };
        // U(n+1) = U_{R F}.
        let np1 = n.add(&one());
        let Some((ukm, _, _)) = lucas_uv(n, &av, &bv, &d, &np1) else {
            return Some((f, r, false));
        };
        if !ukm.is_zero() {
            return Some((f, r, false));
        }
        let mut c = ExactInt::one();
        let mut ok = true;
        for (q, _) in &fac {
            let ex = match np1.div_rem(q) {
                Some((e, _)) => e,
                None => return Some((f, r, false)),
            };
            let Some((uk, _, _)) = lucas_uv(n, &av, &bv, &d, &ex) else {
                return Some((f, r, false));
            };
            if uk.is_zero() {
                ok = false;
                break;
            }
            c = rem(&c.mul(&uk), n);
        }
        if ok {
            let g = c.gcd(n);
            return Some((f, r, g.is_one()));
        }
        // Unused U((n+1)/F) keeps the FLINT structure; a failed q-step
        // tries the next (A, B).
        let _ = um;
        b = b.saturating_add(1);
        if b > 10_000 {
            return Some((f, r, false));
        }
    }
}

fn bls_nminus1_square(n: &ExactInt, f1: &ExactInt) -> Option<bool> {
    // n = c2 F^2 + c1 F + 1; composite if c1^2 - 4 c2 is a square.
    let n1 = n.sub(&one());
    let (n1f, rem1) = n1.div_rem(f1)?;
    if !rem1.is_zero() {
        return None;
    }
    let (c2, c1) = n1f.div_rem(f1)?;
    let disc = c1.mul(&c1).sub(&c2.mul(&ExactInt::from_u64(4)));
    Some(!is_square(&disc))
}

fn combined_or_aprcl(n: &ExactInt, f1: &ExactInt, f2: &ExactInt) -> bool {
    let mut f = f1.mul(f2);
    if is_even(f1) && is_even(f2) {
        f = f.shr_towards_zero(1);
    }
    let fsqr = f.mul(&f);
    if fsqr.soft_cmp(n) == Ordering::Greater {
        let nmod = rem(n, &f);
        if nmod != *n && !nmod.is_one() {
            if let Some((_, r)) = n.div_rem(&nmod) {
                if r.is_zero() {
                    return false;
                }
            }
        }
        return true;
    }
    let fcub = fsqr.mul(&f);
    if fcub.soft_cmp(n) == Ordering::Greater {
        // Lenstra residue-class step is a composite check; if it does not
        // fire we still need a proof, so fall through to APRCL.
        return aprcl_is_prime(n);
    }
    aprcl_is_prime(n)
}

fn trial_smooth_part(m: &ExactInt, primes: &[u64]) -> ExactInt {
    let mut f = ExactInt::one();
    // The cofactor of m; `None` while it is still m.
    let mut rest: Option<ExactInt> = None;
    for &p in primes {
        let pv = ExactInt::from_u64(p);
        loop {
            match rest.as_ref().unwrap_or(m).div_rem(&pv) {
                Some((d, r)) if r.is_zero() => {
                    f = f.mul(&pv);
                    rest = Some(d);
                }
                _ => break,
            }
        }
    }
    f
}

/// Prove an odd Baillie-PSW probable prime, as FLINT 3.6 `fmpz_is_prime` does
/// after the BPSW filter: partial `n-1` / `n+1` BLS, then Jacobi-sum APRCL.
pub(crate) fn prove_prime(n: &ExactInt) -> bool {
    let mut limit = dlog_limit(n);
    for _pass in 0..4 {
        let (pm1, pp1) = trial_pm1_pp1(n, limit);
        // Cheap size check first: only run the Pocklington or Morrison
        // witness search when the trial-smooth part could actually decide.
        // Otherwise FLINT falls through to APRCL, and so do we, without
        // spending a long Lucas-chain search that cannot meet the bound.
        let f1_trial = trial_smooth_part(&n.sub(&one()), &pm1);
        let f2_trial = trial_smooth_part(&n.add(&one()), &pp1);
        let f1_cub = f1_trial.mul(&f1_trial).mul(&f1_trial);
        let f2m1 = if f2_trial.soft_cmp(&one()) == Ordering::Greater {
            f2_trial.sub(&one())
        } else {
            ExactInt::zero()
        };
        let f2_cub = f2m1.mul(&f2m1).mul(&f2m1);
        let mut lcm = f1_trial.mul(&f2_trial);
        if is_even(&f1_trial) && is_even(&f2_trial) {
            lcm = lcm.shr_towards_zero(1);
        }
        let lcm_cub = lcm.mul(&lcm).mul(&lcm);
        let n1_can = f1_cub.soft_cmp(n) != Ordering::Less;
        let n1_sq = f1_trial.mul(&f1_trial).soft_cmp(n) != Ordering::Less;
        let n2_can = f2_cub.soft_cmp(n) == Ordering::Greater
            || f2m1.mul(&f2m1).soft_cmp(n) == Ordering::Greater;
        let comb_can = lcm_cub.soft_cmp(n) == Ordering::Greater
            || lcm.mul(&lcm).soft_cmp(n) == Ordering::Greater;
        if !n1_can && !n2_can && !comb_can {
            return aprcl_is_prime(n);
        }
        let Some((f1, _r1, pok)) = pocklington(n, &pm1) else {
            return false;
        };
        if pok {
            if n1_sq {
                return true;
            }
            if n1_can {
                return bls_nminus1_square(n, &f1).unwrap_or(false);
            }
            if let Some((f2, _r2, mok)) = morrison(n, &pp1) {
                if mok {
                    let fm1 = f2.sub(&one());
                    let fsqr2 = fm1.mul(&fm1);
                    if fsqr2.soft_cmp(n) == Ordering::Greater {
                        return true;
                    }
                    let fcub2 = fsqr2.mul(&fm1);
                    if fcub2.soft_cmp(n) == Ordering::Greater {
                        return aprcl_is_prime(n);
                    }
                    return combined_or_aprcl(n, &f1, &f2);
                }
            }
            return aprcl_is_prime(n);
        }
        // pocklington returned 0 (composite) or we need a larger trial
        // limit. A Fermat failure is a composite; otherwise raise the limit.
        // FLINT only continues the loop while res == -1; a hard 0 stops
        // as composite. Our pocklington returns false both for "no proof"
        // and for a failed Fermat check. If a^{n-1} != 1, n is composite.
        let a = two();
        let nm1 = n.sub(&one());
        match mod_pow(&a, &nm1, n) {
            Some(x) if x.is_one() => {}
            _ => return false,
        }
        match limit.checked_mul(10) {
            Some(v) => limit = v,
            None => break,
        }
    }
    aprcl_is_prime(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn r_value_matches_flint_at_255_bits() {
        let n = ExactInt::from_u64(2).pow(255).sub(&ExactInt::from_u64(19));
        assert_eq!(r_value(n.bit_length()), 2520);
        assert_eq!(n.bit_length(), 255);
    }

    #[test]
    fn limb_add_sub_mul_and_unity_pow() {
        fn to_e(w: &[u64]) -> ExactInt {
            exact_from_limbs(w)
        }
        let mut s = 0x123456789abcdefu64;
        let mut next = || {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
            s
        };
        let n = ExactInt::from_u64(2).pow(255).sub(&ExactInt::from_u64(19));
        let nl = limbs_of(&n);
        for _ in 0..100 {
            let a = rem_limbs(&[next(), next(), next(), next()], nl);
            let b = rem_limbs(&[next(), next(), next(), next()], nl);
            let mut dest = a.clone();
            dest.resize(nl.len(), 0);
            add_mod_limbs(&mut dest, &b, nl);
            let want = rem(&to_e(&a).add(&to_e(&b)), &n);
            assert_eq!(to_e(&dest), want, "add_mod");
            dest = a.clone();
            dest.resize(nl.len(), 0);
            sub_mod_limbs(&mut dest, &b, nl);
            let want = rem(&to_e(&a).sub(&to_e(&b)), &n);
            assert_eq!(to_e(&dest), want, "sub_mod");
            let mut prod = vec![0u64; 10];
            mul_raw(&a, &b, &mut prod);
            let got = rem_limbs(&prod, nl);
            let want = rem(&to_e(&a).mul(&to_e(&b)), &n);
            assert_eq!(to_e(&got), want, "mul_mod");
        }
        let nm1 = n.sub(&ExactInt::one());
        let sq = rem(&nm1.mul(&nm1), &n);
        assert_eq!(sq, ExactInt::one(), "exact (n-1)^2");
        let mut prod = vec![0u64; 10];
        let a = limbs_of(&nm1);
        mul_raw(a, a, &mut prod);
        let got = rem_limbs(&prod, nl);
        assert_eq!(
            exact_from_limbs(&got),
            ExactInt::one(),
            "rem_limbs (n-1)^2 got {got:?} prod={prod:?}"
        );

        // zeta_{32}^{16} == -1 in Z[x]/(x^{16}+1, n)
        // A zero modulus is refused once, at construction (was a debug_assert in rem_limbs).
        assert!(UnityZp::new(2, 5, &ExactInt::zero()).is_none());
        let z = UnityZp::new(2, 5, &n);
        assert!(z.is_some(), "UnityZp::new refused a nonzero modulus");
        let Some(mut z) = z else { return };
        z.set_ui(1, 1);
        let z16 = z.pow_ui(16);
        // expect coeff 0 = n-1
        let c0 = z16
            .slot(0)
            .map(exact_from_limbs)
            .unwrap_or_else(ExactInt::zero);
        assert_eq!(c0, nm1, "zeta^16 coeff0, slots={}", z16.slots);
        for i in 1..z16.slots.max(z16.p_pow()) {
            if let Some(s) = z16.slot(i) {
                assert!(limbs_is_zero(s), "zeta^16 extra coeff {i}");
            }
        }
        let z16s = z16.sqr();
        let s0 = z16s
            .slot(0)
            .map(exact_from_limbs)
            .unwrap_or_else(ExactInt::zero);
        assert_eq!(s0, ExactInt::one(), "(-1)^2 via UnityZp::sqr");
        let z32 = z.pow_ui(32);
        let one_c = z32
            .slot(0)
            .map(exact_from_limbs)
            .unwrap_or_else(ExactInt::zero);
        assert_eq!(one_c, ExactInt::one(), "zeta^32 != 1, slots={}", z32.slots);
        assert_eq!(z32.is_unity(), 0);
        assert_eq!(z16.is_unity(), 16);
    }

    #[test]
    fn rem_limbs_matches_exactint() {
        fn check(u: &ExactInt, v: &ExactInt) {
            let got = rem_limbs(limbs_of(u), limbs_of(v));
            let want = rem(u, v);
            let got_e = exact_from_limbs(&got);
            assert_eq!(got_e, want, "rem({u:?}, {v:?})");
        }
        let n = ExactInt::from_u64(2).pow(255).sub(&ExactInt::from_u64(19));
        check(&ExactInt::zero(), &n);
        check(&n.sub(&ExactInt::one()).mul(&n.sub(&ExactInt::one())), &n);
        check(&n, &n);
        check(&n.sub(&ExactInt::one()), &n);
        check(&n.add(&ExactInt::one()), &n);
        check(&ExactInt::from_u64(1), &n);
        check(&ExactInt::from_u64(u64::MAX), &ExactInt::from_u64(3));
        // random products against 2^255-19 and a few small moduli
        let mut s = 0x9e3779b97f4a7c15u64;
        let mut next = || {
            s = s.wrapping_mul(0x5851f42d4c957f2d).wrapping_add(1);
            s
        };
        for _ in 0..80 {
            let u = ExactInt::from_unsigned_limbs(vec![
                next(),
                next(),
                next(),
                next(),
                next(),
                next(),
                next(),
                next(),
            ]);
            check(&u, &n);
            let v = ExactInt::from_unsigned_limbs(vec![next() | 1, next(), next(), next()]);
            if !v.is_zero() {
                check(&u, &v);
            }
        }
        // one-limb and two-limb divisors
        for d in [2u64, 3, 5, 7, 19, u64::MAX, (1u64 << 32) + 15] {
            check(
                &ExactInt::from_unsigned_limbs(vec![next(), next(), next(), next()]),
                &ExactInt::from_u64(d),
            );
        }
        // 2-limb divisor
        let v2 = ExactInt::from_unsigned_limbs(vec![0x123456789abcdefu64, 0x1111111111111111]);
        check(&ExactInt::from_unsigned_limbs(vec![1, 2, 3, 4]), &v2);

        // Wide numerators (mul accumulators are 2*nlen+2 limbs) and q*n+r identities.
        for _ in 0..200 {
            let q =
                ExactInt::from_unsigned_limbs(vec![next(), next(), next(), next(), next(), next()]);
            let r_lim = vec![next(), next(), next(), next()];
            let r = rem(&ExactInt::from_unsigned_limbs(r_lim), &n);
            let u = q.mul(&n).add(&r);
            check(&u, &n);
            let wide = ExactInt::from_unsigned_limbs(vec![
                next(),
                next(),
                next(),
                next(),
                next(),
                next(),
                next(),
                next(),
                next(),
                next(),
            ]);
            check(&wide, &n);
        }
    }

    #[test]
    fn aprcl_small_primes_and_composites() {
        for p in [5u64, 7, 13, 31, 61, 127, 257, 8191] {
            let v = ExactInt::from_u64(p);
            assert!(aprcl_is_prime(&v), "APRCL must prove {p} prime");
        }
        for c in [9u64, 15, 21, 25, 27, 35, 91, 561, 1105, 2047] {
            let v = ExactInt::from_u64(c);
            assert!(!aprcl_is_prime(&v), "APRCL must reject {c}");
        }
    }

    /// Every odd n in [5, 6000) against trial division, and multi-limb Mersenne primes and their products
    /// (the Jacobi-sum checks run in Z[zeta]/n with n of 1, 2 and 4 limbs).
    #[test]
    fn aprcl_agrees_with_trial_division_and_mersenne() -> Result<(), String> {
        let trial = |n: u64| {
            (2..)
                .take_while(|d| d * d <= n)
                .all(|d| !n.is_multiple_of(d))
        };
        for n in (5u64..6000).step_by(2) {
            if aprcl_is_prime(&ExactInt::from_u64(n)) != trial(n) {
                return Err(format!("n={n}"));
            }
        }
        let mersenne = |e: u32| {
            let mut w = vec![u64::MAX; (e / 64) as usize];
            if !e.is_multiple_of(64) {
                w.push((1u64 << (e % 64)) - 1);
            }
            ExactInt::from_unsigned_limbs(w)
        };
        for e in [61u32, 89, 107, 127, 521] {
            if !aprcl_is_prime(&mersenne(e)) {
                return Err(format!("2^{e}-1 not proved prime"));
            }
        }
        for (a, b) in [(61u32, 31u32), (89, 61), (127, 107), (127, 89)] {
            if aprcl_is_prime(&mersenne(a).mul(&mersenne(b))) {
                return Err(format!("(2^{a}-1)(2^{b}-1) accepted"));
            }
        }
        if aprcl_is_prime(&mersenne(67)) {
            return Err("2^67-1 (composite) accepted".to_string());
        }
        Ok(())
    }
}
