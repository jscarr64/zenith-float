//! Error function, gamma, and Bessel J_n (integer order).

use crate::common::util::bump_prec_retry;
use crate::common::util::round_p;
use crate::defs::Error;
use crate::defs::RoundingMode;
use crate::num::ExactNumNumber;
use crate::ops::consts::Consts;
use crate::Sign;
use crate::WORD_BIT_SIZE;
use crate::Word;

use alloc::vec::Vec;

impl ExactNumNumber {
    /// Error function `erf(self)` at precision `p`.
    pub fn erf(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;

        if self.is_zero() {
            return Self::new2(p, self.sign(), self.inexact());
        }

        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc;

        loop {
            let p_x = p_wrk + WORD_BIT_SIZE;
            let mut ret = self.erf_at(p_x, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    /// Complementary error function `erfc(self) = 1 - erf(self)` at precision `p`.
    pub fn erfc(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;

        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc;

        loop {
            let p_x = p_wrk + WORD_BIT_SIZE;
            let one = Self::from_word(1, p_x)?;
            let erf = self.erf_at(p_x, cc)?;
            let mut ret = one.sub(&erf, p_x, RoundingMode::None)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn erf_at(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let mut x = self.clone()?;
        x.set_inexact(false);
        let neg = x.is_negative();
        x.set_sign(Sign::Pos);

        // |x| large enough that erfc underflows at this precision.
        if (x.exponent() as isize).saturating_mul(2) > p as isize + 4 {
            let mut one = Self::from_word(1, p)?;
            if neg {
                one.set_sign(Sign::Neg);
            }
            return Ok(one);
        }

        let ret = if x.exponent() <= 2 {
            x.erf_series(p, cc)?
        } else {
            let erfc = x.erfc_asymptotic(p, cc)?;
            let one = Self::from_word(1, p)?;
            one.sub(&erfc, p, RoundingMode::None)?
        };

        if neg {
            ret.neg()
        } else {
            Ok(ret)
        }
    }

    fn erf_series(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        // (2/√π) ∑ (-1)^n x^{2n+1} / (n! (2n+1))
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let sqrt_pi = pi.sqrt(p, RoundingMode::None)?;
        let two = Self::from_word(2, p)?;
        let scale = two.div(&sqrt_pi, p, RoundingMode::None)?;

        let x2 = self.mul(self, p, RoundingMode::None)?;
        let mut xpow = self.clone()?;
        let mut nfact = Self::from_word(1, p)?;
        let mut sum = xpow.clone()?;

        for n in 1..=(p + 8) {
            nfact = nfact.mul(&Self::from_word(n as Word, p)?, p, RoundingMode::None)?;
            xpow = xpow.mul(&x2, p, RoundingMode::None)?;
            let two_n_1 = Self::from_word((2 * n + 1) as Word, p)?;
            let den = nfact.mul(&two_n_1, p, RoundingMode::None)?;
            let mut t = xpow.div(&den, p, RoundingMode::None)?;
            if n % 2 == 1 {
                t.set_sign(Sign::Neg);
            }
            sum = sum.add(&t, p, RoundingMode::None)?;
            if t.is_zero() || (t.exponent() as isize) + (p as isize) < 0 {
                break;
            }
        }

        scale.mul(&sum, p, RoundingMode::None)
    }

    fn erfc_asymptotic(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        // e^{-x^2} / (x √π) * ∑ (-1)^k (2k-1)!! / (2x^2)^k
        let x2 = self.mul(self, p, RoundingMode::None)?;
        let mut nx2 = x2.clone()?;
        nx2.set_sign(Sign::Neg);
        let expm = nx2.exp(p, RoundingMode::None, cc)?;
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let sqrt_pi = pi.sqrt(p, RoundingMode::None)?;
        let den = self.mul(&sqrt_pi, p, RoundingMode::None)?;
        let pre = expm.div(&den, p, RoundingMode::None)?;

        let two_x2 = x2.mul(&Self::from_word(2, p)?, p, RoundingMode::None)?;
        let mut s = Self::from_word(1, p)?;
        let mut term = Self::from_word(1, p)?;
        for k in 1..=(p + 8) {
            let odd = Self::from_word((2 * k - 1) as Word, p)?;
            term = term.mul(&odd, p, RoundingMode::None)?;
            term = term.div(&two_x2, p, RoundingMode::None)?;
            term.set_sign(if k % 2 == 1 { Sign::Neg } else { Sign::Pos });
            s = s.add(&term, p, RoundingMode::None)?;
            if term.is_zero() || (term.exponent() as isize) + (p as isize) < 0 {
                break;
            }
            term.set_sign(Sign::Pos);
        }
        pre.mul(&s, p, RoundingMode::None)
    }

    /// Gamma function `Γ(self)` at precision `p`.
    pub fn gamma(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;

        if self.is_zero() {
            return Err(Error::ExponentOverflow(Sign::Pos));
        }

        let fr = self.fract()?;
        if !self.is_positive() && fr.is_zero() {
            return Err(Error::InvalidArgument);
        }

        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc;

        loop {
            let p_x = p_wrk + WORD_BIT_SIZE * 2;
            let mut ret = self.gamma_at(p_x, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    /// `ln Γ(self)` for positive `self` at precision `p`.
    pub fn ln_gamma(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;

        if !self.is_positive() {
            return Err(Error::InvalidArgument);
        }

        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc;

        loop {
            let p_x = p_wrk + WORD_BIT_SIZE * 2;
            let g = self.gamma_at(p_x, cc)?;
            let mut ret = g.ln(p_x, RoundingMode::None, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn gamma_at(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let mut z = self.clone()?;
        z.set_inexact(false);

        if !z.is_positive() {
            // Reflection: Γ(z) = π / (sin(πz) Γ(1-z))
            let pi = cc.pi_num(p, RoundingMode::None)?;
            let piz = pi.mul(&z, p, RoundingMode::None)?;
            let s = piz.sin(p, RoundingMode::None, cc)?;
            let one = Self::from_word(1, p)?;
            let omz = one.sub(&z, p, RoundingMode::None)?;
            let g = omz.gamma_positive(p, cc)?;
            let den = s.mul(&g, p, RoundingMode::None)?;
            return pi.div(&den, p, RoundingMode::None);
        }

        z.gamma_positive(p, cc)
    }

    fn gamma_positive(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        if self.fract()?.is_zero() && self.is_positive() {
            for n in 1u32..=64 {
                let w = Self::from_word(n as Word, p)?;
                if self.cmp(&w) == 0 {
                    return factorial((n - 1) as usize, p);
                }
            }
        }

        let one = Self::from_word(1, p)?;
        let mut z = self.clone()?;
        let mut acc = Self::from_word(1, p)?;
        let mut shifts = 0usize;

        // Raise argument until |z| >= 64 so Stirling terms decay quickly.
        while z.exponent() < 6 && shifts < 512 {
            acc = acc.mul(&z, p, RoundingMode::None)?;
            z = z.add(&one, p, RoundingMode::None)?;
            shifts += 1;
        }

        let lg = z.ln_gamma_stirling(p, cc)?;
        let g = lg.exp(p, RoundingMode::None, cc)?;
        g.div(&acc, p, RoundingMode::None)
    }

    fn ln_gamma_stirling(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        // (z-1/2) ln z - z + (1/2) ln(2π) + ∑ B_{2k} / (2k(2k-1) z^{2k-1})
        let ln_z = self.ln(p, RoundingMode::None, cc)?;
        let two = Self::from_word(2, p)?;
        let half = one_half(p)?;
        let zmh = self.sub(&half, p, RoundingMode::None)?;
        let mut s = zmh.mul(&ln_z, p, RoundingMode::None)?;
        s = s.sub(self, p, RoundingMode::None)?;

        let pi = cc.pi_num(p, RoundingMode::None)?;
        let two_pi = two.mul(&pi, p, RoundingMode::None)?;
        let mut ln_two_pi = two_pi.ln(p, RoundingMode::None, cc)?;
        ln_two_pi.div_by_2(RoundingMode::None);
        s = s.add(&ln_two_pi, p, RoundingMode::None)?;

        let mut zpow = self.clone()?; // z^{1} for k=1 term z^{2k-1}
        for k in 1..=64 {
            let b = bernoulli_even(k, p)?;
            let two_k = Self::from_word((2 * k) as Word, p)?;
            let two_k_m1 = Self::from_word((2 * k - 1) as Word, p)?;
            let den = two_k
                .mul(&two_k_m1, p, RoundingMode::None)?
                .mul(&zpow, p, RoundingMode::None)?;
            let term = b.div(&den, p, RoundingMode::None)?;
            s = s.add(&term, p, RoundingMode::None)?;
            if term.is_zero() || (term.exponent() as isize) + (p as isize) < 0 {
                break;
            }
            zpow = zpow.mul(self, p, RoundingMode::None)?;
            zpow = zpow.mul(self, p, RoundingMode::None)?;
        }
        Ok(s)
    }

    /// Bessel function of the first kind `J_n(self)` for integer order `n`.
    pub fn bessel_j(
        &self,
        n: usize,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let _ = cc;
        let p = round_p(p);
        Self::p_assertion(p)?;

        if n > 1024 {
            return Err(Error::InvalidArgument);
        }

        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc;

        loop {
            let p_x = p_wrk + WORD_BIT_SIZE;
            let mut ret = self.bessel_j_series(n, p_x)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn bessel_j_series(&self, n: usize, p: usize) -> Result<Self, Error> {
        // J_n(x) = (x/2)^n ∑_{k=0} (-1)^k / (k! (k+n)!) (x/2)^{2k}
        let two = Self::from_word(2, p)?;
        let xh = self.div(&two, p, RoundingMode::None)?;
        let mut pow = Self::from_word(1, p)?;
        for _ in 0..n {
            pow = pow.mul(&xh, p, RoundingMode::None)?;
        }

        let mut kfact = Self::from_word(1, p)?;
        let mut knfact = factorial(n, p)?;
        let mut sum = pow.div(&knfact, p, RoundingMode::None)?;
        let xh2 = xh.mul(&xh, p, RoundingMode::None)?;
        let mut num = pow;

        for k in 1..=(p + n + 8) {
            kfact = kfact.mul(&Self::from_word(k as Word, p)?, p, RoundingMode::None)?;
            knfact = knfact.mul(&Self::from_word((k + n) as Word, p)?, p, RoundingMode::None)?;
            num = num.mul(&xh2, p, RoundingMode::None)?;
            let den = kfact.mul(&knfact, p, RoundingMode::None)?;
            let mut t = num.div(&den, p, RoundingMode::None)?;
            if k % 2 == 1 {
                t.set_sign(if t.is_positive() { Sign::Neg } else { Sign::Pos });
            }
            sum = sum.add(&t, p, RoundingMode::None)?;
            if t.is_zero() || (t.exponent() as isize) + (p as isize) < 0 {
                break;
            }
        }
        Ok(sum)
    }
}

fn one_half(p: usize) -> Result<ExactNumNumber, Error> {
    let mut h = ExactNumNumber::from_word(1, p)?;
    h.div_by_2(RoundingMode::None);
    Ok(h)
}

fn factorial(n: usize, p: usize) -> Result<ExactNumNumber, Error> {
    let mut acc = ExactNumNumber::from_word(1, p)?;
    for i in 2..=n {
        acc = acc.mul(
            &ExactNumNumber::from_word(i as Word, p)?,
            p,
            RoundingMode::None,
        )?;
    }
    Ok(acc)
}

/// Even Bernoulli number B_{2k} via Akiyama–Tanigawa.
fn bernoulli_even(k: usize, p: usize) -> Result<ExactNumNumber, Error> {
    let m = 2 * k;
    let mut a: Vec<ExactNumNumber> = Vec::new();
    for i in 0..=m {
        let num = ExactNumNumber::from_word(1, p)?;
        let den = ExactNumNumber::from_word((i + 1) as Word, p)?;
        a.push(num.div(&den, p, RoundingMode::None)?);
    }
    for j in 1..=m {
        for i in 0..=(m - j) {
            let diff = a[i].sub(&a[i + 1], p, RoundingMode::None)?;
            let fac = ExactNumNumber::from_word((i + 1) as Word, p)?;
            a[i] = fac.mul(&diff, p, RoundingMode::None)?;
        }
    }
    a[0].clone()
}

/// Euler–Mascheroni constant via `H_{n-1} - ln n + 1/(2n) + Σ B_{2k}/(2k n^{2k})` with `n = 128`.
pub(crate) fn euler_mascheroni(p: usize, ln2: &ExactNumNumber) -> Result<ExactNumNumber, Error> {
    const N: usize = 128;
    let mut h = ExactNumNumber::from_word(1, p)?;
    for i in 2..N {
        let t = ExactNumNumber::from_word(1, p)?.div(
            &ExactNumNumber::from_word(i as Word, p)?,
            p,
            RoundingMode::None,
        )?;
        h = h.add(&t, p, RoundingMode::None)?;
    }
    let ln_n = ExactNumNumber::from_word(7, p)?.mul(ln2, p, RoundingMode::None)?;
    let mut g = h.sub(&ln_n, p, RoundingMode::None)?;
    let n = ExactNumNumber::from_word(N as Word, p)?;
    let two_n = n.add(&n, p, RoundingMode::None)?;
    let half_n = ExactNumNumber::from_word(1, p)?.div(&two_n, p, RoundingMode::None)?;
    g = g.add(&half_n, p, RoundingMode::None)?;

    let mut npow = n.mul(&n, p, RoundingMode::None)?;
    for k in 1..=40 {
        let b = bernoulli_even(k, p)?;
        let two_k = ExactNumNumber::from_word((2 * k) as Word, p)?;
        let den = two_k.mul(&npow, p, RoundingMode::None)?;
        let term = b.div(&den, p, RoundingMode::None)?;
        g = g.add(&term, p, RoundingMode::None)?;
        if term.is_zero() || (term.exponent() as isize) + (p as isize) < 0 {
            break;
        }
        npow = npow.mul(&n, p, RoundingMode::None)?;
        npow = npow.mul(&n, p, RoundingMode::None)?;
    }
    Ok(g)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_erf_gamma_bessel() {
        let p = 256;
        let mut cc = Consts::new().unwrap();
        let rm = RoundingMode::ToEven;

        let zero = ExactNumNumber::new(p).unwrap();
        assert!(zero.erf(p, rm, &mut cc).unwrap().is_zero());

        let one = ExactNumNumber::from_word(1, p).unwrap();
        let g1 = one.gamma(p, rm, &mut cc).unwrap();
        let d = g1.sub(&one, p, RoundingMode::None).unwrap().abs().unwrap();
        assert!(d.is_zero() || d.exponent() < -((p as i32) / 4));

        let five = ExactNumNumber::from_word(5, p).unwrap();
        let g5 = five.gamma(p, rm, &mut cc).unwrap();
        let tf = ExactNumNumber::from_word(24, p).unwrap();
        let d = g5.sub(&tf, p, RoundingMode::None).unwrap().abs().unwrap();
        assert!(d.is_zero() || d.exponent() < -((p as i32) / 4));

        let j0 = zero.bessel_j(0, p, rm, &mut cc).unwrap();
        let d = j0.sub(&one, p, RoundingMode::None).unwrap().abs().unwrap();
        assert!(d.is_zero() || d.exponent() < -((p as i32) / 4));

        let x = ExactNumNumber::from_word(1, p).unwrap();
        let e = x.erf(p, rm, &mut cc).unwrap();
        let en = x.neg().unwrap().erf(p, rm, &mut cc).unwrap();
        assert!(e.cmp(&en.neg().unwrap()) == 0);
    }
}
