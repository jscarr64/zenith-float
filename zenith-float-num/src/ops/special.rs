//! Error function, gamma, integral specials, and Bessel J_n (integer order).

use crate::common::util::bump_prec_retry;
use crate::common::util::round_p;
use crate::defs::Error;
use crate::defs::RoundingMode;
use crate::num::ExactNumNumber;
use crate::ops::consts::Consts;
use crate::Sign;
use crate::Word;
use crate::WORD_BIT_SIZE;

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
                    let mut f = factorial((n - 1) as usize, p)?;
                    f.set_inexact(false);
                    return Ok(f);
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
            let den =
                two_k
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

    /// Digamma \(\psi(\mathrm{self})=\Gamma'/\Gamma\). Poles at non-positive integers.
    /// For \(z<0\) uses \(\psi(z)=\psi(1-z)-\pi\cot(\pi z)\).
    pub fn digamma(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        if !self.is_positive() && self.fract()?.is_zero() {
            return Err(Error::InvalidArgument);
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE * 2;
            let mut ret = self.digamma_at(p_x, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn digamma_at(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let one = Self::from_word(1, p)?;
        let mut z = self.clone()?;
        z.set_inexact(false);
        if !z.is_positive() {
            // ψ(z) = ψ(1−z) − π cot(πz)
            let omz = one.sub(&z, p, RoundingMode::None)?;
            let psi = omz.digamma_positive(p, cc)?;
            let pi = cc.pi_num(p, RoundingMode::None)?;
            let piz = pi.mul(&z, p, RoundingMode::None)?;
            let cot = piz.cos(p, RoundingMode::None, cc)?.div(
                &piz.sin(p, RoundingMode::None, cc)?,
                p,
                RoundingMode::None,
            )?;
            return psi.sub(&pi.mul(&cot, p, RoundingMode::None)?, p, RoundingMode::None);
        }
        z.digamma_positive(p, cc)
    }

    fn digamma_positive(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let one = Self::from_word(1, p)?;
        let mut z = self.clone()?;
        z.set_inexact(false);
        let mut acc = Self::from_word(0, p)?;
        let mut shifts = 0usize;
        // Raise until |z| >= 128 so Bernoulli terms reach typical p.
        while z.exponent() < 8 && shifts < 512 {
            let rec = one.div(&z, p, RoundingMode::None)?;
            acc = acc.sub(&rec, p, RoundingMode::None)?;
            z = z.add(&one, p, RoundingMode::None)?;
            shifts += 1;
        }
        acc.add(&z.digamma_asymp(p, cc)?, p, RoundingMode::None)
    }

    fn digamma_asymp(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        // ln z − 1/(2z) − Σ_{k≥1} B_{2k} / (2k z^{2k})
        let ln_z = self.ln(p, RoundingMode::None, cc)?;
        let two = Self::from_word(2, p)?;
        let two_z = two.mul(self, p, RoundingMode::None)?;
        let half_inv = Self::from_word(1, p)?.div(&two_z, p, RoundingMode::None)?;
        let mut s = ln_z.sub(&half_inv, p, RoundingMode::None)?;
        let z2 = self.mul(self, p, RoundingMode::None)?;
        let mut zp = Self::from_word(1, p)?;
        let mut prev_e = i32::MIN;
        for k in 1..=64 {
            zp = zp.mul(&z2, p, RoundingMode::None)?;
            let b = bernoulli_even(k, p)?;
            let two_k = Self::from_word((2 * k) as Word, p)?;
            let den = two_k.mul(&zp, p, RoundingMode::None)?;
            let term = b.div(&den, p, RoundingMode::None)?;
            s = s.sub(&term, p, RoundingMode::None)?;
            let e = term.exponent();
            if term.is_zero() || (e as isize) + (p as isize) < 0 {
                break;
            }
            if k > 2 && e > prev_e {
                break;
            }
            prev_e = e;
        }
        Ok(s)
    }

    /// Lower incomplete gamma \(\gamma(s=\mathrm{self}, x)\) for `self > 0`, `x >= 0`.
    pub fn gammainc(
        &self,
        x: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        if !self.is_positive() || x.is_negative() {
            return Err(Error::InvalidArgument);
        }
        if x.is_zero() {
            return Self::new2(p, Sign::Pos, self.inexact() | x.inexact());
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p
            .max(self.mantissa_max_bit_len())
            .max(x.mantissa_max_bit_len())
            + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE * 2;
            let mut ret = self.gammainc_at(x, p_x, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact() | x.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    /// Upper incomplete gamma \(\Gamma(s=\mathrm{self}, x)=\Gamma(s)-\gamma(s,x)\).
    pub fn gammainc_upper(
        &self,
        x: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        if !self.is_positive() || x.is_negative() {
            return Err(Error::InvalidArgument);
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p
            .max(self.mantissa_max_bit_len())
            .max(x.mantissa_max_bit_len())
            + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE * 2;
            let g = self.gamma_at(p_x, cc)?;
            let lo = self.gammainc_at(x, p_x, cc)?;
            let mut ret = g.sub(&lo, p_x, RoundingMode::None)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact() | x.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn gammainc_at(&self, x: &Self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        // e^{-x} x^s underflows → γ(s,x) = Γ(s).
        if (x.exponent() as isize).saturating_mul(2) > p as isize + 8 {
            return self.gamma_at(p, cc);
        }
        // γ(s,x) = x^s e^{-x} Σ_{k=0} x^k / (s)_{k+1}
        let xs = x.pow(self, p, RoundingMode::None, cc)?;
        let mut nx = x.clone()?;
        nx.set_sign(Sign::Neg);
        let exm = nx.exp(p, RoundingMode::None, cc)?;
        let pre = xs.mul(&exm, p, RoundingMode::None)?;
        let mut poch = self.clone()?;
        let mut term = Self::from_word(1, p)?.div(&poch, p, RoundingMode::None)?;
        let mut sum = term.clone()?;
        let one = Self::from_word(1, p)?;
        for _k in 1..=series_n_max(p, x.exponent()) {
            poch = poch.add(&one, p, RoundingMode::None)?;
            term = term.mul(x, p, RoundingMode::None)?;
            term = term.div(&poch, p, RoundingMode::None)?;
            sum = sum.add(&term, p, RoundingMode::None)?;
            if term.is_zero() || (term.exponent() as isize) + (p as isize) < 0 {
                break;
            }
        }
        pre.mul(&sum, p, RoundingMode::None)
    }

    /// Exponential integral `Ei(self)` (Cauchy PV for `self < 0`). `self = 0` is a pole.
    pub fn ei(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        if self.is_zero() {
            return Err(Error::InvalidArgument);
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE;
            let mut ret = self.ei_at(p_x, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn ei_at(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        if self.expint_use_asymptotic(p) {
            self.ei_asymptotic(p, cc)
        } else {
            self.ei_series(p, cc)
        }
    }

    fn ei_series(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        // γ + ln x + Σ_{n=1}^∞ x^n / (n n!)
        let g = cc.euler_gamma_num(p, RoundingMode::None)?;
        let lnx = self.abs()?.ln(p, RoundingMode::None, cc)?;
        let mut term = self.clone()?;
        let mut sum = term.clone()?;
        for n in 2..=series_n_max(p, self.exponent()) {
            let nw = Self::from_word(n as Word, p)?;
            term = term.mul(self, p, RoundingMode::None)?;
            term = term.div(&nw, p, RoundingMode::None)?;
            let piece = term.div(&nw, p, RoundingMode::None)?;
            sum = sum.add(&piece, p, RoundingMode::None)?;
            if piece.is_zero() || (piece.exponent() as isize) + (p as isize) < 0 {
                break;
            }
        }
        g.add(&lnx, p, RoundingMode::None)?
            .add(&sum, p, RoundingMode::None)
    }

    fn ei_asymptotic(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        // e^x / x · Σ_{k=0} k! / x^k (stop when terms grow)
        let ex = self.exp(p, RoundingMode::None, cc)?;
        let pre = ex.div(self, p, RoundingMode::None)?;
        let mut term = Self::from_word(1, p)?;
        let mut s = term.clone()?;
        let mut prev_e = i32::MAX;
        for k in 1..=(p + 8) {
            term = term.mul(&Self::from_word(k as Word, p)?, p, RoundingMode::None)?;
            term = term.div(self, p, RoundingMode::None)?;
            let e = term.exponent();
            if e > prev_e {
                break;
            }
            prev_e = e;
            s = s.add(&term, p, RoundingMode::None)?;
            if term.is_zero() || (e as isize) + (p as isize) < 0 {
                break;
            }
        }
        pre.mul(&s, p, RoundingMode::None)
    }

    /// Sine integral `Si(self)` for all real `self`.
    pub fn si(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        if self.is_zero() {
            return Self::new2(p, self.sign(), self.inexact());
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE;
            let mut ret = self.si_at(p_x, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn si_at(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let mut x = self.clone()?;
        x.set_inexact(false);
        let neg = x.is_negative();
        x.set_sign(Sign::Pos);
        let ret = if x.expint_use_asymptotic(p) {
            x.si_asymptotic(p, cc)?
        } else {
            x.si_series(p)?
        };
        if neg {
            ret.neg()
        } else {
            Ok(ret)
        }
    }

    fn si_series(&self, p: usize) -> Result<Self, Error> {
        let x2 = self.mul(self, p, RoundingMode::None)?;
        let mut t = self.clone()?;
        let mut sum = t.clone()?;
        for n in 1..=series_n_max(p, self.exponent()) {
            let two_n = Self::from_word((2 * n) as Word, p)?;
            let two_n_1 = Self::from_word((2 * n + 1) as Word, p)?;
            t = t.mul(&x2, p, RoundingMode::None)?;
            t = t.div(&two_n, p, RoundingMode::None)?;
            t = t.div(&two_n_1, p, RoundingMode::None)?;
            let mut piece = t.div(&two_n_1, p, RoundingMode::None)?;
            if n % 2 == 1 {
                piece.set_sign(Sign::Neg);
            }
            sum = sum.add(&piece, p, RoundingMode::None)?;
            if piece.is_zero() || (piece.exponent() as isize) + (p as isize) < 0 {
                break;
            }
        }
        Ok(sum)
    }

    fn si_asymptotic(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        // π/2 − f(x) cos x − g(x) sin x
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let mut half_pi = pi.clone()?;
        half_pi.div_by_2(RoundingMode::None);
        let (f, g) = self.si_ci_aux_fg(p)?;
        let c = self.cos(p, RoundingMode::None, cc)?;
        let s = self.sin(p, RoundingMode::None, cc)?;
        let fc = f.mul(&c, p, RoundingMode::None)?;
        let gs = g.mul(&s, p, RoundingMode::None)?;
        half_pi
            .sub(&fc, p, RoundingMode::None)?
            .sub(&gs, p, RoundingMode::None)
    }

    /// Cosine integral `Ci(self)` for `self > 0`.
    pub fn ci(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        if !self.is_positive() {
            return Err(Error::InvalidArgument);
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE;
            let mut ret = self.ci_at(p_x, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn ci_at(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        if self.expint_use_asymptotic(p) {
            self.ci_asymptotic(p, cc)
        } else {
            self.ci_series(p, cc)
        }
    }

    fn ci_series(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        // γ + ln x + Σ_{n=1}^∞ (−1)^n x^{2n} / (2n (2n)!)
        let g = cc.euler_gamma_num(p, RoundingMode::None)?;
        let lnx = self.ln(p, RoundingMode::None, cc)?;
        let x2 = self.mul(self, p, RoundingMode::None)?;
        let two = Self::from_word(2, p)?;
        let mut u = x2.div(&two, p, RoundingMode::None)?;
        u.set_sign(Sign::Neg);
        let mut sum = u.div(&two, p, RoundingMode::None)?;
        for n in 2..=series_n_max(p, self.exponent()) {
            let a = Self::from_word((2 * n - 1) as Word, p)?;
            let b = Self::from_word((2 * n) as Word, p)?;
            u = u.mul(&x2, p, RoundingMode::None)?;
            u = u.div(&a, p, RoundingMode::None)?;
            u = u.div(&b, p, RoundingMode::None)?;
            u = u.neg()?;
            let piece = u.div(&b, p, RoundingMode::None)?;
            sum = sum.add(&piece, p, RoundingMode::None)?;
            if piece.is_zero() || (piece.exponent() as isize) + (p as isize) < 0 {
                break;
            }
        }
        g.add(&lnx, p, RoundingMode::None)?
            .add(&sum, p, RoundingMode::None)
    }

    fn ci_asymptotic(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        // f(x) sin x − g(x) cos x
        let (f, g) = self.si_ci_aux_fg(p)?;
        let s = self.sin(p, RoundingMode::None, cc)?;
        let c = self.cos(p, RoundingMode::None, cc)?;
        let fs = f.mul(&s, p, RoundingMode::None)?;
        let gc = g.mul(&c, p, RoundingMode::None)?;
        fs.sub(&gc, p, RoundingMode::None)
    }

    /// Auxiliary `f,g` for `Si`/`Ci`: `f ∼ Σ (−1)^k (2k)! / x^{2k+1}`, `g ∼ Σ (−1)^k (2k+1)! / x^{2k+2}`.
    fn si_ci_aux_fg(&self, p: usize) -> Result<(Self, Self), Error> {
        let one = Self::from_word(1, p)?;
        let invx = one.div(self, p, RoundingMode::None)?;
        let invx2 = invx.mul(&invx, p, RoundingMode::None)?;
        let mut f_term = invx;
        let mut g_term = invx2.clone()?;
        let mut f = f_term.clone()?;
        let mut g = g_term.clone()?;
        let mut prev_e = f_term.exponent();
        for k in 0..=(p + 8) {
            let a = Self::from_word((2 * k + 1) as Word, p)?;
            let b = Self::from_word((2 * k + 2) as Word, p)?;
            let c = Self::from_word((2 * k + 3) as Word, p)?;
            f_term = f_term
                .mul(&a, p, RoundingMode::None)?
                .mul(&b, p, RoundingMode::None)?
                .mul(&invx2, p, RoundingMode::None)?
                .neg()?;
            g_term = g_term
                .mul(&b, p, RoundingMode::None)?
                .mul(&c, p, RoundingMode::None)?
                .mul(&invx2, p, RoundingMode::None)?
                .neg()?;
            let e = f_term.exponent();
            if e > prev_e {
                break;
            }
            prev_e = e;
            f = f.add(&f_term, p, RoundingMode::None)?;
            g = g.add(&g_term, p, RoundingMode::None)?;
            if f_term.is_zero() || (e as isize) + (p as isize) < 0 {
                break;
            }
        }
        Ok((f, g))
    }

    /// Asymptotic `e^{-|x|}` remainder is below `2^{-p}` when `|x| ≳ 0.7 p`.
    fn expint_use_asymptotic(&self, p: usize) -> bool {
        let e = self.exponent();
        if e <= 6 {
            return false;
        }
        if e >= 20 {
            return true;
        }
        let xmin = 1i32 << (e - 1);
        xmin >= ((p as i32) * 2) / 3
    }

    /// Logarithmic integral `li(self) = Ei(ln self)` for `self > 0`, `self ≠ 1`.
    pub fn li(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        let one = Self::from_word(1, p)?;
        if !self.is_positive() || self.cmp(&one) == 0 {
            return Err(Error::InvalidArgument);
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE;
            let lnx = self.ln(p_x, RoundingMode::None, cc)?;
            let mut ret = lnx.ei(p_x, RoundingMode::None, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    /// Fresnel sine integral `S(self)`.
    pub fn fresnel_s(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        self.fresnel_sc(true, p, rm, cc)
    }

    /// Fresnel cosine integral `C(self)`.
    pub fn fresnel_c(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        self.fresnel_sc(false, p, rm, cc)
    }

    fn fresnel_sc(
        &self,
        sine: bool,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        if self.is_zero() {
            return Self::new2(p, self.sign(), self.inexact());
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE;
            let mut ret = self.fresnel_sc_at(sine, p_x, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn fresnel_sc_at(&self, sine: bool, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let mut x = self.clone()?;
        x.set_inexact(false);
        let neg = x.is_negative();
        x.set_sign(Sign::Pos);
        let ret = if x.fresnel_use_asymptotic(p) {
            x.fresnel_sc_asymptotic(sine, p, cc)?
        } else {
            x.fresnel_sc_series(sine, p, cc)?
        };
        if neg {
            ret.neg()
        } else {
            Ok(ret)
        }
    }

    fn fresnel_use_asymptotic(&self, p: usize) -> bool {
        let e = self.exponent();
        if e <= 3 {
            return false;
        }
        if e >= 8 {
            return true;
        }
        let xmin = 1i64 << (e - 1);
        xmin * xmin * 3 > p as i64
    }

    fn fresnel_sc_asymptotic(&self, sine: bool, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        // S = 1/2 − f cos(πx²/2) − g sin(πx²/2)
        // C = 1/2 + f sin(πx²/2) − g cos(πx²/2)
        let (f, g) = self.fresnel_aux_fg(p, cc)?;
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let two = Self::from_word(2, p)?;
        let pi2 = pi.div(&two, p, RoundingMode::None)?;
        let x2 = self.mul(self, p, RoundingMode::None)?;
        let arg = pi2.mul(&x2, p, RoundingMode::None)?;
        let s = arg.sin(p, RoundingMode::None, cc)?;
        let c = arg.cos(p, RoundingMode::None, cc)?;
        let half = one_half(p)?;
        let fc = f.mul(&c, p, RoundingMode::None)?;
        let gs = g.mul(&s, p, RoundingMode::None)?;
        if sine {
            half.sub(&fc, p, RoundingMode::None)?
                .sub(&gs, p, RoundingMode::None)
        } else {
            let fs = f.mul(&s, p, RoundingMode::None)?;
            let gc = g.mul(&c, p, RoundingMode::None)?;
            half.add(&fs, p, RoundingMode::None)?
                .sub(&gc, p, RoundingMode::None)
        }
    }

    fn fresnel_aux_fg(&self, p: usize, cc: &mut Consts) -> Result<(Self, Self), Error> {
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let pix = pi.mul(self, p, RoundingMode::None)?;
        let one = Self::from_word(1, p)?;
        let inv_pix = one.div(&pix, p, RoundingMode::None)?;
        let pix2 = pix.mul(self, p, RoundingMode::None)?;
        let v = one.div(&pix2, p, RoundingMode::None)?;
        let v2 = v.mul(&v, p, RoundingMode::None)?;
        let mut f_term = Self::from_word(1, p)?;
        let mut g_term = v;
        let mut f_sum = f_term.clone()?;
        let mut g_sum = g_term.clone()?;
        let mut prev_e = f_term.exponent();
        for n in 0..=(p + 8) {
            let a = Self::from_word((4 * n + 1) as Word, p)?;
            let b = Self::from_word((4 * n + 3) as Word, p)?;
            let c = Self::from_word((4 * n + 5) as Word, p)?;
            f_term = f_term
                .mul(&a, p, RoundingMode::None)?
                .mul(&b, p, RoundingMode::None)?
                .mul(&v2, p, RoundingMode::None)?
                .neg()?;
            g_term = g_term
                .mul(&b, p, RoundingMode::None)?
                .mul(&c, p, RoundingMode::None)?
                .mul(&v2, p, RoundingMode::None)?
                .neg()?;
            let e = f_term.exponent();
            if e > prev_e {
                break;
            }
            prev_e = e;
            f_sum = f_sum.add(&f_term, p, RoundingMode::None)?;
            g_sum = g_sum.add(&g_term, p, RoundingMode::None)?;
            if f_term.is_zero() || (e as isize) + (p as isize) < 0 {
                break;
            }
        }
        Ok((
            inv_pix.mul(&f_sum, p, RoundingMode::None)?,
            inv_pix.mul(&g_sum, p, RoundingMode::None)?,
        ))
    }

    fn fresnel_sc_series(&self, sine: bool, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let two = Self::from_word(2, p)?;
        let pi2 = pi.div(&two, p, RoundingMode::None)?;
        let x2 = self.mul(self, p, RoundingMode::None)?;
        let x4 = x2.mul(&x2, p, RoundingMode::None)?;
        let mut pow_x = if sine { self.mul(&x2, p, RoundingMode::None)? } else { self.clone()? };
        let mut fact = Self::from_word(1, p)?;
        let mut pi2_pow = if sine { pi2.clone()? } else { Self::from_word(1, p)? };
        let mut sign_pos = true;
        let mut sum = Self::from_word(1, p)?;
        sum.set_sign(Sign::Pos);
        let mut first = true;
        for n in 0..=series_n_max(p, self.exponent()) {
            let den_i = if sine { 4 * n + 3 } else { 4 * n + 1 };
            let den = Self::from_word(den_i as Word, p)?;
            let mut piece = pi2_pow.mul(&pow_x, p, RoundingMode::None)?;
            piece = piece.div(&fact, p, RoundingMode::None)?;
            piece = piece.div(&den, p, RoundingMode::None)?;
            if !sign_pos {
                piece.set_sign(Sign::Neg);
            }
            if first {
                sum = piece.clone()?;
                first = false;
            } else {
                sum = sum.add(&piece, p, RoundingMode::None)?;
            }
            if piece.is_zero() || (piece.exponent() as isize) + (p as isize) < 0 {
                break;
            }
            let np1 = n + 1;
            if sine {
                fact = fact
                    .mul(
                        &Self::from_word((2 * np1) as Word, p)?,
                        p,
                        RoundingMode::None,
                    )?
                    .mul(
                        &Self::from_word((2 * np1 + 1) as Word, p)?,
                        p,
                        RoundingMode::None,
                    )?;
            } else if np1 >= 1 {
                fact = fact
                    .mul(
                        &Self::from_word((2 * np1 - 1) as Word, p)?,
                        p,
                        RoundingMode::None,
                    )?
                    .mul(
                        &Self::from_word((2 * np1) as Word, p)?,
                        p,
                        RoundingMode::None,
                    )?;
            }
            pi2_pow = pi2_pow
                .mul(&pi2, p, RoundingMode::None)?
                .mul(&pi2, p, RoundingMode::None)?;
            pow_x = pow_x.mul(&x4, p, RoundingMode::None)?;
            sign_pos = !sign_pos;
        }
        Ok(sum)
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

        if self.is_zero() {
            let mut z =
                if n == 0 { Self::from_word(1, p)? } else { Self::new2(p, Sign::Pos, false)? };
            z.set_inexact(self.inexact());
            return Ok(z);
        }

        let mut p_inc = WORD_BIT_SIZE;
        let extra = extra_bits_for_degree(n as u32);
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc + extra;

        loop {
            let p_x = p_wrk + WORD_BIT_SIZE + extra;
            let mut ret = if n > 48 {
                self.bessel_j_miller(n, p_x)?
            } else {
                self.bessel_j_series(n, p_x)?
            };
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

    /// Miller backward recurrence. Normalize with \(J_0+2\sum_{k\ge 1}J_{2k}=1\).
    fn bessel_j_miller(&self, n: usize, p: usize) -> Result<Self, Error> {
        if self.is_zero() {
            return if n == 0 { Self::from_word(1, p) } else { Self::new2(p, Sign::Pos, false) };
        }
        let extra = (p / 2).saturating_add(8).max(16);
        let nstart = n.saturating_add(extra);
        let two = Self::from_word(2, p)?;
        let mut jp1 = Self::from_word(0, p)?;
        let mut jp = Self::from_word(1, p)?;
        let mut jn = Self::from_word(0, p)?;
        let mut even_sum = Self::from_word(0, p)?;
        let mut j0 = Self::from_word(0, p)?;
        for m in (1..=nstart).rev() {
            let two_m = Self::from_word((2 * m) as Word, p)?;
            let jm1 = two_m
                .div(self, p, RoundingMode::None)?
                .mul(&jp, p, RoundingMode::None)?
                .sub(&jp1, p, RoundingMode::None)?;
            if m == n {
                jn = jp.clone()?;
            }
            if m == 1 {
                j0 = jm1.clone()?;
            }
            if m % 2 == 0 {
                even_sum = even_sum.add(&jp, p, RoundingMode::None)?;
            }
            jp1 = jp;
            jp = jm1;
        }
        if n == 0 {
            jn = j0.clone()?;
        }
        let scale = j0.add(
            &two.mul(&even_sum, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
        if scale.is_zero() {
            return Err(Error::InvalidArgument);
        }
        jn.div(&scale, p, RoundingMode::None)
    }

    /// \(J_\nu(\mathrm{self})\) for real order `nu`. Negative `self` only for integer `nu`.
    pub fn bessel_j_nu(
        &self,
        nu: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        if let Some(n) = as_i32_exact(nu, p)? {
            let an = n.unsigned_abs() as usize;
            let j = self.bessel_j(an, p, rm, cc)?;
            if n >= 0 || n % 2 == 0 {
                Ok(j)
            } else {
                j.neg()
            }
        } else if !self.is_positive() {
            Err(Error::InvalidArgument)
        } else {
            self.bessel_series(nu, p, rm, cc, true)
        }
    }

    /// \(Y_\nu(\mathrm{self})\) for `self > 0`.
    pub fn bessel_y(
        &self,
        nu: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        if !self.is_positive() {
            return Err(Error::InvalidArgument);
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p
            .max(self.mantissa_max_bit_len())
            .max(nu.mantissa_max_bit_len())
            + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE;
            let mut ret = self.bessel_y_at(nu, p_x, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact() | nu.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn bessel_y_at(&self, nu: &Self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        if let Some(n) = as_i32_exact(nu, p)? {
            let an = n.unsigned_abs();
            let y = self.bessel_y_int(an, p, cc)?;
            if n >= 0 || n % 2 == 0 {
                Ok(y)
            } else {
                y.neg()
            }
        } else {
            let pi = cc.pi_num(p, RoundingMode::None)?;
            let nupi = nu.mul(&pi, p, RoundingMode::None)?;
            let s = nupi.sin(p, RoundingMode::None, cc)?;
            if s.is_zero() || (s.exponent() as isize) + (p as isize) < 0 {
                return Err(Error::InvalidArgument);
            }
            let jp = self.bessel_jn_series(nu, p, cc, true)?;
            let mut nneg = nu.clone()?;
            nneg.set_sign(if nu.is_positive() { Sign::Neg } else { Sign::Pos });
            let jm = self.bessel_jn_series(&nneg, p, cc, true)?;
            let c = nupi.cos(p, RoundingMode::None, cc)?;
            jp.mul(&c, p, RoundingMode::None)?
                .sub(&jm, p, RoundingMode::None)?
                .div(&s, p, RoundingMode::None)
        }
    }

    /// \(I_\nu(\mathrm{self})\). Integer \(\nu\) allows `self ≤ 0`; otherwise `self > 0`.
    pub fn bessel_i(
        &self,
        nu: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        if self.is_zero() {
            if nu.is_zero() {
                return Self::from_word(1, p);
            }
            if nu.is_positive() {
                return Self::new2(p, Sign::Pos, self.inexact() | nu.inexact());
            }
            return Err(Error::InvalidArgument);
        }
        if let Some(n) = as_i32_exact(nu, p)? {
            let an = Self::from_word(n.unsigned_abs() as Word, p)?;
            let mut ax = self.clone()?;
            ax.set_sign(Sign::Pos);
            let i = ax.bessel_series(&an, p, rm, cc, false)?;
            if self.is_negative() && n.unsigned_abs() % 2 == 1 {
                i.neg()
            } else {
                Ok(i)
            }
        } else if !self.is_positive() {
            Err(Error::InvalidArgument)
        } else {
            self.bessel_series(nu, p, rm, cc, false)
        }
    }

    /// \(K_\nu(\mathrm{self})\) for `self > 0`. \(K_{-ν}=K_ν\).
    pub fn bessel_k(
        &self,
        nu: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        if !self.is_positive() {
            return Err(Error::InvalidArgument);
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p
            .max(self.mantissa_max_bit_len())
            .max(nu.mantissa_max_bit_len())
            + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE * 2;
            let mut ret = self.bessel_k_at(nu, p_x, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact() | nu.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn bessel_k_at(&self, nu: &Self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let mut nu = nu.clone()?;
        nu.set_sign(Sign::Pos);
        if let Some(n) = half_integer_n(&nu, p)? {
            return self.k_half_integer(n, p, cc);
        }
        if self.exponent() >= 5 {
            return self.k_asymptotic(&nu, p, cc);
        }
        if let Some(n) = as_i32_exact(&nu, p)? {
            if n < 0 {
                return Err(Error::InvalidArgument);
            }
            return self.k_integer(n as u32, p, cc);
        }
        self.k_real(&nu, p, cc)
    }

    fn bessel_series(
        &self,
        nu: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
        alternating: bool,
    ) -> Result<Self, Error> {
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p
            .max(self.mantissa_max_bit_len())
            .max(nu.mantissa_max_bit_len())
            + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE;
            let mut ret = self.bessel_jn_series(nu, p_x, cc, alternating)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact() | nu.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn bessel_jn_series(
        &self,
        nu: &Self,
        p: usize,
        cc: &mut Consts,
        alternating: bool,
    ) -> Result<Self, Error> {
        let two = Self::from_word(2, p)?;
        let half = self.div(&two, p, RoundingMode::None)?;
        let one = Self::from_word(1, p)?;
        let g0 = nu.add(&one, p, RoundingMode::None)?.gamma_at(p, cc)?;
        let mut term = half
            .pow(nu, p, RoundingMode::None, cc)?
            .div(&g0, p, RoundingMode::None)?;
        let mut sum = term.clone()?;
        let hh = half.mul(&half, p, RoundingMode::None)?;
        for k in 1..=series_n_max(p, self.exponent()) {
            let kk = Self::from_word(k as Word, p)?;
            let den = kk
                .add(nu, p, RoundingMode::None)?
                .mul(&kk, p, RoundingMode::None)?;
            term = term
                .mul(&hh, p, RoundingMode::None)?
                .div(&den, p, RoundingMode::None)?;
            if alternating {
                term = term.neg()?;
            }
            sum = sum.add(&term, p, RoundingMode::None)?;
            if term.is_zero() || (term.exponent() as isize) + (p as isize) < 0 {
                break;
            }
        }
        Ok(sum)
    }

    fn bessel_y_int(&self, n: u32, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        if n == 0 {
            return self.bessel_y0(p, cc);
        }
        if n == 1 {
            return self.bessel_y1(p, cc);
        }
        let mut ym2 = self.bessel_y0(p, cc)?;
        let mut ym1 = self.bessel_y1(p, cc)?;
        for m in 1..n {
            let two_m = Self::from_word((2 * m) as Word, p)?;
            let ym = two_m
                .div(self, p, RoundingMode::None)?
                .mul(&ym1, p, RoundingMode::None)?
                .sub(&ym2, p, RoundingMode::None)?;
            ym2 = ym1;
            ym1 = ym;
        }
        Ok(ym1)
    }

    fn bessel_y0(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let two = Self::from_word(2, p)?;
        let two_pi = two.div(&pi, p, RoundingMode::None)?;
        let half = self.div(&two, p, RoundingMode::None)?;
        let j0 = self.bessel_j(0, p, RoundingMode::None, cc)?;
        let g = cc.euler_gamma_num(p, RoundingMode::None)?;
        let prefix = g.add(&half.ln(p, RoundingMode::None, cc)?, p, RoundingMode::None)?;
        let z = half.mul(&half, p, RoundingMode::None)?;
        let mut fact = Self::from_word(1, p)?;
        let mut zk = Self::from_word(1, p)?;
        let mut sum = Self::from_word(0, p)?;
        for m in 1..=series_n_max(p, self.exponent()) {
            fact = fact.mul(&Self::from_word(m as Word, p)?, p, RoundingMode::None)?;
            zk = zk.mul(&z, p, RoundingMode::None)?;
            let h = harmonic_u(m, p)?;
            let mut term = h.div(
                &fact.mul(&fact, p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )?;
            term = term.mul(&zk, p, RoundingMode::None)?;
            if m % 2 == 0 {
                term = term.neg()?;
            }
            sum = sum.add(&term, p, RoundingMode::None)?;
            if term.is_zero() || (term.exponent() as isize) + (p as isize) < 0 {
                break;
            }
        }
        two_pi.mul(
            &prefix
                .mul(&j0, p, RoundingMode::None)?
                .add(&sum, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )
    }

    fn bessel_y1(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let two = Self::from_word(2, p)?;
        let two_pi = two.div(&pi, p, RoundingMode::None)?;
        let half = self.div(&two, p, RoundingMode::None)?;
        let j1 = self.bessel_j(1, p, RoundingMode::None, cc)?;
        let g = cc.euler_gamma_num(p, RoundingMode::None)?;
        let prefix = g.add(&half.ln(p, RoundingMode::None, cc)?, p, RoundingMode::None)?;
        let z = half.mul(&half, p, RoundingMode::None)?;
        let mut kfact = Self::from_word(1, p)?;
        let mut kp1fact = Self::from_word(1, p)?;
        let mut zk = Self::from_word(1, p)?;
        let mut sum = Self::from_word(0, p)?;
        for k in 0..=series_n_max(p, self.exponent()) {
            let hk = harmonic_u(k, p)?;
            let hkp1 = hk.add(
                &Self::from_word((k + 1) as Word, p)?.reciprocal(p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )?;
            let mut term = hk
                .add(&hkp1, p, RoundingMode::None)?
                .div(
                    &kfact.mul(&kp1fact, p, RoundingMode::None)?,
                    p,
                    RoundingMode::None,
                )?
                .mul(&zk, p, RoundingMode::None)?;
            if k % 2 == 1 {
                term = term.neg()?;
            }
            sum = sum.add(&term, p, RoundingMode::None)?;
            if (term.is_zero() || (term.exponent() as isize) + (p as isize) < 0) && k > 0 {
                break;
            }
            let kp = k + 1;
            kfact = kfact.mul(&Self::from_word(kp as Word, p)?, p, RoundingMode::None)?;
            kp1fact = kp1fact.mul(
                &Self::from_word((kp + 1) as Word, p)?,
                p,
                RoundingMode::None,
            )?;
            zk = zk.mul(&z, p, RoundingMode::None)?;
        }
        let a = two_pi.mul(
            &prefix.mul(&j1, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
        let b = two.div(&pi.mul(self, p, RoundingMode::None)?, p, RoundingMode::None)?;
        let c = self
            .div(&two.mul(&pi, p, RoundingMode::None)?, p, RoundingMode::None)?
            .mul(&sum, p, RoundingMode::None)?;
        a.sub(&b, p, RoundingMode::None)?
            .sub(&c, p, RoundingMode::None)
    }

    fn k_half_integer(&self, n: u32, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let two = Self::from_word(2, p)?;
        let two_x = two.mul(self, p, RoundingMode::None)?;
        let mut term = Self::from_word(1, p)?;
        let mut sum = Self::from_word(1, p)?;
        for k in 0..n {
            let kk = Self::from_word((k + 1) as Word, p)?;
            let num = Self::from_word((n + k + 1) as Word, p)?.mul(
                &Self::from_word((n - k) as Word, p)?,
                p,
                RoundingMode::None,
            )?;
            term = term.mul(&num, p, RoundingMode::None)?.div(
                &kk.mul(&two_x, p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )?;
            sum = sum.add(&term, p, RoundingMode::None)?;
        }
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let mut nx = self.clone()?;
        nx.set_sign(Sign::Neg);
        let pref = pi
            .div(&two_x, p, RoundingMode::None)?
            .sqrt(p, RoundingMode::None)?
            .mul(&nx.exp(p, RoundingMode::None, cc)?, p, RoundingMode::None)?;
        pref.mul(&sum, p, RoundingMode::None)
    }

    fn k_integer(&self, n: u32, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let k0 = self.k0_series(p, cc)?;
        if n == 0 {
            return Ok(k0);
        }
        let k1 = self.k1_from_wronskian(&k0, p, cc)?;
        if n == 1 {
            return Ok(k1);
        }
        let mut km2 = k0;
        let mut km1 = k1;
        for m in 1..n {
            let two_m = Self::from_word((2 * m) as Word, p)?;
            let km = two_m
                .div(self, p, RoundingMode::None)?
                .mul(&km1, p, RoundingMode::None)?
                .add(&km2, p, RoundingMode::None)?;
            km2 = km1;
            km1 = km;
        }
        Ok(km1)
    }

    fn k0_series(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let zero = Self::from_word(0, p)?;
        let one = Self::from_word(1, p)?;
        let i0 = self.bessel_i(&zero, p, RoundingMode::None, cc)?;
        let two = Self::from_word(2, p)?;
        let half = self.div(&two, p, RoundingMode::None)?;
        let z = half.mul(&half, p, RoundingMode::None)?;
        let mut fact = Self::from_word(1, p)?;
        let mut zk = Self::from_word(1, p)?;
        let mut sum = one.digamma_at(p, cc)?;
        for k in 1..=series_n_max(p, self.exponent()) {
            fact = fact.mul(&Self::from_word(k as Word, p)?, p, RoundingMode::None)?;
            zk = zk.mul(&z, p, RoundingMode::None)?;
            let psi = Self::from_word((k + 1) as Word, p)?.digamma_at(p, cc)?;
            let term = psi
                .div(
                    &fact.mul(&fact, p, RoundingMode::None)?,
                    p,
                    RoundingMode::None,
                )?
                .mul(&zk, p, RoundingMode::None)?;
            sum = sum.add(&term, p, RoundingMode::None)?;
            if term.is_zero() || (term.exponent() as isize) + (p as isize) < 0 {
                break;
            }
        }
        half.ln(p, RoundingMode::None, cc)?
            .neg()?
            .mul(&i0, p, RoundingMode::None)?
            .add(&sum, p, RoundingMode::None)
    }

    fn k1_from_wronskian(&self, k0: &Self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let zero = Self::from_word(0, p)?;
        let one = Self::from_word(1, p)?;
        let i0 = self.bessel_i(&zero, p, RoundingMode::None, cc)?;
        let i1 = self.bessel_i(&one, p, RoundingMode::None, cc)?;
        if i0.is_zero() {
            return Err(Error::InvalidArgument);
        }
        one.div(self, p, RoundingMode::None)?
            .sub(&i1.mul(k0, p, RoundingMode::None)?, p, RoundingMode::None)?
            .div(&i0, p, RoundingMode::None)
    }

    fn k_real(&self, nu: &Self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let n_floor = floor_nonneg_u32(nu, p)?;
        let f = nu.sub(&Self::from_word(n_floor as Word, p)?, p, RoundingMode::None)?;
        if f.is_zero() {
            return self.k_integer(n_floor, p, cc);
        }
        let kf = self.k_connection(&f, p, cc)?;
        let one = Self::from_word(1, p)?;
        let k1mf = self.k_connection(&one.sub(&f, p, RoundingMode::None)?, p, cc)?;
        let mut km1 = k1mf;
        let mut k0 = kf;
        for m in 0..n_floor {
            let nu_m = f.add(&Self::from_word(m as Word, p)?, p, RoundingMode::None)?;
            let kp = two_times(&nu_m, p)?
                .div(self, p, RoundingMode::None)?
                .mul(&k0, p, RoundingMode::None)?
                .add(&km1, p, RoundingMode::None)?;
            km1 = k0;
            k0 = kp;
        }
        Ok(k0)
    }

    fn k_connection(&self, mu: &Self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let s = mu
            .mul(&pi, p, RoundingMode::None)?
            .sin(p, RoundingMode::None, cc)?;
        if s.is_zero() || (s.exponent() as isize) + (p as isize) < 0 {
            return Err(Error::InvalidArgument);
        }
        let ip = self.bessel_i(mu, p, RoundingMode::None, cc)?;
        let mut mneg = mu.clone()?;
        mneg.set_sign(Sign::Neg);
        let im = self.bessel_i(&mneg, p, RoundingMode::None, cc)?;
        let two = Self::from_word(2, p)?;
        pi.div(&two, p, RoundingMode::None)?
            .mul(&im.sub(&ip, p, RoundingMode::None)?, p, RoundingMode::None)?
            .div(&s, p, RoundingMode::None)
    }

    fn k_asymptotic(&self, nu: &Self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let four = Self::from_word(4, p)?;
        let eight = Self::from_word(8, p)?;
        let four_nu2 = four.mul(&nu.mul(nu, p, RoundingMode::None)?, p, RoundingMode::None)?;
        let eight_x = eight.mul(self, p, RoundingMode::None)?;
        let mut term = Self::from_word(1, p)?;
        let mut sum = Self::from_word(1, p)?;
        let mut prev_e = i32::MAX;
        for k in 1..=40 {
            let odd = Self::from_word((2 * k - 1) as Word, p)?;
            let factor = four_nu2
                .sub(
                    &odd.mul(&odd, p, RoundingMode::None)?,
                    p,
                    RoundingMode::None,
                )?
                .div(
                    &Self::from_word(k as Word, p)?.mul(&eight_x, p, RoundingMode::None)?,
                    p,
                    RoundingMode::None,
                )?;
            term = term.mul(&factor, p, RoundingMode::None)?;
            let e = term.exponent();
            if k > 3 && e > prev_e {
                break;
            }
            sum = sum.add(&term, p, RoundingMode::None)?;
            if term.is_zero() || (e as isize) + (p as isize) < 0 {
                break;
            }
            prev_e = e;
        }
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let two = Self::from_word(2, p)?;
        let mut nx = self.clone()?;
        nx.set_sign(Sign::Neg);
        let pref = pi
            .div(
                &two.mul(self, p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )?
            .sqrt(p, RoundingMode::None)?
            .mul(&nx.exp(p, RoundingMode::None, cc)?, p, RoundingMode::None)?;
        pref.mul(&sum, p, RoundingMode::None)
    }

    /// Complete elliptic \(K(m)\). Parameter \(m=k^2\). \(m=1\) is \(+\infty\);
    /// \(m>1\) uses \(K(m)=m^{-1/2}K(1/m)\).
    pub fn elliptic_k(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        let one = Self::from_word(1, p)?;
        if self.cmp(&one) == 0 {
            return Err(Error::ExponentOverflow(Sign::Pos));
        }
        if self.cmp(&one) > 0 {
            if self.is_negative() {
                return Err(Error::InvalidArgument);
            }
            let mut p_inc = WORD_BIT_SIZE;
            let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc;
            loop {
                let p_x = p_wrk + WORD_BIT_SIZE * 2;
                let inv = one.div(self, p_x, RoundingMode::None)?;
                let k = inv.elliptic_k_at(p_x, cc)?;
                let mut ret = k.div(
                    &self.sqrt(p_x, RoundingMode::None)?,
                    p_x,
                    RoundingMode::None,
                )?;
                if ret.try_set_precision(p, rm, p_wrk)? {
                    ret.set_inexact(ret.inexact() | self.inexact());
                    return Ok(ret);
                }
                bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
            }
        }
        if self.is_zero() {
            let mut hp = cc.pi_num(p, rm)?;
            hp.div_by_2(rm);
            return Ok(hp);
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE * 2;
            let mut ret = self.elliptic_k_at(p_x, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn elliptic_k_at(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let zero = Self::from_word(0, p)?;
        let one = Self::from_word(1, p)?;
        let om = one.sub(self, p, RoundingMode::None)?;
        carlson_rf(&zero, &om, &one, p, cc)
    }

    /// Complete elliptic \(E(m)\), \(m\le 1\). \(E(1)=1\).
    pub fn elliptic_e_complete(
        &self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        let one = Self::from_word(1, p)?;
        if self.cmp(&one) > 0 {
            return Err(Error::InvalidArgument);
        }
        if self.is_zero() {
            let mut hp = cc.pi_num(p, rm)?;
            hp.div_by_2(rm);
            return Ok(hp);
        }
        if self.cmp(&one) == 0 {
            return Self::from_word(1, p);
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE * 2;
            let mut ret = self.elliptic_e_complete_at(p_x, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn elliptic_e_complete_at(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let zero = Self::from_word(0, p)?;
        let one = Self::from_word(1, p)?;
        let three = Self::from_word(3, p)?;
        let om = one.sub(self, p, RoundingMode::None)?;
        let rf = carlson_rf(&zero, &om, &one, p, cc)?;
        let rd = carlson_rd(&zero, &om, &one, p, cc)?;
        rf.sub(
            &self
                .div(&three, p, RoundingMode::None)?
                .mul(&rd, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )
    }

    /// Incomplete \(F(x|m)\), \(|x|\le 1\). \(x=\sin\varphi\), \(m=k^2\).
    pub fn elliptic_f(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        let one = Self::from_word(1, p)?;
        if self.abs()?.cmp(&one) > 0 {
            return Err(Error::InvalidArgument);
        }
        if self.is_zero() {
            return Self::new2(p, Sign::Pos, self.inexact() | m.inexact());
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p
            .max(self.mantissa_max_bit_len())
            .max(m.mantissa_max_bit_len())
            + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE * 2;
            let mut ret = self.elliptic_f_at(m, p_x, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact() | m.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn elliptic_f_at(&self, m: &Self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let one = Self::from_word(1, p)?;
        let xa = self.abs()?;
        let sign_neg = self.is_negative();
        if m.is_zero() {
            let mut r = xa.asin(p, RoundingMode::None, cc)?;
            if sign_neg {
                r.set_sign(Sign::Neg);
            }
            return Ok(r);
        }
        if m.cmp(&one) == 0 {
            if xa.cmp(&one) == 0 {
                return Err(Error::InvalidArgument);
            }
            let mut r = xa.atanh(p, RoundingMode::None, cc)?;
            if sign_neg {
                r.set_sign(Sign::Neg);
            }
            return Ok(r);
        }
        let x2 = xa.mul(&xa, p, RoundingMode::None)?;
        let a = one.sub(&x2, p, RoundingMode::None)?;
        let b = one.sub(&m.mul(&x2, p, RoundingMode::None)?, p, RoundingMode::None)?;
        if a.is_negative() || b.is_negative() {
            return Err(Error::InvalidArgument);
        }
        let rf = carlson_rf(&a, &b, &one, p, cc)?;
        let mut r = xa.mul(&rf, p, RoundingMode::None)?;
        if sign_neg {
            r.set_sign(Sign::Neg);
        }
        Ok(r)
    }

    /// Incomplete \(E(x|m)\), \(|x|\le 1\).
    pub fn elliptic_e(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        let one = Self::from_word(1, p)?;
        if self.abs()?.cmp(&one) > 0 {
            return Err(Error::InvalidArgument);
        }
        if self.is_zero() {
            return Self::new2(p, Sign::Pos, self.inexact() | m.inexact());
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p
            .max(self.mantissa_max_bit_len())
            .max(m.mantissa_max_bit_len())
            + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE * 2;
            let mut ret = self.elliptic_e_at(m, p_x, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact() | m.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn elliptic_e_at(&self, m: &Self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let one = Self::from_word(1, p)?;
        let xa = self.abs()?;
        let sign_neg = self.is_negative();
        if m.is_zero() {
            let mut r = xa.asin(p, RoundingMode::None, cc)?;
            if sign_neg {
                r.set_sign(Sign::Neg);
            }
            return Ok(r);
        }
        if m.cmp(&one) == 0 {
            let mut r = xa;
            if sign_neg {
                r.set_sign(Sign::Neg);
            }
            return Ok(r);
        }
        let x2 = xa.mul(&xa, p, RoundingMode::None)?;
        let a = one.sub(&x2, p, RoundingMode::None)?;
        let b = one.sub(&m.mul(&x2, p, RoundingMode::None)?, p, RoundingMode::None)?;
        if a.is_negative() || b.is_negative() {
            return Err(Error::InvalidArgument);
        }
        let three = Self::from_word(3, p)?;
        let rf = carlson_rf(&a, &b, &one, p, cc)?;
        let rd = carlson_rd(&a, &b, &one, p, cc)?;
        let t = xa.mul(&rf, p, RoundingMode::None)?.sub(
            &m.mul(&xa, p, RoundingMode::None)?
                .mul(&x2, p, RoundingMode::None)?
                .div(&three, p, RoundingMode::None)?
                .mul(&rd, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
        if sign_neg {
            t.neg()
        } else {
            Ok(t)
        }
    }

    /// Complete \(\Pi(n,m)\), \(m<1\), \(n<1\). `self` is \(n\).
    pub fn elliptic_pi_complete(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        let one = Self::from_word(1, p)?;
        if m.cmp(&one) >= 0 || self.cmp(&one) >= 0 {
            return Err(Error::InvalidArgument);
        }
        if self.is_zero() {
            return m.elliptic_k(p, rm, cc);
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p
            .max(self.mantissa_max_bit_len())
            .max(m.mantissa_max_bit_len())
            + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE * 2;
            let mut ret = self.elliptic_pi_complete_at(m, p_x, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact() | m.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn elliptic_pi_complete_at(&self, m: &Self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let zero = Self::from_word(0, p)?;
        let one = Self::from_word(1, p)?;
        let three = Self::from_word(3, p)?;
        let om = one.sub(m, p, RoundingMode::None)?;
        let on = one.sub(self, p, RoundingMode::None)?;
        let rf = carlson_rf(&zero, &om, &one, p, cc)?;
        let rj = carlson_rj(&zero, &om, &one, &on, p, cc)?;
        rf.add(
            &self
                .div(&three, p, RoundingMode::None)?
                .mul(&rj, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )
    }

    /// Incomplete \(\Pi(n;x|m)\). `self` is \(n\).
    pub fn elliptic_pi(
        &self,
        x: &Self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        let one = Self::from_word(1, p)?;
        if x.abs()?.cmp(&one) > 0 {
            return Err(Error::InvalidArgument);
        }
        if x.is_zero() {
            return Self::new2(p, Sign::Pos, self.inexact() | x.inexact() | m.inexact());
        }
        if self.is_zero() {
            return x.elliptic_f(m, p, rm, cc);
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p
            .max(self.mantissa_max_bit_len())
            .max(x.mantissa_max_bit_len())
            .max(m.mantissa_max_bit_len())
            + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE * 2;
            let mut ret = self.elliptic_pi_at(x, m, p_x, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact() | x.inexact() | m.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn elliptic_pi_at(&self, x: &Self, m: &Self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let one = Self::from_word(1, p)?;
        let xa = x.abs()?;
        let sign_neg = x.is_negative();
        let x2 = xa.mul(&xa, p, RoundingMode::None)?;
        let a = one.sub(&x2, p, RoundingMode::None)?;
        let b = one.sub(&m.mul(&x2, p, RoundingMode::None)?, p, RoundingMode::None)?;
        let pv = one.sub(
            &self.mul(&x2, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
        if a.is_negative() || b.is_negative() || !pv.is_positive() {
            return Err(Error::InvalidArgument);
        }
        if m.is_zero() && self.cmp(&one) < 0 {
            let root = one
                .sub(self, p, RoundingMode::None)?
                .sqrt(p, RoundingMode::None)?;
            let inner = root
                .mul(&xa, p, RoundingMode::None)?
                .div(&a.sqrt(p, RoundingMode::None)?, p, RoundingMode::None)?
                .atan(p, RoundingMode::None, cc)?;
            let mut r = inner.div(&root, p, RoundingMode::None)?;
            if sign_neg {
                r.set_sign(Sign::Neg);
            }
            return Ok(r);
        }
        let three = Self::from_word(3, p)?;
        let rf = carlson_rf(&a, &b, &one, p, cc)?;
        let rj = carlson_rj(&a, &b, &one, &pv, p, cc)?;
        let t = xa.mul(&rf, p, RoundingMode::None)?.add(
            &self
                .mul(&xa, p, RoundingMode::None)?
                .mul(&x2, p, RoundingMode::None)?
                .div(&three, p, RoundingMode::None)?
                .mul(&rj, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
        if sign_neg {
            t.neg()
        } else {
            Ok(t)
        }
    }

    /// Legendre \(P_n(\mathrm{self})\) for integer \(n\). Recurrence at \(p+O(n)\) bits.
    pub fn legendre_p(&self, n: u32, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        let one = Self::from_word(1, p)?;
        if self.cmp(&one) == 0 {
            let mut r = one;
            r.set_inexact(self.inexact());
            return Ok(r);
        }
        if self.cmp(&one.neg()?) == 0 {
            let mut r = if n % 2 == 0 { one } else { one.neg()? };
            r.set_inexact(self.inexact());
            return Ok(r);
        }
        let extra = extra_bits_for_degree(n);
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc + extra;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE + extra;
            let mut ret = self.legendre_p_at(n, p_x)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn legendre_p_at(&self, n: u32, p: usize) -> Result<Self, Error> {
        if n == 0 {
            return Self::from_word(1, p);
        }
        if n == 1 {
            return self.clone();
        }
        let mut pm2 = Self::from_word(1, p)?;
        let mut pm1 = self.clone()?;
        let mut k = 1u32;
        while k < n {
            let kk = Self::from_word(k as Word, p)?;
            let kp1 = Self::from_word((k + 1) as Word, p)?;
            let two_k1 = Self::from_word((2 * k + 1) as Word, p)?;
            let pnext = two_k1
                .mul(self, p, RoundingMode::None)?
                .mul(&pm1, p, RoundingMode::None)?
                .sub(&kk.mul(&pm2, p, RoundingMode::None)?, p, RoundingMode::None)?
                .div(&kp1, p, RoundingMode::None)?;
            pm2 = pm1;
            pm1 = pnext;
            k += 1;
        }
        Ok(pm1)
    }

    /// Associated \(P_n^m(\mathrm{self})\) with Condon–Shortley phase \((-1)^m\).
    /// \(\lvert m\rvert\le n\); \(\lvert x\rvert\le 1\) when \(m\ne 0\).
    pub fn assoc_legendre_p(
        &self,
        n: u32,
        m: i32,
        p: usize,
        rm: RoundingMode,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        let am = m.unsigned_abs();
        if am > n {
            return Err(Error::InvalidArgument);
        }
        if m != 0 && self.abs()?.cmp(&Self::from_word(1, p)?) > 0 {
            return Err(Error::InvalidArgument);
        }
        let extra = extra_bits_for_degree(n);
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc + extra;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE + extra;
            let mut ret = self.assoc_legendre_p_at(n, m, p_x)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn assoc_legendre_p_at(&self, n: u32, m: i32, p: usize) -> Result<Self, Error> {
        let am = m.unsigned_abs();
        if m == 0 {
            return self.legendre_p_at(n, p);
        }
        let ppos = self.assoc_pos(n, am, p)?;
        if m > 0 {
            return Ok(ppos);
        }
        let ratio = fact_ratio_down(n - am, n + am, p)?;
        let mut out = ratio.mul(&ppos, p, RoundingMode::None)?;
        if am % 2 == 1 {
            out.set_sign(if out.is_positive() { Sign::Neg } else { Sign::Pos });
        }
        Ok(out)
    }

    fn assoc_pos(&self, n: u32, m: u32, p: usize) -> Result<Self, Error> {
        let one = Self::from_word(1, p)?;
        let one_x2 = one.sub(
            &self.mul(self, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
        if one_x2.is_negative() {
            return Err(Error::InvalidArgument);
        }
        let mut pmm = Self::from_word(1, p)?;
        if m > 0 {
            let root = one_x2.sqrt(p, RoundingMode::None)?;
            let mut odd = Self::from_word(1, p)?;
            let two = Self::from_word(2, p)?;
            for _j in 1..=m {
                pmm = pmm.neg()?.mul(&odd, p, RoundingMode::None)?.mul(
                    &root,
                    p,
                    RoundingMode::None,
                )?;
                odd = odd.add(&two, p, RoundingMode::None)?;
            }
        }
        if n == m {
            return Ok(pmm);
        }
        let mut pm1 = self
            .mul(
                &Self::from_word((2 * m + 1) as Word, p)?,
                p,
                RoundingMode::None,
            )?
            .mul(&pmm, p, RoundingMode::None)?;
        if n == m + 1 {
            return Ok(pm1);
        }
        let mut pm2 = pmm;
        let mut k = m + 1;
        while k < n {
            let num1 = Self::from_word((2 * k + 1) as Word, p)?
                .mul(self, p, RoundingMode::None)?
                .mul(&pm1, p, RoundingMode::None)?;
            let num2 = Self::from_word((k + m) as Word, p)?.mul(&pm2, p, RoundingMode::None)?;
            let den = Self::from_word((k - m + 1) as Word, p)?;
            let pn = num1
                .sub(&num2, p, RoundingMode::None)?
                .div(&den, p, RoundingMode::None)?;
            pm2 = pm1;
            pm1 = pn;
            k += 1;
        }
        Ok(pm1)
    }

    /// Gaussian \({}_2F_1(a=\mathrm{self},b;c;z)\).
    /// Series for \(\lvert z\rvert<1\) or terminating \(a\) or \(b\); Gauss at \(z=1\)
    /// when \(c-a-b>0\); Pfaff on \(z\in(1/2,1)\); real continuation for \(z\le -1\).
    /// Non-terminating \(z>1\) (typically not real) is `InvalidArgument`.
    pub fn hypergeom_2f1(
        &self,
        b: &Self,
        c: &Self,
        z: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        if z.is_zero() || self.is_zero() || b.is_zero() {
            let mut one = Self::from_word(1, p)?;
            one.set_inexact(self.inexact() | b.inexact() | c.inexact() | z.inexact());
            return Ok(one);
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p
            .max(self.mantissa_max_bit_len())
            .max(b.mantissa_max_bit_len())
            .max(c.mantissa_max_bit_len())
            .max(z.mantissa_max_bit_len())
            + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE * 2;
            let mut ret = self.hypergeom_2f1_at(b, c, z, p_x, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(
                    ret.inexact() | self.inexact() | b.inexact() | c.inexact() | z.inexact(),
                );
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn hypergeom_2f1_at(
        &self,
        b: &Self,
        c: &Self,
        z: &Self,
        p: usize,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        if pole_c_before_term(c, self, b, p)? {
            return Err(Error::InvalidArgument);
        }
        if terminating_neg_int(self, p)? || terminating_neg_int(b, p)? {
            return hypergeom_series(self, b, c, z, p);
        }
        let one = Self::from_word(1, p)?;
        if z.cmp(&one) == 0 {
            return gauss_z_one(self, b, c, p, cc);
        }
        if z.cmp(&one) > 0 {
            return Err(Error::InvalidArgument);
        }
        if z.cmp(&one.neg()?) <= 0 {
            return hypergeom_z_le_neg_one(self, b, c, z, p, cc);
        }
        let half = one_half(p)?;
        if z.cmp(&half) > 0 && z.cmp(&one) < 0 {
            let zm1 = z.sub(&one, p, RoundingMode::None)?;
            let w = z.div(&zm1, p, RoundingMode::None)?;
            let mut na = self.clone()?;
            na.set_sign(if self.is_positive() { Sign::Neg } else { Sign::Pos });
            let pref = one
                .sub(z, p, RoundingMode::None)?
                .pow(&na, p, RoundingMode::None, cc)?;
            let cb = c.sub(b, p, RoundingMode::None)?;
            let inner = hypergeom_series(self, &cb, c, &w, p)?;
            return pref.mul(&inner, p, RoundingMode::None);
        }
        hypergeom_series(self, b, c, z, p)
    }

    /// Regularized incomplete beta \(I_x(a=\mathrm{self},b)\) for \(a>0\), \(b>0\), \(x\in[0,1]\).
    pub fn betainc(
        &self,
        b: &Self,
        x: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        if !self.is_positive() || !b.is_positive() {
            return Err(Error::InvalidArgument);
        }
        if x.is_negative() || x.cmp(&Self::from_word(1, p)?) > 0 {
            return Err(Error::InvalidArgument);
        }
        if x.is_zero() {
            return Self::new2(p, Sign::Pos, self.inexact() | b.inexact() | x.inexact());
        }
        if x.cmp(&Self::from_word(1, p)?) == 0 {
            let mut one = Self::from_word(1, p)?;
            one.set_inexact(self.inexact() | b.inexact() | x.inexact());
            return Ok(one);
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p
            .max(self.mantissa_max_bit_len())
            .max(b.mantissa_max_bit_len())
            .max(x.mantissa_max_bit_len())
            + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE * 2;
            let mut ret = self.betainc_at(b, x, p_x, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact() | b.inexact() | x.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn betainc_at(&self, b: &Self, x: &Self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let one = Self::from_word(1, p)?;
        let two = Self::from_word(2, p)?;
        let thresh = self.add(&one, p, RoundingMode::None)?.div(
            &self
                .add(b, p, RoundingMode::None)?
                .add(&two, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
        if x.cmp(&thresh) > 0 {
            let ox = one.sub(x, p, RoundingMode::None)?;
            let t = b.betainc_direct(self, &ox, p, cc)?;
            return one.sub(&t, p, RoundingMode::None);
        }
        self.betainc_direct(b, x, p, cc)
    }

    fn betainc_direct(&self, b: &Self, x: &Self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        // I_x(a,b) = x^a / (a B(a,b)) · ₂F₁(a, 1−b; a+1; x)
        let one = Self::from_word(1, p)?;
        let ga = self.gamma_at(p, cc)?;
        let gb = b.gamma_at(p, cc)?;
        let gab = self.add(b, p, RoundingMode::None)?.gamma_at(p, cc)?;
        let beta = ga
            .mul(&gb, p, RoundingMode::None)?
            .div(&gab, p, RoundingMode::None)?;
        let mut nb = b.clone()?;
        nb.set_sign(if b.is_positive() { Sign::Neg } else { Sign::Pos });
        let one_b = one.add(&nb, p, RoundingMode::None)?;
        let ap1 = self.add(&one, p, RoundingMode::None)?;
        let f = self.hypergeom_2f1_at(&one_b, &ap1, x, p, cc)?;
        let xa = x.pow(self, p, RoundingMode::None, cc)?;
        xa.div(self, p, RoundingMode::None)?
            .div(&beta, p, RoundingMode::None)?
            .mul(&f, p, RoundingMode::None)
    }
}

const CARLSON_DUPE_MAX: u32 = 128;
const HYPERGEOM_TERM_MAX: u32 = 10_000;

fn fact_ratio_down(a: u32, b: u32, p: usize) -> Result<ExactNumNumber, Error> {
    let mut acc = ExactNumNumber::from_word(1, p)?;
    let mut k = a + 1;
    while k <= b {
        acc = acc.div(
            &ExactNumNumber::from_word(k as Word, p)?,
            p,
            RoundingMode::None,
        )?;
        k += 1;
    }
    Ok(acc)
}

fn terminating_neg_int(v: &ExactNumNumber, p: usize) -> Result<bool, Error> {
    Ok(as_i32_exact(v, p)?.is_some_and(|n| n <= 0))
}

fn pole_c_before_term(
    c: &ExactNumNumber,
    a: &ExactNumNumber,
    b: &ExactNumNumber,
    p: usize,
) -> Result<bool, Error> {
    let Some(cn) = as_i32_exact(c, p)? else {
        return Ok(false);
    };
    if cn > 0 {
        return Ok(false);
    }
    let stop_a = as_i32_exact(a, p)?.filter(|&n| n <= 0);
    let stop_b = as_i32_exact(b, p)?.filter(|&n| n <= 0);
    let pole_k = -cn;
    Ok(match (stop_a, stop_b) {
        (Some(sa), _) if -sa < pole_k => false,
        (_, Some(sb)) if -sb < pole_k => false,
        _ => true,
    })
}

fn hypergeom_series(
    a: &ExactNumNumber,
    b: &ExactNumNumber,
    c: &ExactNumNumber,
    z: &ExactNumNumber,
    p: usize,
) -> Result<ExactNumNumber, Error> {
    let one = ExactNumNumber::from_word(1, p)?;
    let mut term = ExactNumNumber::from_word(1, p)?;
    let mut sum = ExactNumNumber::from_word(1, p)?;
    let tiny = ExactNumNumber::from_word(1, p)?.ldexp(-((p as i32) - 8), p, RoundingMode::None)?;
    let n_max = series_n_max(p, z.exponent()).min(HYPERGEOM_TERM_MAX as usize);
    for n in 0..n_max {
        if n > 0 && (term.is_zero() || (term.exponent() as isize) + (p as isize) < 0) {
            break;
        }
        let nn = ExactNumNumber::from_word(n as Word, p)?;
        let den = c.add(&nn, p, RoundingMode::None)?.mul(
            &one.add(&nn, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
        if den.abs()?.cmp(&tiny) < 0 {
            return Err(Error::InvalidArgument);
        }
        term = term
            .mul(&a.add(&nn, p, RoundingMode::None)?, p, RoundingMode::None)?
            .mul(&b.add(&nn, p, RoundingMode::None)?, p, RoundingMode::None)?
            .div(&den, p, RoundingMode::None)?
            .mul(z, p, RoundingMode::None)?;
        sum = sum.add(&term, p, RoundingMode::None)?;
        if terminating_after(a, n + 1, p)? || terminating_after(b, n + 1, p)? {
            sum.set_inexact(false);
            return Ok(sum);
        }
    }
    Ok(sum)
}

fn terminating_after(v: &ExactNumNumber, next_n: usize, p: usize) -> Result<bool, Error> {
    Ok(as_i32_exact(v, p)?.is_some_and(|m| m <= 0 && next_n as i32 > -m))
}

fn gauss_z_one(
    a: &ExactNumNumber,
    b: &ExactNumNumber,
    c: &ExactNumNumber,
    p: usize,
    cc: &mut Consts,
) -> Result<ExactNumNumber, Error> {
    let cab = c
        .sub(a, p, RoundingMode::None)?
        .sub(b, p, RoundingMode::None)?;
    if !cab.is_positive() {
        return Err(Error::InvalidArgument);
    }
    let gc = c.gamma_at(p, cc)?;
    let gcab = cab.gamma_at(p, cc)?;
    let gca = c.sub(a, p, RoundingMode::None)?.gamma_at(p, cc)?;
    let gcb = c.sub(b, p, RoundingMode::None)?.gamma_at(p, cc)?;
    gc.mul(&gcab, p, RoundingMode::None)?.div(
        &gca.mul(&gcb, p, RoundingMode::None)?,
        p,
        RoundingMode::None,
    )
}

fn tiny_spread(p: usize) -> Result<ExactNumNumber, Error> {
    ExactNumNumber::from_word(1, p)?.ldexp(-((p as i32) / 3 + 16), p, RoundingMode::None)
}

fn en_max(a: &ExactNumNumber, b: &ExactNumNumber) -> Result<ExactNumNumber, Error> {
    if a.cmp(b) >= 0 {
        a.clone()
    } else {
        b.clone()
    }
}

fn close_enough(dev: &ExactNumNumber, an: &ExactNumNumber, p: usize) -> Result<bool, Error> {
    let one = ExactNumNumber::from_word(1, p)?;
    let scale = en_max(&an.abs()?, &one)?;
    let thresh = tiny_spread(p)?.mul(&scale, p, RoundingMode::None)?;
    Ok(dev.cmp(&thresh) < 0)
}

fn max_dev3(
    an: &ExactNumNumber,
    x: &ExactNumNumber,
    y: &ExactNumNumber,
    z: &ExactNumNumber,
    p: usize,
) -> Result<ExactNumNumber, Error> {
    let dx = an.sub(x, p, RoundingMode::None)?.abs()?;
    let dy = an.sub(y, p, RoundingMode::None)?.abs()?;
    let dz = an.sub(z, p, RoundingMode::None)?.abs()?;
    en_max(&en_max(&dx, &dy)?, &dz)
}

fn nonnegative(x: &ExactNumNumber) -> bool {
    !x.is_negative()
}

/// \(R_C(x,y)=R_F(x,y,y)\) for \(x\ge 0\), \(y>0\).
fn carlson_rc(
    x: &ExactNumNumber,
    y: &ExactNumNumber,
    p: usize,
    cc: &mut Consts,
) -> Result<ExactNumNumber, Error> {
    if x.is_negative() || !y.is_positive() {
        return Err(Error::InvalidArgument);
    }
    if x.sub(y, p, RoundingMode::None)?
        .abs()?
        .cmp(&tiny_spread(p)?)
        < 0
    {
        return ExactNumNumber::from_word(1, p)?.div(
            &x.sqrt(p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        );
    }
    if x.is_zero() {
        let mut hp = cc.pi_num(p, RoundingMode::None)?;
        hp.div_by_2(RoundingMode::None);
        return hp.div(&y.sqrt(p, RoundingMode::None)?, p, RoundingMode::None);
    }
    if y.cmp(x) > 0 {
        let d = y
            .sub(x, p, RoundingMode::None)?
            .sqrt(p, RoundingMode::None)?;
        let arg = y
            .sub(x, p, RoundingMode::None)?
            .div(x, p, RoundingMode::None)?
            .sqrt(p, RoundingMode::None)?;
        arg.atan(p, RoundingMode::None, cc)?
            .div(&d, p, RoundingMode::None)
    } else {
        let d = x
            .sub(y, p, RoundingMode::None)?
            .sqrt(p, RoundingMode::None)?;
        let arg = x
            .sub(y, p, RoundingMode::None)?
            .div(x, p, RoundingMode::None)?
            .sqrt(p, RoundingMode::None)?;
        arg.atanh(p, RoundingMode::None, cc)?
            .div(&d, p, RoundingMode::None)
    }
}

/// Symmetric \(R_F(x,y,z)\). Arguments \(\ge 0\); at most one may be 0.
fn carlson_rf(
    x0: &ExactNumNumber,
    y0: &ExactNumNumber,
    z0: &ExactNumNumber,
    p: usize,
    _cc: &mut Consts,
) -> Result<ExactNumNumber, Error> {
    if !nonnegative(x0) || !nonnegative(y0) || !nonnegative(z0) {
        return Err(Error::InvalidArgument);
    }
    let zeros = usize::from(x0.is_zero()) + usize::from(y0.is_zero()) + usize::from(z0.is_zero());
    if zeros > 1 {
        return Err(Error::InvalidArgument);
    }
    let four = ExactNumNumber::from_word(4, p)?;
    let three = ExactNumNumber::from_word(3, p)?;
    let mut x = x0.clone()?;
    let mut y = y0.clone()?;
    let mut z = z0.clone()?;
    for _n in 0..CARLSON_DUPE_MAX {
        let an = x
            .add(&y, p, RoundingMode::None)?
            .add(&z, p, RoundingMode::None)?
            .div(&three, p, RoundingMode::None)?;
        if close_enough(&max_dev3(&an, &x, &y, &z, p)?, &an, p)? {
            return rf_series(&an, &x, &y, &z, p);
        }
        let sx = x.sqrt(p, RoundingMode::None)?;
        let sy = y.sqrt(p, RoundingMode::None)?;
        let sz = z.sqrt(p, RoundingMode::None)?;
        let lam = sx
            .mul(&sy, p, RoundingMode::None)?
            .add(&sy.mul(&sz, p, RoundingMode::None)?, p, RoundingMode::None)?
            .add(&sz.mul(&sx, p, RoundingMode::None)?, p, RoundingMode::None)?;
        x = x
            .add(&lam, p, RoundingMode::None)?
            .div(&four, p, RoundingMode::None)?;
        y = y
            .add(&lam, p, RoundingMode::None)?
            .div(&four, p, RoundingMode::None)?;
        z = z
            .add(&lam, p, RoundingMode::None)?
            .div(&four, p, RoundingMode::None)?;
    }
    Err(Error::InvalidArgument)
}

fn rf_series(
    an: &ExactNumNumber,
    x: &ExactNumNumber,
    y: &ExactNumNumber,
    z: &ExactNumNumber,
    p: usize,
) -> Result<ExactNumNumber, Error> {
    let xx = an
        .sub(x, p, RoundingMode::None)?
        .div(an, p, RoundingMode::None)?;
    let yy = an
        .sub(y, p, RoundingMode::None)?
        .div(an, p, RoundingMode::None)?;
    let zz = an
        .sub(z, p, RoundingMode::None)?
        .div(an, p, RoundingMode::None)?;
    let e2 = xx.mul(&yy, p, RoundingMode::None)?.sub(
        &zz.mul(&zz, p, RoundingMode::None)?,
        p,
        RoundingMode::None,
    )?;
    let e3 = xx
        .mul(&yy, p, RoundingMode::None)?
        .mul(&zz, p, RoundingMode::None)?;
    let e2s = e2.mul(&e2, p, RoundingMode::None)?;
    let one = ExactNumNumber::from_word(1, p)?;
    let w = |n: Word| ExactNumNumber::from_word(n, p);
    // Carlson 1995 power series once the spread is \(2^{-p/3}\).
    let s = one
        .sub(
            &e2.div(&w(10)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?
        .add(
            &e3.div(&w(14)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?
        .add(
            &e2s.div(&w(24)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?
        .sub(
            &w(3)?
                .mul(&e2, p, RoundingMode::None)?
                .mul(&e3, p, RoundingMode::None)?
                .div(&w(44)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?
        .sub(
            &w(5)?
                .mul(&e2, p, RoundingMode::None)?
                .mul(&e2s, p, RoundingMode::None)?
                .div(&w(208)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?
        .add(
            &w(3)?
                .mul(&e2s, p, RoundingMode::None)?
                .mul(&e3, p, RoundingMode::None)?
                .div(&w(104)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
    s.div(&an.sqrt(p, RoundingMode::None)?, p, RoundingMode::None)
}

/// \(R_D(x,y,z)\). \(z>0\); \(x,y\ge 0\); at most one of \(x,y\) may be 0.
fn carlson_rd(
    x0: &ExactNumNumber,
    y0: &ExactNumNumber,
    z0: &ExactNumNumber,
    p: usize,
    _cc: &mut Consts,
) -> Result<ExactNumNumber, Error> {
    if !nonnegative(x0) || !nonnegative(y0) || !z0.is_positive() {
        return Err(Error::InvalidArgument);
    }
    if x0.is_zero() && y0.is_zero() {
        return Err(Error::InvalidArgument);
    }
    let four = ExactNumNumber::from_word(4, p)?;
    let three = ExactNumNumber::from_word(3, p)?;
    let five = ExactNumNumber::from_word(5, p)?;
    let mut x = x0.clone()?;
    let mut y = y0.clone()?;
    let mut z = z0.clone()?;
    let mut sum = ExactNumNumber::from_word(0, p)?;
    let mut fac = ExactNumNumber::from_word(1, p)?;
    for _n in 0..CARLSON_DUPE_MAX {
        let an = x
            .add(&y, p, RoundingMode::None)?
            .add(
                &three.mul(&z, p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )?
            .div(&five, p, RoundingMode::None)?;
        if close_enough(&max_dev3(&an, &x, &y, &z, p)?, &an, p)? {
            let series = rd_series(&an, &x, &y, &z, p)?;
            return three.mul(&sum, p, RoundingMode::None)?.add(
                &fac.mul(&series, p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            );
        }
        let sx = x.sqrt(p, RoundingMode::None)?;
        let sy = y.sqrt(p, RoundingMode::None)?;
        let sz = z.sqrt(p, RoundingMode::None)?;
        let lam = sx
            .mul(&sy, p, RoundingMode::None)?
            .add(&sy.mul(&sz, p, RoundingMode::None)?, p, RoundingMode::None)?
            .add(&sz.mul(&sx, p, RoundingMode::None)?, p, RoundingMode::None)?;
        sum = sum.add(
            &fac.div(
                &sz.mul(&z.add(&lam, p, RoundingMode::None)?, p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )?,
            p,
            RoundingMode::None,
        )?;
        fac = fac.div(&four, p, RoundingMode::None)?;
        x = x
            .add(&lam, p, RoundingMode::None)?
            .div(&four, p, RoundingMode::None)?;
        y = y
            .add(&lam, p, RoundingMode::None)?
            .div(&four, p, RoundingMode::None)?;
        z = z
            .add(&lam, p, RoundingMode::None)?
            .div(&four, p, RoundingMode::None)?;
    }
    Err(Error::InvalidArgument)
}

fn rd_series(
    an: &ExactNumNumber,
    x: &ExactNumNumber,
    y: &ExactNumNumber,
    z: &ExactNumNumber,
    p: usize,
) -> Result<ExactNumNumber, Error> {
    let xx = an
        .sub(x, p, RoundingMode::None)?
        .div(an, p, RoundingMode::None)?;
    let yy = an
        .sub(y, p, RoundingMode::None)?
        .div(an, p, RoundingMode::None)?;
    let zz = an
        .sub(z, p, RoundingMode::None)?
        .div(an, p, RoundingMode::None)?;
    let e2 = xx.mul(&yy, p, RoundingMode::None)?.sub(
        &zz.mul(&zz, p, RoundingMode::None)?,
        p,
        RoundingMode::None,
    )?;
    let e3 = xx
        .mul(&yy, p, RoundingMode::None)?
        .mul(&zz, p, RoundingMode::None)?;
    let e2s = e2.mul(&e2, p, RoundingMode::None)?;
    let one = ExactNumNumber::from_word(1, p)?;
    let w = |n: Word| ExactNumNumber::from_word(n, p);
    let s = one
        .sub(
            &w(3)?
                .mul(&e2, p, RoundingMode::None)?
                .div(&w(14)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?
        .add(
            &e3.div(&w(6)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?
        .add(
            &w(9)?
                .mul(&e2s, p, RoundingMode::None)?
                .div(&w(88)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?
        .sub(
            &w(3)?
                .mul(&e2, p, RoundingMode::None)?
                .mul(&e3, p, RoundingMode::None)?
                .div(&w(22)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?
        .add(
            &w(9)?
                .mul(&e3, p, RoundingMode::None)?
                .mul(&e3, p, RoundingMode::None)?
                .div(&w(52)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?
        .sub(
            &w(3)?
                .mul(&e2, p, RoundingMode::None)?
                .mul(&e2s, p, RoundingMode::None)?
                .div(&w(26)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
    s.div(
        &an.mul(&an.sqrt(p, RoundingMode::None)?, p, RoundingMode::None)?,
        p,
        RoundingMode::None,
    )
}

/// \(R_J(x,y,z,p)\) for \(x,y,z\ge 0\), \(p>0\); at most one of \(x,y,z\) may be 0.
fn carlson_rj(
    x0: &ExactNumNumber,
    y0: &ExactNumNumber,
    z0: &ExactNumNumber,
    p0: &ExactNumNumber,
    p: usize,
    cc: &mut Consts,
) -> Result<ExactNumNumber, Error> {
    if !nonnegative(x0) || !nonnegative(y0) || !nonnegative(z0) || !p0.is_positive() {
        return Err(Error::InvalidArgument);
    }
    let zeros = usize::from(x0.is_zero()) + usize::from(y0.is_zero()) + usize::from(z0.is_zero());
    if zeros > 1 {
        return Err(Error::InvalidArgument);
    }
    let two = ExactNumNumber::from_word(2, p)?;
    let three = ExactNumNumber::from_word(3, p)?;
    let four = ExactNumNumber::from_word(4, p)?;
    let five = ExactNumNumber::from_word(5, p)?;
    let mut x = x0.clone()?;
    let mut y = y0.clone()?;
    let mut z = z0.clone()?;
    let mut pv = p0.clone()?;
    let mut sum = ExactNumNumber::from_word(0, p)?;
    let mut fac = ExactNumNumber::from_word(1, p)?;
    for _n in 0..CARLSON_DUPE_MAX {
        let an = x
            .add(&y, p, RoundingMode::None)?
            .add(&z, p, RoundingMode::None)?
            .add(&two.mul(&pv, p, RoundingMode::None)?, p, RoundingMode::None)?
            .div(&five, p, RoundingMode::None)?;
        let d4 = max_dev3(&an, &x, &y, &z, p)?;
        let dp = an.sub(&pv, p, RoundingMode::None)?.abs()?;
        if close_enough(&en_max(&d4, &dp)?, &an, p)? {
            let series = rj_series(&an, &x, &y, &z, &pv, p)?;
            return three.mul(&sum, p, RoundingMode::None)?.add(
                &fac.mul(&series, p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            );
        }
        let sx = x.sqrt(p, RoundingMode::None)?;
        let sy = y.sqrt(p, RoundingMode::None)?;
        let sz = z.sqrt(p, RoundingMode::None)?;
        let lam = sx
            .mul(&sy, p, RoundingMode::None)?
            .add(&sy.mul(&sz, p, RoundingMode::None)?, p, RoundingMode::None)?
            .add(&sz.mul(&sx, p, RoundingMode::None)?, p, RoundingMode::None)?;
        let alpha = pv
            .mul(
                &sx.add(&sy, p, RoundingMode::None)?
                    .add(&sz, p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )?
            .add(
                &sx.mul(&sy, p, RoundingMode::None)?
                    .mul(&sz, p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )?;
        let alpha = alpha.mul(&alpha, p, RoundingMode::None)?;
        let pl = pv.add(&lam, p, RoundingMode::None)?;
        let beta = pv
            .mul(&pl, p, RoundingMode::None)?
            .mul(&pl, p, RoundingMode::None)?;
        sum = sum.add(
            &fac.mul(&carlson_rc(&alpha, &beta, p, cc)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
        fac = fac.div(&four, p, RoundingMode::None)?;
        x = x
            .add(&lam, p, RoundingMode::None)?
            .div(&four, p, RoundingMode::None)?;
        y = y
            .add(&lam, p, RoundingMode::None)?
            .div(&four, p, RoundingMode::None)?;
        z = z
            .add(&lam, p, RoundingMode::None)?
            .div(&four, p, RoundingMode::None)?;
        pv = pv
            .add(&lam, p, RoundingMode::None)?
            .div(&four, p, RoundingMode::None)?;
    }
    Err(Error::InvalidArgument)
}

fn rj_series(
    an: &ExactNumNumber,
    x: &ExactNumNumber,
    y: &ExactNumNumber,
    z: &ExactNumNumber,
    pv: &ExactNumNumber,
    p: usize,
) -> Result<ExactNumNumber, Error> {
    let xx = an
        .sub(x, p, RoundingMode::None)?
        .div(an, p, RoundingMode::None)?;
    let yy = an
        .sub(y, p, RoundingMode::None)?
        .div(an, p, RoundingMode::None)?;
    let zz = an
        .sub(z, p, RoundingMode::None)?
        .div(an, p, RoundingMode::None)?;
    let pp = an
        .sub(pv, p, RoundingMode::None)?
        .div(an, p, RoundingMode::None)?;
    let xyz = xx
        .mul(&yy, p, RoundingMode::None)?
        .mul(&zz, p, RoundingMode::None)?;
    let xy_xz_yz = xx
        .mul(&yy, p, RoundingMode::None)?
        .add(&xx.mul(&zz, p, RoundingMode::None)?, p, RoundingMode::None)?
        .add(&yy.mul(&zz, p, RoundingMode::None)?, p, RoundingMode::None)?;
    let p2 = pp.mul(&pp, p, RoundingMode::None)?;
    let p3 = p2.mul(&pp, p, RoundingMode::None)?;
    let two = ExactNumNumber::from_word(2, p)?;
    let three = ExactNumNumber::from_word(3, p)?;
    let e2 = xy_xz_yz.sub(
        &three.mul(&p2, p, RoundingMode::None)?,
        p,
        RoundingMode::None,
    )?;
    let e3 = xyz
        .add(&two.mul(&p3, p, RoundingMode::None)?, p, RoundingMode::None)?
        .sub(
            &pp.mul(&xy_xz_yz, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
    let e2s = e2.mul(&e2, p, RoundingMode::None)?;
    let one = ExactNumNumber::from_word(1, p)?;
    let w = |n: Word| ExactNumNumber::from_word(n, p);
    let s = one
        .sub(
            &w(3)?
                .mul(&e2, p, RoundingMode::None)?
                .div(&w(14)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?
        .add(
            &e3.div(&w(6)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?
        .add(
            &w(9)?
                .mul(&e2s, p, RoundingMode::None)?
                .div(&w(88)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?
        .sub(
            &w(3)?
                .mul(&e2, p, RoundingMode::None)?
                .mul(&e3, p, RoundingMode::None)?
                .div(&w(22)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?
        .add(
            &w(9)?
                .mul(&e3, p, RoundingMode::None)?
                .mul(&e3, p, RoundingMode::None)?
                .div(&w(52)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?
        .sub(
            &w(3)?
                .mul(&e2, p, RoundingMode::None)?
                .mul(&e2s, p, RoundingMode::None)?
                .div(&w(26)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
    s.div(
        &an.mul(&an.sqrt(p, RoundingMode::None)?, p, RoundingMode::None)?,
        p,
        RoundingMode::None,
    )
}

fn extra_bits_for_degree(n: u32) -> usize {
    round_p(n as usize)
}

/// Analytic continuation of \({}_2F_1\) for real \(z\le -1\) when the value is real.
fn hypergeom_z_le_neg_one(
    a: &ExactNumNumber,
    b: &ExactNumNumber,
    c: &ExactNumNumber,
    z: &ExactNumNumber,
    p: usize,
    cc: &mut Consts,
) -> Result<ExactNumNumber, Error> {
    let one = ExactNumNumber::from_word(1, p)?;
    let two = ExactNumNumber::from_word(2, p)?;
    // 2F1(1,1;2;z) = −ln(1−z)/z for z ≠ 0, 1−z > 0.
    if a.cmp(&one) == 0 && b.cmp(&one) == 0 && c.cmp(&two) == 0 {
        let omz = one.sub(z, p, RoundingMode::None)?;
        if !omz.is_positive() {
            return Err(Error::InvalidArgument);
        }
        return omz
            .ln(p, RoundingMode::None, cc)?
            .neg()?
            .div(z, p, RoundingMode::None);
    }
    if z.cmp(&one.neg()?) == 0 {
        return hypergeom_series(a, b, c, z, p);
    }
    // 15.3.7: 1/z ∈ (−1, 0) when z < −1; (−z)^{−a} is real.
    let inv = one.div(z, p, RoundingMode::None)?;
    if inv.abs()?.cmp(&one) >= 0 {
        return Err(Error::InvalidArgument);
    }
    let mz = z.neg()?;
    let t1 = hypergeom_15_3_7_term(a, b, c, &mz, &inv, p, cc)?;
    let t2 = hypergeom_15_3_7_term(b, a, c, &mz, &inv, p, cc)?;
    t1.add(&t2, p, RoundingMode::None)
}

fn hypergeom_15_3_7_term(
    a: &ExactNumNumber,
    b: &ExactNumNumber,
    c: &ExactNumNumber,
    mz: &ExactNumNumber,
    inv: &ExactNumNumber,
    p: usize,
    cc: &mut Consts,
) -> Result<ExactNumNumber, Error> {
    let one = ExactNumNumber::from_word(1, p)?;
    let bma = b.sub(a, p, RoundingMode::None)?;
    let cma = c.sub(a, p, RoundingMode::None)?;
    let gc = c.gamma_at(p, cc)?;
    let gbma = bma.gamma_at(p, cc)?;
    let gb = b.gamma_at(p, cc)?;
    let gcma = cma.gamma_at(p, cc)?;
    let pref = gc.mul(&gbma, p, RoundingMode::None)?.div(
        &gb.mul(&gcma, p, RoundingMode::None)?,
        p,
        RoundingMode::None,
    )?;
    let mut na = a.clone()?;
    na.set_sign(if a.is_positive() { Sign::Neg } else { Sign::Pos });
    let pow = mz.pow(&na, p, RoundingMode::None, cc)?;
    let ac1 = a
        .sub(c, p, RoundingMode::None)?
        .add(&one, p, RoundingMode::None)?;
    let ab1 = a
        .sub(b, p, RoundingMode::None)?
        .add(&one, p, RoundingMode::None)?;
    let f = a.hypergeom_2f1_at(&ac1, &ab1, inv, p, cc)?;
    pref.mul(&pow, p, RoundingMode::None)?
        .mul(&f, p, RoundingMode::None)
}

pub(crate) fn series_n_max(p: usize, exp: i32) -> usize {
    let mag = if exp <= 0 {
        0
    } else if exp >= 16 {
        p
    } else {
        1usize << (exp as usize).min(12)
    };
    p.saturating_add(32).saturating_add(mag)
}

fn as_i32_exact(x: &ExactNumNumber, p: usize) -> Result<Option<i32>, Error> {
    let _ = p;
    if x.is_zero() {
        return Ok(Some(0));
    }
    if !x.fract()?.is_zero() {
        return Ok(None);
    }
    let u = match x.abs()?.int_as_usize() {
        Ok(v) => v,
        Err(_) => return Ok(None),
    };
    if u > i32::MAX as usize {
        return Ok(None);
    }
    let n = u as i32;
    Ok(Some(if x.is_negative() { -n } else { n }))
}

fn harmonic_u(n: usize, p: usize) -> Result<ExactNumNumber, Error> {
    let mut h = ExactNumNumber::from_word(0, p)?;
    for k in 1..=n {
        let t = ExactNumNumber::from_word(1, p)?.div(
            &ExactNumNumber::from_word(k as Word, p)?,
            p,
            RoundingMode::None,
        )?;
        h = h.add(&t, p, RoundingMode::None)?;
    }
    Ok(h)
}

fn half_integer_n(nu: &ExactNumNumber, p: usize) -> Result<Option<u32>, Error> {
    let two = ExactNumNumber::from_word(2, p)?;
    let two_nu = nu.mul(&two, p, RoundingMode::None)?;
    let Some(t) = as_i32_exact(&two_nu, p)? else {
        return Ok(None);
    };
    if t < 0 || t % 2 == 0 {
        return Ok(None);
    }
    Ok(Some((t as u32) / 2))
}

fn floor_nonneg_u32(x: &ExactNumNumber, p: usize) -> Result<u32, Error> {
    let _ = p;
    if x.is_negative() {
        return Err(Error::InvalidArgument);
    }
    if x.is_zero() {
        return Ok(0);
    }
    let u = x.int_as_usize()?;
    if u > u32::MAX as usize {
        return Err(Error::ExponentOverflow(Sign::Pos));
    }
    Ok(u as u32)
}

fn two_times(x: &ExactNumNumber, p: usize) -> Result<ExactNumNumber, Error> {
    x.mul(&ExactNumNumber::from_word(2, p)?, p, RoundingMode::None)
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

    fn bits_agree(a: &ExactNumNumber, b: &ExactNumNumber, p: usize, min_bits: i32, label: &str) {
        let d = a.sub(b, p, RoundingMode::None).unwrap().abs().unwrap();
        if d.is_zero() {
            return;
        }
        let rel = d.exponent() - a.exponent();
        assert!(
            rel < -min_bits,
            "{label}: relative exponent {rel} (a.exp={}, b.exp={})",
            a.exponent(),
            b.exponent()
        );
    }

    #[test]
    fn test_ei_si_li_fresnel() {
        let p = 256;
        let mut cc = Consts::new().unwrap();
        let rm = RoundingMode::ToEven;
        let one = ExactNumNumber::from_word(1, p).unwrap();
        let zero = ExactNumNumber::new(p).unwrap();

        assert!(one.ei(p, rm, &mut cc).is_ok());
        assert!(zero.ei(p, rm, &mut cc).is_err());
        assert!(one.neg().unwrap().ei(p, rm, &mut cc).is_ok());

        let s0 = zero.si(p, rm, &mut cc).unwrap();
        assert!(s0.is_zero());
        let s1 = one.si(p, rm, &mut cc).unwrap();
        let sn = one.neg().unwrap().si(p, rm, &mut cc).unwrap();
        assert!(s1.cmp(&sn.neg().unwrap()) == 0);

        let e = ExactNumNumber::from_word(1, p)
            .unwrap()
            .exp(p, rm, &mut cc)
            .unwrap();
        let li_e = e.li(p, rm, &mut cc).unwrap();
        let ei1 = one.ei(p, rm, &mut cc).unwrap();
        bits_agree(&li_e, &ei1, p, (p as i32) / 4, "li(e)=Ei(1)");

        assert!(one.li(p, rm, &mut cc).is_err());
        assert!(zero.ci(p, rm, &mut cc).is_err());
        assert!(one.ci(p, rm, &mut cc).is_ok());

        assert!(zero.fresnel_s(p, rm, &mut cc).unwrap().is_zero());
        assert!(zero.fresnel_c(p, rm, &mut cc).unwrap().is_zero());
        let fs = one.fresnel_s(p, rm, &mut cc).unwrap();
        let fns = one.neg().unwrap().fresnel_s(p, rm, &mut cc).unwrap();
        assert!(fs.cmp(&fns.neg().unwrap()) == 0);

        let forty = ExactNumNumber::from_word(40, p).unwrap();
        bits_agree(
            &forty.si(p, rm, &mut cc).unwrap(),
            &forty.si(128, rm, &mut cc).unwrap(),
            p,
            80,
            "Si(40)",
        );
        bits_agree(
            &forty.ci(p, rm, &mut cc).unwrap(),
            &forty.ci(128, rm, &mut cc).unwrap(),
            p,
            80,
            "Ci(40)",
        );

        let three_hundred = ExactNumNumber::from_word(300, p).unwrap();
        let si300 = three_hundred.si(p, rm, &mut cc).unwrap();
        let mut half_pi = cc.pi_num(p, RoundingMode::None).unwrap();
        half_pi.div_by_2(RoundingMode::None);
        let gap = half_pi
            .sub(&si300, p, RoundingMode::None)
            .unwrap()
            .abs()
            .unwrap();
        assert!(gap.exponent() < -6, "Si(300) within O(1/x) of π/2");
        bits_agree(
            &si300,
            &three_hundred.si(128, rm, &mut cc).unwrap(),
            p,
            80,
            "Si(300)",
        );

        let sixty_four = ExactNumNumber::from_word(64, p).unwrap();
        bits_agree(
            &sixty_four.si(p, rm, &mut cc).unwrap(),
            &sixty_four.si(64, rm, &mut cc).unwrap(),
            p,
            40,
            "Si(64) series vs f,g",
        );
        bits_agree(
            &sixty_four.ci(p, rm, &mut cc).unwrap(),
            &sixty_four.ci(64, rm, &mut cc).unwrap(),
            p,
            40,
            "Ci(64) series vs f,g",
        );
        let eight = ExactNumNumber::from_word(8, p).unwrap();
        bits_agree(
            &eight.fresnel_s(p, rm, &mut cc).unwrap(),
            &eight.fresnel_s(64, rm, &mut cc).unwrap(),
            p,
            40,
            "S(8) series vs f,g",
        );
        bits_agree(
            &three_hundred.ei(p, rm, &mut cc).unwrap(),
            &three_hundred.ei(128, rm, &mut cc).unwrap(),
            p,
            80,
            "Ei(300)",
        );

        let twenty = ExactNumNumber::from_word(20, p).unwrap();
        let fs20 = twenty.fresnel_s(p, rm, &mut cc).unwrap();
        let fc20 = twenty.fresnel_c(p, rm, &mut cc).unwrap();
        let half = one_half(p).unwrap();
        let ds = half
            .sub(&fs20, p, RoundingMode::None)
            .unwrap()
            .abs()
            .unwrap();
        let dc = half
            .sub(&fc20, p, RoundingMode::None)
            .unwrap()
            .abs()
            .unwrap();
        assert!(ds.exponent() < -4, "S(20) near 1/2");
        assert!(dc.exponent() < -4, "C(20) near 1/2");
        bits_agree(
            &fs20,
            &twenty.fresnel_s(128, rm, &mut cc).unwrap(),
            p,
            80,
            "S(20)",
        );

        let half_ext = crate::ExactNum::from_u8(1, p).div(&crate::ExactNum::from_u8(2, p), p, rm);
        let inf_si = crate::INF_POS.si(p, rm, &mut cc);
        let hp = cc.pi(p, rm).div(&crate::ExactNum::from_u8(2, p), p, rm);
        let d = inf_si.sub(&hp, p, rm).abs();
        assert!(d.is_zero() || d.exponent().unwrap_or(0) < -((p as i32) / 4));
        let d = crate::INF_POS
            .fresnel_s(p, rm, &mut cc)
            .sub(&half_ext, p, rm)
            .abs();
        assert!(d.is_zero() || d.exponent().unwrap_or(0) < -((p as i32) / 4));
        let d = crate::INF_POS
            .fresnel_c(p, rm, &mut cc)
            .sub(&half_ext, p, rm)
            .abs();
        assert!(d.is_zero() || d.exponent().unwrap_or(0) < -((p as i32) / 4));
    }

    #[test]
    fn test_digamma_gammainc() {
        let p = 256;
        let mut cc = Consts::new().unwrap();
        let rm = RoundingMode::ToEven;
        let one = ExactNumNumber::from_word(1, p).unwrap();
        let two = ExactNumNumber::from_word(2, p).unwrap();
        let zero = ExactNumNumber::new(p).unwrap();
        let g = cc.euler_gamma_num(p, RoundingMode::None).unwrap();
        let half = one_half(p).unwrap();

        assert!(zero.digamma(p, rm, &mut cc).is_err());
        let neg_half = half.neg().unwrap();
        let psi_nh = neg_half.digamma(p, rm, &mut cc).unwrap();
        // ψ(1/2) = −γ − 2 ln 2; ψ(−1/2) = ψ(1/2) − π cot(−π/2) wait: 1−(−1/2)=3/2
        // ψ(−1/2) = ψ(3/2) − π cot(−π/2). cot(−π/2)=0, ψ(3/2)=ψ(1/2)+2
        let three_half = ExactNumNumber::from_word(3, p)
            .unwrap()
            .div(&two, p, RoundingMode::None)
            .unwrap();
        let psi_32 = three_half.digamma(p, rm, &mut cc).unwrap();
        bits_agree(&psi_nh, &psi_32, p, 40, "psi(-1/2)=psi(3/2)");
        let psi1 = one.digamma(p, rm, &mut cc).unwrap();
        let r1 = psi1.add(&g, p, RoundingMode::None).unwrap().abs().unwrap();
        assert!(r1.is_zero() || r1.exponent() < -80, "psi(1)+γ");

        let psi2 = two.digamma(p, rm, &mut cc).unwrap();
        let r2 = psi2
            .add(&g, p, RoundingMode::None)
            .unwrap()
            .sub(&one, p, RoundingMode::None)
            .unwrap()
            .abs()
            .unwrap();
        assert!(r2.is_zero() || r2.exponent() < -80, "psi(2)+γ-1");

        let psih = half.digamma(p, rm, &mut cc).unwrap();
        let ln2 = cc.ln_2_num(p, RoundingMode::None).unwrap();
        let two_ln2 = two.mul(&ln2, p, RoundingMode::None).unwrap();
        let want = g
            .neg()
            .unwrap()
            .sub(&two_ln2, p, RoundingMode::None)
            .unwrap();
        bits_agree(&psih, &want, p, 80, "psi(1/2)");

        let five = ExactNumNumber::from_word(5, p).unwrap();
        assert!(five.gammainc(&zero, p, rm, &mut cc).unwrap().is_zero());
        let g11 = one.gammainc(&one, p, rm, &mut cc).unwrap();
        let en = one.neg().unwrap().exp(p, rm, &mut cc).unwrap();
        let want = one.sub(&en, p, RoundingMode::None).unwrap();
        bits_agree(&g11, &want, p, 80, "γ(1,1)");
        let gu = one.gammainc_upper(&one, p, rm, &mut cc).unwrap();
        bits_agree(&gu, &en, p, 80, "Γ(1,1)=e^{-1}");
        assert!(one.gammainc(&one.neg().unwrap(), p, rm, &mut cc).is_err());
    }

    #[test]
    fn test_bessel_family() {
        let p = 256;
        let mut cc = Consts::new().unwrap();
        let rm = RoundingMode::ToEven;
        let one = ExactNumNumber::from_word(1, p).unwrap();
        let zero = ExactNumNumber::new(p).unwrap();
        let half = one_half(p).unwrap();
        let two = ExactNumNumber::from_word(2, p).unwrap();

        let j0 = zero.bessel_j_nu(&zero, p, rm, &mut cc).unwrap();
        bits_agree(&j0, &one, p, 80, "J_0(0)");

        let i0 = zero.bessel_i(&zero, p, rm, &mut cc).unwrap();
        bits_agree(&i0, &one, p, 80, "I_0(0)");

        assert!(zero.bessel_y(&zero, p, rm, &mut cc).is_err());
        assert!(zero.bessel_k(&zero, p, rm, &mut cc).is_err());

        let k_half = one.bessel_k(&half, p, rm, &mut cc).unwrap();
        let pi = cc.pi_num(p, RoundingMode::None).unwrap();
        let mut nx = one.clone().unwrap();
        nx.set_sign(Sign::Neg);
        let want = pi
            .div(&two, p, RoundingMode::None)
            .unwrap()
            .sqrt(p, RoundingMode::None)
            .unwrap()
            .mul(&nx.exp(p, rm, &mut cc).unwrap(), p, RoundingMode::None)
            .unwrap();
        bits_agree(&k_half, &want, p, 80, "K_{1/2}(1)");

        let k_neg = one.bessel_k(&half.neg().unwrap(), p, rm, &mut cc).unwrap();
        bits_agree(&k_half, &k_neg, p, 80, "K_{-1/2}=K_{1/2}");

        let k0 = one.bessel_k(&zero, p, rm, &mut cc).unwrap();
        let k1 = one.bessel_k(&one, p, rm, &mut cc).unwrap();
        let i0 = one.bessel_i(&zero, p, rm, &mut cc).unwrap();
        let i1 = one.bessel_i(&one, p, rm, &mut cc).unwrap();
        let w = i0
            .mul(&k1, p, RoundingMode::None)
            .unwrap()
            .add(
                &i1.mul(&k0, p, RoundingMode::None).unwrap(),
                p,
                RoundingMode::None,
            )
            .unwrap();
        bits_agree(&w, &one, p, 40, "Wronskian at 1");

        let j2000 = zero.bessel_j(2000, p, rm, &mut cc).unwrap();
        assert!(j2000.is_zero(), "J_2000(0)=0");
        let k50 = one
            .bessel_k(&ExactNumNumber::from_word(50, p).unwrap(), p, rm, &mut cc)
            .unwrap();
        let k49 = one
            .bessel_k(&ExactNumNumber::from_word(49, p).unwrap(), p, rm, &mut cc)
            .unwrap();
        let k51 = one
            .bessel_k(&ExactNumNumber::from_word(51, p).unwrap(), p, rm, &mut cc)
            .unwrap();
        // K_{ν+1} = (2ν/x) K_ν + K_{ν−1}
        let rec = ExactNumNumber::from_word(100, p)
            .unwrap()
            .mul(&k50, p, RoundingMode::None)
            .unwrap()
            .add(&k49, p, RoundingMode::None)
            .unwrap();
        bits_agree(&k51, &rec, p, 20, "K recurrence at 50");
    }

    #[test]
    fn test_elliptic_integrals() {
        let p = 256;
        let mut cc = Consts::new().unwrap();
        let rm = RoundingMode::ToEven;
        let one = ExactNumNumber::from_word(1, p).unwrap();
        let zero = ExactNumNumber::new(p).unwrap();
        let two = ExactNumNumber::from_word(2, p).unwrap();
        let three = ExactNumNumber::from_word(3, p).unwrap();
        let four = ExactNumNumber::from_word(4, p).unwrap();
        let half = one_half(p).unwrap();
        let pi = cc.pi_num(p, RoundingMode::None).unwrap();
        let mut half_pi = pi.clone().unwrap();
        half_pi.div_by_2(RoundingMode::None);

        let k0 = zero.elliptic_k(p, rm, &mut cc).unwrap();
        bits_agree(&k0, &half_pi, p, 80, "K(0)");
        let e0 = zero.elliptic_e_complete(p, rm, &mut cc).unwrap();
        bits_agree(&e0, &half_pi, p, 80, "E(0)");
        let e1 = one.elliptic_e_complete(p, rm, &mut cc).unwrap();
        bits_agree(&e1, &one, p, 80, "E(1)");
        assert!(matches!(
            one.elliptic_k(p, rm, &mut cc),
            Err(Error::ExponentOverflow(Sign::Pos))
        ));

        let k_half = half.elliptic_k(p, rm, &mut cc).unwrap();
        let q = ExactNumNumber::from_word(1, p)
            .unwrap()
            .div(&four, p, RoundingMode::None)
            .unwrap();
        let k4 = four.elliptic_k(p, rm, &mut cc).unwrap();
        let kq = q.elliptic_k(p, rm, &mut cc).unwrap();
        let want4 = kq.div(&two, p, RoundingMode::None).unwrap();
        bits_agree(&k4, &want4, p, 40, "K(4)=K(1/4)/2");
        let g14 = q.gamma(p, rm, &mut cc).unwrap();
        let want_k = g14
            .mul(&g14, p, RoundingMode::None)
            .unwrap()
            .div(
                &four
                    .mul(
                        &pi.sqrt(p, RoundingMode::None).unwrap(),
                        p,
                        RoundingMode::None,
                    )
                    .unwrap(),
                p,
                RoundingMode::None,
            )
            .unwrap();
        bits_agree(&k_half, &want_k, p, 80, "K(1/2)");

        let e_half = half.elliptic_e_complete(p, rm, &mut cc).unwrap();
        let want_e = pi
            .div(
                &four.mul(&k_half, p, RoundingMode::None).unwrap(),
                p,
                RoundingMode::None,
            )
            .unwrap()
            .add(
                &k_half.div(&two, p, RoundingMode::None).unwrap(),
                p,
                RoundingMode::None,
            )
            .unwrap();
        bits_agree(&e_half, &want_e, p, 80, "E(1/2)");

        let x = half.clone().unwrap();
        let f0 = x.elliptic_f(&zero, p, rm, &mut cc).unwrap();
        bits_agree(&f0, &x.asin(p, rm, &mut cc).unwrap(), p, 80, "F(x|0)");
        let f1 = x.elliptic_f(&one, p, rm, &mut cc).unwrap();
        bits_agree(&f1, &x.atanh(p, rm, &mut cc).unwrap(), p, 80, "F(x|1)");
        let ei1 = x.elliptic_e(&one, p, rm, &mut cc).unwrap();
        bits_agree(&ei1, &x, p, 80, "E(x|1)");
        assert!(one.elliptic_f(&one, p, rm, &mut cc).is_err());

        let m13 = one.div(&three, p, RoundingMode::None).unwrap();
        let k13 = m13.elliptic_k(p, rm, &mut cc).unwrap();
        let f1 = one.elliptic_f(&m13, p, rm, &mut cc).unwrap();
        bits_agree(&k13, &f1, p, 80, "K=F(1)");
        let e13 = m13.elliptic_e_complete(p, rm, &mut cc).unwrap();
        let ei = one.elliptic_e(&m13, p, rm, &mut cc).unwrap();
        bits_agree(&e13, &ei, p, 80, "E=E(1)");

        let pim = zero.elliptic_pi_complete(&half, p, rm, &mut cc).unwrap();
        bits_agree(&pim, &k_half, p, 80, "Π(0,m)=K");
        let pi0 = zero.elliptic_pi(&x, &half, p, rm, &mut cc).unwrap();
        let fx = x.elliptic_f(&half, p, rm, &mut cc).unwrap();
        bits_agree(&pi0, &fx, p, 80, "Π(0;x|m)=F");

        let n = one.div(&four, p, RoundingMode::None).unwrap();
        let got = n.elliptic_pi(&x, &zero, p, rm, &mut cc).unwrap();
        let root = one
            .sub(&n, p, RoundingMode::None)
            .unwrap()
            .sqrt(p, RoundingMode::None)
            .unwrap();
        let want = root
            .mul(&x, p, RoundingMode::None)
            .unwrap()
            .div(
                &one.sub(
                    &x.mul(&x, p, RoundingMode::None).unwrap(),
                    p,
                    RoundingMode::None,
                )
                .unwrap()
                .sqrt(p, RoundingMode::None)
                .unwrap(),
                p,
                RoundingMode::None,
            )
            .unwrap()
            .atan(p, rm, &mut cc)
            .unwrap()
            .div(&root, p, RoundingMode::None)
            .unwrap();
        bits_agree(&got, &want, p, 80, "Π(n;x|0)");

        let pmm = m13.elliptic_pi_complete(&m13, p, rm, &mut cc).unwrap();
        let want = e13
            .div(
                &one.sub(&m13, p, RoundingMode::None).unwrap(),
                p,
                RoundingMode::None,
            )
            .unwrap();
        bits_agree(&pmm, &want, p, 40, "Π(m,m)");

        assert!(one.elliptic_pi_complete(&half, p, rm, &mut cc).is_err());
        assert!(two.elliptic_f(&zero, p, rm, &mut cc).is_err());
    }

    #[test]
    fn test_legendre_hypergeom_betainc() {
        let p = 256;
        let mut cc = Consts::new().unwrap();
        let rm = RoundingMode::ToEven;
        let one = ExactNumNumber::from_word(1, p).unwrap();
        let zero = ExactNumNumber::new(p).unwrap();
        let two = ExactNumNumber::from_word(2, p).unwrap();
        let three = ExactNumNumber::from_word(3, p).unwrap();
        let half = one_half(p).unwrap();

        for n in 0u32..=8 {
            let p1 = one.legendre_p(n, p, rm).unwrap();
            bits_agree(&p1, &one, p, 80, "P_n(1)");
            let pm = one.neg().unwrap().legendre_p(n, p, rm).unwrap();
            let want = if n % 2 == 0 { one.clone().unwrap() } else { one.neg().unwrap() };
            bits_agree(&pm, &want, p, 80, "P_n(-1)");
        }
        let p2 = zero.legendre_p(2, p, rm).unwrap();
        bits_agree(&p2, &half.neg().unwrap(), p, 80, "P_2(0)");
        let p100 = one.legendre_p(100, p, rm).unwrap();
        bits_agree(&p100, &one, p, 80, "P_100(1)");
        let p128 = one.neg().unwrap().legendre_p(128, p, rm).unwrap();
        bits_agree(&p128, &one, p, 80, "P_128(-1)");
        let xh = half.clone().unwrap();
        let pn = xh.legendre_p(128, p, rm).unwrap();
        let pn1 = xh.legendre_p(127, p, rm).unwrap();
        let pn2 = xh.legendre_p(126, p, rm).unwrap();
        // (n+1) P_{n+1} = (2n+1) x P_n − n P_{n−1}  at n=127
        let n = ExactNumNumber::from_word(127, p).unwrap();
        let np1 = ExactNumNumber::from_word(128, p).unwrap();
        let two_n1 = ExactNumNumber::from_word(255, p).unwrap();
        let lhs = np1.mul(&pn, p, RoundingMode::None).unwrap();
        let rhs = two_n1
            .mul(&xh, p, RoundingMode::None)
            .unwrap()
            .mul(&pn1, p, RoundingMode::None)
            .unwrap()
            .sub(
                &n.mul(&pn2, p, RoundingMode::None).unwrap(),
                p,
                RoundingMode::None,
            )
            .unwrap();
        bits_agree(&lhs, &rhs, p, 20, "Bonnet n=127");

        let p11 = half.assoc_legendre_p(1, 1, p, rm).unwrap();
        let want = one
            .sub(
                &half.mul(&half, p, RoundingMode::None).unwrap(),
                p,
                RoundingMode::None,
            )
            .unwrap()
            .sqrt(p, RoundingMode::None)
            .unwrap()
            .neg()
            .unwrap();
        bits_agree(&p11, &want, p, 80, "P_1^1(1/2)");

        let f0 = half
            .hypergeom_2f1(&half, &one, &zero, p, rm, &mut cc)
            .unwrap();
        bits_agree(&f0, &one, p, 80, "2F1(...,0)");

        let nm2 = two.neg().unwrap();
        let fpoly = nm2
            .hypergeom_2f1(&one, &one, &half, p, rm, &mut cc)
            .unwrap();
        bits_agree(
            &fpoly,
            &one.div(&four_word(p), p, RoundingMode::None).unwrap(),
            p,
            80,
            "2F1(-2,1;1;1/2)",
        );

        let fln = one
            .hypergeom_2f1(&one, &two, &half, p, rm, &mut cc)
            .unwrap();
        let ln2 = cc.ln_2_num(p, RoundingMode::None).unwrap();
        bits_agree(
            &fln,
            &two.mul(&ln2, p, RoundingMode::None).unwrap(),
            p,
            80,
            "2F1(1,1;2;1/2)",
        );

        let fg = nm2
            .hypergeom_2f1(&one, &three, &one, p, rm, &mut cc)
            .unwrap();
        bits_agree(&fg, &half, p, 80, "2F1(-2,1;3;1)");

        let fm1 = one
            .hypergeom_2f1(&one, &two, &one.neg().unwrap(), p, rm, &mut cc)
            .unwrap();
        bits_agree(&fm1, &ln2, p, 40, "2F1(1,1;2;-1)=ln2");
        let fterm2 = nm2
            .hypergeom_2f1(&one, &three, &two, p, rm, &mut cc)
            .unwrap();
        let third = one.div(&three, p, RoundingMode::None).unwrap();
        bits_agree(&fterm2, &third, p, 80, "2F1(-2,1;3;2)=1/3");

        assert!(one
            .hypergeom_2f1(&one, &three, &two, p, rm, &mut cc)
            .is_err());

        let fk = half
            .hypergeom_2f1(&half, &one, &half, p, rm, &mut cc)
            .unwrap();
        let k = half.elliptic_k(p, rm, &mut cc).unwrap();
        let pi = cc.pi_num(p, RoundingMode::None).unwrap();
        let want = two
            .div(&pi, p, RoundingMode::None)
            .unwrap()
            .mul(&k, p, RoundingMode::None)
            .unwrap();
        bits_agree(&fk, &want, p, 40, "2F1 = 2K/π");

        assert!(half
            .betainc(&half, &zero, p, rm, &mut cc)
            .unwrap()
            .is_zero());
        let i1 = half.betainc(&half, &one, p, rm, &mut cc).unwrap();
        bits_agree(&i1, &one, p, 80, "I_1(a,b)");
        // I_{1/2}(1,1) = 1/2
        let i11 = one.betainc(&one, &half, p, rm, &mut cc).unwrap();
        bits_agree(&i11, &half, p, 80, "I_{1/2}(1,1)");
        assert!(zero.betainc(&one, &half, p, rm, &mut cc).is_err());
    }

    fn four_word(p: usize) -> ExactNumNumber {
        ExactNumNumber::from_word(4, p).unwrap()
    }
}
