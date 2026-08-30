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

    /// Exponential integral `Ei(self)` for `self > 0`.
    pub fn ei(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        if !self.is_positive() {
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
        let lnx = self.ln(p, RoundingMode::None, cc)?;
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

    /// Logarithmic integral `li(self) = Ei(ln self)` for `self > 1`.
    pub fn li(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        let one = Self::from_word(1, p)?;
        if self.cmp(&one) <= 0 {
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
        let mut pow_x = if sine {
            self.mul(&x2, p, RoundingMode::None)?
        } else {
            self.clone()?
        };
        let mut fact = Self::from_word(1, p)?;
        let mut pi2_pow = if sine {
            pi2.clone()?
        } else {
            Self::from_word(1, p)?
        };
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
                    .mul(&Self::from_word((2 * np1) as Word, p)?, p, RoundingMode::None)?
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
                    .mul(&Self::from_word((2 * np1) as Word, p)?, p, RoundingMode::None)?;
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

fn series_n_max(p: usize, exp: i32) -> usize {
    let mag = if exp <= 0 {
        0
    } else if exp >= 16 {
        p
    } else {
        1usize << (exp as usize).min(12)
    };
    p.saturating_add(32).saturating_add(mag)
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

    fn bits_agree(a: &ExactNumNumber, b: &ExactNumNumber, _p: usize, min_bits: i32, label: &str) {
        let d = a
            .sub(b, p, RoundingMode::None)
            .unwrap()
            .abs()
            .unwrap();
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
        assert!(one.neg().unwrap().ei(p, rm, &mut cc).is_err());

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
        let ds = half.sub(&fs20, p, RoundingMode::None).unwrap().abs().unwrap();
        let dc = half.sub(&fc20, p, RoundingMode::None).unwrap().abs().unwrap();
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
}
