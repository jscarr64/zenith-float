//! Error function, gamma, integral specials, and Bessel J_n (integer order).

use crate::common::util::bump_prec_retry;
use crate::common::util::round_p;
use crate::defs::Error;
use crate::defs::Exponent;
use crate::defs::RoundingMode;
use crate::defs::EXPONENT_MAX;
use crate::defs::EXPONENT_MIN;
use crate::integer::ExactInt;
use crate::num::ExactNumNumber;
use crate::ops::consts::Consts;
use crate::ops::util::compute_small_exp;
use crate::ExactNum;
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

        // erf(x) = ±(1 - erfc|x|) with erfc|x| below half an ulp of 1: only its sign matters.
        if erf_tail_below(self, p) {
            let mut one = Self::from_word(1, p)?;
            one.set_sign(self.sign());
            return Self::round_offset(&one, self.sign().invert(), p, rm);
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

        // x ≥ 2^31: erfc(x) < e^{-x²} < 2^-(2^62·log₂e) lies below the whole subnormal grid.
        // It is positive, so it never rounds to -0.
        if self.is_positive() && !self.is_zero() && self.exponent() >= ERFC_UNDERFLOW_EXPONENT {
            let mut r = Self::tiny_positive(p, rm)?;
            r.set_inexact(true);
            return Ok(r);
        }

        // erfc(x) = 2 - erfc|x| for very negative x, with the tail below half an ulp of 2.
        if self.is_negative() && erf_tail_below(self, p) {
            let two = Self::from_word(2, p)?;
            return Self::round_offset(&two, Sign::Neg, p, rm);
        }

        // erfc(x) = 1 - 2x/√π + ... for tiny x: the correction is below half an ulp of 1.
        let p_exp = Exponent::try_from(p).map_err(|_| Error::InvalidArgument)?;
        if !self.is_zero() && self.exponent() < -(p_exp + 2) {
            let one = Self::from_word(1, p)?;
            return Self::round_offset(&one, self.sign().invert(), p, rm);
        }

        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc;

        loop {
            let p_x = p_wrk + WORD_BIT_SIZE;
            let mut ret = self.erfc_at(p_x, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    /// A positive value below half the smallest subnormal, rounded at `p` bits with `rm`:
    /// `+0`, or the smallest positive subnormal when `rm` rounds away from zero or to odd.
    pub(crate) fn tiny_positive(p: usize, rm: RoundingMode) -> Result<Self, Error> {
        match rm {
            RoundingMode::Up | RoundingMode::FromZero | RoundingMode::ToOdd => {
                Self::min_positive(p)
            }
            RoundingMode::None
            | RoundingMode::Down
            | RoundingMode::ToZero
            | RoundingMode::ToEven => Self::new2(p, Sign::Pos, true),
        }
    }

    /// Rounds `base + δ` once at `p` bits with `rm`, where `base` is exact and nonzero and `δ`
    /// has sign `toward` and a magnitude strictly between 0 and half an ulp of `base` at `p`
    /// bits. Any such `δ` rounds the same way, so a proxy `δ` 128 bits down stands in for it.
    /// The result is marked inexact.
    pub(crate) fn round_offset(
        base: &Self,
        toward: Sign,
        p: usize,
        rm: RoundingMode,
    ) -> Result<Self, Error> {
        let pw = round_p(p) + 3 * WORD_BIT_SIZE;
        let e0 = base.exponent();
        let mut b = Self::new(pw)?;
        b = b.add(base, pw, RoundingMode::None)?;
        b.set_exponent(0);
        let mut t = Self::from_word(1, pw)?;
        let shift = Exponent::try_from(round_p(p) + 2 * WORD_BIT_SIZE)
            .map_err(|_| Error::InvalidArgument)?;
        t.set_exponent(1 - shift);
        t.set_sign(toward);
        let mut r = b.add(&t, pw, RoundingMode::None)?;
        r.set_precision(p, rm)?;
        let e1 = r.exponent();
        r.set_exponent(e0 + e1);
        r.set_inexact(true);
        Ok(r)
    }

    fn erf_at(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let mut x = self.clone()?;
        x.set_inexact(false);
        let neg = x.is_negative();
        x.set_sign(Sign::Pos);

        // |x| large enough that erfc underflows at this precision.
        if x.exponent().saturating_mul(2) > p as i64 + 4 {
            let mut one = Self::from_word(1, p)?;
            if neg {
                one.set_sign(Sign::Neg);
            }
            return Ok(one);
        }

        let ret = match erf_x2_bits(&x, p)? {
            None => {
                let erfc = x.erfc_asymptotic(p, cc)?;
                let one = Self::from_word(1, p)?;
                one.sub(&erfc, p, RoundingMode::None)?
            }
            // The alternating series peaks near e^{x²}: carry x²·log₂e extra bits.
            Some(guard) => x.erf_series(p + guard, cc)?,
        };

        if neg {
            ret.neg()
        } else {
            Ok(ret)
        }
    }

    fn erfc_at(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let mut x = self.clone()?;
        x.set_inexact(false);
        if x.is_positive() && !x.is_zero() {
            match erf_x2_bits(&x, p)? {
                None => return x.erfc_asymptotic(p, cc),
                Some(guard) => {
                    // erfc(x) ≈ e^{-x²}/(x√π): 1 − erf(x) cancels about x²·log₂e bits.
                    let pw = p + guard;
                    let one = Self::from_word(1, pw)?;
                    return one.sub(&x.erf_at(pw, cc)?, pw, RoundingMode::None);
                }
            }
        }
        let one = Self::from_word(1, p)?;
        one.sub(&x.erf_at(p, cc)?, p, RoundingMode::None)
    }

    fn erf_series(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        // (2/√π) ∑ (-1)^n x^{2n+1} / (n! (2n+1))
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let sqrt_pi = pi.sqrt(p, RoundingMode::None)?;
        let two = Self::from_word(2, p)?;
        let scale = two.div(&sqrt_pi, p, RoundingMode::None)?;

        let x2 = self.mul(self, p, RoundingMode::None)?;
        // x^{2n+1} and the partial sum both start at x itself: None stands for `self` (no copy of it).
        let mut xpow: Option<Self> = None;
        let mut nfact = Self::from_word(1, p)?;
        let mut sum: Option<Self> = None;

        // Terms grow until n ≈ x², so only stop once past that point.
        // The term count is that peak plus the tail bound, not a fixed cap (P9).
        let mut n: usize = 1;
        loop {
            nfact = nfact.mul(&from_usize_p(n, p)?, p, RoundingMode::None)?;
            let xp = xpow
                .as_ref()
                .unwrap_or(self)
                .mul(&x2, p, RoundingMode::None)?;
            let two_n_1 = from_usize_p(2 * n + 1, p)?;
            let den = nfact.mul(&two_n_1, p, RoundingMode::None)?;
            let mut t = xp.div(&den, p, RoundingMode::None)?;
            if n % 2 == 1 {
                t.set_sign(Sign::Neg);
            }
            let sm = sum
                .as_ref()
                .unwrap_or(self)
                .add(&t, p, RoundingMode::None)?;
            let past_peak = from_usize_p(n, p)?.cmp(&x2) > 0;
            if t.is_zero() || (past_peak && t.exponent() + (p as i64) < sm.exponent()) {
                return scale.mul(&sm, p, RoundingMode::None);
            }
            xpow = Some(xp);
            sum = Some(sm);
            n = n.checked_add(1).ok_or(Error::PrecisionRetryExhausted)?;
        }
    }

    /// Asymptotic `erfc`, only used once `x² ≥ p·ln 2` (see [`erf_x2_bits`]), where the
    /// smallest term (≈ e^{-x²}) is already below `2^-p`.
    fn erfc_asymptotic(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        // e^{-x^2} / (x √π) * ∑ (-1)^k (2k-1)!! / (2x^2)^k
        let x2 = self.mul(self, p, RoundingMode::None)?;
        let nx2 = x2.with_sign(Sign::Neg)?;
        let expm = nx2.exp(p, RoundingMode::None, cc)?;
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let sqrt_pi = pi.sqrt(p, RoundingMode::None)?;
        let den = self.mul(&sqrt_pi, p, RoundingMode::None)?;
        let pre = expm.div(&den, p, RoundingMode::None)?;

        let two_x2 = x2.mul(&Self::from_word(2, p)?, p, RoundingMode::None)?;
        let mut s = Self::from_word(1, p)?;
        let mut term = Self::from_word(1, p)?;
        // Cap or divergence without meeting the bound: never use the partial sum.
        let mut done = false;
        for k in 1..=(p + 8) {
            let odd = Self::from_word((2 * k - 1) as Word, p)?;
            let ratio = odd.div(&two_x2, p, RoundingMode::None)?;
            if ratio.cmp(&Self::from_word(1, p)?) >= 0 {
                // Past the smallest term: the series starts to diverge.
                break;
            }
            term = term.mul(&ratio, p, RoundingMode::None)?;
            term.set_sign(if k % 2 == 1 { Sign::Neg } else { Sign::Pos });
            s = s.add(&term, p, RoundingMode::None)?;
            if term.is_zero() || term.exponent() + (p as i64) < 0 {
                done = true;
                break;
            }
            term.set_sign(Sign::Pos);
        }
        if !done {
            return Err(Error::PrecisionRetryExhausted);
        }
        pre.mul(&s, p, RoundingMode::None)
    }

    /// Gamma function `Γ(self)` at precision `p`.
    pub fn gamma(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;

        // Γ(±0) = ±Inf exactly, with the sign of the zero (MPFR).
        if self.is_zero() {
            return Err(Error::ExponentOverflow(self.sign()));
        }

        if !self.is_positive() && self.is_int() {
            return Err(Error::InvalidArgument);
        }

        // Tiny |x|: Γ(x) = 1/x - γ + O(x), and γ·|x| is far below every breakpoint gap of 1/x
        // (an m-bit x keeps 1/x at least 2^-(m+p+1) away from p-bit breakpoints unless 1/x is
        // exact), so the Ziv loop could never separate them. Round once instead.
        let m = self.mantissa_max_bit_len();
        let tiny_exp = -Exponent::try_from(m + p + 8).map_err(|_| Error::InvalidArgument)?;
        if self.exponent() < tiny_exp {
            let pw = round_p(m + p + 3 * WORD_BIT_SIZE);
            let one = Self::from_word(1, pw)?;
            let q = one.div(self, pw, RoundingMode::ToZero)?;
            let mut r = if q.inexact() {
                let mut t = q;
                t.set_precision(p, rm)?;
                t
            } else {
                Self::round_offset(&q, Sign::Neg, p, rm)?
            };
            r.set_inexact(true);
            return Ok(r);
        }

        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(m) + p_inc;

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

    /// Beta function `B(self, b) = Γ(self) Γ(b) / Γ(self + b)` at precision `p`, rounded once.
    ///
    /// Poles of either argument are `InvalidArgument`. Large positive arguments go through
    /// `exp(ln Γ(a) + ln Γ(b) − ln Γ(a + b))` with guard bits for the absolute error of the sum.
    pub fn beta(
        &self,
        b: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        if self.is_beta_pole() || b.is_beta_pole() {
            return Err(Error::InvalidArgument);
        }
        let inexact = self.inexact() | b.inexact();
        let one = Self::from_word(1, WORD_BIT_SIZE)?;

        // B(a, 1) = 1/a and B(1, b) = 1/b: one division, exact when the quotient is.
        for (x, y) in [(self, b), (b, self)] {
            if y.cmp(&one) == 0 {
                let mut r = Self::from_word(1, p)?.div(x, p, rm)?;
                r.set_inexact(r.inexact() | inexact);
                return Ok(r);
            }
        }

        // ZF-37: one positive integer n gives the rational B(a, n) = (n−1)! / (a (a+1) ⋯ (a+n−1))
        // for any non-pole a (DLMF 5.12.1 with Γ(a+n) = (a)_n Γ(a)), e.g. B(−½, 2) = −4, which
        // the Ziv loop over the Γ ratio cannot certify. Exact product, one division; if that width is unallocatable the floating path decides.
        for (x, y) in [(self, b), (b, self)] {
            if !y.is_positive() || !y.is_int() {
                continue;
            }
            let Ok(n) = y.int_as_usize() else { continue };
            // Exact product is a method switch, not a refusal: take it only when
            // the proven bit length of (n−1)! and (x)_n fits in p+2 words.
            // n·per is allocatable under MAX_PRECISION=2^{60}−1 even at n=2^{20},
            // and then the loop of exact muls does not finish.
            if !exact_beta_product_fits(p, n, Some(x)) {
                continue;
            }
            let mut ex = x.clone()?;
            ex.set_inexact(false);
            let num = exact_product(1, n - 1)?;
            // den = ex (ex+1) ⋯ (ex+n−1); `None` = the product so far is `ex` itself.
            let mut den: Option<Self> = None;
            for k in 1..n {
                let t = exact_add_num(&ex, &Self::from_word(k as Word, WORD_BIT_SIZE)?)?;
                den = Some(exact_mul_num(den.as_ref().unwrap_or(&ex), &t)?);
            }
            let mut r = num.div(den.as_ref().unwrap_or(&ex), p, rm)?;
            r.set_inexact(r.inexact() | inexact);
            return Ok(r);
        }

        // Positive integers: B(m, n) = (n−1)! / (m (m+1) ⋯ (m+n−1)), both sides exact.
        // Same method switch as the one-integer path: if the proven factorial
        // width does not fit in p+2 words, the Γ/Ziv path below is used.
        if let (Some(m), Some(n)) = (self.beta_int_index()?, b.beta_int_index()?) {
            if exact_beta_product_fits(p, n, None) && exact_beta_ints_fit(p, m, n) {
                let num = exact_product(1, n - 1)?;
                let den = exact_product(m, m + n - 1)?;
                let mut r = num.div(&den, p, rm)?;
                r.set_inexact(r.inexact() | inexact);
                return Ok(r);
            }
        }

        let large = self.is_positive()
            && b.is_positive()
            && (self.exponent() > BETA_RATIO_EXP || b.exponent() > BETA_RATIO_EXP);
        // |ln Γ(x)| < x·log2(x) for x ≥ 2: that many bits of the log-sum are not fraction bits.
        let guard = if large {
            let e = self.exponent().max(b.exponent()) + 1;
            let lz = Exponent::BITS - e.leading_zeros();
            usize::try_from(e).map_err(|_| Error::InvalidArgument)?
                + usize::try_from(lz).map_err(|_| Error::InvalidArgument)?
                + 4
        } else {
            0
        };
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p
            .max(self.mantissa_max_bit_len())
            .max(b.mantissa_max_bit_len())
            + p_inc;

        loop {
            let p_x = round_p(p_wrk + WORD_BIT_SIZE * 2 + guard);
            let s = self.add(b, p_x, RoundingMode::None)?;
            let mut ret = if large {
                let la = self.ln_gamma_positive_at(p_x, cc)?;
                let lb = b.ln_gamma_positive_at(p_x, cc)?;
                let ls = s.ln_gamma_positive_at(p_x, cc)?;
                let l = la
                    .add(&lb, p_x, RoundingMode::None)?
                    .sub(&ls, p_x, RoundingMode::None)?;
                let e = match l.exp(p_x, RoundingMode::None, cc) {
                    Err(Error::ExponentOverflow(Sign::Neg)) => Self::new(p_x)?,
                    other => other?,
                };
                if e.is_zero() {
                    // B > 0 far below the subnormal range.
                    let mut r = Self::tiny_positive(p, rm)?;
                    r.set_inexact(true);
                    return Ok(r);
                }
                e
            } else {
                if s.is_beta_pole() {
                    // Γ(a + b) has a pole: B(a, b) = 0.
                    return Self::new2(p, Sign::Pos, inexact);
                }
                let ga = self.gamma_at(p_x, cc)?;
                let gb = b.gamma_at(p_x, cc)?;
                let gs = s.gamma_at(p_x, cc)?;
                ga.mul(&gb, p_x, RoundingMode::None)?
                    .div(&gs, p_x, RoundingMode::None)?
            };
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | inexact);
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    /// Zero or a negative integer: a pole of Γ.
    fn is_beta_pole(&self) -> bool {
        self.is_zero() || (self.is_negative() && self.is_int())
    }

    /// `Some(n)` when `self` is a positive integer that fits in `usize`, else `None`.
    /// Larger indices belong on the Γ / ln Γ path, not `InvalidArgument`.
    fn beta_int_index(&self) -> Result<Option<usize>, Error> {
        if !self.is_positive() || self.is_zero() || !self.is_int() {
            return Ok(None);
        }
        match self.int_as_usize() {
            Ok(n) if n >= 1 => Ok(Some(n)),
            Ok(_) => Ok(None),
            Err(Error::InvalidArgument) => Ok(None),
            Err(e) => Err(e),
        }
    }

    /// `ln Γ(self)` for positive `self` at precision `p`.
    pub fn ln_gamma(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;

        // ln|Γ| has a pole at 0 (either sign): an exact +Inf in every mode.
        if self.is_zero() {
            return Err(Error::ExponentOverflow(Sign::Pos));
        }

        // Negative x: reflect with the fractional part g ∈ (0, 1/2], so that
        // |sin(πx)| = sin(πg) is formed without cancellation near the poles.
        let reflect = if self.is_negative() {
            let pw = self.mantissa_max_bit_len() + 2 * WORD_BIT_SIZE;
            let x = Self::new(pw)?.add(self, pw, RoundingMode::None)?;
            let f = x.fract()?.abs()?;
            if f.is_zero() {
                // Non-positive integer: a pole of Γ, so ln|Γ| is +Inf.
                return Err(Error::DivisionByZero);
            }
            let one = Self::from_word(1, pw)?;
            let half = one_half(pw)?;
            let g = if f.cmp(&half) > 0 { one.sub(&f, pw, RoundingMode::None)? } else { f };
            let omx = one.sub(self, pw, RoundingMode::None)?;
            Some((g, omx))
        } else {
            None
        };

        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc;

        loop {
            let p_x = p_wrk + WORD_BIT_SIZE * 2;
            let mut ret = match &reflect {
                // ln|Γ(x)| = ln π - ln sin(πg) - ln Γ(1 - x)
                Some((g, omx)) => {
                    let pi = cc.pi_num(p_x, RoundingMode::None)?;
                    let s = pi
                        .mul(g, p_x, RoundingMode::None)?
                        .sin(p_x, RoundingMode::None, cc)?;
                    let ln_pi = pi.ln(p_x, RoundingMode::None, cc)?;
                    let ln_s = s.ln(p_x, RoundingMode::None, cc)?;
                    ln_pi.sub(&ln_s, p_x, RoundingMode::None)?.sub(
                        &omx.ln_gamma_positive_at(p_x, cc)?,
                        p_x,
                        RoundingMode::None,
                    )?
                }
                None => self.ln_gamma_positive_at(p_x, cc)?,
            };
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    /// Riemann `ζ(self)` for any real `self` at precision `p`, rounded once with `rm`.
    ///
    /// ## Errors
    ///
    ///  - DivisionByZero: `self == 1`, the pole (POLE-08: ζ(1) = +Inf, divideByZero, exact; the caller maps it to +Inf).
    ///  - InvalidArgument: the precision is incorrect.
    ///  - ExponentOverflow: the result overflows (very negative `self`).
    pub(crate) fn zeta_real(
        &self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;

        let one = Self::from_word(1, p)?;
        if self.is_zero() {
            let mut h = one_half(p)?;
            h.set_sign(Sign::Neg);
            return Ok(h);
        }
        if self.cmp(&one) == 0 {
            return Err(Error::DivisionByZero);
        }
        if self.is_negative() && self.is_even_int() {
            // Trivial zeros ζ(-2k) are +0.
            return Self::new2(p, Sign::Pos, false);
        }

        let p_exp = Exponent::try_from(p).map_err(|_| Error::InvalidArgument)?;
        // s > p + 4: ζ(s) = 1 + 2^-s + ... with the tail below half an ulp of 1.
        let lim_bits = usize::BITS - (p + 4).leading_zeros();
        let lim_exp = Exponent::from(lim_bits);
        if self.is_positive() && self.exponent() > lim_exp {
            return Self::round_offset(&one, Sign::Pos, p, rm);
        }
        // |s| tiny: ζ(s) = -1/2 - s·ln(2π)/2 + ... with the correction below half an ulp.
        if self.exponent() < -(p_exp + 4) {
            let mut h = one_half(p)?;
            h.set_sign(Sign::Neg);
            return Self::round_offset(&h, self.sign().invert(), p, rm);
        }

        let half = one_half(p)?;
        let right = self.cmp(&half) >= 0;

        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc;

        loop {
            let p_x = p_wrk + WORD_BIT_SIZE * 2;
            let mut ret = if right { self.zeta_em(p_x, cc)? } else { self.zeta_reflect(p_x, cc)? };
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(true);
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    /// `ζ(s) = 2^s π^(s-1) sin(πs/2) Γ(1-s) ζ(1-s)` for `s < 1/2`, at working precision `p`.
    ///
    /// Every factor is formed at `p + g` bits with `g ≈ log2|s| + log2 log2|s|` (ZF-11): an
    /// error ε in an argument of size `|s|` moves `2^s`, `π^(s−1)` and `Γ(1−s)` by about
    /// `|s| ln|s| ε`. `1 − s` is kept exact, and `sin(πs/2)` is taken of `s mod 4` (exact),
    /// so no factor sees the rounding of a large argument. Past `|s| = 2^20` the product is
    /// formed as `exp` of the summed logarithms, so no factor under- or overflows when the
    /// result does not; for `|s| ≥ 2^64` the result is beyond `EXPONENT_MAX` and overflows
    /// with the sign of `sin(πs/2)`.
    fn zeta_reflect(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let es = usize::try_from(Exponent::max(0, self.exponent()))
            .map_err(|_| Error::InvalidArgument)?;
        let g = es + (usize::BITS - es.leading_zeros()) as usize + 8;
        let pg = round_p(p + g);
        // 1 - s exactly.
        let frac = Exponent::max(0, -self.exponent());
        let frac_bits = usize::try_from(frac).map_err(|_| Error::InvalidArgument)?;
        let pw = pg + self.mantissa_max_bit_len() + frac_bits + WORD_BIT_SIZE;
        let one = Self::from_word(1, pw)?;
        let omz = one.sub(self, pw, RoundingMode::None)?;

        // sin(πs/2) = sin(π t/2), t = s rem 4 exactly (|t| < 4).
        let t = self.rem(&Self::from_word(4, WORD_BIT_SIZE)?)?;
        let pi = cc.pi_num(pg, RoundingMode::None)?;
        let mut arg = pi.mul(&t, pg, RoundingMode::None)?;
        arg.div_by_2(RoundingMode::None);
        let sn = arg.sin(pg, RoundingMode::None, cc)?;
        if es >= 64 {
            // log2|ζ(s)| ≥ |s|(log2|s| − 5) > EXPONENT_MAX = 2^62 − 1.
            return Err(Error::ExponentOverflow(sn.sign()));
        }
        let ln2 = Self::from_word(2, pg)?.ln(pg, RoundingMode::None, cc)?;
        let ln_pi = pi.ln(pg, RoundingMode::None, cc)?;
        let sm1 = omz.neg()?;
        let z = omz.zeta_em(pg, cc)?;
        if es > 20 {
            let mut l = self.mul(&ln2, pg, RoundingMode::None)?;
            l = l.add(
                &sm1.mul(&ln_pi, pg, RoundingMode::None)?,
                pg,
                RoundingMode::None,
            )?;
            l = l.add(
                &sn.abs()?.ln(pg, RoundingMode::None, cc)?,
                pg,
                RoundingMode::None,
            )?;
            l = l.add(
                &omz.ln_gamma(pg, RoundingMode::None, cc)?,
                pg,
                RoundingMode::None,
            )?;
            l = l.add(&z.ln(pg, RoundingMode::None, cc)?, pg, RoundingMode::None)?;
            let mut r = match l.exp(pg, RoundingMode::None, cc) {
                Err(Error::ExponentOverflow(_)) => return Err(Error::ExponentOverflow(sn.sign())),
                other => other?,
            };
            r.set_sign(sn.sign());
            return Ok(r);
        }
        let two_s = self
            .mul(&ln2, pg, RoundingMode::None)?
            .exp(pg, RoundingMode::None, cc)?;
        let pi_s = sm1
            .mul(&ln_pi, pg, RoundingMode::None)?
            .exp(pg, RoundingMode::None, cc)?;
        let gm = omz.gamma(pg, RoundingMode::None, cc)?;
        two_s
            .mul(&pi_s, pg, RoundingMode::None)?
            .mul(&sn, pg, RoundingMode::None)?
            .mul(&gm, pg, RoundingMode::None)?
            .mul(&z, pg, RoundingMode::None)
    }

    /// Euler–Maclaurin `ζ(s)` for `s ≥ 1/2`, `s ≠ 1`, at working precision `p`:
    /// `Σ_{n<N} n^-s + N^(1-s)/(s-1) + N^-s/2 + Σ_k B_2k/(2k)! · s(s+1)…(s+2k-2) · N^(-s-2k+1)`.
    fn zeta_em(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        // ZF-11: the cut below grows as 2s (about 5·10^7 trial-divided terms for the reflected
        // ζ(2.6·10^7)). For s ≥ p + 16 the direct sum needs two terms: the tail after 2^-s is
        // Σ_{n≥3} n^-s ≤ 3^-s (1 + 3/(s-1)) ≤ 2·3^-s < 2^-(p+16).
        let big = Self::from_usize(p + 16)?;
        if self.cmp(&big) >= 0 {
            let two_ms = Self::from_word(2, p)?
                .ln(p, RoundingMode::None, cc)?
                .mul(self, p, RoundingMode::None)?
                .neg()?
                .exp(p, RoundingMode::None, cc)?;
            return Self::from_word(1, p)?.add(&two_ms, p, RoundingMode::None);
        }
        // FLINT 3.6 path: choose_param / hurwitz_em_jet (value component), with a
        // sieved power sum at a = 1. Replaces the old n_cut = p/4 cut.
        let s = ExactNum::from(self.clone()?);
        let a = ExactNum::from_u8(1, p);
        match crate::hurwitz_jet::hurwitz_em_jet(&s, &a, p, cc) {
            Some((z, _)) => z.try_to_number(),
            None => Err(Error::PrecisionRetryExhausted),
        }
    }

    /// `ln Γ(self)` for positive `self` at working precision `p`.
    fn ln_gamma_positive_at(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        // Large arguments go straight to Stirling: Γ itself would overflow.
        let e = self.exponent();
        if e >= stirling_min_exponent(p, 372, 127, 6) {
            return self.ln_gamma_stirling(p, cc);
        }
        // ZF-11: the 64-term threshold grows like p/127, so a wide working precision (a
        // caller's wide-mantissa argument such as the exact `1 − s` of the ζ reflection)
        // sent z ≈ 2^62 to Γ(z), which overflows although ln Γ(z) ≈ 2^67.6 is finite. With
        // more terms the Stirling tail meets 2^-(p+8) at a smaller exponent: use the
        // digamma term count when it suffices, and for z ≥ 2^47 (where Γ can overflow at
        // EXPONENT_MAX) double the count until it does; never fall back to Γ there.
        let pw = p + WORD_BIT_SIZE;
        let mut terms = digamma_terms(p);
        if terms > 64 && e >= asymp_min_exponent(pw, terms) {
            return with_even_bernoulli(terms, p, |bern| self.ln_gamma_stirling_with(p, bern, cc));
        }
        if e >= 48 {
            while e < asymp_min_exponent(pw, terms) {
                terms = terms.checked_mul(2).ok_or(Error::MemoryAllocation)?;
            }
            return with_even_bernoulli(terms, p, |bern| self.ln_gamma_stirling_with(p, bern, cc));
        }
        self.ln_gamma_small_at(p, cc)
    }

    /// `ln Γ(x)` for `0 < x < 2^48` off the Stirling range (ZF-30). `ln(Γ(x))` at `p` bits
    /// has absolute error `2^-p`, which is a large relative error where `ln Γ` is small (near
    /// `x = 1` and `x = 2`), and `Γ(x) ≈ 1/x` overflows for the tiniest `x`. So:
    /// - `x < 2^-(2p+64)`: `ln Γ(x) = −ln x − γx` (next term `ζ(2)x²/2`, below `2^-2p` relative);
    /// - `x < 1/2`: `ln Γ(x) = ln Γ(1 + x) − ln x` with `1 + x` exact (no cancellation);
    /// - `|x − 1|` or `|x − 2| < 2^-(2p+64)`: `−γ(x−1)` resp. `(1−γ)(x−2)`;
    /// - otherwise `ln Γ(x)` from `Γ(x)` at `p + log2(1/|x−1|, 1/|x−2|) + 64` bits.
    fn ln_gamma_small_at(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let nr = RoundingMode::None;
        let e = self.exponent();
        let tiny = -Exponent::try_from(2 * p + 64).map_err(|_| Error::InvalidArgument)?;
        if e < tiny {
            let g = cc.euler_gamma_num(p, nr)?;
            let lx = self.ln(p, nr, cc)?;
            return lx.neg()?.sub(&g.mul(self, p, nr)?, p, nr);
        }
        let mant = self.mantissa_max_bit_len();
        if e <= -1 {
            let w = round_p(mant + usize::try_from(-e).unwrap_or(0) + 2 * WORD_BIT_SIZE);
            let y = Self::from_word(1, WORD_BIT_SIZE)?.add(self, w, nr)?;
            let pw = p + WORD_BIT_SIZE;
            let v = y.ln_gamma_small_at(pw, cc)?;
            return v.sub(&self.ln(pw, nr, cc)?, p, nr);
        }
        let w = round_p(mant + usize::try_from(e.max(0)).unwrap_or(0) + 2 * WORD_BIT_SIZE);
        let t1 = self.sub(&Self::from_word(1, WORD_BIT_SIZE)?, w, nr)?;
        let t2 = self.sub(&Self::from_word(2, WORD_BIT_SIZE)?, w, nr)?;
        if t1.is_zero() || t2.is_zero() {
            return Self::new(p);
        }
        let g = cc.euler_gamma_num(p + WORD_BIT_SIZE, nr)?;
        if t1.exponent() < tiny {
            return g.mul(&t1, p, nr)?.neg();
        }
        if t2.exponent() < tiny {
            let c = Self::from_word(1, WORD_BIT_SIZE)?.sub(&g, p + WORD_BIT_SIZE, nr)?;
            return c.mul(&t2, p, nr);
        }
        let ext = usize::try_from((-t1.exponent()).max(-t2.exponent()).max(0)).unwrap_or(0);
        let pp = round_p(p + ext + WORD_BIT_SIZE);
        let mut v = self.gamma_at(pp, cc)?.ln(pp, nr, cc)?;
        v.set_precision(p, nr)?;
        Ok(v)
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
        // Term count and shift target come from the remainder bound
        // (`gamma_stirling_plan` / `asymp_min_exponent`). The step count is
        // 2^{e_min-1}+2, the walk needed to reach 2^{e_min-1}, not a
        // separate 2^20 cap (P9).
        let e_min = gamma_stirling_plan(p).1;
        self.gamma_positive_capped(p, cc, stirling_shift_count(e_min))
    }

    /// [`Self::gamma_positive`] with the Stirling shift cap as a parameter (tests force it).
    fn gamma_positive_capped(
        &self,
        p: usize,
        cc: &mut Consts,
        shift_max: usize,
    ) -> Result<Self, Error> {
        // Γ(n) = (n−1)! exactly whenever the factorial is at most 4p + 4096 bits wide (ZF-37:
        // the old 1..=64 limit left e.g. Γ(171) = 170!, 852 significant bits and exactly
        // representable at p = 1024, to the Stirling path, whose Ziv loop can never certify an
        // exact value and exhausted after ~5 s). Wider factorials go to Stirling under Ziv.
        if self.is_int() && self.is_positive() && self.exponent() <= 32 {
            let n = self.int_as_usize()?;
            let m = n.saturating_sub(1);
            let bitlen = (usize::BITS - m.leading_zeros()) as usize;
            let bits = m.saturating_mul(bitlen);
            if bits <= p.saturating_mul(4).saturating_add(4096) {
                let w = round_p(bits.max(p) + WORD_BIT_SIZE);
                let mut f = factorial(m, w)?;
                f.set_inexact(false);
                return Ok(f);
            }
        }

        // Half-integers: Γ(n+½) = √π · ∏_{k=1}^{n} (2k−1)/2, with the product exact.
        if let Some(n) = self.half_integer_index(p)? {
            let mut prod = Self::from_word(1, WORD_BIT_SIZE)?;
            for k in 1..=n {
                let pb = round_p(prod.mantissa_max_bit_len() + 2 * WORD_BIT_SIZE);
                prod = prod.mul(&from_usize_p(2 * k - 1, pb)?, pb, RoundingMode::None)?;
            }
            let e = prod.exponent() - Exponent::try_from(n).map_err(|_| Error::InvalidArgument)?;
            prod.set_exponent(e);
            let pw = p + WORD_BIT_SIZE;
            let rpi = cc
                .pi_num(pw, RoundingMode::None)?
                .sqrt(pw, RoundingMode::None)?;
            return rpi.mul(&prod, p, RoundingMode::None);
        }

        // One extra word covers the shift product and exp() turning the absolute error of
        // ln Γ into a relative one.
        let pw = p + WORD_BIT_SIZE;
        let one = Self::from_word(1, pw)?;
        let mut acc = Self::from_word(1, pw)?;

        // Raise the argument until the last Stirling term is below 2^-(p+8).
        // Term count and e_min come from [`gamma_stirling_plan`].
        let (terms, e_min) = gamma_stirling_plan(p);
        // A caller-forced cap below the walk that e_min needs returns
        // PrecisionRetryExhausted at once (tests). The production path
        // passes stirling_shift_count(e_min), not a fixed 2^20 cap.
        let shifted = shift_up(self, e_min, shift_max, pw, |z| {
            acc = acc.mul(z, pw, RoundingMode::None)?;
            z.add(&one, pw, RoundingMode::None)
        })?;
        let z = shifted.as_ref().unwrap_or(self);

        let lg = if terms == 64 {
            z.ln_gamma_stirling(pw, cc)?
        } else {
            with_even_bernoulli(terms, pw, |bern| z.ln_gamma_stirling_with(pw, bern, cc))?
        };
        let g = lg.exp(pw, RoundingMode::None, cc)?;
        g.div(&acc, p, RoundingMode::None)
    }

    /// `Some(n)` when `self = n + 1/2` and the closed-form product or sum is
    /// within the same bit-length / term bound as the integer factorial path
    /// (`n ≤ 4p + 4096`). Larger half-integers use the Stirling / asymptotic kernel.
    fn half_integer_index(&self, p: usize) -> Result<Option<usize>, Error> {
        if !self.is_positive() || self.is_zero() || self.is_int() {
            return Ok(None);
        }
        if !self.is_int_when_doubled() {
            return Ok(None);
        }
        let n = match self.int_as_usize() {
            Ok(n) => n,
            Err(Error::InvalidArgument) => return Ok(None),
            Err(e) => return Err(e),
        };
        let cap = p.saturating_mul(4).saturating_add(4096);
        Ok(if n <= cap { Some(n) } else { None })
    }

    fn ln_gamma_stirling(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let bern = even_bernoulli(64, p)?;
        self.ln_gamma_stirling_with(p, &bern, cc)
    }

    /// Stirling series for `ln Γ(self)` with the even Bernoulli numbers `bern` (B₂, B₄, …):
    /// as many terms as `bern` holds, `PrecisionRetryExhausted` if none falls below 2^-p.
    fn ln_gamma_stirling_with(
        &self,
        p: usize,
        bern: &[Self],
        cc: &mut Consts,
    ) -> Result<Self, Error> {
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

        // z^{2k-1}: None stands for z^1 = `self` itself at k = 1 (no copy).
        let mut zpow: Option<Self> = None;
        // Cap or divergence without meeting the bound: never use the partial sum.
        let mut done = false;
        for (i, b) in bern.iter().enumerate() {
            let k = i + 1;
            let two_k = Self::from_word((2 * k) as Word, p)?;
            let two_k_m1 = Self::from_word((2 * k - 1) as Word, p)?;
            let den = two_k.mul(&two_k_m1, p, RoundingMode::None)?.mul(
                zpow.as_ref().unwrap_or(self),
                p,
                RoundingMode::None,
            )?;
            let term = b.div(&den, p, RoundingMode::None)?;
            s = s.add(&term, p, RoundingMode::None)?;
            if term.is_zero() || term.exponent() + (p as i64) < 0 {
                done = true;
                break;
            }
            zpow = Some(
                zpow.as_ref()
                    .unwrap_or(self)
                    .mul(self, p, RoundingMode::None)?
                    .mul(self, p, RoundingMode::None)?,
            );
        }
        if !done {
            return Err(Error::PrecisionRetryExhausted);
        }
        Ok(s)
    }

    /// Digamma \(\psi(\mathrm{self})=\Gamma'/\Gamma\). Poles at non-positive integers.
    /// For \(z<0\) uses \(\psi(z)=\psi(1-z)-\pi\cot(\pi z)\).
    pub fn digamma(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        if self.is_zero() || (!self.is_positive() && self.is_int()) {
            return Err(Error::InvalidArgument);
        }
        // Tiny |x|: ψ(x) = -1/x - γ + O(x); as for Γ, the Ziv loop cannot separate it from
        // -1/x, so round -1/x once toward -γ.
        let m = self.mantissa_max_bit_len();
        let tiny_exp = -Exponent::try_from(m + p + 8).map_err(|_| Error::InvalidArgument)?;
        if self.exponent() < tiny_exp {
            let pw = round_p(m + p + 3 * WORD_BIT_SIZE);
            let mone = Self::from_word(1, pw)?.neg()?;
            let q = mone.div(self, pw, RoundingMode::ToZero)?;
            let mut r = if q.inexact() {
                let mut t = q;
                t.set_precision(p, rm)?;
                t
            } else {
                Self::round_offset(&q, Sign::Neg, p, rm)?
            };
            r.set_inexact(true);
            return Ok(r);
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(m) + p_inc;
        let mut lost_done = 0;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE * 2;
            let mut ret = self.digamma_at(p_x, cc)?;
            // Near the root x₀ ≈ 1.4616 the shift sum and the asymptotic part cancel: the
            // error is absolute, so a result of exponent -k has k fewer correct bits.
            let lost = if ret.exponent() < 0 {
                usize::try_from(-ret.exponent()).map_err(|_| Error::InvalidArgument)?
            } else {
                0
            };
            if lost > lost_done {
                p_wrk = round_p(p_wrk + lost - lost_done);
                lost_done = lost;
                continue;
            }
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
            // cot(πz) has period 1: reduce with the exact fractional part of z. π·z itself at
            // `p` bits keeps only `p - exponent(z)` fraction bits, a large absolute error in
            // the cot argument once z has many integer bits.
            let fz = z.clone()?.fract()?;
            if fz.is_zero() {
                return Err(Error::InvalidArgument);
            }
            let piz = pi.mul(&fz, p, RoundingMode::None)?;
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
        // Term count and shift target come from the remainder bound
        // (`digamma_terms` / `asymp_min_exponent`). The step count is
        // 2^{e_min-1}+2, the walk needed to reach 2^{e_min-1}, not a
        // separate 2^20 cap (P9).
        let pw = p + WORD_BIT_SIZE;
        let e_min = asymp_min_exponent(pw, digamma_terms(p));
        self.digamma_positive_capped(p, cc, stirling_shift_count(e_min))
    }

    /// [`Self::digamma_positive`] with the shift cap as a parameter (tests force it).
    fn digamma_positive_capped(
        &self,
        p: usize,
        cc: &mut Consts,
        shift_max: usize,
    ) -> Result<Self, Error> {
        // One extra word absorbs the rounding of the shift sum.
        let pw = p + WORD_BIT_SIZE;
        let one = Self::from_word(1, pw)?;
        // Half-integers: ψ(n+½) = −γ − 2 ln 2 + Σ_{k=1}^{n} 2/(2k−1).
        if let Some(n) = self.half_integer_index(p)? {
            let g = cc.euler_gamma_num(pw, RoundingMode::None)?;
            let mut ln4 = cc.ln_2_num(pw, RoundingMode::None)?;
            ln4.set_exponent(ln4.exponent() + 1);
            let mut acc = Self::new(pw)?.sub(&g, pw, RoundingMode::None)?;
            acc = acc.sub(&ln4, pw, RoundingMode::None)?;
            let two = Self::from_word(2, pw)?;
            for k in 1..=n {
                let t = two.div(&from_usize_p(2 * k - 1, pw)?, pw, RoundingMode::None)?;
                acc = acc.add(&t, pw, RoundingMode::None)?;
            }
            return acc.add(&Self::new(p)?, p, RoundingMode::None);
        }
        // Huge z: the series past -1/(2z) is below 2^-(2·pw) relative; z² would overflow.
        if self.exponent() > Exponent::try_from(pw).map_err(|_| Error::InvalidArgument)? {
            let ln_z = self.ln(pw, RoundingMode::None, cc)?;
            let mut half_inv = one.div(self, pw, RoundingMode::None)?;
            half_inv.set_exponent(half_inv.exponent() - 1);
            return ln_z.sub(&half_inv, p, RoundingMode::None);
        }
        let mut z0 = self.clone()?;
        z0.set_inexact(false);
        let mut acc = Self::from_word(0, pw)?;
        // Raise until the last asymptotic term (K = digamma_terms(p) Bernoulli numbers)
        // is below 2^-(p+8). More terms at higher p keep the shift count small.
        let terms = digamma_terms(p);
        // Bound at `pw`, the precision digamma_asymp checks its stop rule against.
        let e_min = asymp_min_exponent(pw, terms);
        let shifted = shift_up(&z0, e_min, shift_max, pw, |z| {
            let rec = one.div(z, pw, RoundingMode::None)?;
            acc = acc.sub(&rec, pw, RoundingMode::None)?;
            z.add(&one, pw, RoundingMode::None)
        })?;
        let z = shifted.as_ref().unwrap_or(&z0);
        acc.add(&z.digamma_asymp(pw, terms, cc)?, p, RoundingMode::None)
    }

    fn digamma_asymp(&self, p: usize, terms: usize, cc: &mut Consts) -> Result<Self, Error> {
        // ln z − 1/(2z) − Σ_{k≥1} B_{2k} / (2k z^{2k})
        let ln_z = self.ln(p, RoundingMode::None, cc)?;
        let two = Self::from_word(2, p)?;
        let two_z = two.mul(self, p, RoundingMode::None)?;
        let half_inv = Self::from_word(1, p)?.div(&two_z, p, RoundingMode::None)?;
        let mut s = ln_z.sub(&half_inv, p, RoundingMode::None)?;
        let z2 = self.mul(self, p, RoundingMode::None)?;
        let mut zp = Self::from_word(1, p)?;
        let mut prev_e = Exponent::MIN;
        let p_exp = Exponent::try_from(p).map_err(|_| Error::InvalidArgument)?;
        // Borrowed from the shared exact-Bernoulli cache (no per-call recurrence).
        with_even_bernoulli(terms, p, |bern| {
            // Cap or divergence without meeting the bound: never use the partial sum.
            let mut done = false;
            for (i, b) in bern.iter().enumerate() {
                let k = i + 1;
                zp = zp.mul(&z2, p, RoundingMode::None)?;
                let two_k = from_usize_p(2 * k, p)?;
                let den = two_k.mul(&zp, p, RoundingMode::None)?;
                let term = b.div(&den, p, RoundingMode::None)?;
                s = s.sub(&term, p, RoundingMode::None)?;
                let e = term.exponent();
                if term.is_zero() || e + p_exp < 0 {
                    done = true;
                    break;
                }
                if k > 2 && e > prev_e {
                    break;
                }
                prev_e = e;
            }
            if !done {
                return Err(Error::PrecisionRetryExhausted);
            }
            Ok(s)
        })
    }

    /// Lambert `W₀(self)` (`minus1 = false`, `self ≥ −1/e`) or `W₋₁(self)` (`minus1 = true`,
    /// `−1/e ≤ self < 0`), correctly rounded (ZF-31).
    ///
    /// Newton on `w e^w = x`, or on `w + ln|w| = ln|x|` away from the branch point and the
    /// origin (no overflow or underflow of `e^w` for huge `|w|`), then a rigorous bracket: the
    /// equation is monotone on each branch, so opposite signs at `w ± δ` (each sign certified
    /// against a bound on its own rounding error) put `W` within `δ` and Ziv rounds. Near
    /// `−1/e` the working precision grows by `log2(1/|1 + w|)`. Dyadic `x ≠ 0` never gives a
    /// dyadic `W` (`e^w` is transcendental for dyadic `w ≠ 0`), so only `W(0) = 0` is exact.
    pub fn lambert_w(
        &self,
        minus1: bool,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        if self.is_zero() {
            if minus1 {
                return Err(Error::InvalidArgument);
            }
            return Self::new2(p, self.sign(), self.inexact());
        }
        if minus1 && self.is_positive() {
            return Err(Error::InvalidArgument);
        }
        let mut x = self.clone()?;
        x.set_inexact(false);
        // Domain: 1 + e x >= 0. x is dyadic and −1/e is not, so the sign is never 0; widen e
        // until the sign of 1 + e x is certain.
        let one = Self::from_word(1, WORD_BIT_SIZE)?;
        if x.is_negative() {
            let mut pp = round_p(p.max(x.mantissa_max_bit_len()) + 2 * WORD_BIT_SIZE);
            loop {
                let e = cc.e_num(pp, RoundingMode::None)?;
                let t = one.add(&e.mul(&x, pp, RoundingMode::None)?, pp, RoundingMode::None)?;
                // |error| <= 2^(2 - pp) (e and the product each within 2^-pp relative, |e x| < 2)
                if !t.is_zero() && t.exponent() > 3 - pp as Exponent {
                    if t.is_negative() {
                        return Err(Error::InvalidArgument);
                    }
                    break;
                }
                pp = pp.checked_mul(2).ok_or(Error::PrecisionRetryExhausted)?;
                if pp > p.saturating_mul(64).max(1 << 16) {
                    return Err(Error::PrecisionRetryExhausted);
                }
            }
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(x.mantissa_max_bit_len()) + p_inc;
        let guess = lambert_w_guess(&x, minus1, cc)?;
        loop {
            if let Some(mut ret) = lambert_w_bracketed(&x, minus1, &guess, p_wrk, cc)? {
                if ret.try_set_precision(p, rm, p_wrk)? {
                    ret.set_inexact(true);
                    return Ok(ret);
                }
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
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
            // ZF-37: for x > s, γ(s, x) = Γ(s) − Γ(s, x) with the asymptotic Γ(s, x) (DLMF
            // 8.11.2) whenever it converges; Γ(s, x) ≤ Γ(s)/2 there, so at most a bit or two
            // cancel. The power series alone needs about x terms (γ(125, 2^20) exhausted its
            // ~2p-term cap). If Γ(s) is exact (integer s) and Γ(s, x) is below half an ulp
            // of it, only the side is left to decide.
            if x.cmp(self) > 0 {
                let up = match self.gammainc_upper_asymptotic(x, p_x, cc) {
                    Ok(v) => v,
                    Err(Error::ExponentOverflow(Sign::Neg)) => Some(Self::new(p_x)?),
                    Err(e) => return Err(e),
                };
                if let Some(up) = up {
                    let g = self.gamma_at(p_x, cc)?;
                    let below = up.is_zero()
                        || up.exponent() < g.exponent().saturating_sub(p_x as Exponent + 2);
                    if !g.inexact() && below {
                        let mut r = Self::round_offset(&g, Sign::Neg, p, rm)?;
                        r.set_inexact(true);
                        return Ok(r);
                    }
                    let mut ret = g.sub(&up, p_x, RoundingMode::None)?;
                    if ret.try_set_precision(p, rm, p_wrk)? {
                        ret.set_inexact(true);
                        return Ok(ret);
                    }
                    bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
                    continue;
                }
            }
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
        // Γ(0, x) = E₁(x) = −Ei(−x); Γ(0) itself is a pole, so Γ(s) − γ(s, x) cannot be used.
        if self.is_zero() {
            if x.is_zero() {
                return Err(Error::ExponentOverflow(Sign::Pos));
            }
            // x > 0 here, so -x is x with the sign set negative.
            let nx = x.neg()?;
            let mut r = nx
                .ei(p, crate::common::util::invert_rm_for_sign(rm), cc)?
                .neg()?;
            r.set_inexact(r.inexact() | self.inexact() | x.inexact());
            return Ok(r);
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p
            .max(self.mantissa_max_bit_len())
            .max(x.mantissa_max_bit_len())
            + p_inc;
        if x.is_zero() {
            let mut r = self.gamma(p, rm, cc)?;
            r.set_inexact(r.inexact() | self.inexact() | x.inexact());
            return Ok(r);
        }
        // Γ(s) − γ(s, x) cancels about log2(Γ(s) / Γ(s, x)) bits (≈ x / ln 2 for large x):
        // that loss is measured and added to the working precision (ZF-25), and for x > s the
        // asymptotic series (DLMF 8.11.2, remainder bounded by the first omitted term once
        // n ≥ s − 1, 8.11(ii)) avoids the cancellation altogether.
        // Γ(n) = (n−1)! for integer 1 <= n <= 64 is exact at 320 bits (63! < 2^296). Then the
        // difference's error is only γ's, which can be far below 2^-p_wrk of Γ(s): the
        // claimed accuracy below uses that, so Γ(n, tiny x) just under (n−1)! still rounds.
        let g_exact = self.is_int() && self.cmp(&Self::from_word(64, WORD_BIT_SIZE)?) <= 0;
        let mut extra: usize = 0;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE * 2;
            let (mut ret, acc) = match self.gammainc_upper_asymptotic(x, p_x, cc)? {
                Some(v) => (v, p_wrk),
                None => loop {
                    let pe = round_p(p_x + extra);
                    let g = self.gamma_at(if g_exact { pe.max(320) } else { pe }, cc)?;
                    let lo = self.gammainc_at(x, pe, cc)?;
                    // Keep γ(s, x) ≪ Γ(s) visible: subtract at a width covering the gap, with
                    // a lower part beyond 2^-limit of Γ(s) replaced by that bound (same side
                    // of Γ(s), inside the error budget), so directed modes see the tail.
                    let limit = pe + WORD_BIT_SIZE;
                    let gap = usize::try_from((g.exponent() - lo.exponent()).max(0))
                        .unwrap_or(usize::MAX);
                    // g is read only for its exponent after v is formed.
                    let g_exp = g.exponent();
                    let v = if lo.is_zero() {
                        g
                    } else if gap > limit {
                        let mut t = Self::from_word(1, WORD_BIT_SIZE)?;
                        t.set_exponent(g.exponent() - limit as Exponent + 1);
                        t.set_sign(lo.sign());
                        g.sub(
                            &t,
                            round_p(g.mantissa_max_bit_len() + limit),
                            RoundingMode::None,
                        )?
                    } else {
                        g.sub(
                            &lo,
                            round_p(g.mantissa_max_bit_len() + gap + pe),
                            RoundingMode::None,
                        )?
                    };
                    let top = g_exp.max(lo.exponent());
                    let lost = if v.is_zero() {
                        pe
                    } else {
                        usize::try_from((top - v.exponent()).max(0)).unwrap_or(usize::MAX)
                    };
                    // p_x carries 128 guard bits over p_wrk; up to 56 of them may absorb the loss.
                    if lost <= extra.saturating_add(56) {
                        let acc = if g_exact && lost == 0 {
                            if gap >= WORD_BIT_SIZE {
                                pe
                            } else {
                                pe - WORD_BIT_SIZE
                            }
                        } else {
                            p_wrk
                        };
                        break (v, acc.max(p_wrk));
                    }
                    // Ziv: bump by the measured loss; stop only when bump_prec_retry
                    // exhausts MAX_PREC_RETRY or the new width is unallocatable.
                    extra = lost.saturating_add(WORD_BIT_SIZE);
                },
            };
            if ret.try_set_precision(p, rm, acc)? {
                ret.set_inexact(ret.inexact() | self.inexact() | x.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    /// `Γ(s, x) = x^(s−1) e^(−x) (Σ_{k<n} u_k / x^k + R_n)`, `u_k = (s−1)(s−2)…(s−k)`
    /// (DLMF 8.11.2). For real `s` and `x > 0`, once `n ≥ s − 1` the remainder `R_n` has the
    /// sign of, and is at most, the first omitted term (DLMF 8.11(ii)); an exact zero `u_n`
    /// (integer `s`) ends the sum exactly. `None` when the terms start growing before they
    /// fall below `2^-(p+8)` of the sum, or when `x <= s`.
    fn gammainc_upper_asymptotic(
        &self,
        x: &Self,
        p: usize,
        cc: &mut Consts,
    ) -> Result<Option<Self>, Error> {
        let nr = RoundingMode::None;
        if x.exponent() < 4 || x.cmp(self) <= 0 {
            return Ok(None);
        }
        let one = Self::from_word(1, WORD_BIT_SIZE)?;
        let sm1 = self.sub(&one, p, nr)?;
        let mut term = Self::from_word(1, p)?;
        let mut sum = Self::from_word(1, p)?;
        let mut a = sm1.clone()?;
        let mut done = false;
        let n_max = 4 * p + 64;
        for k in 1..=n_max {
            // term_k = term_{k-1} (s − k) / x; `a` = s − k.
            let next = term.mul(&a, p, nr)?.div(x, p, nr)?;
            let past = Self::from_word(k as Word, WORD_BIT_SIZE)?.cmp(&sm1) >= 0;
            if next.is_zero() {
                done = true;
                break;
            }
            if past && next.exponent().saturating_add(p as Exponent + 8) < sum.exponent() {
                done = true;
                break;
            }
            if past && next.exponent() > term.exponent() {
                break;
            }
            sum = sum.add(&next, p, nr)?;
            term = next;
            a = a.sub(&one, p, nr)?;
        }
        if !done || sum.is_zero() {
            return Ok(None);
        }
        // x^(s−1) e^(−x) = exp((s−1) ln x − x); the argument's absolute error becomes the
        // result's relative error, so carry its exponent in extra bits.
        let mag = x
            .exponent()
            .max(0)
            .saturating_add(WORD_BIT_SIZE as Exponent);
        let pa = round_p(p + usize::try_from(mag).unwrap_or(0));
        let arg = sm1.mul(&x.ln(pa, nr, cc)?, pa, nr)?.sub(x, pa, nr)?;
        let pre = arg.exp(p, nr, cc)?;
        Ok(Some(pre.mul(&sum, p, nr)?))
    }

    fn gammainc_at(&self, x: &Self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        // ZF-37: the old shortcut here returned Γ(s) whenever 2·exponent(x) > p + 8, which is
        // wrong for s > x (γ(2^300, 2^200) ≪ Γ(2^300)); x > s is now taken by the
        // complement in `gammainc`, so the series is always summed.
        // γ(1, x) = 1 − e^{−x} (DLMF 8.4.7); do not send this through the power series.
        let one = Self::from_word(1, p)?;
        if self.cmp(&one) == 0 {
            let nx = x.with_sign(Sign::Neg)?;
            return one.sub(&nx.exp(p, RoundingMode::None, cc)?, p, RoundingMode::None);
        }
        // γ(s,x) = x^s e^{-x} Σ_{k=0} x^k / (s)_{k+1}
        let xs = x.pow(self, p, RoundingMode::None, cc)?;
        let nx = x.with_sign(Sign::Neg)?;
        let exm = nx.exp(p, RoundingMode::None, cc)?;
        let pre = xs.mul(&exm, p, RoundingMode::None)?;
        // s + k and the partial sum: None stands for s itself and for the first term (no copies).
        let mut poch: Option<Self> = None;
        let mut term = Self::from_word(1, p)?.div(self, p, RoundingMode::None)?;
        let mut sum: Option<Self> = None;
        // All terms are positive and the ratio x/(s+k) falls below 1 once k > x − s,
        // so the relative stop is met after about max(0, x − s) + O(√((x+s)p)) terms;
        // while the terms still grow it cannot fire (sum ≤ k·term). Past the peak, a
        // growing term is a stall, not a truncated value (P9: no 2^28 term cap).
        let mut k: usize = 0;
        let mut prev_e = Exponent::MIN;
        loop {
            k = k.checked_add(1).ok_or(Error::PrecisionRetryExhausted)?;
            let pk = poch
                .as_ref()
                .unwrap_or(self)
                .add(&one, p, RoundingMode::None)?;
            let next = term
                .mul(x, p, RoundingMode::None)?
                .div(&pk, p, RoundingMode::None)?;
            let sm = sum
                .as_ref()
                .unwrap_or(&term)
                .add(&next, p, RoundingMode::None)?;
            term = next;
            poch = Some(pk);
            sum = Some(sm);
            if sum
                .as_ref()
                .is_some_and(|sm| series_negligible(&term, sm, p))
            {
                break;
            }
            let te = term.exponent();
            if k > 1 && te > prev_e {
                // Past the Pochhammer peak the terms must shrink. A rise means the
                // series is the wrong kernel for this (s, x); do not spin.
                let past_peak = from_usize_p(k, p)?.cmp(&x.sub(self, p, RoundingMode::None)?) > 0;
                if past_peak {
                    return Err(Error::PrecisionRetryExhausted);
                }
            }
            prev_e = te;
        }
        pre.mul(sum.as_ref().unwrap_or(&term), p, RoundingMode::None)
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
            // ZF-37: near the switch (|x| ≈ 2p/3) the smallest asymptotic term e^-|x| can exceed
            // 2^-p; then the (uncapped) series below is used instead of giving up.
            match self.ei_asymptotic(p, cc) {
                Err(Error::PrecisionRetryExhausted) => {}
                r => return r,
            }
        }
        // The series cancels: for x < 0 its terms reach e^{|x|} while Ei(x) ≈ e^{-|x|}/|x|, and
        // near the root x₀ ≈ 0.3725 the parts γ + ln x and Σ cancel. Measure the loss and
        // re-evaluate with enough guard bits.
        let mut extra = if self.is_negative() && self.exponent() > 0 {
            (3usize << (self.exponent() as u32).min(40)) + 16
        } else {
            16
        };
        loop {
            let (v, top) = self.ei_series(p + extra, cc)?;
            let lost = if v.is_zero() { usize::MAX } else { (top - v.exponent()).max(0) as usize };
            if lost.saturating_add(8) <= extra {
                let mut v = v;
                v.set_precision(p, RoundingMode::None)?;
                return Ok(v);
            }
            let next = lost + 32;
            if next <= extra {
                return Err(Error::PrecisionRetryExhausted);
            }
            extra = next;
        }
    }

    /// `γ + ln|x| + Σ_{n≥1} x^n / (n n!)` and the largest exponent among the parts.
    fn ei_series(&self, p: usize, cc: &mut Consts) -> Result<(Self, Exponent), Error> {
        let g = cc.euler_gamma_num(p, RoundingMode::None)?;
        let lnx = self.abs()?.ln(p, RoundingMode::None, cc)?;
        // x^n/n! and the partial sum both start at x: None stands for `self` (no copy of it).
        let mut term: Option<Self> = None;
        let mut sum: Option<Self> = None;
        let mut top = max_exp(&[&g, &lnx, self]);
        // ZF-37: no term cap and a stop relative to the largest part. The terms x^n/(n·n!)
        // fall monotonically once n > |x| and reach 2^-p of the largest part after about
        // e|x| + p terms; the old absolute stop 2^-p with a cap of ~2p terms exhausted for
        // Ei(105) at p=128 (sum ~2^144). The caller measures cancellation against `top`.
        let ax = self.abs()?;
        let mut n: Word = 1;
        loop {
            n = n.checked_add(1).ok_or(Error::MemoryAllocation)?;
            let nw = Self::from_word(n, p)?;
            let tn = term
                .as_ref()
                .unwrap_or(self)
                .mul(self, p, RoundingMode::None)?
                .div(&nw, p, RoundingMode::None)?;
            let piece = tn.div(&nw, p, RoundingMode::None)?;
            sum = Some(
                sum.as_ref()
                    .unwrap_or(self)
                    .add(&piece, p, RoundingMode::None)?,
            );
            term = Some(tn);
            top = top.max(max_exp(&[&piece]));
            let past = nw.cmp(&ax) > 0;
            if piece.is_zero() || (past && piece.exponent().saturating_add(p as Exponent) < top) {
                break;
            }
        }
        let v = g.add(&lnx, p, RoundingMode::None)?.add(
            sum.as_ref().unwrap_or(self),
            p,
            RoundingMode::None,
        )?;
        Ok((v, top))
    }

    fn ei_asymptotic(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        // e^x / x · Σ_{k=0} k! / x^k (stop when terms grow)
        let ex = self.exp(p, RoundingMode::None, cc)?;
        let pre = ex.div(self, p, RoundingMode::None)?;
        let mut term = Self::from_word(1, p)?;
        let mut s = Self::from_word(1, p)?;
        let mut prev_e = Exponent::MAX;
        // Cap or divergence without meeting the bound: never use the partial sum.
        let mut done = false;
        for k in 1..=(p + 8) {
            term = term.mul(&Self::from_word(k as Word, p)?, p, RoundingMode::None)?;
            term = term.div(self, p, RoundingMode::None)?;
            let e = term.exponent();
            if e > prev_e {
                break;
            }
            prev_e = e;
            s = s.add(&term, p, RoundingMode::None)?;
            if term.is_zero() || e + (p as i64) < 0 {
                done = true;
                break;
            }
        }
        if !done {
            return Err(Error::PrecisionRetryExhausted);
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
        // ZF-37: Si(x) = x(1 − x²/18 + …); below 2^-(p_wrk+1) relative the correction only
        // decides the side (toward zero). Ziv cannot see it (Si(2^-2^40) is x − 2^-3·2^40/18).
        compute_small_exp!(
            self,
            self.exponent().saturating_mul(2).saturating_sub(4),
            true,
            p_wrk,
            p,
            rm
        );
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
        let asym = if x.expint_use_asymptotic(p) {
            // ZF-37: fall back to the series when the asymptotic terms turn before 2^-p.
            match x.si_asymptotic(p, cc) {
                Err(Error::PrecisionRetryExhausted) => None,
                r => Some(r?),
            }
        } else {
            None
        };
        let ret = if let Some(v) = asym {
            v
        } else {
            let pe = p + Self::si_ci_series_guard(x.exponent());
            let mut r = x.si_series(pe)?;
            r.set_precision(p, RoundingMode::None)?;
            r
        };
        if neg {
            ret.neg()
        } else {
            Ok(ret)
        }
    }

    fn si_series(&self, p: usize) -> Result<Self, Error> {
        let x2 = self.mul(self, p, RoundingMode::None)?;
        // x^{2n+1}/(2n+1)! and the partial sum both start at x: None stands for `self` (no copy of it).
        let mut t: Option<Self> = None;
        let mut sum: Option<Self> = None;
        // Cap or divergence without meeting the bound: never use the partial sum.
        let mut n: usize = 1;
        loop {
            let two_n = Self::from_word(
                n.checked_mul(2)
                    .and_then(|v| Word::try_from(v).ok())
                    .ok_or(Error::PrecisionRetryExhausted)?,
                p,
            )?;
            let two_n_1 = Self::from_word(
                n.checked_mul(2)
                    .and_then(|v| v.checked_add(1))
                    .and_then(|v| Word::try_from(v).ok())
                    .ok_or(Error::PrecisionRetryExhausted)?,
                p,
            )?;
            let tn = t
                .as_ref()
                .unwrap_or(self)
                .mul(&x2, p, RoundingMode::None)?
                .div(&two_n, p, RoundingMode::None)?
                .div(&two_n_1, p, RoundingMode::None)?;
            let mut piece = tn.div(&two_n_1, p, RoundingMode::None)?;
            t = Some(tn);
            if n % 2 == 1 {
                piece.set_sign(Sign::Neg);
            }
            let sm = sum
                .as_ref()
                .unwrap_or(self)
                .add(&piece, p, RoundingMode::None)?;
            if series_negligible(&piece, &sm, p) {
                return Ok(sm);
            }
            sum = Some(sm);
            n = n.checked_add(1).ok_or(Error::PrecisionRetryExhausted)?;
        }
    }

    fn si_asymptotic(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        // π/2 − f(x) cos x − g(x) sin x
        let mut half_pi = cc.pi_num(p, RoundingMode::None)?;
        half_pi.div_by_2(RoundingMode::None);
        // ZF-37: 0 < f < 1/x and 0 < g < 1/x² (DLMF 6.12(ii)), so for x ≥ 2^(p+4) the
        // correction is below 2^-(p+2) of π/2 and cos x, sin x (reduction of a huge x, e.g.
        // the largest finite value) are not needed.
        if self.exponent() > p as Exponent + 4 {
            return Ok(half_pi);
        }
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
            // ZF-37: as in `si_at`, the series is the fallback when the asymptotic one turns.
            match self.ci_asymptotic(p, cc) {
                Err(Error::PrecisionRetryExhausted) => {}
                r => return r,
            }
        }
        {
            let pe = p + Self::si_ci_series_guard(self.exponent());
            let mut r = self.ci_series(pe, cc)?;
            r.set_precision(p, RoundingMode::None)?;
            Ok(r)
        }
    }

    /// Guard bits for the Si/Ci power series (ZF-24): its alternating terms peak near
    /// `e^x / √x` while the result is O(1), so about `x log2 e < 1.5 · 2^exponent(x)` bits
    /// cancel. The series is only used below `2^20`.
    fn si_ci_series_guard(e: Exponent) -> usize {
        let cancel = if e > 0 { (3usize << (e as u32).min(40)) / 2 } else { 0 };
        round_p(cancel + 2 * WORD_BIT_SIZE)
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
        // Cap or divergence without meeting the bound: never use the partial sum.
        let mut n: usize = 2;
        loop {
            let a = Self::from_word(
                n.checked_mul(2)
                    .and_then(|v| v.checked_sub(1))
                    .and_then(|v| Word::try_from(v).ok())
                    .ok_or(Error::PrecisionRetryExhausted)?,
                p,
            )?;
            let b = Self::from_word(
                n.checked_mul(2)
                    .and_then(|v| Word::try_from(v).ok())
                    .ok_or(Error::PrecisionRetryExhausted)?,
                p,
            )?;
            u = u.mul(&x2, p, RoundingMode::None)?;
            u = u.div(&a, p, RoundingMode::None)?;
            u = u.div(&b, p, RoundingMode::None)?;
            u = u.neg()?;
            let piece = u.div(&b, p, RoundingMode::None)?;
            sum = sum.add(&piece, p, RoundingMode::None)?;
            if piece.is_zero() || piece.exponent() + (p as i64) < 0 {
                break;
            }
            n = n.checked_add(1).ok_or(Error::PrecisionRetryExhausted)?;
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
        // The k = 0 terms are 1/x and 1/x². g_term None stands for invx2 itself, and the sums f, g are None until the
        // first k > 0 term is added (they equal the k = 0 terms then): no copies.
        let mut f_term = invx;
        let mut g_term: Option<Self> = None;
        let mut f: Option<Self> = None;
        let mut g: Option<Self> = None;
        let mut prev_e = f_term.exponent();
        // Cap or divergence without meeting the bound: never use the partial sum.
        for k in 0..=(p + 8) {
            let a = Self::from_word((2 * k + 1) as Word, p)?;
            let b = Self::from_word((2 * k + 2) as Word, p)?;
            let c = Self::from_word((2 * k + 3) as Word, p)?;
            let ft = f_term
                .mul(&a, p, RoundingMode::None)?
                .mul(&b, p, RoundingMode::None)?
                .mul(&invx2, p, RoundingMode::None)?
                .neg()?;
            let gt = g_term
                .as_ref()
                .unwrap_or(&invx2)
                .mul(&b, p, RoundingMode::None)?
                .mul(&c, p, RoundingMode::None)?
                .mul(&invx2, p, RoundingMode::None)?
                .neg()?;
            let e = ft.exponent();
            if e > prev_e {
                break;
            }
            prev_e = e;
            let fs = f
                .as_ref()
                .unwrap_or(&f_term)
                .add(&ft, p, RoundingMode::None)?;
            let gs = g
                .as_ref()
                .unwrap_or(&invx2)
                .add(&gt, p, RoundingMode::None)?;
            if ft.is_zero() || e + (p as i64) < 0 {
                return Ok((fs, gs));
            }
            f = Some(fs);
            g = Some(gs);
            f_term = ft;
            g_term = Some(gt);
        }
        Err(Error::PrecisionRetryExhausted)
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
        let xmin: Exponent = 1 << (e - 1);
        xmin >= ((p as Exponent) * 2) / 3
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
        if !sine {
            // ZF-37: C(x) = x(1 − (π²/40)x⁴ + …), π²/40 < 1/4: the same side-only case as Si.
            compute_small_exp!(
                self,
                self.exponent().saturating_mul(4).saturating_sub(2),
                true,
                p_wrk,
                p,
                rm
            );
        }
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
            // The series terms peak near e^{πx²/2} while the result is O(1): add
            // log2(e)·πx²/2 ≈ 2.27x² guard bits (bounded via x² < 2^e).
            let x2 = x.mul(&x, WORD_BIT_SIZE, RoundingMode::Up)?;
            let e2 = x2.exponent();
            let extra = if e2 > 0 { (5usize << (e2 as u32).min(40)) / 2 } else { 0 } + 32;
            let mut r = x.fresnel_sc_series(sine, p + extra, cc)?;
            r.set_precision(p, RoundingMode::None)?;
            r
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
        // The asymptotic series' smallest term is about e^{-πx²/2}, i.e. 2.27x² bits.
        let xmin = 1i64 << (e - 1);
        xmin * xmin * 9 / 4 > p as i64 + 16
    }

    fn fresnel_sc_asymptotic(&self, sine: bool, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        // S = 1/2 − f cos(πx²/2) − g sin(πx²/2)
        // C = 1/2 + f sin(πx²/2) − g cos(πx²/2)
        let (f, g) = self.fresnel_aux_fg(p, cc)?;
        // sin/cos(πx²/2) has period 4 in x², and x² is an exact dyadic: reduce it exactly mod 4
        // first (ZF-32). Rounding x² or πx²/2 to p bits leaves an absolute error of
        // 2^(2·exponent(x) − p) in the angle, i.e. no correct bits once x² ≳ 2^p.
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let two = Self::from_word(2, p)?;
        let pi2 = pi.div(&two, p, RoundingMode::None)?;
        let w2 = round_p(2 * self.mantissa_max_bit_len() + WORD_BIT_SIZE);
        let x2 = self.mul(self, w2, RoundingMode::None)?;
        let t = x2.rem(&Self::from_word(4, WORD_BIT_SIZE)?)?;
        let arg = pi2.mul(&t, p, RoundingMode::None)?;
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
        // The sums start at the n = 0 terms: f_sum is a fresh 1, g_sum is None (it equals g_term) until the first n > 0
        // term is added, so g_term is never copied.
        let mut f_sum = Self::from_word(1, p)?;
        let mut g_sum: Option<Self> = None;
        let mut prev_e = f_term.exponent();
        // Cap or divergence without meeting the bound: never use the partial sum.
        for n in 0..=(p + 8) {
            let a = Self::from_word((4 * n + 1) as Word, p)?;
            let b = Self::from_word((4 * n + 3) as Word, p)?;
            let c = Self::from_word((4 * n + 5) as Word, p)?;
            f_term = f_term
                .mul(&a, p, RoundingMode::None)?
                .mul(&b, p, RoundingMode::None)?
                .mul(&v2, p, RoundingMode::None)?
                .neg()?;
            let gt = g_term
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
            let gs = g_sum
                .as_ref()
                .unwrap_or(&g_term)
                .add(&gt, p, RoundingMode::None)?;
            g_term = gt;
            if f_term.is_zero() || e + (p as i64) < 0 {
                return Ok((
                    inv_pix.mul(&f_sum, p, RoundingMode::None)?,
                    inv_pix.mul(&gs, p, RoundingMode::None)?,
                ));
            }
            g_sum = Some(gs);
        }
        Err(Error::PrecisionRetryExhausted)
    }

    fn fresnel_sc_series(&self, sine: bool, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let two = Self::from_word(2, p)?;
        let pi2 = pi.div(&two, p, RoundingMode::None)?;
        let x2 = self.mul(self, p, RoundingMode::None)?;
        let x4 = x2.mul(&x2, p, RoundingMode::None)?;
        // x^{4n+1} or x^{4n+3}, and (π/2)^{2n} or (π/2)^{2n+1}: None stands for x itself and for π/2 itself (no copies).
        let mut pow_x: Option<Self> =
            if sine { Some(self.mul(&x2, p, RoundingMode::None)?) } else { None };
        let mut fact = Self::from_word(1, p)?;
        let mut pi2_pow: Option<Self> = if sine { None } else { Some(Self::from_word(1, p)?) };
        let mut sign_pos = true;
        // None until the first piece, which becomes the sum itself.
        let mut sum: Option<Self> = None;
        // Cap or divergence without meeting the bound: never use the partial sum.
        let mut n: usize = 0;
        loop {
            let den_i = if sine {
                n.checked_mul(4).and_then(|v| v.checked_add(3))
            } else {
                n.checked_mul(4).and_then(|v| v.checked_add(1))
            }
            .and_then(|v| Word::try_from(v).ok())
            .ok_or(Error::PrecisionRetryExhausted)?;
            let den = Self::from_word(den_i, p)?;
            let mut piece = pi2_pow.as_ref().unwrap_or(&pi2).mul(
                pow_x.as_ref().unwrap_or(self),
                p,
                RoundingMode::None,
            )?;
            piece = piece.div(&fact, p, RoundingMode::None)?;
            piece = piece.div(&den, p, RoundingMode::None)?;
            if !sign_pos {
                piece.set_sign(Sign::Neg);
            }
            let (sm, small) = match sum.take() {
                None => {
                    let small = series_negligible(&piece, &piece, p);
                    (piece, small)
                }
                Some(s) => {
                    let sm = s.add(&piece, p, RoundingMode::None)?;
                    let small = series_negligible(&piece, &sm, p);
                    (sm, small)
                }
            };
            if small {
                return Ok(sm);
            }
            sum = Some(sm);
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
            pi2_pow = Some(
                pi2_pow
                    .as_ref()
                    .unwrap_or(&pi2)
                    .mul(&pi2, p, RoundingMode::None)?
                    .mul(&pi2, p, RoundingMode::None)?,
            );
            pow_x = Some(
                pow_x
                    .as_ref()
                    .unwrap_or(self)
                    .mul(&x4, p, RoundingMode::None)?,
            );
            sign_pos = !sign_pos;
            n = n.checked_add(1).ok_or(Error::PrecisionRetryExhausted)?;
        }
    }

    /// Bessel function of the first kind `J_n(self)` for integer order `n`.
    pub fn bessel_j(
        &self,
        n: usize,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;

        if self.is_zero() {
            let mut z = if n == 0 {
                Self::from_word(1, p)?
            } else if n % 2 == 1 {
                Self::new2(p, self.sign(), false)?
            } else {
                Self::new2(p, Sign::Pos, false)?
            };
            z.set_inexact(self.inexact());
            return Ok(z);
        }

        // Tiny |x|: J0 = 1 - x²/4 + …, J1 = x/2 - x³/16 + …; the correction is far below
        // half an ulp, so round the exact leading term once toward it.
        let tiny_exp = -Exponent::try_from(p / 2 + 4).map_err(|_| Error::InvalidArgument)?;
        if n == 1 && self.exponent() == EXPONENT_MIN {
            // ZF-11: x/2 lies on the subnormal grid or exactly midway between two of its points
            // (one bit is dropped), and the series would widen by the ~2^62-bit exponent.
            return Self::bessel_j1_at_min_exponent(self, p, rm);
        }
        if n <= 1 && self.exponent() < tiny_exp && (n == 0 || self.exponent() > EXPONENT_MIN) {
            let mut r = if n == 0 {
                Self::round_offset(&Self::from_word(1, p)?, Sign::Neg, p, rm)?
            } else {
                // x / 2, exact (x is normal here): the product by one half keeps every bit
                let mut one_half = Self::from_word(1, WORD_BIT_SIZE)?;
                one_half.set_exponent(0);
                let half = self.mul(&one_half, self.mantissa_max_bit_len(), RoundingMode::None)?;
                Self::round_offset(&half, self.sign().invert(), p, rm)?
            };
            r.set_inexact(true);
            return Ok(r);
        }

        let nu_ord = from_usize_p(n, p)?;
        // ZF-37: a leading term (x/2)^n/n! outside the range, or a correction below the precision.
        if let Some(v) = bessel_ji_scaled(self, &nu_ord, p, rm, cc, false)? {
            return Ok(v);
        }
        if use_bessel_hankel_real(self, &nu_ord, p) {
            // Hankel's accuracy is capped by its smallest term (≈ e^{-2|x|}); once the Ziv
            // precision passes it the expansion reports exhaustion and the series takes over.
            match self.bessel_j_hankel_signed(n, p, rm, cc) {
                Err(Error::PrecisionRetryExhausted) => {}
                other => return other,
            }
        }

        let mut p_inc = WORD_BIT_SIZE;
        let extra = extra_bits_for_degree(n as u32);
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc + extra;
        // The power series cancels: terms reach ≈ e^{|x|}
        // while |J_n(x)| ≈ (π|x|/2)^{-1/2}.
        let guard = bessel_cancel_bits(self, 3, 2)?;
        let mut lost_done = 0;

        loop {
            let p_x = p_wrk + WORD_BIT_SIZE + extra + guard;
            // ZF-63: the power series for every n. Miller's recurrence started at n + p/2 + 8
            // whatever x was, below x for x > 4n + 2p + 256 (J_100(1000) at 64 bits), where the
            // backward recurrence has not settled and the normalized value is wrong.
            let mut ret = self.bessel_j_series(n, p_x)?;
            // Near a zero of J_n the series sum cancels: its error is absolute (≈ 2^-p_x),
            // so a result of exponent -k has k fewer correct bits. Widen once by k.
            let lost = if ret.exponent() < 0 {
                usize::try_from(-ret.exponent()).map_err(|_| Error::InvalidArgument)?
            } else {
                0
            };
            if lost > lost_done {
                p_wrk = round_p(p_wrk + lost - lost_done);
                lost_done = lost;
                continue;
            }
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    /// `J_1(x) = x/2 − x³/16 + …` for `x` at the minimum exponent: the result is `x/2 − δ`
    /// (magnitude) with `δ` far below the subnormal spacing, rounded once in `rm`.
    fn bessel_j1_at_min_exponent(x: &Self, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        let neg = x.is_negative();
        let mut xa = x.abs()?;
        xa.set_inexact(false);
        let mut half = Self::from_word(1, WORD_BIT_SIZE)?;
        half.set_exponent(0);
        let h = xa.mul(&half, p, RoundingMode::ToZero)?;
        let two = Self::from_word(2, WORD_BIT_SIZE)?;
        let on_grid = h
            .mul(&two, p.max(xa.mantissa_max_bit_len()), RoundingMode::None)?
            .cmp(&xa)
            == 0;
        let away = match rm {
            RoundingMode::FromZero => Some(true),
            // `None` drops the excess bits (as `round_offset` does for `J_0`): toward zero.
            RoundingMode::ToZero | RoundingMode::None => Some(false),
            RoundingMode::Up => Some(!neg),
            RoundingMode::Down => Some(neg),
            _ => None,
        };
        let mut r = match (away, on_grid) {
            // Toward zero from an on-grid x/2: the next value below it.
            (Some(false), true) => h.next_after(&Self::new(p)?, p, RoundingMode::None)?,
            // Away from zero from a midpoint: the grid point above.
            (Some(true), false) => xa.mul(&half, p, RoundingMode::FromZero)?,
            // On grid and away, or toward zero from a midpoint, or nearest (a midpoint minus δ
            // rounds down; an on-grid x/2 minus δ rounds to x/2).
            _ => h,
        };
        r.set_sign(if neg { Sign::Neg } else { Sign::Pos });
        r.set_inexact(true);
        Ok(r)
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

        // Cap or divergence without meeting the bound: never use the partial sum.
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
            if series_negligible(&t, &sum, p) {
                break;
            }
        }
        Ok(sum)
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
            // J_{−n} = (−1)^n J_n: when the sign flips, round J_n in the mirrored direction
            // so that −J_n is rounded in `rm` (Up of −v is −(Down of v)).
            if n >= 0 || n % 2 == 0 {
                self.bessel_j(an, p, rm, cc)
            } else {
                self.bessel_j(an, p, crate::common::util::invert_rm_for_sign(rm), cc)?
                    .neg()
            }
        } else if self.is_zero() {
            bessel_ji_at_zero(self, nu, p, false)
        } else if self.is_negative() && nu.is_int() {
            // ZF-37: an integer order past i32: J_ν(−x) = (−1)^ν J_ν(x).
            bessel_int_order_neg_x(self, nu, p, rm, cc, ExactNumNumber::bessel_j_nu)
        } else if !self.is_positive() {
            Err(Error::InvalidArgument)
        } else if let Some(v) = bessel_ji_scaled(self, nu, p, rm, cc, false)? {
            Ok(v)
        } else if use_bessel_hankel_real(self, nu, p) {
            self.bessel_jy_hankel(nu, p, rm, cc, true)
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
        if self.is_zero() {
            return bessel_y_at_zero(nu, p);
        }
        if !self.is_positive() {
            return Err(Error::InvalidArgument);
        }
        if let Some(v) = bessel_yk_scaled(self, nu, p, rm, cc, false)? {
            return Ok(v);
        }
        if use_bessel_hankel_real(self, nu, p) {
            return self.bessel_jy_hankel(nu, p, rm, cc, false);
        }
        // Integer |ν| can fail hankel_order_ok (|ν|² ≥ 2|x|) while |x| is still
        // in the Hankel range for Y_0 and Y_1. Form those by Hankel and recur
        // up (Y is dominant in ν, and ν < |x| in the cases this reaches).
        if hankel_mag_ok(self, p) {
            if let Some(n) = as_i32_exact(nu, p)? {
                return self.bessel_y_int_via_hankel01(n, p, rm, cc);
            }
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p
            .max(self.mantissa_max_bit_len())
            .max(nu.mantissa_max_bit_len())
            + p_inc;
        let guard = bessel_cancel_bits(self, 3, 2)?;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE + guard;
            let mut ret = self.bessel_y_at(nu, p_x, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact() | nu.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    /// \(Y_n\) from Hankel \(Y_0,Y_1\) and the forward recurrence, under Ziv.
    fn bessel_y_int_via_hankel01(
        &self,
        n: i32,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE;
            let mut ret = self.bessel_y_int_hankel_at(n, p_x, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn bessel_y_int_hankel_at(&self, n: i32, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let an = n.unsigned_abs();
        let zero = Self::from_word(0, p)?;
        let y0 = self.bessel_hankel_y(&zero, p, cc)?;
        let y = if an == 0 {
            y0
        } else {
            let one = Self::from_word(1, p)?;
            let y1 = self.bessel_hankel_y(&one, p, cc)?;
            if an == 1 {
                y1
            } else {
                let mut ym2 = y0;
                let mut ym1 = y1;
                for m in 1..an {
                    let two_m = Self::from_word((2 * m) as Word, p)?;
                    let ym = two_m
                        .div(self, p, RoundingMode::None)?
                        .mul(&ym1, p, RoundingMode::None)?
                        .sub(&ym2, p, RoundingMode::None)?;
                    ym2 = ym1;
                    ym1 = ym;
                }
                ym1
            }
        };
        if n >= 0 || n % 2 == 0 {
            Ok(y)
        } else {
            y.neg()
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
            if s.is_zero() || s.exponent() + (p as i64) < 0 {
                return Err(Error::InvalidArgument);
            }
            let jp = self.bessel_jn_series(nu, p, cc, true)?;
            let nneg = nu.with_sign(if nu.is_positive() { Sign::Neg } else { Sign::Pos })?;
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
            return bessel_ji_at_zero(self, nu, p, true);
        }
        // ZF-37: integer orders take any sign of x; non-integer ones only x > 0 (below).
        if nu.is_int() || self.is_positive() {
            if let Some(v) = bessel_ji_scaled(self, nu, p, rm, cc, true)? {
                return Ok(v);
            }
        }
        if let Some(n) = as_i32_exact(nu, p)? {
            let an = Self::from_word(n.unsigned_abs() as Word, p)?;
            let ax = self.with_sign(Sign::Pos)?;
            // I_n(−x) = (−1)^n I_n(x): a flipped sign rounds I_n(|x|) in the mirrored direction.
            if self.is_negative() && n.unsigned_abs() % 2 == 1 {
                ax.bessel_series(
                    &an,
                    p,
                    crate::common::util::invert_rm_for_sign(rm),
                    cc,
                    false,
                )?
                .neg()
            } else {
                ax.bessel_series(&an, p, rm, cc, false)
            }
        } else if self.is_negative() && nu.is_int() {
            // ZF-37: an integer order past i32: I_ν(−x) = (−1)^ν I_ν(x).
            bessel_int_order_neg_x(self, nu, p, rm, cc, ExactNumNumber::bessel_i)
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
        if self.is_zero() {
            // K_ν(x) → +∞ as x → 0 for every real ν: a pole (C99 7.12.1), +Inf like Y_ν(0).
            return Err(Error::ExponentOverflow(Sign::Pos));
        }
        if !self.is_positive() {
            return Err(Error::InvalidArgument);
        }
        if let Some(v) = bessel_yk_scaled(self, nu, p, rm, cc, true)? {
            return Ok(v);
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
        let nu = nu.with_sign(Sign::Pos)?;
        if let Some(n) = half_integer_n(&nu, p)? {
            return self.k_half_integer(n, p, cc);
        }
        // The asymptotic series' smallest term is ≈ e^{-2x} (2.885x bits), so it is used only
        // when that covers `p`. Below that the I-based formulas cancel: I ≈ e^{x}, K ≈ e^{-x},
        // so 2.885x guard bits are added. If that width is unallocatable, `|x|` is
        // far past the asymptotic threshold: use the asymptotic series instead.
        let two_x_bits = match bessel_cancel_bits(self, 3, 1) {
            Ok(bits) => bits,
            Err(Error::MemoryAllocation) => return self.k_asymptotic(&nu, p, cc),
            Err(e) => return Err(e),
        };
        if two_x_bits >= p + 16 && self.exponent() >= 3 {
            // ZF-37: the Hankel series turns early when ν is comparable to x (its terms grow
            // like (4ν²)/(8xk) until k ≈ ν²/(2x)); then the I-based formula below, with its
            // 2.885x cancellation guard, gives the value instead of PrecisionRetryExhausted
            // (K_100(100) at p = 64).
            match self.k_asymptotic(&nu, p, cc) {
                Err(Error::PrecisionRetryExhausted) => {}
                r => return r,
            }
        }
        let pe = p + two_x_bits;
        let v = if let Some(n) = as_i32_exact(&nu, pe)? {
            if n < 0 {
                return Err(Error::InvalidArgument);
            }
            self.k_integer(n as u32, pe, cc)?
        } else {
            self.k_real(&nu, pe, cc)?
        };
        let mut v = v;
        v.set_precision(p, RoundingMode::None)?;
        Ok(v)
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
        let guard = if alternating { bessel_cancel_bits(self, 3, 2)? } else { 0 };
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE + guard;
            let (mut ret, peak, terms) = self.bessel_jn_series_peak(nu, p_x, cc, alternating)?;
            // Near a zero of J_ν the sum is far below its largest term. Term k comes from
            // term k−1 by a few roundings, so its relative error is at most (4k+4) ulps at
            // p_x, and the T additions add T ulps of the peak: the absolute error is below
            // (4T² + 5T)·2^(peak − p_x) < 2^(peak − p_x + 2·bitlen(T) + 3). The sum is good
            // to p_x − (peak − e(sum)) − 2·bitlen(T) − 5 bits. Claim no more than that.
            let claim = if ret.is_zero() {
                0
            } else {
                let loss = peak.saturating_sub(ret.exponent()).max(0) as usize;
                let acc = 2 * (usize::BITS - terms.leading_zeros()) as usize + 5;
                // Word-aligned, as `try_set_precision` requires.
                p_x.saturating_sub(loss + acc).min(p_wrk) / WORD_BIT_SIZE * WORD_BIT_SIZE
            };
            if claim > p && ret.try_set_precision(p, rm, claim)? {
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
        Ok(self.bessel_jn_series_peak(nu, p, cc, alternating)?.0)
    }

    /// [`Self::bessel_jn_series`] with the largest term exponent and the number of terms,
    /// for the cancellation-aware claim in [`Self::bessel_series`].
    fn bessel_jn_series_peak(
        &self,
        nu: &Self,
        p: usize,
        cc: &mut Consts,
        alternating: bool,
    ) -> Result<(Self, Exponent, usize), Error> {
        let two = Self::from_word(2, p)?;
        let half = self.div(&two, p, RoundingMode::None)?;
        let one = Self::from_word(1, p)?;
        let g0 = nu.add(&one, p, RoundingMode::None)?.gamma_at(p, cc)?;
        let mut term = half
            .pow(nu, p, RoundingMode::None, cc)?
            .div(&g0, p, RoundingMode::None)?;
        // None until the first k > 0 term is added: the sum is the k = 0 term itself (no copy).
        let mut sum: Option<Self> = None;
        let mut peak = if term.is_zero() { Exponent::MIN } else { term.exponent() };
        let mut terms = 1usize;
        let hh = half.mul(&half, p, RoundingMode::None)?;
        // Terms grow until k ≈ |x|/2 and then decay. Run until the tail is
        // below 2^{-p}; never use a partial sum (the old series_n_max, even
        // with a 2^20 bump, stopped I_ν short of the peak for large |x|).
        let mut k: usize = 1;
        loop {
            let kk = Self::from_word(
                Word::try_from(k).map_err(|_| Error::PrecisionRetryExhausted)?,
                p,
            )?;
            let den = kk
                .add(nu, p, RoundingMode::None)?
                .mul(&kk, p, RoundingMode::None)?;
            let mut tn = term
                .mul(&hh, p, RoundingMode::None)?
                .div(&den, p, RoundingMode::None)?;
            if alternating {
                tn = tn.neg()?;
            }
            let sm = sum
                .as_ref()
                .unwrap_or(&term)
                .add(&tn, p, RoundingMode::None)?;
            term = tn;
            terms += 1;
            if !term.is_zero() {
                peak = peak.max(term.exponent());
            }
            if series_negligible(&term, &sm, p) {
                return Ok((sm, peak, terms));
            }
            sum = Some(sm);
            k = k.checked_add(1).ok_or(Error::PrecisionRetryExhausted)?;
        }
    }

    /// Large-|x| Hankel expansion of `J_n` on `|self|`, then `J_n(-x)=(-1)^n J_n(x)`.
    fn bessel_j_hankel_signed(
        &self,
        n: usize,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let mut p_inc = WORD_BIT_SIZE;
        let extra = extra_bits_for_degree(n.min(u32::MAX as usize) as u32);
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc + extra;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE + extra;
            let nu = from_usize_p(n, p_x)?;
            let ax = self.abs()?;
            let mut ret = ax.bessel_hankel_j(&nu, p_x, cc)?;
            if self.is_negative() && n % 2 == 1 {
                ret = ret.neg()?;
            }
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn bessel_jy_hankel(
        &self,
        nu: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
        first_kind: bool,
    ) -> Result<Self, Error> {
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p
            .max(self.mantissa_max_bit_len())
            .max(nu.mantissa_max_bit_len())
            + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE;
            let mut ret = if first_kind {
                self.bessel_hankel_j(nu, p_x, cc)?
            } else {
                self.bessel_hankel_y(nu, p_x, cc)?
            };
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact() | nu.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    /// χ = x − (2ν+1)π/4 (DLMF 10.17.1).
    fn bessel_hankel_chi(&self, nu: &Self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let two = Self::from_word(2, p)?;
        let four = Self::from_word(4, p)?;
        let one = Self::from_word(1, p)?;
        let two_nu_1 = two
            .mul(nu, p, RoundingMode::None)?
            .add(&one, p, RoundingMode::None)?;
        let phase = two_nu_1
            .mul(&pi, p, RoundingMode::None)?
            .div(&four, p, RoundingMode::None)?;
        self.sub(&phase, p, RoundingMode::None)
    }

    /// ω = √(2/(πx)).
    fn bessel_hankel_omega(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let two = Self::from_word(2, p)?;
        two.div(&pi.mul(self, p, RoundingMode::None)?, p, RoundingMode::None)?
            .sqrt(p, RoundingMode::None)
    }

    /// Hankel amplitude series P, Q (DLMF 10.17.3–4), stopped at the smallest term.
    fn bessel_hankel_pq(&self, nu: &Self, p: usize) -> Result<(Self, Self), Error> {
        let two = Self::from_word(2, p)?;
        let eight = Self::from_word(8, p)?;
        let two_nu = two.mul(nu, p, RoundingMode::None)?;
        let mu = two_nu.mul(&two_nu, p, RoundingMode::None)?;
        let eight_x = eight.mul(self, p, RoundingMode::None)?;
        let mut prod = Self::from_word(1, p)?;
        let mut kf = Self::from_word(1, p)?;
        let mut pz = Self::from_word(1, p)?;
        let mut psum = Self::from_word(1, p)?;
        let mut qsum = Self::from_word(0, p)?;
        let mut prev_e = Exponent::MAX;
        let k_max = (p + 64).min((Word::MAX as usize - 1) / 2);
        // Cap or divergence without meeting the bound: never use the partial sum.
        let mut done = false;
        for k in 1..=k_max {
            let odd = Self::from_word((2 * k - 1) as Word, p)?;
            let odd2 = odd.mul(&odd, p, RoundingMode::None)?;
            prod = prod.mul(
                &mu.sub(&odd2, p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )?;
            kf = kf.mul(&Self::from_word(k as Word, p)?, p, RoundingMode::None)?;
            pz = pz.mul(&eight_x, p, RoundingMode::None)?;
            let term = prod.div(&kf.mul(&pz, p, RoundingMode::None)?, p, RoundingMode::None)?;
            if k % 2 == 0 {
                let negated;
                let signed = if (k / 2) % 2 == 1 {
                    negated = term.neg()?;
                    &negated
                } else {
                    &term
                };
                psum = psum.add(signed, p, RoundingMode::None)?;
            } else {
                let negated;
                let signed = if ((k - 1) / 2) % 2 == 1 {
                    negated = term.neg()?;
                    &negated
                } else {
                    &term
                };
                qsum = qsum.add(signed, p, RoundingMode::None)?;
            }
            let e = term.exponent();
            if term.is_zero() || e + (p as i64) < 0 {
                done = true;
                break;
            }
            if k > 2 && e > prev_e {
                break;
            }
            prev_e = e;
        }
        if !done {
            return Err(Error::PrecisionRetryExhausted);
        }
        Ok((psum, qsum))
    }

    /// `J_ν(self)` for `self > 0` via the Hankel expansion (DLMF 10.17.5).
    fn bessel_hankel_j(&self, nu: &Self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let chi = self.bessel_hankel_chi(nu, p, cc)?;
        let omega = self.bessel_hankel_omega(p, cc)?;
        let (pp, qq) = self.bessel_hankel_pq(nu, p)?;
        let c = chi.cos(p, RoundingMode::None, cc)?;
        let s = chi.sin(p, RoundingMode::None, cc)?;
        omega.mul(
            &pp.mul(&c, p, RoundingMode::None)?.sub(
                &qq.mul(&s, p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )?,
            p,
            RoundingMode::None,
        )
    }

    /// `Y_ν(self)` for `self > 0` via the Hankel expansion (DLMF 10.17.6).
    fn bessel_hankel_y(&self, nu: &Self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let chi = self.bessel_hankel_chi(nu, p, cc)?;
        let omega = self.bessel_hankel_omega(p, cc)?;
        let (pp, qq) = self.bessel_hankel_pq(nu, p)?;
        let c = chi.cos(p, RoundingMode::None, cc)?;
        let s = chi.sin(p, RoundingMode::None, cc)?;
        omega.mul(
            &pp.mul(&s, p, RoundingMode::None)?.add(
                &qq.mul(&c, p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )?,
            p,
            RoundingMode::None,
        )
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
        // Cap or divergence without meeting the bound: never use the partial sum.
        let mut m: usize = 1;
        loop {
            fact = fact.mul(
                &Self::from_word(
                    Word::try_from(m).map_err(|_| Error::PrecisionRetryExhausted)?,
                    p,
                )?,
                p,
                RoundingMode::None,
            )?;
            zk = zk.mul(&z, p, RoundingMode::None)?;
            let h = harmonic_u(m, p)?;
            let mut term = h.div(
                &fact.mul(&fact, p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )?;
            term = term.mul(&zk, p, RoundingMode::None)?;
            if m.is_multiple_of(2) {
                term = term.neg()?;
            }
            sum = sum.add(&term, p, RoundingMode::None)?;
            if term.is_zero() || term.exponent() + (p as i64) < 0 {
                break;
            }
            m = m.checked_add(1).ok_or(Error::PrecisionRetryExhausted)?;
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
        // Cap or divergence without meeting the bound: never use the partial sum.
        let mut k: usize = 0;
        loop {
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
            if (term.is_zero() || term.exponent() + (p as i64) < 0) && k > 0 {
                break;
            }
            let kp = k.checked_add(1).ok_or(Error::PrecisionRetryExhausted)?;
            kfact = kfact.mul(
                &Self::from_word(
                    Word::try_from(kp).map_err(|_| Error::PrecisionRetryExhausted)?,
                    p,
                )?,
                p,
                RoundingMode::None,
            )?;
            kp1fact = kp1fact.mul(
                &Self::from_word(
                    Word::try_from(kp.checked_add(1).ok_or(Error::PrecisionRetryExhausted)?)
                        .map_err(|_| Error::PrecisionRetryExhausted)?,
                    p,
                )?,
                p,
                RoundingMode::None,
            )?;
            zk = zk.mul(&z, p, RoundingMode::None)?;
            k = kp;
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
        let nx = self.with_sign(Sign::Neg)?;
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
        let mut psi = one.digamma_at(p, cc)?;
        // None until the first k > 0 term is added: the sum is ψ(1) itself, held in psi until then (no copy).
        let mut sum: Option<Self> = None;
        // Cap or divergence without meeting the bound: never use the partial sum.
        let mut k: usize = 1;
        loop {
            let kw = Word::try_from(k).map_err(|_| Error::PrecisionRetryExhausted)?;
            fact = fact.mul(&Self::from_word(kw, p)?, p, RoundingMode::None)?;
            zk = zk.mul(&z, p, RoundingMode::None)?;
            // ψ(k+1) = ψ(k) + 1/k (a fresh digamma per term made K_n very slow).
            let psi_k = psi.add(
                &Self::from_word(kw, p)?.reciprocal(p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )?;
            let term = psi_k
                .div(
                    &fact.mul(&fact, p, RoundingMode::None)?,
                    p,
                    RoundingMode::None,
                )?
                .mul(&zk, p, RoundingMode::None)?;
            // while sum is None, psi is still ψ(1), the sum's first part
            sum = Some(
                sum.as_ref()
                    .unwrap_or(&psi)
                    .add(&term, p, RoundingMode::None)?,
            );
            psi = psi_k;
            if term.is_zero() || term.exponent() + (p as i64) < 0 {
                break;
            }
            k = k.checked_add(1).ok_or(Error::PrecisionRetryExhausted)?;
        }
        half.ln(p, RoundingMode::None, cc)?
            .neg()?
            .mul(&i0, p, RoundingMode::None)?
            .add(sum.as_ref().unwrap_or(&psi), p, RoundingMode::None)
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
        if s.is_zero() || s.exponent() + (p as i64) < 0 {
            return Err(Error::InvalidArgument);
        }
        let ip = self.bessel_i(mu, p, RoundingMode::None, cc)?;
        let mneg = mu.with_sign(Sign::Neg)?;
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
        let mut prev_e = Exponent::MAX;
        // The smallest term is near k ≈ 2x; stop there or at 2^-p.
        // Cap or divergence without meeting the bound: never use the partial sum.
        let mut done = false;
        for k in 1..=(p + 64) {
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
            if term.is_zero() || e + (p as i64) < 0 {
                done = true;
                break;
            }
            prev_e = e;
        }
        if !done {
            return Err(Error::PrecisionRetryExhausted);
        }
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let two = Self::from_word(2, p)?;
        let nx = self.with_sign(Sign::Neg)?;
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

    pub(crate) fn elliptic_k_at(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
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
        // F(x|0) = arcsin(x). One call rounds to (p, rm). asin at a wider
        // width with RoundingMode::None, then a second rounding, moves a
        // directed mode by an ulp or fails the roundability test.
        if m.is_zero() {
            let mut r = self.asin(p, rm, cc)?;
            r.set_inexact(r.inexact() | self.inexact() | m.inexact());
            return Ok(r);
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
        Ok(self.elliptic_f_sc(None, m, p, cc)?.0)
    }

    /// Bits a cancelling sum `t = t1 ± t2` loses: \(\max(e_1,e_2)-e_t\), never negative.
    fn cancel_bits(t1: &Self, t2: &Self, t: &Self) -> usize {
        if t.is_zero() {
            return usize::MAX / 4;
        }
        let top = if t1.is_zero() {
            t2.exponent()
        } else if t2.is_zero() {
            t1.exponent()
        } else {
            t1.exponent().max(t2.exponent())
        };
        (top - t.exponent()).max(0) as usize
    }

    /// \((x^2,\ 1-x^2,\ 1-m x^2)\) for \(x=|\sin\varphi|\), each to a few ulps relative.
    ///
    /// `c` is \(\cos\varphi\) when the caller has it. \(1-x^2\) formed from a rounded
    /// sine keeps only the absolute error of the sine, so near \(\varphi=\pi/2\) it
    /// loses the whole cosine; \(c^2\) does not. Without `c`, \(x\) is exact and
    /// \((1-x)(1+x)\) has no cancellation (Sterbenz). For \(0\le m\le 1\),
    /// \(1-mx^2=(1-m)+m(1-x^2)\) is a sum of non-negative terms.
    fn amp_parts(
        xa: &Self,
        c: Option<&Self>,
        m: &Self,
        p: usize,
    ) -> Result<(Self, Self, Self), Error> {
        let one = Self::from_word(1, p)?;
        let x2 = xa.mul(xa, p, RoundingMode::None)?;
        let a = match c {
            Some(c) => c.mul(c, p, RoundingMode::None)?,
            None => one.sub(xa, p, RoundingMode::None)?.mul(
                &one.add(xa, p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )?,
        };
        let b = Self::one_minus_scaled(m, &x2, &a, p)?;
        Ok((x2, a, b))
    }

    /// \(1-\lambda x^2\) given \(x^2\) and \(a=1-x^2\); \((1-\lambda)+\lambda a\) when
    /// \(0\le\lambda\le1\) so no cancellation occurs.
    fn one_minus_scaled(l: &Self, x2: &Self, a: &Self, p: usize) -> Result<Self, Error> {
        let one = Self::from_word(1, p)?;
        if !l.is_negative() && l.cmp(&one) <= 0 {
            one.sub(l, p, RoundingMode::None)?.add(
                &l.mul(a, p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )
        } else {
            one.sub(&l.mul(x2, p, RoundingMode::None)?, p, RoundingMode::None)
        }
    }

    /// [`Self::elliptic_f_at`] with an optional cosine; returns the value and the
    /// bits lost to cancellation (none for F: every term is positive).
    fn elliptic_f_sc(
        &self,
        c: Option<&Self>,
        m: &Self,
        p: usize,
        cc: &mut Consts,
    ) -> Result<(Self, usize), Error> {
        let one = Self::from_word(1, p)?;
        let xa = self.abs()?;
        let sign_neg = self.is_negative();
        if m.is_zero() {
            let mut r = xa.asin(p, RoundingMode::None, cc)?;
            if sign_neg {
                r.set_sign(Sign::Neg);
            }
            return Ok((r, 0));
        }
        if m.cmp(&one) == 0 {
            if xa.cmp(&one) == 0 || c.is_some_and(|c| c.is_zero()) {
                return Err(Error::InvalidArgument);
            }
            // F(φ|1) = atanh(sin φ) = asinh(tan φ); the tangent keeps the cosine.
            let mut r = match c {
                Some(c) => xa
                    .div(c, p, RoundingMode::None)?
                    .asinh(p, RoundingMode::None, cc)?,
                None => xa.atanh(p, RoundingMode::None, cc)?,
            };
            if sign_neg {
                r.set_sign(Sign::Neg);
            }
            return Ok((r, 0));
        }
        let (_, a, b) = Self::amp_parts(&xa, c, m, p)?;
        if a.is_negative() || b.is_negative() {
            return Err(Error::InvalidArgument);
        }
        let rf = carlson_rf(&a, &b, &one, p, cc)?;
        let mut r = xa.mul(&rf, p, RoundingMode::None)?;
        if sign_neg {
            r.set_sign(Sign::Neg);
        }
        Ok((r, 0))
    }

    /// Incomplete \(E(x|m)\), \(|x|\le 1\). \(x=\sin\varphi\), not the angle
    /// \(\varphi\). \(m=k^2\), not the modulus \(k\).
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
        // E(φ|0) = φ, so E(x|0) = arcsin(x). Round that once.
        if m.is_zero() {
            let mut r = self.asin(p, rm, cc)?;
            r.set_inexact(r.inexact() | self.inexact() | m.inexact());
            return Ok(r);
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p
            .max(self.mantissa_max_bit_len())
            .max(m.mantissa_max_bit_len())
            + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE * 2;
            let (mut ret, loss) = self.elliptic_e_sc(None, m, p_x, cc)?;
            let claim = p_wrk.saturating_sub(loss) / WORD_BIT_SIZE * WORD_BIT_SIZE;
            if claim > p && ret.try_set_precision(p, rm, claim)? {
                ret.set_inexact(ret.inexact() | self.inexact() | m.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    /// Sine-form \(E(x|m)\) with an optional cosine (see [`Self::amp_parts`]).
    /// Returns the value and the bits lost where \(xR_F\) and \(mx^3R_D/3\) cancel;
    /// as \(m\to1,\ \varphi\to\pi/2\) both grow like a logarithm while E stays near 1.
    fn elliptic_e_sc(
        &self,
        c: Option<&Self>,
        m: &Self,
        p: usize,
        cc: &mut Consts,
    ) -> Result<(Self, usize), Error> {
        let one = Self::from_word(1, p)?;
        let xa = self.abs()?;
        let sign_neg = self.is_negative();
        if m.is_zero() {
            let mut r = xa.asin(p, RoundingMode::None, cc)?;
            if sign_neg {
                r.set_sign(Sign::Neg);
            }
            return Ok((r, 0));
        }
        if m.cmp(&one) == 0 {
            let mut r = xa;
            if sign_neg {
                r.set_sign(Sign::Neg);
            }
            return Ok((r, 0));
        }
        let (x2, a, b) = Self::amp_parts(&xa, c, m, p)?;
        if a.is_negative() || b.is_negative() {
            return Err(Error::InvalidArgument);
        }
        let three = Self::from_word(3, p)?;
        let rf = carlson_rf(&a, &b, &one, p, cc)?;
        let rd = carlson_rd(&a, &b, &one, p, cc)?;
        let t1 = xa.mul(&rf, p, RoundingMode::None)?;
        let t2 = m
            .mul(&xa, p, RoundingMode::None)?
            .mul(&x2, p, RoundingMode::None)?
            .div(&three, p, RoundingMode::None)?
            .mul(&rd, p, RoundingMode::None)?;
        let t = t1.sub(&t2, p, RoundingMode::None)?;
        let loss = Self::cancel_bits(&t1, &t2, &t);
        let t = if sign_neg { t.neg()? } else { t };
        Ok((t, loss))
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

    /// Incomplete \(\Pi(n;x|m)\). `self` is the characteristic \(n\).
    /// \(x=\sin\varphi\), not the angle \(\varphi\). \(m=k^2\), not the modulus \(k\).
    /// \(|x|>1\) is `InvalidArgument`.
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
        // x = ±1 is φ = ±π/2: the complete integral, odd in x. The incomplete path
        // divides by sqrt(1 - x²) = 0 there.
        if x.abs()?.cmp(&one) == 0 {
            let inexact = x.inexact();
            let mut r = if x.is_negative() {
                self.elliptic_pi_complete(m, p, crate::common::util::invert_rm_for_sign(rm), cc)?
                    .neg()?
            } else {
                self.elliptic_pi_complete(m, p, rm, cc)?
            };
            r.set_inexact(r.inexact() | inexact);
            return Ok(r);
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p
            .max(self.mantissa_max_bit_len())
            .max(x.mantissa_max_bit_len())
            .max(m.mantissa_max_bit_len())
            + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE * 2;
            let (mut ret, loss) = self.elliptic_pi_sc(x, None, m, p_x, cc)?;
            let claim = p_wrk.saturating_sub(loss) / WORD_BIT_SIZE * WORD_BIT_SIZE;
            if claim > p && ret.try_set_precision(p, rm, claim)? {
                ret.set_inexact(ret.inexact() | self.inexact() | x.inexact() | m.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    /// Sine-form \(\Pi(n;x|m)\) with an optional cosine (see [`Self::amp_parts`]).
    /// Returns the value and the bits lost where \(xR_F\) and \(nx^3R_J/3\) cancel
    /// (\(n<0\)).
    fn elliptic_pi_sc(
        &self,
        x: &Self,
        c: Option<&Self>,
        m: &Self,
        p: usize,
        cc: &mut Consts,
    ) -> Result<(Self, usize), Error> {
        let one = Self::from_word(1, p)?;
        let xa = x.abs()?;
        let sign_neg = x.is_negative();
        let (x2, a, b) = Self::amp_parts(&xa, c, m, p)?;
        let pv = Self::one_minus_scaled(self, &x2, &a, p)?;
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
            return Ok((r, 0));
        }
        let three = Self::from_word(3, p)?;
        let rf = carlson_rf(&a, &b, &one, p, cc)?;
        let rj = carlson_rj(&a, &b, &one, &pv, p, cc)?;
        let t1 = xa.mul(&rf, p, RoundingMode::None)?;
        let t2 = self
            .mul(&xa, p, RoundingMode::None)?
            .mul(&x2, p, RoundingMode::None)?
            .div(&three, p, RoundingMode::None)?
            .mul(&rj, p, RoundingMode::None)?;
        let t = t1.add(&t2, p, RoundingMode::None)?;
        let loss = Self::cancel_bits(&t1, &t2, &t);
        let t = if sign_neg { t.neg()? } else { t };
        Ok((t, loss))
    }

    /// Round `v` once to `p` in `rm`. Identity shortcuts (`m = 0`, a zero
    /// amplitude) use this so reduction cannot round a second time.
    fn round_once(
        v: &Self,
        p: usize,
        rm: RoundingMode,
        extra_inexact: bool,
    ) -> Result<Self, Error> {
        let mut v = v.try_rounded_to(p, rm)?;
        v.set_inexact(v.inexact() | extra_inexact);
        Ok(v)
    }

    /// Bits of \(|\varphi|\) above 1. Range reduction keeps this many extra bits.
    fn phi_mag_bits(phi: &Self) -> usize {
        let e = phi.exponent();
        if e <= 0 {
            return 0;
        }
        match usize::try_from(e) {
            Ok(bits) => bits,
            Err(_) => usize::MAX / 4,
        }
    }

    /// \((\sin r, k, \cos r)\) for \(\varphi = r + k\pi\), \(|r|\le\pi/2\).
    ///
    /// `Ok(None)` means the working precision did not keep the remainder inside
    /// that interval. The caller raises the precision and tries again.
    fn reduced_sin(&self, p: usize, cc: &mut Consts) -> Result<Option<(Self, Self, Self)>, Error> {
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let half = pi.div(&Self::from_word(2, WORD_BIT_SIZE)?, p, RoundingMode::None)?; // exact
        let q = self.div(&pi, p, RoundingMode::None)?;
        let k = q.round(0, RoundingMode::ToEven)?;
        if !k.is_int() {
            return Ok(None);
        }
        let r = self.sub(&k.mul(&pi, p, RoundingMode::None)?, p, RoundingMode::None)?;
        if r.abs()?.cmp(&half) > 0 {
            return Ok(None);
        }
        let s = r.sin(p, RoundingMode::None, cc)?;
        let one = Self::from_word(1, p)?;
        if s.abs()?.cmp(&one) > 0 {
            return Ok(None);
        }
        let c = r.cos(p, RoundingMode::None, cc)?.abs()?;
        Ok(Some((s, k, c)))
    }

    /// `part + 2k·complete` and the bits that sum cancels.
    fn add_turns(part: Self, k: &Self, complete: &Self, p: usize) -> Result<(Self, usize), Error> {
        let two = Self::from_word(2, p)?;
        let shift = two
            .mul(k, p, RoundingMode::None)?
            .mul(complete, p, RoundingMode::None)?;
        let t = part.add(&shift, p, RoundingMode::None)?;
        let loss = Self::cancel_bits(&part, &shift, &t);
        Ok((t, loss))
    }

    /// Bits of \(|\delta r|\) amplification when \(k\ne0\): the reduced angle carries
    /// an absolute error near \(2^{-p}|k\pi|\), and \(dF/d\varphi=1/\sqrt{1-m\sin^2}\le
    /// (1-m)^{-1/2}\), \(d\Pi/d\varphi\le(1-n)^{-1}(1-m)^{-1/2}\) for \(0<m,n<1\).
    /// The sum is at least the complete integral, so this bounds the relative loss.
    /// For \(m,n<0\) the derivative is at most 1 while the complete integral decays
    /// like \((1-m)^{-1/2}\) or \((1-n)^{-1/2}\); the same gap bits cover it.
    fn turn_cond(m: &Self, n: Option<&Self>, p: usize) -> Result<usize, Error> {
        let one = Self::from_word(1, p)?;
        // |log2(1 - l)| + 1: the derivative bound for 0 < l < 1, and for l < 0 the
        // complete integral's decay like (1 - l)^{-1/2} against a derivative <= 1.
        let gap_bits = |l: &Self| -> Result<usize, Error> {
            if l.is_zero() || l.cmp(&one) >= 0 {
                return Ok(0);
            }
            let g = one.sub(l, p, RoundingMode::None)?;
            Ok(g.exponent().unsigned_abs() as usize + 1)
        };
        let mut bits = gap_bits(m)?.div_ceil(2) + 2;
        if let Some(n) = n {
            bits += gap_bits(n)?;
        }
        Ok(bits)
    }

    fn complete_k_at(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        if self.is_zero() {
            let mut half = cc.pi_num(p, RoundingMode::None)?;
            half.div_by_2(RoundingMode::None);
            return Ok(half);
        }
        self.elliptic_k_at(p, cc)
    }

    fn complete_e_at(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let one = Self::from_word(1, p)?;
        if self.is_zero() {
            let mut half = cc.pi_num(p, RoundingMode::None)?;
            half.div_by_2(RoundingMode::None);
            return Ok(half);
        }
        if self.cmp(&one) == 0 {
            return Self::from_word(1, p);
        }
        self.elliptic_e_complete_at(p, cc)
    }

    fn elliptic_f_from_reduced(
        s: &Self,
        k: &Self,
        c: &Self,
        m: &Self,
        phi_sign: Sign,
        p: usize,
        cc: &mut Consts,
    ) -> Result<(Self, usize), Error> {
        let one = Self::from_word(1, p)?;
        if k.is_zero() {
            if m.cmp(&one) == 0 && s.abs()?.cmp(&one) == 0 {
                return Err(Error::ExponentOverflow(phi_sign));
            }
            return s.elliptic_f_sc(Some(c), m, p, cc);
        }
        if m.cmp(&one) > 0 {
            return Err(Error::InvalidArgument);
        }
        if m.cmp(&one) == 0 {
            return Err(Error::ExponentOverflow(phi_sign));
        }
        let (part, l0) = s.elliptic_f_sc(Some(c), m, p, cc)?;
        let complete = m.complete_k_at(p, cc)?;
        let (t, l1) = Self::add_turns(part, k, &complete, p)?;
        Ok((t, l0 + l1 + Self::turn_cond(m, None, p)?))
    }

    /// \(F(\varphi|m)\) for an amplitude \(\varphi\). See [`ExactNum::elliptic_f_phi`].
    pub fn elliptic_f_phi(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        if self.is_zero() {
            return Self::new2(p, Sign::Pos, self.inexact() | m.inexact());
        }
        // F(φ|0) = φ. Round the amplitude once. Reduction through sin and
        // arcsin is a second rounding, and directed modes then fail the
        // roundability test.
        if m.is_zero() {
            return Self::round_once(self, p, rm, m.inexact());
        }
        let mag = Self::phi_mag_bits(self);
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = round_p(
            p.max(self.mantissa_max_bit_len())
                .max(m.mantissa_max_bit_len())
                .saturating_add(mag)
                .saturating_add(p_inc),
        );
        loop {
            let p_x = p_wrk.saturating_add(WORD_BIT_SIZE.saturating_mul(2));
            if let Some((s, k, c)) = self.reduced_sin(p_x, cc)? {
                let (mut ret, loss) =
                    Self::elliptic_f_from_reduced(&s, &k, &c, m, self.sign(), p_x, cc)?;
                let claim = p_wrk.saturating_sub(loss) / WORD_BIT_SIZE * WORD_BIT_SIZE;
                if claim > p && ret.try_set_precision(p, rm, claim)? {
                    ret.set_inexact(ret.inexact() | self.inexact() | m.inexact());
                    return Ok(ret);
                }
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn elliptic_e_from_reduced(
        s: &Self,
        k: &Self,
        c: &Self,
        m: &Self,
        p: usize,
        cc: &mut Consts,
    ) -> Result<(Self, usize), Error> {
        let one = Self::from_word(1, p)?;
        if k.is_zero() {
            return s.elliptic_e_sc(Some(c), m, p, cc);
        }
        if m.cmp(&one) > 0 {
            return Err(Error::InvalidArgument);
        }
        let (part, l0) = s.elliptic_e_sc(Some(c), m, p, cc)?;
        let complete = m.complete_e_at(p, cc)?;
        let (t, l1) = Self::add_turns(part, k, &complete, p)?;
        // dE/dφ = sqrt(1 - m sin²φ) <= sqrt(1 - m) (m < 0), and the sum is >= E(m) >= 1.
        let grow = if m.is_negative() {
            let g = one.sub(m, p, RoundingMode::None)?;
            (g.exponent().max(0) as usize).div_ceil(2)
        } else {
            0
        };
        Ok((t, l0 + l1 + grow + 2))
    }

    /// \(E(\varphi|m)\) for an amplitude \(\varphi\). See [`ExactNum::elliptic_e_phi`].
    pub fn elliptic_e_phi(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        if self.is_zero() {
            return Self::new2(p, Sign::Pos, self.inexact() | m.inexact());
        }
        // E(φ|0) = φ. Same single rounding as `elliptic_f_phi`.
        if m.is_zero() {
            return Self::round_once(self, p, rm, m.inexact());
        }
        let mag = Self::phi_mag_bits(self);
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = round_p(
            p.max(self.mantissa_max_bit_len())
                .max(m.mantissa_max_bit_len())
                .saturating_add(mag)
                .saturating_add(p_inc),
        );
        loop {
            let p_x = p_wrk.saturating_add(WORD_BIT_SIZE.saturating_mul(2));
            if let Some((s, k, c)) = self.reduced_sin(p_x, cc)? {
                let (mut ret, loss) = Self::elliptic_e_from_reduced(&s, &k, &c, m, p_x, cc)?;
                let claim = p_wrk.saturating_sub(loss) / WORD_BIT_SIZE * WORD_BIT_SIZE;
                if claim > p && ret.try_set_precision(p, rm, claim)? {
                    ret.set_inexact(ret.inexact() | self.inexact() | m.inexact());
                    return Ok(ret);
                }
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn elliptic_pi_from_reduced(
        n: &Self,
        s: &Self,
        k: &Self,
        c: &Self,
        m: &Self,
        p: usize,
        cc: &mut Consts,
    ) -> Result<(Self, usize), Error> {
        let one = Self::from_word(1, p)?;
        if k.is_zero() {
            return n.elliptic_pi_sc(s, Some(c), m, p, cc);
        }
        if m.cmp(&one) >= 0 || n.cmp(&one) >= 0 {
            return Err(Error::InvalidArgument);
        }
        let (part, l0) = n.elliptic_pi_sc(s, Some(c), m, p, cc)?;
        let complete = if n.is_zero() {
            m.complete_k_at(p, cc)?
        } else {
            n.elliptic_pi_complete_at(m, p, cc)?
        };
        let (t, l1) = Self::add_turns(part, k, &complete, p)?;
        Ok((t, l0 + l1 + Self::turn_cond(m, Some(n), p)?))
    }

    /// \(\Pi(n;\varphi|m)\) for an amplitude \(\varphi\). `self` is \(n\).
    /// See [`ExactNum::elliptic_pi_phi`].
    pub fn elliptic_pi_phi(
        &self,
        phi: &Self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        // Π(0; φ|m) = F(φ|m). Delegate so n = 0 is not a second rounding of
        // the Carlson value that F already rounds.
        if self.is_zero() {
            return phi.elliptic_f_phi(m, p, rm, cc);
        }
        if phi.is_zero() {
            return Self::new2(p, Sign::Pos, self.inexact() | phi.inexact() | m.inexact());
        }
        let mag = Self::phi_mag_bits(phi);
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = round_p(
            p.max(self.mantissa_max_bit_len())
                .max(phi.mantissa_max_bit_len())
                .max(m.mantissa_max_bit_len())
                .saturating_add(mag)
                .saturating_add(p_inc),
        );
        loop {
            let p_x = p_wrk.saturating_add(WORD_BIT_SIZE.saturating_mul(2));
            if let Some((s, k, c)) = phi.reduced_sin(p_x, cc)? {
                let (mut ret, loss) = Self::elliptic_pi_from_reduced(self, &s, &k, &c, m, p_x, cc)?;
                let claim = p_wrk.saturating_sub(loss) / WORD_BIT_SIZE * WORD_BIT_SIZE;
                if claim > p && ret.try_set_precision(p, rm, claim)? {
                    ret.set_inexact(ret.inexact() | self.inexact() | phi.inexact() | m.inexact());
                    return Ok(ret);
                }
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    /// `a + b` with no rounding: the width spans both operands' bits. When that
    /// width would be absurd the exponents are far apart, nothing cancels, and a
    /// rounding at `p` keeps a few-ulp relative error.
    fn sum_exact_or(a: &Self, b: &Self, p: usize) -> Result<Self, Error> {
        if a.is_zero() || b.is_zero() {
            let w = p
                .max(a.mantissa_max_bit_len())
                .max(b.mantissa_max_bit_len());
            return a.add(b, round_p(w), RoundingMode::None);
        }
        let (ea, eb) = (a.exponent(), b.exponent());
        let lo = (ea - a.mantissa_max_bit_len() as i64).min(eb - b.mantissa_max_bit_len() as i64);
        let width = ea.max(eb) + 2 - lo;
        let w = if width > 0 && width < (1 << 22) { (width as usize).max(p) } else { p };
        a.add(b, round_p(w), RoundingMode::None)
    }

    /// `a·b` with no rounding.
    fn mul_exact(a: &Self, b: &Self) -> Result<Self, Error> {
        let w = round_p(a.mantissa_max_bit_len() + b.mantissa_max_bit_len());
        a.mul(b, w, RoundingMode::None)
    }

    /// Inverse Jacobi \(\mathrm{arccn}(x|m)=F(\arccos x|m)\), `self = x`.
    ///
    /// With \(\cos\varphi=x\): \(\sin^2\varphi=1-x^2\), \(1-m\sin^2\varphi=(1-m)+mx^2\), all
    /// formed exactly, so \(F=\sqrt{1-x^2}\,R_F(x^2,(1-m)+mx^2,1)\) loses nothing to
    /// cancellation near \(x=\pm1\) (where \(\sqrt{1-x^2}\) from a rounded square
    /// would). For \(x<0\), \(\mathrm{arccn}=2K-F(\arccos|x|\,|m)\) with \(F\le K\).
    pub fn jacobi_arccn_num(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        let one = Self::from_word(1, p)?;
        let xa = self.abs()?;
        if xa.cmp(&one) > 0 {
            return Err(Error::InvalidArgument);
        }
        let neg = self.is_negative() && !self.is_zero();
        if neg && m.cmp(&one) >= 0 {
            return Err(Error::InvalidArgument);
        }
        if !neg && xa.cmp(&one) == 0 {
            return Self::new2(p, Sign::Pos, self.inexact() | m.inexact());
        }
        let x2 = Self::mul_exact(&xa, &xa)?;
        let s2 = Self::sum_exact_or(&one, &x2.neg()?, p)?;
        let mx2 = Self::mul_exact(m, &x2)?;
        let one_m = Self::sum_exact_or(&one, &m.neg()?, p)?;
        let b = Self::sum_exact_or(&one_m, &mx2, p)?;
        if b.is_negative() {
            return Err(Error::InvalidArgument);
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = round_p(
            p.max(self.mantissa_max_bit_len())
                .max(m.mantissa_max_bit_len())
                .saturating_add(p_inc),
        );
        loop {
            let p_x = p_wrk.saturating_add(WORD_BIT_SIZE * 2);
            let one_x = Self::from_word(1, p_x)?;
            let rf = carlson_rf(&x2, &b, &one_x, p_x, cc)?;
            let f = s2
                .sqrt(p_x, RoundingMode::None)?
                .mul(&rf, p_x, RoundingMode::None)?;
            let (mut ret, loss) = if neg {
                let k2 = m.complete_k_at(p_x, cc)?.mul(
                    &Self::from_word(2, p_x)?,
                    p_x,
                    RoundingMode::None,
                )?;
                let t = k2.sub(&f, p_x, RoundingMode::None)?;
                let l = Self::cancel_bits(&k2, &f, &t);
                (t, l)
            } else {
                (f, 0)
            };
            let claim = p_wrk.saturating_sub(loss) / WORD_BIT_SIZE * WORD_BIT_SIZE;
            if claim > p && ret.try_set_precision(p, rm, claim)? {
                ret.set_inexact(ret.inexact() | self.inexact() | m.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    /// Inverse Jacobi \(\mathrm{arcdn}(x|m)\) for \(m>0\), `self = x`, on \([0,K]\).
    ///
    /// With \(\mathrm{dn}=x\): \(\sin^2\varphi=(1-x^2)/m\), \(\cos^2\varphi=(x^2+m-1)/m\),
    /// \(1-m\sin^2\varphi=x^2\). The numerators are formed exactly, so the
    /// cancellation in \(x^2+m-1\) near \(x=\sqrt{1-m}\) (\(u\to K\)) and in \(1-x^2\)
    /// near \(x=1\) costs nothing; \(F=\sqrt{\sin^2\varphi}\,R_F(\cos^2\varphi,x^2,1)\).
    pub fn jacobi_arcdn_num(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;
        if !m.is_positive() || m.is_zero() {
            return Err(Error::InvalidArgument);
        }
        let one = Self::from_word(1, p)?;
        let xa = self.abs()?;
        if xa.cmp(&one) > 0 {
            return Err(Error::InvalidArgument);
        }
        if xa.cmp(&one) == 0 {
            return Self::new2(p, Sign::Pos, self.inexact() | m.inexact());
        }
        // arcdn(0|1) = K(1) = +∞ (dn(u|1) = sech(u) → 0 as u → +∞). A pole,
        // not a finite overflow: ExponentOverflow maps to +Inf in every mode
        // through result_to_ext. Both signs of zero.
        if xa.is_zero() && m.cmp(&one) == 0 {
            return Err(Error::ExponentOverflow(Sign::Pos));
        }
        let x2 = Self::mul_exact(&xa, &xa)?;
        let s2n = Self::sum_exact_or(&one, &x2.neg()?, p)?;
        let m1 = Self::sum_exact_or(m, &one.neg()?, p)?;
        let c2n = Self::sum_exact_or(&x2, &m1, p)?;
        if c2n.is_negative() {
            return Err(Error::InvalidArgument);
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = round_p(
            p.max(self.mantissa_max_bit_len())
                .max(m.mantissa_max_bit_len())
                .saturating_add(p_inc),
        );
        loop {
            let p_x = p_wrk.saturating_add(WORD_BIT_SIZE * 2);
            let one_x = Self::from_word(1, p_x)?;
            let s2 = s2n.div(m, p_x, RoundingMode::None)?;
            let c2 = c2n.div(m, p_x, RoundingMode::None)?;
            let rf = carlson_rf(&c2, &x2, &one_x, p_x, cc)?;
            let mut ret = s2
                .sqrt(p_x, RoundingMode::None)?
                .mul(&rf, p_x, RoundingMode::None)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact() | m.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
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
            let mut r = if n.is_multiple_of(2) { one } else { one.neg()? };
            r.set_inexact(self.inexact());
            return Ok(r);
        }
        let extra = extra_bits_for_degree(n);
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc + extra;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE + extra;
            let v = self.legendre_p_at(n, p_x)?;
            if let Some(ret) = self.ziv_accept_or_self(v, p, rm, p_wrk)? {
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    /// `None` stands for \(P_1(x) = x\): the caller rounds `self` itself.
    fn legendre_p_at(&self, n: u32, p: usize) -> Result<Option<Self>, Error> {
        if n == 0 {
            return Self::from_word(1, p).map(Some);
        }
        if n == 1 {
            return Ok(None);
        }
        // P_{k-1}, P_k: None stands for P_1 = x itself (no copy).
        let mut pm2: Option<Self> = Some(Self::from_word(1, p)?);
        let mut pm1: Option<Self> = None;
        let mut k = 1u32;
        while k < n {
            let kk = Self::from_word(k as Word, p)?;
            let kp1 = Self::from_word((k + 1) as Word, p)?;
            let two_k1 = Self::from_word((2 * k + 1) as Word, p)?;
            let pnext = two_k1
                .mul(self, p, RoundingMode::None)?
                .mul(pm1.as_ref().unwrap_or(self), p, RoundingMode::None)?
                .sub(
                    &kk.mul(pm2.as_ref().unwrap_or(self), p, RoundingMode::None)?,
                    p,
                    RoundingMode::None,
                )?
                .div(&kp1, p, RoundingMode::None)?;
            pm2 = pm1;
            pm1 = Some(pnext);
            k += 1;
        }
        pm1.ok_or(Error::InvalidArgument).map(Some) // n >= 2 ran the loop at least once
    }

    /// One Ziv pass over a kernel value. `None` stands for `self` itself, which is
    /// rounded from a borrow; the returned value carries `self`'s inexact flag.
    fn ziv_accept_or_self(
        &self,
        v: Option<Self>,
        p: usize,
        rm: RoundingMode,
        p_wrk: usize,
    ) -> Result<Option<Self>, Error> {
        let ret = match v {
            Some(mut r) => {
                if r.try_set_precision(p, rm, p_wrk)? {
                    Some(r)
                } else {
                    None
                }
            }
            None => self.try_rounded_checked(p, rm, p_wrk)?,
        };
        Ok(ret.map(|mut r| {
            r.set_inexact(r.inexact() | self.inexact());
            r
        }))
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
            let v = self.assoc_legendre_p_at(n, m, p_x)?;
            if let Some(ret) = self.ziv_accept_or_self(v, p, rm, p_wrk)? {
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    /// `None` stands for `self` itself (\(m = 0\), \(n = 1\)).
    fn assoc_legendre_p_at(&self, n: u32, m: i32, p: usize) -> Result<Option<Self>, Error> {
        let am = m.unsigned_abs();
        if m == 0 {
            return self.legendre_p_at(n, p);
        }
        let ppos = self.assoc_pos(n, am, p)?;
        if m > 0 {
            return Ok(Some(ppos));
        }
        let ratio = fact_ratio_down(n - am, n + am, p)?;
        let mut out = ratio.mul(&ppos, p, RoundingMode::None)?;
        if am % 2 == 1 {
            out.set_sign(if out.is_positive() { Sign::Neg } else { Sign::Pos });
        }
        Ok(Some(out))
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

    /// Gaussian \({}_2F_1(a=\mathrm{self},b;c;z)\) for real arguments.
    ///
    /// - Terminating (\(a\) or \(b\) a non-positive integer): the polynomial, for any \(z\).
    /// - \(\lvert z\rvert\le 1/2\): the Gauss series.
    /// - \(z\in(1/2,1)\): the \(1-z\) connection formula (DLMF 15.8.4); when \(c-a-b\) is an
    ///   integer, its logarithmic limit (DLMF 15.8.10, A&S 15.3.10–15.3.12).
    /// - \(z=1\): Gauss's sum when \(c-a-b>0\); `InvalidArgument` otherwise (divergent).
    /// - \(z<-1/2\): Pfaff, \((1-z)^{-a}\,{}_2F_1(a,c-b;c;z/(z-1))\) with \(z/(z-1)\in(1/3,1)\).
    ///
    /// `InvalidArgument` for non-terminating \(z>1\) (not real in general) and for an uncancelled
    /// pole at a non-positive integer \(c\). The series runs until the term bound, or until
    /// terms start growing after \(n > \max(\lvert a\rvert,\lvert b\rvert,\lvert c\rvert)\)
    /// (the series is then not in its disk). Cancellation between the connection terms is
    /// measured and the evaluation repeated with that many extra bits.
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
        if let Some(mut v) = self.hyp2f1_terminating_exact(b, c, z, p, rm)? {
            v.set_inexact(v.inexact() | self.inexact() | b.inexact() | c.inexact() | z.inexact());
            return Ok(v);
        }
        if let Some(mut v) = self.hyp2f1_closed_exact(b, c, z, p, rm, cc)? {
            v.set_inexact(v.inexact() | self.inexact() | b.inexact() | c.inexact() | z.inexact());
            return Ok(v);
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

    /// Terminating ₂F₁ (ZF-70): with `a` or `b = −N`, the value is the polynomial
    /// `Σ_{k≤N} (a)_k (b)_k / ((c)_k k!) z^k`, a rational number. It is summed exactly, at
    /// whatever width the exact rational needs, by binary splitting over
    /// `r_k z = (a+k−1)(b+k−1) z / ((c+k−1) k)` (one numerator over one denominator), and rounded
    /// once, so dyadic results such as ₂F₁(−3, −½; 2; 2) = 21/16 are exact in every mode and no
    /// size threshold sends an exactly representable value to the floating path. `None` (the
    /// floating path decides) only for a pole `c + k − 1 = 0`. A real allocation failure
    /// while storing the exact rational is `Err(MemoryAllocation)` from `try_reserve`; there
    /// is no size cap (the 2^46-bit / 8 TiB guard was removed).
    fn hyp2f1_terminating_exact(
        &self,
        b: &Self,
        c: &Self,
        z: &Self,
        p: usize,
        rm: RoundingMode,
    ) -> Result<Option<Self>, Error> {
        let Some((num, q)) = self.hyp2f1_terminating_ratio(b, c, z)? else {
            return Ok(None);
        };
        if num.is_zero() {
            return Ok(Some(Self::new2(p, Sign::Pos, false)?));
        }
        Ok(Some(num.div(&q, p, rm)?))
    }

    /// The exact terminating ₂F₁ of [`Self::hyp2f1_terminating_exact`] as `(numerator,
    /// denominator)`.
    fn hyp2f1_terminating_ratio(
        &self,
        b: &Self,
        c: &Self,
        z: &Self,
    ) -> Result<Option<(Self, Self)>, Error> {
        let degree = |x: &Self| -> Option<usize> {
            if x.is_zero() || x.is_positive() || !x.is_int() {
                return None;
            }
            // A degree past usize is past any memory: report it below as allocation failure.
            Some(x.neg().ok()?.int_as_usize().unwrap_or(usize::MAX))
        };
        let (n_a, n_b) = (degree(self), degree(b));
        let n = match (n_a, n_b) {
            (Some(x), Some(y)) => x.min(y),
            (Some(x), None) | (None, Some(x)) => x,
            (None, None) => return Ok(None),
        };
        // Pole c + k − 1 = 0 for some 1 <= k <= n: c a non-positive integer with −c < n.
        if c.is_int() && !c.is_positive() {
            let mc = if c.is_zero() { 0 } else { c.neg()?.int_as_usize().unwrap_or(usize::MAX) };
            if mc < n {
                return Ok(None);
            }
        }
        if n == 0 {
            let one = Self::from_word(1, WORD_BIT_SIZE)?;
            return Ok(Some((Self::from_word(1, WORD_BIT_SIZE)?, one)));
        }
        // The stored inputs are exact values; their inexact flags are ORed in by the caller.
        let exact = |x: &Self| -> Result<Self, Error> {
            let mut v = x.clone()?;
            v.set_inexact(false);
            Ok(v)
        };
        let (ea, eb, ec, ez) = (
            trim_low_zero_words(exact(self)?)?,
            trim_low_zero_words(exact(b)?)?,
            trim_low_zero_words(exact(c)?)?,
            trim_low_zero_words(exact(z)?)?,
        );
        let ctx = Hyp2f1Split {
            a: &ea,
            b: &eb,
            c: &ec,
            z: &ez,
        };
        let (_, q, t) = ctx.split(1, n + 1)?;
        // F = 1 + T/Q = (Q + T) / Q
        Ok(Some((exact_add_num(&q, &t)?, q)))
    }

    /// ZF-37: closed forms of a non-terminating ₂F₁ whose value is often an exact dyadic that
    /// the Ziv loop can never certify; each is one exact rational rounded once.
    /// - z = 1 with `a = m` (or `b = m`) a positive integer and `c − a − b > 0` (Gauss,
    ///   DLMF 15.4.20): `F = (c−m)_m / (c−b−m)_m`, e.g. ₂F₁(1, 1; 10; 1) = 9/8.
    /// - z < 1 with `c − a` (or `c − b`) a non-positive integer (Euler, DLMF 15.8.1):
    ///   `F = (1−z)^{c−a−b} ₂F₁(c−a, c−b; c; z)`, a terminating polynomial times a power;
    ///   used when the power is exact, e.g. ₂F₁(10, ½; 10; 1 − 2^-20) = 2^10 (c = a is the
    ///   degree-0 case `(1−z)^{−b}`).
    ///
    /// `None` (the floating path decides) otherwise, or when the exact work width is unallocatable (then the series / Ziv path decides).
    fn hyp2f1_closed_exact(
        &self,
        b: &Self,
        c: &Self,
        z: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Option<Self>, Error> {
        let exact = |v: &Self| -> Result<Self, Error> {
            let mut v = v.clone()?;
            v.set_inexact(false);
            Ok(v)
        };
        let (ea, eb, ec, ez) = (exact(self)?, exact(b)?, exact(c)?, exact(z)?);
        let one = Self::from_word(1, WORD_BIT_SIZE)?;
        let sub = |x: &Self, y: &Self| -> Result<Self, Error> { exact_add_num(x, &y.neg()?) };
        let vs_one = ez.cmp(&one);
        if vs_one > 0 {
            return Ok(None);
        }
        if vs_one == 0 {
            let s = sub(&sub(&ec, &ea)?, &eb)?;
            if !s.is_positive() || (ec.is_int() && !ec.is_positive()) {
                return Ok(None);
            }
            for (m_src, other) in [(&ea, &eb), (&eb, &ea)] {
                if !m_src.is_positive() || !m_src.is_int() {
                    continue;
                }
                let Ok(m) = m_src.int_as_usize() else {
                    continue;
                };
                let per = ec
                    .mantissa_max_bit_len()
                    .max(other.mantissa_max_bit_len())
                    .saturating_add(3 * WORD_BIT_SIZE);
                let Some(w) = m.checked_mul(per) else {
                    continue;
                };
                if Self::p_assertion(round_p(w)).is_err() {
                    continue;
                }
                let n0 = sub(&ec, m_src)?;
                let d0 = sub(&n0, other)?;
                let (mut num, mut den) = (
                    Self::from_word(1, WORD_BIT_SIZE)?,
                    Self::from_word(1, WORD_BIT_SIZE)?,
                );
                for k in 0..m {
                    let kw = Self::from_word(k as Word, WORD_BIT_SIZE)?;
                    num = exact_mul_num(&num, &exact_add_num(&n0, &kw)?)?;
                    den = exact_mul_num(&den, &exact_add_num(&d0, &kw)?)?;
                }
                if num.is_zero() {
                    return Ok(Some(Self::new2(p, Sign::Pos, false)?));
                }
                return Ok(Some(num.div(&den, p, rm)?));
            }
            return Ok(None);
        }
        let (ca, cb) = (sub(&ec, &ea)?, sub(&ec, &eb)?);
        let nonpos_int = |v: &Self| v.is_zero() || (v.is_int() && v.is_negative());
        if !nonpos_int(&ca) && !nonpos_int(&cb) {
            return Ok(None);
        }
        let ze = usize::try_from(ez.exponent().unsigned_abs()).unwrap_or(usize::MAX);
        let w_omz = ez
            .mantissa_max_bit_len()
            .saturating_add(ze)
            .saturating_add(WORD_BIT_SIZE);
        if Self::p_assertion(round_p(w_omz)).is_err() {
            return Ok(None);
        }
        let omz = one.sub(&ez, round_p(w_omz), RoundingMode::None)?;
        let deg = |v: &Self| -> usize {
            if nonpos_int(v) {
                v.neg()
                    .ok()
                    .and_then(|n| n.int_as_usize().ok())
                    .unwrap_or(usize::MAX)
            } else {
                usize::MAX
            }
        };
        let n = deg(&ca).min(deg(&cb));
        let per = 3
            * (ca
                .mantissa_max_bit_len()
                .max(cb.mantissa_max_bit_len())
                .max(ec.mantissa_max_bit_len())
                .max(ez.mantissa_max_bit_len())
                + 2 * WORD_BIT_SIZE);
        // n=0 is Euler F(a,b;a;z)=(1-z)^{-b}: no Pochhammer product, bit length 1.
        // After WORK_BITS_MAX, p_assertion(round_p(0*per)) fails (precision 0) and
        // skipped the exact path; Ziv then cannot certify the exact dyadic and
        // returns PrecisionRetryExhausted (NaN at the ExactNum wrapper). That is
        // the zf_37 regression on F(1,-1/2;1;1-2^{-20})=2^{-10}.
        if n > 0 {
            let Some(w) = n.checked_mul(per) else {
                return Ok(None);
            };
            if Self::p_assertion(round_p(w)).is_err() {
                return Ok(None);
            }
        }
        let (num, q) = if n == 0 {
            (
                Self::from_word(1, WORD_BIT_SIZE)?,
                Self::from_word(1, WORD_BIT_SIZE)?,
            )
        } else {
            let Some(r) = ca.hyp2f1_terminating_ratio(&cb, &ec, &ez)? else {
                return Ok(None);
            };
            r
        };
        let d = sub(&sub(&ec, &ea)?, &eb)?;
        let pw = if d.is_zero() {
            Self::from_word(1, WORD_BIT_SIZE)?
        } else {
            let d_up = if d.exponent() <= 0 {
                1
            } else if d.exponent() < 32 {
                1usize << d.exponent()
            } else {
                return Ok(None);
            };
            let wa = omz
                .mantissa_max_bit_len()
                .saturating_mul(d_up)
                .saturating_add(WORD_BIT_SIZE);
            // As in `betainc_exact`: a breakpoint value v = pw·num/q bounds pw's odd part by
            // p + 2 + bits(q); a power inexact there is not such a value.
            let wl = p
                .saturating_add(2)
                .saturating_add(q.mantissa_max_bit_len())
                .saturating_add(2 * WORD_BIT_SIZE)
                .min(wa);
            if Self::p_assertion(round_p(wl)).is_err() {
                return Ok(None);
            }
            let v = omz.pow(&d, round_p(wl), RoundingMode::ToZero, cc)?;
            if v.inexact() {
                return Ok(None);
            }
            v
        };
        let num = exact_mul_num(&pw, &num)?;
        if num.is_zero() {
            return Ok(Some(Self::new2(p, Sign::Pos, false)?));
        }
        Ok(Some(num.div(&q, p, rm)?))
    }

    fn hypergeom_2f1_at(
        &self,
        b: &Self,
        c: &Self,
        z: &Self,
        p: usize,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let (v, lost) = hyp2f1_real(self, b, c, z, p, cc)?;
        if lost <= WORD_BIT_SIZE {
            return Ok(v);
        }
        let p_more = p + round_p(lost) + WORD_BIT_SIZE;
        if Self::p_assertion(p_more).is_err() {
            return Ok(v);
        }
        Ok(hyp2f1_real(self, b, c, z, p_more, cc)?.0)
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
        if !self.is_positive() || !b.is_positive() || self.is_zero() || b.is_zero() {
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
        if self.cmp(b) == 0 && x.cmp(&one_half(p)?) == 0 {
            // I_{1/2}(a,a) = 1/2 exactly; the rounding test cannot certify an exact value.
            let mut half = one_half(p)?;
            half.set_inexact(self.inexact() | b.inexact() | x.inexact());
            return Ok(half);
        }
        if let Some(mut v) = self.betainc_exact(b, x, p, rm, cc)? {
            v.set_inexact(v.inexact() | self.inexact() | b.inexact() | x.inexact());
            return Ok(v);
        }
        // ZF-37: near x = 1, I_x(a,b) = 1 − I_{1−x}(b,a) with the complement below a quarter
        // ulp of 1 (e.g. I_{1−2^-256}(½, 100) = 1 − O(2^-25600)): a value 1 − δ, δ > 0, that
        // Ziv cannot decide; rounded once by `round_offset`. The complement comes from the
        // direct (relative-accuracy) evaluation at 1 − x, which is exact here.
        if x.cmp(&one_half(p)?) > 0 {
            let xe = usize::try_from(x.exponent().unsigned_abs()).unwrap_or(usize::MAX);
            let w = x
                .mantissa_max_bit_len()
                .saturating_add(xe)
                .saturating_add(WORD_BIT_SIZE);
            let one = Self::from_word(1, WORD_BIT_SIZE)?;
            let mut ox = one.sub(x, round_p(w), RoundingMode::None)?;
            ox.set_inexact(false);
            let q = round_p(p + 2 * WORD_BIT_SIZE);
            let t = b.betainc_direct(self, &ox, q, cc)?;
            if t.is_positive() && t.exponent() < -(p as Exponent) - 4 {
                return Self::round_offset(&one, Sign::Neg, p, rm);
            }
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

    /// ZF-37: an integer parameter makes I_x a closed form (DLMF 8.17.5):
    /// `I_x(a, n) = x^a Σ_{k<n} (a)_k/k! (1−x)^k` and `I_x(n, b) = 1 − I_{1−x}(b, n)`. When the
    /// power (`x^a`, resp. `(1−x)^b`) is an exact dyadic, the value is the rational
    /// `x^a W / (n−1)!`, with `W` summed exactly by the division-free Horner recurrence
    /// `P_k = k P_{k+1}`, `W_k = P_k + (a+k−1) y W_{k+1}` (`P_n = W_n = 1`), and rounded once:
    /// exact results such as `I_{1/2}(100, 1) = 2^-100` or `I_{2^-30}(½, 1) = 2^-15`, which no
    /// Ziv loop can certify, come out correctly rounded in every mode. `None` (the floating
    /// path decides; its values are irrational there or Ziv reports exhaustion, never a
    /// wrong value) when the power is inexact or the exact work width is unallocatable.
    fn betainc_exact(
        &self,
        b: &Self,
        x: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Option<Self>, Error> {
        let exact = |v: &Self| -> Result<Self, Error> {
            let mut v = v.clone()?;
            v.set_inexact(false);
            Ok(v)
        };
        let (ea, eb, ex) = (exact(self)?, exact(b)?, exact(x)?);
        let one = Self::from_word(1, WORD_BIT_SIZE)?;
        // 1 − x exactly (x ∈ (0, 1)).
        let xe = usize::try_from(ex.exponent().unsigned_abs()).unwrap_or(usize::MAX);
        let w_omx = ex
            .mantissa_max_bit_len()
            .saturating_add(xe)
            .saturating_add(WORD_BIT_SIZE);
        if Self::p_assertion(round_p(w_omx)).is_err() {
            return Ok(None);
        }
        let omx = one.sub(&ex, round_p(w_omx), RoundingMode::None)?;
        for complement in [false, true] {
            // I_base-form: alpha = the non-integer side's parameter, n the integer one,
            // power base^alpha, series variable y.
            let (alpha, n_src, base, y) =
                if complement { (&eb, &ea, &omx, &ex) } else { (&ea, &eb, &ex, &omx) };
            if !n_src.is_int() {
                continue;
            }
            let Ok(n) = n_src.int_as_usize() else {
                continue;
            };
            let per = alpha
                .mantissa_max_bit_len()
                .saturating_add(y.mantissa_max_bit_len())
                .saturating_add(2 * WORD_BIT_SIZE);
            let Some(w) = n.checked_mul(per) else {
                continue;
            };
            if Self::p_assertion(round_p(w)).is_err() {
                continue;
            }
            let alpha_up = if alpha.exponent() <= 0 {
                1
            } else if alpha.exponent() < 32 {
                1usize << alpha.exponent()
            } else {
                continue;
            };
            let wa = base
                .mantissa_max_bit_len()
                .saturating_mul(alpha_up)
                .saturating_add(WORD_BIT_SIZE);
            let mut pk = Self::from_word(1, WORD_BIT_SIZE)?;
            let mut wk = Self::from_word(1, WORD_BIT_SIZE)?;
            for k in (1..n).rev() {
                let kw = Self::from_word(k as Word, WORD_BIT_SIZE)?;
                pk = exact_mul_num(&pk, &kw)?;
                let akm1 = exact_add_num(alpha, &Self::from_word((k - 1) as Word, WORD_BIT_SIZE)?)?;
                let t = exact_mul_num(&exact_mul_num(&akm1, y)?, &wk)?;
                wk = exact_add_num(&pk, &t)?;
            }
            // Only a value with at most p + 2 significant bits can sit on a rounding
            // breakpoint (the case Ziv cannot certify). With v = base^α W / P that bounds the
            // odd part of base^α by p + 2 + bits(P) bits; for the complement v = 1 − v′ the
            // bound grows by the exponent gap −log₂ v, estimated at 128 bits (+ slack). A
            // power inexact at that width means v is not a breakpoint value: the floating
            // path certifies it.
            let lp = pk.mantissa_max_bit_len();
            let mut gap = 0usize;
            if complement {
                let q = round_p(2 * WORD_BIT_SIZE);
                let approx = base
                    .pow(alpha, q, RoundingMode::None, cc)?
                    .mul(&wk, q, RoundingMode::None)?
                    .div(&pk, q, RoundingMode::None)?;
                let v = one.sub(&approx, q, RoundingMode::None)?;
                gap = if v.is_positive() && v.exponent() > -(WORD_BIT_SIZE as Exponent) {
                    usize::try_from(-v.exponent()).unwrap_or(0)
                } else {
                    usize::MAX
                };
            }
            let wl = p
                .saturating_add(2)
                .saturating_add(lp)
                .saturating_add(gap)
                .saturating_add(2 * WORD_BIT_SIZE)
                .min(wa);
            if Self::p_assertion(round_p(wl)).is_err() {
                continue;
            }
            let xa = base.pow(alpha, round_p(wl), RoundingMode::ToZero, cc)?;
            if xa.inexact() {
                continue;
            }
            let num = exact_mul_num(&xa, &wk)?;
            let num = if complement { exact_add_num(&pk, &num.neg()?)? } else { num };
            if num.is_zero() {
                return Ok(Some(Self::new2(p, Sign::Pos, false)?));
            }
            return Ok(Some(num.div(&pk, p, rm)?));
        }
        Ok(None)
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
        let one = Self::from_word(1, p)?;
        let ga = self.gamma_at(p, cc)?;
        let gb = b.gamma_at(p, cc)?;
        let gab = self.add(b, p, RoundingMode::None)?.gamma_at(p, cc)?;
        let beta = ga
            .mul(&gb, p, RoundingMode::None)?
            .div(&gab, p, RoundingMode::None)?;
        let xa = x.pow(self, p, RoundingMode::None, cc)?;
        let one_b = one.sub(b, p, RoundingMode::None)?;
        if terminating_neg_int(&one_b, p)? {
            // Integer b: I_x(a,b) = x^a / (a B(a,b)) · ₂F₁(a, 1−b; a+1; x), a polynomial in x
            // (exact inputs give an exact result). Its terms alternate, so large b cancels;
            // then the continued fraction below is used instead.
            let ap1 = self.add(&one, p, RoundingMode::None)?;
            let (f, lost) = hypergeom_series(self, &one_b, &ap1, x, p)?;
            if lost <= WORD_BIT_SIZE {
                return xa
                    .div(self, p, RoundingMode::None)?
                    .div(&beta, p, RoundingMode::None)?
                    .mul(&f, p, RoundingMode::None);
            }
        }
        // I_x(a,b) = x^a (1−x)^b / (a B(a,b)) · CF(a,b,x); the continued fraction converges
        // quickly for x < (a+1)/(a+b+2), which `betainc_at` guarantees.
        let omxb = one
            .sub(x, p, RoundingMode::None)?
            .pow(b, p, RoundingMode::None, cc)?;
        let cf = betainc_cf(self, b, x, p)?;
        xa.mul(&omxb, p, RoundingMode::None)?
            .div(self, p, RoundingMode::None)?
            .div(&beta, p, RoundingMode::None)?
            .mul(&cf, p, RoundingMode::None)
    }
}

/// Continued fraction of the incomplete beta (DLMF 8.17.22, modified Lentz):
/// `1/(1+ d₁/(1+ d₂/(1+ …)))` with
/// `d₂ₘ₊₁ = −(a+m)(a+b+m)x / ((a+2m)(a+2m+1))`, `d₂ₘ = m(b−m)x / ((a+2m−1)(a+2m))`.
fn betainc_cf(
    a: &ExactNumNumber,
    b: &ExactNumNumber,
    x: &ExactNumNumber,
    p: usize,
) -> Result<ExactNumNumber, Error> {
    let one = ExactNumNumber::from_word(1, p)?;
    let make_tiny =
        || ExactNumNumber::from_word(1, p)?.ldexp(-(2 * p as Exponent), p, RoundingMode::None);
    let tiny = make_tiny()?;
    // a fresh 2^(-2p) for each floor (exact, the same value as `tiny`): no copy of it
    let floor_tiny = |v: ExactNumNumber| -> Result<ExactNumNumber, Error> {
        if v.abs()?.cmp(&tiny) < 0 {
            make_tiny()
        } else {
            Ok(v)
        }
    };
    let apb = a.add(b, p, RoundingMode::None)?;
    let ap1 = a.add(&one, p, RoundingMode::None)?;
    let mut c = ExactNumNumber::from_word(1, p)?;
    let d0 = floor_tiny(
        one.sub(
            &apb.mul(x, p, RoundingMode::None)?
                .div(&ap1, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?,
    )?;
    // h starts at d; d is None (it equals h) until the first coefficient updates it: no copy.
    let mut h = one.div(&d0, p, RoundingMode::None)?;
    let mut d: Option<ExactNumNumber> = None;
    let mut m: usize = 1;
    loop {
        let mw = from_usize_p(m, p)?;
        let m2 = from_usize_p(2 * m, p)?;
        let a_m2 = a.add(&m2, p, RoundingMode::None)?;
        let even = mw
            .mul(&b.sub(&mw, p, RoundingMode::None)?, p, RoundingMode::None)?
            .mul(x, p, RoundingMode::None)?
            .div(
                &a_m2
                    .sub(&one, p, RoundingMode::None)?
                    .mul(&a_m2, p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )?;
        let odd = a
            .add(&mw, p, RoundingMode::None)?
            .mul(&apb.add(&mw, p, RoundingMode::None)?, p, RoundingMode::None)?
            .mul(x, p, RoundingMode::None)?
            .div(
                &a_m2.mul(
                    &a_m2.add(&one, p, RoundingMode::None)?,
                    p,
                    RoundingMode::None,
                )?,
                p,
                RoundingMode::None,
            )?
            .neg()?;
        let mut delta = ExactNumNumber::from_word(1, p)?;
        for coef in [even, odd] {
            let dn = floor_tiny(one.add(
                &coef.mul(d.as_ref().unwrap_or(&h), p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )?)?;
            c = floor_tiny(one.add(
                &coef.div(&c, p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )?)?;
            let dn = one.div(&dn, p, RoundingMode::None)?;
            delta = dn.mul(&c, p, RoundingMode::None)?;
            d = Some(dn);
            h = h.mul(&delta, p, RoundingMode::None)?;
        }
        // Stop one word short of p: rounding keeps |Δ−1| near 2^-p.
        let step = delta.sub(&one, p, RoundingMode::None)?;
        if step.is_zero() || step.exponent() + ((p - WORD_BIT_SIZE) as i64) < 0 {
            return Ok(h);
        }
        m = m.checked_add(1).ok_or(Error::PrecisionRetryExhausted)?;
    }
}

/// `None` when `x² ≥ 0.7·p` (the asymptotic `erfc` reaches `2^-p` before diverging);
/// otherwise the extra bits (`≈ x²·log₂e` plus a word) the series paths need.
fn erf_x2_bits(x: &ExactNumNumber, p: usize) -> Result<Option<usize>, Error> {
    if x.exponent() >= 16 {
        return Ok(None);
    }
    let x2 = x.mul(x, p, RoundingMode::None)?.int_as_usize()?;
    if x2 >= p * 7 / 10 + 4 {
        return Ok(None);
    }
    Ok(Some(round_p((x2 + 1) * 3 / 2) + WORD_BIT_SIZE))
}

/// True when `x² ≥ p + 64`, so `erfc|x| < 2^-(p+64)` sits below half an ulp of 1 at `p` bits.
fn erf_tail_below(x: &ExactNumNumber, p: usize) -> bool {
    if x.is_zero() {
        return false;
    }
    let e = x.exponent();
    if e < 2 {
        return false;
    }
    // |x| ≥ 2^(e-1), so x² ≥ 2^(2e-2).
    match u32::try_from(2 * e - 2) {
        Ok(k) => 1usize.checked_shl(k).is_none_or(|v| v >= p + 64),
        Err(_) => true,
    }
}

/// Binary exponent from which `erfc` underflows below every subnormal (`x ≥ 2^31`).
const ERFC_UNDERFLOW_EXPONENT: Exponent = 32;

/// Bernoulli terms for the digamma asymptotic series at precision `p`: `p/8`, at least 64.
/// Combined with [`asymp_min_exponent`], this is the remainder bound that sets
/// both the term count and the Stirling shift target (no fixed term or walk cap).
pub(crate) fn digamma_terms(p: usize) -> usize {
    (p / 8).max(64)
}

/// Stirling term count and minimum exponent for Γ / ln Γ at dest `p`.
/// A fixed 64-term tail (`B_128/(128·127·z^127) ≈ 2^363.8/z^127`) needs
/// about `2^(p/127)` shifts and would hit a 2^20 walk cap from p ≈ 2048.
/// When `digamma_terms(p)` lowers the target, use that count; the extra
/// 64 bits in the bound cover the factor `z/(2K-1)` between a Stirling
/// term and the digamma term [`asymp_min_exponent`] measures. Bound at
/// `p + WORD_BIT_SIZE`, the precision `ln_gamma_stirling` checks its
/// stop rule against.
pub(crate) fn gamma_stirling_plan(p: usize) -> (usize, Exponent) {
    let pw = p + WORD_BIT_SIZE;
    let e64 = stirling_min_exponent(pw, 372, 127, 6);
    let k = digamma_terms(p);
    if k > 64 {
        let ek = asymp_min_exponent(pw + WORD_BIT_SIZE, k);
        if ek < e64 {
            return (k, ek);
        }
    }
    (64, e64)
}

/// Smallest binary exponent `e` (`z ≥ 2^(e-1)`) for which the `k`-th term of a
/// `terms`-term Stirling-type tail, `≈ (k/(πe·z))^{2k}`, is below `2^-(p+8)`:
/// `e = ⌈(p+8)/(2k)⌉ + ⌈log2 k⌉ − 2`, never below 8.
pub(crate) fn asymp_min_exponent(p: usize, terms: usize) -> Exponent {
    let mut lg = 0usize;
    while (1usize << lg) < terms {
        lg += 1;
    }
    let need = (p + 8).div_ceil(2 * terms) + lg;
    Exponent::try_from(need).map_or(Exponent::MAX, |e| (e - 2).max(8))
}

/// Largest `n` for which Γ(n+½) and ψ(n+½) use their closed forms.
/// Beta function `B(a, b)` at `p` bits with `rm`; see [`ExactNumNumber::beta`].
pub fn beta(a: &ExactNum, b: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> ExactNum {
    let a_c = a.exact_copy_ext();
    let b_c = b.exact_copy_ext();
    ExactNum::round_exact(
        a.inexact_accumulated() || b.inexact_accumulated(),
        p,
        rm,
        |q, r| a_c.beta_word(&b_c, q, r, cc),
    )
}

/// The exact integer `lo · (lo+1) ⋯ hi` (1 when `lo > hi`).
/// Proven upper bound on `sum_{k=1}^n floor(log2(k))`. If `2^b ≤ n < 2^{b+1}`,
/// that sum is `sum_{i=1}^{b-1} i 2^i + b (n - 2^b + 1)`, and
/// `sum_{i=1}^{k} i 2^i = (k-1) 2^{k+1} + 2`.
fn sum_floor_log2(n: usize) -> Option<usize> {
    if n <= 1 {
        return Some(0);
    }
    let b = (usize::BITS - n.leading_zeros() - 1) as u128;
    let two_b = 1u128 << b;
    let sum_full = if b >= 2 {
        b.checked_sub(2)?.checked_mul(two_b)?.checked_add(2)
    } else {
        Some(0)
    }?;
    let tail = b.checked_mul((n as u128).checked_sub(two_b)?.checked_add(1)?)?;
    usize::try_from(sum_full.checked_add(tail)?).ok()
}

/// Proven upper bound on the bit length of `m!`. From
/// `log2(m!) = sum_{k=2}^m log2(k) ≤ (m-1) + sum_{k=2}^m floor(log2(k))`,
/// the bit length is at most `m + sum_{k=1}^m floor(log2(k))`.
fn factorial_bit_len_upper(m: usize) -> Option<usize> {
    if m <= 1 {
        return Some(1);
    }
    m.checked_add(sum_floor_log2(m)?)
}

/// `| (a)_n | ≤ (|a| + n − 1)^n`, so the bit length is at most
/// `1 + n (max(0, exp(a), floor(log2(n))) + 1)` with `exp(a)` the exponent
/// convention `|a| < 2^{exp(a)}` (0.mantissa · 2^{exp}).
fn pochhammer_bit_len_upper(a: &ExactNumNumber, n: usize) -> Option<usize> {
    if n == 0 {
        return Some(1);
    }
    let ea = if a.is_zero() { 0i64 } else { a.exponent().max(0) };
    let en = if n <= 1 { 0i64 } else { i64::from(usize::BITS - n.leading_zeros() - 1) };
    let logb = ea.max(en).checked_add(1)?;
    let logb = usize::try_from(logb).ok()?;
    n.checked_mul(logb)?.checked_add(1)
}

/// The exact integer-beta product is the right kernel when the proven bit
/// length of `(n−1)!` and, if `a` is given, of `(a)_n` fits in `p + 2` words.
/// Otherwise the caller uses the Ziv-rounded Γ / ln Γ path (a method switch).
fn exact_beta_product_fits(p: usize, n: usize, a: Option<&ExactNumNumber>) -> bool {
    let budget = p.saturating_add(2 * WORD_BIT_SIZE);
    let Some(fb) = factorial_bit_len_upper(n.saturating_sub(1)) else {
        return false;
    };
    if fb > budget {
        return false;
    }
    if let Some(a) = a {
        let Some(pb) = pochhammer_bit_len_upper(a, n) else {
            return false;
        };
        if pb > budget {
            return false;
        }
    }
    ExactNumNumber::p_assertion(round_p(budget)).is_ok()
}

/// Two-integer B(m, n): the denominator `m (m+1) ⋯ (m+n−1)` has bit length
/// at most that of `(m+n−1)!`.
fn exact_beta_ints_fit(p: usize, m: usize, n: usize) -> bool {
    let budget = p.saturating_add(2 * WORD_BIT_SIZE);
    let Some(end) = m.checked_add(n).and_then(|s| s.checked_sub(1)) else {
        return false;
    };
    let Some(db) = factorial_bit_len_upper(end) else {
        return false;
    };
    db <= budget
}

fn exact_product(lo: usize, hi: usize) -> Result<ExactNumNumber, Error> {
    let mut prod = ExactNumNumber::from_word(1, WORD_BIT_SIZE)?;
    for k in lo..=hi {
        let pb = round_p(prod.mantissa_max_bit_len() + 2 * WORD_BIT_SIZE);
        prod = prod.mul(&from_usize_p(k, pb)?, pb, RoundingMode::None)?;
    }
    Ok(prod)
}

/// Beta arguments up to 2^BETA_RATIO_EXP use the Γ ratio; larger ones the ln Γ sum.
const BETA_RATIO_EXP: Exponent = 10;
/// Fixed walk cap used only by tests that force a cap below the
/// remainder-bound step count. Production Γ and ψ use
/// [`stirling_shift_count`] (P9).
#[cfg(test)]
const STIRLING_SHIFT_MAX: usize = 1 << 20;

/// Number of `z → z+1` steps needed to reach `2^{e_min-1}` from a positive
/// argument less than 1. This is the mathematical length of the walk, not
/// a fixed algorithm cap.
pub(crate) fn stirling_shift_count(e_min: Exponent) -> usize {
    let k = e_min.max(1) - 1;
    if k >= Exponent::from(usize::BITS) {
        usize::MAX
    } else {
        (1usize << k).saturating_add(2)
    }
}

/// Applies `step` (in the callers: `z → z+1` plus their accumulator update) while
/// `z.exponent() < e_min`, i.e. until `z ≥ 2^(e_min−1)`; `z > 0`.
///
/// Each step adds 1, so at least `2^(e_min−1) − z` steps are needed. When that is above `shift_max`
/// it returns `PrecisionRetryExhausted` at once, without doing the steps (before, a high precision
/// ran all `shift_max` multiplications at full width and then failed). Otherwise it stops with the
/// same error after `shift_max` steps, so a step that does not move `z` up cannot loop forever.
/// An `Err` from `step` is returned as is.
///
/// `z0` is borrowed and `step` maps the current `z` to the next one, so the caller needs no owned copy of its
/// argument: the result is `None` when no step was taken (`z` is still `z0`), else the shifted `z`.
fn shift_up<S>(
    z0: &ExactNumNumber,
    e_min: Exponent,
    shift_max: usize,
    p: usize,
    mut step: S,
) -> Result<Option<ExactNumNumber>, Error>
where
    S: FnMut(&ExactNumNumber) -> Result<ExactNumNumber, Error>,
{
    let mut target = ExactNumNumber::from_word(1, p)?;
    target.set_exponent(e_min.min(EXPONENT_MAX));
    let need = target.sub(z0, p, RoundingMode::None)?;
    if need.cmp(&from_usize_p(shift_max, p)?) > 0 {
        return Err(Error::PrecisionRetryExhausted);
    }
    let mut z: Option<ExactNumNumber> = None;
    let mut left = shift_max;
    while z.as_ref().unwrap_or(z0).exponent() < e_min {
        if left == 0 {
            return Err(Error::PrecisionRetryExhausted);
        }
        left -= 1;
        z = Some(step(z.as_ref().unwrap_or(z0))?);
    }
    Ok(z)
}

/// Smallest binary exponent `e` (so `z ≥ 2^(e-1)`) for which a 64-term Stirling-type tail
/// `≈ 2^c / z^slope` (with `c` already including the 8 guard bits) is below `2^-p`;
/// never below `floor`.
pub(crate) fn stirling_min_exponent(p: usize, c: usize, slope: usize, floor: Exponent) -> Exponent {
    let need = (p + c).div_ceil(slope) + 1;
    Exponent::try_from(need).map_or(Exponent::MAX, |e| e.max(floor))
}

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

/// A real 2F1 kernel value and the bits lost to cancellation while forming it.
type Lossy = (ExactNumNumber, usize);

/// Largest exponent among the non-zero values (`Exponent::MIN` when all are zero).
fn max_exp(vals: &[&ExactNumNumber]) -> Exponent {
    vals.iter()
        .filter(|v| !v.is_zero())
        .map(|v| v.exponent())
        .max()
        .unwrap_or(Exponent::MIN)
}

/// Bits of `v` lost when it was formed from parts as large as `2^mag`.
fn lost_bits(mag: Exponent, v: &ExactNumNumber) -> usize {
    if v.is_zero() || mag == Exponent::MIN {
        return 0;
    }
    mag.saturating_sub(v.exponent()).max(0) as usize
}

/// `1/Γ(x)`, zero at the poles of Γ.
fn rgamma(x: &ExactNumNumber, p: usize, cc: &mut Consts) -> Result<ExactNumNumber, Error> {
    if terminating_neg_int(x, p)? {
        return ExactNumNumber::new2(p, Sign::Pos, false);
    }
    ExactNumNumber::from_word(1, p)?.div(&x.gamma_at(p, cc)?, p, RoundingMode::None)
}

/// `max(|v|)` over `vals`.
fn max_abs(vals: &[&ExactNumNumber]) -> Result<ExactNumNumber, Error> {
    let mut m = vals[0].abs()?;
    for v in &vals[1..] {
        m = en_max_owned(m, v.abs()?);
    }
    Ok(m)
}

/// Series tail test for terms `t_{n+1} = factor · t_n`: the term is below `2^-(p+4)` of the
/// sum, `n` is past every parameter (so later factors stay near `|z| ≤ 1/2`), and the current
/// ratio is below 7/8, which bounds the tail by 7 terms of this size.
fn series_done(
    term_exp: Exponent,
    sum: &ExactNumNumber,
    factor: &ExactNumNumber,
    n: &ExactNumNumber,
    max_par: &ExactNumNumber,
    p: usize,
) -> Result<bool, Error> {
    if sum.is_zero() || n.cmp(max_par) <= 0 {
        return Ok(false);
    }
    if term_exp >= sum.exponent() - (p as i64) - 4 {
        return Ok(false);
    }
    let seven_eighths = ExactNumNumber::from_word(7, p)?.ldexp(-3, p, RoundingMode::None)?;
    Ok(factor.abs()?.cmp(&seven_eighths) < 0)
}

/// Gauss series \(\sum (a)_n (b)_n / ((c)_n n!) z^n\). Callers keep \(\lvert z\rvert\le 1/2\)
/// unless the series terminates.
fn hypergeom_series(
    a: &ExactNumNumber,
    b: &ExactNumNumber,
    c: &ExactNumNumber,
    z: &ExactNumNumber,
    p: usize,
) -> Result<Lossy, Error> {
    let one = ExactNumNumber::from_word(1, p)?;
    let mut term = ExactNumNumber::from_word(1, p)?;
    let mut sum = ExactNumNumber::from_word(1, p)?;
    let tiny =
        ExactNumNumber::from_word(1, p)?.ldexp(-((p as Exponent) - 8), p, RoundingMode::None)?;
    let max_par = max_abs(&[a, b, c])?;
    let mut mag = one.exponent();
    let mut n: usize = 0;
    loop {
        let nn = from_usize_p(n, p)?;
        let den = c.add(&nn, p, RoundingMode::None)?.mul(
            &one.add(&nn, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
        if den.abs()?.cmp(&tiny) < 0 {
            return Err(Error::InvalidArgument);
        }
        let factor = a
            .add(&nn, p, RoundingMode::None)?
            .mul(&b.add(&nn, p, RoundingMode::None)?, p, RoundingMode::None)?
            .mul(z, p, RoundingMode::None)?
            .div(&den, p, RoundingMode::None)?;
        term = term.mul(&factor, p, RoundingMode::None)?;
        sum = sum.add(&term, p, RoundingMode::None)?;
        mag = mag.max(max_exp(&[&term]));
        if terminating_after(a, n + 1, p)? || terminating_after(b, n + 1, p)? {
            sum.set_inexact(false);
            let lost = lost_bits(mag, &sum);
            return Ok((sum, lost));
        }
        if series_done(max_exp(&[&term]), &sum, &factor, &nn, &max_par, p)? {
            let lost = lost_bits(mag, &sum);
            return Ok((sum, lost));
        }
        // As n → ∞ the factor tends to z, so |z| < 1 still converges even if terms
        // grow for a while after n > max(|a|, |b|, |c|). |z| ≥ 1 and a factor of
        // magnitude ≥ 1 past that peak is not in the disk.
        if nn.cmp(&max_par) > 0 && factor.abs()?.cmp(&one) >= 0 && z.abs()?.cmp(&one) >= 0 {
            return Err(Error::InvalidArgument);
        }
        n = n.checked_add(1).ok_or(Error::PrecisionRetryExhausted)?;
    }
}

/// Real \({}_2F_1\) dispatch; see [`ExactNumNumber::hypergeom_2f1`].
fn hyp2f1_real(
    a: &ExactNumNumber,
    b: &ExactNumNumber,
    c: &ExactNumNumber,
    z: &ExactNumNumber,
    p: usize,
    cc: &mut Consts,
) -> Result<Lossy, Error> {
    if pole_c_before_term(c, a, b, p)? {
        return Err(Error::InvalidArgument);
    }
    if terminating_neg_int(a, p)? || terminating_neg_int(b, p)? {
        return hypergeom_series(a, b, c, z, p);
    }
    let one = ExactNumNumber::from_word(1, p)?;
    let vs_one = z.cmp(&one);
    if vs_one == 0 {
        return Ok((gauss_z_one(a, b, c, p, cc)?, 0));
    }
    if vs_one > 0 {
        return Err(Error::InvalidArgument);
    }
    let half = one_half(p)?;
    if z.cmp(&half.neg()?) < 0 {
        return hyp2f1_pfaff(a, b, c, z, p, cc);
    }
    if z.cmp(&half) <= 0 {
        return hypergeom_series(a, b, c, z, p);
    }
    hyp2f1_one_minus_z(a, b, c, z, p, cc)
}

/// Pfaff (DLMF 15.8.1): \({}_2F_1(a,b;c;z)=(1-z)^{-a}{}_2F_1(a,c-b;c;z/(z-1))\) for
/// \(z<-1/2\), so \(z/(z-1)\in(1/3,1)\). Uses the \(a\leftrightarrow b\) twin when that
/// makes the new series terminate.
fn hyp2f1_pfaff(
    a: &ExactNumNumber,
    b: &ExactNumNumber,
    c: &ExactNumNumber,
    z: &ExactNumNumber,
    p: usize,
    cc: &mut Consts,
) -> Result<Lossy, Error> {
    let one = ExactNumNumber::from_word(1, p)?;
    let omz = one.sub(z, p, RoundingMode::None)?;
    let w = z.div(&z.sub(&one, p, RoundingMode::None)?, p, RoundingMode::None)?;
    let cb = c.sub(b, p, RoundingMode::None)?;
    let ca = c.sub(a, p, RoundingMode::None)?;
    let (s, t) = if !terminating_neg_int(&cb, p)? && terminating_neg_int(&ca, p)? {
        (b, ca)
    } else {
        (a, cb)
    };
    let (f, lost) = hyp2f1_real(s, &t, c, &w, p, cc)?;
    let pref = omz.pow(&s.neg()?, p, RoundingMode::None, cc)?;
    Ok((pref.mul(&f, p, RoundingMode::None)?, lost))
}

/// \(1-z\) connection (DLMF 15.8.4) for \(z\in(1/2,1)\), non-terminating \(a,b\):
/// \(\frac{\Gamma(c)\Gamma(m)}{\Gamma(c-a)\Gamma(c-b)}F(a,b;1-m;\zeta)
/// +\zeta^{m}\frac{\Gamma(c)\Gamma(-m)}{\Gamma(a)\Gamma(b)}F(c-a,c-b;1+m;\zeta)\),
/// \(m=c-a-b\), \(\zeta=1-z\). Integer \(m\) goes to the logarithmic form.
fn hyp2f1_one_minus_z(
    a: &ExactNumNumber,
    b: &ExactNumNumber,
    c: &ExactNumNumber,
    z: &ExactNumNumber,
    p: usize,
    cc: &mut Consts,
) -> Result<Lossy, Error> {
    let one = ExactNumNumber::from_word(1, p)?;
    let zeta = one.sub(z, p, RoundingMode::None)?;
    let ca = c.sub(a, p, RoundingMode::None)?;
    let cb = c.sub(b, p, RoundingMode::None)?;
    let m = ca.sub(b, p, RoundingMode::None)?;
    if let Some(mi) = as_i32_exact(&m, p)? {
        if mi < 0 {
            // Euler (DLMF 15.8.1): F(a,b;c;z) = ζ^m F(c−a,c−b;c;z), whose c−a−b is −m > 0.
            let (f, lost) = hyp2f1_real(&ca, &cb, c, z, p, cc)?;
            let pref = zeta.pow(&m, p, RoundingMode::None, cc)?;
            return Ok((pref.mul(&f, p, RoundingMode::None)?, lost));
        }
        return hyp2f1_log_case(a, b, c, &zeta, mi as u32, p, cc);
    }
    let gc = c.gamma_at(p, cc)?;
    let r1 = rgamma(&ca, p, cc)?.mul(&rgamma(&cb, p, cc)?, p, RoundingMode::None)?;
    let (t1, l1) = if r1.is_zero() {
        (r1, 0)
    } else {
        let (f, l) = hypergeom_series(a, b, &one.sub(&m, p, RoundingMode::None)?, &zeta, p)?;
        let t = gc
            .mul(&m.gamma_at(p, cc)?, p, RoundingMode::None)?
            .mul(&r1, p, RoundingMode::None)?
            .mul(&f, p, RoundingMode::None)?;
        (t, l)
    };
    let r2 = rgamma(a, p, cc)?.mul(&rgamma(b, p, cc)?, p, RoundingMode::None)?;
    let (t2, l2) = if r2.is_zero() {
        (r2, 0)
    } else {
        let (f, l) = hypergeom_series(&ca, &cb, &one.add(&m, p, RoundingMode::None)?, &zeta, p)?;
        let t = gc
            .mul(&m.neg()?.gamma_at(p, cc)?, p, RoundingMode::None)?
            .mul(&r2, p, RoundingMode::None)?
            .mul(
                &zeta.pow(&m, p, RoundingMode::None, cc)?,
                p,
                RoundingMode::None,
            )?
            .mul(&f, p, RoundingMode::None)?;
        (t, l)
    };
    let v = t1.add(&t2, p, RoundingMode::None)?;
    let lost = l1.max(l2).max(lost_bits(max_exp(&[&t1, &t2]), &v));
    Ok((v, lost))
}

/// Logarithmic case of the \(1-z\) connection, integer \(m=c-a-b\ge 0\), \(\zeta=1-z\)
/// (A&S 15.3.10–15.3.11, DLMF 15.8.10):
/// \(F=\frac{\Gamma(m)\Gamma(c)}{\Gamma(a+m)\Gamma(b+m)}\sum_{n<m}\frac{(a)_n(b)_n}{n!(1-m)_n}\zeta^n
/// -(-\zeta)^m\frac{\Gamma(c)}{\Gamma(a)\Gamma(b)}\sum_{n\ge0}\frac{(a+m)_n(b+m)_n}{n!(n+m)!}\zeta^n
/// [\ln\zeta-\psi(n+1)-\psi(n+m+1)+\psi(a+m+n)+\psi(b+m+n)]\).
fn hyp2f1_log_case(
    a: &ExactNumNumber,
    b: &ExactNumNumber,
    c: &ExactNumNumber,
    zeta: &ExactNumNumber,
    m: u32,
    p: usize,
    cc: &mut Consts,
) -> Result<Lossy, Error> {
    let one = ExactNumNumber::from_word(1, p)?;
    let mw = ExactNumNumber::from_word(m as Word, p)?;
    let am = a.add(&mw, p, RoundingMode::None)?;
    let bm = b.add(&mw, p, RoundingMode::None)?;
    let gc = c.gamma_at(p, cc)?;

    let mut s1 = ExactNumNumber::new2(p, Sign::Pos, false)?;
    let mut l1 = 0;
    if m > 0 {
        let mut t = ExactNumNumber::from_word(1, p)?;
        let mut mag = Exponent::MIN;
        for n in 0..m {
            s1 = s1.add(&t, p, RoundingMode::None)?;
            mag = mag.max(max_exp(&[&t]));
            if n + 1 < m {
                let nw = ExactNumNumber::from_word(n as Word, p)?;
                let n1 = ExactNumNumber::from_word((n + 1) as Word, p)?;
                let num = a
                    .add(&nw, p, RoundingMode::None)?
                    .mul(&b.add(&nw, p, RoundingMode::None)?, p, RoundingMode::None)?
                    .mul(zeta, p, RoundingMode::None)?;
                let den = n1.mul(&n1.sub(&mw, p, RoundingMode::None)?, p, RoundingMode::None)?;
                t = t
                    .mul(&num, p, RoundingMode::None)?
                    .div(&den, p, RoundingMode::None)?;
            }
        }
        l1 = lost_bits(mag, &s1);
        let pref = mw
            .gamma_at(p, cc)?
            .mul(&gc, p, RoundingMode::None)?
            .mul(&rgamma(&am, p, cc)?, p, RoundingMode::None)?
            .mul(&rgamma(&bm, p, cc)?, p, RoundingMode::None)?;
        s1 = pref.mul(&s1, p, RoundingMode::None)?;
    }

    let mut coef = zeta
        .powi(m as usize, p, RoundingMode::None)?
        .mul(&gc, p, RoundingMode::None)?
        .mul(&rgamma(a, p, cc)?, p, RoundingMode::None)?
        .mul(&rgamma(b, p, cc)?, p, RoundingMode::None)?;
    if coef.is_zero() {
        return Ok((s1, l1));
    }
    if m % 2 == 1 {
        coef = coef.neg()?;
    }

    let ln_z = zeta.ln(p, RoundingMode::None, cc)?;
    let mut psi_n1 = one.digamma_at(p, cc)?;
    let mut psi_nm1 = mw.add(&one, p, RoundingMode::None)?.digamma_at(p, cc)?;
    let mut psi_a = am.digamma_at(p, cc)?;
    let mut psi_b = bm.digamma_at(p, cc)?;
    let mut t = one.div(&factorial(m as usize, p)?, p, RoundingMode::None)?;
    let mut s2 = ExactNumNumber::new2(p, Sign::Pos, false)?;
    let max_par = max_abs(&[&am, &bm, &mw])?;
    let mut mag = Exponent::MIN;
    let mut n: usize = 0;
    loop {
        let br = ln_z
            .sub(&psi_n1, p, RoundingMode::None)?
            .sub(&psi_nm1, p, RoundingMode::None)?
            .add(&psi_a, p, RoundingMode::None)?
            .add(&psi_b, p, RoundingMode::None)?;
        s2 = s2.add(&t.mul(&br, p, RoundingMode::None)?, p, RoundingMode::None)?;
        let part =
            max_exp(&[&t]).saturating_add(max_exp(&[&ln_z, &psi_n1, &psi_nm1, &psi_a, &psi_b]));
        mag = mag.max(part);
        let nw = from_usize_p(n, p)?;
        let n1 = from_usize_p(n + 1, p)?;
        let nm1 = from_usize_p(n + (m as usize) + 1, p)?;
        let an = am.add(&nw, p, RoundingMode::None)?;
        let bn = bm.add(&nw, p, RoundingMode::None)?;
        let factor = an
            .mul(&bn, p, RoundingMode::None)?
            .mul(zeta, p, RoundingMode::None)?
            .div(&n1.mul(&nm1, p, RoundingMode::None)?, p, RoundingMode::None)?;
        if series_done(part, &s2, &factor, &nw, &max_par, p)? {
            let l2 = lost_bits(mag, &s2);
            let t2 = coef.mul(&s2, p, RoundingMode::None)?;
            let v = s1.sub(&t2, p, RoundingMode::None)?;
            let lost = l1.max(l2).max(lost_bits(max_exp(&[&s1, &t2]), &v));
            return Ok((v, lost));
        }
        if nw.cmp(&max_par) > 0 && factor.abs()?.cmp(&one) >= 0 && zeta.abs()?.cmp(&one) >= 0 {
            return Err(Error::InvalidArgument);
        }
        t = t.mul(&factor, p, RoundingMode::None)?;
        psi_n1 = psi_n1.add(&one.div(&n1, p, RoundingMode::None)?, p, RoundingMode::None)?;
        psi_nm1 = psi_nm1.add(
            &one.div(&nm1, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
        psi_a = psi_a.add(&one.div(&an, p, RoundingMode::None)?, p, RoundingMode::None)?;
        psi_b = psi_b.add(&one.div(&bn, p, RoundingMode::None)?, p, RoundingMode::None)?;
        n = n.checked_add(1).ok_or(Error::PrecisionRetryExhausted)?;
    }
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
    if !cab.is_positive() || cab.is_zero() {
        return Err(Error::InvalidArgument);
    }
    let gc = c.gamma_at(p, cc)?;
    let gcab = cab.gamma_at(p, cc)?;
    let rca = rgamma(&c.sub(a, p, RoundingMode::None)?, p, cc)?;
    let rcb = rgamma(&c.sub(b, p, RoundingMode::None)?, p, cc)?;
    gc.mul(&gcab, p, RoundingMode::None)?
        .mul(&rca, p, RoundingMode::None)?
        .mul(&rcb, p, RoundingMode::None)
}

fn tiny_spread(p: usize) -> Result<ExactNumNumber, Error> {
    ExactNumNumber::from_word(1, p)?.ldexp(-((p as Exponent) / 3 + 16), p, RoundingMode::None)
}

/// The larger of `a` and `b` (`a` on a tie), borrowed.
fn en_max<'a>(a: &'a ExactNumNumber, b: &'a ExactNumNumber) -> &'a ExactNumNumber {
    if a.cmp(b) >= 0 {
        a
    } else {
        b
    }
}

/// The larger of `a` and `b` (`a` on a tie), moved out.
fn en_max_owned(a: ExactNumNumber, b: ExactNumNumber) -> ExactNumNumber {
    if a.cmp(&b) >= 0 {
        a
    } else {
        b
    }
}

fn close_enough(dev: &ExactNumNumber, an: &ExactNumNumber, p: usize) -> Result<bool, Error> {
    let one = ExactNumNumber::from_word(1, p)?;
    let a = an.abs()?;
    let scale = en_max(&a, &one);
    let thresh = tiny_spread(p)?.mul(scale, p, RoundingMode::None)?;
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
    Ok(en_max_owned(en_max_owned(dx, dy), dz))
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
    // RC(x, x(1+t)) = x^{-1/2}(1 − t/3 + …): the shortcut 1/√x is only exact to 2^-p when
    // |x − y| < 2^-(p+8)·max(x, y). (It used the absolute 2^-(p/3+16) before, which capped R_J —
    // and so Π — at about 340 bits for p = 512; the atan/atanh forms below do not cancel.)
    let rel = ExactNumNumber::from_word(1, p)?
        .ldexp(-((p as Exponent) + 8), p, RoundingMode::None)?
        .mul(en_max(x, y), p, RoundingMode::None)?;
    if x.sub(y, p, RoundingMode::None)?.abs()?.cmp(&rel) < 0 {
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
    // The iterates start at the arguments: None stands for x0, y0, z0 themselves (no copies).
    let mut x_s: Option<ExactNumNumber> = None;
    let mut y_s: Option<ExactNumNumber> = None;
    let mut z_s: Option<ExactNumNumber> = None;
    loop {
        let x = x_s.as_ref().unwrap_or(x0);
        let y = y_s.as_ref().unwrap_or(y0);
        let z = z_s.as_ref().unwrap_or(z0);
        let an = x
            .add(y, p, RoundingMode::None)?
            .add(z, p, RoundingMode::None)?
            .div(&three, p, RoundingMode::None)?;
        if close_enough(&max_dev3(&an, x, y, z, p)?, &an, p)? {
            return rf_series(&an, x, y, z, p);
        }
        let sx = x.sqrt(p, RoundingMode::None)?;
        let sy = y.sqrt(p, RoundingMode::None)?;
        let sz = z.sqrt(p, RoundingMode::None)?;
        let lam = sx
            .mul(&sy, p, RoundingMode::None)?
            .add(&sy.mul(&sz, p, RoundingMode::None)?, p, RoundingMode::None)?
            .add(&sz.mul(&sx, p, RoundingMode::None)?, p, RoundingMode::None)?;
        x_s = Some(
            x.add(&lam, p, RoundingMode::None)?
                .div(&four, p, RoundingMode::None)?,
        );
        y_s = Some(
            y.add(&lam, p, RoundingMode::None)?
                .div(&four, p, RoundingMode::None)?,
        );
        z_s = Some(
            z.add(&lam, p, RoundingMode::None)?
                .div(&four, p, RoundingMode::None)?,
        );
    }
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
    // The iterates start at the arguments: None stands for x0, y0, z0 themselves (no copies).
    let mut x_s: Option<ExactNumNumber> = None;
    let mut y_s: Option<ExactNumNumber> = None;
    let mut z_s: Option<ExactNumNumber> = None;
    let mut sum = ExactNumNumber::from_word(0, p)?;
    let mut fac = ExactNumNumber::from_word(1, p)?;
    loop {
        let x = x_s.as_ref().unwrap_or(x0);
        let y = y_s.as_ref().unwrap_or(y0);
        let z = z_s.as_ref().unwrap_or(z0);
        let an = x
            .add(y, p, RoundingMode::None)?
            .add(&three.mul(z, p, RoundingMode::None)?, p, RoundingMode::None)?
            .div(&five, p, RoundingMode::None)?;
        if close_enough(&max_dev3(&an, x, y, z, p)?, &an, p)? {
            let series = rd_series(&an, x, y, z, p)?;
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
        x_s = Some(
            x.add(&lam, p, RoundingMode::None)?
                .div(&four, p, RoundingMode::None)?,
        );
        y_s = Some(
            y.add(&lam, p, RoundingMode::None)?
                .div(&four, p, RoundingMode::None)?,
        );
        z_s = Some(
            z.add(&lam, p, RoundingMode::None)?
                .div(&four, p, RoundingMode::None)?,
        );
    }
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
    // The iterates start at the arguments: None stands for x0, y0, z0, p0 themselves (no copies).
    let mut x_s: Option<ExactNumNumber> = None;
    let mut y_s: Option<ExactNumNumber> = None;
    let mut z_s: Option<ExactNumNumber> = None;
    let mut pv_s: Option<ExactNumNumber> = None;
    let mut sum = ExactNumNumber::from_word(0, p)?;
    let mut fac = ExactNumNumber::from_word(1, p)?;
    loop {
        let x = x_s.as_ref().unwrap_or(x0);
        let y = y_s.as_ref().unwrap_or(y0);
        let z = z_s.as_ref().unwrap_or(z0);
        let pv = pv_s.as_ref().unwrap_or(p0);
        let an = x
            .add(y, p, RoundingMode::None)?
            .add(z, p, RoundingMode::None)?
            .add(&two.mul(pv, p, RoundingMode::None)?, p, RoundingMode::None)?
            .div(&five, p, RoundingMode::None)?;
        let d4 = max_dev3(&an, x, y, z, p)?;
        let dp = an.sub(pv, p, RoundingMode::None)?.abs()?;
        if close_enough(en_max(&d4, &dp), &an, p)? {
            let series = rj_series(&an, x, y, z, pv, p)?;
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
        x_s = Some(
            x.add(&lam, p, RoundingMode::None)?
                .div(&four, p, RoundingMode::None)?,
        );
        y_s = Some(
            y.add(&lam, p, RoundingMode::None)?
                .div(&four, p, RoundingMode::None)?,
        );
        z_s = Some(
            z.add(&lam, p, RoundingMode::None)?
                .div(&four, p, RoundingMode::None)?,
        );
        pv_s = Some(
            pv.add(&lam, p, RoundingMode::None)?
                .div(&four, p, RoundingMode::None)?,
        );
    }
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

/// `|t| < 2^-p |sum|` (or `t = 0`): relative stopping rule for series whose sum can be far
/// below 1 (e.g. `J_ν(x) ≈ (x/2)^ν/Γ(ν+1)` for small `x`, `Si(x) ≈ x`). The absolute rule
/// `|t| < 2^-p` used before truncated such sums early (J_10.5(1e-10) had 367 correct bits at p=512).
fn series_negligible(t: &ExactNumNumber, sum: &ExactNumNumber, p: usize) -> bool {
    if t.is_zero() {
        return true;
    }
    let se = if sum.is_zero() { 0 } else { sum.exponent() };
    t.exponent() + (p as i64) < se
}

// Cancellation guard for a Bessel series whose terms grow like e^{c|x|}:
// ceil(num/den * |x|) + 16, or MemoryAllocation if that width cannot be formed.

/// `|x| ≥ max(16, 0.35(dest_p+112))` — same switch as the complex Hankel path.
fn hankel_mag_ok(x: &ExactNumNumber, dest_p: usize) -> bool {
    if x.is_zero() || x.exponent() < 5 {
        return false;
    }
    let t = dest_p.saturating_add(112).saturating_mul(35) / 100;
    if x.exponent() > 30 {
        return true;
    }
    match x
        .abs()
        .and_then(|a| a.ceil())
        .and_then(|c| c.int_as_usize())
    {
        Ok(c) => c >= t.max(16),
        Err(_) => true,
    }
}

/// First Hankel term decreases when `|ν|² < 2|x|`.
fn hankel_order_ok(x: &ExactNumNumber, nu: &ExactNumNumber) -> bool {
    let p = WORD_BIT_SIZE;
    let Ok(two) = ExactNumNumber::from_word(2, p) else {
        return false;
    };
    let Ok(ax) = x.abs() else {
        return false;
    };
    let Ok(an) = nu.abs() else {
        return false;
    };
    let Ok(nu2) = an.mul(&an, p, RoundingMode::None) else {
        return false;
    };
    let Ok(two_x) = two.mul(&ax, p, RoundingMode::None) else {
        return false;
    };
    nu2.cmp(&two_x) < 0
}

fn use_bessel_hankel_real(x: &ExactNumNumber, nu: &ExactNumNumber, dest_p: usize) -> bool {
    if x.is_zero() || !hankel_order_ok(x, nu) {
        return false;
    }
    hankel_mag_ok(x, dest_p)
}

fn from_usize_p(n: usize, p: usize) -> Result<ExactNumNumber, Error> {
    if n <= Word::MAX as usize {
        return ExactNumNumber::from_word(n as Word, p);
    }
    let pw = p.max(128);
    let mut v = ExactNumNumber::from_u128(n as u128, pw)?;
    if pw != p {
        v.set_precision(p, RoundingMode::None)?;
    }
    Ok(v)
}

/// `⌈num/den · |x|⌉ + 16` guard bits for a series whose terms grow like `e^{c|x|}`
/// (`num/den ≥ c·log2(e)`); 16 for `|x| < 1`.
fn bessel_cancel_bits(x: &ExactNumNumber, num: usize, den: usize) -> Result<usize, Error> {
    if x.is_zero() || x.exponent() <= 0 {
        return Ok(16);
    }
    let c = match x.abs()?.ceil()?.int_as_usize() {
        Ok(c) => c,
        Err(Error::InvalidArgument) => return Err(Error::MemoryAllocation),
        Err(e) => return Err(e),
    };
    c.saturating_mul(num)
        .checked_div(den)
        .and_then(|v| v.checked_add(16))
        .ok_or(Error::MemoryAllocation)
}

/// `f(−|x|) = (−1)^ν f(|x|)` for an integer order `ν` (J and I, DLMF 10.4.1, 10.27.1): a flipped
/// sign rounds `f(|x|)` in the mirrored direction.
fn bessel_int_order_neg_x(
    x: &ExactNumNumber,
    nu: &ExactNumNumber,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
    f: fn(
        &ExactNumNumber,
        &ExactNumNumber,
        usize,
        RoundingMode,
        &mut Consts,
    ) -> Result<ExactNumNumber, Error>,
) -> Result<ExactNumNumber, Error> {
    let pe = round_p(nu.mantissa_max_bit_len().max(WORD_BIT_SIZE) + WORD_BIT_SIZE);
    let two = ExactNumNumber::from_word(2, WORD_BIT_SIZE)?;
    let odd = !nu.is_zero() && !nu.div(&two, pe, RoundingMode::None)?.is_int();
    let ax = x.abs()?;
    if odd {
        match f(&ax, nu, p, crate::common::util::invert_rm_for_sign(rm), cc) {
            Ok(v) => v.neg(),
            Err(Error::ExponentOverflow(s)) => Err(Error::ExponentOverflow(s.invert())),
            Err(e) => Err(e),
        }
    } else {
        f(&ax, nu, p, rm, cc)
    }
}

/// `±max` for a finite overflow when `rm` rounds toward zero for its sign, else `ExponentOverflow`
/// (C99 7.12.1 / IEEE 754 7.4: a genuine overflow, not a pole, follows the rounding direction).
fn overflow_in_rm(s: Sign, p: usize, rm: RoundingMode) -> Result<ExactNumNumber, Error> {
    let toward_zero = matches!(
        (s.is_positive(), rm),
        (_, RoundingMode::ToZero) | (true, RoundingMode::Down) | (false, RoundingMode::Up)
    );
    if !toward_zero {
        return Err(Error::ExponentOverflow(s));
    }
    let mut m = ExactNumNumber::max_value(p)?;
    m.set_precision(p, RoundingMode::ToZero)?;
    m.set_sign(s);
    m.set_inexact(true);
    Ok(m)
}

/// `v · 2^n` rounded once in `rm` for any `n` (saturated far outside the range: the value is
/// then below the subnormal grid by far more than `p` bits, or an overflow).
fn scale_in_rm(
    v: ExactNumNumber,
    n: i128,
    p: usize,
    rm: RoundingMode,
) -> Result<ExactNumNumber, Error> {
    let lo = 2 * (EXPONENT_MIN as i128);
    if n > EXPONENT_MAX as i128 + 8 {
        return overflow_in_rm(v.sign(), p, rm);
    }
    let s = v.sign();
    match v.ldexp(n.max(lo) as Exponent, p, rm) {
        Err(Error::ExponentOverflow(_)) => overflow_in_rm(s, p, rm),
        r => r,
    }
}

/// ZF-37: `J_ν(x)` / `I_ν(x)` (`i`) as `L · S` with the leading term
/// `L = (x/2)^ν / Γ(ν+1)` (DLMF 10.2.2, 10.25.2) held as `g · 2^N` (`N` unbounded) and
/// `S = Σ_k t_k`, `t_0 = 1`, `t_k = t_{k−1} (∓x²/4) / (k (k+ν))`. Used when `|t_k/t_{k−1}|`
/// is below `r ≤ 2^lr ≤ 1/8` and either `L` is near or past the exponent range (the series
/// in `x` then underflows to an exact zero sum, which certifies 0 bits, or overflows with
/// the wrong sign), or `r` is below the precision (the correction is invisible: a dyadic
/// `L` such as `I_0(2^-100000) = 1 + 2^-200002` would never round in Ziv). `None` otherwise.
///
/// Bounds: `k |k+ν| ≥ b` (`b = 1+ν` for `ν ≥ 0`, else the distance of `ν` to the integers),
/// so `|t_k| ≤ r^k`, the tail after term `K` is below `|t_K| / 7`, and `S ∈ [6/7, 8/7]`.
/// Each `t_k` carries at most `4k` ulps; with the additions `|S − Ŝ| ≤ (K+5) 2^-w`.
/// `ln|x|`, `ln|Γ(ν+1)|` and `ln 2` are faithful at `q` bits, so `Λ = log₂|L|` is good to
/// `2^(eΛ − q + 4)`; `q` is raised until that is below `2^-(w+12)`. When `r` is negligible and
/// `L` is exact (`ν ∈ {0, 1, 2}`: `1`, `x/2`, `x²/8`), the value is `L(1+δ)` with `δ` of the
/// sign of `t_1` and `|δ| < 2^-(W+2)`; it rounds as `L ± 2^(e(L)−W−1)` (no rounding
/// boundary of `W ≥ p+2` bits lies strictly between them).
/// Cheap estimate of the exponent `n` in `L = g · 2^n = (x/2)^ν / Γ(ν+1)`.
/// `a` is ν (signed when non-integer). ν and x are first rounded to 4 words so an
/// 8192-bit order does not run `ln_gamma` at working precision just to discover
/// that `L` is a normal (J/I at x = 2^{-30}, p = 8192, ν ≈ 1/3 were 1.36 s vs
/// 0.07 s; K at the same point, which calls `I_ν` and `I_{-ν}`, was 5.8 s vs 0.29 s).
///
/// The estimate is accepted only away from poles of Γ(ν+1). Then, with
/// `|ν| < 2^{16}` and `x` in the representable range, the ToEven 4-word
/// rounding plus a 4-word `ln_gamma` shift `Λ` by well under one unit
/// (`|δν| ≤ 2^{16-256}`, `|dΛ/dν|` is moderate), so `|n_est − n| ≤ 4`.
/// The skip uses a `2^{20}` margin around the already-padded `lo`/`hi`
/// (`EMIN + 2p + 256`, `EMAX − 256`), which cannot shrink enough for a
/// wrong skip. Near a pole, or if `ln_gamma` overflows, the estimate is
/// refused (`None`) and the precise path runs. Rounding a high-bit ν onto
/// a non-positive integer must not leak `ExponentOverflow` / `DivisionByZero`
/// from this helper — those would turn a finite `J_ν`/`I_ν` into a NaN.
fn estimate_log2_leading_ji(
    x: &ExactNumNumber,
    a: &ExactNumNumber,
    cc: &mut Consts,
) -> Result<Option<i128>, Error> {
    type N = ExactNumNumber;
    let q = 4 * WORD_BIT_SIZE;
    let mut aq = a.clone()?;
    aq.set_precision(q, RoundingMode::ToEven)?;
    aq.set_inexact(false);
    let mut xq = x.abs()?;
    xq.set_precision(q, RoundingMode::ToEven)?;
    xq.set_inexact(false);
    let one = N::from_word(1, q)?;
    let ap1 = aq.add(&one, q, RoundingMode::None)?;
    if gamma_arg_too_close_to_nonpos_int(&ap1, q)? {
        return Ok(None);
    }
    let ln2 = cc.ln_2_num(q, RoundingMode::None)?;
    let la = xq.ln(q, RoundingMode::None, cc)?;
    let al = aq.mul(&la, q, RoundingMode::None)?;
    let lg = match ap1.ln_gamma(q, RoundingMode::None, cc) {
        Ok(v) => v,
        Err(Error::ExponentOverflow(_)) | Err(Error::DivisionByZero) => return Ok(None),
        Err(e) => return Err(e),
    };
    if (lg.exponent() as i128).unsigned_abs() > 60 {
        return Ok(None);
    }
    let d = al.sub(&lg, q, RoundingMode::None)?;
    let lam = d
        .div(&ln2, q, RoundingMode::None)?
        .sub(&aq, q, RoundingMode::None)?;
    if lam.exponent() > 62 {
        let n = if lam.is_negative() { -(1i128 << 100) } else { 1i128 << 100 };
        return Ok(Some(n));
    }
    let fl = lam.floor_ref()?;
    let nf = if fl.is_zero() {
        0
    } else {
        match fl.abs()?.int_as_usize() {
            Ok(u) => u as i128,
            Err(_) => return Ok(None),
        }
    };
    Ok(Some(if fl.is_negative() { -nf } else { nf }))
}

/// `true` if `z` is within `2^{-80}` of an integer `≤ 0` (a pole of `Γ`).
fn gamma_arg_too_close_to_nonpos_int(z: &ExactNumNumber, q: usize) -> Result<bool, Error> {
    type N = ExactNumNumber;
    if z.is_zero() {
        return Ok(true);
    }
    let fl = z.floor_ref()?;
    let one = N::from_word(1, q)?;
    let fl_up;
    let ceil: &N = if z.is_int() {
        &fl
    } else {
        fl_up = fl.add(&one, q, RoundingMode::None)?;
        &fl_up
    };
    let d_fl = z.sub(&fl, q, RoundingMode::None)?.abs()?;
    let d_ce = ceil.sub(z, q, RoundingMode::None)?.abs()?;
    let (d, n) = if d_fl.cmp(&d_ce) <= 0 { (d_fl, &fl) } else { (d_ce, ceil) };
    let close = d.is_zero() || d.exponent() <= -80;
    Ok(close && !n.is_positive())
}

fn bessel_ji_scaled(
    x: &ExactNumNumber,
    nu: &ExactNumNumber,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
    i: bool,
) -> Result<Option<ExactNumNumber>, Error> {
    type N = ExactNumNumber;
    let int = nu.is_int();
    let mut a = if int { nu.abs()? } else { nu.clone()? };
    a.set_inexact(false);
    let pa = round_p(a.mantissa_max_bit_len().max(WORD_BIT_SIZE) + 2 * WORD_BIT_SIZE);
    let two = N::from_word(2, WORD_BIT_SIZE)?;
    let one = N::from_word(1, WORD_BIT_SIZE)?;
    let is_odd = |v: &N| -> Result<bool, Error> {
        Ok(!v.is_zero() && !v.div(&two, pa, RoundingMode::None)?.is_int())
    };
    // Sign of L (and of the result: S > 0).
    let neg = if int {
        let odd = is_odd(&a)?;
        (x.is_negative() && odd) != (!i && nu.is_negative() && odd)
    } else if a.is_negative() {
        // sgn Γ(ν+1) = (−1)^m, m = ⌊−ν⌋.
        is_odd(&a.abs()?.floor()?)?
    } else {
        false
    };
    let (mut m, ex) = x.abs()?.frexp()?;
    m.set_inexact(false);
    let b = if a.is_negative() {
        let f = a.sub(&a.floor_ref()?, pa, RoundingMode::None)?;
        let g = one.sub(&f, pa, RoundingMode::None)?;
        if f.cmp(&g) < 0 {
            f
        } else {
            g
        }
    } else {
        a.add(&one, WORD_BIT_SIZE, RoundingMode::Down)?
    };
    let lr = 2 * ex as i128 - 2 - (b.exponent() as i128 - 1);
    if lr > -3 {
        return Ok(None);
    }
    let xb = x.mantissa_max_bit_len();
    let w_big = round_p(p.max(2 * xb + WORD_BIT_SIZE) + 3 * WORD_BIT_SIZE);
    let tiny = lr <= -(w_big as i128 + 8);
    // Integer order of moderate size: L = m^a 2^(a·ex) / (2^a a!) from exact m^a and a!, so a
    // dyadic L (J_3(3·2^-2000) = −(9/16)·2^-6000 (1 − δ)) is recognised as exact.
    let mb = m.mantissa_max_bit_len();
    let small_int = if int && a.exponent() <= 16 {
        let k = if a.is_zero() { 0 } else { a.int_as_usize()? };
        (k.saturating_mul(mb) <= 1 << 16).then_some(k)
    } else {
        None
    };

    // Series path is preferred when L is a comfortable normal and we are not
    // on the tiny-x exact-leading path. Decide from a short estimate of n so
    // we do not compute ln_gamma(ν+1) at p bits first. `margin` is 2^{20};
    // a accepted estimate is accurate to a few units (see
    // [`estimate_log2_leading_ji`]), so this cannot skip when the true `n`
    // is outside `[lo, hi]`.
    if small_int.is_none() && !tiny && a.exponent() <= 16 {
        let lo = EXPONENT_MIN as i128 + 2 * p as i128 + 256;
        let hi = EXPONENT_MAX as i128 - 256;
        let margin = 1i128 << 20;
        if let Some(n_est) = estimate_log2_leading_ji(x, &a, cc)? {
            if n_est > lo + margin && n_est < hi - margin {
                return Ok(None);
            }
        }
    }

    let mut w = round_p(p + 2 * WORD_BIT_SIZE);
    let mut w_inc = WORD_BIT_SIZE;
    let mut first = true;
    loop {
        // g · 2^n = |L|; `g_exact` when g carries no rounding error.
        let (g, n, g_exact) = if let Some(k) = small_int {
            let pm = round_p(k * mb + WORD_BIT_SIZE);
            let mut mk = N::from_word(1, pm)?;
            for _ in 0..k {
                mk = mk.mul(&m, pm, RoundingMode::None)?;
            }
            let bl = (usize::BITS - k.leading_zeros()) as usize;
            let fk = factorial(k, round_p(k * bl + WORD_BIT_SIZE))?;
            let pq = pm.max(w + WORD_BIT_SIZE);
            let mut g = mk.div(&fk, pq, RoundingMode::None)?;
            let exact = !g.inexact();
            g.set_inexact(false);
            let k = k as i128;
            (g, k * ex as i128 - k, exact)
        } else {
            let mut q = w + 4 * WORD_BIT_SIZE;
            let lam = loop {
                let ln2 = cc.ln_2_num(q, RoundingMode::None)?;
                let la = x.abs()?.ln(q, RoundingMode::None, cc)?;
                let al = a.mul(&la, q, RoundingMode::None)?;
                let ap1 = a.add(&one, pa + q, RoundingMode::None)?;
                let lg = ap1.ln_gamma(q, RoundingMode::None, cc)?;
                let d = al.sub(&lg, q, RoundingMode::None)?;
                let lam = d
                    .div(&ln2, q, RoundingMode::None)?
                    .sub(&a, q, RoundingMode::None)?;
                let e_max = [al.exponent(), lg.exponent(), a.exponent(), lam.exponent(), 1]
                    .into_iter()
                    .max()
                    .unwrap_or(1)
                    .max(0) as usize;
                if e_max + 4 + w + 12 <= q {
                    break lam;
                }
                q = round_p(w + e_max + 4 * WORD_BIT_SIZE);
            };
            if lam.exponent() > 62 {
                // |Λ| ≥ 2^62: certainly outside the range whatever g is.
                let n = if lam.is_negative() { -(1i128 << 100) } else { 1i128 << 100 };
                (N::from_word(1, WORD_BIT_SIZE)?, n, false)
            } else {
                let fl = lam.floor_ref()?;
                let fr = lam.sub(&fl, q, RoundingMode::None)?;
                let nf = if fl.is_zero() { 0 } else { fl.abs()?.int_as_usize()? as i128 };
                let nf = if fl.is_negative() { -nf } else { nf };
                let ln2 = cc.ln_2_num(w + WORD_BIT_SIZE, RoundingMode::None)?;
                let g = fr.mul(&ln2, w + WORD_BIT_SIZE, RoundingMode::None)?.exp(
                    w + WORD_BIT_SIZE,
                    RoundingMode::None,
                    cc,
                )?;
                (g, nf, false)
            }
        };
        // Large orders (|ν| ≥ 2^16) always come here: the x-series path forms Γ(ν+1) and
        // (x/2)^ν directly (J_{2^28}(2^-40) at 64 bits ran for minutes).
        if first && !tiny && a.exponent() <= 16 {
            first = false;
            let lo = EXPONENT_MIN as i128 + 2 * p as i128 + 256;
            let hi = EXPONENT_MAX as i128 - 256;
            if (lo..=hi).contains(&n) {
                return Ok(None);
            }
        }
        // S and its relative error 2^s_err: S = 1 while the whole correction (< 2^(lr+1)) is
        // below the working precision, else the series.
        let (s, s_err) = if lr <= -(w as i128) - 16 {
            (N::from_word(1, WORD_BIT_SIZE)?, lr + 1)
        } else {
            let ws = w + WORD_BIT_SIZE;
            let half = m.div(&two, ws, RoundingMode::None)?;
            // (x/2)² = (m/2)² · 2^(2ex), scaled in two steps (2ex can pass the exponent type).
            let mut hh = half.mul(&half, ws, RoundingMode::None)?;
            hh = hh
                .ldexp(ex, ws, RoundingMode::None)?
                .ldexp(ex, ws, RoundingMode::None)?;
            if hh.is_zero() || hh.is_subnormal() {
                // ν within ~2^EXPONENT_MIN of an integer: t_1 would lose bits; refuse.
                return Err(Error::PrecisionRetryExhausted);
            }
            if !i {
                hh = hh.neg()?;
            }
            let mut t = N::from_word(1, ws)?;
            let mut s = N::from_word(1, ws)?;
            let mut k = 0usize;
            let stop = -(ws as i64) - 4;
            loop {
                k += 1;
                let kk = from_usize_p(k, ws)?;
                let den =
                    kk.add(&a, pa + ws, RoundingMode::None)?
                        .mul(&kk, ws, RoundingMode::None)?;
                t = t
                    .mul(&hh, ws, RoundingMode::None)?
                    .div(&den, ws, RoundingMode::None)?;
                s = s.add(&t, ws, RoundingMode::None)?;
                if t.is_zero() || t.exponent() < stop {
                    break;
                }
            }
            let kb = (usize::BITS - (k + 5).leading_zeros()) as i128;
            (s, kb + 1 - ws as i128)
        };
        let sign = if neg { Sign::Neg } else { Sign::Pos };
        let wp = round_p(g.mantissa_max_bit_len().max(p) + WORD_BIT_SIZE);
        if g_exact && lr < -(wp as i128) - 2 {
            // δ has the sign of t_1: − for J, + for I (an integer order ν ≥ 0 here).
            let mut d = N::from_word(1, WORD_BIT_SIZE)?;
            d.set_exponent(g.exponent() - wp as Exponent);
            if !i {
                d = d.neg()?;
            }
            let mut v = g.add(&d, wp + 2 * WORD_BIT_SIZE, RoundingMode::None)?;
            v.set_sign(sign);
            let mut r = scale_in_rm(v, n, p, rm)?;
            r.set_inexact(true);
            return Ok(Some(r));
        }
        let wp = w + 2 * WORD_BIT_SIZE;
        let pv = g.mul(&s, wp, RoundingMode::None)?;
        let c = (s_err.max(-(w as i128) - 10) + 2).max(-(wp as i128) + 4);
        let mut eps = N::from_word(1, WORD_BIT_SIZE)?;
        eps.set_exponent(c as Exponent + 1);
        let up = N::from_word(1, WORD_BIT_SIZE)?.add(&eps, wp, RoundingMode::Up)?;
        let dn = N::from_word(1, WORD_BIT_SIZE)?.sub(&eps, wp, RoundingMode::Down)?;
        let mut lo = pv.mul(&dn, wp, RoundingMode::Down)?;
        let mut hi = pv.mul(&up, wp, RoundingMode::Up)?;
        lo.set_sign(sign);
        hi.set_sign(sign);
        let rl = scale_in_rm(lo, n, p, rm);
        let rh = scale_in_rm(hi, n, p, rm);
        match (rl, rh) {
            (Err(Error::ExponentOverflow(a1)), Err(Error::ExponentOverflow(a2))) if a1 == a2 => {
                return Err(Error::ExponentOverflow(a1));
            }
            (Ok(mut u), Ok(v)) if u.cmp(&v) == 0 && u.sign() == v.sign() => {
                u.set_inexact(true);
                return Ok(Some(u));
            }
            (Err(e), _) | (_, Err(e)) if !matches!(e, Error::ExponentOverflow(_)) => return Err(e),
            _ => {}
        }
        bump_prec_retry(&mut w, &mut w_inc, p)?;
    }
}

/// ZF-37: `Y_ν(x)` (`is_k = false`) or `K_ν(x)` for tiny `x > 0`. The leading
/// factor `Γ(μ) (x/2)^{-μ}`, `μ = |ν| > 0`, is held as `g · 2^N` (`N` unbounded)
/// and rounded once through [`scale_in_rm`] / [`overflow_in_rm`]. `None` when
/// the first-correction ratio is not small (`ν = 0` stays on the existing log
/// series, which is already right).
///
/// - `K_ν = K_μ ∼ (1/2) Γ(μ) (x/2)^{-μ} S` (DLMF 10.30.2 / 10.31).
/// - `Y_μ ∼ -(1/π) Γ(μ) (x/2)^{-μ} S` (DLMF 10.7.4). `Y_{-n} = (-1)^n Y_n`.
/// - `Y_{-μ}` for non-integer `μ` carries an extra `cos(μπ)`, except at
///   half-odd-integer `μ = n + 1/2` where that term vanishes and
///   `Y_{-(n+1/2)} = (-1)^n J_{n+1/2}` (delegated to [`bessel_ji_scaled`]).
fn bessel_yk_scaled(
    x: &ExactNumNumber,
    nu: &ExactNumNumber,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
    is_k: bool,
) -> Result<Option<ExactNumNumber>, Error> {
    type N = ExactNumNumber;
    if !x.is_positive() || nu.is_zero() {
        return Ok(None);
    }
    // Y_{-(n+1/2)}(x) = (-1)^n J_{n+1/2}(x): the huge term is exactly zero.
    if !is_k && nu.is_negative() {
        if let Some(n) = half_integer_n(&nu.abs()?, p)? {
            let mu = nu.abs()?;
            let rm_j = if n % 2 == 1 { crate::common::util::invert_rm_for_sign(rm) } else { rm };
            if let Some(j) = bessel_ji_scaled(x, &mu, p, rm_j, cc, false)? {
                return Ok(Some(if n % 2 == 0 { j } else { j.neg()? }));
            }
            // ji_scaled is None when the leading J term is already in the normal
            // range (e.g. x = 2^{-1000}, p = 2048). The identity still holds:
            // evaluate J_μ itself and apply (−1)^n. Never fall through to Y.
            let j = x.bessel_j_nu(&mu, p, rm_j, cc)?;
            return Ok(Some(if n % 2 == 0 { j } else { j.neg()? }));
        }
    }

    let mut a = nu.abs()?;
    a.set_inexact(false);
    let pa = round_p(a.mantissa_max_bit_len().max(WORD_BIT_SIZE) + 2 * WORD_BIT_SIZE);
    let one = N::from_word(1, WORD_BIT_SIZE)?;
    let two = N::from_word(2, WORD_BIT_SIZE)?;
    let int = a.is_int();
    let den = if int {
        // |1 - μ|; μ ≥ 1 integer: μ - 1.
        if a.cmp(&one) <= 0 {
            // μ = 1: relative next term is O(x^2 log x), still ~ 2^(2 ex).
            N::from_word(1, WORD_BIT_SIZE)?
        } else {
            a.sub(&one, pa, RoundingMode::None)?
        }
    } else {
        one.sub(&a, pa, RoundingMode::None)?.abs()?
    };
    let (m, ex) = x.abs()?.frexp()?;
    let lr = 2 * (ex as i128 - 1) - (den.exponent() as i128 - 1);
    if lr > -3 {
        return Ok(None);
    }
    // Remainder after t_1 is O((x/2)^4) ~ 2^{2 lr}; after the μ=1 log term,
    // O(x^4 ln x) ~ 2^{4(ex-1)}. If that is still on the p-bit grid, leave
    // this path: the regular series is fast (Y_ν(2^{-30}) at p=2048 is
    // 0.007 s there; spinning the retry budget here was ~29 s).
    let rem_est = if int && a.cmp(&one) <= 0 { 4 * (ex as i128 - 1) + 16 } else { 2 * lr + 16 };
    if rem_est > -(p as i128) - 16 {
        return Ok(None);
    }

    // Sign of the leading term.
    let mut sign = if is_k { Sign::Pos } else { Sign::Neg };
    if !is_k && nu.is_negative() {
        if int {
            // Y_{-n} = (-1)^n Y_n, Y_n < 0 for small x: flip when n is odd? Y_n
            // leading is negative; (-1)^n * (negative) is positive iff n is odd.
            let odd = !a.is_zero() && !a.div(&two, pa, RoundingMode::None)?.is_int();
            if odd {
                sign = Sign::Pos;
            }
        } else {
            // Y_{-μ} ∼ -(1/π) Γ(μ) (x/2)^{-μ} cos(μπ): flip when cos(μπ) < 0.
            let pi = cc.pi_num(pa, RoundingMode::None)?;
            let c = a
                .mul(&pi, pa, RoundingMode::None)?
                .cos(pa, RoundingMode::None, cc)?;
            if c.is_negative() {
                sign = Sign::Pos;
            }
            // |cos| extremely small would have been a half-odd-integer, above.
        }
    }

    let mb = m.mantissa_max_bit_len();
    let small_int = if int && a.exponent() <= 16 {
        let k = if a.is_zero() { 0 } else { a.int_as_usize()? };
        (k >= 1 && k.saturating_mul(mb) <= 1 << 16).then_some(k)
    } else {
        None
    };

    let mut w = round_p(p + 2 * WORD_BIT_SIZE);
    let mut w_inc = WORD_BIT_SIZE;
    loop {
        // g · 2^n = Γ(μ) (x/2)^{-μ}.
        let (mut g, n, g_exact) = if let Some(k) = small_int {
            let pm = round_p(k * mb + WORD_BIT_SIZE);
            let mut mk = N::from_word(1, pm)?;
            for _ in 0..k {
                mk = mk.mul(&m, pm, RoundingMode::None)?;
            }
            let bl = (usize::BITS - (k.saturating_sub(1)).leading_zeros()) as usize;
            let fk = factorial(k - 1, round_p((k - 1) * bl.max(1) + WORD_BIT_SIZE))?;
            let pq = pm.max(w + WORD_BIT_SIZE);
            let mut g = fk.div(&mk, pq, RoundingMode::None)?;
            let exact = !g.inexact();
            g.set_inexact(false);
            let k = k as i128;
            (g, -k * (ex as i128 - 1), exact)
        } else {
            let mut q = w + 4 * WORD_BIT_SIZE;
            let lam = loop {
                let ln2 = cc.ln_2_num(q, RoundingMode::None)?;
                let lx = x.ln(q, RoundingMode::None, cc)?;
                let lg = a.ln_gamma(q, RoundingMode::None, cc)?;
                // Λ = log2(Γ(μ)) - μ (log2(x) - 1) = lnΓ(μ)/ln2 - μ log2(x) + μ
                let lam = lg
                    .div(&ln2, q, RoundingMode::None)?
                    .sub(
                        &a.mul(&lx, q, RoundingMode::None)?
                            .div(&ln2, q, RoundingMode::None)?,
                        q,
                        RoundingMode::None,
                    )?
                    .add(&a, q, RoundingMode::None)?;
                let e_max = [lg.exponent(), a.exponent(), lam.exponent(), 1]
                    .into_iter()
                    .max()
                    .unwrap_or(1)
                    .max(0) as usize;
                if e_max + 4 + w + 12 <= q {
                    break lam;
                }
                q = round_p(w + e_max + 4 * WORD_BIT_SIZE);
            };
            if lam.exponent() > 62 {
                let n = if lam.is_negative() { -(1i128 << 100) } else { 1i128 << 100 };
                (N::from_word(1, WORD_BIT_SIZE)?, n, false)
            } else {
                let fl = lam.floor_ref()?;
                let fr = lam.sub(&fl, q, RoundingMode::None)?;
                let nf = if fl.is_zero() { 0 } else { fl.abs()?.int_as_usize()? as i128 };
                let nf = if fl.is_negative() { -nf } else { nf };
                let ln2 = cc.ln_2_num(w + WORD_BIT_SIZE, RoundingMode::None)?;
                let g = fr.mul(&ln2, w + WORD_BIT_SIZE, RoundingMode::None)?.exp(
                    w + WORD_BIT_SIZE,
                    RoundingMode::None,
                    cc,
                )?;
                (g, nf, false)
            }
        };
        // K: × 1/2; Y: × 1/π. |cos(μπ)| for non-integer ν < 0 is absorbed in
        // the lo/hi interval (a factor in [0, 1] cannot un-overflow a minpos
        // leading term, and the half-odd case was already dispatched).
        if is_k {
            g = g.div(&two, w + WORD_BIT_SIZE, RoundingMode::None)?;
        } else {
            let pi = cc.pi_num(w + WORD_BIT_SIZE, RoundingMode::None)?;
            g = g.div(&pi, w + WORD_BIT_SIZE, RoundingMode::None)?;
        }
        let g_exact = g_exact && is_k && !g.inexact();

        let wp = round_p(g.mantissa_max_bit_len().max(p) + WORD_BIT_SIZE);
        if g_exact && lr < -(wp as i128) - 2 {
            // Integer K_n, n ≥ 1: the first correction is negative (K_1's log
            // term, and t_1 for n ≥ 2), so the exact leading term is an upper
            // bound. A downward tick smaller than 2^{-(p+2)} decides directed
            // rounding when the leading term is a p-bit dyadic.
            let mut d = N::from_word(1, WORD_BIT_SIZE)?;
            d.set_exponent(g.exponent() - wp as Exponent);
            let mut v = g.sub(&d, wp + 2 * WORD_BIT_SIZE, RoundingMode::None)?;
            v.set_sign(sign);
            let mut r = scale_in_rm(v, n, p, rm)?;
            r.set_inexact(true);
            return Ok(Some(r));
        }

        let wp = w + 2 * WORD_BIT_SIZE;
        // S = 1 only when every correction (t_1 and, for non-integer μ, the
        // opposite-power term (x/2)^{2μ}) sits below the working precision.
        let opp_maybe = !int && a.exponent() <= 16;
        let (s, s_err) = if lr <= -(w as i128) - 16 && !opp_maybe {
            (N::from_word(1, WORD_BIT_SIZE)?, lr + 1)
        } else if int && a.cmp(&one) <= 0 {
            // μ = 1: K_1 has S = 1 + (x^2/2)(ln(x/2) + γ − 1/2); Y_1 has S = 1 − that
            // (the log term increases |Y_1|). Remainder O(x^4 ln x).
            let ws = w + WORD_BIT_SIZE;
            let e_m1 = ex as i128 - 1;
            if e_m1.unsigned_abs() > Word::MAX as u128 {
                (N::from_word(1, WORD_BIT_SIZE)?, lr + 1)
            } else {
                let ln2 = cc.ln_2_num(ws, RoundingMode::None)?;
                let lm = m.ln(ws, RoundingMode::None, cc)?;
                let mut ev = N::from_word(e_m1.unsigned_abs() as Word, ws)?;
                if e_m1 < 0 {
                    ev = ev.neg()?;
                }
                let lxd2 =
                    ev.mul(&ln2, ws, RoundingMode::None)?
                        .add(&lm, ws, RoundingMode::None)?;
                let gam = cc.euler_gamma_num(ws, RoundingMode::None)?;
                let half = one.div(&two, ws, RoundingMode::None)?;
                let paren =
                    lxd2.add(&gam, ws, RoundingMode::None)?
                        .sub(&half, ws, RoundingMode::None)?;
                let mut xx2 = m.mul(&m, ws, RoundingMode::None)?;
                xx2 = xx2
                    .ldexp(ex, ws, RoundingMode::None)?
                    .ldexp(ex, ws, RoundingMode::None)?
                    .div(&two, ws, RoundingMode::None)?;
                let t = xx2.mul(&paren, ws, RoundingMode::None)?;
                let s = if is_k {
                    one.add(&t, ws, RoundingMode::None)?
                } else {
                    one.sub(&t, ws, RoundingMode::None)?
                };
                let rem = 4 * e_m1 + (e_m1.unsigned_abs().max(1).ilog2() as i128) + 8;
                (s, rem)
            }
        } else {
            // t_1 / t_0 = (x/2)^2 / (μ − 1); K: S = 1 − t_1, Y: S = 1 + t_1.
            // For non-integer μ, also subtract the opposite-power term
            // (x/2)^{2μ} Γ(1−μ)/Γ(1+μ) (times cos(μπ) for Y), which is O(x^{2μ})
            // and sits inside a 2048-bit grid when μ is a fraction and x = 2^{-1000}.
            let ws = w + WORD_BIT_SIZE;
            let half = m.div(&two, ws, RoundingMode::None)?;
            let mut hh = half.mul(&half, ws, RoundingMode::None)?;
            hh = hh
                .ldexp(ex, ws, RoundingMode::None)?
                .ldexp(ex, ws, RoundingMode::None)?;
            let t1 = hh.div(
                &a.sub(&one, pa + ws, RoundingMode::None)?,
                ws,
                RoundingMode::None,
            )?;
            let mut s = if is_k {
                one.sub(&t1, ws, RoundingMode::None)?
            } else {
                one.add(&t1, ws, RoundingMode::None)?
            };
            let mut rem = 2 * lr + 8;
            if opp_maybe {
                let q = ws + WORD_BIT_SIZE;
                let ln2 = cc.ln_2_num(q, RoundingMode::None)?;
                let lm = m.ln(q, RoundingMode::None, cc)?;
                let e_m1 = ex as i128 - 1;
                if e_m1.unsigned_abs() <= Word::MAX as u128 {
                    let mut ev = N::from_word(e_m1.unsigned_abs() as Word, q)?;
                    if e_m1 < 0 {
                        ev = ev.neg()?;
                    }
                    let lxd2 =
                        ev.mul(&ln2, q, RoundingMode::None)?
                            .add(&lm, q, RoundingMode::None)?;
                    let two_mu_lx =
                        a.mul(&lxd2, q, RoundingMode::None)?
                            .mul(&two, q, RoundingMode::None)?;
                    let lg_m =
                        one.sub(&a, q, RoundingMode::None)?
                            .ln_gamma(q, RoundingMode::None, cc)?;
                    let lg_p =
                        one.add(&a, q, RoundingMode::None)?
                            .ln_gamma(q, RoundingMode::None, cc)?;
                    let lam = two_mu_lx.add(&lg_m, q, RoundingMode::None)?.sub(
                        &lg_p,
                        q,
                        RoundingMode::None,
                    )?;
                    let rterm = lam.exp(q, RoundingMode::None, cc)?;
                    if is_k {
                        s = s.sub(&rterm, q, RoundingMode::None)?;
                    } else {
                        let pi = cc.pi_num(q, RoundingMode::None)?;
                        let c =
                            a.mul(&pi, q, RoundingMode::None)?
                                .cos(q, RoundingMode::None, cc)?;
                        s = s.sub(
                            &rterm.mul(&c, q, RoundingMode::None)?,
                            q,
                            RoundingMode::None,
                        )?;
                    }
                    // leftover r · (x/2)^2
                    rem = rem.max(rterm.exponent() as i128 + 2 * (ex as i128 - 1) + 8);
                }
            }
            (s, rem)
        };
        let pv = g.mul(&s, wp, RoundingMode::None)?;
        let c = (s_err.max(-(w as i128) - 10) + 2).max(-(wp as i128) + 4);
        let mut eps = N::from_word(1, WORD_BIT_SIZE)?;
        eps.set_exponent(c as Exponent + 1);
        let up = N::from_word(1, WORD_BIT_SIZE)?.add(&eps, wp, RoundingMode::Up)?;
        let dn = N::from_word(1, WORD_BIT_SIZE)?.sub(&eps, wp, RoundingMode::Down)?;
        let mut lo = pv.mul(&dn, wp, RoundingMode::Down)?;
        let mut hi = pv.mul(&up, wp, RoundingMode::Up)?;
        lo.set_sign(sign);
        hi.set_sign(sign);
        let rl = scale_in_rm(lo, n, p, rm);
        let rh = scale_in_rm(hi, n, p, rm);
        match (rl, rh) {
            (Err(Error::ExponentOverflow(a1)), Err(Error::ExponentOverflow(a2))) if a1 == a2 => {
                return overflow_in_rm(a1, p, rm).map(Some);
            }
            (Ok(mut u), Ok(v)) if u.cmp(&v) == 0 && u.sign() == v.sign() => {
                u.set_inexact(true);
                return Ok(Some(u));
            }
            (Err(e), _) | (_, Err(e)) if !matches!(e, Error::ExponentOverflow(_)) => return Err(e),
            _ => {}
        }
        match bump_prec_retry(&mut w, &mut w_inc, p) {
            Ok(()) => {}
            Err(Error::PrecisionRetryExhausted) => return Ok(None),
            Err(e) => return Err(e),
        }
    }
}

/// \(J_\nu(\pm0)\) (`i = false`) or \(I_\nu(\pm0)\). The power series' only surviving term is
/// \((x/2)^\nu/\Gamma(\nu+1)\) (DLMF 10.2.2, 10.25.2), so its sum is an exact zero or a pole and the
/// cancellation claim of [`ExactNumNumber::bessel_series`] (0 bits for a zero sum) must not see it.
/// - \(\nu=0\): \(1\).
/// - integer \(n\ne0\): a zero, \(J_n(-x)=(-1)^nJ_n(x)\), \(J_{-n}=(-1)^nJ_n\), \(I_{-n}=I_n\).
/// - non-integer \(\nu>0\): \(+0\) (\(\pm0\) alike, as C99 `pow(±0, ν)`).
/// - non-integer \(\nu<0\): a pole, signed \(\mathrm{sgn}\,\Gamma(\nu+1)=(-1)^m\), \(m=\lfloor-\nu\rfloor\),
///   returned as `ExponentOverflow` so the wrapper gives signed Inf (C99 7.12.1, as for \(Y_\nu(0)\)).
fn bessel_ji_at_zero(
    x: &ExactNumNumber,
    nu: &ExactNumNumber,
    p: usize,
    i: bool,
) -> Result<ExactNumNumber, Error> {
    let inexact = x.inexact() | nu.inexact();
    if nu.is_zero() {
        let mut one = ExactNumNumber::from_word(1, p)?;
        one.set_inexact(inexact);
        return Ok(one);
    }
    let a = nu.abs()?;
    let pe = round_p(a.mantissa_max_bit_len().max(WORD_BIT_SIZE) + WORD_BIT_SIZE);
    let m = a.floor()?;
    let two = ExactNumNumber::from_word(2, WORD_BIT_SIZE)?;
    let m_odd = !m.is_zero() && !m.div(&two, pe, RoundingMode::None)?.is_int();
    if nu.is_int() {
        let neg = (x.is_negative() && m_odd) != (!i && nu.is_negative() && m_odd);
        return ExactNumNumber::new2(p, if neg { Sign::Neg } else { Sign::Pos }, inexact);
    }
    if nu.is_positive() {
        return ExactNumNumber::new2(p, Sign::Pos, inexact);
    }
    Err(Error::ExponentOverflow(if m_odd {
        Sign::Neg
    } else {
        Sign::Pos
    }))
}

/// \(Y_\nu(0^+)\): `-Inf` for `nu >= 0`; for `nu = -m` it is `(-1)^m Y_m(0^+) = (-1)^{m+1} Inf`;
/// for `nu = -(m + 1/2)` it is `(-1)^m J_{m+1/2}(0) = (-1)^m 0`; otherwise the
/// `J_nu cot(nu pi)` term dominates, with sign `sgn Gamma(nu+1) sgn cot(nu pi)`
/// `= (-1)^m (+1 if t > 1/2 else -1)`, where `m = floor(-nu)` and `t = -nu - m`.
/// Infinities are returned as `ExponentOverflow` so the wrapper maps them to signed Inf.
fn bessel_y_at_zero(nu: &ExactNumNumber, p: usize) -> Result<ExactNumNumber, Error> {
    let inf = |pos: bool| {
        Err(Error::ExponentOverflow(if pos {
            Sign::Pos
        } else {
            Sign::Neg
        }))
    };
    if nu.is_zero() || nu.is_positive() {
        return inf(false);
    }
    let a = nu.abs()?;
    let pe = round_p(a.mantissa_max_bit_len().max(WORD_BIT_SIZE) + WORD_BIT_SIZE);
    let m = nu.abs()?.floor()?;
    let two = ExactNumNumber::from_word(2, WORD_BIT_SIZE)?;
    let m_odd = !m.is_zero() && !m.div(&two, pe, RoundingMode::None)?.is_int();
    let t = a.sub(&m, pe, RoundingMode::None)?;
    if t.is_zero() {
        return inf(m_odd);
    }
    let half = ExactNumNumber::from_word(1, WORD_BIT_SIZE)?.ldexp(
        -1,
        WORD_BIT_SIZE,
        RoundingMode::None,
    )?;
    let c = t.cmp(&half);
    if c == 0 {
        return ExactNumNumber::new2(p, if m_odd { Sign::Neg } else { Sign::Pos }, false);
    }
    inf((c > 0) != m_odd)
}

fn as_i32_exact(x: &ExactNumNumber, p: usize) -> Result<Option<i32>, Error> {
    let _ = p;
    if x.is_zero() {
        return Ok(Some(0));
    }
    if !x.is_int() {
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
    // Product in ExactInt (one-word factors), then a single rounding to `p`.
    // Multiplying by `from_word(i, p)` used to inflate each factor to `p` bits
    // and made n! at 2*10^5 bits take tens of thousands of full-width muls.
    if n <= 1 {
        return ExactNumNumber::from_word(1, p);
    }
    let mut acc = ExactInt::one();
    for i in 2..=n {
        acc.mul_u64_assign(i as u64);
    }
    let limbs = acc.as_limbs();
    let e = Exponent::try_from(limbs.len().saturating_mul(WORD_BIT_SIZE))
        .map_err(|_| Error::InvalidArgument)?;
    let sign = if acc.is_negative() { Sign::Neg } else { Sign::Pos };
    let mut f = ExactNumNumber::from_words(limbs, sign, e)?;
    let _ = f.set_precision(p, RoundingMode::ToEven);
    Ok(f)
}

/// Even Bernoulli number B_{2k} via Akiyama–Tanigawa.
/// Even Bernoulli numbers \(B_2, B_4, \ldots, B_{2n}\), each rounded once to `p` bits.
///
/// Tangent numbers \(T_k\) are built exactly in integer arithmetic (Brent & Harvey,
/// "Fast computation of Bernoulli, Tangent and Secant numbers", 2011, algorithm
/// *TangentNumbers*), then \(B_{2k} = (-1)^{k-1}\,2k\,T_k / (2^{2k}(2^{2k}-1))\) with one
/// division. (The previous Akiyama–Tanigawa recurrence in floating point cancelled badly: at
/// 192–256 bits Γ lost 60–90 bits, more with 32-bit limbs.)
pub(crate) fn even_bernoulli(n: usize, p: usize) -> Result<Vec<ExactNumNumber>, Error> {
    if n == 0 {
        return Ok(Vec::new());
    }
    // T_k < (2k)! < 2^{2k·log2(2k)}: this precision holds every integer exactly.
    let log2 = usize::BITS as usize - (2 * n).leading_zeros() as usize;
    let pi = round_p(2 * n * (log2 + 1) + 2 * WORD_BIT_SIZE);
    let exact = RoundingMode::None;
    let word = |v: usize| ExactNumNumber::from_word(v as Word, pi);
    let mut t: Vec<ExactNumNumber> = Vec::with_capacity(n + 1);
    t.push(ExactNumNumber::new(pi)?);
    t.push(word(1)?);
    for k in 2..=n {
        let v = t[k - 1].mul(&word(k - 1)?, pi, exact)?;
        t.push(v);
    }
    for k in 2..=n {
        for j in k..=n {
            let a = t[j - 1].mul(&word(j - k)?, pi, exact)?;
            let b = t[j].mul(&word(j - k + 2)?, pi, exact)?;
            t[j] = a.add(&b, pi, exact)?;
        }
    }
    let one = word(1)?;
    let mut out = Vec::with_capacity(n);
    for (k, tk) in t.iter().enumerate().skip(1) {
        let num = tk.mul(&word(2 * k)?, pi, exact)?;
        let mut four_k = word(1)?;
        four_k.set_exponent(four_k.exponent() + 2 * k as Exponent);
        let den = four_k.sub(&one, pi, exact)?;
        let mut b = num.div(&den, p, RoundingMode::ToEven)?;
        b.set_exponent(b.exponent() - 2 * k as Exponent);
        if k % 2 == 0 {
            b = b.neg()?;
        }
        out.push(b);
    }
    Ok(out)
}

/// FLINT `bernoulli_cache_compute` growth, in even-`k` units (`k` so `B_{2k}`).
/// FLINT sizes by the index `n`: `+32` when `n≤128`, else `+128`. That is
/// `+16` / `+64` in the even-`k` count.
fn flint_even_k_len(have: usize, need: usize) -> usize {
    if need <= 64 {
        have.saturating_add(16).max(need)
    } else {
        have.saturating_add(64).max(need)
    }
}

#[cfg(test)]
fn bitlen_usize(n: usize) -> usize {
    (usize::BITS - n.leading_zeros()) as usize
}

/// Proven upper bound on `log2(|B_n|)` for even `n ≥ 2`: `|B_n| < 4 n!/(2π)^n`
/// because `ζ(n) < 2`. `log2(n!)` is bounded by `Σ_{i=2}^n ⌈log2 i⌉`, and
/// `log2(2π) > 53/20`.
#[cfg(test)]
fn bernoulli_abs_log2_ub(n: usize) -> usize {
    let mut log2_fact = 0usize;
    for i in 2..=n {
        log2_fact = log2_fact.saturating_add(bitlen_usize(i));
    }
    let sub = n.saturating_mul(53) / 20;
    2usize.saturating_add(log2_fact).saturating_sub(sub).max(8)
}

fn pow_ui(base: &ExactNumNumber, exp: usize, p: usize) -> Result<ExactNumNumber, Error> {
    let mut acc = ExactNumNumber::from_word(1, p)?;
    if exp == 0 {
        return Ok(acc);
    }
    // base^(2^j); `None` while j = 0 (base read in place).
    let mut b_slot: Option<ExactNumNumber> = None;
    let mut e = exp;
    loop {
        let b = b_slot.as_ref().unwrap_or(base);
        if e & 1 == 1 {
            acc = acc.mul(b, p, RoundingMode::None)?;
        }
        e >>= 1;
        if e == 0 {
            return Ok(acc);
        }
        let sq = b.mul(b, p, RoundingMode::None)?;
        b_slot = Some(sq);
    }
}

fn primes_upto(max: usize) -> Vec<u32> {
    if max < 2 {
        return Vec::new();
    }
    let mut is_p = Vec::new();
    is_p.resize(max + 1, true);
    is_p[0] = false;
    is_p[1] = false;
    let mut p = 2usize;
    while p.saturating_mul(p) <= max {
        if is_p[p] {
            let mut m = p.saturating_mul(p);
            while m <= max {
                is_p[m] = false;
                m = m.saturating_add(p);
            }
        }
        p += 1;
    }
    (2..=max).filter(|&i| is_p[i]).map(|i| i as u32).collect()
}

/// Primes `P` so the Euler-product tail `P^{1-s}/(s-1) < 2^{-(wp+8)}`.
#[cfg(test)]
fn euler_product_limit(s: usize, wp: usize) -> Result<usize, Error> {
    let num = wp.saturating_add(8).saturating_add(bitlen_usize(s));
    let den = s.saturating_sub(1).max(1);
    let e = num.div_ceil(den);
    if e >= 24 {
        return Err(Error::PrecisionRetryExhausted);
    }
    Ok((1usize << e).saturating_add(2))
}

#[cfg(test)]
fn zeta_inv_euler_product(s: usize, wp: usize) -> Result<ExactNumNumber, Error> {
    let pmax = euler_product_limit(s, wp)?;
    let primes = primes_upto(pmax);
    let one = ExactNumNumber::from_word(1, wp)?;
    let mut inv = ExactNumNumber::from_word(1, wp)?;
    for q in primes {
        let qn = pow_ui(&from_usize_p(q as usize, wp)?, s, wp)?;
        let term = one.sub(
            &one.div(&qn, wp, RoundingMode::None)?,
            wp,
            RoundingMode::None,
        )?;
        inv = inv.mul(&term, wp, RoundingMode::None)?;
    }
    Ok(inv)
}

#[cfg(test)]
fn zeta_ui_asymp(s: usize, wp: usize) -> Result<ExactNumNumber, Error> {
    let one = ExactNumNumber::from_word(1, wp)?;
    let mut z = ExactNumNumber::from_word(1, wp)?;
    let lim = Exponent::try_from(wp.saturating_add(8)).unwrap_or(Exponent::MAX);
    let mut j = 2usize;
    loop {
        let jn = pow_ui(&from_usize_p(j, wp)?, s, wp)?;
        let term = one.div(&jn, wp, RoundingMode::None)?;
        if term.is_zero() || term.exponent() < -lim {
            break;
        }
        z = z.add(&term, wp, RoundingMode::None)?;
        j += 1;
        if j > (1 << 20) {
            return Err(Error::PrecisionRetryExhausted);
        }
    }
    Ok(z)
}

/// Working precision for `B_{2k}` at dest `p`. `bernoulli_abs_log2_ub` is a
/// bound on the magnitude, not a significant-bit count. When |B_{2k}| has more
/// than `p` bits, work at p+64. When it is a short exact rational, work at
/// enough bits to recover it (von Staudt–Clausen: `log2(d_n) < n`).
#[cfg(test)]
fn even_bernoulli_zeta_wp(k: usize, p: usize) -> Result<usize, Error> {
    let n = k.checked_mul(2).ok_or(Error::InvalidArgument)?;
    let bits = bernoulli_abs_log2_ub(n);
    let ident = bits.saturating_add(n).saturating_add(256);
    Ok(round_p(if ident > p {
        p.saturating_add(64)
    } else {
        ident
    }))
}

/// `F = 2 (2k)! ζ(2k) / (2π)^{2k}` at `wp`, then the sign of `B_{2k}`.
#[cfg(test)]
fn even_bernoulli_from_f(
    f: ExactNumNumber,
    k: usize,
    wp: usize,
    p: usize,
) -> Result<ExactNumNumber, Error> {
    let n = k.checked_mul(2).ok_or(Error::InvalidArgument)?;
    let mut b = if (n as u128).saturating_mul(10) > (wp as u128).saturating_mul(7) {
        let z = zeta_ui_asymp(n, wp)?;
        f.mul(&z, wp, RoundingMode::None)?
    } else {
        let inv = zeta_inv_euler_product(n, wp)?;
        f.div(&inv, wp, RoundingMode::None)?
    };
    if n.is_multiple_of(4) {
        b = b.neg()?;
    }
    if wp != p {
        b.set_precision(p, RoundingMode::ToEven)?;
    }
    Ok(b)
}

/// `B_{2k}` at `p` bits ToEven, from
/// `B_{2k} = (-1)^{k+1} 2 (2k)! ζ(2k) / (2π)^{2k}`
/// (Arb `arb_bernoulli_ui_zeta`: Euler product, or the short tail when `2k > 0.7 wp`).
#[cfg(test)]
fn even_bernoulli_via_zeta(k: usize, p: usize) -> Result<ExactNumNumber, Error> {
    if k == 0 {
        return Err(Error::InvalidArgument);
    }
    let n = k.checked_mul(2).ok_or(Error::InvalidArgument)?;
    let wp = even_bernoulli_zeta_wp(k, p)?;
    let piwp = round_p(wp.saturating_add(2 * bitlen_usize(n)));
    let mut cc = Consts::new()?;
    let fac = factorial(n, piwp)?;
    let pi = cc.pi_num(piwp, RoundingMode::None)?;
    let two_pi = ExactNumNumber::from_word(2, piwp)?.mul(&pi, piwp, RoundingMode::None)?;
    let pow = pow_ui(&two_pi, n, piwp)?;
    let two = ExactNumNumber::from_word(2, wp)?;
    let f = two
        .mul(&fac, wp, RoundingMode::None)?
        .div(&pow, wp, RoundingMode::None)?;
    even_bernoulli_from_f(f, k, wp, p)
}

/// Performance banding only: split when the zeta-tail / Euler-product length
/// would exceed this, so each band does one factorial and one term table.
/// Not a correctness cap (P9): every `B_{2k}` is still computed at
/// [`even_bernoulli_zeta_wp`]. 8192 covers the mid-`k` gap at p=65536
/// (`e=12`, `J=4098`) in one extra band instead of hundreds of singleton
/// factorials, which is what kept `digamma_half_at_high_p_from_error_bound`
/// from finishing.
#[cfg(test)]
const BERNOULLI_ZETA_BAND_J_MAX: usize = 8192;

/// Fill `B_{2k}` for `k0..=k1` by walking `F(n)=2 n!/(2π)^n` down from
/// `n=2 k1` (`F(n-2) = F(n) (2π)^2 / (n(n-1))`), one factorial at the top
/// of each band. When `n ≤ 0.7 wp` the Euler product is walked (primes
/// only, as in Arb `arb_bernoulli_ui_zeta`); otherwise the short zeta
/// tail over every `j`. A `PrecisionRetryExhausted` from
/// [`euler_product_limit`] ends the band rather than failing the fill.
#[cfg(test)]
fn even_bernoulli_zeta_range(k0: usize, k1: usize, p: usize) -> Result<Vec<ExactNumNumber>, Error> {
    if k0 == 0 || k0 > k1 {
        return Err(Error::InvalidArgument);
    }
    let mut out: Vec<Option<ExactNumNumber>> = (0..=(k1 - k0)).map(|_| None).collect();
    let mut cc = Consts::new()?;
    let mut k_hi = k1;
    while k_hi >= k0 {
        // `wp(k)` is non-decreasing in `k`, so the band works at `wp(k_hi)`.
        let wp = even_bernoulli_zeta_wp(k_hi, p)?;
        let mut k_lo = k_hi;
        let hi_euler =
            (k_hi as u128).saturating_mul(2).saturating_mul(10) <= (wp as u128).saturating_mul(7);
        while k_lo > k0 {
            let nxt = k_lo - 1;
            let n_nxt = 2 * nxt;
            // Keep one ζ method per band: the short tail is only valid when
            // n > 0.7 wp, so crossing that threshold starts a new band.
            let nxt_euler = (n_nxt as u128).saturating_mul(10) <= (wp as u128).saturating_mul(7);
            if nxt_euler != hi_euler {
                break;
            }
            match euler_product_limit(n_nxt, wp) {
                Ok(j) if j <= BERNOULLI_ZETA_BAND_J_MAX => k_lo = nxt,
                _ => break,
            }
        }
        let n_hi = k_hi.checked_mul(2).ok_or(Error::InvalidArgument)?;
        let piwp = round_p(wp.saturating_add(2 * bitlen_usize(n_hi)));
        let fac = factorial(n_hi, piwp)?;
        let pi = cc.pi_num(piwp, RoundingMode::None)?;
        let two_pi = ExactNumNumber::from_word(2, piwp)?.mul(&pi, piwp, RoundingMode::None)?;
        let two_pi_sq = two_pi.mul(&two_pi, piwp, RoundingMode::None)?;
        let pow = pow_ui(&two_pi, n_hi, piwp)?;
        let two = ExactNumNumber::from_word(2, wp)?;
        let one = ExactNumNumber::from_word(1, wp)?;
        let mut f = two
            .mul(&fac, wp, RoundingMode::None)?
            .div(&pow, wp, RoundingMode::None)?;
        // Bloemen reverse iterator (FLINT/Arb): each stored power walks
        // down by a multiply by `q^2`, so ζ(n) / 1/ζ(n) is not rebuilt.
        //
        // Hold only the terms needed at the current n. Pre-building
        // j_max(n_lo) powers at n_hi underflows for the larger primes
        // (`q^{-n_hi}` is 0 at `wp`), and a later multiply by `q^2`
        // cannot recover them, so mid-band B_{2k} would miss Euler
        // factors. It is also the cost leak: a p=65536 band from
        // k=4097 to 10596 used 54 primes on every step, and the
        // k=2732..4096 band used ~564, even though the top of each
        // band needs only a handful. New primes (or j) are added,
        // with a fresh q^{-n}, when the walk reaches an n that needs
        // them.
        let mut j_max = euler_product_limit(n_hi, wp).unwrap_or(2);
        // Same threshold as `even_bernoulli_from_f` / Arb: Euler product
        // when `n ≤ 0.7 wp`, otherwise the short tail over every `j`.
        // Chosen at n_hi; the k_lo split above keeps the method fixed
        // for the band.
        let use_euler = (n_hi as u128).saturating_mul(10) <= (wp as u128).saturating_mul(7);
        let mut primes = if use_euler { primes_upto(j_max) } else { Vec::new() };
        let mut terms: Vec<ExactNumNumber> = Vec::new();
        if use_euler {
            for &q in &primes {
                let qn = pow_ui(&from_usize_p(q as usize, wp)?, n_hi, wp)?;
                terms.push(one.div(&qn, wp, RoundingMode::None)?);
            }
        } else {
            for j in 2..=j_max {
                let jn = pow_ui(&from_usize_p(j, wp)?, n_hi, wp)?;
                terms.push(one.div(&jn, wp, RoundingMode::None)?);
            }
        }
        let tiny = Exponent::try_from(wp.saturating_add(8)).unwrap_or(Exponent::MAX);
        let mut k = k_hi;
        loop {
            let n_cur = 2 * k;
            if let Ok(need) = euler_product_limit(n_cur, wp) {
                if need > j_max {
                    if use_euler {
                        let have = primes.last().copied().unwrap_or(0) as usize;
                        for q in primes_upto(need) {
                            if (q as usize) <= have {
                                continue;
                            }
                            let qn = pow_ui(&from_usize_p(q as usize, wp)?, n_cur, wp)?;
                            terms.push(one.div(&qn, wp, RoundingMode::None)?);
                            primes.push(q);
                        }
                    } else {
                        for j in (j_max + 1)..=need {
                            let jn = pow_ui(&from_usize_p(j, wp)?, n_cur, wp)?;
                            terms.push(one.div(&jn, wp, RoundingMode::None)?);
                        }
                    }
                    j_max = need;
                }
            }
            let mut b = if use_euler {
                let mut inv = ExactNumNumber::from_word(1, wp)?;
                for t in &terms {
                    if t.is_zero() || t.exponent() < -tiny {
                        continue;
                    }
                    let factor = one.sub(t, wp, RoundingMode::None)?;
                    inv = inv.mul(&factor, wp, RoundingMode::None)?;
                }
                f.div(&inv, wp, RoundingMode::None)?
            } else {
                let mut z = ExactNumNumber::from_word(1, wp)?;
                for t in &terms {
                    if t.is_zero() || t.exponent() < -tiny {
                        continue;
                    }
                    z = z.add(t, wp, RoundingMode::None)?;
                }
                f.mul(&z, wp, RoundingMode::None)?
            };
            if (2 * k).is_multiple_of(4) {
                b = b.neg()?;
            }
            if wp != p {
                b.set_precision(p, RoundingMode::ToEven)?;
            }
            out[k - k0] = Some(b);
            if k == k_lo {
                break;
            }
            let n = 2 * k;
            let den =
                from_usize_p(n, wp)?.mul(&from_usize_p(n - 1, wp)?, wp, RoundingMode::None)?;
            f = f
                .mul(&two_pi_sq, wp, RoundingMode::None)?
                .div(&den, wp, RoundingMode::None)?;
            if use_euler {
                for (t, &q) in terms.iter_mut().zip(primes.iter()) {
                    let q2 = from_usize_p((q as usize).saturating_mul(q as usize), wp)?;
                    *t = t.mul(&q2, wp, RoundingMode::None)?;
                }
            } else {
                for (i, t) in terms.iter_mut().enumerate() {
                    let j2 = from_usize_p((i + 2) * (i + 2), wp)?;
                    *t = t.mul(&j2, wp, RoundingMode::None)?;
                }
            }
            k -= 1;
        }
        if k_lo == k0 {
            break;
        }
        k_hi = k_lo - 1;
    }
    out.into_iter()
        .map(|b| b.ok_or(Error::InvalidArgument))
        .collect()
}

/// Exact even Bernoulli numbers, grown on demand, with p-bit roundings keyed by precision.
///
/// Tangent numbers `T_k` are the odd zigzag numbers `A_{2k-1}` of the Seidel–Entringer
/// boustrophedon: `E(n,0) = 0`, `E(n,j) = E(n,j-1) + E(n-1,n-j)`, `A_n = E(n,n)`. Each new row
/// costs `n` exact additions, so `B_2 .. B_2K` cost `O(K²)` additions in total (the former
/// tangent recurrence in `even_bernoulli` restarts its `O(K²)` triangle with multiplications
/// for every call). Then `B_2k = (-1)^(k-1) 2k T_k / (2^2k (2^2k - 1))`.
#[path = "bernoulli_rev.rs"]
mod bernoulli_rev;

struct BernoulliCache {
    /// Boustrophedon row `row_n` (empty before the first growth); exact integers.
    row: Vec<ExactNumNumber>,
    row_n: usize,
    /// `exact[k-1] = (num, den)` with `B_2k = num / den`, both exact.
    exact: Vec<(ExactInt, ExactInt)>,
    /// `(p, [B_2 .. B_2n] rounded to p bits ToEven)`, at most `BERNOULLI_P_SLOTS` entries.
    rounded: Vec<(usize, Vec<ExactNumNumber>)>,
}

/// Distinct precisions kept in the Bernoulli cache; the oldest is dropped beyond this.
const BERNOULLI_P_SLOTS: usize = 8;

impl BernoulliCache {
    const fn new() -> Self {
        Self {
            row: Vec::new(),
            row_n: 0,
            exact: Vec::new(),
            rounded: Vec::new(),
        }
    }

    /// Advances the boustrophedon one row.
    fn next_row(&mut self) -> Result<(), Error> {
        if self.row.is_empty() {
            self.row.push(ExactNumNumber::from_word(1, WORD_BIT_SIZE)?);
            self.row_n = 0;
        }
        let n = self.row_n + 1;
        // A_n ≤ n! < 2^(n·bitlen(n)): every entry of row n fits exactly.
        let bits =
            usize::try_from(usize::BITS - n.leading_zeros()).map_err(|_| Error::InvalidArgument)?;
        let pn = round_p(n * bits + 2 * WORD_BIT_SIZE);
        let mut next = Vec::with_capacity(n + 1);
        let mut acc = ExactNumNumber::new(pn)?;
        for j in 1..=n {
            let prev = self.row.get(n - j).ok_or(Error::InvalidArgument)?;
            let sum = acc.add(prev, pn, RoundingMode::None)?;
            next.push(acc);
            acc = sum;
        }
        next.push(acc);
        self.row = next;
        self.row_n = n;
        Ok(())
    }

    /// Grows the exact table to at least `n` entries, with FLINT cache sizing.
    /// A request that only extends a filled table walks the new high end
    /// (`exact_even_bernoulli_rev_range`) instead of recomputing the prefix.
    fn grow_exact(&mut self, n: usize) -> Result<(), Error> {
        if self.exact.len() >= n {
            return Ok(());
        }
        let n = flint_even_k_len(self.exact.len(), n);
        if n <= 15 {
            return self.grow_exact_boustrophedon(n);
        }
        if self.exact.len() > 15 {
            let k_lo = self.exact.len() + 1;
            let extra = bernoulli_rev::exact_even_bernoulli_rev_range(k_lo, n)?;
            self.exact.extend(extra);
            return Ok(());
        }
        self.exact = bernoulli_rev::exact_even_bernoulli_rev(n)?;
        Ok(())
    }

    /// Seidel–Entringer triangle for small `n` (and tests).
    fn grow_exact_boustrophedon(&mut self, n: usize) -> Result<(), Error> {
        while self.exact.len() < n {
            self.next_row()?;
            if self.row_n.is_multiple_of(2) {
                continue;
            }
            let k = self.exact.len() + 1;
            let t = self.row.last().ok_or(Error::InvalidArgument)?;
            let pb = round_p(t.mantissa_max_bit_len() + 2 * WORD_BIT_SIZE);
            let mut num = t.mul(&from_usize_p(2 * k, pb)?, pb, RoundingMode::None)?;
            if k.is_multiple_of(2) {
                num = num.neg()?;
            }
            let pd = round_p(4 * k + 2 * WORD_BIT_SIZE);
            let one = ExactNumNumber::from_word(1, pd)?;
            let mut pow = ExactNumNumber::from_word(1, pd)?;
            pow.set_exponent(
                pow.exponent() + Exponent::try_from(2 * k).map_err(|_| Error::InvalidArgument)?,
            );
            let den = pow
                .sub(&one, pd, RoundingMode::None)?
                .mul(&pow, pd, RoundingMode::None)?;
            self.exact.push((
                bernoulli_rev::num_to_int(&num)?,
                bernoulli_rev::num_to_int(&den)?,
            ));
        }
        Ok(())
    }

    /// Grows the p-bit roundings to at least `n` entries from the exact
    /// `fmpq`-style cache (FLINT `arb_set_round_fmpz` / `arb_div_fmpz`).
    fn grow_rounded(&mut self, n: usize, p: usize) -> Result<(), Error> {
        self.grow_exact(n)?;
        if !self.rounded.iter().any(|(q, _)| *q == p) {
            if self.rounded.len() >= BERNOULLI_P_SLOTS {
                self.rounded.remove(0);
            }
            self.rounded.push((p, Vec::new()));
        }
        let exact = &self.exact;
        let slot = self
            .rounded
            .iter_mut()
            .find(|(q, _)| *q == p)
            .ok_or(Error::InvalidArgument)?;
        while slot.1.len() < n {
            let (num, den) = exact.get(slot.1.len()).ok_or(Error::InvalidArgument)?;
            // Convert the exact integers at full bit length, then divide at
            // dest `p` (FLINT `arb_set_round_fmpz` / `arb_div_fmpz`).
            let pn = round_p(num.bit_length().max(64).saturating_add(WORD_BIT_SIZE));
            let pd = round_p(den.bit_length().max(64).saturating_add(WORD_BIT_SIZE));
            let nf = num.to_exact_num_number(pn, RoundingMode::None)?;
            let df = den.to_exact_num_number(pd, RoundingMode::None)?;
            slot.1.push(nf.div(&df, p, RoundingMode::ToEven)?);
        }
        Ok(())
    }

    fn rounded_slice(&self, n: usize, p: usize) -> Option<&[ExactNumNumber]> {
        self.rounded
            .iter()
            .find(|(q, _)| *q == p)
            .and_then(|(_, v)| v.get(..n))
    }
}

#[cfg(feature = "std")]
static BERNOULLI: std::sync::RwLock<BernoulliCache> = std::sync::RwLock::new(BernoulliCache::new());

/// Exact even Bernoulli numbers `B_2 .. B_{2n}` by the Seidel–Entringer
/// boustrophedon (used for the tiny table `n ≤ 15` and as a test oracle).
pub(super) fn exact_even_bernoulli_boustrophedon(
    n: usize,
) -> Result<Vec<(ExactInt, ExactInt)>, Error> {
    let mut c = BernoulliCache::new();
    c.grow_exact_boustrophedon(n)?;
    Ok(c.exact)
}

/// Exact even Bernoulli numbers `B_2 .. B_2n` as `(num, den)` `ExactInt` pairs.
#[cfg(test)]
pub(crate) fn with_exact_bernoulli<R>(
    n: usize,
    f: impl FnOnce(&[(ExactInt, ExactInt)]) -> Result<R, Error>,
) -> Result<R, Error> {
    #[cfg(feature = "std")]
    {
        {
            let g = BERNOULLI.read().map_err(|_| Error::InvalidArgument)?;
            if let Some(v) = g.exact.get(..n) {
                return f(v);
            }
        }
        BERNOULLI
            .write()
            .map_err(|_| Error::InvalidArgument)?
            .grow_exact(n)?;
        let g = BERNOULLI.read().map_err(|_| Error::InvalidArgument)?;
        f(g.exact.get(..n).ok_or(Error::InvalidArgument)?)
    }
    #[cfg(not(feature = "std"))]
    {
        let mut c = BernoulliCache::new();
        c.grow_exact(n)?;
        f(c.exact.get(..n).ok_or(Error::InvalidArgument)?)
    }
}

/// `B_n` at `p` bits (ToEven). `B_0 = 1`, `B_1 = -1/2`, odd `n > 1` is 0,
/// and even `n = 2k` uses [`with_even_bernoulli`].
#[cfg(test)]
pub(crate) fn bernoulli_ui(n: usize, p: usize) -> Result<ExactNumNumber, Error> {
    match n {
        0 => ExactNumNumber::from_word(1, p),
        1 => {
            let one = ExactNumNumber::from_word(1, p)?;
            let two = ExactNumNumber::from_word(2, p)?;
            one.div(&two, p, RoundingMode::ToEven)?.neg()
        }
        n if !n.is_multiple_of(2) => ExactNumNumber::from_word(0, p),
        n => {
            let k = n / 2;
            with_even_bernoulli(k, p, |v| {
                v.get(k - 1).ok_or(Error::InvalidArgument)?.clone()
            })
        }
    }
}

/// `B_2 .. B_2n` rounded to `p` bits (ToEven), borrowed from the shared cache for `f`.
/// Round Bernoulli work precision up to a 256-bit band so Ziv retries
/// at `p`, `p+64`, `p+128`, … share one rounded table.
fn bernoulli_prec_bucket(p: usize) -> usize {
    let p = crate::common::util::round_p(p.max(64));
    const STEP: usize = 256;
    let up = p.saturating_add(STEP - 1) / STEP * STEP;
    crate::common::util::round_p(up.max(STEP))
}

pub(crate) fn with_even_bernoulli<R>(
    n: usize,
    p: usize,
    f: impl FnOnce(&[ExactNumNumber]) -> Result<R, Error>,
) -> Result<R, Error> {
    let p = bernoulli_prec_bucket(p);
    #[cfg(feature = "std")]
    {
        {
            let g = BERNOULLI.read().map_err(|_| Error::InvalidArgument)?;
            if let Some(v) = g.rounded_slice(n, p) {
                return f(v);
            }
        }
        // Use the grown slot under the same write guard: between releasing it and taking a
        // read guard another thread may evict this precision (more than
        // `BERNOULLI_P_SLOTS` in use), which returned `InvalidArgument` for a valid request.
        // Only a miss runs `f` under the write lock; hits share the read lock as before.
        let mut g = BERNOULLI.write().map_err(|_| Error::InvalidArgument)?;
        g.grow_rounded(n, p)?;
        f(g.rounded_slice(n, p).ok_or(Error::InvalidArgument)?)
    }
    #[cfg(not(feature = "std"))]
    {
        let mut c = BernoulliCache::new();
        c.grow_rounded(n, p)?;
        f(c.rounded_slice(n, p).ok_or(Error::InvalidArgument)?)
    }
}

/// Binary-splitting state of the Brent–McMillan sums over the terms `k = a+1 ..= b`, all exact
/// integers (`P = n^(2(b−a))` is a power of two and never stored):
/// `D = Π k`, `Q = D²`, `C/D = Σ 1/k`, `T = Q·Σ_k P_{a,k}/Q_{a,k}` and
/// `V = Q·D·Σ_k (P_{a,k}/Q_{a,k})·(H_k − H_a)`, with `P_{a,k}/Q_{a,k} = Π_{i=a+1}^{k} n²/i²`.
struct BmSplit {
    d: ExactNumNumber,
    c: ExactNumNumber,
    t: ExactNumNumber,
    v: ExactNumNumber,
}

/// `x · 2^s` (exact; `x` an integer, `s ≥ 0`).
fn bm_shift(mut x: ExactNumNumber, s: Exponent) -> Result<ExactNumNumber, Error> {
    if !x.is_zero() {
        let e = x.exponent().checked_add(s).ok_or(Error::InvalidArgument)?;
        x.set_exponent(e);
    }
    Ok(x)
}

/// `2^s` as an exact integer.
fn bm_pow2(s: Exponent) -> Result<ExactNumNumber, Error> {
    bm_shift(ExactNumNumber::from_word(1, WORD_BIT_SIZE)?, s)
}

/// [`BmSplit`] over `a+1 ..= b` for `n = 2^j` (`two_j = 2j`, so `n² = 2^two_j`).
fn bm_split(a: u64, b: u64, two_j: Exponent) -> Result<BmSplit, Error> {
    if b == a + 1 {
        return Ok(BmSplit {
            d: ExactNumNumber::from_u64(b, 64)?,
            c: ExactNumNumber::from_word(1, WORD_BIT_SIZE)?,
            t: bm_pow2(two_j)?,
            v: bm_pow2(two_j)?,
        });
    }
    let m = a + (b - a) / 2;
    let l = bm_split(a, m, two_j)?;
    let r = bm_split(m, b, two_j)?;
    // P_L = n^(2(m−a)) = 2^(two_j·(m−a)).
    let len = Exponent::try_from(m - a).map_err(|_| Error::InvalidArgument)?;
    let pl = two_j.checked_mul(len).ok_or(Error::InvalidArgument)?;
    let q_r = r.d.mul_full_prec(&r.d)?;
    let t_r_d_r = r.t.mul_full_prec(&r.d)?;
    // V = V_L·Q_R·D_R + P_L·(V_R·D_L + C_L·T_R·D_R)
    let inner =
        r.v.mul_full_prec(&l.d)?
            .add_full_prec(&l.c.mul_full_prec(&t_r_d_r)?)?;
    let v =
        l.v.mul_full_prec(&q_r)?
            .mul_full_prec(&r.d)?
            .add_full_prec(&bm_shift(inner, pl)?)?;
    // T = T_L·Q_R + P_L·T_R
    let t =
        l.t.mul_full_prec(&q_r)?
            .add_full_prec(&bm_shift(r.t, pl)?)?;
    // C = C_L·D_R + C_R·D_L, D = D_L·D_R
    let c =
        l.c.mul_full_prec(&r.d)?
            .add_full_prec(&r.c.mul_full_prec(&l.d)?)?;
    let d = l.d.mul_full_prec(&r.d)?;
    Ok(BmSplit { d, c, t, v })
}

/// Euler–Mascheroni constant γ to `p` bits (Brent–McMillan B1 by exact binary splitting).
///
/// With `n = 2^j`, `t_k = (n^k/k!)²`, `A = Σ_{k≥0} t_k H_k` and `B = Σ_{k≥0} t_k`,
/// `0 < A/B − ln n − γ < π e^{−4n}` (Brent and McMillan 1980). `j` is the least with
/// `4n log2 e ≥ pw + 16`, so that error is below `2^−(pw+14)`. The sums stop at `N = 5n`
/// (> α n, α ≈ 4.9706 solving α(ln α − 1) = 3): for `k > N`, `t_{k}/t_{k−1} ≤ 1/25`
/// and `t_N ≤ (e n/N)^(2N) = e^{−6.09 n}`, so both tails are below `2^−(pw+16)` relative.
/// Up to `N` everything is exact integers: `A/B = V/(D·(Q+T))` from [`bm_split`] over
/// `1..=N`. The only roundings are that one division and `j·ln 2` (`ln2` within an ulp at
/// the argument precision `p`), each at the working precision `pw = p + 32` bits; with the
/// series truncation below `2^−(pw+14)` = `2^−(p+46)`, the total is a few ulps at `pw`, i.e.
/// under `2^−(p+24)` relative to the argument `p`.
///
/// There is no Ziv retry on this path: `Consts::euler_gamma` goes `exact_p` → `round_exact`,
/// which re-rounds the cached value once (`for_prec`, `set_precision`). Correct rounding
/// relies on a fixed margin instead: `Consts::euler_gamma_num` calls this with
/// `p = round_p(p_req) + 8·WORD_BIT_SIZE` (`+ 512` bits with 64-bit words, `+ 256` with
/// 32-bit), so the error is under `2^−(p_req+536)` (`2^−(p_req+280)` on 32-bit) relative to
/// the requested precision. Rounding could only go wrong if γ sat within that distance of a
/// `p_req`-bit rounding boundary.
/// Cost: O(M(p log p) log p) for M the multiplication cost (the old term-by-term sum did
/// ≈ 0.6 p full `p`-bit multiplications: 44 s at p = 65536).
pub(crate) fn euler_mascheroni(p: usize, ln2: &ExactNumNumber) -> Result<ExactNumNumber, Error> {
    let pw = p.checked_add(32).ok_or(Error::InvalidArgument)?;
    let mut j: u32 = 3;
    // 4·2^j·log2(e) ≥ pw + 16, with 4 log2 e < 5.771 (577/100 rounds it down: one more j at most).
    while (1usize << j)
        .checked_mul(577)
        .ok_or(Error::InvalidArgument)?
        < pw.checked_add(16)
            .and_then(|x| x.checked_mul(100))
            .ok_or(Error::InvalidArgument)?
    {
        j += 1;
    }
    let n_terms = 5u64
        .checked_shl(j)
        .filter(|x| x >> j == 5)
        .ok_or(Error::InvalidArgument)?;
    let two_j = Exponent::from(j) * 2;
    let s = bm_split(0, n_terms, two_j)?;
    let q = s.d.mul_full_prec(&s.d)?;
    // A/B = V / (D·(Q + T)) (the k = 0 term of B is 1, i.e. Q/Q).
    let den = s.d.mul_full_prec(&q.add_full_prec(&s.t)?)?;
    let ab = s.v.div(&den, pw, RoundingMode::None)?;
    let ln_n = ExactNumNumber::from_word(Word::from(j), pw)?.mul(ln2, pw, RoundingMode::None)?;
    ab.sub(&ln_n, pw, RoundingMode::None)
}

fn exact_mul_num(x: &ExactNumNumber, y: &ExactNumNumber) -> Result<ExactNumNumber, Error> {
    let w = round_p(x.mantissa_max_bit_len() + y.mantissa_max_bit_len() + WORD_BIT_SIZE);
    let v = x.mul(y, w, RoundingMode::None)?;
    if v.inexact() {
        return Err(Error::InvalidArgument);
    }
    trim_low_zero_words(v)
}

fn exact_add_num(x: &ExactNumNumber, y: &ExactNumNumber) -> Result<ExactNumNumber, Error> {
    if x.is_zero() {
        return y.clone();
    }
    if y.is_zero() {
        return x.clone();
    }
    let w_of = |v: &ExactNumNumber| Exponent::try_from(v.mantissa_max_bit_len()).unwrap_or(0);
    let ex = x.exponent().max(y.exponent());
    let lo = (x.exponent() - w_of(x)).min(y.exponent() - w_of(y));
    let w = round_p(usize::try_from(ex - lo).unwrap_or(0) + 2 * WORD_BIT_SIZE);
    let v = x.add(y, w, RoundingMode::None)?;
    if v.inexact() {
        return Err(Error::InvalidArgument);
    }
    trim_low_zero_words(v)
}

/// Drop all-zero low mantissa words (exact), so exact products cost by significant words.
fn trim_low_zero_words(mut v: ExactNumNumber) -> Result<ExactNumNumber, Error> {
    if v.is_zero() {
        return Ok(v);
    }
    let words = v.mantissa().digits();
    let zeros = words.iter().take_while(|w| **w == 0).count();
    if zeros == 0 || zeros >= words.len() {
        return Ok(v);
    }
    let q = (words.len() - zeros) * WORD_BIT_SIZE;
    v.set_precision(q, RoundingMode::ToZero)?;
    Ok(v)
}

/// Binary splitting for the terminating ₂F₁ sum `Σ_{k=lo}^{hi−1} Π_{j=lo}^{k} p_j / q_j`,
/// `p_j = (a+j−1)(b+j−1) z`, `q_j = (c+j−1) j`, all exact.
struct Hyp2f1Split<'a> {
    a: &'a ExactNumNumber,
    b: &'a ExactNumNumber,
    c: &'a ExactNumNumber,
    z: &'a ExactNumNumber,
}

impl Hyp2f1Split<'_> {
    /// `(P, Q, T)` with `P = Π p_j`, `Q = Π q_j`, `T / Q` = the partial sum over `[lo, hi)`.
    fn split(
        &self,
        lo: usize,
        hi: usize,
    ) -> Result<(ExactNumNumber, ExactNumNumber, ExactNumNumber), Error> {
        if hi - lo == 1 {
            let jm1 = ExactNumNumber::from_usize(lo - 1)?;
            let j = ExactNumNumber::from_usize(lo)?;
            let fa = exact_add_num(self.a, &jm1)?;
            let fb = exact_add_num(self.b, &jm1)?;
            let fc = exact_add_num(self.c, &jm1)?;
            let pj = exact_mul_num(&exact_mul_num(&fa, &fb)?, self.z)?;
            let qj = exact_mul_num(&fc, &j)?;
            let t = pj.clone()?;
            return Ok((pj, qj, t));
        }
        let mid = lo + (hi - lo) / 2;
        let (p1, q1, t1) = self.split(lo, mid)?;
        let (p2, q2, t2) = self.split(mid, hi)?;
        // T = T1 Q2 + P1 T2
        let t = exact_add_num(&exact_mul_num(&t1, &q2)?, &exact_mul_num(&p1, &t2)?)?;
        Ok((exact_mul_num(&p1, &p2)?, exact_mul_num(&q1, &q2)?, t))
    }
}

/// A 64-bit starting value for Lambert W: branch-point series `−1 ± q − q²/3 ± 11q³/72`,
/// `q = √(2(1 + e x))` (DLMF 4.13.6), for `x < −1/4`; Winitzki's `l (1 − ln(1+l)/(2+l))`,
/// `l = ln(1 + x)`, for `W₀` up to 3; `L₁ − L₂ + L₂/L₁` (DLMF 4.13.10) otherwise.
fn lambert_w_guess(
    x: &ExactNumNumber,
    minus1: bool,
    cc: &mut Consts,
) -> Result<ExactNumNumber, Error> {
    type N = ExactNumNumber;
    let p = 2 * WORD_BIT_SIZE;
    let nr = RoundingMode::None;
    let one = N::from_word(1, p)?;
    let quarter = {
        let mut q = N::from_word(1, p)?;
        q.set_exponent(-1);
        q.neg()?
    };
    if x.cmp(&quarter) < 0 {
        let e = cc.e_num(p, nr)?;
        let mut t = one.add(&e.mul(x, p, nr)?, p, nr)?;
        if t.is_negative() {
            t = N::new(p)?;
        }
        let q = N::from_word(2, p)?.mul(&t, p, nr)?.sqrt(p, nr)?;
        let q2 = q.mul(&q, p, nr)?.div(&N::from_word(3, p)?, p, nr)?;
        let q3 = q
            .mul(&q, p, nr)?
            .mul(&q, p, nr)?
            .mul(&N::from_word(11, p)?, p, nr)?;
        let q3 = q3.div(&N::from_word(72, p)?, p, nr)?;
        let m1 = one.neg()?;
        return if minus1 {
            m1.sub(&q, p, nr)?.sub(&q2, p, nr)?.sub(&q3, p, nr)
        } else {
            m1.add(&q, p, nr)?.sub(&q2, p, nr)?.add(&q3, p, nr)
        };
    }
    if !minus1 && x.cmp(&N::from_word(3, p)?) <= 0 {
        let l = x.log1p(p, nr, cc)?;
        if l.is_zero() {
            return x.try_rounded_to(p, nr);
        }
        let num = l.log1p(p, nr, cc)?;
        let den = N::from_word(2, p)?.add(&l, p, nr)?;
        return l.mul(&one.sub(&num.div(&den, p, nr)?, p, nr)?, p, nr);
    }
    let l1 = x.abs()?.ln(p, nr, cc)?;
    let l2 = l1.abs()?.ln(p, nr, cc)?;
    l1.sub(&l2, p, nr)?.add(&l2.div(&l1, p, nr)?, p, nr)
}

/// Lambert W to `p_wrk + 64` bits, then certify `|W − w| < 2^(exponent(w) − p_wrk − 64)`
/// by a sign change of the (monotone) defining equation; `None` when the bracket cannot be
/// certified at this precision.
fn lambert_w_bracketed(
    x: &ExactNumNumber,
    minus1: bool,
    guess: &ExactNumNumber,
    p_wrk: usize,
    cc: &mut Consts,
) -> Result<Option<ExactNumNumber>, Error> {
    type N = ExactNumNumber;
    let nr = RoundingMode::None;
    let one = N::from_word(1, WORD_BIT_SIZE)?;
    // Log form w + ln|w| = ln|x| for W₀(x > 1) and W₋₁(x >= −1/4): no e^w over/underflow.
    let log_form = if minus1 {
        let mut q = N::from_word(1, WORD_BIT_SIZE)?;
        q.set_exponent(-1);
        x.cmp(&q.neg()?) >= 0
    } else {
        x.cmp(&one) > 0
    };
    let opw = one.add(guess, 2 * WORD_BIT_SIZE, nr)?;
    let cond = if opw.is_zero() {
        0
    } else {
        usize::try_from((-opw.exponent()).max(0)).unwrap_or(0)
    };
    let pp = round_p(p_wrk + 2 * WORD_BIT_SIZE + cond + WORD_BIT_SIZE);
    let lnx = if log_form { Some(x.abs()?.ln(pp, nr, cc)?) } else { None };
    // f(w) and a bound on its evaluation error.
    let eval = |w: &N, cc: &mut Consts| -> Result<(N, N), Error> {
        let mut bound = N::from_word(1, WORD_BIT_SIZE)?;
        if let Some(l) = &lnx {
            let lw = w.abs()?.ln(pp, nr, cc)?;
            let f = w.add(&lw, pp, nr)?.sub(l, pp, nr)?;
            let mag = max_exp(&[w, &lw, l]);
            bound.set_exponent(mag + 4 - pp as Exponent);
            Ok((f, bound))
        } else {
            let ew = w.exp(pp, nr, cc)?;
            let wew = w.mul(&ew, pp, nr)?;
            let f = wew.sub(x, pp, nr)?;
            let mag = max_exp(&[&wew, x]);
            bound.set_exponent(mag + 4 - pp as Exponent);
            Ok((f, bound))
        }
    };
    let mut w = guess.clone()?;
    let mut converged = false;
    for _ in 0..(64 + 4 * usize::BITS as usize) {
        let (f, _) = eval(&w, cc)?;
        if f.is_zero() {
            converged = true;
            break;
        }
        let d = if lnx.is_some() {
            one.add(&one.div(&w, pp, nr)?, pp, nr)?
        } else {
            w.exp(pp, nr, cc)?.mul(&one.add(&w, pp, nr)?, pp, nr)?
        };
        if d.is_zero() {
            break;
        }
        let step = f.div(&d, pp, nr)?;
        w = w.sub(&step, pp, nr)?;
        if step.is_zero() || step.exponent() < w.exponent() - (pp as Exponent) + 8 {
            converged = true;
            break;
        }
    }
    if !converged || w.is_zero() {
        return Ok(None);
    }
    let mut delta = N::from_word(1, WORD_BIT_SIZE)?;
    delta.set_exponent(w.exponent() - (p_wrk + WORD_BIT_SIZE) as Exponent + 1);
    let lo = w.sub(&delta, pp, nr)?;
    let hi = w.add(&delta, pp, nr)?;
    let sign_of = |v: &N, cc: &mut Consts| -> Result<Option<bool>, Error> {
        let (f, b) = eval(v, cc)?;
        if f.is_zero() || f.abs()?.cmp(&b) <= 0 {
            return Ok(None);
        }
        Ok(Some(f.is_positive()))
    };
    let (Some(sl), Some(sh)) = (sign_of(&lo, cc)?, sign_of(&hi, cc)?) else {
        return Ok(None);
    };
    if sl == sh {
        return Ok(None);
    }
    Ok(Some(w))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The term cap the I/J series used before fb2809d (removed from the library: the series now runs to its tail
    /// bound). Kept here only so the regression test below can show the old cap stopped short of the I_1 peak.
    fn old_series_n_max(p: usize, exp: Exponent) -> usize {
        let mag = if exp <= 0 {
            0
        } else if exp >= 16 {
            p
        } else {
            1usize << (exp as usize).min(12)
        };
        p.saturating_add(32).saturating_add(mag)
    }

    /// `x` as an exact number at `p` bits: `2^e` (`e` may be negative).
    fn pow2(e: Exponent, p: usize) -> Result<ExactNumNumber, crate::Error> {
        let mut x = ExactNumNumber::from_word(1, p)?;
        x.set_exponent(e + 1);
        Ok(x)
    }

    /// `n / 4` exactly at `p` bits.
    fn quarters(n: Word, p: usize) -> Result<ExactNumNumber, crate::Error> {
        let mut x = ExactNumNumber::from_word(n, p)?;
        x.set_exponent(x.exponent() - 2);
        Ok(x)
    }

    fn exhausted<T>(r: Result<T, crate::Error>) -> bool {
        matches!(r, Err(crate::Error::PrecisionRetryExhausted))
    }

    /// Every capped loop in this file must end in PrecisionRetryExhausted when its cap
    /// (or divergence) stops it before its error bound holds; none may return a value.
    #[test]
    fn test_capped_loops_shift_caps() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let p = 256;
        // Γ: Stirling shift loop capped at 2 shifts, z = 0.75 stays far below e_min.
        assert!(exhausted(
            quarters(3, p)?.gamma_positive_capped(p, &mut cc, 2)
        ));
        // ψ: shift loop capped at 2 shifts, z = 1.25.
        assert!(exhausted(
            quarters(5, p)?.digamma_positive_capped(p, &mut cc, 2)
        ));
        Ok(())
    }

    /// Runs `shift_up` from `z` (target exponent 6, i.e. z >= 32) with the real z -> z+1 step,
    /// or a step that does not move z when `stuck`; returns the result and the steps taken. Bails
    /// out with `InvalidArgument` after 1000 steps so a missing cap fails instead of hanging.
    fn count_shift_up(z: Word, cap: usize, stuck: bool) -> (Result<(), crate::Error>, usize) {
        let p = 64;
        let mut steps = 0usize;
        let z = match ExactNumNumber::from_word(z, p) {
            Ok(z) => z,
            Err(e) => return (Err(e), 0),
        };
        let r = shift_up(&z, 6, cap, p, |z| {
            steps += 1;
            if steps > 1000 {
                return Err(crate::Error::InvalidArgument);
            }
            // a stuck step returns z unchanged (z + 0)
            let d = ExactNumNumber::from_word(if stuck { 0 } else { 1 }, p)?;
            z.add(&d, p, RoundingMode::None)
        });
        (r.map(|_| ()), steps)
    }

    /// `shift_up` fails at once, with no steps, when 2^(e_min-1) - z is above the cap; runs at
    /// exactly the cap; takes no step once z >= 2^(e_min-1); stops a step that does not move z
    /// after `cap` steps; and returns a step error as is. Before, a high precision ran all
    /// 2^20 full-width shifts before failing.
    #[test]
    fn test_shift_up_bounds() -> Result<(), crate::Error> {
        let px = crate::Error::PrecisionRetryExhausted;
        assert_eq!(count_shift_up(1, 10, false), (Err(px), 0));
        assert_eq!(count_shift_up(22, 10, false), (Ok(()), 10));
        assert_eq!(count_shift_up(22, 9, false), (Err(px), 0));
        assert_eq!(count_shift_up(31, 1, false), (Ok(()), 1));
        assert_eq!(count_shift_up(32, 0, false), (Ok(()), 0));
        assert_eq!(count_shift_up(40, 0, false), (Ok(()), 0));
        assert_eq!(count_shift_up(22, 10, true), (Err(px), 10));
        let z = ExactNumNumber::from_word(1, 64)?;
        let r = shift_up(&z, 6, 100, 64, |_| Err(crate::Error::MemoryAllocation));
        assert_eq!(r.map(|_| ()), Err(crate::Error::MemoryAllocation));
        // A target past the exponent range fails without stepping.
        let z = ExactNumNumber::from_word(1, 64)?;
        let r = shift_up(&z, Exponent::MAX, usize::MAX, 64, |_| {
            Err(crate::Error::InvalidArgument)
        });
        assert_eq!(r.map(|_| ()), Err(px));
        Ok(())
    }

    /// The Γ kernel at p = 1 << 14. With the Stirling term count scaled with p (as ψ scales its
    /// Bernoulli terms, 85b1964), the shift target from z = 0.75 is reachable under
    /// STIRLING_SHIFT_MAX, so Γ answers. With a cap below the shifts that target needs, the
    /// bail-out in `shift_up` refuses before the first step instead of running the loop at
    /// full width: that is what keeps a hopeless call fast.
    #[test]
    fn test_capped_loops_fail_fast_at_high_precision() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let p = 1 << 14;
        let x = quarters(3, p)?;
        // The value, not just Ok: Γ(3/4) at p = 2^14 against the FLINT (Arb) floor of
        // Γ(3/4)·2^16639 (16640 bits, error < 2^-16639; python-flint 0.9.0). A refusal or a wrong
        // value both fail here.
        let g = x.gamma_positive_capped(p, &mut cc, STIRLING_SHIFT_MAX)?;
        let reference = gamma_3_4_flint()?;
        let diff = g.sub(&reference, 2 * p, RoundingMode::None)?;
        let bound = 1 - Exponent::try_from(p).map_err(|_| crate::Error::InvalidArgument)? + 4;
        assert!(
            diff.is_zero() || diff.exponent() <= bound,
            "Γ(3/4) at p = 2^14 off by 2^{} (bound 2^{bound})",
            diff.exponent()
        );
        assert!(exhausted(x.gamma_positive_capped(p, &mut cc, 1)));
        Ok(())
    }

    /// Γ(3/4) at p=2^14 through the public `gamma` path (not the test-only
    /// capped API). A 64-term tail wants ~2^132 shifts and a 2^20 cap would
    /// refuse; the remainder-bound plan (more terms, `stirling_shift_count`)
    /// must finish and match the same FLINT floor as
    /// `test_capped_loops_fail_fast_at_high_precision`.
    #[test]
    fn gamma_three_quarters_at_high_p_from_error_bound() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let p = 1 << 14;
        let x = quarters(3, p)?;
        let t0 = std::time::Instant::now();
        let g = x.gamma(p, RoundingMode::ToEven, &mut cc)?;
        assert!(
            crate::test_time::within(t0.elapsed(), std::time::Duration::from_secs(60)),
            "gamma(3/4) at p=2^14 took {:?}",
            t0.elapsed()
        );
        let reference = gamma_3_4_flint()?;
        let diff = g.sub(&reference, 2 * p, RoundingMode::None)?;
        let bound = 1 - Exponent::try_from(p).map_err(|_| crate::Error::InvalidArgument)? + 4;
        assert!(
            diff.is_zero() || diff.exponent() <= bound,
            "Γ(3/4) at p=2^14 off by 2^{} (bound 2^{bound})",
            diff.exponent()
        );
        Ok(())
    }

    /// FLINT `floor(Γ(3/4)·2^16639)` as hex (16640 bits; Γ(3/4) ∈ [1, 2), so exponent 1).
    fn gamma_3_4_flint() -> Result<ExactNumNumber, crate::Error> {
        const HEX: [&str; 33] = [
        "9cda745a87b16013d7a9256698c602bbb7a79f91e7e526fd3294e130ed0bb0ffa823fffd6e9fcd25e50324be1b46ec4920dd2eac284104f194ca9180f09ebb02",
        "4313dee21938a703ec222d56b9e89d61e6f3d18ca14b9392576f3029bbc750cff48e75a5a96c658b1db227d86afab1cafeeaf49ffc1d6872cef621be9c1f1198",
        "276666d6b33e03c26bb78d76de369837a13034ee3498c7543aaa70576dde701476953e0e343620abdb49c6b128dcfc512a0be1ec75ec79bc0f191c7099ed0ad8",
        "3839c218103575089a9f42567a4c298d5f0329dec90d162e1d2dd2a968e434d4bd03df078f548de971e474f1db48dff571303ff0b49cd5000339bb10c86a7b76",
        "a6e3270e8a33addd82dc1be2dc0ff82dfc8ed76895db5fb2a9be35860f56a7fc8a5f68b09e54756ddee6907f4fe1aad39fae1d3369f2597fe5bba677c0ad72cf",
        "4ab3b6b715c449247ff4b84b9ab9d1466bf4fd09fcba4c91cc17d713697fc3849ffa1ea3ae427212fba7128c2706f99ab86c5df97b4a92cbc3d89e8c797a7b95",
        "669d325342170a7ae4d7049d865c74ba6523df0b87e500dd041e5d9a5f7211e56bd15ebd52336d2e13f6eed11c62690abf3e49d196d653c04ce4214b2fdee37f",
        "6bedafb128dc60ed83bbff264a35d8003ad16ca6eda32a5bc82b90fdf626430e44f141aa576adb0919f88844dc43a746eb1d8ac7132340e3c88ce18568c37532",
        "1e36f4a53359e77b829cab4378e32a36d4eb14088a30762b3e005045d880955f4dfef91b9bee6826d329c55799b07fe0b136c5e4fec81661de012c0341687b20",
        "e3e2c6adea4a7ad206c89cc142429440f26edcbd4f74b8a8b621024ece979821794d5741c07d091ca54eae139c54b061ab3dc5ab460e241356d34355330b252a",
        "7915e9d9d27a3b9a9112264632a1bec429297c586bdc80bb3ef4ce5b14729b5b064ded5f2b2130c9b4c62fb13e68dff3a57bddd195264a91b3e403f3b55ea837",
        "4dd57152b150146e1fe023be4d0d22191ed90fae06c062cb450044a044abd86d57082b5b8744259b48eac3919a8b95bd87af8626cd9babb4b076cc79e3fa22c7",
        "7e6fb365778710392efb11d215c0423d556a3cba415eaddb12de398e98f4a6b9f081e0cbd681eb9a35c6539a54759d7456933ef44142c5e71814c9231cb5dc46",
        "b37e981faf3cb603fe407e16d92fc164a965f48f9ca7da756a68f9f9a5d71714cad01260390f3b8af4587103071cf8a179bd3ee7f0f33768a7300e47d369a137",
        "f7761a4e64728708f44634afbbae9acdb5a90a293480070b6fce07b27dde8b2bf61762125a34ec423fa48715bb3d64488a1765efd149fe033caf52f99102dbe1",
        "2f6b11619ce1535f01ccde7d74c09ba0edc1471922305eb3da713ab45acf4b2e6fb786d461c56460a784b739508893e06fcc0f5e04304376d88dac9de0f4713b",
        "f554a260dd6b063d0ee8a94f032bb2c1d2abd4d7ee8e3c4054e5b2a434d45e1c2e2b1c43fd569921180e14e0c157539fbb55b88ef9ea632dbd854b0cf92ddaf2",
        "9736972845a1a9a094b937178700013de039e2045c9897dc8c0191ddb4c5b5ebed6535a9d82a1e688658b3fb00a7d54bd2189ff7e5ddccde5c78b236f29128db",
        "820602a73dd1f577fa518dd19ed2bb067c3792e0f7e479993efd28c3ddac13586877548e92e7c36701cf7511381641a22cd3fe05a9d517dd090a53cb0d1a071d",
        "76b9d4b7f5882dc665d8588363c10dd7820923bace586fe8b6a79a37464bd4f5842f7f95f3397386cfa229d0856f5ef051100c0809fab65964c6cb110e9b9b79",
        "bc9dbd1012f801c506c2143d1b4d9ff4416e2b09d6ad98aed22d336a4250311d1d0284058f801050b28af932fe102e7e606175baf946393e4dadc8a0cb86e0ee",
        "1d02d65852cae3b294223e8cf0abe0fcc690c580cd94ccd54944887c6f95f90faaf7d2810c5a0708e266a5fb614ce4fdf3656eff495a24158e06f06ce82840d4",
        "7c45e7fd5a59557f35c6c0e003971a4dd647163774d66ba51d71e2de66642be69809bdaa1bf99fd7b7bcfba1b93d307f7a45abbfc1a0491fc77212a7a1f2ec3e",
        "b278f9f2d0caaf3186d3ce98fef71e4053caa24a68d883a3f651fbf059d72ccb2f69c2e70b94332f655eb1a122d91a6f877632801dbd3a0731a882f5b880648b",
        "456b7aa0d94d6068a46c89d48a2f249b136c353ffb85d118bb61a65fc021126fc2528832e5d0b859eabc6c9355353d4b856869f7f3687a87b451b6f56db8a7aa",
        "28cb6911bb2657b8a05b1430727c01a728b929e7d3fb5297cb1e858fe4ed4d490736845a607282f4a87ea0c320c989563c76cb9170c435c85cb4187129be935a",
        "d27000e043bba9e5f5f4e24d5c57169ef3a1f5965dd2d7379d4fcb63aeb6167dde65aa6aa3366facf3b479f4503041d84278f0c4400e8fe4788745d9c2b1bd3b",
        "261438dbe50176c6cd3cb20ec98f6d53200007230b40b702f87443d3dd6cbd0d421e15ec0f487224147c56c8b94b67f223b31fb84c7add0be06856d66609e6d9",
        "7d8248211785982e4c45f577b9f19660786de513ab70b0e2d5d667592129c151750f5b76e9799fc8e37598474240cc1abf7d151ee4a87b4b8e92c0d0854e7978",
        "e5c814a738e58ffb8b500b74efcb061c6a15ece6f5464dbb72a75656af003d32e85103529fee2c9cf96d1845cdda482aa8f2d871375769ad4c9d18c0b1e85644",
        "3b0d1c7c556a29fa0650cb40ad276e9c82107e46169dab16eb7b44191c592bce5eb06090124762c70d41f2f5cfcfb8ed8f198bbe83d8db2169a670070ad0ff59",
        "d40200342319b093a74ab4fb2609b1361c6247571f1adfb90512c957aea9c999089dbdd3be38d6e12febf5b40a7ab66d08db2138746aa9dcd2dc6aeef00b1bae",
        "48bc36315a250833c87f7d1d87a13c60a4aeaef30ffc55eb56a9d91eadc55ce3",
        ];
        let h: alloc::string::String = HEX.concat();
        let mut words: Vec<Word> = Vec::new();
        let mut end = h.len();
        while end > 0 {
            let start = end.saturating_sub(16);
            words.push(
                Word::from_str_radix(&h[start..end], 16)
                    .map_err(|_| crate::Error::InvalidArgument)?,
            );
            end = start;
        }
        ExactNumNumber::from_raw_parts(&words, 16640, crate::Sign::Pos, 1, false)
    }

    #[test]
    fn test_capped_loops_asymptotic_tails() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let p = 4096;
        let one = ExactNumNumber::from_word(1, p)?;
        let nu0 = ExactNumNumber::from_word(0, p)?;
        // Fixed Bernoulli tables evaluated at z = 1, far outside their bound.
        assert!(exhausted(one.ln_gamma_stirling(p, &mut cc)));
        assert!(exhausted(one.digamma_asymp(p, 64, &mut cc)));
        // Divergent asymptotic series at x = 1: the smallest term is far above 2^-p.
        assert!(exhausted(one.erfc_asymptotic(p, &mut cc)));
        assert!(exhausted(one.ei_asymptotic(p, &mut cc)));
        assert!(exhausted(one.si_ci_aux_fg(p)));
        assert!(exhausted(one.fresnel_aux_fg(p, &mut cc)));
        assert!(exhausted(one.bessel_hankel_pq(&nu0, p)));
        assert!(exhausted(one.k_asymptotic(&nu0, p, &mut cc)));
        Ok(())
    }

    #[test]
    fn test_capped_loops_power_series() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let p = 64;
        // ZF-37: the Ei series has no term cap (all terms positive; it stops once k > x and the
        // term is below 2^-p of the largest). At x = 2^9 it needs ~1400 terms, past the old
        // ~2p cap, and now completes and agrees with Ei.
        let mid = pow2(9, p)?;
        let (series, _) = mid.ei_series(p + 64, &mut cc)?;
        let ei = mid.ei(p, RoundingMode::ToEven, &mut cc)?;
        bits_agree(&series, &ei, p, 56, "Ei(2^9) series")?;
        // x = 2^9 is past the old series_n_max at p=64 (nmax = p+32+256).
        // Production Si/Ci add cancellation-guard bits; the uncapped series at
        // that width must agree with the public kernel, not stop early.
        let pe = p + ExactNumNumber::si_ci_series_guard(mid.exponent());
        let si_s = mid.si_series(pe)?;
        bits_agree(
            &si_s,
            &mid.si(p, RoundingMode::ToEven, &mut cc)?,
            p,
            56,
            "Si(2^9) series",
        )?;
        let ci_s = mid.ci_series(pe, &mut cc)?;
        bits_agree(
            &ci_s,
            &mid.ci(p, RoundingMode::ToEven, &mut cc)?,
            p,
            56,
            "Ci(2^9) series",
        )?;
        let s1 = ExactNumNumber::from_word(1, 4096)?;
        // γ(1, 2^30) = 1 − e^{−2^30} used to hit the 2^28 series-term cap.
        let x30 = pow2(30, 4096)?;
        let got = s1.gammainc_at(&x30, 4096, &mut cc)?;
        let nx = x30.with_sign(Sign::Neg)?;
        let want = ExactNumNumber::from_word(1, 4096)?.sub(
            &nx.exp(4096, RoundingMode::None, &mut cc)?,
            4096,
            RoundingMode::None,
        )?;
        bits_agree(&got, &want, 4096, 4000, "γ(1, 2^30)")?;
        Ok(())
    }

    /// The cached exact B_2k equal sympy `bernoulli(2k)` (sympy 1.14.0), checked by an exact
    /// cross-multiplication; no timing involved.
    #[test]
    fn test_bernoulli_rev_vs_boustrophedon_modest() -> Result<(), String> {
        for n in [16, 30, 40] {
            match bernoulli_rev::exact_even_bernoulli_rev(n) {
                Ok(rev) => {
                    let bou = exact_even_bernoulli_boustrophedon(n)?;
                    assert_eq!(rev.len(), bou.len());
                    for k in 1..=n {
                        let (rn, rd) = &rev[k - 1];
                        let (bn, bd) = &bou[k - 1];
                        let lhs = rn.mul(bd);
                        let rhs = rd.mul(bn);
                        assert!(
                            lhs.soft_cmp(&rhs) == core::cmp::Ordering::Equal,
                            "k={k} n={n}: rev differs from boustrophedon"
                        );
                    }
                }
                Err(e) => return Err(format!("exact_even_bernoulli_rev({n}) failed: {e:?}")),
            }
        }
        Ok(())
    }

    #[test]
    fn test_bernoulli_cache_exact_vs_sympy() -> Result<(), crate::Error> {
        let rows: [(usize, &str, &str); 8] = [
        (1, "1", "6"),
        (2, "-1", "30"),
        (3, "1", "42"),
        (6, "-691", "2730"),
        (10, "-174611", "330"),
        (30, "-1215233140483755572040304994079820246041491", "56786730"),
        (64, "-267754707742548082886954405585282394779291459592551740629978686063357792734863530145362663093519862048495908453718017", "510"),
        (100, "-498384049428333414764928632140399662108495887457206674968055822617263669621523687568865802302210999132601412697613279391058654527145340515840099290478026350382802884371712359337984274122861159800280019110197888555893671151", "1366530"),
        ];
        with_exact_bernoulli(100, |rats| {
            for (k, num_s, den_s) in rows {
                let (num, den) = rats.get(k - 1).ok_or(crate::Error::InvalidArgument)?;
                let sn = ExactInt::from_dec_string(num_s).ok_or(crate::Error::InvalidArgument)?;
                let sd = ExactInt::from_dec_string(den_s).ok_or(crate::Error::InvalidArgument)?;
                let lhs = num.mul(&sd);
                let rhs = den.mul(&sn);
                assert!(
                    lhs.soft_cmp(&rhs) == core::cmp::Ordering::Equal,
                    "k={k}: B_2k differs from sympy"
                );
            }
            Ok(())
        })
    }

    /// The zeta path (Arb `arb_bernoulli_ui_zeta`) agrees with the exact cache
    /// at p=256, including indices that later use only the zeta route.
    #[test]
    fn test_bernoulli_zeta_matches_exact_cache() -> Result<(), crate::Error> {
        let p = 256;
        with_even_bernoulli(80, p, |exact| {
            for k in 20..=80 {
                let z = even_bernoulli_via_zeta(k, p)?;
                let e = &exact[k - 1];
                assert!(z.cmp(e) == 0, "k={k}: zeta path differs from exact cache");
            }
            Ok(())
        })
    }

    /// Rounded values past the old `BERNOULLI_EXACT_K = 4096` cap come from
    /// the exact cache, not the working-precision zeta walker, and still
    /// match the zeta-series midpoint at modest precision.
    #[test]
    fn zf_even_bernoulli_past_old_exact_cap_matches_zeta() -> Result<(), crate::Error> {
        let p = 256;
        let k = 4112;
        with_even_bernoulli(k, p, |v| {
            assert_eq!(v.len(), k);
            let z = even_bernoulli_via_zeta(k, p)?;
            assert!(
                v[k - 1].cmp(&z) == 0,
                "B_{} at p={p} from exact cache differs from zeta path",
                2 * k
            );
            Ok(())
        })
    }

    /// FLINT 3.6 exact `B_n` (num/den, including B_0, B_1, B_2 and odd n)
    /// plus working-precision midpoints at dest p = 64, 256, and 4096.
    #[test]
    fn zf_bernoulli_matches_flint_golds() -> Result<(), crate::Error> {
        // Cap at B_4096 so debug unit tests do not fill exact B_21192.
        // Larger indices are checked in the p=65536 gold test (release).
        bernoulli_flint_golds_at_capped(&[64, 256, 4096], 4096)
    }

    /// FLINT 3.6 working-precision midpoints at dest p = 65536, including
    /// the largest even index used by the Barnes G EM tail (`B_21192`).
    #[test]
    fn zf_bernoulli_matches_flint_golds_p65536() -> Result<(), crate::Error> {
        bernoulli_flint_golds_at(&[65536])
    }

    #[test]
    fn zf_bernoulli_edge_cases() -> Result<(), crate::Error> {
        let p = 256usize;
        let b0 = bernoulli_ui(0, p)?;
        let one = ExactNumNumber::from_word(1, p)?;
        assert!(b0.cmp(&one) == 0, "B_0 must be 1");

        let b1 = bernoulli_ui(1, p)?;
        let two = ExactNumNumber::from_word(2, p)?;
        let half = one.div(&two, p, RoundingMode::ToEven)?.neg()?;
        assert!(b1.cmp(&half) == 0, "B_1 must be -1/2");

        let b2 = bernoulli_ui(2, p)?;
        let six = ExactNumNumber::from_word(6, p)?;
        let sixth = one.div(&six, p, RoundingMode::ToEven)?;
        assert!(b2.cmp(&sixth) == 0, "B_2 must be 1/6");

        for n in [3usize, 5, 7, 11, 21, 101] {
            let b = bernoulli_ui(n, p)?;
            assert!(b.is_zero(), "B_{n} for odd n>1 must be 0");
        }
        Ok(())
    }

    fn bernoulli_flint_golds_at(precs: &[usize]) -> Result<(), crate::Error> {
        bernoulli_flint_golds_at_capped(precs, usize::MAX)
    }

    /// `max_n` caps the arb `B_n` rows so the default (debug) gold test does
    /// not first-fill the exact cache through `B_21192`.
    fn bernoulli_flint_golds_at_capped(precs: &[usize], max_n: usize) -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let raw = include_str!("../testdata/bernoulli_flint_golds.txt");
        for line in raw.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut parts = line.splitn(4, ' ');
            let kind = parts.next().ok_or(crate::Error::InvalidArgument)?;
            if kind == "exact" {
                let n: usize = parts
                    .next()
                    .ok_or(crate::Error::InvalidArgument)?
                    .parse()
                    .map_err(|_| crate::Error::InvalidArgument)?;
                let frac = parts.next().ok_or(crate::Error::InvalidArgument)?;
                check_exact_bernoulli_gold(n, frac)?;
            } else if kind == "arb" {
                let p: usize = parts
                    .next()
                    .ok_or(crate::Error::InvalidArgument)?
                    .parse()
                    .map_err(|_| crate::Error::InvalidArgument)?;
                if !precs.contains(&p) {
                    continue;
                }
                let n: usize = parts
                    .next()
                    .ok_or(crate::Error::InvalidArgument)?
                    .parse()
                    .map_err(|_| crate::Error::InvalidArgument)?;
                if n > max_n {
                    continue;
                }
                let dec = parts.next().ok_or(crate::Error::InvalidArgument)?;
                let got = bernoulli_ui(n, p)?;
                let want = ExactNumNumber::parse(
                    dec,
                    crate::defs::Radix::Dec,
                    4 * p,
                    RoundingMode::ToEven,
                    &mut cc,
                )?;
                let min = Exponent::try_from(p.saturating_sub(16)).unwrap_or(16);
                bits_agree(&got, &want, p, min, &format!("B_{n} p={p}"))?;
            } else {
                return Err(crate::Error::InvalidArgument);
            }
        }
        Ok(())
    }

    fn check_exact_bernoulli_gold(n: usize, frac: &str) -> Result<(), crate::Error> {
        if n == 0 {
            let b = bernoulli_ui(0, 64)?;
            let one = ExactNumNumber::from_word(1, 64)?;
            assert!(b.cmp(&one) == 0, "exact B_0");
            return Ok(());
        }
        if n == 1 {
            let b = bernoulli_ui(1, 64)?;
            let one = ExactNumNumber::from_word(1, 64)?;
            let two = ExactNumNumber::from_word(2, 64)?;
            let half = one.div(&two, 64, RoundingMode::ToEven)?.neg()?;
            assert!(b.cmp(&half) == 0, "exact B_1");
            return Ok(());
        }
        if !n.is_multiple_of(2) {
            let b = bernoulli_ui(n, 64)?;
            assert!(b.is_zero(), "exact B_{n} odd");
            return Ok(());
        }
        let (num_s, den_s) =
            if let Some((a, b)) = frac.split_once('/') { (a, b) } else { (frac, "1") };
        let k = n / 2;
        with_exact_bernoulli(k, |rats| {
            let (num, den) = rats.get(k - 1).ok_or(crate::Error::InvalidArgument)?;
            let sn = ExactInt::from_dec_string(num_s).ok_or(crate::Error::InvalidArgument)?;
            let sd = ExactInt::from_dec_string(den_s).ok_or(crate::Error::InvalidArgument)?;
            let lhs = num.mul(&sd);
            let rhs = den.mul(&sn);
            assert!(
                lhs.soft_cmp(&rhs) == core::cmp::Ordering::Equal,
                "exact B_{n} differs from FLINT"
            );
            Ok(())
        })
    }

    /// The Bloemen range walker agrees with the single-`k` zeta path, on
    /// both the Euler-product side (`n ≤ 0.7 wp`) and the short-tail side.
    #[test]
    fn test_bernoulli_zeta_range_matches_single() -> Result<(), crate::Error> {
        let p = 256;
        for &(k0, k1) in &[(40, 80), (120, 140)] {
            let range = even_bernoulli_zeta_range(k0, k1, p)?;
            for k in k0..=k1 {
                let z = even_bernoulli_via_zeta(k, p)?;
                assert!(
                    range[k - k0].cmp(&z) == 0,
                    "k={k}: range walker differs from single-k zeta"
                );
            }
        }
        Ok(())
    }

    /// First-fill of the exact even-Bernoulli table at the improved-Stirling
    /// head size for dest 65536 (M ≈ 3316). Ignored: used to time the
    /// `bernoulli_rev` port in official rustc release.
    #[test]
    fn time_exact_bernoulli_rev_3316() -> Result<(), crate::Error> {
        let t0 = std::time::Instant::now();
        with_exact_bernoulli(3316, |rats| {
            assert_eq!(rats.len(), 3316);
            let (n, d) = &rats[3315];
            assert!(!n.is_zero() && !d.is_zero(), "B_6632 zero");
            Ok(())
        })?;
        eprintln!(
            "exact_even_bernoulli_rev(3316) first fill {:?}",
            t0.elapsed()
        );
        let t1 = std::time::Instant::now();
        with_exact_bernoulli(3316, |_| Ok(()))?;
        eprintln!(
            "exact_even_bernoulli_rev(3316) cache hit {:?}",
            t1.elapsed()
        );
        Ok(())
    }

    /// Profile the exact cache fill at the EM size for dest p=65536
    /// (`M ≈ 10596`, `B_2..B_21192`) and a subsequent rounded lookup.
    /// Ignored: long run; has a progress log from `bernoulli_rev`.
    #[test]
    fn time_exact_bernoulli_cache_em_size() -> Result<(), crate::Error> {
        let k = 10596;
        let p = 65536;
        eprintln!("BERN profile phase=exact_fill k={k} p={p} counter=0/1 ETA=unknown");
        let t0 = std::time::Instant::now();
        with_exact_bernoulli(k, |rats| {
            assert!(rats.len() >= k);
            Ok(())
        })?;
        let fill_ms = t0.elapsed().as_millis();
        eprintln!(
            "BERN profile phase=exact_fill_done k={k} ms={fill_ms} (FLINT B_0..B_21192 exact was 3.81 s)"
        );
        let t1 = std::time::Instant::now();
        with_even_bernoulli(k, p, |v| {
            assert_eq!(v.len(), k);
            Ok(())
        })?;
        let rnd_ms = t1.elapsed().as_millis();
        eprintln!(
            "BERN profile phase=round_done k={k} p={p} ms={rnd_ms} (FLINT arb_bernoulli_ui through 21192 at 65536 bits was 16.0 s standalone)"
        );
        Ok(())
    }

    #[test]
    fn test_erf_gamma_bessel() -> Result<(), crate::Error> {
        let p = 256;
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;

        let zero = ExactNumNumber::new(p)?;
        assert!(zero.erf(p, rm, &mut cc)?.is_zero());

        let one = ExactNumNumber::from_word(1, p)?;
        let g1 = one.gamma(p, rm, &mut cc)?;
        let d = g1.sub(&one, p, RoundingMode::None)?.abs()?;
        assert!(d.is_zero() || d.exponent() < -((p as Exponent) / 4));

        let five = ExactNumNumber::from_word(5, p)?;
        let g5 = five.gamma(p, rm, &mut cc)?;
        let tf = ExactNumNumber::from_word(24, p)?;
        let d = g5.sub(&tf, p, RoundingMode::None)?.abs()?;
        assert!(d.is_zero() || d.exponent() < -((p as Exponent) / 4));

        let j0 = zero.bessel_j(0, p, rm, &mut cc)?;
        let d = j0.sub(&one, p, RoundingMode::None)?.abs()?;
        assert!(d.is_zero() || d.exponent() < -((p as Exponent) / 4));

        let x = ExactNumNumber::from_word(1, p)?;
        let e = x.erf(p, rm, &mut cc)?;
        let en = x.neg()?.erf(p, rm, &mut cc)?;
        assert!(e.cmp(&en.neg()?) == 0);
        Ok(())
    }

    fn bits_agree(
        a: &ExactNumNumber,
        b: &ExactNumNumber,
        p: usize,
        min_bits: Exponent,
        label: &str,
    ) -> Result<(), crate::Error> {
        let d = a.sub(b, p, RoundingMode::None)?.abs()?;
        if d.is_zero() {
            return Ok(());
        }
        let rel = d.exponent() - a.exponent();
        assert!(
            rel < -min_bits,
            "{label}: relative exponent {rel} (a.exp={}, b.exp={})",
            a.exponent(),
            b.exponent()
        );
        Ok(())
    }

    #[test]
    fn test_ei_si_li_fresnel() -> Result<(), crate::Error> {
        let p = 256;
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;
        let one = ExactNumNumber::from_word(1, p)?;
        let zero = ExactNumNumber::new(p)?;

        assert!(one.ei(p, rm, &mut cc).is_ok());
        assert!(zero.ei(p, rm, &mut cc).is_err());
        assert!(one.neg()?.ei(p, rm, &mut cc).is_ok());

        let s0 = zero.si(p, rm, &mut cc)?;
        assert!(s0.is_zero());
        let s1 = one.si(p, rm, &mut cc)?;
        let sn = one.neg()?.si(p, rm, &mut cc)?;
        assert!(s1.cmp(&sn.neg()?) == 0);

        let e = ExactNumNumber::from_word(1, p)?.exp(p, rm, &mut cc)?;
        let li_e = e.li(p, rm, &mut cc)?;
        let ei1 = one.ei(p, rm, &mut cc)?;
        bits_agree(&li_e, &ei1, p, (p as Exponent) / 4, "li(e)=Ei(1)")?;

        assert!(one.li(p, rm, &mut cc).is_err());
        assert!(zero.ci(p, rm, &mut cc).is_err());
        assert!(one.ci(p, rm, &mut cc).is_ok());

        assert!(zero.fresnel_s(p, rm, &mut cc)?.is_zero());
        assert!(zero.fresnel_c(p, rm, &mut cc)?.is_zero());
        let fs = one.fresnel_s(p, rm, &mut cc)?;
        let fns = one.neg()?.fresnel_s(p, rm, &mut cc)?;
        assert!(fs.cmp(&fns.neg()?) == 0);

        let forty = ExactNumNumber::from_word(40, p)?;
        bits_agree(
            &forty.si(p, rm, &mut cc)?,
            &forty.si(128, rm, &mut cc)?,
            p,
            80,
            "Si(40)",
        )?;
        bits_agree(
            &forty.ci(p, rm, &mut cc)?,
            &forty.ci(128, rm, &mut cc)?,
            p,
            80,
            "Ci(40)",
        )?;

        let three_hundred = ExactNumNumber::from_word(300, p)?;
        let si300 = three_hundred.si(p, rm, &mut cc)?;
        let mut half_pi = cc.pi_num(p, RoundingMode::None)?;
        half_pi.div_by_2(RoundingMode::None);
        let gap = half_pi.sub(&si300, p, RoundingMode::None)?.abs()?;
        assert!(gap.exponent() < -6, "Si(300) within O(1/x) of π/2");
        bits_agree(
            &si300,
            &three_hundred.si(128, rm, &mut cc)?,
            p,
            80,
            "Si(300)",
        )?;

        let sixty_four = ExactNumNumber::from_word(64, p)?;
        bits_agree(
            &sixty_four.si(p, rm, &mut cc)?,
            &sixty_four.si(64, rm, &mut cc)?,
            p,
            40,
            "Si(64) series vs f,g",
        )?;
        bits_agree(
            &sixty_four.ci(p, rm, &mut cc)?,
            &sixty_four.ci(64, rm, &mut cc)?,
            p,
            40,
            "Ci(64) series vs f,g",
        )?;
        let eight = ExactNumNumber::from_word(8, p)?;
        bits_agree(
            &eight.fresnel_s(p, rm, &mut cc)?,
            &eight.fresnel_s(64, rm, &mut cc)?,
            p,
            40,
            "S(8) series vs f,g",
        )?;
        bits_agree(
            &three_hundred.ei(p, rm, &mut cc)?,
            &three_hundred.ei(128, rm, &mut cc)?,
            p,
            80,
            "Ei(300)",
        )?;

        let twenty = ExactNumNumber::from_word(20, p)?;
        let fs20 = twenty.fresnel_s(p, rm, &mut cc)?;
        let fc20 = twenty.fresnel_c(p, rm, &mut cc)?;
        let half = one_half(p)?;
        let ds = half.sub(&fs20, p, RoundingMode::None)?.abs()?;
        let dc = half.sub(&fc20, p, RoundingMode::None)?.abs()?;
        assert!(ds.exponent() < -4, "S(20) near 1/2");
        assert!(dc.exponent() < -4, "C(20) near 1/2");
        bits_agree(&fs20, &twenty.fresnel_s(128, rm, &mut cc)?, p, 80, "S(20)")?;

        let half_ext = crate::ExactNum::from_u8(1, p).div(&crate::ExactNum::from_u8(2, p), p, rm);
        let inf_si = crate::INF_POS.si(p, rm, &mut cc);
        let hp = cc.pi(p, rm).div(&crate::ExactNum::from_u8(2, p), p, rm);
        let d = inf_si.sub(&hp, p, rm).abs();
        assert!(d.is_zero() || d.exponent().unwrap_or(0) < -((p as Exponent) / 4));
        let d = crate::INF_POS
            .fresnel_s(p, rm, &mut cc)
            .sub(&half_ext, p, rm)
            .abs();
        assert!(d.is_zero() || d.exponent().unwrap_or(0) < -((p as Exponent) / 4));
        let d = crate::INF_POS
            .fresnel_c(p, rm, &mut cc)
            .sub(&half_ext, p, rm)
            .abs();
        assert!(d.is_zero() || d.exponent().unwrap_or(0) < -((p as Exponent) / 4));
        Ok(())
    }

    #[test]
    fn test_digamma_gammainc() -> Result<(), crate::Error> {
        let p = 256;
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;
        let one = ExactNumNumber::from_word(1, p)?;
        let two = ExactNumNumber::from_word(2, p)?;
        let zero = ExactNumNumber::new(p)?;
        let g = cc.euler_gamma_num(p, RoundingMode::None)?;
        let half = one_half(p)?;

        assert!(zero.digamma(p, rm, &mut cc).is_err());
        let neg_half = half.neg()?;
        let psi_nh = neg_half.digamma(p, rm, &mut cc)?;
        // ψ(1/2) = −γ − 2 ln 2; ψ(−1/2) = ψ(1/2) − π cot(−π/2) wait: 1−(−1/2)=3/2
        // ψ(−1/2) = ψ(3/2) − π cot(−π/2). cot(−π/2)=0, ψ(3/2)=ψ(1/2)+2
        let three_half = ExactNumNumber::from_word(3, p)?.div(&two, p, RoundingMode::None)?;
        let psi_32 = three_half.digamma(p, rm, &mut cc)?;
        bits_agree(&psi_nh, &psi_32, p, 40, "psi(-1/2)=psi(3/2)")?;
        let psi1 = one.digamma(p, rm, &mut cc)?;
        let r1 = psi1.add(&g, p, RoundingMode::None)?.abs()?;
        assert!(r1.is_zero() || r1.exponent() < -80, "psi(1)+γ");

        let psi2 = two.digamma(p, rm, &mut cc)?;
        let r2 = psi2
            .add(&g, p, RoundingMode::None)?
            .sub(&one, p, RoundingMode::None)?
            .abs()?;
        assert!(r2.is_zero() || r2.exponent() < -80, "psi(2)+γ-1");

        let psih = half.digamma(p, rm, &mut cc)?;
        let ln2 = cc.ln_2_num(p, RoundingMode::None)?;
        let two_ln2 = two.mul(&ln2, p, RoundingMode::None)?;
        let want = g.neg()?.sub(&two_ln2, p, RoundingMode::None)?;
        bits_agree(&psih, &want, p, 80, "psi(1/2)")?;

        let five = ExactNumNumber::from_word(5, p)?;
        assert!(five.gammainc(&zero, p, rm, &mut cc)?.is_zero());
        let g11 = one.gammainc(&one, p, rm, &mut cc)?;
        let en = one.neg()?.exp(p, rm, &mut cc)?;
        let want = one.sub(&en, p, RoundingMode::None)?;
        bits_agree(&g11, &want, p, 80, "γ(1,1)")?;
        let gu = one.gammainc_upper(&one, p, rm, &mut cc)?;
        bits_agree(&gu, &en, p, 80, "Γ(1,1)=e^{-1}")?;
        assert!(one.gammainc(&one.neg()?, p, rm, &mut cc).is_err());
        Ok(())
    }

    #[test]
    fn test_bessel_family() -> Result<(), crate::Error> {
        let p = 256;
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;
        let one = ExactNumNumber::from_word(1, p)?;
        let zero = ExactNumNumber::new(p)?;
        let half = one_half(p)?;
        let two = ExactNumNumber::from_word(2, p)?;

        let j0 = zero.bessel_j_nu(&zero, p, rm, &mut cc)?;
        bits_agree(&j0, &one, p, 80, "J_0(0)")?;

        let i0 = zero.bessel_i(&zero, p, rm, &mut cc)?;
        bits_agree(&i0, &one, p, 80, "I_0(0)")?;

        assert!(zero.bessel_y(&zero, p, rm, &mut cc).is_err());
        assert!(zero.bessel_k(&zero, p, rm, &mut cc).is_err());

        let k_half = one.bessel_k(&half, p, rm, &mut cc)?;
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let mut nx = ExactNumNumber::from_word(1, p)?;
        nx.set_sign(Sign::Neg);
        let want = pi
            .div(&two, p, RoundingMode::None)?
            .sqrt(p, RoundingMode::None)?
            .mul(&nx.exp(p, rm, &mut cc)?, p, RoundingMode::None)?;
        bits_agree(&k_half, &want, p, 80, "K_{1/2}(1)")?;

        let k_neg = one.bessel_k(&half.neg()?, p, rm, &mut cc)?;
        bits_agree(&k_half, &k_neg, p, 80, "K_{-1/2}=K_{1/2}")?;

        let k0 = one.bessel_k(&zero, p, rm, &mut cc)?;
        let k1 = one.bessel_k(&one, p, rm, &mut cc)?;
        let i0 = one.bessel_i(&zero, p, rm, &mut cc)?;
        let i1 = one.bessel_i(&one, p, rm, &mut cc)?;
        let w = i0.mul(&k1, p, RoundingMode::None)?.add(
            &i1.mul(&k0, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
        bits_agree(&w, &one, p, 40, "Wronskian at 1")?;

        let j2000 = zero.bessel_j(2000, p, rm, &mut cc)?;
        assert!(j2000.is_zero(), "J_2000(0)=0");
        let k50 = one.bessel_k(&ExactNumNumber::from_word(50, p)?, p, rm, &mut cc)?;
        let k49 = one.bessel_k(&ExactNumNumber::from_word(49, p)?, p, rm, &mut cc)?;
        let k51 = one.bessel_k(&ExactNumNumber::from_word(51, p)?, p, rm, &mut cc)?;
        // K_{ν+1} = (2ν/x) K_ν + K_{ν−1}
        let rec = ExactNumNumber::from_word(100, p)?
            .mul(&k50, p, RoundingMode::None)?
            .add(&k49, p, RoundingMode::None)?;
        bits_agree(&k51, &rec, p, 20, "K recurrence at 50")?;

        // Large-|x| Hankel: J_0(20000) class used to be InvalidArgument.
        let x20k = ExactNumNumber::from_word(20000, p)?;
        let j0_big = x20k.bessel_j(0, p, rm, &mut cc)?;
        let y0_big = x20k.bessel_y(&zero, p, rm, &mut cc)?;
        let j1_big = x20k.bessel_j(1, p, rm, &mut cc)?;
        let y1_big = x20k.bessel_y(&one, p, rm, &mut cc)?;
        let wr = j0_big.mul(&y1_big, p, RoundingMode::None)?.sub(
            &j1_big.mul(&y0_big, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
        let want_w = two
            .div(
                &pi.mul(&x20k, p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )?
            .neg()?;
        bits_agree(&wr, &want_w, p, 40, "J/Y Wronskian at 20000")?;
        Ok(())
    }

    /// The old `HALF_INTEGER_MAX = 2^12` cutoff rejected n = 4097. The index
    /// helper now accepts it whenever n ≤ 4p + 4096. The public Γ closed form
    /// at this n builds a ~n log n-bit product, so the value check runs in
    /// release only (debug was 110 s).
    #[test]
    fn test_gamma_half_integer_beyond_old_4096_cap() -> Result<(), crate::Error> {
        let p = 64;
        let half = one_half(p)?;
        let x = from_usize_p(4097, p)?.add(&half, p, RoundingMode::None)?;
        assert_eq!(x.half_integer_index(p)?, Some(4097));
        Ok(())
    }

    /// `bessel_cancel_bits` used to return `InvalidArgument` once the guard
    /// exceeded `BESSEL_GUARD_MAX = 2^14` (about |x| ≳ 10913 for the 3/2 rule).
    #[test]
    fn test_bessel_cancel_bits_no_fixed_guard_cap() -> Result<(), crate::Error> {
        let x = from_usize_p(20_000, 64)?;
        let g = bessel_cancel_bits(&x, 3, 2)?;
        assert!(g > (1 << 14), "guard for |x|=20000 was {g}, still capped?");
        Ok(())
    }

    /// The old `series_n_max` cap is O(p) once the exponent is ≥ 16, and even at exponent 14
    /// it is only p+32+2^12. For I_1(10000) the peak is at k ≈ 5000, past that
    /// old limit, so a cap there returns PrecisionRetryExhausted (or, before
    /// the `!done` check, a truncated sum). The series must run past the peak
    /// until the tail is below 2^{-p}.
    #[test]
    fn test_bessel_i_series_past_old_n_max() -> Result<(), crate::Error> {
        let p = 64;
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;
        let x = ExactNumNumber::from_word(10_000, p)?;
        let nmax = old_series_n_max(p, x.exponent());
        assert!(
            nmax < 5000,
            "old_series_n_max({p}, exp={}) = {nmax}, no longer short of the I-peak?",
            x.exponent()
        );
        let one = ExactNumNumber::from_word(1, p)?;
        let got = x.bessel_i(&one, p, rm, &mut cc)?;
        let gold = ExactNumNumber::parse(
            "3.51328038895374889521677916179094689012e+4340",
            crate::defs::Radix::Dec,
            p,
            rm,
            &mut cc,
        )?;
        bits_agree(&got, &gold, p, 56, "I_1(10000)")?;
        Ok(())
    }

    /// B(a, 2^{20}) at p=2048 used to take the exact integer-n product
    /// (n·per ≈ 2^{27} bits is still allocatable under MAX_PRECISION=2^{60}−1),
    /// then hang on ~2^{20} growing exact muls. It must finish and agree with
    /// the closed form / gold; B(−½, 2)=−4 stays exact.
    #[test]
    fn test_beta_large_integer_n_uses_ziv_not_huge_product() -> Result<(), crate::Error> {
        let p = 2048;
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;
        let n20 = ExactNumNumber::from_word(1u64 << 20, p)?;
        let two = ExactNumNumber::from_word(2, p)?;
        // B(2, 2^{20}) = 1 / (2^{20} (2^{20}+1)), exact rational.
        let got = two.beta(&n20, p, rm, &mut cc)?;
        let np1 = ExactNumNumber::from_word((1u64 << 20) + 1, p)?;
        let den = n20.mul(&np1, p, RoundingMode::None)?;
        let want = ExactNumNumber::from_word(1, p)?.div(&den, p, rm)?;
        bits_agree(&got, &want, p, 2000, "B(2, 2^20)")?;
        let half = one_half(p)?;
        let got_h = half.beta(&n20, p, rm, &mut cc)?;
        let gold_h = ExactNumNumber::parse(
            "0.00173091217011571577082584230373951918876295553419736",
            crate::defs::Radix::Dec,
            p,
            rm,
            &mut cc,
        )?;
        bits_agree(&got_h, &gold_h, p, 160, "B(1/2, 2^20)")?;
        let mone_half = half.neg()?;
        let bm = mone_half.beta(&two, p, rm, &mut cc)?;
        let mfour = ExactNumNumber::from_word(4, p)?.neg()?;
        assert_eq!(bm.cmp(&mfour), 0, "B(-1/2, 2) = -4");
        Ok(())
    }

    /// ZF-37: `_2F_1(1, -1/2; 1; 1-2^{-20}) = (1-z)^{1/2} = 2^{-10}` exactly
    /// (DLMF 15.4.6, `F(a,b;a;z)=(1-z)^{-b}`). Must not be NaN.
    #[test]
    fn test_hyp2f1_one_minus_pow2_is_exact_sqrt() -> Result<(), crate::Error> {
        let p = 128;
        let mut cc = Consts::new()?;
        let rm = RoundingMode::Up;
        let one = ExactNumNumber::from_word(1, p)?;
        let half = one_half(p)?;
        let mhalf = half.neg()?;
        let two20 = ExactNumNumber::from_word(1u64 << 20, p)?;
        let omz = one.div(&two20, p, RoundingMode::None)?;
        let z = one.sub(&omz, p, RoundingMode::None)?;
        let got = one.hypergeom_2f1(&mhalf, &one, &z, p, rm, &mut cc)?;
        let want = one.div(&ExactNumNumber::from_word(1u64 << 10, p)?, p, rm)?;
        assert_eq!(
            got.cmp(&want),
            0,
            "2F1(1, -1/2; 1; 1-2^-20) p=128 Up must be 2^-10"
        );
        Ok(())
    }

    /// Γ(4097.5) matches the closed form `√π · (2n−1)!! / 2^n` (56 bits).
    #[test]
    fn test_gamma_half_integer_4097_closed_form() -> Result<(), crate::Error> {
        let p = 64;
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;
        let n: usize = 4097;
        let half = one_half(p)?;
        let x = from_usize_p(n, p)?.add(&half, p, RoundingMode::None)?;
        let got = x.gamma(p, rm, &mut cc)?;
        let mut prod = ExactNumNumber::from_word(1, WORD_BIT_SIZE)?;
        for k in 1..=n {
            let pb = round_p(prod.mantissa_max_bit_len() + 2 * WORD_BIT_SIZE);
            prod = prod.mul(&from_usize_p(2 * k - 1, pb)?, pb, RoundingMode::None)?;
        }
        let e =
            prod.exponent() - Exponent::try_from(n).map_err(|_| crate::Error::InvalidArgument)?;
        prod.set_exponent(e);
        let pw = p + WORD_BIT_SIZE;
        let rpi = cc
            .pi_num(pw, RoundingMode::None)?
            .sqrt(pw, RoundingMode::None)?;
        let want = rpi.mul(&prod, p, RoundingMode::None)?;
        bits_agree(&got, &want, p, 56, "Γ(4097.5)")?;
        Ok(())
    }

    #[test]
    fn test_elliptic_integrals() -> Result<(), crate::Error> {
        let p = 256;
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;
        let one = ExactNumNumber::from_word(1, p)?;
        let zero = ExactNumNumber::new(p)?;
        let two = ExactNumNumber::from_word(2, p)?;
        let three = ExactNumNumber::from_word(3, p)?;
        let four = ExactNumNumber::from_word(4, p)?;
        let half = one_half(p)?;
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let mut half_pi = cc.pi_num(p, RoundingMode::None)?;
        half_pi.div_by_2(RoundingMode::None);

        let k0 = zero.elliptic_k(p, rm, &mut cc)?;
        bits_agree(&k0, &half_pi, p, 80, "K(0)")?;
        let e0 = zero.elliptic_e_complete(p, rm, &mut cc)?;
        bits_agree(&e0, &half_pi, p, 80, "E(0)")?;
        let e1 = one.elliptic_e_complete(p, rm, &mut cc)?;
        bits_agree(&e1, &one, p, 80, "E(1)")?;
        assert!(matches!(
            one.elliptic_k(p, rm, &mut cc),
            Err(Error::ExponentOverflow(Sign::Pos))
        ));

        let k_half = half.elliptic_k(p, rm, &mut cc)?;
        let q = ExactNumNumber::from_word(1, p)?.div(&four, p, RoundingMode::None)?;
        let k4 = four.elliptic_k(p, rm, &mut cc)?;
        let kq = q.elliptic_k(p, rm, &mut cc)?;
        let want4 = kq.div(&two, p, RoundingMode::None)?;
        bits_agree(&k4, &want4, p, 40, "K(4)=K(1/4)/2")?;
        let g14 = q.gamma(p, rm, &mut cc)?;
        let want_k = g14.mul(&g14, p, RoundingMode::None)?.div(
            &four.mul(&pi.sqrt(p, RoundingMode::None)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
        bits_agree(&k_half, &want_k, p, 80, "K(1/2)")?;

        let e_half = half.elliptic_e_complete(p, rm, &mut cc)?;
        let want_e = pi
            .div(
                &four.mul(&k_half, p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )?
            .add(
                &k_half.div(&two, p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )?;
        bits_agree(&e_half, &want_e, p, 80, "E(1/2)")?;

        let x = one_half(p)?;
        let f0 = x.elliptic_f(&zero, p, rm, &mut cc)?;
        bits_agree(&f0, &x.asin(p, rm, &mut cc)?, p, 80, "F(x|0)")?;
        let f1 = x.elliptic_f(&one, p, rm, &mut cc)?;
        bits_agree(&f1, &x.atanh(p, rm, &mut cc)?, p, 80, "F(x|1)")?;
        let ei1 = x.elliptic_e(&one, p, rm, &mut cc)?;
        bits_agree(&ei1, &x, p, 80, "E(x|1)")?;
        assert!(one.elliptic_f(&one, p, rm, &mut cc).is_err());

        let m13 = one.div(&three, p, RoundingMode::None)?;
        let k13 = m13.elliptic_k(p, rm, &mut cc)?;
        let f1 = one.elliptic_f(&m13, p, rm, &mut cc)?;
        bits_agree(&k13, &f1, p, 80, "K=F(1)")?;
        let e13 = m13.elliptic_e_complete(p, rm, &mut cc)?;
        let ei = one.elliptic_e(&m13, p, rm, &mut cc)?;
        bits_agree(&e13, &ei, p, 80, "E=E(1)")?;

        let pim = zero.elliptic_pi_complete(&half, p, rm, &mut cc)?;
        bits_agree(&pim, &k_half, p, 80, "Π(0,m)=K")?;
        let pi0 = zero.elliptic_pi(&x, &half, p, rm, &mut cc)?;
        let fx = x.elliptic_f(&half, p, rm, &mut cc)?;
        bits_agree(&pi0, &fx, p, 80, "Π(0;x|m)=F")?;

        let n = one.div(&four, p, RoundingMode::None)?;
        let got = n.elliptic_pi(&x, &zero, p, rm, &mut cc)?;
        let root = one
            .sub(&n, p, RoundingMode::None)?
            .sqrt(p, RoundingMode::None)?;
        let want = root
            .mul(&x, p, RoundingMode::None)?
            .div(
                &one.sub(&x.mul(&x, p, RoundingMode::None)?, p, RoundingMode::None)?
                    .sqrt(p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )?
            .atan(p, rm, &mut cc)?
            .div(&root, p, RoundingMode::None)?;
        bits_agree(&got, &want, p, 80, "Π(n;x|0)")?;

        let pmm = m13.elliptic_pi_complete(&m13, p, rm, &mut cc)?;
        let want = e13.div(
            &one.sub(&m13, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
        bits_agree(&pmm, &want, p, 40, "Π(m,m)")?;

        assert!(one.elliptic_pi_complete(&half, p, rm, &mut cc).is_err());
        assert!(two.elliptic_f(&zero, p, rm, &mut cc).is_err());
        Ok(())
    }

    // References: /workspace/zf/refs/gen_2f1_refs.py (mpmath 1.4.1, mp.dps = 50, each
    // decimal input first rounded to the 128-bit value parsed here). Covers z in (1/2, 1) with
    // non-integer and integer c-a-b (log cases, both signs), near-integer c-a-b, z near 1,
    // |z| <= 1/2, negative z down to -1e6, terminating series, and z = 1.
    #[test]
    fn test_hypergeom_2f1_mpmath() -> Result<(), crate::Error> {
        let p = 128;
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;
        let dec =
            |s: &str, cc: &mut Consts| ExactNumNumber::parse(s, crate::defs::Radix::Dec, p, rm, cc);
        let cases = [
            (
                "1",
                "1",
                "2",
                "0.7",
                "1.71996114903705141803249459680262643278951811",
            ),
            (
                "0.5",
                "0.5",
                "1.5",
                "0.7",
                "1.1846587084327787287338142454272246065570138",
            ),
            (
                "0.3",
                "0.7",
                "1.9",
                "0.9",
                "1.1786297224210022954609654913371386168786381",
            ),
            (
                "1.5",
                "2.5",
                "2",
                "0.8",
                "22.1469514371671108416364442041104122572588905",
            ),
            (
                "0.25",
                "0.5",
                "0.75",
                "0.99",
                "1.839926250774784062953575844257627101164147",
            ),
            (
                "2.5",
                "-1.5",
                "0.5",
                "0.8",
                "-0.626099033699941114994568627244757345917984746",
            ),
            (
                "0.5",
                "0.5",
                "2.5",
                "0.75",
                "1.10459978807807261686469275254738524409468875",
            ),
            (
                "3",
                "1.5",
                "3",
                "0.7",
                "6.08580619450184570507744203112002371056816759",
            ),
            (
                "5",
                "7",
                "13",
                "0.9",
                "81.5980783011695808602321348648454229291473389",
            ),
            (
                "1",
                "2",
                "5",
                "0.8",
                "1.57356005363371137648526484362443864006988493",
            ),
            (
                "2",
                "3",
                "1",
                "0.6",
                "85.937500000000000000000000000000000001102026",
            ),
            (
                "1.5",
                "2.5",
                "1",
                "0.95",
                "13325.0437834819832041243868639672252078542763",
            ),
            (
                "0.5",
                "0.5",
                "1.000000001",
                "0.9",
                "1.64126441320209173433786497074607887944555805",
            ),
            (
                "0.1",
                "0.2",
                "0.3",
                "0.9999",
                "1.62806497725092035401785531448348646268340797",
            ),
            (
                "-0.75",
                "1.25",
                "0.5",
                "0.6",
                "-0.324842102999361368979248665126531883694813935",
            ),
            (
                "0.5",
                "0.5",
                "1",
                "0.5000001",
                "1.18034065295136224658707771685792578406603286",
            ),
            (
                "1",
                "1",
                "2",
                "0.99",
                "4.6516870565536276444807908175441701162413267",
            ),
            (
                "0.5",
                "0.5",
                "1",
                "0.999999",
                "5.2801571547718662756646199176225633673643885",
            ),
            (
                "0.5",
                "0.5",
                "1",
                "0.999999999999999999999999999999",
                "22.8706103669194821708152195250935907947524458",
            ),
            (
                "1.5",
                "0.5",
                "1.75",
                "0.9999999",
                "117.919774460116271662376787995205289131958635",
            ),
            (
                "0.3",
                "0.4",
                "2.2",
                "0.99999999",
                "1.09006069218348579683600983399311026782044901",
            ),
            (
                "1.5",
                "2.5",
                "3",
                "0.5",
                "2.33449175759488983444516920804357125900741473",
            ),
            (
                "20",
                "20",
                "0.5",
                "0.3",
                "2.04754351497133932937461952694010408655766245e+13",
            ),
            (
                "0.3",
                "0.4",
                "2",
                "1e-30",
                "1.00000000000000000000000000000006",
            ),
            (
                "1",
                "1",
                "2",
                "-0.5",
                "0.810930216216328763956026230928698273143980847",
            ),
            (
                "-2.5",
                "3.25",
                "0.125",
                "-0.5",
                "88.071831107804353538989414947295234489045681",
            ),
            (
                "0.5",
                "1",
                "1.5",
                "-0.9",
                "0.800130255087768810416961486533550204101507955",
            ),
            (
                "1",
                "1",
                "2",
                "-1",
                "0.693147180559945309417232121458176568075500134",
            ),
            (
                "0.5",
                "0.5",
                "1",
                "-1",
                "0.834626841674073186281429732799046808993993013",
            ),
            (
                "0.5",
                "0.25",
                "1.5",
                "-3",
                "0.866870889001128608899121186546323506610332879",
            ),
            (
                "0.5",
                "1.5",
                "2",
                "-3",
                "0.570349449920576636067466296715804977195622434",
            ),
            (
                "1",
                "1",
                "3",
                "-7",
                "0.393287034017905609225043710816172964645387887",
            ),
            (
                "1.25",
                "1.25",
                "1.5",
                "-2",
                "0.314907247447028609342530404747985646146595563",
            ),
            (
                "2",
                "2",
                "3",
                "-100",
                "0.000725004301388053870374859455362795851180156715",
            ),
            (
                "0.5",
                "0.75",
                "1.25",
                "-1000000",
                "0.00214371790435314040885647402372034039486027383",
            ),
            (
                "-3",
                "2",
                "1.5",
                "0.9",
                "-0.04502857142857142857142857142857142857165964",
            ),
            (
                "-10",
                "2.5",
                "3.5",
                "-5",
                "14136956.2280906601811165009442031921504259639",
            ),
            ("-2", "1", "1", "3", "4.0"),
            (
                "2.5",
                "-4",
                "-6",
                "0.95",
                "5.77798943684895833333333333333333333332709832",
            ),
            (
                "1.5",
                "0.5",
                "2.5",
                "1",
                "2.35619449019234492884698253745962716314787705",
            ),
            ("0.5", "-0.25", "0.5", "1", "0"),
        ];
        for (a, b, c, z, want) in cases {
            let av = dec(a, &mut cc)?;
            let bv = dec(b, &mut cc)?;
            let cv = dec(c, &mut cc)?;
            let zv = dec(z, &mut cc)?;
            let label = alloc::format!("2F1({a}, {b}; {c}; {z})");
            let got = av.hypergeom_2f1(&bv, &cv, &zv, p, rm, &mut cc);
            assert!(got.is_ok(), "{label}: {got:?}");
            let got = got?;
            let want = dec(want, &mut cc)?;
            if want.is_zero() {
                assert!(got.is_zero(), "{label}: expected exact 0");
            } else {
                bits_agree(&got, &want, p, 120, &label)?;
            }
        }
        let invalid = [
            ("1", "1", "2", "1.5", "non-terminating z > 1"),
            ("1", "1", "1.5", "1", "z = 1 with c-a-b < 0"),
            ("1", "1", "2", "1", "z = 1 with c-a-b = 0"),
            ("0.5", "0.5", "-2", "0.3", "pole at c = -2"),
        ];
        for (a, b, c, z, why) in invalid {
            let av = dec(a, &mut cc)?;
            let bv = dec(b, &mut cc)?;
            let cv = dec(c, &mut cc)?;
            let zv = dec(z, &mut cc)?;
            assert!(
                av.hypergeom_2f1(&bv, &cv, &zv, p, rm, &mut cc).is_err(),
                "{why}"
            );
        }
        Ok(())
    }

    /// `_2F_1(10000, 1; 10001; 1/2)` needs more than the old 10_000-term cap before
    /// `series_done` (n > max(|a|,|b|,|c|)) can fire. Gold: mpmath `hyp2f1`.
    #[test]
    fn test_hypergeom_2f1_param_10000_beyond_old_term_cap() -> Result<(), crate::Error> {
        let p = 64;
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;
        let dec =
            |s: &str, cc: &mut Consts| ExactNumNumber::parse(s, crate::defs::Radix::Dec, p, rm, cc);
        let a = ExactNumNumber::from_word(10_000, p)?;
        let b = ExactNumNumber::from_word(1, p)?;
        let c = ExactNumNumber::from_word(10_001, p)?;
        let z = one_half(p)?;
        let got = a.hypergeom_2f1(&b, &c, &z, p, rm, &mut cc)?;
        let want = dec(
            "1.999800059974014989189356552302545895116026515887677357177912089351346",
            &mut cc,
        )?;
        bits_agree(&got, &want, p, 56, "2F1(10000, 1; 10001; 1/2)")?;
        Ok(())
    }

    // References: /workspace/zf/refs/gen_2f1_refs.py (mpmath betainc(a, b, 0, x,
    // regularized=True), mp.dps = 50, inputs rounded to 128 bits). I_0.7(30.5, 2.25) used to
    // come back as 2.68e240 through the old 2F1 Pfaff branch.
    #[test]
    fn test_betainc_mpmath() -> Result<(), crate::Error> {
        let p = 128;
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;
        let dec =
            |s: &str, cc: &mut Consts| ExactNumNumber::parse(s, crate::defs::Radix::Dec, p, rm, cc);
        let cases = [
            (
                "30.5",
                "2.25",
                "0.7",
                "0.00030323093547038030663005371621383111406105688",
            ),
            (
                "2.25",
                "30.5",
                "0.3",
                "0.999696769064529619693369946283786168885938943",
            ),
            (
                "2.5",
                "1.5",
                "0.3",
                "0.0889437231706655993544962677753584934056995263",
            ),
            (
                "0.5",
                "0.5",
                "0.001",
                "0.0201350416333774909723443944048709219134992636",
            ),
            (
                "0.5",
                "0.5",
                "1e-10",
                "6.36619772378191672615472207370620472228663812e-6",
            ),
            (
                "0.5",
                "0.5",
                "0.999",
                "0.979864958366622509027655605595129078099969071",
            ),
            (
                "10",
                "30",
                "0.05",
                "1.6241391155028577143348453969451423955897253e-5",
            ),
            (
                "5.5",
                "0.75",
                "0.6",
                "0.0378173332267801158863537056009773691191750344",
            ),
            (
                "100",
                "100",
                "0.45",
                "0.0783879327122204842061564596521365858040127826",
            ),
            (
                "50",
                "50",
                "0.2",
                "1.32861996318603359425556163603788268483320429e-11",
            ),
            (
                "50",
                "50",
                "0.8",
                "0.999999999986713800368139664057444383639621173",
            ),
            (
                "0.01",
                "3.5",
                "0.2",
                "0.996055418033318043865926563058446081518333464",
            ),
            (
                "1000.5",
                "2.5",
                "0.999",
                "0.848660344684032250833214132642228229968976541",
            ),
            (
                "0.1",
                "10",
                "0.9",
                "0.999999999998564632108319054323590977126689725",
            ),
            (
                "3.5",
                "4",
                "0.6",
                "0.766962244469439316871787712000428856055679206",
            ),
            (
                "2",
                "5",
                "0.95",
                "0.999998203124999999999999999999999999999999895",
            ),
            ("0.001", "0.001", "0.5", "0.5"),
        ];
        for (a, b, x, want) in cases {
            let av = dec(a, &mut cc)?;
            let bv = dec(b, &mut cc)?;
            let xv = dec(x, &mut cc)?;
            let label = alloc::format!("I_{x}({a}, {b})");
            let got = av.betainc(&bv, &xv, p, rm, &mut cc);
            assert!(got.is_ok(), "{label}: {got:?}");
            let got = got?;
            bits_agree(&got, &dec(want, &mut cc)?, p, 120, &label)?;
        }
        Ok(())
    }

    // References: /workspace/zf/refs/gen_gamma_refs.py (mpmath, mp.dps = 400). Before 1.0.5
    // the Stirling series started at z >= 32 (gamma) / z >= 128 (digamma) with 64 terms, which
    // capped gamma near 274 correct bits and digamma near 530 whatever the precision asked.
    #[test]
    fn test_gamma_digamma_high_precision() -> Result<(), crate::Error> {
        let p = 1024;
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;
        let dec =
            |s: &str, cc: &mut Consts| ExactNumNumber::parse(s, crate::defs::Radix::Dec, p, rm, cc);
        let cases = [
            ("0.3", "2.99156898768759062831251651590491779111280602492171511274411965095638876787632021799025995255765056542142562682614202237570044313685592565039738332500699651731492217445693233669501544707963377073843625467587119914446338609689724358296069475868570923910448073904676197524186031618193782973959519480201743913441846269484020926939138", "-3.50252422220013298896449450737198159953790828840450209566491975126416371903591035905203493019247296315255849330906099086131413050655055954527794120627288788834403942521473264272871185619845792144602545923064374445895433893233464772352890648958999888003674386339797798702859422961831709565727973413316725322338864624067819849202019"),
            ("3.7", "4.1706517837966031653936029986179837279404455809898292945722466324606422685813692411505269066909944153872906767134758267882411321237339947614951553474226244444095217987337371855360796708629477024837169605083452577126872521433265350148023719876668471470669216309124355964013589712118318138521849706134504731375899736002728295617768", "1.16715353936151138587386396614504688117374878787689916455965074078023789550185395821468166859491773967994604749319569991382969422962925178688049288952435893289539968114159009249272573586281515296928168200307293681328111704651702375912858415352758641207254910912020702174151466809993841444669828620914213481044476105894875417630202"),
            ("70.5", "1429158828303246846573838562651630902882496420796728590630912851373230265522054175373279959655886352.17790297183463623799296473556344225544418545904118072813016602611874981916619673887717371489733555185792563204257087407395041763541149993405568290470895024805331342693136816762678344437425087011896968236201817908802084497213971124", "4.248503745147059295990802672685553304359490704691494971797198001646067548177592356909609322645124920315927005996751096556357511162519211891173394099381260107244573276858072927243357078818288905238442575941336918427362361945654183685846907098161387913639662258560146323261645527801019331192688800792418372732602795783524642032161"),
            ("-2.5", "-0.945308720482941881225689324448610764158693043265273135047364154588219351781883830066640350260557154888654305932950704355325948195053025388013383311788153103179753439281656065519960092282107107322097117535609580160093140808187455850942240101432124809912238386397286741500758280280483315416223644000483412484060657854653415164622928", "1.10315664064524318722569033366791109947350706200623255961953941279501169594961256451799294938208254206803225737571821271833512218066386271637715116575159120655171119126698645854837864662628806147732505042720246680802275486308815174019100773441769687963010702067494673598401202657231149432143243715766674712012979544692612931066196"),
        ];
        for (x, g, psi) in cases {
            let xv = dec(x, &mut cc)?;
            let got_g = xv.gamma(p, rm, &mut cc)?;
            bits_agree(
                &got_g,
                &dec(g, &mut cc)?,
                p,
                1016,
                &alloc::format!("gamma({x})"),
            )?;
            let got_psi = xv.digamma(p, rm, &mut cc)?;
            bits_agree(
                &got_psi,
                &dec(psi, &mut cc)?,
                p,
                1016,
                &alloc::format!("digamma({x})"),
            )?;
        }
        Ok(())
    }

    #[test]
    fn test_legendre_hypergeom_betainc() -> Result<(), crate::Error> {
        let p = 256;
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;
        let one = ExactNumNumber::from_word(1, p)?;
        let zero = ExactNumNumber::new(p)?;
        let two = ExactNumNumber::from_word(2, p)?;
        let three = ExactNumNumber::from_word(3, p)?;
        let half = one_half(p)?;

        for n in 0u32..=8 {
            let p1 = one.legendre_p(n, p, rm)?;
            bits_agree(&p1, &one, p, 80, "P_n(1)")?;
            let pm = one.neg()?.legendre_p(n, p, rm)?;
            let want = if n % 2 == 0 { ExactNumNumber::from_word(1, p)? } else { one.neg()? };
            bits_agree(&pm, &want, p, 80, "P_n(-1)")?;
        }
        let p2 = zero.legendre_p(2, p, rm)?;
        bits_agree(&p2, &half.neg()?, p, 80, "P_2(0)")?;
        let p100 = one.legendre_p(100, p, rm)?;
        bits_agree(&p100, &one, p, 80, "P_100(1)")?;
        let p128 = one.neg()?.legendre_p(128, p, rm)?;
        bits_agree(&p128, &one, p, 80, "P_128(-1)")?;
        let xh = one_half(p)?;
        let pn = xh.legendre_p(128, p, rm)?;
        let pn1 = xh.legendre_p(127, p, rm)?;
        let pn2 = xh.legendre_p(126, p, rm)?;
        // (n+1) P_{n+1} = (2n+1) x P_n − n P_{n−1}  at n=127
        let n = ExactNumNumber::from_word(127, p)?;
        let np1 = ExactNumNumber::from_word(128, p)?;
        let two_n1 = ExactNumNumber::from_word(255, p)?;
        let lhs = np1.mul(&pn, p, RoundingMode::None)?;
        let rhs = two_n1
            .mul(&xh, p, RoundingMode::None)?
            .mul(&pn1, p, RoundingMode::None)?
            .sub(&n.mul(&pn2, p, RoundingMode::None)?, p, RoundingMode::None)?;
        bits_agree(&lhs, &rhs, p, 20, "Bonnet n=127")?;

        let p11 = half.assoc_legendre_p(1, 1, p, rm)?;
        let want = one
            .sub(
                &half.mul(&half, p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )?
            .sqrt(p, RoundingMode::None)?
            .neg()?;
        bits_agree(&p11, &want, p, 80, "P_1^1(1/2)")?;

        let f0 = half.hypergeom_2f1(&half, &one, &zero, p, rm, &mut cc)?;
        bits_agree(&f0, &one, p, 80, "2F1(...,0)")?;

        let nm2 = two.neg()?;
        let fpoly = nm2.hypergeom_2f1(&one, &one, &half, p, rm, &mut cc)?;
        bits_agree(
            &fpoly,
            &one.div(&four_word(p)?, p, RoundingMode::None)?,
            p,
            80,
            "2F1(-2,1;1;1/2)",
        )?;

        let fln = one.hypergeom_2f1(&one, &two, &half, p, rm, &mut cc)?;
        let ln2 = cc.ln_2_num(p, RoundingMode::None)?;
        bits_agree(
            &fln,
            &two.mul(&ln2, p, RoundingMode::None)?,
            p,
            80,
            "2F1(1,1;2;1/2)",
        )?;

        let fg = nm2.hypergeom_2f1(&one, &three, &one, p, rm, &mut cc)?;
        bits_agree(&fg, &half, p, 80, "2F1(-2,1;3;1)")?;

        let fm1 = one.hypergeom_2f1(&one, &two, &one.neg()?, p, rm, &mut cc)?;
        bits_agree(&fm1, &ln2, p, 40, "2F1(1,1;2;-1)=ln2")?;
        let fterm2 = nm2.hypergeom_2f1(&one, &three, &two, p, rm, &mut cc)?;
        let third = one.div(&three, p, RoundingMode::None)?;
        bits_agree(&fterm2, &third, p, 80, "2F1(-2,1;3;2)=1/3")?;

        assert!(one
            .hypergeom_2f1(&one, &three, &two, p, rm, &mut cc)
            .is_err());

        let fk = half.hypergeom_2f1(&half, &one, &half, p, rm, &mut cc)?;
        let k = half.elliptic_k(p, rm, &mut cc)?;
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let want = two
            .div(&pi, p, RoundingMode::None)?
            .mul(&k, p, RoundingMode::None)?;
        bits_agree(&fk, &want, p, 40, "2F1 = 2K/π")?;

        assert!(half.betainc(&half, &zero, p, rm, &mut cc)?.is_zero());
        let i1 = half.betainc(&half, &one, p, rm, &mut cc)?;
        bits_agree(&i1, &one, p, 80, "I_1(a,b)")?;
        // I_{1/2}(1,1) = 1/2
        let i11 = one.betainc(&one, &half, p, rm, &mut cc)?;
        bits_agree(&i11, &half, p, 80, "I_{1/2}(1,1)")?;
        assert!(zero.betainc(&one, &half, p, rm, &mut cc).is_err());
        Ok(())
    }

    fn four_word(p: usize) -> Result<ExactNumNumber, crate::Error> {
        ExactNumNumber::from_word(4, p)
    }

    /// More precisions than `BERNOULLI_P_SLOTS` requested from many threads: a slot grown
    /// under the write lock must still be there when the caller reads it. Before the fix
    /// another thread could evict it between the write and the read, and the call returned
    /// `InvalidArgument` for a valid request (seen as a NaN Cl3 in tier1 under load).
    #[cfg(feature = "std")]
    #[test]
    fn bernoulli_cache_slot_survives_concurrent_eviction() {
        let handles: Vec<_> = (0..16usize)
            .map(|t| {
                std::thread::spawn(move || {
                    let mut fails = 0usize;
                    for i in 0..400usize {
                        let p = WORD_BIT_SIZE * (1 + (t + i) % (3 * BERNOULLI_P_SLOTS));
                        if with_even_bernoulli(24, p, |v| Ok(v.len())) != Ok(24) {
                            fails += 1;
                        }
                    }
                    fails
                })
            })
            .collect();
        let fails: usize = handles
            .into_iter()
            .map(|h| h.join().unwrap_or(usize::MAX))
            .sum();
        assert_eq!(fails, 0, "with_even_bernoulli lost its slot {fails} times");
    }

    fn n_at(
        x: &ExactNumNumber,
        a: &ExactNumNumber,
        q: usize,
        cc: &mut Consts,
    ) -> Result<i128, crate::Error> {
        let ln2 = cc.ln_2_num(q, RoundingMode::None)?;
        let la = x.abs()?.ln(q, RoundingMode::None, cc)?;
        let al = a.mul(&la, q, RoundingMode::None)?;
        let one = ExactNumNumber::from_word(1, q)?;
        let ap1 = a.add(&one, q, RoundingMode::None)?;
        let lg = ap1.ln_gamma(q, RoundingMode::None, cc)?;
        let lam = al
            .sub(&lg, q, RoundingMode::None)?
            .div(&ln2, q, RoundingMode::None)?
            .sub(a, q, RoundingMode::None)?;
        if lam.exponent() > 62 {
            return Ok(if lam.is_negative() { -(1i128 << 100) } else { 1i128 << 100 });
        }
        let fl = lam.floor_ref()?;
        let nf = if fl.is_zero() { 0 } else { fl.abs()?.int_as_usize()? as i128 };
        Ok(if fl.is_negative() { -nf } else { nf })
    }

    fn check_skip_is_conservative(
        x: &ExactNumNumber,
        nu: &ExactNumNumber,
        p: usize,
        cc: &mut Consts,
        label: &str,
    ) -> Result<(), crate::Error> {
        let a = nu;
        let n_est = super::estimate_log2_leading_ji(x, a, cc)?;
        let q = (16 * WORD_BIT_SIZE)
            .max(a.mantissa_max_bit_len() + 4 * WORD_BIT_SIZE)
            .max(x.mantissa_max_bit_len() + 4 * WORD_BIT_SIZE);
        let n_true = n_at(x, a, q, cc)?;
        let lo = EXPONENT_MIN as i128 + 2 * p as i128 + 256;
        let hi = EXPONENT_MAX as i128 - 256;
        let margin = 1i128 << 20;
        if let Some(n_est) = n_est {
            if n_est > lo + margin && n_est < hi - margin {
                assert!(
                    (lo..=hi).contains(&n_true),
                    "{label}: estimate skipped (n_est={n_est}) but precise n={n_true} is outside [{lo}, {hi}]"
                );
                let err = (n_est - n_true).abs();
                assert!(
                    err <= 4,
                    "{label}: |n_est-n|={err} (n_est={n_est}, n={n_true}); bound is 4"
                );
            }
        }
        Ok(())
    }

    /// The 4-word skip must never fire when the true leading exponent is outside
    /// the padded normal window, including ν that 256-bit rounding maps onto a pole.
    #[test]
    fn zf_37_ji_scaled_estimate_never_skips_near_edges() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let p64 = 64;
        let x30 = pow2(-30, p64)?;
        let one = ExactNumNumber::from_word(1, p64)?;
        let two = ExactNumNumber::from_word(2, p64)?;
        let three = ExactNumNumber::from_word(3, p64)?;
        let half = one.div(&two, p64, RoundingMode::ToEven)?;
        let third = one.div(&three, p64, RoundingMode::ToEven)?;
        let five_two = ExactNumNumber::from_word(5, p64)?.div(&two, p64, RoundingMode::ToEven)?;

        // ν = -1 + 2^{-300}: 256-bit rounding of ν+1 is 0, a pole of Γ.
        let eps = pow2(-300, 384)?;
        let mut nu_pole = ExactNumNumber::from_word(1, p64)?;
        nu_pole.set_precision(384, RoundingMode::None)?;
        nu_pole = eps.sub(&nu_pole, 384, RoundingMode::None)?; // 2^{-300} - 1

        let third_wide = ExactNumNumber::from_word(1, 8192)?.div(
            &ExactNumNumber::from_word(3, 8192)?,
            8192,
            RoundingMode::ToEven,
        )?;

        let neg_half = half.neg()?;
        let (x_m4, x_m40) = (pow2(-4, p64)?, pow2(-40, p64)?);
        let pos_201 = ExactNumNumber::from_word(201, p64)?.div(&two, p64, RoundingMode::ToEven)?;
        let neg_201 = pos_201.neg()?;
        let cases: [(&str, &ExactNumNumber, &ExactNumNumber, usize); 11] = [
            ("J-shape 1/2 p64", &x30, &half, p64),
            ("-1/2 p64", &x30, &neg_half, p64),
            ("5/2 p64", &x30, &five_two, p64),
            ("1/3 p64", &x30, &third, p64),
            ("1/3 8192-bit p64", &x30, &third_wide, p64),
            ("1/3 8192-bit p8192", &x30, &third_wide, 8192),
            ("pole-adjacent p64", &x30, &nu_pole, p64),
            ("pole-adjacent p2048", &x30, &nu_pole, 2048),
            ("2^{-4} 1/2", &x_m4, &half, p64),
            ("2^{-40} 201/2", &x_m40, &pos_201, p64),
            ("2^{-40} -201/2", &x_m40, &neg_201, p64),
        ];
        for (label, x, nu, p) in cases {
            check_skip_is_conservative(x, nu, p, &mut cc, label)?;
        }

        // The pole-adjacent order must not leak ExponentOverflow / DivisionByZero
        // from the 4-word estimate (that would be a NaN at the ExactNum wrapper).
        let j = x30.bessel_j_nu(&nu_pole, p64, RoundingMode::ToEven, &mut cc)?;
        assert!(
            !j.is_zero(),
            "J of -1 plus 2^-300 at 2^-30 should be finite nonzero"
        );
        let i = x30.bessel_i(&nu_pole, p64, RoundingMode::ToEven, &mut cc)?;
        assert!(
            !i.is_zero(),
            "I of -1 plus 2^-300 at 2^-30 should be finite nonzero"
        );
        // Estimate itself must refuse (ν+1 rounds to 0 at 4 words).
        assert_eq!(
            super::estimate_log2_leading_ji(&x30, &nu_pole, &mut cc)?,
            None,
            "estimate should refuse a pole after 256-bit rounding"
        );
        Ok(())
    }
}
