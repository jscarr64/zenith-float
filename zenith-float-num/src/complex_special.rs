//! Complex `erf` / `erfc` (Faddeeva) and `gamma` / `ln_gamma` / `digamma`.
//!
//! Software-limb [`ExactComplex`] arithmetic only. The error-function path is the
//! Faddeeva function throughout — not the real series evaluated at \(\lvert z\rvert\).

use crate::common::util::bump_prec_retry;
use crate::common::util::round_p;
use crate::defs::Exponent;
use crate::Consts;
use crate::Error;
use crate::ExactComplex;
use crate::ExactNum;
use crate::RoundingMode;
use crate::WORD_BIT_SIZE;
use alloc::vec::Vec;
use core::ops::Neg;

/// L1 radius \(\lvert\mathrm{Re}\,z\rvert+\lvert\mathrm{Im}\,z\rvert\) below which
/// the Faddeeva power series is used (Poppe–Wijers inner region).
const FADDEEVA_SERIES_L1: u32 = 8;

/// Bernoulli terms kept in the Stirling series for \(\ln\Gamma\) and \(\psi\).
const GAMMA_STIRLING_TERMS: usize = 64;

/// Minimum shift target: raise \(z\) by integers until \(\mathrm{exponent}(\lvert z\rvert)\)
/// is at least this (\(\lvert z\rvert\ge 2^{k-1}\)) so Stirling terms decay. The actual target
/// grows with the precision (see [`gamma_shift_exp`]): with a fixed target of 6 the 64-term series
/// capped \(\Gamma\) / \(\ln\Gamma\) at about 275 correct bits whatever `p` was.
const GAMMA_STIRLING_MIN_ABS_EXP: Exponent = 6;

/// Same minimum shift target for digamma (real kernel uses exponent \(< 8\)).
const DIGAMMA_STIRLING_MIN_ABS_EXP: Exponent = 8;

/// Shift target so the last kept Stirling term, \(\approx 2^{364}/\lvert z\rvert^{127}\), is below
/// \(2^{-(p+8)}\) (same bound as the real kernel).
fn gamma_shift_exp(p: usize) -> Exponent {
    crate::ops::special::stirling_min_exponent(p, 372, 127, GAMMA_STIRLING_MIN_ABS_EXP)
}

/// Digamma analogue of [`gamma_shift_exp`].
fn digamma_shift_exp(p: usize) -> Exponent {
    crate::ops::special::stirling_min_exponent(p, 379, 128, DIGAMMA_STIRLING_MIN_ABS_EXP)
}

/// Positive integers \(n\) for which \(\Gamma(n)=(n-1)!\) is evaluated by multiplying
/// \(1\ldots n-1\) instead of Stirling.
const GAMMA_FACTORIAL_MAX: u32 = 64;

pub(crate) fn nan_pair(e: Error) -> ExactComplex {
    ExactComplex::new(ExactNum::nan(Some(e)), ExactNum::nan(Some(e)))
}

pub(crate) fn neg_c(z: ExactComplex) -> ExactComplex {
    z.neg()
}

pub(crate) fn two_c(p: usize) -> ExactComplex {
    ExactComplex::from_real(ExactNum::from_u8(2, p), p)
}

pub(crate) fn half_c(p: usize) -> ExactComplex {
    let h = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, RoundingMode::None);
    ExactComplex::from_real(h, p)
}

pub(crate) fn pi_c(p: usize, cc: &mut Consts) -> ExactComplex {
    ExactComplex::from_real(cc.pi(p, RoundingMode::None), p)
}

/// Guard bits for [`ExactComplex::erf_power_series`]: `1.5·2^e + 16` with \(\lvert u\rvert^2<2^e\).
fn erf_series_guard_bits(u: &ExactComplex) -> usize {
    let a2 = u.re().sqr(64, RoundingMode::None).add(
        &u.im().sqr(64, RoundingMode::None),
        64,
        RoundingMode::None,
    );
    match a2.exponent() {
        Some(e) if e > 0 && !a2.is_zero() => (3usize << (e as u32).min(40)) / 2 + 16,
        _ => 16,
    }
}

pub(crate) fn term_negligible(t: &ExactComplex, p: usize) -> bool {
    let m = t.abs(p, RoundingMode::None);
    m.is_zero() || m.exponent().is_some_and(|e| (e as i64) + (p as i64) < 0)
}

fn l1_below_series_bound(z: &ExactComplex, p: usize) -> bool {
    let (re, im) = (z.re(), z.im());
    let s = match (re.is_negative(), im.is_negative()) {
        (false, false) => ExactNum::add_ref(re, im, p, RoundingMode::None),
        (false, true) => ExactNum::sub_ref(re, im, p, RoundingMode::None),
        (true, false) => ExactNum::sub_ref(im, re, p, RoundingMode::None),
        (true, true) => re.neg_via_mul().sub(im, p, RoundingMode::None),
    };
    let bound = ExactNum::from_u32(FADDEEVA_SERIES_L1, p);
    matches!(s.compare(&bound), Some(c) if c < 0)
}

/// Power series when inside the Poppe–Wijers L1 ball, or when \(\lvert z\rvert^2 < p\)
/// so a truncated asymptotic cannot meet working precision (optimal cut is \(m\sim\lvert z\rvert^2\)).
fn use_faddeeva_series(z: &ExactComplex, p: usize) -> bool {
    if l1_below_series_bound(z, p) {
        return true;
    }
    let az = z.abs(p, RoundingMode::None);
    let az2 = az.sqr(p, RoundingMode::None);
    let thresh = ExactNum::from_u32(p.min(u32::MAX as usize) as u32, p);
    matches!(az2.compare(&thresh), Some(c) if c < 0)
}

fn re_positive(z: &ExactComplex) -> bool {
    z.re().is_positive() && !z.re().is_zero()
}

fn im_negative_or_neg_real(z: &ExactComplex) -> bool {
    z.im().is_negative() || (z.im().is_zero() && z.re().is_negative())
}

pub(crate) fn is_nonpos_integer(z: &ExactComplex) -> bool {
    z.im().is_zero() && z.re().is_int() && (z.re().is_zero() || z.re().is_negative())
}

fn abs_needs_shift(z: &ExactComplex, p: usize, min_exp: Exponent) -> bool {
    match z.abs(p, RoundingMode::None).exponent() {
        Some(e) => e < min_exp,
        None => false,
    }
}

pub(crate) fn series_term_cap(p: usize) -> usize {
    p.saturating_add(WORD_BIT_SIZE)
}

pub(crate) fn ziv_complex<F>(p: usize, rm: RoundingMode, mut compute: F) -> ExactComplex
where
    F: FnMut(usize) -> ExactComplex,
{
    let p = round_p(p);
    let mut p_inc = WORD_BIT_SIZE;
    let Some(mut p_wrk) = p.checked_add(p_inc) else {
        return nan_pair(Error::InvalidArgument);
    };
    p_wrk = round_p(p_wrk);
    loop {
        let p_x = match p_wrk.checked_add(WORD_BIT_SIZE.saturating_mul(2)) {
            Some(v) => v,
            None => return nan_pair(Error::InvalidArgument),
        };
        let z = compute(p_x);
        let mut re = z.re().clone();
        let mut im = z.im().clone();
        let ok_re = re.try_set_precision(p, rm, p_wrk);
        let ok_im = im.try_set_precision(p, rm, p_wrk);
        if ok_re && ok_im {
            return ExactComplex::new(re, im);
        }
        if bump_prec_retry(&mut p_wrk, &mut p_inc, p).is_err() {
            return nan_pair(Error::PrecisionRetryExhausted);
        }
    }
}

/// Even Bernoulli numbers \(B_2,\ldots,B_{2k_{\max}}\) (exact tangent numbers, one rounding each;
/// see `ops::special::even_bernoulli`).
fn even_bernoulli_numbers(kmax: usize, p: usize) -> Vec<ExactNum> {
    match crate::ops::special::even_bernoulli(kmax, p) {
        Ok(v) => v.into_iter().map(ExactNum::from).collect(),
        Err(e) => (0..kmax).map(|_| ExactNum::nan(Some(e))).collect(),
    }
}

fn factorial_um1(n: u32, p: usize) -> ExactComplex {
    let mut acc = ExactNum::from_u8(1, p);
    if n >= 2 {
        for k in 2..n {
            acc = acc.mul(&ExactNum::from_u32(k, p), p, RoundingMode::None);
        }
    }
    ExactComplex::from_real(acc, p)
}

impl ExactComplex {
    /// Error function \(\mathrm{erf}(z)=1-\mathrm{erfc}(z)\). Entire; NaN in → NaN out.
    ///
    /// # Precision
    ///
    /// - Algorithm: Faddeeva `w(z)` series for `|z|` below `FADDEEVA_SERIES_L1 = 8`; continued fraction otherwise.
    /// - Bound: Ziv on each part (`MAX_PREC_RETRY`).
    /// - MPFR oracle: real axis vs `mpfr_erf`. GNU MPC has no `mpc_erf`. Off-axis: `erf` odd, `erfc=1-erf`.
    ///
    /// Accuracy: correctly-rounded
    pub fn erf(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.erf_word(p, rm, cc))
    }

    /// Implementation of [`Self::erf`] on operands without inexact history.
    fn erf_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return ExactComplex::new(self.re().clone(), self.im().clone());
        }
        let dest = round_p(p);
        ziv_complex(dest, rm, |pw| self.erf_at(pw, dest, cc))
    }

    /// Complementary error function via Faddeeva: \(\mathrm{erfc}(z)=e^{-z^2}w(iz)\).
    /// Entire; NaN in → NaN out.
    ///
    /// # Precision
    ///
    /// - Algorithm: same Faddeeva path as [`Self::erf`].
    /// - Bound: Ziv on each part (`MAX_PREC_RETRY`).
    /// - MPFR oracle: real axis vs `mpfr_erfc` (via `1-erf` identity). GNU MPC has no `mpc_erfc`.
    ///
    /// Accuracy: correctly-rounded
    pub fn erfc(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.erfc_word(p, rm, cc))
    }

    /// Implementation of [`Self::erfc`] on operands without inexact history.
    fn erfc_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return ExactComplex::new(self.re().clone(), self.im().clone());
        }
        let dest = round_p(p);
        ziv_complex(dest, rm, |pw| self.erfc_at(pw, dest, cc))
    }

    /// Gamma \(\Gamma(z)\). Poles at non-positive integers → NaN.
    ///
    /// # Precision
    ///
    /// - Algorithm: Stirling (`GAMMA_STIRLING_TERMS = 64`) plus reflection; factorial for small integers (`GAMMA_FACTORIAL_MAX = 64`).
    /// - Bound: Ziv on each part (`MAX_PREC_RETRY`).
    /// - MPFR oracle: real axis vs `mpfr_gamma`. GNU MPC has no `mpc_gamma`. Integers: `Γ(n)=(n-1)!`.
    ///
    /// Accuracy: correctly-rounded
    pub fn gamma(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.gamma_word(p, rm, cc))
    }

    /// Implementation of [`Self::gamma`] on operands without inexact history.
    fn gamma_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return ExactComplex::new(self.re().clone(), self.im().clone());
        }
        if is_nonpos_integer(self) {
            return nan_pair(Error::InvalidArgument);
        }
        ziv_complex(round_p(p), rm, |pw| self.gamma_at(pw, cc))
    }

    /// Principal \(\ln(\Gamma(z))\): the principal logarithm of \(\Gamma(z)\), imaginary part in
    /// \((-\pi,\pi]\). (This differs from the analytic `loggamma` of mpmath by \(2\pi i k\).)
    /// Poles → NaN.
    ///
    /// # Precision
    ///
    /// - Algorithm: Stirling (`GAMMA_STIRLING_TERMS = 64`) plus reflection.
    /// - Bound: Ziv on each part (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn ln_gamma(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.ln_gamma_word(p, rm, cc))
    }

    /// Implementation of [`Self::ln_gamma`] on operands without inexact history.
    fn ln_gamma_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return ExactComplex::new(self.re().clone(), self.im().clone());
        }
        if is_nonpos_integer(self) {
            return nan_pair(Error::InvalidArgument);
        }
        ziv_complex(round_p(p), rm, |pw| {
            // The shifted Stirling sum yields the analytic continuation (mpmath `loggamma`), whose
            // imaginary part can leave (−π, π]; reduce it to the documented principal
            // ln(Γ(z)): subtract 2π·floor((Im + π)/(2π)).
            let v = self.ln_gamma_at(pw, cc);
            let im = v.im();
            if im.is_nan() {
                return v;
            }
            let pi = cc.pi(pw, RoundingMode::None);
            let two_pi = cc
                .pi(pw, RoundingMode::None)
                .ldexp(1, pw, RoundingMode::None);
            let k = pi
                .add(im, pw, RoundingMode::None)
                .div(&two_pi, pw, RoundingMode::None)
                .floor();
            if k.is_zero() {
                return v;
            }
            let (re, v_im) = v.into_parts();
            let im =
                k.mul(&two_pi, pw, RoundingMode::None)
                    .neg()
                    .add(&v_im, pw, RoundingMode::None);
            ExactComplex::new(re, im)
        })
    }

    /// Digamma \(\psi(z)=\Gamma'/\Gamma\). Poles at non-positive integers → NaN.
    ///
    /// # Precision
    ///
    /// - Algorithm: recurrence plus Bernoulli; reflection for \(\operatorname{Re} z < 0\).
    /// - Bound: Ziv on each part (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn digamma(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.digamma_word(p, rm, cc))
    }

    /// Implementation of [`Self::digamma`] on operands without inexact history.
    fn digamma_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return ExactComplex::new(self.re().clone(), self.im().clone());
        }
        if is_nonpos_integer(self) {
            return nan_pair(Error::InvalidArgument);
        }
        ziv_complex(round_p(p), rm, |pw| self.digamma_at(pw, cc))
    }

    fn erf_at(&self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        if use_faddeeva_series(self, dest_p) {
            Self::erf_power_series(self, work_p, cc)
        } else {
            ExactComplex::one(work_p).sub(
                &self.erfc_via_faddeeva(work_p, dest_p, cc),
                work_p,
                RoundingMode::None,
            )
        }
    }

    fn erfc_at(&self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        if use_faddeeva_series(self, dest_p) {
            ExactComplex::one(work_p).sub(
                &Self::erf_power_series(self, work_p, cc),
                work_p,
                RoundingMode::None,
            )
        } else {
            self.erfc_via_faddeeva(work_p, dest_p, cc)
        }
    }

    /// `erfc(z) = e^{-z²} w(iz)`, with the reflection form when `Re(z) ≤ 0`.
    fn erfc_via_faddeeva(&self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        let z2 = self.mul(self, work_p, RoundingMode::None);
        let em = neg_c(z2).exp(work_p, RoundingMode::None, cc);
        let iz = ExactComplex::i(work_p).mul(self, work_p, RoundingMode::None);
        if re_positive(self) {
            em.mul(
                &Self::faddeeva(&iz, work_p, dest_p, cc),
                work_p,
                RoundingMode::None,
            )
        } else {
            two_c(work_p).sub(
                &em.mul(
                    &Self::faddeeva(&neg_c(iz), work_p, dest_p, cc),
                    work_p,
                    RoundingMode::None,
                ),
                work_p,
                RoundingMode::None,
            )
        }
    }

    /// Faddeeva \(w(z)=e^{-z^2}\mathrm{erfc}(-iz)\). Reduces to \(\mathrm{Im}\,z\ge 0\).
    fn faddeeva(z: &Self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        if im_negative_or_neg_real(z) {
            let z2 = z.mul(z, work_p, RoundingMode::None);
            let two_exp = two_c(work_p).mul(
                &neg_c(z2).exp(work_p, RoundingMode::None, cc),
                work_p,
                RoundingMode::None,
            );
            return two_exp.sub(
                &Self::faddeeva_upper(&z.neg_via_mul(), work_p, dest_p, cc),
                work_p,
                RoundingMode::None,
            );
        }
        Self::faddeeva_upper(z, work_p, dest_p, cc)
    }

    fn faddeeva_upper(z: &Self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        if use_faddeeva_series(z, dest_p) {
            Self::faddeeva_series(z, work_p, cc)
        } else {
            Self::faddeeva_asymp(z, work_p, cc)
        }
    }

    /// Poppe–Wijers power series: \(w(z)=e^{-z^2}(1-\mathrm{erf}(-iz))\).
    fn faddeeva_series(z: &Self, p: usize, cc: &mut Consts) -> Self {
        let u = neg_c(ExactComplex::i(p)).mul(z, p, RoundingMode::None);
        let erf_u = Self::erf_power_series(&u, p, cc);
        let erfc_u = ExactComplex::one(p).sub(&erf_u, p, RoundingMode::None);
        let z2 = z.mul(z, p, RoundingMode::None);
        neg_c(z2)
            .exp(p, RoundingMode::None, cc)
            .mul(&erfc_u, p, RoundingMode::None)
    }

    /// Entire power series \(\mathrm{erf}(u)=\frac{2}{\sqrt\pi}\sum(-1)^n u^{2n+1}/(n!(2n+1))\).
    ///
    /// Terms reach \(\approx e^{\lvert u\rvert^2}\) while the sum can be \(O(1)\), so the series
    /// runs with `1.5·bound(|u|²) + 16` guard bits (\(\log_2 e<1.5\)).
    fn erf_power_series(u: &Self, dest_p: usize, cc: &mut Consts) -> Self {
        let p = dest_p + erf_series_guard_bits(u);
        let pi = cc.pi(p, RoundingMode::None);
        let sqrt_pi = pi.sqrt(p, RoundingMode::None);
        let scale = ExactNum::from_u8(2, p).div(&sqrt_pi, p, RoundingMode::None);
        let u2 = u.mul(u, p, RoundingMode::None);
        let mut upow = u.clone();
        let mut nfact = ExactNum::from_u8(1, p);
        let mut sum = u.clone();
        for n in 1..=series_term_cap(p) {
            nfact = nfact.mul(&ExactNum::from_u32(n as u32, p), p, RoundingMode::None);
            upow = upow.mul(&u2, p, RoundingMode::None);
            let two_n_1 = ExactNum::from_u32((2 * n + 1) as u32, p);
            let den = two_n_1.mul(&nfact, p, RoundingMode::None);
            let mut t = upow.div(&ExactComplex::from_real(den, p), p, RoundingMode::None);
            if n % 2 == 1 {
                t = neg_c(t);
            }
            sum = sum.add(&t, p, RoundingMode::None);
            if term_negligible(&t, p) {
                break;
            }
        }
        ExactComplex::from_real(scale, p).mul(&sum, dest_p, RoundingMode::None)
    }

    /// \(w(z)\sim i/(z\sqrt\pi)\sum(2m-1)!!/(2z^2)^m\) for \(\mathrm{Im}\,z\ge 0\).
    fn faddeeva_asymp(z: &Self, p: usize, cc: &mut Consts) -> Self {
        let z2 = z.mul(z, p, RoundingMode::None);
        let two_z2 = two_c(p).mul(&z2, p, RoundingMode::None);
        let mut term = ExactComplex::one(p);
        let mut s = ExactComplex::one(p);
        let mut prev_e = Exponent::MIN;
        for m in 1..=series_term_cap(p) {
            let odd = ExactComplex::from_real(ExactNum::from_u32((2 * m - 1) as u32, p), p);
            term = term
                .mul(&odd, p, RoundingMode::None)
                .div(&two_z2, p, RoundingMode::None);
            s = s.add(&term, p, RoundingMode::None);
            let e = term
                .abs(p, RoundingMode::None)
                .exponent()
                .unwrap_or(Exponent::MIN);
            if term_negligible(&term, p) {
                break;
            }
            if m > 1 && e > prev_e {
                break;
            }
            prev_e = e;
        }
        let sqrt_pi =
            ExactComplex::from_real(cc.pi(p, RoundingMode::None).sqrt(p, RoundingMode::None), p);
        let den = z.mul(&sqrt_pi, p, RoundingMode::None);
        ExactComplex::i(p)
            .div(&den, p, RoundingMode::None)
            .mul(&s, p, RoundingMode::None)
    }

    pub(crate) fn gamma_at(&self, p: usize, cc: &mut Consts) -> Self {
        if let Some(n) = self.small_pos_int(p) {
            return factorial_um1(n, p);
        }
        if !re_positive(self) && !self.re().is_zero() {
            let pi = pi_c(p, cc);
            let piz = pi.mul(self, p, RoundingMode::None);
            let s = piz.sin(p, RoundingMode::None, cc);
            let omz = ExactComplex::one(p).sub(self, p, RoundingMode::None);
            let g = omz.gamma_positive(p, cc);
            return pi.div(&s.mul(&g, p, RoundingMode::None), p, RoundingMode::None);
        }
        self.gamma_positive(p, cc)
    }

    fn small_pos_int(&self, p: usize) -> Option<u32> {
        if !self.im().is_zero() || !self.re().is_int() || !self.re().is_positive() {
            return None;
        }
        for n in 1u32..=GAMMA_FACTORIAL_MAX {
            let w = ExactNum::from_u32(n, p);
            if self.re().compare(&w) == Some(0) {
                return Some(n);
            }
        }
        None
    }

    fn gamma_positive(&self, p: usize, cc: &mut Consts) -> Self {
        let one = ExactComplex::one(p);
        let mut z = self.clone();
        let mut acc = ExactComplex::one(p);
        while abs_needs_shift(&z, p, gamma_shift_exp(p)) {
            acc = acc.mul(&z, p, RoundingMode::None);
            z = z.add(&one, p, RoundingMode::None);
        }
        let lg = z.ln_gamma_stirling(p, cc);
        let g = lg.exp(p, RoundingMode::None, cc);
        g.div(&acc, p, RoundingMode::None)
    }

    fn ln_gamma_at(&self, p: usize, cc: &mut Consts) -> Self {
        if let Some(n) = self.small_pos_int(p) {
            return factorial_um1(n, p).ln(p, RoundingMode::None, cc);
        }
        if !re_positive(self) && !self.re().is_zero() {
            let piz = pi_c(p, cc).mul(self, p, RoundingMode::None);
            let ln_sin = piz
                .sin(p, RoundingMode::None, cc)
                .ln(p, RoundingMode::None, cc);
            let omz = ExactComplex::one(p).sub(self, p, RoundingMode::None);
            return pi_c(p, cc)
                .ln(p, RoundingMode::None, cc)
                .sub(&ln_sin, p, RoundingMode::None)
                .sub(&omz.ln_gamma_positive(p, cc), p, RoundingMode::None);
        }
        self.ln_gamma_positive(p, cc)
    }

    fn ln_gamma_positive(&self, p: usize, cc: &mut Consts) -> Self {
        let one = ExactComplex::one(p);
        let mut z = self.clone();
        let mut ln_acc = ExactComplex::zero(p);
        while abs_needs_shift(&z, p, gamma_shift_exp(p)) {
            ln_acc = ln_acc.add(&z.ln(p, RoundingMode::None, cc), p, RoundingMode::None);
            z = z.add(&one, p, RoundingMode::None);
        }
        z.ln_gamma_stirling(p, cc)
            .sub(&ln_acc, p, RoundingMode::None)
    }

    fn ln_gamma_stirling(&self, p: usize, cc: &mut Consts) -> Self {
        let ln_z = self.ln(p, RoundingMode::None, cc);
        let zmh = half_c(p).neg().add(self, p, RoundingMode::None);
        let mut s = zmh.mul(&ln_z, p, RoundingMode::None);
        s = s.sub(self, p, RoundingMode::None);
        let two_pi = two_c(p).mul(&pi_c(p, cc), p, RoundingMode::None);
        let ln_two_pi = two_pi.ln(p, RoundingMode::None, cc);
        s = s.add(
            &ln_two_pi.mul(&half_c(p), p, RoundingMode::None),
            p,
            RoundingMode::None,
        );
        let bs = even_bernoulli_numbers(GAMMA_STIRLING_TERMS, p);
        let mut zpow = self.clone();
        let mut prev_e = Exponent::MIN;
        for (k, b) in bs.iter().enumerate() {
            let k = k + 1;
            let two_k = ExactNum::from_u32((2 * k) as u32, p);
            let two_k_m1 = ExactNum::from_u32((2 * k - 1) as u32, p);
            let den_r = two_k.mul(&two_k_m1, p, RoundingMode::None);
            let den = ExactComplex::from_real(den_r, p).mul(&zpow, p, RoundingMode::None);
            let term = ExactComplex::from_real(b.clone(), p).div(&den, p, RoundingMode::None);
            s = s.add(&term, p, RoundingMode::None);
            let e = term
                .abs(p, RoundingMode::None)
                .exponent()
                .unwrap_or(Exponent::MIN);
            if term_negligible(&term, p) {
                break;
            }
            if k > 2 && e > prev_e {
                break;
            }
            prev_e = e;
            zpow = zpow
                .mul(self, p, RoundingMode::None)
                .mul(self, p, RoundingMode::None);
        }
        s
    }

    pub(crate) fn digamma_at(&self, p: usize, cc: &mut Consts) -> Self {
        if !re_positive(self) && !self.re().is_zero() {
            let one = ExactComplex::one(p);
            let omz = one.sub(self, p, RoundingMode::None);
            let psi = omz.digamma_positive(p, cc);
            let piz = pi_c(p, cc).mul(self, p, RoundingMode::None);
            let cot = piz.cos(p, RoundingMode::None, cc).div(
                &piz.sin(p, RoundingMode::None, cc),
                p,
                RoundingMode::None,
            );
            return psi.sub(
                &pi_c(p, cc).mul(&cot, p, RoundingMode::None),
                p,
                RoundingMode::None,
            );
        }
        self.digamma_positive(p, cc)
    }

    fn digamma_positive(&self, p: usize, cc: &mut Consts) -> Self {
        let one = ExactComplex::one(p);
        let mut z = self.clone();
        let mut acc = ExactComplex::zero(p);
        while abs_needs_shift(&z, p, digamma_shift_exp(p)) {
            let rec = one.div(&z, p, RoundingMode::None);
            acc = acc.sub(&rec, p, RoundingMode::None);
            z = z.add(&one, p, RoundingMode::None);
        }
        acc.add(&z.digamma_asymp(p, cc), p, RoundingMode::None)
    }

    fn digamma_asymp(&self, p: usize, cc: &mut Consts) -> Self {
        let ln_z = self.ln(p, RoundingMode::None, cc);
        let two_z = two_c(p).mul(self, p, RoundingMode::None);
        let half_inv = ExactComplex::one(p).div(&two_z, p, RoundingMode::None);
        let mut s = ln_z.sub(&half_inv, p, RoundingMode::None);
        let z2 = self.mul(self, p, RoundingMode::None);
        let mut zp = ExactComplex::one(p);
        let bs = even_bernoulli_numbers(GAMMA_STIRLING_TERMS, p);
        let mut prev_e = Exponent::MIN;
        for (k, b) in bs.iter().enumerate() {
            let k = k + 1;
            zp = zp.mul(&z2, p, RoundingMode::None);
            let two_k = ExactComplex::from_real(ExactNum::from_u32((2 * k) as u32, p), p);
            let den = two_k.mul(&zp, p, RoundingMode::None);
            let term = ExactComplex::from_real(b.clone(), p).div(&den, p, RoundingMode::None);
            s = s.sub(&term, p, RoundingMode::None);
            let e = term
                .abs(p, RoundingMode::None)
                .exponent()
                .unwrap_or(Exponent::MIN);
            if term_negligible(&term, p) {
                break;
            }
            if k > 2 && e > prev_e {
                break;
            }
            prev_e = e;
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::ops::Neg;

    fn near_bits(a: &ExactNum, b: &ExactNum, p: usize, slack: Exponent) -> bool {
        let d = ExactNum::sub_ref(a, b, p, RoundingMode::None).abs();
        d.is_zero() || d.exponent().is_some_and(|e| e < -((p as Exponent) / slack))
    }

    fn near(a: &ExactNum, b: &ExactNum, p: usize) -> bool {
        near_bits(a, b, p, 4)
    }

    fn cnear(a: &ExactComplex, b: &ExactComplex, p: usize) -> bool {
        near(a.re(), b.re(), p) && near(a.im(), b.im(), p)
    }

    fn cnear_bits(a: &ExactComplex, b: &ExactComplex, p: usize, slack: Exponent) -> bool {
        near_bits(a.re(), b.re(), p, slack) && near_bits(a.im(), b.im(), p, slack)
    }

    fn tiny(x: &ExactNum, p: usize) -> bool {
        x.is_zero() || x.exponent().is_some_and(|e| e < -((p as Exponent) / 4))
    }

    /// Independent real series \(\mathrm{erfi}(x)=\frac{2}{\sqrt\pi}\sum x^{2n+1}/(n!(2n+1))\).
    fn erfi_real(mk_x: impl Fn() -> ExactNum, p: usize, cc: &mut Consts) -> ExactNum {
        let pi = cc.pi(p, RoundingMode::None);
        let sqrt_pi = pi.sqrt(p, RoundingMode::None);
        let scale = ExactNum::from_u8(2, p).div(&sqrt_pi, p, RoundingMode::None);
        let x2 = mk_x().sqr(p, RoundingMode::None);
        let mut xpow = mk_x();
        let mut nfact = ExactNum::from_u8(1, p);
        let mut sum = mk_x();
        for n in 1..=series_term_cap(p) {
            nfact = nfact.mul(&ExactNum::from_u32(n as u32, p), p, RoundingMode::None);
            xpow = xpow.mul(&x2, p, RoundingMode::None);
            let den = ExactNum::mul_ref(
                &nfact,
                &ExactNum::from_u32((2 * n + 1) as u32, p),
                p,
                RoundingMode::None,
            );
            let t = ExactNum::div_ref(&xpow, &den, p, RoundingMode::None);
            sum = sum.add(&t, p, RoundingMode::None);
            if t.is_zero() || t.exponent().is_some_and(|e| (e as i64) + (p as i64) < 0) {
                break;
            }
        }
        scale.mul(&sum, p, RoundingMode::None)
    }

    #[test]
    fn test_complex_erf_golds() -> Result<(), crate::Error> {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;

        let z0 = ExactComplex::zero(p);
        let e0 = z0.erf(p, rm, &mut cc);
        assert!(tiny(e0.re(), p) && tiny(e0.im(), p));

        let one = ExactComplex::one(p);
        let e1 = one.erf(p, rm, &mut cc);
        let r1 = ExactNum::from_u8(1, p).erf(p, rm, &mut cc);
        assert!(near(e1.re(), &r1, p));
        assert!(tiny(e1.im(), p));

        let z = ExactComplex::new(ExactNum::from_u8(1, p), ExactNum::from_u8(1, p));
        let ez = z.erf(p, rm, &mut cc);
        let em = z.neg_via_mul().erf(p, rm, &mut cc);
        assert!(cnear(&ez, &neg_c(em), p));

        let one_c = ExactComplex::one(p);
        let erfc_z = z.erfc(p, rm, &mut cc);
        let id = one_c.sub(&ez, p, rm);
        assert!(cnear(&erfc_z, &id, p));

        let i = ExactComplex::i(p);
        let ei = i.erf(p, rm, &mut cc);
        let erfi1 = erfi_real(|| ExactNum::from_u8(1, p), p, &mut cc);
        assert!(tiny(ei.re(), p));
        assert!(near(ei.im(), &erfi1, p));

        let h = ExactNum::from_u8(2, p).powsi(-((p as i64) / 8), p, rm);
        let hc = ExactComplex::from_real(h, p);
        let zp = ExactComplex::add_ref(&z, &hc, p, rm);
        let zm = ExactComplex::sub_ref(&z, &hc, p, rm);
        let num = zp.erf(p, rm, &mut cc).sub(&zm.erf(p, rm, &mut cc), p, rm);
        let two_h = hc.mul(&two_c(p), p, rm);
        let deriv = num.div(&two_h, p, rm);
        let z2 = z.mul(&z, p, rm);
        let expm = neg_c(z2).exp(p, rm, &mut cc);
        let two_over = ExactNum::from_u8(2, p).div(&cc.pi(p, rm).sqrt(p, rm), p, rm);
        let expect = ExactComplex::from_real(two_over, p).mul(&expm, p, rm);
        assert!(cnear_bits(&deriv, &expect, p, 8));

        let nan = ExactComplex::new(crate::NAN, ExactNum::new(p));
        assert!(nan.erf(p, rm, &mut cc).is_nan());
        assert!(nan.erfc(p, rm, &mut cc).is_nan());

        let big = ExactComplex::from_real(ExactNum::from_u8(10, p), p);
        let eb = big.erf(p, rm, &mut cc);
        assert!(near(eb.re(), &ExactNum::from_u8(1, p), p));
        assert!(tiny(eb.im(), p));
        assert!(cnear(&eb, &neg_c(neg_c(big).erf(p, rm, &mut cc)), p));
        let far = ExactComplex::from_real(ExactNum::from_u8(20, p), p);
        let ef = far.erf(p, rm, &mut cc);
        assert!(near(ef.re(), &ExactNum::from_u8(1, p), p));
        assert!(tiny(ef.im(), p));
        Ok(())
    }

    #[test]
    fn test_complex_gamma_golds() -> Result<(), crate::Error> {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;

        let one = ExactComplex::one(p);
        let g1 = one.gamma(p, rm, &mut cc);
        assert!(cnear(&g1, &one, p));

        let two = two_c(p);
        let g2 = two.gamma(p, rm, &mut cc);
        assert!(cnear(&g2, &one, p));

        let half = half_c(p);
        let ghalf = half.gamma(p, rm, &mut cc);
        let sqrt_pi = cc.pi(p, rm).sqrt(p, rm);
        assert!(near(ghalf.re(), &sqrt_pi, p));
        assert!(tiny(ghalf.im(), p));

        let five = ExactComplex::from_real(ExactNum::from_u8(5, p), p);
        let g5 = five.gamma(p, rm, &mut cc);
        let tf = ExactComplex::from_real(ExactNum::from_u8(24, p), p);
        assert!(cnear(&g5, &tf, p));

        let z = ExactComplex::new(
            ExactNum::from_u8(1, p).div(&ExactNum::from_u8(3, p), p, rm),
            ExactNum::from_u8(2, p).div(&ExactNum::from_u8(5, p), p, rm),
        );
        let gz = z.gamma(p, rm, &mut cc);
        let omz = ExactComplex::one(p).sub(&z, p, rm);
        let gom = omz.gamma(p, rm, &mut cc);
        let lhs = gz.mul(&gom, p, rm);
        let piz = pi_c(p, &mut cc).mul(&z, p, rm);
        let rhs = pi_c(p, &mut cc).div(&piz.sin(p, rm, &mut cc), p, rm);
        assert!(cnear(&lhs, &rhs, p));

        let lg1 = one.ln_gamma(p, rm, &mut cc);
        assert!(tiny(lg1.re(), p) && tiny(lg1.im(), p));

        let lgh = half.ln_gamma(p, rm, &mut cc);
        let ln_sqrt_pi = sqrt_pi.ln(p, rm, &mut cc);
        assert!(near(lgh.re(), &ln_sqrt_pi, p));
        assert!(tiny(lgh.im(), p));

        let psi1 = one.digamma(p, rm, &mut cc);
        let neg_g = cc.euler_gamma(p, rm).neg();
        assert!(near(psi1.re(), &neg_g, p));
        assert!(tiny(psi1.im(), p));

        let zp1 = ExactComplex::one(p).add(&z, p, rm);
        let dpsi = zp1
            .digamma(p, rm, &mut cc)
            .sub(&z.digamma(p, rm, &mut cc), p, rm);
        let rec = one.div(&z, p, rm);
        assert!(cnear(&dpsi, &rec, p));

        for n in [0i8, -1, -2] {
            let pole = ExactComplex::from_real(ExactNum::from_i8(n, p), p);
            assert!(pole.gamma(p, rm, &mut cc).is_nan(), "gamma pole {n}");
            assert!(pole.ln_gamma(p, rm, &mut cc).is_nan(), "ln_gamma pole {n}");
            assert!(pole.digamma(p, rm, &mut cc).is_nan(), "digamma pole {n}");
        }

        let mk_h = || ExactNum::from_u8(2, p).powsi(-((p as i64) / 8), p, rm);
        let hc = ExactComplex::from_real(mk_h(), p);
        let num = ExactComplex::from_real(mk_h(), p)
            .add(&z, p, rm)
            .ln_gamma(p, rm, &mut cc)
            .sub(
                &ExactComplex::from_real(mk_h(), p)
                    .neg()
                    .add(&z, p, rm)
                    .ln_gamma(p, rm, &mut cc),
                p,
                rm,
            );
        let deriv = num.div(&hc.mul(&two_c(p), p, rm), p, rm);
        let psi = z.digamma(p, rm, &mut cc);
        assert!(cnear_bits(&deriv, &psi, p, 8));

        let nan = ExactComplex::new(crate::NAN, ExactNum::new(p));
        assert!(nan.gamma(p, rm, &mut cc).is_nan());
        Ok(())
    }
}
