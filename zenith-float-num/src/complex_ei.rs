//! Complex `Ei`, `Si`, `Ci`, `li`, and Fresnel integrals.

use crate::common::util::round_p;
use crate::complex_ball::CB;
use crate::complex_ball::RB;
use crate::complex_special::half_c;
use crate::complex_special::nan_pair;
use crate::complex_special::neg_c;
use crate::complex_special::pi_c;
use crate::complex_special::term_negligible;
use crate::complex_special::two_c;
use crate::complex_special::ziv_complex;
use crate::defs::Exponent;
use crate::Consts;
use crate::Error;
use crate::ExactComplex;
use crate::ExactNum;
use crate::RoundingMode;
use crate::WORD_BIT_SIZE;
use core::ops::Neg;

/// Use the power series for \(\mathrm{Ei}\) when \(\lvert z\rvert\) is below this, whatever the precision.
const EI_SERIES_THRESHOLD: u32 = 16;

fn abs_below(z: &ExactComplex, bound: u32, p: usize) -> bool {
    let a = z.abs(p, RoundingMode::None);
    let b = ExactNum::from_u32(bound, p);
    matches!(a.compare(&b), Some(c) if c < 0)
}

/// The asymptotic series' smallest term is about \(e^{-\lvert z\rvert}\), i.e. \(1.44\lvert z\rvert\)
/// bits, so it is used only when \(\lvert z\rvert \ge 0.7(p+96)+8\) (covering the Ziv guard
/// bits). Below that the power series is used with extra precision for its cancellation.
fn use_ei_series(z: &ExactComplex, dest_p: usize) -> bool {
    if abs_below(z, EI_SERIES_THRESHOLD, dest_p) {
        return true;
    }
    let t = (dest_p.saturating_add(96).saturating_mul(7) / 10).saturating_add(8);
    abs_below(z, t.min(u32::MAX as usize) as u32, dest_p)
}

/// Upper bound on the bits cancelled by the \(\mathrm{Ei}\) power series: its largest term is
/// about \(e^{\lvert z\rvert}\) and the result about \(e^{\mathrm{Re}\,z}/\lvert z\rvert\).
fn ei_series_loss_bits(z: &ExactComplex) -> usize {
    let p = 64;
    let az = z.abs(p, RoundingMode::None);
    let d =
        az.sub(z.re(), p, RoundingMode::None)
            .add(&ExactNum::from_u32(1, p), p, RoundingMode::None);
    match d.exponent() {
        Some(e) if e > 0 => {
            // d < 2^e and 1.4427 d < 1.5 * 2^e
            let e = (e as u32).min(40);
            (3usize << e) / 2 + 16
        }
        _ => 16,
    }
}

impl ExactComplex {
    /// Exponential integral \(\mathrm{Ei}(z)\). Cut on \((-\infty,0]\); pole at \(0\) → NaN.
    ///
    /// # Precision
    ///
    /// - Algorithm: power series (with extra precision for cancellation) for `|z| < max(16, 0.7(p+96)+8)`;
    ///   asymptotic `e^z/z Σ k!/z^k + iπ·sgn(Im z)` otherwise.
    /// - Bound: Ziv on each part (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn ei(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.ei_word(p, rm, cc))
    }

    /// Implementation of [`Self::ei`] on operands without inexact history.
    fn ei_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return ExactComplex::new(self.re().clone(), self.im().clone());
        }
        if self.re().is_zero() && self.im().is_zero() {
            return nan_pair(Error::InvalidArgument);
        }
        let dest = round_p(p);
        ziv_complex(dest, rm, |pw| self.ei_at(pw, dest, cc))
    }

    /// Sine integral \(\mathrm{Si}(z)=(E_i(iz)-E_i(-iz))/(2i)-\pi/2\).
    ///
    /// # Precision
    ///
    /// - Algorithm: via [`Self::ei`].
    /// - Bound: Ziv on each part (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn si(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.si_word(p, rm, cc))
    }

    /// Implementation of [`Self::si`] on operands without inexact history.
    fn si_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return ExactComplex::new(self.re().clone(), self.im().clone());
        }
        // Si(0) = 0, odd, with C99 signed zeros. Must not allocate dest.
        if self.re().is_zero() && self.im().is_zero() {
            return crate::complex::compact_complex_zero(
                self.re().is_negative(),
                self.im().is_negative(),
            );
        }
        let dest = round_p(p);
        // Si is odd with real Taylor coefficients: Si(iy) is purely imaginary and Si(x) is real.
        // Pin the exactly-zero part so the Ziv test does not chase rounding noise.
        // The pinned zero carries the sign of the limit: Re Si(x + iy) ≈ x sinh(y)/y, so it is
        // the sign of x; Im Si(x + iy) ≈ y sin(x)/x, so it is sign(y) · sign(sin x / x).
        let re_zero = self.re().is_zero();
        let im_zero = self.im().is_zero();
        let im_zero_neg = im_zero && !re_zero && {
            let sx = self.re().sin(WORD_BIT_SIZE, RoundingMode::ToEven, cc);
            self.im().is_negative() ^ (sx.is_negative() ^ self.re().is_negative())
        };
        let zero = |neg: bool, pw: usize| {
            let mut z = ExactNum::new(pw);
            if neg {
                z = z.neg();
            }
            z
        };
        let re_neg = self.re().is_negative();
        ziv_complex(dest, rm, |pw| {
            let v = self.si_at(pw, dest, cc);
            match (re_zero, im_zero) {
                (true, false) => ExactComplex::new(zero(re_neg, pw), v.im().clone()),
                (false, true) => ExactComplex::new(v.re().clone(), zero(im_zero_neg, pw)),
                _ => v,
            }
        })
    }

    /// Cosine integral \(\mathrm{Ci}(z)=\gamma+\ln z+\int_0^z(\cos t-1)/t\,dt\) (principal branch,
    /// cut on \((-\infty,0]\), matching mpmath `ci`). Pole at \(0\) → NaN.
    ///
    /// # Precision
    ///
    /// - Algorithm: via [`Self::ei`].
    /// - Bound: Ziv on each part (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn ci(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.ci_word(p, rm, cc))
    }

    /// Implementation of [`Self::ci`] on operands without inexact history.
    fn ci_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return ExactComplex::new(self.re().clone(), self.im().clone());
        }
        if self.re().is_zero() && self.im().is_zero() {
            return nan_pair(Error::InvalidArgument);
        }
        let dest = round_p(p);
        // ZF-37: on the positive real axis Ci(x ± 0i) = Ci(x) ± 0i (real Taylor coefficients
        // off the cut); the exactly-zero imaginary part is pinned, as in `si_word`, so Ziv does
        // not chase rounding noise. Its sign is the limit's: Im Ci(x + iy) ≈ y cos(x)/x.
        if self.im().is_zero() && self.re().is_positive() {
            let cx = self.re().cos(WORD_BIT_SIZE, RoundingMode::ToEven, cc);
            // The imaginary part is an exact signed zero; do not allocate dest for it.
            let z = crate::complex::compact_signed_zero(self.im().is_negative() ^ cx.is_negative());
            return ExactComplex::new(self.re().ci(dest, rm, cc), z);
        }
        ziv_complex(dest, rm, |pw| self.ci_at(pw, dest, cc))
    }

    /// Logarithmic integral \(\mathrm{li}(z)=\mathrm{Ei}(\ln z)\). Cut on \((-\infty,1]\); pole at \(1\) → NaN.
    ///
    /// # Precision
    ///
    /// - Algorithm: [`Self::ei`] of `ln z`.
    /// - Bound: Ziv on each part (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn li(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.li_word(p, rm, cc))
    }

    /// Implementation of [`Self::li`] on operands without inexact history.
    fn li_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return ExactComplex::new(self.re().clone(), self.im().clone());
        }
        if self.re().is_zero() && self.im().is_zero() {
            return nan_pair(Error::InvalidArgument);
        }
        let dest = round_p(p);
        ziv_complex(dest, rm, |pw| self.li_at(pw, dest, cc))
    }

    /// Fresnel sine integral \(S(z)\). Entire.
    ///
    /// # Precision
    ///
    /// - Algorithm: via complex `erf`; Ziv on each part (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn fresnel_s(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.fresnel_s_word(p, rm, cc))
    }

    /// Implementation of [`Self::fresnel_s`] on operands without inexact history.
    fn fresnel_s_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return ExactComplex::new(self.re().clone(), self.im().clone());
        }
        // S(0) is an exact signed zero (C99 odd, conjugate-symmetric; Im(S)
        // takes the opposite sign of Im(z)). Must not allocate dest.
        if self.re().is_zero() && self.im().is_zero() {
            return crate::complex::compact_complex_zero(
                self.re().is_negative(),
                !self.im().is_negative(),
            );
        }
        let dest = round_p(p);
        if !finite_parts(self) {
            let pw = dest + WORD_BIT_SIZE;
            let mut v = self.fresnel_pair_plain(pw, cc).0;
            let (mut re, mut im) = (v.re().clone(), v.im().clone());
            let _ = re.set_precision(dest, rm);
            let _ = im.set_precision(dest, rm);
            v = ExactComplex::new(re, im);
            return v;
        }
        if let Some(v) = self.fresnel_axis(true, dest, rm, cc) {
            return v;
        }
        crate::complex_ball::ziv_complex_ball(dest, rm, |q| {
            self.fresnel_pair_ball(q, cc).map(|pair| pair.0)
        })
    }

    /// Fresnel cosine integral \(C(z)\). Entire.
    ///
    /// # Precision
    ///
    /// - Algorithm: via complex `erf`; Ziv on each part (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn fresnel_c(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.fresnel_c_word(p, rm, cc))
    }

    /// Implementation of [`Self::fresnel_c`] on operands without inexact history.
    fn fresnel_c_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return ExactComplex::new(self.re().clone(), self.im().clone());
        }
        // C(0) is an exact signed zero that copies the operand's signed zeros
        // (C99 odd, conjugate-symmetric; C(z) ~ z). Must not allocate dest.
        if self.re().is_zero() && self.im().is_zero() {
            return crate::complex::compact_complex_zero(
                self.re().is_negative(),
                self.im().is_negative(),
            );
        }
        let dest = round_p(p);
        if !finite_parts(self) {
            let pw = dest + WORD_BIT_SIZE;
            let mut v = self.fresnel_pair_plain(pw, cc).1;
            let (mut re, mut im) = (v.re().clone(), v.im().clone());
            let _ = re.set_precision(dest, rm);
            let _ = im.set_precision(dest, rm);
            v = ExactComplex::new(re, im);
            return v;
        }
        if let Some(v) = self.fresnel_axis(false, dest, rm, cc) {
            return v;
        }
        crate::complex_ball::ziv_complex_ball(dest, rm, |q| {
            self.fresnel_pair_ball(q, cc).map(|pair| pair.1)
        })
    }

    fn ei_at(&self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        if use_ei_series(self, dest_p) {
            let extra = ei_series_loss_bits(self);
            self.ei_series(work_p.saturating_add(extra), cc)
        } else {
            // Ei(z) ~ e^z/z Σ k!/z^k + iπ·sgn(Im z) for the principal-Ln branch (the same branch
            // as the series; on the cut Im z = 0, Re z < 0 the upper side, +iπ, is used).
            let a = self.ei_asymptotic(work_p, cc);
            let s = if self.im().is_positive() || (self.im().is_zero() && self.re().is_negative()) {
                1
            } else if self.im().is_negative() && !self.im().is_zero() {
                -1
            } else {
                0
            };
            if s == 0 {
                a
            } else {
                let pi = cc.pi(work_p, RoundingMode::None);
                let pi = if s < 0 { pi.neg() } else { pi };
                let (re, im) = a.into_parts();
                ExactComplex::new(re, im.add(&pi, work_p, RoundingMode::None))
            }
        }
    }

    /// \(\mathrm{Ei}(z)=\gamma+\mathrm{Ln}\,z+\sum z^n/(n\cdot n!)\). Principal \(\mathrm{Ln}\).
    fn ei_series(&self, p: usize, cc: &mut Consts) -> Self {
        let g = ExactComplex::from_real(cc.euler_gamma(p, RoundingMode::None), p);
        let lnz = self.ln(p, RoundingMode::None, cc);
        let mut term = self.clone();
        let mut sum = term.clone();
        // Terms peak near n = |z| and then decay; run until the tail is below 2^{-p}.
        let mut n: u32 = 2;
        loop {
            let nw = ExactComplex::from_real(ExactNum::from_u32(n, p), p);
            term = term
                .mul(self, p, RoundingMode::None)
                .div(&nw, p, RoundingMode::None);
            let piece = term.div(&nw, p, RoundingMode::None);
            sum = sum.add(&piece, p, RoundingMode::None);
            if term_negligible(&piece, p) {
                break;
            }
            let Some(next) = n.checked_add(1) else {
                return nan_pair(Error::PrecisionRetryExhausted);
            };
            n = next;
        }
        g.add(&lnz, p, RoundingMode::None)
            .add(&sum, p, RoundingMode::None)
    }

    /// \(\mathrm{Ei}(z)\sim e^z/z\sum k!/z^k\), stopped at the smallest term.
    fn ei_asymptotic(&self, p: usize, cc: &mut Consts) -> Self {
        let pre = self
            .exp(p, RoundingMode::None, cc)
            .div(self, p, RoundingMode::None);
        let mut term = ExactComplex::one(p);
        let mut s = term.clone();
        let mut prev_e = Exponent::MIN;
        let mut k: u32 = 1;
        loop {
            let kk = ExactComplex::from_real(ExactNum::from_u32(k, p), p);
            term = term
                .mul(&kk, p, RoundingMode::None)
                .div(self, p, RoundingMode::None);
            let e = term
                .abs(p, RoundingMode::None)
                .exponent()
                .unwrap_or(Exponent::MIN);
            if k > 1 && e > prev_e {
                break;
            }
            prev_e = e;
            s = s.add(&term, p, RoundingMode::None);
            if term_negligible(&term, p) {
                break;
            }
            let Some(next) = k.checked_add(1) else {
                return nan_pair(Error::PrecisionRetryExhausted);
            };
            k = next;
        }
        pre.mul(&s, p, RoundingMode::None)
    }

    fn si_at(&self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        if self.re().is_zero() && self.im().is_zero() {
            return ExactComplex::zero(work_p);
        }
        // \(\mathrm{Si}\) is odd; the \(\mathrm{Ei}\) combination minus \(\pi/2\) is the
        // principal value for \(\mathrm{Re}\,z\ge 0\).
        // Si is odd and entire. Decide the half-plane by value, not by the zero's sign: with
        // re = −0 the old test negated to re = −0 again (−1 · +0 is −0 under neg_via_mul) and
        // recursed until the stack overflowed (ZF-72).
        let re_neg = !self.re().is_zero() && self.re().is_negative();
        if re_neg || (self.re().is_zero() && self.im().is_negative()) {
            return neg_c(self.neg_via_mul().si_at(work_p, dest_p, cc));
        }
        let iz = ExactComplex::i(work_p).mul(self, work_p, RoundingMode::None);
        let e_plus = iz.ei_at(work_p, dest_p, cc);
        let e_minus = neg_c(iz).ei_at(work_p, dest_p, cc);
        let two_i = two_c(work_p).mul(&ExactComplex::i(work_p), work_p, RoundingMode::None);
        let half_pi = pi_c(work_p, cc).mul(&half_c(work_p), work_p, RoundingMode::None);
        e_plus
            .sub(&e_minus, work_p, RoundingMode::None)
            .div(&two_i, work_p, RoundingMode::None)
            .sub(&half_pi, work_p, RoundingMode::None)
    }

    fn ci_at(&self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        let iz = ExactComplex::i(work_p).mul(self, work_p, RoundingMode::None);
        let e_plus = iz.ei_at(work_p, dest_p, cc);
        let e_minus = iz.neg_via_mul().ei_at(work_p, dest_p, cc);
        // (Ei(iz) + Ei(-iz)) / 2 = γ + (ln(iz) + ln(-iz)) / 2 + (entire part). The principal
        // Ci(z) = γ + ln z + (entire part) (DLMF 6.2.11), so add i·(arg z − (arg(iz) + arg(−iz)) / 2),
        // which is 0 for Re z > 0 and ±π for Re z < 0.
        let half_sum = e_plus.add(&e_minus, work_p, RoundingMode::None).div(
            &two_c(work_p),
            work_p,
            RoundingMode::None,
        );
        let arg_z = self.ln(work_p, RoundingMode::None, cc).im().clone();
        let arg_p = iz.ln(work_p, RoundingMode::None, cc).im().clone();
        let (_, arg_m) = neg_c(iz).ln(work_p, RoundingMode::None, cc).into_parts();
        let corr = arg_z.sub(
            &arg_p
                .add(&arg_m, work_p, RoundingMode::None)
                .ldexp(-1, work_p, RoundingMode::None),
            work_p,
            RoundingMode::None,
        );
        let (re, im) = half_sum.into_parts();
        ExactComplex::new(re, im.add(&corr, work_p, RoundingMode::None))
    }

    fn li_at(&self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        let one = ExactComplex::one(work_p);
        if self.im().is_zero() && self.re().compare(one.re()) == Some(0) {
            return nan_pair(Error::InvalidArgument);
        }
        self.ln(work_p, RoundingMode::None, cc)
            .ei_at(work_p, dest_p, cc)
    }

    /// One-shot `(S, C)` from complex `erf` (non-finite operands keep its NaN results).
    fn fresnel_pair_plain(&self, p: usize, cc: &mut Consts) -> (Self, Self) {
        let sqrt_pi =
            ExactComplex::from_real(cc.pi(p, RoundingMode::None).sqrt(p, RoundingMode::None), p);
        let scale = sqrt_pi
            .mul(self, p, RoundingMode::None)
            .mul(&half_c(p), p, RoundingMode::None);
        let one = ExactComplex::one(p);
        let i = ExactComplex::i(p);
        let one_p_i = ExactComplex::one(p).add(&i, p, RoundingMode::None);
        let one_m_i = one.sub(&i, p, RoundingMode::None);
        let erf_m = one_m_i
            .mul(&scale, p, RoundingMode::None)
            .erf(p, RoundingMode::None, cc);
        let erf_p = one_p_i
            .mul(&scale, p, RoundingMode::None)
            .erf(p, RoundingMode::None, cc);
        let c_plus_is =
            one_p_i
                .mul(&half_c(p), p, RoundingMode::None)
                .mul(&erf_m, p, RoundingMode::None);
        let c_minus_is =
            one_m_i
                .mul(&half_c(p), p, RoundingMode::None)
                .mul(&erf_p, p, RoundingMode::None);
        let two_i = two_c(p).mul(&i, p, RoundingMode::None);
        // `c+ - c-` is `-c- + c+`: both stay borrowed until the sum below takes `c+`.
        let s = ExactComplex::sub_ref(&c_plus_is, &c_minus_is, p, RoundingMode::None).div(
            &two_i,
            p,
            RoundingMode::None,
        );
        let c = c_plus_is.add(&c_minus_is, p, RoundingMode::None).mul(
            &half_c(p),
            p,
            RoundingMode::None,
        );
        (s, c)
    }

    /// Fresnel `S`/`C` on an axis, from the real functions (DLMF 7.4.7–7.4.8: `S(iy) = −i S(y)`,
    /// `C(iy) = i C(y)`). The zero part is signed as its limit from the side of the operand's
    /// signed zero: `Im S(x ± iε) = ±ε sin(πx²/2) + O(ε³)` (and `∓ε³π cos(πx²/2)/6` where the sine
    /// vanishes), and likewise for `C` with `cos`/`sin` and for the imaginary axis.
    /// `None` off the axes and for non-finite parts.
    fn fresnel_axis(
        &self,
        sine: bool,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Option<Self> {
        let (re, im) = (self.re(), self.im());
        if re.is_nan() || im.is_nan() || re.is_inf() || im.is_inf() {
            return None;
        }
        let real_fn = |x: &ExactNum, rm: RoundingMode, cc: &mut Consts| {
            if sine {
                x.fresnel_s(p, rm, cc)
            } else {
                x.fresnel_c(p, rm, cc)
            }
        };
        let zero = |neg: bool| if neg { ExactNum::new(p).neg() } else { ExactNum::new(p) };
        if im.is_zero() {
            let (s, c) = pi_half_sq_signs(re)?;
            let sigma = if im.is_negative() { -1 } else { 1 };
            let sgn = if sine {
                if s != 0 {
                    sigma * s
                } else {
                    -sigma * c
                }
            } else if c != 0 {
                sigma * c
            } else {
                sigma * s
            };
            return Some(ExactComplex::new(real_fn(re, rm, cc), zero(sgn < 0)));
        }
        if re.is_zero() {
            let (s, c) = pi_half_sq_signs(im)?;
            let sigma = if re.is_negative() { -1 } else { 1 };
            let sgn = if sine {
                if s != 0 {
                    -sigma * s
                } else {
                    sigma * c
                }
            } else if c != 0 {
                sigma * c
            } else {
                sigma * s
            };
            let v = if sine {
                let mirrored = match rm {
                    RoundingMode::Up => RoundingMode::Down,
                    RoundingMode::Down => RoundingMode::Up,
                    other => other,
                };
                real_fn(im, mirrored, cc).neg()
            } else {
                real_fn(im, rm, cc)
            };
            return Some(ExactComplex::new(zero(sgn < 0), v));
        }
        None
    }

    /// `(S, C)` as balls: `S ± iC`-type combinations of `erf((1 ∓ i)√π z/2)` (DLMF 7.5.7),
    /// with the erf arguments' rounding carried by [`ExactComplex::erf_of_ball`].
    fn fresnel_pair_ball(&self, q: usize, cc: &mut Consts) -> Option<(CB, CB)> {
        let sqrt_pi = RB::rounded(cc.pi(q, RoundingMode::None), q)?.sqrt(q)?;
        let scale = CB::exact(self).mul_real(&sqrt_pi, q)?.scale2(-1)?;
        let one = || RB::exact(ExactNum::from_u8(1, WORD_BIT_SIZE));
        let one_p_i = CB::from_parts(one(), one());
        let one_m_i = CB::from_parts(one(), one().neg());
        let erf_m = ExactComplex::erf_of_ball(&one_m_i.mul(&scale, q)?, q, cc)?;
        let erf_p = ExactComplex::erf_of_ball(&one_p_i.mul(&scale, q)?, q, cc)?;
        let c_plus_is = one_p_i.scale2(-1)?.mul(&erf_m, q)?;
        let c_minus_is = one_m_i.scale2(-1)?.mul(&erf_p, q)?;
        // S = (c+ − c−)/(2i) = −i (c+ − c−)/2, C = (c+ + c−)/2.
        let s = c_plus_is.sub(&c_minus_is, q)?.mul_neg_i().scale2(-1)?;
        let c = c_plus_is.add(&c_minus_is, q)?.scale2(-1)?;
        Some((s, c))
    }
}

fn finite_parts(z: &ExactComplex) -> bool {
    !z.re().is_nan() && !z.im().is_nan() && !z.re().is_inf() && !z.im().is_inf()
}

/// Signs (−1, 0, 1) of `sin(πx²/2)` and `cos(πx²/2)`, exactly, from `t = x²/2 mod 2`.
fn pi_half_sq_signs(x: &ExactNum) -> Option<(i32, i32)> {
    let t = x.mul_full_prec(x).ldexp(
        -1,
        x.mantissa_max_bit_len()?.saturating_mul(2) + WORD_BIT_SIZE,
        RoundingMode::None,
    );
    let r = t.rem(&ExactNum::from_u8(2, WORD_BIT_SIZE));
    let at = |k: u8, d: i64| {
        let v = ExactNum::from_u8(k, WORD_BIT_SIZE).ldexp(d, WORD_BIT_SIZE, RoundingMode::None);
        r.compare(&v)
    };
    let c0 = r.compare(&ExactNum::new(WORD_BIT_SIZE))?; // r ≥ 0
    let c_half = at(1, -1)?;
    let c_one = at(1, 0)?;
    let c_3half = at(3, -1)?;
    let sin = if c0 == 0 || c_one == 0 {
        0
    } else if c_one < 0 {
        1
    } else {
        -1
    };
    let cos = if c_half == 0 || c_3half == 0 {
        0
    } else if c_half < 0 || c_3half > 0 {
        1
    } else {
        -1
    };
    Some((sin, cos))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::complex_special::neg_c;
    use crate::complex_special::pi_c;
    use crate::complex_special::two_c;

    fn near(a: &ExactNum, b: &ExactNum, p: usize) -> bool {
        let d = ExactNum::sub_ref(a, b, p, RoundingMode::None).abs();
        d.is_zero() || d.exponent().is_some_and(|e| e < -((p as Exponent) / 4))
    }

    fn cnear(a: &ExactComplex, b: &ExactComplex, p: usize) -> bool {
        near(a.re(), b.re(), p) && near(a.im(), b.im(), p)
    }

    fn cnear_bits(a: &ExactComplex, b: &ExactComplex, p: usize, slack: Exponent) -> bool {
        let dr = ExactNum::sub_ref(a.re(), b.re(), p, RoundingMode::None).abs();
        let di = ExactNum::sub_ref(a.im(), b.im(), p, RoundingMode::None).abs();
        (dr.is_zero()
            || dr
                .exponent()
                .is_some_and(|e| e < -((p as Exponent) / slack)))
            && (di.is_zero()
                || di
                    .exponent()
                    .is_some_and(|e| e < -((p as Exponent) / slack)))
    }

    fn tiny(x: &ExactNum, p: usize) -> bool {
        x.is_zero() || x.exponent().is_some_and(|e| e < -((p as Exponent) / 4))
    }

    #[test]
    fn test_complex_ei_golds() -> Result<(), crate::Error> {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;

        let one = ExactComplex::one(p);
        let ei1 = one.ei(p, rm, &mut cc);
        let r1 = ExactNum::from_u8(1, p).ei(p, rm, &mut cc);
        assert!(near(ei1.re(), &r1, p));
        assert!(tiny(ei1.im(), p));

        let z0 = ExactComplex::zero(p);
        let s0 = z0.si(p, rm, &mut cc);
        assert!(tiny(s0.re(), p) && tiny(s0.im(), p));
        assert!(z0.ci(p, rm, &mut cc).is_nan());

        let z = ExactComplex::new(ExactNum::from_u8(1, p), half_c(p).into_parts().0);
        let sz = z.si(p, rm, &mut cc);
        assert!(cnear(&sz, &neg_c(z.neg_via_mul().si(p, rm, &mut cc)), p));

        let fs0 = z0.fresnel_s(p, rm, &mut cc);
        let fc0 = z0.fresnel_c(p, rm, &mut cc);
        assert!(tiny(fs0.re(), p) && tiny(fs0.im(), p));
        assert!(tiny(fc0.re(), p) && tiny(fc0.im(), p));

        let mk_h = || ExactNum::from_u8(2, p).powsi(-((p as i64) / 8), p, rm);
        let hc = ExactComplex::from_real(mk_h(), p);
        let num_ei = ExactComplex::from_real(mk_h(), p)
            .add(&z, p, rm)
            .ei(p, rm, &mut cc)
            .sub(
                &ExactComplex::from_real(mk_h(), p)
                    .neg()
                    .add(&z, p, rm)
                    .ei(p, rm, &mut cc),
                p,
                rm,
            );
        let deriv_ei = num_ei.div(&hc.mul(&two_c(p), p, rm), p, rm);
        let expect_ei = z.exp(p, rm, &mut cc).div(&z, p, rm);
        assert!(cnear_bits(&deriv_ei, &expect_ei, p, 8));

        let num_si = ExactComplex::from_real(mk_h(), p)
            .add(&z, p, rm)
            .si(p, rm, &mut cc)
            .sub(
                &ExactComplex::from_real(mk_h(), p)
                    .neg()
                    .add(&z, p, rm)
                    .si(p, rm, &mut cc),
                p,
                rm,
            );
        let deriv_si = num_si.div(&hc.mul(&two_c(p), p, rm), p, rm);
        let expect_si = z.sin(p, rm, &mut cc).div(&z, p, rm);
        assert!(cnear_bits(&deriv_si, &expect_si, p, 8));

        let e = ExactComplex::from_real(cc.e(p, rm), p);
        let lie = e.li(p, rm, &mut cc);
        let rli = cc.e(p, rm).li(p, rm, &mut cc);
        assert!(near(lie.re(), &rli, p));
        assert!(tiny(lie.im(), p));

        let above = ExactComplex::new(ExactNum::from_i8(-1, p), ExactNum::new(p));
        let below = ExactComplex::new(ExactNum::from_i8(-1, p), ExactNum::new(p).neg());
        assert!(above.im().is_positive());
        assert!(below.im().is_negative());
        let jump = above
            .ei(p, rm, &mut cc)
            .sub(&below.ei(p, rm, &mut cc), p, rm);
        let two_pi_i = two_c(p)
            .mul(&pi_c(p, &mut cc), p, rm)
            .mul(&ExactComplex::i(p), p, rm);
        assert!(cnear_bits(&jump, &two_pi_i, p, 8));

        let nan = ExactComplex::new(crate::NAN, ExactNum::new(p));
        assert!(nan.ei(p, rm, &mut cc).is_nan());
        Ok(())
    }
    fn compact_signed_zero(neg: bool) -> ExactNum {
        let mut z = ExactNum::from_word(0, crate::DEFAULT_P);
        if neg {
            z.set_sign(crate::Sign::Neg);
        }
        z
    }

    /// Exact signed zeros of S, C and Si at dest \(2^{59}\) must return immediately
    /// with a compact mantissa. C99 oddness and conjugate symmetry fix the signs;
    /// S takes the opposite imaginary sign from the operand (the existing p=64
    /// axis convention).
    #[test]
    fn zf_fresnel_and_si_signed_zero_at_dest_2_59_returns_immediately() {
        let p = 1usize << 59;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new().expect("Consts");
        for (re_neg, im_neg) in [(false, false), (false, true), (true, false), (true, true)] {
            let z = ExactComplex::new(compact_signed_zero(re_neg), compact_signed_zero(im_neg));
            let t0 = std::time::Instant::now();
            let s = z.fresnel_s(p, rm, &mut cc);
            assert!(
                t0.elapsed().as_secs() < 1,
                "fresnel_s at dest 2^59 took {:?}",
                t0.elapsed()
            );
            assert!(
                s.re().is_zero() && s.re().is_negative() == re_neg,
                "S({z:?}) re {:?}",
                s
            );
            assert!(
                s.im().is_zero() && s.im().is_negative() == !im_neg,
                "S({z:?}) im {:?}",
                s
            );

            let t0 = std::time::Instant::now();
            let c = z.fresnel_c(p, rm, &mut cc);
            assert!(
                t0.elapsed().as_secs() < 1,
                "fresnel_c at dest 2^59 took {:?}",
                t0.elapsed()
            );
            assert!(
                c.re().is_zero() && c.re().is_negative() == re_neg,
                "C({z:?}) re {:?}",
                c
            );
            assert!(
                c.im().is_zero() && c.im().is_negative() == im_neg,
                "C({z:?}) im {:?}",
                c
            );

            let t0 = std::time::Instant::now();
            let si = z.si(p, rm, &mut cc);
            assert!(
                t0.elapsed().as_secs() < 1,
                "si at dest 2^59 took {:?}",
                t0.elapsed()
            );
            assert!(
                si.re().is_zero() && si.re().is_negative() == re_neg,
                "Si({z:?}) re {:?}",
                si
            );
            assert!(
                si.im().is_zero() && si.im().is_negative() == im_neg,
                "Si({z:?}) im {:?}",
                si
            );
        }
    }
}
