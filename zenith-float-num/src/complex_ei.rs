//! Complex `Ei`, `Si`, `Ci`, `li`, and Fresnel integrals.

use crate::common::util::round_p;
use crate::complex_special::half_c;
use crate::complex_special::nan_pair;
use crate::complex_special::neg_c;
use crate::complex_special::pi_c;
use crate::complex_special::series_term_cap;
use crate::complex_special::term_negligible;
use crate::complex_special::two_c;
use crate::complex_special::ziv_complex;
use crate::Consts;
use crate::Error;
use crate::ExactComplex;
use crate::ExactNum;
use crate::RoundingMode;
use crate::WORD_BIT_SIZE;

/// Use the power series for \(\mathrm{Ei}\) when \(\lvert z\rvert\) is below this, whatever the precision.
const EI_SERIES_THRESHOLD: u32 = 16;

fn abs_below(z: &ExactComplex, bound: u32, p: usize) -> bool {
    let a = z.abs(p, RoundingMode::None);
    let b = ExactNum::from_u32(bound, p);
    matches!(a.cmp(&b), Some(c) if c < 0)
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
    pub fn ei(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
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
    pub fn si(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return ExactComplex::new(self.re().clone(), self.im().clone());
        }
        let dest = round_p(p);
        // Si is odd with real Taylor coefficients: Si(iy) is purely imaginary and Si(x) is real.
        // Pin the exactly-zero part so the Ziv test does not chase rounding noise.
        let re_zero = self.re().is_zero();
        let im_zero = self.im().is_zero();
        ziv_complex(dest, rm, |pw| {
            let v = self.si_at(pw, dest, cc);
            match (re_zero, im_zero) {
                (true, false) => ExactComplex::new(ExactNum::new(pw), v.im().clone()),
                (false, true) => ExactComplex::new(v.re().clone(), ExactNum::new(pw)),
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
    pub fn ci(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return ExactComplex::new(self.re().clone(), self.im().clone());
        }
        if self.re().is_zero() && self.im().is_zero() {
            return nan_pair(Error::InvalidArgument);
        }
        let dest = round_p(p);
        ziv_complex(dest, rm, |pw| self.ci_at(pw, dest, cc))
    }

    /// Logarithmic integral \(\mathrm{li}(z)=\mathrm{Ei}(\ln z)\). Cut on \((-\infty,1]\); pole at \(1\) → NaN.
    ///
    /// # Precision
    ///
    /// - Algorithm: [`Self::ei`] of `ln z`.
    /// - Bound: Ziv on each part (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    pub fn li(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
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
    pub fn fresnel_s(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return ExactComplex::new(self.re().clone(), self.im().clone());
        }
        let dest = round_p(p);
        ziv_complex(dest, rm, |pw| self.fresnel_s_at(pw, cc))
    }

    /// Fresnel cosine integral \(C(z)\). Entire.
    ///
    /// # Precision
    ///
    /// - Algorithm: via complex `erf`; Ziv on each part (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    pub fn fresnel_c(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return ExactComplex::new(self.re().clone(), self.im().clone());
        }
        let dest = round_p(p);
        ziv_complex(dest, rm, |pw| self.fresnel_c_at(pw, cc))
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
                ExactComplex::new(a.re().clone(), a.im().add(&pi, work_p, RoundingMode::None))
            }
        }
    }

    /// \(\mathrm{Ei}(z)=\gamma+\mathrm{Ln}\,z+\sum z^n/(n\cdot n!)\). Principal \(\mathrm{Ln}\).
    fn ei_series(&self, p: usize, cc: &mut Consts) -> Self {
        let g = ExactComplex::from_real(cc.euler_gamma(p, RoundingMode::None), p);
        let lnz = self.ln(p, RoundingMode::None, cc);
        let mut term = self.clone();
        let mut sum = term.clone();
        // Terms peak near n = |z|; allow e·|z| terms beyond the precision-based cap.
        let az_bound = self
            .abs(64, RoundingMode::None)
            .exponent()
            .map_or(
                0usize,
                |e| if e > 0 { 3usize << (e as u32).min(40) } else { 0 },
            );
        for n in 2..=series_term_cap(p).max(az_bound.saturating_add(WORD_BIT_SIZE)) {
            let nw = ExactComplex::from_real(ExactNum::from_u32(n as u32, p), p);
            term = term
                .mul(self, p, RoundingMode::None)
                .div(&nw, p, RoundingMode::None);
            let piece = term.div(&nw, p, RoundingMode::None);
            sum = sum.add(&piece, p, RoundingMode::None);
            if term_negligible(&piece, p) {
                break;
            }
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
        let mut prev_e = i32::MIN;
        for k in 1..=series_term_cap(p) {
            let kk = ExactComplex::from_real(ExactNum::from_u32(k as u32, p), p);
            term = term
                .mul(&kk, p, RoundingMode::None)
                .div(self, p, RoundingMode::None);
            let e = term
                .abs(p, RoundingMode::None)
                .exponent()
                .unwrap_or(i32::MIN);
            if k > 1 && e > prev_e {
                break;
            }
            prev_e = e;
            s = s.add(&term, p, RoundingMode::None);
            if term_negligible(&term, p) {
                break;
            }
        }
        pre.mul(&s, p, RoundingMode::None)
    }

    fn si_at(&self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        if self.re().is_zero() && self.im().is_zero() {
            return ExactComplex::zero(work_p);
        }
        // \(\mathrm{Si}\) is odd; the \(\mathrm{Ei}\) combination minus \(\pi/2\) is the
        // principal value for \(\mathrm{Re}\,z\ge 0\).
        if self.re().is_negative() || (self.re().is_zero() && self.im().is_negative()) {
            return neg_c(&neg_c(self).si_at(work_p, dest_p, cc));
        }
        let iz = ExactComplex::i(work_p).mul(self, work_p, RoundingMode::None);
        let e_plus = iz.ei_at(work_p, dest_p, cc);
        let e_minus = neg_c(&iz).ei_at(work_p, dest_p, cc);
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
        let e_minus = neg_c(&iz).ei_at(work_p, dest_p, cc);
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
        let arg_m = neg_c(&iz).ln(work_p, RoundingMode::None, cc).im().clone();
        let corr = arg_z.sub(
            &arg_p
                .add(&arg_m, work_p, RoundingMode::None)
                .ldexp(-1, work_p, RoundingMode::None),
            work_p,
            RoundingMode::None,
        );
        ExactComplex::new(
            half_sum.re().clone(),
            half_sum.im().add(&corr, work_p, RoundingMode::None),
        )
    }

    fn li_at(&self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        let one = ExactComplex::one(work_p);
        if self.im().is_zero() && self.re().cmp(one.re()) == Some(0) {
            return nan_pair(Error::InvalidArgument);
        }
        self.ln(work_p, RoundingMode::None, cc)
            .ei_at(work_p, dest_p, cc)
    }

    fn fresnel_pair(&self, p: usize, cc: &mut Consts) -> (Self, Self) {
        let sqrt_pi =
            ExactComplex::from_real(cc.pi(p, RoundingMode::None).sqrt(p, RoundingMode::None), p);
        let scale = sqrt_pi
            .mul(self, p, RoundingMode::None)
            .mul(&half_c(p), p, RoundingMode::None);
        let one = ExactComplex::one(p);
        let i = ExactComplex::i(p);
        let one_p_i = one.add(&i, p, RoundingMode::None);
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
        let c = c_plus_is.add(&c_minus_is, p, RoundingMode::None).mul(
            &half_c(p),
            p,
            RoundingMode::None,
        );
        let two_i = two_c(p).mul(&i, p, RoundingMode::None);
        let s =
            c_plus_is
                .sub(&c_minus_is, p, RoundingMode::None)
                .div(&two_i, p, RoundingMode::None);
        (s, c)
    }

    fn fresnel_s_at(&self, p: usize, cc: &mut Consts) -> Self {
        self.fresnel_pair(p, cc).0
    }

    fn fresnel_c_at(&self, p: usize, cc: &mut Consts) -> Self {
        self.fresnel_pair(p, cc).1
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::complex_special::neg_c;
    use crate::complex_special::pi_c;
    use crate::complex_special::two_c;

    fn near(a: &ExactNum, b: &ExactNum, p: usize) -> bool {
        let d = a.sub(b, p, RoundingMode::None).abs();
        d.is_zero() || d.exponent().is_some_and(|e| e < -((p as i32) / 4))
    }

    fn cnear(a: &ExactComplex, b: &ExactComplex, p: usize) -> bool {
        near(a.re(), b.re(), p) && near(a.im(), b.im(), p)
    }

    fn cnear_bits(a: &ExactComplex, b: &ExactComplex, p: usize, slack: i32) -> bool {
        let dr = a.re().sub(b.re(), p, RoundingMode::None).abs();
        let di = a.im().sub(b.im(), p, RoundingMode::None).abs();
        (dr.is_zero() || dr.exponent().is_some_and(|e| e < -((p as i32) / slack)))
            && (di.is_zero() || di.exponent().is_some_and(|e| e < -((p as i32) / slack)))
    }

    fn tiny(x: &ExactNum, p: usize) -> bool {
        x.is_zero() || x.exponent().is_some_and(|e| e < -((p as i32) / 4))
    }

    #[test]
    fn test_complex_ei_golds() {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new().unwrap();

        let one = ExactComplex::one(p);
        let ei1 = one.ei(p, rm, &mut cc);
        let r1 = ExactNum::from_u8(1, p).ei(p, rm, &mut cc);
        assert!(near(ei1.re(), &r1, p));
        assert!(tiny(ei1.im(), p));

        let z0 = ExactComplex::zero(p);
        let s0 = z0.si(p, rm, &mut cc);
        assert!(tiny(s0.re(), p) && tiny(s0.im(), p));
        assert!(z0.ci(p, rm, &mut cc).is_nan());

        let z = ExactComplex::new(ExactNum::from_u8(1, p), half_c(p).re().clone());
        let sz = z.si(p, rm, &mut cc);
        assert!(cnear(&sz, &neg_c(&neg_c(&z).si(p, rm, &mut cc)), p));

        let fs0 = z0.fresnel_s(p, rm, &mut cc);
        let fc0 = z0.fresnel_c(p, rm, &mut cc);
        assert!(tiny(fs0.re(), p) && tiny(fs0.im(), p));
        assert!(tiny(fc0.re(), p) && tiny(fc0.im(), p));

        let h = ExactNum::from_u8(2, p).powsi(-((p as isize) / 8), p, rm);
        let hc = ExactComplex::from_real(h, p);
        let num_ei =
            z.add(&hc, p, rm)
                .ei(p, rm, &mut cc)
                .sub(&z.sub(&hc, p, rm).ei(p, rm, &mut cc), p, rm);
        let deriv_ei = num_ei.div(&hc.mul(&two_c(p), p, rm), p, rm);
        let expect_ei = z.exp(p, rm, &mut cc).div(&z, p, rm);
        assert!(cnear_bits(&deriv_ei, &expect_ei, p, 8));

        let num_si =
            z.add(&hc, p, rm)
                .si(p, rm, &mut cc)
                .sub(&z.sub(&hc, p, rm).si(p, rm, &mut cc), p, rm);
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

        let nan = ExactComplex::new(crate::NAN.clone(), ExactNum::new(p));
        assert!(nan.ei(p, rm, &mut cc).is_nan());
    }
}
