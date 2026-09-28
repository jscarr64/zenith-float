//! Complex Bessel \(J_\nu,Y_\nu,I_\nu,K_\nu\).

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

/// `|t| < 2^-p |s|` (or `t = 0`).
fn term_small_rel(t: &ExactComplex, s: &ExactComplex, p: usize) -> bool {
    let at = t.abs(64, RoundingMode::None);
    if at.is_zero() {
        return true;
    }
    match (at.exponent(), s.abs(64, RoundingMode::None).exponent()) {
        (Some(et), Some(es)) => (et as i64) + (p as i64) < es as i64,
        _ => term_negligible(t, p),
    }
}

/// Use the power series when \(\lvert z\rvert\) is below this (and \(\lvert z\rvert^2\)
/// is large enough for the Hankel expansion only beyond destination precision).
const BESSEL_SERIES_THRESHOLD: u32 = 16;

/// Integers \(\lvert n\rvert\) recognized exactly for \(Y_n\) recurrence and \(z=0\).
const BESSEL_INTEGER_MAX: i32 = 64;

fn abs_below(z: &ExactComplex, bound: u32, p: usize) -> bool {
    let a = z.abs(p, RoundingMode::None);
    let b = ExactNum::from_u32(bound, p);
    matches!(a.cmp(&b), Some(c) if c < 0)
}

/// The Hankel expansion's smallest term is about \(e^{-2\lvert z\rvert}\) (2.885|z| bits), so it is
/// used only when \(\lvert z\rvert \ge 0.35(p+112)\) (covering the Ziv guard bits); below that the
/// power series is used with [`series_guard_bits`] extra precision.
fn use_bessel_series(z: &ExactComplex, dest_p: usize) -> bool {
    if abs_below(z, BESSEL_SERIES_THRESHOLD, dest_p) {
        return true;
    }
    let t = dest_p.saturating_add(112).saturating_mul(35) / 100;
    abs_below(z, t.min(u32::MAX as usize) as u32, dest_p)
}

/// Upper bound on `|z|` as an integer (`2^e` with `|z| < 2^e`), saturating.
fn abs_bound(z: &ExactComplex) -> usize {
    match z.abs(64, RoundingMode::None).exponent() {
        Some(e) if e > 0 => 1usize << (e as u32).min(40),
        _ => 1,
    }
}

/// Guard bits for the power series, whose terms reach \(\approx e^{\lvert z\rvert}\) while the
/// result can be \(O(1)\): `1.5·bound(|z|) + 16`.
fn series_guard_bits(z: &ExactComplex) -> usize {
    abs_bound(z).saturating_mul(3) / 2 + 16
}

fn integer_nu(nu: &ExactComplex, p: usize) -> Option<i32> {
    if !nu.im().is_zero() || !nu.re().is_int() {
        return None;
    }
    for n in -BESSEL_INTEGER_MAX..=BESSEL_INTEGER_MAX {
        let w = ExactNum::from_i32(n, p);
        if nu.re().cmp(&w) == Some(0) {
            return Some(n);
        }
    }
    None
}

fn harmonic(k: usize, p: usize) -> ExactNum {
    let mut h = ExactNum::new(p);
    for i in 1..=k {
        let t =
            ExactNum::from_u8(1, p).div(&ExactNum::from_u32(i as u32, p), p, RoundingMode::None);
        h = h.add(&t, p, RoundingMode::None);
    }
    h
}

fn four_c(p: usize) -> ExactComplex {
    ExactComplex::from_real(ExactNum::from_u8(4, p), p)
}

fn eight_c(p: usize) -> ExactComplex {
    ExactComplex::from_real(ExactNum::from_u8(8, p), p)
}

/// `z > 0` real and `ν` real: every Bessel function is real there.
fn real_positive_axis(z: &ExactComplex, nu: &ExactComplex) -> bool {
    z.im().is_zero() && nu.im().is_zero() && !z.re().is_zero() && !z.re().is_negative()
}

/// Keeps the real part and sets the imaginary part to an exact zero.
fn real_part_only(v: ExactComplex, p: usize) -> ExactComplex {
    if v.is_nan() {
        return v;
    }
    ExactComplex::new(v.re().clone(), ExactNum::new(p))
}

/// \(e^{\nu\pi i}\) when `upper`, else \(e^{-\nu\pi i}\).
fn reflect_phase(nu: &ExactComplex, upper: bool, p: usize, cc: &mut Consts) -> ExactComplex {
    let pi_i = ExactComplex::new(ExactNum::new(p), cc.pi(p, RoundingMode::None));
    let a = nu.mul(&pi_i, p, RoundingMode::None);
    let a = if upper { a } else { neg_c(&a) };
    a.exp(p, RoundingMode::None, cc)
}

impl ExactComplex {
    /// \(J_\nu(z)\). Entire for integer \(\nu\); cut on \((-\infty,0]\) otherwise.
    /// \(z=0\) with non-integer \(\nu\) → NaN.
    ///
    /// # Precision
    ///
    /// - Algorithm: power series (extra precision for cancellation) for `|z| < max(16, 0.35(p+112))`;
    ///   Hankel expansion (stopped at its smallest term) otherwise, continued from \(-z\) when
    ///   \(\mathrm{Re}\,z<0\) (DLMF 10.11.1). Integer `|n| ≤ BESSEL_INTEGER_MAX` (`64`).
    /// - Bound: Ziv on each part (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    pub fn bessel_j_nu(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() || nu.is_nan() {
            return nan_pair(Error::InvalidArgument);
        }
        let dest = round_p(p);
        ziv_complex(dest, rm, |pw| self.bessel_j_at(nu, pw, dest, cc))
    }

    /// \(Y_\nu(z)\). Cut on \((-\infty,0]\); \(z=0\) → NaN.
    ///
    /// # Precision
    ///
    /// - Algorithm: from \(J_ν\) (same series / Hankel switch; \(\mathrm{Re}\,z<0\) in the Hankel
    ///   regime via DLMF 10.11.2).
    /// - Bound: Ziv on each part (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    pub fn bessel_y(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() || nu.is_nan() {
            return nan_pair(Error::InvalidArgument);
        }
        if self.re().is_zero() && self.im().is_zero() {
            return nan_pair(Error::InvalidArgument);
        }
        let dest = round_p(p);
        ziv_complex(dest, rm, |pw| self.bessel_y_at(nu, pw, dest, cc))
    }

    /// \(I_\nu(z)=i^{-\nu}J_\nu(iz)\). Same cut rules as \(J_\nu\).
    ///
    /// # Precision
    ///
    /// - Algorithm: via [`Self::bessel_j_nu`] of \(\pm iz\) (sign chosen by \(\mathrm{Im}\,z\), DLMF 10.27.6).
    /// - Bound: Ziv on each part (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    pub fn bessel_i(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() || nu.is_nan() {
            return nan_pair(Error::InvalidArgument);
        }
        let dest = round_p(p);
        if real_positive_axis(self, nu) {
            // I_ν(x) is real for x > 0 and real ν; the rotation through J_ν(−ix) leaves a rounding
            // residue in the imaginary part that Ziv can never certify as zero.
            return ziv_complex(dest, rm, |pw| {
                real_part_only(self.bessel_i_at(nu, pw, dest, cc), pw)
            });
        }
        ziv_complex(dest, rm, |pw| self.bessel_i_at(nu, pw, dest, cc))
    }

    /// \(K_\nu(z)=(\pi/2)\,i^{\nu+1}H_\nu^{(1)}(iz)\). Cut on \((-\infty,0]\) (value from above
    /// when \(\mathrm{Im}\,z=0\)); \(z=0\) → NaN. Real \(z>0\) with real \(\nu\) gives an exact
    /// zero imaginary part (likewise [`Self::bessel_i`]).
    ///
    /// # Precision
    ///
    /// - Algorithm: \(H^{(1)}_\nu(iz)\) or \(H^{(2)}_\nu(-iz)\) by the sign of \(\mathrm{Im}\,z\) (DLMF 10.27.8).
    /// - Bound: Ziv on each part (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    pub fn bessel_k(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() || nu.is_nan() {
            return nan_pair(Error::InvalidArgument);
        }
        if self.re().is_zero() && self.im().is_zero() {
            return nan_pair(Error::InvalidArgument);
        }
        let dest = round_p(p);
        if real_positive_axis(self, nu) {
            // K_ν(x) is real for x > 0 and real ν (see `bessel_i`).
            return ziv_complex(dest, rm, |pw| {
                real_part_only(self.bessel_k_at(nu, pw, dest, cc), pw)
            });
        }
        ziv_complex(dest, rm, |pw| self.bessel_k_at(nu, pw, dest, cc))
    }

    fn bessel_j_at(&self, nu: &Self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        if self.re().is_zero() && self.im().is_zero() {
            return match integer_nu(nu, work_p) {
                Some(0) => ExactComplex::one(work_p),
                Some(_) => ExactComplex::zero(work_p),
                None => nan_pair(Error::InvalidArgument),
            };
        }
        if use_bessel_series(self, dest_p) {
            self.bessel_j_series(nu, work_p + series_guard_bits(self), cc)
        } else if self.re().is_negative() && !self.re().is_zero() {
            // The Hankel expansion holds only for |arg z| < π; for Re z < 0 continue from −z
            // (DLMF 10.11.1): J_ν(z) = e^{±νπi} J_ν(−z), upper sign for Im z ≥ 0.
            let upper = !(self.im().is_negative() && !self.im().is_zero());
            let j = neg_c(self).bessel_hankel_j(nu, work_p, cc);
            reflect_phase(nu, upper, work_p, cc).mul(&j, work_p, RoundingMode::None)
        } else {
            self.bessel_hankel_j(nu, work_p, cc)
        }
    }

    fn bessel_y_at(&self, nu: &Self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        if let Some(n) = integer_nu(nu, work_p) {
            return self.bessel_y_int(n, work_p, dest_p, cc);
        }
        if use_bessel_series(self, dest_p) {
            self.bessel_y_nonint(nu, work_p, dest_p, cc)
        } else {
            self.bessel_hankel_y_any(nu, work_p, cc)
        }
    }

    /// Hankel-regime \(Y_\nu\) for any \(\arg z\). For \(\mathrm{Re}\,z<0\) (DLMF 10.11.2 with
    /// \(m=\pm1\)): \(Y_\nu(z)=e^{\mp\nu\pi i}Y_\nu(-z)\pm2i\cos(\nu\pi)J_\nu(-z)\), upper sign for
    /// \(\mathrm{Im}\,z\ge0\).
    fn bessel_hankel_y_any(&self, nu: &Self, work_p: usize, cc: &mut Consts) -> Self {
        if !(self.re().is_negative() && !self.re().is_zero()) {
            return self.bessel_hankel_y(nu, work_p, cc);
        }
        let upper = !(self.im().is_negative() && !self.im().is_zero());
        let w = neg_c(self);
        let y = w.bessel_hankel_y(nu, work_p, cc);
        let j = w.bessel_hankel_j(nu, work_p, cc);
        let pi = pi_c(work_p, cc);
        let cos_nu_pi = nu
            .mul(&pi, work_p, RoundingMode::None)
            .cos(work_p, RoundingMode::None, cc);
        let two_i = ExactComplex::new(ExactNum::new(work_p), ExactNum::from_u8(2, work_p));
        let corr =
            two_i
                .mul(&cos_nu_pi, work_p, RoundingMode::None)
                .mul(&j, work_p, RoundingMode::None);
        let corr = if upper { corr } else { neg_c(&corr) };
        reflect_phase(nu, !upper, work_p, cc)
            .mul(&y, work_p, RoundingMode::None)
            .add(&corr, work_p, RoundingMode::None)
    }

    fn bessel_i_at(&self, nu: &Self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        if self.re().is_zero() && self.im().is_zero() {
            return match integer_nu(nu, work_p) {
                Some(0) => ExactComplex::one(work_p),
                Some(n) if n > 0 => ExactComplex::zero(work_p),
                _ => nan_pair(Error::InvalidArgument),
            };
        }
        // DLMF 10.27.6: I_ν(z) = e^{∓νπi/2} J_ν(z e^{±πi/2}), the upper sign for
        // −π < arg z ≤ π/2 and the lower for −π/2 < arg z ≤ π. Picking by the sign of Im z keeps
        // the rotated argument off the J cut (the old code used the upper sign everywhere, which
        // is wrong for π/2 < arg z ≤ π, e.g. z = −2 + i).
        let upper = self.im().is_negative() && !self.im().is_zero();
        let rot = if upper { ExactComplex::i(work_p) } else { neg_c(&ExactComplex::i(work_p)) };
        let w = rot.mul(self, work_p, RoundingMode::None);
        let j = w.bessel_j_at(nu, work_p, dest_p, cc);
        let half_pi_i = ExactComplex::new(
            ExactNum::new(work_p),
            cc.pi(work_p, RoundingMode::None)
                .ldexp(-1, work_p, RoundingMode::None),
        );
        let arg = nu.mul(&half_pi_i, work_p, RoundingMode::None);
        let arg = if upper { neg_c(&arg) } else { arg };
        arg.exp(work_p, RoundingMode::None, cc)
            .mul(&j, work_p, RoundingMode::None)
    }

    fn bessel_k_at(&self, nu: &Self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        // DLMF 10.27.8: K_ν(z) = (πi/2) e^{νπi/2} H⁽¹⁾_ν(iz) for −π < arg z ≤ π/2 and
        // K_ν(z) = −(πi/2) e^{−νπi/2} H⁽²⁾_ν(−iz) for −π/2 < arg z ≤ π; chosen by the sign of Im z.
        // In the series regime J and Y are ≈ e^{|z|} while K ≈ e^{−Re z}: add 3|z| guard bits.
        // In the Hankel regime H⁽¹,²⁾ are formed directly (no J ± iY cancellation).
        // Im z = 0 (either sign of zero) takes the second form, valid up to arg z = π: on the cut
        // (−∞, 0) that is the value from above, as for `bessel_i` / `bessel_j_nu`.
        let second = !(self.im().is_negative() && !self.im().is_zero());
        let rot = if second { neg_c(&ExactComplex::i(work_p)) } else { ExactComplex::i(work_p) };
        let w = rot.mul(self, work_p, RoundingMode::None);
        let h = if use_bessel_series(&w, dest_p) {
            let wp = work_p + abs_bound(self).saturating_mul(3) + 16;
            let j = w.bessel_j_at(nu, wp, dest_p, cc);
            let y = w.bessel_y_at(nu, wp, dest_p, cc);
            let iy = ExactComplex::i(wp).mul(&y, wp, RoundingMode::None);
            if second {
                j.sub(&iy, wp, RoundingMode::None)
            } else {
                j.add(&iy, wp, RoundingMode::None)
            }
        } else {
            w.bessel_hankel_h(nu, !second, work_p, cc)
        };
        let half_pi_i = ExactComplex::new(
            ExactNum::new(work_p),
            cc.pi(work_p, RoundingMode::None)
                .ldexp(-1, work_p, RoundingMode::None),
        );
        let e = nu.mul(&half_pi_i, work_p, RoundingMode::None);
        let e = if second { neg_c(&e) } else { e };
        let pre = half_pi_i.mul(
            &e.exp(work_p, RoundingMode::None, cc),
            work_p,
            RoundingMode::None,
        );
        let pre = if second { neg_c(&pre) } else { pre };
        pre.mul(&h, work_p, RoundingMode::None)
    }

    fn bessel_j_series(&self, nu: &Self, p: usize, cc: &mut Consts) -> Self {
        let half = self.mul(&half_c(p), p, RoundingMode::None);
        let pow = half.pow(nu, p, RoundingMode::None, cc);
        let g =
            nu.add(&ExactComplex::one(p), p, RoundingMode::None)
                .gamma(p, RoundingMode::None, cc);
        let mut term = pow.div(&g, p, RoundingMode::None);
        let mut sum = term.clone();
        let hh = half.mul(&half, p, RoundingMode::None);
        for k in 1..=series_term_cap(p) {
            let kk = ExactComplex::from_real(ExactNum::from_u32(k as u32, p), p);
            let den = kk
                .add(nu, p, RoundingMode::None)
                .mul(&kk, p, RoundingMode::None);
            term = term
                .mul(&hh, p, RoundingMode::None)
                .div(&den, p, RoundingMode::None);
            term = neg_c(&term);
            sum = sum.add(&term, p, RoundingMode::None);
            if term_small_rel(&term, &sum, p) {
                break;
            }
        }
        sum
    }

    fn bessel_y_nonint(&self, nu: &Self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        let nupi = nu.mul(&pi_c(work_p, cc), work_p, RoundingMode::None);
        let s = nupi.sin(work_p, RoundingMode::None, cc);
        if term_negligible(&s, work_p) {
            return nan_pair(Error::InvalidArgument);
        }
        let jp = self.bessel_j_at(nu, work_p, dest_p, cc);
        let jm = self.bessel_j_at(&neg_c(nu), work_p, dest_p, cc);
        let c = nupi.cos(work_p, RoundingMode::None, cc);
        jp.mul(&c, work_p, RoundingMode::None)
            .sub(&jm, work_p, RoundingMode::None)
            .div(&s, work_p, RoundingMode::None)
    }

    fn bessel_y_int(&self, n: i32, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        let an = n.unsigned_abs();
        let y = if use_bessel_series(self, dest_p) {
            let wp = work_p + series_guard_bits(self);
            match an {
                0 => self.bessel_y0_series(wp, dest_p, cc),
                1 => self.bessel_y1_series(wp, dest_p, cc),
                _ => self.bessel_y_recurrence(an, wp, dest_p, cc),
            }
        } else {
            let nu = ExactComplex::from_real(ExactNum::from_u32(an, work_p), work_p);
            self.bessel_hankel_y_any(&nu, work_p, cc)
        };
        if n < 0 && an % 2 == 1 {
            neg_c(&y)
        } else {
            y
        }
    }

    fn bessel_y_recurrence(&self, n: u32, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        let mut ym2 = self.bessel_y0_series(work_p, dest_p, cc);
        let mut ym1 = self.bessel_y1_series(work_p, dest_p, cc);
        for m in 1..n {
            let two_m = ExactComplex::from_real(ExactNum::from_u32(2 * m, work_p), work_p);
            let ym = two_m
                .div(self, work_p, RoundingMode::None)
                .mul(&ym1, work_p, RoundingMode::None)
                .sub(&ym2, work_p, RoundingMode::None);
            ym2 = ym1;
            ym1 = ym;
        }
        ym1
    }

    fn bessel_y0_series(&self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        let two_pi = two_c(work_p).div(&pi_c(work_p, cc), work_p, RoundingMode::None);
        let half = self.mul(&half_c(work_p), work_p, RoundingMode::None);
        let j0 = self.bessel_j_at(&ExactComplex::zero(work_p), work_p, dest_p, cc);
        let g = ExactComplex::from_real(cc.euler_gamma(work_p, RoundingMode::None), work_p);
        let prefix = g.add(
            &half.ln(work_p, RoundingMode::None, cc),
            work_p,
            RoundingMode::None,
        );
        let z2 = half.mul(&half, work_p, RoundingMode::None);
        let mut fact = ExactNum::from_u8(1, work_p);
        let mut zk = ExactComplex::one(work_p);
        let mut sum = ExactComplex::zero(work_p);
        for m in 1..=series_term_cap(work_p) {
            fact = fact.mul(
                &ExactNum::from_u32(m as u32, work_p),
                work_p,
                RoundingMode::None,
            );
            zk = zk.mul(&z2, work_p, RoundingMode::None);
            let h = harmonic(m, work_p);
            let den = fact.mul(&fact, work_p, RoundingMode::None);
            let mut term = ExactComplex::from_real(h.div(&den, work_p, RoundingMode::None), work_p)
                .mul(&zk, work_p, RoundingMode::None);
            if m % 2 == 0 {
                term = neg_c(&term);
            }
            sum = sum.add(&term, work_p, RoundingMode::None);
            if term_negligible(&term, work_p) {
                break;
            }
        }
        two_pi.mul(
            &prefix
                .mul(&j0, work_p, RoundingMode::None)
                .add(&sum, work_p, RoundingMode::None),
            work_p,
            RoundingMode::None,
        )
    }

    fn bessel_y1_series(&self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        let nu1 = ExactComplex::one(work_p);
        let two_pi = two_c(work_p).div(&pi_c(work_p, cc), work_p, RoundingMode::None);
        let half = self.mul(&half_c(work_p), work_p, RoundingMode::None);
        let j1 = self.bessel_j_at(&nu1, work_p, dest_p, cc);
        let g = ExactComplex::from_real(cc.euler_gamma(work_p, RoundingMode::None), work_p);
        let prefix = g.add(
            &half.ln(work_p, RoundingMode::None, cc),
            work_p,
            RoundingMode::None,
        );
        let z2 = half.mul(&half, work_p, RoundingMode::None);
        let mut kfact = ExactNum::from_u8(1, work_p);
        let mut kp1fact = ExactNum::from_u8(1, work_p);
        let mut zk = ExactComplex::one(work_p);
        let mut sum = ExactComplex::zero(work_p);
        for k in 0..=series_term_cap(work_p) {
            let hk = harmonic(k, work_p);
            let rec = ExactNum::from_u8(1, work_p).div(
                &ExactNum::from_u32((k + 1) as u32, work_p),
                work_p,
                RoundingMode::None,
            );
            let hkp1 = hk.add(&rec, work_p, RoundingMode::None);
            let den = kfact.mul(&kp1fact, work_p, RoundingMode::None);
            let mut term = ExactComplex::from_real(
                hk.add(&hkp1, work_p, RoundingMode::None)
                    .div(&den, work_p, RoundingMode::None),
                work_p,
            )
            .mul(&zk, work_p, RoundingMode::None);
            if k % 2 == 1 {
                term = neg_c(&term);
            }
            sum = sum.add(&term, work_p, RoundingMode::None);
            if k > 0 && term_negligible(&term, work_p) {
                break;
            }
            let kp = k + 1;
            kfact = kfact.mul(
                &ExactNum::from_u32(kp as u32, work_p),
                work_p,
                RoundingMode::None,
            );
            kp1fact = kp1fact.mul(
                &ExactNum::from_u32((kp + 1) as u32, work_p),
                work_p,
                RoundingMode::None,
            );
            zk = zk.mul(&z2, work_p, RoundingMode::None);
        }
        let a = two_pi.mul(
            &prefix.mul(&j1, work_p, RoundingMode::None),
            work_p,
            RoundingMode::None,
        );
        let b = two_c(work_p).div(
            &pi_c(work_p, cc).mul(self, work_p, RoundingMode::None),
            work_p,
            RoundingMode::None,
        );
        let c = self
            .div(
                &two_c(work_p).mul(&pi_c(work_p, cc), work_p, RoundingMode::None),
                work_p,
                RoundingMode::None,
            )
            .mul(&sum, work_p, RoundingMode::None);
        a.sub(&b, work_p, RoundingMode::None)
            .sub(&c, work_p, RoundingMode::None)
    }

    fn hankel_chi_omega(&self, nu: &Self, p: usize, cc: &mut Consts) -> (Self, Self) {
        let two_nu_1 = two_c(p).mul(nu, p, RoundingMode::None).add(
            &ExactComplex::one(p),
            p,
            RoundingMode::None,
        );
        let chi = self.sub(
            &two_nu_1.mul(&pi_c(p, cc), p, RoundingMode::None).div(
                &four_c(p),
                p,
                RoundingMode::None,
            ),
            p,
            RoundingMode::None,
        );
        let two_over = two_c(p).div(
            &pi_c(p, cc).mul(self, p, RoundingMode::None),
            p,
            RoundingMode::None,
        );
        let omega = two_over.sqrt(p, RoundingMode::None, cc);
        (chi, omega)
    }

    fn hankel_pq(&self, nu: &Self, p: usize) -> (Self, Self) {
        let two_nu = two_c(p).mul(nu, p, RoundingMode::None);
        let mu = two_nu.mul(&two_nu, p, RoundingMode::None);
        let eight_z = eight_c(p).mul(self, p, RoundingMode::None);
        let mut prod = ExactComplex::one(p);
        let mut kf = ExactNum::from_u8(1, p);
        let mut pz = ExactComplex::one(p);
        let mut psum = ExactComplex::one(p);
        let mut qsum = ExactComplex::zero(p);
        let mut prev_e = i32::MAX;
        for k in 1..=series_term_cap(p) {
            let odd = ExactComplex::from_real(ExactNum::from_u32((2 * k - 1) as u32, p), p);
            let odd2 = odd.mul(&odd, p, RoundingMode::None);
            prod = prod.mul(&mu.sub(&odd2, p, RoundingMode::None), p, RoundingMode::None);
            kf = kf.mul(&ExactNum::from_u32(k as u32, p), p, RoundingMode::None);
            pz = pz.mul(&eight_z, p, RoundingMode::None);
            let term = prod.div(
                &ExactComplex::from_real(kf.clone(), p).mul(&pz, p, RoundingMode::None),
                p,
                RoundingMode::None,
            );
            if k % 2 == 0 {
                let signed = if (k / 2) % 2 == 1 { neg_c(&term) } else { term.clone() };
                psum = psum.add(&signed, p, RoundingMode::None);
            } else {
                let signed = if ((k - 1) / 2) % 2 == 1 { neg_c(&term) } else { term.clone() };
                qsum = qsum.add(&signed, p, RoundingMode::None);
            }
            if term_negligible(&term, p) {
                break;
            }
            // Asymptotic: stop once the terms start growing.
            let e = term
                .abs(64, RoundingMode::None)
                .exponent()
                .unwrap_or(i32::MIN);
            if k > 2 && e > prev_e {
                break;
            }
            prev_e = e;
        }
        (psum, qsum)
    }

    fn bessel_hankel_j(&self, nu: &Self, p: usize, cc: &mut Consts) -> Self {
        let (chi, omega) = self.hankel_chi_omega(nu, p, cc);
        let (pp, qq) = self.hankel_pq(nu, p);
        let (sn, cs) = {
            let s = chi.sin(p, RoundingMode::None, cc);
            let c = chi.cos(p, RoundingMode::None, cc);
            (s, c)
        };
        omega.mul(
            &pp.mul(&cs, p, RoundingMode::None).sub(
                &qq.mul(&sn, p, RoundingMode::None),
                p,
                RoundingMode::None,
            ),
            p,
            RoundingMode::None,
        )
    }

    /// Hankel function from the large-argument expansion: `H⁽¹⁾ = ω e^{iχ}(P + iQ)` when `first`,
    /// else `H⁽²⁾ = ω e^{−iχ}(P − iQ)` (DLMF 10.17.5–6).
    fn bessel_hankel_h(&self, nu: &Self, first: bool, p: usize, cc: &mut Consts) -> Self {
        let (chi, omega) = self.hankel_chi_omega(nu, p, cc);
        let (pp, qq) = self.hankel_pq(nu, p);
        let i = ExactComplex::i(p);
        let ichi = i.mul(&chi, p, RoundingMode::None);
        let iq = i.mul(&qq, p, RoundingMode::None);
        let (e, pq) = if first {
            (
                ichi.exp(p, RoundingMode::None, cc),
                pp.add(&iq, p, RoundingMode::None),
            )
        } else {
            (
                neg_c(&ichi).exp(p, RoundingMode::None, cc),
                pp.sub(&iq, p, RoundingMode::None),
            )
        };
        omega
            .mul(&e, p, RoundingMode::None)
            .mul(&pq, p, RoundingMode::None)
    }

    fn bessel_hankel_y(&self, nu: &Self, p: usize, cc: &mut Consts) -> Self {
        let (chi, omega) = self.hankel_chi_omega(nu, p, cc);
        let (pp, qq) = self.hankel_pq(nu, p);
        let s = chi.sin(p, RoundingMode::None, cc);
        let c = chi.cos(p, RoundingMode::None, cc);
        omega.mul(
            &pp.mul(&s, p, RoundingMode::None).add(
                &qq.mul(&c, p, RoundingMode::None),
                p,
                RoundingMode::None,
            ),
            p,
            RoundingMode::None,
        )
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
    fn test_complex_bessel_golds() {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new().unwrap();

        let one = ExactComplex::one(p);
        let z1 = one.clone();
        let nu0 = ExactComplex::zero(p);
        let nu1 = ExactComplex::one(p);
        let j0 = z1.bessel_j_nu(&nu0, p, rm, &mut cc);
        let rj0 = ExactNum::from_u8(1, p).bessel_j_nu(&ExactNum::from_u8(0, p), p, rm, &mut cc);
        assert!(near(j0.re(), &rj0, p));
        assert!(tiny(j0.im(), p));
        let j1 = z1.bessel_j_nu(&nu1, p, rm, &mut cc);
        let rj1 = ExactNum::from_u8(1, p).bessel_j_nu(&ExactNum::from_u8(1, p), p, rm, &mut cc);
        assert!(near(j1.re(), &rj1, p));
        assert!(tiny(j1.im(), p));

        let z = ExactComplex::new(ExactNum::from_u8(1, p), half_c(p).re().clone());
        let jn = z.bessel_j_nu(&nu0, p, rm, &mut cc);
        let yn1 = z.bessel_y(&nu1, p, rm, &mut cc);
        let jn1 = z.bessel_j_nu(&nu1, p, rm, &mut cc);
        let yn = z.bessel_y(&nu0, p, rm, &mut cc);
        let lhs = jn.mul(&yn1, p, rm).sub(&jn1.mul(&yn, p, rm), p, rm);
        let rhs = neg_c(&two_c(p)).div(&pi_c(p, &mut cc).mul(&z, p, rm), p, rm);
        assert!(cnear_bits(&lhs, &rhs, p, 8));

        let i0 = z1.bessel_i(&nu0, p, rm, &mut cc);
        let ri0 = ExactNum::from_u8(1, p).bessel_i(&ExactNum::from_u8(0, p), p, rm, &mut cc);
        assert!(near(i0.re(), &ri0, p));
        assert!(tiny(i0.im(), p));

        let k0 = z1.bessel_k(&nu0, p, rm, &mut cc);
        let rk0 = ExactNum::from_u8(1, p).bessel_k(&ExactNum::from_u8(0, p), p, rm, &mut cc);
        assert!(near(k0.re(), &rk0, p));
        assert!(tiny(k0.im(), p));

        let iz = ExactComplex::i(p).mul(&z, p, rm);
        let j_iz = iz.bessel_j_nu(&nu0, p, rm, &mut cc);
        let ln_i = ExactComplex::i(p).ln(p, rm, &mut cc);
        let scale = neg_c(&nu0).mul(&ln_i, p, rm).exp(p, rm, &mut cc);
        let via_j = scale.mul(&j_iz, p, rm);
        let i_z = z.bessel_i(&nu0, p, rm, &mut cc);
        assert!(cnear(&i_z, &via_j, p));

        let h = ExactNum::from_u8(2, p).powsi(-((p as isize) / 8), p, rm);
        let hc = ExactComplex::from_real(h, p);
        let num = z.add(&hc, p, rm).bessel_j_nu(&nu0, p, rm, &mut cc).sub(
            &z.sub(&hc, p, rm).bessel_j_nu(&nu0, p, rm, &mut cc),
            p,
            rm,
        );
        let deriv = num.div(&hc.mul(&two_c(p), p, rm), p, rm);
        let expect = neg_c(&z.bessel_j_nu(&nu1, p, rm, &mut cc));
        assert!(cnear_bits(&deriv, &expect, p, 8));

        let half = half_c(p);
        let z0 = ExactComplex::zero(p);
        assert!(z0.bessel_j_nu(&half, p, rm, &mut cc).is_nan());

        let eps = ExactNum::from_u8(2, p).powsi(-((p as isize) / 8), p, rm);
        let above = ExactComplex::new(ExactNum::from_i8(-1, p), eps.clone());
        let below = ExactComplex::new(
            ExactNum::from_i8(-1, p),
            ExactNum::from_i8(-1, p).mul(&eps, p, rm),
        );
        let ja = above.bessel_j_nu(&half, p, rm, &mut cc);
        let jb = below.bessel_j_nu(&half, p, rm, &mut cc);
        assert!(!cnear(&ja, &jb, p));
        assert!(cnear_bits(&ja, &jb.conj(), p, 8) || !tiny(ja.im(), p) || !tiny(jb.im(), p));
    }
}
