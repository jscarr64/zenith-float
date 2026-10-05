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
use crate::ops::special::BESSEL_GUARD_MAX;
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

/// Guard bits for the alternating \(J\)/\(Y\) power series, whose terms reach
/// \(\approx e^{\lvert z\rvert}\) while the result can be \(O(\lvert z\rvert^{-1/2})\):
/// `1.5·bound(|z|) + 16`.
///
/// Not used for \(I_\nu\). That series does not cancel against itself on the positive-real
/// axis; its guard is [`i_series_guard_bits`].
fn series_guard_bits(z: &ExactComplex) -> usize {
    abs_bound(z).saturating_mul(3) / 2 + 16
}

/// Fixed pad when the \(I\) series does not cancel. The 1.0.6 plan's 16–32 bit range
/// for a direct \({}_0F_1\) on \(\mathrm{Re}\,z>0\). 32 covers term-count rounding on
/// top of that range; it is not `1.5·|z|`.
const I_SERIES_PAD: usize = 32;

/// Real \(z>0\) and real \(\nu>-1\): every factor of \((z/2)^{2k}/(k!(\nu+1)_k)\) is positive,
/// so the running sum dominates each term.
fn i_terms_all_positive(z: &ExactComplex, nu: &ExactComplex) -> bool {
    if !z.im().is_zero() || !nu.im().is_zero() {
        return false;
    }
    if z.re().is_zero() || !z.re().is_positive() {
        return false;
    }
    let p = nu.re().precision().unwrap_or(64).max(64);
    let minus_one = ExactNum::from_i8(-1, p);
    matches!(nu.re().cmp(&minus_one), Some(c) if c > 0)
}

/// Upper bound on \(\lceil(|z|-\mathrm{Re}\,z)/\ln 2\rceil\).
///
/// A term of size \(e^{\lvert z\rvert}\) exceeds a result of size \(e^{\mathrm{Re}\,z}\) by that
/// many bits. It is 0 for real positive \(z\). It is only a lower estimate near a zero of
/// \(I_\nu\) (the sum can be arbitrarily smaller than \(e^{\mathrm{Re}\,z}\)); the series
/// measures peak term versus the sum and raises the guard when this bound is short.
///
/// `diff < 2^e` and \(\log_2 e < 36744/25469\), so the value is strictly above the real quotient.
fn i_series_exp_cancel_bits(z: &ExactComplex) -> usize {
    let az = z.abs(64, RoundingMode::None);
    let diff = az.sub(z.re(), 64, RoundingMode::None);
    if diff.is_zero() || !diff.is_positive() {
        return 0;
    }
    let Some(e) = diff.exponent() else {
        return 0;
    };
    if e <= 0 {
        return 2;
    }
    let e = (e as u32).min(40);
    let bound = ((1u128 << e) * 36744) / 25469 + 1;
    usize::try_from(bound).unwrap_or(usize::MAX / 4)
}

/// Guard for the \({}_0F_1\) series of \(I_\nu\) on \(\mathrm{Re}\,z>0\).
///
/// Where terms do not cancel (real \(z>0\), real \(\nu>-1\)), this is [`I_SERIES_PAD`]
/// only. Where they can, it is that pad plus the \((|z|-\mathrm{Re}\,z)/\ln 2\) cancellation
/// bound above — not a blanket `1.5·|z|` from "terms reach \(e^{\lvert z\rvert}\) while the
/// result can be \(O(1)\)".
fn i_series_guard_bits(z: &ExactComplex, nu: &ExactComplex) -> usize {
    if i_terms_all_positive(z, nu) {
        I_SERIES_PAD
    } else {
        i_series_exp_cancel_bits(z).saturating_add(I_SERIES_PAD)
    }
}

/// Enough terms to pass the peak near \(k\approx\lvert z\rvert/2\) and the Gaussian tail.
/// The cap must not track the guard: a `1.5·|z|` guard made this accidentally large, and a
/// small guard must still be allowed to run \(O(|z|)\) terms.
fn i_series_term_cap(z: &ExactComplex, dest_p: usize) -> usize {
    series_term_cap(dest_p)
        .saturating_add(abs_bound(z).saturating_mul(2))
        .min(u32::MAX as usize)
}

fn term_mag_exp(t: &ExactComplex) -> Option<i32> {
    let at = t.abs(64, RoundingMode::None);
    if at.is_zero() {
        None
    } else {
        at.exponent()
    }
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

fn cpx_tiny(z: &ExactComplex, p: usize) -> bool {
    let a = z.abs(64, RoundingMode::None);
    a.is_zero() || a.exponent().is_some_and(|e| (e as i64) + (p as i64) < 0)
}

fn term_near_one(delta: &ExactComplex, p: usize) -> bool {
    let one = ExactComplex::one(p.max(64));
    let d = delta.sub(&one, 64, RoundingMode::None);
    cpx_tiny(&d, p)
}

/// Modified Lentz for \({}_2F_0(a,b;;w)=1+\alpha_1/(1+\alpha_2/(1+\cdots))\),
/// \(\alpha_k=(a+k-1)(b+k-1)w/k\).
fn hypergeom_2f0_lentz(
    a: &ExactComplex,
    b: &ExactComplex,
    w: &ExactComplex,
    work_p: usize,
    dest_p: usize,
) -> Option<ExactComplex> {
    let one = ExactComplex::one(work_p);
    let tiny = ExactComplex::from_real(
        ExactNum::from_u8(2, work_p).powsi(-(work_p as isize), work_p, RoundingMode::None),
        work_p,
    );
    let mut f = one.clone();
    let mut c = f.clone();
    let mut d = ExactComplex::zero(work_p);
    let cap = series_term_cap(work_p).saturating_mul(4).max(64);
    for k in 1..=cap {
        let km1 = ExactComplex::from_real(ExactNum::from_u32((k - 1) as u32, work_p), work_p);
        let kk = ExactComplex::from_real(ExactNum::from_u32(k as u32, work_p), work_p);
        let ak = a
            .add(&km1, work_p, RoundingMode::None)
            .mul(
                &b.add(&km1, work_p, RoundingMode::None),
                work_p,
                RoundingMode::None,
            )
            .mul(w, work_p, RoundingMode::None)
            .div(&kk, work_p, RoundingMode::None);
        d = one.add(
            &ak.mul(&d, work_p, RoundingMode::None),
            work_p,
            RoundingMode::None,
        );
        if cpx_tiny(&d, dest_p) {
            d = tiny.clone();
        }
        c = one.add(
            &ak.div(&c, work_p, RoundingMode::None),
            work_p,
            RoundingMode::None,
        );
        if cpx_tiny(&c, dest_p) {
            c = tiny.clone();
        }
        d = one.div(&d, work_p, RoundingMode::None);
        let delta = c.mul(&d, work_p, RoundingMode::None);
        f = f.mul(&delta, work_p, RoundingMode::None);
        if term_near_one(&delta, dest_p) {
            return Some(f);
        }
    }
    None
}

fn eight_c(p: usize) -> ExactComplex {
    ExactComplex::from_real(ExactNum::from_u8(8, p), p)
}

fn c_mag_exp(z: &ExactComplex) -> Option<i32> {
    let a = z.abs(64, RoundingMode::None);
    if a.is_zero() {
        None
    } else {
        a.exponent()
    }
}

/// Upper bound on \(2\mathrm{Re}(z)/\ln 2\), the bit cancellation in \(I_{-\nu}-I_\nu\).
///
/// \(\mathrm{Re}\,z<2^e\) and \(2/\ln 2<73488/25469\), so the bound is strictly above the quotient.
fn re_pair_cancel_bits(z: &ExactComplex) -> usize {
    let re = z.re();
    if !re.is_positive() {
        return 32;
    }
    let Some(e) = re.exponent() else {
        return 32;
    };
    if e <= 0 {
        return 32;
    }
    if e > 20 {
        return BESSEL_GUARD_MAX;
    }
    let bound = ((1u128 << e) * 73488) / 25469 + 2;
    usize::try_from(bound).unwrap_or(BESSEL_GUARD_MAX)
}

/// Extra bits when \(\sin(\nu\pi)\) is small (near an integer order).
fn sin_nu_pi_loss(nu: &ExactComplex, cc: &mut Consts) -> Option<usize> {
    let p = 192;
    let mut n = nu.clone();
    n.set_precision(p, RoundingMode::None).ok()?;
    let s = n
        .mul(&pi_c(p, cc), p, RoundingMode::None)
        .sin(p, RoundingMode::None, cc);
    if s.is_nan() {
        return None;
    }
    let a = s.abs(64, RoundingMode::None);
    if a.is_zero() {
        return None;
    }
    let e = a.exponent()?;
    if e < -4096 {
        return None;
    }
    Some((8i64 - i64::from(e)).clamp(0, 4096) as usize)
}

fn i_diff_loss(ip: &ExactComplex, im: &ExactComplex, diff: &ExactComplex) -> usize {
    let ei = c_mag_exp(ip).into_iter().chain(c_mag_exp(im)).max();
    let Some(ei) = ei else {
        return 0;
    };
    let Some(ed) = c_mag_exp(diff) else {
        return BESSEL_GUARD_MAX;
    };
    (i64::from(ei) - i64::from(ed) + 8).max(0) as usize
}

/// `z > 0` real and `ν` real: every Bessel function is real there.
fn real_positive_axis(z: &ExactComplex, nu: &ExactComplex) -> bool {
    z.im().is_zero() && nu.im().is_zero() && !z.re().is_zero() && !z.re().is_negative()
}

/// Cut \((-\infty,0)\): \(\operatorname{Im} z=0\), \(\operatorname{Re} z<0\), real \(\nu\).
fn real_negative_axis(z: &ExactComplex, nu: &ExactComplex) -> bool {
    z.im().is_zero() && nu.im().is_zero() && z.re().is_negative()
}

/// Exact integer in \([-64,64]\).
fn real_small_int(x: &ExactNum, p: usize) -> Option<i32> {
    if !x.is_int() {
        return None;
    }
    (-BESSEL_INTEGER_MAX..=BESSEL_INTEGER_MAX)
        .find(|&n| x.cmp(&ExactNum::from_i32(n, p)) == Some(0))
}

/// \(\cos(\nu\pi),\sin(\nu\pi)\) with exact zeros for integer / half-odd-integer \(\nu\).
fn cos_sin_nu_pi(nu: &ExactNum, p: usize, cc: &mut Consts) -> (ExactNum, ExactNum) {
    if let Some(n) = real_small_int(nu, p) {
        let one = ExactNum::from_u8(1, p);
        return if n % 2 == 0 { (one, ExactNum::new(p)) } else { (one.neg(), ExactNum::new(p)) };
    }
    let two = ExactNum::from_u8(2, p);
    if let Some(t) = real_small_int(&two.mul(nu, p, RoundingMode::None), p) {
        if t % 2 != 0 {
            // ν = t/2 = n+1/2, n = (t-1)/2; sin(νπ)=(-1)^n, cos=0.
            let n = (t - 1) / 2;
            let s = ExactNum::from_u8(1, p);
            return (ExactNum::new(p), if n % 2 == 0 { s } else { s.neg() });
        }
    }
    let pi = cc.pi(p, RoundingMode::None);
    let a = nu.mul(&pi, p, RoundingMode::None);
    (
        a.cos(p, RoundingMode::None, cc),
        a.sin(p, RoundingMode::None, cc),
    )
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
    /// \(z=0\) with non-integer \(\nu\) → NaN. On the cut \(\operatorname{Im} z=0\),
    /// \(\operatorname{Re} z<0\) with real \(\nu\), the principal value is formed from the
    /// real kernels (\(J_\nu(-x+0i)=e^{i\nu\pi}J_\nu(x)\)) so algebraically zero parts stay zero.
    pub fn bessel_j_nu(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() || nu.is_nan() {
            return nan_pair(Error::InvalidArgument);
        }
        let dest = round_p(p);
        if real_negative_axis(self, nu) {
            return ziv_complex(dest, rm, |pw| self.bessel_neg_real_j(nu, pw, cc));
        }
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
        if real_negative_axis(self, nu) {
            return ziv_complex(dest, rm, |pw| self.bessel_neg_real_y(nu, pw, cc));
        }
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
        if real_negative_axis(self, nu) {
            return ziv_complex(dest, rm, |pw| self.bessel_neg_real_i(nu, pw, cc));
        }
        ziv_complex(dest, rm, |pw| self.bessel_i_at(nu, pw, dest, cc))
    }

    /// \(K_\nu(z)=(\pi/2)\,i^{\nu+1}H_\nu^{(1)}(iz)\). Cut on \((-\infty,0]\) (value from above
    /// when \(\mathrm{Im}\,z=0\)); \(z=0\) → NaN. Real \(z>0\) with real \(\nu\) gives an exact
    /// zero imaginary part (likewise [`Self::bessel_i`]). On the cut with real \(\nu\),
    /// \(K_\nu(-x+0i)=e^{-i\nu\pi}K_\nu(x)-i\pi I_\nu(x)\).
    ///
    /// # Precision
    ///
    /// - Algorithm: for non-integer \(\nu\) and \(\mathrm{Re}\,z>0\) inside the series regime,
    ///   DLMF 10.27.4, \(K_\nu=\pi(I_{-\nu}-I_\nu)/(2\sin\nu\pi)\), with the guard set by the
    ///   \(e^{2\mathrm{Re}\,z}\) cancellation. Temme \({}_2F_0\) (DLMF 10.32.10) when
    ///   \(\mathrm{Re}\,z>0\) and \(|z|\) is past the series/Hankel switch. Otherwise
    ///   \(H^{(1)}_\nu(iz)\) or \(H^{(2)}_\nu(-iz)\) by the sign of \(\mathrm{Im}\,z\)
    ///   (DLMF 10.27.8).
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
            // Past the series/Hankel switch the real kernel (asymptotic, or the I-connection
            // with a 2.885|x| guard) is used. Complex Temme 2F0 stays the fallback; its Lentz
            // step compares δ−1 at 64 bits, which is not enough for a 128/256-bit result.
            return ziv_complex(dest, rm, |pw| {
                if !use_bessel_series(self, dest) {
                    let k = self.re().bessel_k(nu.re(), pw, RoundingMode::None, cc);
                    if !k.is_nan() {
                        return ExactComplex::from_real(k, pw);
                    }
                }
                real_part_only(self.bessel_k_at(nu, pw, dest, cc), pw)
            });
        }
        if real_negative_axis(self, nu) {
            return ziv_complex(dest, rm, |pw| self.bessel_neg_real_k(nu, pw, cc));
        }
        ziv_complex(dest, rm, |pw| self.bessel_k_at(nu, pw, dest, cc))
    }

    /// Principal value on \(x<0\): \(J_\nu(-x+0i)=e^{i\nu\pi}J_\nu(x)\).
    fn bessel_neg_real_j(&self, nu: &Self, p: usize, cc: &mut Consts) -> Self {
        let x = self.re().abs();
        let nure = nu.re();
        let jx = x.bessel_j_nu(nure, p, RoundingMode::None, cc);
        let (c, s) = cos_sin_nu_pi(nure, p, cc);
        ExactComplex::new(
            c.mul(&jx, p, RoundingMode::None),
            s.mul(&jx, p, RoundingMode::None),
        )
    }

    /// \(Y_\nu(-x+0i)=e^{-i\nu\pi}Y_\nu(x)+2i\cos(\nu\pi)J_\nu(x)\).
    fn bessel_neg_real_y(&self, nu: &Self, p: usize, cc: &mut Consts) -> Self {
        let x = self.re().abs();
        let nure = nu.re();
        let jx = x.bessel_j_nu(nure, p, RoundingMode::None, cc);
        let yx = x.bessel_y(nure, p, RoundingMode::None, cc);
        let (c, s) = cos_sin_nu_pi(nure, p, cc);
        let two = ExactNum::from_u8(2, p);
        let re = c.mul(&yx, p, RoundingMode::None);
        let im = two
            .mul(&c, p, RoundingMode::None)
            .mul(&jx, p, RoundingMode::None)
            .sub(&s.mul(&yx, p, RoundingMode::None), p, RoundingMode::None);
        ExactComplex::new(re, im)
    }

    /// \(I_\nu(-x+0i)=e^{i\nu\pi}I_\nu(x)\).
    fn bessel_neg_real_i(&self, nu: &Self, p: usize, cc: &mut Consts) -> Self {
        let x = self.re().abs();
        let nure = nu.re();
        let ix = x.bessel_i(nure, p, RoundingMode::None, cc);
        let (c, s) = cos_sin_nu_pi(nure, p, cc);
        ExactComplex::new(
            c.mul(&ix, p, RoundingMode::None),
            s.mul(&ix, p, RoundingMode::None),
        )
    }

    /// \(K_\nu(-x+0i)=e^{-i\nu\pi}K_\nu(x)-i\pi I_\nu(x)\).
    fn bessel_neg_real_k(&self, nu: &Self, p: usize, cc: &mut Consts) -> Self {
        let x = self.re().abs();
        let nure = nu.re();
        let ix = x.bessel_i(nure, p, RoundingMode::None, cc);
        let kx = x.bessel_k(nure, p, RoundingMode::None, cc);
        let (c, s) = cos_sin_nu_pi(nure, p, cc);
        let pi = cc.pi(p, RoundingMode::None);
        let re = c.mul(&kx, p, RoundingMode::None);
        let im = pi
            .mul(&ix, p, RoundingMode::None)
            .add(&s.mul(&kx, p, RoundingMode::None), p, RoundingMode::None)
            .neg();
        ExactComplex::new(re, im)
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
            let upper = !self.im().is_negative() || self.im().is_zero();
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
        if !self.re().is_negative() || self.re().is_zero() {
            return self.bessel_hankel_y(nu, work_p, cc);
        }
        let upper = !self.im().is_negative() || self.im().is_zero();
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
        if !self.re().is_negative() && !self.re().is_zero() {
            return self.bessel_i_series(nu, work_p, dest_p, cc);
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

    /// \(I_\nu(z)=(z/2)^\nu/\Gamma(\nu+1)\,{}_0F_1(;\nu+1;z^2/4)\) for \(\mathrm{Re}\,z>0\).
    ///
    /// The guard is [`i_series_guard_bits`] (a 32-bit pad when every term is positive,
    /// otherwise pad plus cancellation). It is not [`series_guard_bits`]. If the summed
    /// peak still exceeds the result by more than that guard — \(I_\nu\) near a zero, where
    /// \((|z|-\mathrm{Re}\,z)/\ln 2\) is only a lower bound — the sum is repeated with the
    /// measured peak-versus-result gap.
    fn bessel_i_series(&self, nu: &Self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        let positive = i_terms_all_positive(self, nu);
        let mut guard = i_series_guard_bits(self, nu);
        let cap = i_series_term_cap(self, dest_p);
        if positive {
            return self
                .sum_bessel_i_series(nu, work_p.saturating_add(guard), dest_p, cap, cc)
                .0;
        }
        let guard_limit = guard
            .saturating_mul(4)
            .saturating_add(work_p)
            .saturating_add(64);
        for _ in 0..5 {
            let (sum, peak_exp, sum_exp) =
                self.sum_bessel_i_series(nu, work_p.saturating_add(guard), dest_p, cap, cc);
            if sum.is_nan() || guard >= guard_limit {
                return sum;
            }
            let measured = match (peak_exp, sum_exp) {
                (Some(pe), Some(se)) => (pe as i64 - se as i64 + 1).max(0) as usize,
                (Some(_), None) => guard.saturating_add(work_p),
                _ => 0,
            };
            let need = measured.saturating_add(I_SERIES_PAD);
            if need <= guard {
                return sum;
            }
            guard = need.min(guard_limit);
        }
        self.sum_bessel_i_series(nu, work_p.saturating_add(guard), dest_p, cap, cc)
            .0
    }

    /// One pass of the \({}_0F_1\) series at working precision `wp`.
    /// Returns the sum, the binary exponent of the largest term, and the exponent of the sum.
    fn sum_bessel_i_series(
        &self,
        nu: &Self,
        wp: usize,
        dest_p: usize,
        term_cap: usize,
        cc: &mut Consts,
    ) -> (Self, Option<i32>, Option<i32>) {
        let half = self.mul(&half_c(wp), wp, RoundingMode::None);
        let pow = half.pow(nu, wp, RoundingMode::None, cc);
        let g = nu
            .add(&ExactComplex::one(wp), wp, RoundingMode::None)
            .gamma(wp, RoundingMode::None, cc);
        let mut term = pow.div(&g, wp, RoundingMode::None);
        let mut sum = term.clone();
        let mut peak_exp = term_mag_exp(&term);
        let hh = half.mul(&half, wp, RoundingMode::None);
        for k in 1..=term_cap {
            let kk = ExactComplex::from_real(ExactNum::from_u32(k as u32, wp), wp);
            let den = kk
                .add(nu, wp, RoundingMode::None)
                .mul(&kk, wp, RoundingMode::None);
            term = term
                .mul(&hh, wp, RoundingMode::None)
                .div(&den, wp, RoundingMode::None);
            if let Some(e) = term_mag_exp(&term) {
                peak_exp = Some(peak_exp.map_or(e, |p| p.max(e)));
            }
            sum = sum.add(&term, wp, RoundingMode::None);
            if term_small_rel(&term, &sum, dest_p) {
                break;
            }
        }
        let sum_exp = term_mag_exp(&sum);
        (sum, peak_exp, sum_exp)
    }

    fn bessel_k_at(&self, nu: &Self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        if !self.re().is_negative() && !self.re().is_zero() {
            if use_bessel_series(self, dest_p) {
                // J±iY cancels ~3|z| bits here and builds both J_±ν. The I connection
                // is a direct 0F1; only the final subtraction needs the e^{2 Re z} guard.
                if let Some(k) = self.bessel_k_from_i(nu, work_p, cc) {
                    return k;
                }
            } else if let Some(k) = self.bessel_k_temme(nu, work_p, dest_p, cc) {
                return k;
            }
        }
        // DLMF 10.27.8: K_ν(z) = (πi/2) e^{νπi/2} H⁽¹⁾_ν(iz) for −π < arg z ≤ π/2 and
        // K_ν(z) = −(πi/2) e^{−νπi/2} H⁽²⁾_ν(−iz) for −π/2 < arg z ≤ π; chosen by the sign of Im z.
        // In the series regime J and Y are ≈ e^{|z|} while K ≈ e^{−Re z}: add 3|z| guard bits.
        // In the Hankel regime H⁽¹,²⁾ are formed directly (no J ± iY cancellation).
        // Im z = 0 (either sign of zero) takes the second form, valid up to arg z = π: on the cut
        // (−∞, 0) that is the value from above, as for `bessel_i` / `bessel_j_nu`.
        let second = !self.im().is_negative() || self.im().is_zero();
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

    /// DLMF 10.32.10: \(K_\nu(z)=\sqrt{\pi/(2z)}\,e^{-z}\,{}_2F_0(\tfrac12+\nu,\tfrac12-\nu;;-1/(2z))\).
    fn bessel_k_temme(
        &self,
        nu: &Self,
        work_p: usize,
        dest_p: usize,
        cc: &mut Consts,
    ) -> Option<Self> {
        let two = two_c(work_p);
        let half = half_c(work_p);
        let w = neg_c(&ExactComplex::one(work_p).div(
            &two.mul(self, work_p, RoundingMode::None),
            work_p,
            RoundingMode::None,
        ));
        let a = half.add(nu, work_p, RoundingMode::None);
        let b = half.sub(nu, work_p, RoundingMode::None);
        let f = hypergeom_2f0_lentz(&a, &b, &w, work_p, dest_p)?;
        let pi = pi_c(work_p, cc);
        let omega = pi
            .div(
                &two.mul(self, work_p, RoundingMode::None),
                work_p,
                RoundingMode::None,
            )
            .sqrt(work_p, RoundingMode::None, cc);
        let ez = neg_c(self).exp(work_p, RoundingMode::None, cc);
        Some(
            omega
                .mul(&ez, work_p, RoundingMode::None)
                .mul(&f, work_p, RoundingMode::None),
        )
    }

    /// DLMF 10.27.4. `None` for integer \(\nu\) (the \(Y\) recurrence is the stable form)
    /// or when \(\sin(\nu\pi)\) is too small to guard.
    fn bessel_k_from_i(&self, nu: &Self, work_p: usize, cc: &mut Consts) -> Option<Self> {
        if integer_nu(nu, work_p).is_some() {
            return None;
        }
        let sbits = sin_nu_pi_loss(nu, cc)?;
        let mut guard = re_pair_cancel_bits(self)
            .saturating_add(sbits)
            .saturating_add(48);
        if guard > BESSEL_GUARD_MAX {
            return None;
        }
        for _ in 0..3 {
            let acc = work_p.saturating_add(guard);
            if acc > work_p.saturating_add(BESSEL_GUARD_MAX) {
                return None;
            }
            let (k, loss) = self.k_i_diff(nu, acc, cc)?;
            if loss.saturating_add(32) <= guard {
                return Some(k);
            }
            guard = loss.saturating_add(64).max(guard.saturating_add(32));
            if guard > BESSEL_GUARD_MAX {
                return None;
            }
        }
        None
    }

    fn k_i_diff(&self, nu: &Self, acc: usize, cc: &mut Consts) -> Option<(Self, usize)> {
        let ip = self.bessel_i_series(nu, acc, acc, cc);
        let im = self.bessel_i_series(&neg_c(nu), acc, acc, cc);
        if ip.is_nan() || im.is_nan() {
            return None;
        }
        let diff = im.sub(&ip, acc, RoundingMode::None);
        if diff.is_nan() {
            return None;
        }
        let loss = i_diff_loss(&ip, &im, &diff);
        let pi = pi_c(acc, cc);
        let s = nu
            .mul(&pi, acc, RoundingMode::None)
            .sin(acc, RoundingMode::None, cc);
        if s.is_nan() || cpx_tiny(&s, acc.min(4096)) {
            return None;
        }
        let k = pi
            .div(&two_c(acc), acc, RoundingMode::None)
            .mul(&diff, acc, RoundingMode::None)
            .div(&s, acc, RoundingMode::None);
        if k.is_nan() {
            None
        } else {
            Some((k, loss))
        }
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

    #[test]
    fn test_complex_bessel_negative_real_axis() {
        let p = 128;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new().unwrap();
        let x = ExactNum::from_u8(2, p);
        let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, rm);
        let zneg = ExactComplex::from_real(x.neg(), p);
        let nu = ExactComplex::from_real(half.clone(), p);
        let i = zneg.bessel_i(&nu, p, rm, &mut cc);
        assert!(!i.is_nan());
        assert!(
            tiny(i.re(), p),
            "I_{{1/2}}(-2) real part should be exact 0, got {:?}",
            i.re()
        );
        let ix = x.bessel_i(&half, p, rm, &mut cc);
        assert!(near(i.im(), &ix, p));

        let k = zneg.bessel_k(&ExactComplex::zero(p), p, rm, &mut cc);
        let kx = x.bessel_k(&ExactNum::new(p), p, rm, &mut cc);
        let ix0 = x.bessel_i(&ExactNum::new(p), p, rm, &mut cc);
        let pi = cc.pi(p, rm);
        assert!(near(k.re(), &kx, p));
        assert!(near(k.im(), &pi.mul(&ix0, p, rm).neg(), p));

        let j = zneg.bessel_j_nu(&nu, p, rm, &mut cc);
        let jx = x.bessel_j_nu(&half, p, rm, &mut cc);
        assert!(tiny(j.re(), p));
        assert!(near(j.im(), &jx, p));
    }

    #[test]
    fn i_series_guard_sized_to_cancellation() {
        // Not an mpmath value. Locks the guard the series actually requests.
        // Restoring `1.5·bound(|z|)+16` on the positive-real axis is about
        // 49168 bits at |z|=20000 (`abs_bound` is 32768) and must fail here,
        // not merely make `complex_bessel_mpmath` slow.
        let p = 128;
        let z = ExactComplex::from_real(ExactNum::from_u32(20000, p), p);
        let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, RoundingMode::None);
        let nu = ExactComplex::from_real(half, p);
        let g = i_series_guard_bits(&z, &nu);
        assert_eq!(g, I_SERIES_PAD);
        assert!(
            (16..=32).contains(&g),
            "positive-real I guard {g} left the 16–32 bit pad"
        );
        let blanket = series_guard_bits(&z);
        assert!(
            blanket > 10_000,
            "test no longer distinguishes the old 1.5|z| guard ({blanket})"
        );
        assert!(g < blanket);

        // ν = -11/4 is outside the all-positive shortcut. z is real, so
        // (|z|−Re z)/ln 2 contributes 0 and the guard stays the pad.
        // Same source as above: the guard formula, not mpmath.
        let neg = ExactNum::from_i8(-11, p).div(&ExactNum::from_u8(4, p), p, RoundingMode::None);
        let nu_neg = ExactComplex::from_real(neg, p);
        assert_eq!(i_series_guard_bits(&z, &nu_neg), I_SERIES_PAD);

        // Complex z: cancellation bound, not the blanket 1.5|z| rule.
        // The expected width is `i_series_exp_cancel_bits` (a proven upper
        // bound on (|z|−Re z)/ln 2), not an mpmath digit string.
        let zc = ExactComplex::new(ExactNum::from_u8(1, p), ExactNum::from_u8(20, p));
        let cancel = i_series_exp_cancel_bits(&zc);
        assert!(cancel > 0, "1+20i should cancel");
        assert_eq!(i_series_guard_bits(&zc, &nu), cancel + I_SERIES_PAD);
        assert!(cancel + I_SERIES_PAD < series_guard_bits(&z));
    }
}
