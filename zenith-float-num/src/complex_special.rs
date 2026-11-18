//! Complex `erf` / `erfc` (Faddeeva) and `gamma` / `ln_gamma` / `digamma`.
//!
//! Software-limb [`ExactComplex`] arithmetic only. The error-function path is the
//! Faddeeva function throughout — not the real series evaluated at \(\lvert z\rvert\).

use crate::common::util::bump_prec_retry;
use crate::common::util::round_p;
use crate::complex_ball::ziv_complex_ball;
use crate::complex_ball::Mag;
use crate::complex_ball::CB;
use crate::complex_ball::RB;
use crate::defs::Exponent;
use crate::Consts;
use crate::Error;
use crate::ExactComplex;
use crate::ExactNum;
use crate::RoundingMode;
use crate::Word;
use crate::WORD_BIT_SIZE;
use alloc::vec::Vec;
use core::ops::Neg;

/// L1 radius \(\lvert\mathrm{Re}\,z\rvert+\lvert\mathrm{Im}\,z\rvert\) below which
/// the Faddeeva power series is used (Poppe–Wijers inner region).
const FADDEEVA_SERIES_L1: u32 = 8;

/// Minimum shift target: raise \(z\) by integers until \(\mathrm{exponent}(\lvert z\rvert)\)
/// is at least this (\(\lvert z\rvert\ge 2^{k-1}\)) so Stirling terms decay. The actual target
/// grows with the precision (see [`gamma_shift_exp`]): with a fixed target of 6 the 64-term series
/// capped \(\Gamma\) / \(\ln\Gamma\) at about 275 correct bits whatever `p` was.
#[allow(dead_code)]
const GAMMA_STIRLING_MIN_ABS_EXP: Exponent = 6;

/// Same minimum shift target for digamma (real kernel uses exponent \(< 8\)).
const DIGAMMA_STIRLING_MIN_ABS_EXP: Exponent = 8;

/// Shift target so the last kept Stirling term, \(\approx 2^{364}/\lvert z\rvert^{127}\), is below
/// \(2^{-(p+8)}\) (same bound as the real kernel).
#[allow(dead_code)] // kept for the shift-target unit tests
fn gamma_shift_exp(p: usize) -> Exponent {
    crate::ops::special::stirling_min_exponent(p, 372, 127, GAMMA_STIRLING_MIN_ABS_EXP)
}

/// Digamma analogue of [`gamma_shift_exp`]. Kept for the shift-target unit tests;
/// the runtime path uses [`digamma_stirling_plan`].
#[allow(dead_code)]
fn digamma_shift_exp(p: usize) -> Exponent {
    crate::ops::special::stirling_min_exponent(p, 379, 128, DIGAMMA_STIRLING_MIN_ABS_EXP)
}

/// Stirling term count and shift-target exponent at working precision `p`.
///
/// A fixed 64-term tail needs \(|z|\ge 2^{p/127}\), i.e. about a million
/// \(z\to z+1\) steps at `p = 2048`, which is why complex `gamma` / `digamma`
/// (and Bessel / Kelvin, which call them) hung there. The real kernel already
/// raises the term count with `p` (`digamma_terms`) so the target stays near
/// [`crate::ops::special::asymp_min_exponent`]. This is the same plan.
fn gamma_stirling_plan(p: usize) -> (usize, Exponent) {
    let terms = crate::ops::special::digamma_terms(p);
    let e_min = crate::ops::special::asymp_min_exponent(p + WORD_BIT_SIZE, terms);
    (terms, e_min)
}

/// Digamma plan: always `digamma_terms` / `asymp_min_exponent`, as in the real kernel.
fn digamma_stirling_plan(p: usize) -> (usize, Exponent) {
    let terms = crate::ops::special::digamma_terms(p);
    let e_min = crate::ops::special::asymp_min_exponent(p, terms);
    (terms, e_min)
}

// Positive integers n for which Γ(n)=(n-1)! is evaluated by multiplying 1…n-1
// instead of Stirling, when the bit-length bound of the real kernel holds
// ((n-1)·log2(n-1) ≤ 4p + 4096). Larger integers use Stirling under Ziv.

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

fn finite_parts(z: &ExactComplex) -> bool {
    !z.re().is_nan() && !z.im().is_nan() && !z.re().is_inf() && !z.im().is_inf()
}

pub(crate) fn term_negligible(t: &ExactComplex, p: usize) -> bool {
    let m = t.abs(p, RoundingMode::None);
    m.is_zero() || m.exponent().is_some_and(|e| e + (p as i64) < 0)
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

/// Most integer shifts \(z\to z+1\) the gamma-family loops take to reach
/// \(\lvert z\rvert\ge 2^{k-1}\) for the shift target `k` (`min_exp`). The loops start from
/// \(\mathrm{Re}\,z\ge 0\) and each shift adds 1 to \(\mathrm{Re}\,z\), so \(2^{k-1}\) shifts are
/// enough; 2 more cover rounding of \(\lvert z\rvert\). `usize::MAX` once \(2^{k-1}\) does not fit.
fn shift_cap(min_exp: Exponent) -> usize {
    crate::ops::special::stirling_shift_count(min_exp)
}

/// Applies `step` (in the callers: \(z\to z+1\) plus their accumulator update) while
/// [`abs_needs_shift`], at most `cap` times.
///
/// \(\lvert z+n\rvert\le\lvert z\rvert+n\), so at least \(2^{k-1}-\lvert z\rvert\) steps are
/// needed; when that is above `cap` it returns `PrecisionRetryExhausted` at once, without doing
/// the steps. Otherwise it stops with the same error after `cap` steps, so arithmetic that does not
/// move \(z\) up cannot loop forever. An `Err` from `step` is returned as is.
fn shift_while_small<S>(
    z: &mut ExactComplex,
    p: usize,
    min_exp: Exponent,
    cap: usize,
    mut step: S,
) -> Result<(), Error>
where
    S: FnMut(&mut ExactComplex) -> Result<(), Error>,
{
    let target = ExactNum::from_u8(1, p).ldexp(min_exp - 1, p, RoundingMode::None);
    let need = target.sub(&z.abs(p, RoundingMode::None), p, RoundingMode::None);
    let cap_num = ExactNum::from_word(Word::try_from(cap).unwrap_or(Word::MAX), p);
    if matches!(need.compare(&cap_num), Some(c) if c > 0) {
        return Err(Error::PrecisionRetryExhausted);
    }
    let mut left = cap;
    while abs_needs_shift(z, p, min_exp) {
        if left == 0 {
            return Err(Error::PrecisionRetryExhausted);
        }
        left -= 1;
        step(z)?;
    }
    Ok(())
}

#[allow(dead_code)]
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
fn even_bernoulli_numbers(kmax: usize, p: usize) -> Result<Vec<ExactNum>, Error> {
    crate::ops::special::with_even_bernoulli(kmax, p, |v| {
        v.iter()
            .map(|b| b.clone().map(ExactNum::from))
            .collect::<Result<Vec<_>, _>>()
    })
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

/// `(n-1)!` at a width that holds the bit-length bound used by `small_pos_int`.
fn factorial_um1_work(n: u32, p: usize) -> ExactComplex {
    let m = (n as usize).saturating_sub(1);
    let bitlen = (usize::BITS - m.leading_zeros()) as usize;
    let bits = m.saturating_mul(bitlen);
    let w = round_p(bits.max(p) + WORD_BIT_SIZE);
    factorial_um1(n, w)
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
        if !finite_parts(self) {
            return ziv_complex(dest, rm, |pw| self.erf_at_plain(pw, dest, cc));
        }
        // On the real axis the value is the real erf (correctly rounded, with its own
        // treatment of ±1 − δ). The imaginary part is a zero with the sign of the limit
        // Im erf(x + iε) = ε·(2/√π)e^{−x²}: the sign of Im z (erf(z̄) = conj erf(z)).
        if self.im().is_zero() {
            return ExactComplex::new(self.re().erf(dest, rm, cc), self.im().clone());
        }
        if let Some(v) = self.erf_near_unit(dest, rm, cc) {
            return v;
        }
        // The parts are only absolutely accurate (a tiny Im next to Re ≈ ±1, e^{-z²} phase,
        // asymptotic remainders): the ball loop carries an error bound for each part and
        // rounds a part only when the bound certifies it (ZF-68).
        ziv_complex_ball(dest, rm, |q| self.erf_ball(q, dest, cc))
    }

    /// ZF-37: far from the imaginary axis erf(z) = ±1 ∓ erfc(±z) with |erfc(±z)| below half an
    /// ulp of 1 (e.g. z = 100001 + i/2, |erfc| ≈ 2^-1.4e10). The real part is then `±1 + δ`
    /// for a nonzero δ of known sign, which no Ziv loop can decide; it is rounded once by
    /// `round_offset` from the sign of the (correctly rounded) Re erfc. The imaginary part is
    /// ∓Im erfc(±z), a normal value, under the ball loop. `None` everywhere else.
    fn erf_near_unit(&self, dest: usize, rm: RoundingMode, cc: &mut Consts) -> Option<Self> {
        if self.re().is_zero() || self.im().is_zero() || use_faddeeva_series(self, dest) {
            return None;
        }
        let pos = re_positive(self);
        let w = if pos { self.clone() } else { neg_c(self.clone()) };
        let e = w.erfc_word(dest, RoundingMode::ToEven, cc);
        let er = e.re();
        if er.is_nan() || er.is_zero() {
            return None;
        }
        let ex = er.exponent()?;
        if ex >= -(dest as Exponent) - 2 {
            return None;
        }
        // erf = 1 − erfc(z) (Re z > 0) or −1 + erfc(−z) (Re z < 0).
        let toward = if er.is_negative() == pos { crate::Sign::Pos } else { crate::Sign::Neg };
        let base = if pos { 1 } else { -1 };
        let re = match crate::num::ExactNumNumber::from_i8(base, WORD_BIT_SIZE)
            .and_then(|b| crate::num::ExactNumNumber::round_offset(&b, toward, dest, rm))
        {
            Ok(v) => ExactNum::from(v),
            Err(err) => return Some(nan_pair(err)),
        };
        let im = ziv_complex_ball(dest, rm, |q| {
            let v = self.erf_off_axis_ball(q, dest, cc)?;
            Some(CB::from_parts(
                RB::exact(ExactNum::from_u8(1, WORD_BIT_SIZE)),
                v.im,
            ))
        });
        Some(ExactComplex::new(re, im.im().clone()))
    }

    /// For `Re z < 0` past the series switch, `erfc(z) = 2 − erfc(−z)` with `|erfc(−z)|` below
    /// half an ulp of 2: the real part is `2 − δ` with δ of known sign, rounded once by
    /// `round_offset`; the imaginary part is `−Im erfc(−z)`, relatively accurate with no
    /// cancellation, rounded with the mirrored mode. Mirrors [`Self::erf_near_unit`]; without it the loop had to widen to
    /// |log2 δ| bits to see the real part (ZF-68; erfc(−76) took 60 s).
    fn erfc_near_two(&self, dest: usize, rm: RoundingMode, cc: &mut Consts) -> Option<Self> {
        if !self.re().is_negative() || self.re().is_zero() || use_faddeeva_series(self, dest) {
            return None;
        }
        let w = neg_c(self.clone());
        let er = w.erfc_word(dest, RoundingMode::ToEven, cc).re().clone();
        if er.is_nan() || er.is_zero() {
            return None;
        }
        if er.exponent()? >= -(dest as Exponent) - 2 {
            return None;
        }
        let toward = if er.is_negative() { crate::Sign::Pos } else { crate::Sign::Neg };
        let re = match crate::num::ExactNumNumber::from_i8(2, WORD_BIT_SIZE)
            .and_then(|b| crate::num::ExactNumNumber::round_offset(&b, toward, dest, rm))
        {
            Ok(v) => ExactNum::from(v),
            Err(err) => return Some(nan_pair(err)),
        };
        let mirrored = match rm {
            RoundingMode::Up => RoundingMode::Down,
            RoundingMode::Down => RoundingMode::Up,
            other => other,
        };
        let im = w.erfc_word(dest, mirrored, cc).im().clone().neg();
        Some(ExactComplex::new(re, im))
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
        if !finite_parts(self) {
            return ziv_complex(dest, rm, |pw| self.erfc_at_plain(pw, dest, cc));
        }
        // Real axis: the real erfc (correctly rounded, including 2 − δ). The imaginary part
        // is a zero signed by the limit Im erfc(x + iε) = −ε·(2/√π)e^{−x²}: opposite to Im z.
        if self.im().is_zero() {
            return ExactComplex::new(self.re().erfc(dest, rm, cc), self.im().clone().neg());
        }
        if let Some(v) = self.erfc_near_two(dest, rm, cc) {
            return v;
        }
        ziv_complex_ball(dest, rm, |q| self.erfc_ball(q, dest, cc))
    }

    /// Gamma \(\Gamma(z)\). Poles at non-positive integers → NaN.
    ///
    /// # Precision
    ///
    /// - Algorithm: Stirling (term count from `digamma_terms`, shift from `asymp_min_exponent`) plus reflection; factorial when `(n-1)!` fits the bit-length bound.
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
    /// - Algorithm: Stirling (term count from `digamma_terms`, shift from `asymp_min_exponent`) plus reflection.
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

    // One-shot forms kept for non-finite operands only (their NaN kinds are the contract there);
    // finite operands use the ball forms below.
    fn erf_at_plain(&self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        // On the real axis the imaginary part is exactly zero: the complex path only leaves rounding noise
        // there, which never rounds.
        if self.im().is_zero() {
            return ExactComplex::from_real(self.re().erf(work_p, RoundingMode::None, cc), work_p);
        }
        // On the imaginary axis erf(±0 + iy) = ±0 + i·erfi(y) (erf is odd and real on the real
        // axis): the real part is that exact signed zero, not the rounding noise of 1 − erfc.
        if self.re().is_zero() {
            let v = self.erf_off_axis_plain(work_p, dest_p, cc);
            let (_, im) = v.into_parts();
            return ExactComplex::new(self.re().clone(), im);
        }
        self.erf_off_axis_plain(work_p, dest_p, cc)
    }

    fn erf_off_axis_plain(&self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        if use_faddeeva_series(self, dest_p) {
            Self::erf_power_series_plain(self, work_p, cc)
        } else {
            ExactComplex::one(work_p).sub(
                &self.erfc_via_faddeeva_plain(work_p, dest_p, cc),
                work_p,
                RoundingMode::None,
            )
        }
    }

    fn erfc_at_plain(&self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        if self.im().is_zero() {
            return ExactComplex::from_real(self.re().erfc(work_p, RoundingMode::None, cc), work_p);
        }
        // erfc(±0 + iy) = 1 − i·erfi(y): the real part is exactly 1.
        if self.re().is_zero() {
            let (_, im) = self.erf_off_axis_plain(work_p, dest_p, cc).into_parts();
            return ExactComplex::new(ExactNum::from_u8(1, work_p), im.neg());
        }
        if use_faddeeva_series(self, dest_p) {
            ExactComplex::one(work_p).sub(
                &Self::erf_power_series_plain(self, work_p, cc),
                work_p,
                RoundingMode::None,
            )
        } else {
            self.erfc_via_faddeeva_plain(work_p, dest_p, cc)
        }
    }

    /// `erfc(z) = e^{-z²} w(iz)`, with the reflection form when `Re(z) ≤ 0`.
    fn erfc_via_faddeeva_plain(&self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        let z2 = self.mul(self, work_p, RoundingMode::None);
        let em = neg_c(z2).exp(work_p, RoundingMode::None, cc);
        let iz = ExactComplex::i(work_p).mul(self, work_p, RoundingMode::None);
        if re_positive(self) {
            em.mul(
                &Self::faddeeva_plain(&iz, work_p, dest_p, cc),
                work_p,
                RoundingMode::None,
            )
        } else {
            two_c(work_p).sub(
                &em.mul(
                    &Self::faddeeva_plain(&neg_c(iz), work_p, dest_p, cc),
                    work_p,
                    RoundingMode::None,
                ),
                work_p,
                RoundingMode::None,
            )
        }
    }

    /// Faddeeva \(w(z)=e^{-z^2}\mathrm{erfc}(-iz)\). Reduces to \(\mathrm{Im}\,z\ge 0\).
    fn faddeeva_plain(z: &Self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        if im_negative_or_neg_real(z) {
            let z2 = z.mul(z, work_p, RoundingMode::None);
            let two_exp = two_c(work_p).mul(
                &neg_c(z2).exp(work_p, RoundingMode::None, cc),
                work_p,
                RoundingMode::None,
            );
            return two_exp.sub(
                &Self::faddeeva_upper_plain(&z.neg_via_mul(), work_p, dest_p, cc),
                work_p,
                RoundingMode::None,
            );
        }
        Self::faddeeva_upper_plain(z, work_p, dest_p, cc)
    }

    fn faddeeva_upper_plain(z: &Self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        if use_faddeeva_series(z, dest_p) {
            Self::faddeeva_series_plain(z, work_p, cc)
        } else {
            Self::faddeeva_asymp_plain(z, work_p, cc)
        }
    }

    /// Poppe–Wijers power series: \(w(z)=e^{-z^2}(1-\mathrm{erf}(-iz))\).
    fn faddeeva_series_plain(z: &Self, p: usize, cc: &mut Consts) -> Self {
        let u = neg_c(ExactComplex::i(p)).mul(z, p, RoundingMode::None);
        let erf_u = Self::erf_power_series_plain(&u, p, cc);
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
    fn erf_power_series_plain(u: &Self, dest_p: usize, cc: &mut Consts) -> Self {
        let p = dest_p + erf_series_guard_bits(u);
        let pi = cc.pi(p, RoundingMode::None);
        let sqrt_pi = pi.sqrt(p, RoundingMode::None);
        let scale = ExactNum::from_u8(2, p).div(&sqrt_pi, p, RoundingMode::None);
        let u2 = u.mul(u, p, RoundingMode::None);
        let mut upow = u.clone();
        let mut nfact = ExactNum::from_u8(1, p);
        let mut sum = u.clone();
        let mut n: u32 = 1;
        loop {
            nfact = nfact.mul(&ExactNum::from_u32(n, p), p, RoundingMode::None);
            upow = upow.mul(&u2, p, RoundingMode::None);
            let two_n_1 = n
                .checked_mul(2)
                .and_then(|v| v.checked_add(1))
                .map(|v| ExactNum::from_u32(v, p));
            let Some(two_n_1) = two_n_1 else {
                return nan_pair(Error::PrecisionRetryExhausted);
            };
            let den = two_n_1.mul(&nfact, p, RoundingMode::None);
            let mut t = upow.div(&ExactComplex::from_real(den, p), p, RoundingMode::None);
            if n % 2 == 1 {
                t = neg_c(t);
            }
            sum = sum.add(&t, p, RoundingMode::None);
            if term_negligible(&t, p) {
                break;
            }
            let Some(next) = n.checked_add(1) else {
                return nan_pair(Error::PrecisionRetryExhausted);
            };
            n = next;
        }
        ExactComplex::from_real(scale, p).mul(&sum, dest_p, RoundingMode::None)
    }

    /// \(w(z)\sim i/(z\sqrt\pi)\sum(2m-1)!!/(2z^2)^m\) for \(\mathrm{Im}\,z\ge 0\).
    fn faddeeva_asymp_plain(z: &Self, p: usize, cc: &mut Consts) -> Self {
        let z2 = z.mul(z, p, RoundingMode::None);
        let two_z2 = two_c(p).mul(&z2, p, RoundingMode::None);
        let mut term = ExactComplex::one(p);
        let mut s = ExactComplex::one(p);
        let mut prev_e = Exponent::MIN;
        let mut m: u32 = 1;
        loop {
            let odd_n = m.checked_mul(2).and_then(|v| v.checked_sub(1));
            let Some(odd_n) = odd_n else {
                return nan_pair(Error::PrecisionRetryExhausted);
            };
            let odd = ExactComplex::from_real(ExactNum::from_u32(odd_n, p), p);
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
            let Some(next) = m.checked_add(1) else {
                return nan_pair(Error::PrecisionRetryExhausted);
            };
            m = next;
        }
        let sqrt_pi =
            ExactComplex::from_real(cc.pi(p, RoundingMode::None).sqrt(p, RoundingMode::None), p);
        let den = z.mul(&sqrt_pi, p, RoundingMode::None);
        ExactComplex::i(p)
            .div(&den, p, RoundingMode::None)
            .mul(&s, p, RoundingMode::None)
    }

    /// `erf` of a ball `w`: the library `erf` at the center (each part within 4 ulps at `q`)
    /// plus the propagated radius \(\max|\mathrm{erf}'|\cdot|\delta|\) with
    /// \(|\mathrm{erf}'(W)| = \tfrac{2}{\sqrt\pi}|e^{-W^2}| \le \tfrac{2}{\sqrt\pi}e^{-\mathrm{Re}\,w^2}
    /// e^{2|w||\delta|+|\delta|^2} < 4e^{-\mathrm{Re}\,w^2}\) for \(|w||\delta| \le 1/8\),
    /// \(|\delta| \le 1/4\). `None` when the ball is too wide for that bound.
    pub(crate) fn erf_of_ball(w: &CB, q: usize, cc: &mut Consts) -> Option<CB> {
        let c = ExactComplex::new(w.re.v.clone(), w.im.v.clone());
        let v = c.erf(q, RoundingMode::None, cc);
        let re = RB::rounded(v.re().clone(), q)?;
        let im = RB::rounded(v.im().clone(), q)?;
        let d = w.re.r.add(w.im.r);
        if d.is_zero() {
            return Some(CB::from_parts(re, im));
        }
        let de = d.ub_exp()?;
        let we = w.abs_ub().ub_exp().unwrap_or(Exponent::MIN / 2);
        if de > -2 || de.saturating_add(we) > -3 {
            return None;
        }
        // e^{-Re w²} = e^{im² − re²}, bounded from a 64-bit ball.
        let wr = RB::exact(c.re().clone());
        let wi = RB::exact(c.im().clone());
        let g = wi.sqr(64)?.sub(&wr.sqr(64)?, 64)?;
        let eg = g.exp(64, cc)?.abs_ub();
        Some(CB::from_parts(re, im).widen(eg.mul(d).mul_pow2(2)))
    }

    fn erf_ball(&self, q: usize, dest_p: usize, cc: &mut Consts) -> Option<CB> {
        // On the real axis the imaginary part is exactly zero (+0, as `from_real`).
        if self.im().is_zero() {
            return Some(CB::from_real(RB::rounded(
                self.re().erf(q, RoundingMode::None, cc),
                q,
            )?));
        }
        // On the imaginary axis erf(±0 + iy) = ±0 + i·erfi(y) (erf is odd and real on the real
        // axis): the real part is that exact signed zero.
        if self.re().is_zero() {
            let v = self.erf_off_axis_ball(q, dest_p, cc)?;
            return Some(CB::from_parts(RB::exact(self.re().clone()), v.im));
        }
        self.erf_off_axis_ball(q, dest_p, cc)
    }

    fn erf_off_axis_ball(&self, q: usize, dest_p: usize, cc: &mut Consts) -> Option<CB> {
        if use_faddeeva_series(self, dest_p) {
            Self::erf_series_ball(self, q, cc)
        } else {
            CB::int(1).sub(&self.erfc_via_faddeeva_ball(q, dest_p, cc)?, q)
        }
    }

    fn erfc_ball(&self, q: usize, dest_p: usize, cc: &mut Consts) -> Option<CB> {
        if self.im().is_zero() {
            return Some(CB::from_real(RB::rounded(
                self.re().erfc(q, RoundingMode::None, cc),
                q,
            )?));
        }
        // erfc(±0 + iy) = 1 − i·erfi(y): the real part is exactly 1.
        if self.re().is_zero() {
            let v = self.erf_off_axis_ball(q, dest_p, cc)?;
            return Some(CB::from_parts(
                RB::exact(ExactNum::from_u8(1, WORD_BIT_SIZE)),
                v.im.neg(),
            ));
        }
        if use_faddeeva_series(self, dest_p) {
            CB::int(1).sub(&Self::erf_series_ball(self, q, cc)?, q)
        } else {
            self.erfc_via_faddeeva_ball(q, dest_p, cc)
        }
    }

    /// `erfc(z) = e^{-z²} w(iz)`, with the reflection form when `Re(z) ≤ 0`.
    fn erfc_via_faddeeva_ball(&self, q: usize, dest_p: usize, cc: &mut Consts) -> Option<CB> {
        let z = CB::exact(self);
        let em = z.mul(&z, q)?.neg().exp(q, cc)?;
        let iz = ExactComplex::new(self.im().clone().neg(), self.re().clone());
        if re_positive(self) {
            em.mul(&Self::faddeeva_ball(&iz, q, dest_p, cc)?, q)
        } else {
            let w = Self::faddeeva_ball(&neg_c(iz), q, dest_p, cc)?;
            CB::int(2).sub(&em.mul(&w, q)?, q)
        }
    }

    /// Faddeeva \(w(z)=e^{-z^2}\mathrm{erfc}(-iz)\). Reduces to \(\mathrm{Im}\,z\ge 0\) by
    /// \(w(z) = 2e^{-z^2} - w(-z)\). `z` is exact.
    fn faddeeva_ball(z: &Self, q: usize, dest_p: usize, cc: &mut Consts) -> Option<CB> {
        if im_negative_or_neg_real(z) {
            let zb = CB::exact(z);
            let two_exp = zb.mul(&zb, q)?.neg().exp(q, cc)?.scale2(1)?;
            return two_exp.sub(
                &Self::faddeeva_upper_ball(&neg_c(z.clone()), q, dest_p, cc)?,
                q,
            );
        }
        Self::faddeeva_upper_ball(z, q, dest_p, cc)
    }

    fn faddeeva_upper_ball(z: &Self, q: usize, dest_p: usize, cc: &mut Consts) -> Option<CB> {
        if !use_faddeeva_series(z, dest_p) {
            if let Some(w) = Self::faddeeva_asymp_ball(z, q, cc) {
                return Some(w);
            }
        }
        Self::faddeeva_series_ball(z, q, cc)
    }

    /// Poppe–Wijers power series: \(w(z)=e^{-z^2}(1-\mathrm{erf}(-iz))\).
    fn faddeeva_series_ball(z: &Self, q: usize, cc: &mut Consts) -> Option<CB> {
        let u = ExactComplex::new(z.im().clone(), z.re().clone().neg());
        let erfc_u = CB::int(1).sub(&Self::erf_series_ball(&u, q, cc)?, q)?;
        let zb = CB::exact(z);
        zb.mul(&zb, q)?.neg().exp(q, cc)?.mul(&erfc_u, q)
    }

    /// Entire power series \(\mathrm{erf}(u)=\frac{2}{\sqrt\pi}\sum a_n\),
    /// \(a_n=(-1)^n u^{2n+1}/(n!(2n+1))\), `u` exact, as a ball.
    ///
    /// Terms reach \(\approx e^{\lvert u\rvert^2}\) while the sum can be \(O(1)\), so the series
    /// runs with `1.5·bound(|u|²) + 16` guard bits (\(\log_2 e<1.5\)); the radius records the
    /// rounding of every term. It stops at a term below \(2^{-p}\) once \(n+1 \ge 2|u|^2\):
    /// from there \(|a_{k+1}/a_k| \le |u|^2/(k+1) \le 1/2\), so the tail is at most
    /// \(|a_n|\), which is added to the radius. A sum that has not stopped within the term
    /// limit is `None` (never a truncation).
    fn erf_series_ball(u: &Self, dest_p: usize, cc: &mut Consts) -> Option<CB> {
        let guard = erf_series_guard_bits(u);
        let p = dest_p + guard;
        let pi = RB::rounded(cc.pi(p, RoundingMode::None), p)?;
        let scale = RB::exact(ExactNum::from_u8(2, WORD_BIT_SIZE)).div(&pi.sqrt(p)?, p)?;
        let ub = CB::exact(u);
        let mu2 = ub.mul(&ub, p)?.neg();
        let u_abs = ub.abs_ub();
        let two_u2 = u_abs.mul(u_abs).mul_pow2(1);
        let mut term = ub.clone();
        let mut sum = ub;
        let mut n: usize = 1;
        loop {
            let n32 = u32::try_from(n).ok()?;
            let k2 = n32.checked_mul(2)?.checked_add(1)?;
            term = term
                .mul(&mu2, p)?
                .div_real(&RB::exact(ExactNum::from_u32(n32, WORD_BIT_SIZE)), p)?;
            let a = term.div_real(&RB::exact(ExactNum::from_u32(k2, WORD_BIT_SIZE)), p)?;
            sum = sum.add(&a, p)?;
            let am = a.abs_ub();
            let small = am
                .ub_exp()
                .is_none_or(|e| e.saturating_add(p as Exponent) < 0);
            // 2|u|² ≤ n + 1: 2|u|² < 2^k ≤ n + 1 with k = floor(log2(n + 1)).
            let k = Exponent::from(63 - ((n as u64) + 1).leading_zeros());
            let ratio_ok = two_u2.ub_exp().is_none_or(|e| e <= k);
            if small && ratio_ok {
                return sum.widen(am).mul_real(&scale, p);
            }
            n = n.checked_add(1)?;
        }
    }

    /// \(w(z)\sim i/(z\sqrt\pi)\sum_m(2m-1)!!/(2z^2)^m\) for \(\mathrm{Im}\,z > 0\), `z` exact,
    /// as a ball. This is DLMF 7.12.1 for \(\mathrm{erfc}\,\zeta\), \(\zeta=-iz\); the remainder
    /// after the last summed term is bounded by the first neglected term when
    /// \(|\mathrm{ph}\,\zeta|\le\pi/4\) (\(\mathrm{Im}\,z\ge|\mathrm{Re}\,z|\)) and by
    /// \(\csc(2|\mathrm{ph}\,\zeta|)=|z|^2/(2|\mathrm{Re}\,z|\,\mathrm{Im}\,z)\) times it for
    /// \(\pi/4\le|\mathrm{ph}\,\zeta|<\pi/2\) (DLMF 7.12(i)); that bound is added to the radius.
    /// `None` on the real axis (no bound), when the terms grow before one falls below
    /// \(2^{-q}\) (the expansion's smallest term is about \(e^{-|z|^2}\); ZF-68), or when the
    /// term limit is reached.
    fn faddeeva_asymp_ball(z: &Self, q: usize, cc: &mut Consts) -> Option<CB> {
        if z.im().is_zero() || z.im().is_negative() {
            return None;
        }
        let zb = CB::exact(z);
        let two_z2 = zb.mul(&zb, q)?.scale2(1)?;
        let lo2 = match (two_z2.re.lo_exp(), two_z2.im.lo_exp()) {
            (Some(a), Some(b)) => a.max(b),
            (Some(a), None) | (None, Some(a)) => a,
            (None, None) => return None,
        };
        let mut term = CB::int(1);
        let mut s = CB::int(1);
        let mut prev_e: Option<Exponent> = None;
        let mut m: usize = 1;
        loop {
            let odd = u32::try_from(m.checked_mul(2)?.checked_sub(1)?).ok()?;
            term = term
                .mul_real(&RB::exact(ExactNum::from_u32(odd, WORD_BIT_SIZE)), q)?
                .div(&two_z2, q)?;
            s = s.add(&term, q)?;
            let tm = term.abs_ub();
            let e = tm.ub_exp()?;
            if e.saturating_add(q as Exponent) < 0 {
                // First neglected term: |term|·(2m+1)/|2z²|.
                let next = tm.mul_u64(2 * m as u64 + 1).mul_pow2(-lo2);
                let k = Self::asymp_csc_bound(z)?;
                s = s.widen(next.mul(k));
                let sqrt_pi = RB::rounded(cc.pi(q, RoundingMode::None), q)?.sqrt(q)?;
                let den = zb.mul_real(&sqrt_pi, q)?;
                return CB::i().div(&den, q)?.mul(&s, q);
            }
            if prev_e.is_some_and(|pe| e > pe) {
                return None;
            }
            prev_e = Some(e);
            m = m.checked_add(1)?;
        }
    }

    /// Upper bound of the DLMF 7.12(i) remainder factor for `z` (`Im z > 0`): 1 when
    /// `Im z ≥ |Re z|`, else `csc(2|ph ζ|) = |z|²/(2|Re z| Im z) < 2^(2e_max − e_re − e_im + 2)`.
    fn asymp_csc_bound(z: &Self) -> Option<Mag> {
        if z.re().is_zero() || matches!(z.im().abs_cmp(z.re()), Some(c) if c >= 0) {
            return Some(Mag::pow2(0));
        }
        let er = z.re().exponent()?;
        let ei = z.im().exponent()?;
        let emax = er.max(ei);
        Some(Mag::pow2(2 * emax - er - ei + 2))
    }

    pub(crate) fn gamma_at(&self, p: usize, cc: &mut Consts) -> Self {
        if let Some(n) = self.small_pos_int(p) {
            return factorial_um1_work(n, p);
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
        let n = self.re().to_u32(RoundingMode::ToZero).ok()?;
        let m = (n as usize).saturating_sub(1);
        let bitlen = (usize::BITS - m.leading_zeros()) as usize;
        let bits = m.saturating_mul(bitlen);
        if bits <= p.saturating_mul(4).saturating_add(4096) {
            Some(n)
        } else {
            None
        }
    }

    fn gamma_positive(&self, p: usize, cc: &mut Consts) -> Self {
        let one = ExactComplex::one(p);
        let mut z = self.clone();
        let mut acc = ExactComplex::one(p);
        let (terms, k) = gamma_stirling_plan(p);
        let shifted = shift_while_small(&mut z, p, k, shift_cap(k), |z| {
            acc = acc.mul(z, p, RoundingMode::None);
            *z = ExactComplex::add_ref(z, &one, p, RoundingMode::None);
            Ok(())
        });
        if let Err(e) = shifted {
            return nan_pair(e);
        }
        let lg = z.ln_gamma_stirling_with(p, terms, cc);
        let g = lg.exp(p, RoundingMode::None, cc);
        g.div(&acc, p, RoundingMode::None)
    }

    fn ln_gamma_at(&self, p: usize, cc: &mut Consts) -> Self {
        if let Some(n) = self.small_pos_int(p) {
            return factorial_um1_work(n, p).ln(p, RoundingMode::None, cc);
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
        let (terms, k) = gamma_stirling_plan(p);
        let shifted = shift_while_small(&mut z, p, k, shift_cap(k), |z| {
            ln_acc = ExactComplex::add_ref(
                &ln_acc,
                &z.ln(p, RoundingMode::None, cc),
                p,
                RoundingMode::None,
            );
            *z = ExactComplex::add_ref(z, &one, p, RoundingMode::None);
            Ok(())
        });
        if let Err(e) = shifted {
            return nan_pair(e);
        }
        z.ln_gamma_stirling_with(p, terms, cc)
            .sub(&ln_acc, p, RoundingMode::None)
    }

    fn ln_gamma_stirling_with(&self, p: usize, terms: usize, cc: &mut Consts) -> Self {
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
        let bs = match even_bernoulli_numbers(terms, p) {
            Ok(v) => v,
            Err(e) => return nan_pair(e),
        };
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
        // Real-axis closed forms, same as the real kernel. The Stirling shift
        // would otherwise walk from 1 to 2^{e_min−1} (≈ 2^{14} steps at p=65536
        // with the bound-based term count, or 2^{40} with a 1024-term clamp).
        if self.im().is_zero() && self.re().is_positive() {
            // ψ(n) = −γ + H_{n−1} for positive integers.
            if self.re().is_int() {
                if let Ok(n) = self.re().to_u32(RoundingMode::ToZero) {
                    let g = cc.euler_gamma(p, RoundingMode::None);
                    let mut v = -g;
                    if n >= 2 {
                        let mut h = ExactNum::new(p);
                        for k in 1..n {
                            let t = ExactNum::from_u8(1, p).div(
                                &ExactNum::from_u32(k, p),
                                p,
                                RoundingMode::None,
                            );
                            h = h.add(&t, p, RoundingMode::None);
                        }
                        v = v.add(&h, p, RoundingMode::None);
                    }
                    return ExactComplex::from_real(v, p);
                }
            } else {
                // ψ(n+½) = −γ − 2 ln 2 + Σ_{k=1}^{n} 2/(2k−1) (DLMF 5.4.13–14).
                let twice = self.re().clone().ldexp(1, p, RoundingMode::None);
                if twice.is_int() {
                    if let Ok(m) = twice.to_u32(RoundingMode::ToZero) {
                        if m % 2 == 1 {
                            let n = ((m - 1) / 2) as usize;
                            // Same finite-sum cutoff as the real kernel
                            // (n ≤ 4p + 4096); larger half-integers use Stirling.
                            if n <= p.saturating_mul(4).saturating_add(4096) {
                                let g = cc.euler_gamma(p, RoundingMode::None);
                                let ln2 = cc.ln_2(p, RoundingMode::None);
                                let two_ln2 = ln2.ldexp(1, p, RoundingMode::None);
                                let mut v = g.neg().sub(&two_ln2, p, RoundingMode::None);
                                let two = ExactNum::from_u8(2, p);
                                for k in 1..=n {
                                    let t = two.clone().div(
                                        &ExactNum::from_u32((2 * k - 1) as u32, p),
                                        p,
                                        RoundingMode::None,
                                    );
                                    v = v.add(&t, p, RoundingMode::None);
                                }
                                return ExactComplex::from_real(v, p);
                            }
                        }
                    }
                }
            }
        }
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
        let (terms, k) = digamma_stirling_plan(p);
        // Term count and shift target come from the remainder bound
        // (`digamma_terms` / `asymp_min_exponent`). shift_cap(k) is the
        // step count needed to reach 2^{k-1} from Re z ≥ 0, not a
        // separate walk cap.
        let shifted = shift_while_small(&mut z, p, k, shift_cap(k), |z| {
            let rec = one.div(z, p, RoundingMode::None);
            acc = ExactComplex::sub_ref(&acc, &rec, p, RoundingMode::None);
            *z = ExactComplex::add_ref(z, &one, p, RoundingMode::None);
            Ok(())
        });
        if let Err(e) = shifted {
            return nan_pair(e);
        }
        acc.add(&z.digamma_asymp(p, terms, cc), p, RoundingMode::None)
    }

    fn digamma_asymp(&self, p: usize, terms: usize, cc: &mut Consts) -> Self {
        let ln_z = self.ln(p, RoundingMode::None, cc);
        let two_z = two_c(p).mul(self, p, RoundingMode::None);
        let half_inv = ExactComplex::one(p).div(&two_z, p, RoundingMode::None);
        let mut s = ln_z.sub(&half_inv, p, RoundingMode::None);
        let z2 = self.mul(self, p, RoundingMode::None);
        let mut zp = ExactComplex::one(p);
        let bs = match even_bernoulli_numbers(terms, p) {
            Ok(v) => v,
            Err(e) => return nan_pair(e),
        };
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
        let mut n: u32 = 1;
        loop {
            nfact = nfact.mul(&ExactNum::from_u32(n, p), p, RoundingMode::None);
            xpow = xpow.mul(&x2, p, RoundingMode::None);
            let Some(odd) = n.checked_mul(2).and_then(|v| v.checked_add(1)) else {
                break;
            };
            let den = ExactNum::mul_ref(&nfact, &ExactNum::from_u32(odd, p), p, RoundingMode::None);
            let t = ExactNum::div_ref(&xpow, &den, p, RoundingMode::None);
            sum = sum.add(&t, p, RoundingMode::None);
            if t.is_zero() || t.exponent().is_some_and(|e| e + (p as i64) < 0) {
                break;
            }
            let Some(next) = n.checked_add(1) else {
                break;
            };
            n = next;
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

    /// `shift_cap(k)` is 2^(k-1) + 2, with k clamped to at least 1, and `usize::MAX` once 2^(k-1)
    /// does not fit a `usize`. No other limit: from dest p 1984 (k = 22) the cap is past 2^20.
    #[test]
    fn shift_cap_is_two_to_the_k_minus_one_plus_two() -> Result<(), crate::Error> {
        assert_eq!(shift_cap(6), 34);
        assert_eq!(shift_cap(8), 130);
        assert_eq!(shift_cap(2), 4);
        assert_eq!(shift_cap(1), 3);
        assert_eq!(shift_cap(0), 3);
        assert_eq!(shift_cap(-5), 3);
        assert_eq!(shift_cap(20), (1usize << 19) + 2);
        assert_eq!(shift_cap(21), (1usize << 20) + 2);
        assert_eq!(shift_cap(22), (1usize << 21) + 2);
        assert_eq!(shift_cap(23), (1usize << 22) + 2);
        let bits = Exponent::from(usize::BITS);
        assert_eq!(shift_cap(bits), (1usize << (usize::BITS - 1)) + 2);
        assert_eq!(shift_cap(bits + 1), usize::MAX);
        assert_eq!(shift_cap(Exponent::MAX), usize::MAX);
        assert_eq!(shift_cap(gamma_shift_exp(256)), 34);
        assert_eq!(shift_cap(digamma_shift_exp(256)), 130);
        Ok(())
    }

    /// Runs `shift_while_small` from `z` with the real z -> z+1 step (or a step that does not move
    /// z when `stuck`) and returns the result with the number of steps taken. Bails out with
    /// `InvalidArgument` after 2000 steps so a missing cap fails instead of hanging.
    fn count_shifts(
        z: ExactComplex,
        k: Exponent,
        cap: usize,
        stuck: bool,
    ) -> (Result<(), Error>, usize) {
        let p = 64;
        let one = ExactComplex::one(p);
        let mut z = z;
        let mut steps = 0usize;
        let r = shift_while_small(&mut z, p, k, cap, |z| {
            steps += 1;
            if steps > 2000 {
                return Err(Error::InvalidArgument);
            }
            if !stuck {
                *z = ExactComplex::add_ref(z, &one, p, RoundingMode::None);
            }
            Ok(())
        });
        (r, steps)
    }

    /// Dest p 1920 works at p_x = 1920 + 192 = 2112 bits on its first Ziv pass (`ziv_complex`:
    /// round_p(p) + 64, plus 2 words), where the Γ shift target is k = 21; dest p 1984 works at
    /// 2176, where it is k = 22. `gamma_shift_exp(p) = ceil((p + 372) / 127) + 1` steps from 21 to
    /// 22 between p = 2168 and 2169, `digamma_shift_exp(p) = ceil((p + 379) / 128) + 1` between
    /// 2181 and 2182. The shift cap 2^(k-1) + 2 follows: 2^20 + 2 at k = 21, 2^21 + 2 at k = 22.
    /// Oracle: exact integer arithmetic.
    #[test]
    fn shift_target_steps_at_the_dest_p_1920_1984_threshold() -> Result<(), crate::Error> {
        assert_eq!(gamma_shift_exp(1920 + 192), 21);
        assert_eq!(gamma_shift_exp(1984 + 192), 22);
        assert_eq!(gamma_shift_exp(2167), 21);
        assert_eq!(gamma_shift_exp(2168), 21);
        assert_eq!(gamma_shift_exp(2169), 22);
        assert_eq!(gamma_shift_exp(2170), 22);
        assert_eq!(digamma_shift_exp(2181), 21);
        assert_eq!(digamma_shift_exp(2182), 22);
        assert_eq!(shift_cap(gamma_shift_exp(2168)), (1usize << 20) + 2);
        assert_eq!(shift_cap(gamma_shift_exp(2169)), (1usize << 21) + 2);
        assert_eq!(shift_cap(digamma_shift_exp(2182)), (1usize << 21) + 2);
        Ok(())
    }

    /// A fixed 64-term tail at dest p = 2048 still has target k = 21 (~2^20 shifts).
    /// The scaled plan must take more terms and a cap of a few thousand, not a million.
    #[test]
    fn gamma_stirling_plan_keeps_high_p_shifts_small() {
        let (terms, k) = gamma_stirling_plan(2048);
        assert!(
            terms > 64 && shift_cap(k) < (1 << 16),
            "p=2048 terms={terms} k={k} cap={}",
            shift_cap(k)
        );
        let (terms, k) = gamma_stirling_plan(2048 + 192);
        assert!(
            terms > 64 && shift_cap(k) < (1 << 16),
            "p=2240 terms={terms} k={k} cap={}",
            shift_cap(k)
        );
        let (terms, k) = digamma_stirling_plan(2048);
        assert!(
            terms > 64 && shift_cap(k) < (1 << 16),
            "digamma p=2048 terms={terms} k={k} cap={}",
            shift_cap(k)
        );
    }

    /// Γ(171) = 170! has 852 significant bits and is exact at p = 1024. The old
    /// `GAMMA_FACTORIAL_MAX = 64` sent this through Stirling, which cannot certify
    /// an exact integer (same bug the real kernel fixed under ZF-37).
    #[test]
    fn test_complex_gamma_171_uses_factorial_not_stirling() -> Result<(), crate::Error> {
        let p = 1024;
        let z = ExactComplex::from_real(ExactNum::from_u32(171, p), p);
        let mut cc = Consts::new()?;
        let got = z.gamma(p, RoundingMode::ToEven, &mut cc);
        assert!(!got.is_nan(), "Γ(171) was NaN: {got:?}");
        assert!(got.im().is_zero(), "Γ(171) should be real");
        assert!(!got.re().inexact(), "Γ(171) should be an exact integer");
        let w = 1024 + 128;
        let mut acc = ExactNum::from_u8(1, w);
        for k in 2u32..171 {
            acc = acc.mul(&ExactNum::from_u32(k, w), w, RoundingMode::None);
        }
        acc.set_inexact(false);
        assert_eq!(got.re().compare(&acc), Some(0), "Γ(171) != 170!");
        Ok(())
    }

    /// Exact shift counts at the shift bound for k = 11 (target |z| >= 2^10 = 1024), from the first
    /// to the last iteration: Re z = 1023, 1024, 1025 (and with Im z = 1/2) take 1, 0, 0 steps;
    /// Re z = 2^-42 takes 1024 (the most, inside the cap 2^10 + 2 = 1026); a step that does not
    /// move z stops after exactly 1026 steps. The early refusal: from z = 1 the target is 1023
    /// steps away, so cap 1023 just reaches it and cap 1022 just misses (refused, no step).
    /// Oracle: exact integer arithmetic (|n + x + i/2| >= 1024 iff n + x >= 1024 for x + n
    /// integer, since (1023)^2 + 1/4 < 1024^2).
    #[test]
    fn shift_while_small_at_the_shift_bound() -> Result<(), crate::Error> {
        let p = 64;
        let k: Exponent = 11;
        let cap = shift_cap(k);
        assert_eq!(cap, 1026);
        let re = |n: u32| ExactNum::from_word(Word::from(n), p);
        let half = || ExactNum::from_u8(1, p).ldexp(-1, p, RoundingMode::None);
        for (x, steps) in [(1023u32, 1usize), (1024, 0), (1025, 0)] {
            assert_eq!(
                count_shifts(ExactComplex::from_real(re(x), p), k, cap, false),
                (Ok(()), steps),
                "x = {x}"
            );
            assert_eq!(
                count_shifts(ExactComplex::new(re(x), half()), k, cap, false),
                (Ok(()), steps),
                "x = {x} + i/2"
            );
        }
        let tiny =
            ExactComplex::from_real(ExactNum::from_u8(1, p).ldexp(-42, p, RoundingMode::None), p);
        assert_eq!(count_shifts(tiny, k, cap, false), (Ok(()), 1024));
        assert_eq!(
            count_shifts(ExactComplex::one(p), k, cap, true),
            (Err(Error::PrecisionRetryExhausted), cap)
        );
        assert_eq!(
            count_shifts(ExactComplex::one(p), k, 1023, false),
            (Ok(()), 1023)
        );
        assert_eq!(
            count_shifts(ExactComplex::one(p), k, 1022, false),
            (Err(Error::PrecisionRetryExhausted), 0)
        );
        Ok(())
    }

    /// `shift_while_small` steps only while |z| < 2^(k-1): from a tiny Re z the real z -> z+1 step
    /// takes exactly 2^(k-1) steps (inside the cap), a large z takes none, a step that does not move
    /// z stops after exactly `cap` steps with `PrecisionRetryExhausted` (this is what broken
    /// arithmetic, e.g. every result 0, used to turn into an endless loop), and a step error comes back
    /// as is.
    #[test]
    fn shift_while_small_is_bounded() -> Result<(), crate::Error> {
        let p = 64;
        let k: Exponent = 6;

        let tiny =
            ExactComplex::from_real(ExactNum::from_u8(3, p).ldexp(-42, p, RoundingMode::None), p);
        assert_eq!(count_shifts(tiny, k, shift_cap(k), false), (Ok(()), 32));

        let big = ExactComplex::from_real(ExactNum::from_u8(32, p), p);
        assert_eq!(count_shifts(big, k, shift_cap(k), false), (Ok(()), 0));

        let stuck = ExactComplex::one(p);
        assert_eq!(
            count_shifts(stuck, k, shift_cap(k), true),
            (Err(Error::PrecisionRetryExhausted), shift_cap(k))
        );

        let mut failing = ExactComplex::one(p);
        let r = shift_while_small(&mut failing, p, k, shift_cap(k), |_| {
            Err(Error::MemoryAllocation)
        });
        assert_eq!(r, Err(Error::MemoryAllocation));
        Ok(())
    }

    /// When 2^(k-1) - |z| (the fewest steps that can reach the target) is above `cap`,
    /// `shift_while_small` fails at once with no steps; at exactly `cap` it runs. From z = 1 with
    /// k = 6 the target 32 is 31 steps away.
    #[test]
    fn shift_while_small_fails_early_past_the_cap() -> Result<(), crate::Error> {
        let p = 64;
        let one = || ExactComplex::one(p);
        assert_eq!(count_shifts(one(), 6, 31, false), (Ok(()), 31));
        assert_eq!(
            count_shifts(one(), 6, 30, false),
            (Err(Error::PrecisionRetryExhausted), 0)
        );
        assert_eq!(
            count_shifts(one(), 6, 31, true),
            (Err(Error::PrecisionRetryExhausted), 31)
        );
        // The early check uses |z|, not Re z: from z = 8i the target 16 is at least 8 steps away
        // (|n + 8i| >= 16 first at n = 14), so cap 14 runs all 14 and cap 7 fails at once.
        let tilted = || ExactComplex::new(ExactNum::new(p), ExactNum::from_u8(8, p));
        assert_eq!(count_shifts(tilted(), 5, 14, false), (Ok(()), 14));
        assert_eq!(
            count_shifts(tilted(), 5, 7, false),
            (Err(Error::PrecisionRetryExhausted), 0)
        );
        // A target 2^(k-1) too large to represent fails without stepping.
        assert_eq!(
            count_shifts(one(), Exponent::MAX, usize::MAX, false),
            (Err(Error::PrecisionRetryExhausted), 0)
        );
        Ok(())
    }

    /// A tiny `Re z > 0` takes the most integer shifts before Stirling (about 2^(k-1) for the
    /// shift target k), and `Re z = 0` goes straight into the shift loop: both must still finish
    /// with the right value. Checked against the real kernel on the axis (x·Γ(x), ln Γ(x)),
    /// with x·Γ(x) = Γ(1+x), with z·(ψ(z+1) − ψ(z)) = 1, and on the imaginary axis.
    #[test]
    fn gamma_family_tiny_argument_takes_the_most_shifts() -> Result<(), crate::Error> {
        let rm = RoundingMode::ToEven;
        // Each Ziv pass below takes up to 2^(k-1)+2 shifts, and the passes only stop once
        // `ExactNum::try_set_precision` certifies. Check that it certifies an exact, plainly
        // roundable value (1 + 2^-127 to 64 bits is 1) before the loops, so a broken
        // certification fails here at once instead of running the loops to their budget.
        let mut probe = ExactNum::from_words(&[1, 0x8000_0000_0000_0000], crate::Sign::Pos, 1);
        assert!(
            probe.try_set_precision(64, rm, 128),
            "try_set_precision does not certify 1 + 2^-127 at 64 bits: {probe:?}"
        );
        assert_eq!(
            probe.compare(&ExactNum::from_u8(1, 64)),
            Some(0),
            "{probe:?}"
        );
        let mut cc = Consts::new()?;
        for p in [64usize, 256, 1024] {
            let one = ExactNum::from_u8(1, p);
            let mk_x = || ExactNum::from_u8(3, p).ldexp(-42, p, RoundingMode::None);
            let x = mk_x();
            let z = ExactComplex::from_real(mk_x(), p);

            let gz = z.gamma(p, rm, &mut cc);
            let xg_c = ExactNum::mul_ref(gz.re(), &x, p, RoundingMode::None);
            let xg_r = x.gamma(p, rm, &mut cc).mul(&x, p, RoundingMode::None);
            let g1px = one.add(&x, p, RoundingMode::None).gamma(p, rm, &mut cc);
            assert!(
                near(&xg_c, &xg_r, p),
                "p={p}: x*gamma(x) {xg_c:?} vs real {xg_r:?}"
            );
            assert!(
                near(&xg_c, &g1px, p),
                "p={p}: x*gamma(x) {xg_c:?} vs gamma(1+x) {g1px:?}"
            );
            assert!(tiny(gz.im(), p), "p={p}: im gamma(x) {:?}", gz.im());

            let lgz = z.ln_gamma(p, rm, &mut cc);
            let lgx = x.ln_gamma(p, rm, &mut cc);
            assert!(
                near(lgz.re(), &lgx, p),
                "p={p}: ln_gamma {:?} vs real {lgx:?}",
                lgz.re()
            );
            assert!(tiny(lgz.im(), p), "p={p}: im ln_gamma(x) {:?}", lgz.im());

            let zp1 = ExactComplex::one(p).add(&z, p, rm);
            let dpsi = zp1
                .digamma(p, rm, &mut cc)
                .sub(&z.digamma(p, rm, &mut cc), p, rm);
            let unit = dpsi.mul(&z, p, RoundingMode::None);
            assert!(
                cnear(&unit, &ExactComplex::one(p), p),
                "p={p}: z*dpsi {unit:?}"
            );

            let w = ExactComplex::new(ExactNum::new(p), mk_x());
            let wg = w.gamma(p, rm, &mut cc).mul(&w, p, RoundingMode::None);
            let g1pw = ExactComplex::one(p).add(&w, p, rm).gamma(p, rm, &mut cc);
            assert!(
                cnear(&wg, &g1pw, p),
                "p={p}: w*gamma(w) {wg:?} vs gamma(1+w) {g1pw:?}"
            );
        }
        Ok(())
    }

    /// `erf` / `erfc` of a real `z` (an `mpfr_compare_complex_real_axis` case): the imaginary part is exactly
    /// zero, so it must come out zero instead of exhausting the retry budget on rounding noise.
    #[test]
    #[cfg(target_pointer_width = "64")]
    fn erf_real_axis_im_is_zero() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let p = 128;
        let x = ExactNum::from_words(
            &[5566959752443415075, 17370124544999774818],
            crate::Sign::Neg,
            -1,
        );
        let z = ExactComplex::from_real(x, p);
        for rm in [
            RoundingMode::ToEven,
            RoundingMode::ToOdd,
            RoundingMode::Up,
            RoundingMode::Down,
            RoundingMode::ToZero,
            RoundingMode::FromZero,
            RoundingMode::None,
        ] {
            let e = z.erf(p, rm, &mut cc);
            assert!(e.im().is_zero() && !e.re().is_nan(), "erf {rm:?}: {e:?}");
            let c = z.erfc(p, rm, &mut cc);
            assert!(c.im().is_zero() && !c.re().is_nan(), "erfc {rm:?}: {c:?}");
        }
        Ok(())
    }

    /// ψ(1) = −γ must use the integer closed form. The Stirling walk from 1
    /// to 2^{e_min−1} is ~2^{40} steps at p=65536 and does not finish.
    #[test]
    fn digamma_one_is_neg_euler_gamma_at_high_p() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let p = 65536;
        let z = ExactComplex::from_real(ExactNum::from_u8(1, p), p);
        let t0 = std::time::Instant::now();
        let psi = z.digamma_at(p, &mut cc);
        assert!(
            t0.elapsed().as_secs() < 5,
            "digamma(1) at p=65536 took {:?}",
            t0.elapsed()
        );
        assert!(!psi.is_nan(), "digamma(1) NaN: {:?}", psi);
        assert!(psi.im().is_zero());
        let g = cc.euler_gamma(p, RoundingMode::None);
        assert_eq!(psi.re().compare(&g.neg()), Some(0));
        Ok(())
    }

    /// ψ(1/2) = −γ − 2 ln 2 (DLMF 5.4.13 / 5.4.14). At p=65536 the 1024-term
    /// Stirling plan still wants |z|≥2^{39}, and the 2^{20} walk cap refused
    /// at once (PrecisionRetryExhausted). That cap does not come from the
    /// remainder bound (P9). This value now uses the half-integer closed
    /// form (same as the real kernel); term and shift counts on the Stirling
    /// path follow the bound, with no walk cap, so the p=65536 Kelvin path
    /// stays short.
    #[test]
    fn digamma_half_at_high_p_from_error_bound() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let p = 65536;
        let half = ExactNum::from_u8(1, p).ldexp(-1, p, RoundingMode::None);
        let z = ExactComplex::from_real(half, p);
        let t0 = std::time::Instant::now();
        let psi = z.digamma(p, RoundingMode::ToEven, &mut cc);
        assert!(
            t0.elapsed().as_secs() < 60,
            "digamma(1/2) at p=65536 took {:?}",
            t0.elapsed()
        );
        assert!(!psi.is_nan(), "digamma(1/2) NaN: {:?}", psi);
        assert!(psi.im().is_zero(), "digamma(1/2) im {}", psi.im());
        let wrk = p + 2 * crate::WORD_BIT_SIZE;
        let g = cc.euler_gamma(wrk, RoundingMode::None);
        let ln2 = cc.ln_2(wrk, RoundingMode::None);
        let two_ln2 = ln2.ldexp(1, wrk, RoundingMode::None);
        let mut want = g.neg().sub(&two_ln2, wrk, RoundingMode::None);
        assert!(want.try_set_precision(p, RoundingMode::ToEven, wrk));
        assert_eq!(psi.re().compare(&want), Some(0));
        Ok(())
    }
}
