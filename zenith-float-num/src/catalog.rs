//! Additive catalog specials on [`ExactNum`]: Scorer, Kelvin, Struve, Anger–Weber,
//! Clausen, Barnes \(G\), polygamma / Hurwitz, inverse Jacobi, \({}_pF_q\), Lambert \(W\),
//! and the polylogarithm. SoftFloat only (integer limbs; no IEEE arithmetic).
//!
//! Meijer \(G\), Fox \(H\), and Heun live in `heun_fox`.

use crate::defs::WORD_BIT_SIZE;
use crate::quadrature::tanh_sinh;
use crate::Consts;
use crate::Error;
use crate::ExactComplex;
use crate::ExactNum;
use crate::RoundingMode;
use alloc::vec::Vec;

pub(crate) fn work_p(p: usize) -> usize {
    p.saturating_add(WORD_BIT_SIZE.saturating_mul(2))
}

pub(crate) fn cat_nan() -> ExactNum {
    ExactNum::nan(Some(Error::InvalidArgument))
}

pub(crate) fn finite(x: &ExactNum) -> bool {
    !x.is_nan() && !x.is_inf()
}

pub(crate) fn finish(mut y: ExactNum, p: usize, rm: RoundingMode) -> ExactNum {
    let _ = y.set_precision(p, rm);
    y
}

pub(crate) fn series_cap(p: usize) -> usize {
    p.saturating_add(32).max(64)
}

fn even_bernoulli_ext(n: usize, p: usize) -> Vec<ExactNum> {
    match crate::ops::special::even_bernoulli(n, p) {
        Ok(v) => v.into_iter().map(ExactNum::from).collect(),
        Err(e) => (0..n).map(|_| ExactNum::nan(Some(e))).collect(),
    }
}

fn factorial_n(n: usize, p: usize) -> ExactNum {
    let mut acc = ExactNum::from_u8(1, p);
    for k in 2..=n {
        acc = acc.mul(&ExactNum::from_u32(k as u32, p), p, RoundingMode::None);
    }
    acc
}

fn abs_lt(x: &ExactNum, bound: u32, p: usize) -> bool {
    let a = x.abs();
    matches!(a.cmp(&ExactNum::from_u32(bound, p)), Some(c) if c < 0)
}

pub(crate) fn small_int(x: &ExactNum, p: usize, lo: i32, hi: i32) -> Option<i32> {
    if !x.is_int() {
        return None;
    }
    (lo..=hi).find(|&n| x.cmp(&ExactNum::from_i32(n, p)) == Some(0))
}

/// \(\zeta(2k)=(-1)^{k-1}B_{2k}(2\pi)^{2k}/(2\,(2k)!)\).
fn zeta_even(k: usize, p: usize, cc: &mut Consts) -> ExactNum {
    if k == 0 {
        return cat_nan();
    }
    let bs = even_bernoulli_ext(k, p);
    let b = &bs[k - 1];
    let pi = cc.pi(p, RoundingMode::None);
    let two = ExactNum::from_u8(2, p);
    let two_pi = two.mul(&pi, p, RoundingMode::None);
    let mut pow = ExactNum::from_u8(1, p);
    for _ in 0..(2 * k) {
        pow = pow.mul(&two_pi, p, RoundingMode::None);
    }
    let den = two.mul(&factorial_n(2 * k, p), p, RoundingMode::None);
    let mut z = b
        .mul(&pow, p, RoundingMode::None)
        .div(&den, p, RoundingMode::None);
    if k % 2 == 0 {
        z = z.neg();
    }
    z
}

/// Euler–Maclaurin for \(\zeta(s,a)=\sum_{k=0}^\infty(k+a)^{-s}\), integer \(s\ge 2\), \(a>0\).
fn hurwitz_em(s: u32, a: &ExactNum, p: usize, cc: &mut Consts) -> ExactNum {
    if s < 2 || !finite(a) || a.is_negative() || a.is_zero() {
        return cat_nan();
    }
    let n_front = (16 + p / 4).max(24);
    let m_bern = (8 + p / 16).clamp(4, 48);
    let mut sum = ExactNum::new(p);
    for k in 0..n_front {
        let t = a.add(&ExactNum::from_u32(k as u32, p), p, RoundingMode::None);
        let term = ExactNum::from_u8(1, p).div(
            &t.powsi(s as isize, p, RoundingMode::None),
            p,
            RoundingMode::None,
        );
        sum = sum.add(&term, p, RoundingMode::None);
    }
    let n_a = a.add(
        &ExactNum::from_u32(n_front as u32, p),
        p,
        RoundingMode::None,
    );
    let sm1 = ExactNum::from_u32(s - 1, p);
    let tail_int = ExactNum::from_u8(1, p).div(
        &n_a.powsi((s - 1) as isize, p, RoundingMode::None)
            .mul(&sm1, p, RoundingMode::None),
        p,
        RoundingMode::None,
    );
    let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, RoundingMode::None);
    let tail_half = half.div(
        &n_a.powsi(s as isize, p, RoundingMode::None),
        p,
        RoundingMode::None,
    );
    sum = sum
        .add(&tail_int, p, RoundingMode::None)
        .add(&tail_half, p, RoundingMode::None);
    let bs = even_bernoulli_ext(m_bern, p);
    for j in 1..=m_bern {
        let b2j = &bs[j - 1];
        let mut rising = ExactNum::from_u8(1, p);
        for i in 0..(2 * j - 1) {
            rising = rising.mul(&ExactNum::from_u32(s + i as u32, p), p, RoundingMode::None);
        }
        let fact = factorial_n(2 * j, p);
        let pow = n_a.powsi((s + 2 * j as u32 - 1) as isize, p, RoundingMode::None);
        let term = b2j.mul(&rising, p, RoundingMode::None).div(
            &fact.mul(&pow, p, RoundingMode::None),
            p,
            RoundingMode::None,
        );
        sum = sum.add(&term, p, RoundingMode::None);
        let _ = cc;
    }
    sum
}

/// \(m\)-th derivative of \(\ln x/x^2\) at `x` (Leibniz).
fn deriv_ln_over_x2(x: &ExactNum, m: usize, p: usize, cc: &mut Consts) -> ExactNum {
    let mut acc = ExactNum::new(p);
    let lnx = x.ln(p, RoundingMode::None, cc);
    for k in 0..=m {
        let u = if k == 0 {
            lnx.clone()
        } else {
            let mut t = factorial_n(k - 1, p).div(
                &x.powsi(k as isize, p, RoundingMode::None),
                p,
                RoundingMode::None,
            );
            if k % 2 == 0 {
                t = t.neg();
            }
            t
        };
        let j = m - k;
        let mut v = factorial_n(j + 1, p).div(
            &x.powsi((2 + j) as isize, p, RoundingMode::None),
            p,
            RoundingMode::None,
        );
        if j % 2 == 1 {
            v = v.neg();
        }
        let bin = factorial_n(m, p).div(
            &factorial_n(k, p).mul(&factorial_n(m - k, p), p, RoundingMode::None),
            p,
            RoundingMode::None,
        );
        acc = acc.add(
            &bin.mul(&u, p, RoundingMode::None)
                .mul(&v, p, RoundingMode::None),
            p,
            RoundingMode::None,
        );
    }
    acc
}

/// \(\zeta'(2)=-\sum_{n=1}^\infty \ln n/n^2\).
fn zeta_prime_2(p: usize, cc: &mut Consts) -> ExactNum {
    let n_front = (48 + p / 2).clamp(64, 2048);
    let mut sum = ExactNum::new(p);
    for n in 2..=n_front {
        let nf = ExactNum::from_u32(n as u32, p);
        let ln = nf.ln(p, RoundingMode::None, cc);
        let den = nf.mul(&nf, p, RoundingMode::None);
        sum = sum.add(&ln.div(&den, p, RoundingMode::None), p, RoundingMode::None);
    }
    let n0 = ExactNum::from_u32((n_front + 1) as u32, p);
    let ln = n0.ln(p, RoundingMode::None, cc);
    let one = ExactNum::from_u8(1, p);
    let two = ExactNum::from_u8(2, p);
    let half = one.div(&two, p, RoundingMode::None);
    let n2 = n0.mul(&n0, p, RoundingMode::None);
    let integral = ln
        .add(&one, p, RoundingMode::None)
        .div(&n0, p, RoundingMode::None);
    let mid = half
        .mul(&ln, p, RoundingMode::None)
        .div(&n2, p, RoundingMode::None);
    sum = sum
        .add(&integral, p, RoundingMode::None)
        .add(&mid, p, RoundingMode::None);
    let m_bern = (6 + p / 32).clamp(4, 12);
    let bs = even_bernoulli_ext(m_bern, p);
    for j in 1..=m_bern {
        let b = &bs[j - 1];
        let fact = factorial_n(2 * j, p);
        let der = deriv_ln_over_x2(&n0, 2 * j - 1, p, cc);
        let term = b
            .mul(&der, p, RoundingMode::None)
            .div(&fact, p, RoundingMode::None);
        sum = sum.sub(&term, p, RoundingMode::None);
    }
    sum.neg()
}

/// \(\ln A = (\gamma+\ln(2\pi))/12 - \zeta'(2)/(2\pi^2)\).
fn ln_glaisher(p: usize, cc: &mut Consts) -> ExactNum {
    let g = cc.euler_gamma(p, RoundingMode::None);
    let two = ExactNum::from_u8(2, p);
    let pi = cc.pi(p, RoundingMode::None);
    let twopi = two.mul(&pi, p, RoundingMode::None);
    let twelve = ExactNum::from_u32(12, p);
    let a = g
        .add(&twopi.ln(p, RoundingMode::None, cc), p, RoundingMode::None)
        .div(&twelve, p, RoundingMode::None);
    let zp = zeta_prime_2(p, cc);
    let pi2 = pi.mul(&pi, p, RoundingMode::None);
    let b = zp.div(&two.mul(&pi2, p, RoundingMode::None), p, RoundingMode::None);
    a.sub(&b, p, RoundingMode::None)
}

fn reduce_clausen_angle(theta: &ExactNum, p: usize, cc: &mut Consts) -> ExactNum {
    let pi = cc.pi(p, RoundingMode::None);
    let two_pi = ExactNum::from_u8(2, p).mul(&pi, p, RoundingMode::None);
    let mut t = theta.rem(&two_pi);
    if matches!(t.cmp(&pi), Some(c) if c > 0) {
        t = t.sub(&two_pi, p, RoundingMode::None);
    } else if matches!(t.cmp(&pi.neg()), Some(c) if c < 0) {
        t = t.add(&two_pi, p, RoundingMode::None);
    }
    t
}

/// Inhomogeneous Airy series for Scorer \(Gi\) (shared \(a_2=1/(2\pi)\)).
fn scorer_gi_series(x: &ExactNum, p: usize, cc: &mut Consts) -> ExactNum {
    let bi0 = ExactNum::new(p).bi(p, RoundingMode::None, cc);
    let bip0 = ExactNum::new(p).bi_prime(p, RoundingMode::None, cc);
    let three = ExactNum::from_u8(3, p);
    let a0 = bi0.div(&three, p, RoundingMode::None);
    let a1 = bip0.div(&three, p, RoundingMode::None);
    let pi = cc.pi(p, RoundingMode::None);
    let two = ExactNum::from_u8(2, p);
    // Gi'' − z Gi = −1/π ⇒ 2 a₂ = −1/π.
    let a2 = ExactNum::from_u8(1, p)
        .div(&two.mul(&pi, p, RoundingMode::None), p, RoundingMode::None)
        .neg();
    let mut am3 = a0.clone();
    let mut am2 = a1.clone();
    let mut am1 = a2.clone();
    let mut xn = ExactNum::from_u8(1, p);
    let mut sum = a0.clone();
    xn = xn.mul(x, p, RoundingMode::None);
    sum = sum.add(&a1.mul(&xn, p, RoundingMode::None), p, RoundingMode::None);
    xn = xn.mul(x, p, RoundingMode::None);
    sum = sum.add(&a2.mul(&xn, p, RoundingMode::None), p, RoundingMode::None);
    for n in 3..=series_cap(p) {
        let nf = ExactNum::from_u32(n as u32, p);
        let nm1 = ExactNum::from_u32((n - 1) as u32, p);
        let an = am3.div(&nf.mul(&nm1, p, RoundingMode::None), p, RoundingMode::None);
        xn = xn.mul(x, p, RoundingMode::None);
        let term = an.mul(&xn, p, RoundingMode::None);
        sum = sum.add(&term, p, RoundingMode::None);
        if term.is_zero()
            || term
                .abs()
                .exponent()
                .is_some_and(|e| (e as i64) + (p as i64) < sum.abs().exponent().unwrap_or(0) as i64)
        {
            break;
        }
        am3 = am2;
        am2 = am1;
        am1 = an;
    }
    sum
}

/// \(Gi(x)\sim(\pi x)^{-1}\sum (3k+1)^{\overline{\times}}/x^{3k}\) for \(x\to+\infty\).
fn scorer_gi_asymp_plus(x: &ExactNum, p: usize, cc: &mut Consts) -> ExactNum {
    let pi = cc.pi(p, RoundingMode::None);
    let x3 = x
        .mul(x, p, RoundingMode::None)
        .mul(x, p, RoundingMode::None);
    let mut term = ExactNum::from_u8(1, p);
    let mut sum = term.clone();
    for k in 0..=series_cap(p) {
        let fac = ExactNum::from_u32((3 * k + 1) as u32, p);
        term = term
            .mul(&fac, p, RoundingMode::None)
            .div(&x3, p, RoundingMode::None);
        sum = sum.add(&term, p, RoundingMode::None);
        if term.is_zero()
            || term
                .abs()
                .exponent()
                .is_some_and(|e| (e as i64) + (p as i64) < 0)
        {
            break;
        }
    }
    sum.div(&pi.mul(x, p, RoundingMode::None), p, RoundingMode::None)
}

/// \(Hi(x)\sim-(\pi x)^{-1}\sum (-1)^k(3k+1)^{\overline{\times}}/x^{3k}\) as \(x\to-\infty\).
fn scorer_hi_asymp_minus(x: &ExactNum, p: usize, cc: &mut Consts) -> ExactNum {
    let pi = cc.pi(p, RoundingMode::None);
    let x3 = x
        .mul(x, p, RoundingMode::None)
        .mul(x, p, RoundingMode::None);
    let mut term = ExactNum::from_u8(1, p);
    let mut sum = term.clone();
    let mut sign = -1i32;
    for k in 0..=series_cap(p) {
        let fac = ExactNum::from_u32((3 * k + 1) as u32, p);
        term = term
            .mul(&fac, p, RoundingMode::None)
            .div(&x3, p, RoundingMode::None);
        let t = if sign < 0 { term.clone().neg() } else { term.clone() };
        sum = sum.add(&t, p, RoundingMode::None);
        sign = -sign;
        if term.is_zero()
            || term
                .abs()
                .exponent()
                .is_some_and(|e| (e as i64) + (p as i64) < 0)
        {
            break;
        }
    }
    sum.div(&pi.mul(x, p, RoundingMode::None), p, RoundingMode::None)
        .neg()
}

fn kelvin_j_arg(x: &ExactNum, p: usize, _cc: &mut Consts) -> ExactComplex {
    // x exp(3πi/4) = x (−√2/2 + i √2/2)
    let two = ExactNum::from_u8(2, p);
    let s = two.sqrt(p, RoundingMode::None);
    let h = x.div(&s, p, RoundingMode::None);
    ExactComplex::new(h.neg(), h)
}

fn kelvin_k_arg(x: &ExactNum, p: usize, _cc: &mut Consts) -> ExactComplex {
    // x exp(πi/4) = x (√2/2 + i √2/2)
    let two = ExactNum::from_u8(2, p);
    let s = two.sqrt(p, RoundingMode::None);
    let h = x.div(&s, p, RoundingMode::None);
    ExactComplex::new(h.clone(), h)
}

impl ExactNum {
    /// Scorer \(\mathrm{Gi}(\mathrm{self})\). Entire; NaN in → NaN out.
    ///
    /// # Precision
    ///
    /// - Algorithm: Taylor (\(w''-zw=1/\pi\)) for `|x| < 8`; \(1/(\pi x)\) asymptotic for
    ///   large `+x`; \(\mathrm{Hi}_\mathrm{asymp}-\mathrm{Bi}/3\) for large `-x`.
    /// - Bound: extra word of working precision.
    pub fn scorer_gi(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) {
            return if self.is_nan() { self.clone() } else { cat_nan() };
        }
        let wrk = work_p(p);
        let y = if abs_lt(self, 8, wrk) {
            scorer_gi_series(self, wrk, cc)
        } else if self.is_positive() {
            scorer_gi_asymp_plus(self, wrk, cc)
        } else {
            let hi = scorer_hi_asymp_minus(self, wrk, cc);
            let bi = self.bi(wrk, RoundingMode::None, cc);
            bi.sub(&hi, wrk, RoundingMode::None)
        };
        finish(y, p, rm)
    }

    /// Scorer \(\mathrm{Hi}(\mathrm{self})=\mathrm{Bi}(\mathrm{self})-\mathrm{Gi}(\mathrm{self})\).
    pub fn scorer_hi(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) {
            return if self.is_nan() { self.clone() } else { cat_nan() };
        }
        let wrk = work_p(p);
        let y = if self.is_negative() && !abs_lt(self, 8, wrk) {
            scorer_hi_asymp_minus(self, wrk, cc)
        } else {
            let gi = self.scorer_gi(wrk, RoundingMode::None, cc);
            let bi = self.bi(wrk, RoundingMode::None, cc);
            bi.sub(&gi, wrk, RoundingMode::None)
        };
        finish(y, p, rm)
    }

    /// Kelvin \(\mathrm{ber}_\nu(\mathrm{self})=\mathrm{Re}\,J_\nu(x e^{3\pi i/4})\).
    pub fn kelvin_ber_nu(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) || !finite(nu) {
            return cat_nan();
        }
        let wrk = work_p(p);
        let z = kelvin_j_arg(self, wrk, cc);
        let n = ExactComplex::from_real(nu.clone(), wrk);
        finish(
            z.bessel_j_nu(&n, wrk, RoundingMode::None, cc).re().clone(),
            p,
            rm,
        )
    }

    /// Kelvin \(\mathrm{bei}_\nu(\mathrm{self})=\mathrm{Im}\,J_\nu(x e^{3\pi i/4})\).
    pub fn kelvin_bei_nu(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) || !finite(nu) {
            return cat_nan();
        }
        let wrk = work_p(p);
        let z = kelvin_j_arg(self, wrk, cc);
        let n = ExactComplex::from_real(nu.clone(), wrk);
        finish(
            z.bessel_j_nu(&n, wrk, RoundingMode::None, cc).im().clone(),
            p,
            rm,
        )
    }

    /// Kelvin \(\mathrm{ker}_\nu(\mathrm{self})=\mathrm{Re}\,e^{-\nu\pi i/2}K_\nu(x e^{\pi i/4})\).
    pub fn kelvin_ker_nu(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) || !finite(nu) || self.is_zero() {
            return if self.is_zero() { crate::ext::INF_POS } else { cat_nan() };
        }
        let wrk = work_p(p);
        let z = kelvin_k_arg(self, wrk, cc);
        let n = ExactComplex::from_real(nu.clone(), wrk);
        let k = z.bessel_k(&n, wrk, RoundingMode::None, cc);
        let half_pi = ExactComplex::new(
            ExactNum::new(wrk),
            cc.pi(wrk, RoundingMode::None)
                .ldexp(-1, wrk, RoundingMode::None)
                .neg(),
        );
        let phase = n
            .mul(&half_pi, wrk, RoundingMode::None)
            .exp(wrk, RoundingMode::None, cc);
        finish(phase.mul(&k, wrk, RoundingMode::None).re().clone(), p, rm)
    }

    /// Kelvin \(\mathrm{kei}_\nu(\mathrm{self})=\mathrm{Im}\,e^{-\nu\pi i/2}K_\nu(x e^{\pi i/4})\).
    pub fn kelvin_kei_nu(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) || !finite(nu) || self.is_zero() {
            return cat_nan();
        }
        let wrk = work_p(p);
        let z = kelvin_k_arg(self, wrk, cc);
        let n = ExactComplex::from_real(nu.clone(), wrk);
        let k = z.bessel_k(&n, wrk, RoundingMode::None, cc);
        let half_pi = ExactComplex::new(
            ExactNum::new(wrk),
            cc.pi(wrk, RoundingMode::None)
                .ldexp(-1, wrk, RoundingMode::None)
                .neg(),
        );
        let phase = n
            .mul(&half_pi, wrk, RoundingMode::None)
            .exp(wrk, RoundingMode::None, cc);
        finish(phase.mul(&k, wrk, RoundingMode::None).im().clone(), p, rm)
    }

    /// Order-0 Kelvin \(\mathrm{ber}(\mathrm{self})\).
    pub fn kelvin_ber(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.kelvin_ber_nu(&ExactNum::new(p), p, rm, cc)
    }

    /// Order-0 Kelvin \(\mathrm{bei}(\mathrm{self})\).
    pub fn kelvin_bei(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.kelvin_bei_nu(&ExactNum::new(p), p, rm, cc)
    }

    /// Order-0 Kelvin \(\mathrm{ker}(\mathrm{self})\).
    pub fn kelvin_ker(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.kelvin_ker_nu(&ExactNum::new(p), p, rm, cc)
    }

    /// Order-0 Kelvin \(\mathrm{kei}(\mathrm{self})\).
    pub fn kelvin_kei(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.kelvin_kei_nu(&ExactNum::new(p), p, rm, cc)
    }

    /// Struve \(\mathbf{H}_\nu(\mathrm{self})\).
    ///
    /// # Precision
    ///
    /// - Algorithm: power series (DLMF 11.2.1). For large `|x|` uses \(Y_\nu+(2/\pi)\) particular.
    pub fn struve_h(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) || !finite(nu) {
            return cat_nan();
        }
        let wrk = work_p(p);
        if !abs_lt(self, 16, wrk) && self.is_positive() {
            let y = self.bessel_y(nu, wrk, RoundingMode::None, cc);
            let two = ExactNum::from_u8(2, wrk);
            let pi = cc.pi(wrk, RoundingMode::None);
            let half = ExactNum::from_u8(1, wrk).div(&two, wrk, RoundingMode::None);
            let pref = self
                .div(&two, wrk, RoundingMode::None)
                .pow(
                    &nu.sub(&half, wrk, RoundingMode::None),
                    wrk,
                    RoundingMode::None,
                    cc,
                )
                .mul(&two, wrk, RoundingMode::None)
                .div(
                    &pi.mul(
                        &nu.add(&half, wrk, RoundingMode::None)
                            .gamma(wrk, RoundingMode::None, cc),
                        wrk,
                        RoundingMode::None,
                    )
                    .mul(
                        &half.gamma(wrk, RoundingMode::None, cc),
                        wrk,
                        RoundingMode::None,
                    ),
                    wrk,
                    RoundingMode::None,
                );
            return finish(y.add(&pref, wrk, RoundingMode::None), p, rm);
        }
        let two = ExactNum::from_u8(2, wrk);
        let half = self.div(&two, wrk, RoundingMode::None);
        let z2 = half.mul(&half, wrk, RoundingMode::None);
        let three_half = ExactNum::from_u8(3, wrk).div(&two, wrk, RoundingMode::None);
        let mut sum = ExactNum::new(wrk);
        let mut zk = ExactNum::from_u8(1, wrk);
        let mut sign = 1i32;
        for k in 0..=series_cap(wrk) {
            let gk1 = ExactNum::from_u32(k as u32, wrk)
                .add(&three_half, wrk, RoundingMode::None)
                .gamma(wrk, RoundingMode::None, cc);
            let gk2 = ExactNum::from_u32(k as u32, wrk)
                .add(nu, wrk, RoundingMode::None)
                .add(&three_half, wrk, RoundingMode::None)
                .gamma(wrk, RoundingMode::None, cc);
            let mut term = zk.div(
                &gk1.mul(&gk2, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            );
            if sign < 0 {
                term = term.neg();
            }
            sum = sum.add(&term, wrk, RoundingMode::None);
            if term.is_zero()
                || term.abs().exponent().is_some_and(|e| {
                    (e as i64) + (p as i64) < sum.abs().exponent().unwrap_or(0) as i64
                })
            {
                break;
            }
            zk = zk.mul(&z2, wrk, RoundingMode::None);
            sign = -sign;
        }
        let pref = half.pow(
            &nu.add(&ExactNum::from_u8(1, wrk), wrk, RoundingMode::None),
            wrk,
            RoundingMode::None,
            cc,
        );
        finish(pref.mul(&sum, wrk, RoundingMode::None), p, rm)
    }

    /// Anger \(\mathbf{J}_\nu(\mathrm{self})\). Integer \(\nu\) reduces to Bessel \(J\).
    pub fn anger_j(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) || !finite(nu) {
            return cat_nan();
        }
        let wrk = work_p(p);
        if let Some(n) = small_int(nu, wrk, -64, 64) {
            let j = self.bessel_j(n.unsigned_abs() as usize, wrk, RoundingMode::None, cc);
            return finish(if n < 0 && n % 2 != 0 { j.neg() } else { j }, p, rm);
        }
        let z = self.clone();
        let n = nu.clone();
        let pi = cc.pi(wrk, RoundingMode::None);
        let zero = ExactNum::new(wrk);
        let q = tanh_sinh(
            |th, pw, r, c| {
                let inner = n
                    .mul(th, pw, r)
                    .sub(&z.mul(&th.sin(pw, r, c), pw, r), pw, r);
                inner.cos(pw, r, c)
            },
            &zero,
            &pi,
            wrk,
            RoundingMode::None,
            cc,
        );
        match q {
            Some(v) => finish(v.div(&pi, wrk, RoundingMode::None), p, rm),
            None => cat_nan(),
        }
    }

    /// Weber \(\mathbf{E}_\nu(\mathrm{self})\).
    pub fn weber_e(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) || !finite(nu) {
            return cat_nan();
        }
        let wrk = work_p(p);
        let z = self.clone();
        let n = nu.clone();
        let pi = cc.pi(wrk, RoundingMode::None);
        let zero = ExactNum::new(wrk);
        let q = tanh_sinh(
            |th, pw, r, c| {
                let inner = n
                    .mul(th, pw, r)
                    .sub(&z.mul(&th.sin(pw, r, c), pw, r), pw, r);
                inner.sin(pw, r, c)
            },
            &zero,
            &pi,
            wrk,
            RoundingMode::None,
            cc,
        );
        match q {
            Some(v) => finish(v.div(&pi, wrk, RoundingMode::None), p, rm),
            None => cat_nan(),
        }
    }

    /// Clausen \(\mathrm{Cl}_2(\mathrm{self})=\sum_{k\ge1}\sin(k\theta)/k^2=-\int_0^\theta\ln\lvert 2\sin(t/2)\rvert\,dt\).
    pub fn clausen_cl2(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) {
            return cat_nan();
        }
        let wrk = work_p(p);
        let th = reduce_clausen_angle(self, wrk, cc);
        if th.is_zero() {
            return finish(ExactNum::new(wrk), p, rm);
        }
        let pi = cc.pi(wrk, RoundingMode::None);
        if th.cmp(&pi) == Some(0) || th.cmp(&pi.neg()) == Some(0) {
            return finish(ExactNum::new(wrk), p, rm);
        }
        let sign = if th.is_negative() { -1i32 } else { 1 };
        let th_abs = th.abs();
        let two = ExactNum::from_u8(2, wrk);
        let zero = ExactNum::new(wrk);
        let q = tanh_sinh(
            |t, pw, r, c| {
                let s = t.div(&two, pw, r).sin(pw, r, c).abs().mul(&two, pw, r);
                if s.is_zero() {
                    return ExactNum::from_i32(-64, pw);
                }
                s.ln(pw, r, c)
            },
            &zero,
            &th_abs,
            wrk,
            RoundingMode::None,
            cc,
        );
        match q {
            Some(v) => {
                let mut y = v.neg();
                if sign < 0 {
                    y = y.neg();
                }
                finish(y, p, rm)
            }
            None => cat_nan(),
        }
    }

    /// Clausen \(\mathrm{Cl}_3(\mathrm{self})=\sum_{k\ge1}\cos(k\theta)/k^3=\zeta(3)-\int_0^\theta\mathrm{Cl}_2(t)\,dt\).
    pub fn clausen_cl3(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) {
            return cat_nan();
        }
        let wrk = work_p(p);
        let th = reduce_clausen_angle(self, wrk, cc);
        let z3 = ExactNum::from_u8(3, wrk).riemann_zeta(wrk, RoundingMode::None, cc);
        if th.is_zero() {
            return finish(z3, p, rm);
        }
        let sign = if th.is_negative() { -1i32 } else { 1 };
        let th_abs = th.abs();
        let zero = ExactNum::new(wrk);
        let q = tanh_sinh(
            |t, pw, r, c| t.clausen_cl2(pw, r, c),
            &zero,
            &th_abs,
            wrk,
            RoundingMode::None,
            cc,
        );
        match q {
            Some(v) => {
                let y = z3.sub(&v, wrk, RoundingMode::None);
                let _ = sign;
                finish(y, p, rm)
            }
            None => cat_nan(),
        }
    }

    /// Hurwitz \(\zeta(s,a)\) for integer \(s\ge 2\) and \(a>0\). `self` is \(s\).
    pub fn hurwitz_zeta(&self, a: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) || !finite(a) {
            return cat_nan();
        }
        let wrk = work_p(p);
        let s = match small_int(self, wrk, 2, 64) {
            Some(n) => n as u32,
            None => return cat_nan(),
        };
        let mut aa = a.clone();
        let mut extra = ExactNum::new(wrk);
        while matches!(aa.cmp(&ExactNum::from_u8(1, wrk)), Some(c) if c < 0) {
            if aa.is_zero() || (aa.is_int() && !aa.is_positive()) {
                return cat_nan();
            }
            extra = extra.add(
                &ExactNum::from_u8(1, wrk).div(
                    &aa.powsi(s as isize, wrk, RoundingMode::None),
                    wrk,
                    RoundingMode::None,
                ),
                wrk,
                RoundingMode::None,
            );
            aa = aa.add(&ExactNum::from_u8(1, wrk), wrk, RoundingMode::None);
        }
        finish(
            hurwitz_em(s, &aa, wrk, cc).add(&extra, wrk, RoundingMode::None),
            p,
            rm,
        )
    }

    /// Riemann \(\zeta(s)\) for integer \(s\ge 2\). `self` is \(s\).
    pub fn riemann_zeta(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let wrk = work_p(p);
        if let Some(s) = small_int(self, wrk, 2, 64) {
            if s % 2 == 0 {
                return finish(zeta_even((s / 2) as usize, wrk, cc), p, rm);
            }
        }
        self.hurwitz_zeta(&ExactNum::from_u8(1, p), p, rm, cc)
    }

    /// Polygamma \(\psi^{(n)}(\mathrm{self})\). \(n=0\) is digamma; \(n\ge1\) uses Hurwitz.
    pub fn polygamma(&self, n: usize, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) {
            return cat_nan();
        }
        if n == 0 {
            return self.digamma(p, rm, cc);
        }
        if n > 62 {
            return cat_nan();
        }
        let wrk = work_p(p);
        let z =
            ExactNum::from_u32((n + 1) as u32, wrk).hurwitz_zeta(self, wrk, RoundingMode::None, cc);
        let fact = factorial_n(n, wrk);
        let mut y = fact.mul(&z, wrk, RoundingMode::None);
        // ψ^{(n)} = (-1)^{n+1} n! ζ(n+1, z)
        if n % 2 == 0 {
            y = y.neg();
        }
        finish(y, p, rm)
    }

    /// Barnes \(G(\mathrm{self})\). Entire; zeros at non-positive integers.
    pub fn barnes_g(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) {
            return cat_nan();
        }
        let wrk = work_p(p);
        if let Some(n) = small_int(self, wrk, -64, 64) {
            if n <= 0 {
                return finish(ExactNum::new(wrk), p, rm);
            }
            let mut g = ExactNum::from_u8(1, wrk);
            for k in 2..n {
                g = g.mul(&factorial_n((k - 1) as usize, wrk), wrk, RoundingMode::None);
            }
            return finish(g, p, rm);
        }
        if !self.is_positive() {
            // Recur up to a positive argument: G(z) = G(z+1)/Γ(z).
            let zp = self.add(&ExactNum::from_u8(1, wrk), wrk, RoundingMode::None);
            let g1 = zp.barnes_g(wrk, RoundingMode::None, cc);
            let gz = self.gamma(wrk, RoundingMode::None, cc);
            return finish(g1.div(&gz, wrk, RoundingMode::None), p, rm);
        }
        // Shift up until z ≳ 16, Stirling for ln G(z+1), then divide by the Γ product.
        let mut z = self.clone();
        let mut ln_prod = ExactNum::new(wrk);
        let shift_to = ExactNum::from_u32((32 + p / 4) as u32, wrk);
        while matches!(z.cmp(&shift_to), Some(c) if c < 0) {
            ln_prod = ln_prod.add(
                &z.ln_gamma(wrk, RoundingMode::None, cc),
                wrk,
                RoundingMode::None,
            );
            z = z.add(&ExactNum::from_u8(1, wrk), wrk, RoundingMode::None);
        }
        // Now G(self) = G(z) / Π Γ, and z = self + n with n steps. Stirling for ln G(z):
        // ln G(z+1) = z ln Γ(z+1) − z(z+1)/2 + (z/2) ln(2π) − ln A + Σ B_{2k+2}/(4k(k+1) z^{2k})
        // We have G(z) = G((z-1)+1).
        let zm1 = z.sub(&ExactNum::from_u8(1, wrk), wrk, RoundingMode::None);
        let ln_g_z = barnes_ln_g_plus1(&zm1, wrk, cc);
        let ln_g_self = ln_g_z.sub(&ln_prod, wrk, RoundingMode::None);
        finish(ln_g_self.exp(wrk, RoundingMode::None, cc), p, rm)
    }

    /// Inverse Jacobi \(\mathrm{arcsn}(\mathrm{self}|m)=F(\mathrm{self}|m)\).
    pub fn jacobi_arcsn(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.elliptic_f(m, p, rm, cc)
    }

    /// Inverse Jacobi \(\mathrm{arccn}(x|m)=F(\sqrt{1-x^2}|m)\).
    pub fn jacobi_arccn(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) || !finite(m) {
            return cat_nan();
        }
        let wrk = work_p(p);
        let one = ExactNum::from_u8(1, wrk);
        let s = one
            .sub(
                &self.mul(self, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            )
            .sqrt(wrk, RoundingMode::None);
        finish(s.elliptic_f(m, wrk, RoundingMode::None, cc), p, rm)
    }

    /// Inverse Jacobi \(\mathrm{arcdn}(x|m)=F(\sqrt{(1-x^2)/m}|m)\) for \(m>0\).
    pub fn jacobi_arcdn(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) || !finite(m) || !m.is_positive() {
            return cat_nan();
        }
        let wrk = work_p(p);
        let one = ExactNum::from_u8(1, wrk);
        let s = one
            .sub(
                &self.mul(self, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            )
            .div(m, wrk, RoundingMode::None)
            .sqrt(wrk, RoundingMode::None);
        finish(s.elliptic_f(m, wrk, RoundingMode::None, cc), p, rm)
    }

    /// \({}_0F_1(;b;z)\) with `self = z`.
    pub fn hypergeom_0f1(&self, b: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.hypergeom_pfq(&[], core::slice::from_ref(b), p, rm, cc)
    }

    /// Kummer \({}_1F_1(a;b;z)\) with `self = z`.
    pub fn hypergeom_1f1(
        &self,
        a: &Self,
        b: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        self.hypergeom_pfq(
            core::slice::from_ref(a),
            core::slice::from_ref(b),
            p,
            rm,
            cc,
        )
    }

    /// Generalized \({}_pF_q(a_1,\ldots,a_p;b_1,\ldots,b_q;z)\) with `self = z`.
    ///
    /// Direct series. Divergent when \(p>q+1\) unless \(z=0\). Poles in any \(b_j\in\{0,-1,\ldots\}\).
    pub fn hypergeom_pfq(
        &self,
        a: &[Self],
        b: &[Self],
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        let _ = cc;
        if !finite(self) || a.iter().any(|x| !finite(x)) || b.iter().any(|x| !finite(x)) {
            return cat_nan();
        }
        if a.len() > b.len() + 1 && !self.is_zero() {
            return cat_nan();
        }
        let wrk = work_p(p);
        for bj in b {
            if let Some(n) = small_int(bj, wrk, -64, 0) {
                if n <= 0 {
                    return cat_nan();
                }
            }
        }
        let mut term = ExactNum::from_u8(1, wrk);
        let mut sum = term.clone();
        for k in 1..=series_cap(wrk) {
            let kf = ExactNum::from_u32(k as u32, wrk);
            let mut num = ExactNum::from_u8(1, wrk);
            for aj in a {
                num = num.mul(
                    &aj.add(
                        &ExactNum::from_u32((k - 1) as u32, wrk),
                        wrk,
                        RoundingMode::None,
                    ),
                    wrk,
                    RoundingMode::None,
                );
            }
            let mut den = kf;
            for bj in b {
                den = den.mul(
                    &bj.add(
                        &ExactNum::from_u32((k - 1) as u32, wrk),
                        wrk,
                        RoundingMode::None,
                    ),
                    wrk,
                    RoundingMode::None,
                );
            }
            term = term
                .mul(&num, wrk, RoundingMode::None)
                .div(&den, wrk, RoundingMode::None)
                .mul(self, wrk, RoundingMode::None);
            sum = sum.add(&term, wrk, RoundingMode::None);
            if term.is_zero()
                || term.abs().exponent().is_some_and(|e| {
                    (e as i64) + (p as i64) < sum.abs().exponent().unwrap_or(0) as i64
                })
            {
                break;
            }
        }
        finish(sum, p, rm)
    }

    /// Principal Lambert \(W_0(\mathrm{self})\). Domain \([-1/e,\infty)\).
    pub fn lambert_w0(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.lambert_w_branch(false, p, rm, cc)
    }

    /// Lambert \(W_{-1}(\mathrm{self})\). Domain \([-1/e,0)\).
    pub fn lambert_wm1(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.lambert_w_branch(true, p, rm, cc)
    }

    fn lambert_w_branch(&self, minus1: bool, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) {
            return cat_nan();
        }
        let wrk = work_p(p);
        let e = ExactNum::from_u8(1, wrk).exp(wrk, RoundingMode::None, cc);
        let inv_e = ExactNum::from_u8(1, wrk)
            .div(&e, wrk, RoundingMode::None)
            .neg();
        let gap = self.sub(&inv_e, wrk, RoundingMode::None);
        if gap
            .abs()
            .exponent()
            .is_some_and(|e| (e as i64) + ((p as i64) / 2) < 0)
            || self.cmp(&inv_e) == Some(0)
        {
            return finish(ExactNum::from_i8(-1, wrk), p, rm);
        }
        if matches!(self.cmp(&inv_e), Some(c) if c < 0) {
            return cat_nan();
        }
        if minus1 {
            if !self.is_negative() || self.is_zero() {
                return cat_nan();
            }
        }
        if self.is_zero() {
            return finish(ExactNum::new(wrk), p, rm);
        }
        let mut w = if minus1 {
            let two = ExactNum::from_u8(2, wrk);
            let d = two
                .mul(&e, wrk, RoundingMode::None)
                .mul(
                    &self.sub(&inv_e, wrk, RoundingMode::None),
                    wrk,
                    RoundingMode::None,
                )
                .sqrt(wrk, RoundingMode::None);
            ExactNum::from_i8(-1, wrk).sub(&d, wrk, RoundingMode::None)
        } else if matches!(self.cmp(&ExactNum::from_u8(1, wrk)), Some(c) if c > 0) {
            let l = self.ln(wrk, RoundingMode::None, cc);
            l.sub(&l.ln(wrk, RoundingMode::None, cc), wrk, RoundingMode::None)
        } else {
            self.clone()
        };
        let one = ExactNum::from_u8(1, wrk);
        let two = ExactNum::from_u8(2, wrk);
        for _ in 0..(8 + p / 16) {
            let ew = w.exp(wrk, RoundingMode::None, cc);
            let wew = w.mul(&ew, wrk, RoundingMode::None);
            let num = wew.sub(self, wrk, RoundingMode::None);
            let den = ew
                .mul(
                    &w.add(&one, wrk, RoundingMode::None),
                    wrk,
                    RoundingMode::None,
                )
                .sub(
                    &w.add(&two, wrk, RoundingMode::None)
                        .mul(&num, wrk, RoundingMode::None)
                        .div(
                            &two.mul(
                                &w.add(&one, wrk, RoundingMode::None),
                                wrk,
                                RoundingMode::None,
                            ),
                            wrk,
                            RoundingMode::None,
                        ),
                    wrk,
                    RoundingMode::None,
                );
            let nxt = w.sub(
                &num.div(&den, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            );
            let err = nxt.sub(&w, wrk, RoundingMode::None).abs();
            w = nxt;
            if err.is_zero() || err.exponent().is_some_and(|e| (e as i64) + (p as i64) < 0) {
                break;
            }
        }
        finish(w, p, rm)
    }

    /// Polylogarithm \(\mathrm{Li}_n(\mathrm{self})\) for integer \(n\ge 2\) and \(\lvert x\rvert\le 1\).
    pub fn polylog(&self, n: usize, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if !finite(self) || !(2..=64).contains(&n) {
            return cat_nan();
        }
        let wrk = work_p(p);
        let one = ExactNum::from_u8(1, wrk);
        if self.cmp(&one) == Some(0) {
            return ExactNum::from_u32(n as u32, wrk).riemann_zeta(p, rm, cc);
        }
        if matches!(self.abs().cmp(&one), Some(c) if c > 0) {
            return cat_nan();
        }
        if self.cmp(&one.neg()) == Some(0) {
            let z = ExactNum::from_u32(n as u32, wrk).riemann_zeta(wrk, RoundingMode::None, cc);
            let two = ExactNum::from_u8(2, wrk);
            let f = two.powsi(1 - (n as isize), wrk, RoundingMode::None).sub(
                &one,
                wrk,
                RoundingMode::None,
            );
            return finish(f.mul(&z, wrk, RoundingMode::None), p, rm);
        }
        let mut xk = self.clone();
        let mut sum = ExactNum::new(wrk);
        for k in 1..=series_cap(wrk).saturating_mul(4) {
            let den = ExactNum::from_u32(k as u32, wrk).powsi(n as isize, wrk, RoundingMode::None);
            let term = xk.div(&den, wrk, RoundingMode::None);
            sum = sum.add(&term, wrk, RoundingMode::None);
            if term
                .abs()
                .exponent()
                .is_some_and(|e| (e as i64) + (p as i64) < 0)
            {
                break;
            }
            xk = xk.mul(self, wrk, RoundingMode::None);
        }
        let _ = cc;
        finish(sum, p, rm)
    }
}

fn barnes_ln_g_plus1(z: &ExactNum, p: usize, cc: &mut Consts) -> ExactNum {
    // DLMF 5.17.5: ln G(z+1) ∼ (z²/2) ln z − (3/4) z² + (z/2) ln(2π) − (1/12) ln z
    //              + ζ'(−1) + Σ_{k≥1} B_{2k+2} / (4k(k+1) z^{2k}),  ζ'(−1)=1/12−ln A.
    let two = ExactNum::from_u8(2, p);
    let three = ExactNum::from_u8(3, p);
    let four = ExactNum::from_u8(4, p);
    let twelve = ExactNum::from_u32(12, p);
    let z2 = z.mul(z, p, RoundingMode::None);
    let lnz = z.ln(p, RoundingMode::None, cc);
    let pi = cc.pi(p, RoundingMode::None);
    let twopi = two.mul(&pi, p, RoundingMode::None);
    let t1 = z2
        .div(&two, p, RoundingMode::None)
        .mul(&lnz, p, RoundingMode::None);
    let t2 = three
        .mul(&z2, p, RoundingMode::None)
        .div(&four, p, RoundingMode::None);
    let t3 = z.div(&two, p, RoundingMode::None).mul(
        &twopi.ln(p, RoundingMode::None, cc),
        p,
        RoundingMode::None,
    );
    let t4 = lnz.div(&twelve, p, RoundingMode::None);
    let lna = ln_glaisher(p, cc);
    let zeta_pm1 = ExactNum::from_u8(1, p)
        .div(&twelve, p, RoundingMode::None)
        .sub(&lna, p, RoundingMode::None);
    let mut s = t1
        .sub(&t2, p, RoundingMode::None)
        .add(&t3, p, RoundingMode::None)
        .sub(&t4, p, RoundingMode::None)
        .add(&zeta_pm1, p, RoundingMode::None);
    let m = (12 + p / 12).min(48);
    let bs = even_bernoulli_ext(m + 1, p);
    let mut z2k = ExactNum::from_u8(1, p);
    for (k, b) in bs.iter().enumerate().skip(1) {
        z2k = z2k.mul(&z2, p, RoundingMode::None);
        let kf = ExactNum::from_u32(k as u32, p);
        let kp1 = ExactNum::from_u32((k + 1) as u32, p);
        let den = four
            .mul(&kf, p, RoundingMode::None)
            .mul(&kp1, p, RoundingMode::None)
            .mul(&z2k, p, RoundingMode::None);
        s = s.add(&b.div(&den, p, RoundingMode::None), p, RoundingMode::None);
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(a: &ExactNum, b: &ExactNum, p: usize) -> bool {
        let d = a.sub(b, p, RoundingMode::None).abs();
        d.is_zero() || d.exponent().is_some_and(|e| e < -((p as i32) / 4))
    }

    #[test]
    fn test_catalog_identities() {
        let p = 128;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new().unwrap();
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);
        let half = one.div(&ExactNum::from_u8(2, p), p, rm);

        let gi0 = zero.scorer_gi(p, rm, &mut cc);
        let bi0 = zero.bi(p, rm, &mut cc);
        assert!(near(&gi0, &bi0.div(&ExactNum::from_u8(3, p), p, rm), p));
        let hi0 = zero.scorer_hi(p, rm, &mut cc);
        assert!(near(
            &hi0,
            &bi0.mul(&ExactNum::from_u8(2, p), p, rm)
                .div(&ExactNum::from_u8(3, p), p, rm),
            p
        ));
        let x2 = ExactNum::from_u8(2, p);
        let gi2 = x2.scorer_gi(p, rm, &mut cc);
        let hi2 = x2.scorer_hi(p, rm, &mut cc);
        let bi2 = x2.bi(p, rm, &mut cc);
        assert!(near(&hi2.add(&gi2, p, rm), &bi2, p), "Hi+Gi=Bi at 2");

        let ber0 = zero.kelvin_ber(p, rm, &mut cc);
        assert!(near(&ber0, &one, p));
        let bei0 = zero.kelvin_bei(p, rm, &mut cc);
        assert!(bei0.is_zero() || near(&bei0, &zero, p));

        let cl2pi = cc.pi(p, rm).clausen_cl2(p, rm, &mut cc);
        assert!(cl2pi.is_zero() || near(&cl2pi, &zero, p), "Cl2(pi)");
        let z2 = ExactNum::from_u8(2, p).riemann_zeta(p, rm, &mut cc);
        let pi2_6 = cc
            .pi(p, rm)
            .mul(&cc.pi(p, rm), p, rm)
            .div(&ExactNum::from_u32(6, p), p, rm);
        assert!(near(&z2, &pi2_6, p), "zeta(2)");

        let z3 = ExactNum::from_u8(3, p).riemann_zeta(p, rm, &mut cc);
        let cl30 = zero.clausen_cl3(p, rm, &mut cc);
        assert!(near(&cl30, &z3, p));

        assert!(near(
            &ExactNum::from_u8(1, p).barnes_g(p, rm, &mut cc),
            &one,
            p
        ));
        assert!(near(
            &ExactNum::from_u8(2, p).barnes_g(p, rm, &mut cc),
            &one,
            p
        ));
        assert!(near(
            &ExactNum::from_u8(4, p).barnes_g(p, rm, &mut cc),
            &ExactNum::from_u8(2, p),
            p
        ));
        assert!(near(
            &ExactNum::from_u8(5, p).barnes_g(p, rm, &mut cc),
            &ExactNum::from_u32(12, p),
            p
        ));

        let psi0 = one.polygamma(0, p, rm, &mut cc);
        let g = cc.euler_gamma(p, rm).neg();
        assert!(near(&psi0, &g, p));

        let j = one.bessel_j(1, p, rm, &mut cc);
        let aj = one.anger_j(&one, p, rm, &mut cc);
        assert!(near(&aj, &j, p));

        let w0 = zero.lambert_w0(p, rm, &mut cc);
        assert!(w0.is_zero() || near(&w0, &zero, p));
        let inv_e = ExactNum::from_u8(1, p)
            .exp(p, rm, &mut cc)
            .reciprocal(p, rm)
            .neg();
        let wm = inv_e.lambert_w0(p, rm, &mut cc);
        assert!(near(&wm, &ExactNum::from_i8(-1, p), p));

        let e = one.exp(p, rm, &mut cc);
        let m = one.hypergeom_1f1(&half, &half, p, rm, &mut cc);
        assert!(near(&m, &e, p));

        let li2h = half.polylog(2, p, rm, &mut cc);
        // Li_2(1/2) = π²/12 − (ln 2)²/2
        let pi = cc.pi(p, rm);
        let ln2 = cc.ln_2(p, rm);
        let want = pi
            .mul(&pi, p, rm)
            .div(&ExactNum::from_u32(12, p), p, rm)
            .sub(
                &ln2.mul(&ln2, p, rm).div(&ExactNum::from_u8(2, p), p, rm),
                p,
                rm,
            );
        assert!(near(&li2h, &want, p));

        let lna = ln_glaisher(p, &mut cc);
        let lna_ref = ExactNum::parse(
            "0.2487544770337842625472529935761139760974",
            crate::Radix::Dec,
            4 * p,
            rm,
            &mut cc,
        );
        let dl = lna.sub(&lna_ref, 2 * p, RoundingMode::None).abs();
        let lb = if dl.is_zero() {
            i32::MAX
        } else {
            lna_ref.exponent().unwrap_or(0) - dl.exponent().unwrap_or(0)
        };
        assert!(lb >= 80, "ln A only {lb} bits");

        let mpar = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(4, p), p, rm);
        let u = half.elliptic_f(&mpar, p, rm, &mut cc);
        let sn = u.jacobi_sn(&mpar, p, rm, &mut cc);
        let back = sn.jacobi_arcsn(&mpar, p, rm, &mut cc);
        assert!(near(&back, &u, p));
    }
}
