//! Real distribution PDF / CDF / PMF kernels on [`ExactNum`].

use crate::defs::WORD_BIT_SIZE;
use crate::Consts;
use crate::Error;
use crate::ExactNum;
use crate::RoundingMode;
use crate::NAN;

fn dist_nan() -> ExactNum {
    ExactNum::nan(Some(Error::InvalidArgument))
}

fn work_p(p: usize) -> usize {
    p.saturating_add(WORD_BIT_SIZE)
}

fn finite_pos(x: &ExactNum) -> bool {
    x.is_positive() && !x.is_inf()
}

fn nn_int(x: &ExactNum) -> bool {
    x.is_int() && !x.is_negative()
}

impl ExactNum {
    /// Standard-form normal density
    /// \(\varphi(x;\mu,\sigma)=\exp(-(x-\mu)^2/(2\sigma^2))/(\sigma\sqrt{2\pi})\).
    ///
    /// `sigma ≤ 0` or a non-finite parameter is `NaN`.
    ///
    /// # Precision
    ///
    /// - Algorithm: `exp` / `sqrt` / `π` at `p + WORD_BIT_SIZE`, then one rounding to `p`.
    /// - Bound: working precision `p + WORD_BIT_SIZE`.
    /// - MPFR oracle: no (composite of existing leaves).
    pub fn normal_pdf(
        &self,
        mu: &Self,
        sigma: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        if self.is_nan() || mu.is_nan() || sigma.is_nan() {
            return NAN;
        }
        if !finite_pos(sigma) {
            return dist_nan();
        }
        if self.is_inf() {
            return ExactNum::from_u8(0, p);
        }
        let pw = work_p(p);
        let none = RoundingMode::None;
        let z = self.sub(mu, pw, none);
        let two = ExactNum::from_u8(2, pw);
        let sig2 = sigma.mul(sigma, pw, none);
        let expo = z
            .mul(&z, pw, none)
            .div(&two.mul(&sig2, pw, none), pw, none)
            .neg();
        let num = expo.exp(pw, none, cc);
        let two_pi = two.mul(&cc.pi(pw, none), pw, none);
        let den = sigma.mul(&two_pi.sqrt(pw, none), pw, none);
        num.div(&den, p, rm)
    }

    /// Normal CDF \(\Phi(x;\mu,\sigma)=(1+\mathrm{erf}((x-\mu)/(\sigma\sqrt{2})))/2\).
    ///
    /// `sigma ≤ 0` is `NaN`.
    ///
    /// # Precision
    ///
    /// - Algorithm: existing `erf` at `p + WORD_BIT_SIZE`.
    /// - Bound: working precision `p + WORD_BIT_SIZE`.
    /// - MPFR oracle: no (composite of `erf`).
    pub fn normal_cdf(
        &self,
        mu: &Self,
        sigma: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        if self.is_nan() || mu.is_nan() || sigma.is_nan() {
            return NAN;
        }
        if !finite_pos(sigma) {
            return dist_nan();
        }
        let pw = work_p(p);
        let none = RoundingMode::None;
        let two = ExactNum::from_u8(2, pw);
        let z = self
            .sub(mu, pw, none)
            .div(&sigma.mul(&two.sqrt(pw, none), pw, none), pw, none);
        let one = ExactNum::from_u8(1, pw);
        one.add(&z.erf(pw, none, cc), pw, none).div(&two, p, rm)
    }

    /// Gamma density on the scale parameterization
    /// \(x^{\alpha-1}e^{-x/\beta}/(\beta^\alpha\Gamma(\alpha))\).
    ///
    /// Requires `x ≥ 0`, `alpha > 0`, `beta > 0`.
    ///
    /// # Precision
    ///
    /// - Algorithm: `pow` / `exp` / `gamma` at `p + WORD_BIT_SIZE`.
    /// - Bound: working precision `p + WORD_BIT_SIZE`.
    /// - MPFR oracle: no.
    pub fn gamma_pdf(
        &self,
        alpha: &Self,
        beta: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        if self.is_nan() || alpha.is_nan() || beta.is_nan() {
            return NAN;
        }
        if self.is_negative() || !finite_pos(alpha) || !finite_pos(beta) {
            return dist_nan();
        }
        let pw = work_p(p);
        let none = RoundingMode::None;
        let one = ExactNum::from_u8(1, pw);
        let am1 = alpha.sub(&one, pw, none);
        let xb = self.div(beta, pw, none);
        let num = self
            .pow(&am1, pw, none, cc)
            .mul(&xb.neg().exp(pw, none, cc), pw, none);
        let den = beta
            .pow(alpha, pw, none, cc)
            .mul(&alpha.gamma(pw, none, cc), pw, none);
        num.div(&den, p, rm)
    }

    /// Beta density \(x^{\alpha-1}(1-x)^{\beta-1}/B(\alpha,\beta)\) with
    /// \(B(\alpha,\beta)=\Gamma(\alpha)\Gamma(\beta)/\Gamma(\alpha+\beta)\).
    ///
    /// Requires `x ∈ [0, 1]`, `alpha > 0`, `beta > 0`.
    ///
    /// # Precision
    ///
    /// - Algorithm: `pow` / `gamma` at `p + WORD_BIT_SIZE`.
    /// - Bound: working precision `p + WORD_BIT_SIZE`.
    /// - MPFR oracle: no.
    pub fn beta_pdf(
        &self,
        alpha: &Self,
        beta: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        if self.is_nan() || alpha.is_nan() || beta.is_nan() {
            return NAN;
        }
        if self.is_negative()
            || matches!(self.cmp(&ExactNum::from_u8(1, p)), Some(c) if c > 0)
            || !finite_pos(alpha)
            || !finite_pos(beta)
        {
            return dist_nan();
        }
        let pw = work_p(p);
        let none = RoundingMode::None;
        let one = ExactNum::from_u8(1, pw);
        let am1 = alpha.sub(&one, pw, none);
        let bm1 = beta.sub(&one, pw, none);
        let num = self.pow(&am1, pw, none, cc).mul(
            &one.sub(self, pw, none).pow(&bm1, pw, none, cc),
            pw,
            none,
        );
        let bfn = alpha
            .gamma(pw, none, cc)
            .mul(&beta.gamma(pw, none, cc), pw, none)
            .div(&alpha.add(beta, pw, none).gamma(pw, none, cc), pw, none);
        num.div(&bfn, p, rm)
    }

    /// Poisson PMF \(\lambda^k e^{-\lambda}/k!\) for a non-negative integer `self` \(= k\).
    ///
    /// # Precision
    ///
    /// - Algorithm: `pow` / `exp` / \(\Gamma(k+1)\) at `p + WORD_BIT_SIZE`.
    /// - Bound: working precision `p + WORD_BIT_SIZE`.
    /// - MPFR oracle: no.
    pub fn poisson_pmf(&self, lambda: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() || lambda.is_nan() {
            return NAN;
        }
        if !nn_int(self) || lambda.is_negative() {
            return dist_nan();
        }
        let pw = work_p(p);
        let none = RoundingMode::None;
        let one = ExactNum::from_u8(1, pw);
        let kf = self.add(&one, pw, none).gamma(pw, none, cc);
        lambda
            .pow(self, pw, none, cc)
            .mul(&lambda.neg().exp(pw, none, cc), pw, none)
            .div(&kf, p, rm)
    }

    /// Binomial PMF \(\binom{n}{k} \mathrm{prob}^k (1-\mathrm{prob})^{n-k}\).
    ///
    /// `self` is \(k\). Requires non-negative integers `k ≤ n` and `prob ∈ [0, 1]`.
    ///
    /// # Precision
    ///
    /// - Algorithm: multiplicative binomial coefficient, then `pow`.
    /// - Bound: working precision `p + WORD_BIT_SIZE`.
    /// - MPFR oracle: no.
    pub fn binomial_pmf(
        &self,
        n: &Self,
        prob: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        if self.is_nan() || n.is_nan() || prob.is_nan() {
            return NAN;
        }
        if !nn_int(self)
            || !nn_int(n)
            || matches!(self.cmp(n), Some(c) if c > 0)
            || prob.is_negative()
            || matches!(prob.cmp(&ExactNum::from_u8(1, p)), Some(c) if c > 0)
        {
            return dist_nan();
        }
        let pw = work_p(p);
        let none = RoundingMode::None;
        let c = binom_mul(n, self, pw);
        let q = ExactNum::from_u8(1, pw).sub(prob, pw, none);
        let nmk = n.sub(self, pw, none);
        c.mul(&prob.pow(self, pw, none, cc), pw, none)
            .mul(&q.pow(&nmk, pw, none, cc), pw, none)
            .set_prec_val(p, rm)
    }

    /// Chi-squared CDF \(P(k/2, x/2)=\gamma(k/2, x/2)/\Gamma(k/2)\).
    ///
    /// `self` is \(x\). Requires `x ≥ 0` and `k > 0`.
    ///
    /// # Precision
    ///
    /// - Algorithm: lower `gammainc` over `gamma` at `p + WORD_BIT_SIZE`.
    /// - Bound: working precision `p + WORD_BIT_SIZE`.
    /// - MPFR oracle: no.
    pub fn chi_squared_cdf(&self, k: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() || k.is_nan() {
            return NAN;
        }
        if self.is_negative() || !finite_pos(k) {
            return dist_nan();
        }
        let pw = work_p(p);
        let none = RoundingMode::None;
        let two = ExactNum::from_u8(2, pw);
        let s = k.div(&two, pw, none);
        let xh = self.div(&two, pw, none);
        s.gammainc(&xh, pw, none, cc)
            .div(&s.gamma(pw, none, cc), p, rm)
    }

    /// Student-\(t\) density via \(\Gamma\):
    /// \(\Gamma((\nu+1)/2)/(\sqrt{\nu\pi}\,\Gamma(\nu/2))\,(1+x^2/\nu)^{-(\nu+1)/2}\).
    ///
    /// `self` is \(x\). Requires `nu > 0`.
    ///
    /// # Precision
    ///
    /// - Algorithm: `gamma` / `pow` / `sqrt` at `p + WORD_BIT_SIZE`.
    /// - Bound: working precision `p + WORD_BIT_SIZE`.
    /// - MPFR oracle: no.
    pub fn student_t_pdf(&self, nu: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() || nu.is_nan() {
            return NAN;
        }
        if !finite_pos(nu) {
            return dist_nan();
        }
        let pw = work_p(p);
        let none = RoundingMode::None;
        let one = ExactNum::from_u8(1, pw);
        let two = ExactNum::from_u8(2, pw);
        let np1 = nu.add(&one, pw, none);
        let half_np1 = np1.div(&two, pw, none);
        let half_n = nu.div(&two, pw, none);
        let pref = half_np1.gamma(pw, none, cc).div(
            &nu.mul(&cc.pi(pw, none), pw, none).sqrt(pw, none).mul(
                &half_n.gamma(pw, none, cc),
                pw,
                none,
            ),
            pw,
            none,
        );
        let body = one
            .add(&self.mul(self, pw, none).div(nu, pw, none), pw, none)
            .pow(&half_np1.neg(), pw, none, cc);
        pref.mul(&body, p, rm)
    }
}

fn binom_mul(n: &ExactNum, k: &ExactNum, pw: usize) -> ExactNum {
    let none = RoundingMode::None;
    let one = ExactNum::from_u8(1, pw);
    if k.is_zero() {
        return one;
    }
    let mut i = one.clone();
    let mut c = one.clone();
    let nmk = n.sub(k, pw, none);
    loop {
        let term = nmk.add(&i, pw, none);
        c = c.mul(&term, pw, none).div(&i, pw, none);
        if i.cmp(k) == Some(0) {
            return c;
        }
        i = i.add(&one, pw, none);
        if i.cmp(k) == Some(1) {
            return c;
        }
    }
}

impl ExactNum {
    fn set_prec_val(mut self, p: usize, rm: RoundingMode) -> Self {
        let _ = self.set_precision(p, rm);
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dist_plan_golds() {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new().unwrap();
        let zero = ExactNum::from_u8(0, p);
        let one = ExactNum::from_u8(1, p);
        let two = ExactNum::from_u8(2, p);

        let np = zero.normal_pdf(&zero, &one, p, rm, &mut cc);
        let two_pi = two.mul(&cc.pi(p, rm), p, rm);
        let want_np = two_pi.sqrt(p, rm).reciprocal(p, rm);
        assert_eq!(np.cmp(&want_np), Some(0));

        let nc = zero.normal_cdf(&zero, &one, p, rm, &mut cc);
        let half = one.div(&two, p, rm);
        assert_eq!(nc.cmp(&half), Some(0));

        let gp = one.gamma_pdf(&one, &one, p, rm, &mut cc);
        let em1 = one.neg().exp(p, rm, &mut cc);
        assert_eq!(gp.cmp(&em1), Some(0));

        let po = zero.poisson_pmf(&one, p, rm, &mut cc);
        assert_eq!(po.cmp(&em1), Some(0));

        // χ²(2) CDF at the 95% table quantile −2 ln(0.05) = 2 ln 20 is exactly 19/20.
        let twenty = ExactNum::from_u8(20, p);
        let x95 = two.mul(&twenty.ln(p, rm, &mut cc), p, rm);
        let chi = x95.chi_squared_cdf(&two, p, rm, &mut cc);
        let table = ExactNum::from_u8(19, p).div(&ExactNum::from_u8(20, p), p, rm);
        assert_eq!(chi.cmp(&table), Some(0));

        assert!(one.normal_pdf(&zero, &zero, p, rm, &mut cc).is_nan());
        assert!(one.neg().poisson_pmf(&one, p, rm, &mut cc).is_nan());
    }
}
