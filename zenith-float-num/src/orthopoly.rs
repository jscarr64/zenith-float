//! Classical orthogonal polynomials on [`ExactNum`] via three-term recurrences.

use crate::defs::WORD_BIT_SIZE;
use crate::Error;
use crate::ExactNum;
use crate::RoundingMode;

/// Maximum degree for every family in this module. Larger `n` is `NaN`.
pub const ORTHOPOLY_N_MAX: usize = 256;

fn work_p(p: usize) -> usize {
    p.saturating_add(WORD_BIT_SIZE)
}

fn op_nan() -> ExactNum {
    ExactNum::nan(Some(Error::InvalidArgument))
}

fn finite(x: &ExactNum) -> bool {
    !x.is_nan() && !x.is_inf()
}

fn n_ok(n: usize) -> bool {
    n <= ORTHOPOLY_N_MAX
}

fn finish(mut y: ExactNum, p: usize, rm: RoundingMode) -> ExactNum {
    let _ = y.set_precision(p, rm);
    y
}

impl ExactNum {
    /// Probabilist's Hermite polynomial `He_n(self)`.
    ///
    /// `He_0 = 1`, `He_1 = x`, `He_{n+1} = x He_n − n He_{n−1}`.
    /// `n > ORTHOPOLY_N_MAX` or non-finite `self` is `NaN`.
    pub fn hermite_he(&self, n: usize, p: usize, rm: RoundingMode) -> Self {
        if !n_ok(n) || !finite(self) {
            return op_nan();
        }
        let wrk = work_p(p);
        if n == 0 {
            return finish(ExactNum::from_u8(1, wrk), p, rm);
        }
        if n == 1 {
            return finish(self.clone(), p, rm);
        }
        let mut prev = ExactNum::from_u8(1, wrk);
        let mut cur = self.clone();
        let _ = cur.set_precision(wrk, RoundingMode::None);
        for k in 1..n {
            let kf = ExactNum::from_u32(k as u32, wrk);
            let next = self.mul(&cur, wrk, RoundingMode::None).sub(
                &kf.mul(&prev, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            );
            prev = cur;
            cur = next;
        }
        finish(cur, p, rm)
    }

    /// Physicist's Hermite polynomial `H_n(self)`.
    ///
    /// `H_0 = 1`, `H_1 = 2x`, `H_{n+1} = 2x H_n − 2n H_{n−1}`.
    /// `n > ORTHOPOLY_N_MAX` or non-finite `self` is `NaN`.
    pub fn hermite_h(&self, n: usize, p: usize, rm: RoundingMode) -> Self {
        if !n_ok(n) || !finite(self) {
            return op_nan();
        }
        let wrk = work_p(p);
        let two = ExactNum::from_u8(2, wrk);
        if n == 0 {
            return finish(ExactNum::from_u8(1, wrk), p, rm);
        }
        if n == 1 {
            return finish(two.mul(self, wrk, RoundingMode::None), p, rm);
        }
        let mut prev = ExactNum::from_u8(1, wrk);
        let mut cur = two.mul(self, wrk, RoundingMode::None);
        for k in 1..n {
            let two_k = ExactNum::from_u32((2 * k) as u32, wrk);
            let next = two
                .mul(self, wrk, RoundingMode::None)
                .mul(&cur, wrk, RoundingMode::None)
                .sub(
                    &two_k.mul(&prev, wrk, RoundingMode::None),
                    wrk,
                    RoundingMode::None,
                );
            prev = cur;
            cur = next;
        }
        finish(cur, p, rm)
    }

    /// Laguerre polynomial `L_n(self)`.
    ///
    /// `L_0 = 1`, `L_1 = 1 − x`,
    /// `L_{n+1} = ((2n+1−x) L_n − n L_{n−1}) / (n+1)`.
    /// `n > ORTHOPOLY_N_MAX` or non-finite `self` is `NaN`.
    pub fn laguerre(&self, n: usize, p: usize, rm: RoundingMode) -> Self {
        self.gen_laguerre(n, &ExactNum::new(p), p, rm)
    }

    /// Generalized Laguerre `L_n^{(α)}(self)`.
    ///
    /// `L_0^{(α)} = 1`, `L_1^{(α)} = 1+α−x`,
    /// `L_{n+1}^{(α)} = (((2n+1+α−x) L_n − (n+α) L_{n−1}) / (n+1)`.
    /// `n > ORTHOPOLY_N_MAX` or a non-finite argument is `NaN`.
    pub fn gen_laguerre(&self, n: usize, alpha: &Self, p: usize, rm: RoundingMode) -> Self {
        if !n_ok(n) || !finite(self) || !finite(alpha) {
            return op_nan();
        }
        let wrk = work_p(p);
        let one = ExactNum::from_u8(1, wrk);
        if n == 0 {
            return finish(one, p, rm);
        }
        if n == 1 {
            return finish(
                one.add(alpha, wrk, RoundingMode::None)
                    .sub(self, wrk, RoundingMode::None),
                p,
                rm,
            );
        }
        let mut prev = one.clone();
        let mut cur = one
            .add(alpha, wrk, RoundingMode::None)
            .sub(self, wrk, RoundingMode::None);
        for k in 1..n {
            let kf = ExactNum::from_u32(k as u32, wrk);
            let two_k_1 = ExactNum::from_u32((2 * k + 1) as u32, wrk);
            let coeff =
                two_k_1
                    .add(alpha, wrk, RoundingMode::None)
                    .sub(self, wrk, RoundingMode::None);
            let k_a = kf.add(alpha, wrk, RoundingMode::None);
            let num = coeff.mul(&cur, wrk, RoundingMode::None).sub(
                &k_a.mul(&prev, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            );
            let den = ExactNum::from_u32((k + 1) as u32, wrk);
            let next = num.div(&den, wrk, RoundingMode::None);
            prev = cur;
            cur = next;
        }
        finish(cur, p, rm)
    }

    /// Chebyshev polynomial of the first kind `T_n(self)`.
    ///
    /// `T_0 = 1`, `T_1 = x`, `T_{n+1} = 2x T_n − T_{n−1}`.
    /// `n > ORTHOPOLY_N_MAX` or non-finite `self` is `NaN`.
    pub fn chebyshev_t(&self, n: usize, p: usize, rm: RoundingMode) -> Self {
        if !n_ok(n) || !finite(self) {
            return op_nan();
        }
        let wrk = work_p(p);
        if n == 0 {
            return finish(ExactNum::from_u8(1, wrk), p, rm);
        }
        if n == 1 {
            return finish(self.clone(), p, rm);
        }
        let two = ExactNum::from_u8(2, wrk);
        let mut prev = ExactNum::from_u8(1, wrk);
        let mut cur = self.clone();
        let _ = cur.set_precision(wrk, RoundingMode::None);
        for _ in 1..n {
            let next = two
                .mul(self, wrk, RoundingMode::None)
                .mul(&cur, wrk, RoundingMode::None)
                .sub(&prev, wrk, RoundingMode::None);
            prev = cur;
            cur = next;
        }
        finish(cur, p, rm)
    }

    /// Chebyshev polynomial of the second kind `U_n(self)`.
    ///
    /// `U_0 = 1`, `U_1 = 2x`, `U_{n+1} = 2x U_n − U_{n−1}`.
    /// `n > ORTHOPOLY_N_MAX` or non-finite `self` is `NaN`.
    pub fn chebyshev_u(&self, n: usize, p: usize, rm: RoundingMode) -> Self {
        if !n_ok(n) || !finite(self) {
            return op_nan();
        }
        let wrk = work_p(p);
        let two = ExactNum::from_u8(2, wrk);
        if n == 0 {
            return finish(ExactNum::from_u8(1, wrk), p, rm);
        }
        if n == 1 {
            return finish(two.mul(self, wrk, RoundingMode::None), p, rm);
        }
        let mut prev = ExactNum::from_u8(1, wrk);
        let mut cur = two.mul(self, wrk, RoundingMode::None);
        for _ in 1..n {
            let next = two
                .mul(self, wrk, RoundingMode::None)
                .mul(&cur, wrk, RoundingMode::None)
                .sub(&prev, wrk, RoundingMode::None);
            prev = cur;
            cur = next;
        }
        finish(cur, p, rm)
    }

    /// Gegenbauer (ultraspherical) polynomial `C_n^{(λ)}(self)`.
    ///
    /// `C_0 = 1`, `C_1 = 2λ x`,
    /// `C_{n+1} = (2(n+λ) x C_n − (n+2λ−1) C_{n−1}) / (n+1)`.
    /// Standard `C_n^{(λ)}`: `C_2^{(1)} = 4x² − 1 = U_2`, and
    /// `C_n^{(1/2)} = P_n` (Legendre). `n > ORTHOPOLY_N_MAX` or a
    /// non-finite argument is `NaN`.
    pub fn gegenbauer(&self, n: usize, lambda: &Self, p: usize, rm: RoundingMode) -> Self {
        if !n_ok(n) || !finite(self) || !finite(lambda) {
            return op_nan();
        }
        let wrk = work_p(p);
        let one = ExactNum::from_u8(1, wrk);
        if n == 0 {
            return finish(one, p, rm);
        }
        let two = ExactNum::from_u8(2, wrk);
        if n == 1 {
            return finish(
                two.mul(lambda, wrk, RoundingMode::None)
                    .mul(self, wrk, RoundingMode::None),
                p,
                rm,
            );
        }
        let mut prev = one.clone();
        let mut cur = two
            .mul(lambda, wrk, RoundingMode::None)
            .mul(self, wrk, RoundingMode::None);
        for k in 1..n {
            let kf = ExactNum::from_u32(k as u32, wrk);
            let two_n_l = two.mul(
                &kf.add(lambda, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            );
            let n_2l_1 = kf
                .add(
                    &two.mul(lambda, wrk, RoundingMode::None),
                    wrk,
                    RoundingMode::None,
                )
                .sub(&one, wrk, RoundingMode::None);
            let num = two_n_l
                .mul(self, wrk, RoundingMode::None)
                .mul(&cur, wrk, RoundingMode::None)
                .sub(
                    &n_2l_1.mul(&prev, wrk, RoundingMode::None),
                    wrk,
                    RoundingMode::None,
                );
            let den = ExactNum::from_u32((k + 1) as u32, wrk);
            let next = num.div(&den, wrk, RoundingMode::None);
            prev = cur;
            cur = next;
        }
        finish(cur, p, rm)
    }

    /// Jacobi polynomial `P_n^{(α,β)}(self)`.
    ///
    /// `P_0 = 1`, `P_1 = ((α−β)+(α+β+2)x)/2`, and for `n ≥ 1`
    /// `2(n+1)(n+α+β+1)(2n+α+β) P_{n+1}`
    /// `= [(2n+α+β+1)(α²−β²)+(2n+α+β)(2n+α+β+1)(2n+α+β+2)x] P_n`
    /// `− 2(n+α)(n+β)(2n+α+β+2) P_{n−1}`.
    /// `α=β=0` is Legendre `P_n`. `n > ORTHOPOLY_N_MAX` or a non-finite
    /// argument is `NaN`.
    pub fn jacobi_p(
        &self,
        n: usize,
        alpha: &Self,
        beta: &Self,
        p: usize,
        rm: RoundingMode,
    ) -> Self {
        if !n_ok(n) || !finite(self) || !finite(alpha) || !finite(beta) {
            return op_nan();
        }
        let wrk = work_p(p);
        let one = ExactNum::from_u8(1, wrk);
        let two = ExactNum::from_u8(2, wrk);
        if n == 0 {
            return finish(one, p, rm);
        }
        let apb = alpha.add(beta, wrk, RoundingMode::None);
        let amb = alpha.sub(beta, wrk, RoundingMode::None);
        let p1 = amb
            .add(
                &apb.add(&two, wrk, RoundingMode::None)
                    .mul(self, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            )
            .div(&two, wrk, RoundingMode::None);
        if n == 1 {
            return finish(p1, p, rm);
        }
        let mut prev = one.clone();
        let mut cur = p1;
        for k in 1..n {
            let kf = ExactNum::from_u32(k as u32, wrk);
            let k1 = ExactNum::from_u32((k + 1) as u32, wrk);
            let two_k_apb =
                ExactNum::from_u32((2 * k) as u32, wrk).add(&apb, wrk, RoundingMode::None);
            let two_k_apb_1 = two_k_apb.add(&one, wrk, RoundingMode::None);
            let two_k_apb_2 = two_k_apb.add(&two, wrk, RoundingMode::None);
            let a_den = two
                .mul(&k1, wrk, RoundingMode::None)
                .mul(
                    &kf.add(&apb, wrk, RoundingMode::None)
                        .add(&one, wrk, RoundingMode::None),
                    wrk,
                    RoundingMode::None,
                )
                .mul(&two_k_apb, wrk, RoundingMode::None);
            let a2mb2 = alpha.mul(alpha, wrk, RoundingMode::None).sub(
                &beta.mul(beta, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            );
            let b_num = two_k_apb_1.mul(&a2mb2, wrk, RoundingMode::None);
            let c_num = two_k_apb.mul(&two_k_apb_1, wrk, RoundingMode::None).mul(
                &two_k_apb_2,
                wrk,
                RoundingMode::None,
            );
            let d_num = two
                .mul(
                    &kf.add(alpha, wrk, RoundingMode::None),
                    wrk,
                    RoundingMode::None,
                )
                .mul(
                    &kf.add(beta, wrk, RoundingMode::None),
                    wrk,
                    RoundingMode::None,
                )
                .mul(&two_k_apb_2, wrk, RoundingMode::None);
            let lin = b_num.add(
                &c_num.mul(self, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            );
            let num = lin.mul(&cur, wrk, RoundingMode::None).sub(
                &d_num.mul(&prev, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            );
            let next = num.div(&a_den, wrk, RoundingMode::None);
            prev = cur;
            cur = next;
        }
        finish(cur, p, rm)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Consts;

    fn gold_p() -> (usize, RoundingMode) {
        (256, RoundingMode::ToEven)
    }

    #[test]
    fn orthopoly_he4_laguerre_t5_gegenbauer_recurrence() {
        let (p, rm) = gold_p();
        let mut cc = Consts::new().expect("consts");
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);
        let two = ExactNum::from_u8(2, p);
        let three = ExactNum::from_u8(3, p);
        let half = one.div(&two, p, rm);

        assert_eq!(zero.hermite_he(4, p, rm).cmp(&three), Some(0));
        assert_eq!(zero.laguerre(3, p, rm).cmp(&one), Some(0));
        assert_eq!(zero.gen_laguerre(3, &zero, p, rm).cmp(&one), Some(0));

        let pi = cc.pi(p, rm);
        let five = ExactNum::from_u8(5, p);
        let c = pi.div(&five, p, rm).cos(p, rm, &mut cc);
        let t5 = c.chebyshev_t(5, p, rm);
        let neg_one = one.neg();
        assert_eq!(t5.cmp(&neg_one), Some(0));

        // Standard C_2^{(1)} = 4x² − 1 = U_2. The plan's "3x²−1 at λ=1"
        // is twice Legendre P_2 = 2 C_2^{(1/2)}.
        let g1 = half.gegenbauer(2, &one, p, rm);
        let u2 = half.chebyshev_u(2, p, rm);
        let four_x2_m1 = ExactNum::from_u8(4, p)
            .mul(&half, p, rm)
            .mul(&half, p, rm)
            .sub(&one, p, rm);
        assert_eq!(g1.cmp(&four_x2_m1), Some(0));
        assert_eq!(g1.cmp(&u2), Some(0));
        let lam_half = half.clone();
        let g_leg = half.gegenbauer(2, &lam_half, p, rm);
        let three_x2_m1 = three.mul(&half, p, rm).mul(&half, p, rm).sub(&one, p, rm);
        let two_p2 = g_leg.mul(&two, p, rm);
        assert_eq!(two_p2.cmp(&three_x2_m1), Some(0));

        let t6 = half.chebyshev_t(6, p, rm);
        let t5h = half.chebyshev_t(5, p, rm);
        let t4h = half.chebyshev_t(4, p, rm);
        let rec = two.mul(&half, p, rm).mul(&t5h, p, rm).sub(&t4h, p, rm);
        assert_eq!(t6.cmp(&rec), Some(0));

        assert!(zero.hermite_he(ORTHOPOLY_N_MAX + 1, p, rm).is_nan());
        assert_eq!(
            zero.hermite_h(4, p, rm).cmp(&ExactNum::from_u8(12, p)),
            Some(0)
        );
    }

    fn near(a: &ExactNum, b: &ExactNum, p: usize) -> bool {
        let d = a.sub(b, p, RoundingMode::None).abs();
        if d.is_zero() {
            return true;
        }
        let ae = a.exponent().unwrap_or(0);
        d.exponent().is_some_and(|e| e - ae < -((p as i32) / 4))
    }

    #[test]
    fn jacobi_p_legendre_and_identities() {
        let (p, rm) = gold_p();
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);
        let two = ExactNum::from_u8(2, p);
        let three = ExactNum::from_u8(3, p);
        let half = one.div(&two, p, rm);
        let three_tenths = three.div(&ExactNum::from_u8(10, p), p, rm);

        // P_2^{(1,1)}(1/2) = (3/4)(5/4 − 1) = 3/16.
        let p2 = half.jacobi_p(2, &one, &one, p, rm);
        let want = three.div(&ExactNum::from_u8(16, p), p, rm);
        assert_eq!(p2.cmp(&want), Some(0));

        // α = β = 0 is Legendre P_n.
        let j5 = three_tenths.jacobi_p(5, &zero, &zero, p, rm);
        let l5 = three_tenths.legendre_p(5, p, rm);
        assert!(near(&j5, &l5, p), "P_5^{{(0,0)}}(3/10) vs Legendre");

        // P_n^{(α,β)}(−x) = (−1)^n P_n^{(β,α)}(x).
        let alpha = half.clone();
        let beta = one.div(&three, p, rm);
        let x = three_tenths.clone();
        let nx = x.neg();
        let lhs = nx.jacobi_p(6, &alpha, &beta, p, rm);
        let rhs = x.jacobi_p(6, &beta, &alpha, p, rm);
        assert!(near(&lhs, &rhs, p), "even-n reflection");
        let lhs7 = nx.jacobi_p(7, &alpha, &beta, p, rm);
        let rhs7 = x.jacobi_p(7, &beta, &alpha, p, rm).neg();
        assert!(near(&lhs7, &rhs7, p), "odd-n reflection");

        assert!(zero
            .jacobi_p(ORTHOPOLY_N_MAX + 1, &zero, &zero, p, rm)
            .is_nan());
        assert_eq!(zero.jacobi_p(0, &one, &two, p, rm).cmp(&one), Some(0));
    }
}
