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
}
