//! Complex Airy \(\mathrm{Ai}\) / \(\mathrm{Bi}\).
//!
//! Software-limb [`ExactComplex`] arithmetic. The real line uses the real
//! kernels. Off the real axis the Taylor pair is used for small \(\lvert z\rvert\);
//! large \(\lvert z\rvert\) uses the decaying/growing asymptotic in
//! \(\lvert\mathrm{arg}\,z\rvert\le 2\pi/3\) (Ai) or \(\pi/3\) (Bi), otherwise
//! the \(\omega\)-connection. The real series is not evaluated at \(\lvert z\rvert\).

use crate::common::util::round_p;
use crate::complex_special::half_c;
use crate::complex_special::neg_c;
use crate::complex_special::pi_c;
use crate::complex_special::series_term_cap;
use crate::complex_special::term_negligible;
use crate::complex_special::two_c;
use crate::complex_special::ziv_complex;
use crate::Consts;
use crate::ExactComplex;
use crate::ExactNum;
use crate::RoundingMode;

/// \(\lvert z\rvert\) below this uses the Taylor pair \((f,g)\).
const AIRY_SERIES_THRESHOLD: u32 = 8;

fn abs_below(z: &ExactComplex, bound: u32, p: usize) -> bool {
    let a = z.abs(p, RoundingMode::None);
    let b = ExactNum::from_u32(bound, p);
    matches!(a.cmp(&b), Some(c) if c < 0)
}

fn arg_abs_le(z: &ExactComplex, bound: &ExactNum, p: usize, cc: &mut Consts) -> bool {
    let a = z.arg(p, RoundingMode::None, cc).abs();
    matches!(a.cmp(bound), Some(c) if c <= 0)
}

fn omega_pair(p: usize) -> (ExactComplex, ExactComplex) {
    let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, RoundingMode::None);
    let s3h = ExactNum::from_u8(3, p)
        .sqrt(p, RoundingMode::None)
        .mul(&half, p, RoundingMode::None);
    let re = half.neg();
    (
        ExactComplex::new(re.clone(), s3h.clone()),
        ExactComplex::new(re, s3h.neg()),
    )
}

fn exp_i_pi_6(p: usize, cc: &mut Consts) -> (ExactComplex, ExactComplex) {
    let a = cc
        .pi(p, RoundingMode::None)
        .div(&ExactNum::from_u8(6, p), p, RoundingMode::None);
    let (s, c) = a.sin_cos(p, RoundingMode::None, cc);
    (
        ExactComplex::new(c.clone(), s.clone()),
        ExactComplex::new(c, s.neg()),
    )
}

fn two_pi_over_three(p: usize, cc: &mut Consts) -> ExactNum {
    two_c(p)
        .re()
        .mul(&cc.pi(p, RoundingMode::None), p, RoundingMode::None)
        .div(&ExactNum::from_u8(3, p), p, RoundingMode::None)
}

fn pi_over_three(p: usize, cc: &mut Consts) -> ExactNum {
    cc.pi(p, RoundingMode::None)
        .div(&ExactNum::from_u8(3, p), p, RoundingMode::None)
}

fn airy_cs_real(p: usize, cc: &mut Consts) -> (ExactNum, ExactNum, ExactNum) {
    let one = ExactNum::from_u8(1, p);
    let two = ExactNum::from_u8(2, p);
    let three = ExactNum::from_u8(3, p);
    let third = one.div(&three, p, RoundingMode::None);
    let two_third = two.div(&three, p, RoundingMode::None);
    let g23 = two_third.gamma(p, RoundingMode::None, cc);
    let g13 = third.gamma(p, RoundingMode::None, cc);
    let c1 =
        three
            .pow(&two_third.neg(), p, RoundingMode::None, cc)
            .div(&g23, p, RoundingMode::None);
    let c2 = three
        .pow(&third.neg(), p, RoundingMode::None, cc)
        .div(&g13, p, RoundingMode::None);
    let sqrt3 = three.sqrt(p, RoundingMode::None);
    (c1, c2, sqrt3)
}

impl ExactComplex {
    /// Airy \(\mathrm{Ai}(z)\). Entire. NaN in → NaN out.
    ///
    /// # Precision
    ///
    /// - Algorithm: series for `|z| < AIRY_SERIES_THRESHOLD` (`8`); asymptotic otherwise.
    /// - Bound: Ziv on each part (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    pub fn ai(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return ExactComplex::new(self.re().clone(), self.im().clone());
        }
        let dest = round_p(p);
        ziv_complex(dest, rm, |pw| self.ai_at(pw, dest, cc, true))
    }

    /// Airy \(\mathrm{Bi}(z)\). Entire. NaN in → NaN out.
    ///
    /// # Precision
    ///
    /// - Algorithm: same `AIRY_SERIES_THRESHOLD = 8` as [`Self::ai`].
    /// - Bound: Ziv on each part (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    pub fn bi(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return ExactComplex::new(self.re().clone(), self.im().clone());
        }
        let dest = round_p(p);
        ziv_complex(dest, rm, |pw| self.bi_at(pw, dest, cc))
    }

    fn ai_at(&self, work_p: usize, dest_p: usize, cc: &mut Consts, connect: bool) -> Self {
        if self.im().is_zero() {
            return ExactComplex::from_real(self.re().ai(work_p, RoundingMode::None, cc), work_p);
        }
        if abs_below(self, AIRY_SERIES_THRESHOLD, dest_p) {
            return self.airy_series_ai(work_p, cc);
        }
        let sector = two_pi_over_three(work_p, cc);
        if arg_abs_le(self, &sector, work_p, cc) {
            return self.airy_asymp_ai(work_p, cc);
        }
        if connect {
            let (w, w2) = omega_pair(work_p);
            let zw = self.mul(&w, work_p, RoundingMode::None);
            let zw2 = self.mul(&w2, work_p, RoundingMode::None);
            let t1 = w.mul(
                &zw.ai_at(work_p, dest_p, cc, false),
                work_p,
                RoundingMode::None,
            );
            let t2 = w2.mul(
                &zw2.ai_at(work_p, dest_p, cc, false),
                work_p,
                RoundingMode::None,
            );
            return neg_c(&t1.add(&t2, work_p, RoundingMode::None));
        }
        self.airy_asymp_ai(work_p, cc)
    }

    fn bi_at(&self, work_p: usize, dest_p: usize, cc: &mut Consts) -> Self {
        if self.im().is_zero() {
            return ExactComplex::from_real(self.re().bi(work_p, RoundingMode::None, cc), work_p);
        }
        if abs_below(self, AIRY_SERIES_THRESHOLD, dest_p) {
            return self.airy_series_bi(work_p, cc);
        }
        let sector = pi_over_three(work_p, cc);
        if arg_abs_le(self, &sector, work_p, cc) {
            return self.airy_asymp_bi(work_p, cc);
        }
        let (w, w2) = omega_pair(work_p);
        let (ep, em) = exp_i_pi_6(work_p, cc);
        let a1 = self
            .mul(&w2, work_p, RoundingMode::None)
            .ai_at(work_p, dest_p, cc, true);
        let a2 = self
            .mul(&w, work_p, RoundingMode::None)
            .ai_at(work_p, dest_p, cc, true);
        ep.mul(&a1, work_p, RoundingMode::None).add(
            &em.mul(&a2, work_p, RoundingMode::None),
            work_p,
            RoundingMode::None,
        )
    }

    fn airy_series_pair(&self, p: usize) -> (ExactComplex, ExactComplex) {
        let x3 = self
            .mul(self, p, RoundingMode::None)
            .mul(self, p, RoundingMode::None);
        let one = ExactComplex::one(p);
        let mut tf = one.clone();
        let mut f = one;
        let mut tg = self.clone();
        let mut g = self.clone();
        let cap = series_term_cap(p);
        for k in 1..=cap {
            let k3 = ExactComplex::from_real(ExactNum::from_u32((3 * k) as u32, p), p);
            let k3m1 = ExactComplex::from_real(ExactNum::from_u32((3 * k - 1) as u32, p), p);
            let k3p1 = ExactComplex::from_real(ExactNum::from_u32((3 * k + 1) as u32, p), p);
            tf = tf.mul(&x3, p, RoundingMode::None).div(
                &k3.mul(&k3m1, p, RoundingMode::None),
                p,
                RoundingMode::None,
            );
            f = f.add(&tf, p, RoundingMode::None);
            tg = tg.mul(&x3, p, RoundingMode::None).div(
                &k3p1.mul(&k3, p, RoundingMode::None),
                p,
                RoundingMode::None,
            );
            g = g.add(&tg, p, RoundingMode::None);
            if term_negligible(&tf, p) && term_negligible(&tg, p) {
                break;
            }
        }
        (f, g)
    }

    fn airy_series_ai(&self, p: usize, cc: &mut Consts) -> Self {
        let (c1, c2, _) = airy_cs_real(p, cc);
        let (f, g) = self.airy_series_pair(p);
        let c1c = ExactComplex::from_real(c1, p);
        let c2c = ExactComplex::from_real(c2, p);
        c1c.mul(&f, p, RoundingMode::None).sub(
            &c2c.mul(&g, p, RoundingMode::None),
            p,
            RoundingMode::None,
        )
    }

    fn airy_series_bi(&self, p: usize, cc: &mut Consts) -> Self {
        let (c1, c2, sqrt3) = airy_cs_real(p, cc);
        let (f, g) = self.airy_series_pair(p);
        let c1c = ExactComplex::from_real(c1, p);
        let c2c = ExactComplex::from_real(c2, p);
        let s3 = ExactComplex::from_real(sqrt3, p);
        s3.mul(
            &c1c.mul(&f, p, RoundingMode::None).add(
                &c2c.mul(&g, p, RoundingMode::None),
                p,
                RoundingMode::None,
            ),
            p,
            RoundingMode::None,
        )
    }

    fn airy_xi_z(&self, p: usize, cc: &mut Consts) -> (ExactComplex, ExactComplex) {
        let two_thirds = two_c(p).div(
            &ExactComplex::from_real(ExactNum::from_u8(3, p), p),
            p,
            RoundingMode::None,
        );
        let sz = self.sqrt(p, RoundingMode::None, cc);
        let xi = two_thirds.mul(&self.mul(&sz, p, RoundingMode::None), p, RoundingMode::None);
        let z14 = sz.sqrt(p, RoundingMode::None, cc);
        (xi, z14)
    }

    fn airy_u_sum(&self, xi: &ExactComplex, p: usize, alt: bool) -> ExactComplex {
        let one = ExactComplex::one(p);
        let mut u = one.clone();
        let mut sum = one;
        let mut prev_abs = ExactNum::from_u8(1, p);
        let mut xi_pow = xi.clone();
        let cap = series_term_cap(p).min(p.saturating_add(8));
        for k in 1..=cap {
            let num = ExactNum::from_u32((6 * k - 5) as u32, p).mul(
                &ExactNum::from_u32((6 * k - 1) as u32, p),
                p,
                RoundingMode::None,
            );
            let den = ExactNum::from_u32(72, p).mul(
                &ExactNum::from_u32(k as u32, p),
                p,
                RoundingMode::None,
            );
            let ratio = ExactComplex::from_real(num.div(&den, p, RoundingMode::None), p);
            u = u.mul(&ratio, p, RoundingMode::None);
            let mut t = u.div(&xi_pow, p, RoundingMode::None);
            if alt && k % 2 == 1 {
                t = neg_c(&t);
            }
            let ta = t.abs(p, RoundingMode::None);
            if matches!(ta.cmp(&prev_abs), Some(c) if c > 0) {
                break;
            }
            sum = sum.add(&t, p, RoundingMode::None);
            if term_negligible(&t, p) {
                break;
            }
            prev_abs = ta;
            xi_pow = xi_pow.mul(xi, p, RoundingMode::None);
        }
        sum
    }

    fn airy_asymp_ai(&self, p: usize, cc: &mut Consts) -> Self {
        let (xi, z14) = self.airy_xi_z(p, cc);
        let half = half_c(p);
        let pi = pi_c(p, cc);
        let pi_m12 =
            ExactComplex::one(p).div(&pi.sqrt(p, RoundingMode::None, cc), p, RoundingMode::None);
        let z_m14 = ExactComplex::one(p).div(&z14, p, RoundingMode::None);
        let exp_m = neg_c(&xi).exp(p, RoundingMode::None, cc);
        let su = self.airy_u_sum(&xi, p, true);
        half.mul(&pi_m12, p, RoundingMode::None)
            .mul(&z_m14, p, RoundingMode::None)
            .mul(&exp_m, p, RoundingMode::None)
            .mul(&su, p, RoundingMode::None)
    }

    fn airy_asymp_bi(&self, p: usize, cc: &mut Consts) -> Self {
        let (xi, z14) = self.airy_xi_z(p, cc);
        let pi = pi_c(p, cc);
        let pi_m12 =
            ExactComplex::one(p).div(&pi.sqrt(p, RoundingMode::None, cc), p, RoundingMode::None);
        let z_m14 = ExactComplex::one(p).div(&z14, p, RoundingMode::None);
        let exp_p = xi.exp(p, RoundingMode::None, cc);
        let su = self.airy_u_sum(&xi, p, false);
        pi_m12
            .mul(&z_m14, p, RoundingMode::None)
            .mul(&exp_p, p, RoundingMode::None)
            .mul(&su, p, RoundingMode::None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near_bits(a: &ExactNum, b: &ExactNum, p: usize, slack: i32) -> bool {
        let d = a.sub(b, p, RoundingMode::None).abs();
        d.is_zero() || d.exponent().is_some_and(|e| e < -((p as i32) / slack))
    }

    fn tiny(x: &ExactNum, p: usize) -> bool {
        x.is_zero() || x.exponent().is_some_and(|e| e < -((p as i32) / 4))
    }

    #[test]
    fn test_complex_airy_golds() {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new().unwrap();

        let z0 = ExactComplex::zero(p);
        let a0 = z0.ai(p, rm, &mut cc);
        let r0 = ExactNum::new(p).ai(p, rm, &mut cc);
        assert!(near_bits(a0.re(), &r0, p, 4));
        assert!(tiny(a0.im(), p));
        let b0 = z0.bi(p, rm, &mut cc);
        let rb0 = ExactNum::new(p).bi(p, rm, &mut cc);
        assert!(near_bits(b0.re(), &rb0, p, 4));
        assert!(tiny(b0.im(), p));

        let one = ExactComplex::one(p);
        let a1 = one.ai(p, rm, &mut cc);
        let r1 = ExactNum::from_u8(1, p).ai(p, rm, &mut cc);
        assert!(near_bits(a1.re(), &r1, p, 4));
        assert!(tiny(a1.im(), p));

        let i = ExactComplex::i(p);
        let zi = i.ai(p, rm, &mut cc);
        assert!(!zi.is_nan());
        let (w, w2) = omega_pair(p);
        // Ai(z) + ω Ai(ωz) + ω² Ai(ω²z) = 0 at z = i (series region).
        let t0 = zi;
        let t1 = w.mul(&i.mul(&w, p, rm).ai(p, rm, &mut cc), p, rm);
        let t2 = w2.mul(&i.mul(&w2, p, rm).ai(p, rm, &mut cc), p, rm);
        let sum = t0.add(&t1, p, rm).add(&t2, p, rm);
        assert!(
            tiny(sum.re(), p) && tiny(sum.im(), p),
            "connection identity"
        );

        let z = ExactComplex::new(ExactNum::from_u8(1, p), ExactNum::from_u8(1, p));
        let az = z.ai(p, rm, &mut cc);
        let bz = z.bi(p, rm, &mut cc);
        assert!(!az.is_nan() && !bz.is_nan());

        let nan = ExactComplex::new(crate::NAN.clone(), ExactNum::new(p));
        assert!(nan.ai(p, rm, &mut cc).is_nan());
        assert!(nan.bi(p, rm, &mut cc).is_nan());
    }
}
