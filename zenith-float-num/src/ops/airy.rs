//! Airy \(\mathrm{Ai}\) / \(\mathrm{Bi}\) and their first derivatives.

use crate::common::util::bump_prec_retry;
use crate::common::util::round_p;
use crate::defs::Error;
use crate::defs::RoundingMode;
use crate::num::ExactNumNumber;
use crate::ops::consts::Consts;
use crate::Sign;
use crate::Word;
use crate::WORD_BIT_SIZE;

/// \(|x|\) below this uses the Taylor pair \((f,g)\); at or above, the asymptotic.
const AIRY_SERIES_THRESHOLD: Word = 8;

impl ExactNumNumber {
    /// Airy \(\mathrm{Ai}(\mathrm{self})\) at precision `p`.
    pub fn ai(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        self.airy_component(AiryWhich::Ai, p, rm, cc)
    }

    /// Airy \(\mathrm{Bi}(\mathrm{self})\) at precision `p`.
    pub fn bi(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        self.airy_component(AiryWhich::Bi, p, rm, cc)
    }

    /// \(\mathrm{Ai}'(\mathrm{self})\) at precision `p`.
    pub fn ai_prime(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        self.airy_component(AiryWhich::AiPrime, p, rm, cc)
    }

    /// \(\mathrm{Bi}'(\mathrm{self})\) at precision `p`.
    pub fn bi_prime(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        self.airy_component(AiryWhich::BiPrime, p, rm, cc)
    }

    fn airy_component(
        &self,
        which: AiryWhich,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        Self::p_assertion(p)?;

        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc;

        loop {
            let p_x = p_wrk + WORD_BIT_SIZE;
            let mut ret = self.airy_at(which, p_x, cc)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn airy_at(&self, which: AiryWhich, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        let mut x = self.clone()?;
        x.set_inexact(false);
        let (ai, bi, aip, bip) = x.airy_all(p, cc)?;
        Ok(match which {
            AiryWhich::Ai => ai,
            AiryWhich::Bi => bi,
            AiryWhich::AiPrime => aip,
            AiryWhich::BiPrime => bip,
        })
    }

    fn airy_all(&self, p: usize, cc: &mut Consts) -> Result<(Self, Self, Self, Self), Error> {
        let (c1, c2, sqrt3) = airy_cs(p, cc)?;
        if self.is_zero() {
            let bi0 = sqrt3.mul(&c1, p, RoundingMode::None)?;
            let aip0 = c2.clone()?.neg()?;
            let bip0 = sqrt3.mul(&c2, p, RoundingMode::None)?;
            return Ok((c1, bi0, aip0, bip0));
        }
        let thr = Self::from_word(AIRY_SERIES_THRESHOLD, p)?;
        if self.abs_cmp(&thr) < 0 {
            let (f, g, fp, gp) = self.airy_fg_series(p)?;
            return combine_fg(&c1, &c2, &sqrt3, &f, &g, &fp, &gp, p);
        }
        if self.is_positive() {
            self.airy_asymp_pos(p, cc)
        } else {
            self.airy_asymp_neg(p, cc)
        }
    }

    fn airy_fg_series(&self, p: usize) -> Result<(Self, Self, Self, Self), Error> {
        let x3 = self
            .mul(self, p, RoundingMode::None)?
            .mul(self, p, RoundingMode::None)?;
        let one = Self::from_word(1, p)?;
        let mut tf = one.clone()?;
        let mut f = one.clone()?;
        let mut fp = Self::new2(p, Sign::Pos, false)?;
        let mut tg = self.clone()?;
        let mut g = self.clone()?;
        let mut gp = one;
        let cap = p.saturating_add(WORD_BIT_SIZE);
        for k in 1..=cap {
            let k3 = Self::from_word((3 * k) as Word, p)?;
            let k3m1 = Self::from_word((3 * k - 1) as Word, p)?;
            let k3p1 = Self::from_word((3 * k + 1) as Word, p)?;
            tf = tf
                .mul(&x3, p, RoundingMode::None)?
                .div(&k3.mul(&k3m1, p, RoundingMode::None)?, p, RoundingMode::None)?;
            f = f.add(&tf, p, RoundingMode::None)?;
            let dtf = tf.mul(&k3, p, RoundingMode::None)?.div(self, p, RoundingMode::None)?;
            fp = fp.add(&dtf, p, RoundingMode::None)?;
            tg = tg
                .mul(&x3, p, RoundingMode::None)?
                .div(&k3p1.mul(&k3, p, RoundingMode::None)?, p, RoundingMode::None)?;
            g = g.add(&tg, p, RoundingMode::None)?;
            let dtg = tg
                .mul(&k3p1, p, RoundingMode::None)?
                .div(self, p, RoundingMode::None)?;
            gp = gp.add(&dtg, p, RoundingMode::None)?;
            let tf_done = tf.is_zero() || (tf.exponent() as isize) + (p as isize) < 0;
            let tg_done = tg.is_zero() || (tg.exponent() as isize) + (p as isize) < 0;
            if tf_done && tg_done {
                break;
            }
        }
        Ok((f, g, fp, gp))
    }

    fn airy_asymp_pos(&self, p: usize, cc: &mut Consts) -> Result<(Self, Self, Self, Self), Error> {
        let (xi, x14, x_m14, pi_m12) = self.airy_xi_scales(p, cc)?;
        let mut nxi = xi.clone()?;
        nxi.set_sign(Sign::Neg);
        let exp_m = nxi.exp(p, RoundingMode::None, cc)?;
        let exp_p = xi.exp(p, RoundingMode::None, cc)?;
        let (su, sv) = airy_uv_sums(&xi, p, true)?;
        let (tu, tv) = airy_uv_sums(&xi, p, false)?;
        let half = one_half_airy(p)?;
        let pre_ai = half
            .mul(&pi_m12, p, RoundingMode::None)?
            .mul(&x_m14, p, RoundingMode::None)?
            .mul(&exp_m, p, RoundingMode::None)?;
        let pre_bi = pi_m12
            .mul(&x_m14, p, RoundingMode::None)?
            .mul(&exp_p, p, RoundingMode::None)?;
        let pre_aip = half
            .mul(&pi_m12, p, RoundingMode::None)?
            .mul(&x14, p, RoundingMode::None)?
            .mul(&exp_m, p, RoundingMode::None)?
            .neg()?;
        let pre_bip = pi_m12
            .mul(&x14, p, RoundingMode::None)?
            .mul(&exp_p, p, RoundingMode::None)?;
        Ok((
            pre_ai.mul(&su, p, RoundingMode::None)?,
            pre_bi.mul(&tu, p, RoundingMode::None)?,
            pre_aip.mul(&sv, p, RoundingMode::None)?,
            pre_bip.mul(&tv, p, RoundingMode::None)?,
        ))
    }

    fn airy_asymp_neg(&self, p: usize, cc: &mut Consts) -> Result<(Self, Self, Self, Self), Error> {
        let mut z = self.clone()?;
        z.set_sign(Sign::Pos);
        let (zeta, z14, z_m14, pi_m12) = z.airy_xi_scales(p, cc)?;
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let four = Self::from_word(4, p)?;
        let pi4 = pi.div(&four, p, RoundingMode::None)?;
        let chi = zeta.add(&pi4, p, RoundingMode::None)?;
        let (s, c) = chi.sin_cos(p, RoundingMode::None, cc)?;
        let (p_even, q_odd) = airy_pq_sums(&zeta, p, false)?;
        let (r_even, s_odd) = airy_pq_sums(&zeta, p, true)?;
        let amp = pi_m12.mul(&z_m14, p, RoundingMode::None)?;
        let amp_p = pi_m12.mul(&z14, p, RoundingMode::None)?;
        // Ai(-z) ~ π^{-1/2} z^{-1/4} (sin χ P − cos χ Q)
        let ai = amp.mul(
            &s.mul(&p_even, p, RoundingMode::None)?
                .sub(&c.mul(&q_odd, p, RoundingMode::None)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
        // Bi(-z) ~ π^{-1/2} z^{-1/4} (cos χ P + sin χ Q)
        let bi = amp.mul(
            &c.mul(&p_even, p, RoundingMode::None)?
                .add(&s.mul(&q_odd, p, RoundingMode::None)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
        // Ai'(-z) ~ −π^{-1/2} z^{1/4} (cos χ R + sin χ S)
        let aip = amp_p
            .mul(
                &c.mul(&r_even, p, RoundingMode::None)?
                    .add(&s.mul(&s_odd, p, RoundingMode::None)?, p, RoundingMode::None)?,
                p,
                RoundingMode::None,
            )?
            .neg()?;
        // Bi'(-z) ~ π^{-1/2} z^{1/4} (sin χ R − cos χ S)
        let bip = amp_p.mul(
            &s.mul(&r_even, p, RoundingMode::None)?
                .sub(&c.mul(&s_odd, p, RoundingMode::None)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
        Ok((ai, bi, aip, bip))
    }

    /// \(\xi=(2/3)x^{3/2}\), \(x^{1/4}\), \(x^{-1/4}\), \(\pi^{-1/2}\) for \(x>0\).
    fn airy_xi_scales(
        &self,
        p: usize,
        cc: &mut Consts,
    ) -> Result<(Self, Self, Self, Self), Error> {
        let two = Self::from_word(2, p)?;
        let three = Self::from_word(3, p)?;
        let sx = self.sqrt(p, RoundingMode::None)?;
        let x32 = self.mul(&sx, p, RoundingMode::None)?;
        let xi = two
            .div(&three, p, RoundingMode::None)?
            .mul(&x32, p, RoundingMode::None)?;
        let x14 = sx.sqrt(p, RoundingMode::None)?;
        let one = Self::from_word(1, p)?;
        let x_m14 = one.div(&x14, p, RoundingMode::None)?;
        let pi = cc.pi_num(p, RoundingMode::None)?;
        let pi_m12 = one.div(&pi.sqrt(p, RoundingMode::None)?, p, RoundingMode::None)?;
        Ok((xi, x14, x_m14, pi_m12))
    }
}

#[derive(Clone, Copy)]
enum AiryWhich {
    Ai,
    Bi,
    AiPrime,
    BiPrime,
}

fn one_half_airy(p: usize) -> Result<ExactNumNumber, Error> {
    let mut h = ExactNumNumber::from_word(1, p)?;
    h.div_by_2(RoundingMode::None);
    Ok(h)
}

fn airy_cs(
    p: usize,
    cc: &mut Consts,
) -> Result<(ExactNumNumber, ExactNumNumber, ExactNumNumber), Error> {
    let one = ExactNumNumber::from_word(1, p)?;
    let two = ExactNumNumber::from_word(2, p)?;
    let three = ExactNumNumber::from_word(3, p)?;
    let third = one.div(&three, p, RoundingMode::None)?;
    let two_third = two.div(&three, p, RoundingMode::None)?;
    let g23 = two_third.gamma(p, RoundingMode::None, cc)?;
    let g13 = third.gamma(p, RoundingMode::None, cc)?;
    let mut n2t = two_third.clone()?;
    n2t.set_sign(Sign::Neg);
    let mut nt = third.clone()?;
    nt.set_sign(Sign::Neg);
    let c1 = three
        .pow(&n2t, p, RoundingMode::None, cc)?
        .div(&g23, p, RoundingMode::None)?;
    let c2 = three
        .pow(&nt, p, RoundingMode::None, cc)?
        .div(&g13, p, RoundingMode::None)?;
    let sqrt3 = three.sqrt(p, RoundingMode::None)?;
    Ok((c1, c2, sqrt3))
}

fn combine_fg(
    c1: &ExactNumNumber,
    c2: &ExactNumNumber,
    sqrt3: &ExactNumNumber,
    f: &ExactNumNumber,
    g: &ExactNumNumber,
    fp: &ExactNumNumber,
    gp: &ExactNumNumber,
    p: usize,
) -> Result<(ExactNumNumber, ExactNumNumber, ExactNumNumber, ExactNumNumber), Error> {
    let c1f = c1.mul(f, p, RoundingMode::None)?;
    let c2g = c2.mul(g, p, RoundingMode::None)?;
    let c1fp = c1.mul(fp, p, RoundingMode::None)?;
    let c2gp = c2.mul(gp, p, RoundingMode::None)?;
    let ai = c1f.sub(&c2g, p, RoundingMode::None)?;
    let bi = sqrt3.mul(&c1f.add(&c2g, p, RoundingMode::None)?, p, RoundingMode::None)?;
    let aip = c1fp.sub(&c2gp, p, RoundingMode::None)?;
    let bip = sqrt3.mul(&c1fp.add(&c2gp, p, RoundingMode::None)?, p, RoundingMode::None)?;
    Ok((ai, bi, aip, bip))
}

/// \(\sum u_k\xi^{-k}\) and \(\sum v_k\xi^{-k}\). `alt` applies \((-1)^k\).
fn airy_uv_sums(
    xi: &ExactNumNumber,
    p: usize,
    alt: bool,
) -> Result<(ExactNumNumber, ExactNumNumber), Error> {
    let one = ExactNumNumber::from_word(1, p)?;
    let mut u = one.clone()?;
    let mut sum_u = one.clone()?;
    let mut sum_v = one;
    let mut prev_u = sum_u.clone()?;
    let mut xi_pow = xi.clone()?;
    let cap = p.saturating_add(8);
    for k in 1..=cap {
        let num = ExactNumNumber::from_word((6 * k - 5) as Word, p)?
            .mul(
                &ExactNumNumber::from_word((6 * k - 1) as Word, p)?,
                p,
                RoundingMode::None,
            )?;
        let den = ExactNumNumber::from_word(72, p)?.mul(
            &ExactNumNumber::from_word(k as Word, p)?,
            p,
            RoundingMode::None,
        )?;
        u = u.mul(&num, p, RoundingMode::None)?.div(&den, p, RoundingMode::None)?;
        let mut tu = u.div(&xi_pow, p, RoundingMode::None)?;
        let vfac = ExactNumNumber::from_word((6 * k + 1) as Word, p)?.div(
            &ExactNumNumber::from_word((6 * k - 1) as Word, p)?,
            p,
            RoundingMode::None,
        )?;
        let mut tv = tu.mul(&vfac, p, RoundingMode::None)?;
        if alt && k % 2 == 1 {
            tu.set_sign(Sign::Neg);
            tv.set_sign(Sign::Neg);
        }
        if tu.abs_cmp(&prev_u) > 0 {
            break;
        }
        sum_u = sum_u.add(&tu, p, RoundingMode::None)?;
        sum_v = sum_v.add(&tv, p, RoundingMode::None)?;
        if tu.is_zero() || (tu.exponent() as isize) + (p as isize) < 0 {
            break;
        }
        prev_u = tu.abs()?;
        xi_pow = xi_pow.mul(xi, p, RoundingMode::None)?;
    }
    Ok((sum_u, sum_v))
}

/// Even / odd grouped sums for the oscillatory expansion.
/// `use_v`: \(v_k\) instead of \(u_k\).
fn airy_pq_sums(
    zeta: &ExactNumNumber,
    p: usize,
    use_v: bool,
) -> Result<(ExactNumNumber, ExactNumNumber), Error> {
    let one = ExactNumNumber::from_word(1, p)?;
    let mut even = one.clone()?;
    let mut odd = ExactNumNumber::new2(p, Sign::Pos, false)?;
    let mut u = one;
    let mut prev = even.clone()?;
    let mut zpow = ExactNumNumber::from_word(1, p)?;
    let cap = p.saturating_add(8);
    for k in 1..=cap {
        let num = ExactNumNumber::from_word((6 * k - 5) as Word, p)?
            .mul(
                &ExactNumNumber::from_word((6 * k - 1) as Word, p)?,
                p,
                RoundingMode::None,
            )?;
        let den = ExactNumNumber::from_word(72, p)?.mul(
            &ExactNumNumber::from_word(k as Word, p)?,
            p,
            RoundingMode::None,
        )?;
        u = u.mul(&num, p, RoundingMode::None)?.div(&den, p, RoundingMode::None)?;
        zpow = zpow.mul(zeta, p, RoundingMode::None)?;
        let mut t = u.div(&zpow, p, RoundingMode::None)?;
        if use_v {
            let vfac = ExactNumNumber::from_word((6 * k + 1) as Word, p)?.div(
                &ExactNumNumber::from_word((6 * k - 1) as Word, p)?,
                p,
                RoundingMode::None,
            )?;
            t = t.mul(&vfac, p, RoundingMode::None)?;
        }
        // (−1)^{⌊k/2⌋} on u_k ζ^{-k} for both P (even k) and Q (odd k).
        let sign_neg = (k / 2) % 2 == 1;
        if sign_neg {
            t.set_sign(Sign::Neg);
        }
        if t.abs_cmp(&prev) > 0 {
            break;
        }
        if k % 2 == 0 {
            even = even.add(&t, p, RoundingMode::None)?;
        } else {
            odd = odd.add(&t, p, RoundingMode::None)?;
        }
        if t.is_zero() || (t.exponent() as isize) + (p as isize) < 0 {
            break;
        }
        prev = t.abs()?;
    }
    Ok((even, odd))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::defs::RoundingMode;

    fn bits_agree(a: &ExactNumNumber, b: &ExactNumNumber, p: usize, min_bits: i32, label: &str) {
        let d = a.sub(b, p, RoundingMode::None).unwrap().abs().unwrap();
        if d.is_zero() {
            return;
        }
        let rel = d.exponent() - a.exponent();
        assert!(
            rel < -min_bits,
            "{label}: relative exponent {rel} (a.exp={}, b.exp={})",
            a.exponent(),
            b.exponent()
        );
    }

    #[test]
    fn test_airy_golds() {
        let p = 256;
        let mut cc = Consts::new().unwrap();
        let rm = RoundingMode::ToEven;
        let zero = ExactNumNumber::new(p).unwrap();
        let one = ExactNumNumber::from_word(1, p).unwrap();
        let (c1, _c2, sqrt3) = airy_cs(p, &mut cc).unwrap();
        let ai0 = zero.ai(p, rm, &mut cc).unwrap();
        bits_agree(&ai0, &c1, p, 80, "Ai(0)=c1");
        let bi0 = zero.bi(p, rm, &mut cc).unwrap();
        let three = ExactNumNumber::from_word(3, p).unwrap();
        let six = ExactNumNumber::from_word(6, p).unwrap();
        let one_n = ExactNumNumber::from_word(1, p).unwrap();
        let mut m16 = one_n.div(&six, p, RoundingMode::None).unwrap();
        m16.set_sign(Sign::Neg);
        let two_third = ExactNumNumber::from_word(2, p)
            .unwrap()
            .div(&three, p, RoundingMode::None)
            .unwrap();
        let g23 = two_third.gamma(p, RoundingMode::None, &mut cc).unwrap();
        let bi0_want = three
            .pow(&m16, p, RoundingMode::None, &mut cc)
            .unwrap()
            .div(&g23, p, RoundingMode::None)
            .unwrap();
        bits_agree(&bi0, &bi0_want, p, 80, "Bi(0)=3^{-1/6}/Γ(2/3)");
        bits_agree(
            &bi0,
            &sqrt3.mul(&c1, p, RoundingMode::None).unwrap(),
            p,
            80,
            "Bi(0)=√3 c1",
        );

        let x = one.clone().unwrap();
        let ai = x.ai(p, rm, &mut cc).unwrap();
        let bi = x.bi(p, rm, &mut cc).unwrap();
        let aip = x.ai_prime(p, rm, &mut cc).unwrap();
        let bip = x.bi_prime(p, rm, &mut cc).unwrap();
        let wr = ai
            .mul(&bip, p, RoundingMode::None)
            .unwrap()
            .sub(&aip.mul(&bi, p, RoundingMode::None).unwrap(), p, RoundingMode::None)
            .unwrap();
        let pi = cc.pi_num(p, RoundingMode::None).unwrap();
        let wr_want = ExactNumNumber::from_word(1, p)
            .unwrap()
            .div(&pi, p, RoundingMode::None)
            .unwrap();
        bits_agree(&wr, &wr_want, p, 40, "Wronskian at x=1");

        let xm = one.neg().unwrap();
        let wrn = {
            let ai = xm.ai(p, rm, &mut cc).unwrap();
            let bi = xm.bi(p, rm, &mut cc).unwrap();
            let aip = xm.ai_prime(p, rm, &mut cc).unwrap();
            let bip = xm.bi_prime(p, rm, &mut cc).unwrap();
            ai.mul(&bip, p, RoundingMode::None)
                .unwrap()
                .sub(&aip.mul(&bi, p, RoundingMode::None).unwrap(), p, RoundingMode::None)
                .unwrap()
        };
        bits_agree(&wrn, &wr_want, p, 40, "Wronskian at x=-1");

        let ten = ExactNumNumber::from_word(10, p).unwrap();
        assert!(ten.ai(p, rm, &mut cc).unwrap().is_positive());
        assert!(ten.bi(p, rm, &mut cc).unwrap().is_positive());

        // Ai'' ≈ x Ai via centered difference on Ai'.
        let mut h = ExactNumNumber::from_word(1, p).unwrap();
        h.set_exponent(-((p as i32) / 8));
        let xp = x.add(&h, p, RoundingMode::None).unwrap();
        let xm = x.sub(&h, p, RoundingMode::None).unwrap();
        let d2 = xp
            .ai_prime(p, rm, &mut cc)
            .unwrap()
            .sub(&xm.ai_prime(p, rm, &mut cc).unwrap(), p, RoundingMode::None)
            .unwrap()
            .div(
                &h.mul(&ExactNumNumber::from_word(2, p).unwrap(), p, RoundingMode::None)
                    .unwrap(),
                p,
                RoundingMode::None,
            )
            .unwrap();
        let xai = x.mul(&ai, p, RoundingMode::None).unwrap();
        bits_agree(&d2, &xai, p, 8, "Ai''=x Ai at x=1");
    }
}
