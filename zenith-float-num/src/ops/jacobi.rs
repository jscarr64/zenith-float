//! Jacobi elliptic functions. Parameter \(m=k^2\). Real \(m\in[0,1]\).

use alloc::vec::Vec;

use crate::common::util::bump_prec_retry;
use crate::common::util::round_p;
use crate::defs::Error;
use crate::defs::Exponent;
use crate::defs::RoundingMode;
use crate::defs::Sign;
use crate::num::ExactNumNumber;
use crate::ops::consts::Consts;
use crate::WORD_BIT_SIZE;

/// AGM / descending Landen steps for `am`, `sn`, `cn`, `dn`.
pub const JACOBI_AGM_MAX: u32 = 128;

impl ExactNumNumber {
    /// Jacobi amplitude \(\operatorname{am}(u|m)\). Parameter \(m=k^2\in[0,1]\).
    pub fn jacobi_am(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        self.jacobi_ziv(m, p, rm, cc, |am, _sn, _cn, _dn, _px| am.clone())
    }

    /// \(\operatorname{sn}(u|m)=\sin\operatorname{am}(u|m)\).
    pub fn jacobi_sn(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        self.jacobi_ziv(m, p, rm, cc, |_am, sn, _cn, _dn, _px| sn.clone())
    }

    /// \(\operatorname{cn}(u|m)=\cos\operatorname{am}(u|m)\).
    pub fn jacobi_cn(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        self.jacobi_ziv(m, p, rm, cc, |_am, _sn, cn, _dn, _px| cn.clone())
    }

    /// \(\operatorname{dn}(u|m)=\sqrt{1-m\,\operatorname{sn}^2(u|m)}\).
    pub fn jacobi_dn(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        self.jacobi_ziv(m, p, rm, cc, |_am, _sn, _cn, dn, _px| dn.clone())
    }

    /// \(\operatorname{cd}(u|m)=\operatorname{cn}(u|m)/\operatorname{dn}(u|m)\).
    pub fn jacobi_cd(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        self.jacobi_ziv(m, p, rm, cc, |_am, _sn, cn, dn, px| div_nz(cn, dn, px))
    }

    /// \(\operatorname{ns}(u|m)=1/\operatorname{sn}(u|m)\).
    pub fn jacobi_ns(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        self.jacobi_ziv(m, p, rm, cc, |_am, sn, _cn, _dn, px| rec_nz(sn, px))
    }

    /// \(\operatorname{nc}(u|m)=1/\operatorname{cn}(u|m)\).
    pub fn jacobi_nc(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        self.jacobi_ziv(m, p, rm, cc, |_am, _sn, cn, _dn, px| rec_nz(cn, px))
    }

    /// \(\operatorname{nd}(u|m)=1/\operatorname{dn}(u|m)\).
    pub fn jacobi_nd(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        self.jacobi_ziv(m, p, rm, cc, |_am, _sn, _cn, dn, px| rec_nz(dn, px))
    }

    /// \(\operatorname{sc}(u|m)=\operatorname{sn}/\operatorname{cn}\).
    pub fn jacobi_sc(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        self.jacobi_ziv(m, p, rm, cc, |_am, sn, cn, _dn, px| div_nz(sn, cn, px))
    }

    /// \(\operatorname{sd}(u|m)=\operatorname{sn}/\operatorname{dn}\).
    pub fn jacobi_sd(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        self.jacobi_ziv(m, p, rm, cc, |_am, sn, _cn, dn, px| div_nz(sn, dn, px))
    }

    /// \(\operatorname{cs}(u|m)=\operatorname{cn}/\operatorname{sn}\).
    pub fn jacobi_cs(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        self.jacobi_ziv(m, p, rm, cc, |_am, sn, cn, _dn, px| div_nz(cn, sn, px))
    }

    /// \(\operatorname{ds}(u|m)=\operatorname{dn}/\operatorname{sn}\).
    pub fn jacobi_ds(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        self.jacobi_ziv(m, p, rm, cc, |_am, sn, _cn, dn, px| div_nz(dn, sn, px))
    }

    /// \(\operatorname{dc}(u|m)=\operatorname{dn}/\operatorname{cn}\).
    pub fn jacobi_dc(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        self.jacobi_ziv(m, p, rm, cc, |_am, _sn, cn, dn, px| div_nz(dn, cn, px))
    }

    fn jacobi_ziv<F>(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
        pick: F,
    ) -> Result<Self, Error>
    where
        F: Fn(&Self, &Self, &Self, &Self, usize) -> Result<Self, Error>,
    {
        let p = round_p(p);
        Self::p_assertion(p)?;
        jacobi_m_in_unit(m, p)?;
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p
            .max(self.mantissa_max_bit_len())
            .max(m.mantissa_max_bit_len())
            + p_inc;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE * 2;
            let (am, sn, cn, dn) = self.jacobi_sncndn_at(m, p_x, cc)?;
            let mut ret = pick(&am, &sn, &cn, &dn, p_x)?;
            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact() | m.inexact());
                return Ok(ret);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    fn jacobi_sncndn_at(
        &self,
        m: &Self,
        p: usize,
        cc: &mut Consts,
    ) -> Result<(Self, Self, Self, Self), Error> {
        let zero = Self::new2(p, Sign::Pos, false)?;
        let one = Self::from_word(1, p)?;
        if self.is_zero() {
            return Ok((zero.clone()?, zero, one.clone()?, one));
        }
        if m.is_zero() {
            let am = self.clone()?;
            let sn = self.sin(p, RoundingMode::None, cc)?;
            let cn = self.cos(p, RoundingMode::None, cc)?;
            return Ok((am, sn, cn, one));
        }
        if m.cmp(&one) == 0 {
            let sn = self.tanh(p, RoundingMode::None, cc)?;
            let ch = self.cosh(p, RoundingMode::None, cc)?;
            let cn = one.div(&ch, p, RoundingMode::None)?;
            let sh = self.sinh(p, RoundingMode::None, cc)?;
            let am = sh.atan(p, RoundingMode::None, cc)?;
            return Ok((am, sn, cn.clone()?, cn));
        }
        let u = self.jacobi_reduce_period(m, p, cc)?;
        jacobi_agm(&u, m, p, cc)
    }

    fn jacobi_reduce_period(&self, m: &Self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        // \(K(m)\ge\pi/2\), so \(2K\ge\pi\). Skip \(K\) (and a nested Ziv) when \(\lvert u\rvert<\pi\).
        let pi = cc.pi_num(p, RoundingMode::None)?;
        if self.abs()?.cmp(&pi) < 0 {
            return self.clone();
        }
        let k = m.elliptic_k_at(p, cc)?;
        let two = Self::from_word(2, p)?;
        let four = Self::from_word(4, p)?;
        let two_k = two.mul(&k, p, RoundingMode::None)?;
        if self.abs()?.cmp(&two_k) <= 0 {
            return self.clone();
        }
        let four_k = four.mul(&k, p, RoundingMode::None)?;
        let mut r = self.rem(&four_k)?;
        if r.cmp(&two_k) > 0 {
            r = r.sub(&four_k, p, RoundingMode::None)?;
        } else if r.is_negative() && r.abs()?.cmp(&two_k) > 0 {
            r = r.add(&four_k, p, RoundingMode::None)?;
        }
        Ok(r)
    }
}

fn div_nz(n: &ExactNumNumber, d: &ExactNumNumber, p: usize) -> Result<ExactNumNumber, Error> {
    if d.is_zero() {
        return Err(Error::InvalidArgument);
    }
    n.div(d, p, RoundingMode::None)
}

fn rec_nz(d: &ExactNumNumber, p: usize) -> Result<ExactNumNumber, Error> {
    if d.is_zero() {
        return Err(Error::InvalidArgument);
    }
    ExactNumNumber::from_word(1, p)?.div(d, p, RoundingMode::None)
}

fn jacobi_m_in_unit(m: &ExactNumNumber, p: usize) -> Result<(), Error> {
    if m.is_negative() {
        return Err(Error::InvalidArgument);
    }
    let one = ExactNumNumber::from_word(1, p)?;
    if m.cmp(&one) > 0 {
        return Err(Error::InvalidArgument);
    }
    Ok(())
}

fn agm_small(c: &ExactNumNumber, a: &ExactNumNumber, p: usize) -> Result<bool, Error> {
    if c.is_zero() {
        return Ok(true);
    }
    if a.is_zero() {
        return Ok(false);
    }
    let rel = c.abs()?.div(&a.abs()?, p, RoundingMode::None)?;
    let thresh = -((p as i32) / 2 + 16);
    Ok(rel.exponent() < thresh)
}

fn clamp_unit(x: &ExactNumNumber, p: usize) -> Result<ExactNumNumber, Error> {
    let one = ExactNumNumber::from_word(1, p)?;
    let ax = x.abs()?;
    if ax.cmp(&one) <= 0 {
        return x.clone();
    }
    let mut y = one;
    y.set_sign(x.sign());
    Ok(y)
}

fn jacobi_agm(
    u: &ExactNumNumber,
    m: &ExactNumNumber,
    p: usize,
    cc: &mut Consts,
) -> Result<
    (
        ExactNumNumber,
        ExactNumNumber,
        ExactNumNumber,
        ExactNumNumber,
    ),
    Error,
> {
    let one = ExactNumNumber::from_word(1, p)?;
    let two = ExactNumNumber::from_word(2, p)?;
    let om = one.sub(m, p, RoundingMode::None)?;
    if om.is_negative() {
        return Err(Error::InvalidArgument);
    }
    let mut b = om.sqrt(p, RoundingMode::None)?;
    let mut aa = Vec::new();
    let mut cc_seq = Vec::new();
    aa.push(one.clone()?);
    cc_seq.push(m.sqrt(p, RoundingMode::None)?);
    let mut n = 0u32;
    while n < JACOBI_AGM_MAX {
        if agm_small(&cc_seq[n as usize], &aa[n as usize], p)? {
            break;
        }
        let ai = aa[n as usize].clone()?;
        let c_next = ai
            .sub(&b, p, RoundingMode::None)?
            .div(&two, p, RoundingMode::None)?;
        let t = ai
            .mul(&b, p, RoundingMode::None)?
            .sqrt(p, RoundingMode::None)?;
        let a_next = ai
            .add(&b, p, RoundingMode::None)?
            .div(&two, p, RoundingMode::None)?;
        b = t;
        aa.push(a_next);
        cc_seq.push(c_next);
        n += 1;
    }
    if n == JACOBI_AGM_MAX && !agm_small(&cc_seq[n as usize], &aa[n as usize], p)? {
        return Err(Error::InvalidArgument);
    }
    let last = aa.len() - 1;
    let shift = Exponent::try_from(last).map_err(|_| Error::InvalidArgument)?;
    let twon = one.ldexp(shift, p, RoundingMode::None)?;
    let mut phi = twon
        .mul(&aa[last], p, RoundingMode::None)?
        .mul(u, p, RoundingMode::None)?;
    let mut idx = last;
    while idx > 0 {
        let sphi = phi.sin(p, RoundingMode::None, cc)?;
        let t = clamp_unit(
            &cc_seq[idx]
                .mul(&sphi, p, RoundingMode::None)?
                .div(&aa[idx], p, RoundingMode::None)?,
            p,
        )?;
        let asint = t.asin(p, RoundingMode::None, cc)?;
        phi = asint
            .add(&phi, p, RoundingMode::None)?
            .div(&two, p, RoundingMode::None)?;
        idx -= 1;
    }
    let sn = phi.sin(p, RoundingMode::None, cc)?;
    let cn = phi.cos(p, RoundingMode::None, cc)?;
    let msn2 = m.mul(&sn.mul(&sn, p, RoundingMode::None)?, p, RoundingMode::None)?;
    let dnarg = one.sub(&msn2, p, RoundingMode::None)?;
    if dnarg.is_negative() {
        return Err(Error::InvalidArgument);
    }
    let dn = dnarg.sqrt(p, RoundingMode::None)?;
    Ok((phi, sn, cn, dn))
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
    fn test_jacobi_identities() {
        let p = 256;
        let mut cc = Consts::new().unwrap();
        let rm = RoundingMode::ToEven;
        let one = ExactNumNumber::from_word(1, p).unwrap();
        let zero = ExactNumNumber::new(p).unwrap();
        let two = ExactNumNumber::from_word(2, p).unwrap();
        let mut half = one.clone().unwrap();
        half.div_by_2(RoundingMode::None);

        let sn0 = zero.jacobi_sn(&half, p, rm, &mut cc).unwrap();
        assert!(sn0.is_zero(), "sn(0)");
        let cn0 = zero.jacobi_cn(&half, p, rm, &mut cc).unwrap();
        bits_agree(&cn0, &one, p, 80, "cn(0)");
        let dn0 = zero.jacobi_dn(&half, p, rm, &mut cc).unwrap();
        bits_agree(&dn0, &one, p, 80, "dn(0)");
        let am0 = zero.jacobi_am(&half, p, rm, &mut cc).unwrap();
        assert!(am0.is_zero(), "am(0)");

        let u = one.clone().unwrap();
        let sn_m0 = u.jacobi_sn(&zero, p, rm, &mut cc).unwrap();
        let sinu = u.sin(p, rm, &mut cc).unwrap();
        bits_agree(&sn_m0, &sinu, p, 80, "sn(u|0)=sin u");
        let cn_m0 = u.jacobi_cn(&zero, p, rm, &mut cc).unwrap();
        let cosu = u.cos(p, rm, &mut cc).unwrap();
        bits_agree(&cn_m0, &cosu, p, 80, "cn(u|0)=cos u");
        let dn_m0 = u.jacobi_dn(&zero, p, rm, &mut cc).unwrap();
        bits_agree(&dn_m0, &one, p, 80, "dn(u|0)=1");

        let sn_m1 = u.jacobi_sn(&one, p, rm, &mut cc).unwrap();
        let tanhu = u.tanh(p, rm, &mut cc).unwrap();
        bits_agree(&sn_m1, &tanhu, p, 80, "sn(u|1)=tanh u");
        let cn_m1 = u.jacobi_cn(&one, p, rm, &mut cc).unwrap();
        let sech = one
            .div(&u.cosh(p, rm, &mut cc).unwrap(), p, RoundingMode::None)
            .unwrap();
        bits_agree(&cn_m1, &sech, p, 80, "cn(u|1)=sech u");
        let dn_m1 = u.jacobi_dn(&one, p, rm, &mut cc).unwrap();
        bits_agree(&dn_m1, &sech, p, 80, "dn(u|1)=sech u");

        let sn = u.jacobi_sn(&half, p, rm, &mut cc).unwrap();
        let cn = u.jacobi_cn(&half, p, rm, &mut cc).unwrap();
        let dn = u.jacobi_dn(&half, p, rm, &mut cc).unwrap();
        let s2c2 = sn
            .mul(&sn, p, RoundingMode::None)
            .unwrap()
            .add(
                &cn.mul(&cn, p, RoundingMode::None).unwrap(),
                p,
                RoundingMode::None,
            )
            .unwrap();
        bits_agree(&s2c2, &one, p, 80, "sn²+cn²=1");
        let d2ms = dn
            .mul(&dn, p, RoundingMode::None)
            .unwrap()
            .add(
                &half
                    .mul(
                        &sn.mul(&sn, p, RoundingMode::None).unwrap(),
                        p,
                        RoundingMode::None,
                    )
                    .unwrap(),
                p,
                RoundingMode::None,
            )
            .unwrap();
        bits_agree(&d2ms, &one, p, 80, "dn²+m sn²=1");

        let k = half.elliptic_k(p, rm, &mut cc).unwrap();
        let mut kh = k.clone().unwrap();
        kh.div_by_2(RoundingMode::None);
        let snkh = kh.jacobi_sn(&half, p, rm, &mut cc).unwrap();
        let kp = one
            .sub(&half, p, RoundingMode::None)
            .unwrap()
            .sqrt(p, RoundingMode::None)
            .unwrap();
        let want = one
            .div(
                &one.add(&kp, p, RoundingMode::None)
                    .unwrap()
                    .sqrt(p, RoundingMode::None)
                    .unwrap(),
                p,
                RoundingMode::None,
            )
            .unwrap();
        bits_agree(&snkh, &want, p, 40, "sn(K/2|1/2)");

        let four_k = four_word(p).mul(&k, p, RoundingMode::None).unwrap();
        let u4 = u.add(&four_k, p, RoundingMode::None).unwrap();
        let sn_per = u4.jacobi_sn(&half, p, rm, &mut cc).unwrap();
        bits_agree(&sn_per, &sn, p, 40, "sn(u+4K)=sn(u)");

        let mut h = one.clone().unwrap();
        for _ in 0..20 {
            h.div_by_2(RoundingMode::None);
        }
        let up = u.add(&h, p, RoundingMode::None).unwrap();
        let um = u.sub(&h, p, RoundingMode::None).unwrap();
        let snp = up.jacobi_sn(&half, p, rm, &mut cc).unwrap();
        let snm = um.jacobi_sn(&half, p, rm, &mut cc).unwrap();
        let fd = snp
            .sub(&snm, p, RoundingMode::None)
            .unwrap()
            .div(
                &two.mul(&h, p, RoundingMode::None).unwrap(),
                p,
                RoundingMode::None,
            )
            .unwrap();
        let deriv = cn.mul(&dn, p, RoundingMode::None).unwrap();
        bits_agree(&fd, &deriv, p, 20, "d(sn)/du = cn dn");

        let sn_none = u
            .jacobi_sn(&half, 512, RoundingMode::None, &mut cc)
            .unwrap();
        bits_agree(&sn_none, &sn, p, 20, "sn None 512 vs ToEven 256");

        assert!(u.jacobi_sn(&one.neg().unwrap(), p, rm, &mut cc).is_err());
        assert!(u.jacobi_sn(&two, p, rm, &mut cc).is_err());
        let cd = u.jacobi_cd(&half, p, rm, &mut cc).unwrap();
        let want_cd = cn.div(&dn, p, RoundingMode::None).unwrap();
        bits_agree(&cd, &want_cd, p, 80, "cd=cn/dn");
        let ns = u.jacobi_ns(&half, p, rm, &mut cc).unwrap();
        let nssn = ns.mul(&sn, p, RoundingMode::None).unwrap();
        bits_agree(&nssn, &one, p, 80, "ns·sn=1");
        assert!(zero.jacobi_ns(&half, p, rm, &mut cc).is_err());

        // Inverse of sn is incomplete F: F(sn(u|m)|m)=u for |u|<K(m).
        let back = sn.elliptic_f(&half, p, rm, &mut cc).unwrap();
        bits_agree(&back, &u, p, 40, "F(sn(u|m)|m)=u");
    }

    fn four_word(p: usize) -> ExactNumNumber {
        ExactNumNumber::from_word(4, p).unwrap()
    }
}
