//! Jacobi elliptic functions. Parameter \(m=k^2\). Real \(m\in[0,1]\).

use alloc::vec::Vec;

use crate::common::util::bump_prec_retry;
use crate::common::util::prec_retry_exhausted;
use crate::common::util::reduction_width;
use crate::common::util::round_p;
use crate::defs::Error;
use crate::defs::Exponent;
use crate::defs::RoundingMode;
use crate::defs::Sign;
use crate::defs::EXPONENT_MIN;
use crate::num::ExactNumNumber;
use crate::ops::consts::Consts;
use crate::WORD_BIT_SIZE;

/// At `m = 1`, `cn = dn = sech u`, so `cd = dc = 1` exactly. Returning it directly avoids a
/// Ziv loop that can never certify an exactly representable result.
fn jacobi_ratio_m_one(m: &ExactNumNumber, p: usize) -> Result<Option<ExactNumNumber>, Error> {
    let p = round_p(p);
    ExactNumNumber::p_assertion(p)?;
    jacobi_m_in_unit(m, p)?;
    let one = ExactNumNumber::from_word(1, p)?;
    if m.cmp(&one) == 0 {
        let mut r = one;
        r.set_inexact(m.inexact());
        Ok(Some(r))
    } else {
        Ok(None)
    }
}

/// Historical public constant; AGM now runs until `|c|` is below the
/// working-precision bound, or stalls.
#[allow(dead_code)]
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
        self.jacobi_ziv(m, p, rm, cc, |am, _sn, _cn, _dn, _px| Ok(am))
    }

    /// \(\operatorname{sn}(u|m)=\sin\operatorname{am}(u|m)\).
    pub fn jacobi_sn(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        self.jacobi_ziv(m, p, rm, cc, |_am, sn, _cn, _dn, _px| Ok(sn))
    }

    /// \(\operatorname{cn}(u|m)=\cos\operatorname{am}(u|m)\).
    pub fn jacobi_cn(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        self.jacobi_ziv(m, p, rm, cc, |_am, _sn, cn, _dn, _px| Ok(cn))
    }

    /// \(\operatorname{dn}(u|m)=\sqrt{1-m\,\operatorname{sn}^2(u|m)}\).
    pub fn jacobi_dn(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        self.jacobi_ziv(m, p, rm, cc, |_am, _sn, _cn, dn, _px| Ok(dn))
    }

    /// \(\operatorname{cd}(u|m)=\operatorname{cn}(u|m)/\operatorname{dn}(u|m)\).
    pub fn jacobi_cd(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        if let Some(one) = jacobi_ratio_m_one(m, p)? {
            return Ok(one);
        }
        self.jacobi_ziv(m, p, rm, cc, |_am, _sn, cn, dn, px| div_nz(&cn, &dn, px))
    }

    /// \(\operatorname{ns}(u|m)=1/\operatorname{sn}(u|m)\).
    pub fn jacobi_ns(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        self.jacobi_ziv(m, p, rm, cc, |_am, sn, _cn, _dn, px| rec_nz(&sn, px))
    }

    /// \(\operatorname{nc}(u|m)=1/\operatorname{cn}(u|m)\).
    pub fn jacobi_nc(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        self.jacobi_ziv(m, p, rm, cc, |_am, _sn, cn, _dn, px| rec_nz(&cn, px))
    }

    /// \(\operatorname{nd}(u|m)=1/\operatorname{dn}(u|m)\).
    pub fn jacobi_nd(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        self.jacobi_ziv(m, p, rm, cc, |_am, _sn, _cn, dn, px| rec_nz(&dn, px))
    }

    /// \(\operatorname{sc}(u|m)=\operatorname{sn}/\operatorname{cn}\).
    pub fn jacobi_sc(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        self.jacobi_ziv(m, p, rm, cc, |_am, sn, cn, _dn, px| div_nz(&sn, &cn, px))
    }

    /// \(\operatorname{sd}(u|m)=\operatorname{sn}/\operatorname{dn}\).
    pub fn jacobi_sd(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        self.jacobi_ziv(m, p, rm, cc, |_am, sn, _cn, dn, px| div_nz(&sn, &dn, px))
    }

    /// \(\operatorname{cs}(u|m)=\operatorname{cn}/\operatorname{sn}\).
    pub fn jacobi_cs(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        self.jacobi_ziv(m, p, rm, cc, |_am, sn, cn, _dn, px| div_nz(&cn, &sn, px))
    }

    /// \(\operatorname{ds}(u|m)=\operatorname{dn}/\operatorname{sn}\).
    pub fn jacobi_ds(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        self.jacobi_ziv(m, p, rm, cc, |_am, sn, _cn, dn, px| div_nz(&dn, &sn, px))
    }

    /// \(\operatorname{dc}(u|m)=\operatorname{dn}/\operatorname{cn}\).
    pub fn jacobi_dc(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        if let Some(one) = jacobi_ratio_m_one(m, p)? {
            return Ok(one);
        }
        self.jacobi_ziv(m, p, rm, cc, |_am, _sn, cn, dn, px| div_nz(&dn, &cn, px))
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
        F: Fn(Self, Self, Self, Self, usize) -> Result<Self, Error>,
    {
        let p = round_p(p);
        Self::p_assertion(p)?;
        jacobi_m_in_unit(m, p)?;
        if let Some(r) = self.jacobi_m1_far(m, p, rm, cc, &pick)? {
            return Ok(r);
        }
        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p
            .max(self.mantissa_max_bit_len())
            .max(m.mantissa_max_bit_len())
            + p_inc;
        let mut first_lost: Option<usize> = None;
        loop {
            let p_x = p_wrk + WORD_BIT_SIZE * 2;
            let (am, sn, cn, dn, good) = self.jacobi_sncndn_at(m, p_x, cc)?;
            // `good` bounds the relative error of am, sn, cn and dn (2^-good each); a quotient of
            // two of them, rounded once more, has at least good - 2 bits.
            // (Whole words: the roundability check works on word boundaries.)
            let acc = good.saturating_sub(4).min(p_wrk) / WORD_BIT_SIZE * WORD_BIT_SIZE;
            if acc > p {
                let mut ret = pick(am, sn, cn, dn, p_x)?;
                if ret.try_set_precision(p, rm, acc)? {
                    ret.set_inexact(ret.inexact() | self.inexact() | m.inexact());
                    return Ok(ret);
                }
            }
            // Bits lost to cancellation (cn or dn near a zero, asin near ±1 in the Landen
            // recurrence) are restored at once. The retry budget is measured past the loss seen at
            // the first evaluation (fixed), so a bound that does not improve still ends.
            let lost = p_wrk.saturating_sub(acc);
            let allow = *first_lost.get_or_insert(lost);
            if good > 0 && lost > WORD_BIT_SIZE {
                if prec_retry_exhausted(p_wrk + lost, p.saturating_add(allow)) {
                    return Err(Error::PrecisionRetryExhausted);
                }
                p_wrk = round_p(p_wrk + lost);
            } else if prec_retry_exhausted(p_wrk, p.saturating_add(allow)) {
                return Err(Error::PrecisionRetryExhausted);
            } else {
                bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
            }
        }
    }

    /// am, sn, cn, dn at `p` working bits and the number of bits `g` certified for all four
    /// (each has relative error below `2^-g`).
    fn jacobi_sncndn_at(
        &self,
        m: &Self,
        p: usize,
        cc: &mut Consts,
    ) -> Result<(Self, Self, Self, Self, usize), Error> {
        // Correctly rounded building blocks (sin, cos, tanh, ...) at p bits: a few ulps.
        let direct = p.saturating_sub(8);
        let zero = Self::new2(p, Sign::Pos, false)?;
        let one = Self::from_word(1, p)?;
        if self.is_zero() {
            // rebuilt from the literals, one owner each
            return Ok((
                zero,
                Self::new2(p, Sign::Pos, false)?,
                Self::from_word(1, p)?,
                one,
                p,
            ));
        }
        if m.is_zero() {
            // am(u|0) = u, carried at the working precision like the other three results; a narrower value
            // makes the caller's roundability check (which assumes `p` working bits) overflow.
            let am = self.clone()?;
            let mut am = am;
            am.set_precision(p, RoundingMode::None)?;
            let sn = self.sin(p, RoundingMode::None, cc)?;
            let cn = self.cos(p, RoundingMode::None, cc)?;
            return Ok((am, sn, cn, one, direct));
        }
        if m.cmp(&one) == 0 {
            let sn = self.tanh(p, RoundingMode::None, cc)?;
            let ch = self.cosh(p, RoundingMode::None, cc)?;
            let cn = one.div(&ch, p, RoundingMode::None)?;
            // dn(u|1) = cn(u|1): the same quotient formed again, one owner each.
            let dn = one.div(&ch, p, RoundingMode::None)?;
            let sh = self.sinh(p, RoundingMode::None, cc)?;
            let am = sh.atan(p, RoundingMode::None, cc)?;
            return Ok((am, sn, cn, dn, direct));
        }
        let reduced = self.jacobi_reduce_period(m, p, cc)?;
        let u = match &reduced {
            Some((r, _)) => r,
            None => self,
        };
        let (mut am, sn, cn, dn, err) = jacobi_agm(u, reduced.is_some(), m, p, cc)?;
        // |u| < 2^-p: am ≈ sn ≈ u, cn ≈ dn ≈ 1, nothing cancels and every value is good to the
        // working precision (the bound below would be formed near the exponent floor).
        if reduced.is_none() && self.exponent() < -(p as Exponent) {
            return Ok((am, sn, cn, dn, direct));
        }
        // Absolute error bound of am (None: the recurrence could not bound it at p bits).
        let Some(mut err) = err else {
            return Ok((am, sn, cn, dn, 0));
        };
        if let Some((_, Some(n))) = reduced {
            // am(u + 4nK) = am(u) + 2nπ; sn, cn, dn have period 4K.
            let pe = round_p(
                p + usize::try_from(n.exponent().max(0)).map_err(|_| Error::InvalidArgument)?
                    + WORD_BIT_SIZE,
            );
            let mut two_pi = cc.pi_num(pe, RoundingMode::None)?;
            two_pi.set_exponent(two_pi.exponent() + 1);
            let off = n.mul(&two_pi, pe, RoundingMode::None)?;
            // |off - 2nπ| ≤ 2 ulp(off) at pe bits.
            err = err.add(
                &bound_pow2(off.exponent() + 2 - pe as Exponent)?,
                BOUND_P,
                RoundingMode::Up,
            )?;
            am = am.add(&off, p, RoundingMode::None)?;
            err = err.add(&bound_ulp(&am, p)?, BOUND_P, RoundingMode::Up)?;
        }
        let g = [
            rel_bits(&err, &am)?,
            rel_bits(&err, &sn)?,
            rel_bits(&err, &cn)?,
            rel_bits(&err, &dn)?,
        ]
        .into_iter()
        .min()
        .unwrap_or(0)
        .min(direct);
        Ok((am, sn, cn, dn, g))
    }

    /// `m = 1` with `cosh u` above the range: `sech u = 2e^{-|u|}` is below every subnormal,
    /// so `cn = dn` underflow, `sn = tanh u = ±(1 - tiny)` and `am = gd u = ±(π/2 - tiny)`.
    /// The four values are formed with `cn = dn` at the least normal (any quotient by them
    /// overflows, as the true one does) and the picked result is rounded once: an underflowed
    /// pick goes to zero or the least subnormal by the mode, an overflowed one is returned as
    /// `ExponentOverflow` for the caller to map by the mode. `None` when the regular path applies.
    fn jacobi_m1_far<F>(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
        pick: &F,
    ) -> Result<Option<Self>, Error>
    where
        F: Fn(Self, Self, Self, Self, usize) -> Result<Self, Error>,
    {
        if self.is_zero() || m.cmp(&Self::from_word(1, WORD_BIT_SIZE)?) != 0 {
            return Ok(None);
        }
        match self.cosh(WORD_BIT_SIZE, RoundingMode::None, cc) {
            Err(Error::ExponentOverflow(_)) => {}
            Ok(_) => return Ok(None),
            Err(e) => return Err(e),
        }
        let sign = if self.is_negative() { Sign::Neg } else { Sign::Pos };
        let p_x = round_p(p + 2 * WORD_BIT_SIZE);
        let mut am = cc.pi_num(p_x, RoundingMode::None)?;
        am.set_exponent(am.exponent() - 1);
        am.set_sign(sign);
        let one = Self::from_word(1, p_x)?;
        let mut eps = Self::from_word(1, p_x)?;
        // 2^-(p+64): below every ulp at p, yet visible to the p_x-bit quotients 1/sn, sn/dn.
        eps.set_exponent(1 - (p + WORD_BIT_SIZE) as Exponent);
        let mut sn = one.sub(&eps, p_x, RoundingMode::None)?;
        sn.set_sign(sign);
        let mut tiny = Self::from_word(1, p_x)?;
        tiny.set_exponent(EXPONENT_MIN);
        let mut cn = Self::from_word(1, p_x)?;
        cn.set_exponent(EXPONENT_MIN);
        let mut r = pick(am, sn, cn, tiny, p_x)?;
        if !r.is_zero() && r.exponent() <= EXPONENT_MIN + 2 {
            let neg = r.is_negative();
            r = Self::from_word(1, p)?;
            r.set_sign(if neg { Sign::Neg } else { Sign::Pos });
            r.set_exponent(EXPONENT_MIN);
            r.subnormalize(2 * EXPONENT_MIN, rm);
        } else {
            r.set_precision(p, rm)?;
        }
        r.set_inexact(true);
        Ok(Some(r))
    }

    /// `u` reduced into `[-2K, 2K]` and the number of whole periods `n` removed
    /// (`u = r + 4nK`; `n` is `None` when it rounds to zero). `None` when `|u| ≤ 2K` already, so
    /// the caller uses `u` itself. K is computed with `exponent(u)` extra
    /// bits, as the trig reduction does for π, so the reduced argument keeps all `p` bits.
    fn jacobi_reduce_period(
        &self,
        m: &Self,
        p: usize,
        cc: &mut Consts,
    ) -> Result<Option<(Self, Option<Self>)>, Error> {
        // \(K(m)\ge\pi/2\), so \(2K\ge\pi\). Skip \(K\) (and a nested Ziv) when \(\lvert u\rvert<\pi\).
        let pi = cc.pi_num(p, RoundingMode::None)?;
        if self.abs()?.cmp(&pi) < 0 {
            return Ok(None);
        }
        // K carries exponent(u) extra bits (checked, never capped: `MemoryAllocation` past what
        // the machine holds).
        let pk = round_p(reduction_width(p + WORD_BIT_SIZE, self.exponent())?);
        let k = m.elliptic_k_at(pk, cc)?;
        // |u| ≤ 2K as |u|/2 ≤ K (the exponent shift is exact; |u| ≥ π here).
        let mut half_u = self.abs()?;
        half_u.set_exponent(half_u.exponent() - 1);
        if half_u.cmp(&k) <= 0 {
            return Ok(None);
        }
        let mut four_k = k;
        four_k.set_exponent(four_k.exponent() + 2);
        // x > 2K as 2x > 4K (exact exponent shift, undone after the compare).
        let above_two_k = |x: &mut Self| -> bool {
            if x.is_zero() {
                return false;
            }
            x.set_exponent(x.exponent() + 1);
            let above = x.cmp(&four_k) > 0;
            x.set_exponent(x.exponent() - 1);
            above
        };
        let mut r = self.rem(&four_k)?;
        if above_two_k(&mut r) {
            r = r.sub(&four_k, pk, RoundingMode::None)?;
        } else if r.is_negative() && above_two_k(&mut r.abs()?) {
            r = r.add(&four_k, pk, RoundingMode::None)?;
        }
        // n = (u - r) / 4K is an integer up to the rounding error of r at pk bits.
        let n = self
            .sub(&r, pk, RoundingMode::None)?
            .div(&four_k, pk, RoundingMode::None)?
            .round(0, RoundingMode::ToEven)?;
        r.set_precision(p, RoundingMode::None)?;
        Ok(Some((r, if n.is_zero() { None } else { Some(n) })))
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
    let thresh = -((p as Exponent) / 2 + 16);
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

/// Precision of the running error bounds (rounded up, so they stay upper bounds).
const BOUND_P: usize = WORD_BIT_SIZE;

/// `2^e` at [`BOUND_P`] bits.
fn bound_pow2(e: Exponent) -> Result<ExactNumNumber, Error> {
    let mut x = ExactNumNumber::from_word(1, BOUND_P)?;
    x.set_exponent(e.saturating_add(1));
    Ok(x)
}

/// An upper bound of one ulp of `x` at `p` bits (`2^(exponent(x) - p)`; `2^-p` for zero).
fn bound_ulp(x: &ExactNumNumber, p: usize) -> Result<ExactNumNumber, Error> {
    let e = if x.is_zero() { 0 } else { x.exponent() };
    bound_pow2(e - p as Exponent)
}

/// `|x|` rounded up to [`BOUND_P`] bits.
fn bound_abs(x: &ExactNumNumber) -> Result<ExactNumNumber, Error> {
    let mut a = x.abs()?;
    a.set_precision(BOUND_P, RoundingMode::Up)?;
    Ok(a)
}

/// Bits certified by an absolute error bound `err` on a value `x`: `g` with `err/|x| < 2^-g`
/// (0 when `x` is zero or the bound reaches `|x|`).
fn rel_bits(err: &ExactNumNumber, x: &ExactNumNumber) -> Result<usize, Error> {
    if x.is_zero() {
        return Ok(0);
    }
    if err.is_zero() {
        return Ok(usize::MAX);
    }
    let mut ax = x.abs()?;
    ax.set_precision(BOUND_P, RoundingMode::Down)?;
    let r = err.div(&ax, BOUND_P, RoundingMode::Up)?;
    Ok(usize::try_from(-r.exponent()).unwrap_or(0))
}

/// am, sn, cn, dn by the AGM and the descending Landen recurrence (DLMF 22.20(ii)), with an
/// absolute error bound of am carried through the recurrence (`None` when it cannot be bounded
/// at `p` bits). A first-order bound: every working operation is within 2 ulp; the error of
/// `phi_n` reaches `phi_{n-1} = (phi_n + asin t_n)/2` through `t_n = c_n sin phi_n / a_n`,
/// amplified by `1/sqrt(1 - t_n^2)`. sn, cn inherit the absolute bound of am; dn is formed as
/// `sqrt(cn^2 + m' sn^2)` (both terms non-negative), so its relative error is that of cn or sn.
#[allow(clippy::type_complexity)]
fn jacobi_agm(
    u: &ExactNumNumber,
    reduced: bool,
    m: &ExactNumNumber,
    p: usize,
    cc: &mut Consts,
) -> Result<
    (
        ExactNumNumber,
        ExactNumNumber,
        ExactNumNumber,
        ExactNumNumber,
        Option<ExactNumNumber>,
    ),
    Error,
> {
    let one = ExactNumNumber::from_word(1, p)?;
    let two = ExactNumNumber::from_word(2, p)?;
    // m' = 1 - m: exact for m ≥ 1/2 (Sterbenz; p covers m's bits), one rounding otherwise.
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
    loop {
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
        // Quadratic AGM: a stall is a |c| that does not shrink.
        if !cc_seq[n as usize].is_zero() && c_next.abs()?.cmp(&cc_seq[n as usize].abs()?) >= 0 {
            return Err(Error::PrecisionRetryExhausted);
        }
        b = t;
        aa.push(a_next);
        cc_seq.push(c_next);
        n = n.checked_add(1).ok_or(Error::PrecisionRetryExhausted)?;
    }
    let last = aa.len() - 1;
    let shift = Exponent::try_from(last).map_err(|_| Error::InvalidArgument)?;
    let twon = ExactNumNumber::from_word(1, p)?.ldexp(shift, p, RoundingMode::None)?;
    let mut phi = twon
        .mul(&aa[last], p, RoundingMode::None)?
        .mul(u, p, RoundingMode::None)?;
    let wp = p as Exponent;
    let up = RoundingMode::Up;
    let steps = ExactNumNumber::from_word(4 * last as crate::defs::Word + 8, BOUND_P)?;
    let eps = bound_pow2(-wp)?;
    let steps_eps = steps.mul(&eps, BOUND_P, up)?;
    let bone = ExactNumNumber::from_word(1, BOUND_P)?;
    // phi_N = 2^N a_N u: a_N after N AGM steps (≤ (4N+8) ulp relative); a reduced argument
    // also carries an absolute error ≤ 2^(4-p) (the period removed at extra bits), scaled by 2^N a_N.
    // An unreduced u is the exact input.
    let mut err = bound_abs(&phi)?.mul(&steps_eps, BOUND_P, up)?;
    if reduced {
        err = err.add(
            &bound_abs(&twon)?
                .mul(&bound_abs(&aa[last])?, BOUND_P, up)?
                .mul(&bound_pow2(4 - wp)?, BOUND_P, up)?,
            BOUND_P,
            up,
        )?;
    }
    let mut bounded = true;
    let mut idx = last;
    while idx > 0 {
        let sphi = phi.sin(p, RoundingMode::None, cc)?;
        let ca = cc_seq[idx].div(&aa[idx], p, RoundingMode::None)?;
        let t_raw = ca.mul(&sphi, p, RoundingMode::None)?;
        let t = clamp_unit(&t_raw, p)?;
        let asint = t.asin(p, RoundingMode::None, cc)?;
        if bounded {
            // |cos phi_true| ≤ |cos phi_n| + err.
            let c = phi.cos(BOUND_P, up, cc)?;
            let cphi = bound_abs(&c)?
                .add(&bound_ulp(&c, BOUND_P - 1)?, BOUND_P, up)?
                .add(&err, BOUND_P, up)?;
            let dt = bound_abs(&ca)?
                .mul(&cphi, BOUND_P, up)?
                .mul(&err, BOUND_P, up)?
                .add(
                    // c_n carries an absolute error ≤ (4N+8) ulp of a_n (c = (a - b)/2 cancels):
                    // relative to |sin phi_n|, not to |t|.
                    &bound_abs(&t)?
                        .add(&bound_abs(&sphi)?, BOUND_P, up)?
                        .mul(&steps_eps, BOUND_P, up)?,
                    BOUND_P,
                    up,
                )?;
            // sqrt(1 - t^2) ≥ sqrt(1 - |t_true|) with |t_true| ≤ |t_raw| + dt.
            // 1 - |t_raw| from the full-width t_raw (a 64-bit |t_raw| would round to 1).
            let gap = bone.sub(&t_raw.abs()?, BOUND_P, RoundingMode::Down)?.sub(
                &dt,
                BOUND_P,
                RoundingMode::Down,
            )?;
            if gap.is_zero() || gap.is_negative() {
                bounded = false;
            } else {
                let den = gap.sqrt(BOUND_P, RoundingMode::Down)?;
                let da = dt
                    .div(&den, BOUND_P, up)?
                    .add(&bound_ulp(&asint, p - 1)?, BOUND_P, up)?;
                err = err.add(&da, BOUND_P, up)?;
                err.div_by_2(up);
            }
        }
        phi = asint
            .add(&phi, p, RoundingMode::None)?
            .div(&two, p, RoundingMode::None)?;
        if bounded {
            err = err.add(&bound_ulp(&phi, p - 1)?, BOUND_P, up)?;
        }
        idx -= 1;
    }
    let sn = phi.sin(p, RoundingMode::None, cc)?;
    let cn = phi.cos(p, RoundingMode::None, cc)?;
    // dn = sqrt(cn^2 + m' sn^2) (DLMF 22.6.1), no cancellation near the zeros of cn and dn.
    let dn = cn
        .mul(&cn, p, RoundingMode::None)?
        .add(
            &om.mul(&sn.mul(&sn, p, RoundingMode::None)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?
        .sqrt(p, RoundingMode::None)?;
    // The final roundings of sin, cos, dn (relative, a few ulp) are inside the caller's `p - 8`.
    let err = if bounded { Some(err) } else { None };
    Ok((phi, sn, cn, dn, err))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::defs::RoundingMode;

    fn bits_agree(
        a: &ExactNumNumber,
        b: &ExactNumNumber,
        p: usize,
        min_bits: Exponent,
        label: &str,
    ) -> Result<(), crate::Error> {
        let d = a.sub(b, p, RoundingMode::None)?.abs()?;
        if d.is_zero() {
            return Ok(());
        }
        let rel = d.exponent() - a.exponent();
        assert!(
            rel < -min_bits,
            "{label}: relative exponent {rel} (a.exp={}, b.exp={})",
            a.exponent(),
            b.exponent()
        );
        Ok(())
    }

    #[test]
    fn test_jacobi_identities() -> Result<(), crate::Error> {
        let p = 256;
        let mut cc = Consts::new()?;
        let rm = RoundingMode::ToEven;
        let one = ExactNumNumber::from_word(1, p)?;
        let zero = ExactNumNumber::new(p)?;
        let two = ExactNumNumber::from_word(2, p)?;
        let mut half = ExactNumNumber::from_word(1, p)?;
        half.div_by_2(RoundingMode::None);

        let sn0 = zero.jacobi_sn(&half, p, rm, &mut cc)?;
        assert!(sn0.is_zero(), "sn(0)");
        let cn0 = zero.jacobi_cn(&half, p, rm, &mut cc)?;
        bits_agree(&cn0, &one, p, 80, "cn(0)")?;
        let dn0 = zero.jacobi_dn(&half, p, rm, &mut cc)?;
        bits_agree(&dn0, &one, p, 80, "dn(0)")?;
        let am0 = zero.jacobi_am(&half, p, rm, &mut cc)?;
        assert!(am0.is_zero(), "am(0)");

        let u = ExactNumNumber::from_word(1, p)?;
        let sn_m0 = u.jacobi_sn(&zero, p, rm, &mut cc)?;
        let sinu = u.sin(p, rm, &mut cc)?;
        bits_agree(&sn_m0, &sinu, p, 80, "sn(u|0)=sin u")?;
        let cn_m0 = u.jacobi_cn(&zero, p, rm, &mut cc)?;
        let cosu = u.cos(p, rm, &mut cc)?;
        bits_agree(&cn_m0, &cosu, p, 80, "cn(u|0)=cos u")?;
        let dn_m0 = u.jacobi_dn(&zero, p, rm, &mut cc)?;
        bits_agree(&dn_m0, &one, p, 80, "dn(u|0)=1")?;

        let sn_m1 = u.jacobi_sn(&one, p, rm, &mut cc)?;
        let tanhu = u.tanh(p, rm, &mut cc)?;
        bits_agree(&sn_m1, &tanhu, p, 80, "sn(u|1)=tanh u")?;
        let cn_m1 = u.jacobi_cn(&one, p, rm, &mut cc)?;
        let sech = one.div(&u.cosh(p, rm, &mut cc)?, p, RoundingMode::None)?;
        bits_agree(&cn_m1, &sech, p, 80, "cn(u|1)=sech u")?;
        let dn_m1 = u.jacobi_dn(&one, p, rm, &mut cc)?;
        bits_agree(&dn_m1, &sech, p, 80, "dn(u|1)=sech u")?;

        let sn = u.jacobi_sn(&half, p, rm, &mut cc)?;
        let cn = u.jacobi_cn(&half, p, rm, &mut cc)?;
        let dn = u.jacobi_dn(&half, p, rm, &mut cc)?;
        let s2c2 = sn.mul(&sn, p, RoundingMode::None)?.add(
            &cn.mul(&cn, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
        bits_agree(&s2c2, &one, p, 80, "sn²+cn²=1")?;
        let d2ms = dn.mul(&dn, p, RoundingMode::None)?.add(
            &half.mul(&sn.mul(&sn, p, RoundingMode::None)?, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
        bits_agree(&d2ms, &one, p, 80, "dn²+m sn²=1")?;

        let k = half.elliptic_k(p, rm, &mut cc)?;
        let mut kh = half.elliptic_k(p, rm, &mut cc)?;
        kh.div_by_2(RoundingMode::None);
        let snkh = kh.jacobi_sn(&half, p, rm, &mut cc)?;
        let kp = one
            .sub(&half, p, RoundingMode::None)?
            .sqrt(p, RoundingMode::None)?;
        let want = one.div(
            &one.add(&kp, p, RoundingMode::None)?
                .sqrt(p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
        bits_agree(&snkh, &want, p, 40, "sn(K/2|1/2)")?;

        let four_k = four_word(p)?.mul(&k, p, RoundingMode::None)?;
        let u4 = u.add(&four_k, p, RoundingMode::None)?;
        let sn_per = u4.jacobi_sn(&half, p, rm, &mut cc)?;
        bits_agree(&sn_per, &sn, p, 40, "sn(u+4K)=sn(u)")?;

        let mut h = ExactNumNumber::from_word(1, p)?;
        for _ in 0..20 {
            h.div_by_2(RoundingMode::None);
        }
        let up = u.add(&h, p, RoundingMode::None)?;
        let um = u.sub(&h, p, RoundingMode::None)?;
        let snp = up.jacobi_sn(&half, p, rm, &mut cc)?;
        let snm = um.jacobi_sn(&half, p, rm, &mut cc)?;
        let fd = snp.sub(&snm, p, RoundingMode::None)?.div(
            &two.mul(&h, p, RoundingMode::None)?,
            p,
            RoundingMode::None,
        )?;
        let deriv = cn.mul(&dn, p, RoundingMode::None)?;
        bits_agree(&fd, &deriv, p, 20, "d(sn)/du = cn dn")?;

        let sn_none = u.jacobi_sn(&half, 512, RoundingMode::None, &mut cc)?;
        bits_agree(&sn_none, &sn, p, 20, "sn None 512 vs ToEven 256")?;

        assert!(u.jacobi_sn(&one.neg()?, p, rm, &mut cc).is_err());
        assert!(u.jacobi_sn(&two, p, rm, &mut cc).is_err());
        let cd = u.jacobi_cd(&half, p, rm, &mut cc)?;
        let want_cd = cn.div(&dn, p, RoundingMode::None)?;
        bits_agree(&cd, &want_cd, p, 80, "cd=cn/dn")?;
        let ns = u.jacobi_ns(&half, p, rm, &mut cc)?;
        let nssn = ns.mul(&sn, p, RoundingMode::None)?;
        bits_agree(&nssn, &one, p, 80, "ns·sn=1")?;
        assert!(zero.jacobi_ns(&half, p, rm, &mut cc).is_err());

        // Inverse of sn is incomplete F: F(sn(u|m)|m)=u for |u|<K(m).
        let back = sn.elliptic_f(&half, p, rm, &mut cc)?;
        bits_agree(&back, &u, p, 40, "F(sn(u|m)|m)=u")?;
        Ok(())
    }

    fn four_word(p: usize) -> Result<ExactNumNumber, crate::Error> {
        ExactNumNumber::from_word(4, p)
    }
}
