//! Meijer \(G\), Fox \(H\), and local / confluent Heun on [`ExactNum`].
//!
//! SoftFloat only. Residue \(G\) follows DLMF 16.17 / mpmath `hypercomb`; Fox \(H\)
//! with positive rational \(A,B\) reduces through the Gauss multiplication formula
//! (mpmath `foxh`). Local Heun is DLMF 31.3; confluent Heun is the \(z=0\)
//! Frobenius solution of DLMF 31.12.1.

use crate::catalog::{cat_nan, finish, finite, series_cap, small_int, work_p};
use crate::Consts;
use crate::ExactNum;
use crate::RoundingMode;
use alloc::vec::Vec;

/// Max length of any one Meijer / Fox parameter group.
const PARAM_MAX: usize = 8;
/// Max positive integer scale \(A_i,B_j\) after clearing denominators.
const SCALE_MAX: i32 = 12;
/// Max expanded Meijer length after the Fox–Gauss lift.
const TILDE_MAX: usize = 16;

fn gamma_pole(x: &ExactNum) -> bool {
    x.is_int() && !x.is_positive()
}

fn all_finite(xs: &[ExactNum]) -> bool {
    xs.iter().all(finite)
}

fn pairs_finite(xs: &[(ExactNum, ExactNum)]) -> bool {
    xs.iter().all(|(a, a_s)| finite(a) && finite(a_s))
}

fn pow_real(z: &ExactNum, e: &ExactNum, p: usize, cc: &mut Consts) -> ExactNum {
    if z.is_zero() {
        if e.is_positive() {
            return ExactNum::new(p);
        }
        if e.is_zero() {
            return ExactNum::from_u8(1, p);
        }
        return cat_nan();
    }
    if z.is_negative() {
        if let Some(n) = small_int(e, p, -64, 64) {
            return z.powsi(n as isize, p, RoundingMode::None);
        }
        return cat_nan();
    }
    z.pow(e, p, RoundingMode::None, cc)
}

fn gamma_ratio(
    gn: &[ExactNum],
    gd: &[ExactNum],
    p: usize,
    cc: &mut Consts,
) -> Result<Option<ExactNum>, ()> {
    if gn.iter().any(gamma_pole) {
        return Err(());
    }
    if gd.iter().any(gamma_pole) {
        return Ok(None);
    }
    let mut num = ExactNum::from_u8(1, p);
    for x in gn {
        let g = x.gamma(p, RoundingMode::None, cc);
        if !finite(&g) {
            return Err(());
        }
        num = num.mul(&g, p, RoundingMode::None);
    }
    let mut den = ExactNum::from_u8(1, p);
    for x in gd {
        let g = x.gamma(p, RoundingMode::None, cc);
        if !finite(&g) {
            return Ok(None);
        }
        den = den.mul(&g, p, RoundingMode::None);
    }
    Ok(Some(num.div(&den, p, RoundingMode::None)))
}

fn odd_i(n: i32) -> bool {
    n.rem_euclid(2) == 1
}

fn abs_gt_one(z: &ExactNum, p: usize) -> bool {
    matches!(z.abs().cmp(&ExactNum::from_u8(1, p)), Some(c) if c > 0)
}

fn choose_series(an_n: usize, ap_n: usize, bm_n: usize, bq_n: usize, z: &ExactNum, p: usize) -> u8 {
    let pn = an_n + ap_n;
    let qn = bm_n + bq_n;
    // Type 1 is {}_p F_{q-1}; type 2 is {}_q F_{p-1}. When p=q both are
    // disk series, so |z| picks the side (we do not continue {}_r F_{r-1}).
    let mut s = if pn < qn {
        1u8
    } else if pn > qn {
        2
    } else if abs_gt_one(z, p) {
        2
    } else {
        1
    };
    if s == 1 && bm_n == 0 {
        s = 2;
    }
    if s == 2 && an_n == 0 {
        s = 1;
    }
    s
}

fn meijer_g_inner(
    z: &ExactNum,
    an: &[ExactNum],
    ap: &[ExactNum],
    bm: &[ExactNum],
    bq: &[ExactNum],
    r: &ExactNum,
    p: usize,
    cc: &mut Consts,
) -> ExactNum {
    if an.len() + ap.len() > PARAM_MAX || bm.len() + bq.len() > PARAM_MAX {
        return cat_nan();
    }
    if !finite(z)
        || !finite(r)
        || !r.is_positive()
        || !all_finite(an)
        || !all_finite(ap)
        || !all_finite(bm)
        || !all_finite(bq)
    {
        return cat_nan();
    }
    if an.is_empty() && bm.is_empty() {
        return cat_nan();
    }
    let mut a = Vec::with_capacity(an.len() + ap.len());
    a.extend_from_slice(an);
    a.extend_from_slice(ap);
    let mut b = Vec::with_capacity(bm.len() + bq.len());
    b.extend_from_slice(bm);
    b.extend_from_slice(bq);
    let n = an.len();
    let m = bm.len();
    let pn = a.len();
    let qn = b.len();
    let series = choose_series(n, ap.len(), m, bq.len(), z, p);
    let one = ExactNum::from_u8(1, p);
    let r_is_one = r.cmp(&one) == Some(0);
    if !r_is_one && !z.is_positive() && !z.is_zero() {
        return cat_nan();
    }
    let z_to_inv_r = if r_is_one {
        z.clone()
    } else {
        pow_real(z, &one.div(r, p, RoundingMode::None), p, cc)
    };
    if !finite(&z_to_inv_r) {
        return cat_nan();
    }
    let mut sum = ExactNum::new(p);
    if series == 1 {
        let sign_odd = odd_i((pn as i32) - (m as i32) - (n as i32));
        let hz = if sign_odd { z_to_inv_r.neg() } else { z_to_inv_r.clone() };
        for h in 0..m {
            let bh = &b[h];
            let mut gn = Vec::new();
            for (j, bj) in b.iter().enumerate().take(m) {
                if j != h {
                    gn.push(bj.sub(bh, p, RoundingMode::None));
                }
            }
            for aj in a.iter().take(n) {
                gn.push(
                    one.sub(aj, p, RoundingMode::None)
                        .add(bh, p, RoundingMode::None),
                );
            }
            let mut gd = Vec::new();
            for aj in a.iter().skip(n) {
                gd.push(aj.sub(bh, p, RoundingMode::None));
            }
            for bj in b.iter().skip(m) {
                gd.push(
                    one.sub(bj, p, RoundingMode::None)
                        .add(bh, p, RoundingMode::None),
                );
            }
            let ratio = match gamma_ratio(&gn, &gd, p, cc) {
                Err(()) => return cat_nan(),
                Ok(None) => continue,
                Ok(Some(v)) => v,
            };
            let mut hn = Vec::with_capacity(pn);
            for aj in &a {
                hn.push(
                    one.sub(aj, p, RoundingMode::None)
                        .add(bh, p, RoundingMode::None),
                );
            }
            let mut hd = Vec::with_capacity(qn.saturating_sub(1));
            for (j, bj) in b.iter().enumerate() {
                if j != h {
                    hd.push(
                        one.sub(bj, p, RoundingMode::None)
                            .add(bh, p, RoundingMode::None),
                    );
                }
            }
            let f = hz.hypergeom_pfq(&hn, &hd, p, RoundingMode::None, cc);
            if !finite(&f) {
                return cat_nan();
            }
            let expn = if r_is_one { bh.clone() } else { bh.div(r, p, RoundingMode::None) };
            let zp = pow_real(z, &expn, p, cc);
            if !finite(&zp) {
                return cat_nan();
            }
            sum = sum.add(
                &ratio
                    .mul(&zp, p, RoundingMode::None)
                    .mul(&f, p, RoundingMode::None),
                p,
                RoundingMode::None,
            );
        }
    } else {
        let sign_odd = odd_i((qn as i32) - (m as i32) - (n as i32));
        let inv_z = if z_to_inv_r.is_zero() {
            return cat_nan();
        } else {
            one.div(&z_to_inv_r, p, RoundingMode::None)
        };
        let hz = if sign_odd { inv_z.neg() } else { inv_z };
        for k in 0..n {
            let ak = &a[k];
            let mut gn = Vec::new();
            for (j, aj) in a.iter().enumerate().take(n) {
                if j != k {
                    gn.push(ak.sub(aj, p, RoundingMode::None));
                }
            }
            for bj in b.iter().take(m) {
                gn.push(
                    one.sub(ak, p, RoundingMode::None)
                        .add(bj, p, RoundingMode::None),
                );
            }
            let mut gd = Vec::new();
            for bj in b.iter().skip(m) {
                gd.push(ak.sub(bj, p, RoundingMode::None));
            }
            for aj in a.iter().skip(n) {
                gd.push(
                    one.sub(ak, p, RoundingMode::None)
                        .add(aj, p, RoundingMode::None),
                );
            }
            let ratio = match gamma_ratio(&gn, &gd, p, cc) {
                Err(()) => return cat_nan(),
                Ok(None) => continue,
                Ok(Some(v)) => v,
            };
            let mut hn = Vec::with_capacity(qn);
            for bj in &b {
                hn.push(
                    one.sub(ak, p, RoundingMode::None)
                        .add(bj, p, RoundingMode::None),
                );
            }
            let mut hd = Vec::with_capacity(pn.saturating_sub(1));
            for (j, aj) in a.iter().enumerate() {
                if j != k {
                    hd.push(
                        one.add(aj, p, RoundingMode::None)
                            .sub(ak, p, RoundingMode::None),
                    );
                }
            }
            let f = hz.hypergeom_pfq(&hn, &hd, p, RoundingMode::None, cc);
            if !finite(&f) {
                return cat_nan();
            }
            let expn = if r_is_one {
                ak.sub(&one, p, RoundingMode::None)
            } else {
                ak.sub(&one, p, RoundingMode::None)
                    .div(r, p, RoundingMode::None)
            };
            let zp = pow_real(z, &expn, p, cc);
            if !finite(&zp) {
                return cat_nan();
            }
            sum = sum.add(
                &ratio
                    .mul(&zp, p, RoundingMode::None)
                    .mul(&f, p, RoundingMode::None),
                p,
                RoundingMode::None,
            );
        }
    }
    sum
}

fn gcd_i(mut a: i32, mut b: i32) -> i32 {
    a = a.unsigned_abs() as i32;
    b = b.unsigned_abs() as i32;
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

fn lcm_i(a: i32, b: i32) -> Option<i32> {
    if a == 0 || b == 0 {
        return Some(0);
    }
    let g = gcd_i(a, b);
    a.checked_div(g)?.checked_mul(b)?.checked_abs()
}

fn near_small_int(x: &ExactNum, p: usize, lo: i32, hi: i32) -> Option<i32> {
    if let Some(n) = small_int(x, p, lo, hi) {
        return Some(n);
    }
    for n in lo..=hi {
        let d = x
            .sub(&ExactNum::from_i32(n, p), p, RoundingMode::None)
            .abs();
        if d.is_zero()
            || d.exponent()
                .is_some_and(|e| (e as i64) + ((p as i64) / 2) < 0)
        {
            return Some(n);
        }
    }
    None
}

fn as_pos_rational(x: &ExactNum, p: usize) -> Option<(i32, i32)> {
    if !finite(x) || !x.is_positive() {
        return None;
    }
    for den in 1..=16 {
        let prod = x.mul(&ExactNum::from_u32(den, p), p, RoundingMode::None);
        if let Some(num) = near_small_int(&prod, p, 1, 48) {
            let d = den as i32;
            let g = gcd_i(num, d);
            return Some((num / g, d / g));
        }
    }
    None
}

fn all_unit_scales(xs: &[(ExactNum, ExactNum)], p: usize) -> bool {
    let one = ExactNum::from_u8(1, p);
    xs.iter().all(|(_, s)| s.cmp(&one) == Some(0))
}

fn pow_i32(base: &ExactNum, n: i32, p: usize, _cc: &mut Consts) -> ExactNum {
    if n >= 0 {
        base.powsi(n as isize, p, RoundingMode::None)
    } else {
        let y = base.powsi((-n) as isize, p, RoundingMode::None);
        ExactNum::from_u8(1, p).div(&y, p, RoundingMode::None)
    }
}

impl ExactNum {
    /// Meijer \(G^{m,n}_{p,q}(\mathrm{self}\mid a_1,\ldots,a_n;a_{n+1},\ldots,a_p; b_1,\ldots,b_m;b_{m+1},\ldots,b_q)\).
    ///
    /// Groups match mpmath `meijerg([[an],[ap]],[[bm],[bq]], z)` and DLMF 16.17.
    /// Residue series of type 1 (poles \(b_h+\mathbb{N}_0\)) or type 2 (poles \(a_k-1-\mathbb{N}_0\)).
    ///
    /// # Precision
    ///
    /// - Algorithm: DLMF 16.17 residue times [`ExactNum::hypergeom_pfq`]. Type 1 when
    ///   \(p<q\) or \(p=q\) and \(\lvert z\rvert\le 1\); type 2 otherwise. Extra word
    ///   of working precision.
    /// - Bound: working precision `p + 2 WORD_BIT_SIZE`; not Ziv-certified.
    ///
    /// # Limitations
    ///
    /// - Real line only. \(\mathrm{self}<0\) is `NaN` unless every used power is an integer.
    /// - Coincident poles (\(b_j-b_h\) a non-positive integer among the first \(m\)
    ///   parameters, or a numerator \(\Gamma\) pole) return `NaN`. No logarithmic
    ///   residue / `hypercomb` limit is taken.
    /// - At most 8 parameters in \(a\) and in \(b\).
    pub fn meijer_g(
        &self,
        an: &[Self],
        ap: &[Self],
        bm: &[Self],
        bq: &[Self],
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        let wrk = work_p(p);
        let one = ExactNum::from_u8(1, wrk);
        finish(meijer_g_inner(self, an, ap, bm, bq, &one, wrk, cc), p, rm)
    }

    /// Fox \(H^{m,n}_{p,q}(\mathrm{self}\mid (a_i,A_i); (b_j,B_j))\).
    ///
    /// Each slice is pairs \((a,A)\) or \((b,B)\) in the same \(n,p,m,q\) grouping as
    /// [`meijer_g`]. When every \(A_i=B_j=1\) this is Meijer \(G\). Positive rational
    /// scales are lifted by the Gauss multiplication formula (mpmath `foxh`).
    ///
    /// # Precision
    ///
    /// - Algorithm: unit scales call [`ExactNum::meijer_g`]; otherwise Gauss
    ///   multiplication to a larger \(G\) (mpmath `foxh`). Extra word of working
    ///   precision.
    /// - Bound: same as Meijer \(G\) after the lift.
    ///
    /// # Limitations
    ///
    /// - \(A_i,B_j\) must be positive rationals with small numerator/denominator
    ///   (denominators \(\le 16\), cleared scales \(\le 12\)). Otherwise `NaN`.
    /// - Inherits Meijer coincident-pole `NaN` after the lift.
    /// - Real \(\mathrm{self}>0\) when a scale is not 1 (fractional \(z^{1/r}\)).
    pub fn fox_h(
        &self,
        an: &[(Self, Self)],
        ap: &[(Self, Self)],
        bm: &[(Self, Self)],
        bq: &[(Self, Self)],
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        let wrk = work_p(p);
        if an.len() + ap.len() > PARAM_MAX || bm.len() + bq.len() > PARAM_MAX {
            return cat_nan();
        }
        if !finite(self)
            || !pairs_finite(an)
            || !pairs_finite(ap)
            || !pairs_finite(bm)
            || !pairs_finite(bq)
        {
            return cat_nan();
        }
        if an.is_empty() && bm.is_empty() {
            return cat_nan();
        }
        if all_unit_scales(an, wrk)
            && all_unit_scales(ap, wrk)
            && all_unit_scales(bm, wrk)
            && all_unit_scales(bq, wrk)
        {
            let an_a: Vec<ExactNum> = an.iter().map(|(a, _)| a.clone()).collect();
            let ap_a: Vec<ExactNum> = ap.iter().map(|(a, _)| a.clone()).collect();
            let bm_a: Vec<ExactNum> = bm.iter().map(|(a, _)| a.clone()).collect();
            let bq_a: Vec<ExactNum> = bq.iter().map(|(a, _)| a.clone()).collect();
            return finish(
                meijer_g_inner(
                    self,
                    &an_a,
                    &ap_a,
                    &bm_a,
                    &bq_a,
                    &ExactNum::from_u8(1, wrk),
                    wrk,
                    cc,
                ),
                p,
                rm,
            );
        }
        fox_h_gauss(self, an, ap, bm, bq, wrk, cc)
            .map(|y| finish(y, p, rm))
            .unwrap_or_else(cat_nan)
    }

    /// Local Heun \(\mathrm{Hl}(a,q;\alpha,\beta,\gamma,\delta;\mathrm{self})\) (DLMF 31.3).
    ///
    /// Power series about \(z=0\) with \(c_0=1\), \(c_1=q/(a\gamma)\), and the
    /// three-term recurrence 31.3.3–4. \(\varepsilon=\alpha+\beta-\gamma-\delta+1\).
    ///
    /// # Precision
    ///
    /// - Algorithm: DLMF 31.3 power series, stopped at a relative term of \(2^{-p}\)
    ///   or `series_cap`. Extra word of working precision.
    /// - Bound: working precision `p + 2 WORD_BIT_SIZE`; not Ziv-certified.
    ///
    /// # Limitations
    ///
    /// - Converges for \(\lvert z\rvert<1\). \(\lvert z\rvert\ge 1\) is `NaN` (no
    ///   connection formulas).
    /// - \(a=0\) or \(\gamma\in\{0,-1,-2,\ldots\}\) (recurrence pole) is `NaN`.
    /// - When \(q=a\alpha\beta\) and \(\delta=\alpha+\beta-\gamma+1\), this is
    ///   \({}_2F_1(\alpha,\beta;\gamma;z)\).
    pub fn heun_g(
        &self,
        a: &Self,
        q: &Self,
        alpha: &Self,
        beta: &Self,
        gamma: &Self,
        delta: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        let _ = cc;
        if !finite(self)
            || !finite(a)
            || !finite(q)
            || !finite(alpha)
            || !finite(beta)
            || !finite(gamma)
            || !finite(delta)
        {
            return cat_nan();
        }
        let wrk = work_p(p);
        if matches!(self.abs().cmp(&ExactNum::from_u8(1, wrk)), Some(c) if c >= 0) {
            return cat_nan();
        }
        if a.is_zero() || gamma_pole(gamma) {
            return cat_nan();
        }
        if self.is_zero() {
            return finish(ExactNum::from_u8(1, wrk), p, rm);
        }
        let one = ExactNum::from_u8(1, wrk);
        let eps = alpha
            .add(beta, wrk, RoundingMode::None)
            .sub(gamma, wrk, RoundingMode::None)
            .sub(delta, wrk, RoundingMode::None)
            .add(&one, wrk, RoundingMode::None);
        let mut c_prev = ExactNum::from_u8(1, wrk);
        let mut c = q.div(
            &a.mul(gamma, wrk, RoundingMode::None),
            wrk,
            RoundingMode::None,
        );
        let mut zp = self.clone();
        let mut sum = c_prev.add(
            &c.mul(&zp, wrk, RoundingMode::None),
            wrk,
            RoundingMode::None,
        );
        for j in 1..=series_cap(wrk) {
            let jf = ExactNum::from_u32(j as u32, wrk);
            let jm1 = ExactNum::from_u32((j - 1) as u32, wrk);
            let jp1 = ExactNum::from_u32((j + 1) as u32, wrk);
            let pj = jm1.add(alpha, wrk, RoundingMode::None).mul(
                &jm1.add(beta, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            );
            // Q_j = j((j-1+γ)(1+a) + aδ + ε)
            let qj = jf.mul(
                &jm1.add(gamma, wrk, RoundingMode::None)
                    .mul(
                        &one.add(a, wrk, RoundingMode::None),
                        wrk,
                        RoundingMode::None,
                    )
                    .add(
                        &a.mul(delta, wrk, RoundingMode::None),
                        wrk,
                        RoundingMode::None,
                    )
                    .add(&eps, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            );
            let rj = a.mul(&jp1, wrk, RoundingMode::None).mul(
                &jf.add(gamma, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            );
            if rj.is_zero() || !finite(&rj) {
                return cat_nan();
            }
            let c_next = qj
                .add(q, wrk, RoundingMode::None)
                .mul(&c, wrk, RoundingMode::None)
                .sub(
                    &pj.mul(&c_prev, wrk, RoundingMode::None),
                    wrk,
                    RoundingMode::None,
                )
                .div(&rj, wrk, RoundingMode::None);
            zp = zp.mul(self, wrk, RoundingMode::None);
            let term = c_next.mul(&zp, wrk, RoundingMode::None);
            sum = sum.add(&term, wrk, RoundingMode::None);
            if term.is_zero()
                || term.abs().exponent().is_some_and(|e| {
                    (e as i64) + (p as i64) < sum.abs().exponent().unwrap_or(0) as i64
                })
            {
                break;
            }
            c_prev = c;
            c = c_next;
        }
        finish(sum, p, rm)
    }

    /// Confluent Heun: the \(y(0)=1\) Frobenius solution of DLMF 31.12.1,
    ///
    /// \[
    /// y''+\Bigl(\frac{\gamma}{z}+\frac{\delta}{z-1}+\varepsilon\Bigr)y'
    /// +\frac{\alpha z-q}{z(z-1)}y=0.
    /// \]
    ///
    /// Series \(\sum b_n z^n\) with \(b_0=1\), \(b_1=-q/\gamma\), and
    /// \((n+1)(n+\gamma)b_{n+1}=[n(n-1+\gamma+\delta-\varepsilon)-q]b_n+[\alpha+\varepsilon(n-1)]b_{n-1}\).
    ///
    /// Maple `HeunC(α,β,γ,δ,η,z)` is a different parameterization — not this method.
    ///
    /// # Precision
    ///
    /// - Algorithm: Frobenius series of DLMF 31.12.1 about \(z=0\), stopped at a
    ///   relative term of \(2^{-p}\) or `series_cap`. Extra word of working precision.
    /// - Bound: working precision `p + 2 WORD_BIT_SIZE`; not Ziv-certified.
    ///
    /// # Limitations
    ///
    /// - \(\lvert z\rvert<1\). \(\lvert z\rvert\ge 1\) is `NaN`.
    /// - \(\gamma\in\{0,-1,-2,\ldots\}\) is `NaN`.
    /// - When \(\delta=0\), \(\varepsilon=-1\), \(q=\alpha\), this is
    ///   \({}_1F_1(-\alpha;\gamma;z)\).
    /// - Biconfluent / double-confluent / triconfluent Heun (DLMF 31.12.2–4) are
    ///   not implemented: Accumath has no public SoftFloat signature for them,
    ///   and the \(z=0\) series is not unique (triconfluent) or not regular
    ///   (double-confluent).
    pub fn heun_c(
        &self,
        alpha: &Self,
        gamma: &Self,
        delta: &Self,
        eps: &Self,
        q: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        let _ = cc;
        if !finite(self)
            || !finite(alpha)
            || !finite(gamma)
            || !finite(delta)
            || !finite(eps)
            || !finite(q)
        {
            return cat_nan();
        }
        let wrk = work_p(p);
        if matches!(self.abs().cmp(&ExactNum::from_u8(1, wrk)), Some(c) if c >= 0) {
            return cat_nan();
        }
        if gamma_pole(gamma) {
            return cat_nan();
        }
        if self.is_zero() {
            return finish(ExactNum::from_u8(1, wrk), p, rm);
        }
        let mut b_prev = ExactNum::from_u8(1, wrk);
        let mut b = q.neg().div(gamma, wrk, RoundingMode::None);
        let mut zp = self.clone();
        let mut sum = b_prev.add(
            &b.mul(&zp, wrk, RoundingMode::None),
            wrk,
            RoundingMode::None,
        );
        for n in 1..=series_cap(wrk) {
            let nf = ExactNum::from_u32(n as u32, wrk);
            let nm1 = ExactNum::from_u32((n - 1) as u32, wrk);
            let np1 = ExactNum::from_u32((n + 1) as u32, wrk);
            let left = np1.mul(
                &nf.add(gamma, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            );
            if left.is_zero() || !finite(&left) {
                return cat_nan();
            }
            // n(n-1+γ+δ-ε)-q
            let coef_n = nf
                .mul(
                    &nm1.add(gamma, wrk, RoundingMode::None)
                        .add(delta, wrk, RoundingMode::None)
                        .sub(eps, wrk, RoundingMode::None),
                    wrk,
                    RoundingMode::None,
                )
                .sub(q, wrk, RoundingMode::None);
            let coef_nm1 = alpha.add(
                &eps.mul(&nm1, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            );
            let b_next = coef_n
                .mul(&b, wrk, RoundingMode::None)
                .add(
                    &coef_nm1.mul(&b_prev, wrk, RoundingMode::None),
                    wrk,
                    RoundingMode::None,
                )
                .div(&left, wrk, RoundingMode::None);
            zp = zp.mul(self, wrk, RoundingMode::None);
            let term = b_next.mul(&zp, wrk, RoundingMode::None);
            sum = sum.add(&term, wrk, RoundingMode::None);
            if term.is_zero()
                || term.abs().exponent().is_some_and(|e| {
                    (e as i64) + (p as i64) < sum.abs().exponent().unwrap_or(0) as i64
                })
            {
                break;
            }
            b_prev = b;
            b = b_next;
        }
        finish(sum, p, rm)
    }
}

fn fox_h_gauss(
    z: &ExactNum,
    an: &[(ExactNum, ExactNum)],
    ap: &[(ExactNum, ExactNum)],
    bm: &[(ExactNum, ExactNum)],
    bq: &[(ExactNum, ExactNum)],
    p: usize,
    cc: &mut Consts,
) -> Option<ExactNum> {
    let n = an.len();
    let m = bm.len();
    let pn = n + ap.len();
    let qn = m + bq.len();
    let mut a: Vec<ExactNum> = Vec::new();
    let mut a_frac: Vec<(i32, i32)> = Vec::new();
    for (v, s) in an.iter().chain(ap.iter()) {
        a.push(v.clone());
        a_frac.push(as_pos_rational(s, p)?);
    }
    let mut b: Vec<ExactNum> = Vec::new();
    let mut b_frac: Vec<(i32, i32)> = Vec::new();
    for (v, s) in bm.iter().chain(bq.iter()) {
        b.push(v.clone());
        b_frac.push(as_pos_rational(s, p)?);
    }
    let mut d = 1i32;
    for &(_, den) in a_frac.iter().chain(b_frac.iter()) {
        d = lcm_i(d, den)?;
        if d > 48 {
            return None;
        }
    }
    let mut a_int: Vec<i32> = Vec::new();
    for &(num, den) in &a_frac {
        let ai = num.checked_mul(d / den)?;
        if ai <= 0 || ai > SCALE_MAX {
            return None;
        }
        a_int.push(ai);
    }
    let mut b_int: Vec<i32> = Vec::new();
    for &(num, den) in &b_frac {
        let bi = num.checked_mul(d / den)?;
        if bi <= 0 || bi > SCALE_MAX {
            return None;
        }
        b_int.push(bi);
    }
    let mut a_tilde = Vec::new();
    for (ai, &ai_int) in a.iter().zip(a_int.iter()) {
        let aif = ExactNum::from_i32(ai_int, p);
        for k in 0..ai_int {
            let num = ai.add(&ExactNum::from_i32(k, p), p, RoundingMode::None);
            a_tilde.push(num.div(&aif, p, RoundingMode::None));
        }
    }
    let mut b_tilde = Vec::new();
    for (bj, &bj_int) in b.iter().zip(b_int.iter()) {
        let bif = ExactNum::from_i32(bj_int, p);
        for k in 0..bj_int {
            let num = bj.add(&ExactNum::from_i32(k, p), p, RoundingMode::None);
            b_tilde.push(num.div(&bif, p, RoundingMode::None));
        }
    }
    if a_tilde.len() > TILDE_MAX || b_tilde.len() > TILDE_MAX {
        return None;
    }
    let n_tilde: i32 = a_int.iter().take(n).sum();
    let m_tilde: i32 = b_int.iter().take(m).sum();
    if n_tilde < 0 || m_tilde < 0 {
        return None;
    }
    let n_t = n_tilde as usize;
    let m_t = m_tilde as usize;
    if n_t > a_tilde.len() || m_t > b_tilde.len() {
        return None;
    }
    let a_star = a_int.iter().take(n).sum::<i32>() - a_int.iter().skip(n).sum::<i32>()
        + b_int.iter().take(m).sum::<i32>()
        - b_int.iter().skip(m).sum::<i32>();
    // c* = m + n - (p+q)/2  (original counts)
    let c_star_num = 2 * ((m + n) as i32) - (pn as i32) - (qn as i32);
    let one = ExactNum::from_u8(1, p);
    let two = ExactNum::from_u8(2, p);
    let half = one.div(&two, p, RoundingMode::None);
    let mut beta = ExactNum::from_u8(1, p);
    for &ai in &a_int {
        let base = ExactNum::from_i32(ai, p);
        beta = beta.div(&pow_i32(&base, ai, p, cc), p, RoundingMode::None);
    }
    for &bj in &b_int {
        let base = ExactNum::from_i32(bj, p);
        beta = beta.mul(&pow_i32(&base, bj, p, cc), p, RoundingMode::None);
    }
    let mut m_fac = ExactNum::from_u8(1, p);
    for (bj, &bj_int) in b.iter().zip(b_int.iter()) {
        let base = ExactNum::from_i32(bj_int, p);
        let expn = bj.sub(&half, p, RoundingMode::None);
        m_fac = m_fac.mul(&pow_real(&base, &expn, p, cc), p, RoundingMode::None);
    }
    for (ai, &ai_int) in a.iter().zip(a_int.iter()) {
        let base = ExactNum::from_i32(ai_int, p);
        let expn = ai.sub(&half, p, RoundingMode::None);
        m_fac = m_fac.div(&pow_real(&base, &expn, p, cc), p, RoundingMode::None);
    }
    let pi = cc.pi(p, RoundingMode::None);
    let two_pi = two.mul(&pi, p, RoundingMode::None);
    // (2π)^{c* - a*/2} with c* = c_star_num/2
    let expn = ExactNum::from_i32(c_star_num, p)
        .div(&two, p, RoundingMode::None)
        .sub(
            &ExactNum::from_i32(a_star, p).div(&two, p, RoundingMode::None),
            p,
            RoundingMode::None,
        );
    let pre = ExactNum::from_i32(d, p)
        .mul(&pow_real(&two_pi, &expn, p, cc), p, RoundingMode::None)
        .mul(&m_fac, p, RoundingMode::None);
    let r = one.div(&ExactNum::from_i32(d, p), p, RoundingMode::None);
    let beta_r = pow_real(&beta, &r, p, cc);
    if !finite(&beta_r) || beta_r.is_zero() {
        return None;
    }
    let z_arg = z.div(&beta_r, p, RoundingMode::None);
    let an_t = &a_tilde[..n_t];
    let ap_t = &a_tilde[n_t..];
    let bm_t = &b_tilde[..m_t];
    let bq_t = &b_tilde[m_t..];
    let g = meijer_g_inner(&z_arg, an_t, ap_t, bm_t, bq_t, &r, p, cc);
    if !finite(&g) || !finite(&pre) {
        return None;
    }
    Some(pre.mul(&g, p, RoundingMode::None))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(a: &ExactNum, b: &ExactNum, p: usize) -> bool {
        let d = a.sub(b, p, RoundingMode::None).abs();
        d.is_zero() || d.exponent().is_some_and(|e| e < -((p as i32) / 4))
    }

    #[test]
    fn test_meijer_heun_identities() {
        let p = 128;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new().unwrap();
        let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, rm);
        let z = half.clone();
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);

        let g = z.meijer_g(&[], &[], &[zero.clone()], &[], p, rm, &mut cc);
        let e = z.neg().exp(p, rm, &mut cc);
        assert!(near(&g, &e, p), "G^{{1,0}}_{{0,1}}(1/2) = e^{{-1/2}}");

        let fox = z.fox_h(
            &[],
            &[],
            &[(zero.clone(), one.clone())],
            &[],
            p,
            rm,
            &mut cc,
        );
        assert!(near(&fox, &e, p), "H with A=B=1 is G");

        let a03 = ExactNum::from_u8(3, p).div(&ExactNum::from_u8(10, p), p, rm);
        let g11 = z.meijer_g(&[a03.clone()], &[], &[zero.clone()], &[], p, rm, &mut cc);
        let want = ExactNum::from_u8(1, p)
            .sub(&a03, p, rm)
            .gamma(p, rm, &mut cc)
            .mul(
                &ExactNum::from_u8(1, p)
                    .add(&z, p, rm)
                    .pow(&a03.sub(&one, p, rm), p, rm, &mut cc),
                p,
                rm,
            );
        assert!(
            near(&g11, &want, p),
            "G^{{1,1}}_{{1,1}}(z|a;0)=Γ(1-a)(1+z)^{{a-1}}"
        );

        // HeunG → 2F1 when q = a α β and δ = α+β-γ+1.
        let alpha = ExactNum::from_u8(2, p).div(&ExactNum::from_u8(5, p), p, rm);
        let beta = ExactNum::from_u8(7, p).div(&ExactNum::from_u8(10, p), p, rm);
        let gamma = ExactNum::from_u8(6, p).div(&ExactNum::from_u8(5, p), p, rm);
        let zz = ExactNum::from_u8(3, p).div(&ExactNum::from_u8(10, p), p, rm);
        let aa = ExactNum::from_u8(5, p).div(&ExactNum::from_u8(2, p), p, rm);
        let delta = alpha.add(&beta, p, rm).sub(&gamma, p, rm).add(&one, p, rm);
        let qv = aa.mul(&alpha, p, rm).mul(&beta, p, rm);
        let hg = zz.heun_g(&aa, &qv, &alpha, &beta, &gamma, &delta, p, rm, &mut cc);
        let f = alpha.hypergeom_2f1(&beta, &gamma, &zz, p, rm, &mut cc);
        assert!(near(&hg, &f, p), "HeunG 2F1 reduction");

        // HeunC → 1F1 when δ=0, ε=-1, q=α.
        let al = ExactNum::from_u8(2, p).div(&ExactNum::from_u8(5, p), p, rm);
        let gam = ExactNum::from_u8(13, p).div(&ExactNum::from_u8(10, p), p, rm);
        let zc = ExactNum::from_u8(7, p).div(&ExactNum::from_u8(20, p), p, rm);
        let hc = zc.heun_c(&al, &gam, &zero, &one.neg(), &al, p, rm, &mut cc);
        let f1 = zc.hypergeom_1f1(&al.neg(), &gam, p, rm, &mut cc);
        assert!(near(&hc, &f1, p), "HeunC 1F1 reduction");

        let h1 = zz.heun_g(
            &ExactNum::from_u8(2, p),
            &zero,
            &zero,
            &ExactNum::from_u8(1, p),
            &half,
            &half,
            p,
            rm,
            &mut cc,
        );
        assert!(near(&h1, &one, p), "HeunG α=q=0 is 1");

        assert!(ExactNum::from_u8(2, p)
            .heun_g(&aa, &qv, &alpha, &beta, &gamma, &delta, p, rm, &mut cc)
            .is_nan());
        assert!(z.meijer_g(&[], &[], &[], &[], p, rm, &mut cc).is_nan());
    }
}
