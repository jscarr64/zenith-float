//! Meijer \(G\), Fox \(H\), and local / confluent Heun on [`ExactNum`].
//!
//! SoftFloat only. Residue \(G\) follows DLMF 16.17; Fox \(H\)
//! with positive rational \(A,B\) reduces through the Gauss multiplication formula
//! (mpmath `foxh`). Local Heun is DLMF 31.3; confluent Heun is the \(z=0\)
//! Frobenius solution of DLMF 31.12.1.

use crate::catalog::{cat_nan, finish, finite, series_cap, small_int, work_p};
use crate::defs::Exponent;
use crate::Consts;
use crate::ExactNum;
use crate::RoundingMode;
use alloc::vec::Vec;
use core::ops::Neg;

/// Max length of any one Meijer / Fox parameter group.
const PARAM_MAX: usize = 8;
/// Max positive integer scale \(A_i,B_j\) after clearing denominators.
const SCALE_MAX: i32 = 12;
/// Max expanded Meijer length after the Fox–Gauss lift.
const TILDE_MAX: usize = 16;

/// \(x\in\{0,-1,-2,\dots\}\). `is_positive` is true for `+0`, so zero is tested first.
fn gamma_pole(x: &ExactNum) -> bool {
    x.is_int() && (x.is_zero() || !x.is_positive())
}

fn all_finite(xs: &[ExactNum]) -> bool {
    xs.iter().all(finite)
}

fn pairs_finite(xs: &[(ExactNum, ExactNum)]) -> bool {
    xs.iter().all(|(a, a_s)| finite(a) && finite(a_s))
}

fn pow_real(z: &ExactNum, e: &ExactNum, p: usize, cc: &mut Consts) -> ExactNum {
    if z.is_zero() {
        // 0^0 = 1 first: `is_positive` is true for +0, so testing it first made 0^(+0) = 0
        // and dropped the b = 0 residue of every G / H at z = 0 (ZF-40).
        if e.is_zero() {
            return ExactNum::from_u8(1, p);
        }
        if e.is_positive() {
            return ExactNum::new(p);
        }
        return cat_nan();
    }
    if z.is_negative() {
        if let Some(n) = small_int(e, p, -64, 64) {
            return z.powsi(n as i64, p, RoundingMode::None);
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
    matches!(z.abs_cmp(&ExactNum::from_u8(1, p)), Some(c) if c > 0)
}

fn choose_series(an_n: usize, ap_n: usize, bm_n: usize, bq_n: usize, z: &ExactNum, p: usize) -> u8 {
    let pn = an_n + ap_n;
    let qn = bm_n + bq_n;
    // Type 1 is {}_p F_{q-1}; type 2 is {}_q F_{p-1}. When p=q both are
    // disk series, so |z| picks the side (we do not continue {}_r F_{r-1}).
    let mut s = if pn > qn || (pn == qn && abs_gt_one(z, p)) { 2u8 } else { 1 };
    if s == 1 && bm_n == 0 {
        s = 2;
    }
    if s == 2 && an_n == 0 {
        s = 1;
    }
    s
}

enum GEval {
    Val(ExactNum),
    /// Removable coincident pole (logarithmic residue still to take).
    Singular,
    Bad,
}

fn meijer_g_series(
    z: &ExactNum,
    an: &[ExactNum],
    ap: &[ExactNum],
    bm: &[ExactNum],
    bq: &[ExactNum],
    r: &ExactNum,
    p: usize,
    cc: &mut Consts,
) -> GEval {
    if an.len() + ap.len() > PARAM_MAX || bm.len() + bq.len() > PARAM_MAX {
        return GEval::Bad;
    }
    if !finite(z)
        || !finite(r)
        || !r.is_positive()
        || !all_finite(an)
        || !all_finite(ap)
        || !all_finite(bm)
        || !all_finite(bq)
    {
        return GEval::Bad;
    }
    if an.is_empty() && bm.is_empty() {
        return GEval::Bad;
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
    let r_is_one = r.compare(&one) == Some(0);
    if !r_is_one && !z.is_positive() && !z.is_zero() {
        return GEval::Bad;
    }
    let z_to_inv_r = if r_is_one {
        z.clone()
    } else {
        pow_real(
            z,
            &ExactNum::from_u8(1, p).div(r, p, RoundingMode::None),
            p,
            cc,
        )
    };
    if !finite(&z_to_inv_r) {
        return GEval::Bad;
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
                    gn.push(ExactNum::sub_ref(bj, bh, p, RoundingMode::None));
                }
            }
            for aj in a.iter().take(n) {
                gn.push(ExactNum::sub_ref(&one, aj, p, RoundingMode::None).add(
                    bh,
                    p,
                    RoundingMode::None,
                ));
            }
            let mut gd = Vec::new();
            for aj in a.iter().skip(n) {
                gd.push(ExactNum::sub_ref(aj, bh, p, RoundingMode::None));
            }
            for bj in b.iter().skip(m) {
                gd.push(ExactNum::sub_ref(&one, bj, p, RoundingMode::None).add(
                    bh,
                    p,
                    RoundingMode::None,
                ));
            }
            let ratio = match gamma_ratio(&gn, &gd, p, cc) {
                Err(()) => return GEval::Singular,
                Ok(None) => continue,
                Ok(Some(v)) => v,
            };
            let mut hn = Vec::with_capacity(pn);
            for aj in &a {
                hn.push(ExactNum::sub_ref(&one, aj, p, RoundingMode::None).add(
                    bh,
                    p,
                    RoundingMode::None,
                ));
            }
            let mut hd = Vec::with_capacity(qn.saturating_sub(1));
            for (j, bj) in b.iter().enumerate() {
                if j != h {
                    hd.push(ExactNum::sub_ref(&one, bj, p, RoundingMode::None).add(
                        bh,
                        p,
                        RoundingMode::None,
                    ));
                }
            }
            let f = hz.hypergeom_pfq(&hn, &hd, p, RoundingMode::None, cc);
            if f.is_nan() {
                return GEval::Singular;
            }
            let expn_div;
            let expn = if r_is_one {
                bh
            } else {
                expn_div = ExactNum::div_ref(bh, r, p, RoundingMode::None);
                &expn_div
            };
            let zp = pow_real(z, expn, p, cc);
            if !finite(&zp) {
                return GEval::Bad;
            }
            if f.is_inf() {
                // ZF-41: the pFq diverges at |z| = 1 (e.g. ln(1+z) at z = −1 is
                // −z·₂F₁(1,1;2;1)): a pole, so this residue term is a signed infinity.
                // Two opposite infinities are not a value; that goes to the limit path.
                if ratio.is_zero() || zp.is_zero() {
                    return GEval::Singular;
                }
                let neg = ratio.is_negative() ^ zp.is_negative() ^ f.is_negative();
                let term = if neg { crate::ext::INF_NEG } else { crate::ext::INF_POS };
                sum = sum.add(&term, p, RoundingMode::None);
                if sum.is_nan() {
                    return GEval::Singular;
                }
                continue;
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
            return GEval::Bad;
        } else {
            ExactNum::from_u8(1, p).div(&z_to_inv_r, p, RoundingMode::None)
        };
        let hz = if sign_odd { inv_z.neg() } else { inv_z };
        for k in 0..n {
            let ak = &a[k];
            let mut gn = Vec::new();
            for (j, aj) in a.iter().enumerate().take(n) {
                if j != k {
                    gn.push(ExactNum::sub_ref(ak, aj, p, RoundingMode::None));
                }
            }
            for bj in b.iter().take(m) {
                gn.push(ExactNum::sub_ref(&one, ak, p, RoundingMode::None).add(
                    bj,
                    p,
                    RoundingMode::None,
                ));
            }
            let mut gd = Vec::new();
            for bj in b.iter().skip(m) {
                gd.push(ExactNum::sub_ref(ak, bj, p, RoundingMode::None));
            }
            for aj in a.iter().skip(n) {
                gd.push(ExactNum::sub_ref(&one, ak, p, RoundingMode::None).add(
                    aj,
                    p,
                    RoundingMode::None,
                ));
            }
            let ratio = match gamma_ratio(&gn, &gd, p, cc) {
                Err(()) => return GEval::Singular,
                Ok(None) => continue,
                Ok(Some(v)) => v,
            };
            let mut hn = Vec::with_capacity(qn);
            for bj in &b {
                hn.push(ExactNum::sub_ref(&one, ak, p, RoundingMode::None).add(
                    bj,
                    p,
                    RoundingMode::None,
                ));
            }
            let mut hd = Vec::with_capacity(pn.saturating_sub(1));
            for (j, aj) in a.iter().enumerate() {
                if j != k {
                    hd.push(ExactNum::from_u8(1, p).add(aj, p, RoundingMode::None).sub(
                        ak,
                        p,
                        RoundingMode::None,
                    ));
                }
            }
            let f = hz.hypergeom_pfq(&hn, &hd, p, RoundingMode::None, cc);
            if !finite(&f) {
                return GEval::Singular;
            }
            let expn = if r_is_one {
                ExactNum::sub_ref(ak, &one, p, RoundingMode::None)
            } else {
                ExactNum::from_u8(1, p)
                    .neg()
                    .add(ak, p, RoundingMode::None)
                    .div(r, p, RoundingMode::None)
            };
            let zp = pow_real(z, &expn, p, cc);
            if !finite(&zp) {
                return GEval::Bad;
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
    GEval::Val(sum)
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
    match meijer_g_series(z, an, ap, bm, bq, r, p, cc) {
        GEval::Val(v) => v,
        GEval::Bad => cat_nan(),
        GEval::Singular => meijer_removable_limit(z, an, ap, bm, bq, r, p, cc),
    }
}

/// Order window, working precision, and the largest intermediate magnitude
/// (`mag`, an exponent) seen by any [`Jet`] operation of one evaluation.
struct JetW {
    lo: i32,
    hi: i32,
    wp: usize,
    mag: core::cell::Cell<i64>,
    ops: core::cell::Cell<u64>,
}

/// Truncated Laurent series in \(\varepsilon\): `c[i]` is the coefficient of
/// \(\varepsilon^{i+\mathrm{lo}}\).
#[derive(Clone)]
struct Jet {
    c: Vec<ExactNum>,
}

impl JetW {
    fn new(lo: i32, hi: i32, wp: usize) -> Self {
        JetW {
            lo,
            hi,
            wp,
            mag: core::cell::Cell::new(i64::MIN / 4),
            ops: core::cell::Cell::new(0),
        }
    }
    fn len(&self) -> usize {
        (self.hi - self.lo + 1) as usize
    }
    fn idx(&self, o: i32) -> Option<usize> {
        if o < self.lo || o > self.hi {
            None
        } else {
            Some((o - self.lo) as usize)
        }
    }
    fn note(&self, x: &ExactNum) {
        self.ops.set(self.ops.get().saturating_add(1));
        if let (false, Some(e)) = (x.is_zero(), x.exponent()) {
            if e > self.mag.get() {
                self.mag.set(e);
            }
        }
    }
    fn mul_n(&self, x: &ExactNum, y: &ExactNum) -> ExactNum {
        let r = ExactNum::mul_ref(x, y, self.wp, RoundingMode::None);
        self.note(&r);
        r
    }
    fn add_n(&self, x: ExactNum, y: &ExactNum) -> ExactNum {
        self.note(y);
        x.add(y, self.wp, RoundingMode::None)
    }
    fn zero(&self) -> Jet {
        Jet {
            c: (0..self.len()).map(|_| ExactNum::new(self.wp)).collect(),
        }
    }
    fn scalar(&self, x: &ExactNum) -> Jet {
        let mut j = self.zero();
        if let Some(i) = self.idx(0) {
            j.c[i] = x.clone();
        }
        j
    }
    /// \(x + c\varepsilon\).
    fn lin(&self, x: &ExactNum, c: i32) -> Jet {
        let mut j = self.scalar(x);
        if let Some(i) = self.idx(1) {
            j.c[i] = ExactNum::from_i32(c, self.wp);
        }
        j
    }
    fn add(&self, a: &Jet, b: &Jet) -> Jet {
        Jet {
            c: a.c
                .iter()
                .zip(&b.c)
                .map(|(x, y)| self.add_n(x.clone(), y))
                .collect(),
        }
    }
    fn scale(&self, a: &Jet, x: &ExactNum) -> Jet {
        Jet {
            c: a.c.iter().map(|y| self.mul_n(y, x)).collect(),
        }
    }
    fn mul(&self, a: &Jet, b: &Jet) -> Jet {
        let mut r = self.zero();
        for (i, x) in a.c.iter().enumerate() {
            if x.is_zero() {
                continue;
            }
            for (j, y) in b.c.iter().enumerate() {
                if y.is_zero() {
                    continue;
                }
                let o = i as i32 + j as i32 + 2 * self.lo;
                if let Some(k) = self.idx(o) {
                    let t = self.mul_n(x, y);
                    r.c[k] = self.add_n(core::mem::take(&mut r.c[k]), &t);
                }
            }
        }
        r
    }
    /// Lowest order with a nonzero coefficient. Every divisor is built from exact
    /// linear factors or Γ values, so a structural zero is an exact zero.
    fn val(&self, a: &Jet) -> Option<i32> {
        a.c.iter()
            .position(|x| !x.is_zero())
            .map(|i| i as i32 + self.lo)
    }
    fn div(&self, a: &Jet, b: &Jet) -> Option<Jet> {
        let v = self.val(b)?;
        let u: Vec<&ExactNum> = b.c.iter().skip((v - self.lo) as usize).collect();
        let n = self.len();
        let mut inv: Vec<ExactNum> = Vec::with_capacity(n);
        let u0 = u[0];
        inv.push(ExactNum::from_u8(1, self.wp).div(u0, self.wp, RoundingMode::None));
        self.note(&inv[0]);
        for k in 1..n {
            let mut s = ExactNum::new(self.wp);
            for i in 1..=k {
                if let Some(ui) = u.get(i) {
                    let t = self.mul_n(ui, &inv[k - i]);
                    s = self.add_n(s, &t);
                }
            }
            let q = s.neg().div(u0, self.wp, RoundingMode::None);
            self.note(&q);
            inv.push(q);
        }
        let mut r = self.zero();
        for (i, x) in a.c.iter().enumerate() {
            if x.is_zero() {
                continue;
            }
            for (k, y) in inv.iter().enumerate() {
                let o = i as i32 + self.lo + k as i32 - v;
                if let Some(t) = self.idx(o) {
                    let pr = self.mul_n(x, y);
                    r.c[t] = self.add_n(core::mem::take(&mut r.c[t]), &pr);
                }
            }
        }
        Some(r)
    }
    /// \(\exp(s)\) for `s` with no terms below order 1.
    fn exp0(&self, s: &Jet) -> Jet {
        let mut r = self.scalar(&ExactNum::from_u8(1, self.wp));
        let mut t = r.clone();
        for n in 1..=self.hi.max(0) {
            let inv_n = ExactNum::from_u8(1, self.wp).div(
                &ExactNum::from_i32(n, self.wp),
                self.wp,
                RoundingMode::None,
            );
            t = self.scale(&self.mul(&t, s), &inv_n);
            r = self.add(&r, &t);
        }
        r
    }
    /// \(\Gamma(x + c\varepsilon)\). At \(x=-n\) this is
    /// \(\Gamma(1+c\varepsilon)/\bigl(c\varepsilon\prod_{i=1}^{n}(c\varepsilon-i)\bigr)\);
    /// elsewhere \(\Gamma(x)\exp\sum_{j\ge1}\psi^{(j-1)}(x)(c\varepsilon)^j/j!\) (DLMF 5.7.3 / 5.15).
    fn gamma(&self, x: &ExactNum, c: i32, cc: &mut Consts) -> Option<Jet> {
        if gamma_pole(x) {
            let n = x.clone().neg().to_u32(RoundingMode::None).ok()?;
            if n > 4096 {
                return None;
            }
            let g1 = self.gamma(&ExactNum::from_u8(1, self.wp), c, cc)?;
            let mut den = self.lin(&ExactNum::new(self.wp), c);
            for i in 1..=n {
                den = self.mul(&den, &self.lin(&ExactNum::from_u32(i, self.wp).neg(), c));
            }
            return self.div(&g1, &den);
        }
        let g = x.gamma(self.wp, RoundingMode::None, cc);
        if !finite(&g) || g.is_zero() {
            return None;
        }
        let mut s = self.zero();
        let mut ck = ExactNum::from_u8(1, self.wp);
        let mut fact = ExactNum::from_u8(1, self.wp);
        for j in 1..=self.hi.max(0) {
            ck = ck.mul(&ExactNum::from_i32(c, self.wp), self.wp, RoundingMode::None);
            fact = fact.mul(&ExactNum::from_i32(j, self.wp), self.wp, RoundingMode::None);
            let psi = x.polygamma((j - 1) as usize, self.wp, RoundingMode::None, cc);
            if !finite(&psi) {
                return None;
            }
            if let Some(i) = self.idx(j) {
                s.c[i] = psi.mul(&ck, self.wp, RoundingMode::None).div(
                    &fact,
                    self.wp,
                    RoundingMode::None,
                );
                self.note(&s.c[i]);
            }
        }
        Some(self.scale(&self.exp0(&s), &g))
    }
}

/// `|x|` rounded up to 64 bits.
fn up64(x: &ExactNum) -> ExactNum {
    let mut y = x.clone().abs();
    let _ = y.set_precision(64, RoundingMode::Up);
    y
}

/// `|x|` rounded down to 64 bits.
fn dn64(x: &ExactNum) -> ExactNum {
    let mut y = x.clone().abs();
    let _ = y.set_precision(64, RoundingMode::Down);
    y
}

fn lt(a: &ExactNum, b: &ExactNum) -> bool {
    matches!(a.compare(b), Some(c) if c < 0)
}

/// Bits that hold every parameter and every difference or `1 - x` of them exactly.
fn exact_span(xs: &[&ExactNum]) -> Option<usize> {
    let mut hi: i64 = 1;
    let mut lo: i64 = 0;
    for x in xs {
        if x.is_zero() {
            continue;
        }
        let e = x.exponent()?;
        let mb = x.mantissa_max_bit_len()? as i64;
        hi = hi.max(e);
        lo = lo.min(e - mb);
    }
    usize::try_from(hi - lo + 4).ok()
}

/// Smallest nonzero \(|b+k|\) over integers \(k\ge0\), as a 64-bit lower bound.
fn min_shift_gap(b: &ExactNum, wp: usize) -> ExactNum {
    let one = ExactNum::from_u8(1, 64);
    if !b.is_negative() {
        return if b.is_zero() { one } else { dn64(b) };
    }
    let f = ExactNum::sub_ref(b, &b.clone().floor(), wp, RoundingMode::None);
    if f.is_zero() {
        return one;
    }
    let g = ExactNum::sub_ref(&ExactNum::from_u8(1, wp), &f, wp, RoundingMode::None);
    if lt(&f, &g) {
        dn64(&f)
    } else {
        dn64(&g)
    }
}

/// \({}_pF_q(x)\) with parameters \(v_i + c_i\varepsilon\), as a [`Jet`].
///
/// Terms are added until a Cauchy majorant on \(|\varepsilon|=\rho\) bounds the
/// tail's contribution to \(\sum_o W_{-o}F_o\) by \(2^{\mathrm{tol}}\), where
/// `weights[o - lo]` \(\ge|W_{-o}|\). \(\rho\le 1/4\) is a power of two with every
/// root \(-(b+k)/c\) of a lower parameter at distance at least \(4\rho\), so the
/// tail is analytic on \(0<|\varepsilon|\le\rho\) and its order-\(o\) coefficient is
/// at most \(\sup_{|\varepsilon|=\rho}|\text{tail}|\,\rho^{-o}\).
fn pfq_jet(
    hn: &[(ExactNum, i32)],
    hd: &[(ExactNum, i32)],
    x: &ExactNum,
    w: &JetW,
    weights: &[ExactNum],
    tol: i64,
) -> Option<Jet> {
    let wp = w.wp;
    let up = RoundingMode::Up;
    let dn = RoundingMode::Down;
    let one64 = ExactNum::from_u8(1, 64);
    let ax = up64(x);
    if hn.len() > hd.len() + 1 || (hn.len() == hd.len() + 1 && !lt(&ax, &one64)) {
        return None;
    }
    let mut rho = ExactNum::from_u8(1, 64).ldexp(-2, 64, dn);
    for (b, c) in hd {
        let ac = ExactNum::from_i32(c.abs().max(1), 64);
        let r1 = min_shift_gap(b, wp).div(&ac, 64, dn).ldexp(-2, 64, dn);
        if lt(&r1, &rho) {
            rho = r1;
        }
    }
    let re = rho.exponent()? - 1;
    let rho = ExactNum::from_u8(1, 64).ldexp(re, 64, RoundingMode::None);
    let crho = |c: i32| ExactNum::from_i32(c.abs(), 64).mul(&rho, 64, up);
    let amax = hn
        .iter()
        .chain(hd.iter())
        .fold(ExactNum::new(64), |m, (v, _)| {
            let a = up64(v);
            if lt(&m, &a) {
                a
            } else {
                m
            }
        });
    let kmin = amax
        .add(&ExactNum::from_u8(2, 64), 64, up)
        .to_u64(RoundingMode::Up)
        .ok()?;
    let mut t = w.scalar(&ExactNum::from_u8(1, wp));
    let mut sum = t.clone();
    let mut mk = ExactNum::from_u8(1, 64);
    let mut k: u64 = 0;
    loop {
        if k > 20_000_000 {
            return None;
        }
        let kk = ExactNum::from_u64(k, wp);
        let mut num = w.scalar(&ExactNum::from_u8(1, wp));
        let mut nb = ax.clone();
        for (a, c) in hn {
            let v = ExactNum::add_ref(a, &kk, wp, RoundingMode::None);
            nb = nb.mul(&up64(&v).add(&crho(*c), 64, up), 64, up);
            num = w.mul(&num, &w.lin(&v, *c));
        }
        let mut den = w.scalar(&ExactNum::from_u64(k + 1, wp));
        let mut db = ExactNum::from_u64(k + 1, 64);
        for (b, c) in hd {
            let v = ExactNum::add_ref(b, &kk, wp, RoundingMode::None);
            let inf = if v.is_zero() {
                ExactNum::from_i32(c.abs(), 64).mul(&rho, 64, dn)
            } else {
                dn64(&v).sub(&crho(*c), 64, dn)
            };
            if !inf.is_positive() || inf.is_zero() {
                return None;
            }
            db = db.mul(&inf, 64, dn);
            den = w.mul(&den, &w.lin(&v, *c));
        }
        mk = mk.mul(&nb, 64, up).div(&db, 64, up);
        t = w.scale(&w.div(&w.mul(&t, &num), &den)?, x);
        sum = w.add(&sum, &t);
        k += 1;
        if k < kmin {
            continue;
        }
        // Ratio bound for every step j >= k: |x| prod (j + max(|a|+|c|rho, 1)) /
        // ((j+1) prod (j - |b| - |c|rho)); each factor is non-increasing in j.
        let kk = ExactNum::from_u64(k, 64);
        let mut bnd = ax.clone();
        for (a, c) in hn {
            let al = up64(a).add(&crho(*c), 64, up);
            let al = if lt(&al, &one64) { one64.clone() } else { al };
            bnd = bnd.mul(&ExactNum::add_ref(&kk, &al, 64, up), 64, up);
        }
        let mut d = ExactNum::from_u64(k + 1, 64);
        for (b, c) in hd {
            let q = ExactNum::sub_ref(&kk, &up64(b), 64, dn).sub(&crho(*c), 64, dn);
            if !q.is_positive() || q.is_zero() {
                return None;
            }
            d = d.mul(&q, 64, dn);
        }
        let bnd = bnd.div(&d, 64, up);
        if !lt(&bnd, &one64) {
            continue;
        }
        let gap = ExactNum::sub_ref(&one64, &bnd, 64, dn);
        let tail = ExactNum::mul_ref(&mk, &bnd, 64, up).div(&gap, 64, up);
        let mut err = ExactNum::new(64);
        for o in w.lo..=w.hi {
            let wt = &weights[(o - w.lo) as usize];
            if wt.is_zero() {
                continue;
            }
            err = err.add(
                &ExactNum::mul_ref(&tail, wt, 64, up).ldexp(-i64::from(o) * re, 64, up),
                64,
                up,
            );
        }
        if err.is_zero() || err.exponent().is_some_and(|e| e < tol) {
            return Some(sum);
        }
    }
}

/// \(G^{m,n}_{p,q}\) at a removable coincident pole of the residue sum, for
/// \(w = z^{1/r} > 0\) (or \(1/w\) when `invert`, with the parameters already
/// reflected by DLMF 16.19.1).
///
/// Every parameter moves to \(x_i + c_i\varepsilon\) with distinct \(c_i\). The
/// perturbed residue sum (DLMF 16.17.2) is analytic in \(\varepsilon\) with value
/// \(G\) at \(\varepsilon=0\), so \(G\) is the constant term of its Laurent
/// expansion: Γ jets from \(\psi^{(j)}\) (the digamma / derivative-of-Gamma form of
/// the logarithmic residue), \(z^{b_h+c\varepsilon}=z^{b_h}e^{c\varepsilon\ln z}\), and
/// \({}_pF_q\) jets with a rigorous tail. No step size and no agreement test.
/// The working precision grows with the largest intermediate magnitude against
/// \(|G|\); a precision past `16p + 2^{15}` refuses with NaN.
fn meijer_jet_limit(
    z: &ExactNum,
    r: &ExactNum,
    invert: bool,
    an: &[ExactNum],
    ap: &[ExactNum],
    bm: &[ExactNum],
    bq: &[ExactNum],
    p: usize,
    cc: &mut Consts,
) -> ExactNum {
    if !z.is_positive() || z.is_zero() || bm.is_empty() {
        return cat_nan();
    }
    let n = an.len();
    let m = bm.len();
    let a: Vec<&ExactNum> = an.iter().chain(ap.iter()).collect();
    let b: Vec<&ExactNum> = bm.iter().chain(bq.iter()).collect();
    let ca: Vec<i32> = (0..a.len()).map(|i| i as i32 + 1).collect();
    let cb: Vec<i32> = (0..b.len()).map(|i| (a.len() + i) as i32 + 1).collect();
    let all: Vec<&ExactNum> = a.iter().chain(b.iter()).copied().collect();
    let Some(span) = exact_span(&all) else {
        return cat_nan();
    };
    let one_r = ExactNum::from_u8(1, p);
    let r_is_one = r.compare(&one_r) == Some(0);
    let mut extra = 4 * crate::defs::WORD_BIT_SIZE;
    loop {
        let wp = p.max(span).saturating_add(extra);
        if wp > 16 * p + 32_768 {
            return cat_nan();
        }
        let one = ExactNum::from_u8(1, wp);
        // w and ln w at wp.
        let lnz = z.ln(wp, RoundingMode::None, cc);
        let lnw0 = if r_is_one { lnz } else { lnz.div(r, wp, RoundingMode::None) };
        let lnw = if invert { lnw0.neg() } else { lnw0 };
        let w = lnw.exp(wp, RoundingMode::None, cc);
        let w = if r_is_one && !invert { z.clone() } else { w };
        if !finite(&w) || !finite(&lnw) {
            return cat_nan();
        }
        let x = if odd_i(a.len() as i32 - m as i32 - n as i32) {
            w.clone().neg()
        } else {
            w.clone()
        };
        let mut g = ExactNum::new(wp);
        let mut mag = i64::MIN / 4;
        let mut ops: u64 = 0;
        let mut tol_used = i64::MIN / 4;
        let d1 = |u: &ExactNum, v: &ExactNum| ExactNum::sub_ref(u, v, wp, RoundingMode::None);
        let d2 = |u: &ExactNum, v: &ExactNum| {
            ExactNum::sub_ref(&one, u, wp, RoundingMode::None).add(v, wp, RoundingMode::None)
        };
        for h in 0..m {
            let bh = b[h];
            let mut gn: Vec<(ExactNum, i32)> = Vec::new();
            for j in 0..m {
                if j != h {
                    gn.push((d1(b[j], bh), cb[j] - cb[h]));
                }
            }
            for j in 0..n {
                gn.push((d2(a[j], bh), cb[h] - ca[j]));
            }
            let mut gd: Vec<(ExactNum, i32)> = Vec::new();
            for j in n..a.len() {
                gd.push((d1(a[j], bh), ca[j] - cb[h]));
            }
            for j in m..b.len() {
                gd.push((d2(b[j], bh), cb[h] - cb[j]));
            }
            let hn: Vec<(ExactNum, i32)> = (0..a.len())
                .map(|j| (d2(a[j], bh), cb[h] - ca[j]))
                .collect();
            let hd: Vec<(ExactNum, i32)> = (0..b.len())
                .filter(|&j| j != h)
                .map(|j| (d2(b[j], bh), cb[h] - cb[j]))
                .collect();
            let poles = gn
                .iter()
                .chain(gd.iter())
                .chain(hd.iter())
                .filter(|(v, _)| gamma_pole(v))
                .count() as i32;
            let hi = 4 * poles + 4;
            if hi > 63 {
                return cat_nan();
            }
            let jw = JetW::new(-poles - 1, hi, wp);
            let mut num = jw.scalar(&one);
            for (v, c) in &gn {
                let Some(gj) = jw.gamma(v, *c, cc) else {
                    return cat_nan();
                };
                num = jw.mul(&num, &gj);
            }
            let mut den = jw.scalar(&one);
            for (v, c) in &gd {
                let Some(gj) = jw.gamma(v, *c, cc) else {
                    return cat_nan();
                };
                den = jw.mul(&den, &gj);
            }
            let Some(ratio) = jw.div(&num, &den) else {
                return cat_nan();
            };
            let wb = pow_real(&w, bh, wp, cc);
            if !finite(&wb) {
                return cat_nan();
            }
            let mut s = jw.zero();
            if let Some(i) = jw.idx(1) {
                s.c[i] = ExactNum::from_i32(cb[h], wp).mul(&lnw, wp, RoundingMode::None);
            }
            let wj = jw.scale(&jw.mul(&ratio, &jw.exp0(&s)), &wb);
            let weights: Vec<ExactNum> = (jw.lo..=jw.hi)
                .map(|o| {
                    jw.idx(-o)
                        .map_or_else(|| ExactNum::new(64), |i| up64(&wj.c[i]))
                })
                .collect();
            let wmax =
                wj.c.iter()
                    .filter(|c| !c.is_zero())
                    .filter_map(|c| c.exponent())
                    .max()
                    .unwrap_or(0);
            let tol = wmax.min(0) - wp as i64 - 8;
            tol_used = tol_used.max(tol);
            let Some(f) = pfq_jet(&hn, &hd, &x, &jw, &weights, tol) else {
                return cat_nan();
            };
            for o in jw.lo..=jw.hi {
                if let (Some(i), Some(k)) = (jw.idx(-o), jw.idx(o)) {
                    let pr = jw.mul_n(&wj.c[i], &f.c[k]);
                    g = jw.add_n(g, &pr);
                }
            }
            mag = mag.max(jw.mag.get());
            ops = ops.saturating_add(jw.ops.get());
        }
        if !finite(&g) {
            return cat_nan();
        }
        if !g.is_zero() {
            let eg = g.exponent().unwrap_or(0);
            // First-order rounding bound: each of `ops` roundings errs by at most
            // 2^{mag - wp}; ask for p + 64 clean bits plus the tail tolerance.
            let ob = 64 - i64::from(ops.leading_zeros());
            let loss = mag.saturating_sub(eg).max(0).saturating_add(ob);
            if loss.saturating_add(p as i64 + 64) <= wp as i64 && tol_used <= eg - p as i64 - 64 {
                return g;
            }
            extra = extra
                .saturating_add(usize::try_from(loss.clamp(64, 1 << 20)).unwrap_or(64))
                .saturating_add(64);
        } else {
            extra = extra.saturating_mul(2);
        }
    }
}

/// \(\lim_{z\to0^+}G\) at a coincident pole, from the leftmost pole \(s^*\) of the
/// Mellin–Barnes integrand (DLMF 16.17.1), of order \(d+1\ge1\):
/// \(G\sim-Cz^{s^*}(\ln z)^d/d!\), \(C\) the product of the factors' leading Laurent
/// coefficients. \(s^*>0\): signed zero; \(s^*<0\) or \(d\ge1\): signed infinity;
/// \(s^*=0,d=0\): \(-C\). No pole at all: \(+0\).
fn meijer_zero_limit(
    an: &[ExactNum],
    ap: &[ExactNum],
    bm: &[ExactNum],
    bq: &[ExactNum],
    p: usize,
    cc: &mut Consts,
) -> ExactNum {
    let all: Vec<&ExactNum> = an
        .iter()
        .chain(ap.iter())
        .chain(bm.iter())
        .chain(bq.iter())
        .collect();
    let Some(span) = exact_span(&all) else {
        return cat_nan();
    };
    let wp = p.max(span).saturating_add(4 * crate::defs::WORD_BIT_SIZE);
    let one = ExactNum::from_u8(1, wp);
    let int_diff =
        |x: &ExactNum, y: &ExactNum| ExactNum::sub_ref(x, y, wp, RoundingMode::None).is_int();
    let le = |x: &ExactNum, y: &ExactNum| matches!(x.compare(y), Some(c) if c <= 0);
    let order = |s: &ExactNum| -> i64 {
        let nb = bm.iter().filter(|b| int_diff(s, b) && le(b, s)).count() as i64;
        let na = ap.iter().filter(|a| int_diff(s, a) && le(a, s)).count() as i64;
        let nq = bq
            .iter()
            .filter(|b| {
                int_diff(s, b) && le(s, &ExactNum::sub_ref(b, &one, wp, RoundingMode::None))
            })
            .count() as i64;
        nb - na - nq
    };
    let mut best: Option<(ExactNum, i64)> = None;
    for b0 in bm {
        let mut cands: Vec<ExactNum> = bm
            .iter()
            .chain(ap.iter())
            .chain(bq.iter())
            .filter(|x| int_diff(x, b0))
            .cloned()
            .collect();
        cands.sort_by(|x, y| {
            x.compare(y)
                .map_or(core::cmp::Ordering::Equal, |c| c.cmp(&0))
        });
        for s in cands {
            let f = order(&s);
            if f >= 1 {
                if best.as_ref().is_none_or(|(bs, _)| lt(&s, bs)) {
                    best = Some((s, f));
                }
                break;
            }
        }
    }
    let Some((s, f)) = best else {
        return ExactNum::new(p);
    };
    let d = f - 1;
    let fact = |k: &ExactNum, cc: &mut Consts| {
        ExactNum::add_ref(k, &one, wp, RoundingMode::None).gamma(wp, RoundingMode::None, cc)
    };
    let even = |k: &ExactNum| k.clone().ldexp(-1, wp, RoundingMode::None).is_int();
    let mut c = ExactNum::from_u8(1, wp);
    for bj in bm {
        // Γ(b - s) near s*: argument -k gives -(-1)^k/k! per unit of (s - s*).
        let v = ExactNum::sub_ref(bj, &s, wp, RoundingMode::None);
        let t = if gamma_pole(&v) {
            let k = v.clone().neg();
            let q = one.clone().div(&fact(&k, cc), wp, RoundingMode::None);
            if even(&k) {
                q.neg()
            } else {
                q
            }
        } else {
            v.gamma(wp, RoundingMode::None, cc)
        };
        c = c.mul(&t, wp, RoundingMode::None);
    }
    for aj in an {
        let v = ExactNum::sub_ref(&one, aj, wp, RoundingMode::None).add(&s, wp, RoundingMode::None);
        if gamma_pole(&v) {
            return cat_nan();
        }
        c = c.mul(&v.gamma(wp, RoundingMode::None, cc), wp, RoundingMode::None);
    }
    for aj in ap {
        // 1/Γ(a - s): argument -k gives -(-1)^k k!.
        let v = ExactNum::sub_ref(aj, &s, wp, RoundingMode::None);
        let t = if gamma_pole(&v) {
            let k = v.clone().neg();
            let q = fact(&k, cc);
            if even(&k) {
                q.neg()
            } else {
                q
            }
        } else {
            one.clone()
                .div(&v.gamma(wp, RoundingMode::None, cc), wp, RoundingMode::None)
        };
        c = c.mul(&t, wp, RoundingMode::None);
    }
    for bj in bq {
        // 1/Γ(1 - b + s): argument -k gives (-1)^k k!.
        let v = ExactNum::sub_ref(&one, bj, wp, RoundingMode::None).add(&s, wp, RoundingMode::None);
        let t = if gamma_pole(&v) {
            let k = v.clone().neg();
            let q = fact(&k, cc);
            if even(&k) {
                q
            } else {
                q.neg()
            }
        } else {
            one.clone()
                .div(&v.gamma(wp, RoundingMode::None, cc), wp, RoundingMode::None)
        };
        c = c.mul(&t, wp, RoundingMode::None);
    }
    if !finite(&c) || c.is_zero() {
        return cat_nan();
    }
    if s.is_zero() && d == 0 {
        return c.neg();
    }
    // sign of -C (ln z)^d with ln z -> -inf
    let neg = c.is_positive() ^ (d % 2 == 1);
    if s.is_positive() && !s.is_zero() {
        let zr = ExactNum::new(p);
        return if neg { zr.neg() } else { zr };
    }
    if neg {
        crate::ext::INF_NEG
    } else {
        crate::ext::INF_POS
    }
}

/// DLMF 16.17: the Mellin–Barnes contour exists only when no \(a_k-b_h\)
/// (\(k\le n\), \(h\le m\)) is a positive integer.
fn nonremovable(an: &[ExactNum], bm: &[ExactNum], wp: usize) -> bool {
    an.iter().any(|a| {
        bm.iter().any(|b| {
            let d = ExactNum::sub_ref(a, b, wp, RoundingMode::None);
            d.is_int() && d.is_positive() && !d.is_zero()
        })
    })
}

/// A coincident pole of the residue sum. Non-removable parameters are
/// `NaN(InvalidArgument)`; \(z=0\) is the \(z\to0^+\) limit; otherwise the analytic
/// constant term ([`meijer_jet_limit`]) of type 1 at \(w=z^{1/r}\), or of type 2 via
/// \(G^{m,n}_{p,q}(w\mid a;b)=G^{n,m}_{q,p}(1/w\mid1-b;1-a)\) (DLMF 16.19.1).
fn meijer_removable_limit(
    z: &ExactNum,
    an: &[ExactNum],
    ap: &[ExactNum],
    bm: &[ExactNum],
    bq: &[ExactNum],
    r: &ExactNum,
    p: usize,
    cc: &mut Consts,
) -> ExactNum {
    let all: Vec<&ExactNum> = an
        .iter()
        .chain(ap.iter())
        .chain(bm.iter())
        .chain(bq.iter())
        .collect();
    let Some(span) = exact_span(&all) else {
        return cat_nan();
    };
    let wp = p.max(span).saturating_add(64);
    if nonremovable(an, bm, wp) {
        return cat_nan();
    }
    if z.is_zero() {
        return meijer_zero_limit(an, ap, bm, bq, p, cc);
    }
    if !z.is_positive() {
        return cat_nan();
    }
    let one = ExactNum::from_u8(1, p);
    let w = if r.compare(&one) == Some(0) {
        z.clone()
    } else {
        pow_real(
            z,
            &ExactNum::from_u8(1, p).div(r, p, RoundingMode::None),
            p,
            cc,
        )
    };
    if !finite(&w) {
        return cat_nan();
    }
    if choose_series(an.len(), ap.len(), bm.len(), bq.len(), &w, p) == 1 {
        return meijer_jet_limit(z, r, false, an, ap, bm, bq, p, cc);
    }
    let refl = |xs: &[ExactNum]| -> Vec<ExactNum> {
        xs.iter()
            .map(|x| ExactNum::sub_ref(&ExactNum::from_u8(1, wp), x, wp, RoundingMode::None))
            .collect()
    };
    meijer_jet_limit(
        z,
        r,
        true,
        &refl(bm),
        &refl(bq),
        &refl(an),
        &refl(ap),
        p,
        cc,
    )
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
        let d = ExactNum::from_i32(n, p)
            .neg()
            .add(x, p, RoundingMode::None)
            .abs();
        if d.is_zero() || d.exponent().is_some_and(|e| e + ((p as i64) / 2) < 0) {
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
        let prod = ExactNum::from_u32(den, p).mul(x, p, RoundingMode::None);
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
    xs.iter().all(|(_, s)| s.compare(&one) == Some(0))
}

fn pow_i32(base: &ExactNum, n: i32, p: usize, _cc: &mut Consts) -> ExactNum {
    if n >= 0 {
        base.powsi(n as i64, p, RoundingMode::None)
    } else {
        let y = base.powsi((-n) as i64, p, RoundingMode::None);
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
    ///   At a coincident pole a negative argument is always `NaN`.
    /// - Coincident poles (integer differences among the active \(b_1,\ldots,b_m\) or
    ///   \(a_1,\ldots,a_n\), and cancelling numerator \(\Gamma\) poles): if some
    ///   \(a_k-b_h\) (\(k\le n\), \(h\le m\)) is a positive integer the Mellin–Barnes
    ///   contour does not exist and the result is `NaN(InvalidArgument)` (DLMF 16.17).
    ///   Otherwise \(z=0\) is the \(z\to0^+\) limit from the leftmost pole (signed
    ///   zero, finite, or signed infinity for a logarithmic or negative-power pole),
    ///   and \(z>0\) is the constant term of the Laurent expansion in a parameter
    ///   shift: \(\Gamma\) jets from \(\psi^{(j)}\), \({}_pF_q\) jets with a rigorous
    ///   tail bound. Pole orders above 14 are `NaN`.
    /// - At most 8 parameters in \(a\) and in \(b\).
    ///
    /// Accuracy: unspecified
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
        let this = self.exact_copy_ext();
        let an_c: Vec<Self> = an.iter().map(|v| v.exact_copy_ext()).collect();
        let ap_c: Vec<Self> = ap.iter().map(|v| v.exact_copy_ext()).collect();
        let bm_c: Vec<Self> = bm.iter().map(|v| v.exact_copy_ext()).collect();
        let bq_c: Vec<Self> = bq.iter().map(|v| v.exact_copy_ext()).collect();
        Self::round_exact(
            self.inexact_accumulated()
                || an.iter().any(|v| v.inexact_accumulated())
                || ap.iter().any(|v| v.inexact_accumulated())
                || bm.iter().any(|v| v.inexact_accumulated())
                || bq.iter().any(|v| v.inexact_accumulated()),
            p,
            rm,
            |q, r| this.meijer_g_word(&an_c, &ap_c, &bm_c, &bq_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::meijer_g`].
    fn meijer_g_word(
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
    /// - Coincident poles after the lift use the same limit as [`meijer_g`].
    /// - Real \(\mathrm{self}>0\) when a scale is not 1 (fractional \(z^{1/r}\)).
    ///
    /// Accuracy: unspecified
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
        let this = self.exact_copy_ext();
        let an_c: Vec<(Self, Self)> = an
            .iter()
            .map(|(u, v)| (u.exact_copy_ext(), v.exact_copy_ext()))
            .collect();
        let ap_c: Vec<(Self, Self)> = ap
            .iter()
            .map(|(u, v)| (u.exact_copy_ext(), v.exact_copy_ext()))
            .collect();
        let bm_c: Vec<(Self, Self)> = bm
            .iter()
            .map(|(u, v)| (u.exact_copy_ext(), v.exact_copy_ext()))
            .collect();
        let bq_c: Vec<(Self, Self)> = bq
            .iter()
            .map(|(u, v)| (u.exact_copy_ext(), v.exact_copy_ext()))
            .collect();
        Self::round_exact(
            self.inexact_accumulated()
                || an
                    .iter()
                    .any(|(u, v)| u.inexact_accumulated() || v.inexact_accumulated())
                || ap
                    .iter()
                    .any(|(u, v)| u.inexact_accumulated() || v.inexact_accumulated())
                || bm
                    .iter()
                    .any(|(u, v)| u.inexact_accumulated() || v.inexact_accumulated())
                || bq
                    .iter()
                    .any(|(u, v)| u.inexact_accumulated() || v.inexact_accumulated()),
            p,
            rm,
            |q, r| this.fox_h_word(&an_c, &ap_c, &bm_c, &bq_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::fox_h`].
    fn fox_h_word(
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
    /// - Algorithm: DLMF 31.3 power series inside \(\lvert z\rvert<\min(1,\lvert a\rvert)\).
    ///   Outside that disk, Taylor steps of the ODE along the real segment from a
    ///   seed in the disk, each step half the distance to \(\{0,1,a\}\). Extra word
    ///   of working precision.
    /// - Bound: working precision `p + 2 WORD_BIT_SIZE`; not Ziv-certified.
    ///
    /// # Limitations
    ///
    /// - The value outside the disk is the real-analytic continuation along the
    ///   component of \(\mathbb{R}\setminus\{1,a\}\) that contains \(0\). A singular
    ///   point on the segment (including the endpoint) is `NaN(InvalidArgument)`.
    /// - When \(q=a\alpha\beta\) and \(\delta=\alpha+\beta-\gamma+1\) hold as exact
    ///   dyadics, or the two sides agree to one ulp of `p`, this is
    ///   \({}_2F_1(\alpha,\beta;\gamma;z)\) and uses that function's domain
    ///   (including real \(z>1\) when the series terminates). A gap of several
    ///   ulps stays on the general Heun path.
    /// - Two consecutive power-series coefficients that are exact zeros make a
    ///   polynomial, returned for every finite \(z\).
    /// - More than 8192 Taylor steps is `NaN(InvalidArgument)`.
    /// - \(a=0\) or \(\gamma\in\{0,-1,-2,\ldots\}\) (recurrence pole) is `NaN`.
    ///
    /// Accuracy: unspecified
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
        let this = self.exact_copy_ext();
        let a_c = a.exact_copy_ext();
        let q_c = q.exact_copy_ext();
        let alpha_c = alpha.exact_copy_ext();
        let beta_c = beta.exact_copy_ext();
        let gamma_c = gamma.exact_copy_ext();
        let delta_c = delta.exact_copy_ext();
        Self::round_exact(
            self.inexact_accumulated()
                || a.inexact_accumulated()
                || q.inexact_accumulated()
                || alpha.inexact_accumulated()
                || beta.inexact_accumulated()
                || gamma.inexact_accumulated()
                || delta.inexact_accumulated(),
            p,
            rm,
            |q, r| this.heun_g_word(&a_c, &q_c, &alpha_c, &beta_c, &gamma_c, &delta_c, q, r, cc),
        )
    }

    /// Word-precision implementation of [`Self::heun_g`].
    fn heun_g_word(
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
        if a.is_zero() || gamma_pole(gamma) {
            return cat_nan();
        }
        if self.is_zero() {
            return finish(ExactNum::from_u8(1, p), p, rm);
        }
        if heun_g_is_2f1(a, q, alpha, beta, gamma, delta, p, rm) {
            return alpha.hypergeom_2f1(beta, gamma, self, p, rm, cc);
        }
        let wrk = work_p(p);
        let one = ExactNum::from_u8(1, wrk);
        let negated;
        let rad: &ExactNum = if matches!(a.abs_cmp(&one), Some(c) if c < 0) {
            if a.is_negative() {
                negated = a.neg_via_mul();
                &negated
            } else {
                a
            }
        } else {
            &one
        };
        if matches!(self.abs_cmp(rad), Some(c) if c >= 0) {
            return heun_g_outside(self, a, q, alpha, beta, gamma, delta, p, rm);
        }
        let eps = ExactNum::add_ref(alpha, beta, wrk, RoundingMode::None)
            .sub(gamma, wrk, RoundingMode::None)
            .sub(delta, wrk, RoundingMode::None)
            .add(&one, wrk, RoundingMode::None);
        let mut c_prev = ExactNum::from_u8(1, wrk);
        let mut c = ExactNum::div_ref(
            q,
            &ExactNum::mul_ref(a, gamma, wrk, RoundingMode::None),
            wrk,
            RoundingMode::None,
        );
        let mut zp = self.clone();
        let mut sum = ExactNum::from_u8(1, wrk).add(
            &ExactNum::mul_ref(&c, &zp, wrk, RoundingMode::None),
            wrk,
            RoundingMode::None,
        );
        for j in 1..=series_cap(wrk) {
            let jf = ExactNum::from_u32(j as u32, wrk);
            let jm1 = ExactNum::from_u32((j - 1) as u32, wrk);
            let jp1 = ExactNum::from_u32((j + 1) as u32, wrk);
            let pj = ExactNum::from_u32((j - 1) as u32, wrk)
                .add(alpha, wrk, RoundingMode::None)
                .mul(
                    &ExactNum::from_u32((j - 1) as u32, wrk).add(beta, wrk, RoundingMode::None),
                    wrk,
                    RoundingMode::None,
                );
            // Q_j = j((j-1+γ)(1+a) + aδ + ε)
            let qj = ExactNum::from_u32(j as u32, wrk).mul(
                &jm1.add(gamma, wrk, RoundingMode::None)
                    .mul(
                        &ExactNum::from_u8(1, wrk).add(a, wrk, RoundingMode::None),
                        wrk,
                        RoundingMode::None,
                    )
                    .add(
                        &ExactNum::mul_ref(a, delta, wrk, RoundingMode::None),
                        wrk,
                        RoundingMode::None,
                    )
                    .add(&eps, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            );
            let rj = jp1.mul(a, wrk, RoundingMode::None).mul(
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
            let term = ExactNum::mul_ref(&zp, &c_next, wrk, RoundingMode::None);
            sum = sum.add(&term, wrk, RoundingMode::None);
            if term.is_zero()
                || term
                    .abs()
                    .exponent()
                    .is_some_and(|e| e + (p as i64) < sum.exponent().unwrap_or(0))
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
    /// - Algorithm: Frobenius series of DLMF 31.12.1 about \(z=0\) for \(\lvert z\rvert<1\).
    ///   For real \(z<1\) outside that disk, Taylor steps along \((-\infty,1)\), each
    ///   step half the distance to \(\{0,1\}\). Extra word of working precision.
    /// - Bound: working precision `p + 2 WORD_BIT_SIZE`; not Ziv-certified.
    ///
    /// # Limitations
    ///
    /// - \(z\ge 1\) meets the branch point at \(1\) and is `NaN(InvalidArgument)`
    ///   unless two consecutive series coefficients are exact zeros (a polynomial,
    ///   returned for every finite \(z\)) or the \({}_1F_1\) reduction applies.
    /// - More than 8192 Taylor steps is `NaN(InvalidArgument)`.
    /// - \(\gamma\in\{0,-1,-2,\ldots\}\) is `NaN`.
    /// - When \(\delta=0\), \(\varepsilon=-1\), \(q=\alpha\) hold as exact dyadics,
    ///   or agree to one ulp of `p`, this is \({}_1F_1(-\alpha;\gamma;z)\),
    ///   entire, including \(\lvert z\rvert\ge 1\). A gap of several ulps stays
    ///   on the general confluent path.
    /// - Biconfluent / double-confluent / triconfluent Heun (DLMF 31.12.2–4) are
    ///   not implemented: there is no public SoftFloat signature for them,
    ///   and the \(z=0\) series is not unique (triconfluent) or not regular
    ///   (double-confluent).
    ///
    /// Accuracy: unspecified
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
        if !finite(self)
            || !finite(alpha)
            || !finite(gamma)
            || !finite(delta)
            || !finite(eps)
            || !finite(q)
        {
            return cat_nan();
        }
        if gamma_pole(gamma) {
            return cat_nan();
        }
        if self.is_zero() {
            return finish(ExactNum::from_u8(1, p), p, rm);
        }
        if heun_c_is_1f1(alpha, delta, eps, q, p) {
            return self.hypergeom_1f1(&alpha.neg_via_mul(), gamma, p, rm, cc);
        }
        let wrk = work_p(p);
        let one = ExactNum::from_u8(1, wrk);
        if matches!(self.abs_cmp(&one), Some(c) if c >= 0) {
            return heun_c_outside(self, alpha, gamma, delta, eps, q, p, rm);
        }
        let mut b_prev = ExactNum::from_u8(1, wrk);
        let mut b = q.neg_via_mul().div(gamma, wrk, RoundingMode::None);
        let mut zp = self.clone();
        let mut sum = ExactNum::from_u8(1, wrk).add(
            &ExactNum::mul_ref(&b, &zp, wrk, RoundingMode::None),
            wrk,
            RoundingMode::None,
        );
        for n in 1..=series_cap(wrk) {
            let nf = ExactNum::from_u32(n as u32, wrk);
            let nm1 = ExactNum::from_u32((n - 1) as u32, wrk);
            let np1 = ExactNum::from_u32((n + 1) as u32, wrk);
            let left = np1.mul(
                &ExactNum::from_u32(n as u32, wrk).add(gamma, wrk, RoundingMode::None),
                wrk,
                RoundingMode::None,
            );
            if left.is_zero() || !finite(&left) {
                return cat_nan();
            }
            // n(n-1+γ+δ-ε)-q
            let coef_n = nf
                .mul(
                    &ExactNum::from_u32((n - 1) as u32, wrk)
                        .add(gamma, wrk, RoundingMode::None)
                        .add(delta, wrk, RoundingMode::None)
                        .sub(eps, wrk, RoundingMode::None),
                    wrk,
                    RoundingMode::None,
                )
                .sub(q, wrk, RoundingMode::None);
            let coef_nm1 =
                nm1.mul(eps, wrk, RoundingMode::None)
                    .add(alpha, wrk, RoundingMode::None);
            let b_next = coef_n
                .mul(&b, wrk, RoundingMode::None)
                .add(
                    &coef_nm1.mul(&b_prev, wrk, RoundingMode::None),
                    wrk,
                    RoundingMode::None,
                )
                .div(&left, wrk, RoundingMode::None);
            zp = zp.mul(self, wrk, RoundingMode::None);
            let term = ExactNum::mul_ref(&zp, &b_next, wrk, RoundingMode::None);
            sum = sum.add(&term, wrk, RoundingMode::None);
            if term.is_zero()
                || term
                    .abs()
                    .exponent()
                    .is_some_and(|e| e + (p as i64) < sum.exponent().unwrap_or(0))
            {
                break;
            }
            b_prev = b;
            b = b_next;
        }
        finish(sum, p, rm)
    }
}

/// Taylor steps along one real segment. Each step is at most half the distance
/// to the nearest point of `sings` (the local radius of the ODE).
const HEUN_STEP_MAX: usize = 8192;

struct OdePoly {
    a0: ExactNum,
    a1: ExactNum,
    a2: ExactNum,
    a3: ExactNum,
    b0: ExactNum,
    b1: ExactNum,
    b2: ExactNum,
    c0: ExactNum,
    c1: ExactNum,
}

fn un(n: u32, p: usize) -> ExactNum {
    ExactNum::from_u32(n, p)
}

fn rel_below(term: &ExactNum, peak: Exponent, bits: usize) -> bool {
    if !finite(term) {
        return false;
    }
    if term.is_zero() {
        return true;
    }
    match term.exponent() {
        Some(e) => e + (bits as i64) < peak,
        None => false,
    }
}

fn bump_peak(peak: &mut Exponent, v: &ExactNum) {
    if let Some(e) = v.exponent() {
        if e > *peak {
            *peak = e;
        }
    }
}

fn segment_hits(z: &ExactNum, s: &ExactNum) -> bool {
    if !finite(z) || !finite(s) || s.is_zero() || z.is_zero() {
        return false;
    }
    if z.is_negative() != s.is_negative() {
        return false;
    }
    matches!(s.abs_cmp(z), Some(c) if c <= 0)
}

fn min_dist(z: &ExactNum, sings: &[ExactNum], p: usize) -> ExactNum {
    let mut m = ExactNum::sub_ref(z, &sings[0], p, RoundingMode::None).abs();
    for s in &sings[1..] {
        let d = ExactNum::sub_ref(z, s, p, RoundingMode::None).abs();
        if matches!(d.compare(&m), Some(c) if c < 0) {
            m = d;
        }
    }
    m
}

fn horner(coefs: &[ExactNum], z: &ExactNum, p: usize) -> ExactNum {
    let mut acc = ExactNum::new(p);
    for c in coefs.iter().rev() {
        acc = acc
            .mul(z, p, RoundingMode::None)
            .add(c, p, RoundingMode::None);
    }
    acc
}

fn heun_eps(
    alpha: &ExactNum,
    beta: &ExactNum,
    gamma: &ExactNum,
    delta: &ExactNum,
    p: usize,
) -> ExactNum {
    ExactNum::add_ref(alpha, beta, p, RoundingMode::None)
        .sub(gamma, p, RoundingMode::None)
        .sub(delta, p, RoundingMode::None)
        .add(&ExactNum::from_u8(1, p), p, RoundingMode::None)
}

fn params_exact(xs: &[&ExactNum]) -> bool {
    xs.iter().all(|x| finite(x) && !x.inexact_accumulated())
}

fn exactly_equal(a: &ExactNum, b: &ExactNum) -> bool {
    finite(a) && finite(b) && matches!(a.compare(b), Some(0))
}

/// Agreement to one ulp of `p`, measured at precision `2p`.
///
/// A looser test (about eight bits of slack) once treated a
/// near-miss as the exact Kummer / \({}_2F_1\) reduction.
fn agrees_tight(a: &ExactNum, b: &ExactNum, p: usize) -> bool {
    if !finite(a) || !finite(b) {
        return false;
    }
    let tight = p.saturating_mul(2).max(64);
    let d = ExactNum::sub_ref(a, b, tight, RoundingMode::None).abs();
    if d.is_zero() {
        return true;
    }
    // |a| and |b| have the exponent and the zero-ness of a and b.
    let scale = if matches!(a.abs_cmp(b), Some(c) if c > 0) { a } else { b };
    let Some(de) = d.exponent() else {
        return false;
    };
    let mag = if scale.is_zero() { 0 } else { scale.exponent().unwrap_or(0) };
    de + (p as i64) < mag + 1
}

fn heun_g_is_2f1(
    a: &ExactNum,
    q: &ExactNum,
    alpha: &ExactNum,
    beta: &ExactNum,
    gamma: &ExactNum,
    delta: &ExactNum,
    p: usize,
    _rm: RoundingMode,
) -> bool {
    let one = ExactNum::from_u8(1, p);
    if params_exact(&[a, q, alpha, beta, gamma, delta]) {
        let q_star = a.mul_full_prec(alpha).mul_full_prec(beta);
        let d_star = alpha
            .add_full_prec(beta)
            .sub_full_prec(gamma)
            .add_full_prec(&one);
        return exactly_equal(q, &q_star) && exactly_equal(delta, &d_star);
    }
    let q_star =
        ExactNum::mul_ref(a, alpha, p, RoundingMode::None).mul(beta, p, RoundingMode::None);
    let d_star = ExactNum::add_ref(alpha, beta, p, RoundingMode::None)
        .sub(gamma, p, RoundingMode::None)
        .add(&one, p, RoundingMode::None);
    agrees_tight(q, &q_star, p) && agrees_tight(delta, &d_star, p)
}

fn heun_c_is_1f1(
    alpha: &ExactNum,
    delta: &ExactNum,
    eps: &ExactNum,
    q: &ExactNum,
    p: usize,
) -> bool {
    let zero = ExactNum::new(p);
    let neg_one = ExactNum::from_u8(1, p).neg();
    if params_exact(&[alpha, delta, eps, q]) {
        return delta.is_zero() && exactly_equal(eps, &neg_one) && exactly_equal(q, alpha);
    }
    agrees_tight(delta, &zero, p) && agrees_tight(eps, &neg_one, p) && agrees_tight(q, alpha, p)
}

/// Power-series coefficients of local Heun about 0. `true` when two consecutive
/// coefficients are exact zeros (a polynomial; later coefficients stay zero).
fn heun_g_coefs(
    a: &ExactNum,
    q: &ExactNum,
    alpha: &ExactNum,
    beta: &ExactNum,
    gamma: &ExactNum,
    delta: &ExactNum,
    p: usize,
) -> Option<(Vec<ExactNum>, bool)> {
    let one = ExactNum::from_u8(1, p);
    let eps = heun_eps(alpha, beta, gamma, delta, p);
    let c1 = ExactNum::div_ref(
        q,
        &ExactNum::mul_ref(a, gamma, p, RoundingMode::None),
        p,
        RoundingMode::None,
    );
    if !finite(&c1) {
        return None;
    }
    let mut coefs = Vec::new();
    coefs.push(one.clone());
    coefs.push(c1);
    let mut c_prev = one.clone();
    let mut c = coefs[1].clone();
    let mut poly = false;
    for j in 1..=series_cap(p) {
        let jf = un(j as u32, p);
        let jm1 = un((j - 1) as u32, p);
        let jp1 = un((j + 1) as u32, p);
        let pj = un((j - 1) as u32, p).add(alpha, p, RoundingMode::None).mul(
            &un((j - 1) as u32, p).add(beta, p, RoundingMode::None),
            p,
            RoundingMode::None,
        );
        let qj = un(j as u32, p).mul(
            &jm1.add(gamma, p, RoundingMode::None)
                .mul(
                    &ExactNum::from_u8(1, p).add(a, p, RoundingMode::None),
                    p,
                    RoundingMode::None,
                )
                .add(
                    &ExactNum::mul_ref(a, delta, p, RoundingMode::None),
                    p,
                    RoundingMode::None,
                )
                .add(&eps, p, RoundingMode::None),
            p,
            RoundingMode::None,
        );
        let rj = jp1.mul(a, p, RoundingMode::None).mul(
            &jf.add(gamma, p, RoundingMode::None),
            p,
            RoundingMode::None,
        );
        if rj.is_zero() || !finite(&rj) {
            return None;
        }
        let c_next = qj
            .add(q, p, RoundingMode::None)
            .mul(&c, p, RoundingMode::None)
            .sub(
                &pj.mul(&c_prev, p, RoundingMode::None),
                p,
                RoundingMode::None,
            )
            .div(&rj, p, RoundingMode::None);
        if !finite(&c_next) {
            return None;
        }
        if c.is_zero() && c_next.is_zero() {
            poly = true;
            break;
        }
        coefs.push(c_next.clone());
        c_prev = c;
        c = c_next;
    }
    Some((coefs, poly))
}

fn heun_c_coefs(
    alpha: &ExactNum,
    gamma: &ExactNum,
    delta: &ExactNum,
    eps: &ExactNum,
    q: &ExactNum,
    p: usize,
) -> Option<(Vec<ExactNum>, bool)> {
    let b1 = q.neg_via_mul().div(gamma, p, RoundingMode::None);
    if !finite(&b1) {
        return None;
    }
    let mut coefs = Vec::new();
    coefs.push(ExactNum::from_u8(1, p));
    coefs.push(b1);
    let mut b_prev = coefs[0].clone();
    let mut b = coefs[1].clone();
    let mut poly = false;
    for n in 1..=series_cap(p) {
        let nf = un(n as u32, p);
        let nm1 = un((n - 1) as u32, p);
        let np1 = un((n + 1) as u32, p);
        let left = np1.mul(
            &un(n as u32, p).add(gamma, p, RoundingMode::None),
            p,
            RoundingMode::None,
        );
        if left.is_zero() || !finite(&left) {
            return None;
        }
        let coef_n = nf
            .mul(
                &un((n - 1) as u32, p)
                    .add(gamma, p, RoundingMode::None)
                    .add(delta, p, RoundingMode::None)
                    .sub(eps, p, RoundingMode::None),
                p,
                RoundingMode::None,
            )
            .sub(q, p, RoundingMode::None);
        let coef_nm1 = nm1
            .mul(eps, p, RoundingMode::None)
            .add(alpha, p, RoundingMode::None);
        let b_next = coef_n
            .mul(&b, p, RoundingMode::None)
            .add(
                &coef_nm1.mul(&b_prev, p, RoundingMode::None),
                p,
                RoundingMode::None,
            )
            .div(&left, p, RoundingMode::None);
        if !finite(&b_next) {
            return None;
        }
        if b.is_zero() && b_next.is_zero() {
            poly = true;
            break;
        }
        coefs.push(b_next.clone());
        b_prev = b;
        b = b_next;
    }
    Some((coefs, poly))
}

fn eval_jet(
    coefs: &[ExactNum],
    z: &ExactNum,
    p: usize,
    bits: usize,
) -> Option<(ExactNum, ExactNum)> {
    if coefs.is_empty() || z.is_zero() || !finite(z) {
        return None;
    }
    let mut zp = ExactNum::from_u8(1, p);
    let mut y = ExactNum::new(p);
    let mut yp = ExactNum::new(p);
    let mut peak_y = Exponent::MIN;
    let mut peak_d = Exponent::MIN;
    let mut term_ok = false;
    let mut deriv_ok = false;
    for (n, c) in coefs.iter().enumerate() {
        let term = ExactNum::mul_ref(c, &zp, p, RoundingMode::None);
        if !finite(&term) {
            return None;
        }
        y = y.add(&term, p, RoundingMode::None);
        bump_peak(&mut peak_y, &term);
        bump_peak(&mut peak_y, &y);
        term_ok = rel_below(&term, peak_y, bits);
        if n == 0 {
            deriv_ok = true;
        } else {
            let dterm =
                un(n as u32, p).mul(&term.div(z, p, RoundingMode::None), p, RoundingMode::None);
            if !finite(&dterm) {
                return None;
            }
            yp = yp.add(&dterm, p, RoundingMode::None);
            bump_peak(&mut peak_d, &dterm);
            bump_peak(&mut peak_d, &yp);
            deriv_ok = rel_below(&dterm, peak_d, bits);
        }
        zp = zp.mul(z, p, RoundingMode::None);
    }
    if term_ok && deriv_ok && finite(&y) && finite(&yp) {
        Some((y, yp))
    } else {
        None
    }
}

fn heun_g_ode(
    z0: &ExactNum,
    a: &ExactNum,
    q: &ExactNum,
    alpha: &ExactNum,
    beta: &ExactNum,
    gamma: &ExactNum,
    delta: &ExactNum,
    p: usize,
) -> Option<OdePoly> {
    let one = ExactNum::from_u8(1, p);
    let eps = heun_eps(alpha, beta, gamma, delta, p);
    let zm1 = ExactNum::from_u8(1, p).neg().add(z0, p, RoundingMode::None);
    let zma = ExactNum::sub_ref(z0, a, p, RoundingMode::None);
    let a0 = ExactNum::mul_ref(z0, &zm1, p, RoundingMode::None).mul(&zma, p, RoundingMode::None);
    if a0.is_zero() || !finite(&a0) {
        return None;
    }
    let three = ExactNum::from_u8(3, p);
    let a2 = three
        .mul(z0, p, RoundingMode::None)
        .sub(&one, p, RoundingMode::None)
        .sub(a, p, RoundingMode::None);
    let a1 = ExactNum::mul_ref(z0, &zm1, p, RoundingMode::None)
        .add(
            &ExactNum::mul_ref(z0, &zma, p, RoundingMode::None),
            p,
            RoundingMode::None,
        )
        .add(
            &ExactNum::mul_ref(&zm1, &zma, p, RoundingMode::None),
            p,
            RoundingMode::None,
        );
    let b0 = ExactNum::mul_ref(&zm1, &zma, p, RoundingMode::None)
        .mul(gamma, p, RoundingMode::None)
        .add(
            &ExactNum::mul_ref(z0, &zma, p, RoundingMode::None).mul(delta, p, RoundingMode::None),
            p,
            RoundingMode::None,
        )
        .add(
            &ExactNum::mul_ref(z0, &zm1, p, RoundingMode::None).mul(&eps, p, RoundingMode::None),
            p,
            RoundingMode::None,
        );
    // (z0 - 1) + (z0 - a), with z0 - 1 formed afresh
    let b1 = ExactNum::from_u8(1, p)
        .neg()
        .add(z0, p, RoundingMode::None)
        .add(&zma, p, RoundingMode::None)
        .mul(gamma, p, RoundingMode::None)
        .add(
            &zma.add(z0, p, RoundingMode::None)
                .mul(delta, p, RoundingMode::None),
            p,
            RoundingMode::None,
        )
        .add(
            &zm1.add(z0, p, RoundingMode::None)
                .mul(&eps, p, RoundingMode::None),
            p,
            RoundingMode::None,
        );
    let b2 =
        ExactNum::add_ref(gamma, delta, p, RoundingMode::None).add(&eps, p, RoundingMode::None);
    let ab = ExactNum::mul_ref(alpha, beta, p, RoundingMode::None);
    Some(OdePoly {
        a0,
        a1,
        a2,
        a3: one,
        b0,
        b1,
        b2,
        c0: ExactNum::mul_ref(&ab, z0, p, RoundingMode::None).sub(q, p, RoundingMode::None),
        c1: ab,
    })
}

fn heun_c_ode(
    z0: &ExactNum,
    alpha: &ExactNum,
    gamma: &ExactNum,
    delta: &ExactNum,
    eps: &ExactNum,
    q: &ExactNum,
    p: usize,
) -> Option<OdePoly> {
    let one = ExactNum::from_u8(1, p);
    let zm1 = ExactNum::from_u8(1, p).neg().add(z0, p, RoundingMode::None);
    let a0 = ExactNum::mul_ref(z0, &zm1, p, RoundingMode::None);
    if a0.is_zero() || !finite(&a0) {
        return None;
    }
    let two = ExactNum::from_u8(2, p);
    let a1 = two
        .mul(z0, p, RoundingMode::None)
        .sub(&one, p, RoundingMode::None);
    let b0 = ExactNum::mul_ref(gamma, &zm1, p, RoundingMode::None)
        .add(
            &ExactNum::mul_ref(delta, z0, p, RoundingMode::None),
            p,
            RoundingMode::None,
        )
        .add(
            &zm1.mul(z0, p, RoundingMode::None)
                .mul(eps, p, RoundingMode::None),
            p,
            RoundingMode::None,
        );
    let b1 = ExactNum::add_ref(gamma, delta, p, RoundingMode::None).add(
        &ExactNum::mul_ref(eps, &a1, p, RoundingMode::None),
        p,
        RoundingMode::None,
    );
    Some(OdePoly {
        a0,
        a1,
        a2: one,
        a3: ExactNum::new(p),
        b0,
        b1,
        b2: eps.clone(),
        c0: ExactNum::mul_ref(alpha, z0, p, RoundingMode::None).sub(q, p, RoundingMode::None),
        c1: alpha.clone(),
    })
}

fn taylor_step(
    ode: &OdePoly,
    y: &ExactNum,
    yp: &ExactNum,
    h: &ExactNum,
    p: usize,
    bits: usize,
) -> Option<(ExactNum, ExactNum)> {
    if h.is_zero() {
        return Some((y.clone(), yp.clone()));
    }
    let mut c_km1 = ExactNum::new(p);
    let mut c_k = y.clone();
    let mut c_kp1 = yp.clone();
    let mut hpow = h.clone();
    let mut sum = ExactNum::mul_ref(yp, h, p, RoundingMode::None).add(y, p, RoundingMode::None);
    let mut dsum = yp.clone();
    let mut peak_y = Exponent::MIN;
    let mut peak_d = Exponent::MIN;
    bump_peak(&mut peak_y, &sum);
    bump_peak(&mut peak_d, &dsum);
    let cap = series_cap(p);
    for k in 0..cap {
        let ku = k as u32;
        let kp1 = ku + 1;
        let kp2 = ku + 2;
        let left =
            un(kp1, p)
                .mul(&un(kp2, p), p, RoundingMode::None)
                .mul(&ode.a0, p, RoundingMode::None);
        if left.is_zero() || !finite(&left) {
            return None;
        }
        let coef_kp1 = un(kp1, p).mul(
            &un(ku, p)
                .mul(&ode.a1, p, RoundingMode::None)
                .add(&ode.b0, p, RoundingMode::None),
            p,
            RoundingMode::None,
        );
        let (coef_k, coef_km1) = if k == 0 {
            (ode.c0.clone(), ExactNum::new(p))
        } else {
            let km1 = ku - 1;
            let km2 = ku.saturating_sub(2);
            let coef_k = un(ku, p)
                .mul(&un(km1, p), p, RoundingMode::None)
                .mul(&ode.a2, p, RoundingMode::None)
                .add(
                    &un(ku, p).mul(&ode.b1, p, RoundingMode::None),
                    p,
                    RoundingMode::None,
                )
                .add(&ode.c0, p, RoundingMode::None);
            let coef_km1 = un(km1, p)
                .mul(&un(km2, p), p, RoundingMode::None)
                .mul(&ode.a3, p, RoundingMode::None)
                .add(
                    &un(km1, p).mul(&ode.b2, p, RoundingMode::None),
                    p,
                    RoundingMode::None,
                )
                .add(&ode.c1, p, RoundingMode::None);
            (coef_k, coef_km1)
        };
        let rhs = coef_kp1
            .mul(&c_kp1, p, RoundingMode::None)
            .add(
                &coef_k.mul(&c_k, p, RoundingMode::None),
                p,
                RoundingMode::None,
            )
            .add(
                &coef_km1.mul(&c_km1, p, RoundingMode::None),
                p,
                RoundingMode::None,
            );
        let c_kp2 = rhs.neg().div(&left, p, RoundingMode::None);
        if !finite(&c_kp2) {
            return None;
        }
        let h_prev = hpow.clone();
        hpow = hpow.mul(h, p, RoundingMode::None);
        let term = ExactNum::mul_ref(&hpow, &c_kp2, p, RoundingMode::None);
        let dterm = un(kp2, p).mul(
            &h_prev.mul(&c_kp2, p, RoundingMode::None),
            p,
            RoundingMode::None,
        );
        if !finite(&term) || !finite(&dterm) {
            return None;
        }
        sum = sum.add(&term, p, RoundingMode::None);
        dsum = dsum.add(&dterm, p, RoundingMode::None);
        bump_peak(&mut peak_y, &term);
        bump_peak(&mut peak_y, &sum);
        bump_peak(&mut peak_d, &dterm);
        bump_peak(&mut peak_d, &dsum);
        c_km1 = c_k;
        c_k = c_kp1;
        c_kp1 = c_kp2;
        if rel_below(&term, peak_y, bits) && rel_below(&dterm, peak_d, bits) {
            return Some((sum, dsum));
        }
    }
    None
}

fn close_enough(rem: &ExactNum, target: &ExactNum, bits: usize) -> bool {
    if rem.is_zero() {
        return true;
    }
    let mut peak: Exponent = 0;
    bump_peak(&mut peak, target);
    bump_peak(&mut peak, &ExactNum::from_u8(1, bits.max(64)));
    rel_below(rem, peak, bits.saturating_add(16))
}

fn taylor_walk<F>(
    mut zc: ExactNum,
    mut y: ExactNum,
    mut yp: ExactNum,
    target: &ExactNum,
    sings: &[ExactNum],
    p: usize,
    bits: usize,
    mut ode_at: F,
) -> Option<ExactNum>
where
    F: FnMut(&ExactNum) -> Option<OdePoly>,
{
    let acc = bits.saturating_add(32);
    for _ in 0..HEUN_STEP_MAX {
        let rem = ExactNum::sub_ref(target, &zc, p, RoundingMode::None);
        if close_enough(&rem, target, bits) {
            return if finite(&y) { Some(y) } else { None };
        }
        let dist = min_dist(&zc, sings, p);
        if dist.is_zero() || !finite(&dist) {
            return None;
        }
        let step = dist.ldexp(-1, p, RoundingMode::None);
        let h = if matches!(rem.abs_cmp(&step), Some(c) if c <= 0) {
            rem
        } else if rem.is_negative() {
            step.neg()
        } else {
            step
        };
        if h.is_zero() {
            return if finite(&y) { Some(y) } else { None };
        }
        let ode = ode_at(&zc)?;
        let (yn, ypn) = taylor_step(&ode, &y, &yp, &h, p, acc)?;
        let next = h.add(&zc, p, RoundingMode::None);
        if next.compare(&zc) == Some(0) {
            return if finite(&yn) { Some(yn) } else { None };
        }
        zc = next;
        y = yn;
        yp = ypn;
    }
    None
}

fn disk_seed(z: &ExactNum, radius: ExactNum, p: usize) -> ExactNum {
    let half = radius.ldexp(-1, p, RoundingMode::None);
    if z.is_negative() {
        half.neg()
    } else {
        half
    }
}

fn heun_g_outside(
    z: &ExactNum,
    a: &ExactNum,
    q: &ExactNum,
    alpha: &ExactNum,
    beta: &ExactNum,
    gamma: &ExactNum,
    delta: &ExactNum,
    p: usize,
    rm: RoundingMode,
) -> ExactNum {
    let wrk = work_p(p);
    let Some((coefs, poly)) = heun_g_coefs(a, q, alpha, beta, gamma, delta, wrk) else {
        return cat_nan();
    };
    if poly {
        return finish(horner(&coefs, z, wrk), p, rm);
    }
    let one = ExactNum::from_u8(1, wrk);
    if segment_hits(z, &one) || segment_hits(z, a) {
        return cat_nan();
    }
    let seed = if matches!(a.abs_cmp(&one), Some(c) if c < 0) {
        // ±|a|/2 with the sign of z: scale -a (fresh, owned) and flip the sign when it differs.
        let mut half = a.neg_via_mul().ldexp(-1, wrk, RoundingMode::None);
        if half.is_negative() != z.is_negative() {
            half = -half;
        }
        half
    } else {
        disk_seed(z, ExactNum::from_u8(1, wrk), wrk)
    };
    let Some((y, yp)) = eval_jet(&coefs, &seed, wrk, p.saturating_add(32)) else {
        return cat_nan();
    };
    let zero = ExactNum::new(wrk);
    let sings = [zero, one, a.clone()];
    match taylor_walk(seed, y, yp, z, &sings, wrk, p, |z0| {
        heun_g_ode(z0, a, q, alpha, beta, gamma, delta, wrk)
    }) {
        Some(v) if finite(&v) => finish(v, p, rm),
        _ => cat_nan(),
    }
}

fn heun_c_outside(
    z: &ExactNum,
    alpha: &ExactNum,
    gamma: &ExactNum,
    delta: &ExactNum,
    eps: &ExactNum,
    q: &ExactNum,
    p: usize,
    rm: RoundingMode,
) -> ExactNum {
    let wrk = work_p(p);
    let Some((coefs, poly)) = heun_c_coefs(alpha, gamma, delta, eps, q, wrk) else {
        return cat_nan();
    };
    if poly {
        return finish(horner(&coefs, z, wrk), p, rm);
    }
    let one = ExactNum::from_u8(1, wrk);
    if !z.is_negative() {
        return cat_nan();
    }
    let seed = disk_seed(z, ExactNum::from_u8(1, wrk), wrk);
    let Some((y, yp)) = eval_jet(&coefs, &seed, wrk, p.saturating_add(32)) else {
        return cat_nan();
    };
    let sings = [ExactNum::new(wrk), one];
    match taylor_walk(seed, y, yp, z, &sings, wrk, p, |z0| {
        heun_c_ode(z0, alpha, gamma, delta, eps, q, wrk)
    }) {
        Some(v) if finite(&v) => finish(v, p, rm),
        _ => cat_nan(),
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
            let num = ExactNum::from_i32(k, p).add(ai, p, RoundingMode::None);
            a_tilde.push(num.div(&aif, p, RoundingMode::None));
        }
    }
    let mut b_tilde = Vec::new();
    for (bj, &bj_int) in b.iter().zip(b_int.iter()) {
        let bif = ExactNum::from_i32(bj_int, p);
        for k in 0..bj_int {
            let num = ExactNum::from_i32(k, p).add(bj, p, RoundingMode::None);
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
    let half = ExactNum::from_u8(1, p).div(&two, p, RoundingMode::None);
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
        let expn = ExactNum::sub_ref(bj, &half, p, RoundingMode::None);
        m_fac = m_fac.mul(&pow_real(&base, &expn, p, cc), p, RoundingMode::None);
    }
    for (ai, &ai_int) in a.iter().zip(a_int.iter()) {
        let base = ExactNum::from_i32(ai_int, p);
        let expn = ExactNum::from_u8(1, p)
            .div(&two, p, RoundingMode::None)
            .neg()
            .add(ai, p, RoundingMode::None);
        m_fac = m_fac.div(&pow_real(&base, &expn, p, cc), p, RoundingMode::None);
    }
    let pi = cc.pi(p, RoundingMode::None);
    let two_pi = ExactNum::from_u8(2, p).mul(&pi, p, RoundingMode::None);
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
    let z_arg = ExactNum::div_ref(z, &beta_r, p, RoundingMode::None);
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
        let d = ExactNum::sub_ref(a, b, p, RoundingMode::None).abs();
        d.is_zero() || d.exponent().is_some_and(|e| e < -((p as Exponent) / 4))
    }

    #[test]
    fn test_meijer_heun_identities() -> Result<(), crate::Error> {
        let p = 128;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, rm);
        let z = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, rm);
        let zero = ExactNum::new(p);
        let one = ExactNum::from_u8(1, p);

        let g = z.meijer_g(&[], &[], core::slice::from_ref(&zero), &[], p, rm, &mut cc);
        let e = z.neg_via_mul().exp(p, rm, &mut cc);
        assert!(near(&g, &e, p), "G^{{1,0}}_{{0,1}}(1/2) = e^{{-1/2}}");

        let fox = z.fox_h(
            &[],
            &[],
            &[(ExactNum::new(p), ExactNum::from_u8(1, p))],
            &[],
            p,
            rm,
            &mut cc,
        );
        assert!(near(&fox, &e, p), "H with A=B=1 is G");

        let a03 = ExactNum::from_u8(3, p).div(&ExactNum::from_u8(10, p), p, rm);
        let g11 = z.meijer_g(
            core::slice::from_ref(&a03),
            &[],
            core::slice::from_ref(&zero),
            &[],
            p,
            rm,
            &mut cc,
        );
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
        let delta = ExactNum::add_ref(&alpha, &beta, p, rm)
            .sub(&gamma, p, rm)
            .add(&one, p, rm);
        let qv = ExactNum::mul_ref(&aa, &alpha, p, rm).mul(&beta, p, rm);
        let hg = zz.heun_g(&aa, &qv, &alpha, &beta, &gamma, &delta, p, rm, &mut cc);
        let f = alpha.hypergeom_2f1(&beta, &gamma, &zz, p, rm, &mut cc);
        assert!(near(&hg, &f, p), "HeunG 2F1 reduction");

        // HeunC → 1F1 when δ=0, ε=-1, q=α.
        let al = ExactNum::from_u8(2, p).div(&ExactNum::from_u8(5, p), p, rm);
        let gam = ExactNum::from_u8(13, p).div(&ExactNum::from_u8(10, p), p, rm);
        let zc = ExactNum::from_u8(7, p).div(&ExactNum::from_u8(20, p), p, rm);
        let hc = zc.heun_c(&al, &gam, &zero, &one.neg_via_mul(), &al, p, rm, &mut cc);
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
        Ok(())
    }
}
