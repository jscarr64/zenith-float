//! Meijer \(G\), Fox \(H\), and local / confluent Heun on [`ExactNum`].
//!
//! SoftFloat only. Residue \(G\) follows DLMF 16.17. Fox \(H\) uses the Gauss
//! multiplication formula when every scale is an exact fraction whose lift
//! fits in memory, and the Mellin–Barnes residue series at the exact scales
//! otherwise. Local Heun is DLMF 31.3; confluent Heun is the \(z=0\)
//! Frobenius solution of DLMF 31.12.1.

use crate::catalog::{cat_nan, finish, finite, small_int, work_p};
use crate::defs::Exponent;
use crate::Consts;
use crate::ExactNum;
use crate::RoundingMode;
use alloc::vec::Vec;
use core::ops::Neg;

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
    let z_pow;
    let z_to_inv_r: &ExactNum = if r_is_one {
        z
    } else {
        z_pow = pow_real(
            z,
            &ExactNum::from_u8(1, p).div(r, p, RoundingMode::None),
            p,
            cc,
        );
        &z_pow
    };
    if !finite(z_to_inv_r) {
        return GEval::Bad;
    }
    let mut sum = ExactNum::new(p);
    if series == 1 {
        let sign_odd = odd_i((pn as i32) - (m as i32) - (n as i32));
        let hz_neg;
        let hz: &ExactNum = if sign_odd {
            hz_neg = z_to_inv_r.neg_ref();
            &hz_neg
        } else {
            z_to_inv_r
        };
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
            ExactNum::from_u8(1, p).div(z_to_inv_r, p, RoundingMode::None)
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
    fn scalar(&self, x: ExactNum) -> Jet {
        let mut j = self.zero();
        if let Some(i) = self.idx(0) {
            j.c[i] = x;
        }
        j
    }
    /// \(x + c\varepsilon\).
    fn lin(&self, x: ExactNum, c: i32) -> Jet {
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
                .map(|(x, y)| {
                    self.note(y);
                    ExactNum::add_ref(x, y, self.wp, RoundingMode::None)
                })
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
        let mut r = self.scalar(ExactNum::from_u8(1, self.wp));
        let mut t = self.scalar(ExactNum::from_u8(1, self.wp));
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
            let n = x.neg_ref().to_u32(RoundingMode::None).ok()?;
            if n > 4096 {
                return None;
            }
            let g1 = self.gamma(&ExactNum::from_u8(1, self.wp), c, cc)?;
            let mut den = self.lin(ExactNum::new(self.wp), c);
            for i in 1..=n {
                den = self.mul(&den, &self.lin(ExactNum::from_u32(i, self.wp).neg(), c));
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
    let mut y = x.abs_ref();
    let _ = y.set_precision(64, RoundingMode::Up);
    y
}

/// `|x|` rounded down to 64 bits.
fn dn64(x: &ExactNum) -> ExactNum {
    let mut y = x.abs_ref();
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
    let f = ExactNum::sub_ref(b, &b.floor_ref(), wp, RoundingMode::None);
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
    let mut t = w.scalar(ExactNum::from_u8(1, wp));
    let mut sum = w.scalar(ExactNum::from_u8(1, wp));
    let mut mk = ExactNum::from_u8(1, 64);
    let mut k: u64 = 0;
    loop {
        let k_next = k.checked_add(1)?;
        let kk = ExactNum::from_u64(k, wp);
        let mut num = w.scalar(ExactNum::from_u8(1, wp));
        // |x| prod(...); `None` while it is still |x| (read in place).
        let mut nb_slot: Option<ExactNum> = None;
        for (a, c) in hn {
            let v = ExactNum::add_ref(a, &kk, wp, RoundingMode::None);
            let f = up64(&v).add(&crho(*c), 64, up);
            nb_slot = Some(ExactNum::mul_ref(
                nb_slot.as_ref().unwrap_or(&ax),
                &f,
                64,
                up,
            ));
            num = w.mul(&num, &w.lin(v, *c));
        }
        let nb = nb_slot.as_ref().unwrap_or(&ax);
        let mut den = w.scalar(ExactNum::from_u64(k_next, wp));
        let mut db = ExactNum::from_u64(k_next, 64);
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
            den = w.mul(&den, &w.lin(v, *c));
        }
        mk = mk.mul(nb, 64, up).div(&db, 64, up);
        t = w.scale(&w.div(&w.mul(&t, &num), &den)?, x);
        sum = w.add(&sum, &t);
        k = k_next;
        if k < kmin {
            continue;
        }
        // Ratio bound for every step j >= k: |x| prod (j + max(|a|+|c|rho, 1)) /
        // ((j+1) prod (j - |b| - |c|rho)); each factor is non-increasing in j.
        let k_bound = k.checked_add(1)?;
        let kk = ExactNum::from_u64(k, 64);
        // `None` while the bound is still |x| (read in place).
        let mut bnd_slot: Option<ExactNum> = None;
        for (a, c) in hn {
            let al = up64(a).add(&crho(*c), 64, up);
            let al = if lt(&al, &one64) { ExactNum::from_u8(1, 64) } else { al };
            let f = ExactNum::add_ref(&kk, &al, 64, up);
            bnd_slot = Some(ExactNum::mul_ref(
                bnd_slot.as_ref().unwrap_or(&ax),
                &f,
                64,
                up,
            ));
        }
        let mut d = ExactNum::from_u64(k_bound, 64);
        for (b, c) in hd {
            let q = ExactNum::sub_ref(&kk, &up64(b), 64, dn).sub(&crho(*c), 64, dn);
            if !q.is_positive() || q.is_zero() {
                return None;
            }
            d = d.mul(&q, 64, dn);
        }
        let bnd = ExactNum::div_ref(bnd_slot.as_ref().unwrap_or(&ax), &d, 64, up);
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
/// \(|G|\). If the first-order slack `wp - loss` stops improving, the result is
/// NaN. A jet that evaluates to zero is returned as zero once the same bound
/// certifies \(|G| < 2^{-(p+64)}\).
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
    let mut prev_slack: Option<i64> = None;
    loop {
        let wp = p.max(span).saturating_add(extra);
        let one = ExactNum::from_u8(1, wp);
        // w and ln w at wp.
        let lnz = z.ln(wp, RoundingMode::None, cc);
        let lnw0 = if r_is_one { lnz } else { lnz.div(r, wp, RoundingMode::None) };
        let lnw = if invert { lnw0.neg() } else { lnw0 };
        let w = lnw.exp(wp, RoundingMode::None, cc);
        let w_exp = w;
        let w: &ExactNum = if r_is_one && !invert { z } else { &w_exp };
        if !finite(w) || !finite(&lnw) {
            return cat_nan();
        }
        let x_neg;
        let x: &ExactNum = if odd_i(a.len() as i32 - m as i32 - n as i32) {
            x_neg = w.neg_ref();
            &x_neg
        } else {
            w
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
            let jw = JetW::new(-poles - 1, hi, wp);
            let mut num = jw.scalar(ExactNum::from_u8(1, wp));
            for (v, c) in &gn {
                let Some(gj) = jw.gamma(v, *c, cc) else {
                    return cat_nan();
                };
                num = jw.mul(&num, &gj);
            }
            let mut den = jw.scalar(ExactNum::from_u8(1, wp));
            for (v, c) in &gd {
                let Some(gj) = jw.gamma(v, *c, cc) else {
                    return cat_nan();
                };
                den = jw.mul(&den, &gj);
            }
            let Some(ratio) = jw.div(&num, &den) else {
                return cat_nan();
            };
            let wb = pow_real(w, bh, wp, cc);
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
            let Some(f) = pfq_jet(&hn, &hd, x, &jw, &weights, tol) else {
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
            let slack = (wp as i64).saturating_sub(loss);
            if prev_slack.is_some_and(|prev| slack <= prev) {
                return cat_nan();
            }
            prev_slack = Some(slack);
            extra = extra
                .saturating_add(usize::try_from(loss.clamp(64, 1 << 20)).unwrap_or(64))
                .saturating_add(64);
        } else {
            // Certified zero: ops * 2^{mag - wp} < 2^{-(p+64)}.
            let ob = 64 - i64::from(ops.leading_zeros());
            let bound_exp = mag.saturating_add(ob).saturating_sub(wp as i64);
            if bound_exp < -(p as i64) - 64 {
                return g;
            }
            let next = extra.saturating_mul(2);
            if next <= extra {
                return cat_nan();
            }
            extra = next;
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
    let mut best: Option<(&ExactNum, i64)> = None;
    for b0 in bm {
        let mut cands: Vec<&ExactNum> = bm
            .iter()
            .chain(ap.iter())
            .chain(bq.iter())
            .filter(|x| int_diff(x, b0))
            .collect();
        cands.sort_by(|x, y| {
            x.compare(y)
                .map_or(core::cmp::Ordering::Equal, |c| c.cmp(&0))
        });
        for s in cands {
            let f = order(s);
            if f >= 1 {
                if best.as_ref().is_none_or(|(bs, _)| lt(s, bs)) {
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
    let even = |k: &ExactNum| k.ldexp_ref(-1, wp, RoundingMode::None).is_int();
    let mut c = ExactNum::from_u8(1, wp);
    for bj in bm {
        // Γ(b - s) near s*: argument -k gives -(-1)^k/k! per unit of (s - s*).
        let v = ExactNum::sub_ref(bj, s, wp, RoundingMode::None);
        let t = if gamma_pole(&v) {
            let k = v.neg_ref();
            let q = ExactNum::div_ref(&one, &fact(&k, cc), wp, RoundingMode::None);
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
        let v = ExactNum::sub_ref(&one, aj, wp, RoundingMode::None).add(s, wp, RoundingMode::None);
        if gamma_pole(&v) {
            return cat_nan();
        }
        c = c.mul(&v.gamma(wp, RoundingMode::None, cc), wp, RoundingMode::None);
    }
    for aj in ap {
        // 1/Γ(a - s): argument -k gives -(-1)^k k!.
        let v = ExactNum::sub_ref(aj, s, wp, RoundingMode::None);
        let t = if gamma_pole(&v) {
            let k = v.neg_ref();
            let q = fact(&k, cc);
            if even(&k) {
                q.neg()
            } else {
                q
            }
        } else {
            ExactNum::from_u8(1, wp).div(
                &v.gamma(wp, RoundingMode::None, cc),
                wp,
                RoundingMode::None,
            )
        };
        c = c.mul(&t, wp, RoundingMode::None);
    }
    for bj in bq {
        // 1/Γ(1 - b + s): argument -k gives (-1)^k k!.
        let v = ExactNum::sub_ref(&one, bj, wp, RoundingMode::None).add(s, wp, RoundingMode::None);
        let t = if gamma_pole(&v) {
            let k = v.neg_ref();
            let q = fact(&k, cc);
            if even(&k) {
                q
            } else {
                q.neg()
            }
        } else {
            ExactNum::from_u8(1, wp).div(
                &v.gamma(wp, RoundingMode::None, cc),
                wp,
                RoundingMode::None,
            )
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
    let w_pow;
    let w: &ExactNum = if r.compare(&one) == Some(0) {
        z
    } else {
        w_pow = pow_real(
            z,
            &ExactNum::from_u8(1, p).div(r, p, RoundingMode::None),
            p,
            cc,
        );
        &w_pow
    };
    if !finite(w) {
        return cat_nan();
    }
    if choose_series(an.len(), ap.len(), bm.len(), bq.len(), w, p) == 1 {
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

fn gcd_u(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

fn lcm_u(a: u64, b: u64) -> Option<u64> {
    if a == 0 || b == 0 {
        return Some(0);
    }
    let g = gcd_u(a, b);
    a.checked_div(g)?.checked_mul(b)
}

/// The exact positive dyadic \(n/d = x\) in lowest terms, when both \(n\)
/// and \(d\) fit in `u64`. No slack: a nearby rational is never accepted.
/// A generic dest-\(p\) scale has a denominator \(2^{O(p)}\) and returns
/// `None` here, so Fox \(H\) uses the residue series at the exact input.
fn as_pos_rational(x: &ExactNum, _p: usize) -> Option<(u64, u64)> {
    exact_pos_frac(x)
}

/// Exact \(x = n/d\) in lowest terms with \(n,d\in\) `u64`. `None` if \(x\)
/// is not a positive finite number or the reduced numerator or denominator
/// does not fit in `u64`.
fn exact_pos_frac(x: &ExactNum) -> Option<(u64, u64)> {
    if !finite(x) || !x.is_positive() || x.is_zero() {
        return None;
    }
    let words = x.mantissa_digits()?;
    let e = x.exponent()?;
    if words.is_empty() {
        return None;
    }
    let mut tz_words = 0usize;
    let mut first_nz = 0u64;
    for &w in words {
        if w == 0 {
            tz_words += 1;
        } else {
            first_nz = w;
            break;
        }
    }
    if first_nz == 0 {
        return None;
    }
    let tz = tz_words
        .checked_mul(64)?
        .checked_add(first_nz.trailing_zeros() as usize)?;
    let n = shr_words_u64(words, tz)?;
    let exp2 = e
        .checked_sub((64 * words.len()) as i64)?
        .checked_add(tz as i64)?;
    if exp2 >= 0 {
        if exp2 >= 64 {
            return None;
        }
        let n = n.checked_shl(exp2 as u32)?;
        if n == 0 {
            return None;
        }
        Some((n, 1))
    } else {
        let k = (-exp2) as u32;
        if k >= 64 {
            return None;
        }
        let d = 1u64.checked_shl(k)?;
        if n == 0 || d == 0 {
            return None;
        }
        let g = gcd_u(n, d);
        Some((n / g, d / g))
    }
}

/// Integer value of the little-endian word mantissa shifted right by `tz`
/// bits, if it fits in `u64`.
fn shr_words_u64(words: &[u64], tz: usize) -> Option<u64> {
    let woff = tz / 64;
    let b = tz % 64;
    let rest = words.get(woff..)?;
    if rest.is_empty() {
        return None;
    }
    let lo = rest[0] >> b;
    let hi = if b == 0 || rest.len() < 2 { 0 } else { rest[1] << (64 - b) };
    let n = lo | hi;
    if b == 0 {
        if rest[1..].iter().any(|&w| w != 0) {
            return None;
        }
    } else {
        if rest.len() >= 2 && rest[1] >> b != 0 {
            return None;
        }
        if rest.get(2..).is_some_and(|t| t.iter().any(|&w| w != 0)) {
            return None;
        }
    }
    if n == 0 {
        None
    } else {
        Some(n)
    }
}

fn all_unit_scales(xs: &[(ExactNum, ExactNum)], p: usize) -> bool {
    let one = ExactNum::from_u8(1, p);
    xs.iter().all(|(_, s)| s.compare(&one) == Some(0))
}

fn pow_u64(base: &ExactNum, n: u64, p: usize, cc: &mut Consts) -> ExactNum {
    if n <= i64::MAX as u64 {
        base.powsi(n as i64, p, RoundingMode::None)
    } else {
        pow_real(base, &ExactNum::from_u64(n, p), p, cc)
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
    ///   tail bound.
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
    /// [`meijer_g`]. When every \(A_i=B_j=1\) this is Meijer \(G\).
    ///
    /// # Precision
    ///
    /// - Algorithm: unit scales call [`ExactNum::meijer_g`]. If every scale is
    ///   an exact fraction \(n/d\) (the reduced dyadic equals the input, with
    ///   no slack) and the Gauss lift fits in `usize` and in memory, the lift
    ///   is used. Every other case sums the Mellin–Barnes residue series at
    ///   the exact scales (DLMF 16.17 generalizes in the obvious way to
    ///   \(A_i,B_j\neq 1\)). Extra word of working precision.
    /// - Bound: same as Meijer \(G\) after a lift; the residue series stops
    ///   when the term is below a unit in the last place of the working
    ///   precision, or stalls if the terms stop shrinking.
    ///
    /// # Limitations
    ///
    /// - Scales may be any positive reals (exact dyadic inputs). There is no
    ///   approximate or slack-based identification: a dest-\(p\) rounding of
    ///   \(2/17\), or \(2/17+2^{-200}\), is evaluated at that exact value,
    ///   not at \(2/17\). A Gauss lift that does not fit in memory or in a
    ///   `usize` index is not attempted; the residue series is used instead.
    /// - Coincident poles after a lift use the same limit as [`meijer_g`].
    ///   Coincident poles of the residue series are `NaN`.
    /// - Real \(\mathrm{self}>0\) when a scale is not 1 (fractional \(z^{1/B}\)).
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
        if let Some(y) = fox_h_gauss(self, an, ap, bm, bq, wrk, cc) {
            return finish(y, p, rm);
        }
        match fox_h_residue(self, an, ap, bm, bq, wrk, cc) {
            GEval::Val(y) => finish(y, p, rm),
            GEval::Bad | GEval::Singular => cat_nan(),
        }
    }

    /// Local Heun \(\mathrm{Hl}(a,q;\alpha,\beta,\gamma,\delta;\mathrm{self})\) (DLMF 31.3).
    ///
    /// Power series about \(z=0\) with \(c_0=1\), \(c_1=q/(a\gamma)\), and the
    /// three-term recurrence 31.3.3–4. \(\varepsilon=\alpha+\beta-\gamma-\delta+1\).
    ///
    /// # Precision
    ///
    /// - Algorithm: DLMF 31.3 power series inside \(\lvert z\rvert<\min(1,\lvert a\rvert)\),
    ///   summed until a majorant of the three-term recurrence (roots \(1\), \(1/a\))
    ///   bounds the tail; the working precision grows with the largest term.
    ///   Outside that disk, a seed jet in the disk uses the same majorant; then
    ///   Taylor steps of the ODE along the real segment, each step at most half
    ///   the distance to \(\{0,1,a\}\). Each step carries a midpoint and a
    ///   radius from the 2×2 local transition and a majorant remainder. Working
    ///   precision is raised until that enclosure uniquely determines the dest-p
    ///   rounding in the requested mode. A walk that meets a cut is `NaN`.
    /// - Bound (inside the disk): working precision `p + 2 WORD_BIT_SIZE`; not
    ///   Ziv-certified. Bound (outside the disk): the ball enclosure; Ziv on
    ///   working precision until unique rounding.
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
    /// - A Taylor walk that cannot advance (step does not change the
    ///   abscissa, or a remainder bound cannot be formed) is `NaN(InvalidArgument)`.
    /// - \(a=0\) or \(\gamma\in\{0,-1,-2,\ldots\}\) (recurrence pole) is `NaN`.
    ///
    /// Accuracy: correctly-rounded outside the disk; unspecified inside the disk
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
        heun_g_series(self, a, q, alpha, beta, gamma, delta, p, rm)
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
    /// - Algorithm: Frobenius series of DLMF 31.12.1 about \(z=0\) for \(\lvert z\rvert<1\),
    ///   summed until a majorant of the recurrence (roots \(1\), \(0\)) bounds the tail.
    ///   For real \(z<1\) outside that disk, a seed jet in the disk uses the same
    ///   majorant; then Taylor steps along \((-\infty,1)\), each step at most half
    ///   the distance to \(\{0,1\}\). Each step carries a midpoint and a radius
    ///   from the 2×2 local transition and a majorant remainder. Working precision
    ///   is raised until that enclosure uniquely determines the dest-p rounding
    ///   in the requested mode. A walk that meets a cut is `NaN`.
    /// - Bound (inside the disk): working precision `p + 2 WORD_BIT_SIZE`; not
    ///   Ziv-certified. Bound (outside the disk): the ball enclosure; Ziv on
    ///   working precision until unique rounding.
    ///
    /// # Limitations
    ///
    /// - \(z\ge 1\) meets the branch point at \(1\) and is `NaN(InvalidArgument)`
    ///   unless two consecutive series coefficients are exact zeros (a polynomial,
    ///   returned for every finite \(z\)) or the \({}_1F_1\) reduction applies.
    /// - A Taylor walk that cannot advance (step does not change the
    ///   abscissa, or a remainder bound cannot be formed) is `NaN(InvalidArgument)`.
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
    /// Accuracy: correctly-rounded outside the disk; unspecified inside the disk
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
        heun_c_series(self, alpha, gamma, delta, eps, q, p, rm)
    }
}

/// Rigorous tail data for a power series \(\sum c_j z^j\) whose coefficients obey
/// \(c_{j+1}=A_jc_j+B_jc_{j-1}\) with \(A_j\to\lambda_d+\lambda_s\), \(B_j\to-\lambda_d\lambda_s\),
/// \(\lvert\lambda_s\rvert\le\lvert\lambda_d\rvert\), and
/// \(A_j-A_\infty=(\ell_1j+\ell_0)/D_j\), \((A_j-A_\infty)\lambda_s+(B_j-B_\infty)=(k_1j+k_0)/D_j\),
/// \(\lvert D_j\rvert\ge s\,j\,(j-g)\) for \(j>g\). All fields are 64-bit upper bounds
/// (`s` a lower bound).
struct RecTail {
    lam_d: ExactNum,
    lam_s: ExactNum,
    /// \(\lambda_s\) itself at the working precision (signed), for \(d_j\).
    lam_s_v: ExactNum,
    /// Constant part of \(\kappa_j\) from rounding \(\lambda_s\) to the working precision.
    k_c: ExactNum,
    l1: ExactNum,
    l0: ExactNum,
    k1: ExactNum,
    k0: ExactNum,
    s: ExactNum,
    g: ExactNum,
}

/// `|x|` as a 64-bit upper bound with slack for the rounding of `x` itself.
fn ub64(x: &ExactNum) -> ExactNum {
    let up = RoundingMode::Up;
    let mut y = x.abs_ref();
    let _ = y.set_precision(64, up);
    let slack = ExactNum::add_ref(&y, &ExactNum::from_u8(1, 64), 64, up).ldexp(-60, 64, up);
    y.add(&slack, 64, up)
}

/// Upper bound of a (possibly tiny) positive error radius: relative 2^{-60}
/// pad only, so a 2^{-wp} bound is not inflated to 2^{-60}.
fn ub_rel(x: &ExactNum) -> ExactNum {
    let up = RoundingMode::Up;
    let mut y = x.abs_ref();
    let _ = y.set_precision(64, up);
    let slack = y.ldexp_ref(-60, 64, up);
    y.add(&slack, 64, up)
}

impl RecTail {
    /// For \(r\) and every \(j\ge J\): with \(d_j=c_j-\lambda_sc_{j-1}\),
    /// \(\max(\lvert d_j\rvert r^j, w\lvert c_{j-1}\rvert r^{j-1})\) does not increase, where
    /// \(w=(1-\lvert\lambda_s\rvert r)/r\). Returns `w` when that holds.
    ///
    /// Rows: \(\lvert d_{j+1}\rvert r^{j+1}\le(\lvert\lambda_d\rvert+\lvert\alpha_j\rvert)r\cdot\lvert d_j\rvert r^j
    /// +\lvert\kappa_j\rvert r^2\cdot\lvert c_{j-1}\rvert r^{j-1}\) and
    /// \(\lvert c_j\rvert r^j\le r\lvert d_j\rvert r^j+\lvert\lambda_s\rvert r\cdot\lvert c_{j-1}\rvert r^{j-1}\),
    /// with \(\sup_{j\ge J}\lvert\alpha_j\rvert\le(\ell_1+\ell_0/J)/(s(J-g))\) and the same for \(\kappa\).
    fn weight(&self, r: &ExactNum, jj: u64) -> Option<ExactNum> {
        let up = RoundingMode::Up;
        let dn = RoundingMode::Down;
        let j = ExactNum::from_u64(jj, 64);
        let jg = ExactNum::sub_ref(&j, &self.g, 64, dn);
        if !jg.is_positive() || jg.is_zero() {
            return None;
        }
        let den = ExactNum::mul_ref(&self.s, &jg, 64, dn);
        if !den.is_positive() || den.is_zero() {
            return None;
        }
        let sup = |c1: &ExactNum, c0: &ExactNum| {
            ExactNum::div_ref(c0, &j, 64, up)
                .add(c1, 64, up)
                .div(&den, 64, up)
        };
        let alpha = sup(&self.l1, &self.l0);
        let kappa = sup(&self.k1, &self.k0).add(&self.k_c, 64, up);
        let one = ExactNum::from_u8(1, 64);
        let ls_r = ExactNum::mul_ref(&self.lam_s, r, 64, up);
        let wn = ExactNum::sub_ref(&one, &ls_r, 64, dn);
        if !wn.is_positive() || wn.is_zero() {
            return None;
        }
        let w = wn.div(r, 64, dn);
        let row1 = ExactNum::add_ref(&self.lam_d, &alpha, 64, up)
            .mul(r, 64, up)
            .add(&kappa.mul(r, 64, up).mul(r, 64, up).div(&w, 64, up), 64, up);
        if matches!(row1.compare(&one), Some(c) if c <= 0) {
            Some(w)
        } else {
            None
        }
    }
}

/// \(\sum_j c_jz^j\) for the three-term recurrence `step(j) = (A_j, B_j)`, \(c_0=1\),
/// \(c_1\) given, \(\lvert z\rvert<1/\lvert\lambda_d\rvert\). Terms are added until the
/// [`RecTail`] majorant bounds the rest below \(2^{-(\mathrm{wp}+8)}\) of the sum:
/// with \(r=(\lvert z\rvert+1/\lvert\lambda_d\rvert)/2\), \(t=\lvert z\rvert/r\) and \(M_N\) the
/// weighted norm at the last index, \(\sum_{k>N}\lvert c_kz^k\rvert\le(M_N/w)\,t^{N+1}/(1-t)\).
/// There is no term cap. Returns the sum and the rounding loss in bits (largest
/// term against the sum, plus the term count and the majorant's amplification),
/// or `None` when the bound cannot be formed.
struct RecJet {
    val: ExactNum,
    deriv: ExactNum,
    loss: i64,
    poly: bool,
    coefs: Vec<ExactNum>,
}

fn rec_series<F>(
    z: &ExactNum,
    c1: &ExactNum,
    tail: &RecTail,
    wp: usize,
    step: F,
) -> Option<(ExactNum, i64)>
where
    F: FnMut(u64, usize) -> Option<(ExactNum, ExactNum)>,
{
    rec_series_inner(z, c1, tail, wp, false, step).map(|s| (s.val, s.loss))
}

fn rec_series_jet<F>(
    z: &ExactNum,
    c1: &ExactNum,
    tail: &RecTail,
    wp: usize,
    step: F,
) -> Option<RecJet>
where
    F: FnMut(u64, usize) -> Option<(ExactNum, ExactNum)>,
{
    rec_series_inner(z, c1, tail, wp, true, step)
}

fn rec_series_inner<F>(
    z: &ExactNum,
    c1: &ExactNum,
    tail: &RecTail,
    wp: usize,
    jet: bool,
    mut step: F,
) -> Option<RecJet>
where
    F: FnMut(u64, usize) -> Option<(ExactNum, ExactNum)>,
{
    let up = RoundingMode::Up;
    let dn = RoundingMode::Down;
    let nr = RoundingMode::None;
    let one64 = ExactNum::from_u8(1, 64);
    let az = ub64(z);
    let mut az_dn = z.abs_ref();
    let _ = az_dn.set_precision(64, dn);
    let slack = ExactNum::add_ref(&az_dn, &one64, 64, up).ldexp(-60, 64, up);
    az_dn = az_dn.sub(&slack, 64, dn);
    if az_dn.is_negative() || az_dn.is_zero() {
        az_dn = ExactNum::new(64);
    }
    let rho = ExactNum::div_ref(&one64, &tail.lam_d, 64, dn);
    let r = ExactNum::add_ref(&az, &rho, 64, dn).ldexp(-1, 64, dn);
    if !matches!(az.compare(&r), Some(c) if c < 0) {
        return None;
    }
    let t = ExactNum::div_ref(&az, &r, 64, up);
    let geo = ExactNum::sub_ref(&one64, &t, 64, dn);
    if !geo.is_positive() || geo.is_zero() {
        return None;
    }
    let lam_s_wp = &tail.lam_s_v;
    let mut c_prev = ExactNum::from_u8(1, wp);
    let mut c = c1.clone();
    // z^j; `None` while j = 1 (z read in place).
    let mut zp_slot: Option<ExactNum> = None;
    let mut sum = ExactNum::from_u8(1, wp).add(&ExactNum::mul_ref(c1, z, wp, nr), wp, nr);
    let mut dsum = c1.clone();
    let mut peak = sum.exponent().unwrap_or(0).max(0);
    let mut peak_d = dsum.exponent().unwrap_or(0).max(0);
    let mut coefs = Vec::new();
    if jet {
        coefs.push(c_prev.clone());
        coefs.push(c.clone());
    }
    // log2 r^j and t^j tracked through exponents of 64-bit powers.
    // r^j and t^j; `None` while j = 1 (r and t read in place).
    let mut rj_slot: Option<ExactNum> = None;
    let mut tj_slot: Option<ExactNum> = None;
    let mut j: u64 = 1;
    let mut w_ok: Option<ExactNum> = None;
    loop {
        // c holds c_j, c_prev c_{j-1}; rj = r^j, tj = t^j.
        let rj = rj_slot.as_ref().unwrap_or(&r);
        let tj = tj_slot.as_ref().unwrap_or(&t);
        if w_ok.is_none() {
            w_ok = tail.weight(&r, j);
        }
        if let Some(w) = &w_ok {
            // |d_j| from the working-precision difference plus its rounding.
            let ls_c = ExactNum::mul_ref(lam_s_wp, &c_prev, wp, nr);
            let d = ExactNum::sub_ref(&c, &ls_c, wp, nr);
            let rnd = ub64(&c)
                .add(&ub64(&ls_c), 64, up)
                .ldexp(4 - wp as i64, 64, up);
            let dn_ = ub64(&d).add(&rnd, 64, up).mul(rj, 64, up);
            let cp = ub64(&c_prev).mul(rj, 64, up).div(&r, 64, up).mul(w, 64, up);
            let m = if matches!(dn_.compare(&cp), Some(o) if o > 0) { dn_ } else { cp };
            let bound = m
                .div(w, 64, up)
                .mul(tj, 64, up)
                .mul(&t, 64, up)
                .div(&geo, 64, up);
            let es = sum.exponent().unwrap_or(0);
            let val_ok = bound.is_zero()
                || (!sum.is_zero() && bound.exponent().is_some_and(|e| e < es - wp as i64 - 8));
            let deriv_ok = if !jet {
                true
            } else if az_dn.is_zero() || !finite(&dsum) {
                false
            } else {
                // |tail of y'| ≤ bound/|z| · (j + 1/(1-t)); f(t)=j+1/(1-t) grows with t.
                let fac =
                    ExactNum::from_u64(j, 64).add(&ExactNum::div_ref(&one64, &geo, 64, up), 64, up);
                let db = bound.div(&az_dn, 64, up).mul(&fac, 64, up);
                let ed = dsum.exponent().unwrap_or(0);
                db.is_zero()
                    || (!dsum.is_zero() && db.exponent().is_some_and(|e| e < ed - wp as i64 - 8))
            };
            if val_ok && deriv_ok {
                let n_bits = 64 - i64::from(j.leading_zeros());
                let amp = ExactNum::div_ref(&one64, w, 64, up)
                    .add(w, 64, up)
                    .div(&geo, 64, up);
                let mut loss = (peak - es).max(0) + n_bits + amp.exponent().unwrap_or(0).max(0) + 4;
                if jet {
                    let ed = dsum.exponent().unwrap_or(0);
                    let loss_d =
                        (peak_d - ed).max(0) + n_bits + amp.exponent().unwrap_or(0).max(0) + 4;
                    if loss_d > loss {
                        loss = loss_d;
                    }
                }
                return Some(RecJet {
                    val: sum,
                    deriv: dsum,
                    loss,
                    poly: false,
                    coefs,
                });
            }
        }
        let (aj, bj) = step(j, wp)?;
        let c_next = ExactNum::mul_ref(&aj, &c, wp, nr).add(
            &ExactNum::mul_ref(&bj, &c_prev, wp, nr),
            wp,
            nr,
        );
        if !finite(&c_next) {
            return None;
        }
        if jet && c.is_zero() && c_next.is_zero() {
            return Some(RecJet {
                val: sum,
                deriv: dsum,
                loss: 0,
                poly: true,
                coefs,
            });
        }
        let zp = ExactNum::mul_ref(zp_slot.as_ref().unwrap_or(z), z, wp, nr);
        let term = ExactNum::mul_ref(&c_next, &zp, wp, nr);
        zp_slot = Some(zp);
        if let (false, Some(e)) = (term.is_zero(), term.exponent()) {
            peak = peak.max(e);
        }
        sum = sum.add(&term, wp, nr);
        if jet {
            let dterm = ExactNum::from_u64(j + 1, wp).mul(&term.div(z, wp, nr), wp, nr);
            if !finite(&dterm) {
                return None;
            }
            if let (false, Some(e)) = (dterm.is_zero(), dterm.exponent()) {
                peak_d = peak_d.max(e);
            }
            dsum = dsum.add(&dterm, wp, nr);
            coefs.push(c_next.clone());
        }
        c_prev = c;
        c = c_next;
        let (rj, tj) = (
            ExactNum::mul_ref(rj, &r, 64, up),
            ExactNum::mul_ref(tj, &t, 64, up),
        );
        j += 1;
        let stop = rj.is_zero() || tj.is_zero();
        rj_slot = Some(rj);
        tj_slot = Some(tj);
        if stop {
            // 64-bit powers underflow only far past any useful length.
            return None;
        }
    }
}

/// Retry `f(wp)` with more precision until its reported loss leaves `p + 16` bits.
fn with_loss<F>(p: usize, mut f: F) -> Option<ExactNum>
where
    F: FnMut(usize) -> Option<(ExactNum, i64)>,
{
    let mut wp = work_p(p);
    let mut prev_slack: Option<i64> = None;
    loop {
        let (v, loss) = f(wp)?;
        let wp_i = i64::try_from(wp).ok()?;
        let need = i64::try_from(p)
            .ok()?
            .saturating_add(16)
            .saturating_add(loss);
        if need <= wp_i {
            return Some(v);
        }
        let slack = wp_i.saturating_sub(loss);
        if prev_slack.is_some_and(|s| slack <= s) {
            return None;
        }
        prev_slack = Some(slack);
        wp = usize::try_from(need)
            .ok()?
            .saturating_add(crate::defs::WORD_BIT_SIZE);
        if wp == 0 {
            return None;
        }
    }
}

fn with_loss_jet<F>(p: usize, mut f: F) -> Option<(ExactNum, ExactNum, bool, Vec<ExactNum>)>
where
    F: FnMut(usize) -> Option<RecJet>,
{
    let mut wp = work_p(p);
    let mut prev_slack: Option<i64> = None;
    loop {
        let s = f(wp)?;
        if s.poly {
            return Some((s.val, s.deriv, true, s.coefs));
        }
        let wp_i = i64::try_from(wp).ok()?;
        let need = i64::try_from(p)
            .ok()?
            .saturating_add(16)
            .saturating_add(s.loss);
        if need <= wp_i {
            return Some((s.val, s.deriv, false, s.coefs));
        }
        let slack = wp_i.saturating_sub(s.loss);
        if prev_slack.is_some_and(|s| slack <= s) {
            return None;
        }
        prev_slack = Some(slack);
        wp = usize::try_from(need)
            .ok()?
            .saturating_add(crate::defs::WORD_BIT_SIZE);
        if wp == 0 {
            return None;
        }
    }
}

/// DLMF 31.3 series of \(\mathrm{Hl}\) for \(\lvert z\rvert<\min(1,\lvert a\rvert)\) with the
/// [`rec_series`] tail bound (characteristic roots \(1\) and \(1/a\)).
fn heun_g_series(
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
    let nr = RoundingMode::None;
    let v = with_loss(p, |wp| {
        let one = ExactNum::from_u8(1, wp);
        let eps = heun_eps(alpha, beta, gamma, delta, wp);
        let opa = ExactNum::add_ref(&one, a, wp, nr);
        let inv_a = ExactNum::div_ref(&one, a, wp, nr);
        let big_a = matches!(a.abs_cmp(&one), Some(c) if c >= 0);
        let (lam_d, lam_s) = if big_a {
            (ExactNum::from_u8(1, wp), inv_a)
        } else {
            (inv_a, ExactNum::from_u8(1, wp))
        };
        // n1 j + n0 = Q_j + q - (1+a)(j+1)(j+γ); m1 j + m0 = P_j - (j+1)(j+γ).
        let n1 = ExactNum::mul_ref(a, delta, wp, nr).add(&eps, wp, nr).sub(
            &opa.ldexp_ref(1, wp, nr),
            wp,
            nr,
        );
        let n0 = ExactNum::sub_ref(q, &ExactNum::mul_ref(&opa, gamma, wp, nr), wp, nr);
        let m1 = ExactNum::add_ref(alpha, beta, wp, nr)
            .sub(&ExactNum::from_u8(3, wp), wp, nr)
            .sub(gamma, wp, nr);
        let m0 = ExactNum::sub_ref(alpha, &one, wp, nr)
            .mul(&ExactNum::sub_ref(beta, &one, wp, nr), wp, nr)
            .sub(gamma, wp, nr);
        let k1 = ExactNum::mul_ref(&lam_s, &n1, wp, nr).sub(&m1, wp, nr);
        let k0 = ExactNum::mul_ref(&lam_s, &n0, wp, nr).sub(&m0, wp, nr);
        let mut s = a.abs_ref();
        let _ = s.set_precision(64, RoundingMode::Down);
        let tail = RecTail {
            lam_d: ub64(&lam_d),
            lam_s: ub64(&lam_s),
            k_c: ub64(&lam_d)
                .add(&ub64(&lam_s), 64, RoundingMode::Up)
                .add(&ub64(&n1), 64, RoundingMode::Up)
                .add(&ub64(&n0), 64, RoundingMode::Up)
                .ldexp(8 - wp as i64, 64, RoundingMode::Up),
            lam_s_v: lam_s,
            l1: ub64(&n1),
            l0: ub64(&n0),
            k1: ub64(&k1),
            k0: ub64(&k0),
            s,
            g: ub64(gamma),
        };
        let c1 = ExactNum::div_ref(q, &ExactNum::mul_ref(a, gamma, wp, nr), wp, nr);
        if !finite(&c1) {
            return None;
        }
        rec_series(z, &c1, &tail, wp, |j, wp| {
            let jf = ExactNum::from_u64(j, wp);
            let jm1 = ExactNum::from_u64(j - 1, wp);
            let pj = ExactNum::add_ref(&jm1, alpha, wp, nr).mul(
                &ExactNum::add_ref(&jm1, beta, wp, nr),
                wp,
                nr,
            );
            let qj = ExactNum::add_ref(&jm1, gamma, wp, nr)
                .mul(&opa, wp, nr)
                .add(&ExactNum::mul_ref(a, delta, wp, nr), wp, nr)
                .add(&eps, wp, nr)
                .mul(&jf, wp, nr);
            let rj = ExactNum::from_u64(j + 1, wp).mul(a, wp, nr).mul(
                &ExactNum::add_ref(&jf, gamma, wp, nr),
                wp,
                nr,
            );
            if rj.is_zero() || !finite(&rj) {
                return None;
            }
            let aj = qj.add(q, wp, nr).div(&rj, wp, nr);
            let bj = pj.div(&rj, wp, nr).neg();
            Some((aj, bj))
        })
    });
    match v {
        Some(v) => finish(v, p, rm),
        None => cat_nan(),
    }
}

/// DLMF 31.12.1 Frobenius series of confluent Heun for \(\lvert z\rvert<1\) with the
/// [`rec_series`] tail bound (characteristic roots \(1\) and \(0\)).
fn heun_c_series(
    z: &ExactNum,
    alpha: &ExactNum,
    gamma: &ExactNum,
    delta: &ExactNum,
    eps: &ExactNum,
    q: &ExactNum,
    p: usize,
    rm: RoundingMode,
) -> ExactNum {
    let nr = RoundingMode::None;
    let v = with_loss(p, |wp| {
        let one = ExactNum::from_u8(1, wp);
        // l1 n + l0 = n(n-1+γ+δ-ε) - q - (n+1)(n+γ); κ_n = B_n = (ε n + α - ε)/((n+1)(n+γ)).
        let l1 = ExactNum::sub_ref(delta, eps, wp, nr).sub(&ExactNum::from_u8(2, wp), wp, nr);
        let l0 = ExactNum::add_ref(q, gamma, wp, nr);
        let k0 = ExactNum::sub_ref(alpha, eps, wp, nr);
        let tail = RecTail {
            lam_d: ub64(&one),
            lam_s: ExactNum::new(64),
            lam_s_v: ExactNum::new(wp),
            k_c: ExactNum::new(64),
            l1: ub64(&l1),
            l0: ub64(&l0),
            k1: ub64(eps),
            k0: ub64(&k0),
            s: ExactNum::from_u8(1, 64),
            g: ub64(gamma),
        };
        let b1 = q.neg_via_mul().div(gamma, wp, nr);
        if !finite(&b1) {
            return None;
        }
        rec_series(z, &b1, &tail, wp, |n, wp| {
            let nf = ExactNum::from_u64(n, wp);
            let nm1 = ExactNum::from_u64(n - 1, wp);
            let left =
                ExactNum::from_u64(n + 1, wp).mul(&ExactNum::add_ref(&nf, gamma, wp, nr), wp, nr);
            if left.is_zero() || !finite(&left) {
                return None;
            }
            let coef_n = ExactNum::add_ref(&nm1, gamma, wp, nr)
                .add(delta, wp, nr)
                .sub(eps, wp, nr)
                .mul(&nf, wp, nr)
                .sub(q, wp, nr);
            let coef_nm1 = ExactNum::mul_ref(&nm1, eps, wp, nr).add(alpha, wp, nr);
            Some((coef_n.div(&left, wp, nr), coef_nm1.div(&left, wp, nr)))
        })
    });
    match v {
        Some(v) => finish(v, p, rm),
        None => cat_nan(),
    }
}

/// Seed value and derivative of HeunG at a point inside the disk, using the
/// same [`rec_series`] majorant as [`heun_g_series`]. A terminating series
/// is returned as a polynomial.
fn heun_g_jet(
    z: &ExactNum,
    a: &ExactNum,
    q: &ExactNum,
    alpha: &ExactNum,
    beta: &ExactNum,
    gamma: &ExactNum,
    delta: &ExactNum,
    p: usize,
) -> Option<(ExactNum, ExactNum, bool, Vec<ExactNum>)> {
    let nr = RoundingMode::None;
    with_loss_jet(p, |wp| {
        let one = ExactNum::from_u8(1, wp);
        let eps = heun_eps(alpha, beta, gamma, delta, wp);
        let opa = ExactNum::add_ref(&one, a, wp, nr);
        let inv_a = ExactNum::div_ref(&one, a, wp, nr);
        let big_a = matches!(a.abs_cmp(&one), Some(c) if c >= 0);
        let (lam_d, lam_s) = if big_a {
            (ExactNum::from_u8(1, wp), inv_a)
        } else {
            (inv_a, ExactNum::from_u8(1, wp))
        };
        let n1 = ExactNum::mul_ref(a, delta, wp, nr).add(&eps, wp, nr).sub(
            &opa.ldexp_ref(1, wp, nr),
            wp,
            nr,
        );
        let n0 = ExactNum::sub_ref(q, &ExactNum::mul_ref(&opa, gamma, wp, nr), wp, nr);
        let m1 = ExactNum::add_ref(alpha, beta, wp, nr)
            .sub(&ExactNum::from_u8(3, wp), wp, nr)
            .sub(gamma, wp, nr);
        let m0 = ExactNum::sub_ref(alpha, &one, wp, nr)
            .mul(&ExactNum::sub_ref(beta, &one, wp, nr), wp, nr)
            .sub(gamma, wp, nr);
        let k1 = ExactNum::mul_ref(&lam_s, &n1, wp, nr).sub(&m1, wp, nr);
        let k0 = ExactNum::mul_ref(&lam_s, &n0, wp, nr).sub(&m0, wp, nr);
        let mut s = a.abs_ref();
        let _ = s.set_precision(64, RoundingMode::Down);
        let tail = RecTail {
            lam_d: ub64(&lam_d),
            lam_s: ub64(&lam_s),
            k_c: ub64(&lam_d)
                .add(&ub64(&lam_s), 64, RoundingMode::Up)
                .add(&ub64(&n1), 64, RoundingMode::Up)
                .add(&ub64(&n0), 64, RoundingMode::Up)
                .ldexp(8 - wp as i64, 64, RoundingMode::Up),
            lam_s_v: lam_s,
            l1: ub64(&n1),
            l0: ub64(&n0),
            k1: ub64(&k1),
            k0: ub64(&k0),
            s,
            g: ub64(gamma),
        };
        let c1 = ExactNum::div_ref(q, &ExactNum::mul_ref(a, gamma, wp, nr), wp, nr);
        if !finite(&c1) {
            return None;
        }
        rec_series_jet(z, &c1, &tail, wp, |j, wp| {
            let jf = ExactNum::from_u64(j, wp);
            let jm1 = ExactNum::from_u64(j - 1, wp);
            let pj = ExactNum::add_ref(&jm1, alpha, wp, nr).mul(
                &ExactNum::add_ref(&jm1, beta, wp, nr),
                wp,
                nr,
            );
            let qj = ExactNum::add_ref(&jm1, gamma, wp, nr)
                .mul(&opa, wp, nr)
                .add(&ExactNum::mul_ref(a, delta, wp, nr), wp, nr)
                .add(&eps, wp, nr)
                .mul(&jf, wp, nr);
            let rj = ExactNum::from_u64(j + 1, wp).mul(a, wp, nr).mul(
                &ExactNum::add_ref(&jf, gamma, wp, nr),
                wp,
                nr,
            );
            if rj.is_zero() || !finite(&rj) {
                return None;
            }
            let aj = qj.add(q, wp, nr).div(&rj, wp, nr);
            let bj = pj.div(&rj, wp, nr).neg();
            Some((aj, bj))
        })
    })
}

/// Seed value and derivative of HeunC at a point inside the unit disk.
fn heun_c_jet(
    z: &ExactNum,
    alpha: &ExactNum,
    gamma: &ExactNum,
    delta: &ExactNum,
    eps: &ExactNum,
    q: &ExactNum,
    p: usize,
) -> Option<(ExactNum, ExactNum, bool, Vec<ExactNum>)> {
    let nr = RoundingMode::None;
    with_loss_jet(p, |wp| {
        let one = ExactNum::from_u8(1, wp);
        let l1 = ExactNum::sub_ref(delta, eps, wp, nr).sub(&ExactNum::from_u8(2, wp), wp, nr);
        let l0 = ExactNum::add_ref(q, gamma, wp, nr);
        let k0 = ExactNum::sub_ref(alpha, eps, wp, nr);
        let tail = RecTail {
            lam_d: ub64(&one),
            lam_s: ExactNum::new(64),
            lam_s_v: ExactNum::new(wp),
            k_c: ExactNum::new(64),
            l1: ub64(&l1),
            l0: ub64(&l0),
            k1: ub64(eps),
            k0: ub64(&k0),
            s: ExactNum::from_u8(1, 64),
            g: ub64(gamma),
        };
        let b1 = q.neg_via_mul().div(gamma, wp, nr);
        if !finite(&b1) {
            return None;
        }
        rec_series_jet(z, &b1, &tail, wp, |n, wp| {
            let nf = ExactNum::from_u64(n, wp);
            let nm1 = ExactNum::from_u64(n - 1, wp);
            let left =
                ExactNum::from_u64(n + 1, wp).mul(&ExactNum::add_ref(&nf, gamma, wp, nr), wp, nr);
            if left.is_zero() || !finite(&left) {
                return None;
            }
            let coef_n = ExactNum::add_ref(&nm1, gamma, wp, nr)
                .add(delta, wp, nr)
                .sub(eps, wp, nr)
                .mul(&nf, wp, nr)
                .sub(q, wp, nr);
            let coef_nm1 = ExactNum::mul_ref(&nm1, eps, wp, nr).add(alpha, wp, nr);
            Some((coef_n.div(&left, wp, nr), coef_nm1.div(&left, wp, nr)))
        })
    })
}

/// Taylor steps along one real segment. Each step is at most half the distance
/// to the nearest point of `sings` (the local radius of the ODE).
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

/// An upper bound for \(e^x\) when \(x\ge 0\), via argument
/// halving and a 24-term series with a geometric remainder.
fn exp_up_bound(x: &ExactNum) -> Option<ExactNum> {
    let up = RoundingMode::Up;
    if !finite(x) {
        return None;
    }
    if x.is_zero() || x.is_negative() {
        return Some(ExactNum::from_u8(1, 64));
    }
    let one = ExactNum::from_u8(1, 64);
    let mut t = ub64(x);
    let mut halves: u32 = 0;
    if let Some(ex) = t.exponent() {
        if ex >= 0 {
            let s = ex.checked_add(1)?;
            if s > 32 {
                // e^{2^{32}} overflows any 64-bit exponent bound.
                return None;
            }
            t = t.ldexp(-s, 64, up);
            halves = u32::try_from(s).ok()?;
        }
    }
    if matches!(t.compare(&one), Some(c) if c >= 0) {
        t = t.ldexp(-1, 64, up);
        halves = halves.checked_add(1)?;
    }
    if matches!(t.compare(&one), Some(c) if c >= 0) {
        return None;
    }
    // e^t = Σ_{k=0}^{24} t^k/k! + R with R ≤ (t^{24}/24!) t/(25-t) for t<1.
    let mut term = ExactNum::from_u8(1, 64);
    let mut sum = ExactNum::from_u8(1, 64);
    for k in 1..=24u32 {
        term = term.mul(&t, 64, up).div(&ExactNum::from_u32(k, 64), 64, up);
        sum = sum.add(&term, 64, up);
    }
    let den = ExactNum::from_u8(25, 64).sub(&t, 64, RoundingMode::Down);
    if !den.is_positive() {
        return None;
    }
    let rem = term.mul(&t, 64, up).div(&den, 64, up);
    let mut e = sum.add(&rem, 64, up);
    for _ in 0..halves {
        e = ExactNum::mul_ref(&e, &e, 64, up);
        if !finite(&e) {
            return None;
        }
    }
    Some(e)
}

/// Product \(\prod(\lvert z_c-s\rvert-\rho)\) (downward) — a lower bound for
/// \(\lvert A(\zeta)\rvert\) on \(\lvert\zeta\rvert\le\rho\) when the leading
/// coefficient vanishes only at `sings`.
fn amin_on_disk(zc: &ExactNum, sings: &[ExactNum], rho: &ExactNum) -> Option<ExactNum> {
    let dn = RoundingMode::Down;
    let mut prod = ExactNum::from_u8(1, 64);
    for s in sings {
        let mut d = ExactNum::sub_ref(zc, s, 64, RoundingMode::None).abs();
        let _ = d.set_precision(64, dn);
        let gap = ExactNum::sub_ref(&d, rho, 64, dn);
        if !gap.is_positive() || gap.is_zero() || !finite(&gap) {
            return None;
        }
        prod = prod.mul(&gap, 64, dn);
        if !prod.is_positive() || prod.is_zero() {
            return None;
        }
    }
    Some(prod)
}

/// Upper bounds for \(\lvert P\rvert=\lvert B/A\rvert\) and \(\lvert Q\rvert=\lvert C/A\rvert\)
/// on the disk \(\lvert\zeta\rvert\le\rho\).
fn pq_bounds(ode: &OdePoly, amin: &ExactNum, rho: &ExactNum) -> Option<(ExactNum, ExactNum)> {
    let up = RoundingMode::Up;
    let b = ub64(&ode.b0)
        .add(&ub64(&ode.b1).mul(rho, 64, up), 64, up)
        .add(&ub64(&ode.b2).mul(rho, 64, up).mul(rho, 64, up), 64, up);
    let c = ub64(&ode.c0).add(&ub64(&ode.c1).mul(rho, 64, up), 64, up);
    if amin.is_zero() || !amin.is_positive() {
        return None;
    }
    Some((b.div(amin, 64, up), c.div(amin, 64, up)))
}

/// Majorant \(M\ge\max_{\lvert\zeta\rvert=\rho}\lvert y(z_c+\zeta)\rvert\) for the
/// linear majorant \(u''=P_{\max}u'+Q_{\max}u\) with \(u(0)=\lvert y\rvert\),
/// \(u'(0)=\lvert y'\rvert\).
fn majorant_m(
    y: &ExactNum,
    yp: &ExactNum,
    pmax: &ExactNum,
    qmax: &ExactNum,
    rho: &ExactNum,
) -> Option<ExactNum> {
    let up = RoundingMode::Up;
    let y0 = ub_rel(y);
    let yp0 = ub_rel(yp);
    let pm = ub_rel(pmax);
    let qm = ub_rel(qmax);
    if qm.is_zero() && pm.is_zero() {
        return Some(ExactNum::add_ref(
            &y0,
            &ExactNum::mul_ref(&yp0, rho, 64, up),
            64,
            up,
        ));
    }
    if qm.is_zero() {
        // u(t)=y0 + yp0 (e^{P t}-1)/P ≤ y0 + yp0 e^{P t}/P
        let e = exp_up_bound(&ExactNum::mul_ref(&pm, rho, 64, up))?;
        let term = ExactNum::mul_ref(&yp0, &e, 64, up).div(&pm, 64, up);
        return Some(ExactNum::add_ref(&y0, &term, 64, up));
    }
    let disc = ExactNum::mul_ref(&pm, &pm, 64, up).add(&qm.ldexp(2, 64, up), 64, up);
    let sd = disc.sqrt(64, up);
    if !finite(&sd) {
        return None;
    }
    let two = ExactNum::from_u8(2, 64);
    let rp = ExactNum::add_ref(&pm, &sd, 64, up).div(&two, 64, up);
    let rm = ExactNum::sub_ref(&pm, &sd, 64, RoundingMode::Down).div(&two, 64, RoundingMode::Down);
    // r+ - r- = sqrt(disc)
    if sd.is_zero() || !sd.is_positive() {
        // repeated root r=P/2: u(t)=(y0+(yp0-r y0)t) e^{r t}
        let coef = ExactNum::add_ref(&yp0, &ub_rel(&ExactNum::mul_ref(&rp, &y0, 64, up)), 64, up);
        let e = exp_up_bound(&ExactNum::mul_ref(&rp, rho, 64, up))?;
        return Some(
            ExactNum::add_ref(&y0, &ExactNum::mul_ref(&coef, rho, 64, up), 64, up).mul(&e, 64, up),
        );
    }
    let a = ExactNum::add_ref(&yp0, &ExactNum::mul_ref(&ub_rel(&rm), &y0, 64, up), 64, up)
        .div(&sd, 64, up);
    let b = ub_rel(&ExactNum::sub_ref(
        &ExactNum::mul_ref(&rp, &y0, 64, up),
        &yp0,
        64,
        RoundingMode::None,
    ))
    .div(&sd, 64, up);
    let e = exp_up_bound(&ExactNum::mul_ref(&rp, rho, 64, up))?;
    Some(ExactNum::mul_ref(&a, &e, 64, up).add(&b, 64, up))
}

/// Midpoint + radius enclosure used for the outside-disk Taylor walk.
struct Br {
    mid: ExactNum,
    rad: ExactNum,
}

fn br_enclose(mid: ExactNum, rad: ExactNum) -> Br {
    Br {
        mid,
        rad: rad.abs(),
    }
}

fn ball_finish(mid: &ExactNum, rad: &ExactNum, p: usize, rm: RoundingMode) -> Option<ExactNum> {
    if !finite(mid) || !finite(rad) {
        return None;
    }
    let dest = crate::catalog::word_p(p);
    if rad.is_zero() {
        let v = mid.rounded_to(dest, rm);
        return Some(v);
    }
    let wide = dest.saturating_add(64).max(
        mid.mantissa_max_bit_len()
            .unwrap_or(dest)
            .saturating_add(64),
    );
    let lo = ExactNum::sub_ref(mid, rad, wide, RoundingMode::Down);
    let hi = ExactNum::add_ref(mid, rad, wide, RoundingMode::Up);
    if !finite(&lo) || !finite(&hi) {
        return None;
    }
    let mut rlo = lo;
    let mut rhi = hi;
    if rlo.set_precision(dest, rm).is_err() || rhi.set_precision(dest, rm).is_err() {
        return None;
    }
    if rlo.compare(&rhi) == Some(0) {
        Some(rlo)
    } else {
        None
    }
}

/// Redo the walk at higher working precision until the enclosure uniquely
/// determines the dest-p rounding. A zero radius finishes immediately.
/// If the relative radius stops shrinking, the result cannot be certified.
fn ball_ziv<F>(p: usize, rm: RoundingMode, mut walk: F) -> ExactNum
where
    F: FnMut(usize) -> Option<(ExactNum, ExactNum)>,
{
    let mut wp = work_p(p);
    let mut prev_rel: Option<Exponent> = None;
    loop {
        match walk(wp) {
            Some((mid, rad)) => {
                if let Some(v) = ball_finish(&mid, &rad, p, rm) {
                    return v;
                }
                let up = RoundingMode::Up;
                let den = ub_rel(&mid).add(&ExactNum::from_u8(1, 64), 64, up);
                let rel = ub_rel(&rad).div(&den, 64, up);
                let e = rel.exponent();
                if prev_rel.is_some_and(|pe| e.is_none_or(|e| e >= pe)) {
                    return cat_nan();
                }
                prev_rel = e;
            }
            None => {
                // A cut or a step that cannot be formed does not become
                // formable at a wider width after a failed attempt.
                return cat_nan();
            }
        }
        let next = wp.saturating_add(crate::defs::WORD_BIT_SIZE.max(32));
        if next <= wp {
            return cat_nan();
        }
        wp = next;
    }
}

/// Relative rounding unit at precision `p` (4 ulp plus a few guard bits).
fn unit_p(p: usize) -> ExactNum {
    ExactNum::from_u8(1, 64).ldexp(
        i64::try_from(p).map(|n| 16 - n).unwrap_or(0),
        64,
        RoundingMode::Up,
    )
}

/// Majorant for \(v''=P_{\max}v'+Q_{\max}v+\Delta\) with \(v(0)=v'(0)=0\).
/// If \(Q_{\max}>0\), \(w=v+\Delta/Q_{\max}\) solves the homogeneous majorant
/// with \(w(0)=\Delta/Q_{\max}\), \(w'(0)=0\), so \(v\le w\).
fn majorant_forced(
    pmax: &ExactNum,
    qmax: &ExactNum,
    delta: &ExactNum,
    rho: &ExactNum,
) -> Option<ExactNum> {
    let up = RoundingMode::Up;
    let d = ub_rel(delta);
    if d.is_zero() {
        return Some(ExactNum::new(64));
    }
    let pm = ub_rel(pmax);
    let qm = ub_rel(qmax);
    if qm.is_zero() && pm.is_zero() {
        return Some(d.mul(rho, 64, up).mul(rho, 64, up).ldexp(-1, 64, up));
    }
    if qm.is_zero() {
        let e = exp_up_bound(&ExactNum::mul_ref(&pm, rho, 64, up))?;
        return Some(d.mul(&e, 64, up).div(&pm, 64, up).div(&pm, 64, up));
    }
    let w0 = d.div(&qm, 64, up);
    majorant_m(&w0, &ExactNum::new(64), pmax, qmax, rho)
}

fn taylor_step(
    ode: &OdePoly,
    y: &ExactNum,
    yp: &ExactNum,
    h: &ExactNum,
    zc: &ExactNum,
    sings: &[ExactNum],
    p: usize,
    bits: usize,
) -> Option<(ExactNum, ExactNum, ExactNum, ExactNum)> {
    if h.is_zero() {
        return Some((y.clone(), yp.clone(), ExactNum::new(64), ExactNum::new(64)));
    }
    let up = RoundingMode::Up;
    let dn = RoundingMode::Down;
    let nr = RoundingMode::None;
    let mut ah = h.abs_ref();
    let _ = ah.set_precision(64, up);
    let mut radius = min_dist(zc, sings, 64);
    let _ = radius.set_precision(64, dn);
    if !radius.is_positive() || !finite(&radius) {
        return None;
    }
    let mut rho = ExactNum::add_ref(&ah, &radius, 64, up).ldexp(-1, 64, up);
    if !matches!(ah.compare(&rho), Some(c) if c < 0)
        || !matches!(rho.compare(&radius), Some(c) if c < 0)
    {
        rho = ExactNum::add_ref(&ah, &radius, 64, dn)
            .ldexp(-1, 64, dn)
            .add(&radius.ldexp_ref(-4, 64, dn), 64, up);
        if !matches!(ah.compare(&rho), Some(c) if c < 0)
            || !matches!(rho.compare(&radius), Some(c) if c < 0)
        {
            return None;
        }
    }
    let amin = amin_on_disk(zc, sings, &rho)?;
    let (pmax, qmax) = pq_bounds(ode, &amin, &rho)?;
    let m_bound = majorant_m(y, yp, &pmax, &qmax, &rho)?;
    if !finite(&m_bound) {
        return None;
    }
    let theta = ah.div(&rho, 64, up);
    if !matches!(theta.compare(&ExactNum::from_u8(1, 64)), Some(c) if c < 0) {
        return None;
    }
    let geo = ExactNum::from_u8(1, 64).sub(&theta, 64, dn);
    if !geo.is_positive() || geo.is_zero() {
        return None;
    }

    let prec = bits.max(p);
    let unit = unit_p(prec);
    let mut c_km1 = ExactNum::new(prec);
    let mut c_k = y.clone();
    let mut c_kp1 = yp.clone();
    let mut hpow = h.clone();
    let mut sum = ExactNum::add_ref(&ExactNum::mul_ref(yp, h, prec, nr), y, prec, nr);
    let mut dsum = yp.clone();
    let mut peak_y = Exponent::MIN;
    let mut peak_d = Exponent::MIN;
    bump_peak(&mut peak_y, &sum);
    bump_peak(&mut peak_d, &dsum);
    let mut tn = theta.clone();
    let mut k: u32 = 0;
    // A posteriori defect: Σ |residual_k| ρ^k of Aỹ''+Bỹ'+Cỹ, then /amin.
    let mut def = ExactNum::new(64);
    let mut rhok = ExactNum::from_u8(1, 64);
    loop {
        let ku = k;
        let kp1 = ku + 1;
        let kp2 = ku + 2;
        let ukp1 = un(kp1, prec);
        let ukp2 = un(kp2, prec);
        let uku = un(ku, prec);
        let left = ExactNum::mul_ref(
            &ExactNum::mul_ref(&ukp1, &ukp2, prec, nr),
            &ode.a0,
            prec,
            nr,
        );
        if left.is_zero() || !finite(&left) {
            return None;
        }
        let coef_kp1 = ExactNum::mul_ref(
            &ukp1,
            &ExactNum::add_ref(
                &ExactNum::mul_ref(&uku, &ode.a1, prec, nr),
                &ode.b0,
                prec,
                nr,
            ),
            p,
            nr,
        );
        let coef_k_owned;
        let (coef_k, coef_km1): (&ExactNum, ExactNum) = if k == 0 {
            (&ode.c0, ExactNum::new(prec))
        } else {
            let km1 = ku - 1;
            let km2 = ku.saturating_sub(2);
            let ukm1 = un(km1, prec);
            let ukm2 = un(km2, prec);
            let coef_k = ExactNum::add_ref(
                &ExactNum::add_ref(
                    &ExactNum::mul_ref(
                        &ExactNum::mul_ref(&uku, &ukm1, prec, nr),
                        &ode.a2,
                        prec,
                        nr,
                    ),
                    &ExactNum::mul_ref(&uku, &ode.b1, prec, nr),
                    p,
                    nr,
                ),
                &ode.c0,
                p,
                nr,
            );
            let coef_km1 = ExactNum::add_ref(
                &ExactNum::add_ref(
                    &ExactNum::mul_ref(
                        &ExactNum::mul_ref(&ukm1, &ukm2, prec, nr),
                        &ode.a3,
                        prec,
                        nr,
                    ),
                    &ExactNum::mul_ref(&ukm1, &ode.b2, prec, nr),
                    p,
                    nr,
                ),
                &ode.c1,
                p,
                nr,
            );
            coef_k_owned = coef_k;
            (&coef_k_owned, coef_km1)
        };
        let t_kp1 = ExactNum::mul_ref(&coef_kp1, &c_kp1, prec, nr);
        let t_k = ExactNum::mul_ref(coef_k, &c_k, prec, nr);
        let t_km1 = ExactNum::mul_ref(&coef_km1, &c_km1, prec, nr);
        let rhs = ExactNum::add_ref(&ExactNum::add_ref(&t_kp1, &t_k, prec, nr), &t_km1, prec, nr);
        let c_kp2 = ExactNum::div_ref(&rhs.neg_ref(), &left, prec, nr);
        if !finite(&c_kp2) {
            return None;
        }
        let rbound = ExactNum::mul_ref(
            &unit,
            &ub_rel(&rhs).add(&ub_rel(&left).mul(&ub_rel(&c_kp2), 64, up), 64, up),
            64,
            up,
        );
        def = def.add(&rbound.mul(&rhok, 64, up), 64, up);
        rhok = rhok.mul(&rho, 64, up);
        let h_prev = hpow.clone();
        hpow = ExactNum::mul_ref(&hpow, h, prec, nr);
        let term = ExactNum::mul_ref(&hpow, &c_kp2, prec, nr);
        let dterm = ExactNum::mul_ref(
            &ukp2,
            &ExactNum::mul_ref(&h_prev, &c_kp2, prec, nr),
            prec,
            nr,
        );
        if !finite(&term) || !finite(&dterm) {
            return None;
        }
        sum = ExactNum::add_ref(&sum, &term, prec, nr);
        dsum = ExactNum::add_ref(&dsum, &dterm, prec, nr);
        bump_peak(&mut peak_y, &term);
        bump_peak(&mut peak_y, &sum);
        bump_peak(&mut peak_d, &dterm);
        bump_peak(&mut peak_d, &dsum);
        c_km1 = c_k;
        c_k = c_kp1;
        c_kp1 = c_kp2;
        tn = tn.mul(&theta, 64, up);
        let n = (k as u64) + 2;
        let rem_y = ExactNum::mul_ref(&m_bound, &tn, 64, up)
            .mul(&theta, 64, up)
            .div(&geo, 64, up);
        let np1 = ExactNum::from_u64(n + 1, 64);
        let nn = ExactNum::from_u64(n, 64);
        let num = np1.sub(&nn.mul(&theta, 64, dn), 64, up);
        let rem_d = ExactNum::div_ref(&m_bound, &rho, 64, up)
            .mul(&tn, 64, up)
            .mul(&num, 64, up)
            .div(&ExactNum::mul_ref(&geo, &geo, 64, dn), 64, up);
        let y_ok = rem_y.is_zero() || rem_y.exponent().is_some_and(|e| e + (bits as i64) < peak_y);
        let d_ok = rem_d.is_zero() || rem_d.exponent().is_some_and(|e| e + (bits as i64) < peak_d);
        let term_ok = term.is_zero() || term.exponent().is_some_and(|e| e + (bits as i64) < peak_y);
        let dterm_ok =
            dterm.is_zero() || dterm.exponent().is_some_and(|e| e + (bits as i64) < peak_d);
        if (y_ok && d_ok && term_ok && dterm_ok && finite(&sum) && finite(&dsum)) || tn.is_zero() {
            let term_tail = ub_rel(&term).mul(&theta, 64, up).div(&geo, 64, up);
            let dterm_tail = ub_rel(&dterm).div(&geo, 64, up);
            let mut rem_y = if tn.is_zero() {
                ExactNum::new(64)
            } else if matches!(term_tail.compare(&rem_y), Some(c) if c > 0) {
                term_tail
            } else {
                rem_y
            };
            let mut rem_d = if tn.is_zero() {
                ExactNum::new(64)
            } else if matches!(dterm_tail.compare(&rem_d), Some(c) if c > 0) {
                dterm_tail
            } else {
                rem_d
            };
            // Perturbation of the working-precision OdePoly versus the exact
            // polynomial coefficients: |δA ỹ'' + δB ỹ' + δC ỹ| on the step disk.
            let da = ExactNum::mul_ref(
                &unit,
                &ub_rel(&ode.a0)
                    .add(&ub_rel(&ode.a1).mul(&rho, 64, up), 64, up)
                    .add(&ub_rel(&ode.a2).mul(&rho, 64, up).mul(&rho, 64, up), 64, up)
                    .add(
                        &ub_rel(&ode.a3)
                            .mul(&rho, 64, up)
                            .mul(&rho, 64, up)
                            .mul(&rho, 64, up),
                        64,
                        up,
                    ),
                64,
                up,
            );
            let db = ExactNum::mul_ref(
                &unit,
                &ub_rel(&ode.b0)
                    .add(&ub_rel(&ode.b1).mul(&rho, 64, up), 64, up)
                    .add(&ub_rel(&ode.b2).mul(&rho, 64, up).mul(&rho, 64, up), 64, up),
                64,
                up,
            );
            let dc = ExactNum::mul_ref(
                &unit,
                &ub_rel(&ode.c0).add(&ub_rel(&ode.c1).mul(&rho, 64, up), 64, up),
                64,
                up,
            );
            let den_h = ExactNum::mul_ref(&rho, &geo, 64, dn);
            let my = ub_rel(&m_bound);
            // Without a positive den_h both derivative bounds fall back to my itself.
            let (myp_div, mypp_div);
            let (myp, mypp): (&ExactNum, &ExactNum) = if den_h.is_positive() {
                myp_div = ExactNum::div_ref(&my, &den_h, 64, up);
                mypp_div =
                    my.ldexp_ref(1, 64, up)
                        .div(&ExactNum::mul_ref(&den_h, &den_h, 64, dn), 64, up);
                (&myp_div, &mypp_div)
            } else {
                (&my, &my)
            };
            def = def.add(&da.mul(mypp, 64, up), 64, up);
            def = def.add(&db.mul(myp, 64, up), 64, up);
            def = def.add(&dc.mul(&my, 64, up), 64, up);
            let nterm = ExactNum::from_u64((k as u64).saturating_add(4), 64);
            let sum_rnd = ExactNum::mul_ref(&unit, &nterm, 64, up).mul(
                &ub_rel(&sum).add(&ub_rel(&dsum), 64, up).add(&my, 64, up),
                64,
                up,
            );
            rem_y = rem_y.add(&sum_rnd, 64, up);
            rem_d = rem_d.add(&sum_rnd, 64, up);
            let delta = def.div(&amin, 64, up);
            if let Some(ey) = majorant_forced(&pmax, &qmax, &delta, &rho) {
                rem_y = rem_y.add(&ey, 64, up);
                let den = ExactNum::mul_ref(&rho, &geo, 64, dn);
                if den.is_positive() {
                    rem_d = rem_d.add(&ey.div(&den, 64, up), 64, up);
                }
            }
            return Some((sum, dsum, rem_y, rem_d));
        }
        k = k.checked_add(1)?;
    }
}

fn step_with_rad(
    ode: &OdePoly,
    y: &Br,
    yp: &Br,
    h: &ExactNum,
    zc: &ExactNum,
    sings: &[ExactNum],
    p: usize,
    acc: usize,
) -> Option<(Br, Br)> {
    let y_ub = ub_rel(&y.mid).add(&ub_rel(&y.rad), 64, RoundingMode::Up);
    let yp_ub = ub_rel(&yp.mid).add(&ub_rel(&yp.rad), 64, RoundingMode::Up);
    let (_, _, rem_y, rem_d) = taylor_step(ode, &y_ub, &yp_ub, h, zc, sings, p, acc)?;
    let (ym, ypm, _, _) = taylor_step(ode, &y.mid, &yp.mid, h, zc, sings, p, acc)?;
    let one = ExactNum::from_u8(1, p);
    let zero = ExactNum::new(p);
    let (y10, yp10, r10y, r10d) = taylor_step(ode, &one, &zero, h, zc, sings, p, acc)?;
    let (y01, yp01, r01y, r01d) = taylor_step(ode, &zero, &one, h, zc, sings, p, acc)?;
    let up = RoundingMode::Up;
    let ey = ub_rel(&y.rad);
    let eyp = ub_rel(&yp.rad);
    let phi11 = ub_rel(&y10).add(&ub_rel(&r10y), 64, up);
    let phi12 = ub_rel(&y01).add(&ub_rel(&r01y), 64, up);
    let phi21 = ub_rel(&yp10).add(&ub_rel(&r10d), 64, up);
    let phi22 = ub_rel(&yp01).add(&ub_rel(&r01d), 64, up);
    let peak = ub_rel(&ym)
        .add(&ub_rel(&ypm), 64, up)
        .add(&ExactNum::from_u8(1, 64), 64, up);
    let rnd = ExactNum::from_u8(1, 64)
        .ldexp(i64::try_from(p).map(|n| 8 - n).unwrap_or(0), 64, up)
        .mul(&peak, 64, up);
    let ey_new = phi11
        .mul(&ey, 64, up)
        .add(&ExactNum::mul_ref(&phi12, &eyp, 64, up), 64, up)
        .add(&ub_rel(&rem_y), 64, up)
        .add(&rnd, 64, up);
    let eyp_new = phi21
        .mul(&ey, 64, up)
        .add(&ExactNum::mul_ref(&phi22, &eyp, 64, up), 64, up)
        .add(&ub_rel(&rem_d), 64, up)
        .add(&rnd, 64, up);
    Some((br_enclose(ym, ey_new), br_enclose(ypm, eyp_new)))
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

fn err_from_loss(val: &ExactNum, loss: i64, wp: usize) -> ExactNum {
    let up = RoundingMode::Up;
    let bits = i64::try_from(wp).unwrap_or(i64::MAX).saturating_sub(loss);
    let body = ub_rel(val).ldexp(-bits, 64, up);
    let exp = i64::try_from(wp).map(|n| -n).unwrap_or(i64::MIN / 4);
    let unit = ExactNum::from_u8(1, 64).ldexp(exp, 64, up);
    body.add(&unit, 64, up)
}

fn taylor_walk_ball<F>(
    mut zc: ExactNum,
    mut y: Br,
    mut yp: Br,
    target: &ExactNum,
    sings: &[ExactNum],
    p: usize,
    mut ode_at: F,
) -> Option<(ExactNum, ExactNum)>
where
    F: FnMut(&ExactNum, usize) -> Option<OdePoly>,
{
    let acc = p.saturating_add(64);
    let up = RoundingMode::Up;
    loop {
        let rem = ExactNum::sub_ref(target, &zc, p, RoundingMode::None);
        if close_enough(&rem, target, p) || rem.is_zero() {
            if !finite(&y.mid) {
                return None;
            }
            // Include the last-step z-error: |Δz| (|yp| + rad_yp + 1).
            let dz = ub_rel(&rem.abs());
            let slope = ub_rel(&yp.mid).add(&ub_rel(&yp.rad), 64, up).add(
                &ExactNum::from_u8(1, 64),
                64,
                up,
            );
            let extra = ExactNum::mul_ref(&slope, &dz, 64, up);
            let rad = ub_rel(&y.rad).add(&extra, 64, up);
            return Some((y.mid, rad));
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
            if !finite(&y.mid) {
                return None;
            }
            return Some((y.mid, y.rad));
        }
        let ode = ode_at(&zc, acc)?;
        let (yn, ypn) = step_with_rad(&ode, &y, &yp, &h, &zc, sings, p, acc)?;
        let next = ExactNum::add_ref(&h, &zc, p, RoundingMode::None);
        if next.compare(&zc) == Some(0) {
            if !finite(&yn.mid) {
                return None;
            }
            return Some((yn.mid, yn.rad));
        }
        zc = next;
        y = yn;
        yp = ypn;
    }
}

fn disk_seed(z: &ExactNum, radius: ExactNum, p: usize) -> ExactNum {
    let half = radius.ldexp(-1, p, RoundingMode::None);
    if z.is_negative() {
        half.neg()
    } else {
        half
    }
}

fn heun_g_seed(z: &ExactNum, a: &ExactNum, wp: usize) -> ExactNum {
    let one = ExactNum::from_u8(1, wp);
    if matches!(a.abs_cmp(&one), Some(c) if c < 0) {
        let mut half = a.neg_via_mul().ldexp(-1, wp, RoundingMode::None);
        if half.is_negative() != z.is_negative() {
            half = -half;
        }
        half
    } else {
        disk_seed(z, ExactNum::from_u8(1, wp), wp)
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
    ball_ziv(p, rm, |wp| {
        let seed = heun_g_seed(z, a, wp);
        let (y, yp, poly, coefs) = heun_g_jet(&seed, a, q, alpha, beta, gamma, delta, wp)?;
        if poly {
            let v = horner(&coefs, z, wp);
            return Some((v, ExactNum::new(64)));
        }
        let one = ExactNum::from_u8(1, wp);
        if segment_hits(z, &one) || segment_hits(z, a) {
            return None;
        }
        let sings = [ExactNum::new(wp), one, a.clone()];
        let (ye, ype) = (err_from_loss(&y, 16, wp), err_from_loss(&yp, 16, wp));
        let yb = br_enclose(y, ye);
        let ypb = br_enclose(yp, ype);
        taylor_walk_ball(seed, yb, ypb, z, &sings, wp, |z0, pr| {
            heun_g_ode(z0, a, q, alpha, beta, gamma, delta, pr)
        })
    })
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
    if !z.is_negative() {
        return cat_nan();
    }
    ball_ziv(p, rm, |wp| {
        let seed = disk_seed(z, ExactNum::from_u8(1, wp), wp);
        let (y, yp, poly, coefs) = heun_c_jet(&seed, alpha, gamma, delta, eps, q, wp)?;
        if poly {
            let v = horner(&coefs, z, wp);
            return Some((v, ExactNum::new(64)));
        }
        let sings = [ExactNum::new(wp), ExactNum::from_u8(1, wp)];
        let (ye, ype) = (err_from_loss(&y, 16, wp), err_from_loss(&yp, 16, wp));
        let yb = br_enclose(y, ye);
        let ypb = br_enclose(yp, ype);
        taylor_walk_ball(seed, yb, ypb, z, &sings, wp, |z0, pr| {
            heun_c_ode(z0, alpha, gamma, delta, eps, q, pr)
        })
    })
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
    let mut a_frac: Vec<(u64, u64)> = Vec::new();
    for (v, s) in an.iter().chain(ap.iter()) {
        a.push(v.clone());
        a_frac.push(as_pos_rational(s, p)?);
    }
    let mut b: Vec<ExactNum> = Vec::new();
    let mut b_frac: Vec<(u64, u64)> = Vec::new();
    for (v, s) in bm.iter().chain(bq.iter()) {
        b.push(v.clone());
        b_frac.push(as_pos_rational(s, p)?);
    }
    let mut d: u64 = 1;
    for &(_, den) in a_frac.iter().chain(b_frac.iter()) {
        d = lcm_u(d, den)?;
    }
    let mut a_int: Vec<u64> = Vec::new();
    for &(num, den) in &a_frac {
        let ai = num.checked_mul(d / den)?;
        if ai == 0 {
            return None;
        }
        a_int.push(ai);
    }
    let mut b_int: Vec<u64> = Vec::new();
    for &(num, den) in &b_frac {
        let bi = num.checked_mul(d / den)?;
        if bi == 0 {
            return None;
        }
        b_int.push(bi);
    }
    let mut lift_len: usize = 0;
    for &ai in a_int.iter().chain(b_int.iter()) {
        let n = usize::try_from(ai).ok()?;
        lift_len = lift_len.checked_add(n)?;
    }
    // A dest-64 exact dyadic like 1085102592571150095/2^{63} fits in u64
    // but must not attempt a ~2^{63}-parameter lift: fall back to the
    // residue series when the allocation cannot be reserved.
    // Throwaway reserve: without observing the buffer, LLVM may legally elide
    // try_reserve and attempt the lift. black_box keeps the failure observable.
    let mut lift_probe: Vec<ExactNum> = Vec::new();
    if lift_probe.try_reserve(lift_len).is_err() {
        return None;
    }
    core::hint::black_box(&lift_probe);
    let mut a_tilde = Vec::new();
    for (ai, &ai_int) in a.iter().zip(a_int.iter()) {
        let ai_us = usize::try_from(ai_int).ok()?;
        let aif = ExactNum::from_u64(ai_int, p);
        for k in 0..ai_us {
            let num = ExactNum::from_u64(k as u64, p).add(ai, p, RoundingMode::None);
            a_tilde.push(num.div(&aif, p, RoundingMode::None));
        }
    }
    let mut b_tilde = Vec::new();
    for (bj, &bj_int) in b.iter().zip(b_int.iter()) {
        let bj_us = usize::try_from(bj_int).ok()?;
        let bif = ExactNum::from_u64(bj_int, p);
        for k in 0..bj_us {
            let num = ExactNum::from_u64(k as u64, p).add(bj, p, RoundingMode::None);
            b_tilde.push(num.div(&bif, p, RoundingMode::None));
        }
    }
    let n_tilde = a_int
        .iter()
        .take(n)
        .try_fold(0u64, |s, &x| s.checked_add(x))?;
    let m_tilde = b_int
        .iter()
        .take(m)
        .try_fold(0u64, |s, &x| s.checked_add(x))?;
    let n_t = usize::try_from(n_tilde).ok()?;
    let m_t = usize::try_from(m_tilde).ok()?;
    if n_t > a_tilde.len() || m_t > b_tilde.len() {
        return None;
    }
    let sum_u = |xs: &[u64]| xs.iter().try_fold(0i128, |s, &x| s.checked_add(x as i128));
    let a_star = sum_u(&a_int[..n])?
        .checked_sub(sum_u(&a_int[n..])?)?
        .checked_add(sum_u(&b_int[..m])?)?
        .checked_sub(sum_u(&b_int[m..])?)?;
    // c* = m + n - (p+q)/2  (original counts)
    let c_star_num: i128 = 2 * ((m + n) as i128) - (pn as i128) - (qn as i128);
    let one = ExactNum::from_u8(1, p);
    let two = ExactNum::from_u8(2, p);
    let half = ExactNum::from_u8(1, p).div(&two, p, RoundingMode::None);
    let mut beta = ExactNum::from_u8(1, p);
    for &ai in &a_int {
        let base = ExactNum::from_u64(ai, p);
        beta = beta.div(&pow_u64(&base, ai, p, cc), p, RoundingMode::None);
    }
    for &bj in &b_int {
        let base = ExactNum::from_u64(bj, p);
        beta = beta.mul(&pow_u64(&base, bj, p, cc), p, RoundingMode::None);
    }
    let mut m_fac = ExactNum::from_u8(1, p);
    for (bj, &bj_int) in b.iter().zip(b_int.iter()) {
        let base = ExactNum::from_u64(bj_int, p);
        let expn = ExactNum::sub_ref(bj, &half, p, RoundingMode::None);
        m_fac = m_fac.mul(&pow_real(&base, &expn, p, cc), p, RoundingMode::None);
    }
    for (ai, &ai_int) in a.iter().zip(a_int.iter()) {
        let base = ExactNum::from_u64(ai_int, p);
        let expn = ExactNum::from_u8(1, p)
            .div(&two, p, RoundingMode::None)
            .neg()
            .add(ai, p, RoundingMode::None);
        m_fac = m_fac.div(&pow_real(&base, &expn, p, cc), p, RoundingMode::None);
    }
    let pi = cc.pi(p, RoundingMode::None);
    let two_pi = ExactNum::from_u8(2, p).mul(&pi, p, RoundingMode::None);
    // (2π)^{c* - a*/2} with c* = c_star_num/2
    let expn = ExactNum::from_i128(c_star_num, p)
        .div(&two, p, RoundingMode::None)
        .sub(
            &ExactNum::from_i128(a_star, p).div(&two, p, RoundingMode::None),
            p,
            RoundingMode::None,
        );
    let pre = ExactNum::from_u64(d, p)
        .mul(&pow_real(&two_pi, &expn, p, cc), p, RoundingMode::None)
        .mul(&m_fac, p, RoundingMode::None);
    let r = one.div(&ExactNum::from_u64(d, p), p, RoundingMode::None);
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

fn all_scales_positive(
    an: &[(ExactNum, ExactNum)],
    ap: &[(ExactNum, ExactNum)],
    bm: &[(ExactNum, ExactNum)],
    bq: &[(ExactNum, ExactNum)],
) -> bool {
    an.iter()
        .chain(ap.iter())
        .chain(bm.iter())
        .chain(bq.iter())
        .all(|(_, s)| finite(s) && s.is_positive())
}

fn fox_choose_series(
    an: &[(ExactNum, ExactNum)],
    ap: &[(ExactNum, ExactNum)],
    bm: &[(ExactNum, ExactNum)],
    bq: &[(ExactNum, ExactNum)],
    z: &ExactNum,
    p: usize,
    cc: &mut Consts,
) -> u8 {
    let mut delta = ExactNum::new(p);
    for (_, s) in bm.iter().chain(bq.iter()) {
        delta = delta.add(s, p, RoundingMode::None);
    }
    for (_, s) in an.iter().chain(ap.iter()) {
        delta = delta.sub(s, p, RoundingMode::None);
    }
    let mut typ = if delta.is_zero() {
        let mut beta = ExactNum::from_u8(1, p);
        for (_, s) in an.iter().chain(ap.iter()) {
            beta = beta.mul(&pow_real(s, s, p, cc), p, RoundingMode::None);
        }
        for (_, s) in bm.iter().chain(bq.iter()) {
            let den = pow_real(s, s, p, cc);
            if !finite(&den) || den.is_zero() {
                return 1;
            }
            beta = beta.div(&den, p, RoundingMode::None);
        }
        if matches!(z.abs_cmp(&beta), Some(c) if c > 0) {
            2u8
        } else {
            1
        }
    } else if delta.is_positive() {
        1
    } else {
        2
    };
    if typ == 1 && bm.is_empty() {
        typ = 2;
    }
    if typ == 2 && an.is_empty() {
        typ = 1;
    }
    typ
}

fn lin_scale(c: &ExactNum, s: &ExactNum, slope: &ExactNum, p: usize) -> ExactNum {
    ExactNum::add_ref(
        c,
        &ExactNum::mul_ref(s, slope, p, RoundingMode::None),
        p,
        RoundingMode::None,
    )
}

fn fox_term_ratio(
    gn: &[ExactNum],
    gd: &[ExactNum],
    pref: &ExactNum,
    zp: &ExactNum,
    p: usize,
    cc: &mut Consts,
) -> Result<Option<ExactNum>, ()> {
    let ratio = match gamma_ratio(gn, gd, p, cc) {
        Err(()) => return Err(()),
        Ok(None) => return Ok(None),
        Ok(Some(v)) => v,
    };
    if !finite(pref) || !finite(zp) || !finite(&ratio) {
        return Err(());
    }
    Ok(Some(
        ExactNum::mul_ref(pref, &ratio, p, RoundingMode::None).mul(zp, p, RoundingMode::None),
    ))
}

/// Mellin–Barnes residue series at the exact scales. No term cap: the sum
/// stops when a term is below a unit in the last place of `p`, or stalls if
/// the terms stop shrinking. A `u64` term index overflow is an index-type
/// limit.
fn fox_h_residue(
    z: &ExactNum,
    an: &[(ExactNum, ExactNum)],
    ap: &[(ExactNum, ExactNum)],
    bm: &[(ExactNum, ExactNum)],
    bq: &[(ExactNum, ExactNum)],
    p: usize,
    cc: &mut Consts,
) -> GEval {
    if !finite(z)
        || !pairs_finite(an)
        || !pairs_finite(ap)
        || !pairs_finite(bm)
        || !pairs_finite(bq)
        || !all_scales_positive(an, ap, bm, bq)
    {
        return GEval::Bad;
    }
    if an.is_empty() && bm.is_empty() {
        return GEval::Bad;
    }
    if !z.is_positive() && !z.is_zero() {
        return GEval::Bad;
    }
    let typ = fox_choose_series(an, ap, bm, bq, z, p, cc);
    let one = ExactNum::from_u8(1, p);
    let mut sum = ExactNum::new(p);
    let cutoff = -((p as i64).saturating_add(16));
    if typ == 1 {
        for h in 0..bm.len() {
            let (bh, bh_s) = &bm[h];
            if !bh_s.is_positive() {
                return GEval::Bad;
            }
            let inv_b = ExactNum::div_ref(&one, bh_s, p, RoundingMode::None);
            if !finite(&inv_b) {
                return GEval::Bad;
            }
            // s = -(b_h + k)/B_h, z^{-s} = z^{b_h/B_h} (z^{1/B_h})^k
            let mut zp = pow_real(
                z,
                &ExactNum::mul_ref(bh, &inv_b, p, RoundingMode::None),
                p,
                cc,
            );
            let z_step = pow_real(z, &inv_b, p, cc);
            if !finite(&zp) || !finite(&z_step) {
                return GEval::Bad;
            }
            let mut pref = ExactNum::div_ref(&one, bh_s, p, RoundingMode::None);
            let mut k: u64 = 0;
            let mut grew = 0u32;
            let mut prev_exp: Option<i64> = None;
            loop {
                let s = ExactNum::add_ref(bh, &ExactNum::from_u64(k, p), p, RoundingMode::None)
                    .neg()
                    .div(bh_s, p, RoundingMode::None);
                let mut gn = Vec::new();
                let mut gd = Vec::new();
                for (j, (bj, bj_s)) in bm.iter().enumerate() {
                    if j == h {
                        continue;
                    }
                    gn.push(lin_scale(bj, bj_s, &s, p));
                }
                for (ai, ai_s) in an {
                    gn.push(lin_scale(
                        &ExactNum::sub_ref(&one, ai, p, RoundingMode::None),
                        &ai_s.neg_ref(),
                        &s,
                        p,
                    ));
                }
                for (ai, ai_s) in ap {
                    gd.push(lin_scale(ai, ai_s, &s, p));
                }
                for (bj, bj_s) in bq {
                    gd.push(lin_scale(
                        &ExactNum::sub_ref(&one, bj, p, RoundingMode::None),
                        &bj_s.neg_ref(),
                        &s,
                        p,
                    ));
                }
                match fox_term_ratio(&gn, &gd, &pref, &zp, p, cc) {
                    Err(()) => return GEval::Singular,
                    Ok(None) => {}
                    Ok(Some(term)) => {
                        if !finite(&term) {
                            return GEval::Singular;
                        }
                        sum = sum.add(&term, p, RoundingMode::None);
                        if !finite(&sum) {
                            return GEval::Singular;
                        }
                        let te = term.abs().exponent();
                        if k >= 4 {
                            if let Some(e) = te {
                                let se = sum.abs_ref().exponent().unwrap_or(0);
                                if e - se < cutoff {
                                    break;
                                }
                                if prev_exp.is_some_and(|pe| e >= pe) {
                                    grew += 1;
                                    if grew >= 4 {
                                        return GEval::Bad;
                                    }
                                } else {
                                    grew = 0;
                                }
                            }
                        }
                        prev_exp = te;
                    }
                }
                let k_next = match k.checked_add(1) {
                    Some(v) => v,
                    None => return GEval::Bad,
                };
                pref = pref
                    .neg()
                    .div(&ExactNum::from_u64(k_next, p), p, RoundingMode::None);
                zp = zp.mul(&z_step, p, RoundingMode::None);
                k = k_next;
            }
        }
    } else {
        for k_idx in 0..an.len() {
            let (ak, ak_s) = &an[k_idx];
            if !ak_s.is_positive() {
                return GEval::Bad;
            }
            let inv_a = ExactNum::div_ref(&one, ak_s, p, RoundingMode::None);
            if !finite(&inv_a) {
                return GEval::Bad;
            }
            // s = (1 - a_k + ell)/A_k,  residue of Γ(1-a-A s) is -(-1)^ell/(ell! A)
            // Closing to the right, the MB integral picks up + [(-1)^ell/(ell! A)] * rest * z^{-s}
            // (the minus from dw = -A ds cancels the clockwise sign). z^{-s} starts
            // at z^{-(1-a_k)/A_k} and multiplies by z^{-1/A_k} each step.
            let s0 =
                ExactNum::sub_ref(&one, ak, p, RoundingMode::None).div(ak_s, p, RoundingMode::None);
            let mut zp = pow_real(z, &s0.neg(), p, cc);
            let z_step = pow_real(z, &inv_a.neg(), p, cc);
            if !finite(&zp) || !finite(&z_step) {
                return GEval::Bad;
            }
            let mut pref = ExactNum::div_ref(&one, ak_s, p, RoundingMode::None);
            let mut ell: u64 = 0;
            let mut grew = 0u32;
            let mut prev_exp: Option<i64> = None;
            loop {
                let s = ExactNum::add_ref(
                    &ExactNum::sub_ref(&one, ak, p, RoundingMode::None),
                    &ExactNum::from_u64(ell, p),
                    p,
                    RoundingMode::None,
                )
                .div(ak_s, p, RoundingMode::None);
                let mut gn = Vec::new();
                let mut gd = Vec::new();
                for (j, (ai, ai_s)) in an.iter().enumerate() {
                    if j == k_idx {
                        continue;
                    }
                    gn.push(lin_scale(
                        &ExactNum::sub_ref(&one, ai, p, RoundingMode::None),
                        &ai_s.neg_ref(),
                        &s,
                        p,
                    ));
                }
                for (bj, bj_s) in bm {
                    gn.push(lin_scale(bj, bj_s, &s, p));
                }
                for (bj, bj_s) in bq {
                    gd.push(lin_scale(
                        &ExactNum::sub_ref(&one, bj, p, RoundingMode::None),
                        &bj_s.neg_ref(),
                        &s,
                        p,
                    ));
                }
                for (ai, ai_s) in ap {
                    gd.push(lin_scale(ai, ai_s, &s, p));
                }
                match fox_term_ratio(&gn, &gd, &pref, &zp, p, cc) {
                    Err(()) => return GEval::Singular,
                    Ok(None) => {}
                    Ok(Some(term)) => {
                        if !finite(&term) {
                            return GEval::Singular;
                        }
                        sum = sum.add(&term, p, RoundingMode::None);
                        if !finite(&sum) {
                            return GEval::Singular;
                        }
                        let te = term.abs().exponent();
                        if ell >= 4 {
                            if let Some(e) = te {
                                let se = sum.abs_ref().exponent().unwrap_or(0);
                                if e - se < cutoff {
                                    break;
                                }
                                if prev_exp.is_some_and(|pe| e >= pe) {
                                    grew += 1;
                                    if grew >= 4 {
                                        return GEval::Bad;
                                    }
                                } else {
                                    grew = 0;
                                }
                            }
                        }
                        prev_exp = te;
                    }
                }
                let ell_next = match ell.checked_add(1) {
                    Some(v) => v,
                    None => return GEval::Bad,
                };
                pref = pref
                    .neg()
                    .div(&ExactNum::from_u64(ell_next, p), p, RoundingMode::None);
                zp = zp.mul(&z_step, p, RoundingMode::None);
                ell = ell_next;
            }
        }
    }
    GEval::Val(sum)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_must::Must;

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
        let hc = zc.heun_c(
            &al,
            &gam,
            &zero,
            &ExactNum::from_i8(-1, p),
            &al,
            p,
            rm,
            &mut cc,
        );
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

    #[test]
    fn caps_heun_outside_disk() -> Result<(), String> {
        // Independent golds: high-dps seed series + Taylor walk (mpmath, 200
        // decimals). Rows are real-analytic continuations along the component
        // of R\{1,a} that contains 0. All 7 rounding modes, correctly rounded
        // (0 ulps) at p=64 and p=256.
        let mut cc = Consts::new()?;
        let modes = [
            RoundingMode::ToEven,
            RoundingMode::ToOdd,
            RoundingMode::Up,
            RoundingMode::Down,
            RoundingMode::ToZero,
            RoundingMode::FromZero,
            RoundingMode::None,
        ];
        // kind, z, params..., gold
        // G: a,q,alpha,beta,gamma,delta; C: alpha,gamma,delta,eps,q
        type HeunGRow = (
            &'static str,
            &'static str,
            &'static str,
            &'static str,
            &'static str,
            &'static str,
            &'static str,
            &'static str,
        );
        let rows_g: &[HeunGRow] = &[
            (
                "-2",
                "2",
                "1/2",
                "1/3",
                "1/4",
                "1/2",
                "1/5",
                "0.489531139762585335975564627741117594291665579108339509717535235713547166753031876516108640578662697762686563085023431212",
            ),
            (
                "-10",
                "2",
                "1/2",
                "1/3",
                "1/4",
                "1/2",
                "1/5",
                "0.0722608484365945565377386205725345690771598706286165174292745548539479051178652876978056268162763678158612251427065379989",
            ),
            (
                "-1",
                "2",
                "1/2",
                "1/3",
                "1/4",
                "1/2",
                "1/5",
                "0.665318669703791779022287903523602300823525459068884726095333766816505264741601700649182096720874293746452627804307966501",
            ),
            (
                "-100",
                "2",
                "1/2",
                "1/3",
                "1/4",
                "1/2",
                "1/5",
                "-0.213629016605290085507329077652964168387620909199819320666683535541853961748965422270344358336047111974201683020057385713",
            ),
            (
                "3/4",
                "-1/2",
                "1/3",
                "1/3",
                "1/4",
                "4/5",
                "1/5",
                "0.541109480341224921497926456318611001263606675253815533927604013369321946441540084169633512336262543380447760423949067410",
            ),
            (
                "-2",
                "2",
                "1/2",
                "20",
                "20",
                "1/2",
                "1/5",
                "-0.000000566423072895740360266964394373159162709548961058477664999987667492511361009333570878610634962981025672750947742101730383",
            ),
        ];
        let rows_c: &[(&str, &str, &str, &str, &str, &str, &str)] = &[
            (
                "-2",
                "1/3",
                "1/2",
                "1/4",
                "1/5",
                "1/2",
                "3.59125837941479058712814911724782167694979476701259856800909417467855931362641846636120541112206018450666418784947132537",
            ),
            (
                "-10",
                "1/3",
                "1/2",
                "1/4",
                "1/5",
                "1/2",
                "51.8835611416647192527759604682647386371500331599599686107185103964073371484675603452348090453289141671103273769271419770",
            ),
            (
                "-1",
                "1/3",
                "1/2",
                "1/4",
                "1/5",
                "1/2",
                "2.12433591747382292256155153347925090856883937067446492672334734879491714563922844198420459825297571527356755835158096473",
            ),
        ];
        fn dec(s: &str, wp: usize, cc: &mut Consts) -> ExactNum {
            ExactNum::parse(s, crate::Radix::Dec, wp, RoundingMode::ToEven, cc)
        }
        fn rat(s: &str, wp: usize, rm: RoundingMode) -> Result<ExactNum, String> {
            if let Some((ns, ds)) = s.split_once('/') {
                let n: i32 = ns.parse().must("parse")?;
                let d: i32 = ds.parse().must("parse")?;
                let mut x = ExactNum::from_u32(n.unsigned_abs(), wp);
                if n < 0 {
                    x = x.neg();
                }
                if d == 1 {
                    Ok(x)
                } else {
                    Ok(x.div(&ExactNum::from_u32(d as u32, wp), wp, rm))
                }
            } else {
                let n: i32 = s.parse().must("parse")?;
                let mut x = ExactNum::from_u32(n.unsigned_abs(), wp);
                if n < 0 {
                    x = x.neg();
                }
                Ok(x)
            }
        }
        // Parameters are built at a width far above dest p. A dest-p rounding
        // of a non-dyadic parameter (1/3, 1/5, …) is a different exact Heun
        // problem; on a long walk that 2^{-p} perturbation can move the dest-p
        // value by one ulp relative to the exact-rational gold. z in these
        // rows is a small dyadic and is built at dest p. The walk enclosure is
        // for the values actually passed.
        const PARAM_WP: usize = 1024;
        let mut bad = Vec::new();
        for p in [64usize, 256] {
            for rm in modes {
                for (z, a, q, al, be, ga, de, gold) in rows_g {
                    let w = dec(gold, 8960, &mut cc);
                    let v = rat(z, p, rm)?.heun_g(
                        &rat(a, PARAM_WP, RoundingMode::ToEven)?,
                        &rat(q, PARAM_WP, RoundingMode::ToEven)?,
                        &rat(al, PARAM_WP, RoundingMode::ToEven)?,
                        &rat(be, PARAM_WP, RoundingMode::ToEven)?,
                        &rat(ga, PARAM_WP, RoundingMode::ToEven)?,
                        &rat(de, PARAM_WP, RoundingMode::ToEven)?,
                        p,
                        rm,
                        &mut cc,
                    );
                    let mut want = w;
                    let _ = want.set_precision(p, rm);
                    if v.compare(&want) != Some(0) {
                        bad.push(format!(
                            "heun_g({z}; a={a}, q={q}, α={al}, β={be}, γ={ga}, δ={de}) p={p} {rm:?} got {v} want {want}"
                        ));
                    }
                }
                for (z, al, ga, de, eps, q, gold) in rows_c {
                    let w = dec(gold, 8960, &mut cc);
                    let v = rat(z, p, rm)?.heun_c(
                        &rat(al, PARAM_WP, RoundingMode::ToEven)?,
                        &rat(ga, PARAM_WP, RoundingMode::ToEven)?,
                        &rat(de, PARAM_WP, RoundingMode::ToEven)?,
                        &rat(eps, PARAM_WP, RoundingMode::ToEven)?,
                        &rat(q, PARAM_WP, RoundingMode::ToEven)?,
                        p,
                        rm,
                        &mut cc,
                    );
                    let mut want = w;
                    let _ = want.set_precision(p, rm);
                    if v.compare(&want) != Some(0) {
                        bad.push(format!(
                            "heun_c({z}; α={al}, γ={ga}, δ={de}, ε={eps}, q={q}) p={p} {rm:?} got {v} want {want}"
                        ));
                    }
                }
            }
        }
        assert!(bad.is_empty(), "{} wrong:\n{}", bad.len(), bad.join("\n"));
        Ok(())
    }

    /// Public `meijer_g` used to refuse more than 8 parameters in `a` or in `b`.
    /// `G^{9,0}_{0,9}(1/2 | ; 0^{9})` is 9 `b`-parameters; 0 ulp vs mpmath 95 dps.
    #[test]
    fn caps_meijer_param_count_above_8() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let modes = [
            RoundingMode::ToEven,
            RoundingMode::ToOdd,
            RoundingMode::Up,
            RoundingMode::Down,
            RoundingMode::ToZero,
            RoundingMode::FromZero,
            RoundingMode::None,
        ];
        // mpmath 95 dps, agrees with 90 dps to 93 digits.
        let gold_s = "1.1747085653955017064263103776019618650648148969496353924236083981194706914930405850797897626190e-1";
        let gold = ExactNum::parse(
            gold_s,
            crate::Radix::Dec,
            1024,
            RoundingMode::ToEven,
            &mut cc,
        );
        let mut bad = Vec::new();
        for &p in &[64usize, 256] {
            for &rm in &modes {
                let z =
                    ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, RoundingMode::None);
                let bm: Vec<ExactNum> = (0..9).map(|_| ExactNum::new(p)).collect();
                let got = z.meijer_g(&[], &[], &bm, &[], p, rm, &mut cc);
                let want = crate::catalog::finish_ref(&gold, p, rm);
                if got.compare(&want) != Some(0) {
                    bad.push(format!(
                        "G^{{9,0}}_{{0,9}}(1/2) p={p} {rm:?} got {got} want {want}"
                    ));
                }
            }
        }
        assert!(bad.is_empty(), "{} wrong:\n{}", bad.len(), bad.join("\n"));
        Ok(())
    }

    /// `G^{16,0}_{0,16}(1/2 | ; 0^{16})` has 15 counted gamma poles, so
    /// `hi = 4*15 + 4 = 64` used to refuse (`hi > 63`).
    #[test]
    fn caps_meijer_pole_order_above_14() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let modes = [
            RoundingMode::ToEven,
            RoundingMode::ToOdd,
            RoundingMode::Up,
            RoundingMode::Down,
            RoundingMode::ToZero,
            RoundingMode::FromZero,
            RoundingMode::None,
        ];
        // mpmath 95 dps, agrees with 90 dps to 94 digits.
        let gold_s = "3.8637998573893274293136670187237271567980851122233596666791557290041543018069486541266016024764e-2";
        let gold = ExactNum::parse(
            gold_s,
            crate::Radix::Dec,
            1024,
            RoundingMode::ToEven,
            &mut cc,
        );
        let mut bad = Vec::new();
        // One jet at dest width 256 has enough bits for both dest widths.
        let z = ExactNum::from_u8(1, 256).div(&ExactNum::from_u8(2, 256), 256, RoundingMode::None);
        let r = ExactNum::from_u8(1, 256);
        let bm: Vec<ExactNum> = (0..16).map(|_| ExactNum::new(256)).collect();
        let g = meijer_jet_limit(&z, &r, false, &[], &[], &bm, &[], 256, &mut cc);
        for &p in &[64usize, 256] {
            for &rm in &modes {
                let got = crate::catalog::finish_ref(&g, p, rm);
                let want = crate::catalog::finish_ref(&gold, p, rm);
                if got.compare(&want) != Some(0) {
                    bad.push(format!(
                        "G^{{16,0}}_{{0,16}}(1/2) p={p} {rm:?} got {got} want {want}"
                    ));
                }
            }
        }
        assert!(bad.is_empty(), "{} wrong:\n{}", bad.len(), bad.join("\n"));
        Ok(())
    }

    /// Public `fox_h` at dest `p` (not work_p): hybrid of Jeff's 2026-11-17
    /// 8:33 AM CST ruling, "Fox H: hybrid, lift exact fractions and use the
    /// residue series for everything else." Inputs are built at the public
    /// dest precision. Golds are independent and evaluated at the exact
    /// dest-p dyadic: \(H^{1,0}_{0,1}\) is the closed form
    /// \((1/B) z^{b/B} \exp(-z^{1/B})\); \(H^{2,0}_{0,2}\) is the two-family
    /// Mellin–Barnes residue series (mpmath, 160 dps). 0 ulp, 7 modes,
    /// p=64 and 256. Includes dest-64/256 2/17 and 1/3 (must not be NaN),
    /// the dest-64 `RoundingMode::None` 2/17 dyadic (must not attempt a
    /// \(2^{63}\)-parameter lift), and the \(2/17+2^{-200}\) trap at dest 256
    /// (must be H of that dyadic, not H(2/17)).
    #[test]
    fn caps_fox_h_rational_scale_caps() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let modes = [
            RoundingMode::ToEven,
            RoundingMode::ToOdd,
            RoundingMode::Up,
            RoundingMode::Down,
            RoundingMode::ToZero,
            RoundingMode::FromZero,
            RoundingMode::None,
        ];
        fn rat(n: u32, d: u32, p: usize, rm: RoundingMode) -> ExactNum {
            if d == 1 {
                ExactNum::from_u32(n, p)
            } else {
                ExactNum::from_u32(n, p).div(&ExactNum::from_u32(d, p), p, rm)
            }
        }
        fn parse_gold(s: &str, cc: &mut Consts) -> ExactNum {
            ExactNum::parse(s, crate::Radix::Dec, 1024, RoundingMode::ToEven, cc)
        }
        // Independent golds at the exact dest-p input (see /tmp/fox_h_gold.py).
        // name, p, input rm, bm as (b_n,b_d,B_n,B_d) plus optional trap flag, gold
        type FoxRow<'a> = (
            &'a str,
            usize,
            RoundingMode,
            &'a [(u32, u32, u32, u32)],
            bool,
            &'a str,
        );
        let rows: &[FoxRow] = &[
            (
                "H^{1,0}_{0,1}(1/2|(0,2/17)) dest-64 ToEven",
                64,
                RoundingMode::ToEven,
                &[(0, 1, 2, 17)],
                false,
                "8.4765542402492593615599565141120465009571342701261058065558072752501918028265430126284747993886851565936525932",
            ),
            (
                "H^{1,0}_{0,1}(1/2|(0,2/17)) dest-64 None",
                64,
                RoundingMode::None,
                &[(0, 1, 2, 17)],
                false,
                "8.4765542402492593620561365637176756071266455746431292165652507766414790195273454744691924617201878164795150088",
            ),
            (
                "H^{1,0}_{0,1}(1/2|(0,2/17)) dest-256 ToEven",
                256,
                RoundingMode::ToEven,
                &[(0, 1, 2, 17)],
                false,
                "8.4765542402492593615891435758535540954361606206614689367285702748218220092725227772318248601952848693342698784",
            ),
            (
                "H^{1,0}_{0,1}(1/2|(0,2/17)) dest-256 None",
                256,
                RoundingMode::None,
                &[(0, 1, 2, 17)],
                false,
                "8.4765542402492593615891435758535540954361606206614689367285702748218220092726018232741740628599301144177862950",
            ),
            (
                "H^{1,0}_{0,1}(1/2|(0,1/3)) dest-64 ToEven",
                64,
                RoundingMode::ToEven,
                &[(0, 1, 1, 3)],
                false,
                "2.6474907077537862085042633610098426968418977610544813471244904302515056854204447016735489644196780419382603567",
            ),
            (
                "H^{1,0}_{0,1}(1/2|(0,1/3)) dest-64 None",
                64,
                RoundingMode::None,
                &[(0, 1, 1, 3)],
                false,
                "2.6474907077537862087755025670417712323292255651526123029758758310686775405528845871331853462801175748655703590",
            ),
            (
                "H^{1,0}_{0,1}(1/2|(0,1/3)) dest-256 ToEven",
                256,
                RoundingMode::ToEven,
                &[(0, 1, 1, 3)],
                false,
                "2.6474907077537862085946764296871522086660144749719522253109279576243278856339117866652283935228949345180793524",
            ),
            (
                "H^{1,0}_{0,1}(1/2|(0,1/3)) dest-256 None",
                256,
                RoundingMode::None,
                &[(0, 1, 1, 3)],
                false,
                "2.6474907077537862085946764296871522086660144749719522253109279576243278856339549975637395389237010787717529303",
            ),
            (
                "H^{1,0}_{0,1}(1/2|(0,13)) exact integer",
                64,
                RoundingMode::ToEven,
                &[(0, 1, 13, 1)],
                false,
                "0.029806557128670610131633200033113069561999100194288326523783881499742055990273122848946489408177306626787358014",
            ),
            (
                "H^{1,0}_{0,1}(1/2|(0,13)) exact integer dest-256",
                256,
                RoundingMode::ToEven,
                &[(0, 1, 13, 1)],
                false,
                "0.029806557128670610131633200033113069561999100194288326523783881499742055990273122848946489408177306626787358014",
            ),
            (
                "H^{2,0}_{0,2}(1/2|(0,9),(1/2,9)) exact dest-64",
                64,
                RoundingMode::ToEven,
                &[(0, 1, 9, 1), (1, 2, 9, 1)],
                false,
                "0.028744543869538835816376707165363370283649562396805557031716529579208402034590765308567248688101310048939276853",
            ),
            (
                "H^{2,0}_{0,2}(1/2|(0,9),(1/2,9)) exact dest-256",
                256,
                RoundingMode::ToEven,
                &[(0, 1, 9, 1), (1, 2, 9, 1)],
                false,
                "0.028744543869538835816376707165363370283649562396805557031716529579208402034590765308567248688101310048939276853",
            ),
            (
                "H^{2,0}_{0,2}(1/2|(0,1/7),(1/3,1/8)) dest-64 ToEven",
                64,
                RoundingMode::ToEven,
                &[(0, 1, 1, 7), (1, 3, 1, 8)],
                false,
                "14.150861693701019943312221271756651817302063884367624916798395324859134509485433566869870699535918786680547826",
            ),
            (
                "H^{2,0}_{0,2}(1/2|(0,1/7),(1/3,1/8)) dest-256 ToEven",
                256,
                RoundingMode::ToEven,
                &[(0, 1, 1, 7), (1, 3, 1, 8)],
                false,
                "14.150861693701019943330055886306502206671506422108512889717862092656741324465824420587589908295954680378089541",
            ),
            (
                "H^{1,0}_{0,1}(1/2|(0,2/17+2^{-200})) dest-256 trap",
                256,
                RoundingMode::ToEven,
                &[(0, 1, 2, 17)],
                true,
                "8.4765542402492593615891435758535540954361606206614689367285247078807827796319787846512311293365979764113378273",
            ),
        ];
        let mut bad = Vec::new();
        for (name, p, in_rm, bms, trap, gold_s) in rows {
            let gold = parse_gold(gold_s, &mut cc);
            let z = ExactNum::from_u8(1, *p).ldexp(-1, *p, RoundingMode::None);
            let bm: Vec<(ExactNum, ExactNum)> = bms
                .iter()
                .map(|&(bn, bd, sn, sd)| {
                    let b = rat(bn, bd, *p, *in_rm);
                    let mut s = rat(sn, sd, *p, *in_rm);
                    if *trap {
                        s = s.add(
                            &ExactNum::from_u8(1, *p).ldexp(-200, *p, RoundingMode::None),
                            *p,
                            RoundingMode::None,
                        );
                    }
                    (b, s)
                })
                .collect();
            // No slack identification of the trap as 2/17.
            if *trap {
                if let Some((n, d)) = exact_pos_frac(&bm[0].1) {
                    if (n, d) == (2, 17) {
                        bad.push(format!("{name}: exact_pos_frac identified 2/17"));
                    }
                }
            }
            for &rm in &modes {
                let got = z.fox_h(&[], &[], &bm, &[], *p, rm, &mut cc);
                if got.is_nan() {
                    bad.push(format!("{name} p={p} {rm:?} got NaN"));
                    continue;
                }
                let want = crate::catalog::finish_ref(&gold, *p, rm);
                if got.compare(&want) != Some(0) {
                    bad.push(format!("{name} p={p} {rm:?} got {got} want {want}"));
                }
            }
        }
        assert!(bad.is_empty(), "{} wrong:\n{}", bad.len(), bad.join("\n"));
        Ok(())
    }

    fn is_ia(x: &ExactNum) -> bool {
        x.is_nan() && x.err() == Some(crate::Error::InvalidArgument)
    }

    /// \(G^{m,n}_{p,q}\) is undefined when some \(a_k-b_j\) (\(k\le n\), \(j\le m\))
    /// is a positive integer (DLMF 16.17). In particular
    /// \(G^{1,1}_{1,1}(z\mid 1;0)=\Gamma(0)\,z^0(1+z)^{0}\) diverges for every \(z\),
    /// not only \(z=0\). The Hadamard finite part is a regularization, not the value.
    #[test]
    fn zf_meijer_fox_ak_minus_bj_positive_integer_is_invalid() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let modes = [
            RoundingMode::Up,
            RoundingMode::Down,
            RoundingMode::ToZero,
            RoundingMode::FromZero,
            RoundingMode::ToEven,
            RoundingMode::ToOdd,
            RoundingMode::None,
        ];
        for p in [64usize, 256, 1024, 2048] {
            let zero = ExactNum::new(p);
            let one = ExactNum::from_u8(1, p);
            let two = ExactNum::from_u8(2, p);
            let half = ExactNum::from_u8(1, p).div(&two, p, RoundingMode::ToEven);
            let an1 = [ExactNum::from_u8(1, p)];
            let bm0 = [ExactNum::new(p)];
            let zs = [
                ExactNum::new(p),
                ExactNum::from_u8(1, p).div(&two, p, RoundingMode::ToEven),
                ExactNum::from_u8(1, p),
                ExactNum::from_u8(2, p),
            ];
            for z in &zs {
                for rm in modes {
                    let g = z.meijer_g(&an1, &[], &bm0, &[], p, rm, &mut cc);
                    assert!(is_ia(&g), "G(1;0) z={z} p={p} {rm:?} got {g}");
                    let h = z.fox_h(
                        &[(ExactNum::from_u8(1, p), ExactNum::from_u8(1, p))],
                        &[],
                        &[(ExactNum::new(p), ExactNum::from_u8(1, p))],
                        &[],
                        p,
                        rm,
                        &mut cc,
                    );
                    assert!(is_ia(&h), "H(1,1;0,1) z={z} p={p} {rm:?} got {h}");
                }
            }
            if p > 256 {
                continue;
            }
            let rm = RoundingMode::ToEven;
            // a-b = 2 is also a positive integer.
            let g2 = zero.meijer_g(core::slice::from_ref(&two), &[], &bm0, &[], p, rm, &mut cc);
            assert!(is_ia(&g2), "G(2;0) at 0 p={p} got {g2}");
            let g2z = one.meijer_g(core::slice::from_ref(&two), &[], &bm0, &[], p, rm, &mut cc);
            assert!(is_ia(&g2z), "G(2;0) at 1 p={p} got {g2z}");
            // Sibling with no an-bm positive integer stays finite: G^{1,0}_{0,1}(0|;0)=1.
            let sib = zero.meijer_g(&[], &[], core::slice::from_ref(&zero), &[], p, rm, &mut cc);
            assert!(
                sib.compare(&one) == Some(0),
                "sibling G^{{1,0}}_{{0,1}}(0) p={p} got {sib}"
            );
            // a-b = 1/2 is not an integer; the value is finite at z=1/2.
            let ghalf = half.meijer_g(core::slice::from_ref(&half), &[], &bm0, &[], p, rm, &mut cc);
            assert!(
                !ghalf.is_nan() && finite(&ghalf),
                "G(1/2;0) at 1/2 should be finite p={p} got {ghalf}"
            );
        }
        Ok(())
    }
}
