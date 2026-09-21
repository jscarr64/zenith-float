//! Complex elliptic integrals via Carlson \(R_F,R_C,R_D,R_J\).
//!
//! Duplication runs in \(\mathbb{C}\) with the principal square root. Parameter
//! \(m=k^2\); incomplete argument \(x=\sin\varphi\). Not the real series at
//! \(\lvert z\rvert\).

use crate::common::util::round_p;
use crate::complex_special::nan_pair;
use crate::complex_special::pi_c;
use crate::complex_special::ziv_complex;
use crate::Consts;
use crate::Error;
use crate::ExactComplex;
use crate::ExactNum;
use crate::RoundingMode;
use crate::INF_POS;

/// Same cap as the real Carlson kernels.
const CARLSON_DUPE_MAX: u32 = 128;

fn rm() -> RoundingMode {
    RoundingMode::None
}

fn c_u32(n: u32, p: usize) -> ExactComplex {
    ExactComplex::from_real(ExactNum::from_u32(n, p), p)
}

fn is_c_zero(z: &ExactComplex) -> bool {
    z.re().is_zero() && z.im().is_zero()
}

fn is_c_one(z: &ExactComplex, p: usize) -> bool {
    z.im().is_zero() && z.re().cmp(&ExactNum::from_u8(1, p)) == Some(0)
}

fn c_abs(z: &ExactComplex, p: usize) -> ExactNum {
    z.abs(p, rm())
}

fn max_abs(a: &ExactNum, b: &ExactNum, _p: usize) -> ExactNum {
    if matches!(a.cmp(b), Some(c) if c >= 0) {
        a.clone()
    } else {
        b.clone()
    }
}

fn tiny_spread(p: usize) -> ExactNum {
    ExactNum::from_u8(1, p).ldexp(-((p as i32) / 3 + 16), p, rm())
}

fn close_enough(dev: &ExactNum, an: &ExactComplex, p: usize) -> bool {
    let one = ExactNum::from_u8(1, p);
    let scale = max_abs(&c_abs(an, p), &one, p);
    let thresh = tiny_spread(p).mul(&scale, p, rm());
    matches!(dev.cmp(&thresh), Some(c) if c < 0)
}

fn max_dev3(
    an: &ExactComplex,
    x: &ExactComplex,
    y: &ExactComplex,
    z: &ExactComplex,
    p: usize,
) -> ExactNum {
    let dx = c_abs(&an.sub(x, p, rm()), p);
    let dy = c_abs(&an.sub(y, p, rm()), p);
    let dz = c_abs(&an.sub(z, p, rm()), p);
    max_abs(&max_abs(&dx, &dy, p), &dz, p)
}

fn any_nan(args: &[&ExactComplex]) -> bool {
    args.iter().any(|z| z.is_nan())
}

/// \(R_C(x,y)=R_F(x,y,y)\).
fn carlson_rc(x: &ExactComplex, y: &ExactComplex, p: usize, cc: &mut Consts) -> ExactComplex {
    if any_nan(&[x, y]) {
        return ExactComplex::new(x.re().clone(), y.im().clone());
    }
    if is_c_zero(y) {
        return nan_pair(Error::InvalidArgument);
    }
    let dxy = x.sub(y, p, rm());
    if close_enough(&c_abs(&dxy, p), x, p) {
        return c_u32(1, p).div(&x.sqrt(p, rm(), cc), p, rm());
    }
    if is_c_zero(x) {
        let mut hp = pi_c(p, cc);
        hp = hp.div(&c_u32(2, p), p, rm());
        return hp.div(&y.sqrt(p, rm(), cc), p, rm());
    }
    carlson_rf(x, y, y, p, cc)
}

/// Symmetric \(R_F(x,y,z)\). At most one argument may be \(0\).
fn carlson_rf(
    x0: &ExactComplex,
    y0: &ExactComplex,
    z0: &ExactComplex,
    p: usize,
    cc: &mut Consts,
) -> ExactComplex {
    if any_nan(&[x0, y0, z0]) {
        return nan_pair(Error::InvalidArgument);
    }
    let zeros =
        usize::from(is_c_zero(x0)) + usize::from(is_c_zero(y0)) + usize::from(is_c_zero(z0));
    if zeros > 1 {
        return nan_pair(Error::InvalidArgument);
    }
    let four = c_u32(4, p);
    let three = c_u32(3, p);
    let mut x = x0.clone();
    let mut y = y0.clone();
    let mut z = z0.clone();
    for _ in 0..CARLSON_DUPE_MAX {
        let an = x.add(&y, p, rm()).add(&z, p, rm()).div(&three, p, rm());
        if close_enough(&max_dev3(&an, &x, &y, &z, p), &an, p) {
            return rf_series(&an, &x, &y, &z, p, cc);
        }
        let sx = x.sqrt(p, rm(), cc);
        let sy = y.sqrt(p, rm(), cc);
        let sz = z.sqrt(p, rm(), cc);
        let lam = sx
            .mul(&sy, p, rm())
            .add(&sy.mul(&sz, p, rm()), p, rm())
            .add(&sz.mul(&sx, p, rm()), p, rm());
        x = x.add(&lam, p, rm()).div(&four, p, rm());
        y = y.add(&lam, p, rm()).div(&four, p, rm());
        z = z.add(&lam, p, rm()).div(&four, p, rm());
    }
    nan_pair(Error::InvalidArgument)
}

fn rf_series(
    an: &ExactComplex,
    x: &ExactComplex,
    y: &ExactComplex,
    z: &ExactComplex,
    p: usize,
    cc: &mut Consts,
) -> ExactComplex {
    let xx = an.sub(x, p, rm()).div(an, p, rm());
    let yy = an.sub(y, p, rm()).div(an, p, rm());
    let zz = an.sub(z, p, rm()).div(an, p, rm());
    let e2 = xx.mul(&yy, p, rm()).sub(&zz.mul(&zz, p, rm()), p, rm());
    let e3 = xx.mul(&yy, p, rm()).mul(&zz, p, rm());
    let e2s = e2.mul(&e2, p, rm());
    let one = c_u32(1, p);
    let w = |n: u32| c_u32(n, p);
    let s = one
        .sub(&e2.div(&w(10), p, rm()), p, rm())
        .add(&e3.div(&w(14), p, rm()), p, rm())
        .add(&e2s.div(&w(24), p, rm()), p, rm())
        .sub(
            &w(3)
                .mul(&e2, p, rm())
                .mul(&e3, p, rm())
                .div(&w(44), p, rm()),
            p,
            rm(),
        )
        .sub(
            &w(5)
                .mul(&e2, p, rm())
                .mul(&e2s, p, rm())
                .div(&w(208), p, rm()),
            p,
            rm(),
        )
        .add(
            &w(3)
                .mul(&e2s, p, rm())
                .mul(&e3, p, rm())
                .div(&w(104), p, rm()),
            p,
            rm(),
        );
    s.div(&an.sqrt(p, rm(), cc), p, rm())
}

/// \(R_D(x,y,z)\). \(z\neq 0\); at most one of \(x,y\) may be \(0\).
fn carlson_rd(
    x0: &ExactComplex,
    y0: &ExactComplex,
    z0: &ExactComplex,
    p: usize,
    cc: &mut Consts,
) -> ExactComplex {
    if any_nan(&[x0, y0, z0]) {
        return nan_pair(Error::InvalidArgument);
    }
    if is_c_zero(z0) || (is_c_zero(x0) && is_c_zero(y0)) {
        return nan_pair(Error::InvalidArgument);
    }
    let four = c_u32(4, p);
    let three = c_u32(3, p);
    let five = c_u32(5, p);
    let mut x = x0.clone();
    let mut y = y0.clone();
    let mut z = z0.clone();
    let mut sum = ExactComplex::zero(p);
    let mut fac = c_u32(1, p);
    for _ in 0..CARLSON_DUPE_MAX {
        let an = x
            .add(&y, p, rm())
            .add(&three.mul(&z, p, rm()), p, rm())
            .div(&five, p, rm());
        if close_enough(&max_dev3(&an, &x, &y, &z, p), &an, p) {
            let series = rd_series(&an, &x, &y, &z, p, cc);
            return three
                .mul(&sum, p, rm())
                .add(&fac.mul(&series, p, rm()), p, rm());
        }
        let sx = x.sqrt(p, rm(), cc);
        let sy = y.sqrt(p, rm(), cc);
        let sz = z.sqrt(p, rm(), cc);
        let lam = sx
            .mul(&sy, p, rm())
            .add(&sy.mul(&sz, p, rm()), p, rm())
            .add(&sz.mul(&sx, p, rm()), p, rm());
        sum = sum.add(
            &fac.div(&sz.mul(&z.add(&lam, p, rm()), p, rm()), p, rm()),
            p,
            rm(),
        );
        fac = fac.div(&four, p, rm());
        x = x.add(&lam, p, rm()).div(&four, p, rm());
        y = y.add(&lam, p, rm()).div(&four, p, rm());
        z = z.add(&lam, p, rm()).div(&four, p, rm());
    }
    nan_pair(Error::InvalidArgument)
}

fn rd_series(
    an: &ExactComplex,
    x: &ExactComplex,
    y: &ExactComplex,
    z: &ExactComplex,
    p: usize,
    cc: &mut Consts,
) -> ExactComplex {
    let xx = an.sub(x, p, rm()).div(an, p, rm());
    let yy = an.sub(y, p, rm()).div(an, p, rm());
    let zz = an.sub(z, p, rm()).div(an, p, rm());
    let e2 = xx.mul(&yy, p, rm()).sub(&zz.mul(&zz, p, rm()), p, rm());
    let e3 = xx.mul(&yy, p, rm()).mul(&zz, p, rm());
    let e2s = e2.mul(&e2, p, rm());
    let one = c_u32(1, p);
    let w = |n: u32| c_u32(n, p);
    let s = one
        .sub(&w(3).mul(&e2, p, rm()).div(&w(14), p, rm()), p, rm())
        .add(&e3.div(&w(6), p, rm()), p, rm())
        .add(&w(9).mul(&e2s, p, rm()).div(&w(88), p, rm()), p, rm())
        .sub(
            &w(3)
                .mul(&e2, p, rm())
                .mul(&e3, p, rm())
                .div(&w(22), p, rm()),
            p,
            rm(),
        )
        .add(
            &w(9)
                .mul(&e3, p, rm())
                .mul(&e3, p, rm())
                .div(&w(52), p, rm()),
            p,
            rm(),
        )
        .sub(
            &w(3)
                .mul(&e2, p, rm())
                .mul(&e2s, p, rm())
                .div(&w(26), p, rm()),
            p,
            rm(),
        );
    s.div(&an.mul(&an.sqrt(p, rm(), cc), p, rm()), p, rm())
}

/// \(R_J(x,y,z,p)\). Characteristic \(p\neq 0\); at most one of \(x,y,z\) may be \(0\).
fn carlson_rj(
    x0: &ExactComplex,
    y0: &ExactComplex,
    z0: &ExactComplex,
    p0: &ExactComplex,
    p: usize,
    cc: &mut Consts,
) -> ExactComplex {
    if any_nan(&[x0, y0, z0, p0]) {
        return nan_pair(Error::InvalidArgument);
    }
    if is_c_zero(p0) {
        return nan_pair(Error::InvalidArgument);
    }
    let zeros =
        usize::from(is_c_zero(x0)) + usize::from(is_c_zero(y0)) + usize::from(is_c_zero(z0));
    if zeros > 1 {
        return nan_pair(Error::InvalidArgument);
    }
    let two = c_u32(2, p);
    let three = c_u32(3, p);
    let four = c_u32(4, p);
    let five = c_u32(5, p);
    let mut x = x0.clone();
    let mut y = y0.clone();
    let mut z = z0.clone();
    let mut pv = p0.clone();
    let mut sum = ExactComplex::zero(p);
    let mut fac = c_u32(1, p);
    for _ in 0..CARLSON_DUPE_MAX {
        let an = x
            .add(&y, p, rm())
            .add(&z, p, rm())
            .add(&two.mul(&pv, p, rm()), p, rm())
            .div(&five, p, rm());
        let d4 = max_dev3(&an, &x, &y, &z, p);
        let dp = c_abs(&an.sub(&pv, p, rm()), p);
        let dmax = max_abs(&d4, &dp, p);
        if close_enough(&dmax, &an, p) {
            let series = rj_series(&an, &x, &y, &z, &pv, p, cc);
            return three
                .mul(&sum, p, rm())
                .add(&fac.mul(&series, p, rm()), p, rm());
        }
        let sx = x.sqrt(p, rm(), cc);
        let sy = y.sqrt(p, rm(), cc);
        let sz = z.sqrt(p, rm(), cc);
        let lam = sx
            .mul(&sy, p, rm())
            .add(&sy.mul(&sz, p, rm()), p, rm())
            .add(&sz.mul(&sx, p, rm()), p, rm());
        let alpha = pv
            .mul(&sx.add(&sy, p, rm()).add(&sz, p, rm()), p, rm())
            .add(&sx.mul(&sy, p, rm()).mul(&sz, p, rm()), p, rm());
        let alpha = alpha.mul(&alpha, p, rm());
        let pl = pv.add(&lam, p, rm());
        let beta = pv.mul(&pl, p, rm()).mul(&pl, p, rm());
        sum = sum.add(
            &fac.mul(&carlson_rc(&alpha, &beta, p, cc), p, rm()),
            p,
            rm(),
        );
        fac = fac.div(&four, p, rm());
        x = x.add(&lam, p, rm()).div(&four, p, rm());
        y = y.add(&lam, p, rm()).div(&four, p, rm());
        z = z.add(&lam, p, rm()).div(&four, p, rm());
        pv = pv.add(&lam, p, rm()).div(&four, p, rm());
    }
    nan_pair(Error::InvalidArgument)
}

fn rj_series(
    an: &ExactComplex,
    x: &ExactComplex,
    y: &ExactComplex,
    z: &ExactComplex,
    pv: &ExactComplex,
    p: usize,
    cc: &mut Consts,
) -> ExactComplex {
    let xx = an.sub(x, p, rm()).div(an, p, rm());
    let yy = an.sub(y, p, rm()).div(an, p, rm());
    let zz = an.sub(z, p, rm()).div(an, p, rm());
    let pp = an.sub(pv, p, rm()).div(an, p, rm());
    let xyz = xx.mul(&yy, p, rm()).mul(&zz, p, rm());
    let xy_xz_yz = xx
        .mul(&yy, p, rm())
        .add(&xx.mul(&zz, p, rm()), p, rm())
        .add(&yy.mul(&zz, p, rm()), p, rm());
    let p2 = pp.mul(&pp, p, rm());
    let p3 = p2.mul(&pp, p, rm());
    let two = c_u32(2, p);
    let three = c_u32(3, p);
    let e2 = xy_xz_yz.sub(&three.mul(&p2, p, rm()), p, rm());
    let e3 = xyz
        .add(&two.mul(&p3, p, rm()), p, rm())
        .sub(&pp.mul(&xy_xz_yz, p, rm()), p, rm());
    let e2s = e2.mul(&e2, p, rm());
    let one = c_u32(1, p);
    let w = |n: u32| c_u32(n, p);
    let s = one
        .sub(&w(3).mul(&e2, p, rm()).div(&w(14), p, rm()), p, rm())
        .add(&e3.div(&w(6), p, rm()), p, rm())
        .add(&w(9).mul(&e2s, p, rm()).div(&w(88), p, rm()), p, rm())
        .sub(
            &w(3)
                .mul(&e2, p, rm())
                .mul(&e3, p, rm())
                .div(&w(22), p, rm()),
            p,
            rm(),
        )
        .add(
            &w(9)
                .mul(&e3, p, rm())
                .mul(&e3, p, rm())
                .div(&w(52), p, rm()),
            p,
            rm(),
        )
        .sub(
            &w(3)
                .mul(&e2, p, rm())
                .mul(&e2s, p, rm())
                .div(&w(26), p, rm()),
            p,
            rm(),
        );
    s.div(&an.mul(&an.sqrt(p, rm(), cc), p, rm()), p, rm())
}

impl ExactComplex {
    /// Complete elliptic \(K(m)\), \(m=k^2\). Cut on \([1,+\infty)\). \(m=1\) is \(+\infty\).
    ///
    /// # Precision
    ///
    /// - Algorithm: Carlson `R_F` in ℂ; `CARLSON_DUPE_MAX = 128`.
    /// - Bound: identities evaluated outside Ziv (nested Ziv would exhaust `MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    pub fn elliptic_k(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return ExactComplex::new(self.re().clone(), self.im().clone());
        }
        let dest = round_p(p);
        if is_c_one(self, dest) {
            return ExactComplex::from_real(INF_POS.clone(), dest);
        }
        if is_c_zero(self) {
            let hp = pi_c(dest, cc);
            return hp.div(&c_u32(2, dest), dest, rm);
        }
        ziv_complex(dest, rm, |pw| self.elliptic_k_at(pw, cc))
    }

    fn elliptic_k_at(&self, p: usize, cc: &mut Consts) -> Self {
        let zero = ExactComplex::zero(p);
        let one = ExactComplex::one(p);
        let om = one.sub(self, p, rm());
        carlson_rf(&zero, &om, &one, p, cc)
    }

    /// Complete elliptic \(E(m)\). \(E(1)=1\). Cut of \(K\) inherited through \(1-m\).
    ///
    /// # Precision
    ///
    /// - Algorithm: Carlson `R_F` / `R_D`; `CARLSON_DUPE_MAX = 128`.
    /// - Bound: same as [`Self::elliptic_k`].
    /// - MPFR oracle: no.
    pub fn elliptic_e_complete(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return ExactComplex::new(self.re().clone(), self.im().clone());
        }
        let dest = round_p(p);
        if is_c_zero(self) {
            let hp = pi_c(dest, cc);
            return hp.div(&c_u32(2, dest), dest, rm);
        }
        if is_c_one(self, dest) {
            return ExactComplex::one(dest);
        }
        ziv_complex(dest, rm, |pw| self.elliptic_e_complete_at(pw, cc))
    }

    fn elliptic_e_complete_at(&self, p: usize, cc: &mut Consts) -> Self {
        let zero = ExactComplex::zero(p);
        let one = ExactComplex::one(p);
        let three = c_u32(3, p);
        let om = one.sub(self, p, rm());
        let rf = carlson_rf(&zero, &om, &one, p, cc);
        let rd = carlson_rd(&zero, &om, &one, p, cc);
        rf.sub(&self.div(&three, p, rm()).mul(&rd, p, rm()), p, rm())
    }

    /// Incomplete \(F(x|m)\), \(x=\sin\varphi\), \(m=k^2\).
    ///
    /// # Precision
    ///
    /// - Algorithm: Carlson `R_F`; `CARLSON_DUPE_MAX = 128`.
    /// - MPFR oracle: no.
    ///
    /// Carlson \(R_F(1-x^2,1-mx^2,1)\). Cuts when \(1-x^2\) or \(1-mx^2\) lies on
    /// \((-\infty,0]\) (principal square-root cut). \(F(x,0)=\arcsin x\).
    pub fn elliptic_f(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() || m.is_nan() {
            return nan_pair(Error::InvalidArgument);
        }
        let dest = round_p(p);
        if is_c_zero(self) {
            return ExactComplex::zero(dest);
        }
        if is_c_zero(m) {
            return self.asin(dest, rm, cc);
        }
        if is_c_one(m, dest) {
            return self.atanh(dest, rm, cc);
        }
        ziv_complex(dest, rm, |pw| self.elliptic_f_at(m, pw, cc))
    }

    fn elliptic_f_at(&self, m: &Self, p: usize, cc: &mut Consts) -> Self {
        let one = ExactComplex::one(p);
        let x2 = self.mul(self, p, rm());
        let a = one.sub(&x2, p, rm());
        let b = one.sub(&m.mul(&x2, p, rm()), p, rm());
        let rf = carlson_rf(&a, &b, &one, p, cc);
        self.mul(&rf, p, rm())
    }

    /// Incomplete \(E(x|m)\). Same \(x,m\) convention as [`Self::elliptic_f`].
    ///
    /// # Precision
    ///
    /// - Algorithm: Carlson `R_F` / `R_D`; `CARLSON_DUPE_MAX = 128`.
    /// - MPFR oracle: no.
    ///
    /// Cuts as for \(F\). \(E(x,0)=\arcsin x\); \(E(x,1)=x\).
    pub fn elliptic_e(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() || m.is_nan() {
            return nan_pair(Error::InvalidArgument);
        }
        let dest = round_p(p);
        if is_c_zero(self) {
            return ExactComplex::zero(dest);
        }
        if is_c_zero(m) {
            return self.asin(dest, rm, cc);
        }
        if is_c_one(m, dest) {
            return self.clone();
        }
        ziv_complex(dest, rm, |pw| self.elliptic_e_at(m, pw, cc))
    }

    fn elliptic_e_at(&self, m: &Self, p: usize, cc: &mut Consts) -> Self {
        let one = ExactComplex::one(p);
        let three = c_u32(3, p);
        let x2 = self.mul(self, p, rm());
        let a = one.sub(&x2, p, rm());
        let b = one.sub(&m.mul(&x2, p, rm()), p, rm());
        let rf = carlson_rf(&a, &b, &one, p, cc);
        let rd = carlson_rd(&a, &b, &one, p, cc);
        self.mul(&rf, p, rm()).sub(
            &m.mul(self, p, rm())
                .mul(&x2, p, rm())
                .div(&three, p, rm())
                .mul(&rd, p, rm()),
            p,
            rm(),
        )
    }

    /// Complete \(\Pi(n,m)\). `self` is \(n\). \(\Pi(0,m)=K(m)\). Pole at \(n=1\).
    ///
    /// # Precision
    ///
    /// - Algorithm: Carlson `R_J`; `CARLSON_DUPE_MAX = 128`.
    /// - MPFR oracle: no.
    pub fn elliptic_pi_complete(
        &self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        if self.is_nan() || m.is_nan() {
            return nan_pair(Error::InvalidArgument);
        }
        let dest = round_p(p);
        if is_c_one(self, dest) {
            return nan_pair(Error::InvalidArgument);
        }
        if is_c_zero(self) {
            return m.elliptic_k(dest, rm, cc);
        }
        ziv_complex(dest, rm, |pw| self.elliptic_pi_complete_at(m, pw, cc))
    }

    fn elliptic_pi_complete_at(&self, m: &Self, p: usize, cc: &mut Consts) -> Self {
        let zero = ExactComplex::zero(p);
        let one = ExactComplex::one(p);
        let three = c_u32(3, p);
        let om = one.sub(m, p, rm());
        let on = one.sub(self, p, rm());
        let rf = carlson_rf(&zero, &om, &one, p, cc);
        let rj = carlson_rj(&zero, &om, &one, &on, p, cc);
        rf.add(&self.div(&three, p, rm()).mul(&rj, p, rm()), p, rm())
    }

    /// Incomplete \(\Pi(n;x|m)\). `self` is \(n\).
    ///
    /// # Precision
    ///
    /// - Algorithm: Carlson `R_J`; `CARLSON_DUPE_MAX = 128`.
    /// - MPFR oracle: no.
    ///
    /// Cuts when \(1-x^2\), \(1-mx^2\), or \(1-nx^2\) meets the Carlson cut
    /// \((-\infty,0]\). \(\Pi(0;x|m)=F(x|m)\).
    pub fn elliptic_pi(
        &self,
        x: &Self,
        m: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        if self.is_nan() || x.is_nan() || m.is_nan() {
            return nan_pair(Error::InvalidArgument);
        }
        let dest = round_p(p);
        if is_c_zero(x) {
            return ExactComplex::zero(dest);
        }
        if is_c_zero(self) {
            return x.elliptic_f(m, dest, rm, cc);
        }
        ziv_complex(dest, rm, |pw| self.elliptic_pi_at(x, m, pw, cc))
    }

    fn elliptic_pi_at(&self, x: &Self, m: &Self, p: usize, cc: &mut Consts) -> Self {
        let one = ExactComplex::one(p);
        let three = c_u32(3, p);
        let x2 = x.mul(x, p, rm());
        let a = one.sub(&x2, p, rm());
        let b = one.sub(&m.mul(&x2, p, rm()), p, rm());
        let pv = one.sub(&self.mul(&x2, p, rm()), p, rm());
        let rf = carlson_rf(&a, &b, &one, p, cc);
        let rj = carlson_rj(&a, &b, &one, &pv, p, cc);
        x.mul(&rf, p, rm()).add(
            &self
                .mul(x, p, rm())
                .mul(&x2, p, rm())
                .div(&three, p, rm())
                .mul(&rj, p, rm()),
            p,
            rm(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(a: &ExactNum, b: &ExactNum, p: usize, slack: i32) -> bool {
        let d = a.sub(b, p, RoundingMode::None).abs();
        d.is_zero() || d.exponent().unwrap_or(0) < -((p as i32) - slack)
    }

    fn tiny(x: &ExactNum, p: usize) -> bool {
        x.is_zero() || x.exponent().unwrap_or(0) < -((p as i32) / 4)
    }

    #[test]
    fn complex_elliptic_golds() {
        let p = 256;
        let r = RoundingMode::ToEven;
        let mut cc = Consts::new().unwrap();
        let z0 = ExactComplex::zero(p);
        let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, r);
        let mh = ExactComplex::from_real(half.clone(), p);

        let k0 = z0.elliptic_k(p, r, &mut cc);
        let mut hp = cc.pi(p, r);
        hp = hp.div(&ExactNum::from_u8(2, p), p, r);
        assert!(near(k0.re(), &hp, p, 40), "K(0) re");
        assert!(tiny(k0.im(), p), "K(0) im");

        let e0 = z0.elliptic_e_complete(p, r, &mut cc);
        assert!(near(e0.re(), &hp, p, 40), "E(0) re");
        assert!(tiny(e0.im(), p), "E(0) im");

        let k_real = half.elliptic_k(p, r, &mut cc);
        let k_c = mh.elliptic_k(p, r, &mut cc);
        assert!(near(k_c.re(), &k_real, p, 40), "K(1/2) matches real");
        assert!(tiny(k_c.im(), p), "K(1/2) im");

        let e1 = ExactComplex::one(p).elliptic_e_complete(p, r, &mut cc);
        assert!(near(e1.re(), &ExactNum::from_u8(1, p), p, 40), "E(1)");
        assert!(tiny(e1.im(), p));

        let k1 = ExactComplex::one(p).elliptic_k(p, r, &mut cc);
        assert!(k1.re().is_inf_pos(), "K(1)=+∞");

        let x = mh.clone();
        let f0 = x.elliptic_f(&z0, p, r, &mut cc);
        let asin = x.asin(p, r, &mut cc);
        assert!(near(f0.re(), asin.re(), p, 8), "F(x,0)=arcsin");
        assert!(tiny(f0.im(), p) && tiny(asin.im(), p));

        let n0 = z0.elliptic_pi_complete(&mh, p, r, &mut cc);
        assert!(near(n0.re(), k_c.re(), p, 40), "Π(0,m)=K(m)");

        // Legendre: E(m)K(1-m)+E(1-m)K(m)-K(m)K(1-m)=π/2
        let m = ExactComplex::new(
            ExactNum::from_u8(3, p).div(&ExactNum::from_u8(10, p), p, r),
            ExactNum::from_u8(1, p).div(&ExactNum::from_u8(10, p), p, r),
        );
        let mp = ExactComplex::one(p).sub(&m, p, r);
        let km = m.elliptic_k(p, r, &mut cc);
        let kp = mp.elliptic_k(p, r, &mut cc);
        let em = m.elliptic_e_complete(p, r, &mut cc);
        let ep = mp.elliptic_e_complete(p, r, &mut cc);
        let lhs = em
            .mul(&kp, p, r)
            .add(&ep.mul(&km, p, r), p, r)
            .sub(&km.mul(&kp, p, r), p, r);
        let want = ExactComplex::from_real(hp, p);
        assert!(
            near(lhs.re(), want.re(), p, 30),
            "Legendre re {:?}",
            lhs.re()
        );
        assert!(
            tiny(lhs.im(), p) || near(lhs.im(), want.im(), p, 30),
            "Legendre im"
        );

        // Cut of K on [1,+∞): K(2+εi) and K(2-εi) are conjugates, not equal.
        let two = ExactNum::from_u8(2, p);
        let eps = ExactNum::from_u8(1, p).ldexp(-40, p, RoundingMode::None);
        let above = ExactComplex::new(two.clone(), eps.clone());
        let below = ExactComplex::new(two, eps.neg());
        let ka = above.elliptic_k(p, r, &mut cc);
        let kb = below.elliptic_k(p, r, &mut cc);
        assert!(near(ka.re(), kb.re(), p, 20), "K cut Re");
        assert!(ka.im().is_positive() != kb.im().is_positive() || !tiny(ka.im(), p));
        assert!(near(&ka.im().abs(), &kb.im().abs(), p, 20), "K cut |Im|");
        assert_ne!(ka.im().cmp(kb.im()), Some(0));

        // Incomplete F cut: x=3/2, m=1/2, ±ε imag on x.
        let three_half = ExactNum::from_u8(3, p).div(&ExactNum::from_u8(2, p), p, r);
        let xa = ExactComplex::new(three_half.clone(), eps.clone());
        let xb = ExactComplex::new(three_half, eps.neg());
        let fa = xa.elliptic_f(&mh, p, r, &mut cc);
        let fb = xb.elliptic_f(&mh, p, r, &mut cc);
        assert!(!fa.is_nan() && !fb.is_nan(), "F cut defined");
        assert_ne!(fa.im().cmp(fb.im()), Some(0), "F cut Im differs");
    }
}
