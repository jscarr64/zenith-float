//! Complex Gaussian \({}_2F_1(a,b;c;z)\).
//!
//! Series, Euler, Pfaff, and Kummer / linear transformations in \(\mathbb{C}\).
//! Not the real series at \(\lvert z\rvert\).

use crate::common::util::round_p;
use crate::complex_special::is_nonpos_integer;
use crate::complex_special::nan_pair;
use crate::complex_special::neg_c;
use crate::complex_special::term_negligible;
use crate::complex_special::ziv_complex;
use crate::Consts;
use crate::Error;
use crate::ExactComplex;
use crate::ExactNum;
use crate::RoundingMode;
use crate::WORD_BIT_SIZE;

/// Series term cap (named; same order as the real kernel).
const HYPERGEOM_SERIES_MAX_TERMS: u32 = 10_000;

/// Linear-transform attempts before `NaN`.
const HYPERGEOM_TRANSFORM_MAX: u32 = 8;

fn rm() -> RoundingMode {
    RoundingMode::None
}

fn is_c_zero(z: &ExactComplex) -> bool {
    z.re().is_zero() && z.im().is_zero()
}

fn is_c_one(z: &ExactComplex, p: usize) -> bool {
    z.im().is_zero() && z.re().cmp(&ExactNum::from_u8(1, p)) == Some(0)
}

fn c_eq(a: &ExactComplex, b: &ExactComplex) -> bool {
    a.re().cmp(b.re()) == Some(0) && a.im().cmp(b.im()) == Some(0)
}

fn abs_lt_one(z: &ExactComplex, p: usize) -> bool {
    let az = z.abs(p, rm());
    matches!(az.cmp(&ExactNum::from_u8(1, p)), Some(c) if c < 0)
}

fn re_lt_half(z: &ExactComplex, p: usize) -> bool {
    let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, rm());
    matches!(z.re().cmp(&half), Some(c) if c < 0)
}

fn re_positive_strict(z: &ExactComplex) -> bool {
    z.re().is_positive() && !z.re().is_zero()
}

fn any_nan(args: &[&ExactComplex]) -> bool {
    args.iter().any(|z| z.is_nan())
}

fn as_i32_real_int(z: &ExactComplex, p: usize) -> Option<i32> {
    if !z.im().is_zero() || !z.re().is_int() {
        return None;
    }
    if z.re().is_zero() {
        return Some(0);
    }
    for n in 1i32..=10_000 {
        let w = ExactNum::from_u32(n as u32, p);
        if z.re().abs().cmp(&w) == Some(0) {
            return Some(if z.re().is_negative() { -n } else { n });
        }
    }
    None
}

fn terminating_neg_int(v: &ExactComplex, p: usize) -> bool {
    as_i32_real_int(v, p).is_some_and(|n| n <= 0)
}

fn terminating_after(v: &ExactComplex, next_n: u32, p: usize) -> bool {
    as_i32_real_int(v, p).is_some_and(|m| m <= 0 && next_n as i32 > -m)
}

fn pole_c(c: &ExactComplex, a: &ExactComplex, b: &ExactComplex, p: usize) -> bool {
    let Some(cn) = as_i32_real_int(c, p) else {
        return false;
    };
    if cn > 0 {
        return false;
    }
    let pole_k = -cn;
    let stop_a = as_i32_real_int(a, p).filter(|&n| n <= 0);
    let stop_b = as_i32_real_int(b, p).filter(|&n| n <= 0);
    match (stop_a, stop_b) {
        (Some(sa), _) if -sa < pole_k => false,
        (_, Some(sb)) if -sb < pole_k => false,
        _ => true,
    }
}

fn gamma_fixed(z: &ExactComplex, p: usize, cc: &mut Consts) -> ExactComplex {
    if is_nonpos_integer(z) {
        return nan_pair(Error::InvalidArgument);
    }
    z.gamma_at(p, cc)
}

fn hypergeom_series(
    a: &ExactComplex,
    b: &ExactComplex,
    c: &ExactComplex,
    z: &ExactComplex,
    p: usize,
) -> ExactComplex {
    let one = ExactComplex::one(p);
    let mut term = ExactComplex::one(p);
    let mut sum = ExactComplex::one(p);
    let tiny = ExactNum::from_u8(1, p).ldexp(-((p as i32) - 8), p, rm());
    let n_max = (p.saturating_add(WORD_BIT_SIZE).saturating_add(32))
        .min(HYPERGEOM_SERIES_MAX_TERMS as usize);
    for n in 0..n_max {
        if n > 0 && term_negligible(&term, p) {
            break;
        }
        let nn = ExactComplex::from_real(ExactNum::from_u32(n as u32, p), p);
        let den = c.add(&nn, p, rm()).mul(&one.add(&nn, p, rm()), p, rm());
        if matches!(den.abs(p, rm()).cmp(&tiny), Some(c) if c < 0) {
            return nan_pair(Error::InvalidArgument);
        }
        term = term
            .mul(&a.add(&nn, p, rm()), p, rm())
            .mul(&b.add(&nn, p, rm()), p, rm())
            .div(&den, p, rm())
            .mul(z, p, rm());
        sum = sum.add(&term, p, rm());
        if terminating_after(a, (n + 1) as u32, p) || terminating_after(b, (n + 1) as u32, p) {
            sum.set_inexact(false);
            return sum;
        }
    }
    sum
}

/// Kummer: \({}_2F_1(a,b;c;1)=\Gamma(c)\Gamma(c-a-b)/(\Gamma(c-a)\Gamma(c-b))\)
/// when \(\mathrm{Re}(c-a-b)>0\).
fn kummer_z_one(
    a: &ExactComplex,
    b: &ExactComplex,
    c: &ExactComplex,
    p: usize,
    cc: &mut Consts,
) -> ExactComplex {
    let cab = c.sub(a, p, rm()).sub(b, p, rm());
    if !re_positive_strict(&cab) {
        return nan_pair(Error::InvalidArgument);
    }
    let gc = gamma_fixed(c, p, cc);
    let gcab = gamma_fixed(&cab, p, cc);
    let gca = gamma_fixed(&c.sub(a, p, rm()), p, cc);
    let gcb = gamma_fixed(&c.sub(b, p, rm()), p, cc);
    if gc.is_nan() || gcab.is_nan() || gca.is_nan() || gcb.is_nan() {
        return nan_pair(Error::InvalidArgument);
    }
    gc.mul(&gcab, p, rm()).div(&gca.mul(&gcb, p, rm()), p, rm())
}

/// \({}_2F_1(1,1;2;z)=-\ln(1-z)/z\).
fn is_112(a: &ExactComplex, b: &ExactComplex, c: &ExactComplex, p: usize) -> bool {
    let one = ExactComplex::one(p);
    let two = ExactComplex::from_real(ExactNum::from_u8(2, p), p);
    c_eq(a, &one) && c_eq(b, &one) && c_eq(c, &two)
}

fn f112(z: &ExactComplex, p: usize, cc: &mut Consts) -> ExactComplex {
    let one = ExactComplex::one(p);
    let omz = one.sub(z, p, rm());
    neg_c(&omz.ln(p, rm(), cc)).div(z, p, rm())
}

fn pfaff(
    a: &ExactComplex,
    b: &ExactComplex,
    c: &ExactComplex,
    z: &ExactComplex,
    p: usize,
    cc: &mut Consts,
    left: u32,
) -> ExactComplex {
    let one = ExactComplex::one(p);
    let zm1 = z.sub(&one, p, rm());
    if is_c_zero(&zm1) {
        return nan_pair(Error::InvalidArgument);
    }
    let w = z.div(&zm1, p, rm());
    let pref = one.sub(z, p, rm()).pow(&neg_c(a), p, rm(), cc);
    let cb = c.sub(b, p, rm());
    pref.mul(&hypergeom_at(a, &cb, c, &w, p, cc, left), p, rm())
}

/// A&S 15.3.6: argument \(1-z\).
fn transform_1mz(
    a: &ExactComplex,
    b: &ExactComplex,
    c: &ExactComplex,
    z: &ExactComplex,
    p: usize,
    cc: &mut Consts,
    left: u32,
) -> ExactComplex {
    let one = ExactComplex::one(p);
    let omz = one.sub(z, p, rm());
    let cab = c.sub(a, p, rm()).sub(b, p, rm());
    if is_nonpos_integer(&cab) || is_nonpos_integer(&neg_c(&cab)) {
        return nan_pair(Error::InvalidArgument);
    }
    let gc = gamma_fixed(c, p, cc);
    let gcab = gamma_fixed(&cab, p, cc);
    let gca = gamma_fixed(&c.sub(a, p, rm()), p, cc);
    let gcb = gamma_fixed(&c.sub(b, p, rm()), p, cc);
    let abc = a.add(b, p, rm()).sub(c, p, rm());
    let gabc = gamma_fixed(&abc, p, cc);
    let ga = gamma_fixed(a, p, cc);
    let gb = gamma_fixed(b, p, cc);
    if gc.is_nan()
        || gcab.is_nan()
        || gca.is_nan()
        || gcb.is_nan()
        || gabc.is_nan()
        || ga.is_nan()
        || gb.is_nan()
    {
        return nan_pair(Error::InvalidArgument);
    }
    let a_pref = gc.mul(&gcab, p, rm()).div(&gca.mul(&gcb, p, rm()), p, rm());
    let b_pref = gc.mul(&gabc, p, rm()).div(&ga.mul(&gb, p, rm()), p, rm());
    let cap1 = a.add(b, p, rm()).sub(c, p, rm()).add(&one, p, rm());
    let t1 = a_pref.mul(&hypergeom_at(a, b, &cap1, &omz, p, cc, left), p, rm());
    let f2 = hypergeom_at(
        &c.sub(a, p, rm()),
        &c.sub(b, p, rm()),
        &cab.add(&one, p, rm()),
        &omz,
        p,
        cc,
        left,
    );
    let t2 = b_pref
        .mul(&omz.pow(&cab, p, rm(), cc), p, rm())
        .mul(&f2, p, rm());
    t1.add(&t2, p, rm())
}

/// A&S 15.3.7: argument \(1/z\). Degenerate when \(a=b\).
fn transform_inv(
    a: &ExactComplex,
    b: &ExactComplex,
    c: &ExactComplex,
    z: &ExactComplex,
    p: usize,
    cc: &mut Consts,
    left: u32,
) -> ExactComplex {
    if c_eq(a, b) {
        return nan_pair(Error::InvalidArgument);
    }
    let one = ExactComplex::one(p);
    let inv = one.div(z, p, rm());
    let mz = neg_c(z);
    let t1 = inv_term(a, b, c, &mz, &inv, p, cc, left);
    let t2 = inv_term(b, a, c, &mz, &inv, p, cc, left);
    t1.add(&t2, p, rm())
}

fn inv_term(
    a: &ExactComplex,
    b: &ExactComplex,
    c: &ExactComplex,
    mz: &ExactComplex,
    inv: &ExactComplex,
    p: usize,
    cc: &mut Consts,
    left: u32,
) -> ExactComplex {
    let one = ExactComplex::one(p);
    let gc = gamma_fixed(c, p, cc);
    let gbma = gamma_fixed(&b.sub(a, p, rm()), p, cc);
    let gb = gamma_fixed(b, p, cc);
    let gcma = gamma_fixed(&c.sub(a, p, rm()), p, cc);
    if gc.is_nan() || gbma.is_nan() || gb.is_nan() || gcma.is_nan() {
        return nan_pair(Error::InvalidArgument);
    }
    let pref = gc.mul(&gbma, p, rm()).div(&gb.mul(&gcma, p, rm()), p, rm());
    let pow = mz.pow(&neg_c(a), p, rm(), cc);
    let ac1 = a.sub(c, p, rm()).add(&one, p, rm());
    let ab1 = a.sub(b, p, rm()).add(&one, p, rm());
    pref.mul(&pow, p, rm())
        .mul(&hypergeom_at(a, &ac1, &ab1, inv, p, cc, left), p, rm())
}

fn hypergeom_at(
    a: &ExactComplex,
    b: &ExactComplex,
    c: &ExactComplex,
    z: &ExactComplex,
    p: usize,
    cc: &mut Consts,
    left: u32,
) -> ExactComplex {
    if any_nan(&[a, b, c, z]) {
        return nan_pair(Error::InvalidArgument);
    }
    if is_c_zero(z) || is_c_zero(a) || is_c_zero(b) {
        return ExactComplex::one(p);
    }
    if pole_c(c, a, b, p) {
        return nan_pair(Error::InvalidArgument);
    }
    if is_c_one(z, p) {
        return kummer_z_one(a, b, c, p, cc);
    }
    if is_112(a, b, c, p) {
        return f112(z, p, cc);
    }
    if terminating_neg_int(a, p) || terminating_neg_int(b, p) {
        return hypergeom_series(a, b, c, z, p);
    }
    if abs_lt_one(z, p) {
        return hypergeom_series(a, b, c, z, p);
    }
    if left == 0 {
        return nan_pair(Error::InvalidArgument);
    }
    let next = left - 1;
    let one = ExactComplex::one(p);
    let omz = one.sub(z, p, rm());
    let zm1 = z.sub(&one, p, rm());
    if re_lt_half(z, p) && !is_c_zero(&zm1) {
        return pfaff(a, b, c, z, p, cc, next);
    }
    if abs_lt_one(&omz, p) && !is_nonpos_integer(&c.sub(a, p, rm()).sub(b, p, rm())) {
        return transform_1mz(a, b, c, z, p, cc, next);
    }
    let inv = one.div(z, p, rm());
    if abs_lt_one(&inv, p) && !c_eq(a, b) {
        return transform_inv(a, b, c, z, p, cc, next);
    }
    if !is_c_zero(&zm1) {
        let w = z.div(&zm1, p, rm());
        if abs_lt_one(&w, p) {
            return pfaff(a, b, c, z, p, cc, next);
        }
    }
    nan_pair(Error::InvalidArgument)
}

impl ExactComplex {
    /// Gaussian \({}_2F_1(a=\mathrm{self},b;c;z)\) in \(\mathbb{C}\).
    ///
    /// Series when \(\lvert z\rvert<1\); Pfaff when \(\mathrm{Re}(z)<1/2\);
    /// Euler / \(1-z\) and \(1/z\) linear transforms otherwise. Kummer at \(z=1\)
    /// when \(\mathrm{Re}(c-a-b)>0\). Cut on \([1,+\infty)\) in \(z\) (principal
    /// value from above). Non-positive integer \(c\) (uncanceled) → NaN.
    ///
    /// # Precision
    ///
    /// - Algorithm: series / Euler / Pfaff / Kummer. Caps `HYPERGEOM_SERIES_MAX_TERMS = 10_000`, `HYPERGEOM_TRANSFORM_MAX = 8`.
    /// - Bound: Ziv on each part (`MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    pub fn hypergeom_2f1(
        &self,
        b: &Self,
        c: &Self,
        z: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        if any_nan(&[self, b, c, z]) {
            return nan_pair(Error::InvalidArgument);
        }
        let dest = round_p(p);
        if is_c_zero(z) || is_c_zero(self) || is_c_zero(b) {
            let mut one = ExactComplex::one(dest);
            one.set_inexact(self.inexact() | b.inexact() | c.inexact() | z.inexact());
            return one;
        }
        if pole_c(c, self, b, dest) {
            return nan_pair(Error::InvalidArgument);
        }
        ziv_complex(dest, rm, |pw| {
            hypergeom_at(self, b, c, z, pw, cc, HYPERGEOM_TRANSFORM_MAX)
        })
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
    fn complex_hypergeom_2f1_golds() {
        let p = 256;
        let r = RoundingMode::ToEven;
        let mut cc = Consts::new().unwrap();
        let one = ExactComplex::one(p);
        let two = ExactComplex::from_real(ExactNum::from_u8(2, p), p);
        let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, r);
        let zh = ExactComplex::from_real(half.clone(), p);
        let a_h = zh.clone();
        let z0 = ExactComplex::zero(p);

        let f1 = one.hypergeom_2f1(&one, &two, &zh, p, r, &mut cc);
        let ln2 = cc.ln_2(p, r);
        let want = ln2.mul(&ExactNum::from_u8(2, p), p, r);
        assert!(near(f1.re(), &want, p, 40), "2F1(1,1;2;1/2)=2ln2");
        assert!(tiny(f1.im(), p));

        let fk = a_h.hypergeom_2f1(&a_h, &one, &zh, p, r, &mut cc);
        let k = half.elliptic_k(p, r, &mut cc);
        let pi = cc.pi(p, r);
        let want_k = k.mul(&ExactNum::from_u8(2, p), p, r).div(&pi, p, r);
        assert!(near(fk.re(), &want_k, p, 40), "2F1(1/2,1/2;1;1/2)=2K/π");
        assert!(tiny(fk.im(), p));

        let a = ExactComplex::new(
            ExactNum::from_u8(2, p).div(&ExactNum::from_u8(3, p), p, r),
            ExactNum::from_u8(1, p).div(&ExactNum::from_u8(4, p), p, r),
        );
        let b = ExactComplex::new(
            ExactNum::from_u8(1, p).div(&ExactNum::from_u8(5, p), p, r),
            ExactNum::from_u8(1, p).div(&ExactNum::from_u8(7, p), p, r),
        );
        let c = ExactComplex::new(
            ExactNum::from_u8(3, p).div(&ExactNum::from_u8(2, p), p, r),
            ExactNum::from_u8(1, p).div(&ExactNum::from_u8(9, p), p, r),
        );
        let f0 = a.hypergeom_2f1(&b, &c, &z0, p, r, &mut cc);
        assert!(
            near(f0.re(), &ExactNum::from_u8(1, p), p, 8),
            "2F1(*,*,*;0)=1"
        );
        assert!(tiny(f0.im(), p));

        let z = ExactComplex::new(
            ExactNum::from_u8(3, p).div(&ExactNum::from_u8(10, p), p, r),
            ExactNum::from_u8(1, p).div(&ExactNum::from_u8(10, p), p, r),
        );
        let lhs = a.hypergeom_2f1(&b, &c, &z, p, r, &mut cc);
        let cab = c.sub(&a, p, r).sub(&b, p, r);
        let pref = ExactComplex::one(p).sub(&z, p, r).pow(&cab, p, r, &mut cc);
        let rhs = pref.mul(
            &c.sub(&a, p, r)
                .hypergeom_2f1(&c.sub(&b, p, r), &c, &z, p, r, &mut cc),
            p,
            r,
        );
        assert!(near(lhs.re(), rhs.re(), p, 30), "Euler re");
        assert!(near(lhs.im(), rhs.im(), p, 30), "Euler im");

        let c0 = ExactComplex::zero(p);
        let bad = one.hypergeom_2f1(&one, &c0, &zh, p, r, &mut cc);
        assert!(bad.is_nan(), "c=0 → NaN");

        let two_r = ExactNum::from_u8(2, p);
        let eps = ExactNum::from_u8(1, p).ldexp(-40, p, RoundingMode::None);
        let above = ExactComplex::new(two_r.clone(), eps.clone());
        let below = ExactComplex::new(two_r, eps.neg());
        let fa = one.hypergeom_2f1(&one, &two, &above, p, r, &mut cc);
        let fb = one.hypergeom_2f1(&one, &two, &below, p, r, &mut cc);
        assert!(!fa.is_nan() && !fb.is_nan(), "cut defined");
        assert!(near(fa.re(), fb.re(), p, 20), "cut Re");
        assert_ne!(fa.im().cmp(fb.im()), Some(0), "cut Im differs");
        assert!(near(&fa.im().abs(), &fb.im().abs(), p, 20), "cut |Im|");
    }
}
