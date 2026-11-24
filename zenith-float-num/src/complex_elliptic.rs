//! Complex elliptic integrals via Carlson \(R_F,R_C,R_D,R_J\).
//!
//! Duplication runs in \(\mathbb{C}\) with the principal square root. Parameter
//! \(m=k^2\); incomplete argument \(x=\sin\varphi\). Not the real series at
//! \(\lvert z\rvert\).

use crate::common::util::round_p;
use crate::complex_special::nan_pair;
use crate::complex_special::pi_c;
use crate::complex_special::typed_nan_pair;
use crate::complex_special::ziv_complex;
use crate::defs::Exponent;
use crate::Consts;
use crate::Error;
use crate::ExactComplex;
use crate::ExactNum;
use crate::RoundingMode;
use crate::Sign;
use crate::WORD_BIT_SIZE;
use alloc::vec::Vec;

// Carlson duplication stops when the spread is below the series remainder bound
// (about 2^-(p/3+16)); there is no step cap (P9). R_J shift-integral jets and
// panels stop when the path is covered or the remaining distance stops shrinking.

fn rm() -> RoundingMode {
    RoundingMode::None
}

fn c_u32(n: u32, p: usize) -> ExactComplex {
    ExactComplex::from_real(ExactNum::from_u32(n, p), p)
}

fn c_u64(n: u64, p: usize) -> ExactComplex {
    ExactComplex::from_real(ExactNum::from_u64(n, p), p)
}

fn is_c_zero(z: &ExactComplex) -> bool {
    z.re().is_zero() && z.im().is_zero()
}

fn is_c_one(z: &ExactComplex, _p: usize) -> bool {
    // Compare against a compact 1: a dest-width `from_u8(1, dest)` would allocate
    // before the exact specials K(1)=∞ and E(1)=1 can return.
    z.im().is_zero() && z.re().compare(&ExactNum::from_u8(1, WORD_BIT_SIZE)) == Some(0)
}

/// Compact limit of complete Π(n, 1). Re(1/(1−n)) has the sign of
/// (1 − Re n); Im has the sign of Im n. A vanishing coefficient is a
/// compact zero, not an infinity (Re n = 1, Im n ≠ 0 → 0 + i∞ sign(Im n);
/// real n → ∞ sign(1−n)+0i).
fn pi_complete_m1_infinity(n: &ExactComplex) -> ExactComplex {
    let one = ExactNum::from_u8(1, WORD_BIT_SIZE);
    let re = match one.compare(n.re()) {
        Some(1) => crate::INF_POS,
        Some(-1) => crate::INF_NEG,
        Some(0) => crate::complex::compact_signed_zero(false),
        _ => return nan_pair(Error::InvalidArgument),
    };
    let im = if n.im().is_zero() {
        crate::complex::compact_signed_zero(false)
    } else if n.im().is_negative() {
        crate::INF_NEG
    } else {
        crate::INF_POS
    };
    ExactComplex::new(re, im)
}

fn c_abs(z: &ExactComplex, p: usize) -> ExactNum {
    z.abs(p, rm())
}

#[path = "carlson_taylor.rs"]
mod carlson_taylor;

fn max_abs(a: ExactNum, b: ExactNum, _p: usize) -> ExactNum {
    if matches!(a.compare(&b), Some(c) if c >= 0) {
        a
    } else {
        b
    }
}

fn max_dev3(
    an: &ExactComplex,
    x: &ExactComplex,
    y: &ExactComplex,
    z: &ExactComplex,
    p: usize,
) -> ExactNum {
    let dx = c_abs(&ExactComplex::sub_ref(an, x, p, rm()), p);
    let dy = c_abs(&ExactComplex::sub_ref(an, y, p, rm()), p);
    let dz = c_abs(&ExactComplex::sub_ref(an, z, p, rm()), p);
    max_abs(max_abs(dx, dy, p), dz, p)
}

fn any_nan(args: &[&ExactComplex]) -> bool {
    args.iter().any(|z| z.is_nan())
}

/// `log1p(z)` at the working width, no correct-rounding loop.
/// Same parts as `log1p_word`: Re = ½ log1p(x(2+x)+y²), Im = atan2(y, 1+x).
fn rc_log1p_none(z: &ExactComplex, p: usize, cc: &mut Consts) -> ExactComplex {
    let x = z.re();
    let y = z.im();
    let two_x = ExactNum::from_u8(2, p).add(x, p, rm());
    let t = ExactNum::mul_ref(x, &two_x, p, rm()).add(&y.sqr(p, rm()), p, rm());
    let re = t.log1p(p, rm(), cc).ldexp(-1, p, rm());
    let one_x = ExactNum::from_u8(1, p).add(x, p, rm());
    let im = y.atan2(&one_x, p, rm(), cc);
    ExactComplex::new(re, im)
}

/// `atan(z)` at the working width, no correct-rounding loop.
/// `atan(z)=(i/2) log1p(2z/(i−z))`. Forming `ln((i+z)/(i−z))` first loses
/// ~p/6 bits when `z` is small (the late R_J increments, where `1+e_m≈1`).
fn rc_atan_none(z: &ExactComplex, p: usize, cc: &mut Consts) -> ExactComplex {
    let i = ExactComplex::new(ExactNum::new(p), ExactNum::from_u8(1, p));
    let twoz = z.mul(&c_u32(2, p), p, rm());
    let w = twoz.div(&i.sub(z, p, rm()), p, rm());
    let ln = rc_log1p_none(&w, p, cc);
    let half_i = ExactComplex::new(ExactNum::new(p), ExactNum::from_u8(1, p).ldexp(-1, p, rm()));
    half_i.mul(&ln, p, rm())
}

/// \(R_C(x,y)=R_F(x,y,y)\).
fn carlson_rc(x0: &ExactComplex, y0: &ExactComplex, p: usize, cc: &mut Consts) -> ExactComplex {
    if any_nan(&[x0, y0]) {
        return ExactComplex::new(x0.re().clone(), y0.im().clone());
    }
    if is_c_zero(y0) {
        return nan_pair(Error::InvalidArgument);
    }
    // Signed-zero negative reals take the C99 side of the principal-sqrt cut.
    let x = leave_cut(x0, p);
    let y = leave_cut(y0, p);
    // Shortcut 1/√x only when |x − y| < 2^-(p+8)|x| (first-order error (y − x)/(3x)); the
    // spread-based `close_enough` test (2^-(p/3+16)) was too loose for R_J's RC terms.
    let dxy = ExactComplex::sub_ref(&x, &y, p, rm());
    let rel = c_abs(&x, p).ldexp(-((p as Exponent) + 8), p, rm());
    if !is_c_zero(&x) && matches!(c_abs(&dxy, p).compare(&rel), Some(o) if o < 0) {
        return c_u32(1, p).div(&x.sqrt(p, rm(), cc), p, rm());
    }
    if is_c_zero(&x) {
        let mut hp = pi_c(p, cc);
        hp = hp.div(&c_u32(2, p), p, rm());
        return hp.div(&y.sqrt(p, rm(), cc), p, rm());
    }
    // Proven source: the real `carlson_rc` in ops/special.rs (DLMF 19.2.18–19).
    //   y>x: atan(√((y−x)/x)) / √(y−x)
    //   y<x: atanh(√((x−y)/x)) / √(x−y)
    // The previous path called R_F(x,y,y), a full 3-argument duplication, on
    // every R_J increment. Complete and incomplete Π then timed out under
    // ZFS_WALL=120 at dest p=8192 (Harper T4 S2).
    // Evaluate the elementary form at the working width with RoundingMode::None;
    // public complex atan/atanh would start a nested correct-rounding Ziv.
    let both_real_pos =
        is_real_c(&x) && is_real_c(&y) && x.re().is_positive() && y.re().is_positive();
    if both_real_pos {
        match y.re().compare(x.re()) {
            Some(1) => {
                let d = ExactNum::sub_ref(y.re(), x.re(), p, rm());
                let arg = ExactNum::div_ref(&d, x.re(), p, rm()).sqrt(p, rm());
                let val = arg.atan(p, rm(), cc).div(&d.sqrt(p, rm()), p, rm());
                ExactComplex::from_real(val, p)
            }
            Some(-1) => {
                let d = ExactNum::sub_ref(x.re(), y.re(), p, rm());
                let arg = ExactNum::div_ref(&d, x.re(), p, rm()).sqrt(p, rm());
                let val = arg.atanh(p, rm(), cc).div(&d.sqrt(p, rm()), p, rm());
                ExactComplex::from_real(val, p)
            }
            _ => c_u32(1, p).div(&x.sqrt(p, rm(), cc), p, rm()),
        }
    } else {
        let d = ExactComplex::sub_ref(&y, &x, p, rm());
        let arg = ExactComplex::sub_ref(&y, &x, p, rm())
            .div(&x, p, rm())
            .sqrt(p, rm(), cc);
        rc_atan_none(&arg, p, cc).div(&d.sqrt(p, rm(), cc), p, rm())
    }
}

/// Homogeneity normalization for the Carlson integrals (Spence PR #40 R2). `R_F` has degree `-1/2`, `R_D` and
/// `R_J` degree `-3/2` (DLMF 19.16), and scaling every argument by the positive real `4^-t` keeps each on its own
/// side of every cut. `t` is half the largest part exponent (floor), so after scaling the largest part is below
/// `2^2` and `|A|` is far from both ends of the exponent range. Without this, `|A|^2` overflows for arguments near
/// `2^(2^61)` (`c_abs(A) = +Inf` stops the duplication at `k = 0` on an unconverged series, or `X = (A-x)/A`
/// overflows to NaN), and underflows near `2^-(2^61)` (`c_abs(A) = 0`). An argument more than about `2^62`
/// binades below the largest underflows when scaled; it is then below `2^-(2^61)` relative, under every
/// precision's last bit, except that two such arguments make two zeros, which `R_F`/`R_J` refuse.
fn carlson_norm_exp(args: &[&ExactComplex]) -> Exponent {
    args.iter()
        .flat_map(|z| [z.re(), z.im()])
        .filter(|x| !x.is_zero())
        .filter_map(|x| x.exponent())
        .max()
        .map_or(0, |e| e.div_euclid(2))
}

/// `z 4^-t`, exact (the width keeps every bit of `z`).
fn carlson_scaled(z: &ExactComplex, t: Exponent, p: usize) -> ExactComplex {
    let w = [z.re().precision(), z.im().precision()]
        .into_iter()
        .flatten()
        .fold(p, usize::max);
    z.ldexp(t.saturating_mul(-2), w, rm())
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
    let t = carlson_norm_exp(&[x0, y0, z0]);
    if t != 0 {
        // R_F(4^t x, 4^t y, 4^t z) = 2^-t R_F(x, y, z)
        let (xs, ys, zs) = (
            carlson_scaled(x0, t, p),
            carlson_scaled(y0, t, p),
            carlson_scaled(z0, t, p),
        );
        return carlson_rf(&xs, &ys, &zs, p, cc).ldexp(t.saturating_neg(), p, rm());
    }
    let four = c_u32(4, p);
    let three = c_u32(3, p);
    // Signed-zero negative reals take the C99 side of the principal-sqrt cut.
    let mut x = leave_cut(x0, p);
    let mut y = leave_cut(y0, p);
    let mut z = leave_cut(z0, p);
    let real = is_real_c(&x) && is_real_c(&y) && is_real_c(&z);
    let order = carlson_taylor::rf_order(p, real);
    let mut prev_d: Option<ExactNum> = None;
    let mut k = 0usize;
    loop {
        let an = ExactComplex::add_ref(&x, &y, p, rm())
            .add(&z, p, rm())
            .div(&three, p, rm());
        let dmax = max_dev3(&an, &x, &y, &z, p);
        // Quick stop: (dmax/|A|)^order < 2^{-p}. This is zf's rule, not
        // FLINT's: rf.c:241-252 uses max(|x-y|,|x-z|) / (lower bound of |x|)
        // after a duplication step. Both are scale-invariant; a floor of 1 on
        // the denominator stops too early when |A|<1 and drops trailing bits.
        // The rigorous gate is the truncation bound below (rf.c:285-295).
        let scale = c_abs(&an, p);
        let no_improve = k > 2
            && prev_d
                .as_ref()
                .is_some_and(|pd| matches!(dmax.compare(pd), Some(c) if c > 0));
        let forced = no_improve || k >= p;
        if carlson_taylor::rel_pow_below_prec(&dmax, &scale, order, p, p) || forced {
            let xx = ExactComplex::sub_ref(&an, &x, p, rm()).div(&an, p, rm());
            let yy = ExactComplex::sub_ref(&an, &y, p, rm()).div(&an, p, rm());
            let zz = ExactComplex::sub_ref(&an, &z, p, rm()).div(&an, p, rm());
            // FLINT rf.c:285-295 truncation bound; below 2^-p, or keep duplicating.
            let eps = max_abs(max_abs(c_abs(&xx, p), c_abs(&yy, p), p), c_abs(&zz, p), p);
            if forced || carlson_taylor::trunc_bound_below_prec(&eps, order, p) {
                let e2 = xx.mul(&yy, p, rm()).sub(&zz.mul(&zz, p, rm()), p, rm());
                let e3 = xx.mul(&yy, p, rm()).mul(&zz, p, rm());
                let series = carlson_taylor::rf_taylor_sum(&e2, &e3, order, p);
                return series.div(&an.sqrt(p, rm(), cc), p, rm());
            }
        }
        prev_d = Some(dmax);
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
        k += 1;
    }
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
    let t = carlson_norm_exp(&[x0, y0, z0]);
    if t != 0 {
        // R_D(4^t x, 4^t y, 4^t z) = 2^-3t R_D(x, y, z)
        let (xs, ys, zs) = (
            carlson_scaled(x0, t, p),
            carlson_scaled(y0, t, p),
            carlson_scaled(z0, t, p),
        );
        return carlson_rd(&xs, &ys, &zs, p, cc).ldexp(t.saturating_mul(-3), p, rm());
    }
    let four = c_u32(4, p);
    let three = c_u32(3, p);
    let five = c_u32(5, p);
    // Signed-zero negative reals take the C99 side of the principal-sqrt cut.
    let mut x = leave_cut(x0, p);
    let mut y = leave_cut(y0, p);
    let mut z = leave_cut(z0, p);
    let mut sum = ExactComplex::zero(p);
    let mut fac = c_u32(1, p);
    let real = is_real_c(&x) && is_real_c(&y) && is_real_c(&z);
    let order = carlson_taylor::rd_order(p, real);
    let mut prev_d: Option<ExactNum> = None;
    let mut k = 0usize;
    loop {
        let an = ExactComplex::add_ref(&x, &y, p, rm())
            .add(&three.mul(&z, p, rm()), p, rm())
            .div(&five, p, rm());
        let dmax = max_dev3(&an, &x, &y, &z, p);
        // Quick stop: (dmax/|A|)^order < 2^{-p}. This is zf's rule, not
        // FLINT's: rf.c:241-252 uses max(|x-y|,|x-z|) / (lower bound of |x|)
        // after a duplication step. Both are scale-invariant; a floor of 1 on
        // the denominator stops too early when |A|<1 and drops trailing bits.
        // The rigorous gate is the truncation bound below (rf.c:285-295).
        let scale = c_abs(&an, p);
        let no_improve = k > 2
            && prev_d
                .as_ref()
                .is_some_and(|pd| matches!(dmax.compare(pd), Some(c) if c > 0));
        let forced = no_improve || k >= p;
        if (carlson_taylor::rel_pow_below_prec(&dmax, &scale, order, p, p) || forced)
            && (forced || rd_trunc_ok(&an, &x, &y, &z, order, p))
        {
            let series = rd_series_arb(&an, &x, &y, &z, order, p, cc);
            return three
                .mul(&sum, p, rm())
                .add(&fac.mul(&series, p, rm()), p, rm());
        }
        prev_d = Some(dmax);
        let sx = x.sqrt(p, rm(), cc);
        let sy = y.sqrt(p, rm(), cc);
        let sz = z.sqrt(p, rm(), cc);
        let lam = sx
            .mul(&sy, p, rm())
            .add(&sy.mul(&sz, p, rm()), p, rm())
            .add(&sz.mul(&sx, p, rm()), p, rm());
        let zl = z.add(&lam, p, rm());
        sum = sum.add(&fac.div(&sz.mul(&zl, p, rm()), p, rm()), p, rm());
        fac = fac.div(&four, p, rm());
        x = x.add(&lam, p, rm()).div(&four, p, rm());
        y = y.add(&lam, p, rm()).div(&four, p, rm());
        z = zl.div(&four, p, rm());
        k += 1;
    }
}

/// FLINT rj.c:466-477 truncation bound for the R_D-as-R_J sum (`P = Z`): `eps = max(|X|, |Y|, |Z|)` with
/// `X = (A - x)/A` and so on, as [`rd_series_arb`] forms them.
fn rd_trunc_ok(
    an: &ExactComplex,
    x: &ExactComplex,
    y: &ExactComplex,
    z: &ExactComplex,
    order: usize,
    p: usize,
) -> bool {
    let dev = |v: &ExactComplex| c_abs(&ExactComplex::sub_ref(an, v, p, rm()).div(an, p, rm()), p);
    let eps = max_abs(max_abs(dev(x), dev(y), p), dev(z), p);
    carlson_taylor::trunc_bound_below_prec(&eps, order, p)
}

/// Remainder of \(R_D\) as the \(R_J\) Taylor sum with \(p=z\) (Arb).
fn rd_series_arb(
    an: &ExactComplex,
    x: &ExactComplex,
    y: &ExactComplex,
    z: &ExactComplex,
    order: usize,
    p: usize,
    cc: &mut Consts,
) -> ExactComplex {
    let two = c_u32(2, p);
    let four = c_u32(4, p);
    let xx = ExactComplex::sub_ref(an, x, p, rm()).div(an, p, rm());
    let yy = ExactComplex::sub_ref(an, y, p, rm()).div(an, p, rm());
    let zz = ExactComplex::sub_ref(an, z, p, rm()).div(an, p, rm());
    // P = Z here.
    let pp = &zz;
    let e2 = xx
        .mul(&yy, p, rm())
        .add(&xx.mul(&zz, p, rm()), p, rm())
        .add(&yy.mul(&zz, p, rm()), p, rm())
        .sub(&c_u32(3, p).mul(&pp.mul(pp, p, rm()), p, rm()), p, rm());
    let xyz = xx.mul(&yy, p, rm()).mul(&zz, p, rm());
    let e3 = ExactComplex::add_ref(&xyz, &two.mul(&e2, p, rm()).mul(pp, p, rm()), p, rm()).add(
        &four.mul(&pp.mul(pp, p, rm()).mul(pp, p, rm()), p, rm()),
        p,
        rm(),
    );
    let e4 = two
        .mul(&xyz, p, rm())
        .add(&e2.mul(pp, p, rm()), p, rm())
        .add(
            &c_u32(3, p).mul(&pp.mul(pp, p, rm()).mul(pp, p, rm()), p, rm()),
            p,
            rm(),
        )
        .mul(pp, p, rm());
    let e5 = xyz.mul(&pp.mul(pp, p, rm()), p, rm());
    let series = carlson_taylor::rj_taylor_sum(&e2, &e3, &e4, &e5, order, p);
    series.div(&an.mul(&an.sqrt(p, rm(), cc), p, rm()), p, rm())
}

fn is_real_c(z: &ExactComplex) -> bool {
    z.im().is_zero()
}

/// `1 - w`, with `Im(1-w) = -Im(w)` including signed zeros. IEEE 754
/// `(+0) - (+0) = +0`, so a plain subtract would lose the C99 side of a
/// real argument on a cut.
fn one_minus(w: &ExactComplex, p: usize) -> ExactComplex {
    let d = ExactComplex::one(p).sub(w, p, rm());
    if !is_real_c(w) {
        return d;
    }
    // Im w is a signed zero; its negation is the C99 side of 1 - w.
    let (re, _) = d.into_parts();
    ExactComplex::new(re, w.im().neg_ref())
}

/// If `z` lies on the open cut (−∞,0), move it off the cut by
/// \(2^{-(\mathrm{dest}+16)}i\), taking the side given by the sign of
/// `Im z` (C99 / IEEE: `+0` is the upper side, `-0` the lower). At
/// (x,y,z,p)=(0,−1,1,−i) the first \(R_J\) increment has \(e_m\approx-1\)
/// and \(1+e_m=O(\varepsilon)\), so the perturbation is sized from the
/// destination precision, not the Ziv work width.
fn leave_cut(z: &ExactComplex, dest: usize) -> ExactComplex {
    if is_real_c(z) && z.re().is_negative() {
        let mut eps = ExactNum::from_u8(1, dest).ldexp(-((dest as Exponent) + 16), dest, rm());
        if z.im().is_negative() {
            eps.set_sign(Sign::Neg);
        }
        ExactComplex::new(z.re().clone(), eps)
    } else {
        z.clone()
    }
}

/// If `z` is real, move it off the real axis by \(2^{-(\mathrm{dest}+16)}i\),
/// taking the C99 side given by the sign of `Im z`. Products such as
/// \(m x^2\) can otherwise collapse a −0 factor against a +0 factor to +0
/// (IEEE `(−0)+(+0)=+0`), putting \(1-mx^2\) on the wrong side of the cut.
fn bump_real(z: &ExactComplex, dest: usize) -> ExactComplex {
    if is_real_c(z) {
        let mut eps = ExactNum::from_u8(1, dest).ldexp(-((dest as Exponent) + 16), dest, rm());
        if z.im().is_negative() {
            eps.set_sign(Sign::Neg);
        }
        ExactComplex::new(z.re().clone(), eps)
    } else {
        z.clone()
    }
}

fn on_nonpos_real(z: &ExactComplex) -> bool {
    is_real_c(z) && !z.re().is_positive()
}

/// Work precision needed so \(1+e_m=O(2^{-(\mathrm{dest}+16)})\) keeps `dest`
/// good bits through the first \(R_C(1,1+e_m)\) cancellation.
fn rj_cut_work(dest: usize, pw: usize) -> usize {
    let need = dest.saturating_mul(2).saturating_add(128);
    pw.max(need)
}

/// Drop O(ε) leave-cut leakage. The dest-rounded value is the ε→0 limit;
/// a leftover \(2^{-(\mathrm{dest}+16)}\) part is representable at dest
/// precision but is not part of the limit (Π(1+i, 2) is real).
fn flush_cut_bias(z: ExactComplex, dest: usize) -> ExactComplex {
    let thresh = ExactNum::from_u8(1, dest).ldexp(-((dest as Exponent) + 8), dest, rm());
    let flush = |c: ExactNum| {
        if matches!(c.abs_ref().compare(&thresh), Some(o) if o < 0) {
            ExactNum::new(dest)
        } else {
            c
        }
    };
    let (re, im) = z.into_parts();
    ExactComplex::new(flush(re), flush(im))
}

fn min_num(a: ExactNum, b: ExactNum) -> ExactNum {
    if matches!(a.compare(&b), Some(c) if c <= 0) {
        a
    } else {
        b
    }
}

/// [`min_num`] of borrowed values: a reference to the smaller one.
fn min_ref<'a>(a: &'a ExactNum, b: &'a ExactNum) -> &'a ExactNum {
    if matches!(a.compare(b), Some(c) if c <= 0) {
        a
    } else {
        b
    }
}

fn c_neg(z: &ExactComplex, p: usize) -> ExactComplex {
    ExactComplex::zero(p).sub(z, p, rm())
}

/// Distance from `w` to the ray \((-\infty,0]\).
fn dist_to_neg_ray(w: &ExactComplex, p: usize) -> ExactNum {
    if w.re().is_positive() {
        c_abs(w, p)
    } else {
        w.im().abs_ref()
    }
}

/// Distance from `t` to the branch cuts of \(\sqrt{t+x}\) etc. and the pole at \(-p\).
fn rj_sing_dist(
    t: &ExactComplex,
    x: &ExactComplex,
    y: &ExactComplex,
    z: &ExactComplex,
    pv: &ExactComplex,
    p: usize,
) -> ExactNum {
    let tx = ExactComplex::add_ref(t, x, p, rm());
    let ty = ExactComplex::add_ref(t, y, p, rm());
    let tz = ExactComplex::add_ref(t, z, p, rm());
    let tp = ExactComplex::add_ref(t, pv, p, rm());
    min_num(
        min_num(dist_to_neg_ray(&tx, p), dist_to_neg_ray(&ty, p)),
        min_num(dist_to_neg_ray(&tz, p), c_abs(&tp, p)),
    )
}

/// Whether `s = num/den` lies in \([0,1]\).
fn ratio_in_unit(num: &ExactNum, den: &ExactNum) -> bool {
    if num.is_zero() {
        return true;
    }
    if den.is_zero() {
        return false;
    }
    if num.is_negative() != den.is_negative() {
        return false;
    }
    matches!(num.abs_cmp(den), Some(c) if c <= 0)
}

/// The segment \([start, start+step]\) meets the ray \((-\infty,0]\).
fn seg_hits_neg_ray(start: &ExactComplex, step: &ExactComplex, p: usize) -> bool {
    if step.im().is_zero() {
        // Real projection of the segment, including paths that run
        // parallel to the cut at a tiny Im (leave-cut ε). Those must
        // detour; a crawl of step 2^{-(dest+16)} along the cut hangs.
        let a = start.re();
        let b = ExactNum::add_ref(start.re(), step.re(), p, rm());
        let (lo, hi) = if matches!(a.compare(&b), Some(c) if c <= 0) { (a, &b) } else { (&b, a) };
        // A mere touch of 0 at one endpoint (integrable 1/√t) is not a hit
        // provided the rest of the segment has Re > 0.
        if lo.is_zero() && hi.is_positive() {
            return false;
        }
        !lo.is_positive()
    } else {
        // Im(start) + s Im(step) = 0 ⇒ s = −Im(start)/Im(step).
        let s_num = ExactNum::new(p).sub(start.im(), p, rm());
        if s_num.is_zero() || !ratio_in_unit(&s_num, step.im()) {
            return false;
        }
        let s_re = s_num.mul(step.re(), p, rm()).div(step.im(), p, rm());
        let re = ExactNum::add_ref(start.re(), &s_re, p, rm());
        !re.is_positive()
    }
}

fn path_hits_sings(
    from: &ExactComplex,
    step: &ExactComplex,
    x: &ExactComplex,
    y: &ExactComplex,
    z: &ExactComplex,
    pv: &ExactComplex,
    p: usize,
) -> bool {
    // (t + w) walks [w+from, w+from+step] and must miss (−∞,0] for w in {x,y,z};
    // t must miss the pole −p, i.e. (from + s step + p) ≠ 0 for s in (0,1].
    let wx = ExactComplex::add_ref(from, x, p, rm());
    let wy = ExactComplex::add_ref(from, y, p, rm());
    let wz = ExactComplex::add_ref(from, z, p, rm());
    if seg_hits_neg_ray(&wx, step, p)
        || seg_hits_neg_ray(&wy, step, p)
        || seg_hits_neg_ray(&wz, step, p)
    {
        return true;
    }
    // Pole: from + s step + p = 0 ⇒ s step = −(from+p). s = −(from+p)/step in (0,1].
    let fp = ExactComplex::add_ref(from, pv, p, rm());
    if is_c_zero(step) {
        return is_c_zero(&fp);
    }
    let s = c_neg(&fp, p).div(step, p, rm());
    is_real_c(&s) && s.re().is_positive() && ratio_in_unit(s.re(), &ExactNum::from_u8(1, p))
        || is_c_zero(&s)
}

fn shift_n(
    x: &ExactComplex,
    y: &ExactComplex,
    z: &ExactComplex,
    pv: &ExactComplex,
    p: usize,
) -> (ExactComplex, Option<ExactComplex>) {
    // N_re = 1 + max(0, −Re of each argument), so Re(arg+N) ≥ 1.
    let zero = ExactNum::new(p);
    let one = ExactNum::from_u8(1, p);
    let rx = min_ref(min_ref(x.re(), y.re()), min_ref(z.re(), pv.re()));
    let n_re = if rx.is_negative() { one.sub(rx, p, rm()) } else { one };
    let n_real = ExactComplex::from_real(n_re, p);
    let half_im = ExactComplex::new(zero, ExactNum::from_u8(1, p).ldexp(-1, p, rm()));
    let n_up = ExactComplex::add_ref(&n_real, &half_im, p, rm());
    let n_dn = ExactComplex::sub_ref(&n_real, &half_im, p, rm());
    let origin = ExactComplex::zero(p);
    if !path_hits_sings(&origin, &n_real, x, y, z, pv, p) {
        (n_real, None)
    } else if !path_hits_sings(&origin, &n_up, x, y, z, pv, p) {
        (n_up, None)
    } else if !path_hits_sings(&origin, &n_dn, x, y, z, pv, p) {
        (n_dn, None)
    } else {
        // Two segments: 0 → i/2 → N_re+i/2 (or the lower half-plane).
        if !path_hits_sings(&origin, &half_im, x, y, z, pv, p)
            && !path_hits_sings(&half_im, &n_real, x, y, z, pv, p)
        {
            (n_up, Some(half_im))
        } else {
            let mh = c_neg(&half_im, p);
            (n_dn, Some(mh))
        }
    }
}

fn inv_sqrt_series(
    shift: &ExactComplex,
    p: usize,
    cc: &mut Consts,
    kmax: usize,
) -> Vec<ExactComplex> {
    let mut s = Vec::with_capacity(kmax);
    s.push(c_u32(1, p).div(&shift.sqrt(p, rm(), cc), p, rm()));
    for k in 0..kmax.saturating_sub(1) {
        // a_{k+1} = a_k * (−(2k+1)/2) / ((k+1) shift)
        let odd = c_u64((k as u64).saturating_mul(2).saturating_add(1), p);
        let num = c_neg(&odd.div(&c_u32(2, p), p, rm()), p);
        let den = c_u64((k as u64).saturating_add(1), p).mul(shift, p, rm());
        s.push(s[k].mul(&num, p, rm()).div(&den, p, rm()));
    }
    s
}

fn inv_series(shift: &ExactComplex, p: usize, kmax: usize) -> Vec<ExactComplex> {
    let mut s = Vec::with_capacity(kmax);
    s.push(c_u32(1, p).div(shift, p, rm()));
    for k in 0..kmax.saturating_sub(1) {
        s.push(c_neg(&s[k], p).div(shift, p, rm()));
    }
    s
}

fn conv_series(a: &[ExactComplex], b: &[ExactComplex], p: usize) -> Vec<ExactComplex> {
    let n = a.len().min(b.len());
    let mut c = Vec::with_capacity(n);
    for k in 0..n {
        let mut acc = ExactComplex::zero(p);
        for i in 0..=k {
            acc = acc.add(&a[i].mul(&b[k - i], p, rm()), p, rm());
        }
        c.push(acc);
    }
    c
}

/// Proven series length for one RJ panel: want r^k < 2^{-(p+8)} with
/// r = |h/2| / (9/10 ρ). Falls back to the conservative p+32 bound if r
/// is not strictly below 1.
fn rj_panel_kmax(wabs: &ExactNum, rho: &ExactNum, p: usize) -> usize {
    let cap = (p + 32).max(16);
    let den = ExactNum::from_u8(9, p)
        .div(&ExactNum::from_u8(10, p), p, rm())
        .mul(rho, p, rm());
    if den.is_zero() {
        return cap;
    }
    let r = ExactNum::div_ref(wabs, &den, p, rm());
    let one = ExactNum::from_u8(1, p);
    if r.is_zero() || r.is_negative() || matches!(r.compare(&one), Some(c) if c >= 0) {
        return cap;
    }
    let thresh = ExactNum::from_u8(1, p).ldexp(-((p as Exponent) + 8), p, rm());
    // r^k; `None` while k = 1 (r read in place).
    let mut uk_slot: Option<ExactNum> = None;
    let mut k = 1usize;
    loop {
        let uk = uk_slot.as_ref().unwrap_or(&r);
        match uk.compare(&thresh) {
            Some(c) if c >= 0 => {
                if k >= cap {
                    return cap;
                }
                uk_slot = Some(ExactNum::mul_ref(uk, &r, p, rm()));
                k += 1;
            }
            Some(_) => break,
            None => return cap,
        }
    }
    k.saturating_add(4).max(8).min(cap)
}

fn rj_panel(
    a: &ExactComplex,
    h: &ExactComplex,
    [x, y, z]: [&ExactComplex; 3],
    pv: &ExactComplex,
    p: usize,
    cc: &mut Consts,
) -> ExactComplex {
    // Midpoint c = a + h/2. Expand f(c+s), integrate s from −h/2 to h/2.
    let half = h.div(&c_u32(2, p), p, rm());
    let c = ExactComplex::add_ref(a, &half, p, rm());
    let rho = rj_sing_dist(&c, x, y, z, pv, p);
    let wabs = c_abs(&half, p);
    // r = |h/2| / (9/10 ρ); the caller should have |h/2| ≤ ρ/2 so r ≤ 5/9.
    // Take only enough terms that r^k < 2^{-(p+8)} (proven remainder).
    let kmax = rj_panel_kmax(&wabs, &rho, p);
    let sx = inv_sqrt_series(&ExactComplex::add_ref(&c, x, p, rm()), p, cc, kmax);
    let sy = inv_sqrt_series(&ExactComplex::add_ref(&c, y, p, rm()), p, cc, kmax);
    let sz = inv_sqrt_series(&ExactComplex::add_ref(&c, z, p, rm()), p, cc, kmax);
    let ip = inv_series(&ExactComplex::add_ref(&c, pv, p, rm()), p, kmax);
    let f = conv_series(&conv_series(&conv_series(&sx, &sy, p), &sz, p), &ip, p);
    // ∫_{-w}^{w} s^k ds = 0 (k odd), 2 w^{k+1}/(k+1) (k even).
    let mut acc = ExactComplex::zero(p);
    // w^{k+1}; `None` at k = 0 (half read in place).
    let mut wpow_slot: Option<ExactComplex> = None;
    for (k, fk) in f.iter().enumerate() {
        let wpow = wpow_slot.as_ref().unwrap_or(&half);
        if k % 2 == 0 {
            let term = fk.mul(&c_u32(2, p), p, rm()).mul(wpow, p, rm()).div(
                &c_u64((k as u64).saturating_add(1), p),
                p,
                rm(),
            );
            acc = acc.add(&term, p, rm());
        }
        wpow_slot = Some(wpow.mul(&half, p, rm()));
    }
    acc
}

/// `w[k]`, or `zero` past the end of the series.
fn series_coeff<'a>(w: &'a [ExactComplex], k: usize, zero: &'a ExactComplex) -> &'a ExactComplex {
    w.get(k).unwrap_or(zero)
}

/// Coefficients of 1/w(s) up to s^{kmax-1}.
fn inv_of_series(w: &[ExactComplex], p: usize, kmax: usize) -> Vec<ExactComplex> {
    let zero = ExactComplex::zero(p);
    let w0 = series_coeff(w, 0, &zero);
    let mut r = Vec::with_capacity(kmax);
    r.push(c_u32(1, p).div(w0, p, rm()));
    for n in 1..kmax {
        let mut s = ExactComplex::zero(p);
        for (i, ri) in r.iter().enumerate() {
            s = s.add(&ri.mul(series_coeff(w, n - i, &zero), p, rm()), p, rm());
        }
        r.push(c_neg(&s, p).div(w0, p, rm()));
    }
    r
}

/// Coefficients of 1/√w(s) up to s^{kmax-1}. Uses (r*r)*w = 1.
fn inv_sqrt_of_series(
    w: &[ExactComplex],
    p: usize,
    cc: &mut Consts,
    kmax: usize,
) -> Vec<ExactComplex> {
    let zero = ExactComplex::zero(p);
    let w0 = series_coeff(w, 0, &zero);
    let mut r = Vec::with_capacity(kmax);
    r.push(c_u32(1, p).div(&w0.sqrt(p, rm(), cc), p, rm()));
    // Cache r*r so each new coefficient is O(n), not O(n²).
    let mut rr = Vec::with_capacity(kmax);
    rr.push(r[0].mul(&r[0], p, rm()));
    for n in 1..kmax {
        let mut known = ExactComplex::zero(p);
        for (i, rri) in rr.iter().enumerate() {
            known = known.add(&rri.mul(series_coeff(w, n - i, &zero), p, rm()), p, rm());
        }
        // rr_n = 2 r0 r_n + Q, Q = sum_{j=1}^{n-1} r_j r_{n-j}.
        // (rr w)_n = (2 r0 r_n + Q) w0 + known = 0.
        let mut q = ExactComplex::zero(p);
        for j in 1..n {
            q = q.add(&r[j].mul(&r[n - j], p, rm()), p, rm());
        }
        let rhs = known.add(&q.mul(w0, p, rm()), p, rm());
        let coef = c_u32(2, p).mul(&r[0], p, rm()).mul(w0, p, rm());
        r.push(c_neg(&rhs, p).div(&coef, p, rm()));
        let mut rrn = ExactComplex::zero(p);
        for j in 0..=n {
            rrn = rrn.add(&r[j].mul(&r[n - j], p, rm()), p, rm());
        }
        rr.push(rrn);
    }
    r
}

/// val + λ(c+s)² = (val + λ c²) + 2 λ c s + λ s².
fn quad_shift_series(
    val: &ExactComplex,
    lam: &ExactComplex,
    c: &ExactComplex,
    p: usize,
) -> [ExactComplex; 3] {
    let c2 = c.mul(c, p, rm());
    let w0 = ExactComplex::add_ref(val, &lam.mul(&c2, p, rm()), p, rm());
    let w1 = c_u32(2, p).mul(lam, p, rm()).mul(c, p, rm());
    [w0, w1, lam.clone()]
}

/// Distance in the u-plane, t = λ u², from `u` to the cuts of √(y+t),
/// √(z+t) and the pole at t = −p. If |Δu| ≤ r then
/// |Δt| = |λ| |2 u Δu + Δu²| ≤ |λ|(2|u|r + r²), and we take the r
/// that makes this equal to ρ_t: r = −|u| + √(|u|² + ρ_t/|λ|).
fn rj_u2_rho(
    u: &ExactComplex,
    lam: &ExactComplex,
    y: &ExactComplex,
    z: &ExactComplex,
    pv: &ExactComplex,
    p: usize,
    cc: &mut Consts,
) -> ExactNum {
    let t = lam.mul(&u.mul(u, p, rm()), p, rm());
    let rho_t = min_num(
        min_num(
            dist_to_neg_ray(&ExactComplex::add_ref(&t, y, p, rm()), p),
            dist_to_neg_ray(&ExactComplex::add_ref(&t, z, p, rm()), p),
        ),
        c_abs(&ExactComplex::add_ref(&t, pv, p, rm()), p),
    );
    let uabs = c_abs(u, p);
    let labs = c_abs(lam, p);
    if labs.is_zero() || rho_t.is_zero() || rho_t.is_negative() {
        return ExactNum::new(p);
    }
    let disc = ExactNum::mul_ref(&uabs, &uabs, p, rm()).add(&rho_t.div(&labs, p, rm()), p, rm());
    let sqrt_disc = ExactComplex::from_real(disc, p).sqrt(p, rm(), cc);
    ExactNum::sub_ref(sqrt_disc.re(), &uabs, p, rm())
}

fn rj_u2_panel(
    u0: &ExactComplex,
    h: &ExactComplex,
    [lam, two_sl]: [&ExactComplex; 2],
    [y, z]: [&ExactComplex; 2],
    pv: &ExactComplex,
    p: usize,
    cc: &mut Consts,
) -> ExactComplex {
    let half = h.div(&c_u32(2, p), p, rm());
    let c = ExactComplex::add_ref(u0, &half, p, rm());
    let rho = rj_u2_rho(&c, lam, y, z, pv, p, cc);
    let wabs = c_abs(&half, p);
    let kmax = rj_panel_kmax(&wabs, &rho, p);
    let wy = quad_shift_series(y, lam, &c, p);
    let wz = quad_shift_series(z, lam, &c, p);
    let wp = quad_shift_series(pv, lam, &c, p);
    let sy = inv_sqrt_of_series(&wy, p, cc, kmax);
    let sz = inv_sqrt_of_series(&wz, p, cc, kmax);
    let ip = inv_of_series(&wp, p, kmax);
    let f = conv_series(&conv_series(&sy, &sz, p), &ip, p);
    let mut acc = ExactComplex::zero(p);
    // w^{k+1}; `None` at k = 0 (half read in place).
    let mut wpow_slot: Option<ExactComplex> = None;
    for (k, fk) in f.iter().enumerate() {
        let wpow = wpow_slot.as_ref().unwrap_or(&half);
        if k % 2 == 0 {
            let term = two_sl
                .mul(fk, p, rm())
                .mul(&c_u32(2, p), p, rm())
                .mul(wpow, p, rm())
                .div(&c_u64((k as u64).saturating_add(1), p), p, rm());
            acc = acc.add(&term, p, rm());
        }
        wpow_slot = Some(wpow.mul(&half, p, rm()));
    }
    acc
}

/// ∫_0^{to} dt / (√(x+t)√(y+t)√(z+t)(p+t)) when exactly one of x,y,z is 0,
/// via t = λ u², λ = to/|to|, u from 0 to √|to|. The 1/√t singularity
/// is removed (dt/√t = 2√λ du for u ≥ 0), so panels in u are O(1).
fn rj_u2_segment(
    to: &ExactComplex,
    x: &ExactComplex,
    y: &ExactComplex,
    z: &ExactComplex,
    pv: &ExactComplex,
    p: usize,
    cc: &mut Consts,
) -> ExactComplex {
    let zeros = (is_c_zero(x), is_c_zero(y), is_c_zero(z));
    if usize::from(zeros.0) + usize::from(zeros.1) + usize::from(zeros.2) != 1 {
        return nan_pair(Error::InvalidArgument);
    }
    let (aa, bb, cc_p) = if zeros.0 {
        (y, z, pv)
    } else if zeros.1 {
        (x, z, pv)
    } else {
        (x, y, pv)
    };
    let to_abs = c_abs(to, p);
    if is_c_zero(to) || to_abs.is_zero() {
        return ExactComplex::zero(p);
    }
    let lam = to.div(&ExactComplex::from_real(to_abs.clone(), p), p, rm());
    let sqrt_lam = lam.sqrt(p, rm(), cc);
    if sqrt_lam.is_nan() {
        return nan_pair(Error::InvalidArgument);
    }
    let two_sl = c_u32(2, p).mul(&sqrt_lam, p, rm());
    let u_end = ExactComplex::from_real(to_abs, p).sqrt(p, rm(), cc);
    let mut u = ExactComplex::zero(p);
    let mut acc = ExactComplex::zero(p);
    let tiny = ExactNum::from_u8(1, p).ldexp(-((p as Exponent) + 16), p, rm());
    let mut prev_abs: Option<ExactNum> = None;
    loop {
        let rem = ExactComplex::sub_ref(&u_end, &u, p, rm());
        let rem_abs = c_abs(&rem, p);
        if matches!(rem_abs.compare(&tiny), Some(c) if c < 0) {
            return acc;
        }
        if let Some(prev) = &prev_abs {
            if rem_abs.compare(prev) != Some(-1) {
                return nan_pair(Error::InvalidArgument);
            }
        }
        prev_abs = Some(rem_abs.clone());
        let rho = rj_u2_rho(&u, &lam, aa, bb, cc_p, p, cc);
        let half_rho = rho.ldexp_ref(-1, p, rm());
        let h = if rho.is_zero() || rho.is_negative() || rem_abs.is_zero() {
            return nan_pair(Error::InvalidArgument);
        } else if matches!(rem_abs.compare(&half_rho), Some(c) if c <= 0) {
            rem
        } else {
            rem.mul(&ExactComplex::from_real(half_rho, p), p, rm()).div(
                &ExactComplex::from_real(rem_abs, p),
                p,
                rm(),
            )
        };
        acc = acc.add(
            &rj_u2_panel(&u, &h, [&lam, &two_sl], [aa, bb], cc_p, p, cc),
            p,
            rm(),
        );
        u = u.add(&h, p, rm());
    }
}

fn rj_line_integral(
    from: &ExactComplex,
    to: &ExactComplex,
    [x, y, z]: [&ExactComplex; 3],
    pv: &ExactComplex,
    p: usize,
    cc: &mut Consts,
) -> ExactComplex {
    // The panel start; `None` while it is still `from`.
    let mut a_slot: Option<ExactComplex> = None;
    let mut acc = ExactComplex::zero(p);
    // Integrable 1/√t at t=0: do the whole 0→to segment in the u²
    // coordinate so panels are O(1) instead of crawling out of 2^{-Θ(p)}.
    if is_c_zero(from) {
        let zeros = (is_c_zero(x), is_c_zero(y), is_c_zero(z));
        if usize::from(zeros.0) + usize::from(zeros.1) + usize::from(zeros.2) == 1 {
            return rj_u2_segment(to, x, y, z, pv, p, cc);
        }
    }
    let tiny = ExactNum::from_u8(1, p).ldexp(-((p as Exponent) + 16), p, rm());
    let mut prev_abs: Option<ExactNum> = None;
    loop {
        let a = a_slot.as_ref().unwrap_or(from);
        let rem = ExactComplex::sub_ref(to, a, p, rm());
        let rem_abs = c_abs(&rem, p);
        if matches!(rem_abs.compare(&tiny), Some(c) if c < 0) {
            return acc;
        }
        // Stall: the remaining distance is not strictly shrinking.
        if let Some(prev) = &prev_abs {
            if rem_abs.compare(prev) != Some(-1) {
                return nan_pair(Error::InvalidArgument);
            }
        }
        prev_abs = Some(rem_abs.clone());
        // Size the panel from the start `a`. ρ at the mid of the remaining
        // stretch can be large while the panel would pass through a nearby
        // singularity (t ≈ −y when y ≈ −1). |h| ≤ ρ(a)/2 keeps the
        // expansion disk inside the region of analyticity.
        let rho = rj_sing_dist(a, x, y, z, pv, p);
        let half_rho = rho.ldexp_ref(-1, p, rm());
        let h = if rho.is_zero() || rho.is_negative() || rem_abs.is_zero() {
            return nan_pair(Error::InvalidArgument);
        } else if matches!(rem_abs.compare(&half_rho), Some(c) if c <= 0) {
            rem
        } else {
            rem.mul(&ExactComplex::from_real(half_rho, p), p, rm()).div(
                &ExactComplex::from_real(rem_abs, p),
                p,
                rm(),
            )
        };
        acc = acc.add(&rj_panel(a, &h, [x, y, z], pv, p, cc), p, rm());
        a_slot = Some(ExactComplex::add_ref(a, &h, p, rm()));
    }
}

fn carlson_rj_dupe(
    x0: &ExactComplex,
    y0: &ExactComplex,
    z0: &ExactComplex,
    p0: &ExactComplex,
    p: usize,
    cc: &mut Consts,
) -> ExactComplex {
    let mut x = x0.clone();
    let mut y = y0.clone();
    let mut z = z0.clone();
    let mut pv = p0.clone();
    let two = c_u32(2, p);
    let five = c_u32(5, p);
    let four = c_u32(4, p);
    let six = c_u32(6, p);
    let a0 = ExactComplex::add_ref(&x, &y, p, rm())
        .add(&z, p, rm())
        .add(&two.mul(&pv, p, rm()), p, rm())
        .div(&five, p, rm());
    let delta = ExactComplex::sub_ref(&pv, &x, p, rm())
        .mul(&ExactComplex::sub_ref(&pv, &y, p, rm()), p, rm())
        .mul(&ExactComplex::sub_ref(&pv, &z, p, rm()), p, rm());
    let mut am = a0.clone();
    let mut sum = ExactComplex::zero(p);
    let mut fac = c_u32(1, p);
    let mut pow4_cubed = c_u32(1, p); // 4^{-3m}
    let real = is_real_c(&x) && is_real_c(&y) && is_real_c(&z) && is_real_c(&pv);
    let order = carlson_taylor::rj_order(p, real);
    let mut prev_d: Option<ExactNum> = None;
    let mut k = 0usize;
    loop {
        let dmax = max_abs(
            max_dev3(&am, &x, &y, &z, p),
            c_abs(&ExactComplex::sub_ref(&am, &pv, p, rm()), p),
            p,
        );
        // Quick stop: (dmax/|A|)^order < 2^{-p}. This is zf's rule, not
        // FLINT's: rf.c:241-252 uses max(|x-y|,|x-z|) / (lower bound of |x|)
        // after a duplication step. Both are scale-invariant; a floor of 1 on
        // the denominator stops too early when |A|<1 and drops trailing bits.
        // The rigorous gate is the truncation bound below (rf.c:285-295).
        let scale = c_abs(&am, p);
        let no_improve = k > 2
            && prev_d
                .as_ref()
                .is_some_and(|pd| matches!(dmax.compare(pd), Some(c) if c > 0));
        let forced = no_improve || k >= p;
        let stop_hint = carlson_taylor::rel_pow_below_prec(&dmax, &scale, order, p, p) || forced;
        let stop = stop_hint
            .then(|| {
                let t = fac.div(&am, p, rm());
                let xx = ExactComplex::sub_ref(&a0, x0, p, rm()).mul(&t, p, rm());
                let yy = ExactComplex::sub_ref(&a0, y0, p, rm()).mul(&t, p, rm());
                let zz = ExactComplex::sub_ref(&a0, z0, p, rm()).mul(&t, p, rm());
                let pp = c_neg(
                    &ExactComplex::add_ref(&xx, &yy, p, rm()).add(&zz, p, rm()),
                    p,
                )
                .div(&two, p, rm());
                (xx, yy, zz, pp)
            })
            .filter(|(xx, yy, zz, pp)| {
                // FLINT rj.c:466-477 truncation bound; below 2^-p, or keep duplicating.
                let eps = max_abs(
                    max_abs(c_abs(xx, p), c_abs(yy, p), p),
                    max_abs(c_abs(zz, p), c_abs(pp, p), p),
                    p,
                );
                forced || carlson_taylor::trunc_bound_below_prec(&eps, order, p)
            });
        if let Some((xx, yy, zz, pp)) = stop {
            let e2 = xx
                .mul(&yy, p, rm())
                .add(&xx.mul(&zz, p, rm()), p, rm())
                .add(&yy.mul(&zz, p, rm()), p, rm())
                .sub(&c_u32(3, p).mul(&pp.mul(&pp, p, rm()), p, rm()), p, rm());
            let xyz = xx.mul(&yy, p, rm()).mul(&zz, p, rm());
            let e3 = ExactComplex::add_ref(&xyz, &two.mul(&e2, p, rm()).mul(&pp, p, rm()), p, rm())
                .add(
                    &four.mul(&pp.mul(&pp, p, rm()).mul(&pp, p, rm()), p, rm()),
                    p,
                    rm(),
                );
            let e4 = two
                .mul(&xyz, p, rm())
                .add(&e2.mul(&pp, p, rm()), p, rm())
                .add(
                    &c_u32(3, p).mul(&pp.mul(&pp, p, rm()).mul(&pp, p, rm()), p, rm()),
                    p,
                    rm(),
                )
                .mul(&pp, p, rm());
            let e5 = xyz.mul(&pp.mul(&pp, p, rm()), p, rm());
            let series = carlson_taylor::rj_taylor_sum(&e2, &e3, &e4, &e5, order, p);
            // 4^{-n} A_n^{-3/2} series + 6 sum
            let ams = am.mul(&am.sqrt(p, rm(), cc), p, rm());
            return fac.mul(&series, p, rm()).div(&ams, p, rm()).add(
                &six.mul(&sum, p, rm()),
                p,
                rm(),
            );
        }
        prev_d = Some(dmax);
        let sx = x.sqrt(p, rm(), cc);
        let sy = y.sqrt(p, rm(), cc);
        let sz = z.sqrt(p, rm(), cc);
        let sp = pv.sqrt(p, rm(), cc);
        let lam = sx
            .mul(&sy, p, rm())
            .add(&sx.mul(&sz, p, rm()), p, rm())
            .add(&sy.mul(&sz, p, rm()), p, rm());
        let dm = ExactComplex::add_ref(&sp, &sx, p, rm())
            .mul(&ExactComplex::add_ref(&sp, &sy, p, rm()), p, rm())
            .mul(&ExactComplex::add_ref(&sp, &sz, p, rm()), p, rm());
        if is_c_zero(&dm) {
            return nan_pair(Error::InvalidArgument);
        }
        let em = delta
            .mul(&pow4_cubed, p, rm())
            .div(&dm.mul(&dm, p, rm()), p, rm());
        let rc = carlson_rc(&c_u32(1, p), &c_u32(1, p).add(&em, p, rm()), p, cc);
        sum = sum.add(&fac.div(&dm, p, rm()).mul(&rc, p, rm()), p, rm());
        fac = fac.div(&four, p, rm());
        pow4_cubed = pow4_cubed.div(&c_u32(64, p), p, rm()); // 4^{-3(m+1)} / 4^{-3m} = 4^{-3} = 1/64
        let xp = x.add(&lam, p, rm()).div(&four, p, rm());
        let yp = y.add(&lam, p, rm()).div(&four, p, rm());
        let zp = z.add(&lam, p, rm()).div(&four, p, rm());
        let pp = pv.add(&lam, p, rm()).div(&four, p, rm());
        am = am.add(&lam, p, rm()).div(&four, p, rm());
        x = xp;
        y = yp;
        z = zp;
        pv = pp;
        k += 1;
    }
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
    let t = carlson_norm_exp(&[x0, y0, z0, p0]);
    if t != 0 {
        // R_J(4^t x, 4^t y, 4^t z, 4^t p) = 2^-3t R_J(x, y, z, p)
        let (xs, ys) = (carlson_scaled(x0, t, p), carlson_scaled(y0, t, p));
        let (zs, ps) = (carlson_scaled(z0, t, p), carlson_scaled(p0, t, p));
        return carlson_rj(&xs, &ys, &zs, &ps, p, cc).ldexp(t.saturating_mul(-3), p, rm());
    }
    let x = leave_cut(x0, p);
    let y = leave_cut(y0, p);
    let z = leave_cut(z0, p);
    let pv = leave_cut(p0, p);
    let dupe = carlson_rj_dupe(&x, &y, &z, &pv, p, cc);
    // Principal square roots in the duplication theorem are the correct
    // sheet when every argument has nonnegative real part. Otherwise the
    // iterates can jump branch (Π(10−3i, 2+i) was all-bits-wrong) and
    // the value is recovered from the shift-and-integral identity instead.
    // Principal-sqrt duplication is always on-sheet when every argument
    // has Re ≥ 0. It can jump when one of x,y,z is 0 and another has
    // Re < 0 (complete Π at m>1: x=0, y=1-m). Other Re < 0 patterns,
    // including incomplete Π(1/2; −2 | 2) and complete Π(n, 0) with
    // Re(1-n)<0, match the principal-sqrt sheet (ZF-78 / 1a8e251).
    let has_zero_sqrt_arg = is_c_zero(&x) || is_c_zero(&y) || is_c_zero(&z);
    let sqrt_arg_neg_re = x.re().is_negative() || y.re().is_negative() || z.re().is_negative();
    let dupe_off_sheet = has_zero_sqrt_arg && sqrt_arg_neg_re;
    if !dupe_off_sheet && !dupe.is_nan() {
        return dupe;
    }
    let (n_end, mid) = shift_n(&x, &y, &z, &pv, p);
    let three_halves = c_u32(3, p).div(&c_u32(2, p), p, rm());
    let origin = ExactComplex::zero(p);
    let integ = if let Some(m) = mid {
        let i1 = rj_line_integral(&origin, &m, [&x, &y, &z], &pv, p, cc);
        let i2 = rj_line_integral(&m, &n_end, [&x, &y, &z], &pv, p, cc);
        i1.add(&i2, p, rm())
    } else {
        rj_line_integral(&origin, &n_end, [&x, &y, &z], &pv, p, cc)
    };
    let xs = x.add(&n_end, p, rm());
    let ys = y.add(&n_end, p, rm());
    let zs = z.add(&n_end, p, rm());
    let ps = pv.add(&n_end, p, rm());
    let rest = carlson_rj_dupe(&xs, &ys, &zs, &ps, p, cc);
    three_halves.mul(&integ, p, rm()).add(&rest, p, rm())
}

/// Move a real negative argument off the cut by a bump that is visible at
/// work width `p` (\(2^{-(p-16)}\), or half the width when `p` is small).
/// `force_lower` puts real \(m\ge 1\) on the lower side regardless of the
/// sign of \(\operatorname{Im} m = 0\).
fn leave_cut_visible(z: &ExactComplex, p: usize, force_lower: bool) -> ExactComplex {
    if is_real_c(z) && z.re().is_negative() {
        let shift = if p > 32 { (p as Exponent) - 16 } else { (p as Exponent) / 2 + 2 };
        let mut eps = ExactNum::from_u8(1, p).ldexp(-shift, p, rm());
        // Lower side of the m-cut is Im(m)→0−, so Im(1-m)→0+.
        if !force_lower && z.im().is_negative() {
            eps.set_sign(Sign::Neg);
        }
        ExactComplex::new(z.re().clone(), eps)
    } else {
        z.clone()
    }
}

// Design note (the dedicated function was folded into the complete-Π path; see the DLMF 19.21.6 comment there):
// Complete \(\Pi(n,m)\) when \(\operatorname{Re}(m)\ge 1\), via DLMF 19.21.6.
// Principal-sqrt \(R_J\) duplication is off-sheet for \(R_J(0,1-m,1,1-n)\),
// and the shift-integral / defining trig quadrature both miss `ZFS_WALL=120`
// at dest p≥4096. Express that \(R_J\) in \(R_F\) and \(R_D\) only:
// \((\sqrt{r}\sqrt{p}/z) R_J(0,y,z,p)=(r-1)R_F(0,y,z)R_D(p,rz,z)+R_D(0,y,z)R_F(p,rz,z)\),
// \(r=(y-p)/(y-z)\). The product of principal square roots (not \(\sqrt{rp}\))
// is the branch that matches \(R_J\) on the lower side of the \(m\)-cut.
// Real \(m\ge 1\), including \(m=2+0i\), is forced onto that lower side (Jeff / 3b5e5bd).

/// True when `z` is real and `Re(z) ≥ 1` (the Legendre cut of \(K\) and \(\Pi\)).
fn on_legendre_cut(z: &ExactComplex, p: usize) -> bool {
    if !is_real_c(z) || z.re().is_negative() {
        return false;
    }
    matches!(z.re().compare(&ExactNum::from_u8(1, p)), Some(c) if c >= 0)
}

/// Sufficient domain for the principal-sqrt DLMF 19.8.7–19.8.8 AGM,
/// pinned against Arb `acb_elliptic_pi` (python-flint 0.9.0) at 128 bits.
///
/// DLMF states the formula for real \(-\infty<k^2<1\) and \(-\infty<\alpha^2<1\)
/// (`m=k^2`, `n=\alpha^2`). The same recurrences with the principal square
/// root continue \(\Pi\) to complex `(n,m)` and match Arb's
/// \(R_F+(n/3)R_J\) branch when `n` and `m` are not on the cut `[1,+\infty)`,
/// are not the poles `n=1` or `m=1`, and `Re(n) ≤ 1`.
///
/// `Re(n) ≤ 1` is sufficient, not necessary: many points with `Re(n)>1`
/// also agree, but a thin tongue just below the `n`-cut (`Re(n)>1`,
/// `Im(n)` slightly negative) for some `m` (e.g. `m=1/2+i/2`) lands on a
/// different sheet from Arb, so those stay on the \(R_F+R_J\) path.
/// Official complete-Π remainder: the `Π(1+i, 2)` family (`m ∈ [1,+\infty)`),
/// already served by the DLMF 19.21.6 / lower-side \(R_D\) path, and the
/// incomplete row `Π(0;1/2|1/2+i/2)`.
fn pi_complete_agm_applicable(n: &ExactComplex, m: &ExactComplex, p: usize) -> bool {
    if is_c_one(n, p) || is_c_one(m, p) {
        return false;
    }
    if on_legendre_cut(n, p) || on_legendre_cut(m, p) {
        return false;
    }
    matches!(n.re().compare(&ExactNum::from_u8(1, p)), Some(c) if c <= 0)
}

/// DLMF 19.8.7–19.8.8, principal square roots.
/// Π(n, m) = π / (4 M(1, √(1−m))) · (2 + n/(1−n) Σ Q_j),
/// a₀=1, g₀=√(1−m), p₀=√(1−n), Q₀=1,
/// p_{j+1}=(p_j² + a_j g_j)/(2 p_j), ε_j=(p_j² − a_j g_j)/(p_j² + a_j g_j),
/// Q_{j+1}=ε_j Q_j / 2.
fn pi_complete_agm(
    n: &ExactComplex,
    m: &ExactComplex,
    p: usize,
    _dest: usize,
    cc: &mut Consts,
) -> ExactComplex {
    let two = c_u32(2, p);
    let mut a = ExactComplex::one(p);
    let om = one_minus(m, p);
    let on = one_minus(n, p);
    if is_c_zero(&on) {
        return nan_pair(Error::InvalidArgument);
    }
    let mut g = om.sqrt(p, rm(), cc);
    let mut pn = on.sqrt(p, rm(), cc);
    if is_c_zero(&pn) || g.is_nan() || pn.is_nan() {
        return nan_pair(Error::InvalidArgument);
    }
    let mut q = ExactComplex::one(p);
    let mut qsum = ExactComplex::one(p);
    let thresh = ExactNum::from_u8(1, p).ldexp(-((p as Exponent) + 32), p, rm());
    // Quadratic AGM: about log2(p) iterations.
    for _ in 0..=p {
        let ag = a.mul(&g, p, rm());
        let pn2 = pn.mul(&pn, p, rm());
        let den = ExactComplex::add_ref(&pn2, &ag, p, rm());
        if is_c_zero(&den) {
            return nan_pair(Error::InvalidArgument);
        }
        let eps = pn2.sub(&ag, p, rm()).div(&den, p, rm());
        q = eps.mul(&q, p, rm()).div(&two, p, rm());
        qsum = qsum.add(&q, p, rm());
        let an = a.add(&g, p, rm()).div(&two, p, rm());
        g = ag.sqrt(p, rm(), cc);
        if is_c_zero(&pn) {
            return nan_pair(Error::InvalidArgument);
        }
        pn = den.div(&two, p, rm()).div(&pn, p, rm());
        a = an;
        let spread = c_abs(&ExactComplex::sub_ref(&a, &g, p, rm()), p);
        let qabs = c_abs(&q, p);
        if matches!(spread.compare(&thresh), Some(c) if c < 0)
            && matches!(qabs.compare(&thresh), Some(c) if c < 0)
        {
            break;
        }
    }
    // M = a; Π = π/(4M) * (2 + n/(1−n) * qsum)
    let four = c_u32(4, p);
    let pi = ExactComplex::from_real(cc.pi(p, rm()), p);
    let pref = pi.div(&four, p, rm()).div(&a, p, rm());
    let n_over = n.div(&on, p, rm());
    let inner = two.add(&n_over.mul(&qsum, p, rm()), p, rm());
    let val = pref.mul(&inner, p, rm());
    if is_real_c(n) && is_real_c(m) {
        return ExactComplex::from_real(val.into_parts().0, p);
    }
    val
}

fn pi_complete_via_rd(
    n: &ExactComplex,
    m: &ExactComplex,
    p: usize,
    cc: &mut Consts,
) -> ExactComplex {
    let zero = ExactComplex::zero(p);
    let three = c_u32(3, p);
    // leave_cut's bump 2^{-(dest+16)} is below the ulp at work width p≈dest,
    // so principal-sqrt R_F(0,1-m,1) landed on the upper side of K(2)
    // (1.311+1.311i). Use a bump visible at width p; Ziv tightens the rest.
    // Real m≥1, including m=2+0i, is forced onto the lower side (Jeff / 3b5e5bd).
    let y = leave_cut_visible(&one_minus(m, p), p, is_real_c(m) && !m.re().is_negative());
    let z = ExactComplex::one(p);
    let pv = leave_cut_visible(&one_minus(n, p), p, false);
    let yz = ExactComplex::sub_ref(&y, &z, p, rm());
    if is_c_zero(&yz) {
        return nan_pair(Error::InvalidArgument);
    }
    let r = ExactComplex::sub_ref(&y, &pv, p, rm()).div(&yz, p, rm());
    let fac = r
        .sqrt(p, rm(), cc)
        .mul(&pv.sqrt(p, rm(), cc), p, rm())
        .div(&z, p, rm());
    if is_c_zero(&fac) || fac.is_nan() {
        return nan_pair(Error::InvalidArgument);
    }
    let rf0 = carlson_rf(&zero, &y, &z, p, cc);
    let rd0 = carlson_rd(&zero, &y, &z, p, cc);
    let rz = r.mul(&z, p, rm());
    let rd1 = carlson_rd(&pv, &rz, &z, p, cc);
    let rf1 = carlson_rf(&pv, &rz, &z, p, cc);
    let right = r
        .sub(&z, p, rm())
        .mul(&rf0, p, rm())
        .mul(&rd1, p, rm())
        .add(&rd0.mul(&rf1, p, rm()), p, rm());
    let rj = right.div(&fac, p, rm());
    rf0.add(&n.div(&three, p, rm()).mul(&rj, p, rm()), p, rm())
}

fn pi_complete_via_rf_rj(
    n: &ExactComplex,
    m: &ExactComplex,
    p: usize,
    dest: usize,
    cc: &mut Consts,
    on_cut: bool,
) -> ExactComplex {
    let wp = if on_cut { rj_cut_work(dest, p) } else { p };
    let zero = ExactComplex::zero(wp);
    let one = ExactComplex::one(wp);
    let three = c_u32(3, wp);
    let om = leave_cut(&one_minus(m, wp), dest);
    let on = leave_cut(&one_minus(n, wp), dest);
    let rf = carlson_rf(&zero, &om, &one, wp, cc);
    let rj = carlson_rj(&zero, &om, &one, &on, wp, cc);
    let v = rf.add(&n.div(&three, wp, rm()).mul(&rj, wp, rm()), wp, rm());
    if on_cut {
        flush_cut_bias(v, dest)
    } else {
        v
    }
}

impl ExactComplex {
    /// Complete elliptic \(K(m)\), \(m=k^2\). Cut on \([1,+\infty)\). \(m=1\) is \(+\infty\).
    ///
    /// # Precision
    ///
    /// - Algorithm: Carlson `R_F` in ℂ; duplication until the series remainder bound.
    /// - Bound: identities evaluated outside Ziv (nested Ziv would exhaust `MAX_PREC_RETRY`).
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn elliptic_k(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.elliptic_k_word(p, rm, cc))
    }

    /// Implementation of [`Self::elliptic_k`] on operands without inexact history.
    fn elliptic_k_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return typed_nan_pair(self);
        }
        let dest = round_p(p);
        if is_c_one(self, dest) {
            return crate::complex::compact_inf(false, false);
        }
        if let Some(n) = crate::complex_special::nan_if_dest_unallocatable(dest) {
            return n;
        }
        if is_c_zero(self) {
            // π/2 + 0i: π rounded once at (dest, rm), halved exactly. Rounding π at dest and
            // then dividing at rm rounded twice (K(0) Up was 1 ulp low, ZF-18).
            let hp = cc.pi(dest, rm).ldexp(-1, dest, rm);
            return ExactComplex::from_real(hp, dest);
        }
        ziv_complex(dest, rm, |pw| self.elliptic_k_at(pw, dest, cc))
    }

    fn elliptic_k_at(&self, p: usize, dest: usize, cc: &mut Consts) -> Self {
        let om_tmp = one_minus(self, p);
        let on_cut = is_real_c(&om_tmp) && om_tmp.re().is_negative();
        // R_F / R_D only: no R_C(1,1+e_m) cancellation, so no extra work.
        let wp = p;
        let zero = ExactComplex::zero(wp);
        let one = ExactComplex::one(wp);
        let om = leave_cut(&one_minus(self, wp), dest);
        let v = carlson_rf(&zero, &om, &one, wp, cc);
        if on_cut {
            flush_cut_bias(v, dest)
        } else {
            v
        }
    }

    /// Complete elliptic \(E(m)\). \(E(1)=1\). Cut of \(K\) inherited through \(1-m\).
    /// At exact \(+\infty+0i\) the value is compact \(+0-i\infty\) (lower side of
    /// the \(m>1\) cut, same sign as finite real \(m>1\) with imag \(+0\)); at
    /// \(+\infty-0i\) it is \(+0+i\infty\); at \(-\infty+0i\) it is compact \(+\infty+0i\).
    ///
    /// # Precision
    ///
    /// - Algorithm: Carlson `R_F` / `R_D`; duplication until the series remainder bound.
    /// - Bound: same as [`Self::elliptic_k`].
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
    pub fn elliptic_e_complete(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let inputs = self.inexact_accumulated();
        Self::per_op(inputs, this.elliptic_e_complete_word(p, rm, cc))
    }

    /// Implementation of [`Self::elliptic_e_complete`] on operands without inexact history.
    fn elliptic_e_complete_word(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if self.is_nan() {
            return typed_nan_pair(self);
        }
        let dest = round_p(p);
        if self.im().is_zero() && self.re().is_inf() {
            if self.re().is_negative() {
                return crate::complex::compact_inf(false, false);
            }
            // Finite real m>1 already returns +0 · π/(4√m) with Im of opposite
            // sign to the signed zero (lower side at +0, matching the m>1 cut).
            let im = if self.im().is_negative() { crate::INF_POS } else { crate::INF_NEG };
            return ExactComplex::new(crate::complex::compact_signed_zero(false), im);
        }
        if is_c_one(self, dest) {
            return crate::complex::compact_one();
        }
        if let Some(n) = crate::complex_special::nan_if_dest_unallocatable(dest) {
            return n;
        }
        if is_c_zero(self) {
            // π/2 + 0i: π rounded once at (dest, rm), halved exactly. Rounding π at dest and
            // then dividing at rm rounded twice (K(0) Up was 1 ulp low, ZF-18).
            let hp = cc.pi(dest, rm).ldexp(-1, dest, rm);
            return ExactComplex::from_real(hp, dest);
        }
        ziv_complex(dest, rm, |pw| self.elliptic_e_complete_at(pw, dest, cc))
    }

    fn elliptic_e_complete_at(&self, p: usize, dest: usize, cc: &mut Consts) -> Self {
        let om_tmp = one_minus(self, p);
        let on_cut = is_real_c(&om_tmp) && om_tmp.re().is_negative();
        // R_F / R_D only: no R_C(1,1+e_m) cancellation, so no extra work.
        let wp = p;
        let zero = ExactComplex::zero(wp);
        let one = ExactComplex::one(wp);
        let three = c_u32(3, wp);
        let om = leave_cut(&one_minus(self, wp), dest);
        let rf = carlson_rf(&zero, &om, &one, wp, cc);
        let rd = carlson_rd(&zero, &om, &one, wp, cc);
        let v = rf.sub(&self.div(&three, wp, rm()).mul(&rd, wp, rm()), wp, rm());
        if on_cut {
            flush_cut_bias(v, dest)
        } else {
            v
        }
    }

    /// Incomplete \(F(x|m)\). `self` is \(x=\sin\varphi\), not the amplitude
    /// \(\varphi\). \(m=k^2\), not the modulus \(k\).
    ///
    /// # Precision
    ///
    /// - Algorithm: Carlson `R_F`; duplication until the series remainder bound.
    /// - MPFR oracle: no.
    ///
    /// Carlson \(R_F(1-x^2,1-mx^2,1)\). Cuts when \(1-x^2\) or \(1-mx^2\) lies on
    /// \((-\infty,0]\) (principal square-root cut). \(F(x,0)=\arcsin x\).
    ///
    /// Accuracy: correctly-rounded
    pub fn elliptic_f(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let m_c = m.exact_copy();
        let inputs = self.inexact_accumulated() || m.inexact_accumulated();
        Self::per_op(inputs, this.elliptic_f_word(&m_c, p, rm, cc))
    }

    /// Implementation of [`Self::elliptic_f`] on operands without inexact history.
    fn elliptic_f_word(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
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
        ziv_complex(dest, rm, |pw| self.elliptic_f_at(m, pw, dest, cc))
    }

    fn elliptic_f_at(&self, m: &Self, p: usize, dest: usize, cc: &mut Consts) -> Self {
        let x2_tmp = self.mul(self, p, rm());
        let a_tmp = one_minus(&x2_tmp, p);
        let b_tmp = one_minus(&m.mul(&x2_tmp, p, rm()), p);
        let cut_a = on_nonpos_real(&a_tmp);
        let cut_b = on_nonpos_real(&b_tmp);
        let on_cut = cut_a || cut_b;
        // R_F / R_D only: no R_C(1,1+e_m) cancellation, so no extra work.
        let wp = p;
        // If 1−x² is on the cut, only bump x: that first-order imag
        // already places 1−m x² (and 1−n x²) on the same side. Bumping
        // m or n as well can overpower x. If 1−x² is off the cut, bump
        // m (and n) when their product term is on the cut.
        let x_bumped;
        let x: &ExactComplex = if cut_a && is_real_c(self) {
            x_bumped = bump_real(self, dest);
            &x_bumped
        } else {
            self
        };
        let mm_bumped;
        let mm: &ExactComplex = if !cut_a && cut_b && is_real_c(m) {
            mm_bumped = bump_real(m, dest);
            &mm_bumped
        } else {
            m
        };
        let one = ExactComplex::one(wp);
        let x2 = x.mul(x, wp, rm());
        let a = leave_cut(&one_minus(&x2, wp), dest);
        let b = leave_cut(&one_minus(&mm.mul(&x2, wp, rm()), wp), dest);
        let rf = carlson_rf(&a, &b, &one, wp, cc);
        let v = self.mul(&rf, wp, rm());
        if on_cut {
            flush_cut_bias(v, dest)
        } else {
            v
        }
    }

    /// Incomplete \(E(x|m)\). Same convention as [`Self::elliptic_f`]:
    /// `self` is \(x=\sin\varphi\), not \(\varphi\), and \(m=k^2\), not \(k\).
    ///
    /// # Precision
    ///
    /// - Algorithm: Carlson `R_F` / `R_D`; duplication until the series remainder bound.
    /// - MPFR oracle: no.
    ///
    /// Cuts as for \(F\). \(E(x,0)=\arcsin x\); \(E(x,1)=x\).
    ///
    /// Accuracy: correctly-rounded
    pub fn elliptic_e(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        let this = self.exact_copy();
        let m_c = m.exact_copy();
        let inputs = self.inexact_accumulated() || m.inexact_accumulated();
        Self::per_op(inputs, this.elliptic_e_word(&m_c, p, rm, cc))
    }

    /// Implementation of [`Self::elliptic_e`] on operands without inexact history.
    fn elliptic_e_word(&self, m: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
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
        ziv_complex(dest, rm, |pw| self.elliptic_e_at(m, pw, dest, cc))
    }

    fn elliptic_e_at(&self, m: &Self, p: usize, dest: usize, cc: &mut Consts) -> Self {
        let x2_tmp = self.mul(self, p, rm());
        let a_tmp = one_minus(&x2_tmp, p);
        let b_tmp = one_minus(&m.mul(&x2_tmp, p, rm()), p);
        let cut_a = on_nonpos_real(&a_tmp);
        let cut_b = on_nonpos_real(&b_tmp);
        let on_cut = cut_a || cut_b;
        // R_F / R_D only: no R_C(1,1+e_m) cancellation, so no extra work.
        let wp = p;
        let x_bumped;
        let x: &ExactComplex = if cut_a && is_real_c(self) {
            x_bumped = bump_real(self, dest);
            &x_bumped
        } else {
            self
        };
        let mm_bumped;
        let mm: &ExactComplex = if !cut_a && cut_b && is_real_c(m) {
            mm_bumped = bump_real(m, dest);
            &mm_bumped
        } else {
            m
        };
        let one = ExactComplex::one(wp);
        let three = c_u32(3, wp);
        let x2 = x.mul(x, wp, rm());
        let a = leave_cut(&one_minus(&x2, wp), dest);
        let b = leave_cut(&one_minus(&mm.mul(&x2, wp, rm()), wp), dest);
        let rf = carlson_rf(&a, &b, &one, wp, cc);
        let rd = carlson_rd(&a, &b, &one, wp, cc);
        // Prefactors use the original (unbumped) x and m; only the Carlson
        // arguments need the first-order imag.
        let x2_pref = self.mul(self, wp, rm());
        let v = self.mul(&rf, wp, rm()).sub(
            &m.mul(self, wp, rm())
                .mul(&x2_pref, wp, rm())
                .div(&three, wp, rm())
                .mul(&rd, wp, rm()),
            wp,
            rm(),
        );
        if on_cut {
            flush_cut_bias(v, dest)
        } else {
            v
        }
    }

    /// Complete \(\Pi(n,m)\). `self` is \(n\). \(\Pi(0,m)=K(m)\). Pole at \(n=1\).
    ///
    /// # Precision
    ///
    /// - Algorithm: Carlson `R_J`; duplication until the series remainder bound.
    /// - MPFR oracle: no.
    ///
    /// Accuracy: correctly-rounded
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
        // Π(n,m) ~ −1/n ∫ dθ /(sin²θ √(1−m sin²θ)) → 0 as |n|→∞ (mpmath
        // ellippi(±∞, m)=0). Without this, R_J(0,1−m,1,1−∞) hangs past
        // ZFS_WALL at dest p=8192, and the Re(m)≥1 path returns NaN
        // (Harper T4 S2, n=+∞).
        if self.re().is_inf() || self.im().is_inf() {
            return ExactComplex::zero(dest);
        }
        // DLMF 19.6.4: as m→1−, Π(n,m) ∼ [1/(2(1−n))] log 1/(1−m)
        // diverges in the direction of 1/(1−n). After |n|=∞ → 0, so
        // Π(∞, 1) stays 0. n=1 is already InvalidArgument.
        if is_c_one(m, dest) {
            return pi_complete_m1_infinity(self);
        }
        ziv_complex(dest, rm, |pw| self.elliptic_pi_complete_at(m, pw, dest, cc))
    }

    fn elliptic_pi_complete_at(&self, m: &Self, p: usize, dest: usize, cc: &mut Consts) -> Self {
        // DLMF 19.6.3 / 19.2.8: Π(n, 0) = π / (2 √(1−n)).
        if is_c_zero(m) {
            let on = one_minus(self, p);
            let s = on.sqrt(p, rm(), cc);
            if s.is_nan() || is_c_zero(&s) {
                return nan_pair(Error::InvalidArgument);
            }
            let half_pi = ExactComplex::from_real(cc.pi(p, rm()), p).div(&c_u32(2, p), p, rm());
            return half_pi.div(&s, p, rm());
        }
        // DLMF 19.8.7–19.8.8 AGM, principal square roots, in the
        // sufficient complex domain pinned against Arb (see
        // `pi_complete_agm_applicable`). About log2(p) iterations
        // instead of O(p) R_J duplications. Official remainder
        // (m on [1,+∞), incomplete Π) keeps the paths below.
        if pi_complete_agm_applicable(self, m, p) {
            return pi_complete_agm(self, m, p, dest, cc);
        }
        let om_tmp = one_minus(m, p);
        let on_tmp = one_minus(self, p);
        let on_cut = (is_real_c(&om_tmp) && om_tmp.re().is_negative())
            || (is_real_c(&on_tmp) && on_tmp.re().is_negative());
        // When Re(1-m) ≤ 0 (Re(m) ≥ 1), principal-sqrt R_J duplication
        // is off-sheet and the shift-integral / trig quadrature miss
        // ZFS_WALL at dest p≥4096. DLMF 19.21.6 writes that R_J in
        // R_F and R_D, which stay fast at dest p=8192.
        if !om_tmp.re().is_positive() {
            let v = pi_complete_via_rd(self, m, p, cc);
            return if on_cut { flush_cut_bias(v, dest) } else { v };
        }
        pi_complete_via_rf_rj(self, m, p, dest, cc, on_cut)
    }

    /// Incomplete \(\Pi(n;x|m)\). `self` is the characteristic \(n\).
    /// `x` is \(\sin\varphi\), not the amplitude \(\varphi\). \(m=k^2\), not \(k\).
    ///
    /// # Precision
    ///
    /// - Algorithm: Carlson `R_J`; duplication until the series remainder bound.
    /// - MPFR oracle: no.
    ///
    /// Cuts when \(1-x^2\), \(1-mx^2\), or \(1-nx^2\) meets the Carlson cut
    /// \((-\infty,0]\). \(\Pi(0;x|m)=F(x|m)\).
    ///
    /// Accuracy: correctly-rounded
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
        ziv_complex(dest, rm, |pw| self.elliptic_pi_at(x, m, pw, dest, cc))
    }

    fn elliptic_pi_at(&self, x: &Self, m: &Self, p: usize, dest: usize, cc: &mut Consts) -> Self {
        let x2_tmp = x.mul(x, p, rm());
        let a_tmp = one_minus(&x2_tmp, p);
        let b_tmp = one_minus(&m.mul(&x2_tmp, p, rm()), p);
        let pv_tmp = one_minus(&self.mul(&x2_tmp, p, rm()), p);
        let cut_a = on_nonpos_real(&a_tmp);
        let cut_b = on_nonpos_real(&b_tmp);
        let cut_n = on_nonpos_real(&pv_tmp);
        let on_cut = cut_a || cut_b || cut_n;
        let wp = if on_cut { rj_cut_work(dest, p) } else { p };
        let xb_bumped;
        let xb: &ExactComplex = if cut_a && is_real_c(x) {
            xb_bumped = bump_real(x, dest);
            &xb_bumped
        } else {
            x
        };
        let mm_bumped;
        let mm: &ExactComplex = if !cut_a && cut_b && is_real_c(m) {
            mm_bumped = bump_real(m, dest);
            &mm_bumped
        } else {
            m
        };
        let nn_bumped;
        let nn: &ExactComplex = if !cut_a && cut_n && is_real_c(self) {
            nn_bumped = bump_real(self, dest);
            &nn_bumped
        } else {
            self
        };
        let one = ExactComplex::one(wp);
        let three = c_u32(3, wp);
        let x2 = xb.mul(xb, wp, rm());
        let a = leave_cut(&one_minus(&x2, wp), dest);
        let b = leave_cut(&one_minus(&mm.mul(&x2, wp, rm()), wp), dest);
        let pv = leave_cut(&one_minus(&nn.mul(&x2, wp, rm()), wp), dest);
        let rf = carlson_rf(&a, &b, &one, wp, cc);
        let rj = carlson_rj(&a, &b, &one, &pv, wp, cc);
        let x2_pref = x.mul(x, wp, rm());
        let v = x.mul(&rf, wp, rm()).add(
            &self
                .mul(x, wp, rm())
                .mul(&x2_pref, wp, rm())
                .div(&three, wp, rm())
                .mul(&rj, wp, rm()),
            wp,
            rm(),
        );
        if on_cut {
            flush_cut_bias(v, dest)
        } else {
            v
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_must::Must;
    use crate::Sign;
    use crate::INF_NEG;
    use crate::INF_POS;
    use core::ops::Neg;
    use core::str::FromStr;

    fn near(a: &ExactNum, b: &ExactNum, p: usize, slack: Exponent) -> bool {
        let d = ExactNum::sub_ref(a, b, p, RoundingMode::None).abs();
        d.is_zero() || d.exponent().unwrap_or(0) < -((p as Exponent) - slack)
    }

    fn tiny(x: &ExactNum, p: usize) -> bool {
        x.is_zero() || x.exponent().unwrap_or(0) < -((p as Exponent) / 4)
    }

    #[test]
    fn complex_elliptic_golds() -> Result<(), crate::Error> {
        let p = 256;
        let r = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let z0 = ExactComplex::zero(p);
        let mk_half = || ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, r);
        let half = mk_half();
        let mh = ExactComplex::from_real(mk_half(), p);

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

        let x = ExactComplex::from_real(mk_half(), p);
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
        let mk_eps = || ExactNum::from_u8(1, p).ldexp(-40, p, RoundingMode::None);
        let eps = mk_eps();
        let above = ExactComplex::new(ExactNum::from_u8(2, p), mk_eps());
        let below = ExactComplex::new(two, mk_eps().neg());
        let ka = above.elliptic_k(p, r, &mut cc);
        let kb = below.elliptic_k(p, r, &mut cc);
        assert!(near(ka.re(), kb.re(), p, 20), "K cut Re");
        assert!(ka.im().is_positive() != kb.im().is_positive() || !tiny(ka.im(), p));
        assert_ne!(ka.im().compare(kb.im()), Some(0));
        assert!(
            near(&ka.into_parts().1.abs(), &kb.into_parts().1.abs(), p, 20),
            "K cut |Im|"
        );

        // Incomplete F cut: x=3/2, m=1/2, ±ε imag on x.
        let three_half = ExactNum::from_u8(3, p).div(&ExactNum::from_u8(2, p), p, r);
        let xa = ExactComplex::new(
            ExactNum::from_u8(3, p).div(&ExactNum::from_u8(2, p), p, r),
            mk_eps(),
        );
        let xb = ExactComplex::new(three_half, eps.neg());
        let fa = xa.elliptic_f(&mh, p, r, &mut cc);
        let fb = xb.elliptic_f(&mh, p, r, &mut cc);
        assert!(!fa.is_nan() && !fb.is_nan(), "F cut defined");
        assert_ne!(fa.im().compare(fb.im()), Some(0), "F cut Im differs");
        Ok(())
    }

    /// ZF-66: Π(1+i, 0) = π/(2√(1-(1+i))). The old R_J duplication
    /// (R_C of α=[p(√x+√y+√z)+√x√y√z]²) was all-bits-wrong here.
    #[test]
    fn zf66_pi_complete_closed_form() -> Result<(), crate::Error> {
        let p = 256;
        let r = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let n = ExactComplex::new(ExactNum::from_u8(1, p), ExactNum::from_u8(1, p));
        let m = ExactComplex::zero(p);
        let got = n.elliptic_pi_complete(&m, p, r, &mut cc);
        let on = ExactComplex::one(p).sub(&n, p, r);
        let want = ExactComplex::from_real(cc.pi(p, r), p)
            .div(&c_u32(2, p), p, r)
            .div(&on.sqrt(p, r, &mut cc), p, r);
        assert!(
            near(got.re(), want.re(), p, 8),
            "Π(1+i,0) re got {:?} want {:?}",
            got.re(),
            want.re()
        );
        assert!(
            near(got.im(), want.im(), p, 8),
            "Π(1+i,0) im got {:?} want {:?}",
            got.im(),
            want.im()
        );
        Ok(())
    }

    /// Official-family values against Arb midpoints at dest 64, 256, and 1024.
    #[test]
    fn zf_pi_family_matches_arb_midpoints() -> Result<(), crate::Error> {
        let r = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let check = |p: usize,
                     got: ExactComplex,
                     re: &str,
                     im: &str,
                     label: &str|
         -> Result<(), crate::Error> {
            let wr = ExactNum::from_str(re)?;
            let wi = ExactNum::from_str(im)?;
            // Golds are 34-digit midpoints; compare in a 64-bit window.
            assert!(
                near(got.re(), &wr, 64, 8) && near(got.im(), &wi, 64, 8),
                "{label} p={p} got {:?} + {:?}i want {re} + {im}i",
                got.re(),
                got.im()
            );
            Ok(())
        };
        for p in [64usize, 256, 1024] {
            let n = ExactComplex::from_real(ExactNum::from_u8(2, p).neg(), p);
            let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, r);
            let m = ExactComplex::from_real(half, p);
            let got = n.elliptic_pi_complete(&m, p, r, &mut cc);
            check(
                p,
                got,
                "1.0226975086065923912210632881204494",
                "0",
                "Pi(-2,1/2)",
            )?;

            let n = ExactComplex::new(ExactNum::from_u8(1, p), ExactNum::from_u8(1, p));
            let m = ExactComplex::from_real(ExactNum::from_u8(1, p).neg(), p);
            let got = n.elliptic_pi_complete(&m, p, r, &mut cc);
            check(
                p,
                got,
                "0.9894414586482923407135834361683811",
                "0.8576009178956446444392334130940548",
                "Pi(1+i,-1)",
            )?;

            let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, r);
            let n = ExactComplex::new(ExactNum::from_u8(1, p), ExactNum::from_u8(1, p));
            let m = ExactComplex::new(
                ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, r),
                ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, r),
            );
            let got = n.elliptic_pi_complete(&m, p, r, &mut cc);
            check(
                p,
                got,
                "0.8493037644403146015738020272566064",
                "1.3478039001871343637887723094489923",
                "Pi(1+i,1/2+i/2)",
            )?;

            let n = ExactComplex::zero(p);
            let x = ExactComplex::from_real(
                ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, r),
                p,
            );
            let m = ExactComplex::new(
                ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, r),
                half,
            );
            let got = n.elliptic_pi(&x, &m, p, r, &mut cc);
            check(
                p,
                got,
                "0.5348066320602228343305056081235160",
                "0.0127136117099412408682670775585325",
                "inc Pi(0;1/2|1/2+i/2)",
            )?;
        }
        Ok(())
    }

    /// Π(−2, 1/2) is in the real DLMF 19.8.7 domain. The AGM path must
    /// agree with dest 64 and finish quickly at dest 16384.
    #[test]
    fn zf_pi_complete_real_agm_n_m2_m_half() -> Result<(), crate::Error> {
        let r = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let p64 = 64usize;
        let n64 = ExactComplex::from_real(ExactNum::from_u8(2, p64).neg(), p64);
        let half64 = ExactNum::from_u8(1, p64).div(&ExactNum::from_u8(2, p64), p64, r);
        let m64 = ExactComplex::from_real(half64, p64);
        let low = n64.elliptic_pi_complete(&m64, p64, r, &mut cc);
        assert!(!low.is_nan(), "Π(-2,1/2) p=64 NaN");
        for p in [256usize, 1024, 8192] {
            let n = ExactComplex::from_real(ExactNum::from_u8(2, p).neg(), p);
            let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, r);
            let m = ExactComplex::from_real(half, p);
            let got = n.elliptic_pi_complete(&m, p, r, &mut cc);
            assert!(!got.is_nan(), "Π(-2,1/2) p={p} NaN");
            assert!(
                near(got.re(), low.re(), p64, 8) && near(got.im(), low.im(), p64, 8),
                "Π(-2,1/2) p={p} disagree with dest 64"
            );
        }
        let p = 16384usize;
        let t0 = std::time::Instant::now();
        let n = ExactComplex::from_real(ExactNum::from_u8(2, p).neg(), p);
        let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, r);
        let m = ExactComplex::from_real(half, p);
        let got = n.elliptic_pi_complete(&m, p, r, &mut cc);
        let elapsed = t0.elapsed();
        assert!(!got.is_nan(), "Π(-2,1/2) p=16384 NaN");
        assert!(
            near(got.re(), low.re(), p64, 8) && near(got.im(), low.im(), p64, 8),
            "Π(-2,1/2) p=16384 disagree with dest 64"
        );
        assert!(
            crate::test_time::within(elapsed, std::time::Duration::from_secs(5)),
            "Π(-2,1/2) p=16384 took {elapsed:?}"
        );
        Ok(())
    }

    /// Complex DLMF 19.8.7–19.8.8 AGM matches the existing R_F+(n/3)R_J
    /// path at dest 64/256/1024 on the official in-domain complete-Π rows,
    /// including both tight rows.
    #[test]
    fn zf_pi_complete_agm_matches_rf_rj_official_rows() -> Result<(), crate::Error> {
        let r = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        for p in [64usize, 256, 1024] {
            let one = || ExactNum::from_u8(1, p);
            let half = || ExactNum::div_ref(&one(), &ExactNum::from_u8(2, p), p, r);
            let cases: &[(&str, ExactComplex, ExactComplex)] = &[
                (
                    "Pi(-2,1/2)",
                    ExactComplex::from_real(ExactNum::from_u8(2, p).neg(), p),
                    ExactComplex::from_real(half(), p),
                ),
                (
                    "Pi(1+i,-1)",
                    ExactComplex::new(ExactNum::from_u8(1, p), ExactNum::from_u8(1, p)),
                    ExactComplex::from_real(ExactNum::from_u8(1, p).neg(), p),
                ),
                (
                    "Pi(1+i,0)",
                    ExactComplex::new(ExactNum::from_u8(1, p), ExactNum::from_u8(1, p)),
                    ExactComplex::zero(p),
                ),
                (
                    "Pi(1+i,1/2+i/2)",
                    ExactComplex::new(ExactNum::from_u8(1, p), ExactNum::from_u8(1, p)),
                    ExactComplex::new(half(), half()),
                ),
            ];
            for (lab, n, m) in cases {
                assert!(
                    pi_complete_agm_applicable(n, m, p),
                    "{lab} p={p} should be in the AGM domain"
                );
                let agm = pi_complete_agm(n, m, p, p, &mut cc);
                let rj = pi_complete_via_rf_rj(n, m, p, p, &mut cc, false);
                assert!(!agm.is_nan() && !rj.is_nan(), "{lab} p={p} NaN");
                assert!(
                    near(agm.re(), rj.re(), p, 16) && near(agm.im(), rj.im(), p, 16),
                    "{lab} p={p} AGM {:?} RJ {:?}",
                    agm,
                    rj
                );
            }
            // Official remainder: m=2 stays off the AGM path.
            let n = ExactComplex::new(one(), one());
            let m2 = ExactComplex::from_real(ExactNum::from_u8(2, p), p);
            assert!(
                !pi_complete_agm_applicable(&n, &m2, p),
                "Pi(1+i,2) must stay on the R_D / lower-side path"
            );
        }
        Ok(())
    }

    /// A box around both tight official rows at dest 64: every in-domain
    /// point of the AGM path matches R_F+R_J.
    #[test]
    fn zf_pi_complete_agm_grid_around_tight_rows() -> Result<(), crate::Error> {
        let r = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let p = 64usize;
        let one = ExactNum::from_u8(1, p);
        let two = ExactNum::from_u8(2, p);
        let four = ExactNum::from_u8(4, p);
        let half = || ExactNum::div_ref(&one, &two, p, r);
        let qtr = || ExactNum::div_ref(&one, &four, p, r);
        let centers: [(ExactComplex, ExactComplex, &str); 2] = [
            (
                ExactComplex::new(ExactNum::from_u8(1, p), ExactNum::from_u8(1, p)),
                ExactComplex::from_real(ExactNum::from_u8(1, p).neg(), p),
                "Pi(1+i,-1)",
            ),
            (
                ExactComplex::new(ExactNum::from_u8(1, p), ExactNum::from_u8(1, p)),
                ExactComplex::new(half(), half()),
                "Pi(1+i,1/2+i/2)",
            ),
        ];
        // Offsets ±0, ±1/4, ±1/2 in Re and Im of n and (separately) of m.
        let offs = [ExactNum::from_u8(0, p), qtr(), qtr().neg(), half(), half().neg()];
        let mut checked = 0usize;
        for (n0, m0, lab) in &centers {
            for dre in &offs {
                for dim in &offs {
                    let n = ExactComplex::new(
                        ExactNum::add_ref(n0.re(), dre, p, r),
                        ExactNum::add_ref(n0.im(), dim, p, r),
                    );
                    let m = m0;
                    if !pi_complete_agm_applicable(&n, m, p) {
                        continue;
                    }
                    let agm = pi_complete_agm(&n, m, p, p, &mut cc);
                    let rj = pi_complete_via_rf_rj(&n, m, p, p, &mut cc, false);
                    assert!(
                        !agm.is_nan()
                            && !rj.is_nan()
                            && near(agm.re(), rj.re(), p, 16)
                            && near(agm.im(), rj.im(), p, 16),
                        "{lab} n-shift dre/dim AGM {:?} RJ {:?}",
                        agm,
                        rj
                    );
                    checked += 1;
                }
            }
            for dre in &offs {
                for dim in &offs {
                    let n = n0;
                    let m = ExactComplex::new(
                        ExactNum::add_ref(m0.re(), dre, p, r),
                        ExactNum::add_ref(m0.im(), dim, p, r),
                    );
                    if !pi_complete_agm_applicable(n, &m, p) {
                        continue;
                    }
                    let agm = pi_complete_agm(n, &m, p, p, &mut cc);
                    let rj = pi_complete_via_rf_rj(n, &m, p, p, &mut cc, false);
                    assert!(
                        !agm.is_nan()
                            && !rj.is_nan()
                            && near(agm.re(), rj.re(), p, 16)
                            && near(agm.im(), rj.im(), p, 16),
                        "{lab} m-shift dre/dim AGM {:?} RJ {:?}",
                        agm,
                        rj
                    );
                    checked += 1;
                }
            }
        }
        assert!(checked >= 40, "grid too small: checked={checked}");
        Ok(())
    }

    /// After the AGM dispatch, both tight official rows at dest 8192 must
    /// finish in well under a second (the R_J path was a few seconds).
    #[test]
    fn zf_pi_complete_agm_tight_rows_p8192_fast() -> Result<(), crate::Error> {
        let r = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let p = 8192usize;
        let one = ExactNum::from_u8(1, p);
        let n = ExactComplex::new(ExactNum::from_u8(1, p), ExactNum::from_u8(1, p));
        let m1 = ExactComplex::from_real(ExactNum::from_u8(1, p).neg(), p);
        let half = || ExactNum::div_ref(&one, &ExactNum::from_u8(2, p), p, r);
        let m2 = ExactComplex::new(half(), half());
        for (lab, m) in [("Pi(1+i,-1)", &m1), ("Pi(1+i,1/2+i/2)", &m2)] {
            let t0 = std::time::Instant::now();
            let got = n.elliptic_pi_complete(m, p, r, &mut cc);
            let dt = t0.elapsed();
            assert!(!got.is_nan(), "{lab} p=8192 NaN");
            assert!(
                crate::test_time::within(dt, std::time::Duration::from_secs(2)),
                "{lab} p=8192 took {dt:?}"
            );
        }
        Ok(())
    }

    /// Official S5 dest 65536 rows Π(−2, 0) and Π(1+i, 0) must use the
    /// closed form and finish immediately, agreeing with dest 64.
    #[test]
    fn zf_pi_complete_m0_at_p65536_is_instant() -> Result<(), crate::Error> {
        let r = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let p64 = 64usize;
        let p = 65536usize;
        let t0 = std::time::Instant::now();
        let n1 = ExactComplex::from_real(ExactNum::from_u8(2, p).neg(), p);
        let m0 = ExactComplex::zero(p);
        let g1 = n1.elliptic_pi_complete(&m0, p, r, &mut cc);
        let n2 = ExactComplex::new(ExactNum::from_u8(1, p), ExactNum::from_u8(1, p));
        let g2 = n2.elliptic_pi_complete(&m0, p, r, &mut cc);
        let elapsed = t0.elapsed();
        assert!(!g1.is_nan() && !g2.is_nan(), "Π(*,0) p=65536 NaN");
        assert!(
            crate::test_time::within(elapsed, std::time::Duration::from_secs(2)),
            "Π(n,0) p=65536 took {elapsed:?}"
        );
        let l1 = ExactComplex::from_real(ExactNum::from_u8(2, p64).neg(), p64)
            .elliptic_pi_complete(&ExactComplex::zero(p64), p64, r, &mut cc);
        let l2 = ExactComplex::new(ExactNum::from_u8(1, p64), ExactNum::from_u8(1, p64))
            .elliptic_pi_complete(&ExactComplex::zero(p64), p64, r, &mut cc);
        assert!(near(g1.re(), l1.re(), p64, 8) && near(g1.im(), l1.im(), p64, 8));
        assert!(near(g2.re(), l2.re(), p64, 8) && near(g2.im(), l2.im(), p64, 8));
        Ok(())
    }

    fn signed_im(neg: bool, p: usize) -> ExactNum {
        let mut z = ExactNum::new(p);
        if neg {
            z.set_sign(Sign::Neg);
        }
        z
    }

    fn m_on_cut(neg_im: bool, p: usize) -> ExactComplex {
        ExactComplex::new(ExactNum::from_u8(2, p), signed_im(neg_im, p))
    }

    /// ZF-66: m on [1,+∞) is a cut. Complete Π at m=2+0i uses the lower
    /// side (Im(m)→0−), like Arb, the same value as m=2−0i
    /// (≈2.8818+0i for n=1+i). Official R_J at exact y=1-m=-1 is
    /// all-bits-wrong; the trigonometric integral takes the lower branch.
    #[test]
    fn zf66_pi_complete_m_cut_signed_zeros() -> Result<(), crate::Error> {
        let p = 64;
        let r = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let n = ExactComplex::new(ExactNum::from_u8(1, p), ExactNum::from_u8(1, p));
        let plus = n.elliptic_pi_complete(&m_on_cut(false, p), p, r, &mut cc);
        let minus = n.elliptic_pi_complete(&m_on_cut(true, p), p, r, &mut cc);
        let want_re = ExactNum::from_str("2.8818251039409565")?;
        assert!(
            near(minus.re(), &want_re, p, 16),
            "Π(1+i, 2-0i) re got {:?} want {:?}",
            minus.re(),
            want_re
        );
        assert!(
            tiny(minus.im(), p),
            "Π(1+i, 2-0i) im should be ~0, got {:?}",
            minus.im()
        );
        assert!(
            near(plus.re(), &want_re, p, 16),
            "Π(1+i, 2+0i) lower-side re got {:?} want {:?}",
            plus.re(),
            want_re
        );
        assert!(
            tiny(plus.im(), p),
            "Π(1+i, 2+0i) lower-side im should be ~0, got {:?}",
            plus.im()
        );
        Ok(())
    }

    /// Π(10−3i, 2+0i) uses the lower side of the cut, like Arb
    /// (≈ −0.032260744 − 0.375062852i), not the C99 0+ limit. Off-cut
    /// m=2+2^{-10}i is the un-extrapolated Arb value used to isolate
    /// leave-cut vs. R_J.
    #[test]
    fn zf66_pi_complete_n10_minus_3i_m_two_plus_zero() -> Result<(), crate::Error> {
        let p = 64;
        let r = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let n = ExactComplex::new(ExactNum::from_u8(10, p), ExactNum::from_u8(3, p).neg());
        let m_off = ExactComplex::new(
            ExactNum::from_u8(2, p),
            ExactNum::from_u8(1, p).ldexp(-10, p, RoundingMode::None),
        );
        let off = n.elliptic_pi_complete(&m_off, p, r, &mut cc);
        let want_off_re = ExactNum::from_str("0.11066294559569917")?;
        let want_off_im = ExactNum::from_str("-0.77654340189875511")?;
        assert!(!off.is_nan(), "Π(10-3i, 2+2^{{-10}}i) is NaN");
        assert!(
            near(off.re(), &want_off_re, p, 8),
            "off-cut re got {:?} want {:?}",
            off.re(),
            want_off_re
        );
        assert!(
            near(off.im(), &want_off_im, p, 8),
            "off-cut im got {:?} want {:?}",
            off.im(),
            want_off_im
        );

        let got = n.elliptic_pi_complete(&m_on_cut(false, p), p, r, &mut cc);
        let want_re = ExactNum::from_str("-0.03226074424035563377376688819209160641")?;
        let want_im = ExactNum::from_str("-0.37506285238448681531122676224843104364")?;
        assert!(!got.is_nan(), "Π(10-3i, 2+0i) is NaN");
        assert!(
            near(got.re(), &want_re, p, 8),
            "lower-side re got {:?} want {:?}",
            got.re(),
            want_re
        );
        assert!(
            near(got.im(), &want_im, p, 8),
            "lower-side im got {:?} want {:?}",
            got.im(),
            want_im
        );

        let p128 = 128;
        let n128 = ExactComplex::new(
            ExactNum::from_u8(10, p128),
            ExactNum::from_u8(3, p128).neg(),
        );
        let got128 = n128.elliptic_pi_complete(&m_on_cut(false, p128), p128, r, &mut cc);
        assert!(!got128.is_nan(), "Π(10-3i, 2+0i) p=128 is NaN {:?}", got128);
        let want128_re = ExactNum::from_str("-0.0322607442403556337750171956587330382258")?;
        let want128_im = ExactNum::from_str("-0.3750628523844868153150179966420258747532")?;
        assert!(
            near(got128.re(), &want128_re, p128, 8),
            "p=128 lower-side re got {:?} want {:?}",
            got128.re(),
            want128_re
        );
        assert!(
            near(got128.im(), &want128_im, p128, 8),
            "p=128 lower-side im got {:?} want {:?}",
            got128.im(),
            want128_im
        );
        Ok(())
    }

    /// Jeff 2026-11-18 6:31 AM CT: Π(n, 2+0i) uses the lower side of the cut,
    /// like Arb, not the C99 0+ sheet. Π(1+i, 2+0i) must match the 0− value
    /// previously computed at m=2−0i.
    #[test]
    fn zf_pi_complete_m_two_plus_zero_uses_lower_side() -> Result<(), crate::Error> {
        let p = 64;
        let r = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let n = ExactComplex::new(ExactNum::from_u8(1, p), ExactNum::from_u8(1, p));
        let got = n.elliptic_pi_complete(&m_on_cut(false, p), p, r, &mut cc);
        let want_re = ExactNum::from_str("2.8818251039409565")?;
        assert!(!got.is_nan(), "Π(1+i, 2+0i) is NaN");
        assert!(
            near(got.re(), &want_re, p, 16),
            "lower-side re got {:?} want {:?}",
            got.re(),
            want_re
        );
        assert!(
            tiny(got.im(), p),
            "lower-side im should be ~0, got {:?}",
            got.im()
        );
        Ok(())
    }

    /// Harper T4 S2: complete Π(1+i, 2+0i) at dest p=4096 timed out under
    /// ZFS_WALL=120 on the defining trigonometric integral. Must finish and
    /// match the lower-side value (Jeff / 3b5e5bd): ≈ 2.8818 + 0i.
    #[test]
    fn zf_elliptic_pi_complete_m2_at_p4096_finishes() -> Result<(), crate::Error> {
        let r = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let want_re = ExactNum::from_str("2.8818251039409565")?;
        for p in [64usize, 4096, 8192] {
            let n = ExactComplex::new(ExactNum::from_u8(1, p), ExactNum::from_u8(1, p));
            let m = ExactComplex::from_real(ExactNum::from_u8(2, p), p);
            let got = n.elliptic_pi_complete(&m, p, r, &mut cc);
            assert!(!got.is_nan(), "Π(1+i, 2+0i) p={p} is NaN {:?}", got);
            assert!(
                near(got.re(), &want_re, 64, 16),
                "p={p} lower-side re got {:?} want {:?}",
                got.re(),
                want_re
            );
            assert!(
                tiny(got.im(), p),
                "p={p} lower-side im should be ~0, got {:?}",
                got.im()
            );
        }
        Ok(())
    }

    /// Harper T4 S2: complete Π(+∞, 1/2+i/2) at dest p=8192 timed out under
    /// ZFS_WALL=120 in R_J(0,1−m,1,−∞). The same n at m=2 returned NaN.
    /// The defining integral vanishes as |n|→∞ (mpmath ellippi(±∞,m)=0).
    #[test]
    fn zf_elliptic_pi_complete_n_inf_is_zero() -> Result<(), crate::Error> {
        let r = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        for p in [64usize, 8192] {
            let n = ExactComplex::from_real(INF_POS, p);
            let m_half = ExactComplex::new(
                ExactNum::from_u8(1, p).ldexp(-1, p, RoundingMode::None),
                ExactNum::from_u8(1, p).ldexp(-1, p, RoundingMode::None),
            );
            let m2 = ExactComplex::from_real(ExactNum::from_u8(2, p), p);
            for (m, label) in [(&m_half, "1/2+i/2"), (&m2, "2+0i")] {
                let got = n.elliptic_pi_complete(m, p, r, &mut cc);
                assert!(!got.is_nan(), "Π(+∞, {label}) p={p} is NaN {:?}", got);
                assert!(got.re().is_zero(), "Π(+∞, {label}) p={p} re {:?}", got.re());
                assert!(got.im().is_zero(), "Π(+∞, {label}) p={p} im {:?}", got.im());
            }
        }
        Ok(())
    }

    /// Complete Π(n, 1) has the log singularity of DLMF 19.6.4, coefficient
    /// 1/(1−n). The via_rd path used to return InvalidArgument at exact m=1.
    #[test]
    fn zf_elliptic_pi_complete_m1_is_signed_infinity() -> Result<(), crate::Error> {
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
        // Each part: None = compact zero, Some(true)=+∞, Some(false)=−∞.
        type Case = (
            &'static str,
            fn(usize) -> ExactComplex,
            Option<bool>,
            Option<bool>,
        );
        let cases: &[Case] = &[
            // remaining-28 n = 1/2+3/4 i → +∞ + i∞
            (
                "1/2+3/4i",
                |p| {
                    let half = ExactNum::from_u8(1, p).ldexp(-1, p, RoundingMode::None);
                    let three_four = ExactNum::from_u8(3, p).ldexp(-2, p, RoundingMode::None);
                    ExactComplex::new(half, three_four)
                },
                Some(true),
                Some(true),
            ),
            // remaining-28 n = 10−3i → −∞ − i∞
            (
                "10-3i",
                |p| ExactComplex::new(ExactNum::from_u8(10, p), ExactNum::from_u8(3, p).neg()),
                Some(false),
                Some(false),
            ),
            // remaining-28 n = 2^{-60}+i → +∞ + i∞
            (
                "2^{-60}+i",
                |p| {
                    ExactComplex::new(
                        ExactNum::from_u8(1, p).ldexp(-60, p, RoundingMode::None),
                        ExactNum::from_u8(1, p),
                    )
                },
                Some(true),
                Some(true),
            ),
            // real n < 1 → +∞ + 0i (DLMF 19.6.4, sign(1-n)=+)
            (
                "1/2",
                |p| {
                    ExactComplex::from_real(
                        ExactNum::from_u8(1, p).ldexp(-1, p, RoundingMode::None),
                        p,
                    )
                },
                Some(true),
                None,
            ),
            (
                "-1",
                |p| ExactComplex::from_real(ExactNum::from_i8(-1, p), p),
                Some(true),
                None,
            ),
            // real n > 1 → −∞ + 0i
            (
                "2",
                |p| ExactComplex::from_real(ExactNum::from_u8(2, p), p),
                Some(false),
                None,
            ),
            // n=0 is Π(0,1)=K(1)=+∞ + 0i
            ("0", |p| ExactComplex::zero(p), Some(true), None),
            // Re n = 1, Im n ≠ 0: 1/(1-n) is purely imaginary → 0 + i∞ sign(Im n)
            (
                "1+i",
                |p| ExactComplex::new(ExactNum::from_u8(1, p), ExactNum::from_u8(1, p)),
                None,
                Some(true),
            ),
            (
                "1-i",
                |p| ExactComplex::new(ExactNum::from_u8(1, p), ExactNum::from_u8(1, p).neg()),
                None,
                Some(false),
            ),
        ];
        for p in [64usize, 256, 8192] {
            let m1 = ExactComplex::from_real(ExactNum::from_u8(1, p), p);
            for rm in modes {
                for &(label, mk, re_pos, im_pos) in cases {
                    let n = mk(p);
                    let got = n.elliptic_pi_complete(&m1, p, rm, &mut cc);
                    assert!(!got.is_nan(), "Π({label}, 1) p={p} {rm:?} is NaN {:?}", got);
                    match re_pos {
                        None => assert!(
                            got.re().is_zero(),
                            "Π({label}, 1) p={p} {rm:?} re {:?}",
                            got.re()
                        ),
                        Some(true) => assert!(
                            got.re().is_inf_pos(),
                            "Π({label}, 1) p={p} {rm:?} re {:?}",
                            got.re()
                        ),
                        Some(false) => assert!(
                            got.re().is_inf_neg(),
                            "Π({label}, 1) p={p} {rm:?} re {:?}",
                            got.re()
                        ),
                    }
                    match im_pos {
                        None => assert!(
                            got.im().is_zero(),
                            "Π({label}, 1) p={p} {rm:?} im {:?}",
                            got.im()
                        ),
                        Some(true) => assert!(
                            got.im().is_inf_pos(),
                            "Π({label}, 1) p={p} {rm:?} im {:?}",
                            got.im()
                        ),
                        Some(false) => assert!(
                            got.im().is_inf_neg(),
                            "Π({label}, 1) p={p} {rm:?} im {:?}",
                            got.im()
                        ),
                    }
                }
                // n=1 stays a pole
                let n1 = ExactComplex::from_real(ExactNum::from_u8(1, p), p);
                let pole = n1.elliptic_pi_complete(&m1, p, rm, &mut cc);
                assert!(
                    pole.is_nan(),
                    "Π(1, 1) p={p} {rm:?} should stay IA, got {:?}",
                    pole
                );
                // n=+∞ and n=−∞ at m=1 stay 0 (existing |n|→∞ path)
                for (ninf, lab) in [
                    (ExactComplex::from_real(INF_POS, p), "+∞"),
                    (ExactComplex::from_real(INF_NEG, p), "−∞"),
                ] {
                    let got = ninf.elliptic_pi_complete(&m1, p, rm, &mut cc);
                    assert!(!got.is_nan(), "Π({lab}, 1) p={p} {rm:?} is NaN");
                    assert!(
                        got.re().is_zero() && got.im().is_zero(),
                        "Π({lab}, 1) p={p} {rm:?} {:?}",
                        got
                    );
                }
                // neighbor m=1−2^{-(p+8)} is finite, not a compact infinity and not IA
                let eps = ExactNum::from_u8(1, p).ldexp(-(p as i64 + 8), p, RoundingMode::None);
                let m_near = ExactComplex::from_real(
                    ExactNum::from_u8(1, p).sub(&eps, p, RoundingMode::None),
                    p,
                );
                let n_half = ExactComplex::from_real(
                    ExactNum::from_u8(1, p).ldexp(-1, p, RoundingMode::None),
                    p,
                );
                let neigh = n_half.elliptic_pi_complete(&m_near, p, rm, &mut cc);
                assert!(!neigh.is_nan(), "Π(1/2, 1-eps) p={p} {rm:?} is NaN");
                assert!(
                    !neigh.re().is_inf(),
                    "Π(1/2, 1-eps) p={p} {rm:?} should be finite re {:?}",
                    neigh.re()
                );
            }
        }
        Ok(())
    }

    /// Harper T4 S2: complete Π(1+i, 0) at dest p=8192 timed out under ZFS_WALL=120.
    /// Each R_J increment called a full R_F for R_C. Closed form: Π(n,0)=π/(2√(1−n)).
    /// For n=1+i this is π√2 (1+i)/4.
    #[test]
    fn zf_elliptic_pi_complete_modest_n_at_p8192_finishes() -> Result<(), crate::Error> {
        let p = 8192;
        let r = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let n = ExactComplex::new(ExactNum::from_u8(1, p), ExactNum::from_u8(1, p));
        let got = n.elliptic_pi_complete(&ExactComplex::zero(p), p, r, &mut cc);
        assert!(!got.is_nan(), "Π(1+i, 0) is NaN {:?}", got);
        let q = cc
            .pi(p, RoundingMode::None)
            .mul(&ExactNum::from_u8(2, p).sqrt(p, r), p, r)
            .div(&ExactNum::from_u8(4, p), p, r);
        assert!(
            near(got.re(), &q, p, 8),
            "Π(1+i, 0) re got {:?} want {:?}",
            got.re(),
            q
        );
        assert!(
            near(got.im(), &q, p, 8),
            "Π(1+i, 0) im got {:?} want {:?}",
            got.im(),
            q
        );
        Ok(())
    }

    fn assert_dest_oom_fast(name: &str, z: ExactComplex, start: std::time::Instant) {
        assert!(
            crate::test_time::within(start.elapsed(), std::time::Duration::from_secs(2)),
            "{name} at dest 2^59 took {:?}",
            start.elapsed()
        );
        assert!(z.is_nan(), "{name} not NaN: {z:?}");
        assert_eq!(
            z.re().err(),
            Some(Error::MemoryAllocation),
            "{name} re err {:?}",
            z.re().err()
        );
        assert_eq!(
            z.im().err(),
            Some(Error::MemoryAllocation),
            "{name} im err {:?}",
            z.im().err()
        );
    }

    #[test]
    fn zf_elliptic_k_e_dest_2_59_fail_fast_memory_allocation() -> Result<(), String> {
        let p = 1usize << 59;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new().must("Consts")?;
        let z0 = ExactComplex::new(
            crate::complex::compact_signed_zero(false),
            crate::complex::compact_signed_zero(true),
        );
        let t0 = std::time::Instant::now();
        let k = z0.elliptic_k(p, rm, &mut cc);
        assert_dest_oom_fast("K(0-0i)", k, t0);
        let t1 = std::time::Instant::now();
        let e = z0.elliptic_e_complete(p, rm, &mut cc);
        assert_dest_oom_fast("E(0-0i)", e, t1);

        // Exact specials still return first (no dest allocation).
        let one = crate::complex::compact_one();
        let t2 = std::time::Instant::now();
        let k1 = one.elliptic_k(p, rm, &mut cc);
        assert!(
            crate::test_time::within(t2.elapsed(), std::time::Duration::from_secs(1)),
            "K(1) took {:?}",
            t2.elapsed()
        );
        assert!(k1.re().is_inf_pos() && k1.im().is_zero(), "K(1) got {k1:?}");
        let t3 = std::time::Instant::now();
        let e1 = one.elliptic_e_complete(p, rm, &mut cc);
        assert!(
            crate::test_time::within(t3.elapsed(), std::time::Duration::from_secs(1)),
            "E(1) took {:?}",
            t3.elapsed()
        );
        assert_eq!(e1.re().compare(&ExactNum::from_u8(1, 64)), Some(0));
        assert!(e1.im().is_zero(), "E(1) im {:?}", e1.im());
        Ok(())
    }

    /// Compact limits of complete E and of F/E(x,0)=asin(x) at exact ±∞.
    /// E(+∞+0i)=+0−i∞ (lower side of the m>1 cut); E(−∞+0i)=+∞+0i; F/E(+∞+0i, 0)=π/2+i∞ (DLMF 4.23.20, upper side).
    /// The via-Ziv / inv_axis path used to return InvalidArgument at every infinity.
    #[test]
    fn zf_elliptic_e_complete_and_fe_at_infinity() -> Result<(), crate::Error> {
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
        for p in [64usize, 256, 2048, 8192] {
            let z0 = ExactComplex::zero(p);
            let one = ExactComplex::from_real(ExactNum::from_u8(1, p), p);
            for rm in modes {
                let epos = ExactComplex::from_real(INF_POS, p).elliptic_e_complete(p, rm, &mut cc);
                assert!(
                    !epos.is_nan()
                        && epos.re().is_zero()
                        && !epos.re().is_negative()
                        && epos.im().is_inf_neg(),
                    "E(+∞+0i) p={p} {rm:?} {:?}",
                    epos
                );
                let epos_m = ExactComplex::new(INF_POS, crate::complex::compact_signed_zero(true))
                    .elliptic_e_complete(p, rm, &mut cc);
                assert!(
                    !epos_m.is_nan()
                        && epos_m.re().is_zero()
                        && !epos_m.re().is_negative()
                        && epos_m.im().is_inf_pos(),
                    "E(+∞-0i) p={p} {rm:?} {:?}",
                    epos_m
                );
                let eneg = ExactComplex::from_real(INF_NEG, p).elliptic_e_complete(p, rm, &mut cc);
                assert!(
                    !eneg.is_nan() && eneg.re().is_inf_pos() && eneg.im().is_zero(),
                    "E(−∞) p={p} {rm:?} {:?}",
                    eneg
                );
                let hp = cc.pi(p, rm).ldexp(-1, p, RoundingMode::None);
                let hp_m = {
                    let rm_m = match rm {
                        RoundingMode::Up => RoundingMode::Down,
                        RoundingMode::Down => RoundingMode::Up,
                        other => other,
                    };
                    cc.pi(p, rm_m).ldexp(-1, p, RoundingMode::None).neg()
                };
                let cases = [
                    (
                        "+∞+0i",
                        INF_POS,
                        crate::complex::compact_signed_zero(false),
                        &hp,
                        true,
                    ),
                    (
                        "+∞-0i",
                        INF_POS,
                        crate::complex::compact_signed_zero(true),
                        &hp,
                        false,
                    ),
                    (
                        "-∞+0i",
                        INF_NEG,
                        crate::complex::compact_signed_zero(false),
                        &hp_m,
                        true,
                    ),
                    (
                        "-∞-0i",
                        INF_NEG,
                        crate::complex::compact_signed_zero(true),
                        &hp_m,
                        false,
                    ),
                ];
                for (lab, re, im, want_re, im_pos) in cases {
                    let x = ExactComplex::new(re, im);
                    for (fnlab, got) in [
                        ("asin", x.asin(p, rm, &mut cc)),
                        ("F", x.elliptic_f(&z0, p, rm, &mut cc)),
                        ("E", x.elliptic_e(&z0, p, rm, &mut cc)),
                    ] {
                        assert!(!got.is_nan(), "{fnlab}({lab}, 0) p={p} {rm:?} is NaN");
                        assert!(
                            got.re().compare(want_re) == Some(0),
                            "{fnlab}({lab}, 0) p={p} {rm:?} re {:?} want {:?}",
                            got.re(),
                            want_re
                        );
                        if im_pos {
                            assert!(
                                got.im().is_inf_pos(),
                                "{fnlab}({lab}, 0) p={p} {rm:?} im {:?}",
                                got.im()
                            );
                        } else {
                            assert!(
                                got.im().is_inf_neg(),
                                "{fnlab}({lab}, 0) p={p} {rm:?} im {:?}",
                                got.im()
                            );
                        }
                    }
                }
                // W4: E(∞, 1) stays +∞+0i
                let ew4 = ExactComplex::from_real(INF_POS, p).elliptic_e(&one, p, rm, &mut cc);
                assert!(
                    !ew4.is_nan() && ew4.re().is_inf_pos() && ew4.im().is_zero(),
                    "E(+∞, 1) p={p} {rm:?} should stay +∞+0i, got {:?}",
                    ew4
                );
            }
            if p > 256 {
                continue;
            }
            let hundred = ExactComplex::from_real(ExactNum::from_u8(100, p), p);
            let e100 = hundred.elliptic_e_complete(p, RoundingMode::ToEven, &mut cc);
            assert!(
                !e100.is_nan()
                    && !e100.re().is_inf()
                    && !e100.im().is_inf()
                    && e100.re().is_positive()
                    && e100.im().is_negative(),
                "E(100+0i) should be a finite + − i (lower side), got {:?}",
                e100
            );
        }
        Ok(())
    }

    /// R2 (Spence PR #40, 2026-11-22): boundary tests for the complex Carlson `|A|` scale and truncation.
    ///
    /// FLINT 3.6.0 golds (`complex_elliptic_r2_golds.rs`). Each row runs at p = 64, 128, 256, 1024, 8192 and is
    /// checked normwise: `|got - want| <= 2^(e - p + SLACK)` with `e` the exponent of the larger part of `want`.
    /// The `*_b` and `*_d` rows also run with every argument scaled by `4^k`, and are compared with
    /// `2^-k want` (R_F) or `2^-3k want` (R_D, R_J), the homogeneity of DLMF 19.16 (R_F has degree -1/2, R_D
    /// and R_J -3/2). `k = 2^60` makes `|A|^2` overflow (`c_abs(A) = +Inf` while the deviations stay finite);
    /// `k = -2^60` and `k = 1 - 2^61` make `|A|^2` underflow (`c_abs(A) = 0`).
    mod r2_carlson_boundaries {
        use super::*;
        use crate::common::util::round_p;
        use crate::ExactInt;

        include!("complex_elliptic_r2_golds.rs");

        const SLACK: Exponent = 24;
        const PRECS: [usize; 5] = [64, 128, 256, 1024, 8192];

        fn gold_num(m: &str, e: i64) -> Result<ExactNum, String> {
            if m == "0" {
                return Ok(ExactNum::new(64));
            }
            let i = ExactInt::from_dec_string(m)
                .ok_or_else(|| alloc::format!("bad gold mantissa {m}"))?;
            let bits = round_p(i.bit_length().max(64));
            Ok(i.to_exact_num(bits, RoundingMode::None)
                .ldexp(e, bits, RoundingMode::None))
        }

        fn gold(name: &str, p: usize) -> Result<ExactComplex, String> {
            let row = R2_GOLDS
                .iter()
                .find(|r| r.0 == name && r.1 == p)
                .ok_or_else(|| alloc::format!("no gold {name} p={p}"))?;
            Ok(ExactComplex::new(
                gold_num(row.2 .0, row.2 .1)?,
                gold_num(row.3 .0, row.3 .1)?,
            ))
        }

        /// `n 2^e`, exact.
        fn d(n: i64, e: Exponent, w: usize) -> ExactNum {
            let m = ExactNum::from_word(n.unsigned_abs() as crate::Word, w).ldexp(
                e,
                w,
                RoundingMode::None,
            );
            if n < 0 {
                m.neg()
            } else {
                m
            }
        }

        fn c(re: ExactNum, im: ExactNum) -> ExactComplex {
            ExactComplex::new(re, im)
        }

        /// The arguments of each row, exact dyadics at width `w` (the generator's rows, same order).
        fn args(name: &str, w: usize) -> Result<Vec<ExactComplex>, String> {
            let z = || ExactNum::new(w);
            let r = |n: i64, e: Exponent| c(d(n, e, w), ExactNum::new(w));
            let one_m =
                |k: i64| ExactNum::sub_ref(&d(1, 0, w), &d(k, -40, w), w, RoundingMode::None);
            let one_p = |k: i64| d(1, 0, w).add(&d(k, -40, w), w, RoundingMode::None);
            Ok(match name {
                "rf_a0" => {
                    alloc::vec![r(1, 0), c(d(-1, -1, w), d(1, 0, w)), c(d(-1, -1, w), d(-1, 0, w))]
                }
                "rf_atiny" => alloc::vec![
                    r(1, 0),
                    c(d(-1, -1, w), d(1, 0, w)),
                    c(
                        d(-1, -1, w).add(&d(3, -200, w), w, RoundingMode::None),
                        d(-1, 0, w)
                    ),
                ],
                "rf_b" | "rd_b" => alloc::vec![r(2, 0), r(1, 0), r(1, 0)],
                "rf_below" => alloc::vec![r(1, 0), r(1, 0), c(one_m(3), z())],
                "rf_above" => alloc::vec![r(1, 0), r(1, 0), c(one_p(3), z())],
                "rf_d" | "rd_d" => {
                    alloc::vec![c(d(1, 0, w), d(1, 0, w)), c(d(2, 0, w), d(-1, 0, w)), r(3, 0)]
                }
                "rd_a0" => {
                    alloc::vec![c(d(1, 0, w), d(1, 0, w)), c(d(2, 0, w), d(-1, 0, w)), r(-1, 0)]
                }
                "rd_below" => alloc::vec![c(one_m(5), z()), r(1, 0), r(1, 0)],
                "rd_above" => alloc::vec![c(one_p(5), z()), r(1, 0), r(1, 0)],
                "rj_a0" => alloc::vec![
                    r(1, 0),
                    c(d(1, 0, w), d(1, 0, w)),
                    c(d(1, 0, w), d(-1, 0, w)),
                    r(-3, -1)
                ],
                "rj_b" => alloc::vec![r(2, 0), r(1, 0), r(1, 0), r(1, 0)],
                "rj_below" => alloc::vec![c(one_m(5), z()), r(1, 0), r(1, 0), r(1, 0)],
                "rj_above" => alloc::vec![c(one_p(5), z()), r(1, 0), r(1, 0), r(1, 0)],
                "rj_d" => alloc::vec![
                    c(d(1, 0, w), d(1, 0, w)),
                    c(d(2, 0, w), d(-1, 0, w)),
                    r(3, 0),
                    c(d(1, 0, w), d(2, 0, w)),
                ],
                _ => return Err(alloc::format!("no row {name}")),
            })
        }

        fn eval(name: &str, a: &[ExactComplex], p: usize, cc: &mut Consts) -> ExactComplex {
            match &name[..2] {
                "rf" => carlson_rf(&a[0], &a[1], &a[2], p, cc),
                "rd" => carlson_rd(&a[0], &a[1], &a[2], p, cc),
                _ => carlson_rj(&a[0], &a[1], &a[2], &a[3], p, cc),
            }
        }

        /// Normwise `|got - want| <= 2^(e - p + SLACK)`, `e` the larger part exponent of `want`.
        fn close(got: &ExactComplex, want: &ExactComplex, p: usize) -> bool {
            let e = [want.re(), want.im()]
                .iter()
                .filter(|x| !x.is_zero())
                .filter_map(|x| x.exponent())
                .max();
            let Some(e) = e else {
                return got.re().is_zero() && got.im().is_zero();
            };
            let w = round_p(p + 128);
            let lim = e - (p as Exponent) + SLACK;
            [(got.re(), want.re()), (got.im(), want.im())]
                .iter()
                .all(|(g, t)| {
                    let dd = ExactNum::sub_ref(g, t, w, RoundingMode::None);
                    !dd.is_nan() && (dd.is_zero() || dd.exponent().is_some_and(|de| de <= lim))
                })
        }

        const ROWS: [&str; 16] = [
            "rf_a0", "rf_atiny", "rf_b", "rf_below", "rf_above", "rf_d", "rd_a0", "rd_b",
            "rd_below", "rd_above", "rd_d", "rj_a0", "rj_b", "rj_below", "rj_above", "rj_d",
        ];

        #[test]
        fn r2_carlson_rows_match_flint_at_each_tier() -> Result<(), String> {
            let mut cc = Consts::new().map_err(|e| alloc::format!("Consts: {e:?}"))?;
            let mut bad = Vec::new();
            for &p in &PRECS {
                for name in ROWS {
                    let got = eval(name, &args(name, p + 64)?, p, &mut cc);
                    let want = gold(name, p)?;
                    if !close(&got, &want, p) {
                        bad.push(alloc::format!("{name} p={p}: got {got:?} want {want:?}"));
                    }
                }
            }
            assert!(
                bad.is_empty(),
                "{} of 80 rows wrong:\n{}",
                bad.len(),
                bad.join("\n")
            );
            Ok(())
        }

        /// FLINT's truncation bound `2 q^order / (1 - q)`, `q = 9 eps / 8` (rf.c:285-295), against hand values:
        /// `eps = 2^-4`, `order = 10` gives `log2 = 1 + 10 log2(9/128) - log2(119/128) = -37.196...`
        /// (mpmath 1.3.0), so it is below `2^-37` and not below `2^-38`. The quick rule stops at prec 40 for
        /// `eps = 2^-4 (1 - 2^-10)` (`eps < 2^-ceil(40/10)`), where the bound is still about `2^-37.2`: the
        /// bound is the stricter gate.
        #[test]
        fn r2_truncation_bound_matches_flint_formula() {
            let tb = carlson_taylor::trunc_bound_below_prec;
            let e4 = ExactNum::from_u8(1, 64).ldexp(-4, 64, RoundingMode::None);
            assert!(tb(&e4, 10, 37));
            assert!(!tb(&e4, 10, 38));
            let e4m = ExactNum::from_word(1023, 64).ldexp(-14, 64, RoundingMode::None);
            assert!(carlson_taylor::rel_pow_below_prec(
                &e4m,
                &ExactNum::from_u8(1, 64),
                10,
                40,
                64
            ));
            assert!(!tb(&e4m, 10, 40));
            assert!(
                tb(&ExactNum::new(64), 10, 1 << 40),
                "eps = 0: nothing is left out"
            );
            assert!(!tb(&ExactNum::nan(None), 10, 1), "NaN eps never stops");
            assert!(!tb(&crate::INF_POS, 10, 1), "Inf eps never stops");
            let q1 = ExactNum::from_u8(8, 64).div(&ExactNum::from_u8(9, 64), 64, RoundingMode::Up);
            assert!(!tb(&q1, 500, 1), "q >= 1: the geometric series diverges");
            assert!(!tb(&ExactNum::from_u8(1, 64), 500, 1));
        }

        /// `4^k` scaling: `|A|` huge (`k = 2^60`, `c_abs(A)` overflows), tiny (`k = -2^60`, and `1 - 2^61`
        /// for R_F, where the arguments sit two binades above `EXPONENT_MIN`), and ordinary `k`.
        #[test]
        fn r2_carlson_scaled_rows_follow_homogeneity() -> Result<(), String> {
            let mut cc = Consts::new().map_err(|e| alloc::format!("Consts: {e:?}"))?;
            let big = 1i64 << 60;
            let (mut bad, mut n) = (Vec::new(), 0usize);
            for &p in &PRECS {
                for name in ["rf_b", "rf_d", "rd_b", "rd_d", "rj_b", "rj_d"] {
                    let deg: i64 = if name.starts_with("rf") { 1 } else { 3 };
                    let mut ks = alloc::vec![1i64, -1, 1000, -1000, big, -big];
                    if deg == 1 {
                        ks.push(1 - (1i64 << 61));
                    }
                    for k in ks {
                        let a: Vec<ExactComplex> = args(name, p + 64)?
                            .iter()
                            .map(|x| x.ldexp(2 * k, p + 64, RoundingMode::None))
                            .collect();
                        let got = eval(name, &a, p, &mut cc);
                        let want =
                            gold(name, p)?.ldexp(-deg * k, round_p(p + 128), RoundingMode::None);
                        n += 1;
                        if !close(&got, &want, p) {
                            bad.push(alloc::format!(
                                "{name} 4^{k} p={p}: got {got:?} want {want:?}"
                            ));
                        }
                    }
                }
            }
            assert!(
                bad.is_empty(),
                "{} of {n} scaled rows wrong:\n{}",
                bad.len(),
                bad.join("\n")
            );
            Ok(())
        }
    }
}
