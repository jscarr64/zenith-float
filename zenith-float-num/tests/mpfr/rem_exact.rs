//! `rem` is exact: the truncated remainder `x - trunc(x/y)·y` is representable at the
//! precision of `x`, so it needs no rounding mode and is never inexact. Checked against an
//! exact big-integer remainder (2.0.0 audit, round 2 leftovers).

use crate::mpfr::common::{conv_to_mpfr, draw_seed, rebuild_from_seed, reset_test_rng};
use crate::mpfr::exact_bit_functions::{input, PRECISIONS};
use rug::{Float, Integer};
use zenith_float_num::Exponent;
use zenith_float_num::{Consts, ExactNum, RoundingMode, EXPONENT_MIN, WORD_BIT_SIZE};

/// Runs `f` with MPFR's exponent range at its widest (both ends: a lifted operand can pass the default `emax`).
fn widest_range<T>(f: impl FnOnce() -> T) -> T {
    use gmp_mpfr_sys::mpfr;
    let (emin, emax) = unsafe { (mpfr::get_emin(), mpfr::get_emax()) };
    unsafe {
        mpfr::set_emin(mpfr::get_emin_min());
        mpfr::set_emax(mpfr::get_emax_max());
    }
    let r = f();
    unsafe {
        mpfr::set_emin(emin);
        mpfr::set_emax(emax);
    }
    r
}

/// The significand and exponent of a finite `f`, `f = m * 2^e`, with the exponent at full `mpfr_exp_t` width
/// (`rug`'s `to_integer_exp` narrows it to `i32`, which panics with "overflow" once the exponents use the i64 range).
fn integer_exp(f: &Float) -> Result<(Integer, i64), Box<dyn std::error::Error>> {
    if !f.is_finite() {
        return Err("none value".into());
    }
    let mut i = Integer::new();
    let e = unsafe { gmp_mpfr_sys::mpfr::get_z_2exp(i.as_raw_mut(), f.as_raw()) };
    Ok((i, e))
}

/// `|x| mod |y|` with the sign of `x`, from the exact integer significands. Exact math at any exponent spread:
/// for `ex >= ey` the lifted dividend is reduced as `|mx| * (2^(ex - ey) mod |my|) mod |my|` (no `2^spread`
/// integer); for `ey - ex > 4096` the divisor exceeds the dividend and the remainder is `x` itself.
fn oracle(fx: &Float, fy: &Float) -> Result<Float, Box<dyn std::error::Error>> {
    let (mx, ex) = integer_exp(fx)?;
    let (my, ey) = integer_exp(fy)?;
    let (ax, ay) = (Integer::from(mx.abs_ref()), Integer::from(my.abs_ref()));
    let (mut r, e) = if ex >= ey {
        let k = Integer::from(ex - ey);
        let two = Integer::from(2);
        let t = two.pow_mod(&k, &ay).map_err(|_| "pow_mod")?;
        (Integer::from(&ax * &t) % &ay, ey)
    } else if ey - ex > 4096 {
        (ax, ex)
    } else {
        (ax % (ay << u32::try_from(ey - ex)?), ex)
    };
    if fx.is_sign_negative() {
        r = -r;
    }
    let bits = r.significant_bits().max(1);
    let sh = isize::try_from(e)?;
    Ok(widest_range(|| Float::with_val(bits, r) << sh))
}

/// `conv_to_mpfr` of `x * 2^lift` for values near the bottom of the range: `lift` moves the zenith subnormals
/// (below MPFR's `emin_min`) inside MPFR's range. Both operands and the result take the same `lift`, and
/// `rem(x * 2^s, y * 2^s) = rem(x, y) * 2^s` exactly. Takes `x` by value (the scaling consumes it); callers
/// rebuild the operand from its seed.
fn conv(
    w: usize,
    x: ExactNum,
    lift: i64,
    cc: &mut Consts,
) -> Result<Float, Box<dyn std::error::Error>> {
    let e = x.exponent().unwrap_or(0);
    // two steps: -EXPONENT_MIN does not fit in an Exponent
    let s = x
        .ldexp(-(e / 2), w, RoundingMode::None)
        .ldexp(-(e - e / 2), w, RoundingMode::None);
    let f = conv_to_mpfr(w, &s, cc)?;
    let sh = isize::try_from(i128::from(e) + i128::from(lift))?;
    Ok(widest_range(|| f << sh))
}

#[test]
fn rem_is_exact() -> Result<(), Box<dyn std::error::Error>> {
    reset_test_rng();
    let mut cc = Consts::new()?;
    // (x band, y band): y far below x, close to x, above x, and the subnormal band
    let bands: [((Exponent, Exponent), (Exponent, Exponent)); 5] = [
        ((-4, 4), (-4, 4)),
        ((200, 900), (-40, 40)),
        ((-40, 40), (60, 120)),
        (
            (EXPONENT_MIN + 1, EXPONENT_MIN + 300),
            (EXPONENT_MIN + 1, EXPONENT_MIN + 60),
        ),
        ((-100, 100), (EXPONENT_MIN + 1, EXPONENT_MIN + 200)),
    ];
    for px in PRECISIONS {
        for py in [24usize, 53, 64, 130, 255] {
            let wx = px.div_ceil(WORD_BIT_SIZE) * WORD_BIT_SIZE;
            let wy = py.div_ceil(WORD_BIT_SIZE) * WORD_BIT_SIZE;
            for (bx, by) in bands {
                if bx.0 < EXPONENT_MIN + 1000 && by.0 > 0 {
                    continue;
                }
                for _ in 0..4 {
                    // Every operand is drawn from its own seed and rebuilt from it for each consuming use.
                    let (sx, sy) = (draw_seed(), draw_seed());
                    let mk_x = || rebuild_from_seed(sx, || input(px, bx.0, bx.1, true));
                    let mk_y = || rebuild_from_seed(sy, || input(py, by.0, by.1, true));
                    let x = mk_x()?;
                    let y = mk_y()?;
                    // one lift for the pair: above the subnormals when either operand is near the bottom
                    let low = x.exponent().unwrap_or(0).min(y.exponent().unwrap_or(0));
                    let lift = if low < -(1i64 << 61) { 1i64 << 61 } else { 0 };
                    let got = x.rem(&y);
                    let info = format!("rem px={px} py={py} x={x:?} y={y:?}");
                    let fx = conv(wx, mk_x()?, lift, &mut cc)?;
                    let fy = conv(wy, mk_y()?, lift, &mut cc)?;
                    let want = oracle(&fx, &fy)?;
                    assert!(!got.inexact(), "{info}: rem reported inexact");
                    let w = got.mantissa_max_bit_len().unwrap_or(wx).max(1);
                    let fg = if got.is_zero() {
                        Float::with_val(w as u32, 0)
                    } else {
                        conv(w, mk_x()?.rem(&y), lift, &mut cc)?
                    };
                    assert_eq!(
                        fg,
                        want,
                        "{info}: zenith {}, exact {}",
                        crate::mpfr::common::show_mp(&fg),
                        crate::mpfr::common::show_mp(&want)
                    );
                    if want.is_zero() {
                        assert_eq!(got.is_negative(), x.is_negative(), "{info}: sign of zero");
                    }
                }
            }
        }
    }
    Ok(())
}
