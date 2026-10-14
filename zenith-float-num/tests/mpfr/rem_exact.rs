//! `rem` is exact: the truncated remainder `x - trunc(x/y)·y` is representable at the
//! precision of `x`, so it needs no rounding mode and is never inexact. Checked against an
//! exact big-integer remainder (2.0.0 audit, round 2 leftovers).

use crate::mpfr::common::{conv_to_mpfr, reset_test_rng, test_random};
use crate::mpfr::exact_bit_functions::{input, PRECISIONS};
use rug::{Float, Integer};
use zenith_float_num::{Consts, ExactNum, RoundingMode, EXPONENT_MIN, WORD_BIT_SIZE};

/// `|x| mod |y|` with the sign of `x`, from the exact integer significands.
fn oracle(fx: &Float, fy: &Float) -> Float {
    let (mx, ex) = fx.to_integer_exp().unwrap();
    let (my, ey) = fy.to_integer_exp().unwrap();
    let e = ex.min(ey);
    let xi = Integer::from(mx.abs_ref()) << (ex - e) as u32;
    let yi = Integer::from(my.abs_ref()) << (ey - e) as u32;
    let mut r = xi % yi;
    if fx.is_sign_negative() {
        r = -r;
    }
    let bits = r.significant_bits().max(1);
    let old = unsafe { gmp_mpfr_sys::mpfr::get_emin() };
    unsafe { gmp_mpfr_sys::mpfr::set_emin(gmp_mpfr_sys::mpfr::get_emin_min()) };
    let f = Float::with_val(bits, r) << e;
    unsafe { gmp_mpfr_sys::mpfr::set_emin(old) };
    f
}

/// `conv_to_mpfr` for values near the bottom of the range.
fn conv(w: usize, x: &ExactNum, cc: &mut Consts) -> Float {
    let e = x.exponent().unwrap_or(0);
    let s = x.ldexp(-e, w, RoundingMode::None);
    let f = conv_to_mpfr(w, &s, cc);
    let old = unsafe { gmp_mpfr_sys::mpfr::get_emin() };
    unsafe { gmp_mpfr_sys::mpfr::set_emin(gmp_mpfr_sys::mpfr::get_emin_min()) };
    let f = f << e;
    unsafe { gmp_mpfr_sys::mpfr::set_emin(old) };
    f
}

#[test]
fn rem_is_exact() {
    reset_test_rng();
    let mut cc = Consts::new().unwrap();
    // (x band, y band): y far below x, close to x, above x, and the subnormal band
    let bands: [((i32, i32), (i32, i32)); 5] = [
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
                    let x = input(px, bx.0, bx.1, true);
                    let y = input(py, by.0, by.1, true);
                    // the integer oracle needs the two exponents within reach of each other
                    let spread =
                        x.exponent().unwrap_or(0) as i64 - y.exponent().unwrap_or(0) as i64;
                    if spread.abs() > 4096 && test_random::<bool>() {
                        continue;
                    }
                    let got = x.rem(&y);
                    let info = format!("rem px={px} py={py} x={x:?} y={y:?}");
                    let fx = conv(wx, &x, &mut cc);
                    let fy = conv(wy, &y, &mut cc);
                    let want = oracle(&fx, &fy);
                    assert!(!got.inexact(), "{info}: rem reported inexact");
                    let w = got.mantissa_max_bit_len().unwrap_or(wx).max(1);
                    let fg = if got.is_zero() {
                        Float::with_val(w as u32, 0)
                    } else {
                        conv(w, &got, &mut cc)
                    };
                    assert_eq!(fg, want, "{info}: zenith {fg}, exact {want}");
                    if want.is_zero() {
                        assert_eq!(got.is_negative(), x.is_negative(), "{info}: sign of zero");
                    }
                }
            }
        }
    }
}
