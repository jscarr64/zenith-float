//! Results in the subnormal band round correctly in every mode, including `ToOdd`
//! (round half to odd) at exact ties on the subnormal grid, against exact references on the shifted frame (`shifted.rs`) (2.0.0 audit, round 2: item d).

use crate::mpfr::common::{draw_seed, rebuild_from_seed, reset_test_rng};
use crate::mpfr::exact_bit_functions::{input, MODES, PRECISIONS};
use crate::mpfr::shifted::{self as sh, Want};
use zenith_float_num::Exponent;
use zenith_float_num::{ExactNum, EXPONENT_MIN};

#[test]
fn subnormal_band_every_mode() -> Result<(), Box<dyn std::error::Error>> {
    reset_test_rng();
    let lo = EXPONENT_MIN;
    for p in PRECISIONS {
        for i in 0..24 {
            // Every operand is drawn from its own seed and rebuilt from it for each consuming use.
            let (sx, sy, sz) = (draw_seed(), draw_seed(), draw_seed());
            let mk_x = || rebuild_from_seed(sx, || input(p, lo + 1, lo + 2 * p as Exponent, true));
            let mk_y = || rebuild_from_seed(sy, || input(p, -(p as Exponent) - 2, -1, true));
            let mk_z = || rebuild_from_seed(sz, || input(p, lo + 1, lo + p as Exponent, true));
            let x = mk_x()?;
            let y = mk_y()?;
            let z = mk_z()?;
            // shifts that land exactly on and between subnormal grid points
            let k = i % 7 + 1;
            // the exact operands on the shifted frame (MPFR cannot hold the zenith subnormal grid)
            let (tx, ty, tz) = (sh::term_of(&x)?, sh::term_of(&y)?, sh::term_of(&z)?);
            let q = p as u32;
            for rm in MODES {
                let cases: [(&str, ExactNum, Want); 5] = [
                    (
                        "mul",
                        ExactNum::mul_ref(&x, &y, p, rm),
                        sh::ref_sum(&[tx.mul(&ty)], q, rm)?,
                    ),
                    (
                        "add",
                        ExactNum::add_ref(&x, &z, p, rm),
                        sh::ref_sum(&[tx.ldexp(0), tz.ldexp(0)], q, rm)?,
                    ),
                    (
                        "sub",
                        ExactNum::sub_ref(&x, &z, p, rm),
                        sh::ref_sum(&[tx.ldexp(0), tz.neg()], q, rm)?,
                    ),
                    (
                        "fma",
                        x.fma(&y, &z, p, rm),
                        sh::ref_sum(&[tx.mul(&ty), tz.ldexp(0)], q, rm)?,
                    ),
                    (
                        "ldexp",
                        mk_x()?.ldexp(-k, p, rm),
                        sh::ref_sum(&[tx.ldexp(-i128::from(k))], q, rm)?,
                    ),
                ];
                for (name, got, want) in cases.iter() {
                    let info = format!("{name} p={p} {rm:?} x={x:?}");
                    sh::check(got, want, !want.v.is_zero(), &info)?;
                }
            }
        }
    }
    Ok(())
}
