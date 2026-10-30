//! CONTRACT 14.7, 3.2, 3.3; oracle: i128 placement of the power of two (exponent field
//! `pre_round_exp - prec + 1` for `ulp_mag`, `mid_exp - prec` for `half_ulp`), and the
//! dyadic model for in-range values. `prec` is the full `Precision` (usize).

use crate::model::dyadic::{Dy, MagV};
use crate::model::oracle::TWO63;
use crate::support::{dcmp, dy_of, to_magv};
use core::cmp::Ordering;
use zenith_bound::mag::Mag;
use zenith_bound::MAX_PRECISION;

const PRECS: [usize; 10] = [
    0,
    1,
    2,
    31,
    32,
    33,
    1 << 31,
    0xFFFF_FFFF,
    MAX_PRECISION,
    usize::MAX,
];
const EXPS: [i64; 7] = [i64::MIN, i64::MIN + 1, -1, 0, 1, i64::MAX - 1, i64::MAX];

fn placed(mag_exp: i128) -> Result<MagV, String> {
    if mag_exp > i128::from(i64::MAX) {
        Ok(MagV::Infinity)
    } else if mag_exp < i128::from(i64::MIN) {
        Ok(MagV::Finite {
            exp: i64::MIN,
            mant: TWO63,
        })
    } else {
        let exp =
            i64::try_from(mag_exp).map_err(|_| format!("exponent {mag_exp} is not an i64"))?;
        Ok(MagV::Finite { exp, mant: TWO63 })
    }
}

fn ulp_field(pre_round_exp: i64, prec: usize) -> Result<i128, String> {
    i128::from(pre_round_exp)
        .checked_sub(i128::try_from(prec).map_err(|e| format!("prec {}: {}", prec, e))?)
        .and_then(|diff| diff.checked_add(1))
        .ok_or_else(|| format!("ulp field overflow at {pre_round_exp}, {prec}"))
}

fn half_field(mid_exp: i64, prec: usize) -> Result<i128, String> {
    i128::from(mid_exp)
        .checked_sub(i128::try_from(prec).map_err(|e| format!("prec {}: {}", prec, e))?)
        .ok_or_else(|| format!("half-ulp field overflow at {mid_exp}, {prec}"))
}

fn matches_power(got: &MagV, mag_exp: i128) -> Result<(), String> {
    if mag_exp < i128::from(i64::MIN) || mag_exp > i128::from(i64::MAX) {
        return Ok(());
    }
    let power = mag_exp
        .checked_sub(1)
        .ok_or_else(|| format!("power exponent overflow at {mag_exp}"))?;
    if dcmp(&dy_of(got)?, &Dy::pow2(power))? != Ordering::Equal {
        return Err(format!("{got:?} is not 2^{power}"));
    }
    Ok(())
}

#[test]
fn ulp_mag_uses_the_full_prec_and_overflows_only_at_max_exp_prec_zero() -> Result<(), String> {
    for &pre_round_exp in EXPS.iter() {
        for &prec in PRECS.iter() {
            let mag_exp = ulp_field(pre_round_exp, prec)?;
            let want = placed(mag_exp)?;
            let got = to_magv(&Mag::ulp_mag(pre_round_exp, prec))?;
            if got != want {
                return Err(format!(
                    "ulp_mag({pre_round_exp}, {prec}) expected {want:?}, got {got:?}"
                ));
            }
            if pre_round_exp == i64::MAX && prec == 0 && !matches!(got, MagV::Infinity) {
                return Err(String::from("ulp_mag(i64::MAX, 0) must be Infinity"));
            }
            if pre_round_exp == i64::MAX && prec >= 1 && matches!(got, MagV::Infinity) {
                return Err(format!("ulp_mag(i64::MAX, {prec}) must be finite"));
            }
            matches_power(&got, mag_exp)?;
        }
    }
    Ok(())
}

#[test]
fn half_ulp_prec_zero_is_two_to_the_mid_exp_minus_one() -> Result<(), String> {
    for &mid_exp in EXPS.iter() {
        for &prec in PRECS.iter() {
            let mag_exp = half_field(mid_exp, prec)?;
            let want = placed(mag_exp)?;
            let got = to_magv(&Mag::half_ulp(mid_exp, prec))?;
            if got != want {
                return Err(format!(
                    "half_ulp({mid_exp}, {prec}) expected {want:?}, got {got:?}"
                ));
            }
            if prec == 0 {
                let power = i128::from(mid_exp)
                    .checked_sub(1)
                    .ok_or_else(|| format!("half-ulp power overflow at {mid_exp}"))?;
                if dcmp(&dy_of(&got)?, &Dy::pow2(power))? != Ordering::Equal {
                    return Err(format!(
                        "half_ulp({mid_exp}, 0) must be 2^({mid_exp}-1), got {got:?}"
                    ));
                }
            }
            matches_power(&got, mag_exp)?;
        }
    }
    Ok(())
}

#[test]
fn smaller_than_half_ulp_uses_the_full_prec_width() -> Result<(), String> {
    let mants = [TWO63, TWO63 + 1, u64::MAX];
    let precs = [
        0usize,
        1,
        31,
        32,
        1 << 31,
        0xFFFF_FFFF,
        MAX_PRECISION,
        usize::MAX,
    ];
    let mids = [i64::MIN, 0, 10, i64::MAX];
    let exps = [i64::MIN, i64::MIN + 1, 0, 5, i64::MAX];
    for &mid_exp in mids.iter() {
        for &prec in precs.iter() {
            let threshold = match i128::from(mid_exp)
                .checked_sub(i128::try_from(prec).map_err(|e| format!("prec {}: {}", prec, e))?)
                .and_then(|diff| diff.checked_sub(1))
            {
                Some(threshold) => threshold,
                None => {
                    return Err(format!("threshold overflow at {mid_exp}, {prec}"));
                }
            };
            for &exp in exps.iter() {
                for &mant in mants.iter() {
                    let sample = MagV::Finite { exp, mant };
                    let mag = Mag::finite(mant, exp)
                        .map_err(|err| format!("Mag::finite({mant}, {exp}) failed: {err:?}"))?;
                    let want = dcmp(&dy_of(&sample)?, &Dy::pow2(threshold))? == Ordering::Less;
                    let got = mag.smaller_than_half_ulp(mid_exp, prec);
                    if got != want {
                        return Err(format!(
                            "smaller_than_half_ulp(Finite {{ exp: {exp}, mant: {mant} }}, {mid_exp}, {prec}) expected {want}, got {got}"
                        ));
                    }
                }
            }
            if !Mag::zero().smaller_than_half_ulp(mid_exp, prec) {
                return Err(format!(
                    "zero must be smaller than half an ulp at mid_exp {mid_exp}, prec {prec}"
                ));
            }
            if Mag::infinity().smaller_than_half_ulp(mid_exp, prec) {
                return Err(format!(
                    "Infinity must not be smaller than half an ulp at mid_exp {mid_exp}, prec {prec}"
                ));
            }
        }
    }
    Ok(())
}
