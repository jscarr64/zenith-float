//! CONTRACT 14.7, 1.4, 3.2; oracle: dyadic-rational model (tests/model)
//! COVERS: Mag::half_ulp, Mag::infinity, Mag::smaller_than_half_ulp, Mag::ulp_mag, Mag::zero

use crate::model::dyadic::{Dy, MagV};
use crate::model::oracle::{min_pos_mag, TWO63};
use crate::model::sets::mantissas;
use crate::support::{dcmp, dy_of, from_magv, to_magv};
use core::cmp::Ordering;
use zenith_bound::mag::Mag;

/// A `Precision` as `i64` / `i128` for the oracle formulas (every prec in this file fits both).
fn p64(p: usize) -> i64 {
    i64::try_from(p).unwrap_or(i64::MAX)
}

fn p128(p: usize) -> i128 {
    i128::try_from(p).unwrap_or(i128::MAX)
}

const PRECS: [usize; 10] = [1, 2, 3, 4, 24, 53, 64, 65, 113, 127];
const EXPS: [i64; 11] = [-1000, -65, -64, -1, 0, 1, 2, 63, 64, 65, 1000];
/// The three mantissas CONTRACT 1.4 names for its boundary set.
const CONTRACT_MANTS: [u64; 3] = [TWO63, TWO63 + 1, u64::MAX];

#[test]
fn ulp_mag_is_two_to_the_pre_round_exp_minus_prec() -> Result<(), String> {
    for e in EXPS.iter() {
        for p in PRECS.iter() {
            let k = *e - p64(*p);
            let r = to_magv(&Mag::ulp_mag(*e, *p))?;
            if dcmp(&dy_of(&r)?, &Dy::pow2(i128::from(k)))? != Ordering::Equal {
                return Err(format!(
                    "ulp_mag({},{}) value differs from 2^{}, got {:?}",
                    e, p, k, r
                ));
            }
            if r != (MagV::Finite {
                exp: k + 1,
                mant: TWO63,
            }) {
                return Err(format!(
                    "ulp_mag({},{}) expected finite with exp={}, mant={}, got {:?}",
                    e,
                    p,
                    k + 1,
                    TWO63,
                    r
                ));
            }
        }
    }
    Ok(())
}

#[test]
fn ulp_mag_exponent_below_i64_min_rounds_up_to_smallest_positive_mag() -> Result<(), String> {
    let cases = [
        (i64::MIN, 2usize, min_pos_mag()),
        (i64::MIN, 3usize, min_pos_mag()),
        (i64::MIN, 64usize, min_pos_mag()),
        (i64::MIN, 127usize, min_pos_mag()),
        (i64::MIN + 1, 3usize, min_pos_mag()),
        (i64::MIN + 1, 64usize, min_pos_mag()),
        (i64::MIN + 1, 127usize, min_pos_mag()),
    ];
    for (e, p, want) in cases.iter() {
        let r = to_magv(&Mag::ulp_mag(*e, *p))?;
        if r != *want {
            return Err(format!(
                "ulp_mag({},{}) expected {:?}, got {:?}",
                e, p, want, r
            ));
        }
    }

    for p in [1usize, 2, 64, 127] {
        let e = i64::MIN + p64(p) - 1;
        let r = to_magv(&Mag::ulp_mag(e, p))?;
        if r != min_pos_mag() {
            return Err(format!(
                "ulp_mag({},{}) expected {:?}, got {:?}",
                e,
                p,
                min_pos_mag(),
                r
            ));
        }
    }

    for p in [1usize, 2, 64, 127] {
        let e = i64::MIN + p64(p);
        let r = to_magv(&Mag::ulp_mag(e, p))?;
        let want = MagV::Finite {
            exp: i64::MIN + 1,
            mant: TWO63,
        };
        if r != want {
            return Err(format!(
                "ulp_mag({},{}) expected {:?}, got {:?}",
                e, p, want, r
            ));
        }
    }

    Ok(())
}

#[test]
fn half_ulp_is_two_to_the_mid_exp_minus_prec_minus_one() -> Result<(), String> {
    for e in EXPS.iter() {
        for p in PRECS.iter() {
            let k = *e - p64(*p) - 1;
            let r = to_magv(&Mag::half_ulp(*e, *p))?;
            if dcmp(&dy_of(&r)?, &Dy::pow2(i128::from(k)))? != Ordering::Equal {
                return Err(format!(
                    "half_ulp({},{}) value differs from 2^{}, got {:?}",
                    e, p, k, r
                ));
            }
            if r != (MagV::Finite {
                exp: k + 1,
                mant: TWO63,
            }) {
                return Err(format!(
                    "half_ulp({},{}) expected finite with exp={}, mant={}, got {:?}",
                    e,
                    p,
                    k + 1,
                    TWO63,
                    r
                ));
            }
        }
    }
    Ok(())
}

#[test]
fn smaller_than_half_ulp_matches_model_on_boundary_grid() -> Result<(), String> {
    for &mid_exp in [-70i64, 0, 10, 130].iter() {
        for p in PRECS.iter() {
            let h = mid_exp - p64(*p) - 1;
            for exp in (h - 2)..=(h + 2) {
                for &mant in mantissas().iter() {
                    let x = from_magv(&MagV::Finite { exp, mant })?;
                    let want = dcmp(
                        &dy_of(&MagV::Finite { exp, mant })?,
                        &Dy::pow2(i128::from(h)),
                    )? == Ordering::Less;
                    let got = x.smaller_than_half_ulp(mid_exp, *p);
                    if got != want {
                        return Err(format!("smaller_than_half_ulp(Finite {{ exp: {}, mant: {} }}, {}, {}) expected {}, got {}", exp, mant, mid_exp, p, want, got));
                    }
                }
            }
        }
    }
    Ok(())
}

#[test]
fn smaller_than_half_ulp_contract_1_4_boundary_set() -> Result<(), String> {
    let mid_exp = 10;
    let prec = 4;
    for exp in 3..=8 {
        for &mant in CONTRACT_MANTS.iter() {
            let x = from_magv(&MagV::Finite { exp, mant })?;
            let want = exp <= 5;
            let got = x.smaller_than_half_ulp(mid_exp, prec);
            if got != want {
                return Err(format!("smaller_than_half_ulp(Finite {{ exp: {}, mant: {} }}, {}, {}) expected {}, got {}", exp, mant, mid_exp, prec, want, got));
            }
        }
    }
    Ok(())
}

#[test]
fn half_ulp_underflow_near_i64_min_rounds_up_to_smallest_positive_mag() -> Result<(), String> {
    for p in PRECS.iter() {
        let p64 = p64(*p);
        for e in [i64::MIN, i64::MIN + p64 - 1, i64::MIN + p64] {
            let r = to_magv(&Mag::half_ulp(e, *p))?;
            if r != min_pos_mag() {
                return Err(format!(
                    "half_ulp({},{}) expected {:?}, got {:?}",
                    e,
                    p,
                    min_pos_mag(),
                    r
                ));
            }
        }
        let e = i64::MIN + p64 + 1;
        let r = to_magv(&Mag::half_ulp(e, *p))?;
        let want = MagV::Finite {
            exp: i64::MIN + 1,
            mant: TWO63,
        };
        if r != want {
            return Err(format!(
                "half_ulp({},{}) expected {:?}, got {:?}",
                e, p, want, r
            ));
        }
    }
    Ok(())
}

#[test]
fn smaller_than_half_ulp_near_i64_min_matches_model_and_is_false_below_the_range(
) -> Result<(), String> {
    for p in PRECS.iter() {
        let p64 = p64(*p);
        for mid_exp in [
            i64::MIN,
            i64::MIN + 1,
            i64::MIN + p64 - 1,
            i64::MIN + p64,
            i64::MIN + p64 + 1,
            i64::MIN + p64 + 2,
        ] {
            let h = i128::from(mid_exp) - p128(*p) - 1;
            for exp in [
                i64::MIN,
                i64::MIN + 1,
                i64::MIN + 2,
                i64::MIN + 3,
                i64::MIN + 64,
                i64::MIN + 130,
            ] {
                for &mant in mantissas().iter() {
                    let x = from_magv(&MagV::Finite { exp, mant })?;
                    let want =
                        dcmp(&dy_of(&MagV::Finite { exp, mant })?, &Dy::pow2(h))? == Ordering::Less;
                    let got = x.smaller_than_half_ulp(mid_exp, *p);
                    if got != want {
                        return Err(format!("smaller_than_half_ulp(Finite {{ exp: {}, mant: {} }}, {}, {}) expected {}, got {}", exp, mant, mid_exp, p, want, got));
                    }
                    if h < i128::from(i64::MIN) && got {
                        return Err(format!("smaller_than_half_ulp(Finite {{ exp: {}, mant: {} }}, {}, {}) must be false: the threshold 2^{} is below every nonzero Mag", exp, mant, mid_exp, p, h));
                    }
                }
            }
        }
    }
    Ok(())
}

#[test]
fn smaller_than_half_ulp_zero_is_true_and_infinity_is_false() -> Result<(), String> {
    for p in PRECS.iter() {
        for mid_exp in [i64::MIN, -1000, 0, 1000, i64::MAX] {
            if !Mag::zero().smaller_than_half_ulp(mid_exp, *p) {
                return Err(format!(
                    "zero must be smaller than half an ulp: mid_exp {}, prec {}",
                    mid_exp, p
                ));
            }
            if Mag::infinity().smaller_than_half_ulp(mid_exp, *p) {
                return Err(format!(
                    "Infinity must not be smaller than half an ulp: mid_exp {}, prec {}",
                    mid_exp, p
                ));
            }
        }
    }
    Ok(())
}

/// CONTRACT 14.7 with prec = 0 (outside 8.1 but total): ulp_mag(e, 0) is 2^e, i.e. Finite { exp: e + 1, mant: 2^63 }, and a result above the i64 exponent range is Infinity (3.3).
#[test]
fn ulp_mag_prec_zero_is_two_to_the_pre_round_exp_and_overflows_at_i64_max() -> Result<(), String> {
    let finite_cases: [(i64, i64); 7] = [
        (i64::MIN, i64::MIN + 1),
        (i64::MIN + 1, i64::MIN + 2),
        (-1, 0),
        (0, 1),
        (1, 2),
        (i64::MAX - 2, i64::MAX - 1),
        (i64::MAX - 1, i64::MAX),
    ];
    for (e, want_exp) in finite_cases.iter() {
        let r = to_magv(&Mag::ulp_mag(*e, 0))?;
        let want = MagV::Finite {
            exp: *want_exp,
            mant: TWO63,
        };
        if r != want {
            return Err(format!("ulp_mag({},0) expected {:?}, got {:?}", e, want, r));
        }
    }
    let r = to_magv(&Mag::ulp_mag(i64::MAX, 0))?;
    if r != MagV::Infinity {
        return Err(format!(
            "ulp_mag(i64::MAX,0) expected Infinity, got {:?}",
            r
        ));
    }
    Ok(())
}

/// CONTRACT 14.7 with prec = 0: half_ulp(m, 0) is 2^(m - 1), i.e. Finite { exp: m, mant: 2^63 }; at m = i64::MIN that is exactly the smallest positive Mag (no underflow), and no exponent overflows.
#[test]
fn half_ulp_prec_zero_is_two_to_the_mid_exp_minus_one() -> Result<(), String> {
    let mids: [i64; 7] = [i64::MIN, i64::MIN + 1, -1, 0, 1, i64::MAX - 1, i64::MAX];
    for m in mids.iter() {
        let r = to_magv(&Mag::half_ulp(*m, 0))?;
        let want = MagV::Finite {
            exp: *m,
            mant: TWO63,
        };
        if r != want {
            return Err(format!(
                "half_ulp({},0) expected {:?}, got {:?}",
                m, want, r
            ));
        }
    }
    if to_magv(&Mag::half_ulp(i64::MIN, 0))? != min_pos_mag() {
        return Err(String::from(
            "half_ulp(i64::MIN,0) must be the smallest positive Mag",
        ));
    }
    Ok(())
}

/// CONTRACT 14.7: the Mag value 2^k is `Finite { exp: k + 1, mant: 2^63 }` (clause 1.2); an exponent below `i64::MIN` rounds up to the smallest positive Mag (3.2), one above `i64::MAX` is `Infinity` (3.3). Evaluated in `i128`.
fn pow2_mag(k: i128) -> Result<MagV, String> {
    let exp = k
        .checked_add(1)
        .ok_or_else(|| String::from("exponent overflow in the test oracle"))?;
    if exp < i128::from(i64::MIN) {
        return Ok(min_pos_mag());
    }
    match i64::try_from(exp) {
        Ok(e) => Ok(MagV::Finite {
            exp: e,
            mant: TWO63,
        }),
        Err(_) => Ok(MagV::Infinity),
    }
}

/// Precisions around the 31-bit and 32-bit boundaries of `u32`.
const BIG_PRECS: [usize; 4] = [
    (1usize << 31) - 1,
    1usize << 31,
    (1usize << 31) + 1,
    0xFFFF_FFFF,
];

/// CONTRACT 14.7 with the pre-round exponent near `i64::MAX` and prec >= 1: `ulp_mag(e, prec)` is `Finite { exp: e - prec + 1, mant: 2^63 }`, which is still finite (only `e - prec + 1 > i64::MAX` is `Infinity`, clause 3.3).
#[test]
fn ulp_mag_near_i64_max_with_prec_at_least_one_is_finite() -> Result<(), String> {
    let cases: [(i64, usize, i64); 6] = [
        (i64::MAX, 1, i64::MAX),
        (i64::MAX, 2, i64::MAX - 1),
        (i64::MAX, 64, i64::MAX - 63),
        (i64::MAX, 127, i64::MAX - 126),
        (i64::MAX - 1, 1, i64::MAX - 1),
        (i64::MAX - 1, 2, i64::MAX - 2),
    ];
    for (e, p, want_exp) in cases.iter() {
        let r = to_magv(&Mag::ulp_mag(*e, *p))?;
        let want = MagV::Finite {
            exp: *want_exp,
            mant: TWO63,
        };
        if r != want {
            return Err(format!(
                "ulp_mag({},{}) expected {:?}, got {:?}",
                e, p, want, r
            ));
        }
    }
    Ok(())
}

/// CONTRACT 14.7 with the mid exponent near `i64::MAX` and prec >= 1: `half_ulp(m, prec)` is `Finite { exp: m - prec, mant: 2^63 }`.
#[test]
fn half_ulp_near_i64_max_with_prec_at_least_one_is_finite() -> Result<(), String> {
    let cases: [(i64, usize, i64); 6] = [
        (i64::MAX, 1, i64::MAX - 1),
        (i64::MAX, 2, i64::MAX - 2),
        (i64::MAX, 64, i64::MAX - 64),
        (i64::MAX, 127, i64::MAX - 127),
        (i64::MAX - 1, 1, i64::MAX - 2),
        (i64::MAX - 1, 2, i64::MAX - 3),
    ];
    for (m, p, want_exp) in cases.iter() {
        let r = to_magv(&Mag::half_ulp(*m, *p))?;
        let want = MagV::Finite {
            exp: *want_exp,
            mant: TWO63,
        };
        if r != want {
            return Err(format!(
                "half_ulp({},{}) expected {:?}, got {:?}",
                m, p, want, r
            ));
        }
    }
    Ok(())
}

/// CONTRACT 14.7 for every `u32` prec: a prec of 2^31 or more must not be truncated to 31 bits. Expected values come from the formulas in `i128`.
#[test]
fn ulp_mag_and_half_ulp_are_exact_for_prec_at_or_above_2_pow_31() -> Result<(), String> {
    let exps: [i64; 5] = [i64::MIN, -1000, 0, 1000, i64::MAX];
    for p in BIG_PRECS.iter() {
        for e in exps.iter() {
            let ulp = to_magv(&Mag::ulp_mag(*e, *p))?;
            let want_ulp = pow2_mag(i128::from(*e) - p128(*p))?;
            if ulp != want_ulp {
                return Err(format!(
                    "ulp_mag({},{}) expected {:?}, got {:?}",
                    e, p, want_ulp, ulp
                ));
            }
            let half = to_magv(&Mag::half_ulp(*e, *p))?;
            let want_half = pow2_mag(i128::from(*e) - p128(*p) - 1)?;
            if half != want_half {
                return Err(format!(
                    "half_ulp({},{}) expected {:?}, got {:?}",
                    e, p, want_half, half
                ));
            }
        }
    }
    Ok(())
}

/// CONTRACT 1.4 and 14.7 for every `u32` prec: a finite nonzero Mag is smaller than half an ulp iff `exp <= mid_exp - prec - 1`, evaluated in `i128`; zero is smaller and `Infinity` is not.
#[test]
fn smaller_than_half_ulp_is_exact_for_prec_at_or_above_2_pow_31() -> Result<(), String> {
    for p in BIG_PRECS.iter() {
        for mid_exp in [-1000i64, 0, 1000] {
            let h = i128::from(mid_exp) - p128(*p) - 1;
            for d in [-1i128, 0, 1] {
                let exp = i64::try_from(h + d)
                    .map_err(|_| format!("exponent {} out of range in the test oracle", h + d))?;
                for mant in [TWO63, u64::MAX] {
                    let x = from_magv(&MagV::Finite { exp, mant })?;
                    let want = d <= 0;
                    let got = x.smaller_than_half_ulp(mid_exp, *p);
                    if got != want {
                        return Err(format!("smaller_than_half_ulp(Finite {{ exp: {}, mant: {} }}, {}, {}) expected {}, got {}", exp, mant, mid_exp, p, want, got));
                    }
                }
            }
            if !Mag::zero().smaller_than_half_ulp(mid_exp, *p) {
                return Err(format!(
                    "zero must be smaller than half an ulp: mid_exp {}, prec {}",
                    mid_exp, p
                ));
            }
            if Mag::infinity().smaller_than_half_ulp(mid_exp, *p) {
                return Err(format!(
                    "Infinity must not be smaller than half an ulp: mid_exp {}, prec {}",
                    mid_exp, p
                ));
            }
        }
    }
    Ok(())
}
