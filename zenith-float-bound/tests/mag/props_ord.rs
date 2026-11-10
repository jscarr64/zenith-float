//! CONTRACT 14.5, 2.1, 15.4; oracles: the integer order (tests/model/direct.rs) and the constructor rule; seeded splitmix64, replay with ZB_PROP_SEED/ZB_PROP_CASES
//! COVERS: Mag::exp, Mag::finite, Mag::mant, Ord for Mag, PartialEq for Mag, PartialOrd for Mag

use crate::model::direct;
use crate::model::draw::{draw_mag, draw_partner, draw_raw};
use crate::model::dyadic::Dy;
use crate::model::rng::run_cases;
use crate::support::{dcmp, from_magv};
use core::cmp::Ordering;
use zenith_bound::mag::Mag;

const CASES: u64 = 3000;

#[test]
fn prop_cmp_matches_the_integer_order_and_eq_follows_ord() -> Result<(), String> {
    run_cases(
        "prop_cmp_matches_the_integer_order_and_eq_follows_ord",
        CASES,
        &|rng| {
            let a = draw_mag(rng);
            let b = draw_partner(rng, &a);
            let (ma, mb) = (from_magv(&a)?, from_magv(&b)?);
            let got = ma.cmp(&mb);
            let want = direct::order(&a, &b);
            if got != want {
                return Err(format!(
                    "a={:?}, b={:?}: cmp()={:?} but direct::order()={:?}",
                    a, b, got, want
                ));
            }
            if (ma == mb) != (got == Ordering::Equal) {
                return Err(format!(
                    "a={:?}, b={:?}: (ma == mb)={} but got={:?}",
                    a,
                    b,
                    ma == mb,
                    got
                ));
            }
            if (ma < mb) != (got == Ordering::Less) {
                return Err(format!(
                    "a={:?}, b={:?}: (ma < mb)={} but got={:?}",
                    a,
                    b,
                    ma < mb,
                    got
                ));
            }
            if (ma > mb) != (got == Ordering::Greater) {
                return Err(format!(
                    "a={:?}, b={:?}: (ma > mb)={} but got={:?}",
                    a,
                    b,
                    ma > mb,
                    got
                ));
            }
            if ma.partial_cmp(&mb) != Some(got) {
                return Err(format!(
                    "a={:?}, b={:?}: partial_cmp()={:?} but got={:?}",
                    a,
                    b,
                    ma.partial_cmp(&mb),
                    got
                ));
            }
            Ok(())
        },
    )
}

#[test]
fn prop_cmp_is_antisymmetric_and_reflexive() -> Result<(), String> {
    run_cases("prop_cmp_is_antisymmetric_and_reflexive", CASES, &|rng| {
        let a = draw_mag(rng);
        let b = draw_partner(rng, &a);
        let ma = from_magv(&a)?;
        let mb = from_magv(&b)?;
        let cmp_ab = ma.cmp(&mb);
        let cmp_ba = mb.cmp(&ma);
        if cmp_ab != cmp_ba.reverse() {
            return Err(format!(
                "a={:?}, b={:?}: cmp(a,b)={:?} but cmp(b,a)={:?}",
                a, b, cmp_ab, cmp_ba
            ));
        }
        if ma.cmp(&ma) != Ordering::Equal {
            return Err(format!("a={:?}: cmp(a,a)={:?}", a, ma.cmp(&ma)));
        }
        Ok(())
    })
}

#[test]
fn prop_cmp_is_transitive() -> Result<(), String> {
    run_cases("prop_cmp_is_transitive", CASES, &|rng| {
        let a = draw_mag(rng);
        let b = draw_partner(rng, &a);
        let c = draw_partner(rng, &b);
        let ma = from_magv(&a)?;
        let mb = from_magv(&b)?;
        let mc = from_magv(&c)?;
        if (ma <= mb && mb <= mc) && !(ma <= mc) {
            return Err(format!("a={:?}, b={:?}, c={:?}: transitivity violated: ma <= mb && mb <= mc but not ma <= mc", a, b, c));
        }
        if (ma >= mb && mb >= mc) && !(ma >= mc) {
            return Err(format!("a={:?}, b={:?}, c={:?}: transitivity violated: ma >= mb && mb >= mc but not ma >= mc", a, b, c));
        }
        Ok(())
    })
}

#[test]
fn prop_constructor_matches_the_rule() -> Result<(), String> {
    run_cases("prop_constructor_matches_the_rule", CASES, &|rng| {
        let (mant, exp) = draw_raw(rng);
        let expected = if mant == 0 {
            if exp == 0 {
                Some((0i64, 0u64))
            } else {
                None
            }
        } else {
            let lz = mant.leading_zeros();
            let e = i128::from(exp) - i128::from(lz);
            if e < i128::from(i64::MIN) {
                None
            } else {
                Some((
                    i64::try_from(e)
                        .map_err(|_| format!("exp overflow for mant={:#x}, exp={}", mant, exp))?,
                    mant << lz,
                ))
            }
        };
        let got = Mag::finite(mant, exp);
        match (got, expected) {
            (Ok(m), Some((e, mnt))) => {
                if (m.exp(), m.mant()) != (Some(e), Some(mnt)) {
                    return Err(format!(
                        "mant={:#x}, exp={}: expected ({:?}, {:?}), got ({:?}, {:?})",
                        mant,
                        exp,
                        e,
                        mnt,
                        m.exp(),
                        m.mant()
                    ));
                }
            }
            (Err(_), None) => {}
            (Ok(m), None) => {
                return Err(format!(
                    "mant={:#x}, exp={}: expected None, got Ok(Mag {{ exp: {:?}, mant: {:?} }})",
                    mant,
                    exp,
                    m.exp(),
                    m.mant()
                ));
            }
            (Err(_), Some(_)) => {
                return Err(format!(
                    "mant={:#x}, exp={}: expected Some(...), got Err",
                    mant, exp
                ));
            }
        }
        Ok(())
    })
}

/// Independent value oracle: `a_mant * 2^a_exp == b_mant * 2^b_exp` (the common `-64` of clause 1.1 cancels), decided with u128 integers only.
fn same_value(a_mant: u64, a_exp: i64, b_mant: u64, b_exp: i64) -> bool {
    if a_mant == 0 || b_mant == 0 {
        return a_mant == 0 && b_mant == 0;
    }
    let (hi_mant, hi_exp, lo_mant, lo_exp) = if a_exp >= b_exp {
        (a_mant, a_exp, b_mant, b_exp)
    } else {
        (b_mant, b_exp, a_mant, a_exp)
    };
    let gap = i128::from(hi_exp) - i128::from(lo_exp);
    match u32::try_from(gap) {
        // the larger-exponent mantissa is >= 1, so a gap of 64 or more gives a value >= 2^64 > any u64 mantissa
        Ok(g) if g < 64 => (u128::from(hi_mant) << g) == u128::from(lo_mant),
        _ => false,
    }
}

/// Checks that a constructed `Mag` denotes exactly `mant * 2^(exp - 64)`, with the u128 oracle and with the dyadic model.
fn check_value(m: &Mag, mant: u64, exp: i64) -> Result<(), String> {
    let (got_exp, got_mant) = match (m.exp(), m.mant()) {
        (Some(e), Some(n)) => (e, n),
        other => {
            return Err(format!(
                "mant={:#x}, exp={}: constructed Mag has accessors {:?}, expected a finite value",
                mant, exp, other
            ))
        }
    };
    if !same_value(mant, exp, got_mant, got_exp) {
        return Err(format!(
            "mant={:#x}, exp={}: constructed Mag (exp {}, mant {:#x}) is not the same value (u128 oracle)",
            mant, exp, got_exp, got_mant
        ));
    }
    let want = Dy::from_mag(exp, mant);
    let got = Dy::from_mag(got_exp, got_mant);
    if dcmp(&want, &got)? != Ordering::Equal {
        return Err(format!(
            "mant={:#x}, exp={}: constructed Mag (exp {}, mant {:#x}) is not the same value (dyadic model)",
            mant, exp, got_exp, got_mant
        ));
    }
    Ok(())
}

#[test]
fn prop_constructor_is_value_preserving() -> Result<(), String> {
    run_cases("prop_constructor_is_value_preserving", CASES, &|rng| {
        let (mant, exp) = draw_raw(rng);
        if mant == 0 {
            return Ok(());
        }
        let lz = mant.leading_zeros();
        let modulus = match lz.checked_add(1) {
            Some(modulus) => modulus,
            None => return Ok(()),
        };
        let k = u32::from(rng.below(64)) % modulus;
        let shifted_mant = mant << k;
        let shifted_exp = match exp.checked_sub(i64::from(k)) {
            Some(e) => e,
            None => return Ok(()),
        };
        let mag1 = Mag::finite(mant, exp);
        let mag2 = Mag::finite(shifted_mant, shifted_exp);
        match (mag1, mag2) {
            (Ok(m1), Ok(m2)) => {
                check_value(&m1, mant, exp)?;
                check_value(&m2, shifted_mant, shifted_exp)?;
                if m1 != m2 || m1.exp() != m2.exp() || m1.mant() != m2.mant() {
                    return Err(format!(
                        "mant={:#x}, exp={} vs shifted_mant={:#x}, shifted_exp={}: values differ",
                        mant, exp, shifted_mant, shifted_exp
                    ));
                }
            }
            (Ok(_), Err(_)) | (Err(_), Ok(_)) => {
                return Err(format!("mant={:#x}, exp={} vs shifted_mant={:#x}, shifted_exp={}: one is Ok and other is Err", mant, exp, shifted_mant, shifted_exp));
            }
            _ => {}
        }
        Ok(())
    })
}
