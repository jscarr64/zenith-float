//! Harper S2 gate: wide-range edges of the Mag/ExactNum bridge (CONTRACT 16.S2.1 to 16.S2.4). Hand-worked rows, each
//! value in the comment beside it: signed zero, i64 exponent extremes, a 100000-bit precision with the last bit deciding,
//! cancellation in both cmp_exact_diff paths, NaN/Inf for every argument, round trips both ways, directed-rounding
//! boundaries and every to_exact_num error detail, and the Constant::ALL order.
//! COVERS: Mag::cmp_exact_num, Mag::cmp_exact_diff, Mag::from_exact_num, Mag::from_exact_num_lower, Mag::to_exact_num, Constant::Pi, Constant::E, Constant::Ln2, Constant::Ln10, Constant::Sqrt2, Constant::Phi, Constant::EulerGamma

use core::cmp::Ordering;

use zenith_bound::bridge::Constant;
use zenith_bound::error::BoundError;
use zenith_bound::mag::Mag;
use zenith_float_num::{ExactNum, RawKind, RawParts, RoundingMode, EXPONENT_MAX, EXPONENT_MIN};

const TOP: u64 = 0x8000_0000_0000_0000;
/// The huge precision: values with a set bit 100000 binades below their top bit.
const DEEP: i64 = 100_000;

fn mag(mant: u64, exp: i64) -> Result<Mag, String> {
    Mag::finite(mant, exp).map_err(|e| format!("Mag::finite({:#x}, {}): {:?}", mant, exp, e))
}

fn num(v: i64) -> ExactNum {
    ExactNum::from(v)
}

/// `0.mant * 2^exp` exactly at 64 bits.
fn parts(mant: u64, exp: i64) -> Result<ExactNum, String> {
    ExactNum::from_mag_parts(mant, exp, 64)
        .map_err(|e| format!("from_mag_parts({:#x}, {}): {:?}", mant, exp, e))
}

/// `2^k` exactly.
fn pow2(k: i64) -> Result<ExactNum, String> {
    parts(TOP, k.checked_add(1).ok_or("pow2 k")?)
}

fn neg(x: &ExactNum) -> ExactNum {
    num(0).sub_full_prec(x)
}

fn add(a: i64, b: i64) -> Result<i64, String> {
    a.checked_add(b).ok_or_else(|| format!("{} + {}", a, b))
}

fn sub(a: i64, b: i64) -> Result<i64, String> {
    a.checked_sub(b).ok_or_else(|| format!("{} - {}", a, b))
}

fn want(what: &str, got: Option<Ordering>, w: Ordering) -> Result<(), String> {
    if got == Some(w) {
        Ok(())
    } else {
        Err(format!("{}: got {:?}, want {:?}", what, got, w))
    }
}

/// `from_exact_num(x)` is exactly `up` and `from_exact_num_lower(x)` exactly `lo`.
fn want_up_lo(what: &str, x: &ExactNum, up: &Mag, lo: &Mag) -> Result<(), String> {
    let (gu, gl) = (Mag::from_exact_num(x), Mag::from_exact_num_lower(x));
    if gu == *up && gl == *lo {
        Ok(())
    } else {
        Err(format!(
            "{}: got up {:?} lo {:?}, want {:?} {:?}",
            what, gu, gl, up, lo
        ))
    }
}

/// Negative zero, built from raw parts (zf arithmetic need not produce it).
fn neg_zero() -> Result<ExactNum, String> {
    let words = [0u64];
    ExactNum::from_raw_parts(&RawParts {
        kind: RawKind::Zero,
        words: &words,
        precision: 64,
        negative: true,
        exponent: 0,
        inexact: false,
        op_inexact: false,
    })
    .map_err(|e| format!("-0: {:?}", e))
}

#[test]
fn s2_edges_signed_zero() -> Result<(), String> {
    let nz = neg_zero()?;
    if !nz.is_zero() || !nz.as_raw_parts().negative {
        return Err(format!("-0 not built: {:?}", nz.as_raw_parts()));
    }
    let two = mag(TOP, 2)?; // 2
    want("0 vs |-0|", Mag::zero().cmp_exact_num(&nz), Ordering::Equal)?;
    want(
        "1 vs |-0|",
        mag(TOP, 1)?.cmp_exact_num(&nz),
        Ordering::Greater,
    )?;
    want(
        "inf vs |-0|",
        Mag::infinity().cmp_exact_num(&nz),
        Ordering::Greater,
    )?;
    want(
        "0 vs |-0 - 0|",
        Mag::zero().cmp_exact_diff(&nz, &num(0)),
        Ordering::Equal,
    )?;
    want(
        "0 vs |0 - -0|",
        Mag::zero().cmp_exact_diff(&num(0), &nz),
        Ordering::Equal,
    )?;
    want(
        "0 vs |-0 - -0|",
        Mag::zero().cmp_exact_diff(&nz, &nz),
        Ordering::Equal,
    )?;
    want(
        "2 vs |-0 - -2|",
        two.cmp_exact_diff(&nz, &num(-2)),
        Ordering::Equal,
    )?;
    want(
        "2 vs |2 - -0|",
        two.cmp_exact_diff(&num(2), &nz),
        Ordering::Equal,
    )?;
    want(
        "2^-60 vs |-0 - 2^-61|",
        mag(TOP, -59)?.cmp_exact_diff(&nz, &pow2(-61)?),
        Ordering::Greater,
    )?;
    want_up_lo("-0", &nz, &Mag::zero(), &Mag::zero())?;
    want_up_lo("+0", &num(0), &Mag::zero(), &Mag::zero())?;
    let z = Mag::zero()
        .to_exact_num()
        .map_err(|e| format!("to_exact_num(0): {:?}", e))?;
    if !z.is_zero() || z.as_raw_parts().negative {
        return Err(format!("to_exact_num(0) = {:?}", z.as_raw_parts()));
    }
    Ok(())
}

#[test]
fn s2_edges_i64_exponent_extremes() -> Result<(), String> {
    let (hi, lo) = (EXPONENT_MAX, EXPONENT_MIN);
    let big = parts(u64::MAX, hi)?; // 2^EXPONENT_MAX (1 - 2^-64), the largest 64-bit zf value
    let tiny = parts(TOP, lo)?; // 2^(EXPONENT_MIN - 1)
    let top = mag(u64::MAX, i64::MAX)?; // the largest finite Mag
    let bottom = mag(TOP, i64::MIN)?; // the smallest nonzero Mag
    want("max vs big", top.cmp_exact_num(&big), Ordering::Greater)?;
    want(
        "max vs -tiny",
        top.cmp_exact_num(&neg(&tiny)),
        Ordering::Greater,
    )?;
    want("min vs tiny", bottom.cmp_exact_num(&tiny), Ordering::Less)?;
    want(
        "min vs -big",
        bottom.cmp_exact_num(&neg(&big)),
        Ordering::Less,
    )?;
    want("min vs 0", bottom.cmp_exact_num(&num(0)), Ordering::Greater)?;
    // Equal at the top of the zf range: Mag(u64::MAX, EXPONENT_MAX) is exactly `big`.
    want(
        "big vs big",
        mag(u64::MAX, hi)?.cmp_exact_num(&neg(&big)),
        Ordering::Equal,
    )?;
    want(
        "tiny vs tiny",
        mag(TOP, lo)?.cmp_exact_num(&tiny),
        Ordering::Equal,
    )?;
    // cmp_exact_diff with radii at the i64 ends: wide gap decided by exponents, zero operand by cmp_exact_num.
    want(
        "max vs |big - tiny|",
        top.cmp_exact_diff(&big, &tiny),
        Ordering::Greater,
    )?;
    want(
        "min vs |big - tiny|",
        bottom.cmp_exact_diff(&big, &tiny),
        Ordering::Less,
    )?;
    want(
        "min vs |tiny - 0|",
        bottom.cmp_exact_diff(&tiny, &num(0)),
        Ordering::Less,
    )?;
    want(
        "min vs |0 - 0|",
        bottom.cmp_exact_diff(&num(0), &num(0)),
        Ordering::Greater,
    )?;
    // Last bit at EXPONENT_MAX: |big - (big - ulp)| = 2^(EXPONENT_MAX - 64) = 0.1 * 2^(EXPONENT_MAX - 63).
    let below = parts(u64::MAX.checked_sub(1).ok_or("m")?, hi)?;
    let ulp_exp = sub(hi, 63)?;
    want(
        "ulp vs |big - below|",
        mag(TOP, ulp_exp)?.cmp_exact_diff(&big, &below),
        Ordering::Equal,
    )?;
    want(
        "ulp vs |-below - -big|",
        mag(TOP, ulp_exp)?.cmp_exact_diff(&neg(&below), &neg(&big)),
        Ordering::Equal,
    )?;
    want(
        "ulp+ vs |big - below|",
        mag(TOP | 1, ulp_exp)?.cmp_exact_diff(&big, &below),
        Ordering::Greater,
    )?;
    want(
        "ulp- vs |big - below|",
        mag(u64::MAX, sub(ulp_exp, 1)?)?.cmp_exact_diff(&big, &below),
        Ordering::Less,
    )?;
    Ok(())
}

#[test]
fn s2_edges_huge_precision_last_bit() -> Result<(), String> {
    let eps = pow2(sub(0, DEEP)?)?; // 2^-100000
    let one = num(1);
    let x = one.add_full_prec(&eps); // 1 + 2^-100000
    let y = one.sub_full_prec(&eps); // 1 - 2^-100000
    let m1 = mag(TOP, 1)?; // 1
    want("1 vs 1+eps", m1.cmp_exact_num(&x), Ordering::Less)?;
    want("1 vs -(1+eps)", m1.cmp_exact_num(&neg(&x)), Ordering::Less)?;
    want(
        "1+2^-63 vs 1+eps",
        mag(TOP | 1, 1)?.cmp_exact_num(&x),
        Ordering::Greater,
    )?;
    want("1 vs 1-eps", m1.cmp_exact_num(&y), Ordering::Greater)?;
    want(
        "1-2^-64 vs 1-eps",
        mag(u64::MAX, 0)?.cmp_exact_num(&y),
        Ordering::Less,
    )?;
    // |x - 1| = 2^-100000 = 0.1 * 2^-99999: decided on the last bit.
    let e_eps = sub(1, DEEP)?;
    want(
        "eps vs |x - 1|",
        mag(TOP, e_eps)?.cmp_exact_diff(&x, &one),
        Ordering::Equal,
    )?;
    want(
        "eps vs |1 - x|",
        mag(TOP, e_eps)?.cmp_exact_diff(&one, &x),
        Ordering::Equal,
    )?;
    want(
        "eps+ vs |x - 1|",
        mag(TOP | 1, e_eps)?.cmp_exact_diff(&x, &one),
        Ordering::Greater,
    )?;
    want(
        "eps- vs |x - 1|",
        mag(u64::MAX, sub(e_eps, 1)?)?.cmp_exact_diff(&x, &one),
        Ordering::Less,
    )?;
    // |x - y| = 2^-99999 = 0.1 * 2^-99998.
    want(
        "2eps vs |x - y|",
        mag(TOP, add(e_eps, 1)?)?.cmp_exact_diff(&x, &y),
        Ordering::Equal,
    )?;
    // Equal huge-precision operands: |x - x| = 0.
    want(
        "0 vs |x - x|",
        Mag::zero().cmp_exact_diff(&x, &x),
        Ordering::Equal,
    )?;
    want(
        "min vs |x - x|",
        mag(TOP, i64::MIN)?.cmp_exact_diff(&x, &x),
        Ordering::Greater,
    )?;
    // Directed rounding with the sticky bit 100000 binades down.
    want_up_lo("1+eps", &x, &mag(TOP | 1, 1)?, &m1)?;
    want_up_lo("-(1+eps)", &neg(&x), &mag(TOP | 1, 1)?, &m1)?;
    want_up_lo("1-eps", &y, &m1, &mag(u64::MAX, 0)?)?;
    Ok(())
}

#[test]
fn s2_edges_cancellation() -> Result<(), String> {
    // Exact path, gap 0: |(2^200 + 2^-50) - 2^200| = 2^-50 = 0.1 * 2^-49.
    let base = pow2(200)?;
    let a = base.add_full_prec(&pow2(-50)?);
    let m = mag(TOP, -49)?;
    want(
        "2^-50 vs |a - base|",
        m.cmp_exact_diff(&a, &base),
        Ordering::Equal,
    )?;
    want(
        "2^-50 vs |base - a|",
        m.cmp_exact_diff(&base, &a),
        Ordering::Equal,
    )?;
    want(
        "2^-50 vs |-a - -base|",
        m.cmp_exact_diff(&neg(&a), &neg(&base)),
        Ordering::Equal,
    )?;
    want(
        "2^-50 vs |a + base|",
        m.cmp_exact_diff(&a, &neg(&base)),
        Ordering::Less,
    )?;
    // 62-bit integers one apart: |3*2^60 + 1 - 3*2^60| = 1.
    let p = num(0x3000_0000_0000_0001);
    let q = num(0x3000_0000_0000_0000);
    want(
        "1 vs |p - q|",
        mag(TOP, 1)?.cmp_exact_diff(&p, &q),
        Ordering::Equal,
    )?;
    want(
        "1- vs |p - q|",
        mag(u64::MAX, 0)?.cmp_exact_diff(&q, &p),
        Ordering::Less,
    )?;
    // Wide gap, the radius cancels the big operand: m = 2^1000 + 2^937 against |2^1000 - 2^-1000| < 2^1000.
    let big = pow2(1000)?;
    let small = pow2(-1000)?;
    let near = mag(TOP | 1, 1001)?;
    want(
        "2^1000+ vs |big - small|",
        near.cmp_exact_diff(&big, &small),
        Ordering::Greater,
    )?;
    want(
        "2^1000+ vs |-big - -small|",
        near.cmp_exact_diff(&neg(&big), &neg(&small)),
        Ordering::Greater,
    )?;
    // |2^1000 + 2^-1000| < 2^1000 + 2^937 too.
    want(
        "2^1000+ vs |big + small|",
        near.cmp_exact_diff(&big, &neg(&small)),
        Ordering::Greater,
    )?;
    // Wide gap, the radius cancels the big operand to within a factor 4 of the small one (the z = ur ± small branch):
    // big = 2^1000 + 2^-500 + 2^-501, m = 2^1000, small = 2^-500: |big - small| = 2^1000 + 2^-501 > m.
    let deep = big.add_full_prec(&pow2(-500)?).add_full_prec(&pow2(-501)?);
    let at = mag(TOP, 1001)?;
    want(
        "2^1000 vs |deep - 2^-500|",
        at.cmp_exact_diff(&deep, &pow2(-500)?),
        Ordering::Less,
    )?;
    // |deep - 3*2^-501| = 2^1000 exactly.
    let three = pow2(-500)?.add_full_prec(&pow2(-501)?);
    want(
        "2^1000 vs |deep - 3*2^-501|",
        at.cmp_exact_diff(&deep, &three),
        Ordering::Equal,
    )?;
    want(
        "2^1000 vs |3*2^-501 - deep|",
        at.cmp_exact_diff(&three, &deep),
        Ordering::Equal,
    )?;
    // |deep - 2^-499| = 2^1000 - 2^-501 < m.
    want(
        "2^1000 vs |deep - 2^-499|",
        at.cmp_exact_diff(&deep, &pow2(-499)?),
        Ordering::Greater,
    )?;
    Ok(())
}

#[test]
fn s2_edges_nan_inf_every_position() -> Result<(), String> {
    let nan = num(0).div(&num(0), 64, RoundingMode::ToEven);
    let pinf = num(1).div(&num(0), 64, RoundingMode::ToEven);
    let ninf = neg(&pinf);
    if !nan.is_nan() || !pinf.is_inf() || !ninf.is_inf() {
        return Err("NaN / Inf not built".to_string());
    }
    let radii = [
        Mag::zero(),
        mag(TOP, 1)?,
        mag(TOP, i64::MIN)?,
        mag(u64::MAX, i64::MAX)?,
        Mag::infinity(),
    ];
    for r in &radii {
        for (w, v) in [("nan", &nan), ("+inf", &pinf), ("-inf", &ninf)] {
            let got = [
                r.cmp_exact_num(v),
                r.cmp_exact_diff(v, &num(1)),
                r.cmp_exact_diff(&num(1), v),
                r.cmp_exact_diff(v, v),
            ];
            if got.iter().any(Option::is_some) {
                return Err(format!("{:?} with {}: want all None, got {:?}", r, w, got));
            }
        }
    }
    // Upper bound Infinity, lower bound zero, for NaN and both infinities.
    for (w, v) in [("nan", &nan), ("+inf", &pinf), ("-inf", &ninf)] {
        if Mag::from_exact_num(v) != Mag::infinity() || Mag::from_exact_num_lower(v) != Mag::zero()
        {
            return Err(format!(
                "{}: up {:?} lo {:?}",
                w,
                Mag::from_exact_num(v),
                Mag::from_exact_num_lower(v)
            ));
        }
    }
    Ok(())
}

#[test]
fn s2_edges_round_trips_and_rounding() -> Result<(), String> {
    let (hi, lo) = (EXPONENT_MAX, EXPONENT_MIN);
    // Mag -> ExactNum -> Mag is the identity across the zf range, and the image compares Equal.
    for (mant, e) in [
        (TOP, hi),
        (u64::MAX, hi),
        (TOP, lo),
        (u64::MAX, lo),
        (TOP | 1, 0),
        (u64::MAX, 1),
        (TOP, -63),
    ] {
        let m = mag(mant, e)?;
        let x = m
            .to_exact_num()
            .map_err(|err| format!("to_exact_num({:#x}@{}): {:?}", mant, e, err))?;
        want(
            &format!("{:#x}@{} vs image", mant, e),
            m.cmp_exact_num(&x),
            Ordering::Equal,
        )?;
        want_up_lo(&format!("{:#x}@{} round trip", mant, e), &x, &m, &m)?;
    }
    // ExactNum (64 bits) -> Mag -> ExactNum gives the same |value|.
    for (mant, e) in [(TOP, hi), (u64::MAX, lo), (0x9e37_79b9_7f4a_7c15, -12345)] {
        let x = neg(&parts(mant, e)?);
        let back = Mag::from_exact_num(&x)
            .to_exact_num()
            .map_err(|err| format!("back {:#x}@{}: {:?}", mant, e, err))?;
        if back.abs_cmp(&x).ok_or("abs_cmp")? != 0 {
            return Err(format!("{:#x}@{}: round trip changed the value", mant, e));
        }
    }
    // Directed-rounding boundaries around 1 (Mag(TOP, 1)):
    // bit 64 set (1 + 2^-63) fits: exact both ways.
    let fits = num(1).add_full_prec(&pow2(-63)?);
    want_up_lo("1+2^-63", &fits, &mag(TOP | 1, 1)?, &mag(TOP | 1, 1)?)?;
    // bit 65 set (1 + 2^-64) is the first dropped bit: up one ulp, down truncates.
    let half = num(1).add_full_prec(&pow2(-64)?);
    want_up_lo("1+2^-64", &half, &mag(TOP | 1, 1)?, &mag(TOP, 1)?)?;
    want_up_lo("-(1+2^-64)", &neg(&half), &mag(TOP | 1, 1)?, &mag(TOP, 1)?)?;
    // All 64 ones plus a dropped bit: up carries into the next binade, down keeps the ones.
    let carry = parts(u64::MAX, 0)?.add_full_prec(&pow2(-200)?);
    want_up_lo(
        "(1-2^-64)+2^-200",
        &carry,
        &mag(TOP, 1)?,
        &mag(u64::MAX, 0)?,
    )?;
    // The same carry at EXPONENT_MAX lands one binade above the zf range: up has no ExactNum image.
    let carry_hi = parts(u64::MAX, hi)?.add_full_prec(&parts(TOP, sub(hi, 100)?)?);
    let up_hi = mag(TOP, add(hi, 1)?)?;
    want_up_lo(
        "carry at EXPONENT_MAX",
        &carry_hi,
        &up_hi,
        &mag(u64::MAX, hi)?,
    )?;
    match up_hi.to_exact_num() {
        Err(BoundError::Representation {
            op: "Mag::to_exact_num",
            detail: "exponent above zf EXPONENT_MAX",
        }) => Ok(()),
        Err(e) => Err(format!("carry at EXPONENT_MAX: wrong error {:?}", e)),
        Ok(v) => Err(format!(
            "carry at EXPONENT_MAX: want an error, got {:?}",
            v.as_raw_parts()
        )),
    }
}

/// Every `Mag::to_exact_num` error detail (all `Representation`; zf's `from_raw_parts` errors are not reachable from a
/// valid `Mag`, see the report).
#[test]
fn s2_edges_to_exact_num_every_error() -> Result<(), String> {
    let (hi, lo) = (EXPONENT_MAX, EXPONENT_MIN);
    let cases = [
        (Mag::infinity(), "Infinity has no ExactNum value"),
        (mag(TOP, add(hi, 1)?)?, "exponent above zf EXPONENT_MAX"),
        (mag(u64::MAX, i64::MAX)?, "exponent above zf EXPONENT_MAX"),
        (
            mag(TOP | 1, sub(lo, 1)?)?,
            "below zf EXPONENT_MIN and not an exact subnormal",
        ),
        (
            mag(TOP, sub(lo, 64)?)?,
            "below zf EXPONENT_MIN and not an exact subnormal",
        ),
        (
            mag(TOP, i64::MIN)?,
            "below zf EXPONENT_MIN and not an exact subnormal",
        ),
    ];
    for (m, want_detail) in &cases {
        match m.to_exact_num() {
            Err(BoundError::Representation { op, detail }) => {
                if op != "Mag::to_exact_num" || detail != *want_detail {
                    return Err(format!(
                        "{:?}: got {} / {}, want {}",
                        m, op, detail, want_detail
                    ));
                }
            }
            Err(e) => return Err(format!("{:?}: wrong variant {:?}", m, e)),
            Ok(v) => {
                return Err(format!(
                    "{:?}: want an error, got {:?}",
                    m,
                    v.as_raw_parts()
                ))
            }
        }
    }
    // The last exact subnormal just above the error: Mag(TOP, EXPONENT_MIN - 63) = 2^(EXPONENT_MIN - 64) is fine.
    mag(TOP, sub(lo, 63)?)?
        .to_exact_num()
        .map(|_| ())
        .map_err(|e| format!("2^(EXPONENT_MIN - 64): {:?}", e))
}

/// `Constant::ALL` lists every constant once, in zenith-float `ConstantId` order (16.S2.6).
#[test]
fn s2_edges_constant_all_order() -> Result<(), String> {
    let order = [
        Constant::Pi,
        Constant::E,
        Constant::Ln2,
        Constant::Ln10,
        Constant::Sqrt2,
        Constant::Phi,
        Constant::EulerGamma,
    ];
    if Constant::ALL == order {
        Ok(())
    } else {
        Err(format!(
            "Constant::ALL = {:?}, want {:?}",
            Constant::ALL,
            order
        ))
    }
}
