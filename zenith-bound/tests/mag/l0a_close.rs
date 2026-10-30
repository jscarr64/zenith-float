//! L0a close-out (L0a-2, -3, -4, -6); oracle: dyadic-rational model (tests/model) and exact ExactNum arithmetic
//! COVERS: Mag::sqrt_up, Mag::sqrt_lower, Mag::mul_up, Mag::mul_lower, Mag::div_up, Mag::div_lower, Mag::add_up,
//! Mag::sub_lower, Mag::from_exact_num, Mag::from_exact_num_lower
use crate::model::direct;
use crate::model::dyadic::MagV;
use crate::support::{from_magv, to_magv};
use core::cmp::Ordering;
use zenith_bound::mag::Mag;
use zenith_float_num::{ExactNum, RoundingMode};

const TOP: u64 = 1 << 63;
const ONES: u64 = u64::MAX;
// 2^63 * sqrt(2), the correctly rounded-up 64-bit mantissa (mpmath 1.3.0, 200 digits)
const SQRT2_UP: u64 = 0xB504_F333_F9DE_6485;
const SQRT2_LOWER: u64 = 0xB504_F333_F9DE_6484;

fn fin(mant: u64, exp: i64) -> Result<Mag, String> {
    Mag::finite(mant, exp).map_err(|e| format!("finite({:#x}, {}): {:?}", mant, exp, e))
}

fn same(what: &str, got: &Mag, want: &Mag) -> Result<(), String> {
    if got.cmp(want) == Ordering::Equal && got.exp() == want.exp() && got.mant() == want.mant() {
        Ok(())
    } else {
        Err(format!(
            "{}: expected {:?}/{:?}, got {:?}/{:?}",
            what,
            want.mant(),
            want.exp(),
            got.mant(),
            got.exp()
        ))
    }
}

fn model(what: &str, got: &Mag, want: Option<MagV>) -> Result<(), String> {
    let want = want.ok_or_else(|| format!("{}: model gave no answer", what))?;
    let g = to_magv(got)?;
    if g == want {
        Ok(())
    } else {
        Err(format!("{}: model {:?}, got {:?}", what, want, g))
    }
}

// ---------------------------------------------------------------- L0a-2: exact and tie corners

#[test]
fn l0a2_sqrt_of_even_powers_of_two_is_exact() -> Result<(), String> {
    for k in [-1000_i64, -64, -1, 0, 1, 2, 63, 1000, 1 << 40] {
        let twok = k.checked_mul(2).ok_or("2k overflow")?;
        let x = fin(TOP, twok.checked_add(1).ok_or("exp overflow")?)?;
        let want = fin(TOP, k.checked_add(1).ok_or("exp overflow")?)?;
        same(&format!("sqrt_up(2^{})", twok), &x.sqrt_up(), &want)?;
        same(&format!("sqrt_lower(2^{})", twok), &x.sqrt_lower(), &want)?;
    }
    Ok(())
}

#[test]
fn l0a2_sqrt_of_odd_powers_of_two_brackets_sqrt2() -> Result<(), String> {
    for k in [-1000_i64, -1, 0, 1, 31, 1000] {
        let x = fin(
            TOP,
            k.checked_mul(2)
                .and_then(|v| v.checked_add(2))
                .ok_or("exp overflow")?,
        )?;
        let e = k.checked_add(1).ok_or("exp overflow")?;
        same(
            &format!("sqrt_up(2^(2*{}+1))", k),
            &x.sqrt_up(),
            &fin(SQRT2_UP, e)?,
        )?;
        same(
            &format!("sqrt_lower(2^(2*{}+1))", k),
            &x.sqrt_lower(),
            &fin(SQRT2_LOWER, e)?,
        )?;
    }
    Ok(())
}

#[test]
fn l0a2_all_ones_mantissas_match_the_model() -> Result<(), String> {
    for e in [-70_i64, -1, 0, 1, 64, 4000] {
        let a = MagV::Finite { exp: e, mant: ONES };
        let b = MagV::Finite { exp: 3, mant: ONES };
        let ma = from_magv(&a)?;
        let mb = from_magv(&b)?;
        model(
            &format!("sqrt_up(ONES@{})", e),
            &ma.sqrt_up(),
            direct::sqrt_up(&a),
        )?;
        model(
            &format!("sqrt_lower(ONES@{})", e),
            &ma.sqrt_lower(),
            direct::sqrt_lower(&a),
        )?;
        model(
            &format!("mul_up(ONES@{}, ONES@3)", e),
            &ma.mul_up(&mb),
            direct::mul_up(&a, &b),
        )?;
        model(
            &format!("mul_lower(ONES@{}, ONES@3)", e),
            &ma.mul_lower(&mb),
            direct::mul_lower(&a, &b),
        )?;
    }
    Ok(())
}

#[test]
fn l0a2_add_carry_at_i64_max() -> Result<(), String> {
    let top_max = fin(TOP, i64::MAX)?;
    if !top_max.add_up(&top_max).is_infinity() {
        return Err("TOP@MAX + TOP@MAX must carry out to Infinity".into());
    }
    let a = fin(ONES, i64::MAX - 1)?;
    let b = fin(TOP, i64::MAX - 64)?;
    same(
        "ONES@MAX-1 + TOP@MAX-64",
        &a.add_up(&b),
        &fin(TOP, i64::MAX)?,
    )
}

// ---------------------------------------------------------------- L0a-3: W1 translation and monotonicity

fn shift(m: &Mag, k: i64) -> Result<Mag, String> {
    match (m.mant(), m.exp()) {
        (Some(mant), Some(e)) if mant != 0 => fin(mant, e.checked_add(k).ok_or("shift overflow")?),
        _ if m.is_infinity() => Ok(Mag::infinity()),
        _ => Ok(Mag::zero()),
    }
}

fn grid() -> Result<Vec<Mag>, String> {
    let mut v = Vec::new();
    for mant in [
        TOP,
        TOP | 1,
        0xC000_0000_0000_0000,
        SQRT2_UP,
        ONES - 1,
        ONES,
    ] {
        for e in [-3_i64, 0, 1, 5, 64, 70] {
            v.push(fin(mant, e)?);
        }
    }
    Ok(v)
}

#[test]
fn l0a3_w1_translation_bit_for_bit() -> Result<(), String> {
    let g = grid()?;
    for k in [-(1_i64 << 61), -1000, -1, 1, 1000, 1 << 61] {
        for a in g.iter() {
            let ak = shift(a, k)?;
            let a2k = shift(a, k.checked_mul(2).ok_or("2k")?)?;
            same("sqrt_up W1", &a2k.sqrt_up(), &shift(&a.sqrt_up(), k)?)?;
            same(
                "sqrt_lower W1",
                &a2k.sqrt_lower(),
                &shift(&a.sqrt_lower(), k)?,
            )?;
            for b in g.iter() {
                let bk = shift(b, k)?;
                let b2 = shift(b, k.checked_neg().ok_or("neg")?)?;
                same("mul_up W1", &ak.mul_up(b), &shift(&a.mul_up(b), k)?)?;
                same(
                    "mul_lower W1",
                    &ak.mul_lower(b),
                    &shift(&a.mul_lower(b), k)?,
                )?;
                same(
                    "div_up W1",
                    &ak.div_up(&b2),
                    &shift(&shift(&a.div_up(b), k)?, k)?,
                )?;
                same(
                    "div_lower W1",
                    &ak.div_lower(&b2),
                    &shift(&shift(&a.div_lower(b), k)?, k)?,
                )?;
                same("add_up W1", &ak.add_up(&bk), &shift(&a.add_up(b), k)?)?;
                same(
                    "sub_lower W1",
                    &ak.sub_lower(&bk),
                    &shift(&a.sub_lower(b), k)?,
                )?;
            }
        }
    }
    Ok(())
}

#[test]
fn l0a3_monotone_in_each_argument() -> Result<(), String> {
    let mut g = grid()?;
    g.push(Mag::zero());
    g.push(Mag::infinity());
    g.sort();
    for (i, a1) in g.iter().enumerate() {
        for a2 in g.iter().skip(i) {
            for b in g.iter() {
                if a1.add_up(b) > a2.add_up(b) || a1.mul_up(b) > a2.mul_up(b) {
                    return Err(format!(
                        "add_up/mul_up not monotone at {:?} <= {:?}, b {:?}",
                        a1, a2, b
                    ));
                }
                if b.mant().is_some_and(|m| m != 0) {
                    if a1.div_up(b) > a2.div_up(b) {
                        return Err(format!(
                            "div_up not monotone in the dividend at {:?} <= {:?}",
                            a1, a2
                        ));
                    }
                    // a1, a2 as divisors: the quotient is non-increasing in the divisor
                    if a1.mant().is_some_and(|m| m != 0)
                        && a2.mant().is_some_and(|m| m != 0)
                        && b.div_up(a1) < b.div_up(a2)
                    {
                        return Err(format!(
                            "div_up not non-increasing in the divisor at {:?} <= {:?}",
                            a1, a2
                        ));
                    }
                }
            }
        }
    }
    Ok(())
}

// ---------------------------------------------------------------- L0a-4: P9 conversion enclosure

const P: usize = 512;

fn pow2(e: i64) -> Result<ExactNum, String> {
    ExactNum::pow2(e, P).map_err(|err| format!("pow2({}): {:?}", e, err))
}

fn exact(m: &Mag) -> Result<ExactNum, String> {
    m.to_exact_num()
        .map_err(|e| format!("to_exact_num({:?}): {:?}", m, e))
}

fn nonneg(what: &str, d: &ExactNum) -> Result<(), String> {
    if d.is_negative() || d.is_nan() {
        Err(format!("P9 violated: {}", what))
    } else {
        Ok(())
    }
}

/// tail -1: (2^64-1) * 2^e + 2^(e-200), whose round-up carries into the next exponent (mutant 190);
/// tail 0: (2^64-1) * 2^e, exact in a Mag; tail t > 0: 2^e + 2^(e-1) + 2^(e-t), which needs t+1 bits
fn sample(e: i64, tail: i64) -> Result<ExactNum, String> {
    if tail < 0 {
        let ones = ExactNum::from_word(u64::MAX, P).mul(&pow2(e)?, P, RoundingMode::ToEven);
        return Ok(ones.add(
            &pow2(e.checked_sub(200).ok_or("e-200")?)?,
            P,
            RoundingMode::ToEven,
        ));
    }
    if tail == 0 {
        return Ok(ExactNum::from_word(u64::MAX, P).mul(&pow2(e)?, P, RoundingMode::ToEven));
    }
    let x = pow2(e)?.add(
        &pow2(e.checked_sub(1).ok_or("e-1")?)?,
        P,
        RoundingMode::ToEven,
    );
    Ok(x.add(
        &pow2(e.checked_sub(tail).ok_or("tail")?)?,
        P,
        RoundingMode::ToEven,
    ))
}

#[test]
fn l0a4_p9_from_exact_num_encloses_abs() -> Result<(), String> {
    for e in [-300_i64, -65, -1, 0, 1, 63, 64, 200] {
        for tail in [-1_i64, 0, 1, 2, 63, 64, 65, 100, 300] {
            for neg in [false, true] {
                let x = sample(e, tail)?;
                let sx = if neg {
                    -sample(e, tail)?
                } else {
                    sample(e, tail)?
                };
                let up = exact(&Mag::from_exact_num(&sx))?;
                let lo = exact(&Mag::from_exact_num_lower(&sx))?;
                nonneg(
                    &format!(
                        "from_exact_num >= |x| (e {}, tail {}, neg {})",
                        e, tail, neg
                    ),
                    &up.sub(&x, P, RoundingMode::ToEven),
                )?;
                nonneg(
                    &format!(
                        "|x| >= from_exact_num_lower (e {}, tail {}, neg {})",
                        e, tail, neg
                    ),
                    &x.sub(&lo, P, RoundingMode::ToEven),
                )?;
            }
        }
    }
    Ok(())
}

// ---------------------------------------------------------------- L0a-6: trap regressions

#[test]
fn l0a6_second_operand_matters() -> Result<(), String> {
    // the golden runner once fed the first operand twice; these rows differ only in the second operand.
    // Killed mutants: 183 (mul_up ignores b), 184 (div ignores b), 185 (add_up ignores b)
    let a = fin(TOP, 1)?;
    let b = fin(0xC000_0000_0000_0000, 1)?;
    if a.mul_up(&a) == a.mul_up(&b) || a.add_up(&a) == a.add_up(&b) || a.div_up(&a) == a.div_up(&b)
    {
        return Err("an op ignored its second operand".into());
    }
    Ok(())
}

#[test]
fn l0a6_round_down_overflow_rows() -> Result<(), String> {
    // a lower bound must never overflow to Infinity: it clamps to the largest finite Mag
    let big = fin(ONES, i64::MAX)?;
    let two = fin(TOP, 2)?;
    let half = fin(TOP, 0)?;
    for (what, r) in [
        ("mul_lower(MAX, 2)", big.mul_lower(&two)),
        ("div_lower(MAX, 1/2)", big.div_lower(&half)),
    ] {
        if r.is_infinity() {
            return Err(format!("{} rounded down to Infinity", what));
        }
        same(what, &r, &big)?;
    }
    if !big.mul_up(&two).is_infinity() {
        return Err("mul_up(MAX, 2) must overflow to Infinity".into());
    }
    Ok(())
}
