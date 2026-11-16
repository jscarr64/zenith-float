//! ZF-N01: a subnormal result is rounded once, not first to `p` bits and then onto the grid.
use super::ExactNumNumber;
use crate::defs::{Error, RoundingMode, Sign, EXPONENT_MIN, WORD_BIT_SIZE};
use crate::ext::{ExactNum, RawKind, RawParts};

fn half_plus_2_pow_neg(k: u32) -> Result<ExactNumNumber, Error> {
    // 1/2 + 2^{-k} = (2^{k-1} + 1) * 2^{-k}: 128-or-more-bit mantissa at exponent 0.
    assert!(k >= 1);
    let w = (k as usize).div_ceil(64).max(1);
    let mut m = vec![0u64; w];
    let top = k as usize - 1; // bit 2^{k-1} in a k-bit integer, stored in k bits (or a word multiple)
    m[top / 64] |= 1u64 << (top % 64);
    m[0] |= 1;
    let p = w * 64;
    let sh = p - k as usize;
    if sh > 0 {
        // left-justify to a full word width so from_raw_parts sees a normal number
        let mut n = vec![0u64; w];
        let ws = sh / 64;
        let b = sh % 64;
        for (i, &x) in m.iter().enumerate() {
            if b == 0 {
                n[i + ws] = x;
            } else {
                n[i + ws] |= x << b;
                if i + ws + 1 < w {
                    n[i + ws + 1] |= x >> (64 - b);
                }
            }
        }
        m = n;
    }
    ExactNumNumber::from_raw_parts(&m, p, Sign::Pos, 0, false)
}

fn half_minus_2_pow_neg(k: u32) -> Result<ExactNumNumber, Error> {
    // 1/2 - 2^{-k} = (2^{k-1} - 1) * 2^{-k}, in [1/4, 1/2): exponent -1, full mantissa 1...10.
    assert!(k >= 2);
    let w = (k as usize).div_ceil(64).max(1);
    let p = w * 64;
    let mut m = vec![0u64; w];
    // integer 2^{k-1}-1, then left-justify into p bits (shift p-(k-1)), exponent -1.
    for i in 0..(k as usize - 1) {
        m[i / 64] |= 1u64 << (i % 64);
    }
    let sh = p - (k as usize - 1);
    if sh > 0 {
        let mut n = vec![0u64; w];
        let ws = sh / 64;
        let b = sh % 64;
        for (i, &x) in m.iter().enumerate() {
            if b == 0 {
                n[i + ws] = x;
            } else {
                n[i + ws] |= x << b;
                if i + ws + 1 < w {
                    n[i + ws + 1] |= x >> (64 - b);
                }
            }
        }
        m = n;
    }
    ExactNumNumber::from_raw_parts(&m, p, Sign::Pos, -1, false)
}

fn to_ext(v: &ExactNumNumber) -> ExactNum {
    ExactNum::from_raw_parts(&RawParts {
        kind: RawKind::Finite,
        words: v.mantissa().digits(),
        precision: v.mantissa_max_bit_len(),
        negative: v.is_negative(),
        exponent: v.exponent(),
        inexact: false,
        op_inexact: false,
    })
    .expect("finite")
}

fn is_minpos(v: &ExactNumNumber, p: usize) -> bool {
    !v.is_zero()
        && v.exponent() == EXPONENT_MIN
        && v.precision() == 1
        && v.mantissa().digits().first().copied() == Some(1)
        && v.mantissa_max_bit_len() == p
}

fn is_signed_zero(v: &ExactNumNumber, neg: bool) -> bool {
    v.is_zero() && v.is_negative() == neg && v.exponent() == 0
}

/// (1/2 + 2^{-128}) × minpos at 64 bits, ToEven, is minpos (above the midpoint), not 0.
#[test]
fn zfn01_mul_just_above_half_minpos_to_even_is_minpos() -> Result<(), Error> {
    let p = WORD_BIT_SIZE;
    let a = half_plus_2_pow_neg(128)?;
    let b = ExactNumNumber::min_positive(p)?;
    let got = a.mul(&b, p, RoundingMode::ToEven)?;
    assert!(
        is_minpos(&got, p) && got.is_positive(),
        "ExactNumNumber mul ToEven: got exponent {} zero={} bits={}",
        got.exponent(),
        got.is_zero(),
        got.precision()
    );
    assert!(got.inexact() || got.op_inexact());

    let ge = to_ext(&a).mul(&to_ext(&b), p, RoundingMode::ToEven);
    let want = ExactNum::min_positive(p);
    assert_eq!(ge, want, "ExactNum mul ToEven: got {ge:?}");
    Ok(())
}

/// (1/2 − 2^{-128}) × minpos at 64 bits, ToEven, is +0 (below the midpoint).
#[test]
fn zfn01_mul_just_below_half_minpos_to_even_is_zero() -> Result<(), Error> {
    let p = WORD_BIT_SIZE;
    let a = half_minus_2_pow_neg(128)?;
    let b = ExactNumNumber::min_positive(p)?;
    let got = a.mul(&b, p, RoundingMode::ToEven)?;
    assert!(
        is_signed_zero(&got, false),
        "ToEven below half: got exponent {} zero={} neg={}",
        got.exponent(),
        got.is_zero(),
        got.is_negative()
    );
    Ok(())
}

/// Directed modes and the other nearest-mode cases around half of minpos.
#[test]
fn zfn01_mul_half_minpos_nearest_modes() -> Result<(), Error> {
    let p = WORD_BIT_SIZE;
    let b = ExactNumNumber::min_positive(p)?;
    let above = half_plus_2_pow_neg(128)?;
    let below = half_minus_2_pow_neg(128)?;
    assert!(is_minpos(&above.mul(&b, p, RoundingMode::ToOdd)?, p));
    // Below the midpoint is not a tie: both nearest modes go to 0.
    assert!(is_signed_zero(
        &below.mul(&b, p, RoundingMode::ToOdd)?,
        false
    ));
    assert!(is_signed_zero(
        &below.mul(&b, p, RoundingMode::ToEven)?,
        false
    ));
    // Directed: anything in (0, minpos) goes to minpos away from zero / up, else 0.
    assert!(is_minpos(&above.mul(&b, p, RoundingMode::Up)?, p));
    assert!(is_minpos(&above.mul(&b, p, RoundingMode::FromZero)?, p));
    assert!(is_signed_zero(
        &above.mul(&b, p, RoundingMode::Down)?,
        false
    ));
    assert!(is_signed_zero(
        &above.mul(&b, p, RoundingMode::ToZero)?,
        false
    ));
    Ok(())
}

/// Division has the same double-rounding hole: a quotient just above half of minpos
/// must not flush to zero in ToEven.
#[test]
fn zfn01_div_just_above_half_minpos_to_even_is_minpos() -> Result<(), Error> {
    let p = WORD_BIT_SIZE;
    // (2^{64} + 2^{63} + 1) / 3 * 2^{q-64} = (1/2 + 2^{-66}) * 2^q, just above half of minpos.
    // Both operands are normals: N * 2^{q-64+136} / (3 * 2^{136}).
    let n: u128 = (1u128 << 64) + (1u128 << 63) + 1;
    let shifted = n << 63; // left-justify 65 bits in 128
    let m = [shifted as u64, (shifted >> 64) as u64];
    let e = EXPONENT_MIN + 73;
    let num = ExactNumNumber::from_raw_parts(&m, 128, Sign::Pos, e, false)?;
    let mut den = ExactNumNumber::from_word(3, p)?;
    den.set_exponent(138); // 3 * 2^{136}: from_word(3) has exponent 2
    let got = num.div(&den, p, RoundingMode::ToEven)?;
    assert!(
        is_minpos(&got, p),
        "div ToEven just above half: zero={} exp={} bits={}",
        got.is_zero(),
        got.exponent(),
        got.precision()
    );
    Ok(())
}

/// Public ExactNum signs of the original repro: ±(1/2 + 2^{-128}) × ±minpos, ToEven.
#[test]
fn zfn01_public_mul_signs_to_even() -> Result<(), Error> {
    let p = WORD_BIT_SIZE;
    let a = to_ext(&half_plus_2_pow_neg(128)?);
    let b = ExactNum::min_positive(p);
    let want = ExactNum::min_positive(p);
    for (sa, sb, neg) in [
        (false, false, false),
        (true, false, true),
        (false, true, true),
        (true, true, false),
    ] {
        let mut x = a.clone();
        let mut y = b.clone();
        if sa {
            x = -x;
        }
        if sb {
            y = -y;
        }
        let got = x.mul(&y, p, RoundingMode::ToEven);
        let mut w = want.clone();
        if neg {
            w = -w;
        }
        w.set_inexact(true);
        assert_eq!(got, w, "sign a={sa} b={sb}");
        assert!(got.inexact());
    }
    Ok(())
}

/// A flushed underflow in RoundingMode::None is a canonical zero (e = 0), still inexact.
#[test]
fn zfn01_mul_none_underflow_is_canonical_zero() -> Result<(), Error> {
    let p = WORD_BIT_SIZE;
    let a = ExactNumNumber::min_positive(p)?;
    let got = a.mul(&a, p, RoundingMode::None)?;
    assert!(
        is_signed_zero(&got, false),
        "None underflow: zero={} exp={} inexact={}",
        got.is_zero(),
        got.exponent(),
        got.inexact()
    );
    assert!(got.inexact() || got.op_inexact());
    Ok(())
}
