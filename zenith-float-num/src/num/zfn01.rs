//! ZF-N01: a subnormal result is rounded once, not first to `p` bits and then onto the grid.
use super::ExactNumNumber;
use crate::defs::{Error, Exponent, RoundingMode, Sign, EXPONENT_MIN, WORD_BIT_SIZE};
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

fn to_ext(v: &ExactNumNumber) -> Result<ExactNum, Error> {
    ExactNum::from_raw_parts(&RawParts {
        kind: RawKind::Finite,
        words: v.mantissa().digits(),
        precision: v.mantissa_max_bit_len(),
        negative: v.is_negative(),
        exponent: v.exponent(),
        inexact: false,
        op_inexact: false,
    })
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

    let ge = to_ext(&a)?.mul(&to_ext(&b)?, p, RoundingMode::ToEven);
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
    let a = to_ext(&half_plus_2_pow_neg(128)?)?;
    let b = ExactNum::min_positive(p);
    let want = ExactNum::min_positive(p);
    for (sa, sb, neg) in [
        (false, false, false),
        (true, false, true),
        (false, true, true),
        (true, true, false),
    ] {
        let (xn, yn);
        let x: &ExactNum = if sa {
            xn = a.neg_ref();
            &xn
        } else {
            &a
        };
        let y: &ExactNum = if sb {
            yn = b.neg_ref();
            &yn
        } else {
            &b
        };
        let got = ExactNum::mul_ref(x, y, p, RoundingMode::ToEven);
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

fn words_of(sig: u128, wbits: usize, shift: i128) -> Vec<u64> {
    let nwords = wbits.div_ceil(64);
    let mut v = sig;
    if shift > 0 {
        v <<= shift;
    } else if shift < 0 {
        v >>= (-shift) as u32;
    }
    let mut out = vec![0u64; nwords];
    for (i, item) in out.iter_mut().enumerate() {
        *item = (v >> (64 * i)) as u64;
    }
    out
}

fn bit_len_words(w: &[u64]) -> usize {
    for (i, &x) in w.iter().enumerate().rev() {
        if x != 0 {
            return i * 64 + (64 - x.leading_zeros() as usize);
        }
    }
    0
}

fn from_ae(sign: Sign, sig: u128, e: i128, p: usize) -> Result<ExactNumNumber, Error> {
    if sig == 0 {
        return ExactNumNumber::new(p);
    }
    let b = 128 - sig.leading_zeros() as i128;
    let lead = e + b;
    let width = p.max(b as usize).div_ceil(64) * 64;
    if lead >= i128::from(EXPONENT_MIN) {
        let w = words_of(sig, width, width as i128 - b);
        ExactNumNumber::from_raw_parts(&w, width, sign, lead as Exponent, false)
    } else {
        let sh = e - i128::from(EXPONENT_MIN) + width as i128;
        let w = words_of(sig, width, sh);
        let n = bit_len_words(&w);
        ExactNumNumber::from_raw_parts(&w, n, sign, EXPONENT_MIN, false)
    }
}

fn fmt_hypot(v: &ExactNumNumber) -> String {
    let s = if v.is_negative() { "-" } else { "+" };
    if v.is_zero() {
        return format!("Z:{s}");
    }
    let p = v.mantissa_max_bit_len();
    let (words, _n, _sg, e, _ix) = v.as_raw_parts();
    let nibbles = p / 4;
    let nwords = p.div_ceil(64);
    let mut hex = String::new();
    for i in (0..nwords).rev() {
        hex.push_str(&format!("{:016x}", words.get(i).copied().unwrap_or(0)));
    }
    format!("F:{s}:{e}:{}", &hex[hex.len() - nibbles..])
}

/// Independent one-rounding hypot: same scaling as `hypot`, then sqrt toward zero
/// one word wide and `place_scaled` once. The 693_504-pair search below locks
/// `hypot` to this path.
fn hypot_one_round(
    x: &ExactNumNumber,
    y: &ExactNumNumber,
    p: usize,
    mode: RoundingMode,
) -> Result<ExactNumNumber, Error> {
    let p = crate::common::util::round_p(p);
    let mut ax = x.abs()?;
    let mut ay = y.abs()?;
    if ax.abs_cmp(&ay) < 0 {
        core::mem::swap(&mut ax, &mut ay);
    }
    if ax.is_zero() {
        return ExactNumNumber::new2(p, Sign::Pos, x.inexact() | y.inexact());
    }
    if ay.is_zero() {
        ax.set_precision(p, mode)?;
        return Ok(ax);
    }
    let ax_exp = ax.exponent();
    let mut sx = ax;
    let ex = ax_exp - sx.normalize2() as i64;
    let ay_exp = ay.exponent();
    let mut sy = ay;
    let ey = ay_exp - sy.normalize2() as i64;
    let wx = sx.mantissa_max_bit_len() as i64;
    sx.set_exponent(0);
    let lim = wx.max(p as i64 / 2 + 1) + WORD_BIT_SIZE as i64;
    let gap = crate::exp::wide(ex) - crate::exp::wide(ey);
    let dy = if gap > crate::exp::wide(lim + WORD_BIT_SIZE as i64) { -lim } else { ey - ex };
    sy.set_exponent(dy as Exponent);
    let s2 = sx
        .mul_full_prec(&sx)?
        .add_full_prec(&sy.mul_full_prec(&sy)?)?;
    let wide = s2.sqrt(p + WORD_BIT_SIZE, RoundingMode::ToZero)?;
    let at = i128::from(wide.exponent()) + i128::from(ex);
    let dropped = wide.inexact() || wide.op_inexact();
    wide.place_scaled(at, dropped, p, mode)
}

/// 693_504 subnormal-band pairs (p=64, all seven modes, A,B = 0..=128) against the
/// one-rounding reference. This is the regression check for hypot's ZF-N01 fix.
const HYPOT_MODES: [RoundingMode; 7] = [
    RoundingMode::None,
    RoundingMode::Up,
    RoundingMode::Down,
    RoundingMode::ToZero,
    RoundingMode::FromZero,
    RoundingMode::ToEven,
    RoundingMode::ToOdd,
];

#[test]
fn zfn01_hypot_one_rounding_search() -> Result<(), Error> {
    let q = i128::from(EXPONENT_MIN) - 64;
    let mut bad = 0u32;
    let mut checked = 0u32;
    let modes = HYPOT_MODES;
    let p = 64usize;
    for mode in modes {
        for e1 in [q, q + 1, q + 8] {
            for e2 in [q, e1] {
                for a in 1u128..=128 {
                    for b in 0u128..=128 {
                        let x = from_ae(Sign::Pos, a, e1, p)?;
                        let y = from_ae(Sign::Pos, b, e2, p)?;
                        let got = x.hypot(&y, p, mode)?;
                        let want = hypot_one_round(&x, &y, p, mode)?;
                        checked += 1;
                        if fmt_hypot(&got) != fmt_hypot(&want) {
                            bad += 1;
                            if bad <= 8 {
                                eprintln!(
                                    "BAD {mode:?} hypot({a}p{e1},{b}p{e2}) got {} want {}",
                                    fmt_hypot(&got),
                                    fmt_hypot(&want)
                                );
                            }
                        }
                    }
                }
            }
        }
    }
    assert_eq!(checked, 693_504, "search size changed: {checked}");
    assert_eq!(bad, 0, "{bad} hypot mismatches vs one-rounding reference");
    Ok(())
}

/// Exact 3-4-5 at the subnormal scale: hypot(3 minpos, 4 minpos) is 5 minpos.
#[test]
fn zfn01_hypot_345_minpos_is_exact() -> Result<(), Error> {
    let p = WORD_BIT_SIZE;
    let q = i128::from(EXPONENT_MIN) - p as i128;
    let x = from_ae(Sign::Pos, 3, q, p)?;
    let y = from_ae(Sign::Pos, 4, q, p)?;
    let want = from_ae(Sign::Pos, 5, q, p)?;
    for rm in HYPOT_MODES {
        let got = x.hypot(&y, p, rm)?;
        assert_eq!(
            fmt_hypot(&got),
            fmt_hypot(&want),
            "3-4-5 minpos {rm:?}: got {} want {}",
            fmt_hypot(&got),
            fmt_hypot(&want)
        );
        assert!(
            !got.inexact() && !got.op_inexact(),
            "3-4-5 should be exact, {rm:?}"
        );
    }
    Ok(())
}
