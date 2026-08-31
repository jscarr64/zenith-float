//! Widen software IEEE values to `ExactNum` and IEEE-round back.

use super::arith::{pack_finite, pack_inf, pack_qnan, unpack, Class, Format};
use crate::defs::{RoundingMode, Word, WORD_BIT_SIZE};
use crate::ExactNum;

pub(super) fn to_exact(bits: u64, f: Format, p: usize) -> ExactNum {
    let u = unpack(bits, f);
    match u.class {
        Class::Nan => crate::NAN,
        Class::Inf => {
            if u.sign {
                crate::INF_NEG
            } else {
                crate::INF_POS
            }
        }
        Class::Zero => {
            let z = ExactNum::from_u8(0, p);
            if u.sign {
                z.neg()
            } else {
                z
            }
        }
        Class::Norm | Class::Sub => {
            let mut x = ExactNum::from_u64(u.sig, p);
            let sh = u.exp - f.bias - f.frac as i32;
            x = x.ldexp(sh, p, RoundingMode::None);
            if u.sign {
                x.neg()
            } else {
                x
            }
        }
    }
}

pub(super) fn from_exact(x: &ExactNum, f: Format) -> u64 {
    if x.is_nan() {
        return pack_qnan(f);
    }
    if x.is_inf() {
        return pack_inf(x.is_negative(), f);
    }
    if x.is_zero() {
        return if x.is_negative() { f.sign_mask() } else { 0 };
    }
    let words = match x.mantissa_digits() {
        Some(w) if !w.is_empty() => w,
        _ => return 0,
    };
    let (top, sticky) = top64_sticky(words);
    if top == 0 {
        return if x.is_negative() { f.sign_mask() } else { 0 };
    }
    let lz = top.leading_zeros();
    let aligned = top << lz;
    let need = f.frac + 1;
    let frac_part = aligned >> (64 - need);
    let rest = aligned << need;
    let half = 1u64 << 63;
    let mut sig = frac_part as u128;
    let g = rest >= half;
    let st = sticky || (rest << 1) != 0;
    let e0 = x.exponent().unwrap_or(0);
    let unbiased = e0 - 1;
    let mut stored = unbiased + f.bias;
    if g && (st || (sig & 1) == 1) {
        sig += 1;
        if sig >= (1u128 << (f.frac + 1)) {
            sig >>= 1;
            stored += 1;
        }
    }
    pack_finite(x.is_negative(), stored, sig, f)
}

fn top64_sticky(words: &[Word]) -> (u64, bool) {
    let mut acc: u128 = 0;
    let mut got = 0u32;
    for w in words.iter().rev() {
        acc = (acc << WORD_BIT_SIZE) | u128::from(*w);
        got += WORD_BIT_SIZE as u32;
        if got >= 64 {
            break;
        }
    }
    let top = if got >= 64 { (acc >> (got - 64)) as u64 } else { acc as u64 };
    let take = (64u32).div_ceil(WORD_BIT_SIZE as u32) as usize;
    let sticky = words
        .iter()
        .rev()
        .skip(take.min(words.len()))
        .any(|&w| w != 0);
    (top, sticky)
}
