//! Red (i686): `exp_neg_scaled` bounds |x| by `Word::BITS - 1` (pow.rs:115),
//! which is 31 on i686, so x = -0xb17217f7d1cf79fc * 2^-2 (|x| ~ 2^61.5) takes
//! the "|x| >= 2^63, far below EMIN" sticky branch and returns +0 instead of
//! the subnormal golden `f:+:EMIN:0000000432cef8` (r9_exp, p=53, ToEven).
use zenith_float_num::{Consts, ExactNum, RoundingMode, Sign};

#[cfg(target_pointer_width = "64")]
const M: [zenith_float_num::Word; 1] = [0xb172_17f7_d1cf_79fc];
#[cfg(target_pointer_width = "32")]
const M: [zenith_float_num::Word; 2] = [0xd1cf_79fc, 0xb172_17f7];

#[test]
fn exp_into_subnormal_is_not_zero() -> Result<(), String> {
    let mut cc = Consts::new().map_err(|e| format!("{e:?}"))?;
    let x = ExactNum::from_words(&M, Sign::Neg, 62);
    let r = x.exp(53, RoundingMode::ToEven, &mut cc);
    if r.is_zero() || r.is_nan() {
        return Err(format!("exp({x:?}) underflowed to {r:?}; golden is a nonzero subnormal"));
    }
    Ok(())
}

/// Rows straddling |x| = 2^63 (the `Word::BITS - 1` cut in `exp_neg_scaled`), p = 53, ToEven.
/// Derivation: |x| >= 2^63 - 1024, so |x| * log2(e) > 1.44 * 2^62.9 > 2^62 > |EXPONENT_MIN| and
/// EXPONENT_MAX = 2^62 - 1 (defs.rs:49/52). Hence e^-|x| is far below the smallest subnormal and
/// rounds to +0, and e^|x| is above the largest finite value and rounds to +Inf (round to nearest;
/// IEEE 754-2019 7.4/7.5 overflow and underflow rules, the semantics MPFR mpfr_exp follows).
/// Inputs are exact 53-bit integers: 2^63 - 1024 and 2^63 + 2048.
#[test]
fn exp_rows_around_two_pow_63() -> Result<(), String> {
    let mut cc = Consts::new().map_err(|e| format!("{e:?}"))?;
    let mut bad = Vec::new();
    for s in ["9223372036854774784", "9223372036854777856"] {
        let x = ExactNum::parse(s, zenith_float_num::Radix::Dec, 64, RoundingMode::None, &mut cc);
        let neg = ExactNum::parse(&format!("-{s}"), zenith_float_num::Radix::Dec, 64, RoundingMode::None, &mut cc);
        let up = x.exp(53, RoundingMode::ToEven, &mut cc);
        let down = neg.exp(53, RoundingMode::ToEven, &mut cc);
        if !up.is_inf_pos() {
            bad.push(format!("exp({s}) = {up:?}, want +Inf"));
        }
        if !(down.is_zero() && down.is_positive()) {
            bad.push(format!("exp(-{s}) = {down:?}, want +0"));
        }
    }
    if bad.is_empty() { Ok(()) } else { Err(bad.join("\n")) }
}
