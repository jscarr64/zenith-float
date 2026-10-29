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
