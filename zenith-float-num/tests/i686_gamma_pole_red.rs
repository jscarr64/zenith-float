//! Red (i686): gamma near the pole at -1 loses ~k bits to cancellation when
//! x = -(1 - 2^-k), but its guard bits are sized in words (`p_wrk + 2 *
//! WORD_BIT_SIZE` in gamma, `p + WORD_BIT_SIZE` in gamma_positive_capped).
//! That is 128 guard bits on 64-bit and 64 on i686, so at k = 160 the i686
//! result is wrong from bit ~170 and at k = 192 even the exponent is wrong;
//! Ziv accepts it. Same mechanism behind the s9_gamma / s9_digamma /
//! s9_ln_gamma 1-ulp misses at p=256 near -1 and near 1.
use zenith_float_num::{Consts, ExactNum, RoundingMode, Sign, Word};

fn near_minus_one(k: usize) -> ExactNum {
    let w = usize::try_from(Word::BITS).unwrap_or(64);
    let n = 256 / w;
    let mut m = vec![0; n];
    for b in 0..k {
        let pos = n * w - 1 - b;
        m[pos / w] |= Word::MAX ^ (Word::MAX >> 1) >> (w - 1 - pos % w);
    }
    ExactNum::from_words(&m, Sign::Neg, 0)
}

#[test]
fn gamma_near_pole_matches_high_precision() -> Result<(), String> {
    let mut cc = Consts::new().map_err(|e| format!("{e:?}"))?;
    for k in [128, 160, 192] {
        let x = near_minus_one(k);
        let lo = x.gamma(256, RoundingMode::ToEven, &mut cc);
        let hi = x.gamma(1024, RoundingMode::ToEven, &mut cc);
        let (Some(el), Some(eh)) = (lo.exponent(), hi.exponent()) else {
            return Err(format!("k={k}: non-finite {lo:?} / {hi:?}"));
        };
        if el != eh {
            return Err(format!("k={k}: exponent {el} vs reference {eh}"));
        }
        let d = lo.sub(&hi, 1024, RoundingMode::ToEven);
        if let Some(ed) = d.exponent() {
            if ed > eh - 250 {
                return Err(format!("k={k}: |gamma_256 - gamma_1024| ~ 2^{ed}, value ~ 2^{eh}"));
            }
        }
    }
    Ok(())
}
