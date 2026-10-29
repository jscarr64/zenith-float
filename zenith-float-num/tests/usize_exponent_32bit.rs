//! Exponents past `usize` on a 32-bit target (`i686-unknown-linux-musl`): an i64 exponent, or an integer
//! derived from one, must not be truncated or saturated through `usize`. Every row also runs on the host,
//! where it already holds. Exact models: powers of two, and the binary/decimal exponents of `10^N` and
//! `2^K` from `bc -l` at scale 80 (`floor(5e9 * log2 10) = 16609640474`,
//! `floor(2^40 * log10 2) = 330985980541`).

use std::ops::Neg;
use zenith_float_num::{Consts, ExactNum, Radix, RoundingMode, EXPONENT_MIN};

const P: usize = 64;

fn cc() -> Result<Consts, String> {
    Consts::new().map_err(|e| format!("Consts::new: {e:?}"))
}

/// `2^k` exactly: 1 = 0.5 * 2^1, so the zf exponent is `k + 1`.
fn pow2(k: i64) -> ExactNum {
    ExactNum::from_u8(1, P).ldexp(k, P, RoundingMode::None)
}

#[test]
fn powsi_positive_beyond_u32() {
    let two = ExactNum::from_u8(2, P);
    for k in [(1i64 << 32) + 1, (1i64 << 33) + 5, 1i64 << 40] {
        let r = two.powsi(k, P, RoundingMode::Down);
        assert_eq!(r.exponent(), Some(k + 1), "2^{k} exponent");
        assert!(!r.inexact(), "2^{k} is exact");
        assert!(r == pow2(k), "2^{k}");
    }
    let m2 = ExactNum::from_u8(2, P).neg();
    assert!(!m2.powsi(1i64 << 32, P, RoundingMode::Down).is_negative(), "(-2)^(2^32) > 0");
    assert!(m2.powsi((1i64 << 32) + 1, P, RoundingMode::Down).is_negative(), "(-2)^(2^32+1) < 0");
}

#[test]
fn powsi_negative_beyond_u32() {
    let two = ExactNum::from_u8(2, P);
    for k in [(1i64 << 32) + 2, (1i64 << 33) + 6, 1i64 << 40] {
        let r = two.powsi(-k, P, RoundingMode::Down);
        assert_eq!(r.exponent(), Some(1 - k), "2^-{k} exponent");
        assert!(r == pow2(-k), "2^-{k}");
    }
    let m2 = ExactNum::from_u8(2, P).neg();
    assert!(!m2.powsi(-(1i64 << 32), P, RoundingMode::Down).is_negative(), "(-2)^-(2^32) > 0");
    assert!(m2.powsi(-((1i64 << 32) + 1), P, RoundingMode::Down).is_negative(), "(-2)^-(2^32+1) < 0");
}



#[test]
fn parse_decimal_exponent_beyond_u32() -> Result<(), String> {
    let mut cc = cc()?;
    let up = ExactNum::parse("1e5000000000", Radix::Dec, P, RoundingMode::ToEven, &mut cc);
    assert_eq!(up.exponent(), Some(16609640475), "1e5000000000");
    let dn = ExactNum::parse("1e-5000000000", Radix::Dec, P, RoundingMode::ToEven, &mut cc);
    assert_eq!(dn.exponent(), Some(-16609640474), "1e-5000000000");
    Ok(())
}


#[test]
fn to_decimal_exponent_beyond_u32() -> Result<(), String> {
    let mut cc = cc()?;
    for (k, d) in [(1i64 << 40, 330985980541i64), (-(1i64 << 40), -330985980542i64)] {
        let (_, digits, e) = pow2(k)
            .convert_to_radix(Radix::Dec, RoundingMode::ToEven, &mut cc)
            .map_err(|err| format!("2^{k} to dec: {err:?}"))?;
        let len = i64::try_from(digits.len()).map_err(|err| format!("{err}"))?;
        // `0.d0 d1 ... * 10^e` with `d0 != 0`, so `10^(e - 1) <= 2^k < 10^e`; record what came back
        assert_eq!(e - 1, d, "2^{k}: decimal exponent (len {len})");
    }
    Ok(())
}

#[test]
fn format_binary_exponent_beyond_u32() -> Result<(), String> {
    let mut cc = cc()?;
    for k in [1i64 << 40, (1i64 << 32) + 3] {
        let s = pow2(k)
            .format(Radix::Bin, RoundingMode::ToEven, &mut cc)
            .map_err(|err| format!("2^{k}: {err:?}"))?;
        assert!(s.ends_with(&format!("e+{k:b}")), "2^{k}: {s}");
        let s = pow2(-k)
            .format(Radix::Bin, RoundingMode::ToEven, &mut cc)
            .map_err(|err| format!("2^-{k}: {err:?}"))?;
        assert!(s.ends_with(&format!("e-{k:b}")), "2^-{k}: {s}");
    }
    Ok(())
}

#[test]
fn format_decimal_exponent_beyond_u32() -> Result<(), String> {
    let mut cc = cc()?;
    let s = pow2(1i64 << 40)
        .format(Radix::Dec, RoundingMode::ToEven, &mut cc)
        .map_err(|err| format!("{err:?}"))?;
    assert!(s.ends_with("e+330985980541"), "{s}");
    Ok(())
}

#[test]
fn format_subnormal_exponent() -> Result<(), String> {
    let mut cc = cc()?;
    let x = ExactNum::min_positive(P);
    let s = x
        .format(Radix::Bin, RoundingMode::ToEven, &mut cc)
        .map_err(|err| format!("{err:?}"))?;
    let tail = s.rsplit('-').next().ok_or("no exponent")?;
    let back = u64::from_str_radix(tail, 2).map_err(|err| format!("{s}: {err}"))?;
    assert_eq!(back, EXPONENT_MIN.unsigned_abs(), "{s}");
    Ok(())
}
