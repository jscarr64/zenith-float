//! 2.0: the constants round at exactly `p` bits, as MPFR does, in all six modes (not at `p` rounded up to the
//! word size). `from_mag_parts` keeps its exact, never-rounding contract.

use crate::mpfr::common::{assert_exact_eq, mpfr_reference_with};
use crate::mpfr::exact_bit_functions::MODES;
use gmp_mpfr_sys::mpfr;
use rug::Float;
use zenith_float_num::{Consts, Error, ExactNum, RoundingMode, EXPONENT_MIN, WORD_BIT_SIZE};

type Res = Result<(), Box<dyn std::error::Error>>;

const PRECISIONS: [usize; 18] =
    [1, 2, 3, 7, 24, 53, 63, 64, 65, 100, 113, 127, 128, 130, 200, 255, 1000, 4097];

type ZfConst = fn(&mut Consts, usize, RoundingMode) -> ExactNum;

fn mp_pi(r: &mut Float, rnd: mpfr::rnd_t) -> i32 {
    unsafe { mpfr::const_pi(r.as_raw_mut(), rnd) }
}

fn mp_ln_2(r: &mut Float, rnd: mpfr::rnd_t) -> i32 {
    unsafe { mpfr::const_log2(r.as_raw_mut(), rnd) }
}

fn mp_e(r: &mut Float, rnd: mpfr::rnd_t) -> i32 {
    let one = Float::with_val(2, 1);
    unsafe { mpfr::exp(r.as_raw_mut(), one.as_raw(), rnd) }
}

fn mp_ln_10(r: &mut Float, rnd: mpfr::rnd_t) -> i32 {
    let ten = Float::with_val(8, 10);
    unsafe { mpfr::log(r.as_raw_mut(), ten.as_raw(), rnd) }
}

/// An MPFR constant (or constant-like value) at the precision of `r`.
type MpConst = fn(&mut Float, mpfr::rnd_t) -> i32;

#[test]
fn consts_round_at_exactly_p_bits() -> Res {
    let table: [(&str, ZfConst, MpConst); 4] = [
        ("pi", Consts::pi, mp_pi),
        ("ln_2", Consts::ln_2, mp_ln_2),
        ("e", Consts::e, mp_e),
        ("ln_10", Consts::ln_10, mp_ln_10),
    ];
    let mut cc = Consts::new()?;
    let mut wrong = Vec::new();
    for (name, zf, mp) in table {
        for p in PRECISIONS {
            for rm in MODES {
                let got = zf(&mut cc, p, rm);
                if rm == RoundingMode::ToOdd && p < 2 {
                    continue;
                }
                let (want, _) = mpfr_reference_with(p as u32, rm, mp)?;
                let info = format!("{name} p={p} {rm:?}");
                let same = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    let mut c2 = Consts::new().map_err(|e| format!("{e:?}"))?;
                    assert_exact_eq(&got, &want, p, &info, &mut c2).map_err(|e| e.to_string())
                }));
                if !matches!(same, Ok(Ok(()))) {
                    wrong.push(info);
                }
            }
        }
    }
    assert!(
        wrong.is_empty(),
        "{} wrong: {:?}",
        wrong.len(),
        &wrong[..wrong.len().min(12)]
    );
    Ok(())
}

/// `from_mag_parts` is an exact constructor (zenith-bound P1 to P3): it never rounds, stores the word-rounded
/// precision, puts the grid at `2^(EXPONENT_MIN - word_p)`, and its value round-trips exactly against MPFR.
#[test]
fn from_mag_parts_keeps_its_exact_contract() -> Res {
    let mants: [u64; 6] = [
        1 << 63,
        3 << 62,
        0xF000_0000_0000_0000,
        0xFFFF_FFFF_FFFF_FFFF,
        0x8000_0000_0000_0001,
        0xdead_beef_1234_5677,
    ];
    let mut cc = Consts::new()?;
    for p in PRECISIONS {
        let wp = p.max(64).div_ceil(WORD_BIT_SIZE) * WORD_BIT_SIZE;
        for mant in mants {
            for exp in [5i64, -70, 1000] {
                let info = format!("from_mag_parts({mant:#x}, {exp}, {p})");
                // never rounds: every bit of mant is kept, whatever p is
                let x =
                    ExactNum::from_mag_parts(mant, exp, p).map_err(|e| format!("{info}: {e:?}"))?;
                assert!(!x.inexact(), "{info}: inexact");
                assert_eq!(x.mantissa_max_bit_len(), Some(wp), "{info}: stored width");
                let mut want = Float::with_val(64, mant);
                if exp >= 64 {
                    want <<= (exp - 64) as u32;
                } else {
                    want >>= (64 - exp) as u32;
                }
                assert_exact_eq(&x, &want, 64, &info, &mut cc)?;
            }
        }
        // grid: the lowest value is 2^(EXPONENT_MIN - word_p) (exp = EXPONENT_MIN - word_p + 1); below it is rejected
        let edge = i64::from(EXPONENT_MIN) - wp as i64 + 1;
        let at =
            ExactNum::from_mag_parts(1 << 63, edge, p).map_err(|e| format!("grid p={p}: {e:?}"))?;
        assert_eq!(at.leading_bit_exponent(), Some(edge - 1), "grid p={p}");
        assert!(!at.inexact(), "grid p={p}: inexact");
        assert_eq!(
            ExactNum::from_mag_parts(1 << 63, edge - 1, p).err(),
            Some(Error::InvalidArgument),
            "below grid p={p}"
        );
    }
    Ok(())
}
