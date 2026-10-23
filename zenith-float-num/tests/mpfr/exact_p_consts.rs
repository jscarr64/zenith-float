//! 2.0: the constants round at exactly `p` bits, as MPFR does, in all six modes (not at `p` rounded up to the
//! word size), and `from_mag_parts` gives exactly the value `mant * 2^(exp - 64)` when it fits in `p` bits and
//! `InvalidArgument` when it does not.

use crate::mpfr::common::{assert_exact_eq, mpfr_reference_with};
use crate::mpfr::exact_bit_functions::MODES;
use gmp_mpfr_sys::mpfr;
use rug::Float;
use zenith_float_num::{Consts, Error, ExactNum, RoundingMode};

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

#[test]
fn consts_round_at_exactly_p_bits() -> Res {
    let table: [(&str, ZfConst, fn(&mut Float, mpfr::rnd_t) -> i32); 4] = [
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

#[test]
fn from_mag_parts_is_exact_at_p_bits() -> Res {
    // (mant, significant bits)
    let mants: [u64; 5] = [
        1 << 63,
        3 << 62,
        0xF000_0000_0000_0000,
        0xFFFF_FFFF_FFFF_FFFF,
        0x8000_0000_0000_0001,
    ];
    let mut wrong = Vec::new();
    for p in PRECISIONS {
        for mant in mants {
            let bits = 64 - mant.trailing_zeros() as usize;
            let got = ExactNum::from_mag_parts(mant, 5, p);
            if bits <= p {
                let mut want = Float::with_val(p as u32, mant);
                want >>= 59;
                let ok = got.as_ref().ok().is_some_and(|x| {
                    let mut cc = Consts::new().ok();
                    cc.as_mut().is_some_and(|c| {
                        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                            assert_exact_eq(x, &want, p, "from_mag_parts", c).is_ok()
                        }))
                        .unwrap_or(false)
                    }) && x.precision() == Some(p)
                });
                if !ok {
                    wrong.push(format!("p={p} mant={mant:#x}: {got:?}"));
                }
            } else if got.as_ref().err() != Some(&Error::InvalidArgument) {
                wrong.push(format!(
                    "p={p} mant={mant:#x} ({bits} bits): want InvalidArgument, got {got:?}"
                ));
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
