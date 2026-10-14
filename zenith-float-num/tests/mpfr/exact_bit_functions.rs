//! Functions beyond add / sub / mul / div round at exactly `p` bits in every mode, and their
//! `inexact()` reports that one operation (2.0.0 audit, round 2: items a and c).

use crate::mpfr::common::assert_exact_eq;
use crate::mpfr::common::conv_to_mpfr;
use crate::mpfr::common::mpfr_reference_with;
use crate::mpfr::common::reset_test_rng;
use crate::mpfr::common::test_random;
use gmp_mpfr_sys::mpfr::{self, rnd_t};
use rug::Float;
use zenith_float_num::{Consts, ExactNum, Exponent, RoundingMode, WORD_BIT_SIZE};

pub const PRECISIONS: [usize; 11] = [24, 53, 64, 65, 100, 113, 127, 128, 130, 200, 255];

pub const MODES: [RoundingMode; 6] = [
    RoundingMode::ToEven,
    RoundingMode::ToOdd,
    RoundingMode::Up,
    RoundingMode::Down,
    RoundingMode::ToZero,
    RoundingMode::FromZero,
];

/// A `p`-bit input in `2^elo ..= 2^ehi`, stored at the word precision zenith uses for `p`.
/// Its accumulated `inexact` flag is set so that a result cannot pass by inheriting it.
pub fn input(p: usize, elo: Exponent, ehi: Exponent, signed: bool) -> ExactNum {
    let storage = p.div_ceil(WORD_BIT_SIZE) * WORD_BIT_SIZE;
    let mut n = ExactNum::random_normal(storage, elo, ehi);
    n.set_precision(p, RoundingMode::ToZero).unwrap();
    // random_normal picks the sign at random
    if n.is_negative() != (signed && test_random::<bool>()) {
        n = n.neg();
    }
    n.set_inexact(true);
    n
}

/// `p`-bit value of a small integer ratio `num / 2^k`, stored at the word precision.
pub fn exact_input(num: i64, k: i32, p: usize) -> ExactNum {
    let storage = p.div_ceil(WORD_BIT_SIZE) * WORD_BIT_SIZE;
    let mut n = ExactNum::from_i64(num, storage);
    n = n.ldexp(-k, storage, RoundingMode::ToEven);
    n.set_inexact(true);
    n
}

/// Checks one result against MPFR: value bit for bit and the per-operation flag.
pub fn check(
    got: &ExactNum,
    p: usize,
    rm: RoundingMode,
    info: &str,
    cc: &mut Consts,
    f: impl Fn(&mut Float, rnd_t) -> i32,
) {
    let (want, inexact) = mpfr_reference_with(p as u32, rm, f);
    assert_exact_eq(got, &want, p, info, cc);
    if !want.is_infinite() && !want.is_nan() {
        assert_eq!(got.inexact(), inexact, "{info}: inexact flag, MPFR {want}");
    }
}

/// Sweeps a unary function over the precisions and modes on random inputs in each exponent
/// band and on the given exact inputs (`num / 2^k`).
pub fn sweep_unary(
    name: &str,
    bands: &[(Exponent, Exponent, bool)],
    exact: &[(i64, i32)],
    per_band: usize,
    zf: impl Fn(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
    mp: impl Fn(&mut Float, &Float, rnd_t) -> i32,
) {
    reset_test_rng();
    let mut cc = Consts::new().unwrap();
    for p in PRECISIONS {
        let mut inputs = Vec::new();
        for &(elo, ehi, signed) in bands {
            for _ in 0..per_band {
                inputs.push(input(p, elo, ehi, signed));
            }
        }
        for &(num, k) in exact {
            inputs.push(exact_input(num, k, p));
        }
        for x in &inputs {
            let fx = conv_to_mpfr(p.div_ceil(WORD_BIT_SIZE) * WORD_BIT_SIZE, x, &mut cc);
            for rm in MODES {
                let got = zf(x, p, rm, &mut cc);
                let info = format!("{name} p={p} {rm:?} x={x:?}");
                check(&got, p, rm, &info, &mut cc, |r, rnd| mp(r, &fx, rnd));
            }
        }
    }
}

#[test]
fn exact_bit_sqrt() {
    sweep_unary(
        "sqrt",
        &[(-40, 40, false), (-3, 3, false)],
        &[(1, 0), (4, 0), (9, 0), (1, 2), (49, 6), (2, 0), (3, 0)],
        12,
        |x, p, rm, _| x.sqrt(p, rm),
        |r, a, rnd| unsafe { mpfr::sqrt(r.as_raw_mut(), a.as_raw(), rnd) },
    );
}
