//! Short MPFR differential: add / mul / sqrt under every IEEE rounding mode.

use crate::mpfr::common::assert_float_close;
use crate::mpfr::common::get_float_pair;
use crate::mpfr::common::reset_test_rng;
use crate::mpfr::common::test_zf_op_no_cc;
use gmp_mpfr_sys::mpfr::{self, rnd_t};
use rug::float::{exp_max, exp_min};
use rug::Float;
use zenith_float_num::{
    Consts, ExactNum, Exponent, RoundingMode, EXPONENT_MAX, EXPONENT_MIN, WORD_BIT_SIZE,
};

const ROUNDS: [(RoundingMode, rnd_t); 5] = [
    (RoundingMode::ToEven, rnd_t::RNDN),
    (RoundingMode::Up, rnd_t::RNDU),
    (RoundingMode::Down, rnd_t::RNDD),
    (RoundingMode::FromZero, rnd_t::RNDA),
    (RoundingMode::ToZero, rnd_t::RNDZ),
];

#[test]
fn fuzz_mpfr_bit_exact_all_round_modes() {
    reset_test_rng();
    let mut cc = Consts::new().unwrap();

    unsafe {
        mpfr::set_emin(EXPONENT_MIN as gmp_mpfr_sys::gmp::exp_t);
        mpfr::set_emax(EXPONENT_MAX as gmp_mpfr_sys::gmp::exp_t);
    }
    assert_eq!(EXPONENT_MIN, exp_min());
    assert_eq!(EXPONENT_MAX, exp_max());

    let p = 2 * WORD_BIT_SIZE;
    let emin: Exponent = -64;
    let emax: Exponent = 64;

    for i in 0..24 {
        let (n1, f1) = get_float_pair(p, emin, emax, &mut cc);
        let (n2, f2) = get_float_pair(p, emin, emax, &mut cc);
        for (rm, rnd) in ROUNDS {
            test_zf_op_no_cc!(
                true,
                n1,
                n2,
                add,
                f1,
                f2,
                add,
                p,
                rm,
                rnd,
                (i, &n1, &n2, p, rm, "fuzz add"),
                cc
            );
            test_zf_op_no_cc!(
                true,
                n1,
                n2,
                mul,
                f1,
                f2,
                mul,
                p,
                rm,
                rnd,
                (i, &n1, &n2, p, rm, "fuzz mul"),
                cc
            );
            if !n1.is_negative() {
                test_zf_op_no_cc!(
                    true,
                    n1,
                    sqrt,
                    f1,
                    sqrt,
                    p,
                    rm,
                    rnd,
                    (i, &n1, p, rm, "fuzz sqrt"),
                    cc
                );
            }
        }
        let _ = (f1, f2, i);
    }
}
