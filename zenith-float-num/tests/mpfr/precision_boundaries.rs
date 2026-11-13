//! Rusty (2026-11-10): add, sub, mul and div at every precision from 1 to 130 and at the word
//! boundaries 191..=193, 255..=257, 511..=513 and 1023..=1025 (a multiple of the word size, one
//! either side; this covers 32, 64, 128, 256, 512 and 1024), in all 7 rounding modes, both signs,
//! against MPFR bit for bit with the per-operation flag. `None` (no final rounding: the result keeps
//! the storage precision) is checked against MPFR cut toward zero at the storage precision, both in the
//! 7-mode sweep and in its own rows below. These are the
//! branches of the exact-rounding helper (`p` a multiple of the word size: one call at `p`; otherwise a
//! word-wider pass toward zero and a round-to-odd step) whose width computation changed.

use crate::mpfr::common::{conv_to_mpfr, mpfr_reference, reset_test_rng, zenith_binop, BinOp};
use crate::mpfr::exact_bit_functions::{exact_input, input, MODES};
use rug::ops::AssignRound;
use rug::Float;
use zenith_float_num::{Consts, ExactNum, RoundingMode, WORD_BIT_SIZE};

fn storage(p: usize) -> usize {
    p.div_ceil(WORD_BIT_SIZE) * WORD_BIT_SIZE
}

#[test]
fn binary_ops_match_mpfr_at_every_small_precision_and_word_boundary(
) -> Result<(), Box<dyn std::error::Error>> {
    reset_test_rng();
    let mut cc = Consts::new()?;
    let ps = (1..=130)
        .chain(191..=193)
        .chain(255..=257)
        .chain(511..=513)
        .chain(1023..=1025);
    for p in ps {
        let mut v: Vec<(ExactNum, ExactNum)> = Vec::new();
        for _ in 0..3 {
            v.push((input(p, -20, 20, true)?, input(p, -20, 20, true)?));
        }
        for (a, b) in [(12, 4), (-12, 4), (12, -4), (-12, -4), (1, 3), (-1, 3)] {
            v.push((exact_input(a, 0, p), exact_input(b, 0, p)));
        }
        for (x, y) in &v {
            let fx = conv_to_mpfr(storage(p), x, &mut cc)?;
            let fy = conv_to_mpfr(storage(p), y, &mut cc)?;
            for op in [BinOp::Add, BinOp::Sub, BinOp::Mul, BinOp::Div] {
                for rm in MODES {
                    let got = zenith_binop(x, y, op, p, rm);
                    let want = mpfr_reference(&fx, &fy, op, u32::try_from(p)?, rm)?;
                    let info = format!("{op:?} p={p} {rm:?} {x:?} {y:?}");
                    crate::mpfr::common::assert_exact_eq(&got, &want, p, &info, &mut cc)?;
                }
            }
            // `None` skips the final rounding: the result keeps the storage precision
            // (`ExactNum::add` docs), i.e. it is cut toward zero at storage(p) bits.
            let sp = u32::try_from(storage(p))?;
            for op in [BinOp::Add, BinOp::Sub, BinOp::Mul, BinOp::Div] {
                let got = zenith_binop(x, y, op, p, RoundingMode::None);
                let mut want = Float::new(sp);
                let rz = rug::float::Round::Zero;
                match op {
                    BinOp::Add => want.assign_round(&fx + &fy, rz),
                    BinOp::Sub => want.assign_round(&fx - &fy, rz),
                    BinOp::Mul => want.assign_round(&fx * &fy, rz),
                    BinOp::Div => want.assign_round(&fx / &fy, rz),
                };
                let info = format!("{op:?} p={p} None {x:?} {y:?}");
                crate::mpfr::common::assert_exact_eq(&got, &want, p, &info, &mut cc)?;
            }
        }
    }
    Ok(())
}
