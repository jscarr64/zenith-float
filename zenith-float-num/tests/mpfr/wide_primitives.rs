//! The arithmetic primitives beyond `+ - * /` over the whole exponent range, bit for bit against MPFR:
//! `sqrt`, `cbrt`, `reciprocal`, `hypot`, `fma`, `atan2` and `nth_root` at operands from below the subnormal
//! band to the top of the range, at precisions from 1 bit to 200 (and a few at 20 000), in the six
//! rounding modes MPFR can reproduce (`ToOdd` through the `p + 1` bit construction of `mpfr_reference_with`).
//! Subnormal and overflowing results are compared too (MPFR's range is emulated by `mpfr_reference_with`), and so is
//! the per-operation inexact flag. Every mismatch of a run is listed at once.

use crate::mpfr::common::{assert_exact_eq, mpfr_reference_with};
use crate::mpfr::retry_budget::use_full_exponent_range;
use crate::mpfr::retry_budget::{round_of, ternary};
use gmp_mpfr_sys::mpfr::rnd_t;
use rug::float::Round;
use rug::ops::AssignRound;
use rug::Float;
use std::cmp::Ordering;
use std::ops::Neg;
use std::panic::{catch_unwind, AssertUnwindSafe};
use zenith_float_num::{
    Consts, ExactComplex, ExactNum, Exponent, RoundingMode, EXPONENT_MAX, EXPONENT_MIN,
};

type Res = Result<(), Box<dyn std::error::Error>>;

const MODES: [RoundingMode; 6] = [
    RoundingMode::ToEven,
    RoundingMode::ToOdd,
    RoundingMode::Up,
    RoundingMode::Down,
    RoundingMode::ToZero,
    RoundingMode::FromZero,
];

/// Result precisions: tiny (1 to 5 bits), around one and two words, and a few in between.
const PRECISIONS: [usize; 10] = [1, 2, 3, 5, 24, 53, 64, 65, 128, 200];

/// Operands are exact in this many bits (a 64-bit significand times a power of two).
const OPERAND_PRECISION: usize = 192;

const MIN: i64 = EXPONENT_MIN as i64;
const MAX: i64 = EXPONENT_MAX as i64;

/// `(significand, exponent)`: the operand is `significand * 2^exponent`. The exponents straddle the bottom of the
/// range (below `2^(EXPONENT_MIN - 1)` is subnormal) and reach the top (`k + 64 <= EXPONENT_MAX`).
fn grid(significands: &[u64], exponents: &[i64]) -> Vec<(u64, i64)> {
    let mut v = Vec::new();
    for &m in significands {
        for &k in exponents {
            v.push((m, k));
        }
    }
    v
}

fn unary_operands() -> Vec<(u64, i64)> {
    grid(
        &[1, 3, u64::MAX],
        &[MIN - 120, MIN - 21, MIN, MIN + 5, -1, 7, MAX - 130, MAX - 64],
    )
}

fn binary_operands() -> Vec<(u64, i64)> {
    grid(
        &[1, 0xdead_beef_1234_5677],
        &[MIN - 70, MIN + 5, 0, MAX - 130, MAX - 64],
    )
}

fn split_exponent(k: i64) -> (Exponent, Exponent) {
    let hi = k.max(MIN);
    ((hi as Exponent), (k - hi) as Exponent)
}

fn zf(m: u64, k: i64, negative: bool) -> ExactNum {
    zf_at(m, k, negative, OPERAND_PRECISION)
}

/// The same operand held in `p` bits of storage (the result is not changed by the extra bits; the kernels see a
/// longer mantissa).
fn zf_at(m: u64, k: i64, negative: bool, p: usize) -> ExactNum {
    let (hi, lo) = split_exponent(k);
    let x = ExactNum::from_word(m, p)
        .ldexp(hi, p, RoundingMode::ToEven)
        .ldexp(lo, p, RoundingMode::ToEven);
    if negative {
        x.neg()
    } else {
        x
    }
}

fn mp(m: u64, k: i64, negative: bool) -> Float {
    let (hi, lo) = split_exponent(k);
    let f = (Float::with_val(OPERAND_PRECISION as u32, m) << hi) << lo;
    if negative {
        -f
    } else {
        f
    }
}

/// The operand's exponent is inside the range MPFR emulates for a result of `p` bits.
fn in_mpfr_range(m: u64, k: i64, p: usize) -> bool {
    k + 64 - m.leading_zeros() as i64 > MIN - p as i64 + 1
}

/// Collects the mismatches of a run: each check may panic inside `assert_exact_eq`.
struct Run {
    cc: Consts,
    failures: Vec<String>,
    checked: usize,
    /// Compare the per-operation inexact flag as well.
    flags: bool,
}

impl Run {
    fn new(flags: bool) -> Result<Self, Box<dyn std::error::Error>> {
        use_full_exponent_range();
        Ok(Run {
            cc: Consts::new()?,
            failures: Vec::new(),
            checked: 0,
            flags,
        })
    }

    fn check(
        &mut self,
        got: &ExactNum,
        p: usize,
        rm: RoundingMode,
        info: &str,
        f: impl Fn(&mut Float, rnd_t) -> i32,
    ) {
        self.checked += 1;
        let (cc, flags) = (&mut self.cc, self.flags);
        let outcome = catch_unwind(AssertUnwindSafe(|| -> Res {
            let (want, inexact) = mpfr_reference_with(p as u32, rm, f)?;
            assert_exact_eq(got, &want, p, info, cc)?;
            if flags && !want.is_nan() && !want.is_infinite() {
                assert_eq!(got.inexact(), inexact, "{info}: inexact flag, MPFR {want}");
            }
            Ok(())
        }));
        let text = match outcome {
            Ok(Ok(())) => return,
            Ok(Err(e)) => format!("{info} p={p} {rm:?}: error {e}"),
            Err(payload) => match payload.downcast_ref::<String>() {
                Some(m) => format!("p={p} {rm:?} {m}"),
                None => format!("{info} p={p} {rm:?}: panic without a message"),
            },
        };
        self.failures.push(text.chars().take(260).collect());
    }

    fn finish(self, what: &str, at_least: usize) {
        assert!(
            self.checked >= at_least,
            "{what}: grid too small ({})",
            self.checked
        );
        assert!(
            self.failures.is_empty(),
            "{what}: {} of {} results differ from MPFR:\n{}",
            self.failures.len(),
            self.checked,
            self.failures.join("\n")
        );
    }
}

/// What a sweep skips: `(significand, exponent, negative, p)` of the operand and precision.
type Skip = fn(u64, i64, bool, usize) -> bool;

/// The operands, signs and precisions of a unary sweep, whether the inexact flag is compared, what is skipped.
struct Sweep<'a> {
    operands: &'a [(u64, i64)],
    signs: &'a [bool],
    precisions: &'a [usize],
    flags: bool,
    skip: Skip,
}

fn no_skip(_: u64, _: i64, _: bool, _: usize) -> bool {
    false
}

fn unary(
    what: &str,
    sweep: Sweep,
    zf_fn: impl Fn(&ExactNum, usize, RoundingMode) -> ExactNum,
    mp_fn: impl Fn(&mut Float, &Float, Round) -> Ordering,
) -> Res {
    let Sweep {
        operands,
        signs,
        precisions,
        flags,
        skip,
    } = sweep;
    let mut run = Run::new(flags)?;
    for &(m, k) in operands {
        for &neg in signs {
            let x = zf(m, k, neg);
            let fx = mp(m, k, neg);
            for &p in precisions {
                if !in_mpfr_range(m, k, p) || skip(m, k, neg, p) {
                    continue;
                }
                for rm in MODES {
                    let got = zf_fn(&x, p, rm);
                    let info = format!("{what}({m:#x} 2^{k}, neg {neg})");
                    run.check(&got, p, rm, &info, |r, rnd| {
                        ternary(mp_fn(r, &fx, round_of(rnd)))
                    });
                }
            }
        }
    }
    run.finish(what, 20);
    Ok(())
}

#[test]
fn sqrt_over_the_whole_range() -> Res {
    unary(
        "sqrt",
        Sweep {
            operands: &unary_operands(),
            signs: &[false],
            precisions: &PRECISIONS,
            flags: true,
            skip: no_skip,
        },
        |x, p, rm| x.sqrt(p, rm),
        |r, a, rd| r.assign_round(a.sqrt_ref(), rd),
    )
}

#[test]
fn cbrt_over_the_whole_range() -> Res {
    unary(
        "cbrt",
        Sweep {
            operands: &unary_operands(),
            signs: &[false, true],
            precisions: &PRECISIONS,
            flags: true,
            skip: no_skip,
        },
        |x, p, rm| x.cbrt(p, rm),
        |r, a, rd| r.assign_round(a.cbrt_ref(), rd),
    )
}

/// Precisions of the reciprocal: the usual ones and those where the Newton iteration starts (above two words).
const RECIPROCAL_PRECISIONS: [usize; 13] = [1, 2, 3, 5, 24, 53, 64, 65, 128, 129, 192, 200, 256];

/// `reciprocal` of every operand of the grid at every precision, subnormal operands included: above 128 bits
/// the Newton start read the top 64 bits of a subnormal's mantissa, which are zero, so `1 / 2^-2147483768` was
/// `Inf` in the directed modes where the largest finite number is the result (`div_ref(1, x)` was right).
#[test]
fn reciprocal_over_the_whole_range() -> Res {
    unary(
        "reciprocal",
        Sweep {
            operands: &unary_operands(),
            signs: &[false, true],
            precisions: &RECIPROCAL_PRECISIONS,
            flags: true,
            skip: no_skip,
        },
        |x, p, rm| x.reciprocal(p, rm),
        |r, a, rd| r.assign_round(a.recip_ref(), rd),
    )
}

/// Exact results at the ends of the range: perfect squares and cubes (the root is exact in every mode and not
/// inexact) and powers of two.
#[test]
fn exact_roots_and_reciprocals_at_the_ends() -> Res {
    let square = 0xffff_fffe_0000_0001u64; // (2^32 - 1)^2
    let cube = 1_000_003u64.pow(3);
    let sq_ops = grid(
        &[1, square],
        &[MIN - 120, MIN - 20, MIN, -2, 4, MAX - 130, MAX - 64],
    );
    let sq_ops: Vec<(u64, i64)> = sq_ops.into_iter().filter(|&(_, k)| k % 2 == 0).collect();
    unary(
        "sqrt exact",
        Sweep {
            operands: &sq_ops,
            signs: &[false],
            precisions: &[1, 2, 24, 53, 64, 128],
            flags: true,
            skip: no_skip,
        },
        |x, p, rm| x.sqrt(p, rm),
        |r, a, rd| r.assign_round(a.sqrt_ref(), rd),
    )?;
    let cu_ops = grid(
        &[1, cube],
        &[MIN - 120, MIN - 21, MIN + 6, -3, 9, MAX - 129, MAX - 66],
    );
    let cu_ops: Vec<(u64, i64)> = cu_ops.into_iter().filter(|&(_, k)| k % 3 == 0).collect();
    unary(
        "cbrt exact",
        Sweep {
            operands: &cu_ops,
            signs: &[false, true],
            precisions: &[1, 2, 24, 53, 64, 128],
            flags: true,
            skip: no_skip,
        },
        |x, p, rm| x.cbrt(p, rm),
        |r, a, rd| r.assign_round(a.cbrt_ref(), rd),
    )
}

#[test]
fn nth_root_over_the_whole_range() -> Res {
    for n in [3usize, 5, 7, 33] {
        let name = format!("nth_root n={n}");
        unary(
            &name,
            Sweep {
                operands: &unary_operands(),
                signs: &[false],
                precisions: &[1, 3, 24, 64, 65, 128],
                flags: true,
                skip: no_skip,
            },
            |x, p, rm| x.nth_root(n, p, rm),
            |r, a, rd| r.assign_round(a.root_ref(n as u32), rd),
        )?;
    }
    // exact: 3^40 = (3^8)^5
    let pow = 3u64.pow(40);
    let ops: Vec<(u64, i64)> = grid(
        &[1, pow],
        &[MIN - 120, MIN - 20, MIN + 5, 0, 10, MAX - 130, MAX - 65],
    )
    .into_iter()
    .filter(|&(_, k)| k % 5 == 0)
    .collect();
    unary(
        "nth_root 5 exact",
        Sweep {
            operands: &ops,
            signs: &[false],
            precisions: &[1, 2, 24, 64, 128],
            flags: true,
            skip: no_skip,
        },
        |x, p, rm| x.nth_root(5, p, rm),
        |r, a, rd| r.assign_round(a.root_ref(5), rd),
    )
}

/// An exact `nth_root` (`2^7` for `n = 7`, also `2^2147483583`) at 64 bits and more was flagged inexact in `ToEven`
/// and `ToOdd` (the value was right): the Newton value is exact only up to its working noise, so rounding it
/// cannot tell. At 24 bits it was not flagged.
#[test]
fn nth_root_exact_result_is_not_inexact() -> Res {
    let ops = [(1u64, 7i64), (1, 2_147_483_583)];
    unary(
        "nth_root 7 exact flag",
        Sweep {
            operands: &ops,
            signs: &[false],
            precisions: &[24, 64, 65, 128, 200],
            flags: true,
            skip: no_skip,
        },
        |x, p, rm| x.nth_root(7, p, rm),
        |r, a, rd| r.assign_round(a.root_ref(7), rd),
    )
}

fn binary(
    what: &str,
    operands: &[(u64, i64)],
    precisions: &[usize],
    skip_pair: fn((u64, i64), (u64, i64)) -> bool,
    zf_fn: impl Fn(&ExactNum, &ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum,
    mp_fn: impl Fn(&mut Float, &Float, &Float, Round) -> Ordering,
) -> Res {
    let mut run = Run::new(true)?;
    for &(m1, k1) in operands {
        for &(m2, k2) in operands {
            if skip_pair((m1, k1), (m2, k2)) {
                continue;
            }
            for (n1, n2) in [(false, false), (true, false), (false, true), (true, true)] {
                let (a, b) = (zf(m1, k1, n1), zf(m2, k2, n2));
                let (fa, fb) = (mp(m1, k1, n1), mp(m2, k2, n2));
                for &p in precisions {
                    if !in_mpfr_range(m1, k1, p) || !in_mpfr_range(m2, k2, p) {
                        continue;
                    }
                    for rm in MODES {
                        let got = zf_fn(&a, &b, p, rm, &mut run.cc);
                        let info =
                            format!("{what}({m1:#x} 2^{k1} neg {n1}, {m2:#x} 2^{k2} neg {n2})");
                        run.check(&got, p, rm, &info, |r, rnd| {
                            ternary(mp_fn(r, &fa, &fb, round_of(rnd)))
                        });
                    }
                }
            }
        }
    }
    run.finish(what, 20);
    Ok(())
}

fn no_pair_skipped(_: (u64, i64), _: (u64, i64)) -> bool {
    false
}

#[test]
fn hypot_over_the_whole_range() -> Res {
    binary(
        "hypot",
        &binary_operands(),
        &PRECISIONS,
        no_pair_skipped,
        |a, b, p, rm, _| a.hypot(b, p, rm),
        |r, a, b, rd| r.assign_round(a.hypot_ref(b), rd),
    )
}

/// `atan2` over the whole range, including extreme quotients: `y / x` beyond the largest number gives `+-pi/2` for
/// either sign of `x` (`atan2(2^2147483517, -2^-2147483643)`), an exact tiny quotient is rounded below itself by
/// the directed modes (`atan2(2^-2147483718, 1, 128, Up)`), and a quotient below the smallest number rounds to
/// zero or to the smallest number.
#[test]
fn atan2_over_the_whole_range() -> Res {
    binary(
        "atan2",
        &binary_operands(),
        &[1, 3, 24, 64, 128],
        no_pair_skipped,
        |y, x, p, rm, cc| y.atan2(x, p, rm, cc),
        |r, y, x, rd| r.assign_round(y.atan2_ref(x), rd),
    )
}

/// `fma(a, b, c)` with `c` the exact negated product (total cancellation, the sign of the zero by mode), the product
/// plus one part in 2^60 (a cancellation leaving a few bits), and an independent `c`; at the ends of the range.
/// `labels` selects the kinds of `c`, `factor_pairs` the factors.
fn fma_run(labels: &[&str], pairs: &[usize], precisions: &[usize], at_least: usize) -> Res {
    let mut run = Run::new(true)?;
    // factors whose product stays inside the range (the operands themselves reach both ends)
    let all_pairs: [((u64, i64), (u64, i64)); 8] = [
        ((1, MAX - 130), (1, MIN + 5)),
        ((0xdead_beef_1234_5677, MAX - 64), (3, MIN - 40)),
        ((u64::MAX, MIN - 70), (u64::MAX, MAX - 64)),
        ((3, 0), (5, -7)),
        ((0xffff_ffff_ffff_fff1, MIN + 5), (1, 90)),
        ((7, MIN - 21), (1, 21)),
        ((1, MAX - 70), (1, -MAX + 140)),
        ((0xdead_beef_1234_5677, -1), (u64::MAX, 3)),
    ];
    for &i in pairs {
        let ((m1, k1), (m2, k2)) = all_pairs[i];
        let (a, b) = (zf(m1, k1, false), zf(m2, k2, true));
        let (fa, fb) = (mp(m1, k1, false), mp(m2, k2, true));
        // c = -(a * b), exact in 192 bits (a 128-bit product)
        let c0 = ExactNum::mul_ref(&a, &b, OPERAND_PRECISION, RoundingMode::ToEven).neg();
        let c0_shift = ExactNum::mul_ref(&a, &b, OPERAND_PRECISION, RoundingMode::ToEven)
            .neg()
            .ldexp(-60, OPERAND_PRECISION, RoundingMode::ToEven);
        let c1 = c0_shift.add(&c0, OPERAND_PRECISION, RoundingMode::ToEven);
        let fc0 = -Float::with_val(OPERAND_PRECISION as u32, &fa * &fb);
        let fc0_shift = Float::with_val(OPERAND_PRECISION as u32, &fc0 >> 60);
        let fc1 = Float::with_val(OPERAND_PRECISION as u32, &fc0 + &fc0_shift);
        let (fc2, c2) = (mp(m2, k2 - 3, false), zf(m2, k2 - 3, false));
        // the product of the factors is exact in 256 bits
        let product = Float::with_val(256, &fa * &fb);
        for (label, c, fc) in
            [("-ab", &c0, &fc0), ("-ab(1+2^-60)", &c1, &fc1), ("indep", &c2, &fc2)]
        {
            if !labels.contains(&label) {
                continue;
            }
            for &p in precisions {
                if !in_mpfr_range(m1, k1, p) || !in_mpfr_range(m2, k2, p) {
                    continue;
                }
                for rm in MODES {
                    let got = a.fma(&b, c, p, rm);
                    let info = format!("fma({m1:#x} 2^{k1}, -{m2:#x} 2^{k2}, {label})");
                    run.check(&got, p, rm, &info, |r, rnd| {
                        ternary(r.assign_round(&product + fc, round_of(rnd)))
                    });
                }
            }
        }
    }
    run.finish("fma", at_least);
    Ok(())
}

const FMA_PRECISIONS: [usize; 8] = [1, 2, 24, 53, 64, 65, 128, 200];

#[test]
fn fma_over_the_whole_range() -> Res {
    // total and near cancellation, and an independent addend, at every pair
    fma_run(
        &["-ab", "-ab(1+2^-60)", "indep"],
        &[0, 1, 2, 3, 4, 5, 6, 7],
        &FMA_PRECISIONS,
        500,
    )
}

/// `fma` whose addend is far below (or far above) the product and of the opposite sign, so that the sum is just
/// below a power of two (`fma(1, -1, 2^-200)` at 128 bits is `-(1 - 2^-200)`: 128 ones for the directed modes
/// that round toward zero or up, `-1` otherwise). The stand-in for a negligible term used to lie right below the
/// other term's last bit, which is above the rounding position once `p` exceeds that term's width.
#[test]
fn fma_with_a_far_smaller_addend() -> Res {
    let mut run = Run::new(true)?;
    let one = zf(1, 0, false);
    let minus_one = zf(1, 0, true);
    let (f_one, f_minus_one) = (mp(1, 0, false), mp(1, 0, true));
    for gap in [100i64, 200, 1000, 100_000] {
        let c = zf(1, -gap, false);
        let fc = mp(1, -gap, false);
        for p in [24usize, 64, 65, 128, 200] {
            for rm in MODES {
                let got = one.fma(&minus_one, &c, p, rm);
                let info = format!("fma(1, -1, 2^-{gap})");
                run.check(&got, p, rm, &info, |r, rnd| {
                    ternary(r.assign_round(
                        &(Float::with_val(gap as u32 + 128, &f_one * &f_minus_one) + &fc),
                        round_of(rnd),
                    ))
                });
            }
        }
    }
    run.finish("fma(1, -1, tiny)", 20);
    Ok(())
}

/// `fma` with a subnormal product and a subnormal addend thousands of bits below it (the stand-in for the
/// negligible addend must not be placed below the smallest exponent: its exponent wrapped around to a huge
/// positive one and the result was `+3 * 2^2147480571` for `3 * 2^(EXPONENT_MIN - 19) * -7 + 3 * 2^(EXPONENT_MIN - 2997)`).
/// The addend is below MPFR's emulated range for the result, which does not matter to MPFR's addition.
#[test]
fn fma_with_subnormal_terms_far_apart() -> Res {
    let mut run = Run::new(true)?;
    for (ka, kc, neg_c) in [
        (MIN - 19, MIN - 2997, false),
        (MIN - 19, MIN - 2997, true),
        (MIN + 3, MIN - 5000, false),
    ] {
        let (a, b, c) = (zf(3, ka, false), zf(7, 0, true), zf_at(3, kc, neg_c, 8192));
        // MPFR's range here only reaches 1024 below the smallest exponent: form the exact sum 6000 binary places
        // higher (both terms exact, the gap is below 8192 bits), then move it down (exact) and round it once.
        let up = 6000;
        let (fa, fb, fc) = (mp(3, ka + up, false), mp(7, 0, true), mp(3, kc + up, neg_c));
        let prod = Float::with_val(256, &fa * &fb);
        let exact = Float::with_val(8192, &prod + &fc) >> (up as i32);
        for p in [24usize, 64, 128, 200] {
            for rm in MODES {
                let got = a.fma(&b, &c, p, rm);
                let info = format!("fma(3 2^{ka}, -7, +-3 2^{kc}) neg_c {neg_c}");
                run.check(&got, p, rm, &info, |r, rnd| {
                    ternary(r.assign_round(&exact, round_of(rnd)))
                });
            }
        }
    }
    run.finish("fma subnormal terms", 20);
    Ok(())
}

/// A few huge precisions (20 000 bits): the same functions at operands near both ends, a handful of cases each.
#[test]
fn primitives_at_20000_bits() -> Res {
    let ops = grid(&[3, u64::MAX], &[MIN - 20, MIN + 5, 0, MAX - 64]);
    unary(
        "sqrt 20000",
        Sweep {
            operands: &ops,
            signs: &[false],
            precisions: &[20_000],
            flags: true,
            skip: no_skip,
        },
        |x, p, rm| x.sqrt(p, rm),
        |r, a, rd| r.assign_round(a.sqrt_ref(), rd),
    )?;
    unary(
        "reciprocal 20000",
        Sweep {
            operands: &ops,
            signs: &[true],
            precisions: &[20_000],
            flags: true,
            skip: no_skip,
        },
        |x, p, rm| x.reciprocal(p, rm),
        |r, a, rd| r.assign_round(a.recip_ref(), rd),
    )?;
    unary(
        "cbrt 20000",
        Sweep {
            operands: &ops,
            signs: &[false],
            precisions: &[20_000],
            flags: true,
            skip: no_skip,
        },
        |x, p, rm| x.cbrt(p, rm),
        |r, a, rd| r.assign_round(a.cbrt_ref(), rd),
    )
}

// ---- W1b: regressions for four wrong answers found in published 1.0.12 --------------------------------------------
//
// 1.0.12 returned a wrong value for each of these; the branch's exact-`p`-bit kernels were expected to be right, and
// these tests keep them so (each is checked against MPFR in the six modes, flags included).

/// Operands of the composite-`n` root sweep: the end-of-range grid plus ordinary operands (the 1.0.12 error was at
/// ordinary exponents: it rounded after each square or cube root).
fn composite_root_operands() -> Vec<(u64, i64)> {
    let mut v = unary_operands();
    v.extend(grid(
        &[1, 3, 5, 7, 0xdead_beef_1234_5677, u64::MAX],
        &[-70, -9, -2, 0, 1, 3, 10, 1000, -1000],
    ));
    v
}

/// (a) `nth_root` with a composite `n` (4, 6, 8, 9, 10, 12, 15), whose 1.0.12 implementation took a square or cube
/// root and rounded after each step (1 ulp wrong in 8 to 27 percent of cases there). On this branch the root is
/// computed once and rounded once: the branch is already right, this test keeps it so.
#[test]
fn nth_root_with_a_composite_n() -> Res {
    for n in [4usize, 6, 8, 9, 10, 12, 15] {
        let name = format!("nth_root n={n}");
        // negative operands only for odd n (an even root of a negative is NaN, checked below)
        let signs: &[bool] = if n % 2 == 1 { &[false, true] } else { &[false] };
        unary(
            &name,
            Sweep {
                operands: &composite_root_operands(),
                signs,
                precisions: &[1, 2, 3, 5, 24, 53, 64, 65, 128, 200],
                flags: true,
                skip: no_skip,
            },
            |x, p, rm| x.nth_root(n, p, rm),
            |r, a, rd| r.assign_round(a.root_ref(n as u32), rd),
        )?;
        // exact roots: 3^(n*j) with the exponent a multiple of n
        let j = 40 / n;
        let pow = 3u64.pow((n * j) as u32);
        let n_i = n as i64;
        let ops: Vec<(u64, i64)> = grid(
            &[1, pow],
            &[MIN - 120, MIN - 20, MIN + 5, 0, 10 * n_i, MAX - 130, MAX - 65],
        )
        .into_iter()
        .filter(|&(_, k)| k % n_i == 0)
        .collect();
        let name = format!("nth_root n={n} exact");
        unary(
            &name,
            Sweep {
                operands: &ops,
                signs: &[false],
                precisions: &[1, 2, 24, 64, 128],
                flags: true,
                skip: no_skip,
            },
            |x, p, rm| x.nth_root(n, p, rm),
            |r, a, rd| r.assign_round(a.root_ref(n as u32), rd),
        )?;
    }
    Ok(())
}

/// (b) `hypot(1, 2^-k)` for `k >= 9887`: 1.0.12 returned `NaN` there (the scaling squared the small operand into
/// the underflow range). The branch is right; every mode and precision is compared with MPFR and must not be NaN.
#[test]
fn hypot_of_one_and_a_tiny_number_is_not_nan() -> Res {
    let ks: [i64; 10] = [
        9_885,
        9_886,
        9_887,
        9_888,
        12_000,
        100_000,
        5_000_000,
        2_000_000_000,
        2_147_483_000,
        MIN + 1,
    ];
    let mut run = Run::new(true)?;
    for &k in &ks {
        let (one, tiny) = (zf(1, 0, false), zf(1, -k, false));
        let (fone, ftiny) = (mp(1, 0, false), mp(1, -k, false));
        for &p in &[2usize, 24, 64, 128, 200] {
            for rm in MODES {
                let got = one.hypot(&tiny, p, rm);
                assert!(!got.is_nan(), "hypot(1, 2^-{k}, {p}, {rm:?}) is NaN");
                let info = format!("hypot(1, 2^-{k})");
                run.check(&got, p, rm, &info, |r, rnd| {
                    ternary(r.assign_round(fone.hypot_ref(&ftiny), round_of(rnd)))
                });
            }
        }
    }
    run.finish("hypot(1, tiny)", 20);
    Ok(())
}

/// (c) `reciprocal` of a tiny negative operand `-2^-2147483649 .. -2^-2147483647`: the quotient overflows and the
/// sign must be kept (`-Inf`, or the most negative finite number in `Up`/`ToZero`): 1.0.12 returned `+Inf` at
/// `p >= 129` in the nearest modes. The branch is right.
#[test]
fn reciprocal_of_a_tiny_negative_operand_keeps_its_sign() -> Res {
    let ops = [(1u64, MIN - 1), (1, MIN), (1, MIN + 1), (3, MIN - 2), (u64::MAX, MIN - 64)];
    unary(
        "reciprocal tiny negative",
        Sweep {
            operands: &ops,
            signs: &[true],
            precisions: &[24, 64, 128, 129, 192, 200, 256, 512],
            flags: true,
            skip: no_skip,
        },
        |x, p, rm| x.reciprocal(p, rm),
        |r, a, rd| r.assign_round(a.recip_ref(), rd),
    )?;
    // explicit: the nearest modes give -Inf at p >= 129, never +Inf
    let mut seen = 0;
    for k in [MIN - 1, MIN, MIN + 1] {
        let x = zf(1, k, true);
        for p in [129usize, 192, 256] {
            for rm in [
                RoundingMode::ToEven,
                RoundingMode::ToOdd,
                RoundingMode::Down,
                RoundingMode::FromZero,
            ] {
                let got = x.reciprocal(p, rm);
                assert!(got.is_inf_neg(), "1 / -2^{k} at {p} bits {rm:?}: {got:?}");
                seen += 1;
            }
        }
    }
    assert_eq!(seen, 36);
    Ok(())
}

/// (d) `atan2` of a tiny quotient, `atan2(2^-2147483718, 1)`: 1.0.12 flushed it to zero, an earlier build of this
/// branch gave `NaN(PrecisionRetryExhausted)`. The result is the quotient itself (nonzero, not NaN) once `p` is
/// large enough to hold it, and MPFR's value in every mode and precision.
#[test]
fn atan2_of_a_tiny_quotient_is_neither_zero_nor_nan() -> Res {
    let ys = [(1u64, MIN - 70), (3, MIN - 70), (1, MIN - 20), (1, MIN), (1, -1000)];
    let mut run = Run::new(true)?;
    let mut cc = Consts::new()?;
    for &(m, k) in &ys {
        for neg in [false, true] {
            let y = zf(m, k, neg);
            let fy = mp(m, k, neg);
            let (one, fone) = (zf(1, 0, false), mp(1, 0, false));
            for &p in &[24usize, 64, 128, 200] {
                if !in_mpfr_range(m, k, p) {
                    continue;
                }
                for rm in MODES {
                    let got = y.atan2(&one, p, rm, &mut cc);
                    assert!(!got.is_nan(), "atan2({m:#x} 2^{k}, 1, {p}, {rm:?}) is NaN");
                    // the quotient has at most 2 significant bits, so from 128 bits (grid 2^(MIN-128)) on it is representable
                    if k >= MIN - 70 && p >= 128 {
                        assert!(
                            !got.is_zero(),
                            "atan2({m:#x} 2^{k}, 1, {p}, {rm:?}) is zero"
                        );
                    }
                    let info = format!("atan2({m:#x} 2^{k} neg {neg}, 1)");
                    run.check(&got, p, rm, &info, |r, rnd| {
                        ternary(r.assign_round(fy.atan2_ref(&fone), round_of(rnd)))
                    });
                }
            }
        }
    }
    run.finish("atan2 tiny quotient", 20);
    Ok(())
}

/// The quotients users reach through `eval` (Accumath, 256-bit operands): `atan2(1e-3000, 1)` (a quotient of about
/// 2^-9966, `NaN` in 1.0.12), `atan2(exp(-1e9), 1)` and `atan2(exp(1e9), -exp(-1e9))` (a quotient beyond 2^2^31, `pi`
/// instead of `pi/2`). Operands are of the size `exp(+-1e9)` (2^+-1442695040), `1e+-3000` (2^+-9966) and one, in
/// both signs, at 256 bits as well.
#[test]
fn atan2_at_the_scales_users_reach() -> Res {
    const K: i64 = 1_442_695_040;
    binary(
        "atan2 user scales",
        &grid(&[1, 0xdead_beef_1234_5677], &[K, -K, 9_966, -9_966, 0]),
        &[24, 64, 128, 256],
        no_pair_skipped,
        |y, x, p, rm, cc| y.atan2(x, p, rm, cc),
        |r, y, x, rd| r.assign_round(y.atan2_ref(x), rd),
    )
}

/// `ln` of a complex number whose parts differ by a factor `2^(+-2 * 1442695040)` (`ln(1 + exp(-1e9) i)`,
/// `ln(-exp(-1e9) + exp(1e9) i)`): the imaginary part is `atan2` and must be the tiny angle (not `NaN`) and `pi/2`
/// (not `pi`); the real part stays finite.
#[test]
fn complex_ln_with_extreme_parts() -> Res {
    const K: i64 = 1_442_695_040;
    let mut run = Run::new(true)?;
    let mut cc = Consts::new()?;
    // (re, im) as (significand, exponent, negative)
    type Part = (u64, i64, bool);
    let cases: [(Part, Part); 3] = [
        ((1, 0, false), (0xdead_beef_1234_5677, -K, false)),
        ((0xdead_beef_1234_5677, -K, true), (1, K, false)),
        ((3, K, false), (0xdead_beef_1234_5677, -K, true)),
    ];
    for &((rm_, rk, rn), (im_, ik, inn)) in &cases {
        let z = ExactComplex::new(zf(rm_, rk, rn), zf(im_, ik, inn));
        let (fre, fim) = (mp(rm_, rk, rn), mp(im_, ik, inn));
        for &p in &[64usize, 128, 256] {
            for rm in MODES {
                let got = z.ln(p, rm, &mut cc);
                assert!(
                    !got.re().is_nan() && !got.im().is_nan(),
                    "ln({rm_:#x} 2^{rk} + {im_:#x} 2^{ik} i) is NaN"
                );
                assert!(!got.re().is_inf(), "ln(..) real part is infinite");
                let info = format!("ln im({rm_:#x} 2^{rk} neg {rn}, {im_:#x} 2^{ik} neg {inn})");
                run.check(got.im(), p, rm, &info, |r, rnd| {
                    ternary(r.assign_round(fim.atan2_ref(&fre), round_of(rnd)))
                });
                // real part: ln(max) + (1/2) ln(1 + (min/max)^2) at 1500 bits. Compared where its size is within
                // MPFR's exponent range (the first case's is about 2^-2.9e9: it underflows, and is only required to
                // be finite).
                if rk.max(ik) > 0 {
                    let (big, small) = if rk >= ik { (&fre, &fim) } else { (&fim, &fre) };
                    let hi = Float::with_val(1500, big.abs_ref()).ln()
                        + Float::with_val(
                            1500,
                            (Float::with_val(1500, small.abs_ref())
                                / Float::with_val(1500, big.abs_ref()))
                            .square()
                            .ln_1p(),
                        ) / 2;
                    let info = format!("ln re({rm_:#x} 2^{rk}, {im_:#x} 2^{ik} i)");
                    run.check(got.re(), p, rm, &info, |r, rnd| {
                        ternary(r.assign_round(&hi, round_of(rnd)))
                    });
                }
            }
        }
    }
    run.finish("complex ln", 20);
    Ok(())
}
