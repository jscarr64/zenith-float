//! The arithmetic primitives beyond `+ - * /` over the whole exponent range, bit for bit against MPFR (ported from the
//! option-a test tree and adapted to 1.0.x: precision rounded up to a word, overflow is `+-Inf` in every mode):
//! `sqrt`, `cbrt`, `reciprocal`, `hypot`, `fma`, `atan2` and `nth_root` at operands from below the subnormal
//! band to the top of the range, at precisions from 1 bit to 200 (and a few at 20 000), in the six
//! rounding modes MPFR can reproduce (`ToOdd` through the `p + 1` bit construction of `mpfr_reference_with`).
//! Subnormal and overflowing results are compared too (MPFR's range is emulated by `mpfr_reference_with`), and so is
//! the per-operation inexact flag. Every mismatch of a run is listed at once.

use gmp_mpfr_sys::mpfr::{self, rnd_t};
use rug::float::Round;
use rug::ops::AssignRound;
use rug::Float;
use std::cmp::Ordering;
use std::ops::Neg;
use std::panic::{catch_unwind, AssertUnwindSafe};
use zenith_float_num::{
    Consts, ExactNum, Exponent, RoundingMode, Sign, EXPONENT_MAX, EXPONENT_MIN, WORD_BIT_SIZE,
};

type Res = Result<(), Box<dyn std::error::Error>>;

// ---------------------------------------------------------------------------------------------
// Reference helpers (ported from the option-a test tree; main's mpfr/common.rs has no equivalent).
// ---------------------------------------------------------------------------------------------

fn round_of(rnd: rnd_t) -> Round {
    match rnd {
        rnd_t::RNDZ => Round::Zero,
        rnd_t::RNDU => Round::Up,
        rnd_t::RNDD => Round::Down,
        rnd_t::RNDA => Round::AwayZero,
        _ => Round::Nearest,
    }
}

fn ternary(o: Ordering) -> i32 {
    o as i32
}

/// The precision zenith-float 1.0.x computes at: the requested precision rounded up to a whole word.
fn p_eff(p: usize) -> usize {
    p.div_ceil(WORD_BIT_SIZE) * WORD_BIT_SIZE
}

/// MPFR's default exponent range is narrower than zenith-float's: widen it to the library's for this
/// thread so that operands (and expected values) can be built at any exponent.
fn use_full_exponent_range() {
    // SAFETY: only sets MPFR's thread-local exponent bounds, as the other tests of this suite do.
    unsafe {
        mpfr::set_emin(EXPONENT_MIN as mpfr::exp_t - 1024);
        mpfr::set_emax(EXPONENT_MAX as mpfr::exp_t);
    }
}

/// The value MPFR rounds to for `f` at `p` bits in zenith's format (`m * 2^e`, `0.5 <= m < 1`, subnormals
/// on the grid `2^(EXPONENT_MIN - p)`, emulated by `emin` and `mpfr_subnormalize`; `ToOdd` by rounding
/// toward zero at `p + 1` bits). `emax` is the largest exponent MPFR may use. The second value is the
/// inexact flag.
fn reference_in_range<F>(
    p: u32,
    rm: RoundingMode,
    emax: i64,
    f: &F,
) -> Result<(Float, bool), String>
where
    F: Fn(&mut Float, rnd_t) -> i32,
{
    // SAFETY: MPFR's exponent range is thread-local state; it is restored before returning.
    let (old_min, old_max) = unsafe { (mpfr::get_emin(), mpfr::get_emax()) };
    unsafe {
        mpfr::set_emax(emax as mpfr::exp_t);
        mpfr::set_emin(EXPONENT_MIN as mpfr::exp_t - p as mpfr::exp_t + 1);
    }
    let rnd = match rm {
        RoundingMode::ToEven => Some(rnd_t::RNDN),
        RoundingMode::Up => Some(rnd_t::RNDU),
        RoundingMode::Down => Some(rnd_t::RNDD),
        RoundingMode::ToZero => Some(rnd_t::RNDZ),
        RoundingMode::FromZero => Some(rnd_t::RNDA),
        RoundingMode::ToOdd => None,
        RoundingMode::None => {
            unsafe {
                mpfr::set_emin(old_min);
                mpfr::set_emax(old_max);
            }
            return Err(format!("no MPFR reference for {rm:?}"));
        }
    };
    let out = if let Some(rnd) = rnd {
        let mut r = Float::new(p);
        let t = f(&mut r, rnd);
        let t2 = unsafe { mpfr::subnormalize(r.as_raw_mut(), t, rnd) };
        (r, t2 != 0 || t != 0)
    } else {
        // `ToOdd` is round to nearest, ties to odd; MPFR has no such mode. In the subnormal range the result
        // has fewer than `p` bits: find its precision `q` from the exponent of the value, then round to `q` bits.
        // SAFETY: reads MPFR exponent bounds / the exponent of a number; the bounds are restored below.
        let q = unsafe {
            mpfr::set_emin(EXPONENT_MIN as mpfr::exp_t - 2 * p as mpfr::exp_t - 8);
            let mut probe = Float::new(p + 1);
            f(&mut probe, rnd_t::RNDZ);
            if probe.is_zero() || !probe.is_finite() {
                p as i64
            } else {
                p as i64 - (EXPONENT_MIN as i64 - mpfr::get_exp(probe.as_raw()) as i64).max(0)
            }
        };
        if q >= p as i64 {
            unsafe { mpfr::set_emin(EXPONENT_MIN as mpfr::exp_t - p as mpfr::exp_t + 1) };
            to_odd(p, p, &f)
        } else if q >= 1 {
            to_odd(p, q as u32, &f)
        } else {
            // below the grid: at or above half the smallest subnormal rounds to it (a tie goes to the odd one), below to zero
            let mut probe = Float::new(p + 1);
            f(&mut probe, rnd_t::RNDZ);
            let neg = probe.is_sign_negative();
            let mut v = Float::with_val(p, 0);
            if q == 0 {
                // SAFETY: a power of two inside the widest range; it only builds the smallest subnormal.
                unsafe {
                    mpfr::set_si_2exp(
                        v.as_raw_mut(),
                        if neg { -1 } else { 1 },
                        (EXPONENT_MIN as mpfr::exp_t) - p as mpfr::exp_t,
                        rnd_t::RNDN,
                    )
                };
            } else if neg {
                v = -v;
            }
            (v, true)
        }
    };
    unsafe {
        mpfr::set_emin(old_min);
        mpfr::set_emax(old_max);
    }
    Ok(out)
}

/// `f` rounded to nearest at `q` bits with ties to the odd neighbour (the caller has set MPFR's exponent range);
/// the value is returned as a `p`-bit number.
fn to_odd<F>(p: u32, q: u32, f: &F) -> (Float, bool)
where
    F: Fn(&mut Float, rnd_t) -> i32,
{
    let mut w = Float::new(q + 1);
    let t = f(&mut w, rnd_t::RNDZ);
    let mut lo = Float::new(q);
    let mut hi = Float::new(q);
    // SAFETY: plain MPFR roundings of a finite number into numbers of their own precision.
    unsafe {
        mpfr::set(lo.as_raw_mut(), w.as_raw(), rnd_t::RNDZ);
        mpfr::set(hi.as_raw_mut(), w.as_raw(), rnd_t::RNDA);
    }
    let odd = |x: &Float, bits: u32| {
        // SAFETY: `min_prec` only reads a number.
        !x.is_zero() && unsafe { mpfr::min_prec(x.as_raw()) } as u32 == bits
    };
    let v = if lo == hi {
        lo
    } else if t == 0 {
        if odd(&lo, q) {
            lo
        } else {
            hi
        }
    } else if odd(&w, q + 1) {
        hi
    } else {
        lo
    };
    // SAFETY: as above.
    let inexact = t != 0 || unsafe { mpfr::min_prec(w.as_raw()) } as u32 > q;
    (Float::with_val(p, v), inexact)
}

/// What zenith-float 1.0.x answers: the MPFR value, except that an exponent overflow is `+-Inf` in every
/// rounding mode (documented: `ExponentOverflow(Sign)` becomes `+-Inf`; `mul`, `div`, `sqrt`... all follow
/// it, MPFR's directed modes would return the largest finite number), and an exact zero is `+0`. A result is an overflow when
/// the correctly rounded value, with the exponent range unbounded above, reaches `2^EXPONENT_MAX`.
fn mpfr_reference_with<F>(p: u32, rm: RoundingMode, f: F) -> Result<(Float, bool), String>
where
    F: Fn(&mut Float, rnd_t) -> i32,
{
    let wide = reference_in_range(p, rm, EXPONENT_MAX as i64 + 4096, &f)?;
    let overflow = !wide.0.is_zero()
        && !wide.0.is_nan()
        && (wide.0.is_infinite()
            // SAFETY: reads the exponent of a finite non-zero MPFR number (rug's `get_exp` is an `i32`).
            || unsafe { mpfr::get_exp(wide.0.as_raw()) } as i64 > EXPONENT_MAX as i64);
    if overflow {
        let inf = if wide.0.is_sign_negative() {
            Float::with_val(p, infinity(true))
        } else {
            Float::with_val(p, infinity(false))
        };
        return Ok((inf, true));
    }
    let (v, inexact) = reference_in_range(p, rm, EXPONENT_MAX as i64, &f)?;
    if v.is_zero() && !inexact {
        // an exact zero (x - x) is `+0` in every rounding mode (MPFR: `-0` in `Down`); `add` does the same
        return Ok((Float::with_val(p, 0), false));
    }
    Ok((v, inexact))
}

/// An infinity of either sign.
fn infinity(negative: bool) -> rug::float::Special {
    if negative {
        rug::float::Special::NegInfinity
    } else {
        rug::float::Special::Infinity
    }
}

/// Bit-exact comparison of a zenith result with an MPFR value, including Inf, the sign of zero, subnormals.
fn assert_exact_eq(n: &ExactNum, f: &Float, op: &str) -> Res {
    if f.is_nan() {
        assert!(n.is_nan(), "{op}: MPFR NaN, zenith {n:?}");
        return Ok(());
    }
    if f.is_infinite() {
        assert!(n.is_inf(), "{op}: MPFR {f}, zenith {n:?}");
        assert_eq!(n.is_negative(), f.is_sign_negative(), "{op}: sign of Inf");
        return Ok(());
    }
    assert!(!n.is_inf() && !n.is_nan(), "{op}: MPFR {f}, zenith {n:?}");
    if f.is_zero() {
        assert!(n.is_zero(), "{op}: MPFR {f}, zenith {n:?}");
        assert_eq!(n.is_negative(), f.is_sign_negative(), "{op}: sign of zero");
        return Ok(());
    }
    // Exact value of the zenith number, rebuilt from its words (parsing would flush subnormals).
    let (words, _, sign, e, _) = n.as_raw_parts().ok_or("none value")?;
    let (old_min, old_max) = unsafe { (mpfr::get_emin(), mpfr::get_emax()) };
    unsafe {
        mpfr::set_emin(mpfr::get_emin_min());
        mpfr::set_emax(mpfr::get_emax_max());
    }
    let mut z = rug::Integer::new();
    for w in words.iter().rev() {
        z <<= 64;
        z += *w;
    }
    let bits = (64 * words.len()) as u32;
    let mut got = Float::with_val(bits.max(64), &z);
    got >>= (64 * words.len()) as isize - e as isize;
    if sign == Sign::Neg {
        got = -got;
    }
    let mut want = Float::with_val(bits.max(64), 0);
    unsafe { mpfr::set(want.as_raw_mut(), f.as_raw(), rnd_t::RNDN) };
    let equal = got == want;
    unsafe {
        mpfr::set_emin(old_min);
        mpfr::set_emax(old_max);
    }
    assert!(equal, "{op}: MPFR {f}, zenith {n:?} = {got}");
    Ok(())
}

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

/// `m * 2^k` exactly. `ldexp` flushes a result below the smallest exponent to zero in 1.0.x, so a subnormal
/// operand is built as the product of two exact factors (`mul` forms subnormals).
fn zf(m: u64, k: i64, negative: bool) -> ExactNum {
    let p = OPERAND_PRECISION;
    let x = if k >= MIN {
        ExactNum::from_word(m, p).ldexp(k as Exponent, p, RoundingMode::ToEven)
    } else {
        let a = ExactNum::from_word(m, p).ldexp(MIN as Exponent, p, RoundingMode::ToEven);
        let b = ExactNum::from_word(1, p).ldexp((k - MIN) as Exponent, p, RoundingMode::ToEven);
        a.mul(&b, p, RoundingMode::ToEven)
    };
    if negative {
        x.neg()
    } else {
        x
    }
}

fn split_exponent(k: i64) -> (Exponent, Exponent) {
    let hi = k.max(MIN);
    ((hi as Exponent), (k - hi) as Exponent)
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
    k + 64 - m.leading_zeros() as i64 > MIN - p_eff(p) as i64 + 1
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
        let flags = self.flags;
        let outcome = catch_unwind(AssertUnwindSafe(|| -> Res {
            let (want, inexact) = mpfr_reference_with(p_eff(p) as u32, rm, f)?;
            assert_exact_eq(got, &want, info)?;
            if flags && !want.is_nan() && !want.is_infinite() {
                assert_eq!(got.inexact(), inexact, "{info}: inexact flag, MPFR {want}");
            }
            Ok(())
        }));
        let text = match outcome {
            Ok(Ok(())) => return,
            Ok(Err(e)) => format!("{info} p={p} {rm:?}: error {e}"),
            Err(payload) => match (
                payload.downcast_ref::<String>(),
                payload.downcast_ref::<&str>(),
            ) {
                (Some(m), _) => format!("p={p} {rm:?} {m}"),
                (None, Some(m)) => format!("{info} p={p} {rm:?}: {m}"),
                (None, None) => format!("{info} p={p} {rm:?}: panic without a message"),
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

/// The overflow of a reciprocal is `Inf` in every rounding mode on 1.0.x, like every other operation (the tests
/// expect that, see `mpfr_reference_with`); its sign is the operand's.
fn reciprocal_run(skip: Skip, only_known_wrong: bool) -> Res {
    let ops: Vec<(u64, i64)> = unary_operands()
        .into_iter()
        .filter(|&(_, k)| !only_known_wrong || k <= MIN)
        .collect();
    let signs: &[bool] = if only_known_wrong { &[true] } else { &[false, true] };
    let precisions: &[usize] = if only_known_wrong { &[129, 200, 256] } else { &PRECISIONS };
    unary(
        "reciprocal",
        Sweep {
            operands: &ops,
            signs,
            precisions,
            flags: true,
            skip,
        },
        |x, p, rm| x.reciprocal(p, rm),
        |r, a, rd| r.assign_round(a.recip_ref(), rd),
    )
}

#[test]
fn reciprocal_over_the_whole_range() -> Res {
    reciprocal_run(no_skip, false)
}

/// Fixed in 1.0.13: the reciprocal of a negative number at or below `2^EXPONENT_MIN` overflows; at 192 bits or more
/// (the Newton path) the `Inf` was `+Inf` instead of `-Inf` (`div` gives `-Inf`).
#[test]
fn reciprocal_of_a_tiny_negative_number_keeps_its_sign() -> Res {
    reciprocal_run(no_skip, true)
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
    for n in [3usize, 5, 7] {
        let name = format!("nth_root n={n}");
        // the inexact flag of an exact result is checked by the known-wrong test below
        unary(
            &name,
            Sweep {
                operands: &unary_operands(),
                signs: &[false],
                precisions: &[1, 3, 24, 64, 65, 128],
                flags: false,
                skip: no_skip,
            },
            |x, p, rm| x.nth_root(n, p, rm),
            |r, a, rd| r.assign_round(a.root_ref(n as u32), rd),
        )?;
    }
    // exact: 3^40 = (3^8)^5, with the flag checked at the precisions below 64 (see the known-wrong test)
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
            flags: false,
            skip: no_skip,
        },
        |x, p, rm| x.nth_root(5, p, rm),
        |r, a, rd| r.assign_round(a.root_ref(5), rd),
    )
}

/// Fixed in 1.0.13: `nth_root` with a composite `n` rounded after every square or cube root step (`n = 4, 6, 8, 9,
/// 12`, and `n = 10, 15, 33` which continued with Newton from the rounded value): the result was one ulp off in about one
/// case in ten (`ToEven`) and one in five (directed modes). Values only (the flag is checked by the exact-result test).
#[test]
fn nth_root_with_a_composite_n() -> Res {
    let mut ops = unary_operands();
    ops.extend([(607, 57), (0xdead_beef_1234_5677, -1000), (u64::MAX, 7), (5, MIN + 3)]);
    for n in [4usize, 6, 8, 9, 10, 12, 15, 33] {
        let name = format!("nth_root n={n}");
        unary(
            &name,
            Sweep {
                operands: &ops,
                signs: &[false],
                precisions: &[1, 24, 64, 65, 128, 200],
                flags: false,
                skip: no_skip,
            },
            |x, p, rm| x.nth_root(n, p, rm),
            |r, a, rd| r.assign_round(a.root_ref(n as u32), rd),
        )?;
    }
    Ok(())
}

/// Fixed in 1.0.13: an exact `nth_root` (`2^7` for `n = 7`, also `2^2147483583`) at 64 bits and more is
/// was flagged inexact in `ToEven` and `ToOdd` (the value was right; the flag was not). A precision below 64 bits
/// computes at 64.
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

/// The exponents of the operands differ by more than `2^20`: the square of the ratio is far below the smallest
/// number (the 1.0.12 result was the larger operand, exact).
fn far_apart(a: (u64, i64), b: (u64, i64)) -> bool {
    (a.1 - b.1).abs() > 1 << 20
}

fn not_far_apart(a: (u64, i64), b: (u64, i64)) -> bool {
    !far_apart(a, b)
}

fn hypot_run(skip_pair: fn((u64, i64), (u64, i64)) -> bool) -> Res {
    binary(
        "hypot",
        &binary_operands(),
        &PRECISIONS,
        skip_pair,
        |a, b, p, rm, _| a.hypot(b, p, rm),
        |r, a, b, rd| r.assign_round(a.hypot_ref(b), rd),
    )
}

#[test]
fn hypot_over_the_whole_range() -> Res {
    hypot_run(|_, _| false)
}

/// Fixed in 1.0.13: `hypot` of operands whose exponents differ by more than `2^20` (for example
/// `hypot(1, 2^-2147483643)`, or with a subnormal operand) returns the larger one exactly (`inexact` false, and
/// `1` instead of `1 + ulp` in `Up`/`FromZero`): the square of the ratio underflows to zero.
#[test]
fn hypot_with_an_extreme_ratio() -> Res {
    hypot_run(not_far_apart)
}

/// Fixed in 1.0.13: `hypot(1, 2^-g)` for `g` from 2^6 to 2^24 returned `NaN(PrecisionRetryExhausted)` from about
/// `g = 9887` in every rounding mode, and `1` flagged exact below `2^-2^30`.
#[test]
fn hypot_with_a_large_ratio() -> Res {
    let mut run = Run::new(true)?;
    for g in [
        64i64,
        300,
        1000,
        9000,
        9887,
        9966,
        20_000,
        1 << 20,
        1 << 24,
        1_000_000_000,
        1_442_695_040,
    ] {
        let (big, small) = (zf(1, 0, false), zf(3, -g, true));
        let (fbig, fsmall) = (mp(1, 0, false), mp(3, -g, true));
        for p in [24usize, 64, 65, 128, 200] {
            for rm in MODES {
                for (x, y, fx, fy) in
                    [(&big, &small, &fbig, &fsmall), (&small, &big, &fsmall, &fbig)]
                {
                    let got = x.hypot(y, p, rm);
                    let info = format!("hypot(1, 3 2^-{g})");
                    run.check(&got, p, rm, &info, |r, rnd| {
                        ternary(r.assign_round(fx.hypot_ref(fy), round_of(rnd)))
                    });
                }
            }
        }
    }
    run.finish("hypot ratio", 100);
    Ok(())
}

/// Exponent classes of the `binary_operands` grid: `A` subnormal, `B` just above the smallest normal, `C` one,
/// `D` and `E` at the top.
const A: i64 = MIN - 70;
const B: i64 = MIN + 5;
const C: i64 = 0;
const D: i64 = MAX - 130;
const E: i64 = MAX - 64;

/// `(exponent of y, exponent of x)` of the operand pairs whose quotient `y / x` is far from one, so that the result
/// is `+-pi/2` (or `pi/2 + tiny`) or a tiny number (wrong, or `NaN`, in 1.0.12).
const ATAN2_EXTREME: [(i64, i64); 14] = [
    (D, B),
    (E, B),
    (C, A),
    (D, A),
    (E, A),
    (C, B),
    (A, C),
    (B, C),
    (C, D),
    (C, E),
    (B, D),
    (B, E),
    (A, D),
    (A, E),
];

fn atan2_is_not_extreme(a: (u64, i64), b: (u64, i64)) -> bool {
    !ATAN2_EXTREME.contains(&(a.1, b.1))
}

fn atan2_run(skip_pair: fn((u64, i64), (u64, i64)) -> bool) -> Res {
    binary(
        "atan2",
        &binary_operands(),
        &[1, 3, 24, 64, 128],
        skip_pair,
        |y, x, p, rm, cc| y.atan2(x, p, rm, cc),
        |r, y, x, rd| r.assign_round(y.atan2_ref(x), rd),
    )
}

#[test]
fn atan2_over_the_whole_range() -> Res {
    atan2_run(|_, _| false)
}

/// Fixed in 1.0.13: `atan2` was wrong, or `NaN(PrecisionRetryExhausted)`, when the quotient of the arguments is
/// extreme: `atan2(2^2147483517, -2^-2147483643)` was `pi` (rounded) instead of `pi/2`; `atan2(1, 2^2147483517)`
/// and `atan2(2^-2147483718, 1)` were `NaN` in the modes that need the retry loop.
#[test]
fn atan2_with_an_extreme_quotient() -> Res {
    atan2_run(atan2_is_not_extreme)
}

/// Fixed in 1.0.13: `atan2(2^-g, 1)`, `atan2(1, 2^g)` and the same with negative signs, `g` from 2^6 to 2^24, were
/// `NaN(PrecisionRetryExhausted)` in the directed modes and `ToOdd` from about `g = 9888` (the result is `y/x` minus
/// a sticky bit; `ToEven` and the flag were right).
#[test]
fn atan2_with_a_large_quotient() -> Res {
    let mut run = Run::new(true)?;
    // 9966 is `1e3000`, 1_442_695_040 is `exp(1e9)`: the quotients an expression evaluator reaches.
    for g in [
        64i64,
        300,
        1000,
        9000,
        9888,
        9966,
        20_000,
        1 << 20,
        1 << 24,
        1_000_000_000,
        1_442_695_040,
    ] {
        // (y as (significand, exponent), x as (significand, exponent)): a tiny quotient, a huge one
        // from a huge `y` or a tiny `x` (in range up to g = 2^30, beyond it `y / x` overflows)
        for ((my, ky), (mx, kx)) in [
            ((3u64, -g), (1u64, 0i64)),
            ((1, 0), (3, g)),
            ((1, g), (3, -g)),
            ((1, -g), (3, g)),
            ((3, g), (1, 0)),
            ((1, 0), (3, -g)),
        ] {
            for (neg_y, neg_x) in [(false, false), (true, false), (false, true), (true, true)] {
                let (y, x) = (zf(my, ky, neg_y), zf(mx, kx, neg_x));
                let (fy, fx) = (mp(my, ky, neg_y), mp(mx, kx, neg_x));
                for p in [24usize, 64, 65, 128, 200] {
                    for rm in MODES {
                        let got = y.atan2(&x, p, rm, &mut run.cc);
                        let info =
                            format!("atan2({my} 2^{ky} neg {neg_y}, {mx} 2^{kx} neg {neg_x})");
                        run.check(&got, p, rm, &info, |r, rnd| {
                            ternary(r.assign_round(fy.atan2_ref(&fx), round_of(rnd)))
                        });
                    }
                }
            }
        }
    }
    run.finish("atan2 quotient", 100);
    Ok(())
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
        let c0 = a.mul(&b, OPERAND_PRECISION, RoundingMode::ToEven).neg();
        let c0_shift = c0.mul(&zf(1, -60, false), OPERAND_PRECISION, RoundingMode::ToEven);
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
    // total and near cancellation at every pair; an independent addend where it is not far below the product
    fma_run(
        &["-ab", "-ab(1+2^-60)"],
        &[0, 1, 2, 3, 4, 5, 6, 7],
        &FMA_PRECISIONS,
        500,
    )?;
    fma_run(&["indep"], &[0, 1, 2, 3, 4, 5, 6, 7], &FMA_PRECISIONS, 50)
}

/// Fixed in 1.0.13: `fma` whose addend is far below (or far above) the product, and in 1.0.12 more than
/// `p + 128` bits away, was `c` or `a * b` rounded alone with the inexact flag set: the directed modes
/// (`fma(1, 1, 2^-300)` in `Up` was `1`) and, when the other term is a power of two minus the tiny one
/// (`fma(1, -1, 2^-200)`), the nearest modes' flag were wrong. Also with the addend far above a tiny product, and
/// at the ends of the exponent range (the independent addends below).
#[test]
fn fma_with_a_far_smaller_addend() -> Res {
    let mut run = Run::new(true)?;
    let one = zf(1, 0, false);
    let minus_one = zf(1, 0, true);
    let (f_one, f_minus_one) = (mp(1, 0, false), mp(1, 0, true));
    // the product is exact (a one-bit times a one-bit number); the sum is rounded once, by `assign_round`
    let product = Float::with_val(256, &f_one * &f_minus_one);
    for gap in [100i64, 200, 1000, 100_000] {
        let c = zf(1, -gap, false);
        let fc = mp(1, -gap, false);
        for p in [24usize, 64, 65, 128, 200] {
            for rm in MODES {
                let got = one.fma(&minus_one, &c, p, rm);
                let info = format!("fma(1, -1, 2^-{gap})");
                run.check(&got, p, rm, &info, |r, rnd| {
                    ternary(r.assign_round(&product + &fc, round_of(rnd)))
                });
            }
        }
    }
    // the addend far ABOVE a tiny product of the opposite sign: fma(2^-gap, 1, -1)
    for gap in [100i64, 200, 1000, 100_000] {
        let a = zf(1, -gap, false);
        let fa = mp(1, -gap, false);
        let product = Float::with_val(256, &fa * &f_one);
        for p in [24usize, 64, 65, 128, 200] {
            for rm in MODES {
                let got = a.fma(&one, &minus_one, p, rm);
                let info = format!("fma(2^-{gap}, 1, -1)");
                run.check(&got, p, rm, &info, |r, rnd| {
                    ternary(r.assign_round(&product + &f_minus_one, round_of(rnd)))
                });
            }
        }
    }
    run.finish("fma(1, -1, tiny)", 20);
    fma_run(&["indep"], &[0, 1, 2, 4, 5, 6], &FMA_PRECISIONS, 50)
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
