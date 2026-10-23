// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! Components used in MPFR integration tests

use gmp_mpfr_sys::mpfr::{self, rnd_t};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rug::Float;
use std::cell::RefCell;
use zenith_float_num::{
    Consts, ExactNum, Exponent, Radix, RoundingMode, Sign, Word, EXPONENT_MAX, EXPONENT_MIN,
    WORD_BIT_SIZE, WORD_SIGNIFICANT_BIT,
};

/// Which MPFR tier runs: `MPFR_TIER=full` keeps every run count and runs the 100032-bit tests;
/// anything else (the default) is the fast tier, which divides the heavy run counts by 10 and skips
/// the 100032-bit tests. `tools/full_gate.sh` runs the full tier unless `MPFR_TIER` says otherwise.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum MpfrTier {
    Fast,
    Full,
}

pub fn mpfr_tier() -> MpfrTier {
    match std::env::var("MPFR_TIER") {
        Ok(t) if t == "full" => MpfrTier::Full,
        _ => MpfrTier::Fast,
    }
}

/// `full` in the full tier; `full / 10` (at least 1) in the fast tier.
pub fn tier_runs(full: usize) -> usize {
    match mpfr_tier() {
        MpfrTier::Full => full,
        MpfrTier::Fast => (full / 10).max(1),
    }
}

/// True (and a note on stderr) when the fast tier skips a 100032-bit test.
pub fn fast_tier_skips_100032_bits(name: &str) -> bool {
    let skip = mpfr_tier() == MpfrTier::Fast;
    if skip {
        eprintln!("{name}: skipped in the fast MPFR tier (set MPFR_TIER=full to run it)");
    }
    skip
}

thread_local! {
    static TEST_RNG: RefCell<StdRng> = RefCell::new(StdRng::seed_from_u64(0x5EED_CAFE_BADC_0D00));
}

pub fn test_random<T>() -> T
where
    rand::distributions::Standard: rand::distributions::Distribution<T>,
{
    TEST_RNG.with(|rng| rng.borrow_mut().gen())
}

pub fn reset_test_rng() {
    TEST_RNG.with(|rng| {
        *rng.borrow_mut() = StdRng::seed_from_u64(0x5EED_CAFE_BADC_0D00);
    });
    // per thread: a reseed from another test thread must not move this test's operands
    zenith_float_num::reseed_random_this_thread(zenith_float_num::DEFAULT_RANDOM_SEED);
}

/// Draws a seed from the test stream; `rebuild_from_seed` replays a generator from it.
pub fn draw_seed() -> u64 {
    test_random::<u64>()
}

/// Runs `build` with both test streams (this file's and the crate's) restarted from `seed`, so it yields the
/// same operand every time, then resumes both streams from a fresh seed drawn beforehand: a rebuild neither
/// repeats nor shifts the draws that follow.
pub fn rebuild_from_seed<T>(seed: u64, build: impl FnOnce() -> T) -> T {
    let resume = test_random::<u64>();
    TEST_RNG.with(|rng| *rng.borrow_mut() = StdRng::seed_from_u64(seed));
    zenith_float_num::reseed_random_this_thread(seed);
    let built = build();
    TEST_RNG.with(|rng| *rng.borrow_mut() = StdRng::seed_from_u64(resume));
    zenith_float_num::reseed_random_this_thread(resume);
    built
}

/// The value `ExactNum::random_normal(p, emin, emax)` produces from `seed`.
pub fn normal_from_seed(seed: u64, p: usize, emin: Exponent, emax: Exponent) -> ExactNum {
    rebuild_from_seed(seed, || ExactNum::random_normal(p, emin, emax))
}

macro_rules! test_zf_op {
    ($eq:literal, $n1:ident, $n2:ident, $zf_op:ident, $f1:ident, $f2:ident, $mpfr_op:ident, $p:ident, $rm:ident, $rnd:ident, $op_info:expr, $cc:ident) => {
        let n3 = ExactNum::$zf_op(&($n1), &($n2), $p, $rm, &mut $cc);

        let mut f3 = Float::with_val($p as u32, 1);

        unsafe { mpfr::$mpfr_op(f3.as_raw_mut(), ($f1).as_raw(), ($f2).as_raw(), $rnd) };

        //println!("\n{:b}\n{:b}", $n1, $n2);
        //println!("\n{:b}\n{}", $n1, $f1.to_string_radix(2, None));
        //println!("\n{:b}\n{}", n3, f3.to_string_radix(2, None));

        assert_float_close(n3, f3, $p, &format!("{:?}", $op_info), $eq, &mut $cc)?;
    };
    ($eq:literal, $n1:ident, $zf_op:ident, $f1:ident, $mpfr_op:ident, $p:ident, $rm:ident, $rnd:ident, $op_info:expr, $cc:ident) => {
        let n3 = ExactNum::$zf_op(&($n1), $p, $rm, &mut $cc);

        let mut f3 = Float::with_val($p as u32, 1);

        unsafe { mpfr::$mpfr_op(f3.as_raw_mut(), ($f1).as_raw(), $rnd) };

        // println!("\n{:b}\n{}", $n1, $f1.to_string_radix(2, None));
        // println!("\n{:b}\n{}", n3, f3.to_string_radix(2, None));

        assert_float_close(n3, f3, $p, &format!("{:?}", $op_info), $eq, &mut $cc)?;
    };
    ($eq:literal, $n1:ident, $zf_op:ident, $f1:ident, $mpfr_op:ident, $p:ident, $rm:ident, $rnd:ident, $op_info:expr, $cc:ident) => {
        let n3 = ExactNum::$zf_op(&($n1), $p, $rm, &mut $cc);

        let mut f3 = Float::with_val($p as u32, 1);

        unsafe { mpfr::$mpfr_op(f3.as_raw_mut(), ($f1).as_raw(), $rnd) };

        // println!("\n{:b}\n{}", $n1, $f1.to_string_radix(2, None));
        // println!("\n{:b}\n{}", n3, f3.to_string_radix(2, None));

        assert_float_close(n3, f3, $p, &format!("{:?}", $op_info), $eq, &mut $cc)?;
    };
}

// `$n1` is an expression yielding the operand to consume (an `ExactNum` the call may move, or a place
// that the op only borrows); callers that need the operand again rebuild it from its seed.
macro_rules! test_zf_op_no_cc {
    ($eq:literal, $n1:expr, $n2:ident, $zf_op:ident, $f1:ident, $f2:ident, $mpfr_op:ident, $p:ident, $rm:ident, $rnd:ident, $op_info:expr, $cc:ident) => {
        let n3 = ($n1).$zf_op(&($n2), $p, $rm);

        let mut f3 = Float::with_val($p as u32, 1);

        unsafe { mpfr::$mpfr_op(f3.as_raw_mut(), ($f1).as_raw(), ($f2).as_raw(), $rnd) };

        //println!("\n{:b}\n{:b}", $n1, $n2);
        //println!("\n{:b}\n{}", $n1, $f1.to_string_radix(2, None));
        //println!("\n{:b}\n{}", n3, f3.to_string_radix(2, None));

        assert_float_close(n3, f3, $p, &format!("{:?}", $op_info), $eq, &mut $cc)?;
    };
    ($eq:literal, $n1:ident, $zf_op:ident, $f1:ident, $mpfr_op:ident, $p:ident, $rm:ident, $rnd:ident, $op_info:expr, $cc:ident) => {
        let n3 = ExactNum::$zf_op(&($n1), $p, $rm);

        let mut f3 = Float::with_val($p as u32, 1);

        unsafe { mpfr::$mpfr_op(f3.as_raw_mut(), ($f1).as_raw(), $rnd) };

        // println!("\n{:b}\n{}", $n1, $f1.to_string_radix(2, None));
        // println!("\n{:b}\n{}", n3, f3.to_string_radix(2, None));

        assert_float_close(n3, f3, $p, &format!("{:?}", $op_info), $eq, &mut $cc)?;
    };
}

macro_rules! test_zf_fma {
    ($eq:literal, $n1:ident, $n2:ident, $nc:ident, $f1:ident, $f2:ident, $fc:ident, $p:ident, $rm:ident, $rnd:ident, $op_info:expr, $cc:ident) => {
        let n4 = ExactNum::fma(&($n1), &($n2), &($nc), $p, $rm);
        let mut f4 = Float::with_val($p as u32, 1);
        unsafe {
            mpfr::fma(
                f4.as_raw_mut(),
                ($f1).as_raw(),
                ($f2).as_raw(),
                ($fc).as_raw(),
                $rnd,
            )
        };
        assert_float_close(n4, f4, $p, &format!("{:?}", $op_info), $eq, &mut $cc)?;
    };
}

macro_rules! test_zf_rem_pi {
    ($n1:ident, $f1:ident, $p:ident, $rm:ident, $rnd:ident, $op_info:expr, $cc:ident) => {
        let mut n3 = ExactNum::rem_pi(&($n1), $p, $rm, &mut $cc);
        let _ = n3.set_precision($p, $rm);
        let _ = $rnd;

        if ($n1).exponent().unwrap_or(0) <= 2 {
            // Identity path: rem_pi leaves |x| < 4 unchanged aside from working precision.
            assert_float_close(
                n3,
                ($f1).clone(),
                $p,
                &format!("{:?}", $op_info),
                true,
                &mut $cc,
            )?;
        } else if !n3.is_nan() && !n3.is_inf() && !n3.is_zero() {
            // Not IEEE fmod; reduction is validated by sin/cos oracles on the unreduced argument.
            let e = n3.exponent().unwrap_or(0);
            assert!(
                e <= 3,
                "{}",
                format!("{:?} rem_pi exponent {e} (want <= 3)", $op_info)
            );
        }
    };
}

// test constant value match
macro_rules! test_zf_const {
    ($zf_const:ident, $mpfr_const:ident, $p:ident, $rm:ident, $rnd:ident, $op_info:expr, $cc:ident) => {
        let n1: ExactNum = $cc.$zf_const($p, $rm);

        let mut f1 = Float::with_val($p as u32, 1);

        unsafe {
            mpfr::$mpfr_const(f1.as_raw_mut(), $rnd);
        }

        assert_float_close(n1, f1, $p, &format!("{:?}", $op_info), true, &mut $cc)?;
    };
}

pub(crate) use test_zf_const;
pub(crate) use test_zf_fma;
pub(crate) use test_zf_op;
pub(crate) use test_zf_op_no_cc;
pub(crate) use test_zf_rem_pi;

/// Exponent band the random MPFR comparisons draw from: the former 32-bit range. The 62-bit
/// edges are covered by `i64_edges`; drawing from them here would make the decimal `Debug`
/// messages and the exp/pow arguments exabyte-sized.
pub const RAND_EMIN: zenith_float_num::Exponent = i32::MIN as zenith_float_num::Exponent;
/// Upper end of [`RAND_EMIN`]'s band.
pub const RAND_EMAX: zenith_float_num::Exponent = i32::MAX as zenith_float_num::Exponent;

pub const fn get_prec_rng() -> usize {
    157
}

pub fn get_float_pair(
    p: usize,
    emin: Exponent,
    emax: Exponent,
    cc: &mut Consts,
) -> Result<(ExactNum, Float), Box<dyn std::error::Error>> {
    let n = ExactNum::random_normal(p, emin, emax);
    let f = conv_to_mpfr(p, &n, cc)?;
    Ok((n, f))
}

pub fn conv_to_mpfr(
    p: usize,
    n: &ExactNum,
    cc: &mut Consts,
) -> Result<Float, Box<dyn std::error::Error>> {
    // Past the former 32-bit range the binary strings get too large for MPFR's string
    // conversion; build the value from the raw words and the 62-bit exponent instead.
    if let Some((w, _, s, e, _)) = n.as_raw_parts() {
        if e < RAND_EMIN || e > RAND_EMAX {
            let bits = w.len() * WORD_BIT_SIZE;
            let int = rug::Integer::from_digits(w, rug::integer::Order::Lsf);
            let mut f = Float::with_val(u32::try_from(p)?, &int);
            f <<= (e - i64::try_from(bits)?) as isize;
            if s == Sign::Neg {
                f = -f;
            }
            return Ok(f);
        }
    }
    let s1 = conv_str_to_mpfr_compat(n.format(Radix::Bin, RoundingMode::ToEven, cc)?)?;
    let f = Float::with_val(p as u32, Float::parse_radix(s1, 2)?);
    let s2 = conv_str_from_mpfr_compat(f.to_string_radix(2, None))?;
    //println!("\n{}\n{}", s1, s2);
    assert_eq!(
        *n,
        ExactNum::parse(&s2, Radix::Bin, p, RoundingMode::None, cc)
    );
    Ok(f)
}

/// `Some(MPFR exponent)` when either value lies outside [`RAND_EMIN`, `RAND_EMAX`].
fn far_exponent(n: &ExactNum, f: &Float) -> Option<i64> {
    let e_f = if f.is_normal() { unsafe { (*f.as_raw()).exp } } else { 0 };
    let e_n = n.exponent().unwrap_or(0);
    let far = |e: i64| !(RAND_EMIN..=RAND_EMAX).contains(&e);
    (far(e_f) || far(e_n)).then_some(e_f)
}

/// Runs `f` with MPFR's exponent range widened to its maximum, then restores it.
pub fn with_full_range<T>(f: impl FnOnce() -> T) -> T {
    let (old_min, old_max) = unsafe { (mpfr::get_emin(), mpfr::get_emax()) };
    unsafe {
        mpfr::set_emin(mpfr::get_emin_min());
        mpfr::set_emax(mpfr::get_emax_max());
    }
    let r = f();
    unsafe {
        mpfr::set_emin(old_min);
        mpfr::set_emax(old_max);
    }
    r
}

// assert float values are equal or differ not more than 1 ulp.
pub fn assert_float_close(
    n: ExactNum,
    f: Float,
    p: usize,
    op: &str,
    eq: bool,
    cc: &mut Consts,
) -> Result<(), Box<dyn std::error::Error>> {
    if n.is_inf() {
        // Inf only where MPFR is Inf: an overflow rounded toward zero is the largest
        // finite value in MPFR, so the overflow flag alone does not admit Inf.
        assert!(f.is_infinite(), "{} zenith Inf, MPFR {}", op, show_mp(&f));
        if n.is_positive() {
            assert!(f.is_sign_positive(), "{}", op);
        } else {
            assert!(f.is_sign_negative(), "{}", op);
        }
    } else if f.is_infinite() {
        return Err(format!("{} MPFR Inf, zenith {:?}", op, n).into());
    } else if n.is_nan() {
        // nan
        if eq {
            assert!(f.is_nan(), "{}", op);
        }
    } else if (n.is_subnormal() || n.is_zero()) && eq {
        // subnormal (bit-exact mode only — otherwise compare below)
        let unf = unsafe { mpfr::underflow_p() };
        if !f.is_zero() {
            assert!(
                unf != 0,
                "{} zenith {:?}, MPFR exp {} prec {}",
                op,
                n.as_raw_parts(),
                if f.is_normal() { unsafe { (*f.as_raw()).exp } } else { 0 },
                f.prec()
            );
        }
    } else if let Some(e_f) = far_exponent(&n, &f) {
        // Past the former 32-bit range: MPFR's binary string would be enormous, so compare in
        // MPFR space with exact arithmetic over the full 62-bit range.
        let mlen = n.mantissa_max_bit_len().ok_or("none value")?;
        let e1 = n.exponent().ok_or("none value")?
            - i64::try_from(mlen - n.precision().ok_or("none value")?)?;
        let ok = with_full_range(|| -> Result<bool, Box<dyn std::error::Error>> {
            let g = conv_to_mpfr(mlen, &n, cc)?;
            if eq {
                return Ok(g == f);
            }
            let d = Float::with_val(u32::try_from(mlen + f.prec() as usize + 2)?, &g - &f).abs();
            if d.is_zero() {
                return Ok(true);
            }
            let e2 = unsafe { (*d.as_raw()).exp };
            Ok(e1 - e2 >= i64::try_from(mlen)? - 1)
        })?;
        assert!(ok, "{} far-range mismatch: MPFR exp {}, zenith {:?}", op, e_f, n.as_raw_parts());
    } else if eq {
        // n == f
        let s1 = f.to_string_radix(2, None);
        let s2 = conv_str_from_mpfr_compat(s1)?;
        let n2 = ExactNum::parse(&s2, Radix::Bin, p, RoundingMode::None, cc);
        assert_eq!(n, n2, "{}", op);
    } else {
        // at most 1 ulp difference
        let s1 = f.to_string_radix(2, None);
        let s2 = conv_str_from_mpfr_compat(s1)?;
        let n2 = ExactNum::parse(&s2, Radix::Bin, p, RoundingMode::None, cc);
        assert_eq!(
            n.mantissa_max_bit_len(),
            n2.mantissa_max_bit_len(),
            "{}",
            op
        );

        let d = ExactNum::sub_ref(&n, &n2, 1, RoundingMode::None).abs();
        assert!(
            d.exponent().is_some(),
            "{} zenith {:?}, MPFR {}, difference {:?}",
            op,
            n,
            f,
            d
        );

        let e1 = n.exponent().ok_or("none value")? as isize
            - (n.mantissa_max_bit_len().ok_or("none value")? - n.precision().ok_or("none value")?)
                as isize;
        let e2 = d.exponent().ok_or("none value")? as isize
            - (d.mantissa_max_bit_len().ok_or("none value")? - d.precision().ok_or("none value")?)
                as isize;

        //println!("\n{:?}\n{:?}\n{:?}\n{} {}", n, n2, d, e1, e2);

        if !d.is_zero() {
            assert!(
                e1 - e2 >= n.mantissa_max_bit_len().ok_or("none value")? as isize - 1,
                "{}",
                op
            );
        }
    }
    Ok(())
}

pub fn conv_str_to_mpfr_compat(s: String) -> Result<String, Box<dyn std::error::Error>> {
    let (sig, exp) = if let Some(pos) = s.find('e') {
        s.split_at(pos + 1)
    } else {
        (s.as_str(), "0")
    };
    let expn = i64::from_str_radix(exp, 2)?;
    Ok(sig.to_owned() + &expn.to_string())
}

pub fn conv_str_from_mpfr_compat(s: String) -> Result<String, Box<dyn std::error::Error>> {
    if let Some(epos) = s.find('e') {
        let (sig, exp) = s.split_at(epos + 1);
        let expn = exp.parse::<i64>()?;
        if expn < 0 {
            Ok(sig.to_owned() + "-" + &format!("{:b}", -expn))
        } else {
            Ok(sig.to_owned() + &format!("{:b}", expn))
        }
    } else {
        Ok(s)
    }
}

pub fn get_random_rnd_pair() -> (RoundingMode, rnd_t) {
    match test_random::<u8>() % 5 {
        0 => (RoundingMode::ToEven, rnd_t::RNDN),
        1 => (RoundingMode::Up, rnd_t::RNDU),
        2 => (RoundingMode::Down, rnd_t::RNDD),
        3 => (RoundingMode::FromZero, rnd_t::RNDA),
        _ => (RoundingMode::ToZero, rnd_t::RNDZ),
    }
}

// Generates a number with mantissa like 1111111..1110000..00000
pub fn get_oned_zeroed(p: usize, exp_from: Exponent, exp_to: Exponent) -> ExactNum {
    let zero_bits = test_random::<usize>() % (p - 1);
    let mut m1 = vec![Word::MAX; p / WORD_BIT_SIZE];
    let i = zero_bits / WORD_BIT_SIZE;
    m1.iter_mut().take(i).for_each(|v| *v = 0);
    m1[i] <<= zero_bits % WORD_BIT_SIZE;

    bf_from_mantissa_and_exp_rng(&m1, exp_from, exp_to)
}

// Generates a number with mantissa like 111..111000..000111..111
pub fn get_oned_sides(p: usize, exp_from: Exponent, exp_to: Exponent) -> ExactNum {
    let one_bits = test_random::<usize>() % (p / 2).min(WORD_BIT_SIZE * 2) + 1;
    let mut m1 = vec![0; p / WORD_BIT_SIZE];
    let i = one_bits / WORD_BIT_SIZE;
    let m1l = m1.len();
    m1.iter_mut().take(i + 1).for_each(|v| *v = Word::MAX);
    m1[i] >>= one_bits % WORD_BIT_SIZE;
    m1.iter_mut().rev().take(i + 1).for_each(|v| *v = Word::MAX);
    m1[m1l - i - 1] <<= one_bits % WORD_BIT_SIZE;

    bf_from_mantissa_and_exp_rng(&m1, exp_from, exp_to)
}

// Generates a number with periodic mantissa
pub fn get_periodic(p: usize, exp_from: Exponent, exp_to: Exponent) -> ExactNum {
    let nbits = 4;
    let bits = test_random::<Word>() % ((1 << nbits) - 1) + 1;
    let mut w: Word = bits;
    for _ in 1..WORD_BIT_SIZE / nbits {
        w <<= nbits;
        w |= bits;
    }
    let mut m1 = vec![w; p / WORD_BIT_SIZE];
    if let Some(last) = m1.last_mut() {
        *last |= 1 << (WORD_BIT_SIZE - 1);
    }

    bf_from_mantissa_and_exp_rng(&m1, exp_from, exp_to)
}

// Generates a number with oned mantissa, but last bit is zero
pub fn get_last_zero(p: usize, exp_from: Exponent, exp_to: Exponent) -> ExactNum {
    let mut m1 = vec![Word::MAX; p / WORD_BIT_SIZE];
    m1[0] = Word::MAX - 1;

    bf_from_mantissa_and_exp_rng(&m1, exp_from, exp_to)
}

// Generates a number near 1.
pub fn get_near_one(p: usize) -> ExactNum {
    let e = (test_random::<u8>() & 1) as Exponent;

    let random_bits = test_random::<usize>() % (p - 1);
    let i = random_bits / WORD_BIT_SIZE;

    let mut m1;
    if e == 0 {
        m1 = vec![Word::MAX; p / WORD_BIT_SIZE];
    } else {
        m1 = vec![0; p / WORD_BIT_SIZE];
    }
    m1[i] ^= test_random::<Word>() >> (random_bits % WORD_BIT_SIZE);

    m1.iter_mut().take(i).for_each(|v| *v = test_random());

    if let Some(last) = m1.last_mut() {
        *last |= WORD_SIGNIFICANT_BIT;
    }

    bf_from_mantissa_and_exp_rng(&m1, e, e)
}

pub fn bf_from_mantissa_and_exp_rng(m: &[Word], exp_from: Exponent, exp_to: Exponent) -> ExactNum {
    let e = if exp_from < exp_to {
        (test_random::<isize>().abs() % (exp_to as isize - exp_from as isize) + exp_from as isize)
            as Exponent
    } else {
        exp_from
    };

    let s = if test_random::<u8>() & 1 == 0 { Sign::Pos } else { Sign::Neg };

    ExactNum::from_words(m, s, e)
}

/// Binary operations with an exact-precision MPFR reference.
#[derive(Copy, Clone, Debug)]
pub enum BinOp {
    Add,
    Sub,
    Mul,
    Div,
}

/// Zenith result of `n1 op n2` where both operands stay live afterwards (the fuzz loops read them again
/// for the failure message and for the other ops/modes): the borrowed-borrowed forms, no copy.
pub fn zenith_binop(
    n1: &ExactNum,
    n2: &ExactNum,
    op: BinOp,
    p: usize,
    rm: RoundingMode,
) -> ExactNum {
    match op {
        BinOp::Add => ExactNum::add_ref(n1, n2, p, rm),
        BinOp::Sub => ExactNum::sub_ref(n1, n2, p, rm),
        BinOp::Mul => ExactNum::mul_ref(n1, n2, p, rm),
        BinOp::Div => ExactNum::div_ref(n1, n2, p, rm),
    }
}

fn mpfr_binop(r: &mut Float, a: &Float, b: &Float, op: BinOp, rnd: rnd_t) -> i32 {
    match op {
        BinOp::Add => unsafe { mpfr::add(r.as_raw_mut(), a.as_raw(), b.as_raw(), rnd) },
        BinOp::Sub => unsafe { mpfr::sub(r.as_raw_mut(), a.as_raw(), b.as_raw(), rnd) },
        BinOp::Mul => unsafe { mpfr::mul(r.as_raw_mut(), a.as_raw(), b.as_raw(), rnd) },
        BinOp::Div => unsafe { mpfr::div(r.as_raw_mut(), a.as_raw(), b.as_raw(), rnd) },
    }
}

/// `a op b` correctly rounded at exactly `p` bits in the zenith format: `m * 2^e` with
/// `0.5 <= m < 1`, `e` in `EXPONENT_MIN..=EXPONENT_MAX`, and subnormals on the fixed grid
/// `2^(EXPONENT_MIN - p)`. MPFR emulates that grid with `emin = EXPONENT_MIN - p + 1` and
/// `mpfr_subnormalize`. MPFR has no ties-to-odd mode: `ToOdd` rounds toward zero at `p + 1`
/// bits and decides from the extra bit and the ternary value. That emulation is for results
/// in the normal range only.
pub fn mpfr_reference(
    a: &Float,
    b: &Float,
    op: BinOp,
    p: u32,
    rm: RoundingMode,
) -> Result<Float, Box<dyn std::error::Error>> {
    Ok(mpfr_reference_with(p, rm, |r, rnd| mpfr_binop(r, a, b, op, rnd))?.0)
}

/// Like [`mpfr_reference`] for any MPFR function: `f(result, rnd)` computes into `result`
/// (whose precision is already set) and returns the MPFR ternary value. Returns the
/// result and whether it is inexact.
pub fn mpfr_reference_with<F>(
    p: u32,
    rm: RoundingMode,
    f: F,
) -> Result<(Float, bool), Box<dyn std::error::Error>>
where
    F: Fn(&mut Float, rnd_t) -> i32,
{
    mpfr_reference_shifted(p, rm, 0, f)
}

/// Like [`mpfr_reference_with`] in a space scaled by `2^shift`: the operands `f` reads are the zenith
/// operands times `2^shift`, and the zenith exponent range `EXPONENT_MIN..=EXPONENT_MAX` (with its
/// subnormal grid) is emulated as `EXPONENT_MIN + shift ..= EXPONENT_MAX + shift`. The result is in the
/// scaled space too. A positive shift lets MPFR (whose `emin_min` is `EXPONENT_MIN` on LP64 once the
/// exponent is 62 bits) hold the subnormal grid below the zenith minimum. The emulated range is clamped
/// to MPFR's; [`emulated_range_ok`] says whether the clamp was needed.
pub fn mpfr_reference_shifted<F>(
    p: u32,
    rm: RoundingMode,
    shift: i64,
    f: F,
) -> Result<(Float, bool), Box<dyn std::error::Error>>
where
    F: Fn(&mut Float, rnd_t) -> i32,
{
    let (old_min, old_max) = unsafe { (mpfr::get_emin(), mpfr::get_emax()) };
    let (emin, emax) = emulated_range(p, shift);
    unsafe {
        mpfr::set_emax(emax as mpfr::exp_t);
        mpfr::set_emin(emin as mpfr::exp_t);
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
            return Err(format!("no MPFR reference for {rm:?}").into());
        }
    };
    let out = if let Some(rnd) = rnd {
        let mut r = Float::new(p);
        let t = f(&mut r, rnd);
        let t2 = unsafe { mpfr::subnormalize(r.as_raw_mut(), t, rnd) };
        let inexact = t2 != 0 || t != 0;
        (r, inexact)
    } else {
        // Two binades below the smallest subnormal so that results in `[g / 2, g)` survive the truncation.
        unsafe {
            let lo = widen(mpfr::get_emin_min());
            mpfr::set_emin((emin - 2).max(lo) as mpfr::exp_t);
        }
        let mut w = Float::new(p + 1);
        let t = f(&mut w, rnd_t::RNDZ);
        unsafe {
            mpfr::set_emin(emin as mpfr::exp_t);
        }
        // A result under the library's smallest normal exponent has fewer than `p` bits (the subnormal grid is
        // `2^(EXPONENT_MIN - p)`): the nearest-ties-to-odd rounding is done at that many bits, from the
        // `p + 1` bit truncation (a truncation of a truncation is the truncation, and bits lost are a sticky bit).
        // (rug's `get_exp` is an `i32`, which cannot hold an exponent below the library's minimum)
        let exp = if w.is_zero() || !w.is_finite() {
            None
        } else {
            Some(unsafe { mpfr::get_exp(w.as_raw()) })
        };
        // At or below `2^(grid exponent)` no mantissa bit is left: the neighbours are 0 and the smallest
        // subnormal `g`, the midpoint `g / 2` is representable in `w`, so `|w| >= g / 2` (with the tie going
        // to the odd neighbour `g`) picks `g`. The result is always inexact.
        let grid = EXPONENT_MIN as i64 + shift - p as i64;
        if let Some(e) = exp {
            if widen(e) <= grid {
                let neg = w.is_sign_negative();
                // Build 1 (exponent 1) under the widest range: the emulated `emin` can be far above it.
                unsafe {
                    mpfr::set_emin(mpfr::get_emin_min());
                }
                let mut half = Float::with_val(64, 1);
                let mut g = Float::with_val(p, 1);
                unsafe {
                    mpfr::mul_2si(
                        half.as_raw_mut(),
                        half.as_raw(),
                        (grid - 1) as _,
                        rnd_t::RNDN,
                    );
                    mpfr::mul_2si(g.as_raw_mut(), g.as_raw(), grid as _, rnd_t::RNDN);
                }
                let r = if Float::with_val(p + 1, w.abs_ref()) >= half {
                    g
                } else {
                    Float::with_val(p, 0)
                };
                let r = if neg { -r } else { r };
                unsafe {
                    mpfr::set_emin(old_min);
                    mpfr::set_emax(old_max);
                }
                return Ok((r, true));
            }
        }
        let bits = match exp {
            Some(e) if widen(e) < EXPONENT_MIN as i64 + shift => {
                (p as i64 - (EXPONENT_MIN as i64 + shift - widen(e))).max(1) as u32
            }
            _ => p,
        };
        let mut w2 = Float::new(bits + 1);
        let t2 = unsafe { mpfr::set(w2.as_raw_mut(), w.as_raw(), rnd_t::RNDZ) };
        let t = if t2 != 0 { 1 } else { t };
        let mut lo = Float::new(bits);
        let mut hi = Float::new(bits);
        unsafe {
            mpfr::set(lo.as_raw_mut(), w2.as_raw(), rnd_t::RNDZ);
            mpfr::set(hi.as_raw_mut(), w2.as_raw(), rnd_t::RNDA);
        }
        let odd =
            |x: &Float, n: u32| !x.is_zero() && unsafe { mpfr::min_prec(x.as_raw()) } as u32 == n;
        let chosen = if lo == hi {
            lo
        } else if t == 0 {
            // w2 is the exact midpoint: the tie goes to the odd neighbour.
            if odd(&lo, bits) {
                lo
            } else {
                hi
            }
        } else if odd(&w2, bits + 1) {
            // the exact value lies above the midpoint
            hi
        } else {
            lo
        };
        let inexact = t != 0 || unsafe { mpfr::min_prec(w2.as_raw()) } as u32 > bits;
        let mut v = Float::new(p);
        unsafe { mpfr::set(v.as_raw_mut(), chosen.as_raw(), rnd_t::RNDZ) };
        (v, inexact)
    };
    unsafe {
        mpfr::set_emin(old_min);
        mpfr::set_emax(old_max);
    }
    Ok(out)
}

/// MPFR `(emin, emax)` that emulates the zenith range scaled by `2^shift` at precision `p`, clamped to
/// what MPFR allows.
pub fn emulated_range(p: u32, shift: i64) -> (i64, i64) {
    let (lo, hi) = unsafe { (widen(mpfr::get_emin_min()), widen(mpfr::get_emax_max())) };
    let emin = (EXPONENT_MIN as i64 + shift - p as i64 + 1).max(lo);
    let emax = (EXPONENT_MAX as i64).saturating_add(shift).min(hi);
    (emin, emax)
}

/// True if MPFR can emulate the zenith subnormal floor scaled by `2^shift` at precision `p` without
/// clamping. A positive shift may clamp the top of the scaled range: callers using one (the bottom-edge
/// tests) keep their results far below `EXPONENT_MAX`, so only the unshifted top must fit.
pub fn emulated_range_ok(p: u32, shift: i64) -> bool {
    let (lo, hi) = unsafe { (widen(mpfr::get_emin_min()), widen(mpfr::get_emax_max())) };
    EXPONENT_MIN as i64 + shift - p as i64 + 1 >= lo && (EXPONENT_MAX as i64).saturating_add(shift.min(0)) <= hi
}

/// The exact value of a finite zenith number times `2^shift` as an MPFR number (widest exponent range
/// while building; the caller's range must hold the result).
pub fn exact_mpfr_shifted(n: &ExactNum, shift: i64) -> Result<Float, Box<dyn std::error::Error>> {
    if n.is_zero() {
        let z = Float::with_val(64, 0);
        return Ok(if n.is_negative() { -z } else { z });
    }
    let (words, _, sign, e, _) = n.as_raw_parts().ok_or("not a finite number")?;
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
    let sh = e as i64 + shift - 64 * words.len() as i64;
    unsafe {
        mpfr::mul_2si(got.as_raw_mut(), got.as_raw(), sh as _, rnd_t::RNDN);
    }
    if sign == Sign::Neg {
        got = -got;
    }
    unsafe {
        mpfr::set_emin(old_min);
        mpfr::set_emax(old_max);
    }
    Ok(got)
}

/// An MPFR exponent (`c_long`: 32 or 64 bits by target) as an `i64`.
fn widen<T: Into<i64>>(x: T) -> i64 {
    x.into()
}

/// Bit-exact comparison of a zenith result with an MPFR value of the same precision,
/// including Inf, the sign of zero, and subnormals.
pub fn assert_exact_eq(
    n: &ExactNum,
    f: &Float,
    p: usize,
    op: &str,
    cc: &mut Consts,
) -> Result<(), Box<dyn std::error::Error>> {
    if f.is_nan() {
        assert!(n.is_nan(), "{op}: MPFR NaN, zenith {n:?}");
        return Ok(());
    }
    if f.is_infinite() {
        assert!(n.is_inf(), "{op}: MPFR {}, zenith {n:?}", show_mp(f));
        assert_eq!(n.is_negative(), f.is_sign_negative(), "{op}: sign of Inf");
        return Ok(());
    }
    assert!(!n.is_inf() && !n.is_nan(), "{op}: MPFR {}, zenith {n:?}", show_mp(f));
    if f.is_zero() {
        assert!(n.is_zero(), "{op}: MPFR {}, zenith {n:?}", show_mp(f));
        assert_eq!(n.is_negative(), f.is_sign_negative(), "{op}: sign of zero");
        return Ok(());
    }
    let _ = (p, &cc);
    // Exact value of the zenith number, rebuilt from its words (parsing would flush
    // subnormals): sum of word_i * 2^(64 i) * 2^(e - 64 len), with the widest MPFR range.
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
    assert!(equal, "{op}: MPFR {}, zenith {n:?} = {}", show_mp(f), show_mp(&got));
    Ok(())
}

/// `m * 2^e` in hex with MPFR's raw 64-bit exponent: decimal (or long binary) formatting of a
/// number near 2^(+-2^62) would allocate exabytes.
pub fn show_mp(f: &Float) -> String {
    if !f.is_normal() {
        return format!("{f}");
    }
    let e = unsafe { (*f.as_raw()).exp };
    with_full_range(|| {
        let frac = Float::with_val(f.prec(), f >> (e as isize));
        format!("{frac:x} * 2^{e}")
    })
}
