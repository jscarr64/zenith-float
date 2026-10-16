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
    zenith_float_num::reseed_random(zenith_float_num::DEFAULT_RANDOM_SEED);
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

macro_rules! test_zf_op_no_cc {
    ($eq:literal, $n1:ident, $n2:ident, $zf_op:ident, $f1:ident, $f2:ident, $mpfr_op:ident, $p:ident, $rm:ident, $rnd:ident, $op_info:expr, $cc:ident) => {
        let n3 = ($n1).clone().$zf_op(&($n2), $p, $rm);

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
        assert!(f.is_infinite(), "{} zenith Inf, MPFR {}", op, f);
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
            assert!(unf != 0, "{}", op);
        }
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
    let (old_min, old_max) = unsafe { (mpfr::get_emin(), mpfr::get_emax()) };
    unsafe {
        mpfr::set_emax(EXPONENT_MAX as mpfr::exp_t);
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
        let mut w = Float::new(p + 1);
        let t = f(&mut w, rnd_t::RNDZ);
        let mut lo = Float::new(p);
        let mut hi = Float::new(p);
        unsafe {
            mpfr::set(lo.as_raw_mut(), w.as_raw(), rnd_t::RNDZ);
            mpfr::set(hi.as_raw_mut(), w.as_raw(), rnd_t::RNDA);
        }
        let odd = |x: &Float, bits: u32| {
            !x.is_zero() && unsafe { mpfr::min_prec(x.as_raw()) } as u32 == bits
        };
        let v = if lo == hi {
            lo
        } else if t == 0 {
            // w is the exact midpoint: the tie goes to the odd neighbour.
            if odd(&lo, p) {
                lo
            } else {
                hi
            }
        } else if odd(&w, p + 1) {
            // the exact value lies above the midpoint
            hi
        } else {
            lo
        };
        let inexact = t != 0 || unsafe { mpfr::min_prec(w.as_raw()) } as u32 > p;
        (v, inexact)
    };
    unsafe {
        mpfr::set_emin(old_min);
        mpfr::set_emax(old_max);
    }
    Ok(out)
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
    assert!(equal, "{op}: MPFR {f}, zenith {n:?} = {got}");
    Ok(())
}
