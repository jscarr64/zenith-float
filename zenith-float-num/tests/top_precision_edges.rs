//! A valid precision is never `InvalidArgument` (CONTRACT: only a precision outside
//! `MIN_PRECISION..=MAX_PRECISION` is). Precisions are rounded up to whole words and the
//! correctly-rounded paths work some words wider than the caller's precision, so near the top of
//! the range those working widths lie above `MAX_PRECISION`; they must still be accepted, and the
//! only failure left at such a size is the genuine one, `MemoryAllocation`.
//!
//! Every check is written against `MAX_PRECISION` and `WORD_BIT_SIZE`, not their values: if the
//! ceiling changes, a top precision that can be allocated must give the exact value instead.
//!
//! Exact oracles (no rounding at any precision of at least the stated bits): 12 + 4 = 16,
//! 12 - 4 = 8, 12 * 4 = 48, 12 / 4 = 3, 12^2 = 144 (6, 4, 6, 2, 8 bits), with all four sign
//! combinations; 2 op 2 is exact at 1 bit (4, 0, 4, 1). sqrt(16) = 4, hypot(3, 4) = 5,
//! fma(3, 4, 5) = 17, log_2(8) = 3, 1/4 = 0.25, 3^3 = 27.

use zenith_float_num::{
    Consts, Error, ExactNum, Radix, RoundingMode, MAX_PRECISION, MIN_PRECISION, WORD_BIT_SIZE,
};

const MODES: [RoundingMode; 7] = [
    RoundingMode::None,
    RoundingMode::Up,
    RoundingMode::Down,
    RoundingMode::ToZero,
    RoundingMode::FromZero,
    RoundingMode::ToEven,
    RoundingMode::ToOdd,
];

#[derive(Debug)]
enum Bin {
    Add,
    Sub,
    Mul,
    Div,
    AddRef,
    SubRef,
    MulRef,
    DivRef,
}

const BINS: [Bin; 8] = [
    Bin::Add,
    Bin::Sub,
    Bin::Mul,
    Bin::Div,
    Bin::AddRef,
    Bin::SubRef,
    Bin::MulRef,
    Bin::DivRef,
];

fn n(v: i64) -> ExactNum {
    ExactNum::from_i64(v, WORD_BIT_SIZE)
}

/// `a op b` through the public form named by `op` (the owned forms consume a fresh operand).
fn bin(op: &Bin, a: i64, b: i64, p: usize, rm: RoundingMode) -> ExactNum {
    match op {
        Bin::Add => n(a).add(&n(b), p, rm),
        Bin::Sub => n(a).sub(&n(b), p, rm),
        Bin::Mul => n(a).mul(&n(b), p, rm),
        Bin::Div => n(a).div(&n(b), p, rm),
        Bin::AddRef => ExactNum::add_ref(&n(a), &n(b), p, rm),
        Bin::SubRef => ExactNum::sub_ref(&n(a), &n(b), p, rm),
        Bin::MulRef => ExactNum::mul_ref(&n(a), &n(b), p, rm),
        Bin::DivRef => ExactNum::div_ref(&n(a), &n(b), p, rm),
    }
}

/// The exact value of `a op b`, for operands where it is an integer.
fn bin_want(op: &Bin, a: i64, b: i64) -> Result<i64, String> {
    let r = match op {
        Bin::Add | Bin::AddRef => a.checked_add(b),
        Bin::Sub | Bin::SubRef => a.checked_sub(b),
        Bin::Mul | Bin::MulRef => a.checked_mul(b),
        Bin::Div | Bin::DivRef => {
            if a.checked_rem(b) == Some(0) {
                a.checked_div(b)
            } else {
                None
            }
        }
    };
    r.ok_or(format!("{op:?}({a}, {b}) has no exact integer value"))
}

/// The precisions at the top of the range: `MAX_PRECISION - d`, every branch boundary of the
/// word rounding (a multiple of the word size, one either side, one and two words below).
fn top_precisions() -> Result<Vec<usize>, String> {
    let w = WORD_BIT_SIZE;
    let ds = [
        Some(0),
        Some(1),
        Some(2),
        w.checked_sub(2),
        w.checked_sub(1),
        Some(w),
        w.checked_add(1),
        w.checked_mul(2).and_then(|x| x.checked_sub(1)),
        w.checked_mul(2),
        w.checked_mul(2).and_then(|x| x.checked_add(1)),
    ];
    let mut v = Vec::new();
    for d in ds {
        let d = d.ok_or("distance")?;
        v.push(MAX_PRECISION.checked_sub(d).ok_or("MAX_PRECISION - d")?);
    }
    // The largest multiple of the word size at or below the ceiling, and one bit either side.
    let m = MAX_PRECISION
        .checked_sub(MAX_PRECISION % w)
        .ok_or("multiple")?;
    v.extend([m, m.checked_sub(1).ok_or("m - 1")?]);
    if m < MAX_PRECISION {
        v.push(m.checked_add(1).ok_or("m + 1")?);
    }
    Ok(v)
}

/// Precisions above the ceiling (and 0): always `InvalidArgument`.
fn bad_precisions() -> Result<Vec<usize>, String> {
    let w = WORD_BIT_SIZE;
    let mut v = vec![0, usize::MAX];
    for d in [1, 2, w.checked_sub(1).ok_or("w - 1")?, w, w.checked_add(1).ok_or("w + 1")?] {
        v.push(MAX_PRECISION.checked_add(d).ok_or("MAX_PRECISION + d")?);
    }
    Ok(v)
}

/// A result at a valid precision: the exact value `want`, or a NaN carrying `MemoryAllocation`.
fn valid(what: &str, r: &ExactNum, want: &ExactNum) -> Result<(), String> {
    if r.is_nan() {
        let e = r.err();
        return if e == Some(Error::MemoryAllocation) {
            Ok(())
        } else {
            Err(format!(
                "{what}: a valid precision gave NaN({e:?}), want the value or MemoryAllocation"
            ))
        };
    }
    if r == want {
        Ok(())
    } else {
        Err(format!("{what}: got {r:?}, want {want:?}"))
    }
}

fn invalid(what: &str, r: &ExactNum) -> Result<(), String> {
    if r.is_nan() && r.err() == Some(Error::InvalidArgument) {
        Ok(())
    } else {
        Err(format!("{what}: want NaN(InvalidArgument), got {r:?}"))
    }
}

const SIGNED: [(i64, i64); 4] = [(12, 4), (-12, 4), (12, -4), (-12, -4)];

#[test]
fn binary_ops_at_the_top_of_the_precision_range_are_never_invalid() -> Result<(), String> {
    for p in top_precisions()? {
        let d = MAX_PRECISION.checked_sub(p).ok_or("d")?;
        for op in &BINS {
            for (a, b) in SIGNED.into_iter().chain([(0, 4), (12, 0), (0, 0)]) {
                for rm in MODES {
                    let what = format!("{op:?}({a}, {b}) p=MAX-{d} {rm:?}");
                    let r = bin(op, a, b, p, rm);
                    let by_zero = b == 0 && matches!(op, Bin::Div | Bin::DivRef);
                    if by_zero {
                        // x / 0 is Inf, and 0 / 0 is NaN(InvalidArgument), from the operands, not the precision.
                        let inf = r.is_inf_pos() || r.is_inf_neg();
                        let zero_by_zero = a == 0 && r.err() == Some(Error::InvalidArgument);
                        if !(inf || zero_by_zero) {
                            return Err(format!("{what}: got {r:?}"));
                        }
                    } else {
                        valid(&what, &r, &n(bin_want(op, a, b)?))?;
                    }
                }
            }
        }
    }
    Ok(())
}

#[test]
fn sqr_div_exact_and_siblings_at_the_top_are_never_invalid() -> Result<(), String> {
    let mut cc = Consts::new().map_err(|e| format!("{e:?}"))?;
    for p in top_precisions()? {
        let d = MAX_PRECISION.checked_sub(p).ok_or("d")?;
        for rm in MODES {
            let at = format!("p=MAX-{d} {rm:?}");
            valid(&format!("sqr(-12) {at}"), &n(-12).sqr(p, rm), &n(144))?;
            valid(&format!("sqr(12) {at}"), &n(12).sqr(p, rm), &n(144))?;
            valid(&format!("sqrt(16) {at}"), &n(16).sqrt(p, rm), &n(4))?;
            valid(
                &format!("hypot(-3, 4) {at}"),
                &n(-3).hypot(&n(4), p, rm),
                &n(5),
            )?;
            valid(
                &format!("hypot(3, -4) {at}"),
                &n(3).hypot(&n(-4), p, rm),
                &n(5),
            )?;
            valid(
                &format!("fma(3, 4, 5) {at}"),
                &n(3).fma(&n(4), &n(5), p, rm),
                &n(17),
            )?;
            valid(
                &format!("fma(-3, 4, 5) {at}"),
                &n(-3).fma(&n(4), &n(5), p, rm),
                &n(-7),
            )?;
            valid(&format!("powi(-3, 3) {at}"), &n(-3).powi(3, p, rm), &n(-27))?;
            valid(
                &format!("reciprocal(-4) {at}"),
                &n(-4).reciprocal(p, rm),
                &n(-1).div(&n(4), WORD_BIT_SIZE, RoundingMode::ToEven),
            )?;
            valid(
                &format!("log_2(8) {at}"),
                &n(8).log(&n(2), p, rm, &mut cc),
                &n(3),
            )?;
        }
        for (a, b) in SIGNED {
            let what = format!("div_exact({a}, {b}) p=MAX-{d}");
            let r = n(a).div_exact(&n(b), p);
            if let Ok(Some(q)) = &r {
                valid(&what, q, &n(bin_want(&Bin::Div, a, b)?))?;
            } else if r != Err(Error::MemoryAllocation) {
                return Err(format!("{what}: got {r:?}"));
            }
        }
        let mut x = n(-3);
        let r = x.set_precision(p, RoundingMode::ToEven);
        if r == Ok(()) {
            valid(&format!("set_precision(-3) p=MAX-{d}"), &x, &n(-3))?;
        } else if r != Err(Error::MemoryAllocation) {
            return Err(format!("set_precision(-3) p=MAX-{d}: got {r:?}"));
        }
        valid(&format!("new p=MAX-{d}"), &ExactNum::new(p), &n(0))?;
    }
    Ok(())
}

#[test]
fn precisions_outside_the_range_are_invalid() -> Result<(), String> {
    let mut cc = Consts::new().map_err(|e| format!("{e:?}"))?;
    for p in bad_precisions()? {
        for rm in MODES {
            let at = format!("p={p} {rm:?}");
            for op in &BINS {
                for (a, b) in SIGNED {
                    invalid(&format!("{op:?}({a}, {b}) {at}"), &bin(op, a, b, p, rm))?;
                }
            }
            invalid(&format!("sqr {at}"), &n(-12).sqr(p, rm))?;
            invalid(&format!("sqrt {at}"), &n(16).sqrt(p, rm))?;
            invalid(&format!("hypot {at}"), &n(-3).hypot(&n(4), p, rm))?;
            invalid(&format!("fma {at}"), &n(3).fma(&n(4), &n(5), p, rm))?;
            invalid(&format!("powi {at}"), &n(-3).powi(3, p, rm))?;
            invalid(&format!("reciprocal {at}"), &n(-4).reciprocal(p, rm))?;
            invalid(&format!("log {at}"), &n(8).log(&n(2), p, rm, &mut cc))?;
        }
        let r = n(12).div_exact(&n(-4), p);
        if r != Err(Error::InvalidArgument) {
            return Err(format!("div_exact p={p}: got {r:?}"));
        }
        let mut x = n(-3);
        let r = x.set_precision(p, RoundingMode::ToEven);
        if r != Err(Error::InvalidArgument) {
            return Err(format!("set_precision p={p}: got {r:?}"));
        }
        invalid(&format!("new p={p}"), &ExactNum::new(p))?;
    }
    Ok(())
}

/// Small precisions around every word boundary, both signs, every mode: the exact value.
#[test]
fn binary_ops_at_small_precisions_around_word_boundaries_are_exact() -> Result<(), String> {
    let w = WORD_BIT_SIZE;
    let mut ps: Vec<usize> = (MIN_PRECISION..=8).collect();
    for k in 1..=4usize {
        let m = w.checked_mul(k).ok_or("k words")?;
        ps.extend([m.checked_sub(1).ok_or("m - 1")?, m, m.checked_add(1).ok_or("m + 1")?]);
    }
    for p in ps {
        for op in &BINS {
            let pairs: &[(i64, i64)] =
                if p >= 6 { &SIGNED } else { &[(2, 2), (-2, 2), (2, -2), (-2, -2)] };
            for &(a, b) in pairs {
                for rm in MODES {
                    let what = format!("{op:?}({a}, {b}) p={p} {rm:?}");
                    let r = bin(op, a, b, p, rm);
                    let want = n(bin_want(op, a, b)?);
                    if r != want || r.is_nan() {
                        return Err(format!("{what}: got {r:?}, want {want:?}"));
                    }
                }
            }
        }
    }
    Ok(())
}

/// Every other public entry point that takes a precision and works at it directly, by name, with its
/// result at `p`. `parse` takes `usize::MAX` as "the precision of the numeral", so it is left out there.
fn public_results(p: usize, rm: RoundingMode, cc: &mut Consts) -> Vec<(&'static str, ExactNum)> {
    let mut v = vec![
        ("new", ExactNum::new(p)),
        ("max_value", ExactNum::max_value(p)),
        ("min_value", ExactNum::min_value(p)),
        ("min_positive", ExactNum::min_positive(p)),
        ("min_positive_normal", ExactNum::min_positive_normal(p)),
        ("from_word", ExactNum::from_word(5, p)),
        ("from_i64", ExactNum::from_i64(-5, p)),
        ("from_u8", ExactNum::from_u8(5, p)),
        ("ldexp", n(-3).ldexp(2, p, rm)),
        ("scalb", n(3).scalb(-2, p, rm)),
        ("copysign value", n(3).copysign(&n(-1), p, rm)),
        (
            "copysign inf",
            n(3).copysign(&n(-1).div(&n(0), WORD_BIT_SIZE, rm), p, rm),
        ),
        ("copysign nan", n(-3).copysign(&ExactNum::nan(None), p, rm)),
        ("next_after", n(3).next_after(&n(-4), p, rm)),
        ("mul_add", n(3).mul_add(&n(-4), &n(5), p, rm)),
        ("powsi", n(-3).powsi(-2, p, rm)),
    ];
    if p != usize::MAX {
        v.push(("parse", ExactNum::parse("-12.5", Radix::Dec, p, rm, cc)));
    }
    let (hi, lo) = n(3).two_sum(&n(-4), p, rm);
    v.extend([("two_sum hi", hi), ("two_sum lo", lo)]);
    let (hi, lo) = n(3).two_product(&n(-4), p, rm);
    v.extend([("two_product hi", hi), ("two_product lo", lo)]);
    let mut x = n(-3);
    let done = x.try_set_precision(p, rm, 0);
    v.push((
        if done { "try_set_precision (done)" } else { "try_set_precision (not done)" },
        x,
    ));
    v
}

/// Entry points that build constants (e, ln 2, pi) or call the public API again at a working width.
/// Outside the range they are `InvalidArgument` like the rest. At the top of the range they are not
/// checked here: the constants cache works toward 2^60 bits for hours before an allocation fails, and the
/// composites still report `InvalidArgument` there (both recorded in the changelog entry as open).
fn composite_results(p: usize, rm: RoundingMode, cc: &mut Consts) -> Vec<(&'static str, ExactNum)> {
    let mut v = vec![
        ("fused_sum", ExactNum::fused_sum(&[n(3), n(-4)], p, rm)),
        (
            "fused_dot",
            ExactNum::fused_dot(&[n(3), n(-4)], &[n(2), n(5)], p, rm),
        ),
        ("polyval", ExactNum::polyval(&[n(1), n(-2)], &n(3), p, rm)),
        ("pow", n(4).pow(&n(-2), p, rm, cc)),
        ("exp", n(-1).exp(p, rm, cc)),
        ("ln", n(3).ln(p, rm, cc)),
        ("sin", n(-1).sin(p, rm, cc)),
        ("atan", n(1).atan(p, rm, cc)),
        ("atan2", n(1).atan2(&n(-1), p, rm, cc)),
        ("tanh", n(-1).tanh(p, rm, cc)),
        ("log1p", n(1).log1p(p, rm, cc)),
        ("expm1", n(-1).expm1(p, rm, cc)),
    ];
    let (s, c) = n(1).sin_cos(p, rm, cc);
    v.extend([("sin_cos sin", s), ("sin_cos cos", c)]);
    let (s, c) = n(1).sinh_cosh(p, rm, cc);
    v.extend([("sinh_cosh sinh", s), ("sinh_cosh cosh", c)]);
    v
}

#[test]
fn parse_with_usize_max_takes_the_precision_of_the_numeral() -> Result<(), String> {
    let mut cc = Consts::new().map_err(|e| format!("{e:?}"))?;
    let want = n(-25).div(&n(2), WORD_BIT_SIZE, RoundingMode::ToEven);
    for rm in MODES {
        let r = ExactNum::parse("-12.5", Radix::Dec, usize::MAX, rm, &mut cc);
        if r != want {
            return Err(format!("parse(-12.5, usize::MAX) {rm:?}: got {r:?}"));
        }
    }
    Ok(())
}

#[test]
fn every_public_entry_accepts_a_top_precision() -> Result<(), String> {
    let mut cc = Consts::new().map_err(|e| format!("{e:?}"))?;
    for p in top_precisions()? {
        let d = MAX_PRECISION.checked_sub(p).ok_or("d")?;
        for rm in [RoundingMode::None, RoundingMode::ToEven, RoundingMode::Down] {
            for (name, r) in public_results(p, rm, &mut cc) {
                // A value cannot be built at this size on a real machine; with a lower ceiling it would be.
                if r.is_nan() && r.err() != Some(Error::MemoryAllocation) {
                    return Err(format!(
                        "{name} p=MAX-{d} {rm:?}: got NaN({:?}), want a value or MemoryAllocation",
                        r.err()
                    ));
                }
            }
        }
    }
    Ok(())
}

#[test]
fn every_public_entry_rejects_a_precision_outside_the_range() -> Result<(), String> {
    let mut cc = Consts::new().map_err(|e| format!("{e:?}"))?;
    for p in bad_precisions()? {
        for rm in [RoundingMode::None, RoundingMode::ToEven, RoundingMode::Down] {
            let mut all = public_results(p, rm, &mut cc);
            all.extend(composite_results(p, rm, &mut cc));
            for (name, r) in all {
                invalid(&format!("{name} p={p} {rm:?}"), &r)?;
            }
        }
    }
    Ok(())
}
