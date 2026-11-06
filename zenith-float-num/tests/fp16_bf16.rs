//! zenith-float 2.0 Slice 1: `Ieee16` (IEEE 754 binary16) and `Bfloat16` (bfloat16) `u16` bit patterns
//! to and from `ExactNum`, with every `RoundingMode`.
//!
//! Oracles (no hardware float anywhere):
//! A. Decode: every one of the 65,536 patterns of each format is checked against the exact integer
//!    relation `|x| = sig * 2^k` from the IEEE 754-2019 §3.4 field layout (binary16: 5 exponent bits,
//!    bias 15, 10 fraction bits; bfloat16: 8 exponent bits, bias 127, 7 fraction bits), so
//!    `|x| * 2^-k` must compare equal to the integer `sig`. Predicates follow the same field decode.
//! B. Round trip: every pattern widened to `ExactNum` (p = 64 and p = 1024, so the low limbs are zero
//!    words) comes back bit-identical in all seven modes; NaN comes back as the canonical quiet NaN
//!    (0x7E00 / 0x7FC0), -0 keeps its sign.
//! C. Brackets by construction: for every pair of adjacent non-negative patterns (lo, lo + 1), including
//!    (max finite, +Inf) where the upper neighbour is 2^(emax+1), the exact midpoint, the midpoint plus
//!    and minus delta * 2^-900, and lo + delta * 2^-900 (delta = the gap) are rounded in all seven modes
//!    and both signs. The tiny offset sits ~900 bits below the leading bit, i.e. in the low limbs of a
//!    1024-bit operand, so a sticky bit read from only the top word fails here. The expected pattern
//!    follows from IEEE 754-2019 §4.3 directly: ties-to-even / ties-to-odd pick by the pattern's low bit,
//!    directed modes pick by sign, RoundingMode::None refuses an inexact input.
//! D. Binary-search oracle: 4,000 pseudo-random 1024-bit operands per format (fixed LCG seed, exponents
//!    spanning below the subnormal range to above the overflow threshold) are bracketed by searching the
//!    ordered patterns with `ExactNum::compare`, and rounded by comparison with the exact midpoint.
//! E. Extremes: 2^-10000 and 2^10000 (both signs) underflow / overflow per mode; NaN, ±Inf, ±0.
//! F. Golds: 1, 65504, 2^-14, 2^-24 (binary16); 1, (2^8 - 1) * 2^120, 2^-126, 2^-133 (bfloat16):
//!    the format limits from IEEE 754-2019 Table 3.5 and the bfloat16 layout (binary32 truncated to 16
//!    bits).

use zenith_float_num::{Bfloat16, Error, ExactNum, Ieee16, RoundingMode};

const P: usize = 1024;
const MODES: [RoundingMode; 7] = [
    RoundingMode::ToEven,
    RoundingMode::ToOdd,
    RoundingMode::Up,
    RoundingMode::Down,
    RoundingMode::ToZero,
    RoundingMode::FromZero,
    RoundingMode::None,
];

struct Fmt {
    frac: u32,
    bias: i64,
    qnan: u16,
    max: u16,
    to_exact: fn(u16, usize) -> ExactNum,
    from_rm: fn(&ExactNum, RoundingMode) -> Result<u16, Error>,
    from_even: fn(&ExactNum) -> u16,
    preds: fn(u16) -> [bool; 7],
}

fn h_preds(b: u16) -> [bool; 7] {
    let x = Ieee16::from_bits(b);
    [
        x.is_nan(),
        x.is_infinite(),
        x.is_finite(),
        x.is_zero(),
        x.is_subnormal(),
        x.is_sign_negative(),
        x.classify() == core::num::FpCategory::Normal,
    ]
}

fn b_preds(b: u16) -> [bool; 7] {
    let x = Bfloat16::from_bits(b);
    [
        x.is_nan(),
        x.is_infinite(),
        x.is_finite(),
        x.is_zero(),
        x.is_subnormal(),
        x.is_sign_negative(),
        x.classify() == core::num::FpCategory::Normal,
    ]
}

const H: Fmt = Fmt {
    frac: 10,
    bias: 15,
    qnan: 0x7E00,
    max: 0x7BFF,
    to_exact: |b, p| Ieee16::from_bits(b).to_exact(p),
    from_rm: |x, rm| Ieee16::from_exact_rm(x, rm).map(Ieee16::to_bits),
    from_even: |x| Ieee16::from_exact(x).to_bits(),
    preds: h_preds,
};

const B: Fmt = Fmt {
    frac: 7,
    bias: 127,
    qnan: 0x7FC0,
    max: 0x7F7F,
    to_exact: |b, p| Bfloat16::from_bits(b).to_exact(p),
    from_rm: |x, rm| Bfloat16::from_exact_rm(x, rm).map(Bfloat16::to_bits),
    from_even: |x| Bfloat16::from_exact(x).to_bits(),
    preds: b_preds,
};

fn pow2(k: i64) -> ExactNum {
    ExactNum::from_u64(1, P).ldexp(k, P, RoundingMode::None)
}

fn eq(a: &ExactNum, b: &ExactNum) -> bool {
    a.compare(b) == Some(0)
}

/// Exact value of a finite non-negative pattern, or 2^(emax+1) for the +Inf pattern (the upper
/// neighbour of max finite when the exponent range is unbounded).
fn upper_val(f: &Fmt, b: u16) -> ExactNum {
    if b == f.max + 1 {
        pow2(f.bias + 1)
    } else {
        (f.to_exact)(b, P)
    }
}

/// Field decode, independent of the library: (class, negative, sig, k) with |x| = sig * 2^k.
/// class: 0 zero, 1 subnormal, 2 normal, 3 inf, 4 nan.
fn decode(f: &Fmt, b: u16) -> (u8, bool, u64, i64) {
    let neg = b >> 15 == 1;
    let ebits = 15 - f.frac;
    let e = i64::from((b >> f.frac) & ((1u16 << ebits) - 1));
    let m = u64::from(b & ((1u16 << f.frac) - 1));
    let eall = (1i64 << ebits) - 1;
    let fr = i64::from(f.frac);
    if e == eall {
        (if m == 0 { 3 } else { 4 }, neg, m, 0)
    } else if e == 0 {
        (if m == 0 { 0 } else { 1 }, neg, m, 1 - f.bias - fr)
    } else {
        (2, neg, m + (1u64 << f.frac), e - f.bias - fr)
    }
}

fn decode_exhaustive(f: &Fmt) {
    for b in 0..=u16::MAX {
        let (class, neg, sig, k) = decode(f, b);
        let x = (f.to_exact)(b, 64);
        let pr = (f.preds)(b);
        assert_eq!(
            pr,
            [class == 4, class == 3, class < 3, class == 0, class == 1, neg, class == 2],
            "predicates {b:#06x}"
        );
        match class {
            4 => assert!(x.is_nan(), "{b:#06x}"),
            3 => assert!(x.is_inf() && x.is_negative() == neg, "{b:#06x}"),
            0 => assert!(x.is_zero() && x.is_negative() == neg, "{b:#06x}"),
            _ => {
                assert_eq!(x.is_negative(), neg, "{b:#06x}");
                let scaled = x.abs().ldexp(-k, P, RoundingMode::None);
                assert!(eq(&scaled, &ExactNum::from_u64(sig, P)), "{b:#06x}");
            }
        }
    }
}

#[test]
fn ieee16_decode_exhaustive() {
    decode_exhaustive(&H);
}

#[test]
fn bfloat16_decode_exhaustive() {
    decode_exhaustive(&B);
}

fn roundtrip_exhaustive(f: &Fmt) {
    for p in [64usize, P] {
        for b in 0..=u16::MAX {
            let (class, _, _, _) = decode(f, b);
            let want = if class == 4 { f.qnan } else { b };
            let x = (f.to_exact)(b, p);
            for rm in MODES {
                assert_eq!((f.from_rm)(&x, rm), Ok(want), "{b:#06x} {rm:?} p={p}");
            }
            assert_eq!((f.from_even)(&x), want, "{b:#06x} p={p}");
        }
    }
}

#[test]
fn ieee16_roundtrip_exhaustive_all_modes() {
    roundtrip_exhaustive(&H);
}

#[test]
fn bfloat16_roundtrip_exhaustive_all_modes() {
    roundtrip_exhaustive(&B);
}

/// IEEE 754-2019 §4.3 choice between the magnitudes `lo` and `lo + 1` for an inexact operand.
/// `cmp` is the operand against the midpoint (-1 below, 0 tie, 1 above).
fn pick(lo: u16, cmp: i32, neg: bool, rm: RoundingMode) -> Result<u16, Error> {
    let hi = lo + 1;
    let s = if neg { 0x8000 } else { 0 };
    let mag = match rm {
        RoundingMode::None => return Err(Error::InvalidArgument),
        RoundingMode::ToEven | RoundingMode::ToOdd => {
            if cmp < 0 {
                lo
            } else if cmp > 0 {
                hi
            } else if lo.is_multiple_of(2) == (rm == RoundingMode::ToEven) {
                lo
            } else {
                hi
            }
        }
        RoundingMode::Up => {
            if neg {
                lo
            } else {
                hi
            }
        }
        RoundingMode::Down => {
            if neg {
                hi
            } else {
                lo
            }
        }
        RoundingMode::ToZero => lo,
        RoundingMode::FromZero => hi,
    };
    Ok(s | mag)
}

fn check(f: &Fmt, x: &ExactNum, lo: u16, cmp: i32, neg: bool, tag: &str) {
    for rm in MODES {
        let want = pick(lo, cmp, neg, rm);
        assert_eq!(
            (f.from_rm)(x, rm),
            want,
            "{tag} lo={lo:#06x} neg={neg} {rm:?}"
        );
        if rm == RoundingMode::ToEven {
            assert_eq!(Ok((f.from_even)(x)), want, "{tag} lo={lo:#06x} neg={neg}");
        }
    }
}

fn brackets_by_construction(f: &Fmt) {
    for lo in 0..=f.max {
        let delta = upper_val(f, lo + 1).sub(&(f.to_exact)(lo, P), P, RoundingMode::ToEven);
        let tiny = || delta_tiny(f, lo);
        let mid = || {
            upper_val(f, lo + 1)
                .add(&(f.to_exact)(lo, P), P, RoundingMode::ToEven)
                .ldexp(-1, P, RoundingMode::None)
        };
        assert!(!delta.is_zero() && !delta.is_negative());
        for neg in [false, true] {
            let sgn = |x: ExactNum| if neg { -x } else { x };
            check(f, &sgn(mid()), lo, 0, neg, "tie");
            check(
                f,
                &sgn(mid().add(&tiny(), P, RoundingMode::ToEven)),
                lo,
                1,
                neg,
                "above",
            );
            check(
                f,
                &sgn(mid().sub(&tiny(), P, RoundingMode::ToEven)),
                lo,
                -1,
                neg,
                "below",
            );
            check(
                f,
                &sgn((f.to_exact)(lo, P).add(&tiny(), P, RoundingMode::ToEven)),
                lo,
                -1,
                neg,
                "lo+tiny",
            );
        }
    }
}

/// (gap above `lo`) * 2^-900.
fn delta_tiny(f: &Fmt, lo: u16) -> ExactNum {
    upper_val(f, lo + 1)
        .sub(&(f.to_exact)(lo, P), P, RoundingMode::ToEven)
        .ldexp(-900, P, RoundingMode::None)
}

#[test]
fn ieee16_ties_sticky_subnormal_overflow_brackets() {
    brackets_by_construction(&H);
}

#[test]
fn bfloat16_ties_sticky_subnormal_overflow_brackets() {
    brackets_by_construction(&B);
}

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        // Knuth MMIX constants.
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }
}

fn search_oracle(f: &Fmt, seed: u64) {
    let mut r = Lcg(seed);
    let emin = 1 - f.bias;
    let span = u64::try_from(f.bias - emin + i64::from(f.frac) + 24).unwrap_or(1);
    for i in 0..4000 {
        // 1024-bit mantissa in [1, 2): 16 random words, top bit forced.
        let mut m =
            ExactNum::from_u64(r.next() | (1u64 << 63), P).ldexp(-63, P, RoundingMode::None);
        for w in 1..16i64 {
            let piece = ExactNum::from_u64(r.next(), P).ldexp(-63 - 64 * w, P, RoundingMode::None);
            m = m.add(&piece, P, RoundingMode::ToEven);
        }
        // Every 7th operand keeps only its top 40 bits plus a lone low bit, so exact
        // ties and near-ties with a far sticky bit also show up.
        if i % 7 == 0 {
            m = m
                .ldexp(39, P, RoundingMode::None)
                .floor()
                .ldexp(-39, P, RoundingMode::None)
                .add(&pow2(-1000), P, RoundingMode::ToEven);
        }
        let off = i64::try_from(r.next() % span).unwrap_or(0);
        let e = emin - i64::from(f.frac) - 12 + off;
        let x = m.ldexp(e, P, RoundingMode::None);
        // Largest non-negative pattern lo <= max with value(lo) <= x.
        let (mut a, mut b) = (0u16, f.max);
        while a < b {
            let c = a + (b - a).div_ceil(2);
            if (f.to_exact)(c, P).compare(&x).unwrap_or(1) <= 0 {
                a = c;
            } else {
                b = c - 1;
            }
        }
        let lo = a;
        let neg = r.next() % 2 == 1;
        let sx = || {
            let v = m_copy(&x);
            if neg {
                -v
            } else {
                v
            }
        };
        if eq(&(f.to_exact)(lo, P), &x) {
            for rm in MODES {
                let s = if neg { 0x8000 } else { 0 };
                assert_eq!((f.from_rm)(&sx(), rm), Ok(s | lo), "exact {i} {rm:?}");
            }
            continue;
        }
        let mid = upper_val(f, lo + 1)
            .add(&(f.to_exact)(lo, P), P, RoundingMode::ToEven)
            .ldexp(-1, P, RoundingMode::None);
        let cmp = i32::try_from(x.compare(&mid).unwrap_or(0))
            .unwrap_or(0)
            .signum();
        check(f, &sx(), lo, cmp, neg, "search");
    }
}

/// A fresh value equal to `x` (built by an exact zero-shift, not a copy helper).
fn m_copy(x: &ExactNum) -> ExactNum {
    ExactNum::add_ref(x, &ExactNum::from_u64(0, P), P, RoundingMode::None)
}

#[test]
fn ieee16_binary_search_oracle_multi_limb() {
    search_oracle(&H, 0x1616_0000_dead_beef);
}

#[test]
fn bfloat16_binary_search_oracle_multi_limb() {
    search_oracle(&B, 0xbf16_0000_cafe_f00d);
}

fn extremes(f: &Fmt) {
    let minsub = 1u16;
    let inf = f.max + 1;
    for neg in [false, true] {
        let s = if neg { 0x8000u16 } else { 0 };
        let tiny = || if neg { -pow2(-10000) } else { pow2(-10000) };
        let huge = || if neg { -pow2(10000) } else { pow2(10000) };
        for rm in MODES {
            let away_t = match rm {
                RoundingMode::Up => !neg,
                RoundingMode::Down => neg,
                RoundingMode::FromZero => true,
                RoundingMode::ToEven
                | RoundingMode::ToOdd
                | RoundingMode::ToZero
                | RoundingMode::None => false,
            };
            let to_inf = match rm {
                RoundingMode::Up => !neg,
                RoundingMode::Down => neg,
                RoundingMode::FromZero | RoundingMode::ToEven | RoundingMode::ToOdd => true,
                RoundingMode::ToZero | RoundingMode::None => false,
            };
            let (wt, wh) = if rm == RoundingMode::None {
                (Err(Error::InvalidArgument), Err(Error::InvalidArgument))
            } else {
                (
                    Ok(s | if away_t { minsub } else { 0 }),
                    Ok(s | if to_inf { inf } else { f.max }),
                )
            };
            assert_eq!((f.from_rm)(&tiny(), rm), wt, "tiny neg={neg} {rm:?}");
            assert_eq!((f.from_rm)(&huge(), rm), wh, "huge neg={neg} {rm:?}");
            assert_eq!((f.from_rm)(&zenith_float_num::NAN, rm), Ok(f.qnan));
            assert_eq!((f.from_rm)(&zenith_float_num::INF_POS, rm), Ok(inf));
            assert_eq!(
                (f.from_rm)(&zenith_float_num::INF_NEG, rm),
                Ok(0x8000 | inf)
            );
            assert_eq!((f.from_rm)(&-ExactNum::from_u64(0, P), rm), Ok(0x8000));
            assert_eq!((f.from_rm)(&ExactNum::from_u64(0, P), rm), Ok(0));
        }
    }
}

#[test]
fn ieee16_extremes() {
    extremes(&H);
}

#[test]
fn bfloat16_extremes() {
    extremes(&B);
}

#[test]
fn golds_and_constants() {
    let g = |x: ExactNum, k: i64, n: u64| {
        eq(
            &x.ldexp(-k, P, RoundingMode::None),
            &ExactNum::from_u64(n, P),
        )
    };
    assert!(g(Ieee16::from_bits(0x3C00).to_exact(64), 0, 1));
    assert!(g(Ieee16::from_bits(0x7BFF).to_exact(64), 0, 65504));
    assert!(g(Ieee16::from_bits(0x0400).to_exact(64), -14, 1));
    assert!(g(Ieee16::from_bits(0x0001).to_exact(64), -24, 1));
    assert!(g(-Ieee16::from_bits(0xC000).to_exact(64), 0, 2));
    assert!(g(Bfloat16::from_bits(0x3F80).to_exact(64), 0, 1));
    assert!(g(Bfloat16::from_bits(0x7F7F).to_exact(64), 120, 255));
    assert!(g(Bfloat16::from_bits(0x0080).to_exact(64), -126, 1));
    assert!(g(Bfloat16::from_bits(0x0001).to_exact(64), -133, 1));
    assert_eq!(Ieee16::NAN.to_bits(), 0x7E00);
    assert_eq!(Ieee16::INFINITY.to_bits(), 0x7C00);
    assert_eq!(Ieee16::NEG_INFINITY.to_bits(), 0xFC00);
    assert_eq!(Ieee16::NEG_ZERO.to_bits(), 0x8000);
    assert_eq!(Ieee16::ZERO.to_bits(), 0);
    assert_eq!(Ieee16::MAX.to_bits(), 0x7BFF);
    assert_eq!(Ieee16::MIN_POSITIVE.to_bits(), 0x0400);
    assert_eq!(Ieee16::MIN_POSITIVE_SUBNORMAL.to_bits(), 0x0001);
    assert_eq!(Bfloat16::NAN.to_bits(), 0x7FC0);
    assert_eq!(Bfloat16::INFINITY.to_bits(), 0x7F80);
    assert_eq!(Bfloat16::NEG_INFINITY.to_bits(), 0xFF80);
    assert_eq!(Bfloat16::NEG_ZERO.to_bits(), 0x8000);
    assert_eq!(Bfloat16::ZERO.to_bits(), 0);
    assert_eq!(Bfloat16::MAX.to_bits(), 0x7F7F);
    assert_eq!(Bfloat16::MIN_POSITIVE.to_bits(), 0x0080);
    assert_eq!(Bfloat16::MIN_POSITIVE_SUBNORMAL.to_bits(), 0x0001);
    // 1 + 2^-1000 at 1088 bits: the only set bit below the top word sits 1000 bits down.
    let x = || ExactNum::from_u64(1, 1088).add(&pow2(-1000), 1088, RoundingMode::ToEven);
    assert_eq!(
        Ieee16::from_exact_rm(&x(), RoundingMode::ToEven).map(Ieee16::to_bits),
        Ok(0x3C00)
    );
    assert_eq!(
        Ieee16::from_exact_rm(&x(), RoundingMode::Up).map(Ieee16::to_bits),
        Ok(0x3C01)
    );
    assert_eq!(
        Ieee16::from_exact_rm(&x(), RoundingMode::None).map(Ieee16::to_bits),
        Err(Error::InvalidArgument)
    );
    assert_eq!(
        Bfloat16::from_exact_rm(&x(), RoundingMode::FromZero).map(Bfloat16::to_bits),
        Ok(0x3F81)
    );
    assert_eq!(
        Bfloat16::from_exact_rm(&x(), RoundingMode::Down).map(Bfloat16::to_bits),
        Ok(0x3F80)
    );
}
