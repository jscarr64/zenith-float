//! `leading_bit_exponent`, `top_u64_mantissa_and_sticky`, `from_mag_parts` and `pow2` (the prerequisites P1 to P3 of
//! zenith-bound) over the whole exponent range, subnormals included, and widths from one word to 313 words
//! (20 000 bits). The oracle is exact integer arithmetic (`rug::Integer`) on the words of the operand, and MPFR
//! (`Float`) for the powers of two that MPFR can hold.

use crate::mpfr::retry_budget::use_full_exponent_range;
use rug::integer::Order;
use rug::{Float, Integer};
use zenith_float_num::{
    Error, ExactNum, RoundingMode, Sign, EXPONENT_MAX, EXPONENT_MIN, INF_NEG, INF_POS,
};

type Res = Result<(), Box<dyn std::error::Error>>;

const MIN: i64 = EXPONENT_MIN as i64;
const MAX: i64 = EXPONENT_MAX as i64;

/// A deterministic word source (a 64-bit linear congruential step).
fn next_word(state: &mut u64) -> u64 {
    *state = state
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    *state ^ (*state >> 29)
}

/// `len` words (least significant first) whose integer has exactly `bits` significant bits.
fn words_with_bits(len: usize, bits: usize, seed: u64) -> Vec<u64> {
    let mut st = seed;
    let mut w: Vec<u64> = (0..len).map(|_| next_word(&mut st)).collect();
    // clear everything above `bits`, set bit `bits - 1`
    let top = (bits - 1) / 64;
    let off = (bits - 1) % 64;
    for (i, x) in w.iter_mut().enumerate() {
        if i > top {
            *x = 0;
        } else if i == top {
            let mask = if off == 63 { u64::MAX } else { (1u64 << (off + 1)) - 1 };
            *x = (*x & mask) | (1u64 << off);
        }
    }
    w
}

/// What the oracle expects of `(words, exponent)`: the leading-bit exponent, the top 64 bits and the sticky bit, and the
/// `0.m * 2^e` exponent of the normalized value.
struct Expect {
    leading: i64,
    top: u64,
    sticky: bool,
}

fn expect_of(words: &[u64], e: i64) -> Expect {
    let m = Integer::from_digits(words, Order::Lsf);
    let bl = m.significant_bits() as i64;
    let lz = 64 * words.len() as i64 - bl;
    let leading = e - 1 - lz;
    if bl >= 64 {
        let sh = (bl - 64) as u32;
        let t = Integer::from(&m >> sh);
        let back = Integer::from(&t << sh);
        let low = Integer::from(&m - &back);
        Expect {
            leading,
            top: t.to_u64_wrapping(),
            sticky: low != 0,
        }
    } else {
        let t = Integer::from(&m << (64 - bl) as u32);
        Expect {
            leading,
            top: t.to_u64_wrapping(),
            sticky: false,
        }
    }
}

/// An operand from words and an exponent: normal when the top bit of the top word is set, subnormal otherwise.
fn operand(words: &[u64], negative: bool, e: i64) -> Result<ExactNum, Box<dyn std::error::Error>> {
    let m = Integer::from_digits(words, Order::Lsf);
    let n = if words.last().is_some_and(|w| w >> 63 == 1) {
        64 * words.len()
    } else {
        m.significant_bits() as usize
    };
    let s = if negative { Sign::Neg } else { Sign::Pos };
    let x = ExactNum::from_raw_parts(words, n, s, e as i32, false);
    if x.is_nan() {
        return Err(format!("from_raw_parts refused {words:x?} n={n} e={e}").into());
    }
    Ok(x)
}

/// Operands: normal and subnormal, 1 to 313 words, exponents over the whole range, both signs.
fn operands() -> Result<Vec<(ExactNum, Expect, String)>, Box<dyn std::error::Error>> {
    let mut v = Vec::new();
    let es = [MIN, MIN + 1, -1, 0, 1, 77, MAX - 1, MAX];
    let mut seed = 1u64;
    for &len in &[1usize, 2, 3, 5, 313] {
        for &e in &es {
            for neg in [false, true] {
                seed += 1;
                // normal: full width, power of two, all ones, random
                let normals = [
                    words_with_bits(len, 64 * len, seed),
                    {
                        let mut w = vec![0u64; len];
                        if let Some(t) = w.last_mut() {
                            *t = 1 << 63;
                        }
                        w
                    },
                    vec![u64::MAX; len],
                ];
                for w in &normals {
                    let x = operand(w, neg, e)?;
                    v.push((
                        x,
                        expect_of(w, e),
                        format!("normal len {len} e {e} neg {neg}"),
                    ));
                }
            }
        }
        // subnormal: exponent MIN, `lz` leading zero bits
        let max_lz = 64 * len as i64 - 1;
        for lz in [1i64, 2, 63, 64, 65, 100, 127, 128, 129, 200, max_lz] {
            if lz > max_lz {
                continue;
            }
            for neg in [false, true] {
                seed += 1;
                let bits = (64 * len as i64 - lz) as usize;
                let w = words_with_bits(len, bits, seed);
                let x = operand(&w, neg, MIN)?;
                v.push((
                    x,
                    expect_of(&w, MIN),
                    format!("subnormal len {len} lz {lz} neg {neg}"),
                ));
            }
        }
    }
    Ok(v)
}

#[test]
fn leading_bit_exponent_matches_the_integer_model() -> Res {
    for (x, want, what) in operands()? {
        assert_eq!(x.leading_bit_exponent(), Some(want.leading), "{what}");
    }
    // the named values
    for p in [1usize, 24, 64, 65, 128, 1000, 20_000] {
        let pw = (p.div_ceil(64) * 64) as i64;
        assert_eq!(
            ExactNum::min_positive(p).leading_bit_exponent(),
            Some(MIN - pw),
            "min_positive({p})"
        );
        assert_eq!(
            ExactNum::min_positive_normal(p).leading_bit_exponent(),
            Some(MIN - 1),
            "min_positive_normal({p})"
        );
        assert_eq!(
            ExactNum::max_value(p).leading_bit_exponent(),
            Some(MAX - 1),
            "max_value({p})"
        );
        assert_eq!(
            ExactNum::min_value(p).leading_bit_exponent(),
            Some(MAX - 1),
            "min_value({p})"
        );
    }
    // one is 2^0; a subnormal's value is below the exponent field: 2^-2147483718 is MIN - 70
    assert_eq!(ExactNum::from_word(1, 64).leading_bit_exponent(), Some(0));
    assert_eq!(ExactNum::from_word(3, 64).leading_bit_exponent(), Some(1));
    let tiny = ExactNum::from_word(1, 192)
        .ldexp(EXPONENT_MIN, 192, RoundingMode::ToEven)
        .ldexp(-70, 192, RoundingMode::ToEven);
    assert_eq!(tiny.leading_bit_exponent(), Some(MIN - 70));
    // not a finite nonzero number
    assert_eq!(ExactNum::new(64).leading_bit_exponent(), None);
    assert_eq!(INF_POS.leading_bit_exponent(), None);
    assert_eq!(ExactNum::nan(None).leading_bit_exponent(), None);
    Ok(())
}

#[test]
fn leading_bit_exponent_agrees_with_mpfr() -> Res {
    use_full_exponent_range();
    // values m * 2^k that MPFR holds exactly: its exponent is the leading-bit exponent + 1
    for &m in &[1u64, 3, 0xdead_beef_1234_5677, u64::MAX] {
        for &k in &[-2_000_000_000i64, -1000, -1, 0, 5, 1000, 2_000_000_000] {
            let f = Float::with_val(64, m) << (k as i32);
            let Some(e) = f.get_exp() else {
                return Err("no exponent".into());
            };
            let x = ExactNum::from_word(m, 64).ldexp(k as i32, 64, RoundingMode::ToEven);
            assert_eq!(
                x.leading_bit_exponent(),
                Some(i64::from(e) - 1),
                "{m:#x} 2^{k}"
            );
        }
    }
    Ok(())
}

#[test]
fn top_u64_mantissa_and_sticky_matches_the_integer_model() -> Res {
    let mut sticky_seen = 0;
    for (x, want, what) in operands()? {
        let got = x.top_u64_mantissa_and_sticky();
        assert_eq!(
            got,
            Some((want.top, want.sticky, want.leading + 1)),
            "{what}"
        );
        if want.sticky {
            sticky_seen += 1;
        }
        // the top bit is always set
        if let Some((t, _, _)) = got {
            assert_eq!(t >> 63, 1, "{what}");
        }
    }
    assert!(sticky_seen > 100, "the sweep must contain sticky operands");
    // a single low bit far below the top word makes the sticky bit
    let mut w = vec![0u64; 313];
    if let Some(t) = w.last_mut() {
        *t = 1 << 63;
    }
    if let Some(b) = w.first_mut() {
        *b = 1;
    }
    let x = operand(&w, false, 5)?;
    assert_eq!(x.top_u64_mantissa_and_sticky(), Some((1 << 63, true, 5)));
    // the sign does not matter; zero, Inf and NaN have no parts
    let y = operand(&w, true, 5)?;
    assert_eq!(y.top_u64_mantissa_and_sticky(), Some((1 << 63, true, 5)));
    assert_eq!(ExactNum::new(64).top_u64_mantissa_and_sticky(), None);
    assert_eq!(INF_NEG.top_u64_mantissa_and_sticky(), None);
    assert_eq!(ExactNum::nan(None).top_u64_mantissa_and_sticky(), None);
    // 1, 3 and the smallest subnormal
    assert_eq!(
        ExactNum::from_word(1, 64).top_u64_mantissa_and_sticky(),
        Some((1 << 63, false, 1))
    );
    assert_eq!(
        ExactNum::from_word(3, 64).top_u64_mantissa_and_sticky(),
        Some((3 << 62, false, 2))
    );
    assert_eq!(
        ExactNum::min_positive(20_000).top_u64_mantissa_and_sticky(),
        Some((1 << 63, false, MIN - 20_032 + 1))
    );
    Ok(())
}

/// The value of an `ExactNum` as `(odd integer, exponent of its lowest set bit)`, from its raw words.
fn odd_and_exp(x: &ExactNum) -> Result<(Integer, i64), Box<dyn std::error::Error>> {
    let Some((words, _, _, e, _)) = x.as_raw_parts() else {
        return Err("no raw parts".into());
    };
    let m = Integer::from_digits(words, Order::Lsf);
    let Some(tz) = m.find_one(0) else {
        return Err("zero mantissa".into());
    };
    let odd = Integer::from(&m >> tz);
    Ok((odd, i64::from(e) - 64 * words.len() as i64 + i64::from(tz)))
}

/// The expected outcome of `from_mag_parts(mant, exp, p)`: `None` when it must be an error.
fn mag_expect(mant: u64, exp: i64, p: usize) -> Option<(Integer, i64)> {
    if mant == 0 {
        return Some((Integer::new(), 0));
    }
    if exp > MAX {
        return None;
    }
    let tz = i128::from(mant.trailing_zeros());
    let lowest = i128::from(exp) - 64 + tz; // exponent of the lowest set bit
    let grid = i128::from(MIN) - (p.div_ceil(64) * 64) as i128; // exponent of the smallest subnormal
    if lowest < grid {
        return None;
    }
    Some((Integer::from(mant >> tz), lowest as i64))
}

#[test]
fn from_mag_parts_is_exact_or_an_error() -> Res {
    let mants = [1u64 << 63, u64::MAX, 1, 3, 0xdead_beef_1234_5677, (1 << 63) | 1, 0];
    let exps = [
        i64::MAX,
        i64::MIN,
        MAX + 1,
        MAX,
        MAX - 1,
        100,
        1,
        0,
        -5,
        MIN + 64,
        MIN + 65,
        MIN + 1,
        MIN,
        MIN - 1,
        MIN - 63,
        MIN - 64,
        MIN - 65,
        MIN - 127,
        MIN - 128,
        MIN - 129,
        MIN - 200,
        MIN - 20_000,
        MIN - 20_063,
    ];
    let mut ok = 0;
    let mut err = 0;
    for &p in &[1usize, 24, 64, 65, 128, 200, 1000, 20_000] {
        for &mant in &mants {
            for &exp in &exps {
                let got = ExactNum::from_mag_parts(mant, exp, p);
                match (mag_expect(mant, exp, p), got) {
                    (Some((odd, e)), Ok(x)) => {
                        ok += 1;
                        if mant == 0 {
                            assert!(
                                x.is_zero() && !x.is_negative(),
                                "zero mant, exp {exp} p {p}"
                            );
                            assert!(!x.inexact(), "exact");
                            continue;
                        }
                        let (go, ge) = odd_and_exp(&x)?;
                        assert_eq!((go, ge), (odd, e), "mant {mant:#x} exp {exp} p {p}");
                        assert!(!x.inexact(), "mant {mant:#x} exp {exp} p {p}: inexact");
                        assert!(x.is_positive(), "positive");
                        // the stored width holds p bits
                        assert!(
                            x.mantissa_max_bit_len().is_some_and(|b| b >= p),
                            "width at p {p}"
                        );
                    }
                    (None, Err(_)) => err += 1,
                    (want, got) => {
                        return Err(format!(
                            "mant {mant:#x} exp {exp} p {p}: expected {} got {got:?}",
                            if want.is_some() { "Ok" } else { "an error" }
                        )
                        .into());
                    }
                }
            }
        }
    }
    assert!(ok > 300 && err > 100, "ok {ok} err {err}");
    // the overflow error carries the sign of a positive overflow
    assert_eq!(
        ExactNum::from_mag_parts(1 << 63, MAX + 1, 64).err(),
        Some(Error::ExponentOverflow(Sign::Pos))
    );
    Ok(())
}

#[test]
fn pow2_is_exact_over_the_range() -> Res {
    for &p in &[1usize, 24, 64, 65, 128, 1000, 20_000] {
        let pw = (p.div_ceil(64) * 64) as i64;
        for &exp in &[
            i64::MAX,
            MAX + 1,
            MAX,
            MAX - 1,
            1000,
            63,
            64,
            1,
            0,
            -1,
            -63,
            -64,
            -65,
            -1000,
            MIN + 1,
            MIN,
            MIN - 1,
            MIN - 70,
            MIN - pw + 1,
            MIN - pw,
            MIN - pw - 1,
            i64::MIN,
        ] {
            let got = ExactNum::pow2(exp, p);
            // 2^exp has its lowest (only) bit at exp; it fits iff exp <= MAX - 1 and exp >= grid
            let fits = exp < MAX && exp >= MIN - pw;
            match got {
                Ok(x) => {
                    assert!(fits, "pow2({exp}, {p}) must be an error");
                    assert_eq!(x.leading_bit_exponent(), Some(exp), "pow2({exp}, {p})");
                    let (odd, e) = odd_and_exp(&x)?;
                    assert_eq!((odd, e), (Integer::from(1), exp));
                    assert!(!x.inexact());
                    // MPFR holds it when the exponent is in its default range
                    if exp.abs() < 2_000_000_000 {
                        use_full_exponent_range();
                        let f = Float::with_val(2, 1) << (exp as i32);
                        assert_eq!(f.get_exp(), Some(exp as i32 + 1));
                    }
                }
                Err(_) => assert!(!fits, "pow2({exp}, {p}) must be Ok"),
            }
        }
    }
    Ok(())
}
