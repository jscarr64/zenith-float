//! Exact-integer fixed-point `log2` and rational-power helpers for the
//! parameter choices that FLINT/Arb make with C `double`s: Stirling term
//! counts and block lengths, Carlson Taylor orders, and `pfq` term counts.
//!
//! The crate denies `clippy::float_arithmetic` (no hardware floats in the
//! library), so these helpers use integers only. A choice is then the same on
//! every target and in both the `std` and `no_std` builds, and never depends on
//! a libm. A `log2` value is `Q32`: an `i128` `v` stands for `v / 2^32`.
//! `*_down` and `*_up` bracket the true value: `down <= 2^32 log2 x <= up`.
//!
//! These values only choose how much work to do (how many terms, which order,
//! which block length). A returned number never depends on them.

use crate::ExactInt;
use crate::ExactNum;
use core::cmp::Ordering;

/// Fractional bits of a `Q32` value.
pub(crate) const FRAC_BITS: u32 = 32;
/// `1.0` in `Q32`.
pub(crate) const ONE: i128 = 1 << FRAC_BITS;

/// `log2_*_down` is at most this many `Q32` units below the true value (see
/// [`log2_mant_down`]).
const LOG2_SLACK: i128 = 2;

// floor / ceil of `c * 2^32`, from mpmath 1.3.0 at 400 bits:
//   log2(pi)   = 1.651496129472318798...
//   log2(e)    = 1.4426950408889634074...
//   log2(12)   = 3.5849625007211561815...
//   log2(2 pi) = 2.651496129472318798...
//   log2(1e-30)  = -99.657842846620870436...
//   log2(1e-100) = -332.19280948873623479...
/// `floor(2^32 log2 pi)`.
pub(crate) const LOG2_PI_DOWN: i128 = 7_093_121_865;
/// `ceil(2^32 log2 pi)`.
pub(crate) const LOG2_PI_UP: i128 = 7_093_121_866;
/// `floor(2^32 log2 e)`.
pub(crate) const LOG2_E_DOWN: i128 = 6_196_328_018;
/// `floor(2^32 log2 12)`.
pub(crate) const LOG2_12_DOWN: i128 = 15_397_296_697;
/// `floor(2^32 log2 (2 pi))`.
pub(crate) const LOG2_2PI_DOWN: i128 = 11_388_089_161;
/// `floor(2^32 log2 1e-30)`.
pub(crate) const LOG2_1E_M30_DOWN: i128 = -428_027_175_817;
/// `floor(2^32 log2 1e-100)`.
pub(crate) const LOG2_1E_M100_DOWN: i128 = -1_426_757_252_721;

/// `floor(2^32 log2(m / 2^63))` or one less, for `m` in `[2^63, 2^64)`.
///
/// Bit-by-bit: with `y = m / 2^63` in `[1, 2)`, square `y`; the next bit of
/// `log2 y` is 1 exactly when `y^2 >= 2`, and then `y^2 / 2` carries on. `y` is
/// held in `Q63` and each square is truncated, so `y` only shrinks. The bits
/// found are then those of a value no larger than `log2(m / 2^63)`. The
/// truncation errors (below `2^-62` relative per step, weighted by `2^-k` at
/// step `k`) total under one `Q32` unit, so the result is `floor` or
/// `floor - 1`, and the true value is below `result + LOG2_SLACK`.
fn log2_mant_down(m: u64) -> i128 {
    let mut y = u128::from(m);
    let mut r: i128 = 0;
    for i in 1..=FRAC_BITS {
        // `y < 2^64`, so `y * y < 2^128`.
        y = (y * y) >> 63;
        if y >> 64 != 0 {
            y >>= 1;
            r |= 1 << (FRAC_BITS - i);
        }
    }
    r
}

/// `Q32` lower bound on `log2 n`; `None` for `n = 0`. Exact for a power of two.
pub(crate) fn log2_u64_down(n: u64) -> Option<i128> {
    if n == 0 {
        return None;
    }
    let lz = n.leading_zeros();
    Some(i128::from(63 - lz) * ONE + log2_mant_down(n << lz))
}

/// `Q32` upper bound on `log2 n`; `None` for `n = 0`. Exact for a power of two.
pub(crate) fn log2_u64_up(n: u64) -> Option<i128> {
    let down = log2_u64_down(n)?;
    Some(if n.is_power_of_two() { down } else { down + LOG2_SLACK })
}

/// `Q32` lower bound on `log2 |x|`; `None` for zero, Inf and NaN.
pub(crate) fn log2_abs_down(x: &ExactNum) -> Option<i128> {
    let (m, _, e) = x.top_u64_mantissa_and_sticky()?;
    // |x| = (m / 2^63) 2^(e - 1) plus the bits below the top 64.
    Some(i128::from(e).checked_sub(1)?.checked_mul(ONE)? + log2_mant_down(m))
}

/// `Q32` upper bound on `log2 |x|`; `None` for zero, Inf and NaN. Exact for a
/// power of two.
pub(crate) fn log2_abs_up(x: &ExactNum) -> Option<i128> {
    let (m, sticky, _) = x.top_u64_mantissa_and_sticky()?;
    let down = log2_abs_down(x)?;
    // The bits below the top 64 add less than 2^-63 to log2.
    Some(if !sticky && m.is_power_of_two() { down } else { down + LOG2_SLACK })
}

/// `v / 2^32` rounded to nearest, ties up (C `round` for the values used here,
/// which are never negative ties).
pub(crate) fn round_q(v: i128) -> i128 {
    v.saturating_add(ONE >> 1) >> FRAC_BITS
}

/// The largest `k` in `[0, hi]` with `k <= (cn / cd) (xn / xd)^(a / b)`,
/// decided exactly by `(cd k)^b xd^a <= cn^b xn^a` (all exact integers), so it is
/// the floor of the exact real value, capped at `hi`. `b >= 1`.
///
/// FLINT computes these as `double` expressions such as `2.1 * pow(prec, 0.4)`
/// and truncates. The exact floor differs from that only when the real value is
/// within a few `double` ulps of an integer.
pub(crate) fn floor_scaled_root(
    cn: u64,
    cd: u64,
    xn: u64,
    xd: u64,
    a: u64,
    b: u64,
    hi: u64,
) -> u64 {
    let rhs = ExactInt::from_u64(cn)
        .pow(b)
        .mul(&ExactInt::from_u64(xn).pow(a));
    let xd_a = ExactInt::from_u64(xd).pow(a);
    let fits = |k: u64| {
        ExactInt::from_u64(cd)
            .mul_u64(k)
            .pow(b)
            .mul(&xd_a)
            .soft_cmp(&rhs)
            != Ordering::Greater
    };
    if fits(hi) {
        return hi;
    }
    // fits(lo) holds, fits(top) does not (fits(0) holds: 0 <= rhs).
    let (mut lo, mut top) = (0u64, hi);
    while top - lo > 1 {
        let mid = lo + (top - lo) / 2;
        if fits(mid) {
            lo = mid;
        } else {
            top = mid;
        }
    }
    lo
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::RoundingMode;
    use core::ops::Neg;

    #[test]
    fn log2_mant_brackets_and_is_exact_on_powers_of_two() -> Result<(), &'static str> {
        // Golds: floor(2^32 log2 n), mpmath 1.3.0 at 400 bits (E/b7/fixlog_golds.out).
        for &(n, fl) in GOLDS {
            let d = log2_u64_down(n).ok_or("down None")?;
            let u = log2_u64_up(n).ok_or("up None")?;
            assert!(
                d <= fl && fl < u.max(d + 1),
                "n={n} down={d} floor={fl} up={u}"
            );
            assert!(
                fl - d <= 1,
                "n={n} down={d} more than one unit below floor={fl}"
            );
            if n.is_power_of_two() {
                assert_eq!(d, u);
                assert_eq!(d, i128::from(n.trailing_zeros()) * ONE);
            }
        }
        assert_eq!(log2_u64_down(0), None);
        assert_eq!(log2_u64_up(0), None);
        Ok(())
    }

    #[test]
    fn log2_abs_matches_u64_and_handles_specials() -> Result<(), &'static str> {
        for &(n, _) in GOLDS {
            let x = ExactNum::from_word(n as crate::Word, 128);
            assert_eq!(log2_abs_down(&x), log2_u64_down(n), "n={n}");
            assert_eq!(log2_abs_down(&x.neg()), log2_u64_down(n), "-n={n}");
        }
        // 3 * 2^-200 and its exponent shift.
        let x = ExactNum::from_word(3, 64).ldexp(-200, 64, RoundingMode::None);
        assert_eq!(log2_abs_down(&x), log2_u64_down(3).map(|v| v - 200 * ONE));
        // A value with bits below the top 64: 2^64 + 1 at 128 bits.
        let y = ExactNum::from_word(1, 128)
            .ldexp(64, 128, RoundingMode::None)
            .add(&ExactNum::from_word(1, 128), 128, RoundingMode::None);
        let (d, u) = (
            log2_abs_down(&y).ok_or("down None")?,
            log2_abs_up(&y).ok_or("up None")?,
        );
        assert!(d <= 64 * ONE && 64 * ONE < u, "2^64+1: {d} {u}");
        assert_eq!(log2_abs_down(&ExactNum::new(64)), None);
        assert_eq!(log2_abs_down(&ExactNum::nan(None)), None);
        assert_eq!(log2_abs_down(&crate::INF_POS), None);
        Ok(())
    }

    #[test]
    fn round_q_ties_up() {
        assert_eq!(round_q(5 * ONE / 2), 3);
        assert_eq!(round_q(5 * ONE / 2 - 1), 2);
        assert_eq!(round_q(-(ONE / 4)), 0);
    }

    #[test]
    fn floor_scaled_root_is_the_exact_floor() {
        // floor(2.1 p^0.4): p = 2^20 gives 2.1 * 2^8 = 537.6 -> 537 (cap 500 -> 500).
        assert_eq!(floor_scaled_root(21, 10, 1 << 20, 1, 2, 5, 10_000), 537);
        assert_eq!(floor_scaled_root(21, 10, 1 << 20, 1, 2, 5, 500), 500);
        // 2 * 32^0.4 = 2 * 4 = 8 exactly: the floor is 8, not 7.
        assert_eq!(floor_scaled_root(2, 1, 32, 1, 2, 5, 1000), 8);
        // x = 0 gives 0.
        assert_eq!(floor_scaled_root(27, 100, 0, 2, 2, 5, 56), 0);
        // Brute force against the defining inequality on a grid.
        for p in 1u64..300 {
            let k = floor_scaled_root(23, 10, p, 1, 17, 50, 1000);
            let le = |k: u64| {
                ExactInt::from_u64(10 * k).pow(50).soft_cmp(
                    &ExactInt::from_u64(23)
                        .pow(50)
                        .mul(&ExactInt::from_u64(p).pow(17)),
                ) != Ordering::Greater
            };
            assert!(le(k) && !le(k + 1), "p={p} k={k}");
        }
    }

    // (n, floor(2^32 log2 n)), mpmath 1.3.0 at 400 bits (E/b7/fixlog_golds.py).
    const GOLDS: &[(u64, i128)] = &[
        (1, 0),
        (2, 4294967296),
        (3, 6807362105),
        (5, 9972605231),
        (7, 12057497578),
        (10, 14267572527),
        (12, 15397296697),
        (100, 28535145054),
        (255, 34335486563),
        (256, 34359738368),
        (257, 34383895622),
        (1000, 42802717581),
        (4095, 51538094591),
        (4097, 51541120142),
        (65535, 68719382186),
        (65537, 68719571283),
        (1000000, 85605435163),
        (2147483647, 133143986173),
        (4294967297, 137438953473),
        (12157665459056928801, 272294484239),
        (1000000000000000000, 256816305489),
        (4611686018427387905, 266287972352),
        (9223372036854775807, 270582939647),
        (9223372036854775808, 270582939648),
        (9223372036854775809, 270582939648),
        (18446744073709551615, 274877906943),
        (3735928559, 136574890753),
        (11400714819323198485, 271896160629),
    ];
}
