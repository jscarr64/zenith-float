//! Normalise a positive dyadic integer to a canonical `Mag`.
//!
//! `norm(m, scale, sticky, up)` rounds the real value `m * 2^scale`, plus an unresolved
//! positive residue when `sticky` is set. `m` is a `NonZeroU128`, so the state "no bits at
//! all" cannot be expressed. The `Mag` exponent field of a 64-bit mantissa `mant` is
//! `scale + 64` after the significand has been aligned to `[2^63, 2^64)`.
//!
//! Exponent arithmetic uses `checked_*` only (CONTRACT 3.1). A step that cannot fail for any
//! `i64` exponent says so in a comment; the `None` it carries is settled by [`settle`] with the
//! one sound answer: `Infinity` on the way up, zero on the way down.

use super::Mag;
use super::TOP;
use core::num::NonZeroU128;

/// The sound stand-in when no exact result can be formed: `Infinity` on the way up, zero on
/// the way down. It is also the specified answer for an `Infinity` operand or a zero divisor
/// (CONTRACT 14.6).
pub(super) fn loose(up: bool) -> Mag {
    if up {
        Mag::infinity()
    } else {
        Mag::zero()
    }
}

/// `r`, or [`loose`] when a checked step returned `None`. Every operation ends here once.
pub(super) fn settle(r: Option<Mag>, up: bool) -> Mag {
    match r {
        Some(m) => m,
        None => loose(up),
    }
}

/// Low 64 bits of `x`. Callers pass a value below `2^64`, so nothing is lost.
fn low64(x: u128) -> u64 {
    let [_, _, _, _, _, _, _, _, a, b, c, d, e, f, g, h] = x.to_be_bytes();
    u64::from_be_bytes([a, b, c, d, e, f, g, h])
}

/// Place the exponent field `exp` of a mantissa whose top bit is set, clamping to the `Mag`
/// range (CONTRACT 3.2 to 3.4).
fn pack(exp: i128, mant: u64, up: bool) -> Mag {
    match i64::try_from(exp) {
        Ok(e) => Mag::from_canonical(e, mant),
        Err(_) => {
            if exp > i128::from(i64::MAX) {
                if up {
                    Mag::infinity()
                } else {
                    Mag::max_finite()
                }
            } else if up {
                Mag::min_pos()
            } else {
                Mag::zero()
            }
        }
    }
}

/// Round `m * 2^scale` toward `+∞` when `up` is set, otherwise toward zero.
pub(super) fn norm(m: NonZeroU128, scale: i128, sticky: bool, up: bool) -> Option<Mag> {
    // `lz` is in 0..=127. A value with more than 64 bits drops `drop` low bits (at most 64);
    // a shorter one is shifted left by `add` (at most 63). Both arms of the match are reachable.
    let lz = m.leading_zeros();
    let (drop, add) = match 64u32.checked_sub(lz) {
        Some(drop) => (drop, 0),
        None => (0, lz.abs_diff(64)),
    };
    let m = m.get();
    let lost = (m & !(u128::MAX << drop)) != 0;
    // The aligned significand has exactly 64 bits, so its top bit is set and `low64` is exact.
    let mant = low64((m >> drop) << add);
    // |scale| < 2^66 for every caller and drop, add <= 64, so these checked steps cannot fail.
    let exp = scale
        .checked_add(i128::from(drop))?
        .checked_sub(i128::from(add))?
        .checked_add(64)?;
    if up && (sticky || lost) {
        match mant.checked_add(1) {
            Some(next) => Some(pack(exp, next, up)),
            // mant was u64::MAX: the carry gives 2^63 at the next exponent (reachable).
            None => Some(pack(exp.checked_add(1)?, TOP, up)),
        }
    } else {
        Some(pack(exp, mant, up))
    }
}
