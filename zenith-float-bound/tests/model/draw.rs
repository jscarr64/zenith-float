//! Generators over the WHOLE exponent range: about 30 percent of exponents are within 3 ticks of i64::MIN or i64::MAX, and mantissas follow bit patterns (all ones, single bit, alternating, top bit only).

use super::dyadic::MagV;
use super::rng::SplitMix64;

const TOP: u64 = 1 << 63;

/// Generate a random exponent.
pub fn draw_exp(rng: &mut SplitMix64) -> i64 {
    let r = rng.below(100);
    if r < 15 {
        i64::MIN + i64::from(rng.below(4))
    } else if r < 30 {
        i64::MAX - i64::from(rng.below(4))
    } else if r < 45 {
        rng.next_i64() % 401
    } else {
        rng.next_i64()
    }
}

/// Generate a random mantissa.
pub fn draw_mant(rng: &mut SplitMix64) -> u64 {
    let choice = rng.below(8);
    match choice {
        0 => u64::MAX,
        1 => TOP,
        2 => 0xAAAA_AAAA_AAAA_AAAA,
        3 => 0xD555_5555_5555_5555,
        4 => TOP + 1,
        5 => u64::MAX - 1,
        6 => TOP | (1u64 << rng.below(63)),
        _ => rng.next_u64() | TOP,
    }
}

/// Generate a random finite magnitude.
pub fn draw_finite(rng: &mut SplitMix64) -> MagV {
    let exp = draw_exp(rng);
    let mant = draw_mant(rng);
    MagV::Finite { exp, mant }
}

/// Generate a random magnitude (including infinity).
pub fn draw_mag(rng: &mut SplitMix64) -> MagV {
    let r = rng.below(100);
    if r < 4 {
        MagV::Finite { exp: 0, mant: 0 }
    } else if r < 7 {
        MagV::Infinity
    } else {
        draw_finite(rng)
    }
}

/// Generate a partner magnitude for alignment tests.
pub fn draw_partner(rng: &mut SplitMix64, a: &MagV) -> MagV {
    const GAPS: [i64; 12] = [0, 1, 2, 62, 63, 64, 65, 66, 127, 128, 129, 200];
    if let MagV::Finite { exp, mant: a_mant } = a {
        if *a_mant != 0 && rng.below(2) == 0 {
            let gap = match GAPS.get(usize::from(rng.below(12))) {
                Some(&g) => g,
                None => 0,
            };
            let sign = rng.below(2);
            let shifted = if sign == 0 {
                exp.checked_add(gap)
            } else {
                exp.checked_sub(gap)
            };
            if let Some(e) = shifted {
                let m = draw_mant(rng);
                return MagV::Finite { exp: e, mant: m };
            }
        }
    }
    draw_mag(rng)
}

/// Generate a raw (possibly non-canonical) mantissa and exponent.
pub fn draw_raw(rng: &mut SplitMix64) -> (u64, i64) {
    let mant = match rng.below(8) {
        0 => 0,
        1 => 1,
        2 => u64::MAX,
        3 => 1u64 << rng.below(64),
        4 => rng.next_u64() >> rng.below(64),
        5 => TOP | rng.next_u64(),
        6 => (1u64 << rng.below(64)) - 1,
        _ => rng.next_u64(),
    };
    (mant, draw_exp(rng))
}
