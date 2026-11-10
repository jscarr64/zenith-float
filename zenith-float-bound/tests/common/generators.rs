//! S2 generators (TEST_PLAN section 7): exponent draws over the whole range, mantissa patterns, controlled gaps, precision draws,
//! ball descriptors and the six rounding modes. Integer-only, deterministic, no dependency on zenith-float or on the library.
//!
//! FROZEN DRAW ORDER. Every draw consumes a fixed number of generator outputs per branch (documented at each function); changing the
//! tables, the order or the consumption changes every case of every suite and bumps `seeds::GENERATOR_VERSION` (the self-tests pin counts
//! and digests for the fixed seed so drift is caught).
//!
//! The builder's draws in tests/model/draw.rs (draw_exp: 15% near i64::MIN, 15% near i64::MAX, 15% in +-200, 55% uniform; draw_mant: eight
//! patterns; draw_partner gaps 0,1,2,62..66,127..129,200) are a subset of what is built here; they can move onto this module because
//! `draw_exponent` covers the same classes (edge ticks, band around 0, uniform) and the gap rule covers the same gap values at p = 64, 65, 66.
//!
//! Exponents come out as raw `i64`. zenith-float 2.0 mid exponents are `i64` too (Jeff, 2026-10-23 12:36 PM CT), so MID and MAG share the
//! whole `i64` range; the one-tick-outside values of a narrower domain cannot exist there and are left out (`edge_values`, `gap_pair`).
//!
//! Plan steps that need real balls or chains (shrinking radii to 0 then 1 ulp, deleting chain steps) are TODO for later items (L0b, chains);
//! this module only describes the radius classes.

use super::rng::SplitMix64;

/// Gap guard bits (PLAN C6, TEST_PLAN section 7): the exact-add rule is bounded at `p + GUARD`.
pub const GUARD: u32 = 64;
/// Largest precision of the schedule.
pub const MAX_PREC: u32 = 4096;
/// The top bit of a `Mag` mantissa.
pub const TOP: u64 = 1 << 63;

/// zenith-float 2.0 `Exponent` range (`Exponent = i64`, Jeff 2026-10-23 12:36 PM CT; zf 1.x was `i32`): Rusty set it to +-(2^62 - 1), so
/// the sum or difference of two in-range exponents fits in i64. Defined ONCE here; every zb test uses these. Bound to zf's i64
/// `EXPONENT_MIN` / `EXPONENT_MAX` (±(2^62−1), pinned 829f112).
pub const ZF_EXPONENT_MAX: i64 = zenith_float_num::EXPONENT_MAX;
pub const ZF_EXPONENT_MIN: i64 = zenith_float_num::EXPONENT_MIN;
/// Cosmology-scale exponents (2^8 ..= 2^20; Jeff 2026-10-23 12:20 PM CT; 405 is cosmology's 2^+-405), forced into every run, not only
/// reached through the Uniform band. 405 and 1000 are not within one of a power of two, so the PowerOfTwo band does not hit them either.
pub const COSMOLOGY_EXPONENTS: [i64; 8] =
    [256, 405, 1000, 1024, 4097, 65_536, 1_048_575, 1_048_576];

/// The forced exponents, in this fixed order: the four zf bounds `MIN`, `MIN + 1`, `MAX - 1`, `MAX`, then each of `COSMOLOGY_EXPONENTS`
/// followed by its negation (20 values, all inside both `ExpDomain::MAG` and `ExpDomain::MID`). Fuzz suites run them after the random cases.
pub fn forced_exponents() -> Vec<i64> {
    let mut out = vec![
        ZF_EXPONENT_MIN,
        ZF_EXPONENT_MIN + 1,
        ZF_EXPONENT_MAX - 1,
        ZF_EXPONENT_MAX,
    ];
    for e in COSMOLOGY_EXPONENTS {
        out.push(e);
        out.push(-e);
    }
    out
}

/// The forced exponents of every L0a-5 Mag property (tests/harness/fuzz_mag_props.rs), in this fixed order: the four zf bounds `MIN`,
/// `MIN + 1`, `MAX - 1`, `MAX`, then `0`, then each `COSMOLOGY_EXPONENTS` value up to the 1024 scale (256, 405, 1000, 1024) followed by
/// its negation (13 values).
pub fn forced_prop_exponents() -> Vec<i64> {
    let mut out = vec![
        ZF_EXPONENT_MIN,
        ZF_EXPONENT_MIN + 1,
        ZF_EXPONENT_MAX - 1,
        ZF_EXPONENT_MAX,
        0,
    ];
    for e in COSMOLOGY_EXPONENTS {
        if e <= 1024 {
            out.push(e);
            out.push(-e);
        }
    }
    out
}

/// PLACEHOLDER for zf's `MIN_NORMAL_EXPONENT` until the zf pin exists (`seeds::ZF_PIN`): `EXPONENT_MIN + 4096`, the lowest exponent at
/// which a 4096-bit precision still fits above `EXPONENT_MIN`. Replace with zf's constant (it changes the edge table, so bump the generator version).
pub const ZF_MIN_NORMAL_EXPONENT: i64 = ZF_EXPONENT_MIN + 4096;

/// The precision schedule PS (41 values, TEST_PLAN section 9).
pub const PS: [u32; 41] = [
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 20, 24, 31, 32, 33, 48, 53, 63, 64, 65,
    96, 113, 127, 128, 129, 192, 256, 257, 384, 512, 1000, 1024, 2048, 4095, 4096,
];
/// The fast subset PS-fast (18 values) run in the gate.
pub const PS_FAST: [u32; 18] = [
    1, 2, 3, 4, 8, 24, 31, 32, 33, 53, 63, 64, 65, 127, 128, 129, 256, 1024,
];

// ---- exponents ------------------------------------------------------------------------------------------------------------

/// The range exponents are drawn from, with an optional subnormal boundary.
pub struct ExpDomain {
    pub min: i64,
    pub max: i64,
    pub subnormal: Option<i64>,
}

impl ExpDomain {
    /// `Mag` exponents: the whole `i64` range, no subnormal boundary.
    pub const MAG: ExpDomain = ExpDomain {
        min: i64::MIN,
        max: i64::MAX,
        subnormal: None,
    };
    /// zenith-float mid exponents (see `ZF_*`).
    pub const MID: ExpDomain = ExpDomain {
        min: ZF_EXPONENT_MIN,
        max: ZF_EXPONENT_MAX,
        subnormal: Some(ZF_MIN_NORMAL_EXPONENT),
    };

    /// True when `v` is inside the domain or exactly one tick outside it (`min - 1`, `max + 1`): the raw values a draw may return.
    pub fn in_raw_range(&self, v: i64) -> bool {
        self.min.checked_sub(1).is_none_or(|lo| v >= lo)
            && self.max.checked_add(1).is_none_or(|hi| v <= hi)
    }

    /// The edge table, in this fixed order: `min-1 .. min+3`, `max+1 .. max-3`, then `subnormal-1 .. subnormal+1`.
    /// Entries not representable as `i64` (`i64::MIN - 1`, `i64::MAX + 1`) are left out.
    pub fn edge_values(&self) -> Vec<i64> {
        let mut out: Vec<i64> = Vec::new();
        for d in [-1i64, 0, 1, 2, 3] {
            if let Some(v) = self.min.checked_add(d) {
                out.push(v);
            }
        }
        for d in [1i64, 0, -1, -2, -3] {
            if let Some(v) = self.max.checked_add(d) {
                out.push(v);
            }
        }
        if let Some(s) = self.subnormal {
            for d in [-1i64, 0, 1] {
                if let Some(v) = s.checked_add(d) {
                    out.push(v);
                }
            }
        }
        out
    }
}

/// The bands of an exponent draw (TEST_PLAN section 7).
#[derive(Debug, PartialEq, Eq)]
pub enum ExpBand {
    Edge,
    Near0,
    Uniform,
    PowerOfTwo,
    Gap,
}

impl ExpBand {
    pub fn index(&self) -> usize {
        match self {
            ExpBand::Edge => 0,
            ExpBand::Near0 => 1,
            ExpBand::Uniform => 2,
            ExpBand::PowerOfTwo => 3,
            ExpBand::Gap => 4,
        }
    }

    pub fn from_index(i: usize) -> Option<ExpBand> {
        match i {
            0 => Some(ExpBand::Edge),
            1 => Some(ExpBand::Near0),
            2 => Some(ExpBand::Uniform),
            3 => Some(ExpBand::PowerOfTwo),
            4 => Some(ExpBand::Gap),
            _ => None,
        }
    }
}

/// The mix as a constant table, percent per band, in this order. Section 7: 30% edge, ~20% around 0, ~30% uniform, ~10% powers of two, ~10% gap rule.
pub const EXP_MIX: [(ExpBand, u64); 5] = [
    (ExpBand::Edge, 30),
    (ExpBand::Near0, 20),
    (ExpBand::Uniform, 30),
    (ExpBand::PowerOfTwo, 10),
    (ExpBand::Gap, 10),
];

/// The percentages TEST_PLAN section 7 prescribes, kept apart from `EXP_MIX` so a table edited by mistake is caught by `check_exp_mix`.
pub const EXP_MIX_SPEC_PERCENT: [u64; 5] = [30, 20, 30, 10, 10];

/// Checks a mix table against the plan: five bands in the documented order, the documented percentages, a total of 100.
pub fn check_exp_mix(mix: &[(ExpBand, u64)]) -> Result<(), String> {
    if mix.len() != EXP_MIX_SPEC_PERCENT.len() {
        return Err(format!(
            "mix has {} bands, want {}",
            mix.len(),
            EXP_MIX_SPEC_PERCENT.len()
        ));
    }
    let mut total = 0u64;
    for (i, ((band, pct), want)) in mix.iter().zip(EXP_MIX_SPEC_PERCENT.iter()).enumerate() {
        if band.index() != i {
            return Err(format!(
                "band {:?} is at position {}, want {:?}",
                band,
                i,
                ExpBand::from_index(i)
            ));
        }
        if pct != want {
            return Err(format!(
                "band {:?} has {}%, the plan says {}%",
                band, pct, want
            ));
        }
        total = total
            .checked_add(*pct)
            .ok_or_else(|| String::from("mix total overflows"))?;
    }
    if total != 100 {
        return Err(format!("mix totals {}%, want 100%", total));
    }
    Ok(())
}

/// One exponent draw: the band that produced it and the raw value.
pub struct ExpDraw {
    pub band: ExpBand,
    pub value: i64,
}

/// A value in `min..=max` (one draw, widening multiply).
pub fn uniform_in(rng: &mut SplitMix64, min: i64, max: i64) -> Result<i64, String> {
    if min > max {
        return Err(format!("empty range {}..={}", min, max));
    }
    let span = u128::from(max.abs_diff(min)) + 1;
    let offset = (u128::from(rng.next_u64()) * span) >> 64;
    let value = i128::from(min)
        .checked_add(i128::try_from(offset).map_err(|_| String::from("offset conversion"))?)
        .ok_or_else(|| String::from("uniform_in overflow"))?;
    i64::try_from(value).map_err(|_| format!("uniform_in result {} does not fit i64", value))
}

/// The Uniform band's exponent in `min..=max` from ONE output (Jeff, 2026-10-23 12:20 PM CT, cosmology-scale coverage). The top bit picks
/// linear (uniform over the range, as `uniform_in`) or log-scale: a bit length `0..=63` uniform, then a magnitude uniform within that
/// length and a random sign, so every scale from 1 to `2^62` (cosmology's `2^±405` included) is drawn about as often as any other. A
/// log-scale value outside the range falls back to the linear one.
fn wide_uniform_in(rng: &mut SplitMix64, min: i64, max: i64) -> Result<i64, String> {
    if min > max {
        return Err(format!("empty range {}..={}", min, max));
    }
    let x = rng.next_u64();
    let span = u128::from(max.abs_diff(min)) + 1;
    let offset = (u128::from(x << 1) * span) >> 64;
    let linear = i128::from(min)
        .checked_add(i128::try_from(offset).map_err(|_| String::from("offset conversion"))?)
        .ok_or_else(|| String::from("wide_uniform_in overflow"))?;
    let linear = i64::try_from(linear)
        .map_err(|_| format!("wide_uniform_in result {} does not fit i64", linear))?;
    if x >> 63 == 0 {
        return Ok(linear);
    }
    let len = u32::try_from((x >> 57) & 63).map_err(|_| String::from("bit length"))?;
    let negative = (x >> 56) & 1 == 1;
    let low = x & ((1u64 << 56) - 1);
    let mag = match len.checked_sub(1) {
        None => 0,
        Some(k) => {
            let top = 1u64 << k;
            let within = (u128::from(low) * u128::from(top)) >> 56;
            top | u64::try_from(within).map_err(|_| String::from("within"))?
        }
    };
    let magnitude = i64::try_from(mag).map_err(|_| String::from("magnitude"))?;
    let value = if negative { -magnitude } else { magnitude };
    if value < min || value > max {
        Ok(linear)
    } else {
        Ok(value)
    }
}

/// `base`, `base + gap` as exponents. `Err` (never a wrap) when the second does not fit `i64` or lies more than one tick outside `dom`.
pub fn gap_pair(base: i64, gap: i128, dom: &ExpDomain) -> Result<(i64, i64), String> {
    let sum = i128::from(base)
        .checked_add(gap)
        .ok_or_else(|| String::from("gap_pair overflow"))?;
    let second =
        i64::try_from(sum).map_err(|_| format!("base {} + gap {} does not fit i64", base, gap))?;
    if !dom.in_raw_range(second) {
        return Err(format!(
            "base {} + gap {} = {} is outside the domain",
            base, gap, second
        ));
    }
    Ok((base, second))
}

/// The gap classes of the controlled-gap rule, in table order; `p` is the precision, `guard` is `GUARD`.
#[derive(Debug, PartialEq, Eq)]
pub enum GapClass {
    Zero,
    One,
    Two,
    PMinus1,
    P,
    PPlus1,
    PGuardMinus1,
    PGuard,
    PGuardPlus1,
    Random,
}

impl GapClass {
    pub fn index(&self) -> usize {
        match self {
            GapClass::Zero => 0,
            GapClass::One => 1,
            GapClass::Two => 2,
            GapClass::PMinus1 => 3,
            GapClass::P => 4,
            GapClass::PPlus1 => 5,
            GapClass::PGuardMinus1 => 6,
            GapClass::PGuard => 7,
            GapClass::PGuardPlus1 => 8,
            GapClass::Random => 9,
        }
    }

    pub fn from_index(i: usize) -> Option<GapClass> {
        match i {
            0 => Some(GapClass::Zero),
            1 => Some(GapClass::One),
            2 => Some(GapClass::Two),
            3 => Some(GapClass::PMinus1),
            4 => Some(GapClass::P),
            5 => Some(GapClass::PPlus1),
            6 => Some(GapClass::PGuardMinus1),
            7 => Some(GapClass::PGuard),
            8 => Some(GapClass::PGuardPlus1),
            9 => Some(GapClass::Random),
            _ => None,
        }
    }

    /// The magnitude of the gap for precision `prec`; `random` is the magnitude used by `Random` (ignored by the others).
    pub fn magnitude(&self, prec: u32, random: i128) -> Result<i128, String> {
        let p = i128::from(prec);
        let g = i128::from(GUARD);
        let m = match self {
            GapClass::Zero => 0,
            GapClass::One => 1,
            GapClass::Two => 2,
            GapClass::PMinus1 => p - 1,
            GapClass::P => p,
            GapClass::PPlus1 => p + 1,
            GapClass::PGuardMinus1 => p + g - 1,
            GapClass::PGuard => p + g,
            GapClass::PGuardPlus1 => p + g + 1,
            GapClass::Random => random,
        };
        if m < 0 {
            Err(format!("negative gap magnitude for precision {}", prec))
        } else {
            Ok(m)
        }
    }
}

/// One gap draw: the class and the signed gap.
pub struct GapDraw {
    pub class: GapClass,
    pub gap: i128,
}

/// Draws a gap for precision `prec` (>= 1). Always consumes 3 outputs: class (uniform over the 10 classes), sign, random magnitude
/// (uniform in `0..=` the full span of `dom`, used only by `Random`). `Zero` ignores the sign.
pub fn draw_gap(rng: &mut SplitMix64, prec: u32, dom: &ExpDomain) -> Result<GapDraw, String> {
    if prec == 0 {
        return Err(String::from("precision must be >= 1"));
    }
    let class_index = usize::try_from(rng.below(10)).map_err(|_| String::from("class index"))?;
    let negative = rng.below(2) == 1;
    let span = u128::from(dom.max.abs_diff(dom.min));
    let wide = (u128::from(rng.next_u64()) * (span + 1)) >> 64;
    let random = i128::try_from(wide).map_err(|_| String::from("random gap conversion"))?;
    let class =
        GapClass::from_index(class_index).ok_or_else(|| format!("no gap class {}", class_index))?;
    let magnitude = class.magnitude(prec, random)?;
    let gap = if negative { -magnitude } else { magnitude };
    Ok(GapDraw { class, gap })
}

/// A pair of exponents `(base, base + gap)` for a controlled-gap partner draw (3 outputs, see `draw_gap`). `Err` when the partner would leave `dom`.
pub fn draw_exponent_pair(
    rng: &mut SplitMix64,
    base: i64,
    prec: u32,
    dom: &ExpDomain,
) -> Result<((i64, i64), GapDraw), String> {
    let g = draw_gap(rng, prec, dom)?;
    let pair = gap_pair(base, g.gap, dom)?;
    Ok((pair, g))
}

/// Draws one exponent with an explicit mix table. Consumption: 1 output for the band, then Edge 1, Near0 1, Uniform 1, PowerOfTwo 3
/// (bit, sign, delta), Gap 3 (see `draw_gap`). `anchor` and `prec` matter only in the Gap band (value = `anchor + gap`, or `anchor - gap`
/// when that leaves the domain, or `anchor` when both do).
pub fn draw_exponent_with(
    rng: &mut SplitMix64,
    dom: &ExpDomain,
    anchor: i64,
    prec: u32,
    mix: &[(ExpBand, u64)],
) -> Result<ExpDraw, String> {
    let r = rng.below(100);
    let mut acc = 0u64;
    let mut chosen: Option<usize> = None;
    for (band, pct) in mix.iter() {
        acc = acc
            .checked_add(*pct)
            .ok_or_else(|| String::from("mix total overflows"))?;
        if r < acc {
            chosen = Some(band.index());
            break;
        }
    }
    let band = chosen
        .and_then(ExpBand::from_index)
        .ok_or_else(|| format!("mix does not cover the draw {}", r))?;
    let value = match band {
        ExpBand::Edge => {
            let edge = dom.edge_values();
            let n = u64::try_from(edge.len()).map_err(|_| String::from("edge length"))?;
            let i = usize::try_from(rng.below(n)).map_err(|_| String::from("edge index"))?;
            *edge
                .get(i)
                .ok_or_else(|| String::from("empty edge table"))?
        }
        ExpBand::Near0 => {
            let k = i64::try_from(rng.below(129)).map_err(|_| String::from("near0 offset"))?;
            k - 64
        }
        ExpBand::Uniform => wide_uniform_in(rng, dom.min, dom.max)?,
        ExpBand::PowerOfTwo => {
            let bits = u64::from(64 - dom.max.leading_zeros());
            let k = u32::try_from(rng.below(bits)).map_err(|_| String::from("pow2 bit"))?;
            let negative = rng.below(2) == 1;
            let delta = i64::try_from(rng.below(3)).map_err(|_| String::from("pow2 delta"))? - 1;
            let pow = 1i64
                .checked_shl(k)
                .ok_or_else(|| String::from("pow2 shift"))?;
            let signed = if negative { -pow } else { pow };
            signed
                .checked_add(delta)
                .ok_or_else(|| String::from("pow2 overflow"))?
        }
        ExpBand::Gap => {
            let g = draw_gap(rng, prec, dom)?;
            match gap_pair(anchor, g.gap, dom) {
                Ok((_, v)) => v,
                Err(_) => match gap_pair(anchor, -g.gap, dom) {
                    Ok((_, v)) => v,
                    Err(_) => anchor,
                },
            }
        }
    };
    if !dom.in_raw_range(value) {
        return Err(format!(
            "band {:?} produced {} outside the domain",
            band, value
        ));
    }
    Ok(ExpDraw { band, value })
}

/// `draw_exponent_with` and the plan's mix (`EXP_MIX`).
pub fn draw_exponent(
    rng: &mut SplitMix64,
    dom: &ExpDomain,
    anchor: i64,
    prec: u32,
) -> Result<ExpDraw, String> {
    draw_exponent_with(rng, dom, anchor, prec, &EXP_MIX)
}

// ---- mantissas ------------------------------------------------------------------------------------------------------------

/// Mantissa patterns of a `Mag` (top bit set), in table order.
#[derive(Debug, PartialEq, Eq)]
pub enum MantPattern {
    Top,
    TopPlus1,
    AllOnes,
    MaxMinus1,
    C000,
    Aaab,
    Alt5555,
    SingleBit,
    TieAtP,
    TieMinus1,
    TiePlus1,
    Uniform,
}

impl MantPattern {
    pub fn index(&self) -> usize {
        match self {
            MantPattern::Top => 0,
            MantPattern::TopPlus1 => 1,
            MantPattern::AllOnes => 2,
            MantPattern::MaxMinus1 => 3,
            MantPattern::C000 => 4,
            MantPattern::Aaab => 5,
            MantPattern::Alt5555 => 6,
            MantPattern::SingleBit => 7,
            MantPattern::TieAtP => 8,
            MantPattern::TieMinus1 => 9,
            MantPattern::TiePlus1 => 10,
            MantPattern::Uniform => 11,
        }
    }

    pub fn from_index(i: usize) -> Option<MantPattern> {
        match i {
            0 => Some(MantPattern::Top),
            1 => Some(MantPattern::TopPlus1),
            2 => Some(MantPattern::AllOnes),
            3 => Some(MantPattern::MaxMinus1),
            4 => Some(MantPattern::C000),
            5 => Some(MantPattern::Aaab),
            6 => Some(MantPattern::Alt5555),
            7 => Some(MantPattern::SingleBit),
            8 => Some(MantPattern::TieAtP),
            9 => Some(MantPattern::TieMinus1),
            10 => Some(MantPattern::TiePlus1),
            11 => Some(MantPattern::Uniform),
            _ => None,
        }
    }
}

/// Percent per pattern, in table order (total 100).
pub const MANT_MIX: [(MantPattern, u64); 12] = [
    (MantPattern::Top, 8),
    (MantPattern::TopPlus1, 6),
    (MantPattern::AllOnes, 8),
    (MantPattern::MaxMinus1, 6),
    (MantPattern::C000, 6),
    (MantPattern::Aaab, 6),
    (MantPattern::Alt5555, 6),
    (MantPattern::SingleBit, 8),
    (MantPattern::TieAtP, 8),
    (MantPattern::TieMinus1, 6),
    (MantPattern::TiePlus1, 6),
    (MantPattern::Uniform, 26),
];

/// Checks a mantissa mix table: twelve patterns in table order, each with a positive weight, total 100.
pub fn check_mant_mix(mix: &[(MantPattern, u64)]) -> Result<(), String> {
    if mix.len() != 12 {
        return Err(format!("mantissa mix has {} patterns, want 12", mix.len()));
    }
    let mut total = 0u64;
    for (i, (pattern, pct)) in mix.iter().enumerate() {
        if pattern.index() != i {
            return Err(format!("pattern {:?} is at position {}", pattern, i));
        }
        if *pct == 0 {
            return Err(format!("pattern {:?} has weight 0", pattern));
        }
        total = total
            .checked_add(*pct)
            .ok_or_else(|| String::from("mix total overflows"))?;
    }
    if total != 100 {
        return Err(format!("mantissa mix totals {}%, want 100%", total));
    }
    Ok(())
}

/// The halfway-tie mantissa at precision `prec`: the top bit, plus the bit just below the `prec` kept bits (bit `63 - prec`), all lower bits zero.
/// For `prec >= 64` no bit is dropped: the tie degenerates to `TOP`.
pub fn tie_mantissa(prec: u32) -> Result<u64, String> {
    if prec == 0 {
        return Err(String::from("precision must be >= 1"));
    }
    let shift = 63u32.checked_sub(prec);
    match shift {
        Some(s) => Ok(TOP | (1u64 << s)),
        None => Ok(TOP),
    }
}

/// One mantissa draw: the pattern and a value with the top bit set.
pub struct MantDraw {
    pub pattern: MantPattern,
    pub value: u64,
}

/// Draws a mantissa for precision `prec` (>= 1) from `mix`. Always consumes 2 outputs: the pattern, then one auxiliary output
/// (the bit of `SingleBit`, the low bits of `Uniform`; unused by the fixed patterns).
pub fn draw_mantissa_with(
    rng: &mut SplitMix64,
    prec: u32,
    mix: &[(MantPattern, u64)],
) -> Result<MantDraw, String> {
    let r = rng.below(100);
    let aux = rng.next_u64();
    let mut acc = 0u64;
    let mut chosen: Option<usize> = None;
    for (pattern, pct) in mix.iter() {
        acc = acc
            .checked_add(*pct)
            .ok_or_else(|| String::from("mix total overflows"))?;
        if r < acc {
            chosen = Some(pattern.index());
            break;
        }
    }
    let pattern = chosen
        .and_then(MantPattern::from_index)
        .ok_or_else(|| format!("mantissa mix does not cover the draw {}", r))?;
    let tie = tie_mantissa(prec)?;
    let value = match pattern {
        MantPattern::Top => TOP,
        MantPattern::TopPlus1 => TOP + 1,
        MantPattern::AllOnes => u64::MAX,
        MantPattern::MaxMinus1 => u64::MAX - 1,
        MantPattern::C000 => 0xC000_0000_0000_0000,
        MantPattern::Aaab => 0xAAAA_AAAA_AAAA_AAAB,
        MantPattern::Alt5555 => 0x5555_5555_5555_5555 | TOP,
        MantPattern::SingleBit => {
            let bit = u32::try_from((u128::from(aux) * 63) >> 64)
                .map_err(|_| String::from("single bit"))?;
            TOP | (1u64 << bit)
        }
        MantPattern::TieAtP => tie,
        MantPattern::TieMinus1 => tie.checked_sub(1).map_or(TOP, |t| t.max(TOP)),
        MantPattern::TiePlus1 => tie + 1,
        MantPattern::Uniform => aux | TOP,
    };
    Ok(MantDraw { pattern, value })
}

/// `draw_mantissa_with` and the table `MANT_MIX`.
pub fn draw_mantissa(rng: &mut SplitMix64, prec: u32) -> Result<MantDraw, String> {
    draw_mantissa_with(rng, prec, &MANT_MIX)
}

// ---- precisions -----------------------------------------------------------------------------------------------------------

/// Where a precision comes from.
#[derive(Debug, PartialEq, Eq)]
pub enum PrecSource {
    /// The 41-value schedule `PS`.
    Schedule,
    /// The 18-value gate subset `PS_FAST`.
    Fast,
    /// Any `u32` precision `1..=u32::MAX` (CONTRACT 8.1; headroom far past `MAX_PREC`, Jeff 2026-10-23 12:33 PM CT): one draw in 8 is
    /// an edge of `PREC_EDGES`, the rest log-scale (bit length `1..=32` uniform, then uniform within it).
    Uniform,
    /// One of `op - 1`, `op`, `op + 1` for an operand precision `op >= 1`; a result of 0 or past `u32::MAX` becomes `op`.
    Operand(u32),
}

/// The precision edges `PrecSource::Uniform` draws on purpose: the small ones, `MAX_PREC` and both sides of it, and the top of `u32`.
pub const PREC_EDGES: [u32; 12] = [
    1,
    2,
    3,
    MAX_PREC - 1,
    MAX_PREC,
    MAX_PREC + 1,
    (1 << 31) - 1,
    1 << 31,
    (1 << 31) + 1,
    u32::MAX - 2,
    u32::MAX - 1,
    u32::MAX,
];

/// Draws a precision (always 1 output).
pub fn draw_precision(rng: &mut SplitMix64, source: &PrecSource) -> Result<u32, String> {
    match source {
        PrecSource::Schedule => {
            let i = usize::try_from(rng.below(41)).map_err(|_| String::from("schedule index"))?;
            match PS.get(i) {
                Some(&v) => Ok(v),
                None => Err(String::from("schedule index out of range")),
            }
        }
        PrecSource::Fast => {
            let i = usize::try_from(rng.below(18)).map_err(|_| String::from("fast index"))?;
            match PS_FAST.get(i) {
                Some(&v) => Ok(v),
                None => Err(String::from("fast index out of range")),
            }
        }
        PrecSource::Uniform => {
            let x = rng.next_u64();
            let low = u128::from(x & ((1u64 << 56) - 1));
            if x >> 61 == 0 {
                let n = u128::try_from(PREC_EDGES.len()).map_err(|_| String::from("edges"))?;
                let i = usize::try_from((low * n) >> 56).map_err(|_| String::from("edge index"))?;
                return match PREC_EDGES.get(i) {
                    Some(&edge) => Ok(edge),
                    None => Err(String::from("edge index out of range")),
                };
            }
            let k = u32::try_from((x >> 56) & 31).map_err(|_| String::from("bit length"))?;
            let top = 1u64 << k;
            let within =
                u64::try_from((low * u128::from(top)) >> 56).map_err(|_| String::from("within"))?;
            u32::try_from(top | within).map_err(|_| String::from("uniform precision"))
        }
        PrecSource::Operand(op) => {
            if *op == 0 {
                return Err(String::from("operand precision 0"));
            }
            let d = rng.below(3);
            let p = match d {
                0 => op.checked_sub(1),
                1 => Some(*op),
                _ => op.checked_add(1),
            };
            Ok(match p {
                Some(v) if v >= 1 => v,
                _ => *op,
            })
        }
    }
}

// ---- ball descriptors -----------------------------------------------------------------------------------------------------

/// A sign (both zeros are `Plus`/`Minus` with `mid_is_zero`).
#[derive(Debug, PartialEq, Eq)]
pub enum Sign {
    Plus,
    Minus,
}

impl Sign {
    pub fn symbol(&self) -> char {
        match self {
            Sign::Plus => '+',
            Sign::Minus => '-',
        }
    }
}

/// How the radius of a ball is chosen; `Zero`, `OneUlp`, `HalfUlp`, `Huge`, `Infinity` are derived from the mid by whoever builds the real ball,
/// `RandomMag` uses the drawn `rad_exp` and `rad_mant`.
#[derive(Debug, PartialEq, Eq)]
pub enum RadiusClass {
    Zero,
    OneUlp,
    HalfUlp,
    Huge,
    Infinity,
    RandomMag,
}

impl RadiusClass {
    pub fn index(&self) -> usize {
        match self {
            RadiusClass::Zero => 0,
            RadiusClass::OneUlp => 1,
            RadiusClass::HalfUlp => 2,
            RadiusClass::Huge => 3,
            RadiusClass::Infinity => 4,
            RadiusClass::RandomMag => 5,
        }
    }

    pub fn from_index(i: usize) -> Option<RadiusClass> {
        match i {
            0 => Some(RadiusClass::Zero),
            1 => Some(RadiusClass::OneUlp),
            2 => Some(RadiusClass::HalfUlp),
            3 => Some(RadiusClass::Huge),
            4 => Some(RadiusClass::Infinity),
            5 => Some(RadiusClass::RandomMag),
            _ => None,
        }
    }

    pub fn tag(&self) -> &'static str {
        match self {
            RadiusClass::Zero => "zero",
            RadiusClass::OneUlp => "ulp",
            RadiusClass::HalfUlp => "half",
            RadiusClass::Huge => "huge",
            RadiusClass::Infinity => "inf",
            RadiusClass::RandomMag => "mag",
        }
    }
}

/// Percent per radius class, in `RadiusClass::index` order (total 100).
pub const RADIUS_MIX: [u64; 6] = [20, 20, 15, 15, 10, 20];

/// A ball described as plain data (no zenith-float): the real ball is built later by L0b from this description.
pub struct BallDesc {
    pub sign: Sign,
    /// The mid is a (signed) zero; `mid_exp` and `mid_mant` are then drawn but ignored.
    pub mid_is_zero: bool,
    pub mid_exp: i64,
    pub mid_mant: u64,
    pub radius: RadiusClass,
    /// Independent of `mid_exp`; used by `RadiusClass::RandomMag`.
    pub rad_exp: i64,
    pub rad_mant: u64,
}

impl BallDesc {
    /// A whitespace-free text form for failure lines, e.g. `+8000000000000000e-5~r=ulp`, `-0~r=zero`, `+c000000000000000e7~r=mag(aaaaaaaaaaaaaaab,e-3)`.
    pub fn encode(&self) -> String {
        let mid = if self.mid_is_zero {
            format!("{}0", self.sign.symbol())
        } else {
            format!(
                "{}{:016x}e{}",
                self.sign.symbol(),
                self.mid_mant,
                self.mid_exp
            )
        };
        let rad = match self.radius {
            RadiusClass::RandomMag => format!("mag({:016x},e{})", self.rad_mant, self.rad_exp),
            RadiusClass::Zero
            | RadiusClass::OneUlp
            | RadiusClass::HalfUlp
            | RadiusClass::Huge
            | RadiusClass::Infinity => String::from(self.radius.tag()),
        };
        format!("{}~r={}", mid, rad)
    }
}

/// Draws a ball descriptor for precision `prec`. Consumption (frozen): sign 1, zero-mid 1, mid exponent (variable, see `draw_exponent_with`),
/// mid mantissa 2, radius class 1, radius exponent (variable, `Mag` domain, anchored at the mid exponent), radius mantissa 2.
/// A mid is zero with probability 1/20; the radius exponent is drawn independently of the mid exponent.
pub fn draw_ball(rng: &mut SplitMix64, prec: u32, mid_dom: &ExpDomain) -> Result<BallDesc, String> {
    let sign = if rng.below(2) == 1 {
        Sign::Minus
    } else {
        Sign::Plus
    };
    let mid_is_zero = rng.below(20) == 0;
    let mid_exp = draw_exponent(rng, mid_dom, 0, prec)?.value;
    let mid_mant = draw_mantissa(rng, prec)?.value;
    let r = rng.below(100);
    let mut acc = 0u64;
    let mut class_index = 0usize;
    for (i, pct) in RADIUS_MIX.iter().enumerate() {
        acc += pct;
        if r < acc {
            class_index = i;
            break;
        }
    }
    let radius =
        RadiusClass::from_index(class_index).ok_or_else(|| String::from("radius class"))?;
    let rad_exp = draw_exponent(rng, &ExpDomain::MAG, mid_exp, prec)?.value;
    let rad_mant = draw_mantissa(rng, prec)?.value;
    Ok(BallDesc {
        sign,
        mid_is_zero,
        mid_exp,
        mid_mant,
        radius,
        rad_exp,
        rad_mant,
    })
}

// ---- rounding modes -------------------------------------------------------------------------------------------------------

/// The six valid rounding modes (test-side enum; there is no wildcard `match` on it anywhere).
#[derive(Debug, PartialEq, Eq)]
pub enum RoundingMode {
    ToEven,
    ToOdd,
    Up,
    Down,
    ToZero,
    FromZero,
}

impl RoundingMode {
    /// All six modes, fresh values, in the order of `index`.
    pub fn all() -> core::array::IntoIter<RoundingMode, 6> {
        [
            RoundingMode::ToEven,
            RoundingMode::ToOdd,
            RoundingMode::Up,
            RoundingMode::Down,
            RoundingMode::ToZero,
            RoundingMode::FromZero,
        ]
        .into_iter()
    }

    pub fn index(&self) -> usize {
        match self {
            RoundingMode::ToEven => 0,
            RoundingMode::ToOdd => 1,
            RoundingMode::Up => 2,
            RoundingMode::Down => 3,
            RoundingMode::ToZero => 4,
            RoundingMode::FromZero => 5,
        }
    }

    pub fn name(&self) -> &'static str {
        match self {
            RoundingMode::ToEven => "ToEven",
            RoundingMode::ToOdd => "ToOdd",
            RoundingMode::Up => "Up",
            RoundingMode::Down => "Down",
            RoundingMode::ToZero => "ToZero",
            RoundingMode::FromZero => "FromZero",
        }
    }

    /// The mode with this name, `None` for anything else (including the name of the deliberate `None` marker).
    pub fn parse(name: &str) -> Option<RoundingMode> {
        RoundingMode::all().find(|m| m.name() == name)
    }
}

/// A rounding mode or the deliberate `None` marker: kept apart from `RoundingMode` so tests can assert that an op rejects `None`.
#[derive(Debug, PartialEq, Eq)]
pub enum ModeOrNone {
    Valid(RoundingMode),
    None,
}

impl ModeOrNone {
    /// The six valid modes followed by `None`.
    pub fn all_with_none() -> Vec<ModeOrNone> {
        let mut v: Vec<ModeOrNone> = RoundingMode::all().map(ModeOrNone::Valid).collect();
        v.push(ModeOrNone::None);
        v
    }

    pub fn is_valid(&self) -> bool {
        match self {
            ModeOrNone::Valid(_) => true,
            ModeOrNone::None => false,
        }
    }

    pub fn name(&self) -> &'static str {
        match self {
            ModeOrNone::Valid(m) => m.name(),
            ModeOrNone::None => "None",
        }
    }
}

/// Draws one of the six valid modes (1 output).
pub fn draw_mode(rng: &mut SplitMix64) -> Result<RoundingMode, String> {
    let i = usize::try_from(rng.below(6)).map_err(|_| String::from("mode index"))?;
    RoundingMode::all()
        .nth(i)
        .ok_or_else(|| String::from("mode index out of range"))
}
