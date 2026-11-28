//! COVERS: common::Term::new, common::Term::is_zero, common::points::radius_value, common::points::sample_points, common::points::cmp_points, common::SampleKind::Seeded, common::SampleKind::Mid, common::SampleKind::Edge, common::SampleKind::Eighth
//! Self-tests of tests/common/points.rs (L0b P1 grid, TEST_PLAN 7.1). Oracle: hand-computed dyadic values stated in each test; the
//! grid order and the in-ball / monotone properties follow from the definition `mid + r*j/8`. The seeded `j` for seed 7 is pinned from
//! this implementation and only guards determinism of the draw.

use core::cmp::Ordering;

use crate::common::expect_eq;
use crate::common::generators::{BallDesc, RadiusClass, Sign, TOP};
use crate::common::points::{
    cmp_points, radius_value, sample_points, Point, SampleKind, Term, DEFAULT_MAX_SHIFT, MAX_T,
};
use crate::common::rng::SplitMix64;

fn desc(
    sign: Sign,
    zero: bool,
    mid_exp: i64,
    mid_mant: u64,
    radius: RadiusClass,
    rad_exp: i64,
    rad_mant: u64,
) -> BallDesc {
    BallDesc {
        sign,
        mid_is_zero: zero,
        mid_exp,
        mid_mant,
        radius,
        rad_exp,
        rad_mant,
    }
}

fn t(neg: bool, mant: u128, exp: i128) -> Result<Term, String> {
    Term::new(neg, mant, exp)
}

fn pt(neg: bool, mant: u128, exp: i128) -> Result<Point, String> {
    Ok(Point {
        base: t(neg, mant, exp)?,
        off: t(false, 0, 0)?,
    })
}

/// mid = 1 (0.1b * 2^1), r = 0.75 * 2^-2 = 3/16.
fn one_ball() -> BallDesc {
    desc(
        Sign::Plus,
        false,
        1,
        TOP,
        RadiusClass::RandomMag,
        -2,
        TOP | (TOP >> 1),
    )
}

#[test]
fn term_normalises() -> Result<(), String> {
    expect_eq(
        "12*2^0",
        &t(false, 12, 0)?,
        &Term {
            neg: false,
            mant: 3,
            exp: 2,
        },
    )?;
    expect_eq(
        "-0 is +0",
        &t(true, 0, 99)?,
        &Term {
            neg: false,
            mant: 0,
            exp: 0,
        },
    )?;
    expect_eq("zero", &t(true, 0, 5)?.is_zero(), &true)?;
    expect_eq("overflow", &t(false, 2, i128::MAX).is_err(), &true)
}

#[test]
fn grid_order_and_exact_values() -> Result<(), String> {
    let pts = sample_points(&one_ball(), 53, &mut SplitMix64::new(7), 10)?;
    expect_eq("count", &pts.len(), &18)?;
    let kinds: Vec<&SampleKind> = pts.iter().map(|(k, _)| k).collect();
    expect_eq("first", &kinds.first(), &Some(&&SampleKind::Mid))?;
    expect_eq(
        "edge+",
        &kinds.get(1),
        &Some(&&SampleKind::Edge { plus: true }),
    )?;
    expect_eq(
        "eighth -7",
        &kinds.get(3),
        &Some(&&SampleKind::Eighth { j: -7 }),
    )?;
    expect_eq(
        "eighth 7",
        &kinds.get(16),
        &Some(&&SampleKind::Eighth { j: 7 }),
    )?;
    // mid + r = 1 + 3/16 = 19/16; mid - r = 13/16; mid + 5r/8 = 1 + 15/128 = 143/128.
    let at = |i: usize| {
        pts.get(i)
            .map(|(_, p)| p)
            .ok_or_else(|| format!("no point {}", i))
    };
    expect_eq(
        "mid+r",
        &cmp_points(at(1)?, &pt(false, 19, -4)?, DEFAULT_MAX_SHIFT)?,
        &Ordering::Equal,
    )?;
    expect_eq(
        "mid-r",
        &cmp_points(at(2)?, &pt(false, 13, -4)?, DEFAULT_MAX_SHIFT)?,
        &Ordering::Equal,
    )?;
    expect_eq(
        "mid+5r/8",
        &cmp_points(at(14)?, &pt(false, 143, -7)?, DEFAULT_MAX_SHIFT)?,
        &Ordering::Equal,
    )?;
    expect_eq(
        "mid",
        &cmp_points(at(0)?, &pt(false, 1, 0)?, DEFAULT_MAX_SHIFT)?,
        &Ordering::Equal,
    )
}

#[test]
fn grid_is_inside_and_monotone() -> Result<(), String> {
    let b = desc(
        Sign::Minus,
        false,
        -300,
        TOP | 12345,
        RadiusClass::RandomMag,
        -350,
        TOP | 777,
    );
    let pts = sample_points(&b, 200, &mut SplitMix64::new(1), MAX_T)?;
    let get = |i: usize| {
        pts.get(i)
            .map(|(_, p)| p)
            .ok_or_else(|| format!("no point {}", i))
    };
    let (hi, lo) = (get(1)?, get(2)?);
    for (i, (k, p)) in pts.iter().enumerate() {
        if cmp_points(p, lo, DEFAULT_MAX_SHIFT)? == Ordering::Less
            || cmp_points(p, hi, DEFAULT_MAX_SHIFT)? == Ordering::Greater
        {
            return Err(format!("point {} {:?} outside the ball", i, k));
        }
    }
    for i in 3..16 {
        expect_eq(
            "eighths increase",
            &cmp_points(get(i)?, get(i + 1)?, DEFAULT_MAX_SHIFT)?,
            &Ordering::Less,
        )?;
    }
    Ok(())
}

#[test]
fn zero_radius_and_signed_zero_mid() -> Result<(), String> {
    let pts = sample_points(
        &desc(Sign::Plus, false, 3, TOP, RadiusClass::Zero, 0, 0),
        8,
        &mut SplitMix64::new(2),
        4,
    )?;
    for (_, p) in &pts {
        expect_eq(
            "all at mid 4",
            &cmp_points(p, &pt(false, 1, 2)?, DEFAULT_MAX_SHIFT)?,
            &Ordering::Equal,
        )?;
    }
    let plus = sample_points(
        &desc(Sign::Plus, true, 0, 0, RadiusClass::RandomMag, 0, TOP),
        8,
        &mut SplitMix64::new(3),
        4,
    )?;
    let minus = sample_points(
        &desc(Sign::Minus, true, 0, 0, RadiusClass::RandomMag, 0, TOP),
        8,
        &mut SplitMix64::new(3),
        4,
    )?;
    expect_eq("+0 and -0 mids give the same grid", &plus, &minus)?;
    expect_eq(
        "+0 mid, edge = r = 1/2",
        &plus.get(1).map(|(_, p)| {
            cmp_points(
                p,
                &Point {
                    base: Term {
                        neg: false,
                        mant: 0,
                        exp: 0,
                    },
                    off: Term {
                        neg: false,
                        mant: 1,
                        exp: -1,
                    },
                },
                64,
            )
        }),
        &Some(Ok(Ordering::Equal)),
    )
}

#[test]
fn ulp_radii_and_refusals() -> Result<(), String> {
    let b = |r: RadiusClass| desc(Sign::Plus, false, 1, TOP, r, 0, 0);
    expect_eq(
        "one ulp at 53",
        &radius_value(&b(RadiusClass::OneUlp), 53)?,
        &Some((1, -52)),
    )?;
    expect_eq(
        "half ulp at 53",
        &radius_value(&b(RadiusClass::HalfUlp), 53)?,
        &Some((1, -53)),
    )?;
    expect_eq("zero", &radius_value(&b(RadiusClass::Zero), 53)?, &None)?;
    let mut rng = SplitMix64::new(4);
    expect_eq(
        "huge refused",
        &sample_points(&b(RadiusClass::Huge), 53, &mut rng, 4).is_err(),
        &true,
    )?;
    expect_eq(
        "entire line refused",
        &sample_points(&b(RadiusClass::Infinity), 53, &mut rng, 4).is_err(),
        &true,
    )?;
    expect_eq(
        "ulp of zero mid refused",
        &radius_value(&desc(Sign::Plus, true, 0, 0, RadiusClass::OneUlp, 0, 0), 8).is_err(),
        &true,
    )?;
    expect_eq(
        "t too big",
        &sample_points(&one_ball(), 53, &mut rng, MAX_T + 1).is_err(),
        &true,
    )?;
    expect_eq(
        "prec 0",
        &sample_points(&one_ball(), 0, &mut rng, 4).is_err(),
        &true,
    )
}

#[test]
fn seeded_point_is_deterministic_and_in_range() -> Result<(), String> {
    let j_of = |seed: u64, tt: u32| -> Result<(i64, u32), String> {
        match sample_points(&one_ball(), 53, &mut SplitMix64::new(seed), tt)?.pop() {
            Some((SampleKind::Seeded { j, t }, _)) => Ok((j, t)),
            other => Err(format!("last point not seeded: {:?}", other)),
        }
    };
    expect_eq("repeatable", &j_of(7, 10)?, &j_of(7, 10)?)?;
    expect_eq("pinned seed 7", &j_of(7, 10)?, &(PINNED_J, 10))?;
    let (mut lo, mut hi) = (false, false);
    for s in 0..400u64 {
        let (j, _) = j_of(s, 2)?;
        if !(-4..=4).contains(&j) {
            return Err(format!("j = {} outside -4..=4", j));
        }
        lo |= j == -4;
        hi |= j == 4;
    }
    expect_eq("both ends reached", &(lo, hi), &(true, true))
}

const PINNED_J: i64 = -226;

#[test]
fn cmp_cancellation_and_spread() -> Result<(), String> {
    // (2^200 + 1) vs 2^200: exact, needs the low bit.
    let a = Point {
        base: t(false, 1, 200)?,
        off: t(false, 1, 0)?,
    };
    expect_eq(
        "2^200+1 > 2^200",
        &cmp_points(&a, &pt(false, 1, 200)?, DEFAULT_MAX_SHIFT)?,
        &Ordering::Greater,
    )?;
    // 1 + 2^-60 vs 1 - 2^-60.
    let p = Point {
        base: t(false, 1, 0)?,
        off: t(false, 1, -60)?,
    };
    let m = Point {
        base: t(false, 1, 0)?,
        off: t(true, 1, -60)?,
    };
    expect_eq(
        "1+2^-60 > 1-2^-60",
        &cmp_points(&p, &m, DEFAULT_MAX_SHIFT)?,
        &Ordering::Greater,
    )?;
    expect_eq(
        "negatives",
        &cmp_points(&pt(true, 3, 0)?, &pt(true, 1, 1)?, DEFAULT_MAX_SHIFT)?,
        &Ordering::Less,
    )?;
    expect_eq(
        "equal zero",
        &cmp_points(&pt(true, 0, 0)?, &pt(false, 0, 9)?, 0)?,
        &Ordering::Equal,
    )?;
    expect_eq(
        "spread refused",
        &cmp_points(&a, &pt(false, 1, 0)?, 199).is_err(),
        &true,
    )?;
    expect_eq(
        "extreme exponents refused",
        &cmp_points(
            &pt(false, 1, i128::MAX)?,
            &pt(false, 1, i128::MIN)?,
            DEFAULT_MAX_SHIFT,
        )
        .is_err(),
        &true,
    )
}

/// Headroom (Jeff, 2026-10-23 12:33 PM CT): a grid exists for every `i64` mid exponent at every `u32` precision, including the
/// cosmology-scale and extreme ones; the ulp exponent `mid_exp - prec` and the term exponents leave `i64` there, so they must not fail.
#[test]
fn grids_exist_at_the_i64_exponent_extremes_for_every_u32_precision() -> Result<(), String> {
    for mid_exp in [i64::MIN, -405, 405, i64::MAX] {
        for half in [false, true] {
            for prec in [1u32, 53, 4096, u32::MAX] {
                let radius = if half {
                    RadiusClass::HalfUlp
                } else {
                    RadiusClass::OneUlp
                };
                let b = desc(Sign::Minus, false, mid_exp, TOP, radius, 0, 0);
                sample_points(&b, prec, &mut SplitMix64::new(3), MAX_T)?;
            }
        }
    }
    let low = desc(
        Sign::Plus,
        false,
        i64::MIN,
        TOP,
        RadiusClass::RandomMag,
        i64::MIN,
        TOP,
    );
    sample_points(&low, 53, &mut SplitMix64::new(3), MAX_T)?;
    let want = i128::from(i64::MIN)
        .checked_sub(i128::from(u32::MAX))
        .ok_or_else(|| String::from("want"))?;
    expect_eq(
        "one ulp at i64::MIN, prec u32::MAX",
        &radius_value(
            &desc(Sign::Plus, false, i64::MIN, TOP, RadiusClass::OneUlp, 0, 0),
            u32::MAX,
        )?,
        &Some((1, want)),
    )
}
