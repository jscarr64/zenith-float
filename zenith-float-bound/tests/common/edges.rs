//! L0b Ball descriptor edge table (TEST_PLAN 7.1 edge rows, checklist L0b): named `BallDesc` cases built from `generators` (no Ball type, no
//! zenith-float). `Grid::Finite` cases have a P1 sample grid (`points::sample_points` is Ok); `Grid::Refused` ones (entire line, Huge radius
//! whose value is undecided) have none. `placeholder` marks a row that depends on `generators::ZF_MIN_NORMAL_EXPONENT`, which is a
//! placeholder until zf exposes a minimum-normal exponent (zf fe5c15d has only EXPONENT_MIN/EXPONENT_MAX). The control `2^200 + 1` needs
//! 201 mantissa bits and has no `BallDesc` (64-bit mantissas): it lives as an exact `points::Point` in points_selftest.

use super::generators::{
    tie_mantissa, BallDesc, RadiusClass, Sign, TOP, ZF_EXPONENT_MAX, ZF_EXPONENT_MIN,
    ZF_MIN_NORMAL_EXPONENT,
};

/// Whether a row has a finite P1 grid.
#[derive(Debug, PartialEq, Eq)]
pub enum Grid {
    Finite,
    Refused,
}

pub struct Edge {
    pub name: &'static str,
    pub desc: BallDesc,
    pub grid: Grid,
    pub placeholder: bool,
}

fn d(
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

fn e(name: &'static str, desc: BallDesc, grid: Grid) -> Edge {
    Edge {
        name,
        desc,
        grid,
        placeholder: false,
    }
}

/// The edge rows for precision `prec` (>= 1), in a fixed order.
pub fn edge_table(prec: u32) -> Result<Vec<Edge>, String> {
    let tie = tie_mantissa(prec)?;
    // 1 + 2^-60 = 0.(1 0^59 1) * 2^1 and 1 - 2^-60 = 0.(1^60) * 2^0, both exact in 64 bits.
    let one_plus = TOP | (1u64 << 3);
    let one_minus = u64::MAX << 4;
    let mut v = vec![
        e(
            "zero-ball(+0,r=0)",
            d(Sign::Plus, true, 0, 0, RadiusClass::Zero, 0, 0),
            Grid::Finite,
        ),
        e(
            "zero-ball(-0,r=0)",
            d(Sign::Minus, true, 0, 0, RadiusClass::Zero, 0, 0),
            Grid::Finite,
        ),
        e(
            "+0-mid,r=1/2",
            d(Sign::Plus, true, 0, 0, RadiusClass::RandomMag, 0, TOP),
            Grid::Finite,
        ),
        e(
            "-0-mid,r=1/2",
            d(Sign::Minus, true, 0, 0, RadiusClass::RandomMag, 0, TOP),
            Grid::Finite,
        ),
        e(
            "exact-point(1)",
            d(Sign::Plus, false, 1, TOP, RadiusClass::Zero, 0, 0),
            Grid::Finite,
        ),
        e(
            "exact-point(-1)",
            d(Sign::Minus, false, 1, TOP, RadiusClass::Zero, 0, 0),
            Grid::Finite,
        ),
        e(
            "entire-line",
            d(Sign::Plus, true, 0, 0, RadiusClass::Infinity, 0, 0),
            Grid::Refused,
        ),
        e(
            "entire-line(mid=1)",
            d(Sign::Plus, false, 1, TOP, RadiusClass::Infinity, 0, 0),
            Grid::Refused,
        ),
        e(
            "huge-radius",
            d(Sign::Plus, false, 1, TOP, RadiusClass::Huge, 0, 0),
            Grid::Refused,
        ),
        e(
            "max-exponent,1ulp",
            d(
                Sign::Plus,
                false,
                ZF_EXPONENT_MAX,
                u64::MAX,
                RadiusClass::OneUlp,
                0,
                0,
            ),
            Grid::Finite,
        ),
        e(
            "min-exponent,r=0",
            d(
                Sign::Minus,
                false,
                ZF_EXPONENT_MIN,
                TOP,
                RadiusClass::Zero,
                0,
                0,
            ),
            Grid::Finite,
        ),
        e(
            "min-exponent,1ulp",
            d(
                Sign::Plus,
                false,
                ZF_EXPONENT_MIN,
                TOP,
                RadiusClass::OneUlp,
                0,
                0,
            ),
            Grid::Finite,
        ),
        e(
            "tie-at-prec,half-ulp",
            d(Sign::Plus, false, 1, tie, RadiusClass::HalfUlp, 0, 0),
            Grid::Finite,
        ),
        e(
            "tie-at-prec,1ulp",
            d(Sign::Minus, false, 1, tie, RadiusClass::OneUlp, 0, 0),
            Grid::Finite,
        ),
        e(
            "2^200,1ulp",
            d(Sign::Plus, false, 201, TOP, RadiusClass::OneUlp, 0, 0),
            Grid::Finite,
        ),
        e(
            "2^200,r=1",
            d(Sign::Plus, false, 201, TOP, RadiusClass::RandomMag, 1, TOP),
            Grid::Finite,
        ),
        e(
            "1+2^-60,r=0",
            d(Sign::Plus, false, 1, one_plus, RadiusClass::Zero, 0, 0),
            Grid::Finite,
        ),
        e(
            "1-2^-60,r=0",
            d(Sign::Plus, false, 0, one_minus, RadiusClass::Zero, 0, 0),
            Grid::Finite,
        ),
        e(
            "1,r=2^-60",
            d(Sign::Plus, false, 1, TOP, RadiusClass::RandomMag, -59, TOP),
            Grid::Finite,
        ),
        // Forced (generator v5): one tick inside each zf exponent bound, and cosmology scale (2^405 = 0.1b*2^406, 2^-405 = 0.1b*2^-404).
        e(
            "max-exponent-1,1ulp",
            d(
                Sign::Plus,
                false,
                ZF_EXPONENT_MAX - 1,
                u64::MAX,
                RadiusClass::OneUlp,
                0,
                0,
            ),
            Grid::Finite,
        ),
        e(
            "min-exponent+1,1ulp",
            d(
                Sign::Plus,
                false,
                ZF_EXPONENT_MIN + 1,
                TOP,
                RadiusClass::OneUlp,
                0,
                0,
            ),
            Grid::Finite,
        ),
        e(
            "2^405,1ulp",
            d(Sign::Plus, false, 406, TOP, RadiusClass::OneUlp, 0, 0),
            Grid::Finite,
        ),
        e(
            "2^-405,1ulp",
            d(Sign::Plus, false, -404, TOP, RadiusClass::OneUlp, 0, 0),
            Grid::Finite,
        ),
        e(
            "-2^-405,r=0",
            d(Sign::Minus, false, -404, TOP, RadiusClass::Zero, 0, 0),
            Grid::Finite,
        ),
    ];
    v.push(Edge {
        name: "below-min-normal(placeholder),1ulp",
        desc: d(
            Sign::Plus,
            false,
            ZF_MIN_NORMAL_EXPONENT - 1,
            TOP,
            RadiusClass::OneUlp,
            0,
            0,
        ),
        grid: Grid::Finite,
        placeholder: true,
    });
    Ok(v)
}
