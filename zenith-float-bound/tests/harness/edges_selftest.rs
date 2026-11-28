//! COVERS: common::edges::edge_table, common::Grid::Finite, common::Grid::Refused
//! Self-tests of tests/common/edges.rs. Oracle: the values named in each row (hand-derived: 1 = 0.1b*2^1, 2^200 = 0.1b*2^201,
//! 1-2^-60 = 0.(1^60)b*2^0) checked exactly through points::cmp_points; Finite rows must have a grid, Refused rows must not.

use core::cmp::Ordering;

use crate::common::edges::{edge_table, Edge, Grid};
use crate::common::expect_eq;
use crate::common::generators::PS;
use crate::common::points::{cmp_points, sample_points, Point, Term, DEFAULT_MAX_SHIFT};
use crate::common::rng::SplitMix64;

fn find<'a>(rows: &'a [Edge], name: &str) -> Result<&'a Edge, String> {
    rows.iter()
        .find(|r| r.name == name)
        .ok_or_else(|| format!("no row {}", name))
}

fn mid_of(r: &Edge) -> Result<Point, String> {
    let pts = sample_points(&r.desc, 53, &mut SplitMix64::new(0), 4)?;
    pts.into_iter()
        .next()
        .map(|(_, p)| p)
        .ok_or_else(|| String::from("empty grid"))
}

fn dy(neg: bool, mant: u128, exp: i128, off: Term) -> Result<Point, String> {
    Ok(Point {
        base: Term::new(neg, mant, exp)?,
        off,
    })
}

#[test]
fn every_precision_has_unique_rows_and_the_right_grid() -> Result<(), String> {
    for prec in PS {
        let rows = edge_table(prec)?;
        expect_eq("row count", &rows.len(), &25)?;
        for (i, r) in rows.iter().enumerate() {
            for s in rows.iter().skip(i + 1) {
                if r.name == s.name || r.desc.encode() == s.desc.encode() {
                    return Err(format!(
                        "prec {}: rows {} and {} collide",
                        prec, r.name, s.name
                    ));
                }
            }
            let got = sample_points(&r.desc, prec, &mut SplitMix64::new(1), 8);
            match (&r.grid, got) {
                (Grid::Finite, Ok(_)) | (Grid::Refused, Err(_)) => {}
                (Grid::Finite, Err(e)) => {
                    return Err(format!("prec {}: {} has no grid: {}", prec, r.name, e))
                }
                (Grid::Refused, Ok(_)) => {
                    return Err(format!("prec {}: {} has a grid", prec, r.name))
                }
            }
        }
        expect_eq(
            "one placeholder row",
            &rows.iter().filter(|r| r.placeholder).count(),
            &1,
        )?;
    }
    Ok(())
}

#[test]
fn row_values_are_exact() -> Result<(), String> {
    let rows = edge_table(53)?;
    let z = || Term::new(false, 0, 0);
    let eq = |name: &str, want: Point| -> Result<(), String> {
        expect_eq(
            name,
            &cmp_points(&mid_of(find(&rows, name)?)?, &want, DEFAULT_MAX_SHIFT)?,
            &Ordering::Equal,
        )
    };
    eq("exact-point(1)", dy(false, 1, 0, z()?)?)?;
    eq("exact-point(-1)", dy(true, 1, 0, z()?)?)?;
    eq("2^200,1ulp", dy(false, 1, 200, z()?)?)?;
    eq("1+2^-60,r=0", dy(false, 1, 0, Term::new(false, 1, -60)?)?)?;
    eq("1-2^-60,r=0", dy(false, 1, 0, Term::new(true, 1, -60)?)?)?;
    eq("zero-ball(-0,r=0)", dy(false, 0, 0, z()?)?)?;
    // forced cosmology-scale rows: 2^405 = 0.1b*2^406, 2^-405 = 0.1b*2^-404
    eq("2^405,1ulp", dy(false, 1, 405, z()?)?)?;
    eq("2^-405,1ulp", dy(false, 1, -405, z()?)?)?;
    eq("-2^-405,r=0", dy(true, 1, -405, z()?)?)?;
    let edge = |name: &str| -> Result<Point, String> {
        let pts = sample_points(&find(&rows, name)?.desc, 53, &mut SplitMix64::new(0), 4)?;
        pts.into_iter()
            .nth(1)
            .map(|(_, p)| p)
            .ok_or_else(|| String::from("no edge"))
    };
    expect_eq(
        "2^200 + r=1 edge is 2^200+1",
        &cmp_points(
            &edge("2^200,r=1")?,
            &dy(false, 1, 200, Term::new(false, 1, 0)?)?,
            DEFAULT_MAX_SHIFT,
        )?,
        &Ordering::Equal,
    )?;
    expect_eq(
        "1 + 2^-60 edge",
        &cmp_points(
            &edge("1,r=2^-60")?,
            &dy(false, 1, 0, Term::new(false, 1, -60)?)?,
            DEFAULT_MAX_SHIFT,
        )?,
        &Ordering::Equal,
    )?;
    expect_eq(
        "1ulp at 53 above 2^200",
        &cmp_points(
            &edge("2^200,1ulp")?,
            &dy(false, 1, 200, Term::new(false, 1, 148)?)?,
            DEFAULT_MAX_SHIFT,
        )?,
        &Ordering::Equal,
    )
}
