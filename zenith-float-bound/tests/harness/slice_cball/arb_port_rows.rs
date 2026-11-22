//! Slice CB rows for the Arb port of `CBall::mul` / `CBall::div` (16.CB.4 / 16.CB.7; FLINT 3.6.0 `acb_mul`, `acb_div`, `acb_inv`).
//! Gold: tests/goldens/slices/cball/cball_arb_rows.txt (scripts/gen_cball_arb_goldens.py). Each row carries the exact image of the
//! input rectangles (rounded outward) and FLINT's own enclosure at the same prec. A result must contain the image, must be all reals
//! where the image is unbounded or undefined, and must be no looser than FLINT: finite where FLINT is finite, radius at most 2× FLINT's.
//! Classes: exact, inexact (relative radius 2^-8, 2^-30, 2^-(p-2)), an exact-zero part in each position, a zero midpoint with a
//! radius, an infinite radius in each position (alone and with an exact-zero partner), mixed signs, huge and tiny exponents,
//! cancellation, and for div a divisor rectangle containing / touching 0, the exact zero divisor, axis divisors reaching 0, wide
//! balls, a thin strip away from 0 and a short exact divisor at prec > 256. Precs 2, 53, 64, 128, 256, 300, 1024, 2048.

use super::*;

const ARB_ROWS: usize = 646;

enum Part {
    /// The whole line: an unbounded or undefined image, or a non-finite FLINT part.
    All,
    /// Two values: (lo, hi) for an image, (mid, rad) for FLINT.
    Two(ExactNum, ExactNum),
}

struct ArbRow {
    op: String,
    prec: usize,
    class: String,
    inputs: Vec<(ExactNum, Option<ExactNum>)>,
    image: [Part; 2],
    flint: [Part; 2],
}

fn take_rad(f: &[&str], i: &mut usize) -> Result<Option<ExactNum>, String> {
    if f.get(*i) == Some(&"inf") {
        *i += 1;
        return Ok(None);
    }
    take_val(f, i).map(Some)
}

fn take_part(f: &[&str], i: &mut usize, all: &str) -> Result<Part, String> {
    if f.get(*i) == Some(&all) {
        *i += 1;
        return Ok(Part::All);
    }
    let a = take_val(f, i)?;
    let b = take_val(f, i)?;
    Ok(Part::Two(a, b))
}

fn arb_rows() -> Result<Vec<ArbRow>, String> {
    let path = golden_dir()
        .join("slices")
        .join("cball")
        .join("cball_arb_rows.txt");
    let text = fs::read_to_string(&path).map_err(|e| format!("{}: {}", path.display(), e))?;
    let mut rows = Vec::new();
    for line in text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
    {
        let f: Vec<&str> = line.split(' ').collect();
        let op = String::from(*f.first().ok_or("empty row")?);
        let prec: usize = f
            .get(1)
            .ok_or("short row")?
            .parse()
            .map_err(|_| "bad prec")?;
        let class = String::from(*f.get(2).ok_or("short row")?);
        let mut i = 3;
        let mut inputs = Vec::new();
        for _ in 0..4 {
            let m = take_val(&f, &mut i)?;
            let r = take_rad(&f, &mut i)?;
            inputs.push((m, r));
        }
        let image = [take_part(&f, &mut i, "ALL")?, take_part(&f, &mut i, "ALL")?];
        let flint = [take_part(&f, &mut i, "INF")?, take_part(&f, &mut i, "INF")?];
        if i != f.len() {
            return Err(format!("trailing fields in row {} {} {}", op, prec, class));
        }
        rows.push(ArbRow {
            op,
            prec,
            class,
            inputs,
            image,
            flint,
        });
    }
    if rows.len() != ARB_ROWS {
        return Err(format!("want {} Arb rows, got {}", ARB_ROWS, rows.len()));
    }
    Ok(rows)
}

fn arb_part(r: &ArbRow, k: usize) -> Result<Ball, String> {
    let (m, rad) = r.inputs.get(k).ok_or("input index")?;
    let mag = match rad {
        None => Mag::infinity(),
        Some(v) => Mag::from_exact_num(v),
    };
    ok(
        "from_mid_rad",
        Ball::from_mid_rad(m.add_full_prec(&num(0)), mag, r.prec),
    )
}

fn arb_eval(r: &ArbRow) -> Result<CBall, String> {
    let a = ok("CBall::new", CBall::new(arb_part(r, 0)?, arb_part(r, 1)?))?;
    let b = ok("CBall::new", CBall::new(arb_part(r, 2)?, arb_part(r, 3)?))?;
    match r.op.as_str() {
        "mul" => ok("mul", a.mul(&b, r.prec)),
        "div" => ok("div", a.div(&b, r.prec)),
        other => Err(format!("unknown op {:?}", other)),
    }
}

fn want_image(what: &str, b: &Ball, image: &Part) -> Result<(), String> {
    match image {
        Part::All => {
            if b.rad().is_infinity() {
                Ok(())
            } else {
                Err(format!(
                    "{}: image is unbounded or undefined, got finite {:?}",
                    what, b
                ))
            }
        }
        Part::Two(lo, hi) => {
            want_contains(&format!("{} image lo", what), b, lo)?;
            want_contains(&format!("{} image hi", what), b, hi)
        }
    }
}

fn want_vs_flint(what: &str, b: &Ball, flint: &Part) -> Result<(), String> {
    match flint {
        Part::All => Ok(()),
        Part::Two(_, rad) => {
            if b.rad().is_infinity() {
                return Err(format!(
                    "{}: FLINT is finite (rad {:?}), got all reals",
                    what, rad
                ));
            }
            let bound = Mag::from_exact_num(&rad.mul_full_prec(&num(2)));
            if b.rad() > &bound {
                return Err(format!(
                    "{}: rad {:?} > 2× FLINT rad {:?}",
                    what,
                    b.rad(),
                    rad
                ));
            }
            Ok(())
        }
    }
}

fn arb_check(want_flint: bool) -> Result<(), String> {
    let mut bad = Vec::new();
    for r in arb_rows()?.iter() {
        let what = format!("{} p={} {}", r.op, r.prec, r.class);
        let z = arb_eval(r)?;
        for (k, part) in [z.re(), z.im()].into_iter().enumerate() {
            let name = format!("{} {}", what, if k == 0 { "re" } else { "im" });
            let res = if want_flint {
                want_vs_flint(&name, part, r.flint.get(k).ok_or("part index")?)
            } else {
                want_image(&name, part, r.image.get(k).ok_or("part index")?)
            };
            if let Err(e) = res {
                bad.push(e);
            }
        }
    }
    if bad.is_empty() {
        Ok(())
    } else {
        Err(format!("{} failing parts:\n{}", bad.len(), bad.join("\n")))
    }
}

#[test]
fn cball_arb_rows_contain_exact_image() -> Result<(), String> {
    arb_check(false)
}

#[test]
fn cball_arb_rows_no_looser_than_flint() -> Result<(), String> {
    arb_check(true)
}

/// zf exponents stop at EXPONENT_MAX (Arb's do not): a part whose midpoint overflows is all reals (CONTRACT 10.4); the other part
/// is untouched.
#[test]
fn cball_mul_div_overflow_part_is_all_reals() -> Result<(), String> {
    let big = || pow2(zenith_float_num::EXPONENT_MAX - 2);
    let z = ok(
        "new",
        CBall::new(ok("exact", Ball::exact(big()?, 64))?, pt(0, 64)?),
    )?;
    let w = ok(
        "new",
        CBall::new(ok("exact", Ball::exact(big()?, 64))?, pt(0, 64)?),
    )?;
    let p = ok("mul", z.mul(&w, 64))?;
    if p.re().is_informative() {
        return Err(format!(
            "overflowing re part must be all reals, got {:?}",
            p
        ));
    }
    want_point_part("mul im", p.im(), 0)?;
    // 4 / 2^(EXPONENT_MIN - 1) = 2^(EXPONENT_MAX + 4): past the top binade.
    let tiny = pow2(zenith_float_num::EXPONENT_MIN - 1)?;
    let d = ok(
        "new",
        CBall::new(ok("exact", Ball::exact(tiny, 64))?, pt(0, 64)?),
    )?;
    let q = ok("div", cpt(4, 0, 64)?.div(&d, 64))?;
    if q.re().is_informative() {
        return Err(format!(
            "overflowing quotient part must be all reals, got {:?}",
            q
        ));
    }
    want_point_part("div im", q.im(), 0)
}

fn want_point_part(what: &str, b: &Ball, v: i64) -> Result<(), String> {
    if b.mid() != &num(v) || !b.rad().is_zero() {
        return Err(format!("{}: want exact {}, got {:?}", what, v, b));
    }
    Ok(())
}
