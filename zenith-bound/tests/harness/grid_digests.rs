//! COVERS: common::Digest::count, common::Digest::fold_bytes, common::Digest::fold_i64, common::Digest::fold_u64, common::Digest::row, common::Digest::start, common::Digest::value
//! S13 (TEST_PLAN P10, section 13): deterministic digests of the golden and grid runs, one row per suite, compared across targets,
//! profiles and test-thread counts by tools/digests.sh. Suites: `splitmix` (the rng stream), `generators` (S2 draws), `mag_grid` (every Mag
//! op over an edge grid of operands, plus the ulp helpers over PS-fast) and `ball_rows` (the exact-integer recomputation of every row of
//! tests/goldens/ball_exact_ops.txt). The expected digests are pinned in tools/digest_expected.txt; a mismatch prints the actual rows.
//! The digest file of the current run is written to `$CARGO_TARGET_TMPDIR/grid_digests.txt` (rows only matter to the script; the `#` header
//! lines name the arch and pointer width and are NOT compared). Oracle: FNV-1a published vectors; sensitivity and thread-independence tests.

use crate::common::ballrows::{check_row, parse_row, Checked};
use crate::common::digest::Digest;
use crate::common::expect_eq;
use crate::common::generators::{
    draw_ball, draw_exponent, draw_mantissa, draw_precision, ExpDomain, PrecSource, PS_FAST,
};
use crate::common::goldens::golden_dir;
use crate::common::rng::SplitMix64;
use crate::common::seeds::{seed_for, CASE_CLASS_GATE};
use zenith_bound::mag::Mag;

const EXPS: [i64; 25] = [
    i64::MIN,
    i64::MIN + 1,
    i64::MIN + 2,
    i64::MIN + 3,
    -1000,
    -130,
    -128,
    -127,
    -65,
    -64,
    -63,
    -1,
    0,
    1,
    63,
    64,
    65,
    127,
    128,
    130,
    1000,
    i64::MAX - 3,
    i64::MAX - 2,
    i64::MAX - 1,
    i64::MAX,
];

const MANTS: [u64; 8] = [
    0x8000_0000_0000_0000,
    0x8000_0000_0000_0001,
    u64::MAX,
    u64::MAX - 1,
    0xC000_0000_0000_0000,
    0xAAAA_AAAA_AAAA_AAAB,
    1,
    3,
];

enum Binary {
    AddUp,
    SubLower,
    MulUp,
    MulLower,
    DivUp,
    DivLower,
}

impl Binary {
    fn name(&self) -> &'static str {
        match self {
            Binary::AddUp => "add_up",
            Binary::SubLower => "sub_lower",
            Binary::MulUp => "mul_up",
            Binary::MulLower => "mul_lower",
            Binary::DivUp => "div_up",
            Binary::DivLower => "div_lower",
        }
    }

    fn apply(&self, a: &Mag, b: &Mag) -> Mag {
        match self {
            Binary::AddUp => a.add_up(b),
            Binary::SubLower => a.sub_lower(b),
            Binary::MulUp => a.mul_up(b),
            Binary::MulLower => a.mul_lower(b),
            Binary::DivUp => a.div_up(b),
            Binary::DivLower => a.div_lower(b),
        }
    }
}

const BINARY: [Binary; 6] = [
    Binary::AddUp,
    Binary::SubLower,
    Binary::MulUp,
    Binary::MulLower,
    Binary::DivUp,
    Binary::DivLower,
];

fn rng() -> Result<SplitMix64, String> {
    seed_for("common", CASE_CLASS_GATE)
        .map(SplitMix64::new)
        .ok_or_else(|| String::from("suite common is not in the seeds table"))
}

fn fold_mag(d: &mut Digest, m: &Mag) {
    match (m.exp(), m.mant()) {
        (Some(e), Some(mant)) => {
            d.fold_u64(2);
            d.fold_i64(e);
            d.fold_u64(mant);
        }
        (None, None) => d.fold_u64(1),
        (None, Some(_)) | (Some(_), None) => d.fold_u64(3),
    }
}

fn splitmix_suite() -> Result<Digest, String> {
    let mut r = rng()?;
    let mut d = Digest::start();
    for _ in 0..4096 {
        d.fold_u64(r.next_u64());
    }
    for _ in 0..4096 {
        d.fold_u64(r.below(1000));
    }
    Ok(d)
}

fn generators_suite() -> Result<Digest, String> {
    let mut r = rng()?;
    let mut d = Digest::start();
    for _ in 0..2000 {
        let prec = draw_precision(&mut r, &PrecSource::Fast)?;
        d.fold_u64(u64::from(prec));
        d.fold_i64(draw_exponent(&mut r, &ExpDomain::MAG, 0, prec)?.value);
        d.fold_u64(draw_mantissa(&mut r, prec)?.value);
        let b = draw_ball(&mut r, prec, &ExpDomain::MAG)?;
        d.fold_u64(u64::from(u32::from(b.sign.symbol())));
        d.fold_u64(u64::from(b.mid_is_zero));
        d.fold_i64(b.mid_exp);
        d.fold_u64(b.mid_mant);
        d.fold_u64(u64::try_from(b.radius.index()).map_err(|_| String::from("radius index"))?);
        d.fold_i64(b.rad_exp);
        d.fold_u64(b.rad_mant);
    }
    Ok(d)
}

fn operands() -> Result<(Vec<Mag>, u64), String> {
    let mut v = vec![Mag::zero(), Mag::infinity()];
    let mut rejected = 0u64;
    for exp in EXPS.iter() {
        for mant in MANTS.iter() {
            match Mag::finite(*mant, *exp) {
                Ok(m) => v.push(m),
                Err(_) => rejected += 1,
            }
        }
    }
    Ok((v, rejected))
}

fn mag_suite() -> Result<Digest, String> {
    let (ops, rejected) = operands()?;
    let mut d = Digest::start();
    d.fold_u64(rejected);
    d.fold_u64(u64::try_from(ops.len()).map_err(|_| String::from("operand count"))?);
    for op in BINARY.iter() {
        d.fold_bytes(op.name().as_bytes());
        for a in ops.iter() {
            for b in ops.iter() {
                fold_mag(&mut d, &op.apply(a, b));
            }
        }
    }
    d.fold_bytes(b"sqrt_up");
    for a in ops.iter() {
        fold_mag(&mut d, &a.sqrt_up());
    }
    d.fold_bytes(b"sqrt_lower");
    for a in ops.iter() {
        fold_mag(&mut d, &a.sqrt_lower());
    }
    for exp in EXPS.iter() {
        for prec in PS_FAST.iter() {
            fold_mag(
                &mut d,
                &Mag::ulp_mag(*exp, usize::try_from(*prec).unwrap_or(usize::MAX)),
            );
            fold_mag(
                &mut d,
                &Mag::half_ulp(*exp, usize::try_from(*prec).unwrap_or(usize::MAX)),
            );
        }
    }
    Ok(d)
}

fn ball_rows_suite() -> Result<Digest, String> {
    let path = golden_dir().join("ball_exact_ops.txt");
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{:?}: {}", path, e))?;
    let mut d = Digest::start();
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let tag = match parse_row(line).and_then(|row| check_row(&row)) {
            Ok(Checked::Recomputed) => 1,
            Ok(Checked::Structural) => 2,
            Err(_) => 3,
        };
        d.fold_u64(tag);
        d.fold_bytes(line.as_bytes());
    }
    Ok(d)
}

/// The rows of all suites, in the fixed order of the file.
fn rows_sequential() -> Result<String, String> {
    let mut s = String::new();
    s.push_str(&splitmix_suite()?.row("splitmix"));
    s.push_str(&generators_suite()?.row("generators"));
    s.push_str(&mag_suite()?.row("mag_grid"));
    s.push_str(&ball_rows_suite()?.row("ball_rows"));
    Ok(s)
}

fn expected_rows() -> Result<String, String> {
    let path = crate::common::goldens::repo_root()
        .join("tools")
        .join("digest_expected.txt");
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{:?}: {}", path, e))?;
    let mut s = String::new();
    for l in text.lines().filter(|l| !l.starts_with('#')) {
        s.push_str(l);
        s.push('\n');
    }
    Ok(s)
}

#[test]
fn fnv1a_matches_the_published_vectors() -> Result<(), String> {
    let mut d = Digest::start();
    expect_eq("empty", &d.value(), &0xCBF2_9CE4_8422_2325)?;
    d.fold_bytes(b"a");
    expect_eq("a", &d.value(), &0xAF63_DC4C_8601_EC8C)?;
    let mut f = Digest::start();
    f.fold_bytes(b"foobar");
    expect_eq("foobar", &f.value(), &0x8594_4171_F739_67E8)?;
    expect_eq("foobar byte count", &f.count(), &6)?;
    expect_eq("row text", &f.row("x").as_str(), &"x 6 85944171f73967e8\n")
}

#[test]
fn integers_fold_as_little_endian_bytes_and_order_matters() -> Result<(), String> {
    let mut a = Digest::start();
    a.fold_u64(0x0102_0304_0506_0708);
    let mut b = Digest::start();
    b.fold_bytes(&[8, 7, 6, 5, 4, 3, 2, 1]);
    expect_eq("u64 is its LE bytes", &a.value(), &b.value())?;
    expect_eq("word count of a u64", &a.count(), &8)?;
    let mut n = Digest::start();
    n.fold_i64(-2);
    let mut m = Digest::start();
    m.fold_bytes(&[0xFE, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF]);
    expect_eq("i64 is two's complement LE", &n.value(), &m.value())?;
    let (mut x, mut y) = (Digest::start(), Digest::start());
    x.fold_u64(1);
    x.fold_u64(2);
    y.fold_u64(2);
    y.fold_u64(1);
    if x.value() == y.value() {
        return Err(String::from("the digest ignores the order of folded words"));
    }
    Ok(())
}

#[test]
fn a_single_flipped_bit_changes_the_digest_of_every_suite_input() -> Result<(), String> {
    // sensitivity: perturb one operand of the Mag grid (one mantissa bit) and the mag digest must move
    let base = mag_suite()?;
    let (ops, _) = operands()?;
    let mut d = Digest::start();
    for a in ops.iter() {
        fold_mag(&mut d, &Mag::add_up(a, a));
    }
    let mut e = Digest::start();
    for a in ops.iter() {
        fold_mag(&mut e, &Mag::sub_lower(a, a));
    }
    if d.value() == e.value() {
        return Err(String::from(
            "add_up and sub_lower give the same digest over the grid",
        ));
    }
    let mut again = mag_suite()?;
    expect_eq("the mag suite is repeatable", &again.value(), &base.value())?;
    again.fold_u64(0);
    if again.value() == base.value() {
        return Err(String::from(
            "folding one more word left the digest unchanged",
        ));
    }
    Ok(())
}

#[test]
fn suites_are_independent_of_the_thread_that_runs_them() -> Result<(), String> {
    let seq = rows_sequential()?;
    let rows = std::thread::scope(|s| {
        let a = s.spawn(|| splitmix_suite().map(|d| d.row("splitmix")));
        let b = s.spawn(|| generators_suite().map(|d| d.row("generators")));
        let c = s.spawn(|| mag_suite().map(|d| d.row("mag_grid")));
        let d = s.spawn(|| ball_rows_suite().map(|d| d.row("ball_rows")));
        let mut out = String::new();
        for h in [a, b, c, d] {
            match h.join() {
                Ok(r) => out.push_str(&r?),
                Err(_) => return Err(String::from("a suite thread panicked")),
            }
        }
        Ok(out)
    })?;
    expect_eq("4 threads vs 1", &rows, &seq)
}

#[test]
fn grid_digests_match_the_committed_expected_digests() -> Result<(), String> {
    let got = rows_sequential()?;
    let want = expected_rows()?;
    if got == want {
        Ok(())
    } else {
        Err(format!(
            "grid digests drifted (update tools/digest_expected.txt only for a reviewed reason).\nactual rows:\n{}expected rows:\n{}",
            got, want
        ))
    }
}

#[test]
fn grid_digests_write_the_digest_file_for_the_script() -> Result<(), String> {
    let rows = rows_sequential()?;
    let path = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("grid_digests.txt");
    let text = format!(
        "# zenith-bound grid digests (S13); arch={} pointer_width={} (these two lines are not compared)\n# columns: SUITE WORDS FNV1A64\n{}",
        std::env::consts::ARCH,
        usize::BITS,
        rows
    );
    std::fs::write(&path, &text).map_err(|e| format!("{:?}: {}", path, e))?;
    let back = std::fs::read_to_string(&path).map_err(|e| format!("{:?}: {}", path, e))?;
    expect_eq("file text", &back, &text)
}
