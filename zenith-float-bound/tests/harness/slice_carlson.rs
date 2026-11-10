//! Slice S4b harness rows (`src/carlson.rs`). Owned by slice S4b (Cassie) alone.
//! CONTRACT 16.S4b.1 to 16.S4b.3 (red markers LW.1 to LW.3): Lambert W on `Ball`, both real branches.
//! Oracle: tests/goldens/slices/carlson/lambertw_rows.txt (scripts/gen_carlson_goldens.py: python-flint arb brackets,
//! mpmath.lambertw checked inside each), plus hand rows for bad `prec`, out-of-range radii and call-order determinism.
//! COVERS: Ball::lambert_w0, Ball::lambert_wm1

use std::fs;

use crate::common::goldens::golden_dir;
use zenith_float_bound::ball::Ball;
use zenith_float_bound::enclosure::Enclosure;
use zenith_float_bound::error::BoundError;
use zenith_float_bound::mag::Mag;
use zenith_float_bound::policy::DomainPolicy;
use zenith_float_bound::tri::Tri;
use zenith_float_bound::MAX_PRECISION;
use zenith_float_num::{Consts, ExactNum, RoundingMode, EXPONENT_MAX, EXPONENT_MIN};

const TOP: u64 = 0x8000_0000_0000_0000;
const GOLD_ROWS: usize = 680;

fn num(v: i64) -> ExactNum {
    ExactNum::from(v)
}

fn mag(mant: u64, exp: i64) -> Result<Mag, String> {
    Mag::finite(mant, exp).map_err(|e| format!("Mag::finite({:#x}, {}): {:?}", mant, exp, e))
}

fn ball(m: ExactNum, r: Mag, prec: usize) -> Result<Ball, String> {
    Ball::from_mid_rad(m, r, prec).map_err(|e| format!("Ball::from_mid_rad: {:?}", e))
}

/// `SIGN EXP MANT` (scripts/gen_carlson_goldens.py): value = (-1)^SIGN · MANT · 2^(EXP − bitlen(MANT)), MANT hex.
fn gold_val(sign: &str, exp: &str, hex: &str) -> Result<ExactNum, String> {
    let e: i64 = exp.parse().map_err(|_| format!("bad exponent {:?}", exp))?;
    let first = hex
        .chars()
        .next()
        .ok_or("empty mantissa")?
        .to_digit(16)
        .ok_or("bad digit")?;
    if first == 0 {
        return if hex == "0" && sign == "0" && e == 0 {
            Ok(num(0))
        } else {
            Err(format!("non-canonical zero {} {} {}", sign, exp, hex))
        };
    }
    let len = i64::try_from(hex.len()).map_err(|_| "mantissa length")?;
    let bitlen = 4 * (len - 1) + i64::from(32 - first.leading_zeros());
    let base = e - bitlen;
    let mut acc = num(0);
    let digits = hex.as_bytes();
    let mut shift = 0i64;
    // Least significant 16 hex digits (64 bits) first.
    for bytes in digits.rchunks(16) {
        let chunk = std::str::from_utf8(bytes).map_err(|_| "utf8")?;
        let c = u64::from_str_radix(chunk, 16).map_err(|_| format!("bad hex {:?}", chunk))?;
        if c != 0 {
            let lz = i64::from(c.leading_zeros());
            let part = ExactNum::from_mag_parts(c << lz, base + shift + 64 - lz, 64)
                .map_err(|e| format!("from_mag_parts: {:?}", e))?;
            acc = acc.add_full_prec(&part);
        }
        shift += 64;
    }
    if sign == "0" {
        Ok(acc)
    } else if sign == "1" {
        Ok(num(0).sub_full_prec(&acc))
    } else {
        Err(format!("bad sign {:?}", sign))
    }
}

fn take_val(f: &[&str], i: &mut usize) -> Result<ExactNum, String> {
    let v = gold_val(
        f.get(*i).ok_or("short row")?,
        f.get(*i + 1).ok_or("short row")?,
        f.get(*i + 2).ok_or("short row")?,
    )?;
    *i += 3;
    Ok(v)
}

fn want_contains(what: &str, b: &Ball, x: &ExactNum) -> Result<(), String> {
    match b.contains(x) {
        Tri::True => Ok(()),
        t @ (Tri::False | Tri::Unknown) => Err(format!(
            "{}: {:?} does not contain {:?} ({:?})",
            what, b, x, t
        )),
    }
}

#[derive(Debug, PartialEq, Eq)]
enum Branch {
    W0,
    Wm1,
}

fn call(br: &Branch, x: &Ball, prec: usize, policy: &DomainPolicy) -> Result<Ball, BoundError> {
    match br {
        Branch::W0 => x.lambert_w0(prec, policy),
        Branch::Wm1 => x.lambert_wm1(prec, policy),
    }
}

enum Want {
    Ok {
        llo: ExactNum,
        hhi: ExactNum,
        tb: ExactNum,
    },
    Domain,
    AllReals,
}

struct Row {
    line: usize,
    branch: Branch,
    policy: DomainPolicy,
    prec: usize,
    x: Ball,
    want: Want,
}

fn parse_row(line: usize, text: &str) -> Result<Row, String> {
    let f: Vec<&str> = text.split_whitespace().collect();
    let branch = match f.first() {
        Some(&"w0") => Branch::W0,
        Some(&"wm1") => Branch::Wm1,
        other => return Err(format!("line {}: branch {:?}", line, other)),
    };
    let policy = match f.get(1) {
        Some(&"strict") => DomainPolicy::Strict,
        Some(&"clip") => DomainPolicy::Clip,
        other => return Err(format!("line {}: policy {:?}", line, other)),
    };
    let prec: usize = f
        .get(2)
        .ok_or("short row")?
        .parse()
        .map_err(|_| format!("line {}: prec", line))?;
    let mut i = 3;
    let xm = take_val(&f, &mut i)?;
    let x = if f.get(i) == Some(&"inf") {
        i += 1;
        Ball::all_reals(prec).map_err(|e| format!("line {}: all reals {:?}", line, e))?
    } else {
        let xr = take_val(&f, &mut i)?;
        ball(xm, Mag::from_exact_num(&xr), prec)?
    };
    let want = match f.get(i) {
        Some(&"ok") => {
            i += 1;
            let llo = take_val(&f, &mut i)?;
            take_val(&f, &mut i)?; // LHI
            take_val(&f, &mut i)?; // HLO
            let hhi = take_val(&f, &mut i)?;
            let tb = take_val(&f, &mut i)?;
            Want::Ok { llo, hhi, tb }
        }
        Some(&"domain") => Want::Domain,
        Some(&"allreals") => Want::AllReals,
        other => return Err(format!("line {}: outcome {:?}", line, other)),
    };
    Ok(Row {
        line,
        branch,
        policy,
        prec,
        x,
        want,
    })
}

fn gold_rows() -> Result<Vec<Row>, String> {
    let path = golden_dir()
        .join("slices")
        .join("carlson")
        .join("lambertw_rows.txt");
    let text = fs::read_to_string(&path).map_err(|e| format!("{}: {}", path.display(), e))?;
    let mut rows = Vec::new();
    for (n, l) in text.lines().enumerate() {
        if l.starts_with('#') || l.trim().is_empty() {
            continue;
        }
        rows.push(parse_row(n + 1, l)?);
    }
    if rows.len() != GOLD_ROWS {
        return Err(format!("{} gold rows, expected {}", rows.len(), GOLD_ROWS));
    }
    Ok(rows)
}

/// `Err((message, fatal))`; `fatal` when the call hit `Resource`, which no gold row expects.
fn check_row(r: &Row) -> Result<(), (String, bool)> {
    let got = call(&r.branch, &r.x, r.prec, &r.policy);
    let fatal = matches!(got, Err(BoundError::Resource { .. }));
    check_got(r, got).map_err(|e| (e, fatal))
}

fn check_got(r: &Row, got: Result<Ball, BoundError>) -> Result<(), String> {
    let what = format!(
        "line {} {:?} {:?} prec {}",
        r.line, r.branch, r.policy, r.prec
    );
    match (&r.want, got) {
        (Want::Domain, Err(BoundError::Domain { .. })) => Ok(()),
        (Want::AllReals, Ok(b)) if !b.is_informative() && b.prec() == r.prec => Ok(()),
        (Want::Ok { llo, hhi, tb }, Ok(b)) => {
            if b.prec() != r.prec {
                return Err(format!("{}: prec {} kept as {}", what, r.prec, b.prec()));
            }
            want_contains(&format!("{} LLO", what), &b, llo)?;
            want_contains(&format!("{} HHI", what), &b, hhi)?;
            let bound = Mag::from_exact_num(tb);
            if b.rad() > &bound {
                return Err(format!(
                    "{}: rad {:?} above TB {:?} (16.S4b.3)",
                    what,
                    b.rad(),
                    bound
                ));
            }
            Ok(())
        }
        (Want::Domain, got) => Err(format!("{}: want Domain, got {:?}", what, got)),
        (Want::AllReals, got) => Err(format!("{}: want all reals, got {:?}", what, got)),
        (Want::Ok { .. }, got) => Err(format!("{}: want a ball, got {:?}", what, got)),
    }
}

/// LW.3: every gold row (outcome, `prec`, LLO and HHI contained, rad <= TB).
/// Stops at the first `Resource` error: each costs a run to the precision cap, and one is already a failure.
#[test]
fn lambertw_gold_rows() -> Result<(), String> {
    let rows = gold_rows()?;
    let mut bad = Vec::new();
    for r in &rows {
        if let Err((e, fatal)) = check_row(r) {
            bad.push(e);
            if fatal {
                bad.push(format!(
                    "stopped at line {}: Resource is never a gold outcome",
                    r.line
                ));
                break;
            }
        }
    }
    if bad.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "{} of {} rows failed:\n{}",
            bad.len(),
            rows.len(),
            bad.join("\n")
        ))
    }
}

/// LW.1: `InvalidArgument` for `prec` 0 and `MAX_PRECISION + 1`, both branches and policies, before any domain check.
#[test]
fn lambertw_bad_prec() -> Result<(), String> {
    let x = ball(num(-5), Mag::zero(), 53)?;
    for prec in [0, MAX_PRECISION + 1] {
        for br in [Branch::W0, Branch::Wm1] {
            for pol in [DomainPolicy::Strict, DomainPolicy::Clip] {
                match call(&br, &x, prec, &pol) {
                    Err(BoundError::InvalidArgument { .. }) => {}
                    got => return Err(format!("{:?} prec {} {:?}: {:?}", br, prec, pol, got)),
                }
            }
        }
    }
    // The public methods themselves, not only through `call`.
    match (
        x.lambert_w0(0, &DomainPolicy::Strict),
        x.lambert_wm1(0, &DomainPolicy::Clip),
    ) {
        (Err(BoundError::InvalidArgument { .. }), Err(BoundError::InvalidArgument { .. })) => {
            Ok(())
        }
        got => Err(format!("direct prec 0: {:?}", got)),
    }
}

/// `n / d` rounded with `rm` at 128 bits.
fn ratio(n: i64, d: i64, rm: RoundingMode) -> ExactNum {
    num(n).div(&num(d), 128, rm)
}

/// LW.2: a radius above zf's range is `Domain` (Strict) or all reals (Clip); one below it still gives a ball that sits
/// inside a 10-digit bracket of the point value.
#[test]
fn lambertw_out_of_range_radius() -> Result<(), String> {
    for br in [Branch::W0, Branch::Wm1] {
        let x = ball(num(0), mag(TOP, EXPONENT_MAX + 100)?, 53)?;
        match call(&br, &x, 53, &DomainPolicy::Strict) {
            Err(BoundError::Domain { .. }) => {}
            got => return Err(format!("{:?} huge radius Strict: {:?}", br, got)),
        }
        match call(&br, &x, 53, &DomainPolicy::Clip) {
            Ok(b) if !b.is_informative() => {}
            got => return Err(format!("{:?} huge radius Clip: {:?}", br, got)),
        }
    }
    // W0(1) = 0.5671432904097..., W-1(-1/4) = -2.1532923641103...
    let quarter = ExactNum::from_mag_parts(TOP, -1, 64).map_err(|e| format!("{:?}", e))?;
    let cases = [
        (Branch::W0, num(1), (5671432904, 5671432905)),
        (
            Branch::Wm1,
            num(0).sub_full_prec(&quarter),
            (-21532923642, -21532923641),
        ),
    ];
    for (br, m, (lo, hi)) in cases {
        let x = ball(m, mag(TOP, EXPONENT_MIN - 200)?, 53)?;
        let b = call(&br, &x, 53, &DomainPolicy::Strict)
            .map_err(|e| format!("{:?} tiny radius: {:?}", br, e))?;
        let lo = ratio(lo, 10_000_000_000, RoundingMode::Down);
        let hi = ratio(hi, 10_000_000_000, RoundingMode::Up);
        if !b.is_informative() || b.contains(&lo) != Tri::False || b.contains(&hi) != Tri::False {
            return Err(format!(
                "{:?} tiny radius: {:?} not inside ({:?}, {:?})",
                br, b, lo, hi
            ));
        }
        let mid = b.mid();
        if mid.compare(&lo).is_none_or(|c| c <= 0) || mid.compare(&hi).is_none_or(|c| c >= 0) {
            return Err(format!(
                "{:?} tiny radius: mid {:?} outside ({:?}, {:?})",
                br, mid, lo, hi
            ));
        }
    }
    Ok(())
}

/// LW.1: results are a function of the inputs only: identical in forward and reverse call order.
/// No gold row may give `Resource` (it fails at once rather than running every row to the cap).
#[test]
fn lambertw_call_order_determinism() -> Result<(), String> {
    let rows: Vec<Row> = gold_rows()?.into_iter().filter(|r| r.prec <= 256).collect();
    let run = |r: &Row| match call(&r.branch, &r.x, r.prec, &r.policy) {
        Err(e @ BoundError::Resource { .. }) => Err(format!("line {}: {:?}", r.line, e)),
        got => Ok(got
            .map(|b| format!("{:?}", b))
            .map_err(|e| format!("{:?}", e))),
    };
    let mut fwd = Vec::new();
    for r in &rows {
        fwd.push(run(r)?);
    }
    if !fwd.iter().any(|r| r.is_ok()) {
        return Err(format!(
            "no row of {} gave a ball: {:?}",
            rows.len(),
            fwd.first()
        ));
    }
    let mut rev = Vec::new();
    for r in rows.iter().rev() {
        rev.push(run(r)?);
    }
    rev.reverse();
    for ((a, b), r) in fwd.iter().zip(rev.iter()).zip(rows.iter()) {
        if a != b {
            return Err(format!("line {}: forward {:?} reverse {:?}", r.line, a, b));
        }
    }
    Ok(())
}

/// LW.1: an input still undecided against `-1/e` at the working-precision cap counts as outside. `x` is `-1/e` rounded
/// towards zero at 20000 bits (truly inside, but `|x + 1/e| < 2^-19990` is far below the `prec = 53` cap of 4520 bits).
#[test]
fn lambertw_undecided_counts_outside() -> Result<(), String> {
    let mut cc = Consts::new().map_err(|e| format!("Consts::new: {:?}", e))?;
    let e_hi = cc.e(20_100, RoundingMode::Up);
    let x = num(-1).div(&e_hi, 20_000, RoundingMode::Up);
    let b = ball(x, Mag::zero(), 53)?;
    for pol in [DomainPolicy::Strict, DomainPolicy::Clip] {
        match (b.lambert_w0(53, &pol), b.lambert_wm1(53, &pol)) {
            (Err(BoundError::Domain { .. }), Err(BoundError::Domain { .. })) => {}
            got => return Err(format!("{:?}: {:?}", pol, got)),
        }
    }
    Ok(())
}

/// LW.2: `W0` under `Clip` on a ball straddling `-1/e` has lower end `W0(-1/e) = -1`, so the result must contain `-1`,
/// also when `mid + 1` needs more than 64 bits (prec 100 to 300). Mutant-killer for 626-RADBELOW64: it passes on
/// the code at 41df22b.
#[test]
fn lambertw_w0_clip_contains_minus_one() -> Result<(), String> {
    let mut bad = Vec::new();
    // source: W0(-1/e) = -1, DLMF 4.13 (https://dlmf.nist.gov/4.13).
    // no-source: hand-chosen inputs: mid -0.360..-0.367, r = 0.01, so each ball straddles -1/e ≈ -0.367879.
    for prec in [100, 200, 300] {
        for j in 0..8 {
            let r = Mag::from_exact_num(&ratio(1, 100, RoundingMode::Up));
            let x = ball(ratio(-360 - j, 1000, RoundingMode::Down), r, prec)?;
            match x.lambert_w0(prec, &DomainPolicy::Clip) {
                Ok(b) if b.contains(&num(-1)) == Tri::True => {}
                got => bad.push(format!("prec {} m -0.{}: {:?}", prec, 360 + j, got)),
            }
        }
    }
    if bad.is_empty() {
        Ok(())
    } else {
        Err(format!("{} of 24 miss -1:\n{}", bad.len(), bad.join("\n")))
    }
}
