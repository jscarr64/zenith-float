//! Slice LA2 harness rows (`src/lu.rs`). Owned by slice LA2 alone (Cassie).
//! Exact golds: Cassie (owns the zb LA exact golds, 2026-10-24).
//! Goldens: tests/goldens/slices/lu/ (its own SHA256SUMS and GENERATORS). Mutants: tools/mutants/manifest.d/lu.txt, ids 1100..=1199.
//! Golden file headers for this slice are declared here, not in tests/common/goldens.rs.
//! CONTRACT 16.LA2.1 to 16.LA2.7, 16.LA1.4, 16.LA1.6, 16.LA1.7, 16.LA.2 to 16.LA.6, LA4.2. Every replay runs at the precisions
//! {53, 113, 192} (16.LA1.7; the gold files have no prec column). The logdet replay covers every gate row (the `# rows:` header) of
//! tests/goldens/slices/cholesky/la_logdet.txt, not only cosmo_*. (LA4.2 names slice_cholesky.rs for the logdet replay;
//! the logdet API is LA2's, so its replay lives here.)
//! COVERS: BallMat::lu, BallMat::logdet, BallMat::solve, BallMat::solve_mat, BallLu::perm, BallLu::l, BallLu::u, BallLu::perm_sign, BallLu::det, BallLu::solve, LogDet::sign, LogDet::det, LogDet::ln_abs_det

use crate::common::goldens::golden_dir;
use zenith_float_bound::array::{BallMat, BallVec};
use zenith_float_bound::ball::Ball;
use zenith_float_bound::enclosure::Enclosure;
use zenith_float_bound::error::BoundError;
use zenith_float_bound::lu::BallLu;
use zenith_float_bound::mag::Mag;
use zenith_float_bound::tri::Tri;
use zenith_float_num::{ExactNum, EXPONENT_MAX};

/// The test precisions of 16.LA1.7, used by every replay in this file.
const PRECS: [usize; 3] = [53, 113, 192];
/// `Mag`/`ExactNum` mantissa of one half.
const TOP: u64 = 0x8000_0000_0000_0000;
/// 16.LA.4: solve rows certified with exact arithmetic (they must return the named error, never `undecided_*`).
const SOLVE_CERTIFIED: [&str; 4] = [
    "solve.zero.02",
    "solve.zero.05",
    "solve.onebyone.zero",
    "solve.singular.dup_rows3",
];

// ---- helpers ------------------------------------------------------------------------------------------------------------
fn ok<T>(what: &str, r: Result<T, BoundError>) -> Result<T, String> {
    r.map_err(|e| format!("{what}: unexpected {e:?}"))
}

fn num(v: i64) -> ExactNum {
    ExactNum::from(v)
}

fn plus_zero(x: &ExactNum) -> Result<ExactNum, String> {
    Ok(x.add_full_prec(&num(0)))
}

/// `2^k` exactly.
fn pow2(k: i64) -> Result<ExactNum, String> {
    let e = k.checked_add(1).ok_or("pow2 exponent")?;
    ExactNum::from_mag_parts(TOP, e, 64).map_err(|e| format!("pow2({k}): {e:?}"))
}

/// A signed decimal integer of any size, exactly (Horner in 18-digit chunks).
fn big_int(s: &str) -> Result<ExactNum, String> {
    let (neg, digits) = match s.strip_prefix('-') {
        Some(d) => (true, d),
        None => (false, s),
    };
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return Err(format!("bad integer {s:?}"));
    }
    let mut acc = num(0);
    let mut rest = digits;
    while !rest.is_empty() {
        let take = rest.len().min(18);
        let (head, tail) = rest.split_at(take);
        let chunk: i64 = head.parse().map_err(|e| format!("{s:?}: {e}"))?;
        let mut scale = 1_i64;
        for _ in 0..take {
            scale = scale.checked_mul(10).ok_or("scale overflow")?;
        }
        acc = acc.mul_full_prec(&num(scale)).add_full_prec(&num(chunk));
        rest = tail;
    }
    if neg {
        acc = -acc;
    }
    Ok(acc)
}

/// An exact rational `p/q` (q > 0) from a gold token.
struct Rat {
    p: ExactNum,
    q: ExactNum,
}

fn rat(tok: &str) -> Result<Rat, String> {
    let (p, q) = tok
        .split_once('/')
        .ok_or_else(|| format!("not p/q: {tok:?}"))?;
    let q = big_int(q)?;
    if q <= num(0) {
        return Err(format!("non-positive denominator in {tok:?}"));
    }
    Ok(Rat { p: big_int(p)?, q })
}

/// A dyadic `m*2^e` gold token, exactly.
fn dyadic(tok: &str) -> Result<ExactNum, String> {
    let (m, e) = tok
        .split_once("*2^")
        .ok_or_else(|| format!("not m*2^e: {tok:?}"))?;
    let e: i64 = e.parse().map_err(|x| format!("{tok:?}: {x}"))?;
    Ok(big_int(m)?.mul_full_prec(&pow2(e)?))
}

/// The input ball of an exact rational at `prec`: exact for integers, else S1 `p / q`.
fn input_ball(r: &Rat, prec: usize) -> Result<Ball, String> {
    let p = ok("exact p", Ball::exact(plus_zero(&r.p)?, prec))?;
    if r.q == num(1) {
        return Ok(p);
    }
    let q = ok("exact q", Ball::exact(plus_zero(&r.q)?, prec))?;
    ok("p / q", p.div(&q, prec))
}

/// Whether the ball contains the exact rational: `|m·q − p| ≤ r·q` (an `Infinity` radius contains everything).
fn holds(b: &Ball, r: &Rat) -> Result<bool, String> {
    if b.rad().is_infinity() {
        return Ok(true);
    }
    let rad = b.rad().to_exact_num().map_err(|e| format!("{e:?}"))?;
    let d = b.mid().mul_full_prec(&r.q).sub_full_prec(&r.p);
    let rq = rad.mul_full_prec(&r.q);
    Ok(d.abs_cmp(&rq).is_some_and(|c| c <= 0))
}

fn want_holds(what: &str, b: &Ball, r: &Rat) -> Result<(), String> {
    if holds(b, r)? {
        Ok(())
    } else {
        Err(format!(
            "{what}: {b:?} does not contain {:?}/{:?}",
            r.p, r.q
        ))
    }
}

fn is_exact_zero(b: &Ball) -> bool {
    b.mid().is_zero() && b.rad().is_zero()
}

fn has_zero(b: &Ball) -> bool {
    match b.contains(&num(0)) {
        Tri::True | Tri::Unknown => true,
        Tri::False => false,
    }
}

fn read_gold(dir: &str, name: &str) -> Result<String, String> {
    let path = golden_dir().join("slices").join(dir).join(name);
    std::fs::read_to_string(&path).map_err(|e| format!("cannot read {}: {e}", path.display()))
}

fn header_rows(text: &str) -> Result<usize, String> {
    let line = text
        .lines()
        .find_map(|l| l.strip_prefix("# rows: "))
        .ok_or("no '# rows:' header")?;
    line.trim()
        .parse::<usize>()
        .map_err(|e| format!("bad '# rows:' header {line:?}: {e}"))
}

/// Data rows, checked against the `# rows:` header.
fn data_rows(text: &str) -> Result<Vec<Vec<&str>>, String> {
    let rows: Vec<Vec<&str>> = text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| l.split_whitespace().collect())
        .collect();
    let want = header_rows(text)?;
    if rows.len() != want {
        return Err(format!("{} data rows, header says {want}", rows.len()));
    }
    Ok(rows)
}

/// A token cursor over one gold row.
struct Cur<'a> {
    toks: &'a [&'a str],
    at: usize,
}

impl<'a> Cur<'a> {
    fn tok(&mut self) -> Result<&'a str, String> {
        let t = self.toks.get(self.at).ok_or("row too short")?;
        self.at = self.at.checked_add(1).ok_or("cursor overflow")?;
        Ok(t)
    }
    fn rats(&mut self, n: usize) -> Result<Vec<Rat>, String> {
        let mut v = Vec::new();
        for _ in 0..n {
            v.push(rat(self.tok()?)?);
        }
        Ok(v)
    }
    fn rest_is(&self, want: &[&str]) -> bool {
        self.toks.get(self.at..) == Some(want)
    }
}

fn ball_mat(n: usize, m: usize, a: &[Rat], prec: usize) -> Result<BallMat, String> {
    let mut balls = Vec::new();
    for r in a {
        balls.push(input_ball(r, prec)?);
    }
    ok("BallMat::from_balls", BallMat::from_balls(n, m, balls))
}

fn ball_vec(b: &[Rat], prec: usize) -> Result<BallVec, String> {
    let mut balls = Vec::new();
    for r in b {
        balls.push(input_ball(r, prec)?);
    }
    Ok(BallVec::from_balls(balls))
}

fn entry(m: &BallMat, i: usize, j: usize) -> Result<&Ball, String> {
    m.get(i, j).ok_or_else(|| format!("no entry ({i},{j})"))
}

/// `perm_sign · Π sign(mid u_kk)`, or `None` when some `u_kk` contains 0.
fn lu_sign(lu: &BallLu) -> Result<Option<i8>, String> {
    let mut s = lu.perm_sign();
    for k in 0..lu.u().nrows() {
        let u = entry(lu.u(), k, k)?;
        if has_zero(u) {
            return Ok(None);
        }
        if u.mid().is_negative() {
            s = s.checked_neg().ok_or("sign overflow")?;
        }
    }
    Ok(Some(s))
}

fn parse_dims(tok: &str) -> Result<(usize, usize), String> {
    match tok.split_once('x') {
        Some((m, n)) => Ok((
            m.parse().map_err(|e| format!("dims {tok:?}: {e}"))?,
            n.parse().map_err(|e| format!("dims {tok:?}: {e}"))?,
        )),
        None => {
            let n: usize = tok.parse().map_err(|e| format!("dims {tok:?}: {e}"))?;
            Ok((n, n))
        }
    }
}

fn mul_dims(m: usize, n: usize) -> Result<usize, String> {
    m.checked_mul(n).ok_or_else(|| "dims overflow".to_string())
}

fn idx(i: usize, n: usize, j: usize) -> Result<usize, String> {
    mul_dims(i, n)?
        .checked_add(j)
        .ok_or_else(|| "index overflow".to_string())
}

/// The error must be `variant` with `detail` (an empty `detail` accepts any).
fn want_err_detail<T: core::fmt::Debug>(
    what: &str,
    r: Result<T, BoundError>,
    variant: &str,
    detail: &str,
) -> Result<(), String> {
    let (got_v, got_d) = match &r {
        Err(BoundError::Singular { detail, .. }) => ("Singular", *detail),
        Err(BoundError::Resource { detail, .. }) => ("Resource", *detail),
        Err(BoundError::DimensionMismatch { detail, .. }) => ("DimensionMismatch", *detail),
        Err(BoundError::InvalidArgument { detail, .. }) => ("InvalidArgument", *detail),
        Err(e) => return Err(format!("{what}: want {variant} {detail:?}, got {e:?}")),
        Ok(v) => return Err(format!("{what}: want {variant} {detail:?}, got Ok({v:?})")),
    };
    if got_v == variant && (detail.is_empty() || got_d == detail) {
        Ok(())
    } else {
        Err(format!("{what}: want {variant} {detail:?}, got {r:?}"))
    }
}

/// 16.LA.4: a refusal row passes on `Singular "singular"` or (outside the certified list) on `Resource undecided`.
fn want_refusal<T: core::fmt::Debug>(
    what: &str,
    r: Result<T, BoundError>,
    certified: bool,
    undecided: &str,
) -> Result<(), String> {
    match &r {
        Err(BoundError::Singular { detail, .. }) if *detail == "singular" => Ok(()),
        Err(BoundError::Resource { detail, .. }) if !certified && *detail == undecided => Ok(()),
        Err(_) | Ok(_) => Err(format!(
            "{what}: want Singular \"singular\"{}, got {r:?}",
            if certified {
                " (certified row, 16.LA.4)"
            } else {
                " or undecided"
            }
        )),
    }
}

// ---- la_lu.txt replay (16.LA2.6) ------------------------------------------------------------------------------------------
fn replay_lu_row(toks: &[&str], prec: usize) -> Result<(), String> {
    let id = toks.first().ok_or("empty row")?;
    let mut c = Cur { toks, at: 1 };
    let (m, n) = parse_dims(c.tok()?)?;
    let a_r = c.rats(mul_dims(m, n)?)?;
    let a = ball_mat(m, n, &a_r, prec)?;
    if m != n {
        if !c.rest_is(&["ERR", "DimensionMismatch", "rectangular"]) {
            return Err(format!("{id}: non-square row without its ERR tail"));
        }
        return want_err_detail(id, a.lu(prec), "DimensionMismatch", "rectangular");
    }
    let mut perm = Vec::new();
    for _ in 0..n {
        perm.push(
            c.tok()?
                .parse::<usize>()
                .map_err(|e| format!("{id}: perm {e}"))?,
        );
    }
    let l_r = c.rats(mul_dims(n, n)?)?;
    let u_r = c.rats(mul_dims(n, n)?)?;
    let lu = ok(id, a.lu(prec))?;
    if lu.l().nrows() != n || lu.u().ncols() != n || lu.perm().len() != n {
        return Err(format!("{id}: factor shapes wrong"));
    }
    if lu.perm() == perm.as_slice() {
        for i in 0..n {
            for j in 0..n {
                let k = idx(i, n, j)?;
                let (lr, ur) = (l_r.get(k).ok_or("L")?, u_r.get(k).ok_or("U")?);
                want_holds(&format!("{id}@{prec} L[{i},{j}]"), entry(lu.l(), i, j)?, lr)?;
                want_holds(&format!("{id}@{prec} U[{i},{j}]"), entry(lu.u(), i, j)?, ur)?;
            }
        }
        return Ok(());
    }
    // 16.LA2.6: a different tie break; P_zb·A (exact) must lie in the ball product L·U.
    eprintln!(
        "{id}@{prec}: perm {:?} differs from gold {perm:?}; checking P·A in L·U",
        lu.perm()
    );
    let prod = ok("matmul", lu.l().matmul(lu.u(), prec))?;
    for i in 0..n {
        let src = *lu.perm().get(i).ok_or("perm")?;
        for j in 0..n {
            let k = idx(src, n, j)?;
            want_holds(
                &format!("{id}@{prec} (LU)[{i},{j}]"),
                entry(&prod, i, j)?,
                a_r.get(k).ok_or("A")?,
            )?;
        }
    }
    Ok(())
}

#[test]
fn la2_lu_gate_rows_replay_through_zb() -> Result<(), String> {
    let text = read_gold("lu", "la_lu.txt")?;
    for prec in PRECS {
        for row in data_rows(&text)? {
            replay_lu_row(&row, prec)?;
        }
    }
    Ok(())
}

// ---- det, sign and logdet against la_logdet.txt (16.LA1.4, 16.LA1.6, 16.LA1.7, LA4.2) ---------------------------------------
struct LogdetRow<'a> {
    id: &'a str,
    a: Vec<Rat>,
    n: usize,
    /// `(sign, det, lo256, hi256)`, or `None` for `ERR Singular singular`.
    gold: Option<(i8, Rat, ExactNum, ExactNum)>,
}

fn logdet_row<'a>(toks: &'a [&'a str]) -> Result<LogdetRow<'a>, String> {
    let id = toks.first().ok_or("empty row")?;
    let mut c = Cur { toks, at: 1 };
    let (n, _) = parse_dims(c.tok()?)?;
    let a = c.rats(mul_dims(n, n)?)?;
    if c.rest_is(&["ERR", "Singular", "singular"]) {
        return Ok(LogdetRow {
            id,
            a,
            n,
            gold: None,
        });
    }
    let sign: i8 = c.tok()?.parse().map_err(|e| format!("{id}: sign {e}"))?;
    let det = rat(c.tok()?)?;
    let lo = dyadic(c.tok()?)?;
    let hi = dyadic(c.tok()?)?;
    Ok(LogdetRow {
        id,
        a,
        n,
        gold: Some((sign, det, lo, hi)),
    })
}

/// lu → det encloses the exact det; the sign `perm_sign · Π sign(u_kk)` is exact whenever every `u_kk` excludes 0
/// (16.LA2.3, 16.LA2.5, 16.LA1.4). Singular rows: det encloses 0, and is exactly `0 ± 0` when a column was skipped.
#[test]
fn la2_det_and_sign_enclose_every_logdet_gold_row() -> Result<(), String> {
    let text = read_gold("cholesky", "la_logdet.txt")?;
    let rows = data_rows(&text)?;
    let want = header_rows(&text)?;
    if rows.len() != want {
        return Err(format!(
            "want all {want} logdet rows (`# rows:` header), got {}",
            rows.len()
        ));
    }
    for prec in PRECS {
        for toks in &rows {
            let r = logdet_row(toks)?;
            let a = ball_mat(r.n, r.n, &r.a, prec)?;
            let lu = ok(r.id, a.lu(prec))?;
            let det = ok(r.id, lu.det(prec))?;
            let what = format!("{}@{prec} det", r.id);
            match &r.gold {
                Some((sign, d, _, _)) => {
                    want_holds(&what, &det, d)?;
                    if let Some(s) = lu_sign(&lu)? {
                        if s != *sign {
                            return Err(format!("{what}: sign {s}, gold {sign}"));
                        }
                    }
                }
                None => {
                    if !has_zero(&det) {
                        return Err(format!("{what}: singular row, det {det:?} excludes 0"));
                    }
                    let skipped = (0..r.n).any(|k| lu.u().get(k, k).is_some_and(is_exact_zero));
                    if skipped && !is_exact_zero(&det) {
                        return Err(format!("{what}: a column was skipped but det is {det:?}"));
                    }
                }
            }
        }
    }
    Ok(())
}

/// Refusal rows of `logdet` (16.LA1.7, 16.LA.4): `Singular "singular"` or `Resource "undecided_det"`; a value row whose
/// `u_kk` straddles 0 gives `undecided_det`; every other value row is decided and, until S3 lands `ln`, returns
/// `Resource "needs_ln"` (16.LA1.6), or once ln exists a value with the gold sign whose det encloses the gold.
#[test]
fn la2_logdet_refusals_and_needs_ln_on_every_gold_row() -> Result<(), String> {
    let text = read_gold("cholesky", "la_logdet.txt")?;
    for prec in PRECS {
        for toks in &data_rows(&text)? {
            let r = logdet_row(toks)?;
            let what = format!("{}@{prec} logdet", r.id);
            let a = ball_mat(r.n, r.n, &r.a, prec)?;
            let lu = ok(r.id, a.lu(prec))?;
            let got = a.logdet(prec);
            match &r.gold {
                None => want_refusal(&what, got, false, "undecided_det")?,
                Some((sign, d, _, _)) => match (lu_sign(&lu)?, got) {
                    (None, got) => want_err_detail(&what, got, "Resource", "undecided_det")?,
                    (
                        Some(_),
                        Err(BoundError::Resource {
                            detail: "needs_ln", ..
                        }),
                    ) => {}
                    (Some(_), Ok(v)) => {
                        if v.sign() != *sign {
                            return Err(format!("{what}: sign {}, gold {sign}", v.sign()));
                        }
                        want_holds(&what, v.det(), d)?;
                    }
                    (Some(_), Err(e)) => return Err(format!("{what}: decided row gave {e:?}")),
                },
            }
        }
    }
    Ok(())
}

/// Heavy `enc4096`/`enc16384` enclosure of a row from `$ZB_LA_HEAVY_DIR/la_logdet_enc<b>_heavy.txt` (last two tokens lo hi).
fn heavy_enclosure(id: &str, bits: u32) -> Result<Option<(ExactNum, ExactNum)>, String> {
    let Ok(dir) = std::env::var("ZB_LA_HEAVY_DIR") else {
        return Ok(None);
    };
    let path = std::path::Path::new(&dir).join(format!("la_logdet_enc{bits}_heavy.txt"));
    let Ok(text) = std::fs::read_to_string(&path) else {
        return Ok(None);
    };
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.first() == Some(&id) {
            let n = toks.len();
            let lo = dyadic(toks.get(n.checked_sub(2).ok_or("short")?).ok_or("lo")?)?;
            let hi = dyadic(toks.get(n.checked_sub(1).ok_or("short")?).ok_or("hi")?)?;
            return Ok(Some((lo, hi)));
        }
    }
    Ok(None)
}

/// `(covers, overlaps)` of the ball against `[lo, hi]`.
fn vs_interval(b: &Ball, lo: &ExactNum, hi: &ExactNum) -> Result<(bool, bool), String> {
    if b.rad().is_infinity() {
        return Ok((true, true));
    }
    let r = b.rad().to_exact_num().map_err(|e| format!("{e:?}"))?;
    let m = b.mid();
    let d_lo = m.sub_full_prec(lo);
    let d_hi = hi.sub_full_prec(m);
    let covers = d_lo <= r && d_hi <= r;
    let neg_r = -plus_zero(&r)?;
    let overlaps = d_lo >= neg_r && d_hi >= neg_r;
    Ok((covers, overlaps))
}

/// The full LA4.2 replay: sign exact, det encloses det, `ln_abs_det` contains the exact ln|det| (16.LA.5 containment).
#[test]
#[ignore = "red: not implemented: CONTRACT LA4.2: ln_abs_det needs S3 ln (16.LA1.6: logdet returns Resource needs_ln until S3 merges); un-ignore when src/elementary.rs exports ln"]
fn la2_logdet_gate_rows_replay_through_zb() -> Result<(), String> {
    let text = read_gold("cholesky", "la_logdet.txt")?;
    for prec in PRECS {
        for toks in &data_rows(&text)? {
            let r = logdet_row(toks)?;
            let what = format!("{}@{prec} logdet", r.id);
            let a = ball_mat(r.n, r.n, &r.a, prec)?;
            let got = a.logdet(prec);
            let Some((sign, d, lo, hi)) = &r.gold else {
                want_refusal(&what, got, false, "undecided_det")?;
                continue;
            };
            let v = match got {
                Ok(v) => v,
                Err(BoundError::Resource {
                    detail: "undecided_det",
                    ..
                }) => continue,
                Err(e) => return Err(format!("{what}: {e:?}")),
            };
            if v.sign() != *sign {
                return Err(format!("{what}: sign {}, gold {sign}", v.sign()));
            }
            want_holds(&what, v.det(), d)?;
            if !v.ln_abs_det().is_informative() {
                return Err(format!(
                    "{what}: ln_abs_det is all-reals (must stay finite, 16.LA1.7)"
                ));
            }
            let (covers, overlaps) = vs_interval(v.ln_abs_det(), lo, hi)?;
            if covers {
                continue;
            }
            if !overlaps {
                return Err(format!(
                    "{what}: ln_abs_det {:?} disjoint from [lo256, hi256]",
                    v.ln_abs_det()
                ));
            }
            let mut decided = false;
            for bits in [4096_u32, 16384] {
                if let Some((hlo, hhi)) = heavy_enclosure(r.id, bits)? {
                    if vs_interval(v.ln_abs_det(), &hlo, &hhi)?.0 {
                        decided = true;
                    }
                }
            }
            if !decided {
                return Err(format!("{what}: partial overlap not re-decided by enc4096/enc16384 (set ZB_LA_HEAVY_DIR)"));
            }
        }
    }
    Ok(())
}

/// Gold format from the `# format: N` header (format 2 adds the six cond columns of scripts/la_cond_cols.py); none = 1.
fn gold_format(text: &str) -> Result<u32, String> {
    match text.lines().find_map(|l| l.strip_prefix("# format: ")) {
        None => Ok(1),
        Some(rest) => rest
            .split_whitespace()
            .next()
            .ok_or("empty '# format:' header")?
            .parse::<u32>()
            .map_err(|e| format!("bad '# format:' header {rest:?}: {e}")),
    }
}

/// The format-2 cond columns `cond_row k_row scale_row containment_only reason rref` of one non-refusal row (16.LA.5).
struct CondCols {
    k_row: usize,
    scale_row: Rat,
    containment_only: bool,
}

/// Parses the six cond columns at `toks[at..]` and cross-checks `containment_only` against `k_row + 8 >= max(PRECS)`.
fn cond_cols(id: &str, toks: &[&str], at: usize) -> Result<CondCols, String> {
    let rest = toks
        .get(at..)
        .ok_or_else(|| format!("{id}: row too short for cond columns"))?;
    let [cond_row, k_row, scale_row, flag, reason, rref] = rest else {
        return Err(format!("{id}: want 6 cond columns, got {}", rest.len()));
    };
    rat(cond_row)?;
    rat(rref)?;
    let k_row: usize = k_row
        .parse()
        .map_err(|e| format!("{id}: k_row {k_row:?}: {e}"))?;
    let containment_only = match (*flag, *reason) {
        ("0", "tight") => false,
        ("1", "k_row_too_large") => true,
        (f, r) => return Err(format!("{id}: bad containment_only/reason {f} {r}")),
    };
    let top = PRECS.iter().max().ok_or("PRECS")?;
    if containment_only != (k_row.checked_add(8).ok_or("k_row overflows")? >= *top) {
        return Err(format!(
            "{id}: containment_only {flag} disagrees with k_row {k_row}"
        ));
    }
    Ok(CondCols {
        k_row,
        scale_row: rat(scale_row)?,
        containment_only,
    })
}

/// 16.LA.5: `rad <= 2^(k_row - prec)·max(|v|, scale_row)` when `prec > k_row + 8` and the row is not containment-only;
/// `|v|` is read as `|mid|`. Returns whether the bound was enforced.
fn want_tight(what: &str, b: &Ball, c: &CondCols, prec: usize) -> Result<bool, String> {
    if c.containment_only || prec <= c.k_row.checked_add(8).ok_or("k_row overflows")? {
        return Ok(false);
    }
    if b.rad().is_infinity() {
        return Err(format!("{what}: infinite radius under the k_row bound"));
    }
    let k = i64::try_from(c.k_row).map_err(|e| format!("{e}"))?;
    let p = i64::try_from(prec).map_err(|e| format!("{e}"))?;
    let factor = pow2(k.checked_sub(p).ok_or("k_row - prec")?)?;
    let q = &c.scale_row.q;
    let mut m = b.mid().mul_full_prec(q);
    if m.is_negative() {
        m = -m;
    }
    let mut s = plus_zero(&c.scale_row.p)?;
    if s.is_negative() {
        s = -s;
    }
    let big = if m >= s { m } else { s };
    let lhs = b
        .rad()
        .to_exact_num()
        .map_err(|e| format!("{e:?}"))?
        .mul_full_prec(q);
    if lhs <= factor.mul_full_prec(&big) {
        Ok(true)
    } else {
        Err(format!(
            "{what}: rad {:?} over the k_row {} bound at prec {prec}",
            b.rad(),
            c.k_row
        ))
    }
}

/// 16.LA.5 tightness on every format-2 logdet row (`ln_abs_det`, `scale_row = 1`) and solve row (`x`, `scale_row = ‖x‖_∞`).
/// Red until the golds are format 2 (harper/zb-la-golds), S3 `ln` lands (16.LA1.6) and 16.S6.12 trisolve lands.
#[test]
#[ignore = "red: CONTRACT 16.LA.5: needs format-2 la_logdet.txt/la_solve.txt (cond_row k_row scale_row containment_only reason rref, harper/zb-la-golds), S3 ln (16.LA1.6) and 16.S6.12 trisolve; un-ignore when all three are on the branch"]
fn la2_tightness_k_row_bound_on_logdet_and_solve_rows() -> Result<(), String> {
    let ld = read_gold("cholesky", "la_logdet.txt")?;
    let sv = read_gold("lu", "la_solve.txt")?;
    for (name, text) in [("la_logdet.txt", &ld), ("la_solve.txt", &sv)] {
        if gold_format(text)? != 2 {
            return Err(format!(
                "{name} is format {}: no cond columns (16.LA.5)",
                gold_format(text)?
            ));
        }
    }
    let mut enforced = 0_usize;
    for prec in PRECS {
        for toks in &data_rows(&ld)? {
            let r = logdet_row(toks)?;
            if r.gold.is_none() {
                continue;
            }
            let at = mul_dims(r.n, r.n)?.checked_add(6).ok_or("width")?;
            let c = cond_cols(r.id, toks, at)?;
            if c.scale_row.p != num(1) || c.scale_row.q != num(1) {
                return Err(format!("{}: logdet scale_row must be 1/1", r.id));
            }
            let what = format!("{}@{prec} ln_abs_det", r.id);
            let a = ball_mat(r.n, r.n, &r.a, prec)?;
            let v = match a.logdet(prec) {
                Ok(v) => v,
                Err(BoundError::Resource {
                    detail: "undecided_det",
                    ..
                }) => continue,
                Err(e) => return Err(format!("{what}: {e:?}")),
            };
            if want_tight(&what, v.ln_abs_det(), &c, prec)? {
                enforced = enforced.checked_add(1).ok_or("count")?;
            }
        }
        for toks in &data_rows(&sv)? {
            let r = solve_row(toks)?;
            if r.x.is_none() {
                continue;
            }
            let at = mul_dims(r.n, r.n)?
                .checked_add(r.n.checked_mul(2).ok_or("width")?)
                .and_then(|w| w.checked_add(2))
                .ok_or("width")?;
            let c = cond_cols(r.id, toks, at)?;
            let what = format!("{}@{prec} solve", r.id);
            let a = ball_mat(r.n, r.n, &r.a, prec)?;
            let x = match a.solve(&ball_vec(&r.b, prec)?, prec) {
                Ok(x) => x,
                Err(BoundError::Resource { detail, .. }) if detail.starts_with("undecided") => {
                    continue
                }
                Err(e) => return Err(format!("{what}: {e:?}")),
            };
            for i in 0..r.n {
                if want_tight(&format!("{what} x[{i}]"), x.get(i).ok_or("x")?, &c, prec)? {
                    enforced = enforced.checked_add(1).ok_or("count")?;
                }
            }
        }
    }
    if enforced == 0 {
        return Err("no row enforced the k_row bound (16.LA.5)".to_string());
    }
    Ok(())
}

// ---- la_solve.txt (16.LA2.3, 16.LA2.4) ----------------------------------------------------------------------------------------
struct SolveRow<'a> {
    id: &'a str,
    n: usize,
    a: Vec<Rat>,
    b: Vec<Rat>,
    x: Option<Vec<Rat>>,
}

fn solve_row<'a>(toks: &'a [&'a str]) -> Result<SolveRow<'a>, String> {
    let id = toks.first().ok_or("empty row")?;
    let mut c = Cur { toks, at: 1 };
    let (n, _) = parse_dims(c.tok()?)?;
    let a = c.rats(mul_dims(n, n)?)?;
    let b = c.rats(n)?;
    if c.rest_is(&["ERR", "Singular", "singular"]) {
        return Ok(SolveRow {
            id,
            n,
            a,
            b,
            x: None,
        });
    }
    let x = c.rats(n)?;
    Ok(SolveRow {
        id,
        n,
        a,
        b,
        x: Some(x),
    })
}

/// Singular solve rows refuse before any triangular solve (16.LA2.4): certified rows give `Singular "singular"`, others
/// `Singular` or `Resource "undecided_pivot"` (16.LA.4). `BallMat::solve`, `BallLu::solve` and `solve_mat` agree.
#[test]
fn la2_solve_refusal_rows_replay_through_zb() -> Result<(), String> {
    let text = read_gold("lu", "la_solve.txt")?;
    let mut seen = 0_usize;
    let mut cert_seen = 0_usize;
    for prec in PRECS {
        for toks in &data_rows(&text)? {
            let r = solve_row(toks)?;
            if r.x.is_some() {
                continue;
            }
            seen = seen.checked_add(1).ok_or("count")?;
            let cert = SOLVE_CERTIFIED.contains(&r.id);
            if cert {
                cert_seen = cert_seen.checked_add(1).ok_or("count")?;
            }
            let what = format!("{}@{prec} solve", r.id);
            let a = ball_mat(r.n, r.n, &r.a, prec)?;
            let b = ball_vec(&r.b, prec)?;
            want_refusal(&what, a.solve(&b, prec), cert, "undecided_pivot")?;
            let lu = ok(r.id, a.lu(prec))?;
            want_refusal(
                &format!("{what} (BallLu)"),
                lu.solve(&b, prec),
                cert,
                "undecided_pivot",
            )?;
            let bm = ball_mat(r.n, 1, &r.b, prec)?;
            want_refusal(
                &format!("{what} (mat)"),
                a.solve_mat(&bm, prec),
                cert,
                "undecided_pivot",
            )?;
        }
    }
    if seen == 0 || cert_seen == 0 {
        return Err(format!(
            "refusal rows {seen}, certified {cert_seen}: gold ids changed?"
        ));
    }
    Ok(())
}

#[test]
#[ignore = "red: CONTRACT 16.S6.12: BallLu::solve calls solve_lower/solve_upper (spence/zb-la0-trisolve), which are Resource not_implemented stubs on this base; un-ignore when LA0 merges"]
fn la2_solve_gate_value_rows_replay_through_zb() -> Result<(), String> {
    let text = read_gold("lu", "la_solve.txt")?;
    for prec in PRECS {
        for toks in &data_rows(&text)? {
            let r = solve_row(toks)?;
            let Some(x) = &r.x else { continue };
            let what = format!("{}@{prec} solve", r.id);
            let a = ball_mat(r.n, r.n, &r.a, prec)?;
            let b = ball_vec(&r.b, prec)?;
            let got = match a.solve(&b, prec) {
                Ok(v) => v,
                Err(BoundError::Resource { detail, .. }) if detail.starts_with("undecided") => {
                    continue
                }
                Err(e) => return Err(format!("{what}: {e:?}")),
            };
            let bm = ball_mat(r.n, 1, &r.b, prec)?;
            let gm = ok(&what, a.solve_mat(&bm, prec))?;
            for (i, xr) in x.iter().enumerate() {
                want_holds(&format!("{what} x[{i}]"), got.get(i).ok_or("x")?, xr)?;
                want_holds(&format!("{what} X[{i},0]"), entry(&gm, i, 0)?, xr)?;
            }
        }
    }
    Ok(())
}

// ---- hand rows ------------------------------------------------------------------------------------------------------------
fn pt(v: i64) -> Result<Ball, String> {
    ok("exact", Ball::exact(num(v), 64))
}

fn int_mat(n: usize, m: usize, vals: &[i64]) -> Result<BallMat, String> {
    let mut balls = Vec::new();
    for v in vals {
        balls.push(pt(*v)?);
    }
    ok("from_balls", BallMat::from_balls(n, m, balls))
}

fn int_vec(vals: &[i64]) -> Result<BallVec, String> {
    let mut balls = Vec::new();
    for v in vals {
        balls.push(pt(*v)?);
    }
    Ok(BallVec::from_balls(balls))
}

/// `1 ± 2` (straddles 0, nonzero midpoint).
fn straddle() -> Result<Ball, String> {
    ok(
        "straddle",
        Ball::from_mid_rad(num(1), ok("mag", Mag::finite(TOP, 2))?, 64),
    )
}

#[test]
fn la2_rectangular_input_is_dimension_mismatch_everywhere() -> Result<(), String> {
    // no-source: hand-derived: arbitrary 2x3 integer input; only the DimensionMismatch error kind is asserted.
    let a = int_mat(2, 3, &[1, 2, 3, 4, 5, 6])?;
    want_err_detail("lu", a.lu(53), "DimensionMismatch", "rectangular")?;
    want_err_detail("logdet", a.logdet(53), "DimensionMismatch", "rectangular")?;
    want_err_detail(
        "solve",
        a.solve(&int_vec(&[1, 2])?, 53),
        "DimensionMismatch",
        "",
    )?;
    want_err_detail(
        "solve_mat",
        a.solve_mat(&int_mat(2, 1, &[1, 2])?, 53),
        "DimensionMismatch",
        "",
    )
}

/// 16.LA.3: prec is checked first, before shape.
#[test]
fn la2_bad_prec_is_invalid_argument_before_shape() -> Result<(), String> {
    // no-source: hand-derived: arbitrary integer inputs; only the InvalidArgument error kind is asserted (prec 0 is out of range).
    let a = int_mat(2, 3, &[1, 2, 3, 4, 5, 6])?;
    want_err_detail("lu", a.lu(0), "InvalidArgument", "")?;
    want_err_detail("logdet", a.logdet(0), "InvalidArgument", "")?;
    want_err_detail("solve", a.solve(&int_vec(&[1])?, 0), "InvalidArgument", "")?;
    want_err_detail(
        "solve_mat",
        a.solve_mat(&int_mat(1, 1, &[1])?, 0),
        "InvalidArgument",
        "",
    )?;
    let sq = ok("lu", int_mat(1, 1, &[2])?.lu(53))?;
    want_err_detail("BallLu::det", sq.det(0), "InvalidArgument", "")?;
    want_err_detail(
        "BallLu::solve",
        sq.solve(&int_vec(&[1, 2])?, 0),
        "InvalidArgument",
        "",
    )
}

/// 16.LA.3: a right-hand side of the wrong length is `DimensionMismatch "b_len"`.
#[test]
fn la2_solve_b_len_mismatch() -> Result<(), String> {
    // no-source: hand-derived: arbitrary nonsingular 2x2 input; only the b_len error kind is asserted.
    let a = int_mat(2, 2, &[2, 1, 1, 3])?;
    want_err_detail(
        "solve",
        a.solve(&int_vec(&[1, 2, 3])?, 53),
        "DimensionMismatch",
        "b_len",
    )?;
    want_err_detail(
        "solve_mat",
        a.solve_mat(&int_mat(3, 1, &[1, 2, 3])?, 53),
        "DimensionMismatch",
        "b_len",
    )?;
    let lu = ok("lu", a.lu(53))?;
    want_err_detail(
        "BallLu::solve",
        lu.solve(&int_vec(&[1])?, 53),
        "DimensionMismatch",
        "b_len",
    )
}

/// 16.LA2.3: pivot on the largest |mid|; 16.LA2.1 ties go to the lowest row index; perm records the source rows.
#[test]
fn la2_lu_pivot_choice_ties_and_perm() -> Result<(), String> {
    // no-source: hand-derived: perms, signs and det = -1 of these permutation/integer matrices worked out by hand from 16.LA2.1/16.LA2.3.
    // column 0 = (1, -3, 2): pivot row 1 (largest |mid|, negative), perm[0] = 1.
    let lu = ok("lu", int_mat(3, 3, &[1, 0, 0, -3, 1, 0, 2, 0, 1])?.lu(53))?;
    if lu.perm().first() != Some(&1) {
        return Err(format!("largest |mid| pivot: perm {:?}", lu.perm()));
    }
    // column 0 = (2, -2, 2): every |mid| ties, the lowest index (0) wins, no swap.
    let lt = ok("lu", int_mat(3, 3, &[2, 1, 0, -2, 0, 1, 2, 3, 4])?.lu(53))?;
    if lt.perm().first() != Some(&0) {
        return Err(format!("tie must keep row 0: perm {:?}", lt.perm()));
    }
    // cyclic rows (0,1,0) (0,0,1) (1,0,0): P·A = I with perm = [2, 0, 1], an even permutation.
    let lc = ok("lu", int_mat(3, 3, &[0, 1, 0, 0, 0, 1, 1, 0, 0])?.lu(53))?;
    if lc.perm() != [2, 0, 1] || lc.perm_sign() != 1 {
        return Err(format!(
            "cyclic: perm {:?} sign {}",
            lc.perm(),
            lc.perm_sign()
        ));
    }
    // one swap: sign -1 and det = -1 exactly.
    let ls = ok("lu", int_mat(2, 2, &[0, 1, 1, 0])?.lu(53))?;
    let d = ok("det", ls.det(53))?;
    if ls.perm_sign() != -1
        || d.contains(&num(-1)) != Tri::True
        || d.contains(&num(1)) != Tri::False
    {
        return Err(format!("swap2: sign {} det {d:?}", ls.perm_sign()));
    }
    // [[4,3],[2,1]]: no swap, L = [[1,0],[1/2,1]] exact, U = [[4,3],[0,-1/2]] exact.
    // no-source: hand-derived: [[4,3],[2,1]]: l_10 = 2/4 = 1/2, u_11 = 1 - (1/2)·3 = -1/2, all exact dyadics.
    let l2 = ok("lu", int_mat(2, 2, &[4, 3, 2, 1])?.lu(53))?;
    let half = pow2(-1)?;
    let mhalf = -plus_zero(&half)?;
    let checks = [
        (entry(l2.l(), 0, 0)?, num(1)),
        (entry(l2.l(), 0, 1)?, num(0)),
        (entry(l2.l(), 1, 0)?, half),
        (entry(l2.l(), 1, 1)?, num(1)),
        (entry(l2.u(), 0, 0)?, num(4)),
        (entry(l2.u(), 0, 1)?, num(3)),
        (entry(l2.u(), 1, 0)?, num(0)),
        (entry(l2.u(), 1, 1)?, mhalf),
    ];
    for (b, want) in checks {
        if !b.is_exact() || b.contains(&want) != Tri::True {
            return Err(format!("[[4,3],[2,1]]: {b:?} is not exactly {want:?}"));
        }
    }
    Ok(())
}

/// 16.LA2.3/16.LA1.4: an exactly zero column is skipped (no swap, no elimination, `u_kk = 0 ± 0`), `lu` stays `Ok`,
/// `det` is exactly `0 ± 0`, `logdet` and `solve` give the certified `Singular "singular"`.
#[test]
fn la2_zero_column_is_skipped_and_certified_singular() -> Result<(), String> {
    // no-source: hand-derived: column 0 is exactly zero, so the skip rule (16.LA2.3) leaves u_00 = 0 ± 0 and det = 0 exactly.
    let a = int_mat(3, 3, &[0, 1, 2, 0, 3, 4, 0, 5, 7])?;
    let lu = ok("lu", a.lu(53))?;
    if lu.perm().first() != Some(&0) {
        return Err(format!("zero column must not swap: perm {:?}", lu.perm()));
    }
    for i in 1..3 {
        let l = entry(lu.l(), i, 0)?;
        if !is_exact_zero(l) {
            return Err(format!("zero column must not eliminate: L[{i},0] = {l:?}"));
        }
    }
    let u00 = entry(lu.u(), 0, 0)?;
    if !is_exact_zero(u00) {
        return Err(format!("u_00 {u00:?}"));
    }
    let d = ok("det", lu.det(53))?;
    if !is_exact_zero(&d) {
        return Err(format!(
            "det of a skipped column must be exactly 0, got {d:?}"
        ));
    }
    want_err_detail("logdet", a.logdet(53), "Singular", "singular")?;
    want_err_detail(
        "solve",
        a.solve(&int_vec(&[1, 1, 1])?, 53),
        "Singular",
        "singular",
    )?;
    want_err_detail(
        "BallLu::solve",
        lu.solve(&int_vec(&[1, 1, 1])?, 53),
        "Singular",
        "singular",
    )
}

/// A pivot ball straddling 0 (not exactly 0): `lu` is `Ok`, `det` contains 0 but is not exactly 0, `logdet` is
/// `Resource "undecided_det"` (never Domain or a value), `solve` is `Resource "undecided_pivot"` (16.LA2.4).
#[test]
fn la2_straddling_pivot_is_undecided_not_singular() -> Result<(), String> {
    // no-source: hand-derived: 1 ± 2 encloses 0 but its mid is nonzero, so 16.LA2.4 makes it undecided, not singular.
    let a = ok("from_balls", BallMat::from_balls(1, 1, vec![straddle()?]))?;
    let lu = ok("lu", a.lu(53))?;
    let d = ok("det", lu.det(53))?;
    if is_exact_zero(&d) || !has_zero(&d) {
        return Err(format!("straddle det {d:?}"));
    }
    want_err_detail("logdet", a.logdet(53), "Resource", "undecided_det")?;
    want_err_detail(
        "solve",
        a.solve(&int_vec(&[1])?, 53),
        "Resource",
        "undecided_pivot",
    )?;
    want_err_detail(
        "BallLu::solve",
        lu.solve(&int_vec(&[1])?, 53),
        "Resource",
        "undecided_pivot",
    )?;
    // 2x2 with a straddling first pivot: the multiplier 0/(1±2) still encloses 0 and lu is Ok.
    let b = ok(
        "from_balls",
        BallMat::from_balls(2, 2, vec![straddle()?, pt(1)?, pt(0)?, pt(1)?]),
    )?;
    let lb = ok("lu", b.lu(53))?;
    let l10 = entry(lb.l(), 1, 0)?;
    if !has_zero(l10) {
        return Err(format!("l10 {l10:?} must enclose 0"));
    }
    Ok(())
}
/// 16.LA.2 working precision: LU/det run at `prec + 64 + 2·bitlen(n)` and round once, so a det of -1 that cancels
/// 80-bit products comes back tight at prec 53 (kills 1118-WORKPREC: at plain prec 53 the cancellation leaves rad ~2^27).
#[test]
fn la2_det_uses_guarded_work_prec() -> Result<(), String> {
    // no-source: hand-derived: det [[2^40+1, 2^40], [2^40, 2^40-1]] = (2^80 - 1) - 2^80 = -1; 2^-40 is a loose cap (prec-53 rounding of -1 is exact, guard bits keep the rest far below 2^-40).
    let e = |v: i64| ok("exact", Ball::exact(num(v), 53));
    let a = ok(
        "from_balls",
        BallMat::from_balls(
            2,
            2,
            vec![
                e(1_099_511_627_777)?,
                e(1_099_511_627_776)?,
                e(1_099_511_627_776)?,
                e(1_099_511_627_775)?,
            ],
        ),
    )?;
    let d = ok("det", ok("lu", a.lu(53))?.det(53))?;
    if !holds(
        &d,
        &Rat {
            p: num(-1),
            q: num(1),
        },
    )? {
        return Err(format!("det must contain -1, got {d:?}"));
    }
    let rad = d.rad().to_exact_num().map_err(|e| format!("{e:?}"))?;
    if rad > pow2(-40)? {
        return Err(format!("det radius not tight at prec 53: {d:?}"));
    }
    Ok(())
}

/// 16.LA1.4: an overflowing det product gives an all-reals `det` as `Ok`, never an error.
#[test]
fn la2_det_overflow_is_all_reals_ok() -> Result<(), String> {
    // no-source: hand-derived: diag(M, M) with M = TOP·2^EXPONENT_MAX overflows M·M past EXPONENT_MAX, so the det must be all-reals (16.LA1.4).
    let big = ExactNum::from_mag_parts(TOP, EXPONENT_MAX, 64).map_err(|e| format!("{e:?}"))?;
    let b0 = ok("big", Ball::exact(plus_zero(&big)?, 64))?;
    let b1 = ok("big", Ball::exact(big, 64))?;
    let a = ok(
        "from_balls",
        BallMat::from_balls(2, 2, vec![b0, pt(0)?, pt(0)?, b1]),
    )?;
    let lu = ok("lu", a.lu(53))?;
    let d = ok("det", lu.det(53))?;
    if d.is_informative() {
        return Err(format!("overflowing det must be all-reals, got {d:?}"));
    }
    match a.logdet(53) {
        Ok(v) if v.ln_abs_det().is_informative() && v.sign() == 1 => Ok(()),
        Err(BoundError::Resource {
            detail: "needs_ln", ..
        }) => Ok(()),
        other => Err(format!("overflow logdet: {other:?}")),
    }
}

/// 16.LA.6: results are rounded to `prec` once (the balls carry `prec`), and an exact integer system whose factors are
/// dyadic stays exact: det of [[2,1],[4,3]] (pivot 4) is exactly 2 at every test prec, sign from one swap.
#[test]
fn la2_results_carry_prec_and_exact_rows_stay_exact() -> Result<(), String> {
    for prec in PRECS {
        // no-source: hand-derived: P·A = [[4,3],[2,1]] gives U = [[4,3],[0,-1/2]] and det = -(4·(-1/2)) = 2, all dyadic so exact at every prec.
        let lu = ok("lu", int_mat(2, 2, &[2, 1, 4, 3])?.lu(prec))?;
        let d = ok("det", lu.det(prec))?;
        if d.prec() != prec || !d.is_exact() || d.contains(&num(2)) != Tri::True {
            return Err(format!("@{prec}: det {d:?}"));
        }
        for i in 0..2 {
            for j in 0..2 {
                let (l, u) = (entry(lu.l(), i, j)?, entry(lu.u(), i, j)?);
                if l.prec() != prec || u.prec() != prec {
                    return Err(format!(
                        "@{prec}: L/U[{i},{j}] prec {} {}",
                        l.prec(),
                        u.prec()
                    ));
                }
            }
        }
    }
    Ok(())
}

// Held child module (HOLD.md at the repo root): not compiled; the mod line stays here, commented.
// HOLD H09: mod harper_la_fact;
