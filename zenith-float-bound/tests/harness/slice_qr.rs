//! Slice LA3 harness rows (`src/qr.rs`): gold replay of `la_qr.txt` (`BallMat::qr`) and `la_lstsq.txt` (`BallMat::lstsq`).
//! Owner: Rusty (slice LA3, Cassie review 2026-10-29). Exact qr golds: Cassie; lstsq golds: scripts/gen_la_lstsq_goldens.py.
//! Goldens: tests/goldens/slices/qr/ (its own SHA256SUMS and GENERATORS). Mutants: tools/mutants/manifest.d/qr.txt, ids 1200..=1299.
//! COVERS: BallMat::qr, BallMat::lstsq, BallQr::q, BallQr::r, BallQr::into_parts, BallQr::lstsq
//! Golden file headers for this slice are declared here, not in tests/common/goldens.rs.
//!
//! Replay (CONTRACT 16.LA.4, 16.LA.5, 16.LA3.3): every row at each of `PRECS`. A value passes when every returned ball contains
//! the exact value; here, when it covers the gold `[lo256, hi256]` (sufficient, 16.LA.5). A ball that only overlaps it would be
//! re-decided against the heavy enc4096 block, which the gate tier does not carry, so it is reported as a failure naming the entry.
//! lstsq value rows carry the format-2 cond columns (16.LA.5, scripts/la_cond_cols.py la_cond_cols_lstsq); the replay is
//! containment-only for now (the `k_row` tightness bound is a follow-up commit), the columns themselves are checked below.

use crate::common::goldens::golden_dir;
use zenith_float_bound::array::{BallMat, BallVec};
use zenith_float_bound::ball::Ball;
use zenith_float_bound::error::BoundError;
use zenith_float_num::{ExactNum, RoundingMode};

/// Replay precisions (the files have no `prec` column), as the LA1 logdet replay (16.LA1.7).
const PRECS: [usize; 3] = [53, 113, 192];
/// 16.LA.4: qr rows certified with exact arithmetic; they must return the named error, never `Resource "undecided_*"`.
const QR_EXACT_CERTIFIED: [&str; 2] = ["qr.refuse.zero2", "qr.refuse.onebyone_zero"];
/// lstsq rows whose refusal is certified exactly (an exact `0 ± 0` column, 16.LA3.3). Not yet in the 16.LA.4 list (asked of Spence).
const LSTSQ_EXACT_CERTIFIED: [&str; 4] = [
    "lstsq.refuse.zero_col_3x2",
    "lstsq.refuse.blen_3x2_b2",
    "lstsq.refuse.blen_3x2_b4",
    "lstsq.refuse.blen_3x0_b2",
];
/// Trailing format-2 cond columns of an lstsq value row with `n >= 1` (16.LA.5): `cond_row k_row scale_row containment_only reason
/// rref`. `m x 0` rows (empty x) and refusal rows have none.
const LSTSQ_COND_COLS: usize = 6;
/// 16.LA.5 lstsq `reason` tokens, keyed by `containment_only`.
const LSTSQ_REASONS: [(&str, &str); 3] = [
    ("0", "tight"),
    ("1", "k_row_too_large"),
    ("2", "lstsq_x_zero"),
];

/// Value rows whose computed `σ` ball reaches 0 at a replay precision, pinned to the measured `Resource "undecided_rank"`
/// (16.LA3.3, measured 2026-10-29; matched on the id after the `qr.`/`lstsq.` prefix). A pin that starts returning a value fails
/// as stale. `nearsing.2`: the input balls at `prec <= 192` already contain singular matrices (perturbation `2^-200`).
/// `hilbert.09..12` @53: naive ball Householder wraps the trailing radius by about 3 bits per step (`2^-54` to `2^-30` over
/// 8 steps for n = 9) while `σ_n` falls to about `2^-35`. `scale_huge.undecided_rank`: `R22/R11 ≈ 2^-10000`, beyond every working precision.
/// `straddle_pivot.3x2` @53 (mutant 1200): the step-2 pivot ball straddles 0 and so does `σ`.
const UNDECIDED_PINS: [(&str, usize); 11] = [
    ("nearsing.2", 53),
    ("nearsing.2", 113),
    ("nearsing.2", 192),
    ("hilbert.09", 53),
    ("hilbert.10", 53),
    ("hilbert.11", 53),
    ("hilbert.12", 53),
    ("scale_huge.undecided_rank", 53),
    ("scale_huge.undecided_rank", 113),
    ("scale_huge.undecided_rank", 192),
    ("straddle_pivot.3x2", 53),
];

/// The Hilbert rows measured at 256 bits (task d): all factor; each ball must at least overlap the gold `[lo256, hi256]`.
const HILBERT_256_ROWS: usize = 8;

fn is_pinned(id: &str, prec: usize) -> bool {
    let key = id.split_once('.').map_or(id, |(_, s)| s);
    UNDECIDED_PINS.iter().any(|&(k, p)| k == key && p == prec)
}

/// The pinned outcome: `Resource "undecided_rank"`; anything else fails.
fn check_pin(id: &str, prec: usize, e: Option<&BoundError>) -> Result<(), String> {
    match e {
        Some(BoundError::Resource { detail, .. }) if *detail == "undecided_rank" => Ok(()),
        Some(e) => Err(format!("{id} @{prec}: pinned undecided_rank, got {e:?}")),
        None => Err(format!(
            "{id} @{prec}: pinned undecided_rank, got a value (stale pin)"
        )),
    }
}

fn gold_text(name: &str) -> Result<String, String> {
    let path = golden_dir().join("slices").join("qr").join(name);
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

fn data_rows(text: &str) -> impl Iterator<Item = Vec<&str>> {
    text.lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split_whitespace().collect())
}

fn dims(t: &str) -> Result<(usize, usize), String> {
    let (m, n) = t.split_once('x').ok_or_else(|| format!("dims {t:?}"))?;
    let p = |s: &str| s.parse::<usize>().map_err(|e| format!("dims {t:?}: {e}"));
    Ok((p(m)?, p(n)?))
}

fn mul(a: usize, b: usize) -> Result<usize, String> {
    a.checked_mul(b).ok_or_else(|| "size overflow".to_string())
}

fn add(a: usize, b: usize) -> Result<usize, String> {
    a.checked_add(b).ok_or_else(|| "size overflow".to_string())
}

/// A decimal integer, exactly (Horner in 18-digit chunks with full-precision ops).
fn int(t: &str) -> Result<ExactNum, String> {
    let (neg, digits) = match t.strip_prefix('-') {
        Some(d) => (true, d),
        None => (false, t),
    };
    if digits.is_empty() || !digits.bytes().all(|c| c.is_ascii_digit()) {
        return Err(format!("bad integer {t:?}"));
    }
    let mut acc = ExactNum::from(0i64);
    let mut rest = digits;
    while !rest.is_empty() {
        let k = rest.len().min(18);
        let (head, tail) = rest.split_at(k);
        let chunk = head.parse::<i64>().map_err(|e| format!("{t:?}: {e}"))?;
        let k32 = u32::try_from(k).map_err(|e| format!("{e}"))?;
        let scale = 10i64.checked_pow(k32).ok_or("chunk scale")?;
        acc = acc
            .mul_full_prec(&ExactNum::from(scale))
            .add_full_prec(&ExactNum::from(chunk));
        rest = tail;
    }
    let v = if neg {
        ExactNum::from(0i64).sub_full_prec(&acc)
    } else {
        acc
    };
    if v.is_nan() {
        return Err(format!("integer {t:?} did not parse exactly"));
    }
    Ok(v)
}

/// A gold rational `p/q` as a ball at `prec` (exact when `q = 1`, else the S1 `div` enclosure).
fn rational(t: &str, prec: usize) -> Result<Ball, String> {
    let (p, q) = t
        .split_once('/')
        .ok_or_else(|| format!("bad rational {t:?}"))?;
    let num = Ball::exact(int(p)?, prec).map_err(|e| format!("{t}: {e:?}"))?;
    if q == "1" {
        return Ok(num);
    }
    let den = Ball::exact(int(q)?, prec).map_err(|e| format!("{t}: {e:?}"))?;
    num.div(&den, prec).map_err(|e| format!("{t}: {e:?}"))
}

/// A gold dyadic `m*2^e`, exactly.
fn dyadic(t: &str) -> Result<ExactNum, String> {
    let (m, e) = t
        .split_once("*2^")
        .ok_or_else(|| format!("bad dyadic {t:?}"))?;
    let e = e.parse::<i64>().map_err(|err| format!("{t:?}: {err}"))?;
    let bits = m
        .len()
        .checked_mul(4)
        .and_then(|b| b.checked_add(64))
        .ok_or("dyadic width")?;
    let v = int(m)?.ldexp(e, bits, RoundingMode::None);
    if v.is_nan() || v.inexact() {
        return Err(format!("dyadic {t:?} not exact"));
    }
    Ok(v)
}

/// `true` when `b` covers `[lo, hi]` exactly: `mid - rad <= lo` and `hi <= mid + rad`.
fn covers(b: &Ball, lo: &ExactNum, hi: &ExactNum) -> Result<bool, String> {
    if b.rad().is_infinity() {
        return Ok(true);
    }
    let r = b.rad().to_exact_num().map_err(|e| format!("{e:?}"))?;
    let low_gap = lo.sub_full_prec(b.mid()).add_full_prec(&r);
    let high_gap = b.mid().add_full_prec(&r).sub_full_prec(hi);
    if low_gap.is_nan() || high_gap.is_nan() {
        return Err("containment arithmetic failed".to_string());
    }
    Ok(!low_gap.is_negative() && !high_gap.is_negative())
}

/// `true` when `b` meets `[lo, hi]`: `mid - rad <= hi` and `lo <= mid + rad`.
fn overlaps(b: &Ball, lo: &ExactNum, hi: &ExactNum) -> Result<bool, String> {
    if b.rad().is_infinity() {
        return Ok(true);
    }
    let r = b.rad().to_exact_num().map_err(|e| format!("{e:?}"))?;
    let below = hi.sub_full_prec(b.mid()).add_full_prec(&r);
    let above = b.mid().add_full_prec(&r).sub_full_prec(lo);
    if below.is_nan() || above.is_nan() {
        return Err("overlap arithmetic failed".to_string());
    }
    Ok(!below.is_negative() && !above.is_negative())
}

/// 256-bit check: the gate tier cannot re-decide a partial overlap (no enc4096), so only disjointness fails.
fn check_overlap(what: &str, b: Option<&Ball>, pair: &[&str]) -> Result<(), String> {
    let b = b.ok_or_else(|| format!("{what}: missing entry"))?;
    let (lo, hi) = match pair {
        [lo, hi] => (dyadic(lo)?, dyadic(hi)?),
        _ => return Err(format!("{what}: bad lo/hi pair")),
    };
    if overlaps(b, &lo, &hi)? {
        Ok(())
    } else {
        Err(format!("{what}: ball is disjoint from gold [lo256, hi256]"))
    }
}

type Check = fn(&str, Option<&Ball>, &[&str]) -> Result<(), String>;

fn check_cover(what: &str, b: Option<&Ball>, pair: &[&str]) -> Result<(), String> {
    let b = b.ok_or_else(|| format!("{what}: missing entry"))?;
    let (lo, hi) = match pair {
        [lo, hi] => (dyadic(lo)?, dyadic(hi)?),
        _ => return Err(format!("{what}: bad lo/hi pair")),
    };
    if covers(b, &lo, &hi)? {
        Ok(())
    } else {
        Err(format!(
            "{what}: ball {:?} ± {:?} does not cover gold [lo256, hi256] (enc4096 re-decision needs the heavy tier)",
            b.mid(),
            b.rad()
        ))
    }
}

fn variant(e: &BoundError) -> String {
    format!("{e:?}")
        .split(|c: char| !c.is_alphanumeric())
        .next()
        .map_or_else(String::new, str::to_string)
}

/// 16.LA.4: a refusal row passes on the named error, or on `Resource "undecided_*"` unless the row is exact-certified.
fn check_refusal(
    id: &str,
    got: &BoundError,
    want: &[&str],
    certified: &[&str],
) -> Result<(), String> {
    let (v, d) = match want {
        [v, d] => (*v, *d),
        _ => return Err(format!("{id}: bad ERR tail")),
    };
    if variant(got) == v && got.detail() == d {
        return Ok(());
    }
    let shape = v == "DimensionMismatch";
    if !shape
        && !certified.contains(&id)
        && variant(got) == "Resource"
        && got.detail().starts_with("undecided_")
    {
        return Ok(());
    }
    Err(format!("{id}: want ERR {v} {d}, got {got:?}"))
}

fn input_mat<'a>(
    it: &mut impl Iterator<Item = &'a &'a str>,
    m: usize,
    n: usize,
    prec: usize,
) -> Result<BallMat, String> {
    let mut v = Vec::new();
    for _ in 0..mul(m, n)? {
        v.push(rational(it.next().ok_or("short row (A)")?, prec)?);
    }
    BallMat::from_balls(m, n, v).map_err(|e| format!("{e:?}"))
}

/// One qr row at `prec`: `case_id mxn A[m*n] (Qt[m*n] Rt[n*n] d[n] Q_lohi[2mn] R_lohi[2nn] | ERR v d)`.
fn replay_qr(toks: &[&str], prec: usize) -> Result<(), String> {
    replay_qr_with(toks, prec, check_cover)
}

fn replay_qr_with(toks: &[&str], prec: usize, check: Check) -> Result<(), String> {
    let id = toks.first().ok_or("empty row")?;
    let (m, n) = dims(toks.get(1).ok_or("no dims")?)?;
    let mut it = toks.iter().skip(2);
    let a = input_mat(&mut it, m, n, prec)?;
    let rest = toks.get(add(2, mul(m, n)?)?..).ok_or("short row")?;
    let got = a.qr(prec);
    if is_pinned(id, prec) {
        return check_pin(id, prec, got.as_ref().err());
    }
    if rest.first() == Some(&"ERR") {
        return match got {
            Ok(_) => Err(format!("{id} @{prec}: value where the gold refuses")),
            Err(e) => check_refusal(
                id,
                &e,
                rest.get(1..).ok_or("short ERR row")?,
                &QR_EXACT_CERTIFIED,
            ),
        };
    }
    let f = got.map_err(|e| format!("{id} @{prec}: {e:?}"))?;
    let (q, r) = (f.q(), f.r());
    if (q.nrows(), q.ncols(), r.nrows(), r.ncols()) != (m, n, n, n) {
        return Err(format!(
            "{id} @{prec}: Q {}x{}, R {}x{} (want thin {m}x{n}, {n}x{n})",
            q.nrows(),
            q.ncols(),
            r.nrows(),
            r.ncols()
        ));
    }
    let skip = add(add(mul(m, n)?, mul(n, n)?)?, n)?;
    let mut pairs = rest
        .get(skip..)
        .ok_or("short row (enclosures)")?
        .as_chunks::<2>()
        .0
        .iter();
    for i in 0..m {
        for j in 0..n {
            check(
                &format!("{id} @{prec} Q[{i}][{j}]"),
                q.get(i, j),
                &pairs.next().ok_or("short Q")?[..],
            )?;
        }
    }
    for i in 0..n {
        for j in 0..n {
            check(
                &format!("{id} @{prec} R[{i}][{j}]"),
                r.get(i, j),
                &pairs.next().ok_or("short R")?[..],
            )?;
        }
    }
    Ok(())
}

/// One lstsq row at `prec`: `case_id mxn blen A[m*n] b[blen] (x[n] x_lohi[2n] COND_COLS[6] | ERR v d)`.
fn replay_lstsq(toks: &[&str], prec: usize) -> Result<(), String> {
    let id = toks.first().ok_or("empty row")?;
    let (m, n) = dims(toks.get(1).ok_or("no dims")?)?;
    let blen = toks
        .get(2)
        .ok_or("no blen")?
        .parse::<usize>()
        .map_err(|e| format!("{id}: {e}"))?;
    let mut it = toks.iter().skip(3);
    let a = input_mat(&mut it, m, n, prec)?;
    let mut b = Vec::new();
    for _ in 0..blen {
        b.push(rational(it.next().ok_or("short row (b)")?, prec)?);
    }
    let rest = toks
        .get(add(add(3, mul(m, n)?)?, blen)?..)
        .ok_or("short row")?;
    let got = a.lstsq(&BallVec::from_balls(b), prec);
    if is_pinned(id, prec) {
        return check_pin(id, prec, got.as_ref().err());
    }
    if rest.first() == Some(&"ERR") {
        return match got {
            Ok(_) => Err(format!("{id} @{prec}: value where the gold refuses")),
            Err(e) => check_refusal(
                id,
                &e,
                rest.get(1..).ok_or("short ERR row")?,
                &LSTSQ_EXACT_CERTIFIED,
            ),
        };
    }
    let x = got.map_err(|e| format!("{id} @{prec}: {e:?}"))?;
    if x.len() != n {
        return Err(format!("{id} @{prec}: x has {} entries, want {n}", x.len()));
    }
    let mut pairs = rest
        .get(n..)
        .ok_or("short row (x enclosures)")?
        .as_chunks::<2>()
        .0
        .iter();
    for i in 0..n {
        check_cover(
            &format!("{id} @{prec} x[{i}]"),
            x.get(i),
            &pairs.next().ok_or("short x")?[..],
        )?;
    }
    Ok(())
}

/// Replays every row at every `PRECS`, collecting failures so one run names them all.
fn replay_all(text: &str, f: Replay) -> Result<usize, String> {
    replay_some(text, f, |_, _| true)
}

type Replay = fn(&[&str], usize) -> Result<(), String>;

/// As `replay_all`, for the `(row, prec)` pairs `keep` selects.
fn replay_some(text: &str, f: Replay, keep: fn(&[&str], usize) -> bool) -> Result<usize, String> {
    let mut fails = Vec::new();
    let mut n = 0_usize;
    for toks in data_rows(text) {
        for prec in PRECS {
            if !keep(&toks, prec) {
                continue;
            }
            n = add(n, 1)?;
            if let Err(e) = f(&toks, prec) {
                fails.push(e);
            }
        }
    }
    if fails.is_empty() {
        Ok(n)
    } else {
        Err(format!(
            "{} of {n} replays failed:\n{}",
            fails.len(),
            fails.join("\n")
        ))
    }
}

fn qr_width(m: usize, n: usize, refusal: bool) -> Result<usize, String> {
    let tail = if refusal {
        3
    } else {
        add(
            add(add(mul(m, n)?, mul(n, n)?)?, n)?,
            mul(2, add(mul(m, n)?, mul(n, n)?)?)?,
        )?
    };
    add(add(2, mul(m, n)?)?, tail)
}

fn lstsq_width(m: usize, n: usize, blen: usize, refusal: bool) -> Result<usize, String> {
    let tail = if refusal {
        3
    } else if n == 0 {
        0
    } else {
        add(mul(3, n)?, LSTSQ_COND_COLS)?
    };
    add(add(add(3, mul(m, n)?)?, blen)?, tail)
}

fn check_widths(text: &str, lstsq: bool) -> Result<usize, String> {
    let mut seen = 0_usize;
    for toks in data_rows(text) {
        let id = toks.first().ok_or("empty row")?;
        let (m, n) = dims(toks.get(1).ok_or("no dims")?)?;
        let refusal = toks.contains(&"ERR");
        let w = if lstsq {
            let blen = toks
                .get(2)
                .ok_or("no blen")?
                .parse::<usize>()
                .map_err(|e| format!("{id}: {e}"))?;
            lstsq_width(m, n, blen, refusal)?
        } else {
            qr_width(m, n, refusal)?
        };
        if toks.len() != w {
            return Err(format!("{id}: {} fields, want {w}", toks.len()));
        }
        seen = add(seen, 1)?;
    }
    let want = header_rows(text)?;
    if seen != want {
        return Err(format!("{seen} data rows, header says {want}"));
    }
    Ok(seen)
}

#[test]
fn la3_qr_gate_rows_match_header_count_and_field_widths() -> Result<(), String> {
    check_widths(&gold_text("la_qr.txt")?, false).map(|_| ())
}

#[test]
fn la3_lstsq_gate_rows_match_header_count_and_field_widths() -> Result<(), String> {
    check_widths(&gold_text("la_lstsq.txt")?, true).map(|_| ())
}

/// The cond columns of every lstsq value row with `n >= 1` (16.LA.5).
fn lstsq_cond_tail<'a>(toks: &'a [&'a str]) -> Result<&'a [&'a str], String> {
    toks.len()
        .checked_sub(LSTSQ_COND_COLS)
        .and_then(|k| toks.get(k..))
        .ok_or_else(|| format!("{:?}: short row (cond columns)", toks.first()))
}

#[test]
fn la3_lstsq_cond_columns_are_format_2_16_la_5() -> Result<(), String> {
    let text = gold_text("la_lstsq.txt")?;
    let mut seen = 0_usize;
    for toks in data_rows(&text).filter(|t| !t.contains(&"ERR")) {
        if dims(toks.get(1).ok_or("no dims")?)?.1 == 0 {
            continue;
        }
        let id = toks.first().ok_or("empty row")?;
        let t = lstsq_cond_tail(&toks)?;
        let [cond, k_row, scale, co, reason, rref] = t else {
            return Err(format!("{id}: cond columns {t:?}"));
        };
        for (name, v) in [("cond_row", cond), ("scale_row", scale), ("rref", rref)] {
            let (p, q) = v
                .split_once('/')
                .ok_or_else(|| format!("{id}: {name} {v} not p/q"))?;
            let ok = !p.is_empty()
                && p.bytes().all(|c| c.is_ascii_digit())
                && !q.is_empty()
                && q.bytes().all(|c| c.is_ascii_digit())
                && !q.starts_with('0');
            if !ok {
                return Err(format!(
                    "{id}: {name} {v} is not a nonnegative p/q with q > 0"
                ));
            }
        }
        let k = k_row
            .parse::<i64>()
            .map_err(|e| format!("{id}: k_row {e}"))?;
        if !LSTSQ_REASONS.contains(&(*co, *reason)) {
            return Err(format!(
                "{id}: containment_only {co} reason {reason} (16.LA.5)"
            ));
        }
        let max_prec = PRECS.iter().max().ok_or("no PRECS")?;
        let big = i64::try_from(*max_prec).map_err(|e| e.to_string())?;
        let too_large = k.checked_add(8).ok_or("k_row overflow")? >= big;
        if *reason != "lstsq_x_zero" && too_large != (*co == "1") {
            return Err(format!(
                "{id}: k_row {k} vs containment_only {co} (k_row + 8 >= {big})"
            ));
        }
        if (*reason == "lstsq_x_zero") != (*scale == "0/1") {
            return Err(format!(
                "{id}: reason {reason} with scale_row {scale} (x = 0 iff lstsq_x_zero)"
            ));
        }
        seen = add(seen, 1)?;
    }
    if seen == 0 {
        return Err("no lstsq value rows".to_string());
    }
    Ok(())
}

/// 16.LA.5: an exact-fit row (r = 0) reduces cond_row to the first term `|| |A^+||A| ||_inf`. `exact_fit.3x2`:
/// A = [[1,0],[0,1],[1,1]], b = [1,2,3], x = [1,2], |A^+| = [[2,1,1],[1,2,1]]/3, |A^+||A| = [[1,2/3],[2/3,1]]: cond 5/3,
/// k_row = 4 + ceil(log2 3) + ceil(log2 5/3) = 7, scale 2, rref = 3 * max(1 + 4/3, 2/3 + 2) = 8. `x_zero.3x2`: x = 0, r = b != 0 (A^T r = 0, so rref is 0).
#[test]
fn la3_lstsq_exact_fit_and_x_zero_cond_rows_16_la_5() -> Result<(), String> {
    let text = gold_text("la_lstsq.txt")?;
    for (id, want) in [
        (
            "lstsq.exact_fit.3x2",
            ["5/3", "7", "2/1", "0", "tight", "8/1"],
        ),
        (
            "lstsq.x_zero.3x2",
            ["1/1", "6", "0/1", "2", "lstsq_x_zero", "0/1"],
        ),
    ] {
        let toks = data_rows(&text)
            .find(|t| t.first() == Some(&id))
            .ok_or_else(|| format!("{id}: row missing"))?;
        let t = lstsq_cond_tail(&toks)?;
        if t != want {
            return Err(format!("{id}: cond columns {t:?}, want {want:?}"));
        }
    }
    Ok(())
}

#[test]
fn la3_exact_certified_rows_exist_and_refuse() -> Result<(), String> {
    for (file, list) in [
        ("la_qr.txt", &QR_EXACT_CERTIFIED[..]),
        ("la_lstsq.txt", &LSTSQ_EXACT_CERTIFIED[..]),
    ] {
        let text = gold_text(file)?;
        for id in list {
            let row = data_rows(&text)
                .find(|t| t.first() == Some(id))
                .ok_or_else(|| format!("{file}: no row {id}"))?;
            if !row.contains(&"ERR") {
                return Err(format!("{id}: exact-certified row is not a refusal"));
            }
        }
    }
    Ok(())
}

#[test]
fn la3_harness_parses_and_covers() -> Result<(), String> {
    let v = dyadic("-3*2^-2")?;
    let b = Ball::exact(int("-3")?.ldexp(-2, 64, RoundingMode::None), 53)
        .map_err(|e| format!("{e:?}"))?;
    if !covers(&b, &v, &v)? || covers(&b, &dyadic("-1*2^0")?, &v)? {
        return Err("covers wrong on an exact point".to_string());
    }
    let third = rational("1/3", 53)?;
    let lo = dyadic("6148914691236517205*2^-64")?;
    let hi = dyadic("6148914691236517206*2^-64")?;
    if !covers(&third, &lo, &hi)? {
        return Err("1/3 ball does not cover its 64-bit bracket".to_string());
    }
    let big = int("123456789012345678901234567890")?;
    let want = ExactNum::from(123456789012i64)
        .mul_full_prec(&ExactNum::from(1_000_000_000_000_000_000i64))
        .add_full_prec(&ExactNum::from(345678901234567890i64));
    if !big.sub_full_prec(&want).is_zero() {
        return Err("decimal parse wrong".to_string());
    }
    Ok(())
}

/// Every qr row at every `PRECS` (green with the LA3 bodies; pins in `UNDECIDED_PINS`).
#[test]
fn la3_qr_gate_rows_replay_through_zb() -> Result<(), String> {
    replay_all(&gold_text("la_qr.txt")?, replay_qr).map(|_| ())
}

/// lstsq rows that stop before 16.S6.12 `solve_upper`: refusals, `m×0` systems, and the pinned undecided rows (green now).
#[test]
fn la3_lstsq_refusal_empty_and_pinned_rows_replay_through_zb() -> Result<(), String> {
    replay_some(
        &gold_text("la_lstsq.txt")?,
        replay_lstsq,
        stops_before_solve,
    )
    .map(|_| ())
}

fn stops_before_solve(toks: &[&str], prec: usize) -> bool {
    toks.contains(&"ERR")
        || toks.get(1).is_some_and(|d| d.ends_with("x0"))
        || toks.first().is_some_and(|id| is_pinned(id, prec))
}

/// Every lstsq row (16.S6.12 `solve_upper`, LA0, is implemented).
#[test]
fn la3_lstsq_gate_rows_replay_through_zb() -> Result<(), String> {
    replay_all(&gold_text("la_lstsq.txt")?, replay_lstsq).map(|_| ())
}

/// Task (d): Hilbert n = 5..=12 through `qr` at 256 bits. Every row factors (no `Resource`), and every ball overlaps the gold.
#[test]
fn la3_qr_hilbert_rows_factor_at_256_bits() -> Result<(), String> {
    let text = gold_text("la_qr.txt")?;
    let mut n = 0_usize;
    for toks in data_rows(&text) {
        if toks.first().is_some_and(|id| id.starts_with("qr.hilbert.")) {
            replay_qr_with(&toks, 256, check_overlap)?;
            n = add(n, 1)?;
        }
    }
    if n == HILBERT_256_ROWS {
        Ok(())
    } else {
        Err(format!("{n} hilbert rows, want {HILBERT_256_ROWS}"))
    }
}

/// Item 8 (Spence, 2026-10-29, P0): Hilbert n = 5..=12 at 256 bits must CONTAIN the exact Q, R. Each ball that only overlaps
/// `[lo256, hi256]` is re-decided against the heavy enc4096 block, then enc16384; disjoint at any tier fails. `(row, tier)` is the
/// highest tier a row's entries needed (256, 4096, 16384), measured 2026-10-29; a drift fails, naming the row.
const HILBERT_TIERS: [(&str, usize); 8] = [
    ("hilbert.05", 4096),
    ("hilbert.06", 4096),
    ("hilbert.07", 4096),
    ("hilbert.08", 4096),
    ("hilbert.09", 4096),
    ("hilbert.10", 4096),
    ("hilbert.11", 4096),
    ("hilbert.12", 4096),
];
/// Rows still only overlapping at enc16384: `containment_only` with a recorded reason (none: every row resolved at enc4096, 2026-10-29).
const HILBERT_CONTAINMENT_ONLY: [(&str, &str); 0] = [];

fn heavy_rows(dir: &str, name: &str) -> Result<Option<String>, String> {
    let path = std::path::Path::new(dir).join(name);
    if path.exists() {
        std::fs::read_to_string(&path)
            .map(Some)
            .map_err(|e| format!("cannot read {}: {e}", path.display()))
    } else {
        Ok(None)
    }
}

/// The lo/hi pairs of row `id` in a heavy enc file: `case_id dims A[m*n] Q_lohi[2mn] R_lohi[2nn]`.
fn enc_pairs<'a>(text: &'a str, id: &str, mn: usize) -> Result<Vec<&'a str>, String> {
    let mut toks = data_rows(text)
        .find(|t| t.first() == Some(&id))
        .ok_or_else(|| format!("{id}: no row in the heavy enc file"))?;
    let start = add(2, mn)?;
    if toks.len() < start {
        return Err(format!("{id}: short heavy row"));
    }
    // Drop the leading case_id, dims and A entries in place: the row's own Vec is returned, nothing is copied.
    toks.drain(..start);
    Ok(toks)
}

fn pair_at(p: &[&str], k: usize) -> Result<(ExactNum, ExactNum), String> {
    let i = mul(k, 2)?;
    match (p.get(i), p.get(add(i, 1)?)) {
        (Some(lo), Some(hi)) => Ok((dyadic(lo)?, dyadic(hi)?)),
        _ => Err(format!("entry {k}: missing lo/hi pair")),
    }
}

#[test]
#[ignore = "heavy tier (16.LA.5): set ZB_QR_HEAVY_DIR to the dir with la_qr_enc4096_heavy.txt (scripts/gen_la_goldens.py --heavy DIR; sha256 in HEAVY.SHA256) and optionally la_qr_enc16384_heavy.txt"]
fn la3_qr_hilbert_256_contains_by_heavy_tier() -> Result<(), String> {
    // As s6_heavy_rows_replay_through_zb (ZB-P03): with ZB_QR_HEAVY_DIR unset this ignored test skips with a visible
    // SKIP line, so a `--ignored` run without the heavy files is not red; it runs nothing and says so.
    let dir = match std::env::var("ZB_QR_HEAVY_DIR") {
        Ok(d) => d,
        Err(std::env::VarError::NotPresent) => {
            println!("SKIP la3_qr_hilbert_256_contains_by_heavy_tier: ZB_QR_HEAVY_DIR unset, no heavy rows run");
            return Ok(());
        }
        Err(std::env::VarError::NotUnicode(v)) => {
            return Err(format!("ZB_QR_HEAVY_DIR is not unicode: {v:?}"))
        }
    };
    let e4 = heavy_rows(&dir, "la_qr_enc4096_heavy.txt")?.ok_or("no la_qr_enc4096_heavy.txt")?;
    let e16 = heavy_rows(&dir, "la_qr_enc16384_heavy.txt")?;
    let text = gold_text("la_qr.txt")?;
    let mut got: Vec<(String, usize)> = Vec::new();
    for toks in data_rows(&text) {
        let Some(id) = toks.first().filter(|id| id.starts_with("qr.hilbert.")) else {
            continue;
        };
        let (m, n) = dims(toks.get(1).ok_or("no dims")?)?;
        let mn = mul(m, n)?;
        let mut it = toks.iter().skip(2);
        let f = input_mat(&mut it, m, n, 256)?
            .qr(256)
            .map_err(|e| format!("{id} @256: {e:?}"))?;
        let skip = add(add(add(2, mn)?, add(mn, mul(n, n)?)?)?, n)?;
        let g256 = toks.get(skip..).ok_or("short row (enclosures)")?;
        let p4 = enc_pairs(&e4, id, mn)?;
        let p16 = match &e16 {
            Some(t) => Some(enc_pairs(t, id, mn)?),
            None => None,
        };
        let short = id.strip_prefix("qr.").unwrap_or(id);
        let mut tier = 256_usize;
        let total = add(mn, mul(n, n)?)?;
        for k in 0..total {
            let b = if k < mn {
                f.q().get(k / n, k % n)
            } else {
                let r = k.checked_sub(mn).ok_or("index")?;
                f.r().get(r / n, r % n)
            }
            .ok_or_else(|| format!("{id}: missing entry {k}"))?;
            let (lo, hi) = pair_at(g256, k)?;
            if covers(b, &lo, &hi)? {
                continue;
            }
            if !overlaps(b, &lo, &hi)? {
                return Err(format!("{id} entry {k}: disjoint from [lo256, hi256]"));
            }
            let (lo, hi) = pair_at(&p4, k)?;
            if covers(b, &lo, &hi)? {
                tier = tier.max(4096);
                continue;
            }
            if !overlaps(b, &lo, &hi)? {
                return Err(format!("{id} entry {k}: disjoint from [lo4096, hi4096]"));
            }
            let Some(p16) = &p16 else {
                return Err(format!(
                    "{id} entry {k}: overlaps enc4096; needs la_qr_enc16384_heavy.txt"
                ));
            };
            let (lo, hi) = pair_at(p16, k)?;
            if covers(b, &lo, &hi)? {
                tier = tier.max(16384);
            } else if !overlaps(b, &lo, &hi)? {
                return Err(format!("{id} entry {k}: disjoint from [lo16384, hi16384]"));
            } else if HILBERT_CONTAINMENT_ONLY.iter().any(|(c, _)| *c == short) {
                tier = tier.max(16385);
            } else {
                return Err(format!(
                    "{id} entry {k}: still only overlaps at enc16384 and is not containment_only"
                ));
            }
        }
        got.push((short.to_string(), tier));
    }
    let want: Vec<(String, usize)> = HILBERT_TIERS
        .iter()
        .map(|(i, t)| ((*i).to_string(), *t))
        .collect();
    if got == want {
        Ok(())
    } else {
        Err(format!("hilbert tiers measured {got:?}, recorded {want:?}"))
    }
}

/// `BallQr::into_parts` on A = [3; 4]: R is 1x1 with |R_00| = sqrt(3^2 + 4^2) = 5, Q is 2x1; contained at every replay prec.
#[test]
fn la3_qr_into_parts_3_4_5() -> Result<(), String> {
    let e = |n: i64| ExactNum::from(n);
    for prec in PRECS {
        let ball = |v: i64| Ball::exact(ExactNum::from(v), prec).map_err(|er| format!("{er:?}"));
        let a =
            BallMat::from_balls(2, 1, vec![ball(3)?, ball(4)?]).map_err(|er| format!("{er:?}"))?;
        let (q, r) = a.qr(prec).map_err(|er| format!("qr: {er:?}"))?.into_parts();
        let r00 = r.get(0, 0).ok_or("R_00 missing")?;
        if !(covers(r00, &e(5), &e(5))? || covers(r00, &e(-5), &e(-5))?) {
            return Err(format!("prec {prec}: into_parts R_00 does not contain ±5"));
        }
        if q.get(1, 0).is_none()
            || q.get(2, 0).is_some()
            || q.get(0, 1).is_some()
            || r.get(0, 1).is_some()
        {
            return Err(format!(
                "prec {prec}: into_parts shapes are not Q 2x1, R 1x1"
            ));
        }
    }
    Ok(())
}

/// `BallQr::lstsq` and `BallMat::lstsq` on A = [3; 4], b = [3; 4]: the exact fit x = 1 (A·1 = b), contained at every replay prec.
#[test]
fn la3_factored_and_direct_lstsq_on_3_4_5() -> Result<(), String> {
    let e = |n: i64| ExactNum::from(n);
    for prec in PRECS {
        let ball = |v: i64| Ball::exact(ExactNum::from(v), prec).map_err(|er| format!("{er:?}"));
        let a =
            BallMat::from_balls(2, 1, vec![ball(3)?, ball(4)?]).map_err(|er| format!("{er:?}"))?;
        let x_f = a
            .qr(prec)
            .map_err(|er| format!("qr: {er:?}"))?
            .lstsq(&BallVec::from_balls(vec![ball(3)?, ball(4)?]), prec)
            .map_err(|er| format!("BallQr::lstsq: {er:?}"))?;
        let x_m = a
            .lstsq(&BallVec::from_balls(vec![ball(3)?, ball(4)?]), prec)
            .map_err(|er| format!("BallMat::lstsq: {er:?}"))?;
        for (what, x) in [("BallQr::lstsq", &x_f), ("BallMat::lstsq", &x_m)] {
            if !covers(x.get(0).ok_or("x empty")?, &e(1), &e(1))? {
                return Err(format!("prec {prec}: {what} x does not contain 1"));
            }
        }
    }
    Ok(())
}
