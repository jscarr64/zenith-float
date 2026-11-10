//! Slice S6 harness rows (`src/matrix.rs`). Owned by slice S6 alone. RED: the API is a stub until the green unit lands.
//! Goldens: tests/goldens/slices/matrix/ (its own SHA256SUMS and GENERATORS). Mutants: tools/mutants/manifest.d/matrix.txt, ids 900..=999.
//! Golden file headers for this slice are declared here, not in tests/common/goldens.rs.
//! CONTRACT 16.S6.0 to 16.S6.11. Oracle: exact Fraction golds from scripts/gen_matrix_goldens.py (gate file matrix_rows.txt; the heavy
//! file matrix_rows_heavy.txt is replayed by the ignored heavy test when ZB_MATRIX_HEAVY points at it), plus hand rows below.
//! COVERS: MagArray::from_mags, MagArray::get, MagArray::shape, MagArray::rank, MagArray::len, MagArray::is_empty, MagArray::as_slice, MagArray::into_parts, BallArray::add, BallArray::sub, BallArray::mid_rad, BallArray::from_mid_rad, BallVec::dot, BallVec::norm_inf, BallVec::norm_1, BallVec::norm_2_up, BallMat::matvec, BallMat::matmul, BallMat::norm_inf, BallMat::norm_1, BallMat::norm_frob_up, BallVec::add, BallVec::sub, BallMat::add, BallMat::sub

use crate::common::goldens::golden_dir;
use zenith_bound::array::{BallArray, BallMat, BallVec};
use zenith_bound::ball::Ball;
use zenith_bound::enclosure::Enclosure;
use zenith_bound::error::BoundError;
use zenith_bound::mag::Mag;
use zenith_bound::matrix::MagArray;
use zenith_bound::tri::Tri;
use zenith_float_num::{ExactNum, RoundingMode};

const GATE_ROWS: usize = 150;
const HEAVY_ROWS: usize = 30;
const TOP: u64 = 0x8000_0000_0000_0000;

// ---- helpers ------------------------------------------------------------------------------------------------------------
fn num(v: i64) -> ExactNum {
    ExactNum::from(v)
}

fn ok<T>(what: &str, r: Result<T, BoundError>) -> Result<T, String> {
    r.map_err(|e| format!("{what}: unexpected {e:?}"))
}

/// `2^k` exactly.
fn pow2(k: i64) -> Result<ExactNum, String> {
    let e = k.checked_add(1).ok_or("pow2 exponent")?;
    ExactNum::from_mag_parts(TOP, e, 64).map_err(|e| format!("pow2({k}): {e:?}"))
}

fn ball(v: i64, rad: Mag, prec: usize) -> Result<Ball, String> {
    ok("from_mid_rad", Ball::from_mid_rad(num(v), rad, prec))
}

fn pt(v: i64) -> Result<Ball, String> {
    ok("exact", Ball::exact(num(v), 64))
}

fn mat(nrows: usize, ncols: usize, vals: &[i64]) -> Result<BallMat, String> {
    let mut balls = Vec::new();
    for v in vals {
        balls.push(pt(*v)?);
    }
    ok(
        "BallMat::from_balls",
        BallMat::from_balls(nrows, ncols, balls),
    )
}

fn vec_of(vals: &[i64]) -> Result<BallVec, String> {
    let mut balls = Vec::new();
    for v in vals {
        balls.push(pt(*v)?);
    }
    Ok(BallVec::from_balls(balls))
}

fn sign_of(x: &ExactNum, y: &ExactNum) -> i8 {
    match x.partial_cmp(y) {
        Some(core::cmp::Ordering::Greater) => 1,
        Some(core::cmp::Ordering::Less) => -1,
        Some(core::cmp::Ordering::Equal) => 0,
        // incomparable (NaN-like) never satisfies an equality or a bound
        None => 2,
    }
}

fn want_dim_mismatch<T: core::fmt::Debug>(
    what: &str,
    r: Result<T, BoundError>,
) -> Result<(), String> {
    match r {
        Err(BoundError::DimensionMismatch { .. }) => Ok(()),
        other => Err(format!("{what}: want DimensionMismatch, got {other:?}")),
    }
}

fn want_contains(what: &str, b: &Ball, x: &ExactNum) -> Result<(), String> {
    match b.contains(x) {
        Tri::True => Ok(()),
        t @ (Tri::False | Tri::Unknown) => {
            Err(format!("{what}: {b:?} does not contain {x:?} ({t:?})"))
        }
    }
}

/// `x <= bound * (1 + 2^-56)` (the tightness slack of 16.S6.6 and 16.S6.9).
fn within_slack(what: &str, x: &Mag, bound: &ExactNum) -> Result<(), String> {
    let xe = x
        .to_exact_num()
        .map_err(|e| format!("{what}: radius {x:?} not finite: {e:?}"))?;
    let slack = bound.mul_full_prec(&pow2(-56)?);
    let limit = bound.add_full_prec(&slack);
    if sign_of(&xe, &limit) > 0 {
        return Err(format!("{what}: {x:?} exceeds {bound:?} * (1 + 2^-56)"));
    }
    Ok(())
}

/// `x <= 2 * tref` (the tightness factor of 16.S6.9; TREF = R1 + 2^(1-prec) S from the gold row).
fn within_twice(what: &str, x: &Mag, tref: &ExactNum) -> Result<(), String> {
    let xe = x
        .to_exact_num()
        .map_err(|e| format!("{what}: radius {x:?} not finite: {e:?}"))?;
    let limit = tref.add_full_prec(tref);
    if sign_of(&xe, &limit) > 0 {
        return Err(format!("{what}: {x:?} exceeds 2 * TREF {tref:?}"));
    }
    Ok(())
}

// ---- gold rows ------------------------------------------------------------------------------------------------------------
/// `SIGN EXP MANT` (scripts/gen_matrix_goldens.py): value = (-1)^SIGN · MANT · 2^(EXP − bitlen(MANT)), MANT hex.
fn gold_val(sign: &str, exp: &str, hex: &str) -> Result<ExactNum, String> {
    let e: i64 = exp.parse().map_err(|_| format!("bad exponent {exp:?}"))?;
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
            Err(format!("non-canonical zero {sign} {exp} {hex}"))
        };
    }
    let len = i64::try_from(hex.len()).map_err(|_| "mantissa length")?;
    let bitlen = 4 * (len - 1) + i64::from(32 - first.leading_zeros());
    let base = e - bitlen;
    let mut acc = num(0);
    let mut shift = 0i64;
    for bytes in hex.as_bytes().rchunks(16) {
        let chunk = std::str::from_utf8(bytes).map_err(|_| "utf8")?;
        let c = u64::from_str_radix(chunk, 16).map_err(|_| format!("bad hex {chunk:?}"))?;
        if c != 0 {
            let lz = i64::from(c.leading_zeros());
            let part = ExactNum::from_mag_parts(c << lz, base + shift + 64 - lz, 64)
                .map_err(|e| format!("from_mag_parts: {e:?}"))?;
            acc = acc.add_full_prec(&part);
        }
        shift += 64;
    }
    match sign {
        "0" => Ok(acc),
        "1" => Ok(num(0).sub_full_prec(&acc)),
        other => Err(format!("bad sign {other:?}")),
    }
}

struct Cursor<'a> {
    toks: Vec<&'a str>,
    at: usize,
}

impl<'a> Cursor<'a> {
    fn tok(&mut self) -> Result<&'a str, String> {
        let t = self.toks.get(self.at).ok_or("short row")?;
        self.at += 1;
        Ok(t)
    }
    fn usize(&mut self) -> Result<usize, String> {
        let t = self.tok()?;
        t.parse::<usize>()
            .map_err(|e| format!("bad count {t:?}: {e}"))
    }
    fn val(&mut self) -> Result<ExactNum, String> {
        let s = self.tok()?;
        let e = self.tok()?;
        let m = self.tok()?;
        gold_val(s, e, m)
    }
    fn balls(&mut self, n: usize, prec: usize) -> Result<Vec<Ball>, String> {
        let mut out = Vec::new();
        for _ in 0..n {
            let mid = self.val()?;
            let rad = Mag::from_exact_num(&self.val()?);
            out.push(ok("gold ball", Ball::from_mid_rad(mid, rad, prec))?);
        }
        Ok(out)
    }
    fn done(&self) -> Result<(), String> {
        if self.at == self.toks.len() {
            Ok(())
        } else {
            Err(format!("{} trailing tokens", self.toks.len() - self.at))
        }
    }
}

fn read_gold(name: &str) -> Result<String, String> {
    let path = golden_dir().join("slices").join("matrix").join(name);
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

/// Checks one product output entry: contains P0 P1 P2 (16.S6.8) and rad <= 2 TREF (16.S6.9). No midpoint is pinned (16.S6.7).
fn check_entry(what: &str, c: &Ball, cur: &mut Cursor) -> Result<(), String> {
    for k in 0..3 {
        want_contains(&format!("{what} P{k}"), c, &cur.val()?)?;
    }
    let tref = cur.val()?;
    within_twice(&format!("{what} rad"), c.rad(), &tref)
}

/// Replays one gold row through zenith-bound; Ok(()) when every check holds.
fn replay_row(line: &str) -> Result<(), String> {
    let mut cur = Cursor {
        toks: line.split(' ').collect(),
        at: 0,
    };
    let id = cur.tok()?;
    let op = cur.tok()?;
    let prec = cur.usize()?;
    let m = cur.usize()?;
    let k = cur.usize()?;
    let n = cur.usize()?;
    let size = |a: usize, b: usize| {
        a.checked_mul(b)
            .ok_or_else(|| format!("{id}: dims overflow"))
    };
    match op {
        "norm_inf" | "norm_1" | "norm_frob" | "vnorm_inf" | "vnorm_1" | "vnorm_2" => {
            let balls = cur.balls(size(m, n)?, prec)?;
            let s = cur.val()?;
            let hi = cur.val()?;
            let squared = op == "norm_frob" || op == "vnorm_2";
            let r = if op.starts_with('v') {
                if m != 1 {
                    return Err(format!("{id}: vector norm row with M = {m}"));
                }
                let v = BallVec::from_balls(balls);
                match op {
                    "vnorm_inf" => v.norm_inf(),
                    "vnorm_1" => v.norm_1(),
                    other => {
                        if other != "vnorm_2" {
                            return Err(format!("{id}: op {other}"));
                        }
                        v.norm_2_up()
                    }
                }
            } else {
                let a = ok(id, BallMat::from_balls(m, n, balls))?;
                match op {
                    "norm_inf" => a.norm_inf(),
                    "norm_1" => a.norm_1(),
                    other => {
                        if other != "norm_frob" {
                            return Err(format!("{id}: op {other}"));
                        }
                        a.norm_frob_up()
                    }
                }
            };
            let re = r
                .to_exact_num()
                .map_err(|e| format!("{id}: norm {r:?} not finite: {e:?}"))?;
            let lhs = if squared { re.mul_full_prec(&re) } else { re };
            if sign_of(&lhs, &s) < 0 {
                return Err(format!("{id}: {op} {r:?} is below the exact value"));
            }
            within_slack(id, &r, &hi)?;
        }
        "dot" | "matvec" | "matmul" => {
            let av = cur.balls(size(m, k)?, prec)?;
            let bv = cur.balls(size(k, n)?, prec)?;
            match op {
                "dot" => {
                    let c = ok(
                        id,
                        BallVec::from_balls(av).dot(&BallVec::from_balls(bv), prec),
                    )?;
                    check_entry(id, &c, &mut cur)?;
                }
                "matvec" => {
                    let a = ok(id, BallMat::from_balls(m, k, av))?;
                    let y = ok(id, a.matvec(&BallVec::from_balls(bv), prec))?;
                    for i in 0..m {
                        let c = y.get(i).ok_or_else(|| format!("{id}: no y[{i}]"))?;
                        check_entry(&format!("{id} y[{i}]"), c, &mut cur)?;
                    }
                }
                other => {
                    if other != "matmul" {
                        return Err(format!("{id}: op {other}"));
                    }
                    let a = ok(id, BallMat::from_balls(m, k, av))?;
                    let b = ok(id, BallMat::from_balls(k, n, bv))?;
                    let c = ok(id, a.matmul(&b, prec))?;
                    if c.nrows() != m || c.ncols() != n {
                        return Err(format!("{id}: result is {}x{}", c.nrows(), c.ncols()));
                    }
                    for i in 0..m {
                        for j in 0..n {
                            let e = c.get(i, j).ok_or_else(|| format!("{id}: no c[{i},{j}]"))?;
                            check_entry(&format!("{id} c[{i},{j}]"), e, &mut cur)?;
                        }
                    }
                }
            }
        }
        other => return Err(format!("{id}: unknown op {other:?}")),
    }
    cur.done().map_err(|e| format!("{id}: {e}"))
}

/// Replays every row; returns the number of rows and the first few failures.
fn replay_text(text: &str, want_rows: usize) -> Result<(), String> {
    let rows = header_rows(text)?;
    if rows != want_rows {
        return Err(format!("header says {rows} rows, slice pins {want_rows}"));
    }
    let mut seen = 0usize;
    let mut failures = Vec::new();
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        seen += 1;
        if let Err(e) = replay_row(line) {
            failures.push(e);
        }
    }
    if seen != rows {
        return Err(format!("{seen} data rows, header says {rows}"));
    }
    if failures.is_empty() {
        Ok(())
    } else {
        let first: Vec<&String> = failures.iter().take(5).collect();
        Err(format!(
            "{} of {rows} rows fail; first: {first:#?}",
            failures.len()
        ))
    }
}

// ---- gold tests -----------------------------------------------------------------------------------------------------------
#[test]
fn s6_gate_rows_match_header_and_pinned_count() -> Result<(), String> {
    let text = read_gold("matrix_rows.txt")?;
    let rows = header_rows(&text)?;
    let seen = text.lines().filter(|l| !l.starts_with('#')).count();
    if rows == GATE_ROWS && seen == GATE_ROWS {
        Ok(())
    } else {
        Err(format!("header {rows}, data {seen}, pinned {GATE_ROWS}"))
    }
}

#[test]
fn s6_gate_rows_replay_through_zb() -> Result<(), String> {
    replay_text(&read_gold("matrix_rows.txt")?, GATE_ROWS)
}

#[test]
#[ignore = "heavy tier (16.S6.10): set ZB_MATRIX_HEAVY to matrix_rows_heavy.txt (scripts/gen_matrix_goldens.py --heavy DIR; sha256 in HEAVY.SHA256)"]
fn s6_heavy_rows_replay_through_zb() -> Result<(), String> {
    let path = std::env::var("ZB_MATRIX_HEAVY").map_err(|e| format!("ZB_MATRIX_HEAVY: {e}"))?;
    let text = std::fs::read_to_string(&path).map_err(|e| format!("cannot read {path}: {e}"))?;
    replay_text(&text, HEAVY_ROWS)
}

// ---- hand rows ------------------------------------------------------------------------------------------------------------
#[test]
fn s6_matmul_hand_row_2x2() -> Result<(), String> {
    // [1 2; 3 4] * [5 6; 7 8] = [19 22; 43 50], all exact
    let c = ok(
        "matmul",
        mat(2, 2, &[1, 2, 3, 4])?.matmul(&mat(2, 2, &[5, 6, 7, 8])?, 53),
    )?;
    for (i, j, v) in [(0, 0, 19), (0, 1, 22), (1, 0, 43), (1, 1, 50)] {
        let e = c.get(i, j).ok_or("missing entry")?;
        if sign_of(e.mid(), &num(v)) != 0 || !e.rad().is_zero() {
            return Err(format!("c[{i},{j}] = {e:?}, want {v} exactly"));
        }
    }
    Ok(())
}

#[test]
fn s6_shape_mismatch_is_an_error_not_a_panic() -> Result<(), String> {
    let a = mat(2, 3, &[1, 2, 3, 4, 5, 6])?;
    want_dim_mismatch(
        "matmul 2x3 * 2x3",
        a.matmul(&mat(2, 3, &[1, 2, 3, 4, 5, 6])?, 53),
    )?;
    want_dim_mismatch("matvec 2x3 * 2", a.matvec(&vec_of(&[1, 2])?, 53))?;
    want_dim_mismatch("dot 2 . 3", vec_of(&[1, 2])?.dot(&vec_of(&[1, 2, 3])?, 53))?;
    let six = vec_of(&[1, 2, 3, 4, 5, 6])?.into_array();
    want_dim_mismatch("add [6] + [2,3]", six.add(a.as_array(), 53))?;
    let six = vec_of(&[1, 2, 3, 4, 5, 6])?.into_array();
    want_dim_mismatch("sub [6] - [2,3]", six.sub(a.as_array(), 53))?;
    let rads = ok(
        "from_mags",
        MagArray::from_mags(Box::new([2]), vec![Mag::zero(), Mag::zero()]),
    )?;
    want_dim_mismatch(
        "from_mid_rad [6] with [2]",
        BallArray::from_mid_rad(a.into_array(), &rads, 53),
    )
}

#[test]
fn s6_bad_prec_is_invalid_argument() -> Result<(), String> {
    match mat(1, 1, &[1])?.matmul(&mat(1, 1, &[1])?, 0) {
        Err(BoundError::InvalidArgument { .. }) => Ok(()),
        other => Err(format!("prec 0: {other:?}")),
    }
}

#[test]
fn s6_empty_inner_dimension_gives_exact_zeros() -> Result<(), String> {
    let a = mat(2, 0, &[])?;
    let b = mat(0, 3, &[])?;
    let c = ok("matmul 2x0 * 0x3", a.matmul(&b, 53))?;
    if c.nrows() != 2 || c.ncols() != 3 {
        return Err(format!("shape {}x{}", c.nrows(), c.ncols()));
    }
    for e in c.as_array().as_slice() {
        if !e.mid().is_zero() || !e.rad().is_zero() {
            return Err(format!("entry {e:?} is not 0 +- 0"));
        }
    }
    let d = ok("dot of empty", vec_of(&[])?.dot(&vec_of(&[])?, 53))?;
    if d.mid().is_zero() && d.rad().is_zero() {
        Ok(())
    } else {
        Err(format!("empty dot {d:?}"))
    }
}

#[test]
fn s6_infinite_radius_is_all_reals_entry_and_infinite_norm() -> Result<(), String> {
    let inf = ball(1, Mag::infinity(), 53)?;
    let a = ok("from_balls", BallMat::from_balls(2, 1, vec![inf, pt(2)?]))?;
    let b = mat(1, 2, &[3, 4])?;
    let c = ok("matmul with an Infinity radius", a.matmul(&b, 53))?;
    for j in 0..2 {
        let top = c.get(0, j).ok_or("missing")?;
        if !top.rad().is_infinity() {
            return Err(format!("c[0,{j}] {top:?} should be all-reals"));
        }
        let low = c.get(1, j).ok_or("missing")?;
        if low.rad().is_infinity() {
            return Err(format!("c[1,{j}] {low:?} should stay finite (elementwise)"));
        }
    }
    if a.norm_inf().is_infinity() && a.norm_1().is_infinity() && a.norm_frob_up().is_infinity() {
        Ok(())
    } else {
        Err("norms of a matrix with an Infinity radius must be Infinity".to_string())
    }
}

#[test]
fn s6_norms_hand_row() -> Result<(), String> {
    // A = [1 -2; 3 4]: norm_inf = 7, norm_1 = 6, frob^2 = 30
    let a = mat(2, 2, &[1, -2, 3, 4])?;
    let exact = |m: &Mag| m.to_exact_num().map_err(|e| format!("{e:?}"));
    if sign_of(&exact(&a.norm_inf())?, &num(7)) != 0 || sign_of(&exact(&a.norm_1())?, &num(6)) != 0
    {
        return Err(format!(
            "norm_inf {:?}, norm_1 {:?}",
            a.norm_inf(),
            a.norm_1()
        ));
    }
    let f = exact(&a.norm_frob_up())?;
    let f2 = f.mul_full_prec(&f);
    if sign_of(&f2, &num(30)) < 0 {
        return Err("norm_frob_up below sqrt(30)".to_string());
    }
    within_slack("frob", &a.norm_frob_up(), &num(6))
}

#[test]
fn s6_mag_array_constructor_and_get() -> Result<(), String> {
    let mags = vec![
        Mag::zero(),
        Mag::infinity(),
        Mag::zero(),
        Mag::zero(),
        Mag::zero(),
        Mag::zero(),
    ];
    let m = ok("from_mags", MagArray::from_mags(Box::new([2, 3]), mags))?;
    if m.shape() != [2, 3] || m.rank() != 2 || m.len() != 6 {
        return Err(format!("shape {:?}", m.shape()));
    }
    let hit = m.get(&[0, 1]).ok_or("get [0,1]")?;
    if !hit.is_infinity() || m.get(&[1]).is_some() || m.get(&[2, 0]).is_some() {
        return Err(
            "get: wrong element, or a wrong rank / out-of-range index was accepted".to_string(),
        );
    }
    match MagArray::from_mags(Box::new([4]), vec![Mag::zero()]) {
        Err(BoundError::InvalidArgument { .. }) => Ok(()),
        other => Err(format!("length mismatch: {other:?}")),
    }
}

#[test]
fn s6_mid_rad_round_trip_encloses() -> Result<(), String> {
    let r = Mag::finite(TOP, -3).map_err(|e| format!("{e:?}"))?;
    let balls = vec![
        ball(5, Mag::finite(TOP, -3).map_err(|e| format!("{e:?}"))?, 53)?,
        pt(-7)?,
    ];
    let arr = ok("from_balls", BallArray::from_balls(Box::new([2]), balls))?;
    let (mids, rads) = ok("mid_rad", arr.mid_rad())?;
    if rads.shape() != [2] || mids.shape() != [2] {
        return Err(format!(
            "split shapes {:?} {:?}",
            mids.shape(),
            rads.shape()
        ));
    }
    for b in mids.as_slice() {
        if !b.is_exact() {
            return Err(format!("mid {b:?} not exact"));
        }
    }
    let first = rads.as_slice().first().ok_or("no rad")?;
    let want = r.to_exact_num().map_err(|e| format!("{e:?}"))?;
    if first.to_exact_num().ok().map(|x| sign_of(&x, &want)) != Some(0) {
        return Err(format!("rad[0] {first:?}"));
    }
    let back = ok("from_mid_rad", BallArray::from_mid_rad(mids, &rads, 53))?;
    let b0 = back.as_slice().first().ok_or("empty")?;
    want_contains("round trip lower", b0, &num(5).sub_full_prec(&pow2(-4)?))?;
    want_contains("round trip upper", b0, &num(5).add_full_prec(&pow2(-4)?))
}

#[test]
fn s6_elementwise_add_hand_row() -> Result<(), String> {
    let a = vec_of(&[1, 2, 3])?.into_array();
    let b = vec_of(&[10, 20, 30])?.into_array();
    let c = ok("add", a.add(&b, 53))?;
    for (e, v) in c.as_slice().iter().zip([11, 22, 33]) {
        if sign_of(e.mid(), &num(v)) != 0 || !e.rad().is_zero() {
            return Err(format!("{e:?} != {v}"));
        }
    }
    Ok(())
}

#[test]
fn s6_ballvec_norms_hand_row() -> Result<(), String> {
    // v = (3, -4, 12): norm_inf = 12, norm_1 = 19, norm_2 = 13 exactly (3^2 + 4^2 + 12^2 = 169)
    let v = vec_of(&[3, -4, 12])?;
    let exact = |m: &Mag| m.to_exact_num().map_err(|e| format!("{e:?}"));
    if sign_of(&exact(&v.norm_inf())?, &num(12)) != 0
        || sign_of(&exact(&v.norm_1())?, &num(19)) != 0
    {
        return Err(format!(
            "norm_inf {:?}, norm_1 {:?}",
            v.norm_inf(),
            v.norm_1()
        ));
    }
    let f = exact(&v.norm_2_up())?;
    if sign_of(&f.mul_full_prec(&f), &num(169)) < 0 {
        return Err(format!("norm_2_up {f:?} below 13"));
    }
    within_slack("norm_2_up", &v.norm_2_up(), &num(13))?;
    // radii count (mutants 914, 915): w = (-1 +- 2^-4, 1): |m| + r = (17/16, 1); inf 17/16, 1-norm 33/16, 2-norm^2 = 545/256
    let w = BallVec::from_balls(vec![
        ball(-1, Mag::finite(TOP, -3).map_err(|e| format!("{e:?}"))?, 53)?,
        pt(1)?,
    ]);
    let sixteenth = pow2(-4)?;
    let want_inf = num(17).mul_full_prec(&sixteenth);
    let want_1 = num(33).mul_full_prec(&sixteenth);
    if sign_of(&exact(&w.norm_inf())?, &want_inf) != 0
        || sign_of(&exact(&w.norm_1())?, &want_1) != 0
    {
        return Err(format!(
            "radius: norm_inf {:?} (want 17/16), norm_1 {:?} (want 33/16)",
            w.norm_inf(),
            w.norm_1()
        ));
    }
    let g = exact(&w.norm_2_up())?;
    if sign_of(&g.mul_full_prec(&g), &num(545).mul_full_prec(&pow2(-8)?)) < 0 {
        return Err(format!("radius: norm_2_up {g:?} below sqrt(545)/16"));
    }
    // empty vector: every norm is exactly 0; an Infinity radius makes every norm Infinity
    let e = vec_of(&[])?;
    if !e.norm_inf().is_zero() || !e.norm_1().is_zero() || !e.norm_2_up().is_zero() {
        return Err("norms of the empty vector must be 0".to_string());
    }
    let inf = BallVec::from_balls(vec![pt(2)?, ball(1, Mag::infinity(), 53)?]);
    if inf.norm_inf().is_infinity() && inf.norm_1().is_infinity() && inf.norm_2_up().is_infinity() {
        Ok(())
    } else {
        Err("norms of a vector with an Infinity radius must be Infinity".to_string())
    }
}

/// Spence ruling b (2026-10-25): a K = 1 dot is the S1 product plus its own rounding. With `w = prec + 3` the midpoint is
/// `round(round(m(a)·m(b), w), prec)`; when the product fits `w` the dot is identical to `Ball::mul`; otherwise it still contains the
/// exact product box and its radius stays within 2·TREF (16.S6.9).
#[test]
fn s6_k1_dot_is_s1_mul_plus_own_rounding() -> Result<(), String> {
    let half = |e: i64| Mag::finite(TOP, e).map_err(|x| format!("{x:?}"));
    // 15 fits w = 5 bits but not prec = 2: identical to S1 mul (mid and radius)
    let d = ok(
        "dot",
        BallVec::from_balls(vec![ball(3, half(-9)?, 2)?])
            .dot(&BallVec::from_balls(vec![ball(5, half(-11)?, 2)?]), 2),
    )?;
    let s1 = ok(
        "Ball::mul",
        ball(3, half(-9)?, 2)?.mul(&ball(5, half(-11)?, 2)?, 2),
    )?;
    if sign_of(d.mid(), s1.mid()) != 0 || d.rad() != s1.rad() {
        return Err(format!("short K=1 dot {d:?} != S1 mul {s1:?}"));
    }
    // 255 * 253 = 64515 needs 16 bits: rounded at w = 11, then at prec = 8
    let a = ball(255, half(-19)?, 8)?;
    let d = ok(
        "dot",
        BallVec::from_balls(vec![a]).dot(&BallVec::from_balls(vec![pt(253)?]), 8),
    )?;
    let want = num(255)
        .mul(&num(253), 11, RoundingMode::ToEven)
        .ldexp(0, 8, RoundingMode::ToEven);
    if sign_of(d.mid(), &want) != 0 {
        return Err(format!(
            "long K=1 dot mid {:?} != round(round(64515, 11), 8) {want:?}",
            d.mid()
        ));
    }
    let lo = num(255).sub_full_prec(&pow2(-20)?).mul_full_prec(&num(253));
    let hi = num(255).add_full_prec(&pow2(-20)?).mul_full_prec(&num(253));
    want_contains("K=1 P0", &d, &num(64515))?;
    want_contains("K=1 low corner", &d, &lo)?;
    want_contains("K=1 high corner", &d, &hi)?;
    let tref = num(253)
        .mul_full_prec(&pow2(-20)?)
        .add_full_prec(&num(64515).mul_full_prec(&pow2(-7)?));
    within_twice("K=1 radius", d.rad(), &tref)
}

/// 16.S6.5: `from_mid_rad` keeps the midpoint ball's own radius `s_i` and adds `r_i` (`m_i ± (s_i + r_i + ε)`), so a non-exact
/// midpoint is never shrunk to `± r_i`.
#[test]
fn s6_from_mid_rad_keeps_mid_radius() -> Result<(), String> {
    let s = Mag::finite(TOP, -3).map_err(|e| format!("{e:?}"))?;
    let r = Mag::finite(TOP, -5).map_err(|e| format!("{e:?}"))?;
    let mids = ok(
        "from_balls",
        BallArray::from_balls(Box::new([1]), vec![ball(5, s, 53)?]),
    )?;
    let rads = ok("from_mags", MagArray::from_mags(Box::new([1]), vec![r]))?;
    let back = ok("from_mid_rad", BallArray::from_mid_rad(mids, &rads, 53))?;
    let b0 = back.as_slice().first().ok_or("empty")?;
    let reach = pow2(-4)?.add_full_prec(&pow2(-6)?);
    want_contains("s + r lower", b0, &num(5).sub_full_prec(&reach))?;
    want_contains("s + r upper", b0, &num(5).add_full_prec(&reach))
}

// ---- MagArray accessors (shape, rank, len, is_empty, as_slice, into_parts) ------------------------------------------------
/// The distinct finite magnitudes `1, 2, ..., n` (each `v * 2^-64`), in that order.
fn mags_1_to(n: u64) -> Result<Vec<Mag>, String> {
    (1..=n)
        .map(|v| Mag::finite(v, 0).map_err(|e| format!("finite({v}): {e:?}")))
        .collect()
}

/// Checks every accessor of `m` against the expected shape and the expected row-major elements `1..=n`.
fn want_accessors(what: &str, m: &MagArray, shape: &[usize], n: u64) -> Result<(), String> {
    let want = mags_1_to(n)?;
    let product = shape
        .iter()
        .try_fold(1usize, |a, d| a.checked_mul(*d))
        .ok_or("product")?;
    if m.shape() != shape {
        return Err(format!("{what}: shape {:?}, want {shape:?}", m.shape()));
    }
    if m.rank() != shape.len() {
        return Err(format!("{what}: rank {}, want {}", m.rank(), shape.len()));
    }
    if m.len() != product || m.len() != want.len() {
        return Err(format!("{what}: len {}, want {product}", m.len()));
    }
    if m.is_empty() != (product == 0) {
        return Err(format!("{what}: is_empty {}", m.is_empty()));
    }
    if m.as_slice() != want.as_slice() {
        return Err(format!(
            "{what}: as_slice {:?}, want {want:?}",
            m.as_slice()
        ));
    }
    Ok(())
}

#[test]
fn s6_magarray_rank0_scalar_accessors() -> Result<(), String> {
    let m = ok(
        "from_mags",
        MagArray::from_mags(Box::new([]), mags_1_to(1)?),
    )?;
    want_accessors("rank 0", &m, &[], 1)?;
    if m.rank() != 0 || m.len() != 1 || m.is_empty() || !m.shape().is_empty() {
        return Err("rank 0: a scalar has rank 0, one element and an empty shape".to_string());
    }
    Ok(())
}

#[test]
fn s6_magarray_1d_accessors_keep_order() -> Result<(), String> {
    let m = ok(
        "from_mags",
        MagArray::from_mags(Box::new([3]), mags_1_to(3)?),
    )?;
    want_accessors("[3]", &m, &[3], 3)?;
    let last = m.as_slice().last().ok_or("last")?;
    if *last != Mag::finite(3, 0).map_err(|e| format!("{e:?}"))? || m.len() != 3 || m.rank() != 1 {
        return Err(format!("[3]: last {last:?}"));
    }
    Ok(())
}

#[test]
fn s6_magarray_nd_accessors_row_major() -> Result<(), String> {
    let m = ok(
        "from_mags",
        MagArray::from_mags(Box::new([2, 3, 4]), mags_1_to(24)?),
    )?;
    want_accessors("[2,3,4]", &m, &[2, 3, 4], 24)?;
    if m.len() != 24 || m.rank() != 3 || m.shape() != [2, 3, 4] {
        return Err("[2,3,4]: len 24, rank 3".to_string());
    }
    // Row-major: index [1,2,3] is offset 1*12 + 2*4 + 3 = 23, the 24th element; [0,1,0] is offset 4, the 5th.
    let at = |idx: &[usize], v: u64| -> Result<(), String> {
        let got = m.get(idx).ok_or("get")?;
        let flat = m
            .as_slice()
            .get(usize::try_from(v - 1).map_err(|e| format!("{e}"))?)
            .ok_or("nth")?;
        if *got != Mag::finite(v, 0).map_err(|e| format!("{e:?}"))? || got != flat {
            return Err(format!("[2,3,4] {idx:?}: {got:?}, want {v}"));
        }
        Ok(())
    };
    at(&[1, 2, 3], 24)?;
    at(&[0, 1, 0], 5)
}

#[test]
fn s6_magarray_zero_length_dims_are_empty() -> Result<(), String> {
    let a = ok("from_mags", MagArray::from_mags(Box::new([0]), Vec::new()))?;
    want_accessors("[0]", &a, &[0], 0)?;
    let b = ok(
        "from_mags",
        MagArray::from_mags(Box::new([2, 0, 3]), Vec::new()),
    )?;
    want_accessors("[2,0,3]", &b, &[2, 0, 3], 0)?;
    if !a.is_empty() || a.rank() != 1 || !b.is_empty() || b.rank() != 3 {
        return Err("zero-length dims: empty, len 0, rank kept".to_string());
    }
    if !a.as_slice().is_empty() || !b.as_slice().is_empty() {
        return Err("zero-length dims: as_slice not empty".to_string());
    }
    Ok(())
}

#[test]
fn s6_magarray_into_parts_round_trips() -> Result<(), String> {
    let m = ok(
        "from_mags",
        MagArray::from_mags(Box::new([2, 3]), mags_1_to(6)?),
    )?;
    let (shape, data) = m.into_parts();
    if *shape != [2, 3] {
        return Err(format!("into_parts shape {shape:?}"));
    }
    if *data != *mags_1_to(6)?.into_boxed_slice() {
        return Err(format!("into_parts data {data:?}"));
    }
    let again = ok(
        "from_mags again",
        MagArray::from_mags(shape, data.into_vec()),
    )?;
    want_accessors("rebuilt [2,3]", &again, &[2, 3], 6)?;
    let (s0, d0) = ok(
        "from_mags",
        MagArray::from_mags(Box::new([]), mags_1_to(1)?),
    )?
    .into_parts();
    if !s0.is_empty() || *d0 != *mags_1_to(1)?.into_boxed_slice() {
        return Err(format!("rank-0 into_parts {s0:?} {d0:?}"));
    }
    Ok(())
}

#[test]
fn s6_magarray_from_mags_rejects_bad_shapes() -> Result<(), String> {
    let cases: [(&str, Box<[usize]>, u64); 5] = [
        ("rank 0 with no element", Box::new([]), 0),
        ("rank 0 with two elements", Box::new([]), 2),
        ("[2,3] with 5", Box::new([2, 3]), 5),
        ("[0] with 1", Box::new([0]), 1),
        ("[usize::MAX, 2] overflow", Box::new([usize::MAX, 2]), 0),
    ];
    for (what, shape, n) in cases {
        match MagArray::from_mags(shape, mags_1_to(n)?) {
            Err(BoundError::InvalidArgument {
                op: "MagArray::from_mags",
                ..
            }) => {}
            other => return Err(format!("{what}: {other:?}")),
        }
    }
    Ok(())
}

// ---- elementwise add/sub on BallVec and BallMat (16.S6.4; pubfn coverage) ---------------------------------------------
/// A hand entry `numer·2^k ± 2^rad_k` (`None`: radius 0). Every value is dyadic and exact at 53 bits, so the exact sums and
/// differences below are hand values (Fraction-style), not golds.
type Spec = (i64, i64, Option<i64>);

const VA: [Spec; 3] = [(5, -2, Some(-4)), (-3, 0, None), (7, -3, Some(-2))];
const VB: [Spec; 3] = [(1, -1, Some(-4)), (3, -2, Some(-2)), (-7, -3, None)];
const MA: [Spec; 4] = [
    (5, -2, Some(-4)),
    (-3, 0, None),
    (7, -3, Some(-2)),
    (9, 0, Some(-1)),
];
const MB: [Spec; 4] = [
    (1, -1, Some(-4)),
    (3, -2, Some(-2)),
    (-7, -3, None),
    (-1, -3, Some(-3)),
];

fn spec_mid(s: &Spec) -> Result<ExactNum, String> {
    Ok(num(s.0).mul_full_prec(&pow2(s.1)?))
}

fn spec_rad(s: &Spec) -> Result<ExactNum, String> {
    match s.2 {
        Some(k) => pow2(k),
        None => Ok(num(0)),
    }
}

fn spec_balls(specs: &[Spec]) -> Result<Vec<Ball>, String> {
    let mut out = Vec::new();
    for s in specs {
        let r = match s.2 {
            Some(k) => Mag::finite(TOP, k.checked_add(1).ok_or("rad exponent")?)
                .map_err(|e| format!("{e:?}"))?,
            None => Mag::zero(),
        };
        out.push(ok("from_mid_rad", Ball::from_mid_rad(spec_mid(s)?, r, 53))?);
    }
    Ok(out)
}

/// Entry `i` of `got` against `a_i + b_i` (`neg`: `a_i - b_i`): the exact midpoint, a radius no worse than the 16.S1.1 bound
/// `r_a + r_b` (no rounding term: every midpoint sum is exact at 53 bits), and containment of all four corner values.
fn check_elementwise(
    what: &str,
    got: &[Ball],
    a: &[Spec],
    b: &[Spec],
    neg: bool,
) -> Result<(), String> {
    if got.len() != a.len() || a.len() != b.len() {
        return Err(format!("{what}: {} entries, want {}", got.len(), a.len()));
    }
    let op = |x: &ExactNum, y: &ExactNum| {
        if neg {
            x.sub_full_prec(y)
        } else {
            x.add_full_prec(y)
        }
    };
    for (i, ((g, x), y)) in got.iter().zip(a).zip(b).enumerate() {
        let (mx, rx, my, ry) = (spec_mid(x)?, spec_rad(x)?, spec_mid(y)?, spec_rad(y)?);
        let mid = op(&mx, &my);
        if sign_of(g.mid(), &mid) != 0 {
            return Err(format!("{what}[{i}]: midpoint {g:?}, want {mid:?}"));
        }
        let ge = g
            .rad()
            .to_exact_num()
            .map_err(|e| format!("{what}[{i}]: radius not finite: {e:?}"))?;
        let bound = rx.add_full_prec(&ry);
        if sign_of(&ge, &bound) > 0 {
            return Err(format!(
                "{what}[{i}]: radius {ge:?} exceeds r_a + r_b = {bound:?}"
            ));
        }
        let xs = [mx.sub_full_prec(&rx), mx.add_full_prec(&rx)];
        let ys = [my.sub_full_prec(&ry), my.add_full_prec(&ry)];
        for cx in &xs {
            for cy in &ys {
                want_contains(&format!("{what}[{i}] corner"), g, &op(cx, cy))?;
            }
        }
    }
    Ok(())
}

/// An Infinity radius in entry 0 gives the all-reals entry there and nowhere else (elementwise).
fn want_all_reals_first_only(what: &str, got: &[Ball]) -> Result<(), String> {
    let mut it = got.iter();
    let first = it.next().ok_or("empty result")?;
    if !first.rad().is_infinity() {
        return Err(format!("{what}: entry 0 {first:?} should be all-reals"));
    }
    for e in it {
        if e.rad().is_infinity() {
            return Err(format!("{what}: {e:?} should stay finite"));
        }
    }
    Ok(())
}

#[test]
fn s6_ballvec_add_sub_hand_rows() -> Result<(), String> {
    let va = || Ok::<BallVec, String>(BallVec::from_balls(spec_balls(&VA)?));
    let vb = BallVec::from_balls(spec_balls(&VB)?);
    let s = ok("BallVec::add", va()?.add(&vb, 53))?;
    check_elementwise("BallVec::add", s.as_array().as_slice(), &VA, &VB, false)?;
    let d = ok("BallVec::sub", va()?.sub(&vb, 53))?;
    check_elementwise("BallVec::sub", d.as_array().as_slice(), &VA, &VB, true)?;
    want_dim_mismatch(
        "BallVec::add 3 + 2",
        vec_of(&[1, 2, 3])?.add(&vec_of(&[1, 2])?, 53),
    )?;
    want_dim_mismatch(
        "BallVec::sub 3 - 2",
        vec_of(&[1, 2, 3])?.sub(&vec_of(&[1, 2])?, 53),
    )?;
    let inf = || {
        Ok::<BallVec, String>(BallVec::from_balls(vec![
            ball(1, Mag::infinity(), 53)?,
            pt(2)?,
        ]))
    };
    let w = vec_of(&[3, 4])?;
    let s = ok("BallVec::add Infinity", inf()?.add(&w, 53))?;
    want_all_reals_first_only("BallVec::add", s.as_array().as_slice())?;
    let d = ok("BallVec::sub Infinity", inf()?.sub(&w, 53))?;
    want_all_reals_first_only("BallVec::sub", d.as_array().as_slice())
}

#[test]
fn s6_ballmat_add_sub_hand_rows() -> Result<(), String> {
    let ma = || {
        ok(
            "BallMat::from_balls",
            BallMat::from_balls(2, 2, spec_balls(&MA)?),
        )
    };
    let mb = ok(
        "BallMat::from_balls",
        BallMat::from_balls(2, 2, spec_balls(&MB)?),
    )?;
    let s = ok("BallMat::add", ma()?.add(&mb, 53))?;
    if s.nrows() != 2 || s.ncols() != 2 {
        return Err(format!("BallMat::add shape {}x{}", s.nrows(), s.ncols()));
    }
    check_elementwise("BallMat::add", s.as_array().as_slice(), &MA, &MB, false)?;
    let d = ok("BallMat::sub", ma()?.sub(&mb, 53))?;
    if d.nrows() != 2 || d.ncols() != 2 {
        return Err(format!("BallMat::sub shape {}x{}", d.nrows(), d.ncols()));
    }
    check_elementwise("BallMat::sub", d.as_array().as_slice(), &MA, &MB, true)?;
    // equal length, different shape (2x3 vs 3x2), and different length (2x2 vs 2x3)
    let six = [1, 2, 3, 4, 5, 6];
    want_dim_mismatch(
        "BallMat::add 2x3 + 3x2",
        mat(2, 3, &six)?.add(&mat(3, 2, &six)?, 53),
    )?;
    want_dim_mismatch(
        "BallMat::sub 2x3 - 3x2",
        mat(2, 3, &six)?.sub(&mat(3, 2, &six)?, 53),
    )?;
    want_dim_mismatch(
        "BallMat::add 2x2 + 2x3",
        mat(2, 2, &[1, 2, 3, 4])?.add(&mat(2, 3, &six)?, 53),
    )?;
    want_dim_mismatch(
        "BallMat::sub 2x2 - 2x3",
        mat(2, 2, &[1, 2, 3, 4])?.sub(&mat(2, 3, &six)?, 53),
    )?;
    let inf = || {
        let balls = vec![ball(1, Mag::infinity(), 53)?, pt(2)?];
        ok("BallMat::from_balls", BallMat::from_balls(1, 2, balls))
    };
    let w = mat(1, 2, &[3, 4])?;
    let s = ok("BallMat::add Infinity", inf()?.add(&w, 53))?;
    want_all_reals_first_only("BallMat::add", s.as_array().as_slice())?;
    let d = ok("BallMat::sub Infinity", inf()?.sub(&w, 53))?;
    want_all_reals_first_only("BallMat::sub", d.as_array().as_slice())
}
