// HOLD H12: held, not compiled (no `mod` line names this file). Index and use-or-delete rule: HOLD.md at the repo root.
// The original text follows, each line commented out with a leading `// `; uncomment it whole to go live.
// //! Slice LA2 red rows that need `ln` (16.LA1.6: `logdet` returns `Resource "needs_ln"` until S3 exports ln) and the
// //! format-2 LA golds (16.LA.5), with the helpers only they use. Held (H12) because they are red (P12: no ignored reds
// //! in a PR). Moved verbatim from slice_lu.rs; a `mod` line in slice_lu.rs takes them live, without their `#[ignore]`s.
//
// use super::*;
//
// /// Heavy `enc4096`/`enc16384` enclosure of a row from `$ZB_LA_HEAVY_DIR/la_logdet_enc<b>_heavy.txt` (last two tokens lo hi).
// fn heavy_enclosure(id: &str, bits: u32) -> Result<Option<(ExactNum, ExactNum)>, String> {
//     let Ok(dir) = std::env::var("ZB_LA_HEAVY_DIR") else {
//         return Ok(None);
//     };
//     let path = std::path::Path::new(&dir).join(format!("la_logdet_enc{bits}_heavy.txt"));
//     let Ok(text) = std::fs::read_to_string(&path) else {
//         return Ok(None);
//     };
//     for line in text.lines().filter(|l| !l.starts_with('#')) {
//         let toks: Vec<&str> = line.split_whitespace().collect();
//         if toks.first() == Some(&id) {
//             let n = toks.len();
//             let lo = dyadic(toks.get(n.checked_sub(2).ok_or("short")?).ok_or("lo")?)?;
//             let hi = dyadic(toks.get(n.checked_sub(1).ok_or("short")?).ok_or("hi")?)?;
//             return Ok(Some((lo, hi)));
//         }
//     }
//     Ok(None)
// }
//
// /// `(covers, overlaps)` of the ball against `[lo, hi]`.
// fn vs_interval(b: &Ball, lo: &ExactNum, hi: &ExactNum) -> Result<(bool, bool), String> {
//     if b.rad().is_infinity() {
//         return Ok((true, true));
//     }
//     let r = b.rad().to_exact_num().map_err(|e| format!("{e:?}"))?;
//     let m = b.mid();
//     let d_lo = m.sub_full_prec(lo);
//     let d_hi = hi.sub_full_prec(m);
//     let covers = d_lo <= r && d_hi <= r;
//     let neg_r = -plus_zero(&r)?;
//     let overlaps = d_lo >= neg_r && d_hi >= neg_r;
//     Ok((covers, overlaps))
// }
//
// /// The full LA4.2 replay: sign exact, det encloses det, `ln_abs_det` contains the exact ln|det| (16.LA.5 containment).
// #[test]
// #[ignore = "red: not implemented: CONTRACT LA4.2: ln_abs_det needs S3 ln (16.LA1.6: logdet returns Resource needs_ln until S3 merges); un-ignore when src/elementary.rs exports ln"]
// fn la2_logdet_gate_rows_replay_through_zb() -> Result<(), String> {
//     let text = read_gold("cholesky", "la_logdet.txt")?;
//     for prec in PRECS {
//         for toks in &data_rows(&text)? {
//             let r = logdet_row(toks)?;
//             let what = format!("{}@{prec} logdet", r.id);
//             let a = ball_mat(r.n, r.n, &r.a, prec)?;
//             let got = a.logdet(prec);
//             let Some((sign, d, lo, hi)) = &r.gold else {
//                 want_refusal(&what, got, false, "undecided_det")?;
//                 continue;
//             };
//             let v = match got {
//                 Ok(v) => v,
//                 Err(BoundError::Resource {
//                     detail: "undecided_det",
//                     ..
//                 }) => continue,
//                 Err(e) => return Err(format!("{what}: {e:?}")),
//             };
//             if v.sign() != *sign {
//                 return Err(format!("{what}: sign {}, gold {sign}", v.sign()));
//             }
//             want_holds(&what, v.det(), d)?;
//             if !v.ln_abs_det().is_informative() {
//                 return Err(format!(
//                     "{what}: ln_abs_det is all-reals (must stay finite, 16.LA1.7)"
//                 ));
//             }
//             let (covers, overlaps) = vs_interval(v.ln_abs_det(), lo, hi)?;
//             if covers {
//                 continue;
//             }
//             if !overlaps {
//                 return Err(format!(
//                     "{what}: ln_abs_det {:?} disjoint from [lo256, hi256]",
//                     v.ln_abs_det()
//                 ));
//             }
//             let mut decided = false;
//             for bits in [4096_u32, 16384] {
//                 if let Some((hlo, hhi)) = heavy_enclosure(r.id, bits)? {
//                     if vs_interval(v.ln_abs_det(), &hlo, &hhi)?.0 {
//                         decided = true;
//                     }
//                 }
//             }
//             if !decided {
//                 return Err(format!("{what}: partial overlap not re-decided by enc4096/enc16384 (set ZB_LA_HEAVY_DIR)"));
//             }
//         }
//     }
//     Ok(())
// }
//
// /// Gold format from the `# format: N` header (format 2 adds the six cond columns of scripts/la_cond_cols.py); none = 1.
// fn gold_format(text: &str) -> Result<u32, String> {
//     match text.lines().find_map(|l| l.strip_prefix("# format: ")) {
//         None => Ok(1),
//         Some(rest) => rest
//             .split_whitespace()
//             .next()
//             .ok_or("empty '# format:' header")?
//             .parse::<u32>()
//             .map_err(|e| format!("bad '# format:' header {rest:?}: {e}")),
//     }
// }
//
// /// The format-2 cond columns `cond_row k_row scale_row containment_only reason rref` of one non-refusal row (16.LA.5).
// struct CondCols {
//     k_row: usize,
//     scale_row: Rat,
//     containment_only: bool,
// }
//
// /// Parses the six cond columns at `toks[at..]` and cross-checks `containment_only` against `k_row + 8 >= max(PRECS)`.
// fn cond_cols(id: &str, toks: &[&str], at: usize) -> Result<CondCols, String> {
//     let rest = toks
//         .get(at..)
//         .ok_or_else(|| format!("{id}: row too short for cond columns"))?;
//     let [cond_row, k_row, scale_row, flag, reason, rref] = rest else {
//         return Err(format!("{id}: want 6 cond columns, got {}", rest.len()));
//     };
//     rat(cond_row)?;
//     rat(rref)?;
//     let k_row: usize = k_row
//         .parse()
//         .map_err(|e| format!("{id}: k_row {k_row:?}: {e}"))?;
//     let containment_only = match (*flag, *reason) {
//         ("0", "tight") => false,
//         ("1", "k_row_too_large") => true,
//         (f, r) => return Err(format!("{id}: bad containment_only/reason {f} {r}")),
//     };
//     let top = PRECS.iter().max().ok_or("PRECS")?;
//     if containment_only != (k_row.checked_add(8).ok_or("k_row overflows")? >= *top) {
//         return Err(format!(
//             "{id}: containment_only {flag} disagrees with k_row {k_row}"
//         ));
//     }
//     Ok(CondCols {
//         k_row,
//         scale_row: rat(scale_row)?,
//         containment_only,
//     })
// }
//
// /// 16.LA.5: `rad <= 2^(k_row - prec)·max(|v|, scale_row)` when `prec > k_row + 8` and the row is not containment-only;
// /// `|v|` is read as `|mid|`. Returns whether the bound was enforced.
// fn want_tight(what: &str, b: &Ball, c: &CondCols, prec: usize) -> Result<bool, String> {
//     if c.containment_only || prec <= c.k_row.checked_add(8).ok_or("k_row overflows")? {
//         return Ok(false);
//     }
//     if b.rad().is_infinity() {
//         return Err(format!("{what}: infinite radius under the k_row bound"));
//     }
//     let k = i64::try_from(c.k_row).map_err(|e| format!("{e}"))?;
//     let p = i64::try_from(prec).map_err(|e| format!("{e}"))?;
//     let factor = pow2(k.checked_sub(p).ok_or("k_row - prec")?)?;
//     let q = &c.scale_row.q;
//     let mut m = b.mid().mul_full_prec(q);
//     if m.is_negative() {
//         m = -m;
//     }
//     let mut s = plus_zero(&c.scale_row.p)?;
//     if s.is_negative() {
//         s = -s;
//     }
//     let big = if m >= s { m } else { s };
//     let lhs = b
//         .rad()
//         .to_exact_num()
//         .map_err(|e| format!("{e:?}"))?
//         .mul_full_prec(q);
//     if lhs <= factor.mul_full_prec(&big) {
//         Ok(true)
//     } else {
//         Err(format!(
//             "{what}: rad {:?} over the k_row {} bound at prec {prec}",
//             b.rad(),
//             c.k_row
//         ))
//     }
// }
//
// /// 16.LA.5 tightness on every format-2 logdet row (`ln_abs_det`, `scale_row = 1`) and solve row (`x`, `scale_row = ‖x‖_∞`).
// /// Red until the golds are format 2 (harper/zb-la-golds), S3 `ln` lands (16.LA1.6) and 16.S6.12 trisolve lands.
// #[test]
// #[ignore = "red: CONTRACT 16.LA.5: needs format-2 la_logdet.txt/la_solve.txt (cond_row k_row scale_row containment_only reason rref, harper/zb-la-golds), S3 ln (16.LA1.6) and 16.S6.12 trisolve; un-ignore when all three are on the branch"]
// fn la2_tightness_k_row_bound_on_logdet_and_solve_rows() -> Result<(), String> {
//     let ld = read_gold("cholesky", "la_logdet.txt")?;
//     let sv = read_gold("lu", "la_solve.txt")?;
//     for (name, text) in [("la_logdet.txt", &ld), ("la_solve.txt", &sv)] {
//         if gold_format(text)? != 2 {
//             return Err(format!(
//                 "{name} is format {}: no cond columns (16.LA.5)",
//                 gold_format(text)?
//             ));
//         }
//     }
//     let mut enforced = 0_usize;
//     for prec in PRECS {
//         for toks in &data_rows(&ld)? {
//             let r = logdet_row(toks)?;
//             if r.gold.is_none() {
//                 continue;
//             }
//             let at = mul_dims(r.n, r.n)?.checked_add(6).ok_or("width")?;
//             let c = cond_cols(r.id, toks, at)?;
//             if c.scale_row.p != num(1) || c.scale_row.q != num(1) {
//                 return Err(format!("{}: logdet scale_row must be 1/1", r.id));
//             }
//             let what = format!("{}@{prec} ln_abs_det", r.id);
//             let a = ball_mat(r.n, r.n, &r.a, prec)?;
//             let v = match a.logdet(prec) {
//                 Ok(v) => v,
//                 Err(BoundError::Resource {
//                     detail: "undecided_det",
//                     ..
//                 }) => continue,
//                 Err(e) => return Err(format!("{what}: {e:?}")),
//             };
//             if want_tight(&what, v.ln_abs_det(), &c, prec)? {
//                 enforced = enforced.checked_add(1).ok_or("count")?;
//             }
//         }
//         for toks in &data_rows(&sv)? {
//             let r = solve_row(toks)?;
//             if r.x.is_none() {
//                 continue;
//             }
//             let at = mul_dims(r.n, r.n)?
//                 .checked_add(r.n.checked_mul(2).ok_or("width")?)
//                 .and_then(|w| w.checked_add(2))
//                 .ok_or("width")?;
//             let c = cond_cols(r.id, toks, at)?;
//             let what = format!("{}@{prec} solve", r.id);
//             let a = ball_mat(r.n, r.n, &r.a, prec)?;
//             let x = match a.solve(&ball_vec(&r.b, prec)?, prec) {
//                 Ok(x) => x,
//                 Err(BoundError::Resource { detail, .. }) if detail.starts_with("undecided") => {
//                     continue
//                 }
//                 Err(e) => return Err(format!("{what}: {e:?}")),
//             };
//             for i in 0..r.n {
//                 if want_tight(&format!("{what} x[{i}]"), x.get(i).ok_or("x")?, &c, prec)? {
//                     enforced = enforced.checked_add(1).ok_or("count")?;
//                 }
//             }
//         }
//     }
//     if enforced == 0 {
//         return Err("no row enforced the k_row bound (16.LA.5)".to_string());
//     }
//     Ok(())
// }
