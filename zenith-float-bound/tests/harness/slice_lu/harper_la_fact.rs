// HOLD H09: held, not compiled (no `mod` line names this file). Index and use-or-delete rule: HOLD.md at the repo root.
// The original text follows, each line commented out with a leading `// `; uncomment it whole to go live.
// //! Harper wide-range factorization rows (2026-10-29) for CONTRACT 16.LA1 (Cholesky, LDLᵀ, logdet), 16.LA2 (LU), 16.LA3
// //! (QR, least squares) and the shared LA rules 16.LA.2 to 16.LA.5, written test-first against the stub API at 16384 bits:
// //! 32×32 SPD `C = L·D·Lᵀ` with `D` from 2^±203 out to near EXPONENT_MAX, exact `det = 2^k` and `logdet = k·ln 2` (ln 2 from
// //! zb's rigorous `constant_ball`), permutation LU, dyadic Hadamard QR, near-singular 2^-16000 pivots, singular and
// //! straddling inputs, `n = 0`, non-square input and the `"b_len"` rule.
// //! COVERS: BallCholesky::g, BallCholesky::into_g, BallCholesky::solve, BallCholesky::quad_form, BallCholesky::logdet_spd
// //! COVERS: BallLdl::l, BallLdl::d, BallLu::perm, BallLu::l, BallLu::u, BallLu::perm_sign, BallLu::det, BallLu::solve
// //! COVERS: BallQr::q, BallQr::r, BallQr::into_parts, BallQr::lstsq, LogDet::sign, LogDet::det, LogDet::ln_abs_det
// //! COVERS: BallMat::cholesky, BallMat::is_spd, BallMat::ldl, BallMat::logdet, BallMat::lstsq, BallMat::lu, BallMat::qr
// //! COVERS: BallMat::solve, BallMat::solve_mat
// //! COVERS: common::la_rows::hadamard, common::la_rows::i64_of, common::la_rows::sum_i64
//
// use zenith_float_bound::array::{BallMat, BallVec};
// use zenith_float_bound::ball::Ball;
// use zenith_float_bound::enclosure::Enclosure;
// use zenith_float_bound::error::BoundError;
// use zenith_float_bound::gls::gls;
// use zenith_float_bound::mag::Mag;
// use zenith_float_bound::tri::Tri;
// use zenith_float_num::{ExactNum, EXPONENT_MAX};
//
// use crate::common::generators::TOP;
// use crate::common::la_rows::{
//     at, diag, dy, fuzzy, hadamard, i64_of, int, ldlt, mat, mat_with, mids, mul_exact, ok, pow2,
//     pows, sum_i64, transpose, unit_lower, vat, vecb, want_all_reals, want_contains, want_encloses,
//     want_err, want_k_ln2, want_mat, want_vec, P,
// };
//
// const N: usize = 32;
//
// /// `x_i = (i mod 5 − 2)·2^(i mod 9)`.
// fn xs(n: usize) -> Result<Vec<ExactNum>, String> {
//     let mut v = Vec::new();
//     for i in 0..n {
//         let k = i64_of(i)?;
//         let m = k
//             .checked_rem(5)
//             .and_then(|r| r.checked_sub(2))
//             .ok_or("xs")?;
//         v.push(dy(m, k.checked_rem(9).ok_or("xs")?)?);
//     }
//     Ok(v)
// }
//
// /// Exponents `±e` alternating, `n` of them.
// fn alt(n: usize, e: i64) -> Result<Vec<i64>, String> {
//     let ne = e.checked_neg().ok_or("alt")?;
//     let mut v = Vec::new();
//     for i in 0..n {
//         v.push(if i.checked_rem(2) == Some(0) { e } else { ne });
//     }
//     Ok(v)
// }
//
// /// Halves of even exponents.
// fn halves(ks: &[i64]) -> Result<Vec<i64>, String> {
//     let mut v = Vec::new();
//     for k in ks {
//         if k.checked_rem(2) != Some(0) {
//             return Err(format!("odd exponent {k}"));
//         }
//         v.push(k.checked_div(2).ok_or("half")?);
//     }
//     Ok(v)
// }
//
// /// The sign of a permutation (`perm[i]` is the row of `A` that is row `i` of `P·A`) by inversion count.
// fn sign_of(perm: &[usize]) -> i8 {
//     let mut odd = false;
//     for (i, a) in perm.iter().enumerate() {
//         for b in perm.iter().skip(i).skip(1) {
//             if a > b {
//                 odd = !odd;
//             }
//         }
//     }
//     if odd {
//         -1
//     } else {
//         1
//     }
// }
//
// /// `P·A` for the row-major `n×n` grid `a`.
// fn permute(n: usize, perm: &[usize], a: &[ExactNum]) -> Result<Vec<ExactNum>, String> {
//     let mut g = Vec::new();
//     for p in perm {
//         for j in 0..n {
//             g.push(int(0).add_full_prec(crate::common::la_rows::cell(a, n, *p, j)?));
//         }
//     }
//     Ok(g)
// }
//
// /// The 32×32 SPD row: `L = unit_lower(32, 5)` and `D = 2^ks`.
// fn spd(ks: &[i64]) -> Result<(Vec<ExactNum>, Vec<ExactNum>), String> {
//     let l = unit_lower(N, 5)?;
//     let c = ldlt(N, &l, &pows(ks)?)?;
//     Ok((l, c))
// }
//
// fn is_dm(e: &BoundError) -> bool {
//     matches!(e, BoundError::DimensionMismatch { .. })
// }
// fn is_rect(e: &BoundError) -> bool {
//     matches!(
//         e,
//         BoundError::DimensionMismatch {
//             detail: "rectangular",
//             ..
//         }
//     )
// }
// fn is_blen(e: &BoundError) -> bool {
//     matches!(
//         e,
//         BoundError::DimensionMismatch {
//             detail: "b_len",
//             ..
//         }
//     )
// }
// fn is_npd(e: &BoundError) -> bool {
//     matches!(e, BoundError::NotPositiveDefinite { detail: "a", .. })
// }
// fn is_nonsym(e: &BoundError) -> bool {
//     matches!(e, BoundError::NotSymmetric { .. })
// }
// fn is_sing(e: &BoundError) -> bool {
//     matches!(
//         e,
//         BoundError::Singular {
//             detail: "singular",
//             ..
//         }
//     )
// }
// fn is_prec(e: &BoundError) -> bool {
//     matches!(e, BoundError::InvalidArgument { .. })
// }
// fn is_upiv(e: &BoundError) -> bool {
//     matches!(
//         e,
//         BoundError::Resource {
//             detail: "undecided_pivot",
//             ..
//         }
//     )
// }
// fn is_udet(e: &BoundError) -> bool {
//     matches!(
//         e,
//         BoundError::Resource {
//             detail: "undecided_det",
//             ..
//         }
//     )
// }
// fn is_rank(e: &BoundError) -> bool {
//     matches!(
//         e,
//         BoundError::Singular {
//             detail: "rank_deficient",
//             ..
//         }
//     )
// }
// fn is_urank(e: &BoundError) -> bool {
//     matches!(
//         e,
//         BoundError::Resource {
//             detail: "undecided_rank",
//             ..
//         }
//     )
// }
// fn is_wide(e: &BoundError) -> bool {
//     matches!(e, BoundError::DimensionMismatch { detail: "wide", .. })
// }
//
// // ---- 16.LA1 -----------------------------------------------------------------------------------------------------------
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA1.1: cholesky/ldl/is_spd API on a 32x32 SPD L·D·Lᵀ with D = 2^±202, 2^±203"]
// fn la1_cholesky_ldl_api_32_spd_pow2() -> Result<(), String> {
//     let ks = alt(N, 202)?;
//     let (l, c) = spd(&ks)?;
//     let a = mat(N, N, c)?;
//     let g = mul_exact(N, N, N, &l, &diag(&pows(&halves(&ks)?)?))?;
//     let ch = ok("cholesky", a.cholesky(P))?;
//     want_mat("G", ch.g(), N, N, &g)?;
//     let x = xs(N)?;
//     let r = mul_exact(N, N, 1, &mids(&a)?, &x)?;
//     want_vec(
//         "BallCholesky::solve",
//         &ok(
//             "solve",
//             ch.solve(&vecb(mul_exact(N, N, 1, &mids(&a)?, &x)?)?, P),
//         )?,
//         &x,
//     )?;
//     let q = mul_exact(1, N, 1, &x, &r)?;
//     want_encloses(
//         "quad_form rᵀC⁻¹r = xᵀCx",
//         &ok("quad_form", ch.quad_form(&vecb(r)?, P))?,
//         q.first().ok_or("q")?,
//     )?;
//     want_mat("into_g", &ch.into_g(), N, N, &g)?;
//     let ks3 = alt(N, 203)?;
//     let (l3, c3) = spd(&ks3)?;
//     let ld = ok("ldl", mat(N, N, c3)?.ldl(P))?;
//     want_mat("ldl L", ld.l(), N, N, &l3)?;
//     want_vec("ldl d", ld.d(), &pows(&ks3)?)?;
//     match ok("is_spd", a.is_spd(P))? {
//         Tri::True => Ok(()),
//         t @ (Tri::False | Tri::Unknown) => Err(format!("is_spd: {t:?}, want True")),
//     }
// }
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA1.2: lower triangle read, pivot rules NotPositiveDefinite \"a\" / undecided_pivot, is_spd"]
// fn la1_semantics_lower_triangle_and_pivot_rules() -> Result<(), String> {
//     // the upper triangle is all-reals (overlaps every value, so not NotSymmetric) and is never read by the factorization
//     let ks = alt(N, 202)?;
//     let (l, c) = spd(&ks)?;
//     let ar = || ok("all_reals", Ball::all_reals(P));
//     let a = mat_with(
//         N,
//         N,
//         c,
//         vec![(0, 31, ar()?), (7, 9, ar()?), (30, 31, ar()?)],
//     )?;
//     let g = mul_exact(N, N, N, &l, &diag(&pows(&halves(&ks)?)?))?;
//     want_mat(
//         "G with all-reals upper",
//         ok("cholesky", a.cholesky(P))?.g(),
//         N,
//         N,
//         &g,
//     )?;
//     let neg = mat(2, 2, vec![int(4), int(2), int(2), int(-3)])?;
//     want_err(
//         "pivot 2 = -4",
//         neg.cholesky(P),
//         is_npd,
//         "NotPositiveDefinite a",
//     )?;
//     want_err(
//         "ldl pivot 2 = -4",
//         neg.ldl(P),
//         is_npd,
//         "NotPositiveDefinite a",
//     )?;
//     let straddle = mat_with(
//         2,
//         2,
//         vec![int(4), int(2), int(2), int(1)],
//         vec![(1, 1, fuzzy(int(1), 2)?)],
//     )?;
//     want_err(
//         "pivot 0 ± 2",
//         straddle.cholesky(P),
//         is_upiv,
//         "Resource undecided_pivot",
//     )?;
//     let want_tri = |what: &str, got: Tri, want: Tri| {
//         if got == want {
//             Ok(())
//         } else {
//             Err(format!("is_spd {what}: {got:?} want {want:?}"))
//         }
//     };
//     want_tri("negative", ok("is_spd", neg.is_spd(P))?, Tri::False)?;
//     want_tri(
//         "straddling",
//         ok("is_spd", straddle.is_spd(P))?,
//         Tri::Unknown,
//     )?;
//     let nonsym = mat(2, 2, vec![int(2), int(1), int(3), int(2)])?;
//     want_tri("non-symmetric", ok("is_spd", nonsym.is_spd(P))?, Tri::False)
// }
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA1.3: factors of a fuzzy 32x32 SPD ball enclose the midpoint's exact factors"]
// fn la1_enclosure_fuzzy_spd() -> Result<(), String> {
//     let ks = alt(N, 4)?;
//     let (l, c) = spd(&ks)?;
//     let mut repl = Vec::new();
//     for i in 0..N {
//         for j in 0..N {
//             repl.push((
//                 i,
//                 j,
//                 fuzzy(
//                     int(0).add_full_prec(crate::common::la_rows::cell(&c, N, i, j)?),
//                     -60,
//                 )?,
//             ));
//         }
//     }
//     let a = mat_with(N, N, mids(&mat(N, N, c)?)?, repl)?;
//     let g = mul_exact(N, N, N, &l, &diag(&pows(&halves(&ks)?)?))?;
//     want_mat("fuzzy G", ok("cholesky", a.cholesky(P))?.g(), N, N, &g)?;
//     let ld = ok("ldl", a.ldl(P))?;
//     want_mat("fuzzy ldl L", ld.l(), N, N, &l)?;
//     want_vec("fuzzy ldl d", ld.d(), &pows(&ks)?)
// }
//
// /// The 32×32 row-permuted diagonal: row `i` holds `s_i·2^(k_i)` in column `(i + 5) mod 32`.
// fn perm_diag(ks: &[i64]) -> Result<(Vec<ExactNum>, Vec<usize>, i8, i64), String> {
//     let mut g = Vec::new();
//     let mut neg = 0usize;
//     for i in 0..N {
//         let col = i
//             .checked_add(5)
//             .and_then(|x| x.checked_rem(N))
//             .ok_or("col")?;
//         let k = *ks.get(i).ok_or("ks")?;
//         let s = if i.checked_rem(3) == Some(0) { -1 } else { 1 };
//         if s < 0 {
//             neg = neg.checked_add(1).ok_or("neg")?;
//         }
//         for j in 0..N {
//             g.push(if j == col { dy(s, k)? } else { int(0) });
//         }
//     }
//     let mut perm = Vec::new();
//     for k in 0..N {
//         perm.push(
//             k.checked_add(N)
//                 .and_then(|x| x.checked_sub(5))
//                 .and_then(|x| x.checked_rem(N))
//                 .ok_or("perm")?,
//         );
//     }
//     let ps = sign_of(&perm);
//     let sign = if neg.checked_rem(2) == Some(0) {
//         ps
//     } else {
//         ps.checked_neg().ok_or("sign")?
//     };
//     Ok((g, perm, sign, sum_i64(ks)?))
// }
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA1.4: logdet of a 32x32 permuted 2^k diagonal: sign, det = ±2^Σk, ln|det| = Σk·ln2; refusals"]
// fn la1_general_logdet_permuted_pow2() -> Result<(), String> {
//     let mut ks = Vec::new();
//     for i in 0..N {
//         ks.push(
//             i64_of(i)?
//                 .checked_mul(61)
//                 .and_then(|x| x.checked_sub(900))
//                 .ok_or("ks")?,
//         );
//     }
//     let (g, _, sign, k) = perm_diag(&ks)?;
//     let ld = ok("logdet", mat(N, N, g)?.logdet(P))?;
//     if ld.sign() != sign {
//         return Err(format!("logdet sign {} want {sign}", ld.sign()));
//     }
//     want_encloses("det", ld.det(), &dy(i64::from(sign), k)?)?;
//     want_k_ln2("ln|det|", ld.ln_abs_det(), k)?;
//     want_err(
//         "zero column",
//         mat(2, 2, vec![int(0), int(1), int(0), int(2)])?.logdet(P),
//         is_sing,
//         "Singular singular",
//     )?;
//     let s = BallMat::from_balls(1, 1, vec![fuzzy(int(1), 2)?]);
//     want_err(
//         "straddling u00",
//         ok("from_balls", s)?.logdet(P),
//         is_udet,
//         "Resource undecided_det",
//     )
// }
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA1.5: rectangular, prec-first order, n = 0 empty factors and logdet sign 1, det 1, ln 0"]
// fn la1_errors_and_n0() -> Result<(), String> {
//     let r = mat(2, 3, vec![int(1), int(0), int(0), int(0), int(1), int(0)])?;
//     want_err(
//         "cholesky rect",
//         r.cholesky(P),
//         is_rect,
//         "DimensionMismatch rectangular",
//     )?;
//     want_err(
//         "ldl rect",
//         r.ldl(P),
//         is_rect,
//         "DimensionMismatch rectangular",
//     )?;
//     want_err(
//         "logdet rect",
//         r.logdet(P),
//         is_rect,
//         "DimensionMismatch rectangular",
//     )?;
//     want_err("prec 0 first", r.cholesky(0), is_prec, "InvalidArgument")?;
//     let e = mat(0, 0, Vec::new())?;
//     let ch = ok("cholesky n=0", e.cholesky(P))?;
//     let ld = ok("ldl n=0", e.ldl(P))?;
//     if ch.g().nrows() != 0 || ld.l().nrows() != 0 || !ld.d().is_empty() {
//         return Err("n = 0 factors are not empty".to_string());
//     }
//     let lg = ok("logdet n=0", e.logdet(P))?;
//     if lg.sign() != 1 {
//         return Err(format!("n = 0 sign {}", lg.sign()));
//     }
//     want_encloses("det n=0", lg.det(), &int(1))?;
//     want_encloses("ln n=0", lg.ln_abs_det(), &int(0))?;
//     if !lg.det().rad().is_zero() || !lg.ln_abs_det().rad().is_zero() {
//         return Err(format!("n = 0 want 1 ± 0 and 0 ± 0, got {lg:?}"));
//     }
//     Ok(())
// }
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA1.6: post-S3 logdet/logdet_spd values, det beyond EXPONENT_MAX gives all-reals det, finite ln"]
// fn la1_logdet_near_exponent_max() -> Result<(), String> {
//     let big = EXPONENT_MAX
//         .checked_sub(40)
//         .and_then(|x| x.checked_div(2))
//         .and_then(|x| x.checked_mul(2))
//         .ok_or("big")?;
//     let mut ks = alt(N, 2)?;
//     for k in ks.iter_mut().take(2) {
//         *k = big;
//     }
//     let (_, c) = spd(&ks)?;
//     let k = sum_i64(&ks)?;
//     let a = mat(N, N, c)?;
//     let ld = ok("logdet", a.logdet(P))?;
//     if ld.sign() != 1 {
//         return Err(format!("sign {}", ld.sign()));
//     }
//     want_all_reals("det = 2^(2·(EMAX−40)) overflows", ld.det())?;
//     want_k_ln2("ln_abs_det", ld.ln_abs_det(), k)?;
//     want_k_ln2(
//         "logdet_spd",
//         &ok("logdet_spd", ok("cholesky", a.cholesky(P))?.logdet_spd(P))?,
//         k,
//     )?;
//     let ks2 = alt(N, 202)?;
//     let (_, c2) = spd(&ks2)?;
//     let ch = ok("cholesky", mat(N, N, c2)?.cholesky(P))?;
//     want_k_ln2(
//         "logdet_spd 2^±202",
//         &ok("logdet_spd", ch.logdet_spd(P))?,
//         sum_i64(&ks2)?,
//     )
// }
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA1.7: logdet at precs 53, 113, 192 holds k·ln2 under the 16.LA.5 re-decide"]
// fn la1_logdet_at_replay_precisions() -> Result<(), String> {
//     let mut ks = Vec::new();
//     for i in 0..N {
//         ks.push(
//             i64_of(i)?
//                 .checked_mul(2)
//                 .and_then(|x| x.checked_sub(17))
//                 .ok_or("ks")?,
//         );
//     }
//     let (_, c) = spd(&ks)?;
//     let a = mat(N, N, c)?;
//     for prec in [53, 113, 192] {
//         let ld = ok("logdet", a.logdet(prec))?;
//         want_k_ln2(
//             &format!("ln|det| at {prec}"),
//             ld.ln_abs_det(),
//             sum_i64(&ks)?,
//         )?;
//     }
//     Ok(())
// }
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA1.8: mutant killers 1000 1003 1005 1007 (straddle, mid symmetry, radii, d_j)"]
// fn la1_mutant_killers() -> Result<(), String> {
//     // 1000: mid > 0 but lower bound < 0
//     let s = ok(
//         "from_balls",
//         BallMat::from_balls(1, 1, vec![fuzzy(int(1), 2)?]),
//     )?;
//     want_err(
//         "1000 straddling pivot",
//         s.cholesky(P),
//         is_upiv,
//         "Resource undecided_pivot",
//     )?;
//     // 1003: a01 = 1 ± 1 and a10 = 3/2 overlap: symmetric matrices exist in the ball, so no NotSymmetric
//     let o = mat_with(
//         2,
//         2,
//         vec![int(4), int(1), dy(3, -1)?, int(4)],
//         vec![(0, 1, fuzzy(int(1), 1)?)],
//     )?;
//     let g = ok("1003 overlapping pair", o.cholesky(P))?;
//     want_contains("1003 g00", at(g.g(), 0, 0)?, &int(2))?;
//     // 1005: a00 = 4 ± 1/2 so g00 ⊇ [sqrt(3.5), sqrt(4.5)] ∋ 15/8 and 17/8
//     let f = ok(
//         "from_balls",
//         BallMat::from_balls(1, 1, vec![fuzzy(int(4), 0)?]),
//     )?;
//     let gf = ok("1005", f.cholesky(P))?;
//     want_contains("1005 lower", at(gf.g(), 0, 0)?, &dy(15, -3)?)?;
//     want_contains("1005 upper", at(gf.g(), 0, 0)?, &dy(17, -3)?)?;
//     // 1007: d_1 = 2^-203 enters l_21
//     let ld = ok(
//         "1007",
//         mat(3, 3, ldlt(3, &unit_lower(3, 1)?, &pows(&[0, -203, 203])?)?)?.ldl(P),
//     )?;
//     want_mat("1007 L", ld.l(), 3, 3, &unit_lower(3, 1)?)
// }
//
// // ---- 16.LA2 -----------------------------------------------------------------------------------------------------------
//
// fn lu3() -> Result<Vec<ExactNum>, String> {
//     Ok(vec![
//         int(1),
//         int(2),
//         int(0),
//         int(4),
//         int(0),
//         int(1),
//         int(2),
//         int(1),
//         int(3),
//     ])
// }
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA2.3: lu perm/L/U on 3x3 exact, lowest-index tie, 32x32 permuted diagonal, zero column"]
// fn la2_lu_signatures_and_algorithm() -> Result<(), String> {
//     let lu = ok("lu", mat(3, 3, lu3()?)?.lu(P))?;
//     if lu.perm() != [1, 0, 2] || lu.perm_sign() != -1 {
//         return Err(format!("perm {:?} sign {}", lu.perm(), lu.perm_sign()));
//     }
//     let l = vec![
//         int(1),
//         int(0),
//         int(0),
//         dy(1, -2)?,
//         int(1),
//         int(0),
//         dy(1, -1)?,
//         dy(1, -1)?,
//         int(1),
//     ];
//     let u = vec![
//         int(4),
//         int(0),
//         int(1),
//         int(0),
//         int(2),
//         dy(-1, -2)?,
//         int(0),
//         int(0),
//         dy(21, -3)?,
//     ];
//     want_mat("L", lu.l(), 3, 3, &l)?;
//     want_mat("U", lu.u(), 3, 3, &u)?;
//     want_encloses("det", &ok("det", lu.det(P))?, &int(-21))?;
//     let tie = ok(
//         "lu tie",
//         mat(2, 2, vec![int(1), int(0), int(-1), int(1)])?.lu(P),
//     )?;
//     if tie.perm() != [0, 1] {
//         return Err(format!("tie perm {:?}", tie.perm()));
//     }
//     let mut ks = Vec::new();
//     for i in 0..N {
//         ks.push(
//             i64_of(i)?
//                 .checked_mul(13)
//                 .and_then(|x| x.checked_sub(203))
//                 .ok_or("ks")?,
//         );
//     }
//     let (g, perm, sign, k) = perm_diag(&ks)?;
//     let big = ok("lu 32", mat(N, N, ldlt_free(&g))?.lu(P))?;
//     if big.perm() != perm.as_slice() {
//         return Err(format!("32x32 perm {:?} want {perm:?}", big.perm()));
//     }
//     let pa = permute(N, &perm, &g)?;
//     want_mat("32x32 U = P·A", big.u(), N, N, &pa)?;
//     want_mat("32x32 L = I", big.l(), N, N, &diag(&pows(&[0; N])?))?;
//     want_encloses(
//         "32x32 det",
//         &ok("det", big.det(P))?,
//         &dy(i64::from(sign), k)?,
//     )?;
//     let z = ok(
//         "lu zero column",
//         mat(2, 2, vec![int(0), int(1), int(0), int(2)])?.lu(P),
//     )?;
//     if z.perm() != [0, 1] {
//         return Err(format!("zero column swapped: {:?}", z.perm()));
//     }
//     let u00 = at(z.u(), 0, 0)?;
//     if !u00.mid().is_zero() || !u00.rad().is_zero() {
//         return Err(format!("zero column u00 {u00:?} want 0 ± 0"));
//     }
//     let d = ok("det", z.det(P))?;
//     if !d.mid().is_zero() || !d.rad().is_zero() {
//         return Err(format!("zero column det {d:?} want 0 ± 0"));
//     }
//     Ok(())
// }
//
// /// The grid itself (named for the reader: no LDLᵀ structure).
// fn ldlt_free(g: &[ExactNum]) -> Vec<ExactNum> {
//     let mut v = Vec::new();
//     for x in g {
//         v.push(int(0).add_full_prec(x));
//     }
//     v
// }
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA2.4: solve pre-check: exact 0 ± 0 u_kk Singular \"singular\", straddling Resource undecided_pivot"]
// fn la2_solve_refusals() -> Result<(), String> {
//     let b = vecb(vec![int(1), int(1)])?;
//     let z = mat(2, 2, vec![int(0), int(1), int(0), int(2)])?;
//     want_err(
//         "BallMat::solve zero column",
//         z.solve(&b, P),
//         is_sing,
//         "Singular singular",
//     )?;
//     want_err(
//         "BallLu::solve zero column",
//         ok("lu", z.lu(P))?.solve(&b, P),
//         is_sing,
//         "Singular singular",
//     )?;
//     let s = mat_with(
//         2,
//         2,
//         vec![int(0), int(0), int(0), int(1)],
//         vec![(0, 0, fuzzy(int(1), 2)?)],
//     )?;
//     want_err(
//         "straddling u00",
//         s.solve(&b, P),
//         is_upiv,
//         "Resource undecided_pivot",
//     )?;
//     want_err(
//         "straddling u00 (BallLu)",
//         ok("lu", s.lu(P))?.solve(&b, P),
//         is_upiv,
//         "Resource undecided_pivot",
//     )
// }
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA2.5: solve/solve_mat/det of a dense 32x32 L·U with U diagonal 2^-203..2^200 enclose exact"]
// fn la2_enclosure_dense_32() -> Result<(), String> {
//     let mut ks = Vec::new();
//     for i in 0..N {
//         ks.push(
//             i64_of(i)?
//                 .checked_mul(13)
//                 .and_then(|x| x.checked_sub(203))
//                 .ok_or("ks")?,
//         );
//     }
//     let u = mul_exact(
//         N,
//         N,
//         N,
//         &diag(&pows(&ks)?),
//         &transpose(N, N, &unit_lower(N, 7)?)?,
//     )?;
//     let a = mul_exact(N, N, N, &unit_lower(N, 6)?, &u)?;
//     let am = mat(N, N, a)?;
//     let x = xs(N)?;
//     let b = mul_exact(N, N, 1, &mids(&am)?, &x)?;
//     want_vec("solve", &ok("solve", am.solve(&vecb(b)?, P))?, &x)?;
//     let mut x2 = Vec::new();
//     for v in &x {
//         x2.push(int(0).add_full_prec(v));
//         x2.push(v.mul_full_prec(&pow2(-150)?));
//     }
//     let b2 = mul_exact(N, N, 2, &mids(&am)?, &x2)?;
//     want_mat(
//         "solve_mat",
//         &ok("solve_mat", am.solve_mat(&mat(N, 2, b2)?, P))?,
//         N,
//         2,
//         &x2,
//     )?;
//     want_encloses(
//         "det = 2^Σk",
//         &ok("det", ok("lu", am.lu(P))?.det(P))?,
//         &pow2(sum_i64(&ks)?)?,
//     )
// }
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA2.6: P_zb·A lies in L·U entrywise; rectangular lu is DimensionMismatch \"rectangular\""]
// fn la2_replay_rule_pa_in_lu() -> Result<(), String> {
//     for (n, g) in [
//         (3, lu3()?),
//         (
//             3,
//             vec![
//                 int(1),
//                 int(2),
//                 int(3),
//                 int(-1),
//                 int(5),
//                 int(1),
//                 int(1),
//                 int(0),
//                 int(7),
//             ],
//         ),
//     ] {
//         let a = mat(n, n, ldlt_free(&g))?;
//         let lu = ok("lu", a.lu(P))?;
//         let pa = permute(n, lu.perm(), &g)?;
//         let prod = ok("L·U", lu.l().matmul(lu.u(), P))?;
//         want_mat("P·A ∈ L·U", &prod, n, n, &pa)?;
//     }
//     let r = mat(2, 3, vec![int(1), int(0), int(0), int(0), int(1), int(0)])?;
//     want_err("lu rect", r.lu(P), is_rect, "DimensionMismatch rectangular")?;
//     want_err(
//         "logdet rect",
//         r.logdet(P),
//         is_rect,
//         "DimensionMismatch rectangular",
//     )
// }
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA2.7: mutant killers 1105 1106 1109 1110 (perm sign, zero det, P applied to b)"]
// fn la2_mutant_killers() -> Result<(), String> {
//     let a = mat(3, 3, lu3()?)?;
//     // 1105: det must carry perm_sign: −21, not 21
//     let d = ok("det", ok("lu", a.lu(P))?.det(P))?;
//     if d.contains(&int(21)) == Tri::True {
//         return Err(format!("1105: det {d:?} holds +21"));
//     }
//     // 1110: x = (1, 1, 1) needs b permuted: b = A·1 = (3, 5, 6)
//     want_vec(
//         "1110 solve",
//         &ok("solve", a.solve(&vecb(vec![int(3), int(5), int(6)])?, P))?,
//         &[int(1), int(1), int(1)],
//     )?;
//     // 1109: the exact zero column refuses as Singular, not Resource
//     want_err(
//         "1109",
//         mat(1, 1, vec![int(0)])?.solve(&vecb(vec![int(1)])?, P),
//         is_sing,
//         "Singular singular",
//     )
// }
//
// // ---- 16.LA3 -----------------------------------------------------------------------------------------------------------
//
// const M: usize = 16;
// const Q: usize = 8;
//
// /// Row-major `(Q, R, A)` entries.
// type QrRow = (Vec<ExactNum>, Vec<ExactNum>, Vec<ExactNum>);
//
// /// `A = Q·R`: `Q` the first 8 orthonormal Hadamard columns of order 16, `R = diag(2^rs)·unit_upper`.
// fn qr_row(rs: &[i64]) -> Result<QrRow, String> {
//     let q = hadamard(M, &[0, 1, 2, 3, 4, 5, 6, 7])?;
//     let r = mul_exact(
//         Q,
//         Q,
//         Q,
//         &diag(&pows(rs)?),
//         &transpose(Q, Q, &unit_lower(Q, 9)?)?,
//     )?;
//     let a = mul_exact(M, Q, Q, &q, &r)?;
//     Ok((q, r, a))
// }
//
// const RS: [i64; 8] = [203, -203, 0, 5, -5, 100, -100, 1];
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA3.1: qr q/r/into_parts and both lstsq on a 16x8 dyadic Hadamard·R, R diagonal 2^±203"]
// fn la3_qr_signatures_hadamard() -> Result<(), String> {
//     let (q, r, a) = qr_row(&RS)?;
//     let am = mat(M, Q, a)?;
//     let f = ok("qr", am.qr(P))?;
//     want_mat("Q", f.q(), M, Q, &q)?;
//     want_mat("R", f.r(), Q, Q, &r)?;
//     let x = xs(Q)?;
//     let mut b = mul_exact(M, Q, 1, &mids(&am)?, &x)?;
//     // add 3·(Hadamard column 8), orthogonal to range(A): the least-squares solution is still x
//     let h8 = hadamard(M, &[8])?;
//     for (bi, hi) in b.iter_mut().zip(h8.iter()) {
//         *bi = bi.add_full_prec(&hi.mul_full_prec(&int(3)));
//     }
//     let bv = vecb(b)?;
//     want_vec("BallQr::lstsq", &ok("lstsq", f.lstsq(&bv, P))?, &x)?;
//     want_vec("BallMat::lstsq", &ok("lstsq", am.lstsq(&bv, P))?, &x)?;
//     let (q2, r2) = f.into_parts();
//     want_mat("into_parts Q", &q2, M, Q, &q)?;
//     want_mat("into_parts R", &r2, Q, Q, &r)
// }
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA3.2: diag(R) > 0 sign fix: qr(−A) gives −Q and the same R"]
// fn la3_algorithm_sign_fix() -> Result<(), String> {
//     let (q, r, a) = qr_row(&RS)?;
//     let mut na = Vec::new();
//     for v in &a {
//         na.push(v.mul_full_prec(&int(-1)));
//     }
//     let mut nq = Vec::new();
//     for v in &q {
//         nq.push(v.mul_full_prec(&int(-1)));
//     }
//     let f = ok("qr", mat(M, Q, na)?.qr(P))?;
//     want_mat("−Q", f.q(), M, Q, &nq)?;
//     want_mat("R", f.r(), Q, Q, &r)
// }
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA3.3: wide, b_len, exact zero column rank_deficient, straddling undecided_rank, m×0 empty Ok"]
// fn la3_refusals_and_empty() -> Result<(), String> {
//     let w = mat(2, 3, vec![int(1), int(0), int(0), int(0), int(1), int(0)])?;
//     want_err("qr wide", w.qr(P), is_wide, "DimensionMismatch wide")?;
//     want_err(
//         "lstsq wide",
//         w.lstsq(&vecb(vec![int(1), int(1)])?, P),
//         is_wide,
//         "DimensionMismatch wide",
//     )?;
//     let (_, _, a) = qr_row(&RS)?;
//     let am = mat(M, Q, a)?;
//     want_err(
//         "lstsq b_len",
//         am.lstsq(&vecb(xs(Q)?)?, P),
//         is_blen,
//         "DimensionMismatch b_len",
//     )?;
//     let z = mat(3, 2, vec![int(0), int(1), int(0), int(1), int(0), int(0)])?;
//     want_err("zero column", z.qr(P), is_rank, "Singular rank_deficient")?;
//     let s = mat_with(
//         3,
//         1,
//         vec![int(0), int(0), int(0)],
//         vec![(0, 0, fuzzy(int(0), 0)?)],
//     )?;
//     want_err("straddling σ", s.qr(P), is_urank, "Resource undecided_rank")?;
//     let e = ok("from_balls 3x0", BallMat::from_balls(3, 0, Vec::new()))?;
//     let x = ok(
//         "lstsq 3x0",
//         e.lstsq(&vecb(vec![int(1), int(2), int(3)])?, P),
//     )?;
//     if !x.is_empty() {
//         return Err(format!("3x0 lstsq length {}", x.len()));
//     }
//     want_err(
//         "3x0 b_len first",
//         e.lstsq(&vecb(vec![int(1)])?, P),
//         is_blen,
//         "DimensionMismatch b_len",
//     )?;
//     let x0 = ok(
//         "lstsq 0x0",
//         mat(0, 0, Vec::new())?.lstsq(&vecb(Vec::new())?, P),
//     )?;
//     if !x0.is_empty() {
//         return Err(format!("0x0 lstsq length {}", x0.len()));
//     }
//     Ok(())
// }
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA3.4: Q, R, lstsq of a fuzzy (2^-61) Hadamard·R ball enclose the midpoint's factors"]
// fn la3_enclosure_fuzzy() -> Result<(), String> {
//     let rs = [0, 1, 2, 3, -1, -2, -3, 0];
//     let (q, r, a) = qr_row(&rs)?;
//     let mut repl = Vec::new();
//     for i in 0..M {
//         for j in 0..Q {
//             repl.push((
//                 i,
//                 j,
//                 fuzzy(
//                     int(0).add_full_prec(crate::common::la_rows::cell(&a, Q, i, j)?),
//                     -60,
//                 )?,
//             ));
//         }
//     }
//     let am = mat_with(M, Q, ldlt_free(&a), repl)?;
//     let f = ok("qr", am.qr(P))?;
//     want_mat("fuzzy Q", f.q(), M, Q, &q)?;
//     want_mat("fuzzy R", f.r(), Q, Q, &r)?;
//     let x = xs(Q)?;
//     let b = mul_exact(M, Q, 1, &a, &x)?;
//     want_vec("fuzzy lstsq", &ok("lstsq", am.lstsq(&vecb(b)?, P))?, &x)
// }
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA3.5: mutant killers 1207 1209 (thin Q shape, Qᵀb) on a 4x2 dyadic row"]
// fn la3_mutant_killers() -> Result<(), String> {
//     let q = hadamard(4, &[0, 1])?;
//     let r = vec![int(2), int(1), int(0), dy(1, -1)?];
//     let a = mul_exact(4, 2, 2, &q, &r)?;
//     let am = mat(4, 2, ldlt_free(&a))?;
//     let f = ok("qr", am.qr(P))?;
//     if f.q().nrows() != 4 || f.q().ncols() != 2 {
//         return Err(format!("1207: Q is {}x{}", f.q().nrows(), f.q().ncols()));
//     }
//     // b = A·(1, −2): x = R⁻¹Qᵀb
//     let x = vec![int(1), int(-2)];
//     want_vec(
//         "1209",
//         &ok("lstsq", am.lstsq(&vecb(mul_exact(4, 2, 1, &a, &x)?)?, P))?,
//         &x,
//     )
// }
//
// // ---- 16.LA.2 to 16.LA.5 ------------------------------------------------------------------------------------------------
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA.2: refusals only when certified: 2^-16000 pivots both ways, exact psd, overlapping pair"]
// fn la_certified_refusals_near_singular() -> Result<(), String> {
//     let e = pow2(-16000)?;
//     let up = mat(2, 2, vec![int(1), int(1), int(1), int(1).add_full_prec(&e)])?;
//     let g = ok("cholesky 1 + 2^-16000", up.cholesky(P))?;
//     want_encloses("g11 = 2^-8000", at(g.g(), 1, 1)?, &pow2(-8000)?)?;
//     want_encloses("ldl d1 = 2^-16000", vat(ok("ldl", up.ldl(P))?.d(), 1)?, &e)?;
//     let down = mat(2, 2, vec![int(1), int(1), int(1), int(1).sub_full_prec(&e)])?;
//     want_err(
//         "pivot −2^-16000",
//         down.cholesky(P),
//         is_npd,
//         "NotPositiveDefinite a",
//     )?;
//     want_err(
//         "exact psd [[1,1],[1,1]]",
//         mat(2, 2, vec![int(1), int(1), int(1), int(1)])?.cholesky(P),
//         is_npd,
//         "NotPositiveDefinite a",
//     )?;
//     let lg = ok("logdet near-singular", up.logdet(P))?;
//     want_encloses("det 2^-16000", lg.det(), &e)?;
//     want_k_ln2("ln det", lg.ln_abs_det(), -16000)?;
//     // an all-reals entry: undecided (Resource), never Singular
//     let ar = mat_with(
//         2,
//         2,
//         vec![int(1), int(0), int(0), int(1)],
//         vec![(1, 1, ok("all_reals", Ball::all_reals(P))?)],
//     )?;
//     match ar.solve(&vecb(vec![int(1), int(1)])?, P) {
//         Err(BoundError::Resource { .. }) => Ok(()),
//         other => Err(format!("all-reals pivot: want Resource, got {other:?}")),
//     }
// }
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA.3: error order prec > shape > NotSymmetric > value, and b_len in every LA solve"]
// fn la_error_order_and_b_len() -> Result<(), String> {
//     let ns = mat(2, 3, vec![int(-1), int(5), int(0), int(1), int(-1), int(0)])?;
//     want_err(
//         "prec before shape",
//         ns.cholesky(0),
//         is_prec,
//         "InvalidArgument",
//     )?;
//     want_err(
//         "shape before symmetry",
//         ns.cholesky(P),
//         is_dm,
//         "DimensionMismatch",
//     )?;
//     let nsym_neg = mat(2, 2, vec![int(-1), int(5), int(1), int(-1)])?;
//     want_err(
//         "NotSymmetric before NotPositiveDefinite",
//         nsym_neg.cholesky(P),
//         is_nonsym,
//         "NotSymmetric",
//     )?;
//     want_err(
//         "NotSymmetric before NotPositiveDefinite (ldl)",
//         nsym_neg.ldl(P),
//         is_nonsym,
//         "NotSymmetric",
//     )?;
//     let a = mat(2, 2, vec![int(2), int(1), int(1), int(2)])?;
//     let b3 = vecb(vec![int(1), int(1), int(1)])?;
//     want_err(
//         "BallCholesky::solve b_len",
//         ok("cholesky", a.cholesky(P))?.solve(&b3, P),
//         is_blen,
//         "b_len",
//     )?;
//     want_err(
//         "BallCholesky::quad_form b_len",
//         ok("cholesky", a.cholesky(P))?.quad_form(&b3, P),
//         is_dm,
//         "DimensionMismatch",
//     )?;
//     want_err("BallMat::solve b_len", a.solve(&b3, P), is_blen, "b_len")?;
//     want_err(
//         "BallLu::solve b_len",
//         ok("lu", a.lu(P))?.solve(&b3, P),
//         is_blen,
//         "b_len",
//     )?;
//     want_err(
//         "solve_mat b_len",
//         a.solve_mat(&mat(3, 1, vec![int(1), int(1), int(1)])?, P),
//         is_blen,
//         "b_len",
//     )?;
//     want_err(
//         "BallQr::lstsq b_len",
//         ok("qr", a.qr(P))?.lstsq(&b3, P),
//         is_blen,
//         "b_len",
//     )?;
//     // shape errors come before the b_len check of a singular system
//     let z = mat(2, 3, vec![int(0), int(0), int(0), int(0), int(0), int(0)])?;
//     want_err(
//         "solve rect",
//         z.solve(&b3, P),
//         is_rect,
//         "DimensionMismatch rectangular",
//     )
// }
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA.4: exactly-certified refusal rows return the named error (onebyone, zero2, indefinite2, nonsym2, neg_def3)"]
// fn la_certified_refusal_rows() -> Result<(), String> {
//     let rows: [(&str, usize, Vec<ExactNum>); 6] = [
//         ("onebyone_neg", 1, vec![int(-1)]),
//         ("onebyone_zero", 1, vec![int(0)]),
//         ("zero2", 2, vec![int(0), int(0), int(0), int(0)]),
//         ("indefinite2", 2, vec![int(1), int(2), int(2), int(1)]),
//         (
//             "psd_singular3",
//             3,
//             ldlt(3, &unit_lower(3, 2)?, &[int(1), int(4), int(0)])?,
//         ),
//         (
//             "neg_def3",
//             3,
//             ldlt(3, &unit_lower(3, 3)?, &[int(-1), int(-2), int(-4)])?,
//         ),
//     ];
//     for (name, n, g) in rows {
//         want_err(
//             name,
//             mat(n, n, g)?.cholesky(P),
//             is_npd,
//             "NotPositiveDefinite a",
//         )?;
//     }
//     want_err(
//         "nonsym2",
//         mat(2, 2, vec![int(2), int(1), int(0), int(2)])?.cholesky(P),
//         is_nonsym,
//         "NotSymmetric",
//     )
// }
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA.5: containment at 53 bits of 2^±203-scale logdet_spd and solve, re-decided at 2P+64 / 8P"]
// fn la_containment_low_precision() -> Result<(), String> {
//     let ks = alt(N, 202)?;
//     let (_, c) = spd(&ks)?;
//     let a = mat(N, N, c)?;
//     let ch = ok("cholesky at 53", a.cholesky(53))?;
//     want_k_ln2(
//         "logdet_spd at 53",
//         &ok("logdet_spd", ch.logdet_spd(53))?,
//         sum_i64(&ks)?,
//     )?;
//     let x = xs(N)?;
//     let b = mul_exact(N, N, 1, &mids(&a)?, &x)?;
//     let got: BallVec = ok("solve at 53", ch.solve(&vecb(b)?, 53))?;
//     want_vec("solve at 53", &got, &x)
// }
//
// // ---- 16.LA2.1, 16.LA2.2, 16.LA.1, 16.LA.6 ---------------------------------------------------------------------------
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA2.1: exact pivot ties at steps 0 and 1 pick the lowest row index (perm [1, 0, 2])"]
// fn la2_pivot_tie_lowest_index() -> Result<(), String> {
//     // Column 0 ties rows 1 and 2 at |2|; after the swap, column 1 ties positions 1 and 2 at |1|.
//     let a = vec![
//         int(0),
//         int(1),
//         int(0),
//         int(2),
//         int(0),
//         int(1),
//         int(-2),
//         int(1),
//         int(1),
//     ];
//     let lu = ok("lu tie3", mat(3, 3, a)?.lu(P))?;
//     if lu.perm() != [1, 0, 2] || lu.perm_sign() != -1 {
//         return Err(format!("tie perm {:?} sign {}", lu.perm(), lu.perm_sign()));
//     }
//     let l = vec![
//         int(1),
//         int(0),
//         int(0),
//         int(0),
//         int(1),
//         int(0),
//         int(-1),
//         int(1),
//         int(1),
//     ];
//     let u = vec![
//         int(2),
//         int(0),
//         int(1),
//         int(0),
//         int(1),
//         int(0),
//         int(0),
//         int(0),
//         int(2),
//     ];
//     want_mat("tie L", lu.l(), 3, 3, &l)?;
//     want_mat("tie U", lu.u(), 3, 3, &u)?;
//     want_encloses("tie det", &ok("tie det", lu.det(P))?, &int(-4))
// }
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA2.2: lu is square only: 2x3 and 3x2 inputs give DimensionMismatch"]
// fn la2_lu_square_only() -> Result<(), String> {
//     for (n, m) in [(2, 3), (3, 2)] {
//         let r = mat(n, m, [1, 2, 3, 4, 5, 6].into_iter().map(int).collect())?;
//         want_err(&format!("lu {n}x{m}"), r.lu(P), is_dm, "DimensionMismatch")?;
//     }
//     Ok(())
// }
//
// fn is_npd_any(e: &BoundError) -> bool {
//     matches!(e, BoundError::NotPositiveDefinite { .. })
// }
// fn is_sing_any(e: &BoundError) -> bool {
//     matches!(e, BoundError::Singular { .. })
// }
//
// /// One 16.LA.1 refusal kind: the call's outcome, the variant, the op and (where the CONTRACT names it) the detail.
// struct Refusal {
//     kind: &'static str,
//     got: Result<(), BoundError>,
//     variant: fn(&BoundError) -> bool,
//     op: &'static str,
//     detail: Option<&'static str>,
// }
//
// fn la_refusal_rows() -> Result<Vec<Refusal>, String> {
//     let neg = mat(1, 1, vec![int(-1)])?;
//     let nonsym = mat(2, 2, vec![int(1), int(1), int(0), int(1)])?;
//     let zero = mat(2, 2, vec![int(0), int(0), int(0), int(0)])?;
//     let col0 = mat(2, 2, vec![int(1), int(0), int(1), int(0)])?;
//     let wide = mat(1, 2, vec![int(1), int(1)])?;
//     let x = mat(2, 1, vec![int(1), int(1)])?;
//     let c = mat(2, 2, vec![int(-1), int(0), int(0), int(1)])?;
//     let b = vecb(vec![int(1), int(1)])?;
//     Ok(vec![
//         Refusal {
//             kind: "not_spd",
//             got: neg.cholesky(P).map(drop),
//             variant: is_npd_any,
//             op: "BallMat::cholesky",
//             detail: Some("a"),
//         },
//         Refusal {
//             kind: "sigma_not_spd",
//             got: gls(&x, &c, &b, P).map(drop),
//             variant: is_npd_any,
//             op: "gls",
//             detail: Some("sigma"),
//         },
//         Refusal {
//             kind: "singular",
//             got: zero.solve(&b, P).map(drop),
//             variant: is_sing_any,
//             op: "BallMat::solve",
//             detail: Some("singular"),
//         },
//         Refusal {
//             kind: "rank_deficient",
//             got: col0.lstsq(&b, P).map(drop),
//             variant: is_sing_any,
//             op: "BallMat::lstsq",
//             detail: Some("rank_deficient"),
//         },
//         Refusal {
//             kind: "wide",
//             got: wide.qr(P).map(drop),
//             variant: is_dm,
//             op: "BallMat::qr",
//             detail: Some("wide"),
//         },
//         Refusal {
//             kind: "not_symmetric",
//             got: nonsym.cholesky(P).map(drop),
//             variant: is_nonsym,
//             op: "BallMat::cholesky",
//             detail: None,
//         },
//     ])
// }
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA.1: refusal table not_spd, sigma_not_spd, singular, rank_deficient, wide, not_symmetric"]
// fn la_refusal_mapping_table() -> Result<(), String> {
//     for r in la_refusal_rows()? {
//         let e = match r.got {
//             Ok(()) => return Err(format!("{}: want a refusal, got Ok", r.kind)),
//             Err(e) => e,
//         };
//         let detail_ok = r.detail.is_none_or(|d| e.detail() == d);
//         if !(r.variant)(&e) || e.op() != r.op || !detail_ok {
//             return Err(format!(
//                 "{}: want {} {:?}, got {e:?}",
//                 r.kind, r.op, r.detail
//             ));
//         }
//     }
//     Ok(())
// }
//
// /// `A = [[3, 1], [1, s]]` with `s = (2^40 − 1)/3·2^-40`, so `det = 3s − 1 = −2^-40` and `b = (1, 0)` gives
// /// `x = (−(2^40 − 1)/3, 2^40)`. `l_21 = 1/3` is inexact and `u_22 = s − 1/3` cancels 42 bits, so `x_2` carries a relative
// /// error of about `2^(42 − w)`: the bound `rad ≤ 2^(42 − p)` (= `2^(2 − p)·‖x‖_∞`) holds only when `w ≥ p + 42` or so, and
// /// 16.LA.6's `w = p + 64 + 2·⌈log2 3⌉ = p + 68` meets it with room to spare.
// fn la6_system() -> Result<(BallMat, BallVec), String> {
//     let a = mat(
//         2,
//         2,
//         vec![int(3), int(1), int(1), dy(366_503_875_925, -40)?],
//     )?;
//     Ok((a, vecb(vec![int(1), int(0)])?))
// }
//
// #[test]
// #[ignore = "red: CONTRACT 16.LA.6: solve at p in 64..512 is deterministic, at prec p, holds 2^40, rad <= 2^(42-p)"]
// fn la_working_precision_visible() -> Result<(), String> {
//     let (a, b) = la6_system()?;
//     let x1 = pow2(40)?;
//     for prec in [64, 128, 256, 512] {
//         let got = ok("solve", a.solve(&b, prec))?;
//         let again = ok("solve again", a.solve(&b, prec))?;
//         if format!("{got:?}") != format!("{again:?}") {
//             return Err(format!("prec {prec}: two calls differ"));
//         }
//         want_encloses(&format!("x[1] at {prec}"), vat(&got, 1)?, &x1)?;
//         let k = 43_i64.checked_sub(i64_of(prec)?).ok_or("k")?;
//         let bound = Mag::finite(TOP, k).map_err(|e| format!("Mag::finite(2^{k}): {e:?}"))?;
//         for i in 0..2 {
//             let xi = vat(&got, i)?;
//             if xi.prec() != prec || xi.rad() > &bound {
//                 return Err(format!("prec {prec}: x[{i}] = {xi:?} exceeds 2^({k} - 1)"));
//             }
//         }
//     }
//     Ok(())
// }
