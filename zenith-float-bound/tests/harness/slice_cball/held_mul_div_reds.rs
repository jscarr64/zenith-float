// HOLD H11: held, not compiled (no `mod` line names this file). Index and use-or-delete rule: HOLD.md at the repo root.
// The original text follows, each line commented out with a leading `// `; uncomment it whole to go live.
// //! Slice CB red rows for `CBall::mul` (16.CB.4) and `CBall::div` (16.CB.7), and the rows that need them (sqr vs mul,
// //! all-reals propagation, the prec sweep, gold containment and tightness for mul/div). Held (H11) because they are red:
// //! mul and div refuse with `Resource` until the borrowing once-rounded zf ops land (P12: no ignored reds in a PR).
// //! Moved verbatim from slice_cball.rs; a `mod` line in slice_cball.rs takes them live, without their `#[ignore]`s.
//
// use super::*;
//
// #[test]
// #[ignore = "red: not implemented: CONTRACT CB.4 (mul): waits zf borrowing ops"]
// fn cball_mul_rows() -> Result<(), String> {
//     // (1+2i)(3+4i) = (3-8) + (4+6)i = -5+10i.
//     want_point(
//         "mul",
//         &ok("mul", cpt(1, 2, 53)?.mul(&cpt(3, 4, 53)?, 53))?,
//         -5,
//         10,
//     )?;
//     // cancellation: (1+i)(1-i) = 2 + 0i.
//     want_point(
//         "mul cancel",
//         &ok("mul", cpt(1, 1, 53)?.mul(&cpt(1, -1, 53)?, 53))?,
//         2,
//         0,
//     )?;
//     // i·i = -1.
//     want_point(
//         "i*i",
//         &ok("mul", cpt(0, 1, 53)?.mul(&cpt(0, 1, 53)?, 53))?,
//         -1,
//         0,
//     )?;
//     // boxes, parts containing 0 and mixed signs: containment at every sample corner pair.
//     for &prec in PRECS.iter() {
//         let (ar, ai, br, bi) = (0i64, -3i64, 5i64, 0i64);
//         let p = ok(
//             "mul",
//             cbox(ar, ai, 2, prec)?.mul(&cbox(br, bi, 3, prec)?, prec),
//         )?;
//         for x in samples(ar, 2)?.iter() {
//             for y in samples(ai, 2)?.iter() {
//                 for u in samples(br, 3)?.iter() {
//                     for v in samples(bi, 3)?.iter() {
//                         let re = x.mul_full_prec(u).sub_full_prec(&y.mul_full_prec(v));
//                         let im = x.mul_full_prec(v).add_full_prec(&y.mul_full_prec(u));
//                         want_contains("mul re", p.re(), &re)?;
//                         want_contains("mul im", p.im(), &im)?;
//                     }
//                 }
//             }
//         }
//     }
//     Ok(())
// }
//
// #[test]
// #[ignore = "red: the width comparison waits CBall::mul (zf borrowing ops)"]
// fn cball_sqr_not_wider_than_mul() -> Result<(), String> {
//     // Every row: each part of sqr(a) is no wider than the same part of mul(a, a) (16.CB.5).
//     for &prec in PRECS.iter() {
//         for &(re, im, k) in [
//             (0i64, 0i64, 1i64),
//             (0, 3, 2),
//             (1, -1, 4),
//             (-7, 2, 10),
//             (5, 0, 0),
//         ]
//         .iter()
//         {
//             let a = cbox(re, im, k, prec)?;
//             let s = ok("sqr", cbox(re, im, k, prec)?.sqr(prec))?;
//             let m = ok("mul", cbox(re, im, k, prec)?.mul(&a, prec))?;
//             if s.re().rad() > m.re().rad() || s.im().rad() > m.im().rad() {
//                 return Err(format!(
//                     "sqr wider than mul at prec {} for {:?}: {:?} vs {:?}",
//                     prec, a, s, m
//                 ));
//             }
//             for x in samples(re, k)?.iter() {
//                 for y in samples(im, k)?.iter() {
//                     want_contains(
//                         "sqr re",
//                         s.re(),
//                         &x.mul_full_prec(x).sub_full_prec(&y.mul_full_prec(y)),
//                     )?;
//                     let xy = x.mul_full_prec(y);
//                     want_contains("sqr im", s.im(), &xy.add_full_prec(&xy))?;
//                 }
//             }
//         }
//     }
//     Ok(())
// }
//
// #[test]
// #[ignore = "red: not implemented: CONTRACT CB.7 (div): waits zf borrowing ops"]
// fn cball_div_rows() -> Result<(), String> {
//     // (2+4i)/(1+i) = (2+4i)(1-i)/2 = (2+4 + (4-2)i)/2 = 3+i.
//     let q = ok("div", cpt(2, 4, 53)?.div(&cpt(1, 1, 53)?, 53))?;
//     want_contains("div re", q.re(), &num(3))?;
//     want_contains("div im", q.im(), &num(1))?;
//     // by an exact real: (6-9i)/3 = 2-3i.
//     let q = ok("div", cpt(6, -9, 53)?.div(&cpt(3, 0, 53)?, 53))?;
//     want_contains("div real re", q.re(), &num(2))?;
//     want_contains("div real im", q.im(), &num(-3))?;
//     // one axis of the divisor touches 0, the other is clear: finite, not all reals.
//     let q = ok("div", cpt(1, 0, 53)?.div(&cbox(0, 4, 1, 53)?, 53))?;
//     want_finite_both("one axis touches 0", &q)
// }
//
// #[test]
// #[ignore = "red: not implemented: CONTRACT CB.7 (div all reals R1): waits zf borrowing ops"]
// fn cball_div_rectangle_containing_zero_is_all_reals() -> Result<(), String> {
//     // R1: both parts reach 0 → Ok with both parts all reals (never Domain).
//     let exact_zero = cpt(0, 0, 53)?;
//     want_all_reals_both("0/0", &ok("div", cpt(1, 1, 53)?.div(&exact_zero, 53))?)?;
//     want_all_reals_both(
//         "box at 0",
//         &ok("div", cpt(1, 1, 53)?.div(&cbox(0, 0, 3, 53)?, 53))?,
//     )?;
//     // edge |m| = r on both axes: 1 ± 1 touches 0 exactly.
//     let edge = ok(
//         "new",
//         CBall::new(ball(1, mag(TOP, 1)?, 53)?, ball(-1, mag(TOP, 1)?, 53)?),
//     )?;
//     want_all_reals_both("edge", &ok("div", cpt(1, 1, 53)?.div(&edge, 53))?)?;
//     // just clear on one axis: 2 ± 1 excludes 0, so R1 does not hold.
//     let clear = ok(
//         "new",
//         CBall::new(ball(2, mag(TOP, 1)?, 53)?, ball(-1, mag(TOP, 1)?, 53)?),
//     )?;
//     want_finite_both("clear on re", &ok("div", cpt(1, 1, 53)?.div(&clear, 53))?)
// }
//
// #[test]
// #[ignore = "red: not implemented: CONTRACT CB.7 (div all reals R2): waits zf borrowing ops"]
// fn cball_div_mid_rad_sqr_widening_is_all_reals() -> Result<(), String> {
//     // R2 (16.CB.7), kept as a defensive check: R1 does not hold, but d = sqr(b.re) + sqr(b.im) reaches 0 because S1's
//     // sqr is mid-rad (see cball_r2_row_s1_sqr_reaches_zero). Observed behavior: both parts all reals, at every prec.
//     for &prec in PRECS.iter() {
//         let b = ok("new", CBall::new(r2_divisor_re(prec)?, pt(0, prec)?))?;
//         want_all_reals_both(
//             "R2 mid-rad row",
//             &ok("div", cpt(1, 0, prec)?.div(&b, prec))?,
//         )?;
//     }
//     Ok(())
// }
//
// #[test]
// #[ignore = "red: not implemented: CONTRACT CB.7 (div exact scaling): waits zf borrowing ops"]
// fn cball_div_scaling_exact() -> Result<(), String> {
//     // The common 2^-k on both operands leaves the quotient unchanged: (2+4i)·2^K / ((1+i)·2^K) contains 3+i at
//     // K = ±2^40, where |b|^2 unscaled would leave the midpoint range.
//     for &k in [1i64 << 40, -(1i64 << 40)].iter() {
//         let s = pow2(k)?;
//         let mk = |v: i64| ok("exact", Ball::exact(num(v).mul_full_prec(&s), 64));
//         let a = ok("new", CBall::new(mk(2)?, mk(4)?))?;
//         let b = ok("new", CBall::new(mk(1)?, mk(1)?))?;
//         let q = ok("div", a.div(&b, 64))?;
//         want_finite_both("scaled div", &q)?;
//         want_contains("scaled re", q.re(), &num(3))?;
//         want_contains("scaled im", q.im(), &num(1))?;
//     }
//     Ok(())
// }
//
// #[test]
// #[ignore = "red: not implemented: CONTRACT CB.9 (all-reals propagation): waits zf borrowing ops"]
// fn cball_all_reals_propagation_mul() -> Result<(), String> {
//     // mul: both output parts use the all-reals re part.
//     let h = ok(
//         "new",
//         CBall::new(ok("all_reals", Ball::all_reals(53))?, pt(1, 53)?),
//     )?;
//     want_all_reals_both("mul", &ok("mul", h.mul(&cpt(1, 1, 53)?, 53))?)
// }
//
// #[test]
// #[ignore = "red: not implemented: CONTRACT CB.9 (one prec per call): waits zf borrowing ops"]
// fn cball_prec_sweep_mul_sqr_div() -> Result<(), String> {
//     let a = || cbox(7, -5, 3, 1024);
//     for &prec in PRECS.iter() {
//         let b = cbox(3, 2, 5, 1024)?;
//         for z in [
//             ok("mul", a()?.mul(&b, prec))?,
//             ok("sqr", a()?.sqr(prec))?,
//             ok("div", a()?.div(&b, prec))?,
//         ]
//         .iter()
//         {
//             if z.prec() != prec || z.re().prec() != prec || z.im().prec() != prec {
//                 return Err(format!("prec {}: got {:?}", prec, z));
//             }
//         }
//     }
//     Ok(())
// }
//
// #[test]
// #[ignore = "red: not implemented: CONTRACT CB.10 (golds): waits zf borrowing ops"]
// fn cball_golden_rows_contained_mul_div() -> Result<(), String> {
//     gold_contained(&["mul", "div"]).map(|_| ())
// }
//
// #[test]
// #[ignore = "red: not implemented: CONTRACT CB.10 (tightness): waits zf borrowing ops"]
// fn cball_tightness_factors_mul_div() -> Result<(), String> {
//     gold_tight(&["mul", "div"]).map(|_| ())
// }
