//! Binary-splitting bit-burst exponential, after FLINT 3.6
//! `arb_exp_arf_bb` and `_arb_exp_sum_bs_simple` (Fredrik Johansson,
//! LGPL-3.0-or-later).
//!
//! The main wrapper keeps the existing series path below
//! [`ARB_EXP_BB_MIN_PREC`] bits (matching FLINT's `exp_arf` dispatch).
//! Above that, the argument is reduced by \(2^q\), converted to fixed point,
//! and exponentiated by the serial bit-burst loop: each chunk's
//! \(\exp(u/2^r)-1\) is a binary-split Taylor sum, the product is formed,
//! and \(q\) squarings restore \(\exp(x)\). The parallel / powtab variants
//! are not ported; the simple serial path is the one-thread FLINT route.

use crate::defs::Error;
use crate::defs::RoundingMode;
use crate::defs::Sign;
use crate::integer::ExactInt;
use crate::num::ExactNumNumber;
use crate::WORD_BIT_SIZE;
use alloc::borrow::Cow;

/// Minimum destination precision at which [`exp_arf_bb`] is used, after
/// FLINT 3.6 `arb_exp_arf` (`prec >= 10000` and not the tiny-argument
/// rectangular-split path).
pub(crate) const ARB_EXP_BB_MIN_PREC: usize = 10000;

/// Dest widths at or above [`ARB_EXP_BB_MIN_PREC`]: bit-burst is used for
/// every finite payload. A previous carve-out left full-width arguments
/// at dest p=2^20 on the series path because the first burst was slower
/// (~33 s vs ~20 s) when `int_to_num` still accumulated one limb at a
/// time. After the linear pack, bit-burst is the right path at that
/// width too.
pub(crate) fn use_exp_arf_bb(x: &ExactNumNumber, p: usize) -> bool {
    let _ = x;
    p >= ARB_EXP_BB_MIN_PREC
}

const ARGRED_BITS: i64 = 16;
const START_BITS: i64 = 32;

/// Upper bounds on \(\log_2(1/n!)\) for \(n < 256\), copied from FLINT 3.6
/// `exp_taylor_bound.c`.
const REC_FAC_BOUND_2EXP: [i16; 256] = [
    0, 0, -1, -2, -4, -6, -9, -12, -15, -18, -21, -25, -28, -32, -36, -40, -44, -48, -52, -56, -61,
    -65, -69, -74, -79, -83, -88, -93, -97, -102, -107, -112, -117, -122, -127, -132, -138, -143,
    -148, -153, -159, -164, -169, -175, -180, -186, -191, -197, -202, -208, -214, -219, -225, -231,
    -237, -242, -248, -254, -260, -266, -272, -278, -284, -289, -295, -302, -308, -314, -320, -326,
    -332, -338, -344, -350, -357, -363, -369, -375, -382, -388, -394, -401, -407, -413, -420, -426,
    -433, -439, -446, -452, -458, -465, -472, -478, -485, -491, -498, -504, -511, -518, -524, -531,
    -538, -544, -551, -558, -564, -571, -578, -585, -591, -598, -605, -612, -619, -626, -632, -639,
    -646, -653, -660, -667, -674, -681, -688, -695, -702, -709, -716, -723, -730, -737, -744, -751,
    -758, -765, -772, -779, -786, -793, -801, -808, -815, -822, -829, -836, -844, -851, -858, -865,
    -872, -880, -887, -894, -901, -909, -916, -923, -931, -938, -945, -952, -960, -967, -975, -982,
    -989, -997, -1004, -1011, -1019, -1026, -1034, -1041, -1049, -1056, -1064, -1071, -1078, -1086,
    -1093, -1101, -1108, -1116, -1123, -1131, -1139, -1146, -1154, -1161, -1169, -1176, -1184,
    -1192, -1199, -1207, -1214, -1222, -1230, -1237, -1245, -1253, -1260, -1268, -1276, -1283,
    -1291, -1299, -1306, -1314, -1322, -1329, -1337, -1345, -1353, -1360, -1368, -1376, -1384,
    -1391, -1399, -1407, -1415, -1423, -1430, -1438, -1446, -1454, -1462, -1470, -1477, -1485,
    -1493, -1501, -1509, -1517, -1525, -1532, -1540, -1548, -1556, -1564, -1572, -1580, -1588,
    -1596, -1604, -1612, -1620, -1628, -1636, -1644, -1652, -1660, -1668, -1675,
];

fn bit_length_usize(n: usize) -> usize {
    (usize::BITS - n.leading_zeros()) as usize
}

fn rec_fac_bound_2exp_si(n: i64) -> i64 {
    if n < 256 {
        return i64::from(REC_FAC_BOUND_2EXP[n as usize]);
    }
    // Lower bound on log2(n!): n*floor(log2 n) - 3n/2 for n >= 256.
    // rec_fac = -that is an upper bound on log2(1/n!), so the term count
    // is conservative (never too few terms).
    let bitlen = 64 - n.leading_zeros() as i64;
    (3 * n) / 2 - n * (bitlen - 1)
}

/// Fewest \(N\) (exclusive) such that the tail of \(\exp(x)\) after \(x^N/N!\)
/// is below \(2^{-(\mathrm{prec}+1)}\), after `_arb_exp_taylor_bound`.
fn exp_taylor_bound(mag: i64, prec: i64) -> i64 {
    if mag > -2 {
        // After argument reduction |x| < 2^{-16}; this branch is a
        // conservative fallback, not a panic, if a caller slips.
        return prec.saturating_add(16).max(2);
    }
    let mut i: i64 = 1;
    loop {
        let rec = rec_fac_bound_2exp_si(i);
        let lhs = (mag as i128) * (i as i128) + (rec as i128);
        if lhs < -(prec as i128) - 1 {
            return i;
        }
        i = i.saturating_add(1);
        if i > prec.saturating_add(4) {
            return i;
        }
    }
}

fn bs_num_terms(mag: i64, prec: i64) -> i64 {
    let mut n = exp_taylor_bound(mag, prec) - 1;
    if n < 0 {
        n = 0;
    }
    if n > 10000 {
        n = ((n + 127) / 128) * 128;
    } else if n > 1000 {
        n = ((n + 15) / 16) * 16;
    } else if n > 100 {
        n = ((n + 1) / 2) * 2;
    }
    n
}

/// \(T, Q, Q_{\exp}\) such that \((e^{x/2^r}-1) = T/Q \cdot 2^{-Q_{\exp}}\).
fn bsplit(
    x: &ExactInt,
    r: i64,
    a: i64,
    b: i64,
    cont: bool,
) -> Result<(ExactInt, ExactInt, ExactInt, i64), Error> {
    if b - a == 1 {
        let p = x.clone();
        let q = ExactInt::from_u64((a + 1) as u64);
        let t = p.clone();
        return Ok((p, t, q, r));
    }
    let m = a + (b - a) / 2;
    let (p1, t1, q1, qe1) = bsplit(x, r, a, m, true)?;
    let (p2, t2, q2, qe2) = bsplit(x, r, m, b, true)?;
    // T = T1 * Q2 * 2^{Q2exp} + P1 * T2
    let mut t = int_mul(&t1, &q2)?;
    if qe2 > 0 {
        t = t.shl(qe2 as usize);
    }
    t = t.add(&int_mul(&p1, &t2)?);
    let q = int_mul(&q1, &q2)?;
    let qe = qe1 + qe2;
    let p = if cont { int_mul(&p1, &p2)? } else { ExactInt::zero() };
    Ok((p, t, q, qe))
}

fn exp_sum_bs_simple(x: &ExactInt, r: i64, n: i64) -> Result<(ExactInt, ExactInt, i64), Error> {
    let (_p, t, q, qe) = bsplit(x, r, 0, n, false)?;
    Ok((t, q, qe))
}

fn int_to_num(i: &ExactInt, p: usize) -> Result<ExactNumNumber, Error> {
    i.to_exact_num_number(p, RoundingMode::None)
}

fn int_mul(a: &ExactInt, b: &ExactInt) -> Result<ExactInt, Error> {
    a.mul_fast(b)
}

/// Fixed-point integer \(t \approx x \cdot 2^{-\mathrm{exp}}\), truncated
/// toward zero. `exp` is the FLINT `arf_get_fmpz_fixed_si` exponent
/// (typically \(-wp+q\)).
fn to_fixed_si(x: &ExactNumNumber, exp: i64) -> ExactInt {
    let p = x.mantissa_max_bit_len() as i64;
    let shift = x.exponent() - p - exp;
    let sign = if x.is_negative() { Sign::Neg } else { Sign::Pos };
    let mut t = ExactInt::from_le_words(sign, x.mantissa().digits());
    if shift > 0 {
        t = t.shl(shift as usize);
    } else if shift < 0 {
        t = t.shr_towards_zero((-shift) as usize);
    }
    t
}

/// \(\exp(x)\) at working precision `p` by the serial bit-burst path.
/// The caller has already refused huge/tiny inputs that the main `exp`
/// wrapper handles.
pub(crate) fn exp_arf_bb(x: &ExactNumNumber, p: usize) -> Result<ExactNumNumber, Error> {
    if x.is_zero() {
        return ExactNumNumber::from_word(1, p);
    }

    let mag = x.exponent();
    // FLINT refuses mag > 200 or mag < -2p-100; the wrapper already
    // covers those. A defensive clamp keeps q and wp in range.
    if mag > 200 || mag < -2 * (p as i64) - 100 {
        if mag > 200 {
            return Err(Error::ExponentOverflow(x.sign()));
        }
        let mut one = ExactNumNumber::from_word(1, p)?;
        one.set_inexact(true);
        return Ok(one);
    }

    let q = if mag + ARGRED_BITS > 0 { mag + ARGRED_BITS } else { 0 };
    let wp = (p as i64) + 10 + 2 * q + 2 * (bit_length_usize(p) as i64);
    let wp_usize = wp as usize;

    let mut t = to_fixed_si(x, -wp + q);
    let mut z = ExactNumNumber::from_word(1, wp_usize)?;
    let mut bits = START_BITS;

    while !t.is_zero() {
        let r = if bits < wp { bits } else { wp };
        // When wp <= r the chunk u is all of t (borrowed, no copy).
        let u = t.shr_towards_zero_ref(if wp > r { (wp - r) as usize } else { 0 });

        if !u.is_zero() {
            let umag = (u.bit_length() as i64) - r;
            let n = bs_num_terms(umag, wp);
            let one_plus = if n <= 0 {
                ExactNumNumber::from_word(1, wp_usize)?
            } else {
                let (tt, qq, qexp) = exp_sum_bs_simple(&u, r, n)?;
                // (T/Q) * 2^{-Qexp} via the mantissa divider, then + 1.
                let work = wp_usize.saturating_add(WORD_BIT_SIZE);
                let mut tnum = int_to_num(&tt, work)?;
                tnum.set_exponent(tnum.exponent() - qexp);
                let qnum = int_to_num(&qq, work)?;
                let w = tnum.div(&qnum, wp_usize, RoundingMode::None)?;
                let one = ExactNumNumber::from_word(1, wp_usize)?;
                one.add(&w, wp_usize, RoundingMode::None)?
            };
            z = z.mul(&one_plus, wp_usize, RoundingMode::None)?;
        }

        // Remove the used bits: t -= u << (wp - r)
        if !u.is_zero() {
            t = match u {
                Cow::Owned(u) => t.sub(&u.shl((wp - r) as usize)),
                // u was all of t
                Cow::Borrowed(_) => ExactInt::zero(),
            };
        }
        bits = bits.saturating_mul(2);
        if bits > wp && t.is_zero() {
            break;
        }
        // If bits has grown past wp and t is still nonzero, the next
        // iteration uses r = wp and consumes t in one step.
        if bits > wp * 4 && !t.is_zero() {
            // Should be unreachable; consume the remainder as one chunk.
            let n = bs_num_terms((t.bit_length() as i64) - wp, wp);
            if n > 0 {
                let (tt, qq, qexp) = exp_sum_bs_simple(&t, wp, n)?;
                let work = wp_usize.saturating_add(WORD_BIT_SIZE);
                let mut tnum = int_to_num(&tt, work)?;
                tnum.set_exponent(tnum.exponent() - qexp);
                let qnum = int_to_num(&qq, work)?;
                let w = tnum.div(&qnum, wp_usize, RoundingMode::None)?;
                let one = ExactNumNumber::from_word(1, wp_usize)?;
                let one_plus = one.add(&w, wp_usize, RoundingMode::None)?;
                z = z.mul(&one_plus, wp_usize, RoundingMode::None)?;
            }
            break;
        }
    }

    // exp(x) = exp(x/2^q)^{2^q}
    for _ in 0..q {
        z = z.mul(&z, wp_usize, RoundingMode::None)?;
    }

    z.set_precision(p, RoundingMode::None)?;
    // exp(x) for a nonzero x is never a dyadic number; the sum above can still land exactly on a rounding
    // midpoint (one Taylor term for a tiny x: 1 + x), which the caller must not round as an exact tie.
    z.set_inexact(true);
    Ok(z)
}

#[cfg(test)]
mod tests {
    use super::{exp_taylor_bound, rec_fac_bound_2exp_si};

    /// `exp(2^-p) = 1 + 2^-p + 2^(-2p-1) + ...` lies just above the midpoint `1 + 2^-p` of `1` and `1 + 2^(1-p)`.
    /// Its one-term bit-burst sum is exactly that midpoint, so the result must be flagged inexact or the caller
    /// rounds the "exact tie" to even (`1`). Gold: `1 + 2^(1-p)` to nearest / up, `1` down (python-flint `arb`
    /// Ziv loop, tools/mpfr_golden/gen.py `exp 10048 N +1p-10048`; the MPFR test is retry_budget
    /// `exp_just_above_a_tie_at_100032_bits`). p = 10048 is the first word width on the bit-burst path.
    #[test]
    fn exp_of_a_tiny_power_of_two_is_not_taken_for_an_exact_tie() -> Result<(), String> {
        use crate::{Consts, ExactNum, Exponent, RoundingMode};
        let p = 10048usize;
        let mut cc = Consts::new().map_err(|e| format!("{e:?}"))?;
        let x = ExactNum::from_u8(1, p).ldexp(-(p as Exponent), p, RoundingMode::ToEven);
        let one = ExactNum::from_u8(1, p);
        let above = ExactNum::from_u8(1, p).add(
            &ExactNum::from_u8(1, p).ldexp(1 - p as Exponent, p, RoundingMode::ToEven),
            p,
            RoundingMode::ToEven,
        );
        for (rm, want) in [
            (RoundingMode::ToEven, &above),
            (RoundingMode::Up, &above),
            (RoundingMode::Down, &one),
        ] {
            let got = x.exp(p, rm, &mut cc);
            if got != *want {
                return Err(format!("exp(2^-{p}) {rm:?}: got {got:?}"));
            }
        }
        Ok(())
    }

    /// Golds: floor(log2(1/n!)) = floor(-lgamma(n+1)/ln 2) from FLINT (python-flint 0.9.0, `arb(n+1).lgamma()`,
    /// 256 bits), cross-checked with mpmath 1.3.0 `loggamma` (pinned oracle venv); script and output:
    /// ~/pop-logs/paul/2026-11-15/land39/final/pr40-fixes/r2nd/mutkill/recfac_flint.{py,out}.
    const LOG2_REC_FAC_FLOOR: [(i64, i64); 5] =
        [(256, -1684), (300, -2042), (1000, -8530), (4096, -43251), (100000, -1516705)];

    /// n >= 256 (the formula branch): the bound never under-states log2(1/n!) (so the Taylor term count is never too
    /// small) and is within 1.1 n bits of it (n floor(log2 n) - 3n/2 loses at most about n bits). Kills the surviving
    /// mutant exp_bb.rs:73:21 `*` -> `/`, which returns a positive "bound" (352 at n = 256).
    #[test]
    fn rec_fac_bound_formula_branch_is_a_tight_upper_bound() -> Result<(), String> {
        for (n, floor) in LOG2_REC_FAC_FLOOR {
            let rec = rec_fac_bound_2exp_si(n);
            if rec <= floor {
                return Err(format!(
                    "n={n}: bound {rec} is not above log2(1/n!) (floor {floor})"
                ));
            }
            let slack = n + n / 10;
            if rec > floor + 1 + slack {
                return Err(format!("n={n}: bound {rec} is more than 1.1 n = {slack} bits above log2(1/n!) (floor {floor})"));
            }
        }
        Ok(())
    }

    /// Golds: the fewest N with mag*N + log2(1/N!) < -(prec+1), from exact log2(N!) (FLINT `arb.lgamma`, mpmath
    /// cross-check; same script). The returned count may exceed it (the bound is conservative) by at most 12 %.
    const TERMS_OPT: [(i64, i64, i64); 4] =
        [(-16, 6000, 266), (-16, 20000, 825), (-20, 50000, 1707), (-2, 4000, 430)];

    #[test]
    fn exp_taylor_bound_term_count_is_conservative_and_near_optimal() -> Result<(), String> {
        for (mag, prec, opt) in TERMS_OPT {
            let n = exp_taylor_bound(mag, prec);
            if n < opt {
                return Err(format!(
                    "mag={mag} prec={prec}: {n} terms < the {opt} needed"
                ));
            }
            let most = opt + (opt * 12) / 100;
            if n > most {
                return Err(format!(
                    "mag={mag} prec={prec}: {n} terms > {most} (12 % over the optimal {opt})"
                ));
            }
        }
        Ok(())
    }
}
