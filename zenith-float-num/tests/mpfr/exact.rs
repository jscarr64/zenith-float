//! The exact MPFR operations add, sub, mul, div and remainder, computed here with exact big-integer arithmetic instead
//! of being stored in the golden streams (Paul call 2026-11-24 9:13 PM CT, r6 item 6: the stored answers were most of the
//! 1.5 GB). A line-for-line port of `tools/mpfr_golden/gen.py` (`round_mag`, `add_num`, `op_add` ... `op_remainder`, `fmt`):
//! it returns the answer text gen.py would print, and in gen mode (`mpfr-golden-gen`) `gold.rs` checks every such answer
//! against gen.py's, record for record, before the stream is written.

use crate::mpfr::mp::{cmp_nat, Big, Nat};
use std::cmp::Ordering;

/// Exponent beyond which gen.py prints `huge` / `tiny` (`HUGE` in gen.py).
const HUGE: i128 = (1i128 << 62) + (1i128 << 32);

/// gen.py's values: `(-1)^neg * m * 2^e` with `m` odd.
enum V {
    Nan,
    Inf(bool),
    Zero(bool),
    Num(bool, Big, i128),
}

fn num(neg: bool, m: Big, e: i128) -> V {
    assert!(!m.is_zero(), "exact: num of zero");
    let tz = m.trailing_zeros();
    V::Num(neg, m.shr(tz), e + i128::from(tz))
}

fn top(m: &Big, e: i128) -> i128 {
    e + i128::from(m.bits())
}

fn shift(n: i128) -> Result<u64, String> {
    u64::try_from(n).map_err(|e| format!("exact: shift {n}: {e}"))
}

fn sgn(neg: bool, d: i32) -> i32 {
    if neg {
        -d
    } else {
        d
    }
}

/// `round_mag`: `|v| = n 2^e`, plus `sticky` for nonzero bits below `n`; rounded to `prec` bits.
fn round_mag(
    neg: bool,
    n: Big,
    e: i128,
    prec: u64,
    rnd: u8,
    sticky: bool,
) -> Result<(V, i32), String> {
    let l = n.bits();
    assert!(l > 0, "exact: round_mag of zero");
    if sticky {
        assert!(l >= prec + 2, "exact: sticky rounding needs prec+2 bits");
    }
    if l <= prec {
        return Ok((num(neg, n, e), 0));
    }
    let k = l - prec;
    let mut q = n.shr(k);
    let rb = n.bit(k - 1);
    let sb = sticky || n.low_nonzero(k - 1);
    let e2 = e + i128::from(k);
    if !rb && !sb {
        return Ok((num(neg, q, e2), 0));
    }
    let up = match rnd {
        b'N' => rb && (sb || q.bit(0)),
        b'Z' => false,
        b'U' => !neg,
        b'D' => neg,
        b'A' => true,
        _ => return Err(format!("exact: rounding mode {}", char::from(rnd))),
    };
    if up {
        q = q.add(&Big::from_u64(1));
    }
    Ok((num(neg, q, e2), sgn(neg, if up { 1 } else { -1 })))
}

/// `round_signed`: `(neg, |n|) 2^e`; a zero result is `zero(zero_neg)`.
fn round_signed(
    neg: bool,
    n: Big,
    e: i128,
    prec: u64,
    rnd: u8,
    zero_neg: bool,
) -> Result<(V, i32), String> {
    if n.is_zero() {
        return Ok((V::Zero(zero_neg), 0));
    }
    round_mag(neg, n, e, prec, rnd, false)
}

fn round_val(v: V, prec: u64, rnd: u8) -> Result<(V, i32), String> {
    match v {
        V::Num(neg, m, e) => round_mag(neg, m, e, prec, rnd, false),
        other => Ok((other, 0)),
    }
}

/// `add_num`: x + y for two numbers, exact, one rounding.
fn add_num(
    x: (bool, &Big, i128),
    y: (bool, &Big, i128),
    prec: u64,
    rnd: u8,
) -> Result<(V, i32), String> {
    let (big, small) = if top(x.1, x.2) >= top(y.1, y.2) { (x, y) } else { (y, x) };
    let g = (i128::from(prec) + 3 - i128::from(big.1.bits())).max(0);
    if top(small.1, small.2) < big.2 - g {
        // small lies strictly below big's extended last unit: a sticky bit
        let n = big.1.shl(shift(g + 1)?);
        let n = if big.0 == small.0 {
            n.add(&Big::from_u64(1))
        } else {
            n.sub(&Big::from_u64(1))
        };
        return round_mag(big.0, n, big.2 - g - 1, prec, rnd, false);
    }
    let e = x.2.min(y.2);
    let a = x.1.shl(shift(x.2 - e)?);
    let b = y.1.shl(shift(y.2 - e)?);
    let zero_neg = rnd == b'D';
    if x.0 == y.0 {
        return round_signed(x.0, a.add(&b), e, prec, rnd, zero_neg);
    }
    match cmp_nat(&a, &b) {
        Ordering::Greater => round_signed(x.0, a.sub(&b), e, prec, rnd, zero_neg),
        Ordering::Less => round_signed(y.0, b.sub(&a), e, prec, rnd, zero_neg),
        Ordering::Equal => Ok((V::Zero(zero_neg), 0)),
    }
}

fn op_add(prec: u64, rnd: u8, x: V, y: V) -> Result<(V, i32), String> {
    match (x, y) {
        (V::Nan, _) | (_, V::Nan) => Ok((V::Nan, 0)),
        (V::Inf(a), V::Inf(b)) if a != b => Ok((V::Nan, 0)),
        (V::Inf(a), _) | (_, V::Inf(a)) => Ok((V::Inf(a), 0)),
        (V::Zero(a), V::Zero(b)) => Ok((V::Zero((a && b) || (a != b && rnd == b'D')), 0)),
        (V::Zero(_), y) => round_val(y, prec, rnd),
        (x, V::Zero(_)) => round_val(x, prec, rnd),
        (V::Num(xn, xm, xe), V::Num(yn, ym, ye)) => {
            add_num((xn, &xm, xe), (yn, &ym, ye), prec, rnd)
        }
    }
}

fn neg_v(v: V) -> V {
    match v {
        V::Nan => V::Nan,
        V::Inf(n) => V::Inf(!n),
        V::Zero(n) => V::Zero(!n),
        V::Num(n, m, e) => V::Num(!n, m, e),
    }
}

fn sign(v: &V) -> bool {
    match v {
        V::Nan => false,
        V::Inf(n) | V::Zero(n) | V::Num(n, _, _) => *n,
    }
}

fn op_mul(prec: u64, rnd: u8, x: V, y: V) -> Result<(V, i32), String> {
    let neg = sign(&x) != sign(&y);
    match (x, y) {
        (V::Nan, _) | (_, V::Nan) => Ok((V::Nan, 0)),
        (V::Inf(_), V::Zero(_)) | (V::Zero(_), V::Inf(_)) => Ok((V::Nan, 0)),
        (V::Inf(_), _) | (_, V::Inf(_)) => Ok((V::Inf(neg), 0)),
        (V::Zero(_), _) | (_, V::Zero(_)) => Ok((V::Zero(neg), 0)),
        (V::Num(_, xm, xe), V::Num(_, ym, ye)) => {
            round_val(num(neg, xm.mul(&ym), xe + ye), prec, rnd)
        }
    }
}

fn op_div(prec: u64, rnd: u8, x: V, y: V) -> Result<(V, i32), String> {
    let neg = sign(&x) != sign(&y);
    match (x, y) {
        (V::Nan, _) | (_, V::Nan) => Ok((V::Nan, 0)),
        (V::Inf(_), V::Inf(_)) => Ok((V::Nan, 0)),
        (V::Inf(_), _) => Ok((V::Inf(neg), 0)),
        (_, V::Inf(_)) => Ok((V::Zero(neg), 0)),
        (V::Zero(_), V::Zero(_)) => Ok((V::Nan, 0)),
        (_, V::Zero(_)) => Ok((V::Inf(neg), 0)),
        (V::Zero(_), _) => Ok((V::Zero(neg), 0)),
        (V::Num(_, a, xe), V::Num(_, b, ye)) => {
            let s = (i128::from(prec) + 3 + i128::from(b.bits()) - i128::from(a.bits())).max(0);
            let (q, r) = a.shl(shift(s)?).divrem(&b);
            round_mag(neg, q, xe - ye - s, prec, rnd, !r.is_zero())
        }
    }
}

/// `2^k mod m`, `m > 0`.
fn pow2_mod(k: u128, m: &Big) -> Big {
    let mut r = Big::from_u64(1).divrem(m).1;
    let mut b = Big::from_u64(2).divrem(m).1;
    let mut k = k;
    while k > 0 {
        if k & 1 == 1 {
            r = r.mul(&b).divrem(m).1;
        }
        k >>= 1;
        if k > 0 {
            b = b.mul(&b).divrem(m).1;
        }
    }
    r
}

fn op_remainder(prec: u64, rnd: u8, x: V, y: V) -> Result<(V, i32), String> {
    match (x, y) {
        (V::Nan, _) | (_, V::Nan) | (V::Inf(_), _) | (_, V::Zero(_)) => Ok((V::Nan, 0)),
        (x, V::Inf(_)) | (x @ V::Zero(_), _) => round_val(x, prec, rnd),
        (V::Num(xn, mx, ex), V::Num(_, my, ey)) => {
            if top(&mx, ex) < top(&my, ey) - 1 {
                return round_val(V::Num(xn, mx, ex), prec, rnd); // |x| < |y|/2: quotient 0
            }
            let (xx, yy, e) = if ex >= ey {
                let m2 = my.shl(1);
                let k = u128::try_from(ex - ey)
                    .map_err(|e| format!("exact: remainder exponent: {e}"))?;
                (mx.mul(&pow2_mod(k, &m2)).divrem(&m2).1, my, ey)
            } else {
                let yy = my.shl(shift(ey - ex)?);
                (mx.divrem(&yy.shl(1)).1, yy, ex)
            };
            let odd = cmp_nat(&xx, &yy) != Ordering::Less;
            let r = xx.divrem(&yy).1;
            let twice = r.shl(1);
            let (rneg, rmag) = match cmp_nat(&twice, &yy) {
                Ordering::Greater => (true, yy.sub(&r)),
                Ordering::Equal if odd => (true, yy.sub(&r)),
                _ => (false, r),
            };
            if rmag.is_zero() {
                return Ok((V::Zero(xn), 0));
            }
            round_signed(rneg != xn, rmag, e, prec, rnd, false)
        }
    }
}

/// gen.py's `fmt`.
fn fmt(v: &V) -> String {
    match v {
        V::Nan => "nan".to_owned(),
        V::Inf(n) => format!("{}inf", if *n { '-' } else { '+' }),
        V::Zero(n) => format!("{}0", if *n { '-' } else { '+' }),
        V::Num(n, m, e) => {
            let s = if *n { '-' } else { '+' };
            let t = top(m, *e);
            if t > HUGE {
                format!("{s}huge")
            } else if t < -HUGE {
                format!("{s}tiny")
            } else {
                format!("{s}{}p{e}", m.to_hex())
            }
        }
    }
}

/// gen.py's `parse_val` of an operand token.
fn parse_val(s: &str) -> Result<V, String> {
    if s == "nan" {
        return Ok(V::Nan);
    }
    let neg = match s.as_bytes().first() {
        Some(b'-') => true,
        Some(b'+') => false,
        _ => return Err(format!("exact: operand {s}")),
    };
    let b = &s[1..];
    match b {
        "inf" => return Ok(V::Inf(neg)),
        "0" => return Ok(V::Zero(neg)),
        _ => {}
    }
    let (m, e) = b
        .split_once('p')
        .ok_or_else(|| format!("exact: operand {s}"))?;
    if m.is_empty() || !m.bytes().all(|c| c.is_ascii_hexdigit()) {
        return Err(format!("exact: operand {s}"));
    }
    let e: i128 = e.parse().map_err(|x| format!("exact: exponent {e}: {x}"))?;
    Ok(num(neg, Big::from_hex(m), e))
}

/// Ops answered here (the rest are stored in the golden stream).
pub const OPS: [&str; 5] = ["add", "sub", "mul", "div", "remainder"];

/// gen.py's answer to `query` when it is one of `OPS` on two float operands; `None` for every other query.
pub fn answer(query: &str) -> Result<Option<String>, String> {
    let mut it = query.split(' ');
    let (Some(op), Some(prec), Some(rnd), Some(a), Some(b), None) = (
        it.next(),
        it.next(),
        it.next(),
        it.next(),
        it.next(),
        it.next(),
    ) else {
        return Ok(None);
    };
    if !OPS.contains(&op) || [a, b].iter().any(|t| t.starts_with(['i', 'u', 's'])) {
        return Ok(None);
    }
    let prec: u64 = prec
        .parse()
        .map_err(|e| format!("exact: precision {prec}: {e}"))?;
    let rnd = match rnd.as_bytes() {
        [r] => *r,
        _ => return Err(format!("exact: rounding mode {rnd}")),
    };
    let (x, y) = (parse_val(a)?, parse_val(b)?);
    let (v, t) = match op {
        "add" => op_add(prec, rnd, x, y)?,
        "sub" => op_add(prec, rnd, x, neg_v(y))?,
        "mul" => op_mul(prec, rnd, x, y)?,
        "div" => op_div(prec, rnd, x, y)?,
        _ => op_remainder(prec, rnd, x, y)?,
    };
    if let V::Num(_, m, _) = &v {
        assert!(
            m.bits() <= prec,
            "exact: {op} result wider than {prec} bits"
        );
    }
    Ok(Some(format!("{} {t}", fmt(&v))))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ans(q: &str) -> String {
        crate::mpfr::gold::must(answer(q).and_then(|a| a.ok_or_else(|| format!("{q}: not exact"))))
    }

    #[test]
    fn known_answers_from_gen_py() {
        // gen.py --answer on the same queries (oracle-2026-10-30)
        assert_eq!(ans("add 53 N +1p0 +1p-60"), "+1p0 -1");
        assert_eq!(ans("add 53 U +1p0 +1p-60"), "+10000000000001p-52 1");
        assert_eq!(ans("sub 53 D +1p0 +1p0"), "-0 0");
        assert_eq!(ans("sub 53 N +1p0 +1p0"), "+0 0");
        assert_eq!(ans("div 10 N +1p0 +3p0"), "+2abp-11 1");
        assert_eq!(ans("mul 2 Z +3p0 +3p0"), "+1p3 -1");
        assert_eq!(ans("remainder 53 N +7p0 +2p0"), "-1p0 0");
        assert_eq!(ans("remainder 53 N +5p0 +2p0"), "+1p0 0");
        assert_eq!(ans("remainder 53 N +1p100 +3p0"), "+1p0 0");
        assert_eq!(ans("div 53 N +1p0 +0"), "+inf 0");
        assert_eq!(ans("mul 53 N -inf +0"), "nan 0");
        assert_eq!(answer("sqr 53 N +3p0"), Ok(None));
    }

    /// The word-by-word division against its definition: q d + r = n, r < d, on operands that hit the qhat corrections.
    #[test]
    fn divrem_matches_its_definition() {
        let mut s = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            s
        };
        for round in 0..3000u32 {
            let nl = 1 + (next() % 12) as usize;
            let dl = 1 + (next() % 8) as usize;
            let pick = |r: u64, x: u64| match r % 4 {
                0 => u64::MAX,
                1 => 1u64 << 63,
                2 => 0,
                _ => x,
            };
            let n: Vec<u64> = (0..nl)
                .map(|_| {
                    let r = next();
                    pick(r, next())
                })
                .collect();
            let mut d: Vec<u64> = (0..dl)
                .map(|_| {
                    let r = next();
                    pick(r, next())
                })
                .collect();
            if d.iter().all(|&w| w == 0) {
                d[0] = 1 + u64::from(round);
            }
            let (n, d) = (Big::from_vec(n), Big::from_vec(d));
            let (q, r) = n.divrem(&d);
            assert_eq!(
                cmp_nat(&r, &d),
                Ordering::Less,
                "remainder not below divisor: {n:?} / {d:?}"
            );
            assert_eq!(q.mul(&d).add(&r), n, "q d + r != n: {n:?} / {d:?}");
        }
    }
}
