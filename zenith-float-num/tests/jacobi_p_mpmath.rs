//! Jacobi `P_n^{(α,β)}` vs mpmath 1.4.1 (`mpmath.jacobi`) at 50 digits / 2p+64 bits.
//! Covers the Legendre reduction and non-trivial (including non-integer) parameters.

use zenith_float_num::{Consts, ExactNum, Radix, RoundingMode};

const RM: RoundingMode = RoundingMode::ToEven;

/// (n, α, β, x, reference)
const CASES: &[(u32, &str, &str, &str, &str)] = &[
    (5, "0", "0", "0.3", "0.34538625000000000000000000000000000000000000000000"),
    (
        6,
        "0.5",
        "-0.25",
        "0.4",
        "0.41621648451232910156250000000000000000000000000000",
    ),
    (
        8,
        "1.5",
        "2.25",
        "-0.7",
        "-1.2123341660099519649520516395568847656250000000000",
    ),
    (
        7,
        "-0.4",
        "1.2",
        "-0.15",
        "0.29815395762164100000000000000000000000000000000000",
    ),
    (3, "2.5", "0.5", "1.5", "47.031250000000000000000000000000000000000000000000"),
    (
        10,
        "0.75",
        "1.25",
        "0.9",
        "-1.0842854531406208038330078125000000000000000000000",
    ),
    (2, "1", "1", "0.5", "0.18750000000000000000000000000000000000000000000000"),
];

fn bits(got: &ExactNum, want: &ExactNum, p: usize) -> i32 {
    if got.is_nan() {
        return i32::MIN;
    }
    let d = got.sub(want, 2 * p, RoundingMode::None);
    if d.is_zero() {
        return i32::MAX;
    }
    want.exponent().unwrap_or(0) - d.exponent().unwrap_or(0)
}

#[test]
fn jacobi_p_mpmath() {
    let mut cc = Consts::new().unwrap();
    let parse = |s: &str, p: usize, cc: &mut Consts| ExactNum::parse(s, Radix::Dec, p, RM, cc);
    let mut failures = Vec::new();
    for &(n, a, b, x, re) in CASES {
        for &p in &[128usize, 256] {
            let xv = parse(x, p, &mut cc);
            let av = parse(a, p, &mut cc);
            let bv = parse(b, p, &mut cc);
            let want = parse(re, 4 * p, &mut cc);
            let got = xv.jacobi_p(n as usize, &av, &bv, p, RM);
            let bts = bits(&got, &want, p);
            if bts < (p as i32) - 4 {
                failures.push(format!(
                    "P_{n}^({a},{b})({x}) p={p}: {bts} bits {:?}",
                    got.err()
                ));
            }
        }
    }
    // Awkward non-integer α=1/3, β=2/3 built from integers (not a decimal parse).
    for &p in &[128usize, 256] {
        let one = ExactNum::from_u8(1, p);
        let three = ExactNum::from_u8(3, p);
        let four = ExactNum::from_u8(4, p);
        let alpha = one.div(&three, p, RM);
        let beta = ExactNum::from_u8(2, p).div(&three, p, RM);
        let x = one.div(&four, p, RM);
        let want = parse(
            "0.256488313400205761316872427983539094650205761316872427983539094650205761316872427983539095",
            4 * p,
            &mut cc,
        );
        let got = x.jacobi_p(4, &alpha, &beta, p, RM);
        let bts = bits(&got, &want, p);
        if bts < (p as i32) - 4 {
            failures.push(format!("P_4^(1/3,2/3)(1/4) p={p}: {bts} bits"));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
