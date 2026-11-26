//! Harper S2 wide goldens: every row of tests/goldens/slices/bridge/bridge_wide.txt (scripts/gen_bridge_wide_goldens.py:
//! seeded random 2048-, 10_000- and 100_000-bit mantissas at cosmology binades (1e+-61, 1e+-120, 1e+-300), near 0, near
//! EXPONENT_MAX/EXPONENT_MIN and uniform over the range, 10^k rounded to P bits, and cancellation pairs) is asserted exactly
//! (CONTRACT 16.S2.1, 16.S2.2, 16.S2.4). N rows: `Mag::cmp_exact_num`; D rows: `Mag::cmp_exact_diff` (rows with |a - b|
//! above EXPONENT_MAX are the known red of harper_s2_diff_overflow.rs and are not generated); U rows: `Mag::from_exact_num`
//! and `Mag::from_exact_num_lower`. Each kind's row count is pinned.
//! COVERS: Mag::cmp_exact_num, Mag::cmp_exact_diff, Mag::from_exact_num, Mag::from_exact_num_lower

use zenith_float_bound::mag::Mag;

use super::harper_s2_goldens::{exact, mag, ord, rows_in, verdict};

const FILE: &str = "bridge_wide.txt";
const HEADER: &str = "# zenith-float-bound S2 rows above 1028 bits, scripts/gen_bridge_wide_goldens.py v1 (see its docstring)";

#[test]
fn s2_wide_goldens_cmp_exact_num_rows() -> Result<(), String> {
    let mut bad = Vec::new();
    for (line, row) in rows_in(FILE, HEADER, "N", 206)? {
        let f: Vec<&str> = row.split(' ').collect();
        let [r, x, o] = f.as_slice() else {
            return Err(format!("line {}: bad N row", line));
        };
        let got = mag(r)?.cmp_exact_num(&exact(x)?);
        if got != ord(o)? {
            bad.push(format!("line {}: N got {:?}", line, got));
        }
    }
    verdict("N", &bad)
}

#[test]
fn s2_wide_goldens_cmp_exact_diff_rows() -> Result<(), String> {
    let mut bad = Vec::new();
    for (line, row) in rows_in(FILE, HEADER, "D", 206)? {
        let f: Vec<&str> = row.split(' ').collect();
        let [r, a, b, o] = f.as_slice() else {
            return Err(format!("line {}: bad D row", line));
        };
        let got = mag(r)?.cmp_exact_diff(&exact(a)?, &exact(b)?);
        if got != ord(o)? {
            bad.push(format!("line {}: D got {:?}", line, got));
        }
    }
    verdict("D", &bad)
}

#[test]
fn s2_wide_goldens_from_exact_num_rows() -> Result<(), String> {
    let mut bad = Vec::new();
    for (line, row) in rows_in(FILE, HEADER, "U", 103)? {
        let f: Vec<&str> = row.split(' ').collect();
        let [x, up, lo] = f.as_slice() else {
            return Err(format!("line {}: bad U row", line));
        };
        let v = exact(x)?;
        let (gu, gl) = (Mag::from_exact_num(&v), Mag::from_exact_num_lower(&v));
        if gu != mag(up)? || gl != mag(lo)? {
            bad.push(format!("line {}: U got {:?} {:?}", line, gu, gl));
        }
    }
    verdict("U", &bad)
}
