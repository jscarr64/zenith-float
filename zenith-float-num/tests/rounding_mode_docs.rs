//! RoundingMode docs and directed-rounding behaviour (2.0.0 audit, defect 2).

use std::fs;
use zenith_float_num::{ExactNum, RoundingMode};

#[test]
fn rounding_mode_docs_describe_directed_modes() {
    let content = fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/src/defs.rs"))
        .expect("Failed to read defs.rs");

    let start = content
        .find("pub enum RoundingMode {")
        .expect("RoundingMode enum not found");
    let end = content[start..]
        .find('}')
        .expect("RoundingMode enum not closed");

    let enum_content = &content[start..start + end + 1];

    assert!(!enum_content.contains("Round half toward positive infinity"));
    assert!(!enum_content.contains("Round half toward negative infinity"));
    assert!(!enum_content.contains("Round half toward zero"));
    assert!(!enum_content.contains("Round half away from zero"));

    assert!(enum_content.contains("/// Round toward positive infinity (directed rounding)."));
    assert!(enum_content.contains("/// Round toward negative infinity (directed rounding)."));
    assert!(enum_content.contains("/// Round toward zero (truncation; directed rounding)."));
    assert!(enum_content.contains("/// Round away from zero (directed rounding)."));
}

#[test]
fn directed_modes_are_directed_not_half() {
    for &p in &[64, 128, 192, 256] {
        for &n in &[1i8, 2, 5, 11, 127] {
            for &d in &[3i8, 7, 9, 13] {
                for &sign in &[false, true] {
                    let n = if sign { -n } else { n };
                    let a = ExactNum::from_i8(n, p);
                    let b = ExactNum::from_i8(d, p);

                    let q_up = a.div(&b, p, RoundingMode::Up);
                    let q_down = a.div(&b, p, RoundingMode::Down);
                    let q_zero = a.div(&b, p, RoundingMode::ToZero);
                    let q_from = a.div(&b, p, RoundingMode::FromZero);

                    // q_down < q_up strictly
                    assert!(
                        q_down.compare(&q_up).unwrap() < 0,
                        "For p={}, n={}, d={} (sign={}), q_down should be less than q_up",
                        p,
                        n,
                        d,
                        sign
                    );

                    // Check that results are on correct side of true quotient
                    let prod_up = q_up.mul(&b, 2 * p + 64, RoundingMode::None);
                    assert!(
                        prod_up.compare(&a).unwrap() > 0,
                        "For p={}, n={}, d={} (sign={}), q_up*d should be > n",
                        p,
                        n,
                        d,
                        sign
                    );

                    let prod_down = q_down.mul(&b, 2 * p + 64, RoundingMode::None);
                    assert!(
                        prod_down.compare(&a).unwrap() < 0,
                        "For p={}, n={}, d={} (sign={}), q_down*d should be < n",
                        p,
                        n,
                        d,
                        sign
                    );

                    // ToZero behavior
                    let expected_zero = if n > 0 { &q_down } else { &q_up };
                    assert_eq!(
                        q_zero.compare(expected_zero).unwrap(),
                        0,
                        "For p={}, n={}, d={} (sign={}), q_to_zero should equal {}",
                        p,
                        n,
                        d,
                        sign,
                        if n > 0 { "q_down" } else { "q_up" }
                    );

                    // FromZero behavior
                    let expected_from = if n > 0 { &q_up } else { &q_down };
                    assert_eq!(
                        q_from.compare(expected_from).unwrap(),
                        0,
                        "For p={}, n={}, d={} (sign={}), q_from_zero should equal {}",
                        p,
                        n,
                        d,
                        sign,
                        if n > 0 { "q_up" } else { "q_down" }
                    );
                }
            }
        }
    }
}
