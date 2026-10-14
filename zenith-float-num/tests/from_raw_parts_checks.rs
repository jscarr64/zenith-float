//! `ExactNum::from_raw_parts` enforces its documented checks (2.0.0 audit, defect 5).

use zenith_float_num::{Error, ExactNum, Sign, Word, EXPONENT_MAX, EXPONENT_MIN, WORD_BIT_SIZE};

fn assert_invalid(r: &ExactNum, what: &str) {
    assert!(r.is_nan(), "Expected NaN for {}, got {:?}", what, r);
    assert_eq!(
        r.err(),
        Some(Error::InvalidArgument),
        "Expected InvalidArgument error for {}, got {:?}",
        what,
        r.err()
    );
}

#[test]
fn from_raw_parts_round_trips_valid_parts() {
    let ps = [64, 128, 256, 1024];
    let vs = [1u64, 3, 0xFFFF_FFFF_FFFF_FFFF, 0x8000_0000_0000_0001];
    let es = [EXPONENT_MIN, EXPONENT_MIN + 1, -1, 0, 1, EXPONENT_MAX - 1, EXPONENT_MAX];
    let negatives = [false, true];

    for &p in &ps {
        for &v in &vs {
            for &e in &es {
                for &negative in &negatives {
                    let mut x = ExactNum::from_u64(v, p);
                    x.set_exponent(e);
                    if negative {
                        x = x.neg();
                    }
                    let (m, n, s, ex, inexact) = x.as_raw_parts().unwrap();
                    let y = ExactNum::from_raw_parts(m, n, s, ex, inexact);
                    assert_eq!(
                        x, y,
                        "Round-trip failed for p={}, v={}, e={}, negative={}",
                        p, v, e, negative
                    );
                    assert!(
                        !y.is_nan(),
                        "Round-trip produced NaN for p={}, v={}, e={}, negative={}",
                        p,
                        v,
                        e,
                        negative
                    );
                }
            }
        }
    }

    // Test min_positive and its negation
    for &p in &ps {
        let x = ExactNum::min_positive(p);
        let (m, n, s, ex, inexact) = x.as_raw_parts().unwrap();
        let y = ExactNum::from_raw_parts(m, n, s, ex, inexact);
        assert_eq!(x, y, "Round-trip failed for min_positive p={}", p);
        assert!(
            !y.is_nan(),
            "Round-trip produced NaN for min_positive p={}",
            p
        );

        let x_neg = ExactNum::min_positive(p).neg();
        let (m, n, s, ex, inexact) = x_neg.as_raw_parts().unwrap();
        let y = ExactNum::from_raw_parts(m, n, s, ex, inexact);
        assert_eq!(
            x_neg, y,
            "Round-trip failed for negated min_positive p={}",
            p
        );
        assert!(
            !y.is_nan(),
            "Round-trip produced NaN for negated min_positive p={}",
            p
        );
    }
}

#[test]
fn from_raw_parts_rejects_inconsistent_parts() {
    let ps = [64, 128, 256, 1024];
    let signs = [Sign::Pos, Sign::Neg];

    for &p in &ps {
        // (a) normal mantissa words of from_u64(3,p) with n = p + 1 -> invalid.
        {
            let x = ExactNum::from_u64(3, p);
            let (m, _n, s, ex, inexact) = x.as_raw_parts().unwrap();
            let y = ExactNum::from_raw_parts(m, p + 1, s, ex, inexact);
            assert_invalid(&y, &format!("normal mantissa with n=p+1 for p={}", p));
        }

        // (b) normal mantissa (top bit set) with n = p - 1 -> invalid.
        {
            let x = ExactNum::from_u64(3, p);
            let (m, _n, s, ex, inexact) = x.as_raw_parts().unwrap();
            let y = ExactNum::from_raw_parts(m, p - 1, s, ex, inexact);
            assert_invalid(&y, &format!("normal mantissa with n=p-1 for p={}", p));
        }

        // (c) subnormal: mantissa words all zero except the lowest word = 1 (so 1 significant bit, top bit clear), n = 1
        {
            let mut m: Vec<Word> = vec![0; p / WORD_BIT_SIZE];
            m[0] = 1;
            for &e in &[EXPONENT_MIN + 1, -1, 0, 1, EXPONENT_MAX] {
                for &s in &signs {
                    let y = ExactNum::from_raw_parts(&m, 1, s, e, false);
                    assert_invalid(
                        &y,
                        &format!(
                            "subnormal with e != EXPONENT_MIN for p={}, e={}, s={:?}",
                            p, e, s
                        ),
                    );
                }
            }

            // With e = EXPONENT_MIN the same parts must be valid and equal to min_positive(p) (negated when sign is negative)
            for &s in &signs {
                let y = ExactNum::from_raw_parts(&m, 1, s, EXPONENT_MIN, false);
                assert!(
                    !y.is_nan(),
                    "Expected valid subnormal for p={}, e=EXPONENT_MIN",
                    p
                );
                if s == Sign::Neg {
                    assert_eq!(
                        y,
                        ExactNum::min_positive(p).neg(),
                        "Subnormal with negative sign should equal negated min_positive for p={}",
                        p
                    );
                } else {
                    assert_eq!(
                        y,
                        ExactNum::min_positive(p),
                        "Subnormal with positive sign should equal min_positive for p={}",
                        p
                    );
                }
            }
        }

        // (d) mantissa words all zero with n = 1 -> invalid.
        {
            let m: Vec<Word> = vec![0; p / WORD_BIT_SIZE];
            for &s in &signs {
                let y = ExactNum::from_raw_parts(&m, 1, s, EXPONENT_MIN, false);
                assert_invalid(&y, &format!("all-zero mantissa with n=1 for p={}", p));
            }
        }
    }
}
