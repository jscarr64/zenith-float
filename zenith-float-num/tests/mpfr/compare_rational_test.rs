//! GMP `mpq` oracle for [`ExactRational`] add / sub / mul / div.

use crate::mpfr::common::{reset_test_rng, test_random};
use rug::Rational;
use zenith_float_num::ExactNum;
use zenith_float_num::ExactRational;

fn i32_nz() -> i64 {
    let v = (test_random::<i32>() % 10_000) as i64;
    if v == 0 {
        1
    } else {
        v
    }
}

fn i32_pos_den() -> i64 {
    (test_random::<i32>().unsigned_abs() % 10_000 + 1) as i64
}

fn assert_same(got: &ExactRational, want: &Rational) {
    assert!(!got.is_nan(), "ExactRational was NaN; GMP {want}");
    let wn = want
        .numer()
        .to_i64()
        .unwrap_or_else(|| panic!("GMP numer overflow {want}"));
    let wd = want
        .denom()
        .to_i64()
        .unwrap_or_else(|| panic!("GMP denom overflow {want}"));
    let p = got.num().precision().unwrap_or(128);
    assert_eq!(
        got.num().cmp(&ExactNum::from_i64(wn, p)),
        Some(0),
        "num {got:?} vs GMP {want}"
    );
    assert_eq!(
        got.den().cmp(&ExactNum::from_i64(wd, p)),
        Some(0),
        "den {got:?} vs GMP {want}"
    );
}

#[test]
fn mpfr_compare_exact_rational_gmp() {
    reset_test_rng();

    let tenth = ExactRational::parse_exact("0.1").unwrap();
    assert_same(&tenth, &Rational::from((1, 10)));
    assert_same(
        &tenth.mul(&ExactRational::from_i64(10, 1)),
        &Rational::from((1, 1)),
    );

    assert_same(
        &ExactRational::from_i64(1, 3).add(&ExactRational::from_i64(1, 6)),
        &Rational::from((1, 2)),
    );
    assert_same(
        &ExactRational::from_i64(2, 4),
        &Rational::from((1, 2)),
    );

    for _ in 0..48 {
        let n1 = i32_nz();
        let d1 = i32_pos_den();
        let n2 = i32_nz();
        let d2 = i32_pos_den();
        let a = ExactRational::from_i64(n1, d1);
        let b = ExactRational::from_i64(n2, d2);
        let ga = Rational::from((n1, d1));
        let gb = Rational::from((n2, d2));
        assert_same(&a.add(&b), &Rational::from(&ga + &gb));
        assert_same(&a.sub(&b), &Rational::from(&ga - &gb));
        assert_same(&a.mul(&b), &Rational::from(&ga * &gb));
        if !b.num().is_zero() {
            assert_same(&a.div(&b), &Rational::from(&ga / &gb));
        }
    }
}
