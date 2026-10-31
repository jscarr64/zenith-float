//! M32-R: `MAX_PRECISION` is the largest p at which add, sub, mul, div, sqrt and the conversions run
//! inside the 2^31-byte single-buffer budget (trace: m32-mul-trace.md). mul at p = 2^32 - 1 needed a
//! 6,643,780,704-byte FFT buffer; the widest fitting operand is 16,773,120 words, p = 2^30 - 2^18.
//!
//! Heavy rows are `#[ignore]`: run in release under `popslot --mem 10G` with `--ignored`.

use zenith_float_num::{
    check_alloc_bytes, Consts, Error, ExactNum, Radix, RoundingMode, MAX_ALLOC_BYTES, MAX_PRECISION,
};

const P_FIT: usize = (1 << 30) - (1 << 18);

#[test]
fn max_precision_is_largest_fitting_p() {
    assert_eq!(MAX_PRECISION, P_FIT);
    assert_eq!(MAX_PRECISION % 64, 0, "a whole number of words");
}

#[test]
fn over_new_max_is_invalid_argument() {
    assert_eq!(ExactNum::new(P_FIT + 1).err(), Some(Error::InvalidArgument));
    assert_eq!(ExactNum::from_u8(1, P_FIT + 1).err(), Some(Error::InvalidArgument));
}

#[test]
fn single_buffer_boundary() {
    assert_eq!(check_alloc_bytes(MAX_ALLOC_BYTES, 1), Ok(()));
    assert_eq!(check_alloc_bytes(MAX_ALLOC_BYTES + 1, 1), Err(Error::MemoryAllocation));
    assert_eq!(check_alloc_bytes(MAX_ALLOC_BYTES / 8, 8), Ok(()));
    assert_eq!(check_alloc_bytes(MAX_ALLOC_BYTES / 8 + 1, 8), Err(Error::MemoryAllocation));
}

/// 1, u = 2^-p, a = 1 - u (all p bits set), b = 1 - 2u, all exact at p = MAX_PRECISION.
fn parts() -> (usize, ExactNum, ExactNum, ExactNum, ExactNum) {
    let p = MAX_PRECISION;
    let ep = i64::try_from(p).unwrap_or(i64::MAX);
    let one = ExactNum::from_u8(1, p);
    let u = ExactNum::pow2(-ep, p).unwrap_or_else(|e| ExactNum::nan(Some(e)));
    let a = one.clone().sub(&u, p, RoundingMode::ToEven);
    let b = a.clone().sub(&u, p, RoundingMode::ToEven);
    for (n, x) in [("one", &one), ("u", &u), ("a", &a), ("b", &b)] {
        assert_eq!(x.err(), None, "{n} must be built at p = MAX_PRECISION");
        assert!(!x.inexact(), "{n} must be exact");
    }
    (p, one, u, a, b)
}

fn ok(name: &str, got: ExactNum, want: &ExactNum) {
    assert_eq!(got.err(), None, "{name}: refused at MAX_PRECISION");
    assert!(&got == want, "{name}: wrong value at MAX_PRECISION");
}

#[test]
#[ignore = "heavy: release, popslot --mem 10G"]
fn heavy_add_sub_at_max() {
    let (p, one, u, a, _b) = parts();
    ok("a + u", a.clone().add(&u, p, RoundingMode::ToEven), &one);
    ok("1 - a", one.clone().sub(&a, p, RoundingMode::ToEven), &u);
}

#[test]
#[ignore = "heavy: release, popslot --mem 10G"]
fn heavy_mul_at_max() {
    // a^2 = 1 - 2u + u^2: u^2 is below half an ulp (ulp = u in [1/2, 1)).
    let (p, _one, _u, a, b) = parts();
    ok("a*a ToEven", a.clone().mul(&a, p, RoundingMode::ToEven), &b);
    ok("a*a Down", a.clone().mul(&a, p, RoundingMode::Down), &b);
    ok("a*a Up", a.clone().mul(&a, p, RoundingMode::Up), &a);
}

#[test]
#[ignore = "heavy: release, popslot --mem 10G"]
fn heavy_div_at_max() {
    // 1/a = 1 + u + u^2 + ...: ulp is 2u in [1, 2), so the tail is just over half an ulp.
    let (p, one, u, _a, _b) = parts();
    let a = one.clone().sub(&u, p, RoundingMode::ToEven);
    let up = one.clone().add(&u, p, RoundingMode::ToEven).add(&u, p, RoundingMode::ToEven);
    ok("1/a ToEven", one.clone().div(&a, p, RoundingMode::ToEven), &up);
    ok("1/a Down", one.clone().div(&a, p, RoundingMode::Down), &one);
}

#[test]
#[ignore = "heavy: release, popslot --mem 10G"]
fn heavy_sqrt_at_max() {
    // sqrt(1 - 2u) = 1 - u - u^2/2 - ...: just below a, far above a - u.
    let (p, _one, _u, a, b) = parts();
    ok("sqrt(b) ToEven", b.sqrt(p, RoundingMode::ToEven), &a);
    ok("sqrt(b) Down", b.sqrt(p, RoundingMode::Down), &b);
}

#[test]
#[ignore = "heavy: release, popslot --mem 10G"]
fn heavy_conversions_at_max() {
    let (p, _one, _u, a, _b) = parts();
    let mut cc = match Consts::new() {
        Ok(c) => c,
        Err(e) => panic!("consts: {e:?}"),
    };
    for rdx in [Radix::Bin, Radix::Hex, Radix::Dec] {
        let s = a.format(rdx, RoundingMode::ToEven, &mut cc);
        assert!(s.is_ok(), "{rdx:?} format refused at MAX_PRECISION: {:?}", s.as_ref().err());
        let s = s.unwrap_or_default();
        let back = ExactNum::parse(&s, rdx, p, RoundingMode::ToEven, &mut cc);
        ok(&format!("{rdx:?} round trip"), back, &a);
    }
}
