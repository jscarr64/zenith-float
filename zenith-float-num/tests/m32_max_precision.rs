//! M32: `MAX_PRECISION = 2^32 - 1` and the 2 GiB single-buffer size check (edge list m32-edge-list.md).

use zenith_float_num::{
    check_alloc_bytes, Consts, Error, ExactNum, Radix, RoundingMode, MAX_ALLOC_BYTES, MAX_PRECISION,
};

const GIB2: usize = 1 << 31;

#[test]
fn max_precision_is_u32_max() {
    assert_eq!(u64::try_from(MAX_PRECISION).ok(), Some(u64::from(u32::MAX)));
}

#[test]
fn over_max_is_invalid_argument() {
    assert_eq!(
        ExactNum::new(MAX_PRECISION + 1).err(),
        Some(Error::InvalidArgument)
    );
}

#[test]
fn budget_is_2_pow_31() {
    assert_eq!(MAX_ALLOC_BYTES, GIB2);
}

#[test]
fn size_check_bytes() {
    assert_eq!(check_alloc_bytes(GIB2 - 1, 1), Ok(()));
    assert_eq!(check_alloc_bytes(GIB2, 1), Ok(()));
    assert_eq!(check_alloc_bytes(GIB2 + 1, 1), Err(Error::MemoryAllocation));
}

#[test]
fn size_check_words_threshold() {
    // 4 * MAX_PRECISION bits rounds up to exactly 2^28 64-bit words = 2^31 bytes.
    let words_4p = (4 * MAX_PRECISION).div_ceil(64);
    assert_eq!(words_4p, 1 << 28);
    assert_eq!(check_alloc_bytes(words_4p - 1, 8), Ok(()));
    assert_eq!(check_alloc_bytes(words_4p, 8), Ok(()));
    assert_eq!(
        check_alloc_bytes(words_4p + 1, 8),
        Err(Error::MemoryAllocation)
    );
}

#[test]
fn size_check_overflow() {
    assert_eq!(check_alloc_bytes(usize::MAX, 2), Err(Error::MemoryAllocation));
    assert_eq!(check_alloc_bytes(usize::MAX, 0), Ok(()));
}

// Heavy rows: each mantissa at MAX_PRECISION is 2^26 words = 512 MiB. Run under popslot --mem 10G.
#[test]
fn heavy_max_minus_one_and_max_accepted() {
    let a = ExactNum::new(MAX_PRECISION - 1);
    assert!(!a.is_nan(), "p = MAX_PRECISION-1 must be accepted");
    assert!(a.precision().is_some_and(|p| p >= MAX_PRECISION - 1));
    drop(a);
    let three = ExactNum::from_u8(3, MAX_PRECISION);
    assert!(!three.is_nan(), "p = MAX_PRECISION must be accepted");
    assert!(three.precision().is_some_and(|p| p >= MAX_PRECISION));
    assert!(three == ExactNum::from_u8(3, 64), "value must be exactly 3");
}

#[test]
fn heavy_binary_conversion_over_budget_refused() {
    // Radix-2 conversion keeps one byte per mantissa bit: 2^32 bytes at MAX_PRECISION, over the budget.
    let three = ExactNum::from_u8(3, MAX_PRECISION);
    assert!(!three.is_nan());
    let mut cc = Consts::new().expect("consts");
    assert_eq!(
        three
            .convert_to_radix(Radix::Bin, RoundingMode::ToEven, &mut cc)
            .map(|_| ()),
        Err(Error::MemoryAllocation)
    );
}
