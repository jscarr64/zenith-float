//! M32: `MAX_PRECISION` (2^30 - 2^18 since M32-R, was 2^32 - 1) and the 2 GiB single-buffer size check (edge list m32-edge-list.md).

use zenith_float_num::{
    check_alloc_bytes, Consts, Error, ExactNum, Radix, RoundingMode, MAX_ALLOC_BYTES, MAX_PRECISION,
};

const GIB2: usize = 1 << 31;

#[test]
fn max_precision_value() {
    assert_eq!(MAX_PRECISION, (1 << 30) - (1 << 18));
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
    // 2^28 64-bit words = 2^31 bytes exactly (M32-R: no longer tied to 4 * MAX_PRECISION, see m32-mul-trace.md).
    let words_4p: usize = 1 << 28;
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

// Heavy rows: each mantissa at MAX_PRECISION is 16,773,120 words = 128 MiB. Run under popslot --mem 10G.
#[test]
fn heavy_max_minus_one_and_max_accepted() {
    let a = ExactNum::from_u8(1, MAX_PRECISION - 1);
    assert!(!a.is_nan(), "p = MAX_PRECISION-1 must be accepted");
    assert!(a.precision().is_some_and(|p| p >= MAX_PRECISION - 1));
    drop(a);
    let three = ExactNum::from_u8(3, MAX_PRECISION);
    assert!(!three.is_nan(), "p = MAX_PRECISION must be accepted");
    assert!(three.precision().is_some_and(|p| p >= MAX_PRECISION));
    assert!(three == ExactNum::from_u8(3, 64), "value must be exactly 3");
}

#[test]
fn heavy_binary_conversion_fits_at_max() {
    // Radix-2 conversion keeps one byte per mantissa bit: 1,073,479,680 bytes at MAX_PRECISION, inside the
    // budget since M32-R (it was 2^32 bytes and refused at the old MAX).
    let three = ExactNum::from_u8(3, MAX_PRECISION);
    assert!(!three.is_nan());
    let mut cc = Consts::new().expect("consts");
    assert_eq!(
        three
            .convert_to_radix(Radix::Bin, RoundingMode::ToEven, &mut cc)
            .map(|_| ()),
        Ok(())
    );
}
