//! M32 runtime rows (Spence review of 34467a17). They compile against the 1.x API (no new imports),
//! so on the pre-M32 base they fail at runtime, not at compile time.

use zenith_float_num::{Consts, Error, ExactNum, Radix, RoundingMode, MAX_PRECISION};

/// `2^32 - 64`, written out so the row does not depend on the constant it checks.
const CAP: usize = 4_294_967_232;

#[test]
fn max_precision_is_2_pow_32_minus_64_not_the_old_cap() {
    assert_ne!(MAX_PRECISION, (1 << 60) - 1, "old 1.x cap 2^60 - 1");
    assert_eq!(MAX_PRECISION, CAP);
}

// Heavy rows: one mantissa at the cap is 512 MiB. Run under popslot --mem 10G.

#[test]
fn heavy_over_2gib_request_is_memory_allocation_through_the_root_path() {
    // Radix-2 conversion needs one byte per mantissa bit: 2^32 bytes, twice the 2 GiB budget.
    let three = ExactNum::from_u8(3, CAP);
    assert!(!three.is_nan());
    let mut cc = Consts::new().expect("consts");
    assert_eq!(
        three
            .convert_to_radix(Radix::Bin, RoundingMode::ToEven, &mut cc)
            .map(|_| ()),
        Err(Error::MemoryAllocation)
    );
}

#[test]
fn heavy_guard_bits_and_retry_near_the_cap_refuse_with_memory_allocation() {
    // An unaligned precision just under the cap needs a guard word, and the guarded width passes the cap.
    let p = CAP - 1;
    let three = ExactNum::from_u8(3, p);
    assert!(!three.is_nan());
    let sq = ExactNum::from_u8(3, p).mul(&three, p, RoundingMode::ToEven);
    assert_eq!(sq.err(), Some(Error::MemoryAllocation));
    let mut cc = Consts::new().expect("consts");
    let l = three.ln(p, RoundingMode::ToEven, &mut cc);
    assert_eq!(l.err(), Some(Error::MemoryAllocation));
}
