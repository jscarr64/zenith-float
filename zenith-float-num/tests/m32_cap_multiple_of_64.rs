//! M32 (Spence review of d12fa916, item D): `MAX_PRECISION = 2^32 - 64`, a whole number of 64-bit words,
//! so a user precision in `(2^32 - 64, 2^32 - 1]` is out of range (`InvalidArgument`), not a working
//! width past the cap (`MemoryAllocation`).

use zenith_float_num::{Consts, Error, ExactNum, Radix, RoundingMode, MAX_PRECISION};

/// `2^32 - 64`, written out so the rows do not depend on the constant they check.
const CAP: usize = 4_294_967_232;

/// The previous cap `2^32 - 1`, and the precisions between the two caps that are now refused.
const OLD_CAP: usize = 4_294_967_295;
const REFUSED: [usize; 4] = [CAP + 1, CAP + 32, OLD_CAP - 1, OLD_CAP];

#[test]
fn max_precision_is_2_pow_32_minus_64_a_whole_number_of_words() {
    assert_eq!(MAX_PRECISION, CAP);
    assert_eq!(MAX_PRECISION % 64, 0);
}

#[test]
fn round_exact_refuses_p_between_the_new_and_the_old_cap() {
    let mut cc = Consts::new().expect("consts");
    for p in REFUSED {
        let x = ExactNum::parse("3", Radix::Dec, p, RoundingMode::ToEven, &mut cc);
        assert_eq!(x.err(), Some(Error::InvalidArgument), "parse at p = {p}");
        let two = ExactNum::from_u8(2, 64);
        assert_eq!(
            two.powi(3, p, RoundingMode::ToEven).err(),
            Some(Error::InvalidArgument),
            "powi at p = {p}"
        );
    }
}

#[test]
fn arithmetic_and_constructors_refuse_p_between_the_new_and_the_old_cap() {
    let three = ExactNum::from_u8(3, 64);
    for p in REFUSED {
        assert_eq!(
            ExactNum::from_u8(3, p).err(),
            Some(Error::InvalidArgument),
            "from_u8 at p = {p}"
        );
        assert_eq!(
            ExactNum::mul_ref(&three, &three, p, RoundingMode::ToEven).err(),
            Some(Error::InvalidArgument),
            "mul at p = {p}"
        );
    }
}
