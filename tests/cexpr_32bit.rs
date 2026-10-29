//! `cexpr!` with a cancellation count that does not fit a 32-bit `usize`.

use zenith_float::{cexpr, ExactComplex};
use zenith_float_num::{
    ctx::Context, Consts, Error, ExactNum, RoundingMode, EXPONENT_MAX, EXPONENT_MIN,
};

/// `x * x - x * x` with `x = (1 + 2^-63) * 2^(2^39)`: each square needs 127 bits, so at the first working
/// precision both products are inexact, and they cancel to an exact zero (exponent 0) from exponent
/// `2^40 + 1`. The cancelled-bit count `2^40 + 2` is not a 32-bit `usize`, so no retry is possible
/// and `cexpr!` must return NaN with `InvalidArgument` (the error for an unrepresentable precision),
/// never the unretried zero. On 64-bit the retry asks for `2^40` working bits (a 128 GiB mantissa),
/// so the test only runs where the count does not fit.
#[test]
fn cexpr_unrepresentable_cancellation_is_an_error() -> Result<(), Box<dyn std::error::Error>> {
    let k = 1i64 << 39;
    if usize::try_from(2 * k + 2).is_ok() {
        return Ok(());
    }
    let mut ctx = Context::new(
        64,
        RoundingMode::ToEven,
        Consts::new()?,
        EXPONENT_MIN,
        EXPONENT_MAX,
    );
    let one = ExactNum::from_u8(1, 64);
    let ulp = ExactNum::from_u8(1, 64).ldexp(-63, 64, RoundingMode::None);
    let x = ExactComplex::new(
        one.add(&ulp, 64, RoundingMode::None)
            .ldexp(k, 64, RoundingMode::None),
        ExactNum::new(64),
    );
    let z = cexpr!(x * x - x * x, &mut ctx);
    assert!(z.is_nan(), "{z:?}");
    assert_eq!(z.re().err(), Some(Error::InvalidArgument), "{z:?}");
    Ok(())
}
