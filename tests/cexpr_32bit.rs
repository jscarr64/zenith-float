//! `cexpr!` with a cancellation count that does not fit a 32-bit `usize`.

use zenith_float::{cexpr, ExactComplex};
use zenith_float_num::{
    ctx::Context, Consts, Error, ExactNum, RoundingMode, EXPONENT_MAX, EXPONENT_MIN,
};

/// `x / y - x / y` with `x = 2^(2^39)`, `y = 3`: the quotients are inexact at every precision (cexpr!
/// leaves and products are exact), and they cancel to an exact zero (exponent 0) from exponent about
/// `2^39`. The cancelled-bit count is not a 32-bit `usize`, so no retry is possible and `cexpr!` must
/// return NaN with `InvalidArgument` (the error for an unrepresentable precision), never the unretried
/// zero. On 64-bit the retry would ask for about `2^39` working bits, so the test only runs where the
/// count does not fit.
#[test]
fn cexpr_unrepresentable_cancellation_is_an_error() -> Result<(), Box<dyn std::error::Error>> {
    let k = 1i64 << 39;
    if usize::try_from(k + 2).is_ok() {
        return Ok(());
    }
    let mut ctx = Context::new(
        64,
        RoundingMode::ToEven,
        Consts::new()?,
        EXPONENT_MIN,
        EXPONENT_MAX,
    );
    let x = ExactComplex::new(
        ExactNum::from_u8(1, 64).ldexp(k, 64, RoundingMode::None),
        ExactNum::new(64),
    );
    let y = ExactComplex::new(ExactNum::from_u8(3, 64), ExactNum::new(64));
    let z = cexpr!(x / y - x / y, &mut ctx);
    assert!(z.is_nan(), "{z:?}");
    assert_eq!(z.re().err(), Some(Error::InvalidArgument), "{z:?}");
    Ok(())
}
