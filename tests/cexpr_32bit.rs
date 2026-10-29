//! `cexpr!` with a cancellation count that does not fit a 32-bit `usize`.

use zenith_float::{cexpr, ExactComplex};
use zenith_float_num::{ctx::Context, Consts, ExactNum, RoundingMode, EXPONENT_MAX, EXPONENT_MIN};

fn big(k: i64) -> ExactNum {
    ExactNum::from_u8(1, 64).ldexp(k, 64, RoundingMode::None)
}

/// Inexact operands at `2^(2^40)` cancelling to an exact zero: the cancelled-bit count is `2^40 + 2`.
/// On 64-bit the loop retries at that working precision and the zero is exact there. On 32-bit the
/// count does not fit `usize`, no retry is possible, and `cexpr!` must return an error (NaN), never
/// accept the unretried zero.
#[test]
fn cexpr_unrepresentable_cancellation_is_an_error() -> Result<(), Box<dyn std::error::Error>> {
    let mut ctx = Context::new(
        64,
        RoundingMode::ToEven,
        Consts::new()?,
        EXPONENT_MIN,
        EXPONENT_MAX,
    );
    let k = 1i64 << 40;
    let mut a = ExactComplex::new(big(k), ExactNum::new(64));
    a.set_inexact(true);
    let b = ExactComplex::new(-big(k), ExactNum::new(64));
    let z = cexpr!(a + b, &mut ctx);
    if usize::try_from(k + 2).is_ok() {
        assert!(z.re().is_zero() && z.im().is_zero(), "{z:?}");
    } else {
        assert!(z.is_nan(), "{z:?}");
    }
    Ok(())
}
