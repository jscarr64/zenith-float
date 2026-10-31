//! API-freeze unit 5 (Paul/Rusty 2026-10-23): every public prec is `Precision` (= zf's `usize` alias), validated
//! at the constructors against `MIN_PRECISION..=MAX_PRECISION` (the value itself is owned by ZB-TD02 on its branch;
//! this file checks only the refusals, ZB-TD03).
//! COVERS: Ball::exact, Ball::all_reals, Ball::from_mid_rad, BallArray::zeros
//! (Precision, MIN_PRECISION and MAX_PRECISION are a type alias and consts, which pubfn_coverage does not list.)

use zenith_bound::array::BallArray;
use zenith_bound::ball::Ball;
use zenith_bound::error::BoundError;
use zenith_bound::mag::Mag;
use zenith_bound::MAX_PRECISION;
use zenith_float_num::ExactNum;

fn invalid<T>(r: Result<T, BoundError>, what: &str) -> Result<(), String> {
    match r {
        Err(BoundError::InvalidArgument { .. }) => Ok(()),
        Err(e) => Err(format!("{}: want InvalidArgument, got {:?}", what, e)),
        Ok(_) => Err(format!("{}: want InvalidArgument, got Ok", what)),
    }
}

#[test]
fn constructors_reject_out_of_range_precision() -> Result<(), String> {
    for p in [0, MAX_PRECISION + 1, usize::MAX] {
        invalid(
            Ball::exact(ExactNum::from(3i64), p),
            &format!("exact {}", p),
        )?;
        invalid(Ball::all_reals(p), &format!("all_reals {}", p))?;
        invalid(
            Ball::from_mid_rad(ExactNum::from(3i64), Mag::zero(), p),
            &format!("from_mid_rad {}", p),
        )?;
        invalid(
            BallArray::zeros(vec![2].into_boxed_slice(), p),
            &format!("zeros {}", p),
        )?;
    }
    Ok(())
}
