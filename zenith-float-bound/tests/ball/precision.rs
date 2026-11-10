//! API-freeze unit 5 (Paul/Rusty 2026-10-23): every public prec is `Precision` (= zf's `usize` alias), validated
//! at the constructors against `MIN_PRECISION..=MAX_PRECISION` (`min(2^60, usize::MAX/4) - 1`, zf w6 defs.rs).
//! COVERS: Ball::exact, Ball::all_reals, Ball::from_mid_rad, BallArray::zeros
//! (Precision, MIN_PRECISION and MAX_PRECISION are a type alias and consts, which pubfn_coverage does not list.)

use zenith_bound::array::BallArray;
use zenith_bound::ball::Ball;
use zenith_bound::error::BoundError;
use zenith_bound::mag::Mag;
use zenith_bound::{Precision, MAX_PRECISION, MIN_PRECISION};
use zenith_float_num::ExactNum;

fn invalid<T>(r: Result<T, BoundError>, what: &str) -> Result<(), String> {
    match r {
        Err(BoundError::InvalidArgument { .. }) => Ok(()),
        Err(e) => Err(format!("{}: want InvalidArgument, got {:?}", what, e)),
        Ok(_) => Err(format!("{}: want InvalidArgument, got Ok", what)),
    }
}

#[test]
fn precision_limits_match_zf() -> Result<(), String> {
    let one: Precision = 1;
    if MIN_PRECISION != one {
        return Err(format!("MIN_PRECISION {}", MIN_PRECISION));
    }
    let want = if usize::BITS >= 64 {
        (1usize << 60) - 1
    } else {
        usize::MAX / 4 - 1
    };
    if MAX_PRECISION != want {
        return Err(format!("MAX_PRECISION {} want {}", MAX_PRECISION, want));
    }
    Ok(())
}

#[test]
fn constructors_accept_precisions_past_u32() -> Result<(), String> {
    let big: Precision = if usize::BITS >= 64 {
        1usize << 40
    } else {
        1usize << 28
    };
    for p in [1, 53, 16384, 65536, big, MAX_PRECISION] {
        let b =
            Ball::exact(ExactNum::from(3i64), p).map_err(|e| format!("exact {}: {:?}", p, e))?;
        if b.prec() != p {
            return Err(format!("prec {} read back as {}", p, b.prec()));
        }
        let a = Ball::all_reals(p).map_err(|e| format!("all_reals {}: {:?}", p, e))?;
        let r = Ball::from_mid_rad(ExactNum::from(1i64), Mag::zero(), p)
            .map_err(|e| format!("from_mid_rad {}: {:?}", p, e))?;
        if a.prec() != p || r.prec() != p || r.into_parts().2 != p {
            return Err(format!("prec {} lost", p));
        }
    }
    Ok(())
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
