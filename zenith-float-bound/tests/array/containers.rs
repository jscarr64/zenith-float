//! PLAN C5 arrays (rulings c/e), plan addendum L2c-1: BallArray rank-N row-major, private fields,
//! checked constructors; BallVec/BallMat thin newtypes. Oracle: shapes and indices built by hand.
//! COVERS: BallArray::from_balls, BallArray::zeros, BallArray::shape, BallArray::rank, BallArray::len, BallArray::is_empty, BallArray::get, BallArray::as_slice, BallArray::into_parts, BallVec::from_balls, BallVec::from_array, BallVec::len, BallVec::is_empty, BallVec::get, BallVec::as_array, BallVec::into_array, BallMat::from_balls, BallMat::from_array, BallMat::nrows, BallMat::ncols, BallMat::get, BallMat::as_array, BallMat::into_array

use zenith_bound::array::{BallArray, BallMat, BallVec};
use zenith_bound::ball::Ball;
use zenith_bound::error::BoundError;
use zenith_float_num::ExactNum;

fn balls(n: i64) -> Result<Vec<Ball>, String> {
    (0..n)
        .map(|v| Ball::exact(ExactNum::from(v), 64).map_err(|e| format!("{:?}", e)))
        .collect()
}

fn is(b: Option<&Ball>, v: i64) -> bool {
    b.is_some_and(|b| *b.mid() == ExactNum::from(v))
}

fn invalid<T: core::fmt::Debug>(r: Result<T, BoundError>) -> bool {
    matches!(r, Err(BoundError::InvalidArgument { .. }))
}

#[test]
fn from_balls_is_row_major() -> Result<(), String> {
    let a =
        BallArray::from_balls(Box::new([2, 3, 4]), balls(24)?).map_err(|e| format!("{:?}", e))?;
    if a.shape() != [2, 3, 4] || a.rank() != 3 || a.len() != 24 || a.is_empty() {
        return Err(format!(
            "shape {:?} rank {} len {}",
            a.shape(),
            a.rank(),
            a.len()
        ));
    }
    for i in 0usize..2 {
        for j in 0..3 {
            for k in 0..4 {
                let want = i64::try_from(i * 12 + j * 4 + k)
                    .map_err(|e: core::num::TryFromIntError| e.to_string())?;
                if !is(a.get(&[i, j, k]), want) {
                    return Err(format!("get({},{},{}) != {}", i, j, k, want));
                }
            }
        }
    }
    if !is(a.as_slice().get(23), 23) {
        return Err("as_slice is not row-major".into());
    }
    Ok(())
}

#[test]
fn get_rejects_bad_indices() -> Result<(), String> {
    let a = BallArray::from_balls(Box::new([2, 3]), balls(6)?).map_err(|e| format!("{:?}", e))?;
    let bad: [&[usize]; 5] = [&[2, 0], &[0, 3], &[0], &[0, 0, 0], &[usize::MAX, 0]];
    for idx in bad {
        if a.get(idx).is_some() {
            return Err(format!("get({:?}) should be None", idx));
        }
    }
    if !is(a.get(&[1, 2]), 5) {
        return Err("get(1,2) != 5".into());
    }
    Ok(())
}

#[test]
fn from_balls_rejects_shape_mismatch_and_overflow() -> Result<(), String> {
    if !invalid(BallArray::from_balls(Box::new([2, 3]), balls(5)?)) {
        return Err("2x3 from 5 balls accepted".into());
    }
    if !invalid(BallArray::from_balls(Box::new([usize::MAX, 2]), balls(0)?)) {
        return Err("overflowing shape accepted".into());
    }
    if !invalid(BallArray::zeros(Box::new([usize::MAX, usize::MAX]), 64)) {
        return Err("overflowing zeros accepted".into());
    }
    Ok(())
}

#[test]
fn empty_and_rank_zero_shapes() -> Result<(), String> {
    let e = BallArray::from_balls(Box::new([3, 0]), balls(0)?).map_err(|e| format!("{:?}", e))?;
    if !e.is_empty() || e.get(&[0, 0]).is_some() {
        return Err("3x0 array is not empty".into());
    }
    let s = BallArray::from_balls(Box::new([]), balls(1)?).map_err(|e| format!("{:?}", e))?;
    if s.rank() != 0 || s.len() != 1 || !is(s.get(&[]), 0) {
        return Err("rank-0 array must hold one element".into());
    }
    Ok(())
}

#[test]
fn zeros_are_exact_zero_at_prec() -> Result<(), String> {
    let z = BallArray::zeros(Box::new([2, 2]), 128).map_err(|e| format!("{:?}", e))?;
    if z.len() != 4 {
        return Err(format!("len {}", z.len()));
    }
    for b in z.as_slice() {
        if !b.is_exact() || b.prec() != 128 || *b.mid() != ExactNum::from(0i64) {
            return Err(format!("not an exact zero at 128: {:?}", b));
        }
    }
    if !invalid(BallArray::zeros(Box::new([1]), 0)) {
        return Err("zeros at prec 0 accepted".into());
    }
    Ok(())
}

#[test]
fn into_parts_moves_shape_and_data() -> Result<(), String> {
    let a = BallArray::from_balls(Box::new([3]), balls(3)?).map_err(|e| format!("{:?}", e))?;
    let (shape, data) = a.into_parts();
    if *shape != [3] || data.len() != 3 || !is(data.get(2), 2) {
        return Err(format!("into_parts gave {:?} / {}", shape, data.len()));
    }
    Ok(())
}

#[test]
fn ball_vec_is_rank_one() -> Result<(), String> {
    let v = BallVec::from_balls(balls(4)?);
    if v.len() != 4 || v.is_empty() || !is(v.get(3), 3) || v.get(4).is_some() {
        return Err("BallVec accessors".into());
    }
    if v.as_array().rank() != 1 {
        return Err("BallVec rank".into());
    }
    let a = v.into_array();
    let back = BallVec::from_array(a).map_err(|e| format!("{:?}", e))?;
    if back.len() != 4 {
        return Err("round trip".into());
    }
    let m = BallArray::from_balls(Box::new([2, 2]), balls(4)?).map_err(|e| format!("{:?}", e))?;
    if !invalid(BallVec::from_array(m)) {
        return Err("rank-2 array accepted as BallVec".into());
    }
    if !BallVec::from_balls(Vec::new()).is_empty() {
        return Err("empty BallVec".into());
    }
    Ok(())
}

#[test]
fn ball_mat_is_rank_two_row_major() -> Result<(), String> {
    let m = BallMat::from_balls(2, 3, balls(6)?).map_err(|e| format!("{:?}", e))?;
    if m.nrows() != 2 || m.ncols() != 3 || !is(m.get(1, 0), 3) || m.get(0, 3).is_some() {
        return Err("BallMat accessors".into());
    }
    if m.as_array().shape() != [2, 3] {
        return Err("BallMat shape".into());
    }
    let back = BallMat::from_array(m.into_array()).map_err(|e| format!("{:?}", e))?;
    if back.ncols() != 3 {
        return Err("round trip".into());
    }
    if !invalid(BallMat::from_balls(2, 2, balls(3)?)) {
        return Err("2x2 from 3 accepted".into());
    }
    if !invalid(BallMat::from_array(
        BallVec::from_balls(balls(2)?).into_array(),
    )) {
        return Err("rank-1 array accepted as BallMat".into());
    }
    Ok(())
}

#[test]
fn from_balls_rejects_more_balls_than_the_shape() -> Result<(), String> {
    // L0b mutant 148-A-fromlen2: a surplus ball is a mismatch too.
    if !invalid(BallArray::from_balls(Box::new([2, 3]), balls(7)?)) {
        return Err("7 balls for shape [2, 3] must be InvalidArgument".into());
    }
    Ok(())
}

#[test]
fn zeros_too_large_to_allocate_is_resource_not_abort() -> Result<(), String> {
    // Gap F7 (API freeze): allocation failure is BoundError::Resource, never an abort. The
    // shape product fits in usize, so this is not InvalidArgument; the byte size does not.
    for shape in [[1usize << 60], [usize::MAX / 2]] {
        match BallArray::zeros(Box::new(shape), 64) {
            Err(BoundError::Resource {
                op: "BallArray::zeros",
                ..
            }) => {}
            other => {
                return Err(format!(
                    "zeros({:?}): want Resource, got {:?}",
                    shape, other
                ))
            }
        }
    }
    Ok(())
}
