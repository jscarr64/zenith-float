//! TEST_PLAN WR-ND (buildable part, plan addendum): rank 0..6 shapes, empty dims on every
//! axis, shape-product overflow => InvalidArgument. Oracle: row-major offsets computed by hand.
//! COVERS: BallArray::from_balls, BallArray::zeros, BallArray::get, BallArray::rank, BallArray::len

use zenith_bound::array::BallArray;
use zenith_bound::ball::Ball;
use zenith_bound::error::BoundError;
use zenith_float_num::ExactNum;

const DIMS: [usize; 6] = [2, 3, 1, 4, 2, 3];

fn balls(n: usize) -> Result<Vec<Ball>, String> {
    let mut v = Vec::new();
    for k in 0..n {
        let k = i64::try_from(k).map_err(|e| e.to_string())?;
        v.push(Ball::exact(ExactNum::from(k), 64).map_err(|e| format!("{:?}", e))?);
    }
    Ok(v)
}

/// The first `rank` entries of DIMS, with `zero_axis` (if any) set to 0. Built fresh per call.
fn dims(rank: usize, zero_axis: Option<usize>) -> Vec<usize> {
    DIMS.iter()
        .take(rank)
        .enumerate()
        .map(|(i, &d)| if Some(i) == zero_axis { 0 } else { d })
        .collect()
}

/// usize::MAX on `axis`, 2 elsewhere: the product overflows. Built fresh per call.
fn overflow_shape(rank: usize, axis: usize) -> Box<[usize]> {
    (0..rank)
        .map(|i| if i == axis { usize::MAX } else { 2 })
        .collect()
}

/// Every index of `shape` in row-major order.
fn indices(shape: &[usize]) -> Vec<Vec<usize>> {
    let mut out: Vec<Vec<usize>> = vec![Vec::new()];
    for &d in shape {
        let mut next = Vec::new();
        for prefix in &out {
            for i in 0..d {
                let mut ix: Vec<usize> = Vec::with_capacity(prefix.len() + 1);
                ix.extend(prefix);
                ix.push(i);
                next.push(ix);
            }
        }
        out = next;
    }
    out
}

#[test]
fn ranks_zero_to_six_are_row_major() -> Result<(), String> {
    for rank in 0..=6 {
        let shape: Vec<usize> = dims(rank, None);
        let n: usize = shape.iter().product();
        let a = BallArray::from_balls(dims(rank, None).into_boxed_slice(), balls(n)?)
            .map_err(|e| format!("rank {}: {:?}", rank, e))?;
        if a.rank() != rank || a.len() != n {
            return Err(format!("rank {}: rank {} len {}", rank, a.rank(), a.len()));
        }
        for (k, ix) in indices(&shape).iter().enumerate() {
            let want = ExactNum::from(i64::try_from(k).map_err(|e| e.to_string())?);
            match a.get(ix) {
                Some(b) if *b.mid() == want => {}
                other => return Err(format!("rank {} index {:?}: got {:?}", rank, ix, other)),
            }
        }
        let mut long: Vec<usize> = dims(rank, None);
        long.push(0);
        if a.get(&long).is_some() || shape.get(1..).is_some_and(|tail| a.get(tail).is_some()) {
            return Err(format!("rank {}: wrong-rank index accepted", rank));
        }
        let z = BallArray::zeros(shape.into_boxed_slice(), 64)
            .map_err(|e| format!("zeros rank {}: {:?}", rank, e))?;
        if z.len() != n || z.rank() != rank {
            return Err(format!("zeros rank {}: len {}", rank, z.len()));
        }
    }
    Ok(())
}

#[test]
fn an_empty_dim_on_any_axis_empties_the_array() -> Result<(), String> {
    for rank in 1..=6 {
        for axis in 0..rank {
            let shape: Vec<usize> = dims(rank, Some(axis));
            let a = BallArray::from_balls(dims(rank, Some(axis)).into_boxed_slice(), balls(0)?)
                .map_err(|e| format!("{:?}: {:?}", shape, e))?;
            let zero_ix: Vec<usize> = vec![0; rank];
            if !a.is_empty() || a.rank() != rank || a.get(&zero_ix).is_some() {
                return Err(format!("{:?}: not an empty rank-{} array", shape, rank));
            }
            if BallArray::from_balls(dims(rank, Some(axis)).into_boxed_slice(), balls(1)?).is_ok() {
                return Err(format!("{:?}: accepted one ball", shape));
            }
            let z = BallArray::zeros(dims(rank, Some(axis)).into_boxed_slice(), 64)
                .map_err(|e| format!("zeros {:?}: {:?}", shape, e))?;
            if !z.is_empty() {
                return Err(format!("zeros {:?} not empty", shape));
            }
        }
    }
    Ok(())
}

#[test]
fn shape_overflow_is_invalid_argument_at_every_rank_and_axis() -> Result<(), String> {
    for rank in 2..=6 {
        for axis in 0..rank {
            let shape = overflow_shape(rank, axis);
            for (what, r) in [
                (
                    "from_balls",
                    BallArray::from_balls(overflow_shape(rank, axis), balls(0)?).map(|_| ()),
                ),
                (
                    "zeros",
                    BallArray::zeros(overflow_shape(rank, axis), 64).map(|_| ()),
                ),
            ] {
                match r {
                    Err(BoundError::InvalidArgument { .. }) => {}
                    other => return Err(format!("{} {:?}: got {:?}", what, shape, other)),
                }
            }
        }
    }
    Ok(())
}

#[test]
fn a_zero_dim_anywhere_makes_an_empty_array_even_beside_huge_dims() -> Result<(), String> {
    // WR-ND: the element count is 0 whenever any extent is 0, whatever the order; an
    // overflowing partial product must not turn a valid empty shape into InvalidArgument.
    for shape in [
        [0, usize::MAX, 2],
        [usize::MAX, 0, 2],
        [usize::MAX, 2, 0],
        [usize::MAX, usize::MAX, 0],
    ] {
        match BallArray::from_balls(Box::new(shape), balls(0)?) {
            Ok(a) if a.is_empty() && a.rank() == 3 => {}
            other => {
                return Err(format!(
                    "from_balls {:?}: {:?}",
                    shape,
                    other.map(|a| a.len())
                ))
            }
        }
        match BallArray::zeros(Box::new(shape), 64) {
            Ok(a) if a.is_empty() => {}
            other => return Err(format!("zeros {:?}: {:?}", shape, other.map(|a| a.len()))),
        }
    }
    Ok(())
}
