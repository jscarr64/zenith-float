//! Ball containers (PLAN C5 arrays, rulings c/e; plan addendum L2c-1).
//!
//! [`BallArray`] is a dynamic rank-N array of [`Ball`]s in row-major order. [`BallVec`] and
//! [`BallMat`] are thin rank-1 and rank-2 newtypes over it. Fields are private; there is no
//! `Clone` and no `PartialEq` (an enclosure has no unambiguous `==`).
//!
//! Reserved for later units of the 2.0 API (names fixed now, types added with their ops):
//! `ComplexBallArray`, `ComplexBallVec`, `ComplexBallMat` (after `ComplexBall`, L2a),
//! `MagArray` (L2c-2), `BallLu` (L2c-5), `BallCholesky` (L2c-6), `BallQr` (L2c-8),
//! `GlsResult` (generalised least squares). Elementwise ops, `identical`, `contains` and
//! `overlaps` arrive with the matching `Ball` ops.

use zenith_float_num::ExactNum;

use crate::ball::Ball;
use crate::error::{BoundError, BoundResult};
use crate::precision::Precision;

/// Precision of the exact zero midpoint used by [`BallArray::zeros`].
const ZERO_PREC: usize = 64;

/// A rank-N array of balls, row-major.
#[derive(Debug)]
pub struct BallArray {
    shape: Box<[usize]>,
    data: Box<[Ball]>,
}

impl BallArray {
    /// Builds an array of the given shape from `balls` in row-major order.
    ///
    /// # Errors
    /// `InvalidArgument` when the product of `shape` overflows `usize` (a shape with a zero
    /// extent never overflows: it has 0 elements) or differs from `balls.len()`.
    pub fn from_balls(shape: Box<[usize]>, balls: Vec<Ball>) -> BoundResult<BallArray> {
        let len = element_count("BallArray::from_balls", &shape)?;
        if len != balls.len() {
            return Err(BoundError::InvalidArgument {
                op: "BallArray::from_balls",
                detail: "shape does not match the number of balls",
            });
        }
        Ok(BallArray {
            shape,
            data: balls.into_boxed_slice(),
        })
    }

    /// An array of exact zeros (`0 ± 0`) at `prec`.
    ///
    /// # Errors
    /// `InvalidArgument` when `prec == 0` or the shape's product overflows; `Resource` when
    /// the storage cannot be allocated.
    pub fn zeros(shape: Box<[usize]>, prec: Precision) -> BoundResult<BallArray> {
        let len = element_count("BallArray::zeros", &shape)?;
        let mut balls = Vec::new();
        balls
            .try_reserve_exact(len)
            .map_err(|_| BoundError::Resource {
                op: "BallArray::zeros",
                detail: "allocation failed",
            })?;
        for _ in 0..len {
            balls.push(Ball::exact(ExactNum::new(ZERO_PREC), prec)?);
        }
        BallArray::from_balls(shape, balls)
    }

    /// The extent of each axis.
    pub fn shape(&self) -> &[usize] {
        &self.shape
    }

    /// The number of axes.
    pub fn rank(&self) -> usize {
        self.shape.len()
    }

    /// The number of elements.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// True when the array has no elements.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// The element at `index` (one entry per axis), or `None` when out of range.
    pub fn get(&self, index: &[usize]) -> Option<&Ball> {
        if index.len() != self.shape.len() {
            return None;
        }
        let mut offset: usize = 0;
        for (&i, &dim) in index.iter().zip(self.shape.iter()) {
            if i >= dim {
                return None;
            }
            offset = offset.checked_mul(dim)?.checked_add(i)?;
        }
        self.data.get(offset)
    }

    /// The elements in row-major order.
    pub fn as_slice(&self) -> &[Ball] {
        &self.data
    }

    /// Moves the shape and the row-major elements out.
    pub fn into_parts(self) -> (Box<[usize]>, Box<[Ball]>) {
        (self.shape, self.data)
    }
}

/// A rank-1 [`BallArray`].
#[derive(Debug)]
pub struct BallVec(BallArray);

impl BallVec {
    /// A vector of the given balls.
    pub fn from_balls(balls: Vec<Ball>) -> BallVec {
        let shape: Box<[usize]> = Box::new([balls.len()]);
        BallVec(BallArray {
            shape,
            data: balls.into_boxed_slice(),
        })
    }

    /// Checked conversion from an array.
    ///
    /// # Errors
    /// `InvalidArgument` when the array's rank is not 1.
    pub fn from_array(a: BallArray) -> BoundResult<BallVec> {
        match a.rank() {
            1 => Ok(BallVec(a)),
            _ => Err(BoundError::InvalidArgument {
                op: "BallVec::from_array",
                detail: "rank must be 1",
            }),
        }
    }

    /// The number of elements.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// True when the vector has no elements.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Element `i`, or `None` when out of range.
    pub fn get(&self, i: usize) -> Option<&Ball> {
        self.0.get(&[i])
    }

    /// The underlying array.
    pub fn as_array(&self) -> &BallArray {
        &self.0
    }

    /// Moves the underlying array out.
    pub fn into_array(self) -> BallArray {
        self.0
    }
}

/// A rank-2 [`BallArray`], row-major.
#[derive(Debug)]
pub struct BallMat(BallArray);

impl BallMat {
    /// An `nrows × ncols` matrix from row-major `balls`.
    ///
    /// # Errors
    /// `InvalidArgument` when `nrows * ncols` overflows or differs from `balls.len()`.
    pub fn from_balls(nrows: usize, ncols: usize, balls: Vec<Ball>) -> BoundResult<BallMat> {
        BallArray::from_balls(Box::new([nrows, ncols]), balls).map(BallMat)
    }

    /// Checked conversion from an array.
    ///
    /// # Errors
    /// `InvalidArgument` when the array's rank is not 2.
    pub fn from_array(a: BallArray) -> BoundResult<BallMat> {
        match a.rank() {
            2 => Ok(BallMat(a)),
            _ => Err(BoundError::InvalidArgument {
                op: "BallMat::from_array",
                detail: "rank must be 2",
            }),
        }
    }

    /// The number of rows.
    pub fn nrows(&self) -> usize {
        self.0.shape().first().map_or(0, |r| *r)
    }

    /// The number of columns.
    pub fn ncols(&self) -> usize {
        self.0.shape().get(1).map_or(0, |c| *c)
    }

    /// Element `(i, j)`, or `None` when out of range.
    pub fn get(&self, i: usize, j: usize) -> Option<&Ball> {
        self.0.get(&[i, j])
    }

    /// The underlying array.
    pub fn as_array(&self) -> &BallArray {
        &self.0
    }

    /// Moves the underlying array out.
    pub fn into_array(self) -> BallArray {
        self.0
    }
}

/// Product of the axis extents, checked; 0 when any extent is 0.
fn element_count(op: &'static str, shape: &[usize]) -> BoundResult<usize> {
    // Any zero extent empties the array, whatever the other extents (WR-ND).
    if shape.contains(&0) {
        return Ok(0);
    }
    shape
        .iter()
        .try_fold(1usize, |acc, &d| acc.checked_mul(d))
        .ok_or(BoundError::InvalidArgument {
            op,
            detail: "shape product overflows usize",
        })
}
