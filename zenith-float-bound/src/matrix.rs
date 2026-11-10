//! Slice S6: Ball array and matrix operations on `BallArray`, `BallVec` and `BallMat` (the containers stay in `array`): matmul, norms, elementwise ops and N-dimensional enclosures.
//!
//! Owner: slice S6 (Cassie, assigned by Spence 2026-10-25). Only that slice edits this file.
//! CONTRACT 16.S6.0 to 16.S6.11. Elementwise ops and `from_mid_rad` reuse S1 (`Ball::add`, `Ball::sub`, `Ball::mul_2exp`).
//! Products (Spence ruling b, 2026-10-25) do their own rounding on zenith-float ops and their own radius accounting with S1-style
//! `Mag` primitives: sequential accumulation at `w = prec + ceil(log2(K + 1)) + 2`, one final rounding to `prec`, and a rounding
//! term for every inexact step (16.S6.7, 16.S6.9).

use core::cmp::Ordering;
use core::ops::Neg;

use zenith_float_num::{Error, ExactNum, RoundingMode, EXPONENT_MAX, EXPONENT_MIN};

use crate::array::{BallArray, BallMat, BallVec};
use crate::ball::Ball;
use crate::error::{BoundError, BoundResult};
use crate::mag::Mag;
use crate::precision::{Precision, MAX_PRECISION, MIN_PRECISION};

/// Every midpoint rounding is to nearest, ties to even (CONTRACT 15.1, O19), as in S1.
const MID_RM: RoundingMode = RoundingMode::ToEven;
/// Guard bits on top of `ceil(log2(K + 1))` (16.S6.7).
const GUARD_EXTRA: Precision = 2;

/// `prec` outside `MIN_PRECISION..=MAX_PRECISION` is `InvalidArgument` (16.S6.3), checked first in every op.
pub(crate) fn check_prec(op: &'static str, prec: Precision) -> BoundResult<()> {
    if (MIN_PRECISION..=MAX_PRECISION).contains(&prec) {
        Ok(())
    } else {
        Err(BoundError::InvalidArgument {
            op,
            detail: "prec outside MIN_PRECISION..=MAX_PRECISION",
        })
    }
}

pub(crate) fn mismatch(op: &'static str, detail: &'static str) -> BoundError {
    BoundError::DimensionMismatch { op, detail }
}

pub(crate) fn resource(op: &'static str, detail: &'static str) -> BoundError {
    BoundError::Resource { op, detail }
}

/// An empty vector with room for `n` elements, or `Resource` (16.S6.3, clause 15.5).
pub(crate) fn reserve<T>(op: &'static str, n: usize) -> BoundResult<Vec<T>> {
    let mut v = Vec::new();
    v.try_reserve_exact(n)
        .map_err(|_e| resource(op, "result storage cannot be reserved"))?;
    Ok(v)
}

/// The product of `shape`, `None` on `usize` overflow.
fn element_count(shape: &[usize]) -> Option<usize> {
    shape.iter().try_fold(1usize, |acc, d| acc.checked_mul(*d))
}

/// `w = prec + ceil(log2(K + 1)) + 2`, capped at `MAX_PRECISION` (16.S6.7). `ceil(log2(K + 1))` is the bit length of `K`.
fn work_prec(prec: Precision, k: usize) -> Precision {
    let bits = usize::BITS.checked_sub(k.leading_zeros());
    let guard = bits
        .and_then(|b| Precision::try_from(b).ok())
        .and_then(|b| b.checked_add(GUARD_EXTRA));
    match guard.and_then(|g| prec.checked_add(g)) {
        Some(w) => w.min(MAX_PRECISION),
        None => MAX_PRECISION,
    }
}

/// A zenith-float result rounded at some precision: the value and its rounding term, or an overflow (clause 10.4).
enum Step {
    Value(ExactNum, Mag),
    Overflow,
}

/// Classifies a rounded zenith-float result as S1 does (clause 10.4, 16.S1.2): ±Inf, an exponent-overflow NaN and an inexact
/// top-binade value are overflow; the rounding term is zero when exact, `half_ulp(e, prec)` for exponent `e > EXPONENT_MIN`, and
/// `2^EXPONENT_MIN` for a zero or bottom-binade result.
fn rounded(op: &'static str, raw: ExactNum, prec: Precision) -> BoundResult<Step> {
    if raw.is_inf() {
        return Ok(Step::Overflow);
    }
    if raw.is_nan() {
        return match raw.err() {
            Some(Error::ExponentOverflow(_)) => Ok(Step::Overflow),
            Some(Error::MemoryAllocation) => Err(resource(op, "zf memory allocation")),
            Some(Error::PrecisionRetryExhausted) => {
                Err(resource(op, "zf precision retry exhausted"))
            }
            // Mapped as S1 does (`arith::classify`, CONTRACT O2).
            Some(Error::DivisionByZero) => Err(BoundError::Domain {
                op,
                detail: "zf division by zero",
            }),
            Some(Error::InvalidArgument) => Err(BoundError::InvalidArgument {
                op,
                detail: "zf invalid argument",
            }),
            None => Err(BoundError::Representation {
                op,
                detail: "zf NaN without an error kind",
            }),
        };
    }
    if raw.inexact() && raw.exponent() == Some(EXPONENT_MAX) {
        return Ok(Step::Overflow);
    }
    if !raw.inexact() {
        return Ok(Step::Value(raw, Mag::zero()));
    }
    let normal = raw
        .exponent()
        .filter(|e| *e > EXPONENT_MIN && !raw.is_zero());
    let eps = match normal {
        Some(e) => Mag::half_ulp(e, prec),
        None => Mag::ulp_mag(EXPONENT_MIN, 0),
    };
    Ok(Step::Value(raw, eps))
}

/// One product entry `Σ_k a_k·b_k` (16.S6.7): each item is `(m(a), r(a), m(b), r(b))`, `k` the number of items.
fn product_entry<'a, I>(op: &'static str, terms: I, k: usize, prec: Precision) -> BoundResult<Ball>
where
    I: Iterator<Item = (&'a ExactNum, &'a Mag, &'a ExactNum, &'a Mag)>,
{
    entry_core(op, None, terms, false, k, prec)
}

/// The 16.S6.7 accumulation shared by `product_entry` and the LA helpers: `base ± Σ a_k·b_k`, every step rounded at
/// `w = work_prec(prec, k)` with its rounding term added to the radius, then rounded once to `prec`. `base` (when given) is
/// read by borrow and counts as one of the `k` items. `subtract` selects `base − Σ` (16.LA3.2 reflector updates).
fn entry_core<'a, I>(
    op: &'static str,
    base: Option<&Ball>,
    terms: I,
    subtract: bool,
    k: usize,
    prec: Precision,
) -> BoundResult<Ball>
where
    I: Iterator<Item = (&'a ExactNum, &'a Mag, &'a ExactNum, &'a Mag)>,
{
    let w = work_prec(prec, k);
    let mut s = ExactNum::from(0i64);
    let mut rad = Mag::zero();
    if let Some(b) = base {
        rad = rad.add_up(b.rad());
        if rad.is_infinity() {
            return Ball::all_reals(prec);
        }
        s = match rounded(op, ExactNum::from(0i64).add(b.mid(), w, MID_RM), w)? {
            Step::Value(v, e) => {
                rad = rad.add_up(&e);
                v
            }
            Step::Overflow => return Ball::all_reals(prec),
        };
    }
    for (ma, ra, mb, rb) in terms {
        // R1 term: |m(a)|·r(b) + r(a)·(|m(b)| + r(b))
        let mb_up = Mag::from_exact_num(mb);
        let first = Mag::from_exact_num(ma)
            .mul_up(rb)
            .add_up(&ra.mul_up(&mb_up.add_up(rb)));
        rad = rad.add_up(&first);
        if rad.is_infinity() {
            return Ball::all_reals(prec);
        }
        let p = match rounded(op, ma.mul_full_prec(mb).ldexp(0, w, MID_RM), w)? {
            Step::Value(v, e) => {
                rad = rad.add_up(&e);
                v
            }
            Step::Overflow => return Ball::all_reals(prec),
        };
        let next = if subtract {
            s.sub(&p, w, MID_RM)
        } else {
            s.add(&p, w, MID_RM)
        };
        s = match rounded(op, next, w)? {
            Step::Value(v, e) => {
                rad = rad.add_up(&e);
                v
            }
            Step::Overflow => return Ball::all_reals(prec),
        };
    }
    match rounded(op, s.ldexp(0, prec, MID_RM), prec)? {
        Step::Value(m, e) => Ball::from_mid_rad(m, rad.add_up(&e), prec),
        Step::Overflow => Ball::all_reals(prec),
    }
}

/// The parts of a pair of borrowed balls, as `entry_core` reads them.
fn pair_parts<'a>((a, b): (&'a Ball, &'a Ball)) -> (&'a ExactNum, &'a Mag, &'a ExactNum, &'a Mag) {
    (a.mid(), a.rad(), b.mid(), b.rad())
}

/// `Σ a_k·b_k` over borrowed ball pairs as one 16.S6.7 product entry at `prec`; `k` is the number of pairs (LA helper).
pub(crate) fn dot_pairs<'a, I>(
    op: &'static str,
    pairs: I,
    k: usize,
    prec: Precision,
) -> BoundResult<Ball>
where
    I: Iterator<Item = (&'a Ball, &'a Ball)>,
{
    entry_core(op, None, pairs.map(pair_parts), false, k, prec)
}

/// `base − Σ a_k·b_k` over borrowed balls as one 16.S6.7 entry at `prec`; `k` counts the pairs (LA helper, 16.LA3.2).
pub(crate) fn sub_pairs<'a, I>(
    op: &'static str,
    base: &Ball,
    pairs: I,
    k: usize,
    prec: Precision,
) -> BoundResult<Ball>
where
    I: Iterator<Item = (&'a Ball, &'a Ball)>,
{
    let items = k
        .checked_add(1)
        .ok_or_else(|| resource(op, "term count overflows usize"))?;
    entry_core(op, Some(base), pairs.map(pair_parts), true, items, prec)
}

/// `b` rounded outward to `prec` once (16.LA.6): the midpoint to nearest at `prec`, its rounding term added to the radius.
pub(crate) fn round_out(op: &'static str, b: Ball, prec: Precision) -> BoundResult<Ball> {
    let (mid, rad, _) = b.into_parts();
    if rad.is_infinity() {
        return Ball::all_reals(prec);
    }
    match rounded(op, mid.ldexp(0, prec, MID_RM), prec)? {
        Step::Value(m, e) => Ball::from_mid_rad(m, rad.add_up(&e), prec),
        Step::Overflow => Ball::all_reals(prec),
    }
}

/// True when the ball certainly excludes 0: `r < |m|` with `|m|` rounded down (an `Infinity` radius never excludes 0).
pub(crate) fn excludes_zero(b: &Ball) -> bool {
    !b.rad().is_infinity() && Mag::from_exact_num_lower(b.mid()) > *b.rad()
}

/// True when the ball is exactly `0 ± 0` (the certified-refusal test of 16.LA.2).
pub(crate) fn is_exact_zero(b: &Ball) -> bool {
    b.rad().is_zero() && b.mid().is_zero()
}

/// 16.LA.6: `w = prec + 64 + 2·⌈log2(n + 1)⌉`, capped at `MAX_PRECISION`; `⌈log2(n + 1)⌉` is the bit length of `n`.
pub(crate) fn la_work_prec(prec: Precision, n: usize) -> Precision {
    let bits = usize::BITS.checked_sub(n.leading_zeros());
    let guard = bits
        .and_then(|b| Precision::try_from(b).ok())
        .and_then(|b| b.checked_mul(2))
        .and_then(|b| b.checked_add(64));
    match guard.and_then(|g| prec.checked_add(g)) {
        Some(w) => w.min(MAX_PRECISION),
        None => MAX_PRECISION,
    }
}

/// Row `i` of a row-major `? x k` matrix.
fn row_of<'a>(op: &'static str, data: &'a [Ball], i: usize, k: usize) -> BoundResult<&'a [Ball]> {
    let start = i
        .checked_mul(k)
        .ok_or_else(|| resource(op, "row offset overflows usize"))?;
    let end = start
        .checked_add(k)
        .ok_or_else(|| resource(op, "row offset overflows usize"))?;
    data.get(start..end)
        .ok_or_else(|| resource(op, "row out of range"))
}

/// `|m| + r` rounded up (16.S6.6).
fn entry_up(b: &Ball) -> Mag {
    Mag::from_exact_num(b.mid()).add_up(b.rad())
}

/// `Σ (|m| + r)` rounded up.
fn sum_up<'a, I: Iterator<Item = &'a Ball>>(it: I) -> Mag {
    it.fold(Mag::zero(), |acc, b| acc.add_up(&entry_up(b)))
}

/// `Σ (|m| + r)^2` rounded up.
fn sum_sq_up<'a, I: Iterator<Item = &'a Ball>>(it: I) -> Mag {
    it.fold(Mag::zero(), |acc, b| {
        let a = entry_up(b);
        acc.add_up(&a.mul_up(&a))
    })
}

/// A rank-N array of magnitudes, row-major (CONTRACT 16.S6.0, 16.S6.1).
#[derive(Debug)]
pub struct MagArray {
    shape: Box<[usize]>,
    data: Box<[Mag]>,
}

impl MagArray {
    /// Builds an array of the given shape from `mags` in row-major order.
    ///
    /// # Errors
    /// `InvalidArgument` when the product of `shape` overflows `usize` or differs from `mags.len()`.
    pub fn from_mags(shape: Box<[usize]>, mags: Vec<Mag>) -> BoundResult<MagArray> {
        if element_count(&shape) != Some(mags.len()) {
            return Err(BoundError::InvalidArgument {
                op: "MagArray::from_mags",
                detail: "shape product overflows or differs from the number of magnitudes",
            });
        }
        Ok(MagArray {
            shape,
            data: mags.into_boxed_slice(),
        })
    }

    /// The shape.
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

    /// Whether there are no elements.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// The element at a full row-major index; `None` for a wrong rank or an index out of range.
    pub fn get(&self, index: &[usize]) -> Option<&Mag> {
        if index.len() != self.shape.len() {
            return None;
        }
        let mut offset = 0usize;
        for (i, d) in index.iter().zip(self.shape.iter()) {
            if i >= d {
                return None;
            }
            offset = offset.checked_mul(*d)?.checked_add(*i)?;
        }
        self.data.get(offset)
    }

    /// The elements in row-major order.
    pub fn as_slice(&self) -> &[Mag] {
        &self.data
    }

    /// The shape and the elements.
    pub fn into_parts(self) -> (Box<[usize]>, Box<[Mag]>) {
        (self.shape, self.data)
    }
}

impl BallArray {
    /// Elementwise `add`: entry `i` is `Ball::add` of entry `i` (CONTRACT 16.S6.4).
    ///
    /// # Errors
    /// `InvalidArgument` for a bad `prec`, `DimensionMismatch` for unequal shapes, `Resource` on allocation failure, and any
    /// error of `Ball::add`.
    pub fn add(self, other: &BallArray, prec: Precision) -> BoundResult<BallArray> {
        const OP: &str = "BallArray::add";
        check_prec(OP, prec)?;
        if self.shape() != other.shape() {
            return Err(mismatch(OP, "shapes differ"));
        }
        let (shape, data) = self.into_parts();
        let mut out = reserve(OP, data.len())?;
        for (x, y) in data.into_vec().into_iter().zip(other.as_slice()) {
            out.push(x.add(y, prec)?);
        }
        BallArray::from_balls(shape, out)
    }

    /// Elementwise `sub`: entry `i` is `Ball::sub` of entry `i` (CONTRACT 16.S6.4).
    ///
    /// # Errors
    /// As [`BallArray::add`].
    pub fn sub(self, other: &BallArray, prec: Precision) -> BoundResult<BallArray> {
        const OP: &str = "BallArray::sub";
        check_prec(OP, prec)?;
        if self.shape() != other.shape() {
            return Err(mismatch(OP, "shapes differ"));
        }
        let (shape, data) = self.into_parts();
        let mut out = reserve(OP, data.len())?;
        for (x, y) in data.into_vec().into_iter().zip(other.as_slice()) {
            out.push(x.sub(y, prec)?);
        }
        BallArray::from_balls(shape, out)
    }

    /// Splits into exact midpoints `m_i ± 0` at each entry's own precision and radii `r_i` (CONTRACT 16.S6.5; Spence ruling a,
    /// 2026-10-25: fallible, no infallible internal constructor).
    ///
    /// # Errors
    /// `Resource` on allocation failure, and any error of `Ball::exact`.
    pub fn mid_rad(self) -> BoundResult<(BallArray, MagArray)> {
        const OP: &str = "BallArray::mid_rad";
        let (shape, data) = self.into_parts();
        let mut mids = reserve(OP, data.len())?;
        let mut rads = reserve(OP, data.len())?;
        for b in data.into_vec() {
            let (m, r, p) = b.into_parts();
            mids.push(Ball::exact(m, p)?);
            rads.push(r);
        }
        let rad_shape = Box::<[usize]>::from(&*shape);
        Ok((
            BallArray::from_balls(shape, mids)?,
            MagArray::from_mags(rad_shape, rads)?,
        ))
    }

    /// Rebuilds balls `round(m_i) ± (s_i + r_i + ε)` from midpoints and extra radii, enlarging only (CONTRACT 16.S6.5).
    ///
    /// # Errors
    /// `InvalidArgument` for a bad `prec`, `DimensionMismatch` for unequal shapes, `Resource` on allocation failure.
    pub fn from_mid_rad(
        mids: BallArray,
        rads: &MagArray,
        prec: Precision,
    ) -> BoundResult<BallArray> {
        const OP: &str = "BallArray::from_mid_rad";
        check_prec(OP, prec)?;
        if mids.shape() != rads.shape() {
            return Err(mismatch(OP, "shapes differ"));
        }
        let (shape, data) = mids.into_parts();
        let mut out = reserve(OP, data.len())?;
        for (b, r) in data.into_vec().into_iter().zip(rads.as_slice()) {
            let (m, s, _) = b.into_parts();
            out.push(Ball::from_mid_rad(m, s.add_up(r), prec)?.mul_2exp(0, prec)?);
        }
        BallArray::from_balls(shape, out)
    }
}

impl BallVec {
    /// Elementwise `add`, forwarding to [`BallArray::add`] (CONTRACT 16.S6.4).
    ///
    /// # Errors
    /// As [`BallArray::add`].
    pub fn add(self, other: &BallVec, prec: Precision) -> BoundResult<BallVec> {
        BallVec::from_array(self.into_array().add(other.as_array(), prec)?)
    }

    /// Elementwise `sub`, forwarding to [`BallArray::sub`] (CONTRACT 16.S6.4).
    ///
    /// # Errors
    /// As [`BallArray::add`].
    pub fn sub(self, other: &BallVec, prec: Precision) -> BoundResult<BallVec> {
        BallVec::from_array(self.into_array().sub(other.as_array(), prec)?)
    }

    /// Dot product (CONTRACT 16.S6.7 to 16.S6.9).
    ///
    /// # Errors
    /// `InvalidArgument` for a bad `prec`, `DimensionMismatch` for unequal lengths.
    pub fn dot(&self, other: &BallVec, prec: Precision) -> BoundResult<Ball> {
        const OP: &str = "BallVec::dot";
        check_prec(OP, prec)?;
        let k = self.len();
        if other.len() != k {
            return Err(mismatch(OP, "lengths differ"));
        }
        let x = self.as_array().as_slice();
        let y = other.as_array().as_slice();
        let terms = x
            .iter()
            .zip(y)
            .map(|(a, b)| (a.mid(), a.rad(), b.mid(), b.rad()));
        product_entry(OP, terms, k, prec)
    }

    /// Upper bound of the largest `|m| + r` over the entries (CONTRACT 16.S6.6, BallVec norms); 0 for the empty vector.
    pub fn norm_inf(&self) -> Mag {
        self.as_array()
            .as_slice()
            .iter()
            .map(entry_up)
            .fold(Mag::zero(), Ord::max)
    }

    /// Upper bound of the sum of `|m| + r` over the entries (CONTRACT 16.S6.6, BallVec norms); 0 for the empty vector.
    pub fn norm_1(&self) -> Mag {
        sum_up(self.as_array().as_slice().iter())
    }

    /// Upper bound of the Euclidean norm of every vector in the ball, `sqrt(sum (|m| + r)^2)` rounded up (CONTRACT 16.S6.6);
    /// 0 for the empty vector.
    pub fn norm_2_up(&self) -> Mag {
        sum_sq_up(self.as_array().as_slice().iter()).sqrt_up()
    }
}

impl BallMat {
    /// Elementwise `add`, forwarding to [`BallArray::add`] (CONTRACT 16.S6.4).
    ///
    /// # Errors
    /// As [`BallArray::add`].
    pub fn add(self, other: &BallMat, prec: Precision) -> BoundResult<BallMat> {
        BallMat::from_array(self.into_array().add(other.as_array(), prec)?)
    }

    /// Elementwise `sub`, forwarding to [`BallArray::sub`] (CONTRACT 16.S6.4).
    ///
    /// # Errors
    /// As [`BallArray::add`].
    pub fn sub(self, other: &BallMat, prec: Precision) -> BoundResult<BallMat> {
        BallMat::from_array(self.into_array().sub(other.as_array(), prec)?)
    }

    /// Matrix-vector product (CONTRACT 16.S6.7 to 16.S6.9).
    ///
    /// # Errors
    /// `InvalidArgument` for a bad `prec`, `DimensionMismatch` when `ncols != x.len()`, `Resource` on allocation failure.
    pub fn matvec(&self, x: &BallVec, prec: Precision) -> BoundResult<BallVec> {
        const OP: &str = "BallMat::matvec";
        check_prec(OP, prec)?;
        let (m, k) = (self.nrows(), self.ncols());
        if x.len() != k {
            return Err(mismatch(OP, "ncols != x.len()"));
        }
        let a = self.as_array().as_slice();
        let v = x.as_array().as_slice();
        let mut out = reserve(OP, m)?;
        for i in 0..m {
            let row = row_of(OP, a, i, k)?;
            let terms = row
                .iter()
                .zip(v)
                .map(|(p, q)| (p.mid(), p.rad(), q.mid(), q.rad()));
            out.push(product_entry(OP, terms, k, prec)?);
        }
        Ok(BallVec::from_balls(out))
    }

    /// Matrix product (CONTRACT 16.S6.7 to 16.S6.9).
    ///
    /// # Errors
    /// `InvalidArgument` for a bad `prec`, `DimensionMismatch` when `self.ncols() != other.nrows()`, `Resource` on allocation failure.
    pub fn matmul(&self, other: &BallMat, prec: Precision) -> BoundResult<BallMat> {
        const OP: &str = "BallMat::matmul";
        check_prec(OP, prec)?;
        let (m, k, n) = (self.nrows(), self.ncols(), other.ncols());
        if other.nrows() != k {
            return Err(mismatch(OP, "self.ncols() != other.nrows()"));
        }
        let count = m
            .checked_mul(n)
            .ok_or_else(|| resource(OP, "result size overflows usize"))?;
        let a = self.as_array().as_slice();
        let b = other.as_array().as_slice();
        let mut out = reserve(OP, count)?;
        for i in 0..m {
            let row = row_of(OP, a, i, k)?;
            for j in 0..n {
                let col = b.iter().skip(j).step_by(n);
                let terms = row
                    .iter()
                    .zip(col)
                    .map(|(p, q)| (p.mid(), p.rad(), q.mid(), q.rad()));
                out.push(product_entry(OP, terms, k, prec)?);
            }
        }
        BallMat::from_balls(m, n, out)
    }

    /// Upper bound of the max row sum of `|m| + r` (CONTRACT 16.S6.6).
    pub fn norm_inf(&self) -> Mag {
        let n = self.ncols();
        if n == 0 {
            return Mag::zero();
        }
        self.as_array()
            .as_slice()
            .chunks_exact(n)
            .map(|row| sum_up(row.iter()))
            .fold(Mag::zero(), Ord::max)
    }

    /// Upper bound of the max column sum of `|m| + r` (CONTRACT 16.S6.6).
    pub fn norm_1(&self) -> Mag {
        let n = self.ncols();
        if n == 0 {
            return Mag::zero();
        }
        let data = self.as_array().as_slice();
        (0..n)
            .map(|j| sum_up(data.iter().skip(j).step_by(n)))
            .fold(Mag::zero(), Ord::max)
    }

    /// Upper bound of the Frobenius norm of every matrix in the ball (CONTRACT 16.S6.6).
    pub fn norm_frob_up(&self) -> Mag {
        sum_sq_up(self.as_array().as_slice().iter()).sqrt_up()
    }
}

/// The triangle a 16.S6.12 solve reads.
enum Triangle {
    Lower,
    Upper,
}

/// "Contains 0" as the exact test `|m| ≤ r` (16.S6.12.2); an `Infinity` radius contains 0.
fn contains_zero(b: &Ball) -> bool {
    match b.rad().cmp_exact_num(b.mid()) {
        Some(Ordering::Less) => false,
        Some(Ordering::Equal | Ordering::Greater) | None => true,
    }
}

/// The 16.S6.12.2 shape checks: `self` square, then `b` with `n` entries (rows). Returns `n`.
fn tri_shape(op: &'static str, t: &BallMat, b_rows: usize) -> BoundResult<usize> {
    let n = t.nrows();
    if t.ncols() != n {
        return Err(mismatch(op, "rectangular"));
    }
    if b_rows != n {
        return Err(mismatch(op, "b_len"));
    }
    Ok(n)
}

/// The diagonal entry `t_ii` of a row-major `n x n` matrix.
fn diag_of<'a>(op: &'static str, t: &'a [Ball], n: usize, i: usize) -> BoundResult<&'a Ball> {
    row_of(op, t, i, n)?
        .get(i)
        .ok_or_else(|| resource(op, "row out of range"))
}

/// The certified refusal of 16.S6.12.2: with `unit_diag = false`, the first diagonal ball (lowest index) that contains 0 is
/// `Singular "singular"`.
fn tri_diag_check(op: &'static str, t: &[Ball], n: usize, unit_diag: bool) -> BoundResult<()> {
    if unit_diag {
        return Ok(());
    }
    for i in 0..n {
        if contains_zero(diag_of(op, t, n, i)?) {
            return Err(BoundError::Singular {
                op,
                detail: "singular",
            });
        }
    }
    Ok(())
}

/// `Σ_j t_ij·x_j − b_i = −u_i` as one 16.S6.7 product entry at `w` over `k` terms: `(b_i, −1)` then the `(t_ij, x_j)` pairs.
/// This is the 16.S6.12.3 entry `(1, b_i)`, `(−t_ij, x_j)` with every sign flipped: each R1 term and each round-to-nearest-even
/// step is the mirror image, so the ball is exactly `−u_i`.
fn tri_neg_numerator<'a, I>(
    op: &'static str,
    b_i: &'a Ball,
    off: I,
    k: usize,
    w: Precision,
) -> BoundResult<Ball>
where
    I: Iterator<Item = (&'a Ball, &'a Ball)>,
{
    let neg_one = ExactNum::from(-1i64);
    let zero = Mag::zero();
    let terms = core::iter::once((b_i.mid(), b_i.rad(), &neg_one, &zero))
        .chain(off.map(|(t, x)| (t.mid(), t.rad(), x.mid(), x.rad())));
    product_entry(op, terms, k, w)
}

/// One validated triangular system of 16.S6.12: `t` is the row-major `n x n` matrix.
struct TriSys<'t> {
    op: &'static str,
    t: &'t [Ball],
    n: usize,
    side: Triangle,
    unit_diag: bool,
    prec: Precision,
}

/// Solves one column (16.S6.12.3): forward substitution for `Lower`, the mirror with `i` descending for `Upper`. `xw` holds
/// the solutions at the working precision `w` in solve order; each is rounded outward to `prec` once into `out`, in index order.
fn tri_column<'b, F>(
    sys: &TriSys<'_>,
    b_at: F,
    xw: &mut Vec<Ball>,
    out: &mut Vec<Ball>,
) -> BoundResult<()>
where
    F: Fn(usize) -> Option<&'b Ball>,
{
    let TriSys {
        op,
        t,
        n,
        ref side,
        unit_diag,
        prec,
    } = *sys;
    let range = || resource(op, "index out of range");
    let w = la_work_prec(prec, n);
    xw.clear();
    for s in 0..n {
        let i = match side {
            Triangle::Lower => s,
            Triangle::Upper => n
                .checked_sub(s)
                .and_then(|v| v.checked_sub(1))
                .ok_or_else(range)?,
        };
        let row = row_of(op, t, i, n)?;
        let b_i = b_at(i).ok_or_else(range)?;
        let k = s.checked_add(1).ok_or_else(range)?;
        let neg_u = match side {
            Triangle::Lower => {
                let off = row.get(..i).ok_or_else(range)?;
                tri_neg_numerator(op, b_i, off.iter().zip(xw.iter()), k, w)?
            }
            Triangle::Upper => {
                let start = i.checked_add(1).ok_or_else(range)?;
                let off = row.get(start..).ok_or_else(range)?;
                tri_neg_numerator(op, b_i, off.iter().rev().zip(xw.iter()), k, w)?
            }
        };
        let u = Neg::neg(neg_u);
        let x = if unit_diag {
            u
        } else {
            u.div(row.get(i).ok_or_else(range)?, w)?
        };
        xw.push(x);
    }
    match side {
        Triangle::Lower => {
            for x in xw.drain(..) {
                out.push(x.mul_2exp(0, prec)?);
            }
        }
        Triangle::Upper => {
            for x in xw.drain(..).rev() {
                out.push(x.mul_2exp(0, prec)?);
            }
        }
    }
    Ok(())
}

/// The vector form of 16.S6.12 in the error order of 16.S6.12.2.
fn tri_solve_vec(
    op: &'static str,
    t: &BallMat,
    b: &BallVec,
    side: Triangle,
    unit_diag: bool,
    prec: Precision,
) -> BoundResult<BallVec> {
    check_prec(op, prec)?;
    let n = tri_shape(op, t, b.len())?;
    let mut xw = reserve(op, n)?;
    let mut out = reserve(op, n)?;
    let a = t.as_array().as_slice();
    tri_diag_check(op, a, n, unit_diag)?;
    let bs = b.as_array().as_slice();
    let sys = TriSys {
        op,
        t: a,
        n,
        side,
        unit_diag,
        prec,
    };
    tri_column(&sys, |i| bs.get(i), &mut xw, &mut out)?;
    Ok(BallVec::from_balls(out))
}

/// The many-RHS form of 16.S6.12: column `j` of the result is the vector form on column `j` of `b`.
fn tri_solve_mat(
    op: &'static str,
    t: &BallMat,
    b: &BallMat,
    side: Triangle,
    unit_diag: bool,
    prec: Precision,
) -> BoundResult<BallMat> {
    check_prec(op, prec)?;
    let n = tri_shape(op, t, b.nrows())?;
    let k = b.ncols();
    let count = n
        .checked_mul(k)
        .ok_or_else(|| resource(op, "result size overflows usize"))?;
    let mut xw = reserve(op, n)?;
    let mut cols: Vec<Vec<Ball>> = reserve(op, k)?;
    for _ in 0..k {
        cols.push(reserve(op, n)?);
    }
    let mut iters: Vec<std::vec::IntoIter<Ball>> = reserve(op, k)?;
    let mut out = reserve(op, count)?;
    let a = t.as_array().as_slice();
    tri_diag_check(op, a, n, unit_diag)?;
    let bs = b.as_array().as_slice();
    let sys = TriSys {
        op,
        t: a,
        n,
        side,
        unit_diag,
        prec,
    };
    for (j, col) in cols.iter_mut().enumerate() {
        let b_at = |i: usize| {
            i.checked_mul(k)
                .and_then(|p| p.checked_add(j))
                .and_then(|p| bs.get(p))
        };
        tri_column(&sys, b_at, &mut xw, col)?;
    }
    for col in cols {
        iters.push(col.into_iter());
    }
    for _ in 0..n {
        for it in &mut iters {
            out.push(
                it.next()
                    .ok_or_else(|| resource(op, "column shorter than n"))?,
            );
        }
    }
    BallMat::from_balls(n, k, out)
}

impl BallMat {
    /// Forward substitution: `T⁻¹b` for the lower triangle of `self` (CONTRACT 16.S6.12). Only the entries on and below the
    /// diagonal are read; with `unit_diag` the diagonal is taken as exactly 1 and not read. Works at the 16.LA.6 precision and
    /// rounds each entry outward to `prec` once.
    ///
    /// # Errors
    /// In order (16.S6.12.2): `InvalidArgument` for a bad `prec`; `DimensionMismatch` `"rectangular"` for a non-square `self`,
    /// `"b_len"` when `b.len() != n`; `Resource` on allocation failure; with `unit_diag = false`, `Singular` `"singular"` for
    /// the first diagonal ball that contains 0. An `Infinity` term or overflow gives all-reals entries, `Ok`.
    pub fn solve_lower(
        &self,
        b: &BallVec,
        unit_diag: bool,
        prec: Precision,
    ) -> BoundResult<BallVec> {
        tri_solve_vec(
            "BallMat::solve_lower",
            self,
            b,
            Triangle::Lower,
            unit_diag,
            prec,
        )
    }

    /// Backward substitution: `T⁻¹b` for the upper triangle of `self` (CONTRACT 16.S6.12), the mirror of
    /// [`Self::solve_lower`] with `i` descending. Only the entries on and above the diagonal are read.
    ///
    /// # Errors
    /// As [`Self::solve_lower`].
    pub fn solve_upper(
        &self,
        b: &BallVec,
        unit_diag: bool,
        prec: Precision,
    ) -> BoundResult<BallVec> {
        tri_solve_vec(
            "BallMat::solve_upper",
            self,
            b,
            Triangle::Upper,
            unit_diag,
            prec,
        )
    }

    /// Many-RHS forward substitution (CONTRACT 16.S6.12.1): column `j` of the result is [`Self::solve_lower`] on column `j`
    /// of `b`.
    ///
    /// # Errors
    /// As [`Self::solve_lower`], with `"b_len"` when `b.nrows() != n`.
    pub fn solve_lower_mat(
        &self,
        b: &BallMat,
        unit_diag: bool,
        prec: Precision,
    ) -> BoundResult<BallMat> {
        tri_solve_mat(
            "BallMat::solve_lower_mat",
            self,
            b,
            Triangle::Lower,
            unit_diag,
            prec,
        )
    }

    /// Many-RHS backward substitution (CONTRACT 16.S6.12.1): column `j` of the result is [`Self::solve_upper`] on column `j`
    /// of `b`.
    ///
    /// # Errors
    /// As [`Self::solve_lower`], with `"b_len"` when `b.nrows() != n`.
    pub fn solve_upper_mat(
        &self,
        b: &BallMat,
        unit_diag: bool,
        prec: Precision,
    ) -> BoundResult<BallMat> {
        tri_solve_mat(
            "BallMat::solve_upper_mat",
            self,
            b,
            Triangle::Upper,
            unit_diag,
            prec,
        )
    }
}

/// zf NaN error kinds from a product step map as in S1 (`arith::classify`, CONTRACT O2): `DivisionByZero` is `Domain`,
/// `InvalidArgument` is `InvalidArgument`, a NaN without an error kind is `Representation`.
#[cfg(test)]
mod slice_matrix_nan_map {
    use super::{rounded, Step};
    use crate::error::BoundError;
    use zenith_float_num::{Error, ExactNum};

    fn step_error(err: Option<Error>) -> Result<BoundError, String> {
        match rounded("nan_map", ExactNum::nan(err), 53) {
            Err(e) => Ok(e),
            Ok(Step::Value(..)) | Ok(Step::Overflow) => {
                Err(String::from("a NaN step did not error"))
            }
        }
    }

    #[test]
    fn slice_matrix_zf_nan_errors_map_as_s1() -> Result<(), String> {
        let e = step_error(Some(Error::DivisionByZero))?;
        if !(matches!(e, BoundError::Domain { .. }) && e.detail() == "zf division by zero") {
            return Err(format!("DivisionByZero mapped to {e:?}"));
        }
        let e = step_error(Some(Error::InvalidArgument))?;
        if !(matches!(e, BoundError::InvalidArgument { .. }) && e.detail() == "zf invalid argument")
        {
            return Err(format!("InvalidArgument mapped to {e:?}"));
        }
        let e = step_error(None)?;
        if !(matches!(e, BoundError::Representation { .. })
            && e.detail() == "zf NaN without an error kind")
        {
            return Err(format!("kindless NaN mapped to {e:?}"));
        }
        Ok(())
    }
}

/// The pub(crate) LA helpers (16.LA.6, 16.S6.7 split for 16.LA3.2), each checked on exact inputs.
#[cfg(test)]
mod la_helpers {
    use super::{dot_pairs, excludes_zero, is_exact_zero, la_work_prec, round_out, sub_pairs};
    use crate::ball::Ball;
    use crate::mag::Mag;
    use crate::precision::MAX_PRECISION;
    use zenith_float_num::ExactNum;

    fn pt(v: i64) -> Result<Ball, String> {
        Ball::exact(ExactNum::from(v), 64).map_err(|e| format!("{e:?}"))
    }

    fn is_point(b: &Ball, v: i64) -> bool {
        b.rad().is_zero() && b.mid().sub_full_prec(&ExactNum::from(v)).is_zero()
    }

    #[test]
    fn la_work_prec_follows_16_la_6() -> Result<(), String> {
        // n = 0: 0 bits; n = 1: 1; n = 3: 2; n = 4: 3; n = 12: 4.
        for (n, want) in [(0usize, 117usize), (1, 119), (3, 121), (4, 123), (12, 125)] {
            if la_work_prec(53, n) != want {
                return Err(format!("w(53, {n}) = {}, want {want}", la_work_prec(53, n)));
            }
        }
        if la_work_prec(MAX_PRECISION, 5) != MAX_PRECISION {
            return Err("w is not capped at MAX_PRECISION".to_string());
        }
        Ok(())
    }

    #[test]
    fn dot_pairs_is_exact_on_small_integers() -> Result<(), String> {
        let (a, b, c, d) = (pt(3)?, pt(4)?, pt(-2)?, pt(5)?);
        let s = dot_pairs("t", [(&a, &b), (&c, &d)].into_iter(), 2, 53)
            .map_err(|e| format!("{e:?}"))?;
        if !is_point(&s, 2) {
            return Err(format!("3·4 + (−2)·5 = {s:?}"));
        }
        Ok(())
    }

    #[test]
    fn sub_pairs_subtracts_from_the_base_and_keeps_its_radius() -> Result<(), String> {
        let base = Ball::from_mid_rad(
            ExactNum::from(10i64),
            Mag::from_exact_num(&ExactNum::from(1i64)),
            64,
        )
        .map_err(|e| format!("{e:?}"))?;
        let (a, b) = (pt(2)?, pt(3)?);
        let s =
            sub_pairs("t", &base, [(&a, &b)].into_iter(), 1, 53).map_err(|e| format!("{e:?}"))?;
        if !s.mid().sub_full_prec(&ExactNum::from(4i64)).is_zero() || *s.rad() < *base.rad() {
            return Err(format!("10±1 − 2·3 = {s:?}"));
        }
        Ok(())
    }

    #[test]
    fn round_out_encloses_and_keeps_exact_points() -> Result<(), String> {
        let third = pt(1)?.div(&pt(3)?, 200).map_err(|e| format!("{e:?}"))?;
        let r = round_out("t", third, 53).map_err(|e| format!("{e:?}"))?;
        // r must cover 1/3: |3·mid − 1| ≤ 3·rad.
        let err = r
            .mid()
            .mul_full_prec(&ExactNum::from(3i64))
            .sub_full_prec(&ExactNum::from(1i64));
        if Mag::from_exact_num(&err) > Mag::from_exact_num(&ExactNum::from(3i64)).mul_up(r.rad())
            || r.rad().is_zero()
        {
            return Err(format!("round_out(1/3) = {r:?}"));
        }
        let seven = round_out("t", pt(7)?, 53).map_err(|e| format!("{e:?}"))?;
        if !is_point(&seven, 7) {
            return Err(format!("round_out(7) = {seven:?}"));
        }
        Ok(())
    }

    #[test]
    fn zero_tests_separate_certified_and_undecided() -> Result<(), String> {
        let one = Mag::from_exact_num(&ExactNum::from(1i64));
        let straddle =
            Ball::from_mid_rad(ExactNum::from(1i64), one, 64).map_err(|e| format!("{e:?}"))?;
        let zero = pt(0)?;
        let away = pt(-2)?;
        let all = Ball::all_reals(64).map_err(|e| format!("{e:?}"))?;
        let ok = !excludes_zero(&straddle)
            && !excludes_zero(&zero)
            && excludes_zero(&away)
            && !excludes_zero(&all)
            && is_exact_zero(&zero)
            && !is_exact_zero(&straddle)
            && !is_exact_zero(&all);
        if ok {
            Ok(())
        } else {
            Err("excludes_zero / is_exact_zero misclassify".to_string())
        }
    }
}
