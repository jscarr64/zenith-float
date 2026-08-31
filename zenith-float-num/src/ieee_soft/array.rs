//! Dense row-major arrays of software IEEE values and `ExactNum`.
//! A 1-D vector is stored as shape `(1, n)`.

use super::simd::{add_u32_lanes, add_u64_lanes, mul_u32_lanes, mul_u64_lanes};
use super::{Ieee32, Ieee64};
use crate::defs::RoundingMode;
use crate::Consts;
use crate::ExactNum;
use alloc::vec::Vec;

trait LaneBits: Copy {
    fn add_lanes(a: &[Self], b: &[Self]) -> Vec<Self>;
    fn mul_lanes(a: &[Self], b: &[Self]) -> Vec<Self>;
}

impl LaneBits for u32 {
    fn add_lanes(a: &[Self], b: &[Self]) -> Vec<Self> {
        add_u32_lanes(a, b)
    }
    fn mul_lanes(a: &[Self], b: &[Self]) -> Vec<Self> {
        mul_u32_lanes(a, b)
    }
}

impl LaneBits for u64 {
    fn add_lanes(a: &[Self], b: &[Self]) -> Vec<Self> {
        add_u64_lanes(a, b)
    }
    fn mul_lanes(a: &[Self], b: &[Self]) -> Vec<Self> {
        mul_u64_lanes(a, b)
    }
}

/// Contiguous binary32 lanes (`u32` bits), row-major.
#[derive(Clone, Debug)]
pub struct Ieee32Array {
    bits: Vec<u32>,
    rows: usize,
    cols: usize,
}

/// Contiguous binary64 lanes (`u64` bits), row-major.
#[derive(Clone, Debug)]
pub struct Ieee64Array {
    bits: Vec<u64>,
    rows: usize,
    cols: usize,
}

/// Batch of `ExactNum` values at a shared precision `p`, row-major.
#[derive(Clone, Debug)]
pub struct ExactNumArray {
    p: usize,
    vals: Vec<ExactNum>,
    rows: usize,
    cols: usize,
}

macro_rules! impl_ieee_array {
    ($arr:ident, $scalar:ident, $bits:ty) => {
        impl $arr {
            /// Empty 0×0 array.
            pub fn new() -> Self {
                Self {
                    bits: Vec::new(),
                    rows: 0,
                    cols: 0,
                }
            }

            /// Row vector: `n` copies of `fill` (shape `(1, n)`).
            pub fn filled(n: usize, fill: $scalar) -> Self {
                Self {
                    bits: alloc::vec![fill.to_bits(); n],
                    rows: 1,
                    cols: n,
                }
            }

            /// `rows×cols` filled with `fill`.
            pub fn filled_2d(rows: usize, cols: usize, fill: $scalar) -> Option<Self> {
                let n = rows.checked_mul(cols)?;
                Some(Self {
                    bits: alloc::vec![fill.to_bits(); n],
                    rows,
                    cols,
                })
            }

            /// From a slice of IEEE bit patterns as a row vector.
            pub fn from_bits(bits: &[$bits]) -> Self {
                Self {
                    bits: bits.to_vec(),
                    rows: 1,
                    cols: bits.len(),
                }
            }

            /// Row-major `rows×cols`. Length must be `rows*cols`.
            pub fn from_shape(rows: usize, cols: usize, vals: &[$scalar]) -> Option<Self> {
                let n = rows.checked_mul(cols)?;
                if n != vals.len() {
                    return None;
                }
                Some(Self {
                    bits: vals.iter().map(|v| v.to_bits()).collect(),
                    rows,
                    cols,
                })
            }

            /// From scalars as a row vector (shape `(1, n)`).
            pub fn from_values(vals: &[$scalar]) -> Self {
                Self {
                    bits: vals.iter().map(|v| v.to_bits()).collect(),
                    rows: 1,
                    cols: vals.len(),
                }
            }

            /// `(rows, cols)`.
            pub fn shape(&self) -> (usize, usize) {
                (self.rows, self.cols)
            }

            /// Reinterpret the same buffer as `rows×cols` when the product matches.
            pub fn reshape(&self, rows: usize, cols: usize) -> Option<Self> {
                let n = rows.checked_mul(cols)?;
                if n != self.bits.len() {
                    return None;
                }
                Some(Self {
                    bits: self.bits.clone(),
                    rows,
                    cols,
                })
            }

            /// Number of lanes.
            pub fn len(&self) -> usize {
                self.bits.len()
            }

            /// True if there are no lanes.
            pub fn is_empty(&self) -> bool {
                self.bits.is_empty()
            }

            /// IEEE bit patterns.
            pub fn as_bits(&self) -> &[$bits] {
                &self.bits
            }

            /// Lane `i` in storage order, or `None` if out of range.
            pub fn get(&self, i: usize) -> Option<$scalar> {
                self.bits.get(i).copied().map($scalar::from_bits)
            }

            /// Entry `(i, j)`, or `None` if out of range.
            pub fn get2(&self, i: usize, j: usize) -> Option<$scalar> {
                if i >= self.rows || j >= self.cols {
                    return None;
                }
                self.get(i * self.cols + j)
            }

            /// Elementwise add. Shapes must match. Uses integer SIMD when the
            /// architecture provides it (still the software IEEE kernel).
            pub fn add(&self, rhs: &Self) -> Option<Self> {
                if self.rows != rhs.rows || self.cols != rhs.cols {
                    return None;
                }
                Some(Self {
                    bits: <$bits>::add_lanes(&self.bits, &rhs.bits),
                    rows: self.rows,
                    cols: self.cols,
                })
            }

            /// Add a scalar to every lane.
            pub fn add_scalar(&self, s: $scalar) -> Self {
                self.map(|x| x.add(s))
            }

            /// Elementwise sub.
            pub fn sub(&self, rhs: &Self) -> Option<Self> {
                self.zip_op(rhs, $scalar::sub)
            }

            /// Elementwise mul. Shapes must match. Integer SIMD significand
            /// products on the all-normal path.
            pub fn mul(&self, rhs: &Self) -> Option<Self> {
                if self.rows != rhs.rows || self.cols != rhs.cols {
                    return None;
                }
                Some(Self {
                    bits: <$bits>::mul_lanes(&self.bits, &rhs.bits),
                    rows: self.rows,
                    cols: self.cols,
                })
            }

            /// Multiply every lane by a scalar.
            pub fn mul_scalar(&self, s: $scalar) -> Self {
                self.map(|x| x.mul(s))
            }

            /// Elementwise div.
            pub fn div(&self, rhs: &Self) -> Option<Self> {
                self.zip_op(rhs, $scalar::div)
            }

            /// Elementwise sqrt.
            pub fn sqrt(&self) -> Self {
                self.map($scalar::sqrt)
            }

            /// Sequential IEEE sum (one rounding per add).
            pub fn sum(&self) -> $scalar {
                let mut acc = $scalar::ZERO;
                for &b in &self.bits {
                    acc = acc.add($scalar::from_bits(b));
                }
                acc
            }

            /// Sequential dot product (mul then add, two rounds per term).
            pub fn dot(&self, rhs: &Self) -> Option<$scalar> {
                if self.len() != rhs.len() {
                    return None;
                }
                let mut acc = $scalar::ZERO;
                for (a, b) in self.bits.iter().zip(rhs.bits.iter()) {
                    acc = acc.add($scalar::from_bits(*a).mul($scalar::from_bits(*b)));
                }
                Some(acc)
            }

            /// Software matmul: `(m×k)(k×n) → (m×n)`. Sequential IEEE mul-then-add per term.
            pub fn matmul(&self, rhs: &Self) -> Option<Self> {
                if self.cols != rhs.rows {
                    return None;
                }
                let m = self.rows;
                let k = self.cols;
                let n = rhs.cols;
                let mut bits = alloc::vec![$scalar::ZERO.to_bits(); m.checked_mul(n)?];
                for i in 0..m {
                    for j in 0..n {
                        let mut acc = $scalar::ZERO;
                        for t in 0..k {
                            let a = $scalar::from_bits(self.bits[i * k + t]);
                            let b = $scalar::from_bits(rhs.bits[t * n + j]);
                            acc = acc.add(a.mul(b));
                        }
                        bits[i * n + j] = acc.to_bits();
                    }
                }
                Some(Self {
                    bits,
                    rows: m,
                    cols: n,
                })
            }

            fn zip_op(&self, rhs: &Self, op: fn($scalar, $scalar) -> $scalar) -> Option<Self> {
                if self.rows != rhs.rows || self.cols != rhs.cols {
                    return None;
                }
                Some(Self {
                    bits: self
                        .bits
                        .iter()
                        .zip(rhs.bits.iter())
                        .map(|(a, b)| op($scalar::from_bits(*a), $scalar::from_bits(*b)).to_bits())
                        .collect(),
                    rows: self.rows,
                    cols: self.cols,
                })
            }

            fn map(&self, op: impl Fn($scalar) -> $scalar) -> Self {
                Self {
                    bits: self
                        .bits
                        .iter()
                        .map(|&b| op($scalar::from_bits(b)).to_bits())
                        .collect(),
                    rows: self.rows,
                    cols: self.cols,
                }
            }

            /// Widen each lane, apply `op`, IEEE-round back.
            pub fn map_exact<F>(&self, p: usize, mut op: F) -> Self
            where
                F: FnMut(&ExactNum) -> ExactNum,
            {
                Self {
                    bits: self
                        .bits
                        .iter()
                        .map(|&b| {
                            let x = $scalar::from_bits(b).to_exact(p);
                            $scalar::from_exact(&op(&x)).to_bits()
                        })
                        .collect(),
                    rows: self.rows,
                    cols: self.cols,
                }
            }
        }

        impl Default for $arr {
            fn default() -> Self {
                Self::new()
            }
        }
    };
}

impl_ieee_array!(Ieee32Array, Ieee32, u32);
impl_ieee_array!(Ieee64Array, Ieee64, u64);

macro_rules! ieee_unary_exact {
    ($p:expr, $($name:ident),+ $(,)?) => {
        $(
            #[doc = concat!("Elementwise `", stringify!($name), "` via `ExactNum`.")]
            pub fn $name(&self, cc: &mut Consts) -> Self {
                self.map_exact($p, |x| x.$name($p, RoundingMode::ToEven, cc))
            }
        )+
    };
}

macro_rules! ieee_array_specials {
    ($arr:ident, $scalar:ident, $p:expr) => {
        impl $arr {
            ieee_unary_exact!(
                $p, exp, exp2, exp10, expm1, ln, log2, log10, log1p, sin, cos, tan, asin, acos,
                atan, sinh, cosh, tanh, asinh, acosh, atanh, erf, erfc, gamma, ln_gamma, digamma,
                ei, si, ci, li, fresnel_s, fresnel_c, ai, bi, elliptic_k, rem_pi,
            );

            /// Elementwise complete `E(m)` via `ExactNum`.
            pub fn elliptic_e_complete(&self, cc: &mut Consts) -> Self {
                self.map_exact($p, |x| x.elliptic_e_complete($p, RoundingMode::ToEven, cc))
            }

            /// Elementwise `cbrt` via `ExactNum`.
            pub fn cbrt(&self) -> Self {
                self.map_exact($p, |x| x.cbrt($p, RoundingMode::ToEven))
            }

            /// Elementwise `atan2(self, x)` via `ExactNum`.
            pub fn atan2(&self, x: $scalar, cc: &mut Consts) -> Self {
                let xe = x.to_exact($p);
                self.map_exact($p, |y| y.atan2(&xe, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise `hypot(self, other)` via `ExactNum`.
            pub fn hypot(&self, other: $scalar) -> Self {
                let oe = other.to_exact($p);
                self.map_exact($p, |x| x.hypot(&oe, $p, RoundingMode::ToEven))
            }

            /// Elementwise `pow(self, n)` via `ExactNum`.
            pub fn pow(&self, n: $scalar, cc: &mut Consts) -> Self {
                let ne = n.to_exact($p);
                self.map_exact($p, |x| x.pow(&ne, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise `log(self, base)` via `ExactNum`.
            pub fn log(&self, base: $scalar, cc: &mut Consts) -> Self {
                let be = base.to_exact($p);
                self.map_exact($p, |x| x.log(&be, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise `γ(self, x)` via `ExactNum`.
            pub fn gammainc(&self, x: $scalar, cc: &mut Consts) -> Self {
                let xe = x.to_exact($p);
                self.map_exact($p, |s| s.gammainc(&xe, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise `Γ(self, x)` via `ExactNum`.
            pub fn gammainc_upper(&self, x: $scalar, cc: &mut Consts) -> Self {
                let xe = x.to_exact($p);
                self.map_exact($p, |s| s.gammainc_upper(&xe, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise integer-order `J_n(self)`.
            pub fn bessel_j(&self, n: usize, cc: &mut Consts) -> Self {
                self.map_exact($p, |x| x.bessel_j(n, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise `J_ν(self)` for a shared real order.
            pub fn bessel_j_nu(&self, nu: $scalar, cc: &mut Consts) -> Self {
                let n = nu.to_exact($p);
                self.map_exact($p, |x| x.bessel_j_nu(&n, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise `Y_ν(self)`.
            pub fn bessel_y(&self, nu: $scalar, cc: &mut Consts) -> Self {
                let n = nu.to_exact($p);
                self.map_exact($p, |x| x.bessel_y(&n, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise `I_ν(self)`.
            pub fn bessel_i(&self, nu: $scalar, cc: &mut Consts) -> Self {
                let n = nu.to_exact($p);
                self.map_exact($p, |x| x.bessel_i(&n, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise `K_ν(self)`.
            pub fn bessel_k(&self, nu: $scalar, cc: &mut Consts) -> Self {
                let n = nu.to_exact($p);
                self.map_exact($p, |x| x.bessel_k(&n, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise `P_n(self)`.
            pub fn legendre_p(&self, n: u32) -> Self {
                self.map_exact($p, |x| x.legendre_p(n, $p, RoundingMode::ToEven))
            }

            /// Elementwise `P_n^m(self)`.
            pub fn assoc_legendre_p(&self, n: u32, m: i32) -> Self {
                self.map_exact($p, |x| x.assoc_legendre_p(n, m, $p, RoundingMode::ToEven))
            }

            /// Elementwise `F(self | m)`.
            pub fn elliptic_f(&self, m: $scalar, cc: &mut Consts) -> Self {
                let me = m.to_exact($p);
                self.map_exact($p, |x| x.elliptic_f(&me, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise incomplete `E(self | m)`.
            pub fn elliptic_e(&self, m: $scalar, cc: &mut Consts) -> Self {
                let me = m.to_exact($p);
                self.map_exact($p, |x| x.elliptic_e(&me, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise complete `Π(n, m)` with `self = n`.
            pub fn elliptic_pi_complete(&self, m: $scalar, cc: &mut Consts) -> Self {
                let me = m.to_exact($p);
                self.map_exact($p, |n| {
                    n.elliptic_pi_complete(&me, $p, RoundingMode::ToEven, cc)
                })
            }

            /// Elementwise `Π(self; x | m)`.
            pub fn elliptic_pi(&self, x: $scalar, m: $scalar, cc: &mut Consts) -> Self {
                let xe = x.to_exact($p);
                let me = m.to_exact($p);
                self.map_exact($p, |n| {
                    n.elliptic_pi(&xe, &me, $p, RoundingMode::ToEven, cc)
                })
            }

            /// Elementwise `{}_2F_1(self, b; c; z)`.
            pub fn hypergeom_2f1(
                &self,
                b: $scalar,
                c: $scalar,
                z: $scalar,
                cc: &mut Consts,
            ) -> Self {
                let be = b.to_exact($p);
                let ce = c.to_exact($p);
                let ze = z.to_exact($p);
                self.map_exact($p, |a| {
                    a.hypergeom_2f1(&be, &ce, &ze, $p, RoundingMode::ToEven, cc)
                })
            }

            /// Elementwise `I_x(self, b)`.
            pub fn betainc(&self, b: $scalar, x: $scalar, cc: &mut Consts) -> Self {
                let be = b.to_exact($p);
                let xe = x.to_exact($p);
                self.map_exact($p, |a| a.betainc(&be, &xe, $p, RoundingMode::ToEven, cc))
            }
        }
    };
}

ieee_array_specials!(Ieee32Array, Ieee32, 64);
ieee_array_specials!(Ieee64Array, Ieee64, 128);

macro_rules! exact_arr_p_rm_cc {
    ($($name:ident),+ $(,)?) => {
        $(
            #[doc = concat!("Elementwise [`ExactNum::", stringify!($name), "`]. `cc` is the constants cache, not a global.")]
            pub fn $name(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
                self.map_at(p, |x| x.$name(p, rm, cc))
            }
        )+
    };
}

impl ExactNumArray {
    /// Empty 0×0 array at precision `p`.
    pub fn new(p: usize) -> Self {
        Self {
            p,
            vals: Vec::new(),
            rows: 0,
            cols: 0,
        }
    }

    /// Row vector: `n` copies of `fill` (rounded to `p`).
    pub fn filled(p: usize, n: usize, fill: &ExactNum) -> Self {
        let mut v = fill.clone();
        let _ = v.set_precision(p, RoundingMode::ToEven);
        Self {
            p,
            vals: alloc::vec![v; n],
            rows: 1,
            cols: n,
        }
    }

    /// `rows×cols` filled with `fill` (rounded to `p`).
    pub fn filled_2d(p: usize, rows: usize, cols: usize, fill: &ExactNum) -> Option<Self> {
        let n = rows.checked_mul(cols)?;
        let mut v = fill.clone();
        let _ = v.set_precision(p, RoundingMode::ToEven);
        Some(Self {
            p,
            vals: alloc::vec![v; n],
            rows,
            cols,
        })
    }

    /// From values as a row vector; each is rounded to `p`.
    pub fn from_values(p: usize, vals: &[ExactNum]) -> Self {
        Self {
            p,
            vals: vals
                .iter()
                .map(|x| {
                    let mut y = x.clone();
                    let _ = y.set_precision(p, RoundingMode::ToEven);
                    y
                })
                .collect(),
            rows: 1,
            cols: vals.len(),
        }
    }

    /// Row-major `rows×cols`. Length must be `rows*cols`.
    pub fn from_shape(p: usize, rows: usize, cols: usize, vals: &[ExactNum]) -> Option<Self> {
        let n = rows.checked_mul(cols)?;
        if n != vals.len() {
            return None;
        }
        Some(Self {
            p,
            vals: vals
                .iter()
                .map(|x| {
                    let mut y = x.clone();
                    let _ = y.set_precision(p, RoundingMode::ToEven);
                    y
                })
                .collect(),
            rows,
            cols,
        })
    }

    /// Shared precision.
    pub fn precision(&self) -> usize {
        self.p
    }

    /// `(rows, cols)`.
    pub fn shape(&self) -> (usize, usize) {
        (self.rows, self.cols)
    }

    /// Reinterpret the same buffer as `rows×cols` when the product matches.
    pub fn reshape(&self, rows: usize, cols: usize) -> Option<Self> {
        let n = rows.checked_mul(cols)?;
        if n != self.vals.len() {
            return None;
        }
        Some(Self {
            p: self.p,
            vals: self.vals.clone(),
            rows,
            cols,
        })
    }

    /// Number of lanes.
    pub fn len(&self) -> usize {
        self.vals.len()
    }

    /// True if there are no lanes.
    pub fn is_empty(&self) -> bool {
        self.vals.is_empty()
    }

    /// Lane `i` in storage order.
    pub fn get(&self, i: usize) -> Option<&ExactNum> {
        self.vals.get(i)
    }

    /// Entry `(i, j)`, or `None` if out of range.
    pub fn get2(&self, i: usize, j: usize) -> Option<&ExactNum> {
        if i >= self.rows || j >= self.cols {
            return None;
        }
        self.get(i * self.cols + j)
    }

    /// All values in row-major order.
    pub fn as_slice(&self) -> &[ExactNum] {
        &self.vals
    }

    /// Elementwise add at `p`.
    pub fn add(&self, rhs: &Self) -> Option<Self> {
        self.zip(rhs, |a, b| a.add(b, self.p, RoundingMode::ToEven))
    }

    /// Add a scalar to every lane.
    pub fn add_scalar(&self, s: &ExactNum) -> Self {
        Self {
            p: self.p,
            vals: self
                .vals
                .iter()
                .map(|x| x.add(s, self.p, RoundingMode::ToEven))
                .collect(),
            rows: self.rows,
            cols: self.cols,
        }
    }

    /// Elementwise sub.
    pub fn sub(&self, rhs: &Self) -> Option<Self> {
        self.zip(rhs, |a, b| a.sub(b, self.p, RoundingMode::ToEven))
    }

    /// Elementwise mul.
    pub fn mul(&self, rhs: &Self) -> Option<Self> {
        self.zip(rhs, |a, b| a.mul(b, self.p, RoundingMode::ToEven))
    }

    /// Multiply every lane by a scalar.
    pub fn mul_scalar(&self, s: &ExactNum) -> Self {
        Self {
            p: self.p,
            vals: self
                .vals
                .iter()
                .map(|x| x.mul(s, self.p, RoundingMode::ToEven))
                .collect(),
            rows: self.rows,
            cols: self.cols,
        }
    }

    /// Elementwise div.
    pub fn div(&self, rhs: &Self) -> Option<Self> {
        self.zip(rhs, |a, b| a.div(b, self.p, RoundingMode::ToEven))
    }

    /// Sequential sum at `p`.
    pub fn sum(&self) -> ExactNum {
        let mut acc = ExactNum::from_u8(0, self.p);
        for v in &self.vals {
            acc = acc.add(v, self.p, RoundingMode::ToEven);
        }
        acc
    }

    /// Sequential dot product at `p`.
    pub fn dot(&self, rhs: &Self) -> Option<ExactNum> {
        if self.len() != rhs.len() {
            return None;
        }
        let mut acc = ExactNum::from_u8(0, self.p);
        for (a, b) in self.vals.iter().zip(rhs.vals.iter()) {
            let t = a.mul(b, self.p, RoundingMode::ToEven);
            acc = acc.add(&t, self.p, RoundingMode::ToEven);
        }
        Some(acc)
    }

    /// Software matmul: `(m×k)(k×n) → (m×n)`. Sequential mul-then-add at `p`.
    pub fn matmul(&self, rhs: &Self) -> Option<Self> {
        if self.cols != rhs.rows {
            return None;
        }
        let m = self.rows;
        let k = self.cols;
        let n = rhs.cols;
        let mut vals = Vec::with_capacity(m.checked_mul(n)?);
        for i in 0..m {
            for j in 0..n {
                let mut acc = ExactNum::from_u8(0, self.p);
                for t in 0..k {
                    let prod = self.vals[i * k + t].mul(
                        &rhs.vals[t * n + j],
                        self.p,
                        RoundingMode::ToEven,
                    );
                    acc = acc.add(&prod, self.p, RoundingMode::ToEven);
                }
                vals.push(acc);
            }
        }
        Some(Self {
            p: self.p,
            vals,
            rows: m,
            cols: n,
        })
    }

    /// Elementwise integer part.
    pub fn int(&self) -> Self {
        self.map_at(self.p, |x| x.int())
    }
    /// Elementwise fractional part.
    pub fn fract(&self) -> Self {
        self.map_at(self.p, |x| x.fract())
    }
    /// Elementwise `ceil`.
    pub fn ceil(&self) -> Self {
        self.map_at(self.p, |x| x.ceil())
    }
    /// Elementwise `floor`.
    pub fn floor(&self) -> Self {
        self.map_at(self.p, |x| x.floor())
    }
    /// Elementwise `round` with `n` binary fractional bits.
    pub fn round(&self, n: usize, rm: RoundingMode) -> Self {
        self.map_at(self.p, |x| x.round(n, rm))
    }
    /// Elementwise absolute value.
    pub fn abs(&self) -> Self {
        self.map_at(self.p, |x| x.abs())
    }
    /// Elementwise signum.
    pub fn signum(&self) -> Self {
        self.map_at(self.p, |x| x.signum())
    }
    /// Elementwise negation.
    pub fn neg(&self) -> Self {
        self.map_at(self.p, |x| x.neg())
    }
    /// Elementwise reciprocal.
    pub fn reciprocal(&self, p: usize, rm: RoundingMode) -> Self {
        self.map_at(p, |x| x.reciprocal(p, rm))
    }
    /// Elementwise `nth_root`.
    pub fn nth_root(&self, n: usize, p: usize, rm: RoundingMode) -> Self {
        self.map_at(p, |x| x.nth_root(n, p, rm))
    }
    /// Elementwise `powi`.
    pub fn powi(&self, n: usize, p: usize, rm: RoundingMode) -> Self {
        self.map_at(p, |x| x.powi(n, p, rm))
    }
    /// Elementwise `powsi`.
    pub fn powsi(&self, n: isize, p: usize, rm: RoundingMode) -> Self {
        self.map_at(p, |x| x.powsi(n, p, rm))
    }

    exact_arr_p_rm_cc!(
        sin,
        cos,
        tan,
        asin,
        acos,
        atan,
        sinh,
        cosh,
        tanh,
        asinh,
        acosh,
        atanh,
        exp,
        exp2,
        exp10,
        expm1,
        ln,
        log2,
        log10,
        log1p,
        erf,
        erfc,
        gamma,
        ln_gamma,
        digamma,
        ei,
        si,
        ci,
        li,
        fresnel_s,
        fresnel_c,
        ai,
        bi,
        elliptic_k,
        elliptic_e_complete,
        rem_pi,
    );

    /// Elementwise `(sin, cos)` with a shared argument reduction.
    pub fn sin_cos(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> (Self, Self) {
        let mut s = Vec::with_capacity(self.vals.len());
        let mut c = Vec::with_capacity(self.vals.len());
        for x in &self.vals {
            let (sv, cv) = x.sin_cos(p, rm, cc);
            s.push(sv);
            c.push(cv);
        }
        (
            Self {
                p,
                vals: s,
                rows: self.rows,
                cols: self.cols,
            },
            Self {
                p,
                vals: c,
                rows: self.rows,
                cols: self.cols,
            },
        )
    }

    /// Elementwise `(sinh, cosh)` with a shared evaluation.
    pub fn sinh_cosh(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> (Self, Self) {
        let mut s = Vec::with_capacity(self.vals.len());
        let mut c = Vec::with_capacity(self.vals.len());
        for x in &self.vals {
            let (sv, cv) = x.sinh_cosh(p, rm, cc);
            s.push(sv);
            c.push(cv);
        }
        (
            Self {
                p,
                vals: s,
                rows: self.rows,
                cols: self.cols,
            },
            Self {
                p,
                vals: c,
                rows: self.rows,
                cols: self.cols,
            },
        )
    }

    /// Elementwise `sqrt`.
    pub fn sqrt(&self, p: usize, rm: RoundingMode) -> Self {
        self.map_at(p, |x| x.sqrt(p, rm))
    }

    /// Elementwise `cbrt`.
    pub fn cbrt(&self, p: usize, rm: RoundingMode) -> Self {
        self.map_at(p, |x| x.cbrt(p, rm))
    }

    /// Elementwise `P_n(self)`.
    pub fn legendre_p(&self, n: u32, p: usize, rm: RoundingMode) -> Self {
        self.map_at(p, |x| x.legendre_p(n, p, rm))
    }

    /// Elementwise `P_n^m(self)`.
    pub fn assoc_legendre_p(&self, n: u32, m: i32, p: usize, rm: RoundingMode) -> Self {
        self.map_at(p, |x| x.assoc_legendre_p(n, m, p, rm))
    }

    /// Elementwise `hypot(self, other)`.
    pub fn hypot(&self, other: &ExactNum, p: usize, rm: RoundingMode) -> Self {
        self.map_at(p, |x| x.hypot(other, p, rm))
    }

    /// Elementwise `atan2(self, x)`. `cc` is the constants cache, not a global.
    pub fn atan2(&self, x: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |y| y.atan2(x, p, rm, cc))
    }

    /// Elementwise `pow(self, n)`.
    pub fn pow(&self, n: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |x| x.pow(n, p, rm, cc))
    }

    /// Elementwise `log(self, base)`.
    pub fn log(&self, base: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |x| x.log(base, p, rm, cc))
    }

    /// Elementwise `γ(self, x)`.
    pub fn gammainc(&self, x: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |s| s.gammainc(x, p, rm, cc))
    }

    /// Elementwise `Γ(self, x)`.
    pub fn gammainc_upper(
        &self,
        x: &ExactNum,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        self.map_at(p, |s| s.gammainc_upper(x, p, rm, cc))
    }

    /// Elementwise integer-order `J_n(self)`.
    pub fn bessel_j(&self, n: usize, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |x| x.bessel_j(n, p, rm, cc))
    }

    /// Elementwise `J_ν(self)`.
    pub fn bessel_j_nu(&self, nu: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |x| x.bessel_j_nu(nu, p, rm, cc))
    }

    /// Elementwise `Y_ν(self)`.
    pub fn bessel_y(&self, nu: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |x| x.bessel_y(nu, p, rm, cc))
    }

    /// Elementwise `I_ν(self)`.
    pub fn bessel_i(&self, nu: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |x| x.bessel_i(nu, p, rm, cc))
    }

    /// Elementwise `K_ν(self)`.
    pub fn bessel_k(&self, nu: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |x| x.bessel_k(nu, p, rm, cc))
    }

    /// Elementwise `F(self | m)`.
    pub fn elliptic_f(&self, m: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |x| x.elliptic_f(m, p, rm, cc))
    }

    /// Elementwise incomplete `E(self | m)`.
    pub fn elliptic_e(&self, m: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |x| x.elliptic_e(m, p, rm, cc))
    }

    /// Elementwise complete `Π(n, m)` with `self = n`.
    pub fn elliptic_pi_complete(
        &self,
        m: &ExactNum,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        self.map_at(p, |n| n.elliptic_pi_complete(m, p, rm, cc))
    }

    /// Elementwise `Π(self; x | m)`.
    pub fn elliptic_pi(
        &self,
        x: &ExactNum,
        m: &ExactNum,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        self.map_at(p, |n| n.elliptic_pi(x, m, p, rm, cc))
    }

    /// Elementwise `{}_2F_1(self, b; c; z)`.
    pub fn hypergeom_2f1(
        &self,
        b: &ExactNum,
        c: &ExactNum,
        z: &ExactNum,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        self.map_at(p, |a| a.hypergeom_2f1(b, c, z, p, rm, cc))
    }

    /// Elementwise `I_x(self, b)`.
    pub fn betainc(
        &self,
        b: &ExactNum,
        x: &ExactNum,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        self.map_at(p, |a| a.betainc(b, x, p, rm, cc))
    }

    /// LU with partial pivoting: `(L, U, P)` such that row `i` of `P·A` is
    /// original row `P[i]`, and `P·A = L·U` at `(p, rm)`.
    ///
    /// `L` is unit lower (`n×n`). `U` is upper (`n×m`). A zero pivot
    /// (singular) returns `None`.
    pub fn lu_decomp(&self, p: usize, rm: RoundingMode) -> Option<(Self, Self, Vec<usize>)> {
        let n = self.rows;
        let m = self.cols;
        if n == 0 || m == 0 {
            return None;
        }
        let kmax = n.min(m);
        let mut a = self.vals.clone();
        for v in &mut a {
            let _ = v.set_precision(p, rm);
        }
        let mut perm: Vec<usize> = (0..n).collect();
        let mut lvals = alloc::vec![ExactNum::from_u8(0, p); n.checked_mul(n)?];
        for i in 0..n {
            lvals[i * n + i] = ExactNum::from_u8(1, p);
        }
        let ix = |r: usize, c: usize| r * m + c;
        for k in 0..kmax {
            let mut piv = k;
            let mut best = a[ix(k, k)].abs();
            for r in (k + 1)..n {
                let t = a[ix(r, k)].abs();
                if matches!(t.cmp(&best), Some(c) if c > 0) {
                    best = t;
                    piv = r;
                }
            }
            if a[ix(piv, k)].is_zero() {
                return None;
            }
            if piv != k {
                for c in 0..m {
                    a.swap(ix(k, c), ix(piv, c));
                }
                for c in 0..k {
                    lvals.swap(k * n + c, piv * n + c);
                }
                perm.swap(k, piv);
            }
            let akk = a[ix(k, k)].clone();
            for i in (k + 1)..n {
                let lik = a[ix(i, k)].div(&akk, p, rm);
                lvals[i * n + k] = lik.clone();
                a[ix(i, k)] = ExactNum::from_u8(0, p);
                for j in (k + 1)..m {
                    let t = lik.mul(&a[ix(k, j)], p, rm);
                    a[ix(i, j)] = a[ix(i, j)].sub(&t, p, rm);
                }
            }
        }
        Some((
            Self {
                p,
                vals: lvals,
                rows: n,
                cols: n,
            },
            Self {
                p,
                vals: a,
                rows: n,
                cols: m,
            },
            perm,
        ))
    }

    fn zip(&self, rhs: &Self, op: impl Fn(&ExactNum, &ExactNum) -> ExactNum) -> Option<Self> {
        if self.rows != rhs.rows || self.cols != rhs.cols {
            return None;
        }
        Some(Self {
            p: self.p,
            vals: self
                .vals
                .iter()
                .zip(rhs.vals.iter())
                .map(|(a, b)| op(a, b))
                .collect(),
            rows: self.rows,
            cols: self.cols,
        })
    }

    fn map_at(&self, p: usize, mut op: impl FnMut(&ExactNum) -> ExactNum) -> Self {
        Self {
            p,
            vals: self.vals.iter().map(|x| op(x)).collect(),
            rows: self.rows,
            cols: self.cols,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bin64_array_add_dot() {
        let one = Ieee64::from_i32(1);
        let two = Ieee64::from_i32(2);
        let a = Ieee64Array::from_values(&[one, two]);
        let b = Ieee64Array::from_values(&[two, one]);
        let s = a.add(&b).unwrap();
        assert_eq!(s.get(0).unwrap().to_bits(), Ieee64::from_i32(3).to_bits());
        assert_eq!(s.get(1).unwrap().to_bits(), Ieee64::from_i32(3).to_bits());
        let d = a.dot(&b).unwrap();
        assert_eq!(d.to_bits(), Ieee64::from_i32(4).to_bits());
        assert_eq!(a.sum().to_bits(), Ieee64::from_i32(3).to_bits());
        let scaled = a.mul_scalar(two);
        assert_eq!(scaled.get(0).unwrap().to_bits(), two.to_bits());
    }

    #[test]
    fn bin32_array_sqrt() {
        let four = Ieee32::from_i32(4);
        let a = Ieee32Array::filled(3, four);
        let r = a.sqrt();
        assert_eq!(r.get(0).unwrap().to_bits(), Ieee32::from_i32(2).to_bits());
        assert_eq!(r.len(), 3);
    }

    #[test]
    fn bin64_array_exp_sin() {
        let mut cc = Consts::new().unwrap();
        let z = Ieee64Array::from_values(&[Ieee64::ZERO]);
        let e = z.exp(&mut cc);
        assert_eq!(e.get(0).unwrap().to_bits(), Ieee64::from_i32(1).to_bits());
        let s = z.sin(&mut cc);
        assert!(s.get(0).unwrap().is_zero());
    }

    #[test]
    fn exact_array_add() {
        let p = 64;
        let one = ExactNum::from_u8(1, p);
        let two = ExactNum::from_u8(2, p);
        let a = ExactNumArray::from_values(p, &[one.clone(), two.clone()]);
        let b = ExactNumArray::from_values(p, &[two.clone(), one.clone()]);
        let s = a.add(&b).unwrap();
        assert_eq!(s.get(0).unwrap().cmp(&ExactNum::from_u8(3, p)), Some(0));
        let d = a.dot(&b).unwrap();
        assert_eq!(d.cmp(&ExactNum::from_u8(4, p)), Some(0));
    }

    #[test]
    fn bin64_matmul_2x2() {
        let a = Ieee64Array::from_shape(
            2,
            2,
            &[
                Ieee64::from_i32(1),
                Ieee64::from_i32(2),
                Ieee64::from_i32(3),
                Ieee64::from_i32(4),
            ],
        )
        .unwrap();
        let b = Ieee64Array::from_shape(
            2,
            2,
            &[
                Ieee64::from_i32(5),
                Ieee64::from_i32(6),
                Ieee64::from_i32(7),
                Ieee64::from_i32(8),
            ],
        )
        .unwrap();
        let c = a.matmul(&b).unwrap();
        assert_eq!(c.shape(), (2, 2));
        assert_eq!(
            c.get2(0, 0).unwrap().to_bits(),
            Ieee64::from_i32(19).to_bits()
        );
        assert_eq!(
            c.get2(0, 1).unwrap().to_bits(),
            Ieee64::from_i32(22).to_bits()
        );
        assert_eq!(
            c.get2(1, 0).unwrap().to_bits(),
            Ieee64::from_i32(43).to_bits()
        );
        assert_eq!(
            c.get2(1, 1).unwrap().to_bits(),
            Ieee64::from_i32(50).to_bits()
        );
    }

    #[test]
    fn bin64_matmul_identity() {
        let i2 = Ieee64Array::from_shape(
            2,
            2,
            &[Ieee64::from_i32(1), Ieee64::ZERO, Ieee64::ZERO, Ieee64::from_i32(1)],
        )
        .unwrap();
        let a = Ieee64Array::from_shape(
            2,
            2,
            &[
                Ieee64::from_i32(1),
                Ieee64::from_i32(2),
                Ieee64::from_i32(3),
                Ieee64::from_i32(4),
            ],
        )
        .unwrap();
        let c = i2.matmul(&a).unwrap();
        assert_eq!(
            c.get2(0, 0).unwrap().to_bits(),
            Ieee64::from_i32(1).to_bits()
        );
        assert_eq!(
            c.get2(0, 1).unwrap().to_bits(),
            Ieee64::from_i32(2).to_bits()
        );
        assert_eq!(
            c.get2(1, 0).unwrap().to_bits(),
            Ieee64::from_i32(3).to_bits()
        );
        assert_eq!(
            c.get2(1, 1).unwrap().to_bits(),
            Ieee64::from_i32(4).to_bits()
        );
    }

    #[test]
    fn bin64_matmul_shape_mismatch() {
        let a = Ieee64Array::from_shape(2, 2, &[Ieee64::from_i32(1); 4]).unwrap();
        let b = Ieee64Array::from_shape(3, 1, &[Ieee64::from_i32(1); 3]).unwrap();
        assert!(a.matmul(&b).is_none());
        assert!(Ieee64Array::from_shape(2, 2, &[Ieee64::from_i32(1)]).is_none());
        let row = Ieee64Array::from_values(&[Ieee64::from_i32(1); 4]);
        assert_eq!(row.shape(), (1, 4));
        let sq = row.reshape(2, 2).unwrap();
        assert_eq!(sq.shape(), (2, 2));
        assert!(row.add(&sq).is_none());
    }

    #[test]
    fn exact_matmul_2x2() {
        let p = 64;
        let n = |k: u8| ExactNum::from_u8(k, p);
        let a = ExactNumArray::from_shape(p, 2, 2, &[n(1), n(2), n(3), n(4)]).unwrap();
        let b = ExactNumArray::from_shape(p, 2, 2, &[n(5), n(6), n(7), n(8)]).unwrap();
        let c = a.matmul(&b).unwrap();
        assert_eq!(c.shape(), (2, 2));
        assert_eq!(c.get2(0, 0).unwrap().cmp(&n(19)), Some(0));
        assert_eq!(c.get2(0, 1).unwrap().cmp(&n(22)), Some(0));
        assert_eq!(c.get2(1, 0).unwrap().cmp(&n(43)), Some(0));
        assert_eq!(c.get2(1, 1).unwrap().cmp(&n(50)), Some(0));
    }

    #[test]
    fn bin32_simd_add_mul_golds() {
        let one = Ieee32::from_i32(1);
        let two = Ieee32::from_i32(2);
        let four = Ieee32::from_i32(4);
        let a = Ieee32Array::from_values(&[one, two, one, two, one]);
        let b = Ieee32Array::from_values(&[one, two, two, one, one]);
        let s = a.add(&b).unwrap();
        assert_eq!(s.get(0).unwrap().to_bits(), two.to_bits());
        assert_eq!(s.get(1).unwrap().to_bits(), four.to_bits());
        assert_eq!(s.get(4).unwrap().to_bits(), two.to_bits());
        let half = Ieee32::from_bits(0x3F00_0000);
        let t = Ieee32Array::filled(4, two);
        let h = Ieee32Array::filled(4, half);
        let p = t.mul(&h).unwrap();
        assert_eq!(p.get(0).unwrap().to_bits(), one.to_bits());
        assert_eq!(p.get(3).unwrap().to_bits(), one.to_bits());
    }

    #[test]
    fn array_ufunc_identities() {
        let mut cc = Consts::new().unwrap();
        let z32 = Ieee32Array::from_values(&[Ieee32::ZERO]);
        assert!(z32.asin(&mut cc).get(0).unwrap().is_zero());
        assert!(z32.expm1(&mut cc).get(0).unwrap().is_zero());
        assert!(z32.log1p(&mut cc).get(0).unwrap().is_zero());
        assert_eq!(
            z32.bessel_j(0, &mut cc).get(0).unwrap().to_bits(),
            Ieee32::from_i32(1).to_bits()
        );
        let one64 = Ieee64Array::from_values(&[Ieee64::from_i32(1)]);
        assert!(one64.ln_gamma(&mut cc).get(0).unwrap().is_zero());
        let p = 64;
        let rm = RoundingMode::ToEven;
        let z = ExactNumArray::from_values(p, &[ExactNum::from_u8(0, p)]);
        assert!(z.sinh(p, rm, &mut cc).get(0).unwrap().is_zero());
        let one = ExactNumArray::from_values(p, &[ExactNum::from_u8(1, p)]);
        assert_eq!(
            one.ln_gamma(p, rm, &mut cc)
                .get(0)
                .unwrap()
                .cmp(&ExactNum::from_u8(0, p)),
            Some(0)
        );
        let x = ExactNumArray::from_values(p, &[ExactNum::from_u8(3, p)]);
        assert_eq!(
            x.legendre_p(0, p, rm)
                .get(0)
                .unwrap()
                .cmp(&ExactNum::from_u8(1, p)),
            Some(0)
        );
        assert_eq!(
            x.floor().get(0).unwrap().cmp(&ExactNum::from_u8(3, p)),
            Some(0)
        );
    }

    #[test]
    fn exact_array_sin_2x3_matches_scalar() {
        let p = 128;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new().unwrap();
        let n = |k: u8| ExactNum::from_u8(k, p);
        let vals = [n(1), n(2), n(3), n(4), n(5), n(6)];
        let a = ExactNumArray::from_shape(p, 2, 3, &vals).unwrap();
        let s = a.sin(p, rm, &mut cc);
        assert_eq!(s.shape(), (2, 3));
        for i in 0..2 {
            for j in 0..3 {
                let want = vals[i * 3 + j].sin(p, rm, &mut cc);
                assert_eq!(s.get2(i, j).unwrap().cmp(&want), Some(0));
            }
        }
        let row = ExactNumArray::from_values(p, &[n(1), n(2)]);
        assert!(a.add(&row).is_none());
    }

    #[test]
    fn exact_array_bessel_j_nu_matches_scalar() {
        let p = 128;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new().unwrap();
        let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, rm);
        let xs = [ExactNum::from_u8(1, p), ExactNum::from_u8(2, p), ExactNum::from_u8(3, p)];
        let a = ExactNumArray::from_values(p, &xs);
        let out = a.bessel_j_nu(&half, p, rm, &mut cc);
        for (i, x) in xs.iter().enumerate() {
            let want = x.bessel_j_nu(&half, p, rm, &mut cc);
            assert_eq!(out.get(i).unwrap().cmp(&want), Some(0));
        }
        let nan_in = ExactNumArray::from_values(p, &[ExactNum::nan(None)]);
        assert!(nan_in.sin(p, rm, &mut cc).get(0).unwrap().is_nan());
    }

    fn perm_rows(a: &ExactNumArray, perm: &[usize]) -> ExactNumArray {
        let (n, m) = a.shape();
        let mut vals = Vec::with_capacity(n * m);
        for &r in perm {
            for c in 0..m {
                vals.push(a.get2(r, c).unwrap().clone());
            }
        }
        ExactNumArray::from_shape(a.precision(), n, m, &vals).unwrap()
    }

    #[test]
    fn exact_lu_2x2_and_singular() {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let n = |k: u8| ExactNum::from_u8(k, p);
        let a = ExactNumArray::from_shape(p, 2, 2, &[n(2), n(1), n(4), n(3)]).unwrap();
        let (l, u, perm) = a.lu_decomp(p, rm).expect("LU");
        let pa = perm_rows(&a, &perm);
        let lu = l.matmul(&u).expect("L*U");
        assert_eq!(lu.shape(), (2, 2));
        for i in 0..2 {
            for j in 0..2 {
                assert_eq!(
                    lu.get2(i, j).unwrap().cmp(pa.get2(i, j).unwrap()),
                    Some(0),
                    "PA=LU at {i},{j}"
                );
            }
        }
        let sing = ExactNumArray::from_shape(p, 2, 2, &[n(1), n(2), n(2), n(4)]).unwrap();
        assert!(sing.lu_decomp(p, rm).is_none());
    }
}
