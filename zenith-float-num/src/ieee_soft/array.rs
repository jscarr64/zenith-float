//! Dense row-major arrays of software IEEE values and `ExactNum`.
//! A 1-D vector is stored as shape `(1, n)`.

use super::simd::{
    add_u32_lanes, add_u64_lanes, div_u32_lanes, div_u64_lanes, fma_u32_lanes, fma_u64_lanes,
    mul_u32_lanes, mul_u64_lanes, sqrt_u32_lanes, sqrt_u64_lanes, sub_u32_lanes, sub_u64_lanes,
};
use super::{Ieee32, Ieee64};
use crate::defs::RoundingMode;
use crate::Consts;
use crate::Error;
use crate::ExactNum;
use alloc::vec::Vec;

/// QR sweeps allowed per singular value before [`ExactNumArray::svd_decomp`]
/// returns `None`.
const SVD_ITER_MAX: u32 = 64;
/// QR sweeps allowed per eigenvalue before [`ExactNumArray::eigen_decomp`]
/// returns `None`.
const EIGEN_ITER_MAX: u32 = 64;
/// Maximum length of [`ExactNumArray::fft`] / [`ExactNumArray::ifft`].
const FFT_MAX_POINTS: usize = 4096;
/// Bits of working precision reserved when a superdiagonal is compared to
/// its neighboring diagonals.
const SVD_CONV_GUARD_BITS: i32 = 4;

trait LaneBits: Copy {
    fn add_lanes(a: &[Self], b: &[Self]) -> Vec<Self>;
    fn mul_lanes(a: &[Self], b: &[Self]) -> Vec<Self>;
    fn sub_lanes(a: &[Self], b: &[Self]) -> Vec<Self>;
    fn div_lanes(a: &[Self], b: &[Self]) -> Vec<Self>;
    fn sqrt_lanes(a: &[Self]) -> Vec<Self>;
    fn fma_lanes(a: &[Self], b: &[Self], c: &[Self]) -> Vec<Self>;
}

impl LaneBits for u32 {
    fn add_lanes(a: &[Self], b: &[Self]) -> Vec<Self> {
        add_u32_lanes(a, b)
    }
    fn mul_lanes(a: &[Self], b: &[Self]) -> Vec<Self> {
        mul_u32_lanes(a, b)
    }
    fn sub_lanes(a: &[Self], b: &[Self]) -> Vec<Self> {
        sub_u32_lanes(a, b)
    }
    fn div_lanes(a: &[Self], b: &[Self]) -> Vec<Self> {
        div_u32_lanes(a, b)
    }
    fn sqrt_lanes(a: &[Self]) -> Vec<Self> {
        sqrt_u32_lanes(a)
    }
    fn fma_lanes(a: &[Self], b: &[Self], c: &[Self]) -> Vec<Self> {
        fma_u32_lanes(a, b, c)
    }
}

impl LaneBits for u64 {
    fn add_lanes(a: &[Self], b: &[Self]) -> Vec<Self> {
        add_u64_lanes(a, b)
    }
    fn mul_lanes(a: &[Self], b: &[Self]) -> Vec<Self> {
        mul_u64_lanes(a, b)
    }
    fn sub_lanes(a: &[Self], b: &[Self]) -> Vec<Self> {
        sub_u64_lanes(a, b)
    }
    fn div_lanes(a: &[Self], b: &[Self]) -> Vec<Self> {
        div_u64_lanes(a, b)
    }
    fn sqrt_lanes(a: &[Self]) -> Vec<Self> {
        sqrt_u64_lanes(a)
    }
    fn fma_lanes(a: &[Self], b: &[Self], c: &[Self]) -> Vec<Self> {
        fma_u64_lanes(a, b, c)
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

            /// Elementwise sub. Shapes must match. Integer SIMD via sign-bit
            /// flip then add.
            pub fn sub(&self, rhs: &Self) -> Option<Self> {
                if self.rows != rhs.rows || self.cols != rhs.cols {
                    return None;
                }
                Some(Self {
                    bits: <$bits>::sub_lanes(&self.bits, &rhs.bits),
                    rows: self.rows,
                    cols: self.cols,
                })
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

            /// Elementwise div. Shapes must match. Integer SIMD unpack on the
            /// all-normal path; significand quotient is integer `/` per lane.
            pub fn div(&self, rhs: &Self) -> Option<Self> {
                if self.rows != rhs.rows || self.cols != rhs.cols {
                    return None;
                }
                Some(Self {
                    bits: <$bits>::div_lanes(&self.bits, &rhs.bits),
                    rows: self.rows,
                    cols: self.cols,
                })
            }

            /// Elementwise sqrt. Integer SIMD unpack on non-negative normals.
            pub fn sqrt(&self) -> Self {
                Self {
                    bits: <$bits>::sqrt_lanes(&self.bits),
                    rows: self.rows,
                    cols: self.cols,
                }
            }

            /// Elementwise fused multiply-add \(a\cdot b + c\). Shapes must match.
            /// Integer SIMD significand products on the all-normal path.
            pub fn fma(&self, b: &Self, c: &Self) -> Option<Self> {
                if self.rows != b.rows
                    || self.cols != b.cols
                    || self.rows != c.rows
                    || self.cols != c.cols
                {
                    return None;
                }
                Some(Self {
                    bits: <$bits>::fma_lanes(&self.bits, &b.bits, &c.bits),
                    rows: self.rows,
                    cols: self.cols,
                })
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

            /// Elementwise `sn(self | m)`.
            pub fn jacobi_sn(&self, m: $scalar, cc: &mut Consts) -> Self {
                let me = m.to_exact($p);
                self.map_exact($p, |x| x.jacobi_sn(&me, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise `cn(self | m)`.
            pub fn jacobi_cn(&self, m: $scalar, cc: &mut Consts) -> Self {
                let me = m.to_exact($p);
                self.map_exact($p, |x| x.jacobi_cn(&me, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise `dn(self | m)`.
            pub fn jacobi_dn(&self, m: $scalar, cc: &mut Consts) -> Self {
                let me = m.to_exact($p);
                self.map_exact($p, |x| x.jacobi_dn(&me, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise `am(self | m)`.
            pub fn jacobi_am(&self, m: $scalar, cc: &mut Consts) -> Self {
                let me = m.to_exact($p);
                self.map_exact($p, |x| x.jacobi_am(&me, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise `cd(self | m)`.
            pub fn jacobi_cd(&self, m: $scalar, cc: &mut Consts) -> Self {
                let me = m.to_exact($p);
                self.map_exact($p, |x| x.jacobi_cd(&me, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise `ns(self | m)`.
            pub fn jacobi_ns(&self, m: $scalar, cc: &mut Consts) -> Self {
                let me = m.to_exact($p);
                self.map_exact($p, |x| x.jacobi_ns(&me, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise `nc(self | m)`.
            pub fn jacobi_nc(&self, m: $scalar, cc: &mut Consts) -> Self {
                let me = m.to_exact($p);
                self.map_exact($p, |x| x.jacobi_nc(&me, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise `nd(self | m)`.
            pub fn jacobi_nd(&self, m: $scalar, cc: &mut Consts) -> Self {
                let me = m.to_exact($p);
                self.map_exact($p, |x| x.jacobi_nd(&me, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise `sc(self | m)`.
            pub fn jacobi_sc(&self, m: $scalar, cc: &mut Consts) -> Self {
                let me = m.to_exact($p);
                self.map_exact($p, |x| x.jacobi_sc(&me, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise `sd(self | m)`.
            pub fn jacobi_sd(&self, m: $scalar, cc: &mut Consts) -> Self {
                let me = m.to_exact($p);
                self.map_exact($p, |x| x.jacobi_sd(&me, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise `cs(self | m)`.
            pub fn jacobi_cs(&self, m: $scalar, cc: &mut Consts) -> Self {
                let me = m.to_exact($p);
                self.map_exact($p, |x| x.jacobi_cs(&me, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise `ds(self | m)`.
            pub fn jacobi_ds(&self, m: $scalar, cc: &mut Consts) -> Self {
                let me = m.to_exact($p);
                self.map_exact($p, |x| x.jacobi_ds(&me, $p, RoundingMode::ToEven, cc))
            }

            /// Elementwise `dc(self | m)`.
            pub fn jacobi_dc(&self, m: $scalar, cc: &mut Consts) -> Self {
                let me = m.to_exact($p);
                self.map_exact($p, |x| x.jacobi_dc(&me, $p, RoundingMode::ToEven, cc))
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

impl Ieee64Array {
    pub(crate) fn from_parts(rows: usize, cols: usize, bits: Vec<u64>) -> Result<Self, Error> {
        let n = rows.checked_mul(cols).ok_or(Error::InvalidArgument)?;
        if n != bits.len() {
            return Err(Error::InvalidArgument);
        }
        Ok(Self { bits, rows, cols })
    }
}

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

    /// Fill `shape` from `dist` at `(p, rm)`. `None` if the shape product overflows.
    #[cfg(any(test, feature = "random"))]
    pub fn random_fill(
        shape: (usize, usize),
        dist: &crate::RandomDist,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Option<Self> {
        let (rows, cols) = shape;
        let n = rows.checked_mul(cols)?;
        let mut vals = Vec::with_capacity(n);
        for _ in 0..n {
            let v = match dist {
                crate::RandomDist::Uniform(a, b) => ExactNum::random_uniform(a, b, p, rm),
                crate::RandomDist::Normal(mu, sigma) => {
                    ExactNum::random_gaussian(mu, sigma, p, rm, cc)
                }
                crate::RandomDist::Exponential(lambda) => {
                    ExactNum::random_exponential(lambda, p, rm, cc)
                }
            };
            vals.push(v);
        }
        Some(Self {
            p,
            vals,
            rows,
            cols,
        })
    }

    pub(crate) fn from_parts(
        p: usize,
        rows: usize,
        cols: usize,
        vals: Vec<ExactNum>,
    ) -> Result<Self, Error> {
        let n = rows.checked_mul(cols).ok_or(Error::InvalidArgument)?;
        if n != vals.len() {
            return Err(Error::InvalidArgument);
        }
        Ok(Self {
            p,
            vals,
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

    /// Elementwise `sn(self | m)`.
    pub fn jacobi_sn(&self, m: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |x| x.jacobi_sn(m, p, rm, cc))
    }

    /// Elementwise `cn(self | m)`.
    pub fn jacobi_cn(&self, m: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |x| x.jacobi_cn(m, p, rm, cc))
    }

    /// Elementwise `dn(self | m)`.
    pub fn jacobi_dn(&self, m: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |x| x.jacobi_dn(m, p, rm, cc))
    }

    /// Elementwise `am(self | m)`.
    pub fn jacobi_am(&self, m: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |x| x.jacobi_am(m, p, rm, cc))
    }

    /// Elementwise `cd(self | m)`.
    pub fn jacobi_cd(&self, m: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |x| x.jacobi_cd(m, p, rm, cc))
    }

    /// Elementwise `ns(self | m)`.
    pub fn jacobi_ns(&self, m: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |x| x.jacobi_ns(m, p, rm, cc))
    }

    /// Elementwise `nc(self | m)`.
    pub fn jacobi_nc(&self, m: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |x| x.jacobi_nc(m, p, rm, cc))
    }

    /// Elementwise `nd(self | m)`.
    pub fn jacobi_nd(&self, m: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |x| x.jacobi_nd(m, p, rm, cc))
    }

    /// Elementwise `sc(self | m)`.
    pub fn jacobi_sc(&self, m: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |x| x.jacobi_sc(m, p, rm, cc))
    }

    /// Elementwise `sd(self | m)`.
    pub fn jacobi_sd(&self, m: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |x| x.jacobi_sd(m, p, rm, cc))
    }

    /// Elementwise `cs(self | m)`.
    pub fn jacobi_cs(&self, m: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |x| x.jacobi_cs(m, p, rm, cc))
    }

    /// Elementwise `ds(self | m)`.
    pub fn jacobi_ds(&self, m: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |x| x.jacobi_ds(m, p, rm, cc))
    }

    /// Elementwise `dc(self | m)`.
    pub fn jacobi_dc(&self, m: &ExactNum, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        self.map_at(p, |x| x.jacobi_dc(m, p, rm, cc))
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
    /// (singular) or a failed heap reserve (`MemoryAllocation`) returns `None`.
    pub fn lu_decomp(&self, p: usize, rm: RoundingMode) -> Option<(Self, Self, Vec<usize>)> {
        let n = self.rows;
        let m = self.cols;
        if n == 0 || m == 0 {
            return None;
        }
        let kmax = n.min(m);
        let mut a = try_clone_vals(&self.vals)?;
        for v in &mut a {
            let _ = v.set_precision(p, rm);
        }
        let mut perm = try_alloc_vec(n, 0usize)?;
        for (i, slot) in perm.iter_mut().enumerate() {
            *slot = i;
        }
        let mut lvals = try_alloc_vec(n.checked_mul(n)?, ExactNum::from_u8(0, p))?;
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

    /// Row–column transpose.
    pub fn transpose(&self) -> Self {
        let mut vals = Vec::with_capacity(self.vals.len());
        for j in 0..self.cols {
            for i in 0..self.rows {
                vals.push(self.vals[i * self.cols + j].clone());
            }
        }
        Self {
            p: self.p,
            vals,
            rows: self.cols,
            cols: self.rows,
        }
    }

    /// Modified Gram–Schmidt QR at `(p, rm)`.
    ///
    /// Returns `(Q, R)` with `Q` `m×k` having orthonormal columns, `R` `k×n`
    /// upper triangular, `k = min(m, n)`. A rank-deficient column is a zero
    /// column of `Q` and a zero diagonal entry of `R` — not a panic.
    pub fn qr_decomp(&self, p: usize, rm: RoundingMode) -> Option<(Self, Self)> {
        let m = self.rows;
        let n = self.cols;
        if m == 0 || n == 0 {
            return None;
        }
        let k = m.min(n);
        let mut q = try_alloc_vec(m.checked_mul(k)?, ExactNum::from_u8(0, p))?;
        let mut r = try_alloc_vec(k.checked_mul(n)?, ExactNum::from_u8(0, p))?;
        let a = |row: usize, col: usize| -> ExactNum {
            let mut v = self.vals[row * n + col].clone();
            let _ = v.set_precision(p, rm);
            v
        };
        for j in 0..n {
            let mut v: Vec<ExactNum> = (0..m).map(|i| a(i, j)).collect();
            let jlim = j.min(k);
            for i in 0..jlim {
                let mut dot = ExactNum::from_u8(0, p);
                for t in 0..m {
                    let qi = q[t * k + i].clone();
                    dot = dot.add(&qi.mul(&v[t], p, rm), p, rm);
                }
                r[i * n + j] = dot.clone();
                for t in 0..m {
                    let qi = q[t * k + i].clone();
                    v[t] = v[t].sub(&dot.mul(&qi, p, rm), p, rm);
                }
            }
            if j < k {
                let mut nrm = ExactNum::from_u8(0, p);
                for t in 0..m {
                    nrm = nrm.add(&v[t].mul(&v[t], p, rm), p, rm);
                }
                nrm = nrm.sqrt(p, rm);
                r[j * n + j] = nrm.clone();
                if !nrm.is_zero() {
                    for t in 0..m {
                        q[t * k + j] = v[t].div(&nrm, p, rm);
                    }
                }
            }
        }
        Some((
            Self {
                p,
                vals: q,
                rows: m,
                cols: k,
            },
            Self {
                p,
                vals: r,
                rows: k,
                cols: n,
            },
        ))
    }

    /// Golub–Reinsch SVD at `(p, rm)`.
    ///
    /// Returns `(U, Σ, V^T)` with `U` `m×k` (orthonormal columns), `Σ` `k×k`
    /// diagonal (non-negative, descending), `V^T` `k×n` (orthonormal rows),
    /// `k = min(m, n)`, so that `U · Σ · V^T = A` at working precision.
    /// Empty input, a NaN/Inf entry, a failed heap reserve, or failure to
    /// converge within `SVD_ITER_MAX` sweeps per singular value returns `None`.
    pub fn svd_decomp(&self, p: usize, rm: RoundingMode) -> Option<(Self, Self, Self)> {
        let m = self.rows;
        let n = self.cols;
        if m == 0 || n == 0 {
            return None;
        }
        for v in &self.vals {
            if v.is_nan() || v.is_inf() {
                return None;
            }
        }
        if m < n {
            let (ut, s, vtt) = self.transpose().svd_decomp(p, rm)?;
            return Some((vtt.transpose(), s, ut.transpose()));
        }
        svd_decomp_tall(self, p, rm)
    }

    /// Symmetric QR eigendecomposition at `(p, rm)`.
    ///
    /// Returns `(Λ, V)` where `Λ` is a `1×n` row of eigenvalues (descending)
    /// and `V` is `n×n` with orthonormal columns, so `A V = V diag(Λ)` and
    /// `V diag(Λ) V^T = A` at working precision.
    /// Non-square, non-symmetric, empty, non-finite, or failure to converge
    /// within `EIGEN_ITER_MAX` sweeps per value returns `None`.
    pub fn eigen_decomp(&self, p: usize, rm: RoundingMode) -> Option<(Self, Self)> {
        let n = self.rows;
        if n == 0 || n != self.cols {
            return None;
        }
        for v in &self.vals {
            if v.is_nan() || v.is_inf() {
                return None;
            }
        }
        for i in 0..n {
            for j in 0..i {
                let aij = svd_copy_prec(&self.vals[i * n + j], p, rm);
                let aji = svd_copy_prec(&self.vals[j * n + i], p, rm);
                if aij.cmp(&aji) != Some(0) {
                    return None;
                }
            }
        }
        eigen_decomp_sym(self, p, rm)
    }

    /// Radix-2 Cooley–Tukey DFT at `(p, rm, cc)`.
    ///
    /// A `(1, n)` or `(n, 1)` array is real. A `(2, n)` array is complex
    /// (row 0 real, row 1 imaginary). `n` must be a power of two and at most
    /// `FFT_MAX_POINTS`. Returns a `(2, n)` spectrum (unnormalized).
    /// Empty, non-finite, or a bad shape returns `None`.
    pub fn fft(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Option<Self> {
        fft_dit(self, false, p, rm, cc)
    }

    /// Inverse radix-2 DFT at `(p, rm, cc)`. Same layout as [`Self::fft`].
    /// The result is divided by `n` (unitary inverse of the unnormalized DFT).
    pub fn ifft(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Option<Self> {
        fft_dit(self, true, p, rm, cc)
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

fn try_alloc_vec<T: Clone>(n: usize, fill: T) -> Option<Vec<T>> {
    let mut v = Vec::new();
    v.try_reserve_exact(n).ok()?;
    v.resize(n, fill);
    Some(v)
}

fn try_clone_vals(src: &[ExactNum]) -> Option<Vec<ExactNum>> {
    let mut v = Vec::new();
    v.try_reserve_exact(src.len()).ok()?;
    v.extend(src.iter().cloned());
    Some(v)
}

fn svd_zero(p: usize) -> ExactNum {
    ExactNum::from_u8(0, p)
}

fn svd_one(p: usize) -> ExactNum {
    ExactNum::from_u8(1, p)
}

fn svd_copy_prec(x: &ExactNum, p: usize, rm: RoundingMode) -> ExactNum {
    let mut y = x.clone();
    let _ = y.set_precision(p, rm);
    y
}

fn svd_identity(n: usize, p: usize) -> Option<Vec<ExactNum>> {
    let mut v = try_alloc_vec(n.checked_mul(n)?, svd_zero(p))?;
    for i in 0..n {
        v[i * n + i] = svd_one(p);
    }
    Some(v)
}

fn svd_norm(xs: &[ExactNum], p: usize, rm: RoundingMode) -> ExactNum {
    let mut n = svd_zero(p);
    for x in xs {
        n = n.hypot(x, p, rm);
    }
    n
}

fn svd_householder(
    x: &[ExactNum],
    p: usize,
    rm: RoundingMode,
) -> Option<(Vec<ExactNum>, ExactNum)> {
    if x.is_empty() {
        return None;
    }
    let norm = svd_norm(x, p, rm);
    if norm.is_zero() {
        return None;
    }
    let mut v = x.to_vec();
    let signed = if v[0].is_negative() { norm.neg() } else { norm };
    v[0] = v[0].add(&signed, p, rm);
    let mut vtv = svd_zero(p);
    for vi in &v {
        vtv = vtv.add(&vi.mul(vi, p, rm), p, rm);
    }
    if vtv.is_zero() {
        return None;
    }
    let beta = ExactNum::from_u8(2, p).div(&vtv, p, rm);
    Some((v, beta))
}

fn svd_apply_house_left(
    a: &mut [ExactNum],
    cols: usize,
    row0: usize,
    col0: usize,
    v: &[ExactNum],
    beta: &ExactNum,
    p: usize,
    rm: RoundingMode,
) {
    let vlen = v.len();
    for j in col0..cols {
        let mut s = svd_zero(p);
        for i in 0..vlen {
            s = s.add(&v[i].mul(&a[(row0 + i) * cols + j], p, rm), p, rm);
        }
        s = s.mul(beta, p, rm);
        for i in 0..vlen {
            let t = s.mul(&v[i], p, rm);
            let idx = (row0 + i) * cols + j;
            a[idx] = a[idx].sub(&t, p, rm);
        }
    }
}

fn svd_apply_house_right(
    a: &mut [ExactNum],
    rows: usize,
    cols: usize,
    row0: usize,
    col0: usize,
    v: &[ExactNum],
    beta: &ExactNum,
    p: usize,
    rm: RoundingMode,
) {
    let vlen = v.len();
    for i in row0..rows {
        let mut s = svd_zero(p);
        for t in 0..vlen {
            s = s.add(&a[i * cols + col0 + t].mul(&v[t], p, rm), p, rm);
        }
        s = s.mul(beta, p, rm);
        for t in 0..vlen {
            let tt = s.mul(&v[t], p, rm);
            let idx = i * cols + col0 + t;
            a[idx] = a[idx].sub(&tt, p, rm);
        }
    }
}

fn svd_rotg(
    a: &ExactNum,
    b: &ExactNum,
    p: usize,
    rm: RoundingMode,
) -> (ExactNum, ExactNum, ExactNum) {
    let r = a.hypot(b, p, rm);
    if r.is_zero() {
        return (svd_one(p), svd_zero(p), r);
    }
    (a.div(&r, p, rm), b.div(&r, p, rm), r)
}

fn svd_apply_givens_cols(
    mat: &mut [ExactNum],
    rows: usize,
    cols: usize,
    j0: usize,
    j1: usize,
    cs: &ExactNum,
    sn: &ExactNum,
    p: usize,
    rm: RoundingMode,
) {
    for i in 0..rows {
        let a = mat[i * cols + j0].clone();
        let b = mat[i * cols + j1].clone();
        mat[i * cols + j0] = cs.mul(&a, p, rm).add(&sn.mul(&b, p, rm), p, rm);
        mat[i * cols + j1] = cs.mul(&b, p, rm).sub(&sn.mul(&a, p, rm), p, rm);
    }
}

fn svd_negligible(e: &ExactNum, d0: &ExactNum, d1: &ExactNum, p: usize, rm: RoundingMode) -> bool {
    if e.is_zero() {
        return true;
    }
    if e.is_nan() || e.is_inf() {
        return false;
    }
    let scale = d0.abs().add(&d1.abs(), p, rm);
    if scale.is_zero() {
        return e.is_zero();
    }
    match (e.abs().exponent(), scale.exponent()) {
        (Some(ee), Some(se)) => ee < se - (p as i32 - SVD_CONV_GUARD_BITS),
        _ => false,
    }
}

fn svd_wilkinson_shift(
    d_prev: &ExactNum,
    d_last: &ExactNum,
    e_prev: &ExactNum,
    e_last: &ExactNum,
    p: usize,
    rm: RoundingMode,
) -> ExactNum {
    let two = ExactNum::from_u8(2, p);
    let b = d_prev
        .add(d_last, p, rm)
        .mul(&d_prev.sub(d_last, p, rm), p, rm)
        .add(&e_prev.mul(e_prev, p, rm), p, rm)
        .div(&two, p, rm);
    let t = d_last.mul(e_last, p, rm);
    let c = t.mul(&t, p, rm);
    if b.is_zero() && c.is_zero() {
        return svd_zero(p);
    }
    let disc = b.mul(&b, p, rm).add(&c, p, rm).sqrt(p, rm);
    let signed = if b.is_negative() { disc.neg() } else { disc };
    let denom = b.add(&signed, p, rm);
    if denom.is_zero() {
        return svd_zero(p);
    }
    c.div(&denom, p, rm)
}

fn svd_qr_sweep(
    d: &mut [ExactNum],
    e: &mut [ExactNum],
    u: &mut [ExactNum],
    v: &mut [ExactNum],
    m: usize,
    n: usize,
    p_blk: usize,
    q_blk: usize,
    p: usize,
    rm: RoundingMode,
) {
    let last = q_blk - 1;
    let e_prev = if last >= p_blk + 2 { e[last - 2].clone() } else { svd_zero(p) };
    let shift = svd_wilkinson_shift(&d[last - 1], &d[last], &e_prev, &e[last - 1], p, rm);
    let mut f = d[p_blk]
        .add(&d[last], p, rm)
        .mul(&d[p_blk].sub(&d[last], p, rm), p, rm)
        .add(&shift, p, rm);
    let mut g = d[p_blk].mul(&e[p_blk], p, rm);
    for j in p_blk..last {
        let (cs, sn, r) = svd_rotg(&f, &g, p, rm);
        if j > p_blk {
            e[j - 1] = r;
        }
        let dj = d[j].clone();
        let ej = e[j].clone();
        let dj1 = d[j + 1].clone();
        f = cs.mul(&dj, p, rm).add(&sn.mul(&ej, p, rm), p, rm);
        e[j] = cs.mul(&ej, p, rm).sub(&sn.mul(&dj, p, rm), p, rm);
        g = sn.mul(&dj1, p, rm);
        d[j + 1] = cs.mul(&dj1, p, rm);
        svd_apply_givens_cols(v, n, n, j, j + 1, &cs, &sn, p, rm);

        let (cs, sn, r) = svd_rotg(&f, &g, p, rm);
        d[j] = r;
        let ej = e[j].clone();
        let dj1 = d[j + 1].clone();
        f = cs.mul(&ej, p, rm).add(&sn.mul(&dj1, p, rm), p, rm);
        d[j + 1] = cs.mul(&dj1, p, rm).sub(&sn.mul(&ej, p, rm), p, rm);
        if j + 1 < last {
            g = sn.mul(&e[j + 1], p, rm);
            e[j + 1] = cs.mul(&e[j + 1], p, rm);
        }
        svd_apply_givens_cols(u, m, m, j, j + 1, &cs, &sn, p, rm);
    }
    e[last - 1] = f;
}

fn svd_zero_last_super(
    d: &mut [ExactNum],
    e: &mut [ExactNum],
    v: &mut [ExactNum],
    n: usize,
    p_blk: usize,
    q_blk: usize,
    p: usize,
    rm: RoundingMode,
) {
    let k = q_blk - 1;
    let mut f = e[k - 1].clone();
    e[k - 1] = svd_zero(p);
    for j in (p_blk..k).rev() {
        let (cs, sn, t) = svd_rotg(&d[j], &f, p, rm);
        d[j] = t;
        if j > p_blk {
            f = sn.neg().mul(&e[j - 1], p, rm);
            e[j - 1] = cs.mul(&e[j - 1], p, rm);
        }
        svd_apply_givens_cols(v, n, n, j, k, &cs, &sn, p, rm);
    }
}

fn svd_zero_first_super(
    d: &mut [ExactNum],
    e: &mut [ExactNum],
    u: &mut [ExactNum],
    m: usize,
    p_blk: usize,
    q_blk: usize,
    p: usize,
    rm: RoundingMode,
) {
    let mut f = e[p_blk].clone();
    e[p_blk] = svd_zero(p);
    for j in (p_blk + 1)..q_blk {
        let (cs, sn, t) = svd_rotg(&d[j], &f, p, rm);
        d[j] = t;
        if j + 1 < q_blk {
            f = sn.neg().mul(&e[j], p, rm);
            e[j] = cs.mul(&e[j], p, rm);
        }
        svd_apply_givens_cols(u, m, m, p_blk, j, &cs, &sn, p, rm);
    }
}

fn svd_take_cols(vals: &[ExactNum], rows: usize, cols: usize, k: usize, p: usize) -> ExactNumArray {
    let mut out = Vec::with_capacity(rows * k);
    for i in 0..rows {
        for j in 0..k {
            out.push(vals[i * cols + j].clone());
        }
    }
    ExactNumArray {
        p,
        vals: out,
        rows,
        cols: k,
    }
}

fn svd_vt_from_v(v: &[ExactNum], n: usize, k: usize, p: usize) -> ExactNumArray {
    let mut out = Vec::with_capacity(k * n);
    for j in 0..k {
        for i in 0..n {
            out.push(v[i * n + j].clone());
        }
    }
    ExactNumArray {
        p,
        vals: out,
        rows: k,
        cols: n,
    }
}

fn svd_decomp_tall(
    a0: &ExactNumArray,
    p: usize,
    rm: RoundingMode,
) -> Option<(ExactNumArray, ExactNumArray, ExactNumArray)> {
    let m = a0.rows;
    let n = a0.cols;
    let mut a: Vec<ExactNum> = a0.vals.iter().map(|x| svd_copy_prec(x, p, rm)).collect();
    let mut u = svd_identity(m, p)?;
    let mut v = svd_identity(n, p)?;

    for k in 0..n {
        let x: Vec<ExactNum> = (k..m).map(|i| a[i * n + k].clone()).collect();
        if let Some((hv, beta)) = svd_householder(&x, p, rm) {
            svd_apply_house_left(&mut a, n, k, k, &hv, &beta, p, rm);
            svd_apply_house_right(&mut u, m, m, 0, k, &hv, &beta, p, rm);
        }
        if k + 1 < n {
            let x: Vec<ExactNum> = ((k + 1)..n).map(|j| a[k * n + j].clone()).collect();
            if let Some((hv, beta)) = svd_householder(&x, p, rm) {
                svd_apply_house_right(&mut a, m, n, k, k + 1, &hv, &beta, p, rm);
                svd_apply_house_right(&mut v, n, n, 0, k + 1, &hv, &beta, p, rm);
            }
        }
    }

    let mut d: Vec<ExactNum> = (0..n).map(|i| a[i * n + i].clone()).collect();
    let mut e: Vec<ExactNum> = if n >= 2 {
        (0..n - 1).map(|i| a[i * n + i + 1].clone()).collect()
    } else {
        Vec::new()
    };

    let max_sweeps = SVD_ITER_MAX.saturating_mul(n.max(1) as u32);
    let mut sweeps = 0u32;
    loop {
        if n == 1 {
            break;
        }
        for i in 0..n - 1 {
            if svd_negligible(&e[i], &d[i], &d[i + 1], p, rm) {
                e[i] = svd_zero(p);
            }
        }
        if e.iter().all(|x| x.is_zero()) {
            break;
        }
        if sweeps >= max_sweeps {
            return None;
        }
        let mut q = n;
        while q > 1 && e[q - 2].is_zero() {
            q -= 1;
        }
        let mut p_blk = q - 1;
        while p_blk > 0 && !e[p_blk - 1].is_zero() {
            p_blk -= 1;
        }
        if q - p_blk < 2 {
            break;
        }

        let mut did_split = false;
        for i in p_blk..q {
            let el = if i > p_blk { e[i - 1].clone() } else { svd_zero(p) };
            let er = if i + 1 < q { e[i].clone() } else { svd_zero(p) };
            if svd_negligible(&d[i], &el, &er, p, rm) {
                d[i] = svd_zero(p);
                if i == q - 1 && i > p_blk {
                    svd_zero_last_super(&mut d, &mut e, &mut v, n, p_blk, q, p, rm);
                } else if i < q - 1 {
                    svd_zero_first_super(&mut d, &mut e, &mut u, m, i, q, p, rm);
                }
                did_split = true;
                break;
            }
        }
        if did_split {
            sweeps += 1;
            continue;
        }
        for di in d.iter().chain(e.iter()) {
            if di.is_nan() || di.is_inf() {
                return None;
            }
        }
        svd_qr_sweep(&mut d, &mut e, &mut u, &mut v, m, n, p_blk, q, p, rm);
        sweeps += 1;
    }

    for i in 0..n {
        if d[i].is_negative() {
            d[i] = d[i].neg();
            for r in 0..m {
                let idx = r * m + i;
                u[idx] = u[idx].neg();
            }
        }
    }
    for i in 0..n {
        let mut best = i;
        for j in (i + 1)..n {
            if matches!(d[j].cmp(&d[best]), Some(c) if c > 0) {
                best = j;
            }
        }
        if best != i {
            d.swap(i, best);
            for r in 0..m {
                u.swap(r * m + i, r * m + best);
            }
            for r in 0..n {
                v.swap(r * n + i, r * n + best);
            }
        }
    }

    let k = n;
    let mut sigma = try_alloc_vec(k.checked_mul(k)?, svd_zero(p))?;
    for i in 0..k {
        sigma[i * k + i] = d[i].clone();
    }
    Some((
        svd_take_cols(&u, m, m, k, p),
        ExactNumArray {
            p,
            vals: sigma,
            rows: k,
            cols: k,
        },
        svd_vt_from_v(&v, n, k, p),
    ))
}

fn eigen_wilkinson(
    a: &ExactNum,
    b: &ExactNum,
    c: &ExactNum,
    p: usize,
    rm: RoundingMode,
) -> ExactNum {
    let half = svd_one(p).div(&ExactNum::from_u8(2, p), p, rm);
    let delta = a.sub(c, p, rm).mul(&half, p, rm);
    if delta.is_zero() && b.is_zero() {
        return c.clone();
    }
    let h = delta.hypot(b, p, rm);
    let signed = if delta.is_negative() { h.neg() } else { h };
    let denom = delta.add(&signed, p, rm);
    if denom.is_zero() {
        return c.sub(&b.abs(), p, rm);
    }
    c.sub(&b.mul(b, p, rm).div(&denom, p, rm), p, rm)
}

fn eigen_qr_sweep(
    d: &mut [ExactNum],
    e: &mut [ExactNum],
    q: &mut [ExactNum],
    n: usize,
    p_blk: usize,
    q_blk: usize,
    p: usize,
    rm: RoundingMode,
) {
    let last = q_blk - 1;
    let mu = eigen_wilkinson(&d[last - 1], &e[last - 1], &d[last], p, rm);
    let mut f = d[p_blk].sub(&mu, p, rm);
    let mut g = e[p_blk].clone();
    let two = ExactNum::from_u8(2, p);
    for k in p_blk..last {
        let (cs, sn, r) = svd_rotg(&f, &g, p, rm);
        if k > p_blk {
            e[k - 1] = r;
        }
        let d0 = d[k].clone();
        let ee = e[k].clone();
        let d1 = d[k + 1].clone();
        let c2 = cs.mul(&cs, p, rm);
        let s2 = sn.mul(&sn, p, rm);
        let cs2 = cs.mul(&sn, p, rm);
        let two_cse = two.mul(&cs2.mul(&ee, p, rm), p, rm);
        d[k] = c2
            .mul(&d0, p, rm)
            .add(&two_cse, p, rm)
            .add(&s2.mul(&d1, p, rm), p, rm);
        d[k + 1] = s2
            .mul(&d0, p, rm)
            .sub(&two_cse, p, rm)
            .add(&c2.mul(&d1, p, rm), p, rm);
        e[k] = cs2
            .mul(&d1.sub(&d0, p, rm), p, rm)
            .add(&c2.sub(&s2, p, rm).mul(&ee, p, rm), p, rm);
        svd_apply_givens_cols(q, n, n, k, k + 1, &cs, &sn, p, rm);
        if k + 1 < last {
            let ek1 = e[k + 1].clone();
            f = e[k].clone();
            g = sn.mul(&ek1, p, rm);
            e[k + 1] = cs.mul(&ek1, p, rm);
        }
    }
}

fn eigen_decomp_sym(
    a0: &ExactNumArray,
    p: usize,
    rm: RoundingMode,
) -> Option<(ExactNumArray, ExactNumArray)> {
    let n = a0.rows;
    let mut a: Vec<ExactNum> = a0.vals.iter().map(|x| svd_copy_prec(x, p, rm)).collect();
    let mut q = svd_identity(n, p)?;
    for k in 0..n.saturating_sub(2) {
        let x: Vec<ExactNum> = ((k + 1)..n).map(|i| a[i * n + k].clone()).collect();
        if let Some((hv, beta)) = svd_householder(&x, p, rm) {
            svd_apply_house_left(&mut a, n, k + 1, k, &hv, &beta, p, rm);
            svd_apply_house_right(&mut a, n, n, 0, k + 1, &hv, &beta, p, rm);
            svd_apply_house_right(&mut q, n, n, 0, k + 1, &hv, &beta, p, rm);
        }
    }
    let mut d: Vec<ExactNum> = (0..n).map(|i| a[i * n + i].clone()).collect();
    let mut e: Vec<ExactNum> = if n >= 2 {
        (0..n - 1).map(|i| a[i * n + i + 1].clone()).collect()
    } else {
        Vec::new()
    };

    let max_sweeps = EIGEN_ITER_MAX.saturating_mul(n.max(1) as u32);
    let mut sweeps = 0u32;
    loop {
        if n == 1 {
            break;
        }
        for i in 0..n - 1 {
            if svd_negligible(&e[i], &d[i], &d[i + 1], p, rm) {
                e[i] = svd_zero(p);
            }
        }
        if e.iter().all(|x| x.is_zero()) {
            break;
        }
        if sweeps >= max_sweeps {
            return None;
        }
        let mut q_blk = n;
        while q_blk > 1 && e[q_blk - 2].is_zero() {
            q_blk -= 1;
        }
        let mut p_blk = q_blk - 1;
        while p_blk > 0 && !e[p_blk - 1].is_zero() {
            p_blk -= 1;
        }
        if q_blk - p_blk < 2 {
            break;
        }
        for di in d.iter().chain(e.iter()) {
            if di.is_nan() || di.is_inf() {
                return None;
            }
        }
        eigen_qr_sweep(&mut d, &mut e, &mut q, n, p_blk, q_blk, p, rm);
        sweeps += 1;
    }

    for i in 0..n {
        let mut best = i;
        for j in (i + 1)..n {
            if matches!(d[j].cmp(&d[best]), Some(c) if c > 0) {
                best = j;
            }
        }
        if best != i {
            d.swap(i, best);
            for r in 0..n {
                q.swap(r * n + i, r * n + best);
            }
        }
    }
    Some((
        ExactNumArray {
            p,
            vals: d,
            rows: 1,
            cols: n,
        },
        ExactNumArray {
            p,
            vals: q,
            rows: n,
            cols: n,
        },
    ))
}

fn fft_from_len(n: usize, p: usize) -> ExactNum {
    ExactNum::from_word(n as crate::defs::Word, p)
}

fn fft_bitrev(mut i: usize, logn: u32) -> usize {
    let mut r = 0usize;
    for _ in 0..logn {
        r = (r << 1) | (i & 1);
        i >>= 1;
    }
    r
}

fn fft_split(
    a: &ExactNumArray,
    p: usize,
    rm: RoundingMode,
) -> Option<(usize, Vec<ExactNum>, Vec<ExactNum>)> {
    let (rows, cols) = a.shape();
    let pack = |n: usize,
                re: Vec<ExactNum>,
                im: Vec<ExactNum>|
     -> Option<(usize, Vec<ExactNum>, Vec<ExactNum>)> {
        if n == 0 || !n.is_power_of_two() || n > FFT_MAX_POINTS {
            return None;
        }
        Some((n, re, im))
    };
    if rows == 2 && cols > 0 {
        let mut re = Vec::with_capacity(cols);
        let mut im = Vec::with_capacity(cols);
        for j in 0..cols {
            let r = svd_copy_prec(&a.vals[j], p, rm);
            let i = svd_copy_prec(&a.vals[cols + j], p, rm);
            if r.is_nan() || r.is_inf() || i.is_nan() || i.is_inf() {
                return None;
            }
            re.push(r);
            im.push(i);
        }
        pack(cols, re, im)
    } else if (rows == 1 && cols > 0) || (cols == 1 && rows > 0) {
        let n = a.vals.len();
        let mut re = Vec::with_capacity(n);
        for v in &a.vals {
            let r = svd_copy_prec(v, p, rm);
            if r.is_nan() || r.is_inf() {
                return None;
            }
            re.push(r);
        }
        pack(n, re, alloc::vec![svd_zero(p); n])
    } else {
        None
    }
}

fn fft_dit(
    a: &ExactNumArray,
    inverse: bool,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNumArray> {
    let (n, mut re, mut im) = fft_split(a, p, rm)?;
    let logn = n.trailing_zeros();
    for i in 0..n {
        let j = fft_bitrev(i, logn);
        if j > i {
            re.swap(i, j);
            im.swap(i, j);
        }
    }
    let two_pi = ExactNum::from_u8(2, p).mul(&cc.pi(p, rm), p, rm);
    let mut m = 2usize;
    while m <= n {
        let ang = two_pi.div(&fft_from_len(m, p), p, rm);
        let (sn, cs) = ang.sin_cos(p, rm, cc);
        let wm_re = cs;
        let wm_im = if inverse { sn } else { sn.neg() };
        let half = m / 2;
        let mut k = 0usize;
        while k < n {
            let mut w_re = svd_one(p);
            let mut w_im = svd_zero(p);
            for j in 0..half {
                let t = k + j + half;
                let u = k + j;
                let tr = w_re.mul(&re[t], p, rm).sub(&w_im.mul(&im[t], p, rm), p, rm);
                let ti = w_re.mul(&im[t], p, rm).add(&w_im.mul(&re[t], p, rm), p, rm);
                let ur = re[u].clone();
                let ui = im[u].clone();
                re[u] = ur.add(&tr, p, rm);
                im[u] = ui.add(&ti, p, rm);
                re[t] = ur.sub(&tr, p, rm);
                im[t] = ui.sub(&ti, p, rm);
                let nr = w_re.mul(&wm_re, p, rm).sub(&w_im.mul(&wm_im, p, rm), p, rm);
                let ni = w_re.mul(&wm_im, p, rm).add(&w_im.mul(&wm_re, p, rm), p, rm);
                w_re = nr;
                w_im = ni;
            }
            k += m;
        }
        m *= 2;
    }
    if inverse {
        let inv_n = svd_one(p).div(&fft_from_len(n, p), p, rm);
        for i in 0..n {
            re[i] = re[i].mul(&inv_n, p, rm);
            im[i] = im[i].mul(&inv_n, p, rm);
        }
    }
    let mut vals = Vec::with_capacity(n.checked_mul(2)?);
    vals.extend(re);
    vals.extend(im);
    Some(ExactNumArray {
        p,
        vals,
        rows: 2,
        cols: n,
    })
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
    fn ieee64_simd_1000_add_mul_div_sqrt() {
        const N: usize = 1000;
        let p = 128;
        let rm = RoundingMode::ToEven;
        let vals: Vec<Ieee64> = (1..=N as i32).map(Ieee64::from_i32).collect();
        let ones: Vec<Ieee64> = (0..N).map(|_| Ieee64::from_i32(1)).collect();
        let twos: Vec<Ieee64> = (0..N).map(|_| Ieee64::from_i32(2)).collect();
        let a = Ieee64Array::from_values(&vals);
        let one = Ieee64Array::from_values(&ones);
        let two = Ieee64Array::from_values(&twos);
        let add = a.add(&one).unwrap();
        let mul = a.mul(&two).unwrap();
        let div = a.div(&a).unwrap();
        let squares = a.mul(&a).unwrap();
        let sq = squares.sqrt();
        let sub = a.sub(&one).unwrap();
        let fma = a.fma(&one, &one).unwrap();
        assert_eq!(add.len(), N);
        for i in 0..N {
            let ai = a.get(i).unwrap();
            let oi = one.get(i).unwrap();
            let ti = two.get(i).unwrap();
            assert_eq!(add.get(i).unwrap().to_bits(), ai.add(oi).to_bits());
            assert_eq!(mul.get(i).unwrap().to_bits(), ai.mul(ti).to_bits());
            assert_eq!(div.get(i).unwrap().to_bits(), ai.div(ai).to_bits());
            assert_eq!(
                sq.get(i).unwrap().to_bits(),
                ai.mul(ai).sqrt().to_bits()
            );
            assert_eq!(sub.get(i).unwrap().to_bits(), ai.sub(oi).to_bits());
            assert_eq!(
                fma.get(i).unwrap().to_bits(),
                ai.mul_add(oi, oi).to_bits()
            );
            let xa = ai.to_exact(p);
            let x1 = oi.to_exact(p);
            let x2 = ti.to_exact(p);
            assert_eq!(
                add.get(i).unwrap().to_bits(),
                Ieee64::from_exact(&xa.add(&x1, p, rm)).to_bits()
            );
            assert_eq!(
                mul.get(i).unwrap().to_bits(),
                Ieee64::from_exact(&xa.mul(&x2, p, rm)).to_bits()
            );
            assert_eq!(
                div.get(i).unwrap().to_bits(),
                Ieee64::from_exact(&xa.div(&xa, p, rm)).to_bits()
            );
            let sqe = xa.mul(&xa, p, rm).sqrt(p, rm);
            assert_eq!(sq.get(i).unwrap().to_bits(), Ieee64::from_exact(&sqe).to_bits());
        }
        assert!(a.fma(&one, &Ieee64Array::from_values(&vals[..10])).is_none());
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

    #[test]
    fn exact_array_jacobi_sn_matches_scalar() {
        let p = 128;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new().unwrap();
        let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, rm);
        let xs = [
            ExactNum::from_u8(0, p),
            ExactNum::from_u8(1, p),
            ExactNum::from_u8(2, p),
        ];
        let a = ExactNumArray::from_values(p, &xs);
        let out = a.jacobi_sn(&half, p, rm, &mut cc);
        let out_ns = a.jacobi_ns(&half, p, rm, &mut cc);
        for (i, x) in xs.iter().enumerate() {
            let want = x.jacobi_sn(&half, p, rm, &mut cc);
            assert_eq!(out.get(i).unwrap().cmp(&want), Some(0));
            if !x.is_zero() {
                let want_ns = x.jacobi_ns(&half, p, rm, &mut cc);
                assert_eq!(out_ns.get(i).unwrap().cmp(&want_ns), Some(0));
            }
        }
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

    fn near_num(a: &ExactNum, b: &ExactNum, p: usize) -> bool {
        let d = a.sub(b, p, RoundingMode::None).abs();
        d.is_zero() || d.exponent().unwrap_or(0) < -((p as i32) - 40)
    }

    #[test]
    fn exact_qr_recon_orthog_rankdef() {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let n = |k: u8| ExactNum::from_u8(k, p);
        let a = ExactNumArray::from_shape(p, 2, 2, &[n(2), n(1), n(4), n(3)]).unwrap();
        let (q, r) = a.qr_decomp(p, rm).expect("QR");
        let qr = q.matmul(&r).expect("Q*R");
        for i in 0..2 {
            for j in 0..2 {
                assert!(
                    near_num(qr.get2(i, j).unwrap(), a.get2(i, j).unwrap(), p),
                    "QR=A at {i},{j}"
                );
            }
        }
        let qtq = q.transpose().matmul(&q).expect("Q^T Q");
        let one = n(1);
        let zero = n(0);
        assert!(near_num(qtq.get2(0, 0).unwrap(), &one, p));
        assert!(near_num(qtq.get2(1, 1).unwrap(), &one, p));
        assert!(near_num(qtq.get2(0, 1).unwrap(), &zero, p));
        assert!(near_num(qtq.get2(1, 0).unwrap(), &zero, p));

        let def = ExactNumArray::from_shape(p, 2, 2, &[n(1), n(2), n(2), n(4)]).unwrap();
        let (qd, rd) = def.qr_decomp(p, rm).expect("rank-def QR");
        let _ = qd;
        assert!(rd.get2(1, 1).unwrap().is_zero() || near_num(rd.get2(1, 1).unwrap(), &zero, p));
        let recon = qd.matmul(&rd).expect("Qd Rd");
        for i in 0..2 {
            for j in 0..2 {
                assert!(near_num(
                    recon.get2(i, j).unwrap(),
                    def.get2(i, j).unwrap(),
                    p
                ));
            }
        }
    }

    #[test]
    fn exact_svd_diag_3_2() {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let n = |k: u8| ExactNum::from_u8(k, p);
        let a = ExactNumArray::from_shape(p, 2, 2, &[n(3), n(0), n(0), n(2)]).unwrap();
        let (u, s, vt) = a.svd_decomp(p, rm).expect("SVD");
        assert_eq!(s.get2(0, 0).unwrap().cmp(&n(3)), Some(0));
        assert_eq!(s.get2(1, 1).unwrap().cmp(&n(2)), Some(0));
        assert!(near_num(s.get2(0, 1).unwrap(), &n(0), p));
        assert!(near_num(s.get2(1, 0).unwrap(), &n(0), p));
        let us = u.matmul(&s).expect("U Σ");
        let recon = us.matmul(&vt).expect("U Σ V^T");
        for i in 0..2 {
            for j in 0..2 {
                assert!(
                    near_num(recon.get2(i, j).unwrap(), a.get2(i, j).unwrap(), p),
                    "UΣV^T=A at {i},{j}"
                );
            }
        }
    }

    #[test]
    fn exact_svd_recon_orthog() {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let n = |k: u8| ExactNum::from_u8(k, p);
        let a = ExactNumArray::from_shape(p, 2, 2, &[n(2), n(1), n(4), n(3)]).unwrap();
        let (u, s, vt) = a.svd_decomp(p, rm).expect("SVD");
        let us = u.matmul(&s).expect("U Σ");
        let recon = us.matmul(&vt).expect("U Σ V^T");
        for i in 0..2 {
            for j in 0..2 {
                assert!(
                    near_num(recon.get2(i, j).unwrap(), a.get2(i, j).unwrap(), p),
                    "UΣV^T=A at {i},{j}"
                );
            }
        }
        let utu = u.transpose().matmul(&u).expect("U^T U");
        let v = vt.transpose();
        let vtv = vt.matmul(&v).expect("V^T V");
        let one = n(1);
        let zero = n(0);
        for (name, g) in [("U^T U", &utu), ("V^T V", &vtv)] {
            assert!(near_num(g.get2(0, 0).unwrap(), &one, p), "{name}[0,0]");
            assert!(near_num(g.get2(1, 1).unwrap(), &one, p), "{name}[1,1]");
            assert!(near_num(g.get2(0, 1).unwrap(), &zero, p), "{name}[0,1]");
            assert!(near_num(g.get2(1, 0).unwrap(), &zero, p), "{name}[1,0]");
        }
        assert!(ExactNumArray::from_shape(p, 0, 0, &[])
            .unwrap()
            .svd_decomp(p, rm)
            .is_none());

        let wide =
            ExactNumArray::from_shape(p, 2, 3, &[n(1), n(0), n(0), n(0), n(2), n(0)]).unwrap();
        let (uw, sw, vtw) = wide.svd_decomp(p, rm).expect("wide SVD");
        assert_eq!(sw.get2(0, 0).unwrap().cmp(&n(2)), Some(0));
        assert_eq!(sw.get2(1, 1).unwrap().cmp(&n(1)), Some(0));
        let usw = uw.matmul(&sw).expect("Uw Σw");
        let recw = usw.matmul(&vtw).expect("wide recon");
        for i in 0..2 {
            for j in 0..3 {
                assert!(near_num(
                    recw.get2(i, j).unwrap(),
                    wide.get2(i, j).unwrap(),
                    p
                ));
            }
        }
    }

    fn eigen_diag(evals: &ExactNumArray, p: usize) -> ExactNumArray {
        let n = evals.cols;
        let mut vals = Vec::with_capacity(n * n);
        let z = ExactNum::from_u8(0, p);
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    vals.push(evals.get(i).unwrap().clone());
                } else {
                    vals.push(z.clone());
                }
            }
        }
        ExactNumArray::from_shape(p, n, n, &vals).unwrap()
    }

    #[test]
    fn exact_eigen_sym_2x2() {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let n = |k: u8| ExactNum::from_u8(k, p);
        let a = ExactNumArray::from_shape(p, 2, 2, &[n(2), n(1), n(1), n(2)]).unwrap();
        let (evals, v) = a.eigen_decomp(p, rm).expect("eigen");
        assert_eq!(evals.shape(), (1, 2));
        assert!(near_num(evals.get(0).unwrap(), &n(3), p));
        assert!(near_num(evals.get(1).unwrap(), &n(1), p));
        let vtv = v.transpose().matmul(&v).expect("V^T V");
        let one = n(1);
        let zero = n(0);
        assert!(near_num(vtv.get2(0, 0).unwrap(), &one, p));
        assert!(near_num(vtv.get2(1, 1).unwrap(), &one, p));
        assert!(near_num(vtv.get2(0, 1).unwrap(), &zero, p));
        assert!(near_num(vtv.get2(1, 0).unwrap(), &zero, p));
        let av = a.matmul(&v).expect("A V");
        let lam = eigen_diag(&evals, p);
        let vl = v.matmul(&lam).expect("V Λ");
        for i in 0..2 {
            for j in 0..2 {
                assert!(
                    near_num(av.get2(i, j).unwrap(), vl.get2(i, j).unwrap(), p),
                    "Av=λv at {i},{j}"
                );
            }
        }
        let vlvt = vl.matmul(&v.transpose()).expect("V Λ V^T");
        for i in 0..2 {
            for j in 0..2 {
                assert!(near_num(vlvt.get2(i, j).unwrap(), a.get2(i, j).unwrap(), p));
            }
        }
        let nosym = ExactNumArray::from_shape(p, 2, 2, &[n(1), n(2), n(0), n(1)]).unwrap();
        assert!(nosym.eigen_decomp(p, rm).is_none());

        let a3 = ExactNumArray::from_shape(
            p,
            3,
            3,
            &[n(2), n(1), n(0), n(1), n(2), n(1), n(0), n(1), n(2)],
        )
        .unwrap();
        let (w3, v3) = a3.eigen_decomp(p, rm).expect("eigen 3");
        let s2 = n(2).sqrt(p, rm);
        let want = [n(2).add(&s2, p, rm), n(2), n(2).sub(&s2, p, rm)];
        for (i, wi) in want.iter().enumerate() {
            assert!(near_num(w3.get(i).unwrap(), wi, p), "λ[{i}]");
        }
        let av3 = a3.matmul(&v3).expect("A3 V");
        let vl3 = v3.matmul(&eigen_diag(&w3, p)).expect("V3 Λ");
        for i in 0..3 {
            for j in 0..3 {
                assert!(near_num(
                    av3.get2(i, j).unwrap(),
                    vl3.get2(i, j).unwrap(),
                    p
                ));
            }
        }
    }

    #[test]
    fn exact_fft_impulse_cosine_parseval() {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new().unwrap();
        let n = |k: u8| ExactNum::from_u8(k, p);
        let impulse = ExactNumArray::from_values(p, &[n(1), n(0), n(0), n(0)]);
        let spec = impulse.fft(p, rm, &mut cc).expect("FFT impulse");
        assert_eq!(spec.shape(), (2, 4));
        for j in 0..4 {
            assert!(near_num(spec.get2(0, j).unwrap(), &n(1), p), "re[{j}]");
            assert!(near_num(spec.get2(1, j).unwrap(), &n(0), p), "im[{j}]");
        }

        let n8 = ExactNum::from_u8(8, p);
        let two_pi = n(2).mul(&cc.pi(p, rm), p, rm);
        let mut cos_vals = Vec::with_capacity(8);
        for k in 0..8u8 {
            let kn = ExactNum::from_u8(k, p);
            let ang = two_pi.mul(&kn, p, rm).div(&n8, p, rm);
            cos_vals.push(ang.cos(p, rm, &mut cc));
        }
        let cosine = ExactNumArray::from_values(p, &cos_vals);
        let cspec = cosine.fft(p, rm, &mut cc).expect("FFT cos");
        let four = n(4);
        let zero = n(0);
        for j in 0..8 {
            let re = cspec.get2(0, j).unwrap();
            let im = cspec.get2(1, j).unwrap();
            if j == 1 || j == 7 {
                assert!(near_num(re, &four, p), "cos bin {j} re");
            } else {
                assert!(near_num(re, &zero, p), "cos bin {j} re");
            }
            assert!(near_num(im, &zero, p), "cos bin {j} im");
        }

        let back = cspec.ifft(p, rm, &mut cc).expect("IFFT");
        for j in 0..8 {
            assert!(near_num(
                back.get2(0, j).unwrap(),
                cosine.get(j).unwrap(),
                p
            ));
            assert!(near_num(back.get2(1, j).unwrap(), &zero, p));
        }

        let mut e_t = ExactNum::from_u8(0, p);
        let mut e_f = ExactNum::from_u8(0, p);
        for j in 0..8 {
            let x = cosine.get(j).unwrap();
            e_t = e_t.add(&x.mul(x, p, rm), p, rm);
            let xr = cspec.get2(0, j).unwrap();
            let xi = cspec.get2(1, j).unwrap();
            e_f = e_f
                .add(&xr.mul(xr, p, rm), p, rm)
                .add(&xi.mul(xi, p, rm), p, rm);
        }
        let parseval = e_f.div(&n8, p, rm);
        assert!(near_num(&parseval, &e_t, p));
        assert!(ExactNumArray::from_values(p, &[n(1), n(2), n(3)])
            .fft(p, rm, &mut cc)
            .is_none());
    }
}
