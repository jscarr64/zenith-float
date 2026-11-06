//! Discrete cosine / sine transforms, real-input FFT wrappers, and windows.

use crate::defs::WORD_BIT_SIZE;
use crate::ieee_soft::{fft_complex, try_vec, try_zeros};
use crate::Consts;
use crate::ExactNum;
use crate::ExactNumArray;
use crate::RoundingMode;
use alloc::vec::Vec;
use core::ops::Neg;

/// Default maximum real length for [`dct`] / [`idct`] / [`dst`] / [`idst`] and the window
/// generators (16,384; was 2,048 before 2.0). The `*_with_max_points` transforms take the cap
/// per call instead. Every transform runs a `2N`-point FFT internally; that FFT is sized by
/// this cap (or the per-call one), not by `FFT_MAX_POINTS`.
pub const DSP_MAX_POINTS: usize = 16_384;

fn work_p(p: usize) -> usize {
    p.saturating_add(WORD_BIT_SIZE)
}

fn finite(x: &ExactNum) -> bool {
    !x.is_nan() && !x.is_inf()
}

fn copy_p(x: &ExactNum, p: usize, rm: RoundingMode) -> ExactNum {
    let mut y = x.clone();
    let _ = y.set_precision(p, rm);
    y
}

fn to_row(p: usize, rm: RoundingMode, vals: &[ExactNum]) -> Option<ExactNumArray> {
    let rounded: Vec<ExactNum> = vals.iter().map(|x| copy_p(x, p, rm)).collect();
    ExactNumArray::from_shape(p, 1, rounded.len(), &rounded)
}

/// `k` as an `ExactNum` at `p` bits (exact for every `usize`).
fn num(k: usize, p: usize) -> Option<ExactNum> {
    Some(ExactNum::from_u64(u64::try_from(k).ok()?, p))
}

/// Real vector from a `(1, n)` or `(n, 1)` array whose length is a power of two at most
/// `max_points`. The length is checked before anything is copied or allocated.
fn as_real_capped(signal: &ExactNumArray, max_points: usize, p: usize) -> Option<Vec<ExactNum>> {
    let (rows, cols) = signal.shape();
    if !((rows == 1 && cols > 0) || (cols == 1 && rows > 0)) {
        return None;
    }
    let n = signal.len();
    if !n.is_power_of_two() || n > max_points {
        return None;
    }
    let mut out = try_vec(n)?;
    for v in signal.as_slice() {
        if !finite(v) {
            return None;
        }
        out.push(copy_p(v, p, RoundingMode::None));
    }
    Some(out)
}

/// Forward `2N`-point FFT of the real `x` zero-padded to `2N`.
fn pad_fft(
    mut x: Vec<ExactNum>,
    p: usize,
    cc: &mut Consts,
) -> Option<(Vec<ExactNum>, Vec<ExactNum>)> {
    let n2 = x.len().checked_mul(2)?;
    crate::ieee_soft::check_alloc_bytes(n2.checked_mul(2)?)?;
    let extra = n2 - x.len();
    x.try_reserve_exact(extra).ok()?;
    for _ in 0..extra {
        x.push(ExactNum::from_u8(0, p));
    }
    fft_complex(x, try_zeros(n2, p)?, false, p, RoundingMode::None, cc)
}

/// Discrete cosine transform, type II, via a `2N`-point FFT.
///
/// `X_k = Σ_n x_n cos(π(n+1/2)k/N)`. Input is a real row or column. `N` must
/// be a power of two and at most [`DSP_MAX_POINTS`]. Empty, non-finite, or a
/// bad shape returns `None`.
pub fn dct(
    signal: &ExactNumArray,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNumArray> {
    dct_with_max_points(signal, DSP_MAX_POINTS, p, rm, cc)
}

/// [`dct`] with the length cap `max_points` given per call instead of [`DSP_MAX_POINTS`].
/// `N > max_points` returns `None` before anything is copied; `max_points = 0` rejects every
/// input. The internal `2N`-point FFT follows `N`, so a cap above the default works (subject
/// to `check_alloc_bytes` on the working buffers).
pub fn dct_with_max_points(
    signal: &ExactNumArray,
    max_points: usize,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNumArray> {
    dct_dst(signal, max_points, false, p, rm, cc)
}

/// Discrete sine transform, type II, via a `2N`-point FFT.
///
/// `X_k = Σ_n x_n sin(π(n+1/2)(k+1)/N)`. Same length rules as [`dct`].
pub fn dst(
    signal: &ExactNumArray,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNumArray> {
    dst_with_max_points(signal, DSP_MAX_POINTS, p, rm, cc)
}

/// [`dst`] with a per-call length cap (same rules as [`dct_with_max_points`]).
pub fn dst_with_max_points(
    signal: &ExactNumArray,
    max_points: usize,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNumArray> {
    dct_dst(signal, max_points, true, p, rm, cc)
}

/// Type-II DCT (`sine = false`) or DST (`sine = true`) from the `2N`-point FFT `Y` of the
/// zero-padded input: `X_k = Re(Y_k e^{-iπk/2N})`, and for the DST with `m = k + 1`,
/// `X_k = −Im(Y_m e^{−iπm/2N}) = Re(Y_m) sin θ − Im(Y_m) cos θ`.
fn dct_dst(
    signal: &ExactNumArray,
    max_points: usize,
    sine: bool,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNumArray> {
    let wrk = work_p(p);
    let none = RoundingMode::None;
    let x = as_real_capped(signal, max_points, wrk)?;
    let n = x.len();
    let n2 = n.checked_mul(2)?;
    let (yr, yi) = pad_fft(x, wrk, cc)?;
    let pi = cc.pi(wrk, none);
    let two_n = num(n2, wrk)?;
    let mut out = try_vec(n)?;
    for k in 0..n {
        let m = if sine { k + 1 } else { k };
        let ang = num(m, wrk)?.mul(&pi, wrk, none).div(&two_n, wrk, none);
        let (sn, cs) = ang.sin_cos(wrk, none, cc);
        let (r, i) = (yr.get(m)?, yi.get(m)?);
        let v = if sine {
            sn.mul(r, wrk, none).sub(&cs.mul(i, wrk, none), wrk, none)
        } else {
            cs.mul(r, wrk, none).add(&sn.mul(i, wrk, none), wrk, none)
        };
        out.push(v);
    }
    to_row(p, rm, &out)
}

/// Inverse of [`dct`]: type-III DCT scaled by `2/N` so `idct(dct(x)) = x`.
///
/// `x_n = X_0/N + (2/N) Σ_{k=1}^{N-1} X_k cos(π k (n+1/2)/N)`, computed in `O(N log N)`
/// with a `2N`-point FFT.
pub fn idct(
    signal: &ExactNumArray,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNumArray> {
    idct_with_max_points(signal, DSP_MAX_POINTS, p, rm, cc)
}

/// [`idct`] with a per-call length cap (same rules as [`dct_with_max_points`]).
pub fn idct_with_max_points(
    signal: &ExactNumArray,
    max_points: usize,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNumArray> {
    idct_idst(signal, max_points, false, p, rm, cc)
}

/// Inverse of [`dst`]: type-III DST scaled by `2/N` so `idst(dst(x)) = x`.
///
/// `x_n = (2/N)[ (1/2)(−1)^n X_{N−1} + Σ_{k=0}^{N−2} X_k sin(π(k+1)(n+1/2)/N) ]`, computed in
/// `O(N log N)` with a `2N`-point FFT.
pub fn idst(
    signal: &ExactNumArray,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNumArray> {
    idst_with_max_points(signal, DSP_MAX_POINTS, p, rm, cc)
}

/// [`idst`] with a per-call length cap (same rules as [`dct_with_max_points`]).
pub fn idst_with_max_points(
    signal: &ExactNumArray,
    max_points: usize,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNumArray> {
    idct_idst(signal, max_points, true, p, rm, cc)
}

/// Type-III inverses through one `2N`-point forward FFT.
///
/// DCT-III: with `c_0 = 1/2`, `c_k = 1` and `θ_k = πk/2N`, put `Z_k = c_k X_k e^{−iθ_k}` at
/// `k < N` (zero above). Then `Re(FFT(Z))_n = Σ c_k X_k cos(θ_k (2n+1))`, so
/// `x_n = (2/N) Re(FFT(Z))_n`.
///
/// DST-III: with `m = 1..=N`, `d_m = X_{m−1}` (`d_N = X_{N−1}/2`) and `φ_m = πm/2N`, put
/// `Z_m = d_m e^{−iφ_m}` at index `m`. Then `−Im(FFT(Z))_n = Σ d_m sin(φ_m (2n+1))`, and
/// `sin(πN(2n+1)/2N) = (−1)^n` gives the `X_{N−1}` term, so `x_n = −(2/N) Im(FFT(Z))_n`.
fn idct_idst(
    signal: &ExactNumArray,
    max_points: usize,
    sine: bool,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNumArray> {
    let wrk = work_p(p);
    let none = RoundingMode::None;
    let x = as_real_capped(signal, max_points, wrk)?;
    let n = x.len();
    let n2 = n.checked_mul(2)?;
    crate::ieee_soft::check_alloc_bytes(n2.checked_mul(2)?)?;
    let mut zr = try_zeros(n2, wrk)?;
    let mut zi = try_zeros(n2, wrk)?;
    let pi = cc.pi(wrk, none);
    let two_n = num(n2, wrk)?;
    let half = ExactNum::from_u8(1, wrk).div(&ExactNum::from_u8(2, wrk), wrk, none);
    for (k, xk) in x.iter().enumerate() {
        let idx = if sine { k + 1 } else { k };
        let ang = num(idx, wrk)?.mul(&pi, wrk, none).div(&two_n, wrk, none);
        let (sn, cs) = ang.sin_cos(wrk, none, cc);
        let halve = if sine { idx == n } else { k == 0 };
        let d = if halve {
            ExactNum::mul_ref(&half, xk, wrk, none)
        } else {
            ExactNum::mul_ref(&ExactNum::from_u8(1, wrk), xk, wrk, none)
        };
        *zr.get_mut(idx)? = ExactNum::mul_ref(&cs, &d, wrk, none);
        *zi.get_mut(idx)? = sn.mul(&d, wrk, none).neg();
    }
    let (fr, fi) = fft_complex(zr, zi, false, wrk, none, cc)?;
    let scale = ExactNum::from_u8(2, wrk).div(&num(n, wrk)?, wrk, none);
    let mut out = try_vec(n)?;
    for j in 0..n {
        let v = if sine {
            ExactNum::mul_ref(fi.get(j)?, &scale, wrk, none).neg()
        } else {
            ExactNum::mul_ref(fr.get(j)?, &scale, wrk, none)
        };
        if !finite(&v) {
            return None;
        }
        out.push(v);
    }
    to_row(p, rm, &out)
}

/// Unnormalized DFT of a real row or column. Output is `(2, n)` (row 0 real,
/// row 1 imaginary). Same power-of-two rule and default cap as [`ExactNumArray::fft`].
pub fn fft_real(
    signal: &ExactNumArray,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNumArray> {
    fft_real_with_max_points(signal, crate::FFT_MAX_POINTS, p, rm, cc)
}

/// [`fft_real`] with a per-call length cap (see [`ExactNumArray::fft_with_max_points`]).
pub fn fft_real_with_max_points(
    signal: &ExactNumArray,
    max_points: usize,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNumArray> {
    let (rows, cols) = signal.shape();
    if !((rows == 1 && cols > 0) || (cols == 1 && rows > 0)) {
        return None;
    }
    if signal.as_slice().iter().any(|v| !finite(v)) {
        return None;
    }
    signal.fft_with_max_points(max_points, p, rm, cc)
}

fn window_len_ok(n: usize) -> bool {
    n > 0 && n <= DSP_MAX_POINTS
}

fn frac(num: u32, den: u32, p: usize, rm: RoundingMode) -> ExactNum {
    ExactNum::from_u32(num, p).div(&ExactNum::from_u32(den, p), p, rm)
}

/// Symmetric cosine argument `2π k / (n−1)`. `n ≥ 2`.
fn two_pi_k_over_nm1(k: usize, n: usize, p: usize, rm: RoundingMode, cc: &mut Consts) -> ExactNum {
    let two = ExactNum::from_u8(2, p);
    two.mul(&cc.pi(p, rm), p, rm)
        .mul(&ExactNum::from_u32(k as u32, p), p, rm)
        .div(&ExactNum::from_u32((n - 1) as u32, p), p, rm)
}

fn ones_row(n: usize, p: usize, rm: RoundingMode) -> Option<ExactNumArray> {
    if !window_len_ok(n) {
        return None;
    }
    let mut one = ExactNum::from_u8(1, p);
    let _ = one.set_precision(p, rm);
    Some(ExactNumArray::filled(p, n, &one))
}

/// Symmetric Hann: `½(1 − cos(2πk/(n−1)))`. `n = 1` is `[1]`.
///
/// The plan wrote `2πk/N` (periodic). The locked gold `hann_window(4) =
/// [0, 3/4, 3/4, 0]` is the symmetric form.
pub fn hann_window(n: usize, p: usize, rm: RoundingMode, cc: &mut Consts) -> Option<ExactNumArray> {
    if n == 1 {
        return ones_row(n, p, rm);
    }
    if !window_len_ok(n) {
        return None;
    }
    let wrk = work_p(p);
    let none = RoundingMode::None;
    let half = frac(1, 2, wrk, none);
    let mut vals = Vec::with_capacity(n);
    for k in 0..n {
        let one = ExactNum::from_u8(1, wrk);
        let c = two_pi_k_over_nm1(k, n, wrk, none, cc).cos(wrk, none, cc);
        let w = one.sub(&c, wrk, none).mul(&half, wrk, none);
        if !finite(&w) {
            return None;
        }
        vals.push(w);
    }
    to_row(p, rm, &vals)
}

/// Symmetric Hamming: `0.54 − 0.46 cos(2πk/(n−1))`. Endpoints are `0.08`.
pub fn hamming_window(
    n: usize,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNumArray> {
    if n == 1 {
        return ones_row(n, p, rm);
    }
    if !window_len_ok(n) {
        return None;
    }
    let wrk = work_p(p);
    let none = RoundingMode::None;
    let a1 = frac(23, 50, wrk, none);
    let mut vals = Vec::with_capacity(n);
    for k in 0..n {
        let a0 = frac(27, 50, wrk, none);
        let c = two_pi_k_over_nm1(k, n, wrk, none, cc).cos(wrk, none, cc);
        let w = a0.sub(&c.mul(&a1, wrk, none), wrk, none);
        if !finite(&w) {
            return None;
        }
        vals.push(w);
    }
    to_row(p, rm, &vals)
}

/// Symmetric Blackman: `0.42 − 0.5 cos(θ) + 0.08 cos(2θ)`, `θ = 2πk/(n−1)`.
pub fn blackman_window(
    n: usize,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNumArray> {
    if n == 1 {
        return ones_row(n, p, rm);
    }
    if !window_len_ok(n) {
        return None;
    }
    let wrk = work_p(p);
    let none = RoundingMode::None;
    let a1 = frac(1, 2, wrk, none);
    let a2 = frac(2, 25, wrk, none);
    let two = ExactNum::from_u8(2, wrk);
    let mut vals = Vec::with_capacity(n);
    for k in 0..n {
        let a0 = frac(21, 50, wrk, none);
        let th = two_pi_k_over_nm1(k, n, wrk, none, cc);
        let c1 = th.cos(wrk, none, cc);
        let c2 = th.mul(&two, wrk, none).cos(wrk, none, cc);
        let w = a0
            .sub(&c1.mul(&a1, wrk, none), wrk, none)
            .add(&c2.mul(&a2, wrk, none), wrk, none);
        if !finite(&w) {
            return None;
        }
        vals.push(w);
    }
    to_row(p, rm, &vals)
}

/// Kaiser–Bessel: `I_0(β √(1−t_k²)) / I_0(β)` with `t_k = (k−(n−1)/2)/((n−1)/2)`.
///
/// `beta = 0` is the rectangular window. `beta < 0` or non-finite is `None`.
pub fn kaiser_window(
    n: usize,
    beta: &ExactNum,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNumArray> {
    if !window_len_ok(n) || !finite(beta) || beta.is_negative() {
        return None;
    }
    if n == 1 || beta.is_zero() {
        return ones_row(n, p, rm);
    }
    let wrk = work_p(p);
    let none = RoundingMode::None;
    let b = copy_p(beta, wrk, none);
    let nu0 = ExactNum::from_u8(0, wrk);
    let i0b = b.bessel_i(&nu0, wrk, none, cc);
    if !finite(&i0b) || i0b.is_zero() {
        return None;
    }
    let two = ExactNum::from_u8(2, wrk);
    let mid = ExactNum::from_u32((n - 1) as u32, wrk).div(&two, wrk, none);
    let mut vals = Vec::with_capacity(n);
    for k in 0..n {
        let one = ExactNum::from_u8(1, wrk);
        let t = ExactNum::from_u32(k as u32, wrk)
            .sub(&mid, wrk, none)
            .div(&mid, wrk, none);
        let rad = one.sub(&t.sqr(wrk, none), wrk, none);
        if rad.is_negative() {
            return None;
        }
        let arg = rad.sqrt(wrk, none).mul(&b, wrk, none);
        let w = arg.bessel_i(&nu0, wrk, none, cc).div(&i0b, wrk, none);
        if !finite(&w) {
            return None;
        }
        vals.push(w);
    }
    to_row(p, rm, &vals)
}

/// Rectangular window: `n` ones.
pub fn rectangular_window(n: usize, p: usize, rm: RoundingMode) -> Option<ExactNumArray> {
    ones_row(n, p, rm)
}

/// Inverse of [`fft_real`]: [`ExactNumArray::ifft`] then the real row.
///
/// Input must be a `(2, n)` spectrum. Imaginary residuals are dropped.
pub fn ifft_real(
    spectrum: &ExactNumArray,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNumArray> {
    ifft_real_with_max_points(spectrum, crate::FFT_MAX_POINTS, p, rm, cc)
}

/// [`ifft_real`] with a per-call length cap (see [`ExactNumArray::ifft_with_max_points`]).
pub fn ifft_real_with_max_points(
    spectrum: &ExactNumArray,
    max_points: usize,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNumArray> {
    let (rows, cols) = spectrum.shape();
    if rows != 2 || cols == 0 {
        return None;
    }
    let spec = spectrum.ifft_with_max_points(max_points, p, rm, cc)?;
    // Row-major `(2, cols)`: the first `cols` values are the real row; the imaginary
    // residuals are dropped without copying.
    let mut re = spec.into_vals();
    re.truncate(cols);
    to_row(p, rm, &re)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Consts;

    fn gold_p() -> (usize, RoundingMode) {
        (256, RoundingMode::ToEven)
    }

    fn near(a: &ExactNum, b: &ExactNum, p: usize) -> bool {
        let d = ExactNum::sub_ref(a, b, p, RoundingMode::None).abs();
        d.is_zero() || d.exponent().unwrap_or(0) < -((p as crate::Exponent) - 40)
    }

    fn n_at(p: usize, k: u8) -> ExactNum {
        ExactNum::from_u8(k, p)
    }

    #[test]
    fn dsp_dct_dst_fft_real() -> Result<(), crate::Error> {
        let (p, rm) = gold_p();
        let mut cc = Consts::new()?;
        let n = |k: u8| n_at(p, k);
        let zero = n(0);

        let x = ExactNumArray::from_values(p, &[n(1), n(2), n(3), n(4)]);
        let cx = dct(&x, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        let back = idct(&cx, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        assert_eq!(back.len(), 4);
        for j in 0..4 {
            assert!(
                near(
                    back.get(j).ok_or(crate::Error::InvalidArgument)?,
                    x.get(j).ok_or(crate::Error::InvalidArgument)?,
                    p
                ),
                "idct(dct(x))[{j}]"
            );
        }

        let cst = ExactNumArray::from_values(p, &[n(3), n(3), n(3), n(3)]);
        let dc = dct(&cst, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        let twelve = ExactNum::from_u8(12, p);
        assert!(
            near(dc.get(0).ok_or(crate::Error::InvalidArgument)?, &twelve, p),
            "DC bin"
        );
        for j in 1..4 {
            assert!(
                near(dc.get(j).ok_or(crate::Error::InvalidArgument)?, &zero, p),
                "const bin {j}"
            );
        }

        let sx = dst(&x, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        let sback = idst(&sx, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        for j in 0..4 {
            assert!(
                near(
                    sback.get(j).ok_or(crate::Error::InvalidArgument)?,
                    x.get(j).ok_or(crate::Error::InvalidArgument)?,
                    p
                ),
                "idst(dst(x))[{j}]"
            );
        }

        let n8 = ExactNum::from_u8(8, p);
        let two_pi = n(2).mul(&cc.pi(p, rm), p, rm);
        let mut cos_vals = Vec::with_capacity(8);
        for k in 0..8u8 {
            let kn = ExactNum::from_u8(k, p);
            let ang = kn.mul(&two_pi, p, rm).div(&n8, p, rm);
            cos_vals.push(ang.cos(p, rm, &mut cc));
        }
        let cosine = ExactNumArray::from_values(p, &cos_vals);
        let cspec = fft_real(&cosine, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        assert_eq!(cspec.shape(), (2, 8));
        let four = n(4);
        for j in 0..8 {
            let re = cspec.get2(0, j).ok_or(crate::Error::InvalidArgument)?;
            let im = cspec.get2(1, j).ok_or(crate::Error::InvalidArgument)?;
            if j == 1 || j == 7 {
                assert!(near(re, &four, p), "cos bin {j} re");
            } else {
                assert!(near(re, &zero, p), "cos bin {j} re");
            }
            assert!(near(im, &zero, p), "cos bin {j} im");
        }

        let rec = ifft_real(&cspec, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        for j in 0..8 {
            assert!(near(
                rec.get(j).ok_or(crate::Error::InvalidArgument)?,
                cosine.get(j).ok_or(crate::Error::InvalidArgument)?,
                p
            ));
        }

        let mut e_t = ExactNum::from_u8(0, p);
        let mut e_f = ExactNum::from_u8(0, p);
        for j in 0..8 {
            let xv = cosine.get(j).ok_or(crate::Error::InvalidArgument)?;
            e_t = e_t.add(&xv.sqr(p, rm), p, rm);
            let xr = cspec.get2(0, j).ok_or(crate::Error::InvalidArgument)?;
            let xi = cspec.get2(1, j).ok_or(crate::Error::InvalidArgument)?;
            e_f = e_f.add(&xr.sqr(p, rm), p, rm).add(&xi.sqr(p, rm), p, rm);
        }
        let parseval = e_f.div(&n8, p, rm);
        assert!(near(&parseval, &e_t, p));

        assert!(dct(
            &ExactNumArray::from_values(p, &[n(1), n(2), n(3)]),
            p,
            rm,
            &mut cc
        )
        .is_none());
        assert!(fft_real(&cspec, p, rm, &mut cc).is_none());
        Ok(())
    }

    fn window_sum_pos(w: &ExactNumArray, p: usize) -> Result<bool, crate::Error> {
        let mut s = ExactNum::from_u8(0, p);
        for i in 0..w.len() {
            s = s.add(
                w.get(i).ok_or(crate::Error::InvalidArgument)?,
                p,
                RoundingMode::ToEven,
            );
        }
        Ok(s.is_positive() && !s.is_zero())
    }

    #[test]
    fn dsp_windows() -> Result<(), crate::Error> {
        let (p, rm) = gold_p();
        let mut cc = Consts::new()?;
        let zero = ExactNum::from_u8(0, p);
        let three_fourths = ExactNum::from_u8(3, p).div(&ExactNum::from_u8(4, p), p, rm);
        let eight_hundredths = ExactNum::from_u8(2, p).div(&ExactNum::from_u8(25, p), p, rm);

        let hann = hann_window(4, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        assert_eq!(hann.len(), 4);
        assert!(near(
            hann.get(0).ok_or(crate::Error::InvalidArgument)?,
            &zero,
            p
        ));
        assert!(near(
            hann.get(1).ok_or(crate::Error::InvalidArgument)?,
            &three_fourths,
            p
        ));
        assert!(near(
            hann.get(2).ok_or(crate::Error::InvalidArgument)?,
            &three_fourths,
            p
        ));
        assert!(near(
            hann.get(3).ok_or(crate::Error::InvalidArgument)?,
            &zero,
            p
        ));

        let hamm = hamming_window(4, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        assert!(near(
            hamm.get(0).ok_or(crate::Error::InvalidArgument)?,
            &eight_hundredths,
            p
        ));
        assert!(near(
            hamm.get(3).ok_or(crate::Error::InvalidArgument)?,
            &eight_hundredths,
            p
        ));
        assert!(!near(
            hamm.get(0).ok_or(crate::Error::InvalidArgument)?,
            &zero,
            p
        ));

        let rect = rectangular_window(8, p, rm).ok_or(crate::Error::InvalidArgument)?;
        let k0 = kaiser_window(8, &zero, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        assert_eq!(k0.len(), 8);
        for j in 0..8 {
            assert!(near(
                k0.get(j).ok_or(crate::Error::InvalidArgument)?,
                rect.get(j).ok_or(crate::Error::InvalidArgument)?,
                p
            ));
        }

        let blk = blackman_window(8, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        assert!(window_sum_pos(&hann, p)?);
        assert!(window_sum_pos(&hamm, p)?);
        assert!(window_sum_pos(&blk, p)?);
        assert!(window_sum_pos(&rect, p)?);
        assert!(window_sum_pos(&k0, p)?);

        assert!(hann_window(0, p, rm, &mut cc).is_none());
        assert!(rectangular_window(0, p, rm).is_none());
        assert!(kaiser_window(4, &ExactNum::from_i64(-1, p), p, rm, &mut cc).is_none());
        Ok(())
    }

    fn ramp(n: usize, p: usize) -> ExactNumArray {
        // x_j = (j mod 7) − 3: small integers, not symmetric, exact at any p.
        let mut v = Vec::with_capacity(n);
        for j in 0..n {
            let m = u8::try_from(j % 7).unwrap_or(0);
            v.push(ExactNum::from_u8(m, p).sub(&ExactNum::from_u8(3, p), p, RoundingMode::None));
        }
        ExactNumArray::from_values(p, &v)
    }

    fn roundtrip_ok(x: &ExactNumArray, y: &ExactNumArray, p: usize) -> Result<bool, crate::Error> {
        if x.len() != y.len() {
            return Ok(false);
        }
        for j in 0..x.len() {
            let a = x.get(j).ok_or(crate::Error::InvalidArgument)?;
            let b = y.get(j).ok_or(crate::Error::InvalidArgument)?;
            if !near(a, b, p) {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// Slice 2 caps: default 16,384 (was 2,048); per-call caps tighter and looser.
    #[test]
    fn dsp_default_and_per_call_caps() -> Result<(), crate::Error> {
        let p = 64;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        assert_eq!(DSP_MAX_POINTS, 16_384);
        // Default boundary: 32,768 is refused by every transform (16,385 is not a power of two).
        let over = ramp(32_768, p);
        assert!(dct(&over, p, rm, &mut cc).is_none());
        assert!(idct(&over, p, rm, &mut cc).is_none());
        assert!(dst(&over, p, rm, &mut cc).is_none());
        assert!(idst(&over, p, rm, &mut cc).is_none());
        // Tighter per-call cap.
        let x128 = ramp(128, p);
        let x64 = ramp(64, p);
        assert!(dct_with_max_points(&x128, 64, p, rm, &mut cc).is_none());
        assert!(idct_with_max_points(&x128, 64, p, rm, &mut cc).is_none());
        assert!(dst_with_max_points(&x128, 64, p, rm, &mut cc).is_none());
        assert!(idst_with_max_points(&x128, 64, p, rm, &mut cc).is_none());
        assert!(fft_real_with_max_points(&x128, 64, p, rm, &mut cc).is_none());
        assert!(dct_with_max_points(&x64, 0, p, rm, &mut cc).is_none());
        let c =
            dct_with_max_points(&x64, 64, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        let b =
            idct_with_max_points(&c, 64, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        assert!(roundtrip_ok(&x64, &b, p)?);
        let s =
            dst_with_max_points(&x64, 64, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        let b =
            idst_with_max_points(&s, 64, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        assert!(roundtrip_ok(&x64, &b, p)?);
        let f = fft_real_with_max_points(&x64, 64, p, rm, &mut cc)
            .ok_or(crate::Error::InvalidArgument)?;
        assert!(ifft_real_with_max_points(&f, 32, p, rm, &mut cc).is_none());
        let b = ifft_real_with_max_points(&f, 64, p, rm, &mut cc)
            .ok_or(crate::Error::InvalidArgument)?;
        assert!(roundtrip_ok(&x64, &b, p)?);
        // Looser per-call cap: N = 32,768 past the default (internal FFT of 65,536 points).
        let c = dct_with_max_points(&over, 32_768, p, rm, &mut cc)
            .ok_or(crate::Error::InvalidArgument)?;
        assert_eq!(c.len(), 32_768);
        // The DC bin is Σ x_j: 32,768 = 7·4681 + 1, so Σ ((j mod 7) − 3) = 4681·0 + (0 − 3) = −3.
        let m3 = ExactNum::from_u8(3, p).neg();
        assert!(near_mag(
            c.get(0).ok_or(crate::Error::InvalidArgument)?,
            &m3,
            p,
            16
        ));
        Ok(())
    }

    fn near_mag(a: &ExactNum, b: &ExactNum, p: usize, mag: crate::Exponent) -> bool {
        let d = ExactNum::sub_ref(a, b, p, RoundingMode::None).abs();
        d.is_zero() || d.exponent().unwrap_or(0) < mag - ((p as crate::Exponent) - 40)
    }

    /// Past the old 2,048 cap at N = 4,096 (p = 128): DCT-II of cos(π(n+½)k0/N) is N/2 at k0
    /// and 0 elsewhere (orthogonality of the DCT-II basis); idct and idst invert dct and dst.
    #[test]
    fn dsp_accuracy_above_old_cap() -> Result<(), crate::Error> {
        let p = 128;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let n = 4096usize;
        let nn = ExactNum::from_u32(4096, p);
        let pi = cc.pi(p, rm);
        let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, rm);
        let mut xs = Vec::with_capacity(n);
        for j in 0..4096u32 {
            let ang = ExactNum::from_u32(j, p)
                .add(&half, p, rm)
                .mul(&ExactNum::from_u8(9, p), p, rm)
                .mul(&pi, p, rm)
                .div(&nn, p, rm);
            xs.push(ang.cos(p, rm, &mut cc));
        }
        let x = ExactNumArray::from_values(p, &xs);
        let c = dct(&x, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        let half_n = ExactNum::from_u32(2048, p);
        let zero = ExactNum::from_u8(0, p);
        for k in 0..n {
            let want = if k == 9 { &half_n } else { &zero };
            assert!(
                near_mag(c.get(k).ok_or(crate::Error::InvalidArgument)?, want, p, 12),
                "bin {k}"
            );
        }
        let back = idct(&c, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        assert!(roundtrip_ok(&x, &back, p)?);
        let r = ramp(n, p);
        let s = dst(&r, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        let back = idst(&s, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        assert!(roundtrip_ok(&r, &back, p)?);
        Ok(())
    }

    /// The `O(N log N)` inverses agree with the defining sums (the pre-2.0 direct formulas)
    /// at N = 16 for an arbitrary input.
    #[test]
    fn dsp_fast_inverses_match_direct_sums() -> Result<(), crate::Error> {
        let p = 128;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let n = 16usize;
        let x = ramp(n, p);
        let ci = idct(&x, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        let si = idst(&x, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        let nn = ExactNum::from_u8(16, p);
        let pi = cc.pi(p, rm);
        let half = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(2, p), p, rm);
        let two_over_n = ExactNum::from_u8(2, p).div(&nn, p, rm);
        for j in 0..16u8 {
            let jh = ExactNum::from_u8(j, p).add(&half, p, rm);
            let x0 = x.get(0).ok_or(crate::Error::InvalidArgument)?;
            let mut acc_c = ExactNum::mul_ref(x0, &half, p, rm);
            let xl = x.get(n - 1).ok_or(crate::Error::InvalidArgument)?;
            let last = ExactNum::mul_ref(xl, &half, p, rm);
            let mut acc_s = if j % 2 == 0 { last } else { last.neg() };
            for k in 0..16u8 {
                let xk = x.get(usize::from(k)).ok_or(crate::Error::InvalidArgument)?;
                if k > 0 {
                    let a = ExactNum::from_u8(k, p)
                        .mul(&pi, p, rm)
                        .mul(&jh, p, rm)
                        .div(&nn, p, rm);
                    acc_c = acc_c.add(&a.cos(p, rm, &mut cc).mul(xk, p, rm), p, rm);
                }
                if k < 15 {
                    let a = ExactNum::from_u8(k + 1, p)
                        .mul(&pi, p, rm)
                        .mul(&jh, p, rm)
                        .div(&nn, p, rm);
                    acc_s = acc_s.add(&a.sin(p, rm, &mut cc).mul(xk, p, rm), p, rm);
                }
            }
            let jc = usize::from(j);
            let want_c = acc_c.mul(&two_over_n, p, rm);
            let want_s = acc_s.mul(&two_over_n, p, rm);
            assert!(
                near(ci.get(jc).ok_or(crate::Error::InvalidArgument)?, &want_c, p),
                "idct {j}"
            );
            assert!(
                near(si.get(jc).ok_or(crate::Error::InvalidArgument)?, &want_s, p),
                "idst {j}"
            );
        }
        Ok(())
    }

    /// The full default size, N = 16,384 at p = 64: dct/idct and dst/idst round trips.
    #[test]
    fn dsp_at_default_cap_16384() -> Result<(), crate::Error> {
        let p = 64;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new()?;
        let x = ramp(16_384, p);
        let c = dct(&x, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        let b = idct(&c, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        assert!(roundtrip_ok(&x, &b, p)?);
        let s = dst(&x, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        let b = idst(&s, p, rm, &mut cc).ok_or(crate::Error::InvalidArgument)?;
        assert!(roundtrip_ok(&x, &b, p)?);
        Ok(())
    }
}
