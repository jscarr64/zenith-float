//! Discrete cosine / sine transforms, real-input FFT wrappers, and windows.

use crate::defs::WORD_BIT_SIZE;
use crate::Consts;
use crate::ExactNum;
use crate::ExactNumArray;
use crate::RoundingMode;
use alloc::vec::Vec;

/// Maximum real length for [`dct`] / [`idct`] / [`dst`] / [`idst`] and the
/// window generators. Transforms use a `2N`-point FFT, so this is half of
/// `FFT_MAX_POINTS`.
pub const DSP_MAX_POINTS: usize = 2048;

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

/// Real vector from a `(1, n)` or `(n, 1)` array.
fn as_real(signal: &ExactNumArray, p: usize, rm: RoundingMode) -> Option<Vec<ExactNum>> {
    let (rows, cols) = signal.shape();
    if !((rows == 1 && cols > 0) || (cols == 1 && rows > 0)) {
        return None;
    }
    let mut out = Vec::with_capacity(signal.len());
    for v in signal.as_slice() {
        if !finite(v) {
            return None;
        }
        out.push(copy_p(v, p, rm));
    }
    Some(out)
}

fn real_len_ok(n: usize) -> bool {
    n > 0 && n.is_power_of_two() && n <= DSP_MAX_POINTS
}

fn pad_fft(
    x: &[ExactNum],
    n2: usize,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNumArray> {
    let zero = ExactNum::from_u8(0, p);
    let mut y = x.to_vec();
    y.resize(n2, zero);
    ExactNumArray::from_values(p, &y).fft(p, rm, cc)
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
    let wrk = work_p(p);
    let x = as_real(signal, wrk, RoundingMode::None)?;
    let n = x.len();
    let n2 = n.checked_mul(2)?;
    if !real_len_ok(n) {
        return None;
    }
    let spec = pad_fft(&x, n2, wrk, RoundingMode::None, cc)?;
    let pi = cc.pi(wrk, RoundingMode::None);
    let two_n = ExactNum::from_u32(n2 as u32, wrk);
    let mut out = Vec::with_capacity(n);
    for k in 0..n {
        let ang = pi
            .mul(&ExactNum::from_u32(k as u32, wrk), wrk, RoundingMode::None)
            .div(&two_n, wrk, RoundingMode::None);
        let (sn, cs) = ang.sin_cos(wrk, RoundingMode::None, cc);
        let yr = spec.get2(0, k)?;
        let yi = spec.get2(1, k)?;
        let re = yr.mul(&cs, wrk, RoundingMode::None).add(
            &yi.mul(&sn, wrk, RoundingMode::None),
            wrk,
            RoundingMode::None,
        );
        out.push(re);
    }
    to_row(p, rm, &out)
}

/// Inverse of [`dct`]: type-III DCT scaled by `2/N` so `idct(dct(x)) = x`.
///
/// `x_n = X_0/N + (2/N) Σ_{k=1}^{N-1} X_k cos(π k (n+1/2)/N)`.
pub fn idct(
    signal: &ExactNumArray,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNumArray> {
    let wrk = work_p(p);
    let x = as_real(signal, wrk, RoundingMode::None)?;
    let n = x.len();
    if !real_len_ok(n) {
        return None;
    }
    let none = RoundingMode::None;
    let pi = cc.pi(wrk, none);
    let nn = ExactNum::from_u32(n as u32, wrk);
    let two = ExactNum::from_u8(2, wrk);
    let inv_n = ExactNum::from_u8(1, wrk).div(&nn, wrk, none);
    let two_n = two.div(&nn, wrk, none);
    let half = ExactNum::from_u8(1, wrk).div(&two, wrk, none);
    let mut out = Vec::with_capacity(n);
    for ni in 0..n {
        let n_half = ExactNum::from_u32(ni as u32, wrk).add(&half, wrk, none);
        let mut acc = x[0].mul(&inv_n, wrk, none);
        for (k, xk) in x.iter().enumerate().take(n).skip(1) {
            let ang = pi
                .mul(&ExactNum::from_u32(k as u32, wrk), wrk, none)
                .mul(&n_half, wrk, none)
                .div(&nn, wrk, none);
            let c = ang.cos(wrk, none, cc);
            acc = acc.add(&xk.mul(&c, wrk, none).mul(&two_n, wrk, none), wrk, none);
        }
        if !finite(&acc) {
            return None;
        }
        out.push(acc);
    }
    to_row(p, rm, &out)
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
    let wrk = work_p(p);
    let x = as_real(signal, wrk, RoundingMode::None)?;
    let n = x.len();
    let n2 = n.checked_mul(2)?;
    if !real_len_ok(n) {
        return None;
    }
    let spec = pad_fft(&x, n2, wrk, RoundingMode::None, cc)?;
    let pi = cc.pi(wrk, RoundingMode::None);
    let two_n = ExactNum::from_u32(n2 as u32, wrk);
    let mut out = Vec::with_capacity(n);
    for k in 0..n {
        let m = k + 1;
        let ang = pi
            .mul(&ExactNum::from_u32(m as u32, wrk), wrk, RoundingMode::None)
            .div(&two_n, wrk, RoundingMode::None);
        let (sn, cs) = ang.sin_cos(wrk, RoundingMode::None, cc);
        let yr = spec.get2(0, m)?;
        let yi = spec.get2(1, m)?;
        // −Im(Y[m] exp(−iθ)) = yr sinθ − yi cosθ
        let val = yr.mul(&sn, wrk, RoundingMode::None).sub(
            &yi.mul(&cs, wrk, RoundingMode::None),
            wrk,
            RoundingMode::None,
        );
        out.push(val);
    }
    to_row(p, rm, &out)
}

/// Inverse of [`dst`]: type-III DST scaled by `2/N` so `idst(dst(x)) = x`.
///
/// `x_n = (2/N)[ (1/2)(−1)^n X_{N−1} + Σ_{k=0}^{N−2} X_k sin(π(k+1)(n+1/2)/N) ]`.
pub fn idst(
    signal: &ExactNumArray,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNumArray> {
    let wrk = work_p(p);
    let x = as_real(signal, wrk, RoundingMode::None)?;
    let n = x.len();
    if !real_len_ok(n) {
        return None;
    }
    let none = RoundingMode::None;
    let pi = cc.pi(wrk, none);
    let nn = ExactNum::from_u32(n as u32, wrk);
    let two = ExactNum::from_u8(2, wrk);
    let two_n = two.div(&nn, wrk, none);
    let half = ExactNum::from_u8(1, wrk).div(&two, wrk, none);
    let last = x[n - 1].mul(&half, wrk, none);
    let mut out = Vec::with_capacity(n);
    for ni in 0..n {
        let n_half = ExactNum::from_u32(ni as u32, wrk).add(&half, wrk, none);
        let mut acc = if ni % 2 == 0 { last.clone() } else { last.neg() };
        for (k, xk) in x.iter().enumerate().take(n - 1) {
            let ang = pi
                .mul(&ExactNum::from_u32((k + 1) as u32, wrk), wrk, none)
                .mul(&n_half, wrk, none)
                .div(&nn, wrk, none);
            let s = ang.sin(wrk, none, cc);
            acc = acc.add(&xk.mul(&s, wrk, none), wrk, none);
        }
        acc = acc.mul(&two_n, wrk, none);
        if !finite(&acc) {
            return None;
        }
        out.push(acc);
    }
    to_row(p, rm, &out)
}

/// Unnormalized DFT of a real row or column. Output is `(2, n)` (row 0 real,
/// row 1 imaginary). Same power-of-two rule as [`ExactNumArray::fft`].
pub fn fft_real(
    signal: &ExactNumArray,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Option<ExactNumArray> {
    let _ = as_real(signal, p, rm)?;
    signal.fft(p, rm, cc)
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
    let one = ExactNum::from_u8(1, wrk);
    let mut vals = Vec::with_capacity(n);
    for k in 0..n {
        let c = two_pi_k_over_nm1(k, n, wrk, none, cc).cos(wrk, none, cc);
        let w = half.mul(&one.sub(&c, wrk, none), wrk, none);
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
    let a0 = frac(27, 50, wrk, none);
    let a1 = frac(23, 50, wrk, none);
    let mut vals = Vec::with_capacity(n);
    for k in 0..n {
        let c = two_pi_k_over_nm1(k, n, wrk, none, cc).cos(wrk, none, cc);
        let w = a0.sub(&a1.mul(&c, wrk, none), wrk, none);
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
    let a0 = frac(21, 50, wrk, none);
    let a1 = frac(1, 2, wrk, none);
    let a2 = frac(2, 25, wrk, none);
    let two = ExactNum::from_u8(2, wrk);
    let mut vals = Vec::with_capacity(n);
    for k in 0..n {
        let th = two_pi_k_over_nm1(k, n, wrk, none, cc);
        let c1 = th.cos(wrk, none, cc);
        let c2 = th.mul(&two, wrk, none).cos(wrk, none, cc);
        let w = a0
            .sub(&a1.mul(&c1, wrk, none), wrk, none)
            .add(&a2.mul(&c2, wrk, none), wrk, none);
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
    let one = ExactNum::from_u8(1, wrk);
    let mut vals = Vec::with_capacity(n);
    for k in 0..n {
        let t = ExactNum::from_u32(k as u32, wrk)
            .sub(&mid, wrk, none)
            .div(&mid, wrk, none);
        let rad = one.sub(&t.mul(&t, wrk, none), wrk, none);
        if rad.is_negative() {
            return None;
        }
        let arg = b.mul(&rad.sqrt(wrk, none), wrk, none);
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
    let (rows, cols) = spectrum.shape();
    if rows != 2 || cols == 0 {
        return None;
    }
    let spec = spectrum.ifft(p, rm, cc)?;
    let mut re = Vec::with_capacity(cols);
    for j in 0..cols {
        re.push(spec.get2(0, j)?.clone());
    }
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
        let d = a.sub(b, p, RoundingMode::None).abs();
        d.is_zero() || d.exponent().unwrap_or(0) < -((p as i32) - 40)
    }

    fn n_at(p: usize, k: u8) -> ExactNum {
        ExactNum::from_u8(k, p)
    }

    #[test]
    fn dsp_dct_dst_fft_real() {
        let (p, rm) = gold_p();
        let mut cc = Consts::new().expect("consts");
        let n = |k: u8| n_at(p, k);
        let zero = n(0);

        let x = ExactNumArray::from_values(p, &[n(1), n(2), n(3), n(4)]);
        let cx = dct(&x, p, rm, &mut cc).expect("dct");
        let back = idct(&cx, p, rm, &mut cc).expect("idct");
        assert_eq!(back.len(), 4);
        for j in 0..4 {
            assert!(
                near(back.get(j).unwrap(), x.get(j).unwrap(), p),
                "idct(dct(x))[{j}]"
            );
        }

        let cst = ExactNumArray::from_values(p, &[n(3), n(3), n(3), n(3)]);
        let dc = dct(&cst, p, rm, &mut cc).expect("dct const");
        let twelve = ExactNum::from_u8(12, p);
        assert!(near(dc.get(0).unwrap(), &twelve, p), "DC bin");
        for j in 1..4 {
            assert!(near(dc.get(j).unwrap(), &zero, p), "const bin {j}");
        }

        let sx = dst(&x, p, rm, &mut cc).expect("dst");
        let sback = idst(&sx, p, rm, &mut cc).expect("idst");
        for j in 0..4 {
            assert!(
                near(sback.get(j).unwrap(), x.get(j).unwrap(), p),
                "idst(dst(x))[{j}]"
            );
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
        let cspec = fft_real(&cosine, p, rm, &mut cc).expect("fft_real cos");
        assert_eq!(cspec.shape(), (2, 8));
        let four = n(4);
        for j in 0..8 {
            let re = cspec.get2(0, j).unwrap();
            let im = cspec.get2(1, j).unwrap();
            if j == 1 || j == 7 {
                assert!(near(re, &four, p), "cos bin {j} re");
            } else {
                assert!(near(re, &zero, p), "cos bin {j} re");
            }
            assert!(near(im, &zero, p), "cos bin {j} im");
        }

        let rec = ifft_real(&cspec, p, rm, &mut cc).expect("ifft_real");
        for j in 0..8 {
            assert!(near(rec.get(j).unwrap(), cosine.get(j).unwrap(), p));
        }

        let mut e_t = ExactNum::from_u8(0, p);
        let mut e_f = ExactNum::from_u8(0, p);
        for j in 0..8 {
            let xv = cosine.get(j).unwrap();
            e_t = e_t.add(&xv.mul(xv, p, rm), p, rm);
            let xr = cspec.get2(0, j).unwrap();
            let xi = cspec.get2(1, j).unwrap();
            e_f = e_f
                .add(&xr.mul(xr, p, rm), p, rm)
                .add(&xi.mul(xi, p, rm), p, rm);
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
    }

    fn window_sum_pos(w: &ExactNumArray, p: usize) -> bool {
        let mut s = ExactNum::from_u8(0, p);
        for i in 0..w.len() {
            s = s.add(w.get(i).unwrap(), p, RoundingMode::ToEven);
        }
        s.is_positive() && !s.is_zero()
    }

    #[test]
    fn dsp_windows() {
        let (p, rm) = gold_p();
        let mut cc = Consts::new().expect("consts");
        let zero = ExactNum::from_u8(0, p);
        let three_fourths = ExactNum::from_u8(3, p).div(&ExactNum::from_u8(4, p), p, rm);
        let eight_hundredths = ExactNum::from_u8(2, p).div(&ExactNum::from_u8(25, p), p, rm);

        let hann = hann_window(4, p, rm, &mut cc).expect("hann");
        assert_eq!(hann.len(), 4);
        assert!(near(hann.get(0).unwrap(), &zero, p));
        assert!(near(hann.get(1).unwrap(), &three_fourths, p));
        assert!(near(hann.get(2).unwrap(), &three_fourths, p));
        assert!(near(hann.get(3).unwrap(), &zero, p));

        let hamm = hamming_window(4, p, rm, &mut cc).expect("hamming");
        assert!(near(hamm.get(0).unwrap(), &eight_hundredths, p));
        assert!(near(hamm.get(3).unwrap(), &eight_hundredths, p));
        assert!(!near(hamm.get(0).unwrap(), &zero, p));

        let rect = rectangular_window(8, p, rm).expect("rect");
        let k0 = kaiser_window(8, &zero, p, rm, &mut cc).expect("kaiser0");
        assert_eq!(k0.len(), 8);
        for j in 0..8 {
            assert!(near(k0.get(j).unwrap(), rect.get(j).unwrap(), p));
        }

        let blk = blackman_window(8, p, rm, &mut cc).expect("blackman");
        assert!(window_sum_pos(&hann, p));
        assert!(window_sum_pos(&hamm, p));
        assert!(window_sum_pos(&blk, p));
        assert!(window_sum_pos(&rect, p));
        assert!(window_sum_pos(&k0, p));

        assert!(hann_window(0, p, rm, &mut cc).is_none());
        assert!(rectangular_window(0, p, rm).is_none());
        assert!(kaiser_window(4, &ExactNum::from_i64(-1, p), p, rm, &mut cc).is_none());
    }
}
