//! Discrete cosine / sine transforms and real-input FFT wrappers.

use crate::defs::WORD_BIT_SIZE;
use crate::Consts;
use crate::ExactNum;
use crate::ExactNumArray;
use crate::RoundingMode;
use alloc::vec::Vec;

/// Maximum real length for [`dct`] / [`idct`] / [`dst`] / [`idst`].
/// Those transforms use a `2N`-point FFT, so this is half of `FFT_MAX_POINTS`.
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
        for k in 1..n {
            let ang = pi
                .mul(&ExactNum::from_u32(k as u32, wrk), wrk, none)
                .mul(&n_half, wrk, none)
                .div(&nn, wrk, none);
            let c = ang.cos(wrk, none, cc);
            acc = acc.add(&x[k].mul(&c, wrk, none).mul(&two_n, wrk, none), wrk, none);
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
        for k in 0..n - 1 {
            let ang = pi
                .mul(&ExactNum::from_u32((k + 1) as u32, wrk), wrk, none)
                .mul(&n_half, wrk, none)
                .div(&nn, wrk, none);
            let s = ang.sin(wrk, none, cc);
            acc = acc.add(&x[k].mul(&s, wrk, none), wrk, none);
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
}
