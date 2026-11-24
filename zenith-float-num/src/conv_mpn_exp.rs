// Derived from MPFR 4.2.2 `mpn_exp.c` / `get_str.c` (LGPL-3.0-or-later).
//! Truncated `b^e` with a proven ulp bound, after MPFR `mpfr_mpn_exp`.

use crate::common::buf::WordBuf;
use crate::defs::DoubleWord;
use crate::defs::Error;
use crate::defs::Sign;
use crate::defs::Word;
use crate::defs::WORD_BIT_SIZE;
use crate::defs::WORD_SIGNIFICANT_BIT;
use crate::mantissa::Mantissa;

/// A truncated power `a * 2^{exp_r}` of an integer base.
#[derive(Debug)]
pub(crate) struct BasePow {
    /// Normalized `n`-limb integer `a` (MSB set).
    pub a: WordBuf,
    /// `a * 2^{exp_r}` approximates `b^e`.
    pub exp_r: i128,
    /// `None` if `a * 2^{exp_r} == b^e`. `Some(f)` means
    /// `a*2^exp_r <= b^e <= 2^exp_r * (a + 2^f)` (at most `2^f` ulps),
    /// as in `mpfr_mpn_exp`.
    pub err_ulps: Option<u32>,
}

/// Approximation of `x * b^k` as `a * 2^{exp_a}` with a proven error.
#[derive(Debug)]
pub(crate) struct ScaledY {
    pub a: WordBuf,
    pub exp_a: i128,
    /// `None` if exact. `Some(f)` means `|Y - a*2^{exp_a}| < 2^{exp_a + f}`.
    pub err_ulps: Option<u32>,
}

/// Lowest set bit index in a little-endian limb buffer, or `None` if all zero.
pub(crate) fn mpn_scan1(a: &[Word]) -> Option<usize> {
    for (i, &w) in a.iter().enumerate() {
        if w != 0 {
            return Some(i * WORD_BIT_SIZE + w.trailing_zeros() as usize);
        }
    }
    None
}

/// Schoolbook `out = a * b`. `out` must be at least `a.len() + b.len()` limbs, else `Err(InvalidArgument)`.
pub(crate) fn mpn_mul(a: &[Word], b: &[Word], out: &mut [Word]) -> Result<(), Error> {
    if out.len() < a.len() + b.len() {
        return Err(Error::InvalidArgument);
    }
    out.fill(0);
    for (i, &ai) in a.iter().enumerate() {
        if ai == 0 {
            continue;
        }
        let mut k: DoubleWord = 0;
        let ai = ai as DoubleWord;
        for (j, &bj) in b.iter().enumerate() {
            let t = ai * (bj as DoubleWord) + (out[i + j] as DoubleWord) + k;
            out[i + j] = t as Word;
            k = t >> WORD_BIT_SIZE;
        }
        out[i + b.len()] = k as Word;
    }
    Ok(())
}

/// `dst = src * b`; returns the carry limb. `dst` shorter than `src` is `Err(InvalidArgument)`.
fn mpn_mul_1(src: &[Word], b: Word, dst: &mut [Word]) -> Result<Word, Error> {
    if dst.len() < src.len() {
        return Err(Error::InvalidArgument);
    }
    let mut k: DoubleWord = 0;
    let b = b as DoubleWord;
    for (s, d) in src.iter().zip(dst.iter_mut()) {
        let t = (*s as DoubleWord) * b + k;
        *d = t as Word;
        k = t >> WORD_BIT_SIZE;
    }
    Ok(k as Word)
}

/// Left-shift `src` by one bit into `dest`. Returns the bit shifted out of the top. Unequal lengths are
/// `Err(InvalidArgument)`.
fn mpn_lshift1(dest: &mut [Word], src: &[Word]) -> Result<Word, Error> {
    if dest.len() != src.len() {
        return Err(Error::InvalidArgument);
    }
    let mut carry = 0 as Word;
    for (d, &s) in dest.iter_mut().zip(src.iter()) {
        *d = (s << 1) | carry;
        carry = s >> (WORD_BIT_SIZE - 1);
    }
    Ok(carry)
}

fn bit_len(limbs: &[Word]) -> usize {
    for i in (0..limbs.len()).rev() {
        if limbs[i] != 0 {
            return i * WORD_BIT_SIZE + (WORD_BIT_SIZE - limbs[i].leading_zeros() as usize);
        }
    }
    0
}

/// `mpfr_mpn_exp`: truncated `b^e` in `n` limbs, rounded toward zero, with a
/// proven ulp bound. `e > 0`, `2 <= b <= 36`, `n > 0`.
pub(crate) fn mpn_base_pow(base: u8, e: u64, n: usize) -> Result<BasePow, Error> {
    if n == 0 || e == 0 || !(2..=36).contains(&base) {
        return Err(Error::InvalidArgument);
    }

    let mut a = WordBuf::new(n)?;
    let mut c = WordBuf::new(2 * n)?;

    // Normalize the base into the high bit of a limb (`mpn_exp.c` 102–115).
    let mut b_limb: Word = Word::from(base);
    let h0 = b_limb.leading_zeros() as i128;
    b_limb <<= h0;
    let h = -h0;
    a[n - 1] = b_limb;
    // Invariant: A = {a, n} * 2^f.
    let mut f: i128 = h - ((n as i128) - 1) * (WORD_BIT_SIZE as i128);

    let t = 64u32 - e.leading_zeros();
    // `error == t` means no inexact step yet (`mpn_exp.c` 118–121).
    let mut error = t;
    let mut err_s_a2: u32 = 0;
    let mut err_s_ab: u32 = 0;

    for i in (0..t.saturating_sub(1)).rev() {
        // Square the significant limbs of A (`mpn_exp.c` 126–136).
        let bits = match mpn_scan1(&a) {
            Some(s) => n * WORD_BIT_SIZE - s,
            None => return Err(Error::InvalidArgument),
        };
        let n1 = (n * WORD_BIT_SIZE - bits) / WORD_BIT_SIZE;
        c.fill(0);
        mpn_mul(&a[n1..], &a[n1..], &mut c[2 * n1..])?;

        f = f.checked_mul(2).ok_or(Error::ExponentOverflow(Sign::Pos))?;
        f = f
            .checked_add((n * WORD_BIT_SIZE) as i128)
            .ok_or(Error::ExponentOverflow(Sign::Pos))?;

        if c[2 * n - 1] & WORD_SIGNIFICANT_BIT == 0 {
            let _ = mpn_lshift1(&mut a, &c[n..])?;
            a[0] |= c[n - 1] >> (WORD_BIT_SIZE - 1);
            f -= 1;
            if error != t {
                err_s_a2 += 1;
            }
        } else {
            a.copy_from_slice(&c[n..]);
        }

        if error == t && 2 * n1 <= n {
            if let Some(s) = mpn_scan1(&c[2 * n1..]) {
                if s < (n - 2 * n1) * WORD_BIT_SIZE {
                    error = i;
                }
            }
        }

        if ((e >> i) & 1) != 0 {
            c.fill(0);
            c[2 * n - 1] = mpn_mul_1(&a, b_limb, &mut c[n - 1..2 * n - 1])?;
            f = f
                .checked_add(h + i128::from(WORD_BIT_SIZE as u32))
                .ok_or(Error::ExponentOverflow(Sign::Pos))?;
            if c[2 * n - 1] & WORD_SIGNIFICANT_BIT == 0 {
                let _ = mpn_lshift1(&mut a, &c[n..])?;
                a[0] |= c[n - 1] >> (WORD_BIT_SIZE - 1);
                f -= 1;
            } else {
                a.copy_from_slice(&c[n..]);
                if error != t {
                    err_s_ab += 1;
                }
            }
            if error == t && c[n - 1] != 0 {
                error = i;
            }
        }
    }

    let err_ulps = if error == t {
        None
    } else {
        // `error + err_s_ab + err_s_a2/2 + 3`, at most `5t/2 - 1/2` (`mpn_exp.c` 200–205).
        Some(error + err_s_ab + err_s_a2 / 2 + 3)
    };
    Ok(BasePow {
        a,
        exp_r: f,
        err_ulps,
    })
}

fn take_high_limbs(src: &[Word], n: usize) -> &[Word] {
    if src.len() <= n {
        src
    } else {
        &src[src.len() - n..]
    }
}

fn low_limbs_nonzero(src: &[Word], kept_high: usize) -> bool {
    if src.len() <= kept_high {
        return false;
    }
    src[..src.len() - kept_high].iter().any(|&w| w != 0)
}

fn normalize_high(src: &[Word], n: usize, exp: i128) -> Result<(WordBuf, i128), Error> {
    let mut a = WordBuf::new(n)?;
    // Drop high zero limbs so a product that sits in the low half (e.g. `a*1`)
    // is not read as a zero high part.
    let mut last = src.len();
    while last > 0 && src[last - 1] == 0 {
        last -= 1;
    }
    let src = &src[..last];
    let l = bit_len(src);
    if l == 0 {
        return Ok((a, exp));
    }
    let target = n * WORD_BIT_SIZE;
    if l >= target {
        let shift_down = l - target;
        let word_s = shift_down / WORD_BIT_SIZE;
        let bit_s = shift_down % WORD_BIT_SIZE;
        // Keep the high `n` limbs after a right shift of `shift_down` bits:
        // the low `word_s` limbs are discarded, not the high ones.
        // (`src.len() - n - word_s` dropped the high limb of an 11-limb
        // product into 8 limbs, e.g. 1.0 * 3^123 at 192-bit x.)
        let start = word_s;
        if bit_s == 0 {
            a.copy_from_slice(&src[start..start + n]);
        } else {
            let mut carry = if start + n < src.len() { src[start + n] } else { 0 };
            for i in (0..n).rev() {
                let w = src[start + i];
                a[i] = (w >> bit_s) | (carry << (WORD_BIT_SIZE - bit_s));
                carry = w;
            }
        }
        Ok((a, exp + shift_down as i128))
    } else {
        let shift_up = target - l;
        let word_s = shift_up / WORD_BIT_SIZE;
        let bit_s = shift_up % WORD_BIT_SIZE;
        // Room for the shifted integer only (no spare high zero limb that would
        // later be copied in as if it were significant).
        let mut tmp = WordBuf::new(src.len() + word_s + if bit_s != 0 { 1 } else { 0 })?;
        if bit_s == 0 {
            tmp[word_s..word_s + src.len()].copy_from_slice(src);
        } else {
            let mut carry = 0 as Word;
            for i in 0..src.len() {
                let w = src[i];
                tmp[i + word_s] = (w << bit_s) | carry;
                carry = w >> (WORD_BIT_SIZE - bit_s);
            }
            if carry != 0 {
                tmp[src.len() + word_s] = carry;
            }
        }
        let mut used = tmp.len();
        while used > 0 && tmp[used - 1] == 0 {
            used -= 1;
        }
        if used >= n {
            a.copy_from_slice(&tmp[used - n..used]);
        } else {
            a[n - used..].copy_from_slice(&tmp[..used]);
        }
        Ok((a, exp - shift_up as i128))
    }
}

/// `Y = x * b^k` as a truncated limb value with a proven error, after
/// `mpfr_get_str` (get_str.c 2735–2860): `mpfr_mpn_exp` plus mul or div.
pub(crate) fn scale_x_base_pow(
    x_limbs: &[Word],
    x_exp: i64,
    x_prec_bits: usize,
    base: u8,
    k: i64,
    n: usize,
) -> Result<ScaledY, Error> {
    if n == 0 {
        return Err(Error::InvalidArgument);
    }
    let exp_x_unit = i128::from(x_exp) - i128::from(x_prec_bits as i64);

    if k == 0 {
        let (a, exp_a) = normalize_high(x_limbs, n, exp_x_unit)?;
        let err = if low_limbs_nonzero(x_limbs, n) { Some(1u32) } else { None };
        return Ok(ScaledY {
            a,
            exp_a,
            err_ulps: err,
        });
    }

    let exp_abs = k.unsigned_abs();
    let pow = mpn_base_pow(base, exp_abs, n)?;

    if k > 0 {
        // get_str.c 2735–2776: multiply x by b^exp.
        let x1 = take_high_limbs(x_limbs, n);
        let mut exact = pow.err_ulps.is_none() && !low_limbs_nonzero(x_limbs, n);
        let mut prod = WordBuf::new(pow.a.len() + x1.len())?;
        mpn_mul(&pow.a, x1, &mut prod)?;
        if mpn_scan1(&prod).is_some_and(|s| s < x1.len() * WORD_BIT_SIZE) {
            exact = false;
        }
        // prod is (a * x1_int); x = x_int * 2^{exp_x_unit}, x1_int = high limbs
        // of x_int, so if x was truncated, x1_int = high n limbs and the
        // missing low factor is already in `exact`. Value is
        // prod * 2^{exp_r + exp_x_unit} when x1 is the full x_int, and
        // prod * 2^{exp_r + exp_x_unit + (nx-n)*64} when x1 is the high n limbs.
        let low_shift = if x_limbs.len() > n {
            ((x_limbs.len() - n) * WORD_BIT_SIZE) as i128
        } else {
            0
        };
        let (a, exp_a) = normalize_high(&prod, n, pow.exp_r + exp_x_unit + low_shift)?;
        let err = if exact {
            None
        } else {
            // get_str.c 2767–2768: lose one more bit, or two when err <= 0.
            let e0 = pow.err_ulps.unwrap_or(0);
            Some(if e0 == 0 { 2 } else { e0 + 1 })
        };
        Ok(ScaledY {
            a,
            exp_a,
            err_ulps: err,
        })
    } else {
        // get_str.c 2778–2856: divide x by b^exp.
        let two_n = 2 * n;
        let mut x1 = WordBuf::new(two_n)?;
        let mut exact = pow.err_ulps.is_none();
        if x_limbs.len() >= two_n {
            x1.copy_from_slice(&x_limbs[x_limbs.len() - two_n..]);
            if exact && low_limbs_nonzero(x_limbs, two_n) {
                exact = false;
            }
        } else {
            let off = two_n - x_limbs.len();
            x1[off..].copy_from_slice(x_limbs);
        }
        let (q, r) = Mantissa::div_unbalanced(&x1, &pow.a)?;
        if exact && r.iter().any(|&w| w != 0) {
            exact = false;
        }
        // quot * 2^{exp_x_unit - exp_r - 2n*64}  (get_str.c 2833).
        // value = quot * 2^{exp_x_unit - exp_r + (nx - 2n)*64}
        // (MPFR: EXP(x) - exp_r - 2n*64, and EXP(x) = exp_x_unit + nx*64).
        let exp0 = exp_x_unit - pow.exp_r
            + (x_limbs.len() as i128 - two_n as i128) * i128::from(WORD_BIT_SIZE as u32);
        let (a, exp_a) = normalize_high(&q, n, exp0)?;
        let err = if exact {
            None
        } else {
            let e0 = pow.err_ulps.unwrap_or(0);
            Some(if e0 == 0 { 2 } else { e0 + 2 })
        };
        Ok(ScaledY {
            a,
            exp_a,
            err_ulps: err,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_must::Must;

    fn value_as_f_bits(p: &BasePow) -> (u128, i128) {
        // For small n, pack the high limbs into u128 and keep exp_r.
        let n = p.a.len();
        let mut v: u128 = 0;
        let take = n.min(2);
        for &w in p.a[n - take..].iter() {
            v = (v << WORD_BIT_SIZE) | u128::from(w);
        }
        let extra = (n - take) * WORD_BIT_SIZE;
        (v, p.exp_r + extra as i128)
    }

    #[test]
    fn mpn_base_pow_small_exact_powers() -> Result<(), String> {
        // 10^2 = 100, 3^5 = 243, 36^2 = 1296; one limb is enough and exact.
        for (b, e, expect) in
            [(10u8, 2u64, 100u128), (3, 5, 243), (36, 2, 1296), (10, 6, 1_000_000)]
        {
            if e == 0 {
                continue;
            }
            let p = mpn_base_pow(b, e, 1).must("pow")?;
            assert!(p.err_ulps.is_none(), "expected exact {b}^{e}: {p:?} err");
            let (v, exp) = value_as_f_bits(&p);
            // v * 2^exp == expect
            if exp >= 0 {
                assert_eq!(v << exp, expect, "{b}^{e}");
            } else {
                assert_eq!(v >> (-exp), expect, "{b}^{e}");
                assert_eq!(v & ((1u128 << (-exp)) - 1), 0, "{b}^{e} fractional");
            }
        }
        Ok(())
    }

    #[test]
    fn scale_x_base_pow_mul_and_div_small() -> Result<(), String> {
        // 1 * 3^5 = 243, 243 * 3^{-5} = 1, 244 * 3^{-5} slightly above 1.
        let one = [1u64];
        let p = scale_x_base_pow(&one, 64, 64, 3, 5, 2).must("3^5")?;
        // 243 in 2 limbs
        let mut n243 = [0u64; 2];
        n243[0] = 243;
        let q = scale_x_base_pow(&n243, 64 * 2, 64 * 2, 3, -5, 2).must("243/3^5")?;
        let mut n244 = [0u64; 2];
        n244[0] = 244;
        let r = scale_x_base_pow(&n244, 64 * 2, 64 * 2, 3, -5, 2).must("244/3^5")?;
        assert!(p.err_ulps.is_none(), "1*3^5 exact, got {p:?}");
        assert!(q.err_ulps.is_none(), "243/3^5 exact, got {q:?}");
        assert!(r.err_ulps.is_some(), "244/3^5 inexact, got {r:?}");
        // 243 = a * 2^{exp_a}
        let a243 = (u128::from(p.a[1]) << 64) | u128::from(p.a[0]);
        let sh = -p.exp_a;
        assert!((0..128).contains(&sh));
        assert_eq!(a243 >> sh, 243, "1*3^5 value a={a243} exp_a={}", p.exp_a);
        let a1 = (u128::from(q.a[1]) << 64) | u128::from(q.a[0]);
        let sh1 = -q.exp_a;
        assert!((0..128).contains(&sh1));
        assert_eq!(a1 >> sh1, 1, "243/3^5 value a={a1} exp_a={}", q.exp_a);
        Ok(())
    }

    #[test]
    fn scale_x_base_pow_one_times_3_123_at_192_bits() -> Result<(), String> {
        let p = mpn_base_pow(3, 123, 8).must("3^123")?;
        assert!(
            p.err_ulps.is_none(),
            "3^123 should be exact in 8 limbs: {p:?}"
        );
        // 1.0 at 192 bits (3 limbs, MSB of high limb set, exponent 1) times 3^123
        // must equal 3^123. A previous normalize_high used `src.len()-n-word_s`
        // and dropped the high limb of the 11-limb product.
        let mut one192 = [0u64; 3];
        one192[2] = 1u64 << 63;
        let s = scale_x_base_pow(&one192, 1, 192, 3, 123, 8).must("1*3^123")?;
        assert!(s.err_ulps.is_none(), "1*3^123 exact, got {s:?}");
        assert_eq!(&s.a[..], &p.a[..], "limbs of 1*3^123 vs 3^123");
        assert_eq!(s.exp_a, p.exp_r, "exp of 1*3^123 vs 3^123");
        Ok(())
    }
}
