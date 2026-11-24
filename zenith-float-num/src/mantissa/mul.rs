// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! Multiplication algos.

use crate::common::buf::WordBuf;
use crate::common::int::SliceWithSign;
use crate::defs::DoubleWord;
use crate::defs::Error;
use crate::defs::Word;
use crate::defs::WORD_BIT_SIZE;
use crate::mantissa::Mantissa;

impl Mantissa {
    pub(super) fn mul_basic(m1: &[Word], m2: &[Word], m3: &mut [Word]) {
        m3.fill(0);

        for (i, d1mi) in m1.iter().enumerate() {
            let d1mi = *d1mi as DoubleWord;
            if d1mi == 0 {
                continue;
            }

            let mut k = 0;
            for (m2j, m3ij) in m2.iter().zip(m3[i..].iter_mut()) {
                let m = d1mi * (*m2j as DoubleWord) + *m3ij as DoubleWord + k;
                *m3ij = m as Word;
                k = m >> (WORD_BIT_SIZE);
            }

            m3[i + m2.len()] += k as Word;
        }
    }

    fn mul_slices(m1: &[Word], m2: &[Word], m3: &mut [Word]) -> Result<(), Error> {
        debug_assert!(m1.len() <= m2.len());

        if m1.len() <= 32 || m2.len() <= 32 {
            Self::mul_basic(m1, m2, m3);
        } else if m1.len() <= 220 || m2.len() <= 220 {
            Self::toom2(m1, m2, m3)?;
        } else if m1.len() <= 5400 && m2.len() <= 5400 {
            Self::toom3(m1, m2, m3)?;
        } else {
            Mantissa::fft_mul(m1, m2, m3)?;
        }
        Ok(())
    }

    // general case multiplication
    pub(crate) fn mul_unbalanced(m1: &[Word], m2: &[Word], m3: &mut [Word]) -> Result<(), Error> {
        // Zero low words (a short value held at a wide precision) contribute nothing: multiply
        // only the words above them and zero the same number of low product words. Without this a
        // 1-word value at 65568 bits cost a dense 1025-word Toom-3 product (r2nd/zeta root cause).
        let z1 = m1.iter().position(|w| w != &0).unwrap_or(m1.len());
        let z2 = m2.iter().position(|w| w != &0).unwrap_or(m2.len());
        let z = z1 + z2;
        if z > 0 && m3.len() == m1.len() + m2.len() {
            let (lo, hi) = m3.split_at_mut(z);
            lo.fill(0);
            if z1 == m1.len() || z2 == m2.len() {
                hi.fill(0);
                return Ok(());
            }
            return Self::mul_unbalanced(&m1[z1..], &m2[z2..], hi);
        }
        let (sm, lg) = if m1.len() < m2.len() { (m1, m2) } else { (m2, m1) };

        if lg.len() / 2 >= sm.len() && sm.len() > 70 {
            // balancing

            let mut buf = WordBuf::new(2 * sm.len())?;
            let mut even = true;
            let mut lb = 0;
            let mut ub = 0;

            for _ in 0..2 {
                while lb < lg.len() {
                    ub = if lb + sm.len() <= lg.len() { lb + sm.len() } else { lg.len() };

                    Self::mul_slices(&lg[lb..ub], sm, &mut buf)?;

                    let src = SliceWithSign::new(&buf[..ub - lb + sm.len()], 1);
                    let mut dst = SliceWithSign::new_mut(&mut m3[lb..], 1);

                    if even {
                        dst.copy_from(&src)?;
                    } else {
                        dst.add_assign(&src);
                    }

                    lb += sm.len() * 2;
                }

                if even {
                    if ub + sm.len() < m3.len() {
                        m3[ub + sm.len()..].fill(0);
                    }

                    even = false;
                    lb = sm.len();
                }
            }

            Ok(())
        } else {
            Self::mul_slices(sm, lg, m3)
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::common::test_rng::random;

    #[cfg(not(feature = "std"))]
    use alloc::vec::Vec;

    #[test]
    fn test_mul_unbalanced() -> Result<(), crate::Error> {
        let sz1 = random::<usize>() % 10 + 1;
        let sz2 = random::<usize>() % 10 * sz1 + random::<usize>() % sz1 + sz1;
        let f = random_slice(1, sz1);
        let mut ret1 = WordBuf::new(sz1 + sz2)?;
        let mut ret2 = WordBuf::new(sz1 + sz2)?;
        for _ in 0..1000 {
            let v = random_slice(sz1, sz2);
            Mantissa::mul_unbalanced(&f, &v, &mut ret1)?;
            Mantissa::mul_slices(&f, &v, &mut ret2)?;
            assert!(ret1[..] == ret2[..]);
        }
        Ok(())
    }

    fn random_slice(min_len: usize, max_len: usize) -> Vec<Word> {
        let mut s1 = Vec::new();
        let l = if max_len > min_len {
            random::<usize>() % (max_len - min_len) + min_len
        } else {
            min_len
        };
        for _ in 0..l {
            s1.push(random());
        }
        s1
    }

    fn algo_for_limbs(limbs: usize) -> &'static str {
        if limbs <= 32 {
            "schoolbook"
        } else if limbs <= 220 {
            "toom2"
        } else if limbs <= 5400 {
            "toom3"
        } else {
            "fft"
        }
    }

    fn fill_odd(buf: &mut [Word], seed: u64) {
        let mut x = seed;
        for w in buf.iter_mut() {
            x = x.wrapping_mul(0x9E37_79B9_7F4A_7C15).wrapping_add(1);
            *w = x;
        }
        if let Some(top) = buf.last_mut() {
            *top |= 1u64 << 63;
        }
    }

    fn add_limb_slice(a: &[Word], b: &[Word], out: &mut [Word]) {
        let mut carry: crate::defs::DoubleWord = 0;
        for i in 0..a.len() {
            let s = a[i] as crate::defs::DoubleWord + b[i] as crate::defs::DoubleWord + carry;
            out[i] = s as Word;
            carry = s >> crate::defs::WORD_BIT_SIZE;
        }
        if out.len() > a.len() {
            out[a.len()] = carry as Word;
        }
    }

    fn time_ns(
        iters: u32,
        warmup: u32,
        mut body: impl FnMut() -> Result<(), String>,
    ) -> Result<u128, String> {
        for _ in 0..warmup {
            body()?;
        }
        let t0 = std::time::Instant::now();
        for _ in 0..iters {
            body()?;
        }
        Ok(t0.elapsed().as_nanos() / u128::from(iters.max(1)))
    }

    /// Microbench of one full-width multiply, a squaring, a mul-by-one-word,
    /// and an add, at the EM dest width and nearby / algorithm-boundary limb
    /// counts. Paired with `$L/final/timeouts-11/zeta/bench-arb-kernel.py`.
    #[cfg(feature = "std")]
    #[test]
    #[ignore]
    fn time_mul_kernel_vs_flint_widths() -> Result<(), String> {
        use crate::{ExactNum, RoundingMode};
        use std::hint::black_box;

        let rm = RoundingMode::ToEven;
        let word: u64 = 1234567890123456789;

        // Bit widths matching the python-flint arb script, plus the algorithm
        // thresholds (32 / 220 / 5400 limbs).
        let cases: &[(usize, u32, u32, u32)] = &[
            (64, 20000, 50000, 20000),
            (2048, 4000, 12000, 8000), // 32 limbs, schoolbook top
            (2112, 4000, 12000, 8000), // 33 limbs, toom2 bottom
            (4096, 2000, 8000, 4000),
            (14080, 400, 2000, 800), // 220 limbs, toom2 top
            (14144, 400, 2000, 800), // 221 limbs, toom3 bottom
            (16384, 400, 2000, 800),
            (65536, 80, 400, 200),
            (262144, 10, 80, 40),
            (345600, 2, 20, 8), // 5400 limbs, toom3 top
            (345664, 2, 20, 8), // 5401 limbs, fft bottom
        ];

        for &(prec, n_mul, n_add, n_ui) in cases {
            let limbs = prec.div_ceil(crate::defs::WORD_BIT_SIZE);
            eprintln!(
                "==== ZF prec={prec} bits limbs={limbs} algo={} ====",
                algo_for_limbs(limbs)
            );

            let mut ma = WordBuf::new(limbs)?;
            let mut mb = WordBuf::new(limbs)?;
            let mut mc = WordBuf::new(2 * limbs)?;
            let mut md = WordBuf::new(limbs + 1)?;
            let mut me = WordBuf::new(limbs + 1)?;
            fill_odd(&mut ma, 0x2222);
            fill_odd(&mut mb, 0x3333);

            let ns = time_ns(n_mul, n_mul.min(8), || {
                Mantissa::mul_unbalanced(&ma, &mb, &mut mc)?;
                black_box(&mc[0]);
                Ok(())
            })?;
            eprintln!(
                "ZF  mant-mul  p={prec}: n={n_mul}  {ns:.0} ns/op  {} us/op",
                ns / 1_000
            );

            let ns = time_ns(n_mul, n_mul.min(8), || {
                Mantissa::mul_unbalanced(&ma, &ma, &mut mc)?;
                black_box(&mc[0]);
                Ok(())
            })?;
            eprintln!(
                "ZF  mant-sqr  p={prec}: n={n_mul}  {ns:.0} ns/op  {} us/op",
                ns / 1_000
            );

            let ns = time_ns(n_ui, n_ui.min(8), || {
                Mantissa::mul_by_word(&ma, word as crate::defs::DoubleWord, &mut md);
                black_box(&md[0]);
                Ok(())
            })?;
            eprintln!(
                "ZF  mant-mul_ui p={prec}: n={n_ui}  {ns:.0} ns/op  {} us/op",
                ns / 1_000
            );

            let ns = time_ns(n_add, n_add.min(8), || {
                add_limb_slice(&ma, &mb, &mut me);
                black_box(&me[0]);
                Ok(())
            })?;
            eprintln!(
                "ZF  mant-add  p={prec}: n={n_add}  {ns:.0} ns/op  {} us/op",
                ns / 1_000
            );

            let a = ExactNum::from_u8(2, prec).sqrt(prec, rm);
            let b = ExactNum::from_u8(3, prec).sqrt(prec, rm);
            assert!(!a.is_nan() && !b.is_nan() && !a.is_zero() && !b.is_zero());
            let w = ExactNum::from_u64(word, crate::defs::WORD_BIT_SIZE);

            let ns = time_ns(n_mul, n_mul.min(8), || {
                black_box(ExactNum::mul_ref(&a, &b, prec, rm));
                Ok(())
            })?;
            eprintln!(
                "ZF  ExactNum-mul  p={prec}: n={n_mul}  {ns:.0} ns/op  {} us/op",
                ns / 1_000
            );

            let ns = time_ns(n_mul, n_mul.min(8), || {
                black_box(ExactNum::mul_ref(&a, &a, prec, rm));
                Ok(())
            })?;
            eprintln!(
                "ZF  ExactNum-sqr  p={prec}: n={n_mul}  {ns:.0} ns/op  {} us/op",
                ns / 1_000
            );

            let ns = time_ns(n_ui, n_ui.min(8), || {
                black_box(ExactNum::mul_ref(&a, &w, prec, rm));
                Ok(())
            })?;
            eprintln!(
                "ZF  ExactNum-mul_ui p={prec}: n={n_ui}  {ns:.0} ns/op  {} us/op",
                ns / 1_000
            );

            let ns = time_ns(n_add, n_add.min(8), || {
                black_box(ExactNum::add_ref(&a, &b, prec, rm));
                Ok(())
            })?;
            eprintln!(
                "ZF  ExactNum-add  p={prec}: n={n_add}  {ns:.0} ns/op  {} us/op",
                ns / 1_000
            );
        }
        Ok(())
    }

    /// A short value held at a wide precision has zero low words (exact small integers at the
    /// working precision of a binary-splitting tree, `ops::zeta_ui`); its product must be the same
    /// as the schoolbook product of the full slices, for every algorithm range and padding.
    #[test]
    fn zero_low_words_give_the_same_product() -> Result<(), String> {
        for &(l1, l2) in
            &[(1, 1), (3, 40), (33, 33), (40, 300), (221, 230), (1025, 1025), (700, 1500)]
        {
            for &(z1, z2) in &[(0, 1), (1, 0), (l1 - 1, l2 - 1), (l1 / 2, l2 / 3), (l1 - 1, 0)] {
                let mut a = WordBuf::new(l1)?;
                let mut b = WordBuf::new(l2)?;
                fill_odd(&mut a, 0x5151 + l1 as u64);
                fill_odd(&mut b, 0x7373 + l2 as u64);
                a[..z1].fill(0);
                b[..z2].fill(0);
                let mut got = WordBuf::new(l1 + l2)?;
                let mut want = WordBuf::new(l1 + l2)?;
                got.fill(Word::MAX);
                Mantissa::mul_unbalanced(&a, &b, &mut got)?;
                if l1 <= l2 {
                    Mantissa::mul_basic(&a, &b, &mut want);
                } else {
                    Mantissa::mul_basic(&b, &a, &mut want);
                }
                assert!(got[..] == want[..], "l1={l1} l2={l2} z1={z1} z2={z2}");
            }
        }
        Ok(())
    }

    /// Root cause of the debug-build `zf_riemann_zeta_matches_flint_golds` time (r2nd/zeta): the
    /// zero low words of short operands were multiplied as if dense, so a 1-word x 1-word product
    /// at 65568 bits cost a dense 1025 x 1025-word Toom-3. It must cost a small fraction of one.
    #[cfg(feature = "std")]
    #[test]
    fn padded_short_operands_skip_their_zero_words() -> Result<(), String> {
        let limbs = 1025;
        let mut da = WordBuf::new(limbs)?;
        let mut db = WordBuf::new(limbs)?;
        fill_odd(&mut da, 0x2222);
        fill_odd(&mut db, 0x3333);
        let mut pa = WordBuf::new(limbs)?;
        let mut pb = WordBuf::new(limbs)?;
        pa.fill(0);
        pb.fill(0);
        pa[limbs - 1] = 0x8000_0000_075b_cd15;
        pb[limbs - 1] = 0x8000_0000_3ade_68b1;
        let mut mc = WordBuf::new(2 * limbs)?;
        let dense = time_ns(12, 2, || {
            Mantissa::mul_unbalanced(&da, &db, &mut mc).map_err(|e| format!("{e:?}"))
        })?;
        let padded = time_ns(12, 2, || {
            Mantissa::mul_unbalanced(&pa, &pb, &mut mc).map_err(|e| format!("{e:?}"))
        })?;
        assert!(
            crate::test_time::limit_holds(padded.saturating_mul(20) < dense),
            "padded 1-word operands cost {padded} ns vs {dense} ns for a dense 1025-word product"
        );
        Ok(())
    }
}
