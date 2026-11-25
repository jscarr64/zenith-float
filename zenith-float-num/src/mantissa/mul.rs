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
        if !(m1.len() <= m2.len()) {
            return Err(Error::InvalidArgument);
        }
        #[cfg(test)]
        SLICE_WORDS.with(|c| c.set(c.get() + m1.len() * m2.len()));

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
std::thread_local! {
    /// Sum of `len(m1) * len(m2)` over every `mul_slices` call made by this thread: the size of the word products
    /// handed to the multiplication algorithms (unit tests count it instead of timing them).
    static SLICE_WORDS: core::cell::Cell<usize> = const { core::cell::Cell::new(0) };
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::common::test_rng::random;

    #[cfg(not(feature = "std"))]
    use alloc::vec::Vec;

    /// Jeff 2026-11-24 4:40 PM CT item 5: mul_slices needs the shorter operand first; it refuses otherwise, in every
    /// build (was a debug-only assert).
    #[test]
    fn mul_slices_refuses_longer_first_operand() {
        let mut out = [0 as Word; 3];
        assert!(matches!(
            Mantissa::mul_slices(&[1, 2], &[3], &mut out),
            Err(Error::InvalidArgument)
        ));
        assert!(Mantissa::mul_slices(&[3], &[1, 2], &mut out).is_ok());
        assert_eq!(out, [3, 6, 0]);
    }

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
    /// at 65568 bits cost a dense 1025 x 1025-word Toom-3. Counted, not timed (r6: the wall-clock
    /// version was flaky): the padded product hands 1 x 1 words to the algorithms, the dense one at
    /// least 1025 x 1025, and the padded product is right.
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
        let words = |a: &[Word], b: &[Word]| -> Result<(usize, WordBuf), String> {
            let mut mc = WordBuf::new(2 * limbs)?;
            SLICE_WORDS.with(|c| c.set(0));
            Mantissa::mul_unbalanced(a, b, &mut mc).map_err(|e| format!("{e:?}"))?;
            Ok((SLICE_WORDS.with(core::cell::Cell::get), mc))
        };
        let (dense, _) = words(&da, &db)?;
        let (padded, got) = words(&pa, &pb)?;
        // Toom-3 recurses through mul_slices, so the dense count includes its sub-products: at least 1025 x 1025.
        assert!(
            dense >= limbs * limbs,
            "dense 1025-word product counted {dense}"
        );
        assert_eq!(padded, 1, "padded 1-word operands");
        let mut want = WordBuf::new(2 * limbs)?;
        Mantissa::mul_basic(&pa, &pb, &mut want);
        assert!(got[..] == want[..], "padded product");
        Ok(())
    }
}
