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
    pub(super) fn mul_unbalanced(m1: &[Word], m2: &[Word], m3: &mut [Word]) -> Result<(), Error> {
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
                        dst.copy_from(&src);
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
    use rand::random;

    #[cfg(not(feature = "std"))]
    use alloc::vec::Vec;

    #[test]
    fn test_mul_unbalanced() {
        let sz1 = random::<usize>() % 10 + 1;
        let sz2 = random::<usize>() % 10 * sz1 + random::<usize>() % sz1 + sz1;
        let f = random_slice(1, sz1);
        let mut ret1 = WordBuf::new(sz1 + sz2).unwrap();
        let mut ret2 = WordBuf::new(sz1 + sz2).unwrap();
        for _ in 0..1000 {
            let v = random_slice(sz1, sz2);
            Mantissa::mul_unbalanced(&f, &v, &mut ret1).unwrap();
            Mantissa::mul_slices(&f, &v, &mut ret2).unwrap();
            assert!(ret1[..] == ret2[..]);
        }
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
}
