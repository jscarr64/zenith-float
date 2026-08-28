//! Buffer for holding mantissa digits.
//!
//! Lengths of at most [`INLINE_WORDS`] words are stored inline (no heap allocation).
//! Larger buffers promote to `Vec`.

use crate::defs::Error;
use crate::defs::Word;
use crate::defs::WORD_BIT_SIZE;
use core::hash::{Hash, Hasher};
use core::ops::Deref;
use core::ops::DerefMut;
use core::ops::Index;
use core::ops::IndexMut;
use core::slice::SliceIndex;

use crate::common::util::shift_slice_left;
use crate::common::util::shift_slice_right;

use alloc::vec::Vec;

/// Number of words kept on the stack before allocating.
pub const INLINE_WORDS: usize = 2;

/// Buffer for holding mantissa digits.
#[derive(Debug)]
pub struct WordBuf {
    inner: Storage,
}

#[derive(Debug)]
enum Storage {
    Inline { data: [Word; INLINE_WORDS], len: u8 },
    Heap(Vec<Word>),
}

impl WordBuf {
    #[inline]
    pub fn new(sz: usize) -> Result<Self, Error> {
        if sz <= INLINE_WORDS {
            Ok(WordBuf {
                inner: Storage::Inline {
                    data: [0; INLINE_WORDS],
                    len: sz as u8,
                },
            })
        } else {
            let mut inner = Vec::new();
            inner.try_reserve_exact(sz)?;
            inner.resize(sz, 0);
            Ok(WordBuf {
                inner: Storage::Heap(inner),
            })
        }
    }

    #[inline]
    fn as_slice(&self) -> &[Word] {
        match &self.inner {
            Storage::Inline { data, len } => &data[..*len as usize],
            Storage::Heap(v) => v.as_slice(),
        }
    }

    #[inline]
    fn as_mut_slice(&mut self) -> &mut [Word] {
        match &mut self.inner {
            Storage::Inline { data, len } => &mut data[..*len as usize],
            Storage::Heap(v) => v.as_mut_slice(),
        }
    }

    fn resize_len(&mut self, new_len: usize) -> Result<(), Error> {
        if new_len <= self.len() {
            self.truncate_words(new_len);
            return Ok(());
        }
        match &mut self.inner {
            Storage::Inline { len, .. } if new_len <= INLINE_WORDS => {
                *len = new_len as u8;
                Ok(())
            }
            Storage::Heap(v) => {
                v.try_reserve(new_len - v.len())?;
                v.resize(new_len, 0);
                Ok(())
            }
            Storage::Inline { data, len } => {
                let mut v = Vec::new();
                v.try_reserve_exact(new_len)?;
                v.extend_from_slice(&data[..*len as usize]);
                v.resize(new_len, 0);
                self.inner = Storage::Heap(v);
                Ok(())
            }
        }
    }

    #[inline]
    pub fn fill(&mut self, d: Word) {
        self.as_mut_slice().fill(d);
    }

    #[inline]
    pub fn len(&self) -> usize {
        match &self.inner {
            Storage::Inline { len, .. } => *len as usize,
            Storage::Heap(v) => v.len(),
        }
    }

    /// True when the buffer is stored without a heap allocation.
    #[inline]
    pub fn is_inline(&self) -> bool {
        matches!(self.inner, Storage::Inline { .. })
    }

    /// Decrease length of the buffer to l bits. Data is shifted.
    pub fn trunc_to(&mut self, l: usize) {
        let n = (l + WORD_BIT_SIZE - 1) / WORD_BIT_SIZE;
        let sz = self.len();
        if n >= sz {
            return;
        }
        shift_slice_right(self.as_mut_slice(), (sz - n) * WORD_BIT_SIZE);
        self.truncate_words(n);
    }

    /// Decrease length of the buffer to l bits. Data is not moved.
    pub fn trunc_to_2(&mut self, l: usize) {
        let n = (l + WORD_BIT_SIZE - 1) / WORD_BIT_SIZE;
        self.truncate_words(n);
    }

    fn truncate_words(&mut self, n: usize) {
        match &mut self.inner {
            Storage::Inline { len, .. } => {
                *len = (*len).min(n as u8);
            }
            Storage::Heap(v) => v.truncate(n),
        }
    }

    /// Try to extend the size to fit the precision p. Data is shifted to the left.
    pub fn try_extend(&mut self, p: usize) -> Result<(), Error> {
        let n = (p + WORD_BIT_SIZE - 1) / WORD_BIT_SIZE;
        let l = self.len();
        if n > l {
            self.resize_len(n)?;
            shift_slice_left(self.as_mut_slice(), (n - l) * WORD_BIT_SIZE);
        }
        Ok(())
    }

    /// Try to extend the size to fit the precision p. Fill new elements with 0. Data is not moved.
    pub fn try_extend_2(&mut self, p: usize) -> Result<(), Error> {
        let n = (p + WORD_BIT_SIZE - 1) / WORD_BIT_SIZE;
        if n > self.len() {
            self.resize_len(n)?;
        }
        Ok(())
    }

    /// Try to extend the size to fit the precision p. Data is shifted to the left by d bits.
    pub fn try_extend_3(&mut self, p: usize, d: usize) -> Result<(), Error> {
        let n = (p + WORD_BIT_SIZE - 1) / WORD_BIT_SIZE;
        let l = self.len();
        if n > l {
            self.resize_len(n)?;
        }
        shift_slice_left(self.as_mut_slice(), d);
        Ok(())
    }

    // Remove trailing words containing zeroes.
    pub fn trunc_trailing_zeroes(&mut self) {
        let mut n = 0;
        for v in self.as_slice().iter() {
            if *v == 0 {
                n += 1;
            } else {
                break;
            }
        }
        if n > 0 {
            let sz = self.len();
            shift_slice_right(self.as_mut_slice(), n * WORD_BIT_SIZE);
            self.truncate_words(sz - n);
        }
    }

    // Remove leading words containing zeroes.
    pub fn trunc_leading_zeroes(&mut self) {
        let mut n = 0;
        for v in self.as_slice().iter().rev() {
            if *v == 0 {
                n += 1;
            } else {
                break;
            }
        }
        if n > 0 {
            let sz = self.len();
            self.truncate_words(sz - n);
        }
    }
}

impl Hash for WordBuf {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.as_slice().hash(state);
    }
}

impl<I: SliceIndex<[Word]>> IndexMut<I> for WordBuf {
    #[inline]
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        self.as_mut_slice().index_mut(index)
    }
}

impl<I: SliceIndex<[Word]>> Index<I> for WordBuf {
    type Output = I::Output;

    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        self.as_slice().index(index)
    }
}

impl Deref for WordBuf {
    type Target = [Word];

    #[inline]
    fn deref(&self) -> &[Word] {
        self.as_slice()
    }
}

impl DerefMut for WordBuf {
    #[inline]
    fn deref_mut(&mut self) -> &mut [Word] {
        self.as_mut_slice()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_inline_small_and_promote() {
        let mut b = WordBuf::new(1).unwrap();
        assert!(b.is_inline());
        assert_eq!(b.len(), 1);
        b[0] = 7;
        b.try_extend_2(WORD_BIT_SIZE * 8).unwrap();
        assert!(!b.is_inline());
        assert_eq!(b.len(), 8);
        assert_eq!(b[0], 7);
    }

    #[test]
    fn huge_reserve_returns_memory_error() {
        let words = (isize::MAX as usize) / core::mem::size_of::<Word>() + 1;
        assert!(matches!(WordBuf::new(words), Err(Error::MemoryAllocation)));
    }
}
