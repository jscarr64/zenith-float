// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! Lightweigh integer.

use crate::common::util::add_carry;
use crate::common::util::shift_slice_right;
use crate::common::util::{shift_slice_left, sub_borrow};
use crate::defs::{DoubleWord, Error, SignedWord, Word, WORD_BASE, WORD_MAX};
use core::ops::Deref;
use core::ops::DerefMut;
use itertools::izip;

// Internal repr of SliceWithSign.
#[derive(Debug)]
enum SliceWithSignType<'a> {
    Mut(&'a mut [Word]),
    Immut(&'a [Word]),
}

// Slice of words with sign is a lightweight integer number representation.
#[derive(Debug)]
pub struct SliceWithSign<'a> {
    m: SliceWithSignType<'a>,
    sign: i8,
}

impl<'a> SliceWithSign<'a> {
    pub fn new_mut(m: &'a mut [Word], sign: i8) -> Self {
        SliceWithSign {
            m: SliceWithSignType::Mut(m),
            sign,
        }
    }

    pub fn new(m: &'a [Word], sign: i8) -> Self {
        SliceWithSign {
            m: SliceWithSignType::Immut(m),
            sign,
        }
    }

    #[inline]
    pub fn add(&self, s2: &SliceWithSign<'_>, dst: &mut SliceWithSign<'_>) -> Result<(), Error> {
        self.add_sub(s2, dst, 1)
    }

    #[inline]
    pub fn sub(&self, s2: &SliceWithSign<'_>, dst: &mut SliceWithSign<'_>) -> Result<(), Error> {
        self.add_sub(s2, dst, -1)
    }

    #[inline]
    pub fn add_assign(&mut self, s2: &SliceWithSign<'_>) {
        self.add_sub_assign(s2, 1);
    }

    #[inline]
    pub fn sub_assign(&mut self, s2: &SliceWithSign<'_>) {
        self.add_sub_assign(s2, -1);
    }

    fn add_sub(
        &self,
        s2: &SliceWithSign<'_>,
        dst: &mut SliceWithSign<'_>,
        op: i8,
    ) -> Result<(), Error> {
        if (self.sign != s2.sign && op > 0) || (op < 0 && self.sign == s2.sign) {
            // subtract
            let cmp = Self::abs_cmp(self, s2);

            if cmp > 0 {
                dst.sign = self.sign;
                Self::abs_sub(self, s2, dst);
            } else if cmp < 0 {
                dst.sign = s2.sign * op;
                Self::abs_sub(s2, self, dst);
            } else {
                dst.fill(0);
            };
            Ok(())
        } else {
            dst.sign = self.sign;
            Self::abs_add(self, s2, dst)
        }
    }

    fn add_sub_assign(&mut self, s2: &SliceWithSign<'_>, op: i8) {
        if (self.sign != s2.sign && op > 0) || (op < 0 && self.sign == s2.sign) {
            // subtract
            let cmp = Self::abs_cmp(self, s2);
            if cmp > 0 {
                Self::abs_sub_assign_1(self, s2);
            } else if cmp < 0 {
                Self::abs_sub_assign_2(self, s2);
                self.sign = s2.sign * op;
            } else {
                self.fill(0);
            };
        } else {
            Self::abs_add_assign(self, s2);
        }
    }

    #[cfg(test)]
    pub fn mul_assign<'c>(&mut self, s2: &SliceWithSign<'c>, work_buf: &mut [Word]) {
        work_buf.fill(0);
        for (i, d1mi) in self.deref().iter().enumerate() {
            let d1mi = *d1mi as DoubleWord;
            if d1mi == 0 {
                continue;
            }

            let mut k = 0;
            for (m2j, m3ij) in s2.deref().iter().zip(work_buf[i..].iter_mut()) {
                let m = d1mi * (*m2j as DoubleWord) + *m3ij as DoubleWord + k;

                *m3ij = m as Word;
                k = m >> (crate::WORD_BIT_SIZE);
            }
            work_buf[i + s2.len()] += k as Word;
        }
        self.deref_mut().copy_from_slice(work_buf);
        self.sign *= s2.sign;
    }

    #[inline]
    pub fn shift_left(&mut self, shift: usize) {
        shift_slice_left(self, shift);
    }

    #[inline]
    pub fn shift_right(&mut self, shift: usize) {
        shift_slice_right(self, shift);
    }

    pub fn div_by_word(&mut self, d: Word) {
        debug_assert!(d != 0);

        let d = d as DoubleWord;
        let mut rh = 0;
        let m = self.deref_mut();
        let mut iter = m.iter_mut().rev();
        let mut val;

        if let Some(v) = iter.next() {
            val = v;
        } else {
            return;
        }

        if (*val as DoubleWord) < d {
            rh = *val as DoubleWord;
            *val = 0;
            if let Some(v) = iter.next() {
                val = v;
            } else {
                return;
            }
        }

        loop {
            let qh = rh * WORD_BASE as DoubleWord + *val as DoubleWord;

            rh = qh % d;

            *val = (qh / d) as Word;

            if let Some(v) = iter.next() {
                val = v;
            } else {
                break;
            }
        }
    }

    /// Copies the words of `s2` into the front of `self`.
    /// `Err(InvalidArgument)` if `self` is immutable or shorter than `s2`.
    pub fn copy_from(&mut self, s2: &SliceWithSign<'_>) -> Result<(), Error> {
        let SliceWithSignType::Mut(m) = &mut self.m else {
            return Err(Error::InvalidArgument);
        };
        let src: &[Word] = match &s2.m {
            SliceWithSignType::Mut(s) => s,
            SliceWithSignType::Immut(s) => s,
        };
        let Some(dst) = m.get_mut(..src.len()) else {
            return Err(Error::InvalidArgument);
        };
        dst.copy_from_slice(src);
        Ok(())
    }

    fn abs_add(s1: &[Word], s2: &[Word], dst: &mut [Word]) -> Result<(), Error> {
        let mut c = 0;

        let (iter1, mut iter2) = if s1.len() < s2.len() {
            (s1.iter(), s2.iter())
        } else {
            (s2.iter(), s1.iter())
        };

        let mut iter3 = dst.iter_mut();

        for (a, b, x) in izip!(iter1, iter2.by_ref(), iter3.by_ref()) {
            c = add_carry(*a, *b, c, x);
        }

        for (b, x) in iter2.zip(iter3.by_ref()) {
            c = add_carry(0, *b, c, x);
        }

        if c > 0 {
            let Some(carry_slot) = iter3.next() else {
                return Err(Error::InvalidArgument); // dst must be longer than s1 and s2 to hold the carry
            };
            *carry_slot = c as Word;
        }

        for v in iter3 {
            *v = 0;
        }
        Ok(())
    }

    fn abs_add_assign(s1: &mut [Word], s2: &[Word]) {
        let mut c = 0;
        let mut iter1 = s1.iter_mut();
        let iter2 = s2.iter();

        for (b, a) in izip!(iter2, iter1.by_ref()) {
            c = add_carry(*a, *b, c, a);
        }

        for a in iter1 {
            c = add_carry(*a, 0, c, a);
        }
    }

    // prereq: val of s1 >= val of s2
    fn abs_sub_assign_1(s1: &mut [Word], s2: &[Word]) {
        let mut c = 0;
        let mut iter1 = s1.iter_mut();
        let iter2 = s2.iter();

        for (b, a) in izip!(iter2, iter1.by_ref()) {
            c = sub_borrow(*a, *b, c, a);
        }

        for a in iter1 {
            c = sub_borrow(*a, 0, c, a);
        }

        debug_assert!(c == 0);
    }

    // prereq: val of s2 > val of s1
    fn abs_sub_assign_2(s1: &mut [Word], s2: &[Word]) {
        let mut c = 0;

        for (a, b) in s2.iter().zip(s1.iter_mut()) {
            c = sub_borrow(*a, *b, c, b);
        }

        debug_assert!(c == 0);
    }

    fn abs_sub(s1: &[Word], s2: &[Word], dst: &mut [Word]) {
        let mut c = 0;

        let mut iter1 = s1.iter();
        let iter2 = s2.iter();
        let mut iter3 = dst.iter_mut();

        for (b, a, d) in izip!(iter2, iter1.by_ref(), iter3.by_ref()) {
            c = sub_borrow(*a, *b, c, d);
        }

        if c > 0 {
            for (a, d) in iter1.zip(iter3.by_ref()) {
                c = sub_borrow(*a, 0, c, d);
            }
        } else {
            for (a, d) in iter1.zip(iter3.by_ref()) {
                *d = *a;
            }
        }

        for v in iter3 {
            *v = 0;
        }
    }

    // decrement absolute value by 1
    pub fn decrement_abs(&mut self) -> Result<(), Error> {
        for v in self.iter_mut() {
            if *v == 0 {
                *v = WORD_MAX;
            } else {
                *v -= 1;
                return Ok(());
            }
        }
        Err(Error::InvalidArgument) // |self| was zero: nothing to decrement
    }

    pub fn is_zero(&self) -> bool {
        for v in self.iter() {
            if *v != 0 {
                return false;
            }
        }
        true
    }

    fn abs_cmp(s1: &[Word], s2: &[Word]) -> SignedWord {
        let len = s1.len().min(s2.len());

        for v in &s1[len..] {
            if *v != 0 {
                return 1;
            }
        }

        for v in &s2[len..] {
            if *v != 0 {
                return -1;
            }
        }

        for (a, b) in core::iter::zip(s1[..len].iter().rev(), s2[..len].iter().rev()) {
            let diff = *a as SignedWord - *b as SignedWord;
            if diff != 0 {
                return diff;
            }
        }

        0
    }

    #[inline]
    pub fn set_sign(&mut self, sign: i8) {
        self.sign = sign;
    }

    #[inline]
    pub fn sign(&self) -> i8 {
        self.sign
    }

    pub fn cmp(&self, d2: &Self) -> SignedWord {
        if self.is_zero() && d2.is_zero() {
            return 0;
        }

        if self.sign != d2.sign {
            return self.sign as SignedWord;
        }

        Self::abs_cmp(self, d2) * self.sign as SignedWord
    }
}

impl<'a> Deref for SliceWithSign<'a> {
    type Target = [Word];

    #[inline]
    fn deref(&self) -> &[Word] {
        match &self.m {
            SliceWithSignType::Mut(m) => m,
            SliceWithSignType::Immut(m) => m,
        }
    }
}

impl<'a> DerefMut for SliceWithSign<'a> {
    #[inline]
    fn deref_mut(&mut self) -> &mut [Word] {
        match &mut self.m {
            SliceWithSignType::Mut(m) => m,
            _ => panic!(),
        }
    }
}

#[cfg(test)]
mod failure_path_tests {
    use super::*;

    #[test]
    fn add_without_room_for_carry_is_error_not_panic() {
        let a = SliceWithSign::new(&[WORD_MAX, 1], 1);
        let b = SliceWithSign::new(&[1, WORD_MAX], 1);
        let mut buf = [0 as Word; 2];
        let mut dst = SliceWithSign::new_mut(&mut buf, 1);
        assert_eq!(a.add(&b, &mut dst), Err(Error::InvalidArgument));

        let mut buf = [0 as Word; 3];
        let mut dst = SliceWithSign::new_mut(&mut buf, 1);
        assert_eq!(a.add(&b, &mut dst), Ok(()));
        assert_eq!(*dst, [0, 1, 1]);
    }

    #[test]
    fn copy_from_into_immutable_or_short_is_error_not_panic() {
        let src_buf = [1 as Word, 2, 3];
        let src = SliceWithSign::new(&src_buf, 1);
        let target = [0 as Word; 3];
        let mut immut = SliceWithSign::new(&target, 1);
        assert_eq!(immut.copy_from(&src), Err(Error::InvalidArgument));

        let mut short = [0 as Word; 2];
        let mut dst = SliceWithSign::new_mut(&mut short, 1);
        assert_eq!(dst.copy_from(&src), Err(Error::InvalidArgument));

        let mut ok = [0 as Word; 4];
        let mut dst = SliceWithSign::new_mut(&mut ok, 1);
        assert_eq!(dst.copy_from(&src), Ok(()));
        assert_eq!(*dst, [1, 2, 3, 0]);
    }

    #[test]
    fn decrement_of_zero_is_error_not_panic() {
        let mut zero = [0 as Word; 2];
        let mut z = SliceWithSign::new_mut(&mut zero, 1);
        assert_eq!(z.decrement_abs(), Err(Error::InvalidArgument));

        let mut one_carry = [0 as Word, 1];
        let mut v = SliceWithSign::new_mut(&mut one_carry, 1);
        assert_eq!(v.decrement_abs(), Ok(()));
        assert_eq!(*v, [WORD_MAX, 0]);
    }
}
