// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! Auxiliary functions.

use crate::{
    defs::{DoubleWord, Error, Word, WORD_BIT_SIZE, WORD_MAX, WORD_SIGNIFICANT_BIT},
    RoundingMode,
};

#[cfg(test)]
use crate::{num::ExactNumNumber, Sign, EXPONENT_MIN};

#[cfg(test)]
#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

/// integer logarithm base 2 of a number.
pub fn log2_ceil(mut n: usize) -> usize {
    let mut ret = 0;
    let mut sticky = 0;
    while n > 1 {
        if n & 1 != 0 {
            sticky = 1;
        }
        ret += 1;
        n >>= 1;
    }
    ret + sticky
}

/// integer logarithm base 2 of a number.
pub fn log2_floor(mut n: usize) -> usize {
    let mut ret = 0;
    while n > 1 {
        ret += 1;
        n >>= 1;
    }
    ret
}

/// square root integer approximation.
pub fn sqrt_int(a: u32) -> u32 {
    let a = a as u64;
    let mut x = a;
    for _ in 0..20 {
        if x == 0 {
            break;
        }
        x = (a / x + x) >> 1;
    }
    x as u32
}

// cost of multiplication of two numbers with precision p.
pub fn calc_mul_cost(p: usize) -> usize {
    if p < 70 {
        p * p
    } else {
        // toom-3
        if p < 1625 {
            sqrt_int((p * p * p) as u32) as usize
        } else {
            let q = sqrt_int(p as u32) as usize;
            q * q * q
        }
    }
}

// cost of addition/subtraction of two numbers with precision p.
#[inline]
pub fn calc_add_cost(p: usize) -> usize {
    p
}

/// Fixed part of the extra-precision budget of a correct-rounding retry loop, in words.
///
/// The budget is `max(WORD_BIT_SIZE * MAX_PREC_RETRY, 2 * p)` extra bits; see [`prec_retry_exhausted`].
pub const MAX_PREC_RETRY: usize = 256;

#[cfg(test)]
pub const MAX_DEC_SCALE: usize = 1_000_000;

/// Number of working bits a correct-rounding retry loop may add to `p` before it gives up.
///
/// This is `max(WORD_BIT_SIZE * MAX_PREC_RETRY, 2 * p)`, so the loop may work at up to `3p` bits
/// (never fewer than `p + 16384`).
///
/// Why `2p`: when the exact result lies within `2^-(2p+c)` (relative) of a rounding boundary, the
/// working value has to carry about `2p + c` correct bits to fall on one side of it. Two arguments
/// of `p` bits can produce such a result without any special structure: `cos(x)` for `x` the `p`-bit
/// value nearest `4 pi` is `1 - (x - 4 pi)^2 / 2 ~ 1 - 2^-(2p-5)`; `exp(2^-p)` is the midpoint
/// `1 + 2^-p` plus `2^-(2p+1)`; `log1p(2^(1-p))` is a grid point of `p` bits plus a term of relative size `2^-2p`. A budget
/// of a fixed number of extra bits (the old `64 * 256`) is smaller than `p` for `p` above 16384 and
/// ended these in `NaN(PrecisionRetryExhausted)`. A budget of `2p` extra bits reaches `3p`, which is
/// the `2p + c` that such a result needs plus the up-to-20 % overshoot of the last growth step.
///
/// Termination: [`bump_prec_retry`] raises the working precision by at least one word per call
/// (`p_inc` starts at one word and is `round_p(p_wrk / 5)`, never below one word, afterwards), and
/// the limit is a finite number computed from `p` alone, so a loop that never certifies its result
/// (for example an exactly representable result that no finite precision can separate from its
/// boundary) fails with `PrecisionRetryExhausted` after at most `limit / 64` iterations.
#[inline]
pub fn prec_retry_budget(p: usize) -> usize {
    WORD_BIT_SIZE
        .saturating_mul(MAX_PREC_RETRY)
        .max(p.saturating_mul(2))
}

/// True when working precision has grown too far past the requested precision (see [`prec_retry_budget`]).
#[inline]
pub fn prec_retry_exhausted(p_wrk: usize, p: usize) -> bool {
    p_wrk > p.saturating_add(prec_retry_budget(p))
}

/// Bump working precision for a correct-rounding retry, or error instead of spinning.
#[inline]
pub fn bump_prec_retry(p_wrk: &mut usize, p_inc: &mut usize, p: usize) -> Result<(), Error> {
    if prec_retry_exhausted(*p_wrk, p) {
        return Err(Error::PrecisionRetryExhausted);
    }
    *p_wrk += *p_inc;
    *p_inc = round_p(*p_wrk / 5);
    Ok(())
}

#[cfg(test)]
pub const TEST_EXP_BOUND: crate::Exponent = 1024;

#[cfg(test)]
pub fn test_loop_count(full: usize) -> usize {
    full.min(200)
}

// Estimate of sqrt op cost.
#[inline]
pub fn calc_sqrt_cost(p: usize, cost_mul: usize, cost_add: usize) -> usize {
    let log3_estimate = (log2_floor(p) * 41349) >> 16;
    log3_estimate * (5 * cost_mul + 2 * cost_add) / 2
}

/// `a + b + c`. The low limb is written to `r`. The returned carry is the high limb.
///
/// Callers pass a carry of 0 or 1. For that input the result matches the old
/// `_addcarry_u64` / `_addcarry_u32` result on every limb.
#[inline(always)]
pub fn add_carry(a: Word, b: Word, c: Word, r: &mut Word) -> Word {
    let sum = a as DoubleWord + b as DoubleWord + c as DoubleWord;
    *r = sum as Word;
    (sum >> WORD_BIT_SIZE) as Word
}

/// `a - b - c`. The low limb is written to `r`. The returned borrow is 0 or 1.
///
/// Callers pass a borrow of 0 or 1. For that input the result matches the old
/// `_subborrow_u64` / `_subborrow_u32` result on every limb.
#[inline(always)]
pub fn sub_borrow(a: Word, b: Word, c: Word, r: &mut Word) -> Word {
    let v1 = a as DoubleWord;
    let v2 = b as DoubleWord + c as DoubleWord;
    *r = v1.wrapping_sub(v2) as Word;
    Word::from(v1 < v2)
}

// Shift m left by n digits.
pub fn shift_slice_left(m: &mut [Word], n: usize) {
    let idx = n / WORD_BIT_SIZE;
    let shift = n % WORD_BIT_SIZE;
    if idx >= m.len() {
        m.fill(0);
    } else if shift > 0 {
        let mut i = m.len() - 1;
        while i > idx {
            let src = i - idx;
            m[i] = (m[src] << shift) | (m[src - 1] >> (WORD_BIT_SIZE - shift));
            i -= 1;
        }
        m[idx] = m[0] << shift;
        m[..idx].fill(0);
    } else if idx > 0 {
        let r = m.len() - idx;
        m.copy_within(0..r, idx);
        m[..idx].fill(0);
    }
}

// Shift m left by n digits and put result in m2.
pub fn shift_slice_left_copy(m: &[Word], m2: &mut [Word], n: usize) {
    let idx = n / WORD_BIT_SIZE;
    let shift = n % WORD_BIT_SIZE;
    if idx >= m2.len() {
        m2.fill(0);
    } else if shift > 0 {
        m2[..idx].fill(0);
        let mut dst = m2[idx..].iter_mut();
        let src = m.iter();
        let mut prev = 0;
        for (a, b) in src.zip(dst.by_ref()) {
            *b = (prev >> (WORD_BIT_SIZE - shift)) | (*a << shift);
            prev = *a;
        }
        if let Some(b) = dst.next() {
            *b = prev >> (WORD_BIT_SIZE - shift);
        }
        for b in dst {
            *b = 0;
        }
    } else {
        m2[..idx].fill(0);
        let mut dst = m2[idx..].iter_mut();
        for (a, b) in m.iter().zip(dst.by_ref()) {
            *b = *a;
        }
        for b in dst {
            *b = 0;
        }
    }
}

// Shift m right by n digits.
pub fn shift_slice_right(m: &mut [Word], n: usize) {
    let idx = n / WORD_BIT_SIZE;
    let shift = n % WORD_BIT_SIZE;
    if idx >= m.len() {
        m.fill(0);
    } else if shift > 0 {
        let len = m.len();
        let last = len - 1;
        let mut i = 0;
        while idx + i < last {
            let src = idx + i;
            m[i] = (m[src] >> shift) | (m[src + 1] << (WORD_BIT_SIZE - shift));
            i += 1;
        }
        m[i] = m[last] >> shift;
        m[len - idx..].fill(0);
    } else if idx > 0 {
        let r = m.len() - idx;
        m.copy_within(idx.., 0);
        m[r..].fill(0);
    }
}

pub fn count_leading_zeroes_skip_first(m: &[Word]) -> usize {
    let mut iter = m.iter().rev();
    let mut w;
    let mut ret = 0;

    if let Some(v) = iter.next() {
        w = *v & (WORD_MAX >> 1);

        while w == 0 {
            ret += WORD_BIT_SIZE;

            w = match iter.next() {
                Some(v) => *v,
                None => break,
            }
        }

        if w != 0 {
            while w & WORD_SIGNIFICANT_BIT == 0 {
                w <<= 1;
                ret += 1;
            }
        }
    }

    ret
}

pub fn count_leading_ones(m: &[Word]) -> usize {
    let mut ret = 0;

    for &v in m.iter().rev() {
        if v == WORD_MAX {
            ret += WORD_BIT_SIZE;
        } else {
            let mut v = v;

            while v & WORD_SIGNIFICANT_BIT != 0 {
                v <<= 1;
                ret += 1;
            }

            break;
        }
    }

    ret
}

/// `|x|` as `usize`, saturating at `usize::MAX` where `usize` is narrower than 64 bits.
#[inline]
pub fn abs_usize(x: i64) -> usize {
    usize::try_from(x.unsigned_abs()).unwrap_or(usize::MAX)
}

/// Round precision to word bounday.
pub fn round_p(p: usize) -> usize {
    ((p.saturating_add(WORD_BIT_SIZE - 1)) / WORD_BIT_SIZE) * WORD_BIT_SIZE
}

// Convert rounding mode for an opposite sign.
pub fn invert_rm_for_sign(rm: RoundingMode) -> RoundingMode {
    if rm == RoundingMode::Up {
        RoundingMode::Down
    } else if rm == RoundingMode::Down {
        RoundingMode::Up
    } else {
        rm
    }
}

pub fn find_one_from(slice: &[Word], start_pos: usize) -> Option<usize> {
    let start_idx = start_pos / WORD_BIT_SIZE;
    if start_idx >= slice.len() {
        None
    } else {
        let mut iter = slice.iter().rev().skip(start_idx);
        if let Some(v) = iter.next() {
            let mut d = *v;

            let start_bit = start_pos % WORD_BIT_SIZE;
            let mut shift = start_pos;

            d <<= start_bit;

            if d != 0 {
                while d & WORD_SIGNIFICANT_BIT == 0 {
                    d <<= 1;
                    shift += 1;
                }

                return Some(shift);
            }

            shift += WORD_BIT_SIZE - start_bit;

            for v in iter {
                d = *v;

                if d != 0 {
                    while d & WORD_SIGNIFICANT_BIT == 0 {
                        d <<= 1;
                        shift += 1;
                    }

                    return Some(shift);
                }

                shift += WORD_BIT_SIZE;
            }
        }

        None
    }
}

/// Returns random subnormal number.
#[cfg(test)]
pub(crate) fn random_subnormal(p: usize) -> Result<ExactNumNumber, crate::Error> {
    let p = round_p(if p < 3 * WORD_BIT_SIZE { 3 * WORD_BIT_SIZE } else { p });
    let n = p - crate::common::test_rng::random::<usize>() % (2 * WORD_BIT_SIZE) - 1;
    let mut m = Vec::with_capacity(p / WORD_BIT_SIZE);

    for _ in 0..n / WORD_BIT_SIZE {
        m.push(crate::common::test_rng::random::<Word>());
    }

    if !n.is_multiple_of(WORD_BIT_SIZE) {
        let w = (crate::common::test_rng::random::<Word>() | WORD_SIGNIFICANT_BIT)
            >> (WORD_BIT_SIZE - n % WORD_BIT_SIZE);
        m.push(w);
    } else {
        *m.last_mut().ok_or(crate::Error::InvalidArgument)? |= WORD_SIGNIFICANT_BIT;
    }

    m.resize(p / WORD_BIT_SIZE, 0);

    let s = if crate::common::test_rng::random::<u8>() & 1 == 0 {
        Sign::Pos
    } else {
        Sign::Neg
    };

    ExactNumNumber::from_raw_parts(&m, n, s, EXPONENT_MIN, false)
}

#[cfg(test)]
#[inline]
pub fn rand_p() -> usize {
    crate::common::test_rng::random::<usize>() % 1000 + crate::defs::DEFAULT_P
}

#[cfg(test)]
mod carry_tests {
    use super::{add_carry, shift_slice_left, shift_slice_right, sub_borrow};
    use crate::defs::{DoubleWord, Word, WORD_BIT_SIZE, WORD_MAX};

    fn xorshift(state: &mut u64) -> u64 {
        let mut x = *state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        if x == 0 {
            x = 0x9E37_79B9_7F4A_7C15;
        }
        *state = x;
        x
    }

    fn word_from(state: &mut u64) -> Word {
        xorshift(state) as Word
    }

    /// Widening reference: `u128` for 64-bit words, `u64` for 32-bit words.
    /// The sum is the low word of `a + b + c` and the carry is the high part,
    /// so a carry of 2 survives when `c` is a full word.
    fn ref_add(a: Word, b: Word, c: Word) -> (Word, Word) {
        let wide = DoubleWord::from(a) + DoubleWord::from(b) + DoubleWord::from(c);
        ((wide >> WORD_BIT_SIZE) as Word, wide as Word)
    }

    /// Widening reference for subtraction. The difference is the low word of
    /// `a - b - c` taken modulo the double word; any high bit set is a borrow.
    fn ref_sub(a: Word, b: Word, c: Word) -> (Word, Word) {
        let wide = DoubleWord::from(a)
            .wrapping_sub(DoubleWord::from(b))
            .wrapping_sub(DoubleWord::from(c));
        (Word::from(wide >> WORD_BIT_SIZE != 0), wide as Word)
    }

    fn check_add(a: Word, b: Word, c: Word) {
        let mut got = WORD_MAX;
        let carry = add_carry(a, b, c, &mut got);
        let (want_c, want_r) = ref_add(a, b, c);
        assert_eq!(got, want_r, "sum {a}+{b}+{c}");
        assert_eq!(carry, want_c, "carry {a}+{b}+{c}");
        let wide = a as DoubleWord + b as DoubleWord + c as DoubleWord;
        let back = got as DoubleWord + ((carry as DoubleWord) << WORD_BIT_SIZE);
        assert_eq!(back, wide, "wide sum {a}+{b}+{c}");
    }

    fn check_sub(a: Word, b: Word, c: Word) {
        let mut got = 0;
        let borrow = sub_borrow(a, b, c, &mut got);
        let (want_b, want_r) = ref_sub(a, b, c);
        assert_eq!(got, want_r, "diff {a}-{b}-{c}");
        assert_eq!(borrow, want_b, "borrow {a}-{b}-{c}");
    }

    #[test]
    fn carry_and_borrow_edges() {
        let edges = [0, 1, WORD_MAX - 1, WORD_MAX];
        for &a in &edges {
            for &b in &edges {
                for c in [0, 1] {
                    check_add(a, b, c);
                    check_sub(a, b, c);
                }
            }
        }
    }

    #[test]
    fn carry_and_borrow_random_limbs() {
        let mut state = 0xC0FF_EE00_D15E_A5E5_u64;
        for _ in 0..4096 {
            let a = word_from(&mut state);
            let b = word_from(&mut state);
            let full = word_from(&mut state);
            check_add(a, b, full);
            check_sub(a, b, full);
            check_add(a, b, full & 1);
            check_sub(a, b, full & 1);
        }
    }

    fn chain_add(a: &[Word], b: &[Word]) -> (Vec<Word>, Word) {
        chain_add_in(a, b, 0)
    }

    fn chain_add_in(a: &[Word], b: &[Word], carry_in: Word) -> (Vec<Word>, Word) {
        let mut out = vec![0; a.len()];
        let mut carry = carry_in;
        for i in 0..a.len() {
            carry = add_carry(a[i], b[i], carry, &mut out[i]);
        }
        (out, carry)
    }

    fn chain_add_ref(a: &[Word], b: &[Word]) -> (Vec<Word>, Word) {
        chain_add_ref_in(a, b, 0)
    }

    fn chain_add_ref_in(a: &[Word], b: &[Word], carry_in: Word) -> (Vec<Word>, Word) {
        let mut out = vec![0; a.len()];
        let mut carry = carry_in;
        for i in 0..a.len() {
            let (next, limb) = ref_add(a[i], b[i], carry);
            out[i] = limb;
            carry = next;
        }
        (out, carry)
    }

    fn chain_sub(a: &[Word], b: &[Word]) -> (Vec<Word>, Word) {
        chain_sub_in(a, b, 0)
    }

    fn chain_sub_in(a: &[Word], b: &[Word], borrow_in: Word) -> (Vec<Word>, Word) {
        let mut out = vec![0; a.len()];
        let mut borrow = borrow_in;
        for i in 0..a.len() {
            borrow = sub_borrow(a[i], b[i], borrow, &mut out[i]);
        }
        (out, borrow)
    }

    fn chain_sub_ref(a: &[Word], b: &[Word]) -> (Vec<Word>, Word) {
        chain_sub_ref_in(a, b, 0)
    }

    fn chain_sub_ref_in(a: &[Word], b: &[Word], borrow_in: Word) -> (Vec<Word>, Word) {
        let mut out = vec![0; a.len()];
        let mut borrow = borrow_in;
        for i in 0..a.len() {
            let (next, limb) = ref_sub(a[i], b[i], borrow);
            out[i] = limb;
            borrow = next;
        }
        (out, borrow)
    }

    #[test]
    fn carry_chains_on_max_zero_and_random() {
        for len in 1..9 {
            let max = vec![WORD_MAX; len];
            let zero = vec![0; len];
            assert_eq!(chain_add(&max, &max), chain_add_ref(&max, &max));
            assert_eq!(chain_add(&max, &zero), chain_add_ref(&max, &zero));
            assert_eq!(chain_sub(&zero, &max), chain_sub_ref(&zero, &max));
            assert_eq!(chain_sub(&max, &max), chain_sub_ref(&max, &max));
            let (sum, carry) = chain_add(&max, &max);
            assert_eq!(sum[0], WORD_MAX - 1);
            assert!(sum[1..].iter().all(|&w| w == WORD_MAX));
            assert_eq!(carry, 1);
            let (diff, borrow) = chain_sub(&zero, &max);
            assert_eq!(diff[0], 1);
            assert!(diff[1..].iter().all(|&w| w == 0));
            assert_eq!(borrow, 1);
        }

        let mut state = 0x1234_5678_9ABC_DEF0_u64;
        for len in [1usize, 2, 3, 8] {
            for _ in 0..64 {
                let mut a = Vec::new();
                let mut b = Vec::new();
                for _ in 0..len {
                    a.push(word_from(&mut state));
                    b.push(word_from(&mut state));
                }
                assert_eq!(chain_add(&a, &b), chain_add_ref(&a, &b));
                assert_eq!(chain_sub(&a, &b), chain_sub_ref(&a, &b));
            }
        }
    }

    #[test]
    fn carry_chains_eight_limbs_all_max() {
        let max = [WORD_MAX; 8];
        let zero = [0; 8];
        for cin in [0, 1] {
            let sum = chain_add_in(&max, &max, cin);
            assert_eq!(sum, chain_add_ref_in(&max, &max, cin));
            assert_eq!(sum.0[0], WORD_MAX - 1 + cin);
            assert!(sum.0[1..].iter().all(|&w| w == WORD_MAX));
            assert_eq!(sum.1, 1);

            let ripple = chain_add_in(&max, &zero, cin);
            assert_eq!(ripple, chain_add_ref_in(&max, &zero, cin));
            assert_eq!(ripple, (vec![if cin == 1 { 0 } else { WORD_MAX }; 8], cin));

            let diff = chain_sub_in(&zero, &max, cin);
            assert_eq!(diff, chain_sub_ref_in(&zero, &max, cin));
            assert_eq!(diff.0[0], 1 - cin);
            assert!(diff.0[1..].iter().all(|&w| w == 0));
            assert_eq!(diff.1, 1);

            let same = chain_sub_in(&max, &max, cin);
            assert_eq!(same, chain_sub_ref_in(&max, &max, cin));
            assert_eq!(same, (vec![if cin == 1 { WORD_MAX } else { 0 }; 8], cin));
        }
    }

    fn shl1(m: &mut [Word]) {
        let mut carry = 0;
        for w in m.iter_mut() {
            let next = *w >> (WORD_BIT_SIZE - 1);
            *w = (*w << 1) | carry;
            carry = next;
        }
    }

    fn shr1(m: &mut [Word]) {
        let mut carry = 0;
        for w in m.iter_mut().rev() {
            let next = *w & 1;
            *w = (*w >> 1) | (carry << (WORD_BIT_SIZE - 1));
            carry = next;
        }
    }

    fn check_shift(src: &[Word], n: usize) {
        let dup = || {
            let mut v = Vec::with_capacity(src.len());
            v.extend(src.iter().copied());
            v
        };
        let mut got_l = dup();
        let mut want_l = dup();
        shift_slice_left(&mut got_l, n);
        for _ in 0..n {
            shl1(&mut want_l);
        }
        assert_eq!(got_l, want_l, "left {n} {src:?}");

        let mut got_r = dup();
        let mut want_r = dup();
        shift_slice_right(&mut got_r, n);
        for _ in 0..n {
            shr1(&mut want_r);
        }
        assert_eq!(got_r, want_r, "right {n} {src:?}");
    }

    #[test]
    fn shifts_match_one_bit_reference() {
        let patterns: &[&[Word]] =
            &[&[], &[0], &[1], &[WORD_MAX], &[1, 0], &[WORD_MAX, WORD_MAX], &[0, 1, WORD_MAX]];
        let widths = [0, 1, WORD_BIT_SIZE - 1, WORD_BIT_SIZE, WORD_BIT_SIZE + 1, 3 * WORD_BIT_SIZE];
        for src in patterns {
            for n in widths {
                check_shift(src, n);
            }
        }
        let mut state = 0x0F1E_2D3C_4B5A_6978_u64;
        for _ in 0..128 {
            let len = (word_from(&mut state) as usize) % 6;
            let mut src = Vec::new();
            for _ in 0..len {
                src.push(word_from(&mut state));
            }
            let n = (word_from(&mut state) as usize) % (8 * WORD_BIT_SIZE);
            check_shift(&src, n);
        }
    }
}

#[cfg(test)]
mod retry_budget_tests {
    use super::{bump_prec_retry, prec_retry_budget, prec_retry_exhausted, round_p};
    use crate::defs::{Error, WORD_BIT_SIZE};

    #[test]
    fn budget_is_the_larger_of_the_fixed_part_and_twice_p() {
        // 64 * 256 = 16384 extra bits up to p = 8192, then 2p
        assert_eq!(prec_retry_budget(64), 16384);
        assert_eq!(prec_retry_budget(8192), 16384);
        assert_eq!(prec_retry_budget(8256), 16512);
        assert_eq!(prec_retry_budget(100_032), 200_064);
        assert_eq!(prec_retry_budget(usize::MAX), usize::MAX);
    }

    #[test]
    fn exhausted_exactly_above_three_p_for_large_p() {
        let p = 100_032;
        assert!(!prec_retry_exhausted(3 * p, p));
        assert!(prec_retry_exhausted(3 * p + 1, p));
        // the old fixed budget gave up at p + 16384
        assert!(!prec_retry_exhausted(p + 16385, p));
    }

    #[test]
    fn small_p_keeps_the_fixed_budget() {
        let p = 128;
        assert!(!prec_retry_exhausted(p + 16384, p));
        assert!(prec_retry_exhausted(p + 16385, p));
    }

    #[test]
    fn a_loop_that_never_certifies_ends_in_an_error_after_a_bounded_number_of_steps() {
        for p in [64, 128, 1000, 20_032, 100_032] {
            let p = round_p(p);
            let mut p_inc = WORD_BIT_SIZE;
            let mut p_wrk = p + p_inc;
            let mut steps = 0_usize;
            let limit = p + prec_retry_budget(p);
            loop {
                match bump_prec_retry(&mut p_wrk, &mut p_inc, p) {
                    Ok(()) => steps += 1,
                    Err(e) => {
                        assert_eq!(e, Error::PrecisionRetryExhausted);
                        break;
                    }
                }
                // every step adds at least one word, so the loop cannot outlast limit / 64 steps
                assert!(steps <= limit / WORD_BIT_SIZE + 1, "p = {p}: {steps} steps");
            }
            assert!(p_wrk > limit);
            assert!(p_wrk <= limit * 6 / 5 + WORD_BIT_SIZE, "p = {p}: stopped at {p_wrk}");
        }
    }
}
