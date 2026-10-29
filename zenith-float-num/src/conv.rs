// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! Conversion utilities.

use crate::common::consts::TEN;
use crate::common::util::bump_prec_retry;
use crate::common::util::log2_ceil;
use crate::common::util::round_p;
use crate::common::util::MAX_PREC_RETRY;
use crate::defs::DoubleWord;
use crate::defs::Error;
use crate::defs::Exponent;
use crate::defs::Radix;
use crate::defs::RoundingMode;
use crate::defs::Sign;
use crate::defs::Word;
use crate::defs::DEFAULT_P;
use crate::defs::WORD_BIT_SIZE;
use crate::defs::WORD_MAX;
use crate::mantissa::Mantissa;
use crate::num::ExactNumNumber;
use crate::Consts;
use crate::EXPONENT_MAX;
use crate::EXPONENT_MIN;

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

const TEN_PWR_MAX_TO_DEC: usize = EXPONENT_MAX as usize / 4;
/// floor(EXPONENT_MAX * log10(2)) in u64: a 32-bit `usize` would truncate it.
const TEN_PWR_MAX_FROM_DEC: u64 = (EXPONENT_MAX as u128 * 301029995 / 1000000000) as u64;
/// floor(log10(2) * 2^64), a lower bound of log10(2) exact to 2^-64.
const LOG10_2_Q64: i128 = 5553023288523357132;

impl ExactNumNumber {
    /// Converts an array of digits in radix `rdx` to ExactNumNumber with precision `p`.
    /// `digits` represents mantissa and is interpreted as a number smaller than 1 and greater or equal to 1/`rdx`.
    /// The first element in `digits` is the most significant digit.
    /// `e` is the exponent part of the number, such that the number can be represented as `digits` * `rdx` ^ `e`.
    /// Precision is rounded upwards to the word size.
    /// if `p` equals usize::MAX then the precision of the resulting number is determined automatically from the input.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - ExponentOverflow: the resulting exponent becomes greater than the maximum allowed value for the exponent.
    ///  - InvalidArgument: the precision is incorrect, or `digits` contains unacceptable digits for given radix,
    ///    or when `e` is less than EXPONENT_MIN or greater than EXPONENT_MAX.
    pub fn convert_from_radix(
        sign: Sign,
        digits: &[u8],
        e: Exponent,
        rdx: Radix,
        mut p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        if p < usize::MAX {
            p = round_p(p);
            Self::p_assertion(p)?;
        }

        match rdx.value() {
            2 => Self::conv_from_binary(sign, digits, e, p, rm),
            10 => Self::conv_from_dec(sign, digits, e, p, rm, cc),
            b if b.is_power_of_two() => {
                Self::conv_from_commensurable(sign, digits, e, b.trailing_zeros() as usize, p, rm)
            }
            b => Self::conv_from_general(sign, digits, e, b, p, rm),
        }
    }

    fn conv_from_binary(
        sign: Sign,
        digits: &[u8],
        mut e: Exponent,
        p: usize,
        rm: RoundingMode,
    ) -> Result<Self, Error> {
        if digits.is_empty() {
            return Self::zero_with_sign(sign, if p < usize::MAX { p } else { DEFAULT_P });
        }

        Self::p_assertion(round_p(digits.len()))?;

        let mut mantissa = Mantissa::new(digits.len())?;
        let mut d = 0;
        let mut shift = digits.len() % WORD_BIT_SIZE;
        if shift != 0 {
            shift = WORD_BIT_SIZE - shift;
        }

        let mut dst = mantissa.digits_mut().iter_mut();

        for v in digits.iter().rev() {
            if *v > 1 {
                return Err(Error::InvalidArgument);
            }

            d |= (*v as Word) << shift;
            shift += 1;

            if shift == WORD_BIT_SIZE {
                let Some(w) = dst.next() else {
                    return Err(Error::InvalidArgument); // mantissa has precision of digits.len()
                };
                *w = d;
                shift = 0;
                d = 0;
            }
        }

        if shift > 0 {
            let Some(w) = dst.next() else {
                return Err(Error::InvalidArgument); // mantissa has precision of digits.len()
            };
            *w = d;
        }

        mantissa.update_bit_len();

        if p < usize::MAX {
            let mut ret = ExactNumNumber::from_raw_unchecked(mantissa, sign, e, false);

            ret.set_precision(p, rm)?;

            Ok(ret)
        } else {
            let shift = mantissa.max_bit_len() - mantissa.bit_len();

            if shift > 0 && e > EXPONENT_MIN {
                let e_shift = (e - EXPONENT_MIN) as usize;

                if e_shift >= shift {
                    mantissa.normilize2();
                    e -= shift as Exponent;
                } else {
                    mantissa.shift_left(e_shift);
                    mantissa.set_bit_len(mantissa.bit_len() + e_shift);
                    e = EXPONENT_MIN;
                }
            }

            Ok(ExactNumNumber::from_raw_unchecked(mantissa, sign, e, false))
        }
    }

    // radix is power of 2.
    fn conv_from_commensurable(
        sign: Sign,
        digits: &[u8],
        e: Exponent,
        shift: usize,
        p: usize,
        rm: RoundingMode,
    ) -> Result<Self, Error> {
        let significant_bit = 1 << (shift - 1);
        let base = 1 << shift;

        // exponent shift
        let mut e_shift = 0;
        let mut zeroes = 0;
        let mut first_shift = 0;
        for v in digits {
            let mut v = *v;

            if v == 0 {
                e_shift -= shift as i64;
                zeroes += 1;
            } else {
                if v >= base {
                    return Err(Error::InvalidArgument);
                }

                while v & significant_bit == 0 {
                    e_shift -= 1;
                    v <<= 1;
                    first_shift += 1;
                }
                break;
            }
        }

        if zeroes == digits.len() {
            // mantissa is zero
            let m = Mantissa::new(if p < usize::MAX { p } else { DEFAULT_P })?;

            Ok(ExactNumNumber::from_raw_unchecked(m, sign, 0, false))
        } else {
            let msz = round_p((digits.len() - zeroes) * shift + WORD_BIT_SIZE);

            Self::p_assertion(msz)?;

            let mut m = Mantissa::new(msz)?;

            // exponent
            let e = crate::exp::narrow(
                crate::exp::wide(crate::exp::mul(e, shift as i64)) + crate::exp::wide(e_shift),
            );

            if e > EXPONENT_MAX {
                return Err(Error::ExponentOverflow(sign));
            }

            // fill mantissa
            let mut filled = shift - first_shift;
            let mut dst = m.digits_mut().iter_mut().rev();
            let mut d = 0;

            'outer: for &v in digits.iter().skip(zeroes) {
                if v >= base {
                    return Err(Error::InvalidArgument);
                }

                if filled <= WORD_BIT_SIZE {
                    d |= (v as Word) << (WORD_BIT_SIZE - filled);
                } else {
                    d |= (v as Word) >> (filled - WORD_BIT_SIZE);

                    if let Some(w) = dst.next() {
                        *w = d;
                        filled -= WORD_BIT_SIZE;
                        d = if filled > 0 { (v as Word) << (WORD_BIT_SIZE - filled) } else { 0 };
                    } else {
                        break 'outer;
                    }
                }

                filled += shift;
            }

            if d > 0 {
                if let Some(w) = dst.next() {
                    *w = d;
                }
            }

            m.set_bit_len(m.max_bit_len());

            let mut ret = if e < EXPONENT_MIN {
                if p < usize::MAX {
                    let mut num = ExactNumNumber::from_raw_unchecked(m, sign, EXPONENT_MIN, false);

                    if p + WORD_BIT_SIZE > num.mantissa_max_bit_len() {
                        num.set_precision(p + WORD_BIT_SIZE, RoundingMode::None)?;
                    }

                    num.subnormalize(e, RoundingMode::None);

                    if num.inexact() {
                        num.mantissa_mut().digits_mut()[0] |= 1; // sticky for correct rounding when calling set_precision()
                    }

                    num
                } else {
                    let e_shift = (EXPONENT_MIN - e) as usize;
                    let newbitlen = m.bit_len().saturating_add(e_shift);

                    Self::p_assertion(round_p(newbitlen))?;

                    m.extend_subnormal(e_shift)?;

                    let leftshift = m.max_bit_len() - newbitlen;
                    m.shift_left(leftshift);

                    ExactNumNumber::from_raw_unchecked(m, sign, EXPONENT_MIN, false)
                }
            } else {
                ExactNumNumber::from_raw_unchecked(m, sign, e as Exponent, false)
            };

            if p < usize::MAX {
                ret.set_precision(p, rm)?;
            }

            Ok(ret)
        }
    }

    fn conv_from_general(
        sign: Sign,
        digits: &[u8],
        e: Exponent,
        base: u8,
        p: usize,
        rm: RoundingMode,
    ) -> Result<Self, Error> {
        if digits.is_empty() {
            return Self::zero_with_sign(sign, if p < usize::MAX { p } else { DEFAULT_P });
        }

        for &d in digits {
            if d >= base {
                return Err(Error::InvalidArgument);
            }
        }

        let p = if p < usize::MAX {
            round_p(p)
        } else {
            let p = round_p((digits.len() as u64 * 3321928095 / 1000000000) as usize + 1);
            Self::p_assertion(p)?;
            p
        };

        let p_wrk = p + WORD_BIT_SIZE;
        let b = Self::from_usize(base as usize)?;

        let mut acc = Self::new(p_wrk)?;
        for &d in digits {
            acc = acc.mul(&b, p_wrk, RoundingMode::None)?;
            if d > 0 {
                let dv = Self::from_usize(d as usize)?;
                acc = acc.add(&dv, p_wrk, RoundingMode::None)?;
            }
        }

        let exp_adj = e - digits.len() as i64;
        if exp_adj > 0 {
            let pwr = b.powi(exp_adj as usize, p_wrk, RoundingMode::None)?;
            acc = acc.mul(&pwr, p_wrk, RoundingMode::None)?;
        } else if exp_adj < 0 {
            let pwr = b.powi((-exp_adj) as usize, p_wrk, RoundingMode::None)?;
            acc = acc.div(&pwr, p_wrk, RoundingMode::None)?;
        }

        acc.set_sign(sign);
        acc.set_precision(p, rm)?;
        Ok(acc)
    }

    /// The integer formed by the decimal digits (no leading zeroes, non-empty), as a normalised number
    /// with the given sign.
    fn dec_mantissa_number(
        sign: Sign,
        digits: &[u8],
        tenpowers: &[(
            crate::common::buf::WordBuf,
            crate::common::buf::WordBuf,
            usize,
        )],
    ) -> Result<Self, Error> {
        let mut m = Mantissa::conv_from_dec(digits, tenpowers)?;

        if m.bit_len() > EXPONENT_MAX as usize {
            return Err(Error::ExponentOverflow(sign));
        }

        let me = m.bit_len() as Exponent;
        let _ = m.normilize2();

        Ok(ExactNumNumber::from_raw_unchecked(m, sign, me, false))
    }

    fn conv_from_dec(
        sign: Sign,
        digits: &[u8],
        e: Exponent,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        // mantissa part
        let leadzeroes = digits.iter().take_while(|&&x| x == 0).count();

        if digits.len() - leadzeroes == 0 {
            return Self::zero_with_sign(sign, if p < usize::MAX { p } else { DEFAULT_P });
        }

        let k = log2_ceil(digits.len() - leadzeroes);
        let x = Self::dec_mantissa_number(sign, &digits[leadzeroes..], cc.tenpowers(k)?)?;
        // the first pass consumes `x`; a retry rebuilds the operand from the digits (rule R2).
        let mut pending = Some(x);

        // exponent part
        let n = e - digits.len() as i64;

        let p = if p < usize::MAX {
            p
        } else {
            // determine from the input
            let p = round_p((digits.len() as u64 * 3321928095 / 1000000000) as usize + 1);
            Self::p_assertion(p)?;
            p
        };

        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p + p_inc;

        // error estimation
        let mut err = 0;
        // |n| and the chunk in u64: a 32-bit `usize` would saturate or truncate them.
        let n_abs = n.unsigned_abs();
        let npowmax = n_abs / TEN_PWR_MAX_FROM_DEC;
        let tenpowrem = n_abs % TEN_PWR_MAX_FROM_DEC;
        if npowmax != 0 {
            err += usize::try_from(npowmax)
                .ok()
                .and_then(|q| q.checked_mul(3))
                .ok_or(Error::ExponentOverflow(sign))?;
        }
        if tenpowrem != 0 {
            err += 3;
        }

        let mut retries = 0usize;
        loop {
            let p_f = p_wrk + err;

            let mut f = match pending.take() {
                Some(x) => x,
                None => Self::dec_mantissa_number(sign, &digits[leadzeroes..], cc.tenpowers(k)?)?,
            };

            if npowmax != 0 {
                let fpnmax = TEN.powi_u64(TEN_PWR_MAX_FROM_DEC, p_f, RoundingMode::None)?;

                for _ in 0..npowmax {
                    if n < 0 {
                        f = f.div(&fpnmax, p_f, RoundingMode::None)?
                    } else {
                        f = f.mul(&fpnmax, p_f, RoundingMode::None)?
                    }
                }
            };

            if tenpowrem != 0 {
                let fpn = TEN.powi_u64(tenpowrem, p_f, RoundingMode::None)?;
                if n < 0 {
                    f = f.div(&fpn, p_f, RoundingMode::None)?
                } else {
                    f = f.mul(&fpn, p_f, RoundingMode::None)?
                }
            }

            f.set_sign(sign);

            if f.try_set_precision(p, rm, p_wrk)? {
                return Ok(f);
            }

            retries += 1;
            if retries > MAX_PREC_RETRY {
                f.set_precision(p, rm)?;
                return Ok(f);
            }

            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
            if p_wrk > EXPONENT_MAX as usize {
                return Err(Error::InvalidArgument);
            }
        }
    }

    /// Converts `self` to radix `rdx` using rounding mode `rm`.
    /// The function returns sign, mantissa digits in radix `rdx`, and exponent such that the converted number
    /// can be represented as `mantissa digits` * `rdx` ^ `exponent`.
    /// The first element in the mantissa is the most significant digit.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - ExponentOverflow: the resulting exponent becomes greater than the maximum allowed value for the exponent.
    pub fn convert_to_radix(
        &self,
        rdx: Radix,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<(Sign, Vec<u8>, Exponent), Error> {
        match rdx.value() {
            2 => self.conv_to_binary(),
            10 => self.conv_to_dec(rm, cc, 0),
            b if b.is_power_of_two() => self.conv_to_commensurable(b.trailing_zeros() as usize),
            b => self.conv_to_general(b, rm, cc),
        }
    }

    fn zero_with_sign(sign: Sign, p: usize) -> Result<Self, Error> {
        let mut z = Self::new(p)?;
        z.set_sign(sign);
        Ok(z)
    }

    pub(crate) fn convert_to_dec_extra(
        &self,
        rm: RoundingMode,
        cc: &mut Consts,
        extra_digits: usize,
    ) -> Result<(Sign, Vec<u8>, Exponent), Error> {
        self.conv_to_dec(rm, cc, extra_digits)
    }

    fn conv_to_dec(
        &self,
        rm: RoundingMode,
        cc: &mut Consts,
        extra_digits: usize,
    ) -> Result<(Sign, Vec<u8>, Exponent), Error> {
        if self.precision() == 0 {
            return Ok((self.sign(), Vec::new(), 0));
        }

        let p = self.mantissa_max_bit_len();
        let subn_e = p - self.precision();
        let dec_len = (p as u64 * 301029996 / 1000000000) as usize + 1 + extra_digits;

        let mut err = WORD_BIT_SIZE; // speculative
        let mut p_wrk = round_p((dec_len as u64 * 3321928095 / 1000000000) as usize + 1 + err);
        let mut p_inc = WORD_BIT_SIZE;

        let mut retries = 0usize;
        loop {
            let mut x = self.clone()?;
            x.set_inexact(false);

            // i128 with floor(log10(2) * 2^64): with an i64 exponent the old
            // 301029996e-9 factor is off by ~5e8 digits and the product overflows i64.
            let n_wrk = i64::try_from(
                (((p_wrk as i128 - i128::from(self.exponent()) + subn_e as i128) * LOG10_2_Q64)
                    >> 64)
                    + 1,
            )
            .map_err(|_| Error::ExponentOverflow(self.sign()))?;

            let mut err_acc = 0;

            // |n_wrk| in u64: a 32-bit `usize` would saturate it. What is left after the chunk loop is at most
            // the chunk, a `usize`, so its conversion back cannot fail.
            let chunk = u64::try_from(TEN_PWR_MAX_TO_DEC).map_err(|_| Error::InvalidArgument)?;
            let mut pwr = n_wrk.unsigned_abs();
            if pwr > chunk {
                let tp = TEN.powi(TEN_PWR_MAX_TO_DEC, p_wrk, RoundingMode::None)?;
                err_acc += 1;

                while pwr > chunk {
                    if n_wrk < 0 {
                        x = x.div(&tp, p_wrk, RoundingMode::None)?;
                    } else {
                        x = x.mul(&tp, p_wrk, RoundingMode::None)?;
                    }
                    err_acc += 2;
                    pwr -= chunk;
                }
            }

            if pwr != 0 {
                let pwr = usize::try_from(pwr).map_err(|_| Error::InvalidArgument)?;
                let tp = TEN.powi(pwr, p_wrk, RoundingMode::None)?;
                if n_wrk < 0 {
                    x = x.div(&tp, p_wrk, RoundingMode::None)?;
                } else {
                    x = x.mul(&tp, p_wrk, RoundingMode::None)?;
                }
                err_acc += 3;
            }

            if err_acc > err {
                err_acc += err_acc / TEN_PWR_MAX_TO_DEC + 3;
                p_wrk += round_p(err_acc - err);
                err = err_acc;
                continue;
            }

            let (m_slice, _m_bits, _s, e, inexact) = x.as_raw_parts();
            let mut m = Mantissa::from_raw_parts(m_slice, _m_bits)?;

            let shift = e.saturating_sub(p_wrk as i64).max(0) as usize;
            if shift > 0 {
                m.shift_left_resize(shift)?;
            }

            let l = (m.bit_len() as u64 * 301029996 / 1000000000) as usize + 1;

            let k = log2_ceil(l);

            let tenpowers = cc.tenpowers(k)?;
            let mut digits = m.conv_to_dec(1 << k, tenpowers, k.saturating_sub(1), true)?;

            let mut e_out = digits.len() as i64 - n_wrk;

            let mut e_subn = 0;
            if e_out < EXPONENT_MIN {
                e_subn = (EXPONENT_MIN - e_out) as usize;
                e_out = EXPONENT_MIN;
            }

            let mut e_out = e_out as Exponent;

            // cut off digits with error
            let valid =
                digits.len() as i64 - ((shift + err_acc) as i64 * 301029996 / 1000000000) - 1;

            if digits.len() > dec_len && valid > dec_len as i64 {
                // try round
                if Self::try_round_dec(
                    &mut digits[..valid as usize],
                    dec_len,
                    rm,
                    self.sign(),
                    &mut e_out,
                    inexact,
                )? {
                    if e_subn > 0 {
                        if e_out > EXPONENT_MIN {
                            e_subn -= 1;
                            e_out = EXPONENT_MIN;
                        }

                        let rsz = if digits.len() > dec_len { dec_len } else { digits.len() };

                        digits.resize(rsz + e_subn, 0);
                        digits[rsz..].fill(0);
                        digits.rotate_right(e_subn);
                    } else {
                        if digits.len() > dec_len {
                            digits.resize(dec_len, 0);
                        }
                    }

                    // remove trailing zeroes
                    let nzr = digits.iter().rev().take_while(|&&x| x == 0).count();

                    digits.resize(digits.len() - nzr, 0);

                    return Ok((self.sign(), digits, e_out));
                }
            }

            retries += 1;
            if retries > MAX_PREC_RETRY {
                return Err(Error::InvalidArgument);
            }
            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
            if p_wrk > EXPONENT_MAX as usize {
                return Err(Error::InvalidArgument);
            }
        }
    }

    // Try to round a decimal mantissa.
    fn try_round_dec(
        digits: &mut [u8],
        n: usize,
        rm: RoundingMode,
        s: Sign,
        e: &mut Exponent,
        inexact: bool,
    ) -> Result<bool, Error> {
        let mut check_roundable = inexact;

        if n > 0 {
            let ovf = Self::round_dec(digits, n, rm, s.is_positive(), &mut check_roundable);

            if check_roundable {
                return Ok(false);
            }

            if ovf {
                if *e == EXPONENT_MAX {
                    return Err(Error::ExponentOverflow(s));
                }

                *e += 1;
            }
        }

        Ok(true)
    }

    // Round decimal mantissa.
    // The function is similar to Mantissa::round_mantissa.
    fn round_dec(
        digits: &mut [u8],
        n: usize,
        rm: RoundingMode,
        is_positive: bool,
        check_roundable: &mut bool,
    ) -> bool {
        if rm == RoundingMode::None {
            *check_roundable = false;
            return false;
        }

        #[inline]
        fn get_rem(arr: &[u8]) -> (bool, bool) {
            let mut rem9 = true;
            let mut rem0 = true;

            for &d in arr.iter() {
                if d != 9 {
                    rem9 = false;
                }
                if d != 0 {
                    rem0 = false;
                }
            }
            (rem0, rem9)
        }

        if n > 0 && n < digits.len() {
            let mut c = false;

            if rm == RoundingMode::ToEven || rm == RoundingMode::ToOdd {
                let is_even = digits[n - 1].is_multiple_of(2);
                let dn = digits[n];

                let (rem0, rem9) = get_rem(&digits[n + 1..]);

                if *check_roundable && (rem0 || rem9) {
                    return false;
                }

                // need adding 1?
                match rm {
                    RoundingMode::ToEven => {
                        if dn == 5 {
                            if !is_even || !rem0 {
                                c = true;
                            }
                        } else if dn > 5 {
                            c = true;
                        }
                    }
                    // Nearest, ties to odd (not sticky round-to-odd): only the exact tie (`dn == 5`, nothing
                    // after it) looks at parity.
                    RoundingMode::ToOdd => {
                        if dn == 5 {
                            if is_even || !rem0 {
                                c = true;
                            }
                        } else if dn > 5 {
                            c = true;
                        }
                    }
                    RoundingMode::None
                    | RoundingMode::Up
                    | RoundingMode::Down
                    | RoundingMode::ToZero
                    | RoundingMode::FromZero => {}
                };
            } else {
                let (rem0, rem9) = get_rem(&digits[n..]);

                if *check_roundable && (rem0 || rem9) {
                    return false;
                }

                // rounding
                match rm {
                    RoundingMode::ToZero => {}
                    RoundingMode::FromZero => {
                        if !rem0 {
                            // add 1
                            c = true;
                        }
                    }
                    RoundingMode::Up => {
                        if !rem0 && is_positive {
                            // add 1
                            c = true;
                        }
                    }
                    RoundingMode::Down => {
                        if !rem0 && !is_positive {
                            // add 1
                            c = true;
                        }
                    }
                    RoundingMode::None | RoundingMode::ToEven | RoundingMode::ToOdd => {}
                };
            }

            *check_roundable = false; // can round

            digits[n..].fill(0);

            if c {
                for v in digits[..n].iter_mut().rev() {
                    if *v < 9 {
                        *v += 1;
                        return false;
                    } else {
                        *v = 0;
                    }
                }

                digits[0] = 1;

                return true;
            }
        }
        false
    }

    /// Conversion for radixes of power of 2.
    fn conv_to_commensurable(&self, shift: usize) -> Result<(Sign, Vec<u8>, Exponent), Error> {
        let mut e = self.exponent();
        let mut e_shift = e.unsigned_abs() as usize % shift;
        e /= shift as Exponent;
        if e_shift != 0 && self.exponent() > 0 {
            e_shift = shift - e_shift;
            e += 1;
        }

        let mut ret = Vec::new();
        ret.try_reserve_exact(self.mantissa_max_bit_len() / shift)?;

        let mask = (WORD_MAX >> (WORD_BIT_SIZE - shift)) as DoubleWord;
        let mut iter = self.mantissa().digits().iter().rev();

        let mut done = WORD_BIT_SIZE - shift + e_shift;
        let Some(&top) = iter.next() else {
            return Err(Error::InvalidArgument); // a finite nonzero value has at least one word
        };
        let mut d = top as DoubleWord;
        let mut cnt = 0;

        loop {
            let digit = ((d >> done) & mask) as u8;

            if digit == 0 {
                cnt += 1;
            } else {
                cnt = 0;
            }

            ret.push(digit);

            if done < shift {
                d <<= WORD_BIT_SIZE;

                if let Some(v) = iter.next() {
                    d |= *v as DoubleWord;
                } else {
                    break;
                }

                done += WORD_BIT_SIZE;
            }

            done -= shift;
        }

        if done > 0 {
            if done < shift {
                done += WORD_BIT_SIZE - shift;
            }

            let digit = ((d >> done) & mask) as u8;
            if digit > 0 {
                cnt = 0;
                ret.push(digit);
            }
        }

        ret.resize(ret.len() - cnt, 0);

        Ok((self.sign(), ret, e))
    }

    fn conv_to_general(
        &self,
        base: u8,
        _rm: RoundingMode,
        _cc: &mut Consts,
    ) -> Result<(Sign, Vec<u8>, Exponent), Error> {
        if self.is_zero() {
            let mut zero = Vec::new();
            zero.try_reserve_exact(1)?;
            zero.push(0);
            return Ok((self.sign(), zero, 0));
        }

        let p = self.mantissa_max_bit_len() + WORD_BIT_SIZE * 2;
        let b = Self::from_usize(base as usize)?;
        let one = Self::from_usize(1)?;
        let mut x = self.clone()?;
        x.set_precision(p, RoundingMode::None)?;

        let mut e: Exponent = 0;
        while x.cmp(&one) >= 0 {
            x = x.div(&b, p, RoundingMode::None)?;
            e = e.saturating_add(1);
        }

        let mut digits = Vec::new();
        let max_digits = p * 2;
        digits.try_reserve_exact(max_digits)?;

        for _ in 0..max_digits {
            x = x.mul(&b, p, RoundingMode::None)?;
            let mut digit = 0u8;
            for d in 1..base {
                let dv = Self::from_usize(d as usize)?;
                if dv.cmp(&x) > 0 {
                    break;
                }
                digit = d;
            }
            digits.push(digit);
            if digit > 0 {
                let sub = Self::from_usize(digit as usize)?;
                x = x.sub(&sub, p, RoundingMode::None)?;
            }
            if x.is_zero() {
                break;
            }
        }

        while digits.last() == Some(&0) {
            digits.pop();
        }
        if digits.is_empty() {
            digits.push(0);
        }

        Ok((self.sign(), digits, e))
    }

    fn conv_to_binary(&self) -> Result<(Sign, Vec<u8>, Exponent), Error> {
        let mut ret = Vec::new();
        ret.try_reserve_exact(self.mantissa_max_bit_len())?;
        let mut cnt = 0;

        for v in self.mantissa().digits().iter().rev() {
            for i in (0..WORD_BIT_SIZE).rev() {
                let val = ((v >> i) & 1) as u8;
                if val == 1 {
                    cnt = 0;
                } else {
                    cnt += 1;
                }
                ret.push(val);
            }
        }

        // cut off trailing zeroes
        ret.resize(ret.len() - cnt, 0);

        Ok((self.sign(), ret, self.exponent()))
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::common::test_rng::random;
    use crate::common::util::random_subnormal;
    use crate::common::util::test_loop_count;
    use crate::common::util::MAX_DEC_SCALE;
    use crate::common::util::TEST_EXP_BOUND;
    use crate::defs::{Sign, EXPONENT_MAX, EXPONENT_MIN};
    use crate::WORD_SIGNIFICANT_BIT;

    /// The from-decimal exponent chunk is floor(EXPONENT_MAX * 0.301029995) on every target.
    /// A 32-bit `usize` truncated it mod 2^32, so |n| near EXPONENT_MAX / log2(10) split into
    /// ~10^8 chunks: ~10^8 divisions at ~3.5 * 10^8 extra working bits (an effective hang).
    #[test]
    fn dec_exponent_chunk_not_truncated() -> Result<(), crate::Error> {
        let want = u128::try_from(EXPONENT_MAX).map_err(|_| crate::Error::InvalidArgument)?
            * 301029995
            / 1000000000;
        let have = u128::try_from(TEN_PWR_MAX_FROM_DEC).map_err(|_| crate::Error::InvalidArgument)?;
        assert_eq!(have, want);
        Ok(())
    }

    /// The to-decimal exponent chunk is EXPONENT_MAX / 4 on every target. A 32-bit `usize`
    /// truncated it to 2^30 - 1, so an exponent near EXPONENT_MIN / log2(10) took ~1.4 * 10^9
    /// chunk multiplications (an effective hang).
    #[test]
    fn dec_output_exponent_chunk_not_truncated() -> Result<(), crate::Error> {
        let want = u128::from(EXPONENT_MAX.unsigned_abs()) / 4;
        let have = u128::try_from(TEN_PWR_MAX_TO_DEC).map_err(|_| crate::Error::InvalidArgument)?;
        assert_eq!(have, want);
        Ok(())
    }

    #[cfg(not(feature = "std"))]
    use alloc::vec;

    /// mpmath-pinned (g exponent, g words, n exponent, n words) for the decimal vectors of
    /// the "unknown p" round trip in `test_conv_num`, in order. Words are 64-bit, little end first.
    type DecPin = (Exponent, &'static [u64], Exponent, &'static [u64]);
    const CONV_DEC_PINS: [DecPin; 5] = [
        (
            409,
            &[
                9277499933223628938,
                3922711689764733432,
                5698786213872067458,
                13358134941265313276,
            ],
            409,
            &[
                9277499933223628938,
                3922711689764733432,
                5698786213872067458,
                13358134941265313276,
            ],
        ),
        (
            -1531968934941318218,
            &[
                12118989213338987751,
                8129038482058839168,
                957869615584094226,
                14375337014266739808,
            ],
            -1531968934941318218,
            &[
                12118989213338987750,
                8129038482058839168,
                957869615584094226,
                14375337014266739808,
            ],
        ),
        (
            -1531968934941318391,
            &[
                14376713304243956606,
                3717668849003244565,
                4709716642905670367,
                6758455429533423369,
                13654437240512262384,
                863961003625122258,
                17211046775631358814,
            ],
            -1531968934941318391,
            &[
                14376713304243956605,
                3717668849003244565,
                4709716642905670367,
                6758455429533423369,
                13654437240512262384,
                863961003625122258,
                17211046775631358814,
            ],
        ),
        (
            -1531968934941317949,
            &[
                8589000564600802957,
                953034597986637512,
                11784980727772288842,
                1040538886271782467,
                15796102105945313841,
                13627354557981012237,
                15154764830504843509,
            ],
            -1531968934941317949,
            &[
                8589000564600802957,
                953034597986637512,
                11784980727772288842,
                1040538886271782467,
                15796102105945313841,
                13627354557981012237,
                15154764830504843509,
            ],
        ),
        (
            236,
            &[
                5812446391292657664,
                17568215927531755994,
                9442323477031879346,
                15993189243573370911,
            ],
            236,
            &[
                5812446391292657664,
                17568215927531755994,
                9442323477031879346,
                15993189243573370911,
            ],
        ),
    ];

    /// Mantissa words equal the pinned 64-bit words (32-bit targets pair up low/high words).
    fn words_eq_u64(d: &[Word], pin: &[u64]) -> bool {
        #[cfg(not(target_pointer_width = "32"))]
        {
            d == pin
        }
        #[cfg(target_pointer_width = "32")]
        {
            let (pairs, rest) = d.as_chunks::<2>();
            rest.is_empty()
                && pairs.len() == pin.len()
                && pairs
                    .iter()
                    .zip(pin)
                    .all(|([lo, hi], p)| (u64::from(*hi) << 32 | u64::from(*lo)) == *p)
        }
    }

    #[test]
    fn test_conv_num() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;

        // Decimal I/O uses software big-float conversion only (no hardware floats).

        #[cfg(target_pointer_width = "32")]
        {
            let n = ExactNumNumber::from_raw_parts(
                &[2576980377, 2576980377, 2576980377],
                96,
                Sign::Pos,
                -1,
                false,
            )?;
            let (s, m, e) = n.convert_to_radix(Radix::Oct, RoundingMode::None, &mut cc)?;
            let g = ExactNumNumber::convert_from_radix(
                s,
                &m,
                e,
                Radix::Oct,
                160,
                RoundingMode::ToEven,
                &mut cc,
            )?;

            assert!(n.cmp(&g) == 0);

            let n = ExactNumNumber::from_raw_parts(
                &[2576980377, 2576980377, 2576980377],
                96,
                Sign::Pos,
                -0,
                false,
            )?;
            let (s, m, e) = n.convert_to_radix(Radix::Dec, RoundingMode::None, &mut cc)?;

            assert_eq!(
                m,
                [
                    5, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
                    9, 9, 9
                ]
            );
            assert!(s == Sign::Pos);
            assert!(e == 0);

            let g = ExactNumNumber::convert_from_radix(
                s,
                &m,
                e,
                Radix::Dec,
                96,
                RoundingMode::ToEven,
                &mut cc,
            )?;
            assert!(g.cmp(&n) == 0);
        }

        #[cfg(not(target_pointer_width = "32"))]
        {
            let n = ExactNumNumber::from_raw_parts(
                &[0x9999999999999999, 0x9999999999999999, 0x9999999999999999],
                192,
                Sign::Pos,
                -1,
                false,
            )?;
            let (s, m, e) = n.convert_to_radix(Radix::Oct, RoundingMode::None, &mut cc)?;
            let g = ExactNumNumber::convert_from_radix(
                s,
                &m,
                e,
                Radix::Oct,
                192,
                RoundingMode::ToEven,
                &mut cc,
            )?;

            assert!(n.cmp(&g) == 0);

            let n = ExactNumNumber::from_raw_parts(
                &[0x9999999999999999, 0x9999999999999999, 0x9999999999999999],
                192,
                Sign::Pos,
                -0,
                false,
            )?;
            let (s, m, e) = n.convert_to_radix(Radix::Dec, RoundingMode::ToEven, &mut cc)?;

            assert_eq!(
                m,
                [
                    5, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
                    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
                    9, 9, 9, 9, 9, 9
                ]
            );
            assert_eq!(s, Sign::Pos);
            assert_eq!(e, 0);

            let g = ExactNumNumber::convert_from_radix(
                s,
                &m,
                e,
                Radix::Dec,
                192,
                RoundingMode::ToEven,
                &mut cc,
            )?;
            //println!("{:?}", g);
            //println!("{:?}", n);
            assert!(g.cmp(&n) == 0);
        }

        // random normal values
        let mut eps = ExactNumNumber::small_const(1);
        let p_rng = 32;

        for _ in 0..test_loop_count(1000) {
            let p1 = (random::<usize>() % p_rng + 1) * WORD_BIT_SIZE;
            let p2 = (random::<usize>() % p_rng + 1) * WORD_BIT_SIZE;
            let p = p1.min(p2);

            let mut n = ExactNumNumber::random_normal(p1, -TEST_EXP_BOUND, TEST_EXP_BOUND)?;
            let rdx = random_radix();

            let (s1, m1, e1) = n.convert_to_radix(rdx, RoundingMode::ToEven, &mut cc)?;
            let mut g = ExactNumNumber::convert_from_radix(
                s1,
                &m1,
                e1,
                rdx,
                p2,
                RoundingMode::ToEven,
                &mut cc,
            )?;

            //println!("\n{:?}", rdx);
            //println!("{:?} {:?} {}", s1, m1, e1);
            //println!("{:?}\n{:?}", n, g);

            if rdx == Radix::Dec {
                eps.set_exponent(n.exponent() - p as Exponent + 4);
                assert!(n.sub(&g, p, RoundingMode::None)?.abs()?.cmp(&eps) <= 0);
            } else {
                if p2 < p1 {
                    n.set_precision(p, RoundingMode::ToEven)?;
                } else if p2 > p1 {
                    g.set_precision(p, RoundingMode::ToEven)?;
                }

                assert!(n.cmp(&g) == 0);
            }
        }

        // subnormal values
        for _ in 0..test_loop_count(1000) {
            let p1 = (random::<usize>() % p_rng + 3) * WORD_BIT_SIZE;
            let p2 = (random::<usize>() % p_rng + 3) * WORD_BIT_SIZE;
            let p = p1.min(p2);

            let mut n = random_subnormal(p1)?;
            let rdx = random_radix();
            if rdx == Radix::Dec && n.exponent().unsigned_abs() as usize > MAX_DEC_SCALE {
                continue;
            }

            let (s1, m1, e1) = n.convert_to_radix(rdx, RoundingMode::ToEven, &mut cc)?;

            //println!("\n{:?}", rdx);
            //println!("{:?} {:?} {}", s1, m1, e1);

            let mut g = ExactNumNumber::convert_from_radix(
                s1,
                &m1,
                e1,
                rdx,
                p2,
                RoundingMode::ToEven,
                &mut cc,
            )?;

            //println!("\n{:?}", rdx);
            //println!("{:?} {:?} {}", s1, m1, e1);
            //println!("{:?}\n{:?}", n, g);

            if rdx == Radix::Dec {
                let mut eps = ExactNumNumber::min_positive(p)?;
                eps.set_exponent(eps.exponent() + 1);

                assert!(
                    n.sub(&g, p, RoundingMode::None)?.abs()?.cmp(&eps) <= 0,
                    "{:?} {:?}",
                    n,
                    g
                );
            } else {
                if p2 < p1 {
                    n.set_precision(p, RoundingMode::ToEven)?;
                } else if p2 > p1 {
                    g.set_precision(p, RoundingMode::ToEven)?;
                }

                assert!(n.cmp(&g) == 0, "{:?} {:?} {:?}", rdx, n, g);
            }
        }

        // MIN, MAX, min_subnormal
        let p1 = (random::<usize>() % p_rng + 1) * WORD_BIT_SIZE;
        let p2 = (random::<usize>() % p_rng + 1) * WORD_BIT_SIZE;
        let p = p1.min(p2);

        for rdx in [Radix::Bin, Radix::Oct, Radix::Dec, Radix::Hex] {
            if rdx == Radix::Dec {
                continue;
            }
            // min, max
            // for p2 < p1 rounding will cause overflow, for p2 >= p1 no rounding is needed.
            let rm = RoundingMode::None;
            for mut n in [ExactNumNumber::max_value(p1)?, ExactNumNumber::min_value(p1)?] {
                //println!("\n{:?} {} {}", rdx, p1, p2);
                //println!("{:?}", n);

                let (s1, m1, e1) = n.convert_to_radix(rdx, rm, &mut cc)?;

                //println!("{:?} {:?} {}", s1, m1, e1);

                let mut g = ExactNumNumber::convert_from_radix(s1, &m1, e1, rdx, p2, rm, &mut cc)?;

                //println!("{:?}", g);

                if rdx == Radix::Dec {
                    eps.set_exponent(n.exponent() - p as Exponent + 4);
                    assert!(n.sub(&g, p, rm)?.abs()?.cmp(&eps) <= 0);
                } else {
                    if p2 < p1 {
                        n.set_precision(p, rm)?;
                    } else if p2 > p1 {
                        g.set_precision(p, rm)?;
                    }

                    assert!(n.cmp(&g) == 0);
                }
            }

            // min subnormal
            let rm = RoundingMode::ToEven;
            let mut n = ExactNumNumber::min_positive(p1)?;
            //println!("\n{:?} {} {}", rdx, p1, p2);
            //println!("{:?}", n);
            let (s1, m1, e1) = n.convert_to_radix(rdx, rm, &mut cc)?;

            //println!("{:?} {:?} {}", s1, m1, e1);

            let mut g = ExactNumNumber::convert_from_radix(s1, &m1, e1, rdx, p2, rm, &mut cc)?;
            //println!("{:?}", g);

            if rdx == Radix::Dec {
                let mut eps = ExactNumNumber::min_positive(p)?;
                eps.set_exponent(eps.exponent() + 1);
                assert!(n.sub(&g, p, rm)?.abs()?.cmp(&eps) <= 0);
            } else {
                if p2 < p1 {
                    n.set_precision(p, rm)?;
                } else if p2 > p1 {
                    g.set_precision(p, rm)?;
                }
                assert!(n.cmp(&g) == 0);
            }
        }

        // misc/invalid input
        let s1 = Sign::Pos;
        for rdx in [Radix::Bin, Radix::Oct, Radix::Dec, Radix::Hex] {
            for e1 in [123, -123, 0] {
                let m1 = [];
                assert!(ExactNumNumber::convert_from_radix(
                    s1,
                    &m1,
                    e1,
                    rdx,
                    p1,
                    RoundingMode::ToEven,
                    &mut cc
                )?
                .is_zero());
                let m1 = [1, rdx.value(), 0];
                assert!(matches!(
                    ExactNumNumber::convert_from_radix(
                        s1,
                        &m1,
                        e1,
                        rdx,
                        p1,
                        RoundingMode::ToEven,
                        &mut cc
                    ),
                    Err(Error::InvalidArgument)
                ));
                // `p = 0` is rejected (2.0: `MIN_PRECISION` is 1).
                let m1 = [1, rdx.value() - 1, 0];
                assert!(matches!(
                    ExactNumNumber::convert_from_radix(
                        s1,
                        &m1,
                        e1,
                        rdx,
                        0,
                        RoundingMode::ToEven,
                        &mut cc,
                    ),
                    Err(Error::InvalidArgument)
                ));
                let m1 = [0; 256];
                assert!(ExactNumNumber::convert_from_radix(
                    s1,
                    &m1,
                    e1,
                    rdx,
                    p1,
                    RoundingMode::ToEven,
                    &mut cc,
                )?
                .is_zero());
                let m1 = [0];
                assert!(ExactNumNumber::convert_from_radix(
                    s1,
                    &m1,
                    e1,
                    rdx,
                    p1,
                    RoundingMode::ToEven,
                    &mut cc,
                )?
                .is_zero());
            }
        }

        // dec, short digits
        let n = ExactNumNumber::from_words(
            &[1052139549, 0, 0, 0, 0, 0, 0, 0, WORD_SIGNIFICANT_BIT],
            Sign::Pos,
            1,
        )?;
        let (s, d, e) = n.convert_to_radix(Radix::Dec, RoundingMode::ToEven, &mut cc)?;
        let g = ExactNumNumber::convert_from_radix(
            s,
            &d,
            e,
            Radix::Dec,
            n.mantissa_max_bit_len(),
            RoundingMode::None,
            &mut cc,
        )?;
        assert!(
            g.cmp(&n) == 0
                || g.sub(&n, n.mantissa_max_bit_len(), RoundingMode::ToEven)?
                    .exponent()
                    < -(n.mantissa_max_bit_len() as Exponent - 4)
        );

        // unkonwn p: binary, octal, hexadecimal
        let mut dec_ix = 0usize;
        for (s1, exp, rdx) in [
            // Bin
            (
                vec![
                    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1,
                    1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1,
                    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1,
                ],
                123,
                Radix::Bin,
            ),
            (
                vec![
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0,
                    0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1,
                    1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1,
                ],
                EXPONENT_MIN + 3,
                Radix::Bin,
            ),
            (
                vec![
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0,
                    0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1,
                    1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1,
                ],
                EXPONENT_MIN,
                Radix::Bin,
            ),
            (
                vec![
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0,
                    0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1,
                    1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1,
                ],
                EXPONENT_MIN + 256,
                Radix::Bin,
            ),
            // Oct
            (
                vec![
                    7, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1,
                    1, 0, 1, 0, 0, 0, 4, 0, 0, 0, 0, 1, 5, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1,
                    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 6,
                ],
                123,
                Radix::Oct,
            ),
            (
                vec![
                    7, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1,
                    1, 0, 1, 0, 0, 0, 4, 0, 0, 0, 0, 1, 5, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1,
                    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 6,
                ],
                EXPONENT_MIN / 3 - 123,
                Radix::Oct,
            ),
            (
                vec![
                    3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1,
                    1, 0, 1, 0, 0, 0, 4, 0, 0, 0, 0, 1, 5, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1,
                    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 6,
                ],
                EXPONENT_MIN / 3 + 123,
                Radix::Oct,
            ),
            (
                vec![
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 1, 1, 0, 0, 0, 1, 1, 1,
                    1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 4, 0, 0, 0, 0, 1, 5, 1, 0, 0, 0, 0, 1, 1, 1, 0,
                    0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 6,
                ],
                EXPONENT_MIN / 3 - 123,
                Radix::Oct,
            ),
            (
                vec![
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 3, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 1, 1, 0, 0, 0, 1, 1, 1,
                    1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 4, 0, 0, 0, 0, 1, 5, 1, 0, 0, 0, 0, 1, 1, 1, 0,
                    0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 6,
                ],
                EXPONENT_MIN / 3 + 123,
                Radix::Oct,
            ),
            // Hex
            (
                vec![
                    0xF, 0, 2, 0, 0, 0, 0, 0, 0xD, 0, 0, 0, 4, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1,
                    1, 0, 1, 0, 0, 0, 4, 0, 0, 0, 0, 1, 5, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1,
                    1, 0, 0, 0, 0, 0, 0xA, 0, 0, 0, 0, 0, 3, 0, 6,
                ],
                123,
                Radix::Hex,
            ),
            (
                vec![
                    0xF, 0, 2, 0, 0, 0, 0, 0, 0xD, 0, 0, 0, 4, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1,
                    1, 0, 1, 0, 0, 0, 4, 0, 0, 0, 0, 1, 5, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1,
                    1, 0, 0, 0, 0, 0, 0xA, 0, 0, 0, 0, 0, 3, 0, 6,
                ],
                EXPONENT_MIN / 4 - 123,
                Radix::Hex,
            ),
            (
                vec![
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    7, 0, 2, 0, 0, 0, 0, 0, 0xD, 0, 0, 0, 4, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1,
                    0, 1, 0, 0, 0, 4, 0, 0, 0, 0, 1, 5, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1,
                    0, 0, 0, 0, 0, 0xA, 0, 0, 0, 0, 0, 3, 0, 6,
                ],
                EXPONENT_MIN / 4 - 123,
                Radix::Hex,
            ),
            (
                vec![
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    7, 0, 2, 0, 0, 0, 0, 0, 0xD, 0, 0, 0, 4, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1,
                    0, 1, 0, 0, 0, 4, 0, 0, 0, 0, 1, 5, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1,
                    0, 0, 0, 0, 0, 0xA, 0, 0, 0, 0, 0, 3, 0, 6,
                ],
                EXPONENT_MIN / 4 + 123,
                Radix::Hex,
            ),
            // Dec
            (
                vec![
                    9, 5, 7, 4, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 9, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1,
                    1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 9, 1, 1, 7, 0, 0, 0, 0, 1, 1, 1, 5, 0, 0, 1, 1,
                    1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 2, 0, 6,
                ],
                123,
                Radix::Dec,
            ),
            (
                vec![
                    9, 5, 7, 4, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 9, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1,
                    1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 9, 1, 1, 7, 0, 0, 0, 0, 1, 1, 1, 5, 0, 0, 1, 1,
                    1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 2, 0, 5,
                ],
                EXPONENT_MIN / 10 - 123,
                Radix::Dec,
            ),
            (
                vec![
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    9, 5, 7, 4, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 9, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1,
                    1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 9, 1, 1, 7, 0, 0, 0, 0, 1, 1, 1, 5, 0, 0, 1, 1,
                    1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 2, 0, 7,
                ],
                EXPONENT_MIN / 10 - 123,
                Radix::Dec,
            ),
            (
                vec![
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    9, 5, 7, 4, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 9, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1,
                    1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 9, 1, 1, 7, 0, 0, 0, 0, 1, 1, 1, 5, 0, 0, 1, 1,
                    1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 2, 0, 5,
                ],
                EXPONENT_MIN / 10 + 10,
                Radix::Dec,
            ),
            (
                vec![
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    9, 5, 7, 4, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 9, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1,
                    1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 9, 1, 1, 7, 0, 0, 0, 0, 1, 1, 1, 5, 0, 0, 1, 1,
                    1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 2, 0, 5,
                ],
                123,
                Radix::Dec,
            ),
        ] {
            let mut g = ExactNumNumber::convert_from_radix(
                Sign::Pos,
                &s1,
                exp,
                rdx,
                usize::MAX,
                RoundingMode::None,
                &mut cc,
            )?;

            let (s, d, e) = g.convert_to_radix(rdx, RoundingMode::None, &mut cc)?;

            // verify no information loss in mantissa becase of length
            let lead_zeroes_s1 = s1.iter().take_while(|&&x| x == 0).count();
            let lead_zeroes_d = d.iter().take_while(|&&x| x == 0).count();
            assert!(*d.iter().last().ok_or(crate::Error::InvalidArgument)? != 0);
            assert!(*s1.iter().last().ok_or(crate::Error::InvalidArgument)? != 0);
            assert!(d.len() - lead_zeroes_d >= s1.len() - lead_zeroes_s1);

            let mut n = ExactNumNumber::convert_from_radix(
                s,
                &d,
                e,
                rdx,
                usize::MAX,
                RoundingMode::None,
                &mut cc,
            )?;

            if rdx == Radix::Dec {
                if g.mantissa_max_bit_len() < n.mantissa_max_bit_len() {
                    n.set_precision(g.mantissa_max_bit_len(), RoundingMode::ToEven)?;
                } else {
                    g.set_precision(n.mantissa_max_bit_len(), RoundingMode::ToEven)?;
                }

                // Both legs truncate (`RoundingMode::None`). The exact results are pinned from
                // mpmath (~/rusty-work/logs/pinconv.py, two working precisions agreeing):
                // g = trunc_p(V), d = trunc of g to len(d) digits, n = trunc_p(d), then both
                // rounded to even at the common precision. For the EXPONENT_MIN / 10 - 123
                // vectors n is exactly one ulp below g: d < g, so that is the correct answer.
                let (pg_e, pg_w, pn_e, pn_w) = CONV_DEC_PINS
                    .get(dec_ix)
                    .ok_or(crate::Error::InvalidArgument)?;
                dec_ix += 1;
                assert_eq!(g.exponent(), *pg_e);
                assert!(words_eq_u64(g.mantissa().digits(), pg_w), "g {:?}", g);
                assert_eq!(n.exponent(), *pn_e);
                assert!(words_eq_u64(n.mantissa().digits(), pn_w), "n {:?}", n);
            } else {
                assert!(n.cmp(&g) == 0);
            }
        }

        assert_eq!(dec_ix, CONV_DEC_PINS.len());

        // unknown p: decimal
        Ok(())
    }

    fn random_radix() -> Radix {
        match random::<usize>() % 4 {
            0 => Radix::Bin,
            1 => Radix::Oct,
            2 => Radix::Dec,
            _ => Radix::Hex,
        }
    }

    #[test]
    fn test_round_dec() -> Result<(), crate::Error> {
        let mut testset = [
            (
                RoundingMode::ToEven,
                vec![
                    (
                        ([1, 2, 3, 5, 0, 0], 3, true, true), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 5, 0, 0], false, true),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 5, 0, 0], 3, true, false), // input, n, is_positive, check_roundable
                        ([1, 2, 4, 0, 0, 0], false, false),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 5, 0, 1], 3, true, true), // input, n, is_positive, check_roundable
                        ([1, 2, 4, 0, 0, 0], false, false),  // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 4, 0, 0], 3, true, false), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 0, 0, 0], false, false),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 6, 0, 0], 3, true, false), // input, n, is_positive, check_roundable
                        ([1, 2, 4, 0, 0, 0], false, false),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 4, 9, 9], 3, true, true), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 4, 9, 9], false, true),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 4, 5, 0, 0], 3, true, true), // input, n, is_positive, check_roundable
                        ([1, 2, 4, 5, 0, 0], false, true),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 4, 5, 0, 1], 3, true, true), // input, n, is_positive, check_roundable
                        ([1, 2, 5, 0, 0, 0], false, false),  // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 4, 4, 0, 0], 3, true, false), // input, n, is_positive, check_roundable
                        ([1, 2, 4, 0, 0, 0], false, false),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 4, 6, 0, 0], 3, true, false), // input, n, is_positive, check_roundable
                        ([1, 2, 5, 0, 0, 0], false, false),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 4, 4, 9, 9], 3, true, false), // input, n, is_positive, check_roundable
                        ([1, 2, 4, 0, 0, 0], false, false),   // output, overflow, check_roundable
                    ),
                ],
            ),
            (
                RoundingMode::ToOdd,
                vec![
                    (
                        ([1, 2, 3, 5, 0, 0], 3, true, true), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 5, 0, 0], false, true),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 5, 0, 0], 3, true, false), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 0, 0, 0], false, false),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 5, 0, 1], 3, true, true), // input, n, is_positive, check_roundable
                        ([1, 2, 4, 0, 0, 0], false, false),  // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 4, 0, 0], 3, true, false), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 0, 0, 0], false, false),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 6, 0, 0], 3, true, false), // input, n, is_positive, check_roundable
                        ([1, 2, 4, 0, 0, 0], false, false),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 4, 9, 9], 3, true, false), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 0, 0, 0], false, false),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 4, 5, 0, 0], 3, true, false), // input, n, is_positive, check_roundable
                        ([1, 2, 5, 0, 0, 0], false, false),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 4, 5, 0, 1], 3, true, true), // input, n, is_positive, check_roundable
                        ([1, 2, 5, 0, 0, 0], false, false),  // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 4, 4, 0, 0], 3, true, false), // input, n, is_positive, check_roundable
                        ([1, 2, 4, 0, 0, 0], false, false),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 4, 6, 0, 0], 3, true, false), // input, n, is_positive, check_roundable
                        ([1, 2, 5, 0, 0, 0], false, false),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 4, 4, 9, 9], 3, true, true), // input, n, is_positive, check_roundable
                        ([1, 2, 4, 4, 9, 9], false, true),   // output, overflow, check_roundable
                    ),
                ],
            ),
            (
                RoundingMode::FromZero,
                vec![
                    (
                        ([1, 2, 3, 0, 0, 0], 3, true, true), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 0, 0, 0], false, true),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 0, 0, 0], 3, true, false), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 0, 0, 0], false, false),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 0, 0, 1], 3, true, true), // input, n, is_positive, check_roundable
                        ([1, 2, 4, 0, 0, 0], false, false),  // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 9, 9, 9], 3, true, false), // input, n, is_positive, check_roundable
                        ([1, 2, 4, 0, 0, 0], false, false),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 9, 9, 9], 3, true, true), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 9, 9, 9], false, true),   // output, overflow, check_roundable
                    ),
                ],
            ),
            (
                RoundingMode::ToZero,
                vec![
                    (
                        ([1, 2, 3, 0, 0, 0], 3, true, false), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 0, 0, 0], false, false),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 9, 9, 9], 3, true, true), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 9, 9, 9], false, true),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 9, 9, 9], 3, true, false), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 0, 0, 0], false, false),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 1, 1, 1], 3, true, true), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 0, 0, 0], false, false),  // output, overflow, check_roundable
                    ),
                ],
            ),
            (
                RoundingMode::None,
                vec![
                    (
                        ([1, 2, 3, 1, 2, 3], 3, true, true), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 1, 2, 3], false, false),  // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 1, 2, 3], 3, true, false), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 1, 2, 3], false, false),   // output, overflow, check_roundable
                    ),
                ],
            ),
            (
                RoundingMode::Up,
                vec![
                    (
                        ([1, 2, 3, 0, 0, 0], 3, true, true), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 0, 0, 0], false, true),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 0, 0, 0], 3, true, false), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 0, 0, 0], false, false),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 0, 0, 1], 3, true, true), // input, n, is_positive, check_roundable
                        ([1, 2, 4, 0, 0, 0], false, false),  // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 9, 9, 9], 3, true, false), // input, n, is_positive, check_roundable
                        ([1, 2, 4, 0, 0, 0], false, false),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 0, 0, 0], 3, false, false), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 0, 0, 0], false, false),    // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 9, 9, 9], 3, false, true), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 9, 9, 9], false, true),    // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 9, 9, 9], 3, false, false), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 0, 0, 0], false, false),    // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 1, 1, 1], 3, false, true), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 0, 0, 0], false, false),   // output, overflow, check_roundable
                    ),
                ],
            ),
            (
                RoundingMode::Down,
                vec![
                    (
                        ([1, 2, 3, 0, 0, 0], 3, false, true), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 0, 0, 0], false, true),    // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 0, 0, 0], 3, false, false), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 0, 0, 0], false, false),    // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 0, 0, 1], 3, false, true), // input, n, is_positive, check_roundable
                        ([1, 2, 4, 0, 0, 0], false, false),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 9, 9, 9], 3, false, false), // input, n, is_positive, check_roundable
                        ([1, 2, 4, 0, 0, 0], false, false),    // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 0, 0, 0], 3, true, false), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 0, 0, 0], false, false),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 9, 9, 9], 3, true, true), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 9, 9, 9], false, true),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 9, 9, 9], 3, true, false), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 0, 0, 0], false, false),   // output, overflow, check_roundable
                    ),
                    (
                        ([1, 2, 3, 1, 1, 1], 3, true, true), // input, n, is_positive, check_roundable
                        ([1, 2, 3, 0, 0, 0], false, false),  // output, overflow, check_roundable
                    ),
                ],
            ),
        ];

        for (rm, sets) in testset.iter_mut() {
            for (
                (input, n, is_positive, check_roundable),
                (output, overflow, check_roundable_ret),
            ) in sets.iter_mut()
            {
                match rm {
                    RoundingMode::ToEven
                    | RoundingMode::ToOdd
                    | RoundingMode::FromZero
                    | RoundingMode::ToZero
                    | RoundingMode::None => {
                        // indifferent of sign
                        for is_positive in [true, false] {
                            let ovf = ExactNumNumber::round_dec(
                                input,
                                *n,
                                *rm,
                                is_positive,
                                check_roundable,
                            );
                            assert_eq!(*check_roundable, *check_roundable_ret);
                            assert_eq!(*overflow, ovf);
                            assert_eq!(*input, *output);
                        }
                    }
                    RoundingMode::Down | RoundingMode::Up => {
                        let ovf = ExactNumNumber::round_dec(
                            input,
                            *n,
                            *rm,
                            *is_positive,
                            check_roundable,
                        );
                        assert_eq!(*check_roundable, *check_roundable_ret);
                        assert_eq!(*overflow, ovf);
                        assert_eq!(*input, *output);
                    }
                }
            }
        }

        // overflow
        let mut input = [9, 9, 9, 9, 9, 9];
        let mut check_roundable = false;
        let ovf =
            ExactNumNumber::round_dec(&mut input, 3, RoundingMode::Up, true, &mut check_roundable);
        assert!(ovf);
        assert_eq!(input, [1, 0, 0, 0, 0, 0]);

        // n = input.len()
        let mut input = [9, 9, 9, 9, 9, 9];
        let mut check_roundable = false;
        let ovf =
            ExactNumNumber::round_dec(&mut input, 6, RoundingMode::Up, true, &mut check_roundable);
        assert!(!ovf);
        assert_eq!(input, [9, 9, 9, 9, 9, 9]);

        // n > input.len()
        let mut input = [9, 9, 9, 9, 9, 9];
        let mut check_roundable = false;
        let ovf =
            ExactNumNumber::round_dec(&mut input, 7, RoundingMode::Up, true, &mut check_roundable);
        assert!(!ovf);
        assert_eq!(input, [9, 9, 9, 9, 9, 9]);

        // n = 0
        let mut input = [9, 9, 9, 9, 9, 9];
        let mut check_roundable = false;
        let ovf =
            ExactNumNumber::round_dec(&mut input, 0, RoundingMode::Up, true, &mut check_roundable);
        assert!(!ovf);
        assert_eq!(input, [9, 9, 9, 9, 9, 9]);

        // try round
        assert!(ExactNumNumber::try_round_dec(
            &mut [9, 9, 9, 9, 9, 9],
            3,
            RoundingMode::Up,
            Sign::Pos,
            &mut 0,
            false
        )?);
        assert!(!ExactNumNumber::try_round_dec(
            &mut [9, 9, 9, 9, 9, 9],
            3,
            RoundingMode::Up,
            Sign::Pos,
            &mut 0,
            true
        )?);

        let mut e = 0;
        assert!(ExactNumNumber::try_round_dec(
            &mut [9, 9, 9, 9, 9, 0],
            3,
            RoundingMode::Up,
            Sign::Pos,
            &mut e,
            true
        )?);
        assert_eq!(e, 1);

        e = EXPONENT_MAX;
        assert_eq!(
            ExactNumNumber::try_round_dec(
                &mut [9, 9, 9, 9, 9, 0],
                3,
                RoundingMode::FromZero,
                Sign::Neg,
                &mut e,
                true
            ),
            Err(Error::ExponentOverflow(Sign::Neg))
        );
        Ok(())
    }
}
