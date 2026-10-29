// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! ExactNumNumber formatting.

#[cfg(test)]
use crate::common::util::round_p;
use crate::defs::Error;
use crate::defs::Radix;
use crate::defs::RoundingMode;
#[cfg(test)]
use crate::defs::WORD_BIT_SIZE;
use crate::num::ExactNumNumber;
use crate::Consts;
use crate::Exponent;
use crate::Sign;

#[cfg(feature = "std")]
use std::vec::Vec;

#[cfg(not(feature = "std"))]
use {alloc::string::String, alloc::vec::Vec};

const DIGIT_CHARS: [char; 36] = [
    '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I',
    'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
];

/// Rounding mode for a decimal spelling that parses back to the same value under `rm`.
fn decimal_digit_mode(rm: RoundingMode, positive: bool) -> RoundingMode {
    if rm == RoundingMode::ToEven || rm == RoundingMode::ToOdd {
        return rm;
    }
    let not_above = match rm {
        RoundingMode::Up => true,
        RoundingMode::Down => false,
        RoundingMode::FromZero => positive,
        _ => !positive,
    };
    if not_above == positive {
        RoundingMode::ToZero
    } else {
        RoundingMode::FromZero
    }
}

/// `val` is a `u64` so a 32-bit `usize` cannot truncate or saturate an i64 exponent magnitude.
fn format_uint_radix(val: u64, rdx: Radix, out: &mut String) -> Result<(), Error> {
    let base = u64::from(rdx.value());
    if val == 0 {
        out.push('0');
        return Ok(());
    }
    let mut v = val;
    let mut digits = Vec::new();
    digits.try_reserve_exact(32)?;
    while v > 0 {
        digits.push(u8::try_from(v % base).map_err(|_| Error::InvalidArgument)?);
        v /= base;
    }
    digits.reverse();
    for d in digits {
        let ch = DIGIT_CHARS.get(usize::from(d)).copied().ok_or(Error::InvalidArgument)?;
        out.push(ch);
    }
    Ok(())
}

impl ExactNumNumber {
    /// Parses the number from the string `s` using radix `rdx`, precision `p`, and rounding mode `rm`.
    /// Note, since hexadecimal digits include the character "e", the exponent part is separated
    /// from the mantissa by "_".
    /// For example, a number with mantissa `123abcdef` and exponent `123` would be formatted as `123abcdef_e+123`.
    ///
    /// ## Errors
    ///
    ///  - InvalidArgument: failed to parse input or precision is incorrect.
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - ExponentOverflow: the resulting exponent becomes greater than the maximum allowed value for the exponent.
    ///
    /// A binary numeral below the exponent range is not an error: it rounds once to a subnormal or a signed zero.
    #[cfg(test)]
    pub fn parse(
        s: &str,
        rdx: Radix,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let ps = crate::parser::parse(s, rdx)?;

        if ps.is_nan() || ps.is_inf() {
            Err(Error::InvalidArgument)
        } else if ps.out_of_range() < 0 && rdx == Radix::Bin {
            // 0.digits * 2^e with e below the range: exact digits, placed on the subnormal grid and rounded
            // once by `rm` (a subnormal or a signed zero), as the library parse path does.
            let (m, s, e) = ps.raw_parts();
            let mut v = ExactNumNumber::convert_from_radix(s, m, 0, rdx, usize::MAX, rm, cc)?;
            let w = v.mantissa_max_bit_len().max(round_p(p) + WORD_BIT_SIZE);
            v.set_precision(w, RoundingMode::None)?;
            let at = i128::from(e) + i128::from(v.exponent());
            v.place_scaled(at, false, p, rm)
        } else {
            let (m, s, e) = ps.raw_parts();
            ExactNumNumber::convert_from_radix(s, m, e, rdx, p, rm, cc)
        }
    }

    /// Formats the number using radix `rdx` and rounding mode `rm`.
    /// Note, since hexadecimal digits include the character "e", the exponent part is separated
    /// from the mantissa by "_".
    /// For example, a number with mantissa `123abcdef` and exponent `123` would be formatted as `123abcdef_e+123`.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - ExponentOverflow: the resulting exponent becomes greater than the maximum allowed value for the exponent.
    pub fn format(&self, rdx: Radix, rm: RoundingMode, cc: &mut Consts) -> Result<String, Error> {
        if rdx != Radix::Dec {
            let (s, m, e) = self.convert_to_radix(rdx, rm, cc)?;
            return self.render(rdx, s, &m, e);
        }
        let mut extra = 0usize;
        loop {
            // The first spelling uses `rm`. A directed mode can sit on the side of
            // the value that the same mode parses to a neighbor, so later attempts
            // round the decimal the other way and add digits until it lands back.
            let digit_rm = if extra == 0 { rm } else { decimal_digit_mode(rm, self.is_positive()) };
            let (s, m, e) = self.convert_to_dec_extra(digit_rm, cc, extra)?;
            let text = self.render(rdx, s, &m, e)?;
            if self.decimal_round_trips(&text, rm, cc)? {
                return Ok(text);
            }
            if extra >= 16 {
                return Ok(text);
            }
            extra += 1;
        }
    }

    fn decimal_round_trips(
        &self,
        text: &str,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<bool, Error> {
        let ps = match crate::parser::parse(text, Radix::Dec) {
            Ok(ps) => ps,
            Err(_) => return Ok(false),
        };
        if ps.is_nan() || ps.is_inf() || ps.out_of_range() != 0 {
            return Ok(false);
        }
        let (digits, sign, exp) = ps.raw_parts();
        let back = Self::convert_from_radix(
            sign,
            digits,
            exp,
            Radix::Dec,
            self.mantissa_max_bit_len(),
            rm,
            cc,
        )?;
        Ok(self.cmp(&back) == 0 && self.sign() == back.sign())
    }

    fn render(&self, rdx: Radix, s: Sign, m: &[u8], e: Exponent) -> Result<String, Error> {
        let mut mstr = String::new();
        let mstr_sz = 8
            + (self.mantissa_max_bit_len() + core::mem::size_of::<Exponent>() * 8)
                / rdx.bits_per_digit();

        mstr.try_reserve_exact(mstr_sz)?;

        if s == Sign::Neg {
            mstr.push('-');
        }

        if m.is_empty() {
            mstr.push_str("0.0");
        } else {
            let mut iter = m.iter();

            if self.is_subnormal() {
                mstr.push('0');
            } else {
                let Some(lead) = iter.next() else {
                    return Err(Error::InvalidArgument);
                };
                let Some(ch) = DIGIT_CHARS.get(*lead as usize).copied() else {
                    return Err(Error::InvalidArgument);
                };
                mstr.push(ch);
            }

            mstr.push('.');

            for &d in iter {
                let Some(ch) = DIGIT_CHARS.get(d as usize).copied() else {
                    return Err(Error::InvalidArgument);
                };
                mstr.push(ch);
            }

            if rdx.uses_underscore_exponent() {
                mstr.push('_');
            }

            if e < 1 {
                let val = if self.is_subnormal() {
                    e.unsigned_abs()
                } else {
                    (e - 1).unsigned_abs()
                };

                mstr.push('e');
                mstr.push('-');
                format_uint_radix(val, rdx, &mut mstr)?;
            } else {
                mstr.push('e');
                mstr.push('+');
                format_uint_radix((e - 1).unsigned_abs(), rdx, &mut mstr)?;
            };
        }

        Ok(mstr)
    }
}

#[cfg(test)]
mod tests {

    use crate::common::test_rng::random;

    use crate::{
        common::util::random_subnormal, common::util::test_loop_count, common::util::MAX_DEC_SCALE,
        common::util::TEST_EXP_BOUND, Exponent, WORD_BIT_SIZE,
    };

    use super::*;

    /// A binary numeral whose radix exponent is below `EXPONENT_MIN` is an underflow, not an error: the
    /// test-only `ExactNumNumber::parse` gives the same correctly rounded subnormal or signed zero (sign,
    /// exponent, limbs and inexact flag) as the library parse path `ExactNum::parse`, in every mode.
    #[test]
    fn parse_bin_underflow_matches_library() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let emin = crate::EXPONENT_MIN.unsigned_abs();
        let modes = [
            RoundingMode::ToEven,
            RoundingMode::Up,
            RoundingMode::Down,
            RoundingMode::ToZero,
            RoundingMode::FromZero,
        ];
        for pw in [64u64, 128, 256] {
            let p = usize::try_from(pw).map_err(|_| crate::Error::InvalidArgument)?;
            let below = |k: u64| format!("{:b}", emin + k);
            let cases = [
                format!("1e-{}", below(2)),
                format!("1011011e-{}", below(40)),
                format!("1.0110111e-{}", below(60)),
                format!("1e-{}", below(1 + 64)),
                format!("1e-{}", below(pw + 1)),
                format!("11e-{}", below(pw + 2)),
                format!("1e-{}", below(pw + 10)),
            ];
            for body in &cases {
                for sg in ["", "-"] {
                    let s = format!("{sg}{body}");
                    for rm in modes {
                        let want = crate::ExactNum::parse(&s, Radix::Bin, p, rm, &mut cc);
                        let got = crate::ExactNum::from(ExactNumNumber::parse(
                            &s,
                            Radix::Bin,
                            p,
                            rm,
                            &mut cc,
                        )?);
                        assert!(!want.is_nan(), "library parse of {s} is NaN");
                        assert_eq!(got.sign(), want.sign(), "{s} p={p} {rm:?} sign");
                        assert_eq!(
                            got.mantissa_digits(),
                            want.mantissa_digits(),
                            "{s} p={p} {rm:?} limbs"
                        );
                        if !want.is_zero() {
                            assert_eq!(
                                got.exponent(),
                                want.exponent(),
                                "{s} p={p} {rm:?} exponent"
                            );
                        }
                        assert_eq!(got.inexact(), want.inexact(), "{s} p={p} {rm:?} inexact");
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn test_strop() -> Result<(), crate::Error> {
        let mut eps = ExactNumNumber::from_word(1, 192)?;
        let mut cc = Consts::new()?;

        for i in 0..test_loop_count(1000) {
            let p1 = (random::<usize>() % 32 + 3) * WORD_BIT_SIZE;
            let p2 = (random::<usize>() % 32 + 3) * WORD_BIT_SIZE;
            let p = p1.min(p2);

            // n is built twice from one seed (the rounded reference needs its own owner); the stream then resumes
            // from a fresh seed so later draws do not repeat.
            let build_n = |seed: u64| -> Result<ExactNumNumber, crate::Error> {
                let resume = random::<u64>();
                crate::common::test_rng::reseed_random_this_thread(seed);
                let built = if i & 1 == 0 {
                    ExactNumNumber::random_normal(p1, -TEST_EXP_BOUND, TEST_EXP_BOUND)
                } else {
                    Ok(random_subnormal(p1)?)
                };
                crate::common::test_rng::reseed_random_this_thread(resume);
                built
            };
            let n_seed = random::<u64>();
            let n = build_n(n_seed)?;

            let skip_dec = usize::try_from(n.exponent().unsigned_abs()).map_or(true, |m| m > MAX_DEC_SCALE);

            for rdx in [Radix::Bin, Radix::Oct, Radix::Hex] {
                let rm = match random::<u8>() % 7 {
                    0 => RoundingMode::ToEven,
                    1 => RoundingMode::Up,
                    2 => RoundingMode::Down,
                    3 => RoundingMode::FromZero,
                    4 => RoundingMode::ToZero,
                    5 => RoundingMode::ToOdd,
                    _ => RoundingMode::None,
                };

                let s = n.format(rdx, rm, &mut cc)?;
                let d = ExactNumNumber::parse(&s, rdx, p2, rm, &mut cc)?;

                if p2 < p1 {
                    let mut x = build_n(n_seed)?;
                    x.set_precision(p, rm)?;
                    assert!(x.cmp(&d) == 0);
                } else {
                    assert!(n.cmp(&d) == 0);
                }
            }

            if !skip_dec {
                let s = n.format(Radix::Dec, RoundingMode::ToEven, &mut cc)?;
                let d = ExactNumNumber::parse(&s, Radix::Dec, p2, RoundingMode::ToEven, &mut cc)?;

                if i & 1 == 0 {
                    eps.set_exponent(n.exponent() - p as Exponent + 4);
                    assert!(d.sub(&n, p, RoundingMode::None)?.abs()?.cmp(&eps) < 0);
                } else {
                    let mut eps2 = ExactNumNumber::min_positive(p)?;
                    eps2.set_exponent(eps2.exponent() + 2);
                    assert!(d.sub(&n, p, RoundingMode::None)?.abs()?.cmp(&eps2) < 0);
                }
            }
        }

        // MIN, MAX, min_positive
        let p1 = (random::<usize>() % 32 + 1) * WORD_BIT_SIZE;
        let p2 = (random::<usize>() % 32 + 1) * WORD_BIT_SIZE;
        let p = p1.min(p2);
        let rm = RoundingMode::None;

        for rdx in [Radix::Bin, Radix::Oct, Radix::Dec, Radix::Hex] {
            if rdx == Radix::Dec {
                continue;
            }
            // min, max
            // for p2 < p1 rounding will cause overflow, for p2 >= p1 no rounding is needed.
            for mut n in [ExactNumNumber::max_value(p1)?, ExactNumNumber::min_value(p1)?] {
                //println!("\n{:?} {} {}", rdx, p1, p2);
                //println!("{:?}", n);

                let s = n.format(rdx, rm, &mut cc)?;
                let mut g = ExactNumNumber::parse(&s, rdx, p2, rm, &mut cc)?;

                //println!("{:?}", g);

                if rdx == Radix::Dec {
                    eps.set_exponent(n.exponent() - p as Exponent + 4); // 4, since rm is none
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
            let mut n = ExactNumNumber::min_positive(p1)?;
            let s = n.format(rdx, rm, &mut cc)?;
            let mut g = ExactNumNumber::parse(&s, rdx, p2, rm, &mut cc)?;

            if rdx == Radix::Dec {
                let mut eps = ExactNumNumber::min_positive(p)?;
                eps.set_exponent(eps.exponent() + 2);
                assert!(n.sub(&g, p, rm)?.abs()?.cmp(&eps) < 0);
            } else {
                if p2 < p1 {
                    n.set_precision(p, rm)?;
                } else if p2 > p1 {
                    g.set_precision(p, rm)?;
                }
                assert!(n.cmp(&g) == 0);
            }
        }
        Ok(())
    }
}
