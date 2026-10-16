// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! Exponentiation.

use crate::common::consts::{FOUR, THREE};
use crate::common::util::{bump_prec_retry, calc_add_cost, calc_mul_cost, round_p};
use crate::ops::consts::Consts;
use crate::ops::util::compute_small_exp;
use crate::EXPONENT_MIN;
use crate::{
    common::consts::ONE,
    defs::{Error, WORD_BIT_SIZE, WORD_SIGNIFICANT_BIT},
    num::ExactNumNumber,
    RoundingMode, Sign,
};

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

use super::series::{series_cost_optimize, series_run, ArgReductionEstimator, FactPolycoeffGen};

struct SinhArgReductionEstimator {}

impl ArgReductionEstimator for SinhArgReductionEstimator {
    /// Estimates cost of reduction n times for number with precision p.
    fn reduction_cost(n: usize, p: usize) -> u64 {
        let cost_mul = calc_mul_cost(p);
        let cost_add = calc_add_cost(p);
        let cost_mul2 = calc_mul_cost(THREE.mantissa_max_bit_len());

        n as u64 * (2 * cost_mul + 3 * cost_add + cost_mul2) as u64
    }

    /// Given m, the negative power of 2 of a number, returns the negative power of 2 if reduction is applied n times.
    #[inline]
    fn reduction_effect(n: usize, m: i64) -> usize {
        // n*log2(3) + m
        (n as i64 * 1000 / 631 + m) as usize
    }
}

impl ExactNumNumber {
    /// Computes `e` to the power of `self` with precision `p`. The result is rounded using the rounding mode `rm`.
    /// This function requires constants cache `cc` for computing the result.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - ExponentOverflow: the result is too large or too small number.
    ///  - MemoryAllocation: failed to allocate memory.
    ///  - InvalidArgument: the precision is incorrect.
    pub fn exp(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        let p = round_p(p);

        if self.is_zero() {
            let mut ret = Self::from_word(1, p)?;
            ret.set_inexact(self.inexact());
            return Ok(ret);
        }

        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len());

        compute_small_exp!(
            ONE,
            self.exponent() as i64,
            self.is_negative(),
            p_wrk,
            p,
            rm
        );

        p_wrk += p_inc;

        loop {
            let p_x = p_wrk + 2;

            let mut ret = if self.is_negative() {
                let x = match self.exp_positive_arg(p_x, cc) {
                    Ok(v) => Ok(v),
                    Err(e) => match e {
                        Error::ExponentOverflow(sign) => {
                            return Self::new2(p, sign, self.inexact());
                        }
                        Error::DivisionByZero => Err(Error::DivisionByZero),
                        Error::InvalidArgument => Err(Error::InvalidArgument),
                        Error::PrecisionRetryExhausted => Err(Error::PrecisionRetryExhausted),
                        Error::MemoryAllocation => Err(Error::MemoryAllocation),
                    },
                }?;

                ONE.div(&x, p_x, RoundingMode::FromZero)
            } else {
                self.exp_positive_arg(p_x, cc)
            }?;

            if ret.try_set_precision(p, rm, p_wrk)? {
                ret.set_inexact(ret.inexact() | self.inexact());
                return Ok(ret);
            }

            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    /// Computes `2` to the power of `self` with precision `p`. The result is rounded using `rm`.
    ///
    /// ## Errors
    ///
    ///  - ExponentOverflow: the result is too large or too small.
    ///  - MemoryAllocation: failed to allocate memory.
    ///  - InvalidArgument: the precision is incorrect.
    pub fn exp2(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        ExactNumNumber::small_const(2).pow(self, p, rm, cc)
    }

    /// Computes `10` to the power of `self` with precision `p`. The result is rounded using `rm`.
    ///
    /// ## Errors
    ///
    ///  - ExponentOverflow: the result is too large or too small.
    ///  - MemoryAllocation: failed to allocate memory.
    ///  - InvalidArgument: the precision is incorrect.
    pub fn exp10(&self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Result<Self, Error> {
        ExactNumNumber::small_const(10).pow(self, p, rm, cc)
    }

    // exp for positive argument
    fn exp_positive_arg(&self, p: usize, cc: &mut Consts) -> Result<Self, Error> {
        debug_assert!(!self.is_zero());

        if self.exponent() as i64 > WORD_BIT_SIZE as i64 {
            return Err(Error::ExponentOverflow(self.sign()));
        }

        let p_work = p + 4;

        // compute separately for int and fract parts, then combine the results.
        let int = self.int_as_usize().map_err(|e| {
            if matches!(e, Error::InvalidArgument) {
                Error::ExponentOverflow(Sign::Pos)
            } else {
                e
            }
        })?;
        let e_int = if int > 0 {
            let e_const = cc.e_num(p_work, RoundingMode::None)?;

            e_const.powi(int, p_work, RoundingMode::None)
        } else {
            Ok(ExactNumNumber::small_const(1))
        }?;

        // |fract(self)| = |self| - int, computed from the integer part found above.
        let int_num = Self::from_usize(int)?;
        let mut fract = if self.is_negative() {
            self.add_full_prec(&int_num)
        } else {
            self.sub_full_prec(&int_num)
        }?;
        let e_fract = if !fract.is_zero() {
            fract.set_inexact(false);
            fract.set_precision(p_work, RoundingMode::None)?;
            fract.set_sign(Sign::Pos);
            fract.expf()
        } else {
            Ok(ExactNumNumber::small_const(1))
        }?;

        e_int.mul(&e_fract, p, RoundingMode::FromZero)
    }

    /// Compute the power of `self` to the signed integer `n` with precision `p`. The result is rounded using the rounding mode `rm`.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - ExponentOverflow: the result is too large or too small number.
    ///  - MemoryAllocation: failed to allocate memory.
    ///  - InvalidArgument: the precision is incorrect.
    ///  - DivisionByZero: `self` is zero and `n` is negative.
    pub fn powsi(&self, n: i64, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        if n >= 0 {
            self.powi_internal(n as usize, p, rm, true)
        } else {
            let p = round_p(p);

            if self.is_zero() {
                return Err(Error::DivisionByZero);
            }

            if n == -1 {
                return ONE.div(self, p, rm);
            }

            let mut p_inc = WORD_BIT_SIZE;
            let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc;

            loop {
                let p_x = p_wrk + 2;

                let x = self.powi_internal(
                    crate::common::util::abs_usize(n),
                    p_x,
                    RoundingMode::None,
                    false,
                )?;

                if x.inexact() {
                    let mut ret = ONE.div(&x, p_x, RoundingMode::None)?;

                    if ret.try_set_precision(p, rm, p_wrk)? {
                        return Ok(ret);
                    }
                } else {
                    return ONE.div(&x, p, rm);
                }

                bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
            }
        }
    }

    /// Compute the power of `self` to the integer `n` with precision `p`. The result is rounded using the rounding mode `rm`.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - ExponentOverflow: the result is too large or too small number.
    ///  - MemoryAllocation: failed to allocate memory.
    ///  - InvalidArgument: the precision is incorrect.
    pub fn powi(&self, n: usize, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        self.powi_internal(n, p, rm, true)
    }

    fn powi_internal(
        &self,
        n: usize,
        p: usize,
        rm: RoundingMode,
        update_inexact: bool,
    ) -> Result<Self, Error> {
        let p = round_p(p);

        let mut i = n;

        if i == 0 {
            return Self::from_word(1, p);
        }

        if self.is_zero() || i == 1 {
            let mut ret = self.clone()?;
            ret.set_precision(p, rm)?;
            return Ok(ret);
        }

        let mut bit_pos = WORD_BIT_SIZE;
        while bit_pos > 0 {
            if i & WORD_SIGNIFICANT_BIT as usize != 0 {
                bit_pos -= 1;
                i <<= 1;
                break;
            }

            bit_pos -= 1;
            i <<= 1;
        }
        let bit_pos = bit_pos; // make immutable

        let p = round_p(p);

        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len()) + p_inc;

        loop {
            let p_x = p_wrk + bit_pos * 2;

            let mut ret = || -> Result<Self, Error> {
                let mut x = self.clone()?;
                x.set_inexact(false);

                let mut c;
                let arg = if self.inexact() {
                    c = self.clone()?;
                    c.set_inexact(false);
                    &c
                } else {
                    self
                };

                x.set_precision(p_x, RoundingMode::FromZero)?;
                x = Self::powi_sqr_mul(arg, x, bit_pos, i, p_x)?;
                Ok(x)
            }()
            .map_err(|e| -> Error {
                if let Error::ExponentOverflow(_) = e {
                    Error::ExponentOverflow(if self.is_negative() && (n & 1 == 1) {
                        Sign::Neg
                    } else {
                        Sign::Pos
                    })
                } else {
                    e
                }
            })?;

            if ret.exponent() == EXPONENT_MIN && ret.precision() == 1 {
                if update_inexact {
                    ret.set_inexact(ret.inexact() | self.inexact());
                }
                ret.set_precision(p, rm)?;
                return Ok(ret);
            }

            if ret.try_set_precision(p, rm, p_wrk)? {
                if update_inexact {
                    ret.set_inexact(ret.inexact() | self.inexact());
                }
                return Ok(ret);
            }

            bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
        }
    }

    /// Exact integer power by squaring, then a single rounding to `p`.
    fn pow_exact_integer(
        &self,
        mag: usize,
        negative_exponent: bool,
        p: usize,
        rm: RoundingMode,
    ) -> Result<Self, Error> {
        let p = round_p(p);
        let mut square = self.clone()?;
        square.set_sign(Sign::Pos);
        let mut exp = mag;
        let mut acc = ExactNumNumber::small_const(1);
        while exp > 1 {
            if exp & 1 == 1 {
                acc = acc.mul_full_prec(&square)?;
            }
            square = square.mul_full_prec(&square)?;
            exp >>= 1;
        }
        acc = acc.mul_full_prec(&square)?;
        let odd = mag & 1 == 1;
        let ret_sign = if odd { self.sign() } else { Sign::Pos };
        let effective_rm = if ret_sign.is_negative() {
            match rm {
                RoundingMode::Up => RoundingMode::Down,
                RoundingMode::Down => RoundingMode::Up,
                other => other,
            }
        } else {
            rm
        };
        let mut ret = if negative_exponent {
            ONE.div(&acc, p, effective_rm)?
        } else {
            acc.set_precision(p, effective_rm)?;
            acc
        };
        ret.set_sign(ret_sign);
        Ok(ret)
    }

    fn powi_sqr_mul(
        arg: &Self,
        mut x: Self,
        mut bp: usize,
        mut j: usize,
        p_x: usize,
    ) -> Result<Self, Error> {
        const W: usize = 4;
        let rm = RoundingMode::FromZero;
        if bp < 8 {
            while bp > 0 {
                bp -= 1;
                x = x.mul(&x, p_x, rm)?;
                if j & WORD_SIGNIFICANT_BIT as usize != 0 {
                    x = x.mul(arg, p_x, rm)?;
                }
                j <<= 1;
            }
            return Ok(x);
        }

        let a2 = arg.mul(arg, p_x, rm)?;
        let mut odd: Vec<Self> = Vec::with_capacity(8);
        odd.push(arg.clone()?);
        // arg^3, arg^5, ..., arg^15: each power is formed from the previous one first and then moved into the table
        // (capacity 8 was reserved above, so the pushes cannot fail).
        let mut cur = arg.mul(&a2, p_x, rm)?;
        for _ in 2..8 {
            let next = cur.mul(&a2, p_x, rm)?;
            odd.push(cur);
            cur = next;
        }
        odd.push(cur);

        while bp > 0 {
            if j & WORD_SIGNIFICANT_BIT as usize == 0 {
                x = x.mul(&x, p_x, rm)?;
                bp -= 1;
                j <<= 1;
                continue;
            }
            let take = bp.min(W);
            let mut tmp = j;
            let mut w = 0usize;
            for _ in 0..take {
                w <<= 1;
                if tmp & WORD_SIGNIFICANT_BIT as usize != 0 {
                    w |= 1;
                }
                tmp <<= 1;
            }
            let mut used = take;
            while used > 1 && w & 1 == 0 {
                w >>= 1;
                used -= 1;
            }
            for _ in 0..used {
                x = x.mul(&x, p_x, rm)?;
                bp -= 1;
                j <<= 1;
            }
            x = x.mul(&odd[(w - 1) / 2], p_x, rm)?;
        }
        Ok(x)
    }

    // e^self for |self| < 1.
    fn expf(self) -> Result<Self, Error> {
        debug_assert!(!self.is_zero());

        let p = self.mantissa_max_bit_len();

        let sh = self.sinh_series(p, RoundingMode::None)?; // faster convergence than direct series

        // e = sh + sqrt(sh^2 + 1)
        let sq = sh.mul(&sh, p, RoundingMode::None)?;
        let sq2 = sq.add(&ONE, p, RoundingMode::None)?;
        let sq3 = sq2.sqrt(p, RoundingMode::None)?;
        sq3.add(&sh, p, RoundingMode::FromZero)
    }

    /// sinh using series, for |x| < 1
    pub fn sinh_series(mut self, p: usize, rm: RoundingMode) -> Result<Self, Error> {
        // sinh:  x + x^3/3! + x^5/5! + x^7/7! + ...

        let mut polycoeff_gen = FactPolycoeffGen::for_sinh(p)?;
        let (reduction_times, niter, e_eff) = series_cost_optimize::<SinhArgReductionEstimator>(
            p,
            &polycoeff_gen,
            -(self.exponent() as i64),
            2,
            false,
        );

        // Argument reduction gives error 2^(-p+5) per step, and 2^(-p+3) once.
        // First terms of the series give 2^(-p+5).
        // The error of the remaining terms of the series is compensated.
        let add_prec = reduction_times as i64 * 5 + 8 - e_eff as i64;
        let p_arg = p + if add_prec > 0 { add_prec as usize } else { 0 };
        self.set_precision(p_arg, rm)?;

        let arg = if reduction_times > 0 {
            self.sinh_arg_reduce(reduction_times, rm)?
        } else {
            self
        };

        let x_step = arg.mul(&arg, p_arg, rm)?; // x^2
        let x_first = arg.mul(&x_step, p_arg, rm)?; // x^3
        let acc = arg; // x (moved: the powers above are the last borrows of it)

        let ret = series_run(acc, x_first, x_step, niter, &mut polycoeff_gen)?;

        if reduction_times > 0 {
            ret.sinh_arg_restore(reduction_times, rm)
        } else {
            Ok(ret)
        }
    }

    // reduce argument n times.
    // cost: n * O(add)
    fn sinh_arg_reduce(&self, n: usize, rm: RoundingMode) -> Result<Self, Error> {
        // sinh(3*x) = 3*sinh(x) + 4*sinh(x)^3
        let mut d = ExactNumNumber::small_const(3);
        for _ in 1..n {
            d = d.mul_full_prec(&THREE)?;
        }
        self.div(&d, self.mantissa_max_bit_len(), rm)
    }

    // restore value for the argument reduced n times.
    // cost: n * (4*O(mul) + O(add))
    fn sinh_arg_restore(&self, n: usize, rm: RoundingMode) -> Result<Self, Error> {
        // sinh(3*x) = 3*sinh(x) + 4*sinh(x)^3
        let mut sinh = self.clone()?;
        let p = sinh.mantissa_max_bit_len();

        for _ in 0..n {
            let mut sinh_cub = sinh.mul(&sinh, p, rm)?;
            sinh_cub = sinh_cub.mul(&sinh, p, rm)?;
            let p1 = sinh.mul(&THREE, p, rm)?;
            let p2 = sinh_cub.mul(&FOUR, p, rm)?;
            sinh = p1.add(&p2, p, rm)?;
        }

        Ok(sinh)
    }

    /// Compute the power of `self` to the `n` with precision `p`. The result is rounded using the rounding mode `rm`.
    /// This function requires constants cache `cc` for computing the result.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - ExponentOverflow: the result is too large or too small number.
    ///  - MemoryAllocation: failed to allocate memory.
    ///  - InvalidArgument: `self` is negative, and `n` is not an integer number; the precision is incorrect.
    pub fn pow(
        &self,
        n: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        let nint = n.is_int();
        if n.is_zero() {
            let mut ret = Self::from_word(1, p)?;
            if n.inexact() {
                if self.exponent() == 1 && self.abs_cmp(&ONE) == 0 {
                    ret.set_inexact(self.inexact());
                } else {
                    ret.set_inexact(true);
                }
            }
            return Ok(ret);
        } else if self.is_zero() {
            return if n.is_negative() {
                if n.is_odd_int() {
                    Err(Error::ExponentOverflow(self.sign()))
                } else {
                    Err(Error::ExponentOverflow(Sign::Pos))
                }
            } else {
                let s = if n.is_odd_int() { self.sign() } else { Sign::Pos };
                Self::new2(p, s, self.inexact())
            };
        } else if self.exponent() == 1 && self.abs_cmp(&ONE) == 0 {
            // 1^n or (-1)^n
            if n.is_odd_int() {
                let mut ret = Self::from_i8(self.sign().to_int(), p)?;
                ret.set_inexact(self.inexact() || n.inexact());
                return Ok(ret);
            } else {
                if self.is_positive() || n.is_int() {
                    let mut ret = Self::from_word(1, p)?;
                    ret.set_inexact(self.inexact() || n.inexact());
                    return Ok(ret);
                } else {
                    return Err(Error::InvalidArgument);
                }
            };
        } else if self.is_negative() && !nint {
            return Err(Error::InvalidArgument);
        } else if n.exponent() == 1 && n.abs_cmp(&ONE) == 0 {
            if n.is_positive() {
                let mut ret = self.clone()?;
                ret.set_precision(p, rm)?;
                ret.set_inexact(ret.inexact() || n.inexact());
                return Ok(ret);
            } else {
                let mut ret = ONE.div(self, p, rm)?;
                ret.set_inexact(ret.inexact() || n.inexact());
                return Ok(ret);
            }
        }

        // An exact integer power is a product, then one rounding. The logarithm
        // path cannot certify an exact dyadic under a directed mode.
        if n.is_int() && !n.inexact() && !self.inexact() {
            if let Ok(mag) = n.int_as_usize() {
                let width = self.mantissa_max_bit_len().max(1);
                if mag >= 2 && width.saturating_mul(mag) <= 1_048_576 {
                    return self.pow_exact_integer(mag, n.is_negative(), p, rm);
                }
            }
        }

        // A half-integer power of a value with an exact square root is exact or one
        // rounding of an exact product: 9^0.5 = 3, 4^(-1.5) = 1/8.
        if rm != RoundingMode::None
            && !n.inexact()
            && !self.inexact()
            && self.is_positive()
            && !n.is_int()
            && n.exponent() < 64
        {
            let mut twice = n.clone()?;
            twice.set_exponent(n.exponent() + 1);
            if twice.is_int() {
                let w = self.mantissa_max_bit_len() + WORD_BIT_SIZE;
                let root = self.sqrt(w, RoundingMode::ToZero)?;
                if !root.inexact() {
                    if let Ok(mag) = twice.int_as_usize() {
                        if mag == 1 {
                            return if n.is_negative() {
                                ONE.div(&root, p, rm)
                            } else {
                                let mut ret = root;
                                ret.set_precision(p, rm)?;
                                Ok(ret)
                            };
                        }
                        let width = root.mantissa_max_bit_len().max(1);
                        if width.saturating_mul(mag) <= 1_048_576 {
                            return root.pow_exact_integer(mag, n.is_negative(), p, rm);
                        }
                    }
                }
            }
        }

        // actual pow : self^n = e^(n * ln(self))

        let p = round_p(p);

        let mut x = self.clone()?;
        let mut y = n.clone()?;

        x.set_sign(Sign::Pos);
        let ret_sign = if nint && n.is_odd_int() { self.sign() } else { Sign::Pos };
        let effective_rm = if ret_sign.is_negative() {
            match rm {
                RoundingMode::Up => RoundingMode::Down,
                RoundingMode::Down => RoundingMode::Up,
                v => v,
            }
        } else {
            rm
        };

        x.set_inexact(false);
        y.set_inexact(false);

        let mut p_inc = WORD_BIT_SIZE;
        let mut p_wrk = p.max(self.mantissa_max_bit_len().max(y.mantissa_max_bit_len())) + p_inc;

        match || -> Result<Self, Error> {
            loop {
                let p_ext = p_wrk + 2;
                x.set_precision(p_ext, RoundingMode::None)?;

                let ln = x.ln(p_ext, RoundingMode::None, cc)?;

                let mut m = match y.mul(&ln, p_ext, RoundingMode::None) {
                    Ok(v) => Ok(v),
                    Err(e) => match e {
                        Error::ExponentOverflow(Sign::Neg) => {
                            return Self::new2(p, Sign::Pos, n.inexact() || ln.inexact());
                        }
                        Error::ExponentOverflow(Sign::Pos) => {
                            Err(Error::ExponentOverflow(Sign::Pos))
                        }
                        Error::DivisionByZero => Err(Error::DivisionByZero),
                        Error::InvalidArgument => Err(Error::InvalidArgument),
                        Error::PrecisionRetryExhausted => Err(Error::PrecisionRetryExhausted),
                        Error::MemoryAllocation => Err(Error::MemoryAllocation),
                    },
                }?;

                compute_small_exp!(
                    ONE,
                    m.exponent() as i64,
                    m.is_negative(),
                    p_ext,
                    p,
                    effective_rm
                );

                let mut ret = if m.is_negative() {
                    m.set_sign(Sign::Pos);

                    let u = match m.exp(p_ext, RoundingMode::None, cc) {
                        Ok(v) => Ok(v),
                        Err(e) => match e {
                            Error::ExponentOverflow(sign) => {
                                return Self::new2(p, sign, self.inexact());
                            }
                            Error::DivisionByZero => Err(Error::DivisionByZero),
                            Error::InvalidArgument => Err(Error::InvalidArgument),
                            Error::PrecisionRetryExhausted => Err(Error::PrecisionRetryExhausted),
                            Error::MemoryAllocation => Err(Error::MemoryAllocation),
                        },
                    }?;

                    ONE.div(&u, p_ext, RoundingMode::FromZero)
                } else {
                    m.exp(p_ext, RoundingMode::None, cc)
                }?;

                if ret.try_set_precision(p, effective_rm, p_wrk)? {
                    ret.set_inexact(ret.inexact() | self.inexact() | n.inexact());
                    return Ok(ret);
                }

                bump_prec_retry(&mut p_wrk, &mut p_inc, p)?;
            }
        }() {
            Ok(mut v) => {
                v.set_sign(ret_sign);
                Ok(v)
            }
            Err(e) => match e {
                Error::ExponentOverflow(_) => Err(Error::ExponentOverflow(ret_sign)),
                Error::DivisionByZero => Err(Error::DivisionByZero),
                Error::InvalidArgument => Err(Error::InvalidArgument),
                Error::PrecisionRetryExhausted => Err(Error::PrecisionRetryExhausted),
                Error::MemoryAllocation => Err(Error::MemoryAllocation),
            },
        }
    }
}

#[cfg(test)]
mod test {

    use crate::common::{consts::TWO, util::random_subnormal};

    use super::*;

    #[test]
    fn test_power() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;

        // near 1
        let p = 320;
        let d1 = ExactNumNumber::parse(
            "F.FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF2DC85F7E77EC4872DC85F7E77EC487_e-1",
            crate::Radix::Hex,
            p,
            RoundingMode::None,
            &mut cc,
        )?;
        let d2 = d1.exp(p, RoundingMode::ToEven, &mut cc)?;
        let d3 = ExactNumNumber::parse(
            "2.B7E151628AED2A6ABF7158809CF4F3C762E7160F38B4DA56A5496AF9D95160A40F47A2ECF1C6AEA0",
            crate::Radix::Hex,
            p,
            RoundingMode::None,
            &mut cc,
        )?;

        // println!("{:?}", d1.format(crate::Radix::Bin, RoundingMode::None).unwrap());
        // println!("{:?}", d2.format(crate::Radix::Hex, RoundingMode::None).unwrap());

        assert!(d2.cmp(&d3) == 0);

        let d1 = ExactNumNumber::parse(
            "1.00000000000000000000000000000000000000000000000000000000000000002DC85F7E77EC487C",
            crate::Radix::Hex,
            p,
            RoundingMode::None,
            &mut cc,
        )?;
        let d2 = d1.exp(p, RoundingMode::ToEven, &mut cc)?;
        let d3 = ExactNumNumber::parse(
            "2.B7E151628AED2A6ABF7158809CF4F3C762E7160F38B4DA56A784D9045190CFEFAEC1BE22DDEADB48",
            crate::Radix::Hex,
            p,
            RoundingMode::None,
            &mut cc,
        )?;

        // println!("{:?}", d1.format(crate::Radix::Bin, RoundingMode::None).unwrap());
        // println!("{:?}", d2.format(crate::Radix::Hex, RoundingMode::None).unwrap());

        assert!(d2.cmp(&d3) == 0);

        // random subnormal
        let d1 = random_subnormal(p)?;
        assert!(d1.exp(p, RoundingMode::ToEven, &mut cc)?.cmp(&ONE) == 0);

        // pow(small, small)
        let p = 256;
        let mut d1 = ExactNumNumber::parse("1.10010100010011110010001011101010010101111000010101010100010010111001000100011101010111011010110011011111110000010010101000110001001000000101000010101111001100110111100011101000110001001000000101000010101111001100110111100011101010011101101", crate::Radix::Bin, p, RoundingMode::None, &mut cc)?;
        d1.set_exponent(-123456);
        let mut d2 = ExactNumNumber::parse("1.00110001001000010100001010111100110011011110001110101001110110100000101000010010101111000010110010100010011110010001011101010101010001001011100110111111100100010001110101011101101011000001001010100101111001100110111100011101000110001001000", crate::Radix::Bin, p, RoundingMode::None, &mut cc)?;
        d2.set_exponent(-123);
        let d3 = d1.pow(&d2, p, RoundingMode::ToEven, &mut cc)?;
        let d4 = ExactNumNumber::parse("1.111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111100111000110010011110111100111111100111010011001110010010110010100100110010110111110110100011010000110110011010111101110000011100011010100110100000001e-1", crate::Radix::Bin, p, RoundingMode::None, &mut cc)?;

        assert!(d4.cmp(&d3) == 0);

        // pow(large, small)
        let mut d1 = ExactNumNumber::parse("1.10010100010011110010001011101010010101111000010101010100010010111001000100011101010111011010110011011111110000010010101000110001001000000101000010101111001100110111100011101000110001001000000101000010101111001100110111100011101010011101101", crate::Radix::Bin, p, RoundingMode::None, &mut cc)?;
        d1.set_exponent(123456);
        let mut d2 = ExactNumNumber::parse("1.00110001001000010100001010111100110011011110001110101001110110100000101000010010101111000010110010100010011110010001011101010101010001001011100110111111100100010001110101011101101011000001001010100101111001100110111100011101000110001001000", crate::Radix::Bin, p, RoundingMode::None, &mut cc)?;
        d2.set_exponent(-123);
        let d3 = d1.pow(&d2, p, RoundingMode::ToEven, &mut cc)?;
        let d4 = ExactNumNumber::parse("1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000110001110011010111001000101100100101011011111001101101011011010100101111001001000010000010100110110011000001011011111000000000010111010000010110111", crate::Radix::Bin, p, RoundingMode::None, &mut cc)?;

        assert!(d4.cmp(&d3) == 0);

        // pow(>~1, large)
        let d1 = ExactNumNumber::parse("1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000110001110011010111001000101100100101011011111001101101011011010100101111001001000010000010100110110011000001011011111000000000010111010000010110111", crate::Radix::Bin, p, RoundingMode::None, &mut cc)?;
        let mut d2 = ExactNumNumber::parse("1.00110001001000010100001010111100110011011110001110101001110110100000101000010010101111000010110010100010011110010001011101010101010001001011100110111111100100010001110101011101101011000001001010100101111001100110111100011101000110001001000", crate::Radix::Bin, p, RoundingMode::None, &mut cc)?;
        d2.set_exponent(123);
        let d3 = d1.pow(&d2, p, RoundingMode::ToEven, &mut cc)?;
        let d4 = ExactNumNumber::parse("1.000011101011111010010100111000101101011001100110110100101011001110100110011011100100100101011100010100011101011111010011100011010000110001100111001011000001111110001110000110000000101101011001010001000011011001000010001110000100110000110011111100000001e+1010101101000111", crate::Radix::Bin, p, RoundingMode::None, &mut cc)?;

        assert!(d4.cmp(&d3) == 0);

        // pow(<~1, large)
        let d1 = ExactNumNumber::parse("0.1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111110001110011010111001000101100100101011011111001101101011011010100101111001001000010000010100110110011000001011011111000000000010111010000010110111", crate::Radix::Bin, p, RoundingMode::None, &mut cc)?;
        let mut d2 = ExactNumNumber::parse("1.00110001001000010100001010111100110011011110001110101001110110100000101000010010101111000010110010100010011110010001011101010101010001001011100110111111100100010001110101011101101011000001001010100101111001100110111100011101000110001001000", crate::Radix::Bin, p, RoundingMode::None, &mut cc)?;
        d2.set_exponent(123);
        let d3 = d1.pow(&d2, p, RoundingMode::ToEven, &mut cc)?;
        let d4 = ExactNumNumber::parse("1.010011011001100100100000101011101001110011001100011001011011010010101111101100111101011010101101101111100111001001001101000111110001010010111101110010011100001011011111001111001001001100100111010100100010010111011011110010110101110000100111001100100000011e-110000110100111100", crate::Radix::Bin, p, RoundingMode::None, &mut cc)?;

        // println!("{:?}", d1.format(crate::Radix::Bin, RoundingMode::None).unwrap());
        // println!("{:?}", d2.format(crate::Radix::Bin, RoundingMode::None).unwrap());
        // println!("{:?}", d3.format(crate::Radix::Bin, RoundingMode::None).unwrap());

        assert!(d4.cmp(&d3) == 0);

        // max, min, subnormal
        let d1 = ExactNumNumber::max_value(p)?;
        let d2 = ExactNumNumber::min_value(p)?;
        let d3 = ExactNumNumber::min_positive(p)?;

        assert!(matches!(
            d1.exp(p, RoundingMode::ToEven, &mut cc),
            Err(Error::ExponentOverflow(Sign::Pos))
        ));
        assert!(matches!(
            d1.pow(&d1, p, RoundingMode::ToEven, &mut cc),
            Err(Error::ExponentOverflow(Sign::Pos))
        ));
        assert!(d1.pow(&d2, p, RoundingMode::ToEven, &mut cc)?.is_zero());
        assert!(d1.pow(&d3, p, RoundingMode::ToEven, &mut cc)?.cmp(&ONE) == 0);

        assert!(d2.exp(p, RoundingMode::ToEven, &mut cc)?.is_zero());

        assert!(d3.exp(p, RoundingMode::ToEven, &mut cc)?.cmp(&ONE) == 0);

        assert!(d3.pow(&d1, p, RoundingMode::ToEven, &mut cc)?.is_zero());
        assert!(matches!(
            d3.pow(&d2, p, RoundingMode::ToEven, &mut cc),
            Err(Error::ExponentOverflow(Sign::Pos))
        ));
        assert!(d3.pow(&d3, p, RoundingMode::ToEven, &mut cc)?.cmp(&ONE) == 0);

        assert!(d1.pow(&ONE, p, RoundingMode::ToEven, &mut cc)?.cmp(&d1) == 0);
        assert!(d3.pow(&ONE, p, RoundingMode::ToEven, &mut cc)?.cmp(&d3) == 0);
        assert!(ONE.pow(&d1, p, RoundingMode::ToEven, &mut cc)?.cmp(&ONE) == 0);
        assert!(ONE.pow(&d2, p, RoundingMode::ToEven, &mut cc)?.cmp(&ONE) == 0);
        assert!(ONE.pow(&d3, p, RoundingMode::ToEven, &mut cc)?.cmp(&ONE) == 0);

        let rm = RoundingMode::ToOdd;
        let d1 = TWO.pow(&ONE.neg()?, p, rm, &mut cc)?;
        let d2 = ONE.div(&TWO, p, rm)?;

        assert!(d1.cmp(&d2) == 0);
        Ok(())
    }

    #[test]
    fn test_exp2_exp10() -> Result<(), crate::Error> {
        let mut cc = Consts::new()?;
        let p = 320;
        let rm = RoundingMode::ToEven;

        let x = ExactNumNumber::from_word(3, p)?;
        let d2 = x.exp2(p, rm, &mut cc)?;
        let d2_ref = TWO.pow(&x, p, rm, &mut cc)?;
        assert_eq!(d2.cmp(&d2_ref), 0);

        let d10 = x.exp10(p, rm, &mut cc)?;
        let d10_ref = ExactNumNumber::small_const(10).pow(&x, p, rm, &mut cc)?;
        assert_eq!(d10.cmp(&d10_ref), 0);
        Ok(())
    }
}
