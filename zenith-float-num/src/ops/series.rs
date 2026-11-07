// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! Power series computation appliance.

use crate::common::consts::ONE;
use crate::common::util::calc_add_cost;
use crate::common::util::calc_mul_cost;
use crate::common::util::cost_from_u128;
use crate::common::util::cost_product;
use crate::common::util::cost_sum;
use crate::common::util::log2_floor;
use crate::common::util::sqrt_int;
use crate::common::util::wide;
use crate::defs::Error;
use crate::defs::RoundingMode;
use crate::num::ExactNumNumber;
use crate::Sign;

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

const MAX_CACHE: usize = 128;
const RECT_ITER_THRESHOLD: usize = MAX_CACHE / 10 * 9;

//
// Public part
//

/// Generator of polynomial coefficients.
pub(crate) trait PolycoeffGen {
    /// Returns the next polynomial coefficient value.
    fn next(&mut self, rm: RoundingMode) -> Result<&ExactNumNumber, Error>;

    /// Returns the cost of one call to next if numbers have precision p.
    fn iter_cost(&self) -> u64;

    /// Returns true if coefficient is divizor.
    fn is_div(&self) -> bool {
        false
    }
}

/// Shared 1/(k(k+1)) factorial step used by sin, cos, and sinh series.
pub(crate) struct FactPolycoeffGen {
    one_full_p: ExactNumNumber,
    inc: ExactNumNumber,
    fct: ExactNumNumber,
    sign: i8,
    flip_sign: bool,
    iter_cost: u64,
}

impl FactPolycoeffGen {
    pub(crate) fn for_sin(p: usize) -> Result<Self, Error> {
        Self::new(p, true, true)
    }

    pub(crate) fn for_cos(p: usize) -> Result<Self, Error> {
        Self::new(p, false, true)
    }

    pub(crate) fn for_sinh(p: usize) -> Result<Self, Error> {
        Self::new(p, true, false)
    }

    fn new(p: usize, start_from_one: bool, flip_sign: bool) -> Result<Self, Error> {
        let inc = if start_from_one {
            ExactNumNumber::from_word(1, 1)?
        } else {
            ExactNumNumber::new(1)?
        };
        let fct = ExactNumNumber::from_word(1, p)?;
        let one_full_p = ExactNumNumber::from_word(1, p)?;
        let iter_cost = cost_from_u128(cost_product(
            cost_sum(&[
                Some(u128::from(calc_mul_cost(p))),
                Some(u128::from(calc_add_cost(p))),
                Some(u128::from(calc_add_cost(inc.mantissa_max_bit_len()))),
            ]),
            Some(2),
        ));
        Ok(Self {
            one_full_p,
            inc,
            fct,
            sign: 1,
            flip_sign,
            iter_cost,
        })
    }
}

impl PolycoeffGen for FactPolycoeffGen {
    fn next(&mut self, rm: RoundingMode) -> Result<&ExactNumNumber, Error> {
        let p_inc = self.inc.mantissa_max_bit_len();
        let p_one = self.one_full_p.mantissa_max_bit_len();

        self.inc = self.inc.add(&ONE, p_inc, rm)?;
        let inv_inc = self.one_full_p.div(&self.inc, p_one, rm)?;
        self.fct = self.fct.mul(&inv_inc, p_one, rm)?;

        self.inc = self.inc.add(&ONE, p_inc, rm)?;
        let inv_inc = self.one_full_p.div(&self.inc, p_one, rm)?;
        self.fct = self.fct.mul(&inv_inc, p_one, rm)?;

        if self.flip_sign {
            self.sign *= -1;
            if self.sign > 0 {
                self.fct.set_sign(Sign::Pos);
            } else {
                self.fct.set_sign(Sign::Neg);
            }
        }

        Ok(&self.fct)
    }

    #[inline]
    fn iter_cost(&self) -> u64 {
        self.iter_cost
    }
}

/// Estimate how argument reduction influences cost.
pub trait ArgReductionEstimator {
    /// Estimates cost of reduction n times for number with precision p.
    fn reduction_cost(n: usize, p: usize) -> u64;

    /// Given m, the negative power of 2 of a number, returns the negative power of 2 if reduction is applied n times.
    fn reduction_effect(n: usize, m: i64) -> usize;
}

/// Compute the number of reductions required for the best performance.
/// p is the number precision
/// polycoeff_gen is the polynomial coefficient ganerator
/// m is the negative exponent of the number.
/// pwr_step - increment of power of x in each iteration
/// ext - if true use series step cost directly
pub(crate) fn series_cost_optimize<S: ArgReductionEstimator>(
    p: usize,
    polycoeff_gen: &impl PolycoeffGen,
    m: i64,
    pwr_step: usize,
    ext: bool,
) -> (usize, usize, usize) {
    let reduction_num_step = log2_floor(p) / 2;

    if reduction_num_step == 0 {
        let m_eff = S::reduction_effect(0, m);
        let niter = series_niter(p, m_eff) / pwr_step.max(1);
        return (0, niter, m_eff);
    }

    let mut reduction_times = if reduction_num_step as i64 > m {
        (reduction_num_step as i64 - m) as usize
    } else {
        0
    };

    let mut cost1 = u64::MAX;

    loop {
        let m_eff = S::reduction_effect(reduction_times, m);
        let niter = series_niter(p, m_eff) / pwr_step.max(1);
        let series_part = if ext {
            cost_product(
                Some(u128::from(polycoeff_gen.iter_cost())),
                Some(wide(niter)),
            )
        } else {
            Some(u128::from(series_cost(niter, p, polycoeff_gen)))
        };
        let cost2 = cost_from_u128(cost_sum(&[
            series_part,
            Some(u128::from(S::reduction_cost(reduction_times, p))),
        ]));

        if cost2 < cost1 {
            cost1 = cost2;
            reduction_times += reduction_num_step;
            if reduction_times > p.saturating_add(reduction_num_step) {
                return (reduction_times - reduction_num_step, niter, m_eff);
            }
        } else {
            return (reduction_times - reduction_num_step, niter, m_eff);
        }
    }
}

pub(crate) fn series_run<T: PolycoeffGen>(
    acc: ExactNumNumber,
    x_first: ExactNumNumber,
    x_step: ExactNumNumber,
    niter: usize,
    polycoeff_gen: &mut T,
) -> Result<ExactNumNumber, Error> {
    let mut ret = if x_first.is_zero() || x_step.is_zero() {
        series_compute_fast(acc, x_first, polycoeff_gen)
    } else if niter >= RECT_ITER_THRESHOLD {
        series_rectangular(niter, acc, x_first, x_step, polycoeff_gen)
    } else if polycoeff_gen.is_div() {
        series_linear(acc, x_first, x_step, polycoeff_gen)
    } else {
        series_horner(acc, x_first, x_step, polycoeff_gen)
    }?;

    ret.set_inexact(true);

    Ok(ret)
}

//
// Private part
//

/// Estimate of the number of series iterations.
/// p is the precision, m is the negative power of x
/// (i.e. x = f*2^(-m), where 0.5 <= f < 1).
fn series_niter(p: usize, m: usize) -> usize {
    let ln = log2_floor(p);
    let lln = log2_floor(ln);
    p / (ln - lln + m - 2)
}

/// Estimate cost of execution for series.
/// niter is the estimated number of series iterations
/// p is the numbers precision
/// polycoeff_gen is the coefficient generator
fn series_cost<T: PolycoeffGen>(niter: usize, p: usize, polycoeff_gen: &T) -> u64 {
    let cost_mul = u128::from(calc_mul_cost(p));
    let cost_add = u128::from(calc_add_cost(p));
    let iter_cost = u128::from(polycoeff_gen.iter_cost());
    let n = Some(wide(niter));
    // niter * (cost(mul) + cost(add) + cost(polycoeff_gen.next))
    let base = cost_product(
        n,
        cost_sum(&[Some(cost_mul), Some(cost_add), Some(iter_cost)]),
    );

    let total = if niter >= RECT_ITER_THRESHOLD {
        // + isqrt(niter) * cost(mul)
        // + (niter / 10) * (2 * cost(mul) + cost(add) + cost(polycoeff_gen.next))
        let sqrt_term = cost_product(n.map(|v| v.isqrt()), Some(cost_mul));
        let rect_step =
            cost_sum(&[cost_product(Some(2), Some(cost_mul)), Some(cost_add), Some(iter_cost)]);
        let rect_term = cost_product(n.map(|v| v / 10), rect_step);
        cost_sum(&[base, sqrt_term, rect_term])
    } else {
        base
    };
    cost_from_u128(total)
}

// x_first or x_step is zero
fn series_compute_fast<T: PolycoeffGen>(
    acc: ExactNumNumber,
    x_first: ExactNumNumber,
    polycoeff_gen: &mut T,
) -> Result<ExactNumNumber, Error> {
    if x_first.is_zero() {
        Ok(acc)
    } else {
        let p = acc
            .mantissa_max_bit_len()
            .max(x_first.mantissa_max_bit_len());

        let is_div = polycoeff_gen.is_div();

        let coeff = polycoeff_gen.next(RoundingMode::None)?;

        let part = if is_div {
            x_first.div(coeff, p, RoundingMode::None)
        } else {
            x_first.mul(coeff, p, RoundingMode::None)
        }?;

        acc.add(&part, p, RoundingMode::None)
    }
}

// Rectangular series.
// p is the result precision
// niter is the estimated number of iterations
// x_first is the first power of x in the series
// x_step is a multiplication factor for each step
// polycoeff_gen is a generator of polynomial coeeficients for the series
// rm is the rounding mode.
// cost: niter * (O(mul) + O(add) + cost(polcoeff_gen.next)) + sqrt(niter) * O(mul) + cost(series_line(remainder))
fn series_rectangular<T: PolycoeffGen>(
    mut niter: usize,
    add: ExactNumNumber,
    x_first: ExactNumNumber,
    x_step: ExactNumNumber,
    polycoeff_gen: &mut T,
) -> Result<ExactNumNumber, Error> {
    debug_assert!(niter >= 4);

    let p = add
        .mantissa_max_bit_len()
        .max(x_first.mantissa_max_bit_len())
        .max(x_step.mantissa_max_bit_len());

    let mut acc = ExactNumNumber::new(p)?;

    // build cache
    let mut cache = Vec::<ExactNumNumber>::new();
    let sqrt_iter = sqrt_int(niter as u32) as usize;
    let cache_sz = MAX_CACHE.min(sqrt_iter);
    cache.try_reserve_exact(cache_sz)?;
    let mut x_pow = x_step.clone()?;

    for _ in 0..cache_sz {
        // the next power is formed from the current one first; the current one is then moved into the cache
        // (capacity was reserved above, so the push cannot fail).
        let next = x_pow.mul(&x_step, p, RoundingMode::None)?;
        cache.push(x_pow);
        x_pow = next;
    }

    // run computation
    let poly_val = compute_row(p, &cache, polycoeff_gen)?;
    acc = acc.add(&poly_val, p, RoundingMode::None)?;
    // The running terminal power starts equal to `x_pow` (the next power after the cached ones). It is recomputed from
    // the last cached power with the same operands and rounding, which gives the identical value, so that `x_pow` and
    // `terminal_pow` each have one owner.
    let last_cached = cache.last().ok_or(Error::InvalidArgument)?;
    let mut terminal_pow = last_cached.mul(&x_step, p, RoundingMode::None)?;
    niter -= cache_sz;

    loop {
        let poly_val = compute_row(p, &cache, polycoeff_gen)?;
        let part = poly_val.mul(&terminal_pow, p, RoundingMode::None)?;
        acc = acc.add(&part, p, RoundingMode::None)?;
        terminal_pow = terminal_pow.mul(&x_pow, p, RoundingMode::None)?;
        niter -= cache_sz;

        if niter < cache_sz {
            break;
        }
    }
    drop(cache);

    acc = acc.mul(&x_first, p, RoundingMode::None)?;
    terminal_pow = terminal_pow.mul(&x_first, p, RoundingMode::None)?;
    acc = acc.add(&add, p, RoundingMode::None)?;

    acc = if niter < MAX_CACHE * 10 && !polycoeff_gen.is_div() {
        // probably not too many iterations left
        series_horner(acc, terminal_pow, x_step, polycoeff_gen)
    } else {
        series_linear(acc, terminal_pow, x_step, polycoeff_gen)
    }?;

    Ok(acc)
}

// Linear series
// cost: niter * (2 * O(mul) + O(add) + cost(polcoeff_gen.next))
fn series_linear<T: PolycoeffGen>(
    mut acc: ExactNumNumber,
    x_first: ExactNumNumber,
    x_step: ExactNumNumber,
    polycoeff_gen: &mut T,
) -> Result<ExactNumNumber, Error> {
    let p = acc
        .mantissa_max_bit_len()
        .max(x_first.mantissa_max_bit_len())
        .max(x_step.mantissa_max_bit_len());

    let is_div = polycoeff_gen.is_div();
    let mut x_pow = x_first;

    loop {
        let coeff = polycoeff_gen.next(RoundingMode::None)?;
        let part = if is_div {
            x_pow.div(coeff, p, RoundingMode::None)
        } else {
            x_pow.mul(coeff, p, RoundingMode::None)
        }?;

        acc = acc.add(&part, p, RoundingMode::None)?;

        if part.exponent() <= acc.exponent() - acc.mantissa_max_bit_len() as i64 {
            break;
        }

        x_pow = x_pow.mul(&x_step, p, RoundingMode::None)?;
    }

    Ok(acc)
}

// compute row of rectangle series.
fn compute_row<T: PolycoeffGen>(
    p: usize,
    cache: &[ExactNumNumber],
    polycoeff_gen: &mut T,
) -> Result<ExactNumNumber, Error> {
    let is_div = polycoeff_gen.is_div();

    let mut acc = ExactNumNumber::new(p)?;
    let coeff = polycoeff_gen.next(RoundingMode::None)?;
    if is_div {
        let r = coeff.reciprocal(p, RoundingMode::None)?;
        acc = acc.add(&r, p, RoundingMode::None)?;
    } else {
        acc = acc.add(coeff, p, RoundingMode::None)?;
    }

    for x_pow in cache {
        let coeff = polycoeff_gen.next(RoundingMode::None)?;
        let add = if is_div {
            x_pow.div(coeff, p, RoundingMode::None)
        } else {
            x_pow.mul(coeff, p, RoundingMode::None)
        }?;
        acc = acc.add(&add, p, RoundingMode::None)?;
    }

    Ok(acc)
}

/// Horner's method
/// cost: niter*(O(mul) + O(add) + cost(polycoeff_gen.next))
fn series_horner<T: PolycoeffGen>(
    add: ExactNumNumber,
    x_first: ExactNumNumber,
    x_step: ExactNumNumber,
    polycoeff_gen: &mut T,
) -> Result<ExactNumNumber, Error> {
    debug_assert!(x_first.exponent() <= 0);
    debug_assert!(x_step.exponent() <= 0);
    debug_assert!(!polycoeff_gen.is_div());

    let p = add
        .mantissa_max_bit_len()
        .max(x_first.mantissa_max_bit_len())
        .max(x_step.mantissa_max_bit_len());

    // determine number of parts and cache plynomial coeffs.
    let mut cache = Vec::<ExactNumNumber>::new();
    let mut x_p = -x_first.exponent() - x_step.exponent();
    let mut coef_p = 0;

    while x_p + coef_p < p as i64 - add.exponent() {
        let coeff = polycoeff_gen.next(RoundingMode::None)?;
        coef_p = -coeff.exponent();
        x_p += -x_step.exponent();
        cache.push(coeff.clone()?);
    }

    let last_coeff = polycoeff_gen.next(RoundingMode::None)?;
    let mut acc = last_coeff.clone()?;

    for coeff in cache.iter().rev() {
        acc = acc.mul(&x_step, p, RoundingMode::None)?;
        acc = acc.add(coeff, p, RoundingMode::None)?;
    }

    acc = acc.mul(&x_first, p, RoundingMode::None)?;

    add.add(&acc, p, RoundingMode::None)
}

#[cfg(test)]
mod cost_overflow_tests {
    /// Fixed iter cost. `next` is never called by the cost code.
    struct FixedCost(u64);

    impl super::PolycoeffGen for FixedCost {
        fn next(
            &mut self,
            _rm: crate::RoundingMode,
        ) -> Result<&crate::num::ExactNumNumber, crate::defs::Error> {
            Err(crate::defs::Error::InvalidArgument)
        }

        fn iter_cost(&self) -> u64 {
            self.0
        }
    }

    /// Unclamped `f(p)`: `p*p` for `p < 70`, `isqrt(p^3)` for `p < 1625`, `isqrt(p)^3` otherwise.
    fn f(p: u128) -> u128 {
        if p < 70 {
            p * p
        } else if p < 1625 {
            (p * p * p).isqrt()
        } else {
            let q = p.isqrt();
            q * q * q
        }
    }

    /// Series cost for `FixedCost(2 * p)`. `f(p)` stays unclamped; the finished total is clamped.
    fn series_cost_oracle(niter: usize, p: usize) -> u128 {
        let n_u = crate::common::util::wide(niter);
        let p_u = crate::common::util::wide(p);
        let f_p = f(p_u);
        let base = n_u * (f_p + p_u + 2 * p_u);
        let raw = if niter >= super::RECT_ITER_THRESHOLD {
            let sqrt_n = n_u.isqrt();
            base + sqrt_n * f_p + (n_u / 10) * (2 * f_p + p_u + 2 * p_u)
        } else {
            base
        };
        // Clamp the finished cost. f(p) stays unclamped.
        raw.min(u128::from(u64::MAX))
    }

    #[test]
    fn series_cost_matches_integer_oracle() {
        // RECT_ITER_THRESHOLD - 1: MAX_CACHE = 128, threshold = 128/10*9 = 108, so 107 is the last base-only niter.
        // RECT_ITER_THRESHOLD: 108, first niter that adds isqrt(niter)*f(p) + (niter/10)*(2*f(p) + p + 2p).
        // 1_000_000 = 10^6 <= 2^20: rectangular; isqrt = 10^3 and niter/10 = 10^5. f(p) <= 2^90, so the product fits in u128.
        let niters = [super::RECT_ITER_THRESHOLD - 1, super::RECT_ITER_THRESHOLD, 1_000_000_usize];
        // 1625: first p of the isqrt(p)^3 branch; 40^2 = 1_600 <= 1625 < 1_681 = 41^2; <= MAX_PRECISION on every target.
        // 2_643_876 = 1626^2: isqrt = 1626, 1626^3 = 4_298_942_376 > u32::MAX; still <= 2^30 - 2.
        // MAX_PRECISION: 2^30 - 2 on 32-bit (isqrt = 32_767, cube > u32::MAX) and 2^60 - 1 on 64-bit (isqrt = 2^30 - 1, cube > u64::MAX).
        let ps = [1625_usize, 2_643_876_usize, crate::defs::MAX_PRECISION];
        for niter in niters {
            for p in ps {
                // 2p: add cost of two numbers (`calc_add_cost` returns p).
                // 2 * p <= 2^61 always fits u64 on every target.
                let gen = FixedCost(u64::try_from(2 * p).unwrap_or(u64::MAX));
                assert_eq!(
                    u128::from(super::series_cost(niter, p, &gen)),
                    series_cost_oracle(niter, p),
                    "niter {niter} p {p}"
                );
            }
        }
    }

    #[test]
    fn series_cost_exact_below_cap() {
        // RECT_ITER_THRESHOLD = 108. p = 1625 is the first isqrt(p)^3 precision: 40^2 = 1_600 <= 1625 < 41^2, f(p) = 40^3 = 64_000.
        // niter * f(p) and the rectangular terms finish far below u64::MAX, so the oracle clamp does not apply.
        let niter = super::RECT_ITER_THRESHOLD;
        let p = 1625_usize;
        // 2p: add cost of two numbers (`calc_add_cost` returns p).
        // 2 * p <= 2^61 always fits u64 on every target.
        let gen = FixedCost(u64::try_from(2 * p).unwrap_or(u64::MAX));
        assert_eq!(
            u128::from(super::series_cost(niter, p, &gen)),
            series_cost_oracle(niter, p)
        );
    }

    #[cfg(target_pointer_width = "64")]
    #[test]
    fn series_cost_past_u128_is_capped() {
        let p = crate::defs::MAX_PRECISION;
        // 2p: add cost of two numbers (`calc_add_cost` returns p).
        // 2 * p <= 2^61 always fits u64 on every target.
        let gen = FixedCost(u64::try_from(2 * p).unwrap_or(u64::MAX));
        // on 64-bit f(p) = isqrt(2^60 - 1)^3 = (2^30 - 1)^3 ~ 2^90, so niter * f(p) ~ 2^130 does not fit u128
        let niter = 1_usize << 40;
        assert_eq!(super::series_cost(niter, p, &gen), u64::MAX);
    }
}
