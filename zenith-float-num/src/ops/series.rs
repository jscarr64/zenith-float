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
///
/// The returned `(reduction_times, niter, m_eff)` always belong to the same evaluated candidate.
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

    // The first candidate is selected before any comparison. `u64::MAX` is a real capped
    // cost, so it cannot stand in for "no candidate yet".
    let (mut niter, mut m_eff, mut cost) =
        series_candidate_cost::<S>(p, polycoeff_gen, m, pwr_step, ext, reduction_times);

    loop {
        // A sum that does not fit in `usize` is past the limit checked below.
        let Some(next_times) = reduction_times.checked_add(reduction_num_step) else {
            return (reduction_times, niter, m_eff);
        };
        if next_times > p.saturating_add(reduction_num_step) {
            return (reduction_times, niter, m_eff);
        }
        let (next_niter, next_m_eff, next_cost) =
            series_candidate_cost::<S>(p, polycoeff_gen, m, pwr_step, ext, next_times);
        if next_cost < cost {
            reduction_times = next_times;
            niter = next_niter;
            m_eff = next_m_eff;
            cost = next_cost;
        } else {
            return (reduction_times, niter, m_eff);
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

/// `(niter, m_eff, cost)` for one reduction count in `series_cost_optimize`.
fn series_candidate_cost<S: ArgReductionEstimator>(
    p: usize,
    polycoeff_gen: &impl PolycoeffGen,
    m: i64,
    pwr_step: usize,
    ext: bool,
    reduction_times: usize,
) -> (usize, usize, u64) {
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
    let cost = cost_from_u128(cost_sum(&[
        series_part,
        Some(u128::from(S::reduction_cost(reduction_times, p))),
    ]));
    (niter, m_eff, cost)
}

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

    use super::ArgReductionEstimator;
    use super::PolycoeffGen;

    type Triple = (usize, usize, usize);

    /// Floor log2, copied so the model does not call `log2_floor`.
    fn log2_floor_ref(mut n: usize) -> usize {
        let mut ret = 0usize;
        while n > 1 {
            ret += 1;
            n >>= 1;
        }
        ret
    }

    /// Schoolbook / Toom-3 multiply cost, already clamped to `u64` and widened.
    /// `p < 70`: `p^2`. `p < 1625`: `isqrt(p^3)` (`p^3 < 2^32`). Otherwise `isqrt(p)^3`
    /// (`<= 2^96` for a `usize` precision).
    fn mul_clamped(p: usize) -> u128 {
        let p_u = u128::try_from(p).unwrap_or(u128::MAX);
        let raw = if p_u < 70 {
            p_u * p_u
        } else if p_u < 1625 {
            (p_u * p_u * p_u).isqrt()
        } else {
            let q = p_u.isqrt();
            q * q * q
        };
        raw.min(u128::from(u64::MAX))
    }

    fn add_clamped(p: usize) -> u128 {
        u128::from(u64::try_from(p).unwrap_or(u64::MAX))
    }

    fn pwr_div(pwr_step: usize) -> usize {
        if pwr_step == 0 {
            1
        } else {
            pwr_step
        }
    }

    /// `p / (floor(log2(p)) - floor(log2(floor(log2(p)))) + m - 2)`.
    /// A zero or overflowing denominator is outside the cases these tests build;
    /// the fallback is not a series length the production division would return.
    fn series_niter_ref(p: usize, m_eff: usize) -> usize {
        let ln = log2_floor_ref(p);
        let lln = log2_floor_ref(ln);
        match ln
            .checked_sub(lln)
            .and_then(|v| v.checked_add(m_eff))
            .and_then(|v| v.checked_sub(2))
        {
            Some(d) if d > 0 => p / d,
            Some(_) | None => usize::MAX,
        }
    }

    /// Independent `series_cost`: clamped mul and add, then the base or rectangular total, clamped once.
    fn series_cost_ref(niter: usize, p: usize, iter_cost: u64) -> u64 {
        let cost_mul = mul_clamped(p);
        let cost_add = add_clamped(p);
        let iter = u128::from(iter_cost);
        let n_u = u128::try_from(niter).unwrap_or(u128::MAX);
        let per = cost_mul
            .checked_add(cost_add)
            .and_then(|s| s.checked_add(iter));
        let base = per.and_then(|v| n_u.checked_mul(v));
        let total = if niter >= super::RECT_ITER_THRESHOLD {
            let sqrt_term = n_u.isqrt().checked_mul(cost_mul);
            let rect_step = cost_mul
                .checked_mul(2)
                .and_then(|s| s.checked_add(cost_add))
                .and_then(|s| s.checked_add(iter));
            let rect_term = rect_step.and_then(|step| (n_u / 10).checked_mul(step));
            match (base, sqrt_term, rect_term) {
                (Some(b), Some(s), Some(r)) => b.checked_add(s).and_then(|v| v.checked_add(r)),
                (None, _, _) | (_, None, _) | (_, _, None) => None,
            }
        } else {
            base
        };
        match total {
            Some(v) => u64::try_from(v).unwrap_or(u64::MAX),
            None => u64::MAX,
        }
    }

    struct Cand {
        p: usize,
        iter_cost: u64,
        m: i64,
        div: usize,
        ext: bool,
    }

    /// One candidate: `(niter, m_eff, capped cost)`. Does not call `series_cost` or `series_candidate_cost`.
    fn candidate_ref<S: ArgReductionEstimator>(
        inp: &Cand,
        reduction_times: usize,
    ) -> (usize, usize, u64) {
        let m_eff = S::reduction_effect(reduction_times, inp.m);
        let niter = series_niter_ref(inp.p, m_eff) / inp.div;
        let series_part = if inp.ext {
            u128::from(inp.iter_cost).checked_mul(u128::try_from(niter).unwrap_or(u128::MAX))
        } else {
            Some(u128::from(series_cost_ref(niter, inp.p, inp.iter_cost)))
        };
        let red = u128::from(S::reduction_cost(reduction_times, inp.p));
        let cost = match series_part.and_then(|s| s.checked_add(red)) {
            Some(v) => u64::try_from(v).unwrap_or(u64::MAX),
            None => u64::MAX,
        };
        (niter, m_eff, cost)
    }

    fn initial_times(step: usize, m: i64) -> usize {
        let step_i = i64::try_from(step).unwrap_or(i64::MAX);
        if step_i > m {
            let diff = step_i.checked_sub(m).unwrap_or(i64::MAX);
            usize::try_from(diff).unwrap_or(usize::MAX)
        } else {
            0
        }
    }

    /// Reference search. First candidate is held. A later one replaces it only when its
    /// capped cost is strictly lower. The next count is not considered once it passes `p + step`.
    fn model_optimize<S: ArgReductionEstimator>(
        p: usize,
        iter_cost: u64,
        m: i64,
        pwr_step: usize,
        ext: bool,
    ) -> Triple {
        let inp = Cand {
            p,
            iter_cost,
            m,
            div: pwr_div(pwr_step),
            ext,
        };
        let step = log2_floor_ref(p) / 2;
        if step == 0 {
            let m_eff = S::reduction_effect(0, m);
            let niter = series_niter_ref(p, m_eff) / inp.div;
            return (0, niter, m_eff);
        }
        let mut reduction_times = initial_times(step, m);
        let (mut niter, mut m_eff, mut cost) = candidate_ref::<S>(&inp, reduction_times);
        let mut rounds = 0usize;
        loop {
            rounds += 1;
            assert!(
                rounds <= p.checked_add(2).unwrap_or(usize::MAX),
                "model_optimize did not stop p={p}"
            );
            let Some(next_times) = reduction_times.checked_add(step) else {
                return (reduction_times, niter, m_eff);
            };
            let limit = p.checked_add(step).unwrap_or(usize::MAX);
            if next_times > limit {
                return (reduction_times, niter, m_eff);
            }
            let (next_niter, next_m_eff, next_cost) = candidate_ref::<S>(&inp, next_times);
            if next_cost < cost {
                reduction_times = next_times;
                niter = next_niter;
                m_eff = next_m_eff;
                cost = next_cost;
            } else {
                return (reduction_times, niter, m_eff);
            }
        }
    }

    fn optimize_pair<S: ArgReductionEstimator>(
        p: usize,
        iter_cost: u64,
        m: i64,
        pwr_step: usize,
        ext: bool,
    ) -> (Triple, Triple) {
        let gen = FixedCost(iter_cost);
        let got = super::series_cost_optimize::<S>(p, &gen, m, pwr_step, ext);
        let exp = model_optimize::<S>(p, iter_cost, m, pwr_step, ext);
        (got, exp)
    }

    /// Reduction cost grows by at least `2^32` per step. For `p <= 2048` and an iter cost of
    /// 1 or 1000 the whole series term is below `2^32`, so the first candidate is always kept.
    /// `reduction_effect` is `max(3, n + m)`. The floor keeps `series_niter`'s denominator at least 1
    /// for every `p >= 1` (`floor(log2(p)) - floor(log2(floor(log2(p))))` is 0 at `p = 1`).
    struct StepCost;

    impl ArgReductionEstimator for StepCost {
        fn reduction_cost(n: usize, _p: usize) -> u64 {
            let n_u = u64::try_from(n).unwrap_or(u64::MAX);
            n_u.checked_mul(1u64 << 32).unwrap_or(u64::MAX)
        }

        fn reduction_effect(n: usize, m: i64) -> usize {
            let n_i = i64::try_from(n).unwrap_or(i64::MAX);
            let sum = n_i.checked_add(m).unwrap_or(i64::MAX);
            if sum < 3 {
                3
            } else {
                usize::try_from(sum).unwrap_or(usize::MAX)
            }
        }
    }

    /// `m_eff = 64 + n`. At `p = 64` the series denominator is `m_eff + 2 > 64`, so `niter` is 0
    /// and the capped cost is exactly `reduction_cost`. The offset still changes with `n`, so a
    /// swapped `(count, m_eff)` pair fails.
    fn offset_effect(n: usize) -> usize {
        n.checked_add(64).unwrap_or(usize::MAX)
    }

    /// Costs 100, 50, 10 at counts 3, 6, 9, then 40. Two strict improvements, then a rise.
    struct DropTwoThenRise;

    impl ArgReductionEstimator for DropTwoThenRise {
        fn reduction_cost(n: usize, _p: usize) -> u64 {
            if n == 3 {
                100
            } else if n == 6 {
                50
            } else if n == 9 {
                10
            } else {
                40
            }
        }

        fn reduction_effect(n: usize, _m: i64) -> usize {
            offset_effect(n)
        }
    }

    /// `1_000_000_000 - n`, strictly decreasing for every count this search can reach.
    struct DecreaseForever;

    impl ArgReductionEstimator for DecreaseForever {
        fn reduction_cost(n: usize, _p: usize) -> u64 {
            let n_u = u64::try_from(n).unwrap_or(u64::MAX);
            1_000_000_000u64.checked_sub(n_u).unwrap_or(0)
        }

        fn reduction_effect(n: usize, _m: i64) -> usize {
            offset_effect(n)
        }
    }

    /// Iter cost `u64::MAX` with `niter >= 2` caps every candidate. `m_eff = 8 + n`.
    struct AllCapped;

    impl ArgReductionEstimator for AllCapped {
        fn reduction_cost(_n: usize, _p: usize) -> u64 {
            0
        }

        fn reduction_effect(n: usize, _m: i64) -> usize {
            n.checked_add(8).unwrap_or(usize::MAX)
        }
    }

    /// First count costs `u64::MAX`, the next is 50, the one after is 80.
    struct CapThenLower;

    impl ArgReductionEstimator for CapThenLower {
        fn reduction_cost(n: usize, _p: usize) -> u64 {
            if n == 3 {
                u64::MAX
            } else if n == 6 {
                50
            } else {
                80
            }
        }

        fn reduction_effect(n: usize, _m: i64) -> usize {
            offset_effect(n)
        }
    }

    /// 500, then 100, then `u64::MAX` for every later count.
    struct ImproveThenCap;

    impl ArgReductionEstimator for ImproveThenCap {
        fn reduction_cost(n: usize, _p: usize) -> u64 {
            if n == 3 {
                500
            } else if n == 6 {
                100
            } else {
                u64::MAX
            }
        }

        fn reduction_effect(n: usize, _m: i64) -> usize {
            offset_effect(n)
        }
    }

    fn sweep_step_cost(iter_cost: u64) {
        for p in 1..=2048usize {
            for m in -3i64..=20 {
                for pwr_step in [1usize, 2] {
                    for ext in [false, true] {
                        let (got, exp) = optimize_pair::<StepCost>(p, iter_cost, m, pwr_step, ext);
                        assert_eq!(
                            got, exp,
                            "p={p} m={m} pwr={pwr_step} ext={ext} iter={iter_cost}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn series_cost_optimize_matches_model_iter_cost_1() {
        sweep_step_cost(1);
    }

    #[test]
    fn series_cost_optimize_matches_model_iter_cost_1000() {
        sweep_step_cost(1000);
    }

    #[test]
    fn series_cost_optimize_step_zero_returns_zero_reductions() {
        // `floor(log2(p)) / 2` is 0 for `p < 4`. The search returns immediately:
        // count 0, `m_eff = reduction_effect(0, m)`, `niter = series_niter / pwr_step`.
        // `ext` is not read on this path.
        // p = 1, m = 5: effect = 5, denom = 5 - 2 = 3, niter = 1 / 3 = 0.
        // p = 1, m = -3: effect floors at 3, denom = 1, niter = 1.
        // p = 2, m = -3, pwr 1: denom = 1 - 0 + 3 - 2 = 2, niter = 2 / 2 = 1.
        // p = 2, m = -3, pwr 2: niter = 1 / 2 = 0.
        // p = 3, m = 0: effect floors at 3, denom = 1 + 3 - 2 = 2, niter = 3 / 2 = 1.
        // p = 3, m = 20: effect = 20, denom = 19, niter = 0.
        let rows = [
            (1usize, 5i64, 1usize, false, (0usize, 0usize, 5usize)),
            (1, -3, 1, true, (0, 1, 3)),
            (2, -3, 1, false, (0, 1, 3)),
            (2, -3, 2, true, (0, 0, 3)),
            (3, 0, 1, false, (0, 1, 3)),
            (3, 20, 2, true, (0, 0, 20)),
        ];
        for (p, m, pwr, ext, expect) in rows {
            assert_eq!(log2_floor_ref(p) / 2, 0, "p={p}");
            let (got, exp) = optimize_pair::<StepCost>(p, 1, m, pwr, ext);
            assert_eq!(exp, expect, "p={p} m={m} pwr={pwr} ext={ext}");
            assert_eq!(got, expect, "p={p} m={m} pwr={pwr} ext={ext}");
        }
    }

    #[test]
    fn series_cost_optimize_two_advances_keeps_the_last_improvement() {
        // p = 64: floor(log2) = 6, step = 3, m = 0, so the first count is 3.
        // Costs at 3, 6, 9, 12 are 100, 50, 10, 40. The series term is 0
        // (`m_eff = 64 + n >= 67`, denominator `m_eff + 2 > 64`).
        // Two strict improvements select count 9, which is the first count plus two steps.
        let p = 64usize;
        let (got, exp) = optimize_pair::<DropTwoThenRise>(p, 1, 0, 1, false);
        assert_eq!(exp, (9, 0, 73));
        assert_eq!(got, (9, 0, 73));
        let inp = Cand {
            p,
            iter_cost: 1,
            m: 0,
            div: 1,
            ext: false,
        };
        assert_eq!(candidate_ref::<DropTwoThenRise>(&inp, 3).2, 100);
        assert_eq!(candidate_ref::<DropTwoThenRise>(&inp, 6).2, 50);
        assert_eq!(candidate_ref::<DropTwoThenRise>(&inp, 9).2, 10);
        assert_eq!(candidate_ref::<DropTwoThenRise>(&inp, 12).2, 40);
    }

    #[test]
    fn series_cost_optimize_decreasing_cost_stops_at_p_plus_step() {
        // The reduction cost falls by 1 for each extra reduction, and the series term stays 0.
        // p = 64, step = 3, m = 0: first count 3, limit `p + step` = 67.
        // Counts 3, 6, ..., 66 are inside the guard. 66 + 3 = 69 is past it.
        // The rejected count 69 would have been strictly cheaper, so the guard is what stops the loop.
        let p = 64usize;
        let (got, exp) = optimize_pair::<DecreaseForever>(p, 1, 0, 1, false);
        assert_eq!(exp, (66, 0, 130));
        assert_eq!(got, (66, 0, 130));
        assert!(66 > 0);
        let step = 3usize;
        assert!(66 <= p + step);
        assert!(66 + step > p + step);
        let inp = Cand {
            p,
            iter_cost: 1,
            m: 0,
            div: 1,
            ext: false,
        };
        let held = candidate_ref::<DecreaseForever>(&inp, 66).2;
        let rejected = candidate_ref::<DecreaseForever>(&inp, 69).2;
        assert!(rejected < held, "held {held} rejected {rejected}");
    }

    #[test]
    fn series_cost_optimize_first_candidate_capped() {
        // p = 64, step = 3, first count = 3.
        // AllCapped: ext series part is `u64::MAX * niter`. m_eff = 11 gives
        // denom = 6 - 2 + 11 - 2 = 13, niter = 64 / 13 = 4. The same niter at the next
        // count is also 4 (`m_eff = 14`, denom = 16). Both costs are the cap, so the
        // first candidate is kept. Equal caps are not a strict improvement, and the count
        // is not backed up from 3 to 0.
        let (got, exp) = optimize_pair::<AllCapped>(64, u64::MAX, 0, 1, true);
        assert_eq!(exp, (3, 4, 11));
        assert_eq!(got, (3, 4, 11));
        let inp = Cand {
            p: 64,
            iter_cost: u64::MAX,
            m: 0,
            div: 1,
            ext: true,
        };
        assert_eq!(candidate_ref::<AllCapped>(&inp, 3).2, u64::MAX);
        assert_eq!(candidate_ref::<AllCapped>(&inp, 6).2, u64::MAX);

        // The first cost is the cap and the next exact cost is 50, then 80.
        // The cap is a real held cost: the search advances to count 6 and stops.
        // m_eff = 64 + 6 = 70 and niter is 0, so the cost is the reduction cost alone.
        let (got, exp) = optimize_pair::<CapThenLower>(64, 1, 0, 1, false);
        assert_eq!(exp, (6, 0, 70));
        assert_eq!(got, (6, 0, 70));
    }

    #[test]
    fn series_cost_optimize_capped_plateau_after_an_improvement() {
        // Counts 3, 6, 9 cost 500, 100, and `u64::MAX`. The search takes the improvement
        // to 6 and stops on the cap. Every later count is that same cap.
        let (got, exp) = optimize_pair::<ImproveThenCap>(64, 1, 0, 1, false);
        assert_eq!(exp, (6, 0, 70));
        assert_eq!(got, (6, 0, 70));
        let inp = Cand {
            p: 64,
            iter_cost: 1,
            m: 0,
            div: 1,
            ext: false,
        };
        assert_eq!(candidate_ref::<ImproveThenCap>(&inp, 3).2, 500);
        assert_eq!(candidate_ref::<ImproveThenCap>(&inp, 6).2, 100);
        assert_eq!(candidate_ref::<ImproveThenCap>(&inp, 9).2, u64::MAX);
        assert_eq!(candidate_ref::<ImproveThenCap>(&inp, 12).2, u64::MAX);
    }

    #[test]
    fn series_cost_optimize_negative_m_sets_the_first_count() {
        // p = 32: floor(log2) = 5, step = 2. m = -3, so the first count is 2 - (-3) = 5.
        // effect(5, -3) floors at 3. denom = 5 - 2 + 3 - 2 = 4, niter = 32 / 4 = 8.
        // The next count adds at least `2^32` of reduction cost, more than the series term, so the first count is kept.
        let (got, exp) = optimize_pair::<StepCost>(32, 1, -3, 1, false);
        assert_eq!(exp, (5, 8, 3));
        assert_eq!(got, (5, 8, 3));
        // p = 8: floor(log2) = 3, step = 1. First count = 1 - (-3) = 4.
        // effect floors 4 + (-3) = 1 up to 3. denom = 3 - 1 + 3 - 2 = 3, niter = 8 / 3 = 2.
        let (got, exp) = optimize_pair::<StepCost>(8, 1000, -3, 2, true);
        assert_eq!(exp, (4, 1, 3));
        assert_eq!(got, (4, 1, 3));
    }

    /// `2 * (mul(p) + add(p) + add(WORD_BIT_SIZE))`, then the `u64` clamp.
    /// `inc` is built at precision 1. `bit_len_to_word_len(1)` is one word, so
    /// `mantissa_max_bit_len` is `WORD_BIT_SIZE` for both `from_word(1, 1)` (sin, sinh)
    /// and `new(1)` (cos).
    fn fact_iter_oracle(p: usize) -> u128 {
        let word = add_clamped(crate::defs::WORD_BIT_SIZE);
        let per = mul_clamped(p) + add_clamped(p) + word;
        (per * 2).min(u128::from(u64::MAX))
    }

    fn check_fact(p: usize, built: Result<super::FactPolycoeffGen, crate::defs::Error>) {
        let expect = fact_iter_oracle(p);
        match built {
            Ok(gen) => assert_eq!(u128::from(gen.iter_cost()), expect, "p={p}"),
            Err(err) => {
                // `for_sin` / `for_cos` / `for_sinh` allocate two `p`-bit mantissas.
                // `p = 1 << 28` is 32 MiB and must succeed. A wider precision (64-bit
                // `MAX_PRECISION`) is refused with `MemoryAllocation` before `iter_cost` exists.
                assert!(p > (1usize << 28), "p={p} should allocate");
                assert!(matches!(err, crate::defs::Error::MemoryAllocation), "p={p}");
                if mul_clamped(p) == u128::from(u64::MAX) {
                    assert_eq!(expect, u128::from(u64::MAX), "p={p}");
                }
            }
        }
    }

    fn each_fact_precision(mut f: impl FnMut(usize)) {
        for p in [1usize, 2, 69, 70, 1624, 1625, 2_643_875, 2_643_876, 1usize << 28] {
            f(p);
        }
        f(crate::defs::MAX_PRECISION);
    }

    #[test]
    fn fact_iter_cost_matches_integer_oracle() {
        // The `u64` cap of `2 * (mul(p) + p + word)` needs `mul(p)` near `2^63`, which starts
        // around `p > 2^42`. That mantissa cannot be allocated, so the cap is checked on the
        // oracle when construction returns `MemoryAllocation` and `mul` is already clamped.
        each_fact_precision(|p| {
            check_fact(p, super::FactPolycoeffGen::for_sin(p));
            check_fact(p, super::FactPolycoeffGen::for_cos(p));
            check_fact(p, super::FactPolycoeffGen::for_sinh(p));
        });
    }
}
