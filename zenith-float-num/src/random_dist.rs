//! Distribution samplers on the crate RNG (`random` feature / tests).

use crate::common::test_rng::random;
use crate::defs::WORD_BIT_SIZE;
use crate::Consts;
use crate::Error;
use crate::ExactNum;
use crate::RoundingMode;
use crate::Word;
use crate::NAN;

/// Named sample count for the exponential-mean gold.
const RANDOM_EXP_MEAN_SAMPLES: usize = 10_000;

fn dist_nan() -> ExactNum {
    ExactNum::nan(Some(Error::InvalidArgument))
}

fn finite_pos(x: &ExactNum) -> bool {
    !x.is_zero() && x.is_positive() && !x.is_inf()
}

/// Unit interval sample in `[0, 1)` at precision `p` from the crate RNG.
fn random_unit(p: usize, rm: RoundingMode) -> ExactNum {
    let pw = p.max(WORD_BIT_SIZE);
    let nwords = pw.div_ceil(WORD_BIT_SIZE);
    let mut acc = ExactNum::from_u8(0, pw);
    for _ in 0..nwords {
        let w: Word = random();
        acc = acc.ldexp(WORD_BIT_SIZE as i32, pw, RoundingMode::None);
        acc = acc.add(&ExactNum::from_word(w, pw), pw, RoundingMode::None);
    }
    acc.ldexp(-((nwords * WORD_BIT_SIZE) as i32), p, rm)
}

/// `(0, 1]` — redraws a zero so `ln` is defined.
fn random_unit_open(p: usize, rm: RoundingMode) -> ExactNum {
    loop {
        let u = random_unit(p, rm);
        if !u.is_zero() {
            return u;
        }
    }
}

/// Which law [`ExactNumArray::random_fill`] samples.
#[derive(Clone, Debug)]
pub enum RandomDist {
    /// Uniform on `[a, b]`.
    Uniform(ExactNum, ExactNum),
    /// Gaussian \(N(\mu,\sigma^2)\) via Box–Muller. Distinct from
    /// [`ExactNum::random_normal`], which draws a random normalized mantissa.
    Normal(ExactNum, ExactNum),
    /// Exponential with rate `lambda` (mean `1/lambda`).
    Exponential(ExactNum),
}

impl ExactNum {
    /// Uniform sample on `[a, b]` at `(p, rm)`. `a > b` or a non-finite end is `NaN`.
    ///
    /// Uses the crate RNG ([`crate::reseed_random`] when the `random` feature is on).
    pub fn random_uniform(a: &Self, b: &Self, p: usize, rm: RoundingMode) -> Self {
        if a.is_nan() || b.is_nan() || a.is_inf() || b.is_inf() {
            return NAN;
        }
        match a.cmp(b) {
            Some(c) if c > 0 => return dist_nan(),
            Some(0) => {
                let mut v = a.clone();
                let _ = v.set_precision(p, rm);
                return v;
            }
            None => return NAN,
            _ => {}
        }
        let u = random_unit(p, RoundingMode::None);
        let span = b.sub(a, p, RoundingMode::None);
        a.add(&span.mul(&u, p, RoundingMode::None), p, rm)
    }

    /// Gaussian \(N(\mu,\sigma^2)\) via Box–Muller at `(p, rm)`.
    ///
    /// `sigma ≤ 0` is `NaN`. This is not [`Self::random_normal`] (mantissa draw).
    pub fn random_gaussian(
        mu: &Self,
        sigma: &Self,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Self {
        if mu.is_nan() || sigma.is_nan() || !finite_pos(sigma) {
            return dist_nan();
        }
        let pw = p.saturating_add(WORD_BIT_SIZE);
        let none = RoundingMode::None;
        let u1 = random_unit_open(pw, none);
        let u2 = random_unit(pw, none);
        let two = ExactNum::from_u8(2, pw);
        let r = two.neg().mul(&u1.ln(pw, none, cc), pw, none).sqrt(pw, none);
        let theta = two.mul(&cc.pi(pw, none), pw, none).mul(&u2, pw, none);
        let (s, _) = theta.sin_cos(pw, none, cc);
        let z = r.mul(&s, pw, none);
        mu.add(&sigma.mul(&z, pw, none), p, rm)
    }

    /// Exponential sample with rate `lambda` (inverse CDF \(-\ln U/\lambda\)).
    ///
    /// `lambda ≤ 0` is `NaN`.
    pub fn random_exponential(lambda: &Self, p: usize, rm: RoundingMode, cc: &mut Consts) -> Self {
        if lambda.is_nan() || !finite_pos(lambda) {
            return dist_nan();
        }
        let pw = p.saturating_add(WORD_BIT_SIZE);
        let none = RoundingMode::None;
        let u = random_unit_open(pw, none);
        u.ln(pw, none, cc).neg().div(lambda, p, rm)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::test_rng::reseed_random;
    use crate::ExactNumArray;

    #[test]
    fn random_dist_plan_golds() {
        let p = 256;
        let rm = RoundingMode::ToEven;
        let mut cc = Consts::new().unwrap();
        let zero = ExactNum::from_u8(0, p);
        let one = ExactNum::from_u8(1, p);

        reseed_random(1);
        let u = ExactNum::random_uniform(&zero, &one, p, rm);
        assert!(!u.is_negative());
        assert!(matches!(u.cmp(&one), Some(c) if c <= 0));
        assert!(!u.is_nan());

        reseed_random(7);
        let g1 = ExactNum::random_gaussian(&zero, &one, p, rm, &mut cc);
        let g2 = ExactNum::random_gaussian(&zero, &one, p, rm, &mut cc);
        reseed_random(7);
        let h1 = ExactNum::random_gaussian(&zero, &one, p, rm, &mut cc);
        let h2 = ExactNum::random_gaussian(&zero, &one, p, rm, &mut cc);
        assert_eq!(g1.cmp(&h1), Some(0));
        assert_eq!(g2.cmp(&h2), Some(0));

        let p64 = 64;
        let z64 = ExactNum::from_u8(0, p64);
        let o64 = ExactNum::from_u8(1, p64);
        reseed_random(3);
        let n = ExactNum::from_word(RANDOM_EXP_MEAN_SAMPLES as Word, p64);
        let mut sum = ExactNum::from_u8(0, p64);
        for _ in 0..RANDOM_EXP_MEAN_SAMPLES {
            sum = sum.add(
                &ExactNum::random_exponential(&o64, p64, rm, &mut cc),
                p64,
                rm,
            );
        }
        let mean = sum.div(&n, p64, rm);
        let err = mean.sub(&o64, p64, rm).abs();
        let tol = ExactNum::from_u8(1, p64).div(&ExactNum::from_u8(20, p64), p64, rm);
        assert!(matches!(err.cmp(&tol), Some(c) if c <= 0));

        reseed_random(11);
        let fill = ExactNumArray::random_fill(
            (100, 100),
            &RandomDist::Normal(z64.clone(), o64),
            p64,
            rm,
            &mut cc,
        )
        .unwrap();
        assert_eq!(fill.shape(), (100, 100));

        assert!(ExactNum::random_uniform(&one, &zero, p, rm).is_nan());
        assert!(ExactNum::random_gaussian(&zero, &zero, p, rm, &mut cc).is_nan());
        assert!(ExactNum::random_exponential(&z64, p64, rm, &mut cc).is_nan());
    }
}
