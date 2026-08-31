//! Dense univariate polynomials with [`ExactNum`] coefficients.

use crate::defs::RoundingMode;
use crate::defs::WORD_BIT_SIZE;
use crate::Consts;
use crate::ExactInt;
use crate::ExactNum;
use crate::ExactNumArray;
use alloc::vec;
use alloc::vec::Vec;

/// Highest degree whose companion eigenvalues have a closed form on this type.
///
/// Degree 1 is linear. Degree 2 is the characteristic polynomial of the 2×2
/// companion (quadratic formula). Higher-degree companions are not symmetric,
/// so [`ExactNumArray::eigen_decomp`] cannot be used; those degrees return
/// `None` from [`ExactNumPoly::roots_real`].
pub const POLY_COMPANION_CLOSED_DEG: usize = 2;

/// Dense univariate polynomial `c₀ + c₁ x + ⋯ + cₙ xⁿ`.
///
/// Coefficients are stored lowest degree first. Leading zeros are stripped.
/// The zero polynomial has an empty coefficient vector. Arithmetic methods
/// take an explicit precision `p` and rounding mode; stored `(p, rm)` are the
/// values used when the polynomial was constructed.
#[derive(Clone, Debug)]
pub struct ExactNumPoly {
    coeffs: Vec<ExactNum>,
    p: usize,
    rm: RoundingMode,
}

impl ExactNumPoly {
    /// Zero polynomial at `(p, rm)`.
    pub fn zero(p: usize, rm: RoundingMode) -> Self {
        Self {
            coeffs: Vec::new(),
            p,
            rm,
        }
    }

    /// Constant `1` at `(p, rm)`.
    pub fn one(p: usize, rm: RoundingMode) -> Self {
        Self::from_coeffs(p, rm, &[ExactNum::from_u8(1, p)])
    }

    /// Build from coefficients, lowest degree first. Each coefficient is
    /// rounded to `(p, rm)`. Leading zeros are stripped.
    pub fn from_coeffs(p: usize, rm: RoundingMode, coeffs: &[ExactNum]) -> Self {
        let mut out = Vec::with_capacity(coeffs.len());
        for c in coeffs {
            let mut y = c.clone();
            let _ = y.set_precision(p, rm);
            out.push(y);
        }
        let mut s = Self { coeffs: out, p, rm };
        s.strip_leading();
        s
    }

    /// Coefficients from `i64` values, lowest degree first.
    pub fn from_i64_coeffs(p: usize, rm: RoundingMode, coeffs: &[i64]) -> Self {
        let xs: Vec<ExactNum> = coeffs.iter().map(|&k| ExactNum::from_i64(k, p)).collect();
        Self::from_coeffs(p, rm, &xs)
    }

    /// Stored construction precision.
    pub fn precision(&self) -> usize {
        self.p
    }

    /// Stored construction rounding mode.
    pub fn rounding(&self) -> RoundingMode {
        self.rm
    }

    /// Coefficients, lowest degree first. Empty if this is the zero polynomial.
    pub fn coeffs(&self) -> &[ExactNum] {
        &self.coeffs
    }

    /// Degree `n` of a nonzero polynomial. `None` if this is zero.
    pub fn degree(&self) -> Option<usize> {
        if self.coeffs.is_empty() {
            None
        } else {
            Some(self.coeffs.len() - 1)
        }
    }

    /// True if every coefficient is zero (or the vector is empty).
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty() || self.coeffs.iter().all(|c| c.is_zero())
    }

    /// Coefficient of `x^k`, or zero if `k` exceeds the degree.
    pub fn coeff(&self, k: usize) -> ExactNum {
        self.coeffs
            .get(k)
            .cloned()
            .unwrap_or_else(|| ExactNum::new(self.p))
    }

    fn leading(&self) -> ExactNum {
        self.coeffs
            .last()
            .cloned()
            .unwrap_or_else(|| ExactNum::new(self.p))
    }

    fn all_finite(&self) -> bool {
        self.coeffs.iter().all(|c| !c.is_nan() && !c.is_inf())
    }

    fn strip_leading(&mut self) {
        while self.coeffs.last().is_some_and(|c| c.is_zero()) {
            self.coeffs.pop();
        }
    }

    fn work_p(p: usize) -> usize {
        p.saturating_add(WORD_BIT_SIZE)
    }

    /// Horner evaluation `c₀ + x(c₁ + x(c₂ + ⋯))` at `(p, rm)`.
    ///
    /// `cc` is accepted for signature uniformity with other numeric methods;
    /// Horner uses only add and fused multiply-add.
    pub fn eval(&self, x: &ExactNum, p: usize, rm: RoundingMode, _cc: &mut Consts) -> ExactNum {
        ExactNum::polyval(&self.coeffs, x, p, rm)
    }

    /// Coefficient-wise sum at precision `p`.
    pub fn add(&self, rhs: &Self, p: usize, rm: RoundingMode) -> Self {
        let n = self.coeffs.len().max(rhs.coeffs.len());
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let a = self.coeff(i);
            let b = rhs.coeff(i);
            out.push(a.add(&b, p, rm));
        }
        Self::from_coeffs(p, rm, &out)
    }

    /// Coefficient-wise difference at precision `p`.
    pub fn sub(&self, rhs: &Self, p: usize, rm: RoundingMode) -> Self {
        let n = self.coeffs.len().max(rhs.coeffs.len());
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let a = self.coeff(i);
            let b = rhs.coeff(i);
            out.push(a.sub(&b, p, rm));
        }
        Self::from_coeffs(p, rm, &out)
    }

    /// Schoolbook product at precision `p`.
    pub fn mul(&self, rhs: &Self, p: usize, rm: RoundingMode) -> Self {
        if self.is_zero() || rhs.is_zero() {
            return Self::zero(p, rm);
        }
        let n = self.coeffs.len() + rhs.coeffs.len() - 1;
        let mut out = vec![ExactNum::new(p); n];
        for (i, a) in self.coeffs.iter().enumerate() {
            for (j, b) in rhs.coeffs.iter().enumerate() {
                let t = a.mul(b, p, rm);
                out[i + j] = out[i + j].add(&t, p, rm);
            }
        }
        Self::from_coeffs(p, rm, &out)
    }

    /// Polynomial division: `self = q·g + r` with `deg r < deg g`.
    ///
    /// `None` if `g` is the zero polynomial or a coefficient is non-finite.
    pub fn div_rem(&self, g: &Self, p: usize, rm: RoundingMode) -> Option<(Self, Self)> {
        if g.is_zero() || !self.all_finite() || !g.all_finite() {
            return None;
        }
        let deg_g = g.degree()?;
        let lc_g = g.leading();
        if lc_g.is_zero() {
            return None;
        }
        let mut r = self.coeffs.clone();
        for c in &mut r {
            let _ = c.set_precision(p, rm);
        }
        while r.last().is_some_and(|c| c.is_zero()) {
            r.pop();
        }
        if r.is_empty() {
            return Some((Self::zero(p, rm), Self::zero(p, rm)));
        }
        let deg_f = r.len() - 1;
        let mut q = if deg_f >= deg_g {
            vec![ExactNum::new(p); deg_f - deg_g + 1]
        } else {
            Vec::new()
        };
        while r.len() > deg_g {
            let deg_r = r.len() - 1;
            let shift = deg_r - deg_g;
            let t = r[deg_r].div(&lc_g, p, rm);
            q[shift] = q[shift].add(&t, p, rm);
            for i in 0..=deg_g {
                let term = t.mul(&g.coeff(i), p, rm);
                let idx = i + shift;
                r[idx] = r[idx].sub(&term, p, rm);
            }
            while r.last().is_some_and(|c| c.is_zero()) {
                r.pop();
            }
        }
        Some((Self::from_coeffs(p, rm, &q), Self::from_coeffs(p, rm, &r)))
    }

    fn integer_content(&self) -> Option<ExactInt> {
        if self.coeffs.is_empty() {
            return Some(ExactInt::zero());
        }
        let mut g: Option<ExactInt> = None;
        for c in &self.coeffs {
            if c.is_nan() || c.is_inf() || !c.fract().is_zero() {
                return None;
            }
            let Some(i) = ExactInt::from_exact_num(c) else {
                return None;
            };
            g = Some(match g {
                None => i,
                Some(prev) => prev.gcd(&i),
            });
        }
        g
    }

    /// Primitive part with positive leading coefficient, then monic.
    ///
    /// Integer content (GCD of integer coefficients) is divided out when every
    /// coefficient is an integer. Otherwise only the leading coefficient is
    /// removed (float monic form).
    fn primitive_monic(&self, p: usize, rm: RoundingMode) -> Self {
        if self.is_zero() {
            return Self::zero(p, rm);
        }
        let mut s = if let Some(cont) = self.integer_content() {
            if !cont.is_zero() && !cont.is_one() {
                let d = cont.to_exact_num(p, rm);
                let xs: Vec<ExactNum> = self.coeffs.iter().map(|c| c.div(&d, p, rm)).collect();
                Self::from_coeffs(p, rm, &xs)
            } else {
                self.clone()
            }
        } else {
            self.clone()
        };
        let lc = s.leading();
        if lc.is_zero() || lc.is_nan() || lc.is_inf() {
            return s;
        }
        if lc.is_negative() {
            s.coeffs = s.coeffs.into_iter().map(|c| c.neg()).collect();
        }
        let lc = s.leading();
        if lc.is_zero() {
            return s;
        }
        let xs: Vec<ExactNum> = s.coeffs.iter().map(|c| c.div(&lc, p, rm)).collect();
        Self::from_coeffs(p, rm, &xs)
    }

    /// Euclidean GCD with integer content removal, returned monic.
    ///
    /// Both zero → zero. One zero → monic of the other.
    pub fn gcd(&self, other: &Self, p: usize, rm: RoundingMode) -> Self {
        if self.is_zero() && other.is_zero() {
            return Self::zero(p, rm);
        }
        let mut a = self.primitive_monic(p, rm);
        let mut b = other.primitive_monic(p, rm);
        let max_steps = a.degree().unwrap_or(0) + b.degree().unwrap_or(0) + 1;
        let mut steps = 0;
        while !b.is_zero() && steps < max_steps {
            steps += 1;
            let Some((_, r)) = a.div_rem(&b, p, rm) else {
                break;
            };
            a = b;
            b = r.primitive_monic(p, rm);
        }
        a.primitive_monic(p, rm)
    }

    /// Composition `self(g(x))` by Horner at precision `p`.
    pub fn compose(&self, g: &Self, p: usize, rm: RoundingMode, _cc: &mut Consts) -> Self {
        if self.is_zero() {
            return Self::zero(p, rm);
        }
        let mut acc = Self::from_coeffs(p, rm, &[self.leading()]);
        for c in self.coeffs.iter().rev().skip(1) {
            acc = acc.mul(g, p, rm);
            acc = acc.add(&Self::from_coeffs(p, rm, &[c.clone()]), p, rm);
        }
        acc
    }

    /// Formal derivative. The derivative of a constant is zero.
    pub fn derivative(&self, p: usize, rm: RoundingMode) -> Self {
        if self.coeffs.len() <= 1 {
            return Self::zero(p, rm);
        }
        let mut out = Vec::with_capacity(self.coeffs.len() - 1);
        for (k, c) in self.coeffs.iter().enumerate().skip(1) {
            let n = ExactNum::from_u32(k as u32, p);
            out.push(c.mul(&n, p, rm));
        }
        Self::from_coeffs(p, rm, &out)
    }

    /// Indefinite integral with constant term zero.
    pub fn integral(&self, p: usize, rm: RoundingMode) -> Self {
        if self.is_zero() {
            return Self::zero(p, rm);
        }
        let mut out = Vec::with_capacity(self.coeffs.len() + 1);
        out.push(ExactNum::new(p));
        for (k, c) in self.coeffs.iter().enumerate() {
            let den = ExactNum::from_u32((k + 1) as u32, p);
            out.push(c.div(&den, p, rm));
        }
        Self::from_coeffs(p, rm, &out)
    }

    /// Monic companion matrix of this polynomial, or `None` if the degree is
    /// zero or the polynomial is zero / non-finite.
    ///
    /// Last column is `(-c₀/cₙ, …, -cₙ₋₁/cₙ)`. Subdiagonal is ones.
    pub fn companion_matrix(&self, p: usize, rm: RoundingMode) -> Option<ExactNumArray> {
        let n = self.degree()?;
        if n == 0 || !self.all_finite() {
            return None;
        }
        let lc = self.leading();
        if lc.is_zero() {
            return None;
        }
        let mut vals = vec![ExactNum::new(p); n * n];
        let one = ExactNum::from_u8(1, p);
        for i in 1..n {
            vals[i * n + (i - 1)] = one.clone();
        }
        for i in 0..n {
            let a = self.coeff(i).div(&lc, p, rm);
            vals[i * n + (n - 1)] = a.neg();
        }
        ExactNumArray::from_shape(p, n, n, &vals)
    }

    /// Real roots via the companion characteristic equation.
    ///
    /// Degree 1 and 2 use the closed-form eigenvalues of the companion
    /// (linear solve / quadratic formula) at extra working precision, then
    /// one round to `p`. Degree greater than [`POLY_COMPANION_CLOSED_DEG`]
    /// returns `None` — those companions are not symmetric, so
    /// [`ExactNumArray::eigen_decomp`] does not apply.
    ///
    /// The zero polynomial returns `None`. A nonzero constant returns an
    /// empty vector. A negative discriminant returns an empty vector.
    pub fn roots_real(
        &self,
        p: usize,
        rm: RoundingMode,
        _cc: &mut Consts,
    ) -> Option<Vec<ExactNum>> {
        if !self.all_finite() {
            return None;
        }
        let n = match self.degree() {
            None => return None,
            Some(0) => return Some(Vec::new()),
            Some(n) => n,
        };
        if n > POLY_COMPANION_CLOSED_DEG {
            return None;
        }
        let wrk = Self::work_p(p);
        let lc = self.leading();
        if n == 1 {
            let r = self.coeff(0).div(&lc, wrk, RoundingMode::None).neg();
            let mut r = r;
            let _ = r.set_precision(p, rm);
            return Some(vec![r]);
        }
        // Degree 2: a x² + b x + c = 0. Companion eigenvalues.
        let a = self.coeff(2);
        let b = self.coeff(1);
        let c = self.coeff(0);
        let four = ExactNum::from_u8(4, wrk);
        let disc = b.mul(&b, wrk, RoundingMode::None).sub(
            &four
                .mul(&a, wrk, RoundingMode::None)
                .mul(&c, wrk, RoundingMode::None),
            wrk,
            RoundingMode::None,
        );
        if disc.is_negative() {
            return Some(Vec::new());
        }
        let two_a = ExactNum::from_u8(2, wrk).mul(&a, wrk, RoundingMode::None);
        if two_a.is_zero() {
            return None;
        }
        let sqrt_d = disc.sqrt(wrk, RoundingMode::None);
        let mut r0 =
            b.neg()
                .sub(&sqrt_d, wrk, RoundingMode::None)
                .div(&two_a, wrk, RoundingMode::None);
        let mut r1 =
            b.neg()
                .add(&sqrt_d, wrk, RoundingMode::None)
                .div(&two_a, wrk, RoundingMode::None);
        let _ = r0.set_precision(p, rm);
        let _ = r1.set_precision(p, rm);
        if disc.is_zero() {
            Some(vec![r0])
        } else {
            Some(vec![r0, r1])
        }
    }
}

impl PartialEq for ExactNumPoly {
    fn eq(&self, other: &Self) -> bool {
        if self.coeffs.len() != other.coeffs.len() {
            return false;
        }
        self.coeffs
            .iter()
            .zip(other.coeffs.iter())
            .all(|(a, b)| a == b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Consts;

    fn gold_p() -> (usize, RoundingMode) {
        (256, RoundingMode::ToEven)
    }

    #[test]
    fn exact_num_poly_div_rem_gcd_compose_deriv_roots() {
        let (p, rm) = gold_p();
        let mut cc = Consts::new().expect("consts");

        let x2m1 = ExactNumPoly::from_i64_coeffs(p, rm, &[-1, 0, 1]);
        let xm1 = ExactNumPoly::from_i64_coeffs(p, rm, &[-1, 1]);
        let (q, r) = x2m1.div_rem(&xm1, p, rm).expect("div_rem");
        assert_eq!(q, ExactNumPoly::from_i64_coeffs(p, rm, &[1, 1]));
        assert!(r.is_zero());

        let g = x2m1.gcd(&xm1, p, rm);
        assert_eq!(g, xm1);

        let x2 = ExactNumPoly::from_i64_coeffs(p, rm, &[0, 0, 1]);
        let xp1 = ExactNumPoly::from_i64_coeffs(p, rm, &[1, 1]);
        let composed = x2.compose(&xp1, p, rm, &mut cc);
        assert_eq!(composed, ExactNumPoly::from_i64_coeffs(p, rm, &[1, 2, 1]));

        let x3 = ExactNumPoly::from_i64_coeffs(p, rm, &[0, 0, 0, 1]);
        let dx = x3.derivative(p, rm);
        assert_eq!(dx, ExactNumPoly::from_i64_coeffs(p, rm, &[0, 0, 3]));

        let x2m2 = ExactNumPoly::from_i64_coeffs(p, rm, &[-2, 0, 1]);
        let roots = x2m2.roots_real(p, rm, &mut cc).expect("roots");
        assert_eq!(roots.len(), 2);
        let s2 = ExactNum::from_u8(2, p).sqrt(p, rm);
        let ns2 = s2.neg();
        let mut saw_pos = false;
        let mut saw_neg = false;
        for root in &roots {
            if root.cmp(&s2) == Some(0) {
                saw_pos = true;
            }
            if root.cmp(&ns2) == Some(0) {
                saw_neg = true;
            }
        }
        assert!(saw_pos && saw_neg);

        let c = x2m2.companion_matrix(p, rm).expect("companion");
        assert_eq!(c.shape(), (2, 2));
        let two = ExactNum::from_u8(2, p);
        let one = ExactNum::from_u8(1, p);
        let z = ExactNum::new(p);
        assert_eq!(c.get2(0, 0).map(|x| x.cmp(&z)), Some(Some(0)));
        assert_eq!(c.get2(0, 1).map(|x| x.cmp(&two)), Some(Some(0)));
        assert_eq!(c.get2(1, 0).map(|x| x.cmp(&one)), Some(Some(0)));
        assert_eq!(c.get2(1, 1).map(|x| x.cmp(&z)), Some(Some(0)));

        assert!(x2m1.div_rem(&ExactNumPoly::zero(p, rm), p, rm).is_none());
        assert!(x3.roots_real(p, rm, &mut cc).is_none());
    }
}
