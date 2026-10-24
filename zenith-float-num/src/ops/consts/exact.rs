//! Constants correctly rounded at exactly `p` bits: [`Consts::get_exact_prec`].

use crate::common::util::round_p;
use crate::defs::{MAX_PRECISION, MIN_PRECISION};
use crate::num::ExactNumNumber;
use crate::ops::consts::Consts;
use crate::Error;
use crate::ExactNum;
use crate::Exponent;
use crate::RoundingMode;
use crate::WORD_BIT_SIZE;

/// A constant [`Consts::get_exact_prec`] can return.
#[derive(Debug, PartialEq, Eq)]
pub enum ConstantId {
    /// π.
    Pi,
    /// Euler's number e.
    E,
    /// ln 2.
    Ln2,
    /// ln 10.
    Ln10,
    /// √2.
    Sqrt2,
    /// The golden ratio φ = (1+√5)/2.
    Phi,
    /// The Euler–Mascheroni constant γ.
    EulerGamma,
}

/// Ziv rounds before giving up; each doubles the working precision. All seven constants
/// are irrational, so the bracket eventually avoids every p-bit rounding boundary.
const ZIV_ROUNDS: usize = 16;

impl Consts {
    /// Whether constant `c` is exactly representable at some finite precision.
    /// `false` for all seven: each is irrational, so every result of
    /// [`Self::get_exact_prec`] is flagged inexact.
    pub fn is_exact(c: ConstantId) -> bool {
        match c {
            ConstantId::Pi
            | ConstantId::E
            | ConstantId::Ln2
            | ConstantId::Ln10
            | ConstantId::Sqrt2
            | ConstantId::Phi
            | ConstantId::EulerGamma => false,
        }
    }

    /// Constant `c` correctly rounded at exactly `p` bits with `rm`.
    ///
    /// ## Contract
    ///
    /// - The result is the p-bit number `rm` selects for the true value of `c`, for every
    ///   `p ≥ 1` (including 1, 2, 63, 64, 65 and 65536): `RoundingMode::Down ≤ c ≤
    ///   RoundingMode::Up` for EVERY `p ≥ 1`, and since `c` is irrational `Down < c < Up`,
    ///   with `Up` the next p-bit number above `Down`. `ToZero` is `Down` and `FromZero` is
    ///   `Up` (all seven are positive); `ToEven` and `ToOdd` agree (no constant is a tie).
    /// - Bits beyond `p` in the stored words are zero; both inexact flags are set
    ///   ([`Self::is_exact`] is false for every constant).
    /// - Bit-identical whatever the cache holds and in whatever order constants and
    ///   precisions are requested: the result depends only on `(c, p, rm)`. The caches are
    ///   used for speed only; a value from them is never returned without the Ziv check below.
    ///
    /// Method: an approximation `a` at `w + 64` bits (w ≥ p + 64) is trusted only to
    /// `|a − c| ≤ 2^(e(a) − w)`, i.e. 2^64 of its ulps. Both ends of that interval are
    /// computed exactly and rounded at `p` with `rm`; when they agree that is the result,
    /// otherwise `w` doubles (at most 16 times).
    ///
    /// ## Errors
    ///
    ///  - InvalidArgument: `p` is outside `MIN_PRECISION..=MAX_PRECISION`, or `rm` is `RoundingMode::None` (no
    ///    rounding cannot give a value at exactly `p` bits); the bracket did not settle.
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    pub fn get_exact_prec(
        &mut self,
        c: ConstantId,
        p: usize,
        rm: RoundingMode,
    ) -> Result<ExactNum, Error> {
        if !(MIN_PRECISION..=MAX_PRECISION).contains(&p) || rm == RoundingMode::None {
            return Err(Error::InvalidArgument);
        }
        let mut w = round_p(p)
            .checked_add(WORD_BIT_SIZE)
            .ok_or(Error::InvalidArgument)?;
        for _ in 0..ZIV_ROUNDS {
            let q = w.checked_add(WORD_BIT_SIZE).ok_or(Error::InvalidArgument)?;
            let a = self.approx(&c, q)?;
            let lo = Self::bracket_end(&a, w, false, p, rm)?;
            let hi = Self::bracket_end(&a, w, true, p, rm)?;
            if lo.cmp(&hi) == 0 {
                let mut r = lo;
                r.set_inexact(true);
                return Ok(r.into());
            }
            w = w.checked_mul(2).ok_or(Error::InvalidArgument)?;
        }
        Err(Error::InvalidArgument)
    }

    /// `c` at `q` bits (word-rounded), within a few of its ulps.
    fn approx(&mut self, c: &ConstantId, q: usize) -> Result<ExactNumNumber, Error> {
        let rm = RoundingMode::ToEven;
        match c {
            ConstantId::Pi => self.pi_num(q, rm),
            ConstantId::E => self.e_num(q, rm),
            ConstantId::Ln2 => self.ln_2_num(q, rm),
            ConstantId::Ln10 => self.ln_10_num(q, rm),
            ConstantId::Sqrt2 => self.sqrt2_num(q, rm),
            ConstantId::Phi => self.phi_num(q, rm),
            ConstantId::EulerGamma => self.euler_gamma_num(q, rm),
        }
    }

    /// `a ± 2^(e(a) − w)`, computed exactly, then rounded at `p` bits with `rm`.
    fn bracket_end(
        a: &ExactNumNumber,
        w: usize,
        up: bool,
        p: usize,
        rm: RoundingMode,
    ) -> Result<ExactNumNumber, Error> {
        let mut eps = ExactNumNumber::from_word(1, WORD_BIT_SIZE)?;
        // from_word(1) is 0.1b · 2^1; exponent e(a) − w + 1 makes it 2^(e(a) − w).
        let shift = Exponent::try_from(w).map_err(|_| Error::InvalidArgument)?;
        let e = a
            .exponent()
            .checked_sub(shift)
            .and_then(|x| x.checked_add(1))
            .ok_or(Error::InvalidArgument)?;
        eps.set_exponent(e);
        // `a` has at most round_p(w + 64) bits and eps lies on its grid, so one more word
        // holds the sum or difference exactly.
        let wide = round_p(a.mantissa_max_bit_len().max(w + WORD_BIT_SIZE))
            .checked_add(WORD_BIT_SIZE)
            .ok_or(Error::InvalidArgument)?;
        let mut end = if up {
            a.add(&eps, wide, RoundingMode::None)?
        } else {
            a.sub(&eps, wide, RoundingMode::None)?
        };
        end.set_precision(p, rm)?;
        Ok(end)
    }
}
