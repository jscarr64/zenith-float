//! zenith-float implements arbitrary-precision software floating-point numbers.

#![cfg_attr(not(feature = "std"), no_std)]
#![deny(missing_docs)]
#![deny(unused)]
#![deny(clippy::suspicious)]
#![deny(clippy::float_arithmetic)]
#![allow(clippy::manual_div_ceil)]
#![allow(clippy::manual_is_multiple_of)]
#![allow(clippy::neg_cmp_op_on_partial_ord)]
#![allow(clippy::len_zero)]
#![allow(clippy::slow_vector_initialization)]
#![allow(clippy::comparison_chain)]
#![allow(clippy::collapsible_else_if)]
#![allow(clippy::collapsible_if)]

extern crate alloc;

mod ball;
mod chebyshev;
mod common;
mod complex;
mod complex_airy;
mod complex_bessel;
mod complex_ei;
mod complex_elliptic;
mod complex_hypergeom;
mod complex_special;
mod conv;
pub mod ctx;
mod defs;
mod dist;
mod dsp;
mod ext;
mod hash;
mod ieee_soft;
mod integer;
mod mantissa;
mod modular;
mod num;
mod ode;
mod ops;
mod orthopoly;
mod parser;
mod poly;
mod quadrature;
mod radix_float;
#[cfg(any(test, feature = "random"))]
mod random_dist;
mod rational;
mod roots;
mod strop;

#[cfg(feature = "std")]
mod for_3rd;

#[doc(hidden)]
pub mod macro_util;

pub use crate::ball::ziv_round;
pub use crate::ball::ziv_round_vec;
pub use crate::ball::Ball;
pub use crate::ball::ComplexBall;
pub use crate::chebyshev::chebyshev_coeffs;
pub use crate::chebyshev::chebyshev_error_bound;
pub use crate::chebyshev::chebyshev_eval;
pub use crate::chebyshev::clenshaw;
pub use crate::chebyshev::CHEBYSHEV_MAX_DEGREE;
pub use crate::common::buf::INLINE_WORDS;
pub use crate::complex::ExactComplex;
pub use crate::defs::Error;
pub use crate::defs::Exponent;
pub use crate::defs::Radix;
pub use crate::defs::RoundingMode;
pub use crate::defs::Sign;
pub use crate::defs::Word;
pub use crate::dsp::blackman_window;
pub use crate::dsp::dct;
pub use crate::dsp::dst;
pub use crate::dsp::fft_real;
pub use crate::dsp::hamming_window;
pub use crate::dsp::hann_window;
pub use crate::dsp::idct;
pub use crate::dsp::idst;
pub use crate::dsp::ifft_real;
pub use crate::dsp::kaiser_window;
pub use crate::dsp::rectangular_window;
pub use crate::dsp::DSP_MAX_POINTS;
pub use crate::ext::ExactNum;
pub use crate::ext::FromExt;
pub use crate::ext::INF_NEG;
pub use crate::ext::INF_POS;
pub use crate::ext::NAN;
pub use crate::hash::constant_time_eq;
pub use crate::hash::hmac_sha256;
pub use crate::hash::sha256;
pub use crate::hash::sha512;
pub use crate::ieee_soft::{
    ExactNumArray, Ieee32, Ieee32Array, Ieee64, Ieee64Array, IEEE_SIMD_LANE_WIDTH,
};
pub use crate::integer::ExactInt;
pub use crate::modular::miller_rabin;
pub use crate::modular::mod_inv;
pub use crate::modular::mod_pow;
pub use crate::modular::pollard_rho;
pub use crate::modular::POLLARD_RHO_ITER_MAX;
pub use crate::ode::euler;
pub use crate::ode::ode_min_step;
pub use crate::ode::rk4;
pub use crate::ode::rk45_adaptive;
pub use crate::ode::ODE_MAX_STEPS;
pub use crate::ode::ODE_MIN_STEP;
pub use crate::ops::consts::CachedFBig;
pub use crate::ops::consts::ConstCache;
pub use crate::ops::consts::ConstCacheInfo;
pub use crate::ops::consts::Consts;
#[cfg(feature = "std")]
pub use crate::ops::consts::SharedConsts;
pub use crate::orthopoly::ORTHOPOLY_N_MAX;
pub use crate::poly::ExactNumPoly;
pub use crate::poly::POLY_COMPANION_CLOSED_DEG;
pub use crate::quadrature::gauss_hermite;
pub use crate::quadrature::gauss_laguerre;
pub use crate::quadrature::gauss_legendre;
pub use crate::quadrature::tanh_sinh;
pub use crate::quadrature::QUADRATURE_MAX_NODES;
pub use crate::quadrature::TANH_SINH_LEVELS_MAX;
pub use crate::radix_float::RadixFloat;
#[cfg(any(test, feature = "random"))]
pub use crate::random_dist::RandomDist;
pub use crate::rational::ExactRational;
pub use crate::roots::bisect;
pub use crate::roots::brent;
pub use crate::roots::illinois;
pub use crate::roots::newton;
pub use crate::roots::root_default_tol;
pub use crate::roots::ROOT_DEFAULT_TOL;
pub use crate::roots::ROOT_MAX_ITER;

pub use crate::defs::EXPONENT_BIT_SIZE;
pub use crate::defs::EXPONENT_MAX;
pub use crate::defs::EXPONENT_MIN;
pub use crate::defs::WORD_BASE;
pub use crate::defs::WORD_BIT_SIZE;
pub use crate::defs::WORD_MAX;
pub use crate::defs::WORD_SIGNIFICANT_BIT;

pub use crate::common::util::MAX_PREC_RETRY;

#[cfg(feature = "random")]
pub use crate::common::test_rng::{random_seed, reseed_random, seeded_random, DEFAULT_RANDOM_SEED};

#[cfg(test)]
mod tests {

    #[test]
    fn test_bigfloat() {
        use crate::Consts;
        use crate::ExactNum;
        use crate::RoundingMode;

        // Precision with some space for error.
        let p = 1024 + 8;

        // Rounding of all operations
        let rm = RoundingMode::ToEven;

        // Initialize mathematical constants cache
        let mut cc = Consts::new().expect("An error occured when initializing constants");

        // Compute pi: pi = 6*arctan(1/sqrt(3))
        let six = ExactNum::from_word(6, 1);
        let three = ExactNum::from_word(3, p);

        let n = three.sqrt(p, rm);
        let n = n.reciprocal(p, rm);
        let n = n.atan(p, rm, &mut cc);
        let mut pi = six.mul(&n, p, rm);

        // Reduce precision to 1024
        pi.set_precision(1024, rm).expect("Precision updated");

        // Use library's constant for verifying the result
        let pi_lib = cc.pi_num(1024, rm).unwrap().into();

        // Compare computed constant with library's constant
        assert_eq!(pi.cmp(&pi_lib), Some(0));

        // Print using decimal radix.
        //println!("{}", pi);
    }
}
