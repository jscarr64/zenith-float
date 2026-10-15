// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! Static constants.

#[cfg(feature = "std")]
use crate::defs::Error;
use crate::{num::ExactNumNumber, Exponent, WORD_BIT_SIZE};

#[cfg(feature = "std")]
use crate::ops::consts::Consts;
#[cfg(feature = "std")]
use core::cell::RefCell;

macro_rules! small_consts {
    ($($(#[$doc:meta])* $name:ident = $n:expr;)*) => {
        $(
            $(#[$doc])*
            pub(crate) static $name: ExactNumNumber = ExactNumNumber::small_const($n);
        )*
    };
}

small_consts! {
    /// 1
    ONE = 1;
    /// 2
    TWO = 2;
    /// 4
    FOUR = 4;
    /// 3
    THREE = 3;
    /// 10
    TEN = 10;
}

// Thread-local Consts cache used by decimal conversion under std. Created on first use, so a
// failed allocation is an `Err` for the caller instead of a panic in a thread-local initialiser.
#[cfg(feature = "std")]
thread_local! {
    static TENPOWERS: RefCell<Option<Consts>> = const { RefCell::new(None) };
}

/// [`with_tenpowers`] for `fmt` implementations: a cache that cannot be created is `fmt::Error`.
#[cfg(feature = "std")]
pub(crate) fn with_tenpowers_fmt(
    f: impl FnOnce(&mut Consts) -> core::fmt::Result,
) -> core::fmt::Result {
    let mut out = Ok(());
    with_tenpowers(|cc| {
        out = f(cc);
        Ok(())
    })
    .map_err(|_| core::fmt::Error)?;
    out
}

/// Runs `f` with this thread's constants cache, creating it on first use.
///
/// # Errors
///
/// - MemoryAllocation: the cache could not be created.
/// - Any error returned by `f`.
#[cfg(feature = "std")]
pub fn with_tenpowers<R>(f: impl FnOnce(&mut Consts) -> Result<R, Error>) -> Result<R, Error> {
    TENPOWERS.with(|cell| {
        let mut slot = cell.borrow_mut();
        let cc = match slot.as_mut() {
            Some(cc) => cc,
            None => slot.insert(Consts::new()?),
        };
        f(cc)
    })
}

pub const TRIG_EXP_THRES: Exponent = -(WORD_BIT_SIZE as Exponent);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::buf::INLINE_WORDS;
    use crate::defs::{Word, DEFAULT_P};

    // DEFAULT_P bits must fit inline on every target, so building the constants cannot allocate.
    const _: () = assert!(INLINE_WORDS * WORD_BIT_SIZE >= DEFAULT_P);

    // The constants are `const`-constructible: no lazy initialisation, no Result.
    const SEVEN: ExactNumNumber = ExactNumNumber::small_const(7);

    fn same(a: &ExactNumNumber, b: &ExactNumNumber) -> bool {
        a.exponent() == b.exponent()
            && a.sign() == b.sign()
            && a.mantissa().digits() == b.mantissa().digits()
            && a.mantissa().bit_len() == b.mantissa().bit_len()
            && a.mantissa().max_bit_len() == b.mantissa().max_bit_len()
    }

    #[test]
    fn constants_equal_from_word_and_stay_inline() {
        let all: [(&ExactNumNumber, Word); 5] =
            [(&ONE, 1), (&TWO, 2), (&THREE, 3), (&FOUR, 4), (&TEN, 10)];
        for (c, n) in all {
            let expect = ExactNumNumber::from_word(n, DEFAULT_P);
            assert!(matches!(&expect, Ok(e) if same(c, e)), "constant {n}");
            assert!(c.is_inline(), "constant {n} must not allocate");
        }
        let seven = ExactNumNumber::from_word(7, DEFAULT_P);
        assert!(matches!(&seven, Ok(e) if same(&SEVEN, e)));
    }

    #[test]
    fn small_const_zero_is_zero() {
        let z = ExactNumNumber::small_const(0);
        assert!(z.mantissa().is_zero());
        assert_eq!(z.exponent(), 0);
    }

    #[cfg(feature = "std")]
    #[test]
    fn tenpowers_is_fallible_and_reusable() {
        assert_eq!(with_tenpowers(|_| Ok(5usize)), Ok(5));
        assert_eq!(with_tenpowers(|_| Ok(6usize)), Ok(6));
        let e: Result<(), crate::defs::Error> =
            with_tenpowers(|_| Err(crate::defs::Error::InvalidArgument));
        assert_eq!(e, Err(crate::defs::Error::InvalidArgument));
    }
}
