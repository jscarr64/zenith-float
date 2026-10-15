// Derived from astro-float (https://github.com/stencillogic/astro-float),
// Copyright (c) 2022 stencillogic, MIT License.
//! Static constants.

use crate::{defs::DEFAULT_P, num::ExactNumNumber, Exponent, WORD_BIT_SIZE};
use lazy_static::lazy_static;

#[cfg(feature = "std")]
use crate::ops::consts::Consts;
#[cfg(feature = "std")]
use core::cell::RefCell;

lazy_static! {

    /// 1
    pub(crate) static ref ONE: ExactNumNumber = ExactNumNumber::from_word(1, DEFAULT_P).expect("Constant ONE initialization.");

    /// 2
    pub(crate) static ref TWO: ExactNumNumber = ExactNumNumber::from_word(2, DEFAULT_P).expect("Constant TWO initialization.");

    /// 4
    pub(crate) static ref FOUR: ExactNumNumber = ExactNumNumber::from_word(4, DEFAULT_P).expect("Constant FOUR initialization.");

    /// 6
    pub(crate) static ref SIX: ExactNumNumber = ExactNumNumber::from_word(6, DEFAULT_P).expect("Constant SIX initialization.");

    /// 8
    pub(crate) static ref EIGHT: ExactNumNumber = ExactNumNumber::from_word(8, DEFAULT_P).expect("Constant EIGHT initialization.");

    /// 16
    pub(crate) static ref SIXTEEN: ExactNumNumber = ExactNumNumber::from_word(16, DEFAULT_P).expect("Constant SIXTEEN initialization.");

    /// 3
    pub(crate) static ref THREE: ExactNumNumber = ExactNumNumber::from_word(3, DEFAULT_P).expect("Constant THREE initialization.");

    /// 5
    pub(crate) static ref FIVE: ExactNumNumber = ExactNumNumber::from_word(5, DEFAULT_P).expect("Constant FIVE initialization.");

    /// 10
    pub(crate) static ref TEN: ExactNumNumber = ExactNumNumber::from_word(10, DEFAULT_P).expect("Constant TEN initialization.");

    /// 15
    pub(crate) static ref FIFTEEN: ExactNumNumber = ExactNumNumber::from_word(15, DEFAULT_P).expect("Constant FIFTEEN initialization.");

    /// 24
    pub(crate) static ref C24: ExactNumNumber = ExactNumNumber::from_word(24, DEFAULT_P).expect("Constant C24 initialization.");

    /// 40
    pub(crate) static ref FOURTY: ExactNumNumber = ExactNumNumber::from_word(40, DEFAULT_P).expect("Constant FOURTY initialization.");

    /// 120
    pub(crate) static ref C120: ExactNumNumber = ExactNumNumber::from_word(120, DEFAULT_P).expect("Constant C120 initialization.");
}

// Thread-local Consts cache used by decimal conversion under std.
#[cfg(feature = "std")]
thread_local! {
    pub static TENPOWERS: RefCell<Consts> = RefCell::new(Consts::new().expect("Failed to initialize thread-local constants cache"));
}

pub const TRIG_EXP_THRES: Exponent = -(WORD_BIT_SIZE as Exponent);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::buf::INLINE_WORDS;
    use crate::defs::Word;

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
        let all: [(&ExactNumNumber, Word); 13] = [
            (&ONE, 1),
            (&TWO, 2),
            (&THREE, 3),
            (&FOUR, 4),
            (&FIVE, 5),
            (&SIX, 6),
            (&EIGHT, 8),
            (&TEN, 10),
            (&FIFTEEN, 15),
            (&SIXTEEN, 16),
            (&C24, 24),
            (&FOURTY, 40),
            (&C120, 120),
        ];
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
