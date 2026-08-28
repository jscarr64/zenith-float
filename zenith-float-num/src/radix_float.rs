//! Arbitrary-radix float wrapper around [`ExactNum`].
//!
//! Arithmetic uses the binary [`ExactNum`] kernel; the radix affects parse/format only.

use crate::defs::Error;
use crate::ext::ExactNum;
use crate::Consts;
use crate::Radix;
use crate::RoundingMode;
use core::ops::{Deref, DerefMut};

#[cfg(not(feature = "std"))]
use alloc::string::String;

/// A floating-point value with an associated parse/format radix.
#[derive(Debug, Clone)]
pub struct RadixFloat {
    value: ExactNum,
    radix: Radix,
}

impl RadixFloat {
    /// Wraps `value` with radix `radix`.
    ///
    /// ## Errors
    ///
    ///  - InvalidArgument: `radix` is outside 2..=36.
    pub fn new(value: ExactNum, radix: Radix) -> Result<Self, Error> {
        let _ = Radix::try_new(radix.value())?;
        Ok(Self { value, radix })
    }

    /// Wraps `value` with `radix` without validating the radix (for const radix constants).
    pub fn with_radix(value: ExactNum, radix: Radix) -> Self {
        Self { value, radix }
    }

    /// Returns the numeric value.
    pub fn value(&self) -> &ExactNum {
        &self.value
    }

    /// Returns the associated radix.
    pub const fn radix(&self) -> Radix {
        self.radix
    }

    /// Consumes `self` and returns the inner [`ExactNum`].
    pub fn into_inner(self) -> ExactNum {
        self.value
    }

    /// Parses `s` in `radix`.
    pub fn parse(
        s: &str,
        radix: Radix,
        p: usize,
        rm: RoundingMode,
        cc: &mut Consts,
    ) -> Result<Self, Error> {
        Radix::try_new(radix.value())?;
        Ok(Self {
            value: ExactNum::parse(s, radix, p, rm, cc),
            radix,
        })
    }

    /// Formats using the associated radix.
    pub fn format(&self, rm: RoundingMode, cc: &mut Consts) -> Result<String, Error> {
        self.value.format(self.radix, rm, cc)
    }
}

impl Deref for RadixFloat {
    type Target = ExactNum;

    fn deref(&self) -> &ExactNum {
        &self.value
    }
}

impl DerefMut for RadixFloat {
    fn deref_mut(&mut self) -> &mut ExactNum {
        &mut self.value
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_radix_float_roundtrip() {
        let mut cc = Consts::new().unwrap();
        let rdx = Radix::try_new(12).unwrap();
        let p = 128;
        let rm = RoundingMode::ToEven;
        let n = ExactNum::parse("10.5", Radix::Dec, p, rm, &mut cc);
        let rf = RadixFloat::with_radix(n.clone(), rdx);
        let s = rf.format(rm, &mut cc).unwrap();
        let g = RadixFloat::parse(&s, rdx, p, rm, &mut cc).unwrap();
        assert_eq!(n.cmp(g.value()), Some(0));
    }
}
