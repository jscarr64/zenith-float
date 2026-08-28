mod e;
mod ln10;
mod ln2;
mod pi;

use crate::common::buf::WordBuf;
use crate::common::util::round_p;
use crate::mantissa::Mantissa;
use crate::num::ExactNumNumber;
use crate::ops::consts::e::ECache;
use crate::ops::consts::ln10::Ln10Cache;
use crate::ops::consts::ln2::Ln2Cache;
use crate::ops::consts::pi::PiCache;
use crate::Error;
use crate::ExactNum;
use crate::RoundingMode;
use crate::WORD_BIT_SIZE;

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

/// Alias for [`Consts`]: a progressive cache of π, e, ln 2, ln 10, √2, and φ.
pub type ConstCache = Consts;

/// Snapshot of how many mantissa bits of each constant are currently cached.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ConstCacheInfo {
    /// Cached bits of π.
    pub pi: usize,
    /// Cached bits of e.
    pub e: usize,
    /// Cached bits of ln 2.
    pub ln2: usize,
    /// Cached bits of ln 10.
    pub ln10: usize,
    /// Cached bits of √2.
    pub sqrt2: usize,
    /// Cached bits of φ = (1+√5)/2.
    pub phi: usize,
}

/// A float stored at extra working precision so later requests at lower (or equal)
/// precision reuse the cache instead of recomputing.
#[derive(Clone, Debug)]
pub struct CachedFBig {
    inner: ExactNum,
}

impl CachedFBig {
    /// Wrap an already-computed value, retaining its current mantissa width.
    pub fn new(inner: ExactNum) -> Self {
        CachedFBig { inner }
    }

    /// Cached mantissa width in bits (`None` for Inf / NaN).
    pub fn cached_bit_len(&self) -> Option<usize> {
        self.inner.mantissa_max_bit_len()
    }

    /// The stored value.
    pub fn inner(&self) -> &ExactNum {
        &self.inner
    }

    /// Round the cached value to `p` bits.
    pub fn round(&self, p: usize, rm: RoundingMode) -> ExactNum {
        let mut v = self.inner.clone();
        let _ = v.set_precision(p, rm);
        v
    }
}

#[derive(Debug)]
struct ExtraCache {
    bits: usize,
    val: Option<ExactNumNumber>,
}

impl ExtraCache {
    fn new() -> Self {
        ExtraCache {
            bits: 0,
            val: None,
        }
    }

    fn cached_bit_len(&self) -> usize {
        self.bits
    }

    fn for_prec<F>(
        &mut self,
        k: usize,
        rm: RoundingMode,
        mut compute: F,
    ) -> Result<ExactNumNumber, Error>
    where
        F: FnMut(usize) -> Result<ExactNumNumber, Error>,
    {
        let p = round_p(k);
        let p_wrk = p
            .checked_add(WORD_BIT_SIZE)
            .ok_or(Error::InvalidArgument)?;
        if self.bits >= p {
            if let Some(v) = &self.val {
                let mut ret = v.clone()?;
                ret.set_precision(p, rm)?;
                return Ok(ret);
            }
        }
        let computed = compute(p_wrk)?;
        self.bits = computed.mantissa_max_bit_len();
        self.val = Some(computed.clone()?);
        let mut ret = computed;
        ret.set_precision(p, rm)?;
        Ok(ret)
    }
}

/// Constants cache contains arbitrary-precision mathematical constants.
#[derive(Debug)]
pub struct Consts {
    pi: PiCache,
    e: ECache,
    ln2: Ln2Cache,
    ln10: Ln10Cache,
    sqrt2: ExtraCache,
    phi: ExtraCache,
    tenpowers: Vec<(WordBuf, WordBuf, usize)>,
}

/// In an ideal situation, the `Consts` structure is initialized with `Consts::new` only once,
/// and then used where needed.
impl Consts {
    /// Initializes the constants cache.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    pub fn new() -> Result<Self, Error> {
        Ok(Consts {
            pi: PiCache::new()?,
            e: ECache::new()?,
            ln2: Ln2Cache::new()?,
            ln10: Ln10Cache::new()?,
            sqrt2: ExtraCache::new(),
            phi: ExtraCache::new(),
            tenpowers: Vec::new(),
        })
    }

    /// Returns the value of the pi number with precision `p` using rounding mode `rm`.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - InvalidArgument: the precision is incorrect.
    pub(crate) fn pi_num(&mut self, p: usize, rm: RoundingMode) -> Result<ExactNumNumber, Error> {
        let p = round_p(p);
        self.pi.for_prec(p, rm)
    }

    /// Returns the value of the Euler number with precision `p` using rounding mode `rm`.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - InvalidArgument: the precision is incorrect.
    pub(crate) fn e_num(&mut self, p: usize, rm: RoundingMode) -> Result<ExactNumNumber, Error> {
        let p = round_p(p);
        self.e.for_prec(p, rm)
    }

    /// Returns the value of the natural logarithm of 2 with precision `p` using rounding mode `rm`.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - InvalidArgument: the precision is incorrect.
    pub(crate) fn ln_2_num(&mut self, p: usize, rm: RoundingMode) -> Result<ExactNumNumber, Error> {
        let p = round_p(p);
        self.ln2.for_prec(p, rm)
    }

    /// Returns the value of the natural logarithm of 10 with precision `p` using rounding mode `rm`.
    /// Precision is rounded upwards to the word size.
    ///
    /// ## Errors
    ///
    ///  - MemoryAllocation: failed to allocate memory for mantissa.
    ///  - InvalidArgument: the precision is incorrect.
    pub(crate) fn ln_10_num(
        &mut self,
        p: usize,
        rm: RoundingMode,
    ) -> Result<ExactNumNumber, Error> {
        let p = round_p(p);
        self.ln10.for_prec(p, rm)
    }

    /// Returns the value of the pi number with precision `p` using rounding mode `rm`.
    /// Precision is rounded upwards to the word size.
    pub fn pi(&mut self, p: usize, rm: RoundingMode) -> ExactNum {
        match self.pi_num(p, rm) {
            Ok(v) => v.into(),
            Err(e) => ExactNum::nan(Some(e)),
        }
    }

    /// Returns the value of the Euler number with precision `p` using rounding mode `rm`.
    /// Precision is rounded upwards to the word size.
    pub fn e(&mut self, p: usize, rm: RoundingMode) -> ExactNum {
        match self.e_num(p, rm) {
            Ok(v) => v.into(),
            Err(e) => ExactNum::nan(Some(e)),
        }
    }

    /// Returns the value of the natural logarithm of 2 with precision `p` using rounding mode `rm`.
    /// Precision is rounded upwards to the word size.
    pub fn ln_2(&mut self, p: usize, rm: RoundingMode) -> ExactNum {
        match self.ln_2_num(p, rm) {
            Ok(v) => v.into(),
            Err(e) => ExactNum::nan(Some(e)),
        }
    }

    /// Returns the value of the natural logarithm of 10 with precision `p` using rounding mode `rm`.
    /// Precision is rounded upwards to the word size.
    pub fn ln_10(&mut self, p: usize, rm: RoundingMode) -> ExactNum {
        match self.ln_10_num(p, rm) {
            Ok(v) => v.into(),
            Err(e) => ExactNum::nan(Some(e)),
        }
    }

    /// Return powers of 10: 100, 10000, 100000000, ...
    pub(crate) fn tenpowers(&mut self, p: usize) -> Result<&[(WordBuf, WordBuf, usize)], Error> {
        if p >= self.tenpowers.len() {
            Mantissa::compute_tenpowers(&mut self.tenpowers, p)?;
        }

        Ok(&self.tenpowers)
    }

    /// How many bits of each series / extra constant are currently retained.
    pub fn cache_info(&self) -> ConstCacheInfo {
        ConstCacheInfo {
            pi: self.pi.cached_bit_len(),
            e: self.e.cached_bit_len(),
            ln2: self.ln2.cached_bit_len(),
            ln10: self.ln10.cached_bit_len(),
            sqrt2: self.sqrt2.cached_bit_len(),
            phi: self.phi.cached_bit_len(),
        }
    }

    fn sqrt2_num(&mut self, p: usize, rm: RoundingMode) -> Result<ExactNumNumber, Error> {
        self.sqrt2.for_prec(p, rm, |p_wrk| {
            let two = ExactNumNumber::from_word(2, p_wrk)?;
            two.sqrt(p_wrk, RoundingMode::None)
        })
    }

    fn phi_num(&mut self, p: usize, rm: RoundingMode) -> Result<ExactNumNumber, Error> {
        self.phi.for_prec(p, rm, |p_wrk| {
            let five = ExactNumNumber::from_word(5, p_wrk)?;
            let one = ExactNumNumber::from_word(1, p_wrk)?;
            let two = ExactNumNumber::from_word(2, p_wrk)?;
            let s = five.sqrt(p_wrk, RoundingMode::None)?;
            let n = one.add(&s, p_wrk, RoundingMode::None)?;
            n.div(&two, p_wrk, RoundingMode::None)
        })
    }

    /// √2 with precision `p` using rounding mode `rm`.
    /// Higher requests extend the cache; lower requests reuse it.
    pub fn sqrt2(&mut self, p: usize, rm: RoundingMode) -> ExactNum {
        match self.sqrt2_num(p, rm) {
            Ok(v) => v.into(),
            Err(e) => ExactNum::nan(Some(e)),
        }
    }

    /// Golden ratio φ = (1+√5)/2 with precision `p` using rounding mode `rm`.
    /// Higher requests extend the cache; lower requests reuse it.
    pub fn phi(&mut self, p: usize, rm: RoundingMode) -> ExactNum {
        match self.phi_num(p, rm) {
            Ok(v) => v.into(),
            Err(e) => ExactNum::nan(Some(e)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn progressive_constant_cache_extends() {
        let mut cc = Consts::new().expect("constants");
        let rm = RoundingMode::ToEven;
        let _ = cc.pi(128, rm);
        let _ = cc.e(128, rm);
        let _ = cc.ln_2(128, rm);
        let _ = cc.ln_10(128, rm);
        let _ = cc.sqrt2(128, rm);
        let _ = cc.phi(128, rm);
        let before = cc.cache_info();
        assert!(before.sqrt2 >= 128);
        assert!(before.phi >= 128);

        let _ = cc.pi(256, rm);
        let _ = cc.e(256, rm);
        let _ = cc.ln_2(256, rm);
        let _ = cc.ln_10(256, rm);
        let _ = cc.sqrt2(256, rm);
        let _ = cc.phi(256, rm);
        let after = cc.cache_info();
        assert!(after.pi >= before.pi);
        assert!(after.e >= before.e);
        assert!(after.ln2 >= before.ln2);
        assert!(after.ln10 >= before.ln10);
        assert!(after.sqrt2 >= 256);
        assert!(after.phi >= 256);

        let a = cc.sqrt2(128, rm);
        let b = cc.sqrt2(128, rm);
        assert_eq!(a.cmp(&b), Some(0));
        let cached = CachedFBig::new(a);
        assert!(cached.cached_bit_len().unwrap() >= 128);
        let r = cached.round(64, rm);
        assert!(!r.is_nan());
    }
}
