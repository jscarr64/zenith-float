//! Compact binary interchange for [`ExactNum`] and [`ExactNumArray`].
//!
//! The 16-byte inline record is target-independent: two big-endian `u32`
//! mantissa words (64 bits), matching `INLINE_WORDS = 2` on a 32-bit `Word`.
//! Wider mantissas use a heap record with the same `u32` limb unit so a
//! 32-bit reader and a 64-bit writer agree. All multi-byte fields are
//! network byte order.

use crate::ieee_soft::ExactNumArray;
use crate::ExactNum;
use crate::Exponent;
use crate::Sign;
use crate::Word;
use crate::INF_NEG;
use crate::INF_POS;
use crate::WORD_BIT_SIZE;
use crate::{Error, NAN};
use alloc::vec::Vec;

/// First-byte version stored in every record.
pub const BINARY_FORMAT_VERSION: u8 = 1;

/// Bytes in [`InlineBinaryBuffer`] / [`ExactNum::to_inline_bytes`].
pub const BINARY_INLINE_LEN: usize = 16;

/// Mantissa bits that fit in the inline record (`2 × 32`).
pub const BINARY_INLINE_MANT_BITS: usize = 64;

/// `u32` limbs in the inline mantissa field.
pub const BINARY_INLINE_U32_WORDS: usize = 2;

/// Header bytes of a heap or array record (same size as the inline record).
pub const BINARY_HEADER_LEN: usize = 16;

/// Maximum `u32` limbs accepted by [`ExactNum::from_bytes`].
pub const BINARY_MAX_U32: usize = 65_536;

/// Maximum array elements accepted by [`ExactNumArray::from_bytes`].
pub const BINARY_MAX_ELEMS: usize = 1_048_576;

/// Positive finite (inline).
pub const BINARY_FLAG_POS: u8 = 0x00;
/// Negative finite (inline).
pub const BINARY_FLAG_NEG: u8 = 0x01;
/// [`INF_POS`].
pub const BINARY_FLAG_INF_POS: u8 = 0x02;
/// [`INF_NEG`].
pub const BINARY_FLAG_INF_NEG: u8 = 0x03;
/// NaN: [`Error::DivisionByZero`].
pub const BINARY_FLAG_NAN_DIV0: u8 = 0x04;
/// NaN: [`Error::InvalidArgument`].
pub const BINARY_FLAG_NAN_INVALID: u8 = 0x05;
/// NaN: [`Error::PrecisionRetryExhausted`].
pub const BINARY_FLAG_NAN_RETRY: u8 = 0x06;
/// NaN: [`Error::MemoryAllocation`].
pub const BINARY_FLAG_NAN_OOM: u8 = 0x07;
/// NaN: [`Error::ExponentOverflow`] with [`Sign::Pos`].
pub const BINARY_FLAG_NAN_OVF_POS: u8 = 0x08;
/// NaN: [`Error::ExponentOverflow`] with [`Sign::Neg`].
pub const BINARY_FLAG_NAN_OVF_NEG: u8 = 0x09;
/// NaN with no associated [`Error`].
pub const BINARY_FLAG_NAN_BARE: u8 = 0x0A;
/// Heap finite, positive.
pub const BINARY_FLAG_HEAP_POS: u8 = 0x10;
/// Heap finite, negative.
pub const BINARY_FLAG_HEAP_NEG: u8 = 0x11;
/// [`ExactNumArray`] record.
pub const BINARY_FLAG_ARRAY: u8 = 0x20;

const U32_BYTES: usize = 4;
const WORD_U32: usize = WORD_BIT_SIZE / 32;

/// Fixed 16-byte stack buffer for an inlined [`ExactNum`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct InlineBinaryBuffer {
    bytes: [u8; BINARY_INLINE_LEN],
}

impl InlineBinaryBuffer {
    /// The 16 encoded bytes.
    pub fn as_bytes(&self) -> &[u8; BINARY_INLINE_LEN] {
        &self.bytes
    }

    /// Consume and return the 16 encoded bytes.
    pub fn into_bytes(self) -> [u8; BINARY_INLINE_LEN] {
        self.bytes
    }
}

impl ExactNum {
    /// Encode an inlined value into a stack buffer.
    ///
    /// Specials always succeed. A finite mantissa wider than
    /// [`BINARY_INLINE_MANT_BITS`] returns [`Error::MemoryAllocation`].
    pub fn to_inline_bytes(&self) -> Result<InlineBinaryBuffer, Error> {
        let mut bytes = [0u8; BINARY_INLINE_LEN];
        write_inline(self, &mut bytes)?;
        Ok(InlineBinaryBuffer { bytes })
    }

    /// Write the 16-byte inline record into `dest`.
    ///
    /// Returns [`Error::InvalidArgument`] if `dest` is shorter than
    /// [`BINARY_INLINE_LEN`]. Wider finite mantissas return
    /// [`Error::MemoryAllocation`].
    pub fn write_inline_bytes(&self, dest: &mut [u8]) -> Result<usize, Error> {
        if dest.len() < BINARY_INLINE_LEN {
            return Err(Error::InvalidArgument);
        }
        write_inline(self, dest)?;
        Ok(BINARY_INLINE_LEN)
    }

    /// Decode a 16-byte inline record.
    pub fn from_inline_bytes(bytes: &[u8; BINARY_INLINE_LEN]) -> Result<Self, Error> {
        let (v, n) = decode_num(bytes)?;
        if n != BINARY_INLINE_LEN {
            return Err(Error::InvalidArgument);
        }
        Ok(v)
    }

    /// Encode `self` (inline 16 bytes, or a heap record if the mantissa is wider).
    pub fn to_bytes(&self) -> Result<Vec<u8>, Error> {
        let n = encoded_len(self)?;
        let mut out = Vec::new();
        out.try_reserve_exact(n)?;
        out.resize(n, 0);
        let wrote = write_num(self, &mut out)?;
        if wrote != n {
            return Err(Error::InvalidArgument);
        }
        Ok(out)
    }

    /// Write the compact record into `dest` without allocating.
    pub fn write_bytes(&self, dest: &mut [u8]) -> Result<usize, Error> {
        write_num(self, dest)
    }

    /// Decode a compact record. Extra trailing bytes are [`Error::InvalidArgument`].
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
        let (v, n) = decode_num(bytes)?;
        if n != bytes.len() {
            return Err(Error::InvalidArgument);
        }
        Ok(v)
    }
}

impl ExactNumArray {
    /// Encode shape, shared `p`, and every element.
    pub fn to_bytes(&self) -> Result<Vec<u8>, Error> {
        let n = encoded_array_len(self)?;
        let mut out = Vec::new();
        out.try_reserve_exact(n)?;
        out.resize(n, 0);
        let wrote = write_array(self, &mut out)?;
        if wrote != n {
            return Err(Error::InvalidArgument);
        }
        Ok(out)
    }

    /// Write the array record into `dest` without allocating.
    pub fn write_bytes(&self, dest: &mut [u8]) -> Result<usize, Error> {
        write_array(self, dest)
    }

    /// Decode an array record. Shape product must match the element count.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
        decode_array(bytes)
    }
}

fn flag_nan(err: Option<Error>) -> u8 {
    match err {
        Some(Error::DivisionByZero) => BINARY_FLAG_NAN_DIV0,
        Some(Error::InvalidArgument) => BINARY_FLAG_NAN_INVALID,
        Some(Error::PrecisionRetryExhausted) => BINARY_FLAG_NAN_RETRY,
        Some(Error::MemoryAllocation) => BINARY_FLAG_NAN_OOM,
        Some(Error::ExponentOverflow(Sign::Pos)) => BINARY_FLAG_NAN_OVF_POS,
        Some(Error::ExponentOverflow(Sign::Neg)) => BINARY_FLAG_NAN_OVF_NEG,
        None => BINARY_FLAG_NAN_BARE,
    }
}

fn nan_from_flag(flag: u8) -> Result<ExactNum, Error> {
    match flag {
        BINARY_FLAG_NAN_DIV0 => Ok(ExactNum::nan(Some(Error::DivisionByZero))),
        BINARY_FLAG_NAN_INVALID => Ok(ExactNum::nan(Some(Error::InvalidArgument))),
        BINARY_FLAG_NAN_RETRY => Ok(ExactNum::nan(Some(Error::PrecisionRetryExhausted))),
        BINARY_FLAG_NAN_OOM => Ok(ExactNum::nan(Some(Error::MemoryAllocation))),
        BINARY_FLAG_NAN_OVF_POS => Ok(ExactNum::nan(Some(Error::ExponentOverflow(Sign::Pos)))),
        BINARY_FLAG_NAN_OVF_NEG => Ok(ExactNum::nan(Some(Error::ExponentOverflow(Sign::Neg)))),
        BINARY_FLAG_NAN_BARE => Ok(NAN.clone()),
        _ => Err(Error::InvalidArgument),
    }
}

fn is_special_flag(flag: u8) -> bool {
    matches!(
        flag,
        BINARY_FLAG_INF_POS
            | BINARY_FLAG_INF_NEG
            | BINARY_FLAG_NAN_DIV0
            | BINARY_FLAG_NAN_INVALID
            | BINARY_FLAG_NAN_RETRY
            | BINARY_FLAG_NAN_OOM
            | BINARY_FLAG_NAN_OVF_POS
            | BINARY_FLAG_NAN_OVF_NEG
            | BINARY_FLAG_NAN_BARE
    )
}

fn write_special(flag: u8, dest: &mut [u8]) -> Result<usize, Error> {
    if dest.len() < BINARY_INLINE_LEN {
        return Err(Error::InvalidArgument);
    }
    dest[..BINARY_INLINE_LEN].fill(0);
    dest[0] = flag;
    dest[1] = BINARY_FORMAT_VERSION;
    Ok(BINARY_INLINE_LEN)
}

fn write_inline(x: &ExactNum, dest: &mut [u8]) -> Result<usize, Error> {
    if dest.len() < BINARY_INLINE_LEN {
        return Err(Error::InvalidArgument);
    }
    if x.is_nan() {
        return write_special(flag_nan(x.err()), dest);
    }
    if x.is_inf_pos() {
        return write_special(BINARY_FLAG_INF_POS, dest);
    }
    if x.is_inf_neg() {
        return write_special(BINARY_FLAG_INF_NEG, dest);
    }
    let Some((m, n_sig, sign, exp, inexact)) = x.as_raw_parts() else {
        return Err(Error::InvalidArgument);
    };
    let n_u32 = u32_count(m.len());
    if n_u32 > BINARY_INLINE_U32_WORDS {
        return Err(Error::MemoryAllocation);
    }
    dest[..BINARY_INLINE_LEN].fill(0);
    dest[0] = if sign == Sign::Neg { BINARY_FLAG_NEG } else { BINARY_FLAG_POS };
    dest[1] = BINARY_FORMAT_VERSION;
    dest[2] = u8::try_from(n_sig).map_err(|_| Error::InvalidArgument)?;
    dest[3] = u8::from(inexact);
    dest[4..8].copy_from_slice(&exp.to_be_bytes());
    let mut tmp = [0u32; BINARY_INLINE_U32_WORDS];
    words_to_u32(m, &mut tmp)?;
    put_u32_be(&mut dest[8..12], tmp[0]);
    put_u32_be(&mut dest[12..16], tmp[1]);
    Ok(BINARY_INLINE_LEN)
}

fn encoded_len(x: &ExactNum) -> Result<usize, Error> {
    if x.is_nan() || x.is_inf() {
        return Ok(BINARY_INLINE_LEN);
    }
    let Some((m, _, _, _, _)) = x.as_raw_parts() else {
        return Err(Error::InvalidArgument);
    };
    let n_u32 = u32_count(m.len());
    if n_u32 <= BINARY_INLINE_U32_WORDS {
        Ok(BINARY_INLINE_LEN)
    } else {
        heap_len(n_u32)
    }
}

fn heap_len(n_u32: usize) -> Result<usize, Error> {
    if n_u32 > BINARY_MAX_U32 {
        return Err(Error::InvalidArgument);
    }
    n_u32
        .checked_mul(U32_BYTES)
        .and_then(|b| b.checked_add(BINARY_HEADER_LEN))
        .ok_or(Error::InvalidArgument)
}

fn write_num(x: &ExactNum, dest: &mut [u8]) -> Result<usize, Error> {
    if x.is_nan() || x.is_inf() {
        return write_inline(x, dest);
    }
    let Some((m, n_sig, sign, exp, inexact)) = x.as_raw_parts() else {
        return Err(Error::InvalidArgument);
    };
    let n_u32 = u32_count(m.len());
    if n_u32 <= BINARY_INLINE_U32_WORDS {
        return write_inline(x, dest);
    }
    let need = heap_len(n_u32)?;
    if dest.len() < need {
        return Err(Error::InvalidArgument);
    }
    dest[..need].fill(0);
    dest[0] = if sign == Sign::Neg { BINARY_FLAG_HEAP_NEG } else { BINARY_FLAG_HEAP_POS };
    dest[1] = BINARY_FORMAT_VERSION;
    dest[3] = u8::from(inexact);
    dest[4..8].copy_from_slice(&exp.to_be_bytes());
    put_u32_be(&mut dest[8..12], u32_from_usize(n_sig)?);
    put_u32_be(&mut dest[12..16], u32_from_usize(n_u32)?);
    let mut limb = [0u32; 2];
    let mut off = BINARY_HEADER_LEN;
    for &w in m {
        words_to_u32(core::slice::from_ref(&w), &mut limb[..WORD_U32])?;
        for &u in limb.iter().take(WORD_U32) {
            put_u32_be(&mut dest[off..off + U32_BYTES], u);
            off += U32_BYTES;
        }
    }
    Ok(need)
}

fn decode_num(bytes: &[u8]) -> Result<(ExactNum, usize), Error> {
    if bytes.len() < 2 {
        return Err(Error::InvalidArgument);
    }
    let flag = bytes[0];
    if bytes[1] != BINARY_FORMAT_VERSION {
        return Err(Error::InvalidArgument);
    }
    if flag == BINARY_FLAG_ARRAY {
        return Err(Error::InvalidArgument);
    }
    if is_special_flag(flag) {
        if bytes.len() < BINARY_INLINE_LEN {
            return Err(Error::InvalidArgument);
        }
        let v = match flag {
            BINARY_FLAG_INF_POS => INF_POS.clone(),
            BINARY_FLAG_INF_NEG => INF_NEG.clone(),
            _ => nan_from_flag(flag)?,
        };
        return Ok((v, BINARY_INLINE_LEN));
    }
    match flag {
        BINARY_FLAG_POS | BINARY_FLAG_NEG => decode_inline_finite(bytes, flag),
        BINARY_FLAG_HEAP_POS | BINARY_FLAG_HEAP_NEG => decode_heap(bytes, flag),
        _ => Err(Error::InvalidArgument),
    }
}

fn decode_inline_finite(bytes: &[u8], flag: u8) -> Result<(ExactNum, usize), Error> {
    if bytes.len() < BINARY_INLINE_LEN {
        return Err(Error::InvalidArgument);
    }
    let n_sig = bytes[2] as usize;
    let inexact = bytes[3] != 0;
    let exp = i32::from_be_bytes(bytes[4..8].try_into().map_err(|_| Error::InvalidArgument)?);
    let u0 = u32::from_be_bytes(
        bytes[8..12]
            .try_into()
            .map_err(|_| Error::InvalidArgument)?,
    );
    let u1 = u32::from_be_bytes(
        bytes[12..16]
            .try_into()
            .map_err(|_| Error::InvalidArgument)?,
    );
    let sign = if flag == BINARY_FLAG_NEG { Sign::Neg } else { Sign::Pos };
    let x = finite_from_u32(&[u0, u1], n_sig, sign, exp, inexact)?;
    Ok((x, BINARY_INLINE_LEN))
}

fn decode_heap(bytes: &[u8], flag: u8) -> Result<(ExactNum, usize), Error> {
    if bytes.len() < BINARY_HEADER_LEN {
        return Err(Error::InvalidArgument);
    }
    let inexact = bytes[3] != 0;
    let exp = i32::from_be_bytes(bytes[4..8].try_into().map_err(|_| Error::InvalidArgument)?);
    let n_sig = u32::from_be_bytes(
        bytes[8..12]
            .try_into()
            .map_err(|_| Error::InvalidArgument)?,
    ) as usize;
    let n_u32 = u32::from_be_bytes(
        bytes[12..16]
            .try_into()
            .map_err(|_| Error::InvalidArgument)?,
    ) as usize;
    if n_u32 == 0 || n_u32 > BINARY_MAX_U32 {
        return Err(Error::InvalidArgument);
    }
    let need = heap_len(n_u32)?;
    if bytes.len() < need {
        return Err(Error::InvalidArgument);
    }
    let mut u32s = Vec::new();
    u32s.try_reserve_exact(n_u32)?;
    let mut off = BINARY_HEADER_LEN;
    for _ in 0..n_u32 {
        let chunk: [u8; 4] = bytes[off..off + U32_BYTES]
            .try_into()
            .map_err(|_| Error::InvalidArgument)?;
        u32s.push(u32::from_be_bytes(chunk));
        off += U32_BYTES;
    }
    let sign = if flag == BINARY_FLAG_HEAP_NEG { Sign::Neg } else { Sign::Pos };
    let x = finite_from_u32(&u32s, n_sig, sign, exp, inexact)?;
    Ok((x, need))
}

fn finite_from_u32(
    u32s: &[u32],
    n_sig: usize,
    sign: Sign,
    exp: Exponent,
    inexact: bool,
) -> Result<ExactNum, Error> {
    if u32s.len() % WORD_U32 != 0 {
        return Err(Error::InvalidArgument);
    }
    if n_sig > u32s.len() * 32 {
        return Err(Error::InvalidArgument);
    }
    let n_words = u32s.len() / WORD_U32;
    let mut words = Vec::new();
    words.try_reserve_exact(n_words)?;
    for i in 0..n_words {
        words.push(word_from_u32s(&u32s[i * WORD_U32..(i + 1) * WORD_U32]));
    }
    let x = ExactNum::from_raw_parts(&words, n_sig, sign, exp, inexact);
    if x.is_nan() {
        return Err(x.err().unwrap_or(Error::InvalidArgument));
    }
    Ok(x)
}

fn encoded_array_len(a: &ExactNumArray) -> Result<usize, Error> {
    if a.len() > BINARY_MAX_ELEMS {
        return Err(Error::InvalidArgument);
    }
    let mut n = BINARY_HEADER_LEN;
    for v in a.as_slice() {
        n = n
            .checked_add(encoded_len(v)?)
            .ok_or(Error::InvalidArgument)?;
    }
    Ok(n)
}

fn write_array(a: &ExactNumArray, dest: &mut [u8]) -> Result<usize, Error> {
    let need = encoded_array_len(a)?;
    if dest.len() < need {
        return Err(Error::InvalidArgument);
    }
    let (rows, cols) = a.shape();
    dest[..BINARY_HEADER_LEN].fill(0);
    dest[0] = BINARY_FLAG_ARRAY;
    dest[1] = BINARY_FORMAT_VERSION;
    put_u32_be(&mut dest[4..8], u32_from_usize(rows)?);
    put_u32_be(&mut dest[8..12], u32_from_usize(cols)?);
    put_u32_be(&mut dest[12..16], u32_from_usize(a.precision())?);
    let mut off = BINARY_HEADER_LEN;
    for v in a.as_slice() {
        off += write_num(v, &mut dest[off..])?;
    }
    Ok(need)
}

fn decode_array(bytes: &[u8]) -> Result<ExactNumArray, Error> {
    if bytes.len() < BINARY_HEADER_LEN {
        return Err(Error::InvalidArgument);
    }
    if bytes[0] != BINARY_FLAG_ARRAY || bytes[1] != BINARY_FORMAT_VERSION {
        return Err(Error::InvalidArgument);
    }
    let rows =
        u32::from_be_bytes(bytes[4..8].try_into().map_err(|_| Error::InvalidArgument)?) as usize;
    let cols = u32::from_be_bytes(
        bytes[8..12]
            .try_into()
            .map_err(|_| Error::InvalidArgument)?,
    ) as usize;
    let p = u32::from_be_bytes(
        bytes[12..16]
            .try_into()
            .map_err(|_| Error::InvalidArgument)?,
    ) as usize;
    let n = rows.checked_mul(cols).ok_or(Error::InvalidArgument)?;
    if n > BINARY_MAX_ELEMS {
        return Err(Error::InvalidArgument);
    }
    let mut vals = Vec::new();
    vals.try_reserve_exact(n)?;
    let mut off = BINARY_HEADER_LEN;
    for _ in 0..n {
        let (v, used) = decode_num(&bytes[off..])?;
        off += used;
        vals.push(v);
    }
    if off != bytes.len() {
        return Err(Error::InvalidArgument);
    }
    ExactNumArray::from_parts(p, rows, cols, vals)
}

fn u32_count(n_words: usize) -> usize {
    n_words * WORD_U32
}

fn words_to_u32(words: &[Word], out: &mut [u32]) -> Result<(), Error> {
    let need = u32_count(words.len());
    if out.len() < need {
        return Err(Error::InvalidArgument);
    }
    let mut i = 0;
    for &w in words {
        #[cfg(target_pointer_width = "32")]
        {
            out[i] = w;
            i += 1;
        }
        #[cfg(not(target_pointer_width = "32"))]
        {
            out[i] = w as u32;
            out[i + 1] = (w >> 32) as u32;
            i += 2;
        }
    }
    let _ = i;
    Ok(())
}

fn word_from_u32s(part: &[u32]) -> Word {
    #[cfg(target_pointer_width = "32")]
    {
        part[0]
    }
    #[cfg(not(target_pointer_width = "32"))]
    {
        (part[0] as Word) | ((part[1] as Word) << 32)
    }
}

fn put_u32_be(dest: &mut [u8], v: u32) {
    dest[..U32_BYTES].copy_from_slice(&v.to_be_bytes());
}

fn u32_from_usize(v: usize) -> Result<u32, Error> {
    u32::try_from(v).map_err(|_| Error::InvalidArgument)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::RoundingMode;

    fn p64() -> usize {
        64
    }

    fn raw_eq(a: &ExactNum, b: &ExactNum) -> bool {
        a.as_raw_parts() == b.as_raw_parts()
    }

    #[test]
    fn binary_inline_roundtrip_finite() {
        let p = p64();
        let vals = [
            ExactNum::from_u8(0, p),
            ExactNum::from_u8(1, p),
            ExactNum::from_u8(1, p).neg(),
            ExactNum::from_u32(u32::MAX, p),
            ExactNum::from_u8(1, p).div(&ExactNum::from_u8(3, p), p, RoundingMode::ToEven),
        ];
        for x in &vals {
            let buf = x.to_inline_bytes().unwrap();
            assert_eq!(buf.as_bytes().len(), BINARY_INLINE_LEN);
            assert_eq!(buf.as_bytes()[1], BINARY_FORMAT_VERSION);
            let y = ExactNum::from_inline_bytes(buf.as_bytes()).unwrap();
            assert!(raw_eq(x, &y), "inline limbs");
            assert_eq!(x.cmp(&y), Some(0));
            let v = x.to_bytes().unwrap();
            assert_eq!(v.len(), BINARY_INLINE_LEN);
            let z = ExactNum::from_bytes(&v).unwrap();
            assert!(raw_eq(x, &z));
        }
    }

    #[test]
    fn binary_nan_flag_is_recognizable() {
        let buf = NAN.to_inline_bytes().unwrap();
        assert_eq!(buf.as_bytes()[0], BINARY_FLAG_NAN_BARE);
        assert_eq!(buf.as_bytes()[1], BINARY_FORMAT_VERSION);
        let y = ExactNum::from_inline_bytes(buf.as_bytes()).unwrap();
        assert!(y.is_nan());
        let div = ExactNum::nan(Some(Error::DivisionByZero));
        let d = div.to_bytes().unwrap();
        assert_eq!(d[0], BINARY_FLAG_NAN_DIV0);
        let back = ExactNum::from_bytes(&d).unwrap();
        assert!(back.is_nan());
        assert_eq!(back.err(), Some(Error::DivisionByZero));
    }

    #[test]
    fn binary_inf_roundtrip() {
        let p = INF_POS.to_bytes().unwrap();
        let n = INF_NEG.to_bytes().unwrap();
        assert_eq!(p[0], BINARY_FLAG_INF_POS);
        assert_eq!(n[0], BINARY_FLAG_INF_NEG);
        assert!(ExactNum::from_bytes(&p).unwrap().is_inf_pos());
        assert!(ExactNum::from_bytes(&n).unwrap().is_inf_neg());
    }

    #[test]
    fn binary_heap_roundtrip_p256() {
        let p = 256;
        let x = ExactNum::from_u8(1, p).div(&ExactNum::from_u8(3, p), p, RoundingMode::ToEven);
        assert!(x.mantissa_max_bit_len().unwrap() > BINARY_INLINE_MANT_BITS);
        assert!(x.to_inline_bytes().is_err());
        let v = x.to_bytes().unwrap();
        assert!(v.len() > BINARY_INLINE_LEN);
        assert_eq!(v[0], BINARY_FLAG_HEAP_POS);
        let y = ExactNum::from_bytes(&v).unwrap();
        assert!(raw_eq(&x, &y));
        assert_eq!(x.cmp(&y), Some(0));
    }

    #[test]
    fn binary_array_roundtrip_and_shape_err() {
        let p = p64();
        let n = |k: u8| ExactNum::from_u8(k, p);
        let a = ExactNumArray::from_shape(p, 2, 2, &[n(1), n(2), n(3), n(4)]).unwrap();
        let v = a.to_bytes().unwrap();
        assert_eq!(v[0], BINARY_FLAG_ARRAY);
        let b = ExactNumArray::from_bytes(&v).unwrap();
        assert_eq!(b.shape(), (2, 2));
        assert_eq!(b.precision(), p);
        for i in 0..4 {
            assert_eq!(a.get(i).unwrap().cmp(b.get(i).unwrap()), Some(0));
            assert!(raw_eq(a.get(i).unwrap(), b.get(i).unwrap()));
        }
        let mut bad = v.clone();
        // rows = 3, cols = 2, but still four payloads → leftover or short read
        bad[4..8].copy_from_slice(&3u32.to_be_bytes());
        assert!(ExactNumArray::from_bytes(&bad).is_err());
    }

    #[test]
    fn binary_invalid_bytes_are_err() {
        assert!(ExactNum::from_bytes(&[]).is_err());
        assert!(ExactNum::from_bytes(&[0]).is_err());
        assert!(ExactNum::from_bytes(&[BINARY_FLAG_POS, 0]).is_err());
        assert!(ExactNum::from_bytes(&[0x7F, BINARY_FORMAT_VERSION]).is_err());
        let mut short = NAN.to_bytes().unwrap();
        short.pop();
        assert!(ExactNum::from_bytes(&short).is_err());
        let mut extra = ExactNum::from_u8(1, p64()).to_bytes().unwrap();
        extra.push(0);
        assert!(ExactNum::from_bytes(&extra).is_err());
        assert!(ExactNumArray::from_bytes(&[BINARY_FLAG_ARRAY, BINARY_FORMAT_VERSION]).is_err());
    }

    #[test]
    fn binary_write_bytes_zero_alloc_inline() {
        let x = ExactNum::from_u8(2, p64());
        let mut dest = [0u8; BINARY_INLINE_LEN];
        let n = x.write_inline_bytes(&mut dest).unwrap();
        assert_eq!(n, BINARY_INLINE_LEN);
        assert_eq!(dest, x.to_inline_bytes().unwrap().into_bytes());
        assert!(x.write_inline_bytes(&mut dest[..8]).is_err());
    }

    #[test]
    fn u32_word_boundary_add_is_2_pow_32() {
        let p = p64();
        let a = ExactNum::from_u32(u32::MAX, p);
        let b = ExactNum::from_u8(1, p);
        let s = a.add(&b, p, RoundingMode::ToEven);
        let expect = ExactNum::from_u64(1u64 << 32, p);
        assert_eq!(s.cmp(&expect), Some(0));
        assert!(s.err().is_none());
        let buf = s.to_inline_bytes().unwrap();
        assert!(raw_eq(
            &s,
            &ExactNum::from_inline_bytes(buf.as_bytes()).unwrap()
        ));
    }
}
