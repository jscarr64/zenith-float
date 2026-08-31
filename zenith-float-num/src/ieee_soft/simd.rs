//! Integer SIMD for software IEEE arrays. Lanes are `u32` / `u64` bit patterns.
//! Hardware IEEE arithmetic is not used.

/// Number of `u32` lanes in one integer SIMD vector (128-bit register).
/// Binary64 uses `IEEE_SIMD_LANE_WIDTH / 2` `u64` lanes. Scalar fallback
/// uses the same width so wrappers can size buffers without `cfg`.
pub const IEEE_SIMD_LANE_WIDTH: usize = 4;

use super::arith::{
    add_bits, div_bits, fma_add, fma_bits, isqrt, mul_bits, normalize_mul, pack_finite, pack_mag,
    round_rne, sqrt_bits, sub_bits, unpack, Class, Format, BIN32, BIN64,
};

/// Four binary32 bit patterns.
pub type Bin32x4 = [u32; 4];
/// Two binary64 bit patterns.
pub type Bin64x2 = [u64; 2];

#[derive(Clone, Copy)]
struct U32x4([u32; 4]);

impl U32x4 {
    fn load(v: [u32; 4]) -> Self {
        Self(v)
    }

    fn splat(x: u32) -> Self {
        Self([x; 4])
    }

    fn to_array(self) -> [u32; 4] {
        self.0
    }

    fn and(self, o: Self) -> Self {
        self.bop(o, core::ops::BitAnd::bitand)
    }

    fn or(self, o: Self) -> Self {
        self.bop(o, core::ops::BitOr::bitor)
    }

    fn xor(self, o: Self) -> Self {
        self.bop(o, core::ops::BitXor::bitxor)
    }

    fn wrapping_add(self, o: Self) -> Self {
        #[cfg(target_arch = "x86_64")]
        unsafe {
            use core::arch::x86_64::{__m128i, _mm_add_epi32, _mm_loadu_si128, _mm_storeu_si128};
            let a = _mm_loadu_si128(self.0.as_ptr() as *const __m128i);
            let b = _mm_loadu_si128(o.0.as_ptr() as *const __m128i);
            let s = _mm_add_epi32(a, b);
            let mut out = [0u32; 4];
            _mm_storeu_si128(out.as_mut_ptr() as *mut __m128i, s);
            return Self(out);
        }
        #[cfg(target_arch = "aarch64")]
        unsafe {
            use core::arch::aarch64::{uint32x4_t, vaddq_u32, vld1q_u32, vst1q_u32};
            let a = vld1q_u32(self.0.as_ptr());
            let b = vld1q_u32(o.0.as_ptr());
            let s: uint32x4_t = vaddq_u32(a, b);
            let mut out = [0u32; 4];
            vst1q_u32(out.as_mut_ptr(), s);
            return Self(out);
        }
        #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
        {
            self.bop(o, u32::wrapping_add)
        }
    }

    fn srli(self, n: u32) -> Self {
        Self([self.0[0] >> n, self.0[1] >> n, self.0[2] >> n, self.0[3] >> n])
    }

    fn eq_mask(self, o: Self) -> bool {
        self.0 == o.0
    }

    fn lane_eq(self, o: Self) -> [bool; 4] {
        [
            self.0[0] == o.0[0],
            self.0[1] == o.0[1],
            self.0[2] == o.0[2],
            self.0[3] == o.0[3],
        ]
    }

    fn bop(self, o: Self, f: fn(u32, u32) -> u32) -> Self {
        Self([
            f(self.0[0], o.0[0]),
            f(self.0[1], o.0[1]),
            f(self.0[2], o.0[2]),
            f(self.0[3], o.0[3]),
        ])
    }
}

fn all_normal_bin32(bits: U32x4) -> bool {
    let exp = bits.srli(23).and(U32x4::splat(0xff));
    !exp.lane_eq(U32x4::splat(0)).iter().any(|&z| z)
        && !exp.lane_eq(U32x4::splat(0xff)).iter().any(|&z| z)
}

fn mul_u32x4_integer(a: [u32; 4], b: [u32; 4]) -> [u64; 4] {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        use core::arch::x86_64::{
            __m128i, _mm_cvtsi128_si64, _mm_loadu_si128, _mm_mul_epu32, _mm_srli_si128,
        };
        let va = _mm_loadu_si128(a.as_ptr() as *const __m128i);
        let vb = _mm_loadu_si128(b.as_ptr() as *const __m128i);
        let p_even = _mm_mul_epu32(va, vb);
        let p_odd = _mm_mul_epu32(_mm_srli_si128(va, 4), _mm_srli_si128(vb, 4));
        return [
            _mm_cvtsi128_si64(p_even) as u64,
            _mm_cvtsi128_si64(p_odd) as u64,
            _mm_cvtsi128_si64(_mm_srli_si128(p_even, 8)) as u64,
            _mm_cvtsi128_si64(_mm_srli_si128(p_odd, 8)) as u64,
        ];
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        [
            a[0] as u64 * b[0] as u64,
            a[1] as u64 * b[1] as u64,
            a[2] as u64 * b[2] as u64,
            a[3] as u64 * b[3] as u64,
        ]
    }
}

/// Software IEEE add on four binary32 lanes. Integer SIMD on the equal-exp
/// same-sign normal path; otherwise the scalar integer kernel.
pub fn add_bin32_x4(a: Bin32x4, b: Bin32x4) -> Bin32x4 {
    let va = U32x4::load(a);
    let vb = U32x4::load(b);
    if all_normal_bin32(va) && all_normal_bin32(vb) {
        let sign_a = va.srli(31);
        let sign_b = vb.srli(31);
        let exp_a = va.srli(23).and(U32x4::splat(0xff));
        let exp_b = vb.srli(23).and(U32x4::splat(0xff));
        if sign_a.eq_mask(sign_b) && exp_a.eq_mask(exp_b) {
            let hidden = U32x4::splat(1 << 23);
            let frac = U32x4::splat(0x7f_ffff);
            let sa = va.and(frac).or(hidden);
            let sb = vb.and(frac).or(hidden);
            let sa3 = U32x4([sa.0[0] << 3, sa.0[1] << 3, sa.0[2] << 3, sa.0[3] << 3]);
            let sb3 = U32x4([sb.0[0] << 3, sb.0[1] << 3, sb.0[2] << 3, sb.0[3] << 3]);
            let sum = sa3.wrapping_add(sb3);
            let mut out = [0u32; 4];
            for i in 0..4 {
                let (exp, core) = round_rne(sum.0[i] as u128, exp_a.0[i] as i32, false, BIN32);
                out[i] = pack_finite(sign_a.0[i] != 0, exp, core, BIN32) as u32;
            }
            debug_assert_bit_eq32(&out, a, b, add_bits);
            return out;
        }
    }
    [
        add_bits(a[0] as u64, b[0] as u64, BIN32) as u32,
        add_bits(a[1] as u64, b[1] as u64, BIN32) as u32,
        add_bits(a[2] as u64, b[2] as u64, BIN32) as u32,
        add_bits(a[3] as u64, b[3] as u64, BIN32) as u32,
    ]
}

fn debug_assert_bit_eq32(got: &Bin32x4, a: Bin32x4, b: Bin32x4, op: fn(u64, u64, Format) -> u64) {
    for i in 0..4 {
        debug_assert_eq!(
            got[i],
            op(a[i] as u64, b[i] as u64, BIN32) as u32,
            "SIMD lane bits must match the scalar integer kernel"
        );
    }
}

/// Software IEEE mul on four binary32 lanes. Integer SIMD for significand products
/// when every lane is normal.
pub fn mul_bin32_x4(a: Bin32x4, b: Bin32x4) -> Bin32x4 {
    let va = U32x4::load(a);
    let vb = U32x4::load(b);
    if all_normal_bin32(va) && all_normal_bin32(vb) {
        let sign = va.srli(31).xor(vb.srli(31));
        let exp_a = va.srli(23).and(U32x4::splat(0xff));
        let exp_b = vb.srli(23).and(U32x4::splat(0xff));
        let frac = U32x4::splat(0x7f_ffff);
        let hidden = U32x4::splat(1 << 23);
        let sa = va.and(frac).or(hidden).to_array();
        let sb = vb.and(frac).or(hidden).to_array();
        let prod = mul_u32x4_integer(sa, sb);
        let mut out = [0u32; 4];
        for i in 0..4 {
            let e = exp_a.0[i] as i32 + exp_b.0[i] as i32 - BIN32.bias;
            out[i] = normalize_mul(sign.0[i] != 0, e, prod[i] as u128, BIN32) as u32;
        }
        debug_assert_bit_eq32(&out, a, b, mul_bits);
        return out;
    }
    [
        mul_bits(a[0] as u64, b[0] as u64, BIN32) as u32,
        mul_bits(a[1] as u64, b[1] as u64, BIN32) as u32,
        mul_bits(a[2] as u64, b[2] as u64, BIN32) as u32,
        mul_bits(a[3] as u64, b[3] as u64, BIN32) as u32,
    ]
}

fn add_u64x2_integer(a: [u64; 2], b: [u64; 2]) -> [u64; 2] {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        use core::arch::x86_64::{__m128i, _mm_add_epi64, _mm_loadu_si128, _mm_storeu_si128};
        let va = _mm_loadu_si128(a.as_ptr() as *const __m128i);
        let vb = _mm_loadu_si128(b.as_ptr() as *const __m128i);
        let s = _mm_add_epi64(va, vb);
        let mut out = [0u64; 2];
        _mm_storeu_si128(out.as_mut_ptr() as *mut __m128i, s);
        return out;
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        [a[0].wrapping_add(b[0]), a[1].wrapping_add(b[1])]
    }
}

/// Software IEEE add on two binary64 lanes. Integer SIMD on the equal-exp
/// same-sign normal path.
pub fn add_bin64_x2(a: Bin64x2, b: Bin64x2) -> Bin64x2 {
    let ua0 = unpack(a[0], BIN64);
    let ua1 = unpack(a[1], BIN64);
    let ub0 = unpack(b[0], BIN64);
    let ub1 = unpack(b[1], BIN64);
    if ua0.class == Class::Norm
        && ua1.class == Class::Norm
        && ub0.class == Class::Norm
        && ub1.class == Class::Norm
        && ua0.sign == ub0.sign
        && ua1.sign == ub1.sign
        && ua0.exp == ub0.exp
        && ua1.exp == ub1.exp
    {
        let sa = [ua0.sig << 3, ua1.sig << 3];
        let sb = [ub0.sig << 3, ub1.sig << 3];
        let sum = add_u64x2_integer(sa, sb);
        let (e0, c0) = round_rne(sum[0] as u128, ua0.exp, false, BIN64);
        let (e1, c1) = round_rne(sum[1] as u128, ua1.exp, false, BIN64);
        let out = [pack_finite(ua0.sign, e0, c0, BIN64), pack_finite(ua1.sign, e1, c1, BIN64)];
        debug_assert_eq!(out[0], add_bits(a[0], b[0], BIN64));
        debug_assert_eq!(out[1], add_bits(a[1], b[1], BIN64));
        return out;
    }
    [add_bits(a[0], b[0], BIN64), add_bits(a[1], b[1], BIN64)]
}

/// Software IEEE mul on two binary64 lanes.
pub fn mul_bin64_x2(a: Bin64x2, b: Bin64x2) -> Bin64x2 {
    if let (Some(a0), Some(b0)) = (normal_sig64(a[0]), normal_sig64(b[0])) {
        if let (Some(a1), Some(b1)) = (normal_sig64(a[1]), normal_sig64(b[1])) {
            let p0 = mul_u64_integer(a0.sig, b0.sig);
            let p1 = mul_u64_integer(a1.sig, b1.sig);
            let out = [
                normalize_mul(a0.sign ^ b0.sign, a0.exp + b0.exp - BIN64.bias, p0, BIN64),
                normalize_mul(a1.sign ^ b1.sign, a1.exp + b1.exp - BIN64.bias, p1, BIN64),
            ];
            debug_assert_eq!(out[0], mul_bits(a[0], b[0], BIN64));
            debug_assert_eq!(out[1], mul_bits(a[1], b[1], BIN64));
            return out;
        }
    }
    [mul_bits(a[0], b[0], BIN64), mul_bits(a[1], b[1], BIN64)]
}

struct Norm64 {
    sign: bool,
    exp: i32,
    sig: u64,
}

fn normal_sig64(bits: u64) -> Option<Norm64> {
    let u = unpack(bits, BIN64);
    if u.class != Class::Norm {
        return None;
    }
    Some(Norm64 {
        sign: u.sign,
        exp: u.exp,
        sig: u.sig,
    })
}

fn mul_u64_integer(a: u64, b: u64) -> u128 {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        use core::arch::x86_64::{_mm_cvtsi128_si64, _mm_mul_epu32, _mm_set_epi64x};
        // 53-bit significands: split 32+21 and use integer `mul_epu32` on halves.
        let alo = a as u32;
        let ahi = (a >> 32) as u32;
        let blo = b as u32;
        let bhi = (b >> 32) as u32;
        let p0 = _mm_mul_epu32(_mm_set_epi64x(0, alo as i64), _mm_set_epi64x(0, blo as i64));
        let p1 = _mm_mul_epu32(_mm_set_epi64x(0, alo as i64), _mm_set_epi64x(0, bhi as i64));
        let p2 = _mm_mul_epu32(_mm_set_epi64x(0, ahi as i64), _mm_set_epi64x(0, blo as i64));
        let p3 = _mm_mul_epu32(_mm_set_epi64x(0, ahi as i64), _mm_set_epi64x(0, bhi as i64));
        let lo = _mm_cvtsi128_si64(p0) as u128;
        let m1 = _mm_cvtsi128_si64(p1) as u128;
        let m2 = _mm_cvtsi128_si64(p2) as u128;
        let hi = _mm_cvtsi128_si64(p3) as u128;
        return lo + ((m1 + m2) << 32) + (hi << 64);
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        a as u128 * b as u128
    }
}

fn map_pairs_u32(
    a: &[u32],
    b: &[u32],
    chunk: fn(Bin32x4, Bin32x4) -> Bin32x4,
    tail: fn(u64, u64, Format) -> u64,
) -> alloc::vec::Vec<u32> {
    debug_assert_eq!(a.len(), b.len());
    let mut out = alloc::vec::Vec::with_capacity(a.len());
    let mut i = 0;
    while i + 4 <= a.len() {
        let ca = [a[i], a[i + 1], a[i + 2], a[i + 3]];
        let cb = [b[i], b[i + 1], b[i + 2], b[i + 3]];
        out.extend_from_slice(&chunk(ca, cb));
        i += 4;
    }
    while i < a.len() {
        out.push(tail(a[i] as u64, b[i] as u64, BIN32) as u32);
        i += 1;
    }
    out
}

fn map_pairs_u64(
    a: &[u64],
    b: &[u64],
    chunk: fn(Bin64x2, Bin64x2) -> Bin64x2,
    tail: fn(u64, u64, Format) -> u64,
) -> alloc::vec::Vec<u64> {
    debug_assert_eq!(a.len(), b.len());
    let mut out = alloc::vec::Vec::with_capacity(a.len());
    let mut i = 0;
    while i + 2 <= a.len() {
        let r = chunk([a[i], a[i + 1]], [b[i], b[i + 1]]);
        out.extend_from_slice(&r);
        i += 2;
    }
    while i < a.len() {
        out.push(tail(a[i], b[i], BIN64));
        i += 1;
    }
    out
}

/// Elementwise binary32 add over slices of equal length.
pub fn add_u32_lanes(a: &[u32], b: &[u32]) -> alloc::vec::Vec<u32> {
    map_pairs_u32(a, b, add_bin32_x4, add_bits)
}

/// Elementwise binary32 mul over slices of equal length.
pub fn mul_u32_lanes(a: &[u32], b: &[u32]) -> alloc::vec::Vec<u32> {
    map_pairs_u32(a, b, mul_bin32_x4, mul_bits)
}

/// Elementwise binary64 add over slices of equal length.
pub fn add_u64_lanes(a: &[u64], b: &[u64]) -> alloc::vec::Vec<u64> {
    map_pairs_u64(a, b, add_bin64_x2, add_bits)
}

/// Elementwise binary64 mul over slices of equal length.
pub fn mul_u64_lanes(a: &[u64], b: &[u64]) -> alloc::vec::Vec<u64> {
    map_pairs_u64(a, b, mul_bin64_x2, mul_bits)
}

/// Software IEEE sub: flip the sign bit then add.
pub fn sub_bin32_x4(a: Bin32x4, b: Bin32x4) -> Bin32x4 {
    let sign = U32x4::splat(0x8000_0000);
    let nb = U32x4::load(b).xor(sign).to_array();
    add_bin32_x4(a, nb)
}

/// Software IEEE sub on two binary64 lanes.
pub fn sub_bin64_x2(a: Bin64x2, b: Bin64x2) -> Bin64x2 {
    add_bin64_x2(a, [b[0] ^ BIN64.sign_mask(), b[1] ^ BIN64.sign_mask()])
}

/// Software IEEE div on four binary32 lanes. Integer SIMD unpack when every
/// lane is normal; significand quotient is integer `/` per lane (SSE2 has no
/// integer divide). Otherwise the scalar integer kernel.
pub fn div_bin32_x4(a: Bin32x4, b: Bin32x4) -> Bin32x4 {
    let va = U32x4::load(a);
    let vb = U32x4::load(b);
    if all_normal_bin32(va) && all_normal_bin32(vb) {
        let sign = va.srli(31).xor(vb.srli(31));
        let exp_a = va.srli(23).and(U32x4::splat(0xff));
        let exp_b = vb.srli(23).and(U32x4::splat(0xff));
        let frac = U32x4::splat(0x7f_ffff);
        let hidden = U32x4::splat(1 << 23);
        let sa = va.and(frac).or(hidden).to_array();
        let sb = vb.and(frac).or(hidden).to_array();
        let extra = BIN32.frac + 4;
        let mut out = [0u32; 4];
        for i in 0..4 {
            let num = (sa[i] as u128) << extra;
            let den = sb[i] as u128;
            let q = num / den;
            let r = num % den;
            let scale = exp_a.0[i] as i32 - exp_b.0[i] as i32 - extra as i32;
            out[i] = pack_mag(sign.0[i] != 0, scale, q, r != 0, BIN32) as u32;
        }
        debug_assert_bit_eq32(&out, a, b, div_bits);
        return out;
    }
    [
        div_bits(a[0] as u64, b[0] as u64, BIN32) as u32,
        div_bits(a[1] as u64, b[1] as u64, BIN32) as u32,
        div_bits(a[2] as u64, b[2] as u64, BIN32) as u32,
        div_bits(a[3] as u64, b[3] as u64, BIN32) as u32,
    ]
}

/// Software IEEE div on two binary64 lanes.
pub fn div_bin64_x2(a: Bin64x2, b: Bin64x2) -> Bin64x2 {
    if let (Some(a0), Some(b0)) = (normal_sig64(a[0]), normal_sig64(b[0])) {
        if let (Some(a1), Some(b1)) = (normal_sig64(a[1]), normal_sig64(b[1])) {
            let extra = BIN64.frac + 4;
            let out = [
                div_normal_lane(a0, b0, extra, BIN64),
                div_normal_lane(a1, b1, extra, BIN64),
            ];
            debug_assert_eq!(out[0], div_bits(a[0], b[0], BIN64));
            debug_assert_eq!(out[1], div_bits(a[1], b[1], BIN64));
            return out;
        }
    }
    [div_bits(a[0], b[0], BIN64), div_bits(a[1], b[1], BIN64)]
}

fn div_normal_lane(a: Norm64, b: Norm64, extra: u32, f: Format) -> u64 {
    let num = (a.sig as u128) << extra;
    let den = b.sig as u128;
    let q = num / den;
    let r = num % den;
    let scale = a.exp - b.exp - extra as i32;
    pack_mag(a.sign ^ b.sign, scale, q, r != 0, f)
}

/// Software IEEE sqrt on four binary32 lanes. Integer SIMD unpack when every
/// lane is a non-negative normal; `isqrt` per lane. Otherwise the scalar kernel.
pub fn sqrt_bin32_x4(a: Bin32x4) -> Bin32x4 {
    let va = U32x4::load(a);
    if all_normal_bin32(va) && va.srli(31).eq_mask(U32x4::splat(0)) {
        let exp = va.srli(23).and(U32x4::splat(0xff));
        let frac = U32x4::splat(0x7f_ffff);
        let hidden = U32x4::splat(1 << 23);
        let sa = va.and(frac).or(hidden).to_array();
        let mut out = [0u32; 4];
        for i in 0..4 {
            out[i] = sqrt_normal_lane(false, exp.0[i] as i32, sa[i] as u64, BIN32) as u32;
        }
        debug_assert_unary32(&out, a, sqrt_bits);
        return out;
    }
    [
        sqrt_bits(a[0] as u64, BIN32) as u32,
        sqrt_bits(a[1] as u64, BIN32) as u32,
        sqrt_bits(a[2] as u64, BIN32) as u32,
        sqrt_bits(a[3] as u64, BIN32) as u32,
    ]
}

/// Software IEEE sqrt on two binary64 lanes.
pub fn sqrt_bin64_x2(a: Bin64x2) -> Bin64x2 {
    if let (Some(a0), Some(a1)) = (normal_sig64(a[0]), normal_sig64(a[1])) {
        if !a0.sign && !a1.sign {
            let out = [
                sqrt_normal_lane(a0.sign, a0.exp, a0.sig, BIN64),
                sqrt_normal_lane(a1.sign, a1.exp, a1.sig, BIN64),
            ];
            debug_assert_eq!(out[0], sqrt_bits(a[0], BIN64));
            debug_assert_eq!(out[1], sqrt_bits(a[1], BIN64));
            return out;
        }
    }
    [sqrt_bits(a[0], BIN64), sqrt_bits(a[1], BIN64)]
}

fn sqrt_normal_lane(_sign: bool, exp: i32, sig: u64, f: Format) -> u64 {
    let mut s = sig as u128;
    let mut exp2 = exp - f.bias - f.frac as i32;
    if exp2 & 1 != 0 {
        s <<= 1;
        exp2 -= 1;
    }
    let extra = 64u32;
    s <<= extra;
    let root = isqrt(s);
    let rem = s - root * root;
    let scale = exp2 / 2 - extra as i32 / 2;
    pack_mag(false, scale, root, rem != 0, f)
}

fn debug_assert_unary32(got: &Bin32x4, a: Bin32x4, op: fn(u64, Format) -> u64) {
    for i in 0..4 {
        debug_assert_eq!(
            got[i],
            op(a[i] as u64, BIN32) as u32,
            "SIMD lane bits must match the scalar integer kernel"
        );
    }
}

/// Software IEEE FMA \(a\cdot b + c\) on four binary32 lanes. Integer SIMD
/// significand products when every lane is normal; then the scalar `fma_add`.
pub fn fma_bin32_x4(a: Bin32x4, b: Bin32x4, c: Bin32x4) -> Bin32x4 {
    let va = U32x4::load(a);
    let vb = U32x4::load(b);
    let vc = U32x4::load(c);
    if all_normal_bin32(va) && all_normal_bin32(vb) && all_normal_bin32(vc) {
        let sign = va.srli(31).xor(vb.srli(31));
        let exp_a = va.srli(23).and(U32x4::splat(0xff));
        let exp_b = vb.srli(23).and(U32x4::splat(0xff));
        let frac = U32x4::splat(0x7f_ffff);
        let hidden = U32x4::splat(1 << 23);
        let sa = va.and(frac).or(hidden).to_array();
        let sb = vb.and(frac).or(hidden).to_array();
        let prod = mul_u32x4_integer(sa, sb);
        let mut out = [0u32; 4];
        for i in 0..4 {
            let pe = exp_a.0[i] as i32 + exp_b.0[i] as i32 - BIN32.bias;
            let uc = unpack(c[i] as u64, BIN32);
            out[i] = fma_add(sign.0[i] != 0, pe, prod[i] as u128, uc, BIN32) as u32;
        }
        debug_assert_fma32(&out, a, b, c);
        return out;
    }
    [
        fma_bits(a[0] as u64, b[0] as u64, c[0] as u64, BIN32) as u32,
        fma_bits(a[1] as u64, b[1] as u64, c[1] as u64, BIN32) as u32,
        fma_bits(a[2] as u64, b[2] as u64, c[2] as u64, BIN32) as u32,
        fma_bits(a[3] as u64, b[3] as u64, c[3] as u64, BIN32) as u32,
    ]
}

/// Software IEEE FMA on two binary64 lanes.
pub fn fma_bin64_x2(a: Bin64x2, b: Bin64x2, c: Bin64x2) -> Bin64x2 {
    if let (Some(a0), Some(b0)) = (normal_sig64(a[0]), normal_sig64(b[0])) {
        if let (Some(a1), Some(b1)) = (normal_sig64(a[1]), normal_sig64(b[1])) {
            let u0 = unpack(c[0], BIN64);
            let u1 = unpack(c[1], BIN64);
            if u0.class == Class::Norm && u1.class == Class::Norm {
                let p0 = mul_u64_integer(a0.sig, b0.sig);
                let p1 = mul_u64_integer(a1.sig, b1.sig);
                let out = [
                    fma_add(
                        a0.sign ^ b0.sign,
                        a0.exp + b0.exp - BIN64.bias,
                        p0,
                        u0,
                        BIN64,
                    ),
                    fma_add(
                        a1.sign ^ b1.sign,
                        a1.exp + b1.exp - BIN64.bias,
                        p1,
                        u1,
                        BIN64,
                    ),
                ];
                debug_assert_eq!(out[0], fma_bits(a[0], b[0], c[0], BIN64));
                debug_assert_eq!(out[1], fma_bits(a[1], b[1], c[1], BIN64));
                return out;
            }
        }
    }
    [
        fma_bits(a[0], b[0], c[0], BIN64),
        fma_bits(a[1], b[1], c[1], BIN64),
    ]
}

fn debug_assert_fma32(got: &Bin32x4, a: Bin32x4, b: Bin32x4, c: Bin32x4) {
    for i in 0..4 {
        debug_assert_eq!(
            got[i],
            fma_bits(a[i] as u64, b[i] as u64, c[i] as u64, BIN32) as u32,
            "SIMD lane bits must match the scalar integer kernel"
        );
    }
}

fn map_unary_u32(a: &[u32], chunk: fn(Bin32x4) -> Bin32x4, tail: fn(u64, Format) -> u64) -> alloc::vec::Vec<u32> {
    let mut out = alloc::vec::Vec::with_capacity(a.len());
    let mut i = 0;
    while i + 4 <= a.len() {
        out.extend_from_slice(&chunk([a[i], a[i + 1], a[i + 2], a[i + 3]]));
        i += 4;
    }
    while i < a.len() {
        out.push(tail(a[i] as u64, BIN32) as u32);
        i += 1;
    }
    out
}

fn map_unary_u64(a: &[u64], chunk: fn(Bin64x2) -> Bin64x2, tail: fn(u64, Format) -> u64) -> alloc::vec::Vec<u64> {
    let mut out = alloc::vec::Vec::with_capacity(a.len());
    let mut i = 0;
    while i + 2 <= a.len() {
        out.extend_from_slice(&chunk([a[i], a[i + 1]]));
        i += 2;
    }
    while i < a.len() {
        out.push(tail(a[i], BIN64));
        i += 1;
    }
    out
}

fn map_triples_u32(
    a: &[u32],
    b: &[u32],
    c: &[u32],
    chunk: fn(Bin32x4, Bin32x4, Bin32x4) -> Bin32x4,
    tail: fn(u64, u64, u64, Format) -> u64,
) -> alloc::vec::Vec<u32> {
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), c.len());
    let mut out = alloc::vec::Vec::with_capacity(a.len());
    let mut i = 0;
    while i + 4 <= a.len() {
        let ca = [a[i], a[i + 1], a[i + 2], a[i + 3]];
        let cb = [b[i], b[i + 1], b[i + 2], b[i + 3]];
        let cc = [c[i], c[i + 1], c[i + 2], c[i + 3]];
        out.extend_from_slice(&chunk(ca, cb, cc));
        i += 4;
    }
    while i < a.len() {
        out.push(tail(a[i] as u64, b[i] as u64, c[i] as u64, BIN32) as u32);
        i += 1;
    }
    out
}

fn map_triples_u64(
    a: &[u64],
    b: &[u64],
    c: &[u64],
    chunk: fn(Bin64x2, Bin64x2, Bin64x2) -> Bin64x2,
    tail: fn(u64, u64, u64, Format) -> u64,
) -> alloc::vec::Vec<u64> {
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), c.len());
    let mut out = alloc::vec::Vec::with_capacity(a.len());
    let mut i = 0;
    while i + 2 <= a.len() {
        let r = chunk(
            [a[i], a[i + 1]],
            [b[i], b[i + 1]],
            [c[i], c[i + 1]],
        );
        out.extend_from_slice(&r);
        i += 2;
    }
    while i < a.len() {
        out.push(tail(a[i], b[i], c[i], BIN64));
        i += 1;
    }
    out
}

/// Elementwise binary32 sub over slices of equal length.
pub fn sub_u32_lanes(a: &[u32], b: &[u32]) -> alloc::vec::Vec<u32> {
    map_pairs_u32(a, b, sub_bin32_x4, sub_bits)
}

/// Elementwise binary32 div over slices of equal length.
pub fn div_u32_lanes(a: &[u32], b: &[u32]) -> alloc::vec::Vec<u32> {
    map_pairs_u32(a, b, div_bin32_x4, div_bits)
}

/// Elementwise binary32 sqrt.
pub fn sqrt_u32_lanes(a: &[u32]) -> alloc::vec::Vec<u32> {
    map_unary_u32(a, sqrt_bin32_x4, sqrt_bits)
}

/// Elementwise binary32 FMA.
pub fn fma_u32_lanes(a: &[u32], b: &[u32], c: &[u32]) -> alloc::vec::Vec<u32> {
    map_triples_u32(a, b, c, fma_bin32_x4, fma_bits)
}

/// Elementwise binary64 sub over slices of equal length.
pub fn sub_u64_lanes(a: &[u64], b: &[u64]) -> alloc::vec::Vec<u64> {
    map_pairs_u64(a, b, sub_bin64_x2, sub_bits)
}

/// Elementwise binary64 div over slices of equal length.
pub fn div_u64_lanes(a: &[u64], b: &[u64]) -> alloc::vec::Vec<u64> {
    map_pairs_u64(a, b, div_bin64_x2, div_bits)
}

/// Elementwise binary64 sqrt.
pub fn sqrt_u64_lanes(a: &[u64]) -> alloc::vec::Vec<u64> {
    map_unary_u64(a, sqrt_bin64_x2, sqrt_bits)
}

/// Elementwise binary64 FMA.
pub fn fma_u64_lanes(a: &[u64], b: &[u64], c: &[u64]) -> alloc::vec::Vec<u64> {
    map_triples_u64(a, b, c, fma_bin64_x2, fma_bits)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ieee_soft::Ieee32;

    #[test]
    fn bin32_x4_add_matches_scalar() {
        let one = Ieee32::from_i32(1).to_bits();
        let two = Ieee32::from_i32(2).to_bits();
        let a = [one, two, one, two];
        let b = [one, two, two, one];
        let got = add_bin32_x4(a, b);
        for i in 0..4 {
            assert_eq!(got[i], add_bits(a[i] as u64, b[i] as u64, BIN32) as u32);
        }
        assert_eq!(IEEE_SIMD_LANE_WIDTH, 4);
        assert_eq!(got[0], Ieee32::from_i32(2).to_bits());
        assert_eq!(got[1], Ieee32::from_i32(4).to_bits());
    }

    #[test]
    fn bin32_x4_mul_half() {
        let two = Ieee32::from_i32(2).to_bits();
        let half = 0x3F00_0000u32;
        let a = [two, two, two, two];
        let b = [half, half, half, half];
        let got = mul_bin32_x4(a, b);
        let one = Ieee32::from_i32(1).to_bits();
        assert_eq!(got, [one, one, one, one]);
    }

    #[test]
    fn bin32_x4_specials_match_scalar() {
        let one = Ieee32::from_i32(1).to_bits();
        let inf = 0x7F80_0000u32;
        let a = [one, inf, 0, one];
        let b = [one, one, one, 0];
        let got = add_bin32_x4(a, b);
        for i in 0..4 {
            assert_eq!(got[i], add_bits(a[i] as u64, b[i] as u64, BIN32) as u32);
        }
    }

    #[test]
    fn simd_add_mul_bit_identical_to_scalar() {
        let samples = [
            0u32,
            0x0000_0001,
            0x0080_0000,
            0x3F00_0000,
            0x3F80_0000,
            0x3F80_0001,
            0x4000_0000,
            0x7F7F_FFFF,
            0x7F80_0000,
            0x7FC0_0000,
            0x8000_0000,
            0xBF80_0000,
            Ieee32::from_i32(-3).to_bits(),
            Ieee32::from_i32(7).to_bits(),
        ];
        for (i, &ai) in samples.iter().enumerate() {
            for (j, &bj) in samples.iter().enumerate() {
                let a = [
                    ai,
                    samples[(i + 1) % samples.len()],
                    samples[(i + 2) % samples.len()],
                    samples[(i + 3) % samples.len()],
                ];
                let b = [
                    bj,
                    samples[(j + 1) % samples.len()],
                    samples[(j + 2) % samples.len()],
                    samples[(j + 3) % samples.len()],
                ];
                let add = add_bin32_x4(a, b);
                let mul = mul_bin32_x4(a, b);
                for k in 0..4 {
                    assert_eq!(add[k], add_bits(a[k] as u64, b[k] as u64, BIN32) as u32);
                    assert_eq!(mul[k], mul_bits(a[k] as u64, b[k] as u64, BIN32) as u32);
                }
            }
        }
        let one = crate::ieee_soft::Ieee64::from_i32(1).to_bits();
        let two = crate::ieee_soft::Ieee64::from_i32(2).to_bits();
        let add64 = add_bin64_x2([one, two], [two, one]);
        assert_eq!(add64[0], add_bits(one, two, BIN64));
        assert_eq!(add64[1], add_bits(two, one, BIN64));
        let mul64 = mul_bin64_x2([two, two], [one, one]);
        assert_eq!(mul64[0], mul_bits(two, one, BIN64));
        assert_eq!(mul64[1], mul_bits(two, one, BIN64));
        let four = crate::ieee_soft::Ieee64::from_i32(4).to_bits();
        let div64 = div_bin64_x2([four, two], [two, one]);
        assert_eq!(div64[0], div_bits(four, two, BIN64));
        assert_eq!(div64[1], div_bits(two, one, BIN64));
        let sq64 = sqrt_bin64_x2([four, four]);
        assert_eq!(sq64[0], sqrt_bits(four, BIN64));
        assert_eq!(sq64[1], two);
        let fma64 = fma_bin64_x2([two, two], [two, one], [one, one]);
        assert_eq!(fma64[0], fma_bits(two, two, one, BIN64));
        assert_eq!(fma64[1], fma_bits(two, one, one, BIN64));
    }

    #[test]
    fn simd_div_sqrt_fma_match_scalar_samples() {
        let samples = [
            0u32,
            0x0000_0001,
            0x0080_0000,
            0x3F00_0000,
            0x3F80_0000,
            0x3F80_0001,
            0x4000_0000,
            0x7F7F_FFFF,
            0x7F80_0000,
            0x7FC0_0000,
            0x8000_0000,
            0xBF80_0000,
            Ieee32::from_i32(-3).to_bits(),
            Ieee32::from_i32(7).to_bits(),
        ];
        for (i, &ai) in samples.iter().enumerate() {
            for (j, &bj) in samples.iter().enumerate() {
                let a = [
                    ai,
                    samples[(i + 1) % samples.len()],
                    samples[(i + 2) % samples.len()],
                    samples[(i + 3) % samples.len()],
                ];
                let b = [
                    bj,
                    samples[(j + 1) % samples.len()],
                    samples[(j + 2) % samples.len()],
                    samples[(j + 3) % samples.len()],
                ];
                let div = div_bin32_x4(a, b);
                let sub = sub_bin32_x4(a, b);
                let sq = sqrt_bin32_x4(a);
                let fma = fma_bin32_x4(a, b, a);
                for k in 0..4 {
                    assert_eq!(div[k], div_bits(a[k] as u64, b[k] as u64, BIN32) as u32);
                    assert_eq!(sub[k], sub_bits(a[k] as u64, b[k] as u64, BIN32) as u32);
                    assert_eq!(sq[k], sqrt_bits(a[k] as u64, BIN32) as u32);
                    assert_eq!(
                        fma[k],
                        fma_bits(a[k] as u64, b[k] as u64, a[k] as u64, BIN32) as u32
                    );
                }
            }
        }
    }
}
