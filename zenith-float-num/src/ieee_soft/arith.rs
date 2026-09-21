//! Integer IEEE-754 binary32 / binary64 arithmetic (to-nearest, ties to even).

#[derive(Clone, Copy)]
pub(super) struct Format {
    pub frac: u32,
    pub exp_bits: u32,
    pub bias: i32,
}

pub(super) const BIN32: Format = Format {
    frac: 23,
    exp_bits: 8,
    bias: 127,
};

pub(super) const BIN64: Format = Format {
    frac: 52,
    exp_bits: 11,
    bias: 1023,
};

impl Format {
    pub fn sign_shift(self) -> u32 {
        self.frac + self.exp_bits
    }

    pub fn exp_max(self) -> i32 {
        (1 << self.exp_bits) - 1
    }

    pub fn exp_mask(self) -> u64 {
        (1u64 << self.exp_bits) - 1
    }

    pub fn frac_mask(self) -> u64 {
        (1u64 << self.frac) - 1
    }

    pub fn hidden(self) -> u64 {
        1u64 << self.frac
    }

    pub fn sign_mask(self) -> u64 {
        1u64 << self.sign_shift()
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum Class {
    Zero,
    Sub,
    Norm,
    Inf,
    Nan,
}

#[derive(Clone, Copy)]
pub(super) struct Unp {
    pub sign: bool,
    pub exp: i32,
    pub sig: u64,
    pub class: Class,
}

pub(super) fn unpack(bits: u64, f: Format) -> Unp {
    let sign = (bits >> f.sign_shift()) != 0;
    let exp = ((bits >> f.frac) & f.exp_mask()) as i32;
    let frac = bits & f.frac_mask();
    if exp == f.exp_max() {
        if frac != 0 {
            Unp {
                sign,
                exp,
                sig: frac,
                class: Class::Nan,
            }
        } else {
            Unp {
                sign,
                exp,
                sig: 0,
                class: Class::Inf,
            }
        }
    } else if exp == 0 {
        if frac == 0 {
            Unp {
                sign,
                exp: 0,
                sig: 0,
                class: Class::Zero,
            }
        } else {
            Unp {
                sign,
                exp: 1,
                sig: frac,
                class: Class::Sub,
            }
        }
    } else {
        Unp {
            sign,
            exp,
            sig: frac | f.hidden(),
            class: Class::Norm,
        }
    }
}

pub(super) fn pack_finite(sign: bool, mut exp: i32, mut sig: u128, f: Format) -> u64 {
    let hidden = f.hidden() as u128;
    if sig == 0 {
        return if sign { f.sign_mask() } else { 0 };
    }
    if exp >= f.exp_max() {
        return pack_inf(sign, f);
    }
    if exp < 1 || sig < hidden {
        if exp >= 1 && sig < hidden {
            return (if sign { f.sign_mask() } else { 0 }) | (sig as u64 & f.frac_mask());
        }
        let sh = (1 - exp) as u32;
        let sticky = if sh >= 128 {
            let s = sig != 0;
            sig = 0;
            s
        } else {
            let mask = (1u128 << sh) - 1;
            let s = (sig & mask) != 0;
            sig >>= sh;
            s
        };
        if sticky {
            sig |= 1;
        }
        exp = 0;
        if sig >= hidden {
            sig >>= 1;
            exp = 1;
        }
    }
    let frac = (sig as u64) & f.frac_mask();
    let e = (exp as u64) & f.exp_mask();
    (if sign { f.sign_mask() } else { 0 }) | (e << f.frac) | frac
}

pub(super) fn pack_inf(sign: bool, f: Format) -> u64 {
    (if sign { f.sign_mask() } else { 0 }) | ((f.exp_max() as u64) << f.frac)
}

pub(super) fn pack_qnan(f: Format) -> u64 {
    ((f.exp_max() as u64) << f.frac) | (1u64 << (f.frac - 1))
}

pub(super) fn canonical_nan(a: Unp, b: Option<Unp>, f: Format) -> u64 {
    let _ = (a, b);
    pack_qnan(f)
}

const GRS: u32 = 3;

pub(super) fn round_rne(
    mut sig: u128,
    mut exp: i32,
    mut extra_sticky: bool,
    f: Format,
) -> (i32, u128) {
    let target = f.frac + GRS;
    while sig >= (1u128 << (target + 1)) {
        extra_sticky |= (sig & 1) != 0;
        sig >>= 1;
        exp += 1;
    }
    while sig != 0 && sig < (1u128 << target) && exp > 1 {
        sig <<= 1;
        exp -= 1;
    }
    let grs = (sig as u32) & 7;
    let sticky_bit = extra_sticky || (grs & 1) != 0;
    let g = (grs >> 2) & 1;
    let r = (grs >> 1) & 1;
    let mut core = sig >> GRS;
    if g == 1 && (r == 1 || sticky_bit || (core & 1) == 1) {
        core += 1;
        if core >= (1u128 << (f.frac + 1)) {
            core >>= 1;
            exp += 1;
        }
    }
    (exp, core)
}

pub(super) fn add_bits(a: u64, b: u64, f: Format) -> u64 {
    let ua = unpack(a, f);
    let ub = unpack(b, f);
    match (ua.class, ub.class) {
        (Class::Nan, _) | (_, Class::Nan) => canonical_nan(ua, Some(ub), f),
        (Class::Inf, Class::Inf) if ua.sign != ub.sign => pack_qnan(f),
        (Class::Inf, _) => pack_inf(ua.sign, f),
        (_, Class::Inf) => pack_inf(ub.sign, f),
        (Class::Zero, Class::Zero) => {
            if ua.sign && ub.sign {
                f.sign_mask()
            } else {
                0
            }
        }
        (Class::Zero, _) => b,
        (_, Class::Zero) => a,
        _ => add_finite(ua, ub, f),
    }
}

fn add_finite(mut ua: Unp, mut ub: Unp, f: Format) -> u64 {
    if ua.exp < ub.exp || (ua.exp == ub.exp && ua.sig < ub.sig) {
        core::mem::swap(&mut ua, &mut ub);
    }
    let ea = ua.exp;
    let sa = (ua.sig as u128) << GRS;
    let mut sb = (ub.sig as u128) << GRS;
    let diff = (ea - ub.exp) as u32;
    let mut sticky = false;
    if diff > 0 {
        if diff >= 128 {
            sticky = sb != 0;
            sb = 0;
        } else {
            let mask = (1u128 << diff) - 1;
            sticky = (sb & mask) != 0;
            sb >>= diff;
        }
    }
    let (sign, mut sum) = if ua.sign == ub.sign {
        (ua.sign, sa + sb)
    } else if sa > sb {
        (ua.sign, sa - sb)
    } else if sb > sa {
        (ub.sign, sb - sa)
    } else {
        return 0;
    };
    if sticky {
        sum |= 1;
    }
    if sum == 0 {
        return 0;
    }
    let (exp, core) = round_rne(sum, ea, false, f);
    pack_finite(sign, exp, core, f)
}

pub(super) fn sub_bits(a: u64, b: u64, f: Format) -> u64 {
    add_bits(a, b ^ f.sign_mask(), f)
}

pub(super) fn mul_bits(a: u64, b: u64, f: Format) -> u64 {
    let ua = unpack(a, f);
    let ub = unpack(b, f);
    let sign = ua.sign ^ ub.sign;
    match (ua.class, ub.class) {
        (Class::Nan, _) | (_, Class::Nan) => canonical_nan(ua, Some(ub), f),
        (Class::Inf, Class::Zero) | (Class::Zero, Class::Inf) => pack_qnan(f),
        (Class::Inf, _) | (_, Class::Inf) => pack_inf(sign, f),
        (Class::Zero, _) | (_, Class::Zero) => {
            if sign {
                f.sign_mask()
            } else {
                0
            }
        }
        _ => {
            let prod = (ua.sig as u128) * (ub.sig as u128);
            let e = ua.exp + ub.exp - f.bias;
            normalize_mul(sign, e, prod, f)
        }
    }
}

/// Pack `± mag * 2^{scale_exp}` as an IEEE value.
pub(super) fn pack_mag(sign: bool, scale_exp: i32, mag: u128, mut sticky: bool, f: Format) -> u64 {
    if mag == 0 {
        return if sign { f.sign_mask() } else { 0 };
    }
    let hi = 127 - mag.leading_zeros();
    let take = f.frac + 1 + GRS;
    let aligned = if hi + 1 > take {
        let sh = hi + 1 - take;
        let mask = (1u128 << sh) - 1;
        sticky |= (mag & mask) != 0;
        mag >> sh
    } else {
        mag << (take - hi - 1)
    };
    let sh = if hi + 1 > take { (hi + 1 - take) as i32 } else { -((take - hi - 1) as i32) };
    let stored = scale_exp + sh + GRS as i32 + f.bias + f.frac as i32;
    let (exp, core) = round_rne(aligned, stored, sticky, f);
    pack_finite(sign, exp, core, f)
}

pub(super) fn normalize_mul(sign: bool, exp_sum: i32, prod: u128, f: Format) -> u64 {
    let scale = exp_sum - f.bias - 2 * f.frac as i32;
    pack_mag(sign, scale, prod, false, f)
}

pub(super) fn div_bits(a: u64, b: u64, f: Format) -> u64 {
    let ua = unpack(a, f);
    let ub = unpack(b, f);
    let sign = ua.sign ^ ub.sign;
    match (ua.class, ub.class) {
        (Class::Nan, _) | (_, Class::Nan) => canonical_nan(ua, Some(ub), f),
        (Class::Inf, Class::Inf) | (Class::Zero, Class::Zero) => pack_qnan(f),
        (_, Class::Zero) => pack_inf(sign, f),
        (Class::Zero, _) => {
            if sign {
                f.sign_mask()
            } else {
                0
            }
        }
        (Class::Inf, _) => pack_inf(sign, f),
        (_, Class::Inf) => {
            if sign {
                f.sign_mask()
            } else {
                0
            }
        }
        _ => {
            let extra = f.frac + 4;
            let num = (ua.sig as u128) << extra;
            let den = ub.sig as u128;
            let q = num / den;
            let r = num % den;
            let scale = ua.exp - ub.exp - extra as i32;
            pack_mag(sign, scale, q, r != 0, f)
        }
    }
}

pub(super) fn isqrt(n: u128) -> u128 {
    if n < 2 {
        return n;
    }
    let mut x = n;
    let mut y = (x + 1) / 2;
    while y < x {
        x = y;
        y = (x + n / x) / 2;
    }
    x
}

pub(super) fn sqrt_bits(a: u64, f: Format) -> u64 {
    let ua = unpack(a, f);
    match ua.class {
        Class::Nan => pack_qnan(f),
        Class::Zero => a,
        Class::Inf if !ua.sign => pack_inf(false, f),
        _ if ua.sign => pack_qnan(f),
        _ => {
            let mut sig = ua.sig as u128;
            let mut exp2 = ua.exp - f.bias - f.frac as i32;
            if exp2 & 1 != 0 {
                sig <<= 1;
                exp2 -= 1;
            }
            let extra = 64u32;
            sig <<= extra;
            let root = isqrt(sig);
            let rem = sig - root * root;
            let scale = exp2 / 2 - extra as i32 / 2;
            pack_mag(false, scale, root, rem != 0, f)
        }
    }
}

pub(super) fn fma_bits(a: u64, b: u64, c: u64, f: Format) -> u64 {
    let ua = unpack(a, f);
    let ub = unpack(b, f);
    let uc = unpack(c, f);
    if ua.class == Class::Nan || ub.class == Class::Nan || uc.class == Class::Nan {
        return pack_qnan(f);
    }
    let psign = ua.sign ^ ub.sign;
    match (ua.class, ub.class) {
        (Class::Inf, Class::Zero) | (Class::Zero, Class::Inf) => return pack_qnan(f),
        (Class::Inf, _) | (_, Class::Inf) => {
            if uc.class == Class::Inf && uc.sign != psign {
                return pack_qnan(f);
            }
            return pack_inf(psign, f);
        }
        (Class::Zero, _) | (_, Class::Zero) => return c,
        _ => {}
    }
    if uc.class == Class::Inf {
        return pack_inf(uc.sign, f);
    }
    if uc.class == Class::Zero {
        return mul_bits(a, b, f);
    }
    let prod = (ua.sig as u128) * (ub.sig as u128);
    let pe = ua.exp + ub.exp - f.bias;
    fma_add(psign, pe, prod, uc, f)
}

pub(super) fn fma_add(psign: bool, pe: i32, prod: u128, uc: Unp, f: Format) -> u64 {
    if prod == 0 {
        return pack_finite(uc.sign, uc.exp, uc.sig as u128, f);
    }
    let mut mp = prod;
    let mut sp = pe - f.bias - 2 * f.frac as i32;
    let mut mc = uc.sig as u128;
    let sc = uc.exp - f.bias - f.frac as i32;
    let mut sticky = false;
    if sp > sc {
        let d = (sp - sc) as u32;
        if d >= 128 {
            sticky = mc != 0;
            mc = 0;
        } else {
            let mask = (1u128 << d) - 1;
            sticky = (mc & mask) != 0;
            mc >>= d;
        }
    } else if sc > sp {
        let d = (sc - sp) as u32;
        if d >= 128 {
            sticky = mp != 0;
            mp = 0;
        } else {
            let mask = (1u128 << d) - 1;
            sticky = (mp & mask) != 0;
            mp >>= d;
        }
        sp = sc;
    }
    let (sign, mut sum) = if psign == uc.sign {
        (psign, mp + mc)
    } else if mp > mc {
        (psign, mp - mc)
    } else if mc > mp {
        (uc.sign, mc - mp)
    } else {
        return 0;
    };
    if sticky {
        sum |= 1;
    }
    pack_mag(sign, sp, sum, sticky, f)
}

pub(super) fn next_up(bits: u64, f: Format) -> u64 {
    let u = unpack(bits, f);
    if u.class == Class::Nan {
        return pack_qnan(f);
    }
    if u.class == Class::Inf && !u.sign {
        return bits;
    }
    if bits == f.sign_mask() || bits == 0 {
        return 1;
    }
    if u.sign {
        bits.wrapping_sub(1)
    } else {
        bits.wrapping_add(1)
    }
}

pub(super) fn next_down(bits: u64, f: Format) -> u64 {
    let u = unpack(bits, f);
    if u.class == Class::Nan {
        return pack_qnan(f);
    }
    if u.class == Class::Inf && u.sign {
        return bits;
    }
    if bits == 0 {
        return f.sign_mask() | 1;
    }
    if bits == f.sign_mask() {
        return f.sign_mask() | 1;
    }
    if u.sign {
        bits.wrapping_add(1)
    } else {
        bits.wrapping_sub(1)
    }
}

pub(super) fn next_after(x: u64, y: u64, f: Format) -> u64 {
    if unpack(x, f).class == Class::Nan || unpack(y, f).class == Class::Nan {
        return pack_qnan(f);
    }
    match cmp_bits(x, y, f) {
        Some(0) => x,
        Some(d) if d < 0 => next_up(x, f),
        Some(_) => next_down(x, f),
        None => pack_qnan(f),
    }
}

pub(super) fn cmp_bits(a: u64, b: u64, f: Format) -> Option<i32> {
    let ua = unpack(a, f);
    let ub = unpack(b, f);
    if ua.class == Class::Nan || ub.class == Class::Nan {
        return None;
    }
    if ua.class == Class::Zero && ub.class == Class::Zero {
        return Some(0);
    }
    match (ua.sign, ub.sign) {
        (true, false) => Some(-1),
        (false, true) => Some(1),
        (false, false) => Some(ord_mag(a, b, f)),
        (true, true) => Some(-ord_mag(a, b, f)),
    }
}

fn ord_mag(a: u64, b: u64, f: Format) -> i32 {
    let mask = f.sign_mask() - 1;
    let am = a & mask;
    let bm = b & mask;
    am.cmp(&bm) as i32
}

pub(super) fn frexp_bits(bits: u64, f: Format) -> (u64, i32) {
    let u = unpack(bits, f);
    match u.class {
        Class::Nan | Class::Inf | Class::Zero => (bits, 0),
        Class::Norm => {
            let e = u.exp - f.bias + 1;
            let mexp = f.bias - 1;
            let m = (if u.sign { f.sign_mask() } else { 0 })
                | ((mexp as u64) << f.frac)
                | (u.sig & f.frac_mask());
            (m, e)
        }
        Class::Sub => {
            let lz = u.sig.leading_zeros() - (64 - f.frac);
            let e = 1 - f.bias + 1 - lz as i32;
            let sig = u.sig << lz;
            let mexp = f.bias - 1;
            let m = (if u.sign { f.sign_mask() } else { 0 })
                | ((mexp as u64) << f.frac)
                | (sig & f.frac_mask());
            (m, e)
        }
    }
}

pub(super) fn from_i32_bits(n: i32, f: Format) -> u64 {
    if n == 0 {
        return 0;
    }
    let sign = n < 0;
    let mut mag = (n as i64).unsigned_abs();
    let lz = mag.leading_zeros();
    let msb = 63 - lz;
    let exp = msb as i32 + f.bias;
    mag <<= lz;
    let frac = (mag << 1) >> (64 - f.frac);
    let leftover = mag << (1 + f.frac);
    let mut core = (f.hidden() | frac) as u128;
    let g = (leftover >> 63) & 1;
    let sticky = leftover << 1 != 0;
    if g == 1 && (sticky || (core & 1) == 1) {
        core += 1;
    }
    pack_finite(sign, exp, core, f)
}
