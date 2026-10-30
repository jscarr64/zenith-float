//! Independent differential oracle: the Mag operations computed directly with u128/i128 integers (no shared code with the library or with the dyadic model).
use super::dyadic::MagV;
use core::cmp::Ordering;

fn min_pos() -> MagV {
    MagV::Finite {
        exp: i64::MIN,
        mant: 1u64 << 63,
    }
}

fn zero() -> MagV {
    MagV::Finite { exp: 0, mant: 0 }
}

fn pack(exp: i128, mant: u64, up: bool) -> Option<MagV> {
    if exp > i128::from(i64::MAX) {
        if up {
            Some(MagV::Infinity)
        } else {
            Some(MagV::Finite {
                exp: i64::MAX,
                mant: u64::MAX,
            })
        }
    } else if exp < i128::from(i64::MIN) {
        if up {
            Some(min_pos())
        } else {
            Some(zero())
        }
    } else {
        Some(MagV::Finite {
            exp: i64::try_from(exp).ok()?,
            mant,
        })
    }
}

fn norm(m: u128, e: i128, sticky: bool, up: bool) -> Option<MagV> {
    let bl = 128 - i128::from(m.leading_zeros());
    let mut m = m;
    let mut e = e;
    let mut sticky = sticky;
    if bl > 64 {
        let sh = bl - 64;
        let low = m & ((1u128 << sh) - 1);
        m >>= sh;
        e += sh;
        if low != 0 {
            sticky = true;
        }
    } else if bl < 64 {
        let sh = 64 - bl;
        m <<= sh;
        e -= sh;
    }
    let mut mant = u64::try_from(m).ok()?;
    let mut exp = e + 64;
    if up && sticky {
        if mant == u64::MAX {
            mant = 1u64 << 63;
            exp += 1;
        } else {
            mant += 1;
        }
    }
    pack(exp, mant, up)
}

fn isqrt(n: u128) -> u128 {
    let mut r: u128 = 0;
    let mut bit: u128 = 1u128 << 126;
    let mut n = n;
    while bit > n {
        bit >>= 2;
    }
    while bit != 0 {
        if n >= r + bit {
            n -= r + bit;
            r = (r >> 1) + bit;
        } else {
            r >>= 1;
        }
        bit >>= 2;
    }
    r
}

fn parts(x: &MagV) -> Option<(i128, u128)> {
    match x {
        MagV::Finite { exp, mant } if *mant != 0 => Some((i128::from(*exp), u128::from(*mant))),
        MagV::Finite { .. } | MagV::Infinity => None,
    }
}

fn is_zero(x: &MagV) -> bool {
    matches!(x, MagV::Finite { mant: 0, .. })
}

fn succ(x: &MagV) -> Option<MagV> {
    match x {
        MagV::Finite { exp, mant } => {
            if *mant == u64::MAX {
                pack(i128::from(*exp) + 1, 1u64 << 63, true)
            } else {
                Some(MagV::Finite {
                    exp: *exp,
                    mant: *mant + 1,
                })
            }
        }
        MagV::Infinity => Some(MagV::Infinity),
    }
}

/// Smallest Mag at or above a + b.
pub fn add_up(a: &MagV, b: &MagV) -> Option<MagV> {
    if matches!(a, MagV::Infinity) || matches!(b, MagV::Infinity) {
        return Some(MagV::Infinity);
    }
    match (parts(a), parts(b)) {
        (None, Some((e, m))) | (Some((e, m)), None) => Some(MagV::Finite {
            exp: i64::try_from(e).ok()?,
            mant: u64::try_from(m).ok()?,
        }),
        (None, None) => Some(zero()),
        (Some((ea, ma)), Some((eb, mb))) => {
            let (ea, ma, eb, mb) = if ea >= eb {
                (ea, ma, eb, mb)
            } else {
                (eb, mb, ea, ma)
            };
            let gap = ea - eb;
            if gap >= 64 {
                succ(&pack(ea, u64::try_from(ma).ok()?, true)?)
            } else {
                norm((ma << gap) + mb, eb - 64, false, true)
            }
        }
    }
}

/// Smallest Mag at or above a * b.
pub fn mul_up(a: &MagV, b: &MagV) -> Option<MagV> {
    if matches!(a, MagV::Infinity) || matches!(b, MagV::Infinity) {
        return Some(MagV::Infinity);
    }
    match (parts(a), parts(b)) {
        (Some((ea, ma)), Some((eb, mb))) => norm(ma * mb, ea + eb - 128, false, true),
        _ => Some(zero()),
    }
}

/// Largest Mag at or below a * b.
pub fn mul_lower(a: &MagV, b: &MagV) -> Option<MagV> {
    if matches!(a, MagV::Infinity) || matches!(b, MagV::Infinity) {
        return Some(zero());
    }
    match (parts(a), parts(b)) {
        (Some((ea, ma)), Some((eb, mb))) => norm(ma * mb, ea + eb - 128, false, false),
        _ => Some(zero()),
    }
}

/// Smallest Mag at or above a / b (Infinity when b is zero or an operand is Infinity).
pub fn div_up(a: &MagV, b: &MagV) -> Option<MagV> {
    if matches!(a, MagV::Infinity) || matches!(b, MagV::Infinity) || is_zero(b) {
        return Some(MagV::Infinity);
    }
    match (parts(a), parts(b)) {
        (Some((ea, ma)), Some((eb, mb))) => {
            let n = ma << 64;
            norm(n / mb, ea - eb - 64, n % mb != 0, true)
        }
        _ => Some(zero()),
    }
}

/// Largest Mag at or below a / b (zero when b is zero or an operand is Infinity).
pub fn div_lower(a: &MagV, b: &MagV) -> Option<MagV> {
    if matches!(a, MagV::Infinity) || matches!(b, MagV::Infinity) || is_zero(b) {
        return Some(zero());
    }
    match (parts(a), parts(b)) {
        (Some((ea, ma)), Some((eb, mb))) => {
            let n = ma << 64;
            norm(n / mb, ea - eb - 64, n % mb != 0, false)
        }
        _ => Some(zero()),
    }
}

/// Smallest Mag at or above sqrt(a).
pub fn sqrt_up(a: &MagV) -> Option<MagV> {
    if matches!(a, MagV::Infinity) {
        return Some(MagV::Infinity);
    }
    match parts(a) {
        None => Some(zero()),
        Some((ea, ma)) => {
            let (m, e) = if (ea - 128) % 2 == 0 {
                (ma << 64, ea - 128)
            } else {
                (ma << 63, ea - 127)
            };
            let r = isqrt(m);
            norm(r, e / 2, r * r != m, true)
        }
    }
}

/// Largest Mag at or below sqrt(a).
pub fn sqrt_lower(a: &MagV) -> Option<MagV> {
    if matches!(a, MagV::Infinity) {
        return Some(zero());
    }
    match parts(a) {
        None => Some(zero()),
        Some((ea, ma)) => {
            let (m, e) = if (ea - 128) % 2 == 0 {
                (ma << 64, ea - 128)
            } else {
                (ma << 63, ea - 127)
            };
            let r = isqrt(m);
            norm(r, e / 2, r * r != m, false)
        }
    }
}

/// Largest Mag at or below a - b (zero when a <= b).
pub fn sub_lower(a: &MagV, b: &MagV) -> Option<MagV> {
    if matches!(a, MagV::Infinity) || matches!(b, MagV::Infinity) {
        return Some(zero());
    }
    if order(a, b) != Ordering::Greater {
        return Some(zero());
    }
    match (parts(a), parts(b)) {
        (Some((ea, ma)), Some((eb, mb))) => {
            let gap = ea - eb;
            if gap >= 64 {
                let a = ma << 63;
                let s = gap - 63;
                let (b, frac) = if s >= 128 {
                    (0u128, true)
                } else {
                    (mb >> s, (mb & ((1u128 << s) - 1)) != 0)
                };
                let m = a - b - if frac { 1 } else { 0 };
                return norm(m, ea - 64 - 63, false, false);
            }
            norm((ma << gap) - mb, eb - 64, false, false)
        }
        (Some((e, m)), None) => Some(MagV::Finite {
            exp: i64::try_from(e).ok()?,
            mant: u64::try_from(m).ok()?,
        }),
        _ => Some(zero()),
    }
}

/// Total order of bounds: zero least, Infinity greatest and equal to itself.
pub fn order(a: &MagV, b: &MagV) -> Ordering {
    match (a, b) {
        (MagV::Infinity, MagV::Infinity) => Ordering::Equal,
        (MagV::Infinity, _) => Ordering::Greater,
        (_, MagV::Infinity) => Ordering::Less,
        (MagV::Finite { exp: ea, mant: ma }, MagV::Finite { exp: eb, mant: mb }) => {
            if *ma == 0 && *mb == 0 {
                Ordering::Equal
            } else if *ma == 0 {
                Ordering::Less
            } else if *mb == 0 {
                Ordering::Greater
            } else {
                (ea, ma).cmp(&(eb, mb))
            }
        }
    }
}
