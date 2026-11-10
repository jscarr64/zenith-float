//! Slice S2: The `Mag` to `ExactNum` bridge: exact comparisons of a `Mag` radius against zenith-float values that never
//! build the radius as an `ExactNum` unless it lies within a factor 2 of the value (CONTRACT 4.2, 10.2c, 15.2, 16.S2),
//! and the constant bridge (16.S2.6): [`constant_ball`] encloses a zenith-float constant between its two p-bit neighbours
//! from R14 `Consts::get_exact_prec`.
//!
//! Owner: slice S2 (Spence). Only that slice edits this file.

use core::cmp::Ordering;

use zenith_float_num::{ConstantId, Consts, Error, ExactNum, RoundingMode};

use crate::ball::Ball;
use crate::error::{BoundError, BoundResult};
use crate::mag::Mag;
use crate::precision::{Precision, MAX_PRECISION, MIN_PRECISION};

/// Operands whose exponents differ by at most this are subtracted exactly; wider gaps are decided by exponents first
/// (16.S2.3). Any value of at least 66 is correct (the wide-gap cutoffs need `|small| < 2^(e - 65)`); it only trades
/// exact-subtraction cost for case analysis.
const EXACT_GAP: i64 = 128;

/// The normalized `Mag` mantissa of a power of two: `Mag(TOP, k)` is `2^(k-1)`.
const TOP: u64 = 0x8000_0000_0000_0000;

/// The sign of `x` as an ordering against zero (`Less` negative, `Equal` zero, `Greater` positive).
fn sign(x: &ExactNum) -> Ordering {
    if x.is_zero() {
        Ordering::Equal
    } else if x.is_negative() {
        Ordering::Less
    } else {
        Ordering::Greater
    }
}

/// `s` if `keep`, otherwise its opposite.
fn signed(s: Ordering, keep: bool) -> Ordering {
    if keep {
        s
    } else {
        s.reverse()
    }
}

fn not_finite(x: &ExactNum) -> bool {
    x.is_nan() || x.is_inf()
}

/// The exponent `e` of `x` (`|x|` in `[2^(e-1), 2^e)`), as i64; `None` when `x` is zero or not finite.
fn exp_of(x: &ExactNum) -> Option<i64> {
    if x.is_zero() {
        return None;
    }
    x.exponent()
}

impl Mag {
    /// `self` against `|x|`, decided exactly from `x`'s top 64 mantissa bits, its sticky bit and both i64 exponents,
    /// so a radius outside zenith-float's exponent range still compares (16.S2.1). `None` when `x` is NaN or infinite.
    pub fn cmp_exact_num(&self, x: &ExactNum) -> Option<Ordering> {
        if not_finite(x) {
            return None;
        }
        if self.is_infinity() {
            return Some(Ordering::Greater);
        }
        if x.is_zero() {
            return Some(if self.is_zero() {
                Ordering::Equal
            } else {
                Ordering::Greater
            });
        }
        // Zero is stored with exponent 0, so it must not reach the exponent comparison.
        if self.is_zero() {
            return Some(Ordering::Less);
        }
        let (me, mm) = match (self.exp(), self.mant()) {
            (Some(e), Some(m)) => (e, m),
            (None, _) | (_, None) => return Some(Ordering::Less),
        };
        let (xm, sticky, xe) = x.top_u64_mantissa_and_sticky()?;
        let tail = if sticky {
            Ordering::Less
        } else {
            Ordering::Equal
        };
        Some(me.cmp(&xe).then(mm.cmp(&xm)).then(tail))
    }

    /// `self` against `|a - b|`, decided exactly (16.S2.2). Operands within [`EXACT_GAP`] binades are subtracted
    /// exactly; otherwise the exponents decide unless the radius is within a factor 2 of the larger operand, and then
    /// only the exact `|big| - self` and one short exact sum are formed (16.S2.3). `None` when an operand is NaN or
    /// infinite, or when zenith-float reports an exact operation inexact or returns NaN (an allocation failure).
    pub fn cmp_exact_diff(&self, a: &ExactNum, b: &ExactNum) -> Option<Ordering> {
        if not_finite(a) || not_finite(b) {
            return None;
        }
        if self.is_infinity() {
            return Some(Ordering::Greater);
        }
        let (ea, eb) = match (exp_of(a), exp_of(b)) {
            (Some(ea), Some(eb)) => (ea, eb),
            (None, _) => return self.cmp_exact_num(b),
            (_, None) => return self.cmp_exact_num(a),
        };
        let gap = if ea >= eb {
            ea.checked_sub(eb)?
        } else {
            eb.checked_sub(ea)?
        };
        if gap <= EXACT_GAP {
            let d = a.sub_full_prec(b);
            if d.inexact() {
                return None;
            }
            return self.cmp_exact_num(&d);
        }
        self.cmp_far_diff(a, b, ea >= eb)
    }

    /// The wide-gap case of [`Mag::cmp_exact_diff`]: `a_big` says whether `a` is the operand of larger magnitude.
    fn cmp_far_diff(&self, a: &ExactNum, b: &ExactNum, a_big: bool) -> Option<Ordering> {
        let (big, small) = if a_big { (a, b) } else { (b, a) };
        let e = exp_of(big)?;
        let es = exp_of(small)?;
        // |a - b| lies in (2^(e-1) - 2^(e-65), 2^e + 2^(e-65)) since |small| < 2^es <= 2^(e - EXACT_GAP); a radius below
        // 2^(e-1) has mant <= u64::MAX so it is at most 2^(e-1) - 2^(e-65).
        // |a - b| > 0 here; zero is stored with exponent 0, so it must not reach the exponent cutoffs.
        if self.is_zero() {
            return Some(Ordering::Less);
        }
        let me = self.exp()?;
        if me <= e.checked_sub(1)? {
            return Some(Ordering::Less);
        }
        if me >= e.checked_add(2)? {
            return Some(Ordering::Greater);
        }
        // a - b = tb + ts with tb = ±big, ts = ±small; |a - b| = sb (a - b) where sb is the sign of tb.
        let big_neg = big.is_negative();
        let tb_pos = big_neg != a_big;
        let u_keep = a_big;
        // u = tb - sb m, kept as u = (u_pos ? +1 : -1) * ur with ur = big - m (big >= 0) or big + m (big < 0).
        let ur = if me == e.checked_add(1)? {
            // self in [2^e, 2^(e+1)), possibly above EXPONENT_MAX (16.S2.1), so it is never converted here. A mantissa
            // above TOP gives self >= 2^e + 2^(e-63) > 2^e + 2^(e-129) > |a - b|.
            if self.mant()? != TOP {
                return Some(Ordering::Greater);
            }
            // self = 2^e = h + h with h = 2^(e-1) = Mag(TOP, e), in range since big has exponent e.
            let h = Mag::finite(TOP, e).ok()?.to_exact_num().ok()?;
            let once = if big_neg {
                big.add_full_prec(&h)
            } else {
                big.sub_full_prec(&h)
            };
            if once.inexact() || once.is_nan() {
                return None;
            }
            if big_neg {
                once.add_full_prec(&h)
            } else {
                once.sub_full_prec(&h)
            }
        } else {
            // me == e: self < 2^e, in range since big has exponent e.
            let m = self.to_exact_num().ok()?;
            if big_neg {
                big.add_full_prec(&m)
            } else {
                big.sub_full_prec(&m)
            }
        };
        if ur.inexact() || ur.is_nan() {
            return None;
        }
        // Sign of ts: small enters with + when it is a, with - when it is b.
        let ts_pos = small.is_negative() == a_big;
        let ts_sign = if ts_pos {
            Ordering::Greater
        } else {
            Ordering::Less
        };
        // y = u + ts; sign(|a - b| - self) = sb * sign(y).
        let y_sign = match exp_of(&ur) {
            None => ts_sign,
            Some(eu) => {
                let u_sign = signed(sign(&ur), u_keep);
                if es <= eu.checked_sub(2)? {
                    u_sign
                } else if eu <= es.checked_sub(2)? {
                    ts_sign
                } else {
                    // y = fu ur + fs small with fu = u_keep, fs = (sign of small as ±) == ts_pos.
                    let small_pos = !small.is_negative();
                    let fs_plus = small_pos == ts_pos;
                    let z = if u_keep == fs_plus {
                        ur.add_full_prec(small)
                    } else {
                        ur.sub_full_prec(small)
                    };
                    if z.inexact() || z.is_nan() {
                        return None;
                    }
                    signed(sign(&z), u_keep)
                }
            }
        };
        Some(signed(y_sign, tb_pos).reverse())
    }
}

/// A constant [`constant_ball`] encloses: the seven of zenith-float's `ConstantId` (16.S2.6).
#[derive(Debug, PartialEq, Eq)]
pub enum Constant {
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

impl Constant {
    /// Every constant, in `ConstantId` order.
    pub const ALL: [Constant; 7] = [
        Constant::Pi,
        Constant::E,
        Constant::Ln2,
        Constant::Ln10,
        Constant::Sqrt2,
        Constant::Phi,
        Constant::EulerGamma,
    ];

    /// The zenith-float id of `self` (a fresh value: `ConstantId` is not `Clone`).
    fn zf_id(&self) -> ConstantId {
        match self {
            Constant::Pi => ConstantId::Pi,
            Constant::E => ConstantId::E,
            Constant::Ln2 => ConstantId::Ln2,
            Constant::Ln10 => ConstantId::Ln10,
            Constant::Sqrt2 => ConstantId::Sqrt2,
            Constant::Phi => ConstantId::Phi,
            Constant::EulerGamma => ConstantId::EulerGamma,
        }
    }
}

/// A ball at precision `prec` that contains the constant `c` (16.S2.6): its ends are exactly the p-bit neighbours
/// `down < c < up` that `Consts::get_exact_prec` returns with `Down` and `Up`, so its width is one ulp at `prec`. The
/// midpoint `down + ulp/2` has `prec + 1` bits (D15a keeps it exact) and the radius is `ulp/2`. The result depends only
/// on `c` and `prec`, never on what `cc` has cached.
///
/// # Errors
/// `InvalidArgument` when `prec` is outside `MIN_PRECISION..=MAX_PRECISION`; `Resource` when zenith-float cannot settle
/// the constant (`PrecisionRetryExhausted`, `MemoryAllocation`) or its two neighbours are not one ulp apart.
pub fn constant_ball(cc: &mut Consts, c: &Constant, prec: Precision) -> BoundResult<Ball> {
    const OP: &str = "bridge::constant_ball";
    if !(MIN_PRECISION..=MAX_PRECISION).contains(&prec) {
        return Err(BoundError::InvalidArgument {
            op: OP,
            detail: "precision outside MIN_PRECISION..=MAX_PRECISION",
        });
    }
    let down = cc
        .get_exact_prec(c.zf_id(), prec, RoundingMode::Down)
        .map_err(constant_error)?;
    let up = cc
        .get_exact_prec(c.zf_id(), prec, RoundingMode::Up)
        .map_err(constant_error)?;
    let not_one_ulp = BoundError::Resource {
        op: OP,
        detail: "constant neighbours are not one ulp apart",
    };
    // up - down is one ulp of down: a power of two, exact at prec bits.
    let ulp = up.sub(&down, prec, RoundingMode::Up);
    if ulp.inexact() || !ulp.is_positive() || ulp.is_zero() {
        return Err(not_one_ulp);
    }
    let half = ulp.ldexp(-1, prec, RoundingMode::Up);
    if half.inexact() || !half.is_positive() {
        return Err(not_one_ulp);
    }
    let wide = prec.checked_add(1).ok_or(BoundError::InvalidArgument {
        op: OP,
        detail: "precision + 1 overflows",
    })?;
    let mid = down.add(&half, wide, RoundingMode::Down);
    if mid.inexact() || !mid.is_positive() {
        return Err(not_one_ulp);
    }
    let rad = Mag::from_exact_num(&half);
    Ball::from_mid_rad(mid, rad, prec)
}

/// zf errors of `get_exact_prec` and the bracket arithmetic in zb terms (CONTRACT O2 table).
fn constant_error(e: Error) -> BoundError {
    const OP: &str = "bridge::constant_ball";
    match e {
        Error::InvalidArgument => BoundError::InvalidArgument {
            op: OP,
            detail: "zf invalid argument",
        },
        Error::ExponentOverflow(_) => BoundError::Representation {
            op: OP,
            detail: "zf exponent overflow",
        },
        Error::DivisionByZero => BoundError::Domain {
            op: OP,
            detail: "zf division by zero",
        },
        Error::PrecisionRetryExhausted => BoundError::Resource {
            op: OP,
            detail: "zf precision retry exhausted",
        },
        Error::MemoryAllocation => BoundError::Resource {
            op: OP,
            detail: "zf memory allocation",
        },
    }
}
