//! L0a-5 Mag property checkers for the wide-range fuzz (TEST_PLAN 8.1 and section 5; written by Spence, the module scaffold is
//! tests/harness/fuzz_mag.rs). ORACLE: the true result is compared with a candidate `Mag` by exact integer arithmetic only (u128 and i128,
//! no model, no big-integer crate, no use of the library): a `Mag` is `mant * 2^(exp - 64)` with `mant >= 2^63`, and every operation's exact
//! result is one of  `P * 2^s`  (mul, div, sqrt through its square, add and sub with a gap <= 63) or a value squeezed strictly between two
//! neighbouring `Mag`s by a bounded argument (add and sub with a gap >= 64: the smaller operand is below one unit of the larger).
//! The oracle itself is validated against the mpmath goldens (tests/goldens/mag_ops.txt), which it must accept row by row and whose perturbed
//! results it must reject.
//! PROPERTIES (each a `*_case` function over one deterministic draw, run by the scaffold):
//!   soundness and tightness   `up >= exact` and `pred(up) < exact`; `down <= exact` and `succ(down) > exact` (so up is the smallest upper bound,
//!                             down the largest lower bound: within one Mag-ulp of each other), with the range clamps of CONTRACT 3.2 to 3.4, 14.2
//!   bracket                   `up` is `down` or its successor (mul, div, sqrt)
//!   monotonicity              add, mul, sqrt, div (non-increasing in the divisor), sub_lower (non-increasing in the subtrahend), up and lower forms
//!   commutativity             add_up, mul_up, mul_lower
//!   identities                x+0, x*1, x/1, x/x, x-0, x-x, 0-x, 0/x, x/0, x+x, x*2, sqrt of exact squares and of odd-exponent powers of two
//!   canonical form            every result of every operation, with zero and Infinity operands mixed in
//!   ulp helpers               ulp_mag, half_ulp (exact powers of two with the range clamps), smaller_than_half_ulp (exact comparison), their relations

use crate::common::generators::{
    draw_exponent, draw_mantissa, draw_precision, forced_prop_exponents, ExpDomain, PrecSource, TOP,
};
use crate::common::rng::SplitMix64;
use core::cmp::Ordering;
use zenith_bound::mag::Mag;

// ---- the operations under test (swappable, so the self-tests can plant a wrong one) -----------------------------------------

pub struct Ops {
    pub add_up: fn(&Mag, &Mag) -> Mag,
    pub sub_lower: fn(&Mag, &Mag) -> Mag,
    pub mul_up: fn(&Mag, &Mag) -> Mag,
    pub mul_lower: fn(&Mag, &Mag) -> Mag,
    pub div_up: fn(&Mag, &Mag) -> Mag,
    pub div_lower: fn(&Mag, &Mag) -> Mag,
    pub sqrt_up: fn(&Mag) -> Mag,
    pub sqrt_lower: fn(&Mag) -> Mag,
    pub ulp_mag: fn(i64, usize) -> Mag,
    pub half_ulp: fn(i64, usize) -> Mag,
    pub smaller_than_half_ulp: fn(&Mag, i64, usize) -> bool,
}

/// The real library.
pub fn real() -> Ops {
    Ops {
        add_up: Mag::add_up,
        sub_lower: Mag::sub_lower,
        mul_up: Mag::mul_up,
        mul_lower: Mag::mul_lower,
        div_up: Mag::div_up,
        div_lower: Mag::div_lower,
        sqrt_up: Mag::sqrt_up,
        sqrt_lower: Mag::sqrt_lower,
        ulp_mag: Mag::ulp_mag,
        half_ulp: Mag::half_ulp,
        smaller_than_half_ulp: Mag::smaller_than_half_ulp,
    }
}

// ---- results as plain data ----------------------------------------------------------------------------------------------------

#[derive(Debug, PartialEq, Eq)]
pub enum Out {
    Zero,
    Infinity,
    Finite { exp: i64, mant: u64 },
}

/// The result as data. `Err` unless it is canonical: zero is `{0, 0}`, a finite value has the top mantissa bit set, the accessors agree.
pub fn out_of(m: &Mag) -> Result<Out, String> {
    let out = match (m.exp(), m.mant()) {
        (Some(0), Some(0)) => Out::Zero,
        (Some(exp), Some(mant)) if mant >= TOP => Out::Finite { exp, mant },
        (None, None) => Out::Infinity,
        (e, mm) => {
            return Err(format!(
                "non-canonical result: exp() {:?}, mant() {:?}",
                e, mm
            ))
        }
    };
    if m.is_zero() != (out == Out::Zero) || m.is_infinity() != (out == Out::Infinity) {
        return Err(format!(
            "is_zero()/is_infinity() disagree with the fields: {:?}",
            out
        ));
    }
    Ok(out)
}

fn rank(o: &Out) -> u8 {
    match o {
        Out::Zero => 0,
        Out::Finite { .. } => 1,
        Out::Infinity => 2,
    }
}

/// The order of values: zero < every finite value (by exponent, then mantissa) < Infinity.
pub fn out_cmp(a: &Out, b: &Out) -> Ordering {
    match (a, b) {
        (
            Out::Finite {
                exp: ea, mant: ma, ..
            },
            Out::Finite {
                exp: eb, mant: mb, ..
            },
        ) => (ea, ma).cmp(&(eb, mb)),
        _ => rank(a).cmp(&rank(b)),
    }
}

/// The next larger `Mag` value after the finite `(exp, mant)`; `None` means Infinity.
pub fn next_up(exp: i64, mant: u64) -> Option<(i64, u64)> {
    match mant.checked_add(1) {
        Some(m) => Some((exp, m)),
        None => exp.checked_add(1).map(|e| (e, TOP)),
    }
}

/// The next smaller `Mag` value before the finite `(exp, mant)`; `None` means zero.
pub fn next_down(exp: i64, mant: u64) -> Option<(i64, u64)> {
    if mant > TOP {
        Some((exp, mant - 1))
    } else {
        exp.checked_sub(1).map(|e| (e, u64::MAX))
    }
}

// ---- exact comparison ---------------------------------------------------------------------------------------------------------

fn bit_len(x: u128) -> i128 {
    i128::from(128 - x.leading_zeros())
}

fn shl(x: u128, k: i128) -> Result<u128, String> {
    u32::try_from(k)
        .ok()
        .and_then(|k| x.checked_shl(k))
        .ok_or_else(|| format!("shift of {:#x} by {} is out of range", x, k))
}

/// Orders `x * 2^sx` against `y * 2^sy` (both nonnegative).
pub fn cmp_scaled(x: u128, sx: i128, y: u128, sy: i128) -> Result<Ordering, String> {
    match (x == 0, y == 0) {
        (true, true) => return Ok(Ordering::Equal),
        (true, false) => return Ok(Ordering::Less),
        (false, true) => return Ok(Ordering::Greater),
        (false, false) => {}
    }
    let (lx, ly) = (bit_len(x), bit_len(y));
    let (tx, ty) = (sx + lx, sy + ly);
    if tx != ty {
        return Ok(tx.cmp(&ty));
    }
    // same position of the top bit: the scales differ by the difference of the bit lengths (at most 127)
    if lx <= ly {
        Ok(shl(x, ly - lx)?.cmp(&y))
    } else {
        Ok(x.cmp(&shl(y, lx - ly)?))
    }
}

fn wide(m: u64) -> u128 {
    u128::from(m)
}

fn ex(e: i64) -> i128 {
    i128::from(e)
}

/// Order of the exact `a * b` against the candidate `(er, mr)`.
pub fn cmp_mul(
    (ea, ma): (i64, u64),
    (eb, mb): (i64, u64),
    (er, mr): (i64, u64),
) -> Result<Ordering, String> {
    cmp_scaled(
        wide(ma) * wide(mb),
        ex(ea) + ex(eb) - 128,
        wide(mr),
        ex(er) - 64,
    )
}

/// Order of the exact `a / b` against the candidate: `ma * 2^(ea - eb)` against `mr * mb * 2^(er - 64)`.
pub fn cmp_div(
    (ea, ma): (i64, u64),
    (eb, mb): (i64, u64),
    (er, mr): (i64, u64),
) -> Result<Ordering, String> {
    cmp_scaled(wide(ma), ex(ea) - ex(eb), wide(mr) * wide(mb), ex(er) - 64)
}

/// Order of the exact `sqrt(a)` against the candidate: `a` against the candidate's square.
pub fn cmp_sqrt((ea, ma): (i64, u64), (er, mr): (i64, u64)) -> Result<Ordering, String> {
    cmp_scaled(wide(ma), ex(ea) - 64, wide(mr) * wide(mr), 2 * ex(er) - 128)
}

/// Order of the exact `a + b` against the candidate. Gap <= 63: exact in a u128. Gap >= 64: the sum lies strictly between the larger
/// operand and its successor (the smaller one is below one unit of the larger), so it is above every candidate <= the larger one and below
/// every other.
pub fn cmp_add(
    (ea, ma): (i64, u64),
    (eb, mb): (i64, u64),
    (er, mr): (i64, u64),
) -> Result<Ordering, String> {
    let ((eh, mh), (el, ml)) = if (ea, ma) >= (eb, mb) {
        ((ea, ma), (eb, mb))
    } else {
        ((eb, mb), (ea, ma))
    };
    let gap = ex(eh) - ex(el);
    if gap <= 63 {
        let n = shl(wide(mh), gap)? + wide(ml);
        cmp_scaled(n, ex(el) - 64, wide(mr), ex(er) - 64)
    } else if (er, mr) <= (eh, mh) {
        Ok(Ordering::Greater)
    } else {
        Ok(Ordering::Less)
    }
}

/// Order of the exact `a - b` (requires `a > b` as values) against the candidate. Gap <= 63: exact. Gap >= 64: the difference is `a`
/// minus something in `(0, 2^(ea - 64))`, so it is below `a`, above every candidate two or more binades below `a`, and for the two
/// binades next to `a` it is decided by comparing the candidate's distance from `a` with the subtrahend.
pub fn cmp_sub(
    (ea, ma): (i64, u64),
    (eb, mb): (i64, u64),
    (er, mr): (i64, u64),
) -> Result<Ordering, String> {
    let gap = ex(ea) - ex(eb);
    if gap <= 63 {
        let n = shl(wide(ma), gap)? - wide(mb);
        return cmp_scaled(n, ex(eb) - 64, wide(mr), ex(er) - 64);
    }
    if (er, mr) >= (ea, ma) {
        return Ok(Ordering::Less);
    }
    let down = ex(ea) - ex(er);
    if down >= 2 {
        return Ok(Ordering::Greater);
    }
    let a2 = wide(ma) << 1;
    let r2 = if down == 0 { wide(mr) << 1 } else { wide(mr) };
    cmp_scaled(a2 - r2, ex(ea) - 65, wide(mb), ex(eb) - 64)
}

// ---- the checkers of one result ------------------------------------------------------------------------------------------------

pub type Cmp<'a> = dyn Fn(i64, u64) -> Result<Ordering, String> + 'a;

const MAX_EXP: i64 = i64::MAX;
const MIN_EXP: i64 = i64::MIN;

/// A round-up result must be the smallest `Mag` that is `>= exact`: `Infinity` only above the largest Mag, the smallest positive Mag
/// when the exact value is below it (CONTRACT 3.2), else `r >= exact > pred(r)`. `cmp(e, m)` is the order of the exact value against `(e, m)`.
pub fn check_up(out: &Out, cmp: &Cmp) -> Result<(), String> {
    match out {
        Out::Zero => Err(String::from(
            "round-up result is zero for a positive exact value",
        )),
        Out::Infinity => {
            if cmp(MAX_EXP, u64::MAX)? == Ordering::Greater {
                Ok(())
            } else {
                Err(String::from(
                    "Infinity although the exact value fits in the largest Mag",
                ))
            }
        }
        Out::Finite { exp, mant } => {
            if cmp(*exp, *mant)? == Ordering::Greater {
                return Err(String::from("result is below the exact value"));
            }
            if *exp == MIN_EXP && *mant == TOP {
                return Ok(());
            }
            match next_down(*exp, *mant) {
                Some((pe, pm)) if cmp(pe, pm)? != Ordering::Greater => Err(String::from(
                    "not the smallest upper bound: the next smaller Mag is also >= the exact value",
                )),
                _ => Ok(()),
            }
        }
    }
}

/// A round-down result must be the largest `Mag` that is `<= exact`: zero only below the smallest positive Mag, the largest Mag when the
/// exact value is above it (CONTRACT 14.2), else `r <= exact < succ(r)`.
pub fn check_down(out: &Out, cmp: &Cmp) -> Result<(), String> {
    match out {
        Out::Infinity => Err(String::from(
            "round-down result is Infinity for finite operands",
        )),
        Out::Zero => {
            if cmp(MIN_EXP, TOP)? == Ordering::Less {
                Ok(())
            } else {
                Err(String::from(
                    "zero although the exact value reaches the smallest positive Mag",
                ))
            }
        }
        Out::Finite { exp, mant } => {
            if cmp(*exp, *mant)? == Ordering::Less {
                return Err(String::from("result is above the exact value"));
            }
            if *exp == MAX_EXP && *mant == u64::MAX {
                return Ok(());
            }
            match next_up(*exp, *mant) {
                Some((se, sm)) if cmp(se, sm)? != Ordering::Less => Err(String::from(
                    "not the largest lower bound: the next larger Mag is also <= the exact value",
                )),
                _ => Ok(()),
            }
        }
    }
}

/// `up` is `down` or the next larger value (the two bounds of one exact result are at most one Mag-ulp apart).
pub fn check_bracket(up: &Out, down: &Out) -> Result<(), String> {
    if up == down {
        return Ok(());
    }
    let ok = match down {
        Out::Zero => matches!(up, Out::Finite { exp, mant } if *exp == MIN_EXP && *mant == TOP),
        Out::Infinity => false,
        Out::Finite { exp, mant } => match next_up(*exp, *mant) {
            Some((e, m)) => matches!(up, Out::Finite { exp: ue, mant: um } if *ue == e && *um == m),
            None => *up == Out::Infinity,
        },
    };
    if ok {
        Ok(())
    } else {
        Err(format!(
            "round-up {:?} is neither the round-down {:?} nor its successor",
            up, down
        ))
    }
}

// ---- cases ----------------------------------------------------------------------------------------------------------------------

/// What a failed case reports: the operation, the precision of the draw and the operands as text (the replay line fields).
pub struct Bad {
    pub op: &'static str,
    pub prec: usize,
    pub a: String,
    pub b: String,
    pub msg: String,
}

/// One draw: two finite nonzero operands `0.m * 2^e` (`mant` has the top bit set).
pub struct Case {
    pub prec: usize,
    pub ea: i64,
    pub ma: u64,
    pub eb: i64,
    pub mb: u64,
}

/// How the second exponent relates to the first.
pub enum Rel {
    /// Free: the generator's mix, anchored at the first exponent (so the controlled-gap band applies).
    Any,
    /// Half of the cases: `ea + eb` within 3 of one of `0`, `i64::MAX`, `i64::MIN`, `64`, `-64` (mul at the range boundaries).
    Sum,
    /// Half of the cases: `ea - eb` within 3 of one of the same targets (div at the range boundaries).
    Diff,
}

pub fn enc(exp: i64, mant: u64) -> String {
    format!("{:016x}e{}", mant, exp)
}

fn bad(op: &'static str, c: &Case, msg: String) -> Bad {
    Bad {
        op,
        prec: c.prec,
        a: enc(c.ea, c.ma),
        b: enc(c.eb, c.mb),
        msg,
    }
}

fn bare(msg: String) -> Bad {
    Bad {
        op: "draw",
        prec: 0,
        a: String::from("-"),
        b: String::from("-"),
        msg,
    }
}

fn exp_of(src: &mut Src, prec: u32, anchor: i64) -> Result<i64, String> {
    match src.next_forced() {
        Some(e) => Ok(e),
        None => Ok(draw_exponent(src.rng, &ExpDomain::MAG, anchor, prec)?.value),
    }
}

/// Draws a case. The consumption is fixed (frozen draw order): precision, ea, ma, shape, free eb, delta, target, mb.
pub fn draw_case(src: &mut Src, rel: &Rel) -> Result<Case, String> {
    let prec = draw_precision(src.rng, &PrecSource::Fast)?;
    let ea = exp_of(src, prec, 0)?;
    let ma = draw_mantissa(src.rng, prec)?.value;
    let shape = src.below(4);
    let free = exp_of(src, prec, ea)?;
    let delta = i128::from(src.below(7)) - 3;
    let target: i128 = match src.below(5) {
        0 => 0,
        1 => i128::from(i64::MAX),
        2 => i128::from(i64::MIN),
        3 => 64,
        _ => -64,
    };
    let mb = draw_mantissa(src.rng, prec)?.value;
    let related = match rel {
        Rel::Any => None,
        Rel::Sum => Some(target - ex(ea) + delta),
        Rel::Diff => Some(ex(ea) - target - delta),
    };
    // a forced case keeps its forced second exponent (the forced pair is the edge being exercised)
    let eb = match (shape >= 2 && !src.is_forced(), related) {
        (true, Some(v)) => i64::try_from(v).unwrap_or(free),
        _ => free,
    };
    Ok(Case {
        prec: api_prec(prec),
        ea,
        ma,
        eb,
        mb,
    })
}

// ---- where a case's draws come from ---------------------------------------------------------------------------------------------

/// The draws of one case: the random stream, with the case's first exponent draws (in draw order) optionally forced to fixed values.
pub struct Src<'a> {
    rng: &'a mut SplitMix64,
    forced: &'a [i64],
    used: usize,
}

impl<'a> Src<'a> {
    /// Every exponent from the random stream (the random band).
    pub fn random(rng: &'a mut SplitMix64) -> Self {
        Src {
            rng,
            forced: &[],
            used: 0,
        }
    }
    /// The first `forced.len()` exponent draws are `forced`, in order; later ones come from the random stream.
    pub fn forced(rng: &'a mut SplitMix64, forced: &'a [i64]) -> Self {
        Src {
            rng,
            forced,
            used: 0,
        }
    }
    fn is_forced(&self) -> bool {
        !self.forced.is_empty()
    }
    fn next_forced(&mut self) -> Option<i64> {
        let e = self.forced.get(self.used)?;
        self.used += 1;
        Some(*e)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.rng.below(n)
    }
}

/// The forced pairs every property runs on, besides its random band: every ordered pair of `forced_prop_exponents()` (13^2 = 169).
pub const FORCED_PAIRS: u64 = 169;

/// Runs `f` once on every ordered pair `(ea, eb)` of `forced_prop_exponents()`, as the first two exponent draws (mantissas, precisions and
/// shapes from a stream seeded by the pair index). Returns the number of pairs run; `Err` is `forced pair <i> (<ea>, <eb>) op=.. prec=.. a=.. b=.. -- why`.
pub fn run_forced(f: &dyn Fn(&mut Src) -> Result<u8, Bad>) -> Result<u64, String> {
    let es = forced_prop_exponents();
    let mut n = 0u64;
    for ea in es.iter() {
        for eb in es.iter() {
            let pair = [*ea, *eb];
            let mut rng = SplitMix64::new(n);
            if let Err(b) = f(&mut Src::forced(&mut rng, &pair)) {
                return Err(format!(
                    "forced pair {} ({}, {}) op={} prec={} a={} b={} -- {}",
                    n, ea, eb, b.op, b.prec, b.a, b.b, b.msg
                ));
            }
            n += 1;
        }
    }
    Ok(n)
}

fn mag_from(mant: u64, exp: i64) -> Result<Mag, String> {
    Mag::finite(mant, exp).map_err(|e| format!("Mag::finite({:#x}, {}) failed: {:?}", mant, exp, e))
}

fn ap2(f: fn(&Mag, &Mag) -> Mag, a: &Mag, b: &Mag) -> Result<Out, String> {
    out_of(&f(a, b))
}

fn ap1(f: fn(&Mag) -> Mag, a: &Mag) -> Result<Out, String> {
    out_of(&f(a))
}

/// Outcome classes of a case, for the non-vacuity counts of the scaffold.
pub const FIN: u8 = 1;
pub const INF: u8 = 2;
pub const ZERO: u8 = 4;
pub const EXACT: u8 = 8;

fn flag(o: &Out) -> u8 {
    match o {
        Out::Finite { .. } => FIN,
        Out::Infinity => INF,
        Out::Zero => ZERO,
    }
}

/// `add_up`: sound, smallest upper bound.
pub fn add_up_case(ops: &Ops, src: &mut Src) -> Result<u8, Bad> {
    let c = draw_case(src, &Rel::Any).map_err(bare)?;
    let run = || -> Result<u8, String> {
        let (a, b) = (mag_from(c.ma, c.ea)?, mag_from(c.mb, c.eb)?);
        let out = ap2(ops.add_up, &a, &b)?;
        check_up(&out, &|er, mr| {
            cmp_add((c.ea, c.ma), (c.eb, c.mb), (er, mr))
        })?;
        Ok(flag(&out))
    };
    run().map_err(|m| bad("add_up", &c, m))
}

/// `sub_lower`: zero iff `a <= b`, else the largest lower bound.
pub fn sub_lower_case(ops: &Ops, src: &mut Src) -> Result<u8, Bad> {
    let c = draw_case(src, &Rel::Any).map_err(bare)?;
    let run = || -> Result<u8, String> {
        let (a, b) = (mag_from(c.ma, c.ea)?, mag_from(c.mb, c.eb)?);
        let out = ap2(ops.sub_lower, &a, &b)?;
        if (c.ea, c.ma) <= (c.eb, c.mb) {
            if out != Out::Zero {
                return Err(format!("a <= b but the result is {:?}", out));
            }
        } else {
            check_down(&out, &|er, mr| {
                cmp_sub((c.ea, c.ma), (c.eb, c.mb), (er, mr))
            })?;
        }
        Ok(flag(&out))
    };
    run().map_err(|m| bad("sub_lower", &c, m))
}

/// `mul_up` and `mul_lower`: tight bounds, bracketing each other.
pub fn mul_case(ops: &Ops, src: &mut Src) -> Result<u8, Bad> {
    let c = draw_case(src, &Rel::Sum).map_err(bare)?;
    let run = || -> Result<u8, (&'static str, String)> {
        let (a, b) = (
            mag_from(c.ma, c.ea).map_err(|m| ("mul", m))?,
            mag_from(c.mb, c.eb).map_err(|m| ("mul", m))?,
        );
        let cmp = |er, mr| cmp_mul((c.ea, c.ma), (c.eb, c.mb), (er, mr));
        let up = ap2(ops.mul_up, &a, &b).map_err(|m| ("mul_up", m))?;
        check_up(&up, &cmp).map_err(|m| ("mul_up", m))?;
        let lo = ap2(ops.mul_lower, &a, &b).map_err(|m| ("mul_lower", m))?;
        check_down(&lo, &cmp).map_err(|m| ("mul_lower", m))?;
        check_bracket(&up, &lo).map_err(|m| ("mul", m))?;
        Ok(flag(&up) | if up == lo { EXACT } else { 0 } | if lo == Out::Zero { ZERO } else { 0 })
    };
    run().map_err(|(op, m)| bad(op, &c, m))
}

/// `div_up` and `div_lower`.
pub fn div_case(ops: &Ops, src: &mut Src) -> Result<u8, Bad> {
    let c = draw_case(src, &Rel::Diff).map_err(bare)?;
    let run = || -> Result<u8, (&'static str, String)> {
        let (a, b) = (
            mag_from(c.ma, c.ea).map_err(|m| ("div", m))?,
            mag_from(c.mb, c.eb).map_err(|m| ("div", m))?,
        );
        let cmp = |er, mr| cmp_div((c.ea, c.ma), (c.eb, c.mb), (er, mr));
        let up = ap2(ops.div_up, &a, &b).map_err(|m| ("div_up", m))?;
        check_up(&up, &cmp).map_err(|m| ("div_up", m))?;
        let lo = ap2(ops.div_lower, &a, &b).map_err(|m| ("div_lower", m))?;
        check_down(&lo, &cmp).map_err(|m| ("div_lower", m))?;
        check_bracket(&up, &lo).map_err(|m| ("div", m))?;
        Ok(flag(&up) | if up == lo { EXACT } else { 0 } | if lo == Out::Zero { ZERO } else { 0 })
    };
    run().map_err(|(op, m)| bad(op, &c, m))
}

/// `sqrt_up` and `sqrt_lower` (unary: the second operand of the case is unused and printed as `-`).
pub fn sqrt_case(ops: &Ops, src: &mut Src) -> Result<u8, Bad> {
    let c = draw_case(src, &Rel::Any).map_err(bare)?;
    let run = || -> Result<u8, (&'static str, String)> {
        let a = mag_from(c.ma, c.ea).map_err(|m| ("sqrt", m))?;
        let cmp = |er, mr| cmp_sqrt((c.ea, c.ma), (er, mr));
        let up = out_of(&(ops.sqrt_up)(&a)).map_err(|m| ("sqrt_up", m))?;
        check_up(&up, &cmp).map_err(|m| ("sqrt_up", m))?;
        let lo = out_of(&(ops.sqrt_lower)(&a)).map_err(|m| ("sqrt_lower", m))?;
        check_down(&lo, &cmp).map_err(|m| ("sqrt_lower", m))?;
        check_bracket(&up, &lo).map_err(|m| ("sqrt", m))?;
        Ok(flag(&up) | if up == lo { EXACT } else { 0 })
    };
    run().map_err(|(op, m)| bad(op, &c, m)).map_err(unary)
}

fn unary(mut b: Bad) -> Bad {
    b.b = String::from("-");
    b
}

// ---- commutativity ------------------------------------------------------------------------------------------------------------

/// `add_up`, `mul_up`, `mul_lower` give the same result in both argument orders.
pub fn commutative_case(ops: &Ops, src: &mut Src) -> Result<u8, Bad> {
    let c = draw_case(src, &Rel::Sum).map_err(bare)?;
    let run = || -> Result<u8, (&'static str, String)> {
        let (a, b) = (
            mag_from(c.ma, c.ea).map_err(|m| ("commute", m))?,
            mag_from(c.mb, c.eb).map_err(|m| ("commute", m))?,
        );
        for (name, f) in [
            ("add_up", ops.add_up),
            ("mul_up", ops.mul_up),
            ("mul_lower", ops.mul_lower),
        ] {
            let ab = ap2(f, &a, &b).map_err(|m| (name, m))?;
            let ba = ap2(f, &b, &a).map_err(|m| (name, m))?;
            if ab != ba {
                return Err((name, format!("f(a, b) = {:?} but f(b, a) = {:?}", ab, ba)));
            }
        }
        Ok(FIN)
    };
    run().map_err(|(op, m)| bad(op, &c, m))
}

// ---- monotonicity -------------------------------------------------------------------------------------------------------------

/// Which operation and which direction in each argument.
pub enum Mono {
    AddUp,
    MulUp,
    MulLower,
    DivUp,
    DivLower,
    SqrtUp,
    SqrtLower,
    SubLower,
}

impl Mono {
    fn name(&self) -> &'static str {
        match self {
            Mono::AddUp => "add_up",
            Mono::MulUp => "mul_up",
            Mono::MulLower => "mul_lower",
            Mono::DivUp => "div_up",
            Mono::DivLower => "div_lower",
            Mono::SqrtUp => "sqrt_up",
            Mono::SqrtLower => "sqrt_lower",
            Mono::SubLower => "sub_lower",
        }
    }
}

/// One operand as (exponent, mantissa).
type Raw = (i64, u64);

/// Draws a second operand that is at or above `(e, m)`: the next larger value, a fresh draw ordered against it, or the same value.
fn draw_partner(src: &mut Src, prec: u32, e: i64, m: u64) -> Result<(Raw, Raw), String> {
    let pick = src.below(4);
    let fresh_e = exp_of(src, prec, e)?;
    let fresh_m = draw_mantissa(src.rng, prec)?.value;
    let other = match (pick, next_up(e, m)) {
        (0, Some(n)) | (1, Some(n)) => n,
        (2, _) => (fresh_e, fresh_m),
        _ => (e, m),
    };
    if (e, m) <= other {
        Ok(((e, m), other))
    } else {
        Ok((other, (e, m)))
    }
}

fn le(label: &str, lo: &Out, hi: &Out) -> Result<(), String> {
    if out_cmp(lo, hi) == Ordering::Greater {
        Err(format!("{}: {:?} is above {:?}", label, lo, hi))
    } else {
        Ok(())
    }
}

/// `f(x1, y) <= f(x2, y)` for `x1 <= x2` (and the divisor/subtrahend reversed): the monotonicity of CONTRACT 14.1.
pub fn mono_case(ops: &Ops, which: &Mono, src: &mut Src) -> Result<u8, Bad> {
    let rel = match which {
        Mono::MulUp | Mono::MulLower => Rel::Sum,
        Mono::DivUp | Mono::DivLower => Rel::Diff,
        Mono::AddUp | Mono::SqrtUp | Mono::SqrtLower | Mono::SubLower => Rel::Any,
    };
    let c = draw_case(src, &rel).map_err(bare)?;
    let mut run = || -> Result<u8, String> {
        let ((alo_e, alo_m), (ahi_e, ahi_m)) = draw_partner(src, gen_prec(c.prec)?, c.ea, c.ma)?;
        let ((blo_e, blo_m), (bhi_e, bhi_m)) = draw_partner(src, gen_prec(c.prec)?, c.eb, c.mb)?;
        let (alo, ahi) = (mag_from(alo_m, alo_e)?, mag_from(ahi_m, ahi_e)?);
        let (blo, bhi) = (mag_from(blo_m, blo_e)?, mag_from(bhi_m, bhi_e)?);
        // (binary op, second argument non-increasing) or the unary sqrt
        let (f, decreasing) = match which {
            Mono::SqrtUp => {
                return le(
                    "sqrt in the argument",
                    &ap1(ops.sqrt_up, &alo)?,
                    &ap1(ops.sqrt_up, &ahi)?,
                )
                .map(|()| FIN)
            }
            Mono::SqrtLower => {
                return le(
                    "sqrt in the argument",
                    &ap1(ops.sqrt_lower, &alo)?,
                    &ap1(ops.sqrt_lower, &ahi)?,
                )
                .map(|()| FIN)
            }
            Mono::AddUp => (ops.add_up, false),
            Mono::MulUp => (ops.mul_up, false),
            Mono::MulLower => (ops.mul_lower, false),
            Mono::DivUp => (ops.div_up, true),
            Mono::DivLower => (ops.div_lower, true),
            Mono::SubLower => (ops.sub_lower, true),
        };
        // first argument: nondecreasing for every binary operation
        le("first argument", &ap2(f, &alo, &blo)?, &ap2(f, &ahi, &blo)?)?;
        // second argument: nondecreasing for add and mul, non-increasing for div and sub
        let (x, y) = (ap2(f, &alo, &blo)?, ap2(f, &alo, &bhi)?);
        if decreasing {
            le("second argument (non-increasing)", &y, &x)?;
        } else {
            le("second argument", &x, &y)?;
        }
        Ok(FIN)
    };
    run().map_err(|m| bad(which.name(), &c, m))
}

// ---- identities ---------------------------------------------------------------------------------------------------------------

fn require_eq(label: &str, got: &Out, want: &Out) -> Result<(), String> {
    if got == want {
        Ok(())
    } else {
        Err(format!("{}: expected {:?}, got {:?}", label, want, got))
    }
}

fn want_finite(exp: i64, mant: u64) -> Out {
    Out::Finite { exp, mant }
}

/// The identities of CONTRACT 2.3, 2.4, 14.6 on a drawn operand `x`, one and two as `Mag`s.
pub fn identity_case(ops: &Ops, src: &mut Src) -> Result<u8, Bad> {
    let c = draw_case(src, &Rel::Any).map_err(bare)?;
    let run = || -> Result<u8, String> {
        let x = mag_from(c.ma, c.ea)?;
        let (zero, one, two) = (Mag::zero(), mag_from(TOP, 1)?, mag_from(TOP, 2)?);
        let (ex_, mx) = (c.ea, c.ma);
        let me = want_finite(ex_, mx);
        let (one_o, zero_o) = (want_finite(1, TOP), Out::Zero);
        require_eq("x + 0", &ap2(ops.add_up, &x, &zero)?, &me)?;
        require_eq("0 + x", &ap2(ops.add_up, &zero, &x)?, &me)?;
        require_eq("x * 1 up", &ap2(ops.mul_up, &x, &one)?, &me)?;
        require_eq("1 * x up", &ap2(ops.mul_up, &one, &x)?, &me)?;
        require_eq("x * 1 lower", &ap2(ops.mul_lower, &x, &one)?, &me)?;
        require_eq("1 * x lower", &ap2(ops.mul_lower, &one, &x)?, &me)?;
        require_eq("x / 1 up", &ap2(ops.div_up, &x, &one)?, &me)?;
        require_eq("x / 1 lower", &ap2(ops.div_lower, &x, &one)?, &me)?;
        require_eq("x / x up", &ap2(ops.div_up, &x, &x)?, &one_o)?;
        require_eq("x / x lower", &ap2(ops.div_lower, &x, &x)?, &one_o)?;
        require_eq("x - 0", &ap2(ops.sub_lower, &x, &zero)?, &me)?;
        require_eq("x - x", &ap2(ops.sub_lower, &x, &x)?, &zero_o)?;
        require_eq("0 - x", &ap2(ops.sub_lower, &zero, &x)?, &zero_o)?;
        require_eq("x * 0 up", &ap2(ops.mul_up, &x, &zero)?, &zero_o)?;
        require_eq("x * 0 lower", &ap2(ops.mul_lower, &x, &zero)?, &zero_o)?;
        require_eq("0 / x up", &ap2(ops.div_up, &zero, &x)?, &zero_o)?;
        require_eq("0 / x lower", &ap2(ops.div_lower, &zero, &x)?, &zero_o)?;
        require_eq("x / 0 up", &ap2(ops.div_up, &x, &zero)?, &Out::Infinity)?;
        require_eq("x / 0 lower", &ap2(ops.div_lower, &x, &zero)?, &zero_o)?;
        require_eq("sqrt 0 up", &ap1(ops.sqrt_up, &zero)?, &zero_o)?;
        require_eq("sqrt 0 lower", &ap1(ops.sqrt_lower, &zero)?, &zero_o)?;
        // doubling is exact: x + x = x * 2 = (exp + 1, mant), unless the exponent leaves the range
        let (dbl_up, dbl_lo) = match ex_.checked_add(1) {
            Some(e) => (want_finite(e, mx), want_finite(e, mx)),
            None => (Out::Infinity, want_finite(i64::MAX, u64::MAX)),
        };
        require_eq("x + x", &ap2(ops.add_up, &x, &x)?, &dbl_up)?;
        require_eq("x * 2 up", &ap2(ops.mul_up, &x, &two)?, &dbl_up)?;
        require_eq("x * 2 lower", &ap2(ops.mul_lower, &x, &two)?, &dbl_lo)?;
        Ok(FIN)
    };
    run().map_err(|m| bad("identity", &c, m))
}

/// sqrt of exact squares: `s = t * 2^32 * 2^(es - 64)` with `t` a 32-bit integer with the top bit set has the exact square `t^2 * 2^(2 es - 64)`;
/// and of the powers of two `2^(e - 1)` with odd `e`.
pub fn sqrt_exact_case(ops: &Ops, src: &mut Src) -> Result<u8, Bad> {
    let prec = draw_precision(src.rng, &PrecSource::Fast).map_err(bare)?;
    let dom = ExpDomain {
        min: -(1 << 61),
        max: 1 << 61,
        subnormal: None,
    };
    // a forced exponent e puts the square near 2^e: es = e / 2 (inside `dom` for every forced value, which are within +-(2^62 - 1))
    let es = match src.next_forced() {
        Some(e) => e / 2,
        None => draw_exponent(src.rng, &dom, 0, prec).map_err(bare)?.value,
    };
    let pick = src.below(8);
    let t_edge: [u64; 6] = [
        1 << 31,
        (1 << 31) + 1,
        (1 << 32) - 1,
        3_037_000_499,
        3_037_000_500,
        3_037_000_501,
    ];
    let edge_i = usize::try_from(src.below(6)).map_err(|_| bare(String::from("edge index")))?;
    let free_t = (1u64 << 31) + src.below(1 << 31);
    let t = match (pick < 3, t_edge.get(edge_i)) {
        (true, Some(e)) => *e,
        _ => free_t,
    };
    let c = Case {
        prec: api_prec(prec),
        ea: es,
        ma: t,
        eb: 0,
        mb: 0,
    };
    let run = || -> Result<u8, String> {
        let sq = t * t;
        let (sq_exp, sq_mant) = if sq >= TOP {
            (2 * es, sq)
        } else {
            (2 * es - 1, sq << 1)
        };
        let x = mag_from(sq_mant, sq_exp)?;
        let want = want_finite(es, t << 32);
        require_eq("sqrt_up of a square", &ap1(ops.sqrt_up, &x)?, &want)?;
        require_eq("sqrt_lower of a square", &ap1(ops.sqrt_lower, &x)?, &want)?;
        // just above the square the largest lower bound is still the root
        if let Some((ne, nm)) = next_up(sq_exp, sq_mant) {
            let above = mag_from(nm, ne)?;
            require_eq(
                "sqrt_lower just above a square",
                &ap1(ops.sqrt_lower, &above)?,
                &want,
            )?;
            let up = ap1(ops.sqrt_up, &above)?;
            if up == want {
                return Err(String::from(
                    "sqrt_up just above a square is not above the root",
                ));
            }
        }
        // a power of two with an even binary exponent: x = 2^(e - 1), e odd, root 2^((e - 1) / 2) = (((e - 1) / 2) + 1, 2^63)
        let odd = 2 * (es / 2) + 1;
        let p = mag_from(TOP, odd)?;
        let root = want_finite((odd - 1) / 2 + 1, TOP);
        require_eq("sqrt_up of a power of two", &ap1(ops.sqrt_up, &p)?, &root)?;
        require_eq(
            "sqrt_lower of a power of two",
            &ap1(ops.sqrt_lower, &p)?,
            &root,
        )?;
        Ok(FIN | EXACT)
    };
    run().map_err(|m| bad("sqrt_exact", &c, m)).map_err(unary)
}

// ---- canonical form -----------------------------------------------------------------------------------------------------------

fn operand(src: &mut Src, prec: u32, anchor: i64) -> Result<Mag, String> {
    let class = src.below(8);
    let e = exp_of(src, prec, anchor)?;
    let m = draw_mantissa(src.rng, prec)?.value;
    match class {
        0 => Ok(Mag::zero()),
        1 => Ok(Mag::infinity()),
        _ => mag_from(m, e),
    }
}

/// What a zero or Infinity operand must give: a constant, or the first or second operand unchanged.
#[derive(Debug)]
enum Want {
    Zero,
    Infinity,
    First,
    Second,
}

/// The specified result of a binary operation when an operand is zero or Infinity (CONTRACT 2.3, 14.6); `None` when both are finite and nonzero.
fn special2(name: &str, a: &Out, b: &Out) -> Option<Want> {
    let any_inf = *a == Out::Infinity || *b == Out::Infinity;
    let (a0, b0) = (*a == Out::Zero, *b == Out::Zero);
    if !any_inf && !a0 && !b0 {
        return None;
    }
    Some(match name {
        "add_up" if any_inf => Want::Infinity,
        "add_up" if a0 => Want::Second,
        "add_up" => Want::First,
        "sub_lower" if any_inf || a0 => Want::Zero,
        "sub_lower" => Want::First,
        "mul_up" if any_inf => Want::Infinity,
        "div_up" if any_inf || b0 => Want::Infinity,
        // mul_up with a zero operand, mul_lower, div_up of a zero numerator, div_lower: zero
        _ => Want::Zero,
    })
}

/// The specified result of a square root of zero or Infinity; `None` for a finite nonzero operand.
fn special1(name: &str, a: &Out) -> Option<Want> {
    match (name, a) {
        ("sqrt_up", Out::Infinity) => Some(Want::Infinity),
        (_, Out::Infinity) | (_, Out::Zero) => Some(Want::Zero),
        _ => None,
    }
}

fn meets(got: &Out, want: &Want, a: &Out, b: &Out) -> bool {
    match want {
        Want::Zero => *got == Out::Zero,
        Want::Infinity => *got == Out::Infinity,
        Want::First => got == a,
        Want::Second => got == b,
    }
}

/// Every result of every operation is canonical (the check is inside `out_of`), whatever mix of finite, zero and Infinity operands, and a
/// zero or Infinity operand gives the specified value (0 + 0, x + 0, 0 - x, x - Infinity, 0 * Infinity, 0 / 0, x / 0, sqrt(Infinity), ...).
pub fn canonical_case(ops: &Ops, src: &mut Src) -> Result<u8, Bad> {
    let c = draw_case(src, &Rel::Sum).map_err(bare)?;
    let mut run = || -> Result<u8, String> {
        let a = operand(src, gen_prec(c.prec)?, c.ea)?;
        let b = operand(src, gen_prec(c.prec)?, c.eb)?;
        let (oa, ob) = (out_of(&a)?, out_of(&b)?);
        for (name, f) in [
            ("add_up", ops.add_up),
            ("sub_lower", ops.sub_lower),
            ("mul_up", ops.mul_up),
            ("mul_lower", ops.mul_lower),
            ("div_up", ops.div_up),
            ("div_lower", ops.div_lower),
        ] {
            let got = ap2(f, &a, &b).map_err(|m| format!("{}: {}", name, m))?;
            if let Some(want) = special2(name, &oa, &ob) {
                if !meets(&got, &want, &oa, &ob) {
                    return Err(format!(
                        "{}({:?}, {:?}): expected {:?}, got {:?}",
                        name, oa, ob, want, got
                    ));
                }
            }
        }
        for (name, f) in [("sqrt_up", ops.sqrt_up), ("sqrt_lower", ops.sqrt_lower)] {
            let got = ap1(f, &a).map_err(|m| format!("{}: {}", name, m))?;
            if let Some(want) = special1(name, &oa) {
                if !meets(&got, &want, &oa, &oa) {
                    return Err(format!(
                        "{}({:?}): expected {:?}, got {:?}",
                        name, oa, want, got
                    ));
                }
            }
        }
        let mut flags = FIN;
        for o in [&oa, &ob] {
            match o {
                Out::Infinity => flags |= INF,
                Out::Zero => flags |= ZERO,
                Out::Finite { .. } => {}
            }
        }
        Ok(flags)
    };
    run().map_err(|m| bad("canonical", &c, m))
}

/// A generator precision (`u32`) as the API's `Precision` (`usize`, at least 32 bits on every supported target).
fn api_prec(prec: u32) -> usize {
    usize::try_from(prec).unwrap_or(usize::MAX)
}

/// An API precision back to the generator's `u32` (the generators draw at most `MAX_PREC` bits).
fn gen_prec(prec: usize) -> Result<u32, String> {
    u32::try_from(prec).map_err(|_| format!("prec {} is not a generator precision", prec))
}

/// A `Precision` in `i128` (always fits).
fn prec_i128(prec: usize) -> Result<i128, String> {
    i128::try_from(prec).map_err(|_| format!("prec {} does not fit i128", prec))
}

// ---- ulp helpers --------------------------------------------------------------------------------------------------------------

fn draw_prec_wide(src: &mut Src) -> Result<usize, String> {
    match src.below(8) {
        0..=2 => draw_precision(src.rng, &PrecSource::Schedule).map(api_prec),
        3 => Ok(0),
        4 => Ok(usize::MAX),
        5 => Ok(1 << 31),
        6 => Ok((1 << 31) - 1),
        _ => usize::try_from(src.below(1 << 32)).map_err(|_| String::from("precision draw")),
    }
}

/// `2^k` as a `Mag` is `(k + 1, 2^63)`; above `i64::MAX` it is Infinity, below `i64::MIN` the smallest positive Mag (CONTRACT 3.2, 14.7).
fn expect_pow2(label: &str, got: &Out, k: i128) -> Result<u8, String> {
    let field = k + 1;
    let (want, class) = if field > ex(i64::MAX) {
        (Out::Infinity, INF)
    } else if field < ex(i64::MIN) {
        (want_finite(i64::MIN, TOP), ZERO)
    } else {
        let e = i64::try_from(field).map_err(|_| String::from("field conversion"))?;
        (want_finite(e, TOP), FIN)
    };
    require_eq(label, got, &want)?;
    Ok(class)
}

/// `ulp_mag`, `half_ulp` (exact powers of two with the clamps), `smaller_than_half_ulp` (exact comparison with `2^(mid - prec - 1)`), their
/// relations (`ulp_mag(e, p) = half_ulp(e, p - 1)`) and the monotonicity of `ulp_mag` (up in the exponent, down in the precision).
pub fn ulp_case(ops: &Ops, src: &mut Src) -> Result<u8, Bad> {
    let pw = draw_prec_wide(src).map_err(bare)?;
    let anchor_prec = draw_precision(src.rng, &PrecSource::Fast).map_err(bare)?;
    let mid = exp_of(src, anchor_prec, 0).map_err(bare)?;
    let sel = src.below(8);
    let free_e = exp_of(src, anchor_prec, mid).map_err(bare)?;
    let xm = draw_mantissa(src.rng, anchor_prec).map_err(bare)?.value;
    let delta = i128::from(src.below(5)) - 2;
    let thr = ex(mid) - prec_i128(pw).map_err(bare)? - 1;
    let near = i64::try_from(thr + delta).unwrap_or(free_e);
    let c = Case {
        prec: pw,
        ea: mid,
        ma: xm,
        eb: near,
        mb: xm,
    };
    let run = || -> Result<u8, String> {
        let mut class = 0u8;
        let ulp = out_of(&(ops.ulp_mag)(mid, pw))?;
        class |= expect_pow2("ulp_mag", &ulp, ex(mid) - prec_i128(pw)?)?;
        let half = out_of(&(ops.half_ulp)(mid, pw))?;
        class |= expect_pow2("half_ulp", &half, thr)?;
        if let Some(p1) = pw.checked_sub(1) {
            require_eq(
                "ulp_mag(e, p) = half_ulp(e, p - 1)",
                &ulp,
                &out_of(&(ops.half_ulp)(mid, p1))?,
            )?;
        }
        if let Some(m1) = mid.checked_add(1) {
            le(
                "ulp_mag in the exponent",
                &ulp,
                &out_of(&(ops.ulp_mag)(m1, pw))?,
            )?;
        }
        if let Some(p1) = pw.checked_add(1) {
            le(
                "ulp_mag in the precision",
                &out_of(&(ops.ulp_mag)(mid, p1))?,
                &ulp,
            )?;
        }
        let (x, want) = match sel {
            0 => (Mag::zero(), true),
            1 => (Mag::infinity(), false),
            2..=4 => (
                mag_from(xm, near)?,
                cmp_scaled(wide(xm), ex(near) - 64, 1, thr)? == Ordering::Less,
            ),
            _ => (
                mag_from(xm, free_e)?,
                cmp_scaled(wide(xm), ex(free_e) - 64, 1, thr)? == Ordering::Less,
            ),
        };
        let got = (ops.smaller_than_half_ulp)(&x, mid, pw);
        if got != want {
            return Err(format!(
                "smaller_than_half_ulp({:?}, mid {}, prec {}) = {}, exact comparison says {}",
                out_of(&x)?,
                mid,
                pw,
                got,
                want
            ));
        }
        Ok(class | if want { EXACT } else { FIN })
    };
    run().map_err(|m| bad("ulp", &c, m))
}
