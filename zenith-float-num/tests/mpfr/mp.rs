//! Safe stand-in for the parts of MPFR (and rug) the comparison tests use (Paul call 2026-11-24, r6 item 6).
//!
//! Values are held exactly (sign, `prec`-bit mantissa, MPFR exponent: `0.1xxx * 2^exp`). Every
//! value-producing function (arithmetic, roots, transcendentals, constants, string conversion)
//! takes its correctly rounded result from the test's golden stream ([`super::gold`]): python-flint,
//! unbounded exponent range, plus the ternary value. The exponent-range handling MPFR applies to
//! such a result (`mpfr_check_range`, `mpfr_underflow`, `mpfr_overflow`), and the pure
//! representation functions (`set`, `prec_round`, `mul_2si`, `subnormalize`, `nextabove`,
//! `nextbelow`, `get_exp`, `min_prec`) are ported line by line from MPFR 4.2.2
//! (`exceptions.c`, `set.c`, `round_prec.c`, `mul_2si.c`, `subnormal.c`, `next.c`).

use super::gold;
use std::cell::Cell;
use std::cmp::Ordering;
use std::fmt;

/// A harness step that must succeed (a conversion or a golden parse): its failure fails the test with the
/// message (see [`gold::must`]).
trait OrFail<T, E> {
    fn or_fail(self, msg: impl FnOnce(E) -> String) -> T;
}
impl<T: Default, E> OrFail<T, E> for Result<T, E> {
    #[track_caller]
    fn or_fail(self, msg: impl FnOnce(E) -> String) -> T {
        gold::must(self.map_err(msg))
    }
}
impl<T: Default> OrFail<T, ()> for Option<T> {
    #[track_caller]
    fn or_fail(self, msg: impl FnOnce(()) -> String) -> T {
        gold::must(self.ok_or(()).map_err(msg))
    }
}

// ---------------------------------------------------------------- big unsigned integers

/// Unsigned integer, little-endian 64-bit limbs, no high zero limbs.
#[derive(PartialEq, Eq, Debug, Default)]
pub struct Big(Vec<u64>);

impl Big {
    pub fn zero() -> Big {
        Big(Vec::new())
    }
    pub fn from_u64(v: u64) -> Big {
        let mut b = Big(vec![v]);
        b.trim();
        b
    }
    pub fn from_u128(v: u128) -> Big {
        let mut b = Big(vec![v as u64, (v >> 64) as u64]);
        b.trim();
        b
    }
    pub fn from_vec(l: Vec<u64>) -> Big {
        let mut b = Big(l);
        b.trim();
        b
    }
    pub fn limbs(&self) -> &[u64] {
        &self.0
    }
    fn trim(&mut self) {
        while self.0.last() == Some(&0) {
            self.0.pop();
        }
    }
    pub fn from_hex(s: &str) -> Big {
        let b = s.as_bytes();
        let mut limbs = Vec::with_capacity(b.len() / 16 + 1);
        let mut end = b.len();
        while end > 0 {
            let start = end.saturating_sub(16);
            let chunk = std::str::from_utf8(&b[start..end]).or_fail(|e| format!("hex: {e}"));
            limbs.push(u64::from_str_radix(chunk, 16).or_fail(|e| format!("hex {chunk}: {e}")));
            end = start;
        }
        Big::from_vec(limbs)
    }
}

impl std::ops::Deref for Big {
    type Target = [u64];
    fn deref(&self) -> &[u64] {
        &self.0
    }
}

/// Compares two limb slices by value (high zero limbs ignored).
pub fn cmp_nat(a: &[u64], b: &[u64]) -> Ordering {
    let (la, lb) = (a.bits(), b.bits());
    if la != lb {
        return la.cmp(&lb);
    }
    for i in (0..a.len().min(b.len())).rev() {
        match a[i].cmp(&b[i]) {
            Ordering::Equal => {}
            x => return x,
        }
    }
    Ordering::Equal
}

/// Unsigned-integer arithmetic on little-endian 64-bit limbs, read in place (high zero limbs allowed): a
/// [`Big`] derefs to its limbs, and a borrowed word slice (an `ExactNum`'s raw words) is read without a copy.
pub trait Nat {
    fn is_zero(&self) -> bool;
    fn bits(&self) -> u64;
    fn bit(&self, i: u64) -> bool;
    fn trailing_zeros(&self) -> u64;
    fn low_nonzero(&self, n: u64) -> bool;
    fn shl(&self, n: u64) -> Big;
    fn shr(&self, n: u64) -> Big;
    fn add(&self, o: &[u64]) -> Big;
    fn sub(&self, o: &[u64]) -> Big;
    fn mul(&self, o: &[u64]) -> Big;
    fn divrem(&self, d: &[u64]) -> (Big, Big);
    fn to_hex(&self) -> String;
}

impl Nat for [u64] {
    fn is_zero(&self) -> bool {
        self.iter().all(|&x| x == 0)
    }
    fn bits(&self) -> u64 {
        match self.iter().rposition(|&x| x != 0) {
            None => 0,
            Some(i) => 64 * i as u64 + (64 - u64::from(self[i].leading_zeros())),
        }
    }
    fn bit(&self, i: u64) -> bool {
        let (w, b) = ((i / 64) as usize, i % 64);
        self.get(w).is_some_and(|x| (x >> b) & 1 == 1)
    }
    fn trailing_zeros(&self) -> u64 {
        let mut n = 0u64;
        for &x in self {
            if x == 0 {
                n += 64;
            } else {
                return n + u64::from(x.trailing_zeros());
            }
        }
        n
    }
    /// True when any of the low `n` bits is set.
    fn low_nonzero(&self, n: u64) -> bool {
        let (w, b) = ((n / 64) as usize, n % 64);
        self.iter().take(w).any(|&x| x != 0)
            || (b > 0 && self.get(w).is_some_and(|x| x & ((1u64 << b) - 1) != 0))
    }
    fn shl(&self, n: u64) -> Big {
        if self.is_zero() {
            return Big::zero();
        }
        let (w, b) = ((n / 64) as usize, (n % 64) as u32);
        let mut r = vec![0u64; w];
        if b == 0 {
            r.extend_from_slice(self);
        } else {
            let mut carry = 0u64;
            for &x in self {
                r.push((x << b) | carry);
                carry = x >> (64 - b);
            }
            r.push(carry);
        }
        let mut out = Big(r);
        out.trim();
        out
    }
    /// `floor(self / 2^n)`.
    fn shr(&self, n: u64) -> Big {
        let (w, b) = ((n / 64) as usize, (n % 64) as u32);
        if w >= self.len() {
            return Big::zero();
        }
        let src = &self[w..];
        let mut r = Vec::with_capacity(src.len());
        for i in 0..src.len() {
            let lo = src[i] >> b;
            let hi = if b > 0 { src.get(i + 1).map_or(0, |x| x << (64 - b)) } else { 0 };
            r.push(lo | hi);
        }
        let mut out = Big(r);
        out.trim();
        out
    }
    fn add(&self, o: &[u64]) -> Big {
        let n = self.len().max(o.len());
        let mut r = Vec::with_capacity(n + 1);
        let mut c = 0u128;
        for i in 0..n {
            let s =
                u128::from(*self.get(i).unwrap_or(&0)) + u128::from(*o.get(i).unwrap_or(&0)) + c;
            r.push(s as u64);
            c = s >> 64;
        }
        r.push(c as u64);
        let mut out = Big(r);
        out.trim();
        out
    }
    /// `self - o`; `o <= self` is required.
    fn sub(&self, o: &[u64]) -> Big {
        assert!(cmp_nat(o, self) != Ordering::Greater, "Big::sub underflow");
        let mut r = Vec::with_capacity(self.len());
        let mut borrow = 0i128;
        for (i, &a) in self.iter().enumerate() {
            let mut d = i128::from(a) - i128::from(*o.get(i).unwrap_or(&0)) - borrow;
            borrow = 0;
            if d < 0 {
                d += 1i128 << 64;
                borrow = 1;
            }
            r.push(d as u64);
        }
        let mut out = Big(r);
        out.trim();
        out
    }
    fn mul(&self, o: &[u64]) -> Big {
        if self.is_zero() || o.is_zero() {
            return Big::zero();
        }
        let mut r = vec![0u64; self.len() + o.len()];
        for (i, &a) in self.iter().enumerate() {
            let mut c = 0u128;
            for (j, &b) in o.iter().enumerate() {
                let t = u128::from(a) * u128::from(b) + u128::from(r[i + j]) + c;
                r[i + j] = t as u64;
                c = t >> 64;
            }
            let mut k = i + o.len();
            while c != 0 {
                let t = u128::from(r[k]) + c;
                r[k] = t as u64;
                c = t >> 64;
                k += 1;
            }
        }
        let mut out = Big(r);
        out.trim();
        out
    }
    /// `(self / d, self % d)`, word by word (Knuth, TAOCP 4.3.1 algorithm D); `d != 0`.
    fn divrem(&self, d: &[u64]) -> (Big, Big) {
        let n = self.iter().rposition(|&x| x != 0).map_or(0, |i| i + 1);
        let m = d.iter().rposition(|&x| x != 0).map_or(0, |i| i + 1);
        assert!(m > 0, "Big::divrem by zero");
        let (u0, v0) = (&self[..n], &d[..m]);
        if cmp_nat(u0, v0) == Ordering::Less {
            return (Big::zero(), Big::from_vec(u0.to_owned()));
        }
        if m == 1 {
            let dv = u128::from(v0[0]);
            let mut q = vec![0u64; n];
            let mut r = 0u128;
            for i in (0..n).rev() {
                let cur = (r << 64) | u128::from(u0[i]);
                q[i] = (cur / dv) as u64;
                r = cur % dv;
            }
            return (Big::from_vec(q), Big::from_u128(r));
        }
        // normalize: the divisor's top word gets its top bit set
        let s = v0[m - 1].leading_zeros();
        let v = v0.shl(u64::from(s)).0;
        let mut u = u0.shl(u64::from(s)).0;
        u.resize(n + 1, 0);
        let (vt, vs) = (u128::from(v[m - 1]), u128::from(v[m - 2]));
        let b = 1u128 << 64;
        let mut q = vec![0u64; n - m + 1];
        for j in (0..=n - m).rev() {
            let num = (u128::from(u[j + m]) << 64) | u128::from(u[j + m - 1]);
            let mut qhat = num / vt;
            let mut rhat = num % vt;
            while qhat >= b || qhat * vs > ((rhat << 64) | u128::from(u[j + m - 2])) {
                qhat -= 1;
                rhat += vt;
                if rhat >= b {
                    break;
                }
            }
            // u[j..=j+m] -= qhat * v
            let mut k: i128 = 0;
            for i in 0..m {
                let p = qhat * u128::from(v[i]);
                let t = i128::from(u[i + j]) - k - i128::from(p as u64);
                u[i + j] = t as u64;
                k = i128::from((p >> 64) as u64) - (t >> 64);
            }
            let t = i128::from(u[j + m]) - k;
            u[j + m] = t as u64;
            let mut qj = qhat as u64;
            if t < 0 {
                // qhat was one too large: add v back
                qj -= 1;
                let mut c = 0u128;
                for i in 0..m {
                    let sum = u128::from(u[i + j]) + u128::from(v[i]) + c;
                    u[i + j] = sum as u64;
                    c = sum >> 64;
                }
                u[j + m] = u[j + m].wrapping_add(c as u64);
            }
            q[j] = qj;
        }
        u.truncate(m);
        (Big::from_vec(q), Big::from_vec(u).shr(u64::from(s)))
    }
    fn to_hex(&self) -> String {
        match self[..self.iter().rposition(|&x| x != 0).map_or(0, |i| i + 1)].split_last() {
            None => "0".to_owned(),
            Some((top, rest)) => {
                let mut s = format!("{top:x}");
                for x in rest.iter().rev() {
                    s.push_str(&format!("{x:016x}"));
                }
                s
            }
        }
    }
}

impl PartialOrd for Big {
    fn partial_cmp(&self, o: &Big) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Big {
    fn cmp(&self, o: &Big) -> Ordering {
        self.0.len().cmp(&o.0.len()).then_with(|| {
            for (a, b) in self.0.iter().rev().zip(o.0.iter().rev()) {
                match a.cmp(b) {
                    Ordering::Equal => {}
                    x => return x,
                }
            }
            Ordering::Equal
        })
    }
}

// ---------------------------------------------------------------- rounding modes and context

/// MPFR rounding mode (`MPFR_RNDN`, `RNDZ`, `RNDU`, `RNDD`, `RNDA`); also rug's `Round`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Rnd {
    Nearest,
    Zero,
    Up,
    Down,
    AwayZero,
}
pub use Rnd as Round;

impl Rnd {
    fn tag(self) -> &'static str {
        match self {
            Rnd::Nearest => "N",
            Rnd::Zero => "Z",
            Rnd::Up => "U",
            Rnd::Down => "D",
            Rnd::AwayZero => "A",
        }
    }
    /// `MPFR_IS_LIKE_RNDZ(rnd, neg)`.
    fn like_rndz(self, neg: bool) -> bool {
        self == Rnd::Zero || (self == Rnd::Up && neg) || (self == Rnd::Down && !neg)
    }
}

pub const EMIN_MIN: i64 = 1 - (1i64 << 62);
pub const EMAX_MAX: i64 = (1i64 << 62) - 1;
const EMAX_DEFAULT: i64 = (1i64 << 30) - 1;

const FLAG_UNDERFLOW: u32 = 1;
const FLAG_OVERFLOW: u32 = 2;
const FLAG_NAN: u32 = 4;
const FLAG_INEXACT: u32 = 8;

thread_local! {
    static EMIN: Cell<i64> = const { Cell::new(-EMAX_DEFAULT) };
    static EMAX: Cell<i64> = const { Cell::new(EMAX_DEFAULT) };
    static FLAGS: Cell<u32> = const { Cell::new(0) };
}

fn emin() -> i64 {
    EMIN.with(Cell::get)
}
fn emax() -> i64 {
    EMAX.with(Cell::get)
}
fn raise(f: u32) {
    FLAGS.with(|c| c.set(c.get() | f));
}

// ---------------------------------------------------------------- floats

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Kind {
    Nan,
    Inf,
    Zero,
    Reg,
}

/// An MPFR-format float: for a regular value, `m` has exactly `prec` bits and the value is
/// `(-1)^neg * m * 2^(exp - prec)`, i.e. `0.1xxx * 2^exp` (MPFR's exponent convention).
#[derive(Debug)]
pub struct Float {
    prec: u32,
    kind: Kind,
    neg: bool,
    m: Big,
    /// MPFR's exponent (`0.1xxx * 2^exp`); 0 for NaN, Inf and zero.
    pub exp: i64,
}

/// rug's `Special` values.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Special {
    Zero,
    NegZero,
    Infinity,
    NegInfinity,
    Nan,
}

/// Result of rounding `m * 2^e` (exact, plus `sticky` for nonzero bits below) to `prec` bits.
struct Rounded {
    m: Big,
    /// value is `m * 2^e2` with `m` of exactly `prec` bits.
    e2: i64,
    /// sign of (rounded - exact) in magnitude: -1, 0, 1.
    dir: i32,
    /// round-to-nearest tie (MPFR_EVEN_INEX).
    tie: bool,
}

fn round_mag(m: &[u64], e: i64, sticky: bool, prec: u32, rnd: Rnd, neg: bool) -> Rounded {
    let p = u64::from(prec);
    let l = m.bits();
    assert!(l > 0, "round_mag of zero");
    if l <= p {
        assert!(!sticky, "round_mag: sticky bits below a short mantissa");
        return Rounded {
            m: m.shl(p - l),
            e2: e - (p - l) as i64,
            dir: 0,
            tie: false,
        };
    }
    let k = l - p;
    let mut q = m.shr(k);
    let rb = m.bit(k - 1);
    let sb = sticky || m.low_nonzero(k - 1);
    let mut e2 = e + k as i64;
    if !rb && !sb {
        return Rounded {
            m: q,
            e2,
            dir: 0,
            tie: false,
        };
    }
    let up = match rnd {
        Rnd::Nearest => rb && (sb || q.bit(0)),
        Rnd::Zero => false,
        Rnd::Up => !neg,
        Rnd::Down => neg,
        Rnd::AwayZero => true,
    };
    let tie = rnd == Rnd::Nearest && rb && !sb;
    if up {
        q = q.add(&Big::from_u64(1));
        if q.bits() > p {
            q = q.shr(1);
            e2 += 1;
        }
    }
    Rounded {
        m: q,
        e2,
        dir: if up { 1 } else { -1 },
        tie,
    }
}

fn signed(neg: bool, dir: i32) -> i32 {
    if neg {
        -dir
    } else {
        dir
    }
}

impl Float {
    fn singular(prec: u32, kind: Kind, neg: bool) -> Float {
        Float {
            prec,
            kind,
            neg,
            m: Big::zero(),
            exp: 0,
        }
    }
    pub fn new(prec: u32) -> Float {
        Float::singular(prec, Kind::Zero, false)
    }
    pub fn with_val<T: Assign>(prec: u32, v: T) -> Float {
        let mut f = Float::new(prec);
        let _ = v.assign_to(&mut f, Rnd::Nearest);
        f
    }
    pub fn with_val_round<T: Assign>(prec: u32, v: T, r: Rnd) -> (Float, Ordering) {
        let mut f = Float::new(prec);
        let t = v.assign_to(&mut f, r);
        (f, t.cmp(&0))
    }
    pub fn assign_round<T: Assign>(&mut self, v: T, r: Rnd) -> Ordering {
        v.assign_to(self, r).cmp(&0)
    }
    /// Exact value `(-1)^neg * m * 2^e` rounded to `self.prec` (no range check: `mpfr_set`-like).
    fn set_exact(&mut self, neg: bool, m: &[u64], e: i64, sticky: bool, rnd: Rnd) -> i32 {
        if m.is_zero() {
            assert!(!sticky, "set_exact: sticky without bits");
            *self = Float::singular(self.prec, Kind::Zero, neg);
            return 0;
        }
        let r = round_mag(m, e, sticky, self.prec, rnd, neg);
        self.kind = Kind::Reg;
        self.neg = neg;
        self.exp = r.e2 + i64::from(self.prec);
        self.m = r.m;
        signed(neg, r.dir)
    }
    pub fn as_raw(&self) -> &Float {
        self
    }
    pub fn as_raw_mut(&mut self) -> &mut Float {
        self
    }
    pub fn prec(&self) -> u32 {
        self.prec
    }
    pub fn is_nan(&self) -> bool {
        self.kind == Kind::Nan
    }
    pub fn is_infinite(&self) -> bool {
        self.kind == Kind::Inf
    }
    pub fn is_finite(&self) -> bool {
        matches!(self.kind, Kind::Zero | Kind::Reg)
    }
    pub fn is_zero(&self) -> bool {
        self.kind == Kind::Zero
    }
    pub fn is_normal(&self) -> bool {
        self.kind == Kind::Reg
    }
    pub fn is_sign_negative(&self) -> bool {
        self.neg
    }
    pub fn is_sign_positive(&self) -> bool {
        !self.neg
    }
    pub fn is_integer(&self) -> bool {
        match self.kind {
            Kind::Zero => true,
            Kind::Reg => {
                let tz = self.m.trailing_zeros() as i64;
                self.exp - i64::from(self.prec) + tz >= 0
            }
            _ => false,
        }
    }
    /// rug's `get_exp`: the exponent of a regular value (rug's type is `i32`; panics past it).
    pub fn get_exp(&self) -> Option<i32> {
        self.is_normal()
            .then(|| i32::try_from(self.exp).or_fail(|e| format!("get_exp {}: {e}", self.exp)))
    }
    /// rug's `next_up` (`mpfr_nextabove`).
    pub fn next_up(&mut self) {
        mpfr::nextabove(self);
    }
    /// rug's `next_down` (`mpfr_nextbelow`).
    pub fn next_down(&mut self) {
        mpfr::nextbelow(self);
    }
    /// rug's `Float::i_exp(i, e)`: exactly `i * 2^e`.
    pub fn i_exp(i: i32, e: i32) -> Pending<'static> {
        let v = i64::from(i);
        let m = Big::from_u64(v.unsigned_abs());
        if m.is_zero() {
            return Pending::Exact(Val::O(Float::new(2)));
        }
        Pending::Exact(Val::O(Float::from_parts(64, v < 0, &m, i64::from(e))))
    }
    pub fn acos_ref(&self) -> Pending<'_> {
        self.op1("acos")
    }
    pub fn asin_ref(&self) -> Pending<'_> {
        self.op1("asin")
    }
    pub fn acosh_ref(&self) -> Pending<'_> {
        self.op1("acosh")
    }
    pub fn asinh_ref(&self) -> Pending<'_> {
        self.op1("asinh")
    }
    pub fn atanh_ref(&self) -> Pending<'_> {
        self.op1("atanh")
    }
    /// rug's `recip_ref` (`mpfr_ui_div(1, x)`).
    pub fn recip_ref(&self) -> Pending<'_> {
        Pending::op("ui_div", &[Arg::U(1), Arg::F(self)])
    }
    pub fn hypot_ref(&self, y: &Float) -> Pending<'static> {
        Pending::op("hypot", &[Arg::F(self), Arg::F(y)])
    }
    /// rug's `floor_ref` (`mpfr_rint_floor`): the floor, exact, rounded on assignment.
    pub fn floor_ref(&self) -> Pending<'_> {
        match self.kind {
            Kind::Reg => {
                let e = self.exp - i64::from(self.prec);
                if e >= 0 {
                    return Pending::Exact(Val::B(self));
                }
                let k = e.unsigned_abs();
                let mut q = self.m.shr(k);
                if self.neg && self.m.low_nonzero(k) {
                    q = q.add(&Big::from_u64(1));
                }
                if q.is_zero() {
                    return Pending::Exact(Val::O(Float::singular(2, Kind::Zero, self.neg)));
                }
                let bits = u32::try_from(q.bits()).unwrap_or(u32::MAX).max(2);
                Pending::Exact(Val::O(Float::from_parts(bits, self.neg, &q, 0)))
            }
            _ => Pending::Exact(Val::B(self)),
        }
    }
    pub fn abs(mut self) -> Float {
        self.neg = false;
        self
    }
    pub fn abs_ref(&self) -> Pending<'_> {
        // |x|: x itself, negated (exactly) when x is negative
        if self.neg {
            -self
        } else {
            Pending::Exact(Val::B(self))
        }
    }
    pub fn ln_ref(&self) -> Pending<'_> {
        Pending::op("log", &[Arg::F(self)])
    }
    fn op1(&self, op: &'static str) -> Pending<'_> {
        Pending::op(op, &[Arg::F(self)])
    }
    pub fn exp_ref(&self) -> Pending<'_> {
        self.op1("exp")
    }
    pub fn exp2_ref(&self) -> Pending<'_> {
        self.op1("exp2")
    }
    pub fn exp_m1_ref(&self) -> Pending<'_> {
        self.op1("expm1")
    }
    pub fn ln_1p_ref(&self) -> Pending<'_> {
        self.op1("log1p")
    }
    pub fn log2_ref(&self) -> Pending<'_> {
        self.op1("log2")
    }
    pub fn log10_ref(&self) -> Pending<'_> {
        self.op1("log10")
    }
    pub fn sin_ref(&self) -> Pending<'_> {
        self.op1("sin")
    }
    pub fn cos_ref(&self) -> Pending<'_> {
        self.op1("cos")
    }
    pub fn tan_ref(&self) -> Pending<'_> {
        self.op1("tan")
    }
    pub fn atan_ref(&self) -> Pending<'_> {
        self.op1("atan")
    }
    pub fn sinh_ref(&self) -> Pending<'_> {
        self.op1("sinh")
    }
    pub fn cosh_ref(&self) -> Pending<'_> {
        self.op1("cosh")
    }
    pub fn tanh_ref(&self) -> Pending<'_> {
        self.op1("tanh")
    }
    pub fn square_ref(&self) -> Pending<'_> {
        self.op1("sqr")
    }
    pub fn sqrt_ref(&self) -> Pending<'_> {
        self.op1("sqrt")
    }
    /// rug's `root_ref(k)`: `mpfr_rootn_ui`.
    pub fn root_ref(&self, k: u32) -> Pending<'_> {
        Pending::op("rootn_ui", &[Arg::F(self), Arg::U(u64::from(k))])
    }
    /// rug's `Pow`: `mpfr_pow` (float), `mpfr_pow_ui` (u32), `mpfr_pow_si` (i32).
    pub fn pow<Y: PowArg>(&self, y: Y) -> Pending<'static> {
        y.pow_of(self)
    }
    /// rug's `to_integer`: rounded to the nearest integer, ties to even (`None` for NaN and Inf).
    pub fn to_integer(&self) -> Option<Integer> {
        match self.kind {
            Kind::Zero => Some(Integer::new()),
            Kind::Reg => {
                let e = self.exp - i64::from(self.prec);
                if e >= 0 {
                    return Some(Integer::make(self.neg, self.m.shl(e as u64)));
                }
                let s = e.unsigned_abs();
                let mut q = self.m.shr(s);
                // round half to even: the first dropped bit, then the rest and the parity of q
                if self.m.bit(s - 1) && (self.m.low_nonzero(s - 1) || q.bit(0)) {
                    q = q.add(&Big::from_u64(1));
                }
                Some(Integer::make(self.neg, q))
            }
            _ => None,
        }
    }
    /// Truncated toward zero (a regular value with `exp <= 32`).
    fn trunc_integer(&self) -> Integer {
        let e = self.exp - i64::from(self.prec);
        let m = if e >= 0 { self.m.shl(e as u64) } else { self.m.shr(e.unsigned_abs()) };
        Integer::make(self.neg, m)
    }
    /// Truncated toward zero, saturated to `i32` (`None` for NaN).
    pub fn to_i32_saturating(&self) -> Option<i32> {
        match self.kind {
            Kind::Nan => None,
            Kind::Inf => Some(if self.neg { i32::MIN } else { i32::MAX }),
            Kind::Zero => Some(0),
            _ if self.exp > 32 => Some(if self.neg { i32::MIN } else { i32::MAX }),
            _ => {
                let i = self.trunc_integer();
                if i.mag.bits() > 31 {
                    Some(if i.neg { i32::MIN } else { i32::MAX })
                } else {
                    let v = i32::try_from(i.mag.limbs().first().copied().unwrap_or(0))
                        .unwrap_or(i32::MAX);
                    Some(if i.neg { -v } else { v })
                }
            }
        }
    }
    /// Signed integer mantissa and its exponent: `self = m * 2^e` (regular values).
    pub fn mantissa_exp(&self) -> Option<(bool, &Big, i64)> {
        self.is_normal()
            .then(|| (self.neg, &self.m, self.exp - i64::from(self.prec)))
    }
    /// Builds `(-1)^neg * m * 2^e` exactly at `prec` bits (panics when it does not fit).
    pub fn from_parts(prec: u32, neg: bool, m: &[u64], e: i64) -> Float {
        let mut f = Float::new(prec);
        let t = f.set_exact(neg, m, e, false, Rnd::Zero);
        assert_eq!(t, 0, "from_parts: value needs more than {prec} bits");
        f
    }
    fn unary(self, op: &'static str) -> Float {
        let p = self.prec;
        let mut r = Float::new(p);
        let _ = golden_op(&mut r, op, Rnd::Nearest, &[Arg::F(&self)]);
        r
    }
    pub fn square(self) -> Float {
        self.unary("sqr")
    }
    pub fn ln_1p(self) -> Float {
        self.unary("log1p")
    }
    pub fn sqrt(self) -> Float {
        self.unary("sqrt")
    }
    pub fn exp(self) -> Float {
        self.unary("exp")
    }
    pub fn ln(self) -> Float {
        self.unary("log")
    }
    pub fn sin(self) -> Float {
        self.unary("sin")
    }
    pub fn cos(self) -> Float {
        self.unary("cos")
    }
    /// Nearest `f64` (round to nearest even); out-of-range values saturate to +-inf / +-0.
    pub fn to_f64(&self) -> f64 {
        match self.kind {
            Kind::Nan => f64::NAN,
            Kind::Inf => {
                if self.neg {
                    f64::NEG_INFINITY
                } else {
                    f64::INFINITY
                }
            }
            Kind::Zero => {
                if self.neg {
                    -0.0
                } else {
                    0.0
                }
            }
            Kind::Reg => {
                let r = round_mag(
                    &self.m,
                    self.exp - i64::from(self.prec),
                    false,
                    53,
                    Rnd::Nearest,
                    self.neg,
                );
                let mant = r.m.limbs()[0] as f64;
                let v = mant * 2f64.powi(r.e2.clamp(-1100, 1100) as i32);
                if self.neg {
                    -v
                } else {
                    v
                }
            }
        }
    }
    /// MPFR-style binary digits `[-]1.xxxxe<decimal exponent>`, all `prec` digits; `0`/`-0`, `@NaN@`, `[-]@Inf@`.
    pub fn to_string_radix(&self, radix: i32, digits: Option<usize>) -> String {
        assert!(
            radix == 2 && digits.is_none(),
            "to_string_radix: only radix 2 with all digits"
        );
        match self.kind {
            Kind::Nan => "@NaN@".to_owned(),
            Kind::Inf => {
                if self.neg {
                    "-@Inf@".to_owned()
                } else {
                    "@Inf@".to_owned()
                }
            }
            Kind::Zero => {
                if self.neg {
                    "-0".to_owned()
                } else {
                    "0".to_owned()
                }
            }
            Kind::Reg => {
                let mut s = String::with_capacity(self.prec as usize + 24);
                if self.neg {
                    s.push('-');
                }
                let p = u64::from(self.prec);
                for i in (0..p).rev() {
                    s.push(if self.m.bit(i) { '1' } else { '0' });
                    if i == p - 1 && p > 1 {
                        s.push('.');
                    }
                }
                s.push_str(&format!("e{}", self.exp - 1));
                s
            }
        }
    }
    /// Parses MPFR-style binary text (`[-]d[.ddd][e<dec exp>]`, `@NaN@`, `[-]@Inf@`) exactly.
    pub fn parse_radix<S: AsRef<str>>(s: S, radix: i32) -> Result<Pending<'static>, String> {
        let s = s.as_ref();
        if radix != 2 {
            return Err(format!("parse_radix: radix {radix} unsupported"));
        }
        let (neg, body) = match s.strip_prefix('-') {
            Some(b) => (true, b),
            None => (false, s.strip_prefix('+').unwrap_or(s)),
        };
        if body == "@NaN@" {
            return Ok(Pending::Exact(Val::O(Float::singular(2, Kind::Nan, false))));
        }
        if body == "@Inf@" {
            return Ok(Pending::Exact(Val::O(Float::singular(2, Kind::Inf, neg))));
        }
        let (digits, e) = match body.split_once(['e', '@']) {
            Some((d, e)) => (
                d,
                e.parse::<i64>().map_err(|x| format!("exponent {e}: {x}"))?,
            ),
            None => (body, 0),
        };
        let (ip, fp) = digits.split_once('.').unwrap_or((digits, ""));
        let mut m = Big::zero();
        for c in ip.chars().chain(fp.chars()) {
            m = m.shl(1);
            match c {
                '0' => {}
                '1' => m = m.add(&Big::from_u64(1)),
                _ => return Err(format!("parse_radix: bad digit {c:?}")),
            }
        }
        let e2 = e - fp.len() as i64;
        if m.is_zero() {
            return Ok(Pending::Exact(Val::O(Float::singular(2, Kind::Zero, neg))));
        }
        let bits = u32::try_from(m.bits()).map_err(|x| x.to_string())?.max(2);
        Ok(Pending::Exact(Val::O(Float::from_parts(bits, neg, &m, e2))))
    }
    /// Decimal text, rounded by the golden stream (MPFR `mpfr_strtofr` base 10).
    pub fn parse<S: AsRef<str>>(s: S) -> Result<Pending<'static>, String> {
        Ok(Pending::op("parse10", &[Arg::S(s.as_ref())]))
    }
    /// Canonical golden text: `nan`, `[+-]inf`, `[+-]0`, `[+-]<odd hex>p<exp>`.
    pub fn canon(&self) -> String {
        let s = if self.neg { '-' } else { '+' };
        match self.kind {
            Kind::Nan => "nan".to_owned(),
            Kind::Inf => format!("{s}inf"),
            Kind::Zero => format!("{s}0"),
            Kind::Reg => {
                let tz = self.m.trailing_zeros();
                format!(
                    "{s}{}p{}",
                    self.m.shr(tz).to_hex(),
                    self.exp - i64::from(self.prec) + tz as i64
                )
            }
        }
    }
}

impl fmt::Display for Float {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}[{}]", self.canon(), self.prec)
    }
}
impl fmt::LowerHex for Float {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}[{}]", self.canon(), self.prec)
    }
}

/// Exact comparison (`None` with a NaN), as MPFR's `mpfr_cmp`.
fn cmp_val(a: &Float, b: &Float) -> Option<Ordering> {
    if a.is_nan() || b.is_nan() {
        return None;
    }
    let key = |f: &Float| -> i32 {
        match f.kind {
            Kind::Inf => {
                if f.neg {
                    -2
                } else {
                    2
                }
            }
            Kind::Zero => 0,
            _ => {
                if f.neg {
                    -1
                } else {
                    1
                }
            }
        }
    };
    let (ka, kb) = (key(a), key(b));
    if ka != kb || ka.abs() != 1 {
        return Some(ka.cmp(&kb));
    }
    // both regular, same sign
    let mag = a.exp.cmp(&b.exp).then_with(|| {
        let (pa, pb) = (u64::from(a.prec), u64::from(b.prec));
        let w = pa.max(pb);
        a.m.shl(w - pa).cmp(&b.m.shl(w - pb))
    });
    Some(if a.neg { mag.reverse() } else { mag })
}

impl PartialEq for Float {
    fn eq(&self, o: &Float) -> bool {
        cmp_val(self, o) == Some(Ordering::Equal)
    }
}
impl PartialOrd for Float {
    fn partial_cmp(&self, o: &Float) -> Option<Ordering> {
        cmp_val(self, o)
    }
}

impl std::ops::Neg for Float {
    type Output = Float;
    fn neg(mut self) -> Float {
        self.neg = !self.neg;
        self
    }
}
impl<'a> std::ops::Neg for &'a Float {
    type Output = Pending<'a>;
    fn neg(self) -> Pending<'a> {
        Pending::Neg(Box::new(Pending::Exact(Val::B(self))))
    }
}

/// Exponent kinds for [`Float::pow`].
pub trait PowArg {
    fn pow_of(self, x: &Float) -> Pending<'static>;
}
impl PowArg for &Float {
    fn pow_of(self, x: &Float) -> Pending<'static> {
        Pending::op("pow", &[Arg::F(x), Arg::F(self)])
    }
}
impl PowArg for u32 {
    fn pow_of(self, x: &Float) -> Pending<'static> {
        Pending::op("pow_ui", &[Arg::F(x), Arg::U(u64::from(self))])
    }
}
impl PowArg for i32 {
    fn pow_of(self, x: &Float) -> Pending<'static> {
        Pending::op("pow_si", &[Arg::F(x), Arg::I(i64::from(self))])
    }
}

// ---------------------------------------------------------------- assignment sources

/// An operand of a golden query (borrowed: only its query text is used).
#[derive(Debug)]
pub enum Arg<'a> {
    F(&'a Float),
    I(i64),
    U(u64),
    S(&'a str),
}

/// A float operand that is borrowed from the caller or owned by the pending value (a fresh result).
#[derive(Debug)]
pub enum Val<'a> {
    B(&'a Float),
    O(Float),
}

impl Val<'_> {
    fn get(&self) -> &Float {
        match self {
            Val::B(f) => f,
            Val::O(f) => f,
        }
    }
}

/// A value that is rounded into a target precision on assignment: an exact value, or a
/// golden-backed operation (rug's "incomplete computation"). An operation holds its operands as
/// their query text (`op` renders it when the operation is built), so it borrows nothing.
#[derive(Debug)]
pub enum Pending<'a> {
    Exact(Val<'a>),
    Op(&'static str, String),
    Shift(Val<'a>, i64),
    /// `-(inner)`: rounded as the inner value with the mirrored mode, then negated (exact).
    Neg(Box<Pending<'a>>),
}

impl<'a> std::ops::Neg for Pending<'a> {
    type Output = Pending<'a>;
    fn neg(self) -> Pending<'a> {
        match self {
            Pending::Exact(Val::O(f)) => Pending::Exact(Val::O(-f)),
            Pending::Neg(b) => *b,
            p => Pending::Neg(Box::new(p)),
        }
    }
}

impl Pending<'_> {
    /// The golden operation `op` on `args` (rendered to query text now).
    pub fn op(op: &'static str, args: &[Arg]) -> Pending<'static> {
        Pending::Op(op, args_text(args))
    }
    /// rug's `square` on an incomplete value: `self * self`.
    pub fn square(self) -> Pending<'static> {
        let mut x = Float::new(pending_prec(&self).max(64));
        let t = self.assign_to(&mut x, Rnd::Nearest);
        assert_eq!(t, 0, "square of a rounded pending value");
        Pending::op("sqr", &[Arg::F(&x)])
    }
    /// Rounds the value into `f` (at `f.prec()`) with `r`; returns the ternary value.
    fn round_into(&self, f: &mut Float, r: Rnd) -> i32 {
        match self {
            Pending::Exact(x) => set_from(f, x.get(), r),
            Pending::Op(op, args) => golden_op_text(f, op, r, args),
            Pending::Shift(x, n) => mpfr::mul_2si(f, x.get(), *n, r),
            Pending::Neg(inner) => {
                let mr = match r {
                    Rnd::Up => Rnd::Down,
                    Rnd::Down => Rnd::Up,
                    x => x,
                };
                let t = inner.round_into(f, mr);
                f.neg = !f.neg;
                -t
            }
        }
    }
}

impl Assign for &Pending<'_> {
    fn assign_to(self, f: &mut Float, r: Rnd) -> i32 {
        self.round_into(f, r)
    }
}

pub trait Assign {
    /// Rounds `self` into `f` (at `f.prec()`) with `r`; returns the ternary value.
    fn assign_to(self, f: &mut Float, r: Rnd) -> i32;
}

fn set_from(f: &mut Float, x: &Float, r: Rnd) -> i32 {
    match x.kind {
        Kind::Reg => {
            let t = f.set_exact(x.neg, &x.m, x.exp - i64::from(x.prec), false, r);
            if f.exp > emax() && f.exp > x.exp {
                // carry out of the top binade (mpfr_set4's overflow)
                return overflow(f, r, x.neg);
            }
            t
        }
        k => {
            *f = Float::singular(f.prec, k, x.neg);
            if k == Kind::Nan {
                raise(FLAG_NAN);
            }
            0
        }
    }
}

impl Assign for Float {
    fn assign_to(self, f: &mut Float, r: Rnd) -> i32 {
        set_from(f, &self, r)
    }
}
impl Assign for &Float {
    fn assign_to(self, f: &mut Float, r: Rnd) -> i32 {
        set_from(f, self, r)
    }
}
impl Assign for Pending<'_> {
    fn assign_to(self, f: &mut Float, r: Rnd) -> i32 {
        self.round_into(f, r)
    }
}
impl Assign for Special {
    fn assign_to(self, f: &mut Float, _r: Rnd) -> i32 {
        let (k, neg) = match self {
            Special::Zero => (Kind::Zero, false),
            Special::NegZero => (Kind::Zero, true),
            Special::Infinity => (Kind::Inf, false),
            Special::NegInfinity => (Kind::Inf, true),
            Special::Nan => (Kind::Nan, false),
        };
        *f = Float::singular(f.prec, k, neg);
        0
    }
}
macro_rules! assign_int {
    ($($t:ty),*) => {$(
        impl Assign for $t {
            fn assign_to(self, f: &mut Float, r: Rnd) -> i32 {
                let v = i128::from(self);
                f.set_exact(v < 0, &Big::from_u128(v.unsigned_abs()), 0, false, r)
            }
        }
    )*};
}
assign_int!(i8, i16, i32, i64, u8, u16, u32, u64);
impl Assign for usize {
    fn assign_to(self, f: &mut Float, r: Rnd) -> i32 {
        f.set_exact(false, &Big::from_u64(self as u64), 0, false, r)
    }
}
impl Assign for isize {
    fn assign_to(self, f: &mut Float, r: Rnd) -> i32 {
        (self as i64).assign_to(f, r)
    }
}
impl Assign for f64 {
    fn assign_to(self, f: &mut Float, r: Rnd) -> i32 {
        if self.is_nan() {
            *f = Float::singular(f.prec, Kind::Nan, false);
            return 0;
        }
        if self.is_infinite() {
            *f = Float::singular(f.prec, Kind::Inf, self < 0.0);
            return 0;
        }
        let bits = self.to_bits();
        let neg = bits >> 63 == 1;
        let be = ((bits >> 52) & 0x7ff) as i64;
        let frac = bits & ((1u64 << 52) - 1);
        let (m, e) = if be == 0 { (frac, -1074) } else { (frac | (1u64 << 52), be - 1075) };
        f.set_exact(neg, &Big::from_u64(m), e, false, r)
    }
}
/// An unsigned integer given by borrowed little-endian 64-bit words (`mpz_import` of an `ExactNum`'s raw words).
pub struct Words<'a>(pub &'a [u64]);
impl Assign for Words<'_> {
    fn assign_to(self, f: &mut Float, r: Rnd) -> i32 {
        if self.0.is_zero() {
            *f = Float::singular(f.prec, Kind::Zero, false);
            return 0;
        }
        f.set_exact(false, self.0, 0, false, r)
    }
}
impl Assign for &Integer {
    fn assign_to(self, f: &mut Float, r: Rnd) -> i32 {
        f.set_exact(self.neg, &self.mag, 0, false, r)
    }
}
impl Assign for Integer {
    fn assign_to(self, f: &mut Float, r: Rnd) -> i32 {
        (&self).assign_to(f, r)
    }
}

macro_rules! binop {
    ($tr:ident, $m:ident, $op:literal) => {
        impl std::ops::$tr<&Float> for &Float {
            type Output = Pending<'static>;
            fn $m(self, o: &Float) -> Pending<'static> {
                Pending::op($op, &[Arg::F(self), Arg::F(o)])
            }
        }
        impl std::ops::$tr<Float> for &Float {
            type Output = Pending<'static>;
            fn $m(self, o: Float) -> Pending<'static> {
                Pending::op($op, &[Arg::F(self), Arg::F(&o)])
            }
        }
        /// rug: an owned left operand is reused, so the result has its precision (nearest).
        impl std::ops::$tr<&Float> for Float {
            type Output = Float;
            fn $m(self, o: &Float) -> Float {
                let mut r = Float::new(self.prec);
                let _ = golden_op(&mut r, $op, Rnd::Nearest, &[Arg::F(&self), Arg::F(o)]);
                r
            }
        }
        impl std::ops::$tr<Float> for Float {
            type Output = Float;
            fn $m(self, o: Float) -> Float {
                let mut r = Float::new(self.prec);
                let _ = golden_op(&mut r, $op, Rnd::Nearest, &[Arg::F(&self), Arg::F(&o)]);
                r
            }
        }
    };
}
binop!(Add, add, "add");
binop!(Sub, sub, "sub");
binop!(Mul, mul, "mul");
binop!(Div, div, "div");

macro_rules! intop {
    ($tr:ident, $m:ident, $op:literal) => {
        impl std::ops::$tr<&Integer> for &Float {
            type Output = Pending<'static>;
            fn $m(self, o: &Integer) -> Pending<'static> {
                let bits = u32::try_from(o.mag.bits()).unwrap_or(u32::MAX).max(2);
                Pending::op($op, &[Arg::F(self), Arg::F(&Float::with_val(bits, o))])
            }
        }
    };
}
intop!(Add, add, "add");
intop!(Sub, sub, "sub");
intop!(Div, div, "div");

macro_rules! opassign {
    ($tr:ident, $m:ident, $op:literal) => {
        impl std::ops::$tr<Float> for Float {
            fn $m(&mut self, o: Float) {
                let a = take(self);
                let _ = golden_op(self, $op, Rnd::Nearest, &[Arg::F(&a), Arg::F(&o)]);
            }
        }
        impl std::ops::$tr<&Float> for Float {
            fn $m(&mut self, o: &Float) {
                let a = take(self);
                let _ = golden_op(self, $op, Rnd::Nearest, &[Arg::F(&a), Arg::F(o)]);
            }
        }
        impl std::ops::$tr<Pending<'_>> for Float {
            fn $m(&mut self, o: Pending<'_>) {
                let mut b = Float::new(self.prec.max(pending_prec(&o)));
                let _ = o.assign_to(&mut b, Rnd::Nearest);
                let a = take(self);
                let _ = golden_op(self, $op, Rnd::Nearest, &[Arg::F(&a), Arg::F(&b)]);
            }
        }
    };
}
opassign!(AddAssign, add_assign, "add");
opassign!(SubAssign, sub_assign, "sub");
opassign!(MulAssign, mul_assign, "mul");
opassign!(DivAssign, div_assign, "div");

/// Moves `x`'s value out, leaving a zero of the same precision (the target of an in-place op).
fn take(x: &mut Float) -> Float {
    let p = x.prec;
    std::mem::replace(x, Float::new(p))
}

/// Precision that holds an exact pending value (shifts and exact copies); 64 for golden ops.
fn pending_prec(p: &Pending) -> u32 {
    match p {
        Pending::Exact(f) | Pending::Shift(f, _) => f.get().prec,
        Pending::Op(..) => 64,
        Pending::Neg(b) => pending_prec(b),
    }
}

macro_rules! shifts {
    ($($t:ty),*) => {$(
        impl std::ops::Shl<$t> for Float {
            type Output = Float;
            fn shl(mut self, n: $t) -> Float {
                self <<= n;
                self
            }
        }
        impl std::ops::Shr<$t> for Float {
            type Output = Float;
            fn shr(mut self, n: $t) -> Float {
                self >>= n;
                self
            }
        }
        impl<'a> std::ops::Shl<$t> for &'a Float {
            type Output = Pending<'a>;
            fn shl(self, n: $t) -> Pending<'a> {
                Pending::Shift(Val::B(self), shift_amount(n))
            }
        }
        impl<'a> std::ops::Shr<$t> for &'a Float {
            type Output = Pending<'a>;
            fn shr(self, n: $t) -> Pending<'a> {
                Pending::Shift(Val::B(self), -shift_amount(n))
            }
        }
        impl std::ops::ShlAssign<$t> for Float {
            fn shl_assign(&mut self, n: $t) {
                let x = take(self);
                let _ = mpfr::mul_2si(self, &x, shift_amount(n), Rnd::Nearest);
            }
        }
        impl std::ops::ShrAssign<$t> for Float {
            fn shr_assign(&mut self, n: $t) {
                let x = take(self);
                let _ = mpfr::mul_2si(self, &x, -shift_amount(n), Rnd::Nearest);
            }
        }
        impl std::ops::Shr<$t> for &Integer {
            type Output = Integer;
            fn shr(self, n: $t) -> Integer {
                let n = u64::try_from(n).or_fail(|e| format!("Integer shift: {e}"));
                // floor division by 2^n (GMP's mpz_fdiv_q_2exp, as rug's >>)
                let q = self.mag.shr(n);
                if self.neg && self.mag.low_nonzero(n) {
                    Integer::make(true, q.add(&Big::from_u64(1)))
                } else {
                    Integer::make(self.neg, q)
                }
            }
        }
        impl PartialEq<$t> for Float {
            fn eq(&self, o: &$t) -> bool {
                *self == Float::with_val(128, *o)
            }
        }
        impl PartialEq<$t> for Integer {
            fn eq(&self, o: &$t) -> bool {
                let v = i64::try_from(*o).or_fail(|e| format!("Integer eq: {e}"));
                // compared in place: sign, then the magnitude as at most one limb
                self.neg == (v < 0) && self.mag.bits() <= 64 && self.mag.first().copied().unwrap_or(0) == v.unsigned_abs()
            }
        }
        impl std::ops::Div<$t> for Float {
            type Output = Float;
            fn div(mut self, o: $t) -> Float {
                self /= o;
                self
            }
        }
        impl std::ops::DivAssign<$t> for Float {
            fn div_assign(&mut self, o: $t) {
                let d = Float::with_val(128, o);
                let a = take(self);
                let _ = golden_op(self, "div", Rnd::Nearest, &[Arg::F(&a), Arg::F(&d)]);
            }
        }
        impl std::ops::Add<$t> for Integer {
            type Output = Integer;
            fn add(self, o: $t) -> Integer {
                self.add_signed(&Integer::from(i64::try_from(o).or_fail(|e| format!("Integer add: {e}"))), false)
            }
        }
        impl std::ops::Sub<$t> for Integer {
            type Output = Integer;
            fn sub(self, o: $t) -> Integer {
                self.add_signed(&Integer::from(i64::try_from(o).or_fail(|e| format!("Integer sub: {e}"))), true)
            }
        }
    )*};
}
shifts!(i32, i64, isize, u32, u64, usize);

fn shift_amount<T: TryInto<i64>>(n: T) -> i64 {
    n.try_into()
        .or_fail(|_| "shift amount out of i64".to_owned())
}

// ---------------------------------------------------------------- signed big integers (rug::Integer subset)

#[derive(Debug, Default, PartialEq, Eq)]
pub struct Integer {
    neg: bool,
    mag: Big,
}

impl Integer {
    pub fn new() -> Integer {
        Integer::default()
    }
    /// `(-1)^neg * mag` (a zero is never negative).
    pub fn from_big(neg: bool, mag: Big) -> Integer {
        Integer::make(neg, mag)
    }
    fn make(neg: bool, mag: Big) -> Integer {
        let neg = neg && !mag.is_zero();
        Integer { neg, mag }
    }
    pub fn magnitude(&self) -> &Big {
        &self.mag
    }
    pub fn is_negative(&self) -> bool {
        self.neg
    }
    pub fn is_zero(&self) -> bool {
        self.mag.is_zero()
    }
    pub fn abs(self) -> Integer {
        Integer::make(false, self.mag)
    }
    pub fn as_raw_mut(&mut self) -> &mut Integer {
        self
    }
    pub fn significant_bits(&self) -> u32 {
        u32::try_from(self.mag.bits()).unwrap_or(u32::MAX)
    }
    /// Low 64 bits of the two's complement value (rug's `to_u64_wrapping`).
    pub fn to_u64_wrapping(&self) -> u64 {
        let lo = self.mag.limbs().first().copied().unwrap_or(0);
        if self.neg {
            lo.wrapping_neg()
        } else {
            lo
        }
    }
    fn add_signed(&self, o: &Integer, flip: bool) -> Integer {
        let oneg = o.neg != flip;
        if self.neg == oneg {
            return Integer::make(self.neg, self.mag.add(&o.mag));
        }
        if self.mag >= o.mag {
            Integer::make(self.neg, self.mag.sub(&o.mag))
        } else {
            Integer::make(oneg, o.mag.sub(&self.mag))
        }
    }
    /// Index of the lowest set bit at or above `start` (nonnegative values; `None` for zero).
    pub fn find_one(&self, start: u32) -> Option<u32> {
        assert!(!self.neg, "find_one on a negative Integer");
        let m = self.mag.shr(u64::from(start));
        (!m.is_zero()).then(|| start + u32::try_from(m.trailing_zeros()).unwrap_or(u32::MAX))
    }
    /// `self^e mod m` for `e >= 0`, `m > 0` (rug's `pow_mod`; nonnegative result).
    pub fn pow_mod(&self, e: &Integer, m: &Integer) -> Result<Integer, ()> {
        if e.neg || m.mag.is_zero() {
            return Err(());
        }
        let mut r = Big::from_u64(1).divrem(&m.mag).1;
        let base = self.mag.divrem(&m.mag).1;
        let base = if self.neg && !base.is_zero() { m.mag.sub(&base) } else { base };
        for i in (0..e.mag.bits()).rev() {
            r = r.mul(&r).divrem(&m.mag).1;
            if e.mag.bit(i) {
                r = r.mul(&base).divrem(&m.mag).1;
            }
        }
        Ok(Integer::make(false, r))
    }
    pub fn pow<N: Into<u64>>(&self, n: N) -> Integer {
        let n: u64 = n.into();
        // left to right: square, then multiply by the base for each set bit
        let mut r = Big::from_u64(1);
        for i in (0..64 - n.leading_zeros()).rev() {
            r = r.mul(&r);
            if (n >> i) & 1 == 1 {
                r = r.mul(&self.mag);
            }
        }
        Integer::make(self.neg && n % 2 == 1, r)
    }
    /// Truncated remainder (sign of the dividend), as GMP's `mpz_tdiv_r` / rug's `%`.
    pub fn rem_trunc(&self, d: &Integer) -> Integer {
        Integer::make(self.neg, self.mag.divrem(&d.mag).1)
    }
}

impl From<i64> for Integer {
    fn from(v: i64) -> Integer {
        Integer::make(v < 0, Big::from_u64(v.unsigned_abs()))
    }
}
impl From<u64> for Integer {
    fn from(v: u64) -> Integer {
        Integer::make(false, Big::from_u64(v))
    }
}
impl std::ops::Rem<&Integer> for Integer {
    type Output = Integer;
    fn rem(self, d: &Integer) -> Integer {
        self.rem_trunc(d)
    }
}
impl std::ops::Rem<Integer> for Integer {
    type Output = Integer;
    fn rem(self, d: Integer) -> Integer {
        self.rem_trunc(&d)
    }
}
impl From<i32> for Integer {
    fn from(v: i32) -> Integer {
        Integer::from(i64::from(v))
    }
}
impl std::ops::Add<&Integer> for &Integer {
    type Output = Integer;
    fn add(self, o: &Integer) -> Integer {
        self.add_signed(o, false)
    }
}
impl std::ops::Sub<&Integer> for &Integer {
    type Output = Integer;
    fn sub(self, o: &Integer) -> Integer {
        self.add_signed(o, true)
    }
}
impl std::ops::Mul<&Integer> for &Integer {
    type Output = Integer;
    fn mul(self, o: &Integer) -> Integer {
        Integer::make(self.neg != o.neg, self.mag.mul(&o.mag))
    }
}
impl std::ops::Neg for Integer {
    type Output = Integer;
    fn neg(self) -> Integer {
        Integer::make(!self.neg, self.mag)
    }
}
impl std::ops::Shl<u32> for &Integer {
    type Output = Integer;
    fn shl(self, n: u32) -> Integer {
        Integer::make(self.neg, self.mag.shl(u64::from(n)))
    }
}
impl std::ops::Shl<u32> for Integer {
    type Output = Integer;
    fn shl(self, n: u32) -> Integer {
        Integer::make(self.neg, self.mag.shl(u64::from(n)))
    }
}
impl std::ops::Mul<&Integer> for &Float {
    type Output = Pending<'static>;
    fn mul(self, o: &Integer) -> Pending<'static> {
        let bits = u32::try_from(o.mag.bits()).unwrap_or(u32::MAX).max(2);
        let of = Float::with_val(bits, o);
        Pending::op("mul", &[Arg::F(self), Arg::F(&of)])
    }
}
impl std::ops::Mul<&Integer> for Float {
    type Output = Pending<'static>;
    fn mul(self, o: &Integer) -> Pending<'static> {
        &self * o
    }
}
impl std::ops::Rem<&Integer> for &Integer {
    type Output = Integer;
    fn rem(self, d: &Integer) -> Integer {
        self.rem_trunc(d)
    }
}
impl std::ops::Shl<u64> for &Integer {
    type Output = Integer;
    fn shl(self, n: u64) -> Integer {
        Integer::make(self.neg, self.mag.shl(n))
    }
}
impl std::ops::ShlAssign<u64> for Integer {
    fn shl_assign(&mut self, n: u64) {
        self.mag = self.mag.shl(n);
    }
}
impl std::ops::AddAssign<&Integer> for Integer {
    fn add_assign(&mut self, o: &Integer) {
        *self = self.add_signed(o, false);
    }
}
impl std::ops::AddAssign<Integer> for Integer {
    fn add_assign(&mut self, o: Integer) {
        *self = self.add_signed(&o, false);
    }
}
impl std::ops::AddAssign<u64> for Integer {
    fn add_assign(&mut self, o: u64) {
        *self = self.add_signed(&Integer::from(o), false);
    }
}
impl TryFrom<&Integer> for i128 {
    type Error = String;
    fn try_from(v: &Integer) -> Result<i128, String> {
        if v.mag.bits() > 127 {
            return Err(format!("{v} does not fit i128"));
        }
        let l = v.mag.limbs();
        let m = i128::from(l.first().copied().unwrap_or(0))
            | (i128::from(l.get(1).copied().unwrap_or(0)) << 64);
        Ok(if v.neg { -m } else { m })
    }
}
macro_rules! int_ord {
    ($($t:ty),*) => {$(
        impl PartialOrd<$t> for Integer {
            fn partial_cmp(&self, o: &$t) -> Option<Ordering> {
                self.partial_cmp(&Integer::from(i64::try_from(*o).or_fail(|e| format!("Integer cmp: {e}"))))
            }
        }
    )*};
}
int_ord!(i32, i64, u32, u64);
impl PartialOrd for Integer {
    fn partial_cmp(&self, o: &Integer) -> Option<Ordering> {
        Some(match (self.neg, o.neg) {
            (false, true) => Ordering::Greater,
            (true, false) => Ordering::Less,
            (false, false) => self.mag.cmp(&o.mag),
            (true, true) => o.mag.cmp(&self.mag),
        })
    }
}
impl fmt::Display for Integer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}0x{}",
            if self.neg { "-" } else { "" },
            self.mag.to_hex()
        )
    }
}

// ---------------------------------------------------------------- golden-backed operations

fn arg_text(a: &Arg) -> String {
    match a {
        Arg::F(f) => f.canon(),
        Arg::I(i) => format!("i{i}"),
        Arg::U(u) => format!("u{u}"),
        Arg::S(s) => format!("s{}", s.len()) + ":" + s,
    }
}

/// The operands' query text: each operand preceded by a space.
fn args_text(args: &[Arg]) -> String {
    let mut t = String::new();
    for a in args {
        t.push(' ');
        t.push_str(&arg_text(a));
    }
    t
}

/// Parsed golden value at precision `prec`: a float or a beyond-i64 exponent marker.
enum GVal {
    F(Float),
    Huge(bool),
    Tiny(bool),
}

fn parse_val(s: &str, prec: u32) -> GVal {
    match s {
        "nan" => return GVal::F(Float::singular(prec, Kind::Nan, false)),
        "+inf" | "-inf" => return GVal::F(Float::singular(prec, Kind::Inf, s.starts_with('-'))),
        "+0" | "-0" => return GVal::F(Float::singular(prec, Kind::Zero, s.starts_with('-'))),
        "+huge" | "-huge" => return GVal::Huge(s.starts_with('-')),
        "+tiny" | "-tiny" => return GVal::Tiny(s.starts_with('-')),
        _ => {}
    }
    let neg = s.starts_with('-');
    let (m, e) = s[1..]
        .split_once('p')
        .or_fail(|_| format!("golden value {s}"));
    let e: i64 = e.parse().or_fail(|x| format!("golden exponent {e}: {x}"));
    GVal::F(Float::from_parts(prec, neg, &Big::from_hex(m), e))
}

/// One golden query; returns the answers `(value, ternary)` (one pair, two for sin_cos/sinh_cosh).
fn query(op: &str, prec: u32, rnd: Rnd, args: &str) -> Vec<(GVal, i32)> {
    let q = format!("{op} {prec} {}{args}", rnd.tag());
    let ans = gold::answer(&q);
    let parts: Vec<&str> = ans.split(' ').collect();
    assert!(
        parts.len().is_multiple_of(2) && !parts.is_empty(),
        "golden answer {ans} for {op}"
    );
    parts
        .chunks(2)
        .map(|c| {
            let t: i32 = c[1]
                .parse()
                .or_fail(|x| format!("golden ternary {}: {x}", c[1]));
            (parse_val(c[0], prec), t)
        })
        .collect()
}

/// Stores a golden answer in `r` and applies MPFR's exponent-range rules (`mpfr_check_range`).
fn land(r: &mut Float, v: GVal, t: i32, rnd: Rnd) -> i32 {
    match v {
        GVal::F(f) => {
            if f.is_nan() {
                raise(FLAG_NAN);
            }
            *r = f;
            mpfr::check_range(r, t, rnd)
        }
        GVal::Huge(neg) => overflow(r, rnd, neg),
        GVal::Tiny(neg) => {
            let rnd = if rnd == Rnd::Nearest { Rnd::Zero } else { rnd };
            underflow(r, rnd, neg)
        }
    }
}

fn golden_op(r: &mut Float, op: &'static str, rnd: Rnd, args: &[Arg]) -> i32 {
    golden_op_text(r, op, rnd, &args_text(args))
}

fn golden_op_text(r: &mut Float, op: &str, rnd: Rnd, args: &str) -> i32 {
    let mut a = query(op, r.prec, rnd, args);
    assert_eq!(a.len(), 1, "golden {op}: one answer expected");
    let (v, t) = a.remove(0);
    land(r, v, t, rnd)
}

fn setmin(x: &mut Float, e: i64) {
    let p = u64::from(x.prec);
    x.kind = Kind::Reg;
    x.m = Big::from_u64(1).shl(p - 1);
    x.exp = e;
}
fn setmax(x: &mut Float, e: i64) {
    let p = u64::from(x.prec);
    x.kind = Kind::Reg;
    x.m = Big::from_u64(1).shl(p).sub(&Big::from_u64(1));
    x.exp = e;
}
fn powerof2_raw(x: &Float) -> bool {
    x.m.trailing_zeros() == u64::from(x.prec) - 1
}

/// `mpfr_underflow`.
fn underflow(x: &mut Float, rnd: Rnd, neg: bool) -> i32 {
    let inex = if rnd.like_rndz(neg) {
        *x = Float::singular(x.prec, Kind::Zero, neg);
        -1
    } else {
        setmin(x, emin());
        1
    };
    x.neg = neg;
    raise(FLAG_INEXACT | FLAG_UNDERFLOW);
    signed(neg, inex)
}
/// `mpfr_overflow`.
fn overflow(x: &mut Float, rnd: Rnd, neg: bool) -> i32 {
    let inex = if rnd.like_rndz(neg) {
        setmax(x, emax());
        -1
    } else {
        *x = Float::singular(x.prec, Kind::Inf, neg);
        1
    };
    x.neg = neg;
    raise(FLAG_INEXACT | FLAG_OVERFLOW);
    signed(neg, inex)
}

fn next_to_zero(x: &mut Float) {
    match x.kind {
        Kind::Inf => {
            x.kind = Kind::Reg;
            setmax(x, emax());
        }
        Kind::Zero => {
            x.neg = !x.neg;
            setmin(x, emin());
        }
        Kind::Nan => {}
        Kind::Reg => {
            let p = u64::from(x.prec);
            x.m = x.m.sub(&Big::from_u64(1));
            if x.m.bits() < p {
                if x.exp == emin() {
                    *x = Float::singular(x.prec, Kind::Zero, x.neg);
                } else {
                    x.exp -= 1;
                    // was a power of two: one ulp below it is all ones at p bits
                    x.m = Big::from_u64(1).shl(p).sub(&Big::from_u64(1));
                }
            }
        }
    }
}
fn next_to_inf(x: &mut Float) {
    match x.kind {
        Kind::Zero => setmin(x, emin()),
        Kind::Reg => {
            let p = u64::from(x.prec);
            x.m = x.m.add(&Big::from_u64(1));
            if x.m.bits() > p {
                if x.exp == emax() {
                    *x = Float::singular(x.prec, Kind::Inf, x.neg);
                } else {
                    x.exp += 1;
                    x.m = Big::from_u64(1).shl(p - 1);
                }
            }
        }
        _ => {}
    }
}

/// MPFR's API, safe: same names and argument order as `gmp_mpfr_sys::mpfr`, references for pointers.
pub mod mpfr {
    pub use super::Rnd;
    use super::*;

    pub type PrecT = i64;
    pub type ExpT = i64;

    pub fn get_emin() -> i64 {
        emin()
    }
    pub fn get_emax() -> i64 {
        emax()
    }
    pub fn get_emin_min() -> i64 {
        EMIN_MIN
    }
    pub fn get_emax_max() -> i64 {
        EMAX_MAX
    }
    pub fn set_emin(e: i64) -> i32 {
        if (EMIN_MIN..=EMAX_MAX).contains(&e) {
            EMIN.with(|c| c.set(e));
            0
        } else {
            1
        }
    }
    pub fn set_emax(e: i64) -> i32 {
        if (EMIN_MIN..=EMAX_MAX).contains(&e) {
            EMAX.with(|c| c.set(e));
            0
        } else {
            1
        }
    }
    pub fn clear_flags() {
        FLAGS.with(|c| c.set(0));
    }
    pub fn clear_overflow() {
        FLAGS.with(|c| c.set(c.get() & !FLAG_OVERFLOW));
    }
    pub fn overflow_p() -> i32 {
        i32::from(FLAGS.with(Cell::get) & FLAG_OVERFLOW != 0)
    }
    pub fn underflow_p() -> i32 {
        i32::from(FLAGS.with(Cell::get) & FLAG_UNDERFLOW != 0)
    }

    /// `mpfr_check_range`.
    pub fn check_range(x: &mut Float, t: i32, rnd: Rnd) -> i32 {
        if x.kind == Kind::Reg {
            let e = x.exp;
            if e < emin() {
                let rnd = if rnd == Rnd::Nearest
                    && (e + 1 < emin() || (powerof2_raw(x) && if x.neg { t <= 0 } else { t >= 0 }))
                {
                    Rnd::Zero
                } else {
                    rnd
                };
                let neg = x.neg;
                return underflow(x, rnd, neg);
            }
            if e > emax() {
                let neg = x.neg;
                return overflow(x, rnd, neg);
            }
        } else if t != 0 && x.kind == Kind::Inf {
            raise(FLAG_OVERFLOW);
        }
        if t != 0 {
            raise(FLAG_INEXACT);
        }
        t
    }

    /// `mpfr_set`.
    pub fn set(r: &mut Float, x: &Float, rnd: Rnd) -> i32 {
        set_from(r, x, rnd)
    }

    /// `mpfr_prec_round`.
    pub fn prec_round(x: &mut Float, prec: PrecT, rnd: Rnd) -> i32 {
        let p = u32::try_from(prec).or_fail(|e| format!("prec_round {prec}: {e}"));
        let old = std::mem::replace(x, Float::new(p));
        if old.kind != Kind::Reg {
            *x = Float::singular(p, old.kind, old.neg);
            return 0;
        }
        let t = x.set_exact(old.neg, &old.m, old.exp - i64::from(old.prec), false, rnd);
        if x.exp > old.exp && old.exp == emax() {
            let neg = old.neg;
            return overflow(x, rnd, neg);
        }
        t
    }

    /// `mpfr_mul_2si(x, x, n, rnd)`: in place (MPFR's aliased call).
    pub fn mul_2si_self(x: &mut Float, n: i64, rnd: Rnd) -> i32 {
        let a = take(x);
        mul_2si(x, &a, n, rnd)
    }
    /// `mpfr_div_ui(x, x, n, rnd)`: in place.
    pub fn div_ui_self(x: &mut Float, n: u64, rnd: Rnd) -> i32 {
        let a = take(x);
        div_ui(x, &a, n, rnd)
    }
    /// `mpfr_log(x, x, rnd)`: in place.
    pub fn log_self(x: &mut Float, rnd: Rnd) -> i32 {
        let a = take(x);
        log(x, &a, rnd)
    }

    /// `mpfr_mul_2si`.
    pub fn mul_2si(y: &mut Float, x: &Float, n: i64, rnd: Rnd) -> i32 {
        if x.kind != Kind::Reg {
            return set(y, x, rnd);
        }
        let exp = x.exp;
        let inexact = y.set_exact(x.neg, &x.m, -i64::from(x.prec), false, rnd);
        // set_exact put y.exp = (0 or 1 on carry); MPFR_SETRAW keeps exp and adds the carry
        let carry = y.exp;
        let exp = exp + carry;
        let (emn, emx) = (emin(), emax());
        if n >= 0 && (emx < EMIN_MIN + n || exp > emx - n) {
            let neg = y.neg;
            return overflow(y, rnd, neg);
        } else if n < 0 && (emn > EMAX_MAX + n || exp < emn - n) {
            let neg = y.neg;
            let rnd = if rnd == Rnd::Nearest
                && (emn > EMAX_MAX + (n + 1)
                    || exp < emn - (n + 1)
                    || ((if neg { inexact <= 0 } else { inexact >= 0 }) && powerof2_raw(y)))
            {
                Rnd::Zero
            } else {
                rnd
            };
            return underflow(y, rnd, neg);
        }
        y.exp = exp + n;
        inexact
    }

    /// `mpfr_subnormalize`.
    pub fn subnormalize(y: &mut Float, old_inexact: i32, rnd: Rnd) -> i32 {
        let p = i64::from(y.prec);
        if y.kind != Kind::Reg || y.exp >= emin() + p - 1 {
            return old_inexact;
        }
        raise(FLAG_UNDERFLOW);
        let sign = if y.neg { -1 } else { 1 };
        if y.exp == emin() {
            if powerof2_raw(y) {
                return old_inexact;
            }
            let pu = u64::from(y.prec);
            let set_min = if rnd == Rnd::Nearest {
                let rb = y.m.bit(pu - 2);
                if !rb {
                    true
                } else if y.m.low_nonzero(pu - 2) {
                    false
                } else {
                    (old_inexact > 0 && sign > 0) || (old_inexact < 0 && sign < 0)
                }
            } else {
                rnd.like_rndz(y.neg)
            };
            if set_min {
                setmin(y, emin());
                return -sign;
            }
            setmin(y, emin() + 1);
            return sign;
        }
        let q = u32::try_from(y.exp - emin() + 1).or_fail(|e| format!("subnormalize q: {e}"));
        let mut dest = Float::new(q);
        let r = round_mag(&y.m, y.exp - p, false, q, rnd, y.neg);
        dest.kind = Kind::Reg;
        dest.neg = y.neg;
        dest.m = r.m;
        dest.exp = r.e2 + i64::from(q);
        let mut inexact = signed(y.neg, r.dir) * if r.tie { 2 } else { 1 };
        if old_inexact != 0 {
            if rnd == Rnd::Nearest && r.tie {
                if (inexact > 0) == (old_inexact > 0) {
                    if (inexact > 0) == (sign > 0) {
                        next_to_zero(&mut dest);
                    } else {
                        next_to_inf(&mut dest);
                    }
                    inexact = -inexact;
                }
            } else if inexact == 0 {
                inexact = old_inexact;
            }
        }
        let inex2 = set(y, &dest, rnd);
        assert_eq!(inex2, 0, "subnormalize: dest does not fit y");
        inexact
    }

    pub fn nextabove(x: &mut Float) {
        if x.is_nan() {
            raise(FLAG_NAN);
        } else if x.neg {
            next_to_zero(x);
        } else {
            next_to_inf(x);
        }
    }
    pub fn nextbelow(x: &mut Float) {
        if x.is_nan() {
            raise(FLAG_NAN);
        } else if x.neg {
            next_to_inf(x);
        } else {
            next_to_zero(x);
        }
    }
    pub fn get_exp(x: &Float) -> ExpT {
        x.exp
    }
    /// `mpfr_min_prec`: bits needed for the mantissa (0 for singular values).
    pub fn min_prec(x: &Float) -> PrecT {
        if x.kind != Kind::Reg {
            return 0;
        }
        i64::from(x.prec) - x.m.trailing_zeros() as i64
    }
    /// `mpfr_copysign`: `x` with the sign of `y` (exact up to the rounding to `r`'s precision).
    pub fn copysign(r: &mut Float, x: &Float, y: &Float, rnd: Rnd) -> i32 {
        let t = set(r, x, rnd);
        if r.neg == y.neg {
            return t;
        }
        // the magnitude rounds the same whatever its sign, so flip the result and the ternary
        // (with the mirrored mode when the rounding was directed)
        let mr = match rnd {
            Rnd::Up => Rnd::Down,
            Rnd::Down => Rnd::Up,
            m => m,
        };
        let t = set(r, x, mr);
        r.neg = !r.neg;
        -t
    }

    fn f(x: &Float) -> Arg<'_> {
        Arg::F(x)
    }

    macro_rules! unary {
        ($($name:ident),*) => {$(
            pub fn $name(r: &mut Float, x: &Float, rnd: Rnd) -> i32 {
                golden_op(r, stringify!($name), rnd, &[f(x)])
            }
        )*};
    }
    unary!(
        sqrt, cbrt, exp, exp2, exp10, expm1, log, log2, log10, log1p, sin, cos, tan, asin, acos,
        atan, sinh, cosh, tanh, asinh, acosh, atanh, gamma, lngamma, digamma, erf, erfc, eint, ai
    );
    macro_rules! binary {
        ($($name:ident),*) => {$(
            pub fn $name(r: &mut Float, x: &Float, y: &Float, rnd: Rnd) -> i32 {
                golden_op(r, stringify!($name), rnd, &[f(x), f(y)])
            }
        )*};
    }
    binary!(add, sub, mul, div, atan2, hypot, pow, remainder, gamma_inc);
    pub fn fma(r: &mut Float, a: &Float, b: &Float, c: &Float, rnd: Rnd) -> i32 {
        golden_op(r, "fma", rnd, &[f(a), f(b), f(c)])
    }
    pub fn pow_ui(r: &mut Float, x: &Float, n: u64, rnd: Rnd) -> i32 {
        golden_op(r, "pow_ui", rnd, &[f(x), Arg::U(n)])
    }
    pub fn pow_si(r: &mut Float, x: &Float, n: i64, rnd: Rnd) -> i32 {
        golden_op(r, "pow_si", rnd, &[f(x), Arg::I(n)])
    }
    pub fn rootn_ui(r: &mut Float, x: &Float, k: u64, rnd: Rnd) -> i32 {
        golden_op(r, "rootn_ui", rnd, &[f(x), Arg::U(k)])
    }
    pub fn add_ui(r: &mut Float, x: &Float, n: u64, rnd: Rnd) -> i32 {
        golden_op(r, "add_ui", rnd, &[f(x), Arg::U(n)])
    }
    pub fn div_ui(r: &mut Float, x: &Float, n: u64, rnd: Rnd) -> i32 {
        golden_op(r, "div_ui", rnd, &[f(x), Arg::U(n)])
    }
    pub fn ui_div(r: &mut Float, n: u64, x: &Float, rnd: Rnd) -> i32 {
        golden_op(r, "ui_div", rnd, &[Arg::U(n), f(x)])
    }
    pub fn jn(r: &mut Float, n: i64, x: &Float, rnd: Rnd) -> i32 {
        golden_op(r, "jn", rnd, &[Arg::I(n), f(x)])
    }
    pub fn yn(r: &mut Float, n: i64, x: &Float, rnd: Rnd) -> i32 {
        golden_op(r, "yn", rnd, &[Arg::I(n), f(x)])
    }
    pub fn const_pi(r: &mut Float, rnd: Rnd) -> i32 {
        golden_op(r, "const_pi", rnd, &[])
    }
    pub fn const_log2(r: &mut Float, rnd: Rnd) -> i32 {
        golden_op(r, "const_log2", rnd, &[])
    }
    pub fn const_euler(r: &mut Float, rnd: Rnd) -> i32 {
        golden_op(r, "const_euler", rnd, &[])
    }
    /// `mpfr_set_str` (base 2 or 10); returns 0 on success, -1 on a malformed string.
    pub fn set_str(r: &mut Float, s: &str, base: i32, rnd: Rnd) -> i32 {
        let _ = golden_op(
            r,
            if base == 2 { "parse2" } else { "parse10" },
            rnd,
            &[Arg::S(s)],
        );
        0
    }
    /// `mpfr_sin_cos` / `mpfr_sinh_cosh`: two results; MPFR's combined ternary
    /// `INEX(s, c) = s' | c' << 2` with `x' = 0, 1, 2` for exact, positive, negative.
    fn pair(op: &'static str, a: &mut Float, b: &mut Float, x: &Float, rnd: Rnd) -> i32 {
        assert_eq!(a.prec, b.prec, "{op}: the golden pair is at one precision");
        let mut ans = query(op, a.prec, rnd, &args_text(&[f(x)]));
        assert_eq!(ans.len(), 2, "golden {op}: two answers expected");
        let (vb, tb) = ans.remove(1);
        let (va, ta) = ans.remove(0);
        let ia = land(a, va, ta, rnd);
        let ib = land(b, vb, tb, rnd);
        let code = |t: i32| match t.cmp(&0) {
            Ordering::Equal => 0,
            Ordering::Greater => 1,
            Ordering::Less => 2,
        };
        code(ia) | (code(ib) << 2)
    }
    pub fn sin_cos(s: &mut Float, c: &mut Float, x: &Float, rnd: Rnd) -> i32 {
        pair("sin_cos", s, c, x, rnd)
    }
    pub fn sinh_cosh(s: &mut Float, c: &mut Float, x: &Float, rnd: Rnd) -> i32 {
        pair("sinh_cosh", s, c, x, rnd)
    }
}

// ---------------------------------------------------------------- rationals (GMP mpq subset)

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

/// A reduced fraction with `i64` parts; the arithmetic comes from the golden stream.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Rational {
    n: i64,
    d: i64,
}

/// One `Rational` operation waiting for its golden answer.
pub struct RationalOp(&'static str, Rational, Rational);

impl Rational {
    pub fn numer(&self) -> IntPart {
        IntPart(self.n)
    }
    pub fn denom(&self) -> IntPart {
        IntPart(self.d)
    }
}

/// A numerator or denominator (rug's `Integer` view).
pub struct IntPart(i64);
impl IntPart {
    pub fn to_i64(&self) -> Option<i64> {
        Some(self.0)
    }
}

impl From<(i64, i64)> for Rational {
    fn from((n, d): (i64, i64)) -> Rational {
        assert!(d != 0, "Rational with zero denominator");
        let g = gcd(u128::from(n.unsigned_abs()), u128::from(d.unsigned_abs())).max(1);
        let s = if d < 0 { -1 } else { 1 };
        let g = i64::try_from(g).or_fail(|e| format!("gcd: {e}"));
        Rational {
            n: s * n / g,
            d: s * d / g,
        }
    }
}
impl From<(i32, i32)> for Rational {
    fn from((n, d): (i32, i32)) -> Rational {
        Rational::from((i64::from(n), i64::from(d)))
    }
}
impl From<RationalOp> for Rational {
    fn from(op: RationalOp) -> Rational {
        let RationalOp(name, a, b) = op;
        let q = format!("{name} 0 N q{}/{} q{}/{}", a.n, a.d, b.n, b.d);
        let ans = gold::answer(&q);
        let (n, d) = ans
            .strip_prefix('q')
            .and_then(|x| x.split_once('/'))
            .or_fail(|_| format!("golden {name}: {ans}"));
        let n: i64 = n.parse().or_fail(|e| format!("golden numerator {n}: {e}"));
        let d: i64 = d
            .parse()
            .or_fail(|e| format!("golden denominator {d}: {e}"));
        Rational { n, d }
    }
}
impl fmt::Display for Rational {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}/{}", self.n, self.d)
    }
}
macro_rules! qop {
    ($tr:ident, $m:ident, $op:literal) => {
        impl std::ops::$tr<&Rational> for &Rational {
            type Output = RationalOp;
            fn $m(self, o: &Rational) -> RationalOp {
                RationalOp($op, *self, *o)
            }
        }
    };
}
qop!(Add, add, "q_add");
qop!(Sub, sub, "q_sub");
qop!(Mul, mul, "q_mul");
qop!(Div, div, "q_div");
