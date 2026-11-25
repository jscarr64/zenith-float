//! All seven modes for the Jacobi and Airy families (Spence spot-check at cc85211, items 2 and 3). Their grid digests
//! in ext.rs pin ToEven and Down only, and the old pins do not move (Q2). Here every mode is checked against an
//! oracle:
//! - `ai`: MPFR `mpfr_ai`, strict, through `mpfr_reference_with`.
//! - Every function, `ai` included: [`bracket_expect`]. It takes `f` at `PH` = rp(p, 64) + 128 bits in Down and Up,
//!   which bracket the value. When the two are equal and exact, that is the value. Otherwise they are adjacent at
//!   `PH`, and the value lies strictly inside. No point of the p grid or of its midpoints lies strictly inside, so the
//!   midpoint of the bracket rounds like the value in every mode. MPFR rounds that midpoint (`mpfr_reference_with`,
//!   so None is toward zero at rp(p, 64) and ToOdd is nearest with ties to odd).
//!
//! The bracket is `f` itself at a higher precision, so it checks the modes against each other and against `f`'s own
//! error bound. It does not check `f`'s algorithm (the digests are mpmath-verified for that), except for `ai`, which
//! MPFR checks.
//!
//! A bracket that is the largest value on one side and +-Inf on the other is an overflow, checked against the
//! M27-OVF ruling ([`Expect::Overflow`]); a subnormal bracket rounds at `p` on that precision's subnormal grid.
//!
//! Inputs: the digest's operands, which avoid Bi and Bi' overflow. Under None, ToZero and truncating Up/Down, Bi and
//! Bi' overflow are the ZF-M27 overflow rows (F-2), and are not pinned here.

use crate::mpfr::common::{conv_str_to_mpfr_compat, mpfr_reference_with};
use crate::mpfr::mp::mpfr;
use crate::mpfr::mp::Float;
use crate::mpfr::retry_budget::use_full_exponent_range;
use zenith_float_num::{
    Consts, Error, ExactNum, Radix, RoundingMode, EXPONENT_MIN, INF_NEG, INF_POS,
};

type Res = Result<(), Box<dyn std::error::Error>>;
type E = Box<dyn std::error::Error>;

const MODES: [RoundingMode; 7] = [
    RoundingMode::ToEven,
    RoundingMode::ToOdd,
    RoundingMode::ToZero,
    RoundingMode::FromZero,
    RoundingMode::Up,
    RoundingMode::Down,
    RoundingMode::None,
];

/// What `f(p, rm)` must be, for each of [`MODES`].
#[derive(Debug)]
pub enum Expect {
    /// The MPFR value and the inexact flag.
    Strict(Float, bool),
    /// A subnormal bracket (`Bracket::grid`), with no MPFR counterpart: the value rounds at `p` as the bracket's
    /// midpoint does, on the subnormal grid of `p` (`mismatch`).
    Grid,
    /// NaN, or an infinity of this sign.
    Special(Option<bool>),
    /// An overflow of this sign (negative: true): the largest value on the toward-zero side, flagged inexact, and
    /// +-Inf otherwise (ZF-M27 M27-OVF; `None` at the storage precision rp(p, 64)).
    Overflow(bool),
}

/// `x` as an MPFR number holding all its bits; `None` for a subnormal (below MPFR's range).
fn to_mpfr(x: &ExactNum, cc: &mut Consts) -> Result<Option<Float>, E> {
    if x.is_subnormal() {
        return Ok(None);
    }
    let bits = u32::try_from(x.mantissa_max_bit_len().ok_or("not finite")?)?;
    let text = conv_str_to_mpfr_compat(x.format(Radix::Bin, RoundingMode::ToEven, cc)?)?;
    let parsed = Float::parse_radix(text, 2)?;
    let (f, ord) = Float::with_val_round(bits, parsed, crate::mpfr::mp::Round::Nearest);
    if ord != std::cmp::Ordering::Equal {
        return Err(format!("{x:?} is not exact in {bits} bits").into());
    }
    Ok(Some(f))
}

fn set_from(src: &Float) -> impl Fn(&mut Float, mpfr::Rnd) -> i32 + '_ {
    move |r, rnd| mpfr::set(r.as_raw_mut(), src.as_raw(), rnd)
}

/// The expected results of `f` at `p`, one per mode of [`MODES`], and the subnormal bracket when there is one.
pub struct Bracket {
    pub want: Vec<Expect>,
    pub grid: Option<(ExactNum, ExactNum)>,
}

fn all(e: impl Fn() -> Expect) -> Bracket {
    Bracket {
        want: MODES.iter().map(|_| e()).collect(),
        grid: None,
    }
}

/// The expected result of `f` at `p` in each of [`MODES`], from `f` at `PH` bits in Down and Up (see the module doc).
pub fn bracket_expect(
    f: &dyn Fn(usize, RoundingMode, &mut Consts) -> ExactNum,
    p: usize,
    cc: &mut Consts,
) -> Result<Bracket, E> {
    let ph = p.div_ceil(64) * 64 + 128;
    let (lo, hi) = (f(ph, RoundingMode::Down, cc), f(ph, RoundingMode::Up, cc));
    // A finite largest value at `PH` on one side and an infinity of the same sign on the other: an overflow.
    if !lo.is_nan()
        && !hi.is_nan()
        && lo.is_inf() != hi.is_inf()
        && lo.is_negative() == hi.is_negative()
    {
        let neg = lo.is_negative();
        return Ok(all(|| Expect::Overflow(neg)));
    }
    if lo.is_nan() || hi.is_nan() || lo.is_inf() || hi.is_inf() {
        if lo.is_nan() && hi.is_nan() {
            return Ok(all(|| Expect::Special(None)));
        }
        if lo.is_inf() && hi.is_inf() && lo.is_negative() == hi.is_negative() {
            let neg = lo.is_negative();
            return Ok(all(|| Expect::Special(Some(neg))));
        }
        return Err(format!("bracket at {ph} bits is not finite: Down {lo:?}, Up {hi:?}").into());
    }
    if lo.compare(&hi) == Some(0) {
        if lo.inexact() || hi.inexact() {
            return Err(format!("bracket at {ph} bits is one value but inexact: {lo:?}").into());
        }
        let Some(v) = to_mpfr(&lo, cc)? else {
            // an exact subnormal: rounded at `p` on that precision's subnormal grid (`mismatch`)
            return Ok(Bracket {
                grid: Some((lo, hi)),
                ..all(|| Expect::Grid)
            });
        };
        let mut out = Vec::new();
        for rm in MODES {
            let (w, ix) = mpfr_reference_with(u32::try_from(p)?, rm, set_from(&v))?;
            out.push(Expect::Strict(w, ix));
        }
        return Ok(Bracket {
            want: out,
            grid: None,
        });
    }
    if lo.compare(&hi) != Some(-1) {
        return Err(format!("bracket at {ph} bits is reversed: Down {lo:?}, Up {hi:?}").into());
    }
    match (to_mpfr(&lo, cc)?, to_mpfr(&hi, cc)?) {
        (Some(a), Some(b)) => {
            let mut next = Float::with_val(u32::try_from(ph)?, &a);
            next.next_up();
            if next != b {
                return Err(
                    format!("bracket at {ph} bits is not adjacent: Down {a}, Up {b}").into(),
                );
            }
            let mut mid = Float::with_val(u32::try_from(ph + 1)?, &a + &b);
            mid /= 2;
            let mut out = Vec::new();
            for rm in MODES {
                let (w, ix) = mpfr_reference_with(u32::try_from(p)?, rm, set_from(&mid))?;
                if !ix {
                    return Err(format!("the bracket midpoint {mid} is exact at {p} bits").into());
                }
                out.push(Expect::Strict(w, true));
            }
            Ok(Bracket {
                want: out,
                grid: None,
            })
        }
        _ => Ok(Bracket {
            grid: Some((lo, hi)),
            ..all(|| Expect::Grid)
        }),
    }
}

/// `got` against `want` for mode `rm`; `None` if it matches, else what differs.
fn mismatch(
    got: &ExactNum,
    want: &Expect,
    grid: Option<&(ExactNum, ExactNum)>,
    p: usize,
    rm: RoundingMode,
    cc: &mut Consts,
) -> Result<Option<String>, E> {
    Ok(match want {
        Expect::Special(None) => (!got.is_nan()).then(|| format!("want NaN, got {got:?}")),
        Expect::Special(Some(neg)) => (!(got.is_inf() && got.is_negative() == *neg))
            .then(|| format!("want Inf neg={neg}, got {got:?}")),
        Expect::Strict(w, ix) => {
            if got.is_nan() || got.is_inf() {
                Some(format!("want {w}, got {got:?}"))
            } else if got.is_zero() {
                let ok = w.is_zero()
                    && w.is_sign_negative() == got.is_negative()
                    && got.inexact() == *ix;
                (!ok).then(|| {
                    format!(
                        "want {w} (inexact {ix}), got {got:?} (inexact {})",
                        got.inexact()
                    )
                })
            } else {
                match to_mpfr(got, cc)? {
                    Some(g) if g == *w && got.inexact() == *ix => None,
                    g => Some(format!(
                        "want {w} (inexact {ix}), got {g:?} (inexact {})",
                        got.inexact()
                    )),
                }
            }
        }
        Expect::Overflow(neg) => {
            let toward_zero = matches!(rm, RoundingMode::ToZero | RoundingMode::None)
                || (rm == RoundingMode::Down && !*neg)
                || (rm == RoundingMode::Up && *neg);
            let ok = if toward_zero {
                let mut max = ExactNum::max_value(p);
                if rm != RoundingMode::None {
                    max.set_precision(p, RoundingMode::ToZero)?;
                }
                !got.is_inf()
                    && got.is_negative() == *neg
                    && got.abs_cmp(&max) == Some(0)
                    && got.inexact()
            } else {
                got.is_inf() && got.is_negative() == *neg
            };
            (!ok).then(|| format!("want the overflow (negative {neg}) of {rm:?}, got {got:?}"))
        }
        Expect::Grid => {
            // The bracket is adjacent on the `PH`-bit grid (or one exact value). The subnormal grid at `p`,
            // 2^(EXPONENT_MIN - p), is coarser, and neither its points nor its midpoints lie inside the bracket, so
            // the value rounds at `p` as the bracket's midpoint does: one `ldexp` rounding, subnormal grid included.
            let (lo, hi) = grid.ok_or("Grid without a bracket")?;
            let want = lo.add_full_work(hi).ldexp(-1, p, rm);
            let ok = got.compare(&want) == Some(0)
                && got.is_negative() == want.is_negative()
                && got.inexact() == want.inexact();
            (!ok).then(|| format!("want {want:?} (bracket {lo:?} {hi:?}), got {got:?}"))
        }
    })
}

/// Checks `f` at `p` in all seven modes against [`bracket_expect`]; failures go to `fails`.
fn check_all(
    name: &str,
    f: &dyn Fn(usize, RoundingMode, &mut Consts) -> ExactNum,
    p: usize,
    cc: &mut Consts,
    fails: &mut Vec<String>,
) -> Res {
    let b = match bracket_expect(f, p, cc) {
        Ok(w) => w,
        Err(e) => {
            fails.push(format!("{name} p={p}: {e}"));
            return Ok(());
        }
    };
    for (rm, w) in MODES.into_iter().zip(&b.want) {
        let got = f(p, rm, cc);
        if let Some(m) = mismatch(&got, w, b.grid.as_ref(), p, rm, cc)? {
            fails.push(format!("{name} p={p} {rm:?}: {m}"));
        }
    }
    Ok(())
}

fn report(what: &str, n: usize, fails: &[String]) -> Res {
    for f in fails {
        println!("FAIL {what}: {f}");
    }
    assert!(
        fails.is_empty(),
        "{what}: {} of {n} checks failed (printed above)",
        fails.len()
    );
    Ok(())
}

/// The ext.rs digest operands `operand(i)`, for the indices these families use.
fn operand(i: usize) -> ExactNum {
    let third = |p| ExactNum::from_u8(1, p).div(&ExactNum::from_u8(3, p), p, RoundingMode::ToEven);
    let neg = |mut x: ExactNum| {
        x.inv_sign();
        x
    };
    match i {
        0 => ExactNum::new(64),
        1 => neg(ExactNum::new(64)),
        2 => ExactNum::from_u8(1, 64),
        3 => neg(ExactNum::from_u8(7, 64)),
        4 => third(64),
        5 => neg(third(192)),
        6 => third(320),
        7 => INF_POS,
        8 => INF_NEG,
        9 => ExactNum::nan(Some(Error::InvalidArgument)),
        10 => ExactNum::nan(None),
        21 => ExactNum::from_u8(3, 64)
            .ldexp(EXPONENT_MIN + 1, 64, RoundingMode::None)
            .ldexp(-20, 64, RoundingMode::None),
        _ => neg(third(192)
            .ldexp(EXPONENT_MIN + 1, 192, RoundingMode::None)
            .ldexp(-100, 192, RoundingMode::None)),
    }
}

/// `bracket_expect` agrees with MPFR where MPFR has the function: exact results, an exact tie at `p`, inexact
/// results in both signs, and a non-word `p`.
#[test]
fn bracket_expect_matches_mpfr() -> Res {
    use_full_exponent_range();
    let mut cc = Consts::new()?;
    let two_ties = |p: usize| ExactNum::from_u8(1, 64).ldexp(-(p as i64), 256, RoundingMode::None);
    let mut n = 0;
    for p in [64usize, 65, 100, 128] {
        // (name, zenith f, MPFR f)
        let x3 = ExactNum::from_u8(3, 64);
        let x9 = ExactNum::from_u8(9, 64);
        let xm = operand(5);
        let one = ExactNum::from_u8(1, 64);
        let tie = two_ties(p);
        type Zf<'a> = Box<dyn Fn(usize, RoundingMode, &mut Consts) -> ExactNum + 'a>;
        type Mp = Box<dyn Fn(&mut Float, mpfr::Rnd) -> i32>;
        let f3 = Float::with_val(64, 3);
        let f9 = Float::with_val(64, 9);
        let fm = to_mpfr(&xm, &mut cc)?.ok_or("operand 5")?;
        let fm2 = Float::with_val(fm.prec(), &fm);
        let mut ftie = Float::with_val(512, 1);
        ftie += Float::with_val(64, Float::i_exp(1, -(p as i32)));
        let cases: Vec<(&str, Zf, Mp)> = vec![
            (
                "sqrt(9) exact",
                Box::new(|p, rm, _| x9.sqrt(p, rm)),
                Box::new(move |r, d| mpfr::sqrt(r.as_raw_mut(), f9.as_raw(), d)),
            ),
            (
                "sqrt(3)",
                Box::new(|p, rm, _| x3.sqrt(p, rm)),
                Box::new(move |r, d| mpfr::sqrt(r.as_raw_mut(), f3.as_raw(), d)),
            ),
            (
                "exp(-1/3)",
                Box::new(|p, rm, cc| xm.exp(p, rm, cc)),
                Box::new(move |r, d| mpfr::exp(r.as_raw_mut(), fm.as_raw(), d)),
            ),
            (
                "ai(-1/3)",
                Box::new(|p, rm, cc| xm.ai(p, rm, cc)),
                Box::new(move |r, d| mpfr::ai(r.as_raw_mut(), fm2.as_raw(), d)),
            ),
            (
                "1 + 2^-p exact tie",
                Box::new(|p, rm, _| ExactNum::add_ref(&one, &tie, p, rm)),
                Box::new(move |r, d| mpfr::set(r.as_raw_mut(), ftie.as_raw(), d)),
            ),
        ];
        for (name, zf, mp) in &cases {
            let b = bracket_expect(zf.as_ref(), p, &mut cc)?;
            for (rm, w) in MODES.into_iter().zip(&b.want) {
                let (d, dix) = mpfr_reference_with(u32::try_from(p)?, rm, mp.as_ref())?;
                match w {
                    Expect::Strict(b, bix) => assert!(
                        *b == d && *bix == dix,
                        "{name} p={p} {rm:?}: {b} {bix} vs MPFR {d} {dix}"
                    ),
                    other => return Err(format!("{name} p={p} {rm:?}: {other:?}").into()),
                }
                n += 1;
            }
        }
    }
    assert_eq!(n, 4 * 5 * 7);
    Ok(())
}

/// Spence item 2: the thirteen Jacobi functions on the digest's `u` and `m`, at 64, 100
/// and 128 bits, in all seven modes.
#[test]
fn jacobi_all_seven_modes() -> Res {
    use_full_exponent_range();
    type J = fn(&ExactNum, &ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum;
    let fns: [(&str, J); 13] = [
        ("am", ExactNum::jacobi_am),
        ("sn", ExactNum::jacobi_sn),
        ("cn", ExactNum::jacobi_cn),
        ("dn", ExactNum::jacobi_dn),
        ("cd", ExactNum::jacobi_cd),
        ("ns", ExactNum::jacobi_ns),
        ("nc", ExactNum::jacobi_nc),
        ("nd", ExactNum::jacobi_nd),
        ("sc", ExactNum::jacobi_sc),
        ("sd", ExactNum::jacobi_sd),
        ("cs", ExactNum::jacobi_cs),
        ("ds", ExactNum::jacobi_ds),
        ("dc", ExactNum::jacobi_dc),
    ];
    let ms = [
        ExactNum::new(64),
        ExactNum::from_u8(1, 64),
        ExactNum::from_u8(1, 64).ldexp(-1, 64, RoundingMode::None),
        ExactNum::from_u8(3, 64).ldexp(-3, 64, RoundingMode::None),
        ExactNum::from_u8(3, 64),
    ];
    let mut cc = Consts::new()?;
    let (mut fails, mut n) = (Vec::new(), 0usize);
    for (name, f) in fns {
        // Operand 21 is the digest's tiny u = 3 * 2^(EXPONENT_MIN - 19) (F-21/F-22, fixed in b17: the leading term
        // of the series, rounded by the side of the next term); ns, cs and ds overflow there (+-MAX toward zero, F-B).
        for i in [0usize, 1, 2, 3, 5, 7, 21] {
            let u = operand(i);
            for (mi, m) in ms.iter().enumerate() {
                for p in [64usize, 100, 128] {
                    let g = |p, rm, cc: &mut Consts| f(&u, m, p, rm, cc);
                    check_all(&format!("{name}(u{i}|m{mi})"), &g, p, &mut cc, &mut fails)?;
                    n += 7;
                }
            }
        }
    }
    assert_eq!(n, 13 * 7 * 5 * 3 * 7);
    report("jacobi", n, &fails)
}

/// Spence item 3: Ai, Bi, Ai', Bi' on the digest's operands, at 64, 100 and 128 bits, in all seven modes; Ai also
/// strictly against MPFR's `mpfr_ai`.
#[test]
fn airy_all_seven_modes() -> Res {
    use_full_exponent_range();
    type A = fn(&ExactNum, usize, RoundingMode, &mut Consts) -> ExactNum;
    let fns: [(&str, A); 4] = [
        ("ai", ExactNum::ai),
        ("bi", ExactNum::bi),
        ("ai_prime", ExactNum::ai_prime),
        ("bi_prime", ExactNum::bi_prime),
    ];
    let neg = |mut x: ExactNum| {
        x.inv_sign();
        x
    };
    let mut xs: Vec<ExactNum> = (0..11).chain(21..23).map(operand).collect();
    xs.push(ExactNum::from_u8(3, 64).ldexp(-2, 64, RoundingMode::None));
    xs.push(neg(ExactNum::from_u8(5, 64)));
    xs.push(ExactNum::from_u8(30, 64));
    xs.push(neg(ExactNum::from_u8(30, 64)));
    xs.push(ExactNum::from_u8(9, 64).ldexp(-1, 64, RoundingMode::None));
    let mut cc = Consts::new()?;
    let (mut fails, mut n) = (Vec::new(), 0usize);
    for (k, x) in xs.iter().enumerate() {
        for p in [64usize, 100, 128] {
            for (name, f) in fns {
                let g = |p, rm, cc: &mut Consts| f(x, p, rm, cc);
                check_all(&format!("{name}(x{k})"), &g, p, &mut cc, &mut fails)?;
                n += 7;
            }
            // Ai against MPFR, where the argument is finite and in MPFR's range
            if x.is_nan() || x.is_inf() || x.is_subnormal() {
                continue;
            }
            let fx = to_mpfr(x, &mut cc)?.ok_or("argument")?;
            for rm in MODES {
                let (w, ix) = mpfr_reference_with(u32::try_from(p)?, rm, |r, d| {
                    mpfr::ai(r.as_raw_mut(), fx.as_raw(), d)
                })?;
                let got = x.ai(p, rm, &mut cc);
                if let Some(m) = mismatch(&got, &Expect::Strict(w, ix), None, p, rm, &mut cc)? {
                    fails.push(format!("ai(x{k}) vs mpfr_ai p={p} {rm:?}: {m}"));
                }
                n += 1;
            }
        }
    }
    report("airy", n, &fails)
}
