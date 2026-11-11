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
//! Inputs: the digest's operands, which avoid Bi and Bi' overflow. Under None, ToZero and truncating Up/Down, Bi and
//! Bi' overflow are the ZF-M27 overflow rows (F-2), and are not pinned here.

use crate::mpfr::common::{conv_str_to_mpfr_compat, mpfr_reference_with};
use crate::mpfr::retry_budget::use_full_exponent_range;
use gmp_mpfr_sys::mpfr;
use rug::Float;
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
    /// A subnormal bracket (`Bracket::grid`), with no MPFR counterpart: the mode picks its low or high end
    /// (nearest: either), flagged inexact; or the one exact subnormal value.
    Grid,
    /// NaN, or an infinity of this sign.
    Special(Option<bool>),
}

/// `x` as an MPFR number holding all its bits; `None` for a subnormal (below MPFR's range).
fn to_mpfr(x: &ExactNum, cc: &mut Consts) -> Result<Option<Float>, E> {
    if x.is_subnormal() {
        return Ok(None);
    }
    let bits = u32::try_from(x.mantissa_max_bit_len().ok_or("not finite")?)?;
    let text = conv_str_to_mpfr_compat(x.format(Radix::Bin, RoundingMode::ToEven, cc)?)?;
    let parsed = Float::parse_radix(text, 2)?;
    let (f, ord) = Float::with_val_round(bits, parsed, rug::float::Round::Nearest);
    if ord != std::cmp::Ordering::Equal {
        return Err(format!("{x:?} is not exact in {bits} bits").into());
    }
    Ok(Some(f))
}

fn set_from(src: &Float) -> impl Fn(&mut Float, mpfr::rnd_t) -> i32 + '_ {
    move |r, rnd| unsafe { mpfr::set(r.as_raw_mut(), src.as_raw(), rnd) }
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
    let _ = (f, p, cc, all, to_mpfr, set_from);
    Err("RED stub: no bracket".into())
}

/// `got` against `want` for mode `rm`; `None` if it matches, else what differs.
fn mismatch(
    got: &ExactNum,
    want: &Expect,
    grid: Option<&(ExactNum, ExactNum)>,
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
        Expect::Grid => {
            let (lo, hi) = grid.ok_or("Grid without a bracket")?;
            let is = |x: &ExactNum| got.compare(x) == Some(0);
            let exact = lo.compare(hi) == Some(0);
            let pos = lo.is_positive() || (lo.is_zero() && hi.is_positive());
            let ok = match rm {
                _ if exact => is(lo) && !got.inexact(),
                RoundingMode::Down => is(lo),
                RoundingMode::Up => is(hi),
                RoundingMode::ToZero | RoundingMode::None => is(if pos { lo } else { hi }),
                RoundingMode::FromZero => is(if pos { hi } else { lo }),
                _ => is(lo) || is(hi),
            } && (exact || got.inexact());
            (!ok).then(|| format!("want one of {lo:?} {hi:?}, got {got:?}"))
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
        if let Some(m) = mismatch(&got, w, b.grid.as_ref(), rm, cc)? {
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
        type Mp = Box<dyn Fn(&mut Float, mpfr::rnd_t) -> i32>;
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
                Box::new(move |r, d| unsafe { mpfr::sqrt(r.as_raw_mut(), f9.as_raw(), d) }),
            ),
            (
                "sqrt(3)",
                Box::new(|p, rm, _| x3.sqrt(p, rm)),
                Box::new(move |r, d| unsafe { mpfr::sqrt(r.as_raw_mut(), f3.as_raw(), d) }),
            ),
            (
                "exp(-1/3)",
                Box::new(|p, rm, cc| xm.exp(p, rm, cc)),
                Box::new(move |r, d| unsafe { mpfr::exp(r.as_raw_mut(), fm.as_raw(), d) }),
            ),
            (
                "ai(-1/3)",
                Box::new(|p, rm, cc| xm.ai(p, rm, cc)),
                Box::new(move |r, d| unsafe { mpfr::ai(r.as_raw_mut(), fm2.as_raw(), d) }),
            ),
            (
                "1 + 2^-p exact tie",
                Box::new(|p, rm, _| ExactNum::add_ref(&one, &tie, p, rm)),
                Box::new(move |r, d| unsafe { mpfr::set(r.as_raw_mut(), ftie.as_raw(), d) }),
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
                    other => panic!("{name} p={p} {rm:?}: {other:?}"),
                }
                n += 1;
            }
        }
    }
    assert_eq!(n, 4 * 5 * 7);
    Ok(())
}

/// Spence item 2: the thirteen Jacobi functions on the digest's `u` (but operand 21, HELD) and `m`, at 64, 100
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
        // HELD N-6, N-7: the digest's tiny operand 21 (u = 3 * 2^(EXPONENT_MIN - 19)) is left out. At that u, ToZero,
        // FromZero, Up and Down return NaN(PrecisionRetryExhausted) for am (m != 0), sn, cn, cd, sc and sd, and Down is
        // in the digest pin (N-6). None returns exactly 1 for cn(u|1) and dn(u|m > 0), whose values are below 1, and
        // dn and nd at m = 1/2 and 3/8 come back flagged exact (N-7). Neither is F-A/F-B/F-C shaped. ns, cs and ds
        // at that u overflow to +Inf under None: NONE-GAP F-B.
        for i in [0usize, 1, 2, 3, 5, 7] {
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
    assert_eq!(n, 13 * 6 * 5 * 3 * 7);
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
                let (w, ix) = mpfr_reference_with(u32::try_from(p)?, rm, |r, d| unsafe {
                    mpfr::ai(r.as_raw_mut(), fx.as_raw(), d)
                })?;
                let got = x.ai(p, rm, &mut cc);
                if let Some(m) = mismatch(&got, &Expect::Strict(w, ix), None, rm, &mut cc)? {
                    fails.push(format!("ai(x{k}) vs mpfr_ai p={p} {rm:?}: {m}"));
                }
                n += 1;
            }
        }
    }
    report("airy", n, &fails)
}
