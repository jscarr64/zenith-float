//! API-freeze unit 5 (Paul/Rusty 2026-10-23), ZB-TD02 (Eugene 2026-10-30): every public prec is `Precision` (= zf `usize`
//! alias), validated at the constructors against `MIN_PRECISION..=MAX_PRECISION`. P1/P6: `MAX_PRECISION` is `2^32 - 1` (= `u32::MAX`) on every
//! 64-bit target, read from zf `zenith_float_num::MAX_PRECISION` (ZF-M32, decided 2026-10-30 6:37 PM CT). One expectation per case on every target: no cfg and no
//! `usize::BITS` branch. Edge list: ~/eugene-work/zfzb/edge-lists/ZB-TD02.md.
//! P6 (64-bit only): host reds only. RED until the ZF-M32 code lands (host has 2^60 - 1): rows A, B, J, M, N; rows Z1-Z3 (byte budget) are RED until the ZF-M32 size-check API lands.
//! COVERS: Ball::exact, Ball::all_reals, Ball::from_mid_rad, BallArray::zeros, Ball::add, Ball::mul
//! (Precision, MIN_PRECISION and MAX_PRECISION are a type alias and consts, which pubfn_coverage does not list.)

use zenith_bound::array::BallArray;
use zenith_bound::ball::Ball;
use zenith_bound::error::BoundError;
use zenith_bound::mag::Mag;
use zenith_bound::{Precision, MAX_PRECISION, MIN_PRECISION};
use zenith_float_num::ExactNum;

/// `2^32 - 1` (= `u32::MAX`): the one precision limit on every 64-bit target (P6, ZF-M32).
const WANT_MAX: u64 = (1 << 32) - 1;

fn invalid<T>(r: Result<T, BoundError>, what: &str) -> Result<(), String> {
    match r {
        Err(BoundError::InvalidArgument { .. }) => Ok(()),
        Err(e) => Err(format!("{}: want InvalidArgument, got {:?}", what, e)),
        Ok(_) => Err(format!("{}: want InvalidArgument, got Ok", what)),
    }
}

fn three() -> ExactNum {
    ExactNum::from(3i64)
}

/// Every precision-taking entry point refuses `p` with InvalidArgument.
fn all_refuse(p: Precision) -> Result<(), String> {
    invalid(Ball::exact(three(), p), &format!("exact {}", p))?;
    invalid(Ball::all_reals(p), &format!("all_reals {}", p))?;
    invalid(
        Ball::from_mid_rad(three(), Mag::zero(), p),
        &format!("from_mid_rad {}", p),
    )?;
    invalid(
        BallArray::zeros(vec![2].into_boxed_slice(), p),
        &format!("zeros {}", p),
    )?;
    let a = Ball::exact(three(), 53).map_err(|e| format!("exact(3, 53): {:?}", e))?;
    let b = Ball::exact(three(), 53).map_err(|e| format!("exact(3, 53): {:?}", e))?;
    invalid(a.add(&b, p), &format!("add at {}", p))?;
    let c = Ball::exact(three(), 53).map_err(|e| format!("exact(3, 53): {:?}", e))?;
    invalid(c.mul(&b, p), &format!("mul at {}", p))
}


#[test]
fn max_precision_is_2_pow_32_minus_1_on_every_target() -> Result<(), String> {
    // Rows A and D.
    if MIN_PRECISION != 1 {
        return Err(format!("MIN_PRECISION {}", MIN_PRECISION));
    }
    let got = u64::try_from(MAX_PRECISION).map_err(|_| String::from("MAX_PRECISION as u64"))?;
    if got != WANT_MAX {
        return Err(format!(
            "zb MAX_PRECISION {} want {} (waiting on ZF-M32)",
            got, WANT_MAX
        ));
    }
    Ok(())
}

#[test]
fn zf_max_precision_is_2_pow_32_minus_1_on_every_target() -> Result<(), String> {
    // Row B (ZF-M32).
    let got = u64::try_from(zenith_float_num::MAX_PRECISION)
        .map_err(|_| String::from("zf MAX_PRECISION as u64"))?;
    if got != WANT_MAX {
        return Err(format!(
            "zf MAX_PRECISION {} want {} (waiting on ZF-M32)",
            got, WANT_MAX
        ));
    }
    Ok(())
}

#[test]
fn zb_max_precision_reads_zf() -> Result<(), String> {
    // Row C.
    if MAX_PRECISION != zenith_float_num::MAX_PRECISION {
        return Err(format!(
            "zb MAX_PRECISION {} != zf MAX_PRECISION {}",
            MAX_PRECISION,
            zenith_float_num::MAX_PRECISION
        ));
    }
    Ok(())
}

#[test]
fn constructors_accept_every_precision_up_to_2_pow_32_minus_1() -> Result<(), String> {
    // Rows F, G, H, I.
    let ps: [u64; 21] = [
        1,
        2,
        32,
        33,
        53,
        63,
        64,
        65,
        128,
        256,
        512,
        1024,
        4096,
        16384,
        65536,
        1 << 28,
        (1 << 30) - 1,
        1 << 30,
        1 << 31,
        WANT_MAX - 1,
        WANT_MAX,
    ];
    for p in ps {
        let p = usize::try_from(p).map_err(|_| format!("{} as usize", p))?;
        let b = Ball::exact(three(), p).map_err(|e| format!("exact {}: {:?}", p, e))?;
        if b.prec() != p {
            return Err(format!("prec {} read back as {}", p, b.prec()));
        }
        let a = Ball::all_reals(p).map_err(|e| format!("all_reals {}: {:?}", p, e))?;
        let r = Ball::from_mid_rad(ExactNum::from(1i64), Mag::zero(), p)
            .map_err(|e| format!("from_mid_rad {}: {:?}", p, e))?;
        if a.prec() != p || r.prec() != p || r.into_parts().2 != p {
            return Err(format!("prec {} lost", p));
        }
        BallArray::zeros(vec![2].into_boxed_slice(), p)
            .map_err(|e| format!("zeros {}: {:?}", p, e))?;
    }
    Ok(())
}

#[test]
fn ball_at_the_limit_is_exact_three() -> Result<(), String> {
    // Row O: same bits on every target.
    let p = usize::try_from(WANT_MAX).map_err(|_| String::from("2^32 - 1 as usize"))?;
    let b = Ball::exact(three(), p).map_err(|e| format!("exact(3, 2^32 - 1): {:?}", e))?;
    let (m, r, q) = b.into_parts();
    if m != three() || !r.is_zero() || q != p {
        return Err(format!("exact(3, 2^32 - 1): got {:?} +- {:?} at {}", m, r, q));
    }
    Ok(())
}

#[test]
fn precisions_above_2_pow_32_minus_1_are_refused_on_every_target() -> Result<(), String> {
    // Rows E, J, L, M, N.
    for p in [
        0u64,
        WANT_MAX + 1,
        WANT_MAX + 2,
        1 << 33,
        1 << 40,
        (1 << 60) - 1,
        1 << 60,
        (1 << 62) + 1,
    ] {
        all_refuse(usize::try_from(p).map_err(|_| format!("{} as usize", p))?)?;
    }
    all_refuse(usize::MAX)
}

/// ZF-M32 byte budget: a single buffer larger than 2^31 bytes is refused with Error::MemoryAllocation
/// before any allocation (zb maps it to BoundError::Resource, arith.rs:591).
const BUDGET_BYTES: u64 = 1 << 31;

/// Rows Z1-Z3 (budget-1, budget, budget+1 bytes): (bytes, refused?).
const BUDGET_ROWS: [(u64, bool); 3] = [
    (BUDGET_BYTES - 1, false),
    (BUDGET_BYTES, false),
    (BUDGET_BYTES + 1, true),
];

#[test]
fn byte_budget_rows_wait_on_zf_m32_api() -> Result<(), String> {
    // RED with a dependency note: ZF-M32 has not landed a public up-front size check, and no public
    // zb or zf entry point takes a byte count, so these rows cannot be driven yet. Through precision
    // the largest buffer is the 4p work buffer at MAX_PRECISION, exactly 2^31 bytes (row Z2), so
    // budget+1 is unreachable without the M32 check API. When it lands, call it for each row:
    // refused rows want Error::MemoryAllocation (BoundError::Resource in zb), the others want Ok.
    let rows: Vec<String> = BUDGET_ROWS
        .iter()
        .map(|(b, refused)| format!("{} bytes refused={}", b, refused))
        .collect();
    Err(format!(
        "waiting on ZF-M32 byte-budget check API; rows not driven: {}",
        rows.join(", ")
    ))
}
