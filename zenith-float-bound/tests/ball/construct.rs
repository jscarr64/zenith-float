//! CONTRACT 6.1-6.5 (D15a, D15e), 8.1 (O17-W), 15.2-15.3; TEST_PLAN L0b-1, L0b-13 (construction paths), L0b-14
//! (degenerate and all-reals). Oracle: the values are built by hand from small integers and powers of two.
//! COVERS: Ball::exact, Ball::from_mid_rad, Ball::all_reals, Ball::mid, Ball::rad, Ball::prec, Ball::is_informative, Ball::is_exact, Ball::into_parts

use zenith_bound::ball::Ball;
use zenith_bound::error::BoundError;
use zenith_bound::mag::Mag;
use zenith_float_num::{ExactNum, INF_NEG, INF_POS, NAN};

const TOP: u64 = 1u64 << 63;

fn num(v: i64) -> ExactNum {
    ExactNum::from(v)
}

fn mag(mant: u64, exp: i64) -> Result<Mag, String> {
    Mag::finite(mant, exp).map_err(|e| format!("Mag::finite({}, {}): {:?}", mant, exp, e))
}

fn ok(what: &str, r: Result<Ball, BoundError>) -> Result<Ball, String> {
    r.map_err(|e| format!("{}: unexpected {:?}", what, e))
}

fn want_err(what: &str, r: Result<Ball, BoundError>, want: BoundError) -> Result<(), String> {
    match r {
        Err(e) if core::mem::discriminant(&e) == core::mem::discriminant(&want) => Ok(()),
        Err(e) => Err(format!("{}: want {:?}, got {:?}", what, want, e)),
        Ok(b) => Err(format!("{}: want {:?}, got Ok({:?})", what, want, b)),
    }
}

fn is_pos_zero(x: &ExactNum) -> bool {
    x.is_zero() && !x.is_negative()
}

#[test]
fn prec_zero_is_invalid_argument() -> Result<(), String> {
    want_err(
        "exact(3, 0)",
        Ball::exact(num(3), 0),
        BoundError::InvalidArgument { op: "", detail: "" },
    )?;
    want_err(
        "from_mid_rad(3, 0, 0)",
        Ball::from_mid_rad(num(3), Mag::zero(), 0),
        BoundError::InvalidArgument { op: "", detail: "" },
    )?;
    want_err(
        "from_mid_rad(3, inf, 0)",
        Ball::from_mid_rad(num(3), Mag::infinity(), 0),
        BoundError::InvalidArgument { op: "", detail: "" },
    )?;
    want_err(
        "all_reals(0)",
        Ball::all_reals(0),
        BoundError::InvalidArgument { op: "", detail: "" },
    )
}

#[test]
fn nan_and_inf_mid_are_representation() -> Result<(), String> {
    for (name, x) in [("NAN", NAN), ("INF_POS", INF_POS), ("INF_NEG", INF_NEG)] {
        want_err(
            name,
            Ball::exact(x, 64),
            BoundError::Representation { op: "", detail: "" },
        )?;
    }
    for (name, x) in [("NAN", NAN), ("INF_POS", INF_POS), ("INF_NEG", INF_NEG)] {
        want_err(
            name,
            Ball::from_mid_rad(x, mag(TOP, 1)?, 64),
            BoundError::Representation { op: "", detail: "" },
        )?;
    }
    // A NaN midpoint stays an error even with an Infinity radius (6.2: never stored, not silently canonicalised).
    want_err(
        "NAN ± inf",
        Ball::from_mid_rad(NAN, Mag::infinity(), 64),
        BoundError::Representation { op: "", detail: "" },
    )
}

#[test]
fn exact_is_a_point_with_zero_radius() -> Result<(), String> {
    for v in [3i64, -3, 1, -1, i64::MAX, i64::MIN] {
        let b = ok("exact", Ball::exact(num(v), 64))?;
        if b.mid() != &num(v) || !b.rad().is_zero() || b.prec() != 64 {
            return Err(format!("exact({}): got {:?}", v, b));
        }
        if !b.is_exact() || !b.is_informative() {
            return Err(format!("exact({}): is_exact/is_informative wrong", v));
        }
    }
    Ok(())
}

#[test]
fn exact_wider_than_prec_is_not_truncated() -> Result<(), String> {
    // D15a: u64::MAX needs 64 bits; at prec 1 the constructor keeps it exactly with rad = 0.
    let b = ok(
        "exact(u64::MAX, 1)",
        Ball::exact(ExactNum::from(u64::MAX), 1),
    )?;
    if b.mid() != &ExactNum::from(u64::MAX) || !b.rad().is_zero() || b.prec() != 1 {
        return Err(format!("exact(u64::MAX, 1): got {:?}", b));
    }
    let b = ok(
        "from_mid_rad(-u64::MAX, 2^-3, 2)",
        Ball::from_mid_rad(-ExactNum::from(u64::MAX), mag(TOP, -2)?, 2),
    )?;
    if b.mid() != &(-ExactNum::from(u64::MAX)) {
        return Err(format!("wide mid changed: {:?}", b));
    }
    Ok(())
}

#[test]
fn negative_zero_mid_becomes_positive_zero() -> Result<(), String> {
    let nz = -ExactNum::new(64);
    if !(nz.is_zero() && nz.is_negative()) {
        return Err("test setup: -ExactNum::new(64) is not -0".into());
    }
    let b = ok("exact(-0)", Ball::exact(nz, 64))?;
    if !is_pos_zero(b.mid()) || !b.rad().is_zero() {
        return Err(format!("exact(-0): got {:?}", b));
    }
    let b = ok(
        "from_mid_rad(-0, 2^-1)",
        Ball::from_mid_rad(-ExactNum::new(64), mag(TOP, 0)?, 64),
    )?;
    if !is_pos_zero(b.mid()) || b.rad().mant() != Some(TOP) || b.rad().exp() != Some(0) {
        return Err(format!("from_mid_rad(-0, 2^-1): got {:?}", b));
    }
    Ok(())
}

#[test]
fn from_mid_rad_keeps_parts() -> Result<(), String> {
    let b = ok(
        "from_mid_rad(-5, 2^-3, 7)",
        Ball::from_mid_rad(num(-5), mag(TOP | 1, -2)?, 7),
    )?;
    if b.mid() != &num(-5)
        || b.rad().mant() != Some(TOP | 1)
        || b.rad().exp() != Some(-2)
        || b.prec() != 7
    {
        return Err(format!("got {:?}", b));
    }
    if b.is_exact() || !b.is_informative() {
        return Err("finite nonzero radius: is_exact/is_informative wrong".into());
    }
    let b = ok(
        "from_mid_rad(5, 0)",
        Ball::from_mid_rad(num(5), Mag::zero(), 64),
    )?;
    if !b.is_exact() {
        return Err("zero radius is not exact".into());
    }
    Ok(())
}

#[test]
fn infinity_radius_is_canonical_all_reals() -> Result<(), String> {
    for v in [5i64, -5, 0] {
        let b = ok(
            "from_mid_rad(v, inf)",
            Ball::from_mid_rad(num(v), Mag::infinity(), 9),
        )?;
        if !is_pos_zero(b.mid()) || !b.rad().is_infinity() || b.prec() != 9 {
            return Err(format!("from_mid_rad({}, inf): got {:?}", v, b));
        }
        if b.is_informative() || b.is_exact() {
            return Err(format!(
                "from_mid_rad({}, inf): is_informative/is_exact wrong",
                v
            ));
        }
    }
    let b = ok("all_reals(1)", Ball::all_reals(1))?;
    if !is_pos_zero(b.mid()) || !b.rad().is_infinity() || b.prec() != 1 || b.is_informative() {
        return Err(format!("all_reals(1): got {:?}", b));
    }
    Ok(())
}

#[test]
fn prec_is_stored_as_requested() -> Result<(), String> {
    for p in [
        1usize,
        2,
        63,
        64,
        65,
        4096,
        0xFFFF_FFFF,
        1 << 32,
        zenith_bound::MAX_PRECISION,
    ] {
        let b = ok("exact(1, p)", Ball::exact(num(1), p))?;
        if b.prec() != p {
            return Err(format!("prec {}: got {}", p, b.prec()));
        }
    }
    Ok(())
}

#[test]
fn into_parts_moves_the_fields_out() -> Result<(), String> {
    let (m, r, p) = ok(
        "from_mid_rad(7, 2^0, 3)",
        Ball::from_mid_rad(num(7), mag(TOP, 1)?, 3),
    )?
    .into_parts();
    if m != num(7) || r.mant() != Some(TOP) || r.exp() != Some(1) || p != 3 {
        return Err(format!("into_parts: {:?} {:?} {}", m, r, p));
    }
    Ok(())
}

#[test]
fn errors_name_the_constructor_called() -> Result<(), String> {
    // L0b mutant 115-B-exactop: `Ball::exact` reports itself, not the `from_mid_rad` it delegates to.
    let cases: [(&str, Result<Ball, BoundError>); 3] = [
        ("Ball::exact", Ball::exact(num(3), 0)),
        (
            "Ball::from_mid_rad",
            Ball::from_mid_rad(num(3), Mag::zero(), 0),
        ),
        ("Ball::all_reals", Ball::all_reals(0)),
    ];
    for (want, r) in cases {
        match r {
            Err(e) if e.op() == want => {}
            other => return Err(format!("{}: want op {}, got {:?}", want, want, other)),
        }
    }
    Ok(())
}

#[test]
fn positive_zero_mid_keeps_its_precision() -> Result<(), String> {
    // L0b mutant 111-B-poszero: only a -0 midpoint is replaced; a +0 midpoint is kept as given (D15a, CONTRACT 6.3).
    let b = ok("exact(+0 at 200 bits)", Ball::exact(ExactNum::new(200), 64))?;
    // A zero's precision() is Some(0) whatever its buffer, so compare the buffer: 256 bits kept, 64 if rebuilt.
    let want = ExactNum::new(200).mantissa_max_bit_len();
    let got = b.mid().mantissa_max_bit_len();
    if !is_pos_zero(b.mid()) || got != want {
        return Err(format!("+0 midpoint was rebuilt: {:?} vs {:?}", got, want));
    }
    Ok(())
}
