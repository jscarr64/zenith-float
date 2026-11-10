//! CONTRACT 6.4, 15.3, O11-W; PLAN C5 rulings a/b (2026-10-23): struct variants {op, detail},
//! #[non_exhaustive], NotPositiveDefinite and Singular; NotConverged (Paul 2026-10-23: an iteration out of budget, e.g.
//! eig past max_sweeps, returns NotConverged, never a silently loose answer); DimensionMismatch (S0, 2026-10-24: operand shapes
//! disagree, e.g. matmul, solve, GLS; reserved up front so no later slice edits the enum); NotSymmetric (LA plan, 2026-10-24: LA5a).
//! COVERS: BoundError::Representation, BoundError::InvalidArgument, BoundError::Domain, BoundError::Resource, BoundError::NotPositiveDefinite, BoundError::Singular, BoundError::NotConverged, BoundError::DimensionMismatch, BoundError::NotSymmetric, BoundError::op, BoundError::detail, Display for BoundError, Error for BoundError

use zenith_bound::ball::Ball;
use zenith_bound::error::{BoundError, BoundResult};
use zenith_bound::mag::Mag;
use zenith_float_num::ExactNum;

fn all() -> [BoundError; 9] {
    [
        BoundError::Representation {
            op: "op-rep",
            detail: "d-rep",
        },
        BoundError::InvalidArgument {
            op: "op-inv",
            detail: "d-inv",
        },
        BoundError::Domain {
            op: "op-dom",
            detail: "d-dom",
        },
        BoundError::Resource {
            op: "op-res",
            detail: "d-res",
        },
        BoundError::NotPositiveDefinite {
            op: "op-npd",
            detail: "d-npd",
        },
        BoundError::Singular {
            op: "op-sing",
            detail: "d-sing",
        },
        BoundError::NotConverged {
            op: "op-conv",
            detail: "d-conv",
        },
        BoundError::DimensionMismatch {
            op: "op-dim",
            detail: "d-dim",
        },
        BoundError::NotSymmetric {
            op: "op-sym",
            detail: "d-sym",
        },
    ]
}

fn name(e: &BoundError) -> &'static str {
    match e {
        BoundError::Representation { .. } => "Representation",
        BoundError::InvalidArgument { .. } => "InvalidArgument",
        BoundError::Domain { .. } => "Domain",
        BoundError::Resource { .. } => "Resource",
        BoundError::NotPositiveDefinite { .. } => "NotPositiveDefinite",
        BoundError::Singular { .. } => "Singular",
        BoundError::NotConverged { .. } => "NotConverged",
        BoundError::DimensionMismatch { .. } => "DimensionMismatch",
        BoundError::NotSymmetric { .. } => "NotSymmetric",
        _ => "unknown",
    }
}

#[test]
fn every_variant_equals_only_itself() -> Result<(), String> {
    for (i, a) in all().iter().enumerate() {
        for (j, b) in all().iter().enumerate() {
            if (a == b) != (i == j) {
                return Err(format!("{:?} == {:?} gave {}", a, b, a == b));
            }
        }
    }
    Ok(())
}

#[test]
fn op_and_detail_return_the_fields() -> Result<(), String> {
    for e in all() {
        if !e.op().starts_with("op-") || !e.detail().starts_with("d-") {
            return Err(format!(
                "{}: op {:?} detail {:?}",
                name(&e),
                e.op(),
                e.detail()
            ));
        }
        if e.op()[3..] != e.detail()[2..] {
            return Err(format!("{}: op/detail from different variants", name(&e)));
        }
    }
    Ok(())
}

#[test]
fn debug_names_the_variant() -> Result<(), String> {
    for e in all() {
        if !format!("{:?}", e).starts_with(name(&e)) {
            return Err(format!("Debug of {} was {:?}", name(&e), e));
        }
    }
    Ok(())
}

#[test]
fn display_names_kind_op_and_detail() -> Result<(), String> {
    let shown: Vec<String> = all().iter().map(|e| format!("{}", e)).collect();
    for (i, (a, e)) in shown.iter().zip(all().iter()).enumerate() {
        if !a.contains(e.op()) || !a.contains(e.detail()) {
            return Err(format!("Display of {} lacks op/detail: {:?}", name(e), a));
        }
        for (j, b) in shown.iter().enumerate() {
            if j > i && a == b {
                return Err(format!("duplicate Display {:?}", a));
            }
        }
    }
    Ok(())
}

#[test]
fn bound_error_is_a_core_error() -> Result<(), String> {
    let e = BoundError::Domain {
        op: "x",
        detail: "y",
    };
    let shown = format!("{}", e);
    let boxed: Box<dyn core::error::Error> = Box::new(e);
    if boxed.source().is_some() {
        return Err("BoundError has no source".into());
    }
    if boxed.to_string() != shown {
        return Err("Error::to_string differs from Display".into());
    }
    Ok(())
}

#[test]
fn bound_result_is_result_of_bound_error() -> Result<(), String> {
    let err: BoundResult<u32> = Err(BoundError::Resource {
        op: "a",
        detail: "b",
    });
    match err {
        Err(BoundError::Resource {
            op: "a",
            detail: "b",
        }) => Ok(()),
        other => Err(format!("got {:?}", other)),
    }
}

#[test]
fn callers_name_their_op() -> Result<(), String> {
    let cases: [(&str, BoundResult<()>, &str); 4] = [
        (
            "Mag::finite(0, 1)",
            Mag::finite(0, 1).map(|_| ()),
            "Mag::finite",
        ),
        (
            "Mag::finite(1, MIN)",
            Mag::finite(1, i64::MIN).map(|_| ()),
            "Mag::finite",
        ),
        (
            "Mag::infinity().to_exact_num()",
            Mag::infinity().to_exact_num().map(|_| ()),
            "Mag::to_exact_num",
        ),
        (
            "Ball::exact(1, 0)",
            Ball::exact(ExactNum::from(1i64), 0).map(|_| ()),
            "Ball::exact",
        ),
    ];
    for (what, got, op) in cases {
        match got {
            Err(e) if e.op() == op && !e.detail().is_empty() => {}
            other => {
                return Err(format!(
                    "{}: want op {:?} with a detail, got {:?}",
                    what, op, other
                ))
            }
        }
    }
    Ok(())
}

#[test]
fn each_variant_matches_its_own_pattern() -> Result<(), String> {
    let checks = [
        matches!(all()[0], BoundError::Representation { op: "op-rep", .. }),
        matches!(all()[1], BoundError::InvalidArgument { op: "op-inv", .. }),
        matches!(all()[2], BoundError::Domain { op: "op-dom", .. }),
        matches!(all()[3], BoundError::Resource { op: "op-res", .. }),
        matches!(
            all()[4],
            BoundError::NotPositiveDefinite { op: "op-npd", .. }
        ),
        matches!(all()[5], BoundError::Singular { op: "op-sing", .. }),
        matches!(all()[6], BoundError::NotConverged { op: "op-conv", .. }),
        matches!(all()[7], BoundError::DimensionMismatch { op: "op-dim", .. }),
        matches!(all()[8], BoundError::NotSymmetric { op: "op-sym", .. }),
    ];
    match checks.iter().position(|ok| !ok) {
        None => Ok(()),
        Some(i) => Err(format!("variant {} did not match its pattern", i)),
    }
}

#[test]
fn display_kind_text_per_variant() -> Result<(), String> {
    // L0b mutant 104-E-kind-npd: each variant's Display leads with its own kind text.
    let want = [
        "value not representable",
        "invalid argument",
        "input outside the domain",
        "resource limit reached",
        "matrix not positive definite",
        "matrix singular",
        "iteration did not converge",
        "dimensions do not match",
        "matrix not symmetric",
    ];
    for (e, w) in all().iter().zip(want.iter()) {
        if !format!("{}", e).starts_with(w) {
            return Err(format!(
                "Display of {} must start with {:?}: {}",
                name(e),
                w,
                e
            ));
        }
    }
    Ok(())
}

#[test]
fn dimension_mismatch_constructs_and_displays() -> Result<(), String> {
    // S0: the shape variant every array and linear-algebra slice shares.
    let e = BoundError::DimensionMismatch {
        op: "BallMat::matmul",
        detail: "inner dimensions 2 and 3",
    };
    let want = "dimensions do not match in BallMat::matmul: inner dimensions 2 and 3";
    if e.to_string() != want {
        return Err(format!("Display was {:?}, want {:?}", e.to_string(), want));
    }
    if e.op() != "BallMat::matmul" || e.detail() != "inner dimensions 2 and 3" {
        return Err(format!("op/detail wrong: {:?}", e));
    }
    Ok(())
}

#[test]
fn not_symmetric_constructs_and_displays() -> Result<(), String> {
    // LA plan (2026-10-24): LA5a's refusal for a non-symmetric input (CONTRACT 16.LA.1).
    let e = BoundError::NotSymmetric {
        op: "eig_sym",
        detail: "a[0][1] != a[1][0]",
    };
    let want = "matrix not symmetric in eig_sym: a[0][1] != a[1][0]";
    if e.to_string() != want {
        return Err(format!("Display was {:?}, want {:?}", e.to_string(), want));
    }
    if e.op() != "eig_sym" || e.detail() != "a[0][1] != a[1][0]" {
        return Err(format!("op/detail wrong: {:?}", e));
    }
    Ok(())
}
