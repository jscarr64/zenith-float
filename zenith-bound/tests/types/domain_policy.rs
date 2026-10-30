//! CONTRACT 10.1, 15.1, O6; oracle: the clauses (two variants, Strict is Default)
//! COVERS: DomainPolicy::Strict, DomainPolicy::Clip

use zenith_bound::policy::DomainPolicy;

// DomainPolicy is #[non_exhaustive] (options unit), so these tests compare instead of matching exhaustively.
fn name(p: &DomainPolicy) -> &'static str {
    if *p == DomainPolicy::Strict {
        "Strict"
    } else if *p == DomainPolicy::Clip {
        "Clip"
    } else {
        "unknown"
    }
}

#[test]
fn default_is_strict() -> Result<(), String> {
    if DomainPolicy::default() == DomainPolicy::Strict {
        Ok(())
    } else {
        Err(format!("default was {:?}", DomainPolicy::default()))
    }
}

#[test]
fn variants_are_distinct_and_named() -> Result<(), String> {
    let all = [DomainPolicy::Strict, DomainPolicy::Clip];
    for (i, a) in all.iter().enumerate() {
        if format!("{:?}", a) != name(a) {
            return Err(format!("Debug of {} was {:?}", name(a), a));
        }
        for (j, b) in all.iter().enumerate() {
            if (a == b) != (i == j) {
                return Err(format!("{:?} == {:?} gave {}", a, b, a == b));
            }
        }
    }
    if DomainPolicy::default() != DomainPolicy::Strict {
        return Err("default != Strict".into());
    }
    Ok(())
}
