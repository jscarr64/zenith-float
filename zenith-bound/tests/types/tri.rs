//! CONTRACT 15.2, O15-W, D14; oracle: the clauses (three variants; comparison logic is a later unit)
//! COVERS: Tri::True, Tri::False, Tri::Unknown

use zenith_bound::tri::Tri;

fn name(t: &Tri) -> &'static str {
    match t {
        Tri::True => "True",
        Tri::False => "False",
        Tri::Unknown => "Unknown",
    }
}

#[test]
fn tri_variants_are_distinct_and_named() -> Result<(), String> {
    let all = [Tri::True, Tri::False, Tri::Unknown];
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
    Ok(())
}

#[test]
fn tri_is_exhaustive_with_exactly_three_values() -> Result<(), String> {
    // Gap F2 (API freeze): Tri is NOT #[non_exhaustive]; this wildcard-free match in an
    // external crate stops compiling if it ever becomes so (checked by hand: adding it gives E0004, log 1555-zb-l0b-f2-red).
    let all = [Tri::True, Tri::False, Tri::Unknown];
    let mut seen = 0u8;
    for t in &all {
        seen |= match t {
            Tri::True => 1,
            Tri::False => 2,
            Tri::Unknown => 4,
        };
    }
    match seen {
        7 => Ok(()),
        other => Err(format!("variant mask {:b}", other)),
    }
}
