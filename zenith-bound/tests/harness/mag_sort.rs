//! TEST_PLAN 8.1 / L0a-2 "`Ord` transitivity": sorting a set of `Mag` values that includes both zero constructors (zero has one canonical representation), the extreme finite
//! values, a duplicate and `Infinity` gives one and the same sequence under every input permutation, and that sequence is the
//! hand-ordered one (CONTRACT 14.5: zero least, `Infinity` greatest, finite by `(exp, mant)`). A non-transitive or non-total
//! `cmp` makes `sort` depend on the input order, so this is the observable form of transitivity. Permutations: every one of a
//! 6-element subset (720), and seeded Fisher-Yates shuffles of the whole set (SplitMix64, seeds 0..SHUFFLES). The mutant
//! 180-C-cmp-zero-eq (zero compares Equal to a finite value) must be killed here.
//! COVERS: Ord for Mag

use crate::common::expect_eq;
use crate::common::generators::{TOP, ZF_EXPONENT_MAX, ZF_EXPONENT_MIN};
use crate::common::rng::SplitMix64;
use zenith_bound::mag::Mag;

/// How one element is built; listed in ascending order, so the list order is the oracle.
enum Spec {
    Zero,
    ZeroFromFinite,
    Fin(u64, i64),
    Inf,
}

const SET: [Spec; 12] = [
    Spec::Zero,
    Spec::ZeroFromFinite,
    Spec::Fin(TOP, i64::MIN),
    Spec::Fin(TOP, ZF_EXPONENT_MIN),
    Spec::Fin(TOP, -1),
    Spec::Fin(u64::MAX, -1),
    Spec::Fin(TOP, 0),
    Spec::Fin(TOP, 0),
    Spec::Fin(TOP + 1, 0),
    Spec::Fin(TOP, ZF_EXPONENT_MAX),
    Spec::Fin(u64::MAX, i64::MAX),
    Spec::Inf,
];

/// Indices into `SET` of the exhaustive subset: zero, a duplicate pair, the largest finite value, `Infinity`.
const SUBSET: [usize; 6] = [0, 4, 6, 7, 10, 11];

const SHUFFLES: u64 = 500;

/// The observable key of a value: 0 zero, 1 finite with its exponent and mantissa, 2 `Infinity`.
type Key = (u8, i64, u64);

fn build(spec: &Spec) -> Result<Mag, String> {
    match spec {
        Spec::Zero => Ok(Mag::zero()),
        Spec::ZeroFromFinite => {
            Mag::finite(0, 0).map_err(|e| format!("Mag::finite(0, 0) failed: {:?}", e))
        }
        Spec::Fin(m, e) => Mag::finite(*m, *e)
            .map_err(|err| format!("Mag::finite({:#x}, {}) failed: {:?}", m, e, err)),
        Spec::Inf => Ok(Mag::infinity()),
    }
}

fn key(m: &Mag) -> Result<Key, String> {
    match (m.is_zero(), m.is_infinity(), m.exp(), m.mant()) {
        (true, false, _, _) => Ok((0, 0, 0)),
        (false, true, _, _) => Ok((2, 0, 0)),
        (false, false, Some(e), Some(x)) => Ok((1, e, x)),
        // The remaining combinations, listed so the compiler checks the match is exhaustive (no catch-all).
        (true, true, e, x) => Err(format!(
            "inconsistent Mag: both zero and Infinity (exp {:?} mant {:?})",
            e, x
        )),
        (false, false, None, x) => Err(format!(
            "inconsistent Mag: finite with no exp (mant {:?})",
            x
        )),
        (false, false, Some(e), None) => {
            Err(format!("inconsistent Mag: finite with no mant (exp {})", e))
        }
    }
}

/// Builds `SET` in the order `order`, sorts it with `Ord`, and returns the keys.
fn sorted_keys(order: &[usize]) -> Result<Vec<Key>, String> {
    let mut v = Vec::new();
    for i in order {
        v.push(build(
            SET.get(*i).ok_or(format!("index {} outside SET", i))?,
        )?);
    }
    v.sort();
    v.iter().map(key).collect()
}

/// The keys of `order` taken in ascending (`SET`) order: the hand oracle.
fn expected_keys(order: &[usize]) -> Result<Vec<Key>, String> {
    let mut ix = Vec::new();
    ix.extend_from_slice(order);
    ix.sort_unstable();
    let mut out = Vec::new();
    for i in &ix {
        out.push(key(&build(
            SET.get(*i).ok_or(format!("index {} outside SET", i))?,
        )?)?);
    }
    Ok(out)
}

/// Every permutation of `items` (Heap's algorithm, iterative), each passed to `f`; returns the count.
fn for_each_permutation(
    items: &mut [usize],
    f: &mut dyn FnMut(&[usize]) -> Result<(), String>,
) -> Result<u64, String> {
    let n = items.len();
    let mut c = vec![0usize; n];
    let mut count = 1u64;
    f(items)?;
    let mut i = 0usize;
    while i < n {
        let ci = match c.get(i) {
            Some(&x) => x,
            None => return Err(String::from("Heap counter outside its range")),
        };
        if ci < i {
            if i.is_multiple_of(2) {
                items.swap(0, i);
            } else {
                items.swap(ci, i);
            }
            f(items)?;
            count += 1;
            if let Some(slot) = c.get_mut(i) {
                *slot = ci + 1;
            }
            i = 0;
        } else {
            if let Some(slot) = c.get_mut(i) {
                *slot = 0;
            }
            i += 1;
        }
    }
    Ok(count)
}

/// A Fisher-Yates shuffle of `0..SET.len()` from `seed`.
fn shuffled(seed: u64) -> Result<Vec<usize>, String> {
    let mut v: Vec<usize> = (0..SET.len()).collect();
    let mut rng = SplitMix64::new(seed);
    let mut i = v.len();
    while i > 1 {
        let bound = u64::try_from(i).map_err(|e| format!("{:?}", e))?;
        let j = usize::try_from(rng.below(bound)).map_err(|e| format!("{:?}", e))?;
        i -= 1;
        v.swap(i, j);
    }
    Ok(v)
}

#[test]
fn mag_sort_of_the_subset_gives_the_hand_order_under_all_720_permutations() -> Result<(), String> {
    let want = expected_keys(&SUBSET)?;
    let mut items = SUBSET;
    let n = for_each_permutation(&mut items, &mut |p| {
        expect_eq(
            &format!("sort of permutation {:?}", p),
            &sorted_keys(p)?,
            &want,
        )
    })?;
    expect_eq("permutations", &n, &720)
}

#[test]
fn mag_sort_of_seeded_shuffles_of_the_whole_set_gives_the_hand_order() -> Result<(), String> {
    let all: Vec<usize> = (0..SET.len()).collect();
    let want = expected_keys(&all)?;
    for seed in 0..SHUFFLES {
        let p = shuffled(seed)?;
        expect_eq(
            &format!("sort of shuffle seed {} {:?}", seed, p),
            &sorted_keys(&p)?,
            &want,
        )?;
    }
    Ok(())
}

#[test]
fn mag_hand_order_is_strictly_ascending_apart_from_the_two_equal_pairs() -> Result<(), String> {
    let mut v = Vec::new();
    for s in &SET {
        v.push(build(s)?);
    }
    let mut equal = 0u32;
    for (a, b) in v.iter().zip(v.iter().skip(1)) {
        match Ord::cmp(a, b) {
            core::cmp::Ordering::Less => {}
            core::cmp::Ordering::Equal => equal += 1,
            core::cmp::Ordering::Greater => {
                return Err(format!("{:?} > {:?} in the hand order", a, b))
            }
        }
    }
    expect_eq("equal neighbours (two zeros, the duplicate)", &equal, &2)
}

#[test]
fn shuffles_are_deterministic_and_reach_more_than_one_order() -> Result<(), String> {
    expect_eq("seed 7 twice", &shuffled(7)?, &shuffled(7)?)?;
    let first = shuffled(0)?;
    let mut distinct = 0u64;
    for seed in 1..SHUFFLES {
        if shuffled(seed)? != first {
            distinct += 1;
        }
    }
    expect_eq("shuffles that differ from seed 0", &(distinct > 400), &true)
}
