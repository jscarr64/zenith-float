//! CONTRACT 10.1, 15.1 (options unit, Paul's rulings 2026-10-23); oracle: the clauses (SolveOpts { refine_steps },
//! EigOpts { max_sweeps, cluster }, both #[non_exhaustive] + Default; DomainPolicy #[non_exhaustive] { Strict, Clip })
//! COVERS: DomainPolicy::Strict, DomainPolicy::Clip, Default for SolveOpts, Default for EigOpts

use zenith_bound::policy::{DomainPolicy, EigOpts, SolveOpts};

/// The policy source, read to check the `#[non_exhaustive]` markers (a removed marker still compiles in this crate's tests
/// only where they match exhaustively, which these tests deliberately do not).
const POLICY_SRC: &str = include_str!("../../src/policy.rs");

fn marked(item: &str) -> bool {
    POLICY_SRC.split("#[non_exhaustive]").skip(1).any(|after| {
        after
            .lines()
            .take(3)
            .any(|l| l.trim_start().starts_with(item))
    })
}

#[test]
fn solve_opts_default() -> Result<(), String> {
    let o = SolveOpts::default();
    if o.refine_steps != 2 {
        return Err(format!(
            "SolveOpts::default().refine_steps = {}, want 2",
            o.refine_steps
        ));
    }
    Ok(())
}

#[test]
fn eig_opts_default() -> Result<(), String> {
    let o = EigOpts::default();
    if o.max_sweeps != 64 || !o.cluster {
        return Err(format!(
            "EigOpts::default() = {:?}, want max_sweeps 64, cluster true",
            o
        ));
    }
    Ok(())
}

#[test]
fn options_are_set_field_by_field() -> Result<(), String> {
    // #[non_exhaustive] forbids struct literals outside the crate; callers start from Default and set fields.
    let mut s = SolveOpts::default();
    s.refine_steps = 0;
    let mut e = EigOpts::default();
    e.max_sweeps = 1;
    e.cluster = false;
    if s == SolveOpts::default() || e == EigOpts::default() {
        return Err("a changed field must make the options differ from Default".into());
    }
    if s.refine_steps != 0 || e.max_sweeps != 1 || e.cluster {
        return Err(format!("fields did not stick: {:?} {:?}", s, e));
    }
    Ok(())
}

#[test]
fn option_types_are_non_exhaustive() -> Result<(), String> {
    for item in [
        "pub struct SolveOpts",
        "pub struct EigOpts",
        "pub enum DomainPolicy",
    ] {
        if !marked(item) {
            return Err(format!("{} must be #[non_exhaustive]", item));
        }
    }
    Ok(())
}

#[test]
fn domain_policy_keeps_two_variants() -> Result<(), String> {
    let all = [DomainPolicy::Strict, DomainPolicy::Clip];
    if all[0] == all[1] || DomainPolicy::default() != DomainPolicy::Strict {
        return Err("DomainPolicy must be {Strict, Clip} with Strict the Default".into());
    }
    Ok(())
}
