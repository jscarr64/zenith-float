# Contributing to zenith-bound

zenith-bound is pre-release. Every rule below is enforced in review, and the ones marked **(gate)** are also enforced by `tools/full_gate.sh`. The plan is [`docs/PLAN.md`](docs/PLAN.md); the correctness contract is [`docs/CONTRACT.md`](docs/CONTRACT.md) (if they disagree, the contract wins).

## Hard rules

These apply to every line of Rust in the crate: `src/`, `tests/`, `benches/`, `examples/`, `build.rs`, test helpers, `#[cfg(test)]` items and generator output that is Rust source. **There is no test-code exemption.**

1. **100% Rust, no `unsafe`. (gate)** `#![forbid(unsafe_code)]` stays in `src/lib.rs`, and `[lints.rust] unsafe_code = "forbid"` in `Cargo.toml` covers every target. The only `unsafe` the project may ever contain is the FPU-rounding behaviour test in a separate dev-only crate `fpu-check/` (`publish = false`, planned, not yet present), which zenith-bound never depends on.
2. **No IEEE floats. (gate)** No `f32`/`f64` anywhere (not in tests, not in generated Rust), no float literals typed as such, no `as f64`, no `mul_add`/FMA, no FPU rounding-mode control, no `core::arch`/`std::arch`, no `#[target_feature]`, no CPU-feature detection. Directed rounding comes only from zenith-float's `RoundingMode::Up`/`Down`, and radii come from the integer-only `Mag` type. Python generators under `scripts/` may use floats internally only to run the oracle, and their outputs are data files (integers and decimal strings), never Rust source with floats.
3. **Zero clone. (gate)** No `.clone()`, `.clone_from()`, `Clone::clone`, `.to_vec()`, `.cloned()`, `.to_owned()`, no `impl Clone`, and no `#[derive(Clone)]`. **`Copy` is not approved either** (`derive(Copy)` needs `Clone`, and `impl Copy` is flagged): `Mag`, `Tri` and every other type are used through references, and `Mag` methods take `&Mag`. If a design seems to need a copy, consume the value, restructure, or ask. An exception needs a concrete proof and Jeff's approval, and is recorded in `tools/rules_allowlist.txt`, which only Jeff edits.
4. **Zero `unwrap`, `expect`, `panic!`, `unreachable!`, `todo!`, `unimplemented!`. (gate)** Banned in library code **and in tests, test helpers and benches**. Tests are written as `fn name() -> Result<(), String>` functions that return `Err(...)` with a message (a small helper such as `check(cond, msg) -> Result<(), String>` and `?`), never with `unwrap`/`expect`. The only exception is an `unwrap` that is **absolutely necessary**: it is listed in `tools/rules_allowlist.txt` with its reason, for Jeff, who alone approves and edits that file. Library code also has no slice/array indexing that can fail and no unchecked `+ - * <<` on exponents (use `checked_*`); every `Option`/`Result` is propagated into `BoundError`/`Indeterminate`. `assert!`, `assert_eq!`, `assert_ne!`, `assert_matches!` and their `debug_` forms are banned everywhere (a failed check returns `Err`, never panics), as are `unwrap_or_default(`, `unwrap_unchecked(`, `unwrap_err(` and `expect_err(`; `tools/rules_scan.py` and the gate's `no-panic-scan` flag them.
5. **No `#[allow]` and no `#[expect]`, no dead code, no commented-out code. (gate)** Lints are denied for every target in `Cargo.toml` `[lints]` and `clippy.toml`; the gate runs `cargo clippy --all-targets -- -D warnings`.
6. **No new dependencies without Rusty's approval.** Any pin of zenith-float is a git revision in its own commit, filled in by Rusty once the final zenith-float 2.0.0 commit exists (`publish = false` stays). zenith-bound makes **no zenith-float call sites** until that commit exists.
7. **zenith-float ownership forms.** Prefer consuming operands (`self`) and other zero-clone forms. `add_ref`/`sub_ref` are allowed only at sites where both operands are still live. `mul_ref`/`div_ref` are approved **only if absolutely necessary, when everything else fails**, and every use is listed with its reason, for Jeff. `sqr(&self)` is fine. The scan counts every use of `add_ref(`, `sub_ref(`, `mul_ref(`, `div_ref(` against `tools/ref_sum_count_baseline.txt`, which starts at 0 and only goes down (a raise needs Jeff, as a change to a protected file). **(gate)**
8. **MSRV is the frozen pin, Rust 1.99** (`rust-toolchain.toml` with `clippy` and `rustfmt`, `rust-version = "1.99"`). Do not use features newer than 1.99. The pin moves only by a deliberate decision, in a minor (0.x) or major release, named in `CHANGELOG.md`.
9. **License:** contributions are accepted under `MIT OR Apache-2.0`. Files ported from puremp (MIT, KarpelesLab) keep an SPDX header naming the puremp commit and the MIT notice with credit to KarpelesLab (see `NOTICE-puremp`).

## Protected files

Only Jeff changes these; a diff touching them is rejected and sent to him: `.github/CODEOWNERS`, `tools/*` (gate, `rules_scan.*`, self-test), `tools/rules_allowlist.txt`, `tools/clone_count_baseline.txt`, `tools/ref_sum_count_baseline.txt`, `docs/CONTRACT.md` clauses 1 and 5 (LOCKED), `rust-toolchain.toml`, the zenith-float `rev` in `Cargo.toml`, `NOTICE*`, `LICENSE*`. `.github/CODEOWNERS` is `* @jscarr64`.

## The gate

`tools/full_gate.sh TAG` runs, each step logged under `/tmp/zb-gate-logs/`: `cargo test --no-run` (`build-tests`: a nonzero rc here is a compile error, not a red test), `cargo test` (debug, release), the rules-scan self-test and `tools/rules_scan.sh` over the commit range, `cargo fmt --check`, `cargo clippy --all-targets -- -D warnings`, and grep scans for unsafe, floats, panics (including the assert family) and clones (`derive-scan` also sees a multi-line derive). On Pop run it as `flock /tmp/pop-gate.lock tools/full_gate.sh TAG` with `CARGO_TARGET_DIR=/tmp/target-<agent>-<worktree> CARGO_BUILD_JOBS=7`. `CARGO_TARGET_DIR` is required and belongs to one worktree (Paul, 2026-10-23): the `target-dir` step (and `tools/fresh_test_run.sh`, which wraps `build-debug`/`build-release`) records the owning worktree in `$CARGO_TARGET_DIR/.worktree-owner` and fails if the variable is unset or the dir is owned by another worktree. The library crate also denies `missing_docs`, `clippy::allow_attributes` and `clippy::arithmetic_side_effects` (`src/lib.rs`). Run it on **every** commit, not just the last. `tools/rules_scan.sh --working` (after `git add -A`) scans uncommitted work. The scan was derived from zenith-float's (`d12a156`), with the test-code exemption removed.

Every run on Pop gets a run log (`~/pop-logs/`): `tools/*.sh` scripts log themselves through `tools/runlog.sh`; wrap any other verify or test run as `tools/runlog.sh JOB [--tests] -- CMD`.

**Not yet present** (earlier drafts of this file promised them; they are planned, and each lands in its own reviewed commit): `scripts/check_no_float.sh`, a HIR scan, the `thumbv6m-none-eabi` no-FPU build and symbol check (which also needs a decision on `#![no_std]`), `fpu-check/`, `cargo-deny`, `cargo-semver-checks`, CI workflows, `docs/audits/`. Until they exist, nothing here claims to be checked by them.

## Clean-room rule (LGPL oracles)

Arb/FLINT and MPFR are licensed under the **LGPL**. We use them (and mpmath) only as **dev-only test oracles**.

- Implement algorithms **from the published papers and textbooks only** (for example Johansson's Arb papers, Petras, Rump, Lohner, Olver, Carlson). **Never read, copy, translate or paraphrase Arb/FLINT or MPFR source code** while working on zenith-bound, and don't paste their code into prompts for coding models.
- Cite the paper in the module doc of anything that implements a published bound.
- Oracles are never a dependency in any `Cargo.toml` of the published crate, and they are never shipped. `cargo-deny` (planned) will reject LGPL licenses in the dependency graph.
- Goldens derived from oracle output live under `tests/` and `scripts/` and are excluded from the published package (see `NOTICE`).

## Tests

- **Red first.** Every change starts as failing tests that are committed (or shown in review) before the fix. A red state is a test that fails for the right reason against a stub that compiles; stubs never use `todo!()`, `unimplemented!()` or `panic!()`.
- **No single exact-input tests.** Every test is a **property set**, an **oracle set**, an **adversarial set** (a named, enumerated generator) or an **overestimation measurement**. A one-off `assert_eq!(f(3), …)` is rejected. proptest regression files keep failing seeds.
- **Zero is in every set**: zero midpoint, zero radius, and both.
- **Truth must be independent of zenith-float.** Never use zf's `*_full_prec` ops or any other zf result as truth, because that is circular. Use the dev-only rational/bigint model (`tests/model/`) or Python-generated goldens.
- **Every test cites its CONTRACT clause number and names its oracle.**
- **Containment failures are P0** (and security-grade, see `SECURITY.md`). Tightness misses are performance bugs and are tracked as overestimation ratios.
- **Mutation checks:** shrinking any radius by 1 ulp (`cfg(zb_shrink_radius)`) must make the suite fail, and `cargo-mutants` runs nightly on `src/mag`, `src/ball` and `src/encode.rs` (both planned).

## Oracle setup (dev only)

```sh
~/.local/oracle-venv/bin/python scripts/gen_goldens.py    # writes tests/goldens/mag_ops.txt
```

Use the isolated oracle venv `~/.local/oracle-venv/bin/python` (mpmath 1.3.0, sympy 1.14.0, gmpy2 2.3.2, python-flint 0.9.0; `scripts/requirements.lock`). `scripts/goldlib.py` refuses any other interpreter at startup: exact versions, `sys.prefix` is the venv, no `~/.local` on `sys.path`. Never import oracle packages from `~/.local` outside the venv. The output is deterministic, so regenerating must leave `tests/goldens/mag_ops.txt` byte-identical.

MPFR is reached only through zenith-float's own test suite (via `rug`); zenith-bound has no MPFR dependency. Goldens are regenerated only in a commit that says why, and a byte-golden change also needs an `op_ver` bump.

## Review

- Each diff is hunk-sized and gets a bot review. Until 2026-10-30, non-trivial code is written by the local model qwen3-coder:30b (temperature at most 0.1, low top_p) in small one-file steps, and the committer reviews every diff.
- **Every file under `src/mag/` also needs one review by a stronger reviewer (Paul or Cassie)** before it merges.
- Nothing is pushed, tagged, merged or published without Jeff.
- The README status table may call a function **Guaranteed** only after its contract tests pass in CI.
