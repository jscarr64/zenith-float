# Contributing to zenith-float

Thank you for looking at the crate. This is an arbitrary-precision **software** floating-point library. Version **1.0** is the stable public API (`ExactNum`, `ExactComplex`, `Ieee32` / `Ieee64`, `expr!` / `cexpr!`). Breaking changes would be a new major version. Issues and pull requests are welcome when they match that job.

Depend on **`zenith-float`**, not `zenith-float-num`. The kernel crate is an implementation detail.

## What this crate is

- Integer-limb arithmetic only. Hardware IEEE arithmetic is not used in calculations.
- IEEE widths are software types `Ieee32` / `Ieee64` (`from_bits` / `to_bits`). There is no `From` of a hardware float.
- Incomplete paths return `NaN` with a named `Error`, or `None`. They must not panic and must not invent values.

## Before you write code

Open an **issue** first for anything larger than a typo. Say which method you called, the precision and rounding mode, what you expected, and what you got (`Error` variant, bits, or decimal).

## Development

```bash
bash scripts/ci.sh
bash scripts/ci_full.sh
```

Enable the committed pre-commit hook once per clone: `git config core.hooksPath .githooks`. It runs
`tools/rules_scan.sh --working`, `tools/rules_scan.sh --count-only --no-rise HEAD`, `cargo fmt --all -- --check` and
`cargo clippy --workspace --all-targets --all-features -- -D warnings`, and refuses the commit on any FAIL. The same
setting enables `.githooks/commit-msg`, which refuses a commit whose message has no `Agent: <name>` trailer in its last
paragraph (merges too). Never commit with `--no-verify`.

`scripts/ci.sh` is the default gate (workspace tests, `no_std` / thumb, MPFR on Linux x86_64). `scripts/ci_full.sh` is the 12-check pre-publish gate.

MPFR golds need Linux x86_64 and `rug` (`--features mpfr-tests`). Replay a random failure with `ZENITH_TEST_SEED=<seed> cargo test <name> -- --test-threads=1`.

## Per-commit checklist

Every commit, before it is made (the pre-commit hook runs items 3 to 5; run the rest yourself and keep the logs):

1. `cargo check -p zenith-float`, `cargo check -p zenith-float --no-default-features`, `cargo check -p zenith-float-num`,
   `cargo check -p zenith-float-num --no-default-features`.
2. `cargo fmt --all -- --check` and `cargo clippy --workspace --all-targets --all-features -- -D warnings`.
3. `tools/rules_scan.sh --working` (the library scanner, at the repo root).
4. `tools/rules_scan.sh --count-only` and its log, which names the commit. Above baseline the hook allows a commit only
   when no count rises (`--no-rise HEAD`); that is cleanup, not clean. Never raise a baseline.
5. No new clone and no new panic, in library or test code (tools/NO_CLONE_RULES.md; tests return `Result` and use `?`).

Every 'clean' claim cites a log from the library's own `tools/rules_scan.sh` at that exact commit (Jeff 2026-11-23 07:43
CT). The bound layer's `zenith-float-bound/tools/bound_layer_scan.sh` is never evidence for the library.

## Golds

A gold is an expected **bit pattern**, **named constant**, **identity**, or **MPFR 1-ULP**. Tests that only show “did not panic” or `is_ok()` are not enough. Do not weaken existing golds. Do not mark a failing gold `#[ignore]` to go green.

## Pull requests

- Keep the diff to the problem. Do not reformat unrelated files.
- Do not add hardware IEEE types (`f32` / `f64`) in `.rs` sources.
- Do not add CAS, a second float kernel, or a file-format crate.
- Public failures stay named (`Error::…` or `None`), not string-only.
- Run `bash scripts/ci.sh` before you push.

By submitting a change you agree it is licensed under **MIT OR Apache-2.0**, the same as the rest of the crate (`LICENSE-MIT`, `LICENSE-APACHE`).

## Out of scope (will be closed)

- Hardware floating-point arithmetic or `From` of a hardware float
- Computer algebra, formula rewriting, or an expression IR
- HDF5 or any other file-format library inside this crate (CSV and the binary record already here are enough)
- BLAS / blocked matmul
- A total `Ord` / `Hash` that includes NaN

Those belong in a caller, or they contradict the crate contract.
