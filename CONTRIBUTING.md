# Contributing to zenith-float

Thank you for looking at the crate. This is an arbitrary-precision **software** floating-point library. Version **1.0** is the stable public API (`ExactNum`, `ExactComplex`, `Ieee32` / `Ieee64`, `expr!` / `cexpr!`). Breaking changes would be a new major version. Issues and pull requests are welcome when they match that job.

Depend on **`zenith-float`**, not `zenith-float-num`. The kernel crate is an implementation detail.

## What this crate is

- Integer-limb arithmetic only. Hardware IEEE arithmetic is not used in calculations.
- IEEE widths are software types `Ieee32` / `Ieee64` (`from_bits` / `to_bits`). There is no `From` of a hardware float.
- Incomplete paths return `NaN` with a named `Error`, or `None`. They must not panic and must not invent values.
- Feature `hdf5` uses crates.io `hdf5-rust`. There is no `libhdf5`.

## Before you write code

Open an **issue** first for anything larger than a typo. Say which method you called, the precision and rounding mode, what you expected, and what you got (`Error` variant, bits, or decimal).

## Development

```bash
bash scripts/ci.sh
bash scripts/ci_full.sh
```

`scripts/ci.sh` is the default gate (workspace tests, `hdf5`, `no_std` / thumb, MPFR on Linux x86_64). `scripts/ci_full.sh` is the 12-check pre-publish gate.

MPFR golds need Linux x86_64 and `rug` (`--features mpfr-tests`). Replay a random failure with `ZENITH_TEST_SEED=<seed> cargo test <name> -- --test-threads=1`.

## Golds

A gold is an expected **bit pattern**, **named constant**, **identity**, or **MPFR 1-ULP**. Tests that only show “did not panic” or `is_ok()` are not enough. Do not weaken existing golds. Do not mark a failing gold `#[ignore]` to go green.

## Pull requests

- Keep the diff to the problem. Do not reformat unrelated files.
- Do not add hardware IEEE types (`f32` / `f64`) in `.rs` sources.
- Do not add `libhdf5`, CAS, or a second float kernel.
- Public failures stay named (`Error::…` or `None`), not string-only.
- Run `bash scripts/ci.sh` before you push.

By submitting a change you agree it is licensed under **MIT OR Apache-2.0**, the same as the rest of the crate (`LICENSE-MIT`, `LICENSE-APACHE`).

## Out of scope (will be closed)

- Hardware floating-point arithmetic or `From` of a hardware float
- Computer algebra, formula rewriting, or an expression IR
- Linking `libhdf5` or wrapping another HDF5 C library
- BLAS / blocked matmul
- A total `Ord` / `Hash` that includes NaN

Those belong in a caller, or they contradict the crate contract.
