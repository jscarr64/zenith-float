# Reproducibility

This note is for citing a zenith-float result in a paper. It is not a tutorial. See [GETTING_STARTED.md](GETTING_STARTED.md) for a first computation and [LIBRARY.md](LIBRARY.md) for method contracts.

## What determines a result

A finite `ExactNum` or `ExactComplex` value is determined by:

- the source of this crate at a given commit
- the input values
- the requested precision `p`
- the rounding mode `rm`
- for operations that take a cache, the `Consts` (or `SharedConsts`) state used to evaluate named constants

There is no platform-dependent numeric path. Arithmetic is software integer limbs. Hardware binary interchange types are not used. Software `Ieee32` / `Ieee64` are also integer-bit kernels (`from_bits` / `to_bits`); they are not the hardware FPU.

The same inputs at the same commit produce the same limbs on every host that can build the crate. Cross-width (`WORD_BIT_SIZE = 32` vs `64`) and hex-limb CI (`scripts/ci_hex_arm.sh`, `ci_hex_wasm.sh`, `ci_hex_32bit.sh`) lock `to_bytes()` hex in `golds/hex/reference.txt`.

## How to replay

Default unit tests:

```bash
cargo test --workspace
```

Full gold suite, including optional MPFR bit-oracles (Linux x86_64, `rug`):

```bash
cargo test --features mpfr-tests
```

`scripts/ci.sh` / `scripts/ci_full.sh` are the maintainer gates. Debug tests and the MPFR gate have wall-time budgets (`CI_DEBUG_SECS` / `CI_MPFR_SECS`).

Random tests use a seeded software RNG. The default seed is `0x5EED_CAFE_BADC_0D00`. A panic prints the seed. Replay:

```bash
ZENITH_TEST_SEED=<printed value> cargo test <name> -- --test-threads=1
```

## Where the locked values live

Expected values are locked in the crate’s unit tests (`zenith-float-num` lib tests and `tests/`) and in `golds/hex/reference.txt` (canonical `to_bytes()` hex, not a tolerance).

## Citation

Use the crate version from `Cargo.toml` and the git commit you ran:

> Computed with zenith-float 1.0.1 (commit `{hash}`), reproducible by running `cargo test` at that commit.

Replace `{hash}` with `git rev-parse HEAD` of the tree you used. If you enabled `mpfr-tests`, say so.
