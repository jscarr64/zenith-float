## Integration tests

The `mpfr` directory compares zenith-float arithmetic and elementary functions with MPFR at the bit level.

These tests are off by default. They require Linux x86_64 and the `rug` / `gmp-mpfr-sys` stack:

```bash
cargo test -p zenith-float-num --features mpfr-tests -- --test-threads=1
```

`scripts/ci.sh` runs them in **release** on Linux x86_64 automatically.

## Radix round-trip stress

`tests/radix_roundtrip.rs` exercises parse/format for all bases 2–36 (`random` feature):

```bash
cargo test -p zenith-float-num --features random --test radix_roundtrip
```
