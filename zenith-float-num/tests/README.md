## Integration tests

The `mpfr` directory compares zenith-float arithmetic and elementary functions with MPFR at the bit level.

These tests are off by default. They require Linux x86_64 and the `rug` / `gmp-mpfr-sys` stack:

```bash
cargo test -p zenith-float-num --features mpfr-tests -- --test-threads=1
```
