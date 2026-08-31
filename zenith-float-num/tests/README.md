## Integration tests

The `mpfr` directory compares zenith-float arithmetic and specials with MPFR and GMP. GNU MPC 1.3 has no `erf`, `gamma`, Airy, or Bessel — those use a real-axis MPFR restriction or an identity gold. Do not add ARB or link MPC for a function it does not have.

| Function | Oracle |
| --- | --- |
| `ExactNum` arith / elem / `erf` / `Γ` / `ψ` / `Ai` / `Ei` / `J_n` / `Y_n` / `Γ(s,x)` | MPFR |
| Complex `erf` `Γ` `Ai` `J_n` on the real axis | MPFR |
| `ExactRational` `+ − × ÷` | GMP `mpq` |
| `Si` `Ci` `li` Fresnel `Bi` `_2F1` elliptic `K` `I_ν` `K_ν` | identity (no MPFR function) |

These tests are off by default. They require Linux x86_64 and the `rug` / `gmp-mpfr-sys` stack:

```bash
cargo test -p zenith-float-num --features mpfr-tests -- --test-threads=1
```

`scripts/ci.sh` runs them in **release** on Linux x86_64 automatically. That gate includes `fuzz_round_modes` (add / mul / sqrt vs MPFR for every IEEE rounding mode) and `compare_copysign_test` (including ~16k-bit precision). Random tests share seed `0x5EED_CAFE_BADC_0D00` (`ZENITH_TEST_SEED` to override).

## Radix round-trip stress

`tests/radix_roundtrip.rs` exercises parse/format for all bases 2–36 (`random` feature):

```bash
cargo test -p zenith-float-num --features random --test radix_roundtrip
```
