# zenith-float build checklist

Living document for what is **implemented**, **tested**, and **still required** before zenith-float can serve as the numeric kernel for Accumath (~675K equations and formulas).

**Last updated:** 2026-08-28  
**Crate version:** 0.1.0 (+ unreleased changelog items)  
**Reference versions (crates.io):** astro-float 0.9.6, dashu-float 0.6.0  
**Policy:** No hardware floating-point in calculations or identifiers (`f32`/`f64` forbidden in source and docs; enforced in `scripts/ci.sh`).

**Design goal:** As **broad** an API as practical (match or exceed astro-float / dashu-float coverage where Accumath needs it), but **always software limbs** — no machine floating-point registers, literals, or converters inside the library. Callers that need IEEE interchange do conversion outside zenith-float.

---

## How to use this checklist

| Symbol | Meaning |
|--------|---------|
| ✅ | Built, in public API, and covered by default CI |
| 🟡 | Built but partial coverage, manual gate, or known limitations |
| ⬜ | Not implemented or not production-ready |
| 🚫 | Explicitly out of scope (Accumath / caller responsibility) |

**Verify locally**

```bash
./scripts/ci.sh
cargo test -p zenith-float-num --features mpfr-tests -- --test-threads=1   # Linux x86_64 + rug
./scripts/bench.sh                                                          # Criterion (--quick)
```

---

## 1. What is built and running today

### 1.1 Crates

| Crate | Role | Status |
|-------|------|--------|
| `zenith-float` | Public facade, docs, re-exports | ✅ |
| `zenith-float-num` | Limb arithmetic, transcendentals, I/O | ✅ |
| `zenith-float-macro` | `expr!` proc-macro | ✅ |

### 1.2 Core numeric model

| Item | Status | Notes |
|------|--------|-------|
| Limb mantissa (`Word` = 64-bit on 64-bit targets) | ✅ | Word-aligned precision |
| Configurable exponent range (`EXPONENT_MIN` / `EXPONENT_MAX`) | ✅ | Reduced on 32-bit |
| Seven rounding modes | ✅ | `None`, `Up`, `Down`, `ToZero`, `FromZero`, `ToEven`, `ToOdd` |
| `±Inf`, `NaN`, subnormals | ✅ | Software values; errors map to `NaN` at API boundary |
| `inexact` flag + correct-rounding retry loops | ✅ | `try_set_precision` / `bump_prec_retry` |
| `no_std` + global allocator | ✅ | Feature `std` optional |
| No hardware float in kernel | ✅ | `#![deny(clippy::float_arithmetic)]` + CI grep |

### 1.3 Arithmetic (mantissa layer)

| Operation | Algorithm | Unit tests | MPFR oracle |
|-----------|-----------|------------|-------------|
| Add / sub | Limb carry/borrow | ✅ | ✅ |
| Mul | Schoolbook → Toom-2 → Toom-3 → FFT | ✅ | ✅ |
| Div | Knuth-style | ✅ | ✅ |
| `rem` (fmod-like) | ✅ | ✅ | ✅ |
| `sqrt` | Digit algorithm | ✅ | ✅ |
| `cbrt` | Root estimate + refinement | ✅ | ✅ |
| `reciprocal` | Schoolbook (&lt;3 words); Newton (≥3 words) | ✅ | ✅ |
| `powi` / `powsi` | Binary exponentiation | ✅ | ✅ |

### 1.4 Public `ExactNum` API — arithmetic & utility

| API | `expr!` | Status |
|-----|---------|--------|
| `add`, `sub`, `mul`, `div` (+ `_full_prec`) | `+ − * /` | ✅ |
| `rem` | `%` | ✅ |
| `reciprocal` | `recip` | ✅ |
| `neg`, `abs`, `signum` | unary `-` | ✅ |
| `pow`, `powi`, `powsi` | `pow` | ✅ |
| `sqrt`, `cbrt` | same | ✅ |
| `int`, `fract`, `ceil`, `floor`, `round` | — | ✅ |
| `min`, `max`, `clamp` | — | ✅ |
| `cmp`, `abs_cmp`, `PartialEq`/`Ord` | — | ✅ (NaN unordered) |
| `hypot` | `hypot` | ✅ |
| `set_precision`, `try_set_precision` | via `Context` | ✅ |

### 1.5 Transcendentals & constants

| API | `expr!` | MPFR | Notes |
|-----|---------|------|-------|
| `ln`, `log2`, `log10`, `log` (arbitrary base) | ✅ | ✅ | Series + arg reduction |
| `log1p` | ✅ | ✅ | |
| `exp` | ✅ | ✅ | |
| `exp2`, `exp10` | ✅ | 🟡 | Delegate to `pow`; **no dedicated MPFR tests yet** |
| `expm1` | ✅ | ✅ | |
| `sin`, `cos`, `tan` | ✅ | ✅ | `rem_pi` internally |
| `asin`, `acos`, `atan`, `atan2` | ✅ | ✅ | |
| `rem_pi` | ✅ | 🟡 | Public wrapper; **no MPFR tests yet** |
| `sinh`, `cosh`, `tanh` | ✅ | ✅ | |
| `asinh`, `acosh`, `atanh` | ✅ | ✅ | `asinh` fixed for large \|x\| (2\|e\| extra bits) |
| Constants: π, e, ln 2, ln 10 (`Consts`) | `pi`, `e`, `ln_2`, `ln_10` | ✅ | Lazy cache |

### 1.6 I/O and integration

| Item | Status |
|------|--------|
| Parse/format: binary, octal, decimal, hexadecimal | ✅ |
| `convert_from_radix` / `convert_to_radix` | ✅ |
| `Display`, `Binary`, `Octal`, `UpperHex` (`std`) | ✅ |
| `FromStr` (decimal, `std`) | ✅ |
| `serde` (decimal string / integer, feature) | 🟡 |
| `random_normal` (feature `random`) | ✅ |
| `expr!` compile tests (`trybuild`) | ✅ |
| `expr!` cancellation / precision tests (root `tests/`) | ✅ |

### 1.7 Test & CI inventory (default `scripts/ci.sh`)

| Gate | Status |
|------|--------|
| Forbid `f32`/`f64` identifiers in `.rs` / `.md` / `CHANGELOG` | ✅ |
| `cargo test --workspace` (debug) | ✅ ~60 lib tests pass |
| `cargo test -p zenith-float-num --lib --release` | ✅ |
| `cargo test` with `no-default-features --features std` | ✅ |
| `cargo test --features random,serde` | ✅ |
| MPFR bit-oracle tests (`mpfr-tests`) | 🟡 Manual; **not in CI** |
| Criterion / dedicated benches | ✅ `zenith-float-num/benches/` (arithmetic, transcendentals, composite) |
| `proptest` / quickcheck | ⬜ Hand-written random loops (`TEST_ITERS = 256`) |

**Property tests** (`zenith-float-num/src/ops/tests.rs`): inverse pairs (ln↔exp, sin↔asin, log↔pow, etc.) with mathematically derived error bounds; exponent sampling capped at `TEST_EXP_BOUND = 1024` for runtime.

**Error documentation:** `doc/README.md` (ulp / series error bounds used by `expr!` and tests).

---

## 2. Comparison vs astro-float and dashu-float

zenith-float is an **evolution of the astro-float design** (`BigFloat` → `ExactNum`, same `Context` / `expr!` / `Consts` / MPFR harness). [astro-float 0.9.6](https://docs.rs/astro-float) is the direct ancestor. [dashu-float 0.6.0](https://docs.rs/dashu-float) is a **different architecture**: arbitrary base, native operators, Ziv-certified transcendentals, optional complex (`CBig`).

Local reference trees (`dashu-master/`, `astro-float-main/`) are not in this repo; compare via crates.io sources or Accumath `target/` build fingerprints.

### 2.1 Feature matrix (high level)

| Capability | zenith-float | astro-float 0.9.6 | dashu-float 0.6.0 |
|------------|:------------:|:-----------------:|:-----------------:|
| Pure Rust kernel (no MPFR in lib) | ✅ | ✅ | ✅ |
| `no_std` + allocator | ✅ | ✅ | ✅ |
| `expr!` + `Context` | ✅ | ✅ | 🟡 (macros elsewhere) |
| Hardware `f32`/`f64` in library | 🚫 forbidden | `from_f32` / `from_f64` | API-edge literals |
| Elementary transcendentals (core set) | ✅ | ✅ | ✅ |
| `hypot`, `atan2`, `log1p`, `expm1` | ✅ | ⬜ not in astro | ✅ |
| `exp2`, `exp10`, `rem_pi` | ✅ | ⬜ not in astro | 🟡 partial |
| Native `+ − × ÷` operators | ⬜ | ⬜ | ✅ |
| `fbig!` / compile-time float literals | ⬜ | ⬜ | ✅ |
| Parse/format bases 2–36 | ⬜ (bin/oct/dec/hex) | ⬜ (bin/oct/dec/hex) | ✅ |
| Arbitrary-base float type | ⬜ | ⬜ | ✅ |
| Complex (`CBig`) | ⬜ | ⬜ | ✅ |
| `fma` / `mul_add` | ⬜ | ⬜ | 🟡 |
| General `nth_root(n)` | ⬜ | ⬜ | ✅ |
| `sin_cos` / `sinh_cosh` paired APIs | ⬜ | ⬜ | ✅ |
| Special functions (erf, Γ, Bessel, …) | ⬜ | ⬜ | 🟡 / separate |
| MPFR golden tests in repo | ✅ (optional) | ✅ (optional) | fuzz + unit (project policy) |
| Fuzz MPFR bit-exact (all round modes) | ⬜ | ⬜ | ✅ |
| Progressive constant cache | 🟡 series cache | 🟡 | ✅ `ConstCache` / `CachedFBig` |
| Ziv + Ball correct-rounding proof | ⬜ | ⬜ | ✅ on transcendentals |
| `serde` / `random` | ✅ optional | ✅ default-on | ✅ optional |
| Stack-inlined small values | ⬜ | ⬜ | ✅ |
| Published crate + bench history | 🟡 0.1.0 | ✅ 0.9.x | ✅ |

### 2.2 zenith-float vs astro-float (lineage)

zenith-float is already a **superset** of astro-float’s math API. astro-float does **not** exceed zenith on transcendentals.

| zenith has; astro 0.9.6 lacks | Notes |
|------------------------------|-------|
| `exp2`, `exp10` | zenith delegates to `pow` |
| `log1p`, `expm1` | |
| `atan2`, `hypot` | |
| `rem_pi` (public) | astro reduces internally only |
| No `from_f32` / `from_f64` | zenith policy |
| `random` / `serde` off by default | astro enables by default |

Shared weakness vs dashu: no Ziv loop, no guaranteed correctly-rounded transcendentals from `expr!`, binary mantissa only.

### 2.3 Where zenith-float must go beyond both (Accumath scale)

For **675K diverse formulas**, correctness and operability matter more than matching every dashu feature on day one.

| Priority | Gap | Why it matters |
|----------|-----|----------------|
| P0 | MPFR oracle in release CI (or nightly) | Bit-exact reference for every public op |
| P0 | MPFR coverage for `exp2`, `exp10`, `rem_pi` | Recently added; still oracle-blind |
| P0 | Working-precision caps documented + enforced | Prevents pathological hour-long single test cases |
| P1 | `fma` | Exact dot products, compensated summation in long expressions |
| P1 | `nth_root(n)` | Many physics / engineering closed forms |
| P1 | `sinh_cosh` paired evaluation | dashu has this; saves duplicate exp work |
| P1 | Fuzz MPFR differential (all rounding modes) | dashu `fuzz/` model; stronger than property loops alone |
| P1 | Benchmark regression tracking in CI (optional nightly) | Criterion benches exist; not yet gated in `scripts/ci.sh` |
| P1 | `expr!` correct-rounding semantics documented per op | Accumath will lean on macro heavily |
| P2 | Extra constants (φ, √2, γ, …) | Fewer series cold-starts |
| P2 | `LowerExp` / `UpperExp` formatting | Debug / log output at scale |
| P2 | `frexp` / `ldexp` / `scalb` / `logb` | IEEE-style interoperability without hardware floats |
| P3 | Special functions (erf, Gamma, …) | Only where symbolic engine cannot stay exact |
| P3 | Parallel evaluation / thread-safe shared `Consts` | Batch numeric evaluation |

---

## 3. Roadmap checklist

### 3.1 Kernel correctness (P0)

- [x] Software-limb float; no hardware float arithmetic in kernel
- [x] Full elementary + inverse trig and hyperbolic set
- [x] `hypot`, `atan2`, `log1p`, `expm1`
- [x] `exp2`, `exp10`, `rem_pi`
- [x] Newton reciprocal (≥3 words)
- [x] Decimal/radix I/O without hardware floats
- [x] `asinh` large-argument precision (`2|e|` bits + `try_set_precision(p_x)`)
- [x] Property tests without debug-only weakening
- [x] CI: debug + release lib tests, `f32`/`f64` grep
- [ ] Add `exp2`, `exp10`, `rem_pi` to MPFR compare suite
- [ ] MPFR gate in CI (nightly or release-only job; Linux x86_64)
- [ ] Document / bound internal precision growth vs exponent (avoid runaway cost)
- [ ] Fuzz or structured stress harness for parse/format round-trip

### 3.2 API completeness vs dashu/astro (P1)

- [ ] `fma` (`mul_add`) with correct rounding
- [ ] `nth_root(n)` for n ≥ 2 (generalize `sqrt` / `cbrt`)
- [ ] `sinh_cosh` paired evaluation (single exp path)
- [ ] Fuzz harness: MPFR bit-exact under all rounding modes (dashu-style)
- [ ] `copysign`, `next_after` (software, integer-based)
- [ ] `frexp` / `ldexp` / `scalb` / `logb` / `ilogb`
- [ ] `LowerHex` + scientific `Display` options
- [ ] Additional constants in `Consts` (configurable list)

### 3.3 Testing & performance (P1)

- [x] Replace `#[ignore]` `*_perf` with Criterion benches in `benches/`
- [ ] Optional: nightly Criterion regression gate in CI
- [ ] Track CI wall time budget (target: full gate &lt; 10 min debug, &lt; 30 min with MPFR)
- [ ] Seed-controlled random tests for reproducible failures
- [ ] Accumath-driven regression corpus (golden files from real formula subsets)
- [ ] Memory high-water tests for large precisions (OOM → clean `NaN`, not hang)

### 3.4 Accumath integration (P1–P2)

- [ ] Stable ABI surface for engine (`ExactNum` eval from `CanonicalExpr`)
- [ ] Batch evaluator with shared `Consts` / `Context` reuse
- [ ] Precision policy per formula class (exact vs numeric terminal)
- [ ] Wire 675K-formula smoke + deep regression in Accumath CI (separate repo)
- [ ] Purge `f64` / `libm` from Accumath numeric path (🚫 not this crate)

### 3.5 Out of scope for zenith-float

- [ ] 🚫 Hardware `f32`/`f64` converters inside the library
- [ ] 🚫 Arbitrary-complex arithmetic (use separate type or dashu `CBig` if needed)
- [ ] 🚫 Computer algebra (Risch, towers, etc.) — Accumath symbolic layer
- [ ] 🚫 Replacing MPFR/GMP at test time (MPFR remains oracle only)

---

## 4. Release readiness gates

Before calling a version **production-ready for Accumath numeric evaluation**:

1. **All P0 items checked** (including MPFR for every exported transcendental).
2. **`./scripts/ci.sh` green** on Linux x86_64, debug + release.
3. **MPFR suite green** on same platform (`--features mpfr-tests`).
4. **No known hang** in property tests at `TEST_EXP_BOUND` and `TEST_ITERS = 256`.
5. **CHANGELOG + README** match public API (no stale checklist).
6. **Accumath subset**: at least one real formula corpus (1K–10K expressions) evaluated end-to-end through `expr!` or `ExactNum` without `f64`.

---

## 5. Related files

| Path | Purpose |
|------|---------|
| `scripts/ci.sh` | Default public CI |
| `zenith-float-num/tests/README.md` | MPFR test instructions |
| `zenith-float-num/src/ops/tests.rs` | Random inverse property tests |
| `zenith-float-num/tests/mpfr/` | MPFR bit-oracle tests |
| `doc/README.md` | Error bound theory |
| `tests/mod.rs` | `expr!` integration tests |
| `CHANGELOG.md` | Released / unreleased API notes |

---

## 6. Changelog sync (unreleased)

Already implemented but not in a crates.io release:

- Newton reciprocal at three words and up
- `hypot`, `atan2`, `log1p`, `expm1`
- `exp2`, `exp10`, `rem_pi`
- CI hardening and `asinh` precision fix (2026-08-27)
