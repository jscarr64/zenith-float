# zenith-float build checklist

**Maintainer-only.** Not included in the crates.io package. User-facing docs are [GETTING_STARTED.md](GETTING_STARTED.md), [HELP.md](HELP.md), and [LIBRARY.md](LIBRARY.md).

Living document for what is **implemented**, **tested**, and **required** for zenith-float as a public software big-float crate. Application engines (formula corpora, expression ABIs, host hardware-float purge) live in those applications, not here.

**Last updated:** 2026-08-29  
**Crate version:** 0.1.0 (+ unreleased changelog items)  
**Reference versions (crates.io):** astro-float 0.9.6, dashu-float 0.6.0  
**Policy:** No hardware floating-point in calculations. Rust hardware IEEE type tokens are forbidden in `.rs` (`scripts/ci.sh`). Software `Ieee32`/`Ieee64` store binary32/binary64 as integer bits.

**Design goal:** As **broad** an API as practical (match or exceed astro-float / dashu-float coverage), but **always software limbs** — no machine floating-point registers or `libm`. IEEE interchange is `from_bits` / `to_bits` on `Ieee32`/`Ieee64`.

---

## How to use this checklist

| Symbol | Meaning |
| -------- | --------- |
| ✅ | Built, in public API, and covered by default CI |
| 🟡 | Built but partial coverage, manual gate, or known limitations |
| ⬜ | Not implemented or not production-ready |
| 🚫 | Explicitly out of scope (caller / host application) |

### Verify locally

```bash
./scripts/ci.sh
cargo test -p zenith-float-num --features mpfr-tests -- --test-threads=1   # Linux x86_64 + rug
./scripts/bench.sh                                                          # Criterion (--quick)
./scripts/bench-compare.sh                                                  # compare to doc/bench-baselines.tsv
./scripts/compare-bench.sh --quick                                          # zenith vs astro (release compare)
```

---

## 1. What is built and running today

### 1.1 Crates

| Crate | Role | Status |
| ------- | ------ | -------- |
| `zenith-float` | Public facade, docs, re-exports | ✅ |
| `zenith-float-num` | Limb arithmetic, transcendentals, I/O | ✅ |
| `zenith-float-macro` | `expr!` proc-macro | ✅ |

### 1.2 Core numeric model

| Item | Status | Notes |
| ------ | -------- | ------- |
| Limb mantissa (`Word` = 64-bit on 64-bit targets) | ✅ | Word-aligned precision |
| Configurable exponent range (`EXPONENT_MIN` / `EXPONENT_MAX`) | ✅ | Reduced on 32-bit |
| Seven rounding modes | ✅ | `None`, `Up`, `Down`, `ToZero`, `FromZero`, `ToEven`, `ToOdd` |
| `±Inf`, `NaN`, subnormals | ✅ | Software values; errors map to `NaN` at API boundary |
| `inexact` flag + correct-rounding retry loops | ✅ | `try_set_precision` / `bump_prec_retry` |
| `no_std` + global allocator | ✅ | Feature `std` optional |
| No hardware float in kernel | ✅ | `#![deny(clippy::float_arithmetic)]` + CI grep |

### 1.3 Arithmetic (mantissa layer)

| Operation | Algorithm | Unit tests | MPFR oracle |
| ----------- | ----------- | ------------ | ------------- |
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
| ----- | --------- | -------- |
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
| ----- | --------- | ------ | ------- |
| `ln`, `log2`, `log10`, `log` (arbitrary base) | ✅ | ✅ | Series + arg reduction |
| `log1p` | ✅ | ✅ | |
| `exp` | ✅ | ✅ | |
| `exp2`, `exp10` | ✅ | ✅ | Delegate to `pow`; bit-oracle in compare suite |
| `expm1` | ✅ | ✅ | |
| `sin`, `cos`, `tan` | ✅ | ✅ | `rem_pi` internally |
| `asin`, `acos`, `atan`, `atan2` | ✅ | ✅ | |
| `rem_pi` | ✅ | ✅ | Identity for \|x\|&lt;4; large args stay in `(-2π, 2π)`; sin/cos oracles cover reduction |
| `sinh`, `cosh`, `tanh` | ✅ | ✅ | |
| `asinh`, `acosh`, `atanh` | ✅ | ✅ | `asinh` fixed for large \|x\| (2\|e\| extra bits) |
| `erf`, `erfc` | ✅ | ✅ | Series + complementary asymptotic; MPFR 1-ULP on \|x\|≲4 |
| `gamma`, `ln_gamma` | ✅ | ✅ | Stirling + reflection; factorial integers; MPFR 1-ULP |
| `digamma`, `gammainc` | ✅ | ✅ | Recurrence + Bernoulli; lower series; Accumath ψ / γ(s,x) golds |
| `bessel_j` (integer n) | ✅ | ✅ | Power series; `n ≤ 1024`; MPFR `jn` for n=0,1,2 |
| `bessel_j_nu`, `bessel_y`, `bessel_i`, `bessel_k` | ✅ | ✅ | Real order; \(K\): \(x>0\), \(\lvertν\rvert\le 32\) |
| `elliptic_k` / `e` / `f` / `pi` | ✅ | ✅ | Carlson; \(m=k^2\), \(x=\sin\varphi\); Accumath identities |
| `legendre_p`, `assoc_legendre_p` | ✅ | ✅ | \(n\le 48\); Condon–Shortley |
| `hypergeom_2f1`, `betainc` | ✅ | ✅ | Series / Gauss / Pfaff; regularized \(I_x\) |
| `ei`, `si`, `ci`, `li`, `fresnel_s`, `fresnel_c` | ✅ | ✅ | Series + full \(f,g\) / factorial asymptotic; Accumath eval golds; MPFR oracles still open |
| `sin_cos`, `sinh_cosh` | ✅ | ✅ | Tuple methods; `expr!` uses `sin`/`cos` and `sinh`/`cosh` |
| Complex: `ExactComplex` | `cexpr!` | ✅ | Per-part cancel; principal cuts; no `atan2`/`rem_pi`; MPFR add/mul; elliptic Carlson golds |
| Complex elliptic \(K,E,\Pi\) | `elliptic_k` … `elliptic_pi_inc` | ✅ | Carlson in \(\mathbb{C}\); \(K(1)=+\infty\); Legendre + cut golds |
| Complex \({}_2F_1\) | `hypergeom_2f1` | ✅ | Series / Euler / Pfaff / Kummer; cut on \([1,+\infty)\) |
| Constants: π, e, ln 2, ln 10, √2, φ, γ (`Consts`) | `pi`, `e`, `ln_2`, `ln_10`, `sqrt2`, `phi`, `euler_gamma` | ✅ | Progressive cache |

### 1.6 I/O and integration

| Item | Status |
| ------ | -------- |
| Parse/format: binary, octal, decimal, hexadecimal | ✅ |
| `convert_from_radix` / `convert_to_radix` | ✅ |
| `Display`, `LowerExp`, `UpperExp`, `Binary`, `Octal`, `UpperHex`, `LowerHex` (`std`) | ✅ |
| `FromStr` (decimal, `std`) | ✅ |
| `serde` (decimal string / integer, feature; requires `std`) | ✅ |
| `random_normal` (feature `random`) | ✅ |
| `expr!` / `cexpr!` compile and run tests (`trybuild` + `tests/mod.rs`) | ✅ |
| `expr!` cancellation / precision tests (root `tests/`) | ✅ |

### 1.7 Test & CI inventory (default `scripts/ci.sh`)

| Gate | Status |
| ------ | -------- |
| Forbid `f32`/`f64` identifiers in `.rs` / `.md` / `CHANGELOG` | ✅ |
| `cargo test --workspace` (debug) | ✅ (~60 lib tests pass) |
| `cargo test -p zenith-float-num --lib --release` | ✅ |
| `cargo test` with `no-default-features --features std` | ✅ |
| `cargo test --features random,serde` | ✅ |
| MPFR bit-oracle tests (`mpfr-tests`) | ✅ (release gate on Linux x86_64 in `scripts/ci.sh`) |
| CI wall-time budgets | ✅ (debug &lt; 10 min, MPFR &lt; 30 min; `CI_DEBUG_SECS` / `CI_MPFR_SECS`) |
| Seeded random tests | ✅ (default seed `0x5EED_CAFE_BADC_0D00`; `ZENITH_TEST_SEED` to replay) |
| Criterion / dedicated benches | ✅ (`zenith-float-num/benches/`: arithmetic, transcendentals, composite) |
| Cross-library compare (astro / dashu) | 🟡 (`zenith-float-compare/` + `scripts/compare-bench.sh`, release gate) |
| `proptest` / quickcheck | ⬜ (hand-written random loops, `TEST_ITERS = 256`) |

**Property tests** (`zenith-float-num/src/ops/tests.rs`): inverse pairs (ln↔exp, sin↔asin, log↔pow, etc.) with mathematically derived error bounds; exponent sampling capped at `TEST_EXP_BOUND = 1024` for runtime.

**Error documentation:** `doc/README.md` (ulp / series error bounds used by `expr!` and tests). **Precision growth:** `doc/PRECISION.md` (`MAX_PREC_RETRY`, exponent scaling, `expr!` bounds).

---

## 2. Comparison vs astro-float and dashu-float

zenith-float is an **evolution of the astro-float design** (`BigFloat` → `ExactNum`, same `Context` / `expr!` / `Consts` / MPFR harness). [astro-float 0.9.6](https://docs.rs/astro-float) is the direct ancestor. [dashu-float 0.6.0](https://docs.rs/dashu-float) is a **different architecture**: arbitrary base, native operators, Ziv-certified transcendentals, optional complex (`CBig`).

Local reference trees (`dashu-master/`, `astro-float-main/`) are not in this repo; compare via crates.io sources.

### 2.1 Feature matrix (high level)

| Capability | zenith-float | astro-float 0.9.6 | dashu-float 0.6.0 |
| ------------ | :------------: | :-----------------: | :-----------------: |
| Pure Rust kernel (no MPFR in lib) | ✅ | ✅ | ✅ |
| `no_std` + allocator | ✅ | ✅ | ✅ |
| `expr!` + `Context` | ✅ | ✅ | 🟡 (macros elsewhere) |
| Hardware `f32`/`f64` in library | 🚫 forbidden | `from_f32` / `from_f64` | API-edge literals |
| Elementary transcendentals (core set) | ✅ | ✅ | ✅ |
| `hypot`, `atan2`, `log1p`, `expm1` | ✅ | ⬜ not in astro | ✅ |
| `exp2`, `exp10`, `rem_pi` | ✅ | ⬜ not in astro | 🟡 partial |
| Native `+ − × ÷` operators | ✅ | ⬜ | ✅ |
| `fbig!` / compile-time float literals | ✅ | ⬜ | ✅ |
| Parse/format bases 2–36 | ✅ | ⬜ (bin/oct/dec/hex) | ✅ |
| Arbitrary-base float type | ✅ | ⬜ | ✅ |
| Complex (`CBig`) | ✅ `ExactComplex` | ⬜ | ✅ |
| `fma` / `mul_add` | ✅ | ⬜ | 🟡 |
| General `nth_root(n)` | ✅ | ✅ | ✅ |
| `sin_cos` / `sinh_cosh` paired APIs | ✅ | ✅ | ✅ |
| Special functions (erf, Γ, Bessel, …) | ✅ erf, Γ, J_n | ⬜ | 🟡 / separate |
| MPFR golden tests in repo | ✅ (optional) | ✅ (optional) | fuzz + unit (project policy) |
| Fuzz MPFR bit-exact (all round modes) | ✅ `tests/mpfr/fuzz_round_modes.rs` | ⬜ | ✅ |
| Progressive constant cache | ✅ `ConstCache` / `CachedFBig` | 🟡 | ✅ `ConstCache` / `CachedFBig` |
| Ziv + Ball correct-rounding proof | ✅ `ziv_round` / `Ball` | ⬜ | ✅ on transcendentals |
| `serde` / `random` | ✅ optional | ✅ default-on | ✅ optional |
| Stack-inlined small values | ✅ `INLINE_WORDS` | ⬜ | ✅ |
| Published crate + bench history | ✅ 0.1.0 + TSV history | ✅ 0.9.x | ✅ |

### 2.2 zenith-float vs astro-float (lineage)

zenith-float is already a **superset** of astro-float’s math API. astro-float does **not** exceed zenith on transcendentals.

| zenith has; astro 0.9.6 lacks | Notes |
| ------------------------------ | ------- |
| `exp2`, `exp10` | zenith delegates to `pow` |
| `log1p`, `expm1` | |
| `atan2`, `hypot` | |
| `rem_pi` (public) | astro reduces internally only |
| No `from_f32` / `from_f64` | zenith policy |
| `random` / `serde` off by default | astro enables by default |

Shared vs dashu: public `ziv_round` / `Ball` plus the same retry loop already used by transcendentals; binary mantissa (not arbitrary base internally except `RadixFloat`).

### 2.3 Where zenith-float must go beyond both

Correctness and operability at mixed-precision, many-op expressions matter more than matching every dashu feature on day one.

| Priority | Status | Gap | Why it matters |
| ---------- | :------: | ----- | ---------------- |
| P0 | ✅ | MPFR oracle in release CI (or nightly) | `scripts/ci.sh` runs `--features mpfr-tests --release` on Linux x86_64 |
| P0 | ✅ | MPFR coverage for `exp2`, `exp10`, `rem_pi` | `compare_ops` / `compare_special` (`rem_pi`: identity for \|x\| small; sin/cos vs MPFR for large) |
| P0 | ✅ | Working-precision caps documented + enforced | `MAX_PREC_RETRY` / `bump_prec_retry` + `doc/PRECISION.md` |
| P1 | ✅ | `fma` | `ExactNum::fma` / `mul_add` + `expr!`; full-width product; skip when exponents are disjoint |
| P1 | ✅ | `nth_root(n)` | `ExactNum::nth_root` + `expr!` `root`; MPFR 1-ULP for n=4,5 in fuzz harness |
| P1 | ✅ | `sinh_cosh` paired evaluation | `ExactNum::sinh_cosh` + MPFR paired compare |
| P1 | ✅ | Fuzz MPFR differential (all rounding modes) | `tests/mpfr/fuzz_round_modes.rs` |
| P1 | 🟡 | Benchmark regression tracking in CI (optional nightly) | Opt-in: `CI_BENCH=1 ./scripts/ci.sh` → `scripts/bench-compare.sh` |
| P1 | 🟡 | Release compare vs astro-float / dashu-float | `doc/compare-results.tsv` (astro 132-bit); `./scripts/compare-bench.sh --dashu` optional; not default CI |
| P1 | ✅ | `expr!` correct-rounding semantics documented per op | `doc/EXPR.md` (all leaves, including constants and `ldexp`/`logb`) |
| P2 | ✅ | Extra constants (φ, √2, γ) | `ConstCache` √2, φ, `euler_gamma`; `expr!` `sqrt2`/`phi`/`euler_gamma` |
| P2 | ✅ | `LowerExp` / `UpperExp` formatting | `{:e}` / `{:E}` plus `LowerHex` (`std`) |
| P2 | ✅ | `frexp` / `ldexp` / `scalb` / `logb` | Methods + `expr!` `ldexp`/`scalb`/`logb` (`frexp`/`ilogb` are methods) |
| P3 | ✅ | Special functions (erf, Gamma, Bessel J_n) | `erf`/`erfc`, `gamma`/`ln_gamma`, integer-order `bessel_j`; MPFR 1-ULP |
| P3 | ✅ | Parallel evaluation / thread-safe shared `Consts` | `SharedConsts` (`std`, mutex around `Consts`) |

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
- [x] Add `exp2`, `exp10`, `rem_pi` to MPFR compare suite
- [x] MPFR gate in CI (nightly or release-only job; Linux x86_64)
- [x] Document / bound internal precision growth vs exponent (avoid runaway cost)
- [x] Fuzz or structured stress harness for parse/format round-trip

### 3.2 API completeness vs dashu/astro (P1)

**Definition of done (each item):** public `ExactNum` method, `expr!` leaf when the op is macro-expressible, MPFR bit-oracle under every `RoundingMode`, and §2.1 matrix updated to ✅.

**Order:** implement top-to-bottom; later items may call earlier ones.

#### 3.2.1 Arithmetic & roots

- [x] **`fma` / `mul_add`** — `a*b + c` rounded once at precision `p` (no intermediate round of `a*b`). Blocks compensated summation and stable Horner evaluation.
  - [x] `ExactNum::fma` / `mul_add` + `expr!` `fma` / `mul_add` (`mul_full_prec` + `add_full_prec` + retry)
  - [x] Mantissa: full-width product for overlapping magnitudes (required for MPFR-correct sticky bits); skip full product when exponents differ by more than `p + 2` words
  - [x] MPFR oracle (`mpfr_fma`) 1-ULP for all rounding modes
- [x] **`nth_root(n)`** — generalize `sqrt` / `cbrt` (`n = 2` and `n = 3` delegate; composite factors via sqrt/cbrt; prime roots via Newton).
  - [x] `n ≥ 2`; `n = 0` → error; even `n` rejects negative operands
  - [x] `ExactNum::nth_root` + `expr!` `root(x, n)`
  - [x] MPFR oracle for general `nth_root` (n=4,5 at 1-ULP in `fuzz_round_modes`; n=2/3 via sqrt/cbrt; unit-tested)

#### 3.2.2 Transcendentals (paired evaluation)

- [x] **`sinh_cosh`** — single `exp(|x|)` path; each output rounded independently at `p`.
  - [x] `ExactNum::sinh_cosh(x, p, rm, cc) -> (ExactNum, ExactNum)`
  - [x] `ExactNum::sin_cos` paired trig (shared reduction); `sinh_cosh` public paired API
  - [x] MPFR: `sinh`/`cosh`/`sin`/`cos` each match oracle; paired APIs compared via `mpfr_sinh_cosh` / `mpfr_sin_cos`

#### 3.2.3 Sign & successor (software IEEE semantics)

- [x] **`copysign`**, **`next_after`** — integer-limb sign and total-order successor; no hardware floats.
  - [x] `copysign(magnitude, sign)`; `next_after(x, toward)` with explicit direction
  - [x] MPFR oracle for `copysign` at extreme precision (`tests/mpfr/compare_copysign_test.rs`)

### 3.3 Testing & performance (P1)

- [x] Replace `#[ignore]` `*_perf` with Criterion benches in `benches/`
- [x] TSV baseline file + `scripts/bench-compare.sh` (no JSON)
- [x] Cross-library compare harness (`zenith-float-compare`, bigfloat-bench workloads)
- [x] Capture `doc/compare-results.tsv` at release and document vs astro 0.9.6 / dashu 0.6.0
- [x] Fuzz harness: MPFR bit-exact under all rounding modes (dashu-style; complements property tests in `ops/tests.rs`)
- [x] Optional: nightly bench regression gate in CI (`CI_BENCH=1` runs `scripts/bench-compare.sh`)
- [x] Track CI wall time budget (target: full gate &lt; 10 min debug, &lt; 30 min with MPFR) — `scripts/ci.sh` (`CI_DEBUG_SECS` / `CI_MPFR_SECS`; `CI_SKIP_TIME_BUDGET=1` to skip)
- [x] Seed-controlled random tests for reproducible failures — `ZENITH_TEST_SEED` / `reseed_random`; panic hook prints the seed
- [x] Memory high-water tests for large precisions (OOM → clean `NaN`, not hang)

### 3.4 Application integration (not this crate)

zenith-float ships `ExactNum`, `expr!`, `Context`, and `SharedConsts`. Host engines own formula corpora, expression ABIs, and any purge of hardware floats in *their* trees.

- [x] Batch evaluator with shared `Consts` / `Context` reuse (`SharedConsts` for threads; `Context` still owns one `Consts`)
- [ ] 🚫 Stable ABI for a host expression IR (application crate)
- [ ] 🚫 Precision policy per formula class (application)
- [ ] 🚫 Host CI over large formula corpora (application)
- [ ] 🚫 Purge hardware floats from a host numeric path (application)

### 3.5 Out of scope for zenith-float

- [ ] 🚫 Hardware `f32`/`f64` converters inside the library
- [x] Rectangular complex arithmetic (`ExactComplex`); not a dashu `CBig` clone
- [ ] 🚫 Computer algebra (Risch, towers, etc.) — host symbolic layer
- [ ] 🚫 Replacing MPFR/GMP at test time (MPFR remains oracle only)

### 3.6 API polish (P2)

Does not block a public 0.1.x; aligns with §2.3 P2 items.

- [x] `frexp` / `ldexp` / `scalb` / `logb` / `ilogb` — IEEE-style decomposition without hardware floats
- [x] `LowerHex` + scientific `Display` options (`LowerExp` / `UpperExp` formatting)
- [x] Additional `Consts` beyond π, e, ln 2, ln 10 (φ, √2, γ)
- [x] `expr!` correct-rounding semantics documented per op (`doc/EXPR.md`)
- [x] Getting-started narrative (`doc/GETTING_STARTED.md`)
- [x] `Context::with_rounding_mode` scoped combinator
- [x] `two_sum` / `two_product`, `fused_sum` / `fused_dot`, Horner `polyval`
- [x] `ExactComplex` elementary set: `tan`, `sinh`/`cosh`/`tanh`, `sqrt`, `pow`, inverse trig and inverse hyperbolic (principal branches)

### 3.7 Load-bearing exclusions (not a backlog)

See `doc/Additions_to_existing_29-082026.md`. Macros are part of this crate (`expr!`, `cexpr!`, `exact!`, `fbig!`); do not refuse work because Accumath avoids macros.

- [ ] 🚫 Hardware IEEE converters inside the library — CI and §3.5; packing bits belongs in a caller or a separate crate (`frexp` / `ilogb`)
- [x] Complex expressions — `cexpr!` (cancellation on both parts; imaginary unit `I`; leaf set + `tests/mod.rs` / trybuild)
- [x] Bessel \(Y_ν\), \(I_ν\), \(K_ν\), non-integer \(J_ν\) — own numerics (`bessel_y` / `bessel_i` / `bessel_k` / `bessel_j_nu`); MPFR `jn` oracle on integer order

---

## 4. Release readiness gates

Before calling a version **production-ready** as a public crate:

1. **All P0 items checked** (including MPFR for every exported transcendental).
2. **`./scripts/ci.sh` green** on Linux x86_64, debug + release.
3. **MPFR suite green** on same platform (`--features mpfr-tests`).
4. **No known hang** in property tests at `TEST_EXP_BOUND` and `TEST_ITERS = 256`.
5. **CHANGELOG + README** match public API (no stale checklist).

---

## 5. Related files

| Path | Purpose |
| ------ | --------- |
| `scripts/ci.sh` | Default public CI |
| `zenith-float-num/tests/README.md` | MPFR test instructions |
| `zenith-float-num/src/ops/tests.rs` | Random inverse property tests |
| `zenith-float-num/tests/mpfr/` | MPFR bit-oracle tests |
| `doc/GETTING_STARTED.md` | First-use narrative |
| `doc/HELP.md` | User tutorial (ships with the crate) |
| `doc/LIBRARY.md` | Public API inventory |
| `doc/EXPR.md` | `expr!` per-op rounding contract |
| `doc/PRECISION.md` | Working precision vs exponent |
| `tests/mod.rs` | `expr!` integration tests |
| `CHANGELOG.md` | Released / unreleased API notes |

---

## 6. Changelog sync (unreleased)

Already implemented but not in a crates.io release:

- Newton reciprocal at three words and up
- `hypot`, `atan2`, `log1p`, `expm1`
- `exp2`, `exp10`, `rem_pi`
- CI hardening and `asinh` precision fix (2026-08-27)
- `ConstCache` / `CachedFBig`, `Ball` / `ziv_round`, stack-inlined `WordBuf`, MPFR all-round-mode fuzz
- Bench history: `doc/bench-baselines.tsv`, `doc/compare-results.tsv` (astro 0.9.x vs zenith at 132 bits)
- `euler_gamma`, `frexp`/`ldexp`/`scalb`/`logb`/`ilogb`, `LowerExp`/`UpperExp`/`LowerHex`
- Seeded test RNG (`ZENITH_TEST_SEED` / `reseed_random`), CI wall-time budgets, OOM → `NaN` tests
- [x] Getting-started narrative (`doc/GETTING_STARTED.md`); review of extra APIs (`doc/Additions_to_existing_29-082026.md`)
