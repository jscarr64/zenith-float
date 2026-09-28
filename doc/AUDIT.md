# zenith-float 1.0.5 release-candidate audit

Date: 2026-09-27 (CT), updated for the merged release candidate. Baseline: zenith-float 1.0.4 as published on crates.io. The special-function code in 1.0.2 and 1.0.4 is the same: `ops/special.rs` did not change between them, and the erf sweep gives bit-identical wrong values on both.

Paths under `audit/`, `refs/` and `logs/` below are in the maintainer workspace, not in this repository. Files under `tests/`, `zenith-float-num/`, `nostd-tests/` and `scripts/` are in the repository.

Method:

- Every public real and complex special function or math routine is evaluated by a probe (`audit/probe-rc`, built against the release candidate) at P = 128 and P = 512 bits.
- The probe output is compared with mpmath 1.4.1 at 50 digits, using the same binary inputs (each decimal rounded to P bits first).
- Pass: ≥ min(P−8, 158) correct bits. 158 bits is the ceiling of a 50-digit reference.
- Accuracy above that ceiling is covered by the 512/1024-bit cases in `tests/special_audit_mpmath.rs` and by the precision sweeps below.
- Driver: `audit/audit.py`. Results: `audit/rc-128.json`, `audit/rc-512.json`, `audit/v104-128.json`, `audit/v104-512.json`. Table generator: `audit/mk_audit_table.py`.

Totals: 130 functions, 1909 points per precision (the merged candidate adds 40 points: complex Bessel at \(z=2\), \(-500+i\), \(-300-40i\), and complex \({}_2F_1\) near \(z=1\) with integer \(c-a-b\)).

| Build | Points below target at 128 bits | Points below target at 512 bits |
|---|---|---|
| 1.0.4 | 195 (in 37 functions) | 190 (in 36 functions) |
| RC | 3 | 9 |

All 12 RC points are explained in the table: conventions (Γ(+0), `li` outside its domain), 50-digit reference artifacts (near zeros of the function, and one \({}_2F_1\) point where mpmath's own 50-digit evaluation cancels), and one slow point (known issue 1).

## Summary of status

- **fixed**: 45 functions (wrong in 1.0.4, correct in the RC, each with a regression test).
- **verified**: 85 functions (correct in both).
- **broken-not-fixed**: no wrong values at the probe points. Three known issues (one performance, two limitations) are listed at the end with workarounds; a fix is planned for each.
- **untested**: none of the public special functions. Numeric methods are covered separately (below).

## Merged release candidate: follow-up fixes

The first candidate left five open items. Items 3 and 4 are fixed here; items 1 and 2 are known issues (fix planned); item 5 is a convention. Checking the item-2 workaround (complex `bessel_j_nu` for large real arguments) found three more complex Bessel faults, fixed here as well.

- **Item 3: complex `hypergeom_2f1` near \(z=1\) with integer \(m=c-a-b\)** returned NaN.
  - Fix (`zenith-float-num/src/complex_hypergeom.rs`): for \(\lvert 1-z\rvert\le1/2\), or \(\lvert 1-z\rvert<1\) with \(\lvert z\rvert\ge1\), the logarithmic connection (DLMF 15.8.10, A&S 15.3.10–15.3.11) for \(m\ge0\), after Euler's transformation \((1-z)^mF(c-a,c-b;c;z)\) for \(m<0\). Digammas are updated recursively; cancellation is measured and the evaluation repeated with more bits. \(c-a-b\) within \(\delta\) of an integer gets \(\log_2(1/\delta)\) guard bits in the generic \(1-z\) transform.
  - Grid check (`refs/check_c2f1_log.py`): 504 cases (9 \((a,b)\) pairs × \(m\in\{-3,\ldots,3\}\) plus near-integer \(c\) × 9 values of \(z\), including \(0.999+0.001i\), \(0.9999999-10^{-7}i\), \(1+10^{-12}i\), \(1.2+0.3i\), \(1.05-0.4i\), \(1.3+0.001i\)); 0 below \(p-4\) bits at 64, 128, 256 and 512 bits (worst 1–2 bits short of \(p\)).
  - Regression test `zenith-float-num/tests/c2f1_log_mpmath.rs` (generator `refs/gen_c2f1_log_refs.py`): 54 cases at 128 and 256 bits, each to \(p-4\) bits. The first candidate returns NaN for 34 of them.
- **Item 4: 32-bit test failures** were test bugs; the library is correct.
  - `ext::tests::test_ext` asserted `is_inline()` for a 128-bit value, which is 4 limbs on 32-bit (`INLINE_WORDS` = 2). The assertion now depends on the limb width.
  - proptest `add_commutes` used \(p=32\) (one 32-bit word), below `from_i64`'s documented \(p\ge64\), and got `NaN(InvalidArgument)`. It now uses multiples of 64 bits. Arithmetic at \(p\ge64\) was checked separately on i686.
  - Full `cargo test -p zenith-float-num --release --target i686-unknown-linux-musl`: 186 passed, 0 failed (`logs/i686-rc.log`; 1.0.4 and the first candidate fail `test_ext` and `add_commutes` there, and 1.0.4 also `roots_bisect_newton_brent_illinois`).
- **Complex `bessel_j_nu` / `bessel_y` with \(\operatorname{Re} z<0\) in the Hankel regime** (\(\lvert z\rvert\ge0.35(p+112)\)) used the large-argument expansion outside \(\lvert\arg z\rvert<\pi\): `bessel_j_nu(-20000, 0)`, `bessel_j_nu(-500+i, 1)`, `bessel_y(-300+40i, 2.5)` had 0 correct bits, in 1.0.4 and in the first candidate. Fix: continue from \(-z\) (DLMF 10.11.1–10.11.2).
- **Complex `bessel_i` / `bessel_k` at real \(z>0\) with real non-integer \(\nu\)** never converged, because the exactly-zero imaginary part came out as a rounding residue that Ziv's loop cannot certify. `bessel_i(2, 0.5)`: 1.0.4 returned `PrecisionRetryExhausted` after about 3 s; the first candidate took about 310 s and returned NaN. Fix: the imaginary part is an exact 0 there.
- **Complex `bessel_k` on its cut** with \(\operatorname{Im} z=0\) returned the value from below; it now returns the value from above, like `bessel_i`, `bessel_j_nu` and mpmath.
- Bessel sweep (`refs/check_cbessel_sweep.py`; logs `logs/cbessel-sweep-first-rc.txt`, `logs/cbessel-sweep-rc.txt`): 1280 cases (J, Y, I, K × 20 values of \(z\) with \(\lvert z\rvert\) from about 270 to \(10^6\), both half-planes, on and off the axes × 8 orders × 128/256 bits). First candidate: 332 below \(p-4\) bits. RC: 24, all on the negative real axis with half-integer \(\nu\) (known issue 3).
- Regression test `zenith-float-num/tests/complex_bessel_mpmath.rs` (generator `refs/gen_cbessel_refs.py`): 218 cases at 128 and 256 bits, each to \(p-4\) bits norm-wise.
- `no_std`: 4 new cases (two complex \({}_2F_1\) log cases, complex \(J\) at \(-500+i\), complex \(I\) at \(2+0i\)); the harness now has 77.

## erf / erfc / normal_cdf tail (reported 2026-09-27 by another worker)

- **Bug.** For every \(\lvert x\rvert\ge4\), `erf_at` used the asymptotic `erfc` series (switch on binary exponent ≥ 3). That series was summed past its smallest term, and at \(x\approx4\) the smallest term is only about \(e^{-16}\).
- **Confirmed on 1.0.2 and 1.0.4** (identical output), at 256 bits:

  | Call | 1.0.2 / 1.0.4 | True value |
  |---|---|---|
  | `erfc(4.2)` | `1.97610e349` | `2.8554941795921842e-9` |
  | `erfc(10)` | `4.99346e17` | `2.0884875837625448e-45` |
  | `erfc(12.5)` | `7.69497e-56` | `6.2319427819799110e-70` (1.2e14 times too large) |
  | `erf(4)` | `-4.39943e366` | |
  | `erf(-4.2)` | `1.97610e349` | |
  | `normal_cdf(-6)` | `1.52942e433` | `9.8658764503769814e-10` |

- **RC fix** (`zenith-float-num/src/ops/special.rs`, `erf_x2_bits` / `erf_at` / `erfc_at` / `erfc_asymptotic`):
  - asymptotic series only when \(x^2\ge0.7p+4\), stopping at its smallest term;
  - otherwise the Maclaurin series with \(1.5x^2\) guard bits, plus another \(1.5x^2\) for the \(1-\operatorname{erf}\) cancellation in `erfc`;
  - `normal_cdf` via `erfc(-z/√2)/2`.
- **Sweep** (`refs/check_erf_grid.py`): 561 points per function. These are x in [3.5, 30] step 0.1, the reported points, their negatives, and [−3.5, 3.5] step 0.25, for `erf`, `erfc` and `normal_cdf`.

  | Build | Precisions | Points below P−8 bits |
  |---|---|---|
  | RC | 64, 128, 192, 256, 512, 1024 | 0 (worst: 2 bits short of P) |
  | 1.0.4 | 128 | 619 |
  | 1.0.4 | 64 + 256 | 1295 |

- **Regression test**: `zenith-float-num/tests/erf_tail_mpmath.rs` (generator `refs/gen_erf_tail_refs.py`).
  - 221 cases with 50-digit mpmath references: `erf`/`erfc` at 64/128/256 bits for x in {3.5 … 30} (23 values) and negative x (11 values), `erfc` at 1024 bits, and `normal_cdf` lower tail.
  - Plus a symmetry test.
  - Passes on the RC. On 1.0.4, 109 cases fail and the symmetry test fails.
  - The `no_std` harness also carries 8 erf-tail cases.

## Bernoulli numbers (found while making the no_std tests pass)

- **Bug.** `bernoulli_even` (real) and `even_bernoulli_numbers` (complex) used the floating-point Akiyama–Tanigawa recurrence, which cancels catastrophically.
- **64-bit limbs.** Internal Γ at 192 bits had 135 correct bits, and at 256 bits had 197. The outer Ziv loop hid this, and public results passed.
- **32-bit limbs** (`thumbv7em-none-eabihf` under QEMU). Internal Γ at 192 bits had 101 bits, and public `bessel_i(2.25, 1.5)` / `(2.25, 0.5)` / `(2.25, -0.5)` at 128 bits had 101 bits.
- **Fix.** `ops::special::even_bernoulli`: exact tangent numbers (Brent–Harvey), then \(B_{2k}=(-1)^{k-1}2kT_k/(4^k(4^k-1))\) with one rounding. It is used by real and complex Γ / lnΓ / ψ. It is also cheaper: one O(n²) pass per evaluation instead of O(k²) per term.
- **Evidence.**
  - QEMU and host `no_std` runs: 77/77 pass.
  - `refs/check_gamma_sweep.py`: gamma / ln_gamma / digamma at 64–1024 bits, 8 arguments, all within 2 bits. 1.0.4 misses 12 (precision, function) pairs because of the 274/272/526-bit caps.
  - i686 (32-bit limbs): `special_audit_mpmath`, `erf_tail_mpmath` and `numeric_methods_regress` pass.

## no_std

- `zenith-float-num` and `zenith-float` build with `--no-default-features --lib --target thumbv7em-none-eabihf`.
- There is no `alloc` feature (alloc is always required), so there is no separate alloc build.
- Test harness: `nostd-tests/`, a `#![no_std]` crate (core + alloc) with 77 mpmath reference cases (`refs/gen_nostd_cases.py`, references at 2p+64 bits). It is excluded from the workspace and from the published package. It has two runners:
  - `nostd-host`: `#![no_std]` + `#![no_main]` on x86_64 Linux. It links only libc: `posix_memalign` allocator, `write` for output.
  - `nostd-qemu`: bare-metal Cortex-M4F on `thumbv7em-none-eabihf`, run under `qemu-system-arm -machine mps2-an386` with semihosting (`cortex-m-rt`, `embedded-alloc`).
- Script: `scripts/ci_nostd.sh`. Results: host 77 passed / 0 failed; QEMU 77 passed / 0 failed (`logs/ci-nostd-rc.log`).

## Numeric methods (not in the probe table)

`audit/probe/src/bin/spot.rs` checks quadrature, roots, ODE, Chebyshev and polynomial roots against known values. Output: `logs/spot-pristine.txt` (1.0.4) and `logs/spot-fixed.txt` (RC).

| Case | 1.0.4 | RC | Cause |
|---|---|---|---|
| `gauss_legendre` / `tanh_sinh` on [2,3] | `None` | 133 bits | `ExactNum::cmp` returned raw word differences, so `cmp(..) == Some(1)` was false for endpoints with the same exponent |
| `brent` with `root_default_tol` | `None` | 130 bits | Stalled; rewritten as the classic Brent method |
| `rk45_adaptive` with rtol 1e-12 | `None` | 42 bits | `cmp` bug |
| Others | pass | pass | |

Test: `tests/numeric_methods_regress.rs` (5 tests; all fail on 1.0.4).

## Other checks

- **Docs vs behaviour.**
  - `README.md`, `CHANGELOG.md`, `doc/HELP.md`, `doc/LIBRARY.md`, `doc/REPRODUCIBILITY.md` and `doc/ZENITH_FLOAT_CAPABILITIES.md` updated to 1.0.5, including the known issues below.
  - The capabilities `erf` row claimed MPFR 1-ULP only for \(\lvert x\rvert\lesssim4\), and the `normal_cdf` formula row was stale. Both are corrected.
  - `hypergeom_2f1` and Carlson cap docs corrected; the complex `hypergeom_2f1` limitation (NaN near \(z=1\) with integer \(c-a-b\)) is removed from the rustdoc and `doc/LIBRARY.md`.
  - `elliptic_f` now documents \(F(\pm1|1)=\pm\infty\).
- **MSRV.** No `rust-version` is declared and no docs claim one. Empirically the crate needs Rust ≥ 1.93.0: before 1.93, `_addcarry_u64` / `_subborrow_u64` require `unsafe`.
- **32-bit (i686-unknown-linux-musl) full `cargo test -p zenith-float-num --release`.** 186 passed, 0 failed on the merged candidate. The two test-only failures present on 1.0.4 and the first candidate (`ext::tests::test_ext`, proptest `add_commutes`) are fixed in the tests (item 4 above); 1.0.4 also fails `roots_bisect_newton_brent_illinois` there, fixed by the `cmp` normalization.
- **Release checks on the merged candidate** (2026-09-27, CT): `cargo test --workspace` 216 passed; `cargo test -p zenith-float --all-features` 25 and `cargo test -p zenith-float-num --release --all-features` 214 (including 15 MPFR differential tests), 239 in total, 0 failed; `cargo clippy --workspace --all-targets -- -D warnings` with default and with all features clean; `cargo fmt --all --check` clean; `scripts/ci_nostd.sh` host and QEMU 77/77; `cargo publish --workspace --dry-run` succeeds. Logs: `logs/checks-rc-final.log`, `logs/ci-nostd-rc.log`, `logs/i686-rc.log`, `logs/publish-dry-run-rc.log`, `logs/package-list-rc.log`.
- **Package contents.** `doc/` ships in the `zenith-float` crate (only `doc/bench-baselines.tsv` and `doc/compare-results.tsv` are excluded), so this file is in the package; `nostd-tests/` and `scripts/` are excluded.

## Per-function table

The evidence column gives, for the RC, the number of points and the worst correct bits at 128 / 512 bits (512-bit numbers are capped by the 50-digit reference), and the number of 1.0.4 points below target.

| Function | Status | Evidence (RC, correct bits; 1.0.4 failing points) |
|---|---|---|
| `sqrt` | verified | 15 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `cbrt` | verified | 14 pts; worst 128/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `reciprocal` | verified | 12 pts; worst 128/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `nth_root` | fixed | 20 pts; worst 127/128, exact/512; 1.0.4 bad: 2@128, 2@512 |
| `powi` | verified | 20 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `powsi` | verified | 12 pts; worst 128/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `hypot` | verified | 5 pts; worst 129/128, 168/512; 1.0.4 bad: 0@128, 0@512 |
| `fma` | verified | 3 pts; worst 130/128, 170/512; 1.0.4 bad: 0@128, 0@512 |
| `ln` | verified | 16 pts; worst 128/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `log2` | verified | 13 pts; worst 128/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `log10` | verified | 13 pts; worst 130/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `log1p` | verified | 10 pts; worst 128/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `exp` | verified | 11 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `exp2` | verified | 5 pts; worst 129/128, 168/512; 1.0.4 bad: 0@128, 0@512 |
| `exp10` | verified | 5 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `expm1` | verified | 8 pts; worst 128/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `sin` | verified | 12 pts; worst 128/128, 168/512; 1.0.4 bad: 0@128, 0@512 |
| `cos` | verified | 12 pts; worst 128/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `tan` | verified | 12 pts; worst 127/128, 168/512; 1.0.4 bad: 0@128, 0@512 |
| `asin` | verified | 12 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `acos` | verified | 12 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `atan` | verified | 13 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `sinh` | verified | 8 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `cosh` | verified | 7 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `tanh` | verified | 8 pts; worst 129/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `asinh` | verified | 12 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `acosh` | verified | 6 pts; worst 128/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `atanh` | verified | 8 pts; worst 130/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `rem_pi` | verified | 6 pts; worst 169/128, 169/512; 1.0.4 bad: 0@128, 0@512 |
| `pow` | verified | 8 pts; worst 128/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `log` | fixed | 4 pts; worst 127/128, exact/512; 1.0.4 bad: 1@128, 1@512 |
| `atan2` | verified | 8 pts; worst 129/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `erf` | fixed | 19 pts; worst 128/128, exact/512; 1.0.4 bad: 5@128, 8@512 |
| `erfc` | fixed | 21 pts; worst 127/128, exact/512; 1.0.4 bad: 9@128, 11@512 |
| `gamma` | fixed | 17 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512; `0`: Γ(+0) = +∞ (C99 `tgamma` convention; mpmath raises); high precision: Γ capped at ~274 bits in 1.0.4 (gamma(11.5)@1024 had 276); sweep 64-1024 bits all within 2 bits |
| `ln_gamma` | fixed | 9 pts; worst 128/128, exact/512; 1.0.4 bad: 2@128, 2@512; high precision: capped at ~272 bits in 1.0.4; sweep 64-1024 bits all within 2 bits |
| `digamma` | fixed | 11 pts; worst 122/128, 169/512; 1.0.4 bad: 0@128, 0@512; `1.4616321449683623` at 512: reference-limited (50-digit mpmath near a zero); 2000-bit reference gives ≥511 bits; high precision: capped at ~526 bits in 1.0.4; sweep 64-1024 bits all within 2 bits |
| `gammainc` | verified | 10 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `gammainc_upper` | verified | 10 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `betainc` | fixed | 12 pts; worst 127/128, exact/512; 1.0.4 bad: 4@128, 4@512 |
| `hyp2f1` | fixed | 15 pts; worst 127/128, exact/512; 1.0.4 bad: 10@128, 10@512 |
| `ei` | fixed | 15 pts; worst 127/128, exact/512; 1.0.4 bad: 3@128, 1@512; `0.3725074107813666` at 512: reference-limited (50-digit mpmath near a zero); 2000-bit reference gives ≥511 bits; high precision: ei(-1)@512 limited by Euler γ (383 bits) in 1.0.4 |
| `ci` | fixed | 8 pts; worst 128/128, exact/512; 1.0.4 bad: 0@128, 0@512; `0.6165054856207162` at 512: reference-limited (50-digit mpmath near a zero); 2000-bit reference gives ≥511 bits; high precision: ci(0.5)@512 limited by Euler γ in 1.0.4 |
| `si` | verified | 10 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `li` | fixed | 10 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512; `0`: outside documented domain (x > 0, x ≠ 1): NaN(InvalidArgument); `1.4513692348833810` at 512: reference-limited (50-digit mpmath near a zero); 2000-bit reference gives ≥511 bits; `1`: outside documented domain (x > 0, x ≠ 1): NaN(InvalidArgument); high precision: li(10)@512 limited by Euler γ in 1.0.4 |
| `fresnel_s` | fixed | 14 pts; worst 128/128, exact/512; 1.0.4 bad: 3@128, 0@512 |
| `fresnel_c` | fixed | 14 pts; worst 128/128, 162/512; 1.0.4 bad: 3@128, 0@512 |
| `ai` | fixed | 15 pts; worst 128/128, exact/512; 1.0.4 bad: 3@128, 3@512; high precision: ai(-30)@512 had 320 bits in 1.0.4 |
| `bi` | fixed | 15 pts; worst 127/128, exact/512; 1.0.4 bad: 3@128, 3@512 |
| `ai_prime` | fixed | 15 pts; worst 127/128, exact/512; 1.0.4 bad: 9@128, 9@512 |
| `bi_prime` | fixed | 15 pts; worst 127/128, exact/512; 1.0.4 bad: 9@128, 9@512 |
| `bessel_jn` | fixed | 24 pts; worst 127/128, 169/512; 1.0.4 bad: 1@128, 0@512 |
| `bessel_j` | fixed | 28 pts; worst 127/128, 168/512; 1.0.4 bad: 4@128, 0@512; high precision: bessel_j_nu(1e-10, 10.5)@512 had 276 bits in 1.0.4 |
| `bessel_y` | fixed | 35 pts; worst 127/128, exact/512; 1.0.4 bad: 5@128, 0@512; high precision: Euler γ cap (383 bits) at 512 in 1.0.4 |
| `bessel_i` | verified | 32 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `bessel_k` | fixed | 35 pts; worst 127/128, exact/512; 1.0.4 bad: 4@128, 8@512; high precision: bessel_k(20, 3)@512 cases in special_audit_mpmath |
| `elliptic_k` | verified | 9 pts; worst 128/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `elliptic_e_complete` | verified | 9 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `elliptic_f` | fixed | 20 pts; worst 127/128, 168/512; 1.0.4 bad: 1@128, 1@512 |
| `elliptic_e` | verified | 20 pts; worst 127/128, 168/512; 1.0.4 bad: 0@128, 0@512 |
| `elliptic_pi_complete` | fixed | 16 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512; high precision: Π(0.9|0.5)@512 had ~340 bits in 1.0.4 |
| `elliptic_pi` | fixed | 18 pts; worst 127/128, 166/512; 1.0.4 bad: 0@128, 0@512; high precision: Π(0.3,0.5,0.4)@512 had ~340 bits in 1.0.4 |
| `jacobi_sn` | verified | 30 pts; worst 127/128, 168/512; 1.0.4 bad: 0@128, 0@512 |
| `jacobi_cn` | verified | 30 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `jacobi_dn` | verified | 30 pts; worst 128/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `jacobi_cd` | fixed | 30 pts; worst 127/128, exact/512; 1.0.4 bad: 6@128, 6@512 |
| `jacobi_sd` | verified | 30 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `jacobi_nd` | verified | 30 pts; worst 126/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `jacobi_dc` | fixed | 30 pts; worst 126/128, exact/512; 1.0.4 bad: 6@128, 6@512 |
| `jacobi_nc` | verified | 30 pts; worst 126/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `jacobi_sc` | verified | 30 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `jacobi_ns` | verified | 30 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `jacobi_ds` | verified | 30 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `jacobi_cs` | verified | 30 pts; worst 128/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `jacobi_am` | verified | 12 pts; worst 127/128, 168/512; 1.0.4 bad: 0@128, 0@512 |
| `legendre_p` | verified | 20 pts; worst 128/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `assoc_legendre_p` | verified | 12 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `hermite_he` | verified | 9 pts; worst 130/128, 162/512; 1.0.4 bad: 0@128, 0@512 |
| `hermite_h` | verified | 9 pts; worst 128/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `laguerre` | verified | 9 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `gen_laguerre` | verified | 8 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `chebyshev_t` | verified | 12 pts; worst 128/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `chebyshev_u` | verified | 12 pts; worst 128/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `gegenbauer` | verified | 8 pts; worst 128/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `normal_pdf` | verified | 4 pts; worst 128/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `normal_cdf` | fixed | 5 pts; worst 128/128, 168/512; 1.0.4 bad: 2@128, 2@512 |
| `gamma_pdf` | verified | 4 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `beta_pdf` | verified | 3 pts; worst 129/128, 168/512; 1.0.4 bad: 0@128, 0@512 |
| `binomial_pmf` | verified | 3 pts; worst 128/128, 164/512; 1.0.4 bad: 0@128, 0@512 |
| `poisson_pmf` | verified | 4 pts; worst 128/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `chi_squared_cdf` | verified | 4 pts; worst 128/128, 168/512; 1.0.4 bad: 0@128, 0@512 |
| `student_t_pdf` | verified | 4 pts; worst 127/128, 159/512; 1.0.4 bad: 0@128, 0@512 |
| `c.exp` | verified | 13 pts; worst 126/128, 169/512; 1.0.4 bad: 0@128, 0@512 |
| `c.ln` | verified | 13 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `c.sqrt` | verified | 13 pts; worst 127/128, 575/512; 1.0.4 bad: 0@128, 0@512 |
| `c.sin` | verified | 13 pts; worst 126/128, 178/512; 1.0.4 bad: 0@128, 0@512 |
| `c.cos` | verified | 13 pts; worst 126/128, 178/512; 1.0.4 bad: 0@128, 0@512 |
| `c.tan` | verified | 13 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `c.sinh` | verified | 13 pts; worst 127/128, 168/512; 1.0.4 bad: 0@128, 0@512 |
| `c.cosh` | verified | 13 pts; worst 127/128, 178/512; 1.0.4 bad: 0@128, 0@512 |
| `c.tanh` | verified | 13 pts; worst 126/128, 198/512; 1.0.4 bad: 0@128, 0@512 |
| `c.asin` | verified | 13 pts; worst 128/128, 183/512; 1.0.4 bad: 0@128, 0@512 |
| `c.acos` | verified | 13 pts; worst 128/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `c.atan` | verified | 12 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `c.asinh` | verified | 13 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `c.acosh` | fixed | 13 pts; worst 128/128, exact/512; 1.0.4 bad: 1@128, 1@512 |
| `c.atanh` | fixed | 13 pts; worst 127/128, 185/512; 1.0.4 bad: 1@128, 1@512 |
| `c.log1p` | verified | 13 pts; worst 128/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `c.expm1` | verified | 13 pts; worst 127/128, 185/512; 1.0.4 bad: 0@128, 0@512 |
| `c.erf` | fixed | 13 pts; worst 127/128, 174/512; 1.0.4 bad: 0@128, 1@512 |
| `c.erfc` | fixed | 13 pts; worst 128/128, 174/512; 1.0.4 bad: 0@128, 1@512 |
| `c.gamma` | fixed | 13 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512; high precision: capped at ~275 bits in 1.0.4 (c.gamma(0.5+0.5i)@512) |
| `c.ln_gamma` | fixed | 13 pts; worst 127/128, 170/512; 1.0.4 bad: 8@128, 8@512 |
| `c.digamma` | fixed | 13 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512; high precision: capped at ~275 bits in 1.0.4 |
| `c.ai` | fixed | 13 pts; worst 127/128, exact/512; 1.0.4 bad: 2@128, 3@512 |
| `c.bi` | fixed | 13 pts; worst 128/128, 173/512; 1.0.4 bad: 7@128, 7@512 |
| `c.ei` | fixed | 13 pts; worst 127/128, exact/512; 1.0.4 bad: 5@128, 5@512 |
| `c.si` | fixed | 13 pts; worst 127/128, 176/512; 1.0.4 bad: 6@128, 6@512 |
| `c.ci` | fixed | 13 pts; worst 127/128, 179/512; 1.0.4 bad: 13@128, 13@512 |
| `c.li` | verified | 13 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `c.fresnel_s` | verified | 13 pts; worst 127/128, 166/512; 1.0.4 bad: 0@128, 0@512 |
| `c.fresnel_c` | verified | 13 pts; worst 127/128, 166/512; 1.0.4 bad: 0@128, 0@512 |
| `c.elliptic_k` | verified | 13 pts; worst 127/128, exact/512; 1.0.4 bad: 0@128, 0@512 |
| `c.elliptic_e_complete` | verified | 13 pts; worst 127/128, 176/512; 1.0.4 bad: 0@128, 0@512 |
| `c.pow` | verified | 3 pts; worst 126/128, 575/512; 1.0.4 bad: 0@128, 0@512 |
| `c.bessel_j` | fixed | 36 pts; worst 127/128, exact/512; 1.0.4 bad: 8@128, 8@512 |
| `c.bessel_y` | fixed | 36 pts; worst 127/128, exact/512; 1.0.4 bad: 8@128, 8@512 |
| `c.bessel_i` | fixed | 36 pts; worst 127/128, exact/512; 1.0.4 bad: 11@128, 11@512; high precision: c.bessel_i(20+i, 2+0.5i)@512 in special_audit_mpmath |
| `c.bessel_k` | fixed (slow point) | 36 pts; worst 127/128, exact/512; 1.0.4 bad: 17@128, 14@512; `['150,-20', '2,0.5']`: correct (513 bits) but 74 s at 512 bits; audit timeout 60 s; known issue |
| `c.hyp2f1` | fixed | 14 pts; worst 127/128, 161/512; 1.0.4 bad: 9@128, 9@512; `['-0.25,1', '2.5,-0.5', '0.25,0.5', '0.9999999,-1e-7']` at 512: reference-limited (50-digit mpmath hyp2f1 cancels); 1100-bit reference gives 513 bits |

Counts: fixed 45, verified 85

## Open items from the first candidate: resolution

1. **Complex `bessel_k` slow**: known issue 1 (not changed).
2. **Real Bessel J/Y for \(\lvert x\rvert\gtrsim10^4\)**: known issue 2 (not changed).
3. **Complex `hypergeom_2f1` near \(z=1\) with integer \(c-a-b\)**: fixed (above).
4. **32-bit test-only failures**: fixed in the tests (above).
5. **Conventions that differ from mpmath**, left as is and documented: Γ(+0) = +∞ (C99 `tgamma`); `li(x)` is `NaN(InvalidArgument)` outside \(x>0\), \(x\ne1\).

## Known issues (fix planned)

1. **Complex `bessel_k` is slow for non-integer \(\nu\) at high precision when \(\lvert z\rvert\) is just inside the series regime** (\(\lvert z\rvert<0.35(p+112)\)). `bessel_k(150−20i, 2+0.5i)` takes about 74 s at 512 bits and 13 ms at 256 bits (Hankel regime there); the result is correct to 513 bits. 1.0.4 returned a wrong value quickly. Cause: the series path carries about \(3\lvert z\rvert\) guard bits for the \(J\pm iY\) cancellation and needs both \(J_{\pm\nu}\). Planned fix: a uniform asymptotic or recurrence path. Workaround: lower precision, integer \(\nu\) where possible (`bessel_k(150−20i, 2)` at 512 bits: 0.9 s), or allow for the run time.
2. **Real `bessel_j_nu` / `bessel_y` return `InvalidArgument` for \(\lvert x\rvert\gtrsim10^4\)** (`BESSEL_GUARD_MAX`; no large-argument Hankel expansion in the real kernel). Planned fix: the real Hankel expansion. Workaround, verified: `ExactComplex::bessel_j_nu` / `bessel_y` at \(x+0i\) with \(x>0\) agree with mpmath to working precision at \(x=2\cdot10^4\) (J and Y, 128 and 256 bits) and \(x=10^6\) (J, 128 bits); take the real part. For \(x<0\) with integer order use \(J_n(-x)=(-1)^nJ_n(x)\), \(Y\) likewise via DLMF 10.11.2.
3. **Complex Bessel exactly on the negative real axis with real \(\nu\)** (found in this round; present in 1.0.4). A part that is exactly zero, or tiny next to the other part, is not resolved. `bessel_i` / `bessel_k` with half-integer \(\nu\) return `PrecisionRetryExhausted` for large \(\lvert z\rvert\) and take minutes for small \(\lvert z\rvert\) (`bessel_i(−2, 0.5)` at 128 bits: about 5 minutes, then NaN). `bessel_j_nu(−2, 0.5)` returns a real part of about \(4\cdot10^{-136}\) instead of 0, and `bessel_k(−500, 0)` a real part of \(-5\cdot10^{119}\) instead of \(K_0(500)\approx4\cdot10^{-219}\); both are within \(2^{-p}\) of the modulus. Off the axis and on the positive real axis results are correct. Planned fix: evaluate on the axis through the reflection formulas with exact phases. Workaround: evaluate at \(-z\) and apply DLMF 10.11.1 / 10.34.1 (e.g. \(I_\nu(-x+0i)=e^{i\nu\pi}I_\nu(x)\)), or use the real functions.
