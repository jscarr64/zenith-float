# Changelog

## 2.0.0 (unreleased)

Problems in one program led to a forensic audit of the others. This 2.0.0 section is that audit in full: every issue found, the fix, and the migration for each break in the public API. The workspace version is 2.0.0. Nothing here has been published.

### Added

`elliptic_f_phi`, `elliptic_e_phi`, and `elliptic_pi_phi` take the amplitude \(\varphi\). `elliptic_f`, `elliptic_e`, and `elliptic_pi` still take \(\sin\varphi\). Downstream callers passed the angle into the \(\sin\varphi\) slot: `elliptic_f(2, m)` is `NaN(InvalidArgument)`, not \(F(2\mid m)\). The angle form accepts any real \(\varphi\), reduces \(\varphi = r + k\pi\) with \(|r|\le\pi/2\) in extended precision, and adds \(2k\) times the complete integral. \(\varphi=\pm\pi/2\) is \(\pm K\), \(\pm E\), or \(\pm\Pi\). A parameter \(m>1\) whose path is not real, and a characteristic \(n=1\) at a pole, are `NaN(InvalidArgument)` rather than a silent wrong value.

`ExactNum::inexact_accumulated()` keeps the old meaning of `inexact()`: true when this value or anything it was computed from was rounded. See "`inexact` is per operation".

### Changed

The MSRV was raised to 1.99 to match Accumath and latex-rust (breaking). `rust-version` is `"1.99"` in every crate, and `rust-toolchain.toml` pins channel `1.99` with `clippy` and `rustfmt`. A crate that compiled on the previous compiler has to move to Rust 1.99 before it can depend on 2.0.0.

Behaviour changes from the audit, each with its own section below:

- `add`, `sub`, `mul`, and `div` round at exactly `p` bits instead of `p` rounded up to a word, and subnormal results are no longer rounded twice.
- Every other rounded function also rounds at exactly `p` bits: roots, `hypot`, `fma`, `pow`, parsing, the exp/log, trigonometric and hyperbolic families, and every real special function. See "Exact-`p`-bit rounding beyond the basic operations".
- An exponent overflow under `ToZero`, under `Down` for a positive result, and under `Up` for a negative result returns `±max_value(p)` instead of `±Inf`.
- `inexact()` describes the operation that produced the value, for every rounded function and for `ExactComplex`. `inexact_accumulated()` is the old flag; `ExactComplex::inexact_accumulated()` is new.
- `ldexp` and `scalb` below the exponent range round onto the subnormal grid by mode instead of returning zero.
- `sin_cos` and `sinh_cosh` round each result at exactly `p` bits with its own `inexact()`; `sinh_cosh` overflow follows the mode.
- `rem` returns subnormal remainders instead of zero, and its `inexact()` is false: the remainder is always exact.
- `add_full_prec` and `sub_full_prec` with one zero operand return the other operand exactly (with its own `inexact()`). They used to return a value with an empty mantissa flagged inexact. All other `add`/`sub` results are unchanged: the addition kernel now writes the sum of a nonzero operand with a zero mantissa itself instead of copying the operand first (the zero-aware kernel branch), which removes the three `.clone()` sites of `add`/`sub`.
- `from_raw_parts` validates its parts and returns `NaN(InvalidArgument)` for inconsistent ones.
- `EXPONENT_MAX` and `EXPONENT_MIN` are `i32::MAX` and `i32::MIN` on 32-bit targets too.
- The `RoundingMode` documentation describes `Up`, `Down`, `ToZero`, and `FromZero` as directed modes, not tie rules (documentation only).
- `ExactNum::cmp` is now `ExactNum::compare`, and `Context::clone` is removed (breaking). See "No lint exemptions".
- `Radix` is an enum (`Bin`, `Oct`, `Dec`, `Hex`, `Other`); `Radix::Bin` and the others are written the same way (breaking only for code that relied on it being a struct).

### `add_ref` and `sub_ref`: the borrowed-borrowed sum (added)

`ExactNum::add_ref(&a, &b, p, rm)` and `ExactNum::sub_ref(&a, &b, p, rm)` (and the same two on `ExactComplex`) add or subtract two operands that are both borrowed. The result is written to a fresh buffer by the same kernel as `add`/`sub` (the zero-aware branch included); no operand is copied, a subnormal one included (the kernels read it as stored and renormalize the exact result; only `div` shifts a subnormal divisor into its own working buffer, because the division kernel needs a normalized divisor), and there is no `.clone()`, no multiply-by-one and no add-zero. Value, signs of zero in every rounding mode, both inexact flags, Inf/NaN and subnormal results are those of `add`/`sub`, pinned by `add_ref_sub_ref_grid_digest` (the 45030-case `add_sub_grid_digest` digest, taken on the pre-rule-R code) and `complex_add_ref_sub_ref_grid_digest` (the complex grid digest, also pinned on the old code).

They are the sanctioned form for a site where **both operands are still needed afterwards** (a bisection midpoint `(a+b)/2` with `a` and `b` both live, a Carlson duplication step, the `(d_prev+d_last)(d_prev-d_last)` of the SVD shift). Everywhere else the first operand is owned and dead after the call and `add(self, &b)` is used. To keep them from spreading, `tools/rules_scan.sh` counts every use of `add_ref(`/`sub_ref(` in non-test library code against `tools/ref_sum_count_baseline.txt` (a ratchet like the clone count: it may only go down; CODEOWNERS-protected; `--update-baseline` lowers it). The count starts at 4 (the four calls inside the `ExactComplex` wrappers) and is raised once, in the Rule R commit, to the number of both-live sites that conversion needs.

### `mul_ref` and `div_ref`: the borrowed-borrowed product and quotient (added)

`ExactNum::mul_ref(&a, &b, p, rm)` and `ExactNum::div_ref(&a, &b, p, rm)` multiply or divide two operands that are both borrowed, under the same limits as `add_ref`/`sub_ref`: they are the sanctioned form only at a site where **both operands are still needed afterwards** (a square `x * x`, a dot-product or matrix term, a recurrence step, a quotient whose divisor is reused). The result is written to a fresh buffer by the same kernel as `mul`/`div`; no operand is copied, a subnormal one included (the kernels read it as stored and renormalize the exact result; only `div` shifts a subnormal divisor into its own working buffer, because the division kernel needs a normalized divisor), and there is no `.clone()`, no multiply-by-one and no divide-by-one. Value, signs of zero in every rounding mode, both inexact flags, Inf/NaN and subnormal results are those of `mul`/`div`, pinned by `mul_ref_div_ref_grid_digest`, which hashes to the same value as `mul_div_grid_digest`, taken on the code before `mul_ref`/`div_ref` existed (`ExactComplex::mul`/`div` are pinned the same way by `complex_mul_div_grid_digest`).

`tools/rules_scan.sh` counts every use of `mul_ref(`/`div_ref(` in non-test library code against its own baseline, `tools/mul_ref_count_baseline.txt` (CODEOWNERS-protected, only goes down, `--update-baseline` lowers it), separate from the `add_ref`/`sub_ref` count (`tools/ref_sum_count_baseline.txt`, unchanged). The count was 0 before the Rule R mul/div change and is 232 after it.

### `ExactNum::sqr`: the square of a live value (added)

`ExactNum::sqr(&self, p, rm)` returns `self * self` rounded to `p` bits. It is a plain public function, not a copy primitive: it reads `self` twice, builds a fresh result in its own buffer, and does not consume `self` (`self` is not copied, a subnormal one included: the product kernel reads it as stored) (no `.clone()`, no multiply-by-one). Value, sign of zero in every rounding mode, both inexact flags, `op_inexact`, Inf, subnormal results and the error kind of a NaN operand are those of `ExactNum::mul_ref(&x, &x, p, rm)` and of the old `x * x`, pinned by `sqr_grid_digest` (23 operands incl. NaN of both kinds, +-Inf, +-0, the subnormal band and the exponent extremes, 9 precisions, 6 modes), whose literal was recorded on the code before `sqr` existed (`square_via_mul_ref_grid_digest`), by a field-for-field differential against `mul_ref(&x, &x)`, and by an MPFR `mpfr::sqr` oracle in every rounding mode (`fuzz_mpfr_sqr_all_round_modes`).

Squares no longer use `mul_ref`: every `mul_ref(&x, &x, p, rm)` in the library, the tests, the benches and the macro/compare crates now reads `x.sqr(p, rm)`. The library had 42 such squares among its 232 `mul_ref`/`div_ref` uses (the earlier estimate of about 56 counted test-module and non-square lines), so `tools/mul_ref_count_baseline.txt` goes from 232 to 190. `sqr` is not counted by the ratchet (it copies nothing); `tools/rules_scan_selftest.sh` checks that `sqr` uses and definitions leave the `mul_ref`/`div_ref` count unchanged. `div_ref(&x, &x, ..)` is not a square and is unchanged. Migration: `x.mul(&x, p, rm)` on a value you still need is `x.sqr(p, rm)`.

Thirty-three more `mul_ref`/`div_ref` uses are gone by restructuring, with no new primitive, no clone, and no multiply-by-one: the dead operand becomes the receiver of `mul`/`div` (complex `mul_word`/`div_word` now take `self`; `into_iter` over the Bernoulli numbers in `catalog`, the monic coefficients in `poly` and the QR normalisation; `core::mem::replace` on the FFT, LU and SVD-sweep elements that are overwritten right afterwards; a hoisted `y.is_zero()` in `roots`; a rebuilt constant in `gauss_hermite` and `hex_golds`). Every result is bit-identical (all pinned grid digests are unchanged). `tools/mul_ref_count_baseline.txt` goes from 190 to 157. Three sites that the analysis listed stay `mul_ref` because commuting their operands changes which NaN error kind survives when both operands are NaN with different kinds: `ExactComplex::exp`, `ExactComplex::sqrt` and the projection step of `ExactNumArray::qr_decomp`. `tests/nan_error_kinds.rs` now pins those cases.

Test code no longer clones either. The unit-test modules, the integration tests, the root `tests/mod.rs` and the benches had 267 `.clone()`/`.to_vec()`/`.cloned()` sites (the `.clone(` counter `tools/clone_count_baseline.txt` has never counted them: it covers non-test library code only, so it stays at 513); 236 are gone, 31 `.clone()` remain (listed in the commit message; none of the 236 changes a tested value or an oracle). An operand is rebuilt from its literal (`ExactNumNumber::small_const(1)`, `from_word`, a `mk_x = || ...` closure, the `const` NaN/Inf values built per use), a value that is dead afterwards is moved (`into_parts`), a borrowed pair uses `add_ref`/`sub_ref`/`mul_ref` (test code, not counted by the ratchet), and the benches parse their operands from the fixture literals instead of copying parsed ones. Four calls that only check the typed error for a reversed or empty interval (`gauss_kronrod_interval` twice, `gauss_legendre` twice) now pass a constant integrand instead of the identity: the integrand is never evaluated on those paths.

Test code no longer panics either. `tools/rules_scan.py` now applies the `unwrap`/`expect`/`unwrap_err`/`expect_err`/`panic!`/`unreachable!`/`todo!`/`unimplemented!` rules to test code (unit-test modules, integration tests, benches, examples, the compare crate, `nostd-tests`), with a ratchet `tools/test_panic_count_baseline.txt` (CODEOWNERS-protected, only goes down). The measured start was 2220 sites; 6 remain (`macro_util.rs` gen_pair, `ieee_soft/array.rs` perm_rows/eigen_diag: each sits on a line that also carries a `.clone()`, so any edit re-flags the clone; they go with the Option A conversion). Tests return `Result` and use `?`; `assert!`/`assert_eq!`/`matches!` carry the checks, and no oracle value or assertion was dropped. The compare crate reports `BenchError` instead of panicking on an unknown task or a failed constants cache; the benches skip a case whose setup cannot be built instead of panicking.

Layer 0a of the Option A conversion (internal; no public signature changed, every pinned grid digest identical): the `sin`/`cos`/`tan`/`atan` argument reduction and restoration take their value by `self`, the high-precision series scheme moves each cached power instead of copying it, `atan`, `fma` and `hypot` move operands that are dead afterwards, `pow`/`sin` build their small integer constants from literals (`ExactNumNumber::small_const`) instead of copying the shared constants, and decimal parsing rebuilds its operand from the digits on a precision retry instead of keeping a copy (rule R2). The `.clone(` ratchet `tools/clone_count_baseline.txt` went from 541 to 522; the helper-copy ratchet stays at 319 and no ref-form baseline moved. New pinned digests (`trig_hyperbolic`, `reciprocal_nth_root`, `exp_pow`, `series_high_precision`, `fma_stand_in_branches`, `hypot_scaling`, `decimal_parse_grid`) were recorded on the old code before each change.

Found by the layer-0a pin tests: `jacobi_am(u, 0)` (and `jacobi_am` through the `m = 0` shortcut generally) returned `u` at its own width instead of the working precision, so the Ziv roundability check `try_set_precision(p, rm, p_wrk)` subtracted past zero ("attempt to subtract with overflow" in debug builds, a wrapped bit count in release) for an operand with fewer mantissa bits than the working precision, e.g. `u = -0.6666…` (192 bits, inexact) at `p = 64`. The shortcut now carries `u` at the working precision like `sn`, `cn` and `dn`. `jacobi_all_functions_digest` is red on the old code and pins all thirteen functions on the fixed code.

`ExactNum::sqrt(-Inf)` returned NaN with no error kind, while `sqrt(-3)` returned NaN with `InvalidArgument`; both are the square root of a negative number, so `-Inf` now gives `InvalidArgument` too (`sqrt_of_negative_inf_has_the_same_nan_kind_as_a_negative_finite`, red before the fix; `+Inf` is unchanged).

### No cloning in the library (breaking, in progress)

The library no longer copies a value to hand it back. `tools/no_clone.sh` fails on `.clone()`, `.to_vec()`, `.cloned()` and on `Clone` derives and impls on big-number types, with an empty allowlist; `tools/rules_scan.sh` checks every diff. The work goes bottom-up (core, conversions, elementary, special, complex and containers, macro); each layer lists its API changes here with a migration. Counts of `.clone(` in library code: 620 at the start, 594 after layer 1b, 588 after layer 2a, 542 after the compare crate, 532 after `neg`, 530 after `frexp`, 528 after the `abs` restructure, 525 after `ldexp`/`scalb`, 522 after the zero-aware `add`/`sub` kernel branch, 520 after `int`/`fract`/`ceil`/`floor`, 519 after `round`, 514 after Rule R (`add`/`sub` take `self`).

- `ExactNum::ldexp` and `ExactNum::scalb` take `self` (`pub fn ldexp(self, n, p, rm) -> Self`); the exponent is changed in the operand's own storage. Results are bit-identical (value, sign, exponent, precision and both inexact flags; only whether the buffer is inline or on the heap may differ), pinned by a 14112-case grid digest test that gives the same digest before and after. Migration: `x.ldexp(n, p, rm)` on a value you still need afterwards becomes a call on a value you rebuild (for a constant, `ExactNum::from_u8(1, p).ldexp(n, p, rm)` instead of `one.ldexp(..)`), or reorder so the scaling is the last use of `x`; `h.sqrt(p, rm).ldexp(1, p, rm)` is recomputed where the unscaled value is read again. `expr!`'s `ldexp(..)`/`scalb(..)` are unchanged. `ExactComplex::ldexp`/`scalb` still take `&self`. Accumath and latex-rust: no `ldexp`/`scalb` call outside the bundled zenith-float copy.
- `ExactNum::frexp` takes `self` and returns the mantissa built in the operand's own storage (`(Self, Exponent)`; zero, Inf and NaN come back as themselves with exponent 0). Migration: `x.frexp()` on a value you still need becomes a call on a value you rebuild or are done with (`let e = x.exponent();` when only the exponent is wanted). Accumath: no use of `ExactNum::frexp`.
- `ExactNum::int`, `fract`, `ceil` and `floor` take `self` and work in the operand's own storage (the inner `ExactNumNumber` versions too); `ExactNumArray::int`, `fract`, `ceil` and `floor` take `self` and consume the elements. Results are bit-identical (the existing characterization tests and the oracle suites are unchanged). Migration: `x.int()` / `x.fract()` on a value you still need afterwards becomes a call on a value you rebuild or are done with; `x.fract().is_zero()` is `x.is_int()`; a digit loop that compared against `x.int()` compares against `x` itself (an integer is at most `x` exactly when it is at most `int(x)` for `x >= 0`). Inside the library: `ExactInt::from_exact_num` reads the lowest limb from `x.rem(2^w)` and the rest from the truncated quotient instead of copying `x.int()`; `exp`'s fractional part is computed as `|x| - int(x)` from the integer part already in hand instead of `x.fract()` (exact, same bits). Accumath and latex-rust: no use of `ExactNum::{int, fract, ceil, floor}`.
- `ExactNum::round` and `ExactNumArray::round` take `self` and round in the operand's own storage (the inner `ExactNumNumber::round` too). Results are bit-identical, pinned by a 726-case grid digest test (21 operands x 11 positions x 6 modes) taken on the old code. Migration: `x.round(n, rm)` on a value you still need afterwards becomes a call on a value you rebuild, or reorder so the rounding is the last use of `x`. Accumath and latex-rust: no use of `ExactNum::round`.
- `ExactNum::abs` takes `self` (`pub fn abs(self) -> Self`) and flips the sign in place only when it is negative; NaN keeps its error, `-0` becomes `+0`, `-Inf` becomes `+Inf`. `ExactNumArray::abs` takes `self` as well and updates every element in place (it was `abs(&self)`). Migration: `x.abs()` on a value you still need afterwards becomes `x.abs_cmp(&y)` when only a magnitude comparison was wanted, `x.exponent()` / `x.is_zero()` when only the exponent or zero-ness was read (neither depends on the sign), or a branch on the sign that uses the borrowed `x` directly when `x` is not negative and `-1 * x` (`ExactNum::from_i8(-1, p).mul(&x, p, rm)`, exact at the precision of `x`) when it is: `|x| + y` is `if x.is_negative() { y.sub(&x, p, rm) } else { x.add(&y, p, rm) }`, `|x| - y` is `if x.is_negative() { (-1·x).sub(&y, p, rm) } else { x.sub(&y, p, rm) }`. Where the value is mutated afterwards, pass ownership of `x` into the callee (as `gauss_kronrod`'s error estimate does). There is no `±1 * x` form for the non-negative case: nothing in the library multiplies by one to copy. Accumath: `soft.rs:185` (`self.inner.as_ref().abs()`-style call) needs the sign branch above; no other `ExactNum::abs` use outside the bundled zenith-float copy.
- `add` and `sub` take `self` and borrow the second operand (Rule R): `ExactNum::{add, sub}(self, &other, p, rm)`, `ExactComplex::{add, sub}`, `ExactNumArray::add_scalar(self, &s)` and the operators `ExactNumArray + &ExactNumArray` / `- &ExactNumArray` (`Output = Option<ExactNumArray>`, elementwise in the operand's own storage; they replace the inherent `add`/`sub`, which clippy's `should_implement_trait` rejects for a by-value `self`: call `a + &b`, or `Add::add(a, &b)`), `Ball::add`, `ComplexBall::add` and `ExactNum::two_sum`. `Ball::{ln, sqrt, bessel_j1}` and `ComplexBall::{exp, sin, cos}` take `self` as well, so the endpoints are computed once from the ball they consume. The operators `ExactNum + &ExactNum`, `ExactNum + ExactNum`, `&ExactNum + ExactNum` and the matching `-` forms exist; `&a + &b` and `&a - &b` are gone (there is no way to add two borrowed values without a copy except `add_ref`/`sub_ref`). Results are bit-identical to the borrowed forms they replace (value, sign of zero in every rounding mode, flags, Inf/NaN, subnormals), pinned by `add_sub_grid_digest` (45030 cases), `complex_add_sub_grid_digest`, `operator_grid_digest`, `add_ref_sub_ref_grid_digest` and the unchanged `ldexp`/`round` grid digests: every digest was recorded on the pre-Rule-R code and is unchanged. Migration: `a.add(&b, p, rm)` is unchanged when `a` is dead afterwards; when `a` is still needed, use `ExactNum::add_ref(&a, &b, p, rm)` if `b` is too (the sanctioned both-live form, counted by the ratchet) or rebuild `a`; `a - b` with both live is `ExactNum::sub_ref(&a, &b, p, rm)`, never `(-b) + a`. The library has 147 `add_ref`/`sub_ref` uses at both-live sites (bisection and midpoint sums, the 13 former `add_full_prec` sites of `heun_fox`/`dist`, `(a-b)` terms in the special functions, `jacobi_p_word`'s `alpha + beta`, `l1_below_series_bound`'s `re + im`); `tools/ref_sum_count_baseline.txt` pins the count. `neg_via_mul` is not used for subtraction (its count is 73, below the 74 it had before Rule R). Library sites where the owned operand was the second one now add in the other order (`a + b` is computed as `b + a`; about 16 sites, e.g. in `hypergeom`, `heun_fox`, `fft_dit`, `ExactComplex::exp`): the value, the sign of zero, and the flags are those of the original, and only the error code of a result whose two operands are NaN with *different* error kinds can differ (the first operand's error used to win). Pinned by `tests/nan_error_kinds.rs`: the operators, `add`/`sub`, `add_ref`/`sub_ref` and the complex forms are unchanged; the one public function that surfaces a different kind is `ExactComplex::asinh` (its imaginary part now carries the kind of the real part for NaN parts of different kinds); the fft butterfly is back to the original operand order. The operators `&a + b` and `&a - b` keep the old rule for two NaNs explicitly. Accumath: `ExactNum::add`/`sub` call sites that pass a still-live first operand must switch to `add_ref`/`sub_ref` (or clone, outside the library).
- `ExactNum::mul` and `ExactNum::div` take `self` and borrow the second operand (Rule R): `ExactNum::{mul, div}(self, &other, p, rm)`. The operators `ExactNum * &ExactNum`, `ExactNum * ExactNum`, `&ExactNum * ExactNum` and the matching `/` forms exist; `&a * &b` and `&a / &b` are gone (use `mul_ref`/`div_ref` when both operands are still needed). `ExactComplex`, `Ball`, `ExactNumArray` and the other types keep their borrowed `mul`/`div`. Results are bit-identical to the borrowed forms they replace (value, sign of zero in every rounding mode, both inexact flags, Inf/NaN, subnormals): `mul_div_grid_digest`, `complex_mul_div_grid_digest`, `operator_mul_div_grid_digest` and `mul_ref_div_ref_grid_digest` were recorded on the pre-change code and are unchanged, as are the add/sub, `ldexp` and `round` digests. Library sites where the second operand of a `mul` was the owned one now multiply in the other order (`a * b` computed as `b * a`); `div` is never reordered. **Error kind:** a product whose two operands are NaN with different error kinds can carry the other operand's kind at those reordered sites; the operators and `mul`/`div` themselves keep the first operand's kind. No public function surfaces a difference beyond the one already stated for `ExactComplex::asinh`: its imaginary part carries the kind of the real part when the parts are NaN of different kinds. This is accepted and pinned by `tests/nan_error_kinds.rs`. The library has 232 `mul_ref`/`div_ref` uses at both-live sites (42 of them squares `x * x`, converted to `ExactNum::sqr` in the next change, which lowers the baseline to 190); `tools/mul_ref_count_baseline.txt` pins it (raised once, from 0, in this change; it only goes down from here). `tools/ref_sum_count_baseline.txt` is unchanged at 147. Migration: `a.mul(&b, p, rm)` is unchanged when `a` is dead afterwards; for a still-live `a`, use `ExactNum::mul_ref(&a, &b, p, rm)` if `b` is too, or rebuild `a`. `expr!` generated code uses the owned form (its operands are temporaries). Accumath: `ExactNum::mul`/`div` call sites that pass a still-live first operand must switch to `mul_ref`/`div_ref` (or clone, outside the library).
- Negation takes ownership and flips the sign in place: the inherent `ExactNum::neg(&self)` and `ExactNumArray::neg(&self)` are removed, and the `std::ops::Neg` impls for `ExactNum` (existing) and `ExactNumArray` (new) are the only negation. `impl Neg for &ExactNum` is removed. Clippy's `should_implement_trait` forbids an inherent `neg(self)`, so there is none. Migration: `x.neg()` needs `use std::ops::Neg;` (or write `-x`); `x.neg()` / `-&x` on a value you still need afterwards becomes `ExactNum::from_i8(-1, p).mul(&x, p, rm)` (a fresh `-1` times the borrowed operand: exact when `p` is at least the precision of `x`, keeps the sign of zero, Inf, NaN and the accumulated inexact flag; the per-operation flag is false because nothing is rounded), or reorder so the negation is the last use of `x`. For a constant, build it directly (`ExactNum::from_i8(-1, p)`). `expr!`/`cexpr!` unary minus is unchanged. Accumath: `crates/engine_core/src/soft.rs:1532` and `:1539` (`self.inner.as_ref().neg()` on an `Arc<ExactNum>`) need the `from_i8(-1, p).mul(..)` form at the operand's precision; no other `ExactNum::neg` use outside the bundled zenith-float copy.
- `ExactNum::min`, `max`, and `clamp` return a reference to one of their operands (`&'a Self`) instead of a copy. A returned NaN is a shared static NaN. Migration: call methods on the result as before (auto-deref); where an owned value is needed, build it explicitly from the borrowed one (for example `ExactNum::round(..)` to the precision you want) instead of expecting a clone. Accumath: none of its calls use these on `ExactNum` (checked: no `ExactNum::{min,max,clamp}` call in crates/engine_core).
- Added `ExactComplex::into_parts(self) -> (ExactNum, ExactNum)`: moves both parts out without copying. Use it instead of `z.re().clone()` / `z.im().clone()` on a value you own. Accumath: hypergeom.rs 73-74 and elliptic.rs 105-106 can use it on upgrade.
- `binfmt`/`csvfmt` decoders return the `NAN`, `INF_POS` and `INF_NEG` constants directly (no change in behaviour).
- `Clone` is no longer derived for nine types that nothing in the library, tests, benches, examples, macro crate, compare crate, `nostd-tests`, doc examples or Accumath's pinned API surface copied: `Ball`, `ExactRational`, `GkResult`, `RadixFloat`, `RandomDist` (feature `random`), `CachedFBig`, `Ieee32Array`, `Ieee64Array` and `ExactNumArray` (breaking: `.clone()` on any of them no longer compiles). Only the `Clone` derive is removed from each; `Debug` and the other derives are unchanged. Migration: rebuild the value from its inputs, or keep it by reference. Accumath: none of these types is used. Still `Clone`, each pending its own step: `ComplexBall` (needs `ExactComplex: Clone`), `ExactInt`, `ExactNumPoly`, `ExactComplex`, `ExactNum`, and the compare-crate types (`DashuFloat`, `AstroFloat`, `ZenithFloat` and the `*State` types), which the `T: Clone` bounds of the compare crate's `task_for_*` functions and `BenchState`/`BenchFloat` still require. The `derive(Clone)` count in `tools/no_clone.sh` drops by nine; the `.clone(` baseline in `tools/clone_count_baseline.txt` (513) counts `.clone(` calls only, not derives, and is unchanged.
- (Planned, not yet done) `impl Clone for ExactNum` and every `Clone` derive on the big-number types are removed in the last step; migration notes per site follow then. Accumath on upgrade from `=1.0.12`: 4 `.cmp(` to `.compare(` at crates/engine_core/src/soft.rs:174, 179, 1476, 1482, and 23 `.clone()` lines on `ExactNum`/`ExactComplex`/`INF_*` (soft.rs 56, 61, 1208, 1213, 1409; hypergeom.rs 61-64, 73-74; elliptic.rs 105-106, 190-191, 222-223, 292-294; lfunction.rs 86, 94, 179); 9 of them go away with `ExactComplex::into_parts`, the `INF_*` consts used directly, borrowed slice inputs, and an out-of-place precision change; 14 need ownership moved into `ExactNum`/`ExactComplex` constructors.

### No lint exemptions

The library no longer contains any `#[allow]` or `#[expect]`. The underlying code was fixed instead:

- Eight crate-level clippy exemptions are gone from `zenith-float-num` (`manual_div_ceil`, `manual_is_multiple_of`, `neg_cmp_op_on_partial_ord`, `len_zero`, `slow_vector_initialization`, `comparison_chain`, `collapsible_else_if`, `collapsible_if`), and the first two are gone from `zenith-float`. The code uses `div_ceil` and `is_multiple_of`, and nested `if`s are collapsed. The `PartialOrd` tests write `(a < b).not()` instead of `!(a < b)`.
- `Radix` was a struct with associated constants `Bin`, `Oct`, `Dec`, and `Hex`, which needed `non_upper_case_globals`. It is now an enum with those four variants and `Other(OtherRadix)`. `Radix::try_new` returns the named variant for 2, 8, 10, and 16, so equal bases compare equal. `OtherRadix` has a private field, so only `try_new` can build one. Code that writes `Radix::Dec` or calls `try_new`, `value`, or `u8::from` is unchanged.
- `ExactNum::cmp(&self, &ExactNum) -> Option<SignedWord>` had the name of `Ord::cmp` with a different signature, which needed `should_implement_trait`. It is now `ExactNum::compare`, with the same behaviour. Migration: replace `.cmp(` with `.compare(` on `ExactNum`. A deprecated `cmp` alias was not kept, because it would still need the exemption. The `<`, `>`, and `partial_cmp` operators are unchanged.
- `Context::clone(&self) -> Result<Context, Error>` had the name of `Clone::clone` with a different signature, which needed `should_implement_trait`. It is removed, not renamed: nothing in zenith-float, Accumath or any other checkout on the development machine called it, and a copy of a context (its constants cache included) is better made explicitly. Migration: build a new one with `Context::new(p, rm, Consts::new()?, emin, emax)`.
- The internal `Mantissa::abs_add` and `abs_sub` take an `AddSubSpec`, and the internal FFT takes an `FftStep`, instead of eight arguments (`too_many_arguments`).
- The internal `mantissa::mantissa` module is now `mantissa::base` (`module_inception`).
- `zenith-float-macro` documents `expr!` and `cexpr!` and denies `missing_docs` instead of allowing it.

### Invariant violations are errors, not panics

Internal size and shape checks in the arithmetic kernels used to `unwrap()` or `panic!()`. They now return `Err(Error::InvalidArgument)` and the callers propagate it: the FFT multiplier (`fft_forward`, `fft_reverse`, `fft_decompose`, `fft_mul_mod`, `fft_div_mod`), Toom-2 and Toom-3, `div_basic` and its correction step, `SliceWithSign::{add, sub, copy_from, decrement_abs}`, `root_estimate`, `sqrt_rem`, `div_rem_abs` on empty operands, `round`, and the binary/decimal converters. `ExactNum::clamp`, `min`, `max` and the parser no longer unwrap. The rounding `match`es list every `RoundingMode` variant instead of ending in `unreachable!()`. `expr!` and `cexpr!` report a failure to build the constant cache as a `compile_error!` instead of panicking inside the compiler. Results for valid inputs are unchanged. The three library crates now deny `clippy::unwrap_used`, `clippy::expect_used`, and `clippy::unreachable` outside `cfg(test)`, so they cannot come back. Still open: `clippy::panic` (one site, `SliceWithSign`'s `DerefMut` on an immutable view) and unchecked indexing (see PROGRESS.md).

### Constants without start-up panics

The internal constants (`ONE`, `TWO`, `THREE`, `FOUR`, `TEN`) were `lazy_static` values built with `expect`. They are now `static` values built by `const fn`, so there is nothing to allocate and nothing to fail, and the `lazy_static` dependency is gone. That needs 128 bits to fit inline on every target, so `INLINE_WORDS` is `DEFAULT_P / WORD_BIT_SIZE`: 2 with a 64-bit `Word` (unchanged) and 4 with a 32-bit `Word` (was 2). On 32-bit targets a 128-bit `ExactNum` is now stored inline, and `is_inline()` is true for it. The decimal conversion cache behind `FromStr` and `Display` is created on first use; if it cannot be created, `from_str` returns `Err(MemoryAllocation)` and `fmt` returns `fmt::Error`. Eight unused constants were deleted, and the public-module `ONE` and `TWO` in `ext` became test-only.

### Build without default features

`cargo build --no-default-features` failed because the `hex_golds` binary needs `std`. The binary now declares `required-features = ["std"]`, so a `no_std` build skips it (021fa87).

### Special-function oracle under debug builds

The per-row wall-clock budgets in the special-function oracle are for optimised builds. A debug build under load timed out `struve_h` at 53 bits in the prepublish run, though the whole file takes 1.6 s alone. Debug builds now get ten times the budget (6bcf901).

### Exact-`p`-bit rounding beyond the basic operations

Only `add`, `sub`, `mul`, and `div` rounded at exactly `p` bits; every other function rounded at `p` rounded up to a word and then again at `p`, so its last bit could be wrong in any directed mode and at ties. All rounded functions now run toward zero with a guard word, round to odd, and then round once at `p` bits in the requested mode. An exact result (`sqrt(4)`, `exp(0)`, `log(1)`, `log10(1000)`, `hypot(3, 4)`, `9^0.5`, `nth_root(16, 4)`) is exact and reports `inexact() == false`. Fixes found on the way:

- `hypot` never terminated on exact results (`PrecisionRetryExhausted` for `(3, 4)`); it now forms the sum of squares exactly and takes one rounded square root.
- `fma` dropped a far-smaller addend, which gave wrong directed results and flags; an exact cancellation is `-0` under `Down`.
- `nth_root` chained square and cube roots and rounded twice.
- Parsing rounds the decimal value once at `p` bits.

Results can differ from 1.x in the last bit. Code that compared bit patterns with 1.x output has to regenerate them.

### Overflow follows the rounding mode in every function

`gamma` and `Ei` returned `+Inf` on overflow in every mode. They now return `±max_value(p)` under `ToZero`, under `Down` for a positive result, and under `Up` for a negative result, like the other functions. Poles (`gamma(0)`, `Ei(0)`) stay infinite. `ToOdd` rounds to nearest, so it overflows to infinity like `ToEven`.

### `ExactComplex::inexact` is per operation

`ExactComplex::inexact()` was true when either part, or anything it was computed from, was rounded. It now reports whether this operation rounded either part, like `ExactNum::inexact()`; `inexact_accumulated()` keeps the old meaning. `sqrt` of a complex number with an exact root (`sqrt(3+4i) = 2+i`) is exact.

### `ldexp` underflow

`ldexp(x, n)` and `scalb` returned zero whenever the shifted exponent fell below `EXPONENT_MIN`, even when the result was a representable subnormal. They now round onto the subnormal grid in the requested mode.

### 32-bit clippy

`scripts/ci_hex_32bit.sh` runs clippy on `i686-unknown-linux-musl`. Two lints that only fire with 32-bit words were fixed (`modulo_one` in `binfmt.rs`, `needless_return` in `num.rs`).

### `sin_cos` and `sinh_cosh`

The pair functions still rounded at a word multiple, set both flags from the shared computation, and returned `±Inf` on `sinh_cosh` overflow in every mode. Each result now goes through the same exact-`p`-bit rounding as `sin`, `cos`, `sinh`, and `cosh`, with its own per-operation flag and mode-aware overflow. The two values are computed separately, so a call costs about as much as the two single functions.

### `rem` is exact

`rem` takes no precision or rounding mode because the remainder of two floating-point numbers always fits: it is smaller than both operands and its last bit is no finer than the finer operand's. A sweep against an exact big-integer remainder (all test precisions, operands far apart, close, and in the subnormal band) confirms it. The sweep found one bug: a remainder below `EXPONENT_MIN` was always flushed to zero, because the range check compared with `usize::MAX as i64`, which is `-1`. Subnormal remainders are now returned exactly. `inexact()` on a remainder is false; `inexact_accumulated()` keeps the operands' history.

### Full-precision sum with a zero operand (bug fix)

`ExactNum::add_full_prec` and `sub_full_prec` returned a wrong value when exactly one operand was zero: `0.add_full_prec(1)` returned a zero with an empty mantissa (exponent 1) flagged inexact, instead of 1. Both now return the other operand exactly (negated for `sub` when the zero is the first operand), with that operand's own `inexact()`. The fixed-precision `add` and `sub` were not affected and are bit-identical to 1.x. Found by the add/sub grid digest test while removing the operand copies: it differs from the pre-change code in exactly 40 results, all of this kind. Test: `add_sub_zero_tests::full_prec_with_zero_operand_is_the_other_operand`; the grid digest is pinned in `grid_digest_tests::add_sub_grid_digest`.

### Rust 1.99 patch releases in CI

`scripts/use_rust_1_99.sh` installs channel `1.99` and accepts any `1.99.x` compiler (it used to demand exactly 1.99.0, which broke the `hex-32bit` and `prepublish` steps once 1.99.1 was released). Other minor releases are still refused.

### Deterministic `test_compute_added_err`

Under `--features random` the `random_dist` tests reseed the shared test generator from their own threads, which restarted the stream of `macro_util::test_compute_added_err` part-way through and made it fail at random. The test now pins its own stream (`test_rng::pin_thread_seed`, test builds only). This is a test-only change.

### Subnormal underflow

`mul`, `div`, and `add` flushed a result below the minimum positive subnormal to zero whenever the rounding mode was not `Up`, `Down`, or `FromZero`. The discarded significand was ignored, so a value above half an ulp of the minimum positive, and a tie under `ToOdd`, became zero. `ToOdd` of half the minimum positive must be the minimum positive, because zero is even. Directed modes that round away from zero already kept that minimum; `ToEven` of an exact tie correctly stays zero.

The flush now compares the residual with half an ulp and applies the same rounding rule as every other inexact result. Code that depended on flushing a tie, or a value above half an ulp, to zero under `ToEven`, `ToOdd`, or `None` now sees the rounded minimum positive or a kept zero, the same rule the directed modes already used.

### Signed zero

`(+0) + (-0)` returned the sign of the second operand, so the sum was `-0` in every rounding mode. IEEE 754 addition of zeros with opposite signs is `+0`, except for `roundTowardNegative` (`Down`), which is `-0`. The same rule applies to exact cancellation: `1 + (-1)` and `5 - 5` are `+0` except under `Down`, where they are `-0`. `(-0) + (-0)` stays `-0` in every mode. Code that required those sums to be `-0` outside `Down` has to round with `Down` or set the sign itself.

### Invalid operations carry `InvalidArgument`

`0 * Inf`, `Inf * 0`, `Inf - Inf`, `Inf + (-Inf)`, `Inf / Inf`, and remainder of an infinity were `NaN` with an empty error tag. Callers could not tell them from an uninitialized `NaN`. They are now `NaN(InvalidArgument)`. `0 / 0` was already tagged and is unchanged. Code that treated `err() == None` as the result of these operations now sees `Some(InvalidArgument)`. Use `is_nan()` when the tag does not matter.

### Parsing consumes the whole string

`parse` used to accept a valid prefix and ignore the rest. `1.5foo` was `1.5`, `0x10` was `0`, `1.2.3` was `1.2`, and `infinity` was `Inf`. A bare exponent was the same class of bug: `1e`, `1e+`, `1e-`, and `1E` were read as `1` with an exponent of zero. Strings that were not numbers at all, including `""`, `"."`, `"nanx"`, and `++1`, came back as `NaN` with no error tag.

`parse` now requires the entire string to be a finite number, `Inf`, or `NaN`. Anything else is `NaN(InvalidArgument)`. Well-formed `NaN` and `Inf` are unchanged, including their signs for infinities. This is a breaking change. Code that passed a longer buffer and relied on the prefix must slice the number out before calling `parse`. Code that used a bare `e` to mean exponent zero must write `e0` or omit the exponent. A digit outside the radix is the same break: `2.5` in base 3 used to be read as `2` and is now `NaN(InvalidArgument)`. `format` of that `NaN` prints `NaN`, and parsing that text back is still a `NaN` (the error tag is not part of the spelling). A `NaN` from untrusted text now has `err() == Some(InvalidArgument)` instead of `None`; check `is_nan()` if the tag does not matter.

### Decimal round-trip

`format` rounded the decimal digits with the same directed mode that `parse` would apply again. Under `Up` and `ToZero`, `parse(format(x))` moved one ulp. `0.1` at 53 bits printed `1.0000000000000000001e-1` with `Up` and parsed back to the next number up. `ToEven` was unaffected.

`format` in base 10 now retries with extra digits. Those digits are rounded onto the side of the value that the same mode accepts, and the retry stops when `parse` at the same precision returns the same sign, exponent, and limbs. The spelling of a directed value can be longer than in 1.0.12. The value is the one `parse` produced. Code that matched the old spelling has to accept the longer text.

### Integer powers

`exp2` and `exp10` of an exact integer, and `pow` with an exact integer exponent, went through `exp(n · ln x)`. When the power is an exact dyadic, a directed mode sees every discarded bit equal and the Ziv loop gives up. `exp2(2)`, `exp2(20)`, `exp10(2)`, and `exp10(20)` at 24, 53, and 113 bits returned `NaN(PrecisionRetryExhausted)` for `Up`, `Down`, `ToZero`, and `FromZero`. `None` returned the next number down. `ToEven` and `ToOdd` were already the exact power. `pow(0.5, -2)` is the same bug at every width, including 256 and 1024 bits: the value is the integer 4.

An exact integer exponent whose magnitude fits in one word, and whose base is exact, is now exponentiated by squaring at full precision and rounded once. A negative exponent divides that product once. Inexact inputs still use the logarithm. `NaN(PrecisionRetryExhausted)` on these exact powers becomes the integer (or the correctly rounded integer, when `10^20` does not fit). Code that treated that `NaN` as "not an integer power" now receives the integer.

### Constant cache

`sqrt2` and `phi` at 113 bits under `Up` and `FromZero` were one ulp low when a coarser request had been cached first. A 24-, 53-, or 64-bit request is computed at 128 bits. The 113-bit request is also stored in 128 bits, so the cache was reused and `set_precision` did nothing. A fresh cache, or a cache filled at 256 bits, was already the higher neighbor.

The cache is reused only when it still has a guard word below the new precision. Otherwise it is recomputed. `sqrt2(113, Up)` and `phi(113, Up)` after a narrower request now match a fresh cache, and the same holds for `FromZero`. Euler's constant is filled again when the saved value has no guard word left, so a later, wider request does not treat that cache as final. Code that reused a coarse cache and expected the low neighbor at 113 bits under `Up` or `FromZero` now gets the higher neighbor.

### Complex inverse sine on the cut

`asin` of a real number with absolute value greater than 1 ignored the sign of a zero imaginary part. Both sides of the cut came back with the principal imaginary part: negative for `x > 1`, positive for `x < -1`. The value above the cut (`+0`) is unchanged and matches the principal value. The value below the cut (`-0`) is now the conjugate. The parts are rounded to the word-sized precision, same as the real elementary functions. A request for 24, 53, or 113 bits no longer stops at that many bits inside a wider word. Code that read the sign of a zero imaginary part on the cut and expected the principal value on both sides now gets the conjugate below the cut.

### Complex logarithm near the cut

`ln(z)` took the modulus at the destination precision and then the logarithm. For `z = -1 + 10^{-12} i`, `|z| - 1` is about `5·10^{-25}`. At 64 stored bits that difference is below one ulp of 1, so the real part was `+0`. The imaginary part, near `±π`, was already right. The same cancellation hides `2^{-p}` at 113 bits.

The real part is now `ln(max(|re|, |im|)) + ½ ln(1 + (min/max)²)`, computed with extra precision and rounded once. A logarithm just off the cut keeps that small positive real part. The argument is unchanged. Code that treated a zero real part just off the cut as exact now has to accept that small positive part.

### Accuracy labels

Public numeric functions now carry an `Accuracy` label in the rustdoc: correctly rounded, working precision (with a lost-bit bound), an enclosure, or unspecified. `Accuracy`, `accuracy_of`, and `ACCURACY` are new. A function that was already documented as a Ziv leaf stays correctly rounded. Lambert \(W\), Kelvin, Scorer, Clausen, Barnes \(G\), Weber, Anger, Struve, zeta, Hurwitz, \({}_1F_1\), \({}_0F_1\), \(pFq\), the polylogarithm, polygamma, Meijer \(G\), Fox \(H\), and Heun stay unspecified. Orthogonal polynomials and Legendre stay working precision, 32 bits. This does not change any value.

### Riemann zeta at 0 and −1

`riemann_zeta(0)` is \(-1/2\) and `riemann_zeta(-1)` is \(-1/12\). Both used to be `NaN(InvalidArgument)`. This is a breaking change. Callers that treated 0 and −1 as out of domain should test the argument before calling. Other non-integers, the pole \(s=1\), and integers \(s\le -2\) stay `NaN(InvalidArgument)`.

### Clausen

\(\mathrm{Cl}_3\) of a nonzero angle integrated \(\mathrm{Cl}_2\), and \(\mathrm{Cl}_2\) was itself a quadrature, so the pair did not return. Both are now the Bernoulli series on \((-\pi,\pi]\). \(\mathrm{Cl}_2(0)=\mathrm{Cl}_2(\pm\pi)=0\) and \(\mathrm{Cl}_3(0)=\zeta(3)\) are unchanged. There is no timeout argument. The values move by more than an ulp; the class stays unspecified. Callers who pinned the old digits of \(\mathrm{Cl}_2\) or \(\mathrm{Cl}_3\) have to accept the Bernoulli values.

### Weber \(E_0\)

\(E_0(z)\) for \(|z|<8\) is \(-(4/\pi)\sum J_{2k+1}(z)/(2k+1)\). The quadrature at 256 bits had the wrong leading digits (`E_0(0.3)` was about \(-0.043\) instead of \(-0.189\)). Other orders still use the integral. The class stays unspecified. Callers who pinned the old leading digits of \(E_0\) on \(|z|<8\) have to accept the series.

### Barnes \(G\)

Non-integer \(G\) at 256 bits was about 60 bits off. \(\ln A\) uses enough Euler–Maclaurin terms in \(\zeta'(2)\) for the whole word. Integer arguments stay the factorial product. The class stays unspecified. Callers who pinned a non-integer \(G\) at 256 bits have to accept the corrected value.

### Word rounding of catalog and orthogonal polynomials

Catalog and orthogonal-polynomial results were rounded to exactly `p` bits. A request for 53 or 113 bits is stored in a wider word, and the low bits of that word were zero. They are now rounded to the full word (`p` rounded up to a multiple of 64), the same width the elementary functions already use. Word-aligned requests (64, 128, 256) are unchanged when the series was already accurate to that word. Code that compared a 53- or 113-bit catalog value for limb identity with a `p`-bit rounding has to compare at the full word.

### Elliptic argument

`elliptic_f`, `elliptic_e`, and `elliptic_pi` already took the sine of the amplitude and the parameter \(m=k^2\). The rustdoc and the capabilities table now say that in those words. \(|x|>1\), including an angle such as 2 passed in place of \(\sin 2\), stays `NaN(InvalidArgument)`. The argument meaning did not change, and call sites outside this crate were not edited. Callers who passed the angle into that slot should call `elliptic_f_phi`, `elliptic_e_phi`, or `elliptic_pi_phi`.

### Gauss–Kronrod error estimate

`GkResult::abs_err_est` was documented as `|K − G|`, and as the difference of the last two even-order ε entries when that difference stopped the iteration. On a smooth integrand the two rules share the same rounded nodes, so they agree far more closely than either agrees with the integral. At 53 bits, \(\int_0^1 e^x\,dx\) was high by about \(8\cdot10^{-17}\) while the reported error was about \(10^{-19}\), roughly 700 times too small. The same bar missed \(\int_0^1 \sin x\,dx\). On \(\int_0^1 x^{-1/2}\,dx\) the ε difference was about \(10^{-56}\) while the integral was high by about \(5\cdot10^{-21}\). Gauss–Legendre, Gauss–Laguerre, and Gauss–Hermite used the same short working precision, so a low-degree polynomial could miss by many units of the stored word.

The (7, 15) rule and those Gauss rules now run at two guard words above the destination width (that width rounded up to a multiple of 64, plus 128 bits). Tanh–sinh keeps its wider double-exponential guard when that is already more. `abs_err_est` is a conservative error estimate, validated against oracles: at least `|K − G|` on one panel, or the sum of the local `|K − G|` when an ε value is accepted, plus one unit in the last place of the rounded estimate. The estimate is rounded away from zero. It can be larger than `|K − G|`, and it can be larger than the tolerance that stopped the iteration. `estimate_kronrod` and `estimate_gauss` keep their meaning and still use the caller's rounding mode. The ε difference remains the stopping test for an extrapolated value; it is not the published estimate.

This is a breaking change. Code that read `abs_err_est` as `|K − G|` should subtract `estimate_kronrod` and `estimate_gauss`. Build the struct with `GkResult::new`. Do not rely on field order.

### Elliptic shortcuts round once

\(F(\varphi\mid 0)\) and \(E(\varphi\mid 0)\) are \(\varphi\). The angle-form path reduced \(\varphi\) through sine and arcsine at a working precision, then rounded again. `Up`, `Down`, `ToZero`, `FromZero`, and `None` failed the roundability test and returned `NaN` for \(\varphi=0.8\) at 256 and 512 bits. `ToEven` and `ToOdd` happened to land on the input. Those shortcuts now round the amplitude once, in the caller's mode. \(\varphi=0\) stays exact zero. \(\varphi=k\pi/2\) with \(m=0\), including \(\pm\pi/2\), \(\pi\), and \(3\pi/2\), is that same rounding. The sine-form \(m=0\) shortcut is one arcsine. \(\Pi(0;\varphi\mid m)\) calls \(F\) once. \(n=0\) on the sine form already called \(F\).

### Infinite-range Gauss–Kronrod

`integrate_adaptive_gk_infinite` maps a half-line or the whole line onto a finite interval and uses the same (7, 15) rule. `abs_err_est` on those rows is the same conservative error estimate, validated against oracles, and the estimate itself is checked against `est_tol`.

### Platform intrinsics

`add_carry` and `sub_borrow` in the limb kernel called `_addcarry_u64` and `_subborrow_u64` (and the 32-bit pair) through `core::arch`. That code came with `common/util.rs` from astro-float and still carried the MIT notice. On x86-64 the intrinsic was the live path for every addition, subtraction, and product of `ExactNum`. `shift_slice_left` and `shift_slice_right` walked those limbs through raw pointers. Decimal conversion used `Vec::set_len`. `ieee_soft` added and multiplied software binary32 and binary64 lanes with x86-64 and AArch64 integer SIMD (`_mm_add_epi32`, `_mm_add_epi64`, `_mm_mul_epu32`, `vaddq_u32`). None of those paths used `f32` or `f64`. There is no `asm!`, `global_asm!`, or `target_feature`.

The limb operations are now portable integer arithmetic: a widening sum for the carry, a widening difference for the borrow, and index walks for the shifts. Decimal conversion fills the buffer with `resize`. The software IEEE lanes use `wrapping_add` and a widening multiply. For a carry or borrow of 0 or 1 the new limbs match the old intrinsics, including a chain of `WORD_MAX`. `ieee_soft` is not from astro-float; its public bit patterns are unchanged. The astro-float MIT notice on `common/util.rs` and `mantissa/conv.rs` stays.

`scripts/ci_no_intrinsics.sh` rejects `std::arch`, `core::arch`, `asm!`, `global_asm!`, `target_feature`, `f32`, and `f64` in every Rust source, and rejects `unsafe` outside `scripts/no_intrinsics_allowlist.txt`. The allowlist is the `no_std` host and qemu harnesses, plus the optional MPFR test oracle. Those are not linked into the library. The library crates forbid `unsafe_code`. Callers do not see a signature change. A build that depended on the x86 intrinsic being the only implementation of a limb add now gets the same bits from the portable sum.

### Directed rounding modes were documented as ties

The rustdoc on `RoundingMode` called `Up` "round half toward positive infinity", `Down` "round half toward negative infinity", `ToZero` "round half toward zero", and `FromZero` "round half away from zero". Those words describe round-to-nearest with a tie rule. The kernel (`round_mantissa`) has always rounded these four modes in one direction: any discarded bit moves `Up` to the next value above, `Down` to the next value below, `FromZero` away from zero, and `ToZero` truncates. A caller who trusted the rustdoc expected a nearest result and got a directed one.

The rustdoc now says directed rounding toward +∞, −∞, zero, or away from zero. `ToEven` and `ToOdd` are still round-half-to-even and round-half-to-odd. No value changes. `tests/rounding_mode_docs.rs` checks the wording and sweeps non-representable quotients over four precisions and both signs to show each mode lands on its side of the exact value.

### `from_raw_parts` skipped its checks

`ExactNum::from_raw_parts` documents that it returns `NaN` when `n` is larger than the mantissa, when `n` does not match the mantissa's significant bits, when a subnormal mantissa does not carry `EXPONENT_MIN`, when the size is too large, and when `e` is outside `EXPONENT_MIN..=EXPONENT_MAX`. It copied the words into a mantissa and called the unchecked constructor, so none of that happened. `from_raw_parts(&[0xC000…], 65, Pos, 2, false)` on one 64-bit word came back as a finite value with a 65-bit length on a 64-bit mantissa, and a one-bit subnormal mantissa was accepted with any exponent. Later operations on those values read past the significand or misplaced the binary point.

The public constructor now goes through the checked `ExactNumNumber::from_raw_parts`, which was previously compiled only for unit tests. Every documented invalid case is `NaN(InvalidArgument)`; nothing panics. Parts that came from `as_raw_parts` still round-trip unchanged. Code that built numbers from inconsistent parts now gets `NaN(InvalidArgument)` and has to pass the bit count and exponent `as_raw_parts` would report. `tests/from_raw_parts_checks.rs` round-trips 4 precisions × 4 mantissas × 7 exponents including both limits × both signs, and rejects each inconsistent shape at every precision and sign.

### Exponent overflow follows the rounding mode

Every exponent overflow came back as `±Inf`, whatever the rounding mode. `num.rs` reported `ExponentOverflow(sign)` from `mul`, `div`, `rem`, and `add`, and `result_to_ext` turned that into an infinity without looking at `rm`. IEEE 754 and MPFR return the largest finite value when the mode rounds toward zero for the result's sign: `ToZero` always, `Down` for a positive result, `Up` for a negative one. `max_value(p) * 2` under `ToZero` was `+Inf`; it is `max_value(p)`. The oracle hid this twice. `tests/mpfr/common.rs` accepted any `Inf` once MPFR's sticky overflow flag had been raised by any earlier operation, and `scripts/gen_oracle_goldens.py` and `scripts/mpfr_oracle.c` both wrote `Inf` for every mode, so the goldens agreed with the bug.

`add`, `sub`, `mul`, `div`, `fma`, `mul_add`, `pow`, `powi`, `powsi`, `hypot`, `reciprocal`, `ldexp`, `exp`, `exp2`, `exp10`, `sinh`, and `cosh` now return `±max_value(p)`, marked inexact, in those three cases and `±Inf` in every other mode (`ToEven`, `ToOdd`, `FromZero`, `None`, `Up` for a positive result, `Down` for a negative one). Poles are not overflows: `0^(-n)`, `1/0`, and functions that signal a pole with `ExponentOverflow` (`atanh(±1)` and the like) stay `±Inf`. `reciprocal` above the Newton threshold also reported a positive overflow for a negative argument; the sign is now the argument's. Code that tested `is_inf()` to detect overflow under `ToZero`, `Down`, or `Up` has to compare with `max_value(p)` as well.

The MPFR check now requires `Inf` exactly where MPFR is `Inf` and fails when MPFR is `Inf` and zenith is not. Both generators emit the mode-correct value, the goldens were regenerated (`edges.jsonl` gains nine overflow rows per precision and mode, both signs, for `add`, `sub`, `mul`, and `div`; `parse.jsonl` ids shift), and `tests/overflow_modes.rs` sweeps 14 overflowing operations over all seven modes and four precisions. `fuzz_mpfr_overflow_all_round_modes` does the same against MPFR at five precisions.

### `inexact` is per operation

`add`, `sub`, `mul`, and `div` started from `self.inexact || d2.inexact` (`num.rs` 337, 409, 633) and then added whatever this operation rounded. The flag on the result answered "was anything upstream rounded", not "did this operation round". `(1/3) * 2` is exact, and its `inexact()` was `true` because `1/3` had been rounded. A caller could not tell from the flag whether the operation it just performed was exact.

`ExactNum::inexact()` on the result of `add`, `sub`, `mul`, or `div` is now true only when that operation rounded, including rounding onto the subnormal grid. The operands' flags are not inherited. This is a breaking change: code that read `inexact()` to learn whether a chain of operations was exact has to call the new `ExactNum::inexact_accumulated()`.

The accumulated flag is still needed, so it is kept as that separate API. The correct-rounding loops decide whether a value needs a roundability check from it (`try_set_precision`), negative `powi` uses it to know whether repeated squaring was exact, and the `expr!`/`cexpr!` macros use it to know whether a cancellation can carry error. Each number now carries both flags; the kernel's own `add`/`sub`/`mul`/`div` still accumulate, the public wrappers also record the per-operation flag, and any in-place change (`set_precision`, `set_inexact`) drops back to the accumulated flag; subnormalization keeps the per-operation flag and adds its own rounding to it. Other functions (`sqrt`, the transcendentals, `ExactComplex`) report the accumulated flag through `inexact()`, as before. `set_inexact` sets both. The macros, `macro_util`, `ExactComplex::inexact`, and the Heun exactness test read `inexact_accumulated()`, so their behaviour is unchanged. `tests/inexact_per_op.rs` sweeps six modes, four precisions, and both signs over ten exact operations on an inexact operand, two rounding operations, and the accumulated chain.

### `add`, `sub`, `mul`, and `div` round at exactly `p` bits

These four operations rounded at `round_p(p)`, the next multiple of 64, whatever `p` was asked for. `a.add(&b, 24, ToEven)` was a 64-bit result; a caller asking for binary32 or binary64 behaviour got neither, and the requested precision only decided the storage. `set_precision` was already bit-granular, so `chop` was the one place that honoured `p`. A second fault sat under every precision: a result that fell below `2^EXPONENT_MIN` was rounded twice, once at `p` bits and again onto the subnormal grid, so it could be one unit off in the last place. The MPFR sweep found `p = 64` products with `ToEven` landing on the wrong neighbour.

`ExactNumNumber::exact_bits` now evaluates the four operations. With `p` a multiple of 64 and a normal result, it calls the operation as before. Otherwise it computes toward zero at `round_p(p) + 64` bits, folds the sticky bit into the last bit (round to odd), and rounds once more to `p` bits with `rm`. Two bits of headroom make that second rounding exact, so the result is correctly rounded at `p` bits in every mode, including `ToOdd`, and on the subnormal grid `2^(EXPONENT_MIN - p)`. An exact zero is recomputed directly so that its sign follows the mode. An overflow toward zero returns `max_value(p)` truncated to `p` bits. `subnormalize` no longer drops the per-operation `inexact` flag, which the round-to-odd step depends on. Storage is still a whole number of words; the bits below `p` are zero. `RoundingMode::None` skips rounding and keeps the storage precision, as before. Unaligned precisions cost one extra word and a second rounding.

The scope is limited to these four. `sqrt`, `cbrt`, `hypot`, `fma`, `rem`, parsing, and every function built on the kernels still round at `round_p(p)`, and their documentation still says so. Their golden rows are labelled with that effective precision, not with the request. The constant-cache rows keep `p = 113` because they test a 113-bit request after a narrower one. The goldens use `p` bits for the four exact ops (the operand tokens remain the zenith inputs at storage precision) and were regenerated; `fuzz_mpfr_exact_bit_precision` checks the four operations against MPFR. It runs eleven precisions (24, 53, 64, 65, 100, 113, 127, 128, 130, 200, 255) over three exponent bands, including the subnormal range, in six modes plus `ToOdd` in the normal range, with random signs, and compares bit for bit, subnormals included.

### The rounding-mode fuzz covered three operations and one precision

`fuzz_mpfr_bit_exact_all_round_modes` ran `add`, `mul`, and `sqrt` at 128 bits only, in five modes, on uniformly random operands. It never ran `sub` or `div` and never ran `ToOdd`, and random operands almost never hit an exact tie, so most of the tie-breaking code had no MPFR check. With `ToOdd` ties deliberately broken toward even, the old test still passed.

The test now runs `add`, `sub`, `mul`, and `div` bit for bit against `mpfr_reference` in all six rounding modes, `ToOdd` included (MPFR has no ties-to-odd mode; the reference emulates it with a toward-zero result at `p + 1` bits and the ternary value). It uses precisions of 64, 128, 192, 320, and 1024 bits. Half of the operands are runs of ones followed by zeros with random signs, which produces exact ties and carries out of the top word. `sqrt` and `nth_root` keep the five MPFR modes at every precision. Under the same broken `ToOdd`, the new test fails at the 13th case: a 128-bit `Mul`. Unaligned precisions and the subnormal range are covered by `fuzz_mpfr_exact_bit_precision`. `get_random_rnd_pair` still returns only the five MPFR modes, because its 44 callers pass the MPFR mode straight to MPFR; `ToOdd` is covered through `mpfr_reference` instead. No library change.

### The exponent range is the same on 32-bit targets

`EXPONENT_MAX` and `EXPONENT_MIN` were `i32::MAX / 4` and `i32::MIN / 4` on 32-bit pointers (`defs.rs` 49–58). The kernel computed exponents in `isize`, which is 32 bits there, and the quarter range was headroom against `e1 + e2 - shift` overflowing. The same value could therefore be finite on x86_64 and `Inf` on i686 or Cortex-M. That broke the claim that results are bit-identical on every platform, and the hex golds never reached the limits.

Exponent arithmetic is now done in `i64` on every target, and both constants are `i32::MAX` and `i32::MIN` everywhere. Every `isize` in `src` became `i64`. Conversions from an exponent difference back to a count go through `common::util::abs_usize`, which saturates where `usize` is 32 bits. `p_assertion` also rejects `p >= usize::MAX / 4`, so precision arithmetic cannot wrap a 32-bit `usize`; the 64-bit limit is unchanged. The checks and clamps that existed only for 32-bit targets became vacuous and were removed. They were in `from_words`, `from_raw_parts`, `set_exponent`, `random_normal`, the radix conversion, and `ldexp`'s `exponent_in_range`, which took a `#[allow]`. `parser` tests use the same `e+7FFFFFFF` overflow input on both widths.

This is a public change for 32-bit users. Values between `2^(i32::MAX/4)` and `2^i32::MAX` are now finite instead of overflowing, and they underflow at `2^i32::MIN` instead of `2^(i32::MIN/4)`. `doc/LIBRARY.md` and `doc/ZENITH_FLOAT_CAPABILITIES.md` are updated. `tests/exponent_range.rs` checks the constants, `max_value`, `min_positive`, and `min_positive_normal` at four precisions. It sweeps `ldexp` to three exponents at each limit and `mul`, `add`, `sub`, and `div` at the top and bottom of the range, including overflow and subnormal results, in six modes and both signs. It failed on `i686-unknown-linux-musl` before the change. With the change, it passes there in debug, with overflow checks on, together with the library unit tests and the property, edge, overflow, and `inexact` suites. `tests/oracle.rs` does not build for 32-bit words (it passes `u64` limbs), as before.

### Correct rounding within 2^-2p of a boundary at high precision (bug fix)

`cos` of 4π taken at 100032 bits, rounded to 100032 bits with `FromZero`, `Down`, or `Up`, returned `NaN(PrecisionRetryExhausted)`; MPFR returns 1.000…. The result is 1 − 2^-(2p−5), so the rounding has to see about 2p bits before it can tell it from 1. Every correctly rounded function runs a Ziv loop that raises its working precision until the rounding is certain, and every one of them stopped at `p + 64·256` working bits (`PrecisionRetryExhausted`): enough up to about p = 16000, not above. The fixed 16384 extra bits are replaced by `max(64·256, 2p)` extra bits, so the working precision may reach `max(p + 16384, 3p)` (`prec_retry_budget` in `common/util.rs`, shared by all 72 callers of `bump_prec_retry`). The loop still ends: the working precision grows by at least one word per round and the limit is a finite number, so there are at most `limit/64` rounds, and a result that cannot be certified (an exact dyadic that is not recognised) still ends in the error, a little later.

Functions affected, each failing at 100032 bits before the fix and correct against MPFR in all six rounding modes after it: `cos`, `sin`, `sin_cos`, `cosh`, `sinh`, `tanh`, `exp`, `expm1`, and `log1p` (the shapes: `cos(2π)`, `cos(4π)`, `sin(π/2)`, `exp(2^-p)`, `log1p(2^(1-p))`, `cosh` near 1). `tan` (100032 bits), `ln`, `log` with a base, and `nth_root` (40064 bits) go through the same loop and have the same-shape tests; `tan` passed before the fix at the size tried, and `ln`, `log`, and `nth_root` passed before the fix at 20032 bits (the new tests run at 40064 bits, which was not run against the old budget). Results that used to be `NaN(PrecisionRetryExhausted)` are now the correctly rounded number. No number that was returned before changes because of this: the pinned digests `add_sub_grid_digest`, `mul_div_grid_digest`, `sqr_grid_digest`, `jacobi_all_functions_digest`, and the `sqrt` digest are identical (the one pin that changed is explained in the next section).

### `log1p` and `expm1` of a tiny argument rounded the wrong way (bug fix)

For an argument so small that the result is the argument plus one correction bit, `log1p(x)` with `x > 0` and `expm1(x)` with `x < 0` put the correction on the wrong side. Under `Up` (and `Down` for the opposite sign) the result was one unit in the last place too large (too small) compared with MPFR: `log1p(sqrt(2)·2^-2147483648)` at `p = 20032` with `Up` returned the next number above the argument, where the correctly rounded result is the argument itself, and `expm1(-sqrt(2)·2^-20097)` likewise. The wide-range sweep of the new `retry_budget` tests (both signs, from the bottom of the exponent range to 2, all six modes, 128, 20032 and 100032 bits) found it. Fixed. The pinned `ziv_entry_points_digest` recorded the wrong answer for one operand (the positive subnormal operand 21 of its grid, `log1p` with `Down` at `p = 64` and `p = 200` returned the argument itself; the correct result is the argument minus one unit of the subnormal grid), so that pin is re-recorded: 2735357469467868262 becomes 6378838253543692452, and only those two entries differ (checked entry by entry). The other users of the shortcut (`sin`, `cos`, `tan`, `sinh`, `cosh`, `tanh`, `atan`, `exp`) were checked and their correction already goes the right way for both signs.

### The Ziv roundability check accepted numbers it could not decide, and refused numbers it could (bug fix)

`Mantissa::round_mantissa` answers "is this value roundable, given that only its leading `s` bits are correct?". For `ToEven` and `ToOdd` it was wrong in three ways: it looked at the *lowest* words (the unknown ones) instead of the known words next to the round bit, it counted the word that holds the round bit twice (`cc` was one too large when the round bit sat at bit 63 of its word), and its test of the bits just below the round bit was a tautology (`rb & mask == rb`). A value just below a tie whose known bits could not decide the rounding was accepted whenever the unknown words happened to look different from all ones. Example, checked against MPFR: `log1p(-sqrt(2)·2^-20032)` at `p = 20032` with `ToEven` returned the number below the correctly rounded one (the true value is 2^-20034 ulp above the tie and `p_wrk` was 20096); `tests/mpfr/retry_budget.rs::roundability_check_log1p_hard_case_at_20032_bits` is that case. The same check also refused values it could decide: to nearest only a round bit of one followed by known zeros, or a round bit of zero followed by known ones, is close to a tie, but a value just below a representable number (round bit one, then ones) was treated as undecidable too, so `ToEven` returned `NaN(PrecisionRetryExhausted)` for 23 entries (at two precisions) of the Jacobi grid at the subnormal argument `u = 3·2^(-2^31-19)` where the answer is simply `u` or 1.

The check now looks at the known words from the one next to the round bit downwards, counts the round-bit word once, compares the bits below the round bit with all zeros / all ones properly, and for the nearest modes calls a value undecidable only when the pattern is close to the half-way point (`round bit 1` + known zeros, `round bit 0` + known ones). The directed modes were already right apart from a harmless tautology (it can only cost a retry) and are unchanged.

Because the old check accepted exact results by the accident that made it wrong, the loops that can produce an exactly representable result needed their own exactness test. `pow` was the one that failed: a power of two to an integer power (`(±2^k)^n = ±2^(kn)`) is now computed directly, including the overflow, subnormal and underflow ends of the range (it used to return `PrecisionRetryExhausted` under the directed modes for `2^(2^31-1)`, where MPFR returns the mode's overflow value). Every other function with exact results (`ln`, `log2`, `log10`, `log`, `exp`, the trigonometric, inverse and hyperbolic functions at 0 and 1, `nth_root`, `hypot`, `fma`, `reciprocal`) is covered by `exact_results_are_decided_not_retried` against MPFR (all modes, three precisions, up to 4096 bits) and passes with the new check.

Digests: every pinned digest is unchanged except `jacobi_all_functions_digest`, which is re-pinned (13_695_875_014_094_095_401 to 6_259_183_207_121_670_912). It differs in exactly the 46 entries described above (23 function/m pairs at `p` = 64 and 128, all `ToEven`, operand 21): the old value of each was `NaN(PrecisionRetryExhausted)`, the new value is `u` for `am`, `sn`, `sc`, `sd` and 1 for `cn`, `dn`, `cd`, `nc`, `nd`, `dc`. Oracle: the series in `u` (`f(u|m) = t(1 + O(u^2))`, deviation below 2^-(2^32), far under half an ulp) with the constants checked by mpmath at a surrogate argument (`u' = 3·2^-300`, 2000 bits: every deviation is below 8·u'^2), and the new unit test `jacobi_of_a_tiny_argument_rounds_to_the_leading_term`. Not changed: the unit tests, the 80 MPFR tests and the oracle suites pass with the new check.

### A subnormal operand is no longer cloned by `add`, `sub`, `mul`, `div` and `sqr` (internal)

`ExactNumNumber::normalize` used to clone the mantissa of a subnormal operand and shift the clone left before every sum, product and quotient (and so before `add_ref`, `sub_ref`, `mul_ref`, `div_ref`, `sqr` and the owned forms). The kernels now take the operands as stored: the exact sum, difference and product are renormalized (`maximize`) before they are rounded to `p`, and `Mantissa::div` shifts a subnormal dividend while it writes it into its padded working buffer, which exists for every dividend. A subnormal *divisor* is still shifted into a working buffer of its own, because the division kernel needs a normalized divisor; that buffer is written by a shift, it is not a `Mantissa::clone`, and a normal divisor needs none. `rem` still uses `normalize` (and has other copies). Results are bit-identical: every pinned digest is unchanged, and a new MPFR grid (`tests/mpfr/subnormal_ops.rs`: subnormal and tiny operands on either side and of either sign, `add`/`sub`/`mul`/`div` in the borrowed and owned forms and `sqr`, five rounding modes, three result precisions) agrees with MPFR bit for bit. A unit test counts `Mantissa::clone` calls and requires zero for all of these forms.

### `ToOdd` is round to nearest, ties to odd (documentation)

`RoundingMode::ToOdd` always was round to nearest with ties to the odd neighbour, the mirror of `ToEven`. It is not the sticky "round-to-odd" (jamming) mode of the literature, in which every inexact result gets its last bit set to 1; the documentation (`RoundingMode`, README, HELP, the capabilities document) said "round half to odd" and could be read as that mode. It now says so in plain words, with examples that a test pins: at `p = 4`, 17/16 is an exact tie (`ToEven` 1.000, `ToOdd` 1.001) and 33/32 is below the tie (1.000 under both; sticky round-to-odd would give 1.001). No behaviour changed. The sticky mode is used inside the library (`exact_bits`, the division and addition kernels) but is not a public rounding mode; see `/tmp/zf-logs/sticky-round-to-odd-option.md` for the cost of exposing one.

### Tests

- `tests/oracle/goldens/quad.jsonl`: 12 `gk-eps` rows (deep tier) for the Wynn-ε Gauss–Kronrod path on endpoint singularities: 6 integrands (\(x^{-0.75}\), \(x^{-0.9}\), \(x^{-0.99}\), and \(\ln^k x/\sqrt{x}\) for \(k=1,2,3\)) × 256 and 512 bits. 18 `gk-inf` rows cover `integrate_adaptive_gk_infinite`. `abs_err_est` on those rows is checked against the oracle.
- `tests/oracle/goldens/specials/elliptic.jsonl`: \(m=0\) shortcut goldens for \(F\), \(E\), and \(\Pi\) in all seven rounding modes, and angle-form goldens for `elliptic_f_phi`, `elliptic_e_phi`, and `elliptic_pi_phi`. `tests/elliptic_phi.rs` covers the period reduction, \(\varphi=\pm\pi/2\), and the `NaN(InvalidArgument)` cases.
- `carry_tests` in `common/util.rs`: `add_carry` and `sub_borrow` against a widening double-word reference (`u128` for 64-bit words, `u64` for 32-bit words), with carry-in 0 and 1 edges, random limbs, and 8-limb all-`WORD_MAX` chains.
- `tests/oracle/goldens/ode.jsonl` and `balls.jsonl`: 24 ODE and 24 ball regression goldens, recorded before the quadrature change.
- `scripts/ci_no_intrinsics.sh` uses `git grep -nE` and fails hard when `git` is missing or the search itself errors, so a broken search cannot pass as "no match".
- Audit tests, each red before its fix: `tests/rounding_mode_docs.rs`, `tests/from_raw_parts_checks.rs`, `tests/overflow_modes.rs`, `tests/inexact_per_op.rs`, `tests/exponent_range.rs` (red on i686), and in `tests/mpfr/fuzz_round_modes.rs` `fuzz_mpfr_overflow_all_round_modes`, `fuzz_mpfr_exact_bit_precision`, and the widened `fuzz_mpfr_bit_exact_all_round_modes`.
- `.github/workflows/arm-tests.yml`: a test-only ARM64 workflow (`ubuntu-24.04-arm`, Rust 1.99) that runs the workspace tests and the oracle test. It does not publish.

## 1.0.12 — 2026-10-06

Attribution for code derived from [astro-float](https://github.com/stencillogic/astro-float) (MIT, Copyright (c) 2022 stencillogic; crates `astro-float` 0.9.6, `astro-float-num` 0.3.7, `astro-float-macro` 0.4.6, commit `c0d9a26`). `NOTICE` has the full MIT copyright and permission notice. Each derived file names that origin. Original work stays under Apache-2.0 and is also offered as MIT. The file list is `doc/PROVENANCE.md`.

### Heun closed-form detection

`heun_g_is_2f1` and `heun_c_is_1f1` no longer call `values_agree` (about eight bits of slack). Exact dyadic parameters are compared with full-precision arithmetic. Any other parameter must agree to one ulp of `p`, with the residual measured at precision `2p`. A near-miss such as \(q=-20+2^{-120}\) on the terminating \({}_2F_1(-10,1;1;1/2)\) stays on the general Heun series.

Golds are mpmath 1.4.1 `hyp2f1` / `hyp1f1` for the exact reductions, and an mpmath power series for the near-misses. Those near-misses disagree with the closed form.

### Known issues

None are open.

- Complex `bessel_k` for non-integer \(\nu\) inside the series regime: closed in 1.0.11 (DLMF 10.27.4).
- Real `bessel_j_nu` / `bessel_y` for \(\lvert x\rvert\gtrsim 10^4\): closed in 1.0.6 (Hankel expansion).
- Complex Bessel on the negative real axis: closed in 1.0.6 (reflection formulas).

`doc/AUDIT.md` still describes those three as the 1.0.5 audit's open list. That text is the record of 1.0.5, not a list of current defects. A Heun argument whose real segment from 0 meets a singular point is `NaN(InvalidArgument)`, which is the documented domain.

## 1.0.11 — 2026-10-05

Closes the three items that were still listed under Known issues at 1.0.10. No public signature was removed or changed.

### Complex \(K_\nu\) in the series regime

Non-integer \(K_\nu(z)\) with \(\operatorname{Re} z>0\) and \(\lvert z\rvert<0.35(p+112)\) is DLMF 10.27.4,

\[
K_\nu(z)=\frac{\pi}{2\sin(\nu\pi)}\bigl(I_{-\nu}(z)-I_\nu(z)\bigr),
\]

with the guard sized to \(e^{2\operatorname{Re} z}/\ln 2\) plus the loss in \(\sin(\nu\pi)\). The previous path built \(K\) from \(J\pm iY\), which cancels about \(3\lvert z\rvert\) bits and runs both \(Y\) series. Integer \(\nu\) still uses the \(Y\) recurrence. Temme \({}_2F_0\) is unchanged for large \(\lvert z\rvert\).

### Meijer \(G\) at coincident poles

A residue series whose numerator \(\Gamma\) hits a pole (integer spacing in the active \(b_1,\ldots,b_m\) or \(a_1,\ldots,a_n\)) takes the hypercomb limit: every parameter is shifted by a distinct \(O(2^{-\mathrm{hmag}})\) amount and the value is kept only when two step sizes agree. That is the finite part of a logarithmic (or higher) residue. A non-removable pole stays `NaN(InvalidArgument)`.

Golds (mpmath `meijerg`, which is the same limit): \(G^{2,0}_{0,2}(z\mid b,b)=2z^b K_0(2\sqrt{z})\) and \(G^{1,2}_{2,2}(z\mid 1,1;1;0)=\ln(1+z)\).

### Local and confluent Heun past the unit disk

Inside the disk of convergence the Frobenius series is unchanged. Outside it, the value is the real-analytic continuation along the component of the cut plane that contains 0:

- Local Heun: singular points \(1\) and \(a\). Taylor steps of \(z(z-1)(z-a)y''+\cdots=0\), each at most half the distance to \(\{0,1,a\}\).
- Confluent Heun: singular point \(1\). The same stepping on \((-\infty,1)\). \(z\ge 1\) is `NaN(InvalidArgument)`.
- \(q=a\alpha\beta\) and \(\delta=\alpha+\beta-\gamma+1\) reduces to \({}_2F_1(\alpha,\beta;\gamma;z)\).
- \(\delta=0\), \(\varepsilon=-1\), \(q=\alpha\) reduces to \({}_1F_1(-\alpha;\gamma;z)\), entire.
- Two consecutive series coefficients that are exact zeros are a polynomial, returned for every finite \(z\).
- More than 8192 Taylor steps is `NaN(InvalidArgument)`.

Biconfluent, double-confluent, and triconfluent Heun stay unimplemented: there is no separate public signature, and the \(z=0\) series is not the regular solution. That is not a leftover of the \(\lvert z\rvert<1\) limit.

Golds: mpmath `hyp2f1` / `hyp1f1` for the reductions (including \(\lvert z\rvert>1\)), and an independent mpmath ODE Taylor reference for non-reducing accessories at \(z=-1.25\), \(z=0.7\) (local, past \(\lvert a\rvert\)), and \(z=-1.7\) (confluent).

## 1.0.10 — 2026-10-05

### Gauss–Kronrod (7, 15) and Wynn ε

SoftFloat quadrature for `∫_a^b f(x) dx` on a finite interval. No `f64` path.

- `kronrod_pair` builds the (7, 15) nodes and both weight rows at `p + WORD_BIT_SIZE`. The recurrence is Laurie’s Jacobi–Kronrod matrix for the Legendre weight (*Math. Comp.* 66, 1997); weights are Golub–Welsch. Orders other than 7 return `None`.
- `gauss_kronrod_interval` returns the Kronrod sum `K`, the embedded Gauss sum `G`, and `abs_err_est = |K − G|`.
- `wynn_epsilon` is Wynn’s ε table on a sequence of partial sums.
- `integrate_adaptive_gk(f, a, b, tol_rel, max_subintervals, p, rm, cc)` bisects the panel with the largest `|K − G|` and, when that raw sum is still slow, accepts a stable ε value. The absolute stop is `max(tol_rel * max(1, |estimate|), 2^(1−p))`. Failure is [`QuadratureError`], not a panic and not a quiet wrong number.

The local error `|K − G|` and the ε acceleration for endpoint singularities are the QUADPACK design of Piessens, de Doncker-Kapenga, Überhuber, and Kahaner (*QUADPACK*, Springer, 1983; SLATEC’s port of that design is public domain). This tree recomputes the rule in SoftFloat `ExactNum`. It does not copy the ACM / SLATEC Fortran.

Still unsupported: infinite intervals (`QAGI`), oscillatory weights (`QAWO` / `QAWF`), algebraic endpoint weights (`QAWS`), Cauchy principal values (`QAWC`), Kronrod pairs other than (7, 15), and the full `QAGS` heap.

## 1.0.9 — 2026-10-04

The complex \(I_\nu\) power-series guard is sized to cancellation, not \(1.5\lvert z\rvert+16\) (that blanket is what 1.0.6 through 1.0.8 shipped; 1.0.7 and 1.0.8 were not this fix). When every term is positive (real \(z>0\), real \(\nu>-1\)) the guard is a fixed 32-bit pad. Otherwise it is that pad plus an upper bound on \(\lceil(\lvert z\rvert-\mathrm{Re}\,z)/\ln 2\rceil\), raised again if the summed peak still exceeds the result. Real \(I_\nu\) for large \(\lvert x\rvert\) now runs past the peak at \(k\approx\lvert x\rvert/2\). Real \(K_\nu\) uses the asymptotic series when the \(2.885\lvert x\rvert\) guard would exceed `BESSEL_GUARD_MAX`, instead of returning `InvalidArgument`. Complex \(K\) on the positive real axis uses that real kernel past the series switch. Clippy `-D warnings` cleanups in the same tree (manual `find`, `clamp`, `from_ref`, single-iteration test loops) do not change numeric results. The 1.85 `unused_unsafe` allow on `_addcarry_u64` / `_subborrow_u64` stays.

## 1.0.8 — 2026-10-03

README upgrade target corrected to 1.0.8; no math.

## 1.0.7 — 2026-10-02

1.0.7 keeps the 1.85-required unsafe wrappers on `_addcarry_u64`/`_subborrow_u64` and allows `unused_unsafe` so the crate builds on rustc 1.93+ (those intrinsics became safe in 1.93, not 1.87).

(Superseded in 2.0.0: intrinsics and unsafe removed; see Forensic audit.)

## 1.0.6 — 2026-10-02

Additive SoftFloat pack: Jacobi \(P_n^{(\alpha,\beta)}\), real large-\(|x|\) Bessel Hankel, complex Bessel cut/speed fixes, and the catalog specials Accumath was missing. No public signature was removed or changed — **1.0.6, not 2.0.0**. Jeff publishes crates.io — this tag is not published from the PR.

### Version path

Stayed on **1.0.6**. Every new entry is an added `ExactNum` / `ExactComplex` method (or a bugfix of an existing one). No existing signature changed, so a breaking 2.0.0 bump is not required.

### `ExactNum::jacobi_p`

- Jacobi \(P_n^{(\alpha,\beta)}(x)\) via the Bonnet three-term recurrence, same house style as `hermite_h` / `chebyshev_t` in `orthopoly.rs`. Cap `ORTHOPOLY_N_MAX`.
- Golds: \(\alpha=\beta=0\) reduces to Legendre \(P_n\); \(P_2^{(1,1)}(1/2)=3/16\); \(P_n^{(\alpha,\beta)}(-x)=(-1)^n P_n^{(\beta,\alpha)}(x)\).
- Independent mpmath 1.4.1 references for non-integer \((\alpha,\beta)\), including \(\alpha=1/3\), \(\beta=2/3\).

### Real Bessel large \(|x|\)

- `bessel_j` / `bessel_j_nu` / `bessel_y` use the Hankel expansion (DLMF 10.17, \(P,Q\) stopped at the smallest term) when \(|x|\ge\max(16,\,0.35(p+112))\) and \(|\nu|^2<2|x|\), and whenever the series cancellation guard would exceed `BESSEL_GUARD_MAX` under that order condition.
- Closes the 1.0.5 known issue: real \(J_0(20000)\), \(Y_0(20000)\) and the rest of that class no longer return `NaN(InvalidArgument)`.
- Negative \(x\) for integer order uses \(J_n(-x)=(-1)^n J_n(x)\). \(Y\) remains defined for \(x>0\).
- Regression tests vs mpmath at 128 and 256 bits, plus a Wronskian gold at \(x=20000\).

### Complex Bessel (negative real axis + \(K\) speed)

- On the cut \(\operatorname{Im} z=0\), \(\operatorname{Re} z<0\) with real \(\nu\), \(J,Y,I,K\) are assembled from the real kernels (DLMF 10.11 / 10.34) so algebraically zero parts stay exact zeros. Closes the 1.0.5 hang / `PrecisionRetryExhausted` on e.g. `bessel_i(-2, 1/2)`.
- \(I_\nu(z)\) for \(\mathrm{Re}\,z>0\) uses the well-conditioned \({}_0F_1\) series instead of a rotation through \(J\).
- \(K_\nu(z)\) for \(\mathrm{Re}\,z>0\) past the series/Hankel switch uses Temme \({}_2F_0\) (DLMF 10.32.10, modified Lentz) instead of \(J\pm iY\) with \(3|z|\) guard bits. The Hankel / series path remains the fallback. High-precision non-integer \(K\) just inside the series regime can still be slow; the Temme path covers the large-\(|z|\) cases that were the 74 s class.

### Catalog specials (additive)

New `ExactNum` methods, SoftFloat-only, with identity and mpmath golds:

- Scorer \(\mathrm{Gi},\mathrm{Hi}\) (`scorer_gi`, `scorer_hi`); \(\mathrm{Hi}=\mathrm{Bi}-\mathrm{Gi}\).
- Kelvin \(\mathrm{ber},\mathrm{bei},\mathrm{ker},\mathrm{kei}\) (order 0 and `_nu`).
- Struve \(\mathbf{H}_\nu\), Anger \(\mathbf{J}_\nu\), Weber \(\mathbf{E}_\nu\).
- Clausen \(\mathrm{Cl}_2,\mathrm{Cl}_3\); Barnes \(G\); Hurwitz / Riemann \(\zeta\) (integer \(s\ge 2\)); polygamma \(\psi^{(n)}\).
- Inverse Jacobi `jacobi_arcsn` / `arccn` / `arcdn`; \({}_0F_1\), \({}_1F_1\), \({}_pF_q\); Lambert \(W_0,W_{-1}\); \(\mathrm{Li}_n\) for \(|x|\le 1\).

`expr!` leaves for the one- and two-argument catalog names (`scorer_gi`, `kelvin_ber`, `struve_h`, `clausen_cl2`, `barnes_g`, `lambert_w0`, …).

### Meijer \(G\) / Fox \(H\) / Heun (additive, still 1.0.6)

Leftovers after the catalog pack — Accumath still had local SoftFloat clones for these. No existing signature changed.

- `ExactNum::meijer_g(an, ap, bm, bq)` — DLMF 16.17 residue + \({}_pF_q\). Type-1 or type-2 series by \((p,q,\lvert z\rvert)\). Coincident poles / numerator \(\Gamma\) poles → `NaN` (no logarithmic `hypercomb` limit). Negative \(z\) → `NaN` unless every used power is an integer.
- `ExactNum::fox_h` — pairs \((a,A)\), \((b,B)\). \(A=B=1\) is Meijer \(G\); positive rational scales use the Gauss multiplication formula (mpmath `foxh`). Non-rational or large scales → `NaN`.
- `ExactNum::heun_g(a,q,α,β,γ,δ)` — local Heun \(\mathrm{Hl}\), DLMF 31.3, \(\lvert z\rvert<1\).
- `ExactNum::heun_c(α,γ,δ,ε,q)` — confluent Heun, \(y(0)=1\) Frobenius solution of DLMF 31.12.1, \(\lvert z\rvert<1\). Not Maple `HeunC(α,β,γ,δ,η,z)`.
- Left out (no Accumath SoftFloat evidence, or no unique \(z=0\) series): Mathieu, Lamé, Appell \(F_{1..4}\), HeunB / HeunD / HeunT, complex \(G/H\)/Heun. `expr!` skipped (too many arguments, same as `hypergeom_pfq`).
- Golds: mpmath `meijerg` / `foxh` at 50 digits (including \(J_{1/4}\), series-2 \(\lvert z\rvert>1\), rational Fox scale); HeunG vs mpmath \({}_2F_1\); HeunC vs mpmath \({}_1F_1\).

## 1.0.5 — 2026-09-27

Bug-fix release: special-function accuracy. No public API changes (no new public items, no signature changes); results that were wrong now match mpmath to working precision. Every fixed case below is covered by a test with a 50-digit (or 2p+64-bit) mpmath reference.

### `hypergeom_2f1` (real)

- **What was wrong.** For \(1/2<z<1\) the Pfaff transform maps to \(\lvert w\rvert>1\), and the divergent series was summed there and silently truncated at the term cap. Example: \({}_2F_1(1,1;2;0.7)\) returned `6.06e175` (true `1.71996…` = \(-\ln(0.3)/0.7\)). Near \(z=-1\) the series was truncated too. For \(z\le-1\) the A&S 15.3.7 path returned NaN whenever \(b-a\) is an integer (e.g. \({}_2F_1(0.5,1.5;2;-3)\)). The series stop was absolute, not relative.
- **What changed.** New dispatch:
  - poles in \(c\) → NaN;
  - terminating (polynomial) series;
  - Gauss at \(z=1\) when \(c-a-b>0\);
  - \(z>1\) → NaN;
  - Pfaff for \(z<-1/2\);
  - direct series for \(\lvert z\rvert\le1/2\);
  - \(z\to1-z\) connection for \(1/2<z<1\) (DLMF 15.8.4), with the A&S 15.3.10–15.3.11 logarithmic case for integer \(c-a-b\).

  Each path measures cancellation and re-evaluates with that many guard bits. The series returns an error at its term cap instead of truncating.

### `betainc`

- **What was wrong.** The integer-\(b\) polynomial and the series were used at large \(a,b\), where they cancel catastrophically. Examples: \(I_{0.7}(30.5,2.25)\) returned `2.68e240` (true `3.0312…e-4`), \(I_{0.3}(1000,1000)\) returned `1.47e124`, and \(I_{0.5}(5000,5000)\) returned `1.8e2284`.
- **What changed.** A continued fraction (DLMF 8.17.22, modified Lentz), with the symmetry \(x>(a+1)/(a+b+2)\) applied. The polynomial is used only when its measured loss is ≤ 64 bits. \(I_{1/2}(a,a)=1/2\) is returned exactly.

### Other real special functions

- **`erf` / `erfc` / `normal_cdf`.** For every \(\lvert x\rvert\ge4\) (binary exponent ≥ 3), 1.0.2–1.0.4 summed the divergent asymptotic series past its smallest term, so the tail was garbage at every precision. At 256 bits: `erfc(4.2)` = `1.98e349` (true `2.855e-9`), `erfc(10)` = `4.99e17` (true `2.088e-45`), `erfc(12.5)` = `7.69e-56` (true `6.232e-70`, 1.2e14 times too large), `erf(4)` = `-4.40e366`, `erf(-4.2)` = `1.98e349`, `normal_cdf(-6)` = `1.53e433`. At 128 bits `erf(4)` = `-1.22e199` and `erfc(10)` = `1.3e-18`. 1.0.2 and 1.0.4 give identical wrong values (the code did not change between them). Now the asymptotic series is used only when \(x^2\ge0.7p+4\), where its smallest term is below \(2^{-p}\), and it stops at its smallest term. Otherwise the Maclaurin series runs with \(1.5x^2\) guard bits, and `erfc` for \(x>0\) carries another \(1.5x^2\) bits for the cancellation in \(1-\operatorname{erf}\). `normal_cdf` uses `erfc(-(x-μ)/(σ√2))/2`, so the lower tail keeps full relative accuracy (`normal_cdf(-40)` was garbage). A sweep of 561 points (x from −30 to 30, step 0.1 on \(3.5\le\lvert x\rvert\le30\)) for `erf`, `erfc` and `normal_cdf` at 64, 128, 192, 256, 512 and 1024 bits is within 2 bits of mpmath everywhere. 1.0.4 misses at 619 of these points at 128 bits.
- **`gamma` / `digamma`.** The Stirling start was fixed, which capped accuracy at about 274 bits (Γ) and 530 bits (ψ) at any precision. The shift is now precision-dependent. `gamma(11.5)` at 1024 bits had 276 correct bits.
- **Bernoulli numbers (`gamma`, `ln_gamma`, `digamma`, real and complex).** The Stirling coefficients came from the Akiyama–Tanigawa recurrence in floating point, which cancels catastrophically. With 64-bit limbs the Ziv loop usually hid this. With 32-bit limbs (`thumbv7em-none-eabihf`, i686) Γ at an internal 192 bits had about 101 correct bits, so `bessel_i(2.25, 1.5)` at 128 bits returned only 101 correct bits. The \(B_{2k}\) are now built from exact tangent numbers (Brent–Harvey) with one rounding each. They are computed once per evaluation instead of once per term (the old cost was cubic in the number of terms).
- **`ln_gamma`.** `ln_gamma(1e20)` and `ln_gamma(1e300)` returned `inf`; large arguments now use Stirling directly.
- **Euler–Mascheroni γ.** The constant (used by `ei`, `ci`, `li`, `bessel_y`, `bessel_k`, `digamma`) was correct to only about 383 bits at any precision. It is now computed with Brent–McMillan.
- **`ei`.** Series cancellation for \(x<0\) was unguarded: `ei(-50)` had 115 of 128 bits and `ei(-200)` was garbage. Guard bits are now sized from the measured loss.
- **`fresnel_s` / `fresnel_c`.** Series cancellation was unguarded (`fresnel_s(10)` had 36 bits, `fresnel_c(12)` was garbage), and the asymptotic switch now depends on precision.
- **`ai` / `bi` / `ai_prime` / `bi_prime`.**
  - The derivative asymptotic coefficients had the wrong sign/form, so `ai_prime(8)` and `bi_prime(-30)` had 7 correct bits.
  - The series/asymptotic switch now depends on precision: `ai(12)` had 85 bits, and `ai(-30)` had 320 bits at p=512.
- **Real Bessel.**
  - `bessel_j_nu` / `bessel_y` / the integer series had no cancellation guard, so J/Y at \(x=200\) were garbage. They now carry \(1.5\lvert x\rvert\) guard bits. For \(\lvert x\rvert\gtrsim10^4\) they return `InvalidArgument` rather than garbage (no large-argument Hankel expansion yet).
  - `bessel_k` used its asymptotic series below the precision where it converges (`bessel_k(20, 0)` had 62 bits).
  - Series stops are now relative: `bessel_j_nu(1e-10, 10.5)` had 276 of 512 bits.
- **Series stopping rule.** The absolute stopping rule was replaced by a relative one in `gammainc`, `si`, Fresnel and the Bessel series.
- **`elliptic_f(±1 | 1)`** returned `NaN(InvalidArgument)` even though \(\lvert x\rvert\le1\) is the documented domain. It now returns \(\pm\infty\), consistent with `elliptic_k(1)` = \(+\infty\).
- **Carlson \(R_C\), `elliptic_pi`, `elliptic_pi_complete`.** The \(x\approx y\) shortcut used an absolute threshold of \(2^{-(p/3+16)}\), which limited Π to about 340 bits at p=512. The threshold is now relative. The duplication cap scales with precision (`max(128, p/4+64)`), for real and complex.
- **Exactly representable results.** Ziv's loop could not certify these and returned `PrecisionRetryExhausted`: `log(100, 10)`, `log(2, 4)`, `log(0.25, 2)`, `nth_root(1e10, 5)`, `jacobi_cd(u, 1)`, `jacobi_dc(u, 1)`. They are now detected and returned exactly.

### Complex (`ExactComplex`)

- **`erf` / `erfc` / `fresnel_s` / `fresnel_c`.** The power series had no cancellation guard. At 512 bits, `erf(-15-5i)` was wrong from the 45th digit (true value is \(-1+4.4\times10^{-89}+2.3\times10^{-89}i\)) and `fresnel_s(-15-5i)` had 200 bits. The series now carries \(1.5\lvert u\rvert^2\) guard bits.
- **`ci`.** Wrong overall sign/branch: `ci(0.5+0.5i)` and `ci(-2+i)` had 0 correct bits. It is now principal and matches mpmath.
- **`si(i)`.** Returned NaN.
- **`ei`.** The \(i\pi\,\mathrm{sgn}(\operatorname{Im} z)\) term was missing from the asymptotic path, and the switch point was fixed rather than precision-dependent (`ei(20+20i)` had 22 bits).
- **`acosh` / `atanh` on their cuts.** `acosh(-2+0i)` had a negative real part and `acosh(-1+0i)` returned \(-i\pi\). `atanh(±1+0i)` was NaN; it is now \(\pm\infty+0i\).
- **`ai` / `bi`.**
  - The phases in the `bi` connection formula were swapped: `bi(0.001+10i)` had 0 bits.
  - `bi` on the Stokes line used only the dominant asymptotic series (`bi(20+20i)` had 110 bits).
  - The switch point now depends on precision (`ai(0.001+10i)` had 63 bits).
- **`bessel_*`.**
  - The Hankel P/Q sums ran past their smallest term, giving garbage at moderate \(\lvert z\rvert\).
  - `bessel_i` / `bessel_k` used the wrong rotation for \(\operatorname{Re} z<0\): `bessel_i(-2+i, 2+0.5i)` and `bessel_k(-2+i, 0)` were wrong. They now use DLMF 10.27.6 / 10.27.8.
  - The series now has cancellation guards.
  - `bessel_j_nu` / `bessel_y` with \(\operatorname{Re} z<0\) in the Hankel regime (\(\lvert z\rvert\gtrsim 0.35(p+112)\)) summed the large-argument expansion outside its sector \(\lvert\arg z\rvert<\pi\): `bessel_j_nu(-20000, 0)`, `bessel_j_nu(-500+i, 1)` and `bessel_y(-300+40i, 2.5)` had 0 correct bits (1.0.4 too). They now continue from \(-z\) (DLMF 10.11.1–10.11.2).
  - `bessel_i` / `bessel_k` at real \(z>0\) with real non-integer \(\nu\) never converged: the imaginary part is exactly 0 but was computed as a rounding residue that Ziv's loop cannot certify. 1.0.4 returned `PrecisionRetryExhausted` after about 3 s for `bessel_i(2, 0.5)`; the first 1.0.5 candidate took about 5 minutes and returned NaN. The imaginary part is now an exact 0 there.
  - `bessel_k` on its cut \((-\infty,0)\) with \(\operatorname{Im} z=0\) returned the value from below; it now returns the value from above, like `bessel_i`, `bessel_j_nu` and mpmath.
- **`gamma` / `ln_gamma` / `digamma`.** Precision was capped at about 275 bits. `ln_gamma` now returns the documented principal \(\ln(\Gamma(z))\): `ln_gamma(3-4i)` returned a different branch.
- **`hypergeom_2f1`.**
  - The documented 10 000-term cap was actually \(p+96\), and the series stop was absolute. `2F1(0.5,0.5;1;0.999999)` had 1 correct bit.
  - All-real arguments with \(z<1\) now use the real algorithm.
  - Integer \(m=c-a-b\) with \(z\) near 1 (\(\lvert 1-z\rvert\le 1/2\), or \(\lvert 1-z\rvert<1\) with \(\lvert z\rvert\ge1\)) returned NaN: the direct series converges too slowly there (or not at all for \(\lvert z\rvert\ge1\)) and the generic \(1-z\) transform has \(\Gamma(\pm m)\) poles. Such arguments now use the logarithmic connection (DLMF 15.8.10, A&S 15.3.10–15.3.11) for \(m\ge0\), after Euler's transformation \((1-z)^m F(c-a,c-b;c;z)\) for \(m<0\); cancellation is measured and the evaluation repeated with more bits. Example: \({}_2F_1(0.3+0.2i,\,0.7;\,2+0.2i;\,1.2+0.3i)\).
  - \(c-a-b\) within \(\delta\) of an integer (the generic \(1-z\) transform loses about \(\log_2(1/\delta)\) bits) now carries that many guard bits; \(z\) close to 1 inside the unit disc with non-integer \(c-a-b\) uses the transform instead of a slowly converging series.
  - A grid of 504 cases (9 \((a,b)\) pairs, \(m\in\{-3,\ldots,3\}\) plus near-integer \(c\), 9 values of \(z\) around 1 including \(\lvert z\rvert>1\)) agrees with mpmath to within 2 bits at 64, 128, 256 and 512 bits.

### `no_std`: built and tested

- The library already built without `std`. Now the reference tests also run without it.
- `nostd-tests/` is a `#![no_std]` crate (core + alloc only; not published; excluded from the workspace and the package) with 77 mpmath reference cases covering:
  - elementary functions, erf/erfc (including the tail), Γ/lnΓ/ψ at 128–320 bits, 2F1, betainc, Ei/Si/Ci, Fresnel, Airy, Bessel, elliptic, Jacobi, `normal_cdf`;
  - complex Γ/erf/Ci/Ai/J/I/2F1 (including the integer-\(c-a-b\) log case near \(z=1\), \(J\) with \(\operatorname{Re} z<0\) and \(I\) on the real axis);
  - quadrature, Brent and RK45.
- Two runners:
  - `nostd-host`: a `#![no_std]`/`#![no_main]` x86_64 Linux binary linking only libc;
  - `nostd-qemu`: bare-metal `thumbv7em-none-eabihf` (Cortex-M4F, 32-bit limbs) under `qemu-system-arm -machine mps2-an386` with semihosting.
- `scripts/ci_nostd.sh` builds `zenith-float` and `zenith-float-num` with `--no-default-features --target thumbv7em-none-eabihf` and runs both runners. Each passes 77/77.
- There is no separate `alloc` feature; an allocator is always required.

### Comparison, and the numeric methods that depend on it

- **`ExactNum::cmp` / `ExactNum::abs_cmp` return exactly `-1`, `0` or `1`.** The documented contract (positive / negative / zero) is unchanged.
  - For finite values with equal exponents, 1.0.4 returned the raw mantissa word difference: `3.cmp(2)` was `Some(2^62)`.
  - The crate's own numeric routines test `cmp(..) == Some(-1)` / `== Some(1)`, so in 1.0.4:
    - `gauss_legendre`, `tanh_sinh`, `bisect`, `illinois`, `brent`, `rk4`, `euler`, `rk45_adaptive`, `chebyshev_coeffs` and `chebyshev_eval` rejected (`None` / NaN) any interval whose endpoints share a binary exponent, e.g. `[2, 3]`;
    - `rk45_adaptive` with a non-zero `rtol` never accepted a step and returned `None`;
    - several convergence tests misfired.
- **`ExactNum::abs_cmp` compared signed values** for finite operands, contradicting its documentation: `(-2).abs_cmp(1)` was negative. It now compares magnitudes.
- **`brent` rewritten** as Brent's classic zero finder (bracket `[b, c]`, inverse quadratic / secant with a bisection fallback). The old loop had no bisection safeguard: it stalled like regula falsi and returned `None` with `root_default_tol` at 128 bits.
- **Resolution stop in `bisect` / `illinois` / `brent`.** When `tol` is below the working resolution (e.g. `root_default_tol` = \(2^{-256}\) at \(p<256\)), they now return the bracketed root instead of `None` after `ROOT_MAX_ITER`.

### Tests

- `ops::special` unit tests:
  - `test_hypergeom_2f1_mpmath`: 41 values across \(z\in[-10,1]\), including the log case, \(z=1\) and near-1 points, plus 4 error cases;
  - `test_betainc_mpmath`: 17 values;
  - `test_gamma_digamma_high_precision`: 1024 bits.
- New integration test `tests/special_audit_mpmath.rs`: 61 regression points (real and complex, 128/512/1024 bits). Each must agree with mpmath to within 2 ulp; 53 of them fail on 1.0.4.
- Reference generators: `refs/` in the maintainer workspace (mpmath 1.4.1, 50 digits or 2p+64 bits).
- New integration test `tests/numeric_methods_regress.rs`: 5 tests covering `cmp` normalization, `[2, 3]` intervals for quadrature, roots, ODE and Chebyshev, `rk45_adaptive` with `rtol`, and `brent` / `bisect` / `illinois` with `root_default_tol` at 128 bits. All 5 fail on 1.0.4.
- New integration test `tests/erf_tail_mpmath.rs`: 221 `erf` / `erfc` / `normal_cdf` cases with 50-digit mpmath references. They cover x in [3.5, 30] and negative x at 64/128/256 bits, `erfc` at 1024 bits, and the `normal_cdf` lower tail. A second test checks erf odd symmetry and erfc reflection. On 1.0.4, 109 of the 221 cases fail, and so does the symmetry test.
- New integration test `tests/edge_values.rs`: `elliptic_f(±1 | 1)` = ±∞.
- `no_std` harness `nostd-tests/` (77 cases) and `scripts/ci_nostd.sh` (see above).
- New integration test `tests/c2f1_log_mpmath.rs`: 54 complex \({}_2F_1\) cases near \(z=1\) with integer \(c-a-b\in\{-3,\ldots,3\}\) and near-integer \(c-a-b\), at 128 and 256 bits, each to \(p-4\) bits. The first 1.0.5 candidate returned NaN for 34 of them.
- New integration test `tests/complex_bessel_mpmath.rs`: 218 cases (J/Y with \(\operatorname{Re} z<0\) in the Hankel regime, I/K on the positive real axis with real non-integer \(\nu\), K on its cut), at 128 and 256 bits, each to \(p-4\) bits norm-wise.
- 32-bit targets (i686): two tests assumed 64-bit limbs; the library was correct. `ext::test_ext` asserted that a 128-bit value is stored inline, which holds only with 64-bit limbs (4 limbs > `INLINE_WORDS` = 2 on 32-bit); the assertion is now per limb width. The `add_commutes` property test built values at \(p=32\) (one 32-bit word), below the documented \(p\ge64\) minimum of `from_i64`, and got `NaN(InvalidArgument)`; it now uses \(p\) in multiples of 64 bits. `cargo test -p zenith-float-num --target i686-unknown-linux-musl` now passes.
- `cargo fmt --all` applied (1.0.4 was not `fmt --check` clean). Two clippy lints in the MPFR differential tests (`--all-features`) were fixed.

### Known issues (fix planned)

- **Complex `bessel_k` is slow for non-integer \(\nu\) at high precision when \(\lvert z\rvert\) is just inside the series regime** (\(\lvert z\rvert<0.35(p+112)\)). Still open; queued for 1.0.7.
- **Real `bessel_j_nu` / `bessel_y` return `NaN(InvalidArgument)` for \(\lvert x\rvert\gtrsim10^4\)**: closed in 1.0.6 (real Hankel expansion).
- **Complex Bessel on the negative real axis** (\(\operatorname{Im} z=0\), \(\operatorname{Re} z<0\)) with real \(\nu\): still open; queued for 1.0.7.

## 1.0.4 — 2026-09-20

Clippy debt clear under `cargo clippy --workspace --all-targets -- -D warnings`. SoftFloat IEEE inherent ops renamed to `soft_*` where they collided with `std::ops`; needless range loops / doc list / `?` tidy; workspace `clippy.toml` thresholds for SoftFloat kernels.

## 1.0.3 — 2026-09-19

Coordinated patch with latex-rust, hdf5-rust, and redb-view (pure-Rust FOSS family adjacent to Accumath; Accumath itself stays proprietary).

## 1.0.2 — 2026-09-04

This crate is a numeric library. CSV and the 16-byte binary record remain.

- Removed HDF5 I/O and the `hdf5` Cargo feature. This crate does not read or write that file format.
- Build sheets, checklists, and walk lists are gone. Capabilities, contributing, getting started, help, library, and the other user guides stay.

## 1.0.1 — 2026-09-03

Jacobi elliptic functions on the real line.

- `ExactNum` Jacobi family: `am`, `sn`, `cn`, `dn`, and the nine quotients `cd` `ns` `nc` `nd` `sc` `sd` `cs` `ds` `dc`. Parameter \(m=k^2\in[0,1]\). \(m=0\) is trigonometric; \(m=1\) is hyperbolic; otherwise AGM / descending Landen. Cap `JACOBI_AGM_MAX = 128`. Domain errors are `NaN`.
- Period reduction by \(4K(m)\) runs only when \(\lvert u\rvert\ge\pi\) and uses a single Carlson \(K\) pass (no nested Ziv). `RoundingMode::None` at 512 bits is golded.
- `expr!` leaves for every Jacobi method above.
- Elementwise wrappers on `ExactNumArray` and software IEEE arrays for the full Jacobi set (including quotients).
- Hex limb gold `jacobi_sn_1_half` at 64 bits.
- Golds: zeros at \(u=0\); \(m=0\) vs sin/cos; \(m=1\) vs tanh/sech; Pythagorean identities; `sn(K/2)`; `sn(u+4K)=sn(u)`; \(\partial_u sn=cn\,dn\); `ns·sn=1`; `F(sn(u|m)|m)=u`.

## 1.0.0 — 2026-08-31

First crates.io release.

- HDF5 array I/O (feature `hdf5`): crates.io `hdf5-rust` 1.0, no `libhdf5`. `Ieee64Array` / `Ieee32Array` / `ExactNumArray` `to_hdf5` / `from_hdf5`. Golds: 100×3 binary64 bits; 50×50 `ExactNum` at 256 bits; 1-D 1000 binary32 bits; nested `results/data`; append 10×3→20×3; wrong dataset name `Err`. CSV unchanged.
- Pre-publish gate: `scripts/zenith_prepublish.sh` (12 checks) and `scripts/ci_full.sh`. dashu 0.6.0 has no `euler_gamma` and no scoped rounding closure. `LIBRARY.md` `expr!` leaves include `gammainc_upper`, `ai`, `bi`.
- Hex limb CI: `golds/hex/reference.txt` (35 rows) bit-identical on i686 musl, wasm32-wasip1, and aarch64 musl vs x86_64 `to_bytes()`.
- Integer SIMD for software IEEE arrays: add/sub/mul/div/sqrt/fma (`IEEE_SIMD_LANE_WIDTH=4`). Gold: 1000-element `Ieee64Array` bit-identical to the scalar kernel and to `ExactNum` rounded to binary64.
- `lazy_static` uses `spin_no_std` so `thumbv7em-none-eabihf` builds with `default-features = false`.
- Locked `expr!`/`cexpr!` composite golds: `erf+erfc=1` at 256 bits; `J_0²+Y_0²` at working precision; complex `erf` at `p_wrk`.
- Criterion specials (64–1024 bits) and linalg (matmul/FFT/LU). `compare-bench.sh --quick` includes dashu. Baselines refreshed.
- `ExactComplex::mul` rounds each of the four real products at `(p, rm)` before the add/sub, matching the MPFR componentwise gold. `RoundingMode::None` series paths are unchanged.
- `proptest` (`PROPTEST_CASES=1000`): add commutes; directed round-then-coarser; `erf` odd; small-integer 2×2 matmul associativity; `ExactRational` `(a+b)-b=a`.
- `mpfr-tests`: GMP oracle for `ExactRational`; real-axis complex `erf`/`Γ`/`Ai`/`J_n` vs MPFR; `mpfr_gamma_inc`; identity golds where GNU MPFR/MPC have no function. GNU MPC 1.3 has no `erf`/`gamma`/Bessel and is not used.
- `HELP.md` rewritten as the user guide: precision model, rounding, `Consts`, macros vs methods, 30 recipes, mistakes, FAQ.
- `GETTING_STARTED.md` covers software IEEE, arrays, `(p, rm, cc)` specials, `ExactRational`/`ExactInt`, `cexpr!` cuts, and `Ball`.
- CSV for `Ieee64Array` / `ExactNumArray`: bit-pattern cells (empty → `NAN`); `Display@p=` for exact arrays. Golds: 100×3 bit round-trip; missing cell is `NAN`; extra column `Err`.
- Binary interchange: 16-byte big-endian inline record (`to_inline_bytes` / `write_inline_bytes`); heap `u32` limbs for wider mantissas; `ExactNumArray` packed records. Golds: finite / Inf / NaN flag / `p=256` heap / array shape / invalid `Err`; `u32::MAX+1=2^{32}` at `p=64`.
- Serde (`serde` feature): `ExactNum` / `ExactComplex` / `ExactRational` / `ExactInt` / arrays / `Ball`. Decimal strings carry `@p=`. IEEE arrays serialize integer bit patterns. Shape mismatch is `Err`.
- `IEEE_SIMD_LANE_WIDTH=4` public constant. `From<LayoutError>` → `MemoryAllocation`. LU/QR/SVD workspace `try_reserve_exact` returns `None` instead of aborting on reserve failure.
- SHA-256 / SHA-512 / HMAC-SHA-256 / `constant_time_eq`. Golds: FIPS empty and `abc`; RFC 4231 HMAC TC1.
- Modular `ExactInt`: `mod_pow`, `mod_inv`, `miller_rabin`, `pollard_rho` (Brent). Golds: `2^{100} ≡ 976371285 (mod 10^9+7)`; `mod_inv(3,7)=5`; Miller–Rabin on `2^{31}−1`; `8051=83×97`.
- Window functions: symmetric Hann / Hamming / Blackman, Kaiser (`I_0`), rectangular. Golds: `hann(4)=[0,3/4,3/4,0]`; Hamming endpoints `0.08`; Kaiser `β=0` is rectangular.
- Discrete transforms: `dct`/`idct` (type II / III), `dst`/`idst`, `fft_real`/`ifft_real`. Golds: `idct(dct(x))=x`; constant DCT is DC only; cosine energy at bins `k` and `N-k`; Parseval. `DSP_MAX_POINTS=2048`.
- ODE solvers: `rk4`, `rk45_adaptive` (Dormand–Prince 5(4)), `euler`. Golds: `y'=-y` RK4 1000-step error `<10^{-12}`; RK45 meets `atol=10^{-12}`; Euler 1000 vs 2000 is `O(h)`.
- Root finding: `bisect`, `newton`, `brent`, `illinois`. Golds: `bisect(sin,[3,4])=π`; `newton(x²−2)=√2`; Brent fewer iterations than bisection; `bisect(sin,[0,1])=None`.
- Quadrature: `gauss_legendre`, `tanh_sinh`, `gauss_laguerre`, `gauss_hermite`. Golds: \(x^2\) on \([-1,1]\) is \(2/3\); 20-point \(x^{38}\) is \(2/39\); \(1/\sqrt{1-x^2}=\pi\); Laguerre \(x^2=\Gamma(3)=2\); Hermite \(1=\sqrt{\pi}\).
- Orthogonal polynomials: `hermite_he`/`hermite_h`, `laguerre`/`gen_laguerre`, `chebyshev_t`/`chebyshev_u`, `gegenbauer`. Golds: `He_4(0)=3`, `L_3(0)=1`, `T_5(cos(π/5))=-1`, `C_2^{(1)}=4x²−1=U_2`, `2C_2^{(1/2)}=3x²−1`, `T_6=2xT_5−T_4`. `ORTHOPOLY_N_MAX=256`.
- Chebyshev interpolation: `chebyshev_coeffs` / `chebyshev_eval` / `clenshaw` / `chebyshev_error_bound`. `exp` on `[-1,1]` with 20 terms error `< 10^{-15}`; Clenshaw of `[1,2,3]` at `1/2` is `1/2`; `CHEBYSHEV_MAX_DEGREE=256`.
- `ExactNumPoly` — dense univariate (low-to-high coeffs); `div_rem(x²−1, x−1)=(x+1, 0)`; `gcd=x−1`; `compose(x², x+1)=x²+2x+1`; `∂(x³)=3x²`; `roots_real(x²−2)=±√2` (companion closed form through `POLY_COMPANION_CLOSED_DEG=2`).
- RNG: `random_uniform`, `random_gaussian` (Box–Muller), `random_exponential`, `ExactNumArray::random_fill` / `RandomDist`. Existing `random_normal(p, exp_from, exp_to)` is unchanged.
- Distribution kernels: `normal_pdf`/`cdf`, `gamma_pdf`, `beta_pdf`, `poisson_pmf`, `binomial_pmf`, `chi_squared_cdf`, `student_t_pdf`. Golds: \(1/\sqrt{2\pi}\), \(1/2\), \(e^{-1}\), \(\chi^2_2(2\ln 20)=19/20\).
- `ziv_round_vec` — same Ziv loop on a slice; `hypot(3,4)=5`; `atan2(1,1)=π/4` at 256 bits.
- `parse_exact` / `format_exact` — `0.1` is exact `1/10`; dyadic `0.5` / `0.125`; `parse("1.5e3")=1500`.
- `ExactInt` — limb integer; `20!`, `gcd(48,18)=6`, `2^100`, `div_rem(17,5)=(3,2)`.
- `ExactRational` — exact `num/den` with integer GCD reduction; `1/3+1/6=1/2`; `2/4=1/2`; sign onto the numerator; 256-bit `1/3`.
- `doc/REPRODUCIBILITY.md` — what determines a result, how to replay tests, citation line.
- `# Precision` rustdoc on `ExactNum` and `ExactComplex` specials (algorithm, thresholds, Ziv / ULP, MPFR oracle).
- `ExactNumArray::fft` / `ifft` — radix-2 Cooley–Tukey; impulse / cosine / Parseval golds; `FFT_MAX_POINTS=4096`.
- `ExactNumArray::eigen_decomp` — symmetric QR; \(Av=\lambda v\), \(V\Lambda V^T=A\); \(\begin{pmatrix}2&1\\1&2\end{pmatrix}\to(3,1)\); `EIGEN_ITER_MAX=64`.
- `ExactNumArray::svd_decomp` — Golub–Reinsch; \(U\Sigma V^T=A\), \(U^\top U=V^\top V=I\); \(\operatorname{diag}(3,2)\to\sigma=(3,2)\); `SVD_ITER_MAX=64`.
- Living docs are the build plan, capabilities, and build checklist only. `ZENITH_FLOAT_TODO.md` is retired.
- `ExactNumArray::qr_decomp` — modified Gram–Schmidt; \(QR=A\), \(Q^\top Q=I\); rank-deficient zero diagonal.
- `ExactNumArray::lu_decomp` — partial pivoting; \(PA=LU\) gold; singular → `None`.
- `ComplexBall` disk enclosures: `add`/`mul`/`exp`/`ln`/`sin`/`cos`. Unit-disk `exp` and \(\sin^2+\cos^2=1\) golds.
- `Ball::{cos,ln,sqrt,erf,bessel_j0,bessel_j1}` certified enclosures (`BALL_TRANSCENDENTAL_ERROR_TERMS`). `ExactNumArray` `(2×3)` `sin` and `bessel_j_nu` golds; `signum` ufunc.
- `ExactComplex::hypergeom_2f1` — series / Euler / Pfaff / Kummer in \(\mathbb{C}\). `cexpr!` leaf `hypergeom_2f1`. Cut on \([1,+\infty)\); non-positive integer \(c\) → NaN.
- `ExactComplex` elliptic \(K,E,\Pi\) (complete and incomplete) via Carlson \(R_F,R_C,R_D,R_J\) in \(\mathbb{C}\). `cexpr!` leaves `elliptic_k` / `elliptic_e` / `elliptic_f` / `elliptic_e_inc` / `elliptic_pi` / `elliptic_pi_inc`. \(K(1)=+\infty\); cut of \(K\) on \([1,+\infty)\).
- `cexpr!` — per-part `errs[]` for real vs imaginary cancellation; principal branch cuts documented; leaves match the complex-capable `expr!` set (`cbrt`/`root`, logs/exps, `hypot`/`fma`, `abs`/`arg`/`conj`, `ldexp`/`scalb`/`logb`). No `atan2` or `rem_pi` (complex trig uses `x+iy` identities).
- `ExactComplex::{log2, log10, log, log1p, exp2, exp10, expm1, ldexp, scalb, logb, cbrt, nth_root, hypot, fma, mul_add}`.
- Getting-started guide (`doc/GETTING_STARTED.md`) and help tutorial (`doc/HELP.md`).
- Maintainer build checklist is not packaged on crates.io.
- `Context::with_rounding_mode`.
- `ExactNum::{two_sum, two_product, fused_sum, fused_dot, polyval}`.
- `ExactComplex`: `tan`, `sinh`, `cosh`, `tanh`, `sqrt`, `pow`, `asin`, `acos`, `atan`, `asinh`, `acosh`, `atanh`.
- `ExactComplex` rectangular complex type (`+ − × ÷`, `exp`, `ln`, `sin`, `cos`, `abs`, `arg`).
- `ExactNum::sin_cos`; `mul_add` alias of `fma`.
- Special functions: `erf` / `erfc`, `gamma` / `ln_gamma`, integer-order `bessel_j`; `expr!` leaves.
- MPFR compare coverage for `fma`, paired `sin_cos` / `sinh_cosh`.
- MPFR compare coverage for `exp2`, `exp10`, and `rem_pi`; release gate on Linux x86_64 in `scripts/ci.sh`.
- MPFR fuzz: `add` / `mul` / `sqrt` bit-exact under all five IEEE rounding modes (`tests/mpfr/fuzz_round_modes.rs`).
- Progressive `ConstCache` (`cache_info`, √2, φ) and `CachedFBig` extra-precision wrapper.
- Public `Ball` enclosures and `ziv_round` correct-rounding retry helper.
- Stack-inlined mantissas (`INLINE_WORDS` = 2 limbs) before heap promotion.
- `doc/EXPR.md` — `expr!` per-op rounding contract (working precision, cancellation slots, final `set_precision`).
- `fma` skips the full-width product when `|a*b|` and `c` differ by more than `p + 2` words.
- MPFR fuzz covers `nth_root` for n=4,5 at 1-ULP (`mpfr_rootn_ui`).
- Optional nightly bench gate: `CI_BENCH=1 ./scripts/ci.sh`.
- `Consts::euler_gamma` (Euler–Mascheroni γ) on the progressive constant cache.
- `ExactNum::{frexp, ldexp, scalb, logb, ilogb}` — IEEE-style exponent split without hardware floats.
- `LowerExp` / `UpperExp` (`{:e}` / `{:E}`) and `LowerHex` formatting.
- `SharedConsts` — mutex-wrapped constant cache for parallel batch evaluation (`std`).
- MPFR oracles for `erf`/`erfc`, `gamma`/`ln_gamma`, `bessel_j`, extra constants, and `ExactComplex` add.
- `expr!` constants `sqrt2`, `phi`, `euler_gamma`.
- `expr!` leaves `ldexp` / `scalb` / `logb`.
- `serde`: JSON round-trip for `ExactNum` (decimal string / integer, Inf/NaN) and `ExactComplex`; feature implies `std`.
- Structured radix 2–36 parse/format round-trip stress test (`tests/radix_roundtrip.rs`).
- Native `+`, `-`, `*`, `/` operators on `ExactNum` (default precision, round-to-even).
- `exact!` / `fbig!` compile-time decimal float literals.
- Parse/format for radices 2–36 (`Radix::try_new`); `_e` exponent separator for bases > 10.
- `RadixFloat` wrapper and `ExactNum::with_radix` for radix-tagged I/O.
- `ExactNum::nth_root` / `expr!` `root(x, n)` — general n-th root (`sqrt`/`cbrt` delegation, composite factors, Newton for primes).
- `ExactNum::sinh_cosh` — paired hyperbolic evaluation with a single `exp(|x|)` path.
- `ExactNum::copysign`, `ExactNum::next_after` — software IEEE sign and successor.
- MPFR bit-oracle for `copysign` at high precision (`compare_copysign_test`); sign is applied before rounding so directed modes match MPFR.
- Seeded test RNG (`reseed_random`, `ZENITH_TEST_SEED`; panic hook prints the seed).
- `scripts/ci.sh` enforces debug / MPFR wall-time budgets (10 min / 30 min).
- Allocation failure at huge precision returns `NaN` (`Error::MemoryAllocation`), with a large-precision add/mul smoke test.
- `ExactNum::fma` / `expr!` `fma(a, b, c)` — fused multiply-add with single final rounding.
- Criterion benchmarks in `zenith-float-num/benches/` (arithmetic, transcendentals, composite); `./scripts/bench.sh`.
- FFT-scale multiply benchmark tier (`arithmetic/mul_fft`, 346k and 524k bits).
- TSV benchmark baselines (`doc/bench-baselines.tsv`) and `scripts/bench-compare.sh` for regression checks without JSON.
- Cross-library compare harness (`zenith-float-compare/`, `scripts/compare-bench.sh`) for release vs astro/dashu.
- Newton reciprocal at three words and up.
- `hypot`, `atan2`, `log1p`, and `expm1` on `ExactNum` and in `expr!`.
- `exp2`, `exp10`, and public `rem_pi` on `ExactNum` and in `expr!`.

## 0.1.0

First public release.

- Software-limb `ExactNum` (no hardware floating-point in calculations).
- Elementary functions, constants cache, radix parse/format, `expr!`.
- Optional `std`, `random`, and `serde` features.
- Dual license MIT OR Apache-2.0.
