# zenith-float 2.0.0 release checklist

**Release verification: 0/66 applicable items verified; 1 N/A (1.9), 67 items total. Not ready for release.**

Implementation inventory baseline: `origin/main` at `8c7f9729a47bd1af86381f408f648aa7c6cf5d7a`, inspected for this docs update. No whole-library release gate on one exact SHA is established here. Zero verified items does not mean 66 failed tests or 66 missing implementations.

The bound layer is integral to zenith-float 2.0: one library, one version, one release checklist, and one release gate. The legacy standalone bound-layer repository (history now imported under `zenith-float-bound/`) is a historical code source only, not a source of current status, versions, or a separate release/dependency pair. Integration work in P3/P4 is still part of this library's release.

Status keys:
- **MERGED / UNVERIFIED**: implementation is present at the pinned baseline, but whole-library release evidence is missing.
- **BRANCH-ONLY / UNVERIFIED**: work is identified on an unmerged branch, not established as merged.
- **NOT STARTED / UNVERIFIED**: the named integration or task is not present at the baseline; this is not a claim that historical source code is absent.
- **UNKNOWN / UNVERIFIED**: current implementation or execution status has not been established. An old task note alone does not prove work never started.
- **REQUIRED / UNVERIFIED**: a binding acceptance criterion, not an assertion that the implementation exists, passes, or is certified.
- **N/A**: excluded from the applicable denominator, with a reason. Only 1.9 is N/A.

Every applicable item must pass its criterion in the whole-library release gate on **one exact candidate SHA**, with evidence recorded against that SHA, and Jeff must sign off before 2.0 is ready. A merged fix, branch test, earlier gate, or rules approval is not release verification.

Inventory (count each top-level item once; the historical backport subtable is not extra items): P1 10 + P2 7 + P3 7 + P4 8 + P5 7 + P6 15 + full-MPFR requirement 1 + headroom 6 + FP16/FFT/AIR rules 6 = **67 total = 66 applicable + 1 N/A**.

Historical context: the original checklist was written 2026-10-23 around 09:30 CDT against dirty `wip/option-a` at `40c4b5e`. Its scratch paths, logs, branch tips and test results below are historical, non-release-SHA evidence, not live jobs or current completion claims. Those artifacts have not been revalidated in this docs-only update. Later FP16/FFT entries date from 2026-11-06. The 1.0.13 release was cancelled; versions 1.0 through 1.0.12 were yanked and must never be reused.

Required claim wording (a verification requirement, not a current release claim): "Every result from zenith-float's bound layer is a guaranteed enclosure: the true value is always inside; when a tight one can't be proved, the result says so (is_informative). Plain zf operations are correctly rounded in each rounding mode, claimed per function only once verified against MPFR. Never \"exact\"."

## P1 - zf Stage A (owner: Rusty)

| # | Item | Status | Pass criterion | Evidence |
|---|---|---|---|---|
| 1.1 | Stage A gate on fe5c15d | UNKNOWN / UNVERIFIED: historical r2 had all steps rc=0 except MPFR timeout rc=124 (372 passed, 0 failed); r3 was queued. Neither is candidate evidence | all steps rc=0 incl. mpfr | /tmp/zf-gate-logs/fe5c15d-r2.rc, fe5c15d-r2-mpfr.log, fe5c15d-r3-launcher.log |
| 1.2 | Full gate on 689431a (gatefix 75fa2d5 merge) | UNKNOWN / UNVERIFIED: historical gate was queued; final result not established here | all steps rc=0 | /tmp/zf-gate-logs/689431a-launcher.log |
| 1.3 | ln/log2/log10 clone conversion (rework, no 7th shape-A site) | MERGED / UNVERIFIED: ln/log2 and log10 conversions are on main (89adbb8, 3f869a1); historical rework results do not verify the release | 1.1 green, committed, clone 483->482 / combined 804->803 baselines, rules scan 0, digest 15016580156994057137 unchanged | zf-progress.md 07:35-09:18; /tmp/zf-logs/mpfr-rework.log; conv-1-ln-log2.patch, conv-2b-log10.patch |
| 1.4 | W6: pow +-Inf red tests + fixes; powi/powsi of -0 sign | MERGED / UNVERIFIED: pow_specials.rs is tracked and KNOWN_WRONG is empty at the baseline; signed-zero/Inf fixes are merged Scoped evidence (not a release verification): red tests a4d87e6 (17 KNOWN_WRONG); fixes 792cc4f (9 removed), 432852b (2), 71848b8 (4), ee11f4a (KNOWN_WRONG empty). NaN^0 / 1^NaN stay NaN by the documented NaN-operand rule (pow_with_a_nan_operand_is_nan), not a KNOWN_WRONG entry. | KNOWN_WRONG empty, all cases == MPFR mpfr_pow/pow_ui/pow_si, CHANGELOG line per wrong answer | /tmp/zf-logs/w6-pow_specials.rs, w6-mod.diff; signature table "W6" row. Scoped run on Pop 2026-11-07, source equal to 8c7f972 (log has no SHA header): `cargo test -p zenith-float-num --features mpfr-tests --test integration_tests -- pow` 11 passed 0 failed incl. pow_inf_and_signed_zero_powers_match_mpfr, pow_with_a_nan_operand_is_nan (~/pop-logs/rusty/2026-11-07/p1-pow-specials-mpfr-pow-tests.log) |
| 1.5 | 1.0.13 backport fixes present on option-a | MERGED / UNVERIFIED: backport-related fixes and tests are present; the historical 13/14 inventory below is not a current completeness verdict | every row of the table below yes/yes/yes | /tmp/bp-check (scratch crate on git archive of 40c4b5e) |
| 1.5a | NEW: pow of an exact quarter power (16^(k/4), 81^(k/4), 256^(k/4), 625^(k/4), k odd, both signs) | MERGED / UNVERIFIED: pow_of_exact_and_small_cases and its large-precision test are present; historical failures at 40c4b5e are not current test results Scoped evidence (not a release verification): red sweep c1f7330 (port of backport 085bdff); fix 56cea19 (exact roots for N/2^j exponents, j<=6; ports the intent of backport 43db66d, not the patch, whose word-rounded pow(x,1) must not return); follow-up f78f4d6 (underflow through an exact root). | red test (port backport 085bdff pow_of_exact_and_small_cases sweep: 14 bases x int/half/quarter exp x p 1..4096 x 6 modes) then fix (backport 43db66d: integer power of an exact root), all == MPFR | /tmp/bp-check/crate/src/bin/powsweep.rs output. Same scoped Pop 2026-11-07 run (source equal to 8c7f972, log has no SHA header): pow_of_exact_quarter_powers, pow_of_exact_and_small_cases, pow_of_exact_and_small_cases_at_large_precision (zenith-float-num/tests/mpfr/wide_primitives.rs) passed vs MPFR (~/pop-logs/rusty/2026-11-07/p1-pow-specials-mpfr-pow-tests.log) |
| 1.6 | 50-run multithreaded flaky check | UNKNOWN / UNVERIFIED: historical run paused at 17; 16/16 completed runs passed (104 tests each), not 50/50 | 50/50 rc=0 | /tmp/zf-logs/flaky50/summary.txt |
| 1.7 | The 5 W1 wrong answers + the 6th | MERGED / UNVERIFIED: 3cd2b14, 19b5cf8, 5ae5113, 67377c9, 79b9d34 and 4bdef5a are baseline ancestors; gate confirmation still required | each test in wide_primitives.rs green vs MPFR in the gate's mpfr step | git log; CHANGELOG 2.0.0 sections "fma with a far-away addend", "reciprocal of a subnormal", "nth_root of an exact power", "hypot with a subnormal", "atan2 with an extreme quotient", "nth_root of a negative number". Historical red tests b651f27, 4bdef5a W1b and 32c8d54; release gate confirmation pending 1.1/1.2 |
| 1.8 | MPFR emin/emax-restore bug in run_compare_ops | MERGED / UNVERIFIED: MPFR range Drop guard and restore test are present (888e70a); release-SHA execution still required | run_compare_ops (zenith-float-num/tests/mpfr/compare_ops_test.rs:44) restores emin/emax on every exit path (guard), test proves following tests see the default range | - |
| 1.9 | Cancelled 1.0.13 launcher at 61e4c3b | N/A: 1.0.13 was cancelled; this obsolete launcher task is closed. Do not rerun it or reuse a yanked version | - | Historical /tmp/zf-gate-logs/61e4c3b-launcher.log: missing scripts/use_rust_1_99.sh |

### Historical 1.5 backport table (40c4b5e vs 61e4c3b)

Snapshot from 2026-10-23 only. Every yes/no/partial, failure and path in this subtable describes that historical comparison, not the pinned baseline or a release-SHA result. In particular, quarter-power tests and fixes have since landed; see 1.5a.

| backport fix | on option-a | tests | CHANGELOG | evidence |
|---|---|---|---|---|
| 64f0692 reciprocal tiny negative keeps sign | yes | yes | partial | scratch: recip(-2^MIN) p=129/256 ToEven = -Inf; test reciprocal_of_a_tiny_negative_operand_keeps_its_sign (wide_primitives.rs:726, 4bdef5a); CHANGELOG reciprocal section covers subnormal operand, no plain line for the sign |
| c429dcb atan2 +-pi/2 for either sign of x | yes | yes | yes | 79b9d34; scratch atan2(+-2^g, -3*2^-g) = +-1.5707963; CHANGELOG (1) |
| 37c9490 hypot far-smaller operand | yes | yes | partial | scratch hypot(1,2^-2147483643) Up = 1+ulp inexact, ToEven 1 inexact; hypot(1,2^-9887) not NaN; test hypot_of_one_and_a_tiny_number_is_not_nan (:690, k up to 2147483000); CHANGELOG names subnormal case only, not the normal-range NaN case |
| 75b6a27 fma far-smaller term | yes | yes | yes | 3cd2b14/fb5686a; scratch fma(1,1,2^-300,128) Up/FromZero = 1+2^-127, others 1, all inexact; tests fma_with_a_far_smaller_addend (:522) |
| 6af2c9d nth_root exact flag | yes | yes | yes | 5ae5113; scratch nth_root(2^7,7) p=64/128 ToEven/ToOdd = 2, inexact=false; test nth_root_exact_result_is_not_inexact (:363) |
| 9ebef34 nth_root composite n | yes | yes | yes | 7eeee98 (one rounding); scratch nth_root(607*2^57,4,64) Down/Up bracket the root, checked with mpmath; test nth_root_with_a_composite_n (:642); CHANGELOG line 128 |
| 4523e3b atan2 tiny quotient | yes | yes | yes | 79b9d34; scratch atan2(2^(MIN-70),1,128) nonzero not NaN; test :765 |
| 366d963 decimal-reachable gaps | yes | yes | n/a (tests) | 32c8d54 atan2_at_the_scales_users_reach (:805); scratch hypot(1,3*2^-9966) Up p=256 and atan2(3*2^-9966,1) Down not NaN |
| 085bdff pow red tests | partial | no (exact sweep not ported; only exact_results_are_decided_not_retried, powers of two) | n/a | no pow_of_exact_and_small_cases on option-a |
| 43db66d pow exact results | partial: powers of two (f492c2c) and halves ok; exact quarter powers wrong (1.5a) | no | partial (powers of two only) | powsweep: 192 failures, all bases 16/81/256/625 with exponent k/4 |
| 43db66d pow(x,1) precision rounded to a word | n/a by design: option-a rounds at exactly p bits (7eeee98); pow(-3,1,p=1) = -4 ToEven/Down/FromZero, -2 ToOdd/Up/ToZero (= MPFR) | covered by exact-p-bit tests | yes (exact-p-bit section) | backport-only behaviour change; do not port |
| 4a7f472 hex_golds requires std | yes | gate build-nostd/thumbv7em | yes | 021fa87; Cargo.toml:49-52 required-features=["std"]; CHANGELOG line 116 |
| 4d0597d no unwrap in fft.rs/parser.rs | yes | lint deny unwrap_used (585fb6e) | yes | b4e4811 (fft returns Err), 3b6d405 (digit parser); git grep unwrap() in fft/parser at HEAD: 0; CHANGELOG line 108 |
| 5935c2c / 1a9bd76 CHANGELOG 1.0.13 | n/a (1.0.13 cancelled) | - | 2.0.0 section must carry the lines (P6.4) | - |
| 61e4c3b version bump 1.0.13 | n/a (cancelled) | - | - | - |

## P2 - zf API freeze (owner: Rusty)

| # | Item | Status | Pass criterion | Evidence |
|---|---|---|---|---|
| 2.1 | Keep the bound layer's 22 functions and 4 constants exactly as in /tmp/zf-logs/zf-l1-signature-table.md | UNKNOWN / UNVERIFIED: historical signature table was provisional; its final signature and Jeff approval are not established here | Jeff signs the table; a compile probe (zf-l1-probe) against the release SHA passes unchanged | zf-l1-signature-table.md:3, :79-80; zf-l1-probe.fe5c15d.out |
| 2.2 | pub consts for exponent and precision limits | MERGED / UNVERIFIED (exported constants only): MIN_NORMAL_EXPONENT, MIN_PRECISION, MAX_PRECISION, DEFAULT_P and MAX_DEC_SCALE are exported. Cross-target semantics and supported-cap acceptance remain REQUIRED / UNVERIFIED Scoped detail (not a verification, no new test run): defs.rs at 8a79c40 defines Exponent = i64 on every target (EXPONENT_MAX = 2^62 - 1, EXPONENT_MIN = MIN_NORMAL_EXPONENT = -EXPONENT_MAX), pub type Precision = usize, MIN_PRECISION = 1 (p = 0 rejected with InvalidArgument), and a per-target MAX_PRECISION = min(2^60, usize::MAX / 4) - 1 (2^60 - 1 on 64-bit, 2^30 - 2 = 1_073_741_822 on 32-bit); DEFAULT_P = 128; MAX_DEC_SCALE = 1_000_000 and MAX_PREC_RETRY = 256 in common/util.rs; PROPTEST_CASES and TEST_EXP_BOUND doc(hidden). Landed 0de3919 (red) / 0e8a962 and a998ce6 / 02130c8 (2026-10-23). Existing tests present (not re-run here): tests/precision_limits.rs, tests/r14_consts.rs. No 32-bit target executes these tests today. Use by the bound layer's CONTRACT 9.4 is UNVERIFIED (no reference found). | consts added, documented, tested, used by the bound-layer CONTRACT 9.4. REQUIRED / UNVERIFIED: use fixed-width exponent semantics consistently across supported targets, never interpret exponents as usize or derive their meaning from pointer width. Document the supported precision cap for each target from the exact candidate. Requests beyond that cap must return a documented deterministic refusal before allocation or lossy conversion, never wrap, truncate, or silently reduce precision. Test at, below, and above each limit, including internal guard-width and allocation boundaries | git grep at 40c4b5e is historical evidence only; candidate target-specific evidence remains required |
| 2.3 | Wildcard-match fixes (clippy match_wildcard_for_single_variants, enum wildcards) | UNKNOWN / UNVERIFIED: historical wildcard-lint work at c471eca is not a whole-library result | lint denied throughout the unified library, clippy 0 | /tmp/zf-logs/s7b-wildcard-enum.log |
| 2.4 | Remaining Option A loop conversions; Stage B additive APIs | MERGED / UNVERIFIED (partial scope established): P1-P3 APIs and Consts::get_exact_prec exist on main; remaining loop conversions and P8/P13/P14 audits still need evidence | every planned change listed in 2.5 lands before freeze | signature table :59-63, :88 |
| 2.5 | Every API change still planned | UNKNOWN / UNVERIFIED: retain the planned review scope after this table; presence of some APIs is not freeze approval | Jeff decides each before 2.6 | signature table "Could force a later signature change" |
| 2.6 | cargo semver-checks vs 1.0.12 confirms the freeze | UNKNOWN / UNVERIFIED: candidate semver-check result not established | every reported break is listed in CHANGELOG migration; none unintended | - |
| 2.7 | API review for planned breaks (so no 2.0.x API change) | UNKNOWN / UNVERIFIED: final written API approval not established | written review signed by Jeff | - |

Historical API review scope retained for 2.5 (not a claim these implementations are missing): (a) constants decision 5 (signature change possible on Consts::pi/e/ln_2/ln_10/new and the cc parameter of ~20 fns); (b) P11 get_exact_prec (additive); (c) MIN_NORMAL_EXPONENT, MIN/MAX_PRECISION (additive); (d) W6 pow special values (behaviour); (e) bare NaN vs InvalidArgument for Inf args (behaviour); (f) overflow-follows-rm uniformity (behaviour); (g) P14: DECIDED no fallback (true directed rounding, no ToEven+1ulp); (h) exponent width (the original i32/i64 question; baseline now uses i64); (i) zenith_float::bound re-export + unified error/rounding type within this library (P4).

For 2.2, retain the constants/documentation review scope: exponent and precision limits; DEFAULT_P and MAX_DEC_SCALE; MAX_PREC_RETRY as a retry budget, documented as max(p+16384, 3p); hidden test-support constants; and exact-p Consts values in all six rounding modes. Presence of exports alone does not complete that review.

## P3 - bound layers L0a..L4 of zenith-float 2.0 on the frozen API (build: Spence; harness: Harper)

Historical code/design sources: `/tmp/zb-stageb/docs/{PLAN,CONTRACT,TEST_PLAN,REPLAN-2026-10-16,SPLIT-STAGE-B-PROPOSAL}.md` and the legacy standalone bound-layer repository. These locate source material only. Neither a legacy main tip nor its gate is a status source for this release. The baseline workspace has no integrated bound member; P3 implementations must be integrated and gated within zenith-float 2.0, not released as another library.

| # | Layer | Status | Pass criterion | Evidence |
|---|---|---|---|---|
| 3.1 | L0a Mag (integer-only) | UNKNOWN / UNVERIFIED: historical source branch spence/zb-l0a-mag-green* had partial layer evidence; unified integration and gate not established | L0 gate PLAN:302 | REPLAN 3a |
| 3.2 | L0b Ball (construction, add/sub/neg/abs, mul/sqr/powi, div/inv, sqrt/hypot, sum/dot, compare, set ops, set_prec, ...) | UNKNOWN / UNVERIFIED: historical source branch harper/zb-l0b at 69457f4; unified integration and gate not established | L0b-1..L0b-20 green, containment goldens vs mpmath | REPLAN 3b |
| 3.3 | L1 transcendentals | UNKNOWN / UNVERIFIED: historical plan said not started; current source implementation not assessed | PLAN:340-390 gate | REPLAN 3c |
| 3.4 | L2a complex, L2b/c poly/series/matrices, L2d specials | UNKNOWN / UNVERIFIED: historical plan said not started; current source implementation not assessed | PLAN gates | REPLAN 3c |
| 3.5 | L3 quadrature | UNKNOWN / UNVERIFIED: historical plan said not started; current source implementation not assessed | PLAN gate | REPLAN 3c |
| 3.6 | L4 adaptive precision, root isolation, verified solve, eigen, ODE (Lohner) | UNKNOWN / UNVERIFIED: historical plan said not started; current source implementation not assessed | PLAN gate | REPLAN 3c; full scope, no cuts (Jeff 2026-10-16) |
| 3.7 | bound-layer harness: S7 fuzz, S10-S13 tooling | UNKNOWN / UNVERIFIED: historical source branches harper/zb-s7-fuzz at dc382f2 and harper/zb-s10-s13 at c471eca; no unified gate established | Harper's 16 steps green | git log |

## P4 - the fold (owner: Rusty)

At the pinned baseline, `Cargo.toml` lists num, macro and compare; bound integration is not present. This does not establish that the historical source implementations are absent. All eight fold items are release-unverified.

| # | Item | Status | Pass criterion |
|---|---|---|---|
| 4.1 | workspace crate zenith-float-bound depending only on zenith-float-num; re-exported as zenith_float::bound behind default-on `bound` feature | NOT STARTED / UNVERIFIED: unified integration absent at baseline | cargo tree shows only num; facade builds with and without `bound` |
| 4.2 | one type for BoundError / NaN-error and for the rounding-mode duplicates | NOT STARTED / UNVERIFIED: unified integration absent at baseline | one definition, the other a re-export; semver-checks shows it |
| 4.3 | scan: num never imports bound | NOT STARTED / UNVERIFIED: unified integration absent at baseline | rules scan step fails on a planted `use zenith_float_bound` in num (selftest) |
| 4.4 | per-path zero baselines (bound: clone 0, hidden_copy 0, ref-sum 0, sanctioned_rebuild 2) + selftest proving bound can't use zf's allowance (also Harper's per-package baselines ask; bound integration is absent at the pinned baseline) | NOT STARTED / UNVERIFIED: unified integration absent at baseline | selftest red on a planted clone in bound, green otherwise |
| 4.5 | CODEOWNERS paths for bound | NOT STARTED / UNVERIFIED: unified integration absent at baseline | file lists paths, scan's exception check honours it |
| 4.6 | fold historical bound code into this library with --allow-unrelated-histories, keeping its SHAs | NOT STARTED / UNVERIFIED: unified integration absent at baseline | `git log` shows the imported source's original SHAs (including the selected historical harper/spence tips) reachable |
| 4.7 | lints per crate (forbid unsafe, deny unwrap/expect/wildcard) | NOT STARTED / UNVERIFIED: unified integration absent at baseline | clippy 0 per crate |
| 4.8 | no_std builds with and without `bound` on the 64-bit host (thumbv7em and i686 parts VOID: P6, Jeff 2026-11-10 19:15 CT "Yes 64 bit only please") | 64-bit host part: `build-nostd` gate step (bound needs std, so `--no-default-features` builds the facade without it); 32-bit parts VOID | host no_std build green |

## P5 - whole-workspace gate with a bound tier (owner: Harper, Rusty for zf steps)

| # | Item | Status | Pass criterion |
|---|---|---|---|
| 5.1 | Harper's 16 steps over the whole workspace | NOT STARTED / UNVERIFIED: bound tier is not integrated in the baseline workspace; historical layer-only steps are insufficient | all 0 |
| 5.2 | cargo-hack --each-feature (and feature powerset of num/facade) + no_std builds in full_gate.sh | UNKNOWN / UNVERIFIED: no candidate whole-library feature matrix established; host/no_std tooling exists | step rc 0 |
| 5.3 | fuzzing parser and arithmetic plus the bound-layer suite | UNKNOWN / UNVERIFIED: historical bound fuzz work does not establish whole-library fuzz coverage | N CPU-hours without crash/mismatch, corpus committed |
| 5.4 | scoped cargo-mutants on the rounding/arithmetic core (cloud agent) plus bound-layer code | UNKNOWN / UNVERIFIED: candidate mutation results not established | every surviving mutant triaged (killed by a new test or justified) |
| 5.5 | per-package pubfn coverage | UNKNOWN / UNVERIFIED: no whole-library per-package coverage established | 100% of pub fns tested per package |
| 5.6 | MSRV + rust-toolchain file (rust-version 1.99; CI accepts 1.99.x, 194c3c5) | MERGED / UNVERIFIED: rust-version and rust-toolchain.toml are present; exact-MSRV candidate run still required | build on exact MSRV; toolchain file consistent across crates |
| 5.7 | ~~Retained 32-bit targets: i686 and armv7 Linux; existing thumbv7em no_std checks~~ | VOID: P6, Jeff 2026-11-10 19:15 CT "Yes 64 bit only please". The gate no longer builds, lints or tests a 32-bit target (tools/test_full_gate_64bit_only.sh) | - |

## P6 - release checklist over the whole workspace

| # | Item | Owner | Status | Pass criterion |
|---|---|---|---|---|
| 6.1 | Full wide-range sweep vs MPFR: every public fn x 6 modes x precisions (1..23, 24..255, >=20000), subnormal, huge exponents, specials; batches W2 rounding, W3 parse/format, W4 elementary, W5 specials (mpmath) | Rusty | REQUIRED / UNVERIFIED: MPFR-grade numerical validation is the quality target; no candidate result is asserted | zero mismatches; gaps list in zf-wide-range-gaps.md all closed. Require reproducible exact-SHA comparisons across supported targets, rounding modes, exponent extremes, and per-target precision limits. Check values, signed zeros, rounding flags, overflow, underflow, and documented refusal boundaries. Use independently justified reference methods where direct MPFR cannot apply, with oracle and input provenance recorded |
| 6.2 | cargo semver-checks (whole workspace) | Rusty | UNKNOWN / UNVERIFIED: candidate result or approval not established | as 2.6 |
| 6.3 | cargo deny: licenses, bans (no C/Python/LAPACK bindings in normal deps; rug/gmp-mpfr-sys are dev-deps only), advisories | Rusty | UNKNOWN / UNVERIFIED: candidate result or approval not established | deny check 0 |
| 6.4 | Docs + CHANGELOG pass: plain line per wrong answer (add: reciprocal tiny-negative sign, hypot normal-range far operand NaN, pow quarter powers, W6), forensic-audit note (present at CHANGELOG:5, review with Jeff), CONTRACT in docs, puremp MIT notice (zb NOTICE-puremp), claim wording above | Rusty + Harper (bound) | MERGED / UNVERIFIED (partial docs): this update corrects release framing only; full documentation review and Jeff approval remain | Jeff reads and approves |
| 6.5 | Migration guide 1.0 -> 2.0 for Accumath, incl. bound | Rusty | UNKNOWN / UNVERIFIED: candidate result or approval not established | Accumath port follows it with no extra questions |
| 6.6 | docs.rs build: cargo doc --all-features, no warnings; README/doc examples as doctests | Rusty | UNKNOWN / UNVERIFIED: candidate result or approval not established | 0 warnings; every README block compiled |
| 6.7 | Zero clone/unwrap/panic/allow/unsafe/IEEE scans at release SHA | Rusty | MERGED / UNVERIFIED (tooling): historical scans at fe5c15d-r2 are not candidate evidence | all scans 0 / at baseline at release SHA |
| 6.8 | cargo publish --dry-run for the unified release's published workspace packages in dependency order (confirm whether compare is published); no separate bound publication | Rusty | UNKNOWN / UNVERIFIED: historical prepublish rc=0 at fe5c15d-r2 is not candidate evidence | each dry-run 0 |
| 6.9 | latex-rust external consumer pinned to the unified candidate and gated (Accumath is covered by 6.10) | Spence / Rusty | UNKNOWN / UNVERIFIED: candidate result or approval not established | consumer gate green |
| 6.10 | Accumath build + full gate against 2.0 candidate (vendored copy update) | Cassie | UNKNOWN / UNVERIFIED: candidate result or approval not established | Accumath gate green |
| 6.11 | Accumath answer-drift audit (incl. its ball uses) | Cassie, Rusty reviewing | UNKNOWN / UNVERIFIED: candidate result or approval not established | every output diff vs 1.0.12 explained and traced to a zf fix; every gold tuned to an old wrong zf value corrected against an outside oracle |
| 6.12 | Fresh re-audit on the RC with the same probe harness | Rusty (outside agent) | REQUIRED / UNVERIFIED: certification-grade evidence is the quality target; no completed audit or certification is asserted | no new finding. Require reproducible commands, the exact whole-library candidate SHA, target and toolchain versions, runner details, named oracles, input provenance, test counts, and retained result logs. Evidence must cover rounding and refusal boundaries, including per-target supported precision caps. This quality target asserts neither external certification nor market superiority |
| 6.13 | Version and dates (workspace 2.0.0 already; CHANGELOG date, CAPABILITIES history row) | Rusty | MERGED / UNVERIFIED (version only): workspace is 2.0.0; release dates and history still need review | match release day |
| 6.14 | Final read-only re-audit after fixes | outside agent | UNKNOWN / UNVERIFIED: candidate result or approval not established | no finding |
| 6.15 | Jeff's go + pre-publish check shown to him | Jeff | UNKNOWN / UNVERIFIED: candidate result or approval not established | explicit go |
| 6.16 | LAST item, after everything else is 100% done: delete every block still held in HOLD.md, then delete HOLD.md itself (Jeff, 2026-11-10 15:37 HOLD rule addendum) | Rusty | NOT STARTED: HOLD.md lists H01-H08 | no `// HOLD` file, no `// HOLD <id>:` mod line and no HOLD.md in the tree |

### Historical 6.11 read-only evidence: Accumath zf call sites (accumath @ 830e781a, 2026-10-21)
At that historical Accumath SHA, all zf use goes through `SoftFloat` in crates/engine_core/src/soft.rs (const RM = ToEven at soft.rs:24); SoftFloat::inf()/neg_inf() exist (soft.rs:55-61), so Inf, and -0 from negation/underflow, can reach every wrapper; subnormals reach any op whose input came from an underflowing chain.
- reciprocal: soft.rs:224 (reciprocal_rm, caller-chosen rm -> directed modes reach the 19b5cf8 / sign fix), soft.rs:233-234 recip; callers calc/qanalog.rs:141-142, engine_linear/src/sparse.rs:422 (generic F). Inf/-0/subnormal: yes.
- pow: soft.rs:243 (ToEven); exp2 via pow soft.rs:261; callers lfunction.rs:113,268,321,328,486,515,657,800,842,867,877,970,971,1088,1102; soft.rs:1355,1370; precision.rs:58; specials/elliptic.rs:250,311,561; specials/mod.rs:179; specials/cylinder.rs:146 (3^(-2/3)); specials/hypergeom.rs:338,382; calc/qanalog.rs:538; tests/combinatorics_golds.rs:95; engine_linear/src/sparse.rs:31. Most are 2^-k thresholds (exact power-of-two path, fixed); W6 +-Inf / -0 base can reach via precision.rs:58 and lfunction.rs:321/328/800 (user values). Quarter-power exact bug (1.5a) reachable only for user exponents like 0.25/0.75.
- atan2: soft.rs:352-356; precision.rs:131 (complex arg im.atan2(re), -0 and Inf reachable).
- hypot: soft.rs:413-414.
- fma / mul_add, nth_root / root / powi / powsi: no Accumath call site found (cbrt wrapper soft.rs:228 only).

## Historical decisions and wiring notes (not release verification)
- Constants decision 5: historical decision recorded 2026-10-23 (include all constants, list per 2.2).
- Historical One-Ball decision: zf Ball, ComplexBall and ziv_round are removed at the fold (P4); no further work on them.
- P14 has no fallback: every Ziv kernel does true directed rounding; audit in /tmp/zf-logs/p14-remainder-audit.md.
- MPFR test wiring (for Harper): crate rug "~1.20.0" (lock 1.20.1, features float/rand/integer/rational) and gmp-mpfr-sys "~1.6.0" (lock 1.6.8), both dev-dependencies under [target.'cfg(all(target_arch="x86_64", target_os="linux"))'.dev-dependencies] in zenith-float-num/Cargo.toml:45-47; tests behind feature `mpfr-tests`; files zenith-float-num/tests/mpfr/*.rs (common.rs, compare_*_test.rs, exact_bit_*.rs, fuzz_round_modes.rs, mag_parts.rs, overflow_modes.rs, rem_exact.rs, retry_budget.rs, subnormal_*.rs, wide_primitives.rs, + pow_specials.rs, now tracked at the baseline).
- MPFR suite tiers (Paul 2026-10-23): the original plan used a fast tier in ordinary gates; the baseline gate now defaults to full. The FULL MPFR suite (MPFR_TIER=full) is MANDATORY for every RC and release gate (see the release gate requirement below).
- Pinned oracle: /tmp/oracle-lib carries mpmath 1.3.0 (not 1.2.1 as stated in the brief); confirm which is pinned.

## Release gate requirement

Status: UNKNOWN / UNVERIFIED. The baseline full_gate.sh defaults to MPFR_TIER=full; wiring is not a completed run.
- [ ] Full MPFR suite (MPFR_TIER=full, --test-threads=1) green on the RC/release SHA. Mandatory; the fast tier alone does not qualify.

## Headroom beyond cosmology (Jeff 2026-10-23 12:33 CT)
Status: all six headroom items UNKNOWN / UNVERIFIED. Goal: exercise the unified library at its actual limits, not only at cosmology scales. The historical d873c7c i32 limits are superseded. At the pinned baseline, zenith-float-num/src/defs.rs defines Exponent = i64, EXPONENT_MAX = 2^62 - 1, EXPONENT_MIN = MIN_NORMAL_EXPONENT = -EXPONENT_MAX, MIN_PRECISION = 1 and MAX_PRECISION = min(2^60, usize::MAX / 4) - 1, which is per target (2^60 - 1 on 64-bit, 2^30 - 2 = 1_073_741_822 on 32-bit) while the i64 exponent range is the same on every target. Take limits from the exact candidate when gating; branch proposals are not the release contract.
- [ ] (a) Exponent extremes: an MPFR sweep (mpfr emin/emax set to zf's limits, which fit inside MPFR's default range of ±(2^62-1)) for every public function in all supported rounding modes, with inputs at and next to: the overflow edge (largest finite at p, and exponent EXPONENT_MAX-1 / EXPONENT_MAX), the normal floor (min_positive_normal), the subnormal floor (min_positive), and results that land just across overflow or underflow (for example exp/pow/mul/div/sq/hypot near ±EXPONENT_MAX, and log/sqrt/recip of subnormals). Check the value, the Inf/0/-0 result, and the overflow/underflow/inexact flags. Record the commands and test file.
- [ ] (b) High precision: every public function against MPFR at p = 1024, 2048, 4096, 8192, then 16384 and 65536, and near the largest p the box can allocate (MAX_PRECISION itself cannot be allocated; record the largest p tried). Record runtime per function and per p in a TSV (zf ms, MPFR ms, ratio). Any ratio above about 3x MPFR, or growth worse than M(p)·log p, gets a named cause (series term count, missing binary splitting, Ziv retries, constants-cache rebuild, quadratic multiply) before RC.
- [ ] (b2) Retry budget at high p: use the retry_budget test to show that no function returns PrecisionRetryExhausted at p ≥ 8192 on hard-to-round inputs (the budget there is max(p+16384, 3p)).
- [ ] (c) Accumath cosmology against the 2.0 candidate: run the existing FLRW tests in /tmp/acc-or (engine_core/src/geometry.rs: flrw_scale :1547, flrw_flat :1556, static_redshift :743, tests flrw_friedmann_constraint :3427, flrw_dust_and_exponential :3549) with Accumath pinned to the 2.0 candidate SHA. NOTE: these are symbolic (metric) functions. I found no numeric comoving/luminosity distance-redshift function in Accumath (rg 'comoving|luminosity_distance' only hits geometry.rs/lib.rs). Either add a numeric D_C(z)/D_L(z) check through SoftFloat at z up to 10^4 and parameter ratios of 10^±120 (and 10^±400), or record that the numeric cosmology path does not exist yet.
- [ ] (d) Bound-layer wide-range: run the unified zenith-float library's bound-layer suite on the same exact 2.0 candidate SHA, at the exponent extremes and p ladder in (a)/(b). This is part of the whole-library gate, not a separate dependency or legacy-repository gate.
- [ ] (e) Before the API freeze, decide on the limits flagged in the limits audit (MAX_DEC_SCALE is test-only and should be renamed or hidden; the u32 term-count clamps; whether 1M-bit+ precision is a supported claim).

## P-FP16 - FP16/BF16 ExactNum bridges + FFT/DCT caps (Paul / Rusty, 2026-11-06)

Status: **MERGED / UNVERIFIED** for FP16.1-3 and FFT.1/DCT.1: bridges landed in PR #26 (`c432a6a`), caps in PR #27 (`8c7f972`). AIR.1 is UNKNOWN / UNVERIFIED on the Accumath side. Historical rules approval and those merges are not a whole-library release gate. Implementation owner: Rusty (historical Slice 1 then Slice 2 plan at `1f6fc856`). These changes and bound-layer integration belong to the same zenith-float 2.0 release.

Companion: SoftFloat UBO split `#419a–c` / `#437a–d` (Spence `spence/ubo-419-437-native-airgap`)
and historical `~/paul-work/land-seq-2026-11-06/` copies of the split proposal.

### Rules (normative for 2.0)

| # | Rule | Status | Pass criterion |
|---|---|---|---|
| FP16.1 | `Ieee16` / FP16 bit-pattern ↔ `ExactNum` | MERGED / UNVERIFIED: release gate required | Lossless `u16` ↔ `ExactNum` round-trip for every finite/NaN/Inf encoding; test-first |
| FP16.2 | `Bfloat16` / BF16 bit-pattern ↔ `ExactNum` | MERGED / UNVERIFIED: release gate required | Same for BF16 `u16` encodings; test-first |
| FP16.3 | No IEEE host float in conversion path | MERGED / UNVERIFIED: release gate required | No `f32`/`f64` in the bridge; SoftFloat / ExactNum only |
| FFT.1 | Default per-call FFT length cap | MERGED / UNVERIFIED: release gate required | `FFT_MAX_POINTS = 65_536` (was 4,096; larger rejected unless the per-call `fft_with_max_points` / `ifft_with_max_points` cap allows it) |
| DCT.1 | Default per-call DCT length cap | MERGED / UNVERIFIED: release gate required | `DSP_MAX_POINTS = 16_384` (was 2,048; per call: `dct`/`idct`/`dst`/`idst`_`with_max_points`) |
| AIR.1 | Native air-gap mandate for SoftFloat classical ops | UNKNOWN / UNVERIFIED: external integration not established | `#419`/`#437` split: no `ort`, no IEEE-ML, no AGPL YOLO in AccuMath SoftFloat path |

### Ownership (unchanged)

- Rusty Slice 1: FP16.1 + FP16.2 (+ FP16.3) on `zenith-float` off `1f6fc856`.
- Rusty Slice 2: FFT.1 + DCT.1.
- Spence: UBO `#419a–c` / `#437a–d` text + exact-SHA reviews.
- Paul: this checklist section + landing coordination.

### Evidence expected before RC

- Red/green tests naming every encoding class.
- CHANGELOG 2.0.0 bullets for FP16/BF16 bridges and FFT/DCT caps.
- Gate: fmt / clippy / workspace tests / MPFR suite as required by existing Stage A.

