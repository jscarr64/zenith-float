# zf 2.0.0 release checklist (zenith-float + zenith-bound folded in)

Written 2026-10-23 ~09:30 CT by Rusty's worker. Untracked, not committed. Tip: wip/option-a @ 40c4b5e (working tree dirty: ln/log2/log10 rework, W6 pow_specials.rs).
Decisions: no 1.0.13; everything ships as 2.0.0 (Jeff 09:14). zb folded into 2.0, 2.0 waits for it (Jeff 09:19). 2.0 must be fully verified so no 2.0.1/2.0.x API change is needed.
Status keys: DONE (SHA + evidence) / IN PROGRESS / NOT STARTED / BLOCKED.

Doc claim wording (use verbatim): "Every zenith-bound result is a guaranteed enclosure: the true value is always inside; when a tight one can't be proved, the result says so (is_informative). Plain zf operations are correctly rounded in each rounding mode, claimed per function only once verified against MPFR. Never \"exact\"."

## P1 - zf Stage A (owner: Rusty)

| # | Item | Status | Pass criterion | Evidence |
|---|---|---|---|---|
| 1.1 | Stage A gate on fe5c15d | IN PROGRESS. r2 finished 08:20: every step 0 (build debug/release/nostd/thumbv7em, features, oracle fast/deep, rules, fmt, clippy x86_64/i686, audit, no-intrinsics, hex-32bit, unsafe-scan, prepublish) except mpfr=124 (20 min STEP_TIMEOUT; 372 passed, 0 failed). r3 (STEP_TIMEOUT 60m) queued on /tmp/pop-gate.lock since ~08:25 (pid 1711080) | all steps rc=0 incl. mpfr | /tmp/zf-gate-logs/fe5c15d-r2.rc, fe5c15d-r2-mpfr.log, fe5c15d-r3-launcher.log |
| 1.2 | Full gate on 689431a (gatefix 75fa2d5 merge) | IN PROGRESS: queued (pid 1702120, ~54 min waiting at 09:20), ahead of it pr184 driver, gaterun, c471-r199, 2 accumath runs | all steps rc=0 | /tmp/zf-gate-logs/689431a-launcher.log |
| 1.3 | ln/log2/log10 clone conversion (rework, no 7th shape-A site) | IN PROGRESS: uncommitted in /tmp/zf-oa; MPFR suite on reworked tree RC=0 (08:55); wide-range ln/log2/log10 10 tests ok; digests identical; waiting on Stage A re-green, then Jeff review per NO_CLONE_RULES | 1.1 green, committed, clone 483->482 / combined 804->803 baselines, rules scan 0, digest 15016580156994057137 unchanged | zf-progress.md 07:35-09:18; /tmp/zf-logs/mpfr-rework.log; conv-1-ln-log2.patch, conv-2b-log10.patch |
| 1.4 | W6: pow +-Inf red tests + fixes; powi/powsi of -0 sign | IN PROGRESS: red test file written (zenith-float-num/tests/mpfr/pow_specials.rs, untracked; KNOWN_WRONG list: pow(+-0,-Inf), pow(-x,+-Inf), pow(x>1,-Inf)=-Inf, pow(+-Inf,neg odd)=+0, NaN^0/1^NaN), fixes not started | KNOWN_WRONG empty, all cases == MPFR mpfr_pow/pow_ui/pow_si, CHANGELOG line per wrong answer | /tmp/zf-logs/w6-pow_specials.rs, w6-mod.diff; signature table "W6" row |
| 1.5 | 1.0.13 backport fixes present on option-a | IN PROGRESS: 13 of 14 present; **pow exact results only partial** (see 1.5a) | every row of the table below yes/yes/yes | /tmp/bp-check (scratch crate on git archive of 40c4b5e) |
| 1.5a | NEW: pow of an exact quarter power (16^(k/4), 81^(k/4), 256^(k/4), 625^(k/4), k odd, both signs) | NOT STARTED. Wrong on 40c4b5e: ToEven/ToOdd flagged inexact, Up/Down/ToZero/FromZero return NaN(PrecisionRetryExhausted) (192/9792 sweep cases, p=64 and 256). Halves (4^0.5, 9^1.5) and powers of two are right | red test (port backport 085bdff pow_of_exact_and_small_cases sweep: 14 bases x int/half/quarter exp x p 1..4096 x 6 modes) then fix (backport 43db66d: integer power of an exact root), all == MPFR | /tmp/bp-check/crate/src/bin/powsweep.rs output |
| 1.6 | 50-run multithreaded flaky check | IN PROGRESS, paused (SIGSTOP, state T) at run 17; 16/16 pass (104 tests each) | 50/50 rc=0 | /tmp/zf-logs/flaky50/summary.txt |
| 1.7 | The 5 W1 wrong answers + the 6th | DONE: fma 3cd2b14; subnormal reciprocal 19b5cf8; nth_root exact flag 5ae5113; hypot subnormal 67377c9; atan2 extreme 79b9d34 (fix 5/5); 6th negative nth_root odd n>=5 Up/Down 4bdef5a. Red tests b651f27 (5 known_wrong) + 4bdef5a W1b + 32c8d54 (atan2 user scales) | each test in wide_primitives.rs green vs MPFR in the gate's mpfr step | git log; CHANGELOG 2.0.0 sections "fma with a far-away addend", "reciprocal of a subnormal", "nth_root of an exact power", "hypot with a subnormal", "atan2 with an extreme quotient", "nth_root of a negative number". Gate confirmation pending 1.1/1.2 |
| 1.8 | MPFR emin/emax-restore bug in run_compare_ops | NOT STARTED | run_compare_ops (zenith-float-num/tests/mpfr/compare_ops_test.rs:44) restores emin/emax on every exit path (guard), test proves following tests see the default range | - |
| 1.9 | NEW: 61e4c3b gate launcher failed | OPEN: /tmp/zf-gate-logs/61e4c3b-launcher.log "scripts/use_rust_1_99.sh: No such file" (09:14). Moot now that 1.0.13 is cancelled; do not re-run | - | - |

### 1.5 backport table (wip/option-a @ 40c4b5e vs backport/1.0.13-fixes @ 61e4c3b)

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
| 2.1 | Keep zb's 22 functions and 4 constants exactly as in /tmp/zf-logs/zf-l1-signature-table.md | IN PROGRESS: table provisional ("not frozen"); constants decision 5 CLOSED (Jeff 2026-10-23: include all constants, list per 2.2) | Jeff signs the table; a compile probe (zf-l1-probe) against the release SHA passes unchanged | zf-l1-signature-table.md:3, :79-80; zf-l1-probe.fe5c15d.out |
| 2.2 | pub consts for exponent and precision limits | PARTIAL: EXPONENT_MAX (defs.rs:45), EXPONENT_MIN (defs.rs:48), EXPONENT_BIT_SIZE (defs.rs:69) exist, re-exported lib.rs:181-182; MAX_PREC_RETRY (util.rs:78, lib.rs:189) is a retry budget, not a limit. MISSING (planned additive): MIN_NORMAL_EXPONENT (REPLAN 3b: "does not exist"), MIN_PRECISION / MAX_PRECISION Decision 5 CLOSED (Jeff 2026-10-23): include all: MIN_NORMAL_EXPONENT, MIN_PRECISION=1 (p=0 rejected), MAX_PRECISION (num.rs:54-55 limit), pub DEFAULT_P and MAX_DEC_SCALE, MAX_PREC_RETRY doc fixed to max(p+16384, 3p), PROPTEST_CASES/TEST_EXP_BOUND doc(hidden); exact-p Consts (pi, ln_2, ln_10, e, from_mag_parts) correctly rounded at p in all 6 modes | consts added, documented, tested, used by zb CONTRACT 9.4 | git grep at 40c4b5e |
| 2.3 | Wildcard-match fixes (clippy match_wildcard_for_single_variants, enum wildcards) | IN PROGRESS (s7b-wildcard-enum.log; zb side c471eca) | lint denied in zf and zb, clippy 0 | /tmp/zf-logs/s7b-wildcard-enum.log |
| 2.4 | Remaining Option A loop conversions; Stage B additive APIs | IN PROGRESS: layer 0a/1a partial; P1-P3 DONE (leading_bit_exponent, top_u64_mantissa_and_sticky, from_mag_parts/pow2 e876fe9); P8 (ulp/next_up/next_down) re-verify NOT STARTED; P11 get_exact_prec NOT STARTED (additive); P13 re-verify NOT STARTED; P14 scoped Ziv audit NOT STARTED | every planned change listed in 2.5 lands before freeze | signature table :59-63, :88 |
| 2.5 | Every API change still planned | list: (a) constants decision 5 (signature change possible on Consts::pi/e/ln_2/ln_10/new and the cc parameter of ~20 fns); (b) P11 get_exact_prec (additive); (c) MIN_NORMAL_EXPONENT, MIN/MAX_PRECISION (additive); (d) W6 pow special values (behaviour); (e) bare NaN vs InvalidArgument for Inf args (behaviour); (f) overflow-follows-rm uniformity (behaviour); (g) P14: DECIDED no fallback (true directed rounding, no ToEven+1ulp); (h) exponent width i32 vs i64 (none planned); (i) zenith_float::bound re-export + unified error/rounding type (P4) | Jeff decides each before 2.6 | signature table "Could force a later signature change" |
| 2.6 | cargo semver-checks vs 1.0.12 confirms the freeze | NOT STARTED | every reported break is listed in CHANGELOG migration; none unintended | - |
| 2.7 | API review for planned breaks (so no 2.0.x API change) | NOT STARTED | written review signed by Jeff | - |

## P3 - zb L0b..L4 on the frozen API (build: Spence; harness: Harper)

Source: /tmp/zb-stageb/docs/{PLAN,CONTRACT,TEST_PLAN,REPLAN-2026-10-16,SPLIT-STAGE-B-PROPOSAL}.md; repo /home/jsc64/MathWorks/zenith-bound (main @ 3c9f36f; read with git log only).

| # | Layer | Status | Pass criterion | Evidence |
|---|---|---|---|---|
| 3.1 | L0a Mag (integer-only) | IN PROGRESS: green on spence/zb-l0a-mag-green*; L0a-5 fuzz, L0a-7 mutants, L0a-8 exit review open | L0 gate PLAN:302 | REPLAN 3a |
| 3.2 | L0b Ball (construction, add/sub/neg/abs, mul/sqr/powi, div/inv, sqrt/hypot, sum/dot, compare, set ops, set_prec, ...) | IN PROGRESS: harper/zb-l0b @ 69457f4 (07:57 today, L0b-7 enclosure goldens); needs P2/P3 (done) and MIN_NORMAL_EXPONENT (2.2) | L0b-1..L0b-20 green, containment goldens vs mpmath | REPLAN 3b |
| 3.3 | L1 transcendentals | NOT STARTED; gated on P14 and P11 (Stage B) | PLAN:340-390 gate | REPLAN 3c |
| 3.4 | L2a complex, L2b/c poly/series/matrices, L2d specials | NOT STARTED (L2b/c need only L0b) | PLAN gates | REPLAN 3c |
| 3.5 | L3 quadrature | NOT STARTED (after L1/L2a) | PLAN gate | REPLAN 3c |
| 3.6 | L4 adaptive precision, root isolation, verified solve, eigen, ODE (Lohner) | NOT STARTED | PLAN gate | REPLAN 3c; full scope, no cuts (Jeff 2026-10-16) |
| 3.7 | zb harness: S7 fuzz, S10-S13 tooling | IN PROGRESS: harper/zb-s7-fuzz @ dc382f2 (08:56), harper/zb-s10-s13 @ c471eca; gate c471-r199 running | Harper's 16 steps green | git log |

## P4 - the fold (owner: Rusty)

All NOT STARTED.
| # | Item | Pass criterion |
|---|---|---|
| 4.1 | workspace crate zenith-float-bound depending only on zenith-float-num; re-exported as zenith_float::bound behind default-on `bound` feature | cargo tree shows only num; facade builds with and without `bound` |
| 4.2 | one type for BoundError / NaN-error and for the rounding-mode duplicates | one definition, the other a re-export; semver-checks shows it |
| 4.3 | scan: num never imports bound | rules scan step fails on a planted `use zenith_float_bound` in num (selftest) |
| 4.4 | per-path zero baselines (bound: clone 0, hidden_copy 0, ref-sum 0, sanctioned_rebuild 2) + selftest proving bound can't use zf's allowance (also Harper's per-package baselines ask; zb is not a workspace member today: members = num, macro, compare) | selftest red on a planted clone in bound, green otherwise |
| 4.5 | CODEOWNERS paths for bound | file lists paths, scan's exception check honours it |
| 4.6 | merge zb with --allow-unrelated-histories keeping zb SHAs | `git log` shows zb's original SHAs (3c9f36f..., harper/spence branch tips) reachable |
| 4.7 | lints per crate (forbid unsafe, deny unwrap/expect/wildcard) | clippy 0 per crate |
| 4.8 | no_std builds with and without `bound` (host + thumbv7em + i686) | all build |

## P5 - whole-workspace gate with a bound tier (owner: Harper, Rusty for zf steps)

| # | Item | Status | Pass criterion |
|---|---|---|---|
| 5.1 | Harper's 16 steps over the whole workspace | NOT STARTED (zb-only today) | all 0 |
| 5.2 | cargo-hack --each-feature (and feature powerset of num/facade) + no_std builds in full_gate.sh | NOT STARTED: full_gate.sh:33-34 has build-nostd and build-thumbv7em, no cargo-hack | step rc 0 |
| 5.3 | fuzzing parser and arithmetic (zf) + zb fuzz | NOT STARTED for zf; zb S7 in progress | N CPU-hours without crash/mismatch, corpus committed |
| 5.4 | scoped cargo-mutants on the rounding/arithmetic core (cloud agent) + zb | NOT STARTED | every surviving mutant triaged (killed by a new test or justified) |
| 5.5 | per-package pubfn coverage | NOT STARTED for bound (zf has pubfn-coverage step) | 100% of pub fns tested per package |
| 5.6 | MSRV + rust-toolchain file (rust-version 1.99; CI accepts 1.99.x, 194c3c5) | IN PROGRESS | build on exact MSRV; toolchain file consistent across crates |
| 5.7 | i686 + thumbv7em builds | DONE on fe5c15d-r2 (clippy-i686=0, build-thumbv7em=0, hex-32bit=0); re-run on release SHA | rc 0 at release SHA |

## P6 - release checklist over the whole workspace

| # | Item | Owner | Status | Pass criterion |
|---|---|---|---|---|
| 6.1 | Full wide-range sweep vs MPFR: every public fn x 6 modes x precisions (1..23, 24..255, >=20000), subnormal, huge exponents, specials; batches W2 rounding, W3 parse/format, W4 elementary, W5 specials (mpmath) | Rusty | NOT STARTED (W1 done, W6 in progress) | zero mismatches; gaps list in zf-wide-range-gaps.md all closed |
| 6.2 | cargo semver-checks (whole workspace) | Rusty | NOT STARTED | as 2.6 |
| 6.3 | cargo deny: licenses, bans (no C/Python/LAPACK bindings in normal deps; rug/gmp-mpfr-sys are dev-deps only), advisories | Rusty | NOT STARTED (gate has `audit`=0) | deny check 0 |
| 6.4 | Docs + CHANGELOG pass: plain line per wrong answer (add: reciprocal tiny-negative sign, hypot normal-range far operand NaN, pow quarter powers, W6), forensic-audit note (present at CHANGELOG:5, review with Jeff), CONTRACT in docs, puremp MIT notice (zb NOTICE-puremp), claim wording above | Rusty + Harper (bound) | IN PROGRESS | Jeff reads and approves |
| 6.5 | Migration guide 1.0 -> 2.0 for Accumath, incl. bound | Rusty | NOT STARTED | Accumath port follows it with no extra questions |
| 6.6 | docs.rs build: cargo doc --all-features, no warnings; README/doc examples as doctests | Rusty | NOT STARTED (doctests run in mpfr-rework.log, 9 ok, num only) | 0 warnings; every README block compiled |
| 6.7 | Zero clone/unwrap/panic/allow/unsafe/IEEE scans at release SHA | Rusty | NOT STARTED (scans exist and are 0 at fe5c15d-r2 except baselines) | all scans 0 / at baseline at release SHA |
| 6.8 | cargo publish --dry-run in order: num, bound, facade, macro (compare not published?) | Rusty | NOT STARTED (gate prepublish=0 on fe5c15d-r2) | each dry-run 0 |
| 6.9 | zenith-bound and latex-rust consumers pinned to candidate and green | Spence / Rusty | NOT STARTED | their gates green |
| 6.10 | Accumath build + full gate against 2.0 candidate (vendored copy update) | Cassie | NOT STARTED | Accumath gate green |
| 6.11 | Accumath answer-drift audit (incl. its ball uses) | Cassie, Rusty reviewing | NOT STARTED; needs the Accumath port to the 2.0 API | every output diff vs 1.0.12 explained and traced to a zf fix; every gold tuned to an old wrong zf value corrected against an outside oracle |
| 6.12 | Fresh re-audit on the RC with the same probe harness | Rusty (outside agent) | NOT STARTED | no new finding |
| 6.13 | Version and dates (workspace 2.0.0 already; CHANGELOG date, CAPABILITIES history row) | Rusty | NOT STARTED | match release day |
| 6.14 | Final read-only re-audit after fixes | outside agent | NOT STARTED | no finding |
| 6.15 | Jeff's go + pre-publish check shown to him | Jeff | NOT STARTED | explicit go |

### 6.11 read-only part: Accumath zf call sites (accumath @ 830e781a, 2026-10-21)
All zf use goes through `SoftFloat` in crates/engine_core/src/soft.rs (const RM = ToEven at soft.rs:24); SoftFloat::inf()/neg_inf() exist (soft.rs:55-61), so Inf, and -0 from negation/underflow, can reach every wrapper; subnormals reach any op whose input came from an underflowing chain.
- reciprocal: soft.rs:224 (reciprocal_rm, caller-chosen rm -> directed modes reach the 19b5cf8 / sign fix), soft.rs:233-234 recip; callers calc/qanalog.rs:141-142, engine_linear/src/sparse.rs:422 (generic F). Inf/-0/subnormal: yes.
- pow: soft.rs:243 (ToEven); exp2 via pow soft.rs:261; callers lfunction.rs:113,268,321,328,486,515,657,800,842,867,877,970,971,1088,1102; soft.rs:1355,1370; precision.rs:58; specials/elliptic.rs:250,311,561; specials/mod.rs:179; specials/cylinder.rs:146 (3^(-2/3)); specials/hypergeom.rs:338,382; calc/qanalog.rs:538; tests/combinatorics_golds.rs:95; engine_linear/src/sparse.rs:31. Most are 2^-k thresholds (exact power-of-two path, fixed); W6 +-Inf / -0 base can reach via precision.rs:58 and lfunction.rs:321/328/800 (user values). Quarter-power exact bug (1.5a) reachable only for user exponents like 0.25/0.75.
- atan2: soft.rs:352-356; precision.rs:131 (complex arg im.atan2(re), -0 and Inf reachable).
- hypot: soft.rs:413-414.
- fma / mul_add, nth_root / root / powi / powsi: no Accumath call site found (cbrt wrapper soft.rs:228 only).

## Not phase-bound, still open from the logs
- Constants decision 5: CLOSED 2026-10-23 (include all constants, list per 2.2).
- One-Ball plan APPROVED: zf Ball, ComplexBall and ziv_round are removed at the fold (P4); no further work on them.
- P14 has no fallback: every Ziv kernel does true directed rounding; audit in /tmp/zf-logs/p14-remainder-audit.md.
- MPFR test wiring (for Harper): crate rug "~1.20.0" (lock 1.20.1, features float/rand/integer/rational) and gmp-mpfr-sys "~1.6.0" (lock 1.6.8), both dev-dependencies under [target.'cfg(all(target_arch="x86_64", target_os="linux"))'.dev-dependencies] in zenith-float-num/Cargo.toml:45-47; tests behind feature `mpfr-tests`; files zenith-float-num/tests/mpfr/*.rs (common.rs, compare_*_test.rs, exact_bit_*.rs, fuzz_round_modes.rs, mag_parts.rs, overflow_modes.rs, rem_exact.rs, retry_budget.rs, subnormal_*.rs, wide_primitives.rs, + pow_specials.rs untracked).
- MPFR suite tiers (Paul 2026-10-23): a fast tier runs in every full gate; the FULL MPFR suite (MPFR_TIER=full) is MANDATORY for every RC and release gate (see 6.x below).
- Pinned oracle: /tmp/oracle-lib carries mpmath 1.3.0 (not 1.2.1 as stated in the brief); confirm which is pinned.

## Release gate requirement
- [ ] Full MPFR suite (MPFR_TIER=full, --test-threads=1) green on the RC/release SHA. Mandatory; the fast tier alone does not qualify.

## Headroom beyond cosmology (Jeff 2026-10-23 12:33 CT)
Goal: show zf works at its real type limits, not only at the ~10^±120 that cosmology needs. Limits are taken from wip/option-a-w6 (d873c7c): Exponent = i32, EXPONENT_MIN = MIN_NORMAL_EXPONENT = i32::MIN, EXPONENT_MAX = i32::MAX (so |x| goes up to about 2^(2^31) ≈ 10^(±646,456,993)), MIN_PRECISION = 1, MAX_PRECISION = min(i64::MAX/2 + i32::MIN, usize::MAX/4) - 1.
- [ ] (a) Exponent extremes: an MPFR sweep (mpfr emin/emax set to zf's limits, which fit inside MPFR's default range of ±(2^62-1)) for every public function in all 5 rounding modes, with inputs at and next to: the overflow edge (largest finite at p, and exponent EXPONENT_MAX-1 / EXPONENT_MAX), the normal floor (min_positive_normal), the subnormal floor (min_positive), and results that land just across overflow or underflow (for example exp/pow/mul/div/sq/hypot near ±EXPONENT_MAX, and log/sqrt/recip of subnormals). Check the value, the Inf/0/-0 result, and the overflow/underflow/inexact flags. Record the commands and test file.
- [ ] (b) High precision: every public function against MPFR at p = 1024, 2048, 4096, 8192, then 16384 and 65536, and near the largest p the box can allocate (MAX_PRECISION itself cannot be allocated; record the largest p tried). Record runtime per function and per p in a TSV (zf ms, MPFR ms, ratio). Any ratio above about 3x MPFR, or growth worse than M(p)·log p, gets a named cause (series term count, missing binary splitting, Ziv retries, constants-cache rebuild, quadratic multiply) before RC.
- [ ] (b2) Retry budget at high p: use the retry_budget test to show that no function returns PrecisionRetryExhausted at p ≥ 8192 on hard-to-round inputs (the budget there is max(p+16384, 3p)).
- [ ] (c) Accumath cosmology against the 2.0 candidate: run the existing FLRW tests in /tmp/acc-or (engine_core/src/geometry.rs: flrw_scale :1547, flrw_flat :1556, static_redshift :743, tests flrw_friedmann_constraint :3427, flrw_dust_and_exponential :3549) with Accumath pinned to the 2.0 candidate SHA. NOTE: these are symbolic (metric) functions. I found no numeric comoving/luminosity distance-redshift function in Accumath (rg 'comoving|luminosity_distance' only hits geometry.rs/lib.rs). Either add a numeric D_C(z)/D_L(z) check through SoftFloat at z up to 10^4 and parameter ratios of 10^±120 (and 10^±400), or record that the numeric cosmology path does not exist yet.
- [ ] (d) zenith-bound wide-range: run zenith-bound's wide-range suite against the 2.0 candidate at the same exponent extremes and the same p ladder as (a)/(b). Repo path to confirm; not located under /tmp/acc-or.
- [ ] (e) Before the API freeze, decide on the limits flagged in the limits audit (MAX_DEC_SCALE is test-only and should be renamed or hidden; the u32 term-count clamps; whether 1M-bit+ precision is a supported claim).
