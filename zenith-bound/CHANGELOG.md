# Changelog

## Unreleased
- S2 bridge edges (harper, 2026-10-25, gate branch harper/zb-s2-gate = rusty/zb-s2 59c3c82 + harper/zb-pubfn-exempt 59ba053): 8 hand-worked harness tests in tests/harness/slice_bridge/harper_s2_edges.rs (`s2_edges_*`): signed zero, i64::MIN/MAX Mag exponents and the EXPONENT_MAX last bit, 100000-bit last-bit decisions, cancellation on the exact and wide-gap (z = ur ± small) paths of cmp_exact_diff, NaN/±Inf in every argument (also with an Infinity radius), round trips both ways, directed-rounding boundaries (bit 64/65, carry, carry above EXPONENT_MAX), every to_exact_num error detail, Constant::ALL order. Test-only; no src/ change. The approved base failed pubfn-coverage (Constant::* uncovered, slice_bridge.rs COVERS named `constant_ball` instead of `bridge::constant_ball`): COVERS fixed, Constant variants covered by `s2_edges_constant_all_order`. Bridge floor 12 -> 20. No golden rows: the S2 slot holds only bridge_consts.txt (constant_ball); no golden pipeline exists for cmp_exact_num/cmp_exact_diff/conversions. Red-first: "eps vs |x - 1|" flipped to Greater failed s2_edges_huge_precision_last_bit (harper-work/logs/test-red.log).
- pubfn_coverage STALE-EXEMPT (harper, 2026-10-25, Spence): an exempt entry naming no current public item was silently ignored (the list was only a skip filter for UNCOVERED). The check now prints `FAIL STALE-EXEMPT: <name> is in the exempt list but is not a public item` and fails. Selftest +2 cases (24). The shipped exempt list is empty, so no entry was removed. Red-first: pop-logs harper/2026-10-25/1412-pubfn_coverage_selftest.log.
- pubfn_coverage empty exempt list (harper, 2026-10-25, Spence): the exempt load used `grep -v` to drop comments and blanks; on the shipped comment-only list grep selects nothing and exits 1, so under pipefail the ERR trap logged `ERR rc=1 line=119 cmd=sort -u > "$tmp/exempt"` in passing runs. The load is now an awk filter (exit 0 on zero entries, a read error exits 2) and the check prints `exempt list is empty` or `exempt list: N items`. Selftest +7 cases (22): a two-entry list with comments must be read (and the same tree FAILS UNCOVERED with an empty list), and an empty list is reported and logs no ERR (log harper/2026-10-25/1410-full_gate.log). Red-first: pop-logs harper/2026-10-25/1409-pubfn_coverage_selftest.log and 1409-pubfn_coverage_selftest-2580107.log.
- MagArray accessor coverage (harper, 2026-10-25): 6 harness tests in tests/harness/slice_matrix.rs (`s6_magarray_*`) assert exact shape, rank, len, is_empty, as_slice order and into_parts round-trip for rank 0, 1-D, N-D and zero-length dims, plus from_mags InvalidArgument rejections (rank-0 count, length mismatch, [0] with 1, shape overflow). COVERS now names MagArray::shape, rank, len, is_empty, as_slice, into_parts; pubfn-coverage drops them (log harper/2026-10-25/1220-full_gate.log). Red-first: a wrong [1,2,3] expectation (23) failed s6_magarray_nd_accessors_row_major.


- cball follow-up (Spence, 2026-10-25): hand gold row sqr p=53, A = (0 ± 2) + (0 ± 2)i (gen_cball_goldens.py extra_rows; mpmath 1.3.0), where RREF is 0 and the im radius rests on 2·r_re·r_im alone. Not red-first: the new mutant 1524-SQRREFNOXY (sqr_ref without 2·r_re·r_im) is already killed by drawn row 40 (sqr p=2, A = (0 ± 2^-3) + (0 ± 7·2^-23)i; pop-logs cassie/2026-10-25/1142-f2-cb-red); the new row passes and kills it too. 62 gold rows (12 sqr, 42 contained); cball floor 16.

- 16.CB.10 RESOLVED (Spence ruling, 2026-10-25): the sqr tightness bound was mis-specified, not the code. Per part rad(sqr) ≤ 2·(RREF + S2), S2 = r_re² + r_im² (re) and 2·r_re·r_im (im), computed exactly from the gold input radii (test helper `sqr_ref`). `cball_tightness_factors_sqr` is no longer ignored; all 12 gold sqr rows pass (the 10 original rows, including the 3 that failed 2·RREF (pop-logs cassie/2026-10-25/0817-cb10-sqr-try), plus one wide-imaginary sqr row and one zero-mid row with both radii nonzero). Golds extended by that one wide-im sqr row in cball_rows.txt, which kills 1523-SQRREFNOIM; SHA256SUMS updated. CONTRACT 16.CB.10 amended.

- S1 borrowing, option B (2026-10-24, Jeff): `Ball::sqr`, `inv`, `sqrt`, `hypot`, `powi` and `abs` take `&self` (API change; callers
  no longer give up the ball). `sqr` uses `mul_full_prec` and its per-operation `inexact()` flag; `abs` rounds `|m|` as
  `ExactNum::from(±1).mul(m)`. `add`, `sub`, `mul`, `mul_exact`, `div`, `mul_2exp`, `neg` stay by value until zenith-float has
  borrowing correctly-rounded ops (Rusty). No `*_ref`, no new copy; rules_scan baselines unchanged. CONTRACT 15.1, O9, O14-W amended.
  Test: harness `borrowed_ops_reuse_the_same_ball`.
- mutate_slice / full_gate hardening (Spence BLOCK of 9ea44e6, 1742-zb-fgself2/mutself2): the gate runs mutate_slice.sh (--step-timeout and
  --all) under `env -u` MUTATE_FLOOR_ZERO_OK MUTANT_TIMEOUT_S MUTANT_BUDGET_S MUTATE_SETUP_S MUTATE_SH, so a leaked selftest knob cannot
  change a real run. A result is only an UNINDENTED `KILLED|SURVIVED-EQUIVALENT|SURVIVED|INVALID <id>` line for the --only id (a forged
  indented line is suite output; a result for another id fails). A runner SIGKILLed before its limit (exit 137) is `KILLED-OOM <id>` and
  fails, not TIMEOUT. A --step-timeout failure writes its whole output to the step log (`| tail -1` lost the cause of the one-off
  1706-zb-mslice-9ea44e6-selftests failure). full_gate_selftest has 118 cases and logs the count (+6 from main's 25-step gate; +1: a stale MUTATE_REPO/MANIFEST/EQUIVALENT/FLOOR/TIMEOUT/BUILD_TIMEOUT/TMP/PERSIST must not reach mutate.sh, so the gate unsets every env mutate.sh reads except MUTATE_TARGET, which it sets; +19: leaks, duplicate/unknown/forged results, unreadable
  floor/list, runner rc, OOM, --slices cause).

- full_gate mutants step per slice (Spence-approved design, 2026-10-24): the step runs tools/mutate_slice.sh --all, which runs each slice
  (`mutate.sh --slices`: central + manifest.d) with `mutate.sh --slice` and checks the result lines itself instead of trusting the runner's
  exit code. It FAILS on an INVALID mutant, a SURVIVED mutant with no equivalent entry, a slice under its floor (floor.txt for central,
  floor.d/<slice>.txt), a manifest id with no result (`INCOMPLETE: 257 results for 268 ids ...; missing: <ids>`, per slice and over the
  whole manifest), a duplicate/unknown result id, or a slice past its timeout (`TIMEOUT ...; never ran: <ids>`). Per-slice timeout =
  ids x MUTANT_BUDGET_S + MUTATE_SETUP_S; the gate step's timeout is `mutate_slice.sh --step-timeout` (sum + 600s), and the step FAILS if
  it cannot be computed. Budget source: a real tiny-slice run (pop-logs harper/2026-10-24/1644-zb-mslice-measure, f854a9d: 4 mutants,
  one per suite, fresh MUTATE_TARGET) had a max of 14s per mutant and ~1s setup; 1520-zb-mutq-merge-full-gate ran 184 mutants in 4m19s.
  BUDGET_S = 60 (14s x ~4), SETUP_S = 300. mutate.sh gains --slices, --slice SLICE and a `TIME <id> <s>s` line per mutant.
  full_gate_selftest has 82 cases (+10: timeout, missing results, INVALID, unexplained survivor, below floor, each against a runner that
  always exits 0). Spence chose this mutate_slice.sh over fba8ec9's; its ID_REGEX is folded in:
  `mutate_slice.sh SLICE ID_REGEX` runs only the matching ids, prints `PARTIAL run (ID_REGEX=...): floor/completeness not checked` and
  never a full-pass line; a regex matching no id fails; `--all` rejects a regex, so full_gate never runs partial.
- mutate_slice per-mutant timeout 60s (measured max 14s for 178-U-underflow-zero in harper/2026-10-24/1644-zb-mslice-measure, x ~4):
  each id runs as `mutate.sh --slice S --only ID` under that limit; a mutant past it prints `TIMEOUT <id>` and the slice fails (the slice
  continues and lists every timed-out id). Previously a hang ran to mutate.sh's 20m test timeout and counted as KILLED. The per-step
  total is unchanged.

- goldens_regen.sh SHA256SUMS line check quiet/loud (2026-10-24): the malformed-line `grep -v` returns 1 when every line is valid,
  which tripped `ERR rc=1 line=65` on a clean run. rc 1 is now normal (clean rc 0, no ERR); a real malformed line still fails with
  its text, and a grep rc 2 fails with "grep failed (rc 2)". goldens_regen_selftest has 25 cases (+6).

- mutate.sh inrange forces base 10 itself (Spence nit on b5ac93f): `10#` on the id and the range bounds, so 0213 is 213 (not octal
  139) and 0219 is not an arithmetic error. mutate_include_selftest has 36 cases (+1).

- mutate.sh inrange quiet/loud (2026-10-24, Spence zb-s2 ceb4024 log 1422-zb-s2-const-steps: 368 `ERR rc=1 line=92`): inrange now
  returns 0 with empty output on a miss (an id outside every slice range is normal), and rc 2 + a message on a malformed range row.
  Call sites no longer mask it with `|| true`: a nonzero rc trips the ERR trap and fails the run. mutate_include_selftest has 35
  cases (+7). Note: ceb4024 does not contain 7881073; rebasing zb-s2 on main also quiets line 92.

- popguard_selftest no-guard-host goldlib case (2026-10-24, Spence follow-up): asserts the exact WARNING line and exit 0. It was
  rc 1 because the sandbox HOME had no oracle venv (goldlib failed its oracle check, not the guard); the sandbox now links the real
  pinned venv as a stand-in and runs goldlib with it.
- Pop guard bypass closed (2026-10-24, Spence blocked 7881073: POPGUARD=/nonexistent turned the cap off silently): goldens_regen.sh
  and goldlib.py REFUSE (exit 2) when POPGUARD is missing/not executable on a host that has ~/pop-logs/bin/popguard.sh; a host with
  no popguard.sh runs with a loud stderr WARNING. popguard.sh (pop-logs 6d1fc93) now also requires the scope's memory.max to be a
  number <= the 6G cap. popguard_selftest never skips on Pop (fails instead) and has 15 cases (the red commit said 16; 15 is right).
  Red: harper-work/popguard/zb-bypass-red.log; green: zb-bypass-green.log.

- mutate.sh quiet ERR trap (2026-10-24, Spence: 1207-zb-popguard-full-gate1 logged 184 `ERR rc=1 line=92` inrange misses and one
  `ERR rc=101 line=214` per killed mutant): the inrange calls take `|| true` (no owning range is normal), the suite run captures its rc
  explicitly (`rc=0; (...) || rc=$?`), and an unexpected rc (INVALID branch) is re-raised with `(exit "$rc")` so the ERR trap still
  records it. mutate_selftest 34 cases (+2: logged run has no ERR line; unexpected rc 101 still logs ERR and fails).
  Red: harper-work/popguard/zb-errquiet-red.log; green: zb-errquiet-green.log.

- Pop direct-run guard (2026-10-24, 11:49 incident: a bare 58 GB goldens python run global-OOMed Pop): tools/goldens_regen.sh
  (before its first generator) and scripts/goldlib.py (at import, before the oracle check, so every mpmath generator) call
  ~/pop-logs/bin/popguard.sh and exit 3 outside a popslot/poprun memory scope (6G per job, swap 0); hosts without popguard.sh
  run unguarded. New tools/popguard_selftest.sh (8 cases, in the goldens-regen step).

- goldens_regen (2026-10-24): a slice SHA256SUMS holding only '# EMPTY' is skipped cleanly (no ERR trap lines); '# EMPTY' with a *.txt next to it, with entries, or with a malformed line FAILS. goldens_regen_selftest 19 cases. Red: 1200-goldens_regen_selftest, 1200-goldens_regen_selftest-2200915; green: 1201-goldens_regen_selftest
- full_gate (2026-10-24): gate-selftest now runs tools/mutate_selftest.sh (and requires cargo, git, timeout). Red: 0915-zb-red-gate-skips-mutate-selftest
  (step PASS with a planted break); 0915-zb-red2 FAIL with it; 0917-zb-green-gate-mutate-selftest PASS, 28 cases.

- full_gate selftest (2026-10-24): its nested gate runs unset `GATE_STEPS`, so `GATE_STEPS=gate-selftest tools/full_gate.sh` no longer
  fails with "stripped run reported too few steps" (the outer step list leaked into the fixture gates).

- runlog selftest joins the gate (2026-10-24): poplog 729c861 (Rusty) prunes a row only when the pid is gone (ESRCH), so the EPERM
  pid-1 row now survives and `runlog_selftest.sh` passes all 50 cases; full_gate.sh's gate-selftest step now runs it.

- No /tmp in tools (2026-10-24, Jeff: a reboot wiped /tmp and with it a dependency pin and work dirs): new `tools/scratch.sh`
  (`zb_scratch`, sourced by runlog.sh) gives the temp parent, `$ZB_SCRATCH` or `$HOME/.cache/zenith-bound/scratch`. full_gate.sh
  and deep_mag.sh logs (`T`), mutate.sh `MUTATE_TMP`, the fuzz_deep/deep_mag target dirs, runlog.sh, fresh_test_run.sh,
  goldens_regen.sh, digests.sh and every selftest's mktemp now use it. Only the cross-agent locks /tmp/pop-deep.lock and
  /tmp/pop-gate.lock stay. New `tools/no_tmp_selftest.sh` (5 cases, in the gate-selftest step) fails on a hardcoded /tmp in a
  code line of tools/*.sh.

- mutate.sh (2026-10-24): the fixed persistent worktree (`zb-mutate.persist`) now goes under `MUTATE_PERSIST`, by default the parent of
  `MUTATE_TARGET`, instead of `MUTATE_TMP` (default /tmp, wiped by the 10/23 reboot). mutate_selftest gains 2 cases (28): an unwritable
  `MUTATE_TMP` does not matter in persistent mode, with and without `MUTATE_PERSIST`.

- runlog selftest (2026-10-24): the other agent's RUNNING row now uses a live `sleep` pid started by the selftest (it used pid 1, and
  poplog b6e2dcf prunes a row when `kill -0` fails, which for pid 1 is EPERM, so the row was wrongly pruned). A second row with pid 1
  asserts that a live pid this user cannot signal is kept; it FAILS until poplog prunes only on ESRCH (`poplog.sh:65`, reported to
  Rusty), so `runlog_selftest.sh` stays out of the gate. The exec scan no longer flags a printf line writing a fixture script. The
  selftest counts its cases (50).

- Spence nits (2026-10-23): `add_up` gets a gap-65 row (rounds up to TOP+1); the test is renamed `add_up_at_gaps_0_63_64_65_200`
  (still 6 tests in `mag_rows.rs`, test floors stay 379). A mutant floor value of 0 now FAILS BADCOUNT ("must be 1 or more"); only the
  mutate selftest fixture may use 0, via MUTATE_FLOOR_ZERO_OK=1. `tools/mutate_selftest.sh` gains that case (26 cases).

- Selftest case counts (2026-10-23): `tools/selftest_cases.sh` (`st_init`/`ck`/`st_done`) makes 8 more gate selftests count the cases they
  actually run and FAIL on a mismatch: digests 22, fresh_test_run 26, full_gate 72, pubfn_coverage 15, rules_scan 284, tightness 33,
  toolchain_check 14, wildcard_enum 11. `runlog_selftest.sh` is not covered yet: it already FAILS at HEAD ("another agent's row was removed").

- Mutant floor (Paul, 2026-10-23): a full `tools/mutate.sh` run (the gate's mutants step) FAILS if killed or killed+equivalent is below
  `tools/mutants/floor.txt` (159 killed / 176 killed+equivalent now; MUTATE_FLOOR overrides). `--only` skips the floor but still FAILS with
  0 mutants. Missing, garbage or CRLF floor rows FAIL. `tools/mutate_selftest.sh` and `tools/test_floor_selftest.sh` now count their
  cases and FAIL unless the count matches (25 and 13), so a selftest that runs nothing cannot read ok.

- L0a-2 edge rows (TEST_PLAN 8.1): `tests/harness/mag_rows.rs`, 6 tests: sqrt of 1, of powers of two and of all-ones mantissas (odd and
  even exponent), mul of all-ones and of powers of two, and shift distances 0/63/64 in `Mag::finite`, `add_up` and `sub_lower`. New
  mutants 181-S-sqrt-nosticky and 182-F-finite-minexp. Test floors raised 373 -> 379 (build-debug and build-release).

- Test-count floor (Paul, 2026-10-23): new gate step `test-floor` runs `tools/test_floor.sh` on the build-debug and build-release logs.
  It fails if the summed `test result:` passed count is below the step's floor in `tools/test_floor.txt` (373 each), if there is no summary,
  or if fewer summaries than started test binaries appear (a truncated log; the 2052 gate on 7191a1f passed with 599 of 738). Raise a floor
  with the tests that add to it; lowering one needs a CHANGELOG note. Fails closed (Spence on 92ac25f): `set -euo pipefail`, every parsed count must be a plain integer of at most 18 digits
  (else `FAIL BADCOUNT`), the log must be readable (`-r`), and only `passed` counts toward the floor. Self-test: `tools/test_floor_selftest.sh`
  (now also unreadable, garbage, huge-count and passed+failed cases). A CRLF floor row FAILS with a message naming the CR (Spence nit on e61f0ad). poplog is unchanged (Rusty owns its guard).

- `tools/mutate.sh` refuses a dirty tree (modified or untracked files; exit 2, `mutate: FAIL dirty tree`) and prints
  `mutate: testing committed HEAD <sha>` (Spence, 2026-10-23). Mutants run on a worktree of the committed HEAD: 2102-zb-mut180 reported
  180-C-cmp-zero-eq SURVIVED because it tested HEAD 7191a1f while the mag_sort tests were still uncommitted. Self-test: `tools/mutate_selftest.sh`
  (the fixture now commits each manifest edit). `mag_sort.rs` `key()` lists every invalid flag/field combination instead of a catch-all.

- L0a-2 `Ord` permutation sort (TEST_PLAN 8.1, red-first): `tests/harness/mag_sort.rs` sorts a 12-value set (both zero constructors, which build the one canonical zero, `i64::MIN`, the zf bounds, a duplicate, `u64::MAX` at `i64::MAX`, `Infinity`) and checks the hand-ordered result under all 720 permutations of a 6-value subset and 500 seeded SplitMix64 Fisher-Yates shuffles of the whole set. New mutant `180-C-cmp-zero-eq` (zero compares Equal to a finite value, suite harness); red: 3 of the 4 new tests fail with it applied. The library passes.
- mag_carry follow-up (Spence's review of 2d12c51): the inline saturating `add_up` in `tests/harness/mag_carry.rs` is now mutant `179-C-add-nocarry` (`successor` saturates at `u64::MAX` instead of carrying, suite harness), the inline function and its `_ =>` arm are gone, and the header cites CONTRACT 3.3 for the overflow rule.
- L0a-2 carry across the exponent (TEST_PLAN 8.1; red-first): `tests/harness/mag_carry.rs` checks `add_up` (both orders) on hand-derived rows where the mantissa sum, or its round-up, reaches 2^64 (`u64::MAX` at e plus 2^63 at e - 63 or e - 100 gives 2^63 at e + 1) plus the doubling control, at e = 0, zf MAX - 1, zf MAX, i64::MAX - 1 and i64::MAX (where the carry gives Infinity); an `add_up` that saturates instead of carrying is rejected. Red: the planted saturating `add_up` survived with the comparisons stubbed out. The library passes.
- Forced inputs for every L0a-5 Mag property (fuzz_mag_props.rs, handed over by Spence for this change; red-first): besides its random band, every property now runs on every ordered pair of `generators::forced_prop_exponents()` (the zf bounds `MIN`, `MIN + 1`, `MAX - 1`, `MAX`, then `0`, then +-256, +-405, +-1000, +-1024: 13 values, 169 pairs) as its first two exponent draws (`Src::forced`, `run_forced`; a forced `draw_case` keeps its forced second exponent; `sqrt_exact` uses `e / 2` so the square sits near `2^e`). A `mul_up` wrong only at (zf MAX - 1, zf MIN + 1) and (405, -1000) survives 20,000 random cases and is caught by the forced pairs. Red: 18 failed (empty forced list). The real library passes every forced pair; no digest changed. The `forced_exponents` doc now says the forced cases run after the random ones (as they have since 0d30dec).
- Forced zf-bound and cosmology-scale inputs (step 4; red-first): `generators::forced_exponents()` lists the zf bounds `MIN`, `MIN + 1`, `MAX - 1`, `MAX`, then the cosmology-scale exponents `COSMOLOGY_EXPONENTS` (256, 405, 1000, 1024, 4097, 2^16, 2^20 - 1, 2^20) each with both signs. `fuzz_mag` runs every pair of them x mantissas {0.1b, 0.1^64b} at prec 53 (1600 forced cases) after the random cases of every non-replay run; a wrong op that is wrong only at those exponents is now caught. The edge table gains 5 rows (`max-exponent-1,1ulp`, `min-exponent+1,1ulp`, `2^405,1ulp`, `2^-405,1ulp`, `-2^-405,r=0`; 20 -> 25). Red: 6 failed (stubbed list, rows absent). No digest changed.
- zf exponent bounds (Rusty: zf 2.0 `Exponent = i64` with range +-(2^62 - 1); red-first): `ZF_EXPONENT_MIN/MAX` in `tests/common/generators.rs` are now `-(2^62 - 1)` / `2^62 - 1` (were `i64::MIN/MAX`), defined once there with a note that Rusty's zf i64 SHA brings the real zf constant names. The L0b boundary tests moved from `tests/ball/boundary.rs` to the harness (`tests/harness/ball_boundary.rs`, using those constants) and run at `MIN`, `MIN + 1`, `MAX - 1`, `MAX`; the Mag ulp-underflow checks now sit where the underflow really happens inside the range (precisions above 2^62: `ulp_mag` clamps from 2^62 + 3, `half_ulp` from 2^62 + 2), and a new test checks an exponent sum and difference at both edges fit in i64 and give exact `mul_up` / `div_up` powers. Red on the old `i64::MIN/MAX` bounds (4 of 6 failed). Consequences: the MID edge table has 13 entries (min - 1 and max + 1 are representable), a MID gap one tick past a bound is a raw value and two ticks is an error, MID Uniform bit lengths stop at 62; `GENERATOR_VERSION` 4 -> 5 and `BALL_DIGEST` re-pinned (MAG digests unchanged). Mutant 178-U-underflow-zero now runs in the harness suite.
- L0b boundary tests (`tests/ball/boundary.rs`, red-first): at mid exponents `i64::MIN`, `MIN+1`, `MIN+2` and precisions 0, 1, 2, 53, 4096, `u32::MAX`, `ulp_mag` and `half_ulp` round an exponent field below `i64::MIN` outward to the smallest positive `Mag` (Spence, 2026-10-23 7:51 PM CT; CONTRACT 3.2), the clamped radius is never reported below half an ulp, a ball with it is informative and not exact, and `is_informative` turns false exactly at `ulp_mag(i64::MAX, 0)` = Infinity. src already does this; red shown with underflow mapped to zero, kept as mutant 178-U-underflow-zero. Finding for Spence: the pinned zf (fe5c15d) still has an i32 `Exponent`, so a `Ball` midpoint cannot yet sit at exponent `i64::MIN`; the mid-side boundary is tested through the `Mag` ulp terms only until zf 2.0 is pinned.
- `tools/mutate.sh` checks the manifest before running any mutant: every mutant id appears exactly once and every patch file a row names exists; otherwise it fails at once (`mutate: FAIL duplicate manifest id ..` / `.. names missing patch ..`). Prompted by the 092-094 id clash between harper/zb-s7-fuzz and spence/zb-l0b, which git only caught as a text conflict. Red-first in `tools/mutate_selftest.sh`.
- `BallArray` shapes with a zero extent are empty arrays even beside extents whose product would overflow (`[usize::MAX, 2, 0]` was `InvalidArgument`). WR-ND tests: ranks 0..6, empty dims on every axis, overflow at every rank and axis.
- Container/LA op signature rule frozen (API freeze F8; CONTRACT 15.6): explicit `prec: Precision`, no rounding mode, `BoundResult` return.
- Allocation failure is `Resource`, never an abort, for storage zenith-bound sizes from caller input (API freeze F7; CONTRACT 15.5). `BallArray::zeros` with a shape too large to allocate returns `Resource`.
- `Tri` is documented as exhaustive (API freeze F2): `{True, False, Unknown}` is the complete set; no `#[non_exhaustive]`. CONTRACT 15.2.
- `BoundError::NotConverged { op, detail }` (Paul, 2026-10-23): an iteration out of budget (e.g. `EigOpts::max_sweeps`) returns it, never a silently loose answer. CONTRACT O11-W.
- `Enclosure::contains` for `Ball` is tight (Paul's ruling, 2026-10-23): when `rad` is a zenith-float value, False is returned exactly when `|mid - x| > rad` (e.g. `1 ± 1` vs `-2^-300` is now False, was Unknown). Unknown remains only for NaN-producing input, an overflowing working precision, or a `rad` outside zenith-float's exponent range. The lower-bracket fallback is gone. Exact-order property tests in `tests/mag/contains_exact.rs`.
- Options: `SolveOpts { refine_steps }` and `EigOpts { max_sweeps, cluster }`, both `#[non_exhaustive]` with `Default` (2; 64, true); `DomainPolicy` is `#[non_exhaustive]` ({Strict, Clip}). CONTRACT 10.1.
- Precision headroom (Jeff, 2026-10-23 12:33 PM CT, red-first): `PrecSource::Uniform` draws any `u32` precision (CONTRACT 8.1), not just `1..=4096`: one draw in 8 from `PREC_EDGES` (1, 2, 3, `MAX_PREC` and its neighbours, both sides of 2^31, the top three of `u32`), the rest log-scale over bit lengths 1..=32, still one output. `PrecSource::Operand` accepts any `op >= 1` (±1 kept unless it hits 0 or passes `u32::MAX`). Before, the fixed-seed run topped out at 4094. `GENERATOR_VERSION` 4; no pinned digest moved (no suite draws `Uniform` precisions yet).
- Cosmology-scale exponents (Jeff, 2026-10-23 12:20 PM CT, red-first): the Uniform exponent band draws through the new `wide_uniform_in` (still one output per draw): half linear over the domain as before, half log-scale (bit length 0..=63 uniform, then uniform within it), so `2^8 <= |e| < 2^20` (1e±122 is about 2^±405) and every bit length are hit; before, a fixed-seed run of 10,000 drew 0 such non-power-of-two exponents. Band weights and counts are unchanged; the shape inside the Uniform band changed (for Spence, TEST_PLAN section 7). `GENERATOR_VERSION` 3; re-pinned: `EXP_DIGEST_MAG`, the ball-descriptor digest, the W1 random tallies and the `generators` row of `digest_expected.txt`.
- Mid exponents span `i64` (Jeff, 2026-10-23 12:36 PM CT: zenith-float 2.0 `Exponent` is `i64`; red-first): `ZF_EXPONENT_MIN/MAX` are `i64::MIN/MAX`, so `ExpDomain::MID` draws cosmology-scale and extreme mid exponents. The one-tick-outside raw values and `exponent_to_i32` are gone (nothing lies outside `i64`; a MID gap past an extreme is an error). `GENERATOR_VERSION` 2 (`seeds.rs`, `goldlib.py`); the ball-descriptor digest is re-pinned, the MAG band counts and `digest_expected.txt` are unchanged. CONTRACT 8.3 and O3 still say `i32`: for Spence/Rusty.
- Test grid headroom (Jeff, 2026-10-23 12:33 PM CT, red-first): `points::Term` exponents and `radius_value` are `i128`, so the L0b P1 grid of any `i64` mid exponent at any `u32` precision is exact. Before, a 1-ulp radius at `mid_exp = i64::MIN` (or any precision pushing `mid_exp - prec` below `i64::MIN`) failed with "ulp exponent underflow". New self-test `grids_exist_at_the_i64_exponent_extremes_for_every_u32_precision`.
- RUNLOG_PLAIN needs a selftest-only marker (Spence follow-up, non-blocking, 2026-10-23, protected `tools/`, red-first): `runlog.sh` (`runlog_wrap`, so `full_gate.sh` and every wrapped script and one-shot run) exits 2 with a message naming `RUNLOG_SELFTEST` when `RUNLOG_PLAIN=1` is set without `RUNLOG_SELFTEST=1`. Only the self-test fixtures (`full_gate_selftest.sh`, `runlog_selftest.sh`, `fuzz_deep.sh --selftest`) export the marker. New red-first cases in `full_gate_selftest.sh` and `runlog_selftest.sh`.
- One target dir per worktree (Paul, 2026-10-23 1:02 PM CT, protected `tools/`, red-first): `tools/fresh_test_run.sh` fails if `CARGO_TARGET_DIR` is unset or the dir is owned by another worktree (marker `$CARGO_TARGET_DIR/.worktree-owner`, written by the first worktree to use it; `--claim` runs only this check). New gate step `target-dir` runs it before any build; `mutants` no longer falls back to `/tmp/zb-target`. Selftest: two worktrees on one dir fail (SHARED, naming the owner), an unset `CARGO_TARGET_DIR` fails. CONTRIBUTING updated.
- Spence review fixes for c471eca..7b64d57 (2026-10-23, protected `tools/`, red-first):
  - No silent pass (BLOCK): `full_gate.sh` keeps a running `gate_rc`, writes `done` only when every step passed (`failed` otherwise) and ends with `exit $gate_rc`. `runlog.sh`: a missing or unloadable poplog.sh is a hard failure (exit 2); `RUNLOG_PLAIN=1` is the explicit fixture-only plain mode. `runlog_step --tests` enforces 0-tests=FAIL in the zb layer too (with or without poplog). New selftest cases: missing poplog, failed step exits non-zero and writes no `done`, plain-mode zero tests.
  - `fresh_test_run.sh` (Spence's ruling): the age test is gone; a test binary must lie inside `$CARGO_TARGET_DIR`. The 0-binaries check and the manifest-scan gate step stay. The gate may reuse a target dir again; per-job subdirs are optional. Usage is now `tools/fresh_test_run.sh CMD...`.
  - `mutate.sh`: 0 mutants is a failure ("0 mutants ran"), so the gate's mutants step needs at least one mutant; the selftest's empty-manifest case now expects failure.
  - `runlog_wrap` sets `pipefail` (red case in `runlog_selftest`). Pipelines made safe under it: `mutate.sh`'s equivalence lookup is one `awk`; piped `grep -q` became `grep -c > /dev/null` so the writer can't hit SIGPIPE (rc 141); commands whose exit the selftests deliberately ignore (`scan --count-only`, `tightness check/report`) are wrapped in `{ ... || true; }`.
  - `fuzz_deep.sh`: every shard and the replay must report `1 passed` (a renamed or filtered-out target is a failure); `timeout` on the release build (`FUZZ_DEEP_BUILD_TIMEOUT`, default 30m) and each shard (`FUZZ_DEEP_SHARD_TIMEOUT`, default 6h); default deep-lock wait 3600 s, announced; targets come from `tools/fuzz_deep_targets.txt` (new) so Spence's port adds `fuzz_mag_l0a5` there. Shard logs move to `OUT/<target>/shard-N.log`. Selftest fixtures run with `RUNLOG_PLAIN=1`.

- Gate rerun fixes (Jeff, 2026-10-23 10:33 AM CT). These touch protected `tools/` files, so they need Spence's review:
  - `tools/runlog.sh` v3 relies on poplog's FAIL verdict and exit code. The `/dev/null` stdout workaround and the zb-fail markers are gone. A nested full_gate under poprun now exits nonzero, shows FAIL and writes step output into the parent log. There is a red-first selftest case for this.
  - Harness goldens resolve at runtime through `common::goldens::repo_root()`, not compile-time `env!("CARGO_MANIFEST_DIR")`. A binary built in one worktree used to read another worktree's goldens. The new gate step `manifest-scan` bans the `env!` form.
  - A wrapped tool registers its exit cleanup with `runlog_on_exit`, not a raw `trap ... EXIT`, and does not `exec` a command. Both used to skip poplog's finish, which left logs without a footer and RUNNING.tsv rows stuck in RUNNING. Self-tests call `runlog_detach_children` so their deliberately failing fixtures do not FAIL the self-test's own run. A lint in `runlog_selftest.sh` enforces both rules, red first.
  - The gate's build-tests step runs `cargo build --tests`. poplog treats any `cargo test` command as a test step, so the old `cargo test --no-run` read as zero tests and FAILed (first seen in the 11:15 AM CT gate). `full_gate_selftest.sh` now rejects the old form.
  - New `tools/fresh_test_run.sh` is used by build-debug and build-release. It fails if a step runs 0 test binaries, a binary from outside `CARGO_TARGET_DIR`, or a binary older than the gate start. As a result, the gate needs a fresh per-job target dir, `/tmp/target-harper/<job>`.

- Run logs now go through the shared Pop logger (Jeff, 2026-10-23 9:59 AM CT; Paul and Rusty): `tools/runlog.sh` sources
  `~/pop-logs/bin/poplog.sh` instead of writing logs and `RUNNING.tsv` itself, so logs, the canonical 8-column `RUNNING.tsv` row (lock
  `RUNNING.tsv.lock`) and 30-day pruning (`poprotate`) are poplog's. zb keeps a thin layer: step output still on stdout, a failed or zero-test
  step returns nonzero and fails the run, no second log for a wrapped script inside a wrapped run or inside `poprun`, heartbeat kept off stdout.
- Run logs (Jeff, 2026-10-23 9:42 AM CT): new `tools/runlog.sh`, sourced by every `tools/*.sh` (`runlog_wrap JOB "$@"`), and usable one-shot
  for verify and test runs: `tools/runlog.sh JOB [--tests] -- CMD`. Each run logs to `~/pop-logs/<agent>/<date>/<HHMM>-<job>.log` (header,
  per-step timing, exit code and cargo test counts, 5-minute heartbeat, PASS/FAIL footer with warnings) and holds a row in
  `~/pop-logs/RUNNING.tsv` under `flock ~/pop-logs/.running.lock`. A `--tests` step that runs zero tests fails. In `full_gate.sh` every step is
  a run-log step (build-debug and build-release are `--tests`), and a wrapped run now exits 1 when any step failed, even if the script exits 0.
  Self-test: `tools/runlog_selftest.sh`.
- Goldens: the generators now run only in the isolated oracle venv `~/.local/oracle-venv/bin/python` (mpmath 1.3.0, sympy 1.14.0, gmpy2 2.3.2,
  python-flint 0.9.0, the zenith-float `scripts/oracle-requirements.lock` pins), replacing `/usr/bin/python3` with Ubuntu mpmath 1.2.1 (Jeff,
  2026-10-23). `scripts/goldlib.py` checks the four exact versions, that `sys.prefix` is the venv and that no `~/.local` path is on `sys.path`
  at startup. The golden files' provenance header line changed; every value is byte-identical.
