# zf/zb precision matrix: design v4 (Spence third review, items 1-7, on v3)

Agent: Harper. Date: 2026-10-30. Edge table: `precision-matrix-edges.csv` (24620 rows, 17 columns; the first 9 lines start with `#` and carry the pins below; the CSV header is the first line without `#`).

## 0. Pinned inputs (item 11)
The generator `tools/gen_edges.pl` reads only the in-branch copies below. It checks every sha256 before reading and stops with `PIN MISMATCH` if one differs. It also stops if `zenith-float-num/{src,tests}` or `zenith-bound/{src,tests}` differ from the source commit. The originals in ~/paul-work and ~/spence-work were copied unedited.

| input (docs/test-plan/inputs/) | copied from | sha256 |
|---|---|---|
| zf-inventory.csv | ~/paul-work/zf-audit-2026-10-30/inventory.csv | ba79198eb3486944e3b43a898389b81a7ce8fb69e4f67450e49841e8fdac1d59 |
| zb-inventory.csv | ~/paul-work/zb-audit-2026-10-30/inventory.csv (section 1 used) | b65c2bad80431500ab6280eb6ca545194407a35e9ae2f9719b5a6e5a38bdcd65 |
| r4-required-results.csv | Spence R4 (via ~/harper-work/pm-review-inputs) | 2746dc4d05cca03d07db6cd8db600e86369a7f33b450c6b679b1d0c1bf0e8c6d |
| spence-gen_r4.py | ~/spence-work/pm-review/gen_r4.py (G1-G3 lists, lines 6-19) | aa4d614617d086a49a2d625adbef1b516a7cac7ea37de58e391d40df2b51243f |
| FIXLIST.csv | ~/paul-work/zfzb-plan-2026-10-30/FIXLIST.csv | 69391548209c157a3e8cc618abc3b474633b4d9cac89356f88ae28eb30d90887 |
| fixes_after_1.0.12.txt | ~/paul-work/zf-audit-2026-10-30/fixes_after_1.0.12.txt | 03085a0acb3d5d0e51d8e8df13b3c4ad30d7843152d3d32024d6899c9f530c82 |
| zf-audit-SUMMARY.md | ~/paul-work/zf-audit-2026-10-30/SUMMARY.md | 81cf6bb72636dfab03b2a25121724db4b37e3c16e0b896c30269db89be4122f1 |

- Source commit: 7955836a845bd6419dc41d21eda6e739e3b622c5 (goldens, `tests/nan_error_kinds.rs`, `defs.rs`, `error.rs`).
- **Source branch (pmrev3 item 5):** SRC_COMMIT 7955836 is on **paul/zf-bound-merge** (581 commits ahead of zf main), not on zf main (dc6bc87). It **must be re-pinned** (goldens, nan_error_kinds.rs, defs.rs, error.rs, inputs/SHA256SUMS, regenerate, lint, pm-coverage) when that branch merges.
- Oracle env: ~/team-env/oracle-2026-10-30 (mpmath 1.3.0, sympy 1.14.0, gmpy2 2.3.2, python-flint 0.9.0). Rust MPFR 4.2.2 / GMP 6.3.0 via rug 1.20.1 / gmp-mpfr-sys 1.6.8 (Cargo.lock @7955836).
- `inputs/SHA256SUMS` lists the same hashes. The hashes of the generator, lint and CSV are in the Harper log for this commit.

## 1. Verdicts and gating (items 4, 9)
- The only verdicts are PASS, FAIL and NA_PROOF. There are no waivers, and an escalation counts as FAIL.
- Rows: PASS|FAIL zf 18419, zb 2186; NA_PROOF zf 3122, zb 893.
- **Gating plan (item 9):** pm-coverage FAILs a cell when either oracle column of an executed row still reads `sha256:UNGENERATED`. Any such FAIL blocks the release candidate. Today 19742 executed rows carry an UNGENERATED pin (the lint counts them), so no RC can pass until PM-3 writes and pins those golds in SHA256SUMS. pm-coverage also FAILs an NA_PROOF row whose proof file is missing, or whose `Second-agent sign-off:` line reads UNSIGNED.
- **Tie rows (item 4):** every zf tie row now runs (645 rows), except for families in Spence's G1-G3 (accessor, construct-copy, exact-int-modular; spence-gen_r4.py:6-19). Those 74 rows are NA_PROOF. Each has a committed placeholder at `pm-proofs/zf/<key>.md` reading `Second-agent sign-off: UNSIGNED`, so pm-coverage FAILs them until a second agent signs. Families with no RoundingMode are tied at the fixed mode they call internally.

## 2. Family keys (R1)
- Families are keyed by `Type::name` from impl_context. There are 719 zf and 116 zb families, plus `DEFAULT_RANDOM_SEED`.
- `Ball` and `ComplexBall` in zf rows are the zf types in zenith-float-num/src/ball.rs. They are not zenith-bound types.

## 3. Rounding modes
There are 7 modes. ToOdd: the oracle computes RNDN at p+1 and flips exact ties itself. None = RNDZ of the exact value at rp(p,64) (R6).

## 4. Precision axis and MAX_PRECISION (item 2; closes old Q1)
- The base axis is p ∈ {32,64,128,256,512,1024}, with 96 heavy points, and p=1,2,3 run in all 7 modes. p=0 and p=MAX+1 are refused with InvalidArgument.
- `prec-MAX-x86_64` (p=2^60-1) and `prec-MAX-i686` (p=2^30-1) are separate rows, 387 each. Paul ruled that 2^30 bits can be allocated on i686.
- The i686 rows run **only in the heavy tier, in their own shard `pm-max-i686`**. Each row has a 10-minute budget, and a timeout is FAIL.
- The required result is the correct value, or a prompt refusal: `MemoryAllocation` or `PrecisionRetryExhausted` for zf, and `Resource` for zb. A wrong value or a timeout is FAIL, and **every FAIL becomes a FIXLIST row**.

## 5. Error variants (items 8, 10)
- The lint (`tools/lint.pl`) parses the enums from source on every run:
  - zf `Error` (zenith-float-num/src/defs.rs): DivisionByZero, ExponentOverflow, InvalidArgument, MemoryAllocation, PrecisionRetryExhausted.
  - zb `BoundError` (zenith-bound/src/error.rs): DimensionMismatch, Domain, InvalidArgument, NotConverged, NotPositiveDefinite, NotSymmetric, Representation, Resource, Singular.
- The lint rejects:
  - any `Error::X` or `BoundError::X` where X is not in the row's crate enum
  - any bare variant name that exists only in the other crate's enum
  - `BoundError::` in zf rows, and `Error::` in zb rows
  - any "or documented" wording
- A negative test confirmed the lint catches `BoundError::PrecisionRetryExhausted`, `Error::Bogus` and a zf `Resource`.
- **Trig cap in zb (item 8):** zenith-bound has no trig functions (no `fn sin|cos|tan` in zenith-bound/src). If zb reaches the zf cap through the bridge, it surfaces as `Err(BoundError::Resource{..})`: zf `PrecisionRetryExhausted` maps to `Resource` at mag/convert.rs:209, arith.rs:595 and bridge.rs:294. The cap is `trig_reduce_extra` (zenith-float-num/src/ops/util.rs:92-95), `1<<16` with usize saturation, so it is the same on host and i686. The generator writes Resource for any zb trig row. The v2 "zb" trig-huge rows were the zf `Ball::sin/cos` and `ComplexBall::sin/cos` (ball.rs:105,112,352,370). They return `Self` with a NaN midpoint, so `Error::PrecisionRetryExhausted` is the correct zf kind there.

## 6. Thresholds and operand classes (item 6)
- The threshold classes (ovf-max-finite, ovf-max-plus-half-ulp, unf-min-normal, unf-min-subnormal, unf-half-min-subnormal-tie) are kept, with the T1-T12 table and the i64 extremes for scalb, ldexp and ilogb.
- Three operand classes are restored, each with 835 rows (one per family):
  - `huge-exp`: an operand exponent in [EXPONENT_MAX-64, EXPONENT_MAX]
  - `tiny-exp`: an operand exponent in [EXPONENT_MIN, EXPONENT_MIN+64]
  - `subnormal-like`: unnormalised raw parts through from_raw_parts. The input must be normalised exactly and then correctly rounded, or refused with InvalidArgument.
- R4 families take Spence's R4 text for these classes.

## 7. Huge-argument trig (item 7)
- `trig-huge-arg` has 52 rows: 48 from v2, plus 4 for `ExactNum::rem_pi` (FIXLIST:31 ZF-17, FIXLIST:127 ZF-M22).
- Exponents: 2^16-1, 2^16, 2^16+1 and EXPONENT_MAX.
- Required: correct bits (for rem_pi, x mod pi correctly rounded) or PrecisionRetryExhausted, identical on both targets.

## 8. R4 text and special rows (items 1, 5)
- **Ci and li (item 1):**
  - Ci(+0) = -Inf with kind DivisionByZero, in the zero row and the out-of-domain row, for ExactNum and ExactComplex.
  - li(1) = -Inf with DivisionByZero in the out-of-domain row.
  - The li zero row requires li(+0) = +0 exactly, with no error kind, because li(x)=Ei(ln x) and Ei(-Inf)=0 (DLMF 6.2.8). This differs from Spence's wording for that row, and Spence needs to confirm it.
  - "or documented" appears 0 times.
- **R4 text restored (item 5), 78 rows:**
  - as_raw_parts and significant_bits are exact and mode-independent.
  - random_seed, reseed_random and seeded_random are determined by the seed, digest-pinned, and byte-identical on i686 and x86_64.

## 9. NaN propagation (item 12)
Each nan row names the variant it inherits (A=InvalidArgument, B=DivisionByZero, N=no kind).
- **Real binary** rows nan-AB/BA/AN/NA follow Rule R: the first operand's kind wins, including the kind-less N. Pinned EXPECTED lines are cited by file:line. Ops not in EXPECTED (mul, div, ...) cite nan_error_kinds.rs:4-7 and CHANGELOG.md:131-132, and must be pinned at PM-1.
- **Real unary** rows nan-A/B/N carry the operand's kind.
- **ExactComplex** rows cite the pinned line when it exists (cadd/csub/asinh/cexp/csqrt/ln_gamma). Otherwise the row names re and im both carrying the re part's kind, or the first operand's kind for binary ops, and must be pinned at PM-1. A run that disagrees is FAIL.
- zb has 116 `nan` rows, all `Err(BoundError::Representation)`.
- Total: 2548 rows (zf 2432: 1332 unary and 1100 pair rows).

## 10. Oracles and PMCELL (R5)
- Each row pins one oracle per target (`oracle_x86_64`, `oracle_i686`). The i686 golds are cut from the x86_64 run.
- Each cell prints one line:
  `PMCELL crate=<zf|zb> family=<Type::name> class=<c> p=<p> rm=<m> target=<triple> oracle=<id> env=<path|rug-1.20.1/mpfr-4.2.2> verdict=<PASS|FAIL|NA_PROOF> ms=<t>`

## 11. ZF-PM5 spot-checks (item 3; closes old Q5)
- Each release candidate runs N=3 independent spot-checks per family, done by an agent other than the gold cutter. Any family that has ever given a wrong answer gets **N=5**. Any disagreement is FAIL.
- The per-family N and its citations are in `pm5-spotcheck-n.csv`: 212 families have N=5 (zf 202, zb 10), and 623 have N=3.
- The N=5 list is summarised below (`inputs/` files, file:line). `pm5-spotcheck-n.csv` is authoritative: it gives the exact family-by-family N and citations, and pm-coverage reads it:
  - sin: FIXLIST.csv:15 (ZF-01), FIXLIST.csv:127 (ZF-M22), fixes_after_1.0.12.txt:87
  - cos, sin_cos, tan: FIXLIST.csv:127, fixes:87, fixes:98
  - rem_pi: FIXLIST.csv:17 (ZF-03), FIXLIST.csv:31 (ZF-17), FIXLIST.csv:127
  - all ExactComplex::*: FIXLIST.csv:16 (ZF-02), SUMMARY.md:29 (#8), fixes:42-44
  - gamma: FIXLIST.csv:17, SUMMARY.md:24 (#3), fixes:31, 35, 98
  - ln_gamma: SUMMARY.md:24, fixes:38
  - digamma: SUMMARY.md:24, fixes:30, 34
  - exp: FIXLIST.csv:17, fixes:87, 98
  - gamma_pdf, beta_pdf, poisson_pmf, binomial_pmf, chi_squared_cdf, student_t_pdf: FIXLIST.csv:19 (ZF-05)
  - add: FIXLIST.csv:20 (ZF-06), FIXLIST.csv:117 (ZF-M19), SUMMARY.md:25 (#4), fixes:101
  - sub, mul, div: FIXLIST.csv:20, SUMMARY.md:25, fixes:101, fixes:44
  - fma: SUMMARY.md:25, fixes:45, 83
  - reciprocal: SUMMARY.md:25, fixes:19
  - pow, powi, powsi: SUMMARY.md:22 (#1), fixes:8, 21-25, 57-62
  - round: fixes:21-25
  - riemann_zeta: SUMMARY.md:23 (#2), fixes:37, 106
  - expm1: SUMMARY.md:26 (#5), fixes:33, 87
  - log1p: SUMMARY.md:26, fixes:40, 86, 87
  - log10: SUMMARY.md:26, fixes:52
  - log2: SUMMARY.md:26, fixes:15
  - log, ln: SUMMARY.md:26, fixes:18, 98
  - acosh: SUMMARY.md:26, fixes:49
  - atan2: SUMMARY.md:26, fixes:46
  - ldexp: SUMMARY.md:26, fixes:48
  - nth_root: SUMMARY.md:26, fixes:50-51
  - pi, e, ln_2, ln_10: SUMMARY.md:27 (#6), fixes:53
  - erf, erfc: SUMMARY.md:28 (#7), fixes:39, 41, 98
  - sinh, cosh, tanh: fixes:87, 98
  - sqrt: fixes:88, 99
  - cbrt: fixes:98
  - jacobi_am: fixes:89
  - bessel_j0, bessel_j1: fixes:36
  - ei, ci, li: fixes:98
  - elliptic_k, elliptic_e: fixes:98
  - clausen_cl2, clausen_cl3, weber_e, barnes_g: fixes:106
  - zb Ball add, sub, mul_exact, powi, sqrt, div, inv, hypot, abs: FIXLIST.csv:67 (ZB-01), FIXLIST.csv:68 (ZB-02)
- Excluded: FIXLIST.csv:18 (ZF-04, 33 "unspecified" functions) is an unproven accuracy gap, not a recorded wrong answer.

## 12. Tiers
- Fast tier: on every commit.
- Batch gate: fast + the cut heavy tier.
- Full heavy tier: once per RC, overnight, on **4 popslot workers (Paul ruling)**.
- `pm-max-i686` and `pm-max-x86_64` are separate heavy shards.

## 13. Budgets (tools/budget.pl on this CSV, PM-0 cost basis; Spence pmrev3 item 7)
Cost basis: zf 4.1 ms per cell; zb 12 ms per cell (this figure is unmeasured and blocks the RC budget until it is measured); threshold cells 210 per target (PM-0 corrected survey). The i686 factor comes from PM-0 and is passed as `PM0_FAC`: median 1.91, weighted 2.70 (the planning figure), max 5.39. The source is ~/harper-work/pm0-results-2026-10-30.md and its logs: ~/pop-logs/harper/2026-10-30/1638-pm0-factor2-i686.log and 1639-pm0-factor2-x86_64.log (factor), and 1704-pm0-budget2.log :31-35 (1.91), :41-45 (2.70), :51-55 (5.39). The mpfr peak is 724 MiB RSS with a wall time of 20:36 (1642-pm0-mpfr-rss2.log: RSS at :999, wall at :994), so 8 workers at 6G MemoryMax each have room. The old unsourced "4x" factor is gone. The figures below are from `PM0_FAC=<f> PM_WORKERS=8 perl tools/budget.pl precision-matrix-edges.csv`, logged in ~/pop-logs/harper/2026-10-30/1713-pm-rev3-budget-{1.91,2.70,5.39}.log.

| crate | target | fast | cut heavy | full heavy, 1 worker (1.91 / 2.70 / 5.39) |
|---|---|---|---|---|
| zf | x86_64 | 6.7 min | 0.70 h | 11.01 h |
| zf | i686 | 12.8 / 18.1 / 36.1 min | 1.33 / 1.88 / 3.76 h | 21.04 / 29.74 / 59.36 h |
| zb | x86_64 | 2.3 min | 0.27 h | 4.27 h |
| zb | i686 | 4.4 / 6.3 / 12.5 min | 0.51 / 0.72 / 1.44 h | 8.16 / 11.54 / 23.04 h |

- Full heavy, all four crate/target pairs: 44.49 / 56.57 / 97.69 slot-h, which is **5.56 / 7.07 / 12.21 h on 8 workers**.
- MAX shards: 352 zf + 35 zb rows per target, with a 10 min timeout. Under the MAX-timeout skip (section 18) the worst case is 64.5 slot-h, or 8.06 h on 8 workers. Without the skip it would be 129.0 slot-h, or 16.12 h.
- **8-worker RC wall time (full heavy + MAX worst case with skip): 13.62 / 15.13 / 20.27 h.** With no MAX timeouts it is 5.56 / 7.07 / 12.21 h. At the 2.70 planning factor the worst case needs two nights, so the run is resumable (section 19).
- Per-type evidence matching (PM-0 measured the factor per family, not per operand type) is a **PM-1 must**. Until it lands, no "covered" row counts are quoted from these figures.

## 14. Rulings
- R2: test_rng stays under the `random` feature (G13).
- R3: the 10 working-precision families are FIXLIST items.
- R8: Ball::div.
- Primality: the miller_rabin doc row FAILs until its doc is fixed.
- Paul: i686 MAX rule (section 4), ZF-PM5 N=3/N=5 (section 11), and 4 popslot workers (sections 12-13; the RC plan is now 8 workers, section 13). Spence-approved (pmrev3): li(+0)=+0; Ball and ComplexBall as zf types; the MAX-timeout skip (section 18).

## 15. FIXLIST additions
1. The 10 working-precision families.
2. Trig huge-argument reduction (ZF-01, ZF-M22).
3. The miller_rabin doc.
4. Every FAIL in the MAX shards (section 4).

## 16. Open questions
1. The zb cost of 12 ms per cell is unmeasured. Measuring it is a PM-1 must, together with per-type evidence matching (section 13).
2. The UNGENERATED golds are pinned at PM-3, and they block the RC until then (section 1).
3. The 33rd ZF-04 function (section 20) is unidentified.
4. Re-pin SRC_COMMIT when paul/zf-bound-merge merges (section 0).

## 17. pm-coverage (tools/pm_coverage.pl; Spence pmrev3 item 1)
pm_coverage.pl is a separate tool. lint.pl stays the structure and variant lint and must be clean, while pm-coverage is the RC gate. Usage: `perl tools/pm_coverage.pl [--mode rc|design-stage] [--results verdicts.tsv] precision-matrix-edges.csv pm-shards.csv [tree]`. pm-shards.csv is written by gen_edges.pl and gives the row count for each shard and family (fast, heavy, pm-max-x86_64, pm-max-i686, proof).
- a: an NA_PROOF row whose proof file has no `Second-agent sign-off: <agent> <date> <commit>` line (UNSIGNED) is a FAIL.
- b: an NA_PROOF row that cites a proof file that does not exist is a FAIL.
- c: an executed row with an UNGENERATED oracle pin is a FAIL.
- d: a `[PM1-UNPINNED]` row (an unpinned Rule R NaN pair) is a FAIL.
- e: a shard whose row count for a family is below the matrix count is a FAIL.
- f: given `--results`, any executed row whose verdict is not PASS is a FAIL, including SKIPPED_AFTER_TIMEOUT and missing rows. The RC requires a results file.
- Structural rules: S1, every prec-MAX-<target> row has the SKIPPED_AFTER_TIMEOUT verdict space; S2, every NA_PROOF row cites a pm-proofs path; S3, the executed verdict space is PASS|FAIL[|SKIPPED_AFTER_TIMEOUT].
- `--mode rc` (the default) enforces every rule and exits 1 on any count.
- `--mode design-stage` enforces only S1-S3 and e. It lists a-d as GAP lines with their full counts, and its output is marked NOT AN RC VERDICT.
- Each rule was proved red on its own tampered copy (log in the commit message), and the tamper copies were deleted.
- On these tables at design stage, rc is FAIL. The counts, by row, are: a=658 (74 files, all UNSIGNED), b=3357 (435 missing files out of 509 cited), c=19742, d=996, e=0, f not evaluated (no results yet). S1-S3 are 0. design-stage is a structural PASS that lists the same gaps.

## 18. MAX-timeout skip (Paul proposal, Spence pmrev3 item 6)
- Within a family, the x86_64 prec-MAX row runs first.
- If that row hits its 10 min timeout, the family's i686 prec-MAX row is not run. It is recorded as **SKIPPED_AFTER_TIMEOUT**.
- SKIPPED_AFTER_TIMEOUT counts as FAIL in the verdict space (`PASS|FAIL|SKIPPED_AFTER_TIMEOUT` on all 774 prec-MAX-<target> rows) and in pm-coverage rule f. It is never omitted.
- The family's RC stays blocked until its FIXLIST fix lands. After the fix, every skipped row reruns in full.
- The prec-MAXp1-refuse rows are instant refusals, so the skip does not apply to them.

## 19. Resumable RC run (described now, built at PM-1; Spence pmrev3 item 7)
- 8 popslot workers run overnight through poprun.
- Each row's verdict is appended (fsync) to verdicts.tsv as the row finishes. The record holds the row id, verdict, zf/zb SHA, env hash (the team-env manifest plus the target) and the gen_edges.pl hash.
- On restart, the runner skips any row already recorded with the same SHA, env hash and generator hash, and reruns the rest. A change to any of the three invalidates the record.
- If the run is not done by morning (the window is set at PM-1), it pauses: workers finish their current row and stop. It resumes the next night.
- **PM-1 selftest:** kill a run partway with SIGKILL, resume it, and prove that the final verdicts.tsv covers every matrix row exactly once. The check is no missing ids, no duplicate ids, and the count equals the matrix count, verified with pm-coverage rule f. The selftest is a PM-1 exit criterion.

## 20. Round-3 row changes (Spence pmrev3 items 2-5)
- **Item 2, ZF-04 at N=5:** FIXLIST.csv:18 (ZF-04, S1 wrong-answer) claims 33 "Accuracy: unspecified" functions: 29 in catalog.rs:575-1581 and 4 in heun_fox.rs. At 7955836 the source has 25 in catalog.rs plus the 4 in heun_fox.rs, which is 29. The 3 jacobi_arcsn/arccn/arcdn functions in catalog.rs:1255-1301 have no accuracy doc and are included as well. That makes 32 families by name, at N=5 for every type that has them, each cited in pm5-spotcheck-n.csv as FIXLIST.csv:18 plus the source lines. riemann_zeta keeps its recorded wrong answer (zf-audit-SUMMARY.md:23) alongside the ZF-04 citation. **One function of the 33 could not be identified.** It is an open question for Paul/Spence.
- **Item 3:** the 996 unpinned Rule R NaN-pair rows have a required result of `[PM1-UNPINNED] ... FAIL until pinned in tests/nan_error_kinds.rs EXPECTED at PM-1`. A run whose kind matches the guess still FAILs until then (pm-coverage rule d).
- **Item 4:** the signed-zero rows for li, ci, si, nth_root and arg carry exact expected values under C11 Annex G conjugate symmetry and the principal branch (for example ln(-1-0i) = -i*pi). Any mismatch is a red row. There are 10 such rows: ExactComplex arg, nth_root, si, ci and li; ExactNum nth_root, si, ci and li; and ExactNumArray nth_root. The "39 UNCLEAR" in Spence's log is the total of all UNCLEAR rows. The other 29 are R4 rows for decompositions, quadrature, ode and roots, and are unchanged. Paul adds the LIBRARY.md:399 doc row to FIXLIST.
- **Item 5:** see section 0. SRC_COMMIT 7955836 is on paul/zf-bound-merge, not on zf main (dc6bc87).
