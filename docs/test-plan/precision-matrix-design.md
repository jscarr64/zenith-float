# zf/zb precision matrix: design v7 (Jeff ruling 2026-10-30 6:00 PM CT, ZF-M27; v6 was Spence review 5)

Agent: Harper. Date: 2026-10-30. Edge table: `precision-matrix-edges.csv` (24654 rows, 20639 executed, 17 columns; the first 9 lines start with `#` and carry the pins below; the CSV header is the first line without `#`).

## 0. Pinned inputs (item 11)
The generator `tools/gen_edges.pl` reads only the in-branch copies below. It checks every sha256 before reading and stops with `PIN MISMATCH` if one differs. It also stops if `zenith-float-num/{src,tests}` or `zenith-bound/{src,tests}` differ from the source commit. The originals in ~/paul-work and ~/spence-work were copied unedited.

| input (docs/test-plan/inputs/) | copied from | sha256 |
|---|---|---|
| zf-inventory.csv | ~/paul-work/zf-audit-2026-10-30/inventory.csv | ba79198eb3486944e3b43a898389b81a7ce8fb69e4f67450e49841e8fdac1d59 |
| zb-inventory.csv | ~/paul-work/zb-audit-2026-10-30/inventory.csv (section 1 used) | b65c2bad80431500ab6280eb6ca545194407a35e9ae2f9719b5a6e5a38bdcd65 |
| r4-required-results.csv | Spence R4 (via ~/harper-work/pm-review-inputs) | 2746dc4d05cca03d07db6cd8db600e86369a7f33b450c6b679b1d0c1bf0e8c6d |
| spence-gen_r4.py | ~/spence-work/pm-review/gen_r4.py (G1-G3 lists, lines 6-19) | aa4d614617d086a49a2d625adbef1b516a7cac7ea37de58e391d40df2b51243f |
| FIXLIST.csv | ~/paul-work/zfzb-plan-2026-10-30/FIXLIST.csv | 3a8907f41341130e75e8dd2e0871c97e40cee96d8cf63103f7225e5d69f2606b (re-pinned pmrev7: adds ZF-M26, ZF-D02, ZF-D03 and ZF-M27 at :145 after :141; ZF-04=29, ZF-WP01..10 at :131-140, ZF-D01 at :141 unchanged) |
| fixes_after_1.0.12.txt | ~/paul-work/zf-audit-2026-10-30/fixes_after_1.0.12.txt | 03085a0acb3d5d0e51d8e8df13b3c4ad30d7843152d3d32024d6899c9f530c82 |
| zf-audit-SUMMARY.md | ~/paul-work/zf-audit-2026-10-30/SUMMARY.md | a7af866c6b2c6de32f4ed243170a4aa92fff19bc7e3978e86642a3dff1cb75f0 (re-pinned pmrev5) |

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
There are 7 modes. ToOdd: the oracle computes RNDN at p+1 and flips exact ties itself. None = RNDZ of the exact value at rp(p,64). **R6 (extended pmrev7, Jeff ruling 2026-10-30 6:00 PM CT, FIXLIST:145 ZF-M27): all 7 modes round at rp(p,64) and are bit-identical on i686 and x86_64, so each row has one gold for both targets.** Until ZF-M27 lands, cells in the 6 directed modes (every mode except ToEven; FIXLIST:145 says "all 6 directed modes") at a p that is not a multiple of 64 are red (FAIL); each such row carries a ZF-M27 pin note (16373 rows).

## 4. Precision axis and MAX_PRECISION (item 2; closes old Q1)
- The base axis is p ∈ {32,64,128,256,512,1024}, with 96 heavy points, and p=1,2,3 run in all 7 modes. p=0 and p=MAX+1 are refused with InvalidArgument.
- pmrev7 (ZF-M27, one row for both targets): the old `prec-MAX-x86_64`/`prec-MAX-i686` split is gone. MAX_PRECISION is target-dependent at 7955836 (defs.rs:61-90: min(2^60, usize::MAX/4)-1, so 2^60-1 on 64-bit and 2^30-1 on 32-bit targets). Each family now has two shared rows, run on both targets:
  - `prec-2^30-1` (387 rows): p=2^30-1 is MAX_PRECISION on i686 and an ordinary precision on x86_64; Paul ruled 2^30 bits is allocatable, so the required result below is one gold for both targets.
  - `prec-2^60-1` (387 rows): **FLAG, needs a decision from Jeff.** On x86_64 this is MAX_PRECISION and 2^57 bytes cannot be allocated (MemoryAllocation; Resource for zb); on i686 it is above MAX_PRECISION (InvalidArgument). One input cannot give the same result on both targets, so these rows are [PM1-UNPINNED] (rule d FAIL) until the precision limit is the same on both targets. Not guessed.
- The MAX rows run **only in the heavy tier, in their own shard `pm-max`, on both targets**. Each row has a 10-minute budget, and a timeout is FAIL.
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
- **Trig in zb (review 5 item 14):** zenith-bound has no trig families. At 7955836, `git grep -n -E 'fn (sin|cos|tan|sin_cos|rem_pi|sec|csc|cot)\b' -- zenith-bound/src` has 0 hits, `git grep -i -E 'trig|sin\(|cos\('` in zenith-bound/src has 0 hits, and inputs/zb-inventory.csv has no trig entry. So no BoundError::Resource trig rows are added (zb trig rows: 0). The trig-huge rows marked zb in v2 were the zf `Ball::sin/cos` and `ComplexBall::sin/cos` (ball.rs:105,112,352,370) and carry zf kinds.

## 6. Thresholds and operand classes (item 6)
- The threshold classes (ovf-max-finite, ovf-max-plus-half-ulp, unf-min-normal, unf-min-subnormal, unf-half-min-subnormal-tie) are kept, with the T1-T12 table and the i64 extremes for scalb, ldexp and ilogb.
- Three operand classes are restored, each with 835 rows (one per family):
  - `huge-exp`: an operand exponent in [EXPONENT_MAX-64, EXPONENT_MAX]
  - `tiny-exp`: an operand exponent in [EXPONENT_MIN, EXPONENT_MIN+64]
  - `subnormal-like`: unnormalised raw parts through from_raw_parts. The input must be normalised exactly and then correctly rounded, or refused with InvalidArgument.
- R4 families take Spence's R4 text for these classes.

## 7. Huge-argument trig (review 5 items 12-13)
- `trig-huge-arg` has 52 rows (13 families x 4 exponents), plus 26 `trig-pi-budget-*` rows (13 families x under/over).
- e = 2^16-1, 2^16, 2^16+1 (39 rows): pinned to the correctly rounded mpmath value at p in rm, identical on both targets. There is no cap here and no PrecisionRetryExhausted: a refusal is FAIL. Spence's timing: up to 0.07 s per row at p=1025, about 0.002 s at e=2^16+-1.
- e = EXPONENT_MAX (13 rows): pinned to a NaN with Error::MemoryAllocation, refused before any allocation, because the reduction needs more pi bits than the pi-bit budget constant. Any returned value is FAIL.
- **Requirement (ZF-M22, FIXLIST:127):** zf must publish one documented up-front pi-bit budget constant, the same on x86_64 and i686 (Rusty sets the value under ZF-M22). Argument reduction compares the pi bits it needs against the constant before allocating, and refuses with MemoryAllocation when the need is over it. The old usize-saturated `trig_reduce_extra` cap (ops/util.rs:92-95) is not a pass criterion.
- `trig-pi-budget-under` (needs just under the constant: correctly rounded value) and `trig-pi-budget-over` (just over: MemoryAllocation refusal) are `[PM1-UNPINNED]` until the constant is published and the gold is cut.

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

- **Oracle primary rule (review 5 item 8):** every oracle cell starts `primary: MPFR` or `primary: mpmath`. python-flint (arb) appears only as a second check and never generates a gold. The 95 digest rows name mpmath as their generator. Counts: MPFR primary 1228, mpmath primary 19411, flint primary 0, names neither 0.
- Each row has one gold for both targets (pmrev7, ZF-M27): `oracle_x86_64` pins it and `oracle_i686` says "same gold as oracle_x86_64"; the i686 run must be bit-identical. lint v6 rejects any target-specific gold, pin or class (tools/pm_target.pl, shared with gen_edges.pl).
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
- ZF-04 (FIXLIST.csv:18, now 29 = 25 catalog.rs + 4 heun_fox.rs) was excluded in v3. Since v4 every ZF-04 function is at N=5 (section 20).
- ZF-WP01..WP10 (FIXLIST.csv:131-140): the 10 working-precision ExactNum families (hermite_he, hermite_h, laguerre, gen_laguerre, chebyshev_t, chebyshev_u, gegenbauer, jacobi_p, legendre_p, assoc_legendre_p) move from N=3 to N=5, each citing its own ZF-WP row (section 21).

## 12. Tiers
- Fast tier: on every commit.
- Batch gate: fast + the cut heavy tier.
- Full heavy tier: once per RC, overnight, on **4 popslot workers (Paul ruling)**.
- `pm-max` is one heavy shard run on both targets (pmrev7; was `pm-max-i686` and `pm-max-x86_64`).

## 13. Budgets (tools/budget.pl on this CSV, PM-0 cost basis; Spence pmrev3 item 7)
Cost basis: zf 4.1 ms per cell; zb 12 ms per cell (this figure is unmeasured and blocks the RC budget until it is measured); threshold cells 210 per target (PM-0 corrected survey). The i686 factor comes from PM-0 and is passed as `PM0_FAC`: median 1.91, weighted 2.70 (the planning figure), max 5.39. The source is ~/harper-work/pm0-results-2026-10-30.md and its logs: ~/pop-logs/harper/2026-10-30/1638-pm0-factor2-i686.log and 1639-pm0-factor2-x86_64.log (factor), and 1704-pm0-budget2.log :31-35 (1.91), :41-45 (2.70), :51-55 (5.39). The mpfr peak is 724 MiB RSS with a wall time of 20:36 (1642-pm0-mpfr-rss2.log: RSS at :999, wall at :994), so 8 workers at 6G MemoryMax each have room. The old unsourced "4x" factor is gone. The figures below are from `PM0_FAC=<f> PM_WORKERS=8 perl tools/budget.pl precision-matrix-edges.csv`, logged in ~/pop-logs/harper/2026-10-30/1713-pm-rev3-budget-{1.91,2.70,5.39}.log.

| crate | target | fast | cut heavy | full heavy, 1 worker (1.91 / 2.70 / 5.39) |
|---|---|---|---|---|
| zf | x86_64 | 6.7 min | 0.70 h | 11.01 h |
| zf | i686 | 12.8 / 18.1 / 36.1 min | 1.33 / 1.88 / 3.76 h | 21.04 / 29.74 / 59.36 h |
| zb | x86_64 | 2.3 min | 0.27 h | 4.27 h |
| zb | i686 | 4.4 / 6.3 / 12.5 min | 0.51 / 0.72 / 1.44 h | 8.16 / 11.54 / 23.04 h |

- Full heavy, all four crate/target pairs: 44.49 / 56.57 / 97.69 slot-h, which is **5.56 / 7.07 / 12.21 h on 8 workers**.
- MAX shards (pmrev7: 774 shared rows run on each target, so the per-target run count is now 704 zf + 70 zb, not 352 + 35; the budget figures below are the v6 ones and are re-run by tools/budget.pl at PM-0): 352 zf + 35 zb rows per target, with a 10 min timeout. Under the MAX-timeout skip (section 18) the worst case is 64.5 slot-h, or 8.06 h on 8 workers. Without the skip it would be 129.0 slot-h, or 16.12 h.
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
5. API row (review 5 items 2-3): roots::brent and roots::illinois return Option<ExactNum>, so the no-sign-change refusal is pinned as a NaN with Error::InvalidArgument held in the returned Option. rk4 and rk45_adaptive return Option, so their pinned Error kinds FAIL until they can report the kind. Each needs a FIXLIST API row.
6. The pi-bit budget constant (section 7, ZF-M22).

## 16. Open questions
1. The zb cost of 12 ms per cell is unmeasured. Measuring it is a PM-1 must, together with per-type evidence matching (section 13).
2. The UNGENERATED golds are pinned at PM-3, and they block the RC until then (section 1).
3. Closed (pmrev5): Paul corrected ZF-04 to 29 (25 + 4) in FIXLIST.csv:18; re-pinned in section 0.
4. Re-pin SRC_COMMIT when paul/zf-bound-merge merges (section 0).
5. **PM-1 item (review 5 item 11):** the 8 degenerate golds (the 4 `degenerate` rows for lu, qr, svd and eigen, on x86_64 and i686) need a constructed exact input per row: an exactly singular matrix for lu, a zero R_ii for qr, repeated singular values for svd, and repeated eigenvalues for eigen. They are cut at PM-1, and pm-coverage rule c stays red on them until then.
6. Exact poles away from 0 (for example gamma at a negative integer) are excluded from the out-of-domain class and have no rows yet. That is a PM-1 item.

## 17. pm-coverage (tools/pm_coverage.pl; Spence pmrev3 item 1)
pm_coverage.pl is a separate tool. lint.pl stays the structure and variant lint and must be clean, while pm-coverage is the RC gate. Usage: `perl tools/pm_coverage.pl [--mode rc|design-stage] [--results verdicts.tsv --manifest run.manifest] precision-matrix-edges.csv pm-shards.csv [tree]`. pm-shards.csv is written by gen_edges.pl and gives the row count for each shard and family (fast, heavy, pm-max-x86_64, pm-max-i686, proof).
- a: an NA_PROOF row whose proof file has no `Second-agent sign-off: <agent> <date> <commit>` line (UNSIGNED) is a FAIL.
- b: an NA_PROOF row that cites a proof file that does not exist is a FAIL.
- c: an executed row with an UNGENERATED oracle pin is a FAIL.
- d: a `[PM1-UNPINNED]` row is a FAIL. It covers the unpinned Rule R NaN pairs and, since v5, every executed row whose single outcome cannot be pinned yet (section 21).
- e: a shard whose row count for a family is below the matrix count is a FAIL.
- f: given `--results` and `--manifest`, rule f parses every full record (section 17.1) and FAILs on each of: a duplicate row_id (the first record is kept and the duplicate is a FAIL; a later line never overwrites an earlier one, so SKIPPED_AFTER_TIMEOUT then PASS is red), a row_id not in the matrix, a malformed record (fewer than 7 fields), a mismatch in any of the 5 provenance fields against the manifest, an executed row with no record, and any verdict other than PASS (SKIPPED_AFTER_TIMEOUT and FAIL count as FAIL). The RC requires both files. `--results` without `--manifest` is a hard error. Since v3 (review 5 item 9), rule f is enforced in every mode once `--results` is given, so any FAIL record gives a non-zero exit in design-stage too. A duplicate key in the manifest is a hard error (dies).
- Structural rules: S1, every prec-MAX-<target> row has the SKIPPED_AFTER_TIMEOUT verdict space; S2, every NA_PROOF row cites a pm-proofs path; S3, the executed verdict space is PASS|FAIL[|SKIPPED_AFTER_TIMEOUT].
- `--mode rc` (the default) enforces every rule and exits 1 on any count.
- `--mode design-stage` enforces S1-S3 and e, plus f whenever `--results` is given. It lists a-d as GAP lines with their full counts, and its output is marked NOT AN RC VERDICT.
- Each rule was proved red on its own tampered copy (log in the commit message), and the tamper copies were deleted.
- On these tables at design stage (v7), rc is FAIL (exit 1). The counts, by row, are: a=658 (74 files, all UNSIGNED), b=3357 (435 missing files out of 509 cited), c=16022, d=1971, e=0, f not evaluated (no results yet). S1-S3 are 0. design-stage is a structural PASS (exit 0) that lists the same gaps.

### 17.1 Results record and run manifest (pmrev5 item 1)
- verdicts.tsv: UTF-8, one record per line, TAB-separated, `#` lines are comments. Fields in order: `row_id` (`crate|family|class`, exactly the matrix id), `verdict` (PASS, FAIL or SKIPPED_AFTER_TIMEOUT), `zf_sha` (40-hex zenith-float commit), `zb_sha` (40-hex zenith-bound commit), `env_name` (for example oracle-2026-10-30), `env_hash` (sha256 of the team-env manifest plus the target triple), `gen_sha256` (sha256 of tools/gen_edges.pl). Optional trailing fields: `target`, `secs`. Rule f ignores them.
- **Where the expected values come from:** a manifest file passed with `--manifest`, never the results file and never CLI value flags. It has `key=value` lines for exactly `zf_sha`, `zb_sha`, `env_name`, `env_hash`, `gen_sha256`; `#` lines and blank lines are ignored. A missing key or malformed line is a hard error. The PM-1 runner writes the manifest once at run start from the checked-out SHAs, the team-env manifest and the generator hash, and every record is compared to it field by field.
- Tamper proof (pmrev5): a clean control was green (rc 0). Each of dup, unknown id, zf_sha, zb_sha, env_name, env_hash and gen_sha256 was red on its own (rc 1). Logs are listed in the commit message, and the tamper copies were deleted.

## 18. MAX-timeout skip (Paul proposal, Spence pmrev3 item 6)
- Within a family, each MAX row (`prec-2^30-1`, `prec-2^60-1`) runs on the first target first.
- If that run hits its 10 min timeout, the same row's run on the other target is not started. It is recorded as **SKIPPED_AFTER_TIMEOUT**.
- SKIPPED_AFTER_TIMEOUT counts as FAIL in the verdict space (`PASS|FAIL|SKIPPED_AFTER_TIMEOUT` on all 774 MAX rows) and in pm-coverage rule f. It is never omitted.
- The family's RC stays blocked until its FIXLIST fix lands. After the fix, every skipped row reruns in full.
- The prec-MAXp1-refuse rows are instant refusals, so the skip does not apply to them.

## 19. Resumable RC run (described now, built at PM-1; Spence pmrev3 item 7)
- 8 popslot workers run overnight through poprun.
- Each row's verdict is appended (fsync) to verdicts.tsv as the row finishes. The record holds the row id, verdict, zf/zb SHA, env hash (the team-env manifest plus the target) and the gen_edges.pl hash.
- On restart, the runner skips any row already recorded with the same SHA, env hash and generator hash, and reruns the rest. A change to any of the three invalidates the record. **Restart rule (review 5 item 10):** if the run manifest changes in any field, the runner starts a new verdicts.tsv and reruns every row; it never mixes records from two manifests.
- If the run is not done by morning (the window is set at PM-1), it pauses: workers finish their current row and stop. It resumes the next night.
- **PM-1 selftest:** kill a run partway with SIGKILL, resume it, and prove that the final verdicts.tsv covers every executed matrix id exactly once. The check is no missing ids, no duplicate ids, and the count equals the executed-id count (the rows that are not NA_PROOF; 20639 in v6), verified with the fixed pm-coverage rule f (section 17.1: full-record parse, duplicate and unknown ids FAIL, all 5 provenance fields checked against the run manifest). **PM-1 exit criterion:** after SIGKILL and resume, `pm_coverage.pl --mode rc --results verdicts.tsv --manifest run.manifest` shows rule f count=0 (no duplicate, unknown, mismatched or missing record) on the selftest subset, and a second selftest that resumes under a changed manifest field shows rule f red. PM-1 does not exit until both hold.

## 20. Round-3 row changes (Spence pmrev3 items 2-5)
- **Item 2, ZF-04 at N=5 (re-counted in round 4):** FIXLIST.csv:18 and ~/paul-work/zf-audit-2026-10-30/SUMMARY.md:30 claim 29 "Accuracy: unspecified" functions in catalog.rs:575-1581 plus 4 in heun_fox.rs (:576/:640/:771/:967), 33 in total. At SRC_COMMIT 7955836, `git grep -c "Accuracy: unspecified"` finds catalog.rs 25, heun_fox.rs 4 and accuracy.rs 1 (the Accuracy::Unspecified Display string at :42, not a function). catalog.rs:575-1582 has 28 pub fns with an Accuracy line: 25 say unspecified and 3 (jacobi_arcsn :1254, jacobi_arccn :1273, jacobi_arcdn :1300) say "Accuracy: correctly-rounded". Paul's own ~/paul-work/zf-audit-2026-10-30/inventory.csv also has exactly 25 catalog.rs rows with documented_accuracy unspecified. Neither FIXLIST nor SUMMARY names the 29, so 4 cannot be matched; none are invented. The likely explanation, unconfirmed: 25 + 4 = 29 was the true total, and the 4 heun_fox functions were then added again. N=5 covers 25 + 4 + 3 jacobi_arc = **32** families (not 36), each cited in pm5-spotcheck-n.csv with FIXLIST.csv:18 and its own doc and fn source lines. riemann_zeta keeps its recorded wrong answer (zf-audit-SUMMARY.md:23).
- **Item 3:** the 996 unpinned Rule R NaN-pair rows have a required result of `[PM1-UNPINNED] ... FAIL until pinned in tests/nan_error_kinds.rs EXPECTED at PM-1`. A run whose kind matches the guess still FAILs until then (pm-coverage rule d).
- **Item 4:** the signed-zero rows for li, ci, si, nth_root and arg carry exact expected values under C11 Annex G conjugate symmetry and the principal branch (for example ln(-1-0i) = -i*pi). Any mismatch is a red row. There are 10 such rows: ExactComplex arg, nth_root, si, ci and li; ExactNum nth_root, si, ci and li; and ExactNumArray nth_root. The "39 UNCLEAR" in Spence's log is the total of all UNCLEAR rows. The other 29 are R4 rows for decompositions, quadrature, ode and roots, and were pinned in v5 (section 21). The LIBRARY.md:399-410 doc row is now ZF-D01 (FIXLIST.csv:141), and all 10 signed-zero rows cite it.
- **Item 5:** see section 0. SRC_COMMIT 7955836 is on paul/zf-bound-merge, not on zf main (dc6bc87).

## 21. Round-5 changes (Spence fourth review items 3-5)
- **No UNCLEAR, no either/or gold.** Every executed row pins exactly one outcome: one value, or one named refusal. The other branch is a FAIL. lint.pl v4 rejects `UNCLEAR` anywhere, and `or`, `either` or a standalone `|` in required_result or refusal_variant, and exits 1 on any BAD line. Each rule was tamper-tested red against a clean green control.
- Class rulings (Spence): gauss_kronrod and roots illinois/brent are correctly rounded at p (Ziv), with PrecisionRetryExhausted only where the gold pins it; brent/illinois with no sign change give Error::InvalidArgument (new `no-sign-change` rows). rk4 is the exact fixed-step recurrence, correctly rounded per step; InvalidArgument only for a bad step or interval (new `bad-step-interval` row). rk45_adaptive is within tol of mpmath odefun at 2p bits in env oracle-2026-10-30, else PrecisionRetryExhausted; it is pinned per row, and rows not yet pinned are `[PM1-UNPINNED]`. LU/QR/SVD/eigen: each factor is correctly rounded under the pinned normalisation; None only for an exact pivot tie past the lowest-index rule, repeated singular values or eigenvalues, or a zero R_ii (new `degenerate` rows). from_raw_parts is exactly normalized when the parts are valid at p, else InvalidArgument (new `raw-parts-invalid` row). Where a row once said "MPFR or mpmath", the pin now names the one oracle from that row's oracle column.
- MAX rows: x86_64 is pinned to the MemoryAllocation (zb Resource) refusal. i686 is pinned to the correct result or enclosure (Paul: the size is allocatable), and a refusal there is a FAIL.
- Rows that cannot have one outcome yet are `[PM1-UNPINNED]` (rule d): 61 new rows, namely trig-huge-arg e=2^16-1 and e=2^16 (the cap rule is undecided), exponent extremes on helper functions, zero inputs whose doc leaves the result open (for example the zero polynomial and "DivisionByZero or None per doc"), and rk45.
- Row deltas vs 3886199: 24620 -> 24628 (+8 new rows, 0 removed); 3707 rows changed (col10 3637, col11 2157), of which 3646 now carry a single pin and 61 are `[PM1-UNPINNED]`. UNCLEAR 29 -> 0; either/or in cols 10-11 3698 -> 0; `[PM1-UNPINNED]` 996 -> 1057; UNGENERATED pins 19742 -> 19750. The 10 ZF-WP families move N=3 -> N=5; the 4 heun_fox ZF-04 citations now read 29.

## 22. Round-6 changes (Spence review 5 items 1-14)
- **Hedge rule (item 1):** lint.pl v5 and gen_edges.pl share tools/pm_hedge.pl. It rejects in required_result and refusal_variant: TBD, TODO, UNCLEAR, either/or, `, else` and else, unless, otherwise, and `; ` or more than one value in refusal_variant. The `/` rule is scoped to cols 10-11: a `/` is allowed only in a source path, inside parentheses, next to a digit, or in math with `=`, `*` or `^` (for example `4/5`, `cos(pi/n)`, `x=a/b`, a source path `src/x.rs:9`), and any other `/` between alternatives (`a/b`, `Up(+)/FromZero`) is a hedge. Every executed row has exactly one outcome or `[PM1-UNPINNED]`. Tamper proof (pmrev6): against a clean green control (lint rc 0), each of TBD, TODO, slash, `, else`, unless, `; ` in refusal_variant and two values in refusal_variant was red on its own (rc 1, one BAD line each). pm-coverage: a FAIL record is red in rc and design-stage (rc 1), and a duplicate manifest key dies (rc 25) in both modes, each against a green control. The tamper copies were deleted.
- Rulings applied: brent/illinois no-sign-change and rk4 bad-step-interval pin Error::InvalidArgument (section 15 item 5). rk45: all 33 rows cite mpmath odefun at 2p bits, with PrecisionRetryExhausted where tol cannot be met. from_raw_parts: all 20 rows pin one outcome, exact normalised or InvalidArgument. Decompositions: per-class pins on the 132 rows; zero, singular and degenerate are None; no `else ExponentOverflow`. lu: ties go to the lowest index, and None only for an exactly singular matrix. Trig: section 7.
- Row deltas vs 3c0e1e5: 24628 -> 24654 (+26 trig-pi-budget rows); executed 20613 -> 20639. Hedges, executed rows: else 1661 -> 0, slash 5254 -> 0, rv-multi 1006 -> 0. `[PM1-UNPINNED]` executed 1057 -> 1584 (unbounded huge-exp helper families now carry a pin rule). UNGENERATED pins 19750 -> 19776.

## 23. Round-7 changes (Jeff ruling 2026-10-30 6:00 PM CT, FIXLIST:145 ZF-M27)
- **One gold per row:** target-specific rows 18523 before (class 774, gold 18453, pin 1381), 0 after. oracle_i686 now reads "same gold as oracle_x86_64" on every executed row. gen_edges.pl dies on any target-specific row, and lint v6 is BAD on one (tamper: control clean; target gold, target pin and target class each red alone, rc 1).
- **MAX:** 387+387 per-target rows became 387 `prec-2^30-1` shared rows (correct result, one gold) and 387 `prec-2^60-1` shared rows flagged [PM1-UNPINNED] (section 4). Rule d moves 1584 to 1971 (+387).
- **M27 red:** 0 rows carried a ZF-M27 pin note before, 16373 after (= every executed row in a directed mode other than ToEven; the 6 directed modes are read as every mode except ToEven).
- **FIXLIST re-pin:** sha 9fc4a7a5... to 3a8907f4...; existing line citations unchanged.
- **Row deltas:** rows 24654 to 24654, executed 20639 to 20639, NA_PROOF 4015 to 4015; rule c 19776 to 16022 (rows whose x86_64 gold is pinned no longer wait on a separate i686 gold); rule d 1584 to 1971.
