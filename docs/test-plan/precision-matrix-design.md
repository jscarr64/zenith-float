# zf/zb precision matrix: design v10 (Jeff steering 2026-10-30 6:32 PM CT: PRINCIPLES.md P6, drop i686; v9 was Jeff task 2026-10-30 6:24 PM CT: FIXLIST TD rows, tie rows; v8 was Jeff tasks 6:11/6:13/6:15 PM CT: ZF-M28, ZF-M29, ZF-M30, ZF-D04, PRINCIPLES.md P1-P5; v7 was the ZF-M27 ruling)

Agent: Harper. Date: 2026-10-30. Edge table: `precision-matrix-edges.csv` (24314 rows, 20299 executed, 4015 NA_PROOF, 17 columns; the first 9 lines start with `#` and carry the pins below; the CSV header is the first line without `#`). Governing principles: ~/pop-logs/PRINCIPLES.md (read-only; sha256 d39b9a605a02278116f032e65283da376d7d76e904049bffcfdb35c0267061f0 as read 2026-10-30 18:3x CT; was 81780e34... at v9), P1-P5 applied in section 24, P6 in section 26.

## 0. Pinned inputs (item 11)
The generator `tools/gen_edges.pl` reads only the in-branch copies below. It checks every sha256 before reading and stops with `PIN MISMATCH` if one differs. It also stops if `zenith-float-num/{src,tests}` or `zenith-bound/{src,tests}` differ from the source commit. The originals in ~/paul-work and ~/spence-work were copied unedited.

| input (docs/test-plan/inputs/) | copied from | sha256 |
|---|---|---|
| zf-inventory.csv | ~/paul-work/zf-audit-2026-10-30/inventory.csv | ba79198eb3486944e3b43a898389b81a7ce8fb69e4f67450e49841e8fdac1d59 |
| zb-inventory.csv | ~/paul-work/zb-audit-2026-10-30/inventory.csv (section 1 used) | b65c2bad80431500ab6280eb6ca545194407a35e9ae2f9719b5a6e5a38bdcd65 |
| r4-required-results.csv | Spence R4 (via ~/harper-work/pm-review-inputs) | 2746dc4d05cca03d07db6cd8db600e86369a7f33b450c6b679b1d0c1bf0e8c6d |
| spence-gen_r4.py | ~/spence-work/pm-review/gen_r4.py (G1-G3 lists, lines 6-19) | aa4d614617d086a49a2d625adbef1b516a7cac7ea37de58e391d40df2b51243f |
| FIXLIST.csv | ~/paul-work/zfzb-plan-2026-10-30/FIXLIST.csv (first 196 lines) | 56633946854c27c8c702b24157119102fed1e2e7a8044eb18a7111411bce08af (re-pinned pmrev10, 196 lines: adds ZF-M31 :195 and ZF-M32 :196; ZF-M27 :145 .. ZF-D04 :149 and ZF-TD01..ZF-TD45 :150-194 unchanged. DRIFT: at 18:33:43 CT the original grew to 197 lines, sha 02c209a5..., adding ZB-D05 at :197; this pin is its 196-line prefix, per the 6:32 PM CT instruction) |
| tie-resolution.tsv | Harper pmrev9 (one outcome per fixed-mode tie family, read from the source at 7955836; section 25) | 3ee454aefee8ba2521ae63fa94f721f7b1ea219486f93a8583905b54cab6c5fb |
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
There are 7 modes. ToOdd: the oracle computes RNDN at p+1 and flips exact ties itself. None = RNDZ of the exact value at rp(p,64). **R6 (pmrev8, Jeff ruling 2026-10-30 6:00 PM CT, FIXLIST:145 ZF-M27): all 7 modes round at rp(p,64); under P6 every supported target has a 64-bit word, so ZF-M27 is met on every target and each row has one gold.** v9 had 16725 rows red "until ZF-M27" (directed modes at p not a multiple of 64, and zf `prec-2^30-1`); v10 has 0 (section 26).

## 4. Precision axis and MAX_PRECISION (item 2; closes old Q1)
- The base axis is p ∈ {32,64,128,256,512,1024}, with 96 heavy points, and p=1,2,3 run in all 7 modes. p=0 and p=MAX+1 are refused with InvalidArgument.
- pmrev10 (P6; ZF-M32, FIXLIST:196, supersedes ZF-M28 :146): MAX_PRECISION is re-derived for 64-bit targets by Rusty. Today it is min(2^60, usize::MAX/4)-1 (zenith-float-num/src/defs.rs:77) = 2^60-1 on every 64-bit target. Each family keeps two shared rows, one gold (P1):
  - `prec-2^60-1` (387 rows): `[PM1-UNPINNED]`, red until Rusty pins the ZF-M32 constant; the refusal variant is pinned per row in the gold.
  - `prec-MAXp1-refuse` (387 rows): MAX+1 refuses with zf `Error::InvalidArgument` and zb `BoundError::InvalidArgument` (zenith-bound/src/error.rs:22). The v8/v9 `prec-2^30-1` and `prec-2^30` rows (387 + 387) are removed: 2^30-1 was the i686 limit (section 26).
- The MAX rows (`prec-2^60-1`, `m30-p2^27`) run **only in the heavy tier, in the shard `pm-max`, on every 64-bit target with the same gold**. The 10-minute per-run budget is execution config (section 24), not an expected outcome; a timeout is FAIL.
- The required result is the correct value, or the refusal pinned on the row. A wrong value, a different variant, a NaN-as-error, a panic or a timeout is FAIL (P2), and **every FAIL becomes a FIXLIST row**.

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
- e = 2^16-1, 2^16, 2^16+1 (39 rows): pinned to the correctly rounded mpmath value at p in rm, identical on every 64-bit target. There is no cap here and no PrecisionRetryExhausted: a refusal is FAIL. Spence's timing: up to 0.07 s per row at p=1025, about 0.002 s at e=2^16+-1.
- e = EXPONENT_MAX (13 rows): pinned to a NaN with Error::MemoryAllocation, refused before any allocation, because the reduction needs more pi bits than the pi-bit budget constant. Any returned value is FAIL.
- **Requirement (ZF-M22, FIXLIST:127):** zf must publish one documented up-front pi-bit budget constant, the same on every target (Rusty sets the value under ZF-M22). Argument reduction compares the pi bits it needs against the constant before allocating, and refuses with MemoryAllocation when the need is over it. The old usize-saturated `trig_reduce_extra` cap (ops/util.rs:92-95) is not a pass criterion.
- `trig-pi-budget-under` (needs just under the constant: correctly rounded value) and `trig-pi-budget-over` (just over: MemoryAllocation refusal) are `[PM1-UNPINNED]` until the constant is published and the gold is cut.

## 8. R4 text and special rows (items 1, 5)
- **Ci and li (item 1):**
  - Ci(+0) = -Inf with kind DivisionByZero, in the zero row and the out-of-domain row, for ExactNum and ExactComplex.
  - li(1) = -Inf with DivisionByZero in the out-of-domain row.
  - The li zero row requires li(+0) = +0 exactly, with no error kind, because li(x)=Ei(ln x) and Ei(-Inf)=0 (DLMF 6.2.8). This differs from Spence's wording for that row, and Spence needs to confirm it.
  - "or documented" appears 0 times.
- **R4 text restored (item 5), 78 rows:**
  - as_raw_parts and significant_bits are exact and mode-independent.
  - random_seed, reseed_random and seeded_random are determined by the seed, digest-pinned, and byte-identical on every 64-bit target.

## 9. NaN propagation (item 12)
Each nan row names the variant it inherits (A=InvalidArgument, B=DivisionByZero, N=no kind).
- **Real binary** rows nan-AB/BA/AN/NA follow Rule R: the first operand's kind wins, including the kind-less N. Pinned EXPECTED lines are cited by file:line. Ops not in EXPECTED (mul, div, ...) cite nan_error_kinds.rs:4-7 and CHANGELOG.md:131-132, and must be pinned at PM-1.
- **Real unary** rows nan-A/B/N carry the operand's kind.
- **ExactComplex** rows cite the pinned line when it exists (cadd/csub/asinh/cexp/csqrt/ln_gamma). Otherwise the row names re and im both carrying the re part's kind, or the first operand's kind for binary ops, and must be pinned at PM-1. A run that disagrees is FAIL.
- zb has 116 `nan` rows, all `Err(BoundError::Representation)`.
- Total: 2548 rows (zf 2432: 1332 unary and 1100 pair rows).

- **Oracle primary rule (review 5 item 8):** every oracle cell starts `primary: MPFR` or `primary: mpmath`. python-flint (arb) appears only as a second check and never generates a gold. The 95 digest rows name mpmath as their generator. Counts: MPFR primary 1228, mpmath primary 19411, flint primary 0, names neither 0.
- Each row has one gold for every 64-bit target (ZF-M27, P1, P6): `oracle_x86_64` pins it and `oracle_aarch64` (column renamed in v10; was oracle_i686) says "same gold as oracle_x86_64"; the aarch64 run must be bit-identical once TD40 installs a runner. lint v10 FAILs any i686, 32-bit, thumbv*, armv7* or target_pointer_width="32" reference in a row or the header (pm_p6). lint v7 rejects any target-specific gold, pin or class, and (P1, pmrev8) FAILs any row whose limit, refusal variant, budget or result is target-dependent (tools/pm_target.pl `pm_target_dependent`, shared with gen_edges.pl).
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
- `pm-max` is one heavy shard run on every 64-bit target (pmrev7; was split per target before).

## 13. Budgets (tools/budget.pl on this CSV, PM-0 cost basis; re-run pmrev10)
Cost basis: zf 4.1 ms per cell; zb 12 ms per cell (unmeasured; blocks the RC budget until measured); threshold cells 210 (PM-0 corrected survey). x86_64 only: the i686 factor (`PM0_FAC`) is removed under P6. No aarch64 factor: qemu-user is not installed on Pop (log 1833-inv9), so none was measured; one is added only when measured under qemu (TD40 file). The prec axis is the 49-point 64-bit set. Figures from `perl tools/budget.pl precision-matrix-edges.csv`, log ~/pop-logs/harper/2026-10-30/1846-budget10b.log (before: 1846-budget10-before.log).

| crate | target | fast | cut heavy | full heavy, 1 worker |
|---|---|---|---|---|
| zf | x86_64 | 6.7 min | 0.70 h | 5.65 h (v9 11.04 h) |
| zb | x86_64 | 2.3 min | 0.27 h | 2.18 h (v9 4.27 h) |

- Full heavy: 7.83 slot-h, **0.98 h on 8 workers** (v9 at factor 2.70: 56.66 slot-h, 7.08 h).
- MAX shard: 399 zf rows (352 `prec-2^60-1` + 47 `m30-p2^27`) and 35 zb rows, 10 min timeout. Worst case 72.3 slot-h, 9.04 h on 8 workers (v9: 1150 + 105 rows).
- **8-worker RC wall (full heavy + MAX worst case): 10.02 h** (v9 33.23 h); with no MAX timeouts 0.98 h.
- Per-type evidence matching (PM-0 measured the factor per family, not per operand type) is a **PM-1 must**. Until it lands, no "covered" row counts are quoted from these figures.

## 14. Rulings
- R2: test_rng stays under the `random` feature (G13).
- R3: the 10 working-precision families are FIXLIST items.
- R8: Ball::div.
- Primality: the miller_rabin doc row FAILs until its doc is fixed.
- Paul: the i686 MAX rule (section 4; withdrawn under P6, section 26), ZF-PM5 N=3/N=5 (section 11), and 4 popslot workers (sections 12-13; the RC plan is now 8 workers, section 13). Spence-approved (pmrev3): li(+0)=+0; Ball and ComplexBall as zf types; the MAX-timeout skip (section 18).

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
3. Removed (pmrev9, section 25): closed in pmrev5.
4. Re-pin SRC_COMMIT when paul/zf-bound-merge merges (section 0).
5. **PM-1 item (review 5 item 11):** the 4 degenerate golds (the 4 `degenerate` rows for lu, qr, svd and eigen, one gold per row for every 64-bit target) need a constructed exact input per row: an exactly singular matrix for lu, a zero R_ii for qr, repeated singular values for svd, and repeated eigenvalues for eigen. They are cut at PM-1, and pm-coverage rule c stays red on them until then.
6. **PM-1 item (pmrev8): poles away from 0.** Exact poles away from 0 (for example gamma at a negative integer) are excluded from the out-of-domain class and have no rows yet. PM-1 adds one row per pole-bearing family at a pole away from 0, with a P2 refusal or the documented infinity as the one gold.
7. Removed (pmrev9, section 25): the 309 fixed-mode tie rows each have one outcome now (P1, R3, R6, ZF-M27); 2 stay red with a named fact and owner.
No P1-P5 item is an open question; each is applied in section 24.

## 17. pm-coverage (tools/pm_coverage.pl; Spence pmrev3 item 1)
pm_coverage.pl is a separate tool. lint.pl stays the structure and variant lint and must be clean, while pm-coverage is the RC gate. Usage: `perl tools/pm_coverage.pl [--mode rc|design-stage] [--results verdicts.tsv --manifest run.manifest] precision-matrix-edges.csv pm-shards.csv [tree]`. pm-shards.csv is written by gen_edges.pl and gives the row count for each shard and family (fast, heavy, pm-max, proof).
- a: an NA_PROOF row whose proof file has no `Second-agent sign-off: <agent> <date> <commit>` line (UNSIGNED) is a FAIL.
- b: an NA_PROOF row that cites a proof file that does not exist is a FAIL.
- c: an executed row with an UNGENERATED oracle pin is a FAIL.
- d: a `[PM1-UNPINNED]` row is a FAIL. It covers the unpinned Rule R NaN pairs and, since v5, every executed row whose single outcome cannot be pinned yet (section 21).
- e: a shard whose row count for a family is below the matrix count is a FAIL.
- f: given `--results` and `--manifest`, rule f parses every full record (section 17.1) and FAILs on each of: a duplicate row_id (the first record is kept and the duplicate is a FAIL; a later line never overwrites an earlier one, so SKIPPED_AFTER_TIMEOUT then PASS is red), a row_id not in the matrix, a malformed record (fewer than 7 fields), a mismatch in any of the 5 provenance fields against the manifest, an executed row with no record, and any verdict other than PASS (SKIPPED_AFTER_TIMEOUT and FAIL count as FAIL). The RC requires both files. `--results` without `--manifest` is a hard error. Since v3 (review 5 item 9), rule f is enforced in every mode once `--results` is given, so any FAIL record gives a non-zero exit in design-stage too. A duplicate key in the manifest is a hard error (dies).
- Structural rules: S1, every MAX row (`prec-MAX-*`, `prec-2^60-1`, `m30-*`) has the SKIPPED_AFTER_TIMEOUT verdict space; S2, every NA_PROOF row cites a pm-proofs path; S3, the executed verdict space is PASS|FAIL[|SKIPPED_AFTER_TIMEOUT].
- `--mode rc` (the default) enforces every rule and exits 1 on any count.
- `--mode design-stage` enforces S1-S3 and e, plus f whenever `--results` is given. It lists a-d as GAP lines with their full counts, and its output is marked NOT AN RC VERDICT.
- Each rule was proved red on its own tampered copy (log in the commit message), and the tamper copies were deleted.
- On these tables at design stage (v8), rc is FAIL (exit 1; log 1821-cov8rc-final2.log). The counts, by row, are: a=658 (74 files, all UNSIGNED), b=3357 (435 missing files out of 509 cited), c=16346, d=1584, e=0, f not evaluated (no results yet). S1-S3 are 0. design-stage is a structural PASS (exit 0; 1821-cov8design-stage-final2.log) that lists the same gaps.

### 17.1 Results record and run manifest (pmrev5 item 1)
- verdicts.tsv: UTF-8, one record per line, TAB-separated, `#` lines are comments. Fields in order: `row_id` (`crate|family|class`, exactly the matrix id), `verdict` (PASS, FAIL or SKIPPED_AFTER_TIMEOUT), `zf_sha` (40-hex zenith-float commit), `zb_sha` (40-hex zenith-bound commit), `env_name` (for example oracle-2026-10-30), `env_hash` (sha256 of the team-env manifest plus the target triple), `gen_sha256` (sha256 of tools/gen_edges.pl). Optional trailing fields: `target`, `secs`. Rule f ignores them.
- **Where the expected values come from:** a manifest file passed with `--manifest`, never the results file and never CLI value flags. It has `key=value` lines for exactly `zf_sha`, `zb_sha`, `env_name`, `env_hash`, `gen_sha256`; `#` lines and blank lines are ignored. A missing key or malformed line is a hard error. The PM-1 runner writes the manifest once at run start from the checked-out SHAs, the team-env manifest and the generator hash, and every record is compared to it field by field.
- Tamper proof (pmrev5): a clean control was green (rc 0). Each of dup, unknown id, zf_sha, zb_sha, env_name, env_hash and gen_sha256 was red on its own (rc 1). Logs are listed in the commit message, and the tamper copies were deleted.

## 18. MAX-timeout skip (Paul proposal, Spence pmrev3 item 6)
- Within a family, each MAX row (`prec-2^60-1`, `m30-*`) runs on the first target first.
- If that run hits its 10 min timeout, the same row's run on the other target is not started. It is recorded as **SKIPPED_AFTER_TIMEOUT**.
- SKIPPED_AFTER_TIMEOUT counts as FAIL in the verdict space (`PASS|FAIL|SKIPPED_AFTER_TIMEOUT` on all 1255 MAX rows) and in pm-coverage rule f. It is never omitted.
- The family's RC stays blocked until its FIXLIST fix lands. After the fix, every skipped row reruns in full.
- The prec-MAXp1-refuse rows are instant refusals, so the skip does not apply to them.

## 19. Resumable RC run (described now, built at PM-1; Spence pmrev3 item 7)
- 8 popslot workers run overnight through poprun.
- Each row's verdict is appended (fsync) to verdicts.tsv as the row finishes. The record holds the row id, verdict, zf/zb SHA, env hash (the team-env manifest plus the target) and the gen_edges.pl hash.
- On restart, the runner skips any row already recorded with the same SHA, env hash and generator hash, and reruns the rest. A change to any of the three invalidates the record. **Restart rule (review 5 item 10):** if the run manifest changes in any field, the runner starts a new verdicts.tsv and reruns every row; it never mixes records from two manifests.
- If the run is not done by morning (the window is set at PM-1), it pauses: workers finish their current row and stop. It resumes the next night.
- **PM-1 selftest:** kill a run partway with SIGKILL, resume it, and prove that the final verdicts.tsv covers every executed matrix id exactly once. The check is no missing ids, no duplicate ids, and the count equals the executed-id count (the rows that are not NA_PROOF; 20639 in v6), verified with the fixed pm-coverage rule f (section 17.1: full-record parse, duplicate and unknown ids FAIL, all 5 provenance fields checked against the run manifest). **PM-1 exit criterion:** after SIGKILL and resume, `pm_coverage.pl --mode rc --results verdicts.tsv --manifest run.manifest` shows rule f count=0 (no duplicate, unknown, mismatched or missing record) on the selftest subset, and a second selftest that resumes under a changed manifest field shows rule f red. PM-1 does not exit until both hold.

$&
Sections 20-25 are the decision log as written at each round. Where they name i686, 32-bit, `prec-2^30*`, `m30-p2^24` or "both targets", section 26 (P6) supersedes them.
- **Item 2, ZF-04 at N=5 (re-counted in round 4):** FIXLIST.csv:18 and ~/paul-work/zf-audit-2026-10-30/SUMMARY.md:30 claim 29 "Accuracy: unspecified" functions in catalog.rs:575-1581 plus 4 in heun_fox.rs (:576/:640/:771/:967), 33 in total. At SRC_COMMIT 7955836, `git grep -c "Accuracy: unspecified"` finds catalog.rs 25, heun_fox.rs 4 and accuracy.rs 1 (the Accuracy::Unspecified Display string at :42, not a function). catalog.rs:575-1582 has 28 pub fns with an Accuracy line: 25 say unspecified and 3 (jacobi_arcsn :1254, jacobi_arccn :1273, jacobi_arcdn :1300) say "Accuracy: correctly-rounded". Paul's own ~/paul-work/zf-audit-2026-10-30/inventory.csv also has exactly 25 catalog.rs rows with documented_accuracy unspecified. Neither FIXLIST nor SUMMARY names the 29, so 4 cannot be matched; none are invented. The likely explanation, unconfirmed: 25 + 4 = 29 was the true total, and the 4 heun_fox functions were then added again. N=5 covers 25 + 4 + 3 jacobi_arc = **32** families (not 36), each cited in pm5-spotcheck-n.csv with FIXLIST.csv:18 and its own doc and fn source lines. riemann_zeta keeps its recorded wrong answer (zf-audit-SUMMARY.md:23).
- **Item 3:** the 996 unpinned Rule R NaN-pair rows have a required result of `[PM1-UNPINNED] ... FAIL until pinned in tests/nan_error_kinds.rs EXPECTED at PM-1`. A run whose kind matches the guess still FAILs until then (pm-coverage rule d).
- **Item 4:** the signed-zero rows for li, ci, si, nth_root and arg carry exact expected values under C11 Annex G conjugate symmetry and the principal branch (for example ln(-1-0i) = -i*pi). Any mismatch is a red row. There are 10 such rows: ExactComplex arg, nth_root, si, ci and li; ExactNum nth_root, si, ci and li; and ExactNumArray nth_root. The "39 UNCLEAR" in Spence's log is the total of all UNCLEAR rows. The other 29 are R4 rows for decompositions, quadrature, ode and roots, and were pinned in v5 (section 21). The LIBRARY.md:399-410 doc row is now ZF-D01 (FIXLIST.csv:141), and all 10 signed-zero rows cite it.
- **Item 5:** see section 0. SRC_COMMIT 7955836 is on paul/zf-bound-merge, not on zf main (dc6bc87).

## 21. Round-5 changes (Spence fourth review items 3-5)
- **No UNCLEAR, no either/or gold.** Every executed row pins exactly one outcome: one value, or one named refusal. The other branch is a FAIL. lint.pl v4 rejects `UNCLEAR` anywhere, and `or`, `either` or a standalone `|` in required_result or refusal_variant, and exits 1 on any BAD line. Each rule was tamper-tested red against a clean green control.
- Class rulings (Spence): gauss_kronrod and roots illinois/brent are correctly rounded at p (Ziv), with PrecisionRetryExhausted only where the gold pins it; brent/illinois with no sign change give Error::InvalidArgument (new `no-sign-change` rows). rk4 is the exact fixed-step recurrence, correctly rounded per step; InvalidArgument only for a bad step or interval (new `bad-step-interval` row). rk45_adaptive is within tol of mpmath odefun at 2p bits in env oracle-2026-10-30, else PrecisionRetryExhausted; it is pinned per row, and rows not yet pinned are `[PM1-UNPINNED]`. LU/QR/SVD/eigen: each factor is correctly rounded under the pinned normalisation; None only for an exact pivot tie past the lowest-index rule, repeated singular values or eigenvalues, or a zero R_ii (new `degenerate` rows). from_raw_parts is exactly normalized when the parts are valid at p, else InvalidArgument (new `raw-parts-invalid` row). Where a row once said "MPFR or mpmath", the pin now names the one oracle from that row's oracle column.
- MAX rows: x86_64 is pinned to the MemoryAllocation (zb Resource) refusal. i686 is pinned to the correct result or enclosure (Paul: the size is allocatable), and a refusal there is a FAIL. (Superseded by v8, section 4: one gold on both targets.)
- Rows that cannot have one outcome yet are `[PM1-UNPINNED]` (rule d): 61 new rows, namely trig-huge-arg e=2^16-1 and e=2^16 (the cap rule is undecided), exponent extremes on helper functions, zero inputs whose doc leaves the result open (for example the zero polynomial and "DivisionByZero or None per doc"), and rk45.
- Row deltas vs 3886199: 24620 -> 24628 (+8 new rows, 0 removed); 3707 rows changed (col10 3637, col11 2157), of which 3646 now carry a single pin and 61 are `[PM1-UNPINNED]`. UNCLEAR 29 -> 0; either/or in cols 10-11 3698 -> 0; `[PM1-UNPINNED]` 996 -> 1057; UNGENERATED pins 19742 -> 19750. The 10 ZF-WP families move N=3 -> N=5; the 4 heun_fox ZF-04 citations now read 29.

## 22. Round-6 changes (Spence review 5 items 1-14)
- **Hedge rule (item 1):** lint.pl v5 and gen_edges.pl share tools/pm_hedge.pl. It rejects in required_result and refusal_variant: TBD, TODO, UNCLEAR, either/or, `, else` and else, unless, otherwise, and `; ` or more than one value in refusal_variant. The `/` rule is scoped to cols 10-11: a `/` is allowed only in a source path, inside parentheses, next to a digit, or in math with `=`, `*` or `^` (for example `4/5`, `cos(pi/n)`, `x=a/b`, a source path `src/x.rs:9`), and any other `/` between alternatives (`a/b`, `Up(+)/FromZero`) is a hedge. Every executed row has exactly one outcome or `[PM1-UNPINNED]`. Tamper proof (pmrev6): against a clean green control (lint rc 0), each of TBD, TODO, slash, `, else`, unless, `; ` in refusal_variant and two values in refusal_variant was red on its own (rc 1, one BAD line each). pm-coverage: a FAIL record is red in rc and design-stage (rc 1), and a duplicate manifest key dies (rc 25) in both modes, each against a green control. The tamper copies were deleted.
- Rulings applied: brent/illinois no-sign-change and rk4 bad-step-interval pin Error::InvalidArgument (section 15 item 5). rk45: all 33 rows cite mpmath odefun at 2p bits, with PrecisionRetryExhausted where tol cannot be met. from_raw_parts: all 20 rows pin one outcome, exact normalised or InvalidArgument. Decompositions: per-class pins on the 132 rows; zero, singular and degenerate are None; no `else ExponentOverflow`. lu: ties go to the lowest index, and None only for an exactly singular matrix. Trig: section 7.
- Row deltas vs 3c0e1e5: 24628 -> 24654 (+26 trig-pi-budget rows); executed 20613 -> 20639. Hedges, executed rows: else 1661 -> 0, slash 5254 -> 0, rv-multi 1006 -> 0. `[PM1-UNPINNED]` executed 1057 -> 1584 (unbounded huge-exp helper families now carry a pin rule). UNGENERATED pins 19750 -> 19776.

## 23. Round-7 changes (Jeff ruling 2026-10-30 6:00 PM CT, FIXLIST:145 ZF-M27)
- **One gold per row:** target-specific rows 18523 before (class 774, gold 18453, pin 1381), 0 after. oracle_i686 now reads "same gold as oracle_x86_64" on every executed row. gen_edges.pl dies on any target-specific row, and lint v6 is BAD on one (tamper: control clean; target gold, target pin and target class each red alone, rc 1).
- **MAX:** 387+387 split rows became 387 `prec-2^30-1` shared rows and 387 `prec-2^60-1` shared rows flagged [PM1-UNPINNED] (section 4). Rule d moved 1584 to 1971 (+387). Superseded by v8 (section 24).
- **M27 red:** 0 rows carried a ZF-M27 pin note before, 16373 after. **Correction (v8):** v7 read "6 directed modes" as every mode except ToEven; the correct set is every mode except None (16725 rows, section 3).
- **FIXLIST re-pin:** sha 9fc4a7a5... to 3a8907f4...; existing line citations unchanged.
- **Row deltas:** rows 24654 to 24654, executed 20639 to 20639, NA_PROOF 4015 to 4015; rule c 19776 to 16022 (rows whose x86_64 gold is pinned no longer wait on a separate i686 gold); rule d 1584 to 1971.

## 24. Round-8 changes and decision log (pmrev8; Jeff 2026-10-30 6:11, 6:13, 6:15 PM CT)
Decision log (each principle applied directly, none left open; PRINCIPLES.md sha 81780e34...):
- **P1 bit-identical everywhere:** every row has one gold, one limit, one refusal variant and one budget for both targets. MAX_PRECISION is 2^30-1 everywhere (ZF-M28). lint v7 now FAILs (BAD) a target-specific class, gold or pin and any target-dependent limit, refusal, budget or result. Lint BAD rows on the e7ad025 table: 0 under v6 lint, 387 dependent + 387 limit under v7 lint (log 1821-lint8-before.log, the prec-2^60-1 rows); on this table: 0 (1821-tamper8-clean-green.log).
- **P2 correct or refuse:** over-limit rows (`prec-2^30`, `prec-2^60-1`) require InvalidArgument (zf Error, zb BoundError at error.rs:22); NaN-as-error, a panic or another variant is FAIL.
- **P3 32-1024 bits plus boundaries +-1:** the boundary 2^30-1 / 2^30 is covered by twin rows on every family (387 new `prec-2^30` rows); the base axis is unchanged.
- **P4 test-first, exact command and full count:** every count below comes from a poprun log, and each lint rejection class has a saved red tamper log.
- **P5 through Paul:** FIXLIST is Paul's file, pinned at 0ef77110...; no new FIXLIST rows are invented here.
- **Execution config:** the 10-minute per-run budget, the i686 cost factor and the worker count are execution config, not expected outcomes. lint allows them because they are the same on both targets; the budget rule only rejects a target name or target-split budget in the tier or verdict.

Row changes (before e7ad025 -> after):
- Rows 24654 -> 25135 (+387 `prec-2^30`, +94 `m30-*`); executed 20639 -> 21120; NA_PROOF 4015 -> 4015.
- ZF-M27 notes 16373 -> 16725; wrong "except ToEven" text 16373 -> 0.
- ZF-M28 cites 0 -> 1161 (387 x 3 prec classes); unpinned `prec-2^60-1` 387 -> 0; zb over-limit rows 70, all BoundError::InvalidArgument.
- ZF-M29 (FIXLIST:147) cites 0 -> 31, in every Error-kind row of the brent (6), illinois (6), rk4 (6) and rk45 (13) families; the "API row" note 12 -> 0.
- ZF-M30 (FIXLIST:148) rows 0 -> 94: `m30-p2^24` and `m30-p2^27` for each of 47 families, mode ToEven, pm-max shard on both targets, one correct gold, red until ZF-M30 lands.
- ZF-D04 (FIXLIST:149): sinh_series (pow.rs:678) and sin_series (sin.rs:94) are internal, not exported, so no matrix row tests them directly (0 rows); they are reached through the sinh/sin M30 rows.
- Coverage rule c 16022 -> 16346, rule d 1971 -> 1584.

M30 families (series, tan and reduction-cost paths; `series_cost_optimize` at atan.rs:145, cos.rs:101, log.rs:193, pow.rs:682, sin.rs:100, tan.rs:121; `calc_mul_cost` util.rs:55; `reduction_cost` at atan.rs:69, cos.rs:29, log.rs:73, pow.rs:32, series.rs:117,159, sin.rs:27, tan.rs:53):
- ExactComplex (15, per-op `*_word`, e.g. complex.rs:407): exp ln sin cos tan sinh cosh tanh log2 log10 log pow asin acos atan.
- ExactNum (17, gen_wrapper_* in ext.rs, e.g. :2745 exp, :2823 tan): pow log atan atan2 tanh sin_cos ln log2 log10 exp sin cos tan asin acos sinh cosh.
- Ball (4, ball.rs:92/106/113/125): exp sin cos ln. ComplexBall (4, ball.rs:325/342/354/372): exp ln sin cos.
- ieee_soft array.rs $arr (3): atan2 pow log. ExactNumArray (4, map_at/map_exact, array.rs:425, :1133): sin_cos atan2 pow log.
- Excluded: ExactInt::pow (integer square-and-multiply, integer.rs:331).

Target-split scan (CSV, shards, this doc): "32-bit", "64-bit", per-target, pm-max-<target>, -only: 0. "on some targets" appears only inside red-until text that describes today's defect.

## 25. Round-9 changes and decision log (pmrev9; Jeff 2026-10-30 6:24 PM CT)
Decision log:
- **FIXLIST re-pin (P5):** 0ef77110... -> 514193c96f9a088617a19ab578a16e41390b219823b4b33d89d6a3979b962ffa (194 lines). Paul added ZF-TD01..ZF-TD45 at :150-194. ZF-M27..ZF-M30 and ZF-D04 are still at :145-149, so no citation moves. Pinned in gen_edges.pl %PIN, inputs/SHA256SUMS, the CSV header and section 0.
- **Tie rows (P1, R3, R6, ZF-M27):** each of the 309 rows that v8 labelled "fixed mode of the family" now has one outcome, taken from inputs/tie-resolution.tsv (pinned in section 0). gen_edges.pl dies if a family is missing from the table.

| Resolution | Rows | Outcome | Ruling |
|---|---|---|---|
| exact | 178 | no rounding step (integer ops, bytes, CSV and serde, predicates, constructors, abs/neg/min/max/clamp, frexp, ExactRational add/sub/floor/ceil, ExactNum rem/pow2/max_value/min_value, from_raw_parts, IEEE bit accessors, to_exact): the tie input yields its exact value | R3, P1 |
| ieee | 54 | ieee_soft $ty/$arr and Ieee32/Ieee64::from_exact: correctly rounded in roundTiesToEven at the format precision (24 or 53 bits); the tie goes to the even significand; the p axis and ZF-M27 do not apply | R3, IEEE 754-2019 |
| mode | 21 | fixed internal ToEven (18), ToZero (2), FromZero (1) at p, read from the call site: rounded once at rp(p,64); the p=32 cell is red until ZF-M27 lands | R3, ZF-M27 (FIXLIST:145) |
| none-R6 | 6 | fixed internal None: RNDZ at rp(p,64) | R6 |
| rm | 35 | the family takes a RoundingMode (gen_wrapper_*_rm* macro in ext.rs, or ctx.rm for Context::const_*): all 7 modes, one outcome per mode; the M27 pin note applies | ZF-M27, R6 |
| r4 | 13 | Spence's R4 row fixes the outcome (accuracy formulas, the 4 IEEE soft_* ops, random_normal, test_rng seeding) | R4 |
| red | 2 | `[PM1-UNPINNED]`: FAIL until the fact below lands | none |

Red tie rows (owner Rusty for both):
- ExactNum::from_word (ext.rs:1913): Rusty has to state which RoundingMode ExactNumNumber::from_word applies when p is below the bit width of d, and whether d: Word, whose width depends on the build, gets one gold on every 64-bit target (ZF-TD36; under P6 Word is 64 bits on every target).
- ExactNumArray::from_csv_str (csvfmt.rs:92): Rusty has to state which RoundingMode the decimal field parser uses at DEFAULT_P.

Before -> after: unmarked fixed-mode tie rows 309 -> 0; zf executed tie rows 645 -> 645 (rm 336 -> 371); ZF-M27 cites in tie rows 336 -> 446 (+35 rm pin notes, +21 mode red sentences, +54 ieee "does not apply" texts); `[PM1-UNPINNED]` 1584 -> 1586; rows 25135 -> 25135; executed 21120 -> 21120.
- **Generator finding:** for 35 families the inventory signature omits RoundingMode because the rm parameter comes from a gen_wrapper macro (for example ext.rs:2665 cbrt, :2745 exp, :2892 erf) or from ctx.rm. v8 therefore labelled them fixed mode, which was wrong. They now have the all-7-modes tie row.
- **Open questions removed:** Q3 (already closed in pmrev5) and Q7 (answered by P1, R3, R6 and ZF-M27 above). Q1, Q2 and Q4-Q6 stay open: P1-P5 do not answer them, because each one is PM-1 or PM-3 work or depends on a branch merge.
- **ZF-D04 FLAG (owner Rusty):** FIXLIST:149 says "Public sinh_series (pow.rs:686) and sin_series (sin.rs:95)". At 7955836 both are `pub fn` on `pub(crate) struct ExactNumNumber` (num.rs:29) inside the private `mod ops` (lib.rs:53), at pow.rs:678 and sin.rs:94. They are not public, and the cited lines are off by 8 and by 1. The matrix keeps 0 direct rows. Rusty confirms the intended text, and Paul edits FIXLIST.

TD citations (FIXLIST:150-194) for matrix content in this design:
- One limit and the MAX_PRECISION classes (section 4; `prec-2^30`, `prec-2^30-1`, `prec-2^60-1`): ZF-TD01 :150, ZB-TD02 :151, ZB-TD03 :152, ZF-TD13 :162, ZF-TD28 :177, ZF-TD29 :178, ZF-TD30 :179.
- Retry budget and PrecisionRetryExhausted refusals (section 5), and None outputs: ZF-TD04 :153, ZF-TD05 :154, ZF-TD16 :165.
- One gold from equal working precision (numeric-CR rows): ZF-TD06 :155, ZF-TD09..ZF-TD12 :158-161, ZF-TD15 :164, ZF-TD17..ZF-TD21 :166-170, ZF-TD23 :172, ZF-TD24 :173.
- Random rows (G13) and seeds: ZF-TD07 :156, ZF-TD08 :157, ZF-TD43 :192.
- The p=32 cells of the ZF-M27 pin note: ZF-TD14 :163.
- Error variant for x >= 2^32 (section 5): ZF-TD22 :171, ZF-TD34 :183.
- ZF-M30 cost rows (`m30-p2^24`, `m30-p2^27`): ZF-TD25 :174, ZF-TD26 :175, ZF-TD27 :176.
- MemoryAllocation is FAIL (section 5): ZF-TD31 :180, ZF-TD42 :191.
- Golds compared by bytes and not by raw words (bytes/raw-parts rows): ZF-TD32 :181, ZF-TD33 :182, ZF-TD37 :186.
- Digest positions: ZF-TD35 :184. Word-independent inputs on every target: ZF-TD36 :185 (the from_word red row).
- Oracle on every OS: ZF-TD40 :189. CRLF checkout and SHA256SUMS: ZF-TD41 :190. Path joins: ZF-TD44 :193. Timing: ZF-TD45 :194.
- ZF-TD38 :187 and ZB-TD39 :188 are docs only and have no matrix content.

Post-approval PM items (owner Harper; not implemented in v9 or v10):
- PM-a ZF-TD35: digest grid positions in u64, independent of usize, with one pinned digest.
- PM-b ZF-TD36: every matrix test runs on every 64-bit target with Word-independent inputs.
- PM-c ZF-TD37: one byte or hex gold per raw-word row, compared through from_raw bits.
- PM-d ZF-TD40: oracle via the committed gold file on Windows and macOS (proposal: ~/harper-work/td40-cross-os-proposal-2026-10-30.md).
- PM-e ZF-TD41: a `.gitattributes` eol=lf rule for inputs, golds and SHA256SUMS, plus a CRLF tamper check.
- PM-f ZF-TD43: the harness records every seed in the run manifest, and env seeds stay in the harness.
- PM-g ZF-TD44: harness paths use Path joins, with a check that no path is built with a literal separator.
- PM-h ZF-TD45: timing stays in the benchmark crate, and no verdict reads a clock other than the 10-minute budget.

Checks (pmrev9, logs in ~/pop-logs/harper/2026-10-30/): lint clean `*-lint9-clean.log`; the tamper set `*-tamper9-*-red.log`; coverage `*-cov9-rc.log` and `*-cov9-ds.log`; the target-split scan `*-scan9.log`.

## 26. Round-10 changes and decision log (P6; Jeff steering 2026-10-30 6:32 PM CT)
P6 (~/pop-logs/PRINCIPLES.md, sha256 d39b9a605a02278116f032e65283da376d7d76e904049bffcfdb35c0267061f0): targets are x86_64 and aarch64 Linux, plus Windows and macOS where they can run. 32-bit is out of scope.

Decision log:
- **FIXLIST re-pin:** 514193c9... (194 lines) -> 56633946854c27c8c702b24157119102fed1e2e7a8044eb18a7111411bce08af (196 lines; log 1836-fixlist-196). New: ZF-M31 :195 (64-bit minimum, compile_error! on a 32-bit target; owner Rusty) and ZF-M32 :196 (re-derive MAX_PRECISION for 64-bit; supersedes ZF-M28 :146; owner Rusty). Every cited line re-checked: ZF-M27 :145, ZF-M28 :146, ZF-M29 :147, ZF-M30 :148, ZF-D04 :149, ZF-TD01 :150 .. ZF-TD45 :194 unchanged. The original drifted to 197 lines (02c209a5..., ZB-D05 at :197, "CONTRACT.md 8.1/15.1 say prec: u32", owner Eugene); not pinned, flagged for Paul.
- **Removed rows (821), by reason:**
  - `prec-2^30-1` 387: p = 2^30-1 was the i686 MAX_PRECISION (ZF-M28), superseded by ZF-M32.
  - `prec-2^30` 387: the ZF-M28 MAX+1 twin; `prec-MAXp1-refuse` already covers MAX+1.
  - `m30-p2^24` 47: existed only because calc_mul_cost overflowed u32 on i686.
- **Kept `m30-p2^27` (47 rows)** with evidence: series.rs:156, 238, 244 and tan.rs:128 compute the cost in 64-bit arithmetic that still overflows at p = 2^27 (ZF-M30, FIXLIST:148).
- **`prec-2^60-1` (387):** now `[PM1-UNPINNED]`, red until Rusty pins the ZF-M32 constant (rule d 1584 -> 1973).
- **ZF-M27 (:145):** met on every 64-bit target, so the 16725 "red until ZF-M27" sentences are replaced by "ZF-M27 is met on every 64-bit target" and the tie p=32 cells are not red.
- **ZF-M28 (:146):** cited only as "supersedes ZF-M28" on `prec-2^60-1`.
- **prec axis:** the W=32 word points are dropped; the 49-point 64-bit set remains.
- **Header:** `oracle_i686` -> `oracle_aarch64` ("same gold as oracle_x86_64; aarch64 runnable after TD40 installs").
- **Budget:** i686 factor removed; no aarch64 factor (qemu not installed, nothing measured).
- **Lint:** pm_p6 (tools/pm_target.pl) FAILs i686, 32-bit, thumbv*, armv7* and target_pointer_width="32" in any row or the header. Tampers p6target, p6word, p6header are red.
- **ZF-D04 (:149)** mentions |x| >= 2^32 on i686. Design note: under P6 the row content is word-independent and needs no i686 case; the FIXLIST text is Paul's to edit.
- **ZF-M31 test:** rustup keeps i686-unknown-linux-musl installed only so ZF-M31's red check (`cargo check` must fail with compile_error!) can run; it is not a matrix target.

Before -> after (cb744c42 -> this commit): rows 25135 -> 24314; executed 21120 -> 20299; NA_PROOF 4015 -> 4015; pm-max rows 1255 -> 434; fast 19839, heavy 26 and proof 4015 unchanged; families 835 -> 835; target-split scan hits 38916 -> 0. Coverage rc: a 658 -> 658, b 3357 -> 3357, c 16346 -> 15764, d 1584 -> 1973, e 0 -> 0.

Checks (logs in ~/pop-logs/harper/2026-10-30/): gen 1845-gen10c; lint clean 1848-lint10c; tampers 1846/1847-lint10-tamper-* and 1847/1848-lint10b-p6*; budget 1846-budget10b; coverage 1848-cov10-rc and 1848-cov10-design-stage; counts and scan 1848-count10b.
