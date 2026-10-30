# zf/zb precision matrix: design v2 (revised after Spence REQUEST_CHANGES)

Agent: Harper. Date: 2026-10-30. Edge table: `precision-matrix-edges.csv` (20659 rows, 17 columns).
Generator: `gen_edges.pl` (kept in ~/harper-work/pm-rev). Lint: `lint.pl` (it is clean: no BAD lines).

## 0. Inputs
- The zf inventory is at ~/paul-work/zf-audit-2026-10-30/inventory.csv. Rows that are reachable or macro-impl are kept, and trait-decl rows are dropped. That gives 719 families.
- The zb inventory is section 1 of ~/paul-work/zb-audit-2026-10-30/inventory.csv. That gives 116 families.
- The R4 required results come from Spence's r4-required-results.csv, which I copied unedited (sha256 2746dc4d…). All 141 "other" families are mapped, and none are missing.
- Oracle env: ~/team-env/oracle-2026-10-30. It holds mpmath 1.3.0, sympy 1.14.0, gmpy2 2.3.2 and python-flint 0.9.0. Rust uses MPFR 4.2.2 and GMP 6.3.0, through rug 1.20.1 / gmp-mpfr-sys 1.6.8 with feature mpfr-tests, run on x86_64.

## 1. Verdicts (findings 8, 9)
- The only verdicts are PASS, FAIL and NA_PROOF. There are no waivers, and an escalation counts as FAIL.
- An NA_PROOF row needs a committed proof at docs/test-plan/pm-proofs/<crate>/<key>.md, signed by a second agent. pm-coverage FAILs any NA_PROOF row whose proof is missing or unsigned. It also FAILs any row whose gold is `sha256:UNGENERATED`.
- There are 4324 NA_PROOF rows (zf 3431, zb 893) and 16335 PASS|FAIL rows (zf 14532, zb 1803).

## 2. Family keys (R1)
Families are keyed by `Type::name` from impl_context. The CSV adds a `DEFAULT_RANDOM_SEED` row.

## 3. Rounding modes (finding 1)
- There are 7 modes: ToEven, ToZero, Up, Down, FromZero, ToOdd and None.
- For ToOdd, the oracle computes RNDN at p+1 extra and then flips exact ties to the odd neighbour itself, because MPFR has no RNDO.
- None follows R6: it means RNDZ of the exact value at rp(p,64).

## 4. Precision axis (finding 2)
- The base axis is p ∈ {32,64,128,256,512,1024}. The heavy axis has 96 points.
- Rows for p=1, 2 and 3 run in all 7 modes.
- `prec-MAX` is p=MAX_PRECISION, which is 2^60-1 on x86_64 and 2^30-1 on i686. The result is either correct or `Error::MemoryAllocation`.
- The refusal rows p=0 and MAX+1 both give InvalidArgument. APIs that return Option refuse with None.

## 5. Error variants (finding 3)
- zf `Error`: ExponentOverflow(Sign), DivisionByZero, InvalidArgument, PrecisionRetryExhausted, MemoryAllocation.
- zb `BoundError`: Representation, InvalidArgument, Domain, Resource, NotPositiveDefinite, Singular, NotConverged, DimensionMismatch.
- The lint rejects these in the CSV: Domain or Resource in zf rows, bare Overflow or Underflow, and BoundError::Overflow.

## 6. Thresholds (finding 5)
- Classes are pinned at EXPONENT_MIN and EXPONENT_MAX (EXPONENT_MAX = 2^62-1, defs.rs:49): ovf-max-finite, ovf-max-plus-half-ulp, unf-min-normal, unf-min-subnormal and unf-half-min-subnormal-tie.
- The i64 extremes appear only in the scalb, ldexp and ilogb rows (7 rows, class i64-extreme-exponent-arg).
- The T1–T12 table from v1 is kept. Each threshold T gets T-1, T and T+1. That is 140 survey cells on x86_64 and 123 on i686.

## 7. Huge-argument trig (finding 4)
- Class `trig-huge-arg` has 48 rows. The exponents are 2^16-1, 2^16 and 2^16+1, where trig_reduce_extra caps (ops/util.rs:92-95), plus EXPONENT_MAX.
- The required result is the correct bits or PrecisionRetryExhausted, and it must be identical on x86_64 and i686.
- FIXLIST Core: make the huge-argument reduction exact past 2^16, or document the refusal.

## 8. Fast tier (finding 6)
The fast tier crosses the huge, tiny, signed-zero, tie and cancellation classes with all 7 modes at p=64. Each rm-* row is also 1 cell at p=64.

## 9. NaN pairs (finding 11)
Binary zf functions get nan-AB, nan-BA, nan-AN and nan-NA rows (216 each), following Rule R in nan_error_kinds.rs.

## 10. T4 mul operand pairs (finding 12)
Class `mul-pair-T4` has 33 rows. sm is 69, 70 or 71 limbs, and lg = 2·sm, in all 7 modes at p=64.

## 11. Oracles and PMCELL (finding 7, R5)
- Each row pins one oracle per target in `oracle_x86_64` and `oracle_i686`. The oracle is either a version string or a gold file plus SHA256. On i686, golds and digests are cut from the x86_64 run.
- `ORACLE_PY` points only to ~/team-env/oracle-2026-10-30.
- Each cell prints one line in this format:
  `PMCELL crate=<zf|zb> family=<Type::name> class=<c> p=<p> rm=<m> target=<triple> oracle=<id> env=<oracle-env-path|rug-1.20.1/mpfr-4.2.2> verdict=<PASS|FAIL|NA_PROOF> ms=<t>`
- **ZF-PM5:** each release candidate gets N=3 independent spot-checks per family. They are random cells recomputed in mpmath and checked against MPFR or the gold, by an agent other than the one who cut the golds. Any disagreement is a FAIL.

## 12. Tiers (R7)
- Fast tier: runs on every commit.
- Batch gate: the fast tier plus a cut heavy tier, which is the base 6 p × all modes plus the threshold points.
- Full heavy tier: the 96-point axis × all modes. It runs once per release candidate, overnight.

## 13. Budgets (recomputed from the CSV by budget.pl)
- Cost basis (v1, not re-measured): zf 4.1 ms per cell and zb 12 ms per cell, with i686 at 4x. Measured i686 is 1.9x, and those figures are shown in parentheses.
- Cells per row:
  - default rows: fast 6, cut 6×m, full 96×m, where m=7 if the row runs all modes
  - rm-* rows: 1 / 6 / 96
  - refusal rows: 1 / 1 / 1
  - p=1,2,3 rows: 3×m for each tier
  - thresholds: fast 1 cell per point, heavy 7 cells per point
  - prec-MAX rows (zf 352, zb 35) and the doc row are budgeted separately

| crate | target | fast cells | fast | cut cells | cut heavy | full cells | full heavy |
|---|---|---|---|---|---|---|---|
| zf | x86_64 | 77213 | 5.3 min | 475102 | 0.54 h | 7472512 | 8.5 h |
| zf | i686 | 77196 | 21.1 min (10.0) | 474983 | 2.16 h (1.03) | 7472393 | 34.0 h (16.2) |
| zb | x86_64 | 9553 | 1.9 min | 65701 | 0.22 h | 1048591 | 3.5 h |
| zb | i686 | 9553 | 7.6 min (3.6) | 65701 | 0.88 h (0.42) | 1048591 | 14.0 h (6.6) |

The full heavy tier on i686 does not fit in one overnight run on a single worker. It needs about 4 parallel workers at 4x, or about 2 at 1.9x.

## 14. Rulings
- **R2:** test_rng stays under the non-default `random` feature. Its 5 names are at lib.rs:204-207, and it is group G13.
- **R3:** each family is correct to p or refuses with an input error. The 10 working-precision families are FIXLIST items: assoc_legendre_p, chebyshev_t, chebyshev_u, gegenbauer, gen_laguerre, hermite_h, hermite_he, jacobi_p, laguerre and legendre_p. That is 10, not 9, because gegenbauer had been miscategorised.
- **R4:** see section 0.
- **R8:** Ball::div rows 18150 and 18151 use Spence's exact text: "0 in the divisor with rad>0 gives all reals" and "an exact 0 divisor gives Domain". CONTRACT 4.2 and 10.3 and PLAN:287 were updated in commit 5649ad0.
- **Primality:** row `doc-probable-prime` for miller_rabin. Its docs (modular.rs:105-108) say neither "probable prime" nor that it is deterministic, so the row FAILs. FIXLIST: fix the doc.

## 15. FIXLIST additions
1. The 10 working-precision families (R3).
2. The trig huge-argument row (Core).
3. The miller_rabin doc fix.

## 16. Open questions
1. MAX_PRECISION on i686: 2^30 bits can actually be allocated, so the transcendental prec-MAX rows will likely time out, and a timeout is a FAIL. Should those rows get a cap or a separate budget?
2. The per-cell floor and these budgets use v1 costs. They will be reset after the PM-0 cost runs finish; those runs are still in progress and I did not use them.
3. The UNGENERATED golds are pinned at PM-3. Until then, pm-coverage FAILs those rows.
4. Full heavy on i686 needs parallel workers (section 13). We need to confirm how many runners are available.
5. ZF-PM5 uses N=3. Should it be higher?
