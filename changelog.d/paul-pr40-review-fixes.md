### Changed (Spence PR #40 REQUEST_CHANGES at e6f439f, paul, 2026-11-22)

- zenith-float-bound CBall: B3 bottom ε for a midpoint product lost below EXPONENT_MIN; B4 exact 2^-k scaling in
  acb_inv / acb_div so c² + d² stays in range; R1 1,766 FLINT 3.6.0 edge rows (cball_arb_edge_rows.txt, exponent ends,
  lagom, accuracy/fits/hprec cutoffs, exhaustive p = 2..8) with three root-cause fixes (fused_mid past EXPONENT_MAX,
  arb_addmul Inf product, mag_again canonical rebuild); B8 rules_scan helper rename; B2 prep mutate.sh workspace
  layout (selftest case 44). Details in zenith-float-bound/CHANGELOG.md.
- zenith-float-num: B5 complex J_ν(0) / I_ν(0), real non-integer ν < 0, are the signed one-sided limit
  sgn Γ(ν+1)·∞ at all four signed zeros (Jeff 2026-11-22 8:33 PM CT; mpmath hides the sign, land39 RC-14); R2 Carlson
  stop decided by FLINT's rigorous truncation bound, 4^-t homogeneity normalization, boundary tests, and the
  "Arb uses dmax/|A|" wording corrected; B7 no f64 in library parameter choices (common::fixlog), no unwrap/expect,
  workspace clippy -D warnings clean; the three arm-tests --lib failures at e6f439f (stale expectations: struve/anger/weber
  at p = 2^59 is MemoryAllocation, MPFR 40-digit display, J_{-1/2}(0) limit). Details in CHANGELOG.md.
- B6: docs/agent-trailer-map.tsv attributes the 70 post-cutoff commits without an Agent trailer (46731eb … 856f82c,
  including merge 083b930 and the non-bisectable WIP bf82e6b, M2) from Paul's bare-repo reflogs; no history rewrite
  (Jeff 2026-11-22). tools/fold_pr1_selftest.sh row P honors the map read from the checked ref; new plant P-MAP (a map
  row with a non-agent name) is DETECTED.
- M1: changelog.d/paul-tier1-fixes.md no longer calls Γ at p ≥ 2048, ZF-16 and ZF-60 BROKEN (Tier 1 closed, 0 wrong).
- fmt: `cargo fmt --all` over the 16 zenith-float-num files that failed `cargo fmt --all -- --check` at e6f439f (gate step `fmt`); formatting only (Jeff 2026-11-22 ~11:58 PM CT).
