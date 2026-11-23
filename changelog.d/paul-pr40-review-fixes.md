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
- nan_error_kinds: the nine `ln_gamma` pins now follow the P20 contract of f0d0b47 (a NaN in either part of a complex
  log-Gamma input gives a typed NaN in both parts: real part's kind, else imaginary part's, else InvalidArgument). The
  test was stale (it still pinned per-part kinds and typeless `None`); the code is right. Test-only change.
- rules_scan: the borrow call graph no longer dies with IndexError in `_callees` on a macro metavariable as qualified
  self type (`<$bits>::add_lanes(..)`, ieee_soft/array.rs). `<$t>::f` now resolves through the macro's invocations
  (`impl_ieee_array!(.., u32)` -> `u32::f`), else every method `f`. Selftest fixture v-macro-metavar-qself (red: IndexError;
  fixture totals 30 -> 31 / 74 -> 75). The traceback hid the whole borrow check: with it fixed the gate form reports
  52 unlisted copy sites (reach 246 -> 450 fns since 25e2d82); those and the clone / test-panic / neg_via_mul counts
  above baseline are escalated, not baselined.
- Complex 2F1, integer c − a − b = m < 0 near z = 1: Euler's step to c − a′ − b′ = −m > 0 now lets the inner call use
  the logarithmic 1 − z connection again. Since 404a9b5 (first bad commit by bisect be63db9..e6f439f) the connection
  was marked tried before the Euler step, so the inner call ran the plain series at |z| = 1 − 1e-7 (~10^9 terms) and
  c2f1_log_mpmath never finished (>10 min; 0.43 s at be63db9). Now 0.30 s, all 54 mpmath rows pass. New test
  log_case_negative_m_near_one_returns (60 s budget; red: timed out).
- special_audit_mpmath: the two `c.ln_gamma` rows (3−4i at p = 128, 20+20i at p = 512) still held the 1.0.5 principal
  ln(Γ(z)) imaginary parts (1.5405…, −0.6726…). Complex ln Γ is Arb's continuous log-gamma since Jeff's 2026-11-18
  12:14 CT ruling; the code matches FLINT 0.9.0 `acb_lgamma` (−4.7427…, 62.1593…; mpmath `loggamma` agrees), 2πk from
  the old values, so "−2 / −1 correct bits" was a stale reference, not a 1–2 bit error. Test-only change.
- Complex Bessel I_ν / K_ν at |z| = 20000 (complex_bessel_mpmath, 12 rows at p = 128 and 256): the Arb 0F1 term-count
  walk (`choose_n_0f1`, 7b972c3) stopped at Arb's `n_max = 50 + 10 prec` and then returned `n_best = 1`, so the series
  was the first term only (I_{1/2}(20000) = 112.8 instead of 2.19e8683). The walk now runs until the terms have
  decreased 2^-(prec+4) (no cap; NaN PrecisionRetryExhausted only past u32::MAX terms), as the "Caps" rulings require.
  Values match FLINT 0.9.0 `bessel_i` / `bessel_k` and the mpmath refs. New unit test zf_0f1_large_z_runs_to_convergence
  (0F1(;3/2;10^8) = sinh(20000)/20000; red: got 1).
- Unit-test wall-clock budgets scale with `ZF_TIME_FACTOR` (new `#[cfg(test)] mod test_time`, same grammar and [1, 8]
  range as s9_specials_golds): zf_ln_gamma_modest_z_at_p8192_finishes (30 s native) and zf_53_pollard_rho_proof_path
  (2^255−19 in 1 s native) failed only under qemu-aarch64 (5066 ms for the 1 s bound). Native budgets are unchanged;
  under the 6.50 pin they are 195 s and 6.5 s. Green on qemu-aarch64 with ZF_TIME_FACTOR=6.50.
- B5 follow-up (Jeff 2026-11-22 ~11:58 PM CT): J_ν / I_ν at z = +0±0i, real non-integer ν < 0, give the imaginary zero
  −sgn Γ(ν+1)·sgn(Im z)·0 (was: the sign of Im z), matching FLINT 3.6.0 along z = r(1 ± i·r²). So J/I_{−1/2}(+0+0i) is
  +∞−0i. This changes 98 official T2 rows (49 J + 49 I, ν = −1/2, z = +0+0i, 7 p × 7 modes); land39 REF-CORRECTIONS RC-14a
  records the alt-judge basis (official judge files unedited). Tests b5_j_i_negative_order_at_signed_zeros_is_signed_limit,
  zf_j_i_neg_half_at_zero_is_pos_inf_and_small_z_matches_flint, zf_65_… and the dest-2^59 lib test updated first (red).
- zenith-float-bound mutate selftests (Spence: pop-logs conflict). The logged()/zf_special fixture runs of mutate.sh exited 1 under
  the current poplog: they inherited POPLOG_OUTER_* (a "WARNING ignored POPLOG_OUTER_* marker" top-level run), the fixture repo
  ignored rather than committed Cargo.lock (NOT A GATE: tree dirty at start), and mutate_selftest.sh ended inside the fixture dir
  (its own run log: NOT A GATE: cwd changed mid-run). The fixtures now unset POPLOG_OUTER_*, commit a generated Cargo.lock, and
  the selftest returns to its start dir before exit. Red at fcb2ec6 in a clean clone; mutate_selftest 44 cases ok, PASS footer.
