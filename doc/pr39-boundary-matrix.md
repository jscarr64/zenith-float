# PR #39 boundary matrix (ACCEPTANCE-pr39 I9, T2; P3)

One row per math function whose code changed in #39 outside the fold merges, one column per boundary class. Each
cell names the test that covers it (file and test function) or gives a written N/A reason. The 32-bit limit tier is
void (C1/C3: zenith-float 2.0 is 64-bit only). Abbreviations for the test files:

- TPE = `zenith-float-num/tests/top_precision_edges.rs`
- PB = `zenith-float-num/tests/mpfr/precision_boundaries.rs` (MPFR bit for bit, mpfr-tests feature)
- SPE = `zenith-float-num/tests/small_p_exhaustive.rs` (exact integer oracle, every operand at p 1..=8)
- ER = `zenith-float-num/tests/exponent_range.rs`, OM = `zenith-float-num/tests/overflow_modes.rs`
- PL = `zenith-float-num/tests/precision_limits.rs`, NUM = `zenith-float-num/src/num.rs` unit tests
- ENC = `zenith-float-bound/src/enclosure.rs` `mod contains_endpoints`, EB = `zenith-float-bound/tests/enclosure/ball.rs`
- BR = `zenith-float-bound/tests/harness/slice_bridge.rs`, S2E = `.../slice_bridge/spence_s2_diff_overflow_edge.rs`
- BP = `zenith-float-bound/tests/ball/precision.rs`

Changed functions (commit): R1 `ExactNum::add`/`sub`/`mul`/`div` and their `_ref` forms, public precision check
(8da300e); R2 every other public `ExactNum` entry that takes `p` (the same check: `sqr`, `div_exact` siblings,
`sqrt`, `hypot`, `fma`, `log`, `reciprocal`, `powi`, `parse`, conversions; 8da300e); R3 `Ball::contains` (4e82e1b);
R4 `Mag::cmp_exact_diff` above `EXPONENT_MAX` (7163066); R5 `MAX_PRECISION` (d15bb37, a constant). The commits
bdd8ed4 and 5871b06 add tests only (kill MISSED mutants in `Ball::powi`/`Ball::hypot`); no function changed there.

| class | R1 add/sub/mul/div | R2 other public entries | R3 Ball::contains | R4 Mag::cmp_exact_diff | R5 MAX_PRECISION |
|---|---|---|---|---|---|
| 0, ±1, MIN, MAX values | SPE `add_sub_mul_div_every_small_p_operand_pair_in_7_modes` (0, ±1 and every p-bit value); ER `arithmetic_at_the_top_of_the_range` (max_value), `arithmetic_at_the_bottom_of_the_range` (min_positive) | TPE `sqr_div_exact_and_siblings_at_the_top_are_never_invalid` (exact small values 1, 4, 16 ...) | ENC `one_endpoint_verdict` (d = 0, ±below/at/above rad); EB `contains_decides_exact_cases` | BR `cmp_exact_diff_rows` (small values); S2E `s2_diff_overflow_edge_radius_near_emax` (values at EXPONENT_MAX) | PL `max_precision_is_fixed_2_pow_60_minus_1` |
| both signs, signed zero | SPE core window (both signs, +0 and -0 operands, sign of exact-zero results checked in all 7 modes) | TPE `binary_ops_at_small_precisions_around_word_boundaries_are_exact` (four sign combinations) | ENC `one_endpoint_verdict` (both signs of d, 0) | S2E (both signs of a) | N/A: a constant has no sign |
| ±Inf / NaN (or refusal) | OM `overflow_follows_rounding_mode` (Inf results, 7 modes); SPE overflow passes (Inf vs largest finite) | N/A: R2 changed only the precision check; Inf/NaN operand paths are unchanged (CHANGELOG 2.0 entries and their tests) | ENC `one_endpoint_verdict` (Inf of both signs and NaN give False) ; `every_pair_of_endpoint_verdicts` (rad Inf) | N/A: Mag operands are finite or Inf by type; Inf radius in S2E is not reachable (radius exponents E..E+2) — the Inf branch is BR `cmp_exact_diff_rows` | N/A |
| subnormal floor | ER `arithmetic_at_the_bottom_of_the_range` (denormal results, 6 modes) | N/A: precision check only | N/A: contains compares d = x - mid against rad exactly (sub_full_prec), no rounding floor | N/A: exact comparison, no floor | N/A |
| EXPONENT_MAX, EXPONENT_MIN ±1 | ER `arithmetic_at_the_top_of_the_range` (EXPONENT_MAX - 2..EXPONENT_MAX), `ldexp_reaches_both_limits`; SPE `add_sub_at_the_top_of_the_range_overflow_in_7_modes`, `mul_div_at_the_top_of_the_range_overflow_in_7_modes` (EXPONENT_MAX and +1 via overflow) | N/A: precision check only | EB `contains_stays_unknown_when_rad_is_past_zenith_float_range` | S2E (radius exponents E, E+1, E+2) | N/A |
| precisions 32, 64, 128, 256, 512, 1024 | PB `binary_ops_match_mpfr_at_every_small_precision_and_word_boundary` (1..=130, 191..=193, 255..=257, 511..=513, 1023..=1025) | TPE `every_public_entry_accepts_a_top_precision` (top); PB (the four ops) | EB `contains_works_at_the_wider_operand_precision` | BR rows at their gold precisions | N/A |
| word boundaries 63, 64, 65, 127, 128, 129 | PB (1..=130); TPE `binary_ops_at_small_precisions_around_word_boundaries_are_exact` (w·k ± 1, k = 1..4) | TPE `top_precisions` (MAX - w ± 1, MAX - 2w ± 1, the top multiple of w ± 1) | N/A: contains has no precision argument (exact difference) | N/A: exact comparison | N/A |
| method thresholds ±1 | PB comment: the two branches of the exact-rounding helper (p a multiple of the word size or not) are each side of every word boundary above | TPE `top_precisions` (working width = p + 64 crossing MAX_PRECISION) ; NUM `working_width_is_bounded_by_the_exponent_arithmetic` (2^61 bound and 2^61 + 1) | N/A: no method switch | S2E: the overflow branch (`cmp_overflow_diff`) against the finite one (BR) | N/A |
| MIN_PRECISION = 1, MAX_PRECISION, MAX ± 1 | TPE `binary_ops_at_the_top_of_the_precision_range_are_never_invalid` (MAX, MAX - 1 ...), `precisions_outside_the_range_are_invalid` (0, MAX + 1); SPE (p = 1) | TPE `every_public_entry_accepts_a_top_precision`, `every_public_entry_rejects_a_precision_outside_the_range`, `parse_with_usize_max_takes_the_precision_of_the_numeral` | N/A: no precision argument | N/A: no precision argument | PL `max_precision_is_the_largest_accepted_precision`, `precision_zero_is_rejected`; BP `precision_limits_match_zf`, `constructors_reject_out_of_range_precision` |
| all 7 rounding modes | SPE (7 modes); TPE (MODES, 7); PB (6 modes vs MPFR plus None vs MPFR toward zero at storage) | TPE (7 modes) | N/A: contains takes no rounding mode (CONTRACT O19/O20) | N/A: exact, no rounding mode | N/A |
| first and last loop iterations | SPE (the whole p 1..=8 range, first and last operand of each binade) | TPE `every_public_entry_accepts_a_top_precision` (one call per entry, no loop in the check) | ENC `every_pair_of_endpoint_verdicts` (all 9 verdict pairs, the full table) | N/A: no loop | N/A |

Spot check: `git grep -n '<test name>'` finds each named function in the file given.
