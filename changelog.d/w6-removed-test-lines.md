# w6: removed test lines audit (base abb2d5a = merge-base origin/main .. rusty/w6)

Scope: 696 removed lines in tests matching `assert|floor|#[test]`. Data: tn_base.txt, tn_head.txt, tn_gone.txt, rm_by_commit.tsv (this dir).

## Test names that disappeared (2)
- `gk_nonfinite_and_reversed_bounds_are_errors`: renamed to `gk_nonfinite_and_reversed_endpoints_are_errors` (867456d). Coverage kept.
- `test_sinh_cosh`: removed with the internal `ExactNumNumber::sinh_cosh` (ops/sinh_cosh.rs) in b73337a. The public sinh_cosh is still covered by exact_bit_sinh_cosh_pair (tests/mpfr/exact_bit_functions.rs:520), mpfr compare_ops_test.rs:710 and precision_zero_array_api.rs:518. ExponentOverflow cases were NOT covered (overflow_modes inputs are only 2^32..2^33); added tests/sinh_cosh_overflow.rs in 907977c (216 overflow checks: ±3·2^60, ±2^63, ±max at p 64/128/320, six modes, ±Inf or ±max by mode; green first, mutation red 36/216). Its zero and subnormal asserts are restored as sinh_cosh_at_zero_and_at_the_least_subnormal in the same file.

## Large removals that are rewrites (removed about equal to added)
520029b (254/254, unwrap → ?), a3e3388 (245/245, #[allow] removal), 542a0a6, a5a9163/2310060 (i64 exponent), 0e4cdc5, d12a156, 5c5d3a8, a05b9e6 (24 digests re-pinned), f7de24b.

## Net losses, each checked
- db9c84c (21/14): RawParts migration; the asserts are consolidated in from_raw_parts_checks.
- 48e1568 (9/1): abs_via_mul_equals_abs_for_every_kind_of_operand is gone because abs_via_mul was removed.
- ceb15e7/7e1b7d3 (8/5): p=0 is now rejected; the zero-precision asserts became matches!(Err).
- 79b9d34 (2/1): two atan2 wide_primitives tests merged into atan2_over_the_whole_range.
- 40fdbfb (1/0): the 32-bit `!is_inline` assert was dropped because the value is inline on every target.
- b73337a (9/0): dead code (see test_sinh_cosh above).

## w6 tip (this session)
r14_consts golds were not removed. They moved into the fast and full tiers, and the full tier still checks all 1050 rows plus a pinned digest.

Conclusion: no unintentional coverage loss found. The sinh_cosh overflow item is closed by 907977c.

## Run notes
- r14 fix runs, 2026-10-24: 1906-r14_fix (release, full tier) ended with exit 143 (SIGTERM during golds_full_tier_digest_and_brackets, killed, not a test failure); the 1908-r14_fix_dbg rerun passed 7/7 in 2m48s.
