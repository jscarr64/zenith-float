### Changed (PR #40 r7, paul, 2026-11-25; Paul 11:05 PM CT round 7, Jeff 4:42 PM CT rule)

- Clones (P16): 291 -> 204. New `try_rounded_to` / `try_rounded_exact` / `try_rounded_checked` on the number layer
  and `ExactComplex::try_rounded_to` return a rounded owned value (bit-identical to clone + set_precision /
  try_set_precision, including the error and the inexact flag); `ExactNum::rounded_to` is public.
  `meijer_g_inner` is generic over `Borrow<ExactNum>`; the `zf!` / `cplx!` macros convert through new
  `FromExt<&ExactNum>` / `FromExt<&ExactComplex>` impls; `rk45_at_prec` pushes each accepted step when it is
  replaced; complex_airy builds its two rotations fresh; Legendre kernels return `Option` (None = P1(x) = x).
- `tools/rules_scan.py`: a `[clone]` allowlist entry whose line no longer clones fails (stale entries are deleted,
  88 removed); new `ignore` count check fails on any compiled `#[ignore]` / `cfg_attr(.., ignore)` (now 0).
- Oracle tests: the ORACLE_TIER gate is gone, every row runs every time. scripts/gen_specials_goldens.py uses the
  continuous log-gamma (mpmath `loggamma`, Jeff 2026-11-18 12:14 CT ruling) and the two deep `c.ln_gamma` rows of
  goldens/specials/complex.jsonl were regenerated (they held the principal branch and were hidden by the gate).
- Batch 2: clones 204 -> 196. Oracle test rows are moved into their timed thread (specials_oracle gets the row
  back with the value for its verdict); rational gcd borrows num/den (abs_ref); poly div_rem and the trig identity
  path use rounded_to / try_rounded_to (their documented clone + set_precision equivalents); complex_elliptic's
  stall check keeps the remaining distance by move (from the divisor via into_parts).
- Batch 3 (Paul 12:32 AM CT): tools/full_gate.sh has one `oracle` step (the ORACLE_TIER `oracle-deep` duplicate is
  gone; selftests updated). Gate selftests that were red at HEAD are green: time_limits.tsv re-pinned (the two
  removed mantissa timing tests out, the rem_pi row wait in), the rem_pi test waits through `test_time::wait_row`
  and covers all 7 modes, tests/mpfr/mp.rs drops its f64 conversions (fuzz_round_modes checks the MPFR product
  with integer bounds), `ci_no_intrinsics.sh` also catches suffixed float literals (`1.5_f64`), and the
  rules-enforcement fixtures carry the unlisted-clone baseline.
- Batch 4 (Paul 12:32 AM CT items 3, 4): `riemann_zeta(1)` is `+Inf`, exact, divideByZero (POLE-08); the 12 held
  zeta(1) golden rows (POLE08_ZETA_HELD, s9 EXCLUDED_IDS) are gone and run every time. The golden generators give the
  ZF-M27 ruled `None` overflow (signed MAX at the storage precision, inexact; gen_oracle_goldens.py and
  mpfr_oracle.c agree), so edges.jsonl's 60 F-B rows now expect MAX. The library still returns +-Inf there, so
  those 60 rows stay held: they run every time and must still give +-Inf (the test fails when one matches its
  golden). gen_specials_goldens.py resolves each row with a Ziv loop (it used to inherit the previous row's
  precision) and carries the right E(fl(pi/2)|1) note for p=512; a full regeneration of all three generators
  is byte-identical to the committed goldens (checked against FLINT arb and mpmath).
- Batch 5 (Paul 12:32 AM CT item 2, no clones by another name): the helper-copy count also counts the rounding
  helpers that write a new value from a borrowed one (`rounded_to`, `try_rounded_to`, `try_rounded_exact`,
  `try_rounded_checked`, `finish_ref`); a helper's delegation inside its own definition is counted once, at the
  outer call (no_clone selftest covers both, and a mutant that counts delegation fails it). Re-measured: helper
  copies 317 -> 443, combined 513 -> 639 (the same tree measured wider); the per-commit rule still holds on every
  commit of the range. Four owned sites no longer copy (catalog hypergeometric denominator, catalog_aws `u0`, two
  tests). The other 129 rounding-helper uses each read a borrowed source: 66 inputs, 34 `&self` kernels, 12
  goldens rounded per mode, 10 caches, 7 helper equivalence tests.
- Batch 6 (Paul 12:32 AM CT item 1, equal-to-input shortcuts): `ExactInt::shl` / `shr_towards_zero` take `self` by
  value (a zero value or shift returns it, moved); `shl_ref` / `shr_towards_zero_ref` return a `Cow` that borrows
  the input on those shortcuts; bit tests read `magnitude_bit` instead of shifting. modular's `abs_int` returns a
  `Cow` (borrowed when non-negative), owned values use `into_abs`, `gcd` callers drop the `|.|` (gcd reads
  magnitudes), and `egcd` is `bezout_if_coprime` (only mod_inv's coefficients are built; the gcd is only tested).
  exp_bb's last chunk borrows t. Clones 196 -> 191.
- Batch 7 (item 1, part 2): a zero result equal to a zero input is built fresh (`ExactNumNumber::zero_like`: same
  sign, exponent, precision and flags, no payload to copy) in conv's two integer roundings, `rem(0, d)` and
  bernoulli_rev; aprcl's `nth_root`, complex_special's `add_ui_c` and complex_hypergeom's `z_on_lower_cut` return a
  `Cow` that borrows the input on their equal-to-input shortcuts. Clones 191 -> 183.
- Batch 8 (item 1, part 3): complex `elliptic_e_word` takes its caller's working copy by value, so `E(z | 1) = z`
  returns it as is. Clones 183 -> 182. The other equal-to-input shortcuts (ExactInt `add` x2 / `gcd` x2, number
  and ExactNum `rem`, poly `primitive_monic` x2) stay: a by-value probe broke ~45 callers that still need the
  receiver, so the signature change would add more copies than it removes.
- Batch 9: radix_roundtrip's negative mirror is built as `from_i32(-k, p)` (its `[clone] TODO` keep is gone).
  Clones 182 -> 181.
- Batch 10 (Paul 3:36 AM CT): the 8 public-API equal-to-input keeps (ExactInt add x2 / gcd x2, number and ExactNum
  rem, poly primitive_monic x2) state their no-other-way reason: the by-value signature was probed and adds copies.
- Batch 11 (F-B, ZF-M27, Paul 3:36 AM CT): a `None` overflow of finite operands now saturates to the signed largest
  value at rp(p, 64), flagged inexact (rulings M27-FB / M27-OVF); `add/sub/mul_full_prec` saturate at the exact
  length. Internal steps run in the new doc-hidden `RoundingMode::Work`, which keeps +-Inf as the overflow signal
  (`*_full_work` for the full-precision forms). One predicate, `RoundingMode::saturates_overflow`, serves every
  directional overflow site. The 60 held oracle rows and every NONE-GAP F-B skip are asserted directly; six grid
  digests re-pinned with the saturated cells in; seven MPFR golden streams regenerated.
- Batch 12 (F-12, ZF-M27, POW-01): `pow(x, +-0) = 1` and `pow(+1, y) = 1` for any `x`, `y`, NaN included, and
  `powi(NaN, 0) = 1` (IEEE 754). The 84 held grid cells are asserted and hashed; the exp/pow digests re-pinned.
- Batch 13 (F-26/F-27, ZF-M27): `log(x, base)` gives the signed poles of `ln x / ln b`: `log_b(+-0)` is +Inf for
  `0 < b < 1` and -Inf for `b >= 1`; `log_1(x)` is -Inf for `0 < x < 1`, +Inf for `x > 1`.
