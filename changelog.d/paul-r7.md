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
