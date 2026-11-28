### Review fixes (land39, Spence REQUEST_CHANGES on 11377ab)

- Remove the Tier 1 π / Jacobi K reduction exponent limits and the GK top-precision margin: carry the full exponent
  with checked width arithmetic and fallible allocation (`MemoryAllocation`), not a cap or silent truncation (B3/P9).
- Fix exact remainder for gaps ≥ 2^63: Mantissa::pow_mod scanned after shifting and dropped the exponent's top bit;
  square-and-multiply now visits every bit (B4/ZF-55).
- Γ(3/4) at p = 2^14 checks a FLINT/Arb 16640-bit reference value (error bound), not just Ok (R9).
- Fold selftest row P checks every commit in the range, case-insensitive against the ghact agent list (+ Vera's GitHub
  commits, U7), with pre-rule imports exempt. Plants cover missing/unknown/empty/extra-word/not-last/repeated trailers;
  this supersedes 11377ab's weaker any-name regex. Row B includes the deliberately removed .gitignore files; row S
  accepts exactly b0afcb1's rustix 1.1.6 RUSTSEC-2026-0486 lock delta (R6/R7).
- CONTRACT 16.CB.4, 16.CB.7 and 16.CB.11 mark `CBall::mul`/`div` and their rows as held (H11), not exported (Spence A1).
- ln 10 sizes its series in 128-bit arithmetic and probes the working allocation, so a valid top precision is
  `MemoryAllocation`, not `InvalidArgument` (B3 follow-on).
- tests/tier1_regressions.rs: 561 judged rows over 24 Tier 1 families (FLINT first, mpmath where FLINT lacks the function,
  exact elsewhere) plus rem, top-precision, past-the-old-trig-cap and 7-mode/signed-zero/Inf/NaN/subnormal/p 32-1024
  boundary tests; red on 5227f77, green here (B2). Test target count 68 -> 69.
- paul-tier1-fixes.md now marks ZF-16, ZF-60 and Γ p ≥ 2048 BROKEN with their remaining row counts (R8).

## Follow-up fixes under Jeff's 2026-11-15 rulings

Status: 30 FIXED, 47 BROKEN of Harper's 77 on 301b21b; no release-readiness claim.
Only quick tests: no gate, mutants, push, PR update, or messages.

- ZF-05 FIXED (`73f2cdb`): decimal parse below the least subnormal rounds by mode.
- ZF-71 FIXED (`eab202b`): huge negative powsi exponents carry scale through reciprocal factors, rather than returning Inf.
- ZF-06 FIXED (`e5b9bc7`): overflowing decimal parse during format round-trip checking means retry, not an Inf spelling of finite MAX.
- ZF-76 FIXED (`d6c571d`): erf/erfc on the imaginary axis have exactly signed-zero / one real parts. FLINT + mpmath references.
- ZF-10 FIXED (`d7e7013`): polynomial GCD uses an exact primitive pseudo-remainder sequence for finite dyadic coefficients. FLINT fmpz_poly references.
- ZF-47 FIXED (`2882165`): Ball::exp radius reaches certified endpoint bounds instead of using a rounded midpoint as an exact multiplier.
- ZF-46 FIXED (`301b21b`): complex log disks meeting the cut cover both sides of the 2πi jump. All 25 Harper disks exercised.

Each fix has an observed red-first in-repository regression and an Agent: paul commit trailer.
For ZF-46 the generic cut test was red first; the additional all-Harper-disks test was added after the fix.
The seven fixes have not undergone the full release gate (per instruction).
