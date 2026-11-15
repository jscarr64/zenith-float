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
