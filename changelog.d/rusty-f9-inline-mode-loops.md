### Inline mode loops: all 7 modes, or labelled mode-specific

- The loops whose expected value holds in every mode now run all 7 modes: complex div of an exact quotient, complex
  erf on the real axis, round/sqrt of specials, round's own flag, pow near the least normal, fuzz_round_modes'
  fixed MPFR mul/div case, parse/format round trips, and the top-precision entry checks.
- The loops whose expected value holds only for their modes carry `MODE-SUBSET mode-specific: <why>`: the pinned
  digests (Q2), modf's and the integer getters' split by direction, the subnormal-to-integer zero, and the reciprocal
  -Inf rows. The integer getters gain their ToOdd and None rows (tiny values; -0.5 under None).
  `tools/test_rounding_mode_matrices.sh` checks the label.
