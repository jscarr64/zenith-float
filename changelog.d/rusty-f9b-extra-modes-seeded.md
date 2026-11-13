### Seeded ToOdd and None passes for the MPFR call-site families

`get_random_rnd_pair` draws only the five MPFR modes, and its `% 5` streams are unchanged. The new
`tests/mpfr/extra_modes_seeded.rs` adds one seeded pass per call-site family for ToOdd and None through
`extra_mode_reference`, with its own tests: exact, inexact, lsb already 1, lsb already 0, and the binade edge.
- None is MPFR toward zero at rp(p, 64), with the inexact flag.
- ToOdd is nearest with ties to odd, as defs.rs defines it.

Families covered: ops, special_fn (positive arguments; `ln_gamma` is documented for positive arguments only), const,
special values, and complex add. The complex flag is checked as value-wide. copysign is left out because of new
finding N-4: at a precision that is not a word multiple, it rounds at rp(p, 64) in every mode.

Each extra-mode loop is an inline `[ToOdd, None]` loop with a `MODE-SUBSET mode-specific` label, not a named 2-mode
array. The structural check refuses named arrays that are not all 7 modes (fold-pr1-174).
