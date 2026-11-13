### Grid digests: ToOdd joins with digests of its own

- The 16 grid loops in `ext.rs` `grid_digest_tests` now run all 7 modes. The six modes pinned before keep their old
  digests, hashed in their old order. ToOdd is hashed into a new digest of its own (Q2) in 14 pinned tests.
- ToOdd rounds to nearest with ties to odd, so like ToEven it overflows to ±Inf, and its overflow cells are correct
  to pin. Every ToOdd cell was checked against `scripts/verify_grid_digests.py`. The two verifier BAD/FLAG rows are
  verifier artifacts and were checked by hand: `rnd_fixed` rounds a ToOdd tie as ToEven, and x² at p=200 is exact.
- complex_mul_div holds ToOdd back: it would pin new finding N-3, where a complex quotient with a subnormal
  denominator gives NaN in every mode.
- Fixed the text: "7 precisions in all 6 rounding modes" was wrong (the square digests use 9 precisions).
