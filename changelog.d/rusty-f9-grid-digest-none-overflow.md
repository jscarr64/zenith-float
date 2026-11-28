### Grid digests: None overflow cells out of the pins (NONE-GAP F-B)

- Five grid digests in `ext.rs` no longer hash the `None` overflow cells: ldexp/scalb (378 cells), mul/div and
  mul_ref/div_ref (448 each), sqr and square via mul_ref (36 each). Under `None` those cells give ±Inf, and the ruling
  is ±MAX, so the old pins held known-wrong bits. The skipped cells are counted, and the count is asserted.
- The five digests are re-pinned over the remaining cells. A cell-by-cell dump diff shows every remaining cell is
  bit-identical and only the listed cells are gone.
- Old and new digests:
  - ldexp: 1196293272222150608 → 14437472841843974673
  - mul/div: 8827555224262690728 → 6854509046823660977
  - sqr: 17269539713433390012 → 6265526590754367528
- The ZF-M27 fix puts the cells back and re-pins again, with proof.
