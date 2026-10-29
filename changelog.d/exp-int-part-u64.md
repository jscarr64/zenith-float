- `exp` of an argument whose integer part needs more than 32 bits (below `2^64`) no longer reports
  `ExponentOverflow` on 32-bit targets: the integer part is read as a `u64` (`int_as_u64`) instead of a
  `usize`, and the range guard is `2^64` instead of `2^WORD_BIT_SIZE`. 64-bit results are unchanged. This
  fixes `pow_near_min_normal_does_not_underflow`, `test_acosh` and `test_asinh` on i686.
