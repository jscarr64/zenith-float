- `series_niter`: an iteration-count denominator that exceeds `usize` gives 0 iterations instead of overflowing (32-bit).

History: 32-bit targets were dropped by P6; this entry only concerns the `cfg(target_pointer_width = "32")` arms that ZF-M31 removes.
