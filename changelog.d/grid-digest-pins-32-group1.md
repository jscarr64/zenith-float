- Grid digest tests `operator_grid`, `series_high_precision`, `fma_stand_in_branches`, `hypot_scaling` and `powi_windowed` now carry 32-bit pins generated from mpmath at the `Word = u32` layout (64-bit pins unchanged and reproduced by the generator's control run).

History: 32-bit targets were dropped by P6; this entry only concerns the `cfg(target_pointer_width = "32")` arms that ZF-M31 removes.
