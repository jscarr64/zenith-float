- `compute_added_err` (Pow): the exponent bound is compared in `u64`, so it no longer truncates or overflows on 32-bit.

History: 32-bit targets were dropped by P6; this entry only concerns the `cfg(target_pointer_width = "32")` arms that ZF-M31 removes.
