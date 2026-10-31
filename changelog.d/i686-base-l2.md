- `Mantissa::abs_add` / `abs_sub`: the shifted-operand word count is computed without overflow when the shift is clamped to `usize::MAX` (32-bit).

History: 32-bit targets were dropped by P6; this entry only concerns the `cfg(target_pointer_width = "32")` arms that ZF-M31 removes.
