- Tests: `num::test_rounding` keeps its subnormal-grid vectors (`n` reaching `|EXPONENT_MIN|`) on 64-bit targets only; on 32-bit, where `n: usize` cannot reach that range, an exact model at `n = usize::MAX` (`Up` gives `2^-n`, `Down`/`ToEven` give zero) replaces them. Test-only.

History: 32-bit targets were dropped by P6; this entry only concerns the `cfg(target_pointer_width = "32")` arms that ZF-M31 removes.
