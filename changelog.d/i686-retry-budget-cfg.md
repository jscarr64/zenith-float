- `retry_budget_tests`: 32-bit expected values for the fixed retry budget (`WORD_BIT_SIZE * MAX_PREC_RETRY` = 8192 there); 64-bit values unchanged.

History: 32-bit targets were dropped by P6; this entry only concerns the `cfg(target_pointer_width = "32")` arms that ZF-M31 removes.
