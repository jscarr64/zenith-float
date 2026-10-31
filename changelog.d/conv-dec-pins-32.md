- Tests: `test_conv_num` decimal pins for 32-bit targets, computed independently with mpmath at the 32-bit-rounded precisions; 64-bit pins unchanged.

History: 32-bit targets were dropped by P6; this entry only concerns the `cfg(target_pointer_width = "32")` arms that ZF-M31 removes.
