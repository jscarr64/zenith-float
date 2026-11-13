### zenith-float-compare refuses 32-bit targets

- `zenith-float-compare/src/main.rs` now has the same 64-bit `compile_error!` as the facade, num and bound
  crate roots (P6/P8). `tools/test_full_gate_64bit_only.sh` checks it and its removal probe.
