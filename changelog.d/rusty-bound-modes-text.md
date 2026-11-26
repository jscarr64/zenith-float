### zenith-float-bound mode text: six valid modes, None not applicable

- `scripts/README.md` now describes `ball_exact_ops.txt` as exactly rounded add, sub, mul, div and sqrt over the 6 valid modes,
  and says why `RoundingMode::None` does not apply to bound (no public bound fn takes a mode). The stale `ModeOrNone` comment in
  `tests/common/generators.rs` no longer claims a test asserts that an op rejects `None`. Text only; no golden or pin changes.
