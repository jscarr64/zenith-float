### Tests (red): `Consts::get_exact_prec`

- `tests/r14_consts.rs` with 1050 mpmath golds (`scripts/gen_r14_const_goldens.py`,
  `tests/oracle/goldens/r14_consts.jsonl`): π, e, ln 2, ln 10, √2, φ and γ at p in
  1..=65536 (1, 2, 63, 64, 65 and 65536 included) in six modes, settled at p+256 and p+512
  bits; `Down < Up` bracketing each constant at every gold p; results bit-identical for
  1000-then-64 and 64-then-1000 on fresh caches and on a pre-warmed cache;
  `Consts::is_exact` false for all seven; p = 0 and `RoundingMode::None` rejected.
- Red: `ConstantId`, `Consts::get_exact_prec` and `Consts::is_exact` do not exist yet.
