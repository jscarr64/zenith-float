### BORROW rule: closes the copy holes in the add_ref/sub_ref/mul_ref/div_ref call-graph gate

- `tools/rules_scan.py` BORROW rule now scans the whole `zenith-float-num/src` (non-test), strips nested block
  comments, strings and char literals before matching, and resolves calls by impl type, field/local types and arity.
  It follows UFCS/path calls, fns passed as values, operators into op-trait impls and `macro_rules!` bodies, and bans
  `Clone::clone(`, `.clone ()` with spaces, `<[T]>::to_vec(`, `Vec::from(`/`Box::from(`, `|x| *x` deref closures,
  `.extend_from_slice(` and `.clone_from(`. The failure message names the op and the reach chain.
- `tools/rules_scan_selftest.sh`: 12 new BORROW fixtures (11 holes and the reach-chain message, plus a multi-line
  `/* */` comment false-alarm case).
- Library, no behaviour change: `Mantissa::abs_cmp`/`cmp_iters` and `NormalizedView` read words through borrows
  (no deref-copy closure), and the `WordBuf` inline-to-heap spill moves the inline words with `swap_with_slice`
  instead of copying them. `borrow_ops_golds` and every `*_grid_digest` are unchanged (not regenerated).
