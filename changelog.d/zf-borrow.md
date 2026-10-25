### Borrowing arithmetic: `add_ref`/`sub_ref`/`mul_ref`/`div_ref` are sanctioned, copy-free and gated

- `ExactNum::{add_ref, sub_ref, mul_ref, div_ref}(a: &ExactNum, b: &ExactNum, p, rm) -> ExactNum`
  keep their names and signatures and are now documented as the sanctioned borrowing ops (no
  longer a last resort): rounded once, both operands read in place (aliasing `(&x, &x)` allowed),
  only the result and the kernels' own working space allocated, bit-identical to the by-value
  `add`/`sub`/`mul`/`div`. Errors stay NaN values, as for every `ExactNum` op.
- The four `.copied()` calls on their path (`RightShiftedSlice` in `mantissa/util.rs`,
  `write_shifted_left` in `mantissa/base.rs`) now read the word in place. No behaviour change.
- `tools/rules_scan.sh` replaces the `add_ref`/`sub_ref` and `mul_ref`/`div_ref` use-count
  ratchets (and `tools/ref_sum_count_baseline.txt`, `tools/mul_ref_count_baseline.txt`) with the
  BORROW rule: any `.clone(`/`.to_vec(`/`.cloned(`/`.copied(`/`.to_owned(` or helper-copy call
  in a library fn reachable from the four ops fails. `--borrow-reach` lists the reached fns.
- New goldens `tests/oracle/goldens/borrow_ops.jsonl` (3950 rows, mpmath at p+256 and p+512
  agreeing, from `scripts/gen_borrow_goldens.py`) and `tests/borrow_ops_golds.rs`: every rounding
  mode, both signs, zeros, subnormals, underflow, near-`EXPONENT_MAX` and overflow, aliasing,
  identity with the by-value op, and result storage of exactly `ceil(p/64)` words.
