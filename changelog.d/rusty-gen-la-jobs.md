### gen_la_goldens: `--jobs N` computes the logdet cond columns on N cores (byte-identical)

- `scripts/la_cond_cols.py` gains `par_map` and a `jobs=1` argument on `la_cond_cols` / `_inverse`; `lu_pp(A, jobs=1)` in
  `scripts/gen_la_goldens.py` checks `L U == P A` per row. With `jobs > 1`, forked workers compute the |L||U| rows, the
  A^-1 columns and the |A^-1| G rows (exact gmpy2.mpq), and the results are reassembled in index order, so the bytes do
  not depend on N. `--jobs N` applies it to the logdet cond columns; the default 1 is the serial path. The elimination in
  `lu_pp` stays serial (about 4-7% of the `main` time in the fold-pr1-119 profile).
- `tools/gen_la_jobs_selftest.sh` (11 cases) and `tools/gen_la_rowcache_selftest.sh` now run in full_gate's gen-la-heavy step.
