### zenith-float-bound S6 heavy-tier pin follows the crate rename

- `tests/goldens/slices/matrix/HEAVY.SHA256` pinned `matrix_rows_heavy.txt` as generated before the
  fold renamed the crate; the file's first line names the crate, so a fresh
  `gen_matrix_goldens.py --heavy DIR` never matched. The pin is now the regenerated file's sha256
  (same rows; only the header line differs). All 30 heavy rows replay through the crate.
- `tests/goldens/slices/lu/HEAVY.SHA256`: `la_solve_heavy.txt` was still pinned at its format-1
  bytes (zb 3f7f15e); format 2 (cond columns) reached the gate files in b4dee79 without a heavy
  re-pin. Re-pinned to the team-env regeneration (6 rows, format 2).
