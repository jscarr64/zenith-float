### zenith-float-bound LA heavy pin: la_logdet_heavy.txt at format 2

- `tests/goldens/slices/cholesky/HEAVY.SHA256` still pinned `la_logdet_heavy.txt` at its format-1 bytes
  (5386579, 2026-10-24, Cassie's run 1541-la4-heavy3); d9cfacb (2026-10-29) added the format-2 cond columns
  and b4dee79 re-pinned only the gate file. Re-pinned to the regenerated file (11 rows, format 2). Removing
  the format-2 header lines and the 6 cond columns from the new file gives the old pinned bytes exactly.
