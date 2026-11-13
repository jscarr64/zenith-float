### Known-wrong grid-digest cells out of the six-mode pins (Spence delta verdict at e1445390)

Test-only. Every grid digest in ext.rs now dumps its cells (`decimal_parse` and `log_exact_results` moved from
`digest` to `dd`, which hashes the same three fields, so their pins do not move). Cells that are known wrong are
counted and left out of the pinned digests, each under the finding it belongs to, and the count is asserted. The
ZF-M27 fixes put them back, red first. The cell ids (dump lines) are listed in the PR body.
- jacobi_all_functions: 62 Down NaN(PrecisionRetryExhausted) cells (F-21), 24 Down +overflow cells (+Inf, the ruling is
  +MAX) and 16 dn/nd cells with a wrong exact flag (F-22 rider), all at the tiny operand 21.
- trig_hyperbolic: 6 None tanh cells that are exactly +-1 (F-16).
- reciprocal_nth_root: reciprocal(third(320)) at p = 200 under None, exactly 3 (F-17), and 12 None overflow cells
  (NONE-GAP F-B).
- ziv_entry_points (2) and log_exact_results (10): exact log results flagged inexact in ToEven (F-20).
- complex_mul_div: 2970 NaN quotients by a divisor with a subnormal part (F-18), 36 finite quotients given as +-Inf,
  and 1737 toward-zero overflow cells (+-Inf, the ruling is +-MAX).
The six-mode pins otherwise do not move, and no ToOdd pin moves.
