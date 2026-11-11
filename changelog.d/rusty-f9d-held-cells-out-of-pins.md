### Known-wrong grid-digest cells out of the six-mode pins (Spence delta verdict at e1445390)

Test-only. Every grid digest in ext.rs now dumps its cells (`decimal_parse` and `log_exact_results` moved from
`digest` to `dd`, which hashes the same three fields, so their pins do not move). Cells that are known wrong are
counted and left out of the pinned digests, each under the finding it belongs to, and the count is asserted. The
ZF-M27 fixes put them back, red first. The cell ids (dump lines) are listed in the PR body.
