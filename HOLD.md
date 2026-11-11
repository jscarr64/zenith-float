# HOLD index

The one index of held code in zenith-float (Jeff, 2026-11-10 15:36, HOLD method). Held code was folded in from the
bound-layer history but is not live: each held Rust file begins with a `// HOLD <id>:` line and every original line
follows commented out with a leading `// `. No `mod` line, `#[path]` or Cargo target names a held file (its
`mod` line stays in place, commented as `// HOLD <id>: ...`), so held code is not compiled, and held tests do
not run. Held code and held tests are never counted as passing or as ported. `tools/fold_pr1_selftest.sh` rows HOLD
and HOLDIDX check this index against the headers both ways.

Rules (Jeff, 2026-11-10 15:37 addendum): before writing code, check this file first; if a held block fits the work,
take it live, with its tests, instead of writing new code. When zenith-float 2.0 is 100% finished, every block still
held is deleted, and this file with it; that closing step is the LAST item of RELEASE-2.0-CHECKLIST.md (6.16).

To take a block live: uncomment the file whole (drop the two header lines and the leading `// `), restore its
`mod` line, and run the crate tests. To delete one: remove its files, its `mod` line and its row here, in one commit.

Data that belongs to a held slice stays where it was and is listed in the last column: goldens under
`zenith-float-bound/tests/goldens/slices/<slice>/`, mutants entries in `zenith-float-bound/tools/mutants/{manifest,floor}.d/<slice>.txt`
and test floors in `zenith-float-bound/tools/test_floor.d/<slice>.txt`. It is data, not code, and the zf gate does not read the
bound crate's own floors or mutants yet (fold plan M-08), so a held slice's floor is not enforced anywhere. Its
CONTRACT text in `zenith-float-bound/docs/` is unchanged.

| id | path(s) | origin | why held | decides use or delete | held data |
|----|---------|--------|----------|-----------------------|-----------|
| H01 | `zenith-float-bound/src/elementary.rs`, `zenith-float-bound/tests/harness/slice_elementary.rs` | zb 9459b57 (S0 stub, on zb main fcaf135) | S3 stub, doc comments only; no 2.0 function is implemented | Jeff's ruling on the RELEASE-2.0-CHECKLIST 2.1 signature table and 3.x L1 scope: implemented in 2.0, go live, else delete | goldens slices/elementary (GENERATORS, SHA256SUMS only); manifest/floor/test_floor elementary.txt |
| H02 | `zenith-float-bound/src/hypergeom.rs`, `zenith-float-bound/tests/harness/slice_hypergeom.rs` | zb 9459b57 (S0 stub, on zb main fcaf135) | S4a stub, doc comments only | Jeff's ruling on checklist 3.x L2 scope (hypergeometric family in or out of 2.0) | goldens slices/hypergeom (GENERATORS, SHA256SUMS only); manifest/floor/test_floor hypergeom.txt |
| H03 | `zenith-float-bound/src/zf_special.rs`, `zenith-float-bound/tests/harness/slice_zf_special.rs` | zb 9459b57 (S0 stub, on zb main fcaf135) | S4c stub, doc comments only | Jeff's ruling on checklist 3.x L2 scope (zf special-function wrappers) | manifest/floor/test_floor zf_special.txt |
| H04 | `zenith-float-bound/src/quadrature.rs`, `zenith-float-bound/tests/harness/slice_quadrature.rs` | zb 9459b57 (S0 stub, on zb main fcaf135) | S5 stub, doc comments only | Jeff's ruling on checklist 3.x L3 scope (quadrature in or out of 2.0) | manifest/floor/test_floor quadrature.txt |
| H05 | `zenith-float-bound/src/carlson.rs`, `zenith-float-bound/tests/harness/slice_carlson.rs` | zb cassie/zb-s4b-lambertw b25a2dc8 (Lambert W W0/W-1, 889a2c6, 3e81e10), merged at 26e8454 | Implemented and tested on the tip, but Harper's inventory "needs" column says yes only if S4b is still in 2.0 L2; unsure means HOLD | Jeff's ruling on whether S4b (Lambert W, Carlson) is in 2.0; if in, go live and its tests count from then | goldens slices/carlson/lambertw_rows.txt (rehashed for the rename); manifest/floor/test_floor carlson.txt |
| H06 | `zenith-float-bound/src/cholesky.rs`, `zenith-float-bound/tests/harness/slice_cholesky.rs` | zb 4f25b9a (LA API skeleton, via rusty/zb-la3-qr b3c68b62), merged at 0d769b3 | LA1 signatures only: every body returns BoundError::Resource "not_implemented"; harness rows are gold-shape checks for that unimplemented slice | An LA1 implementation landing (Cholesky, LDL, is_spd, logdet_spd) with Jeff keeping LA1 in 2.0; else delete | goldens slices/cholesky (la_cholesky.txt; la_logdet.txt is also read by the live LU slice); manifest/floor/test_floor cholesky.txt |
| H07 | `zenith-float-bound/src/gls.rs`, `zenith-float-bound/tests/harness/slice_gls.rs` | zb 4f25b9a (LA API skeleton, via rusty/zb-la3-qr b3c68b62), merged at 0d769b3 | LA4 signature only: gls() returns Resource "not_implemented"; harness rows are gold-shape checks | An LA4 implementation (needs LA1) with Jeff keeping LA4 in 2.0; else delete | goldens slices/gls; manifest/floor/test_floor gls.txt |
| H08 | `zenith-float-bound/src/eigen.rs`, `zenith-float-bound/tests/harness/slice_eigen.rs` | zb 4f25b9a (LA API skeleton, via rusty/zb-la3-qr b3c68b62), merged at 0d769b3 | LA5 signatures only: eig_sym and eig return Resource "not_implemented" | An LA5 implementation with Jeff keeping LA5 in 2.0; else delete (EigOpts in the live policy.rs goes with it) | goldens slices/eigen; manifest/floor/test_floor eigen.txt |
| H09 | `zenith-float-bound/tests/harness/slice_lu/harper_la_fact.rs` | zb harper/zb-la-rows2 04e9946f (Harper, 2026-10-29: LA own-id reds 16.LA2.1, 16.LA2.2, 16.LA.1, 16.LA.6; found in the local clone /home/jsc64/harper-work/zb-la-rows), not merged; file taken with the rename sed, CHANGELOG hunk not taken | Ignored reds only. The file imports `zenith_float_bound::gls::gls`, and gls is held (H07), so it cannot compile live. Its LU rows would add 4 ignored reds to a live slice | H07 going live (LA4 kept in 2.0): then take it live with its LU reds; if H07 is deleted, drop the gls row and take the rest live, or delete | none (no goldens, no mutants rows) |
