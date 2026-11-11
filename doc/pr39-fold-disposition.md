# PR #39 fold disposition (ACCEPTANCE-pr39 L1, L3)

Source: Harper's inventory `zb-fold-inventory.md` (2026-11-10) section B.3 (48 unmerged tips) and E.3. Ancestor
checks were run on the #39 branch with `git merge-base --is-ancestor <tip> HEAD` (log
`~/pop-logs/rusty/2026-11-10/fold-pr1-96-disposition.log`). Counts: (a) 21, (b) 2, (c) 1,
undisposed 24. The undisposed rows are listed as such; they are not passed silently (L1 and L2 stay FAIL/OPEN
for them).

## B.3 unmerged tips

| where | branch | tip | disposition |
|---|---|---|---|
| local | `harper/zb-cosmo-range` | `ca02c8ad` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| local | `harper/zb-goldlib-teamenv` | `d5a19346` | UNDISPOSED, blocked on access: local-only in a checkout agents may not touch (Jeff's MathWorks checkouts). Needs a Jeff or Paul export, then (b) or (e). |
| local | `harper/zb-grepq` | `a544b204` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| local | `harper/zb-la-golds` | `58b5eb5c` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| local | `harper/zb-la-rows` | `ed3c6c97` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| local | `harper/zb-la-rows2` | `04e9946f` | (c) held: HOLD H09 (d1d4224). Red tests only (P12), so holding avoids E6. |
| local | `harper/zb-nofailfast` | `3995efdc` | UNDISPOSED, blocked on access: as d5a19346. |
| local | `harper/zb-p03` | `003bddfd` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| local | `harper/zb-popguard-quiet` | `a0a46d32` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| local | `harper/zb-s10-s13-backup` | `2c090d28` | UNDISPOSED: historical backup tip (2026-10-15/16, inventory: 'ahead>0 only because the tip is an old line'). The object is not in this repo; the 'nothing unique' diff (E.3) is not done. |
| local | `harper/zb-s10-s13-backup2` | `7af6f64a` | UNDISPOSED: historical backup tip (2026-10-15/16, inventory: 'ahead>0 only because the tip is an old line'). The object is not in this repo; the 'nothing unique' diff (E.3) is not done. |
| local | `harper/zb-s2-f1` | `4a2db745` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| local | `harper/zb-s2-f2` | `8c87b351` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| local | `harper/zb-s2-f3` | `a1f6912f` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| local | `harper/zb-s2-followup` | `fc7b6c0c` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| local | `harper/zb-s7-fuzz-backup` | `d3d7f038` | UNDISPOSED: historical backup tip (2026-10-15/16, inventory: 'ahead>0 only because the tip is an old line'). The object is not in this repo; the 'nothing unique' diff (E.3) is not done. |
| local | `harper/zb-s7-fuzz-backup2` | `6d221c0a` | UNDISPOSED: historical backup tip (2026-10-15/16, inventory: 'ahead>0 only because the tip is an old line'). The object is not in this repo; the 'nothing unique' diff (E.3) is not done. |
| local | `harper/zb-tests-pre-rebase` | `13715bb4` | UNDISPOSED: historical backup tip (2026-10-15/16, inventory: 'ahead>0 only because the tip is an old line'). The object is not in this repo; the 'nothing unique' diff (E.3) is not done. |
| local | `harper/zb-tests-pre-rebase2` | `7bc329eb` | UNDISPOSED: historical backup tip (2026-10-15/16, inventory: 'ahead>0 only because the tip is an old line'). The object is not in this repo; the 'nothing unique' diff (E.3) is not done. |
| local | `harper/zb-tests-pre-rebase3` | `d71af719` | UNDISPOSED: historical backup tip (2026-10-15/16, inventory: 'ahead>0 only because the tip is an old line'). The object is not in this repo; the 'nothing unique' diff (E.3) is not done. |
| local | `harper/zb-tests-pre-rebase4` | `21ed9494` | UNDISPOSED: historical backup tip (2026-10-15/16, inventory: 'ahead>0 only because the tip is an old line'). The object is not in this repo; the 'nothing unique' diff (E.3) is not done. |
| local | `paul/zb-pr13-merge` | `6039d721` | UNDISPOSED: local-only merge-topology commit (1 merge, no own content per the inventory); the object is not in this repo, so the 'nothing unique' diff is not done. |
| local | `paul/zb-pr5-merge` | `35fdfff6` | UNDISPOSED: as 6039d721. |
| local | `rusty/zb-s2-f4-f5` | `3882a6b1` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| local | `spence/zb-l0a-mag-green2-byvalue` | `b2409415` | UNDISPOSED: historical backup tip (2026-10-15/16, inventory: 'ahead>0 only because the tip is an old line'). The object is not in this repo; the 'nothing unique' diff (E.3) is not done. |
| origin | `cassie/zb-la2-lu` | `f0ea86a1` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| origin | `cassie/zb-s4b-lambertw` | `b25a2dc8` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| origin | `cloud/zb-l0a-mag-green` | `993cdbe3` | UNDISPOSED, OPEN: cloud-agent commit (P15 forbids new cloud work). Needs the hunks re-authored by the owner with tests, or a Jeff drop. Not cherry-picked. |
| origin | `docs/accumath-family` | `6f517e3a` | (b) ported by commit c1193a0 (names zb 6f517e3, zb PR #1): the Accumath family note in zenith-float-bound/README.md. Closing zb PR #1 is Paul's, with Jeff's OK (U24). |
| origin | `harper/zb-cosmo-range` | `e6e50b8b` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| origin | `harper/zb-la-golds` | `372c30ca` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| origin | `harper/zb-popguard-merge` | `95910716` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| origin | `rusty/zb-la3-qr` | `b3c68b62` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| origin | `rusty/zb-m31-64bit-only` | `cea810ac` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| origin | `saved/2026-10-25/harper-zb-s10-s13-backup` | `2c090d28` | UNDISPOSED: historical backup tip (2026-10-15/16, inventory: 'ahead>0 only because the tip is an old line'). The object is not in this repo; the 'nothing unique' diff (E.3) is not done. |
| origin | `saved/2026-10-25/harper-zb-s10-s13-backup2` | `7af6f64a` | UNDISPOSED: historical backup tip (2026-10-15/16, inventory: 'ahead>0 only because the tip is an old line'). The object is not in this repo; the 'nothing unique' diff (E.3) is not done. |
| origin | `saved/2026-10-25/harper-zb-s7-fuzz-backup` | `d3d7f038` | UNDISPOSED: historical backup tip (2026-10-15/16, inventory: 'ahead>0 only because the tip is an old line'). The object is not in this repo; the 'nothing unique' diff (E.3) is not done. |
| origin | `saved/2026-10-25/harper-zb-s7-fuzz-backup2` | `6d221c0a` | UNDISPOSED: historical backup tip (2026-10-15/16, inventory: 'ahead>0 only because the tip is an old line'). The object is not in this repo; the 'nothing unique' diff (E.3) is not done. |
| origin | `saved/2026-10-25/harper-zb-tests-pre-rebase` | `13715bb4` | UNDISPOSED: historical backup tip (2026-10-15/16, inventory: 'ahead>0 only because the tip is an old line'). The object is not in this repo; the 'nothing unique' diff (E.3) is not done. |
| origin | `saved/2026-10-25/harper-zb-tests-pre-rebase2` | `7bc329eb` | UNDISPOSED: historical backup tip (2026-10-15/16, inventory: 'ahead>0 only because the tip is an old line'). The object is not in this repo; the 'nothing unique' diff (E.3) is not done. |
| origin | `saved/2026-10-25/harper-zb-tests-pre-rebase3` | `d71af719` | UNDISPOSED: historical backup tip (2026-10-15/16, inventory: 'ahead>0 only because the tip is an old line'). The object is not in this repo; the 'nothing unique' diff (E.3) is not done. |
| origin | `saved/2026-10-25/harper-zb-tests-pre-rebase4` | `21ed9494` | UNDISPOSED: historical backup tip (2026-10-15/16, inventory: 'ahead>0 only because the tip is an old line'). The object is not in this repo; the 'nothing unique' diff (E.3) is not done. |
| origin | `saved/2026-10-25/rusty-zb-s2-contract` | `17e5f1ba` | UNDISPOSED: historical backup tip (2026-10-15/16, inventory: 'ahead>0 only because the tip is an old line'). The object is not in this repo; the 'nothing unique' diff (E.3) is not done. |
| origin | `saved/2026-10-25/spence-zb-l0a-mag-green2-byvalue` | `b2409415` | UNDISPOSED: historical backup tip (2026-10-15/16, inventory: 'ahead>0 only because the tip is an old line'). The object is not in this repo; the 'nothing unique' diff (E.3) is not done. |
| origin | `spence/zb-la-api` | `5f30568b` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| origin | `spence/zb-la-contract` | `ed3c6c97` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| origin | `spence/zb-la0-trisolve` | `e68a221d` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| origin | `spence/zb-s22-overflow` | `13f1122a` | (b) ported by commit 7163066 (names zb 13f1122a): cmp_exact_diff above EXPONENT_MAX, with its tests and mutants 357-362. Unique-hunk diff vs its base 3882a6b1: OPEN for the reviewer (L2). |

## L3: zb main blobs not at the same path under zenith-float-bound/

All 715 blobs of zb main `fcaf135` are present under `zenith-float-bound/` except these, each moved or dropped:

| zb path | disposition |
|---|---|
| `.github/CODEOWNERS` | moved: root `.github/CODEOWNERS`, entries prefixed `zenith-float-bound/` (checklist 4.5) |
| `.gitignore` | moved: root `.gitignore` (`**/target`, `__pycache__/`, `.venv/`; Cargo.lock is the workspace one), 0c23d9f |
| `rust-toolchain.toml` | dropped: the workspace root toolchain applies (one library) |
| `tests/goldens/slices/.gitignore` | moved: root `.gitignore` rule `zenith-float-bound/tests/goldens/slices/**/*_heavy.txt`, 3aa7292 |
| `tools/mutants/117-P-min2.patch` .. `123-P-cmp.patch` (7 files) | dropped in d15bb37: they mutated the bound-local MAX_PRECISION formula, which is gone (one zf definition, re-exported); fold selftest row B lists them as removed |

## E.3 checklist

| E.3 bullet | disposition |
|---|---|
| every zb main blob present | (a) see L3 above |
| public API reachable as `zenith_float::bound::...` | (a) `pub use zenith_float_bound as bound` behind the default-on `bound` feature; facade test `tests/fold_pr1_facade.rs`. Per-item A.3 loop: OPEN (L4) |
| stub slices present; unmerged implementation tips folded | (a)/(c): see B.3; stubs held H01-H08 |
| goldens byte-identical, `-text`, HEAVY.SHA256 pins | (a) golden-sums gate step; heavy regen check: `fold-pr1-83-heavy-gen.log` |
| 374 mutant patches apply | (a) except the 7 retired precision patches (L3); apply-check is part of the edgegate/mutants work (I7) |
| gate tools run from the zf toplevel (M-08) | (c) held: HOLD H10 |
| LICENSE/NOTICE files in the tree | (a) present under zenith-float-bound/ (B6 dedupe: Needs Jeff) |
| CONTRACT, TEST_PLAN, PLAN in the tree | (a) zenith-float-bound/docs/ |
| P6 implemented | (a) compile_error! in all three crate roots (3b9d543), one MAX_PRECISION (d15bb37) |
| zf num `ball.rs` dual Ball | OPEN: C23 / D8, for Jeff |
| no Accumath/latex-rust/hubble source needed | (a) none depends on the crate |
| FIXLIST M-01..M-13 | (d) later: M-08 held (H10); M-13 waits for Jeff |
| zb PR #1 disposed; PR #12 dc51e37 and merge-only commits diffed | PR #1: (b) c1193a0, close is Paul's with Jeff's OK; dc51e37 and the merge-only commits: OPEN |
| historical saved/backup tips diffed | OPEN: see the undisposed B.3 rows |
| tags | none |
| worktrees and target caches not copied | (a) none in the tree |
| zb repo and checkouts not deleted | (a) L5: repo not archived/disabled, 31 branches |
