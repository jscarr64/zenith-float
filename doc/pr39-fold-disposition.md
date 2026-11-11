# PR #39 fold disposition (ACCEPTANCE-pr39 L1, L3)

Source: Harper's inventory `zb-fold-inventory.md` (2026-11-10) section B.3 (48 unmerged tips) and E.3. Ancestor
checks were run on the #39 branch with `git merge-base --is-ancestor <tip> HEAD` (log
`~/pop-logs/rusty/2026-11-10/fold-pr1-96-disposition.log`). Counts: (a) 21, (b) 2, (c) 1, and of the 24 rows that were undisposed: 9 DISPOSED (1) with proof
(patch-id twins or merge-parent ancestry, `~/pop-logs/rusty/2026-11-10/fold-pr1-109-l1-proof.log`), 3 FAIL (2) missing objects,
12 FAIL content not shown in history (objects fetched read-only from the zb clone; added-line containment checked).
The FAIL rows are not closed; L1 stays FAIL for them (Paul's ruling 2026-11-10 23:15 CT).

## B.3 unmerged tips

| where | branch | tip | disposition |
|---|---|---|---|
| local | `harper/zb-cosmo-range` | `ca02c8ad` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| local | `harper/zb-goldlib-teamenv` | `d5a19346` | FAIL (2) missing object: not in this repo nor in the zb clone (`git cat-file -t d5a19346`: not a valid object; fold-pr1-109-l1-proof.log); local-only in a checkout agents may not touch. Paul takes it to Jeff. |
| local | `harper/zb-grepq` | `a544b204` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| local | `harper/zb-la-golds` | `58b5eb5c` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| local | `harper/zb-la-rows` | `ed3c6c97` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| local | `harper/zb-la-rows2` | `04e9946f` | (c) held: HOLD H09 (d1d4224). Red tests only (P12), so holding avoids E6. |
| local | `harper/zb-nofailfast` | `3995efdc` | FAIL (2) missing object: not in this repo nor in the zb clone (`git cat-file -t 3995efdc`: not a valid object; fold-pr1-109-l1-proof.log); local-only in a checkout agents may not touch. Paul takes it to Jeff. |
| local | `harper/zb-p03` | `003bddfd` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| local | `harper/zb-popguard-quiet` | `a0a46d32` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| local | `harper/zb-s10-s13-backup` | `2c090d28` | DISPOSED (1) content in fold history: not an ancestor of H, 5 unique commits, `git cherry H 2c090d28` 0 `+` (every commit has a patch-id twin in H); full SHA and merge-base in fold-pr1-109-l1-proof.log. |
| local | `harper/zb-s10-s13-backup2` | `7af6f64a` | DISPOSED (1) content in fold history: not an ancestor of H, 9 unique commits, `git cherry` 0 `+`; full SHA and merge-base in fold-pr1-109-l1-proof.log. |
| local | `harper/zb-s2-f1` | `4a2db745` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| local | `harper/zb-s2-f2` | `8c87b351` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| local | `harper/zb-s2-f3` | `a1f6912f` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| local | `harper/zb-s2-followup` | `fc7b6c0c` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| local | `harper/zb-s7-fuzz-backup` | `d3d7f038` | DISPOSED (1) content in fold history: not an ancestor of H, 8 unique commits, `git cherry` 0 `+`; full SHA and merge-base in fold-pr1-109-l1-proof.log. |
| local | `harper/zb-s7-fuzz-backup2` | `6d221c0a` | DISPOSED (1) content in fold history: not an ancestor of H, 12 unique commits, `git cherry` 0 `+`; full SHA and merge-base in fold-pr1-109-l1-proof.log. |
| local | `harper/zb-tests-pre-rebase` | `13715bb4` | FAIL content not shown in history: not an ancestor of H; 3 `+` commits (7995ffb5, 7aa6974e, 13715bb4): their added lines are not all in zb main or H (e.g. the 10-15 full_gate.sh pubfn-coverage/goldens-regen step lines, 2 rules_scan_selftest lines) (fold-pr1-109-l1-proof.log). Not closed, not guessed. Paul takes it to Jeff. |
| local | `harper/zb-tests-pre-rebase2` | `7bc329eb` | FAIL content not shown in history: not an ancestor of H; 3 `+` commits (a0b6e0f3, 804b06d3, 3c69ed35), same content as 13715bb4's, not all contained (fold-pr1-109-l1-proof.log). Not closed, not guessed. Paul takes it to Jeff. |
| local | `harper/zb-tests-pre-rebase3` | `d71af719` | FAIL content not shown in history: not an ancestor of H; 5 `+` commits (c5106ad2, f5f8add1, f8bcdce0, 70de9f56, 496df5fe), 43 of 39191 added lines of 70de9f56 and 3 of 5 of 496df5fe not found (fold-pr1-109-l1-proof.log). Not closed, not guessed. Paul takes it to Jeff. |
| local | `harper/zb-tests-pre-rebase4` | `21ed9494` | FAIL content not shown in history: not an ancestor of H; 1 `+` commit 121d1d7a (PYTHONNOUSERSITE README text), 3 of 5 added lines not found (fold-pr1-109-l1-proof.log). Not closed, not guessed. Paul takes it to Jeff. |
| local | `paul/zb-pr13-merge` | `6039d721` | DISPOSED (1) content in fold history: merge commit whose parents d327d398 and a8a4bf8d are both ancestors of H; its only own content is the CHANGELOG conflict resolution (union of both sides), and both lines are in H's zenith-float-bound/CHANGELOG.md (fold-pr1-109-l1-proof.log). |
| local | `paul/zb-pr5-merge` | `35fdfff6` | FAIL (2) missing object: not in this repo nor in the zb clone (`git cat-file -t 35fdfff6`: not a valid object; fold-pr1-109-l1-proof.log); local-only in a checkout agents may not touch. Paul takes it to Jeff. |
| local | `rusty/zb-s2-f4-f5` | `3882a6b1` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| local | `spence/zb-l0a-mag-green2-byvalue` | `b2409415` | FAIL content not shown in history: not an ancestor of H; 2 `+` commits (f498db32, b2409415): by-value add_up/sub_lower; 43 of 47 and 6 of 6 added lines not found (fold-pr1-109-l1-proof.log). Not closed, not guessed. Paul takes it to Jeff. |
| origin | `cassie/zb-la2-lu` | `f0ea86a1` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| origin | `cassie/zb-s4b-lambertw` | `b25a2dc8` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| origin | `cloud/zb-l0a-mag-green` | `993cdbe3` | FAIL content not shown in history: not an ancestor of H; 1 `+` commit 402b9862 (cloud agent, 2026-10-16): 27 of 80 added lines not found; P15 forbids taking cloud work (fold-pr1-109-l1-proof.log). Not closed, not guessed. Paul takes it to Jeff. |
| origin | `docs/accumath-family` | `6f517e3a` | (b) ported by commit c1193a0 (names zb 6f517e3, zb PR #1): the Accumath family note in zenith-float-bound/README.md. Closing zb PR #1 is Paul's, with Jeff's OK (U24). |
| origin | `harper/zb-cosmo-range` | `e6e50b8b` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| origin | `harper/zb-la-golds` | `372c30ca` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| origin | `harper/zb-popguard-merge` | `95910716` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| origin | `rusty/zb-la3-qr` | `b3c68b62` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| origin | `rusty/zb-m31-64bit-only` | `cea810ac` | (a) folded by merge: the tip is an ancestor of the #39 head (`git merge-base --is-ancestor`) |
| origin | `saved/2026-10-25/harper-zb-s10-s13-backup` | `2c090d28` | DISPOSED (1) content in fold history: not an ancestor of H, 5 unique commits, `git cherry H 2c090d28` 0 `+` (every commit has a patch-id twin in H); full SHA and merge-base in fold-pr1-109-l1-proof.log. |
| origin | `saved/2026-10-25/harper-zb-s10-s13-backup2` | `7af6f64a` | DISPOSED (1) content in fold history: not an ancestor of H, 9 unique commits, `git cherry` 0 `+`; full SHA and merge-base in fold-pr1-109-l1-proof.log. |
| origin | `saved/2026-10-25/harper-zb-s7-fuzz-backup` | `d3d7f038` | DISPOSED (1) content in fold history: not an ancestor of H, 8 unique commits, `git cherry` 0 `+`; full SHA and merge-base in fold-pr1-109-l1-proof.log. |
| origin | `saved/2026-10-25/harper-zb-s7-fuzz-backup2` | `6d221c0a` | DISPOSED (1) content in fold history: not an ancestor of H, 12 unique commits, `git cherry` 0 `+`; full SHA and merge-base in fold-pr1-109-l1-proof.log. |
| origin | `saved/2026-10-25/harper-zb-tests-pre-rebase` | `13715bb4` | FAIL content not shown in history: not an ancestor of H; 3 `+` commits (7995ffb5, 7aa6974e, 13715bb4): their added lines are not all in zb main or H (e.g. the 10-15 full_gate.sh pubfn-coverage/goldens-regen step lines, 2 rules_scan_selftest lines) (fold-pr1-109-l1-proof.log). Not closed, not guessed. Paul takes it to Jeff. |
| origin | `saved/2026-10-25/harper-zb-tests-pre-rebase2` | `7bc329eb` | FAIL content not shown in history: not an ancestor of H; 3 `+` commits (a0b6e0f3, 804b06d3, 3c69ed35), same content as 13715bb4's, not all contained (fold-pr1-109-l1-proof.log). Not closed, not guessed. Paul takes it to Jeff. |
| origin | `saved/2026-10-25/harper-zb-tests-pre-rebase3` | `d71af719` | FAIL content not shown in history: not an ancestor of H; 5 `+` commits (c5106ad2, f5f8add1, f8bcdce0, 70de9f56, 496df5fe), 43 of 39191 added lines of 70de9f56 and 3 of 5 of 496df5fe not found (fold-pr1-109-l1-proof.log). Not closed, not guessed. Paul takes it to Jeff. |
| origin | `saved/2026-10-25/harper-zb-tests-pre-rebase4` | `21ed9494` | FAIL content not shown in history: not an ancestor of H; 1 `+` commit 121d1d7a (PYTHONNOUSERSITE README text), 3 of 5 added lines not found (fold-pr1-109-l1-proof.log). Not closed, not guessed. Paul takes it to Jeff. |
| origin | `saved/2026-10-25/rusty-zb-s2-contract` | `17e5f1ba` | FAIL content not shown in history: not an ancestor of H; 4 `+` commits (5c494469, c42b9a92, d390222b, 85d165c1): the 10-25 S2 draft; its CONTRACT/CHANGELOG/test text is not verbatim in H (S2 landed reworded) (fold-pr1-109-l1-proof.log). Not closed, not guessed. Paul takes it to Jeff. |
| origin | `saved/2026-10-25/spence-zb-l0a-mag-green2-byvalue` | `b2409415` | FAIL content not shown in history: not an ancestor of H; 2 `+` commits (f498db32, b2409415): by-value add_up/sub_lower; 43 of 47 and 6 of 6 added lines not found (fold-pr1-109-l1-proof.log). Not closed, not guessed. Paul takes it to Jeff. |
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
