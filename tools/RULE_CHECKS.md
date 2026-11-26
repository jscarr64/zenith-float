# Rule checks: drift, unpushed cap, pre-push, PR done

Plan v2 step 1 (Jeff 2026-11-25: "push and merge as we go", "make sure PRs are closed once merged",
"limit the check to active work trees"). Piled-up local work and forgotten PRs fail loudly instead of growing.

| Check | File | FAILS when | rc |
|---|---|---|---|
| Drift | `tools/drift_check.py` | a local branch or detached HEAD is **more than 1** commit ahead of what is pushed; the oldest unpushed commit is **older than 30 min**; the oldest dirty path (modified, deleted, untracked, not ignored, submodules included) is **older than 30 min**. Exactly 30 min passes. | 0 ok, 1 drift, 2 usage or git error |
| Unpushed cap | `.githooks/pre-commit` -> `drift_check.py --pre-commit` | the current branch (or detached HEAD) already has **1** unpushed commit: push it first. An unborn branch has 0. | 0 / 1 / 2 (not in a work tree, git error) |
| Pre-push | `.githooks/pre-push` | for each pushed ref, in a clean `--shared` clone checked out at the **pushed sha** (tags peeled): `tools/rules_scan.sh "${RATCHET_BASE:-c74e834^}..<sha>"` fails or is missing; `tools/no_suppress_scan.py` (when the pushed tree has it) fails; the clone cannot be made. No tests run, so a red-first commit pushes. | 0 / 1 |
| PR done | `tools/pr_done_check.py` | an **OPEN** PR whose head is already on its base (merged-open); a **landed** PR (MERGED, or CLOSED with head on base) whose head branch still exists (branch-not-deleted). Cross-repo PRs skipped. | 0 ok, 1 not done, 2 any gh error, or a list that reaches `--limit` |

**"Ahead of what is pushed."** Commits on the branch that are not on its upstream, when the upstream is an existing
remote-tracking ref (`refs/remotes/*`). With no upstream, an upstream that is gone, or an upstream that is another
**local** branch, it is the commits on no remote-tracking ref.

**Active mode (`--scan`).** Only active work trees are checked. A tree is active when:
- its newest change is at most **24 h** old (`--active-window-min 1440`; exactly 24 h is active). The changes
  counted are the HEAD commit time, the index mtime, the `logs/HEAD` mtime and the newest dirty file; or
- it has dirty work **under 30 min** old, whatever the window; or
- it has dirty work **over 30 min** old (that is the violation itself, so it is never skipped as idle), whatever the
  window (Cassie R2-B2d); or
- it is listed with `--include DIR` or in `--include-file FILE` (one path per line, `#` comments); or
- its activity cannot be read (fail closed: the tree is checked, and a git error then gives rc 2).

In an active scanned tree the HEAD branch is checked, plus other local branches whose tip is at most 24 h old.
Idle trees and branches print `skip`. `--all` is the full sweep; `--repo` and `--include` trees are always checked
in full. `--scan DIR` checks DIR itself and its direct children. A child that has a `.git` but cannot be read
(missing gitdir, corrupt HEAD), or whose `.git` cannot even be looked up (no permission, e.g. mode 000), is named and
gives rc 2; a child without `.git` (a plain dir, a bare repo) is skipped. An unreadable DIR gives rc 2. Linked worktrees share branches, so each branch
of a repository is reported once.

Options: `drift_check.py --repo DIR | --scan DIR [--all] [--include DIR] [--include-file FILE]
[--active-window-min 1440] [--max-ahead 1] [--max-age-min 30] [--exclude NAME] [--pre-commit] [--now EPOCH]`;
`--exclude NAME` (repeatable) leaves out, and names as `excluded`, every `--scan` entry (DIR or child) with a path
component NAME, as written or after following its symlinks; `--repo` / `--include` trees are never excluded.
Nothing named NAME, and nothing under it, is ever stat'ed, listed or read (Cassie round 4):
- A `--scan` DIR is judged before any isdir/listdir/stat on it.
- Symlinks are followed by name only: the path is walked component by component, each name is checked before it
  is touched, and a symlink's `readlink` string is spliced in and checked the same way (no `realpath`, no stat of
  a target). `..` and absolute targets are walked as the kernel would. More than 40 hops is a symlink loop:
  `--scan DIR is a symlink loop` or `ERROR CHILD is a symlink loop`, rc 2.
- `git status` in a checked repo gets `-- ':(exclude,glob)**/NAME/**'` per NAME, so git prunes untracked or
  ignored NAME folders unread (checked with strace). Pathspec cannot stop git lstat'ing files the index TRACKS
  under NAME, so a repo whose `git ls-files` (index only) lists such a file skips its dirty-work check, said as
  `skip dirty check REPO: it tracks files under NAME`; its branch checks still run (they read `.git` only).
`pr_done_check.py --repo OWNER/NAME [--since YYYY-MM-DD] [--limit 1000]` (GH=<gh binary>; GET calls only).
`--max-age-min` and `--active-window-min` must be finite and >= 0 (`nan`, `inf` or a negative value gives rc 2).
`--include-file` entries are paths as written, resolved against the current directory: use absolute paths.
`pr_done_check.py` fails closed (rc 2) on any malformed gh output, including a list entry that is not an object.
`--since` must be a real `YYYY-MM-DD` date. It is compared with the **UTC** date of mergedAt / closedAt, so a PR
merged 2026-11-24 9 PM CST (03:00Z on 11-25) counts as 2026-11-25.

Known limits (Cassie N4, N5):
- git gives a pre-commit hook no way to tell an amend apart, so `git commit --amend` of the one unpushed commit is
  refused too. Push it first; then the amend is a new commit.
- Dirty-file age is the file mtime, so a `chmod` or `touch` of an old file counts as a change now. This errs
  toward failing.

All git calls use `--no-optional-locks` and drop `GIT_DIR`/`GIT_WORK_TREE`/`GIT_INDEX_FILE`/`GIT_COMMON_DIR`,
so the drift check only reads other agents' worktrees. Scratch goes to `${ZF_SCRATCH:-$HOME/.cache/zenith-float}`,
never `/tmp`. Enable the hooks once per clone with `git config core.hooksPath .githooks`. Fix a refusal; never
bypass it with `--no-verify`.

## Why the pre-push hook was rewritten (2:11 PM CT refusal of d41d938)

Reproduced in `tools/pre_push_selftest.sh` against #40's hook:
1. git exports `GIT_DIR=<gitdir>` to hooks. From a linked worktree it is an absolute path, so #40's
   `git -C <scratch clone> checkout` ran against the pusher's repo. The scratch clone stayed empty, and the hook
   reported "could not check out ... or it has no tools/no_suppress_scan.py". The new hook unsets the git env first.
2. #40's hook ran `rules_scan.sh` in the pusher's worktree. Its count checks read working-tree files, so dirty
   files decided the verdict for a clean commit. The new hook runs every scanner in the clean clone at the pushed sha.

#40 also adds `.githooks/pre-commit` (rules_scan, fmt, clippy). When both land, the merged hook must keep both: the
unpushed cap and #40's steps. `drift_check_selftest.sh` copies only the hook and `drift_check.py` into its
fixture, so update that test in the same merge.

## Tests

`tools/drift_check_selftest.sh`, `tools/pr_done_check_selftest.sh` (mock gh that honours `--state` and `--limit`),
`tools/pre_push_selftest.sh` (real pushes to a local bare origin, stub scanners), `tools/rule_checks_daily_selftest.sh`,
`tools/scoped_gate_selftest.sh` (stub poprun), and `tools/rule_checks_mutation.sh` (231 mutants, including Cassie's
16 round-1, 12 round-2, 25 round-3 and 20 round-4 rows; every one must be killed; a row's `\n` is a newline, and a
mutant that hangs is killed after 15 min). Cassie's round-4 X18 (`--exclude "$x"` -> `--exclude="$x"`) is
equivalent, not a row: argparse reads both forms the same. Rows whose code was rewritten in round 4 were adapted to
the new code with the same intent (D25, N03, Z13-Z15, W05, W15, W16, V10-V12, V16, and Cassie's X01, X04-X08,
X10, X15, X16, X20). Each behavior landed as one red commit, then one fix commit. Test-only commits (56878368,
50229306, bf438dca, b3b7a524, 89d7a4ec and the round-4 walk/git-status cases) add cases for code that was already
right; their red proof is the mutant each one kills, all in the sweep. full_gate.sh runs all six as step `rule-checks`.

## Scoped gate (in the repo)

`tools/scoped_gate.sh rules` is the scoped gate for a PR that changes only these tools (Jeff Q6: a scoped gate per PR,
the full gate at tier boundaries). Steps: `full-gate-syntax` (`bash -n tools/full_gate.sh`), `rules-scan` and
`rule-checks` (their commands are read verbatim from `tools/full_gate.sh`, so the two cannot drift),
`gate-tool-checks` and `full-gate-poplog`. Each step is one `poprun --step`; all steps run; rc 1 names the failed
ones. `--list` / `--dry-run` print the plan. Run it from a clean checkout of the PR head:

    PYTHONDONTWRITEBYTECODE=1 ~/pop-logs/bin/poprun --agent spence --job scoped-gate-rules -- bash tools/scoped_gate.sh rules

The poprun footer must say `RESULT: PASS` and `gate-commit: <PR head> (clean at start and end)`. It is not
full_gate evidence: it covers 2 of full_gate's 25 run steps (no Rust builds or tests).

## Daily full sweep (Jeff 2026-11-25 5:20 PM CT)

"Yes, add a daily full sweep over every worktree that fails on any drift." Active mode hides idle drift (Cassie N1),
so once a day `tools/rule_checks_daily.sh` runs `drift_check.py --all` (every tree and every branch, nothing skipped)
over every root and FAILs (rc 1) on any drift; rc 2 on no roots, a bad root or a git error.
- Roots are dirs or quoted globs of dirs: on Pop `/mnt/hub/{paul,rusty,spence,vera}-work` and
  `/mnt/hub/home-moved/*`. Non-dir matches are ignored, but a glob that matches no directory at all is an error
  (rc 2), so a mistyped root cannot silently drop out of the sweep. Wildcards are allowed in the last component
  only (rc 2 otherwise): expanding `x/*/mw` would stat `x/MathWorks/mw` before any name check.
- The same dir given twice (a literal and a glob, `./`, a symlink) is swept once, first spelling kept.
- Any root with a `MathWorks` path component, as written or after following its symlinks by name (the same
  lexical walk as drift_check: `[ -L ]` and `readlink`, never `readlink -f` or a stat of a target), is left out
  (named `excluded`). The pattern itself is judged before it is expanded, and each glob match before `[ -d ]`. A
  root that is a symlink loop gives `root PATH is a symlink loop`, rc 2. The same names go to `drift_check --exclude`, so a scanned child folder named `MathWorks` is left out
  too. Trees elsewhere whose git data lives in MathWorks are still read (read-only), like in the 15-minute run.
- The summary line counts the FAIL (finding) and ERROR lines, roots and exclusions.

## 15-minute scheduler on Pop

A systemd **user** timer, the same pattern as `poprotate.timer`. No cron, no system units.
- `~/.config/systemd/user/zf-rule-checks.timer`: `OnCalendar=*:0/15`, `Persistent=true`.
- `~/.config/systemd/user/zf-rule-checks.service`: oneshot, `ExecStart=/mnt/hub/spence-work/rule-checks-bin/run.sh`.
- `run.sh` runs copies of `drift_check.py` and `pr_done_check.py` pinned to one commit, kept in
  `/mnt/hub/spence-work/rule-checks-bin/<sha>/` with `SHA256SUMS`. It refuses to run if a hash does not match.
  It runs `drift_check.py --scan` (active mode) over `/mnt/hub/{paul,rusty,spence,vera}-work` with
  `--include-file rule-checks-bin/include.txt` and `--exclude MathWorks`, then `pr_done_check.py` on `jscarr64/zenith-float` and
  `jscarr64/zenith-bound`. Output is appended to `/home/jsc64/pop-logs/spence/rule-checks/YYYY-MM-DD.log`.
- Daily: `zf-rule-checks-daily.timer` (`OnCalendar=*-*-* 06:00`, `Persistent=true`) starts
  `zf-rule-checks-daily.service` (oneshot, `ExecStart=/mnt/hub/spence-work/rule-checks-bin/run-daily.sh`). It checks
  the same pinned `SHA256SUMS` and runs the pinned `rule_checks_daily.sh` over the roots above. Output:
  `/home/jsc64/pop-logs/spence/rule-checks/YYYY-MM-DD-daily.log`.
- To repin after a merge, copy `drift_check.py`, `pr_done_check.py` and `rule_checks_daily.sh` at the new sha into
  `rule-checks-bin/<sha>/`, write `SHA256SUMS`, and change `PIN` in `run.sh` and `run-daily.sh`. To stop them:
  `systemctl --user disable --now zf-rule-checks.timer zf-rule-checks-daily.timer`.
- Linger is on for jsc64 (Jeff 2026-11-25), so both timers run without a session on Pop.
